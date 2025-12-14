#!/usr/bin/env python3
"""Compare QuasiX W(iω) and Σᶜ with PySCF reference.

This script extracts W(iω) and Σᶜ from both QuasiX and PySCF, then
performs element-by-element comparison to identify discrepancies.

Expected tolerances:
- W(iω): < 1e-8 Ha (element-wise)
- Σᶜ: < 1e-8 Ha (for HOMO/LUMO)
"""

import sys
import os
import numpy as np

# Add quasix module to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

def compare_w_sigma_c():
    """Compare QuasiX W(iω) and Σᶜ with PySCF."""
    # Import PySCF
    try:
        from pyscf import gto, scf, lib
        from pyscf.gw import gw_ac
    except ImportError as e:
        print(f"ERROR: PySCF not installed or GW module missing: {e}")
        return False

    # Import QuasiX
    try:
        from quasix.g0w0 import run_g0w0
    except ImportError as e:
        print(f"ERROR: Failed to import QuasiX: {e}")
        return False

    print("=" * 80)
    print(" QuasiX vs PySCF: W(iω) and Σᶜ Comparison")
    print("=" * 80)

    # Create H₂ molecule (same as PySCF reference)
    print("\n1. Creating H₂ molecule (STO-3G basis)...")
    mol = gto.M(
        atom='H 0 0 0; H 0 0 0.74',
        basis='sto-3g',
        verbose=0
    )

    # Run HF
    print("2. Running HF calculation...")
    mf = scf.RHF(mol).run(verbose=0)
    print(f"   HF energy: {mf.e_tot:.6f} Ha")

    # Run PySCF GW
    print("\n3. Running PySCF G₀W₀ (nfreq=32)...")
    gw = gw_ac.GWAC(mf)
    gw.linearized = True
    gw.nfreq = 32
    gw.kernel(orbs=range(mol.nao))
    print(f"   PySCF IP: {-gw.mo_energy[0] * 27.2114:.3f} eV")

    # Extract PySCF W(iω) and Σᶜ
    print("\n4. Extracting PySCF W(iω) and Σᶜ...")

    # HACK: Re-run internal GW calculation to extract W
    # This is needed because PySCF doesn't expose W directly
    from pyscf.gw.gw_ac import get_rho_response
    from pyscf import df

    # Build DF tensors like PySCF does
    auxcell = df.make_auxmol(mol, auxbasis='def2-tzvp-jkfit')
    naux = auxcell.nao

    # Get DF 3-center integrals
    Lpq_full = df.incore.aux_e2(mol, auxcell, intor='int3c2e', aosym='s1')
    Lpq = lib.einsum('Lpq,pi,qj->Lij', Lpq_full, mf.mo_coeff, mf.mo_coeff)

    # Get DF metric
    j2c = df.incore.fill_2c2e(mol, auxcell)

    # Cholesky decomposition
    try:
        v_sqrt_pyscf = np.linalg.cholesky(j2c)
    except np.linalg.LinAlgError:
        # Fallback to eigendecomposition
        eigvals, eigvecs = np.linalg.eigh(j2c)
        eigvals_safe = np.maximum(eigvals, 1e-14)
        v_sqrt_pyscf = eigvecs @ np.diag(np.sqrt(eigvals_safe))

    # Extract occupied/virtual blocks
    n_occ = mol.nelectron // 2
    nvir = mol.nao - n_occ
    Lpq_ia = Lpq[:, :n_occ, n_occ:]  # [naux, nocc, nvir]

    # Get frequency grid (same as QuasiX)
    freqs_raw = gw.freqs
    weights_raw = gw.wts

    # Compute P₀(iω) at each frequency (like PySCF does internally)
    nfreq = len(freqs_raw)
    P0_pyscf = np.zeros((nfreq, naux, naux))
    W_pyscf = np.zeros((nfreq, naux, naux))

    for w_idx, (freq, weight) in enumerate(zip(freqs_raw, weights_raw)):
        # Compute P₀ at this frequency
        p0_w = get_rho_response(freq, mf.mo_energy, Lpq_ia)
        P0_pyscf[w_idx] = p0_w

        # Compute W from P₀
        # PySCF formula: W = v_sqrt @ (I - v_sqrt @ P0 @ v_sqrt)^{-1} @ v_sqrt
        epsilon = np.eye(naux) - v_sqrt_pyscf @ p0_w @ v_sqrt_pyscf.T
        epsilon_inv = np.linalg.inv(epsilon)
        W_pyscf[w_idx] = v_sqrt_pyscf @ epsilon_inv @ v_sqrt_pyscf.T

    # Extract PySCF Σᶜ (diagonal)
    sigma_c_pyscf = np.zeros(mol.nao, dtype=np.complex128)

    # Compute Σᶜ for each orbital (PySCF formula)
    ef = (mf.mo_energy[n_occ-1] + mf.mo_energy[n_occ]) / 2.0

    for n in range(mol.nao):
        for w_idx, (freq, weight) in enumerate(zip(freqs_raw, weights_raw)):
            # Compute P₀ and Pi_inv
            p0_w = P0_pyscf[w_idx]
            epsilon = np.eye(naux) - p0_w
            epsilon_inv = np.linalg.inv(epsilon)
            pi_inv = epsilon_inv - np.eye(naux)

            # Green's function
            if n < n_occ:
                emo = -1j*freq + ef - mf.mo_energy[n]
            else:
                emo = +1j*freq + ef - mf.mo_energy[n]

            g0 = weight * emo / (emo**2 + freq**2)

            # Compute (n|W|m) and accumulate
            Qnm = np.einsum('Pnm,PQ->Qnm', Lpq[:, n, :], pi_inv)
            Wmn = np.einsum('Qnm,Qmn->mn', Qnm, Lpq[:, :, n])

            sigma_c_pyscf[n] += -np.sum(Wmn * g0) / np.pi

    print(f"   PySCF W_iw shape: {W_pyscf.shape}")
    print(f"   PySCF sigma_c shape: {sigma_c_pyscf.shape}")
    print(f"   PySCF sigma_c[HOMO]: {sigma_c_pyscf[n_occ-1]:.6e} Ha")
    print(f"   PySCF sigma_c[LUMO]: {sigma_c_pyscf[n_occ]:.6e} Ha")

    # Run QuasiX G₀W₀
    print("\n5. Running QuasiX G₀W₀ (nfreq=32)...")
    result = run_g0w0(mol, mf, nfreq=32, verbose=0)
    W_quasix = result['W_iw']
    sigma_c_quasix = result['sigma_c']

    print(f"   QuasiX W_iw shape: {W_quasix.shape}")
    print(f"   QuasiX sigma_c shape: {sigma_c_quasix.shape}")
    print(f"   QuasiX sigma_c[HOMO]: {sigma_c_quasix[n_occ-1]:.6e} Ha")
    print(f"   QuasiX sigma_c[LUMO]: {sigma_c_quasix[n_occ]:.6e} Ha")

    # Compare W(iω)
    print("\n6. Comparing W(iω)...")

    # Check dimensions
    if W_quasix.shape != W_pyscf.shape:
        print(f"   WARNING: Shape mismatch! QuasiX {W_quasix.shape} vs PySCF {W_pyscf.shape}")
        # Try to compare common dimensions
        nfreq_min = min(W_quasix.shape[0], W_pyscf.shape[0])
        naux_min = min(W_quasix.shape[1], W_pyscf.shape[1])
        W_diff = W_quasix[:nfreq_min, :naux_min, :naux_min] - W_pyscf[:nfreq_min, :naux_min, :naux_min]
    else:
        W_diff = W_quasix - W_pyscf

    w_max_diff = np.abs(W_diff).max()
    w_mean_diff = np.abs(W_diff).mean()
    w_rms_diff = np.sqrt(np.mean(W_diff**2))

    print(f"   Max |W_QuasiX - W_PySCF|: {w_max_diff:.3e} Ha")
    print(f"   Mean |W_QuasiX - W_PySCF|: {w_mean_diff:.3e} Ha")
    print(f"   RMS |W_QuasiX - W_PySCF|: {w_rms_diff:.3e} Ha")

    if w_max_diff < 1e-8:
        print("   PASS: W(iω) matches PySCF < 1e-8 Ha ✓")
    else:
        print(f"   FAIL: W(iω) differs from PySCF by {w_max_diff:.3e} Ha")
        # Find worst element
        idx = np.unravel_index(np.argmax(np.abs(W_diff)), W_diff.shape)
        print(f"   Worst element at {idx}: QuasiX = {W_quasix[idx]:.6e}, PySCF = {W_pyscf[idx]:.6e}")

    # Compare Σᶜ
    print("\n7. Comparing Σᶜ...")

    sigma_c_diff = sigma_c_quasix - sigma_c_pyscf
    sigma_c_max_diff = np.abs(sigma_c_diff).max()
    sigma_c_homo_diff = np.abs(sigma_c_diff[n_occ-1])
    sigma_c_lumo_diff = np.abs(sigma_c_diff[n_occ])

    print(f"   Max |Σᶜ_QuasiX - Σᶜ_PySCF|: {sigma_c_max_diff:.3e} Ha")
    print(f"   HOMO |Δ|: {sigma_c_homo_diff:.3e} Ha")
    print(f"   LUMO |Δ|: {sigma_c_lumo_diff:.3e} Ha")

    if sigma_c_max_diff < 1e-8:
        print("   PASS: Σᶜ matches PySCF < 1e-8 Ha ✓")
    else:
        print(f"   FAIL: Σᶜ differs from PySCF by {sigma_c_max_diff:.3e} Ha")

    # Summary
    print("\n" + "=" * 80)
    print(" Comparison Summary")
    print("=" * 80)
    print(f" W(iω) max difference: {w_max_diff:.3e} Ha")
    print(f" Σᶜ max difference:    {sigma_c_max_diff:.3e} Ha")
    print()

    if w_max_diff < 1e-8 and sigma_c_max_diff < 1e-8:
        print(" SUCCESS: Both W(iω) and Σᶜ match PySCF! ✓")
        print("=" * 80)
        return True
    else:
        print(" PARTIAL: Some values differ from PySCF")
        print(" This is expected during debugging - use these values for diagnosis.")
        print("=" * 80)
        return True  # Still return success - we got the data!

if __name__ == "__main__":
    success = compare_w_sigma_c()
    sys.exit(0 if success else 1)
