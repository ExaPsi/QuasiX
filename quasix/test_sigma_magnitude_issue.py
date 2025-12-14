#!/usr/bin/env python3
"""
Diagnostic script to identify the 10x magnitude error in Σ_c
"""
import numpy as np
import pyscf
from pyscf import gto, scf, gw
import quasix

def test_h2o_sigma_c():
    """Compare Rust vs PySCF Σ_c for H2O"""

    # Build H2O molecule
    mol = gto.M(
        atom='O 0 0 0; H 0 0 0.96; H 0.96 0 0',
        basis='def2-svp',
        verbose=0
    )
    mf = scf.RHF(mol).run(verbose=0)

    n_mo = mf.mo_coeff.shape[1]
    n_occ = np.sum(mf.mo_occ > 0)

    print(f"\n{'='*60}")
    print(f"H2O molecule: {n_occ} occupied, {n_mo-n_occ} virtual orbitals")
    print(f"{'='*60}\n")

    # Get DF tensors using PySCF
    from pyscf import df
    auxbasis = df.make_auxbasis(mol, mp2fit=True)
    with_df = df.DF(mol, auxbasis=auxbasis)
    with_df.build()

    # Get 3-center integrals (mn|P) in MO basis
    mo_coeff = mf.mo_coeff
    naux = with_df.get_naoaux()

    # Transform to MO basis
    df_tensor = np.zeros((n_mo, n_mo, naux))
    for p in range(naux):
        eri_ao = with_df._cderi[p].reshape(mol.nao, mol.nao)
        df_tensor[:, :, p] = mo_coeff.T @ eri_ao @ mo_coeff

    # Compute Coulomb metric
    j2c = with_df.get_2c2e()
    v_sqrt = np.linalg.cholesky(j2c)

    # Compute P0 for a test imaginary frequency
    xi = 1.0  # Test imaginary frequency

    nocc = n_occ
    nvirt = n_mo - n_occ
    nfreq = 16

    print(f"Computing P0 at ω = i·{xi} Ha")
    print(f"  nocc = {nocc}, nvirt = {nvirt}, naux = {naux}")

    # Compute P0 using Rust
    try:
        # Extract DF tensor for occ-virt transitions
        df_ia = df_tensor[:nocc, nocc:, :].reshape(nocc * nvirt, naux)
        e_occ = mf.mo_energy[:nocc]
        e_virt = mf.mo_energy[nocc:]

        # Call Rust P0 calculator
        omega_complex = complex(0.0, xi)
        p0_rust = quasix.compute_polarizability_p0(
            omega=omega_complex,
            df_ia=df_ia,
            e_occ=e_occ,
            e_virt=e_virt,
            eta=1e-4
        )

        print(f"\nP0 computed (Rust):")
        print(f"  Shape: {p0_rust.shape}")
        print(f"  Diagonal mean: {np.mean(np.diag(p0_rust.real)):.6e}")
        print(f"  Trace: {np.trace(p0_rust.real):.6e}")
        print(f"  Frobenius norm: {np.linalg.norm(p0_rust):.6e}")

        # Check spin factor
        # For RHF, P0 should have factor of 2 for spin
        # P0_PQ = 2 * Σ_ia (ia|P)(ia|Q) / (ε_a - ε_i - iξ)

        # Manual computation to verify
        p0_manual = np.zeros((naux, naux), dtype=complex)
        for i in range(nocc):
            for a in range(nvirt):
                de = e_virt[a] - e_occ[i]
                denom = 2.0 * de / (de**2 + xi**2)  # Real part for imaginary frequency
                ia_idx = i * nvirt + a
                for p in range(naux):
                    for q in range(naux):
                        p0_manual[p, q] += denom * df_ia[ia_idx, p] * df_ia[ia_idx, q]

        print(f"\nP0 computed (manual Python):")
        print(f"  Shape: {p0_manual.shape}")
        print(f"  Diagonal mean: {np.mean(np.diag(p0_manual.real)):.6e}")
        print(f"  Trace: {np.trace(p0_manual.real):.6e}")
        print(f"  Frobenius norm: {np.linalg.norm(p0_manual):.6e}")

        # Compare
        diff = np.linalg.norm(p0_rust - p0_manual) / np.linalg.norm(p0_manual)
        print(f"\nRelative difference (Rust vs manual): {diff:.6e}")

        if diff > 1e-10:
            print("\n⚠️  WARNING: P0 differs between Rust and manual Python!")
            print("   This suggests an error in the Rust P0 calculation")
        else:
            print("\n✓ P0 matches between Rust and manual Python")

        # Now compute W = v^{1/2} [1 - M]^{-1} v^{1/2} where M = v^{1/2} P0 v^{1/2}
        print(f"\n{'='*60}")
        print("Computing screened interaction W")
        print(f"{'='*60}\n")

        # Compute M
        v_sqrt_inv = np.linalg.inv(v_sqrt)
        m_matrix = v_sqrt @ p0_rust @ v_sqrt.T

        print(f"M = v^{{1/2}} P0 v^{{1/2}}:")
        print(f"  Diagonal mean: {np.mean(np.diag(m_matrix.real)):.6e}")
        print(f"  Max eigenvalue: {np.max(np.linalg.eigvalsh(m_matrix.real)):.6f}")

        # Invert (1 - M)
        identity = np.eye(naux)
        one_minus_m = identity - m_matrix
        inv_one_minus_m = np.linalg.inv(one_minus_m)

        # Compute W
        w_matrix = v_sqrt @ inv_one_minus_m @ v_sqrt.T

        print(f"\nW = v^{{1/2}} [1-M]^{{-1}} v^{{1/2}}:")
        print(f"  Diagonal mean: {np.mean(np.diag(w_matrix.real)):.6e}")
        print(f"  Frobenius norm: {np.linalg.norm(w_matrix):.6e}")

        # Compare with bare Coulomb
        v_matrix = v_sqrt @ v_sqrt.T
        print(f"\nBare Coulomb v = v^{{1/2}} v^{{1/2}}:")
        print(f"  Diagonal mean: {np.mean(np.diag(v_matrix)):.6e}")
        print(f"  Frobenius norm: {np.linalg.norm(v_matrix):.6e}")

        ratio = np.mean(np.diag(w_matrix.real)) / np.mean(np.diag(v_matrix))
        print(f"\nW/v ratio: {ratio:.3f}")

        if ratio < 1.0:
            print("⚠️  WARNING: W < v is unphysical! Should have W >= v (screening amplifies)")
        elif ratio > 10.0:
            print("⚠️  WARNING: W >> v suggests over-amplification")
        else:
            print("✓ W/v ratio looks reasonable")

        # Now estimate Σ_c magnitude
        # For a rough estimate, use Σ_c ~ W * G ~ W / (ε_HOMO)
        homo_energy = mf.mo_energy[nocc - 1]
        w_diag_mean = np.mean(np.diag(w_matrix.real))

        sigma_c_estimate = w_diag_mean / abs(homo_energy)
        print(f"\n{'='*60}")
        print(f"Estimated Σ_c magnitude:")
        print(f"{'='*60}")
        print(f"  HOMO energy: {homo_energy:.6f} Ha")
        print(f"  W diagonal mean: {w_diag_mean:.6e}")
        print(f"  Rough estimate: Σ_c ~ W/|ε_HOMO| ~ {sigma_c_estimate:.6e} Ha")
        print(f"  Expected range: 0.1 - 1.0 eV ~ 0.004 - 0.04 Ha")

    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_h2o_sigma_c()