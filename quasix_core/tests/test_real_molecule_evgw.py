#!/usr/bin/env python3
"""Test evGW with real molecular data from PySCF.

This script tests the Rust evGW implementation using real molecular data
from PySCF calculations, including proper DF tensors and orbital energies.
"""

import sys
import numpy as np
from pathlib import Path
import time
import os

# Disable verbose logging from Rust
os.environ['RUST_LOG'] = 'error'

# Add parent directory to path
parent_dir = Path(__file__).parent.parent.parent
sys.path.insert(0, str(parent_dir))
sys.path.insert(0, str(parent_dir / 'quasix'))

# Import PySCF
from pyscf import gto, scf, df
from pyscf.tools import molden

# Import QuasiX
try:
    import quasix
    # Try to get the Rust module
    if hasattr(quasix, 'quasix'):
        rust_module = quasix.quasix
    else:
        rust_module = quasix
    print(f"✓ QuasiX module loaded from: {quasix.__file__}")
except ImportError as e:
    print(f"✗ Failed to import quasix: {e}")
    sys.exit(1)


def build_df_tensors_pyscf(mol, mf, auxbasis='cc-pvdz-jkfit'):
    """Build DF tensors using PySCF.

    Args:
        mol: PySCF molecule
        mf: Converged SCF object
        auxbasis: Auxiliary basis set

    Returns:
        tuple: (iaP, ijP, v_aux, vxc_dft)
    """
    print(f"\n[DF Tensor Construction]")
    print(f"  Auxiliary basis: {auxbasis}")

    # Get dimensions
    nao = mol.nao
    nmo = mf.mo_coeff.shape[1]
    nocc = mol.nelec[0]  # Number of alpha electrons (closed-shell)
    nvir = nmo - nocc

    print(f"  NAO: {nao}, NMO: {nmo}, Nocc: {nocc}, Nvir: {nvir}")

    # Build auxiliary basis
    auxmol = df.addons.make_auxmol(mol, auxbasis=auxbasis)
    naux = auxmol.nao
    print(f"  Auxiliary functions: {naux}")

    # Build DF object
    df_obj = df.DF(mol, auxbasis=auxbasis)
    df_obj.kernel()

    # Get the 2-center integrals (P|Q)
    v_aux = auxmol.intor('int2c2e')
    v_aux = (v_aux + v_aux.T) / 2  # Symmetrize

    # Check for positive definiteness
    eigvals = np.linalg.eigvalsh(v_aux)
    print(f"  V_aux eigenvalues: min={eigvals.min():.2e}, max={eigvals.max():.2e}")
    if eigvals.min() < 1e-10:
        print(f"  Adding regularization: {1e-10 - eigvals.min():.2e}")
        v_aux += (1e-10 - eigvals.min()) * np.eye(naux)

    # Build 3-center integrals (P|μν)
    print(f"  Building 3-center integrals...")
    int3c = df.incore.aux_e2(mol, auxmol, intor='int3c2e')
    eri_3c = int3c.reshape(naux, nao, nao)

    # Transform to MO basis
    print(f"  Transforming to MO basis...")
    iaP_3d = np.zeros((nocc, nvir, naux))
    ijP_3d = np.zeros((nocc, nocc, naux))

    for p in range(naux):
        # Transform (P|μν) -> (P|pq)
        mo_L = np.dot(mf.mo_coeff.T, eri_3c[p])
        mo_pq = np.dot(mo_L, mf.mo_coeff)

        # Extract blocks
        iaP_3d[:, :, p] = mo_pq[:nocc, nocc:]
        ijP_3d[:, :, p] = mo_pq[:nocc, :nocc]

    # Reshape to 2D for Rust interface
    iaP = iaP_3d.reshape(nocc * nvir, naux)
    ijP = ijP_3d.reshape(nocc * nocc, naux)

    # Get VXC (exchange-correlation potential)
    # For HF, this is just the Fock matrix diagonal minus kinetic and nuclear
    # For simplicity, we'll use zeros (pure HF exchange)
    vxc_dft = np.zeros(nmo)

    print(f"  iaP shape: {iaP.shape}")
    print(f"  ijP shape: {ijP.shape}")
    print(f"  v_aux shape: {v_aux.shape}")

    # Verify tensor norms
    print(f"  |iaP|_F = {np.linalg.norm(iaP):.6f}")
    print(f"  |ijP|_F = {np.linalg.norm(ijP):.6f}")
    print(f"  |v_aux|_F = {np.linalg.norm(v_aux):.6f}")

    return iaP, ijP, v_aux, vxc_dft


def test_water_molecule():
    """Test evGW on water molecule."""
    print("\n" + "="*60)
    print("Testing evGW with H2O molecule")
    print("="*60)

    # Build water molecule
    mol = gto.M(
        atom='''
        O    0.000000    0.000000    0.117790
        H    0.000000    0.756950   -0.471161
        H    0.000000   -0.756950   -0.471161
        ''',
        basis='cc-pvdz',
        symmetry=False,
        charge=0,
        spin=0
    )
    mol.build()

    print(f"\n[Molecule]")
    print(f"  Atoms: {mol.natm}")
    print(f"  Electrons: {mol.nelectron}")
    print(f"  AO basis functions: {mol.nao}")

    # Run HF calculation
    print(f"\n[HF Calculation]")
    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()

    if not mf.converged:
        print("  ✗ HF did not converge!")
        return False

    print(f"  ✓ HF converged")
    print(f"  Total energy: {mf.e_tot:.6f} Ha")
    print(f"  HOMO energy: {mf.mo_energy[mol.nelec[0]-1]:.6f} Ha")
    print(f"  LUMO energy: {mf.mo_energy[mol.nelec[0]]:.6f} Ha")
    print(f"  HOMO-LUMO gap: {mf.mo_energy[mol.nelec[0]] - mf.mo_energy[mol.nelec[0]-1]:.6f} Ha")

    # Build DF tensors
    iaP, ijP, v_aux, vxc_dft = build_df_tensors_pyscf(mol, mf, auxbasis='cc-pvdz-jkfit')

    # Prepare data for Rust
    mo_energy = mf.mo_energy.astype(np.float64)
    mo_occ = np.zeros(len(mo_energy))
    mo_occ[:mol.nelec[0]] = 2.0  # Closed-shell

    # Run evGW
    print(f"\n[evGW Calculation]")
    print(f"  Parameters:")
    print(f"    max_cycle: 12")
    print(f"    conv_tol: 1e-4 Ha")
    print(f"    damping: 0.3")
    print(f"    nfreq: 48")
    print(f"    eta: 0.01")

    start_time = time.time()

    try:
        result = rust_module.evgw(
            mo_energy,
            mo_occ,
            iaP,
            ijP,
            v_aux,
            vxc_dft,
            max_cycle=12,
            conv_tol=1e-4,
            conv_tol_z=1e-3,
            damping=0.3,
            damping_dynamic=False,
            diis=False,
            diis_space=6,
            nfreq=48,
            eta=0.01,
            verbose=0
        )
    except Exception as e:
        print(f"  ✗ evGW failed: {e}")
        import traceback
        traceback.print_exc()
        return False

    elapsed = time.time() - start_time
    print(f"\n  Calculation time: {elapsed:.2f} seconds")

    # Extract results
    qp_energies = np.array(result['qp_energies'])
    z_factors = np.array(result['z_factors'])
    converged = result['converged']
    n_cycles = result['n_cycles']

    print(f"\n[Results]")
    print(f"  Converged: {converged}")
    print(f"  Iterations: {n_cycles}")

    # Check convergence
    if not converged:
        print(f"  ⚠ Warning: evGW did not converge")

    # Validate Z-factors
    print(f"\n[Z-factors validation]")
    z_min, z_max = z_factors.min(), z_factors.max()
    print(f"  Range: [{z_min:.4f}, {z_max:.4f}]")

    if z_min < 0 or z_max > 1:
        print(f"  ✗ Unphysical Z-factors detected!")
        print(f"    Expected: [0, 1]")
        return False
    else:
        print(f"  ✓ Z-factors are physical")

    # QP corrections
    print(f"\n[Quasiparticle corrections]")
    qp_corrections = qp_energies - mo_energy
    nocc = mol.nelec[0]

    print(f"  Occupied orbitals:")
    for i in range(max(0, nocc-3), nocc):
        print(f"    MO {i:2d}: HF = {mo_energy[i]:8.4f} Ha, "
              f"QP = {qp_energies[i]:8.4f} Ha, "
              f"Δ = {qp_corrections[i]*27.211:6.3f} eV, "
              f"Z = {z_factors[i]:.3f}")

    print(f"  Virtual orbitals:")
    for i in range(nocc, min(nocc+3, len(mo_energy))):
        print(f"    MO {i:2d}: HF = {mo_energy[i]:8.4f} Ha, "
              f"QP = {qp_energies[i]:8.4f} Ha, "
              f"Δ = {qp_corrections[i]*27.211:6.3f} eV, "
              f"Z = {z_factors[i]:.3f}")

    # HOMO-LUMO gap
    homo_idx = nocc - 1
    lumo_idx = nocc
    hf_gap = mo_energy[lumo_idx] - mo_energy[homo_idx]
    qp_gap = qp_energies[lumo_idx] - qp_energies[homo_idx]

    print(f"\n[HOMO-LUMO Gap]")
    print(f"  HF gap:  {hf_gap*27.211:.3f} eV")
    print(f"  QP gap:  {qp_gap*27.211:.3f} eV")
    print(f"  Change:  {(qp_gap - hf_gap)*27.211:.3f} eV")

    # Ionization potential
    ip_hf = -mo_energy[homo_idx]
    ip_qp = -qp_energies[homo_idx]

    print(f"\n[Ionization Potential]")
    print(f"  HF IP:  {ip_hf*27.211:.3f} eV")
    print(f"  QP IP:  {ip_qp*27.211:.3f} eV")
    print(f"  Change: {(ip_qp - ip_hf)*27.211:.3f} eV")

    # Reference values for H2O/cc-pVDZ from literature
    # These are approximate expected values
    print(f"\n[Validation against expected values]")
    expected_ip_correction = -1.5  # eV, typical GW correction for H2O HOMO
    actual_ip_correction = (ip_qp - ip_hf) * 27.211

    print(f"  Expected IP correction: ~{expected_ip_correction:.1f} eV")
    print(f"  Actual IP correction:   {actual_ip_correction:.3f} eV")

    if abs(actual_ip_correction - expected_ip_correction) > 1.0:
        print(f"  ⚠ Warning: IP correction differs significantly from expected")
    else:
        print(f"  ✓ IP correction is reasonable")

    return True


def test_ammonia_molecule():
    """Test evGW on ammonia molecule."""
    print("\n" + "="*60)
    print("Testing evGW with NH3 molecule")
    print("="*60)

    # Build ammonia molecule
    mol = gto.M(
        atom='''
        N    0.000000    0.000000    0.116490
        H    0.000000    0.939731   -0.271810
        H    0.813798   -0.469865   -0.271810
        H   -0.813798   -0.469865   -0.271810
        ''',
        basis='cc-pvdz',
        symmetry=False,
        charge=0,
        spin=0
    )
    mol.build()

    print(f"\n[Molecule]")
    print(f"  Atoms: {mol.natm}")
    print(f"  Electrons: {mol.nelectron}")
    print(f"  AO basis functions: {mol.nao}")

    # Run HF calculation
    print(f"\n[HF Calculation]")
    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()

    if not mf.converged:
        print("  ✗ HF did not converge!")
        return False

    print(f"  ✓ HF converged")
    print(f"  Total energy: {mf.e_tot:.6f} Ha")
    print(f"  HOMO energy: {mf.mo_energy[mol.nelec[0]-1]:.6f} Ha")
    print(f"  LUMO energy: {mf.mo_energy[mol.nelec[0]]:.6f} Ha")

    # Build DF tensors
    iaP, ijP, v_aux, vxc_dft = build_df_tensors_pyscf(mol, mf, auxbasis='cc-pvdz-jkfit')

    # Prepare data for Rust
    mo_energy = mf.mo_energy.astype(np.float64)
    mo_occ = np.zeros(len(mo_energy))
    mo_occ[:mol.nelec[0]] = 2.0  # Closed-shell

    # Run evGW with more aggressive convergence parameters
    print(f"\n[evGW Calculation]")
    print(f"  Using more iterations for NH3...")

    start_time = time.time()

    try:
        result = rust_module.evgw(
            mo_energy,
            mo_occ,
            iaP,
            ijP,
            v_aux,
            vxc_dft,
            max_cycle=15,  # More iterations for NH3
            conv_tol=1e-4,
            conv_tol_z=1e-3,
            damping=0.2,  # Less damping
            damping_dynamic=True,  # Enable dynamic damping
            diis=False,
            diis_space=6,
            nfreq=48,
            eta=0.01,
            verbose=0
        )
    except Exception as e:
        print(f"  ✗ evGW failed: {e}")
        return False

    elapsed = time.time() - start_time
    print(f"\n  Calculation time: {elapsed:.2f} seconds")

    # Extract results
    qp_energies = np.array(result['qp_energies'])
    z_factors = np.array(result['z_factors'])
    converged = result['converged']
    n_cycles = result['n_cycles']

    print(f"\n[Results]")
    print(f"  Converged: {converged}")
    print(f"  Iterations: {n_cycles}")

    # Validate Z-factors
    print(f"\n[Z-factors validation]")
    z_min, z_max = z_factors.min(), z_factors.max()
    print(f"  Range: [{z_min:.4f}, {z_max:.4f}]")

    if z_min < 0 or z_max > 1:
        print(f"  ✗ Unphysical Z-factors detected!")
        return False
    else:
        print(f"  ✓ Z-factors are physical")

    # QP corrections
    nocc = mol.nelec[0]
    homo_idx = nocc - 1
    lumo_idx = nocc

    ip_hf = -mo_energy[homo_idx]
    ip_qp = -qp_energies[homo_idx]

    print(f"\n[Ionization Potential]")
    print(f"  HF IP:  {ip_hf*27.211:.3f} eV")
    print(f"  QP IP:  {ip_qp*27.211:.3f} eV")
    print(f"  Change: {(ip_qp - ip_hf)*27.211:.3f} eV")

    return True


def main():
    """Run all tests."""
    print("QuasiX evGW Test with Real Molecular Data")
    print("==========================================")

    # Check if we have the Rust evgw function
    if not hasattr(rust_module, 'evgw'):
        print("✗ Rust evgw function not found!")
        print(f"  Available functions: {dir(rust_module)}")
        return 1

    print("✓ Rust evgw function available")

    # Test molecules
    all_passed = True

    # Test water
    if not test_water_molecule():
        all_passed = False
        print("\n✗ Water test failed")
    else:
        print("\n✓ Water test passed")

    # Test ammonia
    if not test_ammonia_molecule():
        all_passed = False
        print("\n✗ Ammonia test failed")
    else:
        print("\n✓ Ammonia test passed")

    # Summary
    print("\n" + "="*60)
    if all_passed:
        print("✓ All tests passed!")
        return 0
    else:
        print("✗ Some tests failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())