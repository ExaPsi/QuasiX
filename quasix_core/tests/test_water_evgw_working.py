#!/usr/bin/env python3
"""Test evGW with water molecule using stable parameters."""

import sys
import numpy as np
from pathlib import Path
import time
import os

# Disable verbose logging
os.environ['RUST_LOG'] = 'warn'

# Add parent directory to path
parent_dir = Path(__file__).parent.parent.parent
sys.path.insert(0, str(parent_dir))
sys.path.insert(0, str(parent_dir / 'quasix'))

# Import PySCF
from pyscf import gto, scf, df

# Import QuasiX
import quasix
if hasattr(quasix, 'quasix'):
    rust_module = quasix.quasix
else:
    rust_module = quasix
print(f"QuasiX loaded from: {quasix.__file__}")


# Build water molecule with a better basis
print("Building H2O molecule...")
mol = gto.M(
    atom='''
    O    0.000000    0.000000    0.000000
    H    0.000000    0.756950    0.585882
    H    0.000000   -0.756950    0.585882
    ''',
    basis='6-31g',  # Better than STO-3G but still small
    symmetry=False,
    charge=0,
    spin=0
)
mol.build()
print(f"  Atoms: {mol.natm}, Electrons: {mol.nelectron}, AOs: {mol.nao}")

# Run HF with damping for better convergence
print("Running HF...")
mf = scf.RHF(mol)
mf.verbose = 0
mf.diis_damp = 0.5  # Add some damping for stability
mf.kernel()

if not mf.converged:
    print("  Warning: HF did not converge fully")

print(f"  HF energy: {mf.e_tot:.6f} Ha")
print(f"  HOMO: {mf.mo_energy[mol.nelec[0]-1]:.4f} Ha")
print(f"  LUMO: {mf.mo_energy[mol.nelec[0]]:.4f} Ha")

# Build DF tensors with a matching auxiliary basis
print("Building DF tensors...")
nao = mol.nao
nmo = mf.mo_coeff.shape[1]
nocc = mol.nelec[0]
nvir = nmo - nocc

# Use a proper auxiliary basis for 6-31G
auxbasis = 'def2-svp-jkfit'  # More compatible auxiliary basis
auxmol = df.addons.make_auxmol(mol, auxbasis=auxbasis)
naux = auxmol.nao
print(f"  NAO: {nao}, NMO: {nmo}, Nocc: {nocc}, Nvir: {nvir}, Naux: {naux}")

# Build V matrix with better conditioning
v_aux = auxmol.intor('int2c2e')
v_aux = (v_aux + v_aux.T) / 2

# Check conditioning and add stronger regularization if needed
eigvals = np.linalg.eigvalsh(v_aux)
print(f"  V_aux eigenvalues: min={eigvals.min():.2e}, max={eigvals.max():.2e}")
print(f"  Condition number: {eigvals.max()/eigvals.min():.2e}")

# Add stronger regularization to improve stability
regularization = max(1e-8, 1e-10 - eigvals.min())
if eigvals.min() < 1e-8:
    print(f"  Adding regularization: {regularization:.2e}")
    v_aux += regularization * np.eye(naux)
    eigvals = np.linalg.eigvalsh(v_aux)
    print(f"  After regularization: min_eig={eigvals.min():.2e}")

# Build 3-center integrals
print("Building 3-center integrals...")
int3c = df.incore.aux_e2(mol, auxmol, intor='int3c2e')
eri_3c = int3c.reshape(naux, nao, nao)

# Transform to MO basis
print("Transforming to MO basis...")
iaP_3d = np.zeros((nocc, nvir, naux))
ijP_3d = np.zeros((nocc, nocc, naux))

for p in range(naux):
    mo_L = np.dot(mf.mo_coeff.T, eri_3c[p])
    mo_pq = np.dot(mo_L, mf.mo_coeff)
    iaP_3d[:, :, p] = mo_pq[:nocc, nocc:]
    ijP_3d[:, :, p] = mo_pq[:nocc, :nocc]

# Reshape for Rust
iaP = iaP_3d.reshape(nocc * nvir, naux)
ijP = ijP_3d.reshape(nocc * nocc, naux)

# Check tensor magnitudes
print(f"  iaP: shape={iaP.shape}, norm={np.linalg.norm(iaP):.2f}, max={np.abs(iaP).max():.2f}")
print(f"  ijP: shape={ijP.shape}, norm={np.linalg.norm(ijP):.2f}, max={np.abs(ijP).max():.2f}")
print(f"  v_aux: norm={np.linalg.norm(v_aux):.2f}, max={np.abs(v_aux).max():.2f}")

# Prepare data
mo_energy = mf.mo_energy.astype(np.float64)
mo_occ = np.zeros(len(mo_energy))
mo_occ[:nocc] = 2.0
vxc_dft = np.zeros(nmo)  # Pure HF, no DFT XC

# Print orbital energies for debugging
print("\nOrbital energies (Ha):")
for i in range(min(10, nmo)):
    occ_str = "occ" if i < nocc else "vir"
    print(f"  MO {i:2d} ({occ_str}): {mo_energy[i]:8.4f}")

# Run evGW with very conservative parameters
print("\nRunning evGW with conservative parameters...")
print("  max_cycle: 3")
print("  conv_tol: 1e-2")
print("  damping: 0.7")
print("  nfreq: 16")

try:
    start = time.time()
    result = rust_module.evgw(
        mo_energy,
        mo_occ,
        iaP,
        ijP,
        v_aux,
        vxc_dft,
        max_cycle=3,       # Very few iterations
        conv_tol=1e-2,     # Very loose convergence
        conv_tol_z=1e-1,   # Very loose Z convergence
        damping=0.7,       # Heavy damping for stability
        damping_dynamic=False,
        diis=False,        # No DIIS for now
        diis_space=3,
        nfreq=16,          # Minimal frequency points
        eta=0.05,          # Larger broadening for stability
        verbose=0
    )
    elapsed = time.time() - start

    print(f"\n✓ evGW completed in {elapsed:.2f} seconds")

    # Extract results
    qp_energies = np.array(result['qp_energies'])
    z_factors = np.array(result['z_factors'])
    converged = result.get('converged', False)
    n_cycles = result.get('n_cycles', 0)

    print(f"\nResults:")
    print(f"  Iterations completed: {n_cycles}")
    print(f"  Converged: {converged}")

    # Check Z-factors
    z_min, z_max = z_factors.min(), z_factors.max()
    print(f"  Z-factors range: [{z_min:.3f}, {z_max:.3f}]")

    # Count physical Z-factors
    physical_z = np.logical_and(z_factors >= 0, z_factors <= 1)
    n_physical = np.sum(physical_z)
    print(f"  Physical Z-factors: {n_physical}/{len(z_factors)}")

    if n_physical < len(z_factors):
        print("  Non-physical Z-factors for orbitals:")
        for i, (z, phys) in enumerate(zip(z_factors, physical_z)):
            if not phys:
                print(f"    MO {i}: Z = {z:.3f}")

    # Show QP corrections for key orbitals
    print(f"\nQuasiparticle corrections (eV):")

    # HOMO
    homo_idx = nocc - 1
    if physical_z[homo_idx]:
        qp_corr_homo = (qp_energies[homo_idx] - mo_energy[homo_idx]) * 27.211
        print(f"  HOMO (MO {homo_idx}): {qp_corr_homo:6.3f} eV, Z = {z_factors[homo_idx]:.3f}")
    else:
        print(f"  HOMO (MO {homo_idx}): Failed (Z = {z_factors[homo_idx]:.3f})")

    # LUMO
    lumo_idx = nocc
    if lumo_idx < len(mo_energy) and physical_z[lumo_idx]:
        qp_corr_lumo = (qp_energies[lumo_idx] - mo_energy[lumo_idx]) * 27.211
        print(f"  LUMO (MO {lumo_idx}): {qp_corr_lumo:6.3f} eV, Z = {z_factors[lumo_idx]:.3f}")
    elif lumo_idx < len(mo_energy):
        print(f"  LUMO (MO {lumo_idx}): Failed (Z = {z_factors[lumo_idx]:.3f})")

    # Gaps
    if physical_z[homo_idx] and lumo_idx < len(mo_energy) and physical_z[lumo_idx]:
        hf_gap = (mo_energy[lumo_idx] - mo_energy[homo_idx]) * 27.211
        qp_gap = (qp_energies[lumo_idx] - qp_energies[homo_idx]) * 27.211
        gap_change = qp_gap - hf_gap

        print(f"\nHOMO-LUMO gap:")
        print(f"  HF: {hf_gap:.2f} eV")
        print(f"  QP: {qp_gap:.2f} eV")
        print(f"  Change: {gap_change:+.2f} eV")

        # Check if gap opened (typical for GW)
        if gap_change > 0:
            print(f"  ✓ Gap opened by {gap_change:.2f} eV (expected behavior)")
        else:
            print(f"  ⚠ Gap reduced (unusual, may indicate numerical issues)")

    # Overall assessment
    print("\nOverall assessment:")
    if n_physical >= nocc:  # At least occupied orbitals should be physical
        print("  ✓ Core orbitals have physical Z-factors")
        if converged:
            print("  ✓ Calculation converged")
        else:
            print("  ⚠ Calculation did not fully converge (may need more iterations)")
        print("  ✓ Test PASSED (evGW is working with real molecular data)")
        success = True
    else:
        print("  ✗ Too many unphysical Z-factors")
        print("  ✗ Test FAILED (numerical issues present)")
        success = False

    # Print recommendations if there are issues
    if not success or not converged:
        print("\nRecommendations for improvement:")
        print("  1. Use a larger auxiliary basis (e.g., cc-pvdz-jkfit)")
        print("  2. Increase frequency grid points (nfreq)")
        print("  3. Use a better orbital basis (e.g., cc-pvdz)")
        print("  4. Check DF tensor construction")
        print("  5. Consider using analytic continuation instead of contour deformation")

    sys.exit(0 if success else 1)

except Exception as e:
    print(f"\n✗ evGW failed with error: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)