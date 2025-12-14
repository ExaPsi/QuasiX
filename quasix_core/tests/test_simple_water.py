#!/usr/bin/env python3
"""Simple test of evGW with water molecule."""

import sys
import numpy as np
from pathlib import Path
import time
import os

# Disable verbose logging
os.environ['RUST_LOG'] = 'error'

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


# Build water molecule
print("Building H2O molecule...")
mol = gto.M(
    atom='''
    O    0.000000    0.000000    0.117790
    H    0.000000    0.756950   -0.471161
    H    0.000000   -0.756950   -0.471161
    ''',
    basis='sto-3g',  # Use minimal basis for speed
    symmetry=False,
    charge=0,
    spin=0
)
mol.build()
print(f"  Atoms: {mol.natm}, Electrons: {mol.nelectron}, AOs: {mol.nao}")

# Run HF
print("Running HF...")
mf = scf.RHF(mol)
mf.verbose = 0
mf.kernel()
print(f"  HF energy: {mf.e_tot:.6f} Ha")
print(f"  HOMO: {mf.mo_energy[mol.nelec[0]-1]:.4f} Ha")
print(f"  LUMO: {mf.mo_energy[mol.nelec[0]]:.4f} Ha")

# Build DF tensors
print("Building DF tensors...")
nao = mol.nao
nmo = mf.mo_coeff.shape[1]
nocc = mol.nelec[0]
nvir = nmo - nocc

# For STO-3G, use a minimal auxiliary basis
auxbasis = 'def2-svp-jkfit'
auxmol = df.addons.make_auxmol(mol, auxbasis=auxbasis)
naux = auxmol.nao
print(f"  NAO: {nao}, NMO: {nmo}, Nocc: {nocc}, Nvir: {nvir}, Naux: {naux}")

# Build V matrix
v_aux = auxmol.intor('int2c2e')
v_aux = (v_aux + v_aux.T) / 2

# Add regularization if needed
eigvals = np.linalg.eigvalsh(v_aux)
if eigvals.min() < 1e-10:
    v_aux += (1e-10 - eigvals.min()) * np.eye(naux)
print(f"  V_aux condition: min_eig={eigvals.min():.2e}, max_eig={eigvals.max():.2e}")

# Build 3-center integrals
int3c = df.incore.aux_e2(mol, auxmol, intor='int3c2e')
eri_3c = int3c.reshape(naux, nao, nao)

# Transform to MO
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

print(f"  iaP shape: {iaP.shape}, norm: {np.linalg.norm(iaP):.2f}")
print(f"  ijP shape: {ijP.shape}, norm: {np.linalg.norm(ijP):.2f}")

# Prepare data
mo_energy = mf.mo_energy.astype(np.float64)
mo_occ = np.zeros(len(mo_energy))
mo_occ[:nocc] = 2.0
vxc_dft = np.zeros(nmo)

# Run evGW with minimal parameters
print("\nRunning evGW...")
print("  max_cycle: 5")
print("  conv_tol: 1e-3")
print("  damping: 0.5")

try:
    start = time.time()
    result = rust_module.evgw(
        mo_energy,
        mo_occ,
        iaP,
        ijP,
        v_aux,
        vxc_dft,
        max_cycle=5,      # Fewer iterations
        conv_tol=1e-3,    # Looser convergence
        conv_tol_z=1e-2,  # Looser Z convergence
        damping=0.5,      # More damping
        damping_dynamic=False,
        diis=False,
        diis_space=6,
        nfreq=24,         # Fewer frequency points
        eta=0.01,
        verbose=0
    )
    elapsed = time.time() - start

    print(f"\n✓ evGW completed in {elapsed:.2f} seconds")

    # Extract results
    qp_energies = np.array(result['qp_energies'])
    z_factors = np.array(result['z_factors'])
    converged = result['converged']
    n_cycles = result['n_cycles']

    print(f"  Converged: {converged}")
    print(f"  Iterations: {n_cycles}")

    # Check Z-factors
    z_min, z_max = z_factors.min(), z_factors.max()
    print(f"  Z-factors: [{z_min:.3f}, {z_max:.3f}]")

    if z_min < 0 or z_max > 1:
        print(f"  ⚠ Warning: Unphysical Z-factors")
    else:
        print(f"  ✓ Z-factors are physical")

    # QP corrections
    qp_corrections = (qp_energies - mo_energy) * 27.211  # to eV

    print(f"\nQP corrections (eV):")
    print(f"  HOMO (MO {nocc-1}): {qp_corrections[nocc-1]:.3f}")
    print(f"  LUMO (MO {nocc}): {qp_corrections[nocc]:.3f}")

    # Gaps
    hf_gap = (mo_energy[nocc] - mo_energy[nocc-1]) * 27.211
    qp_gap = (qp_energies[nocc] - qp_energies[nocc-1]) * 27.211

    print(f"\nHOMO-LUMO gap:")
    print(f"  HF: {hf_gap:.2f} eV")
    print(f"  QP: {qp_gap:.2f} eV")
    print(f"  Change: {qp_gap - hf_gap:.2f} eV")

    print("\n✓ Test passed!")

except Exception as e:
    print(f"\n✗ evGW failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)