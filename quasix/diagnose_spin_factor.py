#!/usr/bin/env python3
"""
Simple diagnostic to check spin factor in P0 calculation
"""
import numpy as np
from pyscf import gto, scf, df
import quasix

def diagnose_p0_spin_factor():
    """Check if P0 has the correct spin factor"""

    # Simple H2 molecule
    mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g', verbose=0)
    mf = scf.RHF(mol).run(verbose=0)

    n_mo = mf.mo_coeff.shape[1]
    n_occ = int(np.sum(mf.mo_occ > 0))
    n_virt = n_mo - n_occ

    print(f"\nH2 molecule (sto-3g):")
    print(f"  n_mo = {n_mo}, n_occ = {n_occ}, n_virt = {n_virt}")

    # Get DF tensors
    auxbasis = df.make_auxbasis(mol, mp2fit=True)
    with_df = df.DF(mol, auxbasis=auxbasis)
    with_df.build()

    naux = with_df.get_naoaux()
    print(f"  naux = {naux}")

    # Transform to MO basis using PySCF helper
    df_mo = with_df.get_mo_eri(mf.mo_coeff).reshape(n_mo, n_mo, naux)

    # Extract occ-virt block
    df_ia = df_mo[:n_occ, n_occ:, :].reshape(n_occ * n_virt, naux)
    e_occ = mf.mo_energy[:n_occ]
    e_virt = mf.mo_energy[n_occ:]

    print(f"\nDF tensor (ia|P):")
    print(f"  Shape: {df_ia.shape}")
    print(f"  Norm: {np.linalg.norm(df_ia):.6e}")

    # Test frequency
    xi = 1.0  # iω = i * 1.0 Ha

    # Compute P0 manually with spin factor = 2
    p0_manual_spin2 = np.zeros((naux, naux), dtype=complex)
    for i in range(n_occ):
        for a in range(n_virt):
            de = e_virt[a] - e_occ[i]
            # For imaginary frequency, denominator is: 2 * de / (de^2 + xi^2)
            denom = 2.0 * de / (de**2 + xi**2)
            ia_idx = i * n_virt + a
            for p in range(naux):
                for q in range(naux):
                    p0_manual_spin2[p, q] += denom * df_ia[ia_idx, p] * df_ia[ia_idx, q]

    # Compute P0 manually WITHOUT spin factor (factor = 1)
    p0_manual_spin1 = np.zeros((naux, naux), dtype=complex)
    for i in range(n_occ):
        for a in range(n_virt):
            de = e_virt[a] - e_occ[i]
            # WITHOUT factor of 2
            denom = de / (de**2 + xi**2)
            ia_idx = i * n_virt + a
            for p in range(naux):
                for q in range(naux):
                    p0_manual_spin1[p, q] += denom * df_ia[ia_idx, p] * df_ia[ia_idx, q]

    # Compute using Rust
    omega_complex = complex(0.0, xi)
    p0_rust = quasix.compute_polarizability_p0(
        omega=omega_complex,
        df_ia=df_ia,
        e_occ=e_occ,
        e_virt=e_virt,
        eta=1e-4
    )

    print(f"\nP0 at ω = i·{xi} Ha:")
    print(f"\n  Manual (spin factor = 2):")
    print(f"    Diagonal mean: {np.mean(np.diag(p0_manual_spin2.real)):.6e}")
    print(f"    Trace: {np.trace(p0_manual_spin2.real):.6e}")

    print(f"\n  Manual (spin factor = 1):")
    print(f"    Diagonal mean: {np.mean(np.diag(p0_manual_spin1.real)):.6e}")
    print(f"    Trace: {np.trace(p0_manual_spin1.real):.6e}")

    print(f"\n  Rust implementation:")
    print(f"    Diagonal mean: {np.mean(np.diag(p0_rust.real)):.6e}")
    print(f"    Trace: {np.trace(p0_rust.real):.6e}")

    # Check which one matches
    diff_spin2 = np.linalg.norm(p0_rust - p0_manual_spin2) / np.linalg.norm(p0_manual_spin2)
    diff_spin1 = np.linalg.norm(p0_rust - p0_manual_spin1) / np.linalg.norm(p0_manual_spin1)

    print(f"\n  Relative difference:")
    print(f"    Rust vs Manual(spin=2): {diff_spin2:.6e}")
    print(f"    Rust vs Manual(spin=1): {diff_spin1:.6e}")

    if diff_spin2 < 1e-10:
        print(f"\n✓ Rust P0 matches Manual with spin factor = 2 (CORRECT for RHF)")
    elif diff_spin1 < 1e-10:
        print(f"\n⚠️  Rust P0 matches Manual with spin factor = 1 (WRONG for RHF!)")
        print(f"    This would cause 2x error in P0, propagating to W and Σ_c")
    else:
        print(f"\n❌ Rust P0 doesn't match either version!")
        print(f"    Check the implementation")

    # Now check the ratio
    ratio = np.trace(p0_manual_spin2.real) / np.trace(p0_manual_spin1.real)
    print(f"\n  P0(spin=2) / P0(spin=1) = {ratio:.3f}")

if __name__ == "__main__":
    diagnose_p0_spin_factor()