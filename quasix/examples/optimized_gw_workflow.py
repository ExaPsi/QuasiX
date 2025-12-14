#!/usr/bin/env python3
"""
Example: Optimized GW workflow using S2-5 PySCF adapter optimizations

This example demonstrates how to use the optimized PySCF adapter
for efficient DF tensor construction in GW calculations.
"""

import numpy as np
import time
from pyscf import gto, scf, df
import sys
import os

# Add quasix to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from quasix.pyscf_adapter import PySCFAdapter, MoleculeFactory
from quasix.df_tensor import DFTensor
from quasix.mo_transform_optimized import (
    build_transition_df_tensor_optimized,
    build_full_df_tensors_optimized
)


def optimized_gw_workflow(mol_name='benzene', basis='def2-svp'):
    """
    Demonstrate optimized GW workflow with performance monitoring
    """

    print(f"\n{'='*60}")
    print(f"Optimized GW Workflow: {mol_name}/{basis}")
    print(f"{'='*60}")

    # Step 1: Create molecule and run SCF
    print("\n1. Setting up molecule and running SCF...")
    factory = MoleculeFactory()

    if mol_name == 'water':
        mol = factory.water(basis)
    elif mol_name == 'benzene':
        mol = factory.benzene(basis)
    elif mol_name == 'ammonia':
        mol = factory.ammonia(basis)
    else:
        raise ValueError(f"Unknown molecule: {mol_name}")

    print(f"   Molecule: {mol.natm} atoms, {mol.nao} AOs, {mol.nelectron} electrons")

    t0 = time.time()
    mf = scf.RHF(mol)
    mf.kernel()
    scf_time = time.time() - t0

    print(f"   SCF converged: {mf.converged}")
    print(f"   SCF energy: {mf.e_tot:.8f} Ha")
    print(f"   SCF time: {scf_time:.2f} s")

    # Step 2: Extract data with PySCF adapter
    print("\n2. Extracting data with PySCF adapter...")
    adapter = PySCFAdapter(tolerance=1e-8)

    # Choose auxiliary basis
    if 'def2' in basis:
        auxbasis = basis.replace('def2-', 'def2-') + '-jkfit'
    else:
        auxbasis = 'def2-svp-jkfit'

    t0 = time.time()
    data = adapter.extract_from_mf(mf, auxbasis=auxbasis, build_df=False)  # Don't build DF yet
    extract_time = time.time() - t0

    print(f"   Dimensions: nocc={data.nocc}, nvir={data.nvir}, naux={data.naux}")
    print(f"   Extraction time: {extract_time:.3f} s")

    # Step 3: Build DF tensors with optimized method
    print("\n3. Building optimized DF tensors...")

    auxmol = df.addons.make_auxmol(mol, auxbasis)

    # Use optimized construction
    t0 = time.time()
    df_tensors = build_full_df_tensors_optimized(
        mol, auxmol, mf.mo_coeff, mf.mo_occ,
        max_memory=4000,  # 4 GB limit
        use_parallel=False  # Can enable if OpenMP available
    )
    df_time = time.time() - t0

    print(f"   DF tensor construction time: {df_time:.3f} s")
    print(f"   Tensor shapes:")
    print(f"     (ia|P): {df_tensors['iaP'].shape}")
    print(f"     (ij|P): {df_tensors['ijP'].shape}")
    print(f"     (ab|P): {df_tensors['abP'].shape}")

    # Calculate memory usage
    total_memory = (df_tensors['iaP'].nbytes + df_tensors['ijP'].nbytes +
                   df_tensors['abP'].nbytes) / 1e9
    print(f"   Total DF tensor memory: {total_memory:.2f} GB")

    # Step 4: Prepare for Rust transfer
    print("\n4. Preparing data for Rust backend...")

    # Attach optimized DF tensors to data
    data.df_tensor = DFTensor(mol, auxbasis)
    data.df_tensor.iaP = df_tensors['iaP']
    data.df_tensor.ijP = df_tensors['ijP']
    data.df_tensor.abP = df_tensors['abP']
    data.df_tensor.metric = df_tensors['metric']
    data.df_tensor.nocc = df_tensors['nocc']
    data.df_tensor.nvir = df_tensors['nvir']

    t0 = time.time()
    rust_data = adapter.prepare_for_rust(data)
    prep_time = time.time() - t0

    print(f"   Rust preparation time: {prep_time:.3f} s")
    print(f"   Number of arrays prepared: {len(rust_data)}")

    # Verify all arrays are C-contiguous
    for key, arr in rust_data.items():
        if isinstance(arr, np.ndarray):
            if not arr.flags['C_CONTIGUOUS']:
                print(f"   WARNING: {key} is not C-contiguous!")

    # Step 5: Simulate GW calculation (would call Rust here)
    print("\n5. Ready for GW calculation...")
    print("   In production, would call: quasix.evgw_diag(...)")

    # Example of how data would be passed to Rust:
    """
    import quasix

    gw_result = quasix.evgw_diag(
        mo_energy=rust_data['mo_energy'],
        mo_coeff=rust_data['mo_coeff'],
        df_ia=rust_data['df_ia'],
        df_ij=rust_data['df_ij'],
        df_ab=rust_data['df_ab'],
        metric=rust_data['metric'],
        n_occ=data.nocc,
        n_vir=data.nvir,
        n_aux=data.naux,
        omega_grid=np.linspace(0, 50, 100),
        eta=0.01,
        convergence_tol=1e-6,
        max_iter=100
    )

    qp_energies = gw_result['qp_energies']
    z_factors = gw_result['z_factors']
    """

    # Step 6: Performance summary
    print(f"\n{'='*60}")
    print("PERFORMANCE SUMMARY")
    print(f"{'='*60}")
    print(f"SCF calculation:     {scf_time:.3f} s")
    print(f"Data extraction:     {extract_time:.3f} s")
    print(f"DF tensor build:     {df_time:.3f} s")
    print(f"Rust preparation:    {prep_time:.3f} s")
    print(f"{'Total setup time:':<20} {scf_time + extract_time + df_time + prep_time:.3f} s")
    print(f"\nMemory usage:        {total_memory:.2f} GB")

    # Return data for further processing
    return {
        'mf': mf,
        'data': data,
        'df_tensors': df_tensors,
        'rust_data': rust_data,
        'timings': {
            'scf': scf_time,
            'extract': extract_time,
            'df_build': df_time,
            'rust_prep': prep_time
        }
    }


def compare_original_vs_optimized(mol_name='water', basis='def2-svp'):
    """
    Compare original vs optimized DF tensor construction
    """

    print(f"\n{'='*60}")
    print(f"Comparing Original vs Optimized: {mol_name}/{basis}")
    print(f"{'='*60}")

    # Setup molecule
    factory = MoleculeFactory()
    if mol_name == 'water':
        mol = factory.water(basis)
    elif mol_name == 'benzene':
        mol = factory.benzene(basis)
    else:
        mol = factory.ammonia(basis)

    # Run SCF
    mf = scf.RHF(mol)
    mf.kernel()

    if not mf.converged:
        print("SCF did not converge!")
        return

    # Setup auxiliary basis
    if 'def2' in basis:
        auxbasis = basis.replace('def2-', 'def2-') + '-jkfit'
    else:
        auxbasis = 'def2-svp-jkfit'

    auxmol = df.addons.make_auxmol(mol, auxbasis)

    print(f"\nSystem size: nao={mol.nao}, naux={auxmol.nao}")
    print(f"MO dimensions: nocc={np.sum(mf.mo_occ > 0)}, nvir={mol.nao - np.sum(mf.mo_occ > 0)}")

    # Import both versions
    from quasix.mo_transform import build_full_df_tensors

    # Time original version
    print("\nTiming original implementation...")
    t0 = time.time()
    df_orig = build_full_df_tensors(mol, auxmol, mf.mo_coeff, mf.mo_occ)
    time_orig = time.time() - t0
    print(f"  Original time: {time_orig:.3f} s")

    # Time optimized version
    print("\nTiming optimized implementation...")
    t0 = time.time()
    df_opt = build_full_df_tensors_optimized(mol, auxmol, mf.mo_coeff, mf.mo_occ)
    time_opt = time.time() - t0
    print(f"  Optimized time: {time_opt:.3f} s")

    # Calculate speedup
    speedup = time_orig / time_opt if time_opt > 0 else float('inf')
    print(f"\nSpeedup: {speedup:.2f}x")

    # Verify correctness
    error_ia = np.max(np.abs(df_orig['iaP'] - df_opt['iaP']))
    error_ij = np.max(np.abs(df_orig['ijP'] - df_opt['ijP']))
    error_ab = np.max(np.abs(df_orig['abP'] - df_opt['abP']))

    print(f"\nCorrectness check:")
    print(f"  Max error (ia|P): {error_ia:.2e}")
    print(f"  Max error (ij|P): {error_ij:.2e}")
    print(f"  Max error (ab|P): {error_ab:.2e}")

    if max(error_ia, error_ij, error_ab) < 1e-10:
        print("  ✓ Results match within numerical precision")
    else:
        print("  ❌ Results differ significantly!")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Optimized GW workflow example")
    parser.add_argument('--molecule', default='water',
                       choices=['water', 'benzene', 'ammonia'],
                       help='Molecule to use')
    parser.add_argument('--basis', default='def2-svp',
                       help='Basis set')
    parser.add_argument('--compare', action='store_true',
                       help='Compare original vs optimized')

    args = parser.parse_args()

    if args.compare:
        compare_original_vs_optimized(args.molecule, args.basis)
    else:
        optimized_gw_workflow(args.molecule, args.basis)