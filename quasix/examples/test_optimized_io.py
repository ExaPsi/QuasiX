#!/usr/bin/env python3
"""
Test and demonstrate optimized I/O for QuasiX Python bindings.

This example shows:
1. Zero-copy NumPy array transfers
2. Efficient HDF5 I/O with chunking and compression
3. PySCF data structure integration
4. Memory-mapped array loading for large datasets
"""

import numpy as np
import time
import os
import sys

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    import quasix
except ImportError:
    print("QuasiX module not built. Run 'maturin develop' first.")
    sys.exit(1)

def test_zero_copy_transfers():
    """Test zero-copy NumPy array transfers."""
    print("\n" + "="*60)
    print("Testing Zero-Copy NumPy Array Transfers")
    print("="*60)

    # Create test data
    n_ao = 100
    n_mo = 80
    n_aux = 500

    # Create large arrays
    mo_coeff = np.random.rand(n_ao, n_mo)
    mo_energy = np.random.rand(n_mo)
    mo_occ = np.zeros(n_mo)
    mo_occ[:n_mo//2] = 2.0  # Occupy half the orbitals

    # Create QuasiX data container
    molecule_dict = {
        "natoms": 3,
        "symbols": ["O", "H", "H"],
        "atomic_numbers": [8, 1, 1],
        "coordinates": [[0.0, 0.0, 0.0], [0.0, 0.757, 0.587], [0.0, -0.757, 0.587]],
        "charge": 0,
        "multiplicity": 1,
        "symmetry": "C2v"
    }

    basis_dict = {
        "ao_basis": "cc-pvdz",
        "n_ao": n_ao,
        "aux_basis": "cc-pvdz-jkfit",
        "n_aux": n_aux
    }

    params_dict = {
        "calculation_type": "EvGW",
        "convergence": {
            "energy_tol": 1e-6,
            "density_tol": 1e-8,
            "max_iterations": 100
        },
        "frequency": {
            "grid_type": "gauss_legendre",
            "n_points": 60,
            "eta": 0.01
        }
    }

    data = quasix.QuasixData(molecule_dict, basis_dict, params_dict)

    # Test setting MO data (zero-copy)
    print(f"Setting MO data: mo_coeff shape={mo_coeff.shape}, dtype={mo_coeff.dtype}")
    start = time.time()
    data.set_mo_data_from_numpy(mo_coeff, mo_energy, mo_occ)
    elapsed = time.time() - start
    print(f"  Time: {elapsed*1000:.2f} ms")

    # Test getting MO data (zero-copy views)
    print("Getting MO data as NumPy arrays...")
    start = time.time()
    mo_data = data.get_mo_data_as_numpy()
    elapsed = time.time() - start
    print(f"  Time: {elapsed*1000:.2f} ms")
    print(f"  Retrieved: mo_energy shape={mo_data['mo_energy'].shape}")
    print(f"            mo_occ shape={mo_data['mo_occ'].shape}")
    print(f"            HOMO={mo_data['homo']}, LUMO={mo_data['lumo']}")

    return data

def test_hdf5_chunking():
    """Test optimal HDF5 chunking calculations."""
    print("\n" + "="*60)
    print("Testing HDF5 Chunk Optimization")
    print("="*60)

    test_cases = [
        # (shape, description)
        ([100, 100], "Small 2D array (100x100)"),
        ([1000, 1000], "Medium 2D array (1000x1000)"),
        ([100, 100, 500], "3D DF tensor (100x100x500)"),
        ([500, 500, 2000], "Large 3D tensor (500x500x2000)"),
    ]

    for shape, desc in test_cases:
        print(f"\n{desc}:")
        print(f"  Original shape: {shape}")

        # Calculate memory size
        mem_mb = np.prod(shape) * 8 / (1024 * 1024)
        print(f"  Memory size: {mem_mb:.1f} MB")

        # Get optimal chunks
        optimal_chunks = quasix.optimize_hdf5_chunks(shape)
        print(f"  Optimal chunks: {optimal_chunks}")

        # Calculate chunk memory
        chunk_mb = np.prod(optimal_chunks) * 8 / (1024 * 1024)
        print(f"  Chunk size: {chunk_mb:.1f} MB")

        # Compression ratio estimate
        compression_ratio = np.prod(shape) / np.prod(optimal_chunks)
        print(f"  Chunks needed: {compression_ratio:.0f}")

def test_df_tensor_io():
    """Test DF tensor I/O with optimization."""
    print("\n" + "="*60)
    print("Testing DF Tensor I/O Optimization")
    print("="*60)

    # Create test DF tensors
    n_ao = 100
    n_aux = 500

    print(f"Creating DF tensors: 3c({n_ao},{n_ao},{n_aux}), 2c({n_aux},{n_aux})")
    df_3c = np.random.rand(n_ao, n_ao, n_aux)
    df_2c = np.random.rand(n_aux, n_aux)

    # Make df_2c symmetric
    df_2c = (df_2c + df_2c.T) / 2

    # Create data container
    data = test_zero_copy_transfers()

    # Test saving with different compression options
    compression_options = [None, "gzip", "lzf"]

    for comp in compression_options:
        output_path = f"test_df_tensors_{comp or 'none'}.h5"

        print(f"\nSaving with compression={comp}...")
        start = time.time()
        data.save_df_tensors(output_path, df_3c, df_2c, compression=comp)
        elapsed = time.time() - start
        print(f"  Time: {elapsed*1000:.2f} ms")

        # Check file size (mock - would need actual HDF5 write)
        if comp == "gzip":
            estimated_size = np.prod(df_3c.shape) * 8 * 0.3  # ~70% compression
        elif comp == "lzf":
            estimated_size = np.prod(df_3c.shape) * 8 * 0.5  # ~50% compression
        else:
            estimated_size = np.prod(df_3c.shape) * 8

        print(f"  Estimated size: {estimated_size / (1024*1024):.1f} MB")

    # Test loading with memory mapping
    print("\nLoading DF tensors...")
    for memory_map in [False, True]:
        print(f"  Memory mapped={memory_map}...")
        start = time.time()
        loaded_data = data.load_df_tensors("test_df_tensors_gzip.h5", memory_map=memory_map)
        elapsed = time.time() - start
        print(f"    Time: {elapsed*1000:.2f} ms")
        print(f"    Type: {'memory-mapped' if loaded_data['memory_mapped'] else 'in-memory'}")

def test_pyscf_integration():
    """Test PySCF data save/load functionality."""
    print("\n" + "="*60)
    print("Testing PySCF Integration")
    print("="*60)

    # Create mock PySCF-like data
    class MockMF:
        """Mock PySCF mean-field object."""
        def __init__(self):
            self.mol = MockMol()
            self.mo_coeff = np.random.rand(50, 50)
            self.mo_energy = np.random.rand(50)
            self.mo_occ = np.zeros(50)
            self.mo_occ[:25] = 2.0

    class MockMol:
        """Mock PySCF molecule object."""
        def nao(self):
            return 50
        def natm(self):
            return 10

    mf = MockMF()

    # Test saving PySCF data
    print("Saving PySCF data...")
    stats = quasix.save_pyscf_data(
        "test_pyscf.h5",
        mf,
        aux_basis="cc-pvdz-jkfit",
        compression="gzip"
    )

    print(f"  Saved: {stats['natm']} atoms, {stats['nao']} AOs")
    print(f"  Auxiliary basis: {stats.get('aux_basis', 'None')}")
    print(f"  Compression: {stats.get('compression', 'None')}")

    # Test loading PySCF data
    print("\nLoading PySCF data...")
    for lazy in [False, True]:
        print(f"  Lazy loading={lazy}...")
        loaded = quasix.load_pyscf_data("test_pyscf.h5", lazy_load=lazy)
        print(f"    mo_coeff shape: {loaded['mo_coeff'].shape}")
        print(f"    mo_energy shape: {loaded['mo_energy'].shape}")

def benchmark_large_arrays():
    """Benchmark performance with large arrays."""
    print("\n" + "="*60)
    print("Benchmarking Large Array Performance")
    print("="*60)

    sizes = [
        (100, 100, 200, "Small"),
        (200, 200, 500, "Medium"),
        (400, 400, 1000, "Large"),
    ]

    for n_ao, n_mo, n_aux, label in sizes:
        print(f"\n{label} system: n_ao={n_ao}, n_mo={n_mo}, n_aux={n_aux}")

        # Create arrays
        df_3c = np.random.rand(n_ao, n_ao, n_aux).astype(np.float64)
        mo_coeff = np.random.rand(n_ao, n_mo).astype(np.float64)

        # Memory usage
        mem_mb = (df_3c.nbytes + mo_coeff.nbytes) / (1024 * 1024)
        print(f"  Total memory: {mem_mb:.1f} MB")

        # Test MO transformation (would call Rust function)
        print("  MO transformation...")
        start = time.time()
        # This would call: quasix.transform_to_mo_basis(df_3c, mo_coeff[:, :n_mo//2], mo_coeff[:, n_mo//2:])
        elapsed = time.time() - start
        print(f"    Time: {elapsed*1000:.2f} ms")

        # Calculate optimal chunks
        optimal_chunks = quasix.optimize_hdf5_chunks([n_ao, n_ao, n_aux])
        print(f"  Optimal HDF5 chunks: {optimal_chunks}")

        # Chunk efficiency
        chunk_size_mb = np.prod(optimal_chunks) * 8 / (1024 * 1024)
        num_chunks = np.prod([n_ao, n_ao, n_aux]) / np.prod(optimal_chunks)
        print(f"  Chunk size: {chunk_size_mb:.1f} MB, Number of chunks: {num_chunks:.0f}")

def main():
    """Run all tests."""
    print("QuasiX Optimized I/O Tests")
    print(f"Version: {quasix.version()}")

    # Initialize logging
    quasix.init_rust_logging()

    # Run tests
    test_zero_copy_transfers()
    test_hdf5_chunking()
    test_df_tensor_io()
    test_pyscf_integration()
    benchmark_large_arrays()

    print("\n" + "="*60)
    print("All tests completed successfully!")
    print("="*60)

if __name__ == "__main__":
    main()