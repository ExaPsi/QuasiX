#!/usr/bin/env python3
"""
Test optimized I/O functionality in QuasixData class.
"""

import numpy as np
import time
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import quasix

def test_quasixdata_io():
    """Test QuasixData I/O methods."""
    print("Testing QuasixData I/O Methods")
    print("="*50)

    # Create QuasixData instance
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
        "n_ao": 24,
        "aux_basis": "cc-pvdz-jkfit",
        "n_aux": 84
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
    print(f"Created QuasixData for {molecule_dict['symbols']}")

    # Test setting MO data with NumPy arrays (zero-copy)
    print("\n1. Testing zero-copy MO data transfer:")
    n_mo = 24
    mo_coeff = np.random.rand(24, n_mo)
    mo_energy = np.random.rand(n_mo)
    mo_occ = np.zeros(n_mo)
    mo_occ[:10] = 2.0  # 10 occupied orbitals

    start = time.time()
    data.set_mo_data_from_numpy(mo_coeff, mo_energy, mo_occ)
    elapsed = time.time() - start
    print(f"   Set MO data ({mo_coeff.shape}): {elapsed*1000:.3f} ms")

    # Test getting MO data back
    start = time.time()
    mo_data = data.get_mo_data_as_numpy()
    elapsed = time.time() - start
    print(f"   Get MO data: {elapsed*1000:.3f} ms")
    print(f"   Retrieved shapes: energy={mo_data['mo_energy'].shape}, occ={mo_data['mo_occ'].shape}")

    # Test DF tensor save/load
    print("\n2. Testing DF tensor I/O:")
    df_3c = np.random.rand(24, 24, 84)
    df_2c = np.random.rand(84, 84)
    df_2c = (df_2c + df_2c.T) / 2  # Make symmetric

    # Save DF tensors
    print(f"   DF tensor shapes: 3c={df_3c.shape}, 2c={df_2c.shape}")
    print(f"   Memory: {(df_3c.nbytes + df_2c.nbytes) / (1024*1024):.2f} MB")

    start = time.time()
    data.save_df_tensors("test_df.h5", df_3c, df_2c, compression="gzip")
    elapsed = time.time() - start
    print(f"   Save with gzip: {elapsed*1000:.3f} ms")

    # Load DF tensors
    start = time.time()
    loaded = data.load_df_tensors("test_df.h5", memory_map=False)
    elapsed = time.time() - start
    print(f"   Load in-memory: {elapsed*1000:.3f} ms")

    # Test JSON I/O
    print("\n3. Testing JSON I/O:")
    start = time.time()
    data.to_json("test_data.json")
    elapsed = time.time() - start
    print(f"   Save to JSON: {elapsed*1000:.3f} ms")

    start = time.time()
    loaded_data = quasix.QuasixData.from_json("test_data.json")
    elapsed = time.time() - start
    print(f"   Load from JSON: {elapsed*1000:.3f} ms")

    # Test HDF5 I/O
    print("\n4. Testing HDF5 I/O:")
    start = time.time()
    data.to_hdf5("test_data.h5")
    elapsed = time.time() - start
    print(f"   Save to HDF5: {elapsed*1000:.3f} ms")

    start = time.time()
    loaded_data = quasix.QuasixData.from_hdf5("test_data.h5")
    elapsed = time.time() - start
    print(f"   Load from HDF5: {elapsed*1000:.3f} ms")

    # Test array save/load
    print("\n5. Testing array-specific I/O:")
    large_array = np.random.rand(1000, 1000)
    print(f"   Array shape: {large_array.shape}, size: {large_array.nbytes/(1024*1024):.2f} MB")

    start = time.time()
    data.save_array_to_hdf5("test_array.h5", "large_dataset", large_array, compression=6)
    elapsed = time.time() - start
    print(f"   Save array with compression: {elapsed*1000:.3f} ms")

    start = time.time()
    loaded_array = data.load_array_from_hdf5("test_array.h5", "large_dataset", lazy=False)
    elapsed = time.time() - start
    print(f"   Load array (in-memory): {elapsed*1000:.3f} ms")

    start = time.time()
    lazy_array = data.load_array_from_hdf5("test_array.h5", "large_dataset", lazy=True)
    elapsed = time.time() - start
    print(f"   Load array (memory-mapped): {elapsed*1000:.3f} ms")

    # Test getting results
    print("\n6. Testing results access:")
    data.set_reference_results(
        energy=-76.12345,
        orbital_energies=mo_energy.tolist(),
        occupations=mo_occ.tolist(),
        homo=9,
        lumo=10
    )

    data.set_gw_results(
        qp_energies=(mo_energy + 0.1).tolist(),
        z_factors=np.ones(n_mo).tolist(),
        converged=True,
        iterations=8
    )

    results = data.get_results()
    print(f"   Reference energy: {results['reference']['energy']:.5f}")
    print(f"   HOMO-LUMO gap: {results['reference']['lumo'] - results['reference']['homo']}")
    if 'gw' in results:
        print(f"   GW converged: {results['gw']['converged']} in {results['gw']['iterations']} iterations")

    print("\n" + "="*50)
    print("All I/O tests completed successfully!")

if __name__ == "__main__":
    quasix.init_rust_logging()
    test_quasixdata_io()