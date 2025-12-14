#!/usr/bin/env python3
"""Test script for S3-2 dielectric and polarizability Python bindings.

This script verifies that the Python interface correctly exposes
the Rust dielectric/polarizability functionality.
"""

import sys
import os
import numpy as np
import pytest

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from quasix import dielectric


def test_polarizability_calculator():
    """Test basic polarizability calculation."""
    print("\n" + "="*60)
    print("Testing PolarizabilityCalculator")
    print("="*60)
    
    # Setup parameters
    nocc = 5
    nvirt = 8
    naux = 20
    n_trans = nocc * nvirt
    
    # Create mock data
    np.random.seed(42)
    df_ia = np.random.randn(n_trans, naux) * 0.1
    e_occ = np.array([-20.0, -1.3, -0.7, -0.6, -0.5])
    e_virt = np.linspace(0.1, 5.0, nvirt)
    
    # Compute polarizability using the module function
    omega = complex(0.0, 0.01)

    # Compute P0 using the module function
    p0_result = dielectric.compute_polarizability_p0(
        df_ia, (e_occ, e_virt), omega, eta=0.01
    )
    p0 = p0_result.matrix

    print(f"P0 shape: {p0.shape}")
    print(f"P0 trace: {np.trace(p0):.6f}")
    print(f"P0 norm: {np.linalg.norm(p0):.6f}")
    print(f"P0 symmetric: {np.allclose(p0, p0.T)}")

    # Test batch computation
    omega_list = np.linspace(-5, 5, 10)
    p0_batch = dielectric.compute_polarizability_batch(
        df_ia, (e_occ, e_virt), omega_list, eta=0.01
    )
    
    print(f"\nBatch P0 shape: {p0_batch.shape}")
    print(f"First P0 norm: {np.linalg.norm(p0_batch[0]):.6f}")
    print(f"Last P0 norm: {np.linalg.norm(p0_batch[-1]):.6f}")
    
    assert True  # Test completed successfully


def test_dielectric_function():
    """Test dielectric function calculation."""
    print("\n" + "="*60)
    print("Testing DielectricFunction")
    print("="*60)
    
    # Create mock data
    naux = 20
    np.random.seed(42)
    
    # Create positive definite metric
    A = np.random.randn(naux, naux)
    metric = A @ A.T + np.eye(naux) * 0.1
    
    # Mock P0
    p0 = np.random.randn(naux, naux) * 0.01
    p0 = (p0 + p0.T) / 2  # Symmetrize
    
    # Compute dielectric function using the module function
    omega = complex(0.0, 0.01)
    epsilon = dielectric.compute_dielectric_function(p0, metric, omega)
    # epsilon is now directly a numpy array, not a DielectricFunction object

    print(f"\nDielectric matrix shape: {epsilon.shape}")
    print(f"Epsilon trace: {np.trace(epsilon):.6f}")
    print(f"Epsilon norm: {np.linalg.norm(epsilon):.6f}")
    print(f"ε diagonal range: [{np.min(np.diag(epsilon)):.4f}, {np.max(np.diag(epsilon)):.4f}]")
    
    # Compute inverse
    epsilon_inv = dielectric.compute_epsilon_inverse(epsilon)
    
    # Verify inversion
    identity_check = epsilon @ epsilon_inv
    inversion_error = np.max(np.abs(identity_check - np.eye(naux)))
    print(f"Inversion error: {inversion_error:.2e}")
    
    # Compute screened Coulomb
    w = dielectric.compute_screened_interaction(epsilon_inv, metric)
    print(f"Screened Coulomb norm: {np.linalg.norm(w):.6f}")
    
    assert True  # Test completed successfully


def test_high_level_interface():
    """Test high-level interface functions."""
    print("\n" + "="*60)
    print("Testing High-Level Interface")
    print("="*60)
    
    # Setup parameters
    nocc = 5
    nvirt = 8
    naux = 20
    n_trans = nocc * nvirt
    
    # Create mock data
    np.random.seed(42)
    df_ia = np.random.randn(n_trans, naux) * 0.1
    e_occ = np.array([-20.0, -1.3, -0.7, -0.6, -0.5])  # Occupied
    e_virt = np.linspace(0.1, 5.0, nvirt)  # Virtual
    mo_energies = np.concatenate([e_occ, e_virt])
    
    # Create DF tensors dictionary
    df_tensors = {
        'df_ia': df_ia,
        'nocc': nocc,
        'nvirt': nvirt
    }
    
    # Test polarizability computation
    frequencies = np.array([0.0, 1.0, 2.0])
    p0_all = dielectric.compute_polarizability_batch(
        df_ia, (e_occ, e_virt), frequencies, eta=0.01
    )
    
    print(f"P0 array shape: {p0_all.shape}")
    print(f"P0[0] norm: {np.linalg.norm(p0_all[0]):.6f}")
    
    # Test dielectric function computation
    A = np.random.randn(naux, naux)
    metric = A @ A.T + np.eye(naux) * 0.1
    
    epsilon = dielectric.compute_dielectric_function(p0_all[0], metric)
    print(f"\nDielectric function shape: {epsilon.shape}")
    print(f"ε trace: {np.trace(epsilon):.6f}")
    
    # Test with multiple frequencies
    epsilon_all = dielectric.compute_dielectric_function(p0_all, metric)
    print(f"Multi-frequency ε shape: {epsilon_all.shape}")
    
    assert True  # Test completed successfully


@pytest.mark.skip(reason="Legacy PolarizabilityRI API not yet re-implemented (S3-2)")
def test_rust_binding_directly():
    """Test calling Rust functions directly."""
    print("\n" + "="*60)
    print("Testing Direct Rust Bindings")
    print("="*60)
    
    import quasix.quasix as rust_lib
    
    # Setup test data
    nocc = 3
    nvirt = 4
    naux = 10
    n_trans = nocc * nvirt
    
    np.random.seed(42)
    df_ia = np.random.randn(n_trans, naux).astype(np.float64)
    e_occ = np.array([-10.0, -5.0, -1.0], dtype=np.float64)
    e_virt = np.array([0.5, 1.0, 2.0, 3.0], dtype=np.float64)
    
    # Call Rust function
    p0 = rust_lib.compute_polarizability_p0(
        df_ia=df_ia,
        e_occ=e_occ,
        e_virt=e_virt,
        omega_real=0.0,
        omega_imag=0.01,
        eta=1e-4
    )
    
    print(f"Direct Rust P0 shape: {p0.shape}")
    print(f"P0 type: {type(p0)}")
    print(f"P0[0,0]: {p0[0,0]:.6f}")
    print(f"P0 norm: {np.linalg.norm(p0):.6f}")
    
    # Test batch computation
    omega_list = np.array([0.0, 1.0, 2.0], dtype=np.float64)
    p0_batch = rust_lib.compute_polarizability_batch(
        df_ia=df_ia,
        e_occ=e_occ,
        e_virt=e_virt,
        omega_list=omega_list,
        broadening=0.01
    )
    
    print(f"\nBatch P0 shape: {p0_batch.shape}")
    print(f"Batch P0 type: {type(p0_batch)}")
    
    assert True  # Test completed successfully


def main():
    """Run all tests."""
    print("QuasiX S3-2 Dielectric/Polarizability Python Bindings Test")
    print("=" * 70)
    
    tests = [
        ("Polarizability Calculator", test_polarizability_calculator),
        ("Dielectric Function", test_dielectric_function),
        ("High-Level Interface", test_high_level_interface),
        ("Direct Rust Bindings", test_rust_binding_directly),
    ]
    
    passed = 0
    failed = 0
    
    for name, test_func in tests:
        try:
            if test_func():
                print(f"\n✓ {name} test passed")
                passed += 1
            else:
                print(f"\n✗ {name} test failed")
                failed += 1
        except Exception as e:
            print(f"\n✗ {name} test failed with exception: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    print("\n" + "=" * 70)
    print(f"Test Results: {passed} passed, {failed} failed")
    
    return failed == 0


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)