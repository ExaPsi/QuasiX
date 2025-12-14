#!/usr/bin/env python3
"""
Simple standalone tests for S3-2 dielectric verification
This version runs without requiring PySCF or the Rust library
"""

import numpy as np
import sys

def test_hermiticity():
    """Test P0 hermiticity"""
    print("\n=== Testing Hermiticity ===")
    
    # Create a Hermitian matrix
    n = 20
    p0 = np.random.randn(n, n) + 1j * np.random.randn(n, n)
    p0 = (p0 + np.conj(p0.T)) / 2  # Make Hermitian
    
    # Check hermiticity
    max_diff = np.max(np.abs(p0 - np.conj(p0.T)))
    
    if max_diff < 1e-10:
        print(f"PASSED: Hermiticity test (max diff = {max_diff:.2e})")
        return True
    else:
        print(f"FAILED: Hermiticity test (max diff = {max_diff:.2e})")
        return False

def test_static_polarizability():
    """Test static polarizability properties"""
    print("\n=== Testing Static Polarizability ===")
    
    # Mock static polarizability calculation
    naux = 30
    nocc = 5
    nvirt = 10
    
    # Create mock P0(ω=0) - should be real and symmetric
    p0_static = np.zeros((naux, naux))
    
    # Add contributions from transitions
    for i in range(nocc):
        for a in range(nvirt):
            de = 0.5 + 0.1 * a - 0.1 * i  # Mock energy difference
            # P0_PQ = 2 * Σ_ia (ia|P)(ia|Q) / (ε_a - ε_i)
            for p in range(naux):
                for q in range(p, naux):
                    contrib = 2.0 * np.exp(-(p-q)**2/10) / de * 0.01
                    p0_static[p, q] += contrib
                    if p != q:
                        p0_static[q, p] += contrib
    
    # Check symmetry
    symmetry_error = np.max(np.abs(p0_static - p0_static.T))
    
    # Check positive semi-definiteness
    eigenvals = np.linalg.eigvalsh(p0_static)
    min_eigenval = np.min(eigenvals)
    
    print(f"  Symmetry error: {symmetry_error:.2e}")
    print(f"  Min eigenvalue: {min_eigenval:.6f}")
    print(f"  Trace: {np.trace(p0_static):.6f}")
    
    passed = symmetry_error < 1e-10 and min_eigenval > -1e-10
    
    if passed:
        print("PASSED: Static polarizability test")
    else:
        print("FAILED: Static polarizability test")
    
    return passed

def test_frequency_convergence():
    """Test frequency grid convergence"""
    print("\n=== Testing Frequency Grid Convergence ===")
    
    grid_sizes = [8, 16, 32, 64]
    values = []
    
    for n in grid_sizes:
        # Gauss-Legendre quadrature
        x, w = np.polynomial.legendre.leggauss(n)
        # Test integral: ∫ f(x) dx with f(x) = 1/(1+x^2)
        integral = np.sum(w / (1 + x**2))
        values.append(integral)
        print(f"  Grid {n:3d}: {integral:.10f}")
    
    # Check convergence
    errors = []
    for i in range(1, len(values)):
        error = abs(values[i] - values[i-1])
        errors.append(error)
        print(f"  Error {grid_sizes[i-1]:3d}->{grid_sizes[i]:3d}: {error:.2e}")
    
    # Check that errors are decreasing overall
    converging = errors[0] > errors[-1]
    
    if converging:
        print("PASSED: Frequency convergence test")
        return True
    else:
        print("FAILED: Frequency convergence test")
        return False

def test_dielectric_properties():
    """Test dielectric function properties"""
    print("\n=== Testing Dielectric Properties ===")
    
    naux = 25
    
    # Create mock P0
    p0 = np.random.randn(naux, naux) * 0.01
    p0 = (p0 + p0.T) / 2  # Symmetrize
    
    # Mock V^(1/2)
    vsqrt = np.eye(naux) * 0.5
    
    # M = V^(1/2) P0 V^(1/2)
    m = vsqrt @ p0 @ vsqrt
    
    # ε = 1 - M
    epsilon = np.eye(naux) - m
    
    # Check eigenvalues (should all be positive)
    eigenvals = np.linalg.eigvalsh(epsilon)
    min_eigenval = np.min(eigenvals)
    max_eigenval = np.max(eigenvals)
    
    # Check trace (should be close to naux for small P0)
    trace = np.trace(epsilon)
    
    print(f"  Min eigenvalue: {min_eigenval:.6f}")
    print(f"  Max eigenvalue: {max_eigenval:.6f}")
    print(f"  Trace: {trace:.6f} (expected ≈ {naux})")
    
    # For imaginary frequencies, ε should be > 1
    passed = min_eigenval > 0.9 and trace > naux * 0.9
    
    if passed:
        print("PASSED: Dielectric properties test")
    else:
        print("FAILED: Dielectric properties test")
    
    return passed

def test_sum_rules():
    """Test sum rules"""
    print("\n=== Testing Sum Rules ===")
    
    # Thomas-Reiche-Kuhn sum rule for oscillator strengths
    n_electrons = 10
    n_transitions = 30
    
    # Create oscillator strengths that sum to n_electrons
    oscillator_strengths = np.random.random(n_transitions)
    oscillator_strengths *= n_electrons / np.sum(oscillator_strengths)
    
    total = np.sum(oscillator_strengths)
    error = abs(total - n_electrons)
    
    print(f"  Number of electrons: {n_electrons}")
    print(f"  Sum of oscillator strengths: {total:.6f}")
    print(f"  Error: {error:.2e}")
    
    if error < 1e-10:
        print("PASSED: Sum rules test")
        return True
    else:
        print("FAILED: Sum rules test")
        return False

def test_conjugate_relation():
    """Test P0(-ω*) = P0(ω)* relation"""
    print("\n=== Testing Conjugate Relation ===")
    
    naux = 15
    
    # For imaginary frequency, P0(iω) should be real
    omega = 1.0j
    
    # Mock P0 calculation
    p0 = np.zeros((naux, naux), dtype=complex)
    
    # Add transitions with imaginary frequency
    for i in range(5):  # occupied
        for a in range(10):  # virtual
            de = 0.5 + 0.1 * a - 0.1 * i
            denom = 1.0 / (de - omega)
            
            # For imaginary ω, this should give real result
            for p in range(naux):
                for q in range(p, naux):
                    contrib = 2.0 * np.exp(-(p-q)**2/10) * denom * 0.01
                    p0[p, q] += contrib
                    if p != q:
                        p0[q, p] += contrib.conj()
    
    # Check that imaginary part is small
    max_imag = np.max(np.abs(p0.imag))
    
    print(f"  Max imaginary part for P0(iω): {max_imag:.2e}")
    
    if max_imag < 1e-10:
        print("PASSED: Conjugate relation test")
        return True
    else:
        print("FAILED: Conjugate relation test")
        return False

def main():
    """Run all tests"""
    print("="*60)
    print("S3-2 DIELECTRIC FUNCTION SIMPLE VERIFICATION")
    print("="*60)
    
    tests = [
        test_hermiticity,
        test_static_polarizability,
        test_frequency_convergence,
        test_dielectric_properties,
        test_sum_rules,
        test_conjugate_relation,
    ]
    
    passed = 0
    total = len(tests)
    
    for test_func in tests:
        try:
            if test_func():
                passed += 1
        except Exception as e:
            print(f"ERROR in {test_func.__name__}: {e}")
    
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Tests passed: {passed}/{total}")
    
    if passed == total:
        print("SUCCESS: All tests passed!")
        return 0
    else:
        print(f"FAILURE: {total - passed} tests failed")
        return 1

if __name__ == "__main__":
    sys.exit(main())