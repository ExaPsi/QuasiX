#!/usr/bin/env python3
"""Test extraction of W(iω) and Σᶜ from QuasiX G₀W₀ calculation.

This script verifies that the new return values (W_iw, P0_iw, sigma_c)
are accessible from the run_g0w0() function.

Expected output:
- W_iw shape: (nfreq, n_aux, n_aux) - e.g., (100, 28, 28) for H₂/STO-3G
- sigma_c shape: (n_mo,) - e.g., (2,) for H₂
- Both values should be finite and non-zero
"""

import sys
import os
import numpy as np

# Add quasix module to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

def test_w_sigma_c_extraction():
    """Test that W(iω) and Σᶜ are returned from run_g0w0()."""
    # Import PySCF
    try:
        from pyscf import gto, scf
    except ImportError:
        print("ERROR: PySCF not installed. Run: pip install pyscf")
        return False

    # Import QuasiX
    try:
        from quasix.g0w0 import run_g0w0
    except ImportError as e:
        print(f"ERROR: Failed to import QuasiX: {e}")
        print("Make sure QuasiX is installed and built.")
        return False

    print("=" * 70)
    print(" Testing W(iω) and Σᶜ Extraction from QuasiX")
    print("=" * 70)

    # Create simple H₂ molecule
    print("\n1. Creating H₂ molecule (STO-3G basis)...")
    mol = gto.M(
        atom='H 0 0 0; H 0 0 0.74',
        basis='sto-3g',
        verbose=0
    )

    # Run HF
    print("2. Running HF calculation...")
    mf = scf.RHF(mol).run(verbose=0)
    print(f"   HF energy: {mf.e_tot:.6f} Ha")

    # Run G₀W₀ with moderate frequency grid
    print("3. Running G₀W₀ calculation (nfreq=32)...")
    try:
        result = run_g0w0(mol, mf, nfreq=32, verbose=0)
    except Exception as e:
        print(f"ERROR: G₀W₀ calculation failed: {e}")
        return False

    print("4. Checking return values...")

    # Check W_iw
    if 'W_iw' not in result:
        print("   FAIL: 'W_iw' not in result dictionary!")
        return False

    W_iw = result['W_iw']
    print(f"   W_iw shape: {W_iw.shape}")
    print(f"   W_iw dtype: {W_iw.dtype}")
    print(f"   W_iw[0, 0, 0] = {W_iw[0, 0, 0]:.6e} Ha")
    print(f"   W_iw min/max: {W_iw.min():.3e} / {W_iw.max():.3e} Ha")

    # Validate W_iw
    if len(W_iw.shape) != 3:
        print(f"   FAIL: W_iw should be 3D, got shape {W_iw.shape}")
        return False

    if not np.all(np.isfinite(W_iw)):
        print("   FAIL: W_iw contains NaN or Inf!")
        return False

    if W_iw.dtype != np.float64:
        print(f"   WARNING: W_iw should be real (float64), got {W_iw.dtype}")

    print("   PASS: W_iw is valid!")

    # Check P0_iw
    if 'P0_iw' not in result:
        print("   FAIL: 'P0_iw' not in result dictionary!")
        return False

    P0_iw = result['P0_iw']
    print(f"\n   P0_iw shape: {P0_iw.shape}")
    print(f"   P0_iw[0, 0, 0] = {P0_iw[0, 0, 0]:.6e} Ha")

    if not np.all(np.isfinite(P0_iw)):
        print("   FAIL: P0_iw contains NaN or Inf!")
        return False

    print("   PASS: P0_iw is valid!")

    # Check sigma_c
    if 'sigma_c' not in result:
        print("   FAIL: 'sigma_c' not in result dictionary!")
        return False

    sigma_c = result['sigma_c']
    print(f"\n   sigma_c shape: {sigma_c.shape}")
    print(f"   sigma_c dtype: {sigma_c.dtype}")
    print(f"   sigma_c[HOMO] = {sigma_c[0]:.6e} Ha")
    print(f"   sigma_c[LUMO] = {sigma_c[1]:.6e} Ha")
    print(f"   |Im(sigma_c)|_max = {np.abs(sigma_c.imag).max():.3e}")

    # Validate sigma_c
    if len(sigma_c.shape) != 1:
        print(f"   FAIL: sigma_c should be 1D, got shape {sigma_c.shape}")
        return False

    if not np.all(np.isfinite(sigma_c)):
        print("   FAIL: sigma_c contains NaN or Inf!")
        return False

    if sigma_c.dtype != np.complex128:
        print(f"   WARNING: sigma_c should be complex128, got {sigma_c.dtype}")

    print("   PASS: sigma_c is valid!")

    # Check freqs and weights
    print("\n5. Checking frequency grid...")
    freqs_iw = result['freqs_iw']
    weights_iw = result['weights_iw']
    print(f"   freqs_iw shape: {freqs_iw.shape}")
    print(f"   freqs_iw range: [{freqs_iw.min():.3e}, {freqs_iw.max():.3e}] Ha")
    print(f"   weights_iw sum: {weights_iw.sum():.6f} (should be close to 2.0)")

    # Summary
    print("\n" + "=" * 70)
    print(" Test Summary")
    print("=" * 70)
    print(f" System: H₂ (STO-3G)")
    print(f" Auxiliary basis: {result['auxbasis']}")
    print(f" Number of frequencies: {result['nfreq']}")
    print(f" Auxiliary basis dimension: {W_iw.shape[1]}")
    print(f" Number of MOs: {sigma_c.shape[0]}")
    print()
    print(" Extracted values:")
    print(f"   W_iw:      [{result['nfreq']}, {W_iw.shape[1]}, {W_iw.shape[2]}] (real)")
    print(f"   P0_iw:     [{result['nfreq']}, {P0_iw.shape[1]}, {P0_iw.shape[2]}] (real)")
    print(f"   sigma_c:   [{sigma_c.shape[0]}] (complex)")
    print(f"   freqs_iw:  [{freqs_iw.shape[0]}]")
    print(f"   weights_iw:[{weights_iw.shape[0]}]")
    print()
    print(" QuasiX results:")
    print(f"   IP:  {result['ip']:.3f} eV")
    print(f"   Gap: {result['gap']:.3f} eV")
    print("=" * 70)
    print(" SUCCESS: All extraction tests passed!")
    print("=" * 70)

    return True

if __name__ == "__main__":
    success = test_w_sigma_c_extraction()
    sys.exit(0 if success else 1)
