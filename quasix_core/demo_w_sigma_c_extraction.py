#!/usr/bin/env python3
"""Demonstrate extraction of W(iω) and Σᶜ from QuasiX.

This script shows that W(iω) and Σᶜ are now accessible from the
run_g0w0() function for validation purposes.

The extracted values can be saved and compared with PySCF reference
data using external comparison scripts.
"""

import sys
import os
import numpy as np

# Add quasix module to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

def demo_extraction():
    """Demonstrate W(iω) and Σᶜ extraction."""
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
        return False

    print("=" * 80)
    print(" QuasiX W(iω) and Σᶜ Extraction Demo")
    print("=" * 80)

    # Create H₂ molecule
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

    # Run G₀W₀
    print("\n3. Running QuasiX G₀W₀ (nfreq=100)...")
    result = run_g0w0(mol, mf, nfreq=100, verbose=1)

    # Extract W(iω)
    print("\n4. Extracted W(iω) screening matrix:")
    W_iw = result['W_iw']
    print(f"   Shape: {W_iw.shape}")
    print(f"   Dtype: {W_iw.dtype}")
    print(f"   Memory: {W_iw.nbytes / 1024**2:.2f} MB")
    print(f"\n   Sample values at first frequency (ω₀ = {result['freqs_iw'][0]:.6f} Ha):")
    print(f"   W[0, 0, 0] = {W_iw[0, 0, 0]:.6e} Ha")
    print(f"   W[0, 0, 1] = {W_iw[0, 0, 1]:.6e} Ha")
    print(f"   W[0, 1, 1] = {W_iw[0, 1, 1]:.6e} Ha")
    print(f"\n   Sample values at mid frequency (ω₅₀ = {result['freqs_iw'][50]:.6f} Ha):")
    print(f"   W[50, 0, 0] = {W_iw[50, 0, 0]:.6e} Ha")
    print(f"\n   Sample values at high frequency (ω₉₉ = {result['freqs_iw'][99]:.6f} Ha):")
    print(f"   W[99, 0, 0] = {W_iw[99, 0, 0]:.6e} Ha")

    # Extract P₀(iω)
    print("\n5. Extracted P₀(iω) polarizability:")
    P0_iw = result['P0_iw']
    print(f"   Shape: {P0_iw.shape}")
    print(f"   P0[0, 0, 0] = {P0_iw[0, 0, 0]:.6e} Ha")

    # Extract Σᶜ
    print("\n6. Extracted Σᶜ correlation self-energy:")
    sigma_c = result['sigma_c']
    n_occ = result['n_occ']
    print(f"   Shape: {sigma_c.shape}")
    print(f"   Dtype: {sigma_c.dtype}")
    print(f"\n   Values:")
    print(f"   Σᶜ[HOMO] = {sigma_c[n_occ-1]:.6e} Ha")
    print(f"           = {sigma_c[n_occ-1].real:.6e} + {sigma_c[n_occ-1].imag:.6e}j Ha")
    print(f"   Σᶜ[LUMO] = {sigma_c[n_occ]:.6e} Ha")
    print(f"           = {sigma_c[n_occ].real:.6e} + {sigma_c[n_occ].imag:.6e}j Ha")
    print(f"\n   Imaginary part:")
    print(f"   |Im(Σᶜ)|_max = {np.abs(sigma_c.imag).max():.3e} Ha")

    # Extract exchange
    print("\n7. Extracted Σˣ exchange self-energy:")
    sigma_x = result['sigma_x']
    print(f"   Σˣ[HOMO] = {sigma_x[n_occ-1]:.6e} Ha")
    print(f"   Σˣ[LUMO] = {sigma_x[n_occ]:.6e} Ha")

    # Show total self-energy
    print("\n8. Total self-energy (Σ = Σˣ + Σᶜ):")
    sigma_total = sigma_x + sigma_c.real
    print(f"   Σ[HOMO] = {sigma_total[n_occ-1]:.6e} Ha")
    print(f"   Σ[LUMO] = {sigma_total[n_occ]:.6e} Ha")

    # Show quasiparticle correction
    print("\n9. Quasiparticle correction (ΔE = E_QP - ε_HF):")
    corrections = result['corrections']
    print(f"   ΔE[HOMO] = {corrections[n_occ-1]:.6e} Ha = {corrections[n_occ-1] * 27.2114:.3f} eV")
    print(f"   ΔE[LUMO] = {corrections[n_occ]:.6e} Ha = {corrections[n_occ] * 27.2114:.3f} eV")

    # Frequency grid info
    print("\n10. Frequency grid information:")
    freqs = result['freqs_iw']
    weights = result['weights_iw']
    print(f"    Number of points: {len(freqs)}")
    print(f"    Range: [{freqs.min():.3e}, {freqs.max():.3e}] Ha")
    print(f"    Weight sum: {weights.sum():.6f}")

    # Summary
    print("\n" + "=" * 80)
    print(" Summary: Validation-Ready Data")
    print("=" * 80)
    print(" The following arrays are now accessible for validation:")
    print()
    print(f"   result['W_iw']       → Screened interaction W(iω)")
    print(f"                           Shape: {W_iw.shape}")
    print(f"                           Use for: Element-wise comparison with PySCF W")
    print()
    print(f"   result['P0_iw']      → RPA polarizability P₀(iω)")
    print(f"                           Shape: {P0_iw.shape}")
    print(f"                           Use for: Verify P₀ computation")
    print()
    print(f"   result['sigma_c']    → Correlation self-energy Σᶜ")
    print(f"                           Shape: {sigma_c.shape}")
    print(f"                           Use for: Compare HOMO/LUMO Σᶜ with PySCF")
    print()
    print(f"   result['sigma_x']    → Exchange self-energy Σˣ")
    print(f"                           Shape: {sigma_x.shape}")
    print(f"                           Use for: Verify Σˣ = V_xc for HF")
    print()
    print(f"   result['freqs_iw']   → Imaginary frequency grid")
    print(f"   result['weights_iw'] → Quadrature weights")
    print()
    print(" Recommended validation workflow:")
    print(" 1. Save QuasiX data: np.save('quasix_W.npy', result['W_iw'])")
    print(" 2. Extract PySCF reference using validation scripts")
    print(" 3. Compare element-wise: |W_QuasiX - W_PySCF|")
    print(" 4. Identify discrepancies and debug")
    print("=" * 80)

    return True

if __name__ == "__main__":
    success = demo_extraction()
    sys.exit(0 if success else 1)
