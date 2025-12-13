#!/usr/bin/env python3
"""Test script for G₀W₀ driver."""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
from pyscf import gto, scf

from python.g0w0_driver import G0W0Driver, G0W0Config


def test_h2o():
    """Test G₀W₀ on water molecule."""
    print("\n" + "="*60)
    print("Testing G₀W₀ on H2O/cc-pVDZ")
    print("="*60)

    # Create water molecule
    mol = gto.M(
        atom='O 0 0 0; H 0 0.757 0.587; H 0 -0.757 0.587',
        basis='cc-pvdz',
        verbose=0
    )

    # Run HF calculation
    mf = scf.RHF(mol)
    mf.kernel()
    print(f"HF Energy: {mf.e_tot:.6f} Ha")

    # Run G₀W₀
    config = G0W0Config(
        auxbasis='cc-pvdz-jkfit',
        nfreq=32,
        verbose=1
    )

    gw_driver = G0W0Driver(mf, config)
    result = gw_driver.kernel()

    # Print summary
    print("\nResults Summary:")
    print(f"  HF HOMO:    {mf.mo_energy[4]*27.211:.3f} eV")
    print(f"  G₀W₀ IP:    {result.ip:.3f} eV")
    print(f"  Z(HOMO):    {result.z_factors[result.homo_index]:.3f}")

    # Reference values
    print("\nReference values:")
    print(f"  PySCF G₀W₀: ~12.2 eV")
    print(f"  Experiment: 12.62 eV")

    return result


def test_nh3():
    """Test G₀W₀ on ammonia molecule."""
    print("\n" + "="*60)
    print("Testing G₀W₀ on NH3/cc-pVDZ")
    print("="*60)

    # Create ammonia molecule
    mol = gto.M(
        atom='''
        N  0.0000  0.0000  0.1173
        H  0.0000  0.9320 -0.2738
        H  0.8068 -0.4660 -0.2738
        H -0.8068 -0.4660 -0.2738
        ''',
        basis='cc-pvdz',
        verbose=0
    )

    # Run HF calculation
    mf = scf.RHF(mol)
    mf.kernel()
    print(f"HF Energy: {mf.e_tot:.6f} Ha")

    # Run G₀W₀
    config = G0W0Config(
        auxbasis='cc-pvdz-jkfit',
        nfreq=32,
        verbose=1
    )

    gw_driver = G0W0Driver(mf, config)
    result = gw_driver.kernel()

    # Print summary
    homo_idx = mol.nelectron // 2 - 1
    print("\nResults Summary:")
    print(f"  HF HOMO:    {mf.mo_energy[homo_idx]*27.211:.3f} eV")
    print(f"  G₀W₀ IP:    {result.ip:.3f} eV")
    print(f"  Z(HOMO):    {result.z_factors[result.homo_index]:.3f}")

    # Reference values
    print("\nReference values:")
    print(f"  Experiment: 10.85 eV")

    return result


def test_co():
    """Test G₀W₀ on carbon monoxide."""
    print("\n" + "="*60)
    print("Testing G₀W₀ on CO/cc-pVDZ")
    print("="*60)

    # Create CO molecule
    mol = gto.M(
        atom='C 0 0 0; O 0 0 1.128',
        basis='cc-pvdz',
        verbose=0
    )

    # Run HF calculation
    mf = scf.RHF(mol)
    mf.kernel()
    print(f"HF Energy: {mf.e_tot:.6f} Ha")

    # Run G₀W₀
    config = G0W0Config(
        auxbasis='cc-pvdz-jkfit',
        nfreq=32,
        verbose=1
    )

    gw_driver = G0W0Driver(mf, config)
    result = gw_driver.kernel()

    # Print summary
    homo_idx = mol.nelectron // 2 - 1
    print("\nResults Summary:")
    print(f"  HF HOMO:    {mf.mo_energy[homo_idx]*27.211:.3f} eV")
    print(f"  G₀W₀ IP:    {result.ip:.3f} eV")
    print(f"  Z(HOMO):    {result.z_factors[result.homo_index]:.3f}")

    # Reference values
    print("\nReference values:")
    print(f"  Experiment: 14.01 eV")

    return result


if __name__ == "__main__":
    # Test on small molecules
    results = {}

    results['h2o'] = test_h2o()
    results['nh3'] = test_nh3()
    results['co'] = test_co()

    # Summary table
    print("\n" + "="*60)
    print("Summary of G₀W₀ Results")
    print("="*60)
    print(f"{'Molecule':<10} {'IP (eV)':<12} {'Z(HOMO)':<10} {'Σˣ(HOMO) eV':<12}")
    print("-"*44)
    for mol_name, result in results.items():
        ip = result.ip
        z_homo = result.z_factors[result.homo_index]
        sigma_x = result.sigma_x[result.homo_index] * 27.211
        print(f"{mol_name.upper():<10} {ip:<12.3f} {z_homo:<10.3f} {sigma_x:<12.3f}")
    print("="*60)

    print("\nNote: Correlation self-energy is approximated.")
    print("Full implementation requires proper W calculation.")