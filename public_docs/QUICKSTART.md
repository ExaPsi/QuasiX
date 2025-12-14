# Quick Start Guide

This guide will walk you through your first QuasiX calculation in under 5 minutes.

## Prerequisites

Ensure QuasiX is installed (see [INSTALL.md](INSTALL.md)).

## Example 1: G₀W₀ Calculation for Water

```python
from quasix import G0W0Driver
from pyscf import gto, scf

# Define water molecule
mol = gto.M(
    atom='''
    O  0.000000  0.000000  0.117369
    H  0.000000  0.756950 -0.469476
    H  0.000000 -0.756950 -0.469476
    ''',
    basis='def2-svp',
    unit='angstrom'
)

# Run Hartree-Fock
mf = scf.RHF(mol)
mf.kernel()

# Run G0W0
gw = G0W0Driver(mf)
result = gw.kernel()

# Print results
print(f"HOMO (DFT):  {result.homo_dft:.3f} eV")
print(f"HOMO (G0W0): {result.homo_qp:.3f} eV")
print(f"LUMO (G0W0): {result.lumo_qp:.3f} eV")
print(f"Gap (G0W0):  {result.gap_qp:.3f} eV")
```

## Example 2: evGW for Better Accuracy

```python
from quasix import evGWDriver
from pyscf import gto, scf, dft

# Use DFT starting point
mol = gto.M(atom='O 0 0 0; H 0 0.757 0.587; H 0 -0.757 0.587', basis='def2-tzvp')
mf = dft.RKS(mol)
mf.xc = 'pbe'
mf.kernel()

# Run eigenvalue self-consistent GW
gw = evGWDriver(mf, max_iter=10, conv_tol=1e-4)
result = gw.kernel()

print(f"evGW HOMO: {result.homo_qp:.3f} eV")
print(f"Converged in {result.n_iter} iterations")
```

## Example 3: BSE Optical Spectrum

```python
from quasix import BSEDriver
from pyscf import gto, scf

# Benzene molecule
mol = gto.M(
    atom='''
    C  1.3970  0.0000  0.0000
    C  0.6985  1.2098  0.0000
    C -0.6985  1.2098  0.0000
    C -1.3970  0.0000  0.0000
    C -0.6985 -1.2098  0.0000
    C  0.6985 -1.2098  0.0000
    H  2.4810  0.0000  0.0000
    H  1.2405  2.1486  0.0000
    H -1.2405  2.1486  0.0000
    H -2.4810  0.0000  0.0000
    H -1.2405 -2.1486  0.0000
    H  1.2405 -2.1486  0.0000
    ''',
    basis='def2-svp'
)
mf = scf.RHF(mol).run()

# Run BSE
bse = BSEDriver(mf, n_states=10, spin='singlet')
result = bse.kernel()

# Print excitation energies
print("Excitation Energies (eV):")
for i, (E, f) in enumerate(zip(result.excitation_energies, result.oscillator_strengths)):
    print(f"  S{i+1}: {E:.3f} eV  (f = {f:.4f})")

# Generate absorption spectrum
spectrum = bse.get_spectrum(broadening=0.1)
```

## Example 4: Custom Frequency Grid

```python
from quasix import G0W0Driver

gw = G0W0Driver(
    mf,
    n_freq=64,              # Number of frequency points
    freq_method='minimax',  # Frequency integration method
    eta=0.01,               # Broadening parameter
)
result = gw.kernel()
```

## Example 5: Parallel Execution

```python
import os
os.environ['RAYON_NUM_THREADS'] = '8'  # Set before import

from quasix import G0W0Driver

gw = G0W0Driver(mf, n_workers=8)
result = gw.kernel()
```

## Output Analysis

### Quasiparticle Energies

```python
# Access all QP energies
print("All QP energies (eV):")
for i, e in enumerate(result.qp_energies):
    print(f"  Orbital {i}: {e:.4f} eV")

# Get frontier orbitals
print(f"\nHOMO: {result.homo_qp:.4f} eV")
print(f"LUMO: {result.lumo_qp:.4f} eV")
print(f"Gap:  {result.gap_qp:.4f} eV")
```

### Self-Energy Components

```python
# Exchange self-energy
print(f"Σx(HOMO): {result.sigma_x[result.homo_idx]:.4f} eV")

# Correlation self-energy
print(f"Σc(HOMO): {result.sigma_c[result.homo_idx]:.4f} eV")

# Renormalization factor
print(f"Z(HOMO): {result.z_factor[result.homo_idx]:.4f}")
```

## Next Steps

- [API Reference](API.md) for detailed function documentation
- [Theory](THEORY.md) for understanding the methods
- [Benchmarks](BENCHMARKS.md) for accuracy and performance data
