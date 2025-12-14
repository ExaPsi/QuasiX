# API Reference

## Core Classes

### G0W0Driver

One-shot G₀W₀ calculations.

```python
class G0W0Driver:
    def __init__(
        self,
        mf,                          # PySCF mean-field object (RHF/RKS)
        basis_aux: str = None,       # Auxiliary basis for RI (default: auto)
        n_freq: int = 32,            # Number of frequency points
        freq_method: str = 'minimax', # 'minimax' or 'gauss-legendre'
        eta: float = 0.001,          # Broadening parameter (Ha)
        frozen_core: int = 0,        # Number of frozen core orbitals
        qp_solver: str = 'linearized', # 'linearized' or 'newton'
    ):
        ...

    def kernel(self) -> G0W0Result:
        """Run G0W0 calculation."""
        ...
```

#### G0W0Result

```python
class G0W0Result:
    qp_energies: np.ndarray    # Quasiparticle energies (eV)
    dft_energies: np.ndarray   # DFT/HF orbital energies (eV)
    sigma_x: np.ndarray        # Exchange self-energy (eV)
    sigma_c: np.ndarray        # Correlation self-energy (eV)
    z_factor: np.ndarray       # Renormalization factors
    homo_idx: int              # HOMO orbital index
    lumo_idx: int              # LUMO orbital index
    homo_qp: float             # HOMO quasiparticle energy (eV)
    lumo_qp: float             # LUMO quasiparticle energy (eV)
    gap_qp: float              # QP band gap (eV)
    homo_dft: float            # HOMO DFT energy (eV)
    converged: bool            # Convergence status
```

---

### evGWDriver

Eigenvalue self-consistent GW calculations.

```python
class evGWDriver:
    def __init__(
        self,
        mf,                          # PySCF mean-field object
        max_iter: int = 30,          # Maximum iterations
        conv_tol: float = 1e-5,      # Convergence tolerance (eV)
        mixing: float = 0.5,         # Mixing parameter for damping
        **g0w0_kwargs                # Additional G0W0 arguments
    ):
        ...

    def kernel(self) -> evGWResult:
        """Run evGW calculation."""
        ...
```

#### evGWResult

Extends `G0W0Result` with:
```python
class evGWResult(G0W0Result):
    n_iter: int                # Number of iterations
    convergence_history: list  # Energy changes per iteration
```

---

### BSEDriver

Bethe-Salpeter Equation for optical excitations.

```python
class BSEDriver:
    def __init__(
        self,
        mf,                          # PySCF mean-field object
        gw_result: G0W0Result = None, # Pre-computed GW result (optional)
        n_states: int = 10,          # Number of excited states
        spin: str = 'singlet',       # 'singlet' or 'triplet'
        tda: bool = False,           # Tamm-Dancoff approximation
        n_occ: int = None,           # Active occupied orbitals
        n_vir: int = None,           # Active virtual orbitals
    ):
        ...

    def kernel(self) -> BSEResult:
        """Solve BSE eigenvalue problem."""
        ...

    def get_spectrum(
        self,
        energy_range: tuple = (0, 20),  # Energy range (eV)
        n_points: int = 1000,            # Number of points
        broadening: float = 0.1,         # Lorentzian broadening (eV)
    ) -> tuple[np.ndarray, np.ndarray]:
        """Generate absorption spectrum."""
        ...
```

#### BSEResult

```python
class BSEResult:
    excitation_energies: np.ndarray  # Excitation energies (eV)
    oscillator_strengths: np.ndarray # Oscillator strengths
    eigenvectors: np.ndarray         # BSE eigenvectors
    transition_dipoles: np.ndarray   # Transition dipole moments
    n_states: int                    # Number of states computed
```

---

## Utility Functions

### Basis Set Handling

```python
from quasix.basis import get_auxiliary_basis, list_available_basis

# Get recommended auxiliary basis
aux_basis = get_auxiliary_basis('def2-svp')  # Returns 'def2-svp-ri'

# List available basis sets
basis_sets = list_available_basis()
```

### Unit Conversion

```python
from quasix.units import eV_to_Ha, Ha_to_eV, Angstrom_to_Bohr

energy_ha = eV_to_Ha(13.6)    # Convert eV to Hartree
energy_ev = Ha_to_eV(0.5)     # Convert Hartree to eV
```

### Analysis Tools

```python
from quasix.analysis import compute_dos, compute_pdos

# Density of states
energies, dos = compute_dos(result.qp_energies, broadening=0.1)

# Projected DOS
energies, pdos = compute_pdos(result, mo_coeffs, ao_labels)
```

---

## Advanced Configuration

### Frequency Integration

```python
from quasix import G0W0Driver

# Minimax grid (default, most efficient)
gw = G0W0Driver(mf, freq_method='minimax', n_freq=32)

# Gauss-Legendre grid (for comparison)
gw = G0W0Driver(mf, freq_method='gauss-legendre', n_freq=64)
```

### QP Solver Options

```python
# Linearized QP equation (fast, usually sufficient)
gw = G0W0Driver(mf, qp_solver='linearized')

# Newton solver (more accurate for difficult cases)
gw = G0W0Driver(mf, qp_solver='newton', newton_tol=1e-6)
```

### Memory Management

```python
# For large systems, use out-of-core mode
gw = G0W0Driver(mf, incore=False, max_memory=4000)  # 4 GB limit
```

---

## Error Handling

```python
from quasix.exceptions import ConvergenceError, BasisError

try:
    result = gw.kernel()
except ConvergenceError as e:
    print(f"GW did not converge: {e}")
except BasisError as e:
    print(f"Basis set error: {e}")
```

---

## Examples

See [QUICKSTART.md](QUICKSTART.md) for complete working examples.
