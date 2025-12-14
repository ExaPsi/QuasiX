#!/usr/bin/env python3
"""G₀W₀ (Single-Shot GW) Driver Implementation.

This module provides a high-level interface to QuasiX G₀W₀ calculations,
which computes quasiparticle energies without self-consistency iteration.

G₀W₀ is mathematically equivalent to evGW with max_cycle=1, providing:
- Fast screening calculations for large molecular datasets
- Standard reference method for GW benchmarks (e.g., GW100)
- Foundation for more advanced methods (BSE, G₀W₀Γ₀)

Mathematical Foundation:
    E_n^QP = ε_n + Z_n × [Σ_n^x + Σ_n^c(ε_n) - v_xc,n]

where:
    - ε_n: HF or DFT eigenvalues (starting point)
    - Σ_n^x: Exchange self-energy
    - Σ_n^c(ω): Correlation self-energy (frequency-dependent)
    - v_xc,n: Exchange-correlation potential
    - Z_n: Renormalization factor (quasiparticle weight)

Physical Constraints:
    - 0 < Z_n < 1 (quasiparticle spectral weight)
    - IP = -E_HOMO > 0 (ionization potential positive)
    - Gap = E_LUMO - E_HOMO > 0 (for molecules)

References:
    - Theory: docs/derivations/g0w0/g0w0_theory.md
    - Algorithms: docs/derivations/g0w0/g0w0_algorithms.md
    - Validation: docs/verification/g0w0_verification_plan.md
"""

import numpy as np
import time
import warnings
from typing import Optional, Dict, Any
from dataclasses import dataclass, field
from pathlib import Path

try:
    from pyscf import gto, scf, df
    HAS_PYSCF = True
except ImportError:
    HAS_PYSCF = False
    warnings.warn("PySCF not available - G₀W₀ calculations require PySCF")

# Import QuasiX evGW driver (G₀W₀ is a special case)
try:
    from .gw_driver import GWDriver, GWConfig, GWResult
    QUASIX_AVAILABLE = True
except ImportError:
    try:
        from gw_driver import GWDriver, GWConfig, GWResult
        QUASIX_AVAILABLE = True
    except ImportError:
        QUASIX_AVAILABLE = False
        warnings.warn("QuasiX GW driver not available")

# Physical constants
HA_TO_EV = 27.211386245988


@dataclass
class G0W0Config:
    """Configuration for G₀W₀ calculations.

    Attributes:
        auxbasis: Auxiliary basis set for RI/DF (e.g., 'cc-pvdz-jkfit')
        nfreq: Number of frequency grid points for contour deformation
        eta: Broadening parameter (Ha) - prevents poles in denominators
        delta_z: Finite difference step for Z-factor derivatives (Ha)
        verbose: Print level (0=quiet, 1=normal, 2=debug)
        validate_z: Check Z-factors are in physical range (0, 1)
        freq_int: Frequency integration method ('cd' for contour deformation)

    Notes:
        - Default parameters are validated for molecules (not solids)
        - For large basis sets, increase nfreq to 64 for better accuracy
        - eta=0.01 Ha provides good compromise between accuracy and stability
    """
    auxbasis: str = 'cc-pvdz-jkfit'
    nfreq: int = 32
    eta: float = 0.01
    delta_z: float = 0.001
    verbose: int = 1
    validate_z: bool = True
    freq_int: str = 'cd'  # Contour deformation (robust, accurate)


@dataclass
class G0W0Result:
    """Results from G₀W₀ calculation.

    Attributes:
        converged: Always True for G₀W₀ (one-shot method)
        qp_energies: Quasiparticle energies (Ha)
        z_factors: Renormalization factors (dimensionless, 0 < Z < 1)
        sigma_x: Exchange self-energy diagonal (Ha)
        sigma_c_real: Real part of correlation self-energy (Ha)
        sigma_c_imag: Imaginary part of correlation self-energy (Ha)
        homo_index: Index of highest occupied orbital
        lumo_index: Index of lowest unoccupied orbital
        ip: Ionization potential (eV) = -E_HOMO
        ea: Electron affinity (eV) = -E_LUMO (None if no virtual orbitals)
        gap: HOMO-LUMO gap (eV)
        elapsed_seconds: Wall-clock time for calculation
        method: Always "G0W0"
        vxc_diag: Exchange-correlation potential diagonal (Ha)
        implementation: "Rust" or "Python" backend used

    Methods:
        get_homo_energy: HOMO energy in eV
        get_lumo_energy: LUMO energy in eV
        get_qp_corrections: QP corrections (QP - HF) in eV
        to_dict: Convert to dictionary for serialization
    """
    converged: bool = True  # G₀W₀ always "converged" (one-shot)
    qp_energies: np.ndarray = None
    z_factors: np.ndarray = None
    sigma_x: np.ndarray = None
    sigma_c_real: np.ndarray = None
    sigma_c_imag: np.ndarray = None
    homo_index: int = None
    lumo_index: int = None
    ip: float = None  # eV
    ea: Optional[float] = None  # eV
    gap: float = None  # eV
    elapsed_seconds: float = None
    method: str = "G0W0"
    vxc_diag: np.ndarray = None
    implementation: str = "Rust"

    def get_homo_energy(self) -> float:
        """Get HOMO energy in eV."""
        if self.qp_energies is None or self.homo_index is None:
            raise ValueError("QP energies not available")
        return self.qp_energies[self.homo_index] * HA_TO_EV

    def get_lumo_energy(self) -> float:
        """Get LUMO energy in eV."""
        if self.qp_energies is None or self.lumo_index is None:
            raise ValueError("QP energies not available")
        return self.qp_energies[self.lumo_index] * HA_TO_EV

    def get_qp_corrections(self, mo_energy: np.ndarray) -> np.ndarray:
        """Get QP corrections (QP - HF) in eV.

        Args:
            mo_energy: HF/DFT orbital energies (Ha)

        Returns:
            QP corrections in eV for each orbital
        """
        if self.qp_energies is None:
            raise ValueError("QP energies not available")
        return (self.qp_energies - mo_energy) * HA_TO_EV

    def to_dict(self) -> Dict[str, Any]:
        """Convert result to dictionary for serialization."""
        return {
            'converged': self.converged,
            'method': self.method,
            'qp_energies_ha': self.qp_energies.tolist() if self.qp_energies is not None else None,
            'qp_energies_ev': (self.qp_energies * HA_TO_EV).tolist() if self.qp_energies is not None else None,
            'z_factors': self.z_factors.tolist() if self.z_factors is not None else None,
            'sigma_x_ha': self.sigma_x.tolist() if self.sigma_x is not None else None,
            'sigma_c_real_ha': self.sigma_c_real.tolist() if self.sigma_c_real is not None else None,
            'sigma_c_imag_ha': self.sigma_c_imag.tolist() if self.sigma_c_imag is not None else None,
            'homo_index': self.homo_index,
            'lumo_index': self.lumo_index,
            'ip_ev': self.ip,
            'ea_ev': self.ea,
            'gap_ev': self.gap,
            'elapsed_seconds': self.elapsed_seconds,
            'implementation': self.implementation,
        }


class G0W0Driver:
    """High-level driver for G₀W₀ (single-shot GW) calculations.

    This class provides a user-friendly interface to QuasiX's production-quality
    G₀W₀ implementation, which uses the same Rust kernel as evGW but with max_cycle=1.

    Mathematical Identity:
        G₀W₀ ≡ evGW(max_cycle=1)

    This identity provides an exact unit test for validation.

    Usage:
        >>> from pyscf import gto, scf
        >>> from quasix.python.g0w0_driver import G0W0Driver, G0W0Config
        >>>
        >>> # Setup molecule
        >>> mol = gto.M(atom='O 0 0 0; H 0 0 0.96; H 0 0.96 0', basis='cc-pvdz')
        >>> mf = scf.RHF(mol).run()
        >>>
        >>> # Run G₀W₀
        >>> config = G0W0Config(auxbasis='cc-pvdz-jkfit', nfreq=32, verbose=1)
        >>> gw = G0W0Driver(mf, config)
        >>> result = gw.kernel()
        >>>
        >>> # Extract results
        >>> print(f"IP: {result.ip:.3f} eV")
        >>> print(f"Gap: {result.gap:.3f} eV")
        >>> print(f"HOMO Z-factor: {result.z_factors[result.homo_index]:.3f}")

    Attributes:
        mf: PySCF mean-field object (RHF or RKS)
        config: G0W0Config object with calculation parameters
        mol: PySCF molecule object (from mf)
        result: G0W0Result object (after kernel() is called)

    Methods:
        kernel: Run G₀W₀ calculation
        analyze: Analyze results and print summary
        save_results: Save results to HDF5 file
    """

    def __init__(self, mf, config: Optional[G0W0Config] = None):
        """Initialize G₀W₀ driver.

        Args:
            mf: PySCF mean-field object (RHF or RKS, must be converged)
            config: G0W0Config object (uses defaults if not provided)

        Raises:
            RuntimeError: If PySCF or QuasiX modules are not available
            ValueError: If mean-field calculation is not converged
        """
        if not HAS_PYSCF:
            raise RuntimeError(
                "PySCF is required for G₀W₀ calculations.\n"
                "Install with: pip install pyscf"
            )

        if not QUASIX_AVAILABLE:
            raise RuntimeError(
                "QuasiX GW driver is required for G₀W₀ calculations.\n"
                "Check that quasix.python.gw_driver is installed."
            )

        # Validate mean-field object
        if not hasattr(mf, 'mol'):
            raise ValueError("Mean-field object must have 'mol' attribute")

        if not hasattr(mf, 'mo_energy') or mf.mo_energy is None:
            raise ValueError("Mean-field calculation must be converged (mo_energy is None)")

        if hasattr(mf, 'converged') and not mf.converged:
            warnings.warn("Mean-field calculation did not converge - G₀W₀ results may be unreliable")

        self.mf = mf
        self.mol = mf.mol
        self.config = config or G0W0Config()
        self.result = None

        # Store MO data
        self.mo_energy = mf.mo_energy
        self.mo_coeff = mf.mo_coeff
        self.mo_occ = mf.mo_occ
        self.n_occ = int(np.sum(self.mo_occ > 0))
        self.n_mo = len(self.mo_energy)

        # Initialize underlying evGW driver
        # G₀W₀ is just evGW with max_cycle=1
        gw_config = GWConfig(
            auxbasis=self.config.auxbasis,
            freq_int=self.config.freq_int,
            nfreq=self.config.nfreq,
            eta=self.config.eta,
            max_cycle=1,  # KEY: Single iteration for G₀W₀
            conv_tol=1e-10,  # Tight tolerance (single iteration anyway)
            conv_tol_z=1e-10,
            damping=1.0,  # No damping needed for single iteration
            damping_dynamic=False,
            diis=False,  # No DIIS needed for single iteration
            validate_z=self.config.validate_z,
            verbose=self.config.verbose,
        )

        self._gw_driver = GWDriver(
            mo_energy=self.mo_energy,
            mo_coeff=self.mo_coeff,
            n_occ=self.n_occ,
            mol=self.mol,
            mf=self.mf,
            auxbasis=self.config.auxbasis,
            config=gw_config
        )

    def kernel(self) -> G0W0Result:
        """Run G₀W₀ calculation using Rust kernel.

        This method delegates to the evGW driver with max_cycle=1, ensuring
        exact equivalence between G₀W₀ and evGW(1).

        Returns:
            G0W0Result: Results including QP energies, Z-factors, IP, EA, gap

        Raises:
            RuntimeError: If calculation fails
            ValueError: If Z-factors are outside physical range

        Notes:
            - Calculation is deterministic (no random elements)
            - Results are reproducible with same inputs
            - Wall-clock time typically < 30s for small molecules (cc-pVDZ)
        """
        start_time = time.time()

        if self.config.verbose >= 1:
            print("=" * 70)
            print("QuasiX G₀W₀ (Single-Shot GW) Calculation")
            print("=" * 70)
            print(f"System: {self.mol.atom_symbol(0)} + {self.mol.natm-1} atoms")
            print(f"Basis: {self.mol.basis}")
            print(f"Auxiliary basis: {self.config.auxbasis}")
            print(f"Orbitals: {self.n_mo} total ({self.n_occ} occupied, {self.n_mo - self.n_occ} virtual)")
            print(f"Frequency grid: {self.config.nfreq} points")
            print(f"Broadening: {self.config.eta:.4f} Ha = {self.config.eta * HA_TO_EV:.3f} eV")
            print("=" * 70)

        # Run evGW with max_cycle=1 (this IS G₀W₀)
        try:
            gw_result = self._gw_driver.kernel()
        except Exception as e:
            raise RuntimeError(
                f"G₀W₀ calculation failed.\n"
                f"Error: {e}\n"
                f"This indicates a problem with the Rust kernel or input data."
            ) from e

        elapsed = time.time() - start_time

        # Convert to G0W0Result format
        result = G0W0Result(
            converged=True,  # G₀W₀ always converged (one-shot)
            qp_energies=gw_result.qp_energies,
            z_factors=gw_result.z_factors,
            sigma_x=gw_result.sigma_x,
            sigma_c_real=np.real(gw_result.sigma_c),
            sigma_c_imag=np.imag(gw_result.sigma_c),
            homo_index=gw_result.homo_index,
            lumo_index=gw_result.lumo_index,
            ip=gw_result.ip,  # Already in eV from GWResult
            ea=gw_result.ea,  # Already in eV from GWResult
            gap=gw_result.gap,  # Already in eV from GWResult
            elapsed_seconds=elapsed,
            implementation=getattr(gw_result, 'implementation', 'Rust'),
        )

        # Store result
        self.result = result

        # Validate Z-factors if requested
        if self.config.validate_z:
            self._validate_z_factors(result.z_factors)

        # Physical sanity checks
        self._check_physical_sanity(result)

        if self.config.verbose >= 1:
            print("\n" + "=" * 70)
            print("G₀W₀ Calculation Complete")
            print("=" * 70)
            print(f"Ionization Potential: {result.ip:.3f} eV")
            if result.ea is not None:
                print(f"Electron Affinity: {result.ea:.3f} eV")
            print(f"HOMO-LUMO Gap: {result.gap:.3f} eV")
            print(f"HOMO Z-factor: {result.z_factors[result.homo_index]:.4f}")
            print(f"LUMO Z-factor: {result.z_factors[result.lumo_index]:.4f}")
            print(f"Elapsed time: {elapsed:.2f} seconds")
            print(f"Implementation: {result.implementation}")
            print("=" * 70)

        return result

    def _validate_z_factors(self, z_factors: np.ndarray):
        """Validate Z-factors are in physical range (0, 1).

        Args:
            z_factors: Array of renormalization factors

        Raises:
            ValueError: If any Z-factors are outside (0, 1)

        Notes:
            - Z < 0: Unphysical (indicates numerical error)
            - Z > 1: Unphysical (quasiparticle weight > 1)
            - Z < 0.3: Warning (weak quasiparticle, strong correlation)
            - Z > 0.95: Warning (nearly free particle, GW may be overkill)
        """
        bad_indices = np.where((z_factors <= 0) | (z_factors > 1))[0]

        if len(bad_indices) > 0:
            error_msg = (
                f"Unphysical Z-factors detected at {len(bad_indices)} orbitals:\n"
            )
            for idx in bad_indices[:5]:  # Show first 5
                error_msg += f"  Orbital {idx}: Z = {z_factors[idx]:.6f}\n"
            if len(bad_indices) > 5:
                error_msg += f"  ... and {len(bad_indices) - 5} more\n"
            error_msg += "\nPossible causes:\n"
            error_msg += "  - Insufficient frequency grid (try nfreq=64)\n"
            error_msg += "  - Numerical instability in derivative calculation\n"
            error_msg += "  - System too metallic for G₀W₀ method\n"

            raise ValueError(error_msg)

        # Warnings for unusual but not necessarily wrong values
        weak_qp = np.where(z_factors < 0.3)[0]
        if len(weak_qp) > 0:
            warnings.warn(
                f"{len(weak_qp)} orbitals have Z < 0.3 (weak quasiparticle character). "
                "This indicates strong correlation effects."
            )

        strong_qp = np.where(z_factors > 0.95)[0]
        if len(strong_qp) > 0 and self.config.verbose >= 2:
            print(f"Note: {len(strong_qp)} orbitals have Z > 0.95 (nearly free particle)")

    def _check_physical_sanity(self, result: G0W0Result):
        """Perform physical sanity checks on results.

        Args:
            result: G0W0Result object to validate

        Warnings are issued for suspicious values, but calculation is not aborted.
        """
        # Check IP is positive
        if result.ip < 0:
            warnings.warn(
                f"Negative ionization potential: {result.ip:.3f} eV. "
                "This is unphysical for molecules."
            )

        # Check gap is positive
        if result.gap < 0:
            warnings.warn(
                f"Negative HOMO-LUMO gap: {result.gap:.3f} eV. "
                "This indicates orbital crossing or numerical error."
            )

        # Check QP corrections are reasonable
        qp_corrections = (result.qp_energies - self.mo_energy) * HA_TO_EV
        max_correction = np.max(np.abs(qp_corrections))

        if max_correction > 10.0:
            warnings.warn(
                f"Large QP correction detected: {max_correction:.2f} eV. "
                "G₀W₀ perturbative treatment may not be valid."
            )

        # Check for NaNs or Infs
        for name, array in [
            ('qp_energies', result.qp_energies),
            ('z_factors', result.z_factors),
            ('sigma_x', result.sigma_x),
        ]:
            if np.any(np.isnan(array)) or np.any(np.isinf(array)):
                raise ValueError(f"{name} contains NaN or Inf values!")

    def analyze(self, print_orbitals: bool = True, n_orbitals: int = 5):
        """Analyze and print G₀W₀ results.

        Args:
            print_orbitals: Print orbital-by-orbital breakdown
            n_orbitals: Number of frontier orbitals to print (around HOMO/LUMO)

        Raises:
            RuntimeError: If no calculation has been performed
        """
        if self.result is None:
            raise RuntimeError("No calculation has been performed. Call kernel() first.")

        r = self.result

        print("\n" + "=" * 70)
        print("G₀W₀ Analysis")
        print("=" * 70)

        # Summary
        print("\nQuasiparticle Energies (Frontier Orbitals):")
        print(f"  Ionization Potential: {r.ip:.4f} eV (= -E_HOMO)")
        if r.ea is not None:
            print(f"  Electron Affinity: {r.ea:.4f} eV (= -E_LUMO)")
        print(f"  HOMO-LUMO Gap: {r.gap:.4f} eV")

        print("\nRenormalization Factors:")
        print(f"  HOMO (orbital {r.homo_index}): Z = {r.z_factors[r.homo_index]:.4f}")
        print(f"  LUMO (orbital {r.lumo_index}): Z = {r.z_factors[r.lumo_index]:.4f}")
        print(f"  Average (all orbitals): Z_avg = {np.mean(r.z_factors):.4f}")
        print(f"  Range: [{np.min(r.z_factors):.4f}, {np.max(r.z_factors):.4f}]")

        # QP corrections
        qp_corrections = r.get_qp_corrections(self.mo_energy)
        print("\nQuasiparticle Corrections (QP - HF):")
        print(f"  HOMO: {qp_corrections[r.homo_index]:+.4f} eV")
        print(f"  LUMO: {qp_corrections[r.lumo_index]:+.4f} eV")
        print(f"  Max correction: {np.max(np.abs(qp_corrections)):.4f} eV")

        # Orbital breakdown
        if print_orbitals:
            print(f"\nFrontier Orbitals (±{n_orbitals} around HOMO/LUMO):")
            print(f"{'Orb':>5} {'HF (eV)':>10} {'QP (eV)':>10} {'Corr (eV)':>11} {'Z':>8} {'Type':>6}")
            print("-" * 70)

            start_idx = max(0, r.homo_index - n_orbitals)
            end_idx = min(self.n_mo, r.lumo_index + n_orbitals + 1)

            for i in range(start_idx, end_idx):
                hf_e = self.mo_energy[i] * HA_TO_EV
                qp_e = r.qp_energies[i] * HA_TO_EV
                corr = qp_corrections[i]
                z = r.z_factors[i]

                if i == r.homo_index:
                    orb_type = "HOMO"
                elif i == r.lumo_index:
                    orb_type = "LUMO"
                elif i < r.homo_index:
                    orb_type = "Occ"
                else:
                    orb_type = "Virt"

                print(f"{i:5d} {hf_e:10.4f} {qp_e:10.4f} {corr:+11.4f} {z:8.4f} {orb_type:>6}")

        print("=" * 70)

    def save_results(self, filepath: str):
        """Save G₀W₀ results to HDF5 file with full provenance.

        Args:
            filepath: Path to HDF5 file (will be created or overwritten)

        Raises:
            RuntimeError: If no calculation has been performed
            ImportError: If h5py is not available
        """
        if self.result is None:
            raise RuntimeError("No calculation has been performed. Call kernel() first.")

        try:
            import h5py
        except ImportError:
            raise ImportError("h5py is required to save results. Install with: pip install h5py")

        import datetime

        r = self.result

        with h5py.File(filepath, 'w') as f:
            # Results
            results_grp = f.create_group('results')
            results_grp.create_dataset('qp_energies_ha', data=r.qp_energies)
            results_grp.create_dataset('qp_energies_ev', data=r.qp_energies * HA_TO_EV)
            results_grp.create_dataset('z_factors', data=r.z_factors)
            results_grp.create_dataset('sigma_x_ha', data=r.sigma_x)
            results_grp.create_dataset('sigma_c_real_ha', data=r.sigma_c_real)
            results_grp.create_dataset('sigma_c_imag_ha', data=r.sigma_c_imag)

            # Metadata
            results_grp.attrs['converged'] = r.converged
            results_grp.attrs['method'] = r.method
            results_grp.attrs['homo_index'] = r.homo_index
            results_grp.attrs['lumo_index'] = r.lumo_index
            results_grp.attrs['ip_ev'] = r.ip
            if r.ea is not None:
                results_grp.attrs['ea_ev'] = r.ea
            results_grp.attrs['gap_ev'] = r.gap
            results_grp.attrs['elapsed_seconds'] = r.elapsed_seconds
            results_grp.attrs['implementation'] = r.implementation

            # Input data
            input_grp = f.create_group('input')
            input_grp.create_dataset('mo_energy_ha', data=self.mo_energy)
            input_grp.create_dataset('mo_coeff', data=self.mo_coeff)
            input_grp.create_dataset('mo_occ', data=self.mo_occ)
            input_grp.attrs['n_occ'] = self.n_occ
            input_grp.attrs['n_mo'] = self.n_mo

            # Configuration
            config_grp = f.create_group('config')
            config_grp.attrs['auxbasis'] = self.config.auxbasis
            config_grp.attrs['nfreq'] = self.config.nfreq
            config_grp.attrs['eta'] = self.config.eta
            config_grp.attrs['delta_z'] = self.config.delta_z
            config_grp.attrs['freq_int'] = self.config.freq_int

            # Provenance
            prov_grp = f.create_group('provenance')
            prov_grp.attrs['timestamp'] = datetime.datetime.now().isoformat()
            prov_grp.attrs['quasix_version'] = "0.2.0"  # Update with actual version
            prov_grp.attrs['method'] = 'G0W0'
            prov_grp.attrs['basis'] = str(self.mol.basis)
            prov_grp.attrs['auxbasis'] = self.config.auxbasis

            # Molecule info
            mol_grp = f.create_group('molecule')
            mol_grp.attrs['natm'] = self.mol.natm
            mol_grp.attrs['nelec'] = self.mol.nelectron
            mol_grp.attrs['spin'] = self.mol.spin
            mol_grp.attrs['charge'] = self.mol.charge

            # Save geometry
            atom_coords = np.array([self.mol.atom_coord(i) for i in range(self.mol.natm)])
            mol_grp.create_dataset('coordinates_bohr', data=atom_coords)
            atom_symbols = [self.mol.atom_symbol(i) for i in range(self.mol.natm)]
            mol_grp.attrs['atom_symbols'] = ','.join(atom_symbols)

        if self.config.verbose >= 1:
            print(f"Results saved to: {filepath}")


def load_gw100_structure(cas_number: str) -> str:
    """Load molecular structure from GW100 dataset.

    Args:
        cas_number: CAS registry number (e.g., '7732-18-5' for water)

    Returns:
        XYZ geometry string in PySCF format

    Raises:
        FileNotFoundError: If structure file not found

    Notes:
        - Structures are from GW100 GOLD STANDARD (experimental from HCP92)
        - Coordinates are in Angstrom
        - Always use these exact geometries for GW100 benchmarks

    Example:
        >>> water_geom = load_gw100_structure('7732-18-5')
        >>> mol = gto.M(atom=water_geom, basis='def2-tzvp', unit='angstrom')
    """
    xyz_file = Path(__file__).parent.parent.parent / 'tests' / 'DataSet' / 'GW100' / 'structures' / f'{cas_number}.xyz'

    if not xyz_file.exists():
        raise FileNotFoundError(
            f"GW100 structure not found: {xyz_file}\n"
            f"Expected CAS number format: XXXX-XX-X\n"
            f"Available structures should be in tests/DataSet/GW100/structures/"
        )

    with open(xyz_file, 'r') as f:
        lines = f.readlines()

    # Skip first 2 lines (atom count and comment)
    return ''.join(lines[2:])


# Convenience function for creating G₀W₀ driver from PySCF
def create_g0w0_driver_from_pyscf(mf, auxbasis: str = 'cc-pvdz-jkfit',
                                  config: Optional[G0W0Config] = None) -> G0W0Driver:
    """Create G0W0Driver from completed PySCF calculation.

    Args:
        mf: Completed PySCF mean-field calculation (RHF/RKS)
        auxbasis: Auxiliary basis set for RI/DF
        config: Optional G0W0Config (uses defaults if not provided)

    Returns:
        G0W0Driver: Initialized driver ready for calculation

    Example:
        >>> from pyscf import gto, scf
        >>> mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='cc-pvdz')
        >>> mf = scf.RHF(mol).run()
        >>> gw = create_g0w0_driver_from_pyscf(mf, auxbasis='cc-pvdz-jkfit')
        >>> result = gw.kernel()
    """
    if config is None:
        config = G0W0Config(auxbasis=auxbasis)

    return G0W0Driver(mf, config)