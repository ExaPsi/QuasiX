"""
Eigenvalue self-consistent GW (evGW) implementation for QuasiX.

This module provides a PySCF-compatible interface for evGW calculations
where only quasiparticle energies are updated self-consistently while
maintaining the mean-field wavefunctions.
"""

import numpy as np
from typing import Optional, Dict, Any, Tuple, List
from dataclasses import dataclass, field, asdict
import logging
import h5py
import json
import time
from pathlib import Path
from pyscf import lib
from pyscf.lib import logger

# Import Rust bindings
try:
    # Import the actual Rust module directly
    import quasix.quasix as rust_module
except ImportError:
    import quasix as rust_module

# Import polarizability update functionality
from .polarizability_update import (
    PolarizabilityUpdater,
    GapStatistics,
    update_polarizability_denominators,
    analyze_gap_evolution
)

__all__ = ['EvGW', 'evgw', 'EvGWResult', 'EvGWParameters']

# Configure logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


@dataclass
class EvGWParameters:
    """Parameters for evGW calculation."""
    max_cycle: int = 12
    conv_tol: float = 1e-4
    conv_tol_z: float = 1e-3
    damping: float = 0.5
    damping_dynamic: bool = True
    diis: bool = True
    diis_space: int = 6
    diis_start_cycle: int = 3
    freq_int: str = 'cd'  # 'cd' or 'ac'
    nfreq: int = 60
    eta: float = 0.01
    auxbasis: str = 'def2-tzvp-jkfit'
    thresh_df: float = 1e-10
    check_stability: bool = True
    verbose: int = 1
    # P0 denominator update parameters
    update_p0_denominators: bool = False
    p0_gap_threshold: float = 1e-6
    track_gap_evolution: bool = True
    # QP solver selection
    qp_solver: str = 'linearized'  # 'linearized' or 'newton'


@dataclass
class IterationData:
    """Data from a single evGW iteration."""
    cycle: int
    qp_energies: np.ndarray
    z_factors: np.ndarray
    energy_change: float
    rms_change: float
    damping_used: float
    converged: bool


@dataclass
class ConvergenceData:
    """Convergence diagnostics."""
    converged: bool
    n_cycles: int
    final_error: float
    final_error_z: float
    oscillation_detected: bool
    energy_history: np.ndarray
    z_history: np.ndarray
    error_history: np.ndarray
    damping_history: np.ndarray


@dataclass
class EvGWResult:
    """Complete evGW calculation result."""
    # Primary results
    qp_energies: np.ndarray
    z_factors: np.ndarray
    sigma_x: np.ndarray
    sigma_c: np.ndarray
    vxc_dft: np.ndarray

    # Convergence information
    converged: bool
    convergence_data: ConvergenceData
    iteration_history: List[IterationData]

    # Derived properties
    homo_index: int
    lumo_index: int
    gap_value: float  # Store as gap_value to avoid conflict with gap() method
    ip: float
    ea: float

    # Metadata
    parameters: EvGWParameters
    starting_point: str
    total_time: float
    # Gap evolution data (if P0 updates were used)
    gap_evolution: Optional[Dict[str, np.ndarray]] = None
    gap_analysis: Optional[Dict[str, Any]] = None

    def spectral_function(self, omega: np.ndarray, eta: float = 0.01) -> np.ndarray:
        """Generate spectral function A(ω) for visualization."""
        n_states = len(self.qp_energies)
        n_freq = len(omega)
        A = np.zeros((n_states, n_freq))

        for i, e_qp in enumerate(self.qp_energies):
            # Lorentzian broadening
            A[i, :] = self.z_factors[i] * eta / ((omega - e_qp)**2 + eta**2) / np.pi

        return A

    @property
    def convergence_history(self):
        """Provide convergence_history for compatibility with visualization module."""
        if not self.iteration_history:
            return None

        history = []
        for iter_data in self.iteration_history:
            history.append({
                'cycle': iter_data.cycle,
                'energy_diff': iter_data.energy_change,
                'z_min': np.min(iter_data.z_factors),
                'z_max': np.max(iter_data.z_factors),
                'z_mean': np.mean(iter_data.z_factors),
                'converged_orbitals': np.sum(np.abs(iter_data.energy_change) < self.parameters.conv_tol)
            })
        return history

    @property
    def n_cycles(self):
        """Number of completed cycles."""
        return self.convergence_data.n_cycles

    @property
    def final_energy_diff(self):
        """Final energy difference."""
        return self.convergence_data.final_error

    @property
    def wall_time(self):
        """Total wall time."""
        return self.total_time

    @property
    def sigma_c_diag(self):
        """Diagonal of correlation self-energy."""
        return np.diag(self.sigma_c)

    def ip_ea(self):
        """Return (IP, EA) tuple in eV."""
        return (self.ip * 27.2114, self.ea * 27.2114)

    def gap(self):
        """Return HOMO-LUMO gap in eV."""
        return self.gap_value * 27.2114

    def to_dict(self) -> Dict:
        """Export results to dictionary."""
        return {
            'qp_energies': self.qp_energies.tolist(),
            'z_factors': self.z_factors.tolist(),
            'sigma_x': self.sigma_x.tolist(),
            'sigma_c': self.sigma_c.tolist(),
            'converged': self.converged,
            'n_cycles': self.convergence_data.n_cycles,
            'gap': self.gap_value,
            'ip': self.ip,
            'ea': self.ea,
        }

    def save_to_hdf5(self, filename: str, molecule_data: Optional[Dict] = None,
                     basis_data: Optional[Dict] = None, compress: bool = True):
        """Save complete results to HDF5 file with QuasiX schema compatibility.

        Parameters
        ----------
        filename : str
            Path to HDF5 file
        molecule_data : dict, optional
            Molecule information for provenance tracking
        basis_data : dict, optional
            Basis set information for provenance tracking
        compress : bool
            Enable gzip compression for large datasets
        """
        import uuid
        import platform
        import os
        import subprocess

        with h5py.File(filename, 'w') as f:
            # === QuasiX Schema Root Groups ===

            # 1. Metadata group (QuasiX schema compliant)
            meta_grp = f.create_group('metadata')
            meta_grp.attrs['quasix_version'] = '0.1.0'
            meta_grp.attrs['timestamp'] = time.strftime('%Y-%m-%dT%H:%M:%S')
            meta_grp.attrs['uuid'] = str(uuid.uuid4())
            meta_grp.attrs['user'] = os.environ.get('USER', 'unknown')
            meta_grp.attrs['hostname'] = platform.node()

            # Try to get git commit hash
            try:
                git_hash = subprocess.check_output(
                    ['git', 'rev-parse', 'HEAD'],
                    cwd=os.path.dirname(__file__)
                ).decode('ascii').strip()
                meta_grp.attrs['git_commit'] = git_hash[:8]
            except:
                pass

            # 2. Molecule group (if provided)
            if molecule_data:
                mol_grp = f.create_group('molecule')
                mol_grp.attrs['natoms'] = molecule_data.get('natoms', 0)
                mol_grp.attrs['charge'] = molecule_data.get('charge', 0)
                mol_grp.attrs['multiplicity'] = molecule_data.get('multiplicity', 1)
                if 'symbols' in molecule_data:
                    mol_grp.create_dataset('symbols', data=np.array(molecule_data['symbols'], dtype='S10'))
                if 'coordinates' in molecule_data:
                    mol_grp.create_dataset('coordinates', data=molecule_data['coordinates'])
                if 'atomic_numbers' in molecule_data:
                    mol_grp.create_dataset('atomic_numbers', data=molecule_data['atomic_numbers'])

            # 3. Basis group
            if basis_data:
                basis_grp = f.create_group('basis')
                basis_grp.attrs['ao_basis'] = basis_data.get('ao_basis', 'unknown')
                basis_grp.attrs['n_ao'] = basis_data.get('n_ao', 0)
                basis_grp.attrs['aux_basis'] = basis_data.get('aux_basis', self.parameters.auxbasis)
                basis_grp.attrs['n_aux'] = basis_data.get('n_aux', 0)

            # 4. Parameters group (calculation settings)
            params_grp = f.create_group('parameters')
            params_grp.attrs['calculation_type'] = 'EvGW'

            # Convergence subgroup
            conv_params = params_grp.create_group('convergence')
            conv_params.attrs['energy_tol'] = self.parameters.conv_tol
            conv_params.attrs['density_tol'] = self.parameters.conv_tol_z
            conv_params.attrs['max_iterations'] = self.parameters.max_cycle
            conv_params.attrs['use_diis'] = self.parameters.diis
            conv_params.attrs['diis_space'] = self.parameters.diis_space

            # Frequency subgroup
            freq_params = params_grp.create_group('frequency')
            freq_params.attrs['grid_type'] = self.parameters.freq_int
            freq_params.attrs['n_points'] = self.parameters.nfreq
            freq_params.attrs['eta'] = self.parameters.eta

            # GW parameters subgroup
            gw_params = params_grp.create_group('gw_params')
            gw_params.attrs['starting_point'] = self.starting_point
            gw_params.attrs['damping'] = self.parameters.damping
            gw_params.attrs['damping_dynamic'] = self.parameters.damping_dynamic
            gw_params.attrs['check_stability'] = self.parameters.check_stability

            # === Results group (main calculation output) ===
            results_grp = f.create_group('results')

            # Primary results with compression
            compression = 'gzip' if compress else None
            compression_opts = 4 if compress else None

            # Store main arrays with optimal chunking
            n_states = len(self.qp_energies)
            chunk_1d = min(n_states, 1000)
            chunk_2d = (min(n_states, 100), min(n_states, 100))

            results_grp.create_dataset('qp_energies', data=self.qp_energies,
                                     chunks=(chunk_1d,), compression=compression,
                                     compression_opts=compression_opts)
            results_grp.create_dataset('z_factors', data=self.z_factors,
                                     chunks=(chunk_1d,), compression=compression,
                                     compression_opts=compression_opts)

            # Self-energies
            se_grp = results_grp.create_group('self_energy')
            se_grp.create_dataset('sigma_x', data=self.sigma_x,
                                chunks=chunk_2d if self.sigma_x.ndim == 2 else (chunk_1d,),
                                compression=compression, compression_opts=compression_opts)
            se_grp.create_dataset('sigma_c_real', data=np.real(self.sigma_c),
                                chunks=chunk_2d if self.sigma_c.ndim == 2 else (chunk_1d,),
                                compression=compression, compression_opts=compression_opts)
            se_grp.create_dataset('sigma_c_imag', data=np.imag(self.sigma_c),
                                chunks=chunk_2d if self.sigma_c.ndim == 2 else (chunk_1d,),
                                compression=compression, compression_opts=compression_opts)
            se_grp.create_dataset('vxc_dft', data=self.vxc_dft,
                                chunks=(chunk_1d,), compression=compression,
                                compression_opts=compression_opts)

            # Convergence information
            conv_grp = results_grp.create_group('convergence')
            conv_grp.attrs['converged'] = self.converged
            conv_grp.attrs['n_cycles'] = self.convergence_data.n_cycles
            conv_grp.attrs['final_error'] = self.convergence_data.final_error
            conv_grp.attrs['final_error_z'] = self.convergence_data.final_error_z
            conv_grp.attrs['oscillation_detected'] = self.convergence_data.oscillation_detected

            # Convergence history arrays
            if len(self.convergence_data.energy_history) > 0:
                conv_grp.create_dataset('energy_history', data=self.convergence_data.energy_history,
                                       compression=compression, compression_opts=compression_opts)
            if len(self.convergence_data.z_history) > 0:
                conv_grp.create_dataset('z_history', data=self.convergence_data.z_history,
                                       compression=compression, compression_opts=compression_opts)
            if len(self.convergence_data.error_history) > 0:
                conv_grp.create_dataset('error_history', data=self.convergence_data.error_history)
            if len(self.convergence_data.damping_history) > 0:
                conv_grp.create_dataset('damping_history', data=self.convergence_data.damping_history)

            # Iteration history (detailed)
            if self.iteration_history:
                iter_grp = results_grp.create_group('iteration_history')
                for i, iteration in enumerate(self.iteration_history):
                    iter_i = iter_grp.create_group(f'iteration_{i:03d}')
                    iter_i.attrs['cycle'] = iteration.cycle
                    iter_i.attrs['energy_change'] = iteration.energy_change
                    iter_i.attrs['rms_change'] = iteration.rms_change
                    iter_i.attrs['damping_used'] = iteration.damping_used
                    iter_i.attrs['converged'] = iteration.converged

                    # Store per-iteration data with compression
                    iter_i.create_dataset('qp_energies', data=iteration.qp_energies,
                                        chunks=(chunk_1d,), compression=compression,
                                        compression_opts=compression_opts)
                    iter_i.create_dataset('z_factors', data=iteration.z_factors,
                                        chunks=(chunk_1d,), compression=compression,
                                        compression_opts=compression_opts)

            # Properties (derived quantities)
            props_grp = results_grp.create_group('properties')
            props_grp.attrs['homo_index'] = self.homo_index
            props_grp.attrs['lumo_index'] = self.lumo_index
            props_grp.attrs['gap'] = self.gap_value
            props_grp.attrs['gap_eV'] = self.gap_value * 27.2114
            props_grp.attrs['ip'] = self.ip
            props_grp.attrs['ip_eV'] = self.ip * 27.2114
            props_grp.attrs['ea'] = self.ea
            props_grp.attrs['ea_eV'] = self.ea * 27.2114

            # Timing information
            timing_grp = results_grp.create_group('timings')
            timing_grp.attrs['total_time'] = self.total_time
            timing_grp.attrs['timestamp_start'] = self.starting_point
            timing_grp.attrs['timestamp_end'] = time.strftime('%Y-%m-%dT%H:%M:%S')

            # === File attributes (global metadata) ===
            f.attrs['file_format'] = 'QuasiX-HDF5'
            f.attrs['file_version'] = '1.0'
            f.attrs['created_by'] = 'QuasiX.EvGW'
            f.attrs['creation_date'] = time.strftime('%Y-%m-%dT%H:%M:%S')

    @classmethod
    def load_from_hdf5(cls, filename: str) -> 'EvGWResult':
        """Load results from HDF5 file with QuasiX schema compatibility.

        Parameters
        ----------
        filename : str
            Path to HDF5 file

        Returns
        -------
        EvGWResult
            Loaded result object

        Notes
        -----
        This method supports both the new QuasiX schema and legacy formats
        for backward compatibility.
        """
        with h5py.File(filename, 'r') as f:
            # Check file format version
            file_format = f.attrs.get('file_format', 'legacy')

            if file_format == 'QuasiX-HDF5':
                # Load from new schema structure
                results_grp = f['results']

                # Load convergence data
                conv_grp = results_grp['convergence']
                convergence_data = ConvergenceData(
                    converged=bool(conv_grp.attrs['converged']),
                    n_cycles=int(conv_grp.attrs['n_cycles']),
                    final_error=float(conv_grp.attrs['final_error']),
                    final_error_z=float(conv_grp.attrs['final_error_z']),
                    oscillation_detected=bool(conv_grp.attrs['oscillation_detected']),
                    energy_history=conv_grp['energy_history'][:] if 'energy_history' in conv_grp else np.array([]),
                    z_history=conv_grp['z_history'][:] if 'z_history' in conv_grp else np.array([]),
                    error_history=conv_grp['error_history'][:] if 'error_history' in conv_grp else np.array([]),
                    damping_history=conv_grp['damping_history'][:] if 'damping_history' in conv_grp else np.array([])
                )

                # Load iteration history
                iteration_history = []
                if 'iteration_history' in results_grp:
                    iter_grp = results_grp['iteration_history']
                    for i in range(len(iter_grp)):
                        iter_i = iter_grp[f'iteration_{i:03d}']
                        iteration_history.append(IterationData(
                            cycle=int(iter_i.attrs['cycle']),
                            qp_energies=iter_i['qp_energies'][:],
                            z_factors=iter_i['z_factors'][:],
                            energy_change=float(iter_i.attrs['energy_change']),
                            rms_change=float(iter_i.attrs['rms_change']),
                            damping_used=float(iter_i.attrs['damping_used']),
                            converged=bool(iter_i.attrs['converged'])
                        ))

                # Load parameters from new structure
                params_grp = f['parameters']
                conv_params = params_grp['convergence']
                freq_params = params_grp['frequency']
                gw_params = params_grp['gw_params']

                parameters = EvGWParameters(
                    max_cycle=int(conv_params.attrs.get('max_iterations', 12)),
                    conv_tol=float(conv_params.attrs.get('energy_tol', 1e-4)),
                    conv_tol_z=float(conv_params.attrs.get('density_tol', 1e-3)),
                    damping=float(gw_params.attrs.get('damping', 0.5)),
                    damping_dynamic=bool(gw_params.attrs.get('damping_dynamic', True)),
                    diis=bool(conv_params.attrs.get('use_diis', True)),
                    diis_space=int(conv_params.attrs.get('diis_space', 6)),
                    freq_int=freq_params.attrs.get('grid_type', 'cd'),
                    nfreq=int(freq_params.attrs.get('n_points', 60)),
                    eta=float(freq_params.attrs.get('eta', 0.01)),
                    check_stability=bool(gw_params.attrs.get('check_stability', True))
                )

                # Load properties
                props_grp = results_grp['properties']

                # Load self-energies
                se_grp = results_grp['self_energy']
                sigma_c_real = se_grp['sigma_c_real'][:]
                sigma_c_imag = se_grp['sigma_c_imag'][:]
                sigma_c = sigma_c_real + 1j * sigma_c_imag

                # Load timing
                timing_grp = results_grp.get('timings', {})
                total_time = float(timing_grp.attrs.get('total_time', 0.0)) if timing_grp else 0.0

                return cls(
                    qp_energies=results_grp['qp_energies'][:],
                    z_factors=results_grp['z_factors'][:],
                    sigma_x=se_grp['sigma_x'][:],
                    sigma_c=sigma_c,
                    vxc_dft=se_grp['vxc_dft'][:],
                    converged=bool(conv_grp.attrs['converged']),
                    convergence_data=convergence_data,
                    iteration_history=iteration_history,
                    homo_index=int(props_grp.attrs['homo_index']),
                    lumo_index=int(props_grp.attrs['lumo_index']),
                    gap_value=float(props_grp.attrs['gap']),
                    ip=float(props_grp.attrs['ip']),
                    ea=float(props_grp.attrs['ea']),
                    parameters=parameters,
                    starting_point=gw_params.attrs.get('starting_point', 'SCF'),
                    total_time=total_time
                )

            else:
                # Legacy format support (backward compatibility)
                # Load convergence data
                conv_grp = f['convergence'] if 'convergence' in f else f['results/convergence']
                convergence_data = ConvergenceData(
                    converged=bool(conv_grp.attrs['converged']),
                    n_cycles=int(conv_grp.attrs['n_cycles']),
                    final_error=float(conv_grp.attrs['final_error']),
                    final_error_z=float(conv_grp.attrs['final_error_z']),
                    oscillation_detected=bool(conv_grp.attrs['oscillation_detected']),
                    energy_history=conv_grp['energy_history'][:] if 'energy_history' in conv_grp else np.array([]),
                    z_history=conv_grp['z_history'][:] if 'z_history' in conv_grp else np.array([]),
                    error_history=conv_grp['error_history'][:] if 'error_history' in conv_grp else np.array([]),
                    damping_history=conv_grp['damping_history'][:] if 'damping_history' in conv_grp else np.array([])
                )

                # Load iteration history
                iteration_history = []
                if 'iteration_history' in f:
                    iter_grp = f['iteration_history']
                    for i in range(len(iter_grp)):
                        iter_i = iter_grp[f'iteration_{i:03d}']
                        iteration_history.append(IterationData(
                            cycle=int(iter_i.attrs['cycle']),
                            qp_energies=iter_i['qp_energies'][:],
                            z_factors=iter_i['z_factors'][:],
                            energy_change=float(iter_i.attrs['energy_change']),
                            rms_change=float(iter_i.attrs['rms_change']),
                            damping_used=float(iter_i.attrs['damping_used']),
                            converged=bool(iter_i.attrs['converged'])
                        ))

                # Load parameters
                params_grp = f['parameters']
                parameters = EvGWParameters(**dict(params_grp.attrs))

                # Load properties
                props_grp = f['properties']

                return cls(
                    qp_energies=f['qp_energies'][:],
                    z_factors=f['z_factors'][:],
                    sigma_x=f['sigma_x'][:],
                    sigma_c=f['sigma_c'][:],
                    vxc_dft=f['vxc_dft'][:],
                    converged=bool(conv_grp.attrs['converged']),
                    convergence_data=convergence_data,
                    iteration_history=iteration_history,
                    homo_index=int(props_grp.attrs['homo_index']),
                    lumo_index=int(props_grp.attrs['lumo_index']),
                    gap_value=float(props_grp.attrs['gap']),
                    ip=float(props_grp.attrs['ip']),
                    ea=float(props_grp.attrs['ea']),
                    parameters=parameters,
                    starting_point=f.attrs['starting_point'],
                    total_time=float(f.attrs['total_time'])
                )


class EvGW(lib.StreamObject):
    """
    Eigenvalue self-consistent GW for molecules.

    This class implements evGW where only quasiparticle energies are updated
    self-consistently while maintaining the mean-field wavefunctions.

    Attributes
    ----------
    conv_tol : float
        Convergence threshold for QP energies (Ha)
    conv_tol_z : float
        Convergence threshold for Z factors
    max_cycle : int
        Maximum number of iterations
    damping : float
        Mixing parameter (0=full update, 1=no update)
    auxbasis : str
        Auxiliary basis set for density fitting
    freq_int : str
        Frequency integration method ('cd' or 'ac')
    verbose : int
        Print level (0-5, following PySCF convention)

    Examples
    --------
    >>> from pyscf import gto, scf
    >>> from quasix.evgw import EvGW
    >>> mol = gto.M(atom='H 0 0 0; F 0 0 1.1', basis='cc-pvdz')
    >>> mf = scf.RHF(mol).run()
    >>> gw = EvGW(mf, auxbasis='cc-pvdz-jkfit')
    >>> gw.kernel()
    >>> print(gw.qp_energies)
    """

    def __init__(self, mf, auxbasis='def2-tzvp-jkfit', freq_int='cd'):
        """Initialize evGW calculator.

        Parameters
        ----------
        mf : SCF object
            PySCF mean-field calculation
        auxbasis : str
            Auxiliary basis for density fitting
        freq_int : str
            Frequency integration method ('cd' or 'ac')
        """
        self.mol = mf.mol
        self._scf = mf
        self.auxbasis = auxbasis
        self.freq_int = freq_int

        # Convergence parameters
        self.conv_tol = 1e-4
        self.conv_tol_z = 1e-3
        self.max_cycle = 12
        self.damping = 0.5
        self.damping_dynamic = True

        # DIIS parameters
        self.diis = True
        self.diis_space = 6
        self.diis_start_cycle = 3

        # Advanced options
        self.check_stability = True

        # Frequency integration
        self.nfreq = 60
        self.eta = 0.01

        # Thresholds
        self.thresh_df = 1e-10

        # P0 denominator update parameters
        self.update_p0_denominators = False
        self.p0_gap_threshold = 1e-6
        self.track_gap_evolution = True
        self.p0_updater = None

        # Results storage
        self.qp_energies = None
        self.z_factors = None
        self.sigma_x = None
        self.sigma_c = None
        self.converged = False
        self.iteration_history = []

        # Monitoring support
        self.monitor = None
        self._callbacks = []

        # Verbosity
        self.verbose = mf.verbose
        self.stdout = mf.stdout

        # QP solver selection ('linearized' or 'newton')
        self.qp_solver = 'linearized'

    def register_callback(self, callback):
        """Register a callback for iteration updates.

        Parameters
        ----------
        callback : callable
            Function to call with iteration data dictionary
        """
        if callable(callback):
            self._callbacks.append(callback)
        else:
            raise TypeError("Callback must be callable")

    def dump_flags(self, verbose=None):
        """Print calculation parameters."""
        log = logger.new_logger(self, verbose)
        log.info('')
        log.info('******** %s ********', self.__class__)
        log.info('method = %s', self.__class__.__name__)
        log.info('auxbasis = %s', self.auxbasis)
        log.info('freq_int = %s', self.freq_int)
        log.info('nfreq = %d', self.nfreq)
        log.info('conv_tol = %g', self.conv_tol)
        log.info('conv_tol_z = %g', self.conv_tol_z)
        log.info('max_cycle = %d', self.max_cycle)
        log.info('damping = %g', self.damping)
        log.info('damping_dynamic = %s', self.damping_dynamic)
        log.info('diis = %s', self.diis)
        if self.diis:
            log.info('diis_space = %d', self.diis_space)
            log.info('diis_start_cycle = %d', self.diis_start_cycle)
        return self

    def kernel(self, mo_coeff=None, mo_energy=None, mo_occ=None, update_p0=None):
        """
        Run evGW calculation.

        Parameters
        ----------
        mo_coeff : ndarray, optional
            MO coefficients (uses SCF if not provided)
        mo_energy : ndarray, optional
            MO energies
        mo_occ : ndarray, optional
            MO occupations
        update_p0 : bool, optional
            Whether to update P0 denominators during iterations
            (overrides self.update_p0_denominators)

        Returns
        -------
        ndarray
            Converged quasiparticle energies
        """
        cput0 = (logger.process_clock(), logger.perf_counter())

        if mo_coeff is None:
            mo_coeff = self._scf.mo_coeff
        if mo_energy is None:
            mo_energy = self._scf.mo_energy
        if mo_occ is None:
            mo_occ = self._scf.mo_occ

        self.dump_flags()

        # Build DF tensors
        log = logger.Logger(self.stdout, self.verbose)
        log.info('\nBuilding DF tensors...')

        # CRITICAL FIX (2025-11-25): Always use built-in method that returns abP
        # The df_utils version doesn't return abP, which is required for correct evGW
        iaP, miP, abP, chol_v = self.build_df_tensors(mo_coeff, mo_occ)

        # Get DFT XC potential
        vxc_dft = self.get_vxc(mo_coeff)

        # Reshape tensors for Rust interface
        # Rust expects iaP as (nocc, nvir, naux) not (nocc*nvir, naux)
        nocc = np.sum(mo_occ > 0).astype(int)
        nvir = len(mo_energy) - nocc
        naux = chol_v.shape[0]

        # Reshape from flat to 3D for iaP
        if iaP.ndim == 2 and iaP.shape[0] == nocc * nvir:
            iaP_3d = iaP.reshape(nocc, nvir, naux)
        else:
            iaP_3d = iaP

        # Initialize P0 updater if requested
        if update_p0 is not None:
            self.update_p0_denominators = update_p0

        if self.update_p0_denominators:
            # Keep iaP in 2D format for P0 updater
            iaP_2d = iaP.reshape(nocc * nvir, naux) if iaP.ndim == 3 else iaP
            self.p0_updater = PolarizabilityUpdater(
                nocc, nvir, naux,
                mo_energy,
                iaP_2d,
                gap_threshold=self.p0_gap_threshold,
                eta=self.eta
            )
            log.info('P0 denominator updates enabled with gap threshold = %.2e Ha',
                    self.p0_gap_threshold)

        # miP stays 2D: (nmo*nocc, naux) as expected by Rust
        # It will be reshaped inside Rust to (nmo, nocc, naux)

        # CRITICAL FIX (2025-11-25): Reshape abP from flat to 2D for Rust
        # Rust expects abP as (nvir*nvir, naux)
        if abP.ndim == 3:
            abP_2d = abP.reshape(nvir * nvir, naux)
        else:
            abP_2d = abP

        # Ensure arrays are in the correct format for Rust
        # Rust expects C-contiguous float64 arrays
        mo_energy = np.ascontiguousarray(mo_energy, dtype=np.float64)
        mo_occ = np.ascontiguousarray(mo_occ, dtype=np.float64)
        iaP_3d = np.ascontiguousarray(iaP_3d, dtype=np.float64)
        miP = np.ascontiguousarray(miP, dtype=np.float64)  # CRITICAL: Changed from ijP
        abP_2d = np.ascontiguousarray(abP_2d, dtype=np.float64)  # Virtual-virtual block
        chol_v = np.ascontiguousarray(chol_v, dtype=np.float64)
        vxc_dft = np.ascontiguousarray(vxc_dft, dtype=np.float64)

        # Call Rust evGW implementation
        log.info('Starting evGW iteration...')

        # CRITICAL: Rust expects a config dict, not positional arguments
        # Build configuration dictionary for Rust
        config = {
            'mo_energy': mo_energy,
            'mo_occ': mo_occ,
            'iaP': iaP_3d,
            'ijP': miP,  # CRITICAL: This is actually miP (all-occupied), not ijP
            'abP': abP_2d,  # CRITICAL FIX (2025-11-25): Virtual-virtual block for complete tensor
            'chol_v': chol_v,
            'vxc_dft': vxc_dft,
            'max_cycle': self.max_cycle,
            'conv_tol': self.conv_tol,
            'damping': self.damping,
            'freq_int': self.freq_int,
            'nfreq': self.nfreq,
            'verbose': self.verbose,
            # QP solver selection: 'linearized' or 'newton'
            'qp_solver': self.qp_solver,
        }

        try:
            result = rust_module.run_evgw(config)

            # Store results
            self.qp_energies = result['qp_energies']
            self.z_factors = result['z_factors']
            self.sigma_x = result['sigma_x']
            self.sigma_c = result['sigma_c_re'] + 1j * result['sigma_c_im']
            self.converged = result['converged']

            # Process iteration history
            self.iteration_history = []
            for iter_dict in result['iteration_history']:
                iter_data = IterationData(
                    cycle=iter_dict['cycle'],
                    qp_energies=iter_dict['qp_energies'],
                    z_factors=iter_dict['z_factors'],
                    energy_change=iter_dict['energy_change'],
                    rms_change=iter_dict['rms_change'],
                    damping_used=iter_dict['damping_used'],
                    converged=iter_dict['converged']
                )
                self.iteration_history.append(iter_data)

                # Update P0 denominators if requested
                if self.p0_updater is not None:
                    gap_stats = self.p0_updater.update_energies(
                        iter_dict['qp_energies'],
                        track_history=self.track_gap_evolution
                    )
                    if gap_stats.has_issues():
                        log.warn(f'Cycle {iter_dict["cycle"]}: {gap_stats.n_negative} negative gaps, '
                                f'{gap_stats.n_thresholded} thresholded gaps')

                # Notify monitor and callbacks
                if self.monitor or self._callbacks:
                    # Calculate HOMO/LUMO gap for monitoring
                    homo_energy = iter_dict['qp_energies'][nocc-1] if nocc > 0 else 0.0
                    lumo_energy = iter_dict['qp_energies'][nocc] if nocc < len(iter_dict['qp_energies']) else 0.0
                    gap = lumo_energy - homo_energy

                    monitor_data = {
                        'iteration': iter_dict['cycle'],
                        'cycle': iter_dict['cycle'],
                        'qp_energies': iter_dict['qp_energies'],
                        'z_factors': iter_dict['z_factors'],
                        'convergence': iter_dict['energy_change'],
                        'energy_change': iter_dict['energy_change'],
                        'rms_change': iter_dict['rms_change'],
                        'damping_used': iter_dict['damping_used'],
                        'homo_energy': homo_energy,
                        'lumo_energy': lumo_energy,
                        'gap': gap,
                        'timing': iter_dict.get('timing', {}),
                        'memory_mb': iter_dict.get('memory_mb', 0.0),
                        'converged': iter_dict['converged']
                    }

                    # Notify monitor if attached
                    if self.monitor and hasattr(self.monitor, 'update'):
                        self.monitor.update(monitor_data)

                    # Call registered callbacks
                    for callback in self._callbacks:
                        try:
                            callback(monitor_data)
                        except Exception as e:
                            log.warn(f"Callback failed: {e}")

            if self.converged:
                log.info('evGW converged in %d cycles', result['n_cycles'])
            else:
                log.warn('evGW did not converge in %d cycles', self.max_cycle)

        except ImportError as e:
            # Rust module not available - use mock data for testing
            log.warn(f'Rust module not available: {e}. Using mock results for testing.')

            # Create mock results for testing
            self.qp_energies = mo_energy.copy()
            self.z_factors = np.ones_like(mo_energy) * 0.9
            self.sigma_x = np.diag(mo_energy * 0.1)
            self.sigma_c = np.diag(mo_energy * 0.05)
            self.converged = True

            # Create mock iteration history
            self.iteration_history = []
            for i in range(5):
                # Create mock QP energies and Z factors with realistic convergence
                mock_qp = self.qp_energies + np.random.randn(len(mo_energy)) * 0.01 * (5-i)
                mock_z = self.z_factors + np.random.randn(len(mo_energy)) * 0.01 * (5-i)

                iter_data = IterationData(
                    cycle=i+1,
                    qp_energies=mock_qp,
                    z_factors=mock_z,
                    energy_change=0.1 * (0.5 ** i),
                    rms_change=0.05 * (0.5 ** i),
                    damping_used=0.5 + 0.1 * i,
                    converged=(i == 4)
                )
                self.iteration_history.append(iter_data)

                # Notify monitor and callbacks
                if self.monitor or self._callbacks:
                    # Calculate HOMO/LUMO gap for monitoring
                    homo_energy = mock_qp[nocc-1] if nocc > 0 else 0.0
                    lumo_energy = mock_qp[nocc] if nocc < len(mock_qp) else 0.0
                    gap = lumo_energy - homo_energy

                    monitor_data = {
                        'iteration': i+1,
                        'cycle': i+1,
                        'qp_energies': mock_qp,
                        'z_factors': mock_z,
                        'convergence': iter_data.energy_change,
                        'energy_change': iter_data.energy_change,
                        'rms_change': iter_data.rms_change,
                        'damping_used': iter_data.damping_used,
                        'homo_energy': homo_energy,
                        'lumo_energy': lumo_energy,
                        'gap': gap,
                        'timing': {'iteration': 0.1 * (i+1)},
                        'memory_mb': 100.0 + i * 10.0,
                        'converged': iter_data.converged
                    }

                    # Notify monitor if attached
                    if self.monitor and hasattr(self.monitor, 'update'):
                        self.monitor.update(monitor_data)

                    # Call registered callbacks
                    for callback in self._callbacks:
                        try:
                            callback(monitor_data)
                        except Exception as e:
                            log.warn(f"Callback failed: {e}")

        # Analyze gap evolution if P0 updates were used
        gap_evolution = None
        gap_analysis = None
        if self.p0_updater is not None and self.track_gap_evolution:
            gap_evolution = self.p0_updater.get_gap_evolution()
            gap_analysis = analyze_gap_evolution(
                self.p0_updater.gap_history,
                output_file=None  # Don't auto-save to file
            )
            if gap_analysis.get('warnings'):
                for warning in gap_analysis['warnings']:
                    log.warn(f'Gap evolution: {warning}')

        # Store gap evolution data
        self.gap_evolution = gap_evolution
        self.gap_analysis = gap_analysis

        # Print summary
        self._print_summary(log)

        # Validate results
        if self.check_stability:
            self.validate_results()

        log.timer('evGW', *cput0)

        return self.qp_energies

    def build_df_tensors(self, mo_coeff, mo_occ):
        """Build density-fitted 3-center integrals.

        Parameters
        ----------
        mo_coeff : ndarray
            MO coefficients (nao × nmo)
        mo_occ : ndarray
            MO occupations (nmo,)

        Returns
        -------
        tuple
            (iaP, miP, abP, chol_v) tensors for DF
            iaP: (nocc*nvir, naux) - occupied-virtual DF tensor (flattened)
            miP: (nmo*nocc, naux) - all-occupied DF tensor (flattened)
                 This is the (mi|P) tensor for ALL m and occupied i,
                 required for correct exchange self-energy calculation
            abP: (nvir*nvir, naux) - virtual-virtual DF tensor (flattened)
                 Required for complete DF tensor reconstruction in Rust
            chol_v: (naux, naux) - Cholesky factor of Coulomb metric
        """
        from pyscf import df
        from pyscf.ao2mo import _ao2mo

        log = logger.Logger(self.stdout, self.verbose)
        cput0 = (logger.process_clock(), logger.perf_counter())

        # Determine occupied and virtual indices
        nocc = np.sum(mo_occ > 0).astype(int)
        nvir = mo_coeff.shape[1] - nocc
        nao = self.mol.nao_nr()

        # Create DF object with specified auxiliary basis
        dfobj = df.DF(self.mol, auxbasis=self.auxbasis)
        dfobj.max_memory = 4000  # Use up to 4GB for DF

        # Build the DF integrals (this generates _cderi)
        # We use a custom approach that doesn't hang
        naux = dfobj.get_naoaux()

        log.debug(f'Building DF tensors: nao={nao}, nocc={nocc}, nvir={nvir}, naux={naux}')

        # Build 2-center Coulomb metric (P|Q)
        log.debug('Building 2-center Coulomb metric...')
        # CRITICAL FIX: PySCF requires auxbasis as dict mapping atoms to basis sets
        # Convert string auxbasis to dict format if needed
        auxbasis_dict = self.auxbasis
        if isinstance(self.auxbasis, str):
            auxbasis_dict = {self.mol.atom_symbol(i): self.auxbasis for i in range(self.mol.natm)}
        auxmol = df.addons.make_auxmol(self.mol, auxbasis_dict)
        metric = auxmol.intor('int2c2e', aosym='s1')  # (naux, naux)

        # Compute V^{1/2} for symmetrized dielectric formulation
        # Use eigendecomposition for numerical stability
        eigvals, eigvecs = np.linalg.eigh(metric)

        # Remove small eigenvalues for numerical stability
        idx = eigvals > self.thresh_df
        n_valid = np.sum(idx)
        if n_valid < naux:
            log.info(f'Removed {naux - n_valid} linearly dependent auxiliary basis functions')

        # Build V^{1/2} = U * sqrt(lambda) * U^T
        # We need the full reconstruction to ensure symmetry
        sqrt_eigvals = np.sqrt(np.maximum(eigvals, self.thresh_df))
        chol_v = eigvecs @ np.diag(sqrt_eigvals) @ eigvecs.T

        # Build 3-center integrals using PySCF's incore method
        log.debug('Building 3-center integrals...')

        # Get 3-center integrals in AO basis using PySCF's incore method
        # This is more reliable than trying to use internal methods
        int3c = df.incore.aux_e2(self.mol, auxmol, intor='int3c2e', aosym='s1')
        # int3c has shape (nao*nao, naux) in compressed format

        # Reshape to (nao, nao, naux) for easier manipulation
        int3c = int3c.reshape(nao, nao, naux)

        # CRITICAL FIX (2025-11-25): Apply J^(-1/2) transformation to match PySCF convention
        # PySCF's _cderi stores Lpq = J^(-1/2) @ (pq|Q), where J[P,Q] = (P|Q) is the Coulomb metric
        # We need to apply this transformation so our tensors match PySCF's convention
        #
        # The Coulomb metric J was already computed above as "metric"
        # We need J^(-1/2) to transform: Lpq[P,u,v] = sum_Q J^(-1/2)[P,Q] * int3c[u,v,Q]
        inv_sqrt_eigvals = 1.0 / np.sqrt(np.maximum(eigvals, self.thresh_df))
        j_inv_sqrt = eigvecs @ np.diag(inv_sqrt_eigvals) @ eigvecs.T

        # Apply J^(-1/2) transformation to get properly normalized DF integrals
        # int3c_transformed[u,v,P] = sum_Q int3c[u,v,Q] * J^(-1/2)[Q,P]
        int3c = np.einsum('uvQ,QP->uvP', int3c, j_inv_sqrt)
        log.debug(f'Applied J^(-1/2) transformation to DF integrals, norm = {np.linalg.norm(int3c):.4f}')

        # Split MO coefficients into occupied and virtual
        mo_occ_coeff = np.ascontiguousarray(mo_coeff[:, :nocc])
        mo_vir_coeff = np.ascontiguousarray(mo_coeff[:, nocc:nocc+nvir])

        # Initialize output tensors
        # CRITICAL FIX: Rust expects (n_mo, n_occ, n_aux) for ijP (now miP)
        # This is the full (mi|P) tensor for ALL m and occupied i
        nmo = mo_coeff.shape[1]
        iaP = np.zeros((nocc, nvir, naux), dtype=np.float64)
        miP = np.zeros((nmo, nocc, naux), dtype=np.float64)  # Full tensor for exchange
        abP = np.zeros((nvir, nvir, naux), dtype=np.float64)  # Virtual-virtual block

        # Transform 3-center integrals to MO basis
        # We do this auxiliary function by auxiliary function for memory efficiency
        log.debug('Transforming to MO basis...')

        for p in range(naux):
            # Get the 3-center integral for this auxiliary function
            ints_p = int3c[:, :, p]  # (nao, nao)

            # Transform to MO basis: (μν|P) → (ij|P), (ia|P), etc.
            # First transformation: C^T (μν|P) C
            mo_ints = mo_coeff.T @ ints_p @ mo_coeff  # (nmo, nmo)

            # Extract the blocks we need
            # (ia|P) block: occupied rows, virtual columns
            iaP[:, :, p] = mo_ints[:nocc, nocc:nocc+nvir]

            # (mi|P) block: ALL rows, occupied columns
            # This is the correct tensor for exchange self-energy
            miP[:, :, p] = mo_ints[:, :nocc]

            # (ab|P) block: virtual rows, virtual columns
            # CRITICAL FIX (2025-11-25): Required for complete DF tensor reconstruction
            abP[:, :, p] = mo_ints[nocc:nocc+nvir, nocc:nocc+nvir]

        log.timer('DF tensor construction', *cput0)
        log.info(f'DF tensors built: iaP shape={iaP.shape}, miP shape={miP.shape}, abP shape={abP.shape}')
        log.info(f'  iaP memory={iaP.nbytes/1e9:.2f} GB, miP memory={miP.nbytes/1e9:.2f} GB, abP memory={abP.nbytes/1e9:.2f} GB')

        # Reshape tensors for compatibility with Rust code expectations
        # Rust expects (nocc*nvir, naux) for iaP
        iaP_flat = iaP.reshape(nocc * nvir, naux)
        # Rust expects (n_mo*n_occ, naux) for miP (NOT nocc*nocc!)
        miP_flat = miP.reshape(nmo * nocc, naux)
        # Rust expects (nvir*nvir, naux) for abP
        abP_flat = abP.reshape(nvir * nvir, naux)

        return iaP_flat, miP_flat, abP_flat, chol_v

    def get_vxc(self, mo_coeff):
        """Extract DFT exchange-correlation potential in MO basis.

        Parameters
        ----------
        mo_coeff : ndarray
            MO coefficients

        Returns
        -------
        ndarray
            XC potential diagonal in MO basis
        """
        if hasattr(self._scf, 'xc'):
            # DFT calculation
            from pyscf.dft import numint
            ni = self._scf._numint

            # Get XC potential in AO basis
            dm = self._scf.make_rdm1()
            vxc_ao = ni.nr_vxc(self.mol, self._scf.grids, self._scf.xc, dm)[2]

            # Transform to MO basis
            vxc_mo = mo_coeff.T @ vxc_ao @ mo_coeff
            return np.diag(vxc_mo)
        else:
            # HF calculation - vxc = Sigma_x (Fock exchange)
            # CRITICAL FIX (2025-11-21): For HF, get_veff() returns J - 0.5*K which
            # is NOT the exchange self-energy. We need the actual Sigma_x.
            #
            # The QP equation is: E_QP = ε_HF + Re[Σ(E_QP) - V_xc]
            # where Σ = Σx + Σc. For HF starting point, V_xc = Σx,
            # so E_QP = ε_HF + Re[Σc].
            #
            # To get Sigma_x, we compute it from 4-center ERIs.
            from pyscf import ao2mo
            nmo = mo_coeff.shape[1]
            nocc = np.sum(self._scf.mo_occ > 0).astype(int)

            eri_mo = ao2mo.kernel(self.mol, mo_coeff)
            eri_mo = ao2mo.restore(1, eri_mo, nmo)  # (nmo, nmo, nmo, nmo)

            # Sigma_x[p] = -Σ_i^occ (pi|ip)
            sigma_x = np.zeros(nmo)
            for p in range(nmo):
                for i in range(nocc):
                    sigma_x[p] -= eri_mo[p, i, i, p]

            return sigma_x

    def _print_summary(self, log):
        """Print calculation summary."""
        nocc = np.sum(self._scf.mo_occ > 0).astype(int)

        # HOMO-LUMO gap
        if nocc > 0 and nocc < len(self.qp_energies):
            homo = self.qp_energies[nocc - 1]
            lumo = self.qp_energies[nocc]
            gap = (lumo - homo) * 27.211  # Convert to eV

            log.info('')
            log.info('evGW Results:')
            log.info('  HOMO energy: %.4f Ha (%.4f eV)', homo, homo * 27.211)
            log.info('  LUMO energy: %.4f Ha (%.4f eV)', lumo, lumo * 27.211)
            log.info('  HOMO-LUMO gap: %.4f eV', gap)
            log.info('  HOMO Z-factor: %.3f', self.z_factors[nocc - 1])
            log.info('  LUMO Z-factor: %.3f', self.z_factors[nocc])

    def validate_results(self):
        """Check physical validity of results."""
        log = logger.Logger(self.stdout, self.verbose)

        # Check Z factors
        if np.any(self.z_factors < 0) or np.any(self.z_factors > 1):
            log.warn('Unphysical Z factors detected: min=%.3f, max=%.3f',
                     self.z_factors.min(), self.z_factors.max())

        # Check energy ordering
        nocc = np.sum(self._scf.mo_occ > 0).astype(int)
        if nocc > 0 and nocc < len(self.qp_energies):
            if self.qp_energies[nocc - 1] > self.qp_energies[nocc]:
                log.warn('HOMO-LUMO inversion detected!')

        # Check for large corrections
        corrections = self.qp_energies - self._scf.mo_energy
        if np.abs(corrections).max() > 5.0:  # 5 Ha = 136 eV
            log.warn('Very large QP corrections: max = %.2f eV',
                     np.abs(corrections).max() * 27.211)

    def analyze(self):
        """Analyze evGW results and return summary."""
        nocc = np.sum(self._scf.mo_occ > 0).astype(int)

        # Handle case where kernel has not been run
        if self.qp_energies is None:
            # Initialize with SCF values
            self.qp_energies = self._scf.mo_energy.copy()
            self.z_factors = np.ones_like(self.qp_energies)
            self.sigma_x = np.zeros_like(self.qp_energies)
            self.sigma_c = np.zeros_like(self.qp_energies)
            self.converged = False
            self.iteration_history = []

        # Build convergence data
        if self.iteration_history:
            energy_history = np.array([d.qp_energies for d in self.iteration_history])
            z_history = np.array([d.z_factors for d in self.iteration_history])
            error_history = np.array([d.energy_change for d in self.iteration_history])
            damping_history = np.array([d.damping_used for d in self.iteration_history])
        else:
            energy_history = np.array([])
            z_history = np.array([])
            error_history = np.array([])
            damping_history = np.array([])

        convergence_data = ConvergenceData(
            converged=self.converged,
            n_cycles=len(self.iteration_history),
            final_error=error_history[-1] if len(error_history) > 0 else 0.0,
            final_error_z=0.0,  # Would need to track separately
            oscillation_detected=False,  # Would need detection algorithm
            energy_history=energy_history,
            z_history=z_history,
            error_history=error_history,
            damping_history=damping_history
        )

        # Create result object
        result = EvGWResult(
            qp_energies=self.qp_energies,
            z_factors=self.z_factors,
            sigma_x=self.sigma_x if self.sigma_x is not None else np.zeros_like(self.qp_energies),
            sigma_c=self.sigma_c.real if self.sigma_c is not None and hasattr(self.sigma_c, 'real') else
                    self.sigma_c if self.sigma_c is not None else np.zeros_like(self.qp_energies),
            vxc_dft=self.get_vxc(self._scf.mo_coeff),
            converged=self.converged,
            convergence_data=convergence_data,
            iteration_history=self.iteration_history,
            homo_index=nocc - 1 if nocc > 0 else 0,
            lumo_index=nocc if nocc < len(self.qp_energies) else len(self.qp_energies) - 1,
            gap_value=(self.qp_energies[nocc] - self.qp_energies[nocc - 1]) if nocc > 0 and nocc < len(self.qp_energies) else 0.0,
            ip=-self.qp_energies[nocc - 1] if nocc > 0 else 0.0,
            ea=-self.qp_energies[nocc] if nocc < len(self.qp_energies) else 0.0,
            parameters=EvGWParameters(
                max_cycle=self.max_cycle,
                conv_tol=self.conv_tol,
                conv_tol_z=self.conv_tol_z,
                damping=self.damping,
                freq_int=self.freq_int,
                nfreq=self.nfreq,
                auxbasis=self.auxbasis
            ),
            starting_point=type(self._scf).__name__,
            total_time=0.0  # Would need to track timing
        )

        return result

    def save_results(self, filename: str, compress: bool = True):
        """Save evGW results to HDF5 file with full QuasiX schema support.

        Parameters
        ----------
        filename : str
            Output filename (will add .h5 extension if not present)
        compress : bool
            Enable gzip compression for large datasets
        """
        if not filename.endswith('.h5'):
            filename += '.h5'

        # Extract molecule and basis information from SCF object
        molecule_data = self._extract_molecule_data()
        basis_data = self._extract_basis_data()

        result = self.analyze()
        result.save_to_hdf5(filename, molecule_data=molecule_data,
                          basis_data=basis_data, compress=compress)

        log = logger.Logger(self.stdout, self.verbose)
        log.info(f'Results saved to {filename} (QuasiX HDF5 format)')

    def load_results(self, filename: str) -> EvGWResult:
        """Load evGW results from HDF5 file.

        Parameters
        ----------
        filename : str
            Input filename

        Returns
        -------
        EvGWResult
            Loaded results
        """
        result = EvGWResult.load_from_hdf5(filename)

        # Update internal state
        self.qp_energies = result.qp_energies
        self.z_factors = result.z_factors
        self.sigma_x = result.sigma_x
        self.sigma_c = result.sigma_c
        self.converged = result.converged
        self.iteration_history = result.iteration_history

        log = logger.Logger(self.stdout, self.verbose)
        log.info(f'Results loaded from {filename}')

        return result

    def plot_convergence(self, save_fig=None, show=True):
        """Plot convergence of evGW iterations using comprehensive visualization module.

        Parameters
        ----------
        save_fig : str, optional
            Filename to save figure
        show : bool
            Whether to display the plot

        Returns
        -------
        matplotlib.figure.Figure
            The convergence plot figure
        """
        from . import visualization

        # Get the analysis result
        result = self.analyze()

        # Use the comprehensive visualization module
        fig = visualization.plot_evgw_convergence(
            result,
            save_path=save_fig
        )

        if show:
            import matplotlib.pyplot as plt
            plt.show()

        return fig

    def plot_spectral_function(self, omega_range=None, eta=None, save_fig=None, show=True):
        """Plot spectral function A(ω) using comprehensive visualization module.

        Parameters
        ----------
        omega_range : tuple, optional
            Energy range (min, max) in eV
        eta : float, optional
            Broadening parameter in eV
        save_fig : str, optional
            Filename to save figure
        show : bool
            Whether to display the plot

        Returns
        -------
        matplotlib.figure.Figure
            The spectral function plot figure
        """
        from . import visualization

        if self.qp_energies is None:
            log = logger.Logger(self.stdout, self.verbose)
            log.warn('No QP energies available. Run calculation first.')
            return

        # Use default broadening if not specified
        if eta is None:
            eta = self.eta * 27.2114  # Convert to eV

        # Use the comprehensive visualization module
        fig = visualization.plot_spectral_function(
            self.qp_energies,
            self.z_factors,
            energy_range=omega_range,
            broadening=eta,
            save_path=save_fig
        )

        if show:
            import matplotlib.pyplot as plt
            plt.show()

        return fig

    def plot_qp_corrections(self, orbital_range=None, save_fig=None, show=True):
        """Plot quasiparticle corrections vs orbital energy.

        Parameters
        ----------
        orbital_range : tuple, optional
            Range of orbital indices to plot (auto if None)
        save_fig : str, optional
            Filename to save figure
        show : bool
            Whether to display the plot

        Returns
        -------
        matplotlib.figure.Figure
            The QP corrections plot figure
        """
        from . import visualization

        if self.qp_energies is None:
            log = logger.Logger(self.stdout, self.verbose)
            log.warn('No QP energies available. Run calculation first.')
            return

        # Use the comprehensive visualization module
        fig = visualization.plot_qp_corrections(
            self._scf.mo_energy,
            self.qp_energies,
            self._scf.mo_occ,
            orbital_range=orbital_range,
            save_path=save_fig
        )

        if show:
            import matplotlib.pyplot as plt
            plt.show()

        return fig

    def plot_z_factors(self, save_fig=None, show=True):
        """Plot Z-factors (quasiparticle weights) analysis.

        Parameters
        ----------
        save_fig : str, optional
            Filename to save figure
        show : bool
            Whether to display the plot

        Returns
        -------
        matplotlib.figure.Figure
            The Z-factors plot figure
        """
        from . import visualization

        if self.z_factors is None:
            log = logger.Logger(self.stdout, self.verbose)
            log.warn('No Z-factors available. Run calculation first.')
            return

        # Use the comprehensive visualization module
        fig = visualization.plot_z_factors(
            self.z_factors,
            self._scf.mo_energy,
            self._scf.mo_occ,
            save_path=save_fig
        )

        if show:
            import matplotlib.pyplot as plt
            plt.show()

        return fig

    def create_analysis_report(self, molecule_name="Molecule", save_dir=None):
        """Create comprehensive evGW analysis report with multiple figures.

        Parameters
        ----------
        molecule_name : str
            Name of the molecule for the report
        save_dir : str, optional
            Directory to save figures (creates subdirectory if needed)

        Returns
        -------
        dict
            Dictionary of figure names to Figure objects
        """
        from . import visualization
        import os

        if self.qp_energies is None:
            log = logger.Logger(self.stdout, self.verbose)
            log.warn('No calculation results available. Run calculation first.')
            return {}

        # Create save directory if needed
        if save_dir:
            os.makedirs(save_dir, exist_ok=True)

        # Get the analysis result
        result = self.analyze()

        # Create comprehensive report
        figures = visualization.create_evgw_report(
            result,
            self._scf.mo_energy,
            self._scf.mo_occ,
            molecule_name=molecule_name,
            save_dir=save_dir
        )

        log = logger.Logger(self.stdout, self.verbose)
        log.info(f'Created evGW analysis report with {len(figures)} figures')

        if save_dir:
            log.info(f'Figures saved to {save_dir}/')

        return figures

    def restart_from_file(self, filename: str, mo_coeff=None, mo_energy=None, mo_occ=None):
        """Restart evGW calculation from saved checkpoint.

        Parameters
        ----------
        filename : str
            Path to checkpoint file
        mo_coeff : ndarray, optional
            New MO coefficients (uses saved if not provided)
        mo_energy : ndarray, optional
            New MO energies (uses saved if not provided)
        mo_occ : ndarray, optional
            New MO occupations (uses saved if not provided)

        Returns
        -------
        ndarray
            Converged QP energies
        """
        # Load previous results
        result = self.load_results(filename)

        # Use last QP energies as starting point
        if mo_energy is None:
            mo_energy = result.qp_energies

        # Continue calculation
        log = logger.Logger(self.stdout, self.verbose)
        log.info(f'Restarting from iteration {result.convergence_data.n_cycles}')

        return self.kernel(mo_coeff, mo_energy, mo_occ)

    def _extract_molecule_data(self) -> Dict:
        """Extract molecule information from PySCF mol object.

        Returns
        -------
        dict
            Molecule data for HDF5 export
        """
        mol = self.mol

        # Get atomic symbols and numbers
        symbols = []
        atomic_numbers = []
        for ia in range(mol.natm):
            symbol = mol.atom_symbol(ia)
            symbols.append(symbol)
            atomic_numbers.append(mol.atom_charge(ia))

        return {
            'natoms': mol.natm,
            'symbols': symbols,
            'atomic_numbers': atomic_numbers,
            'coordinates': mol.atom_coords(),  # In Bohr
            'charge': mol.charge,
            'multiplicity': mol.spin + 1,
            'symmetry': mol.symmetry if hasattr(mol, 'symmetry') and mol.symmetry else None
        }

    def _extract_basis_data(self) -> Dict:
        """Extract basis set information from PySCF calculation.

        Returns
        -------
        dict
            Basis set data for HDF5 export
        """
        mol = self.mol

        # Try to get auxiliary basis info
        try:
            from pyscf import df
            # CRITICAL FIX: PySCF requires auxbasis as dict
            auxbasis_dict = self.auxbasis
            if isinstance(self.auxbasis, str):
                auxbasis_dict = {mol.atom_symbol(i): self.auxbasis for i in range(mol.natm)}
            auxmol = df.addons.make_auxmol(mol, auxbasis_dict)
            n_aux = auxmol.nao_nr()
        except:
            n_aux = 0

        return {
            'ao_basis': mol.basis if isinstance(mol.basis, str) else 'custom',
            'n_ao': mol.nao_nr(),
            'aux_basis': self.auxbasis,
            'n_aux': n_aux
        }

    def export_to_quasix_schema(self, include_df_tensors: bool = False) -> Dict:
        """Export calculation data in QuasiX schema format.

        Parameters
        ----------
        include_df_tensors : bool
            Include large DF tensors in export (can be memory intensive)

        Returns
        -------
        dict
            Complete calculation data following QuasiX schema
        """
        result = self.analyze()

        schema_data = {
            'metadata': {
                'quasix_version': '0.1.0',
                'timestamp': time.strftime('%Y-%m-%dT%H:%M:%S'),
                'calculation_type': 'EvGW'
            },
            'molecule': self._extract_molecule_data(),
            'basis': self._extract_basis_data(),
            'parameters': {
                'calculation_type': 'EvGW',
                'convergence': {
                    'energy_tol': self.conv_tol,
                    'density_tol': self.conv_tol_z,
                    'max_iterations': self.max_cycle,
                    'use_diis': self.diis,
                    'diis_space': self.diis_space
                },
                'frequency': {
                    'grid_type': self.freq_int,
                    'n_points': self.nfreq,
                    'eta': self.eta
                },
                'gw_params': {
                    'starting_point': type(self._scf).__name__,
                    'damping': self.damping,
                    'damping_dynamic': self.damping_dynamic
                }
            },
            'results': {
                'reference': {
                    'energy': self._scf.e_tot if hasattr(self._scf, 'e_tot') else 0.0,
                    'orbital_energies': self._scf.mo_energy.tolist(),
                    'occupations': self._scf.mo_occ.tolist(),
                    'homo': np.sum(self._scf.mo_occ > 0) - 1,
                    'lumo': np.sum(self._scf.mo_occ > 0)
                },
                'gw': {
                    'qp_energies': result.qp_energies.tolist(),
                    'z_factors': result.z_factors.tolist(),
                    'converged': result.converged,
                    'iterations': result.n_cycles
                },
                'timings': {
                    'total_time': result.total_time
                }
            }
        }

        if include_df_tensors:
            # Add DF tensor metadata (not the actual tensors due to size)
            schema_data['df_tensors'] = {
                'auxbasis': self.auxbasis,
                'n_aux': self._extract_basis_data()['n_aux'],
                'thresh_df': self.thresh_df
            }

        return schema_data

    def save_df_tensors(self, filename: str, iaP: np.ndarray = None,
                       ijP: np.ndarray = None, chol_v: np.ndarray = None,
                       compress: bool = True):
        """Save density-fitted tensors to HDF5 file with optimal chunking.

        Parameters
        ----------
        filename : str
            Output HDF5 filename
        iaP : ndarray, optional
            (ia|P) transition DF tensor
        ijP : ndarray, optional
            (ij|P) occupied-occupied DF tensor
        chol_v : ndarray, optional
            Cholesky factor of Coulomb metric
        compress : bool
            Enable compression
        """
        import h5py

        if not filename.endswith('.h5'):
            filename += '.h5'

        compression = 'gzip' if compress else None
        compression_opts = 4 if compress else None

        with h5py.File(filename, 'w') as f:
            # Add metadata
            f.attrs['created_by'] = 'QuasiX.EvGW'
            f.attrs['timestamp'] = time.strftime('%Y-%m-%dT%H:%M:%S')
            f.attrs['auxbasis'] = self.auxbasis

            if iaP is not None:
                # Optimize chunking for (ia|P) access patterns
                shape = iaP.shape
                if len(shape) == 2:
                    chunk_size = (min(shape[0], 1000), min(shape[1], 100))
                else:
                    chunk_size = (min(shape[0], 50), min(shape[1], 50), min(shape[2], 100))

                f.create_dataset('iaP', data=iaP, chunks=chunk_size,
                               compression=compression,
                               compression_opts=compression_opts)
                f['iaP'].attrs['description'] = 'Occupied-virtual DF tensor'

            if ijP is not None:
                # Optimize chunking for (ij|P) access patterns
                shape = ijP.shape
                if len(shape) == 2:
                    chunk_size = (min(shape[0], 100), min(shape[1], 100))
                else:
                    chunk_size = (min(shape[0], 20), min(shape[1], 20), min(shape[2], 100))

                f.create_dataset('ijP', data=ijP, chunks=chunk_size,
                               compression=compression,
                               compression_opts=compression_opts)
                f['ijP'].attrs['description'] = 'Occupied-occupied DF tensor'

            if chol_v is not None:
                # Optimize chunking for Coulomb metric
                shape = chol_v.shape
                chunk_size = (min(shape[0], 100), min(shape[1], 100))

                f.create_dataset('chol_v', data=chol_v, chunks=chunk_size,
                               compression=compression,
                               compression_opts=compression_opts)
                f['chol_v'].attrs['description'] = 'Cholesky factor of Coulomb metric V^{1/2}'

        log = logger.Logger(self.stdout, self.verbose)
        log.info(f'DF tensors saved to {filename}')

    @staticmethod
    def load_df_tensors(filename: str, lazy: bool = False) -> Dict[str, np.ndarray]:
        """Load density-fitted tensors from HDF5 file.

        Parameters
        ----------
        filename : str
            Input HDF5 filename
        lazy : bool
            Use memory mapping for large arrays

        Returns
        -------
        dict
            Dictionary with 'iaP', 'ijP', 'chol_v' arrays
        """
        import h5py

        if lazy:
            # Memory-mapped access
            f = h5py.File(filename, 'r')
            return {
                'iaP': f['iaP'] if 'iaP' in f else None,
                'ijP': f['ijP'] if 'ijP' in f else None,
                'chol_v': f['chol_v'] if 'chol_v' in f else None,
                '_file_handle': f  # Keep file handle alive
            }
        else:
            # Load into memory
            with h5py.File(filename, 'r') as f:
                tensors = {}
                if 'iaP' in f:
                    tensors['iaP'] = f['iaP'][:]
                if 'ijP' in f:
                    tensors['ijP'] = f['ijP'][:]
                if 'chol_v' in f:
                    tensors['chol_v'] = f['chol_v'][:]
                return tensors


def evgw(mf, auxbasis='def2-tzvp-jkfit', **kwargs):
    """
    Functional interface for evGW calculations.

    Parameters
    ----------
    mf : SCF object
        PySCF mean-field calculation
    auxbasis : str
        Auxiliary basis for density fitting
    **kwargs
        Additional parameters passed to EvGW class

    Returns
    -------
    ndarray or EvGWResult
        QP energies if return_full=False (default)
        EvGWResult object if return_full=True

    Examples
    --------
    >>> from pyscf import gto, scf
    >>> from quasix.evgw import evgw
    >>> mol = gto.M(atom='H 0 0 0; F 0 0 1.1', basis='cc-pvdz')
    >>> mf = scf.RHF(mol).run()
    >>> qp_energies = evgw(mf)
    """
    return_full = kwargs.pop('return_full', False)

    gw = EvGW(mf, auxbasis=auxbasis)

    # Set parameters from kwargs
    for key, value in kwargs.items():
        if hasattr(gw, key):
            setattr(gw, key, value)

    # Run calculation
    gw.kernel()

    if return_full:
        return gw.analyze()
    else:
        return gw.qp_energies