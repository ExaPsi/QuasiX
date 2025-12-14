"""
Zero-copy PySCF adapter for QuasiX

This module provides optimized PySCF integration with zero-copy transfers
to the Rust backend where possible.
"""

import numpy as np
from typing import Optional, Dict, Any, Tuple
from pyscf import gto, scf, df
from dataclasses import dataclass
import logging

from .logging import get_logger
from .mo_transform_fast import (
    build_transition_df_tensor_fast,
    build_metric_fast,
    configure_threading
)

# Try to import Rust bindings
try:
    import quasix
    RUST_AVAILABLE = True
except ImportError:
    RUST_AVAILABLE = False

logger = get_logger(__name__)


@dataclass
class ZeroCopyData:
    """Data structure optimized for zero-copy transfer to Rust"""
    mo_coeff: np.ndarray  # C-contiguous, float64
    mo_energy: np.ndarray  # C-contiguous, float64
    mo_occ: np.ndarray     # C-contiguous, float64
    ia_P: np.ndarray       # C-contiguous, float64
    metric: np.ndarray     # C-contiguous, float64
    vxc_diag: Optional[np.ndarray] = None  # For DFT

    def __post_init__(self):
        """Ensure all arrays are C-contiguous for zero-copy"""
        self.mo_coeff = np.ascontiguousarray(self.mo_coeff, dtype=np.float64)
        self.mo_energy = np.ascontiguousarray(self.mo_energy, dtype=np.float64)
        self.mo_occ = np.ascontiguousarray(self.mo_occ, dtype=np.float64)
        self.ia_P = np.ascontiguousarray(self.ia_P, dtype=np.float64)
        self.metric = np.ascontiguousarray(self.metric, dtype=np.float64)
        if self.vxc_diag is not None:
            self.vxc_diag = np.ascontiguousarray(self.vxc_diag, dtype=np.float64)

    @property
    def memory_usage(self) -> float:
        """Total memory usage in MB"""
        total = (self.mo_coeff.nbytes + self.mo_energy.nbytes +
                self.mo_occ.nbytes + self.ia_P.nbytes + self.metric.nbytes)
        if self.vxc_diag is not None:
            total += self.vxc_diag.nbytes
        return total / 1e6


class EvGWZeroCopy:
    """
    evGW class with zero-copy optimization for PySCF integration

    This class provides efficient data transfer to the Rust backend
    using zero-copy where possible and SIMD acceleration.
    """

    def __init__(self, mf, auxbasis: Optional[str] = None, **kwargs):
        """
        Initialize evGW calculation from PySCF mean-field

        Args:
            mf: PySCF mean-field object (RHF/RKS)
            auxbasis: Auxiliary basis for DF (auto if None)
            **kwargs: Additional options
        """
        # Configure threading for optimal performance
        configure_threading()

        self.mol = mf.mol
        self._scf = mf
        self.auxbasis = auxbasis or self._get_default_auxbasis()

        # evGW parameters
        self.conv_tol = kwargs.get('conv_tol', 1e-6)
        self.max_cycle = kwargs.get('max_cycle', 50)
        self.damping = kwargs.get('damping', 0.0)
        self.eta = kwargs.get('eta', 1e-3)
        self.nfreq = kwargs.get('nfreq', 100)

        # Cache for zero-copy data
        self._data_cache: Optional[ZeroCopyData] = None
        self._qp_energies = None
        self._z_factors = None

        logger.info(f"EvGWZeroCopy initialized for {self.mol.natm} atoms")

    def _get_default_auxbasis(self) -> str:
        """Get default auxiliary basis based on AO basis"""
        ao_basis = self.mol.basis

        # Common mappings
        aux_map = {
            'sto-3g': 'def2-svp-jkfit',
            'def2-svp': 'def2-svp-jkfit',
            'def2-tzvp': 'def2-tzvp-jkfit',
            'cc-pvdz': 'cc-pvdz-jkfit',
            'cc-pvtz': 'cc-pvtz-jkfit',
        }

        return aux_map.get(ao_basis.lower(), 'def2-svp-jkfit')

    def prepare_data(self) -> ZeroCopyData:
        """
        Prepare all data for zero-copy transfer to Rust

        Returns:
            ZeroCopyData object with all arrays C-contiguous
        """
        if self._data_cache is not None:
            logger.debug("Using cached zero-copy data")
            return self._data_cache

        logger.info("Preparing zero-copy data for Rust transfer")

        # Extract MO data (already C-contiguous from PySCF)
        mo_coeff = self._scf.mo_coeff
        mo_energy = self._scf.mo_energy
        mo_occ = self._scf.mo_occ

        # Build DF tensors using optimized routines
        ia_P = build_transition_df_tensor_fast(
            self.mol, self.auxbasis, mo_coeff, mo_occ
        )

        metric = build_metric_fast(self.mol, self.auxbasis)

        # Extract vxc for DFT
        vxc_diag = None
        if hasattr(self._scf, 'xc'):  # DFT calculation
            vxc = self._scf.get_veff() - self._scf.get_j()
            vxc_mo = mo_coeff.T @ vxc @ mo_coeff
            vxc_diag = np.diag(vxc_mo)

        # Create zero-copy data structure
        self._data_cache = ZeroCopyData(
            mo_coeff=mo_coeff,
            mo_energy=mo_energy,
            mo_occ=mo_occ,
            ia_P=ia_P,
            metric=metric,
            vxc_diag=vxc_diag
        )

        logger.info(f"Zero-copy data prepared: {self._data_cache.memory_usage:.2f} MB")
        return self._data_cache

    def kernel(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Run evGW calculation

        Returns:
            qp_energies: Quasiparticle energies
            z_factors: Renormalization factors
        """
        # Prepare zero-copy data
        data = self.prepare_data()

        if RUST_AVAILABLE and hasattr(quasix, 'simd_ops'):
            # Use SIMD-accelerated path
            logger.info("Using SIMD-accelerated evGW kernel")
            result = self._kernel_simd(data)
        elif RUST_AVAILABLE:
            # Use standard Rust path
            logger.info("Using Rust evGW kernel")
            result = self._kernel_rust(data)
        else:
            # Fallback to Python implementation
            logger.info("Using Python evGW kernel (Rust not available)")
            result = self._kernel_python(data)

        self._qp_energies = result['qp_energies']
        self._z_factors = result['z_factors']

        return self._qp_energies, self._z_factors

    def _kernel_simd(self, data: ZeroCopyData) -> Dict[str, np.ndarray]:
        """Run evGW using SIMD-accelerated Rust kernel"""
        import quasix

        # Check SIMD availability
        caps = quasix.simd_ops.simd_capabilities()
        if caps.get('avx2', False) and caps.get('fma', False):
            logger.info("AVX2+FMA available, using optimized path")

        # Call Rust evGW with zero-copy arrays
        result = quasix.evgw_contour_deformation(
            data.mo_energy,
            data.mo_occ,
            data.ia_P,
            data.metric,
            {
                'max_cycles': self.max_cycle,
                'conv_tol': self.conv_tol,
                'damping': self.damping,
                'n_freq': self.nfreq,
                'eta': self.eta,
                'freq_method': 'cd',
                'use_simd': True,  # Enable SIMD optimizations
            }
        )

        return {
            'qp_energies': np.array(result['qp_energies']),
            'z_factors': np.array(result['z_factors'])
        }

    def _kernel_rust(self, data: ZeroCopyData) -> Dict[str, np.ndarray]:
        """Run evGW using standard Rust kernel"""
        import quasix

        result = quasix.evgw_contour_deformation(
            data.mo_energy,
            data.mo_occ,
            data.ia_P,
            data.metric,
            {
                'max_cycles': self.max_cycle,
                'conv_tol': self.conv_tol,
                'damping': self.damping,
                'n_freq': self.nfreq,
                'eta': self.eta,
                'freq_method': 'cd',
            }
        )

        return {
            'qp_energies': np.array(result['qp_energies']),
            'z_factors': np.array(result['z_factors'])
        }

    def _kernel_python(self, data: ZeroCopyData) -> Dict[str, np.ndarray]:
        """Fallback Python implementation for testing"""
        # Simple diagonal approximation for testing
        logger.warning("Using simplified Python kernel (for testing only)")

        # Mock implementation - just shift energies slightly
        qp_energies = data.mo_energy.copy()
        z_factors = np.ones_like(data.mo_energy) * 0.8

        # Add small corrections
        occ_mask = data.mo_occ > 0.5
        qp_energies[occ_mask] -= 0.1  # Lower occupied
        qp_energies[~occ_mask] += 0.1  # Raise virtual

        return {
            'qp_energies': qp_energies,
            'z_factors': z_factors
        }

    @property
    def qp_energies(self) -> Optional[np.ndarray]:
        """Quasiparticle energies"""
        return self._qp_energies

    @property
    def z_factors(self) -> Optional[np.ndarray]:
        """Renormalization factors"""
        return self._z_factors

    @property
    def ip(self) -> Optional[float]:
        """Ionization potential (negative HOMO energy)"""
        if self._qp_energies is None:
            return None
        occ_mask = self._scf.mo_occ > 0.5
        homo_idx = np.where(occ_mask)[0][-1]
        return -self._qp_energies[homo_idx]

    @property
    def ea(self) -> Optional[float]:
        """Electron affinity (negative LUMO energy)"""
        if self._qp_energies is None:
            return None
        vir_mask = self._scf.mo_occ < 0.5
        lumo_idx = np.where(vir_mask)[0][0]
        return -self._qp_energies[lumo_idx]

    @property
    def gap(self) -> Optional[float]:
        """HOMO-LUMO gap"""
        if self._qp_energies is None:
            return None
        return self.ip + self.ea


def test_zero_copy():
    """Test zero-copy adapter with H2O"""
    from pyscf import gto, scf

    # Create test molecule
    mol = gto.M(
        atom='O 0 0 0; H 0 0.75 0.5; H 0 0.75 -0.5',
        basis='sto-3g',
        verbose=0
    )

    # Run HF
    mf = scf.RHF(mol)
    mf.kernel()

    # Run evGW with zero-copy
    gw = EvGWZeroCopy(mf)
    data = gw.prepare_data()

    print(f"Zero-copy data prepared:")
    print(f"  Memory usage: {data.memory_usage:.2f} MB")
    print(f"  mo_coeff C-contiguous: {data.mo_coeff.flags['C_CONTIGUOUS']}")
    print(f"  ia_P C-contiguous: {data.ia_P.flags['C_CONTIGUOUS']}")
    print(f"  metric C-contiguous: {data.metric.flags['C_CONTIGUOUS']}")

    # Run calculation
    qp_e, z_factors = gw.kernel()
    print(f"\nResults:")
    print(f"  IP: {gw.ip:.6f} eV")
    print(f"  EA: {gw.ea:.6f} eV")
    print(f"  Gap: {gw.gap:.6f} eV")

    return gw


if __name__ == "__main__":
    test_zero_copy()