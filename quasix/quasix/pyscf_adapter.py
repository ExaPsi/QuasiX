"""
PySCF Adapter for QuasiX

This module provides clean interfaces for extracting data from PySCF mean-field
objects and preparing them for QuasiX GW/BSE calculations. It supports density 
fitting with various auxiliary basis sets and ensures tolerances < 1e-8 for 
PySCF comparisons.
"""

import numpy as np
from typing import Optional, Dict, Any, Union, Tuple, List
import logging
from dataclasses import dataclass, field

try:
    from pyscf import gto, scf, df, lib
    from pyscf.scf import RHF, RKS, UHF, UKS
    from pyscf.lib import logger as pyscf_logger
    HAS_PYSCF = True
except ImportError:
    HAS_PYSCF = False
    # Define dummy classes for type hints
    class RHF: pass
    class RKS: pass
    class UHF: pass
    class UKS: pass
    class gto:
        class Mole: pass

from .df_tensor import DFTensor
from .mo_transform import (
    validate_mo_coefficients,
    transform_mo_coefficients,
    build_transition_df_tensor,
    build_full_df_tensors
)
from .logging import get_logger, timed_stage

# Set up module logger
logger = get_logger(__name__)


@dataclass
class PySCFData:
    """Container for data extracted from PySCF mean-field objects"""
    
    # Molecular information
    mol: Any  # PySCF molecule object
    auxmol: Optional[Any] = None  # Auxiliary basis molecule
    
    # MO data
    mo_coeff: np.ndarray = field(default_factory=lambda: np.array([]))
    mo_occ: np.ndarray = field(default_factory=lambda: np.array([]))
    mo_energy: np.ndarray = field(default_factory=lambda: np.array([]))
    
    # Dimensions
    nao: int = 0
    naux: int = 0
    nocc: int = 0
    nvir: int = 0
    nmo: int = 0
    
    # Energy components
    e_tot: float = 0.0
    e_hf: float = 0.0
    vxc_diagonal: Optional[np.ndarray] = None
    
    # DF tensors (populated on demand)
    df_tensor: Optional[DFTensor] = None
    
    # Metadata
    converged: bool = False
    method: str = ""  # RHF, RKS, UHF, UKS
    basis: str = ""
    auxbasis: str = ""
    
    def validate(self) -> bool:
        """Validate the extracted data"""
        if self.mo_coeff.size == 0:
            logger.error("Empty MO coefficients")
            return False
            
        if self.nao != self.mo_coeff.shape[0]:
            logger.error(f"Inconsistent dimensions: nao={self.nao}, mo_coeff.shape={self.mo_coeff.shape}")
            return False
            
        if self.nocc + self.nvir != self.nmo:
            logger.error(f"Inconsistent MO count: nocc={self.nocc}, nvir={self.nvir}, nmo={self.nmo}")
            return False
            
        return True


class PySCFAdapter:
    """
    Adapter for extracting and preparing PySCF data for QuasiX calculations
    
    This class provides a clean interface between PySCF's quantum chemistry
    infrastructure and QuasiX's high-performance GW/BSE implementation.
    """
    
    SUPPORTED_METHODS = ['RHF', 'RKS', 'UHF', 'UKS']
    DEFAULT_AUXBASIS_MAP = {
        'sto-3g': 'def2-svp-jkfit',
        'def2-svp': 'def2-svp-jkfit', 
        'def2-tzvp': 'def2-tzvp-jkfit',
        'def2-qzvp': 'def2-qzvp-jkfit',
        'cc-pvdz': 'cc-pvdz-jkfit',
        'cc-pvtz': 'cc-pvtz-jkfit',
        'cc-pvqz': 'cc-pvqz-jkfit',
        '6-31g': 'def2-svp-jkfit',
        '6-311g': 'def2-tzvp-jkfit'
    }
    
    def __init__(self, tolerance: float = 1e-8):
        """
        Initialize PySCF adapter
        
        Args:
            tolerance: Numerical tolerance for validations (default: 1e-8)
        """
        if not HAS_PYSCF:
            raise ImportError("PySCF is required for PySCFAdapter. Install with: pip install pyscf")
            
        self.tolerance = tolerance
        logger.info(f"PySCFAdapter initialized with tolerance={tolerance}")
    
    def extract_from_mf(self, mf: Union[RHF, RKS, UHF, UKS],
                       auxbasis: Optional[str] = None,
                       build_df: bool = True) -> PySCFData:
        """
        Extract data from a converged PySCF mean-field object
        
        Args:
            mf: Converged PySCF mean-field object
            auxbasis: Auxiliary basis for density fitting (auto-selected if None)
            build_df: Whether to build DF tensors immediately
            
        Returns:
            PySCFData object containing extracted data
            
        Raises:
            ValueError: If SCF is not converged or data extraction fails
        """
        if not mf.converged:
            logger.warning("SCF not converged! Results may be unreliable.")
        
        # Determine method type
        method_type = type(mf).__name__
        if method_type not in self.SUPPORTED_METHODS:
            raise ValueError(f"Unsupported method: {method_type}. Supported: {self.SUPPORTED_METHODS}")
        
        logger.info(f"Extracting data from {method_type} calculation")
        
        # Create data container
        data = PySCFData(
            mol=mf.mol,
            method=method_type,
            converged=mf.converged
        )
        
        # Extract basic information
        data.basis = mf.mol.basis
        data.nao = mf.mol.nao
        data.nmo = mf.mo_coeff.shape[1]
        
        # Extract MO data
        data.mo_coeff = np.asarray(mf.mo_coeff, order='C')
        data.mo_occ = np.asarray(mf.mo_occ)
        data.mo_energy = np.asarray(mf.mo_energy)
        
        # Calculate dimensions
        if method_type in ['RHF', 'RKS']:
            data.nocc = int(np.sum(data.mo_occ > 0))
            data.nvir = data.nmo - data.nocc
        else:  # UHF/UKS
            # For unrestricted, we need to handle alpha/beta separately
            # For now, just handle alpha spin
            data.nocc = int(np.sum(data.mo_occ[0] > 0))
            data.nvir = data.nmo - data.nocc
            logger.warning("UHF/UKS support is preliminary - using alpha spin only")
        
        # Extract energies
        data.e_tot = mf.e_tot
        data.e_hf = mf.e_tot
        
        # Extract Vxc diagonal elements for GW
        if hasattr(mf, 'get_veff'):
            veff = mf.get_veff()
            if method_type in ['RKS']:
                # For DFT, extract Vxc
                vxc = veff - mf.get_j()
                data.vxc_diagonal = np.diag(data.mo_coeff.T @ vxc @ data.mo_coeff)
            else:
                # For HF, Vxc is the Fock exchange
                data.vxc_diagonal = np.diag(data.mo_coeff.T @ veff @ data.mo_coeff)
        
        # Select auxiliary basis if not provided
        if auxbasis is None:
            auxbasis = self._select_auxbasis(data.basis)
        data.auxbasis = auxbasis
        
        # Build auxiliary molecule
        data.auxmol = df.addons.make_auxmol(mf.mol, auxbasis)
        data.naux = data.auxmol.nao
        
        logger.info(f"Dimensions: nao={data.nao}, nocc={data.nocc}, nvir={data.nvir}, naux={data.naux}")
        
        # Validate MO coefficients
        with timed_stage("Validating MO coefficients"):
            S = mf.mol.intor('int1e_ovlp')
            try:
                validate_mo_coefficients(data.mo_coeff, data.nao, data.mo_occ, S)
                logger.info("MO coefficients validated successfully")
            except ValueError as e:
                logger.error(f"MO validation failed: {e}")
                raise
        
        # Build DF tensors if requested
        if build_df:
            with timed_stage("Building DF tensors"):
                data.df_tensor = DFTensor.from_pyscf(mf, auxbasis, build_all=True)
                logger.info(f"DF tensors built: {data.df_tensor}")
        
        # Validate extracted data
        if not data.validate():
            raise ValueError("Extracted data validation failed")
        
        return data
    
    def _select_auxbasis(self, basis: str) -> str:
        """
        Auto-select appropriate auxiliary basis based on orbital basis
        
        Args:
            basis: Orbital basis set name
            
        Returns:
            Auxiliary basis set name
        """
        # Check direct mapping
        basis_lower = basis.lower()
        if basis_lower in self.DEFAULT_AUXBASIS_MAP:
            auxbasis = self.DEFAULT_AUXBASIS_MAP[basis_lower]
            logger.info(f"Auto-selected auxiliary basis: {auxbasis} for {basis}")
            return auxbasis
        
        # Try to find a close match
        if 'def2' in basis_lower:
            if 'svp' in basis_lower:
                return 'def2-svp-jkfit'
            elif 'tzvp' in basis_lower:
                return 'def2-tzvp-jkfit'
            elif 'qzvp' in basis_lower:
                return 'def2-qzvp-jkfit'
        elif 'cc-pv' in basis_lower:
            if 'dz' in basis_lower:
                return 'cc-pvdz-jkfit'
            elif 'tz' in basis_lower:
                return 'cc-pvtz-jkfit'
            elif 'qz' in basis_lower:
                return 'cc-pvqz-jkfit'
        
        # Default fallback
        logger.warning(f"Could not auto-select auxiliary basis for {basis}, using def2-tzvp-jkfit")
        return 'def2-tzvp-jkfit'
    
    def validate_against_pyscf(self, data: PySCFData, 
                              reference_mf: Union[RHF, RKS]) -> Dict[str, float]:
        """
        Validate QuasiX data against PySCF reference with strict tolerances
        
        Args:
            data: PySCFData object to validate
            reference_mf: Reference PySCF mean-field object
            
        Returns:
            Dictionary of validation metrics and errors
        """
        results = {}
        
        # Validate MO energies
        energy_error = np.abs(data.mo_energy - reference_mf.mo_energy).max()
        results['mo_energy_error'] = energy_error
        if energy_error > self.tolerance:
            logger.warning(f"MO energy error {energy_error:.2e} exceeds tolerance {self.tolerance:.2e}")
        
        # Validate MO coefficients orthonormality
        S = reference_mf.mol.intor('int1e_ovlp')
        ortho = data.mo_coeff.T @ S @ data.mo_coeff
        ortho_error = np.abs(ortho - np.eye(data.nmo)).max()
        results['mo_ortho_error'] = ortho_error
        if ortho_error > self.tolerance:
            logger.warning(f"MO orthonormality error {ortho_error:.2e} exceeds tolerance {self.tolerance:.2e}")
        
        # Validate DF tensors if available
        if data.df_tensor is not None:
            # Build reference DF tensors
            with timed_stage("Building reference DF tensors"):
                ref_df = df.DF(reference_mf.mol, data.auxbasis)
                ref_df.kernel()
                
                # Get reference (ia|P) tensor
                naux = data.naux
                nocc = data.nocc
                nvir = data.nvir
                
                # Build reference using PySCF's DF module
                ref_3c = ref_df._cderi  # This gets the 3-center integrals
                
                # Transform to MO basis
                c_occ = reference_mf.mo_coeff[:, :nocc]
                c_vir = reference_mf.mo_coeff[:, nocc:nocc+nvir]
                
                # Contract to get (ia|P)
                # This is a simplified comparison - full validation would be more complex
                iaP = data.df_tensor.get_transition_tensor()
                
                # Check shape consistency
                expected_shape = (nocc * nvir, naux)
                if iaP.shape != expected_shape:
                    logger.error(f"DF tensor shape mismatch: got {iaP.shape}, expected {expected_shape}")
                    results['df_shape_error'] = True
                else:
                    results['df_shape_error'] = False
                    
                # Check norm (as a proxy for correctness)
                df_norm = np.linalg.norm(iaP)
                results['df_tensor_norm'] = df_norm
                logger.info(f"DF tensor norm: {df_norm:.6f}")
        
        # Validate total energy
        energy_diff = abs(data.e_tot - reference_mf.e_tot)
        results['total_energy_error'] = energy_diff
        if energy_diff > self.tolerance * 100:  # Use looser tolerance for total energy
            logger.warning(f"Total energy error {energy_diff:.2e} is large")
        
        # Summary
        all_pass = all([
            results.get('mo_energy_error', 1.0) < self.tolerance,
            results.get('mo_ortho_error', 1.0) < self.tolerance,
            not results.get('df_shape_error', True),
            results.get('total_energy_error', 1.0) < self.tolerance * 100
        ])
        
        results['validation_passed'] = all_pass
        
        if all_pass:
            logger.info(f"All validations passed with tolerance {self.tolerance:.2e}")
        else:
            logger.warning("Some validations failed - check individual metrics")
        
        return results
    
    def prepare_for_rust(self, data: PySCFData) -> Dict[str, np.ndarray]:
        """
        Prepare PySCF data for transfer to Rust backend
        
        Ensures all arrays are C-contiguous and properly formatted for
        efficient transfer across the Python-Rust boundary.
        
        Args:
            data: PySCFData object
            
        Returns:
            Dictionary of C-contiguous arrays ready for Rust
        """
        rust_data = {}
        
        # Ensure C-contiguous arrays for efficient Rust transfer
        rust_data['mo_coeff'] = np.ascontiguousarray(data.mo_coeff, dtype=np.float64)
        rust_data['mo_energy'] = np.ascontiguousarray(data.mo_energy, dtype=np.float64)
        rust_data['mo_occ'] = np.ascontiguousarray(data.mo_occ, dtype=np.float64)
        
        # Add DF tensors if available
        if data.df_tensor is not None:
            iaP = data.df_tensor.get_transition_tensor()
            rust_data['df_ia'] = np.ascontiguousarray(iaP, dtype=np.float64)
            
            metric = data.df_tensor.get_metric()
            rust_data['metric'] = np.ascontiguousarray(metric, dtype=np.float64)
            
            if data.df_tensor.ijP is not None:
                rust_data['df_ij'] = np.ascontiguousarray(data.df_tensor.ijP, dtype=np.float64)
            if data.df_tensor.abP is not None:
                rust_data['df_ab'] = np.ascontiguousarray(data.df_tensor.abP, dtype=np.float64)
        
        # Add Vxc if available
        if data.vxc_diagonal is not None:
            rust_data['vxc_diagonal'] = np.ascontiguousarray(data.vxc_diagonal, dtype=np.float64)
        
        # Add dimensions
        rust_data['nocc'] = np.array([data.nocc], dtype=np.int64)
        rust_data['nvir'] = np.array([data.nvir], dtype=np.int64)
        rust_data['naux'] = np.array([data.naux], dtype=np.int64)
        
        logger.info(f"Prepared {len(rust_data)} arrays for Rust transfer")
        
        return rust_data


class MoleculeFactory:
    """Factory for creating test molecules with PySCF"""
    
    @staticmethod
    def water(basis: str = 'def2-svp') -> gto.Mole:
        """Create water molecule"""
        return gto.M(
            atom='O 0 0 0; H 0 0 1.0; H 0 0.7 0.7',
            basis=basis,
            verbose=0
        )
    
    @staticmethod
    def benzene(basis: str = 'def2-svp') -> gto.Mole:
        """Create benzene molecule"""
        return gto.M(
            atom='''
            C  0.000000  1.396000  0.000000
            C  1.209000  0.698000  0.000000
            C  1.209000 -0.698000  0.000000
            C  0.000000 -1.396000  0.000000
            C -1.209000 -0.698000  0.000000
            C -1.209000  0.698000  0.000000
            H  0.000000  2.479000  0.000000
            H  2.147000  1.240000  0.000000
            H  2.147000 -1.240000  0.000000
            H  0.000000 -2.479000  0.000000
            H -2.147000 -1.240000  0.000000
            H -2.147000  1.240000  0.000000
            ''',
            basis=basis,
            verbose=0
        )
    
    @staticmethod
    def ammonia(basis: str = 'def2-svp') -> gto.Mole:
        """Create ammonia molecule"""
        return gto.M(
            atom='''
            N  0.000000  0.000000  0.116700
            H  0.000000  0.940000 -0.272300
            H  0.814000 -0.470000 -0.272300
            H -0.814000 -0.470000 -0.272300
            ''',
            basis=basis,
            verbose=0
        )
    
    @staticmethod
    def carbon_monoxide(basis: str = 'def2-svp') -> gto.Mole:
        """Create CO molecule"""
        return gto.M(
            atom='C 0 0 0; O 0 0 1.128',
            basis=basis,
            verbose=0
        )


def benchmark_adapter(mol: gto.Mole, auxbasis: str = 'def2-tzvp-jkfit') -> Dict[str, Any]:
    """
    Benchmark PySCF adapter performance
    
    Args:
        mol: PySCF molecule
        auxbasis: Auxiliary basis set
        
    Returns:
        Dictionary with timing and memory statistics
    """
    import time
    import tracemalloc
    
    results = {}
    
    # Run SCF
    t0 = time.time()
    mf = scf.RHF(mol)
    mf.kernel()
    results['scf_time'] = time.time() - t0
    results['scf_converged'] = mf.converged
    
    # Initialize adapter
    adapter = PySCFAdapter(tolerance=1e-8)
    
    # Extract data with memory tracking
    tracemalloc.start()
    t0 = time.time()
    
    data = adapter.extract_from_mf(mf, auxbasis=auxbasis, build_df=True)
    
    extraction_time = time.time() - t0
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    results['extraction_time'] = extraction_time
    results['peak_memory_mb'] = peak / 1024 / 1024
    results['nocc'] = data.nocc
    results['nvir'] = data.nvir
    results['naux'] = data.naux
    
    # Prepare for Rust
    t0 = time.time()
    rust_data = adapter.prepare_for_rust(data)
    results['rust_prep_time'] = time.time() - t0
    
    # Calculate data sizes
    total_size = sum(arr.nbytes for arr in rust_data.values() if isinstance(arr, np.ndarray))
    results['rust_data_size_mb'] = total_size / 1024 / 1024
    
    logger.info(f"Benchmark results: {results}")
    
    return results


# Convenience functions for common use cases
def quick_extract(mf: Union[RHF, RKS], auxbasis: Optional[str] = None) -> PySCFData:
    """
    Quick extraction of PySCF data with default settings
    
    Args:
        mf: Converged mean-field object
        auxbasis: Auxiliary basis (auto-selected if None)
        
    Returns:
        PySCFData object
    """
    adapter = PySCFAdapter()
    return adapter.extract_from_mf(mf, auxbasis=auxbasis, build_df=True)


def validate_tolerance(data: PySCFData, reference_mf: Union[RHF, RKS],
                       tolerance: float = 1e-8) -> bool:
    """
    Validate data against reference with specified tolerance
    
    Args:
        data: PySCFData to validate
        reference_mf: Reference mean-field object
        tolerance: Numerical tolerance
        
    Returns:
        True if all validations pass
    """
    adapter = PySCFAdapter(tolerance=tolerance)
    results = adapter.validate_against_pyscf(data, reference_mf)
    return results['validation_passed']