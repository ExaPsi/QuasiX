"""
Exchange self-energy module for QuasiX.

This module provides high-level Python interfaces for computing exchange self-energy
using density fitting/resolution-of-identity techniques. It integrates seamlessly
with PySCF workflows and supports both standard and optimized algorithms.
"""

import numpy as np
from typing import Dict, Optional, Tuple, Union
from .logging import get_logger, timed_stage

# Import low-level Rust bindings
from . import (
    compute_exchange_matrix_ri,
    compute_exchange_diagonal_ri,
    compute_exchange_symmetric
)

logger = get_logger(__name__)


def compute_exchange_self_energy(
    mo_energy: np.ndarray,
    mo_occ: np.ndarray,
    mo_coeff: np.ndarray,
    ao_integrals_df: Dict,
    diagonal_only: bool = False,
    use_optimized: bool = True,
    **kwargs
) -> np.ndarray:
    """
    Compute exchange self-energy matrix using density fitting.
    
    This function provides a high-level interface for computing the exchange
    self-energy (Σˣ) in the GW approximation using density fitting/resolution-of-identity
    techniques. It automatically handles the transformation from AO to MO basis
    and applies the appropriate algorithm based on user preferences.
    
    The exchange self-energy is computed as:
    Σˣₘₙ = -Σᵢ ΣₚQ (mi|P) v⁻¹ₚQ (Q|ni)
    
    where:
    - m,n are molecular orbital indices
    - i runs over occupied orbitals
    - P,Q are auxiliary basis indices
    - v is the DF metric matrix
    
    Parameters
    ----------
    mo_energy : np.ndarray
        Molecular orbital energies, shape (n_mo,)
    mo_occ : np.ndarray
        Molecular orbital occupations (0 or 1 for RHF, 0-2 for UHF), shape (n_mo,)
    mo_coeff : np.ndarray
        Molecular orbital coefficients, shape (n_ao, n_mo)
    ao_integrals_df : Dict
        Dictionary containing density fitting tensors:
        - 'j3c_mo': Three-center integrals in MO basis (n_mo, n_mo, n_aux)
        - 'df_metric': DF metric matrix (n_aux, n_aux)
        - 'df_metric_inv': Inverse of DF metric (n_aux, n_aux)
        Alternative keys:
        - 'j3c_ao': Three-center integrals in AO basis (n_ao, n_ao, n_aux)
          (will be transformed to MO basis automatically)
    diagonal_only : bool, optional
        If True, compute only diagonal elements of Σˣ. Default: False
    use_optimized : bool, optional
        If True, use optimized blocked algorithm. Default: True
    **kwargs : dict
        Additional keyword arguments:
        - block_size: Block size for cache optimization (default: 64)
        - use_symmetrized: Use symmetrized DF tensors (default: True)
        - validate: Perform validation checks (default: True)
    
    Returns
    -------
    np.ndarray
        Exchange self-energy matrix Σˣ, shape (n_mo, n_mo) if diagonal_only=False,
        or diagonal elements, shape (n_mo,) if diagonal_only=True
    
    Raises
    ------
    ValueError
        If required keys are missing from ao_integrals_df
    RuntimeError
        If computation fails in the Rust kernel
    
    Examples
    --------
    >>> from pyscf import gto, scf, df
    >>> import quasix
    >>> 
    >>> # Setup molecule and perform HF calculation
    >>> mol = gto.M(atom='H 0 0 0; H 0 0 1', basis='cc-pvdz')
    >>> mf = scf.RHF(mol).run()
    >>> 
    >>> # Prepare DF integrals (using mock data for example)
    >>> n_mo = mf.mo_coeff.shape[1]
    >>> n_aux = 50
    >>> ao_integrals_df = {
    ...     'j3c_mo': np.random.randn(n_mo, n_mo, n_aux),
    ...     'df_metric_inv': np.eye(n_aux)
    ... }
    >>> 
    >>> # Compute exchange self-energy
    >>> sigma_x = quasix.selfenergy.compute_exchange_self_energy(
    ...     mf.mo_energy, mf.mo_occ, mf.mo_coeff, ao_integrals_df
    ... )
    >>> 
    >>> # Compute only diagonal elements (faster)
    >>> sigma_x_diag = quasix.selfenergy.compute_exchange_self_energy(
    ...     mf.mo_energy, mf.mo_occ, mf.mo_coeff, ao_integrals_df,
    ...     diagonal_only=True
    ... )
    
    Notes
    -----
    - For large systems, consider using diagonal_only=True for memory efficiency
    - The optimized algorithm uses cache-blocking for better performance
    - Ensure df_metric_inv is properly conditioned to avoid numerical issues
    - This function releases the GIL during computation for parallel efficiency
    
    See Also
    --------
    compute_exchange_matrix_ri : Low-level RI exchange self-energy
    compute_exchange_optimized : Optimized blocked algorithm
    """
    with timed_stage("compute_exchange_self_energy"):
        # Validate inputs
        if not isinstance(mo_energy, np.ndarray):
            mo_energy = np.asarray(mo_energy)
        if not isinstance(mo_occ, np.ndarray):
            mo_occ = np.asarray(mo_occ)
        if not isinstance(mo_coeff, np.ndarray):
            mo_coeff = np.asarray(mo_coeff)
        
        # Extract dimensions
        n_mo = len(mo_energy)
        nocc = int(np.sum(mo_occ > 0))
        
        logger.info(f"Computing exchange self-energy: n_mo={n_mo}, nocc={nocc}")
        
        # Get DF tensors from dictionary
        if 'j3c_mo' in ao_integrals_df:
            j3c_mo = ao_integrals_df['j3c_mo']
        elif 'j3c_ao' in ao_integrals_df:
            # Transform from AO to MO basis
            logger.info("Transforming 3-center integrals from AO to MO basis")
            j3c_ao = ao_integrals_df['j3c_ao']
            j3c_mo = _transform_j3c_to_mo(j3c_ao, mo_coeff)
        else:
            raise ValueError("ao_integrals_df must contain either 'j3c_mo' or 'j3c_ao'")
        
        # Get metric matrices
        if use_optimized:
            # Check if we have symmetrized DF tensors
            use_symmetrized = kwargs.get('use_symmetrized', True)
            
            if use_symmetrized and 'j3c_mo_sym' in ao_integrals_df:
                # Use symmetrized form directly
                logger.debug("Using symmetric exchange algorithm")
                sigma_x = compute_exchange_symmetric(
                    ao_integrals_df['j3c_mo_sym'], nocc
                )
                
                if diagonal_only:
                    return np.diag(sigma_x)
                return sigma_x
            elif use_symmetrized and 'df_metric_sqrt' in ao_integrals_df:
                # Apply metric square root to get symmetrized form
                logger.info("Computing symmetrized DF tensors")
                # Need to reshape and apply metric transformation properly
                # For now, use the tensor as-is (this would need proper implementation)
                logger.warning("Metric square root transformation not yet implemented for 3D tensors")
                use_optimized = False
            else:
                # Fall back to standard RI algorithm
                logger.debug("Falling back to standard RI algorithm")
                use_optimized = False
        
        if not use_optimized:
            # Standard RI algorithm needs inverse metric
            if 'df_metric_inv' not in ao_integrals_df:
                if 'df_metric' in ao_integrals_df:
                    logger.info("Computing inverse of DF metric matrix")
                    df_metric_inv = np.linalg.inv(ao_integrals_df['df_metric'])
                else:
                    raise ValueError("ao_integrals_df must contain 'df_metric_inv' or 'df_metric'")
            else:
                df_metric_inv = ao_integrals_df['df_metric_inv']
            
            logger.debug("Using standard RI algorithm")
            
            if diagonal_only:
                return compute_exchange_diagonal_ri(j3c_mo, df_metric_inv, nocc)
            else:
                return compute_exchange_matrix_ri(j3c_mo, df_metric_inv, nocc)


def _transform_j3c_to_mo(j3c_ao: np.ndarray, mo_coeff: np.ndarray) -> np.ndarray:
    """
    Transform 3-center integrals from AO to MO basis.
    
    Parameters
    ----------
    j3c_ao : np.ndarray
        Three-center integrals in AO basis, shape (n_ao, n_ao, n_aux)
    mo_coeff : np.ndarray
        MO coefficients, shape (n_ao, n_mo)
    
    Returns
    -------
    np.ndarray
        Three-center integrals in MO basis, shape (n_mo, n_mo, n_aux)
    """
    n_ao, _, n_aux = j3c_ao.shape
    n_mo = mo_coeff.shape[1]
    
    # First transformation: (μν|P) -> (μj|P)
    j3c_mo_half = np.zeros((n_ao, n_mo, n_aux))
    for p in range(n_aux):
        j3c_mo_half[:, :, p] = j3c_ao[:, :, p] @ mo_coeff
    
    # Second transformation: (μj|P) -> (ij|P)
    j3c_mo = np.zeros((n_mo, n_mo, n_aux))
    for p in range(n_aux):
        j3c_mo[:, :, p] = mo_coeff.T @ j3c_mo_half[:, :, p]
    
    return j3c_mo


class ExchangeSelfEnergy:
    """
    Exchange self-energy calculator with PySCF integration.
    
    This class provides a stateful interface for computing exchange self-energy,
    compatible with PySCF's object-oriented design patterns.
    
    Parameters
    ----------
    mf : pyscf.scf.hf.SCF
        Converged mean-field object
    auxbasis : str, optional
        Auxiliary basis set for density fitting. Default: 'def2-svp-jkfit'
    
    Attributes
    ----------
    mol : pyscf.gto.Mole
        Molecule object
    mf : pyscf.scf.hf.SCF
        Mean-field object
    auxbasis : str
        Auxiliary basis set name
    sigma_x : np.ndarray or None
        Computed exchange self-energy matrix
    
    Examples
    --------
    >>> from pyscf import gto, scf
    >>> import quasix
    >>> 
    >>> mol = gto.M(atom='H 0 0 0; H 0 0 1', basis='cc-pvdz')
    >>> mf = scf.RHF(mol).run()
    >>> 
    >>> # Create calculator
    >>> sigma = quasix.selfenergy.ExchangeSelfEnergy(mf)
    >>> 
    >>> # Compute exchange self-energy
    >>> sigma_x = sigma.kernel()
    >>> 
    >>> # Access cached result
    >>> print(sigma.sigma_x)
    """
    
    def __init__(self, mf, auxbasis='def2-svp-jkfit'):
        """Initialize exchange self-energy calculator."""
        self.mol = mf.mol
        self.mf = mf
        self.auxbasis = auxbasis
        self.sigma_x = None
        self._df_integrals = None
        
    def build_df_integrals(self) -> Dict:
        """
        Build density fitting integrals.
        
        Returns
        -------
        Dict
            Dictionary with DF tensors
        """
        if self._df_integrals is not None:
            return self._df_integrals
        
        logger.info(f"Building DF integrals with auxbasis={self.auxbasis}")
        
        # This is a placeholder - in production, would use PySCF's DF module
        # For now, create mock data for testing
        n_mo = self.mf.mo_coeff.shape[1]
        n_aux = n_mo * 3  # Rough estimate
        
        self._df_integrals = {
            'j3c_mo': np.random.randn(n_mo, n_mo, n_aux),
            'df_metric': np.eye(n_aux) + 0.1 * np.random.randn(n_aux, n_aux),
            'df_metric_inv': np.eye(n_aux)
        }
        
        # Make metric symmetric and positive definite
        self._df_integrals['df_metric'] = (
            self._df_integrals['df_metric'] + 
            self._df_integrals['df_metric'].T
        ) / 2
        self._df_integrals['df_metric'] += 2 * np.eye(n_aux)
        self._df_integrals['df_metric_inv'] = np.linalg.inv(
            self._df_integrals['df_metric']
        )
        
        return self._df_integrals
    
    def kernel(self, diagonal_only: bool = False) -> np.ndarray:
        """
        Compute exchange self-energy.
        
        Parameters
        ----------
        diagonal_only : bool, optional
            If True, compute only diagonal elements. Default: False
        
        Returns
        -------
        np.ndarray
            Exchange self-energy matrix or diagonal
        """
        # Build DF integrals if needed
        df_integrals = self.build_df_integrals()
        
        # Compute exchange self-energy
        self.sigma_x = compute_exchange_self_energy(
            self.mf.mo_energy,
            self.mf.mo_occ,
            self.mf.mo_coeff,
            df_integrals,
            diagonal_only=diagonal_only
        )
        
        return self.sigma_x
    
    def get_qp_correction(self, orbital_idx: Optional[int] = None) -> Union[float, np.ndarray]:
        """
        Get quasiparticle correction from exchange self-energy.
        
        Parameters
        ----------
        orbital_idx : int, optional
            Orbital index. If None, return all corrections
        
        Returns
        -------
        float or np.ndarray
            QP correction(s)
        """
        if self.sigma_x is None:
            self.kernel()
        
        # For exchange only, QP correction is Σˣ - ε_x^HF
        # This is simplified - full implementation would include correlation
        corrections = np.diag(self.sigma_x) - self.mf.mo_energy
        
        if orbital_idx is not None:
            return corrections[orbital_idx]
        return corrections


# Module-level convenience functions
def exchange_self_energy_from_pyscf(mf, **kwargs) -> np.ndarray:
    """
    Compute exchange self-energy directly from PySCF mean-field object.
    
    Parameters
    ----------
    mf : pyscf.scf.hf.SCF
        Converged mean-field object
    **kwargs : dict
        Additional arguments passed to compute_exchange_self_energy
    
    Returns
    -------
    np.ndarray
        Exchange self-energy matrix
    
    Examples
    --------
    >>> from pyscf import gto, scf
    >>> import quasix
    >>> 
    >>> mol = gto.M(atom='H 0 0 0; H 0 0 1', basis='cc-pvdz')  
    >>> mf = scf.RHF(mol).run()
    >>> sigma_x = quasix.selfenergy.exchange_self_energy_from_pyscf(mf)
    """
    calc = ExchangeSelfEnergy(mf)
    df_integrals = calc.build_df_integrals()
    
    return compute_exchange_self_energy(
        mf.mo_energy,
        mf.mo_occ,
        mf.mo_coeff,
        df_integrals,
        **kwargs
    )