"""
Density Fitting Tensor Management for QuasiX

This module provides the DFTensor class for extracting and managing
density fitting tensors from PySCF mean-field objects.
"""

import numpy as np
import h5py
from typing import Optional, Dict, Any, Union, Tuple
from pyscf import gto, scf, df, lib
from pyscf.scf import RHF, RKS
from pyscf.lib import logger as pyscf_logger
import logging

from .mo_transform import (
    build_transition_df_tensor,
    build_full_df_tensors,
    transform_mo_coefficients,
    validate_mo_coefficients
)
from .logging import get_logger, timed_stage

# Set up module logger
logger = get_logger(__name__)


class DFTensor:
    """
    Density Fitting Tensor container for QuasiX
    
    This class manages the extraction, storage, and access to density fitting
    tensors from PySCF calculations. It provides a clean interface between
    PySCF's density fitting infrastructure and QuasiX's Rust backend.
    
    Attributes:
        mol: PySCF molecule object
        auxmol: Auxiliary molecule for density fitting
        mo_coeff: Molecular orbital coefficients
        mo_occ: Molecular orbital occupations
        mo_energy: Molecular orbital energies
        nocc: Number of occupied orbitals
        nvir: Number of virtual orbitals
        naux: Number of auxiliary basis functions
        iaP: (ia|P) transition density fitting tensor
        ijP: (ij|P) occupied-occupied DF tensor
        abP: (ab|P) virtual-virtual DF tensor
        metric: Two-center metric matrix (P|Q)
        chol_metric: Cholesky factor of metric matrix
    """
    
    def __init__(self, mol: gto.Mole, auxbasis: str = 'def2-tzvp-jkfit'):
        """
        Initialize DFTensor with molecule and auxiliary basis
        
        Args:
            mol: PySCF molecule object
            auxbasis: Auxiliary basis set for density fitting
        """
        self.mol = mol
        self.auxbasis = auxbasis
        self.auxmol = df.addons.make_auxmol(mol, auxbasis)
        
        # MO data (to be filled by from_pyscf)
        self.mo_coeff = None
        self.mo_occ = None
        self.mo_energy = None
        
        # Dimensions
        self.nao = mol.nao
        self.naux = self.auxmol.nao
        self.nocc = None
        self.nvir = None
        
        # DF tensors
        self.iaP = None  # (ia|P) transition tensor
        self.ijP = None  # (ij|P) occupied-occupied tensor [DEPRECATED - use miP]
        self.miP = None  # (mi|P) all-occupied tensor [NEW - for Σˣ]
        self.abP = None  # (ab|P) virtual-virtual tensor
        self.cderi = None  # Cholesky-decomposed 3-center integrals (naux, nao*(nao+1)//2)

        # Metric tensors
        self.metric = None  # (P|Q) two-center integrals
        self.chol_metric = None  # Cholesky factor of metric
        
        # Memory tracking
        self._memory_usage = 0
        
        logger.info(f"DFTensor initialized: nao={self.nao}, naux={self.naux}")
    
    @classmethod
    def from_pyscf(cls, mf: Union[RHF, RKS], 
                   auxbasis: str = 'def2-tzvp-jkfit',
                   build_all: bool = True) -> 'DFTensor':
        """
        Extract DF tensors from a converged PySCF mean-field object
        
        This is the primary constructor for creating DFTensor objects from
        PySCF calculations. It extracts MO data and builds the necessary
        density fitting tensors.
        
        Args:
            mf: Converged PySCF mean-field object (RHF or RKS)
            auxbasis: Auxiliary basis set for density fitting
            build_all: If True, build all DF tensors (iaP, ijP, abP)
                      If False, only build iaP (for memory efficiency)
        
        Returns:
            DFTensor object with populated DF tensors
        
        Raises:
            ValueError: If SCF is not converged or MO data is invalid
        """
        if not mf.converged:
            logger.warning("SCF not converged! Results may be unreliable.")
        
        # Create instance
        df_tensor = cls(mf.mol, auxbasis)
        
        # Extract MO data
        df_tensor.mo_coeff = mf.mo_coeff.copy()
        df_tensor.mo_occ = mf.mo_occ.copy()
        df_tensor.mo_energy = mf.mo_energy.copy()
        
        # Validate MO coefficients
        try:
            S = mf.mol.intor('int1e_ovlp')
            validate_mo_coefficients(df_tensor.mo_coeff, df_tensor.nao, 
                                    df_tensor.mo_occ, overlap_matrix=S)
        except ValueError as e:
            logger.error(f"MO validation failed: {e}")
            raise
        
        # Set dimensions
        df_tensor.nocc = int(np.sum(df_tensor.mo_occ > 0))
        df_tensor.nvir = df_tensor.nao - df_tensor.nocc
        
        logger.info(f"Extracted MO data: nocc={df_tensor.nocc}, nvir={df_tensor.nvir}")
        
        # Build DF tensors
        with timed_stage("Building DF tensors"):
            if build_all:
                # Build all tensors
                tensors = build_full_df_tensors(
                    df_tensor.mol,
                    df_tensor.auxmol,
                    df_tensor.mo_coeff,
                    df_tensor.mo_occ
                )
                df_tensor.iaP = tensors['iaP']
                df_tensor.miP = tensors['miP']  # NEW: Full (mi|P) tensor
                df_tensor.ijP = tensors['ijP']  # Backward compatibility (subset of miP)
                df_tensor.abP = tensors['abP']

                # Extract metric if available
                if 'metric' in tensors:
                    df_tensor.metric = tensors['metric']

                # Extract Cholesky vectors if available
                if 'cderi' in tensors:
                    df_tensor.cderi = tensors['cderi']

                logger.info(f"Built all DF tensors: iaP {df_tensor.iaP.shape}, "
                           f"miP {df_tensor.miP.shape}, ijP {df_tensor.ijP.shape}, "
                           f"abP {df_tensor.abP.shape}")
            else:
                # Build only transition tensor for memory efficiency
                df_tensor.iaP = build_transition_df_tensor(
                    df_tensor.mol,
                    df_tensor.auxmol,
                    df_tensor.mo_coeff,
                    df_tensor.mo_occ
                )
                logger.info(f"Built transition tensor: iaP {df_tensor.iaP.shape}")
        
        # Build metric tensor
        df_tensor._build_metric()
        
        # Calculate memory usage
        df_tensor._calculate_memory_usage()
        
        return df_tensor
    
    def _build_metric(self):
        """Build two-center metric matrix (P|Q) and its Cholesky factor"""
        with timed_stage("Building metric tensor"):
            # Compute (P|Q) two-center integrals
            self.metric = self.auxmol.intor('int2c2e')
            
            # Compute Cholesky factorization for metric^(1/2)
            try:
                self.chol_metric = np.linalg.cholesky(self.metric)
                logger.info(f"Built metric tensor: {self.metric.shape}")
            except np.linalg.LinAlgError:
                logger.warning("Metric not positive definite, using eigendecomposition")
                # Use eigendecomposition for metric^(1/2)
                w, v = np.linalg.eigh(self.metric)
                # Remove negative eigenvalues (numerical noise)
                idx = w > 1e-10
                w = w[idx]
                v = v[:, idx]
                self.chol_metric = v @ np.diag(np.sqrt(w))
    
    def get_transition_tensor(self) -> np.ndarray:
        """
        Get the (ia|P) transition density fitting tensor
        
        Returns:
            ndarray of shape (nocc*nvir, naux)
        """
        if self.iaP is None:
            raise ValueError("Transition tensor not built. Call from_pyscf() first.")
        return self.iaP
    
    def get_occupied_tensor(self) -> Optional[np.ndarray]:
        """
        Get the (ij|P) occupied-occupied density fitting tensor
        
        Returns:
            ndarray of shape (nocc*nocc, naux) or None if not built
        """
        return self.ijP
    
    def get_virtual_tensor(self) -> Optional[np.ndarray]:
        """
        Get the (ab|P) virtual-virtual density fitting tensor

        Returns:
            ndarray of shape (nvir*nvir, naux) or None if not built
        """
        return self.abP

    def get_cholesky_vectors(self) -> Optional[np.ndarray]:
        """
        Get the Cholesky-decomposed 3-center integrals (AO basis).

        These are the L vectors such that:
            (μν|ρσ) ≈ Σ_P L_μν,P * L_ρσ,P

        Returns:
            ndarray of shape (naux, nao*(nao+1)//2) in packed triangular format,
            or None if not built
        """
        return self.cderi

    def get_metric(self) -> np.ndarray:
        """
        Get the two-center metric matrix (P|Q)
        
        Returns:
            ndarray of shape (naux, naux)
        """
        if self.metric is None:
            self._build_metric()
        return self.metric
    
    def get_cholesky_metric(self) -> np.ndarray:
        """
        Get the Cholesky factor of the metric matrix
        
        Returns:
            ndarray of shape (naux, naux) such that L @ L.T = metric
        """
        if self.chol_metric is None:
            self._build_metric()
        return self.chol_metric
    
    def apply_metric_sqrt(self, x: np.ndarray) -> np.ndarray:
        """
        Apply metric^(1/2) to a vector or matrix
        
        Args:
            x: Input array of shape (..., naux)
        
        Returns:
            metric^(1/2) @ x
        """
        L = self.get_cholesky_metric()
        return x @ L.T
    
    def apply_metric_inv_sqrt(self, x: np.ndarray) -> np.ndarray:
        """
        Apply metric^(-1/2) to a vector or matrix
        
        Args:
            x: Input array of shape (..., naux)
        
        Returns:
            metric^(-1/2) @ x
        """
        L = self.get_cholesky_metric()
        # Solve L.T @ y = x.T for y
        y = np.linalg.solve(L.T, x.T)
        return y.T
    
    def save_to_hdf5(self, filename: str, use_rust: bool = True):
        """
        Save DF tensors to HDF5 file

        Args:
            filename: Path to output HDF5 file
            use_rust: Use Rust HDF5 implementation if available (default: True)
        """
        # Try Rust backend first if requested
        if use_rust:
            try:
                from .io import save_df_tensors_hdf5, HAS_RUST_HDF5
                if HAS_RUST_HDF5:
                    save_df_tensors_hdf5(self, filename)
                    logger.info(f"Saved DF tensors to {filename} using Rust backend")
                    return
            except (ImportError, AttributeError) as e:
                logger.warning(f"Rust HDF5 backend not available: {e}. Falling back to Python.")

        # Fallback to Python h5py implementation
        with h5py.File(filename, 'w') as f:
            # Save metadata
            f.attrs['auxbasis'] = self.auxbasis
            f.attrs['nao'] = self.nao
            f.attrs['naux'] = self.naux
            f.attrs['nocc'] = self.nocc
            f.attrs['nvir'] = self.nvir

            # Save MO data
            f.create_dataset('mo_coeff', data=self.mo_coeff)
            f.create_dataset('mo_occ', data=self.mo_occ)
            f.create_dataset('mo_energy', data=self.mo_energy)

            # Save DF tensors
            if self.iaP is not None:
                f.create_dataset('iaP', data=self.iaP, chunks=(min(1000, self.iaP.shape[0]), self.naux),
                               compression='gzip', compression_opts=4)
            if self.ijP is not None:
                f.create_dataset('ijP', data=self.ijP, chunks=(min(100, self.ijP.shape[0]), self.naux),
                               compression='gzip', compression_opts=4)
            if self.abP is not None:
                f.create_dataset('abP', data=self.abP, chunks=(min(1000, self.abP.shape[0]), self.naux),
                               compression='gzip', compression_opts=4)

            # Save metric
            if self.metric is not None:
                f.create_dataset('metric', data=self.metric)
            if self.chol_metric is not None:
                f.create_dataset('chol_metric', data=self.chol_metric)

            logger.info(f"Saved DF tensors to {filename} using Python backend")
    
    @classmethod
    def load_from_hdf5(cls, filename: str, mol: gto.Mole, use_rust: bool = True) -> 'DFTensor':
        """
        Load DF tensors from HDF5 file

        Args:
            filename: Path to HDF5 file
            mol: PySCF molecule object
            use_rust: Use Rust HDF5 implementation if available (default: True)

        Returns:
            DFTensor object with loaded data
        """
        # Try Rust backend first if requested
        if use_rust:
            try:
                from .io import load_df_tensors_hdf5, convert_to_dftensor, HAS_RUST_HDF5
                if HAS_RUST_HDF5:
                    loaded_data = load_df_tensors_hdf5(filename)
                    df_tensor = convert_to_dftensor(loaded_data, mol)
                    logger.info(f"Loaded DF tensors from {filename} using Rust backend")
                    df_tensor._calculate_memory_usage()
                    return df_tensor
            except (ImportError, AttributeError) as e:
                logger.warning(f"Rust HDF5 backend not available: {e}. Falling back to Python.")

        # Fallback to Python h5py implementation
        with h5py.File(filename, 'r') as f:
            # Create instance
            auxbasis = f.attrs['auxbasis']
            df_tensor = cls(mol, auxbasis)

            # Load dimensions
            df_tensor.nocc = f.attrs['nocc']
            df_tensor.nvir = f.attrs['nvir']

            # Load MO data
            df_tensor.mo_coeff = f['mo_coeff'][:]
            df_tensor.mo_occ = f['mo_occ'][:]
            df_tensor.mo_energy = f['mo_energy'][:]

            # Load DF tensors
            if 'iaP' in f:
                df_tensor.iaP = f['iaP'][:]
            if 'ijP' in f:
                df_tensor.ijP = f['ijP'][:]
            if 'abP' in f:
                df_tensor.abP = f['abP'][:]

            # Load metric
            if 'metric' in f:
                df_tensor.metric = f['metric'][:]
            if 'chol_metric' in f:
                df_tensor.chol_metric = f['chol_metric'][:]

            logger.info(f"Loaded DF tensors from {filename} using Python backend")

        df_tensor._calculate_memory_usage()
        return df_tensor
    
    def _calculate_memory_usage(self):
        """Calculate memory usage of stored tensors"""
        self._memory_usage = 0
        
        # Count tensor memory
        for tensor_name in ['iaP', 'ijP', 'abP', 'metric', 'chol_metric', 
                           'mo_coeff', 'mo_occ', 'mo_energy']:
            tensor = getattr(self, tensor_name)
            if tensor is not None:
                self._memory_usage += tensor.nbytes
        
        # Convert to MB
        self._memory_usage = self._memory_usage / (1024 * 1024)
        logger.info(f"Total memory usage: {self._memory_usage:.2f} MB")
    
    @property
    def memory_usage_mb(self) -> float:
        """Get current memory usage in MB"""
        return self._memory_usage
    
    def __repr__(self) -> str:
        """String representation"""
        return (f"DFTensor(auxbasis='{self.auxbasis}', "
                f"nocc={self.nocc}, nvir={self.nvir}, naux={self.naux}, "
                f"memory={self._memory_usage:.2f} MB)")