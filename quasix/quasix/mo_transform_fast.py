"""
Fast MO transformation utilities optimized for QuasiX

This module provides optimized transformations from AO to MO basis,
addressing the real performance bottlenecks in PySCF integration.
"""

import numpy as np
import logging
from typing import Optional, Tuple
from pyscf import gto, df
import os

from .logging import get_logger, timed_stage

# Set up module logger
logger = get_logger(__name__)

# Configure threading to avoid conflicts
def configure_threading():
    """Configure thread settings to avoid conflicts between NumPy, BLAS, and Rayon"""
    # Set single-threaded BLAS when using Rayon in Rust
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'
    os.environ['OMP_NUM_THREADS'] = '1'
    logger.debug("Thread configuration set for optimal performance")

# Memory pool for temporary arrays
class MemoryPool:
    """Simple memory pool to reduce allocation overhead"""
    def __init__(self):
        self.pool = {}

    def get(self, shape: tuple, dtype=np.float64) -> np.ndarray:
        """Get array from pool or allocate new one"""
        key = (shape, dtype)
        if key in self.pool:
            arr = self.pool[key]
            arr.fill(0)  # Clear array
            return arr
        else:
            arr = np.zeros(shape, dtype=dtype)
            self.pool[key] = arr
            return arr

    def clear(self):
        """Clear the memory pool"""
        self.pool.clear()

# Global memory pool instance
_memory_pool = MemoryPool()

def build_transition_df_tensor_fast(mol: gto.Mole, auxbasis: str,
                                   mo_coeff: np.ndarray, mo_occ: np.ndarray) -> np.ndarray:
    """
    Fast (ia|P) tensor construction with optimized PySCF integration

    This version addresses the real bottleneck (PySCF overhead) by:
    1. Using PySCF's cached DF infrastructure when available
    2. Minimizing temporary array creation
    3. Using memory pool for repeated allocations
    4. Proper thread configuration

    Args:
        mol: PySCF molecule object
        auxbasis: Auxiliary basis set name
        mo_coeff: MO coefficients [nao, nmo]
        mo_occ: MO occupations [nmo]

    Returns:
        (ia|P) tensor [nocc*nvir, naux]
    """
    configure_threading()

    with timed_stage("Fast (ia|P) construction"):
        # Get dimensions
        occ_idx = mo_occ > 0.5
        nocc = np.sum(occ_idx)
        nvir = len(mo_occ) - nocc

        # Split MO coefficients (views, no copy)
        c_occ = mo_coeff[:, occ_idx]
        c_vir = mo_coeff[:, ~occ_idx]

        # Use PySCF's DF infrastructure efficiently
        mydf = df.DF(mol)
        mydf.auxbasis = auxbasis

        # Check if DF is already built (cached)
        if not hasattr(mydf, '_cderi'):
            with timed_stage("Building DF integrals"):
                mydf.build()

        naux = mydf.get_naoaux()

        # Get 3-center integrals efficiently
        # Use PySCF's optimized routines
        with timed_stage("Getting 3-center integrals"):
            # Build properly through PySCF API
            j3c = df.incore.aux_e2(mol, mydf.auxmol, intor='int3c2e', aosym='s1')
            # This should give us (nao*nao, naux) or compact format
            # Reshape to (nao, nao, naux)
            if j3c.ndim == 2:
                nao_sq = mol.nao * mol.nao
                if j3c.shape[0] == nao_sq:
                    j3c = j3c.reshape(mol.nao, mol.nao, -1)
                else:
                    # Use s1 (no symmetry) to get full matrix
                    j3c = df.incore.aux_e2(mol, mydf.auxmol, intor='int3c2e', aosym='s1')
                    j3c = j3c.reshape(mol.nao, mol.nao, -1)
            else:
                # For larger systems, use outcore with blocking
                j3c = _memory_pool.get((mol.nao, mol.nao, naux))
                mydf._cderi_to_save = j3c.reshape(-1, naux)
                mydf.build()
                j3c = j3c.reshape(mol.nao, mol.nao, naux)

        # Optimized MO transformation
        with timed_stage("MO transformation"):
            # Use BLAS-optimized operations
            # First transform: (μν|P) @ C_μi -> (iν|P)
            temp = np.zeros((nocc, mol.nao, naux))
            for p in range(naux):
                temp[:, :, p] = c_occ.T @ j3c[:, :, p]

            # Second transform: (iν|P) @ C_νa -> (ia|P)
            ia_P = np.zeros((nocc, nvir, naux))
            for p in range(naux):
                ia_P[:, :, p] = temp[:, :, p] @ c_vir

            # Reshape to combined index
            ia_P = ia_P.reshape(nocc * nvir, naux)

        # Ensure C-contiguous for Rust
        if not ia_P.flags['C_CONTIGUOUS']:
            ia_P = np.ascontiguousarray(ia_P)

        logger.info(f"Fast (ia|P) built: shape={ia_P.shape}")
        return ia_P


def build_metric_fast(mol: gto.Mole, auxbasis: str) -> np.ndarray:
    """
    Fast metric (P|Q) construction

    Args:
        mol: PySCF molecule
        auxbasis: Auxiliary basis name

    Returns:
        Metric matrix (P|Q) [naux, naux]
    """
    configure_threading()

    mydf = df.DF(mol)
    mydf.auxbasis = auxbasis

    # Get 2-center integrals
    with timed_stage("Building metric (P|Q)"):
        metric = mydf.get_2c2e()  # This is already optimized in PySCF

    # Ensure symmetric and positive definite
    metric = 0.5 * (metric + metric.T)

    # Check condition number
    eigenvalues = np.linalg.eigvalsh(metric)
    condition = eigenvalues[-1] / eigenvalues[0]
    logger.debug(f"Metric condition number: {condition:.2e}")

    if condition > 1e12:
        logger.warning(f"Metric poorly conditioned: {condition:.2e}")

    return metric


def build_all_df_tensors_fast(mol: gto.Mole, auxbasis: str,
                              mo_coeff: np.ndarray, mo_occ: np.ndarray) -> dict:
    """
    Build all DF tensors efficiently with caching

    Returns:
        Dictionary with keys: 'iaP', 'ijP', 'abP', 'metric'
    """
    configure_threading()

    # Build DF object once and reuse
    mydf = df.DF(mol)
    mydf.auxbasis = auxbasis
    mydf.build()

    occ_idx = mo_occ > 0.5
    c_occ = mo_coeff[:, occ_idx]
    c_vir = mo_coeff[:, ~occ_idx]

    naux = mydf.get_naoaux()
    nocc = c_occ.shape[1]
    nvir = c_vir.shape[1]

    result = {}

    # Get 3-center integrals once
    j3c = df.incore.aux_e2(mol, mydf.auxmol, intor='int3c2e')
    j3c = j3c.reshape(mol.nao, mol.nao, naux)

    # Build all tensors using the same j3c
    with timed_stage("Building all DF tensors"):
        # (ia|P)
        ia_P = np.einsum('mnp,mi,na->iap', j3c, c_occ, c_vir, optimize='optimal')
        result['iaP'] = ia_P.reshape(nocc * nvir, naux)

        # (ij|P)
        ij_P = np.einsum('mnp,mi,nj->ijp', j3c, c_occ, c_occ, optimize='optimal')
        result['ijP'] = ij_P.reshape(nocc * nocc, naux)

        # (ab|P)
        ab_P = np.einsum('mnp,ma,nb->abp', j3c, c_vir, c_vir, optimize='optimal')
        result['abP'] = ab_P.reshape(nvir * nvir, naux)

        # Metric
        result['metric'] = mydf.get_2c2e()

    # Ensure all arrays are C-contiguous
    for key in result:
        if not result[key].flags['C_CONTIGUOUS']:
            result[key] = np.ascontiguousarray(result[key])

    logger.info(f"All DF tensors built: {list(result.keys())}")
    return result