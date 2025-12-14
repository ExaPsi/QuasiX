"""
MO transformation utilities for QuasiX

This module provides efficient transformations from AO to MO basis,
specifically optimized for building (ia|P) tensors for GW calculations.
"""

import numpy as np
import logging
from typing import Optional, Tuple, Union
from pyscf import gto, df, lib
from pyscf.lib import unpack_tril
import h5py
import scipy.linalg

from .logging import get_logger, timed_stage

# Set up module logger
logger = get_logger(__name__)


def validate_mo_coefficients(mo_coeff: np.ndarray, nao: int, 
                            mo_occ: np.ndarray, overlap_matrix: np.ndarray = None) -> None:
    """
    Validate MO coefficient matrix
    
    Args:
        mo_coeff: MO coefficient matrix
        nao: Number of atomic orbitals
        mo_occ: MO occupations
        overlap_matrix: AO overlap matrix S (if None, assumes orthonormal AOs)
        
    Raises:
        ValueError: If MO coefficients are invalid
    """
    if mo_coeff.ndim != 2:
        raise ValueError(f"MO coefficients must be 2D, got shape {mo_coeff.shape}")
    
    if mo_coeff.shape[0] != nao:
        raise ValueError(f"MO coefficient shape {mo_coeff.shape} inconsistent with nao={nao}")
    
    if mo_coeff.shape[1] != len(mo_occ):
        raise ValueError(f"MO coefficient shape {mo_coeff.shape} inconsistent with mo_occ length {len(mo_occ)}")
    
    # Check for orthonormality (within tolerance)
    # MO coefficients should satisfy: C^T @ S @ C = I
    if overlap_matrix is not None:
        overlap = mo_coeff.T @ overlap_matrix @ mo_coeff
    else:
        overlap = mo_coeff.T @ mo_coeff
    identity = np.eye(mo_coeff.shape[1])
    ortho_error = np.max(np.abs(overlap - identity))
    
    if ortho_error > 1e-6:
        logger.warning(f"MO coefficients not perfectly orthonormal (max error: {ortho_error:.2e})")
        if ortho_error > 1e-3:
            raise ValueError(f"MO coefficients poorly orthonormal (max error: {ortho_error:.2e})")


def transform_mo_coefficients(mo_coeff: np.ndarray, mo_occ: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Split MO coefficients into occupied and virtual blocks

    Args:
        mo_coeff: Full MO coefficient matrix (nao, nmo)
        mo_occ: MO occupations

    Returns:
        c_occ: Occupied MO coefficients (nao, nocc)
        c_vir: Virtual MO coefficients (nao, nvir)
    """
    occ_mask = mo_occ > 0
    vir_mask = ~occ_mask

    c_occ = mo_coeff[:, occ_mask]
    c_vir = mo_coeff[:, vir_mask]

    logger.debug(f"Split MOs: nocc={c_occ.shape[1]}, nvir={c_vir.shape[1]}")

    return c_occ, c_vir


def apply_cholesky_decomposition(mol: gto.Mole, auxmol: gto.Mole) -> np.ndarray:
    """
    Compute Cholesky-decomposed 3-center integrals following PySCF convention.

    PySCF uses Cholesky decomposition for numerical stability and efficiency:
        L = low^(-1) @ j3c
    where:
        low @ low.T = j2c (2-center metric)
        j3c = raw 3-center integrals (μν|P)
        L = Cholesky-decomposed integrals

    This is equivalent to: L = j2c^(-1/2) @ j3c

    The Cholesky form satisfies:
        (μν|ρσ) ≈ Σ_P L_μν,P L_ρσ,P

    Args:
        mol: PySCF molecule object
        auxmol: Auxiliary basis molecule

    Returns:
        cderi: Cholesky-decomposed integrals, shape (naux, nao*(nao+1)//2)

    References:
        - PySCF implementation: pyscf/df/incore.py:cholesky_eri()
        - Theory: Aquilante et al., J. Chem. Phys. 127, 114107 (2007)
    """
    with timed_stage("apply_cholesky_decomposition"):
        nao = mol.nao
        naux = auxmol.nao

        # Step 1: Compute 2-center metric (P|Q)
        j2c = auxmol.intor('int2c2e')
        logger.debug(f"2-center metric (P|Q): shape={j2c.shape}, norm={np.linalg.norm(j2c):.6f}")

        # Step 2: Cholesky decomposition of metric
        # j2c = low @ low.T
        try:
            low = scipy.linalg.cholesky(j2c, lower=True)
        except scipy.linalg.LinAlgError:
            # Fallback to eigendecomposition if Cholesky fails
            logger.warning("Cholesky decomposition failed, using eigendecomposition")
            w, v = np.linalg.eigh(j2c)
            # Remove small/negative eigenvalues
            idx = w > 1e-10
            w = w[idx]
            v = v[:, idx]
            low = v @ np.diag(np.sqrt(w))

        logger.debug(f"Cholesky factor: shape={low.shape}")

        # Step 3: Compute raw 3-center integrals in packed format
        # Use aosym='s2ij' for packed triangular storage (matches PySCF)
        j3c_packed = df.incore.aux_e2(mol, auxmol, intor='int3c2e', aosym='s2ij')
        # Shape: (nao*(nao+1)//2, naux)

        logger.debug(f"Raw j3c (packed): shape={j3c_packed.shape}, norm={np.linalg.norm(j3c_packed):.6f}")

        # Step 4: Apply Cholesky transformation
        # Solve: low @ cderi = j3c.T for cderi
        # This gives: cderi = low^(-1) @ j3c.T
        cderi = scipy.linalg.solve_triangular(low, j3c_packed.T, lower=True,
                                              check_finite=False, overwrite_b=False)
        # Shape: (naux, nao_pair)

        logger.debug(f"Cholesky-decomposed L: shape={cderi.shape}, norm={np.linalg.norm(cderi):.6f}")

        # Validate result
        if np.any(np.isnan(cderi)) or np.any(np.isinf(cderi)):
            raise RuntimeError("Cholesky decomposition produced NaN or Inf values")

        return cderi


def build_transition_df_tensor(mol: gto.Mole, auxmol: gto.Mole,
                              mo_coeff: np.ndarray, mo_occ: np.ndarray,
                              max_memory: float = 4000,
                              block_size: Optional[int] = None) -> np.ndarray:
    """
    Build (ia|P) transition density fitting tensor with Cholesky decomposition.

    This function constructs Cholesky-decomposed 3-center integrals in the
    occupied-virtual MO basis, following PySCF's DF convention for numerical
    stability and compatibility.

    The returned tensor satisfies:
        (ia|jb) ≈ Σ_P iaP[ia,P] * iaP[jb,P]

    Args:
        mol: PySCF molecule object
        auxmol: Auxiliary basis molecule
        mo_coeff: MO coefficients
        mo_occ: MO occupations
        max_memory: Maximum memory in MB
        block_size: Block size for auxiliary index (auto-determined if None)

    Returns:
        iaP: Cholesky-decomposed transition DF tensor of shape (nocc*nvir, naux)

    Note:
        Uses Cholesky decomposition of 2-center metric for numerical stability.
        Matches PySCF convention: L = j2c^(-1/2) @ j3c
    """
    with timed_stage("build_transition_df"):
        nao = mol.nao
        naux = auxmol.nao
        nocc = np.sum(mo_occ > 0)
        nvir = nao - nocc

        logger.info(f"Building (ia|P) tensor: nocc={nocc}, nvir={nvir}, naux={naux}")

        # Step 1: Get Cholesky-decomposed 3-center integrals
        # This matches PySCF's _cderi convention
        cderi = apply_cholesky_decomposition(mol, auxmol)
        # Shape: (naux, nao*(nao+1)//2)

        # Step 2: Transform to MO basis
        # For each auxiliary function P, transform L_μν,P to L_pq,P
        iaP = np.zeros((nocc, nvir, naux), dtype=np.float64)

        for P in range(naux):
            # Unpack triangular Cholesky vector to full (nao, nao)
            Lp_ao = unpack_tril(cderi[P])

            # Transform to MO basis: L_pq = C^T @ L_μν @ C
            Lp_mo = mo_coeff.T @ Lp_ao @ mo_coeff

            # Extract occupied-virtual block
            iaP[:, :, P] = Lp_mo[:nocc, nocc:]

        # Reshape to (nocc*nvir, naux)
        iaP_flat = iaP.reshape(nocc * nvir, naux)

        # Verify tensor properties
        logger.info(f"(ia|P) tensor built: shape={iaP_flat.shape}, memory={iaP_flat.nbytes/1e9:.2f} GB")
        logger.debug(f"(ia|P) norm: {np.linalg.norm(iaP_flat):.6f}")

        return iaP_flat


def build_full_df_tensors(mol: gto.Mole, auxmol: gto.Mole,
                         mo_coeff: np.ndarray, mo_occ: np.ndarray,
                         max_memory: float = 4000) -> dict:
    """
    Build all Cholesky-decomposed DF tensors needed for GW: (ia|P), (mi|P), (ab|P).

    CRITICAL FIX (2025-11-18): Changed from (ij|P) to (mi|P) for exchange self-energy.
    The exchange self-energy Σˣ requires summing over ALL orbitals m (occupied + virtual),
    not just occupied orbitals i,j. The tensor (mi|P) has first index m ∈ [0,nmo) and
    second index i ∈ [0,nocc).

    All tensors use PySCF's Cholesky DF convention for numerical stability:
        L = j2c^(-1/2) @ j3c

    This ensures compatibility with PySCF and proper DF approximation:
        (pq|rs) ≈ Σ_P L_pq,P * L_rs,P

    Args:
        mol: PySCF molecule object
        auxmol: Auxiliary basis molecule
        mo_coeff: MO coefficients
        mo_occ: MO occupations
        max_memory: Maximum memory in MB

    Returns:
        Dictionary containing:
            'iaP': Cholesky-decomposed (ia|P) tensor, shape (nocc*nvir, naux)
            'miP': Cholesky-decomposed (mi|P) tensor, shape (nmo*nocc, naux) [FULL tensor for Σˣ]
            'abP': Cholesky-decomposed (ab|P) tensor, shape (nvir*nvir, naux)
            'metric': 2-center metric (P|Q), shape (naux, naux)
            'nocc': Number of occupied orbitals
            'nvir': Number of virtual orbitals
            'naux': Number of auxiliary functions
    """
    with timed_stage("build_full_df_tensors"):
        nao = mol.nao
        naux = auxmol.nao
        nocc = np.sum(mo_occ > 0)
        nvir = nao - nocc
        nmo = nao  # Total number of MOs

        logger.info(f"Building full DF tensor set for GW (FIXED: using (mi|P) instead of (ij|P))")

        # Estimate memory requirements
        mem_ia = nocc * nvir * naux * 8 / 1e6  # MB
        mem_mi = nmo * nocc * naux * 8 / 1e6   # MB - FULL tensor for Σˣ
        mem_ab = nvir * nvir * naux * 8 / 1e6  # MB
        total_mem = mem_ia + mem_mi + mem_ab

        logger.info(f"Memory estimate: (ia|P)={mem_ia:.1f} MB, (mi|P)={mem_mi:.1f} MB, (ab|P)={mem_ab:.1f} MB")

        if total_mem > max_memory:
            logger.warning(f"Total memory ({total_mem:.1f} MB) exceeds limit ({max_memory:.1f} MB)")
            logger.warning("Consider using out-of-core storage or reducing basis size")

        # Step 1: Get Cholesky-decomposed 3-center integrals
        cderi = apply_cholesky_decomposition(mol, auxmol)
        # Shape: (naux, nao*(nao+1)//2)

        # Step 2: Transform to MO basis for all blocks
        iaP = np.zeros((nocc, nvir, naux), dtype=np.float64)
        miP = np.zeros((nmo, nocc, naux), dtype=np.float64)  # FULL tensor: ALL orbitals × occupied
        abP = np.zeros((nvir, nvir, naux), dtype=np.float64)

        for P in range(naux):
            # Unpack Cholesky vector for auxiliary function P
            Lp_ao = unpack_tril(cderi[P])  # (nao, nao)

            # Transform to MO basis: L_pq = C^T @ L_μν @ C
            Lp_mo = mo_coeff.T @ Lp_ao @ mo_coeff  # (nmo, nmo)

            # Extract blocks
            iaP[:, :, P] = Lp_mo[:nocc, nocc:]     # occupied-virtual
            miP[:, :, P] = Lp_mo[:, :nocc]         # ALL-occupied (for Σˣ)
            abP[:, :, P] = Lp_mo[nocc:, nocc:]     # virtual-virtual

        # Reshape to flattened format
        iaP_flat = iaP.reshape(nocc * nvir, naux)
        miP_flat = miP.reshape(nmo * nocc, naux)  # Shape: (nmo*nocc, naux)
        abP_flat = abP.reshape(nvir * nvir, naux)

        # Get 2-center metric
        metric = auxmol.intor('int2c2e')

        result = {
            'iaP': iaP_flat,
            'miP': miP_flat,  # CHANGED: (mi|P) instead of (ij|P)
            'abP': abP_flat,
            'metric': metric,
            'cderi': cderi,  # Store Cholesky vectors for validation
            'nocc': nocc,
            'nvir': nvir,
            'naux': naux,
            'nmo': nmo       # Add nmo for clarity
        }

        logger.info("DF tensor construction completed")
        logger.debug(f"iaP norm: {np.linalg.norm(iaP_flat):.6f}")
        logger.debug(f"miP norm: {np.linalg.norm(miP_flat):.6f}")
        logger.debug(f"abP norm: {np.linalg.norm(abP_flat):.6f}")

        # Backward compatibility: also provide ijP as subset of miP for old code
        ijP_flat = miP.reshape(nmo, nocc, naux)[:nocc, :, :].reshape(nocc * nocc, naux)
        result['ijP'] = ijP_flat
        logger.debug(f"ijP (subset) norm: {np.linalg.norm(ijP_flat):.6f}")

        return result


def _transform_3c_to_mo_pairs(j3c: np.ndarray, c_left: np.ndarray, 
                              c_right: np.ndarray) -> np.ndarray:
    """
    Transform 3-center AO integrals to MO pair basis
    
    Args:
        j3c: 3-center integrals (naux_block, nao, nao)
        c_left: Left MO coefficients (nao, nmo_left)
        c_right: Right MO coefficients (nao, nmo_right)
        
    Returns:
        j3c_mo: Transformed integrals (naux_block, nmo_left, nmo_right)
    """
    naux_block = j3c.shape[0]
    nmo_left = c_left.shape[1]
    nmo_right = c_right.shape[1]
    
    # Use einsum for efficient transformation
    # (P,μν) @ (μ,i) @ (ν,j) = (P,ij)
    j3c_mo = np.zeros((naux_block, nmo_left, nmo_right))
    
    for iaux in range(naux_block):
        # Two-step transformation for better performance
        temp = j3c[iaux] @ c_right  # (nao, nmo_right)
        j3c_mo[iaux] = c_left.T @ temp  # (nmo_left, nmo_right)
    
    return j3c_mo


def save_df_tensors_hdf5(filename: str, df_tensors: dict, 
                         compression: str = 'gzip', compression_opts: int = 4):
    """
    Save DF tensors to HDF5 file for later use
    
    Args:
        filename: Output HDF5 file
        df_tensors: Dictionary from build_full_df_tensors
        compression: HDF5 compression algorithm
        compression_opts: Compression level (1-9)
    """
    with h5py.File(filename, 'w') as f:
        # Save tensors with compression
        for key in ['iaP', 'ijP', 'abP']:
            if key in df_tensors:
                f.create_dataset(key, data=df_tensors[key],
                               compression=compression,
                               compression_opts=compression_opts)
        
        # Save dimensions as attributes
        f.attrs['nocc'] = df_tensors['nocc']
        f.attrs['nvir'] = df_tensors['nvir']
        f.attrs['naux'] = df_tensors['naux']
        
    logger.info(f"DF tensors saved to {filename}")


def load_df_tensors_hdf5(filename: str) -> dict:
    """
    Load DF tensors from HDF5 file
    
    Args:
        filename: Input HDF5 file
        
    Returns:
        Dictionary with DF tensors and dimensions
    """
    result = {}
    
    with h5py.File(filename, 'r') as f:
        # Load tensors
        for key in ['iaP', 'ijP', 'abP']:
            if key in f:
                result[key] = f[key][:]
        
        # Load dimensions
        result['nocc'] = f.attrs['nocc']
        result['nvir'] = f.attrs['nvir']
        result['naux'] = f.attrs['naux']
    
    logger.info(f"DF tensors loaded from {filename}")
    return result


def benchmark_mo_transform(mol: gto.Mole, mo_coeff: np.ndarray, 
                          mo_occ: np.ndarray, auxbasis: str = 'def2-svp-jkfit'):
    """
    Benchmark MO transformation performance
    
    Args:
        mol: PySCF molecule
        mo_coeff: MO coefficients
        mo_occ: MO occupations
        auxbasis: Auxiliary basis set
        
    Returns:
        Dictionary with timing results
    """
    import time
    
    auxmol = df.addons.make_auxmol(mol, auxbasis)
    
    results = {}
    
    # Time (ia|P) construction
    t0 = time.time()
    iaP = build_transition_df_tensor(mol, auxmol, mo_coeff, mo_occ)
    dt_ia = time.time() - t0
    
    results['time_iaP'] = dt_ia
    results['size_iaP_GB'] = iaP.nbytes / 1e9
    results['shape_iaP'] = iaP.shape
    
    # Time full tensor construction
    t0 = time.time()
    df_tensors = build_full_df_tensors(mol, auxmol, mo_coeff, mo_occ)
    dt_full = time.time() - t0
    
    results['time_full'] = dt_full
    results['nocc'] = df_tensors['nocc']
    results['nvir'] = df_tensors['nvir']
    results['naux'] = df_tensors['naux']
    
    # Calculate throughput
    n_elements = (df_tensors['nocc'] * df_tensors['nvir'] + 
                 df_tensors['nocc'] * df_tensors['nocc'] +
                 df_tensors['nvir'] * df_tensors['nvir']) * df_tensors['naux']
    results['throughput_Melem_per_s'] = n_elements / dt_full / 1e6
    
    logger.info("MO transformation benchmark results:")
    logger.info(f"  (ia|P) time: {dt_ia:.2f} s")
    logger.info(f"  Full DF time: {dt_full:.2f} s")
    logger.info(f"  Throughput: {results['throughput_Melem_per_s']:.1f} Melem/s")
    
    return results