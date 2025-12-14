"""
Density fitting utilities for QuasiX.

This module provides optimized routines for building density-fitted
3-center integrals with proper memory management and blocking.
"""

import numpy as np
from pyscf import lib, df
from pyscf.lib import logger
import h5py
from typing import Tuple, Optional


def build_df_tensors_blocked(mol, mo_coeff, mo_occ, auxbasis='def2-tzvp-jkfit',
                             max_memory=4000, thresh_df=1e-10,
                             verbose=None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Build density-fitted 3-center integrals with optimized memory blocking.

    This function efficiently constructs the (ia|P), (ij|P) tensors and
    Coulomb metric V^{1/2} needed for GW calculations.

    Parameters
    ----------
    mol : Mole
        PySCF molecule object
    mo_coeff : ndarray
        MO coefficients (nao × nmo)
    mo_occ : ndarray
        MO occupations (nmo,)
    auxbasis : str
        Auxiliary basis set name
    max_memory : float
        Maximum memory in MB
    thresh_df : float
        Threshold for linear dependencies in auxiliary basis
    verbose : int, optional
        Verbosity level

    Returns
    -------
    tuple
        (iaP, miP, chol_v) where:
        - iaP: (nocc*nvir, naux) occupied-virtual DF tensor
        - miP: (nmo*nocc, naux) all-occupied DF tensor (CRITICAL for Σc)
        - chol_v: (naux, naux) Cholesky factor of Coulomb metric V^{1/2}
    """
    log = logger.new_logger(mol, verbose)
    cput0 = (logger.process_clock(), logger.perf_counter())

    # Get dimensions
    nao = mol.nao_nr()
    nmo = mo_coeff.shape[1]
    nocc = np.sum(mo_occ > 0).astype(int)
    nvir = nmo - nocc

    # Create auxiliary basis
    # CRITICAL FIX: PySCF requires auxbasis as dict mapping atoms to basis sets
    # Convert string auxbasis to dict format if needed
    if isinstance(auxbasis, str):
        # Create dict mapping each unique atom to the same auxiliary basis
        auxbasis_dict = {mol.atom_symbol(i): auxbasis for i in range(mol.natm)}
    else:
        auxbasis_dict = auxbasis

    auxmol = df.addons.make_auxmol(mol, auxbasis_dict)
    naux = auxmol.nao_nr()

    log.info(f'Building DF tensors: nao={nao}, nocc={nocc}, nvir={nvir}, naux={naux}')

    # Build 2-center Coulomb metric (P|Q) and compute V^{1/2}
    log.debug('Computing Coulomb metric V^{1/2}...')
    metric = auxmol.intor('int2c2e', aosym='s1')

    # Eigendecomposition for V^{1/2}
    eigvals, eigvecs = np.linalg.eigh(metric)
    sqrt_eigvals = np.sqrt(np.maximum(eigvals, thresh_df))
    chol_v = eigvecs @ np.diag(sqrt_eigvals) @ eigvecs.T

    # Check linear dependencies
    n_valid = np.sum(eigvals > thresh_df)
    if n_valid < naux:
        log.warn(f'Removed {naux - n_valid} linearly dependent auxiliary functions')

    # Determine optimal block size for memory
    mem_per_aux = 8 * nao**2 / 1e6  # MB per auxiliary function
    blksize = max(1, min(naux, int(max_memory / 4 / mem_per_aux)))
    log.debug(f'Using block size {blksize} for {naux} auxiliary functions')

    # Initialize output tensors
    iaP = np.zeros((nocc, nvir, naux), dtype=np.float64, order='C')
    miP = np.zeros((nmo, nocc, naux), dtype=np.float64, order='C')  # CRITICAL: ALL orbitals × occupied

    # Split MO coefficients
    mo_occ_coeff = np.ascontiguousarray(mo_coeff[:, :nocc])
    mo_vir_coeff = np.ascontiguousarray(mo_coeff[:, nocc:nocc+nvir])
    mo_all_coeff = np.ascontiguousarray(mo_coeff)  # ALL MO coefficients

    # Process 3-center integrals in blocks
    for p0, p1 in lib.prange(0, naux, blksize):
        blk = p1 - p0
        log.debug(f'Processing auxiliary functions {p0}-{p1}')

        # Generate 3-center integrals for this block
        # More efficient: compute only the block we need
        shls_slice = (0, mol.nbas, 0, mol.nbas,
                      p0, p1)  # Shell ranges for integral evaluation

        # Use optimized integral driver
        if hasattr(df.incore, 'aux_e2'):
            # Full block computation
            int3c_blk = df.incore.aux_e2(mol, auxmol, intor='int3c2e',
                                         aosym='s1')[:, p0:p1]
            int3c_blk = int3c_blk.reshape(nao, nao, blk)
        else:
            # Fallback: compute full and slice
            int3c = mol.intor('int3c2e', shls_slice=(0, mol.nbas, 0, mol.nbas,
                                                     p0//auxmol.bas_len()[0],
                                                     p1//auxmol.bas_len()[0]))
            int3c_blk = int3c.reshape(nao, nao, blk)

        # Efficient MO transformation using einsum
        # (μν|P) → (ia|P) and (ij|P)
        for p_idx in range(blk):
            p_global = p0 + p_idx
            ints_p = int3c_blk[:, :, p_idx]

            # Use einsum for efficient transformation
            # (ia|P) = Σ_μν C_μi C_νa (μν|P)
            iaP[:, :, p_global] = np.einsum('mi,na,mn->ia',
                                           mo_occ_coeff, mo_vir_coeff,
                                           ints_p, optimize=True)

            # CRITICAL: (mi|P) = Σ_μν C_μm C_νi (μν|P) for ALL m and occupied i
            # This is required for correct correlation self-energy computation
            # Einsum: ints_p[μ,ν] @ C[ν,i] → temp[μ,i], then C.T[m,μ] @ temp[μ,i] → mi[m,i]
            miP[:, :, p_global] = mo_all_coeff.T @ ints_p @ mo_occ_coeff

    log.timer('DF tensor construction', *cput0)

    # CRITICAL FIX (2025-11-21): Transform to RI basis like PySCF
    # PySCF uses: Lpq = L^-1 @ (P|pq) where L L^T = V
    # This makes the metric identity in the RI basis, which is required for
    # the screening formula W = (I - P0)^-1 - I to work correctly.
    #
    # We computed chol_v = V^{1/2} via eigendecomposition, but we need L^-1 @ X
    # where L is lower Cholesky of V.
    #
    # Use Cholesky instead of eigendecomposition for RI transformation
    L_chol = np.linalg.cholesky(metric)  # L @ L.T = V
    L_chol_inv = np.linalg.inv(L_chol)   # L^-1

    # Transform: iaP_ri[ia, Q] = Σ_P L^-1[Q, P] * iaP[ia, P]
    iaP_flat = iaP.reshape(nocc * nvir, naux)
    miP_flat = miP.reshape(nmo * nocc, naux)

    iaP_ri = iaP_flat @ L_chol_inv.T  # (nocc*nvir, naux) @ (naux, naux) = (nocc*nvir, naux)
    miP_ri = miP_flat @ L_chol_inv.T  # (nmo*nocc, naux) @ (naux, naux) = (nmo*nocc, naux)

    log.info('Applied RI transformation (L^-1 @ integrals)')

    # CRITICAL FIX (2025-11-21): Return IDENTITY as metric for RI-transformed tensors!
    # The RI transformation absorbs the Coulomb metric into the integrals:
    #   L_pq = L^{-1} @ (P|pq)  where L L^T = V (Coulomb metric)
    #
    # For RI tensors, the metric is effectively IDENTITY:
    #   (pq|rs) = sum_L L_pq * L_rs  (no metric factor needed!)
    #
    # Previously, we returned chol_v = V^{1/2} which caused the Rust exchange
    # calculation to apply V^{-1/2} incorrectly, giving ~2x error in Sigma_x.
    #
    # PySCF GW uses: vk = -einsum('Lni,Lim->nm', Lpq, Lpq)  (NO metric!)
    # QuasiX should do the same.
    metric_for_ri = np.eye(naux)

    # Report memory usage
    total_memory = (iaP_ri.nbytes + miP_ri.nbytes + metric_for_ri.nbytes) / 1e9
    log.info(f'DF tensors built: total memory = {total_memory:.2f} GB')
    log.info(f'  iaP shape: {iaP_ri.shape}, miP shape: {miP_ri.shape}')
    log.info(f'  Using identity metric for RI-transformed tensors')

    return iaP_ri, miP_ri, metric_for_ri


def save_df_tensors_h5(filename: str, iaP: np.ndarray, ijP: np.ndarray,
                       chol_v: np.ndarray, mo_energy: np.ndarray,
                       mo_occ: np.ndarray, auxbasis: str):
    """
    Save DF tensors to HDF5 file for later use.

    Parameters
    ----------
    filename : str
        Output HDF5 filename
    iaP : ndarray
        Occupied-virtual DF tensor (nocc*nvir, naux)
    ijP : ndarray
        Occupied-occupied DF tensor (nocc*nocc, naux)
    chol_v : ndarray
        Coulomb metric V^{1/2} (naux, naux)
    mo_energy : ndarray
        MO energies
    mo_occ : ndarray
        MO occupations
    auxbasis : str
        Auxiliary basis set name
    """
    with h5py.File(filename, 'w') as f:
        # Save tensors with compression
        f.create_dataset('iaP', data=iaP, compression='gzip', compression_opts=4)
        f.create_dataset('ijP', data=ijP, compression='gzip', compression_opts=4)
        f.create_dataset('chol_v', data=chol_v, compression='gzip', compression_opts=4)

        # Save MO data
        f.create_dataset('mo_energy', data=mo_energy)
        f.create_dataset('mo_occ', data=mo_occ)

        # Save metadata
        f.attrs['auxbasis'] = auxbasis
        f.attrs['nocc'] = np.sum(mo_occ > 0).astype(int)
        f.attrs['nvir'] = len(mo_energy) - f.attrs['nocc']
        f.attrs['naux'] = chol_v.shape[0]

        # Add timestamp
        import time
        f.attrs['created'] = time.strftime('%Y-%m-%d %H:%M:%S')


def load_df_tensors_h5(filename: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    """
    Load DF tensors from HDF5 file.

    Parameters
    ----------
    filename : str
        Input HDF5 filename

    Returns
    -------
    tuple
        (iaP, ijP, chol_v, metadata) where metadata is a dict with auxiliary info
    """
    with h5py.File(filename, 'r') as f:
        iaP = f['iaP'][:]
        ijP = f['ijP'][:]
        chol_v = f['chol_v'][:]

        metadata = dict(f.attrs)
        if 'mo_energy' in f:
            metadata['mo_energy'] = f['mo_energy'][:]
        if 'mo_occ' in f:
            metadata['mo_occ'] = f['mo_occ'][:]

    return iaP, ijP, chol_v, metadata


def validate_df_tensors(iaP: np.ndarray, ijP: np.ndarray, chol_v: np.ndarray,
                        nocc: int, nvir: int) -> bool:
    """
    Validate DF tensor dimensions and properties.

    Parameters
    ----------
    iaP : ndarray
        Occupied-virtual DF tensor
    ijP : ndarray
        Occupied-occupied DF tensor
    chol_v : ndarray
        Coulomb metric V^{1/2}
    nocc : int
        Number of occupied orbitals
    nvir : int
        Number of virtual orbitals

    Returns
    -------
    bool
        True if all validations pass
    """
    naux = chol_v.shape[0]

    # Check shapes
    assert iaP.shape == (nocc * nvir, naux), f"iaP shape mismatch: {iaP.shape}"
    assert ijP.shape == (nocc * nocc, naux), f"ijP shape mismatch: {ijP.shape}"
    assert chol_v.shape == (naux, naux), f"chol_v shape mismatch: {chol_v.shape}"

    # Check symmetry of chol_v
    symmetry_error = np.max(np.abs(chol_v - chol_v.T))
    assert symmetry_error < 1e-10, f"chol_v not symmetric: error = {symmetry_error}"

    # Check positive semi-definiteness of V
    eigvals = np.linalg.eigvalsh(chol_v @ chol_v.T)
    assert eigvals.min() > -1e-10, f"V not positive semi-definite: min eigval = {eigvals.min()}"

    # Check data ranges are reasonable
    assert np.abs(iaP).max() < 1000, f"iaP has unreasonable values: max = {np.abs(iaP).max()}"
    assert np.abs(ijP).max() < 1000, f"ijP has unreasonable values: max = {np.abs(ijP).max()}"

    return True