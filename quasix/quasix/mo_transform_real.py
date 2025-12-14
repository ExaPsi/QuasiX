"""
Real MO transformation utilities for PySCF adapter
Following QuasiX golden rule: real libraries only, results comparable to PySCF
"""

import numpy as np
from pyscf import lib, df, ao2mo
from pyscf.df import incore, addons
from quasix.logging import get_logger
from contextlib import contextmanager
import time

# Set up module logger
logger = get_logger(__name__)

@contextmanager
def timed_stage(name):
    """Context manager for timing stages"""
    start = time.time()
    logger.debug(f"Starting {name}")
    try:
        yield
    finally:
        elapsed = time.time() - start
        logger.debug(f"Completed {name} in {elapsed:.2f}s")

def validate_mo_coefficients(mo_coeff: np.ndarray, nao: int, 
                            mo_occ: np.ndarray, overlap_matrix: np.ndarray = None) -> None:
    """
    Validate MO coefficient matrix with proper orthonormality check
    
    Args:
        mo_coeff: MO coefficient matrix [nao, nmo]
        nao: Number of atomic orbitals
        mo_occ: MO occupations [nmo]
        overlap_matrix: AO overlap matrix S [nao, nao]
        
    Raises:
        ValueError: If MO coefficients are invalid
    """
    if mo_coeff.ndim != 2:
        raise ValueError(f"MO coefficients must be 2D, got shape {mo_coeff.shape}")
    
    if mo_coeff.shape[0] != nao:
        raise ValueError(f"MO coefficient shape {mo_coeff.shape} inconsistent with nao={nao}")
    
    if mo_coeff.shape[1] != len(mo_occ):
        raise ValueError(f"MO coefficient shape {mo_coeff.shape} inconsistent with mo_occ length {len(mo_occ)}")
    
    # Check for NaN or Inf values
    if np.any(np.isnan(mo_coeff)):
        raise ValueError("MO coefficients contain NaN values")
    if np.any(np.isinf(mo_coeff)):
        raise ValueError("MO coefficients contain Inf values")
    
    # Check for orthonormality: C^T @ S @ C = I
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

def transform_mo_coefficients(mo_coeff: np.ndarray, mo_occ: np.ndarray):
    """
    Split MO coefficients into occupied and virtual
    
    Args:
        mo_coeff: Full MO coefficient matrix [nao, nmo]
        mo_occ: MO occupations [nmo]
        
    Returns:
        c_occ: Occupied MO coefficients [nao, nocc]
        c_vir: Virtual MO coefficients [nao, nvir]
    """
    nocc = int(np.sum(mo_occ > 0))
    c_occ = mo_coeff[:, :nocc]
    c_vir = mo_coeff[:, nocc:]
    return c_occ, c_vir

def build_transition_df_tensor(mol, auxbasis, mo_coeff: np.ndarray,
                              mo_occ: np.ndarray, max_memory: int = 2000) -> np.ndarray:
    """
    Build density-fitted transition integrals (ia|P) using real PySCF DF
    
    This properly constructs the 3-center integrals in the occupied-virtual
    MO basis, which are needed for GW calculations.
    
    Args:
        mol: PySCF molecule object
        auxbasis: Auxiliary basis name (e.g., 'def2-svp-jkfit')
        mo_coeff: MO coefficients [nao, nmo]
        mo_occ: MO occupations [nmo]
        max_memory: Maximum memory in MB
        
    Returns:
        iaP: Transition DF tensor of shape [nocc*nvir, naux]
    """
    with timed_stage("build_transition_df"):
        nao = mol.nao
        nmo = mo_coeff.shape[1]
        nocc = int(np.sum(mo_occ > 0))
        nvir = nmo - nocc
        
        logger.info(f"Building (ia|P) tensor: nocc={nocc}, nvir={nvir}, nao={nao}")
        
        # Split MO coefficients
        c_occ, c_vir = transform_mo_coefficients(mo_coeff, mo_occ)
        
        # Build auxiliary molecule
        if isinstance(auxbasis, str):
            auxmol = addons.make_auxmol(mol, auxbasis)
        else:
            auxmol = auxbasis  # Assume it's already an auxmol
        
        naux = auxmol.nao
        logger.info(f"Auxiliary basis size: {naux}")
        
        # Get 3-center integrals (μν|P) using incore module
        with timed_stage("3center_integrals"):
            # Get 3-center integrals in AO basis
            j3c = incore.aux_e2(mol, auxmol, intor='int3c2e', aosym='s1')
            
            # j3c shape is (nao, nao, naux) for s1 symmetry
            # or could be compressed form - let's ensure proper shape
            if j3c.ndim == 2:
                # Compressed form (nao*(nao+1)/2, naux)
                # Need to unpack for transformation
                nao = mol.nao
                j3c_full = np.zeros((nao, nao, naux))
                idx = 0
                for i in range(nao):
                    for j in range(i+1):
                        j3c_full[i,j] = j3c[idx]
                        j3c_full[j,i] = j3c[idx]
                        idx += 1
                j3c = j3c_full
            elif j3c.shape[0] * j3c.shape[1] == mol.nao * mol.nao:
                # Reshape from (nao*nao, naux) to (nao, nao, naux)
                j3c = j3c.reshape(mol.nao, mol.nao, naux)
        
        # Transform to MO basis: (μν|P) -> (ia|P)
        with timed_stage("mo_transform"):
            # (ia|P) = Σ_μν C_μi C_νa (μν|P)
            iaP = np.einsum('pi,qa,pqr->iar', c_occ, c_vir, j3c)
            
            # Reshape to (nocc*nvir, naux) for GW calculations
            iaP = iaP.reshape(nocc * nvir, naux)
        
        logger.info(f"Built (ia|P) tensor with shape {iaP.shape}")
        
        # Validate the result
        if np.any(np.isnan(iaP)):
            raise ValueError("NaN values in (ia|P) tensor")
        if np.any(np.isinf(iaP)):
            raise ValueError("Inf values in (ia|P) tensor")
            
        return iaP

def build_full_df_tensors(mol, auxbasis, mo_coeff: np.ndarray, 
                         mo_occ: np.ndarray, max_memory: int = 2000):
    """
    Build all DF tensors needed for GW calculation
    
    Returns:
        dict with keys:
            'iaP': (ia|P) transition integrals [nocc*nvir, naux]
            'ijP': (ij|P) occupied-occupied [nocc*nocc, naux]
            'abP': (ab|P) virtual-virtual [nvir*nvir, naux]
            'v_PQ': (P|Q) metric matrix [naux, naux]
    """
    with timed_stage("build_full_df_tensors"):
        # Check for valid inputs
        if mo_coeff is None or mo_occ is None:
            raise ValueError("MO coefficients and occupations must be provided. Run SCF first.")
        
        nocc = int(np.sum(mo_occ > 0))
        nvir = mo_coeff.shape[1] - nocc
        
        # Split MO coefficients
        c_occ, c_vir = transform_mo_coefficients(mo_coeff, mo_occ)
        
        # Build auxiliary molecule
        if isinstance(auxbasis, str):
            auxmol = addons.make_auxmol(mol, auxbasis)
        else:
            auxmol = auxbasis
        
        naux = auxmol.nao
        logger.info(f"Building full DF tensors: nocc={nocc}, nvir={nvir}, naux={naux}")
        
        # Get 3-center integrals (μν|P) using incore module
        with timed_stage("3center_integrals"):
            j3c = incore.aux_e2(mol, auxmol, intor='int3c2e', aosym='s1')
            
            # Ensure proper shape
            if j3c.ndim == 2:
                nao = mol.nao
                j3c_full = np.zeros((nao, nao, naux))
                idx = 0
                for i in range(nao):
                    for j in range(i+1):
                        j3c_full[i,j] = j3c[idx]
                        j3c_full[j,i] = j3c[idx]
                        idx += 1
                j3c = j3c_full
            elif j3c.shape[0] * j3c.shape[1] == mol.nao * mol.nao:
                j3c = j3c.reshape(mol.nao, mol.nao, naux)
        
        tensors = {}
        
        # Build (ia|P)
        with timed_stage("iaP"):
            tensors['iaP'] = np.einsum('pi,qa,pqr->iar', c_occ, c_vir, j3c).reshape(nocc * nvir, naux)
        
        # Build (ij|P)
        with timed_stage("ijP"):
            tensors['ijP'] = np.einsum('pi,qj,pqr->ijr', c_occ, c_occ, j3c).reshape(nocc * nocc, naux)
        
        # Build (ab|P)
        with timed_stage("abP"):
            tensors['abP'] = np.einsum('pa,qb,pqr->abr', c_vir, c_vir, j3c).reshape(nvir * nvir, naux)
        
        # Build (P|Q) metric using 2-center integrals
        with timed_stage("v_PQ"):
            tensors['v_PQ'] = incore.fill_2c2e(mol, auxmol)
        
        return tensors

def validate_df_tensors(tensors: dict, nocc: int, nvir: int, naux: int,
                        tolerance: float = 1e-8):
    """
    Validate DF tensors for correctness
    
    Args:
        tensors: Dictionary of DF tensors
        nocc: Number of occupied orbitals
        nvir: Number of virtual orbitals
        naux: Number of auxiliary functions
        tolerance: Numerical tolerance for validation
    """
    # Check shapes
    expected_shapes = {
        'iaP': (nocc * nvir, naux),
        'ijP': (nocc * nocc, naux),
        'abP': (nvir * nvir, naux),
        'v_PQ': (naux, naux)
    }
    
    for key, expected_shape in expected_shapes.items():
        if key in tensors:
            actual_shape = tensors[key].shape
            if actual_shape != expected_shape:
                raise ValueError(f"{key} has wrong shape: {actual_shape} != {expected_shape}")
    
    # Check for NaN/Inf
    for key, tensor in tensors.items():
        if np.any(np.isnan(tensor)):
            raise ValueError(f"NaN values in {key}")
        if np.any(np.isinf(tensor)):
            raise ValueError(f"Inf values in {key}")
    
    # Check metric positive definiteness
    if 'v_PQ' in tensors:
        eigvals = np.linalg.eigvalsh(tensors['v_PQ'])
        if np.min(eigvals) < -tolerance:
            raise ValueError(f"Metric not positive definite: min eigenvalue = {np.min(eigvals)}")
    
    logger.info("DF tensor validation passed")

# For backward compatibility
def build_metric_tensor(mol, auxbasis):
    """Build (P|Q) metric tensor"""
    if isinstance(auxbasis, str):
        auxmol = addons.make_auxmol(mol, auxbasis)
    else:
        auxmol = auxbasis
    return incore.fill_2c2e(mol, auxmol)