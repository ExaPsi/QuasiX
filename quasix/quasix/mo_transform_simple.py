"""
Simplified MO transformation for PySCF adapter - production quality
"""

import numpy as np
from pyscf import lib, df
from quasix.logging import get_logger

logger = get_logger(__name__)

def validate_mo_coefficients(mo_coeff, nao, mo_occ, overlap_matrix=None):
    """Validate MO coefficients with proper orthonormality check"""
    if mo_coeff.shape[0] != nao:
        raise ValueError(f"MO shape {mo_coeff.shape} inconsistent with nao={nao}")
    
    # Check C^T @ S @ C = I
    if overlap_matrix is not None:
        ortho = mo_coeff.T @ overlap_matrix @ mo_coeff
    else:
        ortho = mo_coeff.T @ mo_coeff
    
    ortho_error = np.max(np.abs(ortho - np.eye(ortho.shape[0])))
    if ortho_error > 1e-6:
        logger.warning(f"MO not perfectly orthonormal (error: {ortho_error:.2e})")
    if ortho_error > 1e-3:
        raise ValueError(f"MO poorly orthonormal (error: {ortho_error:.2e})")
    
    return True

def build_transition_df_tensor(mol, auxbasis, mo_coeff, mo_occ):
    """Build (ia|P) tensor using PySCF DF module"""
    # Get dimensions
    nocc = int(np.sum(mo_occ > 0))
    nvir = mol.nao - nocc
    
    # Split MO coefficients
    c_occ = mo_coeff[:, :nocc]
    c_vir = mo_coeff[:, nocc:]
    
    # Build DF object
    mydf = df.DF(mol)
    if auxbasis:
        mydf.auxbasis = auxbasis
    mydf.build()
    
    naux = mydf.get_naoaux()
    logger.info(f"Building (ia|P): nocc={nocc}, nvir={nvir}, naux={naux}")
    
    # Get 3-center integrals transformed to MO basis
    # For DF, we need to use ao2mo with proper syntax
    # The function expects (mo_i, mo_j, mo_k, mo_l) for 4-index, but for DF we use special syntax
    from pyscf.ao2mo import _ao2mo
    
    # Get 3-center integrals (μν|P) and transform to (ia|P)
    # This is more direct: get the DF integrals and transform
    iaP = np.zeros((nocc * nvir, naux))
    
    # Use the DF module's internal method to get transformed integrals
    # This gets (ij|P) where i,j are MO indices
    mo_coeffs = (c_occ, c_vir)
    eri_mo = mydf.ao2mo(mo_coeffs, compact=False)
    
    # Reshape to (ia|P) format
    iaP = eri_mo.reshape(nocc, nvir, naux).transpose(0,1,2).reshape(nocc*nvir, naux)
    
    return iaP

def build_metric_tensor(mol, auxbasis):
    """Build (P|Q) metric tensor"""
    mydf = df.DF(mol)
    if auxbasis:
        mydf.auxbasis = auxbasis
    mydf.build()
    
    # Get 2-center integrals
    metric = mydf.get_2c2e()
    
    return metric