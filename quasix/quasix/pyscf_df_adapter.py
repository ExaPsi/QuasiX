"""
PySCF adapter for density fitting tensor operations.

This module provides seamless integration between PySCF's DF module
and QuasiX's optimized tensor transformations.
"""

import numpy as np
from typing import Optional, Tuple, Union, Dict, Any
import logging

logger = logging.getLogger(__name__)

# Import QuasiX components
from .df_tensors import DFTensorTransformer

# Try to import PySCF
try:
    from pyscf import gto, scf, df, lib
    from pyscf.df import incore
    PYSCF_AVAILABLE = True
except ImportError:
    PYSCF_AVAILABLE = False
    logger.warning("PySCF not available. PySCF adapter functions will be limited.")


class PySCFDFAdapter:
    """
    Adapter for using QuasiX DF tensor operations with PySCF.

    This class provides methods to extract DF tensors from PySCF calculations
    and transform them using QuasiX's optimized routines.
    """

    def __init__(self, mol: 'gto.Mole', auxbasis: str = 'def2-tzvp-jkfit'):
        """
        Initialize PySCF DF adapter.

        Args:
            mol: PySCF molecule object
            auxbasis: Auxiliary basis set name
        """
        if not PYSCF_AVAILABLE:
            raise RuntimeError("PySCF is required for this adapter")

        self.mol = mol
        self.auxbasis = auxbasis
        self.df_obj = df.DF(mol, auxbasis)

        # Get dimensions
        self.n_ao = mol.nao
        self.n_aux = self.df_obj.get_naoaux()

        logger.info(f"Initialized PySCF DF adapter: n_ao={self.n_ao}, n_aux={self.n_aux}, auxbasis={auxbasis}")

    def get_3center_ao(self) -> np.ndarray:
        """
        Get 3-center integrals in AO basis from PySCF.

        Returns:
            3-center integrals (μν|P) of shape (n_ao, n_ao, n_aux)
        """
        # Build the DF integrals if not already built
        if not hasattr(self.df_obj, '_cderi'):
            self.df_obj.build()

        # Get integrals from PySCF - they're stored in _cderi
        # Shape is typically (n_aux, n_ao*(n_ao+1)//2) for compact storage
        # We need to unpack to full (n_ao, n_ao, n_aux)
        cderi = self.df_obj._cderi

        if cderi is None:
            # Build on the fly
            self.df_obj.build()
            cderi = self.df_obj._cderi

        n_ao = self.n_ao
        n_aux = self.n_aux

        # Unpack the compact storage to full matrix
        j3c = np.zeros((n_ao, n_ao, n_aux))
        for p in range(n_aux):
            # Get the p-th auxiliary function's integrals
            # This is stored as a lower triangular matrix in compact form
            j3c_p = lib.unpack_tril(cderi[p])
            j3c[:, :, p] = j3c_p

        logger.debug(f"Retrieved 3-center integrals: shape={j3c.shape}, min={j3c.min():.6f}, max={j3c.max():.6f}")

        return j3c

    def get_2center_metric(self) -> np.ndarray:
        """
        Get 2-center metric matrix from PySCF.

        Returns:
            Metric matrix (P|Q) of shape (n_aux, n_aux)
        """
        # Build if needed
        if not hasattr(self.df_obj, '_cderi'):
            self.df_obj.build()

        # Get the 2-center integrals (P|Q)
        # This is the metric matrix for the auxiliary basis
        metric = self.df_obj.auxmol.intor('int2c2e')

        logger.debug(f"Retrieved 2-center metric: shape={metric.shape}, condition={np.linalg.cond(metric):.2e}")

        return metric

    def transform_to_mo(
        self,
        mf: 'scf.hf.SCF',
        mo_coeff: Optional[np.ndarray] = None
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Transform 3-center integrals to MO basis using QuasiX.

        Args:
            mf: PySCF SCF object (RHF/UHF)
            mo_coeff: Optional MO coefficients (uses mf.mo_coeff if None)

        Returns:
            Tuple of:
                - j3c_ia: (ia|P) tensor of shape (n_occ*n_vir, n_aux)
                - j3c_ij: (ij|P) tensor of shape (n_occ*n_occ, n_aux)
                - j3c_ab: (ab|P) tensor of shape (n_vir*n_vir, n_aux)
        """
        # Get MO coefficients
        if mo_coeff is None:
            mo_coeff = mf.mo_coeff

        # Get occupation
        mo_occ = mf.mo_occ
        n_occ = int(np.sum(mo_occ > 0))
        n_vir = self.n_ao - n_occ

        # Split coefficients
        c_occ = mo_coeff[:, :n_occ]
        c_vir = mo_coeff[:, n_occ:]

        logger.info(f"MO transformation: n_occ={n_occ}, n_vir={n_vir}")

        # Get 3-center integrals
        j3c_ao = self.get_3center_ao()

        # Create transformer
        transformer = DFTensorTransformer(self.n_ao, n_occ, n_vir, self.n_aux)

        # Transform (ia|P)
        j3c_ia = transformer.transform_mo_3center(j3c_ao, c_occ, c_vir)

        # Transform (ij|P)
        transformer_ij = DFTensorTransformer(self.n_ao, n_occ, n_occ, self.n_aux)
        j3c_ij = transformer_ij.transform_mo_3center(j3c_ao, c_occ, c_occ)

        # Transform (ab|P)
        transformer_ab = DFTensorTransformer(self.n_ao, n_vir, n_vir, self.n_aux)
        j3c_ab = transformer_ab.transform_mo_3center(j3c_ao, c_vir, c_vir)

        return j3c_ia, j3c_ij, j3c_ab

    def compute_cholesky_metric(self) -> Dict[str, Any]:
        """
        Compute Cholesky factorization of the metric matrix.

        Returns:
            Dictionary with Cholesky factor and condition number
        """
        metric = self.get_2center_metric()

        # Create transformer (dimensions don't matter for Cholesky)
        transformer = DFTensorTransformer(self.n_ao, 1, 1, self.n_aux)

        return transformer.compute_cholesky_v(metric)

    def build_symmetrized_df(
        self,
        mf: 'scf.hf.SCF'
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Build symmetrized DF tensors for GW calculations.

        Computes B' = V^(1/2) @ B where B are the MO-transformed tensors.

        Args:
            mf: PySCF SCF object

        Returns:
            Tuple of:
                - b_ia_sym: Symmetrized (ia|P) tensor
                - vsqrt: V^(1/2) matrix
        """
        # Get Cholesky factor of metric
        cholesky_result = self.compute_cholesky_metric()
        L = cholesky_result['L']

        # For symmetric DF, we want V^(1/2) not L
        # Use eigendecomposition for V^(1/2)
        metric = self.get_2center_metric()
        eigvals, eigvecs = np.linalg.eigh(metric)

        # Check for negative eigenvalues
        if eigvals.min() < 0:
            logger.warning(f"Negative eigenvalues in metric: min={eigvals.min():.2e}")
            eigvals = np.maximum(eigvals, 1e-14)

        # Compute V^(1/2)
        vsqrt = eigvecs @ np.diag(np.sqrt(eigvals)) @ eigvecs.T

        # Get MO-transformed tensors
        j3c_ia, _, _ = self.transform_to_mo(mf)

        # Symmetrize: B' = B @ V^(1/2)
        b_ia_sym = j3c_ia @ vsqrt

        logger.info(f"Built symmetrized DF tensors: shape={b_ia_sym.shape}")

        return b_ia_sym, vsqrt


def create_df_from_pyscf(
    mol: 'gto.Mole',
    auxbasis: str = 'def2-tzvp-jkfit'
) -> PySCFDFAdapter:
    """
    Convenience function to create DF adapter from PySCF molecule.

    Args:
        mol: PySCF molecule object
        auxbasis: Auxiliary basis set name

    Returns:
        PySCFDFAdapter instance
    """
    if not PYSCF_AVAILABLE:
        raise RuntimeError("PySCF is required for this function")

    return PySCFDFAdapter(mol, auxbasis)


def validate_df_transformation(
    mol: 'gto.Mole',
    mf: 'scf.hf.SCF',
    auxbasis: str = 'def2-tzvp-jkfit',
    tolerance: float = 1e-8
) -> bool:
    """
    Validate QuasiX DF transformation against PySCF.

    Args:
        mol: PySCF molecule
        mf: Converged SCF object
        auxbasis: Auxiliary basis set
        tolerance: Maximum allowed difference

    Returns:
        True if validation passes
    """
    if not PYSCF_AVAILABLE:
        logger.warning("PySCF not available for validation")
        return True

    # Create adapter
    adapter = PySCFDFAdapter(mol, auxbasis)

    # Get dimensions
    mo_occ = mf.mo_occ
    n_occ = int(np.sum(mo_occ > 0))
    n_vir = mol.nao - n_occ

    # Get MO coefficients
    c_occ = mf.mo_coeff[:, :n_occ]
    c_vir = mf.mo_coeff[:, n_occ:]

    # Get 3-center integrals
    j3c_ao = adapter.get_3center_ao()

    # Transform with PySCF
    j3c_ia_pyscf = lib.einsum('mnP,mi,na->iaP', j3c_ao, c_occ, c_vir)
    j3c_ia_pyscf = j3c_ia_pyscf.reshape(n_occ * n_vir, adapter.n_aux)

    # Transform with QuasiX
    transformer = DFTensorTransformer(mol.nao, n_occ, n_vir, adapter.n_aux)
    j3c_ia_quasix = transformer.transform_mo_3center(j3c_ao, c_occ, c_vir)

    # Compare
    max_diff = np.abs(j3c_ia_quasix - j3c_ia_pyscf).max()
    passed = max_diff < tolerance

    if passed:
        logger.info(f"✓ DF transformation validation passed (max diff: {max_diff:.2e})")
    else:
        logger.error(f"✗ DF transformation validation failed (max diff: {max_diff:.2e} > {tolerance:.2e})")

    return passed


def benchmark_pyscf_vs_quasix(
    mol: 'gto.Mole',
    mf: 'scf.hf.SCF',
    auxbasis: str = 'def2-tzvp-jkfit'
) -> Dict[str, float]:
    """
    Benchmark QuasiX against PySCF for DF transformations.

    Args:
        mol: PySCF molecule
        mf: Converged SCF object
        auxbasis: Auxiliary basis set

    Returns:
        Dictionary with timing results
    """
    import time

    if not PYSCF_AVAILABLE:
        logger.warning("PySCF not available for benchmarking")
        return {}

    # Create adapter
    adapter = PySCFDFAdapter(mol, auxbasis)

    # Get dimensions
    mo_occ = mf.mo_occ
    n_occ = int(np.sum(mo_occ > 0))
    n_vir = mol.nao - n_occ

    # Get data
    j3c_ao = adapter.get_3center_ao()
    c_occ = mf.mo_coeff[:, :n_occ]
    c_vir = mf.mo_coeff[:, n_occ:]

    results = {}

    # Benchmark PySCF
    t0 = time.time()
    j3c_ia_pyscf = lib.einsum('mnP,mi,na->iaP', j3c_ao, c_occ, c_vir)
    j3c_ia_pyscf = j3c_ia_pyscf.reshape(n_occ * n_vir, adapter.n_aux)
    results['pyscf_time'] = time.time() - t0

    # Benchmark QuasiX
    transformer = DFTensorTransformer(mol.nao, n_occ, n_vir, adapter.n_aux)
    t0 = time.time()
    j3c_ia_quasix = transformer.transform_mo_3center(j3c_ao, c_occ, c_vir)
    results['quasix_time'] = time.time() - t0

    # Compare
    results['max_diff'] = np.abs(j3c_ia_quasix - j3c_ia_pyscf).max()
    results['speedup'] = results['pyscf_time'] / results['quasix_time']

    logger.info(f"Benchmark results:")
    logger.info(f"  PySCF time:  {results['pyscf_time']:.3f} s")
    logger.info(f"  QuasiX time: {results['quasix_time']:.3f} s")
    logger.info(f"  Speedup:     {results['speedup']:.2f}x")
    logger.info(f"  Max diff:    {results['max_diff']:.2e}")

    return results


# Example usage
if __name__ == "__main__" and PYSCF_AVAILABLE:
    logging.basicConfig(level=logging.INFO)

    # Create a test molecule
    mol = gto.M(
        atom='O 0 0 0; H 0 0 1; H 0 1 0',
        basis='cc-pvdz',
        verbose=0
    )

    # Run HF
    mf = scf.RHF(mol)
    mf.kernel()

    # Create DF adapter
    adapter = create_df_from_pyscf(mol, 'cc-pvdz-jkfit')

    # Transform to MO basis
    j3c_ia, j3c_ij, j3c_ab = adapter.transform_to_mo(mf)
    print(f"Transformed tensors: ia={j3c_ia.shape}, ij={j3c_ij.shape}, ab={j3c_ab.shape}")

    # Compute Cholesky
    cholesky = adapter.compute_cholesky_metric()
    print(f"Cholesky condition number: {cholesky['condition_number']:.2e}")

    # Validate
    is_valid = validate_df_transformation(mol, mf, 'cc-pvdz-jkfit')
    print(f"Validation: {'PASSED' if is_valid else 'FAILED'}")

    # Benchmark
    results = benchmark_pyscf_vs_quasix(mol, mf, 'cc-pvdz-jkfit')
    print(f"Performance: {results.get('speedup', 1.0):.2f}x speedup")