"""
Optimized PySCF adapter for density fitting tensor operations.

This module provides high-performance integration between PySCF's DF module
and QuasiX's optimized tensor transformations.

Key optimizations:
1. Lazy loading of DF integrals
2. Memory-mapped file access for large systems
3. Batch processing of multiple frequencies
4. Efficient unpacking of compact storage
5. Parallel transformation support
"""

import numpy as np
from typing import Optional, Tuple, Union, Dict, Any, List, Callable
import logging
from functools import lru_cache, wraps
from concurrent.futures import ThreadPoolExecutor, as_completed
import h5py
import tempfile
import os

logger = logging.getLogger(__name__)

# Import QuasiX optimized components
from .df_tensors_optimized import OptimizedDFTensorTransformer, ArrayCache

# Try to import PySCF
try:
    from pyscf import gto, scf, df, lib
    from pyscf.df import incore, outcore
    PYSCF_AVAILABLE = True
except ImportError:
    PYSCF_AVAILABLE = False
    logger.warning("PySCF not available. PySCF adapter functions will be limited.")


class OptimizedPySCFDFAdapter:
    """
    High-performance adapter for using QuasiX DF tensor operations with PySCF.

    Key features:
    - Lazy loading of integrals
    - Memory-mapped file support for large systems
    - Efficient unpacking of compact storage
    - Batch transformation support
    - Parallel processing capabilities
    """

    def __init__(
        self,
        mol: 'gto.Mole',
        auxbasis: str = 'def2-tzvp-jkfit',
        use_mmap: bool = False,
        max_memory: int = 4000,  # MB
        parallel: bool = True
    ):
        """
        Initialize optimized PySCF DF adapter.

        Args:
            mol: PySCF molecule object
            auxbasis: Auxiliary basis set name
            use_mmap: Use memory-mapped files for large systems
            max_memory: Maximum memory in MB
            parallel: Enable parallel processing
        """
        if not PYSCF_AVAILABLE:
            raise RuntimeError("PySCF is required for this adapter")

        self.mol = mol
        self.auxbasis = auxbasis
        self.use_mmap = use_mmap
        self.max_memory = max_memory
        self.parallel = parallel and (mol.nao * mol.nao * self._estimate_naux() * 8 / 1e6 > 1000)

        # Create DF object with optimized settings
        self.df_obj = df.DF(mol, auxbasis)
        self.df_obj.max_memory = max_memory

        # Get dimensions
        self.n_ao = mol.nao
        self.n_aux = None  # Lazy evaluation
        self._cderi = None  # Lazy loading
        self._cderi_file = None  # For memory-mapped access

        # Cache for transformed tensors
        self.cache = ArrayCache(max_size=20)
        self._metric_cache = None
        self._vsqrt_cache = None

        # Thread pool for parallel operations
        self.executor = ThreadPoolExecutor(max_workers=4) if parallel else None

        logger.info(f"Initialized OptimizedPySCFDFAdapter: n_ao={self.n_ao}, "
                   f"auxbasis={auxbasis}, use_mmap={use_mmap}, parallel={parallel}")

    def _estimate_naux(self) -> int:
        """Estimate number of auxiliary functions."""
        # Rough estimate based on basis set
        if 'def2' in self.auxbasis.lower():
            return self.n_ao * 3
        elif 'cc-pv' in self.auxbasis.lower():
            return self.n_ao * 4
        else:
            return self.n_ao * 3

    @property
    def naux(self) -> int:
        """Get number of auxiliary basis functions (lazy evaluation)."""
        if self.n_aux is None:
            self.n_aux = self.df_obj.get_naoaux()
        return self.n_aux

    def _build_df_lazy(self):
        """Build DF integrals lazily when first needed."""
        if self._cderi is not None:
            return

        logger.debug("Building DF integrals (lazy loading)...")

        if self.use_mmap:
            # Use memory-mapped file for large systems
            with tempfile.NamedTemporaryFile(delete=False, suffix='.h5') as tmpfile:
                self._cderi_file = tmpfile.name

            self.df_obj._cderi_to_save = self._cderi_file
            self.df_obj.build()

            # Open for memory-mapped access
            self._cderi = h5py.File(self._cderi_file, 'r')['eri_mo']
        else:
            # In-memory storage
            self.df_obj.build()
            self._cderi = self.df_obj._cderi

        logger.debug(f"DF integrals built: shape={self._cderi.shape if hasattr(self._cderi, 'shape') else 'file'}")

    @lru_cache(maxsize=1)
    def get_3center_ao(self) -> np.ndarray:
        """
        Get 3-center integrals in AO basis with optimized unpacking.

        Returns:
            3-center integrals (μν|P) of shape (n_ao, n_ao, n_aux)
        """
        self._build_df_lazy()

        n_ao = self.n_ao
        n_aux = self.naux

        # Check memory requirements
        mem_required = n_ao * n_ao * n_aux * 8 / 1e6  # MB
        if mem_required > self.max_memory:
            logger.warning(f"3-center integrals require {mem_required:.1f} MB, "
                          f"exceeding limit {self.max_memory} MB")

        # Optimized unpacking based on storage format
        if isinstance(self._cderi, h5py.Dataset):
            # Memory-mapped file access
            j3c = np.zeros((n_ao, n_ao, n_aux))
            chunk_size = min(100, n_aux)  # Process in chunks

            for p0 in range(0, n_aux, chunk_size):
                p1 = min(p0 + chunk_size, n_aux)
                chunk = self._cderi[p0:p1]
                for i, p in enumerate(range(p0, p1)):
                    j3c[:, :, p] = lib.unpack_tril(chunk[i])

        elif self._cderi.ndim == 2:
            # Compact storage format
            j3c = np.zeros((n_ao, n_ao, n_aux))

            if self.parallel and self.executor:
                # Parallel unpacking
                futures = []
                for p in range(n_aux):
                    future = self.executor.submit(lib.unpack_tril, self._cderi[p])
                    futures.append((p, future))

                for p, future in futures:
                    j3c[:, :, p] = future.result()
            else:
                # Serial unpacking with optimized loop
                for p in range(n_aux):
                    j3c[:, :, p] = lib.unpack_tril(self._cderi[p])

        else:
            # Already in full format
            j3c = self._cderi

        logger.debug(f"Retrieved 3-center integrals: shape={j3c.shape}, "
                    f"memory={j3c.nbytes/1e6:.1f} MB")

        return j3c

    @lru_cache(maxsize=1)
    def get_2center_metric(self) -> np.ndarray:
        """
        Get 2-center metric matrix with caching.

        Returns:
            Metric matrix (P|Q) of shape (n_aux, n_aux)
        """
        if self._metric_cache is not None:
            return self._metric_cache

        # Build if needed
        self._build_df_lazy()

        # Compute 2-center integrals
        metric = self.df_obj.auxmol.intor('int2c2e')

        # Cache result
        self._metric_cache = metric

        logger.debug(f"Retrieved 2-center metric: shape={metric.shape}, "
                    f"condition={np.linalg.cond(metric):.2e}")

        return metric

    @lru_cache(maxsize=1)
    def get_vsqrt(self) -> np.ndarray:
        """
        Get V^(1/2) matrix with caching.

        Returns:
            Square root of metric matrix
        """
        if self._vsqrt_cache is not None:
            return self._vsqrt_cache

        metric = self.get_2center_metric()

        # Eigendecomposition for V^(1/2)
        eigvals, eigvecs = np.linalg.eigh(metric)

        # Check for negative eigenvalues
        min_eigval = eigvals.min()
        if min_eigval < 0:
            logger.warning(f"Negative eigenvalues in metric: min={min_eigval:.2e}")
            eigvals = np.maximum(eigvals, 1e-14)

        # Compute V^(1/2)
        vsqrt = eigvecs @ np.diag(np.sqrt(eigvals)) @ eigvecs.T

        # Cache result
        self._vsqrt_cache = vsqrt

        return vsqrt

    def transform_to_mo_batch(
        self,
        mf: 'scf.hf.SCF',
        mo_coeff_list: Optional[List[np.ndarray]] = None
    ) -> List[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """
        Transform 3-center integrals to MO basis for multiple MO coefficients.

        Args:
            mf: PySCF SCF object
            mo_coeff_list: List of MO coefficient matrices

        Returns:
            List of (j3c_ia, j3c_ij, j3c_ab) tuples
        """
        if mo_coeff_list is None:
            mo_coeff_list = [mf.mo_coeff]

        # Get 3-center integrals once
        j3c_ao = self.get_3center_ao()

        results = []

        for mo_coeff in mo_coeff_list:
            # Get occupation
            mo_occ = mf.mo_occ if len(mo_coeff_list) == 1 else np.ones(mo_coeff.shape[1])
            n_occ = int(np.sum(mo_occ > 0))
            n_vir = self.n_ao - n_occ

            # Split coefficients
            c_occ = mo_coeff[:, :n_occ]
            c_vir = mo_coeff[:, n_occ:]

            # Create optimized transformer
            transformer = OptimizedDFTensorTransformer(
                self.n_ao, n_occ, n_vir, self.naux, enable_cache=True
            )

            # Transform with optimizations
            j3c_ia = transformer.transform_mo_3center(j3c_ao, c_occ, c_vir, use_blas=True)

            # Also compute ij and ab if needed
            j3c_ij = transformer.transform_mo_3center(j3c_ao, c_occ, c_occ, use_blas=True)
            j3c_ab = transformer.transform_mo_3center(j3c_ao, c_vir, c_vir, use_blas=True)

            results.append((j3c_ia, j3c_ij, j3c_ab))

        return results

    def compute_polarizability_batch(
        self,
        frequencies: np.ndarray,
        mf: 'scf.hf.SCF',
        eta: float = 0.01
    ) -> np.ndarray:
        """
        Compute polarizability P0(ω) for multiple frequencies in parallel.

        Args:
            frequencies: Array of frequencies
            mf: Converged SCF object
            eta: Broadening parameter

        Returns:
            Array of shape (n_freq, n_aux, n_aux) with P0 matrices
        """
        # Get MO-transformed tensors
        mo_occ = mf.mo_occ
        n_occ = int(np.sum(mo_occ > 0))
        n_vir = self.n_ao - n_occ

        # Get energies
        e_occ = mf.mo_energy[:n_occ]
        e_virt = mf.mo_energy[n_occ:]

        # Get DF tensor (ia|P)
        j3c_ao = self.get_3center_ao()
        c_occ = mf.mo_coeff[:, :n_occ]
        c_vir = mf.mo_coeff[:, n_occ:]

        transformer = OptimizedDFTensorTransformer(
            self.n_ao, n_occ, n_vir, self.naux, enable_cache=True
        )
        df_ia = transformer.transform_mo_3center(j3c_ao, c_occ, c_vir, use_blas=True)
        df_ia = df_ia.reshape(n_occ, n_vir, self.naux)

        # Compute P0 for each frequency
        n_freq = len(frequencies)
        p0_all = np.zeros((n_freq, self.naux, self.naux), dtype=np.complex128)

        def compute_p0_single(omega):
            """Compute P0 for single frequency."""
            p0 = np.zeros((self.naux, self.naux), dtype=np.complex128)

            for i in range(n_occ):
                for a in range(n_vir):
                    de = e_virt[a] - e_occ[i]
                    # Resonant term
                    denom1 = omega - de + 1j * eta
                    # Anti-resonant term
                    denom2 = omega + de + 1j * eta

                    # P0(P,Q) = 2 * sum_ia (ia|P) * (ia|Q) * [1/(ω-Δε) - 1/(ω+Δε)]
                    vec = df_ia[i, a, :]
                    contribution = 2.0 * np.outer(vec, vec) * (1.0/denom1 - 1.0/denom2)
                    p0 += contribution

            return p0

        if self.parallel and self.executor and n_freq > 4:
            # Parallel computation
            futures = {self.executor.submit(compute_p0_single, omega): i
                      for i, omega in enumerate(frequencies)}

            for future in as_completed(futures):
                idx = futures[future]
                p0_all[idx] = future.result()
        else:
            # Serial computation
            for i, omega in enumerate(frequencies):
                p0_all[i] = compute_p0_single(omega)

        return p0_all

    def build_symmetrized_df_optimized(
        self,
        mf: 'scf.hf.SCF'
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Build symmetrized DF tensors with optimizations.

        Args:
            mf: PySCF SCF object

        Returns:
            Tuple of (b_ia_sym, vsqrt)
        """
        # Get V^(1/2) from cache
        vsqrt = self.get_vsqrt()

        # Get MO-transformed tensors
        j3c_ia, _, _ = self.transform_to_mo_batch(mf)[0]

        # Symmetrize using optimized matrix multiplication
        # B' = B @ V^(1/2)
        # Use BLAS for efficiency
        b_ia_sym = np.dot(j3c_ia, vsqrt)

        logger.info(f"Built symmetrized DF tensors: shape={b_ia_sym.shape}")

        return b_ia_sym, vsqrt

    def cleanup(self):
        """Clean up resources."""
        # Clear caches
        self.cache.clear()
        self._metric_cache = None
        self._vsqrt_cache = None
        self.get_3center_ao.cache_clear()
        self.get_2center_metric.cache_clear()
        self.get_vsqrt.cache_clear()

        # Close memory-mapped file if used
        if self._cderi_file and os.path.exists(self._cderi_file):
            if hasattr(self._cderi, 'close'):
                self._cderi.close()
            os.unlink(self._cderi_file)

        # Shutdown thread pool
        if self.executor:
            self.executor.shutdown(wait=False)

        logger.debug("Cleaned up OptimizedPySCFDFAdapter resources")

    def __del__(self):
        """Destructor to ensure cleanup."""
        self.cleanup()


def create_optimized_df_adapter(
    mol: 'gto.Mole',
    auxbasis: str = 'def2-tzvp-jkfit',
    **kwargs
) -> OptimizedPySCFDFAdapter:
    """
    Factory function to create optimized DF adapter.

    Args:
        mol: PySCF molecule object
        auxbasis: Auxiliary basis set name
        **kwargs: Additional arguments for OptimizedPySCFDFAdapter

    Returns:
        OptimizedPySCFDFAdapter instance
    """
    if not PYSCF_AVAILABLE:
        raise RuntimeError("PySCF is required for this function")

    # Auto-detect if memory mapping should be used
    estimated_memory = mol.nao ** 2 * mol.nao * 3 * 8 / 1e6  # Rough estimate in MB
    if 'use_mmap' not in kwargs and estimated_memory > 8000:
        kwargs['use_mmap'] = True
        logger.info(f"Auto-enabling memory mapping for large system (estimated {estimated_memory:.0f} MB)")

    return OptimizedPySCFDFAdapter(mol, auxbasis, **kwargs)


def benchmark_pyscf_optimizations(
    mol: 'gto.Mole',
    mf: 'scf.hf.SCF',
    auxbasis: str = 'def2-tzvp-jkfit'
) -> Dict[str, float]:
    """
    Benchmark optimized vs standard PySCF adapter.

    Returns:
        Dictionary with timing comparisons
    """
    import time
    from .pyscf_df_adapter import PySCFDFAdapter

    results = {}

    # Benchmark standard adapter
    logger.info("Benchmarking standard adapter...")
    standard_adapter = PySCFDFAdapter(mol, auxbasis)

    start = time.perf_counter()
    j3c_std = standard_adapter.get_3center_ao()
    results['standard_3center_time'] = time.perf_counter() - start

    start = time.perf_counter()
    j3c_ia_std, _, _ = standard_adapter.transform_to_mo(mf)
    results['standard_transform_time'] = time.perf_counter() - start

    # Benchmark optimized adapter
    logger.info("Benchmarking optimized adapter...")
    opt_adapter = OptimizedPySCFDFAdapter(mol, auxbasis, use_mmap=False, parallel=True)

    start = time.perf_counter()
    j3c_opt = opt_adapter.get_3center_ao()
    results['optimized_3center_time'] = time.perf_counter() - start

    start = time.perf_counter()
    j3c_ia_opt, _, _ = opt_adapter.transform_to_mo_batch(mf)[0]
    results['optimized_transform_time'] = time.perf_counter() - start

    # Verify correctness
    error_3c = np.abs(j3c_std - j3c_opt).max()
    error_ia = np.abs(j3c_ia_std - j3c_ia_opt).max()

    results['error_3center'] = error_3c
    results['error_transform'] = error_ia

    # Calculate speedups
    results['speedup_3center'] = results['standard_3center_time'] / results['optimized_3center_time']
    results['speedup_transform'] = results['standard_transform_time'] / results['optimized_transform_time']

    # Cleanup
    opt_adapter.cleanup()

    logger.info(f"Benchmark results:")
    logger.info(f"  3-center retrieval: {results['speedup_3center']:.2f}x speedup")
    logger.info(f"  MO transformation: {results['speedup_transform']:.2f}x speedup")
    logger.info(f"  Max errors: 3-center={error_3c:.2e}, transform={error_ia:.2e}")

    return results


# Example usage
if __name__ == "__main__" and PYSCF_AVAILABLE:
    logging.basicConfig(level=logging.INFO)

    # Create test molecule
    mol = gto.M(
        atom='O 0 0 0; H 0 0 1; H 0 1 0',
        basis='cc-pvdz',
        verbose=0
    )

    # Run HF
    mf = scf.RHF(mol)
    mf.kernel()

    # Create optimized adapter
    adapter = create_optimized_df_adapter(mol, 'cc-pvdz-jkfit')

    # Test batch transformation
    print("Testing batch MO transformation...")
    results = adapter.transform_to_mo_batch(mf)
    j3c_ia, j3c_ij, j3c_ab = results[0]
    print(f"Transformed tensors: ia={j3c_ia.shape}, ij={j3c_ij.shape}, ab={j3c_ab.shape}")

    # Test polarizability batch computation
    print("\nTesting batch polarizability computation...")
    frequencies = np.linspace(0, 10, 5)
    p0_batch = adapter.compute_polarizability_batch(frequencies, mf)
    print(f"P0 batch shape: {p0_batch.shape}")

    # Run benchmarks
    print("\nRunning benchmarks...")
    benchmark_results = benchmark_pyscf_optimizations(mol, mf, 'cc-pvdz-jkfit')

    # Cleanup
    adapter.cleanup()