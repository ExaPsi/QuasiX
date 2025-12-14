"""Batch processing utilities for correlation self-energy calculations.

This module provides optimized batch processing capabilities for computing
correlation self-energy at multiple evaluation points efficiently.
"""

from typing import Tuple, Dict, List, Optional, Union
import numpy as np
from numpy.typing import NDArray
import threading
from concurrent.futures import ThreadPoolExecutor
import warnings

from .correlation import CDParams, compute_correlation_self_energy_cd


class BatchProcessor:
    """Efficient batch processor for correlation self-energy calculations.

    This class optimizes the computation of correlation self-energy for
    multiple evaluation points by:
    1. Reusing preprocessed data across evaluations
    2. Chunking evaluation points for optimal cache usage
    3. Optional parallel processing
    """

    def __init__(
        self,
        mo_energy: NDArray[np.float64],
        mo_occ: NDArray[np.float64],
        w_matrices_imag: NDArray[np.complex128],
        omega_imag: NDArray[np.float64],
        df_tensors: NDArray[np.float64],
        cd_params: Optional[CDParams] = None,
    ):
        """Initialize batch processor with system data.

        Parameters
        ----------
        mo_energy : ndarray
            Molecular orbital energies
        mo_occ : ndarray
            Orbital occupations
        w_matrices_imag : ndarray
            Screened interaction on imaginary axis
        omega_imag : ndarray
            Imaginary frequency grid
        df_tensors : ndarray
            Density fitting tensors
        cd_params : CDParams, optional
            Contour deformation parameters
        """
        # Store preprocessed data
        self.mo_energy = self._ensure_contiguous(mo_energy, np.float64)
        self.mo_occ = self._ensure_contiguous(mo_occ, np.float64)
        self.omega_imag = self._ensure_contiguous(omega_imag, np.float64)
        self.df_tensors = self._ensure_contiguous(df_tensors, np.float64)

        # Preprocess W matrices once
        if w_matrices_imag.flags['C_CONTIGUOUS']:
            self.w_matrices = w_matrices_imag
        else:
            self.w_matrices = np.ascontiguousarray(w_matrices_imag, dtype=np.complex128)

        self.cd_params = cd_params or CDParams()

        # Validate once
        self._validate_inputs()

        # Cache for reusable computations
        self._cache = {}

    @staticmethod
    def _ensure_contiguous(arr: np.ndarray, dtype: np.dtype) -> np.ndarray:
        """Ensure array is C-contiguous with correct dtype."""
        if arr.flags['C_CONTIGUOUS'] and arr.dtype == dtype:
            return arr
        return np.ascontiguousarray(arr, dtype=dtype)

    def _validate_inputs(self):
        """Validate input dimensions once."""
        n_mo = len(self.mo_energy)
        n_freq, n_aux, n_aux2 = self.w_matrices.shape

        if n_aux != n_aux2:
            raise ValueError(f"W matrices must be square: got ({n_aux}, {n_aux2})")
        if len(self.mo_occ) != n_mo:
            raise ValueError("mo_energy and mo_occ must have same length")
        if self.df_tensors.shape != (n_mo, n_mo, n_aux):
            raise ValueError(f"df_tensors shape mismatch")
        if len(self.omega_imag) != n_freq:
            raise ValueError("omega_imag length mismatch")

    def compute_single(
        self,
        eval_omega: NDArray[np.float64]
    ) -> Tuple[NDArray[np.complex128], Dict]:
        """Compute correlation self-energy for a single set of evaluation points.

        Parameters
        ----------
        eval_omega : ndarray
            Evaluation frequencies

        Returns
        -------
        sigma_c : ndarray
            Correlation self-energy
        metadata : dict
            Calculation metadata
        """
        eval_omega = self._ensure_contiguous(eval_omega, np.float64)

        return compute_correlation_self_energy_cd(
            self.mo_energy,
            self.mo_occ,
            self.w_matrices,
            self.omega_imag,
            self.df_tensors,
            eval_omega,
            self.cd_params,
        )

    def compute_batch(
        self,
        eval_points_list: List[NDArray[np.float64]],
        parallel: bool = False,
        max_workers: Optional[int] = None,
    ) -> List[Tuple[NDArray[np.complex128], Dict]]:
        """Compute correlation self-energy for multiple evaluation point sets.

        Parameters
        ----------
        eval_points_list : list of ndarray
            List of evaluation frequency arrays
        parallel : bool
            Use parallel processing
        max_workers : int, optional
            Maximum number of parallel workers

        Returns
        -------
        results : list of (sigma_c, metadata) tuples
            Results for each evaluation point set
        """
        if not parallel:
            # Sequential processing
            results = []
            for eval_omega in eval_points_list:
                results.append(self.compute_single(eval_omega))
            return results

        # Parallel processing
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [
                executor.submit(self.compute_single, eval_omega)
                for eval_omega in eval_points_list
            ]
            results = [future.result() for future in futures]

        return results

    def compute_orbital_resolved(
        self,
        orbital_indices: Optional[List[int]] = None,
        energy_window: Optional[Tuple[float, float]] = None,
    ) -> Dict[int, Tuple[NDArray[np.complex128], Dict]]:
        """Compute correlation self-energy resolved by orbital.

        Parameters
        ----------
        orbital_indices : list of int, optional
            Specific orbital indices to compute (default: all)
        energy_window : tuple of float, optional
            Energy window (E_min, E_max) for evaluation

        Returns
        -------
        results : dict
            Mapping from orbital index to (sigma_c, metadata)
        """
        if orbital_indices is None:
            orbital_indices = list(range(len(self.mo_energy)))

        results = {}

        for idx in orbital_indices:
            # Evaluate at orbital energy
            eval_omega = np.array([self.mo_energy[idx]])

            # Apply energy window if specified
            if energy_window is not None:
                e_min, e_max = energy_window
                if not (e_min <= eval_omega[0] <= e_max):
                    continue

            sigma_c, metadata = self.compute_single(eval_omega)

            # Extract only the relevant orbital component
            results[idx] = (sigma_c[:, idx], metadata)

        return results

    def compute_on_grid(
        self,
        omega_min: float,
        omega_max: float,
        n_points: int,
        grid_type: str = 'linear',
    ) -> Tuple[NDArray[np.float64], NDArray[np.complex128], Dict]:
        """Compute correlation self-energy on a frequency grid.

        Parameters
        ----------
        omega_min : float
            Minimum frequency
        omega_max : float
            Maximum frequency
        n_points : int
            Number of grid points
        grid_type : str
            Grid type: 'linear', 'log', or 'sinh'

        Returns
        -------
        omega_grid : ndarray
            Frequency grid
        sigma_c : ndarray
            Correlation self-energy on grid
        metadata : dict
            Calculation metadata
        """
        # Create frequency grid
        if grid_type == 'linear':
            omega_grid = np.linspace(omega_min, omega_max, n_points)
        elif grid_type == 'log':
            if omega_min <= 0:
                raise ValueError("Log grid requires positive frequencies")
            omega_grid = np.logspace(np.log10(omega_min), np.log10(omega_max), n_points)
        elif grid_type == 'sinh':
            # Sinh grid for better resolution near zero
            x = np.linspace(-1, 1, n_points)
            omega_grid = omega_min + (omega_max - omega_min) * (np.sinh(x) / np.sinh(1) + 1) / 2
        else:
            raise ValueError(f"Unknown grid type: {grid_type}")

        # Compute in chunks for better cache usage
        chunk_size = min(10, n_points)  # Process 10 points at a time
        all_sigma = []
        all_metadata = []

        for i in range(0, n_points, chunk_size):
            chunk_omega = omega_grid[i:i+chunk_size]
            sigma_c, metadata = self.compute_single(chunk_omega)
            all_sigma.append(sigma_c)
            all_metadata.append(metadata)

        # Concatenate results
        sigma_c_full = np.vstack(all_sigma)

        # Merge metadata (use last chunk's metadata as representative)
        merged_metadata = all_metadata[-1] if all_metadata else {}

        return omega_grid, sigma_c_full, merged_metadata


def compute_correlation_batch(
    mo_energy: NDArray[np.float64],
    mo_occ: NDArray[np.float64],
    w_matrices_imag: NDArray[np.complex128],
    omega_imag: NDArray[np.float64],
    df_tensors: NDArray[np.float64],
    eval_points_list: Union[List[NDArray[np.float64]], NDArray[np.float64]],
    cd_params: Optional[CDParams] = None,
    parallel: bool = False,
) -> Union[List[Tuple[NDArray[np.complex128], Dict]], Tuple[NDArray[np.complex128], Dict]]:
    """Convenience function for batch correlation self-energy calculation.

    Parameters
    ----------
    mo_energy : ndarray
        Molecular orbital energies
    mo_occ : ndarray
        Orbital occupations
    w_matrices_imag : ndarray
        Screened interaction matrices
    omega_imag : ndarray
        Imaginary frequency grid
    df_tensors : ndarray
        Density fitting tensors
    eval_points_list : list of ndarray or ndarray
        Evaluation points (single array or list of arrays)
    cd_params : CDParams, optional
        Contour deformation parameters
    parallel : bool
        Use parallel processing

    Returns
    -------
    results : list of tuples or single tuple
        Correlation self-energy and metadata
    """
    processor = BatchProcessor(
        mo_energy, mo_occ, w_matrices_imag,
        omega_imag, df_tensors, cd_params
    )

    # Handle single array input
    if isinstance(eval_points_list, np.ndarray):
        return processor.compute_single(eval_points_list)

    # Handle list input
    return processor.compute_batch(eval_points_list, parallel=parallel)


def benchmark_batch_processing():
    """Benchmark batch processing performance."""
    import time

    print("Benchmarking Batch Processing for Correlation Self-Energy")
    print("=" * 60)

    # Create test system
    n_mo, n_aux, n_freq = 30, 100, 20
    mo_energy = np.random.randn(n_mo).astype(np.float64) * 10
    mo_occ = np.zeros(n_mo, dtype=np.float64)
    mo_occ[:n_mo//2] = 1.0

    w_matrices = np.random.randn(n_freq, n_aux, n_aux) + \
                 1j * np.random.randn(n_freq, n_aux, n_aux)
    w_matrices = 0.5 * (w_matrices + np.conj(np.transpose(w_matrices, (0, 2, 1))))

    omega_imag = np.logspace(-1, 1, n_freq, dtype=np.float64)
    df_tensors = np.random.randn(n_mo, n_mo, n_aux).astype(np.float64) * 0.1

    # Test different batch sizes
    batch_sizes = [1, 5, 10, 20]

    print(f"System: {n_mo} MOs, {n_aux} aux basis, {n_freq} frequencies\n")

    for batch_size in batch_sizes:
        # Create evaluation points
        eval_points_list = [
            np.random.randn(5).astype(np.float64)
            for _ in range(batch_size)
        ]

        # Time sequential processing
        processor = BatchProcessor(
            mo_energy, mo_occ, w_matrices,
            omega_imag, df_tensors
        )

        start = time.perf_counter()
        results_seq = processor.compute_batch(eval_points_list, parallel=False)
        seq_time = time.perf_counter() - start

        # Time parallel processing
        start = time.perf_counter()
        results_par = processor.compute_batch(eval_points_list, parallel=True, max_workers=4)
        par_time = time.perf_counter() - start

        # Verify results match
        for r_seq, r_par in zip(results_seq, results_par):
            np.testing.assert_allclose(r_seq[0], r_par[0], rtol=1e-10)

        print(f"Batch size {batch_size:2d}: Sequential {seq_time*1000:6.2f} ms, "
              f"Parallel {par_time*1000:6.2f} ms, Speedup {seq_time/par_time:.2f}x")

    # Test grid computation
    print("\n" + "-"*40)
    print("Testing grid computation:")

    processor = BatchProcessor(
        mo_energy, mo_occ, w_matrices,
        omega_imag, df_tensors
    )

    start = time.perf_counter()
    omega_grid, sigma_grid, metadata = processor.compute_on_grid(
        omega_min=-10.0, omega_max=10.0, n_points=50, grid_type='linear'
    )
    grid_time = time.perf_counter() - start

    print(f"Grid with 50 points: {grid_time*1000:.2f} ms")
    print(f"Grid shape: {sigma_grid.shape}")
    print(f"Throughput: {sigma_grid.size / grid_time:.0f} elements/sec")


if __name__ == "__main__":
    benchmark_batch_processing()