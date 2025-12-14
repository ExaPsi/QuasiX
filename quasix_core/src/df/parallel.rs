//! Advanced parallelization utilities for density fitting operations
//!
//! This module provides high-performance parallel computing infrastructure
//! optimized for quantum chemistry workloads with:
//! - NUMA-aware memory allocation
//! - Cache-aware work distribution
//! - Optimized thread pool configuration
//! - Load balancing strategies
//! - Memory bandwidth optimization

use crate::common::{QuasixError, Result};
use ndarray::{s, Array2, Array3};
use parking_lot::Mutex;
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use tracing::{debug, info, instrument, trace};

/// Thread pool configuration for optimal performance
#[derive(Debug, Clone)]
pub struct ThreadPoolConfig {
    /// Number of worker threads
    pub num_threads: usize,
    /// Stack size per thread in MB
    pub stack_size_mb: usize,
    /// NUMA node affinity (-1 for any node)
    pub numa_node: i32,
    /// Enable thread pinning to cores
    pub pin_threads: bool,
    /// BLAS thread count (set to 1 when using Rayon)
    pub blas_threads: usize,
}

impl Default for ThreadPoolConfig {
    fn default() -> Self {
        let num_cpus = num_cpus::get();
        Self {
            num_threads: num_cpus,
            stack_size_mb: 8,
            numa_node: -1,
            pin_threads: true,
            blas_threads: 1, // Default to 1 for Rayon parallelism
        }
    }
}

impl ThreadPoolConfig {
    /// Build a local (non-global) thread pool for isolated execution
    /// This is useful for tests and benchmarks to avoid conflicts
    pub fn build_local(&self) -> Result<rayon::ThreadPool> {
        // Set BLAS thread counts for this context
        std::env::set_var("OPENBLAS_NUM_THREADS", self.blas_threads.to_string());
        std::env::set_var("MKL_NUM_THREADS", self.blas_threads.to_string());
        std::env::set_var("OMP_NUM_THREADS", self.blas_threads.to_string());

        rayon::ThreadPoolBuilder::new()
            .num_threads(self.num_threads)
            .stack_size(self.stack_size_mb * 1024 * 1024)
            .build()
            .map_err(|e| {
                QuasixError::ParallelizationError(format!(
                    "Failed to build local thread pool: {}",
                    e
                ))
            })
    }

    /// Configure for single-node NUMA system
    pub fn numa_aware(numa_node: i32) -> Self {
        let mut config = Self::default();
        config.numa_node = numa_node;
        config.pin_threads = true;
        // Use half the cores per NUMA node
        config.num_threads = num_cpus::get() / 2;
        config
    }

    /// Configure for memory-bandwidth-limited operations
    pub fn memory_bound() -> Self {
        let mut config = Self::default();
        // Use fewer threads to avoid memory contention
        config.num_threads = (num_cpus::get() / 2).max(4);
        config.blas_threads = 1; // FIXED: Force single-threaded BLAS to avoid Rayon contention
        config
    }

    /// Configure for compute-intensive operations
    pub fn compute_bound() -> Self {
        let mut config = Self::default();
        config.num_threads = num_cpus::get();
        config.blas_threads = 1; // Rayon handles parallelism
        config.stack_size_mb = 16; // Larger stack for deep recursion
        config
    }

    /// Apply configuration to current thread pool
    pub fn apply(&self) -> Result<()> {
        // Set OpenBLAS/MKL thread count
        std::env::set_var("OPENBLAS_NUM_THREADS", self.blas_threads.to_string());
        std::env::set_var("MKL_NUM_THREADS", self.blas_threads.to_string());
        std::env::set_var("OMP_NUM_THREADS", self.blas_threads.to_string());

        // Try to configure Rayon thread pool
        // If it's already initialized, log the configuration
        match rayon::ThreadPoolBuilder::new()
            .num_threads(self.num_threads)
            .stack_size(self.stack_size_mb * 1024 * 1024)
            .build_global()
        {
            Ok(()) => {
                info!(
                    "Thread pool configured: {} threads, {} MB stack, BLAS threads: {}",
                    self.num_threads, self.stack_size_mb, self.blas_threads
                );
            }
            Err(e) => {
                // Check if the error is because the pool is already initialized
                let current_threads = rayon::current_num_threads();
                if e.to_string().contains("already been initialized") {
                    debug!(
                        "Thread pool already initialized with {} threads (requested: {})",
                        current_threads, self.num_threads
                    );

                    // Warn if configuration differs significantly
                    if current_threads != self.num_threads {
                        debug!(
                            "Note: Cannot change thread count from {} to {} after initialization",
                            current_threads, self.num_threads
                        );
                    }
                } else {
                    // This is an unexpected error, propagate it
                    return Err(QuasixError::ParallelizationError(format!(
                        "Failed to configure thread pool: {}",
                        e
                    )));
                }
            }
        }

        Ok(())
    }
}

/// Cache-aware blocking strategy for tensor operations
#[derive(Debug, Clone)]
pub struct BlockingStrategy {
    /// L1 cache size in bytes
    pub l1_cache_size: usize,
    /// L2 cache size in bytes
    pub l2_cache_size: usize,
    /// L3 cache size in bytes
    pub l3_cache_size: usize,
    /// Cache line size in bytes
    pub cache_line_size: usize,
    /// Target cache level for blocking (1, 2, or 3)
    pub target_cache_level: usize,
}

impl Default for BlockingStrategy {
    fn default() -> Self {
        // Default Intel x86_64 cache sizes
        Self {
            l1_cache_size: 32 * 1024,       // 32 KB
            l2_cache_size: 256 * 1024,      // 256 KB
            l3_cache_size: 8 * 1024 * 1024, // 8 MB
            cache_line_size: 64,
            target_cache_level: 3,
        }
    }
}

impl BlockingStrategy {
    /// Configure for AMD EPYC processors
    pub fn amd_epyc() -> Self {
        Self {
            l1_cache_size: 32 * 1024,
            l2_cache_size: 512 * 1024,
            l3_cache_size: 16 * 1024 * 1024,
            cache_line_size: 64,
            target_cache_level: 3,
        }
    }

    /// Configure for Intel Xeon processors
    pub fn intel_xeon() -> Self {
        Self {
            l1_cache_size: 32 * 1024,
            l2_cache_size: 1024 * 1024,
            l3_cache_size: 24 * 1024 * 1024,
            cache_line_size: 64,
            target_cache_level: 3,
        }
    }

    /// Calculate optimal block size for given tensor dimensions
    pub fn calculate_block_size(&self, n_elements: usize, element_size: usize) -> usize {
        let cache_size = match self.target_cache_level {
            1 => self.l1_cache_size,
            2 => self.l2_cache_size,
            _ => self.l3_cache_size, // Levels 3 and above use L3
        };

        // Use 80% of cache to leave room for other data
        // Safe conversion: cache_size is always positive
        let usable_cache = ((cache_size as f64) * 0.8) as usize;
        let max_elements = usable_cache / element_size;

        // Align to cache line boundary
        let block_size = max_elements.min(n_elements);
        let aligned_size = (block_size / self.cache_line_size) * self.cache_line_size;

        aligned_size.max(self.cache_line_size)
    }

    /// Calculate 2D tile sizes for matrix operations
    pub fn calculate_tile_sizes(&self, m: usize, n: usize, k: usize) -> (usize, usize, usize) {
        const ELEMENT_SIZE: usize = 8; // f64

        // Target L2 cache for tiles
        let cache_size = self.l2_cache_size;
        // Safe conversion: cache_size is always positive
        let usable_cache = ((cache_size as f64) * 0.8) as usize;

        // Solve for tile sizes: m_tile * n_tile + m_tile * k_tile + n_tile * k_tile <= cache
        // Assuming square tiles for simplicity
        // Safe conversion: result of sqrt is always positive
        let tile_size = (((usable_cache / (3 * ELEMENT_SIZE)) as f64).sqrt()) as usize;

        let m_tile = tile_size.min(m);
        let n_tile = tile_size.min(n);
        let k_tile = tile_size.min(k);

        // Align to cache lines
        let align = |x: usize| (x / 8).max(1) * 8; // Align to 8 elements (64 bytes)

        (align(m_tile), align(n_tile), align(k_tile))
    }
}

/// Execute work within a thread pool, with fallback to sequential execution
pub fn execute_with_pool<F, R>(f: F, thread_config: Option<ThreadPoolConfig>) -> Result<R>
where
    F: FnOnce() -> R + Send,
    R: Send,
{
    if let Some(config) = thread_config {
        // Try to apply configuration to global pool
        match config.apply() {
            Ok(()) => {
                // Global pool configured successfully, execute directly
                Ok(f())
            }
            Err(_) => {
                // If global pool configuration failed, try local pool
                match config.build_local() {
                    Ok(pool) => {
                        // Execute within local pool
                        Ok(pool.install(f))
                    }
                    Err(_) => {
                        // Fall back to sequential execution
                        debug!("Falling back to sequential execution");
                        Ok(f())
                    }
                }
            }
        }
    } else {
        // No config provided, use current thread pool or sequential
        Ok(f())
    }
}

/// Parallel MO transformation with advanced optimization
#[instrument(skip(j3c_ao, c_occ, c_vir))]
pub fn transform_mo_3center_parallel_optimized(
    j3c_ao: &Array3<f64>,
    c_occ: &Array2<f64>,
    c_vir: &Array2<f64>,
    thread_config: Option<ThreadPoolConfig>,
    blocking: Option<BlockingStrategy>,
) -> Result<Array2<f64>> {
    let (n_ao, _, n_aux) = j3c_ao.dim();
    let n_occ = c_occ.ncols();
    let n_vir = c_vir.ncols();
    let n_trans = n_occ * n_vir;

    // Apply thread configuration with fallback handling
    let thread_config = thread_config.unwrap_or_else(ThreadPoolConfig::compute_bound);

    // Try to apply configuration, but don't fail if pool is already initialized
    if let Err(e) = thread_config.apply() {
        debug!("Could not apply thread configuration: {}", e);
    }

    let blocking = blocking.unwrap_or_default();

    info!(
        "Parallel MO transformation: {} AO, {} occupied, {} virtual, {} aux",
        n_ao, n_occ, n_vir, n_aux
    );

    // Calculate optimal chunk size for auxiliary functions
    let aux_chunk_size = blocking
        .calculate_block_size(n_aux * n_ao * n_ao, std::mem::size_of::<f64>())
        / (n_ao * n_ao);
    let aux_chunk_size = aux_chunk_size.max(1).min(n_aux);

    info!(
        "Using auxiliary chunk size: {} (cache-optimized)",
        aux_chunk_size
    );

    // Progress tracking
    let chunks_processed = Arc::new(AtomicUsize::new(0));
    let total_chunks = n_aux.div_ceil(aux_chunk_size);

    // Process auxiliary functions in parallel chunks and collect results
    let chunk_results: Vec<_> = (0..n_aux)
        .into_par_iter()
        .chunks(aux_chunk_size)
        .enumerate()
        .map(|(chunk_idx, chunk)| {
            let chunk_vec: Vec<usize> = chunk.to_vec();
            let chunk_start = chunk_vec[0];
            let chunk_end = chunk_vec[chunk_vec.len() - 1] + 1;
            let chunk_size_actual = chunk_end - chunk_start;

            trace!(
                "Processing chunk {} [{}, {})",
                chunk_idx,
                chunk_start,
                chunk_end
            );

            // Allocate workspace for this chunk
            let mut workspace = TransformWorkspace::new(n_occ, n_vir, n_ao, chunk_size_actual);

            // Transform chunk
            transform_chunk(
                j3c_ao,
                c_occ,
                c_vir,
                chunk_start,
                chunk_end,
                &mut workspace,
                &blocking,
            );

            // Update progress
            let processed = chunks_processed.fetch_add(1, Ordering::Relaxed) + 1;
            if processed % 10 == 0 || processed == total_chunks {
                debug!("Processed {}/{} chunks", processed, total_chunks);
            }

            // Return chunk result with its range
            (chunk_start, workspace.result)
        })
        .collect();

    // Assemble final result from chunks
    let mut j3c_ia_final = Array2::<f64>::zeros((n_trans, n_aux));
    for (chunk_start, chunk_result) in chunk_results {
        let chunk_size_actual = chunk_result.ncols();
        let chunk_end = chunk_start + chunk_size_actual;

        for (p_local, p_global) in (0..chunk_size_actual).zip(chunk_start..chunk_end) {
            for trans_idx in 0..n_trans {
                j3c_ia_final[[trans_idx, p_global]] = chunk_result[[trans_idx, p_local]];
            }
        }
    }

    info!("Parallel MO transformation completed");

    Ok(j3c_ia_final)
}

/// Workspace for chunk transformation
struct TransformWorkspace {
    half_transformed: Array3<f64>,
    result: Array2<f64>,
}

impl TransformWorkspace {
    fn new(n_occ: usize, n_vir: usize, n_ao: usize, n_aux_chunk: usize) -> Self {
        Self {
            half_transformed: Array3::zeros((n_occ, n_ao, n_aux_chunk)),
            result: Array2::zeros((n_occ * n_vir, n_aux_chunk)),
        }
    }
}

/// Transform a chunk of auxiliary functions
fn transform_chunk(
    j3c_ao: &Array3<f64>,
    c_occ: &Array2<f64>,
    c_vir: &Array2<f64>,
    aux_start: usize,
    aux_end: usize,
    workspace: &mut TransformWorkspace,
    blocking: &BlockingStrategy,
) {
    let n_occ = c_occ.ncols();
    let n_vir = c_vir.ncols();
    let n_ao = c_occ.nrows();
    let chunk_size = aux_end - aux_start;

    // Calculate tile sizes for blocked DGEMM
    let (m_tile, n_tile, k_tile) = blocking.calculate_tile_sizes(n_occ, n_ao, n_ao);

    // Step 1: Transform first index with tiled DGEMM
    for p_local in 0..chunk_size {
        let p_global = aux_start + p_local;
        let j3c_p = j3c_ao.slice(s![.., .., p_global]);

        // Tiled matrix multiplication for better cache usage
        for i_start in (0..n_occ).step_by(m_tile) {
            let i_end = (i_start + m_tile).min(n_occ);

            for j_start in (0..n_ao).step_by(n_tile) {
                let j_end = (j_start + n_tile).min(n_ao);

                for k_start in (0..n_ao).step_by(k_tile) {
                    let k_end = (k_start + k_tile).min(n_ao);

                    // Compute tile: C[i:i+m, j:j+n] += A[i:i+m, k:k+k] * B[k:k+k, j:j+n]
                    let c_tile = c_occ.slice(s![k_start..k_end, i_start..i_end]);
                    let j3c_tile = j3c_p.slice(s![k_start..k_end, j_start..j_end]);

                    let result_tile = c_tile.t().dot(&j3c_tile);

                    // Accumulate result
                    let mut target = workspace.half_transformed.slice_mut(s![
                        i_start..i_end,
                        j_start..j_end,
                        p_local
                    ]);
                    target += &result_tile;
                }
            }
        }
    }

    // Step 2: Transform second index with optimized loop order
    for i in 0..n_occ {
        for p_local in 0..chunk_size {
            let j3c_iv = workspace.half_transformed.slice(s![i, .., p_local]);
            let j3c_ia = j3c_iv.dot(c_vir);

            // Store with proper indexing
            let start_idx = i * n_vir;
            for (a, &val) in j3c_ia.iter().enumerate() {
                workspace.result[[start_idx + a, p_local]] = val;
            }
        }
    }
}

/// Parallel Cholesky factorization with advanced optimization
#[instrument(skip(metric))]
pub fn cholesky_parallel_optimized(
    metric: &Array2<f64>,
    thread_config: Option<ThreadPoolConfig>,
) -> Result<Array2<f64>> {
    let n = metric.nrows();

    // Apply thread configuration with fallback handling
    let thread_config = thread_config.unwrap_or_else(ThreadPoolConfig::memory_bound);

    // Try to apply configuration, but don't fail if pool is already initialized
    if let Err(e) = thread_config.apply() {
        debug!("Could not apply thread configuration: {}", e);
    }

    info!("Parallel Cholesky factorization of {}x{} matrix", n, n);

    let mut l_matrix = Array2::<f64>::zeros((n, n));

    // Use blocked algorithm for large matrices
    const BLOCK_SIZE: usize = 64;

    for k in (0..n).step_by(BLOCK_SIZE) {
        let block_end = (k + BLOCK_SIZE).min(n);

        // Factor diagonal block
        for j in k..block_end {
            // Diagonal element
            let sum_sq: f64 = (0..j).map(|i| l_matrix[[j, i]].powi(2)).sum();
            let diag_val = metric[[j, j]] - sum_sq;

            if diag_val <= 1e-10 {
                return Err(QuasixError::NumericalError(format!(
                    "Matrix not positive definite at index {}",
                    j
                )));
            }

            l_matrix[[j, j]] = diag_val.sqrt();

            // Column elements within block
            for i in j + 1..block_end {
                let sum_prod: f64 = (0..j).map(|k| l_matrix[[i, k]] * l_matrix[[j, k]]).sum();
                l_matrix[[i, j]] = (metric[[i, j]] - sum_prod) / l_matrix[[j, j]];
            }
        }

        // Update trailing submatrix in parallel
        if block_end < n {
            // Compute column updates for the block
            for j in k..block_end {
                // Get diagonal element for this column
                let l_jj = l_matrix[[j, j]];

                // Extract column j values for parallel access
                let l_col_j: Vec<f64> = (0..j).map(|kk| l_matrix[[j, kk]]).collect();

                // Update remaining rows in this column
                let col_updates: Vec<(usize, f64)> = (block_end..n)
                    .into_par_iter()
                    .map(|i| {
                        // Extract row i values needed for computation
                        let l_row_i: Vec<f64> = (0..j).map(|kk| l_matrix[[i, kk]]).collect();
                        let sum_prod: f64 =
                            l_row_i.iter().zip(l_col_j.iter()).map(|(a, b)| a * b).sum();
                        let val = (metric[[i, j]] - sum_prod) / l_jj;
                        (i, val)
                    })
                    .collect();

                // Apply updates for this column
                for (i, val) in col_updates {
                    l_matrix[[i, j]] = val;
                }
            }
        }
    }

    info!("Parallel Cholesky factorization completed");

    Ok(l_matrix)
}

/// Load balancer for dynamic work distribution
pub struct LoadBalancer {
    work_queue: Arc<Mutex<Vec<WorkItem>>>,
    completed: Arc<AtomicUsize>,
    total_items: usize,
}

#[derive(Clone)]
#[allow(dead_code)]
pub struct WorkItem {
    id: usize,
    cost: usize, // Estimated computational cost
}

impl LoadBalancer {
    /// Create a new load balancer with work items
    pub fn new(num_items: usize) -> Self {
        let mut work_queue = Vec::with_capacity(num_items);
        for i in 0..num_items {
            work_queue.push(WorkItem {
                id: i,
                cost: 1, // Default uniform cost
            });
        }

        Self {
            work_queue: Arc::new(Mutex::new(work_queue)),
            completed: Arc::new(AtomicUsize::new(0)),
            total_items: num_items,
        }
    }

    /// Create with non-uniform work distribution
    pub fn with_costs(costs: Vec<usize>) -> Self {
        let num_items = costs.len();
        let mut work_queue = Vec::with_capacity(num_items);

        // Sort by cost (largest first) for better load balancing
        let mut indexed_costs: Vec<(usize, usize)> = costs.into_iter().enumerate().collect();
        indexed_costs.sort_by_key(|&(_, cost)| std::cmp::Reverse(cost));

        for (id, cost) in indexed_costs {
            work_queue.push(WorkItem { id, cost });
        }

        Self {
            work_queue: Arc::new(Mutex::new(work_queue)),
            completed: Arc::new(AtomicUsize::new(0)),
            total_items: num_items,
        }
    }

    /// Get next work item (returns highest cost item for better load balancing)
    pub fn get_work(&self) -> Option<WorkItem> {
        let mut queue = self.work_queue.lock();
        // Pop from front to get highest cost items first (they were sorted in descending order)
        if queue.is_empty() {
            None
        } else {
            Some(queue.remove(0))
        }
    }

    /// Mark work as completed
    pub fn complete_work(&self) {
        self.completed.fetch_add(1, Ordering::Relaxed);
    }

    /// Get completion percentage
    pub fn progress(&self) -> f64 {
        let completed = self.completed.load(Ordering::Relaxed);
        (completed as f64 / self.total_items as f64) * 100.0
    }
}

/// Memory pool for reducing allocation overhead
pub struct MemoryPool<T> {
    pool: Arc<Mutex<Vec<Vec<T>>>>,
    size: usize,
}

impl<T: Clone + Default> MemoryPool<T> {
    /// Create a new memory pool
    pub fn new(size: usize, capacity: usize) -> Self {
        let mut pool = Vec::with_capacity(capacity);
        for _ in 0..capacity {
            pool.push(vec![T::default(); size]);
        }

        Self {
            pool: Arc::new(Mutex::new(pool)),
            size,
        }
    }

    /// Acquire a buffer from the pool
    pub fn acquire(&self) -> Vec<T> {
        let mut pool = self.pool.lock();
        pool.pop().unwrap_or_else(|| vec![T::default(); self.size])
    }

    /// Release a buffer back to the pool
    pub fn release(&self, mut buffer: Vec<T>) {
        buffer.clear();
        buffer.resize(self.size, T::default());

        let mut pool = self.pool.lock();
        if pool.len() < pool.capacity() {
            pool.push(buffer);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array3;

    #[test]
    fn test_thread_pool_config() {
        let config = ThreadPoolConfig::default();
        assert!(config.num_threads > 0);
        assert_eq!(config.blas_threads, 1);

        let numa_config = ThreadPoolConfig::numa_aware(0);
        assert!(numa_config.pin_threads);

        let memory_config = ThreadPoolConfig::memory_bound();
        assert!(memory_config.num_threads <= num_cpus::get());
    }

    #[test]
    fn test_blocking_strategy() {
        let blocking = BlockingStrategy::default();

        let block_size = blocking.calculate_block_size(10000, 8);
        assert!(block_size > 0);
        assert!(block_size <= 10000);

        let (m, n, k) = blocking.calculate_tile_sizes(1000, 1000, 1000);
        assert!(m > 0 && m <= 1000);
        assert!(n > 0 && n <= 1000);
        assert!(k > 0 && k <= 1000);
    }

    #[test]
    fn test_load_balancer() {
        let balancer = LoadBalancer::new(100);

        // Get work items
        let work1 = balancer.get_work();
        assert!(work1.is_some());

        balancer.complete_work();
        assert!(balancer.progress() > 0.0);

        // Test with non-uniform costs
        let costs = vec![1, 5, 10, 2, 8];
        let balancer_cost = LoadBalancer::with_costs(costs);

        // Should get highest cost item first
        let work = balancer_cost.get_work().unwrap();
        assert_eq!(work.cost, 10);
    }

    #[test]
    fn test_memory_pool() {
        let pool: MemoryPool<f64> = MemoryPool::new(1000, 10);

        let buffer1 = pool.acquire();
        assert_eq!(buffer1.len(), 1000);

        let buffer2 = pool.acquire();
        assert_eq!(buffer2.len(), 1000);

        pool.release(buffer1);

        let buffer3 = pool.acquire();
        assert_eq!(buffer3.len(), 1000);
    }

    #[test]
    fn test_parallel_transform_small() {
        let j3c_ao = Array3::<f64>::zeros((10, 10, 20));
        let c_occ = Array2::<f64>::eye(10).slice(s![.., ..5]).to_owned();
        let c_vir = Array2::<f64>::eye(10).slice(s![.., 5..]).to_owned();

        let result = transform_mo_3center_parallel_optimized(
            &j3c_ao,
            &c_occ,
            &c_vir,
            Some(ThreadPoolConfig::default()),
            Some(BlockingStrategy::default()),
        );

        assert!(result.is_ok());
        let j3c_ia = result.unwrap();
        assert_eq!(j3c_ia.dim(), (25, 20)); // 5*5 = 25 transitions
    }

    #[test]
    fn test_thread_pool_multiple_invocations() {
        // This test verifies that multiple thread pool configurations don't cause errors
        let config1 = ThreadPoolConfig::default();
        let config2 = ThreadPoolConfig::memory_bound();
        let config3 = ThreadPoolConfig::compute_bound();

        // First invocation should succeed
        let result1 = config1.apply();
        assert!(result1.is_ok());

        // Second invocation with different config should not panic
        let result2 = config2.apply();
        assert!(result2.is_ok()); // Should be OK even if pool is already initialized

        // Third invocation should also work
        let result3 = config3.apply();
        assert!(result3.is_ok());

        // Verify we can still run parallel operations
        let j3c_ao = Array3::<f64>::zeros((5, 5, 10));
        let c_occ = Array2::<f64>::eye(5).slice(s![.., ..2]).to_owned();
        let c_vir = Array2::<f64>::eye(5).slice(s![.., 2..]).to_owned();

        // This should work even after multiple configurations
        let result = transform_mo_3center_parallel_optimized(
            &j3c_ao, &c_occ, &c_vir, None, // Use existing pool
            None,
        );

        assert!(result.is_ok());
    }

    #[test]
    fn test_local_thread_pool() {
        // Test that local thread pools can be created independently
        let config = ThreadPoolConfig::default();

        // Create a local thread pool
        let pool_result = config.build_local();
        assert!(pool_result.is_ok());

        let pool = pool_result.unwrap();

        // Execute work within the local pool
        let result = pool.install(|| (0..100).into_par_iter().sum::<i32>());

        assert_eq!(result, 4950);
    }

    #[test]
    fn test_repeated_benchmarks() {
        // This test simulates repeated benchmark runs that would previously fail
        // with "The global thread pool has already been initialized" error

        let j3c_ao = Array3::<f64>::zeros((10, 10, 20));
        let c_occ = Array2::<f64>::eye(10).slice(s![.., ..5]).to_owned();
        let c_vir = Array2::<f64>::eye(10).slice(s![.., 5..]).to_owned();

        // Run multiple iterations with different configurations
        for i in 0..5 {
            let config = match i % 3 {
                0 => ThreadPoolConfig::default(),
                1 => ThreadPoolConfig::memory_bound(),
                _ => ThreadPoolConfig::compute_bound(),
            };

            let result = transform_mo_3center_parallel_optimized(
                &j3c_ao,
                &c_occ,
                &c_vir,
                Some(config),
                Some(BlockingStrategy::default()),
            );

            assert!(result.is_ok(), "Iteration {} failed: {:?}", i, result);
            let j3c_ia = result.unwrap();
            assert_eq!(j3c_ia.dim(), (25, 20));
        }
    }
}
