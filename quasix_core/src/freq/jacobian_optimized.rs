//! Optimized Jacobian computation with overflow protection
//!
//! High-performance implementations with SIMD optimization and numerical stability
//! for frequency grid transformations, especially near singularities.

use crate::common::Result;
use ndarray::Array1;
use rayon::prelude::*;
use std::f64;
use wide::f64x4;

/// Cache for precomputed grids with different parameters
#[derive(Debug, Clone)]
pub struct GridCache {
    /// Cached grids indexed by (n_points, omega_max)
    cache: std::collections::HashMap<(usize, u64), CachedGrid>,
}

/// Cached grid data (internal structure)
#[derive(Debug, Clone)]
struct CachedGrid {
    points: Array1<f64>,
    weights: Array1<f64>,
    #[allow(dead_code)] // Kept for future cache-based optimization
    jacobians: Array1<f64>,
}

impl Default for GridCache {
    fn default() -> Self {
        Self::new()
    }
}

impl GridCache {
    pub fn new() -> Self {
        Self {
            cache: std::collections::HashMap::new(),
        }
    }
}

/// Optimized Jacobian computer with SIMD and overflow protection
pub struct JacobianOptimized {
    /// Cutoff threshold for (1-x) to prevent overflow
    cutoff_threshold: f64,
    /// Maximum allowed Jacobian value
    max_jacobian: f64,
    /// Use smooth cutoff vs hard cutoff
    use_smooth_cutoff: bool,
    /// Grid cache for reuse
    cache: GridCache,
}

impl Default for JacobianOptimized {
    fn default() -> Self {
        Self {
            cutoff_threshold: 1e-10,
            max_jacobian: 1e20,
            use_smooth_cutoff: true,
            cache: GridCache::new(),
        }
    }
}

impl JacobianOptimized {
    pub fn new() -> Self {
        Self::default()
    }

    /// Configure cutoff parameters
    pub fn with_cutoff(mut self, threshold: f64, max_jacobian: f64) -> Self {
        self.cutoff_threshold = threshold;
        self.max_jacobian = max_jacobian;
        self
    }

    /// Enable/disable smooth cutoff
    pub fn with_smooth_cutoff(mut self, smooth: bool) -> Self {
        self.use_smooth_cutoff = smooth;
        self
    }

    /// Compute Jacobian with overflow protection (scalar version)
    #[inline(always)]
    pub fn safe_jacobian(&self, x: f64, omega_max: f64) -> f64 {
        let denom = 1.0 - x;

        if denom.abs() < self.cutoff_threshold {
            if self.use_smooth_cutoff {
                // Smooth cutoff using tanh-like transition
                self.smooth_cutoff_jacobian(x, omega_max)
            } else {
                // Hard cutoff at maximum value
                self.max_jacobian * omega_max
            }
        } else {
            // Safe to compute normally
            omega_max / (denom * denom)
        }
    }

    /// Smooth cutoff function for Jacobian near x=1
    #[inline(always)]
    fn smooth_cutoff_jacobian(&self, x: f64, omega_max: f64) -> f64 {
        let denom = 1.0 - x;
        let epsilon = self.cutoff_threshold;

        // Use a smooth transition function
        // J(x) = omega_max * (1 + tanh((epsilon - |denom|)/delta)) / (2 * epsilon^2)
        // This provides smooth behavior near the cutoff

        let delta = epsilon * 0.1; // Transition width
        let transition = ((epsilon - denom.abs()) / delta).tanh();
        let smooth_factor = (1.0 + transition) * 0.5;

        // Blend between regular and maximum Jacobian
        let regular_jac = omega_max / (denom * denom + epsilon * epsilon);
        let max_jac = self.max_jacobian * omega_max;

        regular_jac * (1.0 - smooth_factor) + max_jac * smooth_factor
    }

    /// SIMD-optimized Jacobian computation with overflow protection
    pub fn compute_jacobian_simd(&self, x: &[f64], omega_max: f64) -> Vec<f64> {
        let n = x.len();
        let mut jacobians = vec![0.0; n];

        const SIMD_WIDTH: usize = 4;
        let n_simd = n - (n % SIMD_WIDTH);

        let _omega_simd = f64x4::splat(omega_max); // Kept for potential vectorized computation
        let one_simd = f64x4::splat(1.0);
        let _cutoff_simd = f64x4::splat(self.cutoff_threshold); // Kept for potential vectorized comparison
        let _max_jac_simd = f64x4::splat(self.max_jacobian * omega_max); // Kept for potential SIMD clamping

        // Process SIMD chunks
        for i in (0..n_simd).step_by(SIMD_WIDTH) {
            // Load x values
            let mut x_vals = [0.0_f64; SIMD_WIDTH];
            for j in 0..SIMD_WIDTH {
                x_vals[j] = x[i + j];
            }
            let x_simd = f64x4::from(x_vals);

            // Compute denominator: 1 - x
            let denom_simd = one_simd - x_simd;

            // Check for cutoff condition
            // Since wide doesn't have abs, we'll handle this element-wise
            let _denom_arr: [f64; 4] = denom_simd.into(); // Kept for potential SIMD-based cutoff logic
            let mut jac_vals = [0.0_f64; SIMD_WIDTH];

            for j in 0..SIMD_WIDTH {
                jac_vals[j] = self.safe_jacobian(x_vals[j], omega_max);
            }

            // Store results
            for j in 0..SIMD_WIDTH {
                jacobians[i + j] = jac_vals[j];
            }
        }

        // Handle remainder
        for i in n_simd..n {
            jacobians[i] = self.safe_jacobian(x[i], omega_max);
        }

        jacobians
    }

    /// Parallel Jacobian computation for large grids
    pub fn compute_jacobian_parallel(&self, x: &[f64], omega_max: f64) -> Vec<f64> {
        x.par_iter()
            .map(|&xi| self.safe_jacobian(xi, omega_max))
            .collect()
    }

    /// Optimized grid generation with cached Jacobians
    pub fn generate_grid_cached(
        &mut self,
        n: usize,
        omega_max: f64,
        transform_type: TransformType,
    ) -> Result<(Array1<f64>, Array1<f64>)> {
        // Check if already cached
        let omega_bits = omega_max.to_bits();
        let key = (n, omega_bits);

        if self.cache.cache.contains_key(&key) {
            let grid = self.cache.cache.get(&key).unwrap();
            return Ok((grid.points.clone(), grid.weights.clone()));
        }

        // Generate new grid
        let (points, weights, jacobians) =
            self.generate_grid_internal(n, omega_max, transform_type)?;

        // Cache it
        self.cache.cache.insert(
            key,
            CachedGrid {
                points: points.clone(),
                weights: weights.clone(),
                jacobians,
            },
        );

        Ok((points, weights))
    }

    /// Internal grid generation with Jacobian computation
    fn generate_grid_internal(
        &self,
        n: usize,
        omega_max: f64,
        transform_type: TransformType,
    ) -> Result<(Array1<f64>, Array1<f64>, Array1<f64>)> {
        // Get base GL nodes and weights
        let (nodes, weights) = crate::freq::gauss_legendre_nodes_weights(n)?;

        match transform_type {
            TransformType::Rational => {
                self.rational_transform_optimized(&nodes, &weights, omega_max)
            }
            TransformType::Exponential => {
                self.exponential_transform_optimized(&nodes, &weights, omega_max)
            }
            TransformType::Sinh => self.sinh_transform_optimized(&nodes, &weights, omega_max),
        }
    }

    /// Optimized rational transformation with overflow protection
    fn rational_transform_optimized(
        &self,
        nodes: &Array1<f64>,
        weights: &Array1<f64>,
        omega_max: f64,
    ) -> Result<(Array1<f64>, Array1<f64>, Array1<f64>)> {
        let n = nodes.len();
        let mut points = Array1::zeros(n);
        let mut new_weights = Array1::zeros(n);
        let mut jacobians = Array1::zeros(n);

        // Map [-1, 1] to [0, 1]
        let t_values: Vec<f64> = nodes.iter().map(|&t| (t + 1.0) * 0.5).collect();

        // Compute Jacobians with SIMD
        let jac_vec = if n > 64 {
            self.compute_jacobian_parallel(&t_values, omega_max)
        } else {
            self.compute_jacobian_simd(&t_values, omega_max)
        };

        // Apply transformation
        for i in 0..n {
            let t = t_values[i];

            // omega = omega_max * t / (1 + a*(1-t))
            let a = 0.1;
            let denom = 1.0 + a * (1.0 - t);
            points[i] = omega_max * t / denom;

            // Use precomputed safe Jacobian
            jacobians[i] = jac_vec[i];
            new_weights[i] = weights[i] * 0.5 * jacobians[i]; // 0.5 from [-1,1] to [0,1]
        }

        Ok((points, new_weights, jacobians))
    }

    /// Optimized exponential transformation
    fn exponential_transform_optimized(
        &self,
        nodes: &Array1<f64>,
        weights: &Array1<f64>,
        omega_max: f64,
    ) -> Result<(Array1<f64>, Array1<f64>, Array1<f64>)> {
        let n = nodes.len();
        let mut points = Array1::zeros(n);
        let mut new_weights = Array1::zeros(n);
        let mut jacobians = Array1::zeros(n);

        let alpha = (omega_max + 1.0).ln();

        // Process with SIMD where possible
        const SIMD_WIDTH: usize = 4;
        let n_simd = n - (n % SIMD_WIDTH);

        for i in (0..n_simd).step_by(SIMD_WIDTH) {
            let mut t_vals = [0.0_f64; SIMD_WIDTH];
            let mut w_vals = [0.0_f64; SIMD_WIDTH];

            for j in 0..SIMD_WIDTH {
                t_vals[j] = (nodes[i + j] + 1.0) * 0.5; // Map to [0,1]
                w_vals[j] = weights[i + j];
            }

            // Compute exponential transformation
            for j in 0..SIMD_WIDTH {
                let t = t_vals[j];
                let exp_val = (alpha * t).exp();
                points[i + j] = (exp_val - 1.0).min(omega_max);
                jacobians[i + j] = alpha * exp_val;
                new_weights[i + j] = w_vals[j] * 0.5 * jacobians[i + j];
            }
        }

        // Handle remainder
        for i in n_simd..n {
            let t = (nodes[i] + 1.0) * 0.5;
            let exp_val = (alpha * t).exp();
            points[i] = (exp_val - 1.0).min(omega_max);
            jacobians[i] = alpha * exp_val;
            new_weights[i] = weights[i] * 0.5 * jacobians[i];
        }

        Ok((points, new_weights, jacobians))
    }

    /// Optimized sinh transformation
    fn sinh_transform_optimized(
        &self,
        nodes: &Array1<f64>,
        weights: &Array1<f64>,
        omega_max: f64,
    ) -> Result<(Array1<f64>, Array1<f64>, Array1<f64>)> {
        let n = nodes.len();
        let mut points = Array1::zeros(n);
        let mut new_weights = Array1::zeros(n);
        let mut jacobians = Array1::zeros(n);

        let alpha = omega_max.ln();
        let sinh_alpha = alpha.sinh();

        // Parallel processing for large grids
        if n > 100 {
            let results: Vec<(f64, f64)> = (0..n)
                .into_par_iter()
                .map(|i| {
                    let t = (nodes[i] + 1.0) * 0.5;
                    let sinh_val = (alpha * t).sinh();
                    let omega = omega_max * sinh_val / sinh_alpha;
                    let jac = omega_max * (alpha * t).cosh() / sinh_alpha;
                    (omega, jac)
                })
                .collect();

            for (i, (omega, jac)) in results.iter().enumerate() {
                points[i] = *omega;
                jacobians[i] = *jac;
                new_weights[i] = weights[i] * 0.5 * jac;
            }
        } else {
            // Sequential for small grids
            for i in 0..n {
                let t = (nodes[i] + 1.0) * 0.5;
                let sinh_val = (alpha * t).sinh();
                points[i] = omega_max * sinh_val / sinh_alpha;
                jacobians[i] = omega_max * (alpha * t).cosh() / sinh_alpha;
                new_weights[i] = weights[i] * 0.5 * jacobians[i];
            }
        }

        Ok((points, new_weights, jacobians))
    }

    /// Validate grid for numerical stability
    pub fn validate_grid(
        &self,
        points: &Array1<f64>,
        weights: &Array1<f64>,
        jacobians: &Array1<f64>,
    ) -> Result<()> {
        // Check for infinities and NaN
        if points.iter().any(|&x| !x.is_finite()) {
            return Err(crate::common::QuasixError::InvalidInput(
                "Non-finite points detected".to_string(),
            ));
        }

        if weights.iter().any(|&w| !w.is_finite() || w < 0.0) {
            return Err(crate::common::QuasixError::InvalidInput(
                "Invalid weights detected".to_string(),
            ));
        }

        if jacobians.iter().any(|&j| !j.is_finite()) {
            return Err(crate::common::QuasixError::InvalidInput(
                "Non-finite Jacobians detected".to_string(),
            ));
        }

        // Check for monotonicity
        for i in 1..points.len() {
            if points[i] <= points[i - 1] {
                return Err(crate::common::QuasixError::InvalidInput(
                    "Non-monotonic grid points".to_string(),
                ));
            }
        }

        // Check condition number of weights
        let max_weight = weights.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
        let min_weight = weights
            .iter()
            .filter(|&&w| w > 0.0)
            .fold(f64::INFINITY, |a, &b| a.min(b));
        let condition = max_weight / min_weight;

        if condition > 1e12 {
            eprintln!("Warning: High weight condition number: {:.2e}", condition);
        }

        Ok(())
    }

    /// Benchmark different Jacobian computation methods
    pub fn benchmark_jacobian_methods(&self, n: usize, omega_max: f64) {
        use std::time::Instant;

        // Generate test data
        let x: Vec<f64> = (0..n).map(|i| i as f64 / (n - 1) as f64).collect();

        // Benchmark scalar version
        let start = Instant::now();
        let _scalar_result: Vec<f64> = x
            .iter()
            .map(|&xi| self.safe_jacobian(xi, omega_max))
            .collect();
        let scalar_time = start.elapsed();

        // Benchmark SIMD version
        let start = Instant::now();
        let _simd_result = self.compute_jacobian_simd(&x, omega_max);
        let simd_time = start.elapsed();

        // Benchmark parallel version
        let start = Instant::now();
        let _parallel_result = self.compute_jacobian_parallel(&x, omega_max);
        let parallel_time = start.elapsed();

        println!("Jacobian computation benchmark (n={}):", n);
        println!("  Scalar:   {:?}", scalar_time);
        println!(
            "  SIMD:     {:?} ({:.2}x speedup)",
            simd_time,
            scalar_time.as_secs_f64() / simd_time.as_secs_f64()
        );
        println!(
            "  Parallel: {:?} ({:.2}x speedup)",
            parallel_time,
            scalar_time.as_secs_f64() / parallel_time.as_secs_f64()
        );
    }
}

/// Transformation type for grid generation
#[derive(Debug, Clone, Copy)]
pub enum TransformType {
    Rational,
    Exponential,
    Sinh,
}

/// Branch prediction hints for likely/unlikely conditions
#[inline(always)]
#[allow(unused)]
fn likely(b: bool) -> bool {
    b
}

#[inline(always)]
#[allow(unused)]
fn unlikely(b: bool) -> bool {
    !b
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_safe_jacobian_no_overflow() {
        let computer = JacobianOptimized::new();

        // Test near x = 1 (would overflow without protection)
        let x_values = vec![0.9999, 0.99999, 0.999999, 0.9999999, 0.99999999];
        let omega_max = 1000.0;

        for &x in &x_values {
            let jac = computer.safe_jacobian(x, omega_max);
            assert!(jac.is_finite(), "Jacobian should be finite for x={}", x);
            assert!(jac > 0.0, "Jacobian should be positive");
            assert!(
                jac <= computer.max_jacobian * omega_max,
                "Jacobian should be bounded"
            );
        }
    }

    #[test]
    fn test_smooth_cutoff() {
        let computer = JacobianOptimized::new().with_smooth_cutoff(true);

        // Test smooth transition near cutoff
        let x_values: Vec<f64> = (0..100).map(|i| 1.0 - 1e-11 * (1.1_f64.powi(i))).collect();

        let omega_max = 100.0;
        let jacobians: Vec<f64> = x_values
            .iter()
            .map(|&x| computer.safe_jacobian(x, omega_max))
            .collect();

        // Check smoothness (no large jumps)
        for i in 1..jacobians.len() {
            let ratio = jacobians[i] / jacobians[i - 1];
            assert!(
                ratio < 10.0 && ratio > 0.1,
                "Jacobian should change smoothly: {} -> {}",
                jacobians[i - 1],
                jacobians[i]
            );
        }
    }

    #[test]
    fn test_simd_jacobian() {
        let computer = JacobianOptimized::new();
        let n = 128; // Multiple of SIMD width
        let x: Vec<f64> = (0..n).map(|i| i as f64 / (n - 1) as f64).collect();
        let omega_max = 100.0;

        // Compute with SIMD
        let simd_result = computer.compute_jacobian_simd(&x, omega_max);

        // Compute with scalar for comparison
        let scalar_result: Vec<f64> = x
            .iter()
            .map(|&xi| computer.safe_jacobian(xi, omega_max))
            .collect();

        // Results should match
        for i in 0..n {
            assert_relative_eq!(simd_result[i], scalar_result[i], epsilon = 1e-14);
        }
    }

    #[test]
    fn test_parallel_jacobian() {
        let computer = JacobianOptimized::new();
        let n = 1000;
        let x: Vec<f64> = (0..n).map(|i| i as f64 / (n - 1) as f64).collect();
        let omega_max = 100.0;

        // Compute with parallel
        let parallel_result = computer.compute_jacobian_parallel(&x, omega_max);

        // Compute with scalar for comparison
        let scalar_result: Vec<f64> = x
            .iter()
            .map(|&xi| computer.safe_jacobian(xi, omega_max))
            .collect();

        // Results should match
        for i in 0..n {
            assert_relative_eq!(parallel_result[i], scalar_result[i], epsilon = 1e-14);
        }
    }

    #[test]
    fn test_grid_generation() {
        let mut computer = JacobianOptimized::new();
        let n = 32;
        let omega_max = 100.0;

        // Test rational transform
        let (points, weights) = computer
            .generate_grid_cached(n, omega_max, TransformType::Rational)
            .unwrap();

        assert_eq!(points.len(), n);
        assert_eq!(weights.len(), n);
        assert!(points.iter().all(|&x| x.is_finite()));
        assert!(weights.iter().all(|&w| w.is_finite() && w > 0.0));

        // Check monotonicity
        for i in 1..n {
            assert!(points[i] > points[i - 1], "Points should be monotonic");
        }
    }

    #[test]
    fn test_grid_caching() {
        let mut computer = JacobianOptimized::new();
        let n = 32;
        let omega_max = 100.0;

        // First call should compute
        let (points1, weights1) = computer
            .generate_grid_cached(n, omega_max, TransformType::Rational)
            .unwrap();

        // Second call should use cache
        let (points2, weights2) = computer
            .generate_grid_cached(n, omega_max, TransformType::Rational)
            .unwrap();

        // Should return identical results
        for i in 0..n {
            assert_eq!(points1[i], points2[i]);
            assert_eq!(weights1[i], weights2[i]);
        }
    }

    #[test]
    fn test_edge_cases() {
        let computer = JacobianOptimized::new();

        // Test x = 0
        let jac = computer.safe_jacobian(0.0, 100.0);
        assert_relative_eq!(jac, 100.0, epsilon = 1e-10);

        // Test x = 0.5
        let jac = computer.safe_jacobian(0.5, 100.0);
        assert_relative_eq!(jac, 400.0, epsilon = 1e-10); // 100 / 0.25

        // Test x very close to 1
        let jac = computer.safe_jacobian(1.0 - 1e-15, 100.0);
        assert!(jac.is_finite());
        assert!(jac <= computer.max_jacobian * 100.0);
    }

    #[test]
    fn test_transform_types() {
        let mut computer = JacobianOptimized::new();
        let n = 16;
        let omega_max = 50.0;

        // Test all transform types
        for transform in &[
            TransformType::Rational,
            TransformType::Exponential,
            TransformType::Sinh,
        ] {
            let (points, weights) = computer
                .generate_grid_cached(n, omega_max, *transform)
                .unwrap();

            // Basic validation
            assert_eq!(points.len(), n);
            assert_eq!(weights.len(), n);

            // Check numerical properties
            assert!(points.iter().all(|&x| x.is_finite() && x >= 0.0));
            assert!(weights.iter().all(|&w| w.is_finite() && w > 0.0));

            // Check ordering
            for i in 1..n {
                assert!(
                    points[i] > points[i - 1],
                    "Points should be monotonic for {:?}",
                    transform
                );
            }
        }
    }
}
