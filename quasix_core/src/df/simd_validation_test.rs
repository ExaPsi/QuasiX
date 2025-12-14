//! Validation tests for SIMD optimizations
//!
//! Ensures SIMD implementations produce identical results to scalar versions

#[cfg(test)]
mod tests {
    use crate::df::mo_transform::{
        generate_mock_mo_coefficients, transform_mo_3center, TransformConfig,
    };
    use crate::df::simd_ops::{check_simd_support, simd_matmul_transpose};
    use ndarray::{Array2, Array3};

    /// Generate small test data
    fn generate_test_data(
        n_ao: usize,
        n_aux: usize,
        n_occ: usize,
        n_vir: usize,
    ) -> (Array3<f64>, Array2<f64>, Array2<f64>) {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(12345);

        // Generate symmetric 3-center integrals
        let mut j3c = Array3::<f64>::zeros((n_ao, n_ao, n_aux));
        for i in 0..n_ao {
            for j in 0..=i {
                for p in 0..n_aux {
                    let val = rng.random::<f64>() * 0.1;
                    j3c[[i, j, p]] = val;
                    j3c[[j, i, p]] = val;
                }
            }
        }

        let c_occ = generate_mock_mo_coefficients(n_ao, n_occ, 100);
        let c_vir = generate_mock_mo_coefficients(n_ao, n_vir, 101);

        (j3c, c_occ, c_vir)
    }

    #[test]
    fn test_simd_matmul_vs_ndarray() {
        let simd_support = check_simd_support();
        println!("SIMD Support: {:?}", simd_support);

        // Test various matrix sizes
        for &size in &[8, 16, 32, 64] {
            let m = size;
            let k = size;
            let n = size;

            let a = Array2::from_shape_fn((m, k), |(i, j)| (i + j) as f64 * 0.01);
            let b = Array2::from_shape_fn((k, n), |(i, j)| (i * j) as f64 * 0.01);

            // SIMD version
            let mut c_simd = Array2::zeros((m, n));
            simd_matmul_transpose(a.view(), b.view(), c_simd.view_mut());

            // Reference: ndarray dot
            let c_ref = a.dot(&b);

            // Compare
            let max_error = c_simd
                .iter()
                .zip(c_ref.iter())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max);

            assert!(
                max_error < 1e-14,
                "SIMD matmul error too large for size {}: {:.2e}",
                size,
                max_error
            );

            if simd_support.has_high_perf_simd() {
                println!(
                    "✓ Size {} passed with AVX2/FMA (max error: {:.2e})",
                    size, max_error
                );
            } else {
                println!(
                    "✓ Size {} passed with scalar fallback (max error: {:.2e})",
                    size, max_error
                );
            }
        }
    }

    #[test]
    fn test_mo_transform_simd_correctness() {
        // Test on realistic GW100 molecule sizes
        let test_cases = vec![
            ("H2O", 7, 76, 5, 2),
            ("NH3", 8, 90, 5, 3),
            ("CH4", 9, 104, 5, 4),
            ("CO", 10, 112, 7, 3),
        ];

        for (name, n_ao, n_aux, n_occ, n_vir) in test_cases {
            let (j3c_ao, c_occ, c_vir) = generate_test_data(n_ao, n_aux, n_occ, n_vir);

            // Transform with SIMD (when available)
            let result_simd = transform_mo_3center(&j3c_ao, &c_occ, &c_vir).unwrap();

            // Transform with explicit non-SIMD config (small system forces scalar path)
            let config = TransformConfig::small_system();
            let result_scalar = crate::df::mo_transform::transform_mo_3center_with_config(
                &j3c_ao, &c_occ, &c_vir, &config,
            )
            .unwrap();

            // Compare results
            let max_error = result_simd
                .iter()
                .zip(result_scalar.iter())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max);

            assert!(
                max_error < 1e-12,
                "{}: MO transform SIMD/scalar mismatch: {:.2e}",
                name,
                max_error
            );

            println!(
                "✓ {} MO transform validated (max error: {:.2e})",
                name, max_error
            );
        }
    }

    #[test]
    fn test_lock_free_parallel_correctness() {
        // Test that lock-free parallel implementation matches sequential
        let (j3c_ao, c_occ, c_vir) = generate_test_data(10, 112, 7, 3); // N2 size

        // Sequential reference
        let result_seq = transform_mo_3center(&j3c_ao, &c_occ, &c_vir).unwrap();

        // Parallel with lock-free
        use crate::df::parallel::{
            transform_mo_3center_parallel_optimized, BlockingStrategy, ThreadPoolConfig,
        };
        let result_par = transform_mo_3center_parallel_optimized(
            &j3c_ao,
            &c_occ,
            &c_vir,
            Some(ThreadPoolConfig::default()),
            Some(BlockingStrategy::default()),
        )
        .unwrap();

        // Compare
        let max_error = result_seq
            .iter()
            .zip(result_par.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);

        assert!(
            max_error < 1e-12,
            "Lock-free parallel result differs from sequential: {:.2e}",
            max_error
        );

        println!(
            "✓ Lock-free parallel transformation validated (max error: {:.2e})",
            max_error
        );
    }

    #[test]
    fn test_thread_coordination_no_deadlock() {
        // Ensure thread coordination doesn't cause deadlocks
        use crate::df::parallel::{
            transform_mo_3center_parallel_optimized, BlockingStrategy, ThreadPoolConfig,
        };

        let (j3c_ao, c_occ, c_vir) = generate_test_data(12, 132, 8, 4); // H2CO size

        // Test different configurations sequentially
        for &blas_threads in &[1, 2] {
            let mut config = ThreadPoolConfig::default();
            config.blas_threads = blas_threads;

            let result = transform_mo_3center_parallel_optimized(
                &j3c_ao,
                &c_occ,
                &c_vir,
                Some(config),
                Some(BlockingStrategy::default()),
            );

            assert!(
                result.is_ok(),
                "Thread config with blas_threads={} failed",
                blas_threads
            );
            println!(
                "✓ Thread coordination with blas_threads={} succeeded",
                blas_threads
            );
        }
    }

    #[test]
    fn test_simd_alignment() {
        // Test that SIMD handles non-aligned sizes correctly
        for &size in &[7, 13, 17, 23] {
            // Non-multiples of 4 (AVX2 width)
            let m = size;
            let k = size;
            let n = size;

            let a = Array2::from_shape_fn((m, k), |(i, j)| (i + j) as f64 * 0.01);
            let b = Array2::from_shape_fn((k, n), |(i, j)| (i * j) as f64 * 0.01);

            let mut c_simd = Array2::zeros((m, n));
            simd_matmul_transpose(a.view(), b.view(), c_simd.view_mut());

            let c_ref = a.dot(&b);

            let max_error = c_simd
                .iter()
                .zip(c_ref.iter())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max);

            assert!(
                max_error < 1e-14,
                "SIMD non-aligned size {} failed: {:.2e}",
                size,
                max_error
            );

            println!("✓ Non-aligned size {} handled correctly", size);
        }
    }
}
