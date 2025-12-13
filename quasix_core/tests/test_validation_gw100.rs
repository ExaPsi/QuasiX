//! Integration tests for GW100 validation module

#[cfg(test)]
mod tests {
    use quasix_core::validation::{
        BenchmarkResult, SimdStatistics, ValidationStats, compute_validation_stats,
        BenchmarkExecutor, BenchmarkTask,
    };
    use approx::assert_relative_eq;

    #[test]
    fn test_benchmark_result_creation() {
        let result = BenchmarkResult::new(
            "H2O".to_string(),
            "def2-tzvp".to_string(),
            "def2-tzvp-jkfit".to_string(),
            vec![-15.0, -10.0, -5.0, 1.0, 2.0],
            vec![0.8, 0.85, 0.9, 0.92, 0.95],
            vec![-15.1, -10.05, -5.02, 1.05, 2.1],
            2, // HOMO index
            3, // LUMO index
            10.5,
        );

        // Check IP error (HOMO)
        assert_relative_eq!(result.ip_error, 0.02, epsilon = 1e-6);

        // Check EA error (LUMO)
        assert_relative_eq!(result.ea_error, 0.05, epsilon = 1e-6);

        // Check MAD
        let expected_mad = (0.1 + 0.05 + 0.02 + 0.05 + 0.1) / 5.0;
        assert_relative_eq!(result.mad(), expected_mad, epsilon = 1e-6);

        // Check threshold
        assert!(result.passes_threshold(0.1));
        assert!(!result.passes_threshold(0.01));
    }

    #[test]
    fn test_simd_mad_computation() {
        let calculated = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let reference = vec![1.1, 1.9, 3.2, 3.8, 5.1, 5.9, 7.2, 7.8];

        let mad = SimdStatistics::compute_mad_simd(&calculated, &reference).unwrap();

        // Expected MAD
        let expected = vec![0.1, 0.1, 0.2, 0.2, 0.1, 0.1, 0.2, 0.2]
            .iter()
            .sum::<f64>() / 8.0;

        assert_relative_eq!(mad, expected, epsilon = 1e-10);
    }

    #[test]
    fn test_correlation_computation() {
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let y = vec![1.05, 2.05, 3.05, 4.05, 5.05];

        let r2 = SimdStatistics::compute_correlation(&x, &y).unwrap();

        // Near perfect correlation
        assert!(r2 > 0.999);
    }

    #[test]
    fn test_rmse_computation() {
        let calculated = vec![1.0, 2.0, 3.0, 4.0];
        let reference = vec![1.1, 2.1, 3.1, 4.1];

        let rmse = SimdStatistics::compute_rmse(&calculated, &reference).unwrap();

        // RMSE should be 0.1
        assert_relative_eq!(rmse, 0.1, epsilon = 1e-10);
    }

    #[test]
    fn test_outlier_detection() {
        let values = vec![
            0.1, 0.12, 0.09, 0.11, 0.10, // Normal values around 0.1
            1.5,                          // Clear outlier
            0.08, 0.13, 0.095,           // More normal values
        ];

        let outliers = SimdStatistics::detect_outliers(&values, 3.5);

        // Should detect index 5 as outlier
        assert!(outliers.contains(&5));
        assert_eq!(outliers.len(), 1);
    }

    #[test]
    fn test_validation_stats_computation() {
        let results = vec![
            BenchmarkResult::new(
                "H2O".to_string(),
                "def2-tzvp".to_string(),
                "def2-tzvp-jkfit".to_string(),
                vec![-15.0, -10.0, -5.0],
                vec![0.8, 0.85, 0.9],
                vec![-15.1, -10.1, -5.05],
                1, // HOMO
                2, // LUMO
                10.0,
            ),
            BenchmarkResult::new(
                "NH3".to_string(),
                "def2-tzvp".to_string(),
                "def2-tzvp-jkfit".to_string(),
                vec![-12.0, -8.0, -3.0],
                vec![0.82, 0.88, 0.91],
                vec![-12.15, -8.1, -3.08],
                1, // HOMO
                2, // LUMO
                12.0,
            ),
        ];

        let stats = compute_validation_stats(&results).unwrap();

        // Check that stats are reasonable
        assert!(stats.mad > 0.0 && stats.mad < 0.2);
        assert!(stats.rmse > 0.0 && stats.rmse < 0.2);
        assert!(stats.correlation > 0.9);
        assert_eq!(stats.n_samples, 6);

        // Check validation passing
        assert!(stats.passes_validation(0.2, 0.9));
    }

    #[test]
    fn test_benchmark_executor() {
        let executor = BenchmarkExecutor::new(2);

        let tasks = vec![
            BenchmarkTask {
                molecule: "H2O".to_string(),
                basis_set: "def2-svp".to_string(),
                aux_basis: "def2-svp-jkfit".to_string(),
                conv_tol: 1e-5,
                max_iterations: 10,
                freq_integration: "cd".to_string(),
                n_frequencies: 16,
            },
            BenchmarkTask {
                molecule: "NH3".to_string(),
                basis_set: "def2-svp".to_string(),
                aux_basis: "def2-svp-jkfit".to_string(),
                conv_tol: 1e-5,
                max_iterations: 10,
                freq_integration: "cd".to_string(),
                n_frequencies: 16,
            },
        ];

        // Mock executor function
        let executor_fn = |task: &BenchmarkTask| -> quasix_core::common::Result<BenchmarkResult> {
            Ok(BenchmarkResult::new(
                task.molecule.clone(),
                task.basis_set.clone(),
                task.aux_basis.clone(),
                vec![-10.0, -5.0, 1.0],
                vec![0.9, 0.85, 0.92],
                vec![-10.1, -5.05, 1.08],
                1,
                2,
                1.0,
            ))
        };

        let results = executor.run_parallel(tasks, executor_fn);

        assert_eq!(results.len(), 2);
        assert_eq!(results[0].molecule, "H2O");
        assert_eq!(results[1].molecule, "NH3");
    }

    #[test]
    fn test_validation_summary() {
        let stats = ValidationStats {
            mad: 0.15,
            rmse: 0.18,
            max_error: 0.35,
            correlation: 0.98,
            n_samples: 100,
            mad_ci_lower: 0.12,
            mad_ci_upper: 0.18,
            n_outliers: 2,
            outlier_indices: vec![15, 72],
            per_molecule_mad: [
                ("H2O".to_string(), 0.12),
                ("NH3".to_string(), 0.18),
            ].iter().cloned().collect(),
        };

        let summary = stats.summary();
        assert!(summary.contains("MAD: 0.150"));
        assert!(summary.contains("R²: 0.980"));
        assert!(summary.contains("N: 100"));
        assert!(summary.contains("Outliers: 2"));
    }
}