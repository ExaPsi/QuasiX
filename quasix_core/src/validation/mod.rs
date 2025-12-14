//! Validation module for comparing CD and AC methods with statistical analysis
//!
//! This module provides a comprehensive framework for validating the accuracy
//! and performance of contour deformation (CD) versus analytical continuation (AC)
//! methods in GW calculations. Key features include:
//!
//! - Parallel execution of CD and AC calculations
//! - Statistical comparison with MAD < 0.05 eV criterion
//! - Bootstrap confidence intervals for uncertainty quantification
//! - Memory-efficient shared intermediate storage
//! - Performance profiling with minimal overhead
//!
//! # Architecture
//!
//! The validation framework uses a thread pool with 40/40/20 allocation:
//! - 40% threads for CD calculations
//! - 40% threads for AC calculations
//! - 20% threads for statistical analysis
//!
//! # Example
//!
//! ```rust,ignore
//! use quasix_core::validation::{ComparisonHarness, ValidationConfig};
//!
//! let config = ValidationConfig::default();
//! let harness = ComparisonHarness::new(config)?;
//!
//! let results = harness.run_comparison(
//!     &cd_data,
//!     &ac_data,
//!     &reference_data
//! )?;
//!
//! assert!(results.mad < 0.05); // MAD < 0.05 eV criterion
//! ```

pub mod algorithms;
pub mod comparison_harness;
pub mod gw100;
pub mod parameters;
pub mod statistics;

#[cfg(feature = "python")]
pub mod python_bindings;

// Re-export main types from algorithms module
pub use algorithms::{
    BootstrapCI, CacheEntry, CacheMetrics, CachePriority, CachingStrategy, CalculationStep,
    ConvergenceData, ConvergenceMonitor, ErrorEstimate, ErrorPropagator, MADResult, ModifiedZScore,
    OutlierClassification, OutlierClassifier, PerformanceProfiler, RobustMADCalculator,
    TimingStats, UnifiedParameterMapper, ZScoreResult,
};

// Re-export from other modules
pub use comparison_harness::{ComparisonHarness, ComparisonResult, MethodResult, ValidationConfig};
pub use gw100::{
    compute_validation_stats, BenchmarkExecutor, BenchmarkResult, BenchmarkTask, OptimizedArray,
    SimdStatistics, ValidationStats,
};
pub use parameters::{ACParameters, CDParameters, UnifiedParameters};
pub use statistics::{
    bootstrap_confidence_interval, compute_mad_simd, detect_outliers, StatisticalMetrics,
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mad_criterion() {
        // Test that our MAD calculation meets the 0.05 eV criterion
        let cd_values = vec![1.0, 1.02, 0.98, 1.01, 0.99];
        let ac_values = vec![1.01, 1.03, 0.97, 1.02, 0.98];

        let mad = compute_mad_simd(&cd_values, &ac_values).unwrap();
        assert!(mad < 0.05, "MAD should be less than 0.05 eV");
    }

    #[test]
    fn test_module_imports() {
        // Verify all submodules are accessible
        let _ = ValidationConfig::default();
        let _ = UnifiedParameters::default();
        let _ = StatisticalMetrics::default();
    }
}
