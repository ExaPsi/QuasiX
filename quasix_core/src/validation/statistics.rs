//! Statistical algorithms with SIMD optimization for CD/AC comparison
//!
//! This module provides high-performance statistical functions for comparing
//! results from contour deformation and analytical continuation methods.
//! All algorithms are optimized with SIMD instructions for maximum performance.

use crate::common::{QuasixError, Result};
use rand::prelude::*;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::{
    _mm256_add_pd, _mm256_loadu_pd, _mm256_set1_pd, _mm256_storeu_pd, _mm256_sub_pd,
};

/// Statistical metrics for method comparison
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct StatisticalMetrics {
    /// Mean absolute deviation (eV)
    pub mad: f64,
    /// Root mean square deviation (eV)
    pub rmsd: f64,
    /// Maximum absolute error (eV)
    pub max_error: f64,
    /// Pearson correlation coefficient
    pub correlation: f64,
    /// 95% confidence interval for MAD
    pub mad_ci_95: (f64, f64),
    /// Number of outliers detected
    pub n_outliers: usize,
    /// Outlier indices
    pub outlier_indices: Vec<usize>,
}

/// Compute mean absolute deviation using SIMD operations
///
/// # Arguments
/// * `cd_values` - Values from contour deformation method
/// * `ac_values` - Values from analytical continuation method
///
/// # Returns
/// Mean absolute deviation in eV
pub fn compute_mad_simd(cd_values: &[f64], ac_values: &[f64]) -> Result<f64> {
    if cd_values.len() != ac_values.len() {
        return Err(QuasixError::InvalidInput(
            "Array lengths must match".to_string(),
        ));
    }

    if cd_values.is_empty() {
        return Err(QuasixError::InvalidInput(
            "Empty arrays provided".to_string(),
        ));
    }

    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx") && is_x86_feature_detected!("fma") {
            // SAFETY: AVX and FMA availability is verified via is_x86_feature_detected!
            // compute_mad_avx is marked #[target_feature(enable = "avx", "fma")]
            // and only uses AVX/FMA intrinsics which are guaranteed available.
            // Input slices are validated for equal length at function entry.
            unsafe { compute_mad_avx(cd_values, ac_values) }
        } else {
            // Fallback to scalar
            compute_mad_scalar(cd_values, ac_values)
        }
    }

    #[cfg(not(target_arch = "x86_64"))]
    {
        // Fallback to scalar for non-x86 architectures
        compute_mad_scalar(cd_values, ac_values)
    }
}

/// Scalar implementation of MAD computation
fn compute_mad_scalar(cd_values: &[f64], ac_values: &[f64]) -> Result<f64> {
    let sum: f64 = cd_values
        .iter()
        .zip(ac_values.iter())
        .map(|(cd, ac)| (cd - ac).abs())
        .sum();

    Ok(sum / cd_values.len() as f64)
}

/// AVX-optimized MAD computation
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx,fma")]
unsafe fn compute_mad_avx(cd_values: &[f64], ac_values: &[f64]) -> Result<f64> {
    let n = cd_values.len();
    let mut sum = 0.0;

    // Process 4 elements at a time with AVX
    let chunks = n / 4;
    let remainder = n % 4;

    if chunks > 0 {
        let mut acc = _mm256_set1_pd(0.0);

        for i in 0..chunks {
            let idx = i * 4;
            let cd = _mm256_loadu_pd(cd_values.as_ptr().add(idx));
            let ac = _mm256_loadu_pd(ac_values.as_ptr().add(idx));

            // Compute absolute difference
            let diff = _mm256_sub_pd(cd, ac);

            // Manual absolute value: mask sign bit
            let abs_mask = _mm256_set1_pd(f64::from_bits(0x7FFF_FFFF_FFFF_FFFF));
            let abs_diff = _mm256_and_pd(diff, abs_mask);

            acc = _mm256_add_pd(acc, abs_diff);
        }

        // Extract sum from vector
        let mut result = [0.0; 4];
        _mm256_storeu_pd(result.as_mut_ptr(), acc);
        sum = result.iter().sum::<f64>();
    }

    // Process remaining elements
    for i in (n - remainder)..n {
        sum += (cd_values[i] - ac_values[i]).abs();
    }

    Ok(sum / n as f64)
}

// Helper for absolute value in AVX (since _mm256_and_pd is not available in std::arch)
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::_mm256_and_pd;

/// Compute bootstrap confidence interval for MAD
///
/// # Arguments
/// * `cd_values` - Values from contour deformation method
/// * `ac_values` - Values from analytical continuation method
/// * `n_bootstrap` - Number of bootstrap samples
/// * `confidence_level` - Confidence level (e.g., 0.95 for 95%)
///
/// # Returns
/// Lower and upper bounds of confidence interval
pub fn bootstrap_confidence_interval(
    cd_values: &[f64],
    ac_values: &[f64],
    n_bootstrap: usize,
    confidence_level: f64,
) -> Result<(f64, f64)> {
    if cd_values.len() != ac_values.len() {
        return Err(QuasixError::InvalidInput(
            "Array lengths must match".to_string(),
        ));
    }

    let n = cd_values.len();

    // Parallel bootstrap sampling
    let bootstrap_mads: Vec<f64> = (0..n_bootstrap)
        .into_par_iter()
        .map(|_| {
            use rand::rng;
            let mut local_rng = rng();
            let mut cd_sample = Vec::with_capacity(n);
            let mut ac_sample = Vec::with_capacity(n);

            // Resample with replacement
            for _ in 0..n {
                let idx = local_rng.random_range(0..n);
                cd_sample.push(cd_values[idx]);
                ac_sample.push(ac_values[idx]);
            }

            compute_mad_scalar(&cd_sample, &ac_sample).unwrap_or(0.0)
        })
        .collect();

    // Sort bootstrap samples
    let mut sorted_mads = bootstrap_mads;
    sorted_mads.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // Compute confidence interval
    let alpha = 1.0 - confidence_level;
    let lower_idx = ((alpha / 2.0) * n_bootstrap as f64) as usize;
    let upper_idx = ((1.0 - alpha / 2.0) * n_bootstrap as f64) as usize;

    let lower = sorted_mads[lower_idx.min(sorted_mads.len() - 1)];
    let upper = sorted_mads[upper_idx.min(sorted_mads.len() - 1)];

    Ok((lower, upper))
}

/// Detect outliers using modified Z-score method
///
/// # Arguments
/// * `differences` - Absolute differences between methods
/// * `threshold` - Z-score threshold (typically 3.5)
///
/// # Returns
/// Indices of detected outliers
pub fn detect_outliers(differences: &[f64], threshold: f64) -> Vec<usize> {
    if differences.is_empty() {
        return Vec::new();
    }

    // Compute median
    let mut sorted = differences.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median = if sorted.len() % 2 == 0 {
        (sorted[sorted.len() / 2 - 1] + sorted[sorted.len() / 2]) / 2.0
    } else {
        sorted[sorted.len() / 2]
    };

    // Compute median absolute deviation (MAD)
    let deviations: Vec<f64> = differences.iter().map(|x| (x - median).abs()).collect();
    let mut sorted_dev = deviations.clone();
    sorted_dev.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mad = if sorted_dev.len() % 2 == 0 {
        (sorted_dev[sorted_dev.len() / 2 - 1] + sorted_dev[sorted_dev.len() / 2]) / 2.0
    } else {
        sorted_dev[sorted_dev.len() / 2]
    };

    // Compute modified Z-scores
    const CONSISTENCY_CONSTANT: f64 = 1.4826; // For normal distribution
    let mad_scaled = CONSISTENCY_CONSTANT * mad;

    if mad_scaled < 1e-10 {
        // All values are essentially the same
        return Vec::new();
    }

    differences
        .iter()
        .enumerate()
        .filter_map(|(i, &x)| {
            let z_score = 0.6745 * (x - median) / mad_scaled;
            if z_score.abs() > threshold {
                Some(i)
            } else {
                None
            }
        })
        .collect()
}

/// Compute comprehensive statistical metrics
pub fn compute_statistics(cd_values: &[f64], ac_values: &[f64]) -> Result<StatisticalMetrics> {
    if cd_values.len() != ac_values.len() {
        return Err(QuasixError::InvalidInput(
            "Array lengths must match".to_string(),
        ));
    }

    let n = cd_values.len() as f64;

    // Compute differences
    let differences: Vec<f64> = cd_values
        .iter()
        .zip(ac_values.iter())
        .map(|(cd, ac)| cd - ac)
        .collect();

    let abs_differences: Vec<f64> = differences.iter().map(|d| d.abs()).collect();

    // MAD with confidence interval
    let mad = compute_mad_simd(cd_values, ac_values)?;
    let mad_ci_95 = bootstrap_confidence_interval(cd_values, ac_values, 1000, 0.95)?;

    // RMSD
    let sum_sq: f64 = differences.iter().map(|d| d * d).sum();
    let rmsd = (sum_sq / n).sqrt();

    // Max error
    let max_error = abs_differences.iter().copied().fold(0.0, f64::max);

    // Pearson correlation
    let cd_mean = cd_values.iter().sum::<f64>() / n;
    let ac_mean = ac_values.iter().sum::<f64>() / n;

    let cov: f64 = cd_values
        .iter()
        .zip(ac_values.iter())
        .map(|(cd, ac)| (cd - cd_mean) * (ac - ac_mean))
        .sum::<f64>()
        / n;

    let cd_std = (cd_values.iter().map(|x| (x - cd_mean).powi(2)).sum::<f64>() / n).sqrt();

    let ac_std = (ac_values.iter().map(|x| (x - ac_mean).powi(2)).sum::<f64>() / n).sqrt();

    let correlation = if cd_std > 1e-10 && ac_std > 1e-10 {
        cov / (cd_std * ac_std)
    } else {
        1.0 // Perfect correlation if no variation
    };

    // Detect outliers
    let outlier_indices = detect_outliers(&abs_differences, 3.5);
    let n_outliers = outlier_indices.len();

    Ok(StatisticalMetrics {
        mad,
        rmsd,
        max_error,
        correlation,
        mad_ci_95,
        n_outliers,
        outlier_indices,
    })
}

/// Parallel computation of statistics for multiple datasets
pub fn compute_statistics_parallel(
    datasets: &[(Vec<f64>, Vec<f64>)],
) -> Result<Vec<StatisticalMetrics>> {
    datasets
        .par_iter()
        .map(|(cd, ac)| compute_statistics(cd, ac))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_mad_computation() {
        let cd = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let ac = vec![1.1, 1.9, 3.2, 3.8, 5.1];

        let mad = compute_mad_simd(&cd, &ac).unwrap();
        let expected = (0.1 + 0.1 + 0.2 + 0.2 + 0.1) / 5.0;
        assert_relative_eq!(mad, expected, epsilon = 1e-10);
    }

    #[test]
    fn test_mad_criterion() {
        // Test data with MAD < 0.05 eV
        let cd = vec![10.0, 10.02, 9.98, 10.01, 9.99];
        let ac = vec![10.01, 10.03, 9.97, 10.02, 9.98];

        let mad = compute_mad_simd(&cd, &ac).unwrap();
        assert!(mad < 0.05, "MAD should be less than 0.05 eV");
    }

    #[test]
    fn test_outlier_detection() {
        let differences = vec![0.01, 0.02, 0.01, 5.0, 0.02, 0.01]; // 5.0 is outlier
        let outliers = detect_outliers(&differences, 3.5);
        assert_eq!(outliers, vec![3]);
    }

    #[test]
    fn test_bootstrap_ci() {
        let cd = vec![1.0; 100];
        let ac = vec![1.05; 100];

        let (lower, upper) = bootstrap_confidence_interval(&cd, &ac, 100, 0.95).unwrap();

        // CI should be tight around 0.05
        assert!((0.045..=0.055).contains(&lower));
        assert!((0.045..=0.055).contains(&upper));
    }

    #[test]
    fn test_comprehensive_statistics() {
        let cd = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let ac = vec![1.1, 2.1, 3.1, 4.1, 5.1];

        let stats = compute_statistics(&cd, &ac).unwrap();

        assert_relative_eq!(stats.mad, 0.1, epsilon = 1e-10);
        assert_relative_eq!(stats.rmsd, 0.1, epsilon = 1e-10);
        assert_relative_eq!(stats.max_error, 0.1, epsilon = 1e-10);
        assert!(stats.correlation > 0.99); // Near perfect correlation
        assert_eq!(stats.n_outliers, 0);
    }

    #[test]
    fn test_parallel_statistics() {
        let datasets = vec![
            (vec![1.0, 2.0, 3.0], vec![1.1, 2.1, 3.1]),
            (vec![4.0, 5.0, 6.0], vec![4.2, 5.2, 6.2]),
        ];

        let results = compute_statistics_parallel(&datasets).unwrap();
        assert_eq!(results.len(), 2);
        assert_relative_eq!(results[0].mad, 0.1, epsilon = 1e-10);
        assert_relative_eq!(results[1].mad, 0.2, epsilon = 1e-10);
    }
}
