//! Advanced algorithms for CD/AC comparison and validation
//!
//! This module implements sophisticated algorithms for comparing contour deformation
//! and analytical continuation methods, including robust statistical measures,
//! performance profiling, and convergence monitoring.

use crate::common::{QuasixError, Result};
use parking_lot::RwLock;
use rand::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::sync::Arc;
use std::time::{Duration, Instant};

/// Unified parameter mapper for consistent configurations
#[derive(Debug, Clone)]
pub struct UnifiedParameterMapper {
    /// Parameter mapping table
    mappings: HashMap<String, f64>,
}

impl UnifiedParameterMapper {
    /// Create new mapper
    pub fn new() -> Self {
        Self {
            mappings: HashMap::new(),
        }
    }

    /// Add parameter mapping
    pub fn add_mapping(&mut self, key: String, value: f64) {
        self.mappings.insert(key, value);
    }

    /// Get mapped value
    pub fn get(&self, key: &str) -> Option<f64> {
        self.mappings.get(key).copied()
    }
}

impl Default for UnifiedParameterMapper {
    fn default() -> Self {
        Self::new()
    }
}

/// Result of MAD calculation with diagnostics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MADResult {
    /// Mean absolute deviation value
    pub mad: f64,
    /// Standard deviation of differences
    pub std_dev: f64,
    /// Number of data points
    pub n_points: usize,
    /// Confidence interval
    pub confidence_interval: (f64, f64),
}

/// Robust MAD calculator with outlier resistance
pub struct RobustMADCalculator {
    /// Trimming percentage for outliers
    trim_percent: f64,
}

impl RobustMADCalculator {
    /// Create new calculator
    pub fn new(trim_percent: f64) -> Self {
        Self { trim_percent }
    }

    /// Calculate robust MAD
    pub fn calculate(&self, cd_values: &[f64], ac_values: &[f64]) -> Result<MADResult> {
        if cd_values.len() != ac_values.len() {
            return Err(QuasixError::InvalidInput(
                "Array lengths must match".to_string(),
            ));
        }

        let mut differences: Vec<f64> = cd_values
            .iter()
            .zip(ac_values.iter())
            .map(|(cd, ac)| (cd - ac).abs())
            .collect();

        // Sort for trimming
        differences.sort_by(|a, b| a.partial_cmp(b).unwrap());

        // Trim outliers
        let trim_count = (differences.len() as f64 * self.trim_percent / 100.0) as usize;
        let trimmed = &differences[trim_count..differences.len() - trim_count];

        // Calculate MAD
        let mad = trimmed.iter().sum::<f64>() / trimmed.len() as f64;

        // Calculate standard deviation
        let mean = trimmed.iter().sum::<f64>() / trimmed.len() as f64;
        let variance =
            trimmed.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / trimmed.len() as f64;
        let std_dev = variance.sqrt();

        // Simple confidence interval
        let ci_margin = 1.96 * std_dev / (trimmed.len() as f64).sqrt();
        let confidence_interval = (mad - ci_margin, mad + ci_margin);

        Ok(MADResult {
            mad,
            std_dev,
            n_points: trimmed.len(),
            confidence_interval,
        })
    }
}

/// Bootstrap confidence interval calculator
pub struct BootstrapCI {
    /// Number of bootstrap samples
    n_samples: usize,
    /// Confidence level
    confidence_level: f64,
}

impl BootstrapCI {
    /// Create new bootstrap calculator
    pub fn new(n_samples: usize, confidence_level: f64) -> Self {
        Self {
            n_samples,
            confidence_level,
        }
    }

    /// Calculate confidence interval
    pub fn calculate(&self, data: &[f64]) -> (f64, f64) {
        use rand::rng;
        let mut rng = rng();
        let mut bootstrap_means = Vec::with_capacity(self.n_samples);

        for _ in 0..self.n_samples {
            let mut sample_sum = 0.0;
            for _ in 0..data.len() {
                let idx = rng.random_range(0..data.len());
                sample_sum += data[idx];
            }
            bootstrap_means.push(sample_sum / data.len() as f64);
        }

        bootstrap_means.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let alpha = 1.0 - self.confidence_level;
        let lower_idx = ((alpha / 2.0) * self.n_samples as f64) as usize;
        let upper_idx = ((1.0 - alpha / 2.0) * self.n_samples as f64) as usize;

        (
            bootstrap_means[lower_idx],
            bootstrap_means[upper_idx.min(self.n_samples - 1)],
        )
    }
}

/// Modified Z-score for outlier detection
pub struct ModifiedZScore {
    /// Threshold for outlier detection
    threshold: f64,
}

impl ModifiedZScore {
    /// Create new detector
    pub fn new(threshold: f64) -> Self {
        Self { threshold }
    }

    /// Detect outliers
    pub fn detect(&self, values: &[f64]) -> Vec<usize> {
        if values.is_empty() {
            return Vec::new();
        }

        // Calculate median
        let mut sorted = values.to_vec();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let median = if sorted.len() % 2 == 0 {
            (sorted[sorted.len() / 2 - 1] + sorted[sorted.len() / 2]) / 2.0
        } else {
            sorted[sorted.len() / 2]
        };

        // Calculate MAD
        let deviations: Vec<f64> = values.iter().map(|x| (x - median).abs()).collect();
        let mut sorted_dev = deviations.clone();
        sorted_dev.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mad = if sorted_dev.len() % 2 == 0 {
            (sorted_dev[sorted_dev.len() / 2 - 1] + sorted_dev[sorted_dev.len() / 2]) / 2.0
        } else {
            sorted_dev[sorted_dev.len() / 2]
        };

        // Detect outliers
        const CONSISTENCY_CONSTANT: f64 = 1.4826;
        let mad_scaled = CONSISTENCY_CONSTANT * mad;

        if mad_scaled < 1e-10 {
            return Vec::new();
        }

        values
            .iter()
            .enumerate()
            .filter_map(|(i, &x)| {
                let z_score = 0.6745 * (x - median) / mad_scaled;
                if z_score.abs() > self.threshold {
                    Some(i)
                } else {
                    None
                }
            })
            .collect()
    }
}

/// Z-score result
#[derive(Debug, Clone)]
pub struct ZScoreResult {
    /// Z-scores for each data point
    pub z_scores: Vec<f64>,
    /// Outlier indices
    pub outliers: Vec<usize>,
}

/// Performance profiler for method comparison
#[derive(Debug, Clone)]
pub struct PerformanceProfiler {
    /// Timing records
    timings: Arc<RwLock<HashMap<String, Duration>>>,
    /// Memory usage records
    memory: Arc<RwLock<HashMap<String, f64>>>,
}

impl PerformanceProfiler {
    /// Create new profiler
    pub fn new() -> Self {
        Self {
            timings: Arc::new(RwLock::new(HashMap::new())),
            memory: Arc::new(RwLock::new(HashMap::new())),
        }
    }

    /// Start timing a section
    pub fn start_timer(&self, name: &str) -> TimingGuard {
        TimingGuard {
            name: name.to_string(),
            start: Instant::now(),
            profiler: self.clone(),
        }
    }

    /// Record memory usage
    pub fn record_memory(&self, name: &str, memory_mb: f64) {
        self.memory.write().insert(name.to_string(), memory_mb);
    }

    /// Get timing statistics
    pub fn get_timings(&self) -> TimingStats {
        let timings = self.timings.read();
        let total: Duration = timings.values().sum();
        let mut breakdown = HashMap::new();

        for (name, duration) in timings.iter() {
            breakdown.insert(name.clone(), duration.as_millis() as u64);
        }

        TimingStats {
            total_ms: total.as_millis() as u64,
            breakdown,
        }
    }
}

impl Default for PerformanceProfiler {
    fn default() -> Self {
        Self::new()
    }
}

/// Timing guard for automatic duration recording
pub struct TimingGuard {
    name: String,
    start: Instant,
    profiler: PerformanceProfiler,
}

impl Drop for TimingGuard {
    fn drop(&mut self) {
        let duration = self.start.elapsed();
        self.profiler
            .timings
            .write()
            .insert(self.name.clone(), duration);
    }
}

/// Timing statistics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TimingStats {
    /// Total time in milliseconds
    pub total_ms: u64,
    /// Breakdown by section
    pub breakdown: HashMap<String, u64>,
}

/// Convergence monitor for iterative methods
#[derive(Debug, Clone)]
pub struct ConvergenceMonitor {
    /// Convergence history
    history: Arc<RwLock<Vec<ConvergenceData>>>,
    /// Tolerance for convergence
    tolerance: f64,
}

impl ConvergenceMonitor {
    /// Create new monitor
    pub fn new(tolerance: f64) -> Self {
        Self {
            history: Arc::new(RwLock::new(Vec::new())),
            tolerance,
        }
    }

    /// Add iteration data
    pub fn add_iteration(&self, iteration: usize, error: f64, value: f64) {
        self.history.write().push(ConvergenceData {
            iteration,
            error,
            value,
        });
    }

    /// Check if converged
    pub fn is_converged(&self) -> bool {
        let history = self.history.read();
        if let Some(last) = history.last() {
            last.error < self.tolerance
        } else {
            false
        }
    }

    /// Get convergence rate
    pub fn convergence_rate(&self) -> Option<f64> {
        let history = self.history.read();
        if history.len() < 3 {
            return None;
        }

        let n = history.len();
        let e_n = history[n - 1].error;
        let e_n_minus_1 = history[n - 2].error;
        let e_n_minus_2 = history[n - 3].error;

        if e_n_minus_1 > 1e-14 && e_n_minus_2 > 1e-14 {
            Some((e_n / e_n_minus_1).ln() / (e_n_minus_1 / e_n_minus_2).ln())
        } else {
            None
        }
    }
}

/// Convergence data point
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConvergenceData {
    /// Iteration number
    pub iteration: usize,
    /// Error estimate
    pub error: f64,
    /// Current value
    pub value: f64,
}

/// Caching strategy for intermediate results
#[derive(Debug, Clone)]
pub struct CachingStrategy {
    /// Cache storage
    cache: Arc<RwLock<HashMap<String, CacheEntry>>>,
    /// Maximum cache size (MB)
    max_size_mb: f64,
}

impl CachingStrategy {
    /// Create new caching strategy
    pub fn new(max_size_mb: f64) -> Self {
        Self {
            cache: Arc::new(RwLock::new(HashMap::new())),
            max_size_mb,
        }
    }

    /// Store in cache
    pub fn store(&self, key: String, data: Vec<f64>, priority: CachePriority) {
        let size_mb = (data.len() * 8) as f64 / 1_048_576.0;

        // Check if adding would exceed limit
        let current_size = self.total_cache_size();
        if current_size + size_mb > self.max_size_mb {
            self.evict_lru();
        }

        self.cache.write().insert(
            key,
            CacheEntry {
                data,
                priority,
                last_access: Instant::now(),
                size_mb,
            },
        );
    }

    /// Retrieve from cache
    pub fn get(&self, key: &str) -> Option<Vec<f64>> {
        let mut cache = self.cache.write();
        cache.get_mut(key).map(|entry| {
            entry.last_access = Instant::now();
            entry.data.clone()
        })
    }

    /// Get total cache size
    fn total_cache_size(&self) -> f64 {
        self.cache.read().values().map(|e| e.size_mb).sum()
    }

    /// Evict least recently used
    fn evict_lru(&self) {
        let mut cache = self.cache.write();
        if let Some(lru_key) = cache
            .iter()
            .filter(|(_, e)| e.priority != CachePriority::Critical)
            .min_by_key(|(_, e)| e.last_access)
            .map(|(k, _)| k.clone())
        {
            cache.remove(&lru_key);
        }
    }

    /// Get cache metrics
    pub fn get_metrics(&self) -> CacheMetrics {
        let cache = self.cache.read();
        CacheMetrics {
            n_entries: cache.len(),
            total_size_mb: self.total_cache_size(),
            max_size_mb: self.max_size_mb,
        }
    }
}

/// Cache entry
#[derive(Debug, Clone)]
pub struct CacheEntry {
    /// Cached data
    pub data: Vec<f64>,
    /// Priority level
    pub priority: CachePriority,
    /// Last access time
    pub last_access: Instant,
    /// Size in MB
    pub size_mb: f64,
}

/// Cache priority levels
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CachePriority {
    /// Critical - never evict
    Critical,
    /// High priority
    High,
    /// Normal priority
    Normal,
    /// Low priority
    Low,
}

/// Cache metrics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CacheMetrics {
    /// Number of cached entries
    pub n_entries: usize,
    /// Total cache size (MB)
    pub total_size_mb: f64,
    /// Maximum cache size (MB)
    pub max_size_mb: f64,
}

/// Outlier classifier
pub struct OutlierClassifier {
    /// Classification thresholds
    thresholds: Vec<f64>,
}

impl OutlierClassifier {
    /// Create new classifier
    pub fn new(thresholds: Vec<f64>) -> Self {
        Self { thresholds }
    }

    /// Classify outliers
    pub fn classify(&self, values: &[f64]) -> Vec<OutlierClassification> {
        values
            .iter()
            .map(|&v| {
                let abs_v = v.abs();
                if abs_v < self.thresholds[0] {
                    OutlierClassification::Normal
                } else if abs_v < self.thresholds[1] {
                    OutlierClassification::Mild
                } else if abs_v < self.thresholds[2] {
                    OutlierClassification::Moderate
                } else {
                    OutlierClassification::Severe
                }
            })
            .collect()
    }
}

/// Outlier classification levels
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum OutlierClassification {
    /// Normal value
    Normal,
    /// Mild outlier
    Mild,
    /// Moderate outlier
    Moderate,
    /// Severe outlier
    Severe,
}

/// Error propagator for uncertainty quantification
pub struct ErrorPropagator {
    /// Error estimates
    errors: HashMap<String, ErrorEstimate>,
}

impl ErrorPropagator {
    /// Create new propagator
    pub fn new() -> Self {
        Self {
            errors: HashMap::new(),
        }
    }

    /// Add error estimate
    pub fn add_error(&mut self, name: String, estimate: ErrorEstimate) {
        self.errors.insert(name, estimate);
    }

    /// Propagate errors through calculation
    pub fn propagate(&self, steps: &[CalculationStep]) -> ErrorEstimate {
        let mut total_variance = 0.0;

        for step in steps {
            if let Some(error) = self.errors.get(&step.name) {
                let contribution = error.std_dev.powi(2) * step.sensitivity.powi(2);
                total_variance += contribution;
            }
        }

        ErrorEstimate {
            mean: 0.0,
            std_dev: total_variance.sqrt(),
            confidence_level: 0.95,
        }
    }
}

impl Default for ErrorPropagator {
    fn default() -> Self {
        Self::new()
    }
}

/// Calculation step for error propagation
#[derive(Debug, Clone)]
pub struct CalculationStep {
    /// Step name
    pub name: String,
    /// Sensitivity coefficient
    pub sensitivity: f64,
}

/// Error estimate
#[derive(Debug, Clone)]
pub struct ErrorEstimate {
    /// Mean value
    pub mean: f64,
    /// Standard deviation
    pub std_dev: f64,
    /// Confidence level
    pub confidence_level: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_robust_mad_calculator() {
        let calculator = RobustMADCalculator::new(10.0);
        let cd = vec![1.0, 1.1, 0.9, 1.05, 0.95, 10.0]; // 10.0 is outlier
        let ac = vec![1.05, 1.15, 0.85, 1.0, 0.9, 10.5]; // 10.5 is outlier

        let result = calculator.calculate(&cd, &ac).unwrap();
        assert!(result.mad < 0.2); // Should be small after trimming outliers
    }

    #[test]
    fn test_modified_z_score() {
        let detector = ModifiedZScore::new(3.5);
        let values = vec![1.0, 1.1, 0.9, 1.05, 0.95, 10.0]; // 10.0 is outlier
        let outliers = detector.detect(&values);
        assert_eq!(outliers, vec![5]);
    }

    #[test]
    fn test_convergence_monitor() {
        let monitor = ConvergenceMonitor::new(1e-6);

        monitor.add_iteration(1, 0.1, 1.0);
        assert!(!monitor.is_converged());

        monitor.add_iteration(2, 0.01, 1.01);
        assert!(!monitor.is_converged());

        monitor.add_iteration(3, 1e-7, 1.0001);
        assert!(monitor.is_converged());
    }

    #[test]
    fn test_caching_strategy() {
        let cache = CachingStrategy::new(1.0); // 1 MB limit

        let data1 = vec![1.0; 10000];
        let data2 = vec![2.0; 10000];

        cache.store("key1".to_string(), data1.clone(), CachePriority::Normal);
        assert!(cache.get("key1").is_some());

        cache.store("key2".to_string(), data2.clone(), CachePriority::High);
        assert!(cache.get("key2").is_some());

        let metrics = cache.get_metrics();
        assert!(metrics.total_size_mb > 0.0);
    }
}
