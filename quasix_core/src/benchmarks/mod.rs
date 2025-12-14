//! High-performance GW100 benchmark validation framework
//!
//! This module provides efficient data structures and algorithms for
//! systematic validation against the GW100 reference database.
//!
//! ## Core Components
//!
//! - [`BenchmarkRunner`]: Parallel execution engine for molecule sets
//! - [`StatisticalAnalyzer`]: Numerically stable statistical algorithms
//! - [`ValidationPipeline`]: Automated validation and reporting
//! - [`ReferenceDatabase`]: Efficient reference data management
//! - [`parallel`]: High-performance parallel execution infrastructure

pub mod parallel;

use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::common::{QuasixError, Result};
// NOTE: QualityMetrics removed during G₀W₀ cleanup - will be re-implemented in S3-6
// use crate::gw::QualityMetrics;
use crate::io::{GWResults, MolecularData};

/// Placeholder for QualityMetrics (to be re-implemented in S3-6)
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct QualityMetrics {
    /// Convergence information (placeholder)
    pub convergence: ConvergenceInfo,
}

/// Placeholder for convergence information
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ConvergenceInfo {
    /// Oscillation metric (placeholder)
    pub oscillation: f64,
}

/// Benchmark configuration parameters
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkConfig {
    /// Maximum parallel molecules to compute
    pub max_parallel: usize,
    /// Enable result caching
    pub enable_cache: bool,
    /// Cache directory path
    pub cache_dir: String,
    /// Statistical confidence level (e.g., 0.95 for 95%)
    pub confidence_level: f64,
    /// Maximum allowed deviation in eV
    pub max_deviation_ev: f64,
    /// Outlier detection threshold (standard deviations)
    pub outlier_threshold: f64,
    /// Enable incremental computation
    pub incremental: bool,
    /// Enable work-stealing scheduler
    pub enable_work_stealing: bool,
    /// Memory limit in GB (optional)
    pub memory_limit_gb: Option<f64>,
    /// Batch strategy: "balanced", "memory", "time"
    pub batch_strategy: String,
}

impl Default for BenchmarkConfig {
    fn default() -> Self {
        Self {
            max_parallel: 4,
            enable_cache: true,
            cache_dir: ".quasix_benchmark_cache".to_string(),
            confidence_level: 0.95,
            max_deviation_ev: 0.2,
            outlier_threshold: 3.0,
            incremental: true,
            enable_work_stealing: true,
            memory_limit_gb: None,
            batch_strategy: "balanced".to_string(),
        }
    }
}

/// Individual molecule benchmark result
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MoleculeResult {
    /// Molecule identifier
    pub molecule: String,
    /// Calculated ionization potential (eV)
    pub ip_calc: f64,
    /// Reference ionization potential (eV)
    pub ip_ref: f64,
    /// Calculated electron affinity (eV)
    pub ea_calc: Option<f64>,
    /// Reference electron affinity (eV)
    pub ea_ref: Option<f64>,
    /// Quasiparticle energies (Ha)
    pub qp_energies: Vec<f64>,
    /// Z-factors (renormalization)
    pub z_factors: Vec<f64>,
    /// Computation time
    pub elapsed: Duration,
    /// Convergence metrics
    pub quality: QualityMetrics,
    /// Individual orbital deviations
    pub orbital_deviations: Vec<f64>,
}

impl MoleculeResult {
    /// Compute IP deviation in eV
    pub fn ip_deviation(&self) -> f64 {
        self.ip_calc - self.ip_ref
    }

    /// Compute EA deviation in eV if available
    pub fn ea_deviation(&self) -> Option<f64> {
        match (self.ea_calc, self.ea_ref) {
            (Some(calc), Some(ref_val)) => Some(calc - ref_val),
            _ => None,
        }
    }

    /// Check if result passes accuracy threshold
    pub fn passes_threshold(&self, max_dev: f64) -> bool {
        self.ip_deviation().abs() <= max_dev
            && self
                .ea_deviation()
                .is_none_or(|dev| dev.abs() <= max_dev * 1.5)
    }
}

/// Aggregated benchmark statistics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkStatistics {
    /// Number of molecules computed
    pub n_molecules: usize,
    /// Mean absolute deviation (eV)
    pub mad: f64,
    /// Root mean square deviation (eV)
    pub rmsd: f64,
    /// Maximum absolute deviation (eV)
    pub max_dev: f64,
    /// Standard deviation (eV)
    pub std_dev: f64,
    /// Pearson correlation coefficient
    pub correlation: f64,
    /// R-squared value
    pub r_squared: f64,
    /// Linear regression slope
    pub slope: f64,
    /// Linear regression intercept (eV)
    pub intercept: f64,
    /// 95% confidence interval (eV)
    pub confidence_interval: (f64, f64),
    /// Number of outliers detected
    pub n_outliers: usize,
    /// Pass rate (fraction)
    pub pass_rate: f64,
}

/// Statistical analyzer with numerically stable algorithms
pub struct StatisticalAnalyzer {
    config: BenchmarkConfig,
}

impl StatisticalAnalyzer {
    /// Create new analyzer with configuration
    pub fn new(config: BenchmarkConfig) -> Self {
        Self { config }
    }

    /// Compute comprehensive statistics from results
    pub fn analyze(&self, results: &[MoleculeResult]) -> BenchmarkStatistics {
        let n = results.len();
        if n == 0 {
            return self.empty_statistics();
        }

        // Extract deviations
        let deviations: Vec<f64> = results.iter().map(|r| r.ip_deviation()).collect();

        // Compute basic statistics with numerical stability
        let mad = self.compute_mad(&deviations);
        let rmsd = self.compute_rmsd(&deviations);
        let (_mean, std_dev) = self.compute_mean_std(&deviations);
        let max_dev = deviations.iter().map(|d| d.abs()).fold(0.0, f64::max);

        // Linear regression with numerical stability
        let (slope, intercept, correlation, r_squared) = self.compute_regression(results);

        // Confidence interval using t-distribution approximation
        let confidence_interval = self.compute_confidence_interval(&deviations);

        // Outlier detection using modified Z-score
        let n_outliers = self.detect_outliers(&deviations);

        // Pass rate calculation
        let n_passed = results
            .iter()
            .filter(|r| r.passes_threshold(self.config.max_deviation_ev))
            .count();
        let pass_rate = n_passed as f64 / n as f64;

        BenchmarkStatistics {
            n_molecules: n,
            mad,
            rmsd,
            max_dev,
            std_dev,
            correlation,
            r_squared,
            slope,
            intercept,
            confidence_interval,
            n_outliers,
            pass_rate,
        }
    }

    /// Compute Mean Absolute Deviation with numerical stability
    fn compute_mad(&self, deviations: &[f64]) -> f64 {
        if deviations.is_empty() {
            return 0.0;
        }

        // Use Kahan summation for numerical stability
        let mut sum = 0.0;
        let mut c = 0.0;

        for &dev in deviations {
            let y = dev.abs() - c;
            let t = sum + y;
            c = (t - sum) - y;
            sum = t;
        }

        sum / deviations.len() as f64
    }

    /// Compute Root Mean Square Deviation
    fn compute_rmsd(&self, deviations: &[f64]) -> f64 {
        if deviations.is_empty() {
            return 0.0;
        }

        // Numerically stable computation
        let mut sum_sq = 0.0;
        let mut c = 0.0;

        for &dev in deviations {
            let sq = dev * dev;
            let y = sq - c;
            let t = sum_sq + y;
            c = (t - sum_sq) - y;
            sum_sq = t;
        }

        (sum_sq / deviations.len() as f64).sqrt()
    }

    /// Compute mean and standard deviation with Welford's algorithm
    fn compute_mean_std(&self, values: &[f64]) -> (f64, f64) {
        if values.is_empty() {
            return (0.0, 0.0);
        }
        if values.len() == 1 {
            return (values[0], 0.0);
        }

        // Welford's online algorithm for numerical stability
        let mut mean = 0.0;
        let mut m2 = 0.0;
        let mut n = 0.0;

        for &value in values {
            n += 1.0;
            let delta = value - mean;
            mean += delta / n;
            let delta2 = value - mean;
            m2 += delta * delta2;
        }

        let variance = m2 / (n - 1.0);
        (mean, variance.sqrt())
    }

    /// Compute linear regression with numerical stability
    fn compute_regression(&self, results: &[MoleculeResult]) -> (f64, f64, f64, f64) {
        let n = results.len() as f64;
        if n < 2.0 {
            return (1.0, 0.0, 1.0, 1.0);
        }

        // Extract x (reference) and y (calculated) values
        let x: Vec<f64> = results.iter().map(|r| r.ip_ref).collect();
        let y: Vec<f64> = results.iter().map(|r| r.ip_calc).collect();

        // Compute means
        let x_mean = x.iter().sum::<f64>() / n;
        let y_mean = y.iter().sum::<f64>() / n;

        // Compute covariance and variance with numerical stability
        let mut cov = 0.0;
        let mut var_x = 0.0;
        let mut var_y = 0.0;

        for i in 0..results.len() {
            let dx = x[i] - x_mean;
            let dy = y[i] - y_mean;
            cov += dx * dy;
            var_x += dx * dx;
            var_y += dy * dy;
        }

        cov /= n - 1.0;
        var_x /= n - 1.0;
        var_y /= n - 1.0;

        // Compute regression parameters
        let slope = if var_x > 1e-10 { cov / var_x } else { 1.0 };
        let intercept = y_mean - slope * x_mean;

        // Correlation and R-squared
        let correlation = if var_x > 1e-10 && var_y > 1e-10 {
            cov / (var_x.sqrt() * var_y.sqrt())
        } else {
            1.0
        };
        let r_squared = correlation * correlation;

        (slope, intercept, correlation, r_squared)
    }

    /// Compute confidence interval using t-distribution
    fn compute_confidence_interval(&self, deviations: &[f64]) -> (f64, f64) {
        let (mean, std_dev) = self.compute_mean_std(deviations);
        let n = deviations.len() as f64;

        if n < 2.0 {
            return (mean, mean);
        }

        // T-distribution critical value approximation
        // For 95% confidence and large n, approaches 1.96
        let t_critical = if n > 30.0 {
            1.96
        } else {
            // Approximation for smaller samples
            2.0 + 4.0 / n
        };

        let margin = t_critical * std_dev / n.sqrt();
        (mean - margin, mean + margin)
    }

    /// Detect outliers using modified Z-score method
    fn detect_outliers(&self, values: &[f64]) -> usize {
        if values.len() < 4 {
            return 0;
        }

        // Compute median absolute deviation (MAD)
        let mut sorted = values.to_vec();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let median = sorted[sorted.len() / 2];

        let mad = sorted.iter().map(|v| (v - median).abs()).sum::<f64>() / sorted.len() as f64;

        if mad < 1e-10 {
            return 0;
        }

        // Modified Z-score
        let threshold = self.config.outlier_threshold;
        values
            .iter()
            .filter(|&&v| {
                let modified_z = 0.6745 * (v - median) / mad;
                modified_z.abs() > threshold
            })
            .count()
    }

    /// Create empty statistics structure
    fn empty_statistics(&self) -> BenchmarkStatistics {
        BenchmarkStatistics {
            n_molecules: 0,
            mad: 0.0,
            rmsd: 0.0,
            max_dev: 0.0,
            std_dev: 0.0,
            correlation: 1.0,
            r_squared: 1.0,
            slope: 1.0,
            intercept: 0.0,
            confidence_interval: (0.0, 0.0),
            n_outliers: 0,
            pass_rate: 0.0,
        }
    }
}

/// Parallel benchmark execution engine
pub struct BenchmarkRunner {
    config: BenchmarkConfig,
    results: Arc<Mutex<Vec<MoleculeResult>>>,
    cache: Option<BenchmarkCache>,
}

impl BenchmarkRunner {
    /// Create new benchmark runner
    pub fn new(config: BenchmarkConfig) -> Self {
        let cache = if config.enable_cache {
            Some(BenchmarkCache::new(&config.cache_dir))
        } else {
            None
        };

        Self {
            config,
            results: Arc::new(Mutex::new(Vec::new())),
            cache,
        }
    }

    /// Execute benchmarks for molecule set in parallel
    pub fn run_benchmarks(
        &mut self,
        molecules: Vec<MolecularData>,
        references: &ReferenceDatabase,
    ) -> Result<Vec<MoleculeResult>> {
        // Set up parallel execution pool
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.config.max_parallel)
            .build()
            .map_err(|e| QuasixError::NumericalError(e.to_string()))?;

        let results_arc = Arc::clone(&self.results);
        let cache = self.cache.clone();
        let config = self.config.clone();

        let results = pool.install(|| {
            molecules
                .par_iter()
                .map(|mol_data| {
                    // Check cache first
                    if let Some(ref cache) = cache {
                        if let Some(cached) = cache.get(&mol_data.name) {
                            if config.incremental {
                                return Ok(cached);
                            }
                        }
                    }

                    // Run GW calculation
                    let start = Instant::now();
                    let gw_result = self.run_single_gw(mol_data)?;
                    let elapsed = start.elapsed();

                    // Get reference values
                    let ref_data = references.get(&mol_data.name).ok_or_else(|| {
                        QuasixError::ValidationError(format!("No reference for {}", mol_data.name))
                    })?;

                    // Extract IP/EA
                    let ip_calc = self.extract_ip(&gw_result);
                    let ea_calc = self.extract_ea(&gw_result);

                    // Build result
                    let result = MoleculeResult {
                        molecule: mol_data.name.clone(),
                        ip_calc,
                        ip_ref: ref_data.ip,
                        ea_calc,
                        ea_ref: ref_data.ea,
                        qp_energies: gw_result.qp_energies.to_vec(),
                        z_factors: gw_result.z_factors.to_vec(),
                        elapsed,
                        quality: gw_result.quality_metrics.clone(),
                        orbital_deviations: self.compute_orbital_deviations(&gw_result, ref_data),
                    };

                    // Cache result
                    if let Some(ref cache) = cache {
                        cache.store(&mol_data.name, &result);
                    }

                    // Store in results
                    results_arc.lock().unwrap().push(result.clone());

                    Ok(result)
                })
                .collect::<Result<Vec<_>>>()
        })?;

        Ok(results)
    }

    /// Run single GW calculation
    fn run_single_gw(&self, _mol_data: &MolecularData) -> Result<GWResults> {
        // This would call the actual GW implementation
        // For now, placeholder
        todo!("Integrate with actual GW driver")
    }

    /// Extract ionization potential (negative HOMO energy)
    fn extract_ip(&self, gw_result: &GWResults) -> f64 {
        let n_occ = gw_result.n_occ;
        if n_occ > 0 {
            -gw_result.qp_energies[n_occ - 1] * 27.211_386 // Ha to eV
        } else {
            0.0
        }
    }

    /// Extract electron affinity (negative LUMO energy)
    fn extract_ea(&self, gw_result: &GWResults) -> Option<f64> {
        let n_occ = gw_result.n_occ;
        if n_occ < gw_result.qp_energies.len() {
            Some(-gw_result.qp_energies[n_occ] * 27.211_386) // Ha to eV
        } else {
            None
        }
    }

    /// Compute per-orbital deviations
    fn compute_orbital_deviations(
        &self,
        gw_result: &GWResults,
        ref_data: &ReferenceData,
    ) -> Vec<f64> {
        // Match orbital energies and compute deviations
        gw_result
            .qp_energies
            .iter()
            .zip(ref_data.orbital_energies.iter())
            .map(|(&calc, &ref_val)| (calc - ref_val) * 27.211_386)
            .collect()
    }

    /// Get collected results
    pub fn get_results(&self) -> Vec<MoleculeResult> {
        self.results.lock().unwrap().clone()
    }
}

/// Reference data for a single molecule
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReferenceData {
    /// Molecule name
    pub molecule: String,
    /// Reference ionization potential (eV)
    pub ip: f64,
    /// Reference electron affinity (eV)
    pub ea: Option<f64>,
    /// Reference orbital energies (Ha)
    pub orbital_energies: Vec<f64>,
    /// Method used for reference
    pub method: String,
    /// Basis set used
    pub basis: String,
}

/// Reference database manager
#[derive(Debug, Clone)]
pub struct ReferenceDatabase {
    data: HashMap<String, ReferenceData>,
}

impl ReferenceDatabase {
    /// Load GW100 reference data
    pub fn load_gw100() -> Result<Self> {
        // Load from embedded data or file
        let mut data = HashMap::new();

        // H2O reference
        data.insert(
            "H2O".to_string(),
            ReferenceData {
                molecule: "H2O".to_string(),
                ip: 12.62, // eV
                ea: None,
                orbital_energies: vec![-0.4634, -0.3903, -0.2569, -0.1822],
                method: "evGW@PBE".to_string(),
                basis: "cc-pVDZ".to_string(),
            },
        );

        // NH3 reference
        data.insert(
            "NH3".to_string(),
            ReferenceData {
                molecule: "NH3".to_string(),
                ip: 10.82, // eV
                ea: None,
                orbital_energies: vec![-0.3974, -0.3149, -0.2188],
                method: "evGW@PBE".to_string(),
                basis: "cc-pVDZ".to_string(),
            },
        );

        // CO reference
        data.insert(
            "CO".to_string(),
            ReferenceData {
                molecule: "CO".to_string(),
                ip: 14.01, // eV
                ea: Some(1.33),
                orbital_energies: vec![-0.5149, -0.4321, -0.3894],
                method: "evGW@PBE".to_string(),
                basis: "cc-pVDZ".to_string(),
            },
        );

        // Benzene reference
        data.insert(
            "benzene".to_string(),
            ReferenceData {
                molecule: "benzene".to_string(),
                ip: 9.24, // eV
                ea: Some(-1.12),
                orbital_energies: vec![-0.3397, -0.2867, -0.2422],
                method: "evGW@PBE".to_string(),
                basis: "cc-pVDZ".to_string(),
            },
        );

        Ok(Self { data })
    }

    /// Get reference data for molecule
    pub fn get(&self, molecule: &str) -> Option<&ReferenceData> {
        self.data.get(molecule)
    }

    /// Add custom reference data
    pub fn add(&mut self, ref_data: ReferenceData) {
        self.data.insert(ref_data.molecule.clone(), ref_data);
    }
}

/// Benchmark result cache for incremental computation
#[derive(Debug, Clone)]
struct BenchmarkCache {
    cache_dir: String,
    data: Arc<Mutex<HashMap<String, MoleculeResult>>>,
}

impl BenchmarkCache {
    /// Create new cache
    fn new(cache_dir: &str) -> Self {
        // Create cache directory if needed
        std::fs::create_dir_all(cache_dir).ok();

        // Load existing cache
        let data = Self::load_cache(cache_dir).unwrap_or_default();

        Self {
            cache_dir: cache_dir.to_string(),
            data: Arc::new(Mutex::new(data)),
        }
    }

    /// Get cached result
    fn get(&self, molecule: &str) -> Option<MoleculeResult> {
        self.data.lock().unwrap().get(molecule).cloned()
    }

    /// Store result in cache
    fn store(&self, molecule: &str, result: &MoleculeResult) {
        self.data
            .lock()
            .unwrap()
            .insert(molecule.to_string(), result.clone());
        self.save_cache();
    }

    /// Load cache from disk
    fn load_cache(cache_dir: &str) -> Result<HashMap<String, MoleculeResult>> {
        let path = format!("{}/benchmark_cache.json", cache_dir);
        if !std::path::Path::new(&path).exists() {
            return Ok(HashMap::new());
        }

        let data =
            std::fs::read_to_string(&path).map_err(|e| QuasixError::IoError(e.to_string()))?;

        serde_json::from_str(&data).map_err(|e| QuasixError::SerializationError(e.to_string()))
    }

    /// Save cache to disk
    fn save_cache(&self) {
        let data = self.data.lock().unwrap();
        let json = serde_json::to_string_pretty(&*data).ok();

        if let Some(json) = json {
            let path = format!("{}/benchmark_cache.json", self.cache_dir);
            std::fs::write(&path, json).ok();
        }
    }
}

/// Validation pipeline for automated pass/fail determination
pub struct ValidationPipeline {
    config: BenchmarkConfig,
    analyzer: StatisticalAnalyzer,
}

impl ValidationPipeline {
    /// Create new validation pipeline
    pub fn new(config: BenchmarkConfig) -> Self {
        let analyzer = StatisticalAnalyzer::new(config.clone());
        Self { config, analyzer }
    }

    /// Validate results against thresholds
    pub fn validate(&self, results: &[MoleculeResult]) -> ValidationReport {
        let stats = self.analyzer.analyze(results);

        // Check individual pass criteria
        let mad_pass = stats.mad <= self.config.max_deviation_ev;
        let rmsd_pass = stats.rmsd <= self.config.max_deviation_ev * 1.2;
        let correlation_pass = stats.r_squared >= 0.95;
        let pass_rate_pass = stats.pass_rate >= 0.75;

        let overall_pass = mad_pass && rmsd_pass && correlation_pass && pass_rate_pass;

        // Identify problem molecules
        let failed_molecules: Vec<String> = results
            .iter()
            .filter(|r| !r.passes_threshold(self.config.max_deviation_ev))
            .map(|r| r.molecule.clone())
            .collect();

        // Generate diagnostics
        let diagnostics = self.generate_diagnostics(results, &stats);

        ValidationReport {
            statistics: stats,
            mad_pass,
            rmsd_pass,
            correlation_pass,
            pass_rate_pass,
            overall_pass,
            failed_molecules,
            diagnostics,
            timestamp: chrono::Utc::now(),
        }
    }

    /// Generate detailed diagnostics
    fn generate_diagnostics(
        &self,
        results: &[MoleculeResult],
        stats: &BenchmarkStatistics,
    ) -> Vec<String> {
        let mut diagnostics = Vec::new();

        // Overall performance
        diagnostics.push(format!(
            "Benchmark completed for {} molecules",
            stats.n_molecules
        ));
        diagnostics.push(format!("Mean Absolute Deviation: {:.3} eV", stats.mad));
        diagnostics.push(format!("Pass rate: {:.1}%", stats.pass_rate * 100.0));

        // Outlier analysis
        if stats.n_outliers > 0 {
            diagnostics.push(format!("Warning: {} outliers detected", stats.n_outliers));
        }

        // Systematic bias check
        if (stats.intercept).abs() > 0.1 {
            diagnostics.push(format!(
                "Possible systematic bias: intercept = {:.3} eV",
                stats.intercept
            ));
        }

        // Convergence issues
        let convergence_issues: Vec<_> = results
            .iter()
            .filter(|r| r.quality.convergence.oscillation > 0.5)
            .map(|r| r.molecule.clone())
            .collect();

        if !convergence_issues.is_empty() {
            diagnostics.push(format!(
                "Convergence warnings for: {}",
                convergence_issues.join(", ")
            ));
        }

        diagnostics
    }
}

/// Validation report with pass/fail determination
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValidationReport {
    /// Statistical analysis results
    pub statistics: BenchmarkStatistics,
    /// MAD criterion passed
    pub mad_pass: bool,
    /// RMSD criterion passed
    pub rmsd_pass: bool,
    /// Correlation criterion passed
    pub correlation_pass: bool,
    /// Pass rate criterion passed
    pub pass_rate_pass: bool,
    /// Overall validation passed
    pub overall_pass: bool,
    /// List of failed molecules
    pub failed_molecules: Vec<String>,
    /// Diagnostic messages
    pub diagnostics: Vec<String>,
    /// Report timestamp
    pub timestamp: chrono::DateTime<chrono::Utc>,
}

impl ValidationReport {
    /// Generate formatted report string
    pub fn format_report(&self) -> String {
        let mut report = String::new();

        report.push_str("GW100 Mini-Validation Report\n");
        report.push_str(&format!(
            "Generated: {}\n\n",
            self.timestamp.format("%Y-%m-%d %H:%M:%S UTC")
        ));

        report.push_str("Statistical Summary:\n");
        report.push_str(&format!(
            "  Molecules tested: {}\n",
            self.statistics.n_molecules
        ));
        report.push_str(&format!("  MAD: {:.3} eV\n", self.statistics.mad));
        report.push_str(&format!("  RMSD: {:.3} eV\n", self.statistics.rmsd));
        report.push_str(&format!(
            "  Max deviation: {:.3} eV\n",
            self.statistics.max_dev
        ));
        report.push_str(&format!("  R²: {:.4}\n", self.statistics.r_squared));
        report.push_str(&format!(
            "  Pass rate: {:.1}%\n\n",
            self.statistics.pass_rate * 100.0
        ));

        report.push_str("Validation Results:\n");
        report.push_str(&format!(
            "  MAD ≤ 0.2 eV: {}\n",
            if self.mad_pass {
                "✓ PASS"
            } else {
                "✗ FAIL"
            }
        ));
        report.push_str(&format!(
            "  RMSD ≤ 0.24 eV: {}\n",
            if self.rmsd_pass {
                "✓ PASS"
            } else {
                "✗ FAIL"
            }
        ));
        report.push_str(&format!(
            "  R² ≥ 0.95: {}\n",
            if self.correlation_pass {
                "✓ PASS"
            } else {
                "✗ FAIL"
            }
        ));
        report.push_str(&format!(
            "  Pass rate ≥ 75%: {}\n\n",
            if self.pass_rate_pass {
                "✓ PASS"
            } else {
                "✗ FAIL"
            }
        ));

        report.push_str(&format!(
            "Overall Status: {}\n\n",
            if self.overall_pass {
                "✓ PASS"
            } else {
                "✗ FAIL"
            }
        ));

        if !self.failed_molecules.is_empty() {
            report.push_str("Failed Molecules:\n");
            for mol in &self.failed_molecules {
                report.push_str(&format!("  - {}\n", mol));
            }
            report.push('\n');
        }

        if !self.diagnostics.is_empty() {
            report.push_str("Diagnostics:\n");
            for diag in &self.diagnostics {
                report.push_str(&format!("  • {}\n", diag));
            }
        }

        report
    }

    /// Export report to JSON
    pub fn to_json(&self) -> Result<String> {
        serde_json::to_string_pretty(self)
            .map_err(|e| QuasixError::SerializationError(e.to_string()))
    }

    /// Export report to CSV
    pub fn to_csv(&self, results: &[MoleculeResult]) -> String {
        let mut csv = String::new();
        csv.push_str("Molecule,IP_Calc(eV),IP_Ref(eV),IP_Dev(eV),Pass\n");

        for result in results {
            csv.push_str(&format!(
                "{},{:.3},{:.3},{:.3},{}\n",
                result.molecule,
                result.ip_calc,
                result.ip_ref,
                result.ip_deviation(),
                if result.passes_threshold(0.2) {
                    "Yes"
                } else {
                    "No"
                }
            ));
        }

        csv
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_statistical_analyzer() {
        let config = BenchmarkConfig::default();
        let analyzer = StatisticalAnalyzer::new(config);

        // Create test results
        let results = vec![
            MoleculeResult {
                molecule: "H2O".to_string(),
                ip_calc: 12.5,
                ip_ref: 12.62,
                ea_calc: None,
                ea_ref: None,
                qp_energies: vec![-0.46, -0.39, -0.26],
                z_factors: vec![0.92, 0.91, 0.93],
                elapsed: Duration::from_secs(10),
                quality: QualityMetrics::default(),
                orbital_deviations: vec![-0.12, 0.08, -0.05],
            },
            MoleculeResult {
                molecule: "NH3".to_string(),
                ip_calc: 10.9,
                ip_ref: 10.82,
                ea_calc: None,
                ea_ref: None,
                qp_energies: vec![-0.40, -0.31],
                z_factors: vec![0.90, 0.91],
                elapsed: Duration::from_secs(12),
                quality: QualityMetrics::default(),
                orbital_deviations: vec![0.08, 0.05],
            },
        ];

        let stats = analyzer.analyze(&results);

        assert_eq!(stats.n_molecules, 2);
        assert!(stats.mad < 0.2);
        assert!(stats.r_squared > 0.9);
    }

    #[test]
    fn test_validation_pipeline() {
        let config = BenchmarkConfig::default();
        let pipeline = ValidationPipeline::new(config);

        // Create passing results
        let results = vec![MoleculeResult {
            molecule: "H2O".to_string(),
            ip_calc: 12.62,
            ip_ref: 12.62,
            ea_calc: None,
            ea_ref: None,
            qp_energies: vec![],
            z_factors: vec![],
            elapsed: Duration::from_secs(10),
            quality: QualityMetrics::default(),
            orbital_deviations: vec![],
        }];

        let report = pipeline.validate(&results);
        assert!(report.overall_pass);
        assert!(report.mad_pass);
        assert!(report.correlation_pass);
    }

    #[test]
    fn test_reference_database() {
        let db = ReferenceDatabase::load_gw100().unwrap();

        assert!(db.get("H2O").is_some());
        assert!(db.get("NH3").is_some());
        assert!(db.get("CO").is_some());
        assert!(db.get("benzene").is_some());

        let h2o = db.get("H2O").unwrap();
        assert!((h2o.ip - 12.62).abs() < 0.01);
    }

    #[test]
    fn test_outlier_detection() {
        let config = BenchmarkConfig::default();
        let analyzer = StatisticalAnalyzer::new(config);

        let values = vec![1.0, 1.1, 0.9, 1.2, 10.0]; // 10.0 is outlier
        let n_outliers = analyzer.detect_outliers(&values);
        assert_eq!(n_outliers, 1);
    }

    #[test]
    fn test_confidence_interval() {
        let config = BenchmarkConfig::default();
        let analyzer = StatisticalAnalyzer::new(config);

        let values = vec![0.1, -0.1, 0.15, -0.05, 0.0];
        let (lower, upper) = analyzer.compute_confidence_interval(&values);

        assert!(lower < 0.0);
        assert!(upper > 0.0);
        assert!((upper - lower) < 0.5);
    }
}
