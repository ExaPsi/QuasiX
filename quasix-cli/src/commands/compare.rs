//! CLI command for CD vs AC method comparison

use clap::{Args, ValueEnum};
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

#[derive(Args, Debug)]
pub struct CompareCommand {
    /// Molecules to compare (H2O, NH3, CO, etc.)
    #[clap(required = true, value_delimiter = ',')]
    pub molecules: Vec<String>,

    /// Basis set to use
    #[clap(short, long, default_value = "cc-pVDZ")]
    pub basis: String,

    /// Auxiliary basis for DF (auto-select if not specified)
    #[clap(short = 'a', long)]
    pub aux_basis: Option<String>,

    /// Number of frequency grid points
    #[clap(short = 'n', long, default_value = "40")]
    pub n_points: usize,

    /// Maximum imaginary frequency
    #[clap(short = 'x', long, default_value = "50.0")]
    pub xi_max: f64,

    /// Broadening parameter eta
    #[clap(short = 'e', long, default_value = "0.01")]
    pub eta: f64,

    /// MAD threshold in eV for validation
    #[clap(short = 'm', long, default_value = "0.05")]
    pub mad_threshold: f64,

    /// Number of bootstrap samples for confidence intervals
    #[clap(short = 's', long, default_value = "1000")]
    pub bootstrap_samples: usize,

    /// Parallel execution of CD and AC
    #[clap(short = 'p', long, default_value = "true")]
    pub parallel: bool,

    /// Number of threads to use (auto-detect if not specified)
    #[clap(short = 't', long)]
    pub threads: Option<usize>,

    /// Output format
    #[clap(short = 'f', long, value_enum, default_value = "html")]
    pub format: OutputFormat,

    /// Output file path
    #[clap(short = 'o', long, default_value = "cd_vs_ac_report")]
    pub output: PathBuf,

    /// Verbose output
    #[clap(short = 'v', long, action = clap::ArgAction::Count)]
    pub verbose: u8,

    /// Save intermediate results
    #[clap(long)]
    pub save_intermediates: bool,

    /// Generate plots
    #[clap(long, default_value = "true")]
    pub plot: bool,

    /// JSON configuration file (overrides CLI arguments)
    #[clap(short = 'c', long)]
    pub config: Option<PathBuf>,
}

#[derive(Debug, Clone, Copy, ValueEnum, Serialize, Deserialize)]
pub enum OutputFormat {
    Html,
    Json,
    Text,
    All,
}

impl CompareCommand {
    pub fn execute(&self) -> anyhow::Result<()> {
        use anyhow::Context;
        use quasix_core::validation::{
            ACParameters, CDParameters, ComparisonHarness, UnifiedParameters, ValidationConfig,
        };
        use std::time::Instant;

        // Print banner
        println!("\n═══════════════════════════════════════════════════════");
        println!("       CD vs AC Comparison Harness");
        println!("═══════════════════════════════════════════════════════\n");

        // Load configuration if provided
        let config = if let Some(config_path) = &self.config {
            println!("Loading configuration from: {}", config_path.display());
            let config_str = std::fs::read_to_string(config_path)
                .context("Failed to read configuration file")?;
            serde_json::from_str::<ValidationConfig>(&config_str)
                .context("Failed to parse configuration file")?
        } else {
            // Build configuration from CLI arguments
            let mut unified_params = UnifiedParameters {
                n_freq_common: self.n_points,
                energy_range: (-30.0, 30.0),
                eta_common: self.eta,
                temperature: 0.025, // ~300K
                chemical_potential: 0.0,
                cd_params: CDParameters::default(),
                ac_params: ACParameters::default(),
            };

            let cd_params = CDParameters {
                eta: self.eta,
                n_freq: self.n_points,
                omega_max: 100.0,
                tolerance: 1e-6,
                adaptive_grid: true,
                n_segments: 4,
                rotation_angle: 0.0,
            };

            let ac_params = ACParameters {
                n_imag: self.n_points,
                xi_max: self.xi_max,
                n_poles: 20,
                regularization: 1e-8,
                model_type: "pade".to_string(),
                convergence_tol: 1e-6,
                max_iter: 100,
            };

            // Update unified params with custom CD and AC params
            unified_params.cd_params = cd_params;
            unified_params.ac_params = ac_params;

            ValidationConfig {
                parameters: unified_params,
                n_threads: self.threads.unwrap_or(num_cpus::get()),
                thread_allocation: [0.4, 0.4, 0.2],
                enable_profiling: true,
                save_intermediates: self.save_intermediates,
                max_memory_gb: 16.0,
            }
        };

        println!("Configuration:");
        println!("  Molecules: {:?}", self.molecules);
        println!("  Basis set: {}", self.basis);
        println!("  Frequency points: {}", config.parameters.n_freq_common);
        println!("  MAD threshold: {} eV", self.mad_threshold);
        println!("  Parallel execution: {}", self.parallel);
        if config.n_threads > 0 {
            println!("  Threads: {}", config.n_threads);
        }
        println!();

        // Create comparison harness
        let harness = ComparisonHarness::new(config.clone())?;

        // Process each molecule
        let mut all_results = Vec::new();
        let start_time = Instant::now();

        for molecule_name in &self.molecules {
            println!("Processing molecule: {}", molecule_name);

            // Load molecule data (mock for now - would load from file or database)
            let (df_tensors, w_screened, mo_energies, mo_occ) =
                load_molecule_data(molecule_name, &self.basis, &self.aux_basis)?;

            // Run comparison
            let result = harness.run_comparison(df_tensors, w_screened, mo_energies, mo_occ)?;

            // Print immediate results
            let cd_qp = &result.cd_results.qp_energies;
            let ac_qp = &result.ac_results.qp_energies;
            println!(
                "  CD QP energies: {:?}",
                cd_qp.iter().take(5).collect::<Vec<_>>()
            );
            println!(
                "  AC QP energies: {:?}",
                ac_qp.iter().take(5).collect::<Vec<_>>()
            );
            println!("  MAD: {:.4} eV", result.statistics.mad);
            println!("  RMSD: {:.4} eV", result.statistics.rmsd);
            println!("  Max deviation: {:.4} eV", result.statistics.max_error);
            println!("  Correlation: {:.4}", result.statistics.correlation);
            println!(
                "  CD timing: {:.2} s",
                result.cd_results.wall_time.as_secs_f64()
            );
            println!(
                "  AC timing: {:.2} s",
                result.ac_results.wall_time.as_secs_f64()
            );
            println!(
                "  Speedup: {:.1}x",
                result.cd_results.wall_time.as_secs_f64()
                    / result.ac_results.wall_time.as_secs_f64()
            );

            if result.validation_passed {
                println!(
                    "  ✓ Methods agree: MAD = {:.4} eV < {:.4} eV",
                    result.statistics.mad, self.mad_threshold
                );
            } else {
                println!(
                    "  ✗ Methods disagree: MAD = {:.4} eV > {:.4} eV",
                    result.statistics.mad, self.mad_threshold
                );
            }

            if !result.statistics.outlier_indices.is_empty() {
                println!(
                    "  ⚠ {} outliers detected: {:?}",
                    result.statistics.outlier_indices.len(),
                    result.statistics.outlier_indices
                );
            }
            println!();

            all_results.push((molecule_name.clone(), result));
        }

        let total_time = start_time.elapsed();

        // Generate reports based on output format
        match self.format {
            OutputFormat::All => {
                // Generate all formats
                let html_path = self.output.with_extension("html");
                generate_html_report(&all_results, &config, &html_path)?;
                println!("HTML report saved to: {}", html_path.display());

                let json_path = self.output.with_extension("json");
                save_json_results(&all_results, &config, &json_path)?;
                println!("JSON results saved to: {}", json_path.display());

                let text_path = self.output.with_extension("txt");
                generate_text_report(&all_results, &config, &text_path)?;
                println!("Text report saved to: {}", text_path.display());
            }
            OutputFormat::Html => {
                let html_path = self.output.with_extension("html");
                generate_html_report(&all_results, &config, &html_path)?;
                println!("HTML report saved to: {}", html_path.display());
            }
            OutputFormat::Json => {
                let json_path = self.output.with_extension("json");
                save_json_results(&all_results, &config, &json_path)?;
                println!("JSON results saved to: {}", json_path.display());
            }
            OutputFormat::Text => {
                let text_path = self.output.with_extension("txt");
                generate_text_report(&all_results, &config, &text_path)?;
                println!("Text report saved to: {}", text_path.display());
            }
        }

        // Generate plots if requested
        if self.plot {
            let plots_dir = self.output.with_extension("plots");
            std::fs::create_dir_all(&plots_dir)?;
            // Note: Actual plotting would require Python integration
            println!("Plots would be saved to: {}", plots_dir.display());
        }

        // Summary statistics
        println!("\n=== Summary ===");
        let total_passed = all_results
            .iter()
            .filter(|(_, r)| r.validation_passed)
            .count();
        let total_molecules = all_results.len();
        println!("Molecules tested: {}", total_molecules);
        println!("Passed validation: {}/{}", total_passed, total_molecules);

        let overall_mad: f64 = all_results
            .iter()
            .map(|(_, r)| r.statistics.mad)
            .sum::<f64>()
            / total_molecules as f64;
        println!("Overall MAD: {:.4} eV", overall_mad);

        let avg_speedup: f64 = all_results
            .iter()
            .map(|(_, r)| {
                r.cd_results.wall_time.as_secs_f64() / r.ac_results.wall_time.as_secs_f64()
            })
            .sum::<f64>()
            / total_molecules as f64;
        println!("Average speedup (AC vs CD): {:.1}x", avg_speedup);
        println!("Total time: {:.2} s", total_time.as_secs_f64());

        if overall_mad < self.mad_threshold {
            println!(
                "\n✓ Overall validation PASSED: MAD = {:.4} eV < {:.4} eV",
                overall_mad, self.mad_threshold
            );
        } else {
            println!(
                "✗ Overall validation FAILED: MAD = {:.4} eV > {:.4} eV",
                overall_mad, self.mad_threshold
            );
        }

        Ok(())
    }
}

fn load_molecule_data(
    molecule: &str,
    _basis: &str,
    _aux_basis: &Option<String>,
) -> anyhow::Result<(
    ndarray::Array3<f64>,
    ndarray::Array3<num_complex::Complex64>,
    ndarray::Array1<f64>,
    ndarray::Array1<f64>,
)> {
    use ndarray::{Array1, Array3};
    use num_complex::Complex64;

    // Mock implementation - in reality would load from HDF5 or compute
    let (n_basis, n_aux, n_occ, n_virt) = match molecule {
        "H2O" => (13, 38, 5, 8),  // cc-pVDZ
        "NH3" => (18, 50, 5, 13), // cc-pVDZ
        "CO" => (20, 60, 7, 13),  // cc-pVDZ
        _ => (10, 30, 5, 5),      // Default
    };

    let n_mo = n_occ + n_virt;

    // Create mock tensors with physically reasonable values
    let mut df_tensors = Array3::<f64>::zeros((n_mo, n_mo, n_aux));
    let mut w_screened = Array3::<Complex64>::zeros((40, n_aux, n_aux)); // (n_freq, n_aux, n_aux)
    let mut mo_energies = Array1::<f64>::zeros(n_mo);
    let mut mo_occ = Array1::<f64>::zeros(n_mo);

    // Fill with mock data
    for i in 0..n_mo {
        for j in 0..n_mo {
            for p in 0..n_aux {
                df_tensors[[i, j, p]] =
                    (-((i as f64 - j as f64).powi(2) + (p as f64 / 10.0).powi(2)) / 10.0).exp();
            }
        }

        // HOMO-LUMO gap around 8-10 eV
        if i < n_occ {
            mo_energies[i] = -15.0 + (i as f64) * 2.0; // Occupied
            mo_occ[i] = 2.0;
        } else {
            mo_energies[i] = 2.0 + (i as f64 - n_occ as f64) * 1.5; // Virtual
            mo_occ[i] = 0.0;
        }
    }

    // Mock screened Coulomb - shape is (n_freq, n_aux, n_aux)
    for w in 0..40 {
        for p in 0..n_aux {
            for q in 0..n_aux {
                let freq = w as f64 * 2.0;
                let real_part =
                    1.0 / (1.0 + ((p as f64 - q as f64).powi(2) + freq.powi(2)) / 100.0);
                w_screened[[w, p, q]] = Complex64::new(real_part, 0.0);
            }
        }
    }

    println!(
        "  Loaded {} with {} basis functions, {} auxiliary functions",
        molecule, n_basis, n_aux
    );

    Ok((df_tensors, w_screened, mo_energies, mo_occ))
}

fn generate_html_report(
    results: &[(String, quasix_core::validation::ComparisonResult)],
    config: &quasix_core::validation::ValidationConfig,
    path: &PathBuf,
) -> anyhow::Result<()> {
    use std::fs::File;
    use std::io::Write;

    let mut html = String::new();
    let freq_points = config.parameters.n_freq_common;
    html.push_str(&format!(
        r#"<!DOCTYPE html>
<html>
<head>
    <title>CD vs AC Comparison Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        h1 {{ color: #333; }}
        table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #f2f2f2; }}
        .passed {{ color: green; font-weight: bold; }}
        .failed {{ color: red; font-weight: bold; }}
        .summary {{ background-color: #f9f9f9; padding: 15px; border-radius: 5px; }}
    </style>
</head>
<body>
    <h1>CD vs AC Method Comparison Report</h1>
    <div class="summary">
        <h2>Configuration</h2>
        <ul>
            <li>Frequency points: {}</li>
            <li>MAD threshold: {} eV</li>
            <li>Bootstrap samples: {}</li>
            <li>Parallel execution: {}</li>
        </ul>
    </div>
"#,
        freq_points,
        0.05, // MAD threshold hard-coded for now
        1000, // Bootstrap samples hard-coded
        true  // Parallel execution hard-coded
    ));

    html.push_str(
        r#"
    <h2>Results Summary</h2>
    <table>
        <tr>
            <th>Molecule</th>
            <th>MAD (eV)</th>
            <th>RMSD (eV)</th>
            <th>Max Dev (eV)</th>
            <th>Correlation</th>
            <th>CD Time (s)</th>
            <th>AC Time (s)</th>
            <th>Speedup</th>
            <th>Status</th>
        </tr>
"#,
    );

    for (molecule, result) in results {
        let status_class = if result.validation_passed {
            "passed"
        } else {
            "failed"
        };
        let status_text = if result.validation_passed {
            "PASSED"
        } else {
            "FAILED"
        };

        html.push_str(&format!(
            r#"
        <tr>
            <td>{}</td>
            <td>{:.4}</td>
            <td>{:.4}</td>
            <td>{:.4}</td>
            <td>{:.4}</td>
            <td>{:.2}</td>
            <td>{:.2}</td>
            <td>{:.1}x</td>
            <td class="{}">{}</td>
        </tr>
"#,
            molecule,
            result.statistics.mad,
            result.statistics.rmsd,
            result.statistics.max_error,
            result.statistics.correlation,
            result.cd_results.wall_time.as_secs_f64(),
            result.ac_results.wall_time.as_secs_f64(),
            result.cd_results.wall_time.as_secs_f64() / result.ac_results.wall_time.as_secs_f64(),
            status_class,
            status_text
        ));
    }

    html.push_str("    </table>\n");

    // Overall statistics
    let total_passed = results.iter().filter(|(_, r)| r.validation_passed).count();
    let overall_mad: f64 =
        results.iter().map(|(_, r)| r.statistics.mad).sum::<f64>() / results.len() as f64;

    html.push_str(&format!(
        r#"
    <div class="summary">
        <h2>Overall Statistics</h2>
        <ul>
            <li>Molecules tested: {}</li>
            <li>Passed validation: {}/{}</li>
            <li>Overall MAD: {:.4} eV</li>
            <li>Pass rate: {:.1}%</li>
        </ul>
    </div>
</body>
</html>
"#,
        results.len(),
        total_passed,
        results.len(),
        overall_mad,
        (total_passed as f64 / results.len() as f64) * 100.0
    ));

    let mut file = File::create(path)?;
    file.write_all(html.as_bytes())?;

    Ok(())
}

fn save_json_results(
    results: &[(String, quasix_core::validation::ComparisonResult)],
    config: &quasix_core::validation::ValidationConfig,
    path: &PathBuf,
) -> anyhow::Result<()> {
    use serde_json::json;
    use std::fs::File;

    let json_data = json!({
        "configuration": config,
        "results": results.iter().map(|(mol, res)| {
            json!({
                "molecule": mol,
                "mad": res.statistics.mad,
                "rmsd": res.statistics.rmsd,
                "max_deviation": res.statistics.max_error,
                "correlation": res.statistics.correlation,
                "cd_timing": res.cd_results.wall_time.as_secs_f64(),
                "ac_timing": res.ac_results.wall_time.as_secs_f64(),
                "validation_passed": res.validation_passed,
                "outliers": res.statistics.outlier_indices,
            })
        }).collect::<Vec<_>>(),
        "summary": {
            "total_molecules": results.len(),
            "passed": results.iter().filter(|(_, r)| r.validation_passed).count(),
            "overall_mad": results.iter().map(|(_, r)| r.statistics.mad).sum::<f64>() / results.len() as f64,
        }
    });

    let file = File::create(path)?;
    serde_json::to_writer_pretty(file, &json_data)?;

    Ok(())
}

fn generate_text_report(
    results: &[(String, quasix_core::validation::ComparisonResult)],
    config: &quasix_core::validation::ValidationConfig,
    path: &PathBuf,
) -> anyhow::Result<()> {
    use std::fs::File;
    use std::io::Write;

    let mut file = File::create(path)?;

    writeln!(file, "CD vs AC Method Comparison Report")?;
    writeln!(file, "==================================")?;
    writeln!(file)?;
    writeln!(file, "Configuration:")?;
    writeln!(
        file,
        "  Frequency points: {}",
        config.parameters.n_freq_common
    )?;
    writeln!(file, "  MAD threshold: {} eV", 0.05)?; // Hard-coded for now
    writeln!(file)?;

    for (molecule, result) in results {
        writeln!(file, "Molecule: {}", molecule)?;
        writeln!(file, "  MAD: {:.4} eV", result.statistics.mad)?;
        writeln!(file, "  RMSD: {:.4} eV", result.statistics.rmsd)?;
        writeln!(
            file,
            "  Max deviation: {:.4} eV",
            result.statistics.max_error
        )?;
        writeln!(file, "  Correlation: {:.4}", result.statistics.correlation)?;
        writeln!(
            file,
            "  Status: {}",
            if result.validation_passed {
                "PASSED"
            } else {
                "FAILED"
            }
        )?;
        writeln!(file)?;
    }

    Ok(())
}
