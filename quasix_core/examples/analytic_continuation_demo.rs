//! Example demonstrating analytic continuation in QuasiX
//!
//! This example shows how to use the high-performance Rust implementation
//! for transforming GW self-energy from imaginary to real frequency axis.

use ndarray::Array1;
use num_complex::Complex64;
use quasix_core::freq::{ACConfig, AnalyticContinuationFitter, ModelType};

fn main() {
    println!("QuasiX Analytic Continuation Demo");
    println!("==================================\n");

    // Generate synthetic self-energy data on imaginary axis
    // In real GW calculations, this would come from imaginary-axis calculations
    let n_freq = 100;
    let xi = Array1::linspace(0.1, 20.0, n_freq);

    // Create a test self-energy with known pole structure
    // Σ(iξ) with poles at -1-0.5i and -3-1.0i
    let sigma_imag = xi.mapv(|x| {
        let z = Complex64::new(0.0, x);
        let pole1 = Complex64::new(-1.0, -0.5);
        let pole2 = Complex64::new(-3.0, -1.0);
        let residue1 = Complex64::new(0.5, 0.1);
        let residue2 = Complex64::new(0.3, -0.1);

        residue1 / (z - pole1) + residue2 / (z - pole2)
    });

    // Configure analytic continuation
    let config = ACConfig {
        max_poles: 10,
        max_pade_order: (10, 10),
        cv_fraction: 0.3,
        cv_iterations: 5,
        eta: 0.01,
        parallel: true,
        ..Default::default()
    };

    // Perform fitting with automatic model selection
    println!("Fitting analytic continuation models...");
    let mut fitter = AnalyticContinuationFitter::new(config);

    match fitter.fit(&xi, &sigma_imag) {
        Ok(()) => println!("✓ Fitting successful"),
        Err(e) => {
            eprintln!("✗ Fitting failed: {}", e);
            return;
        }
    }

    // Get the best model
    let result = fitter.get_result().unwrap();
    println!("\nSelected model: {:?}", result.model_type);
    println!("Cross-validation error: {:.2e}", result.cv_error);
    println!("Stability score: {:.2}", result.stability_score);

    // Extract poles and residues if available
    if let Some((poles, residues)) = result.model.get_poles_residues() {
        println!("\nFitted poles and residues:");
        for (i, (pole, residue)) in poles.iter().zip(residues.iter()).enumerate() {
            if residue.norm() > 1e-6 {
                // Only show significant poles
                println!(
                    "  Pole {}: {:.3} + {:.3}i, Residue: {:.3} + {:.3}i",
                    i + 1,
                    pole.re,
                    pole.im,
                    residue.re,
                    residue.im
                );
            }
        }
    }

    // Evaluate on real frequency axis
    println!("\nEvaluating on real frequency axis...");
    let omega = Array1::linspace(-10.0, 10.0, 200);
    let _sigma_real = fitter.evaluate_real_axis(&omega, Some(0.01)).unwrap();

    // Compute spectral function A(ω) = -Im[Σ(ω)]/π
    let spectral = fitter.spectral_function(&omega, Some(0.01)).unwrap();

    // Find peaks in spectral function
    let mut peaks = Vec::new();
    for i in 1..spectral.len() - 1 {
        if spectral[i] > spectral[i - 1]
            && spectral[i] > spectral[i + 1]
            && spectral[i] > 0.1 * spectral.iter().cloned().fold(0.0, f64::max)
        {
            peaks.push((omega[i], spectral[i]));
        }
    }

    println!("\nSpectral function peaks:");
    for (omega_peak, intensity) in peaks.iter() {
        println!("  ω = {:.2} eV, intensity = {:.3}", omega_peak, intensity);
    }

    // Validate causality
    println!("\nValidating causality and Kramers-Kronig relations...");
    let metrics = fitter.validate_causality(&omega).unwrap();
    println!(
        "  Causality violations: {} ({:.1}%)",
        metrics.violations,
        metrics.violation_fraction * 100.0
    );
    println!(
        "  Kramers-Kronig error: {:.2e} (relative: {:.2e})",
        metrics.kk_error, metrics.kk_relative_error
    );

    // Performance comparison between models
    println!("\nPerformance comparison:");
    if result.model_type == ModelType::Multipole {
        println!("  Multipole model selected for best accuracy");
    } else {
        println!("  Padé model selected for stability");
    }

    println!("\n✓ Analytic continuation completed successfully!");
}
