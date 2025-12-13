//! Performance benchmarks for spectral analysis module
//!
//! Tests SIMD optimization effectiveness for spectral function computation

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use ndarray::Array1;
use quasix_core::analysis::{BroadeningParams, SpectralAnalyzer};

/// Benchmark spectral function computation with varying system sizes
fn bench_spectral_function(c: &mut Criterion) {
    let mut group = c.benchmark_group("spectral_function");

    // Test different system sizes
    for n_states in [10, 50, 100, 500].iter() {
        let energies = Array1::linspace(-5.0, 5.0, *n_states);
        let z_factors = Array1::from_elem(*n_states, 0.8);
        let fermi_energy = 0.0;

        let analyzer = SpectralAnalyzer::new(energies, z_factors, fermi_energy)
            .expect("Failed to create analyzer");

        // Benchmark with 1000-point omega grid
        let omega_grid = Array1::linspace(-10.0, 10.0, 1000);

        group.bench_function(
            BenchmarkId::new("compute", n_states),
            |b| {
                b.iter(|| {
                    analyzer
                        .compute_spectral_function(black_box(&omega_grid))
                        .expect("Failed to compute")
                });
            },
        );
    }

    group.finish();
}

/// Benchmark different broadening schemes
fn bench_broadening_schemes(c: &mut Criterion) {
    let mut group = c.benchmark_group("broadening_schemes");

    // Setup: 100-state system
    let n_states = 100;
    let energies = Array1::linspace(-5.0, 5.0, n_states);
    let z_factors = Array1::from_elem(n_states, 0.8);
    let omega_grid = Array1::linspace(-10.0, 10.0, 500);

    // Pure Gaussian
    {
        let mut analyzer = SpectralAnalyzer::new(
            energies.clone(),
            z_factors.clone(),
            0.0,
        ).expect("Failed to create analyzer");

        let mut params = BroadeningParams::default();
        params.voigt_mixing = 0.0; // Pure Gaussian
        analyzer.set_broadening(params).expect("Failed to set broadening");

        group.bench_function("gaussian", |b| {
            b.iter(|| {
                analyzer
                    .compute_spectral_function(black_box(&omega_grid))
                    .expect("Failed to compute")
            });
        });
    }

    // Pure Lorentzian
    {
        let mut analyzer = SpectralAnalyzer::new(
            energies.clone(),
            z_factors.clone(),
            0.0,
        ).expect("Failed to create analyzer");

        let mut params = BroadeningParams::default();
        params.voigt_mixing = 1.0; // Pure Lorentzian
        analyzer.set_broadening(params).expect("Failed to set broadening");

        group.bench_function("lorentzian", |b| {
            b.iter(|| {
                analyzer
                    .compute_spectral_function(black_box(&omega_grid))
                    .expect("Failed to compute")
            });
        });
    }

    // Voigt profile (mixed)
    {
        let mut analyzer = SpectralAnalyzer::new(
            energies.clone(),
            z_factors.clone(),
            0.0,
        ).expect("Failed to create analyzer");

        let params = BroadeningParams::default(); // Default is 0.5 mixing
        analyzer.set_broadening(params).expect("Failed to set broadening");

        group.bench_function("voigt", |b| {
            b.iter(|| {
                analyzer
                    .compute_spectral_function(black_box(&omega_grid))
                    .expect("Failed to compute")
            });
        });
    }

    group.finish();
}

/// Benchmark PES/IPES spectrum generation
fn bench_pes_ipes(c: &mut Criterion) {
    let mut group = c.benchmark_group("pes_ipes");

    // Setup: 200-state system with half occupied
    let n_states = 200;
    let energies = Array1::linspace(-5.0, 5.0, n_states);
    let z_factors = Array1::from_elem(n_states, 0.8);
    let fermi_energy = 0.0;

    let analyzer = SpectralAnalyzer::new(energies, z_factors, fermi_energy)
        .expect("Failed to create analyzer");

    group.bench_function("pes_generation", |b| {
        b.iter(|| {
            analyzer
                .generate_pes_spectrum(black_box(0.0), black_box(5.0), black_box(500))
                .expect("Failed to generate PES")
        });
    });

    group.bench_function("ipes_generation", |b| {
        b.iter(|| {
            analyzer
                .generate_ipes_spectrum(black_box(0.0), black_box(5.0), black_box(500))
                .expect("Failed to generate IPES")
        });
    });

    group.finish();
}

/// Benchmark SIMD vs scalar implementation
#[cfg(target_arch = "x86_64")]
fn bench_simd_effectiveness(c: &mut Criterion) {
    let mut group = c.benchmark_group("simd_effectiveness");

    // Large system to show SIMD benefits
    let n_states = 1000;
    let energies = Array1::linspace(-10.0, 10.0, n_states);
    let z_factors = Array1::from_elem(n_states, 0.8);
    let omega_grid = Array1::linspace(-15.0, 15.0, 2000);

    let mut analyzer = SpectralAnalyzer::new(energies, z_factors, 0.0)
        .expect("Failed to create analyzer");

    // Use Gaussian broadening to test SIMD path
    let mut params = BroadeningParams::default();
    params.voigt_mixing = 0.0; // Pure Gaussian for SIMD
    params.gaussian_width = 0.05; // Narrow width for more computation
    analyzer.set_broadening(params).expect("Failed to set broadening");

    group.bench_function("simd_gaussian_large", |b| {
        b.iter(|| {
            analyzer
                .compute_spectral_function(black_box(&omega_grid))
                .expect("Failed to compute")
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_spectral_function,
    bench_broadening_schemes,
    bench_pes_ipes,
);

#[cfg(target_arch = "x86_64")]
criterion_group!(
    simd_benches,
    bench_simd_effectiveness,
);

#[cfg(not(target_arch = "x86_64"))]
criterion_main!(benches);

#[cfg(target_arch = "x86_64")]
criterion_main!(benches, simd_benches);