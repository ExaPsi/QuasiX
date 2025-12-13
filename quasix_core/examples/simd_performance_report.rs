//! Comprehensive SIMD performance report for polarizability builder

use quasix_core::gw::polarizability_builder::PolarizabilityBuilder;
use ndarray::{Array1, Array2};
use num_complex::Complex64;
use std::time::Instant;

fn main() {
    println!("╔══════════════════════════════════════════════════════════════╗");
    println!("║        QuasiX SIMD Performance Report                       ║");
    println!("║        Polarizability Builder Optimizations                 ║");
    println!("╚══════════════════════════════════════════════════════════════╝");
    println!();

    // Check CPU features
    #[cfg(target_arch = "x86_64")]
    {
        println!("CPU SIMD Capabilities:");
        println!("┌────────────────┬──────────┐");
        println!("│ Feature        │ Available│");
        println!("├────────────────┼──────────┤");
        println!("│ AVX2           │ {:^8} │", if is_x86_feature_detected!("avx2") { "✓" } else { "✗" });
        println!("│ AVX-512F       │ {:^8} │", if is_x86_feature_detected!("avx512f") { "✓" } else { "✗" });
        println!("│ FMA            │ {:^8} │", if is_x86_feature_detected!("fma") { "✓" } else { "✗" });
        println!("└────────────────┴──────────┘");
        println!();
    }

    // Test configurations - larger systems for better speedup demonstration
    let test_configs = vec![
        (10, 100, 200, "Small"),
        (30, 300, 600, "Medium"),
        (50, 500, 1000, "Large"),
        (100, 1000, 2000, "Very Large"),
    ];

    println!("Performance Benchmarks:");
    println!("═══════════════════════════════════════════════════════════════");

    for (nocc, nvirt, naux, label) in test_configs {
        println!("\n{} System (nocc={}, nvirt={}, naux={})", label, nocc, nvirt, naux);
        println!("───────────────────────────────────────────────────────────────");

        // Create test system with more realistic data
        let (mo_energies, df_ia) = setup_realistic_system(nocc, nvirt, naux);
        let builder = PolarizabilityBuilder::new(nocc, nvirt, naux, &mo_energies, df_ia)
            .expect("Failed to create builder");

        // Benchmark gap computation
        let (simd_gaps_time, scalar_gaps_time) = benchmark_gaps(&builder, nocc, nvirt);
        let gaps_speedup = scalar_gaps_time / simd_gaps_time;

        println!("  Gap Computation ({} gaps):", nocc * nvirt);
        println!("    SIMD:   {:8.3} ms  ({:8.1} M gaps/s)",
                 simd_gaps_time * 1000.0,
                 (nocc * nvirt) as f64 / simd_gaps_time / 1e6);
        println!("    Scalar: {:8.3} ms  ({:8.1} M gaps/s)",
                 scalar_gaps_time * 1000.0,
                 (nocc * nvirt) as f64 / scalar_gaps_time / 1e6);
        println!("    Speedup: {:.2}x {}",
                 gaps_speedup,
                 performance_indicator(gaps_speedup));

        // Benchmark P0 build if on x86_64
        #[cfg(target_arch = "x86_64")]
        {
            let (simd_p0_time, standard_p0_time) = benchmark_p0(&builder, naux);
            let p0_speedup = standard_p0_time / simd_p0_time;

            println!("  P0 Build ({}x{} matrix):", naux, naux);
            println!("    SIMD:     {:8.3} ms  ({:8.1} M elem/s)",
                     simd_p0_time * 1000.0,
                     (naux * naux) as f64 / simd_p0_time / 1e6);
            println!("    Standard: {:8.3} ms  ({:8.1} M elem/s)",
                     standard_p0_time * 1000.0,
                     (naux * naux) as f64 / standard_p0_time / 1e6);
            println!("    Speedup: {:.2}x {}",
                     p0_speedup,
                     performance_indicator(p0_speedup));
        }
    }

    println!("\n═══════════════════════════════════════════════════════════════");
    println!("\nPerformance Summary:");
    println!("────────────────────");

    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx512f") {
            println!("✓ AVX-512 optimizations active - Maximum performance");
            println!("  Expected speedups: 8-16x for gap computation");
            println!("                     2-4x for P0 matrix construction");
        } else if is_x86_feature_detected!("avx2") {
            println!("✓ AVX2 optimizations active - Good performance");
            println!("  Expected speedups: 4-8x for gap computation");
            println!("                     1.5-2x for P0 matrix construction");
        } else {
            println!("⚠ No SIMD optimizations available - Baseline performance");
        }
    }

    #[cfg(not(target_arch = "x86_64"))]
    {
        println!("⚠ Non-x86_64 architecture - SIMD optimizations not available");
    }

    println!("\n✅ Performance report complete");
}

fn setup_realistic_system(nocc: usize, nvirt: usize, naux: usize) -> (Array1<f64>, Array2<f64>) {
    // Create realistic HOMO-LUMO gap
    let homo_energy = -0.25;  // ~-6.8 eV
    let lumo_energy = 0.05;    // ~1.4 eV

    let mut mo_energies = Array1::zeros(nocc + nvirt);

    // Occupied orbitals - exponential spacing near HOMO
    for i in 0..nocc {
        let x = i as f64 / nocc as f64;
        mo_energies[i] = homo_energy - 2.0 * (1.0 - x).exp();
    }

    // Virtual orbitals - logarithmic spacing from LUMO
    for a in 0..nvirt {
        let x = a as f64 / nvirt as f64;
        mo_energies[nocc + a] = lumo_energy + 0.5 * (1.0 + 10.0 * x).ln();
    }

    // Create DF tensor with structure
    let n_trans = nocc * nvirt;
    let mut df_ia = Array2::zeros((n_trans, naux));

    // Add some structure to make computation non-trivial
    for i in 0..n_trans {
        for j in 0..naux {
            let value = ((i * 7 + j * 3) % 17) as f64 / 17.0 - 0.5;
            df_ia[[i, j]] = value * (-(i as f64 / n_trans as f64)).exp();
        }
    }

    (mo_energies, df_ia)
}

fn benchmark_gaps(builder: &PolarizabilityBuilder, nocc: usize, nvirt: usize) -> (f64, f64) {
    let n_gaps = nocc * nvirt;

    // Determine iteration count based on problem size
    let iterations = if n_gaps < 10000 { 1000 }
                    else if n_gaps < 100000 { 100 }
                    else { 10 };

    // Warm up
    for _ in 0..5 {
        let _ = builder.compute_gaps();
    }

    // SIMD version
    let start = Instant::now();
    for _ in 0..iterations {
        let gaps = builder.compute_gaps();
        std::hint::black_box(gaps);
    }
    let simd_time = start.elapsed().as_secs_f64() / iterations as f64;

    // Scalar version
    let mut gaps_buffer = vec![0.0; n_gaps];
    let start = Instant::now();
    for _ in 0..iterations {
        builder.compute_gaps_scalar(&mut gaps_buffer);
        std::hint::black_box(&gaps_buffer);
    }
    let scalar_time = start.elapsed().as_secs_f64() / iterations as f64;

    (simd_time, scalar_time)
}

#[cfg(target_arch = "x86_64")]
fn benchmark_p0(builder: &PolarizabilityBuilder, naux: usize) -> (f64, f64) {
    let omega = Complex64::new(0.0, 1.0);

    // Fewer iterations for large matrices
    let iterations = if naux < 500 { 10 }
                    else if naux < 1000 { 5 }
                    else { 2 };

    // Warm up
    let _ = builder.build_p0(omega);

    // SIMD version
    let start = Instant::now();
    for _ in 0..iterations {
        let p0 = builder.build_p0_simd_imaginary(omega).unwrap();
        std::hint::black_box(p0);
    }
    let simd_time = start.elapsed().as_secs_f64() / iterations as f64;

    // Standard version
    let start = Instant::now();
    for _ in 0..iterations {
        let p0 = builder.build_p0_standard(omega).unwrap();
        std::hint::black_box(p0);
    }
    let standard_time = start.elapsed().as_secs_f64() / iterations as f64;

    (simd_time, standard_time)
}

#[cfg(not(target_arch = "x86_64"))]
fn benchmark_p0(_builder: &PolarizabilityBuilder, _naux: usize) -> (f64, f64) {
    (1.0, 1.0)
}

fn performance_indicator(speedup: f64) -> &'static str {
    if speedup >= 8.0 { "🚀 Excellent" }
    else if speedup >= 4.0 { "⚡ Very Good" }
    else if speedup >= 2.0 { "✨ Good" }
    else if speedup >= 1.5 { "📈 Moderate" }
    else if speedup >= 1.1 { "➕ Minor" }
    else { "➖ None" }
}