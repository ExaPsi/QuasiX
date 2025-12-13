//! Example demonstrating SIMD optimizations in polarizability builder

use quasix_core::gw::polarizability_builder::PolarizabilityBuilder;
use ndarray::{Array1, Array2};
use num_complex::Complex64;
use std::time::Instant;

fn main() {
    println!("QuasiX Polarizability Builder SIMD Performance Test");
    println!("===================================================\n");

    // Check CPU features
    #[cfg(target_arch = "x86_64")]
    {
        println!("Detected CPU SIMD Features:");
        println!("  AVX2:     {}", is_x86_feature_detected!("avx2"));
        println!("  AVX-512F: {}", is_x86_feature_detected!("avx512f"));
        println!("  FMA:      {}", is_x86_feature_detected!("fma"));
        println!();
    }

    // Test configuration
    let nocc = 20;
    let nvirt = 100;
    let naux = 200;

    println!("Test System:");
    println!("  Occupied orbitals: {}", nocc);
    println!("  Virtual orbitals:  {}", nvirt);
    println!("  Auxiliary basis:   {}", naux);
    println!("  Total gaps:        {}", nocc * nvirt);
    println!();

    // Create test system
    let mut mo_energies = Array1::zeros(nocc + nvirt);
    for i in 0..nocc {
        mo_energies[i] = -2.0 + 0.1 * i as f64;
    }
    for a in 0..nvirt {
        mo_energies[nocc + a] = 0.1 + 0.05 * a as f64;
    }

    // Create simple DF tensor
    let n_trans = nocc * nvirt;
    let mut df_ia = Array2::zeros((n_trans, naux));
    for i in 0..n_trans {
        for j in 0..naux {
            df_ia[[i, j]] = 0.1 * ((i + j) % 7) as f64;
        }
    }

    let builder = PolarizabilityBuilder::new(nocc, nvirt, naux, &mo_energies, df_ia)
        .expect("Failed to create builder");

    // Benchmark gap computation
    println!("Benchmarking Gap Computation:");
    println!("------------------------------");

    let iterations = 1000;

    // SIMD version
    let start = Instant::now();
    for _ in 0..iterations {
        let gaps = builder.compute_gaps();
        std::hint::black_box(gaps);
    }
    let simd_duration = start.elapsed();

    // Scalar version
    let mut gaps_buffer = vec![0.0; nocc * nvirt];
    let start = Instant::now();
    for _ in 0..iterations {
        builder.compute_gaps_scalar(&mut gaps_buffer);
        std::hint::black_box(&gaps_buffer);
    }
    let scalar_duration = start.elapsed();

    let speedup = scalar_duration.as_secs_f64() / simd_duration.as_secs_f64();

    println!("  Iterations: {}", iterations);
    println!("  SIMD time:   {:.3} ms", simd_duration.as_millis());
    println!("  Scalar time: {:.3} ms", scalar_duration.as_millis());
    println!("  Speedup:     {:.2}x", speedup);
    println!();

    // Verify correctness
    println!("Verifying Correctness:");
    println!("----------------------");

    let gaps_simd = builder.compute_gaps();
    let mut gaps_scalar = vec![0.0; nocc * nvirt];
    builder.compute_gaps_scalar(&mut gaps_scalar);

    let mut max_diff: f64 = 0.0;
    for (i, (&simd_val, &scalar_val)) in gaps_simd.iter().zip(gaps_scalar.iter()).enumerate() {
        let diff = (simd_val - scalar_val).abs();
        max_diff = max_diff.max(diff);
        if diff > 1e-14 {
            println!("  WARNING: Gap {} differs: SIMD={:.6e}, Scalar={:.6e}, diff={:.6e}",
                     i, simd_val, scalar_val, diff);
        }
    }

    if max_diff < 1e-14 {
        println!("  ✓ All gaps match within tolerance (max diff: {:.6e})", max_diff);
    } else {
        println!("  ✗ Some gaps differ (max diff: {:.6e})", max_diff);
    }
    println!();

    // Test P0 build
    #[cfg(target_arch = "x86_64")]
    {
        println!("Benchmarking P0 Build (Imaginary Frequency):");
        println!("---------------------------------------------");

        let omega = Complex64::new(0.0, 1.0);
        let iterations = 10;

        // SIMD version
        let start = Instant::now();
        for _ in 0..iterations {
            let p0 = builder.build_p0_simd_imaginary(omega).unwrap();
            std::hint::black_box(p0);
        }
        let simd_duration = start.elapsed();

        // Standard version
        let start = Instant::now();
        for _ in 0..iterations {
            let p0 = builder.build_p0_standard(omega).unwrap();
            std::hint::black_box(p0);
        }
        let standard_duration = start.elapsed();

        let speedup = standard_duration.as_secs_f64() / simd_duration.as_secs_f64();

        println!("  Iterations:     {}", iterations);
        println!("  SIMD time:      {:.3} ms", simd_duration.as_millis());
        println!("  Standard time:  {:.3} ms", standard_duration.as_millis());
        println!("  Speedup:        {:.2}x", speedup);
        println!();

        // Verify correctness
        println!("Verifying P0 Correctness:");
        println!("-------------------------");

        let p0_simd = builder.build_p0_simd_imaginary(omega).unwrap();
        let p0_standard = builder.build_p0_standard(omega).unwrap();

        let mut max_diff_re: f64 = 0.0;
        let mut max_diff_im: f64 = 0.0;

        for i in 0..naux {
            for j in 0..naux {
                let diff_re = (p0_simd[[i, j]].re - p0_standard[[i, j]].re).abs();
                let diff_im = (p0_simd[[i, j]].im - p0_standard[[i, j]].im).abs();
                max_diff_re = max_diff_re.max(diff_re);
                max_diff_im = max_diff_im.max(diff_im);
            }
        }

        if max_diff_re < 1e-12 && max_diff_im < 1e-12 {
            println!("  ✓ P0 matrices match within tolerance");
            println!("    Max diff (real): {:.6e}", max_diff_re);
            println!("    Max diff (imag): {:.6e}", max_diff_im);
        } else {
            println!("  ✗ P0 matrices differ");
            println!("    Max diff (real): {:.6e}", max_diff_re);
            println!("    Max diff (imag): {:.6e}", max_diff_im);
        }
    }

    println!("\n✅ SIMD polarizability builder test complete!");
}