#!/usr/bin/env rust-script
//! Test evGW performance directly

use ndarray::{Array1, Array2, Array3};
use quasix_core::freq::{FrequencyGrid, GridType};
use quasix_core::gw::evgw::{EvGWConfig, EvGWDriver};
use std::time::Instant;

fn main() {
    println!("=== evGW Performance Verification ===\n");

    // Test with different system sizes
    let test_cases = vec![
        ("H2O-like", 13, 7, 41),    // Small
        ("NH3-like", 16, 10, 60),   // Medium
        ("Benzene-like", 42, 21, 168), // Large
    ];

    for (name, nbasis, nocc, naux) in test_cases {
        println!("\n--- {} System ---", name);
        println!("  nbasis={}, nocc={}, nvirt={}, naux={}",
                 nbasis, nocc, nbasis - nocc, naux);

        // Create test data
        let mo_energies = Array1::linspace(-1.0, 1.0, nbasis);
        let mo_occ = {
            let mut occ = Array1::zeros(nbasis);
            for i in 0..nocc {
                occ[i] = 2.0;
            }
            occ
        };

        // Create mock DF tensors
        let nvirt = nbasis - nocc;
        let ia_p = Array3::from_shape_fn((nocc, nvirt, naux), |(i, a, p)| {
            0.01 * ((i + 1) * (a + 1) * (p + 1)) as f64
        });
        let chol_v = Array2::eye(naux);
        let vxc_dft = Array1::zeros(nbasis);
        let freq_grid = FrequencyGrid::new(32, GridType::GaussLegendre).unwrap();

        // Test with different thread counts
        for n_threads in [1, 2, 4, 8] {
            if n_threads > num_cpus::get() {
                continue;
            }

            let config = EvGWConfig {
                n_threads,
                parallel_freq: true,
                cache_align: true,
                block_size: 64,
                max_cycle: 2,  // Just 2 iterations for testing
                verbose: 0,
                ..Default::default()
            };

            let mut driver = EvGWDriver::new(nbasis, nocc, naux, config);

            let start = Instant::now();
            let result = driver.run_evgw_loop(
                &mo_energies,
                &mo_occ,
                &ia_p,
                &chol_v,
                &vxc_dft,
                &freq_grid,
            );
            let elapsed = start.elapsed();

            if result.is_ok() {
                println!("  {} threads: {:.3}s", n_threads, elapsed.as_secs_f64());
            } else {
                println!("  {} threads: FAILED - {:?}", n_threads, result.err());
            }
        }
    }

    println!("\n=== Performance Analysis ===");

    // Calculate theoretical complexity
    let benzene_ops = 42_usize * 21 * 168 * 32;  // N_basis * N_virt * N_aux * N_freq
    let flops_per_op = 100_f64; // Estimate
    let total_gflops = (benzene_ops as f64 * flops_per_op) / 1e9;

    println!("Theoretical complexity for benzene-like:");
    println!("  O(N⁴) operations: ~{:.1} GFLOP", total_gflops);
    println!("  Target per iteration: <1s on 8 cores");
    println!("  Memory requirement: ~{:.1} MB",
             (42 * 168 * 8 * 10) as f64 / (1024.0 * 1024.0));
}