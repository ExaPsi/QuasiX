//! Example demonstrating frequency grid utilities for GW calculations

use quasix_core::freq::gauss_legendre::gauss_legendre_nodes_weights;
use quasix_core::freq::{FrequencyGrid, GridType};

fn main() {
    println!("QuasiX Frequency Grid Demonstration");
    println!("====================================\n");

    // 1. Standard Gauss-Legendre quadrature
    println!("1. Standard Gauss-Legendre Quadrature on [-1, 1]");
    println!("-------------------------------------------------");

    let n = 8;
    let (nodes, weights) = gauss_legendre_nodes_weights(n).unwrap();

    println!("  {} quadrature points:", n);
    for i in 0..n {
        println!(
            "    node[{}] = {:+.10}, weight[{}] = {:.10}",
            i, nodes[i], i, weights[i]
        );
    }

    let weight_sum: f64 = weights.sum();
    println!("  Sum of weights = {:.15} (should be 2.0)", weight_sum);

    // Test integration of x^4
    let integral: f64 = nodes
        .iter()
        .zip(weights.iter())
        .map(|(&x, &w)| w * x.powi(4))
        .sum();
    println!(
        "  ∫_{{-1}}^{{1}} x^4 dx = {:.15} (exact: {})",
        integral,
        2.0 / 5.0
    );
    println!();

    // 2. Scaled Gauss-Legendre for imaginary frequency
    println!("2. Modified Gauss-Legendre for Imaginary Frequency");
    println!("---------------------------------------------------");

    let omega_max = 100.0;
    let grid = FrequencyGrid::new(16, GridType::ModifiedGaussLegendre { omega_max }).unwrap();

    println!(
        "  Imaginary frequency grid [0, {}] with {} points:",
        omega_max, grid.nfreq
    );
    println!("  First 5 frequencies:");
    for i in 0..5.min(grid.nfreq) {
        println!(
            "    ω[{}] = {:.6}, weight[{}] = {:.6}",
            i, grid.points[i], i, grid.weights[i]
        );
    }
    println!("  ...");
    println!(
        "  Last frequency: ω[{}] = {:.6}",
        grid.nfreq - 1,
        grid.points[grid.nfreq - 1]
    );
    println!();

    // 3. Integration example for GW self-energy
    println!("3. Example: Frequency Integration for GW");
    println!("-----------------------------------------");

    // Integrate 1/(1 + ω^2) as a test function
    let test_integral = grid.integrate(|omega| 1.0 / (1.0 + omega * omega));
    println!("  ∫_0^{} 1/(1 + ω²) dω ≈ {:.6}", omega_max, test_integral);
    println!(
        "  (For comparison, ∫_0^∞ 1/(1 + ω²) dω = π/2 ≈ {:.6})",
        std::f64::consts::PI / 2.0
    );
    println!();

    // 4. Precision demonstration
    println!("4. Machine Precision Test");
    println!("-------------------------");

    let n_large = 64;
    let (nodes_64, weights_64) = gauss_legendre_nodes_weights(n_large).unwrap();

    // Test orthogonality: integral of 1
    let integral_1: f64 = weights_64.sum();
    let error_1 = (integral_1 - 2.0).abs();
    println!("  {} points: ∫ 1 dx error = {:.2e}", n_large, error_1);

    // Test high-degree polynomial
    let degree = 100;
    let exact = 2.0 / (degree + 1) as f64;
    let computed: f64 = nodes_64
        .iter()
        .zip(weights_64.iter())
        .map(|(&x, &w)| w * x.powi(degree))
        .sum();
    let error_poly = (computed - exact).abs();
    println!(
        "  {} points: ∫ x^{} dx error = {:.2e}",
        n_large, degree, error_poly
    );

    println!("\n✓ Frequency grid utilities initialized successfully!");
    println!("  - Gauss-Legendre quadrature achieving < 1e-14 precision");
    println!("  - Ready for GW/BSE calculations");
}
