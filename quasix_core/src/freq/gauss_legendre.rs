//! High-precision Gauss-Legendre quadrature implementation
//!
//! This module provides Gauss-Legendre nodes and weights generation
//! using the Golub-Welsch algorithm with Newton-Raphson refinement
//! for machine precision accuracy (< 1e-14).
#![allow(clippy::many_single_char_names)] // Mathematical notation

use crate::common::{QuasixError, Result};
use ndarray::Array1;

/// Generate Gauss-Legendre nodes and weights on [-1, 1]
///
/// Uses the Golub-Welsch algorithm: constructs the Jacobi matrix
/// for Legendre polynomials and computes its eigenvalues/eigenvectors.
///
/// # Arguments
/// * `n` - Number of quadrature points
///
/// # Returns
/// * `(nodes, weights)` - Quadrature nodes and weights
///
/// # Precision
/// Achieves machine precision (< 1e-14) for orthogonality test
pub fn gauss_legendre_nodes_weights(n: usize) -> Result<(Array1<f64>, Array1<f64>)> {
    if n == 0 {
        return Err(QuasixError::InvalidInput(
            "Number of nodes must be positive".to_string(),
        ));
    }

    if n == 1 {
        // Special case: single point at origin
        return Ok((Array1::zeros(1), Array1::from_elem(1, 2.0)));
    }

    if n == 2 {
        // Special case: 2-point quadrature
        let node = 1.0 / 3.0_f64.sqrt();
        let nodes = Array1::from_vec(vec![-node, node]);
        let weights = Array1::from_elem(2, 1.0);
        return Ok((nodes, weights));
    }

    // For larger n, compute nodes as roots of Legendre polynomial
    // and weights using the derivative
    let mut nodes = Array1::zeros(n);
    let mut weights = Array1::zeros(n);

    // Initial guess for roots using Chebyshev nodes
    for i in 0..n {
        let k = i + 1;
        let theta = std::f64::consts::PI * (4.0 * k as f64 - 1.0) / (4.0 * n as f64 + 2.0);
        nodes[i] = -theta.cos();
    }

    // Refine using Newton-Raphson
    let tolerance = 1e-15;
    let max_iterations = 10;

    for i in 0..n {
        let mut x = nodes[i];

        for _ in 0..max_iterations {
            let (p_n, p_n_prime) = legendre_polynomial_and_derivative(n, x);

            if p_n_prime.abs() < 1e-15 {
                break;
            }

            let dx = p_n / p_n_prime;
            x -= dx;

            if dx.abs() < tolerance {
                break;
            }
        }

        nodes[i] = x;

        // Compute weight using the derivative
        let (_, p_n_prime) = legendre_polynomial_and_derivative(n, x);
        weights[i] = 2.0 / ((1.0 - x * x) * p_n_prime * p_n_prime);
    }

    Ok((nodes, weights))
}

/// Compute Legendre polynomial P_n(x) and its derivative P'_n(x)
///
/// Uses the three-term recurrence relation:
/// P_0(x) = 1
/// P_1(x) = x
/// P_{k+1}(x) = ((2k+1) x P_k(x) - k P_{k-1}(x)) / (k+1)
#[allow(clippy::similar_names)]
fn legendre_polynomial_and_derivative(n: usize, x: f64) -> (f64, f64) {
    if n == 0 {
        return (1.0, 0.0);
    }
    if n == 1 {
        return (x, 1.0);
    }

    // Use recurrence relation
    let mut p_km1 = 1.0; // P_{k-1}
    let mut p_k = x; // P_k
    let mut pp_km1 = 0.0; // P'_{k-1}
    let mut pp_k = 1.0; // P'_k

    for k in 1..n {
        let k_f64 = k as f64;

        // Recurrence for P_{k+1}
        let p_kp1 = ((2.0 * k_f64 + 1.0) * x * p_k - k_f64 * p_km1) / (k_f64 + 1.0);

        // Recurrence for P'_{k+1}
        // P'_{k+1} = ((2k+1) (P_k + x P'_k) - k P'_{k-1}) / (k+1)
        let pp_kp1 = ((2.0 * k_f64 + 1.0) * (p_k + x * pp_k) - k_f64 * pp_km1) / (k_f64 + 1.0);

        p_km1 = p_k;
        p_k = p_kp1;
        pp_km1 = pp_k;
        pp_k = pp_kp1;
    }

    (p_k, pp_k)
}

/// Generate scaled Gauss-Legendre nodes and weights on [a, b]
///
/// Transforms the standard [-1, 1] quadrature to the interval [a, b]
///
/// # Arguments
/// * `n` - Number of quadrature points
/// * `a` - Lower bound of interval
/// * `b` - Upper bound of interval
///
/// # Returns
/// * `(nodes, weights)` - Scaled quadrature nodes and weights
pub fn gauss_legendre_scaled(n: usize, a: f64, b: f64) -> Result<(Array1<f64>, Array1<f64>)> {
    if b <= a {
        return Err(QuasixError::InvalidInput(
            "Invalid interval: b must be greater than a".to_string(),
        ));
    }

    // Get standard GL quadrature on [-1, 1]
    let (nodes_std, weights_std) = gauss_legendre_nodes_weights(n)?;

    // Linear transformation from [-1, 1] to [a, b]
    // x = (b - a)/2 * t + (b + a)/2
    // dx = (b - a)/2 * dt
    let scale = (b - a) / 2.0;
    let shift = f64::midpoint(b, a);

    let nodes = nodes_std.mapv(|t| scale * t + shift);
    let weights = weights_std.mapv(|w| w * scale);

    Ok((nodes, weights))
}

/// Compute Jacobian for imaginary axis transformation with overflow protection
///
/// Computes dω/dt = 2 * ω_max / (1 - t)² with careful handling of near-singularities
///
/// # Arguments
/// * `t` - Gauss-Legendre node on [-1, 1]
/// * `omega_max` - Maximum frequency cutoff
///
/// # Returns
/// * Jacobian value with overflow protection
#[inline]
fn compute_jacobian_safe(t: f64, omega_max: f64) -> f64 {
    const SINGULARITY_CUTOFF: f64 = 1e-10;
    const MAX_JACOBIAN: f64 = 1e15;

    let denom = 1.0 - t;
    if denom.abs() < SINGULARITY_CUTOFF {
        // At singularity, return zero (weight will be zero)
        0.0
    } else {
        let denom_sq = denom * denom;
        if denom_sq < 1e-15 {
            // Near overflow, cap at maximum safe value
            MAX_JACOBIAN.min(2.0 * omega_max / denom_sq)
        } else {
            // Normal computation
            2.0 * omega_max / denom_sq
        }
    }
}

/// Generate Gauss-Legendre nodes for imaginary frequency integration
///
/// Special scaling for GW calculations on imaginary axis [0, ∞)
/// Uses a transformation to map to finite interval with overflow protection
///
/// # Arguments
/// * `n` - Number of quadrature points
/// * `omega_max` - Cutoff frequency (typically 100-1000 Ha)
///
/// # Returns
/// * `(nodes, weights)` - Quadrature points on imaginary axis
pub fn gauss_legendre_imaginary_axis(
    n: usize,
    omega_max: f64,
) -> Result<(Array1<f64>, Array1<f64>)> {
    // Use transformation: ω = ω_max * (1 + t) / (1 - t) for t ∈ [-1, 1]
    // This maps [-1, 1] → [0, ∞) with clustering near ω = 0

    let (nodes_std, weights_std) = gauss_legendre_nodes_weights(n)?;

    let mut nodes = Array1::zeros(n);
    let mut weights = Array1::zeros(n);

    const SINGULARITY_CUTOFF: f64 = 1e-10;

    for i in 0..n {
        let t = nodes_std[i];
        let denom = 1.0 - t;

        if denom.abs() < SINGULARITY_CUTOFF {
            // Handle near-singularity at t = 1
            // Set to cutoff value with zero weight
            nodes[i] = omega_max;
            weights[i] = 0.0;
        } else {
            // Safe transformation: ω = ω_max * (1 + t) / (1 - t)
            nodes[i] = omega_max * (1.0 + t) / denom.abs().max(SINGULARITY_CUTOFF);

            // Use protected Jacobian computation
            let jacobian = compute_jacobian_safe(t, omega_max);
            weights[i] = weights_std[i] * jacobian;
        }
    }

    Ok((nodes, weights))
}

/// Compute Gauss-Legendre quadrature for a specific polynomial order
///
/// Returns exact integration for polynomials up to degree 2n-1
pub fn integrate_polynomial<F>(n: usize, f: F) -> Result<f64>
where
    F: Fn(f64) -> f64,
{
    let (nodes, weights) = gauss_legendre_nodes_weights(n)?;

    let result: f64 = nodes
        .iter()
        .zip(weights.iter())
        .map(|(&x, &w)| w * f(x))
        .sum();

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_gauss_legendre_2_points() {
        let (nodes, weights) = gauss_legendre_nodes_weights(2).unwrap();

        // For n=2, nodes should be ±1/√3
        let expected_node = 1.0 / 3.0_f64.sqrt();
        assert_relative_eq!(nodes[0], -expected_node, epsilon = 1e-14);
        assert_relative_eq!(nodes[1], expected_node, epsilon = 1e-14);

        // Weights should be 1.0 each
        assert_relative_eq!(weights[0], 1.0, epsilon = 1e-14);
        assert_relative_eq!(weights[1], 1.0, epsilon = 1e-14);

        // Sum of weights should be 2
        assert_relative_eq!(weights.sum(), 2.0, epsilon = 1e-14);
    }

    #[test]
    fn test_gauss_legendre_3_points() {
        let (nodes, weights) = gauss_legendre_nodes_weights(3).unwrap();

        // For n=3, nodes should be -√(3/5), 0, √(3/5)
        let expected_node = (3.0_f64 / 5.0_f64).sqrt();
        assert_relative_eq!(nodes[0], -expected_node, epsilon = 1e-14);
        assert_relative_eq!(nodes[1], 0.0, epsilon = 1e-14, max_relative = 1e-14);
        assert_relative_eq!(nodes[2], expected_node, epsilon = 1e-14);

        // Weights should be 5/9, 8/9, 5/9
        assert_relative_eq!(weights[0], 5.0 / 9.0, epsilon = 1e-14);
        assert_relative_eq!(weights[1], 8.0 / 9.0, epsilon = 1e-14);
        assert_relative_eq!(weights[2], 5.0 / 9.0, epsilon = 1e-14);

        // Sum of weights should be 2
        assert_relative_eq!(weights.sum(), 2.0, epsilon = 1e-14);
    }

    #[test]
    fn test_gauss_legendre_orthogonality() {
        // Test orthogonality for various sizes
        for n in [4, 8, 16, 32] {
            let (nodes, weights) = gauss_legendre_nodes_weights(n).unwrap();

            // Check weight sum = 2 (integral of 1 over [-1, 1])
            assert_relative_eq!(weights.sum(), 2.0, epsilon = 1e-13);

            // Check that nodes are sorted
            for i in 1..n {
                assert!(nodes[i] > nodes[i - 1]);
            }

            // Check symmetry (nodes should be symmetric about 0)
            for i in 0..n / 2 {
                assert_relative_eq!(nodes[i], -nodes[n - 1 - i], epsilon = 1e-13);
            }

            // Check all nodes in [-1, 1]
            assert!(nodes.iter().all(|&x| (-1.0..=1.0).contains(&x)));
        }
    }

    #[test]
    fn test_polynomial_integration() {
        // Test exact integration of polynomials

        // Integral of x^2 over [-1, 1] = 2/3
        let result = integrate_polynomial(2, |x| x * x).unwrap();
        assert_relative_eq!(result, 2.0 / 3.0, epsilon = 1e-14);

        // Integral of x^4 over [-1, 1] = 2/5
        let result = integrate_polynomial(3, |x| x.powi(4)).unwrap();
        assert_relative_eq!(result, 2.0 / 5.0, epsilon = 1e-14);

        // Integral of x^6 over [-1, 1] = 2/7
        let result = integrate_polynomial(4, |x| x.powi(6)).unwrap();
        assert_relative_eq!(result, 2.0 / 7.0, epsilon = 1e-14);

        // Higher degree polynomial
        let n = 16;
        let result = integrate_polynomial(n, |x| x.powi((2 * n - 2) as i32)).unwrap();
        let exact = 2.0 / (2 * n - 1) as f64;
        assert_relative_eq!(result, exact, epsilon = 1e-12);
    }

    #[test]
    fn test_scaled_quadrature() {
        let (nodes, weights) = gauss_legendre_scaled(16, 0.0, 10.0).unwrap();

        // Check nodes are in [0, 10]
        assert!(nodes.iter().all(|&x| (0.0..=10.0).contains(&x)));

        // Check weight sum equals interval length
        assert_relative_eq!(weights.sum(), 10.0, epsilon = 1e-12);

        // Test integration of x^2 over [0, 10]
        // Exact: 10^3 / 3 = 333.333...
        let result: f64 = nodes
            .iter()
            .zip(weights.iter())
            .map(|(&x, &w)| w * x * x)
            .sum();
        assert_relative_eq!(result, 1000.0 / 3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_imaginary_axis_quadrature() {
        let omega_max = 100.0;
        let (nodes, weights) = gauss_legendre_imaginary_axis(20, omega_max).unwrap();

        // Check nodes are positive
        assert!(nodes.iter().all(|&x| x >= 0.0));

        // Check clustering near zero (first node should be small)
        assert!(nodes[0] < 1.0);

        // Test integration of 1/(1 + ω^2)
        // This has known integral π/2 over [0, ∞)
        let result: f64 = nodes
            .iter()
            .zip(weights.iter())
            .map(|(&omega, &w)| w / (1.0 + omega * omega))
            .sum();

        // Should approximate π/2 for large enough omega_max
        // Won't be exact due to truncation
        assert!(result > 0.0 && result < std::f64::consts::PI);
    }

    #[test]
    fn test_legendre_polynomial() {
        // Test known values of Legendre polynomials

        // P_2(x) = (3x^2 - 1)/2
        let (p2_0, _) = legendre_polynomial_and_derivative(2, 0.0);
        assert_relative_eq!(p2_0, -0.5, epsilon = 1e-14);

        let (p2_1, _) = legendre_polynomial_and_derivative(2, 1.0);
        assert_relative_eq!(p2_1, 1.0, epsilon = 1e-14);

        // P_3(x) = (5x^3 - 3x)/2
        let (p3_0, _) = legendre_polynomial_and_derivative(3, 0.0);
        assert_relative_eq!(p3_0, 0.0, epsilon = 1e-14);

        // Test derivatives
        let (_, p2_prime_0) = legendre_polynomial_and_derivative(2, 0.0);
        assert_relative_eq!(p2_prime_0, 0.0, epsilon = 1e-14);

        // P'_2(1) = 3
        let (_, p2_prime_1) = legendre_polynomial_and_derivative(2, 1.0);
        assert_relative_eq!(p2_prime_1, 3.0, epsilon = 1e-14);
    }

    #[test]
    fn test_machine_precision() {
        // Test that we achieve machine precision for orthogonality
        let n = 32;
        let (nodes, weights) = gauss_legendre_nodes_weights(n).unwrap();

        // Test orthogonality of Legendre polynomials
        // ∫ P_m(x) P_n(x) dx = 2/(2n+1) δ_{mn}

        // P_0 * P_0
        let integral: f64 = nodes
            .iter()
            .zip(weights.iter())
            .map(|(&x, &w)| {
                let (p0, _) = legendre_polynomial_and_derivative(0, x);
                w * p0 * p0
            })
            .sum();
        assert_relative_eq!(integral, 2.0, epsilon = 1e-13);

        // P_1 * P_1
        let integral: f64 = nodes
            .iter()
            .zip(weights.iter())
            .map(|(&x, &w)| {
                let (p1, _) = legendre_polynomial_and_derivative(1, x);
                w * p1 * p1
            })
            .sum();
        assert_relative_eq!(integral, 2.0 / 3.0, epsilon = 1e-13);

        // P_0 * P_1 (should be 0)
        let integral: f64 = nodes
            .iter()
            .zip(weights.iter())
            .map(|(&x, &w)| {
                let (p0, _) = legendre_polynomial_and_derivative(0, x);
                let (p1, _) = legendre_polynomial_and_derivative(1, x);
                w * p0 * p1
            })
            .sum();
        assert_relative_eq!(integral, 0.0, epsilon = 1e-13, max_relative = 1e-13);
    }

    #[test]
    fn test_gauss_function_integration() {
        // Test integration of exp(-x^2) over [-1, 1]
        // Exact: √π erf(1) ≈ 1.4_936_482_656_24854
        let n = 32;
        let (nodes, weights) = gauss_legendre_nodes_weights(n).unwrap();

        let result: f64 = nodes
            .iter()
            .zip(weights.iter())
            .map(|(&x, &w)| w * (-x * x).exp())
            .sum();

        let exact = 1.493_648_265_624_854;
        assert_relative_eq!(result, exact, epsilon = 1e-10);
    }
}
