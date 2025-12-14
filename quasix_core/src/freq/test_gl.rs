// Test Gauss-Legendre with known values
use super::gauss_legendre::GaussLegendre;

#[cfg(test)]
mod known_values {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_gl_n2() {
        // For n=2, the roots are ±1/√3 and weights are both 1
        let gl = GaussLegendre::new(2);

        let expected_nodes = [-1.0 / 3.0_f64.sqrt(), 1.0 / 3.0_f64.sqrt()];
        let expected_weights = [1.0, 1.0];

        for i in 0..2 {
            println!(
                "Node {}: {} (expected {})",
                i, gl.nodes[i], expected_nodes[i]
            );
            println!(
                "Weight {}: {} (expected {})",
                i, gl.weights[i], expected_weights[i]
            );
            assert_relative_eq!(gl.nodes[i], expected_nodes[i], epsilon = 1e-14);
            assert_relative_eq!(gl.weights[i], expected_weights[i], epsilon = 1e-14);
        }
    }

    #[test]
    fn test_gl_n3() {
        // For n=3, roots are -√(3/5), 0, √(3/5) and weights are 5/9, 8/9, 5/9
        let gl = GaussLegendre::new(3);

        let sqrt_3_5 = (3.0 / 5.0_f64).sqrt();
        let expected_nodes = [-sqrt_3_5, 0.0, sqrt_3_5];
        let expected_weights = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0];

        for i in 0..3 {
            println!(
                "Node {}: {} (expected {})",
                i, gl.nodes[i], expected_nodes[i]
            );
            println!(
                "Weight {}: {} (expected {})",
                i, gl.weights[i], expected_weights[i]
            );
            assert_relative_eq!(gl.nodes[i], expected_nodes[i], epsilon = 1e-14);
            assert_relative_eq!(gl.weights[i], expected_weights[i], epsilon = 1e-14);
        }
    }

    #[test]
    fn debug_weight_sum() {
        for n in [2, 3, 4, 5, 10] {
            println!("\n=== n = {} ===", n);
            let gl = GaussLegendre::new(n);

            // Print nodes and weights
            for i in 0..n {
                println!(
                    "  node[{}] = {:.15}, weight[{}] = {:.15}",
                    i, gl.nodes[i], i, gl.weights[i]
                );
            }

            // Check weight sum
            let weight_sum: f64 = gl.weights.iter().sum();
            println!("  Weight sum = {:.15} (should be 2.0)", weight_sum);
            println!("  Error = {:.2e}", (weight_sum - 2.0).abs());
        }
    }
}
