//! Integration tests for improved QP solver bracketing strategies

use ndarray::{Array1, Array2};
use num_complex::Complex64 as Complex;
use quasix_core::qp::solver::{QPEquationSolver, QPSolverConfig};
use quasix_core::common::Result;

/// Test with a simple pole-like self-energy
#[test]
fn test_simple_pole_bracketing() {
    let config = QPSolverConfig {
        energy_tolerance: 1e-6,
        residual_tolerance: 1e-7,
        max_newton_iterations: 50,
        use_bisection_fallback: true,
        use_richardson: true,
        ..Default::default()
    };

    let solver = QPEquationSolver::new(config);

    // Test orbital energies
    let mo_energies = Array1::linspace(-10.0, 5.0, 10);
    let mo_occ = Array1::zeros(10);
    let vxc_diagonal = Array1::from_elem(10, -0.3);

    // Simple self-energy with pole structure
    let mo_energies_clone = mo_energies.clone();
    let sigma_func = move |orbital_idx: usize, energy: f64| -> Result<Complex> {
        let epsilon = mo_energies_clone[orbital_idx];
        // Σ(ω) = -1/(ω - ε + 0.1i)
        Ok(Complex::new(-1.0, 0.0) / Complex::new(energy - epsilon, 0.1))
    };

    let solution = solver.solve_all_orbitals(
        &mo_energies,
        &mo_occ,
        sigma_func,
        &vxc_diagonal,
    );

    assert!(solution.is_ok(), "QP solver should succeed");

    let sol = solution.unwrap();
    assert!(sol.convergence_flags.iter().filter(|&&x| x).count() >= 8,
            "At least 80% of orbitals should converge");

    println!("Simple pole test: {} of {} orbitals converged",
             sol.convergence_flags.iter().filter(|&&x| x).count(),
             mo_energies.len());
}

/// Test with more realistic GW-like self-energy
#[test]
fn test_realistic_gw_bracketing() {
    let config = QPSolverConfig {
        energy_tolerance: 1e-5,
        residual_tolerance: 1e-6,
        max_newton_iterations: 50,
        use_bisection_fallback: true,
        core_threshold: Some(20.0), // Freeze deep core
        ..Default::default()
    };

    let solver = QPEquationSolver::new(config);

    // Mix of core, valence, and virtual orbitals
    let mo_energies = Array1::from_vec(vec![
        -25.0, -12.0, -5.0, -2.0, -0.5, 0.8, 3.0, 8.0
    ]);
    let mo_occ = Array1::from_vec(vec![2.0, 2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0]);
    let vxc_diagonal = Array1::from_elem(8, -0.4);

    // More complex self-energy with multiple poles
    let mo_energies_clone = mo_energies.clone();
    let mo_occ_clone = mo_occ.clone();
    let sigma_func = move |orbital_idx: usize, energy: f64| -> Result<Complex> {
        let epsilon_n = mo_energies_clone[orbital_idx];
        let mut sigma = Complex::new(0.0, 0.0);

        // Add pole contributions (simplified GW structure)
        for i in 0..mo_occ_clone.len() {
            if mo_occ_clone[i] > 0.5 { // Occupied
                for a in 0..mo_occ_clone.len() {
                    if mo_occ_clone[a] < 0.5 { // Virtual
                        let pole_energy = mo_energies_clone[i] - mo_energies_clone[a];
                        let weight = 0.1 * (-((orbital_idx as f64 - i as f64).abs())).exp();
                        sigma += weight / Complex::new(energy - pole_energy, 0.05);
                    }
                }
            }
        }

        // Add smooth background
        sigma += Complex::new(-0.2, 0.0) / Complex::new(energy - epsilon_n + 1.0, 0.3);

        Ok(sigma)
    };

    let solution = solver.solve_all_orbitals(
        &mo_energies,
        &mo_occ,
        sigma_func,
        &vxc_diagonal,
    );

    assert!(solution.is_ok(), "QP solver should succeed for realistic case");

    let sol = solution.unwrap();
    println!("\nRealistic GW test diagnostics:");
    println!("{}", sol.diagnostics);

    // Check physical Z-factors
    for (i, z) in sol.z_factors.iter().enumerate() {
        if mo_energies[i].abs() < 20.0 { // Not frozen core
            assert!(*z > 0.0 && *z <= 1.0,
                    "Z-factor for orbital {} should be physical: {}", i, z);
        }
    }

    assert!(sol.diagnostics.converged_count >= 6,
            "At least 6 of 8 orbitals should converge");
}

/// Test pathological cases that are difficult to bracket
#[test]
fn test_pathological_cases() {
    let config = QPSolverConfig {
        energy_tolerance: 1e-4,
        residual_tolerance: 1e-5,
        max_newton_iterations: 30,
        use_bisection_fallback: true,
        ..Default::default()
    };

    let solver = QPEquationSolver::new(config);

    // Case 1: Nearly monotonic self-energy
    let epsilon = -1.0;
    let vxc = -0.3;

    // This self-energy makes bracketing very difficult
    let sigma_func_mono = |_: usize, energy: f64| -> Result<Complex> {
        // Very weak non-monotonicity
        Ok(Complex::new(-0.99 * energy - 2.0 + 0.001 * energy.sin(), 0.0))
    };

    let mo_energies = Array1::from_vec(vec![epsilon]);
    let mo_occ = Array1::zeros(1);
    let vxc_diagonal = Array1::from_vec(vec![vxc]);

    let solution = solver.solve_all_orbitals(
        &mo_energies,
        &mo_occ,
        sigma_func_mono,
        &vxc_diagonal,
    );

    // Should either succeed with linearized approximation or fail gracefully
    match solution {
        Ok(sol) => {
            println!("\nPathological case 1: Solver used {} method",
                     if sol.orbital_solutions[0].method == quasix_core::qp::solver::SolverMethod::Linearized {
                         "linearized"
                     } else {
                         "iterative"
                     });
        }
        Err(e) => {
            println!("\nPathological case 1: Expected failure: {}", e);
            assert!(e.to_string().contains("bracket") || e.to_string().contains("linearized"),
                    "Error should mention bracketing or linearization");
        }
    }

    // Case 2: Self-energy with multiple crossings
    let sigma_func_multi = |_: usize, energy: f64| -> Result<Complex> {
        Ok(Complex::new(0.5 * (3.0 * energy).sin(), 0.0))
    };

    let solution = solver.solve_all_orbitals(
        &mo_energies,
        &mo_occ,
        sigma_func_multi,
        &vxc_diagonal,
    );

    match solution {
        Ok(sol) => {
            println!("Pathological case 2: Found solution with QP energy = {:.4}",
                     sol.qp_energies[0]);
            // Should find one of the crossings
            assert!((sol.qp_energies[0] - epsilon).abs() < 2.0,
                    "QP energy should be within reasonable range");
        }
        Err(e) => {
            println!("Pathological case 2: {}", e);
        }
    }
}

/// Test frozen core approximation
#[test]
fn test_frozen_core() {
    let config = QPSolverConfig {
        core_threshold: Some(10.0), // Freeze orbitals with |ε| > 10 Ha
        ..Default::default()
    };

    let solver = QPEquationSolver::new(config);

    // Include deep core orbitals
    let mo_energies = Array1::from_vec(vec![
        -50.0, -25.0, -15.0, -8.0, -2.0, 0.5, 2.0
    ]);
    let mo_occ = Array1::from_vec(vec![2.0, 2.0, 2.0, 2.0, 2.0, 0.0, 0.0]);
    let vxc_diagonal = Array1::from_elem(7, -0.3);

    // Dummy self-energy
    let sigma_func = |_: usize, energy: f64| -> Result<Complex> {
        Ok(Complex::new(-0.1 / (energy + 1.0), 0.0))
    };

    let solution = solver.solve_all_orbitals(
        &mo_energies,
        &mo_occ,
        sigma_func,
        &vxc_diagonal,
    ).expect("Frozen core solve should succeed");

    // Check that deep core orbitals are frozen
    for i in 0..3 {
        assert!((solution.qp_energies[i] - mo_energies[i]).abs() < 1e-10,
                "Deep core orbital {} should be frozen: ε = {}, E_QP = {}",
                i, mo_energies[i], solution.qp_energies[i]);
        assert!((solution.z_factors[i] - 1.0).abs() < 1e-10,
                "Frozen core orbital {} should have Z = 1.0, got {}",
                i, solution.z_factors[i]);
    }

    // Check that valence orbitals are solved
    for i in 3..7 {
        assert!(solution.convergence_flags[i],
                "Valence orbital {} should converge", i);
    }

    println!("\nFrozen core test passed: {} frozen, {} solved",
             3, solution.diagnostics.converged_count - 3);
}

/// Test adaptive initial guess generation
#[test]
fn test_initial_guess_strategies() {
    let config = QPSolverConfig::default();
    let solver = QPEquationSolver::new(config);

    // Different orbital types
    let test_energies = vec![
        -30.0, // Deep core
        -10.0, // Core
        -3.0,  // Valence
        -0.5,  // Near Fermi
        0.5,   // Near Fermi
        3.0,   // Virtual
        10.0,  // High virtual
    ];

    for epsilon in test_energies {
        let mo_energies = Array1::from_vec(vec![epsilon]);
        let mo_occ = Array1::from_vec(vec![if epsilon < 0.0 { 2.0 } else { 0.0 }]);
        let vxc_diagonal = Array1::from_vec(vec![-0.3]);

        // Simple test self-energy
        let epsilon_copy = epsilon;
        let sigma_func = move |_: usize, energy: f64| -> Result<Complex> {
            Ok(Complex::new(-0.2 / (energy - epsilon_copy + 1.0), 0.0))
        };

        let solution = solver.solve_all_orbitals(
            &mo_energies,
            &mo_occ,
            sigma_func,
            &vxc_diagonal,
        );

        match solution {
            Ok(sol) => {
                println!("ε = {:6.2} Ha: converged with {} iterations, E_QP = {:.4}, Z = {:.3}",
                         epsilon, sol.orbital_solutions[0].iterations,
                         sol.qp_energies[0], sol.z_factors[0]);
            }
            Err(e) => {
                println!("ε = {:6.2} Ha: failed - {}", epsilon, e);
            }
        }
    }
}