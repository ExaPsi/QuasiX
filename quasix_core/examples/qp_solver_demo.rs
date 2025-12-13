//! Demonstration of the S3-6 Quasiparticle Solver
//!
//! This example shows how to use the Newton-Raphson solver with bisection fallback
//! to solve the quasiparticle equation: E_n = ε_n + Re[Σ_nn(E_n) - v_xc_nn]

use ndarray::array;
use num_complex::Complex64 as Complex;
use quasix_core::common::Result;
use quasix_core::qp::solver::{QPEquationSolver, QPSolverConfig};

/// Simple model self-energy for demonstration
/// Σ(ω) = Σ_x + α·ω/(ω² + β²)
fn model_self_energy(_orbital: usize, energy: f64) -> Result<Complex> {
    let sigma_x = -0.4; // Exchange part (static)
    let alpha = 0.1; // Correlation strength
    let beta = 1.0; // Energy scale

    // Frequency-dependent correlation
    let sigma_c = alpha * energy / (energy * energy + beta * beta);

    Ok(Complex::new(sigma_x + sigma_c, 0.0))
}

fn main() -> Result<()> {
    println!("=== S3-6 Quasiparticle Solver Demo ===\n");

    // Setup: 4 molecular orbitals (2 occupied, 2 virtual)
    let mo_energies = array![-0.6, -0.3, 0.4, 0.8]; // Hartree
    let mo_occ = array![2.0, 2.0, 0.0, 0.0]; // Occupations
    let vxc_diagonal = array![-0.35, -0.25, -0.15, -0.10]; // DFT XC

    println!("System Setup:");
    println!("  MO energies (Ha): {:?}", mo_energies);
    println!("  Occupations:      {:?}", mo_occ);
    println!("  V_xc diagonal:    {:?}", vxc_diagonal);
    println!();

    // Configure solver
    let config = QPSolverConfig {
        energy_tolerance: 1e-6,
        residual_tolerance: 1e-7,
        max_newton_iterations: 30,
        max_bisection_iterations: 50,
        use_richardson: true,  // Use Richardson extrapolation
        use_line_search: true, // Enable line search
        use_bisection_fallback: true,
        n_threads: Some(2), // Use 2 threads
        ..Default::default()
    };

    println!("Solver Configuration:");
    println!("  Method: Newton-Raphson with bisection fallback");
    println!("  Richardson extrapolation: {}", config.use_richardson);
    println!("  Line search: {}", config.use_line_search);
    println!("  Tolerance: {:.1e} Ha", config.energy_tolerance);
    println!();

    // Create solver
    let solver = QPEquationSolver::new(config);

    // Solve QP equations
    println!("Solving quasiparticle equations...\n");

    let solution =
        solver.solve_all_orbitals(&mo_energies, &mo_occ, model_self_energy, &vxc_diagonal)?;

    // Display results
    println!("Results:");
    println!("┌─────────┬──────────┬──────────┬──────────┬─────────┬────────────┬──────────┐");
    println!("│ Orbital │ ε_n (Ha) │ E_QP (Ha)│ Shift(eV)│ Z-factor│ Iterations │  Method  │");
    println!("├─────────┼──────────┼──────────┼──────────┼─────────┼────────────┼──────────┤");

    for sol in &solution.orbital_solutions {
        let shift_ev = sol.energy_shift * 27.2114; // Convert to eV
        let method = if sol.used_bisection {
            "Bisection"
        } else {
            "Newton-R"
        };

        println!(
            "│    {}    │  {:7.4} │  {:7.4} │  {:7.3} │  {:6.3} │     {:2}     │ {:8} │",
            sol.orbital_idx,
            mo_energies[sol.orbital_idx],
            sol.qp_energy,
            shift_ev,
            sol.z_factor,
            sol.iterations,
            method
        );
    }
    println!("└─────────┴──────────┴──────────┴──────────┴─────────┴────────────┴──────────┘");
    println!();

    // Compute and display gap
    let homo_idx = 1; // Last occupied
    let lumo_idx = 2; // First virtual

    let ks_gap = mo_energies[lumo_idx] - mo_energies[homo_idx];
    let qp_gap = solution.qp_energies[lumo_idx] - solution.qp_energies[homo_idx];
    let gap_correction = (qp_gap - ks_gap) * 27.2114; // Convert to eV

    println!("Band Gap Analysis:");
    println!("  Kohn-Sham gap:     {:.3} eV", ks_gap * 27.2114);
    println!("  Quasiparticle gap: {:.3} eV", qp_gap * 27.2114);
    println!("  GW correction:     {:+.3} eV", gap_correction);
    println!();

    // Display diagnostics
    println!("{}", solution.diagnostics);
    println!();

    println!("Wall time: {:.3} seconds", solution.wall_time);

    // Validate physical constraints
    println!("\nPhysical Validation:");

    let all_z_physical = solution.z_factors.iter().all(|&z| z > 0.0 && z < 1.0);
    println!(
        "  All Z-factors in (0,1): {}",
        if all_z_physical { "✓" } else { "✗" }
    );

    let all_converged = solution.convergence_flags.iter().all(|&c| c);
    println!(
        "  All orbitals converged: {}",
        if all_converged { "✓" } else { "✗" }
    );

    if solution.diagnostics.problematic_z_count > 0 {
        println!(
            "  ⚠ Warning: {} orbitals have Z-factors near boundaries",
            solution.diagnostics.problematic_z_count
        );
    }

    Ok(())
}
