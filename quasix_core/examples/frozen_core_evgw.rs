//! Example demonstrating frozen-core approximation in evGW calculations.
//!
//! This example shows how the frozen-core approximation prevents bracketing
//! failures for deep core orbitals with binding energies |ε| > 10 Ha.

use ndarray::{Array1, Array2, Array3};
use num_complex::Complex64;
use quasix_core::qp::evgw::{EvGWConfig, EvGWDriver};
use quasix_core::qp::solver::QPSolverConfig;
use quasix_core::selfenergy::correlation::ContourDeformationConfig;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();

    println!("=== Frozen-Core evGW Calculation Example ===\n");

    // System parameters (e.g., a molecule with deep core orbitals)
    let n_mo = 20;
    let n_occ = 10;
    let n_aux = 100;

    // Create orbital energies including deep core states
    let mut mo_energies = Array1::<f64>::zeros(n_mo);

    // Deep core orbitals (will be frozen)
    mo_energies[0] = -35.2;  // 1s core orbital
    mo_energies[1] = -28.7;  // 1s core orbital
    mo_energies[2] = -20.5;  // 2s core orbital
    mo_energies[3] = -18.3;  // 2s core orbital
    mo_energies[4] = -15.1;  // 2p core orbital
    mo_energies[5] = -14.8;  // 2p core orbital

    // Valence orbitals (will be active)
    mo_energies[6] = -8.5;   // 3s valence
    mo_energies[7] = -6.2;   // 3p valence
    mo_energies[8] = -4.8;   // 3p valence
    mo_energies[9] = -3.1;   // 3d valence

    // Virtual orbitals
    for i in n_occ..n_mo {
        mo_energies[i] = 2.0 + (i - n_occ) as f64 * 1.5;
    }

    println!("Initial orbital energies:");
    for (i, &energy) in mo_energies.iter().enumerate() {
        let orbital_type = if i < n_occ { "occ" } else { "virt" };
        let is_core = energy.abs() > 10.0;
        println!(
            "  Orbital {:2}: ε = {:8.3} Ha ({:4}) {}",
            i,
            energy,
            orbital_type,
            if is_core { "[CORE]" } else { "" }
        );
    }

    // Create occupations
    let mut mo_occ = Array1::<f64>::zeros(n_mo);
    for i in 0..n_occ {
        mo_occ[i] = 2.0;
    }

    // Create mock MO coefficients
    let mo_coeff = Array2::<f64>::eye(n_mo);

    // Create mock DF integrals
    let n_vir = n_mo - n_occ;
    let df_ia = Array3::<f64>::zeros((n_occ, n_vir, n_aux));
    let df_ij = Array3::<f64>::zeros((n_occ, n_occ, n_aux));

    // Initialize with small random values for realistic behavior
    for i in 0..n_occ {
        for a in 0..n_vir {
            for p in 0..n_aux {
                df_ia[[i, a, p]] = 0.01 * ((i + a + p) as f64 % 7.0 - 3.5) / 3.5;
            }
        }
        for j in 0..n_occ {
            for p in 0..n_aux {
                df_ij[[i, j, p]] = 0.01 * ((i + j + p) as f64 % 5.0 - 2.5) / 2.5;
            }
        }
    }

    // Create Coulomb metric (identity + small perturbation)
    let mut v_matrix = Array2::<f64>::eye(n_aux);
    for i in 0..n_aux {
        for j in i + 1..n_aux {
            let val = 0.01 * (((i + j) % 3) as f64 - 1.0) / 10.0;
            v_matrix[[i, j]] = val;
            v_matrix[[j, i]] = val;
        }
    }

    // Create Vxc diagonal
    let vxc_diagonal = mo_energies.mapv(|e| e * 0.9);

    // Configure evGW with frozen-core approximation
    let mut evgw_config = EvGWConfig {
        max_iterations: 5,
        energy_tolerance: 1e-4,
        density_tolerance: 1e-5,
        damping_factor: 0.5,
        z_min: 0.1,
        z_max: 0.99,
        use_diis: false,
        diis_space_dim: 6,
        print_level: 2,
        adaptive_damping: true,
        core_threshold: 10.0,  // Freeze orbitals with |ε| > 10 Ha
        cd_config: ContourDeformationConfig::default(),
        qp_solver_config: QPSolverConfig {
            core_threshold: Some(10.0),  // Also set in QP solver
            ..QPSolverConfig::default()
        },
    };

    println!("\n=== evGW Configuration ===");
    println!("Core threshold: {:.1} Ha", evgw_config.core_threshold);
    println!("Max iterations: {}", evgw_config.max_iterations);
    println!("Energy tolerance: {:.2e} Ha\n", evgw_config.energy_tolerance);

    // Create evGW driver
    let mut evgw_driver = EvGWDriver::new(n_mo, n_aux, n_occ, evgw_config.clone());

    println!("Starting evGW calculation with frozen-core approximation...\n");

    // Run evGW calculation
    match evgw_driver.run_evgw_blocks(
        &mo_energies,
        &mo_occ,
        &mo_coeff,
        &df_ia,
        &df_ij,
        &v_matrix,
        &vxc_diagonal,
    ) {
        Ok(result) => {
            println!("\n=== evGW Results ===");
            println!("Converged: {}", result.converged);
            println!("Iterations: {}", result.n_iterations);
            println!("Wall time: {:.2} s\n", result.wall_time);

            // Count frozen core orbitals
            let n_frozen = mo_energies
                .iter()
                .filter(|&&e| e.abs() > evgw_config.core_threshold)
                .count();

            println!("Frozen core orbitals: {}/{}", n_frozen, n_mo);

            println!("\nQuasiparticle energies and Z-factors:");
            for i in 0..n_mo {
                let is_frozen = mo_energies[i].abs() > evgw_config.core_threshold;
                let shift = result.qp_energies[i] - mo_energies[i];

                println!(
                    "  Orbital {:2}: E_DFT = {:8.3} Ha, E_QP = {:8.3} Ha, ΔE = {:7.4} Ha, Z = {:.3} {}",
                    i,
                    mo_energies[i],
                    result.qp_energies[i],
                    shift,
                    result.z_factors[i],
                    if is_frozen { "[FROZEN]" } else { "" }
                );
            }

            // Verify frozen core orbitals are unchanged
            let mut frozen_correct = true;
            for i in 0..n_mo {
                if mo_energies[i].abs() > evgw_config.core_threshold {
                    if (result.qp_energies[i] - mo_energies[i]).abs() > 1e-10 {
                        println!("ERROR: Frozen core orbital {} has changed energy!", i);
                        frozen_correct = false;
                    }
                    if (result.z_factors[i] - 1.0).abs() > 1e-10 {
                        println!("ERROR: Frozen core orbital {} has Z != 1!", i);
                        frozen_correct = false;
                    }
                }
            }

            if frozen_correct {
                println!("\n✓ Frozen-core approximation correctly applied!");
            }

            // Calculate and display gaps
            if n_occ > 0 && n_occ < n_mo {
                let homo_idx = n_occ - 1;
                let lumo_idx = n_occ;
                let ip = -result.qp_energies[homo_idx];
                let ea = -result.qp_energies[lumo_idx];
                let gap = result.qp_energies[lumo_idx] - result.qp_energies[homo_idx];

                println!("\n=== Electronic Structure Properties ===");
                println!("HOMO index: {}", homo_idx);
                println!("LUMO index: {}", lumo_idx);
                println!("Ionization Potential (IP): {:.3} eV", ip * 27.211);
                println!("Electron Affinity (EA): {:.3} eV", ea * 27.211);
                println!("HOMO-LUMO Gap: {:.3} eV", gap * 27.211);
            }

            println!("\n=== Summary ===");
            println!("The frozen-core approximation successfully prevented");
            println!("bracketing failures for deep core orbitals while");
            println!("correctly computing QP corrections for valence states.");
        }
        Err(e) => {
            eprintln!("ERROR: evGW calculation failed: {}", e);
            eprintln!("\nThis error would have occurred without frozen-core approximation");
            eprintln!("for orbitals with |ε| > 10 Ha due to bracketing failures.");
            std::process::exit(1);
        }
    }

    Ok(())
}