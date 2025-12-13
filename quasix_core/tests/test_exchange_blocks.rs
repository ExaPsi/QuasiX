//! Test exchange self-energy with block-structured DF tensors

use ndarray::{Array1, Array2, Array3};
use quasix_core::qp::evgw::{EvGWConfig, EvGWDriver};

#[test]
fn test_exchange_with_blocks() {
    // Create a small test system
    let n_mo = 10;
    let n_occ = 4;
    let n_vir = n_mo - n_occ;
    let n_aux = 20;

    // Create test molecular orbital energies
    let mut mo_energy = Array1::zeros(n_mo);
    for i in 0..n_occ {
        mo_energy[i] = -1.0 - (i as f64) * 0.5;  // Occupied orbitals (negative)
    }
    for i in n_occ..n_mo {
        mo_energy[i] = 0.5 + (i - n_occ) as f64 * 0.3;  // Virtual orbitals (positive)
    }

    // Create occupation array
    let mut mo_occ = Array1::zeros(n_mo);
    for i in 0..n_occ {
        mo_occ[i] = 2.0;  // Doubly occupied for RHF
    }

    // Create MO coefficients (identity for simplicity)
    let mo_coeff = Array2::eye(n_mo);

    // Create DF tensors with non-zero values
    // Occupied-virtual block (ia|P)
    let mut df_ia = Array3::zeros((n_occ, n_vir, n_aux));
    for i in 0..n_occ {
        for a in 0..n_vir {
            for p in 0..n_aux {
                df_ia[[i, a, p]] = ((i + 1) * (a + 1)) as f64 * 0.1 / (p + 1) as f64;
            }
        }
    }

    // Occupied-occupied block (ij|P) - CRUCIAL for exchange
    let mut df_ij = Array3::zeros((n_occ, n_occ, n_aux));
    for i in 0..n_occ {
        for j in 0..n_occ {
            for p in 0..n_aux {
                // Ensure non-zero values for all occupied pairs
                df_ij[[i, j, p]] = ((i + 1) + (j + 1)) as f64 * 0.2 / (p + 1) as f64;
            }
        }
    }

    // Create a positive-definite metric matrix
    let mut v_matrix = Array2::eye(n_aux);
    for i in 0..n_aux {
        for j in 0..i {
            let val = 0.05 / ((i - j) as f64 + 1.0);
            v_matrix[[i, j]] = val;
            v_matrix[[j, i]] = val;
        }
        v_matrix[[i, i]] += 2.0;  // Ensure positive definiteness
    }

    // Create Vxc diagonal (from DFT)
    let vxc_diagonal = Array1::from_vec(vec![-0.5; n_mo]);

    // Create evGW configuration
    let config = EvGWConfig {
        max_iterations: 1,  // Just one iteration to test exchange
        energy_tolerance: 1e-4,
        density_tolerance: 1e-5,
        damping_factor: 0.5,
        z_min: 0.1,
        z_max: 0.999,
        use_diis: false,
        diis_space_dim: 6,
        print_level: 2,  // Verbose output for debugging
        adaptive_damping: false,
        core_threshold: 10.0,
        cd_config: Default::default(),
        qp_solver_config: Default::default(),
    };

    // Create driver and run calculation
    let mut driver = EvGWDriver::new(n_mo, n_aux, n_occ, config);

    // Run evGW with block-structured tensors
    let result = driver.run_evgw_blocks(
        &mo_energy,
        &mo_occ,
        &mo_coeff,
        &df_ia,
        &df_ij,
        &v_matrix,
        &vxc_diagonal,
    );

    // Check that the calculation succeeded
    assert!(result.is_ok(), "evGW calculation should succeed");

    let evgw_result = result.unwrap();

    // Verify exchange self-energy has non-zero values for ALL occupied orbitals
    let sigma_x = &evgw_result.sigma_x;

    println!("Exchange self-energy diagonal elements:");
    let mut n_nonzero_occ = 0;
    for i in 0..n_mo {
        let sigma_x_ii = sigma_x[[i, i]];
        println!("  Orbital {}: Σˣ_ii = {:.6}", i, sigma_x_ii);

        if i < n_occ {
            // Occupied orbitals MUST have negative (non-zero) exchange
            assert!(
                sigma_x_ii < 0.0,
                "Occupied orbital {} should have negative exchange self-energy, got {}",
                i, sigma_x_ii
            );

            if sigma_x_ii.abs() > 1e-10 {
                n_nonzero_occ += 1;
            }
        }
    }

    // ALL occupied orbitals should have non-zero exchange
    assert_eq!(
        n_nonzero_occ, n_occ,
        "All {} occupied orbitals should have non-zero exchange self-energy, but only {} do",
        n_occ, n_nonzero_occ
    );

    // Check symmetry of exchange matrix
    for i in 0..n_mo {
        for j in 0..i {
            assert!(
                (sigma_x[[i, j]] - sigma_x[[j, i]]).abs() < 1e-10,
                "Exchange matrix should be symmetric"
            );
        }
    }

    println!("✓ Exchange self-energy computed correctly for all {} occupied orbitals", n_occ);
}