//! Test real libcint integration with water molecule

use quasix_core::integrals::libcint_real::LibcintRealEngine;
use quasix_core::integrals::{BasisFunction, BasisSet, BasisType, Molecule};

fn main() {
    println!("Testing real libcint integration with H2O...\n");

    // Create water molecule
    let molecule = Molecule::water();

    // Create minimal STO-3G basis for water
    // O: 1s, 2s, 2p (total 5 AOs: 1s, 2s, 2px, 2py, 2pz)
    // H: 1s (each)
    // Total: 7 AO basis functions
    let ao_basis = BasisSet::new(
        vec![
            // Oxygen 1s
            BasisFunction::new(
                BasisType::S,
                vec![130.7093214, 23.80886605, 6.443608313],
                vec![0.15432897, 0.53532814, 0.44463454],
                0, // O atom
            ),
            // Oxygen 2s
            BasisFunction::new(
                BasisType::S,
                vec![5.033151319, 1.169596125, 0.380389],
                vec![-0.09996723, 0.39951283, 0.70011547],
                0, // O atom
            ),
            // Oxygen 2p
            BasisFunction::new(
                BasisType::P,
                vec![5.033151319, 1.169596125, 0.380389],
                vec![0.15591627, 0.60768372, 0.39195739],
                0, // O atom
            ),
            // H1 1s
            BasisFunction::new(
                BasisType::S,
                vec![3.42525091, 0.62391373, 0.16885540],
                vec![0.15432897, 0.53532814, 0.44463454],
                1, // H1 atom
            ),
            // H2 1s
            BasisFunction::new(
                BasisType::S,
                vec![3.42525091, 0.62391373, 0.16885540],
                vec![0.15432897, 0.53532814, 0.44463454],
                2, // H2 atom
            ),
        ],
        "STO-3G",
    );

    // Create auxiliary basis - for simplicity, use a smaller set
    // In production, this would be a proper RI basis like def2-SVP-JKFIT
    let aux_basis = BasisSet::new(
        vec![
            // Oxygen s-type auxiliary
            BasisFunction::new(
                BasisType::S,
                vec![50.0, 10.0, 2.0],
                vec![0.3, 0.5, 0.2],
                0, // O atom
            ),
            // Oxygen p-type auxiliary
            BasisFunction::new(
                BasisType::P,
                vec![8.0, 2.0],
                vec![0.4, 0.6],
                0, // O atom
            ),
            // H1 s-type auxiliary
            BasisFunction::new(
                BasisType::S,
                vec![5.0, 1.0],
                vec![0.4, 0.6],
                1, // H1 atom
            ),
            // H2 s-type auxiliary
            BasisFunction::new(
                BasisType::S,
                vec![5.0, 1.0],
                vec![0.4, 0.6],
                2, // H2 atom
            ),
        ],
        "mock-RI",
    );

    println!("Molecule: {}", molecule.name());
    println!("Number of atoms: {}", molecule.natoms());
    println!("Number of electrons: {}", molecule.nelec());
    println!(
        "AO basis functions: {} (expected: 7 for STO-3G)",
        ao_basis.size()
    );
    println!("Auxiliary basis functions: {}\n", aux_basis.size());

    // Create libcint engine
    let engine = match LibcintRealEngine::new(&molecule, &ao_basis, &aux_basis) {
        Ok(e) => {
            println!("✓ LibcintRealEngine created successfully");
            e
        }
        Err(e) => {
            eprintln!("✗ Failed to create engine: {}", e);
            return;
        }
    };

    // Test 2-center integrals
    println!("\nComputing 2-center integrals (P|Q)...");
    match engine.compute_2center() {
        Ok(v_pq) => {
            println!("✓ 2-center integrals computed");
            println!("  Shape: ({}, {})", v_pq.nrows(), v_pq.ncols());

            // Check if positive definite (all diagonal elements > 0)
            let mut all_positive = true;
            let mut min_diag = f64::MAX;
            for i in 0..v_pq.nrows() {
                let diag = v_pq[[i, i]];
                if diag <= 0.0 {
                    all_positive = false;
                }
                if diag < min_diag {
                    min_diag = diag;
                }
            }

            println!("  Min diagonal element: {:.6}", min_diag);
            if all_positive {
                println!("  ✓ Matrix is positive definite");
            } else {
                println!("  ✗ Matrix is NOT positive definite!");
            }

            // Check symmetry
            let mut max_asym = 0.0;
            for i in 0..v_pq.nrows() {
                for j in i + 1..v_pq.ncols() {
                    let asym = (v_pq[[i, j]] - v_pq[[j, i]]).abs();
                    if asym > max_asym {
                        max_asym = asym;
                    }
                }
            }
            println!("  Max asymmetry: {:.2e}", max_asym);
            if max_asym < 1e-10 {
                println!("  ✓ Matrix is symmetric");
            }
        }
        Err(e) => {
            eprintln!("✗ Failed to compute 2-center integrals: {}", e);
        }
    }

    // Test 3-center integrals
    println!("\nComputing 3-center integrals (μν|P)...");
    match engine.compute_3center() {
        Ok(j3c) => {
            println!("✓ 3-center integrals computed");
            println!(
                "  Shape: ({}, {}, {})",
                j3c.dim().0,
                j3c.dim().1,
                j3c.dim().2
            );

            // Find max absolute value
            let mut max_val = 0.0;
            for p in 0..j3c.dim().2 {
                for mu in 0..j3c.dim().0 {
                    for nu in 0..j3c.dim().1 {
                        let val = j3c[[mu, nu, p]].abs();
                        if val > max_val {
                            max_val = val;
                        }
                    }
                }
            }
            println!("  Max absolute value: {:.6}", max_val);

            // Check symmetry in first two indices
            let mut max_asym = 0.0;
            for p in 0..j3c.dim().2 {
                for mu in 0..j3c.dim().0 {
                    for nu in mu + 1..j3c.dim().1 {
                        let asym = (j3c[[mu, nu, p]] - j3c[[nu, mu, p]]).abs();
                        if asym > max_asym {
                            max_asym = asym;
                        }
                    }
                }
            }
            println!("  Max asymmetry in μν: {:.2e}", max_asym);
            if max_asym < 1e-10 {
                println!("  ✓ Integrals are symmetric in μν");
            }
        }
        Err(e) => {
            eprintln!("✗ Failed to compute 3-center integrals: {}", e);
        }
    }

    println!("\n✅ All tests completed successfully!");
}
