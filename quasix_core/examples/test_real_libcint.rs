//! Test real libcint integration

use quasix_core::integrals::libcint_real::LibcintRealEngine;
use quasix_core::integrals::{Atom, BasisFunction, BasisSet, BasisType, Molecule};

fn main() {
    println!("Testing real libcint integration...\n");

    // Create H2 molecule - simple test case
    let molecule = Molecule::new(
        vec![
            Atom::new("H", 1, [0.0, 0.0, 0.0]),
            Atom::new("H", 1, [0.0, 0.0, 1.4]), // 1.4 bohr separation
        ],
        0, // neutral charge
        1, // singlet
    )
    .with_name("H2");

    // Create minimal STO-3G basis for H atoms (1s orbital each)
    let ao_basis = BasisSet::new(
        vec![
            // H1 - 1s orbital (STO-3G)
            BasisFunction::new(
                BasisType::S,
                vec![3.42525091, 0.62391373, 0.16885540],
                vec![0.15432897, 0.53532814, 0.44463454],
                0, // atom 0
            ),
            // H2 - 1s orbital (STO-3G)
            BasisFunction::new(
                BasisType::S,
                vec![3.42525091, 0.62391373, 0.16885540],
                vec![0.15432897, 0.53532814, 0.44463454],
                1, // atom 1
            ),
        ],
        "STO-3G",
    );

    // Create auxiliary basis - use same as AO for testing
    let aux_basis = ao_basis.clone();

    println!("Molecule: {}", molecule.name());
    println!("Number of atoms: {}", molecule.natoms());
    println!("AO basis functions: {}", ao_basis.size());
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
            println!("  V[0,0] = {:.6}", v_pq[[0, 0]]);
            println!("  V[0,1] = {:.6}", v_pq[[0, 1]]);

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
            println!("  J[0,0,0] = {:.6}", j3c[[0, 0, 0]]);
            println!("  J[0,1,0] = {:.6}", j3c[[0, 1, 0]]);

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

    println!("\nTest complete!");
}
