// Test program for the correlation self-energy fix
use ndarray::{Array1, Array2, Array3};
use quasix_core::selfenergy::correlation_proper::{
    ProperCorrelationConfig, ProperCorrelationSelfEnergy,
};

fn main() {
    println!("Testing Correlation Self-Energy Fix");
    println!("====================================");

    // Set up a simple test system (like H2O)
    let n_mo = 24;
    let n_occ = 5;
    let n_vir = n_mo - n_occ;
    let n_aux = 50;

    // Create MO energies (in Hartree)
    let mut mo_energy = Array1::<f64>::zeros(n_mo);
    // Occupied orbitals (negative energies)
    mo_energy[0] = -20.0;  // Core
    mo_energy[1] = -1.3;   // Valence
    mo_energy[2] = -0.7;
    mo_energy[3] = -0.6;
    mo_energy[4] = -0.493; // HOMO (~-13.4 eV)
    // Virtual orbitals (positive energies)
    for i in n_occ..n_mo {
        mo_energy[i] = 0.185 + 0.1 * (i - n_occ) as f64; // LUMO starts at ~5 eV
    }

    // Create occupation array
    let mut mo_occ = Array1::<f64>::zeros(n_mo);
    for i in 0..n_occ {
        mo_occ[i] = 2.0;
    }

    // Create a simple DF tensor (ia|P)
    let mut df_ia = Array3::<f64>::zeros((n_occ, n_vir, n_aux));
    for i in 0..n_occ {
        for a in 0..n_vir {
            for p in 0..n_aux {
                // Simple Gaussian-like structure
                let r = ((i as f64 - p as f64 / 3.0).powi(2) +
                        (a as f64 - p as f64 / 2.0).powi(2)) / 20.0;
                df_ia[[i, a, p]] = (-r).exp() * 0.5;
            }
        }
    }

    // Create Coulomb matrix in auxiliary basis
    let mut v_aux = Array2::<f64>::zeros((n_aux, n_aux));
    for p in 0..n_aux {
        for q in 0..n_aux {
            let diff = (p as f64 - q as f64).abs();
            v_aux[[p, q]] = 2.0 / (1.0 + diff);
        }
        v_aux[[p, p]] = 4.0; // Diagonal dominance
    }

    // Configuration for correlation calculation
    let config = ProperCorrelationConfig {
        eta: 0.001,
        n_omega: 32,
        omega_max: 30.0,
        use_gl_grid: true,
        ..Default::default()
    };

    // Create correlation self-energy calculator
    let calc = ProperCorrelationSelfEnergy::new(n_mo, n_aux, n_occ, config);

    // Evaluation points (at MO energies)
    let eval_points = mo_energy.clone();

    // Compute correlation self-energy
    println!("\nComputing correlation self-energy...");
    match calc.compute_sigma_c(&mo_energy, &mo_occ, &df_ia, &v_aux, &eval_points) {
        Ok(sigma_c) => {
            println!("Success! Correlation self-energy computed.");

            // Extract HOMO and LUMO corrections
            let homo_idx = n_occ - 1;
            let lumo_idx = n_occ;

            let sigma_c_homo = sigma_c[[homo_idx, homo_idx]];
            let sigma_c_lumo = sigma_c[[lumo_idx, lumo_idx]];

            println!("\nResults:");
            println!("--------");
            println!("HOMO Σᶜ: {:.3} eV (real), {:.3} eV (imag)",
                     sigma_c_homo.re * 27.211,
                     sigma_c_homo.im * 27.211);
            println!("LUMO Σᶜ: {:.3} eV (real), {:.3} eV (imag)",
                     sigma_c_lumo.re * 27.211,
                     sigma_c_lumo.im * 27.211);

            // Check sign convention (PySCF-like)
            println!("\nSign Convention Check:");
            println!("----------------------");
            if sigma_c_homo.re > 0.0 {
                println!("✓ HOMO Σᶜ is positive (reduces binding) - CORRECT");
            } else {
                println!("✗ HOMO Σᶜ is negative - WRONG (should be positive)");
            }

            if sigma_c_lumo.re < 0.0 {
                println!("✓ LUMO Σᶜ is negative (increases binding) - CORRECT");
            } else {
                println!("✗ LUMO Σᶜ is positive - WRONG (should be negative)");
            }

            // Check magnitude
            println!("\nMagnitude Check:");
            println!("----------------");
            let homo_magnitude = sigma_c_homo.re.abs() * 27.211;
            let lumo_magnitude = sigma_c_lumo.re.abs() * 27.211;

            if homo_magnitude > 0.5 && homo_magnitude < 3.0 {
                println!("✓ HOMO |Σᶜ| = {:.2} eV is in expected range (0.5-3.0 eV)", homo_magnitude);
            } else {
                println!("✗ HOMO |Σᶜ| = {:.2} eV is outside expected range", homo_magnitude);
            }

            if lumo_magnitude > 0.1 && lumo_magnitude < 2.0 {
                println!("✓ LUMO |Σᶜ| = {:.2} eV is in expected range (0.1-2.0 eV)", lumo_magnitude);
            } else {
                println!("✗ LUMO |Σᶜ| = {:.2} eV is outside expected range", lumo_magnitude);
            }

            // Print all orbital corrections
            println!("\nAll Orbital Corrections:");
            println!("------------------------");
            for i in 0..n_mo {
                let sigma = sigma_c[[i, i]];
                let label = if i < n_occ { "occ" } else { "vir" };
                println!("Orbital {:2} ({}): Σᶜ = {:+.3} eV",
                         i, label, sigma.re * 27.211);
            }
        }
        Err(e) => {
            println!("Error computing correlation self-energy: {}", e);
        }
    }
}