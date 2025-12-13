// Test the proper correlation implementation directly in Rust
use ndarray::{Array1, Array2, Array3};
use quasix_core::gw::evgw_proper::run_g0w0;
use quasix_core::selfenergy::correlation_proper::compute_sigma_c_diagonal_approximation;

fn main() {
    println!("Testing Proper Correlation Implementation in Rust");
    println!("{}", "=".repeat(60));

    // Create mock data for H2O
    let n_mo = 7;
    let n_occ = 5;
    let n_aux = 20;

    // Mock DF tensor
    let df_tensor = Array3::<f64>::zeros((n_mo, n_mo, n_aux));

    // Mock MO energies (Hartree)
    let mo_energy = Array1::from_vec(vec![
        -20.5, -1.3, -0.7, -0.55, -0.45,  // Occupied
        0.15, 0.20,  // Virtual
    ]);

    // Mock MO coefficients
    let mo_coeff = Array2::<f64>::eye(n_mo);

    // Mock v_sqrt matrix
    let v_sqrt = Array2::<f64>::eye(n_aux);

    println!("Test configuration:");
    println!("  n_mo = {}", n_mo);
    println!("  n_occ = {}", n_occ);
    println!("  n_aux = {}", n_aux);

    // Test diagonal approximation
    println!("\nTesting compute_sigma_c_diagonal_approximation...");
    let omega_eval = mo_energy.clone();
    let sigma_c_diag = compute_sigma_c_diagonal_approximation(
        &df_tensor,
        &mo_energy,
        n_occ,
        &omega_eval,
    );

    println!("Correlation self-energy diagonal (eV):");
    for (i, &sigma) in sigma_c_diag.iter().enumerate() {
        let label = if i < n_occ { "occ" } else { "vir" };
        println!("  MO {} ({}): Σc = {:.3} eV", i, label, sigma * 27.211);
    }

    // Test G0W0
    println!("\nTesting run_g0w0...");
    match run_g0w0(&df_tensor, &mo_energy, &mo_coeff, &v_sqrt, n_occ) {
        Ok(result) => {
            println!("✓ G0W0 calculation completed");
            println!("  Converged: {}", result.converged);
            println!("  Iterations: {}", result.n_iterations);

            println!("\nQuasiparticle energies (eV):");
            for (i, (&e_dft, &e_qp)) in mo_energy.iter().zip(result.qp_energies.iter()).enumerate() {
                let delta = e_qp - e_dft;
                let label = if i < n_occ { "occ" } else { "vir" };
                println!("  MO {} ({}): DFT = {:.3}, QP = {:.3}, Δ = {:.3}",
                         i, label, e_dft * 27.211, e_qp * 27.211, delta * 27.211);
            }

            // Check HOMO IP
            let homo_idx = n_occ - 1;
            let ip_dft = -mo_energy[homo_idx] * 27.211;
            let ip_qp = -result.qp_energies[homo_idx] * 27.211;
            println!("\nIonization Potential:");
            println!("  DFT: {:.2} eV", ip_dft);
            println!("  GW:  {:.2} eV", ip_qp);

            if 10.0 < ip_qp && ip_qp < 15.0 {
                println!("✓ IP is in reasonable range for H2O (expected ~12.6 eV)");
            } else {
                println!("⚠ IP = {:.2} eV is outside expected range (10-15 eV)", ip_qp);
            }
        }
        Err(e) => {
            println!("✗ G0W0 calculation failed: {:?}", e);
        }
    }
}