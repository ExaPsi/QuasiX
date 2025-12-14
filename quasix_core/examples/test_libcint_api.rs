//! Test program to understand libcint API

use libcint::prelude::*;

fn main() {
    println!("Testing libcint API...");

    // Use the predefined H2O molecule with def2-TZVP basis
    let cint_data = init_h2o_def2_tzvp();

    // Test 2-center integrals (overlap)
    println!("Computing 1e overlap integrals...");
    let (overlap, shape) = cint_data.integrate("int1e_ovlp", None, None).into();
    println!("Overlap shape: {:?}", shape);
    println!(
        "First few overlap values: {:?}",
        &overlap[..5.min(overlap.len())]
    );

    // Test 2-center Coulomb integrals - these are the (P|Q) integrals we need
    println!("\nComputing 2-center Coulomb integrals (P|Q)...");
    let (_coulomb2c, shape2c) = cint_data.integrate("int2c2e", None, None).into();
    println!("2-center Coulomb shape: {:?}", shape2c);

    // Test 3-center integrals - these are the (μν|P) integrals we need
    println!("\nComputing 3-center 2e integrals (μν|P)...");
    let (_int3c, shape3c) = cint_data.integrate("int3c2e", None, None).into();
    println!("3-center shape: {:?}", shape3c);

    println!("\nAPI exploration complete");
}
