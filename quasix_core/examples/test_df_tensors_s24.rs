//! Test program for S2.4 DF Tensor HDF5 implementation
//!
//! Demonstrates the production-ready HDF5 I/O for density-fitted tensors.

use anyhow::Result;
use ndarray::Array2;
use quasix_core::io::df_tensors::{DFMetadata, DFTensorsS24};
use std::path::Path;
use std::time::Instant;

fn main() -> Result<()> {
    println!("=== S2.4 DF Tensor HDF5 I/O Test ===\n");

    // Create test data matching H2/sto-3g system
    let naux = 28;
    let nao = 2;
    let nocc = 1;
    let nvir = 1;

    println!("Creating test DF tensors:");
    println!("  naux = {}", naux);
    println!("  nao = {}", nao);
    println!("  nocc = {}", nocc);
    println!("  nvir = {}\n", nvir);

    let df_tensors = DFTensorsS24 {
        metric_2c: Array2::eye(naux),
        cderi_3c: Array2::from_elem((naux, nao * (nao + 1) / 2), 0.1),
        ia_p: Array2::from_elem((nocc * nvir, naux), 0.01),
        ij_p: Array2::from_elem((nocc * nocc, naux), 0.02),
        metadata: DFMetadata {
            naux,
            nao,
            nocc,
            nvir,
            auxbasis: "def2-svp-jkfit".to_string(),
            basis: "sto-3g".to_string(),
            version: "1.0".to_string(),
            timestamp: Some(chrono::Utc::now().to_rfc3339()),
            git_commit: Some(env!("CARGO_PKG_VERSION").to_string()),
        },
    };

    // Test write
    let output_path = Path::new("/tmp/test_df_tensors_s24.h5");
    println!("Writing to HDF5: {}", output_path.display());

    let start = Instant::now();
    df_tensors.write_hdf5(output_path)?;
    let write_time = start.elapsed();

    // Get file size
    let file_size = std::fs::metadata(output_path)?.len();
    let data_size =
        (naux * naux + naux * (nao * (nao + 1) / 2) + (nocc * nvir) * naux + (nocc * nocc) * naux)
            * 8; // f64 = 8 bytes
    let compression_ratio = data_size as f64 / file_size as f64;

    println!("  Write time: {:.3} ms", write_time.as_secs_f64() * 1000.0);
    println!("  File size: {:.2} KB", file_size as f64 / 1024.0);
    println!(
        "  Data size: {:.2} KB (uncompressed)",
        data_size as f64 / 1024.0
    );
    println!("  Compression ratio: {:.2}x\n", compression_ratio);

    // Test read
    println!("Reading from HDF5...");
    let start = Instant::now();
    let df_loaded = DFTensorsS24::read_hdf5(output_path)?;
    let read_time = start.elapsed();

    println!("  Read time: {:.3} ms", read_time.as_secs_f64() * 1000.0);

    // Verify metadata
    println!("\nVerifying metadata:");
    assert_eq!(df_loaded.metadata.naux, naux);
    assert_eq!(df_loaded.metadata.nao, nao);
    assert_eq!(df_loaded.metadata.nocc, nocc);
    assert_eq!(df_loaded.metadata.nvir, nvir);
    assert_eq!(df_loaded.metadata.basis, "sto-3g");
    assert_eq!(df_loaded.metadata.auxbasis, "def2-svp-jkfit");
    println!("  All metadata fields match!");

    // Verify tensor shapes
    println!("\nVerifying tensor shapes:");
    assert_eq!(df_loaded.metric_2c.shape(), &[naux, naux]);
    assert_eq!(df_loaded.cderi_3c.shape(), &[naux, nao * (nao + 1) / 2]);
    assert_eq!(df_loaded.ia_p.shape(), &[nocc * nvir, naux]);
    assert_eq!(df_loaded.ij_p.shape(), &[nocc * nocc, naux]);
    println!("  All tensor shapes match!");

    // Verify values
    println!("\nVerifying tensor values:");
    let max_error_metric = compare_arrays(&df_tensors.metric_2c, &df_loaded.metric_2c);
    let max_error_cderi = compare_arrays(&df_tensors.cderi_3c, &df_loaded.cderi_3c);
    let max_error_ia = compare_arrays(&df_tensors.ia_p, &df_loaded.ia_p);
    let max_error_ij = compare_arrays(&df_tensors.ij_p, &df_loaded.ij_p);

    println!("  metric_2c max error: {:.2e}", max_error_metric);
    println!("  cderi_3c max error: {:.2e}", max_error_cderi);
    println!("  ia_p max error: {:.2e}", max_error_ia);
    println!("  ij_p max error: {:.2e}", max_error_ij);

    assert!(max_error_metric < 1e-14, "metric_2c error too large");
    assert!(max_error_cderi < 1e-14, "cderi_3c error too large");
    assert!(max_error_ia < 1e-14, "ia_p error too large");
    assert!(max_error_ij < 1e-14, "ij_p error too large");

    println!("\n=== All tests passed! ===");

    // Performance summary
    println!("\nPerformance Summary:");
    println!(
        "  Write speed: {:.1} MB/s",
        (data_size as f64 / 1024.0 / 1024.0) / write_time.as_secs_f64()
    );
    println!(
        "  Read speed: {:.1} MB/s",
        (data_size as f64 / 1024.0 / 1024.0) / read_time.as_secs_f64()
    );

    // Clean up
    std::fs::remove_file(output_path)?;

    Ok(())
}

/// Compare two arrays and return maximum absolute difference
fn compare_arrays(a: &Array2<f64>, b: &Array2<f64>) -> f64 {
    assert_eq!(a.shape(), b.shape());
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}
