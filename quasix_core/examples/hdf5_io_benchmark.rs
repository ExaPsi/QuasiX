//! Benchmark HDF5 I/O optimizations for GW/BSE datasets
//!
//! This example demonstrates the performance improvements from:
//! - Optimal chunking strategies
//! - Compression with gzip/lz4
//! - Streaming writes for large arrays
//! - Zero-copy operations

use anyhow::Result;
use quasix_core::io::hdf5_io::{
    read_array_hdf5, write_dielectric_matrix, write_frequency_batch, write_transition_tensor,
};
use std::path::Path;
use std::time::Instant;

fn main() -> Result<()> {
    println!("QuasiX HDF5 I/O Optimization Benchmark");
    println!("======================================\n");

    // Test different dataset sizes
    let sizes = vec![
        ("Small", 100, 500),   // 100 transitions, 500 aux
        ("Medium", 500, 2000), // 500 transitions, 2000 aux
        ("Large", 2000, 5000), // 2000 transitions, 5000 aux
    ];

    for (label, n_trans, n_aux) in sizes {
        println!(
            "Testing {} dataset: {} transitions x {} aux functions",
            label, n_trans, n_aux
        );
        benchmark_transition_tensor(n_trans, n_aux)?;
        println!();
    }

    // Benchmark dielectric matrix I/O
    println!("Testing dielectric matrix I/O");
    benchmark_dielectric_matrix(1000)?;
    println!();

    // Benchmark frequency batch writes
    println!("Testing frequency batch writes");
    benchmark_frequency_batch(500, 32)?;

    Ok(())
}

fn benchmark_transition_tensor(n_transitions: usize, n_aux: usize) -> Result<()> {
    let path = Path::new("/tmp/quasix_benchmark_tensor.h5");

    // Generate test data
    let data: Vec<f64> = (0..n_transitions * n_aux)
        .map(|i| (i as f64).sin() * 0.1)
        .collect();

    // Benchmark write
    let start = Instant::now();
    write_transition_tensor(path, "ia_P", &data, n_transitions, n_aux)?;
    let write_time = start.elapsed();

    // Benchmark read
    let start = Instant::now();
    let _data_read = read_array_hdf5(path, "/df_tensors", "ia_P")?;
    let read_time = start.elapsed();

    // Calculate throughput
    let data_size_mb = (n_transitions * n_aux * 8) as f64 / (1024.0 * 1024.0);
    let write_throughput = data_size_mb / write_time.as_secs_f64();
    let read_throughput = data_size_mb / read_time.as_secs_f64();

    println!("  Data size: {:.2} MB", data_size_mb);
    println!(
        "  Write time: {:.3} s ({:.1} MB/s)",
        write_time.as_secs_f64(),
        write_throughput
    );
    println!(
        "  Read time:  {:.3} s ({:.1} MB/s)",
        read_time.as_secs_f64(),
        read_throughput
    );

    // Get file size to show compression ratio
    let metadata = std::fs::metadata(path)?;
    let file_size_mb = metadata.len() as f64 / (1024.0 * 1024.0);
    let compression_ratio = data_size_mb / file_size_mb;
    println!(
        "  File size:  {:.2} MB (compression ratio: {:.2}x)",
        file_size_mb, compression_ratio
    );

    // Clean up
    std::fs::remove_file(path).ok();

    Ok(())
}

fn benchmark_dielectric_matrix(n_aux: usize) -> Result<()> {
    let path = Path::new("/tmp/quasix_benchmark_dielectric.h5");

    // Generate test matrix (symmetric)
    let matrix: Vec<f64> = (0..n_aux * n_aux)
        .map(|i| {
            let row = i / n_aux;
            let col = i % n_aux;
            if row <= col {
                (i as f64).cos() * 0.1
            } else {
                // Make it symmetric
                let idx = col * n_aux + row;
                (idx as f64).cos() * 0.1
            }
        })
        .collect();

    let omega = 0.5;

    // Benchmark write
    let start = Instant::now();
    write_dielectric_matrix(path, 0, omega, &matrix, n_aux)?;
    let write_time = start.elapsed();

    // Benchmark read
    let start = Instant::now();
    let _matrix_read = read_array_hdf5(path, "/dielectric", "W_0000")?;
    let read_time = start.elapsed();

    let data_size_mb = (n_aux * n_aux * 8) as f64 / (1024.0 * 1024.0);
    println!(
        "  Matrix size: {}x{} ({:.2} MB)",
        n_aux, n_aux, data_size_mb
    );
    println!("  Write time: {:.3} s", write_time.as_secs_f64());
    println!("  Read time:  {:.3} s", read_time.as_secs_f64());

    // Clean up
    std::fs::remove_file(path).ok();

    Ok(())
}

fn benchmark_frequency_batch(n_aux: usize, n_frequencies: usize) -> Result<()> {
    let path = Path::new("/tmp/quasix_benchmark_batch.h5");

    // Generate test data for multiple frequencies
    let omega_values: Vec<f64> = (0..n_frequencies).map(|i| i as f64 * 0.1).collect();

    let matrices: Vec<Vec<f64>> = omega_values
        .iter()
        .map(|&omega| {
            (0..n_aux * n_aux)
                .map(|i| (i as f64 * omega).sin() * 0.1)
                .collect()
        })
        .collect();

    // Benchmark batch write
    let start = Instant::now();
    write_frequency_batch(path, &omega_values, &matrices, n_aux)?;
    let batch_time = start.elapsed();

    let total_size_mb = (n_aux * n_aux * n_frequencies * 8) as f64 / (1024.0 * 1024.0);
    let throughput = total_size_mb / batch_time.as_secs_f64();

    println!(
        "  Batch size: {} frequencies x {}x{} matrices",
        n_frequencies, n_aux, n_aux
    );
    println!("  Total data: {:.2} MB", total_size_mb);
    println!(
        "  Batch write time: {:.3} s ({:.1} MB/s)",
        batch_time.as_secs_f64(),
        throughput
    );

    // Clean up
    std::fs::remove_file(path).ok();

    Ok(())
}
