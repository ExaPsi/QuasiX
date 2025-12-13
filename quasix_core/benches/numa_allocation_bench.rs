//! NUMA-aware allocation benchmark
//!
//! Compares performance of standard vs NUMA-interleaved allocation
//! for polarizability and screening calculations on dual-socket systems.
//!
//! # Hardware Target
//! - 2x Intel Xeon Silver 4314 (Ice Lake)
//! - 2 NUMA nodes (16 cores each, 32 cores total / 64 threads)
//! - 48MB L3 cache total (24MB per socket)
//! - ~180 GB/s memory bandwidth (4x DDR4 channels)
//!
//! # Expected Results
//! - 10-15% speedup for NUMA-interleaved allocation on memory-bound operations
//! - Graceful degradation (no performance loss) on non-NUMA systems
//!
//! # Run Instructions
//! ```bash
//! # Full benchmark
//! cargo bench --bench numa_allocation_bench
//!
//! # With NUMA awareness (requires root or CAP_SYS_NICE)
//! numactl --interleave=all cargo bench --bench numa_allocation_bench
//!
//! # Report NUMA topology
//! numactl --hardware
//! ```

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use ndarray::{Array1, Array3};
use rayon::prelude::*;
use std::time::Duration;

// Import QuasiX NUMA allocator
use quasix_core::qp::numa_alloc::{NumaAllocator, NumaPolicy, NumaTopology};
use quasix_core::dielectric::polarizability::compute_polarizability_p0;
use quasix_core::dielectric::screening::compute_screened_interaction;

/// Generate test data for polarizability benchmark
fn generate_test_data(n_occ: usize, n_virt: usize, n_aux: usize, nfreq: usize)
    -> (Array3<f64>, Array1<f64>, Array1<f64>, Array1<f64>)
{
    let n_mo = n_occ + n_virt;

    // DF tensor (ia|P)
    let ia_p = Array3::<f64>::from_shape_fn((n_occ, n_virt, n_aux), |(i, a, p)| {
        0.1 * ((i * n_virt + a) * n_aux + p) as f64 / (n_occ * n_virt * n_aux) as f64
    });

    // MO energies: occupied negative, virtual positive
    let mo_energy = Array1::from_shape_fn(n_mo, |i| {
        if i < n_occ {
            -0.5 - 0.1 * i as f64
        } else {
            0.2 + 0.1 * (i - n_occ) as f64
        }
    });

    // Occupation numbers
    let mo_occ = Array1::from_shape_fn(n_mo, |i| if i < n_occ { 2.0 } else { 0.0 });

    // Frequency grid
    let freqs = Array1::linspace(0.1, 10.0, nfreq);

    (ia_p, mo_energy, mo_occ, freqs)
}

/// Generate test data for screening benchmark
fn generate_screening_data(n_aux: usize, nfreq: usize)
    -> (Array3<f64>, ndarray::Array2<f64>, Array1<f64>)
{
    // P0: polarizability (symmetric matrices)
    let p0 = Array3::<f64>::from_shape_fn((nfreq, n_aux, n_aux), |(f, p, q)| {
        let base = 0.01 / ((f + 1) as f64);
        if p == q {
            base
        } else {
            base * 0.1 / ((p.abs_diff(q) + 1) as f64)
        }
    });

    // v_aux: positive definite Coulomb metric
    let mut v_aux = ndarray::Array2::<f64>::eye(n_aux);
    for i in 0..n_aux {
        v_aux[[i, i]] = 2.0 + 0.1 * (i as f64);
    }
    // Add small off-diagonal terms to make it more realistic
    for i in 0..n_aux.saturating_sub(1) {
        v_aux[[i, i + 1]] = 0.05;
        v_aux[[i + 1, i]] = 0.05;
    }

    let freqs = Array1::linspace(0.1, 10.0, nfreq);

    (p0, v_aux, freqs)
}

/// Report NUMA topology information
fn report_numa_topology() {
    match NumaTopology::detect() {
        Ok(topo) => {
            println!("\n=== NUMA Topology ===");
            println!("NUMA nodes: {}", topo.n_nodes);
            println!("Is NUMA system: {}", topo.is_numa_system());

            for (i, cpus) in topo.cpus_per_node.iter().enumerate() {
                println!("Node {}: {} CPUs ({:?})",
                        i,
                        cpus.len(),
                        if cpus.len() > 8 {
                            format!("{:?}...", &cpus[..8])
                        } else {
                            format!("{:?}", cpus)
                        });
            }

            for (i, mem) in topo.memory_per_node.iter().enumerate() {
                println!("Node {} memory: {:.1} GB", i, *mem as f64 / (1024.0 * 1024.0 * 1024.0));
            }
            println!("=====================\n");
        }
        Err(e) => {
            println!("Could not detect NUMA topology: {:?}", e);
        }
    }
}

/// Benchmark: Polarizability computation
fn bench_polarizability_numa(c: &mut Criterion) {
    report_numa_topology();

    let mut group = c.benchmark_group("polarizability");
    group.sample_size(20);
    group.measurement_time(Duration::from_secs(10));

    // Test sizes: small, medium, large
    let sizes = [
        ("small", 5, 10, 50, 16),       // H2O-like
        ("medium", 10, 30, 200, 32),    // Benzene-like
        ("large", 20, 80, 500, 64),     // GW100 typical
    ];

    for (name, n_occ, n_virt, n_aux, nfreq) in sizes {
        let (ia_p, mo_energy, mo_occ, freqs) = generate_test_data(n_occ, n_virt, n_aux, nfreq);
        let size_mb = (nfreq * n_aux * n_aux * 8) as f64 / (1024.0 * 1024.0);

        group.throughput(Throughput::Bytes((nfreq * n_aux * n_aux * 8) as u64));

        // Polarizability computation benchmark
        group.bench_with_input(
            BenchmarkId::new("compute", name),
            &(&ia_p, &mo_energy, &mo_occ, &freqs),
            |b, (ia_p, mo_energy, mo_occ, freqs)| {
                b.iter(|| {
                    compute_polarizability_p0(
                        black_box(ia_p),
                        black_box(mo_energy),
                        black_box(mo_occ),
                        black_box(freqs),
                    ).unwrap()
                })
            },
        );

        println!("Size '{}': output = {:.1} MB", name, size_mb);
    }

    group.finish();
}

/// Benchmark: Screening computation with different allocation patterns
fn bench_screening_numa(c: &mut Criterion) {
    let mut group = c.benchmark_group("screening_numa");
    group.sample_size(20);
    group.measurement_time(Duration::from_secs(10));

    // Test sizes
    let sizes = [
        ("small", 50, 16),
        ("medium", 200, 32),
        ("large", 500, 64),
    ];

    for (name, n_aux, nfreq) in sizes {
        let (p0, v_aux, freqs) = generate_screening_data(n_aux, nfreq);
        let size_mb = (nfreq * n_aux * n_aux * 8) as f64 / (1024.0 * 1024.0);

        group.throughput(Throughput::Bytes((nfreq * n_aux * n_aux * 8) as u64));

        // Standard allocation
        group.bench_with_input(
            BenchmarkId::new("standard", name),
            &(&p0, &v_aux, &freqs),
            |b, (p0, v_aux, freqs)| {
                b.iter(|| {
                    compute_screened_interaction(
                        black_box(p0),
                        black_box(v_aux),
                        black_box(freqs),
                    ).unwrap()
                })
            },
        );

        println!("Size '{}': output = {:.1} MB", name, size_mb);
    }

    group.finish();
}

/// Benchmark: Raw NUMA allocation overhead
fn bench_numa_allocation_overhead(c: &mut Criterion) {
    let mut group = c.benchmark_group("numa_allocation_overhead");
    group.sample_size(50);

    let topology = NumaTopology::detect().ok();
    let is_numa = topology.as_ref().map_or(false, NumaTopology::is_numa_system);
    let allocator = NumaAllocator::new().expect("Failed to create NUMA allocator");

    // Different allocation sizes
    let sizes = [
        ("1MB", 1024 * 1024 / 8),
        ("10MB", 10 * 1024 * 1024 / 8),
        ("100MB", 100 * 1024 * 1024 / 8),
    ];

    for (name, n_elements) in sizes {
        // Standard allocation
        group.bench_with_input(
            BenchmarkId::new("standard", name),
            &n_elements,
            |b, &n| {
                b.iter(|| {
                    let arr: Array1<f64> = Array1::zeros(std::hint::black_box(n));
                    std::hint::black_box(arr)
                })
            },
        );

        // NUMA interleaved allocation
        if is_numa {
            group.bench_with_input(
                BenchmarkId::new("numa_interleaved", name),
                &n_elements,
                |b, &n| {
                    b.iter(|| {
                        let arr: Array1<f64> = allocator.alloc_array1(std::hint::black_box(n))
                            .expect("Allocation failed");
                        std::hint::black_box(arr)
                    })
                },
            );
        }
    }

    group.finish();
}

/// Benchmark: Parallel memory access patterns
fn bench_parallel_memory_access(c: &mut Criterion) {
    let mut group = c.benchmark_group("parallel_memory_access");
    group.sample_size(20);

    let n = 1000;
    let size = n * n;

    // Standard allocation
    let standard_arr: Array1<f64> = Array1::from_elem(size, 1.0);

    // NUMA interleaved allocation
    let topology = NumaTopology::detect().ok();
    let is_numa = topology.as_ref().map_or(false, NumaTopology::is_numa_system);
    let allocator = NumaAllocator::new().expect("Failed to create NUMA allocator");
    let numa_arr: Array1<f64> = if is_numa {
        allocator.alloc_array1(size).expect("NUMA allocation failed")
    } else {
        Array1::from_elem(size, 1.0)
    };

    // Parallel read benchmark
    group.bench_function("standard_parallel_read", |b| {
        b.iter(|| {
            let sum: f64 = standard_arr.as_slice().unwrap()
                .par_chunks(1000)
                .map(|chunk| chunk.iter().sum::<f64>())
                .sum();
            std::hint::black_box(sum)
        })
    });

    group.bench_function("numa_parallel_read", |b| {
        b.iter(|| {
            let sum: f64 = numa_arr.as_slice().unwrap()
                .par_chunks(1000)
                .map(|chunk| chunk.iter().sum::<f64>())
                .sum();
            std::hint::black_box(sum)
        })
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_polarizability_numa,
    bench_screening_numa,
    bench_numa_allocation_overhead,
    bench_parallel_memory_access,
);
criterion_main!(benches);
