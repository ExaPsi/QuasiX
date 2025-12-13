//! Comprehensive performance profiling for DF tensor operations
//!
//! This benchmark identifies bottlenecks in:
//! - SIMD utilization (detect dead code)
//! - Thread coordination (BLAS vs Rayon)
//! - Lock contention
//! - Memory access patterns
//!
//! Run with: cargo bench --bench df_performance_profile

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, black_box};
use ndarray::{Array2, Array3};
use quasix_core::df::mo_transform::{generate_mock_mo_coefficients, TransformConfig};
use quasix_core::df::parallel::{BlockingStrategy, ThreadPoolConfig};
use quasix_core::df::simd_ops::{check_simd_support, simd_matmul_transpose, simd_transform_3center};
use std::time::Duration;

/// GW100 test molecules with realistic dimensions
struct GW100Molecule {
    name: &'static str,
    nao: usize,
    naux: usize,
    nocc: usize,
    nvir: usize,
}

const GW100_SUBSET: &[GW100Molecule] = &[
    GW100Molecule { name: "H2",   nao: 2,  naux: 28,  nocc: 1,  nvir: 1  },
    GW100Molecule { name: "H2O",  nao: 7,  naux: 76,  nocc: 5,  nvir: 2  },
    GW100Molecule { name: "NH3",  nao: 8,  naux: 90,  nocc: 5,  nvir: 3  },
    GW100Molecule { name: "CH4",  nao: 9,  naux: 104, nocc: 5,  nvir: 4  },
    GW100Molecule { name: "CO",   nao: 10, naux: 112, nocc: 7,  nvir: 3  },
    GW100Molecule { name: "N2",   nao: 10, naux: 112, nocc: 7,  nvir: 3  },
    GW100Molecule { name: "C2H2", nao: 13, naux: 140, nocc: 7,  nvir: 6  },
    GW100Molecule { name: "H2CO", nao: 12, naux: 132, nocc: 8,  nvir: 4  },
];

fn generate_mock_j3c(n_ao: usize, n_aux: usize) -> Array3<f64> {
    use rand::{Rng, SeedableRng};
    let mut rng = rand::rngs::StdRng::seed_from_u64(42);

    let mut j3c = Array3::<f64>::zeros((n_ao, n_ao, n_aux));
    for i in 0..n_ao {
        for j in 0..=i {
            for p in 0..n_aux {
                let val = rng.random::<f64>() * 0.1;
                j3c[[i, j, p]] = val;
                j3c[[j, i, p]] = val; // Symmetry
            }
        }
    }
    j3c
}

/// Benchmark 1: SIMD vs Scalar Matrix Multiplication
fn bench_simd_matmul(c: &mut Criterion) {
    let mut group = c.benchmark_group("simd_matmul");
    group.warm_up_time(Duration::from_secs(1));
    group.measurement_time(Duration::from_secs(3));

    let simd_support = check_simd_support();
    println!("SIMD Support: {:?}", simd_support);

    for size in [32, 64, 128, 256].iter() {
        let m = *size;
        let k = *size;
        let n = *size;

        let a = Array2::from_shape_fn((m, k), |(i, j)| (i + j) as f64 * 0.01);
        let b = Array2::from_shape_fn((k, n), |(i, j)| (i * j) as f64 * 0.01);
        let mut c_simd = Array2::zeros((m, n));

        group.bench_with_input(
            BenchmarkId::new("simd", size),
            size,
            |bench, _| {
                bench.iter(|| {
                    simd_matmul_transpose(
                        black_box(a.view()),
                        black_box(b.view()),
                        black_box(c_simd.view_mut())
                    );
                })
            }
        );

        // Baseline: ndarray's native dot product
        group.bench_with_input(
            BenchmarkId::new("ndarray_dot", size),
            size,
            |bench, _| {
                bench.iter(|| {
                    let _ = black_box(a.dot(&b));
                })
            }
        );
    }

    group.finish();
}

/// Benchmark 2: SIMD 3-center transformation
fn bench_simd_3center(c: &mut Criterion) {
    let mut group = c.benchmark_group("simd_3center_transform");
    group.warm_up_time(Duration::from_secs(1));
    group.measurement_time(Duration::from_secs(5));

    for mol in GW100_SUBSET.iter().take(4) {  // Test on smaller molecules
        let j3c_ao = generate_mock_j3c(mol.nao, mol.naux);
        let j3c_flat = j3c_ao.into_shape((mol.nao * mol.nao, mol.naux)).unwrap();

        let c_occ = generate_mock_mo_coefficients(mol.nao, mol.nocc, 100);
        let c_vir = generate_mock_mo_coefficients(mol.nao, mol.nvir, 101);

        group.bench_with_input(
            BenchmarkId::new("simd", mol.name),
            mol.name,
            |bench, _| {
                bench.iter(|| {
                    simd_transform_3center(
                        black_box(j3c_flat.view()),
                        black_box(c_occ.view()),
                        black_box(c_vir.view())
                    )
                })
            }
        );
    }

    group.finish();
}

/// Benchmark 3: Thread coordination (BLAS vs Rayon)
fn bench_thread_coordination(c: &mut Criterion) {
    let mut group = c.benchmark_group("thread_coordination");
    group.warm_up_time(Duration::from_secs(2));
    group.measurement_time(Duration::from_secs(5));
    group.sample_size(15);

    // Test H2CO (medium-small system)
    let mol = &GW100_SUBSET[7];
    let j3c_ao = generate_mock_j3c(mol.nao, mol.naux);
    let c_occ = generate_mock_mo_coefficients(mol.nao, mol.nocc, 100);
    let c_vir = generate_mock_mo_coefficients(mol.nao, mol.nvir, 101);

    // Test different BLAS thread counts with Rayon
    for blas_threads in [1, 2, 4].iter() {
        let mut config = ThreadPoolConfig::default();
        config.blas_threads = *blas_threads;

        group.bench_with_input(
            BenchmarkId::new("blas_threads", blas_threads),
            blas_threads,
            |bench, _| {
                bench.iter(|| {
                    use quasix_core::df::parallel::transform_mo_3center_parallel_optimized;
                    transform_mo_3center_parallel_optimized(
                        black_box(&j3c_ao),
                        black_box(&c_occ),
                        black_box(&c_vir),
                        Some(config.clone()),
                        Some(BlockingStrategy::default())
                    )
                })
            }
        );
    }

    group.finish();
}

/// Benchmark 4: Lock contention analysis
fn bench_lock_contention(c: &mut Criterion) {
    let mut group = c.benchmark_group("lock_contention");
    group.warm_up_time(Duration::from_secs(1));
    group.measurement_time(Duration::from_secs(5));
    group.sample_size(15);

    // Test on C2H2 (larger system)
    let mol = &GW100_SUBSET[6];
    let j3c_ao = generate_mock_j3c(mol.nao, mol.naux);
    let c_occ = generate_mock_mo_coefficients(mol.nao, mol.nocc, 100);
    let c_vir = generate_mock_mo_coefficients(mol.nao, mol.nvir, 101);

    // Current implementation uses RwLock
    group.bench_function("current_rwlock", |bench| {
        bench.iter(|| {
            use quasix_core::df::parallel::transform_mo_3center_parallel_optimized;
            transform_mo_3center_parallel_optimized(
                black_box(&j3c_ao),
                black_box(&c_occ),
                black_box(&c_vir),
                Some(ThreadPoolConfig::default()),
                Some(BlockingStrategy::default())
            )
        })
    });

    // Sequential for comparison
    group.bench_function("sequential_baseline", |bench| {
        bench.iter(|| {
            use quasix_core::df::mo_transform::transform_mo_3center;
            transform_mo_3center(
                black_box(&j3c_ao),
                black_box(&c_occ),
                black_box(&c_vir)
            )
        })
    });

    group.finish();
}

/// Benchmark 5: Full GW100 subset performance
fn bench_gw100_subset(c: &mut Criterion) {
    let mut group = c.benchmark_group("gw100_subset");
    group.warm_up_time(Duration::from_secs(1));
    group.measurement_time(Duration::from_secs(5));
    group.sample_size(10);

    for mol in GW100_SUBSET.iter() {
        let j3c_ao = generate_mock_j3c(mol.nao, mol.naux);
        let c_occ = generate_mock_mo_coefficients(mol.nao, mol.nocc, 100);
        let c_vir = generate_mock_mo_coefficients(mol.nao, mol.nvir, 101);

        group.bench_with_input(
            BenchmarkId::new("parallel", mol.name),
            &mol.name,
            |bench, _| {
                bench.iter(|| {
                    use quasix_core::df::parallel::transform_mo_3center_parallel_optimized;
                    transform_mo_3center_parallel_optimized(
                        black_box(&j3c_ao),
                        black_box(&c_occ),
                        black_box(&c_vir),
                        Some(ThreadPoolConfig::default()),
                        Some(BlockingStrategy::default())
                    )
                })
            }
        );
    }

    group.finish();
}

/// Benchmark 6: Memory access patterns
fn bench_memory_patterns(c: &mut Criterion) {
    let mut group = c.benchmark_group("memory_patterns");
    group.warm_up_time(Duration::from_secs(1));
    group.measurement_time(Duration::from_secs(3));

    // Test cache-friendly vs cache-unfriendly access
    let mol = &GW100_SUBSET[5];  // N2
    let j3c_ao = generate_mock_j3c(mol.nao, mol.naux);
    let c_occ = generate_mock_mo_coefficients(mol.nao, mol.nocc, 100);
    let c_vir = generate_mock_mo_coefficients(mol.nao, mol.nvir, 101);

    // Default blocking strategy
    group.bench_function("default_blocking", |bench| {
        bench.iter(|| {
            use quasix_core::df::parallel::transform_mo_3center_parallel_optimized;
            transform_mo_3center_parallel_optimized(
                black_box(&j3c_ao),
                black_box(&c_occ),
                black_box(&c_vir),
                Some(ThreadPoolConfig::default()),
                Some(BlockingStrategy::default())
            )
        })
    });

    // Intel Xeon optimized
    group.bench_function("intel_xeon_blocking", |bench| {
        bench.iter(|| {
            use quasix_core::df::parallel::transform_mo_3center_parallel_optimized;
            transform_mo_3center_parallel_optimized(
                black_box(&j3c_ao),
                black_box(&c_occ),
                black_box(&c_vir),
                Some(ThreadPoolConfig::default()),
                Some(BlockingStrategy::intel_xeon())
            )
        })
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_simd_matmul,
    bench_simd_3center,
    bench_thread_coordination,
    bench_lock_contention,
    bench_gw100_subset,
    bench_memory_patterns
);
criterion_main!(benches);
