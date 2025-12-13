use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use std::hint::black_box;
use ndarray::{Array1, Array2, Array3};
use num_complex::Complex64;
use quasix_core::selfenergy::{
    WScreenedStorage, WScreenedStorageSimd, InterpolationMethod,
    EnergyDependentSelfEnergy, EnergyDependentSelfEnergySimd,
    ContourDeformationConfig, DFTensorBlocked,
};
use rand::prelude::*;
use rand::thread_rng;

/// Generate test W matrices on imaginary frequency grid
fn generate_test_w_matrices(n_freq: usize, n_aux: usize) -> Vec<Array2<Complex64>> {
    let mut rng = thread_rng();

    (0..n_freq)
        .map(|_| {
            let mut mat = Array2::<Complex64>::zeros((n_aux, n_aux));
            for i in 0..n_aux {
                for j in i..n_aux {
                    let val = Complex64::new(
                        rng.gen_range(-1.0..1.0),
                        rng.gen_range(-0.1..0.1),
                    );
                    mat[[i, j]] = val;
                    mat[[j, i]] = val.conj(); // Hermitian
                }
            }
            mat
        })
        .collect()
}

/// Generate test molecular orbital data
fn generate_test_mo_data(n_mo: usize, n_occ: usize) -> (Array1<f64>, Array1<f64>) {
    let mut rng = thread_rng();

    let mo_energy: Array1<f64> = Array1::from_iter(
        (0..n_mo).map(|i| {
            if i < n_occ {
                -20.0 + i as f64 * 2.0 // Occupied
            } else {
                1.0 + (i - n_occ) as f64 * 3.0 // Virtual
            }
        })
    );

    let mo_occ: Array1<f64> = Array1::from_iter(
        (0..n_mo).map(|i| if i < n_occ { 2.0 } else { 0.0 })
    );

    (mo_energy, mo_occ)
}

/// Generate test DF tensor
fn generate_test_df_tensor(n_mo: usize, n_aux: usize) -> Array3<f64> {
    let mut rng = thread_rng();
    let mut tensor = Array3::<f64>::zeros((n_mo, n_mo, n_aux));

    for i in 0..n_mo {
        for j in 0..n_mo {
            for p in 0..n_aux {
                tensor[[i, j, p]] = rng.gen_range(-0.1..0.1);
            }
        }
    }

    tensor
}

fn benchmark_linear_interpolation(c: &mut Criterion) {
    let mut group = c.benchmark_group("w_linear_interpolation");

    for n_aux in [50, 100, 200].iter() {
        let n_freq = 40;
        let w_matrices = generate_test_w_matrices(n_freq, *n_aux);
        let xi_grid = Array1::linspace(0.1, 100.0, n_freq);

        // Original implementation
        let w_storage = WScreenedStorage::new(
            w_matrices.clone(),
            xi_grid.clone(),
            InterpolationMethod::Linear,
        ).unwrap();

        // SIMD implementation
        let w_storage_simd = WScreenedStorageSimd::from_matrices(
            w_matrices,
            xi_grid,
            None,
        ).unwrap();

        let test_omega = Complex64::new(0.0, 5.0);

        group.throughput(Throughput::Elements((*n_aux * n_aux) as u64));

        group.bench_function(BenchmarkId::new("original", n_aux), |b| {
            b.iter(|| {
                black_box(w_storage.evaluate(test_omega))
            })
        });

        group.bench_function(BenchmarkId::new("simd", n_aux), |b| {
            b.iter(|| {
                black_box(w_storage_simd.linear_interpolation_simd(test_omega))
            })
        });
    }

    group.finish();
}

fn benchmark_cubic_interpolation(c: &mut Criterion) {
    let mut group = c.benchmark_group("w_cubic_interpolation");

    for n_aux in [50, 100, 200].iter() {
        let n_freq = 40;
        let w_matrices = generate_test_w_matrices(n_freq, *n_aux);
        let xi_grid = Array1::linspace(0.1, 100.0, n_freq);

        // Original implementation
        let w_storage = WScreenedStorage::new(
            w_matrices.clone(),
            xi_grid.clone(),
            InterpolationMethod::Cubic,
        ).unwrap();

        // SIMD implementation
        let w_storage_simd = WScreenedStorageSimd::from_matrices(
            w_matrices,
            xi_grid,
            None,
        ).unwrap();

        let test_omega = Complex64::new(0.0, 5.0);

        group.throughput(Throughput::Elements((*n_aux * n_aux) as u64));

        group.bench_function(BenchmarkId::new("original", n_aux), |b| {
            b.iter(|| {
                black_box(w_storage.evaluate(test_omega))
            })
        });

        group.bench_function(BenchmarkId::new("simd", n_aux), |b| {
            b.iter(|| {
                black_box(w_storage_simd.cubic_interpolation_simd(test_omega))
            })
        });
    }

    group.finish();
}

fn benchmark_mo_transformation(c: &mut Criterion) {
    let mut group = c.benchmark_group("w_mo_transformation");

    for n_aux in [50, 100, 200].iter() {
        let n_mo = 100;
        let n_occ = 50;
        let n_freq = 20;

        let w_matrices = generate_test_w_matrices(n_freq, *n_aux);
        let xi_grid = Array1::linspace(0.1, 100.0, n_freq);
        let (mo_energy, mo_occ) = generate_test_mo_data(n_mo, n_occ);
        let df_tensor = generate_test_df_tensor(n_mo, *n_aux);

        let config = ContourDeformationConfig {
            n_imag_points: 40,
            xi_max: 100.0,
            eta: 1e-3,
            pole_threshold: 0.1,
            convergence_tol: 1e-6,
            use_simd: false,
            n_threads: None,
            compute_spectral: false,
        };

        // Original implementation
        let w_storage = WScreenedStorage::new(
            w_matrices.clone(),
            xi_grid.clone(),
            InterpolationMethod::Cubic,
        ).unwrap();

        let evaluator = EnergyDependentSelfEnergy::new(
            w_storage,
            mo_energy.clone(),
            mo_occ.clone(),
            df_tensor.clone(),
            config.clone(),
        ).unwrap();

        // SIMD implementation
        let w_storage_simd = WScreenedStorageSimd::from_matrices(
            w_matrices,
            xi_grid,
            None,
        ).unwrap();

        let evaluator_simd = EnergyDependentSelfEnergySimd::new(
            w_storage_simd,
            mo_energy,
            mo_occ,
            df_tensor,
            config.eta,
        ).unwrap();

        let test_orbital = 25;
        let test_energy = 0.0;

        group.throughput(Throughput::Elements((*n_aux * n_aux) as u64));

        group.bench_function(BenchmarkId::new("original", n_aux), |b| {
            b.iter(|| {
                black_box(evaluator.evaluate_orbital(test_orbital, test_energy))
            })
        });

        group.bench_function(BenchmarkId::new("simd", n_aux), |b| {
            b.iter(|| {
                black_box(evaluator_simd.evaluate_orbital_simd(test_orbital, test_energy))
            })
        });
    }

    group.finish();
}

fn benchmark_parallel_orbital_evaluation(c: &mut Criterion) {
    let mut group = c.benchmark_group("parallel_orbital_evaluation");

    let n_aux = 100;
    let n_mo = 100;
    let n_occ = 50;
    let n_freq = 20;

    let w_matrices = generate_test_w_matrices(n_freq, n_aux);
    let xi_grid = Array1::linspace(0.1, 100.0, n_freq);
    let (mo_energy, mo_occ) = generate_test_mo_data(n_mo, n_occ);
    let df_tensor = generate_test_df_tensor(n_mo, n_aux);

    // SIMD implementation
    let w_storage_simd = WScreenedStorageSimd::from_matrices(
        w_matrices,
        xi_grid,
        None,
    ).unwrap();

    let evaluator_simd = EnergyDependentSelfEnergySimd::new(
        w_storage_simd,
        mo_energy,
        mo_occ,
        df_tensor,
        1e-3,
    ).unwrap();

    for n_orbitals in [1, 4, 8, 16].iter() {
        let orbital_indices: Vec<usize> = (0..*n_orbitals).collect();
        let test_energy = 0.0;

        group.throughput(Throughput::Elements(*n_orbitals as u64));

        group.bench_function(BenchmarkId::new("parallel", n_orbitals), |b| {
            b.iter(|| {
                black_box(evaluator_simd.evaluate_orbitals_parallel(&orbital_indices, test_energy))
            })
        });

        group.bench_function(BenchmarkId::new("sequential", n_orbitals), |b| {
            b.iter(|| {
                let results: Vec<_> = orbital_indices
                    .iter()
                    .map(|&idx| evaluator_simd.evaluate_orbital_simd(idx, test_energy))
                    .collect();
                black_box(results)
            })
        });
    }

    group.finish();
}

fn benchmark_df_tensor_blocking(c: &mut Criterion) {
    let mut group = c.benchmark_group("df_tensor_blocking");

    for n_aux in [100, 200, 400].iter() {
        let n_mo = 100;
        let df_tensor = generate_test_df_tensor(n_mo, *n_aux);

        // Test different block sizes
        for block_size in [16, 32, 64].iter() {
            let blocked = DFTensorBlocked::from_tensor(
                df_tensor.clone(),
                Some(*block_size),
                Some(*block_size),
            );

            let mut rng = thread_rng();
            let test_indices: Vec<(usize, usize, usize)> = (0..1000)
                .map(|_| (
                    rng.gen_range(0..n_mo),
                    rng.gen_range(0..n_mo),
                    rng.gen_range(0..*n_aux),
                ))
                .collect();

            group.throughput(Throughput::Elements(1000));

            group.bench_function(
                BenchmarkId::new(format!("blocked_{}", block_size), n_aux),
                |b| {
                    b.iter(|| {
                        for &(i, j, p) in &test_indices {
                            black_box(blocked.get(i, j, p));
                        }
                    })
                },
            );

            group.bench_function(
                BenchmarkId::new("direct", n_aux),
                |b| {
                    b.iter(|| {
                        for &(i, j, p) in &test_indices {
                            black_box(df_tensor[[i, j, p]]);
                        }
                    })
                },
            );
        }
    }

    group.finish();
}

fn benchmark_cache_efficiency(c: &mut Criterion) {
    let mut group = c.benchmark_group("cache_efficiency");

    let n_aux = 200;
    let n_freq = 40;
    let w_matrices = generate_test_w_matrices(n_freq, n_aux);
    let xi_grid = Array1::linspace(0.1, 100.0, n_freq);

    // Test different block sizes for cache efficiency
    for block_size in [32, 64, 128].iter() {
        let w_storage_simd = WScreenedStorageSimd::from_matrices(
            w_matrices.clone(),
            xi_grid.clone(),
            Some(*block_size),
        ).unwrap();

        // Generate multiple test frequencies
        let test_frequencies: Vec<Complex64> = (0..100)
            .map(|i| Complex64::new(0.0, 0.5 + i as f64 * 0.5))
            .collect();

        group.throughput(Throughput::Elements(100));

        group.bench_function(BenchmarkId::new("block_size", block_size), |b| {
            b.iter(|| {
                for &omega in &test_frequencies {
                    black_box(w_storage_simd.cubic_interpolation_simd(omega).unwrap());
                }
            })
        });
    }

    group.finish();
}

criterion_group!(
    benches,
    benchmark_linear_interpolation,
    benchmark_cubic_interpolation,
    benchmark_mo_transformation,
    benchmark_parallel_orbital_evaluation,
    benchmark_df_tensor_blocking,
    benchmark_cache_efficiency
);
criterion_main!(benches);