//! Benchmark for unified screening implementation optimization
//!
//! This benchmark profiles the performance of the unified screening algorithm
//! across different frequency grid sizes to identify optimization opportunities.

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use ndarray::{Array1, Array2, Array3};
use num_complex::Complex64;
use quasix_core::dielectric::screening_enhanced::{
    EnhancedDielectricSolver, EnhancedSolverConfig,
};
use quasix_core::dielectric::screening_frequency_dependent::{
    FreqDependentScreening, FreqDependentScreeningConfig,
};

/// Generate test data for benchmarking
fn generate_test_data(naux: usize, nocc: usize, nvirt: usize) -> (Array1<f64>, Array3<f64>, Array2<f64>) {
    // MO energies: occupied negative, virtual positive
    let mut mo_energy = Array1::<f64>::zeros(nocc + nvirt);
    for i in 0..nocc {
        mo_energy[i] = -1.0 - (i as f64) * 0.1;
    }
    for i in 0..nvirt {
        mo_energy[nocc + i] = 0.5 + (i as f64) * 0.1;
    }

    // DF tensor (ia|P) with realistic sparsity pattern
    let mut df_ia = Array3::<f64>::zeros((nocc, nvirt, naux));
    for i in 0..nocc {
        for a in 0..nvirt {
            for p in 0..naux {
                // Create sparse pattern: most elements small, few large
                let value = if (i + a + p) % 10 == 0 {
                    0.5 + (i + a + p) as f64 * 0.01
                } else {
                    1e-4 * ((i + a + p) as f64).sin()
                };
                df_ia[[i, a, p]] = value;
            }
        }
    }

    // Coulomb matrix (positive definite)
    let mut v_aux = Array2::<f64>::zeros((naux, naux));
    for i in 0..naux {
        v_aux[[i, i]] = 2.0 + (i as f64) * 0.1;
        for j in 0..i {
            let value = 0.1 * ((i + j) as f64).sin().abs();
            v_aux[[i, j]] = value;
            v_aux[[j, i]] = value;
        }
    }

    (mo_energy, df_ia, v_aux)
}

/// Benchmark unified screening with different frequency counts
fn bench_unified_screening(c: &mut Criterion) {
    let mut group = c.benchmark_group("unified_screening");

    // Test different frequency grid sizes
    let freq_counts = vec![8, 16, 30, 60, 128];
    let naux = 100;
    let nocc = 10;
    let nvirt = 20;

    // Generate test data once
    let (mo_energy, df_ia, v_aux) = generate_test_data(naux, nocc, nvirt);

    for nfreq in freq_counts {
        group.bench_with_input(
            BenchmarkId::new("freq_dependent", nfreq),
            &nfreq,
            |b, &nfreq| {
                let config = FreqDependentScreeningConfig {
                    n_omega: nfreq,
                    omega_max: 20.0,
                    use_gl_grid: true,
                    parallel: false, // Single-threaded for consistent benchmarking
                    ..Default::default()
                };

                b.iter(|| {
                    let mut screening = FreqDependentScreening::new(
                        naux, nocc, nvirt, config.clone()
                    ).unwrap();

                    let _w_matrices = screening.compute_w_imag_axis(
                        &mo_energy,
                        &df_ia,
                        &v_aux,
                    ).unwrap();

                    black_box(_w_matrices);
                });
            },
        );
    }

    group.finish();
}

/// Benchmark P0 construction at different frequencies
fn bench_p0_construction(c: &mut Criterion) {
    let mut group = c.benchmark_group("p0_construction");

    let naux = 100;
    let nocc = 10;
    let nvirt = 20;
    let omega_values = vec![0.1, 1.0, 5.0, 10.0, 20.0];

    let (_mo_energy, df_ia, _v_aux) = generate_test_data(naux, nocc, nvirt);
    let e_occ = Array1::from_vec(vec![-1.0; nocc]);
    let e_virt = Array1::from_vec(vec![0.5; nvirt]);

    for omega in omega_values {
        group.bench_with_input(
            BenchmarkId::new("omega", format!("{:.1}", omega)),
            &omega,
            |b, &omega| {
                b.iter(|| {
                    let mut p0 = Array2::<f64>::zeros((naux, naux));

                    // Inner loop that needs optimization
                    for i in 0..nocc {
                        let e_i = e_occ[i];
                        for a in 0..nvirt {
                            let e_a = e_virt[a];
                            let de = e_i - e_a;
                            let denom = -2.0 * de / (de * de + omega * omega + 1e-14);

                            // This is the hot loop - needs SIMD/blocking
                            for p in 0..naux {
                                let df_iap = df_ia[[i, a, p]];
                                for q in 0..naux {
                                    p0[[p, q]] += denom * df_iap * df_ia[[i, a, q]];
                                }
                            }
                        }
                    }

                    black_box(p0);
                });
            },
        );
    }

    group.finish();
}

/// Benchmark dielectric matrix inversion with different sizes
fn bench_dielectric_inversion(c: &mut Criterion) {
    let mut group = c.benchmark_group("dielectric_inversion");

    let aux_sizes = vec![50, 100, 200, 400];

    for naux in aux_sizes {
        group.bench_with_input(
            BenchmarkId::new("enhanced_solver", naux),
            &naux,
            |b, &naux| {
                let config = EnhancedSolverConfig::default();
                let mut solver = EnhancedDielectricSolver::with_config(naux, config);

                // Create test matrices
                let mut p0 = Array2::<Complex64>::zeros((naux, naux));
                for i in 0..naux {
                    for j in 0..naux {
                        let val = if i == j {
                            0.5
                        } else {
                            0.01 * ((i + j) as f64).sin()
                        };
                        p0[[i, j]] = Complex64::new(val, 0.0);
                    }
                }

                let mut v_sqrt = Array2::<f64>::zeros((naux, naux));
                for i in 0..naux {
                    v_sqrt[[i, i]] = 1.414; // sqrt(2)
                }

                b.iter(|| {
                    let (_w, _diagnostics) = solver
                        .compute_screened_interaction_adaptive(&p0, &v_sqrt)
                        .unwrap();
                    black_box(_w);
                });
            },
        );
    }

    group.finish();
}

/// Benchmark memory access patterns
fn bench_memory_patterns(c: &mut Criterion) {
    let mut group = c.benchmark_group("memory_patterns");

    let naux = 200;
    let nocc = 20;
    let nvirt = 40;

    let (_mo_energy, df_ia, _v_aux) = generate_test_data(naux, nocc, nvirt);

    // Benchmark different loop orderings for P0 construction
    group.bench_function("loop_order_pqia", |b| {
        b.iter(|| {
            let mut p0 = Array2::<f64>::zeros((naux, naux));

            // P-Q outer, I-A inner (potentially better cache usage)
            for p in 0..naux {
                for q in 0..naux {
                    let mut sum = 0.0;
                    for i in 0..nocc {
                        for a in 0..nvirt {
                            let de = -0.5; // Dummy value
                            sum += de * df_ia[[i, a, p]] * df_ia[[i, a, q]];
                        }
                    }
                    p0[[p, q]] = sum;
                }
            }

            black_box(p0);
        });
    });

    group.bench_function("loop_order_iapq", |b| {
        b.iter(|| {
            let mut p0 = Array2::<f64>::zeros((naux, naux));

            // I-A outer, P-Q inner (original ordering)
            for i in 0..nocc {
                for a in 0..nvirt {
                    let de = -0.5; // Dummy value
                    for p in 0..naux {
                        let df_iap = df_ia[[i, a, p]];
                        for q in 0..naux {
                            p0[[p, q]] += de * df_iap * df_ia[[i, a, q]];
                        }
                    }
                }
            }

            black_box(p0);
        });
    });

    group.bench_function("loop_order_blocked", |b| {
        b.iter(|| {
            let mut p0 = Array2::<f64>::zeros((naux, naux));
            let block_size = 32;

            // Blocked version for better cache usage
            for p_block in (0..naux).step_by(block_size) {
                let p_end = (p_block + block_size).min(naux);
                for q_block in (0..naux).step_by(block_size) {
                    let q_end = (q_block + block_size).min(naux);

                    for i in 0..nocc {
                        for a in 0..nvirt {
                            let de = -0.5; // Dummy value
                            for p in p_block..p_end {
                                let df_iap = df_ia[[i, a, p]];
                                for q in q_block..q_end {
                                    p0[[p, q]] += de * df_iap * df_ia[[i, a, q]];
                                }
                            }
                        }
                    }
                }
            }

            black_box(p0);
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_unified_screening,
    bench_p0_construction,
    bench_dielectric_inversion,
    bench_memory_patterns
);
criterion_main!(benches);