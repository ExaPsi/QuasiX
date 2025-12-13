use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use ndarray::{Array2, Array3};
use quasix_core::df::mo_transform::{transform_mo_3center, TransformConfig, transform_mo_3center_with_config};

/// Generate mock orthonormal MO coefficients
fn generate_mock_mo_coefficients(n_ao: usize, n_mo: usize, seed: u64) -> Array2<f64> {
    use rand::{Rng, SeedableRng};
    use rand_distr::StandardNormal;

    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    let mut coeffs = Array2::<f64>::zeros((n_ao, n_mo));

    for i in 0..n_ao {
        for j in 0..n_mo {
            coeffs[[i, j]] = rng.sample(StandardNormal);
        }
    }

    // Gram-Schmidt orthogonalization
    for j in 0..n_mo {
        let mut v = coeffs.column(j).to_owned();
        for k in 0..j {
            let u_k = coeffs.column(k);
            let proj = v.dot(&u_k);
            v = v - proj * &u_k;
        }
        let norm = v.dot(&v).sqrt();
        if norm > 1e-10 {
            v /= norm;
        }
        coeffs.column_mut(j).assign(&v);
    }

    coeffs
}

fn bench_mo_transform_various_sizes(c: &mut Criterion) {
    let mut group = c.benchmark_group("mo_transform_blas");

    // Test various molecule sizes from GW100
    let test_cases = vec![
        ("H2", 2, 28, 1, 1),      // Smallest
        ("H2O", 7, 76, 5, 2),     // Small
        ("NH3", 8, 90, 5, 3),     // Small
        ("CH4", 9, 104, 5, 4),    // Medium-small
        ("CO", 10, 112, 7, 3),    // Medium-small
    ];

    for (name, n_ao, n_aux, nocc, nvir) in test_cases {
        let j3c = Array3::<f64>::zeros((n_ao, n_ao, n_aux));
        let c_occ = generate_mock_mo_coefficients(n_ao, nocc, 42);
        let c_vir = generate_mock_mo_coefficients(n_ao, nvir, 43);

        group.bench_with_input(
            BenchmarkId::from_parameter(name),
            &(j3c, c_occ, c_vir),
            |b, (j3c, c_occ, c_vir)| {
                b.iter(|| {
                    transform_mo_3center(
                        black_box(j3c),
                        black_box(c_occ),
                        black_box(c_vir)
                    )
                })
            }
        );
    }

    group.finish();
}

fn bench_sequential_vs_blocked(c: &mut Criterion) {
    let mut group = c.benchmark_group("sequential_vs_blocked");

    // Medium-sized system where both algorithms apply
    let n_ao = 50;
    let n_aux = 100;
    let nocc = 10;
    let nvir = 40;

    let j3c = Array3::<f64>::zeros((n_ao, n_ao, n_aux));
    let c_occ = generate_mock_mo_coefficients(n_ao, nocc, 100);
    let c_vir = generate_mock_mo_coefficients(n_ao, nvir, 101);

    group.bench_function("sequential", |b| {
        let config = TransformConfig {
            use_blocking: false,
            use_advanced_parallel: false,
            ..Default::default()
        };
        b.iter(|| {
            transform_mo_3center_with_config(
                black_box(&j3c),
                black_box(&c_occ),
                black_box(&c_vir),
                &config
            )
        })
    });

    group.bench_function("blocked", |b| {
        let config = TransformConfig {
            use_blocking: true,
            use_advanced_parallel: false,
            max_memory_gb: 1.0,
            ..Default::default()
        };
        b.iter(|| {
            transform_mo_3center_with_config(
                black_box(&j3c),
                black_box(&c_occ),
                black_box(&c_vir),
                &config
            )
        })
    });

    group.finish();
}

criterion_group!(benches, bench_mo_transform_various_sizes, bench_sequential_vs_blocked);
criterion_main!(benches);
