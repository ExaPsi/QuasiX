use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use ndarray::Array2;
use num_complex::Complex64;
use quasix_core::dielectric::screening_enhanced::EnhancedDielectricSolver;
use std::hint::black_box;

fn benchmark_screening(c: &mut Criterion) {
    let mut group = c.benchmark_group("screening_algorithms");

    for naux in [10, 50, 100].iter() {
        group.bench_with_input(
            BenchmarkId::new("enhanced_adaptive", naux),
            naux,
            |b, &naux| {
                let mut solver = EnhancedDielectricSolver::new(naux);
                let p0 = Array2::from_elem((naux, naux), Complex64::new(0.1, 0.0));
                let v_sqrt = Array2::eye(naux);

                b.iter(|| {
                    solver.compute_screened_interaction_adaptive(black_box(&p0), black_box(&v_sqrt))
                });
            },
        );
    }

    group.finish();
}

criterion_group!(benches, benchmark_screening);
criterion_main!(benches);
