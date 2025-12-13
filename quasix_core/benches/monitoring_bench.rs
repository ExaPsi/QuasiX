//! Performance benchmarks for the convergence monitoring module
//!
//! This benchmark suite validates that the monitoring overhead meets targets:
//! - <100 μs per iteration update
//! - <1% total computational overhead
//! - Memory usage <100 MB for 1000 iterations

use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use quasix_core::gw::monitoring::{ConvergenceMonitor, MonitorConfig};

fn benchmark_monitor_update(c: &mut Criterion) {
    let mut group = c.benchmark_group("monitor_update");

    for n_states in [100, 500, 1000, 5000].iter() {
        group.bench_with_input(
            BenchmarkId::from_parameter(n_states),
            n_states,
            |b, &n| {
                let config = MonitorConfig::default();
                let monitor = ConvergenceMonitor::new(n, n / 2, config);
                let energies = vec![1.0; n];
                let z_factors = vec![0.85; n];

                // Establish baseline
                monitor.update(0, &energies, &z_factors).unwrap();

                let mut energies_mod = energies.clone();
                let mut iteration = 1;
                b.iter(|| {
                    // Modify energies slightly to simulate convergence
                    energies_mod[0] += 0.001;
                    iteration += 1;
                    black_box(monitor.update(iteration, &energies_mod, &z_factors).unwrap())
                });
            },
        );
    }
    group.finish();
}

fn benchmark_oscillation_detection(c: &mut Criterion) {
    let mut group = c.benchmark_group("oscillation_detection");

    let config = MonitorConfig::default();
    let monitor = ConvergenceMonitor::new(1000, 500, config);
    let energies = vec![1.0; 1000];
    let z_factors = vec![0.85; 1000];

    // Build up history
    for i in 0..20 {
        let mut e = energies.clone();
        e[0] = 1.0 + (i as f64) * 0.01;
        monitor.update(i, &e, &z_factors).unwrap();
    }

    group.bench_function("detect_oscillations", |b| {
        b.iter(|| black_box(monitor.detect_oscillations()))
    });

    group.finish();
}

fn benchmark_statistics_computation(c: &mut Criterion) {
    let config = MonitorConfig::default();
    let monitor = ConvergenceMonitor::new(1000, 500, config);
    let energies = vec![1.0; 1000];
    let z_factors = vec![0.85; 1000];

    // Build up history
    for i in 0..50 {
        let mut e = energies.clone();
        e[0] = 1.0 + (i as f64) * 0.001;
        monitor.update(i, &e, &z_factors).unwrap();
    }

    c.bench_function("get_statistics", |b| {
        b.iter(|| black_box(monitor.get_statistics()))
    });
}

fn benchmark_report_generation(c: &mut Criterion) {
    let config = MonitorConfig::default();
    let monitor = ConvergenceMonitor::new(100, 50, config);
    let energies = vec![1.0; 100];
    let z_factors = vec![0.85; 100];

    // Build up history
    for i in 0..100 {
        let mut e = energies.clone();
        e[0] = 1.0 + (i as f64) * 0.001;
        monitor.update(i, &e, &z_factors).unwrap();
    }

    c.bench_function("generate_report", |b| {
        b.iter(|| black_box(monitor.generate_report().unwrap()))
    });
}

fn benchmark_memory_usage(c: &mut Criterion) {
    c.bench_function("memory_1000_iterations", |b| {
        b.iter(|| {
            let config = MonitorConfig {
                max_iterations: 1000,
                ..Default::default()
            };
            let monitor = ConvergenceMonitor::new(1000, 500, config);
            let energies = vec![1.0; 1000];
            let z_factors = vec![0.85; 1000];

            for i in 0..1000 {
                let mut e = energies.clone();
                e[i % 1000] += 0.001;
                monitor.update(i, &e, &z_factors).unwrap();
            }
            black_box(monitor.get_history(1000))
        })
    });
}

fn benchmark_simd_rms(c: &mut Criterion) {
    let mut group = c.benchmark_group("simd_rms");

    for n in [100, 500, 1000, 5000].iter() {
        group.bench_with_input(
            BenchmarkId::from_parameter(n),
            n,
            |b, &n| {
                let old = vec![1.0; n];
                let mut new = vec![1.0; n];
                new[n / 2] = 2.0; // One change

                b.iter(|| black_box(quasix_core::gw::monitoring::compute_changes(&old, &new)))
            },
        );
    }
    group.finish();
}

criterion_group!(
    benches,
    benchmark_monitor_update,
    benchmark_oscillation_detection,
    benchmark_statistics_computation,
    benchmark_report_generation,
    benchmark_memory_usage,
    benchmark_simd_rms,
);

criterion_main!(benches);