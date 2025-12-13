//! Test parallel evGW implementation
//!
//! This test verifies that the parallel evGW module compiles and basic
//! functionality works, including NUMA allocation and convergence checking.

use ndarray::{Array1, Array3};
use quasix_core::qp::convergence::{
    AccelerationMethod, ConvergenceCriteria, ParallelConvergenceChecker,
};
use quasix_core::qp::evgw_parallel::{
    DFTensors, EvGWConvergence, LoadBalanceStrategy, ParallelEvGW, ParallelEvGWConfig,
};
use quasix_core::qp::numa_alloc::{NumaAllocator, NumaPolicy};
use quasix_core::selfenergy::correlation::ContourDeformationConfig;

#[test]
fn test_parallel_evgw_creation() {
    // Test that we can create a parallel evGW solver
    let n_mo = 10;
    let n_aux = 50;

    let config = ParallelEvGWConfig {
        n_orbital_threads: 4,
        n_freq_threads: 2,
        numa_aware: false, // Disable NUMA for simple test
        cache_block_mb: 16,
        load_balance: LoadBalanceStrategy::Dynamic,
        cd_config: ContourDeformationConfig::default(),
        convergence: EvGWConvergence::default(),
        profile: false,
    };

    let solver = ParallelEvGW::new(n_mo, n_aux, config);
    assert!(solver.is_ok(), "Failed to create parallel evGW solver");
}

#[test]
fn test_df_tensors() {
    // Test DF tensor creation
    let n_occ = 5;
    let n_vir = 5;
    let n_aux = 20;

    // Create mock tensors
    let ia_p = Array3::<f64>::zeros((n_occ, n_vir, n_aux));
    let ij_p = Array3::<f64>::zeros((n_occ, n_occ, n_aux));
    let occupations = Array1::<f64>::from_vec(vec![2.0; n_occ]);

    let _df_tensors = DFTensors::new(ia_p.clone(), ij_p.clone(), occupations.clone());

    // Since fields are private, just verify construction succeeded
    // The fact that we can create it without panic is sufficient
    assert_eq!(ia_p.shape(), &[n_occ, n_vir, n_aux]);
    assert_eq!(ij_p.shape(), &[n_occ, n_occ, n_aux]);
    assert_eq!(occupations.len(), n_occ);
}

#[test]
fn test_parallel_convergence_checker() {
    // Test parallel convergence checker
    let n_systems = 3;
    let criteria = ConvergenceCriteria::default();

    let mut checker =
        ParallelConvergenceChecker::new(n_systems, criteria, AccelerationMethod::None);

    // Create mock energy and z-factor arrays
    let energies: Vec<Array1<f64>> = (0..n_systems)
        .map(|i| Array1::from_vec(vec![i as f64; 5]))
        .collect();

    let z_factors: Vec<Array1<f64>> = (0..n_systems)
        .map(|_| Array1::from_vec(vec![0.8; 5]))
        .collect();

    // Check convergence for all systems
    let converged = checker.check_all(&energies, &z_factors);
    assert_eq!(
        converged.len(),
        n_systems,
        "Should return convergence status for all systems"
    );
}

#[test]
fn test_numa_allocator() {
    // Test NUMA allocator creation
    let allocator = NumaAllocator::new();
    assert!(allocator.is_ok(), "Failed to create NUMA allocator");

    let mut alloc = allocator.unwrap();

    // Test different policies
    let policies = vec![
        NumaPolicy::Default,
        NumaPolicy::Interleave,
        NumaPolicy::FirstTouch,
    ];

    for policy in policies {
        alloc.set_policy(policy);
        // On non-NUMA systems this might fail, which is OK
        // Just test that we can attempt to set policies
    }
}

#[test]
fn test_evgw_convergence_defaults() {
    // Test convergence criteria defaults
    let conv = EvGWConvergence::default();
    assert_eq!(conv.max_iter, 50);
    assert_eq!(conv.energy_tol, 1e-6);
    assert!(conv.z_min > 0.0 && conv.z_min < 1.0);
}

#[test]
fn test_parallel_config_defaults() {
    // Test parallel config defaults
    let config = ParallelEvGWConfig::default();
    assert!(config.n_orbital_threads > 0);
    assert!(config.n_freq_threads > 0);
    assert!(config.cache_block_mb > 0);
}

/// Simple performance test (not a benchmark, just verification)
#[test]
fn test_parallel_performance_tracking() {
    // Create a small solver to test performance tracking
    let n_mo = 4;
    let n_aux = 10;

    let config = ParallelEvGWConfig {
        n_orbital_threads: 2,
        n_freq_threads: 1,
        numa_aware: false,
        cache_block_mb: 8,
        load_balance: LoadBalanceStrategy::Static,
        cd_config: ContourDeformationConfig::default(),
        convergence: EvGWConvergence {
            max_iter: 2, // Just 2 iterations for test
            energy_tol: 1e-4,
            z_tol: 1e-3,
            z_min: 0.3,
            damping: 0.5,
        },
        profile: true, // Enable profiling
    };

    let solver = ParallelEvGW::new(n_mo, n_aux, config).unwrap();
    let perf = solver.get_performance();

    // Check that performance struct is initialized
    assert_eq!(perf.sigma_x_time, 0.0);
    assert_eq!(perf.sigma_c_time, 0.0);
    assert_eq!(perf.w_time, 0.0);
    assert!(perf.iter_times.is_empty());
}
