//! Integration tests for spectral analysis module

use ndarray::Array1;
use quasix_core::analysis::{
    BroadeningParams, SpectralAnalyzer, SpectralData, SpectralFunction,
    gaussian_broadening, lorentzian_broadening, voigt_profile,
};
use approx::assert_relative_eq;

#[test]
fn test_spectral_analyzer_creation() {
    // Test data: simple 5-state system
    let energies = Array1::from_vec(vec![-2.0, -1.0, 0.0, 1.0, 2.0]);
    let z_factors = Array1::from_vec(vec![0.8, 0.75, 0.7, 0.75, 0.8]);
    let fermi_energy = 0.0;

    let analyzer = SpectralAnalyzer::new(energies, z_factors, fermi_energy);
    assert!(analyzer.is_ok(), "Failed to create spectral analyzer");
}

#[test]
fn test_invalid_z_factors() {
    let energies = Array1::from_vec(vec![-1.0, 0.0, 1.0]);
    let invalid_z = Array1::from_vec(vec![0.5, 1.5, 0.3]); // 1.5 is invalid

    let analyzer = SpectralAnalyzer::new(energies, invalid_z, 0.0);
    assert!(analyzer.is_err(), "Should reject invalid Z-factors > 1");
}

#[test]
fn test_dimension_mismatch() {
    let energies = Array1::from_vec(vec![-1.0, 0.0, 1.0]);
    let z_factors = Array1::from_vec(vec![0.5, 0.8]); // Wrong size

    let analyzer = SpectralAnalyzer::new(energies, z_factors, 0.0);
    assert!(analyzer.is_err(), "Should reject dimension mismatch");
}

#[test]
fn test_broadening_parameters() {
    let mut params = BroadeningParams::default();
    assert!(params.validate().is_ok());

    // Test invalid parameters
    params.gaussian_width = -1.0;
    assert!(params.validate().is_err(), "Should reject negative width");

    params.gaussian_width = 0.1;
    params.voigt_mixing = 2.0; // Out of [0,1] range
    assert!(params.validate().is_err(), "Should reject invalid mixing");
}

#[test]
fn test_spectral_function_generation() {
    let energies = Array1::from_vec(vec![-2.0, -1.0, 0.0, 1.0, 2.0]);
    let z_factors = Array1::from_vec(vec![0.8, 0.75, 0.7, 0.75, 0.8]);
    let fermi_energy = 0.0;

    let analyzer = SpectralAnalyzer::new(energies, z_factors, fermi_energy).unwrap();
    let omega_grid = Array1::linspace(-3.0, 3.0, 100);

    let spectral_fn = analyzer.compute_spectral_function(&omega_grid);
    assert!(spectral_fn.is_ok(), "Failed to compute spectral function");

    let spectral_fn = spectral_fn.unwrap();
    assert_eq!(spectral_fn.omega_grid.len(), 100);
    assert_eq!(spectral_fn.spectral_weights.len(), 100);
    assert!(spectral_fn.metadata.total_weight > 0.0);
}

#[test]
fn test_pes_ipes_separation() {
    // System with clear occupied/unoccupied separation
    let energies = Array1::from_vec(vec![-2.0, -1.0, 0.5, 1.0, 2.0]);
    let z_factors = Array1::from_vec(vec![0.8, 0.7, 0.6, 0.7, 0.8]);
    let fermi_energy = 0.0;

    let analyzer = SpectralAnalyzer::new(energies, z_factors, fermi_energy).unwrap();

    // Test PES (photoemission - occupied states)
    let pes = analyzer.generate_pes_spectrum(0.0, 3.0, 50);
    assert!(pes.is_ok(), "Failed to generate PES spectrum");
    let pes = pes.unwrap();
    assert_eq!(pes.metadata.n_states, 2, "PES should have 2 occupied states");

    // Test IPES (inverse photoemission - unoccupied states)
    let ipes = analyzer.generate_ipes_spectrum(0.0, 3.0, 50);
    assert!(ipes.is_ok(), "Failed to generate IPES spectrum");
    let ipes = ipes.unwrap();
    assert_eq!(ipes.metadata.n_states, 3, "IPES should have 3 unoccupied states");
}

#[test]
fn test_broadening_functions() {
    let omega = 1.0;
    let energy = 0.0;
    let sigma = 0.1;
    let gamma = 0.05;

    // Test Gaussian broadening
    let gauss = gaussian_broadening(omega, energy, sigma);
    assert!(gauss > 0.0, "Gaussian broadening should be positive");

    // Test Lorentzian broadening
    let lorentz = lorentzian_broadening(omega, energy, gamma);
    assert!(lorentz > 0.0, "Lorentzian broadening should be positive");

    // Test Voigt profile limits
    let pure_gauss = voigt_profile(omega, energy, sigma, gamma, 0.0);
    let expected_gauss = gaussian_broadening(omega, energy, sigma);
    assert_relative_eq!(pure_gauss, expected_gauss, epsilon = 1e-10);

    let pure_lorentz = voigt_profile(omega, energy, sigma, gamma, 1.0);
    let expected_lorentz = lorentzian_broadening(omega, energy, gamma);
    assert_relative_eq!(pure_lorentz, expected_lorentz, epsilon = 1e-10);

    // Test mixture
    let mixed = voigt_profile(omega, energy, sigma, gamma, 0.5);
    assert!(mixed > 0.0, "Mixed Voigt profile should be positive");
    assert!(mixed < pure_gauss + pure_lorentz, "Mixed should be less than sum");
}

#[test]
fn test_broadening_normalization() {
    // Test that broadening functions are normalized
    let sigma = 0.1;
    let gamma = 0.05;
    let omega_grid = Array1::linspace(-10.0, 10.0, 2000);
    let dw = omega_grid[1] - omega_grid[0];

    // Gaussian normalization
    let mut gauss_sum = 0.0;
    for &omega in omega_grid.iter() {
        gauss_sum += gaussian_broadening(omega, 0.0, sigma) * dw;
    }
    assert_relative_eq!(gauss_sum, 1.0, epsilon = 0.01);

    // Lorentzian normalization
    let mut lorentz_sum = 0.0;
    for &omega in omega_grid.iter() {
        lorentz_sum += lorentzian_broadening(omega, 0.0, gamma) * dw;
    }
    assert_relative_eq!(lorentz_sum, 1.0, epsilon = 0.01);
}

#[test]
fn test_json_export() {
    use tempfile::NamedTempFile;
    use std::fs;

    let energies = Array1::from_vec(vec![-1.0, 0.0, 1.0]);
    let z_factors = Array1::from_vec(vec![0.5, 0.8, 0.3]);
    let analyzer = SpectralAnalyzer::new(energies, z_factors, 0.0).unwrap();

    let data = SpectralData {
        pes_spectrum: None,
        ipes_spectrum: None,
        combined_spectrum: None,
        timestamp: chrono::Utc::now().to_rfc3339(),
        method: "evGW".to_string(),
    };

    let temp_file = NamedTempFile::new().unwrap();
    let result = analyzer.export_json(&data, temp_file.path());
    assert!(result.is_ok(), "Failed to export JSON");

    // Verify file contents
    let contents = fs::read_to_string(temp_file.path()).unwrap();
    assert!(contents.contains("\"method\":\"evGW\""));
    assert!(contents.contains("timestamp"));
}

#[cfg(feature = "hdf5_support")]
#[test]
fn test_hdf5_export() {
    use tempfile::NamedTempFile;

    let energies = Array1::from_vec(vec![-1.0, 0.0, 1.0]);
    let z_factors = Array1::from_vec(vec![0.5, 0.8, 0.3]);
    let mut analyzer = SpectralAnalyzer::new(energies, z_factors, 0.0).unwrap();

    // Set custom broadening
    let broadening = BroadeningParams {
        gaussian_width: 0.15,
        lorentzian_width: 0.08,
        temperature: 300.0,
        voigt_mixing: 0.3,
    };
    analyzer.set_broadening(broadening).unwrap();

    // Generate spectra
    let pes = analyzer.generate_pes_spectrum(0.0, 2.0, 50).ok();
    let ipes = analyzer.generate_ipes_spectrum(0.0, 2.0, 50).ok();

    let data = SpectralData {
        pes_spectrum: pes,
        ipes_spectrum: ipes,
        combined_spectrum: None,
        timestamp: chrono::Utc::now().to_rfc3339(),
        method: "evGW".to_string(),
    };

    let temp_file = NamedTempFile::new().unwrap();
    let result = analyzer.export_hdf5(&data, temp_file.path());
    assert!(result.is_ok(), "Failed to export HDF5: {:?}", result);

    // Verify file exists and has non-zero size
    let metadata = std::fs::metadata(temp_file.path()).unwrap();
    assert!(metadata.len() > 0, "HDF5 file should not be empty");
}

#[test]
fn test_empty_omega_grid() {
    let energies = Array1::from_vec(vec![-1.0, 0.0, 1.0]);
    let z_factors = Array1::from_vec(vec![0.5, 0.8, 0.3]);
    let analyzer = SpectralAnalyzer::new(energies, z_factors, 0.0).unwrap();

    let empty_grid = Array1::zeros(0);
    let result = analyzer.compute_spectral_function(&empty_grid);
    assert!(result.is_err(), "Should reject empty omega grid");
}

#[test]
fn test_different_broadening_schemes() {
    let energies = Array1::from_vec(vec![-1.0, 0.0, 1.0]);
    let z_factors = Array1::from_vec(vec![0.5, 0.8, 0.3]);
    let mut analyzer = SpectralAnalyzer::new(energies, z_factors, 0.0).unwrap();
    let omega_grid = Array1::linspace(-2.0, 2.0, 100);

    // Test pure Gaussian
    let mut params = BroadeningParams::default();
    params.voigt_mixing = 0.0;
    analyzer.set_broadening(params).unwrap();
    let gauss_spectrum = analyzer.compute_spectral_function(&omega_grid).unwrap();

    // Test pure Lorentzian
    params.voigt_mixing = 1.0;
    analyzer.set_broadening(params).unwrap();
    let lorentz_spectrum = analyzer.compute_spectral_function(&omega_grid).unwrap();

    // Test mixed Voigt
    params.voigt_mixing = 0.5;
    analyzer.set_broadening(params).unwrap();
    let voigt_spectrum = analyzer.compute_spectral_function(&omega_grid).unwrap();

    // All should have same grid size
    assert_eq!(gauss_spectrum.omega_grid.len(), 100);
    assert_eq!(lorentz_spectrum.omega_grid.len(), 100);
    assert_eq!(voigt_spectrum.omega_grid.len(), 100);

    // All should have positive weights
    assert!(gauss_spectrum.spectral_weights.iter().all(|&w| w >= 0.0));
    assert!(lorentz_spectrum.spectral_weights.iter().all(|&w| w >= 0.0));
    assert!(voigt_spectrum.spectral_weights.iter().all(|&w| w >= 0.0));
}