//! Analysis module for post-processing and spectral analysis
//!
//! This module provides tools for analyzing quasiparticle calculations,
//! including spectral function generation, broadening schemes, and
//! data export functionality.

pub mod spectral;

// Re-export main types for convenience
pub use spectral::{
    BroadeningParams, SpectralAnalyzer, SpectralData, SpectralError, SpectralFunction,
    SpectralMetadata,
};

// Re-export broadening functions
pub use spectral::{gaussian_broadening, lorentzian_broadening, voigt_profile};
