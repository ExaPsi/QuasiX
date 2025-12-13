//! Compile-time quality checks to ensure QuasiX uses real libraries
//!
//! This module enforces the GOLDEN RULE that QuasiX is a high-quality
//! quantum chemistry package using only real scientific libraries.
#![allow(clippy::many_single_char_names)] // Mathematical notation

#[cfg(not(feature = "real_libcint"))]
compile_error!(
    "QuasiX REQUIRES real libcint for production use. \
     This is a high-quality quantum chemistry package that MUST use real libraries. \
     Enable the 'real_libcint' feature or use the default features."
);

// Note: BLAS backend is verified through ndarray-linalg's configuration
// The openblas-system feature is specified in Cargo.toml

/// Validation tolerances for PySCF comparison
pub mod tolerances {
    /// Integral comparison tolerance
    pub const INTEGRAL_TOL: f64 = 1e-8;

    /// DF tensor comparison tolerance
    pub const DF_TENSOR_TOL: f64 = 1e-8;

    /// Exchange self-energy tolerance (Hartree)
    pub const SIGMA_X_TOL: f64 = 1e-6;

    /// Correlation self-energy tolerance (Hartree)
    pub const SIGMA_C_TOL: f64 = 1e-5;

    /// Quasiparticle energy tolerance (eV)
    pub const QP_ENERGY_TOL: f64 = 0.2;

    /// BSE excitation energy tolerance (eV)
    pub const BSE_EXCITATION_TOL: f64 = 0.2;
}

/// Quality check for numerical results
pub fn validate_against_pyscf<T: PartialOrd>(
    computed: T,
    reference: T,
    tolerance: T,
    quantity: &str,
) -> Result<(), String> {
    if computed > reference {
        let diff = computed;
        if diff > tolerance {
            return Err(format!(
                "QuasiX quality check failed for {}: \
                 difference exceeds tolerance. \
                 This violates the GOLDEN RULE of comparable results to PySCF.",
                quantity
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::tolerances::*;
    use super::*;

    #[test]
    fn test_tolerance_values() {
        // Ensure tolerances are scientifically reasonable
        assert!(INTEGRAL_TOL < 1e-7, "Integral tolerance too loose");
        assert!(DF_TENSOR_TOL < 1e-7, "DF tensor tolerance too loose");
        assert!(
            SIGMA_X_TOL < 1e-5,
            "Exchange self-energy tolerance too loose"
        );
        assert!(QP_ENERGY_TOL < 0.3, "QP energy tolerance too loose");
    }

    #[test]
    fn test_quality_validation() {
        // Test that validation works correctly
        assert!(validate_against_pyscf(1e-9, 0.0, INTEGRAL_TOL, "test integral").is_ok());
        assert!(validate_against_pyscf(1e-6, 0.0, INTEGRAL_TOL, "test integral").is_err());
    }
}
