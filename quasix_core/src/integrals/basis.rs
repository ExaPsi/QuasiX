//! Basis set representation for integral calculations
#![allow(clippy::many_single_char_names)] // Mathematical notation

use crate::common::Result;

/// Basis function type
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BasisType {
    /// s-type orbital
    S,
    /// p-type orbital (px, py, pz)
    P,
    /// d-type orbital (6 components)
    D,
    /// f-type orbital (10 components)
    F,
}

impl BasisType {
    /// Get the number of components for this basis type
    pub fn n_components(&self) -> usize {
        match self {
            BasisType::S => 1,
            BasisType::P => 3,
            BasisType::D => 6,
            BasisType::F => 10,
        }
    }
}

/// Gaussian basis function
#[derive(Debug, Clone)]
pub struct BasisFunction {
    /// Angular momentum type
    pub basis_type: BasisType,
    /// Gaussian exponents
    pub exponents: Vec<f64>,
    /// Contraction coefficients
    pub coefficients: Vec<f64>,
    /// Center position (atom index)
    pub center: usize,
}

impl BasisFunction {
    /// Create a new basis function
    #[must_use]
    pub fn new(
        basis_type: BasisType,
        exponents: Vec<f64>,
        coefficients: Vec<f64>,
        center: usize,
    ) -> Self {
        assert_eq!(
            exponents.len(),
            coefficients.len(),
            "Exponents and coefficients must have same length"
        );

        Self {
            basis_type,
            exponents,
            coefficients,
            center,
        }
    }
}

/// Basis set representation
#[derive(Debug, Clone)]
pub struct BasisSet {
    /// List of basis functions
    pub functions: Vec<BasisFunction>,
    /// Total number of basis functions (including components)
    pub n_functions: usize,
    /// Total number of basis functions (alias for n_functions)
    size: usize,
    /// Basis set name
    name: String,
}

impl BasisSet {
    /// Create a new basis set
    #[must_use]
    pub fn new(functions: Vec<BasisFunction>, name: impl Into<String>) -> Self {
        let size = functions.iter().map(|f| f.basis_type.n_components()).sum();

        Self {
            functions,
            n_functions: size,
            size,
            name: name.into(),
        }
    }

    /// Get the total number of basis functions
    pub fn size(&self) -> usize {
        self.size
    }

    /// Get the basis set name
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Create a mock STO-3G basis set with specified size
    pub fn mock_sto3g(size: usize) -> Self {
        // Create mock basis functions
        // For simplicity, we'll create mostly s-type functions
        let mut functions = Vec::new();
        let mut current_size = 0;

        while current_size < size {
            let basis_type = if current_size == 0 || (size - current_size) < 3 {
                BasisType::S
            } else if (size - current_size) >= 3 && functions.len() % 3 == 0 {
                BasisType::P
            } else {
                BasisType::S
            };

            // Mock STO-3G exponents and coefficients
            let exponents = vec![3.425_250_91, 0.623_913_73, 0.168_855_40];
            let coefficients = vec![0.154_328_97, 0.535_328_14, 0.444_634_54];

            functions.push(BasisFunction::new(
                basis_type,
                exponents,
                coefficients,
                functions.len() % 3, // Mock center assignment
            ));

            current_size += basis_type.n_components();
        }

        // Trim to exact size if we overshot
        if current_size > size {
            functions.pop();
            current_size = functions.iter().map(|f| f.basis_type.n_components()).sum();
        }

        Self {
            functions,
            n_functions: current_size,
            size: current_size,
            name: "STO-3G".to_string(),
        }
    }

    /// Create a mock RI basis set with specified size
    pub fn mock_ri(size: usize) -> Self {
        // Create mock auxiliary basis functions
        // RI basis sets typically have more functions
        let mut functions = Vec::new();
        let mut current_size = 0;

        while current_size < size {
            let basis_type = match functions.len() % 4 {
                1 if (size - current_size) >= 3 => BasisType::P,
                2 if (size - current_size) >= 6 => BasisType::D,
                _ => BasisType::S,
            };

            // Mock RI exponents (typically more diffuse and more contracted)
            let exponents = vec![10.0, 4.0, 1.5, 0.5, 0.15];
            let coefficients = vec![0.1, 0.2, 0.3, 0.25, 0.15];

            functions.push(BasisFunction::new(
                basis_type,
                exponents,
                coefficients,
                functions.len() % 3,
            ));

            current_size += basis_type.n_components();
        }

        // Trim to exact size if we overshot
        while current_size > size && !functions.is_empty() {
            let last = functions.pop().unwrap();
            current_size -= last.basis_type.n_components();
        }

        Self {
            functions,
            n_functions: current_size,
            size: current_size,
            name: "RI".to_string(),
        }
    }

    /// Create mock basis for specific molecules
    pub fn for_molecule(molecule_name: &str, is_auxiliary: bool) -> Self {
        match (molecule_name, is_auxiliary) {
            ("H2O", false) => Self::mock_sto3g(8),
            ("H2O", true) => Self::mock_ri(18),
            ("NH3", false) => Self::mock_sto3g(11),
            ("NH3", true) => Self::mock_ri(23),
            ("CO", true) => Self::mock_ri(20),
            ("C6H6", false) => Self::mock_sto3g(36),
            ("C6H6", true) => Self::mock_ri(72),
            _ => Self::mock_sto3g(10), // Default (covers CO regular basis and others)
        }
    }

    /// Validate the basis set
    pub fn validate(&self) -> Result<()> {
        if self.functions.is_empty() {
            return Err(super::IntegralError::InvalidBasis(
                "Basis set has no functions".to_string(),
            )
            .into());
        }

        for (i, func) in self.functions.iter().enumerate() {
            if func.exponents.is_empty() {
                return Err(super::IntegralError::InvalidBasis(format!(
                    "Function {} has no exponents",
                    i
                ))
                .into());
            }

            if func.exponents.len() != func.coefficients.len() {
                return Err(super::IntegralError::InvalidBasis(format!(
                    "Function {} has mismatched exponents and coefficients",
                    i
                ))
                .into());
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basis_type_components() {
        assert_eq!(BasisType::S.n_components(), 1);
        assert_eq!(BasisType::P.n_components(), 3);
        assert_eq!(BasisType::D.n_components(), 6);
        assert_eq!(BasisType::F.n_components(), 10);
    }

    #[test]
    fn test_mock_sto3g() {
        let basis = BasisSet::mock_sto3g(8);
        assert_eq!(basis.size(), 8);
        assert_eq!(basis.name(), "STO-3G");
        assert!(basis.validate().is_ok());
    }

    #[test]
    fn test_mock_ri() {
        let basis = BasisSet::mock_ri(18);
        assert_eq!(basis.size(), 18);
        assert_eq!(basis.name(), "RI");
        assert!(basis.validate().is_ok());
    }

    #[test]
    fn test_molecule_specific_basis() {
        let h2o_basis = BasisSet::for_molecule("H2O", false);
        assert_eq!(h2o_basis.size(), 8);

        let h2o_aux = BasisSet::for_molecule("H2O", true);
        assert_eq!(h2o_aux.size(), 18);

        let nh3_basis = BasisSet::for_molecule("NH3", false);
        assert_eq!(nh3_basis.size(), 11);

        let nh3_aux = BasisSet::for_molecule("NH3", true);
        assert_eq!(nh3_aux.size(), 23);
    }
}
