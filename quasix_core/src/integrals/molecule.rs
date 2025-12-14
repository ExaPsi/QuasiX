//! Molecule representation for integral calculations
#![allow(clippy::many_single_char_names)] // Mathematical notation

use crate::common::Result;

/// Atomic element representation
#[derive(Debug, Clone)]
pub struct Atom {
    /// Atomic symbol (e.g., "H", "O", "N")
    pub symbol: String,
    /// Atomic number
    pub atomic_number: u8,
    /// Cartesian coordinates in Bohr
    pub coords: [f64; 3],
    /// Position (alias for coords)
    pub position: [f64; 3],
}

impl Atom {
    /// Create a new atom
    #[must_use]
    pub fn new(symbol: impl Into<String>, atomic_number: u8, coords: [f64; 3]) -> Self {
        Self {
            symbol: symbol.into(),
            atomic_number,
            coords,
            position: coords,
        }
    }
}

/// Molecule representation
#[derive(Debug, Clone)]
pub struct Molecule {
    /// List of atoms
    pub atoms: Vec<Atom>,
    /// Molecular charge
    pub charge: i32,
    /// Spin multiplicity (2S+1)
    pub multiplicity: u32,
    /// Optional molecule name
    name: Option<String>,
}

impl Molecule {
    /// Create a new molecule
    #[must_use]
    pub fn new(atoms: Vec<Atom>, charge: i32, multiplicity: u32) -> Self {
        Self {
            atoms,
            charge,
            multiplicity,
            name: None,
        }
    }

    /// Create a new molecule with a name
    #[must_use]
    pub fn with_name(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    /// Get the molecule name
    pub fn name(&self) -> &str {
        self.name.as_deref().unwrap_or("Unknown")
    }

    /// Get the number of atoms
    pub fn natoms(&self) -> usize {
        self.atoms.len()
    }

    /// Get the total number of electrons
    pub fn nelec(&self) -> i32 {
        let nuclear_charge: i32 = self.atoms.iter().map(|a| i32::from(a.atomic_number)).sum();
        nuclear_charge - self.charge
    }

    /// Mock water molecule for testing
    pub fn water() -> Self {
        let atoms = vec![
            Atom::new("O", 8, [0.0, 0.0, 0.0]),
            Atom::new("H", 1, [0.0, 1.4, 1.1]),
            Atom::new("H", 1, [0.0, -1.4, 1.1]),
        ];
        Self::new(atoms, 0, 1).with_name("H2O")
    }

    /// Mock ammonia molecule for testing
    pub fn ammonia() -> Self {
        let atoms = vec![
            Atom::new("N", 7, [0.0, 0.0, 0.0]),
            Atom::new("H", 1, [0.0, 1.91, 0.56]),
            Atom::new("H", 1, [1.65, -0.96, 0.56]),
            Atom::new("H", 1, [-1.65, -0.96, 0.56]),
        ];
        Self::new(atoms, 0, 1).with_name("NH3")
    }

    /// Mock carbon monoxide molecule for testing
    pub fn carbon_monoxide() -> Self {
        let atoms = vec![
            Atom::new("C", 6, [0.0, 0.0, 0.0]),
            Atom::new("O", 8, [0.0, 0.0, 2.13]),
        ];
        Self::new(atoms, 0, 1).with_name("CO")
    }

    /// Mock benzene molecule for testing
    pub fn benzene() -> Self {
        let r = 2.63; // C-C bond length in Bohr
        let atoms = vec![
            Atom::new("C", 6, [r, 0.0, 0.0]),
            Atom::new("C", 6, [r / 2.0, r * 0.866, 0.0]),
            Atom::new("C", 6, [-r / 2.0, r * 0.866, 0.0]),
            Atom::new("C", 6, [-r, 0.0, 0.0]),
            Atom::new("C", 6, [-r / 2.0, -r * 0.866, 0.0]),
            Atom::new("C", 6, [r / 2.0, -r * 0.866, 0.0]),
            // Hydrogens
            Atom::new("H", 1, [2.0 * r, 0.0, 0.0]),
            Atom::new("H", 1, [r, r * 1.732, 0.0]),
            Atom::new("H", 1, [-r, r * 1.732, 0.0]),
            Atom::new("H", 1, [-2.0 * r, 0.0, 0.0]),
            Atom::new("H", 1, [-r, -r * 1.732, 0.0]),
            Atom::new("H", 1, [r, -r * 1.732, 0.0]),
        ];
        Self::new(atoms, 0, 1).with_name("C6H6")
    }

    /// Validate the molecule
    pub fn validate(&self) -> Result<()> {
        if self.atoms.is_empty() {
            return Err(
                super::IntegralError::InvalidMolecule("Molecule has no atoms".to_string()).into(),
            );
        }

        if self.nelec() < 0 {
            return Err(super::IntegralError::InvalidMolecule(
                "Molecule has negative number of electrons".to_string(),
            )
            .into());
        }

        if self.multiplicity < 1 {
            return Err(super::IntegralError::InvalidMolecule(
                "Invalid spin multiplicity".to_string(),
            )
            .into());
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_water_molecule() {
        let mol = Molecule::water();
        assert_eq!(mol.natoms(), 3);
        assert_eq!(mol.nelec(), 10);
        assert_eq!(mol.name(), "H2O");
        assert!(mol.validate().is_ok());
    }

    #[test]
    fn test_ammonia_molecule() {
        let mol = Molecule::ammonia();
        assert_eq!(mol.natoms(), 4);
        assert_eq!(mol.nelec(), 10);
        assert_eq!(mol.name(), "NH3");
        assert!(mol.validate().is_ok());
    }

    #[test]
    fn test_co_molecule() {
        let mol = Molecule::carbon_monoxide();
        assert_eq!(mol.natoms(), 2);
        assert_eq!(mol.nelec(), 14);
        assert_eq!(mol.name(), "CO");
        assert!(mol.validate().is_ok());
    }

    #[test]
    fn test_benzene_molecule() {
        let mol = Molecule::benzene();
        assert_eq!(mol.natoms(), 12);
        assert_eq!(mol.nelec(), 42);
        assert_eq!(mol.name(), "C6H6");
        assert!(mol.validate().is_ok());
    }
}
