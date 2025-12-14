//! QuasiX data schema definitions for JSON and HDF5 I/O
//!
//! This module defines the complete data structures for storing
//! GW/BSE calculations with full provenance tracking.
#![warn(clippy::all, clippy::pedantic, clippy::perf)]
#![warn(missing_docs)]
#![allow(clippy::module_name_repetitions)] // Allow GWParams, BSEParams etc.
#![allow(clippy::many_single_char_names)] // Mathematical notation in comments

use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Complete QuasiX calculation data with provenance
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QuasixData {
    /// Metadata and provenance information
    pub metadata: Metadata,
    /// Input molecular system
    pub molecule: MoleculeInput,
    /// Basis set information
    pub basis: BasisSetData,
    /// Calculation parameters
    pub parameters: CalculationParams,
    /// Results from the calculation
    pub results: CalculationResults,
}

/// Metadata for provenance tracking
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Metadata {
    /// QuasiX version
    pub quasix_version: String,
    /// Calculation timestamp
    pub timestamp: DateTime<Utc>,
    /// User who ran the calculation
    pub user: Option<String>,
    /// Host machine
    pub hostname: Option<String>,
    /// Git commit hash (if available)
    pub git_commit: Option<String>,
    /// Calculation UUID for tracking
    pub uuid: String,
    /// Parent calculation UUID (for restarts)
    pub parent_uuid: Option<String>,
    /// Custom metadata
    pub custom: HashMap<String, serde_json::Value>,
}

/// Input molecular structure
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MoleculeInput {
    /// Number of atoms
    pub natoms: usize,
    /// Atomic symbols (e.g., `["H", "O", "H"]`)
    pub symbols: Vec<String>,
    /// Atomic numbers
    pub atomic_numbers: Vec<u32>,
    /// Cartesian coordinates in Bohr
    pub coordinates: Vec<[f64; 3]>,
    /// Total charge
    pub charge: i32,
    /// Spin multiplicity (2S+1)
    pub multiplicity: u32,
    /// Point group symmetry
    pub symmetry: Option<String>,
}

/// Basis set information
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BasisSetData {
    /// AO basis set name (e.g., "def2-SVP")
    pub ao_basis: String,
    /// Number of AO basis functions
    pub n_ao: usize,
    /// Auxiliary basis set name (e.g., "def2-SVP-JKFIT")
    pub aux_basis: String,
    /// Number of auxiliary basis functions
    pub n_aux: usize,
    /// Detailed basis function info (optional)
    pub functions: Option<Vec<BasisFunction>>,
}

/// Individual basis function
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BasisFunction {
    /// Atom index (0-based)
    pub center: usize,
    /// Angular momentum (S=0, P=1, D=2, etc.)
    pub angular_momentum: u32,
    /// Gaussian exponents
    pub exponents: Vec<f64>,
    /// Contraction coefficients
    pub coefficients: Vec<f64>,
}

/// Calculation parameters
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CalculationParams {
    /// Type of calculation
    pub calculation_type: CalculationType,
    /// GW-specific parameters
    pub gw_params: Option<GWParams>,
    /// BSE-specific parameters
    pub bse_params: Option<BSEParams>,
    /// Convergence criteria
    pub convergence: ConvergenceParams,
    /// Frequency grid parameters
    pub frequency: FrequencyParams,
}

/// Type of calculation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum CalculationType {
    /// One-shot G0W0
    G0W0,
    /// Eigenvalue self-consistent GW
    EvGW,
    /// Self-consistent GW
    ScGW,
    /// BSE on top of GW
    BSE,
    /// Combined GW+BSE
    GWBSE,
}

/// GW calculation parameters
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GWParams {
    /// Starting point (HF, DFT functional)
    pub starting_point: String,
    /// Number of states to compute
    pub n_states: usize,
    /// Use contour deformation
    pub contour_deformation: bool,
    /// Use analytic continuation
    pub analytic_continuation: bool,
    /// Linearization scheme
    pub linearization: String,
    /// Maximum iterations for self-consistency
    pub max_iter: usize,
}

/// BSE calculation parameters
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BSEParams {
    /// Use Tamm-Dancoff approximation
    pub tda: bool,
    /// Number of excitations to compute
    pub n_excitations: usize,
    /// Include triplet excitations
    pub triplets: bool,
    /// Kernel approximation
    pub kernel: String,
}

/// Convergence parameters
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConvergenceParams {
    /// Energy convergence threshold
    pub energy_tol: f64,
    /// Density convergence threshold
    pub density_tol: f64,
    /// Maximum iterations
    pub max_iterations: usize,
    /// Use DIIS acceleration
    pub use_diis: bool,
    /// DIIS space size
    pub diis_space: usize,
}

/// Frequency grid parameters
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FrequencyParams {
    /// Grid type (`GaussLegendre`, `GaussLaguerre`, etc.)
    pub grid_type: String,
    /// Number of frequency points
    pub n_points: usize,
    /// Maximum frequency (for real axis)
    pub omega_max: Option<f64>,
    /// Broadening parameter
    pub eta: f64,
}

/// Complete calculation results
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CalculationResults {
    /// Reference (HF/DFT) results
    pub reference: ReferenceResults,
    /// GW results (if computed)
    pub gw: Option<GWResults>,
    /// BSE results (if computed)
    pub bse: Option<BSEResults>,
    /// Timing information
    pub timings: TimingInfo,
    /// Convergence history
    pub convergence_history: Vec<ConvergenceStep>,
}

/// Reference calculation results
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReferenceResults {
    /// Total energy
    pub energy: f64,
    /// Orbital energies
    pub orbital_energies: Vec<f64>,
    /// Occupation numbers
    pub occupations: Vec<f64>,
    /// HOMO index
    pub homo: usize,
    /// LUMO index
    pub lumo: usize,
    /// MO coefficients (optional, can be large)
    pub mo_coefficients: Option<Vec<Vec<f64>>>,
}

/// GW calculation results
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GWResults {
    /// Quasiparticle energies
    pub qp_energies: Vec<f64>,
    /// Z-factors (renormalization)
    pub z_factors: Vec<f64>,
    /// Self-energy values (optional)
    pub self_energy: Option<SelfEnergyData>,
    /// Spectral functions (optional)
    pub spectral_functions: Option<Vec<Vec<f64>>>,
    /// Converged flag
    pub converged: bool,
    /// Number of iterations
    pub iterations: usize,
}

/// Self-energy data
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SelfEnergyData {
    /// Exchange self-energy (diagonal)
    pub sigma_x: Vec<f64>,
    /// Correlation self-energy (diagonal, real part)
    pub sigma_c_real: Vec<f64>,
    /// Correlation self-energy (diagonal, imaginary part)
    pub sigma_c_imag: Vec<f64>,
    /// Full matrix elements (optional)
    pub matrix_elements: Option<Vec<Vec<[f64; 2]>>>,
}

/// BSE calculation results
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BSEResults {
    /// Excitation energies
    pub excitation_energies: Vec<f64>,
    /// Oscillator strengths
    pub oscillator_strengths: Vec<f64>,
    /// Transition dipole moments
    pub dipole_moments: Vec<[f64; 3]>,
    /// Exciton wavefunctions (optional)
    pub wavefunctions: Option<Vec<ExcitonWavefunction>>,
    /// Absorption spectrum
    pub absorption_spectrum: Option<SpectrumData>,
}

/// Exciton wavefunction
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExcitonWavefunction {
    /// Excitation index
    pub index: usize,
    /// Dominant configurations [(hole, particle, coefficient)]
    pub configurations: Vec<(usize, usize, f64)>,
    /// Participation ratio
    pub participation_ratio: f64,
    /// Exciton binding energy
    pub binding_energy: f64,
}

/// Spectrum data
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpectrumData {
    /// Frequency points
    pub omega: Vec<f64>,
    /// Absorption intensity
    pub absorption: Vec<f64>,
    /// Real part of dielectric function
    pub epsilon_real: Option<Vec<f64>>,
    /// Imaginary part of dielectric function
    pub epsilon_imag: Option<Vec<f64>>,
}

/// Timing information
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TimingInfo {
    /// Total wall time in seconds
    pub total_time: f64,
    /// Breakdown by component
    pub components: HashMap<String, f64>,
    /// Peak memory usage in MB
    pub peak_memory: Option<f64>,
}

/// Convergence step information
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConvergenceStep {
    /// Iteration number
    pub iteration: usize,
    /// Energy change
    pub energy_change: f64,
    /// Density change
    pub density_change: f64,
    /// Maximum residual
    pub max_residual: f64,
    /// Wall time for this iteration
    pub time: f64,
}

// HDF5 group structure constants

/// HDF5 root group
pub const HDF5_ROOT: &str = "/";
/// HDF5 metadata group
pub const HDF5_METADATA: &str = "/metadata";
/// HDF5 molecule group
pub const HDF5_MOLECULE: &str = "/molecule";
/// HDF5 basis set group
pub const HDF5_BASIS: &str = "/basis";
/// HDF5 parameters group
pub const HDF5_PARAMS: &str = "/parameters";
/// HDF5 results group
pub const HDF5_RESULTS: &str = "/results";
/// HDF5 reference results group
pub const HDF5_REFERENCE: &str = "/results/reference";
/// HDF5 GW results group
pub const HDF5_GW: &str = "/results/gw";
/// HDF5 BSE results group
pub const HDF5_BSE: &str = "/results/bse";
/// HDF5 timing information group
pub const HDF5_TIMINGS: &str = "/results/timings";
/// HDF5 convergence history group
pub const HDF5_CONVERGENCE: &str = "/results/convergence";

impl QuasixData {
    /// Create a new `QuasixData` instance with default metadata
    #[must_use]
    pub fn new(molecule: MoleculeInput, basis: BasisSetData, params: CalculationParams) -> Self {
        use uuid::Uuid;

        Self {
            metadata: Metadata {
                quasix_version: env!("CARGO_PKG_VERSION").to_string(),
                timestamp: Utc::now(),
                user: std::env::var("USER").ok(),
                hostname: hostname::get().ok().and_then(|h| h.into_string().ok()),
                git_commit: None, // Could use git2 crate to get this
                uuid: Uuid::new_v4().to_string(),
                parent_uuid: None,
                custom: HashMap::new(),
            },
            molecule,
            basis,
            parameters: params,
            results: CalculationResults {
                reference: ReferenceResults {
                    energy: 0.0,
                    orbital_energies: Vec::new(),
                    occupations: Vec::new(),
                    homo: 0,
                    lumo: 0,
                    mo_coefficients: None,
                },
                gw: None,
                bse: None,
                timings: TimingInfo {
                    total_time: 0.0,
                    components: HashMap::new(),
                    peak_memory: None,
                },
                convergence_history: Vec::new(),
            },
        }
    }

    /// Save to JSON file.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be written or serialization fails
    pub fn to_json(&self, path: &std::path::Path) -> Result<(), Box<dyn std::error::Error>> {
        let json = serde_json::to_string_pretty(self)?;
        std::fs::write(path, json)?;
        Ok(())
    }

    /// Load from JSON file.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be read or deserialization fails
    pub fn from_json(path: &std::path::Path) -> Result<Self, Box<dyn std::error::Error>> {
        let json = std::fs::read_to_string(path)?;
        let data = serde_json::from_str(&json)?;
        Ok(data)
    }
}
