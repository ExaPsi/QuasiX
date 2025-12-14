//! QuasiX Core - High-performance quantum chemistry library
//!
//! This library provides production-quality implementations of:
//! - evGW and scGW quasiparticle methods
//! - Bethe-Salpeter equation (BSE) for excitations
//! - Resolution-of-identity (RI) / density fitting
//! - Optimized integral evaluations with SIMD and parallelization

// ============================================================================
// STRICT CLIPPY CONFIGURATION FOR ZERO-WARNING COMPILATION
// ============================================================================

// Enable all clippy lints for maximum code quality
#![warn(clippy::all, clippy::pedantic, clippy::perf)]
// Scientific computing specific allowances
#![allow(clippy::many_single_char_names)] // Mathematical notation (i, j, k, p, q, r, s common in QC)
#![allow(clippy::similar_names)] // p, q, pp, qq are standard in quantum chemistry
#![allow(clippy::cast_precision_loss)] // Array indexing to f64 is standard in numerical code
#![allow(clippy::cast_possible_truncation)] // Safe when bounds are known (e.g., timing to u64)
#![allow(clippy::cast_sign_loss)] // Safe when values are guaranteed positive
#![allow(clippy::cast_possible_wrap)] // Safe in controlled numerical contexts

// Module and documentation allowances
#![allow(clippy::module_name_repetitions)] // Scientific code often has repetitive names for clarity
#![allow(clippy::must_use_candidate)] // Not all functions need #[must_use]
#![allow(clippy::missing_errors_doc)] // TODO: Add comprehensive error documentation
#![allow(clippy::doc_markdown)] // Scientific terms don't need backticks
#![allow(clippy::uninlined_format_args)] // Format string readability preference
#![allow(clippy::missing_panics_doc)] // Would require documenting all potential panics
#![allow(clippy::missing_docs_in_private_items)] // Private implementation details
#![allow(missing_docs)] // Temporarily allow for rapid development

// Code style allowances for scientific computing
#![allow(clippy::items_after_statements)] // Constants in functions improve readability
#![allow(clippy::wildcard_imports)] // Common for schema imports
#![allow(clippy::unused_self)] // Trait implementations may not use self
#![allow(clippy::unnecessary_wraps)] // Sometimes used for API consistency
#![allow(clippy::assertions_on_constants)] // Test assertions may check constant conditions
#![allow(clippy::field_reassign_with_default)] // Sometimes clearer to reassign individual fields
#![allow(clippy::redundant_closure_for_method_calls)] // Sometimes clearer with explicit closure
#![allow(clippy::manual_midpoint)] // Scientific code often needs explicit formulas
#![allow(clippy::single_match_else)] // Sometimes clearer than if let
#![allow(clippy::implicit_clone)] // Clone is acceptable for small types
#![allow(clippy::map_unwrap_or)] // Standard pattern in option/result handling
#![allow(clippy::struct_excessive_bools)] // Configuration structs may have many bools
#![allow(clippy::ptr_cast_constness)] // FFI code needs const casting
#![allow(clippy::used_underscore_binding)] // Sometimes needed for clarity
#![allow(clippy::as_ptr_cast_mut)] // Needed for FFI and unsafe optimizations
#![allow(clippy::needless_pass_by_value)] // Sometimes clearer API design
#![allow(clippy::ptr_as_ptr)] // Raw pointer casts needed for FFI
#![allow(clippy::no_effect_underscore_binding)] // Sometimes used for explicit type inference
#![allow(clippy::manual_let_else)] // Sometimes explicit match is clearer
#![allow(clippy::manual_clamp)] // Explicit min/max chains are clearer in numerical code
#![allow(clippy::manual_is_multiple_of)] // Explicit modulo is clearer in performance code
#![allow(clippy::manual_range_contains)] // Explicit range checks are clearer with tolerances
#![allow(clippy::bool_to_int_with_if)] // Explicit conditional is clearer than From conversion
#![allow(clippy::cloned_instead_of_copied)] // Iterator methods work with both
#![allow(clippy::assigning_clones)] // Assignment to existing variable is clear
#![allow(clippy::needless_range_loop)] // Sometimes index-based loops are clearer for FFI
#![allow(clippy::default_trait_access)] // Default::default() is acceptable
#![allow(clippy::match_same_arms)] // Sometimes clearer to have explicit matching arms
#![allow(clippy::diverging_sub_expression)] // Used intentionally in unimplemented! cases
#![allow(clippy::struct_field_names)] // Scientific structs often have prefixed field names
#![allow(clippy::doc_overindented_list_items)] // Documentation style preference
#![allow(clippy::manual_memcpy)] // Explicit loops may be clearer for SIMD-friendly code
#![allow(clippy::cast_lossless)] // as casts are acceptable when types are known
#![allow(clippy::op_ref)] // Reference operand is clearer in some matrix operations
#![allow(clippy::assign_op_pattern)] // Explicit assignment is clearer in some scientific contexts
#![allow(clippy::elidable_lifetime_names)] // Explicit lifetimes are clearer in PyO3 bindings

// Performance-related allowances
#![allow(clippy::inline_always)] // We manually optimize hot paths
#![allow(clippy::too_many_arguments)] // Complex scientific functions may need many parameters
#![allow(clippy::type_complexity)] // Complex types are sometimes necessary in generic code
#![allow(clippy::too_many_lines)] // Complex scientific algorithms require long functions
#![allow(clippy::format_push_string)] // Format strings for report generation are acceptable
#![allow(clippy::unreadable_literal)] // Physical constants like 27.211386 should be readable as-is

// Test-specific allowances
#![allow(clippy::float_cmp)] // Exact comparison is intentional in tests
#![allow(clippy::return_self_not_must_use)]
// Builder pattern doesn't always need must_use

// Raw pointer allowances (for FFI and performance-critical code)
#![allow(clippy::borrow_as_ptr)] // Sometimes needed for FFI or unsafe optimizations
//!
//! ## Module Organization
//!
//! - [`df`]: Density fitting and resolution-of-identity
//! - [`freq`]: Frequency grids and transformations
//! - [`dielectric`]: Polarizability and screening calculations
//! - [`selfenergy`]: Exchange and correlation self-energy
//! - [`qp`]: Quasiparticle equation solvers
//! - [`bse`]: Bethe-Salpeter equation solver
//! - [`io`]: Input/output and data schemas
//! - [`linalg`]: Linear algebra utilities
//! - [`pbc`]: Periodic boundary conditions

/// Version of the QuasiX core library
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

// Module declarations
pub mod analysis;
pub mod benchmarks;
pub mod blas_utils;
pub mod bse;
pub mod common;
pub mod df;
pub mod dielectric;
pub mod freq;
pub mod gw;
pub mod integrals;
pub mod io;
pub mod linalg;
pub mod logging;
pub mod pbc;
pub mod qp;
pub mod quality_check;
pub mod selfenergy;
pub mod validation;

// NOTE: Fallback module removed during G₀W₀ re-implementation cleanup
// See: docs/G0W0/README.md

// Python bindings module (conditional compilation)
#[cfg(feature = "python")]
pub mod python;

// Re-export commonly used types
pub use analysis::{BroadeningParams, SpectralAnalyzer, SpectralData, SpectralFunction};
pub use bse::{
    BSEKernelConfig, BSETDADriver, BSETDAKernel, DavidsonConfig, OpticalProperties, SpinType,
};
pub use common::{version, QuasixError, Result};
pub use df::{DFBuilder, RIBasis};
pub use dielectric::{DielectricFunction, PolarizabilityRI, ScreenedCoulomb};
pub use freq::{ACFitter, ContourDeformation, FrequencyGrid, GridType};
// NOTE: Fallback and QualityMetrics types removed during G₀W₀ cleanup
// See: docs/G0W0/README.md for re-implementation plan
pub use integrals::{
    compute_2center_integrals, compute_3center_integrals, BasisSet, IntegralEngine, Molecule,
};
pub use io::{BSEResults, CheckpointManager, GWResults, MolecularData};
pub use linalg::{LinalgError, MetricSqrt};
pub use pbc::{CoulombTruncation, KMesh, PeriodicBSE, PeriodicGW};
pub use qp::numa_alloc::{NumaAllocator, NumaPolicy, NumaTopology};
pub use qp::{ConvergenceCriteria, EvGWDriver, QPSolver, ScGWDriver};
pub use selfenergy::{CorrelationSelfEnergy, GWSelfEnergy};
// Re-export the new RI-based exchange self-energy as the primary interface
pub use selfenergy::ExchangeSelfEnergyRI;

/// Library metadata
pub const AUTHORS: &str = env!("CARGO_PKG_AUTHORS");

// PyO3 module entry point (only when Python feature is enabled)
#[cfg(feature = "python")]
use pyo3::prelude::*;

#[cfg(feature = "python")]
#[pymodule]
fn quasix_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();
    python::init_python_module(py, m)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_version() {
        assert_eq!(VERSION, "0.6.0");
        assert_eq!(version(), "0.6.0");
    }

    #[test]
    fn test_module_imports() {
        // Test that all modules are accessible
        let _ = common::constants::HARTREE_TO_EV;
        // DFBuilder requires Molecule + BasisSet - tested in df module
        let _ = FrequencyGrid::new(16, GridType::GaussLegendre);
        let _ = PolarizabilityRI::new(5, 10, 50);
        let _ = ExchangeSelfEnergyRI::new(20, 5, 60); // naux = 3*nbasis
        let _ = QPSolver::new(10);
        // BSETDAKernel requires full input tensors - tested in bse module
        let _ = MolecularData::new(3, 10, 30);
        // MatrixOps is a module, not a value
    }
}
