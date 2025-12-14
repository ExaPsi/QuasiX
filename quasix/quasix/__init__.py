"""
QuasiX: High-performance GW/BSE implementation for molecules and periodic systems

This package provides Python bindings to the QuasiX computational kernel,
enabling efficient quasiparticle and exciton calculations integrated with PySCF.
"""

# Import the Rust extension module
try:
    from .quasix import (
        # Core functions (always available)
        version, metadata, noop_kernel, __version__,
        init_rust_logging,
        # Integral functions
        compute_3center_integrals, compute_2center_integrals,
        get_integral_dimensions,
        # MO transformation and Cholesky functions
        transform_to_mo_basis, compute_cholesky_metric,
        generate_mock_mo_coefficients,
        # Exchange self-energy functions
        compute_exchange_matrix_ri, compute_exchange_diagonal_ri,
        compute_exchange_symmetric,
        # evGW functions (S5-1) - renamed from evgw to evgw_contour_deformation
        evgw as evgw_contour_deformation,
        run_evgw,
        # Dielectric and polarizability functions
        compute_polarizability_p0, compute_polarizability_batch,
        compute_epsilon_inverse,
        # Frequency grid functions (S3-1)
        create_frequency_grid, gauss_legendre_quadrature,
        create_optimized_gl_grid, create_imaginary_axis_grid, create_minimax_grid,
        # Linalg functions
        compute_metric_sqrt, apply_metric_sqrt,
        apply_metric_inv_sqrt, verify_metric_identity,
        # Schema/IO classes
        QuasixData,
        # Optimized I/O functions
        save_pyscf_data, load_pyscf_data, optimize_hdf5_chunks,
        # Analytic continuation classes
        ACConfig, MultipoleModel, PadeModel, AnalyticContinuationFitter
    )
    RUST_EXTENSION_AVAILABLE = True

    # These functions are DISABLED in current Rust bindings - set to None
    compute_correlation_self_energy_cd = None  # DISABLED: API changed
    solve_quasiparticle_equations = None  # DISABLED: QPEquationSolver removed
    compute_z_factors = None  # DISABLED: GW module not registered
    compute_sigma_diag = None  # DISABLED: GW module not registered
    update_polarizability_denominators = None  # DISABLED: GW module not registered
    PolarizabilityBuilder = None  # DISABLED: GW module not registered
    GapStatistics = None  # DISABLED: GW module not registered
    compute_dielectric_function = None  # DISABLED: uses removed DielectricMatrix
    compute_screened_interaction = None  # DISABLED: legacy API
    compute_screened_interaction_batch = None  # DISABLED: legacy API

except ImportError as e:
    # Provide mock implementations for testing
    import warnings
    warnings.warn(f"QuasiX Rust extension not available: {e}. Using mock implementations for testing.")
    RUST_EXTENSION_AVAILABLE = False

    # Mock implementations
    __version__ = "0.1.0"
    
    def version():
        return __version__
    
    def metadata():
        return {"version": __version__, "name": "quasix", "mode": "mock"}
    
    def noop_kernel():
        return {"status": "success", "message": "Mock kernel"}
    
    def init_rust_logging():
        pass
    
    # Provide None for functions that aren't critical
    compute_3center_integrals = None
    compute_2center_integrals = None
    get_integral_dimensions = None
    transform_to_mo_basis = None
    compute_cholesky_metric = None
    generate_mock_mo_coefficients = None
    compute_exchange_matrix_ri = None
    compute_exchange_diagonal_ri = None
    compute_exchange_symmetric = None
    compute_correlation_self_energy_cd = None
    solve_quasiparticle_equations = None
    evgw_contour_deformation = None
    run_evgw = None
    compute_z_factors = None
    compute_sigma_diag = None
    update_polarizability_denominators = None
    PolarizabilityBuilder = None
    GapStatistics = None
    compute_polarizability_p0 = None
    compute_polarizability_batch = None
    compute_dielectric_function = None
    compute_epsilon_inverse = None
    compute_screened_interaction = None
    compute_screened_interaction_batch = None
    create_frequency_grid = None
    gauss_legendre_quadrature = None
    create_optimized_gl_grid = None
    create_imaginary_axis_grid = None
    create_minimax_grid = None
    compute_metric_sqrt = None
    apply_metric_sqrt = None
    apply_metric_inv_sqrt = None
    verify_metric_identity = None
    QuasixData = None
    save_pyscf_data = None
    load_pyscf_data = None
    optimize_hdf5_chunks = None
    ACConfig = None
    MultipoleModel = None
    PadeModel = None
    AnalyticContinuationFitter = None

# Import logging utilities
from .logging import setup_logging, timed_stage, StageTimer, get_logger

# Import monitoring module (S5-3)
try:
    from . import monitoring
    from .monitoring import (
        MonitorConfig,
        EvGWMonitor,
        DataBuffer,
        PlotManager,
        create_monitor,
        monitor_context
    )
except ImportError as e:
    import warnings
    warnings.warn(f"Monitoring module not available: {e}")
    monitoring = None
    MonitorConfig = None
    EvGWMonitor = None
    DataBuffer = None
    PlotManager = None
    create_monitor = None
    monitor_context = None

# Import high-level modules
from . import selfenergy
from . import freq
from . import dielectric
from . import screening
from . import qp
from . import correlation

# Import DF tensor modules for PySCF integration
from .df_tensor import DFTensor
from .df_tensors import (
    DFTensorTransformer,
    transform_mo_3center_numpy,
    validate_against_pyscf,
    benchmark_mo_transformation
)

# Import frequency grid components (S3-1)
from .freq import (
    FrequencyGrid,
    ContourDeformation,
    ACFitter,
    create_frequency_grid as create_frequency_grid_py,
    gauss_legendre_quadrature as gauss_legendre_quadrature_py,
    create_gw_frequency_grid,
    transform_to_imaginary,
    transform_to_matsubara,
)

# Import dielectric components (S3-2, S3-3)
from .dielectric import (
    Polarizability,
    DielectricFunction,
    ScreenedInteraction,
    compute_polarizability_p0 as compute_polarizability_p0_py,
    compute_polarizability_batch as compute_polarizability_batch_py,
    compute_dielectric_function as compute_dielectric_function_py,
    compute_epsilon_inverse as compute_epsilon_inverse_py,
    compute_screened_interaction as compute_screened_interaction_py,
    create_static_dielectric,
    test_kramers_kronig,
)

# Import screening components (S3-3)
from .screening import (
    ScreeningCalculator,
    compute_screened_interaction_symmetric,
    compute_plasmon_pole_model,
)

# Import QP solver components (S3-6)
from .qp import (
    QuasiparticleState,
    QPSolverConfig,
    QPSolver,  # Main S3-6 solver with Rust backend
    LinearizedQPSolver,
    GraphicalQPSolver,
    solve_quasiparticle_equations as solve_quasiparticle_equations_py,
    compute_qp_gap,
    compute_ionization_potential,
    compute_electron_affinity,
    validate_z_factors,
)

# Import Correlation self-energy components (S3-5)
from .correlation import (
    CDParams,
    compute_correlation_self_energy_cd,
    compute_correlation_self_energy_simple,
)

# Import optimized batch processing for correlation self-energy
try:
    from .correlation_batch import (
        BatchProcessor,
        compute_correlation_batch,
    )
except ImportError:
    # Batch processing is optional
    BatchProcessor = None
    compute_correlation_batch = None

# Use Python implementations when Rust ones aren't available
if not RUST_EXTENSION_AVAILABLE:
    create_frequency_grid = create_frequency_grid_py
    gauss_legendre_quadrature = gauss_legendre_quadrature_py
    compute_polarizability_p0 = compute_polarizability_p0_py
    compute_polarizability_batch = compute_polarizability_batch_py
    compute_dielectric_function = compute_dielectric_function_py
    compute_epsilon_inverse = compute_epsilon_inverse_py
    compute_screened_interaction = compute_screened_interaction_py
    solve_quasiparticle_equations = solve_quasiparticle_equations_py

# Import evGW module (S4-1)
try:
    from . import evgw
    from .evgw import EvGW, EvGWParameters as EvGWConfig, EvGWResult
    # IMPORTANT: Use Rust run_evgw that returns dict, not Python evgw.evgw that returns ndarray
    from .quasix import run_evgw as run_evgw_rust
    run_evgw = run_evgw_rust  # Use Rust implementation that returns dict
except ImportError as e:
    # evGW requires additional dependencies
    import warnings
    warnings.warn(f"evGW module not available: {e}")
    evgw = None
    EvGW = None
    EvGWConfig = None
    EvGWResult = None
    run_evgw = None

# Import P0 denominator update module (S5-1)
try:
    from . import polarizability_update
    from .polarizability_update import (
        PolarizabilityUpdater,
        PolarizabilityBuilder as PyPolarizabilityBuilder,
        GapStatistics as PyGapStatistics,
        update_polarizability_denominators as update_p0_denominators_py,
        analyze_gap_evolution,
    )
    # If we don't have the Rust PolarizabilityBuilder, use the one from polarizability_update
    if not RUST_EXTENSION_AVAILABLE or PolarizabilityBuilder is None:
        PolarizabilityBuilder = PyPolarizabilityBuilder
except ImportError as e:
    import warnings
    warnings.warn(f"Polarizability update module not available: {e}")
    polarizability_update = None
    PolarizabilityUpdater = None
    PyPolarizabilityBuilder = None
    PyGapStatistics = None
    update_p0_denominators_py = None
    analyze_gap_evolution = None

# Import spectral analysis module (S5-4)
try:
    from . import spectral
    from .spectral import (
        SpectralAnalyzer,
        SpectralConfig,
        SpectralPeak,
        BroadeningFunctions,
        plot_dos_spectrum,
        compare_spectra,
        analyze_evgw_spectrum,
    )
except ImportError as e:
    import warnings
    warnings.warn(f"Spectral analysis module not available: {e}")
    spectral = None
    SpectralAnalyzer = None
    SpectralConfig = None
    SpectralPeak = None
    BroadeningFunctions = None
    plot_dos_spectrum = None
    compare_spectra = None
    analyze_evgw_spectrum = None

# Import BSE Python driver module (S6)
# Note: 'bse_driver' is used because the Rust extension exports 'bse' as a submodule
try:
    from . import bse_driver
    from .bse_driver import (
        run_bse_tda,
        BSETDAResult,
        BSETDAConfig,
        compute_transition_dipoles,
        compute_absorption_spectrum as compute_bse_absorption_spectrum,
        analyze_exciton as analyze_bse_exciton,
    )
    # Also expose Rust BSE functions
    from .quasix import bse as bse_rust
except ImportError as e:
    import warnings
    warnings.warn(f"BSE module not available: {e}")
    bse_driver = None
    run_bse_tda = None
    BSETDAResult = None
    BSETDAConfig = None
    compute_transition_dipoles = None
    compute_bse_absorption_spectrum = None
    analyze_bse_exciton = None
    bse_rust = None

# Create linalg namespace for cleaner API
class LinalgNamespace:
    """Linear algebra functions for metric tensor operations."""
    compute_metric_sqrt = compute_metric_sqrt
    apply_metric_sqrt = apply_metric_sqrt
    apply_metric_inv_sqrt = apply_metric_inv_sqrt
    verify_metric_identity = verify_metric_identity

linalg = LinalgNamespace()

# Initialize logging if environment variables are set
import os
if os.environ.get("QUASIX_LOG_LEVEL") or os.environ.get("QUASIX_LOG"):
    setup_logging()

# Public API
__all__ = [
    "__version__",
    "version",
    "metadata",
    "noop_kernel",
    "compute_3center_integrals",
    "compute_2center_integrals",
    "get_integral_dimensions",
    "init_rust_logging",
    "transform_to_mo_basis",
    "compute_cholesky_metric",
    "generate_mock_mo_coefficients",
    # Exchange self-energy functions
    "compute_exchange_matrix_ri",
    "compute_exchange_diagonal_ri",
    "compute_exchange_symmetric",
    # Correlation self-energy functions
    "compute_correlation_self_energy_cd",
    # Quasiparticle solver functions
    "solve_quasiparticle_equations",
    # evGW functions (S4-1)
    "evgw_contour_deformation",
    # GW module functions
    "run_evgw",
    "compute_z_factors",
    "compute_sigma_diag",
    "update_polarizability_denominators",
    "PolarizabilityBuilder",
    "GapStatistics",
    # Dielectric and polarizability functions
    "compute_polarizability_p0",
    "compute_polarizability_batch",
    "compute_dielectric_function",
    "compute_epsilon_inverse",
    # Screened interaction functions
    "compute_screened_interaction",
    "compute_screened_interaction_batch",
    # Frequency grid functions (S3-1)
    "create_frequency_grid",
    "gauss_legendre_quadrature",
    "create_optimized_gl_grid",
    "create_imaginary_axis_grid",
    "create_minimax_grid",
    "FrequencyGrid",
    "ContourDeformation",
    "ACFitter",
    # Analytic continuation classes
    "ACConfig",
    "MultipoleModel",
    "PadeModel",
    "AnalyticContinuationFitter",
    "create_gw_frequency_grid",
    "transform_to_imaginary",
    "transform_to_matsubara",
    # Dielectric components (S3-2, S3-3)
    "dielectric",
    "Polarizability",
    "DielectricFunction",
    "ScreenedInteraction",
    "create_static_dielectric",
    "test_kramers_kronig",
    # Screening module (S3-3)
    "screening",
    "ScreeningCalculator",
    "compute_screened_interaction_symmetric",
    "compute_plasmon_pole_model",
    # QP solver module (S3-6)
    "qp",
    "QuasiparticleState",
    "QPSolverConfig",
    "QPSolver",
    "LinearizedQPSolver",
    "GraphicalQPSolver",
    "compute_qp_gap",
    "compute_ionization_potential",
    "compute_electron_affinity",
    "validate_z_factors",
    "setup_logging",
    "timed_stage",
    "StageTimer",
    "get_logger",
    "linalg",  # Add linalg namespace to public API
    "QuasixData",  # Add schema/IO class to public API
    "save_pyscf_data",  # S2-4: Optimized PySCF data save
    "load_pyscf_data",  # S2-4: Optimized PySCF data load
    "optimize_hdf5_chunks",  # S2-4: HDF5 chunk optimization
    "selfenergy",  # Add selfenergy module to public API
    "freq",  # Add frequency module to public API
    "DFTensor",  # Add DFTensor class for PySCF integration
    "DFTensorTransformer",  # S2-2: DF tensor transformer
    "transform_mo_3center_numpy",  # S2-2: NumPy reference implementation
    "validate_against_pyscf",  # S2-2: PySCF validation utility
    "benchmark_mo_transformation",  # S2-2: Performance benchmarking
    # evGW module (S4-1)
    "evgw",
    "EvGW",
    "EvGWConfig",
    "EvGWResult",
    "run_evgw",
    # Correlation self-energy module (S3-5)
    "correlation",
    "CDParams",
    "compute_correlation_self_energy_cd",
    "compute_correlation_self_energy_simple",
    "BatchProcessor",
    "compute_correlation_batch",
    # P0 denominator update module (S5-1)
    "polarizability_update",
    "PolarizabilityUpdater",
    "PyGapStatistics",
    "update_p0_denominators_py",
    "analyze_gap_evolution",
    # Monitoring module (S5-3)
    "monitoring",
    "MonitorConfig",
    "EvGWMonitor",
    "DataBuffer",
    "PlotManager",
    "create_monitor",
    "monitor_context",
    # Spectral analysis module (S5-4)
    "spectral",
    "SpectralAnalyzer",
    "SpectralConfig",
    "SpectralPeak",
    "BroadeningFunctions",
    "plot_dos_spectrum",
    "compare_spectra",
    "analyze_evgw_spectrum",
    # BSE module (S6)
    "bse_driver",
    "bse_rust",
    "run_bse_tda",
    "BSETDAResult",
    "BSETDAConfig",
    "compute_transition_dipoles",
    "compute_bse_absorption_spectrum",
    "analyze_bse_exciton",
]