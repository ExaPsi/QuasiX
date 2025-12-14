#!/usr/bin/env python3
"""GW Driver implementation for benchmark calculations.

This module provides a high-level interface to QuasiX evGW calculations,
specifically designed for GW100 benchmarks and validation.
"""

import numpy as np
import time
from typing import Dict, Optional, Any, Tuple
import warnings
from dataclasses import dataclass, field

try:
    from pyscf import gto, scf, df
    HAS_PYSCF = True
except ImportError:
    HAS_PYSCF = False
    warnings.warn("PySCF not available")

# Import QuasiX modules
try:
    # Try importing QuasiX module
    import sys
    from pathlib import Path
    # Add parent directory to path for the quasix package
    parent_dir = Path(__file__).parent.parent
    if str(parent_dir) not in sys.path:
        sys.path.insert(0, str(parent_dir))

    # Import the main QuasiX module
    import quasix
    from quasix import EvGW, EvGWConfig, EvGWResult
    from quasix import monitoring, spectral

    # Check if we have the run_evgw function (Rust extension)
    # First check in the quasix.quasix submodule (where Rust functions are)
    if hasattr(quasix, 'quasix') and hasattr(quasix.quasix, 'evgw'):
        RUST_AVAILABLE = True
        rust_module = quasix.quasix
    elif hasattr(quasix, 'run_evgw'):
        RUST_AVAILABLE = True
        rust_module = quasix
    else:
        RUST_AVAILABLE = False
        rust_module = None
        warnings.warn("QuasiX Rust extension not available")

    QUASIX_AVAILABLE = True
    # Use EvGWConfig instead of EvGWParameters
    EvGWParameters = EvGWConfig

except ImportError as e:
    QUASIX_AVAILABLE = False
    RUST_AVAILABLE = False
    rust_module = None
    warnings.warn(f"QuasiX modules not available: {e}, using mock implementation")
    # Define mock classes
    EvGW = None
    EvGWParameters = None
    EvGWConfig = None
    EvGWResult = None
    monitoring = None
    spectral = None

# Physical constants
HA_TO_EV = 27.211386245988


@dataclass
class GWConfig:
    """Configuration for GW calculations."""
    auxbasis: str = 'cc-pvdz-jkfit'
    freq_int: str = 'cd'  # Contour deformation
    nfreq: int = 48
    eta: float = 0.01
    max_cycle: int = 10
    conv_tol: float = 1e-5
    conv_tol_z: float = 1e-4
    damping: float = 0.3
    damping_dynamic: bool = True
    diis: bool = True
    diis_space: int = 5
    validate_z: bool = True
    monitor: bool = False
    verbose: int = 0
    # QP solver settings
    qp_solver: str = 'linearized'  # 'linearized' or 'newton'
    newton_max_iterations: int = 50
    newton_energy_tolerance: float = 1e-8


@dataclass
class GWResult:
    """Results from GW calculation."""
    converged: bool
    n_iterations: int
    qp_energies: np.ndarray  # Hartree
    z_factors: np.ndarray
    sigma_x: np.ndarray  # Exchange self-energy (diagonal)
    sigma_c: np.ndarray  # Correlation self-energy (diagonal)
    homo_index: int
    lumo_index: int
    ip: float  # eV
    ea: Optional[float]  # eV
    gap: float  # eV
    elapsed_seconds: float
    convergence_data: Optional[Dict] = None
    spectral_data: Optional[Dict] = None

    def get_homo_energy(self) -> float:
        """Get HOMO energy in eV."""
        return self.qp_energies[self.homo_index] * HA_TO_EV

    def get_lumo_energy(self) -> float:
        """Get LUMO energy in eV."""
        return self.qp_energies[self.lumo_index] * HA_TO_EV

    def get_qp_corrections(self, mo_energy: np.ndarray) -> np.ndarray:
        """Get QP corrections (QP - HF) in eV."""
        return (self.qp_energies - mo_energy) * HA_TO_EV


class GWDriver:
    """High-level driver for GW calculations.

    This class provides a unified interface to QuasiX evGW calculations,
    specifically designed for benchmarks and validation.
    """

    def __init__(self,
                 mo_energy: np.ndarray,
                 mo_coeff: np.ndarray,
                 n_occ: int,
                 mol: Optional[Any] = None,
                 mf: Optional[Any] = None,  # Add optional mf parameter
                 auxbasis: str = 'cc-pvdz-jkfit',
                 config: Optional[GWConfig] = None):
        """Initialize GW driver.

        Args:
            mo_energy: Molecular orbital energies (Hartree)
            mo_coeff: Molecular orbital coefficients
            n_occ: Number of occupied orbitals
            mol: PySCF molecule object (optional)
            mf: PySCF mean-field object (optional, for getting vxc)
            auxbasis: Auxiliary basis set for RI
            config: GW configuration parameters
        """
        self.mo_energy = np.asarray(mo_energy)
        self.mo_coeff = np.asarray(mo_coeff)
        self.n_occ = n_occ
        self.mol = mol
        self.mf_external = mf  # Store external mf if provided
        self.auxbasis = auxbasis
        self.config = config or GWConfig(auxbasis=auxbasis)

        # Validate input
        self._validate_input()

        # Initialize GW object if we have PySCF
        self.gw_obj = None
        if HAS_PYSCF and mol is not None:
            self._initialize_gw()

    def _validate_input(self):
        """Validate input parameters."""
        n_mo = len(self.mo_energy)

        if self.mo_coeff.shape[1] != n_mo:
            raise ValueError(f"Inconsistent MO dimensions: energy {n_mo}, coeff {self.mo_coeff.shape}")

        if self.n_occ <= 0 or self.n_occ >= n_mo:
            raise ValueError(f"Invalid n_occ: {self.n_occ} (must be in [1, {n_mo-1}])")

        # Check orbital energies are ordered (occupied should be increasing, virtual should be increasing)
        # Allow for small numerical variations (0.1 Ha tolerance)
        occ_diff = np.diff(self.mo_energy[:self.n_occ])
        virt_diff = np.diff(self.mo_energy[self.n_occ:])

        # Only warn if there's a significant disorder (> 0.1 Ha)
        if len(occ_diff) > 0 and np.any(occ_diff < -0.1):
            warnings.warn("Occupied orbital energies may not be properly ordered")
        if len(virt_diff) > 0 and np.any(virt_diff < -0.1):
            warnings.warn("Virtual orbital energies may not be properly ordered")

    def _initialize_gw(self):
        """Initialize evGW object with PySCF."""
        if not HAS_PYSCF:
            warnings.warn("PySCF is required for full GW calculations")
            return

        # Create a mock SCF object if needed
        class MockSCF:
            def __init__(self, mol, mo_energy, mo_coeff, mo_occ):
                self.mol = mol
                self.mo_energy = mo_energy
                self.mo_coeff = mo_coeff
                self.mo_occ = mo_occ
                self.e_tot = 0.0
                self.converged = True
                self.verbose = 0
                self.stdout = None  # Add stdout attribute for PySCF compatibility

        # Use external mf if provided, otherwise create mock
        if hasattr(self, 'mf_external') and self.mf_external is not None:
            self.mf = self.mf_external
            # Get mo_occ from the external mf
            if hasattr(self.mf, 'mo_occ'):
                self.mo_occ = self.mf.mo_occ
            else:
                mo_occ = np.zeros(len(self.mo_energy))
                mo_occ[:self.n_occ] = 2.0  # Closed-shell assumption
                self.mo_occ = mo_occ
        else:
            # Create mo_occ array
            mo_occ = np.zeros(len(self.mo_energy))
            mo_occ[:self.n_occ] = 2.0  # Closed-shell assumption
            self.mf = MockSCF(self.mol, self.mo_energy, self.mo_coeff, mo_occ)
            self.mo_occ = mo_occ

        # Initialize evGW object if available
        if QUASIX_AVAILABLE and EvGWConfig and EvGW:
            # Initialize EvGW with correct signature: EvGW(mf, auxbasis, freq_int)
            self.gw_obj = EvGW(self.mf, auxbasis=self.config.auxbasis, freq_int=self.config.freq_int)

            # Set parameters after initialization
            self.gw_obj.max_cycle = self.config.max_cycle
            self.gw_obj.conv_tol = self.config.conv_tol
            self.gw_obj.conv_tol_z = self.config.conv_tol_z
            self.gw_obj.damping = self.config.damping
            self.gw_obj.damping_dynamic = self.config.damping_dynamic
            self.gw_obj.diis = self.config.diis
            self.gw_obj.diis_space = self.config.diis_space
            self.gw_obj.nfreq = self.config.nfreq
            self.gw_obj.eta = self.config.eta
            self.gw_obj.verbose = self.config.verbose
        else:
            self.gw_obj = None

        # Store DF tensors for later use
        self.df_tensors = None

    def extract_vxc_exact(self):
        """Extract exact exchange-correlation potential without DF approximation.

        This is CRITICAL for accurate G0W0 calculations. Using DF-approximated vxc
        causes exact cancellation with DF-approximated exchange self-energy, leading
        to errors > 1 eV in quasiparticle energies.

        Returns:
            np.ndarray: Diagonal of v_xc in MO basis (Ha)

        Raises:
            RuntimeError: If v_xc cannot be computed exactly
        """
        # Import scf module at the top of the method
        from pyscf import scf

        # For HF calculations
        if hasattr(self, 'mf') and self.mf is not None:
            mf = self.mf

            # Check if this is HF or DFT
            # UHF and ROHF are in scf module, not scf.hf
            is_hf = False
            try:
                hf_types = []
                if hasattr(scf, 'hf') and hasattr(scf.hf, 'RHF'):
                    hf_types.append(scf.hf.RHF)
                if hasattr(scf, 'RHF'):
                    hf_types.append(scf.RHF)
                if hasattr(scf, 'UHF'):
                    hf_types.append(scf.UHF)
                if hasattr(scf, 'ROHF'):
                    hf_types.append(scf.ROHF)
                if hf_types:
                    is_hf = isinstance(mf, tuple(hf_types))
            except (TypeError, AttributeError):
                is_hf = False
            # DFT classes are in dft.rks submodule
            is_dft = False
            try:
                from pyscf.dft import rks, uks, roks
                is_dft = isinstance(mf, (rks.RKS, uks.UKS, roks.ROKS))
            except (ImportError, AttributeError):
                # Try alternate import
                try:
                    from pyscf import dft
                    if hasattr(dft, 'RKS'):
                        is_dft = isinstance(mf, dft.RKS)
                except ImportError:
                    pass

            # CRITICAL: Check DFT FIRST since dft.RKS inherits from scf.RHF
            # (is_hf will be True for DFT objects due to inheritance)
            if is_dft:
                # DFT case: v_xc from XC functional
                try:
                    if self.config.verbose >= 1:
                        print("Extracting v_xc from DFT XC functional...")

                    # For DFT, we need the XC potential from the functional
                    dm = mf.make_rdm1()

                    # Get the XC potential
                    ni = mf._numint
                    xc_code = mf.xc

                    # Evaluate XC potential on grid
                    # nr_vxc returns (exc, ecoul, vxc) - we need index [2] for vxc matrix!
                    vxc = ni.nr_vxc(mf.mol, mf.grids, xc_code, dm)[2]

                    # Transform to MO basis
                    vxc_mo = np.dot(mf.mo_coeff.T, np.dot(vxc, mf.mo_coeff))
                    vxc_diag = np.diag(vxc_mo)

                    if self.config.verbose >= 1:
                        print(f"✓ Extracted DFT v_xc from {xc_code} functional")
                        print(f"  HOMO v_xc = {vxc_diag[self.n_occ-1]:.6f} Ha")
                        print(f"  v_xc range: [{vxc_diag.min():.6f}, {vxc_diag.max():.6f}] Ha")

                    return vxc_diag

                except Exception as e:
                    raise RuntimeError(f"Failed to extract DFT v_xc: {e}")

            elif is_hf or hasattr(mf, 'get_k'):
                # HF case: v_xc = -0.5 * K (exact exchange)
                try:
                    # CRITICAL: Must NOT use DF approximation!
                    # Check if mf was created with DF
                    if hasattr(mf, 'with_df') and mf.with_df is not None:
                        # mf uses DF - we need to recompute K exactly
                        if self.config.verbose >= 1:
                            print("WARNING: Mean-field uses DF. Recomputing K matrix with exact integrals...")

                        # Build density matrix
                        dm = mf.make_rdm1()

                        # Create a new HF object without DF for exact K
                        mf_exact = scf.RHF(mf.mol)
                        mf_exact.with_df = None  # Force exact integrals
                        vj, vk = mf_exact.get_jk(mf.mol, dm)
                    else:
                        # mf already uses exact integrals
                        if self.config.verbose >= 1:
                            print("Computing v_xc from exact HF exchange (no DF)...")

                        # Get the exact exchange matrix
                        dm = mf.make_rdm1()
                        vj, vk = mf.get_jk(mf.mol, dm)

                    # Transform to MO basis
                    vxc_ao = -0.5 * vk  # Negative for closed-shell
                    vxc_mo = np.dot(mf.mo_coeff.T, np.dot(vxc_ao, mf.mo_coeff))
                    vxc_diag = np.diag(vxc_mo)

                    if self.config.verbose >= 1:
                        print(f"✓ Extracted exact HF v_xc (no DF approximation)")
                        print(f"  HOMO v_xc = {vxc_diag[self.n_occ-1]:.6f} Ha = {vxc_diag[self.n_occ-1]*27.211:.3f} eV")
                        print(f"  v_xc range: [{vxc_diag.min():.6f}, {vxc_diag.max():.6f}] Ha")

                    # Sanity check for HF
                    if np.any(vxc_diag[:self.n_occ] > 0):
                        warnings.warn("ERROR: HF v_xc has positive values for occupied orbitals!")

                    return vxc_diag

                except Exception as e:
                    raise RuntimeError(f"Failed to extract exact HF v_xc: {e}")

            else:
                # Unknown SCF type - try generic approach
                warnings.warn("Unknown SCF type - attempting generic v_xc extraction")

                # Try to get veff and subtract Coulomb
                if hasattr(mf, 'get_veff'):
                    dm = mf.make_rdm1()
                    veff = mf.get_veff(dm=dm)

                    # Get Coulomb part
                    if hasattr(mf, 'get_j'):
                        vj = mf.get_j(dm=dm)
                    else:
                        vj, _ = mf.get_jk(mf.mol, dm)

                    # v_xc = veff - vj
                    vxc_ao = veff - vj
                    vxc_mo = np.dot(mf.mo_coeff.T, np.dot(vxc_ao, mf.mo_coeff))
                    vxc_diag = np.diag(vxc_mo)

                    return vxc_diag

        # Last resort: Use mol object directly if available
        elif self.mol is not None:
            try:
                if self.config.verbose >= 1:
                    print("Computing v_xc from mol object with exact integrals...")

                # Build density matrix
                dm = 2.0 * np.dot(self.mo_coeff[:, :self.n_occ],
                                self.mo_coeff[:, :self.n_occ].T)

                # Compute exact exchange
                vj, vk = scf.hf.get_jk(self.mol, dm, hermi=1)

                vxc_ao = -0.5 * vk
                vxc_mo = np.dot(self.mo_coeff.T, np.dot(vxc_ao, self.mo_coeff))
                vxc_diag = np.diag(vxc_mo)

                return vxc_diag

            except Exception as e:
                raise RuntimeError(f"Failed to compute v_xc from mol object: {e}")

        # If we reach here, we cannot compute exact v_xc
        raise RuntimeError(
            "CRITICAL: Cannot compute exact v_xc!\n"
            "v_xc must be computed with exact integrals (not DF) to avoid cancellation.\n"
            "Please provide either:\n"
            "  1. A PySCF mf object (RHF/RKS) with mol attribute\n"
            "  2. A PySCF mol object\n"
            "Without exact v_xc, G0W0 results will have errors > 1 eV!"
        )

    def _build_df_tensors(self):
        """Build density-fitted 3-center integrals for Rust interface.

        Uses PySCF's Cholesky-decomposed integrals from df._cderi which are
        properly normalized such that (μν|ρσ) ≈ Σ_P L_μν,P L_ρσ,P.

        Returns:
            tuple: (iaP, ijP, v_aux) tensors
                - iaP: (n_occ*n_vir, n_aux) occupied-virtual DF tensors
                - ijP: (n_occ*n_occ, n_aux) occupied-occupied DF tensors
                - v_aux: (n_aux, n_aux) Identity matrix (no metric needed for Cholesky)

        Raises:
            RuntimeError: If PySCF or molecule object is not available
            RuntimeError: If DF construction fails
            ValueError: If DF tensors have invalid values
        """
        if not HAS_PYSCF:
            raise RuntimeError(
                "PySCF is required for DF tensor construction.\n"
                "Install with: pip install pyscf"
            )

        if self.mol is None:
            raise RuntimeError(
                "Molecule object is required for DF tensor construction.\n"
                "Please provide a PySCF molecule when initializing GWDriver."
            )

        try:
            from pyscf.lib import unpack_tril
        except ImportError:
            # Fallback implementation of unpack_tril
            def unpack_tril(x, axis=-1):
                """Unpack lower triangular matrix from compact form."""
                if axis == -1:
                    n = int((2 * x.shape[-1])**0.5)
                    idx = np.tril_indices(n)
                    mat = np.zeros((n, n))
                    mat[idx] = x
                    mat[idx[::-1]] = x
                    return mat
                else:
                    raise NotImplementedError("Only axis=-1 supported")

        # Step 1: Build auxiliary basis
        if self.config.verbose >= 1:
            print(f"  Step 1/4: Building auxiliary basis '{self.config.auxbasis}'...")

        try:
            # CRITICAL FIX: PySCF requires auxbasis as dict mapping atoms to basis sets
            # Convert string auxbasis to dict format if needed
            auxbasis = self.config.auxbasis
            if isinstance(auxbasis, str):
                auxbasis = {self.mol.atom_symbol(i): auxbasis for i in range(self.mol.natm)}
            auxmol = df.addons.make_auxmol(self.mol, auxbasis=auxbasis)
        except Exception as e:
            raise RuntimeError(
                f"Failed to create auxiliary basis '{self.config.auxbasis}'.\n"
                f"Error: {e}\n"
                f"The auxiliary basis may not be available for your atomic basis set.\n"
                f"Try: 'def2-tzvp-jkfit', 'weigend', or 'weigend+etb'"
            ) from e

        self.aux_mol = auxmol  # Store for later use
        naux = auxmol.nao_nr()
        nmo = len(self.mo_energy)
        nocc = self.n_occ
        nvir = nmo - nocc
        nao = self.mol.nao_nr()

        if self.config.verbose >= 1:
            print(f"    ✓ Auxiliary basis: {naux} functions")
            print(f"    ✓ AO basis: {nao} functions")
            print(f"    ✓ MO basis: {nmo} orbitals ({nocc} occ, {nvir} vir)")

        # Step 2: Build 3-center integrals using PySCF's DF module
        if self.config.verbose >= 1:
            print(f"  Step 2/4: Computing Cholesky-decomposed 3-center integrals...")

        try:
            df_obj = df.DF(self.mol, auxbasis=self.config.auxbasis)
            df_obj.kernel()
        except Exception as e:
            raise RuntimeError(
                f"Failed to compute DF integrals with PySCF.\n"
                f"Error: {e}\n"
                f"This may indicate:\n"
                f"  - Linear dependency in auxiliary basis\n"
                f"  - Numerical instability in Cholesky decomposition\n"
                f"  - Insufficient memory"
            ) from e

        # Step 3: Get Cholesky-decomposed integrals
        if self.config.verbose >= 1:
            print(f"  Step 3/4: Extracting Cholesky vectors...")

        if not hasattr(df_obj, '_cderi'):
            raise RuntimeError(
                "PySCF DF object has no '_cderi' attribute.\n"
                "This indicates the DF.kernel() did not run successfully."
            )

        if df_obj._cderi is None:
            raise RuntimeError(
                "PySCF DF object has '_cderi' = None.\n"
                "This indicates the Cholesky decomposition failed silently.\n"
                "Try increasing basis set quality or checking molecule geometry."
            )

        cderi = df_obj._cderi

        # Validate Cholesky array structure
        if cderi.ndim != 2:
            raise ValueError(
                f"Invalid Cholesky array: expected 2D, got shape {cderi.shape}.\n"
                f"This indicates a fundamental issue with PySCF's DF module."
            )

        expected_compact_size = nao * (nao + 1) // 2
        if cderi.shape[1] != expected_compact_size:
            raise ValueError(
                f"Cholesky array has wrong compact dimension.\n"
                f"Expected: {expected_compact_size} (from {nao} AO basis functions)\n"
                f"Got: {cderi.shape[1]}\n"
                f"This indicates a mismatch between molecule and DF computation."
            )

        if cderi.shape[0] != naux:
            raise ValueError(
                f"Cholesky array has wrong auxiliary dimension.\n"
                f"Expected: {naux} (from auxiliary basis)\n"
                f"Got: {cderi.shape[0]}"
            )

        # Check for pathological values
        cderi_norm = np.linalg.norm(cderi)
        if cderi_norm < 1e-10:
            raise ValueError(
                f"Cholesky integrals have suspiciously small norm: {cderi_norm:.3e}\n"
                f"This indicates numerical issues in DF construction."
            )

        if np.any(np.isnan(cderi)) or np.any(np.isinf(cderi)):
            raise ValueError(
                "Cholesky integrals contain NaN or Inf values!\n"
                "This indicates a serious numerical problem."
            )

        if self.config.verbose >= 1:
            print(f"    ✓ Cholesky array shape: {cderi.shape}")
            print(f"    ✓ Frobenius norm: {cderi_norm:.3e}")

        # Step 4: Transform Cholesky vectors to MO basis
        if self.config.verbose >= 1:
            print(f"  Step 4/4: Transforming to MO basis...")

        iaP_3d = np.zeros((nocc, nvir, naux))
        # CRITICAL: For exchange self-energy, we need (m,i,P) for ALL m and occupied i
        miP_3d = np.zeros((nmo, nocc, naux))

        try:
            for p in range(naux):
                # Unpack the Cholesky vector for auxiliary function p
                Lp = unpack_tril(cderi[p])

                # Validate unpacked matrix
                if np.any(np.isnan(Lp)) or np.any(np.isinf(Lp)):
                    raise ValueError(f"Unpacked Cholesky vector {p} contains NaN/Inf")

                # Transform to MO basis: L_pq,P = C_μp L_μν,P C_νq
                Lp_half = np.dot(Lp, self.mo_coeff)
                Lp_mo = np.dot(self.mo_coeff.T, Lp_half)

                # Validate transformed matrix
                if np.any(np.isnan(Lp_mo)) or np.any(np.isinf(Lp_mo)):
                    raise ValueError(f"MO-transformed Cholesky vector {p} contains NaN/Inf")

                # Extract blocks
                iaP_3d[:, :, p] = Lp_mo[:nocc, nocc:]
                miP_3d[:, :, p] = Lp_mo[:, :nocc]

                # Progress indicator for large systems
                if self.config.verbose >= 2 and (p + 1) % 50 == 0:
                    print(f"    Processed {p+1}/{naux} auxiliary functions")

        except Exception as e:
            raise RuntimeError(
                f"Failed during MO transformation of Cholesky vectors.\n"
                f"Error at auxiliary function {p}/{naux}: {e}\n"
                f"This may indicate:\n"
                f"  - Numerical instability in MO coefficients\n"
                f"  - Memory issues with large tensors"
            ) from e

        # Reshape to 2D for Rust interface
        iaP = iaP_3d.reshape(nocc * nvir, naux)
        miP = miP_3d.reshape(nmo * nocc, naux)

        # Validate final tensors
        iaP_norm = np.linalg.norm(iaP)
        miP_norm = np.linalg.norm(miP)

        if iaP_norm < 1e-10:
            raise ValueError(
                f"iaP tensor has suspiciously small norm: {iaP_norm:.3e}\n"
                f"This indicates the DF tensors are essentially zero."
            )

        if miP_norm < 1e-10:
            raise ValueError(
                f"miP tensor has suspiciously small norm: {miP_norm:.3e}\n"
                f"This indicates the DF tensors are essentially zero."
            )

        if np.any(np.isnan(iaP)) or np.any(np.isinf(iaP)):
            raise ValueError("iaP tensor contains NaN or Inf values!")

        if np.any(np.isnan(miP)) or np.any(np.isinf(miP)):
            raise ValueError("miP tensor contains NaN or Inf values!")

        # CRITICAL FIX: Use j2c (FORWARD Coulomb metric) for screening, NOT j2c_inv!
        #
        # PySCF convention (from pyscf/gw/gw_cd.py, get_W):
        #   v_sqrt = cholesky(j2c)  ← j2c is the FORWARD metric (P|Q)!
        #   epsilon = I - v_sqrt @ P0 @ v_sqrt.T
        #   W = v_sqrt @ epsilon_inv @ v_sqrt.T
        #
        # The formula W = v^{1/2} [I - v^{1/2} P0 v^{1/2}]^{-1} v^{1/2}
        # uses v = j2c (forward), NOT j2c_inv (inverse)!
        #
        # Reference: pyscf/df/incore.py::fill_2c2e(), pyscf/gw/gw_cd.py lines ~180
        if self.config.verbose >= 1:
            print("  Computing Coulomb metric j2c...")

        # Build 2-center Coulomb integrals (P|Q) - this is the forward metric!
        v_aux = df.incore.fill_2c2e(self.mol, auxmol)

        # Verify metric is positive definite
        j2c_eigvals = np.linalg.eigvalsh(v_aux)
        if np.any(j2c_eigvals < -1e-10):
            raise RuntimeError(f"j2c has negative eigenvalues: min = {np.min(j2c_eigvals)}")

        if self.config.verbose >= 1:
            print(f"    ✓ j2c computed: shape {v_aux.shape}")
            print(f"    ✓ Eigenvalue range: [{np.min(j2c_eigvals):.3e}, {np.max(j2c_eigvals):.3e}]")

        # Success: Log completion
        if self.config.verbose >= 1:
            print(f"\n✓ DF tensor construction successful!")
            print(f"  iaP: {iaP.shape}, norm = {iaP_norm:.3e}")
            print(f"  miP: {miP.shape}, norm = {miP_norm:.3e}")
            print(f"  v_aux: {v_aux.shape} (j2c forward Coulomb metric)")
            print(f"  Memory usage: {(iaP.nbytes + miP.nbytes + v_aux.nbytes) / 1024**2:.1f} MB")

        # Store both 2D (for Rust) and 3D (for Python) versions
        self.df_tensors = {
            'iaP': iaP,  # 2D for Rust interface (occupied-virtual)
            'ijP': miP,  # 2D for Rust interface (ALL-occupied, NOT occ-occ!)
            'iaP_3d': iaP_3d,  # 3D for Python use
            'ijP_3d': miP_3d,  # 3D for Python use (renamed from ijP_3d)
            'v_aux': v_aux,  # j2c (forward Coulomb metric) for screening
            'naux': naux,
            'nocc': nocc,
            'nvir': nvir,
            'using_cholesky': True  # Flag to indicate Cholesky vectors
        }

        # Build full DF tensor for correlation self-energy
        if 'mnP' not in self.df_tensors:
            try:
                self.df_tensors['mnP'] = self._build_full_df_tensor()

                if self.config.verbose >= 1:
                    print("✓ Full DF tensor available for correlation")
            except Exception as e:
                import warnings
                warnings.warn(f"Failed to build full DF tensor: {e}")
                warnings.warn("Correlation self-energy will use approximate method")
                # Don't fail - just warn

        return iaP, miP, v_aux  # Return j2c (forward Coulomb metric)

    def _build_full_df_tensor(self):
        """Build full (n_mo, n_mo, n_aux) DF tensor for correlation.

        This is required for the correct GW correlation self-energy formula,
        which needs diagonal elements ⟨mm|W|mm⟩ that cannot be computed
        from the block (ia|P) structure.

        Returns:
            np.ndarray: Full DF tensor with shape (n_mo, n_mo, n_aux)

        Raises:
            RuntimeError: If DF construction fails
            ValueError: If tensor has invalid values
        """
        from pyscf import df
        import numpy as np

        if self.config.verbose >= 1:
            print("\n=== Building full DF tensor for correlation self-energy ===")
            print("  (This is required for accurate Σ_c calculation)")

        nmo = len(self.mo_energy)

        # Get auxiliary molecule if not already stored
        if not hasattr(self, 'aux_mol'):
            # CRITICAL FIX: PySCF requires auxbasis as dict
            auxbasis = self.config.auxbasis
            if isinstance(auxbasis, str):
                auxbasis = {self.mol.atom_symbol(i): auxbasis for i in range(self.mol.natm)}
            self.aux_mol = df.addons.make_auxmol(self.mol, auxbasis=auxbasis)

        naux = self.aux_mol.nao_nr()
        nao = self.mol.nao_nr()

        # Build 3-center integrals using PySCF's DF module
        df_obj = df.DF(self.mol, auxbasis=self.config.auxbasis)
        df_obj.kernel()

        # Get Cholesky-decomposed integrals
        if not hasattr(df_obj, '_cderi') or df_obj._cderi is None:
            raise RuntimeError("Failed to obtain Cholesky integrals from PySCF DF object")

        cderi = df_obj._cderi

        # Helper function for unpacking triangular matrices
        try:
            from pyscf.lib import unpack_tril
        except ImportError:
            def unpack_tril(x, axis=-1):
                """Unpack lower triangular matrix from compact form."""
                if axis == -1:
                    n = int((2 * x.shape[-1])**0.5)
                    idx = np.tril_indices(n)
                    mat = np.zeros((n, n))
                    mat[idx] = x
                    mat[idx[::-1]] = x
                    return mat
                else:
                    raise NotImplementedError("Only axis=-1 supported")

        # Build full tensor: (n_mo, n_mo, n_aux)
        # NOTE: We're using Cholesky-decomposed integrals which are already properly
        # normalized, so no metric correction is needed
        mnP_3d = np.zeros((nmo, nmo, naux))

        if self.config.verbose >= 1:
            print(f"  Transforming {naux} auxiliary basis functions to MO basis...")

        for P in range(naux):
            # Unpack the Cholesky vector for auxiliary function P
            # This gives us L_μν,P in full matrix form (symmetric)
            Lp_ao = unpack_tril(cderi[P])

            # Transform to MO basis: L_mn,P = C_μm^T L_μν,P C_νn
            # This is the standard two-index transformation
            # First transform columns: L_μn,P = L_μν,P @ C
            Lp_half = Lp_ao @ self.mo_coeff  # Shape: (nao, nmo)
            # Then transform rows: L_mn,P = C^T @ L_μn,P
            Lp_mo = self.mo_coeff.T @ Lp_half  # Shape: (nmo, nmo)

            # Store in 3D tensor
            mnP_3d[:, :, P] = Lp_mo

            if self.config.verbose >= 2 and P % 20 == 0:
                print(f"    Processed {P+1}/{naux} auxiliary functions")

        # Validate final tensor
        norm = np.linalg.norm(mnP_3d)

        if norm < 1e-10:
            raise ValueError(
                f"Full DF tensor has suspiciously small norm: {norm:.3e}\n"
                f"This indicates numerical issues in DF construction."
            )

        if np.any(np.isnan(mnP_3d)):
            raise ValueError(
                "Full DF tensor contains NaN values!\n"
                "This indicates a serious numerical problem."
            )

        if np.any(np.isinf(mnP_3d)):
            raise ValueError(
                "Full DF tensor contains Inf values!\n"
                "This indicates a serious numerical problem."
            )

        # Check that diagonal elements are non-zero (critical for correlation)
        diag_norms = [np.linalg.norm(mnP_3d[i, i, :]) for i in range(min(5, nmo))]
        if all(d < 1e-10 for d in diag_norms):
            raise ValueError(
                "All diagonal DF elements are zero!\n"
                "This will cause incorrect correlation self-energy."
            )

        if self.config.verbose >= 1:
            print(f"\n✓ Full DF tensor construction successful!")
            print(f"  Shape: {mnP_3d.shape}")
            print(f"  Frobenius norm: {norm:.3e}")
            print(f"  Memory usage: {mnP_3d.nbytes / 1024**2:.1f} MB")
            print(f"  Sample diagonal norms: {[f'{d:.3e}' for d in diag_norms[:3]]}")

        return mnP_3d

    def _dict_to_gwresult(self, result_dict: Dict[str, Any]) -> GWResult:
        """Convert result dictionary to GWResult object.

        Args:
            result_dict: Dictionary with calculation results

        Returns:
            GWResult object
        """
        # Ensure all required fields are present
        homo_idx = result_dict.get('homo', self.n_occ - 1)
        lumo_idx = result_dict.get('lumo', self.n_occ)

        # Calculate IP, EA, and gap if not provided
        qp_energies = result_dict['qp_energies']

        # DEBUG: Log raw values from Rust
        if 'ip' in result_dict:
            print(f"DEBUG _dict_to_gwresult: Raw IP from Rust = {result_dict['ip']:.6f}")
        if 'ea' in result_dict:
            print(f"DEBUG _dict_to_gwresult: Raw EA from Rust = {result_dict.get('ea')}")
        if 'gap' in result_dict:
            print(f"DEBUG _dict_to_gwresult: Raw gap from Rust = {result_dict['gap']:.6f}")
        print(f"DEBUG _dict_to_gwresult: HOMO QP energy = {qp_energies[homo_idx]:.6f} Ha")
        print(f"DEBUG _dict_to_gwresult: LUMO QP energy = {qp_energies[lumo_idx]:.6f} Ha")

        if 'ip' not in result_dict:
            result_dict['ip'] = -qp_energies[homo_idx] * HA_TO_EV if homo_idx >= 0 else 0.0
        else:
            # Ensure IP is in eV
            print(f"DEBUG _dict_to_gwresult: Checking IP conversion: abs(IP)={abs(result_dict['ip']):.6f}, threshold=1.0")
            if abs(result_dict['ip']) < 1.0:  # Likely in Hartree
                print(f"DEBUG _dict_to_gwresult: Converting IP from Ha to eV: {result_dict['ip']:.6f} -> {result_dict['ip'] * HA_TO_EV:.6f}")
                result_dict['ip'] *= HA_TO_EV
            else:
                print(f"DEBUG _dict_to_gwresult: IP already in eV, not converting")

        if 'ea' not in result_dict:
            result_dict['ea'] = -qp_energies[lumo_idx] * HA_TO_EV if lumo_idx < len(qp_energies) else None
        elif result_dict['ea'] is not None and abs(result_dict['ea']) < 1.0:  # Likely in Hartree
            result_dict['ea'] *= HA_TO_EV

        if 'gap' not in result_dict:
            if lumo_idx < len(qp_energies) and homo_idx >= 0:
                result_dict['gap'] = (qp_energies[lumo_idx] - qp_energies[homo_idx]) * HA_TO_EV
            else:
                result_dict['gap'] = 0.0
        elif abs(result_dict['gap']) < 1.0:  # Likely in Hartree
            result_dict['gap'] *= HA_TO_EV

        # DEBUG: Log final values after conversion
        print(f"DEBUG _dict_to_gwresult: Final IP = {result_dict['ip']:.6f} eV")
        print(f"DEBUG _dict_to_gwresult: Final EA = {result_dict.get('ea')}")
        print(f"DEBUG _dict_to_gwresult: Final gap = {result_dict['gap']:.6f} eV")

        # Handle sigma_c - Rust returns sigma_c_re and sigma_c_im separately
        if 'sigma_c' in result_dict:
            sigma_c = result_dict['sigma_c']
        elif 'sigma_c_re' in result_dict and 'sigma_c_im' in result_dict:
            # Combine real and imaginary parts
            sigma_c = result_dict['sigma_c_re'] + 1j * result_dict['sigma_c_im']

            # Check if sigma_c_re is all zeros (bug in Rust?)
            if np.all(np.abs(result_dict['sigma_c_re']) < 1e-10):
                import warnings
                warnings.warn(
                    "WARNING: Correlation self-energy (sigma_c_re) is all zeros!\n"
                    "This indicates a bug in the Rust implementation.\n"
                    "The GW calculation will be incorrect."
                )
        else:
            sigma_c = np.zeros_like(qp_energies)

        # Create GWResult object
        return GWResult(
            converged=result_dict.get('converged', True),
            n_iterations=result_dict.get('n_iterations', 1),
            qp_energies=qp_energies,
            z_factors=result_dict.get('z_factors', np.ones_like(qp_energies)),
            sigma_x=result_dict.get('sigma_x', np.zeros_like(qp_energies)),
            sigma_c=sigma_c,
            homo_index=homo_idx,
            lumo_index=lumo_idx,
            ip=result_dict['ip'],
            ea=result_dict['ea'],
            gap=result_dict['gap'],
            elapsed_seconds=result_dict.get('elapsed_seconds', 0.0),
            convergence_data=result_dict.get('convergence_data'),
            spectral_data=result_dict.get('spectral_data')
        )

    def _verify_df_tensors(self):
        """Verify DF tensors are correctly constructed."""
        if self.config.verbose >= 2:
            print("\nVerifying DF tensors...")

            # Check shapes
            if 'mnP' in self.df_tensors:
                mnP = self.df_tensors['mnP']
                nmo = len(self.mo_energy)
                naux = self.aux_mol.nao_nr()

                expected_shape = (nmo, nmo, naux)
                if mnP.shape != expected_shape:
                    raise ValueError(f"Full DF tensor shape {mnP.shape} != expected {expected_shape}")

                # Check diagonal elements exist
                diag_elements = np.array([mnP[i, i, :] for i in range(min(5, nmo))])
                if np.all(np.abs(diag_elements) < 1e-10):
                    import warnings
                    warnings.warn("Diagonal DF elements are zero - this is suspicious!")

                print(f"  ✓ Full tensor shape: {mnP.shape}")
                print(f"  ✓ Diagonal elements present: {np.max(np.abs(diag_elements)):.6e}")

            # Check block tensors
            if 'iaP' in self.df_tensors:
                iaP = self.df_tensors['iaP']
                print(f"  ✓ Block iaP shape: {iaP.shape}")

            if 'ijP' in self.df_tensors:
                ijP = self.df_tensors['ijP']
                print(f"  ✓ Block miP shape: {ijP.shape}")

    def kernel(self, method='evgw', qp_solver=None, **kwargs) -> GWResult:
        """Run GW calculation.

        Args:
            method: GW method to use ('g0w0', 'evgw'). Default: 'evgw'
            qp_solver: QP equation solver ('linearized' or 'newton').
                       If None, uses config default. 'newton' solves the full
                       QP equation via Newton-Raphson for higher accuracy.
            **kwargs: Additional keyword arguments

        Returns:
            GWResult object with attributes:
            - qp_energies: Quasiparticle energies (Hartree)
            - z_factors: Renormalization factors
            - converged: Convergence status
            - n_iterations: Number of iterations
            - ip: Ionization potential (eV)
            - ea: Electron affinity (eV)
            - gap: HOMO-LUMO gap (eV)

        Raises:
            RuntimeError: If PySCF or Rust implementation is not available
            RuntimeError: If DF tensor construction fails
        """
        # Override config qp_solver if explicitly provided
        if qp_solver is not None:
            self.config.qp_solver = qp_solver

        start_time = time.time()

        # CRITICAL: NO MOCK IMPLEMENTATIONS - QuasiX Golden Rule
        if not HAS_PYSCF:
            raise RuntimeError(
                "PySCF is not available. QuasiX requires real scientific computing libraries.\n"
                "Please install PySCF: pip install pyscf"
            )

        if self.mol is None:
            raise RuntimeError(
                "Molecule object is required for GW calculations.\n"
                "Please provide a PySCF molecule object when initializing GWDriver."
            )

        # Priority 1: Try Rust evGW (fastest, production-quality)
        if RUST_AVAILABLE and rust_module is not None:
            try:
                if self.config.verbose >= 1:
                    print("Using Rust evGW implementation (production path)...")
                result = self._run_full_gw(**kwargs)
                result['elapsed_seconds'] = time.time() - start_time
                result['implementation'] = 'Rust'
                return self._dict_to_gwresult(result)
            except Exception as e:
                if self.config.verbose >= 0:
                    print(f"\nWARNING: Rust evGW failed: {e}")
                    print("Falling back to Python EvGW implementation...")
                # Fall through to Python fallback
        else:
            if self.config.verbose >= 0:
                print("WARNING: Rust extension not available")
                print("  Build with: cd quasix && maturin develop --release")
                print("Attempting Python EvGW fallback...")

        # Priority 2: Fall back to Python EvGW (slower, but still functional)
        if QUASIX_AVAILABLE and EvGW is not None:
            try:
                if self.gw_obj is None:
                    raise RuntimeError("Python EvGW object not initialized")

                if self.config.verbose >= 1:
                    print("Using Python EvGW implementation (fallback path)...")

                # Run Python EvGW
                gw_result = self.gw_obj.kernel(**kwargs)

                # Convert to standard result format
                result = {
                    'qp_energies': self.gw_obj.qp_energies,
                    'z_factors': self.gw_obj.z_factors,
                    'converged': self.gw_obj.converged,
                    'n_iterations': len(self.gw_obj.iteration_history),
                    'max_error': 0.0,
                    'homo': self.n_occ - 1,
                    'lumo': self.n_occ,
                    'gap': 0.0,
                    'ip': -self.gw_obj.qp_energies[self.n_occ - 1] if self.n_occ > 0 else 0.0,
                    'ea': -self.gw_obj.qp_energies[self.n_occ] if self.n_occ < len(self.gw_obj.qp_energies) else None,
                    'sigma_x': self.gw_obj.sigma_x if self.gw_obj.sigma_x is not None else np.zeros(len(self.mo_energy)),
                    'sigma_c': self.gw_obj.sigma_c if self.gw_obj.sigma_c is not None else np.zeros(len(self.mo_energy)),
                    'elapsed_seconds': time.time() - start_time,
                    'implementation': 'Python'
                }

                # Calculate gap
                if self.n_occ > 0 and self.n_occ < len(result['qp_energies']):
                    result['gap'] = result['qp_energies'][self.n_occ] - result['qp_energies'][self.n_occ - 1]

                return self._dict_to_gwresult(result)

            except Exception as e:
                if self.config.verbose >= 0:
                    print(f"\nERROR: Python EvGW also failed: {e}")
                raise RuntimeError(
                    f"Both Rust and Python evGW implementations failed.\n"
                    f"Rust error: {e}\n"
                    f"This is a critical failure - no working implementation available."
                ) from e
        else:
            raise RuntimeError(
                "No GW implementation available.\n"
                "- Rust extension not built: cd quasix && maturin develop --release\n"
                "- Python EvGW not available: Check QuasiX installation"
            )

    def _run_full_gw(self, **kwargs) -> Dict[str, Any]:
        """Run full evGW calculation using Rust implementation.

        Raises:
            RuntimeError: If DF tensor construction fails
            RuntimeError: If Rust calculation fails
        """

        # Build DF tensors if not already done
        if self.df_tensors is None:
            if self.config.verbose >= 1:
                print("Building DF tensors from PySCF...")

            try:
                self._build_df_tensors()
            except Exception as e:
                import traceback
                error_msg = (
                    f"CRITICAL ERROR: Failed to build DF tensors from PySCF.\n"
                    f"Error: {e}\n"
                    f"Traceback:\n{traceback.format_exc()}\n\n"
                    f"This is a fatal error. QuasiX requires real DF tensors for accurate GW calculations.\n"
                    f"Possible causes:\n"
                    f"  1. PySCF DF module failed to compute Cholesky integrals\n"
                    f"  2. Auxiliary basis set '{self.config.auxbasis}' is incompatible with molecule\n"
                    f"  3. Basis set or molecule structure issues\n\n"
                    f"Suggested fixes:\n"
                    f"  - Try a different auxiliary basis (e.g., 'def2-tzvp-jkfit' or 'weigend')\n"
                    f"  - Verify molecule geometry and basis set are valid\n"
                    f"  - Check PySCF installation: python -c 'from pyscf import df; print(df)'\n"
                )
                raise RuntimeError(error_msg) from e

        # Verify DF tensors if in debug mode
        if self.config.verbose >= 2:
            self._verify_df_tensors()

        # Get VXC from exact HF exchange (NOT DF approximation!)
        # CRITICAL: For RHF, vxc = -0.5*vk (exchange only)
        # Note: PySCF's get_k() returns POSITIVE exchange matrix K where K[i,i] = 2*Σⱼ (ij|ij) > 0
        # The exchange potential for closed-shell RHF is: v_xc = -K/2 = -0.5*vk < 0 (negative)
        # IMPORTANT: We must use EXACT 4-center integrals, NOT DF approximation!
        # If vxc and Σx both come from DF, they cancel identically in the QP equation!

        # Use the new exact extraction method
        try:
            vxc_dft = self.extract_vxc_exact()

            # Validate v_xc values
            if np.any(np.isnan(vxc_dft)) or np.any(np.isinf(vxc_dft)):
                raise ValueError("v_xc contains NaN or Inf values")

            if self.config.verbose >= 1:
                print(f"\n✓ Successfully extracted exact v_xc")
                print(f"  HOMO v_xc = {vxc_dft[self.n_occ-1]:.6f} Ha = {vxc_dft[self.n_occ-1]*27.211:.3f} eV")
                print(f"  v_xc range: [{vxc_dft.min():.6f}, {vxc_dft.max():.6f}] Ha")

                # For HF, check that occupied orbitals have negative v_xc
                if hasattr(self, 'mf') and self.mf is not None:
                    from pyscf import scf
                    if isinstance(self.mf, (scf.hf.RHF, scf.UHF, scf.ROHF)) if hasattr(scf, 'hf') else False:
                        if np.any(vxc_dft[:self.n_occ] > 0):
                            warnings.warn("WARNING: HF v_xc has positive values for occupied orbitals!")

        except (RuntimeError, ValueError) as e:
            # CRITICAL: Without exact v_xc, G0W0 will be wrong!
            error_msg = (
                "\n" + "=" * 70 + "\n"
                "CRITICAL ERROR: Cannot extract exact v_xc!\n"
                f"Error: {e}\n"
                "\n"
                "v_xc MUST be computed with exact integrals to avoid cancellation.\n"
                "Please ensure you provide either:\n"
                "  1. A PySCF mf object (RHF/RKS) with mol attribute\n"
                "  2. A PySCF mol object with mo_coeff data\n"
                "\n"
                "Without exact v_xc, G0W0 results will have errors > 1 eV!\n"
                "=" * 70
            )
            raise RuntimeError(error_msg)

        # Debug: Print vxc values AND mo_energy being passed to Rust
        if self.config.verbose >= 1:
            print(f"\nDEBUG: Passing data to Rust:")
            print(f"  mo_energy shape: {self.mo_energy.shape}")
            print(f"  mo_energy HOMO (idx={self.n_occ-1}): {self.mo_energy[self.n_occ-1]:.6f} Ha = {self.mo_energy[self.n_occ-1]*27.211:.2f} eV")
            print(f"  mo_energy LUMO (idx={self.n_occ}): {self.mo_energy[self.n_occ]:.6f} Ha = {self.mo_energy[self.n_occ]*27.211:.2f} eV")
            print(f"  mo_energy range: [{self.mo_energy.min():.6f}, {self.mo_energy.max():.6f}] Ha")
            print(f"\nPassing vxc to Rust:")
            print(f"  Shape: {vxc_dft.shape}")
            print(f"  HOMO (orbital {self.n_occ-1}): vxc = {vxc_dft[self.n_occ-1]:.4f} Ha = {vxc_dft[self.n_occ-1]*27.211:.2f} eV")
            print(f"  Min: {vxc_dft.min():.4f} Ha, Max: {vxc_dft.max():.4f} Ha")

        # Call Rust evGW implementation using evgw function
        # CRITICAL: Use evgw() with positional arguments (not run_evgw with dict)
        # The evgw function is properly declared with #[pyfunction] and works correctly
        try:
            if not hasattr(rust_module, 'evgw'):
                raise RuntimeError(
                    "evgw function not found in Rust module.\n"
                    "Please rebuild QuasiX: cd quasix && maturin develop --release"
                )

            # Prepare tensors - evgw expects 2D arrays, not 3D!
            # iaP: (nocc*nvir, naux)
            # ijP (miP): (nmo*nocc, naux)
            # abP: (nvir*nvir, naux) - may not be available
            iaP_2d = np.ascontiguousarray(self.df_tensors['iaP'], dtype=np.float64)
            miP_2d = np.ascontiguousarray(self.df_tensors['ijP'], dtype=np.float64)
            chol_v = np.ascontiguousarray(self.df_tensors['v_aux'], dtype=np.float64)
            mo_energy_contig = np.ascontiguousarray(self.mo_energy, dtype=np.float64)
            mo_occ_contig = np.ascontiguousarray(self.mo_occ, dtype=np.float64)
            vxc_dft_contig = np.ascontiguousarray(vxc_dft, dtype=np.float64)

            # Build abP if not available (backward compatibility)
            nocc = self.n_occ
            nvir = len(self.mo_energy) - nocc
            naux = self.df_tensors['naux']

            # Create empty abP (will be filled by Rust if needed)
            abP_2d = np.zeros((nvir * nvir, naux), dtype=np.float64, order='C')

            if self.config.verbose >= 1:
                print(f"\nCalling Rust evgw function...")
                print(f"  mo_energy shape: {mo_energy_contig.shape}")
                print(f"  mo_occ shape: {mo_occ_contig.shape}")
                print(f"  iaP shape: {iaP_2d.shape} (2D)")
                print(f"  ijP shape: {miP_2d.shape} (2D)")
                print(f"  abP shape: {abP_2d.shape} (2D)")
                print(f"  chol_v shape: {chol_v.shape}")
                print(f"  vxc_dft shape: {vxc_dft_contig.shape}")

            # Call evgw with all positional arguments
            # Use the signature from quasix/src/evgw.rs line 43
            result_dict = rust_module.evgw(
                mo_energy_contig,
                mo_occ_contig,
                iaP_2d,
                miP_2d,
                abP_2d,
                chol_v,
                vxc_dft_contig,
                max_cycle=self.config.max_cycle,
                conv_tol=self.config.conv_tol,
                conv_tol_z=self.config.conv_tol_z,
                damping=self.config.damping,
                damping_dynamic=self.config.damping_dynamic,
                diis=self.config.diis,
                diis_space=self.config.diis_space,
                nfreq=self.config.nfreq,
                eta=self.config.eta,
                verbose=self.config.verbose,
                # QP solver parameters
                qp_solver=self.config.qp_solver,
                newton_max_iterations=self.config.newton_max_iterations,
                newton_energy_tolerance=self.config.newton_energy_tolerance
            )

            # Process Rust result
            result = {
                'qp_energies': np.array(result_dict['qp_energies']),
                'z_factors': np.array(result_dict['z_factors']),
                'converged': result_dict['converged'],
                'n_iterations': result_dict['n_cycles'],
                'max_error': result_dict['final_error'],
                'homo': self.n_occ - 1,
                'lumo': self.n_occ,
                'gap': 0.0,  # Will calculate below
                'ip': 0.0,   # Will calculate below
                'ea': None,  # Will calculate below
                'sigma_x': np.array(result_dict['sigma_x']),
                'sigma_c': np.array(result_dict['sigma_c_re']) + 1j * np.array(result_dict['sigma_c_im']),
                'iteration_history': result_dict.get('iteration_history', [])
            }

            # Calculate derived properties
            qp_energies = result['qp_energies']
            if self.n_occ > 0:
                result['ip'] = -qp_energies[self.n_occ - 1]  # Negative of HOMO
            if self.n_occ < len(qp_energies):
                result['ea'] = -qp_energies[self.n_occ]  # Negative of LUMO
                if self.n_occ > 0:
                    result['gap'] = qp_energies[self.n_occ] - qp_energies[self.n_occ - 1]

            # Store for later analysis
            self._last_result = result

            # Validate Z factors if requested
            if self.config.validate_z and 'z_factors' in result:
                self._validate_z_factors(result['z_factors'])

            return result

        except Exception as e:
            import traceback
            error_msg = (
                f"CRITICAL ERROR: Rust evGW calculation failed.\n"
                f"Error: {e}\n"
                f"Traceback:\n{traceback.format_exc()}\n\n"
                f"This is a fatal error. QuasiX requires successful Rust kernel execution.\n"
                f"Possible causes:\n"
                f"  1. Invalid input data (check DF tensors, MO energies, occupations)\n"
                f"  2. Numerical instability in GW iterations\n"
                f"  3. Bug in Rust implementation\n\n"
                f"Debug information:\n"
                f"  - Number of orbitals: {len(self.mo_energy)}\n"
                f"  - Occupied orbitals: {self.n_occ}\n"
                f"  - DF tensor shapes: iaP={self.df_tensors['iaP'].shape}, ijP={self.df_tensors['ijP'].shape}\n"
                f"  - Auxiliary basis size: {self.df_tensors['naux']}\n"
            )
            raise RuntimeError(error_msg) from e

        # Fall back to Python EvGW if available
        if self.gw_obj is not None:
            try:
                # Run kernel - returns qp_energies (ndarray)
                qp_energies = self.gw_obj.kernel(**kwargs)
                # Get full result using analyze()
                evgw_result = self.gw_obj.analyze()
            except Exception as e:
                import traceback
                error_msg = (
                    f"CRITICAL ERROR: Python EvGW calculation failed.\n"
                    f"Error: {e}\n"
                    f"Traceback:\n{traceback.format_exc()}\n\n"
                    f"This is a fatal error. No fallback implementation available.\n"
                )
                raise RuntimeError(error_msg) from e

            # Handle both old and new result formats
            if isinstance(evgw_result, dict):
                # Old format - direct dictionary
                result = evgw_result
            elif hasattr(evgw_result, 'qp_energies'):
                # New format - EvGWResult dataclass
                result = {
                    'qp_energies': evgw_result.qp_energies,
                    'z_factors': evgw_result.z_factors,
                    'converged': evgw_result.converged,
                    'n_iterations': evgw_result.convergence_data.n_cycles if hasattr(evgw_result, 'convergence_data') else 0,
                    'max_error': evgw_result.convergence_data.final_error if hasattr(evgw_result, 'convergence_data') else 0.0,
                    'homo': evgw_result.homo_index if hasattr(evgw_result, 'homo_index') else self.n_occ - 1,
                    'lumo': evgw_result.lumo_index if hasattr(evgw_result, 'lumo_index') else self.n_occ,
                    'gap': evgw_result.gap_value if hasattr(evgw_result, 'gap_value') else 0.0,
                    'ip': evgw_result.ip if hasattr(evgw_result, 'ip') else 0.0,
                    'ea': evgw_result.ea if hasattr(evgw_result, 'ea') else None,
                    'sigma_x': evgw_result.sigma_x if hasattr(evgw_result, 'sigma_x') else np.zeros(len(self.mo_energy)),
                    'sigma_c': evgw_result.sigma_c if hasattr(evgw_result, 'sigma_c') else np.zeros(len(self.mo_energy)),
                }
            else:
                # Fallback for unexpected format
                error_msg = (
                    f"CRITICAL ERROR: Unexpected evGW result format.\n"
                    f"Result type: {type(evgw_result)}\n"
                    f"Result attributes: {dir(evgw_result) if hasattr(evgw_result, '__dict__') else 'N/A'}\n"
                    f"This indicates a mismatch between Python interface and Rust implementation.\n"
                )
                raise RuntimeError(error_msg)
        else:
            # No GW implementation available
            raise RuntimeError(
                "No GW implementation available. This should not happen - "
                "kernel() should have caught this earlier."
            )

        # Store for later analysis
        self._last_result = result

        # Validate Z factors if requested
        if self.config.validate_z and 'z_factors' in result:
            self._validate_z_factors(result['z_factors'])

        return result

    def _validate_z_factors(self, z_factors: np.ndarray):
        """Validate Z factors are physical."""
        if np.any(z_factors <= 0) or np.any(z_factors > 1):
            bad_indices = np.where((z_factors <= 0) | (z_factors > 1))[0]
            warnings.warn(
                f"Unphysical Z factors at orbitals {bad_indices}: "
                f"{z_factors[bad_indices]}"
            )

    def analyze(self) -> GWResult:
        """Analyze GW results and return structured output."""
        if not hasattr(self, '_last_result'):
            raise RuntimeError("No calculation has been performed yet")

        r = self._last_result

        return GWResult(
            converged=r['converged'],
            n_iterations=r['n_iterations'],
            qp_energies=r['qp_energies'],
            z_factors=r['z_factors'],
            sigma_x=r['sigma_x'],
            sigma_c=r['sigma_c'],
            homo_index=r['homo'],
            lumo_index=r['lumo'],
            ip=r['ip'] * HA_TO_EV,
            ea=r['ea'] * HA_TO_EV if r['ea'] is not None else None,
            gap=r['gap'] * HA_TO_EV,
            elapsed_seconds=r['elapsed_seconds'],
            convergence_data=r.get('convergence_data'),
            spectral_data=r.get('spectral_data'),
        )

    def compute_spectral_function(self, omega_grid: np.ndarray,
                                 orbital: int,
                                 broadening: float = 0.01) -> np.ndarray:
        """Compute spectral function A(ω) for given orbital.

        Args:
            omega_grid: Frequency grid (Hartree)
            orbital: Orbital index
            broadening: Broadening parameter (Hartree)

        Returns:
            Spectral function values
        """
        if self.gw_obj is not None and hasattr(self.gw_obj, 'compute_spectral_function'):
            return self.gw_obj.compute_spectral_function(
                omega_grid, orbital, broadening
            )
        else:
            # Simple Lorentzian approximation
            qp_e = self._last_result['qp_energies'][orbital]
            z = self._last_result['z_factors'][orbital]

            # Quasiparticle peak
            spectral = z * broadening / (
                (omega_grid - qp_e)**2 + broadening**2
            ) / np.pi

            return spectral

    def save_results(self, filepath: str):
        """Save GW results to HDF5 file."""
        import h5py

        if not hasattr(self, '_last_result'):
            raise RuntimeError("No calculation has been performed yet")

        r = self._last_result

        with h5py.File(filepath, 'w') as f:
            # Save arrays
            f.create_dataset('qp_energies', data=r['qp_energies'])
            f.create_dataset('z_factors', data=r['z_factors'])
            f.create_dataset('sigma_x', data=r['sigma_x'])
            f.create_dataset('sigma_c', data=r['sigma_c'])
            f.create_dataset('mo_energy', data=self.mo_energy)
            f.create_dataset('mo_coeff', data=self.mo_coeff)

            # Save metadata
            f.attrs['converged'] = r['converged']
            f.attrs['n_iterations'] = r['n_iterations']
            f.attrs['max_error'] = r.get('max_error', 0.0)
            f.attrs['homo'] = r['homo']
            f.attrs['lumo'] = r['lumo']
            f.attrs['gap'] = r['gap']
            f.attrs['ip'] = r['ip']
            if r['ea'] is not None:
                f.attrs['ea'] = r['ea']
            f.attrs['elapsed_seconds'] = r['elapsed_seconds']
            f.attrs['n_occ'] = self.n_occ
            f.attrs['auxbasis'] = self.config.auxbasis
            f.attrs['freq_int'] = self.config.freq_int
            f.attrs['nfreq'] = self.config.nfreq


def create_gw_driver_from_pyscf(mf, auxbasis='def2-tzvp-jkfit', config=None):
    """Create a GWDriver from a completed PySCF calculation.

    Args:
        mf: Completed PySCF mean-field calculation (RHF/RKS)
        auxbasis: Auxiliary basis set for RI/DF
        config: Optional GWConfig object

    Returns:
        GWDriver: Initialized driver ready for GW calculation
    """
    if not HAS_PYSCF:
        raise RuntimeError("PySCF is required to use this function")

    # Get MO data from PySCF
    mo_energy = mf.mo_energy
    mo_coeff = mf.mo_coeff
    mo_occ = mf.mo_occ
    n_occ = int(np.sum(mo_occ > 0))

    # Create GW driver
    driver = GWDriver(
        mo_energy=mo_energy,
        mo_coeff=mo_coeff,
        n_occ=n_occ,
        mol=mf.mol,
        auxbasis=auxbasis,
        config=config
    )

    return driver