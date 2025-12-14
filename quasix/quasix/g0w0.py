"""High-level Python interface for QuasiX G₀W₀ calculations.

This module provides a simple, PySCF-compatible interface for running
single-shot GW (G₀W₀) calculations using the QuasiX Rust kernel.

Example:
    >>> from pyscf import gto, scf
    >>> from quasix.gw import run_g0w0
    >>>
    >>> mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g')
    >>> mf = scf.RHF(mol).run()
    >>> result = run_g0w0(mol, mf)
    >>> print(f"IP = {result['ip']:.3f} eV")
"""

import numpy as np
from typing import Dict, Optional, Any
import warnings

# Physical constants
HA_TO_EV = 27.211386245988  # Hartree to eV conversion


def run_g0w0(
    mol: Any,
    mf: Any,
    auxbasis: str = 'def2-tzvp-jkfit',
    nfreq: int = 32,
    eta: float = 0.01,
    verbose: int = 1
) -> Dict[str, Any]:
    """
    Run G₀W₀ calculation on a PySCF molecule.

    This function performs a single-shot GW calculation using the QuasiX
    Rust kernel. It extracts DF tensors from PySCF, computes the self-energy
    once at the HF/DFT eigenvalues, and returns quasiparticle energies.

    Args:
        mol: PySCF Mole object
        mf: Completed PySCF mean-field calculation (RHF/RKS/UHF)
        auxbasis: Auxiliary basis set for RI/DF (default: 'def2-tzvp-jkfit')
        nfreq: Number of frequency points for integration (default: 32)
        eta: Broadening parameter in Hartree (default: 0.01 Ha ≈ 0.27 eV)
        verbose: Verbosity level (0=quiet, 1=normal, 2=debug)

    Returns:
        Dictionary containing:
        - 'qp_energies': Quasiparticle energies in Hartree [n_mo]
        - 'z_factors': Renormalization factors [n_mo]
        - 'sigma_x': Exchange self-energy in Hartree [n_mo]
        - 'sigma_c_real': Real part of correlation self-energy in Hartree [n_mo]
        - 'sigma_c_imag': Imaginary part of correlation self-energy in Hartree [n_mo]
        - 'sigma_c': Complex correlation self-energy in Hartree [n_mo] (for validation)
        - 'corrections': QP corrections (E_QP - ε_HF) in Hartree [n_mo]
        - 'vxc_diag': Exchange-correlation potential in Hartree [n_mo]
        - 'ip': Ionization potential in eV
        - 'ea': Electron affinity in eV
        - 'gap': Fundamental gap in eV
        - 'homo_index': Index of HOMO orbital
        - 'lumo_index': Index of LUMO orbital
        - 'W_iw': Screened interaction W(iω) in Hartree [nfreq, n_aux, n_aux] (for validation)
        - 'P0_iw': RPA polarizability P₀(iω) in Hartree [nfreq, n_aux, n_aux] (for validation)
        - 'freqs_iw': Imaginary frequency grid in Hartree [nfreq]
        - 'weights_iw': Quadrature weights [nfreq]

    Raises:
        ImportError: If PySCF is not installed
        RuntimeError: If DF tensor construction fails
        RuntimeError: If G₀W₀ calculation fails

    Example:
        >>> mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g')
        >>> mf = scf.RHF(mol).run()
        >>> result = run_g0w0(mol, mf)
        >>> print(f"HOMO IP: {result['ip']:.3f} eV")
        >>> print(f"Gap: {result['gap']:.3f} eV")
    """
    # Import PySCF (required)
    try:
        from pyscf import gto, scf, df as pyscf_df
    except ImportError as e:
        raise ImportError(
            "PySCF is required for G₀W₀ calculations. "
            "Install with: pip install pyscf"
        ) from e

    # Import GWDriver for DF tensor construction
    try:
        import sys
        import os
        # Add quasix/python directory to path
        # __file__ is in quasix/quasix/g0w0.py, so go up one level then to python/
        quasix_python_path = os.path.abspath(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'python'))
        if quasix_python_path not in sys.path:
            sys.path.insert(0, quasix_python_path)
        from gw_driver import GWDriver, GWConfig
    except ImportError as e:
        raise ImportError(
            f"QuasiX GWDriver not found at {quasix_python_path}. "
            "Please ensure QuasiX is properly installed."
        ) from e

    # Import quasix_core functions (individual S3-1 through S3-6 components)
    try:
        import quasix_core
    except ImportError as e:
        raise ImportError(
            "QuasiX core library not available. "
            "Please rebuild: cd quasix_core && maturin develop --release"
        ) from e

    if verbose > 0:
        print("=" * 70)
        print(" QuasiX G₀W₀ Calculation")
        print("=" * 70)
        print(f"  System: {mol.atom}")
        print(f"  Basis: {mol.basis}")
        print(f"  Auxiliary basis: {auxbasis}")
        print(f"  Number of orbitals: {mf.mo_energy.size}")
        print(f"  Number of electrons: {mol.nelectron}")
        print("=" * 70)

    # Extract MO data
    mo_energy = np.asarray(mf.mo_energy, dtype=np.float64)
    mo_coeff = np.asarray(mf.mo_coeff, dtype=np.float64)
    mo_occ = np.asarray(mf.mo_occ, dtype=np.float64)
    n_occ = int(np.sum(mo_occ > 0))
    n_mo = len(mo_energy)

    # Create GW driver for DF tensor construction
    config = GWConfig(
        auxbasis=auxbasis,
        freq_int='cd',
        nfreq=nfreq,
        eta=eta,
        verbose=verbose
    )

    driver = GWDriver(
        mo_energy=mo_energy,
        mo_coeff=mo_coeff,
        n_occ=n_occ,
        mol=mol,
        mf=mf,
        auxbasis=auxbasis,
        config=config
    )

    # Build DF tensors using PySCF
    if verbose > 0:
        print("\nStep 1/3: Building DF tensors from PySCF...")

    try:
        iaP, miP, v_aux = driver._build_df_tensors()
    except Exception as e:
        raise RuntimeError(
            f"Failed to build DF tensors: {e}\n"
            "This may indicate:\n"
            "  - Incompatible auxiliary basis for your system\n"
            "  - Linear dependency in basis sets\n"
            "  - Numerical instability\n"
            "Try a different auxiliary basis (e.g., 'weigend' or 'weigend+etb')"
        ) from e

    # Extract exact vxc (CRITICAL: must use exact integrals, NOT DF!)
    if verbose > 0:
        print("\nStep 2/3: Extracting exact v_xc (no DF approximation)...")

    try:
        vxc_diag = driver.extract_vxc_exact()
    except Exception as e:
        raise RuntimeError(
            f"Failed to extract v_xc: {e}\n"
            "v_xc MUST be computed with exact integrals to avoid cancellation.\n"
            "Without exact v_xc, G₀W₀ results will have errors > 1 eV!"
        ) from e

    # Prepare data for Rust
    mo_energy_contig = np.ascontiguousarray(mo_energy, dtype=np.float64)
    mo_occ_contig = np.ascontiguousarray(mo_occ, dtype=np.float64)
    vxc_diag_contig = np.ascontiguousarray(vxc_diag, dtype=np.float64)

    # Build full DF tensor (n_mo, n_mo, n_aux) from blocks
    naux = driver.df_tensors['naux']
    nvir = n_mo - n_occ

    df_3c_mo = np.zeros((n_mo, n_mo, naux), dtype=np.float64, order='C')

    # Fill from miP (all orbitals × occupied)
    miP_3d = miP.reshape(n_mo, n_occ, naux)
    df_3c_mo[:, :n_occ, :] = miP_3d

    # Fill from iaP (occupied × virtual)
    iaP_3d = iaP.reshape(n_occ, nvir, naux)
    df_3c_mo[:n_occ, n_occ:, :] = iaP_3d

    # v_aux is identity matrix for Cholesky (metric inverse)
    df_metric_inv = np.ascontiguousarray(v_aux, dtype=np.float64)

    if verbose > 0:
        print("\nStep 3/7: Running G₀W₀ calculation using quasix_core components...")
        print(f"  DF tensor shape: {df_3c_mo.shape}")
        print(f"  Metric shape: {df_metric_inv.shape}")

    # Run G₀W₀ using individual quasix_core functions (S3-1 through S3-6)
    try:
        # S3-1: Create imaginary frequency grid (PySCF convention)
        if verbose > 0:
            print(f"\n  S3-1: Creating {nfreq}-point imaginary frequency grid (PySCF rational mapping)...")

        # PySCF uses rational transformation to [0, ∞) with x0 = 0.5
        # This is CRITICAL - linear transformation gives wrong frequencies!
        # Reference: pyscf/pbc/gw/gw_ac.py::_get_scaled_legendre_roots()
        x0 = 0.5  # PySCF scaling parameter

        if verbose > 0:
            print(f"  Using PySCF rational transformation: x0 = {x0}")

        grid = quasix_core.gauss_legendre_grid(nfreq)
        freqs_gl = grid['nodes']  # GL nodes on [-1, 1]
        weights_gl = grid['weights']  # GL weights

        # Transform to imaginary axis [0, ∞) using PySCF rational mapping
        # freqs = x0 * (1 + x) / (1 - x)  maps [-1, 1] → [0, ∞)
        # weights *= 2 * x0 / (1 - x)^2  (Jacobian)
        # PySCF formula: gw/gw_ac.py::_get_scaled_legendre_roots()
        freqs_iw = x0 * (1.0 + freqs_gl) / (1.0 - freqs_gl)  # [nfreq]
        weights_iw = weights_gl * 2.0 * x0 / ((1.0 - freqs_gl)**2)  # [nfreq]

        # S3-2: Compute RPA polarizability P₀(iω)
        if verbose > 0:
            print("  S3-2: Computing RPA polarizability P₀(iω)...")
        # Extract (ia|P) block from full DF tensor
        ia_P = df_3c_mo[:n_occ, n_occ:, :].reshape(n_occ * (n_mo - n_occ), naux)
        p0_iw_complex = quasix_core.compute_polarizability_p0(
            ia_P.reshape(n_occ, n_mo - n_occ, naux),
            mo_energy_contig,
            mo_occ_contig,
            freqs_iw
        )

        # CRITICAL: P₀ is now complex-valued from Rust
        # On imaginary axis, imaginary part should be < 1e-10
        max_imag_p0 = np.max(np.abs(np.imag(p0_iw_complex)))
        if verbose > 1 or max_imag_p0 > 1e-10:
            print(f"  DEBUG: max|Im(P₀)| = {max_imag_p0:.3e} Ha")
            if max_imag_p0 > 1e-10:
                import warnings
                warnings.warn(f"P₀ has large imaginary part: {max_imag_p0:.3e} Ha")

        # Extract real part for further processing (should be nearly identical)
        p0_iw = np.real(p0_iw_complex)

        # S3-3: Compute screened interaction W(iω)
        if verbose > 0:
            print("  S3-3: Computing screened interaction W(iω)...")
        w_iw_complex = quasix_core.compute_screening_w(
            p0_iw_complex,  # Pass complex P₀ to Rust
            df_metric_inv,
            freqs_iw
        )

        # CRITICAL: W is now complex-valued from Rust
        # On imaginary axis, imaginary part should be < 1e-10
        max_imag_w = np.max(np.abs(np.imag(w_iw_complex)))
        if verbose > 1 or max_imag_w > 1e-10:
            print(f"  DEBUG: max|Im(W)| = {max_imag_w:.3e} Ha")
            if max_imag_w > 1e-10:
                import warnings
                warnings.warn(f"W has large imaginary part: {max_imag_w:.3e} Ha")

        # Extract real part for further processing (should be nearly identical)
        w_iw = np.real(w_iw_complex)

        if verbose > 1:
            print(f"  DEBUG: W_iw shape = {w_iw.shape}")
            print(f"  DEBUG: W_iw[0, 0, 0] = {w_iw[0, 0, 0]:.6e} Ha")

        # S3-4: Compute exchange self-energy Σˣ
        if verbose > 0:
            print("  S3-4: Computing exchange self-energy Σˣ...")

        # CRITICAL: For G₀W₀ with HF starting point, Σˣ = V_xc
        # This is because HF exchange IS the exchange self-energy.
        # PySCF does the same thing for G₀W₀ calculations.
        #
        # General formula: Σˣ_nn = -Σ_i Σ_PQ (ni|P) v⁻¹_PQ (Q|ni)
        # But for HF: Σˣ = V_xc (already computed exactly from HF)
        #
        # This avoids the DF approximation error in exchange.
        sigma_x = vxc_diag  # Use exact HF exchange (no DF error!)

        # S3-5: Compute correlation self-energy Σᶜ
        if verbose > 0:
            print("  S3-5: Computing correlation self-energy Σᶜ using PySCF CD formula...")

        # CRITICAL: PySCF CD formula computes Green's function PER TARGET ORBITAL
        # This is because the sign of the shift depends on target orbital's position.
        #
        # Reference: pyscf/gw/gw_cd.py lines 174-176:
        #   emo = omega - 1j*eta*sign - mo_energy
        #   g0 = wts[None,:]*emo[:,None] / ((emo**2)[:,None]+(freqs**2)[None,:])
        #   sigma = -einsum('mw,mw',g0,Wmn)/np.pi
        #
        # where 'sign' depends on target orbital position relative to Fermi level

        # PYTHON-ONLY PROTOTYPE (Phase 1)
        # This implements the PySCF CD formula exactly in pure Python
        # for fast iteration and validation. Rust implementation in Phase 2.

        # Prepare containers
        sigma_c_complex = np.zeros(n_mo, dtype=np.complex128)
        lpq = df_3c_mo.transpose(2, 0, 1)  # [naux, nmo, nmo]
        naux = lpq.shape[0]

        # Chemical potential (Fermi level) - midpoint between HOMO and LUMO
        ef = (mo_energy[n_occ-1] + mo_energy[n_occ]) / 2.0

        # PySCF CD integration parameters
        omega = 0.0  # Static self-energy (ω=0)

        if verbose > 1:
            print(f"    Fermi level (μ): {ef:.6f} Ha")
            print(f"    Broadening (η): {eta:.6f} Ha")
            print(f"    Integration frequency: ω = {omega:.6f} Ha")
            print(f"    Computing target-dependent Green's function for each orbital...")

        # Compute Σᶜ for each orbital separately (TARGET-DEPENDENT Green's function)
        for n_target in range(n_mo):
            # Determine sign based on target orbital position relative to Fermi level
            # PySCF convention: sign = +1 for virtual (n >= n_occ), -1 for occupied (n < n_occ)
            # This matches pyscf/gw/gw_cd.py line 232-234
            sign_target = 1.0 if n_target >= n_occ else -1.0

            # PySCF CD formula: emo = omega - 1j*eta*sign - mo_energy
            # For ALL orbitals m, computed with target orbital n's sign
            emo = omega - 1j*eta*sign_target - mo_energy  # [n_mo] array (complex)

            # Green's function: g0 = wts * emo / (emo^2 + freqs^2)
            # Shape: [nfreq, n_mo] where each row is G(iω_k) for all orbitals
            g0 = np.zeros((nfreq, n_mo), dtype=np.complex128)
            for ifreq in range(nfreq):
                weight = weights_iw[ifreq]
                freq = freqs_iw[ifreq]
                # g0[ifreq, :] = weight * emo / (emo^2 + freq^2)
                denominator = emo**2 + freq**2  # Complex + real
                g0[ifreq, :] = weight * emo / denominator

            # Compute W(iω) × G(iω) for this target orbital
            # W is computed from P₀(iω) at each frequency
            # Formula: Σᶜ_n = -1/π Σ_w Σ_m W_nm(iω) G_m(iω)

            # Accumulate Σᶜ_n over frequencies
            sigma_c_n = 0.0 + 0.0j

            for ifreq in range(nfreq):
                # Get P₀ at this frequency: [naux, naux]
                p0_w = p0_iw[ifreq, :, :]

                # Build W in MO basis using PySCF convention
                # CRITICAL: PySCF computes Pi = [I - P₀]^{-1} - I (screened Coulomb HOLE)
                # Reference: pyscf/gw/gw_cd.py line 163
                #
                # Step 1: Compute epsilon^{-1} = [I - P₀]^{-1}
                epsilon_inv = np.linalg.inv(np.eye(naux) - p0_w)  # [naux, naux]

                # Step 2: Pi = epsilon^{-1} - I (screened Coulomb hole, NOT full W!)
                # This is the INTERACTION part: W - v = [ε^{-1} - I] v
                pi = epsilon_inv - np.eye(naux)  # [naux, naux]

                # Step 3: Transform to MO basis
                # W_nm = Σ_PQ Lpq[P, n, m] * Pi[P, Q] * Lpq[Q, n, m]
                # where Pi is the screened Coulomb hole (NOT bare Coulomb!)
                #
                # We only need W_nm for fixed n=n_target:
                # W_m = Σ_PQ Lpq[P, n_target, m] * Pi[P, Q] * Lpq[Q, n_target, m]

                # Extract Lpq[:, n_target, :] -> [naux, nmo]
                lpq_n = lpq[:, n_target, :]  # [naux, nmo]

                # W_m = einsum('Pm,PQ,Qm->m', lpq_n, Pi, lpq_n)
                # Efficient implementation: tmp = Pi @ lpq_n -> [naux, nmo]
                tmp = np.dot(pi, lpq_n)  # [naux, nmo]
                # Then: w_nm = Σ_P lpq_n[P, m] * tmp[P, m]
                w_nm = np.sum(lpq_n * tmp, axis=0)  # [nmo] (element-wise multiply then sum over P)

                # Get Green's function for all orbitals at this frequency
                g_m = g0[ifreq, :]  # [nmo]

                # Accumulate: Σᶜ_n += Σ_m W_nm(iω) G_m(iω)
                sigma_c_n += np.sum(w_nm * g_m)  # Scalar

            # Apply -1/π normalization (PySCF convention)
            sigma_c_n *= -1.0 / np.pi

            sigma_c_complex[n_target] = sigma_c_n

            if verbose > 1 and (n_target == n_occ-1 or n_target == n_occ):
                orbital_label = "HOMO" if n_target == n_occ-1 else "LUMO"
                print(f"    Σᶜ_{orbital_label} (n={n_target}): {sigma_c_n:.6e} Ha (sign={sign_target:+.0f})")

        if verbose > 0:
            print(f"  ✓ Σᶜ computed using PySCF CD formula (target-dependent, Python prototype)")
        sigma_c_real = sigma_c_complex.real
        sigma_c_imag = sigma_c_complex.imag

        # For diagnostics: save last Green's function (not used in calculation)
        # This is the Green's function for the last orbital computed
        green_iw = g0.T  # [nfreq, nmo] - last computed Green's function

        # DEBUG: Print Σᶜ values
        if verbose > 0:
            print(f"  DEBUG: Σᶜ_HOMO = {sigma_c_complex[n_occ-1]:.6e} Ha")
            print(f"  DEBUG: Σᶜ_LUMO = {sigma_c_complex[n_occ]:.6e} Ha")
            print(f"  DEBUG: Σˣ_HOMO = {sigma_x[n_occ-1]:.6e} Ha")
            print(f"  DEBUG: V_xc_HOMO = {vxc_diag_contig[n_occ-1]:.6e} Ha")

        # S3-6: Solve quasiparticle equation
        if verbose > 0:
            print("  S3-6: Solving quasiparticle equations...")
        qp_result = quasix_core.solve_qp_equation(
            mo_energy_contig,
            sigma_x,
            sigma_c_real,
            vxc_diag_contig
        )
        qp_energies = qp_result['qp_energies']
        z_factors = qp_result['z_factors']

    except Exception as e:
        raise RuntimeError(
            f"G₀W₀ calculation failed: {e}\n"
            "This may indicate:\n"
            "  - Numerical instability in self-energy evaluation\n"
            "  - Invalid input data\n"
            "  - Bug in implementation"
        ) from e

    # Compute corrections (QP - HF)
    corrections = qp_energies - mo_energy

    # Compute observables from QP energies
    homo_idx = n_occ - 1
    lumo_idx = n_occ

    ip = -qp_energies[homo_idx] * HA_TO_EV  # eV (negative of HOMO QP energy)
    ea = -qp_energies[lumo_idx] * HA_TO_EV  # eV (negative of LUMO QP energy)
    gap = (qp_energies[lumo_idx] - qp_energies[homo_idx]) * HA_TO_EV  # eV

    # Print summary
    if verbose > 0:
        print("\n" + "=" * 70)
        print(" G₀W₀ Results Summary")
        print("=" * 70)
        print(f"  HOMO energy:")
        print(f"    HF:  {mo_energy[homo_idx]:.6f} Ha = {mo_energy[homo_idx]*HA_TO_EV:>8.3f} eV")
        print(f"    QP:  {qp_energies[homo_idx]:.6f} Ha = {qp_energies[homo_idx]*HA_TO_EV:>8.3f} eV")
        print(f"    ΔE:  {corrections[homo_idx]:.6f} Ha = {corrections[homo_idx]*HA_TO_EV:>8.3f} eV")
        print(f"  LUMO energy:")
        print(f"    HF:  {mo_energy[lumo_idx]:.6f} Ha = {mo_energy[lumo_idx]*HA_TO_EV:>8.3f} eV")
        print(f"    QP:  {qp_energies[lumo_idx]:.6f} Ha = {qp_energies[lumo_idx]*HA_TO_EV:>8.3f} eV")
        print(f"    ΔE:  {corrections[lumo_idx]:.6f} Ha = {corrections[lumo_idx]*HA_TO_EV:>8.3f} eV")
        print(f"  Ionization potential: {ip:.3f} eV")
        print(f"  Electron affinity:    {ea:.3f} eV")
        print(f"  Fundamental gap:      {gap:.3f} eV")
        print(f"  Z-factor (HOMO): {z_factors[homo_idx]:.4f}")
        print(f"  Z-factor (LUMO): {z_factors[lumo_idx]:.4f}")
        print("=" * 70)

    # Return complete result dictionary
    return {
        'qp_energies': qp_energies,  # Ha
        'z_factors': z_factors,
        'sigma_x': sigma_x,  # Ha
        'sigma_c_real': sigma_c_real,  # Ha
        'sigma_c_imag': sigma_c_imag,  # Ha
        'sigma_c': sigma_c_complex,  # Complex Σᶜ [nmo] (NEW - for validation)
        'corrections': corrections,  # Ha
        'vxc_diag': vxc_diag,  # Ha
        'ip': ip,  # eV
        'ea': ea,  # eV
        'gap': gap,  # eV
        'homo_index': homo_idx,
        'lumo_index': lumo_idx,
        'homo_idx': homo_idx,  # Backward compatibility
        'lumo_idx': lumo_idx,  # Backward compatibility
        'mo_energy_hf': mo_energy,  # Ha (for reference)
        'n_occ': n_occ,
        'nfreq': nfreq,
        'eta': eta,
        'auxbasis': auxbasis,
        # NEW: Intermediate values for validation
        'W_iw': w_iw,  # Screened interaction [nfreq, naux, naux] (real, Ha)
        'W_iw_complex': w_iw_complex,  # Complex W(iω) for validation [nfreq, naux, naux]
        'P0_iw': p0_iw,  # RPA polarizability [nfreq, naux, naux] (real, Ha)
        'P0_iw_complex': p0_iw_complex,  # Complex P₀(iω) for validation [nfreq, naux, naux]
        'green_iw': green_iw,  # Green's function G₀(iω) [nfreq, nmo] (complex, Ha)
        'freqs_iw': freqs_iw,  # Imaginary frequency grid [nfreq] (Ha)
        'weights_iw': weights_iw,  # Quadrature weights [nfreq] (dimensionless)
        'df_3c_mo': df_3c_mo,  # Full DF tensor [nmo, nmo, naux] (real, Ha)
    }
