"""Diagnose W matrix and correlation self-energy magnitude."""
import sys
sys.path.insert(0, '/home/vyv/Working/QuasiX/quasix/python')

import numpy as np
from pyscf import gto, scf, df
import quasix_core

# H2O setup
mol = gto.M(atom='O 0 0 0; H 0 1 0; H 0 0 1', basis='sto-3g', verbose=0)
mf = scf.RHF(mol).density_fit(auxbasis='weigend').run()

# Get data
mo_coeff = mf.mo_coeff
mo_energy = mf.mo_energy
mo_occ = mf.mo_occ
nocc = int(np.sum(mo_occ > 0))
nmo = len(mo_energy)
naux = mf.with_df.get_naoaux()

print("=" * 80)
print("DIAGNOSING W MATRIX AND CORRELATION MAGNITUDE")
print("=" * 80)
print(f"\nSystem: H2O / sto-3g")
print(f"  nmo = {nmo}, nocc = {nocc}, naux = {naux}")
print(f"  HF HOMO = {mo_energy[nocc-1]:.6f} Ha = {mo_energy[nocc-1] * 27.211:.3f} eV")
print(f"  HF LUMO = {mo_energy[nocc]:.6f} Ha = {mo_energy[nocc] * 27.211:.3f} eV")
print(f"  HF gap = {(mo_energy[nocc] - mo_energy[nocc-1]) * 27.211:.3f} eV")

# DF tensors
df_3c2e = mf.with_df._cderi
df_tensor = np.zeros((nmo, nmo, naux))
for P in range(naux):
    df_slice = df_3c2e[P]
    df_ao = np.zeros((mol.nao, mol.nao))
    idx = 0
    for i in range(mol.nao):
        for j in range(i+1):
            df_ao[i,j] = df_ao[j,i] = df_slice[idx]
            idx += 1
    df_tensor[:,:,P] = mo_coeff.T @ df_ao @ mo_coeff

# Metric
metric_ao = mf.with_df.get_2c2e()
V_chol = np.linalg.cholesky(metric_ao)

print(f"\nDF tensor stats:")
print(f"  Max |(mn|P)|: {np.abs(df_tensor).max():.6e}")
print(f"  Mean |(mn|P)|: {np.abs(df_tensor).mean():.6e}")
print(f"  Max diagonal |(mm|P)|: {np.abs([df_tensor[m,m,:] for m in range(nmo)]).max():.6e}")

print(f"\nMetric stats:")
print(f"  Max |V|: {np.abs(metric_ao).max():.6e}")
print(f"  Max |V_chol|: {np.abs(V_chol).max():.6e}")
print(f"  V_chol diagonal: {np.diag(V_chol)[:5]}")

# Simple frequency grid
n_freq = 4
omega_grid = np.linspace(0.0, 2.0, n_freq)

print(f"\n{'='*80}")
print("COMPUTING W(ω) SCREENING")
print(f"{'='*80}")
w_screened, stats = quasix_core.compute_screened_interaction_batch(
    mo_energy, mo_occ, V_chol, df_tensor, omega_grid
)

print(f"\nW(ω=0) statistics:")
W0 = w_screened[0]
print(f"  Shape: {W0.shape}")
print(f"  Max |W(0)|: {np.abs(W0).max():.6e}")
print(f"  Max diagonal |W(0)|: {np.abs(np.diag(W0)).max():.6e}")
print(f"  Mean diagonal W(0).real: {np.diag(W0).real.mean():.6e}")
print(f"  W(0) diagonal[:5]: {np.diag(W0).real[:5]}")

# Compare with bare Coulomb
print(f"\nW vs bare Coulomb V:")
print(f"  Max |V|: {np.abs(metric_ao).max():.6e}")
print(f"  Max |W|: {np.abs(W0).max():.6e}")
print(f"  |W|/|V| ratio: {np.abs(W0).max() / np.abs(metric_ao).max():.3f}")
if np.abs(W0).max() / np.abs(metric_ao).max() > 2.0:
    print(f"  ⚠️  WARNING: W is {np.abs(W0).max() / np.abs(metric_ao).max():.1f}x larger than V!")
    print(f"  This is UNPHYSICAL - W should be ~ V in magnitude")

print(f"\n{'='*80}")
print("COMPUTING CORRELATION SELF-ENERGY")
print(f"{'='*80}")
sigma_c, sigma_c_result = quasix_core.compute_correlation_self_energy(
    mo_energy, mo_occ, w_screened, omega_grid, mo_energy, df_tensor
)

print(f"\nΣ_c(ε) at MO energies:")
print(f"  Shape: {sigma_c.shape}")
print(f"  Σ_c[HOMO] = {sigma_c[0, nocc-1]:.6f} Ha = {sigma_c[0, nocc-1].real * 27.211:.3f} eV")
print(f"  Σ_c[LUMO] = {sigma_c[0, nocc]:.6f} Ha = {sigma_c[0, nocc].real * 27.211:.3f} eV")
print(f"  Max |Σ_c|: {np.abs(sigma_c).max():.6e} Ha = {np.abs(sigma_c).max() * 27.211:.3f} eV")

print(f"\nResidue vs Integral contributions [HOMO]:")
res_homo = sigma_c_result["residue_part"][0, nocc-1]
int_homo = sigma_c_result["integral_part"][0, nocc-1]
print(f"  Residue: {res_homo:.6f} Ha = {res_homo.real * 27.211:.3f} eV")
print(f"  Integral: {int_homo:.6f} Ha = {int_homo.real * 27.211:.3f} eV")
print(f"  Total: {sigma_c[0, nocc-1]:.6f} Ha = {sigma_c[0, nocc-1].real * 27.211:.3f} eV")
print(f"  Residue/Total: {res_homo.real / sigma_c[0, nocc-1].real:.2%}")

# Check for magnitude issue
expected_correction_eV = 4.0  # Typical GW correction for occupied orbital
actual_correction_eV = abs(sigma_c[0, nocc-1].real * 27.211)
if actual_correction_eV > 2 * expected_correction_eV:
    print(f"\n⚠️  MAGNITUDE ERROR DETECTED:")
    print(f"  Expected |Σ_c| ~ {expected_correction_eV} eV for occupied orbitals")
    print(f"  Actual |Σ_c| = {actual_correction_eV:.3f} eV")
    print(f"  Factor: {actual_correction_eV / expected_correction_eV:.1f}x TOO LARGE")

print(f"\n{'='*80}")
print("DIAGNOSIS SUMMARY")
print(f"{'='*80}")
print("\nPotential issues:")
if np.abs(W0).max() / np.abs(metric_ao).max() > 1.5:
    print("❌ W matrix is too large compared to bare Coulomb V")
    print("   → Check dielectric matrix inversion normalization")
if actual_correction_eV > 2 * expected_correction_eV:
    print("❌ Correlation self-energy magnitude is too large")
    print("   → Likely caused by W normalization error")
    print("   → Or frequency integration weights incorrect")
else:
    print("✅ Magnitudes appear reasonable")