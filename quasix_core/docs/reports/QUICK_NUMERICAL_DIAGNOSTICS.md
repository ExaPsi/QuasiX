# Quick Numerical Diagnostics for Σᶜ Debugging

## 🚨 Emergency Checklist

When Σᶜ values look wrong, check these in order:

### 1. Catastrophic Cancellation Test (30 seconds)

```python
# Extract components from your calculation
residue = result['residue_part']  # Pole contributions
integral = result['integral_part']  # Imaginary axis integral
total = result['sigma_c']

# Check HOMO
homo_idx = n_occ - 1
res_mag = abs(residue[homo_idx])
int_mag = abs(integral[homo_idx])
tot_mag = abs(total[homo_idx])

CF = max(res_mag, int_mag) / tot_mag

print(f"Cancellation Factor: {CF:.1f}×")

# Thresholds
if CF < 100:
    print("✅ OK - Normal cancellation")
elif CF < 1000:
    print("⚠️ WARNING - Losing 2-3 digits of precision")
else:
    print("🔴 CRITICAL - Severe cancellation! Reformulate!")
```

**Typical values:**
- CF ~ 10×: ✅ Expected for GW
- CF ~ 100×: ⚠️ Acceptable but monitor
- CF > 1000×: 🔴 Problem - check implementation

---

### 2. Denominator Sanity Check (1 minute)

```python
import numpy as np

# Check all denominators in GW calculation
denoms = []
for i, ei in enumerate(mo_energy):
    for j, ej in enumerate(mo_energy):
        if i != j:
            denom = abs(ei - ej) + eta
            denoms.append(denom)

denoms = np.array(denoms)
min_denom = np.min(denoms)
max_denom = np.max(denoms)
conditioning = max_denom / min_denom

print(f"Min denominator: {min_denom:.6f} Ha = {min_denom*27.211:.3f} eV")
print(f"Max denominator: {max_denom:.6f} Ha = {max_denom*27.211:.3f} eV")
print(f"Conditioning:    {conditioning:.2e}")

# Warnings
if min_denom < 1e-4:
    print("⚠️ Very small denominator - check for degeneracies")
if conditioning > 1e6:
    print("⚠️ High conditioning - results may be sensitive to eta")
```

**Typical values (H₂O):**
- Min: ~0.001 Ha (set by η)
- Max: ~25 Ha (HOMO-LUMO span)
- Conditioning: 10⁴ - 10⁵ (acceptable)

---

### 3. Quick Convergence Test (5 minutes)

```python
# Run with different grid sizes
n_freqs = [8, 16, 32, 64]
sigma_c_homo = []

for n in n_freqs:
    result = compute_gw(n_freq=n)
    sigma_c_homo.append(result['sigma_c'][homo_idx])
    print(f"n={n:2d}: Σᶜ(HOMO) = {sigma_c_homo[-1]:.6f} Ha")

# Check convergence
diffs = np.abs(np.diff(sigma_c_homo))
print("\nConvergence:")
for i in range(len(diffs)):
    print(f"  Δ({n_freqs[i]}→{n_freqs[i+1]}): {diffs[i]:.2e} Ha")

# Should decrease exponentially
if len(diffs) >= 2:
    ratio = diffs[0] / diffs[1]
    print(f"\nConvergence rate: {ratio:.1f}×")
    if ratio < 10:
        print("🔴 NOT CONVERGING - Check quadrature implementation!")
    else:
        print("✅ Exponential convergence detected")
```

**Expected:** Each doubling of n_freq should reduce error by 10-1000×

---

### 4. Eta Independence Test (5 minutes)

```python
# Test different eta values
etas = [0.001, 0.01, 0.1]
sigma_c_homo = []

for eta in etas:
    result = compute_gw(eta=eta)
    sigma_c_homo.append(result['sigma_c'][homo_idx])
    print(f"η={eta:.3f}: Σᶜ(HOMO) = {sigma_c_homo[-1]:.6f} Ha")

# Check variation
sigma_range = max(sigma_c_homo) - min(sigma_c_homo)
sigma_mean = np.mean(sigma_c_homo)
sensitivity = sigma_range / abs(sigma_mean)

print(f"\nSensitivity: {sensitivity:.1%}")

if sensitivity < 0.1:
    print("✅ Low eta dependence")
elif sensitivity < 0.5:
    print("⚠️ Moderate eta dependence")
else:
    print("🔴 TOO SENSITIVE - Check pole treatment!")
```

**Expected:** Variation < 10% for η ∈ [0.001, 0.1] Ha

---

### 5. Physical Sanity Checks (instant)

```python
# 1. Imaginary part should be small
im_part = np.abs(sigma_c.imag)
if np.max(im_part) > 0.1:
    print(f"⚠️ Large Im[Σᶜ]: {np.max(im_part):.6f} Ha")
    print("   Expected: < 0.01 Ha for molecules")

# 2. HOMO correlation should be negative
if sigma_c[homo_idx].real > 0:
    print("🔴 WRONG SIGN: Σᶜ(HOMO) > 0 (should be negative)")
    print("   Check: sign of W matrix elements")

# 3. QP shifts should be ~0.5-2 eV
qp_shift_ev = sigma_c[homo_idx].real * 27.211
if abs(qp_shift_ev) > 5.0:
    print(f"⚠️ Large QP shift: {qp_shift_ev:.2f} eV")
    print("   Typical: 0.5-2 eV for molecules")

# 4. Z-factor check (if available)
if 'z_factor' in result:
    z = result['z_factor'][homo_idx]
    if not (0.0 < z < 1.0):
        print(f"🔴 Unphysical Z-factor: {z:.3f}")
        print("   Must be in (0, 1)")
```

---

## 🔧 Common Issues & Fixes

### Issue: Cancellation Factor > 1000

**Symptoms:**
- Very large residue and integral, tiny total
- Results change dramatically with small parameter changes
- Numerical noise dominates

**Diagnosis:**
```python
print(f"Residue:  {residue:.6f} Ha")
print(f"Integral: {integral:.6f} Ha")
print(f"Total:    {total:.6f} Ha")
print(f"CF:       {CF:.0f}×")
```

**Fixes:**
1. Increase η (try 0.01 instead of 0.001)
2. Use Kahan summation (compensated arithmetic)
3. Consider extended precision (f128)
4. Reformulate to avoid subtraction

---

### Issue: Not Converging with n_freq

**Symptoms:**
- Σᶜ changes by > 0.01 Ha between n_freq=32 and 64
- No clear exponential decay pattern

**Diagnosis:**
```python
# Check actual integrand behavior
import matplotlib.pyplot as plt

omega_grid = np.linspace(0, 30, 100)
integrand = [compute_integrand(omega) for omega in omega_grid]

plt.plot(omega_grid, integrand)
plt.xlabel('ω (Ha)')
plt.ylabel('Integrand')
plt.savefig('integrand.png')
```

**Fixes:**
1. Check quadrature transformation (freq/mod.rs)
2. Verify Jacobian factor included
3. Increase ξ_max (try 50 Ha instead of 30)
4. Check for discontinuities in W(iω)

---

### Issue: Strong η Dependence

**Symptoms:**
- Σᶜ varies > 20% for η ∈ [0.001, 0.01]
- Results unstable

**Diagnosis:**
```python
# Plot Σᶜ vs η
etas = np.logspace(-4, -1, 20)
sigmas = [compute_gw(eta=eta)['sigma_c'][homo_idx] for eta in etas]

plt.semilogx(etas, np.real(sigmas))
plt.xlabel('η (Ha)')
plt.ylabel('Re[Σᶜ(HOMO)] (Ha)')
plt.savefig('eta_dependence.png')
```

**Fixes:**
1. Check pole treatment (residue calculation)
2. Verify ±iη signs for occupied/virtual
3. Check for near-degeneracies (Test 2 above)
4. Ensure contour doesn't cross poles

---

### Issue: Wrong Sign for Σᶜ(HOMO)

**Symptoms:**
- Σᶜ(HOMO) > 0 (positive)
- Expected: negative for correlation

**Diagnosis:**
```python
# Check components separately
print(f"Σˣ(HOMO):  {sigma_x[homo_idx]:.6f} Ha")
print(f"Σᶜ(HOMO):  {sigma_c[homo_idx]:.6f} Ha")
print(f"Vxc(HOMO): {vxc[homo_idx]:.6f} Ha")
print(f"QP shift:  {(sigma_x[homo_idx] + sigma_c[homo_idx] - vxc[homo_idx]):.6f} Ha")
```

**Fixes:**
1. Check W matrix sign (should have W_PQ > 0 on diagonal)
2. Verify factor of -2 in P⁰ (RPA formula)
3. Check contraction: Σᶜ = ∫ G W (both should be symmetric)
4. Verify spin factor (2× for closed-shell)

---

## 📊 Reference Values (H₂O, def2-SVP)

For sanity checks:

```
Orbital energies:
  HOMO (4): -0.48 Ha = -13.1 eV
  LUMO (5): +0.06 Ha = +1.6 eV
  Gap:       0.54 Ha = 14.7 eV

Exchange self-energy:
  Σˣ(HOMO):  -0.52 Ha = -14.2 eV
  Σˣ(LUMO):  -0.02 Ha = -0.5 eV

Correlation self-energy:
  Σᶜ(HOMO):  -0.05 to -0.15 Ha (typical range)
  Σᶜ(LUMO):  -0.01 to -0.05 Ha

Quasiparticle energies:
  IP(HOMO):  ~12.5 eV (experiment: 12.6 eV)
  EA(LUMO):  ~1.0 eV
```

---

## 🎯 Quick Decision Tree

```
Is Σᶜ(HOMO) positive?
├─ YES → 🔴 Sign error - check W matrix, P⁰ formula
└─ NO  → Continue

Is cancellation factor > 100?
├─ YES → 🔴 Numerical instability - increase η or use Kahan sum
└─ NO  → Continue

Does Σᶜ change > 1% from n_freq=32 to 64?
├─ YES → ⚠️ Not converged - increase n_freq or check quadrature
└─ NO  → Continue

Does Σᶜ vary > 20% for η ∈ [0.001, 0.01]?
├─ YES → 🔴 Pole treatment bug - check residue calculation
└─ NO  → Continue

Is |Σᶜ(HOMO)| > 0.5 Ha?
├─ YES → ⚠️ Unusually large - verify with PySCF
└─ NO  → ✅ Likely OK - proceed to validation
```

---

**Quick reference for:** QuasiX correlation self-energy debugging
**Full analysis:** Run `tests/validation/test_numerical_stability.py`
**Report:** See `docs/reports/NUMERICAL_STABILITY_ANALYSIS.md`
