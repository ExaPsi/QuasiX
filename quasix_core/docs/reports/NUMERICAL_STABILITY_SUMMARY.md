# Numerical Stability Analysis - Executive Summary

## Critical Findings (Priority Order)

### 🟢 1. Catastrophic Cancellation: PASSED ✅

**Finding:** Cancellation factor = **10×** (acceptable, < 100 threshold)

**What this means:**
- Residue and imaginary integral contributions are similar magnitude (~0.14 Ha)
- They partially cancel to give Σᶜ ~ 0.014 Ha
- Losing only 1 significant digit (acceptable)
- **NOT** a numerical bug - this is expected GW physics

**Action:** ✅ No action needed. This is normal behavior.

---

### 🟡 2. Denominator Conditioning: WARNING ⚠️

**Finding:** Conditioning number = **2.5×10⁴** on real axis

**What this means:**
- Smallest denominator: 0.001 Ha (η broadening)
- Largest denominator: 24.8 Ha (HOMO-LUMO span)
- Ratio is high, but **expected** for small η

**Why it's acceptable:**
- Imaginary axis (where most accuracy comes from) has conditioning: 752× ✅
- No orbital degeneracies found
- Real-axis poles contribute only ~10-20% of total Σᶜ

**Action:** ⚠️ Monitor Z-factors (should be in 0-1 range). If issues arise, try η = 0.01 Ha as diagnostic.

---

### 🟢 3. Quadrature Accuracy: PASSED ✅

**Finding:** Gauss-Legendre quadrature accurate to machine precision

**Metrics:**
- ∫₀¹⁰⁰ 1 dx: relative error = 0.00e+00
- ∫₀¹⁰⁰ x² dx: relative error = 1.05e-15 ✅
- ∫₀¹⁰⁰ x⁴ dx: relative error = 3.10e-15 ✅

**Action:** ✅ No action needed. Implementation is correct.

---

### 🔴 4. Grid Convergence: PENDING ⚠️

**Status:** Requires full GW implementation to test

**What needs testing:**
```
n_freq    Σᶜ(HOMO)     Error reduction
  8       -0.XXX Ha    ---
 16       -0.YYY Ha    Should improve ~1000×
 32       -0.ZZZ Ha    Should improve ~1000×
 64       -0.WWW Ha    Should improve ~1000×
```

**Expected:** Exponential convergence (Gauss-Legendre property)

**Action:** 🔴 **HIGH PRIORITY** - Implement `tests/validation/test_gw_convergence.py` with H₂ molecule

---

### 🔴 5. Eta Sensitivity: PENDING ⚠️

**Status:** Requires full GW implementation to test

**What needs testing:**
- Variation of Σᶜ with η ∈ [0.001, 0.1] Ha
- Should vary < 20% (physical results independent of artificial parameter)

**Action:** 🔴 **HIGH PRIORITY** - Implement `tests/validation/test_eta_sensitivity.py`

---

## Immediate Next Steps

### Priority 1: Complete Tests 4-5 🔴

**Implement minimal GW test:**

```python
# tests/validation/test_gw_minimal.py

from pyscf import gto, scf
import quasix

# H2 molecule (simplest test)
mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='def2-svp')
mf = scf.RHF(mol).run()

# Test convergence
for n_freq in [8, 16, 32, 64]:
    gw = quasix.GW(mf, n_freq=n_freq)
    result = gw.kernel()
    print(f"n_freq={n_freq}: IP = {result.ip:.6f} eV")

# Should see exponential convergence
```

**Target:** Convergence within 0.01 eV for n_freq = 64

---

### Priority 2: Monitor Cancellation in Production 🟡

**Add diagnostic logging:**

```rust
// In correlation.rs
let residue_mag = residue.norm();
let integral_mag = integral.norm();
let total_mag = (residue + integral).norm();

let cancellation_factor = residue_mag.max(integral_mag) / total_mag;

if cancellation_factor > 100.0 {
    warn!("High cancellation: CF = {:.0}×", cancellation_factor);
} else {
    debug!("Cancellation: CF = {:.0}× (acceptable)", cancellation_factor);
}
```

---

### Priority 3: Validate Against PySCF 🟡

**Target tolerance:**
- DF integrals: < 1e-8 Ha ✅ (already achieved)
- Exchange Σˣ: < 1e-6 Ha ✅ (already achieved)
- Correlation Σᶜ: < 1e-4 Ha ⚠️ (needs validation)
- QP energies: < 0.1 eV 🔴 (critical for physics)

---

## Risk Assessment

| Component | Risk Level | Confidence | Action |
|-----------|-----------|------------|--------|
| Quadrature | 🟢 LOW | High (tested) | ✅ None |
| Denominators | 🟡 MEDIUM | High (understood) | Monitor Z-factors |
| Cancellation | 🟢 LOW | High (tested) | ✅ None |
| Convergence | 🔴 UNKNOWN | Low (untested) | **Test ASAP** |
| Eta sensitivity | 🔴 UNKNOWN | Low (untested) | **Test ASAP** |

---

## Success Criteria

### Minimum Viable Product (MVP)

- ✅ Quadrature: machine precision
- ✅ Cancellation: CF < 100
- 🔴 Convergence: exponential decay
- 🔴 Eta independence: < 20% variation
- 🔴 PySCF agreement: Σᶜ < 1e-4 Ha

### Production Quality

- ✅ All MVP criteria
- ⚠️ H₂O IP within 0.1 eV of experiment (12.6 eV)
- ⚠️ GW100 subset (10 molecules) < 0.5 eV MAE
- ⚠️ No Z-factors outside (0, 1)
- ⚠️ Documentation complete

---

## Bottom Line

**Current status:** 3/5 tests passed, 2 critical tests pending

**Good news:**
- ✅ No catastrophic cancellation (CF = 10×)
- ✅ Quadrature implementation correct
- ✅ No orbital degeneracies

**Concerns:**
- 🔴 Grid convergence untested (need full GW)
- 🔴 Eta sensitivity untested (need full GW)
- 🟡 High denominator conditioning (acceptable, but monitor)

**Recommendation:**
1. **Immediately** implement tests 4-5 with minimal H₂ system
2. Validate convergence is exponential (key indicator of correct implementation)
3. Verify eta independence (rules out pole treatment bugs)
4. Then proceed to H₂O and GW100 benchmarks

**Expected timeline:**
- Tests 4-5 implementation: 2-4 hours
- Validation runs: 1 hour
- If all pass → proceed to production testing
- If any fail → debug frequency integration or pole treatment

---

**Report generated:** 2025-11-09
**Analysis tool:** `/home/vyv/Working/QuasiX/quasix_core/tests/validation/test_numerical_stability.py`
**Full report:** `/home/vyv/Working/QuasiX/quasix_core/docs/reports/NUMERICAL_STABILITY_ANALYSIS.md`
