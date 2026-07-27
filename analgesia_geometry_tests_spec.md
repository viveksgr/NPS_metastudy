# Analgesia geometry: implementation spec for the coding agent

Goal: test whether diverse analgesic interventions **converge behaviorally but diverge neurally**, and specifically whether the neural variance that distinguishes interventions lives in the **null space of the neural→behavioral readout** (the directions that don't change predicted relief). Every test below is a concrete, checkable operation on one subject-per-row table.

---

## 0. Data contract

One row = one subject (no subject ID column needed). Columns:

- `study_id` — **required.** This is the confounding unit *and* the resampling unit. Nothing below is valid without it.
- **Intervention dummies** — the one-hot group columns as designed for the LMM; use them as-is. For the classifiers in Tests 0–3, the multiclass target is simply which dummy is 1 in a row (argmax over the one-hot set) — no separate categorical column or reconstruction needed.
- `y` — subjective analgesia as an **effect size** (high–low pain contrast). Already in standardized units; confirm the sign so that **larger = more relief** and flip if not.
- Neural coordinate columns forming the vector `x` (`k = 3`): `C1`, `C2`, `C3`, with `C3` the positive-activation axis. Optional interpretable set for the robustness run: `NPS`, `SIIPS` (and, if wanted, all five together).

Sign convention: orient every coordinate so that larger = more analgesia before fitting, and record the sign applied to each. `C3` is the positive/activation direction; keep its polarity explicit since it enters the readout with the opposite intuition to a deactivation axis.

---

## 1. Configuration (lock these at the top of the script)

- `COORDS` — list of neural coordinate columns. Primary: `{C1, C2, C3}`. Robustness: `{NPS, SIIPS}` and, optionally, `{NPS, SIIPS, C1, C2, C3}`.
- `Y_TARGET` — `"effect_size"` (use `y` as given) or `"study_z"` (further z-score `y` within study). Default `"effect_size"`, since `y` is already a standardized effect size; keep `"study_z"` as a sensitivity check in case studies still differ in analgesia scale (affects how comparable the predicted-relief coordinate `a` is across studies).
- `N_BOOT` — bootstrap resamples over studies (e.g. 2000).
- `RNG_SEED`.

---

## 2. Preprocessing

1. **Harmonization.** If scanner/paradigm differences are a concern, either ComBat-harmonize `x` across `study_id`, or rely on study-level random effects downstream (§3). Do not skip this silently; log which was used.
2. **Standardize** each coordinate in `COORDS` (mean 0, unit variance) *using training-study statistics only* inside the CV loop (§7) — never fit scaling on held-out studies.
3. **Residual (within-group, within-study) covariance `Σ_resid`.** Remove `group × study` cell means from `x`, then compute the covariance of the residuals. This is the metric used for whitening. Fit on training studies only.
4. **Whitening.** Define `z = Σ_resid^(-1/2) · x`. All length/variance/angle comparisons (Tests 1 & 2) happen in `z`-space so that "orthogonal" is meaningful despite correlated, differently-scaled axes. Tests 0 & 3 don't need this (they condition on a scalar, then decode).
5. **Circularity.** The coordinate/ROI definitions are treated as fixed inputs — they've been shown robust to leave-one-study-out lesions, so no per-fold redefinition is required. (Standardization and `Σ_resid` are still fit train-only per §2.2–2.3; that's about scaling, not axis definition.)

---

## 3. Step A — estimate the readout map `w`

`w` must be the **within-subject, within-study** neural→behavior slope, not the cross-study ecological correlation. Fit a linear mixed model, pooled across all subjects (not per group):

```
y = w · x + (1 | study_id) + ε        # random intercept by study, minimum
```

- The fixed-effect slope vector is `w`.
- If a study has enough subjects, also fit random slopes; a `w` that varies by study is itself a reportable finding (readout not universal).
- Do this in whitened space (fit on `z`), or fit on `x` and transform — be consistent, since the projections below assume the metric.
- Expectation from prior results: the readout is low-rank — `y` tracks an NPS-like axis and is largely blind to the rest. So most of `C1`/`C2`/`C3` (and SIIPS, in the interpretable run) should fall largely in the kernel. That is the mechanism, not a nuisance. With `k = 3` the null space is 2-dimensional.

### Derived geometric objects
- Readout axis: `ŵ = w / ‖w‖`.
- Null space (kernel): `{v : wᵀv = 0}`, dimension `k − 1`.
- Iso-analgesia level sets: hyperplanes `wᵀx = c` (affine translates of the kernel).

---

## 4. Test 0 — model-free degeneracy  *(run and report this first)*

Uses **no** `w`; assumption-light, so it can't be argued with on geometric grounds.

1. Bin subjects into narrow bands of **actual** `y` (equal-width or equal-count; report both).
2. Within each band (subjects matched on *experienced* relief), train a classifier to predict intervention from `x` (and, if feasible, from the full brain map).
3. Two references:
   - **Floor (significance):** permutation null — shuffle intervention labels within band. Establishes that relief-matched decoding beats chance.
   - **Ceiling (benchmark):** the pooled decoder trained/evaluated across *all* pain levels, which is free to exploit between-group relief differences. To keep the comparison fair, subsample the pooled training set to the per-band sample size so the gap reflects lost relief information, not lost N.
4. Aggregate balanced accuracy across bands and read off the **within-band vs pooled gap**:
   - within-band ≈ pooled ⇒ intervention identity is relief-invariant (**same relief, different neural state** — the empirical heart of the claim).
   - within-band ≫ chance but ≪ pooled ⇒ part of decodability rode on relief magnitude; the residual is relief-invariant.
   - within-band ≈ chance ⇒ identity is magnitude; Test-0 falsifier.

The null-space geometry (Tests 1–3) is the *explanation* for a positive Test 0, not the claim itself.

---

## 5. Test 1 — variance-fraction / angle headline

For every intervention pair `(g, g′)` in whitened space:

1. Mean-difference vector `Δ = z̄_g − z̄_g′`.
2. Split:
   - on-axis `a = (ŵᵀΔ)²`
   - null-space `b = ‖Δ − (ŵᵀΔ)ŵ‖²`
3. Report **null-space fraction** `b / (a + b)`, averaged over pairs → the one-number thesis ("X% of between-intervention neural separation lies in the analgesia null space").
4. Also report the distribution of **angles** between `Δ` and `w`; predicted to pile up near 90°.

**Claim is `b ≫ a`, not `a = 0`.** Interventions genuinely differ in mean relief, so some on-axis separation is expected — state this so no one "refutes" a claim you didn't make.

---

## 6. Test 2 — decode from the null space

Per subject, in whitened space:

1. Readout coordinate `a_i = ŵᵀ z_i`.
2. Null-space residual `z_i^⊥ = z_i − a_i ŵ`.
3. Train three `group` classifiers (same estimator, same CV folds):
   - `A_full` — from full `z`
   - `A_readout` — from `a` alone
   - `A_null` — from `z^⊥` alone
4. **Prediction: `A_null ≈ A_full ≫ A_readout`** — intervention identity is recoverable from the directions that don't change relief, and poorly from the relief axis itself. With whitening this is effectively LDA-vs-readout.

---

## 7. Test 3 — level-set classification  *(the figure's engine)*

Same as Test 0, but condition on the **neurally-predicted** relief `a` instead of actual `y`:

1. Bin subjects by `a`.
2. Within each `a`-band (a level set), decode `group` vs permutation null.
3. This is what the money figure draws.

---

## 8. Cross-validation & inference wrapper (wraps Steps A and Tests 1–3)

- **Leave-one-study-out (LOSO):** estimate `w`, `Σ_resid`, standardization, and any circularity-sensitive axis on training studies; compute all projections/decoders on the held-out study. Study is the unit of confounding, so it must be the fold boundary.
- **Bootstrap over studies** (`N_BOOT`) for CIs on: null-space fraction (Test 1), `A_null / A_readout` ratio (Test 2), and Test 0/3 accuracies.
- Rationale to preserve: sampling noise in `w` leaks on-axis signal into the null projection; only out-of-study estimation keeps the null-vs-readout split honest.
- Test 0 still cross-validates its decoder over studies even though it uses no `w`.

**Robustness across coordinate sets:** run the entire pipeline on `COORDS = {C1, C2, C3}` and again on `{NPS, SIIPS}` (optionally the combined five). Null-space dominance surviving both ⇒ not an artifact of a particular set of scores.

---

## 9. The figure

Scatter of subjects, colored by intervention:
- **x-axis** = readout coordinate `a` (neurally predicted relief).
- **y-axis** = the single most group-discriminating direction *within* the null space — run LDA on `z^⊥`, take the top discriminant, project onto it.
- Vertical dashed lines = iso-analgesia level sets (fixed `a`).
- Target visual: clusters strung out **vertically within a level set** (e.g. Remi and Conditioning at the same `a`, pulled apart on `y`). Caption: *same relief, different coordinates.*
- Overlay each intervention's mean vector from the origin: magnitude = engagement, direction = mechanistic signature (the complex-number instinct generalized off the 2D special case).

---

## 10. Decision rules (report these explicitly)

- `b ≫ a` **and** `A_null ≈ A_full` **and** Test 0 positive → thesis holds: analgesia is a degenerate endpoint; behavior is a low-rank readout; intervention identity sits in its null space.
- `a ≈ b` **or** `A_readout ≈ A_full` → interventions differ mostly in **magnitude** along the relief axis; the geometry collapses to "some work harder." This is the falsifier — surface it early rather than after building more analyses on top.

---

## Decision to lock before running

**`Y_TARGET`** — effect-size `y` as-is vs additionally study-standardized. Since `y` is already an effect size, default to as-is; run study-standardized as a sensitivity check to confirm the predicted-relief coordinate `a` is comparable across studies. If the two disagree, the study-standardized version is the safer headline.
