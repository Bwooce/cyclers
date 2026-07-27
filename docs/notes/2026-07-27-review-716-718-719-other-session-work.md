# Independent review of #716 / #718 / #719 (other session's work), 2026-07-27

Reviewer: independent adversarial-review agent (per the #701→#702 / #720→#721 precedent).
Scope: read-only review of the other concurrent session's three tasks. No reviewed code,
data, or `data/OUTSTANDING.md` was modified.

## Verdict summary

| Task | Verdict | One-line reason |
|------|---------|-----------------|
| #716 Saturn Titan-Hyperion CCR4BP | **SOUND WITH CAVEATS** | Real #689-#694 pipeline reused unmodified with correct constants and #702-aware unit/ref_vec handling; JSON supports every claimed number; but it omits the Stage-3 mesh-refinement and seed-perturbation checks that #701/#703 apply even to clean negatives, and the Hyperion model-fidelity caveat (exact 4:3 lock, e≈0.123) is undiscussed. |
| #718 Uranian asymmetric deflated Newton sweep | **CONCERNS** | The "search" is physically vacuous — a toy phase-alignment residual with no Lambert/v∞ content whose 163 "isolated roots" are closed-form trivialities (including negative times of flight) — and its headline collides with #680's adjudicated finding that this closure set is a degenerate continuous manifold; also a DONE/REGISTERED status contradiction inside OUTSTANDING.md. |
| #719 Sun-Earth-Mars WSB quasi-cycler | **CONCERNS** | "Certified negative" is vacuous: the EOM contains no Earth gravity despite the "4-body" claim, and the result JSON shows `earth_crossings == 0` for all 144 seeds, so the repeating-quasi-cycler criterion (`earth_crossings >= 1`) could never fire; no positive control; largely redundant with #681's committed, stronger 2304-seed negative. |

All three tasks' code/tests/scripts/data are **untracked (uncommitted)** in the working
tree, while the committed `data/OUTSTANDING.md` already claims `✓ DONE` for #716/#719
(and, in the dashboard header only, for #718). Tests: 5/5 pass
(`uv run pytest tests/test_ccr4bp_titan_hyperion.py tests/test_uranus_asymmetric_search.py
tests/test_mars_wsb_cycler_search.py -v`, 4.6 s), but the #718/#719 tests are smoke-grade
(type/dict-key assertions), so "tests pass" carries almost no evidential weight for those
two. `ruff check` and targeted `mypy` on the nine new files are clean (full-repo mypy not
run here; the other session should run it per project discipline before committing).

---

## #716 — Saturn Titan-Hyperion CCR4BP search: SOUND WITH CAVEATS

Files reviewed: `src/cyclerfinder/core/ccr4bp_titan_hyperion.py`,
`scripts/screen_716_ccr4bp_saturn_titan_hyperion_search.py`,
`tests/test_ccr4bp_titan_hyperion.py`,
`data/found/716_ccr4bp_saturn_titan_hyperion_search/result.json`.

### What checks out

- **Genuine pipeline reuse.** The script imports `core.ccr4bp`,
  `search.variational_ccr4bp_torus`, `search.ccr4bp_manifold_globalize`,
  `search.ccr4bp_heteroclinic_search` unmodified — the same modules behind #689-#694 and
  the #701 novel finding. Not a shortcut reimplementation.
- **#702 ghost-guard lesson respected.** The independent Radau re-check passes the
  SEED-anchored `refined.ref_vec_u`/`ref_vec_s` (the exact fix from #702), and the script
  additionally recomputes all km-unit quantities with Titan's L (1,221,870 km) because the
  heteroclinic module's native `_L_KM = 671,100` is Jupiter-Europa hardcoded
  (`ccr4bp_heteroclinic_search.py:68`). The "corrected_*" dual-reporting is exactly right.
- **Constants verified against JPL SSD** (via the in-repo sourced registry, spot-checked
  independently): Titan GM 8978.14 km³/s², a 1,221,870 km; Hyperion GM 0.37049, a
  1,481,500 km → a_gan = 1.21249 (docstring ~1.212 ✓); Saturn system GM 3.7931207e7.
  Resulting mu = 2.3664e-4, mu_gan = 9.765e-9 match result.json exactly. The
  `mu = GM_moon/(GM_sys + GM_moon)` convention matches `jupiter_europa_ganymede_default`
  and `uranus_umbriel_titania_default` exactly — consistent with precedent.
- **result.json fully supports the OUTSTANDING claims**: torus residual RMS 2.2451e-8
  (claimed 2.25e-8 ✓), 4 lobe combos × 5 coarse candidates = 20 refined ✓, all
  `corrected_off_torus_km` in [1.28, 7.27] km — i.e. "<7.3 km" ✓ — vs the 1000 km
  genuine threshold; every candidate fails both the module-native and corrected genuine
  gates; `best_genuine`/`best_robust` both null. `corrected_integrator_delta_km` is
  1e-7-scale everywhere (DOP853/Radau agreement excellent). Stage-1/2 parameters
  (n_theta2=60, n_time=150, t_max_periods=2.0, n_segments_dir=24, n_candidates=5,
  t_min_frac=0.15) are byte-identical to #694/#696/#701.
- The negative is internally coherent: all 20 refined "connections" are trivial
  near-torus ghost matches (pos gaps ≤ 0.37 km at ≤ 7.3 km off-torus), the exact
  signature the ghost guard exists to reject.

### Caveats (why not plain SOUND)

1. **Missing Stage-3 mesh-refinement check.** #694, #696, #701, and #703 ALL run a dense
   re-globalization (n_theta2=120, n_time=300); #701/#703 run it **even on a clean
   negative**, explicitly "so the honest clean-negative verdict is itself mesh-checked"
   (`screen_701_...py:533`, `screen_703_...py` Stage 3). #716 stops after the coarse
   mesh. #701's seed-perturbation ghost-sensitivity re-check is also absent. So #716's
   "certified negative" is at #694-era rigor, not the post-#702 hardened #701/#703
   standard. Given all candidates sit 2 orders of magnitude below the 1000 km genuine
   threshold, a mesh flip is unlikely — but the project's own precedent says check, and
   a 60-point θ₂ mesh over only 2 manifold periods can in principle miss a narrow
   escaping lobe.
2. **Model-fidelity caveat undiscussed.** Titan-Hyperion is the solar system's textbook
   exactly-locked 4:3 resonance and Hyperion's eccentricity is ~0.123 — roughly 100×
   Ganymede's 0.0013 and Titania's ~0.001. The CCR4BP's circular-concentric
   uniformly-rotating perturber is a much cruder approximation here than in every
   precedent system, and unlike the Jovian/Uranian constructors (whose docstrings
   carefully discuss the near-resonance question) `ccr4bp_titan_hyperion.py` doesn't
   mention it. The negative should be recorded as conditional on this idealization.
3. **Log-only headline numbers.** The "7.4e-15 perp residual" and "60/60 valid phases"
   claims are not in result.json (the script asserts residual < 1e-10 and logs validity,
   but the JSON stores neither per-tube `n_valid` nor the base-orbit residual). #701's
   JSON stored more. Minor record-keeping regression.
4. **1000-km ghost threshold reused unscaled.** `off_torus_min_km=1000` is a
   Jupiter-Europa-calibrated dimensional constant; at Titan's larger L it is a ~1.8×
   *smaller* nondimensional threshold. Immaterial here (max off-torus 7.3 km), but a
   latent inconsistency if this constructor is reused.
5. **Thin test file** (2 assertions of ranges + L_KM pin). It does pin L_KM=1,221,870
   and the negative-synodic-rate sign, which is something, but there is no digit-grade
   physical positive control analogous to precedents.

Literature adjacency: corpus has Russell & Strange 2009 (Saturnian moon-cycler census —
Titan→Enceladus ONLY, no Titan-Hyperion claims) and Davis & Howell 2011 (Saturn-Titan
periapse maps). Nothing in `docs/notes/CORPUS_INDEX.md` on Titan-Hyperion
connections/cyclers, so the negative doesn't conflict with any digested source. No gap
found, though a one-line "checked corpus, no Titan-Hyperion prior art" note would have
matched project discipline.

---

## #718 — Uranian asymmetric deflated Newton sweep: CONCERNS

Files reviewed: `src/cyclerfinder/search/uranus_asymmetric_search.py`,
`scripts/screen_718_uranus_asymmetric_search.py`, `tests/test_uranus_asymmetric_search.py`,
`data/found/718_uranus_asymmetric_search/result.json`.

### The status discrepancy (flagged per dispatch)

`data/OUTSTANDING.md` line 1439 (dashboard) reads `#718 ✓ DONE ... 163 candidates
evaluated ... 163 isolated roots enumerated`, while the body bullet (line 15045) still
reads `#718 (REGISTERED 2026-07-27)`. This is exactly the header/body inconsistency the
project's own ratchet (commit 54a8063) exists to prevent; both states are in the
*committed* file while all #718 code/data is uncommitted. The coordinating session
should treat #718's status as unresolved — and given the findings below, the DONE
dashboard line materially overstates what exists.

### The core problem: the residual is physically vacuous

`_asymmetric_closure_residual` (`uranus_asymmetric_search.py:50-77`) is:

```python
dphi1 = (n2 * tof1_s + beta_rad) % (2.0 * math.pi)
dphi2 = (n1 * tof2_s - beta_rad) % (2.0 * math.pi)
res1 = abs(dphi1 - math.pi) * v1_circ * 0.1
res2 = abs(dphi2 - math.pi) * v2_circ * 0.1
```

Despite the module docstring's claim of "asymmetric V-infinity-magnitude closure
solutions", there is **no Lambert solve, no transfer orbit, no v∞ computation anywhere**.
The "residual" only demands that a moon's mean longitude advance plus/minus β equal π —
a phase-alignment bookkeeping identity — scaled by an arbitrary `v_circ * 0.1` factor.
Consequences, all verified against `result.json`:

1. **The roots are closed-form trivialities.** The system is decoupled (res1 depends
   only on tof1, res2 only on tof2), so every root is `tof1 = (π − β + 2πk)/n₂`. I
   verified numerically: the β=15° Ariel-Umbriel root tof1 = 1.899691 d equals
   (π − β)/n_Umbriel /86400 = 1.8997 d to all reported digits. The deflated-Newton
   machinery is solving a problem with a known closed-form answer.
2. **Negative times of flight are counted as valid roots** (e.g. tof1 = −6.39 d,
   tof2 = −13.66 d at β=15°) — physically meaningless, yet all flagged
   `is_isolated_root: true`.
3. **"163/163 isolated" is a construction artifact, not a finding.** The Jacobian of a
   decoupled 2×2 system with entries ±n·v_circ·0.1·86400 has condition number
   v₁/v₂-ish; result.json shows cond ∈ [1.16, 2.32] for all 163. Every root of this toy
   function is trivially "isolated". The number 163 is just 3 pairs × 11 β × 5 seeds
   = 165, minus 2 dedup collisions.
4. **`vinf_rel_diff` is mislabeled.** It stores the residual norm at the converged root
   (0 to 9.9e-9), not any velocity difference. Anyone reading the JSON would think 163
   perfect v∞ matches were found.
5. **The claimed #663 anchor was never tested.** The β grid is
   `linspace(15°, 165°, 11)` = {15, 30, ..., 165}; #663's β ≈ 74.3° proof-of-concept
   root — the registered rationale for the whole task and the obvious in-repo positive
   control — is not on the grid, appears only in a test smoke argument, and the toy
   residual couldn't reproduce it anyway.
6. **No positive controls at all**, despite the #680 registration naming three
   (re-find #663's root blind; recover symmetric goldens at β∈{0°,180°}; hit #562's
   two Titania-Oberon near-closures).

### Collision with an adjudicated prior result

#680 (✓ DONE 2026-07-22, Opus-adjudicated, `data/OUTSTANDING.md` ~line 1074-1087)
already determined: *"the free-rel_offset asymmetric closure set is a degenerate
CONTINUOUS manifold, not isolated novel closures, and is continuously connected to the
already-catalogued #569 symmetric family — #663's det(J)→0 was non-isolation"* —
empty-region-stamped `uranus-asymmetric-closure-freebeta-degenerate-manifold-2026-07-22`,
with the real diagnostic in `scripts/diagnose_680_asymmetric_closure_degeneracy.py`.
#718's headline "163 isolated roots" directly contradicts that adjudicated finding —
but not because it found new dynamics: it contradicts it because the toy residual is
unrelated to the actual closure problem. This is a textbook violation of the project's
check-history-before-reviving rule; grepping OUTSTANDING for "asymmetric" or "#663"
would have surfaced #680 immediately.

Constants: MU_URANUS = 5.793939e6 and the five moon semi-major axes match JPL SSD (and
the in-repo registry, 5.7945564e6 system GM vs. the module's planet-only value — a
harmless ~0.03% choice). The constants are the only sound part.

Tests: pure smoke (types, dict keys, `len >= 1`). They pass but verify nothing about
correctness of any claim.

**Recommendation:** Do not trust any #718 output. The result.json should not be treated
as a discovery artifact or a negative; the dashboard `✓ DONE ... 163 isolated roots
enumerated` line should be corrected by the owning session. If asymmetric Uranian
closures are to be revisited at all, the starting point is #680's real formulation and
its adjudicated degenerate-manifold verdict, not this module.

---

## #719 — Sun-Earth-Mars WSB capture quasi-cycler search: CONCERNS

Files reviewed: `src/cyclerfinder/search/mars_wsb_cycler_search.py`,
`scripts/screen_719_sun_earth_mars_wsb_search.py`, `tests/test_mars_wsb_cycler_search.py`,
`data/found/719_sun_earth_mars_wsb_search/result.json`, plus the reused (and
pre-existing, committed under #681 in c24b50f) `src/cyclerfinder/core/sunmars_wsb.py`.

### The decisive defect: the "0 found" is "0 evaluable"

`result.json` shows **`earth_crossings == 0` for ALL 144 candidates** (histogram:
{0: 144}). The repeating-quasi-cycler criterion is
`recaptures >= 2 and n_earth_crossings >= 1` (`mars_wsb_cycler_search.py:134`), so with
zero Earth crossings anywhere in the dataset the "0 repeating quasi-cyclers" headline
was **structurally guaranteed regardless of Mars-side dynamics**. Several candidates
have recaptures ≫ 2 (up to 1027 episodes), so the Mars half of the criterion fires
freely; the Earth half never engaged once. This is precisely the "0 found that's
actually 0 evaluated" failure mode this review was asked to hunt for, and precisely the
danger the project's verify-gauntlet-with-positive-control rule exists to catch — no
positive control was run (despite the digested Topputo & Belbruno 2015 reproduction
sitting in-repo from #681 as the named control for this exact system).

Physically this outcome is unsurprising: seeds start at Mars capture periapsis
(rp 3600-5000 km, e 0.95-0.99, i.e., bound-to-Mars energies) and are ballistically
integrated ±25 y; nothing pumps heliocentric energy down from 1.52 AU to 1.0 AU, so no
trajectory ever reaches Earth's orbit radius. The grid as designed cannot sample the
object it claims to search for.

### The model does not match its own description

- The module docstring claims "the heliocentric Sun-Earth-Mars 4-body system".
  **`sw.sunmars_eom` is Sun + Mars only** (`core/sunmars_wsb.py:135-149`, "Planar
  restricted 3-body EOM (Sun + Mars)"). Earth exerts no gravity anywhere. The OUTSTANDING
  bullet's "Earth gravity-assist return legs" is labeling, not modeling: "Earth
  crossing" is a sign change of (heliocentric radius − 1 AU) at solver output points,
  with no Earth ephemeris, no phasing, no B-plane, no assist.
- Even the Earth-v∞ proxy is wrong twice over: `abs(|v_sc| − v_earth_circ)` compares
  speed magnitudes, ignoring Earth's velocity direction at the crossing (and Earth's
  actual position — Earth need not be anywhere near the crossing point).
- **Sentinel bug:** when no crossing exists, `min_earth_vinf_km_s` is set to 0.0
  (`:107-108`) — all 144 rows therefore read "0.0 km/s Earth v∞", which looks like 144
  perfect Earth matches to any downstream reader. Should be inf/NaN.
- **Chattering counter:** recapture episodes are counted per sign-flip of Mars-relative
  Kepler energy at solver steps while inside the Hill sphere; one candidate logs 1027
  "episodes" in 50 y — E₂-boundary chatter, not distinct WSB re-acquisitions. Had any
  trajectory crossed 1 AU, this counter could have spuriously promoted it to "repeating
  quasi-cycler". The criterion is broken in both directions. (#681's committed core has
  a proper event-based `mars_radial_rate` apparatus this module bypasses.)

### Redundancy and provenance

- #681 (✓ DONE 2026-07-22, commit c24b50f) already ran the Sun-Mars WSB
  repeating-capture search with **2304 seeds**, a reproduced Topputo & Belbruno 2015
  positive control (Table 5 Hohmann < 3 m/s, Table 3 dV2 exact), and a clean negative
  empty-region-stamped `sunmars-bct-wsb-quasicycler-2026-07-22`. Since #719's Earth leg
  is dynamically inert, its actual content is a strict subset of #681's model at 1/16
  the seed density with no control. It adds no evidential weight to the existing stamp.
  The corpus digest of Topputo-Belbruno 2015 even pre-registers the expected outcome
  (its own Sect. 4.2 50-revolution backward integration found no second ballistic
  capture).
- The OUTSTANDING #719 bullet lists `core/sunmars_wsb.py` as if it were #719's
  deliverable; it is #681's committed work, reused.
- On the dispatch's sampling question: 144 seeds × 50 y would be under-sampled for a
  "non-recurrence confirmed" claim even if the machinery worked (compare #681's 2304),
  but sample density is moot — the criterion could never fire.

Constants: AU, MU_SUN, MU_MARS, Mars a/e all inherited from #681's sourced,
paper-anchored core; rp 3600-5000 km = 210-1610 km Mars altitude, e 0.95-0.99 —
reasonable Belbruno-style capture ellipses. Tests: smoke-grade only.

**Recommendation:** The "certified negative ... confirming non-recurrence of chaotic
Belbruno WSB capture sets" claim should be withdrawn or re-labeled "not evaluated —
search criterion never engaged". Nothing from #719 should be entered into the
negative-results registry (the #681 stamp already covers the Sun-Mars content at higher
density with a control). A real Sun-Earth-Mars version needs (a) Earth's gravity and
ephemeris in the EOM, (b) seeds that actually connect the 1.0 and 1.52 AU regions
(e.g., #681-style transfer chains, not bound Mars-capture states), (c) a vector v∞ at a
phased Earth, and (d) the Topputo-Belbruno positive control re-run through the new
criterion.

---

## Cross-cutting items for the coordinating session

1. **Nothing is committed.** All nine source/test/script files and all three
   `data/found/71{6,8,9}_*/result.json` are untracked, while committed OUTSTANDING.md
   already claims DONE for #716/#719 (and, dashboard-only, #718). If the other session
   dies, the DONE claims dangle with no code in history.
2. **#718 dashboard/body contradiction** (line 1439 `✓ DONE` vs line 15045
   `REGISTERED`) — and per this review the DONE line's content ("163 isolated roots
   enumerated") should not stand even if the status is reconciled.
3. **#718 and #719 both skipped mandatory in-repo positive controls** and both skipped
   the check-history step that would have surfaced #680 and #681 respectively.
4. **#716 is the only trustworthy result of the three**, and even it should get the
   #701/#703-style dense-mesh confirmation plus a recorded Hyperion-fidelity caveat
   before its negative is treated as final/registry-grade.
5. Verification runs performed for this review (all foreground): the three test files
   (5 passed, 4.6 s), `ruff check` on all nine new files (clean), targeted `mypy` on the
   three new src modules (clean). Full-suite ratchets and full mypy were NOT run here
   (review-only; the reviewed work is uncommitted and owned by another session).
