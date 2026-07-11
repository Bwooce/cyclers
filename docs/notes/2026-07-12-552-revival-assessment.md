# #552 revival assessment — does the #571 Titan–Iapetus result reopen the 3D/inclined-releg genome build?

**Date:** 2026-07-12
**Task:** adjudicate whether task #552 (a 3D/inclined-releg extension to the moontour
discovery genome, previously *SCOPED, NOT GREENLIT* on 2026-07-10) should be revived on the
strength of #571's 187 gate-passing Titan–Iapetus candidates.
**Scope constraint:** read/analysis/calculation only. No capability built, no
`catalogue.yaml` / `empty_regions.jsonl` edit, no new sweep. Analysis reuses the existing
#571 data plus `core/flyby.py::max_bend` and `core/satellites.py`'s sourced constants.

---

## 1. The question in one line

#552 was killed because its payoff claim (reopen Amalthea / Neptune-Triton empty regions)
failed a back-of-envelope flyby-bend check: Amalthea is mass-limited regardless of geometry,
and Triton modeled *correctly* (retrograde) makes its empty verdict **stronger**, not
reopened — the coplanar-prograde model was being *generous* to Triton and hiding real
infeasibility. #571 then produced 187 coplanar Titan–Iapetus candidates, every one carrying
an explicit `iapetus_inclination_caveat` (Iapetus's real ~15.5° inclination gives a
~0.85 km/s out-of-plane relative velocity, comparable to the entire ~0.93 km/s coplanar V∞
budget). **Does that ~0.85 km/s penalty behave like Triton (kills everything) or is it
genuinely different (signal survives)?** The whole revival decision turns on this.

## 2. The same back-of-envelope check that killed Triton, applied to Titan–Iapetus

Constants (sourced, `core/satellites.py`): Iapetus GM = 120.51511 km³/s², mean R = 734.3 km,
safe_alt = 100 km → r_p,min = 834.3 km, a = 3,561,700 km. Saturn GM = 3.7931207e7 km³/s².

- **Iapetus orbital speed:** v = √(GM_Sat/a) = **3.263 km/s**.
- **Out-of-plane penalty at i = 15.5°:** v·sin(15.5°) = **0.872 km/s** (confirms the caveat's
  ~0.85). This is the velocity mismatch at Iapetus's node — the only point a near-coplanar
  (Titan-plane) trajectory naturally meets Iapetus, where positions match (both at z ≈ 0) but
  the out-of-plane velocity mismatch is **maximal**. So adding the full 0.872 km/s in
  quadrature is the *worst-case* / most-conservative correction, not an average one.
- **Gate ceiling:** the two-sided #324 gate (`physical_sanity.py`, applied to every encounter)
  requires ≥5° bend at Iapetus. `max_bend` gives 5° at Iapetus V∞ = **1.780 km/s** (matches
  Fable's 1.78). So the max *coplanar* Iapetus V∞ that can survive the penalty at the node is
  √(1.780² − 0.872²) = **1.551 km/s**.

**Real-corrected survivorship over all 187 flagged candidates** (real V∞ = √(V∞_coplanar² +
0.872²), re-run through `max_bend`, keep only those ≥5° at every Iapetus encounter):

| quantity | value |
|---|---|
| coplanar Iapetus V∞ across the 187 | 0.982 – 1.756 km/s (median 1.440) |
| **survive 5° gate after worst-case penalty** | **149 / 187 (80%)** |
| real-corrected Iapetus bend | 4.15 – 8.87° (median 5.41°) |
| best cluster (coplanar V∞ ≈ 0.98) real bend | **~8.9°** |

**Inclination sensitivity** (worst-case node penalty in quadrature):

| relative inclination | penalty | survive 5° gate |
|---|---|---|
| 8.0° (local Laplace plane) | 0.454 km/s | 182 / 187 |
| 12.0° | 0.678 km/s | 166 / 187 |
| **15.5° (Iapetus-to-Titan-plane, central estimate)** | **0.872 km/s** | **149 / 187** |
| 17.3° (to ecliptic — pessimistic) | 0.970 km/s | 127 / 187 |

## 3. This is genuinely DIFFERENT from Triton — not the same pattern

| | Triton (#554, killed) | Titan–Iapetus (#571) |
|---|---|---|
| nature of the correction | **retrograde** — a ~180° direction flip | **prograde, inclined ~15.5°** — a plane change only |
| effect on encounter V∞ | additive: forces V∞ to ~8.5 km/s (~2× orbital speed) | quadrature: √(V∞² + 0.87²), a ~15–40% increase |
| resulting bend at flyby body | collapses to ~1.5° (below the 5° floor) | 4.15–8.87°, median 5.41° (mostly above the floor) |
| survivors after correction | **0** | **149 / 187 (127–182 across 8–17°)** |
| what the coplanar model was doing | *generous* — hiding real infeasibility | *mildly* generous — a modest, absorbable penalty |

The logic that killed #552's Triton payoff — "coplanar model generous, real geometry kills
it" — **does not transfer to Titan–Iapetus.** Iapetus is prograde with real mass (GM = 120,
comparable to Rhea/Oberon) and the inclination is a bounded plane-change cost, not a velocity
doubling. The correction is fundamentally milder, and ~80% of the candidates absorb even the
worst-case node penalty at the physical bend gate.

## 4. Distribution of margin — not a fragile cluster at the floor

The Triton-failure signature would be all candidates piled right at the 5° minimum, collapsing
en masse under any real perturbation. That is **not** what these 187 look like. *Before* the
penalty their coplanar Iapetus bends span 5.13–14.97° (median 7.47°): only 27 sit in [5,6°),
105 in [6,8°), and 55 at ≥8°. There is a genuine population with comfortable margin, and the
best-margin cluster sits at the *analytically-predicted minimum-achievable* Iapetus V∞
(≈0.98 km/s, the 0.93 floor Fable derived), which is exactly where a real family would be
expected to concentrate. This is closer to #566's robust Uranian representatives (15–45° bend)
than to a fragile floor-hugging artifact.

## 5. What the survivorship does and does NOT prove

**Does prove:** the max-bend physical-feasibility gate — a *necessary* condition — survives
the real inclination for the large majority of candidates. Unlike Triton, real 3D geometry
does not analytically annihilate this pair.

**Does NOT prove:** that a *closed* 3D resonant trajectory exists. Passing the two-sided bend
gate is necessary-not-sufficient ([[feedback_verify_gauntlet_with_positive_control]] /
[[feedback_orbit_closure_discipline]] — "it matched!" is the danger signal). Two unmodelled
gaps remain:
1. **Closure under real 3D geometry.** The coplanar Lambert solution is a *seed*, not a
   solution, once Iapetus is placed on its real inclined orbit. Whether a real 3D transfer
   closes near that seed with residual comparable to the coplanar value is an empirical
   question the bend check cannot answer.
2. **Eccentricity gap (Fable flag #4, pre-registered on #571).** Titan e ≈ 0.0288 → ±0.16 km/s
   velocity modulation, ~3× the 0.05 km/s residual floor; Iapetus e ≈ 0.028 similarly. The
   idealized→real gap and #568-style duty cycle will be materially worse than the Uranian
   family. This is a downstream V-gauntlet concern, not a capability-build blocker, but it
   caps the realistic upside.

## 6. Cost of the build now — narrow Titan–Iapetus, not the general capability

The original #552 kill correctly rejected the *general* n-body-arbitrary-inclination build
(multi-week, and only one pair now motivates it). #552's own entry already identified the
narrow re-scope: "re-scope around Iapetus/Saturn … genuinely small (3–5 days) since
`core/lambert.py` already supports 3D/retrograde branch selection."

**Verified in code:** `core/lambert.py::lambert` takes full 3D `r1`, `r2` position vectors and
derives the transfer plane from their cross product (`cross_z` branch logic, lines ~688–691),
with an explicit `prograde: bool` branch selector. The Lambert machinery is **already fully
3D-capable** — it is not the bottleneck. The only coplanar assumption lives upstream in state
*generation*: `search/discovery_campaign.py::_moon_state()` hardcodes `z = 0`, and
`core/satellites.py` carries no inclination/node/direction field. So a narrow Titan–Iapetus 3D
re-evaluation is exactly: give Iapetus a real inclined orbit (i, node) in a one-off state
generator, feed the resulting real 3D positions to the already-3D Lambert, recompute V∞
vectors and re-gate. That confirms #552's "small lift" claim — and it means an even cheaper
**single-candidate hand-check is a day-scale probe, not a build.**

## 7. Recommendation: CHEAP-CHECK-FIRST (allocate #572), then narrow build only if it passes

**Do NOT revive the general multi-week #552.** Nothing here motivates an
n-body-arbitrary-inclination capability; exactly one pair is in question.

**Do NOT declare Titan–Iapetus a family either** — the survivorship is a necessary-condition
pass, not a closure proof, and the eccentricity gap caps the upside.

**Do the cheap intermediate first (new task #572):** a bounded, day-scale 3D **closure** probe
on the top 2–3 Titan–Iapetus candidates by real-corrected margin (coplanar Iapetus V∞ ≈ 0.98,
real bend ≈ 8.9°, e.g. rel_offset 18° / tof_scale 1.15 / residual 5.6e-5). Place Iapetus on
its real inclined orbit (i ≈ 15.5°, real node) in a throwaway state generator, feed real 3D
positions to `core/lambert.py` (already 3D + prograde/retrograde branch-aware), and attempt an
actual 3D Lambert closure near each coplanar seed. **Decision gate:** does a real 3D solution
exist with residual near the coplanar value *and* every encounter still clearing the bend
gate? This is the decisive question the bend check alone cannot settle, and it is the same
necessary-condition-then-closure discipline the project already runs.

- **If #572 finds a genuine 3D closure:** *then* scope the narrow Titan–Iapetus 3D corrector
  (#552's own re-scoped ~3–5 day estimate, reusing `core/lambert.py`'s 3D support — **not** the
  general build), followed by the standard Opus adjudication + Fable second opinion + V-gauntlet
  before any writeback, with the eccentricity/duty-cycle caveat honored throughout.
- **If #572 finds nothing closes under real 3D geometry:** stamp Titan–Iapetus as
  conditionally-empty-pending-nothing (the bend gate passed but no closure survives) and close
  the pair — a clean, cheap negative, exactly as [[feedback_never_give_up_reproducing_papers]]
  vs. [[project_negative_results_registry]] intend for a *novel* (not published) search.

This spends ~1 day to answer the only load-bearing open question before committing to any
3–5 day build, mirroring the discipline that has repeatedly served this project: the bend gate
is the cheap necessary filter (done here — **PASSES**, unlike Triton), and the 3D closure probe
is the next cheap increment before a capability investment.

## 8. Uncertainty flags for the Fable second-opinion pass

1. **Relative inclination 15.5°:** I used Iapetus-to-Saturn-equator (~15.47°, well-established)
   as a proxy for Iapetus-to-Titan-plane (Titan is only ~0.35° to the equator, so the relative
   inclination is ~15.1–15.5° depending on node longitudes). Confident within ~1°, and the
   conclusion is robust across the whole 8–17° band (127–182 survivors). To the *local Laplace
   plane* Iapetus is only ~8° — which would be even more favorable — so 15.5° is the
   conservative choice, not a cherry-pick.
2. **Worst-case node penalty:** I add the **full** 0.872 km/s out-of-plane component in
   quadrature (encounter at Iapetus's node). A real 3D-optimized trajectory could tilt the
   spacecraft's own plane to *split* the Titan↔Iapetus plane change (Titan bends hugely, so it
   can absorb a small tilt cheaply), reducing the Iapetus penalty below 0.872. So actual
   survivorship under a proper 3D correction is likely **≥ 149**, i.e. my estimate is a floor.
3. **Necessary-not-sufficient:** the entire §2–4 result is a bend-*feasibility* pass, not a
   closure proof. This is deliberately why the recommendation is a #572 closure probe rather
   than "the family is real." If Fable disagrees anywhere, it is most likely here — is the bend
   gate survivorship strong enough evidence to justify even the day-scale #572 probe? My
   judgment: yes, because it cleanly separates this case from Triton (0 survivors) and the cost
   is a day, but this is the judgment call worth a second opinion.
