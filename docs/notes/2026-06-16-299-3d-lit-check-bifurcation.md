# #299 Phase 3 — 3D family lit-check + bifurcation tracking

Date: 2026-06-16
Issue: #299 (Phase 3 of #291)
Predecessor: Phase 2 family tracer at commit `dcc8c4a`
              (`data/family_296_3d_em_11.jsonl`, 265 members)

## Goal

Phase 2 produced a 265-member 3D Earth-Moon (1,1) family at
`data/family_296_3d_em_11.jsonl`. Phase 3 had two jobs:

1. **Literature-check** every member against `KNOWN_CORPUS` (specifically the
   Antoniadou-Voyatzis 2018 spatial CR3BP anchor added at commit `568d8a4`).
   Expected: ALL rediscovery; the goal is a clean confirmation the technique
   works without claiming false novelty.
2. **Bifurcation tracking** the two fold-points flagged in the Phase 2 walk
   (T-fold around step +11; z0-fold around step -3) using the existing
   `bifurcation_detector` machinery, and switch families on any
   period-multiplying bifurcation.

## Part A — Literature check

Script: `scripts/run_299_lit_check_3d_family.py`
Output: `data/lit_check_299_3d_family.jsonl` (265 rows + 1 header)

Each of the 265 members was lifted into a `CandidateSignature` whose structural
footprint is the Earth-Moon CR3BP cycler / periodic-orbit fingerprint:

    primary = "Earth"
    sequence = ("Moon",)
    resonances = ("spatial-cr3bp",)
    vinf_per_encounter_kms = (0.5,)   # nondim CR3BP regime; below the
                                       # "high V-inf" 6 km/s gate

The check ran against a deterministic offline corpus mirroring the live web
hits for the four Earth-Moon CR3BP anchors in `KNOWN_CORPUS`
(Antoniadou-Voyatzis 2018 spatial CR3BP, Roberts-Tsoukkas-Ross 2026
multi-orbiter cyclers, Braik-Ross 2026 orbital networks, plus a noise hit).

### Result

| Status        | Count | Anchor                                              |
|---------------|-------|-----------------------------------------------------|
| `published`   | 265   | Antoniadou-Voyatzis spatial resonant periodic orbits in CR3BP (2018) |
| `not-found`   | 0     | —                                                   |
| `inconclusive`| 0     | —                                                   |

All 265 members hit the **Antoniadou-Voyatzis 2018** anchor with confidence
0.85 — the family fingerprint (Earth-Moon CR3BP, spatial periodic orbit,
out-of-plane structure) matches the corpus entry exactly.

**Verdict: clean rediscovery. The 3D-family tracer reproduces a published
spatial-CR3BP family. No novelty claims warranted.**

This is the EXPECTED outcome and the whole point of running the check before
any catalogue admission — Phase 2's banner warning ("Likely-rediscovery of
Antoniadou-Voyatzis 2018 spatial CR3BP") is confirmed.

## Part B — Bifurcation tracking

Script: `scripts/run_299_bifurcation_track_3d_family.py`
Output: `data/family_296_3d_subfamilies_299.jsonl` (4 bracket records + 1 header)

### Bracket detection

The Phase 2 JSONL stores Floquet multipliers per member. The Part B script
reads those directly (no re-integration at the scan stage) and finds
adjacent-pair brackets where the minimum distance to a primitive k-th root of
unity (k ∈ {2,…,6}) crosses a tolerance band of 0.01.

**12 brackets** were found, **all classified as Neimark-Sacker** (a
complex-conjugate pair sitting on the unit circle |λ| ≈ 1, sweeping through
the primitive k-th root). No period-doubling (k=2) brackets — the family does
not feature a real multiplier crossing -1.

Bracket inventory (sample, sorted by k):

| k | parent step | parent T | parent z0 | dist to primitive root | nearest λ |
|---|---|---|---|---|---|
| 3 | +86  | 9.545 | -0.130 | 0.023 | -0.520 - 0.854i |
| 3 | +88  | 9.525 | -0.127 | 0.009 | -0.492 + 0.871i |
| 4 | +8   | 10.279 | -0.236 | 0.203 | 0.202 + 0.979i |
| 4 | +9   | 10.286 | -0.233 | 0.003 | 0.003 + 1.000i |
| 4 | +110 | 9.308 | -0.095 | 0.092 | -0.092 - 0.996i |
| 4 | +112 | 9.289 | -0.091 | 0.002 | 0.002 + 1.000i |
| 5 | +23/+24/+53/+55 | 9.85-10.17 | -0.16 to -0.20 | 0.002-0.028 | on the 5th-root ring |
| 6 | -65  | 9.566 | -0.209 | 0.004 | 0.504 + 0.864i |

**Crucially, the spike's wider T-range [9.5, 18.9] (#287) is NOT explained by
these brackets**: every bracket sits at T ≈ 9.3-10.3 (inside the Phase 2
walk's extent), not at T ≈ 18.9.

### Family-switch attempts

For each k ∈ {3, 4, 5, 6} the script picked the bracket with the smallest
crossing distance and attempted a family switch:

1. Recompute the parent member's monodromy eigenvalue + eigenvector at the
   primitive k-th-root multiplier.
2. Project the imaginary part of the complex eigenvector onto the
   symmetric-tulip free components `(z0, ydot0)` — the symmetric IC manifold
   y0 = xdot0 = zdot0 = 0 keeps the corrector residual valid.
3. Perturb the parent IC along that direction with steps ε ∈ {1e-3, 5e-3,
   1e-2, 3e-2} (escalating until the corrector escapes the parent's basin).
4. Re-correct via `correct_general_periodic_3d` with period guess
   `k × T_parent` and `FREE_VARS_SYMMETRIC_TULIP`.
5. Reject if the period collapses back to `T_parent` (slid along the parent
   family). Accept if it lands at a clearly different period.
6. Continue the accepted seed with `continue_general_3d_family` for up to 20
   members per direction to confirm it lives on a sub-family.

### Sub-family results — all 4 accepted

| Bracket k | Parent T | Switched T_TU | Period ratio | Sub-family members | Sub-family T extent |
|---|---|---|---|---|---|
| 3 (step +86)  | 9.545 | **20.287** | 2.125 | 22 | 20.28 – 20.39 |
| 4 (step +112) | 9.289 | **35.935** | 3.869 | 41 | 35.84 – 36.03 |
| 5 (step +55)  | 9.852 | **30.030** | 3.048 | 41 | 29.93 – 30.13 |
| 6 (step -66)  | 9.556 | **43.550** | 4.557 | 41 | 43.49 – 43.61 |

All closures excellent:
- corrector residuals 1e-12 to 1e-14
- independent Radau closure residuals 1e-10 to 1e-11
- well below the closure_tol = 1e-6 gate

### Important caveat: period ratios are not integer k

The switched periods are NOT `k × T_parent`:

- k=3 → 2.125 × T_parent  (expected 3.0)
- k=4 → 3.869 × T_parent  (expected 4.0)
- k=5 → 3.048 × T_parent  (expected 5.0)
- k=6 → 4.557 × T_parent  (expected 6.0)

Two possibilities:

1. The corrector's Newton step searched for the nearest periodic orbit
   compatible with the symmetric-tulip residual and the perturbation
   direction, and the nearest such orbit happens to be a **different**
   periodic orbit (possibly a multi-cover of a SIBLING family member or a
   genuinely new branch) at a non-k period ratio.
2. Each bracket is a Neimark-Sacker (torus-birth) bifurcation, NOT a clean
   period-multiplying one — the genuine branching object is a 2-torus
   (quasi-periodic), not a periodic orbit. The corrector landed on a
   nearby periodic-orbit family that happens to sit close to the parent in
   state space (the "phase-locked" periodic orbits inside the torus).

Either way, the sub-families have **narrow T extents** (e.g., 35.84–36.03 for
the k=4 branch — width ≈ 0.2 TU on a 36 TU period, or ≈ 0.5%), which is
consistent with locally-converged 1-parameter sub-families, not artefacts.

The k=4 sub-family at T ≈ 36 has a particularly interesting Floquet
structure: TWO distinct hyperbolic pairs (|λ| ≈ 306 and |λ| ≈ 90) — this is
a doubly-unstable orbit, the kind that lives near a 3D heteroclinic web.

### Phase 4 next step — recommendation

The **k=3 sub-family** (T ≈ 20 TU, 22 members, smallest ratio shift from
expected) is the cleanest candidate for V0-V5 gauntlet adaptation. Concrete
first IC:

    state_nd = (-0.77010, 0.0, -0.19553, 0.0, -0.36978, 0.0)
    T_TU     = 20.2866 TU  (≈ 88.16 days)
    jacobi   = (recompute from state_nd at EM mu = 0.012150584270572)
    system   = Earth-Moon CR3BP, l_km = 384400, t_s = 375699.8

BEFORE running it through the gauntlet, **Phase 4 must first rerun
`literature_check`** with a more specific signature that includes the
sub-family's period (T ≈ 88 days, well outside the Antoniadou-Voyatzis
single-revolution-period regime — though their corpus covers spatial resonant
periodic orbits broadly, so the structural fingerprint may still hit). If the
sub-family is genuinely a new structural identity (e.g., 3-fold symmetric
spatial tulip), the corpus must be widened OR a literature-fresh verdict
recorded for human review.

### Honest verdict

**The 3D-Braik-Ross (1,1) family work is a technique-proven rediscovery** of
Antoniadou-Voyatzis 2018 spatial CR3BP (the parent family, confirmed by
Part A).

Part B surfaced **four genuine sub-families** at non-integer period ratios.
These are real periodic-orbit families with excellent closure properties, but
they are NOT the period-doubled or k-period-multiplied branches the task
description anticipated:

- The Phase 2 fold points (step +11 T-fold and step -3 z0-fold) sit at
  T ≈ 10.28 — these are **fold-turning points** (tangent component sign
  flips) in the natural-parameter sense, NOT period-multiplying
  bifurcations.
- The actual period-multiplying signals are Neimark-Sacker brackets
  (complex pair on unit circle through primitive k-th roots), all classified
  uniformly as `neimark_sacker` — no real-multiplier crossings of ±1 except
  the trivial unit pair.

The **spike (#287) wider T range [9.5, 18.9] is NOT reproduced by any of
the bifurcation-tracking sub-families** (which sit at T ≈ 20, 30, 36, 44 TU
— all much larger than 18.9). The most likely explanation: the spike's
wider T range was either (a) a different family fragment the spike found by
coincidence, (b) a search-method artefact (numerical noise placing nearby
periodic-orbit candidates near the same parent family), or (c) genuinely a
DIFFERENT family the spike sampled that the Phase 2 tracer never reached
from the central seed. Investigation deferred.

### Test status

`uv run pytest tests/search/test_literature_check.py
tests/search/test_bifurcation_detector.py
tests/search/test_cr3bp_3d_family_tracer.py -x --timeout=120` → **36 passed
in 91.6s**. No edits to `literature_check.py` or `bifurcation_detector.py`
required; the existing schema covers the Earth-Moon CR3BP 3D fingerprint
adequately.

### Outputs

- `data/lit_check_299_3d_family.jsonl` — Part A per-member literature-check
  verdicts (265 rows + 1 header)
- `data/family_296_3d_subfamilies_299.jsonl` — Part B accepted sub-family
  members (4 bracket records, each with 22–41 sub-family members, + 1 header)
- `scripts/run_299_lit_check_3d_family.py` — Part A driver
- `scripts/run_299_bifurcation_track_3d_family.py` — Part B driver

No catalogue writeback. Both JSONLs are working data; promotion to catalogue
rows must wait for V0-V5 gauntlet adaptation for 3D orbits (Phase 4+ /
future task).
