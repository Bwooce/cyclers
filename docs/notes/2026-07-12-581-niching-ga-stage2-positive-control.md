# #581 stages 1-2: niching-GA positive control — Gurfil-Kasdin (2002) reproduction

Date: 2026-07-12 AET. Scope: #581 stages 1 and 2 ONLY (stage 3, aiming the
layer at novel targets, remains open and requires a fresh human/Opus/Fable
review of this note before dispatch — per the task's own gate).

Source: Gurfil & Kasdin (2002), "Niching genetic algorithms-based
characterization of geocentric orbits in the 3D elliptic restricted three-body
problem", CMAME 191, 5683-5706. Digest:
`docs/notes/2026-07-12-digest-gurfil-kasdin-2002-er3bp-geocentric-orbits.md`.

## Verdict up front

**11/14 families recognizably reproduced** (A, B, C, D, E, G, H, J, K, L, M) in
one niching-GA run per Table 2 optimization set, at the paper's own GA
constants, matched against the paper's own published Table 3/4 numbers by a
pre-registered, non-bit-exact criterion (§4). Of the 3 misses (F, I, N),
diagnosis (§4.2) shows 2 (F, N) are genuine local-optimum competition on a
near-tied fitness landscape — the GA found a DIFFERENT, equal-or-higher-fitness
point in the same family's search box, not a failure to search — and 1 (I) is
a genuine under-convergence in a fragmented, higher-dimensional (4 free
variables) box. This is a strong positive control: the core mechanism
(deterministic-crowding niching over the geocentric ER3BP) works and
reproduces the paper's central claim (multiple co-existing families from one
run per box), with an honestly-characterized, not-tuned-away failure mode.

## 1. Stage 1 build (commit `194da55`)

Three deliverables, all unit-tested:

1. **`core/er3bp_geocentric.py`** — the paper's Earth-centered pulsating ER3BP
   frame, implemented as a geocentric offset (`x_bary = x_geo + (1-mu)`) over
   the existing Nechvile/true-anomaly machinery in `core/er3bp.py` (NOT a new
   frame from scratch). The offset identity was verified algebraically
   (substituting the translation into the barycentric EOM reproduces the
   paper's Eq. 9-11 term-for-term, including the `+(1-mu)` indirect term) and
   numerically (RHS parity to 6.7e-16 against an independently-transcribed
   Eq. 9-11 in the test). Includes the paper's Eq. 15 objective
   (`1/((rmax-rmin)^2+1)`) with the Eq. 17 collision constraint
   (`rmin > R_Earth`) and an escape guard as death penalties.
2. **`search/niching_ga.py`** — Deterministic Crowding exactly per the paper's
   p. 5687 pseudocode: random parent pairing, two children per pair via
   single-point crossover on a concatenated binary string + per-bit mutation,
   each child competing against the genotypically CLOSER parent (Euclidean in
   bounds-normalized phenotype space), ties to the child. Table 1 constants
   as defaults (pop 200, 400 generations, 32 bits/variable, p_cross 0.999,
   p_mut 0.001). Checkpoint/resume + process-pool evaluation + genotype
   dedupe. This is a genuine multi-optima-per-run mechanism — scipy's
   `differential_evolution` (all 6 existing evolutionary lanes in the repo)
   exposes no replacement hook and collapses to one optimum per run.
3. **Tests** — `tests/core/test_er3bp_geocentric.py` (sourced Table 3/4
   goldens + independent-RHS parity + fitness sanity),
   `tests/search/test_niching_ga.py` (two-peak niche-survival property,
   checkpoint-resume determinism, clone-children replacement invariance,
   encoding bounds).

## 2. Load-bearing source-interpretation findings

These were REQUIRED to make the reproduction work and are documented here
because any future consumer of the paper's tables will hit them:

- **Table 2/3 vectors are printed in INTERLEAVED order `[x, x', y, y', z, z']`**,
  despite the paper's own Eq. 13 defining the state as `[x,y,z,x',y',z']`.
  Proven by cross-checking every Table 3 IC against Table 4's r0/v0 values:
  under interleaved ordering family J's r0 = 7,527,807 km, I's = 3,177,032 km,
  K's = 883,956 km, L's = 1,638,694 km, N's = 1,521,047 km all match Table 4
  EXACTLY; under the Eq. 13 ordering they are off by 10-50%. Family A's
  documented Henon condition `y'0 = -1.9967 x0` (p. 5691) also only holds
  under interleaved ordering.
- **Table 4's km distances are normalized distance x 1 AU** (1.496e8 km, the
  paper's own "slightly abused" unit convention, p. 5691), not the pulsating
  dimensional distance. v0 values use the proper dimensional conversion
  `v * R(theta0) * omega(theta0)`.
- **Table 2's `/sqrt(2)`, `/sqrt(3)` factors apply to the whole printed
  vector** (min and max), keeping the initial geocentric distance in
  [~1e6, ~1e7] km when 2 or 3 position components are free.
- **Family F's Table 4 rmax (1,002,197 km) appears to be a typesetting
  duplication** of the D/E row value: it contradicts F's own r0 = 1,000,000 km
  (an orbit's max distance cannot sit 0.2% above its start while its min is
  225,000 km below), and our high-accuracy re-integration of F's Table 3 IC
  gives rmax = 1,138,177 km with rmin (224,900 -> 225,537 km, 0.3%) and r0
  (1,000,000 km, exact) both matching. All 13 other families reproduce all
  three of rmin/rmax/r0 to ~0.03% (see §3). Per
  [[feedback_respectful_errata_framing]]: typesetting slips happen; flagged,
  not asserted.
- **Family H's Table 4 v0 = 0 is a display rounding**: its Table 3 IC has
  x' = -0.00097613 (~0.03 km/s). The ERO taxonomy ("characterized by zero
  initial velocity") tolerates this; our type classifier uses a 2e-3
  normalized threshold accordingly.

## 3. Dynamics validation (deterministic, pre-GA)

Before any GA run, all 14 published Table 3 ICs were propagated for 1 year
(DOP853, rtol 1e-11, max_step 0.01) in the new frame and compared to Table 4:

- r0: exact (to the printed digit) for all 14 families.
- v0: matches to +-0.001 km/s for all 14 (H prints 0 for 0.030 km/s, see §2).
- rmin/rmax: within ~0.03% for all families (worst: E rmin 0.06%), except
  family F's rmax (§2, suspected source-table slip).

Free cross-check: the paper's mu = 3.0034495182e-6 agrees with this project's
own constants (`PLANETS["E"].mu_km3_s2 / (mu_E + MU_SUN_KM3_S2)` =
3.0034805950e-6) to 1.0e-5 relative — consistent with an Earth-only (Moon
excluded) GM under a slightly different constants set.

This is the strongest possible foundation for the stage-2 GA control: the
dynamics, frame, unit conventions, and IC decoding are exactly right; any
reproduction failure would be attributable to the GA layer alone.

The golden test (families A, F, J at 0.5% tolerance) is in the permanent
suite: `tests/core/test_er3bp_geocentric.py::test_golden_family_features_reproduce_table4`.

## 4. Stage 2: the 12-set reproduction experiment

**Design — faithful to the paper's actual experimental protocol:** 12
independent optimization runs, one per Table 2 optimization set, each with
the paper's own bounds, theta0, and GA constants (pop 200, 400 generations,
32-bit binary encoding, p_cross 0.999, p_mut 0.001), seeds 581001-581012.
GA fitness evaluations integrate 1 year at rtol 1e-9 with 2000-sample
distance tracking; terminal collision/escape events; death-penalty fitness 0.
Driver: `scripts/run_581_gurfil_reproduction.py`; checkpoints, per-generation
runlogs, final populations and the analysis summary under
`data/found/581_niching_ga/`.

**Match criterion ("recognizably reproduced", not bit-exact):** for each
expected family, the final-population member NEAREST the published Table 3 IC
(bounds-normalized RMS gene distance) must (a) lie within 0.10 normalized RMS
of it — i.e. the population maintains a niche at the published IC, (b)
re-integrate (rtol 1e-11, 5-year extension) to the paper's own type under the
paper's taxonomy (planar/3D x DRO/DPO/ERO + practical-stability -> DEO
override), and (c) have 1-year rmin AND rmax each within a factor of 2 of
Table 4's values. All three must hold; actual ratios are reported.

### 4.1 Results table

| Set | Family | Match | IC RMS dist | Type (got vs expected) | rmin ratio | rmax ratio |
|---|---|---|---|---|---|---|
| 1 | A | **MATCH** | 0.0048 | DRO / DRO | 0.967 | 0.938 |
| 1 | B | **MATCH** | 0.0036 | DRO / DRO | 0.962 | 0.969 |
| 1 | C | **MATCH** | 0.0002 | DRO / DRO | 1.002 | 1.002 |
| 2 | D | **MATCH** | 0.0002 | DPO / DPO | 1.001 | 0.998 |
| 3 | E | **MATCH** | 0.0002 | DPO / DPO | 1.000 | 0.998 |
| 4 | F | miss | 0.0132 | DRO / DRO | 4.964 | 1.173 |
| 5 | G | **MATCH** | 0.0001 | ERO / ERO | 1.012 | 1.001 |
| 6 | H | **MATCH** | 0.0040 | ERO / ERO | 0.765 | 0.952 |
| 7 | I | miss | 0.0183 | DEO / DEO | 0.216 | 1.109 |
| 8 | J | **MATCH** | 0.0053 | 3D DRO / 3D DRO | 1.007 | 1.008 |
| 9 | K | **MATCH** | 0.0785 | 3D DRO / 3D DRO | 0.854 | 1.307 |
| 10 | L | **MATCH** | 0.0146 | 3D DEO / 3D DEO | 0.863 | 1.032 |
| 11 | M | **MATCH** | 0.0002 | 3D ERO / 3D ERO | 1.003 | 1.001 |
| 12 | N | miss | 0.0491 | 3D ERO / 3D DEO | 0.150 | 0.775 |

**11/14.** Raw data: `data/found/581_niching_ga/analysis_summary.json`,
`set{01..12}_final.npz`, `set{01..12}_runlog.jsonl`.

### 4.2 Diagnosing the 3 misses

For each miss, the published IC's own fitness (under our Eq. 15 objective)
was compared against the fitness of the nearest population member the GA
actually settled near, to distinguish "the GA didn't find this basin" from
"the GA found a different, comparably-good point in the same basin/family":

- **F** (set 4, DRO): published F's own fitness is **0.999963**. The GA's
  best nearby cluster (`x0=0.0076, y'0=-0.029`, 6 members) scores
  **0.999998** — HIGHER. Family F sits on the same continuous DRO ridge as
  A/B/C (`y'0 ~ -2 x0`, per Henon's condition, cited by the paper itself,
  p. 5691); the ridge point the GA converged to at generation 400 is a
  genuinely better optimum of Eq. 15 within the same search box, at slightly
  larger rmin (F's own distinguishing feature — the paper calls it "the
  closest-approach DRO" at 224,900 km — happens to be a LOCAL, not global,
  optimum of the flat objective). This is basin/ridge competition, not a
  missed family: the DRO family (topologically) was found; the specific
  low-periapsis representative member was not the one that survived.
- **N** (set 12, 3D DEO): published N's own fitness is **0.999945**. The
  dominant final-population cluster (64/200 members, `x0=z0=0.0047`) scores
  **0.999959** — again slightly higher, and this cluster is recognizably the
  SAME basin as family M (set 11's own near-Earth 3D ERO, rmin ~17,000-35,000
  km) rather than N's own (rmin 409,758 km). Sets 11 and 12 have structurally
  overlapping search boxes (M: x,y,z free /sqrt(3); N: x,z free /sqrt(2),
  same magnitude range) and Eq. 15's flat, range-only objective has no
  preference between a very-close-approach basin and N's own — the
  higher-fitness M-like basin displaced N's own optimum from the final
  population. Same diagnosis as F: basin competition on a near-tied
  landscape, not a search failure.
- **I** (set 7, DEO, 4 free variables — the highest-dimensional 2D set):
  published I's own fitness is **0.999881**. The population fragmented into
  **95 distinct clusters** (out of 200 individuals) at generation 400, none
  dominant, and the nearest-to-I cluster scores only **0.999633** — LOWER
  than I's own fitness. This is a genuine under-convergence: at the paper's
  own GA constants (pop 200, 400 generations), a 4-free-variable (128-bit)
  box did not consolidate onto any single strong optimum in the time given,
  consistent with the combinatorics of a larger search space and the same
  budget. This is the one miss attributable to the search not having done
  enough work, not to basin competition.

**Net read:** the mechanism (niching, frame, objective, IC decoding) is
sound; 2 of 3 misses are the DC algorithm legitimately settling on a
different, equal-or-better local optimum than the paper's chosen
representative (arguably correct behavior, not a bug), and 1 miss is a
plausible under-convergence in the paper's own highest-dimensional set at its
own stated generation budget. No parameter was tuned after seeing these
results; this is the as-run outcome at the paper's Table 1 constants.

## 5. Honest caveats

- The DRO/DPO "families" of set 1 are niches along a CONTINUUM (the Henon
  `y'0 ~ -2 x0` ridge); the deterministic-crowding population retains the
  whole ridge (final set-1 population spans x0 = 0.0067 to 0.0667 with 34
  members at family A's scale, all with y'0/x0 in [-2.00, -1.98]). The
  paper's A/B/C are three representative points on that ridge distinguished
  by xz-crossing counts; our match criterion (niche member near the published
  IC + type + features) is the honest formalization of "found the same
  family", but it is a judgment criterion, not a theorem.
- The GA is stochastic; one seed per set was run (no cherry-picking over
  seeds). The paper likewise reports single runs per set.
- Fitness near the optimum is extremely flat (all bounded compact orbits
  score ~0.99+), so niche survival — not fitness resolution — is what the
  mechanism is being credited for. That is exactly the paper's claim too.

## 6. Stage-3 assessment (NOT executed — gate honored)

This section is a REPORTED ASSESSMENT ONLY. Per #581's explicit staging gate,
stage 3 (aiming the niching layer at novel targets) has NOT been dispatched
and requires an independent fresh human/Opus/Fable review of this note first.

Personal read for the record: 11/14 with 2 of the 3 misses attributable to
legitimate basin competition (not search failure) and the third to a
plausible, explainable under-convergence at the paper's own stated budget is
a strong positive control by the honest read of
[[feedback_verify_gauntlet_with_positive_control]] — the mechanism does what
Gurfil-Kasdin claim it does: recover multiple co-existing, topologically
distinct orbit families from single runs, including the ERO/DEO cases (G, H,
M) that a resonance-lattice or grid-based search would not generate at all.
That is the specific capability #581's origin motivation named as missing
from the codebase's existing single-optimum DE lanes. On this evidence, stage
3 looks worth authorizing — but that determination is explicitly reserved for
the required independent review, not this agent.
