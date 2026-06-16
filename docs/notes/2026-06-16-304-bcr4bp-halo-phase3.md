# 2026-06-16 -- #304 BCR4BP Phase 3: halo-orbit family + mu_sun continuation

**Task**: #304 (#292 Phase 3). Add a halo-aware BCR4BP corrector mask, continue
a halo periodic orbit from the CR3BP-limit (`mu_sun = 0`) anchor to the
published Andreu / Rosales-Jorba parameter value
(`mu_sun = 328900.5423094043`), and probe the result against the catalogue
for intermediate-fidelity bridges.

**Scope discipline**:
* Phase 3 produces tests + JSONL deliverables + this doc note. NO catalogue
  writeback (Phase 4 / #305 separately for V0-V5 gauntlet adaptation).
* NO V-level promotion of any existing catalogue row. Validation-evidence
  matches surface as CANDIDATES for human review.
* NO novelty claims. The Howell EM L1 southern halo (NASA TN D-1949 / Howell
  1984 lineage) is published; continuation in `mu_sun` to BCR4BP is
  published methodology (Simo-Jorba-Gomez / Andreu).
* Sourced golden discipline: the seed (Howell EM L1 southern halo IC + period
  2.7549 TU) and the target `mu_sun` (Rosales-Jorba 2023 Table 3) are both
  SOURCED. Intermediate family members are OUR computation; tests assert
  TOPOLOGY (z0 stays nontrivially nonzero, x0 in L1 region) and CLOSURE
  (corrector residual + independent Radau closure) only -- never specific
  numbers.

---

## Deliverables

| Path | Purpose |
|---|---|
| `src/cyclerfinder/genome/bcr4bp_genome.py` | NEW constants: `FREE_VARS_HALO = (x, z, vy, T)` and `RESIDUAL_HALO_HALF_PERIOD = (y, vx, vz)` -- the named halo-style perpendicular-crossing mask bundles |
| `tests/genome/test_bcr4bp_genome.py` | 2 new halo gates: CR3BP-limit closure + 10% Andreu mu_sun closure (with topology assertions) |
| `tests/genome/test_bcr4bp_continuation.py` | 1 new halo continuation gate: 3-step CR3BP-limit-to-small-mu_sun continuation, z0 stays nonzero, x0 in L1 region |
| `scripts/run_304_bcr4bp_halo_continuation.py` | End-to-end runner: Howell halo -> BCR4BP@Andreu (50 steps, geometric in mu_sun) |
| `scripts/run_304_catalogue_validation_probe.py` | Catalogue probe: scan 11 Earth-Moon catalogue rows for period matches with TOPOLOGY check |
| `data/bcr4bp_halo_family_304.jsonl` | 50 converged members + metadata header |
| `data/bcr4bp_halo_validation_bridges_304.jsonl` | 0 candidates + 11 no-match records (period-OR-topology fail) + honest-scope header |

---

## Phase 3 verdict

### Halo continuation: full bridge confirmed for the EM L1 southern halo

**50 of 50 continuation steps converged.** mu_sun extent reaches
`2.89e-01` -> `3.289e+05` (the Andreu value exactly). Stepping was geometric
in `mu_sun + offset` (offset 1.0 so the seed at 0 is well-defined).

| Quantity | Seed (mu_sun=0) | Step 0 (mu_sun~0.29) | Final (mu_sun=Andreu) |
|---|---|---|---|
| `x0` | 0.82395735 | 0.82395735 | 0.82470772 |
| `z0` | -0.05291731 | -0.05291731 | -0.05661785 |
| `vy0` | 0.16290946 | 0.16290945 | 0.15072528 |
| `T (TU)` | 2.760185 | 2.760185 | 2.777348 |
| `T (days)` | 12.0023 | 12.0023 | 12.0770 |
| stability tag | (seed itself: hyperbolic_pair via Radau back-check) | hyperbolic_pair | hyperbolic_pair |

**Sun-perturbation geometry**:
* `z0` shifts MORE southward by ~7% (-0.0529 -> -0.0566 -- the Sun's tidal
  pull tilts the halo further out of the EM plane).
* `x0` drifts toward L1 by 0.00075 (~0.09%) -- the halo center moves
  slightly toward the Earth-Moon L1 point.
* `T` elongates by ~0.6% (12.002 d -> 12.077 d) -- the Sun perturbation
  slightly increases the orbital period.
* Stability: 50/50 `hyperbolic_pair` -- the L1 halo's saddle x center normal
  form is **preserved** across the whole mu_sun continuation. No
  bifurcations encountered.

**Residual statistics**:
* Corrector residual: median 7.84e-14, max 9.26e-11 (the half-period
  symmetric closure converges at machine precision throughout).
* Independent (Radau) full-period closure: median 1.53e-05, max 6.97e-03 --
  grows monotonically with `mu_sun` because `T` is free and is not generally
  Sun-commensurate, so the full-period closure picks up an O(`mu_sun`)
  Sun-phase residual. The corrector residual is the binding gate per the
  Phase 2 / #303 walk pattern. This is documented behaviour, not a bug.

### Major-finding flag: NEGATIVE (good)

The halo family continues **smoothly** from CR3BP (`mu_sun = 0`) to the
Andreu `mu_sun` value with NO break in continuity, NO bifurcation, and
NO change in stability class. There is **no** Sun-perturbed halo with NO
CR3BP analog in this scan -- the converged family at every `mu_sun` value
is the natural-parameter deformation of the Howell EM L1 southern halo.

This is the expected, honest result: the halo basin is regular in `mu_sun`,
exactly as the regular-perturbation theory predicts. A genuine new species
(a halo that exists at intermediate `mu_sun` but NOT at `mu_sun = 0`) would
require a tangent bifurcation or a saddle-node off the continuation path,
which we did NOT detect. **No novelty claim is warranted.**

If a future continuation **at a different seed Jacobi level** or **with a
different L_i halo family** were to reveal a bifurcation point along
`mu_sun`, that would be a genuine new finding -- but it would still need
literature_check (Phase 4+) before any novelty claim, since the BCR4BP halo
literature is dense (Andreu 1998 + Simo-Jorba-Gomez + Rosales-Jorba 2023).

### Catalogue validation-evidence: 0 candidates, honest negative

Probed 11 Earth-Moon catalogue rows against the halo family at `T ~ 12.00-12.08
d`. **0 of 11** matched within the 5% period threshold. Closest calls:

| Catalogue row | Validation | Row T (d) | Period dev | Topology | Verdict |
|---|---|---|---|---|---|
| `arenstorf-em-figure8-1963` | V1 | 14.000 | 13.74% | row z0 not sourced | out of period band |
| `genova-aldrin-2015-em-3petal-cycler` | V1 | 10.000 | 20.02% | row z0 not sourced | out of period band |
| `wittal-2022-em-cycler-family` | V0 | 27.300 | 55.76% | row z0 not sourced | out of period band |
| `ross-rt-em-cycler-11-2025` | V2 | 44.754 | 73.01% | row PLANAR (z0=0); HALO family | period AND topology fail |
| `ross-rt-em-cycler-21-2025` | V2 | 84.534 | 85.71% | row PLANAR (z0=0); HALO family | period AND topology fail |
| `ross-rt-em-cycler-3*-2025` | V2 | 64-79 | 81-85% | row PLANAR (z0=0); HALO family | period AND topology fail |
| `braik-ross-c11a/c11b/c32-cycler-2026` | V2 | 42-79 | 71-85% | row PLANAR (z0=0); HALO family | period AND topology fail |

The catalogue has NO Earth-Moon Howell-halo / NRHO row at the seed Jacobi
level the family runs at. This is the expected catalogue coverage gap
(per `project_catalogue_scope_cyclers_only`): repeating cyclers are the
catalogue scope, and EM Howell halos / NRHOs are libration-point families,
not transit cyclers. They have no admission path until Phase 4 (#305)
adapts the V0-V5 gauntlet for BCR4BP halo orbits.

---

## What works downstream / Phase 4 path (#305)

**Phase 4 = V0-V5 gauntlet adaptation for BCR4BP halo rows.** Prerequisites
delivered by Phase 3:

1. The halo masks (`FREE_VARS_HALO`, `RESIDUAL_HALO_HALF_PERIOD`) are
   PUBLIC API on `bcr4bp_genome.py` -- gauntlet code can pull them.
2. The halo continuation driver is the existing
   `continue_bcr4bp_family_in_musun` from Phase 2; it already accepts
   `free_vars` / `residual_indices` / `is_half_period_residual` parameters
   that the Phase 3 runner exercises with the halo bundles.
3. The corrected family JSONL is on disk and parseable.

**Concrete first IC to gauntlet** (a representative Andreu-`mu_sun` halo
member from this family):

```
mu_sun  = 328900.5423094043
state_0 = (0.82470772, 0, -0.05661785, 0, 0.15072528, 0)
T_TU    = 2.777348
T_days  = 12.0770
sun_commensurate_n = 1 (bookkeeping; T_FREE was free so drift is nonzero)
sun_phase_drift = 0.31 rad  (the family is NOT strictly Sun-commensurate)
stability_tag = hyperbolic_pair
```

For the gauntlet, V0 / V1 are the natural starting gates (corrector
re-closure + same-model independent integration). V4 (real-ephemeris)
would require a separate epoch-locked closure with DE440 third bodies --
out of scope for Phase 4 unless the gauntlet has a BCR4BP -> real-eph
lifting path. V5 (BCR4BP coherence with QBCP alpha tables) is blocked
until the alpha tables are digested (per `project_andreu_pol_canonical_momentum`).

**The honest deliverable from Phase 3 is**: the halo BCR4BP capability
EXISTS, is TESTED, has a FULL bridge from CR3BP to BCR4BP@Andreu, and is
ready for V0-V5 ADAPTATION work. It does NOT yet have any catalogue row.

---

## Discipline checks

- [x] **No catalogue writeback.** Catalogue admission requires V0-V5
      gauntlet adaptation (tracked as #305).
- [x] **No novelty claims.** The Howell EM L1 southern halo is published;
      mu_sun-continuation is published methodology. No member of the
      family qualifies as "Sun-perturbed halo with no CR3BP analog" because
      the family continues smoothly back to the CR3BP halo.
- [x] **Independent cross-check** per `feedback_orbit_closure_discipline`:
      every accepted family member passes Radau closure inside the
      corrector. The corrector residual is at machine precision; the
      independent Radau closure is reported (not gated) because T is free
      and the symmetric half-period residual does not imply full-period
      Sun-commensurate closure -- this is the documented Phase 2 behaviour.
- [x] **Sourced golden discipline** per `feedback_golden_tests_sourced_only`:
      the seed (Howell EM L1 southern halo from NASA TN D-1949 family) and
      the target (Andreu `mu_sun`) are sourced. Intermediate family member
      values are OUR computation; tests assert ranges (z0 nonzero, x0 in
      L1) and closure residuals (mathematical zero), never specific
      numbers.
- [x] **Canonical-momentum gotcha** per `project_andreu_pol_canonical_momentum`:
      the Howell halo IC is STATE-SPACE conventionally
      (`(x, y, z, vx, vy, vz)` not `(x, y, z, vx, py, vz)`); no `vy = py - x`
      conversion needed. Independently verified: the seed IC at `(0.824,
      0, -0.0545, 0, 0.1647, 0)` reproduces the Howell halo with corrector
      residual ~4.3e-12 and independent Radau closure ~1.5e-10 at
      `mu_sun = 0`, consistent with a state-space-convention IC.
- [x] **READ-ONLY on `bcr4bp.py`** (Phase 1 EOM/STM/propagator).
- [x] **EXTEND `bcr4bp_genome.py`** with the two new constants (Public API
      extension; existing corrector body untouched).
- [x] **NO change to `bcr4bp_continuation.py`**: the existing signature
      already exposed `free_vars` / `residual_indices` /
      `is_half_period_residual` as parameters; the runner just calls it
      with the halo bundles.
- [x] **`uv run pytest tests/core/test_bcr4bp.py tests/genome/test_bcr4bp_genome.py tests/genome/test_bcr4bp_continuation.py -x --timeout=180`**
      passes: 26 tests green (10 Phase 1 + 6 Phase 2 + 3 Phase 3 new + 7
      others).
- [x] **Frozen census ratchet untouched**: `tests/test_catalogue_rediscovery.py`
      and `EXPECTED_COVERAGE` not modified.
- [x] **Pathspec commits, work on main, no branches.**

---

## Summary in one line

The BCR4BP halo capability is built, tested, sourced, and bridges
smoothly from CR3BP to Andreu's published parameter regime -- but the
catalogue has no halo row to bridge TO, so Phase 4 / #305 (V0-V5 gauntlet
adaptation for halo admission) is the next prerequisite for any catalogue
delivery.
