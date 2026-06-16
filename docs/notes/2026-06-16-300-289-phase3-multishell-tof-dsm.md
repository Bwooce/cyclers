# Phase 3 of #289 — multi-shell BFS, TOF optimisation, DSM extension (#300)

Phase 3 of task #289 (Track-A epoch-locked trajectory substrate) closes the
three coupled deficiencies the Phase-2 doc named at the Galileo VEEGA
calibration probe:

  1. **Multi-shell BFS** — `find_mga_chains(..., multi_shell=True)` admits
     per-flyby V_inf shifts within the Strange-Longuski 2002 §12 pump
     envelope.
  2. **Per-leg TOF optimisation** — `optimise_chain_tofs` runs Nelder-Mead
     over (launch_epoch, leg_tofs) with loss
     `closure_residual_kms + 2*flyby_continuity_max_dv_kms`.
  3. **DSM extension** — optional `dsm_specs` field on
     `EpochLockedTrajectory`. The closure driver splits each DSM-equipped
     leg at the DSM fraction, propagates ballistically to the split,
     applies the user-prescribed ΔV, propagates the post-impulse state
     forward to the leg's arrival epoch, and computes the V_inf at the
     arrival body from the actual (post-DSM) velocity.

The independent Kepler cross-check is preserved per
`feedback_orbit_closure_discipline`; the ballistic case (`dsm_specs=None`
or empty) is bit-identical to the Phase-1 closure.

## What sourced anchors gate which test

### Multi-shell BFS

`delta_vinf_max_kms(body, vinf_in_kms)` is a closed-form expression for the
maximum per-flyby V_inf change at the next body, derived from the
Strange-Longuski 2002 §12 leveraging identity
`V_inf * (1 - cos δ_max)` where δ_max follows from `max_bend`. The sourced
calibration is from the corpus rev `568d8a4` mining note for
Strange-Longuski 2002 §12 (the value cited in
`docs/notes/2026-06-16-298-289-phase2-tisserand-enumerator.md`): Earth at
V_inf = 9 km/s admits a per-flyby V_inf change of ~2.7 km/s. Our
geometric formula returns 3.24 km/s — a factor of 1.20 above the sourced
2.7 figure. This is the GEOMETRIC bound (the rotation may not always
reach the maximum cone half-angle because the rotated v_inf must point at
the next body, not anywhere on the V_inf globe). The realised-vs-geometric
factor is exposed via `pump_envelope_factor` on `find_mga_chains` for
callers that want a tighter envelope; the default 1.0 uses the geometric
maximum and is the conservative envelope (never UNDER-bounds the
published number — that would cause silent search misses).

### TOF optimisation

The test seeds a deliberately-off Tito Earth-Mars-Earth chain (TOFs off
by ~25 days each, epoch off by 4 days) and asserts that the optimised
loss strictly drops. On DE440 the optimisation reduces the loss from
~0.5 km/s to ~0.05 km/s — a 10× reduction — recovering the published
Tito V_inf tuple within the 1.5% sourced-side residual.

### DSM extension

Tested with a 10 m/s DSM at the midpoint of Tito's outbound leg
(E → M). The DSM is wired correctly when:

* `dsm_delta_v_kms_per_leg` records the user-supplied ΔV magnitude;
* `per_encounter_vinf_kms` for the MARS flyby changes from the no-DSM
  value (the leg-0 arrival velocity is now the post-DSM-propagated
  velocity, not the ballistic-Lambert v2);
* the independent-cross-check residual changes (the Kepler propagator
  re-traces the split-leg geometry and sees the DSM).

A 10 m/s impulse produces ~3 m/s shift in Mars V_inf — sane proportional
behaviour. The test enforces a 1 km/s ceiling on the Mars V_inf shift.

## Galileo VEEGA Phase 3 verdict

`data/scan_300_galileo_veega_multishell.jsonl` records the (V,E,E,J)
chain reproduction with all three Phase-3 capabilities:

| Stage | Closure residual (km/s) | Worst flyby ΔV (km/s) | Comment |
|---|---|---|---|
| Multi-shell BFS seed (ballistic) | 25.5 | 71.7 | Geometric proposal, no optimisation |
| + TOF / epoch optimisation | 11.5 | 33.4 | Nelder-Mead, ±15 d epoch, ±40% TOF |
| + Hand-placed E2→J DSM (62 m/s) | 11.5 | 33.4 | Mid-leg-2; minor effect |

The TOF optimisation roughly halves the residual; the hand-placed DSM is
essentially ineffective at the Galileo scale because the residual gap is
structurally not on the E2→J leg alone — it's on the Lambert closure of
the whole sequence at V_inf bins that don't match Galileo's actual
geometry.

Honest honest scope (per `feedback_orbit_closure_discipline`):

* The multi-shell BFS surfaces (V,E,E,J) chains with the right Venus / E1
  / E2 V_inf SHELL pattern (~5 → 9 → 9 km/s) — the sourced Galileo
  pattern is admitted by the predicate. The Jupiter V_inf is pinned at
  the predicate threshold (`linkable_3d(E,J)` opens at V_inf ≥ 11 km/s in
  the inclined Tisserand model), structurally NOT Galileo's actual 5.6
  km/s. This is the predicate's leg-by-leg ballistic limitation, not a
  framework bug.
* Galileo's actual mission LAUNCHES FROM EARTH on 1989-10-18, swings BY
  VENUS in Feb 1990. The Phase-3 `(V,E,E,J)` sequence starts at the Venus
  flyby (Feb 1990) — the Earth liftoff and the Earth→Venus injection leg
  are NOT part of the chain. The Lambert closure-residual numbers
  measure the (V→E1→E2→J) interior chain only.
* The 11.5 km/s residual is dominated by Lambert geometry mismatch at
  the V_inf=11 Jupiter pinning. With a single Tisserand-graph hop from
  Earth-2 (V_inf=9) to Jupiter (V_inf=11), the heliocentric energy
  step is too large for a 997-day Hohmann-class arc; the real Galileo
  trajectory takes ~3 years on this leg and uses a perihelion ΔV at
  ~62 m/s with multiple intermediate state corrections that a single
  hand-placed DSM cannot represent.
* **5 km/s gate not reached**; the residual sits in [10, 20] km/s band.

## Phase 4 path

The structural gap is now sharply scoped:

1. **Automated DSM placement.** The hand-placed DSM in #300 fires at a
   chosen fraction with a chosen ΔV; the spacecraft pays the ΔV but the
   DSM doesn't redirect the trajectory toward the next body unless the
   ΔV is chosen specifically. Phase 4 should add DSM-optimisation: at
   each leg, optimise (fraction, ΔV vector) to minimise the closure
   residual. This is the Vasile-Conway 2006 MGA-DSM optimiser around the
   surviving multi-shell seeds.

2. **Precursor-MGA insertion matcher.** Resolve `inserts_into` against
   the live catalogue. The Aldrin S1L1 (Earth-Mars cycler) is the
   high-value first target — its outbound seed conditions
   (V_inf = 4.99 km/s at Earth, M = 5.10 km/s at Mars; see
   `docs/notes/2026-06-08-continuation-chain-s1l1-results.md`) admit a
   Phase-3-style `precursor_mga` rendezvous trajectory:

   ```yaml
   sequence: ("E", "S/C", "S/C", "M")  # depart Earth, two boost burns, intercept the cycler
   vinf_kms_at_encounters: (~3.0, ?, ?, ~5.10)  # the cycler-side V_inf
   inserts_into: "russell-ross-2025-s1l1-stable-A"
   orbit_class: "precursor_mga"
   ```

   The `match_precursor_to_cyclers(candidate, catalogue)` function should
   gate the Phase-3 (V_inf, body, epoch) tuple against the cycler's seed
   conditions.

3. **Multi-rev Lambert.** A single-rev prograde Lambert is the Phase-1
   default; long Earth → Jupiter legs admit multi-revolution geometries
   that the seed enumerator cannot select today.

4. **Eccentric-Earth Tisserand graph.** The `linkable_3d` predicate uses
   the coplanar / inclined Tisserand identity with circular Earth/Jupiter
   orbits. Real-Earth eccentricity (~0.0167) shifts the V_inf-shell
   boundaries; for the high-V_inf E→J case this matters at the 1 km/s
   level.

## Files changed

* `src/cyclerfinder/search/tisserand_mga_window.py` — multi-shell adjacency,
  envelope formula, TOF optimiser.
* `src/cyclerfinder/genome/epoch_aware_genome.py` — `DSMSpec`,
  `dsm_specs` field, DSM-aware closure driver.
* `tests/search/test_tisserand_mga_window.py` — envelope sourced golden,
  multi-shell tests, TOF-optimiser tests.
* `tests/genome/test_epoch_aware_genome.py` — DSM invariants, Tito + tiny
  DSM smoke test, ballistic-still-empty regression test.
* `scripts/scan_300_galileo_veega_multishell.py` — Galileo VEEGA Phase-3
  verdict driver.
* `data/scan_300_galileo_veega_multishell.jsonl` — three-stage scan
  output.

## References

* Strange, N. J. & Longuski, J. M., "Graphical Method for Gravity-Assist
  Trajectory Design", *J. Spacecraft and Rockets*, 39(1):9-16, 2002. The
  §12 pump envelope formula is the sourced anchor for
  `delta_vinf_max_kms`. KNOWN_CORPUS at corpus rev `568d8a4`.
* Petropoulos, A. E. & Longuski, J. M., "Shape-Based Algorithm for the
  Automated Design of Low-Thrust, Gravity Assist Trajectories",
  AAS/AIAA Astrodynamics Specialist Conf., 2000. Pump-tour
  combinatorics.
* Diehl, R. E., Belbruno, E. A., Roberts, R. M., "Galileo Mission Plan
  for the (V, E, E, J) Trajectory", AAS 1986. KNOWN_CORPUS at corpus rev
  `568d8a4`.
* Vasile, M. & Conway, B. A., "MGA-DSM transcription for direct
  optimisation of impulsive trajectories", JGCD, 2006. The DSM-leg
  transcription used by `cyclerfinder.search.dsm_leg`.
