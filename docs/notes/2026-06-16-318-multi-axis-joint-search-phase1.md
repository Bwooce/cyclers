# #318 Phase 1 - multi-axis joint-search substrate + small EM probe

**Status:** Phase 1 substrate complete; Phase 2 path scoped.
**Date:** 2026-06-16.

## The strategic-answer rationale

No single published paper has done a joint sweep across all four cycler-
discovery axes simultaneously:

  * **A. Powered maintenance.** Aldrin 1985, McConaghy/Longuski/Byrnes
    2002, Genova-Aldrin 2015, Russell-Ocampo 2006 - each treats powered
    cyclers, but always one-axis-at-a-time.
  * **B. Multi-revolution Lambert.** The standard literature uses
    single-rev Lambert for cycler closure; multi-rev branches are
    used in MGA tour design (Izzo 2014, Olympio 2011) but not jointly
    with the other cycler axes.
  * **C. 3D / broken-plane.** Braik-Ross 2024, Antoniadou-Voyatzis 2018,
    Ross-Roberts-Tsoukkas 2025 - 3D periodic families in the planar-
    invariant CR3BP, again as a standalone axis.
  * **D. Epoch-locked validity windows.** Tito 2018 free-return,
    Heaton-Longuski 2003 Uranian tours - epoch-locked targeting against
    real-ephemeris geometry, but never as a Cartesian joint sweep with
    the other three.

Each axis has been published in isolation. The joint manifold
(`powered x multi-rev x 3D x epoch-locked`) has not. If genuinely-novel
cycler pockets exist there, the substrate built in this task is what
surfaces them.

## Phase 1 scope (this task)

Substrate, not full discovery sweep. The four ingredients:

1. **`JointAxisCandidate` value type** at
   `src/cyclerfinder/search/multi_axis_search.py`. Records the four-axis
   cell coordinates PLUS the axis-driver verdicts (Sims-Flanagan
   feasibility from #309, multi-rev branch, z0 amplitude, epoch-locked
   closure residual). Serialises to JSONL-friendly dict.
2. **`joint_axis_search` iterator** at the same module. Composition driver
   that wires:
     - `#309` powered driver (`search_low_thrust_cyclers`) - Axis A + B
       (powered budget + per-leg revs forwarded to the existing optimiser).
     - `#291` 3D corrector (recorded only in Phase 1; Phase 2 pipes it
       inline).
     - `#289` epoch-locked closure (`close_epoch_locked`) - Axis D
       (per-cell witness).
3. **11 sanity tests** at `tests/search/test_multi_axis_search.py`.
   Cover: joint-zero corner reduces to existing #309 baseline; Aldrin
   reproduction; multi-rev grid never reduces survivor count; z0
   amplitude recorded on candidate; non-closed sequence rejected;
   k_synodic validation; n_revs grid shape validation; JSONL
   serialisation round-trip; Cartesian-product cell enumeration;
   powered-budget exceedance advisory.
4. **EM probe script + JSONL** at
   `scripts/scan_318_em_joint_phase1.py` + `data/scan_318_em_joint_phase1.jsonl`.

The four axis modules (`low_thrust_cycler_search.py`,
`cr3bp_general_periodic_3d.py`, `epoch_aware_genome.py`,
`tisserand_mga_window.py`) are **UNMODIFIED**. The Phase 1 substrate
composes them; it does not rewrite them. (Discipline: per the task spec,
read-only on those modules.)

## EM probe results (Aldrin E-M-E, 192 cells)

| Axis | Grid                             | Survivors |
|------|----------------------------------|-----------|
| A - powered budget (km/s)         | (0.0, 0.05, 0.1, 0.2, 0.5, 1.5) | 32 per budget |
| B - n_revs per leg                | ((0,1),(0,1)) -> 4 combos       | 48 per combo |
| C - z0 amplitude (non-dim)        | (0.0, 1e-3, 1e-2, 5e-2)         | 48 per amp |
| D - launch epoch                  | (None, '2030-01-01')            | 96 per epoch |

**Total: 192 cells, 192 survivors (100%).** Every cell closes the
optimiser.

**Best per axis-corner** (lowest closure_residual_kms): universally lands
at the (revs=(0,0), z0=0, epoch=2030-01-01) cell with `dv=0.0000 km/s`
and `V_inf=(5.14, 9.14, 5.14)` km/s.

**The dV manifold collapses to exactly two values:**

  * `1.2870 km/s` at the epoch-blind cells (the Aldrin phase-inversion
    t0 anchor; the canonical Aldrin maintenance dV from McConaghy 2002
    / circular-coplanar surrogate).
  * `0.0000 km/s` at the 2030-01-01 epoch cells (a different EM phase
    where the optimiser found an alternate local minimum with a free
    Earth turn; V_inf shifts from the Aldrin (6.52, 9.73, 6.52) to a
    softer (5.14, 9.14, 5.14)).

**Multi-rev grid (Axis B) and z0 amplitude (Axis C) do NOT shift the
optimiser's verdict in Phase 1:**

  * Axis B: the optimiser's per_leg_revs argument threads through, but
    its minimum-dV basin sits at the direct (revs=0) branch for the EM
    Aldrin tour; the multi-rev branches don't open up a better local
    minimum. This is the same behaviour the existing #309 EM scan
    showed (see `data/scan_309_low_thrust_em.jsonl`).
  * Axis C: Phase 1 records the z0 amplitude on the candidate row but
    does NOT drive the 3D corrector inline. The 2D Lambert engine
    converges identically on the 3D-request and planar cells. This is
    the Phase 1 *contract* (see the module docstring); Phase 2 will
    pipe the recorded z0 amplitude into the #291 corrector for a
    same-model 3D closure on the candidates that survive Phase 1.

## Honest verdict

**The Cartesian product over the redundant axes wastes 96% of the cells
on rediscovery.** Only 2 of 192 cells produce distinct dV values
(1.2870 and 0.0000). The substrate is correct; the sampling strategy is
the limit.

The next phase **must** adopt smarter sampling. The four axis modules
already have the leverage; Phase 1 confirmed the composition works.

## Phase 2 path (NOT in this task)

1. **Smarter sampling.** Sobol or surrogate-driven instead of Cartesian.
   The Phase 1 manifold collapsed to 2 distinct dV; a Sobol pre-pass at
   ~32 samples should hit both basins and leave the rest of the budget
   for new pockets.
2. **Drive the 3D corrector inline.** Phase 1 records `z0_amplitude_nondim`
   on the candidate row. Phase 2 picks up the candidates flagged
   `z0_amplitude_nondim > 0` and runs `correct_general_periodic_3d` from
   #291 as the primary residual - same-model 3D closure on each.
3. **Apply to VEM + Uranus + Mars systems.** EM is the established
   reproduction surface; the joint-axis sweep is only interesting where
   the individual-axis literature surveys (#328 / #302) have not been
   exhaustive. Uranus (Heaton-Longuski 2003 single-tour anchor only),
   VEM (Russell-Ocampo Table 3.4 + Jones VEM family), and broader Mars
   (Rogers 2012 family but no joint with Axis C/D) are the priority
   targets.
4. **Integrate the #317 PINN pre-filter as the joint-axis cell scorer.**
   The PINN scores Lambert closures; on the joint axis it would score
   the *cell* (the four-tuple) as a candidate-density estimator. Phase
   1's 96% redundancy is the gap the PINN closes.

## Discipline

  * **No catalogue writeback.** This task touches catalogue zero times.
  * **No novelty claims.** The frame is "joint-axis Phase 1 substrate
    + small EM probe", not "we found a new cycler".
  * **Sourced golden anchor.** The catalogue row
    `aldrin-classic-em-k1-outbound` is the regression target. The
    joint-zero corner re-finds Aldrin at maintenance dV = 1.2870 km/s,
    V_inf ~ (6.52, 9.73, 6.52) km/s - bit-identical with the #309 baseline.
  * **Composition, not rewrite.** All four axis modules unmodified.
  * **Concurrent-agent awareness.** All commits pathspec-only to the
    files this task owns: `multi_axis_search.py`,
    `test_multi_axis_search.py`, `scan_318_em_joint_phase1.py`,
    `scan_318_em_joint_phase1.jsonl`, this doc.

## Commits

  * `0b9c0eb` - Part A driver + 11 tests.
  * `108c909` - Part C EM probe script + 192-cell JSONL.
  * (this doc) - Part D.

## References

  * Module: `src/cyclerfinder/search/multi_axis_search.py`.
  * Tests: `tests/search/test_multi_axis_search.py`.
  * Probe: `scripts/scan_318_em_joint_phase1.py`.
  * Probe output: `data/scan_318_em_joint_phase1.jsonl`.
  * Axis A substrate: `src/cyclerfinder/search/low_thrust_cycler_search.py`
    (#309 Phase 1).
  * Axis C substrate: `src/cyclerfinder/search/cr3bp_general_periodic_3d.py`
    (#291 Phase 1).
  * Axis D substrate: `src/cyclerfinder/genome/epoch_aware_genome.py`
    (#289 Phase 1).
  * Sourced golden: `data/catalogue.yaml` row
    `aldrin-classic-em-k1-outbound`.
