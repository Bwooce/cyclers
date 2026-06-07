# One-DSM-per-leg genome (the Takao eta-coordinate) — task #150

**Date:** 2026-06-07
**Code:** `src/cyclerfinder/search/dsm_leg.py`, `tests/search/test_dsm_leg.py`
**Context read:** `docs/notes/2026-06-07-takao-2025-mpga-1dsm-mining.md` (the
transcription), `src/cyclerfinder/search/free_return.py` (#137 architectural
sibling), `docs/notes/2026-06-07-hughes-2014-fast-mars-free-returns-mining.md`
(WHY low-V∞ E-V-M needs broken-plane/DSM), `docs/notes/multi-arc-classification.md`
§5/§7/§12 (6.44Gg3 and S1L1 are multi-arc, not single-ellipse),
`docs/notes/2026-06-07-mbh-wrapper.md` (the optimiser + Gate-3 negative).

## What this is

The multi-arc-return primitive the S1L1 / Jones closure blocker needs: an
**interior impulse on each leg** (a deep-space maneuver, DSM). Where
`free_return.py` models an E→M→E transfer as a *single* heliocentric ellipse
(one shape `(a, e)` shared by both legs), this module lets a leg follow a
different ballistic arc on its front fraction than on its back fraction — i.e. a
leg that is NOT a piece of one repeating ellipse. The multi-arc rows
(`russell-ch4-*`, S1L1) are exactly the cases where the sourced geometry is *two
generic-return arcs*, so a single-ellipse genome has no sourced-anchor basin
(the #137 / MBH Gate-3 negative result). The one-DSM-per-leg genome can
represent that topology.

## The transcription (Takao 2025, Eq.6-7)

`dsm_leg(r0, v0, tof, eta, target_r, mu)`:

1. propagate `(r0, v0)` ballistically (heliocentric 2-body, `core/kepler.propagate`)
   for `eta·tof` → DSM point `(r_12, v_12)`;
2. Lambert (`core/lambert.lambert`, single-rev) from `r_12` to `target_r` over
   `(1−eta)·tof` → post-impulse departure `v_21`, arrival `v_22`;
3. DSM impulse `ΔV_DSM = ‖v_21 − v_12‖` (Eq.6);
4. incoming hyperbolic velocity at the next body = `v_22 − v_planet` (Eq.7).

`eta ∈ [0, 1]` is the **only genuinely new genome coordinate** vs the existing
free-return / ballistic correctors; everything else is reused from `core/`.

## Conventions mirrored from `free_return.py`

- Frozen-dataclass results carrying both the CONSTRAINED quantity (the residual /
  total ΔV the corrector drives) and the FREE / EMERGED evidence (per-leg ΔV
  breakdown, per-body emerged V∞). **The sourced V∞ is never imposed** — it
  emerges from the converged genome and is comparable as evidence (golden-rule
  separation).
- `converged` decided by residual MAGNITUDE alone (`total_dv < tol_kms`), BY
  DESIGN; the `least_squares` `success` flag is audit-only.
- Full audit trail on the result (per-leg ΔV, DSM states, emerged V∞ in/out,
  solver nfev).
- RNG-free / deterministic everywhere; MBH supplies the hops, this supplies the
  refinement (the corrector is the local-solve closure `make_dsm_chain_step`).

## Module surface

- `dsm_leg(...)` — the leg primitive (Eq.6-7).
- `evaluate_dsm_chain(...)` — chained N-leg evaluator (Eq.3 epoch chaining, Eq.5
  departure parameterisation, Eq.14-15 objective `Σ ΔV_DSM + ΔV_arr`, P-FB term
  deliberately left to `core/flyby.py` at scoring time).
- `sequence_keyed_bounds(...)` — Takao Appendix A.1-A.3 automatic bounds
  (`alpha∈[-π,π]`, `beta∈[-π/2,π/2]`, `eta∈[0,1]`; inner-pair ToF `[30 d, P_s+P_H]`,
  outer-pair ToF `[0.3 P_H, 1.3 P_H]`; departure V∞ `[1, 5.1]` km/s default).
- `dsm_chain_decision_vector` / `dsm_chain_correct` — flat genome
  `[t0, vinf_out0, alpha0, beta0, *tof_days, *eta]` + bounded `least_squares`
  corrector (clips an out-of-box MBH-hopped seed into the interior so `trf`
  never rejects it).
- `make_dsm_chain_step(...)` — MBH adapter (imports + calls `mbh`; never edits it).

## Mechanics gates (constructed, label = mechanics) — VERBATIM RESULTS

`uv run pytest tests/search/test_dsm_leg.py -m 'slow or not slow'` →
**7 passed in 17.79s** (6 mechanics + 1 slow probe).
`ruff check` + `ruff format --check` + `mypy` on both files: **clean**.

Each gate's "expected" value is defined BY CONSTRUCTION by a known transfer
(a Lambert arc, or a hand-applied impulse), NOT by a DSM-evaluator output.

### Gate 1 — η-degeneracy regression
- `test_eta_consistency_zero_dsm_when_v0_is_lambert_velocity`: when `v0` IS the
  Lambert arc's start velocity, `ΔV_DSM < 1e-6 km/s` for **every** η in
  {1e-4, 0.25, 0.5, 0.75, 0.9999} (measured ~3e-14 to ~5e-8), and `v_arrive`
  equals the pure-Lambert arrival — the one-DSM leg reduces to a pure Lambert
  leg whenever the front arc already follows the Lambert solution. **PASS.**
- `test_eta_to_zero_reduces_to_start_endpoint_discontinuity`: with η~0 and `v0 =
  v_lambert + known_offset`, `ΔV_DSM` recovers `‖known_offset‖` to < 1e-3 km/s
  (the start-endpoint discontinuity). **PASS.**
- `test_eta_endpoints_are_rejected_as_singular`: exact η∈{0,1} and out-of-range η
  raise (zero-duration ballistic/Lambert arc). **PASS.**

### Gate 2 — constructed broken-plane two-impulse transfer recovered exactly
- `test_constructed_broken_plane_transfer_recovered`: from an out-of-plane
  departure, propagate `eta·tof`, apply a KNOWN out-of-plane impulse
  `dv_known=[0.8,-0.5,0.3]`, propagate the rest to a target. Feeding the endpoints
  back to `dsm_leg` recovers `‖dv_known‖` (< 1e-6 km/s), the DSM position (< 1e-3
  km), the post-impulse departure (< 1e-6 km/s) and the arrival velocity (< 1e-6
  km/s). The non-zero z-impulse is exactly the broken-plane case a single-ellipse
  genome cannot represent and this primitive can. **PASS.**

### Gate 3 — determinism / audit trail
- `test_chain_is_deterministic_and_carries_full_audit_trail`: two identical chain
  evaluations are byte-identical (rng-free); the result carries the per-leg ΔV
  breakdown, emerged per-body V∞ (in/out), η per leg, and DSM states;
  `total_dv == Σ ΔV_DSM` and `max_residual_kms` aliases the objective. **PASS.**
- `test_sequence_keyed_bounds_layout_and_eta_box`: A.1-A.3 bounds have correct
  arity (8 for E-M-E), `eta∈[0,1]`, `alpha∈[-π,π]`, `beta∈[-π/2,π/2]`, inner-leg
  ToF floor 30 d. **PASS.**

## THE PROBE — 6.44Gg3 sourced-anchor gate (slow, wall 16.9 s ≪ 10 min)

`russell-ch4-6.44Gg3` (Russell 2004 Table 4.13) is a **3-synodic MULTI-ARC**
cycler: two generic Earth-Earth free-return arcs (g(2.087 yr) + G(4.3191 yr))
bracketing the E→M outbound (262 d) and M→E return (262 d). **SOURCED anchors:**
V∞ E = 6.44, M = 3.74 km/s, aphelion 1.54 AU. The single-ellipse free-return
genome PROVABLY cannot reach this (MBH Gate-3: emerged V∞ E=3.01 / M=3.06, far
off-anchor).

Probe: an explicit `E→M→E` chain with an interior DSM on each leg, seeded at the
sourced 262-day transits and the sourced Earth departure V∞ (6.44 km/s,
tangential azimuth), driven by the MBH wrapper (cauchy, seed 6, ≤120 hops,
stall 60). The MBH hop perturbs **every** gene that enters the ΔV objective —
t0, departure V∞ magnitude/azimuth/elevation, per-leg ToF, per-leg η — and the
corrector moves all of them.

### Emerged-vs-sourced (rng_seed=6) — verbatim

```
feasible=False   total_dV = 9.4012 km/s
hops attempted/accepted = 61/0   (stopped on stall)
emerged V_inf_in   {leg1(E→M arrival at M): 9.04,  leg2(M→E arrival at E): 7.09}   sourced E=6.44, M=3.74
emerged V_inf_out  {dep E: 6.65,  M: 9.04,  E: 7.09}
per-leg dV_DSM = (3.43, 5.97) km/s
eta per leg    = (0.633, 0.435)
tof days/leg   = (264.3, 264.8)
wall = 16.9 s
```

### Verdict — a clean NEGATIVE that bounds the frontier honestly

The departure V∞ converged near the sourced Earth anchor (6.65 vs 6.44 km/s),
but the **flyby V∞ continuity never closes**: best total ΔV is 9.40 km/s (vs the
0.1 km/s tol), 0 of 61 hops accepted, and the emerged encounter V∞ (9.04 / 7.09)
sit far from the multi-arc anchors (6.44 / 3.74). This is the **expected and
honest result**: the direct `E→M→E` one-DSM chain with 262-day transits is *not*
the full 3-synodic multi-arc topology — 6.44Gg3 has **two intermediate
Earth-Earth loop arcs** (g + G, the 2.087 yr + 4.319 yr free-returns) between the
outbound and return transits, which a 2-leg chain does not model. The DSM
primitive is mechanically correct (Gates 1-3 prove it); the negative says the
*minimal* E-M-E chaining does not reproduce a topology whose defining feature is
the pair of resonant loop arcs.

**What this DOES establish (the positive in the negative):** the one-DSM-per-leg
primitive is built, validated, and MBH-drivable; the missing piece for 6.44Gg3
is not the interior-impulse mechanic but the **multi-leg sequence** — chaining
the two Earth-Earth loop arcs (e.g. `E→E→E→M→M` per the row's
`sequence_canonical`, each loop arc its own DSM leg). The evaluator already
strings arbitrary `len(sequence)-1` legs, so that is a sequence/seed change, not
new mechanics. Recorded as the follow-up.

**No catalogue writeback** (per task constraint; the row stays multi-arc with its
sourced anchors).

## Follow-ups

1. **Full 6.44Gg3 topology:** drive the evaluator on the row's `E-E-E-M-M`
   sequence (the two generic loop arcs as their own DSM legs) rather than the
   minimal `E-M-E`, seeded from the `free_return_arcs` ToFs (2.087 yr / 4.3191 yr).
   This is the topology whose two-arc nature §12 predicted is the real closure
   key; the minimal chain here is the control.
2. **Oberth-aware P-FB cost (Takao Eq.11/12/13):** the transcription flags Eq.11
   as a strict improvement over `core/flyby.py::dv_from_turn_deficit`. Out of
   scope for #150 (this module charges only the DSM impulses; the flyby cost is
   the scorer's job) but the natural next adoption.
3. **s1l1-2syn-em-cpom two-arc probe (5.65 / 3.05):** the optional second
   application; same machinery, the 2-synodic two-arc geometry.
