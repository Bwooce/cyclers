# Multi-arc basin-selection — explicit flyby-continuity residual (task #162)

> **For agentic workers:** REQUIRED SUB-SKILL: `superpowers:executing-plans` or
> `superpowers:subagent-driven-development`. Checkbox steps; strict TDD (write
> failing test → run **red** → minimal impl → run **green** → commit). Work on
> `main` — **do NOT branch** (project rule). uv-managed venv (no pip). Lint/type
> gate before **every** commit:
> `uv run ruff check .` · `uv run ruff format --check .` · `uv run mypy src tests`.
> Fast suite: `uv run pytest -m "not slow"`; probes that hit the corrector/MBH are
> `slow`.
>
> This plan is the task-level expansion of the **approved** design
> `docs/notes/2026-06-08-multiarc-genome-design.md` — read it first, especially
> §0.1 (the diagnosed mechanism), §3 (the gate), and §4 (the empty-set
> assessment). Where this plan and the design could diverge, **the design wins**,
> except where a design claim no longer holds against live code, which an
> implementer must flag and override with the verified state.
>
> **GOLDEN/HONESTY (binding).** The EXPECTED side of every probe is the row's
> SOURCED V∞ anchor (6.44Gg3 → E 6.44 / M 3.74; from catalogue
> `russell-ch4-6.44Gg3`, Russell 2004 Table 4.13). Emerged V∞ is EVIDENCE, never
> imposed as a residual. Mechanics-gate "expected" values must be defined BY
> CONSTRUCTION (a hand-applied impulse / a known feasible flyby), never by a
> DSM-evaluator output. **No catalogue writeback.**

---

## Goal

Add the intermediate-flyby **V∞-continuity-AND-bend** residual that
`evaluate_dsm_chain` currently omits (it forces heliocentric continuity at
`dsm_leg.py:448` and charges only `Σ ΔV_DSM` at `:474`, with NO turn-feasibility
term — design §0.1). Turn the chain corrector from a single dominated scalar into a
multi-term root-find whose minimum **is** the sourced bend-feasible low-V∞ cycler,
then run the decisive 6.44Gg3 probe against the three-way go/no-go gate (design §3).

Outcome is decisive either way:
- **closes** → first multi-arc closure at sourced anchors;
- **floors with irreducible `flyby_dv` > 0** → quantified empty-set negative
  (publishable);
- **changes basin but no close** → proceed to the hybrid (out of scope here).

All code is **additive in `src/cyclerfinder/search/dsm_leg.py` only**, behind a
default-off `charge_flyby_continuity` flag that preserves every existing caller and
test bit-for-bit. `core/flyby.py` is imported read-only (its `flyby_dv` already
implements the Russell Eq.5.5 powered-SOI surrogate that returns 0 when
ballistic-feasible). **No `core/` edit. No other `src/` file edited.**

---

## VERIFY-FIRST (live state to re-read before any edit)

The design cites these against live code on 2026-06-08; re-read before editing in
case of drift:

1. `evaluate_dsm_chain` (`src/cyclerfinder/search/dsm_leg.py:357`) — the forced
   continuity is `v_depart = leg.v_arrive` (≈ line 448); the scalar objective is
   `total_dv = sum(dv_dsm) + dv_arrive` (≈ line 474). Confirm both still hold.
2. `dsm_chain_decision_vector` / `_unpack` (`dsm_leg.py:582` / `:602`) — layout
   `[t0_sec, vinf_out0, alpha0, beta0, *tof_days, *eta]`. The new direction coords
   append AFTER `*eta`.
3. `dsm_chain_correct._res` (`dsm_leg.py:641`) — today returns `np.array([val])`
   (length-1). `least_squares` natively handles a vector residual; this is the one
   shape change.
4. `core/flyby.py` — `flyby_dv(vin_vec, vout_vec, mu_planet, rp_min)` at `:356`
   (returns 0.0 when `is_ballistic_feasible`, else a positive surrogate Δv);
   `is_ballistic_feasible:318`; `max_bend:40`. Confirm signatures. Per-planet
   `rp_min` and `mu_planet`: source from the registry the way `flyby_dv_for`
   (`:408`) does — re-read it to reuse the planet→(μ, rp_min) resolution rather
   than hardcoding.
5. `russell-ch4-6.44Gg3` in `data/catalogue.yaml` — confirm `free_return_arcs[]`
   (g 2.087 yr / G 4.3191 yr), `vinf_kms_at_encounters` (E 6.44, M 3.74), transit
   262 d. These are the SOURCED anchors + the seed source.

---

## Phase 1 — Mechanics gates for the new residual (RED first, CONSTRUCTED)

> These prove the residual term is correct BY CONSTRUCTION before any probe. Each
> "expected" value is a hand-built feasible/infeasible flyby, never an evaluator
> output. Tests in `tests/search/test_dsm_leg.py` (additive; the existing 14 stay
> green).

- [ ] **Task 1.1** — `test_flyby_continuity_residual_zero_on_feasible_flyby`
  (RED→GREEN). CONSTRUCT a `v∞⁻`/`v∞⁺` pair with equal magnitude and a turn angle
  strictly inside `max_bend(μ_mars, rp_min, |v∞|)`. Assert the new per-flyby
  residual term (magnitude-diff + `flyby_dv`) is `< 1e-9`. Reference = a known
  feasible geometry, NOT the evaluator. Label: mechanics.
- [ ] **Task 1.2** — `test_flyby_continuity_residual_charges_infeasible_turn`
  (RED→GREEN). CONSTRUCT a pair with equal magnitude but a turn angle BEYOND
  `max_bend` (e.g. near-reversal at high V∞). Assert the residual equals the
  hand-computed `flyby_dv` (> 0, within 1e-6). This is the term the old objective
  could not see. Label: mechanics.
- [ ] **Task 1.3** — `test_charge_flyby_continuity_default_off_is_bit_identical`
  (RED→GREEN). With `charge_flyby_continuity=False` (default), `evaluate_dsm_chain`
  and `dsm_chain_correct` return bit-identical results to the current code on the
  existing E-M-E mechanics fixture (the residual is still the length-1 scalar). The
  regression that guarantees no existing caller/test changes. Label: regression.

## Phase 2 — Implement the residual (additive, default-off)

- [ ] **Task 2.1** — Extend the genome layout. In `dsm_chain_decision_vector`,
  `_unpack`, `sequence_keyed_bounds`, `DsmBounds`: append `2·(n_legs−1)` coords
  `[*alpha_int, *beta_int]` (the intermediate-leg departure-V∞ direction for legs
  1..n_legs−1; leg 0 already has `alpha0/beta0`). Bounds: `alpha∈[-π,π]`,
  `beta∈[-π/2,π/2]` (mirror the leg-0 box). Guard arity. Run the existing
  `test_sequence_keyed_bounds_*` — they must still pass (default path) or be
  extended to cover the new arity.
- [ ] **Task 2.2** — In `evaluate_dsm_chain`, add `charge_flyby_continuity: bool =
  False` and the intermediate direction coords. When True: at each intermediate
  flyby, the next leg departs on a V∞ vector built from `(|v∞⁻|, alpha_int,
  beta_int)` (NOT forced to `leg.v_arrive`); compute the per-flyby residual term
  `flyby_dv(v∞⁻_vec, v∞⁺_vec, μ_target, rp_min_target)` using the registry
  resolution from `flyby_dv_for`. Magnitude continuity is automatic if the genome
  sets `|v∞⁺| = |v∞⁻|` — enforce by construction (depart direction free, magnitude
  inherited), so the only flyby residual is `flyby_dv` (the bend term). Return a
  VECTOR residual `[*ΔV_DSM, *flyby_dv_per_flyby, arrival]` via a new
  `residual_vector` field on `DsmChainResult` (additive; `total_dv_kms` stays the
  scalar sum for back-compat / audit).
- [ ] **Task 2.3** — In `dsm_chain_correct._res`, when `charge_flyby_continuity`,
  return `result.residual_vector` instead of `[total_dv]`. `converged` =
  `max(|residual_vector|) < tol_kms` in that mode (every DSM AND every flyby_dv
  below tol), still residual-magnitude only. Thread `charge_flyby_continuity`
  through `dsm_chain_correct` and `make_dsm_chain_step`.
- [ ] **Task 2.4** — Lint/type/regression gate: `ruff check` + `ruff format
  --check` + `mypy src tests` clean; `uv run pytest tests/search/test_dsm_leg.py
  -m "slow or not slow"` — the existing 14 stay green, the 3 new mechanics gates
  pass. **Commit:** `dsm_leg: explicit flyby V∞-continuity+bend residual (default-off)`.

## Phase 3 — Symmetric descriptor seed (thin slice of design #2)

- [ ] **Task 3.1** — Add a probe-side helper (in the test module, NOT new src
  surface unless reused) `symmetric_arc_seed_644gg3()` that maps
  `russell-ch4-6.44Gg3` `free_return_arcs[]` → the E-M-E-M-E genome:
  ToF (262, 500.4, 262, 1315.5) d (the #153 sourced decomposition), η=0.5 per leg,
  departure V∞ = 6.44 km/s tangential, intermediate direction coords initialised to
  continue the arrival direction (so the seed = the old forced-continuity chain at
  t=0). Long-leg branch forced `(1,low)` (the #157 winner).

## Phase 4 — THE DECISIVE PROBE (slow, qualitative-assert + gate verdict)

- [ ] **Task 4.1** — `test_dsm_644gg3_flyby_continuity_probe` (slow). Single
  `dsm_chain_correct(..., charge_flyby_continuity=True, max_revs=3,
  rev_branch_per_leg=...)` from the symmetric seed on `Ephemeris("circular")`
  first (cheapest; the continuation seed if it closes), then DE440. Assert only
  that the search RUNS and returns a finite, fully-audited result (the scientific
  verdict is REPORTED, not gated). PRINT the verbatim block: `total ΔV`, per-flyby
  `flyby_dv`, per-leg `ΔV_DSM`, emerged V∞ in/out vs sourced 6.44/3.74, branches,
  η, ToF, hops if MBH.
- [ ] **Task 4.2** — If the single run shows movement toward the anchors, wrap in
  MBH (`make_dsm_chain_step(..., charge_flyby_continuity=True)`, cauchy, rng_seed=6,
  ≤120 hops, stall 60) and re-probe. Deterministic; seed echoed.
- [ ] **Task 4.3** — Apply the **three-way gate** (design §3) and write the verdict
  to a results note `docs/notes/2026-06-08-multiarc-continuity-probe.md`:
  - **CLOSE:** `max(residual_vector) < 0.1` AND emerged V∞ within tol (E 6.44 ±1.0,
    M 3.74 ±0.5) → first multi-arc closure; record, then (separate task) hybrid +
    S1L1 + continuation-to-DE440.
  - **EMPTY-SET:** solver converges but irreducible `min max(residual_vector)`
    stays > ~1 km/s across a small seed sweep AND emerged V∞ never approaches
    anchors → record the quantified `flyby_dv` floor as the publishable negative;
    STOP (do not assemble hybrid).
  - **AMBIGUOUS:** basin moves (emerged V∞ materially toward anchors, floor below
    #157's 26.9 km/s) but no close → the objective fix works; next experiment =
    full hybrid (design #5). Record and hand off.
- [ ] **Task 4.4** — Lint/type gate clean. **Commit (code+test):**
  `dsm_leg: 6.44Gg3 flyby-continuity probe (multi-arc basin-selection #162)`.
  **Commit (results note):** `docs: 6.44Gg3 flyby-continuity probe verdict (#162)`.

---

## Decisive gate (one-line summary)

On 6.44Gg3, with the explicit per-flyby `flyby_dv` residual charged and the long
leg forced `(1,low)`, does `dsm_chain_correct` reach `max(residual_vector) <
0.1 km/s` with emerged V∞ at E 6.44 / M 3.74? **Yes → first multi-arc closure.
No, with an irreducible flyby_dv floor → quantified empty-set (publishable). Basin
moves but no close → assemble the hybrid.** No path returns "inconclusive": a floor
number + emerged V∞ is always decisive evidence.

## Out of scope (explicit)

- The hybrid (design #5) and S1L1 4.991gG2 re-probe — only if Phase 4 fires CLOSE
  or AMBIGUOUS.
- Multi-arc continuation to DE440 (design #4) — only if a circular-coplanar
  multi-arc closure emerges in Task 4.1.
- Any `core/` edit, any catalogue writeback, any change to the default
  (`charge_flyby_continuity=False`) path.
