# Multi-rev resonant-leveraging releg — IMPLEMENTATION PLAN (#465)

**Date:** 2026-06-26
**Design draft:** `docs/superpowers/specs/2026-06-26-465-multirev-leveraging-releg-design-draft.md`
**Status:** PLAN — bite-sized TDD tasks. No code written by this doc.

Each task: a failing test first, the minimal implementation, a verify step, a
single pathspec commit. **Reuse over rebuild** — the per-hop leverage leg
(`search/leveraging_leg.py`, #179), the quadrature floor/golden (`search/vilm.py`),
the `Releg` protocol + driver + prefilter (`search/releg_solver.py`,
`releg_moontour.py`, `moon_prune.py`) all exist and are tested. The new work is the
**chain orchestration** (walk V∞_H→T via N resonant hops, sum δ) — ONE new backend
class + its tests, behind the existing protocol. No new optimiser, no new cost
model, no driver rewrite.

Conventions: `uv run pytest <path> -q`; `uv run ruff check . && uv run ruff format
--check .` before every commit (`feedback_run_ruff_before_commit`). **Pathspec
commits only**, never `git add -A` (a sibling agent edits
`genome/qp_tori_arclength.py` concurrently — `feedback_concurrent_agent_git_rules`).
No `Co-Authored-By`. Work on `main`, never branch. Prefix multi-step bash with
`date -Iseconds`. If the pre-commit mypy hook fails ONLY on the sibling's
in-progress file, commit `--no-verify` for markdown-only commits; for code commits,
ensure your own files pass and use a pathspec commit (the hook runs on staged
files).

---

## Task 0 — (gate, read-only) confirm the chain primitive + seam contract

**Why first:** the whole plan rests on `evaluate_leveraging_leg` (#179) being a
usable per-hop V∞-shaping primitive and the `Releg`/driver seam taking a swappable
backend. Confirm before building.

- **Verify (read-only, no commit):**
  - `grep -n 'def evaluate_leveraging_leg\|gamma_floor_kms\|vout_t\|near root' src/cyclerfinder/search/leveraging_leg.py` — confirm the hop primitive takes `(moon, n_moon_revs, m_sc_revs, vinf_in_kms, vinf_out_target_kms, exterior)` and returns `(dv_dsm_kms, vinf_out_kms, converged, gamma_floor_ok)`.
  - `grep -n 'class Releg\|isinstance(releg, BallisticReleg)\|vinf_target_in\|vinf_depart_mag' src/cyclerfinder/search/releg_solver.py src/cyclerfinder/search/releg_moontour.py` — confirm the driver selects "powered iff not BallisticReleg" and pins/retargets via the two kwargs.
- **Output:** one-line note confirming both. If `evaluate_leveraging_leg` cannot
  hit an arbitrary `vinf_out_target` (only specific resonances), Task 2's descent
  must search resonances — note it.
- **No commit.**

---

## Task 1 — chain orchestration helper (pure, the V∞ descent)

**Goal:** the resonant-hop descent that walks V∞_H → T at a moon and sums δ —
isolated and unit-tested BEFORE wrapping it in a `Releg`.

- **Failing test:** `tests/search/test_leveraging_chain.py::test_chain_walks_vinf_down_at_floor`
  - At Europa, walk V∞ from 1.8 km/s down to 0.77 km/s (the sourced A6 endgame
    bounds, mining note 436-438). Assert the chain: converges, ends within tol of
    0.77, and its summed ΔV is `≥ vilm._leverage_dv_kms`-style continuous floor
    AND `≤ 1.25 × floor` (the +20% finite-chain band the published 154 vs 128 m/s
    sets — EXPECTED edges both sourced: floor=Eq.13, ceiling=published discrete).
  - `test_chain_zero_walk_is_zero_dv`: target == natural V∞ ⇒ no hops ⇒ ΔV ≈ 0.
  - `test_chain_infeasible_returns_none`: an unreachable target (deeper than any
    feasible hop sequence) returns a sentinel "infeasible", never a fabricated ΔV.
- **Minimal impl:** `src/cyclerfinder/search/leveraging_chain.py`
  - `walk_vinf_down(moon, vinf_hi_kms, vinf_lo_kms, *, exterior=True, max_hops, max_revs) -> ChainResult` (frozen dataclass: `total_dv_kms`, `hops: tuple`, `vinf_end_kms`, `total_revs`, `converged`).
  - Greedy "zigzag low" descent (mining note p.11): at each step pick the cheapest
    FEASIBLE `evaluate_leveraging_leg(...)` (small n:m grid) whose `vinf_out` steps
    toward `vinf_lo`; stop at `vinf_lo` (tol), at `max_hops`/`max_revs` (kill), or
    when no feasible hop advances (infeasible). Sum `dv_dsm_kms`; honour each hop's
    `gamma_floor_ok`.
- **Verify:** `uv run pytest tests/search/test_leveraging_chain.py -q`; ruff.
- **Commit (pathspec):** `git add src/cyclerfinder/search/leveraging_chain.py tests/search/test_leveraging_chain.py && git commit -m "search: #465 resonant-hop V-inf descent chain (leveraging_chain)"`

---

## Task 2 — `MultiRevLeveragingReleg` backend (behind the `Releg` protocol)

**Goal:** wrap the Task-1 chain as a third powered `Releg`, swappable into the
driver with no driver change.

- **Failing test:** `tests/search/test_releg_solver.py::test_multirev_leveraging_releg_*`
  - `test_multirev_releg_zero_retarget_matches_ballistic`: with `vinf_target_in
    is None`, the backend returns `dv_kms ≈ 0` and a V∞ chain matching the
    lowest-energy ballistic branch (the regression limit, like the other backends).
  - `test_multirev_releg_retargets_via_chain`: on a Ganymede→Europa leg with the
    arrival retargeted down to a common `T`, assert `feasible`, arrival V∞ within
    tol of `T`, `dv_kms > 0`, AND `dv_kms` is **below** the single-impulse
    `DsmReleg` cost on the same leg (the whole point: the chain is cheaper than one
    impulse — the multi-rev win, asserted as `chain < dsm`).
  - `test_multirev_releg_dv_geq_leverage_floor` (golden, sourced): `dv_kms ≥
    vilm._leverage_dv_kms`-equivalent continuous floor − tol (a finite chain cannot
    beat the continuous minimum; EXPECTED traces to Eq.13, `feedback_golden_tests_sourced_only`).
- **Minimal impl:** add `class MultiRevLeveragingReleg` to
  `src/cyclerfinder/search/releg_solver.py`
  - `solve(...)`: lambert seed → natural arrival V∞_H (same seed logic as
    `DsmReleg`/`LowThrustReleg`, lines 230-247) → if no retarget, `dv=0` ballistic
    limit → else `leveraging_chain.walk_vinf_down(moon_B, V∞_H, vinf_target_in)` →
    `RelegResult(vinf_out=vinf_depart_mag, vinf_in=vinf_target_in, dv_kms=chain.total_dv_kms, feasible=chain.converged)`.
  - Moon identity at the arrival side: the backend needs the arrival MOON name to
    pick the leverage body. Thread it in (the driver knows `sequence[k+1]`); add an
    optional constructor arg or a per-call moon hint via a thin driver tweak in
    Task 3. Keep `solve`'s signature protocol-compatible.
- **Verify:** `uv run pytest tests/search/test_releg_solver.py -q`; ruff.
- **Commit (pathspec):** `git add src/cyclerfinder/search/releg_solver.py tests/search/test_releg_solver.py && git commit -m "search: #465 MultiRevLeveragingReleg backend — chained VILM endgame releg"`

---

## Task 3 — driver wiring: pass the leverage moon, keep the contract

**Goal:** let `close_powered_cycle` give the new backend the arrival-moon name so
the chain leverages at the right body — WITHOUT breaking the other backends.

- **Failing test:** `tests/search/test_releg_moontour.py::test_multirev_galilean_positive_control_in_band`
  - Io-Europa-Ganymede-Io skeleton (registry positive control, links at vinf=4):
    with `MultiRevLeveragingReleg`, assert the cycle closes (continuity below the
    0.05 km/s gate, by construction) AND total `dv_kms` is **inside the powered
    dv-band** (`verify/dv_band_acceptance`) — the in-band closure #449/#464 missed.
  - `test_multirev_galilean_cheaper_than_dsm`: same skeleton, assert the
    multi-rev total ΔV is well below the `DsmReleg` total (13.18 km/s) — the gate
    evidence, reproduced end-to-end.
  - Regression: `test_uranus_disjoint_prefiltered_empty` still EMPTY with the new
    backend (prefilter skips before any chain solve — chaining can't bridge
    disjoint contours).
- **Minimal impl:** in `src/cyclerfinder/search/releg_moontour.py`
  - `_close_at_target` passes the arrival moon `sequence[k+1]` to the backend.
    Add an optional `arrival_moon: str | None = None` kwarg to the `Releg.solve`
    protocol (default `None`; the ballistic/DSM/SF backends ignore it — they read
    V∞ off the Lambert seed and don't need the body name). The leveraging backend
    uses it to pick the leverage body. This is a protocol ADD (backwards-compatible
    default), not a rewrite.
  - Confirm `DsmReleg`/`LowThrustReleg`/`BallisticReleg` still pass (they ignore
    the new kwarg).
- **Verify:** `uv run pytest tests/search/test_releg_moontour.py tests/search/test_releg_solver.py -q`; ruff.
- **Commit (pathspec):** `git add src/cyclerfinder/search/releg_moontour.py src/cyclerfinder/search/releg_solver.py tests/search/test_releg_moontour.py && git commit -m "search: #465 wire arrival-moon hint into releg driver for the leveraging chain"`

---

## Task 4 — Saturnian positive control + the in-band golden table

**Goal:** lock the §6 in-band numbers as a sourced golden and confirm a second
system closes in-band.

- **Failing test:** `tests/data/test_golden_multirev_leveraging.py::test_*`
  - `test_leverage_only_costs_match_sourced`: assert `vilm._leverage_dv_kms`-derived
    per-transfer leveraging-only costs reproduce the design-draft §6 table
    (Ganymede↔Europa 0.305, Europa↔Io 0.409, Ganymede↔Io 0.735, Titan↔Rhea 0.334,
    Rhea↔Dione 0.196 km/s) to a documented tol. EXPECTED side: these are the Eq.13
    begingame+endgame quadratures (sourced), NOT a fresh self-computed number — the
    test asserts the DECOMPOSITION (full = lev-only + escape + capture) against the
    sourced Table-1 A2 totals (`feedback_golden_tests_sourced_only`).
  - `tests/search/test_releg_moontour.py::test_saturnian_positive_control_in_band`:
    Titan-Rhea-Dione-Titan skeleton closes in-band with `MultiRevLeveragingReleg`.
- **Minimal impl:** none beyond Tasks 1-3 (these are assertions over existing
  code + a small golden constant table). If a `data/golden/` YAML is preferred,
  add `data/golden/multirev_leveraging_inband.yaml` with the sourced lev-only costs
  + their A2 Table-1 provenance lines; else inline the sourced constants in the
  test with `# source:` comments tracing to A2/A6.
- **Verify:** `uv run pytest tests/data/test_golden_multirev_leveraging.py tests/search/test_releg_moontour.py -q`; ruff.
- **Commit (pathspec):** `git add tests/data/test_golden_multirev_leveraging.py [data/golden/multirev_leveraging_inband.yaml] tests/search/test_releg_moontour.py && git commit -m "data: #465 sourced in-band leveraging golden + Saturnian positive control"`

---

## Task 5 — capability tags + powered-empty re-stamp (Uranus/Neptune, stronger)

**Goal:** register the multi-rev capability in the subsumption partial order and
confirm the stronger powered-empty re-stamp for the disjoint-contour systems.

- **Failing test:** `tests/search/test_releg_moontour.py::test_multirev_capability_subsumes_dsm`
  - Assert a `multirev_leveraging_method_capability()` carries tags that subsume the
    DSM releg tags (`multi-rev-leveraging` ⊐ `one-dsm-per-leg`, `leveraging`) via
    `data/method_capability` edges.
  - `test_multirev_uranus_powered_empty_restamp`: on a prefiltered-empty Uranian
    skeleton, `build_powered_empty_restamp` produces a `validate_empty_region`-valid
    record naming the multi-rev method+version (a powered CHAIN that ALSO can't
    bridge is a stronger negative than the single-DSM one).
- **Minimal impl:** in `releg_moontour.py`, add
  `multirev_leveraging_method_capability(*, git_sha)` parallel to the existing
  `powered_releg_method_capability`, with the new tag set; add the
  `multi-rev-leveraging` capability edge(s) to `data/method_capability` if absent.
- **Verify:** `uv run pytest tests/search/test_releg_moontour.py tests/data -q`; ruff.
- **Commit (pathspec):** `git add src/cyclerfinder/search/releg_moontour.py src/cyclerfinder/data/method_capability* tests/search/test_releg_moontour.py && git commit -m "search: #465 multi-rev-leveraging capability tags + stronger powered-empty re-stamp"`

---

## Task 6 — ratchets + capability note

**Goal:** keep frozen-census ratchets green; document the capability.

- **Run ALL ratchets** (`feedback_catalogue_edits_run_all_ratchets`):
  `uv run pytest tests/data tests/search -q`. This plan adds NO catalogue.yaml row
  (the discovery campaign is a separate issue), so census ratchets stay untouched —
  confirm green.
- **Docs:** add a short capability note to `data/OUTSTANDING.md` next to the #449
  note, referencing this design draft + plan: "Multi-rev leveraging releg
  (#465) — CAPABILITY: chained VILM endgame brings the Galilean/Saturnian moon-tour
  IN-BAND (≈1.5 / ≈0.9 km/s/cycle, vs the 13.18 km/s single-DSM close); Uranus/Neptune
  stay a stronger powered-empty (disjoint contours). Golden: Eq.13 leveraging
  quadrature + Europa 154/128 m/s." (`feedback_update_docs_proactively`)
- **Verify:** full `uv run pytest tests/data tests/search -q`; ruff.
- **Commit (pathspec):** `git add data/OUTSTANDING.md && git commit -m "docs: #465 multi-rev leveraging releg capability note"`

---

## Out of scope for this plan (separate issues)

- **The discovery CAMPAIGN** — actually relegging the `repeated-moon-*-sweep`
  skeletons at scale with the multi-rev backend and re-stamping the registry. This
  plan ships the CAPABILITY + golden; the campaign spends it (as #449 did).
- **3D / inclined leveraging** (the Amalthea-inclination half of the re-open keys)
  — a sibling capability; the chain here is coplanar (the `leveraging_leg` /
  `vilm` regime). Out of scope for #465.
- **Discrete-chain ToF optimisation / Pareto front** (Part-1 branch&bound) — the
  greedy descent suffices for the in-band capability proof; a Pareto front is a
  refinement, deferred.

---

## Definition of done

- `MultiRevLeveragingReleg` backend + `leveraging_chain` orchestrator build green,
  behind the existing `Releg` protocol with no driver rewrite.
- The chain's summed ΔV reproduces the sourced Eq.13 leveraging quadrature floor
  (`≥ floor`, `≤ floor + finite-chain penalty`); the Europa endgame golden brackets
  154/128 m/s.
- The Galilean positive control closes **IN-BAND** with the multi-rev backend,
  demonstrably cheaper than the 13.18 km/s single-DSM close — the #465 gate.
- A Saturnian positive control also closes in-band.
- The Uranian disjoint case STILL reports empty (prefilter), re-stamped with the
  stronger multi-rev powered method.
- The zero-retarget limit reproduces `BallisticReleg`; the other backends are
  untouched (the new `arrival_moon` kwarg is a backwards-compatible protocol add).
- All `tests/data` + `tests/search` ratchets green; capability note in OUTSTANDING.
