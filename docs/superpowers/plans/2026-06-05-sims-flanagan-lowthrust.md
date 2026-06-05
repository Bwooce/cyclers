# Sims-Flanagan low-thrust leg model (Phased Plan)

> **For agentic workers:** This is a milestone roadmap (like the `docs/phases/mN`
> plans and `2026-06-03-the-forge-pipeline.md`). Each phase is executed as its own
> superpowers:subagent-driven-development run with a detailed task breakdown written
> at execution time. Work on `main` — do NOT branch. Golden-discipline throughout
> (no fabricated literature numbers; see Phase 5 on the golden-test problem).

**Goal:** Add a low-thrust trajectory capability to cyclerfinder via the
Sims-Flanagan transcription, so the project can model and (eventually) optimise
`trajectory_regime: powered` / low-thrust cyclers — the long-deferred v2 scope.
This plan opens that scope; it does **not** declare it complete.

**Reference:** Yam, Di Lorenzo & Izzo, "Constrained global optimization of
low-thrust interplanetary trajectories," IEEE CEC 2010,
DOI [10.1109/cec.2010.5586019](https://doi.org/10.1109/cec.2010.5586019). The
project read it in full 2026-06-01; what was extracted is recorded in
`docs/v2-future-references.md` §1. **The PDF is not in the local mirror** — work
only from those in-repo notes; never fabricate values from memory of the paper.

**Scope discipline.** The v1 scope (spec §2) is ballistic patched-conic +
impulsive flybys. Low-thrust is explicitly a *stretch / v2* goal (spec §2
"low-thrust is a stretch goal"; spec §3 stretch list "GMAT bridge & low-thrust").
This plan stands up the leg-model substrate first and keeps it fully decoupled
from the v1 ballistic path until Phase 4, so nothing v1 regresses.

**Module placement note (deviation from task brief).** The task brief names the
target directory `src/cyclerfinder/compute/`. That directory does not exist; the
existing compute primitives this work extends — `lambert.py`, `kepler.py`,
`flyby.py` — all live in `src/cyclerfinder/core/`. To match the real repo layout
and import conventions, the Sims-Flanagan module lands at
`src/cyclerfinder/core/sims_flanagan.py` with tests at
`tests/core/test_sims_flanagan.py`. This is the only deviation from the brief.

---

## The Sims-Flanagan transcription (what we are building)

From the Yam 2010 framing recorded in `docs/v2-future-references.md` §1: a
trajectory is divided into **legs** that begin and end with a planet. The
continuous low-thrust arc on each leg is modelled as a sequence of **impulsive
ΔV manoeuvres connected by conic (ballistic Kepler) arcs** — one ΔV per segment,
applied at the segment midpoint. Each leg is propagated **forward from its start
state** and **backward from its end state** to a **match point** (usually the
halfway point in time); the position/velocity mismatch `S_mf − S_mb` (the
**defect**) must fall below tolerance for the leg to be feasible.

Per-segment thrust bound (Yam Eq. 1, as recorded):
`ΔV_max = (T_max / m) · (t_f − t_0) / N` for `N` segments. We apply the bound
**per segment** with the instantaneous segment mass `m_i` (a tighter, physically
honest variant of the constant-mass Eq. 1 form). Mass updates segment to segment
via the rocket equation (Yam Eq. 5, as recorded):
`m_{i+1} = m_i · exp(−ΔV_i / (g0·Isp))`.

---

## Phase 1 — Sims-Flanagan leg model (THIS COMMIT)

**What.** A pure, typed `core/sims_flanagan.py` implementing the leg model only —
no optimiser coupling.

- A frozen `SimsFlanaganLeg` config dataclass: start state `(r0, v0)`, end state
  `(rf, vf)`, total time of flight, segment count `N`, spacecraft mass at leg
  start `m0`, propulsion params `T_max` (thrust), `Isp`, central-body `mu`, and
  the match-point segment index (defaults to the temporal midpoint).
- The per-segment ΔV vectors as the manoeuvre decision variables (a length-`N`
  array of 3-vectors; all-zero ⇒ pure ballistic / coast).
- `propagate_forward(...)`: from `(r0, v0, m0)`, for each segment up to the match
  point, coast a half-segment with `core.kepler.propagate`, apply the segment ΔV
  impulse (updating mass via Tsiolkovsky), coast the remaining half-segment.
  Returns the match-point state + mass.
- `propagate_backward(...)`: the mirror from `(rf, vf)` with negative `dt` coasts
  and ΔV applied with the time-reversed sign convention; mass is propagated
  *backward* (it increases going back in time, since forward burns reduce mass).
  Returns the match-point state + (back-propagated) mass.
- `match_point_defect(...)`: the 7-vector `[Δr (3), Δv (3), Δm (1)]` between the
  forward and backward match-point states (`neq = 7`, the 3-D + mass problem
  dimension recorded from Yam §1). At convergence this → 0.
- `segment_dv_bounds(...)`: the per-segment `ΔV_max` bound array used later as the
  optimiser inequality constraint and asserted as an invariant in tests.

**Conventions to match exactly** (per `core/kepler.py`, `core/flyby.py`,
`core/lambert.py`): `from __future__ import annotations`; `Vec3` alias
`NDArray[np.float64]`; module/`@dataclass(frozen=True)` numpydoc-style docstrings
with a citation block; explicit `mu` defaulting to `MU_SUN_KM3_S2`; pure
functions (no I/O, no globals); km / km/s / s / kg units throughout; `g0` and
the unit conversions sourced from `core/constants.py` (add a sourced
`STANDARD_GRAVITY_M_S2` there if absent — it is, so add it with a citation).

**Exit:** `ruff check` + `ruff format --check` + `mypy src` + the full fast
pytest suite green; new physics-invariant tests (Phase 5 list) green.

## Phase 2 — feasibility / defect constraints

**What.** Turn the raw defect vector into the constraint surface an optimiser
consumes.

- `leg_feasible(leg, dvs, tol) -> bool`: defect norm (position + velocity,
  appropriately scaled) below tolerance.
- A constraint-vector assembler over a *chain* of legs sharing flyby boundaries:
  the flyby turn-angle constraint reuses `core.flyby.max_bend` /
  `is_ballistic_feasible` (the same `sin(δ/2) = 1/(1 + r_p V∞²/μ)` bend formula
  already in `flyby.py`, matching Yam Eq. 3 as recorded).
- NLP-dimension bookkeeping `(8 + 3N)·M` and `neq·M` (Yam §1) so the assembled
  problem can be validated structurally against the recorded paper dimensions.

**Exit:** chain-level defect + flyby constraints assembled and unit-tested on
invariants (no literature numbers).

## Phase 3 — integration with the existing optimiser stack

**What.** Wire the leg model into the M5 optimiser pattern (`search/optimize.py`:
scipy `differential_evolution` global + SLSQP local polish), matching that file's
structure.

- A low-thrust mode whose decision vector is the per-segment ΔV (plus leg ToFs
  and boundary V∞), with the defect as an equality constraint and the per-segment
  `ΔV_max` + flyby bend as inequalities.
- Two-phase solve mirroring Yam §1: Phase-1 minimise total ΔV at constant mass;
  Phase-2 re-optimise to maximise final mass via the rocket equation. (Yam used
  SNOPT; we use SLSQP to match the existing stack — the *fidelity ladder /
  transcription* is what Yam validates, not the specific solver.)

**Exit:** optimiser converges the defect to tol on a self-consistent synthetic
leg; gated behind `@pytest.mark.slow` like the other optimiser tests.

## Phase 4 — low-thrust cycler-maintenance application (the v2 goal)

**What.** Use the Phase 1–3 substrate to model real low-thrust cycler
*maintenance* and emit `trajectory_regime: powered` catalogue rows.

- Replace the single-impulse powered-flyby surrogate (`flyby.flyby_dv`) for
  applicable entries with a distributed low-thrust maintenance arc, reporting a
  `propellant_mass_fraction` so low-thrust candidates are comparable to ballistic
  ones (the `model/score.py` extension flagged in `v2-future-references.md` §1).
- Schema: populate the spec §16.7 powered / `maneuvers[].type = flyby-powered`
  fields and `model_assumption` appropriately. Partition matching by
  `model_assumption` (spec §16.7 / M7 rule) so low-thrust rows never
  cross-compare with circular-coplanar literature.

**Exit:** at least one source-attested low-thrust cycler modelled end-to-end with
its propellant fraction; catalogue + validation gates green.

## Phase 5 — validation gates (and the golden-test problem)

**The golden-test problem (read before writing any test).** The project rule is
*golden/validation EXPECTED values must trace to a published source, never a
value our own code computed* (`docs/.../feedback_golden_tests_sourced_only.md`).
For Sims-Flanagan we have **no usable literature anchor yet**:

- Yam 2010's worked examples are **Jupiter rendezvous** transfers (arrival V∞
  constrained to zero), explicitly **out of scope** for cyclers and not the
  ballistic resonant geometry our schema encodes — recorded as such in
  `v2-future-references.md` §1. They are not cycler goldens.
- The Vasile & Campagnola JBIS tables (the only low-thrust MGA numbers we have
  transcribed — `docs/notes/2026-06-05-vasile-tables-retranscription.md`) are
  **transcription-blocked** as golden EXPECTED: they are a raster OCR of a
  font-broken PDF, which is not acceptable EXPECTED provenance even at 300 DPI.
  Usable as *reference/candidate* anchors only.

**Therefore Phases 1–2 tests rest on PHYSICS INVARIANTS, not literature numbers:**

1. **Zero-thrust reduces to Kepler** (the key, source-free, rigorous regression
   anchor): with all segment ΔV = 0, the forward propagation of a leg must equal
   `core.kepler.propagate(r0, v0, dt)` to tight tolerance (km / µm·s⁻¹ scale),
   and the forward+backward match-point states must coincide (zero defect).
2. **Zero-thrust = Lambert closure:** seed a leg's endpoints from a `core.lambert`
   solution; with zero ΔV the match-point defect is ~0 (the leg already closes
   ballistically).
3. **Mass monotonic non-increasing forward** (strictly decreasing on any nonzero
   burn), via Tsiolkovsky; and the backward mass propagation inverts it.
4. **ΔV ≤ thrust-capability bound:** every per-segment ΔV magnitude respects
   `ΔV_max = (T_max/m_i)·dt_seg`; `segment_dv_bounds` matches the applied bound.
5. **Match-point defect → 0 at convergence** (Phase 3, slow): the optimiser
   drives `‖S_mf − S_mb‖` below tol.
6. **Energy consistency:** on a coast (zero-ΔV) segment, specific orbital energy
   is conserved to `kepler.py`'s relative tolerance; each impulse changes energy
   by exactly the work-energy bookkeeping of the applied ΔV.
7. **Reversibility:** forward-then-backward over a full leg with a fixed ΔV
   schedule round-trips the start state.

**Trigger for adding literature gates.** Add literature-anchored golden tests
**only** when *either*: (a) a clean, non-raster JBIS-typeset copy of Vasile &
Campagnola is obtained (lifting the block in
`2026-06-05-vasile-tables-retranscription.md`), **or** (b) the Yam 2010 PDF (or
another source publishing a *cycler* — not rendezvous — low-thrust trajectory
with full state) is obtained and verified. Until then, invariants are the only
gates; this is intentional and documented, not a coverage gap to paper over with
fabricated numbers.

**Gates (every phase ends green):**
`uv run ruff check .` · `uv run ruff format --check .` · `uv run mypy src` ·
`uv run pytest tests/ -q -m "not slow"`.

---

## Execution deviation (Phase 4, task #37, 2026-06-05)

Phase 4 as written implies catalogue rows (`trajectory_regime: powered`) and
schema edits (spec §16.7 powered / `maneuvers[].type = flyby-powered` fields).
**This was deliberately scoped down to machinery-only.** No published source
supplies a powered low-thrust cycler row we hold with extractable numbers, so
fabricating a catalogue row is forbidden by the project's golden-test discipline
(and by the task brief's binding disciplines). Likewise `data/catalogue.yaml`,
`data/catalogue.schema.json`, `src/cyclerfinder/data/validate.py`, and
`docs/spec.md` were left untouched.

What landed instead: a powered-maintenance **evaluator**
(`search/lowthrust_maintenance.py`) that models a cycler's per-synodic
maintenance manoeuvre as a thrust-bounded low-thrust arc (reusing the Phase 1–3
Sims-Flanagan substrate) and reports a `propellant_mass_fraction`, so future
*sourced* powered rows (or the Forge) can consume it. It is demonstrated in
tests on the existing ballistic rows' geometry (the Aldrin E-M maintenance ΔV
recomputed under a thrust-bounded model and compared to the impulsive
`optimise_maintenance_dv` result as an internal-consistency check — **not** a
golden). A schema hook for `propellant_mass_fraction` / powered fields is a
follow-up for whenever a sourced powered row is obtained.

## Sequencing rule

Each phase ends green (full fast suite + ruff + mypy), commits locally per task,
golden-discipline throughout. Phases are independently valuable — Phase 1 (the
leg model + invariant tests) is shippable on its own and unblocks all later work.
