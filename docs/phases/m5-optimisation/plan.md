# M5 — Optimisation

**Spec reference:** spec.md §4 (architecture — `search/optimize.py`), §5 step 5 (optimise pipeline step), §6 (`find_cyclers` top-level interface sketch), §7 (tech stack — pygmo is *optional*), §8 (M5 milestone definition — "rediscover a published **2-synodic E–M ballistic cycler from scratch**. Gate: matches published low-V∞ values (~5.65 km/s Earth, ~3.05 km/s Mars) within tolerance."), §9 (validation anchors — the 5.65 / 3.05 km/s targets plus the V∞ > 11 km/s degenerate-solution guard), §10 (search-landscape brutality — "**never rely on blind DE over a guessed sequence**"), §12(a) (two optimisation modes — *idealized* strict closure vs *ephemeris* finite-horizon TCM), §12(d) (**hard inequality constraints**, not soft regularisers — supersedes the §5 step 5 original wording), §13.1 (the search atomic unit = a `Cell`), §13.4 (**making the inner timing search near-deterministic** — free-return / resonance construction + fixed multi-start grid + local polish), §13.7 (prioritisation), §13.8 (`reproducibility:` block in the catalogue record, populated from this module's `seed`).

**Purpose:** stand up the **continuous optimisation layer**. M0–M3 built the primitives; M4 produced the discrete enumerator + scoring. M5 closes the gap by **searching the timing DOF inside a single M4 `Cell`** to produce a fully-closed cycler with minimised maintenance ΔV under hard inequality constraints. Compose with M4 to deliver the top-level `find_cyclers(...)` interface from spec §6, i.e. the v1 discovery pipeline end-to-end (Tisserand-pruned enumerator → constructor → per-cell optimiser → ranker). M5 is the **conceptually hardest milestone in the project** (per spec §10's risk list); the plan is sized accordingly.

**Gate (definition of done):**
1. `tests/search/test_optimize.py::test_2syn_em_rediscovers_5_65_kms_earth` asserts the published 2-synodic E-M-E ballistic cycler is **rediscovered from scratch** — i.e. starting from the M4 cell `("E","M","E"), period_k=2, per_leg_revs=(0,0), per_leg_branch=("single","single")` plus `vinf_cap = 7.0`, the optimiser converges to a closed cycler whose Earth-encounter V∞ magnitude is within ±0.2 km/s of **5.65 km/s** and whose Mars-encounter V∞ magnitude is within ±0.2 km/s of **3.05 km/s**, with `hard_constraints_pass=True`. The known-answer V∞ values are the **assertion target**, never an input to the optimiser (no catalogue look-up, no signature short-circuit).
2. `tests/search/test_optimize.py::test_2syn_em_rejects_high_vinf_degenerate` asserts the spec §9 / §10 degenerate-solution guard: with `vinf_cap = 7.0`, the optimiser never returns a "closure" relying on V∞ > 11 km/s. Implementation strategy: feed the optimiser a deliberately bad initial guess that would slide into the high-V∞ basin under an unconstrained run; assert the constraint formulation rejects it (`constraints_satisfied=False` if the cap is breached; the SLSQP polish refuses to exit the V∞ ≤ cap feasible region).
3. `tests/search/test_optimize.py::test_aldrin_regression_anchor` asserts the Aldrin geometry (M4 hand-off: `composite_score = 4.239371`, `max_vinf = 9.743359`, `taxi_cost = 6.530070`) survives an `optimise_cell_idealized` round-trip without score regression beyond ≤ 1e-3 relative tolerance. Catches accidental optimiser-driven regressions of the M4 scoring layer.
4. `tests/search/test_optimize.py::test_find_cyclers_em_top_level` asserts the spec §6 top-level `find_cyclers(bodies=("E","M"), k_synodic=2, vinf_cap=7.0, n_keep=5)` returns a non-empty ranked list whose top entry's cycler matches the same 5.65 / 3.05 V∞ signature (within ±0.2 km/s) as the per-cell gate, with the published cell at or near rank 1.
5. `tests/search/test_optimize.py::test_optimisation_result_frozen_and_seeded` asserts `OptimisationResult` raises `FrozenInstanceError` on assignment, and that two runs with the same `seed` and same `Cell` produce bitwise-identical `best_score` floats. Pins reproducibility per spec §13.8.
6. `tests/search/test_optimize.py::test_ephemeris_mode_stubbed_until_m6` asserts that calling `optimise_cell_ephemeris(...)` from M5 raises `NotImplementedError("requires M6 ephemeris backend")`. The API surface is present; the working implementation arrives in M6b.
7. `uv run pytest`, `uv run ruff check .`, `uv run ruff format --check .`, `uv run mypy src tests` all clean. CI green on the M5 commit.

The §8 published-V∞ values are tolerance ±0.2 km/s rather than tighter because: (a) the literature reports the figures to 0.05 km/s, (b) the published cycler is from a real-ephemeris constrained search (McConaghy et al. 2006) whereas M5 closes in the circular-coplanar model — the small mismatch is the bridge §12(a) calls out — and (c) the discrete `per_leg_revs` / `per_leg_branch` choices in the M4 cell may differ slightly from the literature parameterisation. ±0.2 is loose enough to absorb these legitimate differences and tight enough to *exclude* the degenerate high-V∞ basin (~11 km/s) by an order of magnitude.

---

## 1. What this milestone delivers

One new source module and its test file. M5 is **purely additive** at the module level; no edits to any M0–M4 source file.

### 1.1 `src/cyclerfinder/search/optimize.py` — the inner-timing optimiser

Public surface, in dependency order:

- `OptimisationResult` — frozen dataclass carrying the optimised cycler, its `Score`, closure residual, convergence flags, and the multi-start history (for diagnostics + reproducibility).
- `optimise_cell_idealized(cell, ephem, *, vinf_cap, rp_factors=None, n_starts=5, seed=0, target_period_sec=None) -> OptimisationResult` — spec §12(a) *idealized mode*: strict periodic closure over one period T in the circular-coplanar model. The discovery-mode optimiser; this is what the M5 gate exercises.
- `optimise_cell_ephemeris(cell, ephem, *, vinf_cap, n_laps=5, ...) -> OptimisationResult` — spec §12(a) *ephemeris mode*: minimise summed TCM ΔV over a finite 3–5 lap horizon on a real ephemeris. **API present, body raises `NotImplementedError`** — the working implementation requires M6's astropy backend + the multi-lap propagator from `verify/propagate.py`. M5 ships the signature so M6/M7 don't need to reshape the public surface.
- `find_cyclers(bodies, k_synodic, vinf_cap, *, n_keep=20, ephem=None, L_max=4, N_max=0, branch_set=("single",), n_starts=5, seed=0) -> list[OptimisationResult]` — the **complete v1 discovery pipeline** per spec §6: enumerate Tisserand-pruned cells from M4, run `optimise_cell_idealized` on each, filter by `hard_constraints_pass`, rank by `composite_score`, return the top N.

Module-internal helpers (private):

- `_free_return_seed(cell, ephem, target_period_sec)` — produces the structured initial guess for one start of the multi-start grid per spec §13.4 (free-return / resonance construction fixes most timing parameters).
- `_multi_start_grid(cell, ephem, n_starts, seed, target_period_sec)` — yields `n_starts` initial guesses by perturbing the free-return seed on a fixed reproducible grid.
- `_objective(x, cell, ephem, target_period_sec)` — closure-residual + sum-of-flyby-ΔV objective evaluated on a flat parameter vector `x` of free encounter offsets.
- `_constraints(x, cell, ephem, vinf_cap, rp_factors)` — per-encounter `V∞ ≤ vinf_cap` and `r_p ≥ r_p_min` inequality constraints in SLSQP's `{type: 'ineq', fun: ...}` format.
- `_polish(x0, cell, ephem, vinf_cap, rp_factors, target_period_sec)` — single SLSQP local solve.
- `_de_pass(cell, ephem, vinf_cap, rp_factors, seed, target_period_sec)` — single scipy `differential_evolution` global pass over the bounded parameter box; the result is fed into `_polish`. Optional (controlled by `use_de=True` kwarg on the public function); the structured multi-start grid carries most of the global coverage per spec §13.4, with DE as a defence-in-depth wrapper.

### 1.2 Test file

`tests/search/test_optimize.py` — the gate tests above plus the supporting unit tests on the helpers (objective sign, constraint formulation, multi-start determinism, ephemeris-mode stub).

### 1.3 Explicit non-goals (M5 boundaries)

These belong to adjacent milestones; **do not stub or partially implement them in M5** beyond the `optimise_cell_ephemeris` signature itself:

| Out of M5 | Where it lands | Why deferred |
|---|---|---|
| Real ephemeris backend (astropy, JPL DE) | **M6a** | M5 closes in `Ephemeris("circular")` only. The `optimise_cell_ephemeris` mode's *body* needs M6's astropy backend + multi-lap propagator; the *signature* lives in M5 so M6 only fills the body. |
| Phase-matching to real launch windows (`search/phase_match.py`) | **M6b** | Spec §12.1 — the idealized→ephemeris bridge is its own module, separate from the inner-timing optimiser. |
| Catalogue / canonical-signature matching, ledger writes | **M7** | M5's gate test asserts the V∞ values **independently**, never via a catalogue look-up. The signature-distance matcher (spec §16.3) is M7's concern. |
| CLI / viz | **M8** | M5's output is a `list[OptimisationResult]`, not a serialised artefact. |
| pygmo as the global optimiser | **Stretch (post-M8)** | Decided below (§3.1 → §6): commit to scipy DE; pygmo deferred. |
| The §13.6/§13.8 work-queue + ledger + parallel execution | **M7** | M4 deferred this; M5 inherits the deferral. `find_cyclers` is single-process. |
| Optimising the discrete cell structure (sequence, revs, branches) | **N/A — design intent** | Per spec §13.1 the discrete structure is the M4 enumerator's job. M5 optimises *within* a fixed cell. Conflating the two is the spec §10 failure mode. |
| Powered flyby ΔV minimisation as an objective for closure | **N/A — design intent** | The objective minimises **maintenance** ΔV (the sum of `flyby_dv` at all encounters); a feasible closure has this at or near zero per M4's hard-constraints semantics. The cycler closes when the residual vanishes; ΔV is what we minimise on the path there. |

---

## 2. File tree after M5

```
cyclers/
├── … (M0/M1/M2/M3/M4 layout preserved unchanged)
├── src/cyclerfinder/
│   ├── core/                          # unchanged
│   ├── search/
│   │   ├── __init__.py                # unchanged (or +1 re-export, see §3.6)
│   │   ├── tisserand.py               # M2 — unchanged
│   │   ├── resonance.py               # M2 — unchanged
│   │   ├── sequence.py                # M4 — unchanged
│   │   ├── construct.py               # M3 — unchanged
│   │   └── optimize.py                # NEW (M5)
│   └── model/                         # unchanged
└── tests/
    ├── … (M0/M1/M2/M3/M4 tests preserved)
    └── search/
        ├── __init__.py                # unchanged (M4 created)
        ├── test_sequence.py           # M4 — unchanged
        └── test_optimize.py           # NEW (M5 — includes the gate tests)
```

Subpackages `verify/`, `data/`, `viz/` still remain uncreated — M6/M7/M8 territory. No edits to any M0–M4 source. The only optional touch is a re-export of the M5 public names in `src/cyclerfinder/search/__init__.py` (currently empty); see §3.6.

---

## 3. Module design

This section is rich because the entire module is one design decision repeated at five levels of granularity (mode → algorithm → constraints → result → wrapper). The spec §10 failure mode — blind DE on a guessed sequence at ~26 km/s — is what each subsection is built to avoid.

### 3.1 Two optimisation modes (spec §12(a)) — the binding choice

Spec §12(a) supersedes the original §5 step 5 wording. There are two distinct optimisation targets, and they are not interchangeable:

| Mode | Objective | When used | Backend status in M5 |
|---|---|---|---|
| **Idealized** | Strict periodic closure over one period `T` in the circular-coplanar model. Minimise `closure_residual(cycler, omega_target) + Σ flyby_dv(enc)`, both in km/s. Goal: discover the geometry. | Discovery — the M5 gate. | **Implemented in M5**. |
| **Ephemeris** | Minimise summed TCM ΔV over a finite horizon of 3–5 laps (~20–30 yr) on the real astropy/JPL ephemeris. Closure is *not* exactly required; bounded lap-to-lap drift is. | Realising a candidate to a launch window. | **API present, body `NotImplementedError("requires M6 ephemeris backend")`** until M6b. |

**Binding decisions baked into M5:**

- The two modes are **separate functions**, not a `mode="…"` parameter on one function. The argument lists are genuinely different (`n_laps` only applies to ephemeris mode; `target_period_sec` only applies to idealized mode), and a future reader is better served by two short docstrings than one long one with conditional sections.
- Idealized mode is the one that gates M5. The published 2-synodic E-M cycler's V∞ values *are* a circular-coplanar quantity at this fidelity — the small real-ephemeris differences (eccentricity, inclination) are what §12(a) hands to ephemeris mode to absorb via TCMs and b-plane steering. M5's tolerance (±0.2 km/s) is loose enough to absorb this circular-vs-real gap.
- The ephemeris-mode signature is locked **in M5** so that M6 fills only the body, not the public API. This preserves caller code (`find_cyclers` doesn't change between M5 and M6; the realisation pipeline calls `optimise_cell_ephemeris` directly).

#### 3.1.1 Why scipy, not pygmo — committing the deferred decision

Per `docs/overview.md` §2 the global-optimiser choice was deferred to M5. **Decision: `scipy.optimize.differential_evolution` (DE) as the wrapping global pass; pygmo deferred to stretch.**

Rationale:

1. **Already a dependency.** `scipy>=1.13` is in `pyproject.toml` (M1) and consumed throughout M1–M4. Adding `pygmo` (compiled C++ via pybind11, plus the optional `pagmo2` archipelago dependency) is a substantial install-side burden — Python-3.11-only wheels, larger CI download — for a tool that the §13.4 structured-grid strategy makes mostly redundant.
2. **The §13.4 structured grid carries most of the global coverage.** The spec deliberately specifies "free-return / resonance construction so most timing parameters are fixed by the cell's structure, leaving a **low-dimensional root-find / bounded grid**." Once the multi-start grid is in place, the *global* optimiser's job is small: cover the residual continuous DOF the grid missed and break ties between near-equal basins. A 5-start grid + DE wrapper + SLSQP polish per start is the spec-honest topology; pygmo's specialised global solvers (CMAES, basin hopping, NSGA-II) are overkill for the residual.
3. **pygmo's strength — archipelago parallelism — is a property of the runner, not the optimiser.** When the M7 ledger + parallel work-queue land, parallelism is achieved by running *cells* in parallel across workers (spec §13.6), not by running one cell's optimisation across an archipelago. The per-cell budget is small enough to fit in a single SLSQP polish.
4. **scipy DE's constraint support is sufficient.** scipy ≥ 1.7 supports nonlinear constraints in `differential_evolution` via the `constraints=` argument. Per spec §12(d) we want hard inequalities; SLSQP (the local polish) supports them natively; DE supports them via penalty under the hood. The constraint formulation is the same shape for both.
5. **Reversibility.** If a future user wants pygmo (e.g. for the VEM campaign at M8 where the parameter space is genuinely higher-dimensional), the optimiser is a single module — swap the backend, keep the contract. Documented in §5 risk #6.

Pygmo is therefore **deferred, not rejected** — it stays a spec §7 *optional* dependency. M5 ships scipy-only.

#### 3.1.2 Why SLSQP for local polish

SLSQP (Sequential Least-Squares Quadratic Programming) is scipy's only constrained-NLP solver that handles **nonlinear inequality constraints** without a barrier or penalty reformulation. The alternatives:

- **Nelder-Mead**: unconstrained. Would require a soft penalty for V∞ caps — exactly what spec §12(d) prohibits.
- **L-BFGS-B**: bound constraints only (per-parameter box). The V∞ ≤ cap and r_p ≥ r_p_min constraints are nonlinear functions of the parameter vector, not parameter bounds.
- **trust-constr**: supports nonlinear constraints; heavier than needed for our 3–6 parameter problems and slower in practice. Documented as a fallback if SLSQP misbehaves.
- **COBYLA**: derivative-free constrained. Slower convergence than SLSQP on smooth problems like ours. Documented as a fallback.

SLSQP is the right size of tool; trust-constr / COBYLA are the documented escape hatches.

### 3.2 Inner-search algorithm (spec §13.4)

The spec §13.4 prescription is binding:

> Within a surviving cell, do **not** rely on a single blind optimiser run (the proven failure mode). Instead:
> - Use free-return / resonance construction so most timing parameters are fixed by the cell's structure, leaving a **low-dimensional root-find / bounded grid**.
> - Cover the remaining continuous DOF with a **fixed, reproducible multi-start grid** plus local polish, so coverage within a cell is systematic, not stochastic.

M5 implements this exactly:

```
optimise_cell_idealized(cell, ephem, *, vinf_cap, n_starts=5, seed=0, ...):
    target_period_sec = cell.period_k * synodic_period_days(...) * SECONDS_PER_DAY
    seeds = _multi_start_grid(cell, ephem, n_starts, seed, target_period_sec)
    polish_results = []
    for x0 in seeds:
        x_local = _polish(x0, cell, ephem, vinf_cap, rp_factors, target_period_sec)
        polish_results.append(x_local)
    # Optional global wrapper for defence-in-depth:
    if use_de:
        x_de = _de_pass(cell, ephem, vinf_cap, rp_factors, seed, target_period_sec)
        x_de_polished = _polish(x_de, cell, ephem, vinf_cap, rp_factors, target_period_sec)
        polish_results.append(x_de_polished)
    best = min(polish_results, key=lambda r: r.composite_with_constraints)
    cyc = _build_cycler_from_x(best.x, cell, ephem, target_period_sec)
    s = score(cyc, ephem, vinf_cap, target_period_sec=target_period_sec, rp_factors=rp_factors)
    residual = cyc.closure_residual(omega_rad_per_s=2*pi/target_period_sec)
    return OptimisationResult(cell=cell, best_cycler=cyc, best_score=s,
                              closure_residual_kms=residual,
                              optimiser_history=tuple(polish_results),
                              converged=best.success,
                              constraints_satisfied=s.hard_constraints_pass)
```

#### 3.2.1 Free parameter vector `x`

The continuous DOF per spec §13.4 are the encounter epochs offset from the structural anchor (the period `T` and the equispaced encounter spacing implied by the resonance). For a cell of `N` encounters:

- Encounter 0 is pinned to `t = 0` (the cycler's origin is a free choice; pinning kills a translational symmetry that the optimiser otherwise wastes effort on).
- Encounter `N-1` is pinned to `t = T = period_k * T_syn` (the period is the cell's discrete structure, not a free DOF).
- The interior encounters `t_1 … t_{N-2}` are the free parameters — `N - 2` of them.

For the M5 gate's 2-synodic E-M-E cell (`N=3`), this is **one free parameter**: `t_1`, the Mars encounter time. The free-return seed places it at `t_1 = T/2` (symmetric); the multi-start grid perturbs around this.

For longer sequences (M8's VEM cells at `N=6`–`8`) the free-parameter dimension grows to 4–6 — still small enough that 5–10 starts + SLSQP cover it densely. The module documents this scaling.

#### 3.2.2 `_free_return_seed` — the spec §13.4 structural anchor

The free-return seed places interior encounters at equispaced fractions of `T`, then nudges them onto the resonance lattice if the cell's `period_k > 1`. For a 2-synodic E-M-E cell the seed is `t = [0, T/2, T]`. For a 3-synodic E-V-M-E-M-E cell the seed is `t = [0, T/5, 2T/5, 3T/5, 4T/5, T]`. The seed is intentionally simple; the multi-start grid does the actual coverage.

This step replaces what a blind global optimiser would have to *discover*: that interior encounters live near rational fractions of the period. Without this anchor, DE wastes its first 200 generations finding the resonance lattice; with it, the residual continuous DOF is the small detuning from the lattice.

#### 3.2.3 `_multi_start_grid` — fixed reproducible coverage

`n_starts=5` deterministic perturbations of the free-return seed, generated with `numpy.random.default_rng(seed)`:

- Start 0: the free-return seed exactly.
- Starts 1–(n_starts-1): each interior `t_i` perturbed by `± k * T / (4N)` for `k ∈ {1, 2, …}` chosen from a fixed table (not random) so two runs with the same `seed` produce bitwise-identical starts.

The seed only affects how the *table* is shuffled across starts (so `n_starts < table_size` still picks a varied subset). With a small `n_starts` the table covers ≥ ¼ of the period in each direction, which is wider than the basin of attraction empirically observed in the M3 Aldrin closure work.

#### 3.2.4 `_polish` — SLSQP local solve per start

Each start hands its `x0` to `scipy.optimize.minimize(method='SLSQP', constraints=…)`. The constraints are passed as a `list[dict]` per the SLSQP convention. The polish budget is `maxiter=200` (SLSQP usually converges in 30–80); the tolerance is `ftol=1e-6` on the objective.

#### 3.2.5 `_de_pass` — optional global wrapper

`scipy.optimize.differential_evolution(_objective, bounds, constraints=NonlinearConstraint(...), seed=seed, polish=False)` with `polish=False` because we want the explicit SLSQP polish to follow (DE's built-in polish is L-BFGS-B which doesn't honour our nonlinear constraints). Bounds: each free `t_i` is bounded to `(0, T)` strictly inside the period. Constraints: same `_constraints` as SLSQP, wrapped in `NonlinearConstraint`.

The DE pass is **opt-in via `use_de=True` (default True)** because the cost is bounded (~5–10 seconds per cell on the M5 test machines) and the defence-in-depth value against missed basins is worth it. A user running M8's VEM enumeration over thousands of cells may pass `use_de=False` to skip it; the multi-start grid + polish is still spec-compliant on its own per §13.4.

#### 3.2.6 Selecting the "best" across starts

```python
def _composite_with_constraints(result):
    if not result.constraints_satisfied:
        return math.inf
    return result.objective_value
```

Identical pattern to M4's `composite_score`: infeasibles sort last by being `+inf`; feasibles sort by objective ascending; ties broken deterministically by start index (the table order). This pin makes "best of `n_starts`" reproducible across runs.

### 3.3 Hard-constraints formulation (spec §12(d)) — the second binding choice

Spec §12(d) supersedes the original §5 step 5's "low-V∞ regulariser" wording. The constraints are **hard inequalities**, formulated as:

```python
def _constraints(x, cell, ephem, vinf_cap, rp_factors):
    """Return list of SLSQP-style {'type':'ineq', 'fun': callable} dicts.

    SLSQP convention: fun(x) >= 0 means the constraint is satisfied.
    """
    cyc = _build_cycler_from_x(x, cell, ephem, target_period_sec)
    out = []
    for i, enc in enumerate(cyc.encounters):
        # V∞ cap: vinf_cap - max(||vinf_in||, ||vinf_out||) >= 0
        out.append({"type":"ineq", "fun": lambda x_, i=i:
                    vinf_cap - max(np.linalg.norm(_enc(x_, i).vinf_in),
                                   np.linalg.norm(_enc(x_, i).vinf_out))})
        # r_p floor: r_p_at(vinf, mu, delta) - rp_min(enc.body) >= 0
        rp_min = rp_factors.get(enc.body, 1.0) * SAFE_PERIHELION_KM[enc.body]
        out.append({"type":"ineq", "fun": lambda x_, i=i, rp=rp_min:
                    _r_p_at(_enc(x_, i)) - rp})
    return out
```

**Binding properties:**

- **Per-encounter, not aggregate.** Each flyby gets its own V∞ and r_p constraint. An aggregate "max V∞ ≤ cap" formulation is harder for SLSQP (non-smooth at the max) and less informative when it fails (don't know which encounter violated).
- **Inequality form `g(x) ≥ 0`.** SLSQP's native format. Slack at the bound is allowed; equality at the bound is fine.
- **No soft regulariser.** The objective is exclusively `closure_residual + Σ flyby_dv`; there is no `+ λ * (vinf_max)^2` term anywhere. Per spec §12(d) this is structural.
- **V∞ > 11 km/s degenerate-solution guard (§9) is structural.** Because `vinf_cap=7.0` by default in M5 testing, the V∞ ≤ vinf_cap constraint *already* rules out the degenerate basin. No separate "high-V∞ rejection" code path exists; the same constraint that bounds V∞ at 7 km/s also bounds it at 11 km/s. The §9 anchor is the *test target*; the §12(d) constraint formulation is the *mechanism*.

The `constraints_satisfied` field on `OptimisationResult` records whether SLSQP exited with all constraints satisfied. SLSQP's convergence flag distinguishes "converged but infeasible" (a constraint violation at the local min) from "converged and feasible" (all constraints satisfied) — both states are surfaced separately on the result.

### 3.4 `OptimisationResult` dataclass

```python
from __future__ import annotations
from dataclasses import dataclass

from cyclerfinder.search.sequence import Cell
from cyclerfinder.model import Cycler, Score


@dataclass(frozen=True)
class OptimisationResult:
    """The output of one cell's per-cell optimisation.

    All fields are immutable (frozen). The tuple-valued ``optimiser_history``
    is intentionally a ``tuple`` rather than a ``list`` so the frozen invariant
    is structural, not just by-convention.

    Spec references: §12(a) (two-mode optimisation produces a result of this
    shape), §12(d) (``constraints_satisfied`` records the hard-inequality
    outcome), §13.8 (``cell`` and the implicit ``seed`` populate the catalogue
    record's ``reproducibility`` block).
    """

    cell: Cell                                  # the M4 cell that was optimised
    best_cycler: Cycler                         # the closed/optimised trajectory
    best_score: Score                           # from model.score
    closure_residual_kms: float                 # rotating-frame closure ΔV, km/s
    optimiser_history: tuple[_StartRecord, ...]  # diagnostic per-start outcomes
    converged: bool                             # global success indicator
    constraints_satisfied: bool                 # all hard constraints met


@dataclass(frozen=True)
class _StartRecord:
    """Diagnostic record for one start in the multi-start grid.

    Module-private. Surfaced on ``OptimisationResult.optimiser_history`` for
    debug/inspection; not part of the public API contract beyond ``len()``
    and indexing.
    """

    start_index: int                # position in the multi-start grid
    x0: tuple[float, ...]           # initial parameter vector
    x_final: tuple[float, ...]      # converged parameter vector
    objective_value: float          # objective at x_final
    constraints_satisfied: bool     # at x_final
    nit: int                        # SLSQP iteration count
    success: bool                   # SLSQP convergence flag
```

**Design notes:**

- `cell` is carried so a downstream consumer (the M7 ledger; the M8 reporter) can populate `cell_id` directly from `result.cell.id` without re-deriving it.
- `best_cycler` is a `Cycler` from M3, so `result.best_cycler.maintenance_dv()` etc. all just work — M5 doesn't shadow `Cycler`'s API.
- `best_score` is a `Score` from M4; the consumer reads `taxi_cost_kms`, `max_vinf_kms`, etc. from it directly.
- `closure_residual_kms` duplicates information available via `best_cycler.closure_residual(omega)` but the M5 module computed it with a known `omega_target`, so caching the result avoids recomputing the rotating-frame transform in consumers.
- `optimiser_history` is a `tuple[_StartRecord, ...]`. The leading `_` marks `_StartRecord` as module-private — consumers should treat it as opaque diagnostic data.
- `converged ∧ constraints_satisfied` is the "trustworthy" predicate. Either alone is insufficient: `converged ∧ ¬constraints_satisfied` means SLSQP found a local min outside the feasible region (degenerate basin); `¬converged ∧ constraints_satisfied` means iteration cap hit but the current point is feasible (re-run with more `maxiter` may help).

### 3.5 `find_cyclers` — the spec §6 top-level wrapper

```python
def find_cyclers(
    bodies: tuple[str, ...],
    k_synodic: int,
    vinf_cap: float,
    *,
    n_keep: int = 20,
    ephem: Ephemeris | None = None,
    L_max: int = 4,
    N_max: int = 0,
    branch_set: tuple[str, ...] = ("single",),
    n_starts: int = 5,
    seed: int = 0,
    use_de: bool = True,
    rp_factors: dict[str, float] | None = None,
) -> list[OptimisationResult]:
    """Spec §6 top-level discovery interface.

    Pipeline (composes M4 + M5):
        1. ephem = ephem or Ephemeris("circular")
        2. cells = feasible_cells(bodies, L_max, k_synodic, N_max,
                                  vinf_cap, ephem, branch_set)
        3. results = [optimise_cell_idealized(c, ephem, vinf_cap=vinf_cap,
                                              n_starts=n_starts, seed=seed,
                                              use_de=use_de, rp_factors=rp_factors,
                                              target_period_sec=…)
                      for c in cells]
        4. feasible = [r for r in results if r.constraints_satisfied
                                          and r.best_score.hard_constraints_pass]
        5. feasible.sort(key=lambda r: composite_score(r.best_score))
        6. return feasible[:n_keep]

    Parameters mostly mirror :func:`optimise_cell_idealized`. ``L_max``,
    ``N_max``, ``branch_set`` control the M4 enumerator's cap; defaults match
    M4's gate test (L=4, N=0, single-branch only).

    ``k_synodic`` is interpreted as ``k_max`` for the M4 enumerator AND as
    ``period_k`` for the per-cell ``target_period_sec`` resolution — i.e.
    if you pass ``k_synodic=2`` you get cells with ``period_k ∈ {1, 2}``,
    and the target period for each cell is ``cell.period_k * T_syn``.

    Notes
    -----
    Single-process (per spec §13.6 deferred to M7). Single-pair only —
    the ``k_synodic`` interpretation as a single-pair multiplier is the
    M5/M8 contract; multi-body beats (the §3.4 VEM 6.4-yr case) require
    the §13.4 multi-pair resonance computation which M8's VEM campaign
    is the right place to add.
    """
```

**Why this lives in `optimize.py` and not its own module:**

- The composition is mechanical: enumerator → optimiser → rank. Splitting it into a separate `discover.py` would invert the dependency graph (the discover module would depend on optimize.py).
- The M5 gate test exercises `find_cyclers` directly; keeping it in the same module keeps the gate's import surface small.
- M8 will likely grow a `cli.py` that calls `find_cyclers`; that's the right layer to add a serialisation step.

### 3.6 Re-exports

`src/cyclerfinder/search/__init__.py` may optionally re-export the new public symbols (`OptimisationResult`, `optimise_cell_idealized`, `optimise_cell_ephemeris`, `find_cyclers`) so callers can write `from cyclerfinder.search import find_cyclers`. M4 did not re-export `Cell` etc.; M5 follows the same convention by default (callers import from the specific submodule). The re-export is a nice-to-have; tests use the fully-qualified path. Not load-bearing.

### 3.7 Imports / dependency graph after M5

```
constants.py            (M0)
ephemeris.py            (M1)
lambert.py              (M1)
kepler.py               (M1)
flyby.py                (M2)
frames.py               (M3)
tisserand.py            (M2)
resonance.py            (M2)
model/cycler.py         (M3)
search/construct.py     (M3) ← lambert, ephemeris, model/cycler
search/sequence.py      (M4) ← tisserand
model/score.py          (M4) ← flyby, model/cycler
search/optimize.py      (M5) ← sequence, construct, score, cycler,         [new]
                              ephemeris, resonance, flyby, constants,
                              scipy.optimize
```

No cycles. M5's optimiser depends on essentially everything below it — this is expected: it's the first module that wires the entire pipeline.

### 3.8 API summary

| Symbol | Purpose | Notes |
|---|---|---|
| `OptimisationResult` | Frozen dataclass: per-cell optimisation outcome | Consumers read `.best_score`, `.best_cycler`, `.constraints_satisfied` |
| `optimise_cell_idealized(cell, ephem, *, vinf_cap, …)` | Spec §12(a) idealized mode — strict closure | **M5 implemented** |
| `optimise_cell_ephemeris(cell, ephem, *, vinf_cap, n_laps=5, …)` | Spec §12(a) ephemeris mode — finite-horizon TCM | **M5 raises `NotImplementedError("requires M6 ephemeris backend")`** — signature locked |
| `find_cyclers(bodies, k_synodic, vinf_cap, …)` | Spec §6 top-level discovery pipeline | Composes M4 enumerator + M5 optimiser + M4 ranker |

---

## 4. Tests + gate

Tests live under `tests/search/test_optimize.py`. Tolerances are named at module level.

### 4.1 Gate tests (spec §8 binding)

| Test | Assertion | Tolerance |
|---|---|---|
| `test_2syn_em_rediscovers_5_65_kms_earth` (**gate**) | Cell `("E","M","E"), period_k=2, per_leg_revs=(0,0), per_leg_branch=("single","single")`, `vinf_cap=7.0` → result's `best_cycler.encounters` has an Earth encounter with `||vinf_in|| ≈ 5.65 km/s` and a Mars encounter with `||vinf|| ≈ 3.05 km/s`, AND `constraints_satisfied=True`, AND `closure_residual_kms < 0.5 km/s`. | V∞ ±0.2 km/s; residual hard < 0.5 km/s |
| `test_2syn_em_rejects_high_vinf_degenerate` (**gate**) | Same cell, but seed the optimiser with an `x0` deliberately near the high-V∞ basin (`t_1 = 0.05 * T`, putting Mars very close to Earth in time → demands V∞ ≈ 12 km/s for closure). Assert the polish returns `constraints_satisfied=False` for that start (the SLSQP polish cannot exit the cap-violating region) AND that the *overall result* picks a feasible start. | Boolean asserts |
| `test_aldrin_regression_anchor` | `cell = (("E","M"), ("E","M","E"), 1, (0,), ("single",))` (Aldrin's 1-synodic shape); run `optimise_cell_idealized`; assert `composite_score(result.best_score) ≤ 4.239371 + ε` (no regression worse than the M4 anchor by more than 1e-3 relative), `result.best_score.max_vinf_kms ≈ 9.74 ±0.1`, `taxi_cost_kms ≈ 6.53 ±0.1`. | rel ≤ 1e-3 on composite; ±0.1 km/s on metrics |
| `test_find_cyclers_em_top_level` (**gate**) | `results = find_cyclers(bodies=("E","M"), k_synodic=2, vinf_cap=7.0, n_keep=5, seed=0)`; `len(results) >= 1`; the top result's Earth/Mars V∞ match the M5 published anchors within ±0.2 km/s. | V∞ ±0.2 km/s |
| `test_optimisation_result_frozen_and_seeded` | `result.best_cycler = …` raises `FrozenInstanceError`; two calls with `seed=42` produce bitwise-identical `best_score.max_vinf_kms`. | Bitwise; exception |
| `test_ephemeris_mode_stubbed_until_m6` | `optimise_cell_ephemeris(cell, ephem, vinf_cap=7.0)` raises `NotImplementedError` with the documented message. | Exception type + message substring |

### 4.2 Unit tests (helpers)

| Test | Assertion |
|---|---|
| `test_objective_zero_at_closed_aldrin` | Hand-build the parameter vector for Aldrin's known geometry; `_objective(x, …)` returns a small finite number (objective is the closure residual + flyby ΔV sum; both are small for Aldrin in circular-coplanar). |
| `test_constraints_ineq_form` | `_constraints(x_aldrin, …)` returns a `list[dict]` where each `dict["fun"](x)` is a positive float (constraints satisfied). |
| `test_constraints_reject_at_cap` | A hand-built `x` that produces V∞ = 12 km/s → the corresponding V∞-cap constraint's `fun(x)` is negative. |
| `test_free_return_seed_em_2syn` | `_free_return_seed(cell, ephem, target_period_sec)` for the 2-syn E-M-E cell returns `t = (0, T/2, T)` exactly. |
| `test_multi_start_grid_deterministic` | Two calls with the same `seed` produce identical start vectors. |
| `test_multi_start_grid_distinct` | `n_starts=5` produces 5 distinct start vectors (no accidental duplicates). |
| `test_polish_returns_start_record` | `_polish(x0, …)` returns an object with `x_final`, `objective_value`, `constraints_satisfied`, `success`, `nit` populated. |
| `test_de_pass_obeys_seed` | Two calls to `_de_pass(…, seed=7)` produce bitwise-identical results. |
| `test_de_pass_can_be_disabled` | `optimise_cell_idealized(…, use_de=False)` runs faster than `use_de=True` and still satisfies the M5 gate's V∞ tolerance. (Soft test — runtime measurement is captured as a `pytest-benchmark` aside, not a hard assert.) |

### 4.3 Composition tests (`find_cyclers`)

| Test | Assertion |
|---|---|
| `test_find_cyclers_empty_when_caps_too_low` | `find_cyclers(("E","M"), k_synodic=1, vinf_cap=1.0)` → empty list (no feasible cell at vinf_cap=1.0). |
| `test_find_cyclers_n_keep_truncation` | `n_keep=2` ⇒ at most 2 results. |
| `test_find_cyclers_results_sorted` | Output is sorted ascending by `composite_score(r.best_score)`. |
| `test_find_cyclers_all_results_feasible` | Every `r` in output has `r.constraints_satisfied=True` and `r.best_score.hard_constraints_pass=True`. |

### 4.4 Tolerance summary

| Layer | Quantity | Tolerance |
|---|---|---|
| Earth V∞ at the rediscovered cycler | km/s | ±0.2 |
| Mars V∞ at the rediscovered cycler | km/s | ±0.2 |
| Closure residual at the rediscovered cycler | km/s | < 0.5 (hard) |
| Aldrin composite_score reproducibility | relative | ≤ 1e-3 |
| Aldrin `max_vinf_kms` reproducibility | km/s | ±0.1 |
| Reproducibility with `seed` | float | bitwise |

### 4.5 Why ±0.2 km/s (not ±0.05)

The published 5.65 / 3.05 km/s values are reported to 0.05 km/s precision in McConaghy et al. 2006, but for **real-ephemeris** trajectories with specific launch epochs. Our M5 circular-coplanar idealized closure cannot match real-ephemeris values to that precision — and shouldn't try; that's M6b's job. ±0.2 km/s is loose enough to:

- absorb the circular-coplanar → real-ephemeris gap (~0.1 km/s empirically),
- absorb the per-leg branch and revs choice (single-branch may differ slightly from the literature parameterisation),
- absorb numerical noise in SLSQP convergence,

and tight enough to:

- exclude the degenerate-solution basin (~11 km/s) by an order of magnitude,
- distinguish 2-synodic from 1-synodic geometries (the Aldrin V∞_E ≈ 6.53 is 0.88 away from 5.65 — well outside ±0.2),
- reject any "closure" relying on the trivial null trajectory (V∞ = 0 is 5.45 away from 5.65).

---

## 5. Risks

| # | Risk | Likelihood | Impact | Mitigation in M5 |
|---|---|---|---|---|
| 1 | **Search landscape is brutal — the M5 gate test fails to converge.** This is the spec §10 headline risk. The proven failure mode is blind DE on a guessed sequence at ~26 km/s. | medium | **high — milestone-blocking** | The §13.4 structural strategy is the entire defence: cell-fixed structure + free-return seed + multi-start grid + SLSQP polish + (optional) DE wrapper. The M5 gate test uses a known-good cell from M4's enumerator (not a hand-guessed one). If the test still fails, escalate: (a) raise `n_starts` from 5 to 25, (b) tighten the multi-start grid spacing near `T/2`, (c) add `trust-constr` as a fallback polish, (d) re-examine `_free_return_seed` for the 2-synodic case — the spec literature explicitly documents the half-period anchor; a wrong anchor here means the seed is in a different basin entirely. |
| 2 | **SLSQP returns "converged" but the solution is infeasible.** SLSQP's exit flag distinguishes converged-feasible from converged-infeasible, but a careless caller would conflate them. | medium | high — would let degenerate solutions through | `constraints_satisfied` and `converged` are **separate** fields on `OptimisationResult`. The gate test asserts `constraints_satisfied=True` explicitly, not `converged=True`. The "trustworthy" predicate is `converged ∧ constraints_satisfied`. The `_composite_with_constraints` helper that picks "best across starts" sends infeasibles to `+inf` so they never become the best result. |
| 3 | **Scipy DE wastes the time budget.** DE is slow on small problems where its diversity strategy doesn't help. | low | medium — slow tests | `use_de=True` by default but disablable via `use_de=False`. The gate test uses default; the longer M8 enumeration will likely pass `use_de=False`. The DE budget is capped at `maxiter=50, popsize=8` (small but enough to explore a 1–6 dim parameter box). Document that DE is the "defence-in-depth" pass; the multi-start grid + polish is the primary mechanism. |
| 4 | **Aldrin regression test breaks when the optimiser "improves" Aldrin.** The M4 anchor `composite_score = 4.239` is the M3 hand-built Aldrin seed's score; an optimiser run on Aldrin's cell may find a different (possibly better, possibly different) local min and the regression test fails. | medium | medium | The regression test asserts `composite_score ≤ M4_anchor + ε` (one-sided — improvement is OK) and `max_vinf ≈ 9.74 ±0.1` (the Aldrin family is identifiable by its high Mars V∞; even a re-optimised Aldrin lives near 9.74). If an optimiser run produces e.g. `max_vinf=5.65` for the Aldrin cell, that's diagnostic — it means the optimiser jumped to the **2-synodic basin** and the test should fail loudly so we investigate (likely a `period_k` confusion). |
| 5 | **`optimise_cell_ephemeris` stub leaks into `find_cyclers`.** A caller passing `ephem=Ephemeris("astropy")` (M6) might expect M5's `find_cyclers` to use ephemeris mode, since the parameter is there. | low | medium — confusing error | `find_cyclers` *always* calls `optimise_cell_idealized` in M5 regardless of `ephem`'s `model` attribute. Documented in the docstring. M6's `find_cyclers` extension (if any) will be a parameter `mode="idealized"` (default) or `mode="ephemeris"`, but the M5 contract is idealized-only. The `ephem` parameter exists in M5 because the inner construct() calls take an `Ephemeris`, which can still be circular-coplanar with the M6 backend installed. |
| 6 | **scipy DE proves inadequate at M8 scale (VEM campaign).** The VEM cells at 6–8 encounters with multi-rev branches may genuinely need a heavier global optimiser. | medium | low (in M5) | M5's choice is **reversible**: the global optimiser is a single private helper (`_de_pass`). Swapping it for pygmo (or anything else) is a one-module change with the same input/output contract. Documented in §3.1.1 as the reversibility argument. |
| 7 | **Multi-start determinism breaks under scipy version drift.** scipy's RNG conventions have shifted over versions; pinning `seed` may not be sufficient if `differential_evolution`'s internal sampler changes. | low | medium — bitwise tests fail | Pin `scipy>=1.13` (already done in `pyproject.toml`). The `seed` kwarg on `differential_evolution` is documented stable since scipy 1.9. The `_multi_start_grid` uses `numpy.random.default_rng(seed)` which is bitwise-stable across numpy 2.x. CI's lockfile (uv.lock) makes the version pin enforceable. |
| 8 | **`_build_cycler_from_x` is internally a wrapper around `construct_cycler` — failure modes leak.** The optimiser's objective evaluates `construct_cycler`, which raises `ValueError` on Lambert pathologies (negative TOF, NaN solutions). | medium | medium — optimiser crashes during DE/SLSQP | Wrap `construct_cycler` in a try/except inside `_objective` and `_constraints`: on `ValueError`, return a large finite penalty (1e6) for the objective and a large constraint violation. Documented in `_objective`'s docstring. Tests: parameterised `_objective` over pathological `x` values asserts the function never raises. |
| 9 | **Pinning encounter 0 to `t=0` kills a symmetry but creates a phase-non-uniqueness across cells.** Two cells that differ only by a cyclic rotation of `sequence` will produce different `t_i` values that describe the same physical cycler. | low (in M5) | low (in M5) | This is the spec §16.2 canonicalisation problem; M7 owns it. M5's contract is: a `Cell` in some rotation produces a `Cycler` in that rotation. M7 collapses rotations at signature time. The `find_cyclers` output may therefore contain near-duplicate cyclers that differ only by rotation — documented in the docstring. M7 dedupes at the catalogue layer. |
| 10 | **The optimiser's "best" composite (`_composite_with_constraints`) uses a different objective than `composite_score`.** The optimiser minimises closure-residual + flyby-ΔV; the ranker (M4 `composite_score`) sorts by a weighted sum that also includes `max_vinf`, `taxi_cost`, etc. So the optimiser's "best" start within a cell may not be the cell's best `composite_score` for ranking purposes. | low | low | Documented as intentional. The optimiser minimises a *physically motivated* objective (close the orbit, minimise propellant); the ranker scores *for downstream selection* using a broader notion. These are different jobs. In practice the difference is small because both formulations punish hard-constraint violators identically. The gate test asserts the M5 result has good `composite_score` AND good V∞ values — both layers agree on the published cycler. |
| 11 | **The 2-synodic E-M cell's `period_k=2` interpretation assumes single-pair semantics.** A future VEM cell will use the multi-body beat from `resonance.find_beats`, not `k * synodic_period(pair)`. | low (in M5) | medium (in M8) | M5's `target_period_sec` is computed via `synodic_period_days(bodies[0], bodies[1]) * cell.period_k * SECONDS_PER_DAY`. This is single-pair-correct (the M5 gate uses 2 bodies). When M8 extends to VEM (3 bodies), the helper needs to dispatch to `resonance.find_beats(...)` for `len(bodies) >= 3`. Documented in the `_target_period_sec(cell)` helper's docstring as M8 follow-up; M5's tests don't exercise it. |
| 12 | **Bend-feasibility constraint not exactly an inequality on `x`.** `is_ballistic_feasible(vin, vout, mu, rp_min)` returns a bool, not a smooth function — feeding it to SLSQP as a constraint will fail. | medium | high — constraint formulation invalid | The M5 constraint formulation uses the **r_p floor** (a smooth function of `vin, vout, delta`) rather than the boolean feasibility, exactly because SLSQP needs smooth constraints. The r_p formula: given V∞ and the required bend `delta = angle(vin, vout)`, solve `sin(delta/2) = 1 / (1 + r_p * V∞^2 / mu)` for `r_p`; the constraint is `r_p >= r_p_min`. This is smooth and monotone in V∞ and delta. `_r_p_required(vin, vout, mu)` is a new helper in `optimize.py` (not in `flyby.py` — it's an optimisation-only auxiliary). |

---

## 6. Dependency additions

**None.** M5 uses only `scipy.optimize.{minimize, differential_evolution, NonlinearConstraint}` (already available via the existing `scipy>=1.13` runtime dep), `numpy`, the standard library, and the in-house M0–M4 modules. No edits to `pyproject.toml`; no `uv.lock` regeneration needed.

**`pygmo` deferred to stretch.** Per spec §7 it's *optional*. Per §3.1.1 above scipy is the binding M5 choice; pygmo can be added later as a swap-in backend for `_de_pass` if M8's VEM campaign reveals a need. The M5 commit explicitly notes "pygmo deferred to stretch" so the deferred-decisions table in `docs/overview.md` §2 can be updated accordingly at M5 closeout.

---

## 7. Order of work

The `todo.md` mirrors this with checkboxes.

1. **Re-read predecessor docs.** Confirm M4 hand-off note's Aldrin regression anchor (composite=4.239, max_vinf=9.743, taxi_cost=6.530); confirm M4 `feasible_cells` / `Cell` API; confirm M3 `construct_cycler` signature and the open-sequence convention at the first/last encounters; confirm M2 `synodic_period_days` and `constants.SECONDS_PER_DAY` / `DAYS_PER_JULIAN_YEAR`.
2. **Write `search/optimize.py` skeleton.** Module docstring referencing spec §12(a), §12(d), §13.4. `OptimisationResult` + `_StartRecord` dataclasses. All four public functions with full docstrings and `NotImplementedError` bodies. Run `mypy --strict` — confirm the type signatures are clean before any logic lands.
3. **Implement helpers in dependency order.** Each helper gets a paired unit test before the next one starts.
   - `_target_period_sec(cell)` — single-pair `k * T_syn(pair)`; doc M8 multi-body follow-up.
   - `_free_return_seed(cell, ephem, T)` — equispaced anchor; **unit test against the 2-syn E-M-E case** asserting `t = (0, T/2, T)` exactly.
   - `_multi_start_grid(cell, ephem, n_starts, seed, T)` — deterministic perturbation table; **unit tests** for determinism (same seed → identical), distinctness (n_starts → n distinct vectors).
   - `_build_cycler_from_x(x, cell, ephem, T)` — wrap `construct_cycler` with the parameter-vector unpack; **unit test** that hand-built `x` for Aldrin reproduces the M3 `build_aldrin_seed` cycler within float tolerance.
   - `_r_p_required(vin, vout, mu)` — smooth r_p constraint helper; **unit test** that calling it on a known feasible flyby gives `r_p > rp_min`.
   - `_objective(x, cell, ephem, T)` — closure-residual + sum-of-flyby-ΔV; **unit test** `test_objective_zero_at_closed_aldrin`.
   - `_constraints(x, cell, ephem, vinf_cap, rp_factors)` — list-of-dicts; **unit tests** `test_constraints_ineq_form` (positive at feasible) and `test_constraints_reject_at_cap` (negative at infeasible).
   - `_polish(x0, cell, ephem, vinf_cap, rp_factors, T)` — SLSQP call returning a `_StartRecord`; **unit test** `test_polish_returns_start_record`.
   - `_de_pass(cell, ephem, vinf_cap, rp_factors, seed, T)` — DE call; **unit tests** for seed determinism and disablability.
4. **Implement `optimise_cell_idealized`.** Compose the helpers per the §3.2 pseudocode. Land:
   - `test_2syn_em_rediscovers_5_65_kms_earth` (**gate**) — the M5 binding gate. If it fails, escalate per risk #1.
   - `test_2syn_em_rejects_high_vinf_degenerate` (**gate**) — the §10 degenerate-solution guard test.
   - `test_aldrin_regression_anchor` — the M4 hand-off anchor reproducibility.
   - `test_optimisation_result_frozen_and_seeded`.
5. **Implement `optimise_cell_ephemeris` as a stub.** Body `raise NotImplementedError("requires M6 ephemeris backend")`. Land `test_ephemeris_mode_stubbed_until_m6`.
6. **Implement `find_cyclers`.** Compose M4 enumerator + per-cell optimiser + filter + rank. Land:
   - `test_find_cyclers_em_top_level` (**gate**) — the spec §6 interface working end-to-end.
   - `test_find_cyclers_empty_when_caps_too_low`.
   - `test_find_cyclers_n_keep_truncation`.
   - `test_find_cyclers_results_sorted`.
   - `test_find_cyclers_all_results_feasible`.
7. **Optional re-exports.** Add `from cyclerfinder.search.optimize import OptimisationResult, optimise_cell_idealized, optimise_cell_ephemeris, find_cyclers` to `src/cyclerfinder/search/__init__.py` if the import path becomes verbose in tests. Decided at implementation time; not load-bearing.
8. **Run the full local quality gate:** `uv run pytest`, `uv run ruff check .`, `uv run ruff format --check .`, `uv run mypy src tests`.
9. **Commit** as `m5: per-cell inner-timing optimiser + find_cyclers pipeline (rediscovers 2-syn E-M)`. Push; confirm CI green.
10. **Update `docs/overview.md`.** §2 deferred-decisions table: cross out "Global optimiser — scipy DE vs pygmo (decided in M5)" and replace with a row in the kept-decisions table noting "scipy.optimize.differential_evolution (M5) — pygmo deferred to stretch." §4 milestone table: M5 → completed; M6a → planned.
11. **Hand-off note** appended to `todo.md` under `## Hand-off to M6a` covering: the actual V∞ values reproduced (vs the 5.65 / 3.05 anchors), the closure residual reached, any cells in the L=4 / k=2 enumeration that failed to converge (and why), the empirical SLSQP iteration count / DE generation count for the 2-syn case, and any escalation actions taken from risk #1.

The order is "helpers → idealized mode → ephemeris stub → top-level wrapper" deliberately: each step depends only on what came before; the M5 gate (rediscovering 5.65 / 3.05) lands as soon as `optimise_cell_idealized` is in place; `find_cyclers` is the composition that proves the v1 pipeline end-to-end.

---

## 8. Exit checklist (the gate, restated)

Before declaring M5 done:

- [ ] `uv run pytest tests/search/test_optimize.py` green; all six gate tests in §4.1 pass at the documented tolerances.
- [ ] `uv run pytest` green overall (no regression of M0–M4 tests, including the M4 Aldrin anchor `composite_score=4.239371`).
- [ ] `uv run ruff check .` clean.
- [ ] `uv run ruff format --check .` clean.
- [ ] `uv run mypy src tests` clean under `strict = true` — including `OptimisationResult`, `_StartRecord`, the `Iterator`-or-`list` return types, and the scipy `OptimizeResult` interop (which mypy treats as `Any` — explicit casts at the boundary).
- [ ] CI green on the M5 commit.
- [ ] `docs/overview.md` updated: M5 status = `completed`; M6a row marked `planned`; deferred-decision "Global optimiser — scipy vs pygmo" moved to kept-decisions with the chosen value.
- [ ] `## Hand-off to M6a` section appended to `phases/m5-optimisation/todo.md` covering:
  - The exact V∞ values reproduced for the 2-syn E-M-E gate (vs the 5.65 / 3.05 anchors), and the gap (which M6a/M6b absorbs).
  - The closure residual reached at convergence (the budget for M6a's "verify periodic over ≥3 laps").
  - The empirical SLSQP iteration count and (if `use_de=True`) the DE generation count for the M5 gate — informs M8's compute-budget planning.
  - Any cells in `find_cyclers(("E","M"), k_synodic=2, vinf_cap=7.0, L_max=4, N_max=0)` that failed to converge (and why) — informs M6a's verification ordering.
  - Whether `_build_cycler_from_x`'s try/except trap fired during the test pass — informs M6a's expectation for ephemeris-mode robustness.
  - Whether `use_de=True` was needed or whether the multi-start grid + SLSQP alone hit the V∞ anchors — informs the M8 compute-budget default for VEM.
  - Any spec ambiguities resolved during M5 implementation that M6a should know about (notably the single-pair `period_k` semantics that M8 will need to extend for VEM).

(Writing the M6a plan doc is the first task of M6a, not an M5 exit criterion.)
