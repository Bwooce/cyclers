# M-ED — Real-ephemeris multi-arc cycler discovery (ballistic differential corrector)

> **For agentic workers:** REQUIRED SUB-SKILL: `superpowers:subagent-driven-development`
> or `superpowers:executing-plans`. Checkbox steps; strict TDD (write failing
> test → run **red** → minimal impl → run **green** → commit). Work on `main` —
> **do NOT branch** (project rule). **Do NOT commit without the user's review**
> when this plan is being authored under a docs-only mandate; the commit
> messages below are the messages to *use* when the implementation phase runs.
> uv-managed venv (no pip). Lint/type gate before **every** commit:
> `uv run ruff check .` · `uv run ruff format --check .` · `uv run mypy src tests`.
> Fast suite: `uv run pytest -m "not slow"`; the rediscovery/corrector geometry
> tests are `slow`.
>
> This plan is the task-level expansion of the **approved** design
> `docs/superpowers/specs/2026-06-05-m-ed-realeph-multiarc-discovery-design.md`
> (read it, **including its Approval section** — all six recommendations were
> accepted). Where this plan and the spec could diverge, **the spec + its
> Approval win.** This plan implements **Option B** (corrector-first ballistic
> MVP); the headline gate is the Jones AAS 17-577 VEM multisets; S1L1 stays a
> non-gating diagnostic; the TCM budget and Option C refinement are deferred.

---

## Goal

Make the catalogue's multi-arc rows (the Jones VEM members first; the broader
203 multi-arc rows by descriptor/fallback later) solvable as **closed N-arc
ballistic chains on the real DE440 ephemeris**, by generalising the proven
`scripts/correct_s1l1_twoarc.py` two-arc differential corrector to an arbitrary
closed encounter sequence and wiring it as a new, default-inert `mode="ballistic"`
on `optimise_cell_ephemeris`.

The genuinely new content (per spec §0) is **not** an optimiser — the
epoch-resolution / multi-rev-Lambert / body-agnostic chain-build infrastructure
is already built and general (`optimise_cell_ephemeris`,
`optimize.py:1309-1516`; `_resolve_t0_multi_seed`, `optimize.py:1207-1306`;
`_build_chain`, `maintain.py:303-350`). The new content is:

1. a **ballistic-closure solver mode** (V∞-magnitude continuity residuals driven
   to zero with `least_squares`, period pinned, bend feasibility post-hoc) —
   `search/correct.py`;
2. the **descriptor → topology + ToF-seed** parser — `search/descriptor.py`;
3. the **seeding ladder** (descriptor → sourced anchors → coplanar → epoch scan);
4. **honest validation gates** — the Jones VEM xfail flip as the headline gate,
   S1L1 as a non-gating diagnostic, and a real ballistic-closure gate replacing
   the stale `data/discover.py:99-106` V3 branch.

**Honesty boundary (binding, carried from spec §1 and §7 risk note).** The S1L1
prototype's *closed family floors Mars V∞ ≈ 6.4 km/s* (project memory
`project_s1l1_realeph_closure_blocker.md`; `data/OUTSTANDING.md:45`), which is a
*different* family than the Jones Mars-V∞ targets of 2.50/2.79 km/s (EMEVVE) and
~3.12/3.85 km/s (MEEVEM). **This plan does not assume the VEM gate converges.**
Phase 5 defines an explicit STOP/report branch (§Phase 5, Task 5.4) for the case
where the corrector cannot reach the sourced VEM V∞ multiset within tolerance:
the xfail is left `xfail` with a documented finding, and the plan ships Phases
1–4 (the corrector + descriptor + mode + ladder) as the deliverable, with the
VEM convergence recorded as open research — exactly as M8-Core handed it over.

---

## Architecture

### The N-arc ballistic model (spec §2.1)

Generalise the 3-variable prototype to a closed encounter sequence
`B0-B1-…-Bn` with `B0 == Bn` (`E-M-E-E` for S1L1; `E-M-E-V-V-E` for the Jones
EMEVVE outbound; `M-E-E-V-E-M` for the MEEVEM inbound). The chain is exactly the
`Cell.sequence` already loaded (`sequence.py:106-111`; `Cell` carries
`bodies, sequence, period_k, per_leg_revs, per_leg_branch, period_basis`).

- **Free variables:** `x = [t0, T_1, …, T_{n-1}]` — launch epoch (seconds since
  J2000) + per-leg ToFs (days). This is the *same shape* `_build_chain` already
  consumes (`maintain.py:311` "`x = [t0, tof_1, …, tof_{n-1}]`"), so the entire
  forward map is reused.
- **Period constraint by elimination (spec §2.1(a), recommended):** pin one
  "slack" leg ToF to `T − Σ(others)` exactly as the prototype pins `T_EE`
  (`correct_s1l1_twoarc.py:77`), where `T = _target_period_sec(cell)`
  (`optimize.py:1405`). Dimension `n−1`. (The constraint-residual form, §2.1(b),
  is implemented only if no single leg is a clean slack leg — not needed for the
  VEM members.)
- **Residuals (the ballistic-closure system):** for every *intermediate*
  encounter `Bi` (`1 ≤ i ≤ n−1`), `|V∞_in(Bi)| − |V∞_out(Bi)|` (the flyby
  conserves V∞ magnitude); plus the closure term `|V∞_in(Bn)| − |V∞_out(B0)|`.
  This is the prototype's residual set (`correct_s1l1_twoarc.py:102-106`)
  generalised. Solve with `scipy.optimize.least_squares` (the prototype uses
  `method="lm"`, `correct_s1l1_twoarc.py:125-133`) so over/under-determination is
  irrelevant.
- **Bend feasibility (post-hoc, never in the residual, spec §2.1):** at each
  flyby check `required_turn ≤ max_turn(V∞, r_p_safe)` using the prototype's
  `_max_bend_deg` / `_bend_deg` (`correct_s1l1_twoarc.py:109-120`). A
  ballistically-closed chain whose flyby exceeds `max_turn` is *powered*, not
  ballistic — surfaced honestly, never fitted toward.

### Why a corrector, not the shipped optimiser (spec §0 finding 2, §2.2)

The shipped `optimise_cell_ephemeris` objective is summed flyby turn-deficit ΔV
(`_maintenance_dv_chain`, `maintain.py:353-368`, via `_objective`,
`maintain.py:371-391`). Turn-deficit ΔV is ≈0 across a wide feasible plateau, so
minimising it does **not** pin the ballistic family — it drifts onto a degenerate
high-V∞ basin (`close_s1l1_realeph.py`; project memory
`project_s1l1_realeph_closure_blocker.md`). The prototype's V∞-continuity
*residual* pins the ballistic node conditions. M-ED keeps **both** modes:
`mode="maintenance"` (existing, byte-identical default) answers the M7 TCM
question; `mode="ballistic"` (new) answers the M-ED discovery question.

### Frame-freedom (spec §6.1 — why M-3D is not a dependency)

The ballistic node-closure residual is **V∞-magnitude continuity at nodes**,
which needs no rotating frame. The prototype reads full 3D DE440 states
(`ephem.state` returns 3D `r,v`) and Lambert is 3D-native; `V∞ = v_sc − v_planet`
is a genuine 3D vector and its magnitude + bend angle are correct in 3D without
any frame change. Per the Approval (Q3), **M-ED does not block on M-3D**. M-L
(multi-rev Lambert) **is** confirmed landed (Approval Q4; `lambert(..., max_revs=N)`
exists, `core/lambert.py:518,565-566,652-658`), and the prototype already calls
`lambert(max_revs=1/2)` (`correct_s1l1_twoarc.py:65-67`).

---

## Tech stack

Python 3.11, numpy, `scipy.optimize.least_squares`, pytest + pyyaml. uv-managed
venv. New production modules: `src/cyclerfinder/search/correct.py`,
`src/cyclerfinder/search/descriptor.py`. Extended:
`src/cyclerfinder/search/optimize.py` (`optimise_cell_ephemeris` `mode=` kwarg),
`src/cyclerfinder/data/discover.py` (V3 branch). New tests under `tests/search/`
and `tests/data/`; the existing `tests/test_vem_rediscovery.py` xfail is flipped
(or left xfail per the STOP branch). Dependencies of the corrector are limited to
`core/lambert`, `core/ephemeris`, `core/constants` (the prototype's imports) — it
stays catalogue-agnostic.

---

## Spec references (verified against live code, 2026-06-05)

- spec §0 finding 1 — `optimise_cell_ephemeris` is fully implemented, not a stub
  (`optimize.py:1309-1516`). **Verified.**
- spec §0 finding 2 — shipped objective is `_maintenance_dv_chain`
  (`maintain.py:353-368` via `_objective`, `maintain.py:371-391`). **Verified.**
- spec §0 finding 3 — `data/discover.py:99-106` still believes the stub raises
  (`# raises until M6b lands`, `discover.py:105`). **Verified verbatim.**
- spec §1 / prototype anatomy — `correct_s1l1_twoarc.py` free vars
  `[t0_off, T_EM, T_ME]` (`:127`), `T_EE` pinned (`:77`), residuals (`:102-106`),
  branch-pinned legs (`:65-67`), bend feasibility (`:109-120,151`), epoch×branch
  scan (`:162-181`). **Verified.**
- spec §2.1 forward map shape — `_build_chain(x=[t0, tof_1…])`
  (`maintain.py:303-350`). **Verified.**
- spec §3.1 / §16.7.7 — `free_return_arcs[]` descriptor fields `arc_type`
  (g/G generic, h/H half-rev, f/F full-rev), `tof_years` (g/h only),
  `resonance` (M:N, f/F only), `raw_descriptor` (uppercase = designated transit
  leg) — `docs/spec.md:994-1016`; 12 entries with explicit descriptors
  (`docs/spec.md:1018`). **Verified.**
- spec §3.2 / §5 headline gate — Jones VEM rows
  `jones-2017-vem-emevve-outbound` (`catalogue.yaml:2011`) and
  `jones-2017-vem-meevem-inbound` (`catalogue.yaml:2288`) carry full sourced V∞
  multisets + sourced transit ToFs (309/259, 268/223). **Verified.** Mining
  golden anchors at `docs/notes/2026-06-05-jones-aas17-577-vem-mining.md:393-416`.
- spec §5 xfail target — `tests/test_vem_rediscovery.py:220-258`
  `test_emeeve_idealized_optimiser_converges_feasible`, `strict=False`,
  "Flipped by M-ED". **Verified.**
- spec §5 S1L1 non-gate — `EXPECTED_SKIPS["s1l1-2syn-em-cpom"]`
  (`real_closure.py:222-236`); the M5 xfail is
  `tests/test_catalogue_rediscovery.py::test_2syn_em_rediscovers_5_65_kms_earth`
  (`test_catalogue_rediscovery.py:12-13`) and the discover end-to-end variant
  (`tests/data/test_discover.py:110-126`, both unverified-provenance, flagged
  2026-06-04). **Verified — note the spec wrote "M5 `test_2syn…`"; the live name
  is exactly that.**
- spec §4 / discover wiring — `discover(..., optimiser="ephemeris")` routes to
  `optimise_cell_ephemeris` (`discover.py:179-193`); signs into the
  `analytic-ephemeris` signature pool (`discover.py:218`). **Verified.**

---

## Phasing (independently shippable)

| Phase | Theme | Tasks | M-L dep | discover() dep |
|---|---|---|---|---|
| **1** | Corrector core (`search/correct.py`) | 1.0–1.5 (6) | yes (already landed) | no |
| **2** | Descriptor parser (`search/descriptor.py`) | 2.0–2.3 (4) | no | no |
| **3** | `mode="ballistic"` wiring on `optimise_cell_ephemeris` | 3.0–3.3 (4) | no | no |
| **4** | Seeding ladder (descriptor → anchor → coplanar → scan) | 4.0–4.3 (4) | no | no |
| **5** | Discover V3 gate + Jones VEM headline gate + S1L1 diagnostic + census | 5.0–5.5 (6) | no | yes |

Phases 1–4 are each shippable on their own (the corrector + parser + mode are
useful without `discover()` wiring). Phase 5 wires the consumption surface and
carries the **STOP/report branch** (Task 5.4) if the headline gate does not
converge. Total: **24 tasks.**

---

## Phase 1 — Corrector core (`search/correct.py`)

Lift the prototype's `_solve` (`correct_s1l1_twoarc.py:123-152`) to N arcs as a
pure, body/length-agnostic module. Public surface (spec §4):

```python
ballistic_correct(
    sequence: tuple[str, ...],
    per_leg_revs: tuple[int, ...],
    per_leg_branch: tuple[str, ...],
    t0_seed_sec: float,
    tof_seed_days: Sequence[float],
    period_sec: float,
    ephem: Ephemeris,
    *,
    vinf_cap: float,
    rp_factors: dict[str, float] | None = None,
    slack_leg: int | None = None,
) -> BallisticClosureResult
```

### Task 1.0 — `BallisticClosureResult` dataclass + per-leg V∞ extraction

**Files:** create `src/cyclerfinder/search/correct.py`; test
`tests/search/test_correct_result.py`.

#### Failing test — `tests/search/test_correct_result.py`

```python
"""M-ED Phase 1: BallisticClosureResult shape (plan Phase 1)."""
from __future__ import annotations

from cyclerfinder.search.correct import BallisticClosureResult


def test_result_fields_present() -> None:
    r = BallisticClosureResult(
        t0_sec=0.0,
        tof_days=(154.0, 379.0, 1027.0),
        max_residual_kms=0.04,
        vinf_per_encounter_kms=(5.6, 6.4, 5.6, 5.6),
        converged=True,
        bend_feasible=True,
    )
    assert r.converged is True
    assert r.bend_feasible is True
    assert r.max_residual_kms == 0.04
    assert len(r.tof_days) == 3
    assert len(r.vinf_per_encounter_kms) == 4


def test_constraints_satisfied_is_converged_and_feasible() -> None:
    """constraints_satisfied = converged AND bend_feasible AND vinf-cap met
    (spec §2.2)."""
    r = BallisticClosureResult(
        t0_sec=0.0,
        tof_days=(154.0,),
        max_residual_kms=0.04,
        vinf_per_encounter_kms=(5.6, 6.4),
        converged=True,
        bend_feasible=True,
        vinf_cap_ok=True,
    )
    assert r.constraints_satisfied is True
    r2 = BallisticClosureResult(
        t0_sec=0.0,
        tof_days=(154.0,),
        max_residual_kms=0.04,
        vinf_per_encounter_kms=(5.6, 6.4),
        converged=True,
        bend_feasible=False,
        vinf_cap_ok=True,
    )
    assert r2.constraints_satisfied is False
```

Run: `uv run pytest tests/search/test_correct_result.py -q` → **red** (no module).

#### Minimal impl — `search/correct.py` (dataclass only, this task)

```python
"""N-arc ballistic differential corrector on the real ephemeris (spec §2.1).

Generalises scripts/correct_s1l1_twoarc.py: free vars x = [t0, leg ToFs] with
one leg pinned by the sourced period; residuals = flyby V_inf-magnitude
continuity + periodicity closure, driven to zero with least_squares; bend
feasibility checked post-hoc, never in the residual. Pure: depends only on
core/lambert, core/ephemeris, core/constants.
"""
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass


@dataclass(frozen=True)
class BallisticClosureResult:
    t0_sec: float
    tof_days: tuple[float, ...]
    max_residual_kms: float
    vinf_per_encounter_kms: tuple[float, ...]
    converged: bool
    bend_feasible: bool
    vinf_cap_ok: bool = True

    @property
    def constraints_satisfied(self) -> bool:
        return self.converged and self.bend_feasible and self.vinf_cap_ok
```

Run → **green**. Then lint/type. Commit:

```
search/correct: BallisticClosureResult dataclass (M-ED Phase 1)
```

### Task 1.1 — N-arc leg builder + V∞ node dictionary

**Files:** `search/correct.py`; test `tests/search/test_correct_legs.py`.

Lift `_legs` / `_state_vinf` (`correct_s1l1_twoarc.py:55-93`) to arbitrary N.
The builder takes `x = [t0_sec, tof_1, …, tof_{slack-1}, tof_{slack+1}, …]` with
the slack leg eliminated, reconstructs the slack leg as
`period_days − Σ(free legs)`, walks the cumulative encounter epochs, solves each
leg's Lambert with its `(n_revs, branch)` (via a `_pick` helper lifted from
`correct_s1l1_twoarc.py:48-52`), and returns per-encounter `V∞_in` / `V∞_out`
3D vectors.

#### Failing test — `tests/search/test_correct_legs.py`

```python
"""M-ED Phase 1: N-arc leg/V_inf builder reproduces the prototype on S1L1.

NON-GOLDEN: the V_inf values here are OUR computation (spec §5 / project memory
golden-tests-sourced-only). This is a non-regression fixture for the SOLVER,
not a published-anchor assertion. The numbers are pinned from the live
scripts/correct_s1l1_twoarc.py prototype output, not from any source.
"""
from __future__ import annotations

import numpy as np

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.correct import _vinf_nodes


def test_s1l1_two_arc_nodes_have_in_out_per_encounter() -> None:
    ephem = Ephemeris("astropy")
    # S1L1 E-M-E-E: t0 ~2030-03-22, T_EM 154 d, T_ME 379 d, slack leg = E-E.
    # period = (1.4612 + 2.8096) yr (Russell 4.991gG arcs), days.
    period_days = (1.4612 + 2.8096) * 365.25
    seq = ("E", "M", "E", "E")
    t0_sec = (np.datetime64("2030-03-22T00:00:00") - np.datetime64("2000-01-01T12:00:00")) / np.timedelta64(1, "s")
    nodes = _vinf_nodes(
        sequence=seq,
        per_leg_revs=(0, 1, 2),
        per_leg_branch=("single", "single", "low"),
        t0_sec=float(t0_sec),
        free_tof_days=(154.0, 379.0),  # slack = E-E leg (index 2), eliminated
        slack_leg=2,
        period_days=period_days,
        ephem=ephem,
    )
    # One node per encounter; intermediates carry both in/out; ends carry the
    # closure pair (B0 out vs Bn in).
    assert set(nodes) >= {"m_in", "m_out", "e0", "e1_in", "e1_out", "e2_in"}
    for k in ("m_in", "m_out", "e0", "e1_in", "e1_out", "e2_in"):
        assert np.asarray(nodes[k]).shape == (3,)
```

Run → **red**.

#### Minimal impl

Add `_pick`, `_vinf_nodes` to `correct.py`. `_vinf_nodes` reconstructs the slack
leg ToF (`period_days − sum(free_tof_days)`), inserts it at `slack_leg`, walks
encounter epochs, solves each leg with `lambert(r_i, r_{i+1}, tof*DAY_S,
max_revs=per_leg_revs[i])`, `_pick`s the requested `(n_revs, branch)`, and
returns the `V∞ = v_sc − v_planet` dict keyed `b{i}_in` / `b{i}_out` plus the
prototype-compatible aliases for S1L1 (`m_in`, `m_out`, `e0`, `e1_in`, …). Use
the prototype's `vi(vsc, vpl)` (`correct_s1l1_twoarc.py:83-84`). On any Lambert
pathology raise nothing here — return is only called inside a guarded residual
(Task 1.2).

Run → **green**. Lint/type. Commit:

```
search/correct: N-arc leg builder + per-encounter V_inf nodes (M-ED Phase 1)
```

### Task 1.2 — residual function (V∞-continuity + closure)

**Files:** `search/correct.py`; test `tests/search/test_correct_residuals.py`.

Lift `_residuals` (`correct_s1l1_twoarc.py:96-106`) generalised: for each
intermediate encounter `|V∞_in| − |V∞_out|`; plus closure `|V∞_in(Bn)| −
|V∞_out(B0)|`. Guard Lambert pathologies → large finite penalty (the prototype's
`[1e3, 1e3, 1e3]`, `correct_s1l1_twoarc.py:100`).

#### Failing test

```python
"""M-ED Phase 1: ballistic-closure residual vector (plan Phase 1 Task 1.2)."""
from __future__ import annotations

from cyclerfinder.search.correct import _residual_vector


def test_residual_length_is_n_minus_one_intermediates_plus_closure() -> None:
    # E-M-E-E: encounters B0..B3, intermediates B1,B2 -> 2 residuals + 1 closure.
    fake = {
        "b1_in": (3.0, 0.0, 0.0), "b1_out": (3.0, 0.0, 0.0),
        "b2_in": (5.0, 0.0, 0.0), "b2_out": (5.0, 0.0, 0.0),
        "b3_in": (5.6, 0.0, 0.0), "b0_out": (5.6, 0.0, 0.0),
    }
    res = _residual_vector(fake, n_encounters=4)
    assert len(res) == 3
    assert max(abs(r) for r in res) < 1e-12  # perfectly continuous fixture
```

Run → **red**, then impl `_residual_vector(nodes, n_encounters)` reading the
`b{i}_in/out` keys, then the top-level `_residuals(x, ...)` that calls
`_vinf_nodes` inside the Lambert guard. → **green**. Commit:

```
search/correct: V_inf-continuity + closure residual vector (M-ED Phase 1)
```

### Task 1.3 — bend feasibility (post-hoc)

**Files:** `search/correct.py`; test `tests/search/test_correct_bend.py`.

Lift `_max_bend_deg` / `_bend_deg` verbatim (`correct_s1l1_twoarc.py:109-120`),
applying `rp_factors` to `safe_alt_km` if supplied (spec §2.1
`r_p_safe`). `bend_feasible` = at every intermediate flyby
`required_turn ≤ max_turn(V∞_in, body)`.

#### Failing test

```python
"""M-ED Phase 1: post-hoc bend feasibility (plan Phase 1 Task 1.3)."""
from __future__ import annotations

from cyclerfinder.search.correct import _bend_deg, _max_bend_deg


def test_max_bend_decreases_with_vinf() -> None:
    # Higher V_inf -> tighter max turn at the same body (Mars).
    assert _max_bend_deg(3.0, "M") > _max_bend_deg(8.0, "M")


def test_bend_zero_for_parallel_vectors() -> None:
    assert _bend_deg((1.0, 0.0, 0.0), (2.0, 0.0, 0.0)) == 0.0
```

Run → **red** → impl → **green**. Commit:

```
search/correct: post-hoc bend feasibility (M-ED Phase 1)
```

### Task 1.4 — `ballistic_correct` (the public N-arc solver)

**Files:** `search/correct.py`; test `tests/search/test_correct_s1l1.py` (slow).

Lift `_solve` (`correct_s1l1_twoarc.py:123-152`) to the public signature:
`least_squares(_residuals, x0=[t0_sec, *free_tof_days], method="lm", max_nfev=80,
xtol=1e-9, ftol=1e-9)`; converged iff `max_residual_kms < tol` (default `0.1`
km/s, the prototype's threshold `correct_s1l1_twoarc.py:169`); fill
`vinf_per_encounter_kms`, `bend_feasible`, `vinf_cap_ok` (`max V∞ ≤ vinf_cap`).

#### Failing test (slow — non-regression of the solver, NON-GOLDEN)

```python
"""M-ED Phase 1: ballistic_correct closes the S1L1 two-arc chain on DE440.

NON-GOLDEN non-regression fixture (spec §5, project memory): the asserted V_inf
floor (Mars ~6.4 km/s) is OUR prior computation, pinned from the live prototype,
NOT a published anchor. This guards the SOLVER against regression; it is NOT a
rediscovery of any sourced number. See project memory
project_s1l1_realeph_closure_blocker.md for why this family floors at ~6.4.
"""
from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.correct import ballistic_correct


@pytest.mark.slow
def test_ballistic_correct_closes_s1l1_two_arc() -> None:
    ephem = Ephemeris("astropy")
    period_days = (1.4612 + 2.8096) * 365.25
    t0_seed = (np.datetime64("2030-03-22T00:00:00") - np.datetime64("2000-01-01T12:00:00")) / np.timedelta64(1, "s")
    r = ballistic_correct(
        sequence=("E", "M", "E", "E"),
        per_leg_revs=(0, 1, 2),
        per_leg_branch=("single", "single", "low"),
        t0_seed_sec=float(t0_seed),
        tof_seed_days=(154.0, 379.0),  # the two free legs; E-E is slack
        period_sec=period_days * 86400.0,
        ephem=ephem,
        vinf_cap=9.0,
        slack_leg=2,
    )
    assert r.converged  # ballistic closure reached (residual < 0.1 km/s)
    # The closed S1L1 family floors Mars V_inf ~6.4 km/s (OUR computation, not
    # a sourced anchor). Assert the regime, not a sourced value.
    vinf_mars = r.vinf_per_encounter_kms[1]  # encounter index 1 = Mars
    assert vinf_mars > 5.5
```

Run → **red** → impl `ballistic_correct` → **green** (run with `-m slow`).
Commit:

```
search/correct: ballistic_correct N-arc solver (M-ED Phase 1)
```

### Task 1.5 — full Phase-1 suite + lint/type gate

- [ ] `uv run pytest tests/search/test_correct_*.py -q` (and `-m slow` for the
  S1L1 closure) all green.
- [ ] `uv run ruff check . && uv run ruff format --check . && uv run mypy src tests`.
- [ ] Confirm `correct.py` imports only `core/lambert`, `core/ephemeris`,
  `core/constants` (catalogue-agnostic, spec §4). Commit if any cleanup:

```
search/correct: Phase 1 lint/type gate (M-ED)
```

---

## Phase 2 — Descriptor parser (`search/descriptor.py`)

Russell `free_return_arcs[]` (spec §16.7.7, `docs/spec.md:994-1016`) → per-leg
`n_revs`/`branch` + asymmetric `tof_seed_days`. Isolated so the corrector stays
catalogue-agnostic (spec §4).

### Task 2.0 — arc-type → leg topology mapping

**Files:** create `src/cyclerfinder/search/descriptor.py`; test
`tests/search/test_descriptor_arctype.py`.

Mapping (sourced from `docs/spec.md:998-1000`):
`generic` (g/G) → `n_revs=0`, `branch="single"` (a generic return is a direct
E-E leg); `half-rev` (h/H) → `n_revs=0`, `branch="single"`; `full-rev` (f/F) →
`n_revs = M` from the `M:N` resonance, `branch="low"` (the resonant E-E loop is
multi-rev). The uppercase letter marks the designated transit leg
(`docs/spec.md:1001`) — recorded but not needed for topology.

#### Failing test

```python
"""M-ED Phase 2: arc_type -> leg topology (plan Phase 2; spec §16.7.7)."""
from __future__ import annotations

from cyclerfinder.search.descriptor import arc_to_leg_topology


def test_generic_arc_is_direct_single() -> None:
    assert arc_to_leg_topology("generic", resonance=None) == (0, "single")


def test_half_rev_arc_is_direct_single() -> None:
    assert arc_to_leg_topology("half-rev", resonance=None) == (0, "single")


def test_full_rev_arc_uses_resonance_revs() -> None:
    # "3:2" -> spacecraft does 3 revs (M:N, M = spacecraft revs, spec §16.7.7).
    assert arc_to_leg_topology("full-rev", resonance="3:2") == (3, "low")
```

Run → **red** → impl `arc_to_leg_topology(arc_type, *, resonance)` → **green**.
Commit:

```
search/descriptor: arc_type -> per-leg revs/branch (M-ED Phase 2)
```

### Task 2.1 — `tof_years` → `tof_seed_days` (g/h arcs)

**Files:** `descriptor.py`; test `tests/search/test_descriptor_tof.py`.

For g/h arcs the leg ToF seed is `tof_years * DAYS_PER_JULIAN_YEAR`
(`docs/spec.md:999`). For f/F arcs `tof_years` is `null`; the seed is derived
from the M:N resonance period (`M:N` ⇒ `M` spacecraft revs over `N` Earth years
≈ `N * DAYS_PER_JULIAN_YEAR` for the E-E resonant interval; document the
approximation in the docstring — it is a *seed*, refined by the corrector).

#### Failing test

```python
"""M-ED Phase 2: arc ToF seed in days (plan Phase 2; spec §16.7.7)."""
from __future__ import annotations

import pytest

from cyclerfinder.core.constants import DAYS_PER_JULIAN_YEAR
from cyclerfinder.search.descriptor import arc_tof_seed_days


def test_generic_tof_years_to_days() -> None:
    # g(1.4612,...) -> 1.4612 yr (spec §16.7.7; NOT the 526.02 deg psi param).
    assert arc_tof_seed_days("generic", tof_years=1.4612, resonance=None) == pytest.approx(
        1.4612 * DAYS_PER_JULIAN_YEAR
    )


def test_full_rev_tof_from_resonance() -> None:
    # F(3:2,...) -> resonant E-E interval ~ N=2 Earth years (seed only).
    assert arc_tof_seed_days("full-rev", tof_years=None, resonance="3:2") == pytest.approx(
        2.0 * DAYS_PER_JULIAN_YEAR
    )
```

Run → **red** → impl → **green**. Commit:

```
search/descriptor: tof_years/resonance -> ToF seed days (M-ED Phase 2)
```

### Task 2.2 — full descriptor list → `(per_leg_revs, per_leg_branch, tof_seed_days)`

**Files:** `descriptor.py`; test `tests/search/test_descriptor_parse.py`.

`parse_free_return_arcs(arcs: list[dict]) -> tuple[tuple[int,...],
tuple[str,...], tuple[float,...]]` maps the catalogue `free_return_arcs[]` list
(one arc per E-E leg) onto the three per-leg tuples. Use the S1L1 descriptor
`g(1.4612,…) G(2.8096,…)` → two generic arcs (revs `(0,0)`, branches
`("single","single")`, seeds `[1.4612 yr, 2.8096 yr]` in days), matching the
prototype's pinned arcs (`correct_s1l1_twoarc.py:40`).

#### Failing test

```python
"""M-ED Phase 2: free_return_arcs[] -> per-leg tuples (plan Phase 2)."""
from __future__ import annotations

import pytest

from cyclerfinder.core.constants import DAYS_PER_JULIAN_YEAR
from cyclerfinder.search.descriptor import parse_free_return_arcs


def test_s1l1_two_generic_arcs() -> None:
    arcs = [
        {"arc_type": "generic", "tof_years": 1.4612, "resonance": None,
         "raw_descriptor": "g(1.4612,526.02,Ll)"},
        {"arc_type": "generic", "tof_years": 2.8096, "resonance": None,
         "raw_descriptor": "G(2.8096,...)"},
    ]
    revs, branches, seeds = parse_free_return_arcs(arcs)
    assert revs == (0, 0)
    assert branches == ("single", "single")
    assert seeds[0] == pytest.approx(1.4612 * DAYS_PER_JULIAN_YEAR)
    assert seeds[1] == pytest.approx(2.8096 * DAYS_PER_JULIAN_YEAR)
```

Run → **red** → impl → **green**. Commit:

```
search/descriptor: parse_free_return_arcs -> per-leg topology + seeds (M-ED Phase 2)
```

### Task 2.3 — Phase-2 gate over the 12 sourced-descriptor rows + lint/type

**Files:** test `tests/search/test_descriptor_catalogue.py`.

Parametrise over the catalogue rows that carry a non-null `free_return_arcs[]`
(`docs/spec.md:1018`: 12 rows). Assert `parse_free_return_arcs` returns
well-formed tuples (lengths consistent, revs ≥ 0, branches in
`{"single","low","high"}`, seeds > 0). **Golden discipline:** assert only
*shape/structure* parsed from the descriptor strings, never a V∞ our code
computed. Skip the 3 gapped rows (`russell-ch4-8.165Gfh-f2`,
`russell-ch4-3.77Gh3`, `russell-ch4-5.66Gfh3`, `docs/spec.md:1020-1021`).

- [ ] Run green; `uv run ruff check . && uv run ruff format --check . && uv run mypy src tests`.

Commit:

```
test: descriptor parser gate over 12 sourced free_return_arcs rows (M-ED Phase 2)
```

---

## Phase 3 — `mode="ballistic"` wiring on `optimise_cell_ephemeris`

Extend, do not fork (spec §4). Add `mode: str = "maintenance"` to
`optimise_cell_ephemeris` (`optimize.py:1309-1321`). `"maintenance"` is the
existing path **byte-identical**; `"ballistic"` reuses `_target_period_sec`,
`_ephemeris_tof_seed_and_bounds`, `_resolve_t0_multi_seed`, calls
`ballistic_correct`, and maps to `OptimisationResult`.

### Task 3.0 — characterisation: default mode is byte-identical

**Files:** test `tests/search/test_optimise_ephemeris_mode.py`.

#### Failing test

```python
"""M-ED Phase 3: mode= kwarg on optimise_cell_ephemeris (plan Phase 3).

test_default_mode_unchanged is a CHARACTERISATION test: the default
mode="maintenance" path must be byte-identical to pre-M-ED. If it passes
immediately after the signature change (kwarg added, default branch untouched),
that is the desired outcome.
"""
from __future__ import annotations

import inspect

from cyclerfinder.search.optimize import optimise_cell_ephemeris


def test_mode_kwarg_exists_and_defaults_to_maintenance() -> None:
    sig = inspect.signature(optimise_cell_ephemeris)
    assert "mode" in sig.parameters
    assert sig.parameters["mode"].default == "maintenance"


def test_unknown_mode_raises() -> None:
    import pytest
    from cyclerfinder.core.ephemeris import Ephemeris
    from cyclerfinder.search.sequence import Cell

    cell = Cell(
        bodies=("E", "M"),
        sequence=("E", "M", "E"),
        period_k=2,
        per_leg_revs=(0, 0),
        per_leg_branch=("single", "single"),
    )
    with pytest.raises(ValueError, match="mode"):
        optimise_cell_ephemeris(cell, Ephemeris(model="circular"), vinf_cap=9.0, mode="bogus")
```

Run → **red** (no `mode` param).

#### Minimal impl — `optimize.py`

Add `mode: str = "maintenance"` to the signature. At the top of the body, after
the closed-loop check (`optimize.py:1398-1403`), validate
`mode in {"maintenance", "ballistic"}` else `ValueError`. Leave the entire
existing body as the `mode == "maintenance"` path (default → byte-identical).
Add the `mode == "ballistic"` branch in Task 3.1.

Run `test_mode_kwarg_exists…` + `test_unknown_mode_raises` → **green**. Run the
**whole existing optimise suite** to confirm the default path is unchanged:
`uv run pytest tests/search/test_optimize.py -q`. Commit:

```
search/optimize: add mode= kwarg to optimise_cell_ephemeris (default-inert) (M-ED Phase 3)
```

### Task 3.1 — ballistic branch: resolve epoch, run corrector, map result

**Files:** `optimize.py`; test `tests/search/test_optimise_ephemeris_mode.py` (slow).

In the `mode == "ballistic"` branch: reuse `target_period_sec =
_target_period_sec(cell)` (`optimize.py:1405`) and `seed_days, bounds =
_ephemeris_tof_seed_and_bounds(cell, target_period_sec)` (`optimize.py:1406`),
honour `tof_seed_days` override (`optimize.py:1414-1425`), resolve `t0_sec` via
`_resolve_t0_multi_seed` exactly as the maintenance path
(`optimize.py:1433-1444`). If `t0_sec is None`, return the same non-converged
sentinel (`optimize.py:1446-1464`). Else call `ballistic_correct(...)` with the
cell topology, the resolved `t0_sec`, `seed_days`, `target_period_sec`,
`vinf_cap`, `rp_factors`, and a `slack_leg` (default = the longest seed leg).
Map to `OptimisationResult`: build the closed `Cycler` via `_build_chain`-style
construction at the corrected `[t0, tofs]`, then `score(...)`; set `converged =
result.converged`, `constraints_satisfied = result.constraints_satisfied`,
`closure_residual_kms = result.max_residual_kms` (**a real residual**, unlike the
maintenance-ΔV proxy at `optimize.py:1504-1506`, spec §2.2).

#### Failing test (slow)

```python
@pytest.mark.slow
def test_ballistic_mode_closes_s1l1_returns_real_residual() -> None:
    """ballistic mode returns a REAL closure residual (V_inf-continuity), not
    the maintenance-dV proxy. NON-GOLDEN: closure is asserted, the V_inf value
    is OUR computation (S1L1 floors Mars ~6.4 — see project memory)."""
    from cyclerfinder.core.ephemeris import Ephemeris
    from cyclerfinder.search.optimize import optimise_cell_ephemeris
    from cyclerfinder.search.sequence import Cell

    cell = Cell(
        bodies=("E", "M"),
        sequence=("E", "M", "E", "E"),
        period_k=2,
        per_leg_revs=(0, 1, 2),
        per_leg_branch=("single", "single", "low"),
    )
    result = optimise_cell_ephemeris(
        cell,
        Ephemeris("astropy"),
        vinf_cap=9.0,
        priority_date_iso="2030-03-22",
        vinf_targets_kms={"E": 5.6, "M": 6.4},
        tof_seed_days=[154.0, 379.0, (1.4612 + 2.8096) * 365.25 - 154.0 - 379.0],
        mode="ballistic",
    )
    assert result.converged
    assert result.closure_residual_kms < 0.1  # real V_inf-continuity residual
```

Run → **red** → impl the branch → **green** (`-m slow`). Commit:

```
search/optimize: ballistic mode runs N-arc corrector, real closure residual (M-ED Phase 3)
```

### Task 3.2 — guard: maintenance mode result unchanged for the Aldrin cell

**Files:** test `tests/search/test_optimise_ephemeris_mode.py`.

Pin that `mode="maintenance"` on a known cell produces the same
`OptimisationResult` fields as before the `mode` kwarg (compare against the
pre-existing maintenance behaviour — reuse whatever cell the existing
`tests/search/test_optimize.py` maintenance test uses, or the Aldrin E-M-E
cell). This is the byte-identical-default contract (spec §4).

Run → green. Commit (may fold into 3.1 if trivial):

```
test: pin maintenance-mode result unchanged under mode= kwarg (M-ED Phase 3)
```

### Task 3.3 — Phase-3 lint/type gate

- [ ] `uv run pytest tests/search/test_optimise_ephemeris_mode.py tests/search/test_optimize.py -q`
  (+ `-m slow`), then ruff + mypy. Commit if cleanup needed.

---

## Phase 4 — Seeding ladder

Wire the ladder (spec §3): **descriptor → sourced anchor → coplanar → epoch
scan**, each rung used only when the rung above is absent for that row.

### Task 4.0 — ladder resolver skeleton (descriptor rung)

**Files:** create `src/cyclerfinder/search/seed_ladder.py`; test
`tests/search/test_seed_ladder.py`.

`resolve_seed(cell, *, free_return_arcs=None, anchor_vinf=None,
anchor_tofs=None, coplanar_tofs=None, ephem=None) -> SeedPlan` where `SeedPlan`
carries `per_leg_revs, per_leg_branch, tof_seed_days, source: Literal[
"descriptor","anchor","coplanar","scan"]`. Rung 1: if `free_return_arcs` is
present, use `parse_free_return_arcs` (Phase 2).

#### Failing test

```python
"""M-ED Phase 4: seeding ladder (plan Phase 4; spec §3)."""
from __future__ import annotations

from cyclerfinder.search.seed_ladder import resolve_seed
from cyclerfinder.search.sequence import Cell


def _s1l1_cell() -> Cell:
    return Cell(
        bodies=("E", "M"),
        sequence=("E", "M", "E", "E"),
        period_k=2,
        per_leg_revs=(0, 0, 0),
        per_leg_branch=("single", "single", "single"),
    )


def test_descriptor_rung_used_when_arcs_present() -> None:
    arcs = [
        {"arc_type": "generic", "tof_years": 1.4612, "resonance": None, "raw_descriptor": "g(...)"},
        {"arc_type": "generic", "tof_years": 2.8096, "resonance": None, "raw_descriptor": "G(...)"},
    ]
    plan = resolve_seed(_s1l1_cell(), free_return_arcs=arcs)
    assert plan.source == "descriptor"
    assert len(plan.tof_seed_days) == 2  # two E-E arcs
```

Run → **red** → impl rung 1 → **green**. Commit:

```
search/seed_ladder: descriptor rung (M-ED Phase 4)
```

### Task 4.1 — anchor rung (sourced V∞/ToF multisets)

**Files:** `seed_ladder.py`; test `tests/search/test_seed_ladder.py`.

Rung 2: when no descriptor but `anchor_tofs` (sourced transit ToFs) and
`anchor_vinf` (the `vinf_targets_kms` the epoch resolver needs,
`optimize.py:1355-1358`) are present, use them as the asymmetric seed. Sources
(spec §3.2): Jones VEM transit ToFs 309/259 & 268/223
(`docs/notes/2026-06-05-jones-aas17-577-vem-mining.md:393-416`); Sanchez Net
ME ~136 d / ~1026 d; Russell ch4 arc ToFs. **NOT** the S1L1 5.65/3.05 pair
(unverified provenance, spec §3.2 — diagnostic only).

#### Failing test

```python
def test_anchor_rung_used_when_no_descriptor() -> None:
    plan = resolve_seed(
        _emevve_cell(),
        free_return_arcs=None,
        anchor_tofs=(309.0, 259.0),       # Jones EMEVVE transit legs (sourced)
        anchor_vinf={"E": 4.72, "M": 2.50},  # Jones Table 2 (sourced)
    )
    assert plan.source == "anchor"
    assert plan.tof_seed_days[0] == 309.0
```

(`_emevve_cell` builds the `E-M-E-V-V-E` cell with `period_basis=("E","M")`.)
Run → **red** → impl rung 2 → **green**. Commit:

```
search/seed_ladder: sourced-anchor rung (Jones/Sanchez/Russell ToFs) (M-ED Phase 4)
```

### Task 4.2 — coplanar + epoch-scan rungs

**Files:** `seed_ladder.py`; test `tests/search/test_seed_ladder.py`.

Rung 3 (coplanar warm start): when neither descriptor nor anchor, use
`optimise_cell_idealized` (`optimize.py:955`) leg ToFs as the asymmetric seed
(spec §3.3). Rung 4 (epoch scan, last resort): equispaced seed + the
`_resolve_t0_multi_seed` ±10 yr window default (`optimize.py:1293`); `source =
"scan"`. Test that the ladder degrades in order (descriptor > anchor > coplanar >
scan) and that `source` reflects the rung actually used.

Run → **red** → impl → **green**. Commit:

```
search/seed_ladder: coplanar + epoch-scan fallback rungs (M-ED Phase 4)
```

### Task 4.3 — Phase-4 lint/type gate

- [ ] `uv run pytest tests/search/test_seed_ladder.py -q`, ruff + mypy. Commit
  if cleanup.

---

## Phase 5 — discover() V3 gate, Jones VEM headline gate, S1L1 diagnostic, census

### Task 5.0 — replace the stale `discover.py:99-106` V3 branch (resolves tracked task #109)

**Files:** `src/cyclerfinder/data/discover.py`; test
`tests/data/test_discover_v3_gate.py`.

> **NOTE (task #109):** the design names this as "resolves tracked task #109".
> Task #109 could **not** be located in the live repo (`data/OUTSTANDING.md`,
> `docs/`) — see Self-Review "Unverifiable claims". Record the task id in the
> commit message as the design instructs, but do not fabricate a tracker link.

The current branch (`discover.py:99-106`) calls `optimise_cell_ephemeris(...)`
expecting `NotImplementedError` (`# raises until M6b lands`, `discover.py:105`)
and unconditionally sets `level="V3"`. It no longer raises. Replace it with a
**real ballistic-closure gate**: call `optimise_cell_ephemeris(result.cell,
ephem, vinf_cap=vinf_cap, mode="ballistic", priority_date_iso=…,
vinf_targets_kms=…)` and set `level="V3"` **only if**
`v3_result.constraints_satisfied` (real V∞-continuity closure + bend-feasible +
V∞-cap), else leave `level="V2"`.

#### Failing test

```python
"""M-ED Phase 5: discover V3 branch is a real ballistic-closure gate, not a
stub-raise assumption (plan Phase 5 Task 5.0; spec §0 finding 3)."""
from __future__ import annotations

import inspect

from cyclerfinder.data import discover as discover_mod


def test_v3_branch_no_longer_assumes_notimplemented() -> None:
    src = inspect.getsource(discover_mod._auto_validate)
    # The stale comment / assumption must be gone.
    assert "raises until M6b lands" not in src
    assert "mode=\"ballistic\"" in src or "mode='ballistic'" in src
```

Run → **red** → impl the real gate → **green**. Then the existing
`tests/data/test_discover.py` suite must stay green (the `enable_v3=False`
default path is unchanged). Commit:

```
data/discover: real ballistic-closure V3 gate replacing stale stub branch (resolves #109)

The V3 branch assumed optimise_cell_ephemeris still raised NotImplementedError
(discover.py:99-106, "raises until M6b lands"); it no longer does. V3 now runs
the M-ED ballistic mode and gates level="V3" on real V_inf-continuity closure
(constraints_satisfied), not on an exception. Spec §0 finding 3.
```

### Task 5.1 — Jones VEM headline gate: flip the xfail (HEADLINE GATE)

**Files:** `tests/test_vem_rediscovery.py` (flip the existing xfail,
`test_vem_rediscovery.py:220-258`); helper fixtures.

This is the **headline rediscovery gate** (Approval Q2). Re-target the existing
aspirational xfail from the circular-coplanar idealized optimiser to the **M-ED
ballistic mode on the real ephemeris**, asserting convergence to the **sourced
Jones VEM V∞ multiset** within tolerance.

**Headline gate definition (verbatim — this is the binding criterion):**

> **For each Jones VEM member row (`jones-2017-vem-emevve-outbound`,
> `jones-2017-vem-meevem-inbound`), the M-ED ballistic corrector
> (`optimise_cell_ephemeris(..., mode="ballistic")`), seeded via the sourced-
> anchor rung (the row's sourced transit ToFs and per-encounter V∞ targets from
> AAS 17-577 Tables 2/3), converges to a closed ballistic chain on DE440 whose
> per-encounter V∞ magnitudes match the row's sourced `vinf_kms_at_encounters`
> multiset (compared as a sorted multiset, EXPECTED = the catalogue's sourced
> values) within `VEM_VINF_TOL_KMS = 0.5 km/s` per encounter, AND whose closure
> residual is `< 0.1 km/s`, AND all flybys are bend-feasible. The EXPECTED side
> is the published Jones multiset only; no value our own code computed is ever
> the EXPECTED side (golden discipline, project memory
> `feedback_golden_tests_sourced_only.md`).**

**Flip criteria (explicit, Approval-mandated):** the xfail flips to a green
(non-xfail) assertion when, and only when, the criterion above holds for **at
least one** of the two member rows at `VEM_VINF_TOL_KMS = 0.5`. The tolerance is
tied to the sourced values (Jones rounds V∞ to 0.01 km/s; 0.5 km/s absorbs the
zero-SOI patched-conic vs. our DE440 Lambert model difference, spec §5). If only
one row flips, the other stays xfail with a documented per-row finding.

#### Test (xfail-first, flipped on convergence)

```python
"""M-ED HEADLINE GATE: Jones VEM ballistic rediscovery (plan Phase 5 Task 5.1).

GOLDEN DISCIPLINE: EXPECTED = the catalogue's SOURCED vinf_kms_at_encounters
(AAS 17-577 Tables 2/3). The corrector output is the side under test. No
self-computed value is ever the EXPECTED side.

RISK (spec §7): the S1L1 corrector family floors Mars V_inf ~6.4 km/s; the Jones
Mars targets are 2.50/2.79 (EMEVVE) and ~3.12/3.85 (MEEVEM). Convergence is NOT
assumed. Until it converges within VEM_VINF_TOL_KMS this stays xfail with the
finding recorded; the STOP/report branch (Task 5.4) governs that outcome.
"""
from __future__ import annotations

import pytest
import yaml

from tests._catalogue_loader import CATALOGUE_PATH

VEM_VINF_TOL_KMS = 0.5  # tied to sourced Jones rounding + model difference (spec §5)
_VEM_MEMBERS = ("jones-2017-vem-emevve-outbound", "jones-2017-vem-meevem-inbound")


def _row(entry_id: str) -> dict:
    for row in yaml.safe_load(CATALOGUE_PATH.read_text()):
        if row["id"] == entry_id:
            return row
    raise AssertionError(f"catalogue row {entry_id!r} not found")


def _sourced_vinf_multiset(entry_id: str) -> list[float]:
    return sorted(float(e["vinf_kms"]) for e in _row(entry_id)["vinf_kms_at_encounters"])


@pytest.mark.slow
@pytest.mark.xfail(
    reason="M-ED HEADLINE GATE: ballistic VEM rediscovery to the sourced Jones "
    "multiset within 0.5 km/s. Open until the corrector converges (spec §7 risk: "
    "S1L1 family floors Mars V_inf ~6.4 vs Jones 2.5-3.9). Flip when at least one "
    "member row converges; see plan Phase 5 Task 5.1 flip criteria.",
    strict=False,
)
@pytest.mark.parametrize("entry_id", _VEM_MEMBERS)
def test_jones_vem_ballistic_rediscovers_sourced_multiset(entry_id: str) -> None:
    from cyclerfinder.core.ephemeris import Ephemeris
    from cyclerfinder.search.optimize import optimise_cell_ephemeris
    # Build the cell from the row (sequence_canonical, period.pair beat-token ->
    # period_basis via the loader's _anchor_pair_from_period_pair), seed via the
    # anchor rung from the row's sourced transit ToFs + V_inf targets, run
    # mode="ballistic". (Helper construction omitted here; see Task 5.1 fixtures.)
    expected = _sourced_vinf_multiset(entry_id)
    result = ...  # optimise_cell_ephemeris(cell, Ephemeris("astropy"), mode="ballistic", ...)
    assert result.converged and result.constraints_satisfied
    got = sorted(
        max(float(np.linalg.norm(e.vinf_in)), float(np.linalg.norm(e.vinf_out)))
        for e in result.best_cycler.encounters
    )
    for g, x in zip(got, expected, strict=True):
        assert g == pytest.approx(x, abs=VEM_VINF_TOL_KMS)
```

Run `-m slow`: registers `xfail` (or `xpass` → a reviewable finding that the gate
converged → then Task 5.4 confirms and removes the xfail). Commit:

```
test: M-ED headline gate — Jones VEM ballistic rediscovery to sourced multiset (xfail)

Re-targets the M8-Core aspirational xfail from the circular-coplanar optimiser
to the M-ED ballistic mode on DE440, asserting per-encounter V_inf against the
SOURCED AAS 17-577 Tables 2/3 multiset within 0.5 km/s. xfail until convergence;
flip criteria + tolerance documented (Approval Q2). Golden discipline preserved.
```

### Task 5.2 — Sanchez Net regime cross-validation gate (promote prototype)

**Files:** test `tests/search/test_correct_sanchez_regime.py` (slow).

Promote the prototype's already-demonstrated DE440 2030-2034 near-ballistic
E-M-E-E closure (spec §5) to a test asserting `ballistic_correct` closes ≥1
chain in that regime with V∞ ≤ cap and bend-feasible. **NON-GOLDEN** for the V∞
value (our computation); the gate is *closure exists in the near-ballistic
regime*, not a sourced V∞ assertion (the Sanchez event V∞ are sourced but the
S1L1 family floors elsewhere — assert regime, cap, feasibility only).

Run `-m slow` → green (the prototype proves this closes). Commit:

```
test: Sanchez-regime near-ballistic closure gate for the corrector (M-ED Phase 5)
```

### Task 5.3 — S1L1 stays non-gating diagnostic (no green 5.65/3.05 assertion)

**Files:** confirm-only; possibly a docstring note on the existing xfails.

Per spec §5 + Approval Q2, S1L1 is demoted to a non-gating diagnostic:
- Keep `EXPECTED_SKIPS["s1l1-2syn-em-cpom"]` (`real_closure.py:222-236`) and
  `EXPECTED_SKIPS["mcconaghy-2006-em-k2"]` (`real_closure.py:237-243`) **as-is**.
- Keep the M5 xfail `test_2syn_em_rediscovers_5_65_kms_earth`
  (`test_catalogue_rediscovery.py:12-13`) and the discover end-to-end variant
  (`tests/data/test_discover.py:110-126`) **xfail** — no green test may assert
  5.65/3.05 (unverified provenance).
- [ ] Add (only if not already present) a one-line cross-reference in those
  xfail reasons pointing at the M-ED headline gate as the replacement criterion,
  per spec §5 ("re-target it to the Jones VEM gate"). Do **not** loosen or flip
  them.

Commit (only if a docstring note is added):

```
test: cross-reference S1L1 xfails to the M-ED Jones VEM headline gate (M-ED Phase 5)
```

### Task 5.4 — STOP/report branch (honest risk register) + census ratchet

**Files:** test runs; possibly `tests/test_catalogue_rediscovery.py`
`EXPECTED_COVERAGE`.

**STOP/report branch (binding, spec §7 risk register).** After Tasks 5.1-5.2,
evaluate the headline gate outcome:

- **If the Jones VEM gate converges** (Task 5.1 xpasses for ≥1 row): remove the
  `xfail` marker on the converging row(s), promote to a strict green gate at
  `VEM_VINF_TOL_KMS`, and — if the row's loader class changes from
  `CONSTRUCTIBLE_MULTIBODY` to constructible — update the
  `EXPECTED_COVERAGE` census ratchet (`test_catalogue_rediscovery.py`) **in the
  same commit** (the ratchet rule), keeping the coverage-audit invariant intact
  (spec §5; "no entry vanishes silently").
- **If the gate does NOT converge** (the expected risk: the corrector cannot
  reach the Jones Mars-V∞ 2.5-3.9 from the family that floors ~6.4): **STOP — do
  not loosen the tolerance, do not assert a self-computed value, do not flip the
  xfail.** Leave the xfail in place with a per-row finding documenting the best
  achieved residual and the V∞ gap, and record the result in `data/OUTSTANDING.md`
  as open research (the corrector + descriptor + mode + ladder of Phases 1-4 are
  the shipped deliverable; the VEM convergence is the open item). The plan ships
  either way; the headline gate's *existence* (xfail-first, sourced, with flip
  criteria) is the Phase-5 deliverable, not its green status.

- [ ] Whichever branch: census ratchet (`EXPECTED_COVERAGE`) re-derived from the
  live loader (do not copy stale numbers — see M8-Core Revision R1 delta 4), and
  the coverage-audit census invariant holds.

Commit:

```
test: M-ED headline-gate outcome + census ratchet (STOP/report branch) (M-ED Phase 5)
```

### Task 5.5 — full-suite + lint/type gate; OUTSTANDING note

- [ ] `uv run pytest -m "not slow"` green; `uv run pytest -m slow` evaluated
  (xfail/xpass per Task 5.4).
- [ ] `uv run ruff check . && uv run ruff format --check . && uv run mypy src tests` clean.
- [ ] Update `data/OUTSTANDING.md`: the corrector/descriptor/mode/ladder shipped;
  the VEM headline-gate status (converged-and-green, or xfail-with-finding per
  the STOP branch); the S1L1 ~6.4 floor still open vs the (future) McConaghy
  4.7/5.0 anchor. *(Defer this edit if shared-doc concurrency forbids it at run
  time; flag in the report.)*

Commit:

```
docs: M-ED OUTSTANDING refresh — corrector shipped, VEM gate status, S1L1 floor (M-ED Phase 5)
```

---

## Out of scope (explicit carve-outs, spec §6.2)

- **Powered cyclers / Aldrin BVP** (`bvp.py:92` `solve_powered_periodic_cycler`
  pins E-M-E via `optimise_aldrin_maintenance_dv`, `bvp.py:61,149`). M-ED reports
  bend-infeasible chains honestly but does not solve the powered TCM problem.
  Generalising `solve_powered_periodic_cycler` off the Aldrin lock is deferred to
  M7 (spec §4 "recommend deferring").
- **M7 per-family TCM budget** — the maintenance-ΔV mode already computes
  per-cycle ΔV (`_maintenance_dv_chain`); a per-family ceiling is a reporting
  layer, deferred (Approval Q5).
- **Option C basin-locked maintenance-ΔV refinement** — build B first; layer C
  only if a consumer (Forge/M7) needs the TCM cost (spec §7).
- **Rotating-frame multi-lap drift** (M-3D's `REAL_DRIFT_TOLERANCE_KM`) — a
  separate validation axis (the existing V2 gate, `discover.py:91-97`),
  unchanged.
- **`find_cyclers`** — circular-only by contract (`optimize.py:1592-1598`);
  untouched.
- **Blind (anchorless) discovery** — epoch resolution requires V∞ targets
  (`optimize.py:1356-1358`); M-ED rediscovers *sourced* families. Blind search is
  Forge Phase 4.

---

## Definition of done

0. `search/correct.py` — `ballistic_correct` closes the S1L1 two-arc chain on
   DE440 (slow, non-golden regression). Pure (core-only imports). **Green.**
1. `search/descriptor.py` — `parse_free_return_arcs` maps the 12 sourced-
   descriptor rows to well-formed `(revs, branch, tof_seed_days)`. **Green.**
2. `optimise_cell_ephemeris(mode="ballistic")` runs the corrector and reports a
   **real** closure residual; `mode="maintenance"` default byte-identical. **Green.**
3. `search/seed_ladder.py` — descriptor → anchor → coplanar → scan, `source`
   reflects the rung used. **Green.**
4. `data/discover.py` V3 branch is a real ballistic-closure gate (no
   "raises until M6b lands"); resolves the design's task #109 reference. **Green.**
5. The Jones VEM **headline gate** exists (xfail-first, sourced multiset EXPECTED,
   documented flip criteria + `VEM_VINF_TOL_KMS=0.5`), flipped green iff it
   converges (Task 5.4 STOP/report branch governs the no-converge outcome).
6. S1L1 5.65/3.05 stays **non-gating** — no green test asserts it; the M5 xfail
   and `EXPECTED_SKIPS` entries are unchanged.
7. No assertion anywhere uses a V∞/`(a,e)` value our own code computed as the
   EXPECTED side — only sourced Jones/Sanchez/Russell multisets, sourced ToFs,
   and feasibility/closure predicates.
8. `uv run pytest -m "not slow"`, ruff, ruff format, mypy all clean. Census ratchet
   intact.

---

## Self-Review

### Spec coverage (every spec contract → plan task)

| Spec item | Plan location |
|---|---|
| §2.1 N-arc model, free vars `[t0, ToFs]`, period elimination, residuals, bend post-hoc | Phase 1 (Tasks 1.0-1.4) |
| §2.2 corrector new mode, both modes kept, real residual | Phase 3 (Tasks 3.0-3.1) |
| §3.1 descriptor genome → revs/branch/tof_seed | Phase 2 (Tasks 2.0-2.2) |
| §3.2 sourced anchors (Jones/Sanchez/Russell), NOT S1L1 5.65/3.05 | Phase 4 Task 4.1; Phase 5 Task 5.3 |
| §3.3 coplanar warm start; §3.4 epoch scan | Phase 4 Task 4.2 |
| §4 `search/correct.py`, `search/descriptor.py`, `mode=` extend-not-fork, discover wiring | Phases 1-3, 5 Task 5.0 |
| §4 defer `solve_powered_periodic_cycler`; `find_cyclers` untouched | Out of scope |
| §5 Jones VEM headline gate (xfail flip) | Phase 5 Task 5.1 (HEADLINE) |
| §5 Sanchez regime gate | Phase 5 Task 5.2 |
| §5 S1L1 non-gating; keep `EXPECTED_SKIPS` + M5 xfail | Phase 5 Task 5.3 |
| §5 census ratchet intact | Phase 5 Task 5.4 |
| §0 finding 3 / discover V3 real gate (#109) | Phase 5 Task 5.0 |
| §7 honest risk register, STOP/report branch | Goal honesty boundary; Phase 5 Task 5.4 |
| §6.1 M-3D not a dependency; M-L landed | Architecture (frame-freedom) |
| §6.2 non-goals | Out of scope |
| Approval Q1-Q6 (ballistic objective, Jones gate, no M-3D, M-L landed, Option B MVP, descriptor-first) | Goal + phasing + Out of scope |

### Placeholder scan
No `TODO`/`FIXME`/`...`-as-impl placeholders in production-code instructions.
The single `...` appears inside the Task 5.1 *test sketch* (the cell-construction
helper), explicitly marked "Helper construction omitted here; see Task 5.1
fixtures" — the construction recipe (build cell from `sequence_canonical` +
`_anchor_pair_from_period_pair`, seed via anchor rung, `mode="ballistic"`) is
spelled out in prose in the same task. All file paths are concrete and absolute-
within-repo. All commit messages are exact.

### Type consistency
`BallisticClosureResult` is a frozen dataclass with a `constraints_satisfied`
property (matches `OptimisationResult`'s field semantics). `ballistic_correct`
signature matches spec §4. `mode: str = "maintenance"` is additive-keyword on
`optimise_cell_ephemeris` (existing signature `optimize.py:1309-1321`
unchanged otherwise → default byte-identical). `SeedPlan.source` is a `Literal`.
`parse_free_return_arcs` returns `tuple[tuple[int,...], tuple[str,...],
tuple[float,...]]`. `Cell` fields used (`bodies, sequence, period_k,
per_leg_revs, per_leg_branch, period_basis`) match `sequence.py:106-111`.

### Unverifiable claims (flagged, not papered over)
1. **Task #109.** The design calls the discover V3 fix "resolves tracked task
   #109", but **task #109 is not findable** anywhere in the live repo
   (`data/OUTSTANDING.md` has no #109; no `docs/` match for a discover/V3/stub
   task #109). The "roadmap:109" tokens in the spec are *line-number* citations,
   not a task id. The plan records #109 in the commit message as instructed but
   does not invent a tracker entry. **If the user has an external tracker, the
   #109 link must be confirmed there before the commit lands.**
2. **MEEVEM cycle-2 Mars V∞ "~3.12".** The mining note (Task 5.1 reason) gives
   MEEVEM transit V∞ as "Mars 3.85 / Earth 3.48 … Earth 2.98 / Mars 3.12 region
   cycle 2 (see Table 3 ordering)"
   (`docs/notes/2026-06-05-jones-aas17-577-vem-mining.md:412-414`) — the 3.12 is
   the note author's read of Table 3 ordering, not a clean cell. The headline
   gate uses the **full sourced `vinf_kms_at_encounters` multiset from the
   catalogue row** (a sorted-multiset compare), which sidesteps any per-leg
   ordering ambiguity, so this does not affect the gate; flagged for honesty.
3. **MEEVEM catalogue row V∞ values** were not exhaustively read line-by-line
   (the EMEVVE row was, `catalogue.yaml:2070-2100+`); the MEEVEM row
   (`catalogue.yaml:2288`) is asserted to carry a full sourced multiset on the
   basis of the mining note and the row header. The gate reads the multiset from
   the live YAML at test time, so a partial multiset would surface as a length
   mismatch in the `zip(..., strict=True)` — a visible failure, not a silent pass.
