# Ellison FBS Path-B Re-transcription Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Execute ONE phase at a time; the human reviews between phases. This is a MULTI-DAY effort (#226).

**Goal:** Bring the Ellison 2018 forward-backward-shooting (FBS) MGAnDSMs match-point formulation into the DSM lane so a DSM leg can be corrected with an *analytic* Jacobian — no Lambert solver, no finite differencing.

**Architecture:** A new pure `core/` module assembles the mid-phase match-point defect of a single two-body leg carrying one interior impulse (Δv as an explicit decision variable, the MGAnDSMs n=1 case) by *forward* and *backward* Keplerian propagation to a common match point, together with its analytic Jacobian built from the analytic two-body STM (`core/kepler_stm.shepperd_stm`). The massless cycler-leg simplification collapses Ellison's maneuver transition matrix (MTM, Eq. 29) to identity plus the Eq. 42 velocity-slot connection, so for a massless leg the only non-trivial matrices in the Eq. 31-32 chain are the STMs. Later phases extend to phase-TOF partials (Pitkin Eqs. 43-44), the multi-arc chain assembler, and wiring an analytic-Jacobian option into the DSM corrector.

**Tech Stack:** Python 3.11, numpy, scipy (existing); uv-managed venv; pytest; ruff + mypy (pre-commit hook). Validation is the FD-vs-analytic CONSISTENCY pattern (identical discipline to `core/kepler_stm` and `nbody/flyby_gradients`).

**Source (cite by publication only in committed content):** D. H. Ellison, B. A. Conway, J. A. Englander, M. T. Ozimek, "Analytic Gradient Computation for Bounded-Impulse Trajectory Models Using Two-Sided Shooting," *Journal of Guidance, Control, and Dynamics*, Vol. 41, No. 7, 2018, pp. 1449-1462, doi:10.2514/1.G003077. Mining note: `docs/notes/2026-06-10-ellison-2018-analytic-gradients-mining.md`. The A1-A6 NLP-scaling caveat (see `nbody/flyby_gradients.py` docstring) carries over: where a printed Ellison expression carries an inconsistent nondimensionalisation factor, implement the dimensionally-consistent form and validate by FD, noting it.

---

## Discipline (applies to every phase)

- **TDD, bite-sized commits.** Failing test first; minimal implementation; pass; commit. One logical change per commit, single-line commit message, pathspec-only `git add`.
- **Consistency validation, never goldens.** Ellison prints no unit-level numeric gradient (mining note §6). Every analytic-Jacobian test is FD-vs-analytic agreement; label the test docstring as a consistency test. Do NOT invent sourced numeric gradient values.
- **Additive only.** New `core/` module + new functions. Do NOT change the behaviour of `core/sims_flanagan.py`, `search/dsm_leg.py`, or any other lane — other lanes depend on them. Wire analytic Jacobians as an OPTION.
- **Pre-commit:** `uv run ruff check .`, `uv run ruff format --check .`, `uv run mypy` clean before every commit (the hook enforces it). Re-check `git status` before committing. Never `git reset`/`git stash`. No AI-attribution in commits. Never branch — work on `main`.
- **Units/conventions:** km, km/s, s, km^3/s^2 throughout, matching `core/kepler_stm` (state ordered `[r; v]`, position rows/cols first; `dr/dv0` block carries seconds, `dv/dr0` block 1/seconds). Heliocentric `MU_SUN_KM3_S2` default.

---

## Why a new module (not editing `sims_flanagan.py`)

`core/sims_flanagan.py` is the MGALT (Sims-Flanagan low-thrust) transcription: `N` equal-time segments, one bounded impulse per segment midpoint, a 7-element defect *including a mass row*, and a per-segment Δv schedule that is NOT a free decision vector in the cycler sense. The DSM lane's leg is the MGAnDSMs case: a small number (here one) of interior impulses whose Δv is *directly* a decision variable, and our cycler legs are **massless** (no mass row — mining note §7, Path B). These are different transcriptions sharing only the FBS skeleton. Editing `sims_flanagan.py` would either (a) bolt a foreign massless/DSM mode onto an MGALT class other lanes' tests pin, risking silent behaviour change, or (b) duplicate its skeleton anyway. A clean sibling module `core/fbs_match_point.py` is cheaper, isolates the massless simplification, and keeps `sims_flanagan` untouched. (Justified per the task's "new sibling module if cleaner — your call" clause.)

The existing `sims_flanagan` skeleton needs **no refactor first**: its forward/backward propagate-to-match-point pattern is the conceptual template, but its mass-bearing 7-vector and per-segment-schedule API do not fit the massless single-DSM leg. Phase 1 is therefore the new module, not a refactor.

---

## Phase Breakdown (titles + one-line each)

- **Phase 1 — Single massless two-body DSM leg: analytic match-point defect + its Jacobian (FD-validated).** The foundational, lowest-risk piece: one leg, one interior impulse, forward+backward Kepler propagation to the match point, 6-element `[Δr; Δv]` defect, and its analytic Jacobian w.r.t. the leg's decision variables via `shepperd_stm` and the Eq. 42 connection — validated FD-vs-analytic on several random states.
- **Phase 2 — Phase-TOF partials (Pitkin Eqs. 43-44) for the single leg.** Add the `∂(defect)/∂Δt_p` columns using the Lagrange-coefficient time derivatives, with `∂Δt_k/∂Δt_p = α_k` (Eq. 58), extending Phase 1's Jacobian to the time decision variables; FD-validated.
- **Phase 3 — Boundary v∞ / epoch columns + ephemeris-moving endpoints.** Add the `∂(defect)/∂v∞_out`, `∂(defect)/∂v∞_in` columns and the moving-boundary terms (Eqs. 57, 59-61) so a leg pinned to ephemeris bodies has a complete analytic Jacobian; FD-validated.
- **Phase 4 — Multi-arc chain assembler (Eqs. 31-32) over M legs.** Compose per-leg STM/connection products into the stacked chain defect + block-sparse Jacobian for a sequence of legs joined at flybys; FD-validated against the existing chain machinery's FD.
- **Phase 5 — Flyby-continuity coupling between consecutive legs.** Join the chain with the already-implemented `nbody/flyby_gradients` (Eqs. 3-4, A1-A6) so inter-leg v∞ continuity/turn feasibility enter the assembled Jacobian analytically.
- **Phase 6 — Wire the analytic Jacobian as an OPTION into the DSM corrector.** Add an opt-in path to `search/dsm_leg.py`'s `dsm_chain_correct` that supplies the analytic Jacobian to `least_squares` (`jac=`), default-off (FD path unchanged); validate convergence parity + speedup on a known case.

Each phase leaves the DSM lane working (additive, default-off) and is independently testable.

---

## PHASE 1 (execute now)

### Task 1: New module skeleton + single-leg match-point defect

**Files:**
- Create: `src/cyclerfinder/core/fbs_match_point.py`
- Test: `tests/core/test_fbs_match_point.py`

**Concept.** A *massless* DSM leg, MGAnDSMs n=1 (mining note §3-§7):
- Leg start state `(r0, v0)` at the left boundary; the leg flies for total time `tof_s`.
- One interior impulse Δv (a 3-vector decision variable) applied at burn-fraction `alpha in (0,1)` of the TOF, i.e. at `t_burn = alpha * tof_s` after start. Forward dynamics: coast `t_burn` from `(r0, v0)`, apply `v <- v + Δv`, this is state `X_k^+`.
- Right boundary state `(rf, vf)` at the right boundary; backward dynamics coast `-(tof_s - t_burn)` from `(rf, vf)` to the match point (which we place *at the impulse*, matching Ellison's mid-phase match point sitting at a maneuver point for n=1).
- **Match point = the impulse point.** Forward match state: `X^F = (r_fwd, v_fwd + Δv)` where `(r_fwd, v_fwd) = propagate(r0, v0, t_burn)`. Backward match state: `X^B = propagate(rf, vf, -(tof_s - t_burn))`. (Massless: no mass row; MTM Eq. 29 = identity, the only thing the impulse does to the forward state is `v += Δv` — the Eq. 42 connection.)
- **Defect (Eq. 2, mass row dropped):** `c_mp = X^B - X^F = [r^B - r^F; v^B - v^F]`, a 6-vector. Zero ⇒ the leg is dynamically consistent: the coasted-and-burned forward arc meets the back-propagated right boundary.

This is Lambert-free: the post-impulse velocity is `v_fwd + Δv` with Δv a free variable, not solved from a Lambert mismatch.

**Step 1: Write the failing test (defect is zero on a self-consistent leg)**

Construct a self-consistent leg by *forward* construction so the defect must vanish: pick `(r0, v0)`, `tof_s`, `alpha`; coast to the burn point; choose any Δv; the post-impulse state continued to `tof_s` defines `(rf, vf)`. Feeding those back must give a ~0 defect.

```python
"""Tests for cyclerfinder.core.fbs_match_point (Ellison 2018 FBS, massless DSM leg).

CONSISTENCY-TEST PATTERN (same discipline as tests/core/test_kepler_stm.py and
tests/nbody for flyby_gradients): Ellison 2018 prints no unit-level numeric
gradient (mining note 2026-06-10 §6), so the analytic match-point-defect Jacobian
is validated against central differences of the defect itself, never a sourced
golden. The zero-defect test below uses a forward-CONSTRUCTED self-consistent leg
(no external numbers), not a published value.
"""
from __future__ import annotations

from math import pi, sqrt

import numpy as np
import pytest

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, SECONDS_PER_DAY
from cyclerfinder.core.kepler import coe_to_rv, propagate
from cyclerfinder.core.fbs_match_point import (
    FbsLeg,
    match_point_defect,
    match_point_defect_jacobian,
)

_V_CIRC_1AU = sqrt(MU_SUN_KM3_S2 / AU_KM)


def _self_consistent_leg(r0, v0, tof_s, alpha, dv):
    """Build (rf, vf) so the leg with this (Δv, alpha) has ~zero defect."""
    t_burn = alpha * tof_s
    r_b, v_b = propagate(r0, v0, t_burn)
    v_b_post = v_b + dv
    rf, vf = propagate(r_b, v_b_post, tof_s - t_burn)
    return rf, vf


def test_defect_zero_on_self_consistent_leg() -> None:
    r0 = np.array([AU_KM, 0.0, 0.0])
    v0 = np.array([0.0, _V_CIRC_1AU, 0.03 * _V_CIRC_1AU])
    tof_s = 250.0 * SECONDS_PER_DAY
    alpha = 0.4
    dv = np.array([0.2, -0.1, 0.05])
    rf, vf = _self_consistent_leg(r0, v0, tof_s, alpha, dv)
    leg = FbsLeg(r0=r0, v0=v0, rf=rf, vf=vf, tof_s=tof_s, alpha=alpha)
    c = match_point_defect(leg, dv)
    assert c.shape == (6,)
    assert float(np.linalg.norm(c[0:3])) / AU_KM < 1.0e-10
    assert float(np.linalg.norm(c[3:6])) / _V_CIRC_1AU < 1.0e-10
```

**Step 2: Run it to verify it fails**

Run: `uv run pytest tests/core/test_fbs_match_point.py::test_defect_zero_on_self_consistent_leg -v`
Expected: FAIL — `ImportError` / module not found.

**Step 3: Write minimal implementation (module + dataclass + defect)**

```python
r"""Ellison 2018 forward-backward-shooting (FBS) match point — massless DSM leg.

Path-B re-transcription of the DSM lane (#226; mining note
``docs/notes/2026-06-10-ellison-2018-analytic-gradients-mining.md`` §3, §7). A
single two-body leg carries ONE interior impulse Δv (the MGAnDSMs n=1 case,
Ellison Sec. II.D). The leg is propagated FORWARD from its left boundary and
BACKWARD from its right boundary to a common match point placed AT the impulse;
the position/velocity mismatch there (the defect, Ellison Eq. 2 with the mass row
dropped — our cycler legs are massless) must vanish for the leg to be dynamically
consistent. Unlike the Takao/Lambert DSM leg (``search/dsm_leg.py``) the
post-impulse velocity is ``v_fwd + Δv`` with Δv an explicit decision variable —
there is NO Lambert solve anywhere (the structural payoff of two-sided shooting,
mining note §1).

Massless simplification (mining note §7, Path B): Ellison's maneuver transition
matrix (Eq. 29) collapses to the identity, and the impulse enters the forward
match state only through the Eq. 42 velocity-slot connection ``∂X_k^+/∂Δv =
[0; I]``. The only non-trivial matrices in the Eq. 31-32 match-point chain are
therefore the analytic two-body STMs (``core/kepler_stm.shepperd_stm``).

Phase 1 scope: the single-leg defect and its analytic Jacobian w.r.t. Δv and the
boundary states. Phase-TOF partials (Pitkin Eqs. 43-44), the multi-leg chain
(Eqs. 31-32), and the corrector wiring are later phases (see the plan).

Source: D. H. Ellison, B. A. Conway, J. A. Englander, M. T. Ozimek, "Analytic
Gradient Computation for Bounded-Impulse Trajectory Models Using Two-Sided
Shooting," *Journal of Guidance, Control, and Dynamics*, 41(7), 2018,
pp. 1449-1462, doi:10.2514/1.G003077 (Eqs. 2, 17, 29, 31-32, 42).

CONSISTENCY DISCIPLINE: Ellison prints no unit-level numeric gradient (mining
note §6); the Jacobian is validated against central differences of the defect,
never a sourced golden. The A1-A6 NLP-scaling caveat noted in
``nbody/flyby_gradients.py`` carries over.

Units: km, km/s, s, km^3/s^2. State ordered ``[r; v]`` (position first), matching
``core/kepler_stm``.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import MU_SUN_KM3_S2
from cyclerfinder.core.kepler_stm import shepperd_stm

Vec3 = NDArray[np.float64]   # (3,)
Vec6 = NDArray[np.float64]   # (6,)
Mat6 = NDArray[np.float64]   # (6, 6)


class FbsMatchPointError(Exception):
    """Errors raised by the FBS match-point model."""


@dataclass(frozen=True)
class FbsLeg:
    r"""A single massless two-body leg with one interior impulse (MGAnDSMs n=1).

    Attributes
    ----------
    r0, v0:
        Left-boundary inertial state, ``(3,)`` km / km/s.
    rf, vf:
        Right-boundary inertial state, ``(3,)`` km / km/s.
    tof_s:
        Total leg time of flight, s (> 0).
    alpha:
        Burn-index fraction in ``(0, 1)`` (Ellison p. 6 ``α_1``; our Takao ``eta``
        is exactly this for n=1, mining note §2). The impulse is applied at
        ``t_burn = alpha * tof_s`` after the leg start; the match point sits
        there.
    mu:
        Central-body gravitational parameter, km^3/s^2 (heliocentric default).
    """

    r0: Vec3
    v0: Vec3
    rf: Vec3
    vf: Vec3
    tof_s: float
    alpha: float
    mu: float = MU_SUN_KM3_S2

    def __post_init__(self) -> None:
        if self.tof_s <= 0.0:
            raise FbsMatchPointError(f"tof_s must be positive, got {self.tof_s}")
        if not (0.0 < self.alpha < 1.0):
            raise FbsMatchPointError(f"alpha must lie in (0, 1), got {self.alpha}")

    @property
    def t_burn_s(self) -> float:
        """Forward time from the leg start to the impulse, s (``alpha * tof_s``)."""
        return self.alpha * self.tof_s


def _forward_match_state(leg: FbsLeg, dv: Vec3) -> tuple[Vec6, Vec3, Vec3, Mat6]:
    """Forward state at the match point (just AFTER the impulse) + STM pieces.

    Returns ``(X_fwd_post (6,), r_fwd, v_fwd_pre, phi_fwd)`` where ``phi_fwd`` is
    the 6x6 STM of the forward coast ``propagate(r0, v0, t_burn)``. The impulse
    adds ``Δv`` to the velocity slot (massless MTM = identity + Eq. 42), which is
    a state shift only — it does not alter ``phi_fwd``.
    """
    r_fwd, v_fwd_pre, phi_fwd = shepperd_stm(leg.r0, leg.v0, leg.t_burn_s, leg.mu)
    x_post = np.empty(6, dtype=np.float64)
    x_post[0:3] = r_fwd
    x_post[3:6] = v_fwd_pre + np.asarray(dv, dtype=np.float64)
    return x_post, r_fwd, v_fwd_pre, phi_fwd


def _backward_match_state(leg: FbsLeg) -> tuple[Vec6, Mat6]:
    """Backward state at the match point + the backward-coast STM.

    Coasts ``(rf, vf)`` by ``-(tof_s - t_burn)`` to the impulse point. Returns
    ``(X_bwd (6,), phi_bwd)`` with ``phi_bwd = ∂X_bwd/∂(rf, vf)``.
    """
    dt_back = -(leg.tof_s - leg.t_burn_s)
    r_bwd, v_bwd, phi_bwd = shepperd_stm(leg.rf, leg.vf, dt_back, leg.mu)
    x_bwd = np.empty(6, dtype=np.float64)
    x_bwd[0:3] = r_bwd
    x_bwd[3:6] = v_bwd
    return x_bwd, phi_bwd


def match_point_defect(leg: FbsLeg, dv: Vec3) -> Vec6:
    r"""Match-point defect 6-vector ``c_mp = X^B - X^F`` (Ellison Eq. 2, massless).

    ``[r^B - r^F; v^B - v^F]`` — position mismatch (km, 3) and velocity mismatch
    (km/s, 3) between the backward- and forward-propagated states at the impulse
    point. Zero ⇒ the leg is dynamically consistent. No mass row (massless leg),
    no Lambert solve (Δv is the explicit decision variable).
    """
    arr = np.asarray(dv, dtype=np.float64)
    if arr.shape != (3,):
        raise FbsMatchPointError(f"dv must have shape (3,), got {arr.shape}")
    x_fwd, _, _, _ = _forward_match_state(leg, arr)
    x_bwd, _ = _backward_match_state(leg)
    return x_bwd - x_fwd
```

**Step 4: Run it to verify it passes**

Run: `uv run pytest tests/core/test_fbs_match_point.py::test_defect_zero_on_self_consistent_leg -v`
Expected: PASS.

**Step 5: Pre-commit checks + commit**

```bash
uv run ruff check src/cyclerfinder/core/fbs_match_point.py tests/core/test_fbs_match_point.py
uv run ruff format --check src/cyclerfinder/core/fbs_match_point.py tests/core/test_fbs_match_point.py
uv run mypy src/cyclerfinder/core/fbs_match_point.py
git status
git add src/cyclerfinder/core/fbs_match_point.py tests/core/test_fbs_match_point.py
git commit -m "core(fbs): single massless two-body DSM leg match-point defect (Ellison 2018 Eq.2, Path B #226)"
```

---

### Task 2: Analytic match-point defect Jacobian (FD-validated)

**Files:**
- Modify: `src/cyclerfinder/core/fbs_match_point.py` (add `match_point_defect_jacobian`)
- Test: `tests/core/test_fbs_match_point.py` (add the FD-vs-analytic consistency test)

**Concept — the Jacobian columns.** Decision variables for one leg, in this order, a 9-vector `x = [Δv (3); v0 (3); vf (3)]` (Phase 1 holds `r0, rf, tof_s, alpha` fixed; `v0`/`vf` stand in for the boundary v∞ slots — coordinate-free here, the v∞ → state map is an additive constant so `∂/∂v0 = ∂/∂v∞_out`). The defect is `c = X^B - X^F`. Then (Eq. 16 = `∂X^B/∂x - ∂X^F/∂x`):

- **`∂c/∂Δv`** (6x3): `X^F` depends on Δv only through the velocity slot (Eq. 42 connection `∂X^F/∂Δv = [0; I]`); `X^B` does not depend on Δv. So `∂c/∂Δv = -[0_{3x3}; I_{3x3}]` (shape 6x3, top block zero, bottom block `-I`).
- **`∂c/∂v0`** (6x3): only `X^F` depends on `v0`, through the forward coast STM. `∂X^F/∂v0 = phi_fwd[:, 3:6]` (the `[R; V]` columns = `dr/dv0`, `dv/dv0`). The impulse adds a constant, leaving this unchanged. So `∂c/∂v0 = -phi_fwd[:, 3:6]`.
- **`∂c/∂vf`** (6x3): only `X^B` depends on `vf`, through the backward coast STM. `∂X^B/∂vf = phi_bwd[:, 3:6]`. So `∂c/∂vf = +phi_bwd[:, 3:6]`.

Assemble into a 6x9 `J = [∂c/∂Δv | ∂c/∂v0 | ∂c/∂vf]`. This is the Eq. 31-32 chain for the trivial single-arc-per-side case: the chain is just one STM on each side, the connection matrices are `[0; I]` (Δv) / the STM velocity columns (boundary velocities).

**Step 1: Write the failing test (FD-vs-analytic consistency, several random states)**

```python
# central-difference step + tolerance, mirroring test_kepler_stm.py
_FD_STEP = 1.0e-6
_FD_RTOL = 1.0e-6


def _fd_defect_jacobian(leg, dv):
    """Central-difference 6x9 Jacobian of the defect w.r.t. [Δv; v0; vf]."""
    base_dv = np.asarray(dv, dtype=np.float64)
    cols = []
    # Δv columns
    for j in range(3):
        h = _FD_STEP * max(abs(base_dv[j]), 1.0)
        e = np.zeros(3); e[j] = h
        cp = match_point_defect(leg, base_dv + e)
        cm = match_point_defect(leg, base_dv - e)
        cols.append((cp - cm) / (2.0 * h))
    # v0 columns
    v0 = np.asarray(leg.v0, dtype=np.float64)
    for j in range(3):
        h = _FD_STEP * max(abs(v0[j]), float(np.linalg.norm(v0)))
        e = np.zeros(3); e[j] = h
        lp = FbsLeg(r0=leg.r0, v0=v0 + e, rf=leg.rf, vf=leg.vf, tof_s=leg.tof_s, alpha=leg.alpha, mu=leg.mu)
        lm = FbsLeg(r0=leg.r0, v0=v0 - e, rf=leg.rf, vf=leg.vf, tof_s=leg.tof_s, alpha=leg.alpha, mu=leg.mu)
        cols.append((match_point_defect(lp, base_dv) - match_point_defect(lm, base_dv)) / (2.0 * h))
    # vf columns
    vf = np.asarray(leg.vf, dtype=np.float64)
    for j in range(3):
        h = _FD_STEP * max(abs(vf[j]), float(np.linalg.norm(vf)))
        e = np.zeros(3); e[j] = h
        lp = FbsLeg(r0=leg.r0, v0=leg.v0, rf=leg.rf, vf=vf + e, tof_s=leg.tof_s, alpha=leg.alpha, mu=leg.mu)
        lm = FbsLeg(r0=leg.r0, v0=leg.v0, rf=leg.rf, vf=vf - e, tof_s=leg.tof_s, alpha=leg.alpha, mu=leg.mu)
        cols.append((match_point_defect(lp, base_dv) - match_point_defect(lm, base_dv)) / (2.0 * h))
    return np.column_stack(cols)


def _block_rel_err(a, b):
    """Block-wise (r/v rows) Frobenius-relative error; units differ by block."""
    errs = []
    for r0 in (0, 3):
        aa, bb = a[r0:r0 + 3, :], b[r0:r0 + 3, :]
        denom = float(np.linalg.norm(bb))
        if denom > 0.0:
            errs.append(float(np.linalg.norm(aa - bb)) / denom)
    return max(errs)


def test_jacobian_fd_consistency_random() -> None:
    """Analytic defect Jacobian matches central differences (consistency test)."""
    rng = np.random.default_rng(20260613)
    for _ in range(6):
        a_km = float(rng.uniform(0.6, 3.0)) * AU_KM
        e = float(rng.uniform(0.0, 0.6))
        nu = float(rng.uniform(0.0, 2.0 * pi))
        argp = float(rng.uniform(0.0, 2.0 * pi))
        r0, v0 = coe_to_rv(a_km, e, nu, arg_peri_rad=argp)
        period = 2.0 * pi * sqrt(a_km**3 / MU_SUN_KM3_S2)
        tof_s = float(rng.uniform(0.15, 0.6)) * period
        alpha = float(rng.uniform(0.25, 0.75))
        dv = rng.uniform(-0.3, 0.3, size=3)
        # arbitrary (NOT self-consistent) right boundary -> generic non-zero defect
        rf, vf = coe_to_rv(
            float(rng.uniform(0.6, 3.0)) * AU_KM, float(rng.uniform(0.0, 0.6)),
            float(rng.uniform(0.0, 2.0 * pi)), arg_peri_rad=float(rng.uniform(0.0, 2.0 * pi)),
        )
        leg = FbsLeg(r0=r0, v0=v0, rf=rf, vf=vf, tof_s=tof_s, alpha=alpha)
        j_an = match_point_defect_jacobian(leg, dv)
        j_fd = _fd_defect_jacobian(leg, dv)
        assert j_an.shape == (6, 9)
        assert _block_rel_err(j_an, j_fd) < _FD_RTOL, (a_km, e, nu, alpha, tof_s)
```

**Step 2: Run it to verify it fails**

Run: `uv run pytest tests/core/test_fbs_match_point.py::test_jacobian_fd_consistency_random -v`
Expected: FAIL — `match_point_defect_jacobian` not defined (the import at top of the test file already references it; this whole-file failure also confirms the symbol is missing).

**Step 3: Write minimal implementation**

```python
def match_point_defect_jacobian(leg: FbsLeg, dv: Vec3) -> NDArray[np.float64]:
    r"""Analytic 6x9 Jacobian of the match-point defect (Ellison Eqs. 16, 31-32).

    Columns are ordered ``[∂c/∂Δv (3) | ∂c/∂v0 (3) | ∂c/∂vf (3)]`` for the
    decision sub-vector ``x = [Δv; v0; vf]`` (Phase 1 holds ``r0, rf, tof_s,
    alpha`` fixed; ``v0``/``vf`` are the boundary-velocity / v∞ slots — the
    v∞ → state map is an additive constant, so ``∂c/∂v0 = ∂c/∂v∞_out`` etc.).

    For the massless single-DSM leg the Eq. 31-32 chain is one STM per side and
    the maneuver connection is the Eq. 42 velocity slot ``[0; I]``:

    * ``∂c/∂Δv  = -[0; I]``           (only ``X^F`` sees Δv, via Eq. 42)
    * ``∂c/∂v0  = -Φ_fwd[:, 3:6]``    (only ``X^F`` sees v0, via the fwd-coast STM)
    * ``∂c/∂vf  = +Φ_bwd[:, 3:6]``    (only ``X^B`` sees vf, via the bwd-coast STM)

    Validated FD-vs-analytic (consistency; Ellison publishes no numeric gradient).
    """
    arr = np.asarray(dv, dtype=np.float64)
    if arr.shape != (3,):
        raise FbsMatchPointError(f"dv must have shape (3,), got {arr.shape}")
    _, _, _, phi_fwd = _forward_match_state(leg, arr)
    _, phi_bwd = _backward_match_state(leg)

    jac = np.zeros((6, 9), dtype=np.float64)
    # ∂c/∂Δv = -[0; I]
    jac[3:6, 0:3] = -np.eye(3, dtype=np.float64)
    # ∂c/∂v0 = -Φ_fwd[:, 3:6]
    jac[:, 3:6] = -phi_fwd[:, 3:6]
    # ∂c/∂vf = +Φ_bwd[:, 3:6]
    jac[:, 6:9] = phi_bwd[:, 3:6]
    return jac
```

**Step 4: Run it to verify it passes**

Run: `uv run pytest tests/core/test_fbs_match_point.py -v`
Expected: PASS (both tests). If the random FD test trips on a near-degenerate draw, the issue is FD step scaling, not the analytic form — confirm the analytic vs FD blocks agree to ~1e-6 and only tighten/loosen `_FD_STEP` (NOT the analytic code) within the `cbrt(eps)`-justified band, mirroring `test_kepler_stm`.

**Step 5: Pre-commit checks + commit**

```bash
uv run ruff check src/cyclerfinder/core/fbs_match_point.py tests/core/test_fbs_match_point.py
uv run ruff format --check src/cyclerfinder/core/fbs_match_point.py tests/core/test_fbs_match_point.py
uv run mypy src/cyclerfinder/core/fbs_match_point.py
git status
git add src/cyclerfinder/core/fbs_match_point.py tests/core/test_fbs_match_point.py
git commit -m "core(fbs): analytic match-point defect Jacobian for single DSM leg, FD-validated (Ellison 2018 Eqs.31-32/42, #226)"
```

---

### Phase 1 exit criteria

- `tests/core/test_fbs_match_point.py` passes (zero-defect + FD-vs-analytic consistency, ~6 random states, rel tol 1e-6).
- `core/fbs_match_point.py` is additive: no other module imports change, no other lane touched. `uv run pytest` (full suite) green; ruff + mypy clean.
- **STOP. Phase 2 (Pitkin Eqs. 43-44 phase-TOF partials) is the next phase, started only after human review.**
