# #314 Heteroclinic-Cycle Framework Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a planar-CR3BP corrector that discovers and certifies *closed heteroclinic cycles* — chains O₁→O₂→…→O₁ of transversal invariant-manifold connections among equal-energy unstable orbits — as periodic-up-to-rotation transport loops, validated against the Wilczak-Zgliczyński Sun-Jupiter-Oterma L1↔L2 Lyapunov closed cycle.

**Architecture:** One new module `src/cyclerfinder/genome/heteroclinic_cycle.py`. It reuses (no edits) `core/cr3bp` (propagate + STM + Jacobi) and `search/cr3bp_periodic` (Lyapunov-orbit corrector + Floquet stability), and replicates the focused manifold-seeding pattern from `search/resonance_network` (Floquet eigenvector → ε-perturb → integrate to a Poincaré section). The connection primitive is a 2-D Newton on the {y=0} section gap between the unstable manifold of A and the stable manifold of B, with a finite-differenced Jacobian over the two manifold phases (τ_u, τ_s). The cycle assembler certifies a connection per consecutive pair and checks the chain returns to O₁. Closure = periodic-up-to-rotation (recurrence to the same orbit, phase free), NOT strict state(T)=state(0).

**Tech Stack:** Python 3, numpy, scipy.integrate.solve_ivp (DOP853 primary, Radau independent cross-check), pytest, uv, ruff, mypy. Spec: `docs/superpowers/specs/2026-06-19-314-heteroclinic-cycle-framework-design.md`.

---

## Conventions every task must follow

- **Run from the repo root** `/home/bruce/dev/cyclers` with `uv run` (uv-managed venv; never pip).
- **Pre-commit is mandatory.** Never `git commit --no-verify`. The hooks run ruff-check, ruff-format, mypy (`uv run mypy src tests`), and jsonschema. If a hook fails on a file you don't own, STOP and surface it — do not bypass.
- **Commit with explicit pathspecs** (`git add <path1> <path2>`), never `git add -A` (concurrent agents share this tree).
- **Do not push.** Commits stay local unless Bruce asks.
- **Sourced-golden discipline:** every EXPECTED value in a test must trace to the Wilczak-Zgliczyński publication (cited inline), never to a value our own code produced.
- **Type annotations required** (mypy is strict here): annotate every function signature and dataclass field.

## File structure

- **Create** `src/cyclerfinder/genome/heteroclinic_cycle.py` — the entire framework (dataclasses, Floquet eigenpair helper, manifold-section machinery, connection corrector, cycle assembler, independent cross-check). One module, one responsibility: certify heteroclinic connections and cycles in the planar CR3BP.
- **Create** `tests/genome/test_heteroclinic_cycle.py` — sourced-golden + self-consistency tests.

The W-Z golden constants are captured inline in the test module (cited to W-Z Part I, arXiv:math/0201278) so this plan is executable now; task #403 produces the fuller `data/golden/wz_oterma_heteroclinic.yaml` in parallel, and Task 8 wires the test to it when it lands.

## Reused signatures (verified — do not re-implement)

```python
# core/cr3bp.py
cr3bp.CR3BPSystem(mu, primary, secondary, l_km, t_s)            # frozen dataclass
cr3bp.propagate(system, state6, t, *, with_stm=False, rtol=1e-12, atol=1e-12, stm_mode="variable") -> CR3BPArc
#   CR3BPArc has .state_f (6,), .stm (6x6 or None), .t
cr3bp.cr3bp_eom(t, state6, mu) -> ndarray(6)
cr3bp.jacobi_constant(state6, mu) -> float

# search/cr3bp_periodic.py
correct_symmetric_fixed_jacobi(system, x0_guess, jacobi, period_guess, *, ydot0_sign=1.0,
    half_crossings=None, tol=1e-8, max_iter=30, rtol=1e-12, atol=1e-12,
    x0_bounds=(-2.0, 2.0)) -> SymmetricOrbit
#   SymmetricOrbit has .x0 .ydot0 .jacobi .t_half .period .converged .crossing_residual .n_iter
```

---

## Task 1: Module skeleton — dataclasses + Floquet eigenpair helper

**Files:**
- Create: `src/cyclerfinder/genome/heteroclinic_cycle.py`
- Test: `tests/genome/test_heteroclinic_cycle.py`

- [ ] **Step 1: Write the failing test**

Create `tests/genome/test_heteroclinic_cycle.py`:

```python
"""Tests for the #314 heteroclinic-cycle framework (planar CR3BP).

Sourced-golden discipline (feedback_golden_tests_sourced_only): EXPECTED values
trace to Wilczak & Zgliczyński, "Heteroclinic Connections between Periodic Orbits
in the Planar Restricted Three-Body Problem" Part I (arXiv:math/0201278, Comm.
Math. Phys.) — the computer-assisted proof of the closed L1<->L2 Lyapunov cycle
in the Sun-Jupiter-Oterma PCR3BP. Self-consistency checks (FD-Jacobian, empty-path)
need no external source, mirroring existing corrector tests.
"""

from __future__ import annotations

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.heteroclinic_cycle import (
    LyapunovNode,
    _planar_floquet_pair,
)

# --- W-Z Sun-Jupiter-Oterma golden (arXiv:math/0201278) ---------------------
WZ_MU = 0.0009537            # W-Z fixed Sun-Jupiter mass ratio (published exactly)
WZ_C = 3.03                  # Oterma Jacobi constant
# Lyapunov fixed points on the section {y=0}, params (x, xdot); xdot=0 at the
# perpendicular crossing. W-Z Part I, interval-enclosed centres:
WZ_X_L1 = 0.9208034913207400196
WZ_X_L2 = 1.081929486841799903


def _sun_jupiter() -> cr3bp.CR3BPSystem:
    # l_km / t_s are not used by the corrector math (all dynamics use mu only);
    # plausible Sun-Jupiter values for completeness.
    return cr3bp.CR3BPSystem(
        mu=WZ_MU, primary="sun", secondary="jupiter", l_km=778.57e6, t_s=5.957e8
    )


def test_floquet_pair_gives_unstable_and_stable_reciprocal() -> None:
    """A libration Lyapunov orbit has a real saddle Floquet pair (lambda, 1/lambda)."""
    system = _sun_jupiter()
    # Generate the L1 Lyapunov orbit at the Oterma energy (Task 2 wires the real
    # corrector; here we lean on the same primitive directly).
    node = LyapunovNode.from_libration(system, x0_guess=WZ_X_L1, jacobi=WZ_C,
                                       period_guess=3.0, label="L1")
    lam_u, v_u, lam_s, v_s = _planar_floquet_pair(system, node.state0, node.period)
    assert lam_u > 1.0 + 1e-3, f"unstable multiplier must exceed 1, got {lam_u}"
    assert lam_s < 1.0 - 1e-3, f"stable multiplier must be < 1, got {lam_s}"
    # Reciprocal saddle pair: lam_u * lam_s ~ 1.
    assert abs(lam_u * lam_s - 1.0) < 1e-2, f"not a reciprocal pair: {lam_u}*{lam_s}"
    assert v_u.shape == (4,) and v_s.shape == (4,)
    assert np.isclose(np.linalg.norm(v_u), 1.0) and np.isclose(np.linalg.norm(v_s), 1.0)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/genome/test_heteroclinic_cycle.py -x -q`
Expected: FAIL — `ImportError: cannot import name 'LyapunovNode'` (module doesn't exist).

- [ ] **Step 3: Write minimal implementation**

Create `src/cyclerfinder/genome/heteroclinic_cycle.py`:

```python
"""#314 Heteroclinic-cycle framework (planar CR3BP).

Discovers and certifies CLOSED HETEROCLINIC CYCLES: chains O_1 -> O_2 -> ... -> O_1
of transversal invariant-manifold connections among equal-energy unstable orbits.
This is a new closure definition -- periodic-up-to-rotation (recurrence to the same
orbit, phase along it free) -- distinct from the strict state(T)=state(0) periodicity
every other genome assumes.

Validated against Wilczak & Zgliczynski's computer-assisted proof of the closed
L1<->L2 Lyapunov cycle in the Sun-Jupiter-Oterma PCR3BP (arXiv:math/0201278). See
docs/superpowers/specs/2026-06-19-314-heteroclinic-cycle-framework-design.md.

Reuses core/cr3bp (propagate + STM + Jacobi) and search/cr3bp_periodic (Lyapunov
orbit corrector); replicates the focused Floquet manifold-seeding pattern from
search/resonance_network.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.cr3bp_periodic import correct_symmetric_fixed_jacobi

_PLANAR_IDX = [0, 1, 3, 4]  # (x, y, xdot, ydot) block of the 6x6 STM


@dataclass(frozen=True)
class LyapunovNode:
    """A libration-point Lyapunov orbit serving as a cycle node.

    ``state0`` is the full 6-vector IC ``(x0, 0, 0, 0, ydot0, 0)`` on the section
    {y=0}; ``unstable_eigvec`` / ``stable_eigvec`` are the planar 4-vectors
    (x, y, xdot, ydot) of the Floquet saddle pair, used to seed the manifolds.
    """

    label: str
    state0: NDArray[np.float64]
    period: float
    jacobi: float
    unstable_eigvec: NDArray[np.float64]
    stable_eigvec: NDArray[np.float64]
    converged: bool

    @classmethod
    def from_libration(
        cls,
        system: cr3bp.CR3BPSystem,
        *,
        x0_guess: float,
        jacobi: float,
        period_guess: float,
        label: str,
        ydot0_sign: float = 1.0,
        tol: float = 1e-10,
    ) -> LyapunovNode:
        """Correct a Lyapunov orbit at fixed Jacobi and extract its Floquet pair."""
        orbit = correct_symmetric_fixed_jacobi(
            system, x0_guess, jacobi, period_guess, ydot0_sign=ydot0_sign, tol=tol
        )
        state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0], dtype=np.float64)
        jac = cr3bp.jacobi_constant(state0, system.mu)
        _lu, v_u, _ls, v_s = _planar_floquet_pair(system, state0, orbit.period)
        return cls(
            label=label,
            state0=state0,
            period=orbit.period,
            jacobi=jac,
            unstable_eigvec=v_u,
            stable_eigvec=v_s,
            converged=orbit.converged,
        )


def _real_unit(v: NDArray[np.complex128]) -> NDArray[np.float64]:
    """Real-cast + unit-normalise a (possibly complex) eigenvector (4-vector)."""
    vr = np.real(v)
    n = float(np.linalg.norm(vr))
    if n < 1e-14:
        vr = np.real(v) + np.imag(v)
        n = float(np.linalg.norm(vr))
    return (vr / n).astype(np.float64) if n > 0.0 else vr.astype(np.float64)


def _planar_floquet_pair(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    period: float,
    *,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> tuple[float, NDArray[np.float64], float, NDArray[np.float64]]:
    """Return ``(|lam_u|, v_u, |lam_s|, v_s)`` for the planar monodromy saddle pair.

    Integrates the 6x6 STM over one period, slices the planar (x, y, xdot, ydot)
    block, and returns the largest- and smallest-magnitude eigenvalues with their
    real-normalised eigenvectors (the unstable and stable Floquet directions).
    Mirrors ``search/resonance_network._planar_floquet`` but returns BOTH ends of
    the reciprocal pair (the connection needs unstable-of-A and stable-of-B).
    """
    arc = cr3bp.propagate(system, state0, period, with_stm=True, rtol=rtol, atol=atol)
    assert arc.stm is not None
    phi4 = arc.stm[np.ix_(_PLANAR_IDX, _PLANAR_IDX)]
    eigvals, eigvecs = np.linalg.eig(phi4)
    mags = np.abs(eigvals)
    i_u = int(np.argmax(mags))
    i_s = int(np.argmin(mags))
    return (
        float(mags[i_u]),
        _real_unit(eigvecs[:, i_u]),
        float(mags[i_s]),
        _real_unit(eigvecs[:, i_s]),
    )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/genome/test_heteroclinic_cycle.py -x -q`
Expected: PASS (1 test).

- [ ] **Step 5: Lint + commit**

```bash
uv run ruff check src/cyclerfinder/genome/heteroclinic_cycle.py tests/genome/test_heteroclinic_cycle.py
uv run ruff format src/cyclerfinder/genome/heteroclinic_cycle.py tests/genome/test_heteroclinic_cycle.py
git add src/cyclerfinder/genome/heteroclinic_cycle.py tests/genome/test_heteroclinic_cycle.py
git commit -m "genome: #314 heteroclinic-cycle skeleton — LyapunovNode + Floquet pair"
```

---

## Task 2: Conventions self-check — reproduce the W-Z Lyapunov fixed points (validation-ladder item 1)

**Files:**
- Modify: `tests/genome/test_heteroclinic_cycle.py` (add test)

This is the cheap "do our μ/C/section conventions match the paper?" gate. No new source code — it exercises `LyapunovNode.from_libration` and asserts the corrected x0 equals W-Z's published fixed point.

- [ ] **Step 1: Write the failing test**

Append to `tests/genome/test_heteroclinic_cycle.py`:

```python
def test_lyapunov_fixed_points_match_wz() -> None:
    """Corrected L1/L2 Lyapunov x0 reproduce W-Z's section fixed points at C=3.03.

    EXPECTED = W-Z Part I interval-enclosed centres (arXiv:math/0201278); confirms
    our mu/Jacobi/section conventions agree with the paper before any connection.
    """
    system = _sun_jupiter()
    l1 = LyapunovNode.from_libration(system, x0_guess=WZ_X_L1, jacobi=WZ_C,
                                     period_guess=3.0, label="L1")
    l2 = LyapunovNode.from_libration(system, x0_guess=WZ_X_L2, jacobi=WZ_C,
                                     period_guess=3.0, label="L2")
    assert l1.converged and l2.converged
    # W-Z enclosures are ~1e-13; our corrector tol is 1e-10, so allow 1e-6.
    assert abs(l1.state0[0] - WZ_X_L1) < 1e-6, f"L1 x0={l1.state0[0]} vs {WZ_X_L1}"
    assert abs(l2.state0[0] - WZ_X_L2) < 1e-6, f"L2 x0={l2.state0[0]} vs {WZ_X_L2}"
    # Both nodes sit at the Oterma energy.
    assert abs(l1.jacobi - WZ_C) < 1e-6 and abs(l2.jacobi - WZ_C) < 1e-6
```

- [ ] **Step 2: Run test to verify it fails-or-passes meaningfully**

Run: `uv run pytest tests/genome/test_heteroclinic_cycle.py::test_lyapunov_fixed_points_match_wz -x -q`
Expected: This may PASS immediately (the implementation from Task 1 already supports it). If it FAILS on the x0 tolerance, the seed/period_guess needs adjustment: try `period_guess=3.4` and `ydot0_sign=-1.0` for L2. The corrector converges to the perpendicular-crossing Lyapunov orbit nearest the seed; W-Z's x* is that crossing. Document the working seeds in the test if you change them.

- [ ] **Step 3: (only if Step 2 failed) adjust seeds**

If convergence lands off the W-Z value, the cause is the `half_crossings` index. Pin it explicitly by passing through a `half_crossings=1` in `from_libration` (add the kwarg, default `None`, forward to `correct_symmetric_fixed_jacobi`). Show the change:

```python
    @classmethod
    def from_libration(
        cls,
        system: cr3bp.CR3BPSystem,
        *,
        x0_guess: float,
        jacobi: float,
        period_guess: float,
        label: str,
        ydot0_sign: float = 1.0,
        half_crossings: int | None = None,
        tol: float = 1e-10,
    ) -> LyapunovNode:
        orbit = correct_symmetric_fixed_jacobi(
            system, x0_guess, jacobi, period_guess, ydot0_sign=ydot0_sign,
            half_crossings=half_crossings, tol=tol,
        )
        # ... rest unchanged ...
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/genome/test_heteroclinic_cycle.py -x -q`
Expected: PASS (3 tests).

- [ ] **Step 5: Lint + commit**

```bash
uv run ruff check src/cyclerfinder/genome/heteroclinic_cycle.py tests/genome/test_heteroclinic_cycle.py
uv run ruff format src/cyclerfinder/genome/heteroclinic_cycle.py tests/genome/test_heteroclinic_cycle.py
git add src/cyclerfinder/genome/heteroclinic_cycle.py tests/genome/test_heteroclinic_cycle.py
git commit -m "genome: #314 validation-ladder item 1 — Lyapunov fixed points vs W-Z"
```

---

## Task 3: Manifold → section machinery

Globalize a manifold and record its {y=0} crossings. `_seed_on_manifold` transports the Floquet eigenvector to phase τ via the STM and ε-perturbs the orbit state; `_section_crossing` integrates the seed (forward = unstable, backward = stable) to the k-th qualifying section crossing and returns the section point (x, ẋ), or `None` if it never reaches the section (bounded — no hang, the #380 lesson).

**Files:**
- Modify: `src/cyclerfinder/genome/heteroclinic_cycle.py`
- Modify: `tests/genome/test_heteroclinic_cycle.py`

- [ ] **Step 1: Write the failing test**

Append to the test file:

```python
from cyclerfinder.genome.heteroclinic_cycle import (  # noqa: E402  (grouped import)
    _seed_on_manifold,
    _section_crossing,
)


def test_unstable_manifold_reaches_section() -> None:
    """The L1 unstable manifold crosses {y=0} within a bounded horizon."""
    system = _sun_jupiter()
    l1 = LyapunovNode.from_libration(system, x0_guess=WZ_X_L1, jacobi=WZ_C,
                                     period_guess=3.0, label="L1")
    seed = _seed_on_manifold(system, l1, tau=0.0, direction="unstable",
                             branch=+1, epsilon=1e-6)
    assert seed.shape == (6,)
    pt = _section_crossing(system, seed, direction="unstable", k=1,
                           max_time=8.0 * l1.period)
    assert pt is not None, "manifold must reach the {y=0} section"
    assert pt.shape == (2,)  # (x, xdot)


def test_section_miss_returns_none() -> None:
    """A horizon too short to reach the section yields None (no hang, no fabrication)."""
    system = _sun_jupiter()
    l1 = LyapunovNode.from_libration(system, x0_guess=WZ_X_L1, jacobi=WZ_C,
                                     period_guess=3.0, label="L1")
    seed = _seed_on_manifold(system, l1, tau=0.0, direction="unstable",
                             branch=+1, epsilon=1e-6)
    # Ask for the 9999th crossing — unreachable in this horizon.
    pt = _section_crossing(system, seed, direction="unstable", k=9999,
                           max_time=2.0 * l1.period)
    assert pt is None
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/genome/test_heteroclinic_cycle.py::test_unstable_manifold_reaches_section -x -q`
Expected: FAIL — `ImportError: cannot import name '_seed_on_manifold'`.

- [ ] **Step 3: Write the implementation**

Append to `src/cyclerfinder/genome/heteroclinic_cycle.py`:

```python
def _seed_on_manifold(
    system: cr3bp.CR3BPSystem,
    node: LyapunovNode,
    *,
    tau: float,
    direction: str,
    branch: int,
    epsilon: float,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> NDArray[np.float64]:
    """Manifold seed at phase ``tau`` along ``node``.

    Propagates the orbit to phase ``tau`` (state + STM), transports the chosen
    Floquet eigenvector with ``Phi(tau)``, renormalises, and ε-perturbs:
    ``state(tau) + branch * epsilon * v_hat(tau)``. ``direction`` selects the
    unstable (forward) or stable (backward) eigenvector.
    """
    if direction not in {"stable", "unstable"}:
        raise ValueError(f"direction must be 'stable' or 'unstable'; got {direction!r}")
    if branch not in (+1, -1):
        raise ValueError(f"branch must be +1 or -1; got {branch!r}")
    tau = float(tau) % float(node.period)
    if tau <= 0.0:
        state_tau = np.asarray(node.state0, float)
        phi = np.eye(6)
    else:
        arc = cr3bp.propagate(system, node.state0, tau, with_stm=True, rtol=rtol, atol=atol)
        assert arc.stm is not None
        state_tau = arc.state_f
        phi = arc.stm
    v4 = node.unstable_eigvec if direction == "unstable" else node.stable_eigvec
    v6 = np.array([v4[0], v4[1], 0.0, v4[2], v4[3], 0.0], dtype=np.float64)
    v_tau = phi @ v6
    n = float(np.linalg.norm(v_tau))
    if n > 0.0:
        v_tau = v_tau / n
    return (state_tau + float(branch) * float(epsilon) * v_tau).astype(np.float64)


def _section_crossing(
    system: cr3bp.CR3BPSystem,
    seed: NDArray[np.float64],
    *,
    direction: str,
    k: int,
    max_time: float,
    section_y: float = 0.0,
    ydot_sign: int | None = None,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    method: str = "DOP853",
) -> NDArray[np.float64] | None:
    """Integrate ``seed`` to the k-th qualifying crossing of {y=section_y}.

    Forward in time for ``direction="unstable"``, backward for ``"stable"``.
    ``ydot_sign`` (if given) restricts to the Theta+/Theta- half of the section
    (sign of ydot at the crossing), matching W-Z's split. Returns the section
    point ``(x, xdot)`` at the k-th crossing (1-based), or ``None`` if fewer than
    ``k`` qualifying crossings occur within ``max_time`` (bounded — never hangs,
    never fabricates a crossing).
    """
    if direction not in {"stable", "unstable"}:
        raise ValueError(f"direction must be 'stable' or 'unstable'; got {direction!r}")
    horizon = abs(float(max_time))
    t_span = (0.0, horizon) if direction == "unstable" else (0.0, -horizon)

    def _y_event(t: float, y: NDArray[np.float64], _mu: float) -> float:
        return float(y[1] - section_y)

    _y_event.terminal = False  # type: ignore[attr-defined]
    _y_event.direction = 0.0  # type: ignore[attr-defined]

    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        t_span,
        np.asarray(seed, float),
        args=(system.mu,),  # type: ignore[call-overload]
        method=method,
        rtol=rtol,
        atol=atol,
        events=_y_event,
        max_step=horizon / 500.0,  # fine enough that no {y=0} crossing is stepped over
    )
    t_events = sol.t_events[0] if sol.t_events is not None else np.array([])
    y_events = sol.y_events[0] if sol.y_events is not None else []
    t_floor = 1e-6 * horizon
    count = 0
    for t_ev, y_ev in zip(t_events, y_events, strict=False):
        if abs(float(t_ev)) <= t_floor:
            continue  # skip the t~0 root at the seed
        if ydot_sign is not None and int(np.sign(float(y_ev[4]))) != ydot_sign:
            continue
        count += 1
        if count == k:
            return np.array([float(y_ev[0]), float(y_ev[3])], dtype=np.float64)
    return None
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/genome/test_heteroclinic_cycle.py -x -q`
Expected: PASS (5 tests).

- [ ] **Step 5: Lint + commit**

```bash
uv run ruff check src/cyclerfinder/genome/heteroclinic_cycle.py tests/genome/test_heteroclinic_cycle.py
uv run ruff format src/cyclerfinder/genome/heteroclinic_cycle.py tests/genome/test_heteroclinic_cycle.py
git add src/cyclerfinder/genome/heteroclinic_cycle.py tests/genome/test_heteroclinic_cycle.py
git commit -m "genome: #314 manifold-to-section machinery (seed + bounded {y=0} crossing)"
```

---

## Task 4: The connection corrector — 2-D Newton on the section gap

`HeteroclinicConnection` + `correct_connection`. Free variables (τ_u, τ_s); residual = (x, ẋ) gap on {y=0} between the k_u-th crossing of Wu(A) and the k_s-th crossing of Ws(B); FD-Jacobian Newton with backtracking; equal-energy guard.

**Files:**
- Modify: `src/cyclerfinder/genome/heteroclinic_cycle.py`
- Modify: `tests/genome/test_heteroclinic_cycle.py`

- [ ] **Step 1: Write the failing test**

Append to the test file:

```python
from cyclerfinder.genome.heteroclinic_cycle import (  # noqa: E402
    HeteroclinicConnection,
    correct_connection,
)

import pytest  # noqa: E402


def test_connection_energy_mismatch_raises() -> None:
    """Connections require equal Jacobi; a mismatch is a hard error."""
    system = _sun_jupiter()
    l1 = LyapunovNode.from_libration(system, x0_guess=WZ_X_L1, jacobi=WZ_C,
                                     period_guess=3.0, label="L1")
    l2 = LyapunovNode.from_libration(system, x0_guess=WZ_X_L2, jacobi=WZ_C + 0.05,
                                     period_guess=3.0, label="L2")
    with pytest.raises(ValueError, match="equal Jacobi|energy"):
        correct_connection(system, l1, l2)


def test_connection_l1_to_l2_converges() -> None:
    """Wu(L1) meets Ws(L2) on {y=0}: a transversal heteroclinic connection.

    W-Z Part I proves this connection exists; we certify the section-gap residual
    closes. (The exact crossing-coordinate cross-check is Task 8, gated on #403's
    full golden; here we assert convergence + that the meeting point lies between
    the two libration x positions, i.e. in the L1-L2 neck.)
    """
    system = _sun_jupiter()
    l1 = LyapunovNode.from_libration(system, x0_guess=WZ_X_L1, jacobi=WZ_C,
                                     period_guess=3.0, label="L1")
    l2 = LyapunovNode.from_libration(system, x0_guess=WZ_X_L2, jacobi=WZ_C,
                                     period_guess=3.0, label="L2")
    conn = correct_connection(system, l1, l2, tol=1e-7)
    assert isinstance(conn, HeteroclinicConnection)
    assert conn.converged, f"residual={conn.residual:.3e}, n_iter={conn.n_iter}"
    assert conn.residual < 1e-6
    assert WZ_X_L1 - 0.1 < conn.crossing_xv[0] < WZ_X_L2 + 0.1
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/genome/test_heteroclinic_cycle.py::test_connection_l1_to_l2_converges -x -q`
Expected: FAIL — `ImportError: cannot import name 'correct_connection'`.

- [ ] **Step 3: Write the implementation**

Append to `src/cyclerfinder/genome/heteroclinic_cycle.py`:

```python
@dataclass(frozen=True)
class HeteroclinicConnection:
    """A certified (or attempted) Wu(A) -> Ws(B) connection on {y=0}."""

    orbit_from: str
    orbit_to: str
    jacobi: float
    tau_u: float
    tau_s: float
    k_u: int
    k_s: int
    crossing_xv: NDArray[np.float64]  # (x, xdot) at the matched crossing
    residual: float
    converged: bool
    n_iter: int
    notes: str = ""


def _connection_residual(
    system: cr3bp.CR3BPSystem,
    a: LyapunovNode,
    b: LyapunovNode,
    *,
    tau_u: float,
    tau_s: float,
    k_u: int,
    k_s: int,
    epsilon: float,
    branch_u: int,
    branch_s: int,
    max_time: float,
    ydot_sign_u: int | None,
    ydot_sign_s: int | None,
    method: str = "DOP853",
) -> tuple[NDArray[np.float64] | None, NDArray[np.float64] | None]:
    """Return ``(residual2, crossing_xv)`` or ``(None, None)`` if a leg misses Σ."""
    seed_u = _seed_on_manifold(system, a, tau=tau_u, direction="unstable",
                               branch=branch_u, epsilon=epsilon)
    p_u = _section_crossing(system, seed_u, direction="unstable", k=k_u,
                            max_time=max_time, ydot_sign=ydot_sign_u, method=method)
    seed_s = _seed_on_manifold(system, b, tau=tau_s, direction="stable",
                               branch=branch_s, epsilon=epsilon)
    p_s = _section_crossing(system, seed_s, direction="stable", k=k_s,
                            max_time=max_time, ydot_sign=ydot_sign_s, method=method)
    if p_u is None or p_s is None:
        return None, None
    return (p_u - p_s).astype(np.float64), p_u


def correct_connection(
    system: cr3bp.CR3BPSystem,
    orbit_from: LyapunovNode,
    orbit_to: LyapunovNode,
    *,
    k_u: int = 1,
    k_s: int = 1,
    epsilon: float = 1e-6,
    branch_u: int = +1,
    branch_s: int = -1,
    tau_u0: float | None = None,
    tau_s0: float | None = None,
    ydot_sign_u: int | None = None,
    ydot_sign_s: int | None = None,
    tol: float = 1e-7,
    max_iter: int = 40,
    fd_step: float = 1e-6,
    max_time_factor: float = 8.0,
    jacobi_tol: float = 1e-6,
) -> HeteroclinicConnection:
    """Certify Wu(orbit_from) ∩ Ws(orbit_to) on {y=0} by 2-D Newton on (tau_u, tau_s).

    Residual = section-plane gap (Δx, Δxdot) between the k_u-th crossing of the
    unstable manifold of ``orbit_from`` and the k_s-th crossing of the stable
    manifold of ``orbit_to``. Free variables = the manifold phases (tau_u, tau_s).
    Jacobian is finite-differenced (2x2; cheap, robust); Newton step damped by
    backtracking on ||residual||. Raises ``ValueError`` on an energy mismatch
    (a connection cannot exist off equal Jacobi). A leg that never reaches Σ ->
    ``converged=False`` with a diagnostic note (never a fabricated closure).
    """
    if abs(orbit_from.jacobi - orbit_to.jacobi) > jacobi_tol:
        raise ValueError(
            f"connection requires equal Jacobi (energy): "
            f"{orbit_from.label} C={orbit_from.jacobi:.6f} vs "
            f"{orbit_to.label} C={orbit_to.jacobi:.6f}"
        )
    max_time = max_time_factor * max(orbit_from.period, orbit_to.period)
    tau_u = 0.25 * orbit_from.period if tau_u0 is None else float(tau_u0)
    tau_s = 0.25 * orbit_to.period if tau_s0 is None else float(tau_s0)

    def _resid(tu: float, ts: float) -> tuple[NDArray[np.float64] | None, NDArray[np.float64] | None]:
        return _connection_residual(
            system, orbit_from, orbit_to, tau_u=tu, tau_s=ts, k_u=k_u, k_s=k_s,
            epsilon=epsilon, branch_u=branch_u, branch_s=branch_s, max_time=max_time,
            ydot_sign_u=ydot_sign_u, ydot_sign_s=ydot_sign_s,
        )

    res, xv = _resid(tau_u, tau_s)
    n_iter = 0
    for n_iter in range(1, max_iter + 1):  # noqa: B007 -- returned as iteration count
        if res is None:
            return HeteroclinicConnection(
                orbit_from.label, orbit_to.label, orbit_from.jacobi, tau_u, tau_s,
                k_u, k_s, np.zeros(2), float("inf"), False, n_iter,
                notes="manifold leg did not reach the section",
            )
        rn = float(np.linalg.norm(res))
        if rn < tol:
            break
        # 2x2 finite-difference Jacobian d(residual)/d(tau_u, tau_s).
        jac = np.zeros((2, 2), dtype=np.float64)
        ok = True
        for j, (du, ds) in enumerate([(fd_step, 0.0), (0.0, fd_step)]):
            rp, _ = _resid(tau_u + du, tau_s + ds)
            if rp is None:
                ok = False
                break
            jac[:, j] = (rp - res) / fd_step
        if not ok:
            return HeteroclinicConnection(
                orbit_from.label, orbit_to.label, orbit_from.jacobi, tau_u, tau_s,
                k_u, k_s, xv if xv is not None else np.zeros(2), rn, False, n_iter,
                notes="FD-Jacobian probe left the manifold's section branch",
            )
        try:
            step = np.linalg.solve(jac, -res)
        except np.linalg.LinAlgError:
            step, *_ = np.linalg.lstsq(jac, -res, rcond=None)
        # Backtracking line search on ||residual||.
        alpha = 1.0
        improved = False
        for _ in range(20):
            tu_t = (tau_u + alpha * float(step[0])) % orbit_from.period
            ts_t = (tau_s + alpha * float(step[1])) % orbit_to.period
            res_t, xv_t = _resid(tu_t, ts_t)
            if res_t is not None and float(np.linalg.norm(res_t)) < rn:
                tau_u, tau_s, res, xv = tu_t, ts_t, res_t, xv_t
                improved = True
                break
            alpha *= 0.5
        if not improved:
            break  # cannot make progress; report the best so far
    final_rn = float(np.linalg.norm(res)) if res is not None else float("inf")
    converged = res is not None and final_rn < tol
    return HeteroclinicConnection(
        orbit_from=orbit_from.label,
        orbit_to=orbit_to.label,
        jacobi=orbit_from.jacobi,
        tau_u=tau_u,
        tau_s=tau_s,
        k_u=k_u,
        k_s=k_s,
        crossing_xv=xv if xv is not None else np.zeros(2),
        residual=final_rn,
        converged=converged,
        n_iter=n_iter,
        notes="" if converged else "did not reach tol",
    )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/genome/test_heteroclinic_cycle.py -x -q`
Expected: PASS (7 tests). If `test_connection_l1_to_l2_converges` fails to converge, widen the search: try `k_u=2`/`k_s=2`, or set `tau_u0`/`tau_s0` to scan a few starts (the connection lives at a specific manifold phase; the Newton needs a start in its basin). Record the working `(k_u, k_s, tau_u0, tau_s0, branch_u, branch_s)` as defaults/constants in the test. The W-Z neck crossing is at x≈0.921 / x≈1.0817 with small ẋ≈5e-4 — use that to sanity-check you're on the right branch.

- [ ] **Step 5: Lint + commit**

```bash
uv run ruff check src/cyclerfinder/genome/heteroclinic_cycle.py tests/genome/test_heteroclinic_cycle.py
uv run ruff format src/cyclerfinder/genome/heteroclinic_cycle.py tests/genome/test_heteroclinic_cycle.py
git add src/cyclerfinder/genome/heteroclinic_cycle.py tests/genome/test_heteroclinic_cycle.py
git commit -m "genome: #314 connection corrector — 2-D Newton on the {y=0} section gap"
```

---

## Task 5: Cycle assembly + the periodic-up-to-rotation closure definition

`HeteroclinicCycle` + `assemble_cycle`. Certify a connection per consecutive pair (chain wraps back to O₁) and set `closed` iff every leg converged and every node is at the shared Jacobi.

**Files:**
- Modify: `src/cyclerfinder/genome/heteroclinic_cycle.py`
- Modify: `tests/genome/test_heteroclinic_cycle.py`

- [ ] **Step 1: Write the failing test**

Append to the test file:

```python
from cyclerfinder.genome.heteroclinic_cycle import (  # noqa: E402
    HeteroclinicCycle,
    assemble_cycle,
)


def test_assemble_l1_l2_two_cycle_closes() -> None:
    """The L1->L2->L1 chain forms a closed heteroclinic cycle (W-Z, both directions)."""
    system = _sun_jupiter()
    l1 = LyapunovNode.from_libration(system, x0_guess=WZ_X_L1, jacobi=WZ_C,
                                     period_guess=3.0, label="L1")
    l2 = LyapunovNode.from_libration(system, x0_guess=WZ_X_L2, jacobi=WZ_C,
                                     period_guess=3.0, label="L2")
    cycle = assemble_cycle(system, [l1, l2], tol=1e-7)
    assert isinstance(cycle, HeteroclinicCycle)
    assert cycle.closed, (
        f"max_leg_residual={cycle.max_leg_residual:.3e}, "
        f"symbols={cycle.symbol_sequence}"
    )
    assert len(cycle.connections) == 2  # L1->L2 and L2->L1
    assert cycle.symbol_sequence == ["L1", "L2", "L1"]
    assert abs(cycle.jacobi - WZ_C) < 1e-6


def test_assemble_energy_mismatch_raises() -> None:
    """A node off the shared energy is rejected before any leg is attempted."""
    system = _sun_jupiter()
    l1 = LyapunovNode.from_libration(system, x0_guess=WZ_X_L1, jacobi=WZ_C,
                                     period_guess=3.0, label="L1")
    l2 = LyapunovNode.from_libration(system, x0_guess=WZ_X_L2, jacobi=WZ_C + 0.05,
                                     period_guess=3.0, label="L2")
    with pytest.raises(ValueError, match="equal Jacobi|energy"):
        assemble_cycle(system, [l1, l2])
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/genome/test_heteroclinic_cycle.py::test_assemble_l1_l2_two_cycle_closes -x -q`
Expected: FAIL — `ImportError: cannot import name 'assemble_cycle'`.

- [ ] **Step 3: Write the implementation**

Append to `src/cyclerfinder/genome/heteroclinic_cycle.py`:

```python
@dataclass(frozen=True)
class HeteroclinicCycle:
    """A closed chain of heteroclinic connections (periodic-up-to-rotation).

    ``closed`` means the itinerary returns to the same first orbit O_1 with every
    leg certified at the shared Jacobi -- a recurrence, NOT strict state(T)=state(0)
    periodicity. ``symbol_sequence`` is the node-label itinerary including the
    return to O_1 (e.g. ["L1", "L2", "L1"]).
    """

    orbits: list[str]
    connections: list[HeteroclinicConnection]
    jacobi: float
    closed: bool
    max_leg_residual: float
    independent_residual: float
    symbol_sequence: list[str]
    notes: str = ""


def assemble_cycle(
    system: cr3bp.CR3BPSystem,
    nodes: list[LyapunovNode],
    *,
    tol: float = 1e-7,
    jacobi_tol: float = 1e-6,
    connection_kwargs: dict[str, object] | None = None,
) -> HeteroclinicCycle:
    """Certify a closed heteroclinic cycle over the ordered ``nodes`` (wraps to O_1).

    Requires >= 1 node; a single node is a (degenerate) homoclinic cycle O->O.
    Raises ``ValueError`` if the nodes are not all at the same Jacobi. ``closed``
    iff every consecutive-pair connection converged. ``independent_residual`` is
    filled by Task 6 (set to nan here as a clear "not yet cross-checked" sentinel).
    """
    if not nodes:
        raise ValueError("assemble_cycle requires at least one node")
    c0 = nodes[0].jacobi
    for nd in nodes[1:]:
        if abs(nd.jacobi - c0) > jacobi_tol:
            raise ValueError(
                f"all cycle nodes must share Jacobi (energy): {nd.label} "
                f"C={nd.jacobi:.6f} vs {nodes[0].label} C={c0:.6f}"
            )
    kw = dict(connection_kwargs or {})
    kw.setdefault("tol", tol)
    n = len(nodes)
    connections: list[HeteroclinicConnection] = []
    for i in range(n):
        a = nodes[i]
        b = nodes[(i + 1) % n]
        conn = correct_connection(system, a, b, **kw)  # type: ignore[arg-type]
        connections.append(conn)
    max_leg = max((c.residual for c in connections), default=float("inf"))
    closed = all(c.converged for c in connections)
    symbols = [nd.label for nd in nodes] + [nodes[0].label]
    return HeteroclinicCycle(
        orbits=[nd.label for nd in nodes],
        connections=connections,
        jacobi=c0,
        closed=closed,
        max_leg_residual=max_leg,
        independent_residual=float("nan"),
        symbol_sequence=symbols,
        notes="" if closed else "one or more legs did not converge",
    )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/genome/test_heteroclinic_cycle.py -x -q`
Expected: PASS (9 tests). If the L2→L1 leg needs different `(k, branch, tau0)` than L1→L2, pass them per-leg: extend `assemble_cycle` to accept a `per_leg_kwargs: list[dict] | None` and use `per_leg_kwargs[i]` when present. Add that only if a single global `connection_kwargs` cannot close both legs.

- [ ] **Step 5: Lint + commit**

```bash
uv run ruff check src/cyclerfinder/genome/heteroclinic_cycle.py tests/genome/test_heteroclinic_cycle.py
uv run ruff format src/cyclerfinder/genome/heteroclinic_cycle.py tests/genome/test_heteroclinic_cycle.py
git add src/cyclerfinder/genome/heteroclinic_cycle.py tests/genome/test_heteroclinic_cycle.py
git commit -m "genome: #314 cycle assembler + periodic-up-to-rotation closure"
```

---

## Task 6: Independent Radau cross-check (validation-ladder item 4)

Re-derive each converged leg's section crossing with the Radau integrator (vs the corrector's DOP853) and store the max disagreement in `HeteroclinicCycle.independent_residual`. Mandatory before any "closed" claim is trusted.

**Files:**
- Modify: `src/cyclerfinder/genome/heteroclinic_cycle.py`
- Modify: `tests/genome/test_heteroclinic_cycle.py`

- [ ] **Step 1: Write the failing test**

Append to the test file:

```python
from cyclerfinder.genome.heteroclinic_cycle import (  # noqa: E402
    crosscheck_cycle,
)


def test_cycle_independent_crosscheck() -> None:
    """Radau re-propagation reproduces each leg's section crossing (DOP853)."""
    system = _sun_jupiter()
    l1 = LyapunovNode.from_libration(system, x0_guess=WZ_X_L1, jacobi=WZ_C,
                                     period_guess=3.0, label="L1")
    l2 = LyapunovNode.from_libration(system, x0_guess=WZ_X_L2, jacobi=WZ_C,
                                     period_guess=3.0, label="L2")
    cycle = assemble_cycle(system, [l1, l2], tol=1e-7)
    assert cycle.closed
    checked = crosscheck_cycle(system, [l1, l2], cycle)
    assert checked.independent_residual < 1e-5, (
        f"Radau vs DOP853 disagreement {checked.independent_residual:.3e}"
    )
    # The cross-check returns a NEW cycle object with the residual filled.
    assert not np.isnan(checked.independent_residual)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/genome/test_heteroclinic_cycle.py::test_cycle_independent_crosscheck -x -q`
Expected: FAIL — `ImportError: cannot import name 'crosscheck_cycle'`.

- [ ] **Step 3: Write the implementation**

Append to `src/cyclerfinder/genome/heteroclinic_cycle.py`:

```python
def crosscheck_cycle(
    system: cr3bp.CR3BPSystem,
    nodes: list[LyapunovNode],
    cycle: HeteroclinicCycle,
    *,
    method: str = "Radau",
    rtol: float = 1e-11,
    atol: float = 1e-11,
    max_time_factor: float = 8.0,
) -> HeteroclinicCycle:
    """Re-derive each converged leg's crossing with an independent integrator.

    Returns a copy of ``cycle`` with ``independent_residual`` set to the maximum
    section-point disagreement between the corrector's stored crossing (DOP853) and
    a fresh ``method`` (default Radau) re-propagation at the same (tau_u, tau_s).
    A leg that fails to reproduce its crossing contributes ``inf`` (a real failure,
    surfaced — never silently dropped).
    """
    n = len(nodes)
    label_to_node = {nd.label: nd for nd in nodes}
    worst = 0.0
    for i, conn in enumerate(cycle.connections):
        if not conn.converged:
            worst = float("inf")
            continue
        a = label_to_node[conn.orbit_from]
        b = label_to_node[conn.orbit_to]
        max_time = max_time_factor * max(a.period, b.period)
        seed_u = _seed_on_manifold(system, a, tau=conn.tau_u, direction="unstable",
                                   branch=+1, epsilon=1e-6)
        seed_s = _seed_on_manifold(system, b, tau=conn.tau_s, direction="stable",
                                   branch=-1, epsilon=1e-6)
        p_u = _section_crossing(system, seed_u, direction="unstable", k=conn.k_u,
                                max_time=max_time, method=method, rtol=rtol, atol=atol)
        p_s = _section_crossing(system, seed_s, direction="stable", k=conn.k_s,
                                max_time=max_time, method=method, rtol=rtol, atol=atol)
        if p_u is None or p_s is None:
            worst = float("inf")
            continue
        gap = float(np.linalg.norm(p_u - p_s))  # independent section-gap residual
        # Compare the independently-derived meeting point to the stored one.
        recovered = 0.5 * (p_u + p_s)
        worst = max(worst, float(np.linalg.norm(recovered - conn.crossing_xv)), gap)
    return HeteroclinicCycle(
        orbits=cycle.orbits,
        connections=cycle.connections,
        jacobi=cycle.jacobi,
        closed=cycle.closed,
        max_leg_residual=cycle.max_leg_residual,
        independent_residual=worst,
        symbol_sequence=cycle.symbol_sequence,
        notes=cycle.notes,
    )
```

> Note: `crosscheck_cycle` uses the default `branch_u=+1`/`branch_s=-1`. If Task 4/5 settled on different branches to close the legs, thread the per-leg branches through `HeteroclinicConnection` (add `branch_u: int` / `branch_s: int` fields, populate them in `correct_connection`, and read them here) so the cross-check seeds the same manifold side. Do this if the test disagreement is large (~O(1)) rather than ~1e-6.

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/genome/test_heteroclinic_cycle.py -x -q`
Expected: PASS (10 tests).

- [ ] **Step 5: Lint + commit**

```bash
uv run ruff check src/cyclerfinder/genome/heteroclinic_cycle.py tests/genome/test_heteroclinic_cycle.py
uv run ruff format src/cyclerfinder/genome/heteroclinic_cycle.py tests/genome/test_heteroclinic_cycle.py
git add src/cyclerfinder/genome/heteroclinic_cycle.py tests/genome/test_heteroclinic_cycle.py
git commit -m "genome: #314 independent Radau cross-check (validation-ladder item 4)"
```

---

## Task 7: Edge cases — non-transversal connection reports a clean negative

Confirm the framework reports `converged=False` (never fabricates closure) when no connection exists, and that a non-closing chain yields `closed=False`. This locks in the "a clean negative is an acceptable outcome" discipline.

**Files:**
- Modify: `tests/genome/test_heteroclinic_cycle.py`

- [ ] **Step 1: Write the test**

Append to the test file:

```python
def test_no_connection_reports_clean_negative() -> None:
    """A too-short horizon -> the legs never meet -> converged=False, no exception."""
    system = _sun_jupiter()
    l1 = LyapunovNode.from_libration(system, x0_guess=WZ_X_L1, jacobi=WZ_C,
                                     period_guess=3.0, label="L1")
    l2 = LyapunovNode.from_libration(system, x0_guess=WZ_X_L2, jacobi=WZ_C,
                                     period_guess=3.0, label="L2")
    # Starve the integration horizon so no qualifying crossing is reached.
    conn = correct_connection(system, l1, l2, max_time_factor=0.05, max_iter=5)
    assert not conn.converged
    assert conn.residual == float("inf") or conn.residual > 1e-6
    assert conn.notes  # a diagnostic is recorded


def test_nonclosing_chain_is_not_closed() -> None:
    """If a leg cannot close, the cycle is reported open (not silently 'closed')."""
    system = _sun_jupiter()
    l1 = LyapunovNode.from_libration(system, x0_guess=WZ_X_L1, jacobi=WZ_C,
                                     period_guess=3.0, label="L1")
    l2 = LyapunovNode.from_libration(system, x0_guess=WZ_X_L2, jacobi=WZ_C,
                                     period_guess=3.0, label="L2")
    cycle = assemble_cycle(system, [l1, l2],
                           connection_kwargs={"max_time_factor": 0.05, "max_iter": 5})
    assert not cycle.closed
    assert cycle.notes
```

- [ ] **Step 2: Run tests to verify they pass**

Run: `uv run pytest tests/genome/test_heteroclinic_cycle.py -x -q`
Expected: PASS (12 tests). These should pass against the existing implementation (the not-reached path already returns `converged=False`). If `correct_connection` raises instead of returning a non-converged result on a starved horizon, fix the implementation so the no-crossing path returns the diagnostic dataclass (it must never raise for a legitimate "no connection here").

- [ ] **Step 3: Full-module test + type-check + commit**

```bash
uv run pytest tests/genome/test_heteroclinic_cycle.py -q
uv run mypy src/cyclerfinder/genome/heteroclinic_cycle.py
uv run ruff check src/cyclerfinder/genome/heteroclinic_cycle.py tests/genome/test_heteroclinic_cycle.py
uv run ruff format --check src/cyclerfinder/genome/heteroclinic_cycle.py tests/genome/test_heteroclinic_cycle.py
git add tests/genome/test_heteroclinic_cycle.py
git commit -m "genome: #314 edge cases — clean negatives never fabricate closure"
```

Expected: all 12 tests PASS; mypy clean; ruff clean.

---

## Task 8: Wire the test to #403's golden + per-leg crossing assertions

When task #403 lands `data/golden/wz_oterma_heteroclinic.yaml` (the full per-leg crossing coordinates), strengthen the connection test from "lands in the neck" to "matches the published crossing (x, ẋ)". **If #403 is not yet merged when you reach this task, SKIP it and leave a note in the final summary — do not block the build.**

**Files:**
- Modify: `tests/genome/test_heteroclinic_cycle.py`

- [ ] **Step 1: Check the golden exists**

Run: `ls data/golden/wz_oterma_heteroclinic.yaml`
If absent: skip this task; the inline W-Z constants already provide a sourced item-1 + neck check. Note the skip in the handoff.

- [ ] **Step 2: Write the golden-backed test**

The golden schema (verified — file already exists, committed in `cd659cd`):

```yaml
system: {mu: 0.0009537, C: 3.03}
lyapunov_fixed_points:
  L1_star: {point: [0.9208034913207400196, 0.0], eigenvector_unstable: [1.0, 2.5733011], ...}
  L2_star: {point: [1.081929486841799903, 0.0], ...}
crossings:
  heteroclinic_L1_to_L2:
    source: "Part I, Section 6.1"
    sequence:           # the connecting orbit's successive {y=0} crossings, each [x, xdot]
      - [0.9522928423486199945, 1.23e-05]
      - [0.921005737890425169, 0.0005205932817646883714]
      # ... 8 points total
  homoclinic_L2_exterior_4_symbol: {...}
  homoclinic_L1_interior_4_symbol: {...}
  homoclinic_L2_2_3_resonance_6_symbol: {...}
  homoclinic_L1_5_3_resonance_6_symbol: {...}
```

The corrector certifies ONE matched crossing per leg (`conn.crossing_xv`), which must coincide with one of the published sequence points for that connection. Append:

```python
import pathlib  # noqa: E402

import yaml  # noqa: E402

_GOLDEN = pathlib.Path("data/golden/wz_oterma_heteroclinic.yaml")


def test_l1_to_l2_crossing_matches_wz_golden() -> None:
    """The certified L1->L2 crossing (x, xdot) matches one published W-Z crossing."""
    data = yaml.safe_load(_GOLDEN.read_text())
    system = _sun_jupiter()
    l1 = LyapunovNode.from_libration(system, x0_guess=WZ_X_L1, jacobi=WZ_C,
                                     period_guess=3.0, label="L1")
    l2 = LyapunovNode.from_libration(system, x0_guess=WZ_X_L2, jacobi=WZ_C,
                                     period_guess=3.0, label="L2")
    conn = correct_connection(system, l1, l2, tol=1e-8)
    assert conn.converged
    seq = np.array(data["crossings"]["heteroclinic_L1_to_L2"]["sequence"], dtype=np.float64)
    # The corrector's matched crossing must be one of W-Z's tabulated section points.
    dists = np.linalg.norm(seq - conn.crossing_xv[None, :], axis=1)
    assert float(dists.min()) < 1e-4, (
        f"crossing {conn.crossing_xv} not among W-Z sequence (min dist {dists.min():.3e})"
    )
```

Also verify the inline fixed-point constants `WZ_X_L1`/`WZ_X_L2` equal `data["lyapunov_fixed_points"]["L1_star"]["point"][0]` / `L2_star` — they should match exactly (same source).

- [ ] **Step 3: Run + commit**

```bash
uv run pytest tests/genome/test_heteroclinic_cycle.py -q
uv run ruff check tests/genome/test_heteroclinic_cycle.py
uv run ruff format tests/genome/test_heteroclinic_cycle.py
git add tests/genome/test_heteroclinic_cycle.py
git commit -m "genome: #314 per-leg crossing vs W-Z golden (#403)"
```

---

## Final verification (after all tasks)

```bash
uv run pytest tests/genome/test_heteroclinic_cycle.py -q     # all green
uv run pytest tests/genome tests/search -q                   # no regressions in siblings
uv run mypy src/cyclerfinder/genome/heteroclinic_cycle.py
uv run ruff check . && uv run ruff format --check .
```

Then update `data/OUTSTANDING.md` and mark task #314 completed. The cross-system novel search (#405) is unblocked once this lands and the W-Z validation (item 2, Task 8) is green.

## Notes for the implementer

- **The hard part is Task 4** (the connection Newton finding its basin). The connection is a codimension-1 object; a blind start may not converge. If the default `(k_u=1, k_s=1, tau0=0.25*T, branch_u=+1, branch_s=-1)` doesn't close, do a coarse 2-D scan of `(tau_u, tau_s)` over `[0, period)²` computing `||_resid||`, find the minimum, and start Newton there. Add that scan as a private `_scan_starts` helper if needed — it's legitimate and keeps the corrector robust. Record the working configuration.
- **Energy is the master constraint:** every node and every leg lives at one Jacobi C. The guards enforce this; don't relax them.
- **Never fabricate:** if the W-Z connection genuinely won't reproduce at our tolerances, that's a finding to surface (escalate per `feedback_never_give_up_reproducing_papers` — published result, so escalate the corrector, don't accept a silent miss), not a test to weaken.
