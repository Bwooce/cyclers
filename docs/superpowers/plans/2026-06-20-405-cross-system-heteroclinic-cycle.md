# #405 Cross-System Heteroclinic-Cycle Search (Phase A) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a patched-CR3BP Sun-Earth↔Earth-Moon connection matcher, validate it, and run a bounded search for a closed cross-system heteroclinic cycle (a closed cycle OR a documented clean negative).

**Architecture:** New module `src/cyclerfinder/genome/cross_system_cycle.py` treating the SE and EM systems as two `core/cr3bp` instances, reusing #314's `genome/heteroclinic_cycle.py` machinery (`LyapunovNode`, `_seed_on_manifold`, `_section_crossing`, the Newton+coarse-scan pattern). The new layer is (a) an inter-system frame transform via a common Earth-centered inertial frame, (b) a cross-system connection corrector matching inertial position at a patch section, (c) a bounded closure search.

**Tech Stack:** Python 3, numpy, scipy.integrate.solve_ivp (DOP853 + Radau cross-check), pytest, uv, ruff, mypy. Spec: `docs/superpowers/specs/2026-06-20-405-cross-system-heteroclinic-cycle-design.md`.

---

## Conventions every task must follow
- Run from repo root `/home/bruce/dev/cyclers` with `uv run` (uv venv; never pip).
- **Pre-commit mandatory; never `git commit --no-verify`.** Commit with explicit pathspecs (concurrent agents share the tree; never `git add -A`). Do not push.
- **Sourced-golden discipline:** EXPECTED values trace to a publication (Canalias) or are analytical/identity checks; never to our own computed output.
- Type-annotate everything (mypy strict).
- Expensive convergence/search tests get `@pytest.mark.slow` (the repo skips slow by default; run with `-m "slow or not slow"`).

## File structure
- **Create** `src/cyclerfinder/genome/cross_system_cycle.py` — entire Phase-A framework.
- **Create** `tests/genome/test_cross_system_cycle.py` — two-tier validation.
- **(later)** append to `data/negative_results.yaml` if the search is a clean negative.
- Golden (from #407, may not exist yet): `data/golden/canalias_se_em_connection.yaml`. Golden-gated tests use the #314 Task-8 skip pattern (`pytest.mark.skipif(not path.exists())`).

## Reused signatures (verified — do not re-implement)
```python
# core/cr3bp.py
cr3bp.cr3bp_system(primary: str, secondary: str) -> CR3BPSystem   # "Sun","Earth" / "Earth","Moon"
#   CR3BPSystem(mu, primary, secondary, l_km, t_s); t_s = sqrt(l_km^3 / gm_pair) [seconds]
cr3bp.propagate(system, state6, t, *, with_stm=False, rtol=1e-12, atol=1e-12) -> CR3BPArc  # .state_f .stm
cr3bp.cr3bp_eom(t, state6, mu) -> ndarray(6)
cr3bp.jacobi_constant(state6, mu) -> float
# genome/heteroclinic_cycle.py
LyapunovNode.from_libration(system, *, x0_guess, jacobi, period_guess, label, ydot0_sign=1.0, tol=1e-10)
_seed_on_manifold(system, node, *, tau, direction, branch, epsilon, ...) -> ndarray(6)
_section_crossing(system, seed, *, direction, k, max_time, section_y=0.0, ydot_sign=None, method="DOP853") -> ndarray(2)|None
```

---

## Task 1: SE/EM systems + the inter-system frame transform

The crux. Map a 6-state between SE-rotating, a common **Earth-centered inertial** frame (km, km/s), and EM-rotating, parameterized by the relative frame phase `theta` (angle between the Sun-Earth line and the Earth-Moon line). Correctness is gated by a round-trip identity AND a physical-position test.

**Files:** Create `src/cyclerfinder/genome/cross_system_cycle.py`, `tests/genome/test_cross_system_cycle.py`.

- [ ] **Step 1: Write the failing tests**

Create `tests/genome/test_cross_system_cycle.py`:
```python
"""Tests for the #405 cross-system (SE<->EM) heteroclinic-cycle framework (patched CR3BP).

Two-tier validation (per spec): (1) sourced numeric — SE Lyapunov / connection vs
Canalias 2007 C-values (golden, #407); (2) internal-consistency identities for the new
inter-system frame transform (round-trip, ballistic-dV continuity, energy bookkeeping).
"""

from __future__ import annotations

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.cross_system_cycle import FrameBridge


def _bridge() -> FrameBridge:
    se = cr3bp.cr3bp_system("Sun", "Earth")
    em = cr3bp.cr3bp_system("Earth", "Moon")
    return FrameBridge(se=se, em=em)


def test_frame_round_trip_identity() -> None:
    """SE-rot -> inertial -> EM-rot -> inertial -> SE-rot returns the state to <1e-9."""
    bridge = _bridge()
    rng = np.random.default_rng(0)  # seed: deterministic; NOT a sourced value
    for _ in range(20):
        s0 = rng.normal(scale=0.3, size=6)
        s0[0] += 1.0  # place near the Earth region in SE-rot
        theta = float(rng.uniform(0.0, 2.0 * np.pi))
        inert = bridge.se_rot_to_inertial(s0, theta=theta)
        em = bridge.inertial_to_em_rot(inert, theta=theta)
        inert2 = bridge.em_rot_to_inertial(em, theta=theta)
        s1 = bridge.inertial_to_se_rot(inert2, theta=theta)
        assert np.allclose(s0, s1, atol=1e-9), f"round-trip drift {np.abs(s0 - s1).max():.2e}"


def test_moon_maps_to_em_secondary_position() -> None:
    """The EM secondary (Moon) at rest in EM-rot maps to a consistent inertial point.

    Physical anchor (not just an inverse-consistency check): the Moon sits at
    (1-mu_em, 0,0,0,0,0) in EM-rot. Transforming EM-rot -> inertial -> EM-rot must
    return it, and its inertial distance from Earth must equal the EM length scale
    (Earth-Moon distance) to ~1e-6 relative.
    """
    bridge = _bridge()
    moon_em = np.array([1.0 - bridge.em.mu, 0.0, 0.0, 0.0, 0.0, 0.0])
    theta = 0.7
    inert = bridge.em_rot_to_inertial(moon_em, theta=theta)
    # Earth is the EM primary at (-mu_em,0,0); inertial frame is Earth-centered, so the
    # Moon's inertial position magnitude is the Earth-Moon distance = em.l_km.
    r_km = float(np.linalg.norm(inert[:3]))
    assert abs(r_km - bridge.em.l_km) / bridge.em.l_km < 1e-6, f"{r_km} vs {bridge.em.l_km}"
    back = bridge.inertial_to_em_rot(inert, theta=theta)
    assert np.allclose(moon_em, back, atol=1e-9)
```

- [ ] **Step 2: Run** `uv run pytest tests/genome/test_cross_system_cycle.py -x -q` → expect ImportError (`FrameBridge`).

- [ ] **Step 3: Implement the frame bridge.** Create `src/cyclerfinder/genome/cross_system_cycle.py`.

Formulation (Earth-centered inertial in km, km/s; z planar-through but kept 3-D):
- SE-rot state is nondim in SE units (length `se.l_km`, time `se.t_s`, so velocity unit `se.l_km/se.t_s`). Earth is the SE **secondary** at `(1-mu_se, 0, 0)`.
- EM-rot state is nondim in EM units. Earth is the EM **primary** at `(-mu_em, 0, 0)`.
- Common frame: **Earth-centered inertial**, km & km/s. `omega_se = 1/se.t_s`, `omega_em = 1/em.t_s` (rad/s; the CR3BP mean motion is 1 per nondim time).
- `theta` is the inertial angle of the EM-rot x-axis relative to the SE-rot x-axis at the patch instant (the free relative-phase parameter). We place the SE-rot x-axis at inertial angle 0 and the EM-rot x-axis at inertial angle `theta` (a snapshot; the corrector solves over `theta`).

```python
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp


def _rot_z(angle: float) -> NDArray[np.float64]:
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], dtype=np.float64)


@dataclass(frozen=True)
class FrameBridge:
    """Transforms 6-states between SE-rot, Earth-centered inertial (km, km/s), EM-rot.

    SE-rot and EM-rot are nondimensional CR3BP rotating frames; the inertial frame is
    Earth-centered, dimensional. ``theta`` is the inertial angle of the EM-rot x-axis
    minus that of the SE-rot x-axis (the SE-rot x-axis is taken at inertial angle 0).
    """

    se: cr3bp.CR3BPSystem
    em: cr3bp.CR3BPSystem

    # --- SE-rot <-> Earth-centered inertial -------------------------------------
    def se_rot_to_inertial(self, s: NDArray[np.float64], *, theta: float) -> NDArray[np.float64]:
        r = np.asarray(s[:3], float).copy()
        v = np.asarray(s[3:], float).copy()
        r = r - np.array([1.0 - self.se.mu, 0.0, 0.0])  # Earth-centered, SE-rot, nondim
        lkm, vunit = self.se.l_km, self.se.l_km / self.se.t_s
        R = r * lkm            # km, still in SE-rot axes
        Vrot = v * vunit       # km/s, rotating-frame velocity
        omega = np.array([0.0, 0.0, 1.0 / self.se.t_s])
        Vin_serot = Vrot + np.cross(omega, R)          # inertial velocity in SE-rot axes
        # SE-rot x-axis is at inertial angle 0 -> no rotation into inertial axes.
        return np.concatenate([R, Vin_serot])

    def inertial_to_se_rot(self, x: NDArray[np.float64], *, theta: float) -> NDArray[np.float64]:
        R = np.asarray(x[:3], float).copy()
        Vin = np.asarray(x[3:], float).copy()
        omega = np.array([0.0, 0.0, 1.0 / self.se.t_s])
        Vrot = Vin - np.cross(omega, R)
        lkm, vunit = self.se.l_km, self.se.l_km / self.se.t_s
        r = R / lkm + np.array([1.0 - self.se.mu, 0.0, 0.0])
        v = Vrot / vunit
        return np.concatenate([r, v])

    # --- EM-rot <-> Earth-centered inertial -------------------------------------
    def em_rot_to_inertial(self, s: NDArray[np.float64], *, theta: float) -> NDArray[np.float64]:
        r = np.asarray(s[:3], float).copy()
        v = np.asarray(s[3:], float).copy()
        r = r - np.array([-self.em.mu, 0.0, 0.0])      # Earth-centered, EM-rot, nondim
        lkm, vunit = self.em.l_km, self.em.l_km / self.em.t_s
        R_emrot = r * lkm
        Vrot = v * vunit
        omega = np.array([0.0, 0.0, 1.0 / self.em.t_s])
        Vin_emrot = Vrot + np.cross(omega, R_emrot)    # inertial velocity in EM-rot axes
        # EM-rot axes are rotated by +theta relative to inertial -> rotate into inertial.
        Q = _rot_z(theta)
        return np.concatenate([Q @ R_emrot, Q @ Vin_emrot])

    def inertial_to_em_rot(self, x: NDArray[np.float64], *, theta: float) -> NDArray[np.float64]:
        Qi = _rot_z(-theta)
        R_emrot = Qi @ np.asarray(x[:3], float)
        Vin_emrot = Qi @ np.asarray(x[3:], float)
        omega = np.array([0.0, 0.0, 1.0 / self.em.t_s])
        Vrot = Vin_emrot - np.cross(omega, R_emrot)
        lkm, vunit = self.em.l_km, self.em.l_km / self.em.t_s
        r = R_emrot / lkm + np.array([-self.em.mu, 0.0, 0.0])
        v = Vrot / vunit
        return np.concatenate([r, v])
```

- [ ] **Step 4: Run** `uv run pytest tests/genome/test_cross_system_cycle.py -x -q` → expect PASS (2 tests).

  If `test_moon_maps_to_em_secondary_position` fails on the distance check, the centering sign is wrong — re-derive the Earth-offset (EM primary at `-mu_em`, SE secondary at `1-mu_se`). The round-trip test passing but the physical test failing is the signal that forward-transform signs (not just inverse-consistency) need fixing. Do NOT weaken either tolerance.

- [ ] **Step 5: Lint + commit**
```bash
uv run ruff check src/cyclerfinder/genome/cross_system_cycle.py tests/genome/test_cross_system_cycle.py
uv run ruff format src/cyclerfinder/genome/cross_system_cycle.py tests/genome/test_cross_system_cycle.py
uv run mypy src/cyclerfinder/genome/cross_system_cycle.py
git add src/cyclerfinder/genome/cross_system_cycle.py tests/genome/test_cross_system_cycle.py
git commit -m "genome: #405 SE<->EM frame bridge (round-trip + physical-position tests)"
```

---

## Task 2: SE Lyapunov node + Canalias C-value check (sourced tier, golden-gated)

Reuse #314's `LyapunovNode` to build an SE-L1/L2 Lyapunov orbit and confirm it reproduces the Canalias family Jacobi value — pinning the SE μ/Jacobi/section convention. Golden-gated (skip if #407 hasn't landed yet).

**Files:** Modify `tests/genome/test_cross_system_cycle.py`.

- [ ] **Step 1: Write the test**

Append:
```python
import pathlib  # noqa: E402

import pytest  # noqa: E402

from cyclerfinder.genome.heteroclinic_cycle import LyapunovNode  # noqa: E402

_GOLDEN = pathlib.Path("data/golden/canalias_se_em_connection.yaml")

# Canalias-Gómez-Marcote-Masdemont SE L1->L2 heteroclinic-family bifurcation Jacobi,
# confirmed by the Phase-0 scout (C in the SE rotating CR3BP). Inline fallback until the
# #407 golden lands; the golden-backed assertion (below) supersedes it when present.
CANALIAS_C_SE = 3.000863625


@pytest.mark.slow
def test_se_lyapunov_reproduces_canalias_c() -> None:
    """An SE Lyapunov orbit corrected at the Canalias family C closes; jacobi matches.

    EXPECTED = Canalias C=3.000863625 (sourced). Confirms our SE mu/Jacobi/section
    conventions match the Barcelona-school value (modulo a recorded convention offset).
    """
    se = cr3bp.cr3bp_system("Sun", "Earth")
    c_target = CANALIAS_C_SE
    if _GOLDEN.exists():
        import yaml  # noqa: E402  # type: ignore[import-untyped]
        data = yaml.safe_load(_GOLDEN.read_text())
        c_target = float(data["se"]["lyapunov_family_C"]) - float(data["se"].get("c_offset", 0.0))
    # L1 Lyapunov seed: x0 slightly sunward of SE-L1 (~0.9900 in SE-rot). period ~3.06.
    node = LyapunovNode.from_libration(se, x0_guess=0.9899, jacobi=c_target,
                                       period_guess=3.06, label="SE-L1")
    assert node.converged
    assert abs(node.jacobi - c_target) < 1e-6, f"SE jacobi {node.jacobi} vs {c_target}"
```

- [ ] **Step 2: Run** `uv run pytest tests/genome/test_cross_system_cycle.py::test_se_lyapunov_reproduces_canalias_c -q -m "slow or not slow"`.

  Expected PASS. If the L1 seed `x0_guess=0.9899` does not converge to a Lyapunov orbit at this C, scan a small bracket around the SE-L1 libration point (the SE-L1 x is ≈ 1 − (mu_se/3)^(1/3) ≈ 0.99003); document the working seed. If the jacobi misses by a fixed offset, that's the Canalias convention difference — record `c_offset` (do not loosen the 1e-6 check); when #407 lands it carries the offset.

- [ ] **Step 3: Commit**
```bash
uv run ruff check tests/genome/test_cross_system_cycle.py && uv run ruff format tests/genome/test_cross_system_cycle.py
git add tests/genome/test_cross_system_cycle.py
git commit -m "genome: #405 SE Lyapunov vs Canalias C-value (sourced tier, golden-gated)"
```

---

## Task 3: Cross-system connection corrector + ballistic-ΔV continuity

`CrossConnection` + `correct_cross_connection`: propagate the `orbit_from` (EM) unstable manifold to a patch section in the common inertial frame, transform, match the `orbit_to` (SE) stable manifold there. Residual = inertial **position** gap; inertial **velocity** gap = patch ΔV. Free vars `(tau_u, tau_s, theta)`, FD-Jacobian Newton + coarse scan.

**Files:** Modify `src/cyclerfinder/genome/cross_system_cycle.py`, `tests/genome/test_cross_system_cycle.py`.

- [ ] **Step 1: Write the test** (internal-consistency tier — exercises the matcher end-to-end and asserts the patch ΔV is physically low; reproduces a KLMR-class one-way connection)

Append:
```python
from cyclerfinder.genome.cross_system_cycle import (  # noqa: E402
    CrossConnection,
    correct_cross_connection,
)


@pytest.mark.slow
def test_em_to_se_connection_is_low_energy() -> None:
    """An EM-L2 unstable -> SE-L2 stable connection closes its inertial position gap,
    with a low patch ΔV (KLMR 'Shoot the Moon' is a low-energy transfer, ~tens of m/s..
    a few hundred m/s). We assert the position gap closes AND ΔV is bounded (< 1 km/s),
    the physical content of a near-ballistic cross-system connection."""
    se = cr3bp.cr3bp_system("Sun", "Earth")
    em = cr3bp.cr3bp_system("Earth", "Moon")
    bridge = FrameBridge(se=se, em=em)
    # EM-L2 and SE-L2 Lyapunov nodes at their own (small) amplitudes.
    em_l2 = LyapunovNode.from_libration(em, x0_guess=1.155, jacobi=3.15,
                                        period_guess=3.4, label="EM-L2")
    se_l2 = LyapunovNode.from_libration(se, x0_guess=1.0101, jacobi=CANALIAS_C_SE,
                                        period_guess=3.06, label="SE-L2", ydot0_sign=-1.0)
    conn = correct_cross_connection(bridge, em_l2, se_l2, label_from="EM-L2", label_to="SE-L2")
    assert isinstance(conn, CrossConnection)
    assert conn.converged, f"pos residual {conn.residual:.3e} km, n_iter {conn.n_iter}"
    assert conn.residual < 1e2  # inertial position match to <100 km on the patch section
    assert conn.patch_dv_kms < 1.0  # near-ballistic (low-energy transfer)
```

- [ ] **Step 2: Run** → expect ImportError (`correct_cross_connection`).

- [ ] **Step 3: Implement.** Append to the module the `CrossConnection` dataclass and the corrector. The patch section is taken in the **common inertial frame** as a plane (default: the inertial x=0 plane through Earth, i.e. the Sun-Earth perpendicular near the EM region — the implementer picks the section that both manifolds reach; a robust default is "first crossing of the inertial plane normal to the Earth→Sun direction at the EM-L region"). Use the #314 `_seed_on_manifold` to seed each manifold in its OWN system, propagate with `cr3bp.propagate` capturing the trajectory, transform sampled states to inertial via the bridge, and detect the section crossing in inertial coords. Residual = inertial (x,y,z) gap; ΔV = inertial velocity-gap magnitude (km/s). Free vars `(tau_u, tau_s, theta)`; reuse the #314 FD-Jacobian Newton + coarse-scan structure (now 3-D). NEVER raise for "no connection" → `converged=False` + note.

```python
@dataclass(frozen=True)
class CrossConnection:
    label_from: str
    label_to: str
    c_em: float
    c_se: float
    theta: float
    tau_u: float
    tau_s: float
    patch_state_inertial: NDArray[np.float64]  # (6,) km, km/s at the matched patch point
    patch_dv_kms: float
    residual: float           # inertial position gap, km
    converged: bool
    n_iter: int
    notes: str = ""
```

Implementation guidance (the implementer writes the body following this contract + the #314
`correct_connection` pattern):
- `_manifold_inertial_section(bridge, system, node, *, direction, branch, tau, theta, ...)` →
  inertial 6-state at the manifold's first crossing of the inertial patch plane, or `None`.
  (Propagate in-system with `cr3bp.propagate` over a bounded horizon, sample densely, transform
  each sample through the bridge, find the inertial-plane sign change, refine.)
- residual2 (3-vector position gap) and the velocity gap from the two inertial section states.
- Newton over `(tau_u, tau_s, theta)` (3 unknowns, 3 position residuals — square system),
  FD-Jacobian, backtracking, coarse scan over `(tau_u, tau_s, theta)` for the start. Report ΔV.

- [ ] **Step 4: Run** `uv run pytest tests/genome/test_cross_system_cycle.py -q -m "slow or not slow"`.

  CONVERGENCE IS HARD (this is the #405 analog of #314's Task 4). If it does not converge:
  coarse-scan `theta` over `[0, 2π)` first (the relative phase is the most sensitive unknown),
  then `(tau_u, tau_s)`; try EM-L1↔SE-L1 as well as L2↔L2; widen the manifold horizon. Sanity:
  the patch should sit near the EM region ~1.5 Mkm sunward of Earth (the #316 OOM probe). If
  after a genuine scan it will not close to a low ΔV, report DONE_WITH_CONCERNS with the best
  ΔV + position gap + configs tried — a real non-closure is a finding (and feeds the closure
  search's clean-negative path), NOT something to force by loosening tolerances.

- [ ] **Step 5: Lint + commit**
```bash
uv run ruff check src/cyclerfinder/genome/cross_system_cycle.py tests/genome/test_cross_system_cycle.py
uv run ruff format src/cyclerfinder/genome/cross_system_cycle.py tests/genome/test_cross_system_cycle.py
uv run mypy src/cyclerfinder/genome/cross_system_cycle.py
git add src/cyclerfinder/genome/cross_system_cycle.py tests/genome/test_cross_system_cycle.py
git commit -m "genome: #405 cross-system connection corrector (inertial-patch match + ΔV)"
```

---

## Task 4: Bounded closure search + clean-negative registry

`CrossCycle` + `search_cross_cycle`: chain an EM→SE connection and a SE→EM return into a loop that returns to the start orbit (periodic-up-to-rotation) with `theta` advancing commensurately, over a bounded energy×resonance grid. A non-closing grid → registered clean negative.

**Files:** Modify `src/cyclerfinder/genome/cross_system_cycle.py`, `tests/genome/test_cross_system_cycle.py`.

- [ ] **Step 1: Write the test**

Append:
```python
from cyclerfinder.genome.cross_system_cycle import (  # noqa: E402
    CrossCycle,
    search_cross_cycle,
)


@pytest.mark.slow
def test_closure_search_returns_results_or_clean_negative() -> None:
    """The bounded search returns a list of CrossCycle results; each is either closed
    (with bounded independent residual) or honestly open (closed=False, notes set).
    Whatever the outcome, the search NEVER fabricates closure."""
    se = cr3bp.cr3bp_system("Sun", "Earth")
    em = cr3bp.cr3bp_system("Earth", "Moon")
    bridge = FrameBridge(se=se, em=em)
    results = search_cross_cycle(
        bridge,
        c_em_grid=(3.15,),
        c_se_grid=(CANALIAS_C_SE,),
        libration_pairs=(("EM-L2", "SE-L2"),),
        max_attempts=2,
    )
    assert isinstance(results, list)
    for cyc in results:
        assert isinstance(cyc, CrossCycle)
        if cyc.closed:
            assert cyc.max_leg_residual < 1e2
            assert cyc.theta_closure_residual < 1e-2
        else:
            assert cyc.notes  # an honest diagnostic
```

- [ ] **Step 2: Run** → expect ImportError (`search_cross_cycle`).

- [ ] **Step 3: Implement** `CrossCycle` + `search_cross_cycle`. For each grid point: build the EM and SE nodes, `correct_cross_connection` for EM→SE and SE→EM, and check the chain returns to the start orbit with `theta` advancing commensurately (`theta_closure_residual = |theta_return − theta_start mod 2π|`). `closed` iff both legs converged AND `theta_closure_residual < tol`. Return all results. The caller (a run-script, not this task) registers a clean negative if nothing closes.

```python
@dataclass(frozen=True)
class CrossCycle:
    connections: list[CrossConnection]
    c_em: float
    c_se: float
    libration_pair: tuple[str, str]
    theta_closure_residual: float
    closed: bool
    max_leg_residual: float
    independent_residual: float
    notes: str = ""
```

- [ ] **Step 4: Run** `uv run pytest tests/genome/test_cross_system_cycle.py -q -m "slow or not slow"` → expect PASS (results list well-formed; closed or honest-open).

- [ ] **Step 5: Lint + commit**
```bash
uv run ruff check src/cyclerfinder/genome/cross_system_cycle.py tests/genome/test_cross_system_cycle.py
uv run ruff format src/cyclerfinder/genome/cross_system_cycle.py tests/genome/test_cross_system_cycle.py
uv run mypy src/cyclerfinder/genome/cross_system_cycle.py
git add src/cyclerfinder/genome/cross_system_cycle.py tests/genome/test_cross_system_cycle.py
git commit -m "genome: #405 bounded closure search + CrossCycle (periodic-up-to-rotation)"
```

---

## Task 5: Independent Radau cross-check + run-script with clean-negative registry

Reuse #314's discipline: re-derive any found cycle's legs with Radau and fill `independent_residual`. Add a run-script that executes the bounded search and, if nothing closes, registers a method-versioned clean negative in `data/negative_results.yaml`.

**Files:** Modify `src/cyclerfinder/genome/cross_system_cycle.py`, `tests/genome/test_cross_system_cycle.py`; Create `scripts/run_405_cross_system_search.py`.

- [ ] **Step 1: Write the cross-check test**

Append:
```python
from cyclerfinder.genome.cross_system_cycle import crosscheck_cross_cycle  # noqa: E402


@pytest.mark.slow
def test_crosscheck_is_recorded_or_inf() -> None:
    """crosscheck_cross_cycle returns a cycle whose independent_residual is set
    (finite if every converged leg re-derives under Radau, inf if any fails)."""
    se = cr3bp.cr3bp_system("Sun", "Earth")
    em = cr3bp.cr3bp_system("Earth", "Moon")
    bridge = FrameBridge(se=se, em=em)
    results = search_cross_cycle(bridge, c_em_grid=(3.15,), c_se_grid=(CANALIAS_C_SE,),
                                 libration_pairs=(("EM-L2", "SE-L2"),), max_attempts=2)
    for cyc in results:
        checked = crosscheck_cross_cycle(bridge, cyc)
        assert not np.isnan(checked.independent_residual)
```

- [ ] **Step 2: Run** → expect ImportError (`crosscheck_cross_cycle`).

- [ ] **Step 3: Implement** `crosscheck_cross_cycle` (re-derive each converged leg's inertial patch state with `method="Radau"`, compare to the stored `patch_state_inertial`, take the max disagreement; non-reproducing leg → `inf`). Then create `scripts/run_405_cross_system_search.py` that runs `search_cross_cycle` over the planned bounded grid, prints per-grid-point outcomes with timestamps, and — if NOTHING closes — appends a clean-negative entry to `data/negative_results.yaml`:
```python
# scripts/run_405_cross_system_search.py — outline the engineer fills in:
#   - build SE/EM systems + FrameBridge
#   - run search_cross_cycle over the bounded grid (a few C_em, C_se, L1/L2 pairs,
#     the 235:19 Metonic commensurability first)
#   - for any closed cycle: crosscheck_cross_cycle, print residuals
#   - if zero closed: append to data/negative_results.yaml an entry:
#       id: cross_system_se_em_L2_patched_cr3bp
#       issue: 405
#       method: "patched-CR3BP SE<->EM connection matcher (#405 Phase A)"
#       regime: "Sun-Earth + Earth-Moon coupled CR3BP"
#       failed_rung: "closure-search"
#       physical_reason: "<observed: e.g. position gap floors at X km / ΔV floors at Y km/s>"
#       resweep_condition: "BCR4BP Phase B refinement OR wider energy/resonance grid"
```

- [ ] **Step 4: Run** the test + the script:
```bash
uv run pytest tests/genome/test_cross_system_cycle.py -q -m "slow or not slow"
uv run python scripts/run_405_cross_system_search.py
```
Expected: tests pass; the script prints outcomes and (if no closure) writes the registry entry.

- [ ] **Step 5: Lint + commit** (commit code + script; commit the registry entry separately only if the search ran to a genuine clean negative)
```bash
uv run ruff check src/cyclerfinder/genome/cross_system_cycle.py tests/genome/test_cross_system_cycle.py scripts/run_405_cross_system_search.py
uv run ruff format src/cyclerfinder/genome/cross_system_cycle.py tests/genome/test_cross_system_cycle.py scripts/run_405_cross_system_search.py
uv run mypy src/cyclerfinder/genome/cross_system_cycle.py
git add src/cyclerfinder/genome/cross_system_cycle.py tests/genome/test_cross_system_cycle.py scripts/run_405_cross_system_search.py
git commit -m "genome: #405 Radau cross-check + bounded-search run-script (clean-negative registry)"
```

---

## Final verification (after all tasks)
```bash
uv run pytest tests/genome/test_cross_system_cycle.py -q -m "slow or not slow"   # all green
uv run pytest tests/genome tests/search -q                                       # no sibling regressions
uv run mypy src/cyclerfinder/genome/cross_system_cycle.py
uv run ruff check . && uv run ruff format --check .
```
Then update `data/OUTSTANDING.md` (mark #405 Phase A delivered + the closure result: cycle found OR clean negative) and decide on Phase B (#292 BCR4BP refinement) per the outcome. A FOUND closed cycle is a NOVEL discovery → run `search/literature_check.py` novelty gate + the full V0–V5 gauntlet before ANY catalogue claim (downstream, not this plan).

## Notes for the implementer
- **The two hard tasks are 1 (frame transform) and 3 (cross-system Newton).** Task 1's correctness is gated by the round-trip AND the physical Moon-position test — both must pass; round-trip alone can hide a forward-sign error. Task 3 is the #314-Task-4 analog: expect to coarse-scan `theta` first.
- **`theta` is the most sensitive unknown** — it sets the relative SE/EM geometry. Scan it densely.
- **Never fabricate closure.** A clean negative (no closed cross-system cycle on the bounded grid) is a legitimate, registered Phase-A result — this is a novel search, not a reproduction.
- The frame bridge keeps z (3-D capable) though the orbits are planar; do not assume z=0 in the transforms.
