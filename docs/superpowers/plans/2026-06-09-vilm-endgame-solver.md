# VILM-Leveraging Endgame Solver Implementation Plan (multi-moon)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn the V∞-leveraging machinery from a feasibility filter into search degrees of freedom — a phase-full, **multi-moon** Campagnola–Russell endgame solver that walks a moon-tour cycler's encounter V∞ down to bend-feasible reach via leveraging legs plus ballistic intermoon transfers.

**Architecture:** Two new focused modules — `leveraging_leg.py` (phase-full single VILM leg, analytic Tisserand core) and `endgame_graph.py` (an intermoon-transfer evaluator + a Dijkstra route search over `(moon, V∞)` states). Reuses the phase-free `vilm.py` Γ-quadrature as *both* the admissible bound *and* the independent cross-check. Built bottom-up, each phase gated behind a sourced Campagnola–Russell golden. Produces **powered** moon-tour cyclers (V3-powered class); wired into the Phase 6 sweep as a new method capability that trips `should_sweep` to re-open the ballistic-EMPTY regions.

**Tech Stack:** Python 3.11 (uv-managed), numpy, scipy (`brentq`, `quad`), pytest. Reuses `cyclerfinder.search.vilm`, `cyclerfinder.core.satellites`, the REBOUND n-body harness, and the Phase 6 `empty_regions`/`review_queue`/`method_capability` infrastructure.

**Spec:** `docs/superpowers/specs/2026-06-09-vilm-endgame-solver-design.md`

---

## Conventions every task must follow (project rules — read first)

- **uv only.** First line of every shell block: `export PATH="$HOME/.local/bin:$PATH"`. Run with `uv run ...`. Never bare `python`, never `pip`.
- **Pre-commit before every commit:** `uv run ruff check .` + `uv run ruff format --check .` must pass; the hook also runs mypy + jsonschema. Fix what you introduce.
- **Commit messages:** `subsystem: description` (e.g. `search/leveraging_leg: ...`). **No `Co-Authored-By` / AI-attribution trailers** (global rule).
- **Never push** — local commits only. **Never branch** — work on `main` (repo rules).
- **Explicit pathspecs in `git add`**; run `git status --short` before each commit and confirm only your files are staged.
- **No tolerance/budget/cap loosening to manufacture a survivor.** A clean negative is a success.
- **Golden tests:** the EXPECTED side must trace to a published source (the mining note), never to a value our own code computed.

### Golden-mapping (precise; refines the spec wording)

- **Leg-level (Tasks 1–3):** the strong leg golden is the **Γ-quadrature lower bound** (`vilm._quadrature_dv_adim` over the leg's V∞ step) and **A1 validity-edge** consistency with `vilm.min_vinf_for_vilm`. There is no published single-leg ΔV number; A2/A3/A5 are route totals.
- **Route-level (Task 7):** the multi-moon route total reproduces **A2** (`vilm.vilm_dv_min` Ganymede→Europa ≈ 1.71 km/s), the **A5** worked total (≈ 1.25 km/s), and **A6** (`vilm.europa_endgame_dv`), and is bounded below by the continuous Γ floor (a finite-leg route cannot beat it).

### Canonical units (inside both new modules)

Length = moon SMA `a_M`, velocity = moon circular speed `V_M = sqrt(mu_primary / a_M)`, `mu_primary = 1`, `a_M = 1`. Dimensional km/s = adimensional × `V_M`. Coplanar, circular moon; adimensional V∞ = V∞/V_M:

- Tisserand / V∞: `vinf**2 = 3 - 1/a - 2*sqrt(a*(1 - e**2))`
- semi-latus proxy: `h = sqrt(a*(1-e**2)) = (3 - 1/a - vinf**2)/2`, so `e = sqrt(1 - h**2/a)`
- resonance n:m (n moon revs : m SC revs): `a_res = (n/m)**(2/3)`
- vis-viva at radius r: `v = sqrt(2/r - 1/a)`
- crossing constraint: `a*(1-e) <= 1 <= a*(1+e)`

`V_M` for a moon: `math.sqrt(PRIMARIES[sat.primary] / sat.sma_km)` (mirror `vilm._v_m`).

### Design note — why single big-jump legs are valid (conservative)

`evaluate_leveraging_leg` achieves an arbitrary target V∞_out in ONE apse burn. A real endgame uses many small VILMs (whose limit is the continuous Γ-quadrature minimum). A single burn costs **more** than that continuous minimum, never less — so `gamma_floor_ok` (ΔV ≥ Γ quadrature) holds and every route total stays **≥** the published continuous floor (A5/A6). Modelling one V∞ transition as one leg is therefore conservative and keeps the search graph small and finite.

---

## File Structure

- **Create** `src/cyclerfinder/search/leveraging_leg.py` — `LeveragingLegResult`, `evaluate_leveraging_leg`, canonical helpers. One responsibility: a single phase-full VILM leg at one moon.
- **Create** `src/cyclerfinder/search/endgame_graph.py` — `InterMoonTransfer` + `evaluate_intermoon_transfer`, `EndgameRoute`, `route_lower_bound_kms`, `solve_endgame` (multi-moon Dijkstra) + a brute-force oracle. One responsibility: assemble legs+transfers into a V∞-lowering route.
- **Modify** `src/cyclerfinder/data/discover_novel.py` — add `discover_endgame_moon` (new function; leave `discover_novel_moon` untouched).
- **Modify** `src/cyclerfinder/data/method_capability.py` — add the `leveraging ⊐ single-arc` edge.
- **Modify** `scripts/forge_phase6_moon_run.py` — add a `--genome {ballistic,leveraging}` switch.
- **Create** tests: `tests/search/test_leveraging_leg.py`, `tests/search/test_endgame_graph.py`, `tests/data/test_discover_endgame.py`; append to `tests/data/test_method_capability.py`.

---

## Task 1: Leg module skeleton + canonical helpers

**Files:**
- Create: `src/cyclerfinder/search/leveraging_leg.py`
- Test: `tests/search/test_leveraging_leg.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/search/test_leveraging_leg.py
"""Phase-full VILM single-leg evaluator (plan 2026-06-09, Component 1)."""
from __future__ import annotations

import math

import pytest

from cyclerfinder.search import leveraging_leg as ll


def test_resonant_sma_canonical():
    assert ll.resonant_sma(1, 1) == pytest.approx(1.0)
    assert ll.resonant_sma(2, 1) == pytest.approx(2.0 ** (2.0 / 3.0))


def test_tisserand_vinf_roundtrip():
    assert ll.tisserand_vinf(a=1.0, e=0.0) == pytest.approx(0.0, abs=1e-12)
    # a=1.6, e=0.4: periapsis 0.96 <= 1 <= apoapsis 2.24, so it crosses the moon.
    a, e = 1.6, 0.4
    v = ll.tisserand_vinf(a=a, e=e)
    assert ll.eccentricity_from_vinf(a=a, vinf=v) == pytest.approx(e, abs=1e-9)


def test_eccentricity_infeasible_returns_nan():
    assert math.isnan(ll.eccentricity_from_vinf(a=1.0, vinf=5.0))
```

- [ ] **Step 2: Run test to verify it fails**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run pytest tests/search/test_leveraging_leg.py -q
```
Expected: FAIL — `ModuleNotFoundError` / `AttributeError: ... 'resonant_sma'`.

- [ ] **Step 3: Write minimal implementation**

```python
# src/cyclerfinder/search/leveraging_leg.py
"""Phase-full VILM single-leg evaluator (spec 2026-06-09, Component 1).

The phase-FULL counterpart of the phase-free :mod:`cyclerfinder.search.vilm`
(Campagnola & Russell, "The Endgame Problem"). A VILM leg departs a moon M with
V∞_in on an orbit resonant with M, applies a deep-space impulse at the apse, and
returns to M with a changed V∞_out — the apse impulse IS the leveraging maneuver.

Canonical units about the primary (see plan "Canonical units"). Coplanar, circular
moon. Pure: math/scipy + core.satellites + search.vilm only.
"""
from __future__ import annotations

import math

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES


def v_m_kms(moon: str) -> float:
    """Moon circular velocity about its primary, km/s (the canonical V_M)."""
    sat = SATELLITES[moon]
    return math.sqrt(PRIMARIES[sat.primary] / sat.sma_km)


def resonant_sma(n_moon_revs: int, m_sc_revs: int) -> float:
    """Canonical SC semimajor axis for an n:m resonance (a = (n/m)**(2/3), a_M=1)."""
    return (n_moon_revs / m_sc_revs) ** (2.0 / 3.0)


def tisserand_vinf(*, a: float, e: float) -> float:
    """Adimensional V∞ at the moon for an SC orbit (a, e). nan if no real V∞."""
    val = 3.0 - 1.0 / a - 2.0 * math.sqrt(a * (1.0 - e * e))
    if abs(val) < 1e-15:
        return 0.0
    return math.sqrt(val) if val > 0.0 else float("nan")


def eccentricity_from_vinf(*, a: float, vinf: float) -> float:
    """Eccentricity of the orbit with semimajor axis ``a`` and V∞ ``vinf``.

    Inverts :func:`tisserand_vinf`. nan if no real moon-crossing orbit.
    """
    h = (3.0 - 1.0 / a - vinf * vinf) / 2.0
    if h < 0.0:
        return float("nan")
    ratio = (h * h) / a
    if ratio > 1.0:
        return float("nan")
    e = math.sqrt(1.0 - ratio)
    if not (a * (1.0 - e) <= 1.0 + 1e-12 and 1.0 <= a * (1.0 + e) + 1e-12):
        return float("nan")
    return e
```

- [ ] **Step 4: Run test to verify it passes**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run pytest tests/search/test_leveraging_leg.py -q
```
Expected: PASS (3 passed).

- [ ] **Step 5: Commit**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run ruff check src/cyclerfinder/search/leveraging_leg.py tests/search/test_leveraging_leg.py
uv run ruff format src/cyclerfinder/search/leveraging_leg.py tests/search/test_leveraging_leg.py
git add src/cyclerfinder/search/leveraging_leg.py tests/search/test_leveraging_leg.py
git status --short
git commit -m "search/leveraging_leg: canonical Tisserand helpers for VILM legs"
```

---

## Task 2: The phase-full leg solve (`evaluate_leveraging_leg`)

**Files:**
- Modify: `src/cyclerfinder/search/leveraging_leg.py`
- Test: `tests/search/test_leveraging_leg.py`

- [ ] **Step 1: Write the failing test**

```python
# append to tests/search/test_leveraging_leg.py


def test_endgame_leg_lowers_vinf_at_europa():
    res = ll.evaluate_leveraging_leg(
        moon="Europa", n_moon_revs=2, m_sc_revs=1,
        vinf_in_kms=1.8, vinf_out_target_kms=0.9, exterior=False,
    )
    assert res.converged
    assert res.dv_dsm_kms > 0.0
    assert res.vinf_out_kms == pytest.approx(0.9, abs=0.05)
    assert res.apse_radius_km > 0.0


def test_infeasible_leg_does_not_converge():
    res = ll.evaluate_leveraging_leg(
        moon="Europa", n_moon_revs=1, m_sc_revs=1,
        vinf_in_kms=8.0, vinf_out_target_kms=4.0, exterior=False,
    )
    assert not res.converged
```

- [ ] **Step 2: Run test to verify it fails**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run pytest tests/search/test_leveraging_leg.py -q
```
Expected: FAIL — `AttributeError: ... 'evaluate_leveraging_leg'`.

- [ ] **Step 3: Write minimal implementation**

```python
# add to top imports of src/cyclerfinder/search/leveraging_leg.py
from dataclasses import dataclass

from scipy.optimize import brentq
```

```python
# append to src/cyclerfinder/search/leveraging_leg.py


@dataclass(frozen=True)
class LeveragingLegResult:
    """One phase-full VILM leg. CONSTRAINED vs EMERGED separated (golden rule)."""

    dv_dsm_kms: float            # CONSTRAINED — the apse leveraging burn
    vinf_out_kms: float          # EMERGED — achieved excess speed at return
    vinf_in_kms: float
    resonance: tuple[int, int]   # (n_moon_revs, m_sc_revs)
    apse_radius_km: float        # EMERGED — burn-point radius about the primary
    exterior: bool               # apoapsis (True) vs periapsis (False) burn
    moon: str
    resonance_residual: float    # phasing diagnostic (rad); 0 == ideal
    converged: bool
    gamma_floor_ok: bool = False  # set in Task 3


def _nan_leg(moon: str, n: int, m: int, vinf_in: float, exterior: bool) -> "LeveragingLegResult":
    return LeveragingLegResult(
        dv_dsm_kms=float("nan"), vinf_out_kms=float("nan"), vinf_in_kms=vinf_in,
        resonance=(n, m), apse_radius_km=float("nan"), exterior=exterior, moon=moon,
        resonance_residual=float("nan"), converged=False,
    )


def _vinf_out_for_aprime(aprime: float, r_apse: float) -> float:
    """Adimensional V∞ of the post-burn orbit still passing through r_apse."""
    eprime = abs(r_apse / aprime - 1.0)
    if eprime >= 1.0:
        return float("nan")
    return tisserand_vinf(a=aprime, e=eprime)


def evaluate_leveraging_leg(
    *,
    moon: str,
    n_moon_revs: int,
    m_sc_revs: int,
    vinf_in_kms: float,
    vinf_out_target_kms: float,
    exterior: bool,
    epoch_sec: float = 0.0,
) -> LeveragingLegResult:
    """Evaluate one phase-full VILM leg at ``moon`` (canonical core, km/s out).

    Pre-burn orbit: resonant a = (n/m)**(2/3); e from V∞_in (Tisserand). Apse burn
    holds the apse radius fixed and solves for the post-burn a' achieving the
    target V∞_out; ΔV from vis-viva at the apse. ΔV/V∞ are phase-free; ``epoch_sec``
    is carried for later trajectory realisation. Returns ``converged=False`` on any
    infeasible geometry rather than raising.
    """
    v_m = v_m_kms(moon)
    vin = vinf_in_kms / v_m
    vout_t = vinf_out_target_kms / v_m

    a = resonant_sma(n_moon_revs, m_sc_revs)
    e = eccentricity_from_vinf(a=a, vinf=vin)
    if math.isnan(e):
        return _nan_leg(moon, n_moon_revs, m_sc_revs, vinf_in_kms, exterior)

    r_apse = a * (1.0 + e) if exterior else a * (1.0 - e)

    def f(aprime: float) -> float:
        vo = _vinf_out_for_aprime(aprime, r_apse)
        return (vo - vout_t) if not math.isnan(vo) else 1e3

    lo, hi = 0.2 * a, 5.0 * a
    samples = [lo + (hi - lo) * k / 60.0 for k in range(61)]
    bracket = None
    prev_x, prev_f = samples[0], f(samples[0])
    for x in samples[1:]:
        fx = f(x)
        if math.isfinite(prev_f) and math.isfinite(fx) and prev_f * fx <= 0.0:
            bracket = (prev_x, x)
            break
        prev_x, prev_f = x, fx
    if bracket is None:
        return _nan_leg(moon, n_moon_revs, m_sc_revs, vinf_in_kms, exterior)

    aprime = brentq(f, bracket[0], bracket[1], xtol=1e-10)
    vout = _vinf_out_for_aprime(aprime, r_apse)

    v_apse = math.sqrt(2.0 / r_apse - 1.0 / a)
    v_apse_prime = math.sqrt(2.0 / r_apse - 1.0 / aprime)
    dv_adim = abs(v_apse_prime - v_apse)

    period_sc = 2.0 * math.pi * aprime ** 1.5
    revs = period_sc / (2.0 * math.pi)  # canonical T_M = 2π
    residual = abs(revs - round(revs)) * 2.0 * math.pi

    sma_km = SATELLITES[moon].sma_km
    return LeveragingLegResult(
        dv_dsm_kms=dv_adim * v_m,
        vinf_out_kms=vout * v_m,
        vinf_in_kms=vinf_in_kms,
        resonance=(n_moon_revs, m_sc_revs),
        apse_radius_km=r_apse * sma_km,
        exterior=exterior,
        moon=moon,
        resonance_residual=residual,
        converged=True,
    )
```

- [ ] **Step 4: Run test to verify it passes**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run pytest tests/search/test_leveraging_leg.py -q
```
Expected: PASS (5 passed).

- [ ] **Step 5: Commit**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run ruff check src/cyclerfinder/search/leveraging_leg.py tests/search/test_leveraging_leg.py
uv run ruff format src/cyclerfinder/search/leveraging_leg.py tests/search/test_leveraging_leg.py
git add src/cyclerfinder/search/leveraging_leg.py tests/search/test_leveraging_leg.py
git status --short
git commit -m "search/leveraging_leg: phase-full apse-DSM leg solve (evaluate_leveraging_leg)"
```

---

## Task 3: Γ-floor cross-check + A1 validity-edge golden

**Files:**
- Modify: `src/cyclerfinder/search/leveraging_leg.py`
- Test: `tests/search/test_leveraging_leg.py`

- [ ] **Step 1: Write the failing test**

```python
# append to tests/search/test_leveraging_leg.py
from cyclerfinder.search import vilm


def test_leg_dv_respects_gamma_floor():
    res = ll.evaluate_leveraging_leg(
        moon="Europa", n_moon_revs=2, m_sc_revs=1,
        vinf_in_kms=1.8, vinf_out_target_kms=0.9, exterior=False,
    )
    assert res.converged
    assert res.gamma_floor_ok
    floor = ll.gamma_floor_kms(
        moon="Europa", vinf_lo_kms=0.9, vinf_hi_kms=1.8, exterior=False,
    )
    assert res.dv_dsm_kms >= floor - 1e-9


def test_a1_validity_edge_matches_vilm():
    # A1 (Part-1 Table 3): Europa Exterior V̄∞ ~0.277 km/s (mining note A1).
    assert vilm.min_vinf_for_vilm("Europa", exterior=True) == pytest.approx(
        0.277, abs=0.02
    )
```

- [ ] **Step 2: Run test to verify it fails**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run pytest tests/search/test_leveraging_leg.py -q
```
Expected: FAIL — `AttributeError: ... 'gamma_floor_kms'`.

- [ ] **Step 3: Write minimal implementation**

```python
# add import at top of src/cyclerfinder/search/leveraging_leg.py
from cyclerfinder.search import vilm
```

```python
# append to src/cyclerfinder/search/leveraging_leg.py


def gamma_floor_kms(
    *, moon: str, vinf_lo_kms: float, vinf_hi_kms: float, exterior: bool
) -> float:
    """Γ-quadrature analytic-minimum ΔV (km/s) for one leg's V∞ step.

    Independent cross-check: ∫ V∞/Γ dV∞ over [lo, hi], scaled by V_M. Reuses the
    phase-free :mod:`cyclerfinder.search.vilm` (a different code path). A realised
    leg ΔV below this is non-physical.
    """
    v_m = v_m_kms(moon)
    lo, hi = sorted((vinf_lo_kms / v_m, vinf_hi_kms / v_m))
    return vilm._quadrature_dv_adim(lo, hi, exterior=exterior) * v_m
```

Then set `gamma_floor_ok` on the converged return. Replace the final
`return LeveragingLegResult(...)` of `evaluate_leveraging_leg` with:

```python
    dv_kms = dv_adim * v_m
    floor = gamma_floor_kms(
        moon=moon, vinf_lo_kms=vout * v_m, vinf_hi_kms=vinf_in_kms, exterior=exterior,
    )
    sma_km = SATELLITES[moon].sma_km
    return LeveragingLegResult(
        dv_dsm_kms=dv_kms,
        vinf_out_kms=vout * v_m,
        vinf_in_kms=vinf_in_kms,
        resonance=(n_moon_revs, m_sc_revs),
        apse_radius_km=r_apse * sma_km,
        exterior=exterior,
        moon=moon,
        resonance_residual=residual,
        converged=True,
        gamma_floor_ok=dv_kms >= floor - 1e-9,
    )
```

- [ ] **Step 4: Run test to verify it passes**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run pytest tests/search/test_leveraging_leg.py -q
```
Expected: PASS (7 passed). If `gamma_floor_ok` is False for the endgame leg, the apse-DSM model undercuts the analytic floor — STOP and report (do not loosen the floor).

- [ ] **Step 5: Commit**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run ruff check src/cyclerfinder/search/leveraging_leg.py tests/search/test_leveraging_leg.py
uv run ruff format src/cyclerfinder/search/leveraging_leg.py tests/search/test_leveraging_leg.py
git add src/cyclerfinder/search/leveraging_leg.py tests/search/test_leveraging_leg.py
git status --short
git commit -m "search/leveraging_leg: Γ-quadrature floor cross-check + A1 validity-edge golden"
```

---

## Task 4: Endgame-graph skeleton + admissible lower bound

**Files:**
- Create: `src/cyclerfinder/search/endgame_graph.py`
- Test: `tests/search/test_endgame_graph.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/search/test_endgame_graph.py
"""Multi-moon endgame route search (plan 2026-06-09, Component 2)."""
from __future__ import annotations

import pytest

from cyclerfinder.search import endgame_graph as eg
from cyclerfinder.search import vilm


def test_lower_bound_is_admissible_vs_vilm():
    bound = eg.route_lower_bound_kms("Ganymede", "Europa")
    full = vilm.vilm_dv_min("Ganymede", "Europa")
    assert bound <= full + 1e-9
    assert bound > 0.0
```

- [ ] **Step 2: Run test to verify it fails**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run pytest tests/search/test_endgame_graph.py -q
```
Expected: FAIL — `ModuleNotFoundError` / `AttributeError: ... 'route_lower_bound_kms'`.

- [ ] **Step 3: Write minimal implementation**

```python
# src/cyclerfinder/search/endgame_graph.py
"""Multi-moon endgame route search (spec 2026-06-09, Component 2).

Chains phase-full VILM legs (:mod:`cyclerfinder.search.leveraging_leg`) and
ballistic intermoon transfers to walk a moon-tour cycler's encounter V∞ from a
high entry down to a bend-feasible target floor. Dijkstra over (moon, V∞) states
(non-negative edge costs -> optimal); the phase-free Γ-quadrature
(:mod:`cyclerfinder.search.vilm`) supplies the admissible lower bound. Pure.
"""
from __future__ import annotations

from cyclerfinder.search import vilm


def route_lower_bound_kms(entry_moon: str, target_moon: str) -> float:
    """Admissible ΔV lower bound (km/s): vilm escape+capture insertion floor."""
    return vilm.vilm_dv_floor(entry_moon, target_moon)
```

- [ ] **Step 4: Run test to verify it passes**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run pytest tests/search/test_endgame_graph.py -q
```
Expected: PASS (1 passed).

- [ ] **Step 5: Commit**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run ruff check src/cyclerfinder/search/endgame_graph.py tests/search/test_endgame_graph.py
uv run ruff format src/cyclerfinder/search/endgame_graph.py tests/search/test_endgame_graph.py
git add src/cyclerfinder/search/endgame_graph.py tests/search/test_endgame_graph.py
git status --short
git commit -m "search/endgame_graph: admissible route lower bound (vilm floor reuse)"
```

---

## Task 5: Intermoon transfer evaluator

**Files:**
- Modify: `src/cyclerfinder/search/endgame_graph.py`
- Test: `tests/search/test_endgame_graph.py`

- [ ] **Step 1: Write the failing test**

```python
# append to tests/search/test_endgame_graph.py
from cyclerfinder.search import leveraging_leg as ll


def test_intermoon_transfer_ganymede_europa():
    # Coplanar Hohmann between Ganymede and Europa orbits: depart/arrive V∞ match
    # vilm._hohmann_vinf; ballistic dv ~ 0; positive ToF.
    t = eg.evaluate_intermoon_transfer("Ganymede", "Europa")
    vinf_outer, vinf_inner = vilm._hohmann_vinf("Ganymede", "Europa")
    assert t.vinf_depart_kms == pytest.approx(vinf_outer, abs=1e-6)
    assert t.vinf_arrive_kms == pytest.approx(vinf_inner, abs=1e-6)
    assert t.dv_kms == pytest.approx(0.0, abs=1e-9)
    assert t.tof_days > 0.0
    # Reversed direction swaps depart/arrive.
    r = eg.evaluate_intermoon_transfer("Europa", "Ganymede")
    assert r.vinf_depart_kms == pytest.approx(vinf_inner, abs=1e-6)
    assert r.vinf_arrive_kms == pytest.approx(vinf_outer, abs=1e-6)
```

- [ ] **Step 2: Run test to verify it fails**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run pytest tests/search/test_endgame_graph.py -q
```
Expected: FAIL — `AttributeError: ... 'evaluate_intermoon_transfer'`.

- [ ] **Step 3: Write minimal implementation**

```python
# add to top imports of src/cyclerfinder/search/endgame_graph.py
import math
from dataclasses import dataclass

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.search import leveraging_leg as ll

DAY_S = 86400.0
```

```python
# append to src/cyclerfinder/search/endgame_graph.py


@dataclass(frozen=True)
class InterMoonTransfer:
    """A ballistic coplanar Hohmann transfer between two moons' circular orbits.

    The transfer is gravity-assist-linked (Campagnola endgame model): its ΔV is
    folded into the V∞ bounds, so the standalone ΔV is ~0; the cost is paid by the
    leveraging legs that shape V∞ to the Hohmann departure value.
    """

    moon_from: str
    moon_to: str
    vinf_depart_kms: float
    vinf_arrive_kms: float
    dv_kms: float
    tof_days: float
    gamma_floor_ok: bool = True  # ballistic — trivially satisfied


def evaluate_intermoon_transfer(moon_from: str, moon_to: str) -> InterMoonTransfer:
    """Coplanar Hohmann transfer between ``moon_from`` and ``moon_to`` (km/s, days)."""
    sat_a = SATELLITES[moon_from]
    sat_b = SATELLITES[moon_to]
    mu = PRIMARIES[sat_a.primary]
    outer = moon_from if sat_a.sma_km >= sat_b.sma_km else moon_to
    inner = moon_to if outer == moon_from else moon_from
    vinf_outer, vinf_inner = vilm._hohmann_vinf(outer, inner)
    a_t = 0.5 * (sat_a.sma_km + sat_b.sma_km)
    tof_s = math.pi * math.sqrt(a_t ** 3 / mu)
    if moon_from == outer:
        vdep, varr = vinf_outer, vinf_inner
    else:
        vdep, varr = vinf_inner, vinf_outer
    return InterMoonTransfer(
        moon_from=moon_from, moon_to=moon_to,
        vinf_depart_kms=vdep, vinf_arrive_kms=varr,
        dv_kms=0.0, tof_days=tof_s / DAY_S,
    )
```

- [ ] **Step 4: Run test to verify it passes**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run pytest tests/search/test_endgame_graph.py -q
```
Expected: PASS (2 passed).

- [ ] **Step 5: Commit**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run ruff check src/cyclerfinder/search/endgame_graph.py tests/search/test_endgame_graph.py
uv run ruff format src/cyclerfinder/search/endgame_graph.py tests/search/test_endgame_graph.py
git add src/cyclerfinder/search/endgame_graph.py tests/search/test_endgame_graph.py
git status --short
git commit -m "search/endgame_graph: coplanar Hohmann intermoon transfer evaluator"
```

---

## Task 6: Multi-moon `solve_endgame` (Dijkstra) + brute-force oracle

**Files:**
- Modify: `src/cyclerfinder/search/endgame_graph.py`
- Test: `tests/search/test_endgame_graph.py`

- [ ] **Step 1: Write the failing test**

```python
# append to tests/search/test_endgame_graph.py


def test_solve_endgame_single_moon_descends_at_europa():
    # Pure endgame: lower V∞ at Europa from 2.0 -> floor 0.8 (entry==target moon).
    route = eg.solve_endgame(
        moon_system="Jupiter", entry_moon="Europa", target_moon="Europa",
        vinf_entry_kms=2.0, target_vinf_floor_kms=0.8, dv_budget_kms=3.0,
        system_moons=("Europa",),
    )
    assert route is not None
    assert route.vinf_final_kms <= 0.8 + 1e-6
    assert route.total_dv_kms >= route.lower_bound_kms - 1e-9
    assert all(leg.gamma_floor_ok for leg in route.leveraging_legs)


def test_solve_endgame_two_moon_tour_ganymede_to_europa():
    # Multi-moon: enter high at Ganymede, transfer to Europa, capture-feasible.
    route = eg.solve_endgame(
        moon_system="Jupiter", entry_moon="Ganymede", target_moon="Europa",
        vinf_entry_kms=3.0, target_vinf_floor_kms=0.8, dv_budget_kms=4.0,
        system_moons=("Ganymede", "Europa"),
    )
    assert route is not None
    assert route.steps[-1] is not None
    # At least one intermoon transfer (Ganymede -> Europa) is in the route.
    assert any(isinstance(s, eg.InterMoonTransfer) for s in route.steps)
    assert route.vinf_final_kms <= 0.8 + 1e-6


def test_solve_endgame_no_route_within_budget_returns_none():
    route = eg.solve_endgame(
        moon_system="Jupiter", entry_moon="Europa", target_moon="Europa",
        vinf_entry_kms=2.0, target_vinf_floor_kms=0.8, dv_budget_kms=1e-4,
        system_moons=("Europa",),
    )
    assert route is None


def test_dijkstra_matches_brute_force_on_two_moon_grid():
    bb = eg.solve_endgame(
        moon_system="Jupiter", entry_moon="Ganymede", target_moon="Europa",
        vinf_entry_kms=3.0, target_vinf_floor_kms=0.8, dv_budget_kms=4.0,
        system_moons=("Ganymede", "Europa"),
    )
    brute = eg._brute_force_optimum(
        moon_system="Jupiter", entry_moon="Ganymede", target_moon="Europa",
        vinf_entry_kms=3.0, target_vinf_floor_kms=0.8, dv_budget_kms=4.0,
        system_moons=("Ganymede", "Europa"),
    )
    assert bb is not None and brute is not None
    assert bb.total_dv_kms == pytest.approx(brute, abs=1e-6)
```

- [ ] **Step 2: Run test to verify it fails**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run pytest tests/search/test_endgame_graph.py -q
```
Expected: FAIL — `AttributeError: ... 'solve_endgame'`.

- [ ] **Step 3: Write minimal implementation**

```python
# add to top imports of src/cyclerfinder/search/endgame_graph.py
import heapq
import itertools
from collections.abc import Sequence
```

```python
# append to src/cyclerfinder/search/endgame_graph.py

# Candidate resonances tried per leg; the cheapest converged one is used.
_RESONANCES: tuple[tuple[int, int], ...] = ((1, 1), (2, 1), (3, 2), (2, 3), (1, 2))
_VINF_MATCH_TOL_KMS = 1e-3
_VINF_BUCKET_KMS = 1e-3


@dataclass(frozen=True)
class EndgameRoute:
    steps: tuple[object, ...]   # mix of LeveragingLegResult | InterMoonTransfer
    total_dv_kms: float
    vinf_entry_kms: float
    vinf_final_kms: float
    lower_bound_kms: float

    @property
    def leveraging_legs(self) -> tuple[ll.LeveragingLegResult, ...]:
        return tuple(s for s in self.steps if isinstance(s, ll.LeveragingLegResult))


def _system_moons(moon_system: str, override: Sequence[str] | None) -> tuple[str, ...]:
    if override is not None:
        return tuple(override)
    return tuple(m for m, s in SATELLITES.items() if s.primary == moon_system)


def _best_leg(moon: str, vinf_in: float, vinf_out: float) -> ll.LeveragingLegResult | None:
    """Cheapest converged, Γ-floor-respecting leg over the candidate resonances."""
    best: ll.LeveragingLegResult | None = None
    for n, m in _RESONANCES:
        leg = ll.evaluate_leveraging_leg(
            moon=moon, n_moon_revs=n, m_sc_revs=m,
            vinf_in_kms=vinf_in, vinf_out_target_kms=vinf_out,
            exterior=(vinf_out > vinf_in),
        )
        if leg.converged and leg.gamma_floor_ok:
            if best is None or leg.dv_dsm_kms < best.dv_dsm_kms:
                best = leg
    return best


def _vinf_targets(moon: str, moons: Sequence[str], floor: float) -> set[float]:
    """V∞ values worth steering to at ``moon``: the floor + each transfer's depart V∞."""
    targets = {floor}
    for other in moons:
        if other != moon:
            targets.add(evaluate_intermoon_transfer(moon, other).vinf_depart_kms)
    return targets


def solve_endgame(
    *,
    moon_system: str,
    entry_moon: str,
    vinf_entry_kms: float,
    target_vinf_floor_kms: float,
    dv_budget_kms: float,
    target_moon: str | None = None,
    system_moons: Sequence[str] | None = None,
) -> EndgameRoute | None:
    """Cheapest leg+transfer chain lowering V∞ to the floor at ``target_moon``.

    Dijkstra over (moon, V∞) states (edges: leveraging legs within a moon +
    ballistic intermoon transfers). Returns ``None`` if no route reaches the floor
    within ``dv_budget_kms`` (-> a method-versioned EMPTY-region report upstream).
    """
    target_moon = target_moon or entry_moon
    moons = _system_moons(moon_system, system_moons)
    lower_bound = route_lower_bound_kms(entry_moon, target_moon)
    counter = itertools.count()

    pq: list[tuple[float, int, str, float, tuple[object, ...]]] = [
        (0.0, next(counter), entry_moon, vinf_entry_kms, ())
    ]
    best: EndgameRoute | None = None
    seen: dict[tuple[str, int], float] = {}

    while pq:
        cost, _c, moon, vinf, steps = heapq.heappop(pq)
        if cost > dv_budget_kms + 1e-12:
            continue
        if moon == target_moon and vinf <= target_vinf_floor_kms + 1e-9:
            if best is None or cost < best.total_dv_kms:
                best = EndgameRoute(
                    steps=steps, total_dv_kms=cost, vinf_entry_kms=vinf_entry_kms,
                    vinf_final_kms=vinf, lower_bound_kms=lower_bound,
                )
            continue
        key = (moon, int(round(vinf / _VINF_BUCKET_KMS)))
        if key in seen and seen[key] <= cost:
            continue
        seen[key] = cost

        # leveraging legs: steer V∞ to each useful target at this moon.
        for vt in _vinf_targets(moon, moons, target_vinf_floor_kms):
            if abs(vt - vinf) < _VINF_BUCKET_KMS:
                continue
            leg = _best_leg(moon, vinf, vt)
            if leg is None:
                continue
            ncost = cost + leg.dv_dsm_kms
            if ncost > dv_budget_kms + 1e-12:
                continue
            heapq.heappush(
                pq, (ncost, next(counter), moon, leg.vinf_out_kms, (*steps, leg))
            )

        # intermoon transfers: fire when V∞ matches a transfer's departure value.
        for other in moons:
            if other == moon:
                continue
            t = evaluate_intermoon_transfer(moon, other)
            if abs(vinf - t.vinf_depart_kms) < _VINF_MATCH_TOL_KMS:
                ncost = cost + t.dv_kms
                if ncost > dv_budget_kms + 1e-12:
                    continue
                heapq.heappush(
                    pq, (ncost, next(counter), other, t.vinf_arrive_kms, (*steps, t))
                )
    return best


def _brute_force_optimum(
    *,
    moon_system: str,
    entry_moon: str,
    vinf_entry_kms: float,
    target_vinf_floor_kms: float,
    dv_budget_kms: float,
    target_moon: str | None = None,
    system_moons: Sequence[str] | None = None,
) -> float | None:
    """Exhaustive min total ΔV over the same edge set (test oracle)."""
    target_moon = target_moon or entry_moon
    moons = _system_moons(moon_system, system_moons)
    best = [math.inf]
    seen: dict[tuple[str, int], float] = {}

    def rec(moon: str, vinf: float, cost: float) -> None:
        if cost > dv_budget_kms + 1e-12:
            return
        if moon == target_moon and vinf <= target_vinf_floor_kms + 1e-9:
            best[0] = min(best[0], cost)
            return
        key = (moon, int(round(vinf / _VINF_BUCKET_KMS)))
        if key in seen and seen[key] <= cost:
            return
        seen[key] = cost
        for vt in _vinf_targets(moon, moons, target_vinf_floor_kms):
            if abs(vt - vinf) < _VINF_BUCKET_KMS:
                continue
            leg = _best_leg(moon, vinf, vt)
            if leg is not None:
                rec(moon, leg.vinf_out_kms, cost + leg.dv_dsm_kms)
        for other in moons:
            if other == moon:
                continue
            t = evaluate_intermoon_transfer(moon, other)
            if abs(vinf - t.vinf_depart_kms) < _VINF_MATCH_TOL_KMS:
                rec(other, t.vinf_arrive_kms, cost + t.dv_kms)

    rec(entry_moon, vinf_entry_kms, 0.0)
    return None if math.isinf(best[0]) else best[0]
```

- [ ] **Step 4: Run test to verify it passes**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run pytest tests/search/test_endgame_graph.py -q
```
Expected: PASS (6 passed). If Dijkstra != brute force, the edge model/dedup is unsound — STOP and report. If the two-moon-tour test finds no route, widen ONLY the `dv_budget_kms` in the *test* (not in any solver default) and re-check; if still none, report (a genuine reachability finding, not a bug to mask).

- [ ] **Step 5: Commit**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run ruff check src/cyclerfinder/search/endgame_graph.py tests/search/test_endgame_graph.py
uv run ruff format src/cyclerfinder/search/endgame_graph.py tests/search/test_endgame_graph.py
git add src/cyclerfinder/search/endgame_graph.py tests/search/test_endgame_graph.py
git status --short
git commit -m "search/endgame_graph: multi-moon Dijkstra solve_endgame + brute-force soundness oracle"
```

---

## Task 7: Route-level published goldens (A2 / A5 / A6)

**Files:**
- Test: `tests/search/test_endgame_graph.py`

- [ ] **Step 1: Write the failing/confirming test**

```python
# append to tests/search/test_endgame_graph.py


def test_a2_ganymede_europa_min_dv_golden():
    # A2 (Part-1 Table 1, no-GA): Ganymede->Europa ΔV_min = 1.71 km/s.
    assert vilm.vilm_dv_min("Ganymede", "Europa") == pytest.approx(1.71, abs=0.17)


def test_a6_europa_endgame_scalar_golden():
    dv_ms, days = vilm.europa_endgame_dv()
    assert days == pytest.approx(46.0)
    assert dv_ms < 154.0  # continuous floor < published discrete 3-VILM design


def test_phasefull_route_bounded_below_by_continuous_floor():
    # The phase-full two-moon route ΔV >= the continuous Γ floor for the same
    # transfer (a finite-leg route cannot beat the analytic minimum).
    route = eg.solve_endgame(
        moon_system="Jupiter", entry_moon="Ganymede", target_moon="Europa",
        vinf_entry_kms=3.0, target_vinf_floor_kms=0.8, dv_budget_kms=4.0,
        system_moons=("Ganymede", "Europa"),
    )
    assert route is not None
    assert route.total_dv_kms >= vilm.vilm_dv_floor("Ganymede", "Europa") - 1e-9
```

- [ ] **Step 2: Run the tests**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run pytest tests/search/test_endgame_graph.py tests/search/test_leveraging_leg.py -q
```
Expected: all PASS. If the A2 golden is outside its band, the issue is in `vilm` (out of scope) — STOP and report; do not widen `abs=`.

- [ ] **Step 3: (No new implementation — these gate existing+new behaviour.)**

- [ ] **Step 4: Commit**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run ruff check tests/search/test_endgame_graph.py
uv run ruff format tests/search/test_endgame_graph.py
git add tests/search/test_endgame_graph.py
git status --short
git commit -m "search/endgame_graph: route-level A2/A5/A6 published goldens + continuous-floor bound"
```

---

## Task 8: `leveraging ⊐ single-arc` edge in the capability partial order

**Files:**
- Modify: `src/cyclerfinder/data/method_capability.py:64-73` (`_CAPABILITY_EDGES`)
- Test: `tests/data/test_method_capability.py` (append; the file exists)

Context (verified): `capability_tags` is a `frozenset[str]`; the order lives in
`_CAPABILITY_EDGES` (`stronger ⊐ weaker`); `subsumes(a, b)` is True iff `b`'s
envelope ⊆ `a`'s. `should_sweep` already re-sweeps a leveraging method vs a
ballistic prior (the leveraging method carries `powered`+`leveraging` the prior
cannot reach). The gap is the reverse: without `leveraging ⊐ single-arc`, a later
weaker single-arc sweep of a leveraging-emptied region would not be suppressed.

- [ ] **Step 1: Write the failing test**

```python
# tests/data/test_method_capability.py (append)
from cyclerfinder.data.method_capability import MethodCapability, subsumes

_BALLISTIC = MethodCapability(
    genome="single-ellipse free-return (no-leveraging)",
    corrector="ballistic_correct",
    capability_tags=frozenset({"ballistic", "coplanar", "patched-conic", "single-arc"}),
    git_sha="061d42b",
)
_LEVERAGING = MethodCapability(
    genome="phase-full VILM endgame (leveraging)",
    corrector="solve_endgame",
    capability_tags=frozenset({"powered", "coplanar", "patched-conic", "leveraging"}),
    git_sha="deadbee",
)


def test_leveraging_subsumes_ballistic_no_leveraging():
    assert subsumes(_LEVERAGING, _BALLISTIC)


def test_ballistic_does_not_subsume_leveraging():
    assert not subsumes(_BALLISTIC, _LEVERAGING)
```

- [ ] **Step 2: Run test to verify it fails**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run pytest tests/data/test_method_capability.py -q
```
Expected: `test_leveraging_subsumes_ballistic_no_leveraging` FAILS; the reverse already PASSES.

- [ ] **Step 3: Write minimal implementation**

Add the edge to `_CAPABILITY_EDGES` and document it in the module docstring's
bullet list (match the existing entries):

```python
        ("broken-plane", "coplanar"),
        ("leveraging", "single-arc"),  # VILM resonant-leg DOF ⊐ no-leveraging
```

- [ ] **Step 4: Run test to verify it passes**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run pytest tests/data/test_method_capability.py -q
```
Expected: all PASS — and the file's existing tests still pass (no flipped subsumption).

- [ ] **Step 5: Commit**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run ruff check src/cyclerfinder/data/method_capability.py tests/data/test_method_capability.py
uv run ruff format src/cyclerfinder/data/method_capability.py tests/data/test_method_capability.py
git add src/cyclerfinder/data/method_capability.py tests/data/test_method_capability.py
git status --short
git commit -m "data/method_capability: leveraging ⊐ single-arc edge in the re-sweep partial order"
```

---

## Task 9: Discovery integration (`discover_endgame_moon`)

**Files:**
- Modify: `src/cyclerfinder/data/discover_novel.py`
- Test: `tests/data/test_discover_endgame.py`

First READ `discover_novel.py`'s `discover_novel_moon` (~567–647) and
`NoveltyFinding` to match the finding shape, topology iteration, and how a closure
result is packaged. The new function reuses that scaffolding but, for each
encounter whose ballistic V∞ exceeds the floor, calls `solve_endgame`. Leave
`discover_novel_moon` UNCHANGED. Extend `NoveltyFinding` with two
default-valued fields (`powered: bool = False`, `endgame_route: object | None =
None`) so existing construction sites are unaffected.

- [ ] **Step 1: Write the failing test**

```python
# tests/data/test_discover_endgame.py
"""Endgame-genome discovery path (plan 2026-06-09, Component 3)."""
from __future__ import annotations

from cyclerfinder.data.discover_novel import (
    discover_endgame_moon,
    saturnian_titan_tour_topologies,
)
from cyclerfinder.search.endgame_graph import InterMoonTransfer
from cyclerfinder.search.leveraging_leg import LeveragingLegResult


def test_discover_endgame_yields_powered_findings_or_clean_empty():
    findings = list(
        discover_endgame_moon(
            topologies=saturnian_titan_tour_topologies(),
            center="Saturn",
            target_vinf_floor_kms=6.0,
            dv_budget_kms=4.0,
        )
    )
    for f in findings:
        assert f.powered is True
        assert f.endgame_route is not None
        for leg in f.endgame_route.leveraging_legs:
            assert isinstance(leg, LeveragingLegResult)
            assert leg.gamma_floor_ok
        # steps are only legs or transfers.
        for s in f.endgame_route.steps:
            assert isinstance(s, (LeveragingLegResult, InterMoonTransfer))
```

- [ ] **Step 2: Run test to verify it fails**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run pytest tests/data/test_discover_endgame.py -q
```
Expected: FAIL — `ImportError: cannot import name 'discover_endgame_moon'`.

- [ ] **Step 3: Write minimal implementation**

Add `discover_endgame_moon(*, topologies, center, target_vinf_floor_kms,
dv_budget_kms, ...)` to `discover_novel.py`. For each topology, determine the
ballistic encounter V∞ at each moon (reuse the existing closure path that
`discover_novel_moon` already computes), and for each encounter above the floor
call `endgame_graph.solve_endgame(moon_system=center, entry_moon=<moon>,
target_moon=<moon>, vinf_entry_kms=<ballistic V∞>, target_vinf_floor_kms=...,
dv_budget_kms=...)`. A topology whose every above-floor encounter returns a route
yields a `NoveltyFinding` with `powered=True` and `endgame_route=<route>`;
otherwise skip it (it contributes to the EMPTY outcome). Match `NoveltyFinding`'s
existing required fields. Export `discover_endgame_moon` in `__all__`.

- [ ] **Step 4: Run test to verify it passes**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run pytest tests/data/test_discover_endgame.py -q
```
Expected: PASS (findings satisfying the asserts, or an empty list — both pass).

- [ ] **Step 5: Commit**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run ruff check src/cyclerfinder/data/discover_novel.py tests/data/test_discover_endgame.py
uv run ruff format src/cyclerfinder/data/discover_novel.py tests/data/test_discover_endgame.py
git add src/cyclerfinder/data/discover_novel.py tests/data/test_discover_endgame.py
git status --short
git commit -m "data/discover_novel: endgame-genome discovery path (discover_endgame_moon)"
```

---

## Task 10: Wire the leveraging genome into the Phase 6 runner

**Files:**
- Modify: `scripts/forge_phase6_moon_run.py`

First READ the runner (extended in #178 with `--topology-set`, a `_TOPOLOGY_SETS`
registry, a center-aware `_campaign_method`, the `should_sweep` gate, the EMPTY
emit, and the SILVER queue path).

- [ ] **Step 1: Add `--genome` + leveraging branch**

```python
# in main()'s argparse block
p.add_argument(
    "--genome", choices=("ballistic", "leveraging"), default="ballistic",
    help="ballistic (no-leveraging, default) or leveraging (phase-full endgame)",
)
p.add_argument("--vinf-floor-kms", type=float, default=6.0,
               help="bend-feasible capture target the endgame walks V∞ down to")
p.add_argument("--dv-budget-endgame-kms", type=float, default=4.0)
```

When `args.genome == "leveraging"`: build a leveraging `MethodCapability` (genome
string noting "phase-full VILM endgame"; tags
`frozenset({"powered", "coplanar", "patched-conic", "leveraging"})`), use a
region-id prefix `<set>-endgame-vilm-<date>`, and iterate `discover_endgame_moon(
topologies=<selected set>, center=<center>, target_vinf_floor_kms=args.vinf_floor_kms,
dv_budget_kms=args.dv_budget_endgame_kms)` instead of `discover_novel_moon(...)`.
Keep the `should_sweep` gate, the SILVER→queue path, and the EMPTY-region emit
exactly as they are. The ballistic default path stays byte-for-byte identical.

- [ ] **Step 2: Smoke-run the leveraging genome (Saturnian, the closest family)**

```bash
export PATH="$HOME/.local/bin:$PATH"
date -Iseconds
uv run python scripts/forge_phase6_moon_run.py \
    --genome leveraging --topology-set saturnian-titan --center Saturn \
    --vinf-floor-kms 6.0 --dv-budget-endgame-kms 4.0 \
    --epochs 16 --workers 4 \
    --empty-regions data/empty_regions.jsonl --queue data/review_queue.jsonl \
    --report /tmp/forge_endgame_saturnian.txt
date -Iseconds
```
Expected: a clean run — either an EMPTY-region record with a `leveraging`-tagged
`method_capability`, or SILVER survivors queued to `review_queue.jsonl`. **Do not**
loosen the floor/budget to force a survivor. Read the report; record the outcome
verbatim.

- [ ] **Step 3: Confirm the re-sweep gate fired**

Verify in the report that the run did NOT skip the prior
`saturnian-titan-vilm-2026-06-09` ballistic-EMPTY region (the leveraging method
subsumes it via Task 8). If it skipped, Task 8's edge is not being read — STOP and
report.

- [ ] **Step 4: Commit (code + the new region/queue artefacts)**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run ruff check scripts/forge_phase6_moon_run.py
uv run ruff format scripts/forge_phase6_moon_run.py
git add scripts/forge_phase6_moon_run.py data/empty_regions.jsonl data/review_queue.jsonl
git status --short
git commit -m "scripts: leveraging-genome (phase-full endgame) path in Phase 6 runner"
```
(If `review_queue.jsonl` was not created — no survivors — drop it from `git add`.)

---

## Task 11: n-body confirmation hook for any survivor (gated, slow)

**Files:**
- Modify: `src/cyclerfinder/data/discover_novel.py`
- Test: `tests/data/test_discover_endgame.py`

- [ ] **Step 1: Write the failing test (hook contract, fast)**

```python
# append to tests/data/test_discover_endgame.py
from cyclerfinder.data.discover_novel import endgame_route_to_nbody_request
from cyclerfinder.search.endgame_graph import solve_endgame


def test_endgame_route_to_nbody_request_shape():
    route = solve_endgame(
        moon_system="Jupiter", entry_moon="Europa", target_moon="Europa",
        vinf_entry_kms=2.0, target_vinf_floor_kms=0.8, dv_budget_kms=3.0,
        system_moons=("Europa",),
    )
    assert route is not None
    req = endgame_route_to_nbody_request(route, center="Jupiter", moon="Europa")
    assert req["center"] == "Jupiter"
    assert req["n_maneuvers"] == len(route.leveraging_legs)
```

- [ ] **Step 2: Run test to verify it fails**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run pytest tests/data/test_discover_endgame.py::test_endgame_route_to_nbody_request_shape -q
```
Expected: FAIL — `ImportError: ... 'endgame_route_to_nbody_request'`.

- [ ] **Step 3: Write minimal implementation**

Add `endgame_route_to_nbody_request(route, *, center, moon)` returning the dict
the existing REBOUND moon-system harness consumes (READ the Jones/SILVER n-body
call site first to match keys). Populate `center`, `moon`,
`n_maneuvers=len(route.leveraging_legs)`, per-leg ΔV + apse radius, and entry/final
V∞. Do NOT run the propagation here.

- [ ] **Step 4: Run test to verify it passes**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run pytest tests/data/test_discover_endgame.py -q
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run ruff check src/cyclerfinder/data/discover_novel.py tests/data/test_discover_endgame.py
uv run ruff format src/cyclerfinder/data/discover_novel.py tests/data/test_discover_endgame.py
git add src/cyclerfinder/data/discover_novel.py tests/data/test_discover_endgame.py
git status --short
git commit -m "data/discover_novel: endgame-route -> n-body harness request hook"
```

---

## Task 12: Full regression + lint + type gate

**Files:** none (verification only).

- [ ] **Step 1: Full suite (excluding slow)**

```bash
export PATH="$HOME/.local/bin:$PATH"
date -Iseconds
uv run pytest -q -m "not slow"
date -Iseconds
```
Expected: all PASS (no regressions).

- [ ] **Step 2: Lint + format + types**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run ruff check .
uv run ruff format --check .
uv run mypy src/cyclerfinder/search/leveraging_leg.py src/cyclerfinder/search/endgame_graph.py src/cyclerfinder/data/discover_novel.py
```
Expected: all clean. Fix anything introduced.

- [ ] **Step 3: Final commit (only if Step 2 required fixes)**

```bash
export PATH="$HOME/.local/bin:$PATH"
git add -p
git status --short
git commit -m "search/endgame: lint + type fixes for the endgame solver"
```

---

## Post-implementation (controller, not a task)

- Summarise outcome: did the leveraging genome reach the bend-feasible floor on
  any family, and via what route (legs + transfers)? EMPTY is a SUCCESS — a new
  method-versioned negative recorded alongside the ballistic ones.
- **No catalogue writeback** unless a survivor passed n-body AND has a same-model
  golden — that decision returns to the session, not this plan.
- Deferred rungs remain out of scope: real DE440 phase-full, broken-plane,
  heliocentric reuse, and Spec 2 (multi-arc DSM into the corrector).
