# ER3BP Discovery Campaign — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Run the first ER3BP discovery campaign — continue rotating-frame cycler families from e=0 into target eccentricity, Floquet-monitoring each for survival/death and bifurcations (candidate e>0-only families), report-only.

**Architecture:** New `src/cyclerfinder/search/er3bp_discovery.py` (+ a small Floquet helper) orchestrates the existing, unchanged #293 ER3BP genome (`core/er3bp.py`, `genome/er3bp_{periodic,continuation}.py`). Seeds come from already-encoded CR3BP ICs (Broucke golden, Koblick NRHO table) as the guaranteed floor, plus best-effort catalogue-cycler IC recovery. No catalogue writeback.

**Tech Stack:** Python 3.11, numpy, scipy (`solve_ivp` via the genome, `numpy.linalg.eig` for Floquet), the cyclerfinder core/genome/search packages, pytest, uv, ruff, mypy.

**Spec:** `docs/superpowers/specs/2026-06-24-er3bp-discovery-campaign-design.md`

**Conventions every task:** run `uv run ruff check <files> && uv run ruff format <files> && uv run mypy <files>` before each commit; never add Co-Authored-By trailers; work on `main` (no branches); pathspec `git add` (never `-A`); implementers must NOT spawn their own review subagents; imports at file top (ruff E402). Commit prefix `search/#432:` or `tests/#432:`.

**Confirmed upstream API (use verbatim; do not assume otherwise):**
- `ER3BPSystem(mu: float, e: float, primary_name: str, secondary_name: str)`; classmethod `ER3BPSystem.from_cr3bp(cr3bp: CR3BPSystem, e: float)`. In `cyclerfinder.core.er3bp`.
- `propagate_er3bp(state6, f_span: tuple[float,float], sys: ER3BPSystem, rtol=1e-12, atol=1e-12, with_stm=False, stm_mode="variable", method="DOP853") -> (f_eval, state_history, STM_final)`. With `with_stm=True`, `STM_final` is the 6×6 STM over `f_span`. In `cyclerfinder.core.er3bp`.
- `continue_er3bp_family_in_e(sys_base: ER3BPSystem, seed_state, period_f: float, e_target: float, n_steps: int, *, is_half_period_residual=True, tol=1e-10) -> list[ER3BPPeriodicOrbit]`. In `cyclerfinder.genome.er3bp_continuation`.
- `ER3BPPeriodicOrbit(state0, period_f, mu, e, corrector_residual, independent_residual, iterations, notes)`. In `cyclerfinder.genome.er3bp_periodic`.
- `cr3bp_system(primary: str, secondary: str) -> CR3BPSystem` (fields incl. `.mu`). In `cyclerfinder.core.cr3bp`.
- Encoded CR3BP ICs: `cyclerfinder.genome.tulip.KOBLICK_2023_TABLE4` (Earth-Moon NRHO by petal count); the Broucke-1969 Earth-Moon golden IC used in `tests/genome/test_er3bp_genome.py` (read that test for the exact `mu=0.0121550`, Orbit-1 `x0=0.1520965, ydot0=3.1608994`, `e=0.0001`, half-period `period_f=π`).

---

## File Structure

| File | Responsibility | Status |
|------|----------------|--------|
| `src/cyclerfinder/search/er3bp_floquet.py` | NEW. `er3bp_monodromy` + `floquet_classify`. | create |
| `src/cyclerfinder/search/er3bp_discovery.py` | NEW. `Er3bpSeed`, seed providers, `continue_and_monitor`, outcome classification, adjudication. | create |
| `tests/search/test_er3bp_floquet.py` | NEW. Floquet golden + classification. | create |
| `tests/search/test_er3bp_discovery.py` | NEW. Seeds, driver, classification, determinism. | create |
| `scripts/run_432_er3bp_discovery.py` | NEW. Phase-A/B runner. | create |

---

## Task 1: Floquet monitor (`er3bp_monodromy` + `floquet_classify`)

**Files:** Create `src/cyclerfinder/search/er3bp_floquet.py`, `tests/search/test_er3bp_floquet.py`

**Context:** The monodromy is the full-period (`f ∈ [0, period_f]`) STM. `propagate_er3bp(..., with_stm=True)` returns it as `STM_final`. Floquet eigenvalues classify stability: a periodic orbit is **stable** iff all eigenvalues lie on/inside the unit circle (|λ| ≤ 1 + tol); **unstable** if any |λ| > 1 + tol; an eigenvalue **on** the unit circle (||λ| − 1| ≤ tol) is the bifurcation signal.

- [ ] **Step 1: Write the failing monodromy-shape test**

```python
# tests/search/test_er3bp_floquet.py
from __future__ import annotations

import numpy as np

from cyclerfinder.core.er3bp import ER3BPSystem
from cyclerfinder.search.er3bp_floquet import er3bp_monodromy


def test_monodromy_is_6x6_real() -> None:
    # Broucke-1969 Earth-Moon Orbit 1 (mu=0.0121550, e=0.0001), half-period pi.
    sys = ER3BPSystem(mu=0.0121550, e=0.0001, primary_name="E", secondary_name="M")
    state0 = np.array([0.1520965, 0.0, 0.0, 0.0, 3.1608994, 0.0])
    m = er3bp_monodromy(state0, period_f=2.0 * np.pi, system=sys)
    assert m.shape == (6, 6)
    assert np.isfinite(m).all()
```

- [ ] **Step 2: Run → FAIL** (`cannot import name 'er3bp_monodromy'`). `uv run pytest tests/search/test_er3bp_floquet.py::test_monodromy_is_6x6_real -v`

- [ ] **Step 3: Implement `er3bp_monodromy` + `floquet_classify`**

```python
# src/cyclerfinder/search/er3bp_floquet.py
"""#432 Floquet monitor for the ER3BP discovery campaign.

The monodromy is the full-period (f in [0, period_f]) STM from propagate_er3bp;
its eigenvalues classify stability and flag bifurcations (eigenvalue on the unit
circle). Conventions mirror the #347 Floquet framework.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.er3bp import ER3BPSystem, propagate_er3bp

_UNIT_CIRCLE_TOL = 1.0e-3  # |λ| within this of 1.0 counts as "on the unit circle"


def er3bp_monodromy(
    state0: NDArray[np.float64], period_f: float, system: ER3BPSystem
) -> NDArray[np.float64]:
    """Full-period monodromy (6x6 STM over f in [0, period_f]) via propagate_er3bp."""
    _f, _hist, stm = propagate_er3bp(
        np.asarray(state0, dtype=np.float64),
        (0.0, period_f),
        system,
        with_stm=True,
    )
    return np.asarray(stm, dtype=np.float64)


@dataclass(frozen=True)
class FloquetResult:
    eigenvalues: tuple[complex, ...]
    stability_tag: str  # "stable" | "unstable" | "marginal"
    on_unit_circle: bool  # a non-trivial eigenvalue sits on the unit circle


def floquet_classify(
    monodromy: NDArray[np.float64], *, unit_circle_tol: float = _UNIT_CIRCLE_TOL
) -> FloquetResult:
    """Classify stability from the monodromy eigenvalues.

    stable: all |λ| <= 1 + tol. unstable: some |λ| > 1 + tol. marginal: max |λ|
    is within tol of 1 (no eigenvalue clearly outside). on_unit_circle flags a
    non-trivial eigenvalue (|λ-1| > tol away from the trivial pair) sitting on
    the unit circle — the bifurcation signal.
    """
    eig = np.linalg.eigvals(np.asarray(monodromy, dtype=np.float64))
    mags = np.abs(eig)
    max_mag = float(mags.max())
    if max_mag > 1.0 + unit_circle_tol:
        tag = "unstable"
    elif max_mag < 1.0 - unit_circle_tol:
        tag = "stable"
    else:
        tag = "marginal"
    # bifurcation: any eigenvalue ON the unit circle that is not the trivial λ=1 pair.
    on_uc = bool(
        np.any(
            (np.abs(mags - 1.0) <= unit_circle_tol) & (np.abs(eig - 1.0) > unit_circle_tol)
        )
    )
    return FloquetResult(
        eigenvalues=tuple(complex(v) for v in eig), stability_tag=tag, on_unit_circle=on_uc
    )
```

- [ ] **Step 4: Run → PASS.** `uv run pytest tests/search/test_er3bp_floquet.py::test_monodromy_is_6x6_real -v`

- [ ] **Step 5: Write the classification golden test**

```python
# tests/search/test_er3bp_floquet.py  (append)
from cyclerfinder.search.er3bp_floquet import floquet_classify


def test_floquet_classify_on_broucke_orbit() -> None:
    """The Broucke Orbit-1 monodromy classifies to a definite stability tag and
    its eigenvalues include the trivial λ≈1 pair (Hamiltonian monodromy)."""
    sys = ER3BPSystem(mu=0.0121550, e=0.0001, primary_name="E", secondary_name="M")
    state0 = np.array([0.1520965, 0.0, 0.0, 0.0, 3.1608994, 0.0])
    m = er3bp_monodromy(state0, period_f=2.0 * np.pi, system=sys)
    res = floquet_classify(m)
    assert res.stability_tag in {"stable", "unstable", "marginal"}
    # Hamiltonian/symplectic monodromy: eigenvalues come in reciprocal pairs, so
    # the product of magnitudes is ~1; at least one eigenvalue near the unit circle.
    assert np.isclose(float(np.prod(np.abs(res.eigenvalues))), 1.0, atol=1e-2)
```

- [ ] **Step 6: Run → PASS.** If the `prod≈1` symplectic check fails by more than atol (integration tolerance), tighten `propagate_er3bp` rtol/atol in `er3bp_monodromy` to 1e-13 and re-run; do NOT loosen the test below atol=1e-2 without noting why.

- [ ] **Step 7: Commit**

```bash
uv run ruff check src/cyclerfinder/search/er3bp_floquet.py tests/search/test_er3bp_floquet.py
uv run ruff format src/cyclerfinder/search/er3bp_floquet.py tests/search/test_er3bp_floquet.py
uv run mypy src/cyclerfinder/search/er3bp_floquet.py
git add src/cyclerfinder/search/er3bp_floquet.py tests/search/test_er3bp_floquet.py
git commit -m "search/#432: ER3BP Floquet monitor (monodromy + stability classification)"
```

---

## Task 2: Seed model + standard-family seed provider (the floor)

**Files:** Create `src/cyclerfinder/search/er3bp_discovery.py`, `tests/search/test_er3bp_discovery.py`

**Context:** The guaranteed seed floor is ICs already encoded in the repo: the Broucke-1969 Earth-Moon family (golden) and the Koblick NRHO table (`tulip.KOBLICK_2023_TABLE4`, Earth-Moon). These are CR3BP ICs at e=0 that the ER3BP continuation lifts into e>0. `Er3bpSeed` is the unit the campaign iterates over.

- [ ] **Step 1: Write the failing seed-floor test**

```python
# tests/search/test_er3bp_discovery.py
from __future__ import annotations

from cyclerfinder.search.er3bp_discovery import Er3bpSeed, standard_family_seeds


def test_standard_family_seeds_returns_usable_floor() -> None:
    seeds = standard_family_seeds(target_e=0.0549)
    assert len(seeds) >= 1
    for s in seeds:
        assert isinstance(s, Er3bpSeed)
        assert s.state0.shape == (6,)
        assert s.period_f > 0.0
        assert 0.0 < s.target_e < 1.0
        assert s.system.primary_name and s.system.secondary_name
        assert s.source  # non-empty provenance
```

- [ ] **Step 2: Run → FAIL** (`cannot import name`). `uv run pytest tests/search/test_er3bp_discovery.py::test_standard_family_seeds_returns_usable_floor -v`

- [ ] **Step 3: Implement `Er3bpSeed` + `standard_family_seeds`**

```python
# src/cyclerfinder/search/er3bp_discovery.py
"""#432 ER3BP discovery campaign: continue rotating-frame cycler families into
e>0, Floquet-monitor for survival/death and bifurcations. Report-only; no
catalogue writeback. See docs/superpowers/specs/2026-06-24-er3bp-discovery-campaign-design.md.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.er3bp import ER3BPSystem


@dataclass(frozen=True)
class Er3bpSeed:
    label: str
    system: ER3BPSystem  # the e=0 CR3BP μ lives in system.mu; system.e is the target
    state0: NDArray[np.float64]  # rotating-frame IC at e=0, shape (6,)
    period_f: float  # true-anomaly period (multiple of 2π for full period)
    is_half_period_residual: bool
    target_e: float
    source: str  # provenance string


# Broucke (1969) TR 32-1360 Table 12, Family 7P, Earth-Moon mu=0.0121550, Orbit 1.
_BROUCKE_EM_MU = 0.0121550
_BROUCKE_EM_ORBIT1 = np.array([0.1520965, 0.0, 0.0, 0.0, 3.1608994, 0.0])


def standard_family_seeds(*, target_e: float = 0.0549) -> list[Er3bpSeed]:
    """Guaranteed seed floor from ICs already encoded in the repo (Earth-Moon).

    Currently the sourced Broucke-1969 Earth-Moon family. The Koblick NRHO table
    (tulip.KOBLICK_2023_TABLE4) can be added here as additional Earth-Moon seeds
    once their (x0, z0, ydot0, period) are mapped to the (state0, period_f) form.
    """
    sys = ER3BPSystem(
        mu=_BROUCKE_EM_MU, e=target_e, primary_name="E", secondary_name="M"
    )
    return [
        Er3bpSeed(
            label="broucke-1969-em-7P-orbit1",
            system=sys,
            state0=_BROUCKE_EM_ORBIT1.copy(),
            period_f=2.0 * np.pi,
            is_half_period_residual=True,
            target_e=target_e,
            source="Broucke 1969 TR 32-1360 Table 12 Family 7P Orbit 1 (mu=0.0121550)",
        )
    ]
```

- [ ] **Step 4: Run → PASS.** `uv run pytest tests/search/test_er3bp_discovery.py::test_standard_family_seeds_returns_usable_floor -v`

- [ ] **Step 5: Commit**

```bash
uv run ruff check src/cyclerfinder/search/er3bp_discovery.py tests/search/test_er3bp_discovery.py
uv run ruff format src/cyclerfinder/search/er3bp_discovery.py tests/search/test_er3bp_discovery.py
uv run mypy src/cyclerfinder/search/er3bp_discovery.py
git add src/cyclerfinder/search/er3bp_discovery.py tests/search/test_er3bp_discovery.py
git commit -m "search/#432: Er3bpSeed model + standard-family seed floor (Broucke EM golden)"
```

---

## Task 3: Continuation driver + outcome classification

**Files:** Modify `src/cyclerfinder/search/er3bp_discovery.py`, `tests/search/test_er3bp_discovery.py`

**Context:** `continue_and_monitor` continues a seed e=0→target_e via `continue_er3bp_family_in_e` (which returns one `ER3BPPeriodicOrbit` per step), computes the monodromy+Floquet at each returned step, and classifies the seed: **survives** (reaches target_e), **dies** (continuation returns fewer than `n_steps` orbits / raises before target_e), **bifurcates** (a step's `floquet_classify` reports `on_unit_circle`, recording the first such `e`).

- [ ] **Step 1: Write the failing survives-classification test**

```python
# tests/search/test_er3bp_discovery.py  (append)
from cyclerfinder.search.er3bp_discovery import continue_and_monitor, standard_family_seeds


def test_continue_and_monitor_classifies_survival() -> None:
    """The Broucke EM floor seed continues to a small target e and is classified
    (survives|dies|bifurcates) with a per-step trace."""
    seed = standard_family_seeds(target_e=0.01)[0]
    trace = continue_and_monitor(seed, n_steps=5)
    assert trace.outcome in {"survives", "dies", "bifurcates"}
    assert len(trace.steps) >= 1
    # each step records (e, residual, stability_tag)
    s0 = trace.steps[0]
    assert 0.0 <= s0.e <= 0.01 + 1e-9
    assert s0.stability_tag in {"stable", "unstable", "marginal", "unknown"}
```

- [ ] **Step 2: Run → FAIL** (`cannot import name 'continue_and_monitor'`).

- [ ] **Step 3: Implement the trace types + driver**

```python
# src/cyclerfinder/search/er3bp_discovery.py  (append; imports at top)
from cyclerfinder.genome.er3bp_continuation import continue_er3bp_family_in_e
from cyclerfinder.search.er3bp_floquet import er3bp_monodromy, floquet_classify


@dataclass(frozen=True)
class Er3bpStep:
    e: float
    corrector_residual: float
    stability_tag: str
    on_unit_circle: bool


@dataclass(frozen=True)
class Er3bpContinuationTrace:
    seed_label: str
    outcome: str  # "survives" | "dies" | "bifurcates"
    steps: tuple[Er3bpStep, ...]
    e_max_reached: float
    e_star: float | None  # first bifurcation eccentricity, if any
    target_e: float


def continue_and_monitor(seed: Er3bpSeed, *, n_steps: int = 20) -> Er3bpContinuationTrace:
    """Continue ``seed`` e=0→target_e, Floquet-monitor each step, classify."""
    try:
        family = continue_er3bp_family_in_e(
            ER3BPSystem(
                mu=seed.system.mu, e=0.0,
                primary_name=seed.system.primary_name,
                secondary_name=seed.system.secondary_name,
            ),
            seed.state0,
            seed.period_f,
            seed.target_e,
            n_steps,
            is_half_period_residual=seed.is_half_period_residual,
        )
    except Exception:
        family = []
    steps: list[Er3bpStep] = []
    e_star: float | None = None
    for orb in family:
        sys_e = ER3BPSystem(
            mu=orb.mu, e=orb.e,
            primary_name=seed.system.primary_name,
            secondary_name=seed.system.secondary_name,
        )
        try:
            mono = er3bp_monodromy(orb.state0, orb.period_f, sys_e)
            fl = floquet_classify(mono)
            tag, on_uc = fl.stability_tag, fl.on_unit_circle
        except Exception:
            tag, on_uc = "unknown", False
        steps.append(
            Er3bpStep(e=orb.e, corrector_residual=orb.corrector_residual,
                      stability_tag=tag, on_unit_circle=on_uc)
        )
        if on_uc and e_star is None:
            e_star = orb.e
    e_max = steps[-1].e if steps else 0.0
    reached_target = bool(steps) and abs(e_max - seed.target_e) <= 1e-6
    if e_star is not None:
        outcome = "bifurcates"
    elif reached_target:
        outcome = "survives"
    else:
        outcome = "dies"
    return Er3bpContinuationTrace(
        seed_label=seed.label, outcome=outcome, steps=tuple(steps),
        e_max_reached=e_max, e_star=e_star, target_e=seed.target_e,
    )
```

- [ ] **Step 4: Run → PASS.** `uv run pytest tests/search/test_er3bp_discovery.py::test_continue_and_monitor_classifies_survival -v` (allow ~30s — integrates the ER3BP STM per step).

- [ ] **Step 5: Write the dies-classification test**

```python
# tests/search/test_er3bp_discovery.py  (append)
import numpy as np

from cyclerfinder.core.er3bp import ER3BPSystem
from cyclerfinder.search.er3bp_discovery import Er3bpSeed, continue_and_monitor


def test_continue_and_monitor_dies_on_infeasible_seed() -> None:
    """A garbage IC cannot continue and is classified 'dies' without crashing."""
    seed = Er3bpSeed(
        label="garbage", system=ER3BPSystem(mu=0.0121550, e=0.5, primary_name="E", secondary_name="M"),
        state0=np.array([5.0, 5.0, 0.0, 9.0, 9.0, 0.0]), period_f=2.0 * np.pi,
        is_half_period_residual=True, target_e=0.5, source="test-garbage",
    )
    trace = continue_and_monitor(seed, n_steps=3)
    assert trace.outcome == "dies"
    assert trace.e_max_reached < 0.5
```

- [ ] **Step 6: Run → PASS.** If the garbage IC happens to converge (unlikely), pick a more clearly-divergent IC (e.g. very large velocities) so the corrector provably fails; the point is the no-crash death path.

- [ ] **Step 7: Commit**

```bash
uv run ruff check src/cyclerfinder/search/er3bp_discovery.py tests/search/test_er3bp_discovery.py
uv run ruff format src/cyclerfinder/search/er3bp_discovery.py tests/search/test_er3bp_discovery.py
uv run mypy src/cyclerfinder/search/er3bp_discovery.py
git add src/cyclerfinder/search/er3bp_discovery.py tests/search/test_er3bp_discovery.py
git commit -m "search/#432: continuation driver + survives/dies/bifurcates classification"
```

---

## Task 4: Catalogue-cycler best-effort seed provider

**Files:** Modify `src/cyclerfinder/search/er3bp_discovery.py`, `tests/search/test_er3bp_discovery.py`

**Context:** The 12 CR3BP-class catalogue rows store NO rotating-frame ICs (verified). `catalogue_cr3bp_seeds()` recovers ICs only where a code constant / sourced value exists, and SKIPS the rest with a logged count (no silent drop). For this first pass the recoverable set may be just the Koblick NRHO table (Earth-Moon) mapped from `tulip.KOBLICK_2023_TABLE4`; the function must return a (possibly empty) list and a skip count, and the campaign floor must not depend on it.

- [ ] **Step 1: Write the failing best-effort test**

```python
# tests/search/test_er3bp_discovery.py  (append)
from cyclerfinder.search.er3bp_discovery import catalogue_cr3bp_seeds


def test_catalogue_cr3bp_seeds_returns_list_and_skip_count() -> None:
    seeds, skipped = catalogue_cr3bp_seeds(target_e=0.0549)
    assert isinstance(seeds, list)
    assert isinstance(skipped, int)
    assert skipped >= 0
    for s in seeds:
        assert s.state0.shape == (6,)
        assert s.source
```

- [ ] **Step 2: Run → FAIL** (`cannot import name`).

- [ ] **Step 3: Implement `catalogue_cr3bp_seeds`**

```python
# src/cyclerfinder/search/er3bp_discovery.py  (append; imports at top)
from cyclerfinder.genome.tulip import EARTH_MOON_MU, KOBLICK_2023_TABLE4


def catalogue_cr3bp_seeds(*, target_e: float = 0.0549) -> tuple[list[Er3bpSeed], int]:
    """Best-effort ICs for CR3BP-class catalogue cyclers.

    The catalogue rows store no rotating-frame ICs, so we recover only where an
    encoded constant exists. Currently: the Koblick NRHO table (Earth-Moon).
    Returns (seeds, n_skipped). The campaign floor (standard_family_seeds) does
    NOT depend on this.
    """
    seeds: list[Er3bpSeed] = []
    skipped = 0
    for np_count, row in KOBLICK_2023_TABLE4.items():
        x0 = row.get("x0")
        z0 = row.get("z0")
        ydot0 = row.get("ydot0") or row.get("vy0")
        period = row.get("period") or row.get("period_tu")
        if x0 is None or ydot0 is None or period is None:
            skipped += 1
            continue
        state0 = np.array([float(x0), 0.0, float(z0 or 0.0), 0.0, float(ydot0), 0.0])
        sys = ER3BPSystem(
            mu=float(EARTH_MOON_MU), e=target_e, primary_name="E", secondary_name="M"
        )
        seeds.append(
            Er3bpSeed(
                label=f"koblick-nrho-np{np_count}", system=sys, state0=state0,
                period_f=2.0 * np.pi, is_half_period_residual=True, target_e=target_e,
                source=f"Koblick 2023 Table 4 NRHO Np={np_count}",
            )
        )
    return seeds, skipped
```

NOTE: confirm the real key names in `KOBLICK_2023_TABLE4` (`x0`/`z0`/`ydot0`/`period`) by reading `tulip.py:134`; adapt the `.get(...)` keys to the actual schema, and confirm `EARTH_MOON_MU` is the exported constant name (else read the system's μ via `cr3bp_system("E","M").mu`). The skip-on-missing-key logic is the load-bearing part. If the NRHO period is in time units not true-anomaly, set `period_f = 2π` (one synodic-frame revolution maps to 2π in the ER3BP true-anomaly clock for the continuation's purposes) and note it.

- [ ] **Step 4: Run → PASS.**

- [ ] **Step 5: Commit**

```bash
uv run ruff check src/cyclerfinder/search/er3bp_discovery.py tests/search/test_er3bp_discovery.py
uv run ruff format src/cyclerfinder/search/er3bp_discovery.py tests/search/test_er3bp_discovery.py
uv run mypy src/cyclerfinder/search/er3bp_discovery.py
git add src/cyclerfinder/search/er3bp_discovery.py tests/search/test_er3bp_discovery.py
git commit -m "search/#432: best-effort catalogue-cycler seed recovery (skip-logged)"
```

---

## Task 5: Survivor adjudication (literature + ML flag)

**Files:** Modify `src/cyclerfinder/search/er3bp_discovery.py`, `tests/search/test_er3bp_discovery.py`

**Context:** Survivors and bifurcation-flagged seeds get a literature verdict. Reuse `check_literature`. The structural fingerprint for an ER3BP orbit is its (system, μ, e, period) — build a `CandidateSignature` consistent with how `literature_check` expects it (READ `cyclerfinder.search.literature_check` for the `CandidateSignature` fields + `check_literature(sig, search=...)` signature, and mirror the precursor_matcher usage). When no search fn is supplied, return the deferred/inconclusive sentinel (mirror the local matcher path).

- [ ] **Step 1: Write the failing adjudication test**

```python
# tests/search/test_er3bp_discovery.py  (append)
from cyclerfinder.search.er3bp_discovery import adjudicate_trace


def test_adjudicate_trace_deferred_without_search_fn() -> None:
    seed = standard_family_seeds(target_e=0.01)[0]
    trace = continue_and_monitor(seed, n_steps=3)
    verdict = adjudicate_trace(trace, seed, literature_check_search=None)
    assert verdict.status in {"inconclusive", "not-found", "published"}
```

- [ ] **Step 2: Run → FAIL.**

- [ ] **Step 3: Implement `adjudicate_trace`** mirroring the precursor_matcher literature handling: build a `CandidateSignature` from `(seed.system.primary_name, seed.system.secondary_name, seed.system.mu, trace.target_e, seed.period_f)`; if `literature_check_search is None` return a `LiteratureCheckResult(status="inconclusive", ..., notes="No literature_check_search; verdict deferred")`; else `return check_literature(sig, search=literature_check_search)`. Return type is `LiteratureCheckResult`. Match the exact `CandidateSignature` constructor fields from `literature_check.py`.

- [ ] **Step 4: Run → PASS.**

- [ ] **Step 5: Commit**

```bash
uv run ruff check src/cyclerfinder/search/er3bp_discovery.py tests/search/test_er3bp_discovery.py
uv run ruff format src/cyclerfinder/search/er3bp_discovery.py tests/search/test_er3bp_discovery.py
uv run mypy src/cyclerfinder/search/er3bp_discovery.py
git add src/cyclerfinder/search/er3bp_discovery.py tests/search/test_er3bp_discovery.py
git commit -m "search/#432: survivor literature adjudication"
```

---

## Task 6: Campaign runner script

**Files:** Create `scripts/run_432_er3bp_discovery.py`

**Context:** Phase A (standard floor + best-effort catalogue seeds at native real e) then Phase B (standard families at high-e systems — for the first cut, re-run the Earth-Moon floor seeds at the high-e targets Sun-Mercury μ / Sun-Mars μ by building `ER3BPSystem` from `cr3bp_system("Sun","Me")` / `("Sun","M")` and a generic Lyapunov/DRO seed; if no high-e CR3BP seed IC is available yet, log that Phase B is seed-limited and run whatever standard seeds exist). Incremental timestamped logging; JSONL output per seed.

- [ ] **Step 1: Write the runner** mirroring `scripts/run_430_global_precursor.py` structure: `_print_progress` timestamped logging; iterate seeds from `standard_family_seeds()` + `catalogue_cr3bp_seeds()` (Phase A) then high-e variants (Phase B); for each call `continue_and_monitor` + `adjudicate_trace`; write `data/er3bp_discovery_{phaseA,phaseB}.jsonl` (one JSON record per seed: label, system, outcome, e_max_reached, e_star, per-step stability tags, literature status). Print per-phase outcome breakdown (survives/dies/bifurcates counts) + any bifurcation-flagged seeds.

- [ ] **Step 2: Smoke-test** (small n_steps, the floor seed only):

```bash
timeout 600 uv run python -c "
from cyclerfinder.search.er3bp_discovery import standard_family_seeds, continue_and_monitor, adjudicate_trace
s=standard_family_seeds(target_e=0.0549)[0]
t=continue_and_monitor(s, n_steps=8)
print('SMOKE outcome=', t.outcome, 'e_max=', round(t.e_max_reached,4), 'e_star=', t.e_star)
print('  steps:', [(round(st.e,3), st.stability_tag) for st in t.steps])
print('  lit=', adjudicate_trace(t, s, literature_check_search=None).status)
"
```
Confirm it classifies the Broucke EM floor seed to a definite outcome with a per-step stability trace and no crash.

- [ ] **Step 3: Commit**

```bash
uv run ruff check scripts/run_432_er3bp_discovery.py
uv run ruff format scripts/run_432_er3bp_discovery.py
git add scripts/run_432_er3bp_discovery.py
git commit -m "search/#432: ER3BP discovery campaign runner (Phase A/B)"
```

---

## Task 7: Deliverable run + verdict + registry

**Files:** Create `data/er3bp_discovery_{phaseA,phaseB}.jsonl`, `docs/superpowers/plans/2026-06-24-432-er3bp-discovery-verdict.md`; modify `data/empty_regions.jsonl`

**Context:** The controller (not an implementer subagent) runs the detached campaign, harvests, writes the verdict + registry. This task documents that hand-off.

- [ ] **Step 1: Launch detached** `nohup uv run python scripts/run_432_er3bp_discovery.py > data/runlogs/432_er3bp.log 2>&1 &` and arm a harness-tracked waiter on the PID.
- [ ] **Step 2: Harvest** the two JSONLs → outcome breakdown per system (survives/dies/bifurcates), the bifurcation-flagged seeds (`e_star`), and literature verdicts.
- [ ] **Step 3: Write the verdict note** — what survives to real e, where bifurcations are flagged (candidate e>0-only families), literature status. Honest framing: most likely a structural map; a clean "no bifurcation flagged" is a registry-grade negative; any flagged bifurcation → follow-on branch-switching + ER3BP gauntlet (new tasks, not built here).
- [ ] **Step 4: Register** the method-versioned result in `data/empty_regions.jsonl` (first ER3BP discovery capability; record the seed set, systems, n_steps, git sha). Run `uv run pytest tests/data -q -k "empty_region or registry"` before committing.
- [ ] **Step 5: Commit** the script outputs, verdict note, registry entry (pathspec).

---

## Self-Review

**Spec coverage:** Floquet monitor → Task 1 ✓; seed registry (floor + catalogue best-effort) → Tasks 2+4 ✓; continuation driver + survives/dies/bifurcates → Task 3 ✓; adjudication → Task 5 ✓; runner → Task 6 ✓; deliverable run + verdict + registry → Task 7 ✓; report-only / no catalogue writeback → Task 7 (registry + verdict only). ✓

**Placeholder scan:** Tasks 4 and 5 flag exact upstream schemas to confirm before writing (`KOBLICK_2023_TABLE4` keys, `EARTH_MOON_MU` name, `CandidateSignature` fields) — these are verify-then-adapt directives, not placeholders; the load-bearing logic (skip-on-missing, deferred-sentinel) is fully specified. No "TBD"/"handle edge cases".

**Type consistency:** `Er3bpSeed` (Task 2) consumed by `continue_and_monitor` (Task 3), `catalogue_cr3bp_seeds` (Task 4), `adjudicate_trace` (Task 5), runner (Task 6). `Er3bpContinuationTrace`/`Er3bpStep` (Task 3) consumed by `adjudicate_trace` + runner. `FloquetResult` (Task 1: `.stability_tag`, `.on_unit_circle`) consumed by Task 3's driver. `er3bp_monodromy(state0, period_f, system)` signature consistent between Task 1 def and Task 3 call. Consistent.

**Discipline:** every golden sourced (Broucke 1969); report-only, no catalogue writeback; skip-counts logged (no silent truncation); ER3BP genome reused unedited.
