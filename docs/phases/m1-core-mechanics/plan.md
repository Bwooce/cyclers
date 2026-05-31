# M1 — Core mechanics

**Spec reference:** spec.md §4 (architecture), §6 (interfaces sketch), §7 (tech stack), §8 (M1 milestone definition), §9 (validation anchors), §10 (risks), §12(b) (multi-rev Lambert list-return interface — supersedes §6).

**Purpose:** deliver the three pieces of pure orbital mechanics that every subsequent milestone consumes — a planet-state provider (`ephemeris`), a Lambert two-point boundary-value solver (`lambert`), and a universal-variable Kepler propagator (`kepler`). M0 produced the ground; M1 produces the *primitives*. No flyby logic, no Tisserand graphs, no `Cycler` model — those land in M2/M3.

**Gate (definition of done):** `uv run pytest` shows the Lambert solver agrees with `lamberthub`'s `izzo2015` and `gooding1990` to `max_diff_mps < 1e-3` on **at least three distinct test legs** (short/medium/long ToF regimes), AND a Lambert solution's `(v1, v2)` re-propagated forward by `tof` using `kepler.propagate(r1, v1, tof)` lands within < 1 km of `r2`. `uv run ruff check`, `uv run ruff format --check`, and `uv run mypy src tests` remain clean. CI green.

---

## 1. What this milestone delivers

A package that, in addition to M0's scaffold and constants:

1. **`cyclerfinder.core.ephemeris.Ephemeris`** — a class taking `model: str = "circular"` whose `state(body, t_sec)` returns `(r_km, v_km_s)` heliocentric inertial vectors. Only the `"circular"` backend has a working implementation in M1; `"astropy"` raises `NotImplementedError` with an explicit pointer to M6. The class is structured so M6 swaps in the JPL backend by overriding a single internal method.
2. **`cyclerfinder.core.lambert.lambert`** — a self-contained universal-variable single-revolution Lambert solver returning `list[LambertSolution]` per spec §12(b). The multi-rev branches (`n_revs >= 1`) return an empty list in M1; this is documented as the M4 extension point. Edge cases (180° transfer, retrograde) have explicit, tested behaviour.
3. **`cyclerfinder.core.lambert.lambert_crosscheck`** — a developer-facing convenience that solves a leg with the in-house solver and with `lamberthub.izzo2015` + `lamberthub.gooding1990`, returning a dict with all three answers and the `max_diff_mps`. Used by the gate tests now and by M3's V0/V1 closure verification later.
4. **`cyclerfinder.core.kepler.propagate`** — a universal-variable two-body propagator using Stumpff functions and Newton iteration on the universal anomaly χ. Used by tests to re-propagate Lambert solutions, and by later phases to sample trajectories and verify closure.
5. A `tests/core/` test suite that exercises each module independently and the Lambert↔Kepler self-consistency described in §4.

Explicitly **out of scope** for M1 (these belong to later milestones, do not stub):

- Multi-rev Lambert branch selection (M4).
- `astropy`/JPL ephemeris backend implementation (M6).
- Rotating-frame transforms (M3, `frames.py`).
- Flyby, Tisserand, resonance (M2).
- `Cycler`/`Leg`/`Encounter` dataclasses, closure residual (M3).
- Any CLI or viz code.

---

## 2. File tree after M1

```
cyclers/
├── .github/workflows/ci.yml
├── .gitignore
├── .python-version
├── README.md
├── pyproject.toml                          # updated: +scipy (runtime), +lamberthub (dev)
├── uv.lock                                 # regenerated
├── docs/
│   └── phases/
│       ├── m0-scaffold/{plan,todo}.md
│       └── m1-core-mechanics/
│           ├── plan.md                     # this file
│           └── todo.md
├── src/
│   └── cyclerfinder/
│       ├── __init__.py
│       └── core/
│           ├── __init__.py
│           ├── constants.py                # unchanged from M0
│           ├── ephemeris.py                # NEW
│           ├── lambert.py                  # NEW
│           └── kepler.py                   # NEW
└── tests/
    ├── __init__.py
    ├── test_constants.py                   # unchanged
    └── core/
        ├── __init__.py
        ├── test_ephemeris.py               # NEW
        ├── test_kepler.py                  # NEW
        ├── test_lambert.py                 # NEW — includes the gate tests
        └── conftest.py                     # NEW — shared planet-state fixtures
```

`search/`, `model/`, `verify/`, `data/`, `viz/` directories still do not exist — they appear when their first module lands (per the M0-established *what exists = what works* convention).

---

## 3. Module designs

### 3.1 `core/ephemeris.py`

#### 3.1.1 Public API

```python
from typing import Protocol
import numpy as np
from numpy.typing import NDArray

Vec3 = NDArray[np.float64]   # shape (3,), dtype float64

class Ephemeris:
    def __init__(self, model: str = "circular") -> None: ...
    def state(self, body: str, t_sec: float) -> tuple[Vec3, Vec3]: ...
```

Returns `(r_km, v_km_s)`, both `numpy.ndarray` of shape `(3,)`, `dtype=float64`, in a heliocentric inertial frame with `+x` toward the J2000 vernal equinox and `+z` along the ecliptic north pole. The circular model places every planet in the `z=0` plane, so the returned `r[2]` and `v[2]` are exactly zero for the circular backend.

#### 3.1.2 Construction logic

`__init__` validates `model in {"circular", "astropy"}` and dispatches via an internal strategy pattern: `self._backend = _CircularBackend()` or raises `NotImplementedError("astropy backend lands in M6; use model='circular' for now")`. Public `state()` simply delegates to `self._backend.state(body, t_sec)`. This keeps M6's surgery to a single new class plus a one-line dispatch update.

#### 3.1.3 Circular model

For body code `b` (one of `"V"`, `"E"`, `"M"` from `constants.PLANETS`):

- Semi-major axis `a = PLANETS[b].sma_au * AU_KM`.
- Mean motion `n = PLANETS[b].mean_motion_deg_day * (π/180) / SECONDS_PER_DAY` (rad/s). Use the tabulated value, not `sqrt(μ/a³)` derived from it, to keep `constants.py` the single source of truth.
- True anomaly = mean anomaly = `θ = n * t_sec` (circular, so M = E = ν). Reference epoch is `t_sec = 0` → body at `θ = 0` (i.e. `+x` direction). M1 does not need a J2000 phase calibration; M3's Aldrin reproduction will choose its own epoch, and M6's astropy backend brings real phases.
- Position `r = a * [cos θ, sin θ, 0]`.
- Velocity `v = a * n * [-sin θ, cos θ, 0]`. (Circular orbital speed = `a · n`.)

All planets share the prograde sense (ccw viewed from `+z`). Unknown body codes raise `KeyError(body)` (let the dict do the talking).

#### 3.1.4 What this gets us

A deterministic, fast planet-state source the Lambert and Kepler tests can lean on without needing a real ephemeris. The M2 Tisserand maps and M3 Aldrin reproduction work entirely on this backend; the spec's §12.1 phase-matching only needs the real ephemeris at M6.

### 3.2 `core/lambert.py`

#### 3.2.1 Public API (per spec §12(b))

```python
from dataclasses import dataclass
import numpy as np
from numpy.typing import NDArray

Vec3 = NDArray[np.float64]

@dataclass(frozen=True)
class LambertSolution:
    n_revs: int
    branch: str          # "single" for n_revs == 0; "low" | "high" for n_revs >= 1
    v1: Vec3             # heliocentric velocity at r1, km/s
    v2: Vec3             # heliocentric velocity at r2, km/s

def lambert(
    r1: Vec3,
    r2: Vec3,
    tof: float,
    *,
    mu: float = MU_SUN_KM3_S2,
    prograde: bool = True,
    max_revs: int = 0,
) -> list[LambertSolution]: ...

def lambert_crosscheck(
    r1: Vec3,
    r2: Vec3,
    tof: float,
    *,
    mu: float = MU_SUN_KM3_S2,
    prograde: bool = True,
) -> dict[str, object]:
    """Returns {"mine": LambertSolution, "izzo": (v1, v2), "gooding": (v1, v2),
                "max_diff_mps": float}. Single-rev only in M1."""
```

| Parameter | Type | Meaning |
|---|---|---|
| `r1`, `r2` | `Vec3` (km) | departure / arrival position vectors, same inertial frame |
| `tof` | float (s) | time of flight, must be > 0 |
| `mu` | float (km³/s²) | central body gravitational parameter; defaults to `MU_SUN_KM3_S2` |
| `prograde` | bool | True ⇒ short-way for non-collinear `r1, r2` with positive `(r1 × r2) · ẑ`; False ⇒ retrograde |
| `max_revs` | int ≥ 0 | maximum revolution count searched |

**M1 behaviour:** the returned list contains **at most one** `LambertSolution` with `n_revs=0, branch="single"`. For `max_revs >= 1` the multi-rev branches are not yet computed — the list contains only the single-rev solution. This is documented in the function's docstring and tested. M4 fills in the multi-rev branches; the return type does not change.

#### 3.2.2 Algorithm: universal-variable Lambert (Bate/Mueller/White, Vallado Algorithm 5.2)

Self-contained, no scipy in the inner loop. Sketch:

1. Compute `r1_n = ||r1||`, `r2_n = ||r2||`, transfer angle `Δν` from `cos Δν = (r1·r2)/(r1_n·r2_n)` with sign from `prograde` and `(r1×r2)·ẑ`.
2. Compute `A = sin Δν · sqrt(r1_n · r2_n / (1 − cos Δν))`. If `A == 0` (i.e. 180° transfer) raise `LambertGeometryError("180° transfer is singular in the planar universal-variable form; use multi-rev or perturb the geometry")`. (See §3.2.4 below for the rationale.)
3. Newton iteration on universal anomaly `z` (= `χ² · α`):
   - `S(z), C(z)` are Stumpff functions implemented inline (series expansion near `z=0`, closed form otherwise) — share these with `kepler.py` by importing from a single private module `_stumpff.py`.
   - `y(z) = r1_n + r2_n + A · (z · S(z) − 1) / sqrt(C(z))`.
   - `t(z) = (y(z) / C(z))^(3/2) · S(z) / sqrt(μ) + A · sqrt(y(z) / μ)`.
   - Iterate `z` to make `t(z) == tof`; use Newton with the analytic derivative `dt/dz` (Vallado eq. 7-15).
4. Compute Lagrange coefficients `f, g, ġ` from converged `z, y`. Then `v1 = (r2 − f·r1) / g`, `v2 = (ġ·r2 − r1) / g`.
5. Return `[LambertSolution(n_revs=0, branch="single", v1=v1, v2=v2)]`.

Convergence tolerance: `|t(z) − tof| / tof < 1e-12` OR `|Δz| < 1e-12`. Iteration cap 60; failure raises `LambertConvergenceError(z, residual)`.

A small bracketing/bootstrap step is needed because `y(z) < 0` for some initial `z`: if a Newton step would land in the forbidden region, bisect the previous bracket instead. Standard Vallado workaround — documented inline.

#### 3.2.3 Cross-check helper

`lambert_crosscheck` calls:

- `lamberthub.izzo2015(mu, r1, r2, tof, M=0, prograde=prograde)` → `(v1_izzo, v2_izzo)`.
- `lamberthub.gooding1990(mu, r1, r2, tof, M=0, prograde=prograde)` → `(v1_g, v2_g)`.
- `mine = lambert(r1, r2, tof, mu=mu, prograde=prograde)[0]`.

`max_diff_mps = 1000 * max(||mine.v1 − v1_izzo||, ||mine.v2 − v2_izzo||, ||mine.v1 − v1_g||, ||mine.v2 − v2_g||)`.

This is the single source of truth that the gate test consumes. `lamberthub` is a **dev** dep — the production code path doesn't import it.

#### 3.2.4 Edge cases & errors

| Case | M1 behaviour |
|---|---|
| `tof <= 0` | `ValueError("tof must be positive, got ...")` |
| `||r1|| == 0` or `||r2|| == 0` | `ValueError("r1, r2 must be non-zero position vectors")` |
| `r1` and `r2` collinear (Δν ≈ 0 or 180°) | Raise `LambertGeometryError`. The 180° case is genuinely singular for the universal-variable single-rev form; the multi-rev `low/high` branches resolve it but those are M4. The 0° case implies the orbit is the same point — no transfer makes sense. |
| Newton fails to converge in 60 steps | `LambertConvergenceError` with the last `z` and residual; the caller can retry with a different `prograde` or perturbed `tof`. |
| `n_revs >= 1` requested but `max_revs >= 1` | M1 returns only the single-rev solution (list of length 1). Logged in the docstring as a known M4 stub. |
| `prograde=False` | Flips the sign convention on the transfer angle; tested via a known retrograde leg. |

All three error types subclass a common `LambertError(Exception)` defined in `lambert.py`.

#### 3.2.5 What `scipy` is for (and why it's a runtime dep regardless)

The Lambert inner loop is numpy-only — no scipy. But `scipy` is added in M1 because (a) the Stumpff implementation is easier to validate against `scipy.special` series, (b) M2's Tisserand work needs `scipy.optimize.brentq` for V∞-contour root finding and pulling that dep forward keeps M2's plan smaller. Adding `scipy>=1.13` here, used only at module test time in M1, is honest scope creep.

### 3.3 `core/kepler.py`

#### 3.3.1 Public API

```python
def propagate(
    r0: Vec3,
    v0: Vec3,
    dt: float,
    mu: float = MU_SUN_KM3_S2,
) -> tuple[Vec3, Vec3]:
    """Two-body universal-variable propagation of (r0, v0) by dt seconds.
    Returns (r, v) in the same inertial frame. dt may be positive or negative."""
```

#### 3.3.2 Algorithm: Vallado Algorithm 3.4 (universal anomaly, Stumpff)

1. `r0_n = ||r0||`, `v0_n = ||v0||`, `α = 2/r0_n − v0_n²/μ` (reciprocal semi-major axis; sign distinguishes elliptic/parabolic/hyperbolic).
2. Initial χ guess from Vallado: elliptic uses `χ₀ = sqrt(μ)·|α|·dt`; parabolic uses a special form; hyperbolic uses a log expression.
3. Newton iteration on χ with Stumpff `C(z), S(z)` where `z = χ²·α`:
   - `f(χ) = (r0·v0/sqrt(μ))·χ²·C(z) + (1 − α·r0_n)·χ³·S(z) + r0_n·χ − sqrt(μ)·dt`.
   - `f'(χ) = (r0·v0/sqrt(μ))·χ·(1 − z·S(z)) + (1 − α·r0_n)·χ²·C(z) + r0_n`.
   - Tolerance `|f| < 1e-10`, cap 50 iterations. Failure raises `KeplerConvergenceError`.
4. Compute Lagrange `f, g, ḟ, ġ` and `r = f·r0 + g·v0`, `v = ḟ·r0 + ġ·v0`.

#### 3.3.3 Stumpff functions

Implemented in a private module `_stumpff.py` shared by `lambert.py` and `kepler.py`:

```python
def stumpff_c(z: float) -> float: ...
def stumpff_s(z: float) -> float: ...
```

Series expansion for `|z| < 1e-3` (avoid catastrophic cancellation); closed-form `(1 − cos sqrt z)/z` etc. otherwise. Single source of truth so any future numerical improvement lands once.

#### 3.3.4 Edge cases

| Case | Behaviour |
|---|---|
| `dt == 0` | Return `(r0.copy(), v0.copy())` immediately (skip the iteration). |
| `dt < 0` | Supported — the universal formulation handles backward propagation natively. |
| Hyperbolic orbit | Supported, no special path; the χ iteration just lands in a different regime. |
| Parabolic orbit (α ≈ 0) | Initial guess switches to the parabolic form; tested with a fabricated `v0 = sqrt(2μ/r0)` case. |

### 3.4 API summary table

| Symbol | Module | Purpose | Notes |
|---|---|---|---|
| `Ephemeris(model="circular")` | `ephemeris` | Planet-state provider | `"astropy"` raises `NotImplementedError` (M6) |
| `Ephemeris.state(body, t_sec)` | `ephemeris` | `(r_km, v_km_s)` heliocentric inertial | shape `(3,)` float64 |
| `LambertSolution` | `lambert` | Frozen dataclass | `n_revs, branch, v1, v2` |
| `lambert(r1, r2, tof, *, mu, prograde, max_revs)` | `lambert` | UV solver, list-return | M1: single-rev only |
| `lambert_crosscheck(r1, r2, tof, ...)` | `lambert` | Mine vs izzo vs gooding | Gate-test consumer |
| `LambertError` / `LambertGeometryError` / `LambertConvergenceError` | `lambert` | Error hierarchy | |
| `propagate(r0, v0, dt, mu)` | `kepler` | UV Kepler propagation | `dt` may be negative |
| `KeplerConvergenceError` | `kepler` | Solver failure | |
| `_stumpff.stumpff_c/s(z)` | private | Shared utility | not re-exported |

---

## 4. Tests

Tests live under `tests/core/`. Fixtures in `tests/core/conftest.py` build the three canonical test legs once for reuse.

### 4.1 Test legs (the gate)

Three legs chosen to span short / medium / long ToF regimes. All use `Ephemeris(model="circular")` placed at deliberately chosen epochs so the geometry is non-degenerate (no 180° transfers).

| Leg | From | To | ToF | Regime | Why |
|---|---|---|---|---|---|
| **A — Aldrin E→M** | Earth at `t=0` | Mars at `t=146 d`, with Mars phase-offset so the transfer is well-posed | 146 d | medium | The headline anchor (spec §9). |
| **B — Earth-to-Earth short arc** | Earth at `t=0` | Earth at `t=50 d` (Earth advances ~49.3° in its orbit) | 50 d | short | Stresses small Δν, small `A`; catches sign bugs. |
| **C — Earth-to-Mars-and-back long arc** | Earth at `t=0` | Mars at `t=500 d` (Mars advances ~262°) | 500 d | long | Stresses large transfer angle and slow Newton convergence. |

For each leg, the gate test asserts `lambert_crosscheck(...)["max_diff_mps"] < 1e-3`.

### 4.2 Lambert ↔ Kepler self-consistency

For each of the three legs above: take the Lambert `(v1, v2)` solution, call `kepler.propagate(r1, v1, tof)` → `(r_end, v_end)`, assert:

- `||r_end − r2|| < 1.0` km (positional re-closure).
- `||v_end − v2|| < 1e-4` km/s (velocity re-closure, ~0.1 m/s — generous).

This catches the failure mode where the Lambert solver and `lamberthub` share a numerical bias against the universe — if Kepler integrated from the Lambert state lands on `r2`, the physics is internally consistent.

### 4.3 Ephemeris tests

- `test_circular_earth_period`: `Ephemeris.state("E", 0)` and `Ephemeris.state("E", 365.25 * 86400)` agree to < 1 km in position. Confirms mean-motion-based period closes.
- `test_circular_speeds`: `||v||` matches `sqrt(μ_sun / r)` for each planet within 1 m/s.
- `test_planar`: `r[2] == 0` and `v[2] == 0` exactly for every planet at multiple epochs.
- `test_astropy_not_implemented`: `Ephemeris(model="astropy")` raises `NotImplementedError` whose message contains `"M6"`.
- `test_unknown_body`: `Ephemeris(model="circular").state("Pluto", 0.0)` raises `KeyError`.

### 4.4 Lambert standalone tests

- `test_lambert_returns_list_singleton`: `len(lambert(...)) == 1`, `.n_revs == 0`, `.branch == "single"`.
- `test_lambert_max_revs_stub`: `max_revs=2` still returns a length-1 list (documented M4 stub); a warning isn't necessary but the docstring is checked via `pytest --doctest-modules` if cheap.
- `test_lambert_retrograde`: a known retrograde leg matches lamberthub with `prograde=False`.
- `test_lambert_zero_tof_raises`: `ValueError` on `tof <= 0`.
- `test_lambert_180_deg_raises`: a transfer set up to be exactly 180° raises `LambertGeometryError`. (Construct with `r2 = -k * r1` for some `k > 0`.)
- `test_lambert_solution_dataclass_frozen`: assigning to a field raises `FrozenInstanceError`.

### 4.5 Kepler standalone tests

- `test_kepler_zero_dt_identity`: `propagate(r0, v0, 0.0)` returns `(r0, v0)`.
- `test_kepler_circular_period`: propagating Earth's circular state by one circular period returns the start state within tolerance.
- `test_kepler_reversibility`: `propagate(*propagate(r0, v0, +dt), -dt)` returns `(r0, v0)` within tolerance.
- `test_kepler_energy_conservation`: specific orbital energy `v²/2 − μ/r` constant to < 1e-8 relative across propagation.
- `test_kepler_hyperbolic`: a fabricated hyperbolic state (`v0 > sqrt(2μ/r0)`) propagates forward then back to identity.

### 4.6 Tolerance summary (documented in module docstrings)

| Layer | Quantity | Tolerance |
|---|---|---|
| Stumpff series cutoff | `|z|` | `1e-3` |
| Kepler χ convergence | `|f(χ)|` | `1e-10` |
| Lambert z convergence | `|t(z) − tof| / tof` OR `|Δz|` | `1e-12` |
| Lambert iteration cap | iterations | 60 |
| Kepler iteration cap | iterations | 50 |
| Test: Lambert vs lamberthub | max Δ\|v\| | `< 1e-3` m/s (gate) |
| Test: Lambert→Kepler re-closure | `||r_end − r2||` | `< 1.0` km |
| Test: ephemeris circular period | position drift | `< 1.0` km |

---

## 5. Risks

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Lambert Newton stalls on long-arc leg C | medium | gate fails | Use Vallado's bisection fallback when `y(z) < 0`; cap iterations at 60 with informative error; test leg C explicitly. |
| Stumpff cancellation near `z = 0` | low | silent accuracy loss | Series expansion below `|z| < 1e-3`; unit tests against `scipy.special` reference values at boundary. |
| `lamberthub` API drift between versions | low | dev-only friction | Pin `lamberthub>=1.0,<2.0` in dev-extras; import via thin wrapper so a future API change touches one place. |
| `mypy --strict` friction on numpy 2.0 generics | medium | dev velocity | Use `numpy.typing.NDArray[np.float64]` aliased to `Vec3`; avoid `np.ndarray[Any, Any]`. If `mypy` chokes on a specific call, the documented escape hatch is a targeted `# type: ignore[<code>]` with a comment — *not* turning off strict mode. |
| 180° transfer arises in legitimate later use | medium | confusing UX | Raise `LambertGeometryError` with a message that names the M4 multi-rev workaround; M3's Aldrin leg is not 180° so M1 is fine. |
| Ephemeris epoch convention bites M3 | low | replanning M3 | `t=0` ⇒ all planets at `θ=0`; documented in module docstring. M3 chooses its own epoch shift; it doesn't depend on M1 phases being calibrated to a real date. |
| scipy added as runtime dep but unused in M1 module code | low | wasted install weight | Acceptable — M2 needs it anyway and pulling it forward is honest. README mentions it. |
| Choosing the "single" branch label for `n_revs=0` differs from any future convention | low | minor refactor | `branch` is a string; the dataclass is frozen but the M4 author can decide if `"single"` stays or becomes `"prograde"/"retrograde"`. Documented in the dataclass docstring. |

---

## 6. Dependency additions

Per the M0 plan §7 dependency-policy ("add a runtime dep only when the first module that needs it is being written"), M1 adds:

| Dep | Version | Where | Why |
|---|---|---|---|
| `scipy` | `>=1.13` | `[project.dependencies]` | Stumpff validation reference; M2's `brentq` for V∞ contours. Not in Lambert's inner loop. |
| `lamberthub` | `>=1.0,<2.0` | `[project.optional-dependencies] dev` | Cross-check oracle only. Never imported by production code. |

`uv add scipy` and `uv add --dev lamberthub` regenerate `uv.lock`; the lockfile is committed. CI's `uv sync --frozen --all-extras` picks both up.

No removals. M0's `numpy>=2.0` stays.

---

## 7. Order of work

`todo.md` mirrors this with checkboxes. High-level phases:

1. **Plan/docs first** — this file + `todo.md` committed before any code. (You are here.)
2. **Dependencies** — `uv add scipy` (runtime), `uv add --dev lamberthub` (dev); confirm `uv sync` and CI's `--frozen` mode still pass.
3. **Source skeleton** — create `core/ephemeris.py`, `core/lambert.py`, `core/kepler.py`, `core/_stumpff.py` with module docstrings and signatures only (stubs raise `NotImplementedError`). Verify `mypy --strict` accepts the type signatures before any logic.
4. **Stumpff** — implement and unit-test `_stumpff` first; everything else depends on it.
5. **Kepler** — implement `propagate`; test in isolation (period closure, reversibility, energy). Doing Kepler before Lambert means Lambert's self-consistency check has a trusted propagator from day one.
6. **Ephemeris (circular)** — implement and test. Used as the source of `r1, r2` for Lambert legs.
7. **Lambert** — implement single-rev solver; unit-test standalone behaviours (180° error, retrograde, zero-tof) first; *then* the gate (cross-check + Lambert-Kepler self-consistency) on the three test legs.
8. **Tests/CI** — confirm `uv run pytest`, `ruff check`, `ruff format --check`, `mypy src tests` all clean locally, then push and watch CI.
9. **Closeout** — update `docs/overview.md` milestone table (M1 → completed; M2 → planned); append a `## Hand-off to M2` section to this milestone's `todo.md` with anything learned.

---

## 8. Exit checklist (the gate, restated)

Before declaring M1 done:

- [ ] `uv run pytest` green locally, including:
  - [ ] The three Lambert cross-check tests (legs A, B, C) all assert `max_diff_mps < 1e-3`.
  - [ ] The three Lambert→Kepler self-consistency tests all assert positional re-closure < 1.0 km.
  - [ ] All standalone ephemeris, lambert, kepler tests pass.
- [ ] `uv run ruff check .` clean.
- [ ] `uv run ruff format --check .` clean.
- [ ] `uv run mypy src tests` clean under `strict = true`.
- [ ] CI green on push.
- [ ] `pyproject.toml` lists `scipy` (runtime) and `lamberthub` (dev); `uv.lock` regenerated and committed.
- [ ] `docs/overview.md` §4 milestone table: M1 status `planned` → `completed`; M2 row `not yet planned` → `planned`.
- [ ] `## Hand-off to M2` appended to `phases/m1-core-mechanics/todo.md` covering:
  - Confirmed module APIs (any signature drift from this plan).
  - Numerical surprises (anywhere a tolerance had to be loosened).
  - Anything that should inform M2's flyby/Tisserand work — e.g. observed Lambert performance at long ToF, since M2's enumerator will call it heavily.

(Writing the M2 plan doc is the first task of M2, not an M1 exit criterion.)
