# N-body harness — restricted-n-body propagator, SILVER rungs, Jones shooter, GMAT V4 stamp

> **For agentic workers:** REQUIRED SUB-SKILL: `superpowers:subagent-driven-development`
> or `superpowers:executing-plans`. Checkbox steps; strict TDD (write failing
> test → run **red** → minimal impl → run **green** → commit). Work on `main` —
> **do NOT branch** (project rule). **Do NOT commit while this plan is authored
> under the docs-only mandate**; the commit messages below are the messages to
> *use* when the implementation phase runs. uv-managed venv (**no pip**).
> Lint/type gate before **every** commit:
> `uv run ruff check .` · `uv run ruff format --check .` · `uv run mypy src tests`.
> Fast suite: `uv run pytest -m "not slow"`; every propagation/closure test is
> `slow` (DE440-backed, see `docs/notes/2026-06-06-performance-profile.md`).
>
> This plan is the task-level expansion of the **approved** design
> `docs/superpowers/specs/2026-06-06-nbody-harness-design.md` (read it **including
> its Approval section** — all seven open questions resolved there are binding).
> Where this plan and the design could diverge, **the design + its Approval win.**
> The approved resolutions this plan implements: (Q1) REBOUND-first, finite-diff
> STMs, Tudat only on demonstrated need; (Q2) hand-rolled ephemeris force callback
> against the shared DE440 BSP — **no `assist`**; (Q3) node-impulse correction-ΔV,
> *comparable* (documented mapping) but not numerically identical to
> `maintain.py`'s per-synodic value; (Q4) the §5.3 sensitivity test governs body
> inclusion unless AAS 17-577 states its set; (Q5) V4 tolerance **5%**, declared
> up front, never back-fit; (Q6) rungs run on shared DE440 as **cross-check only**,
> never recorded V4; (Q7) GMAT is a manual, out-of-CI lane.

---

## Goal

Build **one** restricted-n-body harness, sized at once for three consumers
(design §0), behind a single propagator interface so REBOUND↔Tudat is a backend
swap not a rewrite (mirroring how `core/ephemeris.py` hides `circular`/`astropy`/
`inclined` behind one `state()`):

1. **SILVER validation rungs 2–3** (#131) — propagate the two USER-HELD SILVER
   candidates (`data/review_queue.jsonl`, the Forge Phase-4 E-M-E-E pair) through
   full n-body, measure terminal closure degradation + the **node-impulse
   correction ΔV**, and **record the result to the review-queue audit trail —
   never auto-promote** (`review_queue.py:62` `is_catalogue_source()` returns
   `False` by contract; the human decides).
2. **The Jones shooter** (#133, deferred Phase 3 of
   `2026-06-06-jones-family-corrector-variant-design.md`) — multiple-shooting
   differential correction in n-body over encounter nodes, seeded from the
   near-miss conic survey (#133 Phase 3a), gated against the **sourced** Jones
   VEM multisets (fold-corrected comparison, unchanged criteria).
3. **Aldrin V4** (#134, spec §14 V4, `docs/spec.md:398`) — GMAT (independent
   codebase + ephemeris) reproduces the Aldrin powered-periodic solution and our
   computed **2.9138 km/s** per-synodic maintenance ΔV within the **declared 5%**
   tolerance, as a manual out-of-CI stamp.

The genuinely new content is the **substrate + its §0 self-validation** (Phase A):
a propagator interface, a hand-rolled rails force callback, and the time/frame
conversion that **reuses #129's kernel machinery** (`verify/spice_kernels.py`)
rather than reinventing it. Everything downstream (rungs, shooter, stamp)
consumes that substrate.

**Honesty boundary (binding, carried from design §5 + the live #135 verdict).**
The #135 like-for-like diagnostic
(`docs/notes/2026-06-06-russell12-likeforlike.md`) ran coplanar-vs-coplanar on
known-solvable Russell/McConaghy instances and found **0 CLOSE-AND-MATCH**: the
corrector closes geometrically but lands **off-anchor** (our V∞ ≈ 9–28 km/s vs
sourced 3–10 km/s) even with the model mismatch removed — a **seeding/basin**
problem, not pure model mismatch. The two SILVER candidates themselves floor at
**Earth V∞ 9.62–9.75 / Mars V∞ 12.06–13.01 km/s**
(`docs/notes/2026-06-06-silver-candidates-russell-diligence.md`), i.e. the same
high-V∞ basin. **This plan does not assume the Jones shooter reaches the 2.5–3.9
km/s Jones regime.** Phase C carries an explicit **#135-verdict checkpoint**
(Task C.0) and a STOP/report branch (Task C.5): the harness + rungs + shooter
infrastructure are the deliverable; reaching the Jones multiset is recorded as
open research if the basin problem is not overcome by the near-miss seeding.

---

## Architecture

### The propagator interface (design §1 "single interface", §3 backend swap)

A `Propagator` protocol over the **TDB-J2000 seconds axis + heliocentric
J2000-ecliptic frame** that `core/ephemeris.py` already defines (`ephemeris.py:47`
`_J2000_EPOCH`, `:53` `_J2000_OBLIQUITY_RAD = radians(23.4392911)`, `:262-263`
the `Time(..., scale="tdb")` build). One method:

```python
propagate(r0_km, v0_km_s, t0_sec, t1_sec, *, bodies, accuracy) -> NBodyArc
```

returning a frozen `NBodyArc` (final state + energy/anchor diagnostics +
integrator settings + converged/diverged). Backends: `RestrictedNBody("rebound")`
(Phase A baseline) and a deferred `RestrictedNBody("tudat")` slot (Phase C, only
on demonstrated STM need). The **glue owns the §0 conversions** for REBOUND
(design §1 "REBOUND alone has no SPICE time/frame machinery").

### The rails force callback (design §2, Q2 hand-rolled)

**Restricted n-body, heliocentric, ephemeris-driven point-mass perturbers.**
Spacecraft is a massless test particle; planets are **on rails** from the shared
DE440 BSP (`spice_kernels.astropy_de440_bsp_path()`), no back-reaction (exact,
since the spacecraft is massless). The callback at time `t_sec`:

- reads each perturber's heliocentric J2000-ecliptic state via
  `Ephemeris("astropy").state(body, t_sec)` (the **same reader** the rungs
  cross-check against — that is the point of the shared-kernel rung, design §4);
- accumulates `a = -μ_sun r̂/r² + Σ_p μ_p (Δ_p/|Δ_p|³ − r_p/|r_p|³)` (the
  standard third-body indirect term), `μ` from `core/constants` (`MU_SUN_KM3_S2`,
  `PLANETS[code].mu_km3_s2`, `constants.py:85,164`).
- **Velocity-independent external force** (no drag/SRP, design §2) — IAS15 with a
  velocity-independent additional force is well-behaved; this is asserted in the
  energy-conservation gate (Task A.4), not assumed (risk register).

Per-consumer body sets (design §2): **rungs** Sun+E+M+J; **Jones** Sun+V+E+M(+J
per §5.3 sensitivity); **V4** the external tool's standard set. SRP / GR /
non-spherical / moons are **OUT** with the design-§2 justifications.

### The multiple-shooting corrector (design §3)

Nodes = encounters (the patch points, matching `correct.py`'s `b{i}_in/_out`
vocabulary, `correct.py:86` `_vinf_nodes`). Free vars = node states + node epochs
+ ToFs, with the slack-leg/period pin carried from `correct.py`'s
`_reconstruct_tofs` (`correct.py:77-86`). Defects = full-state continuity in real
dynamics (strictly stronger than `correct.py`'s magnitude-continuity residual).
Solver: `least_squares(method="lm")` (mirrors `correct.py`), Jacobians by
finite-difference over REBOUND (Q1 baseline). Divergence is a first-class
non-converged record, not an exception (mirrors `correct.py:354-378`'s honest
non-converged `BallisticClosureResult`).

### The correction-ΔV convention (design §3, Q3)

**Node-impulse**: the sum of the impulsive ΔV at the encounter nodes that the
flybys could *not* absorb gravitationally and that a real mission would burn to
restore n-body periodicity from the conic seed. **Comparable** to `maintain.py`'s
per-synodic maintenance ΔV (`maintain.py:24,41` "the per-synodic maintenance ΔV";
the 2.9138 km/s the V4 consumer reproduces) — same physical meaning, **documented
mapping**, not forced numerically identical (Q3). Sign- and node-explicit.

### Independence accounting (design §4, Q6)

The rungs run on the **shared** astropy DE440 BSP: independent *integrator*
(REBOUND, not `core/kepler.py`) over *identical data* — a disagreement is a
reader/frame/time bug, which is exactly what the rung catches. **The rung result
is NEVER recorded as a V4 pass** (Q6); it is a cross-check rung feeding the
review-queue audit trail (`review_queue.py`) / the gauntlet
(`verify/gauntlet.py`). Only the GMAT lane (independent codebase + ephemeris)
earns V4 (`docs/spec.md:402`).

---

## Tech stack

Python 3.11, numpy, `scipy.optimize.least_squares`, **rebound** (new optional
dep), pytest. uv-managed venv. New production package:
`src/cyclerfinder/nbody/` — `propagator.py` (the protocol + `NBodyArc`),
`forces.py` (the rails callback), `convert.py` (the §0 time/frame glue, reusing
`verify/spice_kernels`), `shooter.py` (multiple shooting), `correction_dv.py`
(the node-impulse metric). New tests under `tests/nbody/`. New manual lane:
`scripts/gmat_v4_aldrin.py` (generator) + `scripts/parse_gmat_report.py` (parser),
out of CI. The rung consumer touches `data/review_queue.py` (audit-trail append
only). Dependencies of `nbody/` are limited to `core/ephemeris`, `core/constants`,
`core/kepler` (for the golden gate), `verify/spice_kernels`, numpy, scipy, rebound.

### Dependency convention (design Q2, pyproject)

The `validation` extra (`pyproject.toml:36-38`, currently `spiceypy>=6.0`) is the
established home for independent-toolchain validation that "reads the SAME DE440
kernel astropy already cached … to catch a frame/time/reader bug". REBOUND on the
shared BSP is **exactly this class** (design §4: cross-check rung, shared
ephemeris). Add `rebound` to the existing `validation` extra rather than minting a
new `nbody` extra — Task A.0 below. Tudat, if ever needed (Q1), is a separate
decision with a separate extra (its conda/binary weight, design §1, is why it is
NOT pulled in now).

---

## Design references (verified against live code, 2026-06-06)

- §0 epoch — `_J2000_EPOCH` = 2000-01-01T12:00:00, `Time(..., scale="tdb")`
  (`ephemeris.py:47,262-263`). **Verified.** TDB-seconds-since-J2000 axis.
- §0 frame — `_J2000_OBLIQUITY_RAD = radians(23.4392911)`, ICRS→ecliptic rotation
  about +x by −obliquity (`ephemeris.py:53,241-242,271-275`). **Verified.**
- §0 same-kernel — `astropy_de440_bsp_path()` + `ensure_leapseconds_kernel()`
  (`spice_kernels.py:44,71`), `naif0012.tls` LSK (`spice_kernels.py:40`).
  **Verified** — the harness MUST reuse these, not re-fetch.
- §3 seed object — `BallisticClosureResult` fields `t0_sec, tof_days,
  vinf_per_encounter_kms, max_residual_kms, converged, bend_feasible,
  vinf_cap_ok`, `constraints_satisfied` property (`correct.py:39-50`).
  **Verified — DRIFT:** the design §3 sketch refers to
  `result.best_cycler.encounters[].vinf_in/vinf_out` (per-node *vectors*); the
  live `BallisticClosureResult` carries **scalar** `vinf_per_encounter_kms` and
  the per-node *vectors* live only inside `_vinf_nodes`'s `b{i}_in/_out` dict
  (`correct.py:86`). The shooter seed (Task C.1) must therefore read node vectors
  from `_vinf_nodes`, not from a `best_cycler` attribute that does not exist on
  the result. Flagged.
- SILVER candidates — `review_queue.py` `ReviewQueueEntry` carries
  `vinf_per_encounter_kms, tof_days, sequence, period_k, bend_feasible,
  verdict_audit, panel` (`review_queue.py:44-59`); only SILVER/GOLD queueable
  (`review_queue.py:28,82`); the two held E-M-E-E candidates float at E∞
  9.62–9.75 / M∞ 12.06–13.01 km/s
  (`docs/notes/2026-06-06-silver-candidates-russell-diligence.md`). **Verified.**
- §3 maintenance convention — `maintain.py:24,41` per-synodic maintenance ΔV;
  `optimise_aldrin_maintenance_dv` (`maintain.py:717`) →
  `solve_powered_periodic_cycler` (`bvp.py:92,149`). 2.9138 km/s is **our
  computed** Aldrin value (design §0/§5.5), recorded in
  `docs/notes/2026-06-06-performance-profile.md`; it is **not** a published Aldrin
  number — golden caveat applies (in-pipeline reproduction = consistency check;
  GMAT = genuine external validation). **Verified.**
- §4 gauntlet tiers — `VerdictTier{GOLD,SILVER,BRONZE,REJECTED}`
  (`gauntlet.py:91-94`); SILVER = "machine-confirmed but unsourced, capped
  pending human review, never auto-promoted to GOLD". **Verified.**
- #135 verdict — `docs/notes/2026-06-06-russell12-likeforlike.md`: 0
  CLOSE-AND-MATCH coplanar-vs-coplanar ⇒ seeding/basin, not model mismatch.
  **Verified verbatim** (the per-row table; ach V∞ 9–28 vs src 3–10).
- two-body golden — `core/kepler.py::propagate` (`kepler.py:141`),
  `KeplerConvergenceError` (`kepler.py:45`). **Verified.**
- pyproject — `validation` extra exists (`pyproject.toml:36`, `spiceypy>=6.0`).
  **Verified.**

---

## Phasing (independently shippable)

| Phase | Theme | Tasks | rebound dep | shared-DE440 | external tool |
|---|---|---|---|---|---|
| **A** | Substrate + §0 conversions + the golden gates | A.0–A.6 (7) | yes | yes | no |
| **B** | SILVER rungs (ingest queue, propagate, correction ΔV, record) | B.0–B.4 (5) | yes | yes | no |
| **C** | Jones shooter (multiple shooting, #135 checkpoint, near-miss seeds) | C.0–C.6 (7) | yes | yes | Tudat (conditional) |
| **D** | GMAT V4 stamp (manual lane, 5% tol, run-book) | D.0–D.3 (4) | no | no | GMAT (manual) |

Phase A is foundational and **must land first** (the §0 proofs gate everything,
design §6). Phase B ships the immediately-useful rung verdicts on the two held
candidates. Phase C is the hardest dynamics and is gated on the #135 verdict +
the #133 near-miss survey. Phase D is a manual out-of-CI stamp. **Total: 23
tasks.**

---

## Phase A — substrate + the §0 conversions + the golden gates

The §0 time/frame proofs come **first** (design Approval: "with §0 time/frame
conversion proofs FIRST"). Build the conversion glue and the propagator, then
prove them with the five golden gates *before* any science number is trusted.

### Task A.0 — add `rebound` to the `validation` extra + `nbody` package skeleton

**Files:** `pyproject.toml`; create `src/cyclerfinder/nbody/__init__.py`; test
`tests/nbody/test_package_imports.py`.

#### Failing test — `tests/nbody/test_package_imports.py`

```python
"""N-body harness Phase A: package + optional-dep wiring (plan Phase A)."""
from __future__ import annotations

import importlib

import pytest


def test_nbody_package_imports() -> None:
    mod = importlib.import_module("cyclerfinder.nbody")
    assert mod is not None


def test_rebound_is_in_validation_extra() -> None:
    import tomllib
    from pathlib import Path

    root = Path(__file__).resolve().parents[2]
    data = tomllib.loads((root / "pyproject.toml").read_text())
    extras = data["project"]["optional-dependencies"]
    assert any(d.startswith("rebound") for d in extras["validation"]), (
        "rebound must join the existing 'validation' extra (design Q2: shared-"
        "DE440 cross-check class), not a new extra"
    )


def test_rebound_skips_cleanly_when_absent() -> None:
    """The fast suite must not hard-require rebound (mirror spiceypy skip)."""
    rebound = pytest.importorskip("rebound")
    assert rebound is not None
```

Run: `uv run pytest tests/nbody/test_package_imports.py -q` → **red** (no package).

#### Minimal impl

Create `src/cyclerfinder/nbody/__init__.py` with a module docstring stating the
package is **validation infrastructure only** (seeds-not-tracks intact: nothing
here is consumed by construct/score/verify of catalogue rows — mirror
`spice_kernels.py:5-8`). Add `"rebound>=4.0"` to the `validation` extra in
`pyproject.toml:36-38`. Run `uv sync --extra validation`.

Run → **green**. Lint/type. Commit:

```
nbody: package skeleton + rebound in the validation extra (n-body harness Phase A)
```

### Task A.1 — §0 time conversion: TDB-J2000-seconds ↔ SPICE ET (reuse #129)

**Files:** create `src/cyclerfinder/nbody/convert.py`; test
`tests/nbody/test_convert_time.py` (slow — needs the LSK).

The harness `t_sec` axis is **TDB seconds since J2000** (`ephemeris.py:262-263`).
SPICE-native tools need ET; the TDB↔ET identity (within ~ms) is the safe path,
TT/UTC is the trap (design §0). Provide `t_sec_to_et(t_sec)` that loads #129's LSK
(`spice_kernels.ensure_leapseconds_kernel()`) and returns ET = J2000-relative
TDB seconds (the identity, asserted against spiceypy `str2et` of the J2000 epoch +
offset). **Reuse, do not reinvent** (`spice_kernels.py:71`).

#### Failing test — `tests/nbody/test_convert_time.py`

```python
"""N-body Phase A: t_sec (TDB-J2000) -> ET reuses #129's LSK (plan Phase A; design §0).

GOLDEN DISCIPLINE: EXPECTED side = SPICE's own str2et over the J2000 epoch (the
LSK-defined time scale), not a value our code computed. This catches the TDB vs
UTC/TT trap (design §0: the "~64.184 s with drift" failure).
"""
from __future__ import annotations

import pytest

spiceypy = pytest.importorskip("spiceypy")

from cyclerfinder.nbody.convert import t_sec_to_et
from cyclerfinder.verify.spice_kernels import ensure_leapseconds_kernel


@pytest.mark.slow
def test_t_sec_zero_is_j2000_et() -> None:
    spiceypy.furnsh(ensure_leapseconds_kernel())
    try:
        et_j2000 = spiceypy.str2et("2000-01-01T12:00:00 TDB")
        assert t_sec_to_et(0.0) == pytest.approx(et_j2000, abs=1e-3)
    finally:
        spiceypy.kclear()


@pytest.mark.slow
def test_t_sec_is_tdb_not_utc() -> None:
    """A one-year offset stays a clean TDB second count (no leap-second jump)."""
    spiceypy.furnsh(ensure_leapseconds_kernel())
    try:
        one_year = 365.25 * 86400.0
        et = spiceypy.str2et("2000-01-01T12:00:00 TDB") + one_year
        assert t_sec_to_et(one_year) == pytest.approx(et, abs=1e-3)
    finally:
        spiceypy.kclear()
```

Run → **red** → impl `t_sec_to_et` (ET = `str2et("J2000 TDB") + t_sec`, the
identity; document that we never route through UTC/TT). → **green** (`-m slow`).
Commit:

```
nbody/convert: TDB-J2000-seconds -> ET reusing #129's LSK (n-body harness Phase A)
```

### Task A.2 — §0 frame conversion: ICRS/barycentric ↔ heliocentric J2000-ecliptic

**Files:** `nbody/convert.py`; test `tests/nbody/test_convert_frame.py`.

Any n-body tool returns barycentric/ICRS-equatorial by default; the harness MUST
(a) Sun-subtract to heliocentric and (b) apply the **same obliquity rotation** as
`ephemeris.py:241-242,271-275` (`R_x(−23.4392911°)`) before any vector is compared
to a `BallisticClosureResult` V∞ or a catalogue value (design §0). Provide
`icrs_eq_to_ecliptic(vec3)` / `ecliptic_to_icrs_eq(vec3)` using the **same**
`_J2000_OBLIQUITY_RAD` constant (import it, do not re-literal it).

#### Failing test — `tests/nbody/test_convert_frame.py`

```python
"""N-body Phase A: ICRS-equatorial <-> J2000-ecliptic reuses ephemeris.py's obliquity.

NON-GOLDEN cross-implementation check: both sides are OUR rotation; this asserts
the harness uses the SAME obliquity constant as core/ephemeris (design §0: a wrong
obliquity reads as a fake out-of-plane V_inf component).
"""
from __future__ import annotations

import numpy as np

from cyclerfinder.nbody.convert import ecliptic_to_icrs_eq, icrs_eq_to_ecliptic


def test_roundtrip_identity() -> None:
    v = np.array([1.0, 2.0, 3.0])
    back = ecliptic_to_icrs_eq(icrs_eq_to_ecliptic(v))
    assert np.allclose(back, v, atol=1e-12)


def test_uses_ephemeris_obliquity_constant() -> None:
    from cyclerfinder.core.ephemeris import _J2000_OBLIQUITY_RAD

    # +z equatorial maps to (0, +sin eps, +cos eps) in ecliptic about +x.
    z_eq = np.array([0.0, 0.0, 1.0])
    got = icrs_eq_to_ecliptic(z_eq)
    eps = _J2000_OBLIQUITY_RAD
    assert np.allclose(got, [0.0, np.sin(eps), np.cos(eps)], atol=1e-12)
```

Run → **red** → impl both with `R_x(-eps)` / `R_x(+eps)`, `eps =
_J2000_OBLIQUITY_RAD` imported from `core.ephemeris`. → **green**. Commit:

```
nbody/convert: ICRS-eq <-> J2000-ecliptic via the shared obliquity (n-body harness Phase A)
```

### Task A.3 — the rails force callback + the `Propagator` interface

**Files:** create `src/cyclerfinder/nbody/forces.py`,
`src/cyclerfinder/nbody/propagator.py`; test `tests/nbody/test_propagator_api.py`.

`propagator.py`: a `Propagator` protocol with `propagate(...)` and a frozen
`NBodyArc` (`r_km, v_km_s, t1_sec, energy_rel_drift, anchor_err_km,
integrator_accuracy, bodies, converged`). `RestrictedNBody("rebound")`
implements it. `forces.py`: `rails_acceleration(r_km, t_sec, bodies, ephem)`
returning the Sun + third-body indirect acceleration, `μ` from `core/constants`.
This task wires the API + the REBOUND `additional_forces` callback shell; the
physics is proven in A.4–A.6.

#### Failing test — `tests/nbody/test_propagator_api.py`

```python
"""N-body Phase A: propagator interface + rails-force shape (plan Phase A)."""
from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.constants import MU_SUN_KM3_S2
from cyclerfinder.nbody.forces import rails_acceleration


def test_sun_only_acceleration_is_central_inverse_square() -> None:
    # No perturbers -> pure -mu_sun r_hat / r^2.
    r = np.array([1.5e8, 0.0, 0.0])  # ~1 AU on +x
    a = rails_acceleration(r, t_sec=0.0, bodies=(), ephem=None)
    expected = -MU_SUN_KM3_S2 / (1.5e8**2)
    assert a[0] == pytest.approx(expected, rel=1e-12)
    assert a[1] == 0.0 and a[2] == 0.0


def test_propagator_backend_selectable() -> None:
    pytest.importorskip("rebound")
    from cyclerfinder.nbody.propagator import RestrictedNBody

    prop = RestrictedNBody("rebound")
    assert prop.backend == "rebound"
    with pytest.raises(ValueError, match="backend"):
        RestrictedNBody("tudat")  # deferred slot: not wired yet (design Q1)
```

Run → **red** → impl `rails_acceleration` (Sun term first; perturber loop reads
`ephem.state(body, t_sec)` and adds the indirect term) + `RestrictedNBody` with a
`backend` attr that accepts only `"rebound"` for now (Tudat raises `ValueError`,
deferred per Q1). → **green**. Commit:

```
nbody: Propagator interface + rails third-body force callback (n-body harness Phase A)
```

### Task A.4 — GOLDEN GATE 1: two-body limit reduction (the keystone)

**Files:** test `tests/nbody/test_golden_twobody.py` (slow).

> **GOLDEN GATE 1 (verbatim, the keystone — design §5.1).** Configure the
> n-body harness with **Sun only** (no perturbers) and propagate a state; it MUST
> equal `core/kepler.py::propagate` to tight tolerance (both are then the same
> two-body problem). This is a *cross-implementation* golden test — neither side
> is the "sourced" side because both must agree with the analytic Kepler solution
> — and it directly catches a frame/time/units slip in the §0 glue.

Tolerance: position agreement `< 1 km` over a multi-month arc (tighten until it
stops moving, Task A.6). Source-free: both sides are independent solvers of the
identical two-body problem.

#### Failing test — `tests/nbody/test_golden_twobody.py`

```python
"""N-body GOLDEN GATE 1: Sun-only n-body == core/kepler.propagate (design §5.1).

SOURCE-FREE cross-implementation golden: REBOUND/IAS15 and core/kepler must both
reduce to the SAME analytic two-body solution. A disagreement is a §0 frame/time/
units slip in the glue, not a physics choice. Neither side is EXPECTED; the
analytic Kepler problem is.
"""
from __future__ import annotations

import numpy as np
import pytest

rebound = pytest.importorskip("rebound")

from cyclerfinder.core import kepler
from cyclerfinder.nbody.propagator import RestrictedNBody


@pytest.mark.slow
def test_sun_only_matches_kepler_propagate() -> None:
    r0 = np.array([1.496e8, 0.0, 0.0])      # ~1 AU
    v0 = np.array([0.0, 29.78, 0.0])        # ~circular
    dt = 120.0 * 86400.0                     # 120 days
    kep_r, kep_v = kepler.propagate(r0, v0, dt)
    arc = RestrictedNBody("rebound").propagate(
        r0, v0, t0_sec=0.0, t1_sec=dt, bodies=(), accuracy=1e-12
    )
    assert np.linalg.norm(arc.r_km - kep_r) < 1.0   # < 1 km
    assert np.linalg.norm(arc.v_km_s - kep_v) < 1e-5
```

Run → **red** → flesh out `RestrictedNBody.propagate` (build a REBOUND sim with a
central Sun mass + one massless particle; for `bodies=()` no `additional_forces`
needed; IAS15; convert back through §0 glue). → **green** (`-m slow`). Commit:

```
nbody: GOLDEN GATE 1 — Sun-only n-body reduces to core/kepler (n-body harness Phase A)
```

### Task A.5 — GOLDEN GATE 2 (energy/integral conservation) + GATE 3 (Earth-anchor)

**Files:** test `tests/nbody/test_golden_conservation.py` (slow),
`tests/nbody/test_golden_anchor.py` (slow).

> **GOLDEN GATE 2 (verbatim — design §5.2).** With Sun-only, the Keplerian energy
> and angular momentum are conserved; the harness's relative energy drift over a
> full propagation MUST stay below an integrator-accuracy floor (IAS15 should give
> ~machine-precision energy conservation on the two-body sub-problem — a direct
> integrator-health check). **With perturbers, there is no exact conserved
> scalar** (the planets-on-rails system is time-dependent and non-autonomous, so a
> true Jacobi integral does NOT exist); the honest invariant is the **Sun-only
> energy floor** (assert ~machine precision) plus a **bounded** total-energy
> bookkeeping with perturbers ON (assert it stays below the integrator floor scaled
> by the perturber work, reported not asserted-tight). Design the right invariant
> honestly: do not claim a conserved Jacobi constant the rails model does not have.

> **GOLDEN GATE 3 (verbatim — design §5.4, the Earth-anchor).** Propagate a
> massless particle co-located with a planet (Earth) — i.e. assert the harness's
> planet-state ingestion (after the §0 conversions) reproduces
> `Ephemeris("astropy").state(body, t_sec)` to numerical precision at sampled
> epochs; the sourced side is DE440 itself. This proves the §0 time/frame
> conversion is correct *before* any spacecraft number is trusted.

#### Failing tests

```python
# tests/nbody/test_golden_conservation.py
"""N-body GOLDEN GATE 2: energy conservation (design §5.2).

Sun-only: relative energy drift ~ machine precision (integrator health).
Perturbers ON: NO exact Jacobi constant exists (rails = time-dependent); we
assert a BOUNDED energy budget, reported, not a fake conserved scalar.
"""
from __future__ import annotations

import numpy as np
import pytest

rebound = pytest.importorskip("rebound")

from cyclerfinder.nbody.propagator import RestrictedNBody


@pytest.mark.slow
def test_sun_only_energy_drift_below_floor() -> None:
    r0 = np.array([1.496e8, 0.0, 0.0])
    v0 = np.array([0.0, 29.78, 0.0])
    arc = RestrictedNBody("rebound").propagate(
        r0, v0, t0_sec=0.0, t1_sec=365.25 * 86400.0, bodies=(), accuracy=1e-12
    )
    assert abs(arc.energy_rel_drift) < 1e-10   # IAS15 two-body ~ machine precision
```

```python
# tests/nbody/test_golden_anchor.py
"""N-body GOLDEN GATE 3: planet-state ingestion == Ephemeris('astropy') (design §5.4).

GOLDEN: EXPECTED = DE440 itself (Ephemeris('astropy').state). This is the §0
anchor: it proves the time/frame conversion before any spacecraft number is
trusted. The sourced side is the ephemeris.
"""
from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.nbody.forces import ingest_planet_state


@pytest.mark.slow
@pytest.mark.parametrize("body", ["E", "M", "J"])
@pytest.mark.parametrize("t_sec", [0.0, 365.25 * 86400.0, 10 * 365.25 * 86400.0])
def test_ingested_planet_state_matches_ephemeris(body: str, t_sec: float) -> None:
    ephem = Ephemeris("astropy")
    r_ref, v_ref = ephem.state(body, t_sec)
    r_got, v_got = ingest_planet_state(body, t_sec, ephem)
    assert np.allclose(r_got, r_ref, atol=1e-6)   # numerical precision
    assert np.allclose(v_got, v_ref, atol=1e-9)
```

Run → **red** → impl `arc.energy_rel_drift` (REBOUND `sim.energy()` deltas) and a
thin `ingest_planet_state` (the §0-converted reader the callback uses; for the
shared-DE440 path it is the identity over `ephem.state`, which is *why* the anchor
passes — the rung's independence is the *integrator*, not the reader, design §4).
→ **green** (`-m slow`). Commit:

```
nbody: GOLDEN GATES 2-3 — energy floor + DE440 ingestion anchor (n-body harness Phase A)
```

### Task A.6 — GOLDEN GATE 4 (timestep/tolerance convergence) + Phase-A gate

**Files:** test `tests/nbody/test_golden_convergence.py` (slow); lint/type.

> **GOLDEN GATE 4 (verbatim — design §5.3).** Tighten the IAS15 accuracy
> parameter until the final state and the reported correction ΔV stop moving below
> the consumer's tolerance; report the converged setting per consumer. This same
> sweep **doubles as the §2 body-inclusion sensitivity test** (re-run with/without
> Jupiter; if the metric moves < tolerance, exclusion is justified — evidence, not
> assertion).

The convergence test asserts that halving the accuracy parameter changes the
final state by less than the looser run's own tolerance band (Cauchy-style
convergence), and **records** the converged accuracy per consumer. The
Jupiter-sensitivity arm is wired here but its *verdict* (include/exclude J for
each consumer) is recorded in Phase B/C against the real candidate baselines (the
sweep needs a real trajectory).

```python
"""N-body GOLDEN GATE 4: IAS15 accuracy convergence (design §5.3).

Asserts the final state is converged wrt the integrator accuracy parameter
(Cauchy: tighter accuracy moves the answer less than the consumer tolerance).
Doubles as the body-inclusion sensitivity harness (Jupiter on/off), whose verdict
is recorded against the real candidate baselines in Phase B/C.
"""
from __future__ import annotations

import numpy as np
import pytest

rebound = pytest.importorskip("rebound")

from cyclerfinder.nbody.propagator import RestrictedNBody


@pytest.mark.slow
def test_state_converges_with_accuracy() -> None:
    r0 = np.array([1.496e8, 0.0, 0.0])
    v0 = np.array([0.0, 24.0, 3.0])
    prop = RestrictedNBody("rebound")
    coarse = prop.propagate(r0, v0, 0.0, 200 * 86400.0, bodies=("E", "M"), accuracy=1e-7)
    fine = prop.propagate(r0, v0, 0.0, 200 * 86400.0, bodies=("E", "M"), accuracy=1e-10)
    assert np.linalg.norm(fine.r_km - coarse.r_km) < 50.0   # converged < 50 km
```

Run → **red** → confirm/extend `accuracy` plumbing → **green** (`-m slow`).

- [ ] `uv run pytest tests/nbody/ -q` (+ `-m slow`) all green/skip-clean.
- [ ] `uv run ruff check . && uv run ruff format --check . && uv run mypy src tests`.

Commit:

```
nbody: GOLDEN GATE 4 — accuracy convergence + Phase A lint/type gate (n-body harness Phase A)
```

---

## Phase B — SILVER rungs

Ingest the two USER-HELD SILVER candidates from the review queue, propagate one
full period in n-body (Sun+E+M+J, design §2), measure terminal closure
degradation + the node-impulse correction ΔV, and **record to the review-queue
audit trail — NEVER auto-promote** (`review_queue.py:62`; the candidates are
user-held, the rung runs and records, the human decides).

### Rung verdict thresholds (human-declared, design §1 / §5)

> **The rung verdict thresholds (verbatim, human-declared per the Jones <200 m/s
> analogue, design §1 gate "small correction ⇒ robust; large ⇒ patched-conic
> artifact"):**
> - **ROBUST** — node-impulse correction ΔV **< 200 m/s** total (the Jones
>   analogue): the conic seed is a faithful shadow of a real n-body-ballistic
>   trajectory.
> - **MARGINAL** — **200 m/s ≤ correction ΔV < 1000 m/s**: real but not free; a
>   powered-cycler-class result, flagged for human judgement.
> - **ARTIFACT** — **correction ΔV ≥ 1000 m/s**, OR the shoot **diverges**: the
>   seed lives only in patched-conic land; strong REJECTED-style evidence.
>
> These thresholds are **recorded, not gating**: the rung writes the tier + the
> ΔV + the terminal closure error to the review-queue audit trail; the human makes
> the promotion/rejection call (golden discipline: no auto-promotion,
> `review_queue.py:82-94`).

### Task B.0 — load the two SILVER candidates from the review queue

**Files:** create `src/cyclerfinder/nbody/rung.py`; test
`tests/nbody/test_rung_ingest.py`.

`load_silver_candidates(path) -> list[ReviewQueueEntry]` via
`review_queue.load_review_queue` (`review_queue.py:121`), filtered to
`verdict_tier == VerdictTier.SILVER.value`. Reconstruct each candidate's seed
(sequence, period from `period_k`, `tof_days`, `vinf_per_encounter_kms`) into the
shooter's seed form (Task C.1 reuse). **Golden discipline:** the candidate
numerics are OUR prior computation (the SILVER tier is unsourced by definition,
`gauntlet.py` SILVER), so the rung asserts *closure/ΔV regime*, never a sourced
value.

#### Failing test

```python
"""N-body Phase B: ingest the two held SILVER candidates (plan Phase B).

NON-GOLDEN: SILVER candidate numerics are OUR computation (unsourced by tier
definition). The rung asserts closure/correction-dV regime, never a sourced V_inf.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from cyclerfinder.nbody.rung import load_silver_candidates


def test_loads_only_silver_entries(tmp_path: Path) -> None:
    qfile = tmp_path / "q.jsonl"
    rows = [
        {"candidate_id": "c1", "verdict_tier": "silver", "sequence": ["E", "M", "E", "E"],
         "period_k": 2, "tof_days": [154.0, 379.0, 1027.0],
         "vinf_per_encounter_kms": [9.75, 13.01, 9.76, 9.75], "bend_feasible": True,
         "signature_hash": "h", "match_outcome": "novel", "known_id": None,
         "superseded_by": [], "max_vinf_kms": 13.01, "model_assumption": "circular",
         "verdict_audit": {}, "panel": {}, "t_added": "2026-06-06T00:00:00Z"},
        {"candidate_id": "c2", "verdict_tier": "bronze", "sequence": ["E", "M", "E"],
         "period_k": 1, "tof_days": [200.0], "vinf_per_encounter_kms": [5.0, 5.0],
         "bend_feasible": True, "signature_hash": "h2", "match_outcome": "novel",
         "known_id": None, "superseded_by": [], "max_vinf_kms": 5.0,
         "model_assumption": "circular", "verdict_audit": {}, "panel": {},
         "t_added": "2026-06-06T00:00:00Z"},
    ]
    qfile.write_text("\n".join(json.dumps(r) for r in rows) + "\n")
    cands = load_silver_candidates(qfile)
    assert [c.candidate_id for c in cands] == ["c1"]  # bronze filtered out
```

Run → **red** → impl `load_silver_candidates` → **green**. Commit:

```
nbody/rung: ingest SILVER candidates from the review queue (n-body harness Phase B)
```

### Task B.1 — single-period n-body propagation of a candidate + terminal closure error

**Files:** `nbody/rung.py`; test `tests/nbody/test_rung_propagate.py` (slow).

`propagate_one_period(entry, ephem, *, bodies=("E","M","J"), accuracy) ->
RungArc`: seed the spacecraft at the home-Earth node (`v_sc = v_planet + vinf`,
the `verify/propagate.py:434`-style reconstruction the design §3 cites), propagate
the full period in n-body, and measure the **terminal closure error** = the
heliocentric state gap between the propagated wrap node and the seeded start node.
NON-GOLDEN (the candidate is unsourced); assert the closure error is *finite and
recorded*, not a value.

#### Failing test (slow)

```python
"""N-body Phase B: full-period propagation + terminal closure error (plan Phase B).

NON-GOLDEN: asserts a finite, recorded closure error for a held SILVER candidate;
the candidate V_inf is OUR computation, never an EXPECTED.
"""
from __future__ import annotations

import math

import pytest

rebound = pytest.importorskip("rebound")

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.nbody.rung import RungArc, propagate_one_period


@pytest.mark.slow
def test_propagate_one_period_reports_finite_closure(silver_fixture) -> None:
    arc = propagate_one_period(silver_fixture, Ephemeris("astropy"), accuracy=1e-10)
    assert isinstance(arc, RungArc)
    assert math.isfinite(arc.terminal_closure_km)
    assert arc.terminal_closure_km >= 0.0
```

(`silver_fixture` builds a `ReviewQueueEntry` from the candidate-1 numerics in
`docs/notes/2026-06-06-silver-candidates-russell-diligence.md`.) Run → **red** →
impl → **green** (`-m slow`). Commit:

```
nbody/rung: single-period n-body propagation + terminal closure error (n-body harness Phase B)
```

### Task B.2 — node-impulse correction ΔV metric (design §3, Q3)

**Files:** create `src/cyclerfinder/nbody/correction_dv.py`; test
`tests/nbody/test_correction_dv.py`.

`node_impulse_correction_dv(seed_nodes, corrected_nodes) -> CorrectionDV`: the
sum of per-node velocity-discontinuity changes between the raw conic seed and the
n-body-corrected solution (design §3, Q3 node-impulse convention). Carries a
`compare_to_maintenance(maintenance_dv_kms)` note documenting the **comparable**
(not identical) mapping to `maintain.py`'s per-synodic value (`maintain.py:24,41`).
Sign- and node-explicit.

#### Failing test

```python
"""N-body Phase B: node-impulse correction-dV convention (plan Phase B; design §3 Q3)."""
from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.nbody.correction_dv import node_impulse_correction_dv


def test_zero_correction_when_seed_equals_corrected() -> None:
    nodes = {"e0": np.zeros(3), "m1": np.array([1.0, 0.0, 0.0])}
    dv = node_impulse_correction_dv(nodes, nodes)
    assert dv.total_kms == pytest.approx(0.0)


def test_correction_is_sum_of_per_node_discontinuity_changes() -> None:
    seed = {"e0": np.zeros(3), "m1": np.zeros(3)}
    corr = {"e0": np.array([0.1, 0.0, 0.0]), "m1": np.array([0.0, 0.05, 0.0])}
    dv = node_impulse_correction_dv(seed, corr)
    assert dv.total_kms == pytest.approx(0.15)
    assert set(dv.per_node_kms) == {"e0", "m1"}
```

Run → **red** → impl → **green**. Commit:

```
nbody/correction_dv: node-impulse correction-dV (comparable to maintain.py) (n-body harness Phase B)
```

### Task B.3 — rung verdict + record to the review-queue audit trail (NEVER promote)

**Files:** `nbody/rung.py`; test `tests/nbody/test_rung_verdict.py`.

`rung_verdict(correction_dv, terminal_closure_km, converged) -> RungVerdict`
applying the human-declared thresholds (ROBUST <200 m/s, MARGINAL <1000 m/s,
ARTIFACT ≥1000 m/s or diverged). `record_rung_result(entry, verdict, path)`
**appends a rung-audit record to the candidate's review-queue audit trail**
(extending `verdict_audit`, `review_queue.py:57`) — it must NOT change
`verdict_tier`, must NOT promote, must NOT write a catalogue row
(`review_queue.py:62-68` `is_catalogue_source()` False). Assert the recording is
non-promoting.

#### Failing test

```python
"""N-body Phase B: rung verdict thresholds + non-promoting audit record (plan Phase B).

The rung RECORDS, the human DECIDES (review_queue golden discipline: no auto-
promotion). This test pins both the thresholds and the no-promotion invariant.
"""
from __future__ import annotations

from cyclerfinder.nbody.rung import RungVerdict, record_rung_result, rung_verdict


def test_thresholds() -> None:
    assert rung_verdict(0.15, terminal_closure_km=10.0, converged=True).tier == "ROBUST"
    assert rung_verdict(0.5, terminal_closure_km=10.0, converged=True).tier == "MARGINAL"
    assert rung_verdict(2.0, terminal_closure_km=10.0, converged=True).tier == "ARTIFACT"
    assert rung_verdict(0.05, terminal_closure_km=1e9, converged=False).tier == "ARTIFACT"


def test_record_does_not_promote(tmp_path, silver_fixture) -> None:
    v = rung_verdict(2.0, terminal_closure_km=1e6, converged=True)
    out = record_rung_result(silver_fixture, v, tmp_path / "audit.jsonl")
    # Tier on the queue entry is untouched; the rung wrote an AUDIT note only.
    assert silver_fixture.verdict_tier == "silver"
    assert out["rung_verdict"] == "ARTIFACT"
    assert out["promoted"] is False
```

Run → **red** → impl → **green**. Commit:

```
nbody/rung: rung verdict + non-promoting review-queue audit record (n-body harness Phase B)
```

### Task B.4 — run the two held candidates + record findings; Phase-B gate

**Files:** test `tests/nbody/test_rung_silver_pair.py` (slow); lint/type;
**a findings note** `docs/notes/2026-06-06-nbody-silver-rungs.md` *(defer the note
if shared-doc concurrency forbids it at run time; flag in the report)*.

Run the rung end-to-end on **both** held SILVER candidates (Sun+E+M+J), record
each verdict + correction ΔV + terminal closure to the audit trail, and write the
findings note. **Expected outcome per the honesty boundary:** both candidates
float at E∞ ~9.7 / M∞ ~12–13 km/s (the high-V∞ basin); the rung likely returns
MARGINAL/ARTIFACT, which is itself a valid recorded finding — the rung's *job* is
to record, not to make them pass. Also record the **Jupiter on/off sensitivity**
(Gate 4 arm) against these real baselines (design §2 standing rule).

- [ ] `uv run pytest tests/nbody/test_rung_silver_pair.py -q -m slow` runs;
  verdicts recorded.
- [ ] `uv run ruff check . && uv run ruff format --check . && uv run mypy src tests`.

Commit:

```
nbody/rung: SILVER rung run on the two held candidates + Jupiter sensitivity (n-body harness Phase B)
```

---

## Phase C — the Jones shooter

Multiple shooting over encounter nodes from a near-miss conic seed. **Gated on
the #135 verdict** (Task C.0 first) and the #133 Phase-3a near-miss survey.

### Task C.0 — INCORPORATE THE #135 VERDICT (checkpoint, FIRST in this phase)

**Files:** read-only review; test `tests/nbody/test_c0_135_checkpoint.py`
(a documentation/guard test); **a checkpoint note**
`docs/notes/2026-06-06-nbody-shooter-135-checkpoint.md` *(defer if concurrency
forbids; flag in report)*.

> **#135-verdict checkpoint (binding).** Before writing any shooter code, read the
> landed #135 like-for-like verdict
> (`docs/notes/2026-06-06-russell12-likeforlike.md`) and decide the seeding
> strategy from it. The verdict at plan-authoring time was **seeding/basin, not
> solver deficiency**: 0 CLOSE-AND-MATCH coplanar-vs-coplanar, the corrector lands
> off-anchor (ach V∞ 9–28 vs src 3–10). **Implication for the shooter:** a
> single-shoot from a naive equispaced/coplanar seed will fall into the same
> high-V∞ basin; the shooter MUST be seeded from the **#133 near-miss conic
> survey** (real near-encounter geometry close to the Jones anchors), NOT from a
> blind scan. If #135's final verdict differs from the authoring-time read
> (e.g. it concludes a solver deficiency), re-derive the seeding strategy here and
> record the change before C.1.

The guard test pins that the shooter module documents its seeding source as the
near-miss survey (not a blind scan), so the #135 lesson cannot be silently
dropped.

```python
"""N-body Phase C.0: the shooter honours the #135 seeding verdict (plan Phase C).

#135 verdict (russell12-likeforlike): the basin problem is SEEDING, not solver.
The shooter must be seeded from the #133 near-miss survey, never a blind scan.
This guard pins that intent in the module docstring so it cannot be lost.
"""
from __future__ import annotations

import inspect


def test_shooter_documents_nearmiss_seeding() -> None:
    from cyclerfinder.nbody import shooter

    doc = inspect.getdoc(shooter) or ""
    assert "near-miss" in doc.lower()
    assert "135" in doc  # cross-reference the verdict that mandated it
```

Run → **red** (no `shooter` module) → create `shooter.py` with the mandated
docstring → **green**. Commit:

```
nbody/shooter: #135-verdict checkpoint — near-miss seeding mandated (n-body harness Phase C)
```

### Task C.1 — node/defect structure + the conic-seed mapping

**Files:** `nbody/shooter.py`; test `tests/nbody/test_shooter_nodes.py`.

Build the multiple-shooting variable vector `x = [{node states}, {node epochs},
{ToFs}]` with the slack-leg/period pin from `correct.py:_reconstruct_tofs`
(`correct.py:77-86`). Map a conic seed in (the **DRIFT** flagged above: read node
V∞ **vectors** from `correct._vinf_nodes`'s `b{i}_in/_out` dict, NOT from a
non-existent `best_cycler.encounters` attribute). Node epochs on the **TDB-J2000
axis** (§0), node states `ephem.state(body, epoch)` + seeded outgoing V∞ (the
`v_sc = v_planet + vinf` reconstruction). The defect for leg i = `(n-body
propagate node_i → node_{i+1}) − node_{i+1} state`.

#### Failing test

```python
"""N-body Phase C: multiple-shooting node/defect structure + conic seed map (plan Phase C).

Reads node V_inf VECTORS from correct._vinf_nodes (b{i}_in/_out), NOT from a
best_cycler attribute (design §3 drift: BallisticClosureResult carries only scalar
vinf_per_encounter_kms; the vectors live in _vinf_nodes).
"""
from __future__ import annotations

import numpy as np

from cyclerfinder.nbody.shooter import build_shooting_vector, defect_count


def test_defect_count_is_legs_times_state_dim() -> None:
    # E-M-E-V-V-E: 5 legs -> 5 full-state defects (6 components each).
    assert defect_count(n_encounters=6) == 5 * 6


def test_shooting_vector_packs_states_epochs_tofs() -> None:
    nodes = {f"b{i}": np.arange(6.0) for i in range(4)}
    epochs = [0.0, 1.0, 2.0, 3.0]
    tofs = [1.0, 1.0, 1.0]
    x = build_shooting_vector(nodes, epochs, tofs, slack_leg=2, period_days=4.0)
    # slack leg eliminated from the free vector.
    assert len(x) == 4 * 6 + 4 + (3 - 1)
```

Run → **red** → impl → **green**. Commit:

```
nbody/shooter: node/defect structure + conic-seed mapping via _vinf_nodes (n-body harness Phase C)
```

### Task C.2 — defect residual + flyby/periodicity constraints (finite-diff Jacobian)

**Files:** `nbody/shooter.py`; test `tests/nbody/test_shooter_residual.py` (slow).

The residual = full-state defects + flyby constraints (`r_p ≥ r_p_safe` via
`PLANETS[body].safe_alt_km`, `constants.py:168`, the `correct.py` bend lessons
carry in) + periodicity wrap. Jacobian by finite difference over REBOUND (Q1
baseline). Divergence is a first-class non-converged record (mirror
`correct.py:354-378`), never an exception.

#### Failing test (slow)

```python
"""N-body Phase C: defect residual is zero on a continuous fixture (plan Phase C)."""
from __future__ import annotations

import numpy as np
import pytest

rebound = pytest.importorskip("rebound")

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.nbody.shooter import defect_residual


@pytest.mark.slow
def test_defect_zero_for_self_consistent_nodes() -> None:
    # Nodes placed ON a real n-body arc -> defects ~ 0 by construction.
    ephem = Ephemeris("astropy")
    res = defect_residual(_arc_consistent_x(ephem), ephem=ephem, bodies=("E", "M"),
                          accuracy=1e-10)
    assert np.max(np.abs(res)) < 1e-3  # km / (km/s) defect floor
```

(`_arc_consistent_x` builds nodes by sampling a single n-body propagation, so the
defects vanish by construction — a self-consistency check, not a sourced anchor.)
Run → **red** → impl → **green** (`-m slow`). Commit:

```
nbody/shooter: defect residual + flyby/periodicity constraints (n-body harness Phase C)
```

### Task C.3 — the multiple-shooting solve (`shoot`) + honest divergence

**Files:** `nbody/shooter.py`; test `tests/nbody/test_shooter_solve.py` (slow).

`shoot(seed, ephem, *, bodies, accuracy) -> ShootResult` running
`least_squares(defect_residual, x0, method="lm")` (mirror `correct.py`), filling
a frozen `ShootResult` (corrected nodes, defect norm, per-node correction ΔV via
Task B.2's `node_impulse_correction_dv`, converged/diverged, integrator settings).
A diverged shoot returns an honest non-converged record (design §3).

#### Failing test (slow)

```python
"""N-body Phase C: multiple-shooting solve drives defects down (plan Phase C).

NON-GOLDEN convergence check: asserts the solver reduces the defect from a
slightly-perturbed self-consistent seed (a solver-health test, not a sourced
rediscovery). The Jones-multiset gate is Task C.4.
"""
from __future__ import annotations

import pytest

rebound = pytest.importorskip("rebound")

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.nbody.shooter import shoot


@pytest.mark.slow
def test_shoot_reduces_defect_from_perturbed_seed() -> None:
    ephem = Ephemeris("astropy")
    result = shoot(_perturbed_consistent_seed(ephem), ephem=ephem,
                   bodies=("E", "M"), accuracy=1e-10)
    assert result.converged
    assert result.defect_norm < result.seed_defect_norm
```

Run → **red** → impl → **green** (`-m slow`). Commit:

```
nbody/shooter: multiple-shooting solve with honest divergence record (n-body harness Phase C)
```

### Task C.4 — the Jones VEM gate (sourced multisets, fold-corrected, unchanged criteria)

**Files:** test `tests/nbody/test_shooter_jones_gate.py` (slow, xfail-first).

> **The Jones gate criteria are UNCHANGED from the conic layer (design §3):**
> sourced multisets only on the EXPECTED side (the catalogue rows'
> `vinf_kms_at_encounters` from AAS 17-577 Tables 2/3), compared as a
> **fold-corrected sorted multiset**, tolerance `VEM_VINF_TOL_KMS = 0.5 km/s` per
> encounter. The shooter output (n-body-ballistic V∞ at the converged nodes) is the
> side under test. **No self-computed value is ever EXPECTED** (golden discipline,
> `feedback_golden_tests_sourced_only`). The gate is **xfail-first**: per the
> honesty boundary + the #135 basin verdict, convergence to the Jones 2.5–3.9 km/s
> regime is NOT assumed. Flip only when the near-miss-seeded shoot reaches the
> sourced multiset within tolerance for at least one Jones member row.

Seeded from the **#133 near-miss survey** (Task C.0 mandate). Mirrors the M-ED
headline-gate pattern (`tests/test_vem_rediscovery.py:1058-1082`).

#### Test (xfail-first)

```python
"""N-body Phase C HEADLINE: Jones VEM ballistic rediscovery in n-body (plan Phase C).

GOLDEN: EXPECTED = catalogue SOURCED vinf_kms_at_encounters (AAS 17-577). The
shooter output is the side under test. Fold-corrected sorted-multiset compare,
VEM_VINF_TOL_KMS = 0.5. xfail until the near-miss-seeded shoot reaches the Jones
basin (honesty boundary + #135 verdict: convergence NOT assumed).
"""
from __future__ import annotations

import pytest
import yaml

from tests._catalogue_loader import CATALOGUE_PATH

VEM_VINF_TOL_KMS = 0.5
_VEM_MEMBERS = ("jones-2017-vem-emevve-outbound", "jones-2017-vem-meevem-inbound")


def _sourced_multiset(entry_id: str) -> list[float]:
    for row in yaml.safe_load(CATALOGUE_PATH.read_text()):
        if row["id"] == entry_id:
            return sorted(float(e["vinf_kms"]) for e in row["vinf_kms_at_encounters"])
    raise AssertionError(f"row {entry_id!r} not found")


@pytest.mark.slow
@pytest.mark.xfail(
    reason="n-body Jones shooter: rediscovery to the sourced AAS 17-577 multiset "
    "within 0.5 km/s. Open until the near-miss-seeded shoot reaches the Jones "
    "2.5-3.9 basin (#135 verdict: seeding/basin problem; honesty boundary). Flip "
    "when >=1 member converges; see plan Phase C Task C.4 + C.5.",
    strict=False,
)
@pytest.mark.parametrize("entry_id", _VEM_MEMBERS)
def test_jones_vem_nbody_rediscovers_sourced_multiset(entry_id: str) -> None:
    expected = _sourced_multiset(entry_id)
    result = ...  # shoot(nearmiss_seed_for(entry_id), ...); helper in Task C.4 fixtures
    assert result.converged
    got = sorted(result.vinf_per_encounter_kms)  # fold-corrected in the helper
    for g, x in zip(got, expected, strict=True):
        assert g == pytest.approx(x, abs=VEM_VINF_TOL_KMS)
```

Run `-m slow`: registers `xfail` (or `xpass` → reviewable finding the gate
converged → Task C.5). The `...` is the near-miss-seed helper, spelled out in
prose: build the cell from the row, pull the #133 near-miss survey seed nearest
the row's sourced anchors, `shoot(...)`. Commit:

```
test: n-body Jones VEM gate — sourced-multiset rediscovery (xfail, near-miss seeded)
```

### Task C.5 — STOP/report branch + Tudat-trigger decision; Phase-C gate

**Files:** test runs; a findings note
`docs/notes/2026-06-06-nbody-jones-shooter.md` *(defer if concurrency forbids)*;
lint/type.

**STOP/report branch (binding, honesty boundary + design §7 risk register).**
After C.4:

- **If the gate converges** (xpass for ≥1 row): remove the `xfail` on that row,
  promote to a strict green gate at `VEM_VINF_TOL_KMS`, record the n-body-ballistic
  Jones solution.
- **If it does NOT converge** (the expected basin risk): **STOP — do not loosen
  the tolerance, do not assert a self-computed value, do not flip the xfail.**
  Leave xfail with a per-row finding (best defect norm, V∞ gap to the Jones basin),
  and record in the findings note as open research. The harness + shooter
  infrastructure is the shipped deliverable; reaching the Jones multiset is the
  open item.

**Tudat-trigger decision (Q1, recorded here).** Record whether finite-diff STMs
over REBOUND were too slow/noisy (the §7 risk: "finite-diff STM noise"). The
trigger to promote consumer 2 to Tudat is **demonstrated**, not anticipated:
finite-diff Jacobian condition number or `least_squares` step failure traced to
Jacobian noise, OR per-iteration cost making the solve impractical (cite the perf
profile). If neither is observed, REBOUND stays (Q1). Record the decision; do NOT
add the Tudat dependency speculatively.

- [ ] `uv run pytest -m "not slow"` green; `uv run pytest tests/nbody/ -m slow`
  evaluated (xfail/xpass per the branch).
- [ ] `uv run ruff check . && uv run ruff format --check . && uv run mypy src tests`.

Commit:

```
test: n-body Jones shooter outcome + Tudat-trigger decision (STOP/report branch) (n-body harness Phase C)
```

### Task C.6 — Phase-C lint/type + full-suite gate

- [ ] `uv run pytest tests/nbody/ -q` (+ `-m slow`) green/xfail-as-designed.
- [ ] ruff + ruff format + mypy clean. Commit if cleanup needed:

```
nbody/shooter: Phase C lint/type gate (n-body harness Phase C)
```

---

## Phase D — GMAT V4 stamp (manual, out-of-CI lane)

A `scripts/gmat_v4_aldrin.py` generator + `scripts/parse_gmat_report.py` parser,
run **manually out of CI** (like the perf-profile reproduction scripts, design
§6 / Q7). Produces a one-line "GMAT reproduces 2.9138 km/s within 5%" record for
the Aldrin row. **No CI test gates on GMAT** (Q7).

### V4 tolerance (declared up front, design Q5)

> **V4 tolerance = 5% (verbatim, human-declared BEFORE any run, design Q5).** The
> GMAT reproduction of our computed **2.9138 km/s** per-synodic maintenance ΔV
> passes V4 iff the GMAT-reported maintenance ΔV is within **±5%** (i.e.
> 2.768–3.060 km/s). Declared up front, **never back-fit**. Golden caveat
> (design §5.5): 2.9138 is OUR computed value, so this is genuine *external*
> validation (independent codebase + ephemeris, `docs/spec.md:402`), unlike the
> in-pipeline consistency check.

### Task D.0 — GMAT script generator for the powered Aldrin reproduction

**Files:** create `scripts/gmat_v4_aldrin.py`; test
`tests/nbody/test_gmat_script_gen.py` (NOT slow — string generation only, no GMAT).

`generate_aldrin_script(out_path, *, epoch_iso, ...)` emits a GMAT `.script`
configuring the Aldrin E-M-E powered-periodic cycler against GMAT's own DE
ephemeris + its standard high-fidelity planetary set (Sun..Jupiter at least,
design §2 V4: use the external tool's standard set). The generator is pure string
templating — testable without GMAT installed.

#### Failing test

```python
"""N-body Phase D: GMAT Aldrin script generation (plan Phase D; no GMAT needed)."""
from __future__ import annotations

from pathlib import Path

from scripts.gmat_v4_aldrin import generate_aldrin_script


def test_script_has_aldrin_force_model_and_epoch(tmp_path: Path) -> None:
    out = tmp_path / "aldrin.script"
    generate_aldrin_script(out, epoch_iso="2030-01-01T00:00:00")
    text = out.read_text()
    assert "ForceModel" in text
    assert "Sun" in text and "Earth" in text and "Mars" in text and "Jupiter" in text
    assert "2030-01-01" in text
```

Run → **red** → impl the template → **green**. Commit:

```
scripts/gmat_v4_aldrin: GMAT script generator for the Aldrin V4 stamp (n-body harness Phase D)
```

### Task D.1 — GMAT report parser + the 5% V4 predicate

**Files:** create `scripts/parse_gmat_report.py`; test
`tests/nbody/test_gmat_parse.py`.

`parse_maintenance_dv(report_text) -> float` extracts the per-synodic maintenance
ΔV from a GMAT report; `v4_pass(gmat_dv, our_dv=2.9138, tol_frac=0.05) -> bool`
applies the declared 5% predicate. Test against a **fixture report string** (no
GMAT). Golden discipline: `our_dv=2.9138` is OUR value and is the *reference being
externally checked*, NOT an EXPECTED-from-source assertion — documented in the
docstring.

#### Failing test

```python
"""N-body Phase D: GMAT report parse + 5% V4 predicate (plan Phase D; design Q5).

2.9138 is OUR computed value; GMAT is the INDEPENDENT external check of it (V4).
The 5% band (2.768-3.060) is declared up front, never back-fit.
"""
from __future__ import annotations

from scripts.parse_gmat_report import parse_maintenance_dv, v4_pass


def test_parse_and_pass_within_5pct() -> None:
    report = "...\nManeuver.TotalDV = 2.95 km/s\n..."
    dv = parse_maintenance_dv(report)
    assert dv == 2.95
    assert v4_pass(dv) is True            # within 5% of 2.9138
    assert v4_pass(3.20) is False         # outside the 2.768-3.060 band
```

Run → **red** → impl → **green**. Commit:

```
scripts/parse_gmat_report: parse + 5% V4 predicate for the Aldrin stamp (n-body harness Phase D)
```

### Task D.2 — run-book documentation (manual lane, NOT a CI test)

**Files:** create `docs/notes/2026-06-06-nbody-gmat-v4-runbook.md` *(defer if
shared-doc concurrency forbids it at run time; flag in the report)*.

Document the **manual** procedure (Q7): install GMAT (version pinned), generate
the Aldrin script (`scripts/gmat_v4_aldrin.py`), run GMAT batch, parse the report
(`scripts/parse_gmat_report.py`), record the one-line V4 result (pass/fail at 5%)
for the Aldrin row. State explicitly: **this is NOT in CI; it is a
once-per-promising-candidate stamp** (Q7). Note the golden caveat (2.9138 is our
value; GMAT is the external check). Run-book documentation **rather than CI
tests** (design Q7).

Commit:

```
docs: GMAT V4 stamp run-book — manual Aldrin reproduction at 5% (n-body harness Phase D)
```

### Task D.3 — Phase-D gate

- [ ] `uv run pytest tests/nbody/test_gmat_*.py -q` green (string-only, no GMAT).
- [ ] ruff + ruff format + mypy clean. Confirm **no CI test invokes GMAT** (Q7).
  Commit if cleanup:

```
scripts: Phase D lint/type gate; confirm GMAT stays out of CI (n-body harness Phase D)
```

---

## Risk register (honest, binding)

1. **Rails-model validity for the correction-ΔV metric.** The planets-on-rails
   restricted n-body has **no spacecraft back-reaction** (exact for a massless
   particle) but also **no true Jacobi/energy invariant** (time-dependent,
   non-autonomous). The correction-ΔV is therefore a *model-relative* number; its
   honesty depends on the §0 anchor (Gate 3) being exact and the body set being
   sensitivity-justified (Gate 4 arm). **Mitigation:** Gates 3 + 4 gate every
   science number; the rung records the body set used; the V4 GMAT stamp is the
   independent backstop for the headline (Aldrin) number.
2. **Finite-diff STM noise (the Tudat trigger).** Finite-difference Jacobians over
   an expensive DE440 propagation are slow and noisy (design §3, perf profile:
   propagation, not the optimiser, is the cost). **Trigger to promote to Tudat
   (Q1, demonstrated not anticipated):** Jacobian condition number / `least_squares`
   step failure traced to Jacobian noise, OR impractical per-iteration cost. If
   neither is observed, REBOUND stays. Recorded in Task C.5; the Tudat dep is NOT
   added speculatively.
3. **REBOUND + rails integration subtleties.** IAS15 is adaptive high-order; the
   external rails force is **velocity-independent** (no drag/SRP), which IAS15
   handles well, but (a) the `additional_forces` callback must read planet states
   at the integrator's *sub-step* times (not just node times) — the §0 reader is
   called per force-eval, which is the perf cost; (b) units (REBOUND is unit-
   agnostic; we fix km, km/s, s and the matching `μ` from `core/constants`); (c)
   energy bookkeeping with an external force needs REBOUND's
   `sim.update_acceleration` discipline. **Mitigation:** Gate 1 (two-body
   reduction) + Gate 2 (Sun-only machine-precision energy) catch a units/force
   wiring bug before any perturber is trusted.
4. **The basin problem (carried from #135).** The whole shooter may fail to reach
   the Jones 2.5–3.9 km/s regime from any seed if the near-miss survey does not
   bracket it. **Mitigation:** Task C.0 mandates near-miss seeding; Task C.5's
   STOP/report branch ships the infrastructure honestly if the basin is not
   reached. The plan does not promise the Jones gate flips green.
5. **GMAT operational brittleness (Q7).** GMAT is a GUI-era app; scripted batch is
   brittle/slow/hard-to-pin. **Mitigation:** it is a manual out-of-CI lane (D.2
   run-book), never a gating check; the generator/parser are unit-tested without
   GMAT installed so the logic is covered even when the app is absent.

---

## Out of scope (explicit carve-outs)

- **`assist` dependency** — not taken (Q2); the ephemeris force is hand-rolled.
  Re-evaluate only if hand-rolled interpolation proves a correctness problem.
- **Tudat** — deferred (Q1); a conditional Phase-C upgrade on demonstrated
  finite-diff STM failure, with its own extra. Not installed now.
- **Independent-kernel rung upgrade** (Q6) — running the rungs on a separately
  downloaded DE440 to push a strong rung toward V4 — a later nice-to-have, not
  this plan.
- **SRP / GR / non-spherical gravity / moons** — OUT (design §2 justifications).
- **GMAT in CI** — OUT (Q7); manual lane only.
- **Powered-cycler BVP generalisation** off the Aldrin lock
  (`bvp.py:92,149`) — the V4 stamp *reproduces* the Aldrin solution; it does not
  generalise the powered solver (that is M7).

---

## Definition of done

0. **Phase A:** `nbody/` package; rebound in the `validation` extra; the §0
   time + frame conversions reuse #129's LSK + `core/ephemeris`'s obliquity; the
   `Propagator` interface + rails callback; **all five golden gates green**
   (two-body reduction, energy floor, DE440 anchor, accuracy convergence, with the
   conservation invariant designed honestly — no fake Jacobi constant). **Green.**
1. **Phase B:** the two held SILVER candidates ingested from the review queue,
   propagated one full period (Sun+E+M+J), correction ΔV + terminal closure
   measured, verdict (ROBUST/MARGINAL/ARTIFACT) **recorded to the review-queue
   audit trail with `promoted=False`** — never auto-promoted. Jupiter sensitivity
   recorded. **Green.**
2. **Phase C:** multiple-shooting solver behind the §3 interface; the
   #135-verdict checkpoint honoured (near-miss seeding mandated); the Jones VEM
   gate exists (xfail-first, sourced multiset EXPECTED, fold-corrected, 0.5 km/s),
   flipped green iff it converges (C.5 STOP/report branch governs the no-converge
   outcome); the Tudat-trigger decision recorded.
3. **Phase D:** GMAT generator + parser (unit-tested without GMAT); the **5%** V4
   predicate declared up front; a run-book documenting the manual out-of-CI lane.
   **No CI test invokes GMAT.**
4. No assertion anywhere uses a value our own code computed as the EXPECTED side —
   only the analytic two-body problem (source-free), DE440 itself (Gate 3), and
   the sourced Jones multisets (Gate C.4). The SILVER candidates are unsourced by
   tier definition, so the rung asserts *regime*, never a value.
5. `uv run pytest -m "not slow"`, ruff, ruff format, mypy all clean.

---

## Self-Review

### Design coverage (every design contract → plan task)

| Design item | Plan location |
|---|---|
| §0 epoch/TDB-J2000 axis, reuse #129 LSK | Phase A Task A.1 |
| §0 frame/obliquity, reuse `core/ephemeris` constant | Phase A Task A.2 |
| §0 same-kernel discipline (shared DE440) | Architecture; Task A.5 anchor |
| §1 single propagator interface (REBOUND-first, Tudat slot) | Phase A Task A.3 |
| §1 REBOUND in venv/CI (Q2 no assist) | Task A.0 (validation extra) |
| §2 rails force model, per-consumer body sets, SRP/GR/moons OUT | Architecture; Task A.3; Risk 1 |
| §3 multiple-shooting node/defect structure | Phase C Tasks C.1-C.3 |
| §3 conic-seed mapping (DRIFT: via `_vinf_nodes`) | Task C.1 (drift flagged) |
| §3 correction-ΔV node-impulse convention (Q3) | Phase B Task B.2 |
| §3 honest divergence (non-converged record) | Tasks B.1, C.3 |
| §4 independence accounting (rungs shared-DE440, never V4, Q6) | Architecture; Task B.3 |
| §5.1 GOLDEN GATE 1 two-body reduction | Task A.4 |
| §5.2 GOLDEN GATE 2 energy/integral (honest invariant) | Task A.5 |
| §5.3 GOLDEN GATE 4 convergence + body sensitivity | Tasks A.6, B.4, C.5 |
| §5.4 GOLDEN GATE 3 Earth/DE440 anchor | Task A.5 |
| §5.5 Aldrin reproduction (consistency vs external) | Phase D (external/GMAT) |
| §6 phasing rungs → shooter → GMAT, §0 first | Phase order; A first |
| Q1 REBOUND-first, Tudat on need | Task A.3; C.5 trigger; Risk 2 |
| Q2 hand-rolled force, no assist | Task A.3; Out of scope |
| Q3 comparable not identical to maintain.py | Task B.2 |
| Q4 §5.3 sensitivity governs body set | Tasks A.6, B.4 |
| Q5 V4 tol 5% up front | Phase D V4-tolerance block; Task D.1 |
| Q6 rungs cross-check only, never V4 | Architecture; Task B.3 |
| Q7 GMAT manual out-of-CI | Phase D; Tasks D.2, D.3 |
| #135 verdict checkpoint at shooter start | Phase C Task C.0 |
| #133 near-miss survey feeds seeds | Tasks C.0, C.4 |
| Jones gate criteria unchanged (sourced, fold-corrected) | Task C.4 |

### Placeholder scan

No `TODO`/`FIXME`/`...`-as-impl placeholders in production-code instructions. The
single `...` is inside the Task C.4 *test sketch* (the near-miss-seed helper),
explicitly marked "helper in Task C.4 fixtures" with the construction recipe
spelled out in prose (build cell from row → nearest #133 near-miss seed →
`shoot`). All paths are concrete; all commit messages exact; all golden gates and
the rung thresholds quoted verbatim from the design / declared inline.

### Type consistency

`NBodyArc`, `RungArc`, `CorrectionDV`, `RungVerdict`, `ShootResult` are frozen
dataclasses (mirror `BallisticClosureResult`'s frozen-result discipline,
`correct.py:39`). `Propagator.propagate` signature is stable across backends
(`bodies`, `accuracy`). `load_silver_candidates` returns `list[ReviewQueueEntry]`
(the live type, `review_queue.py:36`). The shooter reads node vectors from
`correct._vinf_nodes` (scalars-only `BallisticClosureResult` drift handled). The
GMAT lane functions are pure (string in/out), testable without GMAT.

### Design claims that did NOT hold against live code (flagged)

1. **`best_cycler.encounters[].vinf_in/vinf_out` (design §3).** The live
   `BallisticClosureResult` (`correct.py:39-50`) carries **scalar**
   `vinf_per_encounter_kms`, not per-node vectors, and has **no** `best_cycler`
   attribute. The per-node V∞ *vectors* exist only inside `_vinf_nodes`'s
   `b{i}_in/_out` dict (`correct.py:86`). Task C.1 reads from `_vinf_nodes`
   accordingly. The same drift appears in the M-ED headline-gate test sketch
   (`2026-06-05-m-ed-ballistic-corrector.md:1077-1080`), which is itself a sketch
   — the live `correct.py` is the authority.
2. **SILVER candidate regime vs the Jones target (design §1 rung framing).** The
   design frames the rungs as the "Jones <200 m/s analogue" expecting possibly-
   robust candidates. The live candidates float at **E∞ 9.62–9.75 / M∞
   12.06–13.01 km/s** (the high-V∞ basin,
   `2026-06-06-silver-candidates-russell-diligence.md`), and #135 shows this
   family lands off-anchor by construction. So the *expected* rung verdict is
   MARGINAL/ARTIFACT, not ROBUST — the plan's honesty boundary + the rung
   thresholds make that a valid recorded outcome rather than a failure. Not a code
   contradiction, but a science-expectation flag the design's optimistic framing
   understates.
3. **Aldrin 2.9138 km/s provenance.** Confirmed it is **our computed** value
   (recorded only in `docs/notes/2026-06-06-performance-profile.md` + the design),
   **not** a published Aldrin maintenance ΔV (the catalogue records only a
   turn-angle test for classic Aldrin, design Q5). The plan treats it as the
   *reference under external check*, never an EXPECTED-from-source — Phase D
   docstrings state this explicitly. No drift, but the golden caveat is load-
   bearing and is preserved.
