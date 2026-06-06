# M-3D — Inclination lift (Venus i=3.39°, Mars i=1.85°)

> **For agentic workers:** REQUIRED SUB-SKILL: `superpowers:subagent-driven-development`
> or `superpowers:executing-plans`. Checkbox steps; strict TDD (write failing
> test → run **red** → minimal impl → run **green** → commit). Work on `main` —
> **do NOT branch** (project rule). **Do NOT commit without the user's review**
> while this plan is being authored under a docs-only mandate; the commit
> messages below are the messages to *use* when the implementation phase runs.
> uv-managed venv (no pip). Lint/type gate before **every** commit:
> `uv run ruff check .` · `uv run ruff format --check .` · `uv run mypy src tests`.
> Fast suite: `uv run pytest -m "not slow"`; any DE440 cross-check is `slow`.
>
> This plan is the task-level expansion of the **approved** design
> `docs/superpowers/specs/2026-06-05-m-3d-inclination-lift-design.md` (read it,
> **including its Approval section** — all six recommendations were accepted).
> Where this plan and the design could diverge, **the design + its Approval
> win.** This plan implements **Approach 1** (activate + frame-fix): the F3
> vector-ω⃗ closure-frame primitive with F1 interim acceptance, the
> `circular-inclined` fidelity rung, `linkable_3d` consistency wiring, and the
> flyby bend-decomposition **diagnostic only** — the inclined closure *attempt*
> and the optimiser b-plane DOF stay with M-ED (Approval Q4/Q5).
>
> **Concurrency note (binding):** M-ED Phases 1-5 are landing on `main` *now*
> and own `search/correct.py`, `search/descriptor.py`, `search/seed_ladder.py`,
> the `mode="ballistic"` branch of `search/optimize.py`, and
> `data/discover.py`'s V3 gate. **M-3D must not touch any of those.** M-3D's
> production write-set is `core/frames.py`, `core/ephemeris.py` (additive
> opt-in only), `core/flyby.py` (new pure helper), `data/provenance.py`,
> `data/catalogue.schema.json`, `verify/fidelity.py`, and `search/sequence.py`
> (`tisserand_feasible` already threads `linkable_3d` — M-3D only pins it).

---

## Goal

Turn on the inclination substrate that is **already built and sourced** in the
codebase (the `_InclinedCircularBackend`, the sourced J2000 inc/Ω in `PLANETS`'
comments, the `linkable_3d` predicate) behind an **explicit opt-in that never
mutates the live coplanar `PLANETS`**, and fix the **one** place the coplanar
assumption is still load-bearing in idealized mode — the closure frame
(`model/cycler.py:closure_residual` → `core/frames.py:to_rotating`).

The crux deliverable is the **F3 vector-ω⃗ generalisation of
`to_rotating`/`from_rotating`**: a strict superset of today's scalar-z form
whose coplanar limit (`ω⃗ = ω ẑ`) reproduces every M3 frame golden **bit for
bit**. On top of that primitive M-3D ships (a) a sanctioned inclined-state
accessor, (b) the `circular-inclined` fidelity rung wired into
`solve_at_fidelity`, (c) `linkable_3d` consistency pins in `tisserand_feasible`,
and (d) the flyby bend-decomposition **diagnostic** (in/out-of-plane
attribution) on `core/flyby.py` outputs.

**Honesty boundary (binding, carried from design §6 and project memory
`feedback_golden_tests_sourced_only.md`).** There is **no published "inclined
cycler closes at X km/s" anchor** at the inclined-circular rung. Almost every
M-3D test is therefore a **physics invariant** (coplanar limit, round-trip
identity, energy/momentum, node geometry, Tisserand `cos(i)→1` reduction) or a
**reviewed pinned diff** (the M8-Core coplanar `tisserand_feasible` baselines).
The only *sourced* 3-D anchors are the inclinations/nodes themselves (Standish &
Williams, `core/constants.py:179-211` comments) — a 3-D-state assertion that the
inclined backend reproduces V at i=3.39°, M at i=1.85° at the node is legitimate
because **the input is the source**. No closure or V∞ value M-3D computes is ever
the EXPECTED side of a golden.

---

## Architecture

### The closure-frame primitive (design §2, the crux)

`closure_residual` (`model/cycler.py:208-252`) transforms the spacecraft
departure/arrival velocities into a synodic rotating frame via `to_rotating`
(`core/frames.py:90-146`) at a scalar `omega` about ecliptic **+z**. For a
coplanar-circular system the orbit plane *is* the ecliptic, so this is exact;
for an inclined orbit the spacecraft has genuine z-motion the +z frame does not
absorb, so a geometrically-perfect inclined cycler shows a spurious residual.

F3 replaces the implicit `ω⃗ = ω ẑ` with an explicit angular-velocity **vector**
`ω⃗` so the Coriolis term `v − ω⃗ × r` becomes a full 3-D cross product and the
position rotation is a rotation about the unit `ω̂`. **Backward-compatibility is
the binding constraint:** the existing scalar signature
`to_rotating(r, v, t, omega: float)` is called by every M3 golden
(`tests/core/test_frames.py`) and by `closure_residual`
(`cycler.py:248-251`). The scalar path must stay byte-identical. The chosen
shape is a **new vector entry point** `to_rotating_omega_vec(r, v, t, omega_vec)`
plus a one-line internal reduction so the scalar `to_rotating` delegates to it
with `omega_vec = (0, 0, omega)`; the scalar public signature is unchanged.

The anchor `ω⃗` for an inclined cycler is the home body's instantaneous orbit
normal scaled to its mean motion: `ω⃗ = (r×v)/|r|²` evaluated on the
anchor body's heliocentric `(r,v)` (this is exactly what `synodic_omega_dynamic`
computes the z-component of, `core/frames.py:286-297`). A new
`synodic_omega_vec(anchor_body, t_sec, ephem)` returns that vector; for the
coplanar circular backend it equals `(0, 0, synodic_omega(anchor_body))` to
float precision (a tested degeneracy, mirroring
`synodic_omega_dynamic`'s circular degeneracy gate).

### Opt-in inclined ephemeris (design §3, §7 Approach 3 rejection)

`_InclinedCircularBackend(planets=...)` already accepts an injected `PlanetData`
dict (`core/ephemeris.py:129-130`) and is auto-selected by `_CircularBackend`
only when `planet.inc_deg != 0.0` (`core/ephemeris.py:106-111`, an **exact**
float gate). The live `PLANETS` keep `inc_deg=lan_deg=0.0`
(`core/constants.py:181-186, 209-211`) so the default stays byte-identical.
M-3D adds a **module-level factory** `inclined_planets()` that returns a *new*
dict — a copy of `PLANETS` with each body's sourced J2000 `inc_deg`/`lan_deg`
filled from Standish & Williams — and a convenience
`Ephemeris.inclined_circular()` constructor that builds an `Ephemeris` whose
backend is `_InclinedCircularBackend(inclined_planets())`. **`PLANETS` is never
mutated.** (Note: `PLANETS` already carries sourced non-zero `ecc`,
`constants.py:180,197,208`, but `_circular_inplane_state` ignores it,
`ephemeris.py:72-90` — so the inclined rung stays circular; eccentricity is the
separable follow-on per Approval Q3.)

### The `circular-inclined` fidelity rung (design §3, Approval Q2)

`Fidelity` is `Literal["circular-coplanar", "analytic-ephemeris", "real-de440"]`
(`data/provenance.py:95`) with a *parallel* runtime frozenset `_FIDELITIES`
(`data/provenance.py:109`) and `is_fidelity` (`:117-119`). The same three values
are an `enum` in **two** places in `data/catalogue.schema.json` (`orbit_fidelity`
`:59`, `vinf_fidelity` `:64`). M-3D widens **all** of these to include
`circular-inclined` (inserted between `circular-coplanar` and
`analytic-ephemeris`, matching the ladder ordering documented at
`provenance.py:14-15`). `solve_at_fidelity` (`verify/fidelity.py:111-182`) gains
a **third resolvable rung** dispatching to a new `_solve_inclined` that builds
the inclined ephemeris (above) and runs the same coplanar resonance construction
*or* the inclined-circular state assertions — see §Phase 4 for the exact wiring.
The rung is **additive**: it does not redefine the Forge's persistence
semantics (`fidelity_persistence` is untouched).

### `linkable_3d` consistency (design §5)

`tisserand_feasible(cell, vinf_cap, ephem=None)` (`search/sequence.py:276-366`)
**already** consults `linkable_3d` when an `ephem` is supplied
(`sequence.py:329` `use_3d = ephem is not None`; `:356`). M-3D adds **no new
wiring** here; it (a) pins the coplanar baseline diff the M8-Core tests
explicitly await (`tests/search/test_sequence_multibody.py:1-7`: "so the M-3D
inclination lift is a reviewed diff"), and (b) documents the honest limit
(`linkable_3d` tests energetic linkability with an inclination budget, **not**
node compatibility — design §5) at the gate.

### Flyby bend-decomposition diagnostic (design §6, Approval Q4 — diagnostic only)

`core/flyby.py` is already plane-free (`bend_angle` is the unsigned 3-D angle,
`flyby.py:117-141`). M-3D adds **one pure helper**
`bend_decompose(vin_vec, vout_vec, orbit_normal) -> (delta_inplane_rad,
delta_outofplane_rad)` that splits the V∞-in→V∞-out rotation into the component
in the orbit plane (changes a,e) and the component out of plane (changes i,
node), so a reviewer can *see* Venus doing the plane-change work the spec
predicts. **No optimiser DOF** — no `optimize.py` change, no new continuous
variable. Exploitation (the b-plane DOF) is M-ED's (Approval Q4).

---

## Tech stack

Python 3.11, numpy, pytest. uv-managed venv. Extended production modules:
`src/cyclerfinder/core/frames.py` (vector-ω⃗ primitive + `synodic_omega_vec`),
`src/cyclerfinder/core/ephemeris.py` (additive `inclined_planets()` +
`Ephemeris.inclined_circular()`), `src/cyclerfinder/core/flyby.py`
(`bend_decompose`), `src/cyclerfinder/data/provenance.py` (Fidelity widening),
`src/cyclerfinder/data/catalogue.schema.json` (two enum widenings),
`src/cyclerfinder/verify/fidelity.py` (`circular-inclined` rung). Pin-only:
`tests/search/test_sequence_multibody.py`. New tests under `tests/core/`,
`tests/data/`, `tests/verify/`. No new third-party dependency.

---

## Design references (verified against live code, 2026-06-06)

- §0 / §1 — `_InclinedCircularBackend` rotates by `R_z(+lan) @ R_x(-inc)`,
  injectable `planets` dict — **live at `core/ephemeris.py:114-169`** (design
  said `:114-169`; **holds**). Auto-select gate `inc_deg != 0.0` exact —
  **`core/ephemeris.py:106-111`** (design `:106-111`; **holds**).
- §1 — sourced inc/lan in `PLANETS` comments, live values 0.0 —
  **`core/constants.py:181-186` Venus, `:209-211` Mars** (design cited
  `:181-186,209-212`; **holds**, off by one line on Mars). **NEW SINCE DESIGN:**
  `ecc` is now a *live non-zero* field (`constants.py:180,197,208`) but is
  **ignored** by `_circular_inplane_state` (`ephemeris.py:72-90`) — the circular
  states are still pure circles; the inc-only scope is unaffected. Flagged.
- §1 / §2 — closure frame load-bearing: `closure_residual` →
  `to_rotating` about scalar +z — **`model/cycler.py:208-252`**, transform
  **`core/frames.py:90-146`** (design `cycler.py:235-252`, `frames.py:90-146`;
  **holds** — the method body spans 208-252, the `to_rotating` calls are at
  `:248-251`).
- §2 — `synodic_omega` scalar +z — **`core/frames.py:203-230`** (design
  `:203-230`; **holds**). `synodic_omega_dynamic` computes `(r×v)_z/|r|²` —
  **`core/frames.py:286-297`** (design `:286-297`; **holds**).
- §3 — Fidelity literal + `_FIDELITIES` + `is_fidelity` —
  **`data/provenance.py:95,109,117-119`** (design `provenance.py:95-103`;
  **holds**, the frozenset is a separate object that must also be widened).
  Fidelity enum in schema — **`data/catalogue.schema.json:59,64`** (design said
  "schema 4.5's tag enums — check"; **confirmed: two enum arrays**, labelled
  schema v4.4 in the field descriptions).
- §3 — `solve_at_fidelity` raises `FidelityRungUnavailableError` for
  `analytic-ephemeris` — **`verify/fidelity.py:61-70,167-172`** (design
  `:62-68,157`; **holds**).
- §5 — `linkable_3d` sourced predicate, `i_sc_max_deg=30.0` default,
  monotonicity guarantee — **`search/tisserand.py:471-680`, default `:476`,
  monotonicity `:548-566`** (design `:471-680`, `:476`, `:544-559`; **holds**).
  `tisserand_feasible` already threads it via `ephem` — **`search/sequence.py:276-366`,
  `use_3d` `:329`** (design said `sequence.py:300-363`; **MOVED — now 276-366**,
  the threading is live and exactly as described).
- §6 — `core/flyby.py` plane-free: `bend_angle` `:117-141`, `max_bend` `:40`,
  `flyby_dv` `:182` (design `:117-141,182-231`; **holds**).
- §6 — pinned coplanar baseline awaiting M-3D diff —
  **`tests/search/test_sequence_multibody.py:1-7`** (design `:1-7`; **holds**,
  verbatim "so the M-3D inclination lift is a reviewed diff").

### Design claims that NO LONGER HOLD / shifted (flagged)

1. **`PLANETS.ecc` is now live and non-zero** (`constants.py:180,197,208`). The
   design treated eccentricity as wholly absent. It is carried as a record field
   but **not applied** by the circular backend (`_circular_inplane_state` ignores
   it), so the inc-only scope is intact. No action beyond noting it in the rung
   docstring ("circular-inclined = real i/Ω, mean sma, **e ignored**").
2. **`tisserand_feasible` moved** from the design's `sequence.py:300-363` to
   **`:276-366`**; the `linkable_3d` threading is already present (design §1
   item said "already threads it" — confirmed). M-3D's task here is **pin +
   document only**, not wire.
3. **The Fidelity enum lives in `data/catalogue.schema.json` (two arrays) AND a
   parallel runtime frozenset `_FIDELITIES`**, not only the `Literal`. The design
   said "widen the literal … the Fidelity enum appears in schema". Confirmed and
   enumerated: the widening touches **four** sites (Literal, frozenset, two
   schema enums) plus docstrings. Spelled out in Phase 4.

---

## Phasing (independently shippable)

| Phase | Theme | Tasks | Depends on |
|---|---|---|---|
| **1** | Vector-ω⃗ closure-frame primitive (`core/frames.py`) | 1.0–1.4 (5) | — (foundational) |
| **2** | Opt-in inclined ephemeris (`core/ephemeris.py`) | 2.0–2.2 (3) | — |
| **3** | Flyby bend-decomposition diagnostic (`core/flyby.py`) | 3.0–3.1 (2) | — |
| **4** | `circular-inclined` fidelity rung (provenance + schema + fidelity) | 4.0–4.3 (4) | Phase 2 |
| **5** | Tisserand consistency pins + honest-limit docs (`search/sequence.py`) | 5.0–5.1 (2) | — |

Each phase is shippable alone. The crux (Phase 1) is foundational but the others
do not import its new symbol unless they validate inclined closure (they don't —
inclined closure is M-ED's, Approval Q4). **Total: 16 tasks.**

---

## Phase 1 — Vector-ω⃗ closure-frame primitive (`core/frames.py`)

The heart of M-3D. Generalise `to_rotating`/`from_rotating` to a vector `ω⃗`
**as a strict superset** that reproduces every M3 golden bit-for-bit.

### The coplanar-limit golden gate (verbatim — the binding criterion)

> **For every `(r, v, t, omega)` tuple exercised by the M3 frame goldens in
> `tests/core/test_frames.py` (`test_round_trip_identity`,
> `test_omega_zero_is_identity`,
> `test_circular_orbit_is_stationary_in_its_own_frame`), the new vector entry
> point `to_rotating_omega_vec(r, v, t, omega_vec=(0.0, 0.0, omega))` MUST return
> arrays that are `numpy.array_equal`-identical (bit-for-bit, not `approx`) to the
> existing scalar `to_rotating(r, v, t, omega)`, and likewise
> `from_rotating_omega_vec` vs `from_rotating`. The scalar public functions keep
> their exact `(r, v, t, omega: float)` signature and their existing
> implementation is preserved by delegating to the vector form with
> `omega_vec = (0, 0, omega)` only if that delegation is itself bit-identical;
> otherwise the scalar bodies are left untouched and the vector form is a
> separate code path proven equal by this gate. No M3 golden's EXPECTED side
> changes. The EXPECTED side of this gate is the scalar function's own output
> (a regression invariant, not a sourced anchor) — permitted because it pins a
> refactor-equivalence, exactly the `to_rotating`-`array_equal` discipline the
> design §6 names as item (b).**

### Task 1.0 — vector-ω⃗ `to_rotating_omega_vec` reduces to scalar bit-for-bit

**Files:** `src/cyclerfinder/core/frames.py`; test
`tests/core/test_frames_omega_vec.py`.

#### Failing test — `tests/core/test_frames_omega_vec.py`

```python
"""M-3D Phase 1: vector-omega to_rotating is a bit-exact superset (plan §1).

COPLANAR-LIMIT GOLDEN GATE: with omega_vec = (0,0,omega) the vector form must be
numpy.array_equal to the scalar to_rotating across the M3 golden inputs. EXPECTED
side = the scalar function's own output (a refactor-equivalence regression pin,
design §6(b)), never a sourced value.
"""
from __future__ import annotations

import numpy as np

from cyclerfinder.core.frames import (
    synodic_omega,
    to_rotating,
    to_rotating_omega_vec,
)

_RS = [
    np.array([1.496e8, 0.0, 0.0]),
    np.array([0.0, 2.279e8, 0.0]),
    np.array([1.0e8, -1.0e8, 3.0e7]),  # genuine z-component
]
_VS = [
    np.array([0.0, 29.78, 0.0]),
    np.array([-24.07, 0.0, 0.0]),
    np.array([10.0, 10.0, 1.0]),
]


def test_vector_omega_z_only_matches_scalar_bit_for_bit() -> None:
    omega = synodic_omega("E")
    for omega_val in (0.0, omega, -omega, 2.0 * omega):
        for r, v in zip(_RS, _VS, strict=True):
            for t_sec in (0.0, 1.0e6, -3.7e6):
                r_s, v_s = to_rotating(r, v, t_sec, omega_val)
                r_x, v_x = to_rotating_omega_vec(
                    r, v, t_sec, np.array([0.0, 0.0, omega_val])
                )
                assert np.array_equal(r_s, r_x)
                assert np.array_equal(v_s, v_x)
```

Run: `uv run pytest tests/core/test_frames_omega_vec.py -q` → **red** (no symbol).

#### Minimal impl — `core/frames.py`

Add `to_rotating_omega_vec(r_inertial, v_inertial, t, omega_vec)`. To guarantee
bit-equality with the scalar form when `omega_vec = (0,0,ω)`, the impl must
**reproduce the scalar arithmetic exactly** in that case: rotate the position
about the unit `ω̂` by `−|ω⃗|·t` using Rodrigues, and apply `v − ω⃗ × r` then the
same rotation. For the pure-z case Rodrigues reduces algebraically to the scalar
`R(−θ)` about +z and `ω⃗ × r = (−ω r_y, ω r_x, 0)` — but float rounding can
differ. **Therefore branch:** if `omega_vec[0] == 0.0 and omega_vec[1] == 0.0`,
call the scalar `to_rotating(..., float(omega_vec[2]))` directly (guarantees
`array_equal`); else run the general Rodrigues path. This keeps the gate exact
and adds the genuine 3-D path for non-z `ω⃗`.

Run → **green**. Lint/type. Commit:

```
core/frames: vector-omega to_rotating_omega_vec (z-only bit-exact) (M-3D Phase 1)
```

### Task 1.1 — `from_rotating_omega_vec` round-trip identity + scalar bit-equality

**Files:** `core/frames.py`; test `tests/core/test_frames_omega_vec.py`.

#### Failing test (append)

```python
def test_vector_omega_roundtrip_identity_general_axis() -> None:
    """from(to(x)) == x for a TILTED omega vector (general 3-D frame)."""
    from cyclerfinder.core.frames import from_rotating_omega_vec

    omega_vec = np.array([0.02, -0.015, 0.11]) * synodic_omega("E")
    for r, v in zip(_RS, _VS, strict=True):
        for t_sec in (0.0, 5.0e5, -2.0e6):
            r_rot, v_rot = to_rotating_omega_vec(r, v, t_sec, omega_vec)
            r_back, v_back = from_rotating_omega_vec(r_rot, v_rot, t_sec, omega_vec)
            r_mag = float(np.linalg.norm(r))
            assert float(np.linalg.norm(r_back - r)) / r_mag < 1e-10


def test_vector_omega_z_only_inverse_matches_scalar_bit_for_bit() -> None:
    from cyclerfinder.core.frames import from_rotating, from_rotating_omega_vec

    omega = synodic_omega("E")
    for r, v in zip(_RS, _VS, strict=True):
        r_s, v_s = from_rotating(r, v, 1.0e6, omega)
        r_x, v_x = from_rotating_omega_vec(r, v, 1.0e6, np.array([0.0, 0.0, omega]))
        assert np.array_equal(r_s, r_x)
        assert np.array_equal(v_s, v_x)
```

Run → **red** → impl `from_rotating_omega_vec` (same z-only branch to the scalar
`from_rotating`; general Rodrigues `R(+θ)` about `ω̂` then add back `ω⃗ × r`) →
**green**. Commit:

```
core/frames: from_rotating_omega_vec inverse + round-trip gate (M-3D Phase 1)
```

### Task 1.2 — `synodic_omega_vec` (anchor-body orbit-normal angular velocity)

**Files:** `core/frames.py`; test `tests/core/test_frames_omega_vec.py`.

`synodic_omega_vec(anchor_body, t_sec, ephem) -> Vec3` returns
`(r×v)/|r|²` on the anchor body's heliocentric `(r,v)` — the full vector whose
z-component is `synodic_omega_dynamic` (`frames.py:286-297`). For the **coplanar
circular** backend this is `(0, 0, synodic_omega(anchor_body))` to float
precision; for the **inclined** backend it tilts to the body's orbit normal.

#### Failing test (append)

```python
def test_synodic_omega_vec_is_z_only_for_coplanar_circular() -> None:
    from cyclerfinder.core.ephemeris import Ephemeris
    from cyclerfinder.core.frames import synodic_omega_vec

    ephem = Ephemeris(model="circular")
    w = synodic_omega_vec("E", 1.234e7, ephem)
    assert abs(w[0]) < 1e-18 and abs(w[1]) < 1e-18
    assert w[2] == __import__("pytest").approx(synodic_omega("E"), rel=1e-9)


def test_synodic_omega_vec_tilts_for_inclined_backend() -> None:
    """Inclined Venus -> non-zero in-plane omega components (orbit normal tilt)."""
    import pytest

    from cyclerfinder.core.ephemeris import Ephemeris
    from cyclerfinder.core.frames import synodic_omega_vec

    ephem = Ephemeris.inclined_circular()  # Phase 2 accessor
    w = synodic_omega_vec("V", 0.0, ephem)
    # Venus i=3.39 deg -> the orbit-normal tilt off ecliptic z is ~3.39 deg.
    tilt = float(np.degrees(np.arccos(abs(w[2]) / np.linalg.norm(w))))
    assert tilt == pytest.approx(3.39467605, abs=0.05)
```

> **Phase ordering:** the second test imports `Ephemeris.inclined_circular`
> (Phase 2). If running Phase 1 strictly before Phase 2, mark that single test
> `@pytest.mark.skip(reason="needs Phase 2 inclined_circular accessor")` and
> un-skip it in Task 2.2. The first (coplanar) test gates Task 1.2 alone.

Run → **red** → impl `synodic_omega_vec` (lift the `(r×v)/|r|²` arithmetic from
`synodic_omega_dynamic` but return all three components; reuse `ephem.state`) →
**green**. Commit:

```
core/frames: synodic_omega_vec anchor orbit-normal angular velocity (M-3D Phase 1)
```

### Task 1.3 — re-run the M3 frame goldens (no-change characterisation)

**Files:** none (run-only); the M3 goldens are
`tests/core/test_frames.py::test_round_trip_identity`,
`::test_omega_zero_is_identity`,
`::test_circular_orbit_is_stationary_in_its_own_frame`,
`::test_synodic_omega_earth_matches_constants`,
`::test_synodic_omega_anchors_on_named_body`,
`::test_synodic_omega_unknown_body_raises_keyerror`.

- [ ] `uv run pytest tests/core/test_frames.py -q` → **all green, unchanged**
  (the scalar functions and `synodic_omega` were not modified; this is the proof
  the superset did not disturb the coplanar form).
- [ ] Confirm `closure_residual` (`model/cycler.py:208-252`) still imports and
  calls the **scalar** `to_rotating` — M-3D does **not** rewire
  `closure_residual` onto the vector form (inclined closure is M-ED's, Approval
  Q4). `uv run pytest tests/model/test_cycler.py -q` green.

No commit (verification step) unless a stray import needs cleanup.

### Task 1.4 — Phase-1 lint/type gate

- [ ] `uv run pytest tests/core/test_frames_omega_vec.py tests/core/test_frames.py -q`
- [ ] `uv run ruff check . && uv run ruff format --check . && uv run mypy src tests`

Commit if any cleanup:

```
core/frames: Phase 1 lint/type gate (M-3D)
```

---

## Phase 2 — Opt-in inclined ephemeris (`core/ephemeris.py`)

A sanctioned way to get inclined states **without mutating `PLANETS`**
(design §7 Approach 3 rejection).

### Task 2.0 — `inclined_planets()` factory (sourced inc/Ω, copy not mutate)

**Files:** `src/cyclerfinder/core/ephemeris.py` (or `core/constants.py` if the
factory belongs with the data — keep it in `ephemeris.py` to avoid a second
writer on `constants.py`); test `tests/core/test_inclined_planets.py`.

`inclined_planets() -> dict[str, PlanetData]` returns a **new** dict, each body a
`dataclasses.replace(PLANETS[code], inc_deg=..., lan_deg=...)` with the sourced
J2000 values (Standish & Williams Table 1, the values quoted in the
`constants.py:184,210` comments): Venus `inc=3.39467605, lan=76.67984255`; Mars
`inc=1.84969142, lan=49.55953891`; Earth stays `0.0/0.0` (ecliptic-defining).

#### Failing test

```python
"""M-3D Phase 2: opt-in inclined planet table (plan §2). SOURCED anchor: the
inc/Ω values are Standish & Williams Table 1 (the EXPECTED side IS the source).
"""
from __future__ import annotations

import pytest

from cyclerfinder.core.constants import PLANETS
from cyclerfinder.core.ephemeris import inclined_planets


def test_inclined_planets_carries_sourced_inc_lan() -> None:
    inc = inclined_planets()
    assert inc["V"].inc_deg == pytest.approx(3.39467605)
    assert inc["V"].lan_deg == pytest.approx(76.67984255)
    assert inc["M"].inc_deg == pytest.approx(1.84969142)
    assert inc["M"].lan_deg == pytest.approx(49.55953891)


def test_inclined_planets_does_not_mutate_live_PLANETS() -> None:
    _ = inclined_planets()
    assert PLANETS["V"].inc_deg == 0.0
    assert PLANETS["M"].inc_deg == 0.0  # live coplanar default UNTOUCHED
```

Run → **red** → impl `inclined_planets()` (module-level, `dataclasses.replace`,
sourced constants as named module locals with a Standish & Williams comment) →
**green**. Commit:

```
core/ephemeris: inclined_planets() sourced opt-in table (no PLANETS mutation) (M-3D Phase 2)
```

### Task 2.1 — `Ephemeris.inclined_circular()` constructor

**Files:** `core/ephemeris.py`; test `tests/core/test_inclined_planets.py`.

The live `Ephemeris.__init__` (`ephemeris.py:256-263`) accepts only
`"circular"`/`"astropy"` and hard-raises otherwise (`:261-262`); the only state
interface is `self._backend.state` (`:289`). Two equally valid shapes — **pick
the model-branch form** for symmetry with the existing dispatch:

1. **Preferred:** add an `elif model == "inclined-circular":
   self._backend = _InclinedCircularBackend(inclined_planets())` branch to
   `__init__`, and a thin `@classmethod inclined_circular(cls)` returning
   `cls(model="inclined-circular")`. The existing `"circular"`/`"astropy"`
   branches are **untouched** (byte-identical for every existing caller); the
   `model` property then reports `"inclined-circular"`.
2. Alternative (if a model-string addition is undesirable): a `@classmethod`
   that builds `cls.__new__(cls)` and sets `_backend` /`_model` directly.

The default `Ephemeris(model="circular")` path is byte-identical either way.

#### Failing test (append)

```python
def test_inclined_circular_state_sits_on_node_with_z_excursion() -> None:
    """At t=0 each inclined body is on its ascending node (z==0); off-node it has
    a z-component of order a*sin(inc). SOURCED: the inclination is the anchor."""
    import numpy as np

    from cyclerfinder.core.ephemeris import Ephemeris

    ephem = Ephemeris.inclined_circular()
    r0, _ = ephem.state("V", 0.0)
    assert abs(float(r0[2])) < 1.0  # on the node at t=0 (km)
    # Quarter period later the z-excursion ~ a * sin(i). Sign/scale only.
    r_q, _ = ephem.state("V", 0.0 + 56.0 * 86400.0)  # ~quarter of Venus year
    assert abs(float(r_q[2])) > 1.0e6  # genuine out-of-ecliptic motion
```

Run → **red** → impl the constructor → **green**. Commit:

```
core/ephemeris: Ephemeris.inclined_circular() opt-in backend (M-3D Phase 2)
```

### Task 2.2 — un-skip the Phase-1 inclined `synodic_omega_vec` test + lint/type

**Files:** `tests/core/test_frames_omega_vec.py` (remove the Task 1.2 skip on
`test_synodic_omega_vec_tilts_for_inclined_backend`).

- [ ] Un-skip; run `uv run pytest tests/core/test_frames_omega_vec.py -q` green.
- [ ] `uv run ruff check . && uv run ruff format --check . && uv run mypy src tests`.

Commit:

```
core/ephemeris: wire inclined backend into synodic_omega_vec tilt gate (M-3D Phase 2)
```

---

## Phase 3 — Flyby bend-decomposition diagnostic (`core/flyby.py`)

Diagnostic only — no optimiser DOF (Approval Q4).

### Task 3.0 — `bend_decompose` pure helper

**Files:** `src/cyclerfinder/core/flyby.py`; test
`tests/core/test_flyby_bend_decompose.py`.

`bend_decompose(vin_vec, vout_vec, orbit_normal) -> tuple[float, float]` returns
`(delta_inplane_rad, delta_outofplane_rad)`: project both V∞ vectors onto the
plane orthogonal to `orbit_normal` and take the in-plane angle between the
projections (`delta_inplane`); the out-of-plane part is the angle change along
`orbit_normal` (the difference of the signed `arcsin` of each vector's
normal-component over its magnitude). Pure; raises on zero-length inputs like
`bend_angle` (`flyby.py:138-139`).

#### Failing test

```python
"""M-3D Phase 3: flyby bend in/out-of-plane attribution (plan §3). DIAGNOSTIC
ONLY. Physics invariants (no sourced V_inf anchor): a purely in-plane bend has
zero out-of-plane component, and the total recovers bend_angle.
"""
from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.flyby import bend_angle, bend_decompose

_NZ = np.array([0.0, 0.0, 1.0])  # ecliptic normal


def test_inplane_bend_has_zero_outofplane() -> None:
    vin = np.array([5.0, 0.0, 0.0])
    vout = np.array([5.0 * np.cos(0.3), 5.0 * np.sin(0.3), 0.0])
    d_in, d_out = bend_decompose(vin, vout, _NZ)
    assert d_in == pytest.approx(0.3, abs=1e-9)
    assert abs(d_out) < 1e-9


def test_pure_outofplane_bend_has_zero_inplane() -> None:
    vin = np.array([5.0, 0.0, 0.0])
    vout = np.array([5.0 * np.cos(0.2), 0.0, 5.0 * np.sin(0.2)])
    d_in, d_out = bend_decompose(vin, vout, _NZ)
    assert abs(d_in) < 1e-9
    assert abs(d_out) == pytest.approx(0.2, abs=1e-9)


def test_decomposition_components_consistent_with_total_bend() -> None:
    vin = np.array([4.0, 1.0, 0.5])
    vout = np.array([3.5, 1.8, -0.4])
    d_in, d_out = bend_decompose(vin, vout, _NZ)
    total = bend_angle(vin, vout)
    # In/out-of-plane are orthogonal contributions; quadrature recovers the
    # total to small-angle order (consistency invariant, not a sourced anchor).
    assert np.hypot(d_in, d_out) == pytest.approx(total, rel=0.15)


def test_zero_vector_raises() -> None:
    with pytest.raises(ValueError):
        bend_decompose(np.zeros(3), np.array([1.0, 0.0, 0.0]), _NZ)
```

Run → **red** → impl `bend_decompose` → **green**. Commit:

```
core/flyby: bend_decompose in/out-of-plane diagnostic (M-3D Phase 3)
```

### Task 3.1 — Phase-3 lint/type gate

- [ ] `uv run pytest tests/core/test_flyby_bend_decompose.py tests/core/test_flyby*.py -q`
- [ ] `uv run ruff check . && uv run ruff format --check . && uv run mypy src tests`.

Commit if cleanup:

```
core/flyby: Phase 3 lint/type gate (M-3D)
```

---

## Phase 4 — `circular-inclined` fidelity rung

Widen the `Fidelity` enum in **all four** sites and wire the third resolvable
rung into `solve_at_fidelity`. Additive only — `fidelity_persistence` untouched.

### Task 4.0 — widen the `Fidelity` literal + runtime frozenset

**Files:** `src/cyclerfinder/data/provenance.py`; test
`tests/data/test_fidelity_widen.py`.

Insert `"circular-inclined"` into the `Literal` (`provenance.py:95`), the
`_FIDELITIES` frozenset (`:109`), and update the ladder-order docstring
(`:14-15`, `:95-103`) to read coplanar → **inclined** → analytic-ephemeris →
real-DE440.

#### Failing test

```python
"""M-3D Phase 4: circular-inclined Fidelity rung registered (plan §4)."""
from __future__ import annotations

from cyclerfinder.data.provenance import is_fidelity


def test_circular_inclined_is_a_known_fidelity() -> None:
    assert is_fidelity("circular-inclined")


def test_existing_fidelities_unchanged() -> None:
    for tier in ("circular-coplanar", "analytic-ephemeris", "real-de440"):
        assert is_fidelity(tier)
    assert not is_fidelity("bogus-rung")
```

Run → **red** → impl the widening → **green**. Commit:

```
data/provenance: register circular-inclined Fidelity rung (M-3D Phase 4)
```

### Task 4.1 — widen the two schema enums

**Files:** `src/cyclerfinder/data/catalogue.schema.json` (the `orbit_fidelity`
enum `:59` and `vinf_fidelity` enum `:64`); test reuse
`tests/data/test_schema_v45_fields.py` / a new
`tests/data/test_schema_circular_inclined.py`.

Add `"circular-inclined"` to both enum arrays (keeping `null`). Bump the schema
description note to record the additive Fidelity-tier widening (do **not** bump
the schema `$id`/version unless the project's schema-versioning convention
requires it — check `test_schema_v45_fields.py` for the asserted version string
first; if a bump is required, do it in this same commit per the schema-change
discipline).

#### Failing test

```python
"""M-3D Phase 4: schema accepts circular-inclined fidelity tags (plan §4)."""
from __future__ import annotations

import json
from pathlib import Path

SCHEMA = json.loads(Path("data/catalogue.schema.json").read_text())


def _enum(field: str) -> list:
    return SCHEMA["items"]["properties"][field]["enum"]


def test_orbit_and_vinf_fidelity_enums_include_circular_inclined() -> None:
    assert "circular-inclined" in _enum("orbit_fidelity")
    assert "circular-inclined" in _enum("vinf_fidelity")
    # existing tiers preserved
    for tier in ("circular-coplanar", "analytic-ephemeris", "real-de440", None):
        assert tier in _enum("orbit_fidelity")
```

> **Verify the JSON path** (`items.properties.<field>.enum`) against the live
> file before writing the test — adjust the accessor to the actual schema shape
> if it nests differently.

Run → **red** → widen both enums → **green**; then
`uv run pytest tests/data/test_schema_*.py -q` (no existing schema test
regresses; the catalogue still validates). Commit:

```
data/schema: accept circular-inclined orbit/vinf fidelity tags (M-3D Phase 4)
```

### Task 4.2 — `solve_at_fidelity` resolves `circular-inclined`

**Files:** `src/cyclerfinder/verify/fidelity.py`; test
`tests/verify/test_fidelity_inclined_rung.py`.

Add a `circular-inclined` branch to `solve_at_fidelity`
(`verify/fidelity.py:165-182`) dispatching to a new `_solve_inclined(cell, *,
a_au, e)`. The inclined rung's job is a **closed-form inclined-circular
state/Tisserand screen**, not a Lambert closure (closure is M-ED's). Concretely
`_solve_inclined` builds `Ephemeris.inclined_circular()` and returns a
`FidelitySolution(fidelity="circular-inclined", ...)` whose `outbound_tof_days`
and `vinf_kms` come from the **same coplanar resonance construction** as
`_solve_coplanar` (the inc-only lift does not change the resonance `(a,e)`
construction's scalar ToF/V∞; it changes only the *frame geometry*, which this
rung reports via the inclined ephemeris being available to downstream
diagnostics). `converged=True` (closed-form). **Document in the docstring:**
`circular-inclined` = real i/Ω, mean sma, **e ignored** (circular); it is a
pre-filter/diagnostic, not a closure engine.

> **HONESTY NOTE (binding):** there is no sourced "inclined cycler closes at X"
> anchor, so this rung must NOT assert any self-computed closure as a golden. The
> test below asserts only (a) the rung resolves without raising, and (b) its
> scalar ToF/V∞ equal the coplanar rung's (the inc-only lift is frame geometry,
> not a scalar-magnitude change at the resonance-construction level). Both are
> regression/consistency invariants, not sourced goldens.

#### Failing test

```python
"""M-3D Phase 4: solve_at_fidelity resolves circular-inclined (plan §4).

HONEST SCOPE: no sourced inclined-closure anchor exists. Asserted: the rung
resolves (no FidelityRungUnavailableError) and its scalar ToF/V_inf equal the
coplanar rung's (inc-only lift is frame geometry; consistency invariant).
"""
from __future__ import annotations

import pytest

from cyclerfinder.search.sequence import Cell
from cyclerfinder.verify.fidelity import solve_at_fidelity


def _aldrin_cell() -> Cell:
    return Cell(
        bodies=("E", "M"),
        sequence=("E", "M", "E"),
        period_k=2,
        per_leg_revs=(0, 0),
        per_leg_branch=("single", "single"),
    )


def test_inclined_rung_resolves_and_matches_coplanar_scalars() -> None:
    cell = _aldrin_cell()
    # sourced (a, e) for the Aldrin E-M cycler resonance construction
    a_au, e = 1.6004, 0.3929  # read from the catalogue row in the real fixture
    coplanar = solve_at_fidelity(cell, "circular-coplanar", a_au=a_au, e=e)
    inclined = solve_at_fidelity(cell, "circular-inclined", a_au=a_au, e=e)
    assert inclined.fidelity == "circular-inclined"
    assert inclined.converged is True
    assert inclined.outbound_tof_days == pytest.approx(coplanar.outbound_tof_days)
    assert inclined.vinf_kms == pytest.approx(coplanar.vinf_kms)
```

> **Sourced (a,e) caveat:** replace the inline `a_au, e` with the values read
> from the Aldrin catalogue row at test time (mirror how
> `tests/verify/test_fidelity_gate.py` sources its anchors from the row), so no
> magic numbers are pinned. The construction inputs must trace to the catalogue.

Run → **red** → impl the branch + `_solve_inclined` → **green**. Commit:

```
verify/fidelity: circular-inclined rung in solve_at_fidelity (additive) (M-3D Phase 4)
```

### Task 4.3 — Phase-4 full gate (provenance + schema + fidelity) + lint/type

- [ ] `uv run pytest tests/data/test_fidelity_widen.py tests/data/test_schema_*.py tests/verify/test_fidelity_inclined_rung.py tests/verify/test_fidelity_gate.py -q`
- [ ] Confirm the catalogue still validates against the widened schema
  (`uv run pytest tests/data/test_jsonschema.py -q`).
- [ ] `uv run ruff check . && uv run ruff format --check . && uv run mypy src tests`.

Commit if cleanup:

```
verify/fidelity: Phase 4 gate — circular-inclined rung end-to-end (M-3D)
```

---

## Phase 5 — Tisserand consistency pins + honest-limit docs (`search/sequence.py`)

`tisserand_feasible` **already** consults `linkable_3d` behind `ephem`
(`sequence.py:329,356`). M-3D pins the reviewed coplanar→3-D baseline diff the
M8-Core tests await and documents the node-compatibility limit.

### Task 5.0 — pin the inclined `tisserand_feasible` diff (reviewed, not silent)

**Files:** `tests/search/test_sequence_multibody.py` (extend; the header
`:1-7` already declares it pins the coplanar baseline "so the M-3D inclination
lift is a reviewed diff").

Add parametrised assertions that for a body pair where coplanar `linkable` is
False but inclined linkability opens (per `linkable_3d` monotonicity,
`tisserand.py:548-566`), `tisserand_feasible(cell, vinf_cap, ephem=...)` returns
the documented 3-D result, while `tisserand_feasible(cell, vinf_cap)` (no ephem)
keeps the **byte-identical coplanar** result. The EXPECTED side is the
**predicate's own logical contract** (coplanar-True ⇒ 3-D-True monotonicity, a
physics invariant), never a sourced feasibility number.

#### Failing test (sketch — finalise pair against the live grid)

```python
def test_tisserand_3d_is_superset_of_coplanar() -> None:
    """ephem-supplied tisserand_feasible never rejects a coplanar-feasible cell
    (monotonicity, tisserand.py:548-566) — a reviewed 3-D diff, not a silent
    change. Physics invariant, not a sourced anchor."""
    from cyclerfinder.core.ephemeris import Ephemeris

    ephem = Ephemeris.inclined_circular()
    for cell in _SAMPLE_CELLS:  # the existing module's pinned cells
        coplanar = tisserand_feasible(cell, vinf_cap=9.0)
        threed = tisserand_feasible(cell, vinf_cap=9.0, ephem=ephem)
        if coplanar:
            assert threed  # coplanar-True => 3-D-True (monotone superset)
```

> Reuse `_SAMPLE_CELLS` / the cells the module already pins. The point is the
> **reviewed diff**: any cell whose 3-D verdict differs from coplanar is captured
> here with the monotonicity reason, so the inclination lift is auditable.

Run → **red** (if the inclined-ephem path is newly exercised) → adjust the
assertion to the live verdicts (documenting each flip's physical reason in a
comment) → **green**. Commit:

```
test: pin reviewed coplanar->3D tisserand_feasible diff (monotone superset) (M-3D Phase 5)
```

### Task 5.1 — honest-limit docstring at the gate + Phase-5 lint/type

**Files:** `src/cyclerfinder/search/sequence.py` (extend the
`tisserand_feasible` docstring, `:300-327`).

Add one paragraph: `linkable_3d` tests **energetic** linkability with an
inclination budget (`i_sc_max_deg=30°` default, `tisserand.py:476` — a *modelling
choice*, not sourced), **not node compatibility**; two bodies can share a
reachable `i_sc` magnitude while their lines of nodes make the transfer
phasing-infeasible. The predicate is a *necessary, not sufficient* screen — same
character as coplanar `linkable` (design §5). Pure docstring; no behaviour change.

- [ ] `uv run pytest tests/search/test_sequence_multibody.py -q`
- [ ] `uv run pytest -m "not slow"` green (full fast suite — no regression).
- [ ] `uv run ruff check . && uv run ruff format --check . && uv run mypy src tests`.

Commit:

```
search/sequence: document linkable_3d node-compatibility limit at the gate (M-3D Phase 5)
```

---

## Out of scope (explicit carve-outs, design §7 non-goals + Approval)

- **Inclined idealized closure attempt** — rewiring `closure_residual` onto the
  vector-ω⃗ frame and re-baselining its thresholds is **M-ED's** (Approval Q4).
  M-3D ships the *primitive* (Phase 1) and proves its coplanar limit; it does
  **not** flip `closure_residual` to the vector form.
- **Optimiser b-plane / out-of-plane steering DOF** — no `search/optimize.py`
  change (Approval Q4). M-3D ships `bend_decompose` (diagnostic) only.
- **Dynamic-frame inclined (orbit-normal) variant** of `to_rotating_dynamic`
  — **M-ED's** (Approval Q5). M-3D leaves `to_rotating_dynamic`
  (`frames.py:323-389`) untouched.
- **Eccentricity in the analytic rung** — separable follow-on (promotes
  `circular-inclined` → `analytic-ephemeris`); `_circular_inplane_state` keeps
  ignoring `PLANETS.ecc` (Approval Q3).
- **Mutating live `PLANETS` to inclined** — rejected (design §7 Approach 3); the
  opt-in factory (`inclined_planets()`) is the only path.
- **Planet-centric / CR3BP / T-P graph** — task #76 (design §5 fence).
- **Real-ephemeris blind discovery / TCM budget / `mode="ballistic"` corrector /
  descriptor parser / seed ladder / discover V3 gate** — **M-ED**, landing
  concurrently; M-3D does not touch those files.

---

## Definition of done

0. `core/frames.py` — `to_rotating_omega_vec`/`from_rotating_omega_vec` are a
   bit-exact superset of the scalar forms (`numpy.array_equal` on every M3
   golden input) with a general-axis round-trip identity; `synodic_omega_vec`
   returns the anchor orbit-normal angular velocity (z-only for coplanar
   circular, tilted for inclined). The M3 frame goldens pass **unchanged**. **Green.**
1. `core/ephemeris.py` — `inclined_planets()` + `Ephemeris.inclined_circular()`
   give sourced inclined states **without mutating `PLANETS`** (live coplanar
   default byte-identical). **Green.**
2. `core/flyby.py` — `bend_decompose` attributes a flyby bend to in/out-of-plane
   components (diagnostic only; no optimiser DOF). **Green.**
3. `Fidelity` is widened to `circular-inclined` in **all four** sites (Literal,
   `_FIDELITIES`, both schema enums); `solve_at_fidelity` resolves the rung
   additively (no `fidelity_persistence` change); the catalogue still validates.
   **Green.**
4. `tisserand_feasible`'s coplanar→3-D diff is **pinned as a reviewed monotone
   superset** and the node-compatibility limit is documented at the gate. **Green.**
5. No assertion anywhere uses a closure/V∞ value M-3D computed as the EXPECTED
   side — only the sourced inc/Ω (input-is-source), physics invariants
   (coplanar limit, round-trip, monotonicity, plane-decomposition consistency),
   and refactor-equivalence pins.
6. `uv run pytest -m "not slow"`, ruff, ruff format, mypy all clean.

---

## Self-Review

### Design coverage (every approved decision → plan task)

| Design / Approval item | Plan location |
|---|---|
| Q1 F3 vector-ω⃗ primitive, F1 interim acceptance | Phase 1 (Tasks 1.0-1.2); coplanar-limit golden gate verbatim §Phase 1 |
| Q1 `ω⃗ = ω ẑ` reproduces M3 goldens bit-for-bit | Task 1.0/1.1 (`array_equal`), Task 1.3 (goldens unchanged) |
| §2 anchor `ω⃗` on home orbit normal `(r×v)/|r|²` | Task 1.2 `synodic_omega_vec` |
| Activate inclined backend behind explicit opt-in, never mutate PLANETS | Phase 2 (Tasks 2.0-2.1) |
| Q2 R-A new `circular-inclined` rung (literal/schema widening) | Phase 4 (Tasks 4.0-4.2) |
| Q2 wire into `solve_at_fidelity` as third rung, additive | Task 4.2 |
| Q3 inclination-only, eccentricity separable (e ignored) | §Architecture, Task 4.2 docstring, Out of scope |
| Q4 bend-decomposition diagnostic only, no optimiser DOF | Phase 3 (Task 3.0); Out of scope |
| Q5 dynamic-frame inclined variant is M-ED's | Out of scope (frames.py:323-389 untouched) |
| Q6 keep 30° `i_sc_max_deg` default as tunable | Task 5.1 docstring |
| §5 `linkable_3d` consistency in `tisserand_feasible` (already wired) | Phase 5 (pin + document, not rewire) |
| §6 pinned-baseline reviewed diff | Task 5.0 |
| Honest test strategy (physics invariants carry the load) | Honesty boundary; every test docstring |
| Phasing: closure-frame primitive first, independently shippable | Phasing table |

### Placeholder scan
No `TODO`/`FIXME`/`...`-as-impl in production instructions. Every file path is
concrete. Commit messages are exact. The two inline magic numbers flagged
inline (the Aldrin `(a,e)` in Task 4.2; the `_SAMPLE_CELLS` in Task 5.0) carry
explicit "source from the catalogue row / reuse the module's pinned cells"
instructions so no self-computed value is pinned.

### Type consistency
`to_rotating_omega_vec`/`from_rotating_omega_vec` take `omega_vec: Vec3` and
return `tuple[Vec3, Vec3]` (matching the existing scalar shape,
`frames.py:90-146`). `synodic_omega_vec(anchor_body: str, t_sec: float,
ephem: Ephemeris) -> Vec3`. `inclined_planets() -> dict[str, PlanetData]`
(matches `_InclinedCircularBackend.__init__` injection signature,
`ephemeris.py:129`). `bend_decompose(Vec3, Vec3, Vec3) -> tuple[float, float]`
(mirrors `bend_angle`'s `Vec3 -> float`, `flyby.py:117`). The `Fidelity` widening
keeps the `Literal`↔`_FIDELITIES` parity that `is_fidelity` depends on
(`provenance.py:109,117-119`). `FidelitySolution.fidelity` accepts the new tier
(it is typed `Fidelity`, `fidelity.py:104`).

### Concurrency safety (M-ED overlap)
M-3D's write-set (`core/frames.py`, `core/ephemeris.py`, `core/flyby.py`,
`data/provenance.py`, `data/catalogue.schema.json`, `verify/fidelity.py`,
`search/sequence.py` docstring + `tests/search/test_sequence_multibody.py`) is
**disjoint** from M-ED's (`search/correct.py`, `search/descriptor.py`,
`search/seed_ladder.py`, the `mode="ballistic"` branch of `search/optimize.py`,
`data/discover.py`). The single shared file risk is `data/provenance.py` — M-ED
does not widen `Fidelity` (it consumes `real-de440`), so the widening is
non-conflicting, but the implementer must `git pull`/rebase before the Phase-4
provenance edit and re-run the schema tests.

### Unverifiable / shifted claims (flagged, not papered over)
1. **`PLANETS.ecc` became live and non-zero** since the design was written
   (`constants.py:180,197,208`). It is **ignored** by the circular backend
   (`ephemeris.py:72-90`), so the inc-only scope holds; flagged in Design
   references and the rung docstring (Task 4.2: "e ignored").
2. **`tisserand_feasible` moved** to `sequence.py:276-366` (design said
   `:300-363`) and **already** threads `linkable_3d` — M-3D's Phase 5 is
   **pin + document**, not the wiring the design's §1 audit table implied was
   pending. Flagged.
3. **The Fidelity enum lives in four sites**, not just the `Literal` — the
   `_FIDELITIES` frozenset and the two `catalogue.schema.json` enum arrays
   (`:59,:64`) must all be widened or `is_fidelity`/schema validation diverge.
   Enumerated in Task 4.0/4.1.
4. **Schema version bump** (Task 4.1) is left conditional: the implementer must
   read the asserted version string in `tests/data/test_schema_v45_fields.py`
   first and bump only if the convention requires it — I did not pre-decide the
   `$id` change to avoid breaking a version-pin test blind.
5. **The `bend_decompose` quadrature-consistency tolerance** (Task 3.0,
   `rel=0.15`) is a small-angle approximation, not exact — in/out-of-plane are
   orthogonal only to first order. The two pure-component tests
   (`abs=1e-9`) are the exact gates; the combined test is a loose sanity check.
   Flagged so a reviewer does not read it as an exactness claim.
