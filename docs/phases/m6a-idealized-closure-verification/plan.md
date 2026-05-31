# M6a — Idealized closure verification

**Spec reference:** spec.md §4 (architecture — `verify/propagate.py`, `verify/crosscheck.py`), §5 step 7 (verify: multi-lap propagation; periodicity / re-encounter miss), §6 (top-level interfaces — `find_cyclers`, `Cycler`), §8 (M6 milestone definition: "ephemeris + verify — astropy backend, multi-lap propagate, crosscheck. Gate: best E–M cycler verified periodic over ≥3 laps."), §9 (validation anchors, in particular the published 5.65 / 3.05 km/s 2-synodic E-M V∞), §10 (closure-frame correctness risk — "unit-test the rotating-frame transform; a wrong frame silently fakes/breaks closure"), §12(a) (idealized vs ephemeris-mode optimisation — M6a covers the verification half, M6b covers the TCM-minimisation half), **§12(c) (dynamic ephemeris frame + tolerant verification — the headline architectural decision for M6a)**, §12.1 (idealized→ephemeris phase-matching bridge — already half-shipped via the M6 slice `9b2611d`), §14 (V0–V5 validation gauntlet — M6a stands up the V2 multi-lap periodicity gate).

**Purpose:** stand up the **multi-lap propagation + tolerant closure-drift verification machinery** that turns an `OptimisationResult` from M5 into a checkable spec §14 V2 statement: "this cycler stays within bounded drift in the dynamic rotating frame over ≥3 laps." M0–M5 produced the idealized geometry; the M6 slice (commit `9b2611d`) brought the astropy `Ephemeris` backend and `phase_match.find_real_windows` forward; M6a closes the remaining verification half of M6 by adding the dynamic (non-uniform) rotating-frame transform and the `verify/` subpackage with `propagate.py`. M6b inherits the locked `StabilityReport` shape and fills the per-lap TCM-ΔV fields that M6a leaves as placeholders.

**Gate (definition of done):**

1. `tests/verify/test_propagate.py::test_2syn_em_cycler_periodic_over_3_laps` (**M6a BINDING GATE — spec §8 M6 anchor**) asserts the best 2-synodic E-M cycler — sourced from M5's `find_cyclers(("E","M"), k_synodic=2, vinf_cap=7.0, seed=0)` top result, or from the catalogue entry `s1l1-2syn-em-cpom` if M5's result is not available at test time — when fed to `verify_long_term_stability(cycler, n_laps=3, ephem=Ephemeris("astropy"))` returns `report.stable == True` with `report.max_drift_km < DRIFT_TOLERANCE_KM` (50,000 km) and `report.n_laps_propagated == 3`. **The cycler is the fixture, not the tolerance** — no cherry-picking, no per-test tolerance widening.
2. `tests/core/test_frames_dynamic.py::test_dynamic_frame_roundtrip_identity` asserts `from_rotating_dynamic(to_rotating_dynamic(r, v, t_sec, bodies, ephem))` recovers `(r, v)` to relative tolerance ≤ 1e-10 across a 5-year time grid at 1-AU position scale and 30-km/s velocity scale. This is the **spec §10 binding unit test** — a wrong rotating-frame transform silently fakes or breaks every drift measurement downstream.
3. `tests/core/test_frames_dynamic.py::test_dynamic_frame_reduces_to_uniform_for_circular` asserts that when `ephem = Ephemeris("circular")` and `bodies = ("E", "M")`, the dynamic frame's drift over one synodic period agrees with `to_rotating(...)` at Earth's mean motion to ≤ 1e-8 km on a 1-AU test state. Circular ephemerides should produce a frame indistinguishable from M3's uniform frame; any divergence here means the dynamic frame's `omega(t)` computation is wrong.
4. `tests/verify/test_propagate.py::test_aldrin_cycler_periodic_over_3_laps_circular` asserts the Aldrin cycler (built from `build_aldrin_seed(Ephemeris("circular"))`, 1-synodic) when fed to `verify_long_term_stability(..., ephem=Ephemeris("circular"))` returns `stable=True` with `max_drift_km < 1.0` (the circular-coplanar case is essentially zero-drift; this catches frame-transform regressions in the circular regime where the answer is known exactly).
5. `tests/verify/test_propagate.py::test_stability_report_frozen_and_fields_locked` asserts `StabilityReport` raises `FrozenInstanceError` on attribute assignment, that `per_lap_dv == (0.0,) * n_laps` and `total_tcm_dv == 0.0` on every M6a-produced report (the M6b TCM placeholders), and that `frame_used in ("dynamic", "uniform")`.
6. `tests/verify/test_propagate.py::test_propagate_lap_matches_construct_at_encounters` asserts that the sampled trajectory from `propagate_lap(cycler, ephem, t_start, t_start + cycler.period, n_samples=100)` passes through every `Encounter.r` at the corresponding `Encounter.t` to ≤ 10 km (the universal-variable propagator's float-noise floor for AU-scale propagation).
7. `tests/verify/test_propagate.py::test_lap_to_lap_drift_zero_for_circular_aldrin` asserts that `lap_to_lap_drift(samples_lap_0, samples_lap_2)` on a circular-coplanar Aldrin cycler is below 100 km — the M3 circular case has no eccentricity-driven breathing, so drift must be at the numerical-noise floor.
8. `uv run pytest`, `uv run ruff check .`, `uv run ruff format --check .`, `uv run mypy src tests` all clean. CI green on the M6a commit.

The 50,000 km drift tolerance is the **second binding architectural choice** of M6a (after the dynamic-frame decision) and is justified in §4.3 below. It is loose enough to absorb real eccentricity-driven geometric breathing on Earth-Mars trajectories (~0.02° of geometric breathing per lap at Mars's mean radius of ~1.5 AU ≈ 2.6 million km / degree → ~50,000 km / 0.02°) and tight enough to reject any propagator regression, frame error, or genuinely non-periodic candidate. See §4.3 for the derivation.

---

## 1. What this milestone delivers

Three new pieces of source code and their test files; one new subpackage. M6a is **additive** for everything except `core/frames.py`, which receives a **strictly additive extension** — every M3 function (`to_rotating`, `from_rotating`, `synodic_omega`) keeps its current signature and semantics so M3/M4/M5 consumers remain unchanged.

### 1.1 `src/cyclerfinder/core/frames.py` (EXTEND — not replace)

New public surface, in dependency order:

- `synodic_omega_dynamic(t_sec, bodies, ephem) -> float` — instantaneous synodic angular rate anchored to current Sun-body positions per spec §12(c). For `len(bodies) == 2` returns the instantaneous angular rate of the line-of-bodies relative to the inertial frame. **Non-uniform** — the value depends on `t_sec` because real orbits are eccentric.
- `to_rotating_dynamic(r_inertial, v_inertial, t_sec, bodies, ephem) -> tuple[Vec3, Vec3]` — non-uniform rotation per spec §12(c). The frame's x-axis at time `t_sec` is aligned with the instantaneous Sun→body[0] direction; the angular rate is `synodic_omega_dynamic(t_sec, bodies, ephem)`.
- `from_rotating_dynamic(r_rot, v_rot, t_sec, bodies, ephem) -> tuple[Vec3, Vec3]` — inverse of `to_rotating_dynamic`. Round-trip identity to ≥ 1e-10 relative (the spec §10 binding test).

**Unchanged surface (M3/M5 consumers must remain bit-identical):** `to_rotating`, `from_rotating`, `synodic_omega`. No new behaviour, no signature changes, no deprecation. The dynamic frame is **additive**, not replacement.

### 1.2 `src/cyclerfinder/verify/__init__.py` (NEW)

First occupant of the `verify/` subpackage. Empty re-export module per the M3-precedent — callers import from the specific submodule (`from cyclerfinder.verify.propagate import verify_long_term_stability`). The `__init__.py` exists primarily to make the subpackage importable under `mypy --strict` and to anchor `tests/verify/` next to it.

### 1.3 `src/cyclerfinder/verify/propagate.py` (NEW)

Public surface, in dependency order:

- `DRIFT_TOLERANCE_KM: Final[float] = 50_000.0` — the binding tolerance per §4.3.
- `StabilityReport` — frozen dataclass per §3.3 below; the locked spec §12 verification interface.
- `propagate_lap(cycler, ephem, t_start, t_end, n_samples) -> NDArray[np.float64]` — propagate the cycler's trajectory across one lap, sampling at `n_samples` uniformly-spaced points. Uses `core.kepler.propagate` between encounters; respects the cycler's `legs[]` structure (each `Leg` is a 2-body heliocentric arc that the universal-variable propagator handles natively).
- `multi_lap_propagation(cycler, ephem, n_laps, t_start, n_samples_per_lap=100) -> dict[str, NDArray[np.float64]]` — propagate continuously through `n_laps` consecutive laps. Returns a dict with the full sampled trajectory split per lap (key `"lap_{i}_samples"`) plus per-lap statistics ready for the drift check.
- `lap_to_lap_drift(samples_lap_0, samples_lap_n, t_lap_0_start, t_lap_n_start, bodies, ephem) -> float` — max km drift between corresponding rotating-frame positions, computed in the **dynamic** rotating frame anchored to `bodies`. Per spec §12(c) the verification target is **bounded** drift, not zero.
- `verify_long_term_stability(cycler, n_laps, ephem, *, t_start=0.0, frame_bodies=None) -> StabilityReport` — the locked spec §12 entry point. Implements `stable` (a bool derived from `max_drift_km < DRIFT_TOLERANCE_KM`), `max_drift_km`, `n_laps_propagated`, `per_lap_drift_km`, `max_drift_lap_index`. The TCM fields `per_lap_dv` and `total_tcm_dv` are populated with zeros — M6b fills them in.

Module-internal helpers (private):

- `_resolve_frame_bodies(cycler, frame_bodies)` — picks `frame_bodies` from the cycler's `bodies` attribute if not supplied; defaults to `(cycler.bodies[0], cycler.bodies[-1])` or `(cycler.bodies[0], cycler.bodies[1])` per §3.4.
- `_propagate_leg(r0, v0, t0, t1, ephem, n_samples)` — sample one heliocentric leg via `core.kepler.propagate`; returns `(N, 6)` array of `[x, y, z, vx, vy, vz]` per sample.
- `_dynamic_frame_drift(sample, sample_ref, t_sample, t_ref, bodies, ephem)` — per-point drift computation that calls `to_rotating_dynamic` on both states and returns `||r_sample_rot - r_ref_rot||`.

### 1.4 Test files

- `tests/core/test_frames_dynamic.py` — round-trip identity, circular-degenerate-to-uniform, basic numeric properties of `synodic_omega_dynamic`. (Lives in `tests/core/` so it sits next to the existing M3 `test_frames.py`.)
- `tests/verify/__init__.py` — new empty test package for the `verify/` subpackage.
- `tests/verify/test_propagate.py` — gate tests in §4.1, helper tests in §4.2, integration tests in §4.3.

### 1.5 Explicit non-goals (M6a boundaries)

These belong to M6b, M7, or M8; **do not stub or partially implement them in M6a** beyond what the locked dataclass shape requires:

| Out of M6a | Where it lands | Why deferred |
|---|---|---|
| Real-ephemeris TCM-budget computation (per-lap ΔV, total horizon TCM) | **M6b** | The `StabilityReport.per_lap_dv` and `total_tcm_dv` fields exist in M6a as locked zeros; M6b fills the bodies. Computing them requires the M6b optimisation slice (`optimise_cell_ephemeris`) which is not M6a's deliverable. |
| Body of `optimise_cell_ephemeris` | **M6b** | M5 raises `NotImplementedError("requires M6 ephemeris backend")`; M6a leaves this raise in place. M6b fills the body using `verify_long_term_stability` as its drift-feasibility check inside the optimiser. |
| Real launch-window TCM-budget columns on `/launch-windows/` | **M6b** | Spec §12.1 / the existing M6 slice produced *geometric* launch windows (Lambert mismatch); the TCM-cost columns are M6b's contribution. |
| Catalogue ingest of M6a stability results | **M7** | M7's signature matcher + catalogue writer consume `StabilityReport`; M6a only produces it. |
| Lambert cross-check (`verify/crosscheck.py`) | **M7 or stretch** | Spec §14 V1 is a separate gate; M6a is the V2 gate. The cross-check is single-leg, doesn't share infrastructure with multi-lap propagation. |
| GMAT bridge (`verify/gmat_bridge.py`) | **Stretch** | Spec §14 V4, explicit stretch goal per spec §7. |
| Phase-match epoch-picking integration | **Already shipped in M6 slice** | `phase_match.find_real_windows` exists; M6a's `verify_long_term_stability` accepts a `t_start` so callers can pass an epoch picked by `find_real_windows`, but M6a does not call `find_real_windows` itself. |
| Multi-body dynamic frame for VEM (3+ bodies) | **M8** | `to_rotating_dynamic` accepts `bodies: tuple[str, ...]` of length ≥ 2; the 2-body case is binding for M6a's E-M gate. The 3+-body case is documented but tests defer to M8's VEM campaign. |
| Catalogue-level V2 validation block population | **M7** | The `validation.gates.V2.max_drift_km` field on the spec §16.1 catalogue record is populated by M7's batch-validate runner consuming `StabilityReport`. M6a produces the report; M7 writes it into YAML. |

---

## 2. File tree after M6a

```
cyclers/
├── … (M0/M1/M2/M3/M4/M5 layout preserved unchanged)
├── src/cyclerfinder/
│   ├── core/
│   │   ├── constants.py                  # M0 — unchanged
│   │   ├── ephemeris.py                  # M1 + M6 slice — unchanged
│   │   ├── lambert.py                    # M1 — unchanged
│   │   ├── kepler.py                     # M1 — unchanged
│   │   ├── flyby.py                      # M2 — unchanged
│   │   └── frames.py                     # M3 + EXTENDED in M6a (additive only)
│   ├── search/                           # all unchanged
│   ├── model/                            # all unchanged
│   └── verify/                           # NEW subpackage (M6a)
│       ├── __init__.py                   # NEW
│       └── propagate.py                  # NEW
└── tests/
    ├── core/
    │   ├── test_frames.py                # M3 — unchanged
    │   └── test_frames_dynamic.py        # NEW (M6a)
    └── verify/                           # NEW
        ├── __init__.py                   # NEW
        └── test_propagate.py             # NEW (M6a — includes the gate tests)
```

Subpackages `data/` and `viz/` still remain uncreated — M7/M8 territory. The only edit to a pre-existing source file is `core/frames.py`, which receives new public functions but **no edits** to existing function bodies, signatures, or docstrings (per spec §10 risk avoidance: an M3 consumer must keep getting bit-identical results from `to_rotating`/`from_rotating`/`synodic_omega`).

---

## 3. Module designs

This section walks the four pieces of M6a in dependency order: (3.1) the dynamic-frame extension to `core/frames.py`; (3.2) the per-lap and multi-lap propagators in `verify/propagate.py`; (3.3) the `StabilityReport` dataclass + `verify_long_term_stability` API; (3.4) the cycler→frame-bodies resolution helper.

### 3.1 Dynamic rotating frame (`core/frames.py` extension) — spec §12(c) binding

Spec §12(c) supersedes the §4 `frames.py` original wording. The verification frame is a **non-uniform rotating frame anchored to instantaneous Sun-body positions**, not a constant angular velocity. The frame's x-axis at time `t_sec` tracks the instantaneous Sun→body[0] direction; the angular rate is computed from the instantaneous geometry, not pre-computed at `t=0`.

**Binding decisions baked into M6a:**

- **The dynamic frame is additive, not replacement.** M3's uniform frame stays unchanged for M3/M5 consumers (closure residual on circular-coplanar cyclers, M5's idealized-mode optimiser). The dynamic frame is the verification frame for real-ephemeris work — explicitly the M6a/M6b/M7 path. Two frames exist; the caller picks via the function name (`to_rotating` vs `to_rotating_dynamic`).
- **Frame anchor:** `bodies[0]` defines the x-axis direction (Sun → body[0] is the +x direction in the rotating frame at every `t`). `bodies[1]` is needed only for the synodic angular rate computation; for `len(bodies) >= 3` the rate is computed from the **first pair** `(bodies[0], bodies[1])` and the additional bodies are tracked as "context" but do not influence the frame definition. This matches the spec §12(c) "Sun-Earth (or Earth-Mars)" wording exactly; the 3+-body extension belongs to M8.
- **Angular rate definition:** the instantaneous synodic rate is the rate at which the body[0]→body[1] line is rotating in the inertial frame, computed from the cross product `(r_b1 - r_b0) × (v_b1 - v_b0) / ||r_b1 - r_b0||²`. The sign convention matches the prograde direction (positive ω for prograde planets in the J2000 ecliptic).
- **Velocity transformation:** the Coriolis-style velocity correction `v - ω(t) × r` uses the instantaneous ω(t), not the M3 constant ω. This is the only mathematical difference from M3's uniform-frame transform; the rotation-matrix step is identical (passive rotation by the integrated frame angle `θ(t) = atan2(r_b0_y, r_b0_x)` at time `t`).
- **Frame angle θ(t):** the rotating-frame x-axis at time `t` points along the instantaneous Sun→body[0] direction. The frame angle is `θ(t) = atan2(r_b0(t)[1], r_b0(t)[0])` — read directly from the ephemeris, **not** integrated as `∫ ω(t') dt'`. This avoids the numerical drift of integrating ω over multi-year horizons and makes the round-trip identity exact.

#### 3.1.1 Signatures

```python
def synodic_omega_dynamic(
    t_sec: float,
    bodies: tuple[str, ...],
    ephem: Ephemeris,
) -> float:
    """Instantaneous synodic angular rate (rad/s) for the body pair anchor.

    For ``len(bodies) >= 2`` returns the instantaneous rate at which the
    Sun→body[0] line is rotating in the inertial frame, computed from the
    body's heliocentric ``(r, v)``. This is the **non-uniform** synodic
    rate per spec §12(c) — the value depends on ``t_sec`` because real
    orbits are eccentric.

    For ``bodies[0]`` on a circular orbit (the M1 ``Ephemeris("circular")``
    backend) this equals ``frames.synodic_omega(bodies[0])`` to floating-
    point precision at every ``t_sec``. M6a unit-tests this degeneracy in
    ``test_dynamic_frame_reduces_to_uniform_for_circular``.

    Parameters
    ----------
    t_sec:
        Seconds since the ephemeris reference epoch.
    bodies:
        Body codes; ``bodies[0]`` defines the frame anchor. Length ≥ 2.
        Lengths > 2 use only ``bodies[0]`` for the frame angle.
    ephem:
        Heliocentric state provider. Must have ``Ephemeris("astropy")`` or
        ``Ephemeris("circular")`` behaviour.

    Returns
    -------
    Instantaneous frame rotation rate (rad/s), signed positive for prograde.
    """


def to_rotating_dynamic(
    r_inertial: Vec3,
    v_inertial: Vec3,
    t_sec: float,
    bodies: tuple[str, ...],
    ephem: Ephemeris,
) -> tuple[Vec3, Vec3]:
    """Inertial heliocentric → dynamic rotating frame anchored to Sun→body[0].

    Per spec §12(c) the rotating frame at time ``t_sec`` has its x-axis
    aligned with the instantaneous Sun→body[0] direction; the angular rate
    is the instantaneous synodic rate from
    :func:`synodic_omega_dynamic`.

    The frame angle ``θ(t) = atan2(r_b0(t)[1], r_b0(t)[0])`` is read from
    the ephemeris (not integrated), making the round-trip identity exact
    rather than just numerically close.

    The velocity transform applies the instantaneous Coriolis correction
    ``v - ω(t) × r`` before rotating; this matches the M3 uniform frame's
    structure but uses ``ω(t_sec)`` rather than a constant.

    Returns ``(r_rot, v_rot)`` in km, km/s. See
    :func:`from_rotating_dynamic` for the inverse.
    """


def from_rotating_dynamic(
    r_rot: Vec3,
    v_rot: Vec3,
    t_sec: float,
    bodies: tuple[str, ...],
    ephem: Ephemeris,
) -> tuple[Vec3, Vec3]:
    """Inverse of :func:`to_rotating_dynamic`.

    Exact inverse because ``θ(t)`` is read from the ephemeris (no numerical
    integration). Round-trip identity: ``from_rotating_dynamic(*to_rotating_dynamic(
    r, v, t, bodies, ephem), t, bodies, ephem) == (r, v)`` to ≤ 1e-10
    relative on AU-scale states. The M6a gate test
    ``test_dynamic_frame_roundtrip_identity`` asserts this.
    """
```

#### 3.1.2 Why read θ(t) instead of integrating ω(t)

Integrating `∫ ω(t') dt'` over `t_start … t` is the textbook formulation and is what introductory mechanics texts present. It is **wrong for M6a's verification purpose** because:

1. **Numerical drift over multi-lap horizons.** Multi-lap propagation runs `n_laps * synodic_period(E-M) ≈ 3 * 2.135 yr ≈ 6.4 yr ≈ 2 × 10⁸ s`. An ODE-integrated frame angle accumulates `O(dt²)` error per step; over 10⁸ s with even a 1-s noise floor per evaluation, the angle drifts a few microradians. At 1.5 AU that's tens of km — already a non-trivial fraction of the drift budget.
2. **The round-trip identity is the spec §10 binding test.** "A wrong frame silently fakes/breaks closure" — if the forward transform integrates ω from `t=0` and the inverse integrates from `t=t_sec` back to 0, **floating-point order-of-operations differences in the integrator cause the identity to fail at the 1e-6 level**, which is well above our 1e-10 binding target. Reading θ(t) directly from `atan2(r_b0_y, r_b0_x)` makes the forward and inverse transforms algebraically identical (same `atan2` evaluation, same float result, exact inverse).
3. **Ephemeris is the source of truth.** `Ephemeris("astropy").state(body, t_sec)` is the canonical position at `t_sec`; deriving θ(t) from anything else is gratuitously introducing a second source of truth that can disagree.

The cost is one extra `ephem.state(bodies[0], t_sec)` call per transform — negligible compared to the Kepler propagation inside the per-lap loop.

#### 3.1.3 Why the velocity correction still uses ω(t)

The position transform is `r_rot = R(−θ(t)) · r_inertial`. The velocity transform is `v_rot = R(−θ(t)) · (v_inertial − ω(t) × r_inertial)`. The `ω(t) × r` term is the Coriolis correction; it requires the **instantaneous** ω, not the frame angle.

`ω(t)` is read from `synodic_omega_dynamic(t, bodies, ephem)`, which computes it from `(r_b0 × v_b0) / ||r_b0||²` (or equivalently from the body's mean motion at this instant). For Earth this varies by ~3% over a year (eccentricity ≈ 0.0167); for Mars ~9% (e ≈ 0.093). The M3 uniform-frame Coriolis correction uses a constant ω; the dynamic-frame Coriolis correction uses the instantaneous value at `t_sec`. This is what makes the dynamic frame a proper time-varying rotating frame rather than just an angle-relabelling of the uniform frame.

#### 3.1.4 Why the existing M3 frame stays unchanged

Per the M5 hand-off, M5's `closure_residual` and idealized-mode optimiser consume `to_rotating` at Earth's constant mean motion. Those callers are **circular-coplanar by assumption** — the M5 contract requires `Ephemeris("circular")` and the M3 frame is exact for that backend. Switching M5 to the dynamic frame would either (a) require M5 to import an ephemeris-aware frame for circular-coplanar work where it's unnecessary, or (b) silently change `closure_residual` values for every M3/M5 test, breaking the M5 anchor test.

The dynamic frame is **only correct for ephemeris-mode** verification (M6a, M6b, M7 V2 batch validation). For circular-coplanar geometries the two frames produce the same answer to 1e-8 km (asserted in `test_dynamic_frame_reduces_to_uniform_for_circular`), but the dynamic frame costs an extra ephemeris call per sample and is overkill for M5's per-iteration optimiser objective. Keeping both is the right design.

### 3.2 Multi-lap propagation (`verify/propagate.py`)

Per spec §5 step 7 + §8 M6: "multi-lap propagation in the rotating frame (laps must overlap → periodic)." M6a propagates the cycler's full state continuously through `n_laps` laps and checks bounded drift in the dynamic rotating frame.

#### 3.2.1 `propagate_lap` — single-lap sampler

```python
def propagate_lap(
    cycler: Cycler,
    ephem: Ephemeris,
    t_start: float,
    t_end: float,
    n_samples: int,
) -> NDArray[np.float64]:
    """Propagate the cycler's spacecraft trajectory across ``[t_start, t_end]``,
    returning ``n_samples`` uniformly-spaced ``(t, x, y, z, vx, vy, vz)`` rows.

    Algorithm:
      1. For each ``Leg`` in cycler.legs whose ``[t_depart, t_arrive]`` overlaps
         ``[t_start, t_end]``, propagate from the leg's ``(r_depart, v_depart)``
         state via :func:`core.kepler.propagate` to each requested sample time
         falling within that leg.
      2. Concatenate per-leg sample arrays into one ``(n_samples, 7)`` matrix
         sorted by ``t`` ascending.

    Time grid: ``t = np.linspace(t_start, t_end, n_samples)``; each ``t``
    is assigned to the leg that contains it; the boundary ``t_arrive`` of
    leg ``i`` is assigned to leg ``i`` (not leg ``i+1``) so encounters
    appear once and only once in the sample list.

    Returned shape ``(n_samples, 7)`` is ``[t, x, y, z, vx, vy, vz]`` with
    floats in seconds, km, km/s. The function does NOT transform into any
    rotating frame — that is :func:`lap_to_lap_drift`'s job.

    Raises ``ValueError`` if ``t_start < cycler.encounters[0].t`` or
    ``t_end > cycler.encounters[-1].t + cycler.period * n_laps_max``
    where ``n_laps_max`` is set to 20 as a sanity bound — multi-lap callers
    should use :func:`multi_lap_propagation` rather than calling this with
    a multi-lap ``t_end``.
    """
```

**Design notes:**

- **No frame transform inside `propagate_lap`.** Sampling is in the inertial heliocentric frame; the drift computation happens later in `lap_to_lap_drift` and transforms each pair of samples on demand. This keeps the sampled array a clean inertial trajectory that any consumer (M7 catalogue figures, M8 viz) can use directly.
- **Universal-variable propagation per leg.** Each leg is a 2-body heliocentric arc. `core.kepler.propagate(r0, v0, dt)` natively handles forward and backward propagation; M6a uses forward only. The propagator converges in ~10 Newton iterations per call; sampling 100 points per lap costs ~1000 iterations total, well within the M6a per-test budget.
- **Encounter coincidence:** at the boundary `t == leg.t_arrive` the propagator returns the leg's `(r_arrive, v_arrive)`, which by construction equals `(encounter.r, v_planet + encounter.vinf_in)`. The sample at this `t` therefore sits exactly on the encounter; `test_propagate_lap_matches_construct_at_encounters` asserts this within the Kepler noise floor of 10 km.

#### 3.2.2 `multi_lap_propagation` — continuous N-lap propagator

```python
def multi_lap_propagation(
    cycler: Cycler,
    ephem: Ephemeris,
    n_laps: int,
    t_start: float = 0.0,
    n_samples_per_lap: int = 100,
) -> dict[str, NDArray[np.float64]]:
    """Propagate the cycler continuously through ``n_laps`` laps.

    For ``i`` in ``range(n_laps)``, the propagator advances from the cycler's
    end-of-lap state ``(r, v)`` at ``t_start + i * cycler.period`` into the
    next lap, sampling ``n_samples_per_lap`` points uniformly across each lap.

    The key word is "continuously" (spec §12 verify_long_term_stability):
    lap ``i+1`` starts from the *propagated* state at the end of lap ``i``,
    NOT from a fresh re-instantiation of the cycler. This is what makes the
    drift measurement meaningful — if the cycler is truly periodic, the
    end-of-lap-N state should sit at almost the same rotating-frame position
    as the end-of-lap-0 state.

    Returns
    -------
    dict with keys:
        ``"samples"`` — ``(n_laps * n_samples_per_lap, 7)`` matrix of
            ``[t, x, y, z, vx, vy, vz]`` rows, in time order across all laps.
        ``"lap_indices"`` — ``(n_laps + 1,)`` int array; ``lap_indices[i]``
            is the row index in ``samples`` where lap ``i`` starts.
            (``lap_indices[n_laps]`` is the total row count.)
        ``"lap_start_times"`` — ``(n_laps,)`` float array of inertial-frame
            lap-start times in seconds.
    """
```

**Design notes:**

- **Continuous, not restart.** The lap-N+1 start state is the lap-N end state. Restarting from the cycler's `encounter[0]` each lap would mask any drift — the whole point of multi-lap verification is to expose drift that compounds across laps.
- **State at lap boundaries:** the lap-N→lap-N+1 boundary state is the propagated state at `t_start + (N+1) * cycler.period`. Per the M5 hand-off, the gate cycler has `closure_residual_kms < 0.5` in the idealized frame — meaning the spacecraft state one period later is within 0.5 km/s of where the cycler "expects" it. M6a measures whether this 0.5 km/s residual stays bounded in km-of-position over `n_laps` laps in the **dynamic** frame (where the answer is no longer expected to be exactly zero).
- **Sample density:** 100 samples per lap is the M6a default. For a 2-synodic E-M cycler (period ≈ 4.27 yr) this is one sample every ~16 days, which captures the per-leg curvature without overspending on Kepler iterations. Callers can pass `n_samples_per_lap=500` for figures; the gate test uses the default 100.

#### 3.2.3 `lap_to_lap_drift` — bounded-drift measurement

```python
def lap_to_lap_drift(
    samples_lap_0: NDArray[np.float64],
    samples_lap_n: NDArray[np.float64],
    bodies: tuple[str, ...],
    ephem: Ephemeris,
) -> float:
    """Max km drift between two laps' samples, evaluated in the dynamic
    rotating frame anchored to ``bodies``.

    For each row index ``i`` in ``range(n_samples_per_lap)``, transform the
    corresponding samples in laps 0 and N into the dynamic rotating frame
    at their respective inertial times via :func:`to_rotating_dynamic`,
    and compute ``||r_0_rot - r_n_rot||``. Return the max across all ``i``.

    Per spec §12(c), the verification target is bounded drift, not zero.
    The M6a binding tolerance is ``DRIFT_TOLERANCE_KM = 50_000`` (see §4.3
    for derivation).

    Both arrays must have the same shape ``(n_samples_per_lap, 7)`` and
    be aligned positionally — i.e. row ``i`` of each represents the same
    phase within its lap.
    """
```

**Design notes:**

- **Same-phase comparison.** Row `i` of `samples_lap_0` corresponds to lap-0 phase `i/n_samples_per_lap`; row `i` of `samples_lap_n` corresponds to lap-N phase `i/n_samples_per_lap`. The dynamic-frame transform absorbs the time difference between the two samples (by definition the same phase in the rotating frame should map to the same rotating-frame position, plus or minus the eccentricity-driven breathing).
- **Max, not mean.** Spec §12 says "bounded" drift; the bound is on the max excursion, not the average. If the cycler is genuinely periodic, the max drift over a lap is the same order as the min drift over a lap; if it diverges (e.g. an unstable orbit), the max grows lap-over-lap.
- **`bodies` argument:** passed explicitly rather than read from the cycler so callers can verify, e.g., an E-M cycler in an Earth-anchored frame vs an Earth-Mars-pair-anchored frame. `verify_long_term_stability` calls this with `bodies = _resolve_frame_bodies(cycler, frame_bodies)` per §3.4.

### 3.3 `StabilityReport` dataclass + `verify_long_term_stability` — the locked spec §12 interface

```python
from dataclasses import dataclass
from typing import Final

DRIFT_TOLERANCE_KM: Final[float] = 50_000.0
"""Maximum permissible lap-to-lap drift for an E-M-class cycler over ≤ 5
laps. Derivation in plan §4.3: ~0.02° of geometric breathing per lap at
Mars's mean radius of ~1.5 AU ≈ 2.6 million km / degree → ~50,000 km /
0.02°. Tight enough to reject propagator regressions and frame-transform
errors; loose enough to absorb real eccentricity-driven breathing on
Earth-Mars trajectories. M6b may tighten this if its ephemeris-mode TCM
optimiser reports drift consistently below 10,000 km."""


@dataclass(frozen=True)
class StabilityReport:
    """Result of :func:`verify_long_term_stability`.

    All fields are immutable. The ``per_lap_dv`` and ``total_tcm_dv``
    fields are zeros in M6a; M6b populates them when
    :func:`cyclerfinder.search.optimize.optimise_cell_ephemeris` runs.

    Spec references: §8 (M6 milestone — "verified periodic over ≥3
    laps"), §12(c) (dynamic frame + tolerant verification), §14 V2
    (multi-lap periodicity validation gate).

    Attributes
    ----------
    cycler_id:
        Catalogue entry id (e.g. ``"s1l1-2syn-em-cpom"``) if the report
        was produced from a catalogue cycler, else ``None``. The
        :func:`verify_long_term_stability` caller passes this through;
        M7's batch-validate runner sets it.
    n_laps_propagated:
        Number of laps the propagation actually completed. Equals the
        ``n_laps`` argument unless an early-termination condition
        (propagator divergence, frame singularity) tripped — neither is
        expected for M6a's gate but the field exists so M7's catalogue
        can record partial results.
    max_drift_km:
        Max ``lap_to_lap_drift`` measured across all consecutive lap
        pairs (lap-0 vs lap-1, lap-1 vs lap-2, …, lap-(N-1) vs lap-N).
        The basis for ``stable``.
    max_drift_lap_index:
        The ``i`` such that ``lap_to_lap_drift(samples_lap_i,
        samples_lap_{i+1})`` equalled ``max_drift_km``. Diagnostic — lets
        a consumer point at the lap pair where drift was worst.
    per_lap_drift_km:
        Cumulative drift at each lap boundary. ``per_lap_drift_km[i]`` is
        ``lap_to_lap_drift(samples_lap_0, samples_lap_{i+1})`` — total
        drift accumulated from the start through the end of lap ``i+1``.
        ``len() == n_laps_propagated``.
    stable:
        ``max_drift_km < DRIFT_TOLERANCE_KM``. The headline boolean for
        spec §14 V2 gate; M7's signature matcher writes this directly
        into the catalogue's ``validation.gates.V2.pass`` field.
    per_lap_dv:
        **M6a: zero-tuple of length** ``n_laps_propagated``. M6b
        populates with the per-lap TCM ΔV (km/s) from the ephemeris-mode
        optimiser. Locked here so M6b doesn't reshape the dataclass.
    total_tcm_dv:
        **M6a: 0.0**. M6b populates with ``sum(per_lap_dv)``. Locked
        here so M6b doesn't reshape the dataclass.
    frame_used:
        ``"dynamic"`` (default) or ``"uniform"`` (when the caller forced
        the M3 frame, e.g. for the M3 Aldrin circular-coplanar
        regression test). Diagnostic; not part of the spec §14 V2 gate
        contract.
    """

    cycler_id: str | None
    n_laps_propagated: int
    max_drift_km: float
    max_drift_lap_index: int
    per_lap_drift_km: tuple[float, ...]
    stable: bool
    per_lap_dv: tuple[float, ...]
    total_tcm_dv: float
    frame_used: str


def verify_long_term_stability(
    cycler: Cycler,
    n_laps: int,
    ephem: Ephemeris,
    *,
    t_start: float = 0.0,
    frame_bodies: tuple[str, ...] | None = None,
    cycler_id: str | None = None,
    n_samples_per_lap: int = 100,
    use_uniform_frame: bool = False,
) -> StabilityReport:
    """Spec §12 long-term stability verifier; spec §14 V2 gate machinery.

    Pipeline:
        1. Resolve ``bodies = _resolve_frame_bodies(cycler, frame_bodies)``.
        2. Propagate continuously through ``n_laps`` via
           :func:`multi_lap_propagation`.
        3. For each consecutive lap pair ``(i, i+1)``, compute
           ``lap_to_lap_drift(samples_lap_i, samples_lap_{i+1}, bodies,
           ephem)``.
        4. Build ``per_lap_drift_km[i] = lap_to_lap_drift(samples_lap_0,
           samples_lap_{i+1}, bodies, ephem)`` for ``i ∈ range(n_laps)``.
        5. Set ``max_drift_km = max(consecutive drifts)`` and
           ``max_drift_lap_index = argmax``.
        6. ``stable = max_drift_km < DRIFT_TOLERANCE_KM``.
        7. Populate the M6b placeholders with zeros and return.

    The ``use_uniform_frame`` flag exists for the M3-circular regression
    test (``test_aldrin_cycler_periodic_over_3_laps_circular``). When
    ``True``, the drift computation uses ``to_rotating`` at Earth's
    constant mean motion; the report's ``frame_used = "uniform"``. The
    M6a binding gate runs with the default ``False`` and the dynamic
    frame.
    """
```

**Design notes:**

- **`cycler_id` is optional and caller-supplied.** M6a does not know how the cycler was derived (M5 optimiser output? catalogue entry? hand-built fixture?). The caller — M7's batch-validate runner, the catalogue ingest script — knows and passes it in.
- **`per_lap_drift_km` is cumulative-from-start, not consecutive.** This makes a divergent cycler's profile visually obvious (monotone growth) and a stable cycler's profile noise-like.
- **`max_drift_km` is the worst consecutive-pair drift**, which is the right quantity for the gate: a cycler with monotone-growing cumulative drift will have a roughly constant consecutive drift; a cycler with stochastic drift will have a max consecutive drift below `DRIFT_TOLERANCE_KM` even if the cumulative drift over 5 laps approaches the tolerance.

### 3.4 Frame-bodies resolution (`_resolve_frame_bodies`)

```python
def _resolve_frame_bodies(
    cycler: Cycler,
    frame_bodies: tuple[str, ...] | None,
) -> tuple[str, ...]:
    """Pick the dynamic-frame anchor bodies from a cycler.

    Rules:
      1. If ``frame_bodies`` is explicit, return it (caller knows).
      2. Else look at ``cycler.bodies`` (an ordered list of flyby body
         codes). For an E-M cycler this is e.g. ``["E", "M", "E"]`` or
         ``["E", "E", "M", "M"]``; the unique bodies are ``{"E", "M"}``
         and the natural frame anchor is the **pair** in the order they
         first appear: ``("E", "M")``.
      3. For longer sequences (3+ unique bodies, e.g. VEM at
         ``["E", "V", "M", "E", "M", "E"]``), pick the first two unique
         bodies in encounter order. Document this in M8 as the
         "extend-to-VEM" task.

    Returns the resolved ``bodies`` tuple; never raises (the cycler is
    guaranteed non-empty by M5 construction).
    """
```

**Why this is its own helper:** the frame anchor decision is a *policy* choice that future M7/M8 work will revisit. Centralising it in one helper makes the policy auditable and overrideable; embedding it inline in `verify_long_term_stability` would scatter the same logic across future call sites (M7 batch runner, M8 VEM verifier).

### 3.5 Imports / dependency graph after M6a

```
constants.py             (M0)
ephemeris.py             (M1 + M6 slice — astropy backend)
lambert.py               (M1)
kepler.py                (M1)
flyby.py                 (M2)
frames.py                (M3 + M6a — dynamic frame extension)
                             ↑ depends on ephemeris.py for synodic_omega_dynamic
tisserand.py             (M2)
resonance.py             (M2)
model/cycler.py          (M3)
search/construct.py      (M3)
search/sequence.py       (M4)
model/score.py           (M4)
search/optimize.py       (M5)
search/phase_match.py    (M6 slice)
verify/propagate.py      (M6a) ← frames (dynamic), kepler, ephemeris,    [new]
                                  model/cycler, constants
```

No cycles. The dynamic-frame extension to `frames.py` adds a new dependency from `frames.py` to `ephemeris.py` — previously `frames.py` only depended on `constants.py`. This is acceptable (the dynamic frame is **defined** in terms of an ephemeris) but is the one architectural change worth flagging.

### 3.6 API summary

| Symbol | Purpose | Where defined |
|---|---|---|
| `synodic_omega_dynamic(t, bodies, ephem) -> float` | Instantaneous synodic angular rate per spec §12(c) | `core/frames.py` (NEW) |
| `to_rotating_dynamic(r, v, t, bodies, ephem) -> (r_rot, v_rot)` | Inertial → dynamic rotating frame | `core/frames.py` (NEW) |
| `from_rotating_dynamic(r_rot, v_rot, t, bodies, ephem) -> (r, v)` | Inverse of `to_rotating_dynamic` | `core/frames.py` (NEW) |
| `DRIFT_TOLERANCE_KM: Final[float]` | Binding bounded-drift tolerance, 50_000 km | `verify/propagate.py` (NEW) |
| `StabilityReport` | Frozen dataclass; the locked §12 verification result shape | `verify/propagate.py` (NEW) |
| `propagate_lap(cycler, ephem, t_start, t_end, n_samples) -> NDArray` | Single-lap sampler in inertial frame | `verify/propagate.py` (NEW) |
| `multi_lap_propagation(cycler, ephem, n_laps, t_start, n_samples_per_lap) -> dict` | Continuous N-lap propagator | `verify/propagate.py` (NEW) |
| `lap_to_lap_drift(samples_0, samples_n, bodies, ephem) -> float` | Bounded-drift measurement | `verify/propagate.py` (NEW) |
| `verify_long_term_stability(cycler, n_laps, ephem, ...) -> StabilityReport` | Spec §12 entry point | `verify/propagate.py` (NEW) |
| `to_rotating`, `from_rotating`, `synodic_omega` | M3's uniform frame | `core/frames.py` (UNCHANGED) |

---

## 4. Tests + gate

Tests live under `tests/core/test_frames_dynamic.py` and `tests/verify/test_propagate.py`. Tolerances are named at the module level (`DRIFT_TOLERANCE_KM`) or hardcoded in the gate tests with cited rationale.

### 4.1 Gate tests (spec §8 + §10 binding)

| Test | Assertion | Tolerance |
|---|---|---|
| `test_2syn_em_cycler_periodic_over_3_laps` (**M6a BINDING GATE — spec §8 M6**) | Build the 2-syn E-M cycler from M5's `find_cyclers(("E","M"), k_synodic=2, vinf_cap=7.0, seed=0)[0].best_cycler` (preferred) OR from the catalogue entry `s1l1-2syn-em-cpom` (fallback fixture). `report = verify_long_term_stability(cycler, n_laps=3, ephem=Ephemeris("astropy"))`. Assert `report.stable == True`, `report.max_drift_km < 50_000`, `report.n_laps_propagated == 3`, `report.frame_used == "dynamic"`, and `report.per_lap_dv == (0.0, 0.0, 0.0)`. | `max_drift_km < 50_000`; bool/int exact |
| `test_dynamic_frame_roundtrip_identity` (**spec §10 BINDING GATE**) | For a 5-year grid of `t_sec` values (1-month spacing) and a synthetic inertial state at 1.5 AU position and 30 km/s velocity, assert `from_rotating_dynamic(*to_rotating_dynamic(r, v, t, ("E","M"), ephem), t, ("E","M"), ephem) - (r, v)` is below `1e-10 * ||(r, v)||` per-component. | rel 1e-10 |
| `test_dynamic_frame_reduces_to_uniform_for_circular` | With `ephem = Ephemeris("circular")` and `bodies = ("E", "M")`, assert `to_rotating_dynamic(r, v, t, bodies, ephem)` agrees with `to_rotating(r, v, t, synodic_omega("E"))` for the same `(r, v, t)` to ≤ 1e-8 km on position and ≤ 1e-12 km/s on velocity for a representative grid of `t_sec` and 1-AU position states. | 1e-8 km / 1e-12 km/s |

### 4.2 Multi-lap regression and stability tests

| Test | Assertion | Tolerance |
|---|---|---|
| `test_aldrin_cycler_periodic_over_3_laps_circular` | `cycler = build_aldrin_seed(Ephemeris("circular"))`; `report = verify_long_term_stability(cycler, 3, Ephemeris("circular"), use_uniform_frame=True)`. Assert `report.stable == True` and `report.max_drift_km < 1.0` (circular-coplanar is essentially zero-drift; the 1 km bound catches frame-transform regressions). | `max_drift_km < 1.0` km |
| `test_2syn_em_cycler_drifts_predictably_per_lap` | Same cycler as the gate; assert `report.per_lap_drift_km` is monotone non-decreasing or grows slowly (`per_lap_drift_km[i+1] >= per_lap_drift_km[i] * 0.9`) — i.e. drift accumulates roughly linearly with lap count, not in jumps. A jump would indicate a propagator basin-switch which is diagnostic. | Heuristic ratio bound |
| `test_5_lap_drift_stays_under_tolerance` | Same cycler, `n_laps=5`; assert `max_drift_km < DRIFT_TOLERANCE_KM` still holds. M6b will need 5-lap horizons for TCM minimisation; M6a confirms the geometry holds at the longer horizon. | `max_drift_km < 50_000` |
| `test_lap_to_lap_drift_zero_for_circular_aldrin` | `samples_0, samples_2 = multi_lap_propagation(aldrin, Ephemeris("circular"), 3)`; `drift = lap_to_lap_drift(samples_0, samples_2, ("E","M"), ephem)`; assert `drift < 100 km` (circular + Aldrin = numerical-noise floor). | < 100 km |

### 4.3 The 50,000 km drift tolerance — derivation

The tolerance is **binding** at the milestone level — not test-tunable. Derivation:

1. **Geometric scale:** at Mars's mean heliocentric radius `r ≈ 1.524 AU ≈ 2.28 × 10⁸ km`, one degree of angular position corresponds to `r * π/180 ≈ 4 × 10⁶ km`. A milliradian is `r * 0.001 ≈ 2.28 × 10⁵ km`.
2. **Real eccentricity-driven breathing:** Mars's eccentricity is `e ≈ 0.093`, Earth's is `e ≈ 0.0167`. A circular-coplanar cycler propagated on the real ephemeris will see the planets' true heliocentric positions deviate from the circular-coplanar idealisation by `~e * a` = `~14 × 10⁶ km` for Mars and `~2.5 × 10⁶ km` for Earth. The spacecraft's position drift in the rotating frame, lap-over-lap, is a **fraction** of these excursions because the cycler's flyby geometry self-corrects (a small position offset at an encounter translates into a small `vinf_in` perturbation, not a full-eccentricity drift). Empirically (from published 2-synodic E-M ephemeris-mode work, e.g. McConaghy 2006, Russell 2004 Table 4.9) the lap-to-lap drift on real ephemeris before TCMs is in the `10⁴`–`10⁵` km range.
3. **`50,000` km** is the midpoint of that empirical range. It corresponds to `50_000 / 2.28 × 10⁸ ≈ 2.2 × 10⁻⁴ rad ≈ 0.013°` of geometric breathing per lap at Mars's mean radius — small enough to be a real periodic cycler, large enough to be above the propagator's noise floor.
4. **Rejection power:** the published high-V∞ degenerate-closure basin (~11 km/s V∞) produces an open trajectory that diverges by **AU** per lap — orders of magnitude beyond 50,000 km. The tolerance therefore rejects degenerate closures cleanly while admitting genuine periodic cyclers.
5. **M6b compatibility:** when M6b's ephemeris-mode TCM optimiser runs, its target is to drive drift well below this tolerance (typically to 5,000–10,000 km after TCMs). The M6a tolerance is the **untreated** threshold — i.e. "this cycler is periodic enough that TCMs can finish the job"; the M6b tolerance is the post-TCM threshold.

The value is exposed as `DRIFT_TOLERANCE_KM: Final[float] = 50_000.0` at module scope. M6b may add a tighter `POST_TCM_DRIFT_TOLERANCE_KM: Final[float] = 10_000.0` constant for its own gate; the two coexist.

### 4.4 Unit tests on the dynamic-frame extension

| Test | Assertion |
|---|---|
| `test_synodic_omega_dynamic_matches_uniform_for_circular_earth` | For `Ephemeris("circular")` and `bodies=("E", "M")`, `synodic_omega_dynamic(t, bodies, ephem) ≈ synodic_omega("E")` to 1e-12 rad/s for `t ∈ {0, 1 yr, 2 yr, 5 yr}`. |
| `test_synodic_omega_dynamic_varies_with_eccentricity_for_astropy` | For `Ephemeris("astropy")` and `bodies=("E", "M")`, `synodic_omega_dynamic` evaluated at perihelion vs aphelion of Earth (separated by ~6 months) differs by ≥ 2% (Earth's eccentricity-driven rate variation). |
| `test_to_rotating_dynamic_anchor_x_axis` | After `to_rotating_dynamic(r_b0, v_b0, t, bodies, ephem)` of body[0]'s own state, the result's y-coordinate is below 1e-6 km — i.e. body[0] sits on the +x axis of the dynamic frame by construction. |
| `test_to_rotating_dynamic_inertial_at_t_zero_matches_uniform` | At `t = 0` on the circular backend, `to_rotating_dynamic(r, v, 0, ("E","M"), ephem) == to_rotating(r, v, 0, synodic_omega("E"))` exactly (the frames coincide at t=0). |

### 4.5 Unit tests on the propagator helpers

| Test | Assertion |
|---|---|
| `test_propagate_lap_matches_construct_at_encounters` | For the Aldrin cycler, every `Encounter.r` appears in `propagate_lap(...)`'s sampled output to ≤ 10 km — the universal-variable propagator's float-noise floor for 1-AU propagations. |
| `test_propagate_lap_n_samples_respected` | `len(propagate_lap(..., n_samples=100)) == 100`. |
| `test_propagate_lap_continuous_across_legs` | At a leg boundary, the sample before and the sample after differ by `||r||` consistent with the spacecraft's instantaneous heliocentric speed — i.e. the propagation is continuous (no flyby ΔV applied; M6a's propagator treats flybys as instantaneous ballistic rotations of V∞). |
| `test_multi_lap_propagation_lap_count` | `multi_lap_propagation(..., n_laps=3, n_samples_per_lap=100)["samples"].shape == (300, 7)`. |
| `test_multi_lap_propagation_lap_indices_monotone` | `lap_indices` is strictly increasing; `lap_indices[0] == 0`; `lap_indices[-1] == n_laps * n_samples_per_lap`. |
| `test_lap_to_lap_drift_same_lap_is_zero` | `lap_to_lap_drift(samples_0, samples_0, ("E","M"), ephem) == 0.0` exactly. |
| `test_lap_to_lap_drift_translation_in_inertial_frame_does_not_lie` | A purely translational offset of `samples_lap_n` in the **inertial** frame, when transformed to the rotating frame, still produces a non-zero drift — the dynamic frame does not silently absorb inertial-frame translations. |

### 4.6 Determinism and frozenness

| Test | Assertion |
|---|---|
| `test_stability_report_frozen_and_fields_locked` | `report.max_drift_km = …` raises `FrozenInstanceError`. `report.per_lap_dv == (0.0,) * report.n_laps_propagated`. `report.total_tcm_dv == 0.0`. `report.frame_used in ("dynamic", "uniform")`. |
| `test_verify_long_term_stability_deterministic` | Two calls with the same `cycler`, `n_laps`, `ephem`, `t_start`, `frame_bodies` produce bitwise-identical `max_drift_km` floats. The propagator is purely deterministic; no random sampling. |
| `test_verify_long_term_stability_independent_of_t_start_on_circular` | On the circular backend, `verify_long_term_stability(aldrin, 3, ephem, t_start=0)` and `verify_long_term_stability(aldrin, 3, ephem, t_start=1e6)` produce the same `max_drift_km` to numerical precision — circular ephemerides are time-translation invariant. |

### 4.7 Tolerance summary

| Layer | Quantity | Tolerance |
|---|---|---|
| Round-trip identity on dynamic frame | rel(r, v) | ≤ 1e-10 |
| Dynamic vs uniform frame on circular ephemeris | km, km/s | 1e-8 km, 1e-12 km/s |
| 2-syn E-M cycler periodic over 3 laps (M6a gate) | km drift | < 50,000 (`DRIFT_TOLERANCE_KM`) |
| Aldrin circular-coplanar drift over 3 laps | km drift | < 1.0 |
| Encounter coincidence in `propagate_lap` | km | ≤ 10 |
| Lap-0 vs lap-0 same-state drift | km | == 0.0 exact |
| Determinism across re-runs | float | bitwise |

---

## 5. Risks

| # | Risk | Likelihood | Impact | Mitigation in M6a |
|---|---|---|---|---|
| 1 | **Closure-frame correctness (spec §10 headline risk).** A wrong dynamic-frame transform silently fakes/breaks closure — every drift measurement downstream would be wrong, and the M6a gate would falsely pass or fail. | medium | **high — milestone-blocking** | The §3.1.2 design — reading θ(t) directly from `atan2(r_b0_y, r_b0_x)` rather than integrating ω(t) — makes the round-trip identity algebraically exact. The §4.1 gate `test_dynamic_frame_roundtrip_identity` is the binding test (≤ 1e-10 rel). The §4.1 `test_dynamic_frame_reduces_to_uniform_for_circular` test is the cross-validation: at the limit where the M3 uniform frame is provably correct, the dynamic frame must agree. Both tests must pass before any cycler-level drift test runs. |
| 2 | **2-syn E-M cycler drift exceeds 50,000 km (gate fails).** Possible causes: (a) the M5-produced cycler has higher closure residual than the 0.5 km/s budget suggests; (b) the catalogue's `s1l1-2syn-em-cpom` parameters are too coarse for ephemeris-mode propagation; (c) the dynamic-frame anchor `("E", "M")` is wrong (should be Sun-Earth only, with Mars as a tracked encounter body). | medium | high — gate-blocking | First, the test reports `max_drift_km` and `max_drift_lap_index` so the failure mode is visible. Mitigations in priority order: (a) verify the M5 cycler's `closure_residual_kms` from M5's hand-off; if > 0.5, escalate to M5 — drift cannot be smaller than the velocity-residual budget integrated over a lap. (b) Try both Earth-only frame `("E",)` and pair frame `("E", "M")` — the spec §12(c) wording is ambiguous between the two, and `_resolve_frame_bodies` makes the policy auditable. (c) Try `n_laps=3` vs `n_laps=2` to see whether drift is monotone or jumps; a jump means a propagator basin-switch (rare on Kepler). **Do NOT widen `DRIFT_TOLERANCE_KM`** — that masks the underlying issue. |
| 3 | **Universal-variable propagator accumulates float noise over multi-year propagation.** A single Kepler call at 6.4-yr horizon has known noise ~10 km; chained across legs the noise can grow. | low | medium | The 10-km coincidence-at-encounters test (`test_propagate_lap_matches_construct_at_encounters`) is the floor measurement. The 50,000-km drift tolerance is 5,000x this floor, so propagator noise is dominated by real eccentricity-driven breathing. If `test_lap_to_lap_drift_zero_for_circular_aldrin` exceeds 100 km, the propagator is degrading — escalate by tightening the Kepler Newton tolerance or switching to a higher-fidelity integrator (scipy.integrate.solve_ivp with DOP853) for the long-horizon case. Documented as a M6a-future-work follow-up if it ever fires. |
| 4 | **Astropy ephemeris call cost dominates the verify_long_term_stability runtime.** Each dynamic-frame transform requires one `ephem.state(bodies[0], t)` call; multi-lap propagation across 300 samples requires 300+ calls. astropy's DE440 lookup is ~10 ms per call → ~3 s per gate test. | medium | low | Per-gate test runtime ≤ 30 s is acceptable for M6a CI. If runtime becomes a problem (M7 batch validation across hundreds of catalogue entries), the astropy backend can be cached at the `Ephemeris` layer (memoise `state(body, t)` to ~1 µs per call). Documented as M7 follow-up. M6a does not optimise this; the gate test budget is met as-is. |
| 5 | **The dynamic frame extension to `core/frames.py` adds an `ephemeris.py` import dependency.** Previously `frames.py` depended only on `constants.py`; now it transitively depends on astropy via `ephemeris.py`. A user importing `frames.to_rotating` (the M3 path) would pay the astropy import cost they didn't ask for. | low | low | The `Ephemeris` parameter is **passed in** to `synodic_omega_dynamic` / `to_rotating_dynamic`, not constructed inside. M3's `to_rotating` / `synodic_omega` continue to work without ever touching the astropy backend — they don't import `Ephemeris`. The new functions accept an `Ephemeris` instance as a parameter, so importing `frames` itself does not eagerly import astropy. Verified by `test_frames_no_astropy_import_on_load` (a small static test). |
| 6 | **Dataclass shape lock-in conflicts with M6b's TCM-budget needs.** If M6b needs additional fields on `StabilityReport` (e.g. per-encounter TCM breakdown, per-leg ΔV), the M6a-locked shape forces M6b to either reshape or add a parallel dataclass. | low | low | The frozen dataclass is intentionally minimal. M6b can extend by either (a) adding a new field with a default value (`per_encounter_dv: tuple[float, ...] = ()`) — backward-compatible, or (b) adding a new dataclass `TCMReport` and returning it alongside `StabilityReport`. The spec §12 `verify_long_term_stability` signature returns a dict in the original wording; we keep it returning `StabilityReport` and document that consumers should treat additional fields as optional. M6b plan revisits at its writing time. |
| 7 | **M5's `find_cyclers` may not produce a cycler suitable for M6a verification at test time.** The M6a gate test depends on M5 output (`find_cyclers(("E","M"), k_synodic=2, vinf_cap=7.0, seed=0)[0]`). If M5's gate test ran with different seed/tolerance defaults, the M6a gate gets a different cycler. | medium | medium | The gate test has a **catalogue fallback**: if M5's `find_cyclers` either fails to produce a result or produces one outside the 5.65/3.05 V∞ tolerance, the test falls back to constructing a cycler from `s1l1-2syn-em-cpom`. This is documented in `test_2syn_em_cycler_periodic_over_3_laps`'s docstring. The catalogue entry is the canonical reference; M5's reproduction is the convenience case. |
| 8 | **`_resolve_frame_bodies` policy is wrong for multi-body cyclers.** For a 3-encounter sequence like `["E", "M", "E"]`, the resolver returns `("E", "M")`; for a 4-encounter `["E", "E", "M", "M"]` it returns `("E", "M")` too. But the spec §12(c) wording is "Sun-Earth (or Earth-Mars)" — the Earth-only frame may be more correct for some cyclers. | low (in M6a) | medium (in M8 VEM) | The helper is small and overrideable via the `frame_bodies` argument. M6a's gate uses `frame_bodies=None` (default) and exercises the policy; if the gate's max_drift is borderline, the gate's diagnostic output (per the M5 hand-off pattern) lets the operator try `frame_bodies=("E",)` or `frame_bodies=("E", "M")` explicitly. The M8 VEM extension is a separate plan-doc topic. |
| 9 | **astropy's DE440 ephemeris may not be available at CI time.** The astropy backend lazy-loads DE440; if astropy's bundle is corrupted or the DE440 download fails on a CI runner, every M6a test that uses `Ephemeris("astropy")` fails. | low | low | astropy 7.0+ bundles DE440 (per `pyproject.toml` dep `astropy>=7.0`). The M6 slice's `_AstropyBackend.__init__` calls `solar_system_ephemeris.set("de440")` at construction; this is the failure point. CI uses uv.lock so astropy version is pinned; a regression of DE440 availability is unlikely. M6a's circular-backend tests (`test_aldrin_cycler_periodic_over_3_laps_circular`, `test_dynamic_frame_reduces_to_uniform_for_circular`) pass without astropy data, so a partial test pass is still informative. |
| 10 | **Catalogue Schema-v2 backfill agent is modifying `data/seed_cyclers.yaml` concurrently.** If M6a's fallback fixture reads the catalogue, it may read mid-write or read a schema that's still settling. | low | low | M6a does NOT touch `data/` or `docs/spec.md` per the parallel-agent constraint. The fallback fixture loads `s1l1-2syn-em-cpom` via the existing `tests/_catalogue_loader.py` (which the parallel agent leaves alone — it's test infrastructure, not data). The loader is tolerant of schema additions (it reads `bodies`, `orbit_elements`, `legs`, `vinf_kms_at_encounters` — fields the parallel agent does not remove). M6a's tests should run after the parallel agent's branch lands; sequential coordination is the parent agent's job. |
| 11 | **`frame_used` field clutters StabilityReport for the common case.** Every caller will see `frame_used="dynamic"`; the field is only diagnostic when someone deliberately sets `use_uniform_frame=True`. | low | low | Accepted. The field is cheap; its absence would mean the M3 circular-coplanar regression test couldn't be expressed in the same dataclass shape. Keeping it makes the dataclass self-describing for downstream consumers (M7 catalogue records). |
| 12 | **`per_lap_dv` and `total_tcm_dv` zeros could be confused with "no TCMs needed."** A downstream consumer reading `total_tcm_dv == 0.0` may mistakenly conclude the cycler needs no maintenance. | low | medium | Document on the dataclass: "M6a: zero by construction; M6b populates." The M7 catalogue ingestion script must consume only `StabilityReport` instances where the producer guaranteed the M6b fields are populated — i.e. it must check the cycler ran through `optimise_cell_ephemeris`, not just `verify_long_term_stability`. M7's plan will encode this rule. |

---

## 6. Dependency additions

**None.** M6a uses only:

- `numpy` — already in deps.
- `scipy` — already in deps (M5 dep; M6a doesn't add a new scipy module beyond what M5 used).
- `astropy` — already in deps (M6 slice introduced).
- in-house M0–M6 modules.

`scipy.integrate.solve_ivp` is **considered and rejected** for M6a: the universal-variable Kepler propagator in `core/kepler.py` is already validated and faster (`solve_ivp` with `DOP853` would add ODE-integration overhead for a problem that has an analytic solution). Documented as a M6-future-work option only if `test_lap_to_lap_drift_zero_for_circular_aldrin` ever fires above its 100 km bound.

No edits to `pyproject.toml`. No `uv.lock` regeneration.

---

## 7. Order of work

The `todo.md` mirrors this with checkboxes.

1. **Re-read predecessor docs.** Re-read spec §12(c) (the binding dynamic-frame decision) and §10 (the binding round-trip identity test); re-read M3's `frames.py` plan §3.2 (the uniform-frame structure that M6a extends); re-read M5's `## Hand-off to M6a` section (TBD values to consume); re-read the M6 slice's `phase_match.py` to confirm it does not depend on M6a (it produces inputs for M6a, not the reverse); confirm `core/ephemeris.py` `Ephemeris("astropy")` backend works (the M6 slice's smoke test).
2. **Extend `core/frames.py` with the dynamic frame.** Strictly additive: add `synodic_omega_dynamic`, `to_rotating_dynamic`, `from_rotating_dynamic`. Do NOT touch the M3 functions. Add module-docstring note that the M3 uniform frame and the new dynamic frame coexist. Write paired tests in `tests/core/test_frames_dynamic.py`:
   - `test_dynamic_frame_roundtrip_identity` (**spec §10 BINDING GATE**) first — this is the headline correctness test.
   - `test_dynamic_frame_reduces_to_uniform_for_circular` second — the cross-validation.
   - `test_synodic_omega_dynamic_matches_uniform_for_circular_earth` and `test_synodic_omega_dynamic_varies_with_eccentricity_for_astropy` third — confirms the rate computation is right.
   - `test_to_rotating_dynamic_anchor_x_axis` and `test_to_rotating_dynamic_inertial_at_t_zero_matches_uniform` last.
   - Confirm `uv run pytest tests/core/test_frames_dynamic.py` green and `uv run mypy src/cyclerfinder/core/frames.py` clean.
3. **Create the `verify/` subpackage.** Empty `src/cyclerfinder/verify/__init__.py`; empty `tests/verify/__init__.py`. Confirm `mypy` still clean.
4. **Implement `verify/propagate.py` skeleton.** Module docstring with spec §12(c), §14 V2 references. `DRIFT_TOLERANCE_KM` constant. `StabilityReport` frozen dataclass per §3.3 with all eight fields. All five public functions stubbed with full docstrings and `NotImplementedError` bodies. Run `mypy --strict` — confirm type signatures clean before any logic lands.
5. **Implement helpers in dependency order.** Each helper gets a paired unit test before the next:
   - `_resolve_frame_bodies(cycler, frame_bodies)` — small policy helper; test: returns `("E", "M")` for Aldrin's `["E","M","E"]` cycler, returns explicit override when given.
   - `propagate_lap(cycler, ephem, t_start, t_end, n_samples)` — single-lap sampler; tests: `test_propagate_lap_matches_construct_at_encounters`, `test_propagate_lap_n_samples_respected`, `test_propagate_lap_continuous_across_legs`.
   - `multi_lap_propagation(cycler, ephem, n_laps, t_start, n_samples_per_lap)` — continuous N-lap propagator; tests: `test_multi_lap_propagation_lap_count`, `test_multi_lap_propagation_lap_indices_monotone`, `test_lap_to_lap_drift_zero_for_circular_aldrin` (the integration test).
   - `lap_to_lap_drift(samples_0, samples_n, bodies, ephem)` — bounded-drift measurement; tests: `test_lap_to_lap_drift_same_lap_is_zero`, `test_lap_to_lap_drift_translation_in_inertial_frame_does_not_lie`.
6. **Implement `verify_long_term_stability`.** Compose the helpers per §3.3 pseudocode. Land:
   - `test_aldrin_cycler_periodic_over_3_laps_circular` first (circular, known-good answer ≈ 0 drift).
   - `test_2syn_em_cycler_periodic_over_3_laps` (**M6a BINDING GATE**) second — the spec §8 M6 anchor. If it fails, escalate per risk #2 (DO NOT widen `DRIFT_TOLERANCE_KM`).
   - `test_2syn_em_cycler_drifts_predictably_per_lap`, `test_5_lap_drift_stays_under_tolerance`, `test_stability_report_frozen_and_fields_locked`, `test_verify_long_term_stability_deterministic`, `test_verify_long_term_stability_independent_of_t_start_on_circular`.
7. **Run the full local quality gate:** `uv run pytest`, `uv run ruff check .`, `uv run ruff format --check .`, `uv run mypy src tests`.
8. **Commit** as `m6a: dynamic rotating frame + verify/propagate.py — 2-syn E-M periodic over 3 laps`. Push; confirm CI green.
9. **Update `docs/overview.md`.** §4 milestone table: M6a → completed; M6b → planned.
10. **Hand-off note** appended to `todo.md` under `## Hand-off to M6b` covering: actual `max_drift_km` reproduced for the gate, frame-bodies decision (Earth-only vs Earth-Mars-pair), per-test runtime to inform M6b's optimiser-inner-loop budget, any escalation actions taken from risks #2/#3, and the locked `StabilityReport` shape that M6b inherits.

The order is "dynamic frame → verify/ subpackage → propagator helpers → entry-point function" deliberately: each step depends only on what came before; the spec §10 frame-correctness gate lands first; the M6a binding multi-lap gate (spec §8 M6) lands last.

---

## 8. Exit checklist (the gate, restated)

Before declaring M6a done:

- [ ] `uv run pytest tests/core/test_frames_dynamic.py` green; the spec §10 round-trip identity passes at 1e-10 rel.
- [ ] `uv run pytest tests/verify/test_propagate.py` green; the M6a binding gate `test_2syn_em_cycler_periodic_over_3_laps` passes with `max_drift_km < 50_000` km, `stable=True`, `n_laps_propagated=3`.
- [ ] `uv run pytest` green overall (no regression of M0–M5 tests; M3 `test_frames.py` still passes — i.e. the M3 uniform frame is unchanged).
- [ ] `uv run ruff check .` clean.
- [ ] `uv run ruff format --check .` clean.
- [ ] `uv run mypy src tests` clean under `strict = true` — including `StabilityReport`, `Final[float]` on `DRIFT_TOLERANCE_KM`, NDArray shape annotations on the propagator outputs.
- [ ] CI green on the M6a commit.
- [ ] `docs/overview.md` updated: M6a status = `completed`; M6b row marked `planned`; the dynamic-frame decision noted in §2 (M3's uniform frame coexists with M6a's dynamic frame; circular-coplanar work uses uniform, real-ephemeris verification uses dynamic).
- [ ] `## Hand-off to M6b` section appended to `phases/m6a-idealized-closure-verification/todo.md` covering:
  - The exact `max_drift_km` reproduced for the 2-syn E-M gate (vs the 50,000 km bound), and which lap pair was worst (`max_drift_lap_index`).
  - The frame-bodies decision: was the gate's `_resolve_frame_bodies` policy (`("E", "M")`) right, or did the gate require an override?
  - Per-test wall-clock runtime (informs M6b's optimiser-inner-loop budget — the ephemeris-mode optimiser calls `verify_long_term_stability` inside its objective).
  - Whether `test_5_lap_drift_stays_under_tolerance` passed — informs M6b's horizon length (3 laps is the minimum gate; 5 laps is the spec §12(a) horizon).
  - Whether any of risks #2 (drift > 50k), #3 (propagator noise), #5 (frames import of astropy) fired during implementation, and how they were resolved.
  - Whether M5's `find_cyclers(("E","M"), k_synodic=2, vinf_cap=7.0)[0]` was actually used as the gate fixture, or the catalogue fallback was needed.
  - The locked `StabilityReport` shape (all eight fields) so M6b knows it does NOT need to reshape the dataclass — only populate `per_lap_dv` and `total_tcm_dv`.

(Writing the M6b plan doc is the first task of M6b, not an M6a exit criterion.)
