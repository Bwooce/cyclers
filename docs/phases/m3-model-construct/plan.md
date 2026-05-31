# M3 — Model + construction

**Spec reference:** spec.md §3 (closure definition), §4 (architecture), §5 step 3 (seed construction), §6 (interfaces — `Cycler`/`Leg`/`Encounter`), §8 (M3 milestone definition), §9 (Aldrin validation anchors), §10 (closure-frame correctness risk), §12(c) (dynamic ephemeris frame — explicitly deferred to M6).

**Purpose:** stand up the dataclass model for a cycler, the synodic rotating-frame transforms used by closure, and a patched-conic leg constructor — then prove the foundation by **reproducing the Aldrin cycler** (a ≈ 1.659 AU, e ≈ 0.41, E→M ≈ 146 d) end-to-end from the M1 ephemeris + Lambert primitives. M3 is the first milestone that produces a *named object the rest of the project will recognise* — `Cycler` — and the first that gates on a published number rather than an internal solver check.

**Gate (definition of done):** `tests/model/test_aldrin.py` reproduces the Aldrin orbit's `(a, e, perihelion, aphelion)` within the tolerances of §4 below, `tests/core/test_frames.py` passes a sub-1e-10 km rotating-frame round-trip, and `tests/model/test_cycler.py::test_closure_residual_aldrin` shows a small (numerically-limited) closure residual for the constructed Aldrin cycler.

---

## 1. What this milestone delivers

A package layer that, given the M0–M2 primitives, can:

1. **Represent** a cycler trajectory as a frozen dataclass tree (`Cycler` → `[Encounter]` + `[Leg]`) carrying enough state for every downstream metric (maintenance ΔV, closure residual, radial span, max V∞).
2. **Transform** heliocentric inertial states into a synodic rotating frame at the body pair's mean motion, and back. The frame is the **uniform** rotating frame anchored to Earth's mean motion — exact for the M1 circular-coplanar ephemeris.
3. **Construct** a patched-conic cycler from a flyby sequence and a list of encounter times by Lambert-solving each leg against the M1 `Ephemeris`. The constructor consumes encounter times as inputs; it does not search over them (that is M5).
4. **Reproduce the Aldrin cycler** from `build_aldrin_seed(ephem)` and pass spec §9's numeric checks.
5. **Reproduce a 2-synodic E–M cycler** (period ≈ 4.27 yr) as an `E-M-E` sequence and demonstrate that `closure_residual()` is small and computable, not as an optimisation target but as a sanity check that the geometry, frames, and residual math line up.

### Out of scope (M3)

- **No optimiser.** `construct_cycler` consumes encounter times as inputs. Searching encounter times to minimise the residual is M5; searching sequences is M4.
- **No real ephemeris.** Aldrin's elements come from Earth–Mars *circular* semi-major axes (1.0 AU and 1.524 AU); with the M6 astropy backend the numbers shift slightly (reported in M6b).
- **No dynamic (ephemeris-anchored) rotating frame.** Spec §12(c)'s non-uniform frame and the "geometric breathing" tolerance live in M6/M6b. M3's frame is exact for the M1 ephemeris; calling `closure_residual()` against a real-ephemeris cycler would silently mis-measure, so this scope boundary must be enforced at the call site (M6 plan introduces the dynamic frame and rewires the residual).
- **No catalogue lookup.** Aldrin's name appears only as a fixture/seed builder; signature-matching is M7.
- **No multi-lap propagation.** The closure test rotates the *endpoint state* by ω·T and compares to start; it does not propagate the trajectory through multiple laps (that is M6 V2).
- **No bend feasibility, no powered-flyby ΔV.** `maintenance_dv()` sums `flyby_dv` at each encounter, but for the idealized Aldrin cycler this is 0 by construction; the sum is plumbing for M5+.

---

## 2. File tree after M3

```
cyclers/
├── src/cyclerfinder/
│   ├── core/
│   │   ├── constants.py            # M0
│   │   ├── ephemeris.py            # M1
│   │   ├── lambert.py              # M1
│   │   ├── kepler.py               # M1
│   │   ├── flyby.py                # M2
│   │   └── frames.py               # NEW (M3)
│   ├── search/
│   │   ├── tisserand.py            # M2
│   │   ├── resonance.py            # M2
│   │   └── construct.py            # NEW (M3)
│   └── model/
│       ├── __init__.py             # NEW (M3) — re-export Cycler, Leg, Encounter
│       └── cycler.py               # NEW (M3)
└── tests/
    ├── core/
    │   └── test_frames.py          # NEW (M3)
    └── model/
        ├── __init__.py             # NEW (M3)
        ├── test_aldrin.py          # NEW (M3) — gate test
        ├── test_cycler.py          # NEW (M3) — dataclass + closure
        └── test_construct.py       # NEW (M3) — 2-synodic E-M-E
```

Directories `verify/`, `data/`, `viz/` remain uncreated — they belong to M6/M7/M8. `search/sequence.py`, `search/optimize.py`, `model/score.py` likewise wait for M4/M5.

`model/__init__.py` re-exports the three dataclasses so callers write `from cyclerfinder.model import Cycler, Leg, Encounter` instead of reaching into `cycler.py` directly. This is the only convenience aggregator in M3.

---

## 3. Module designs

### 3.1 `model/cycler.py`

Three `frozen=True` dataclasses, plus methods on `Cycler`. Frozen because a cycler is a *result object* — modifying its encounters in place would corrupt every metric already derived from it. Frozen also makes `Cycler` hashable for the M7 catalogue, but only structurally; we will not implement `__hash__` here (numpy arrays are not hashable and the signature-based hash from spec §16.2 is the right identity, not Python object identity).

#### Type aliases (top of module)

```python
from numpy.typing import NDArray
import numpy as np

Vec3 = NDArray[np.float64]   # shape (3,), heliocentric km or km/s
```

`Vec3` is module-private (single underscore would be conventional but mypy treats `NDArray` aliases as fine without). Every array field below carries this annotation; mypy strict will enforce it. Where a method needs a 3×3 rotation matrix we use `NDArray[np.float64]` directly with a shape comment — there is no clean shape generic in numpy 2.x that survives mypy strict without `typing.cast`.

#### `Encounter`

```python
@dataclass(frozen=True)
class Encounter:
    body: str               # one-letter code matching constants.PLANETS keys ("V", "E", "M")
    t: float                # seconds from cycler epoch (t=0 at first encounter)
    r: Vec3                 # heliocentric position of the planet (km)
    v_planet: Vec3          # heliocentric velocity of the planet (km/s)
    vinf_in: Vec3           # spacecraft V∞ approach vector  = v_sc_arrive - v_planet
    vinf_out: Vec3          # spacecraft V∞ departure vector = v_sc_depart - v_planet
```

Notes:

- `vinf_in` at the **first** encounter and `vinf_out` at the **last** encounter of a non-closed sequence are by convention the same as their counterpart (the cycler closes, so the "first arrival" V∞ equals the "last departure" V∞ rotated by the synodic frame's advance). The constructor populates both for the closed (cyclic) case; tests assert `|vinf_in| ≈ |vinf_out|` per the ballistic-flyby invariant.
- `body` is a one-letter code, not the full name, for symmetry with `PLANETS["E"]` in `constants.py`. The dataclass does no validation of the code; that responsibility sits in the constructor.
- All arrays are expected to be C-contiguous, 1-D, length 3, float64. We do not enforce this in `__post_init__` (frozen + mypy strict gives us most of the type guarantee; a runtime shape check is one M3 nice-to-have if a test ever fires on bad input).

#### `Leg`

```python
@dataclass(frozen=True)
class Leg:
    from_body: str
    to_body: str
    t_depart: float          # seconds from cycler epoch
    t_arrive: float
    v_depart: Vec3           # spacecraft heliocentric velocity at departure (km/s)
    v_arrive: Vec3           # spacecraft heliocentric velocity at arrival
    n_revs: int = 0          # heliocentric revolutions on this leg (0 = direct)
    branch: str = "single"   # "single" for n_revs=0; "low" | "high" for n_revs ≥ 1
```

Notes:

- `v_depart` and `v_arrive` are *spacecraft* velocities, not V∞. The encounter's `vinf_out = leg.v_depart - encounter.v_planet`; the constructor enforces this.
- `branch` follows spec §12(b)'s `"low" | "high"` vocabulary for multi-rev Lambert. Default `"single"` for the 0-rev case keeps the default sensible and matches the M1 Lambert solver's `branch="low"` convention for `revs=0` (a 0-rev Lambert has a single branch; "single" reads better than the misleading "low").
- `t_arrive > t_depart` is asserted in the constructor, not in `__post_init__` (consistency with `Encounter`).

#### `Cycler`

```python
@dataclass(frozen=True)
class Cycler:
    bodies: list[str]                     # ordered flyby sequence, e.g. ["E", "M", "E"]
    period: float                         # total period (s) — typically k * synodic
    encounters: list[Encounter]           # length == len(bodies); ordered by time
    legs: list[Leg]                       # length == len(bodies) - 1 for open; len(bodies) for closed-loop reps

    def maintenance_dv(self) -> float: ...
    def closure_residual(self, omega_rad_per_s: float | None = None) -> float: ...
    def radial_span(self) -> tuple[float, float]: ...
    def max_vinf(self) -> float: ...
```

Method semantics:

| Method | Returns | Definition |
|---|---|---|
| `maintenance_dv()` | km/s | `sum(np.linalg.norm(enc.vinf_out - enc.vinf_in) for enc in encounters)` if a flyby cannot achieve the bend; the M3 implementation returns the **velocity discontinuity sum** without checking bend feasibility (no `flyby_dv` call in M3). M4+ replaces this with `core.flyby.flyby_dv(...)`. For an idealized cycler this is 0 by construction; for the test we assert it is < 1e-3 km/s. |
| `closure_residual(omega)` | km/s | After period T, the spacecraft state `(r, v)` at the *last* encounter (rotated by the synodic frame's advance over T) must equal the state at the *first* encounter. Residual = `‖v_first - rot(v_last, ω·T)‖` in km/s. If `omega` is `None`, the method uses `frames.synodic_omega("E")` (Earth's mean motion is the default frame for E–M cyclers per spec §3). For VEM cyclers the caller passes a different ω. |
| `radial_span()` | `(min_au, max_au)` | The per-leg analytic `(perihelion, aphelion)` computed from each leg's departure state `(r, v_depart)` via the standard `a = 1/(2/r - v²/μ)`, `e = ‖e_vec‖`. M3 returns the min perihelion / max aphelion across legs, in **AU** (the unit the spec §9 anchors use). Leg-level `(a, e)` is recomputed each call — no caching on the frozen dataclass. |
| `max_vinf()` | km/s | `max(np.linalg.norm(enc.vinf_in) for enc in encounters)`. Uses `vinf_in` because in steady cycling `‖vinf_in‖ == ‖vinf_out‖`; taking only one avoids double-counting. |

Implementation notes:

- All methods are pure reads of state — no mutation, no caching attribute (frozen dataclass cannot hold one cleanly). If profiling later shows them hot, the right fix is `functools.cached_property` on a non-frozen subclass or a separate `score.py` helper that memoizes externally; M3 does neither.
- `closure_residual` is the **velocity** residual only, not full state. Position closure is automatic for the geometric construction (the last encounter happens at body B at time T; by construction the next encounter would be at body A at time T, which is exactly where body A is in the rotating frame). What can fail is the velocity match — that is the maintenance ΔV.
- `period` is set by the constructor (typically `k * synodic_period(pair)` from M2). `Cycler` itself does nothing with it beyond storing it for closure and metrics.

---

### 3.2 `core/frames.py`

Three functions; no classes. The frame is fully specified by `omega_rad_per_s` — there is no `Frame` object because in M3 there is exactly one frame (the uniform synodic frame at one ω). The dynamic frame from spec §12(c) will introduce a `Frame` abstraction in M6.

#### Frame definition (state this in the module docstring)

> **Uniform synodic rotating frame.** Rotation about the heliocentric z-axis at constant angular rate ω, with the frame's x-axis aligned with the inertial x-axis at t=0. A state `(r_i, v_i)` in the inertial heliocentric frame maps to the rotating frame by `r_r = R(−ωt)·r_i`, `v_r = R(−ωt)·(v_i − ω×r_i)`, where `R(θ)` is the right-handed rotation matrix about z. This frame is exact for the M1 circular-coplanar ephemeris. **It is NOT correct for the real (eccentric, inclined) ephemeris** — see spec §12(c) and M6b for the dynamic frame.

#### Signatures

```python
def to_rotating(
    r_inertial: NDArray[np.float64],
    v_inertial: NDArray[np.float64],
    t: float,
    omega_rad_per_s: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Inertial heliocentric -> uniform rotating frame at angular rate omega about z.
    Returns (r_rot, v_rot) in km, km/s."""

def from_rotating(
    r_rot: NDArray[np.float64],
    v_rot: NDArray[np.float64],
    t: float,
    omega_rad_per_s: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Inverse of to_rotating. (from_rotating ∘ to_rotating) is identity to ~1e-13 rel."""

def synodic_omega(body: str) -> float:
    """Angular rate (rad/s) of the synodic rotating frame for an Earth-<body> pair.

    For Earth-Mars cyclers this is Earth's mean motion (a heliocentric circular orbit
    at 1 AU); the spacecraft's repeating geometry is described in the frame that follows
    Earth's longitude. For Earth-Venus or VEM cyclers the caller chooses which body
    anchors the frame (typically the slowest member of the pair set).

    M3 returns Earth's mean motion for body == "E" or body == "M" (the E-M case); raises
    NotImplementedError for "V" (Venus-anchored frames arrive with the VEM campaign in M8).
    """
```

The `synodic_omega("E")` value is derived from `constants.PLANETS["E"].mean_motion_deg_day` converted to rad/s. This makes the constants module the single source of truth, and avoids hard-coding a second copy of Earth's mean motion.

#### Implementation notes

- The rotation matrix is built inline as a 3×3 `np.array` from `cos(ωt)`, `sin(ωt)`. No `scipy.spatial.transform.Rotation` dependency — overkill for a single z-axis rotation and would force an additional scipy import for M3 (M1 already pulled scipy in for Lambert root-finding, so this is just a style preference; we keep the M3 footprint minimal).
- The frame's z-axis is the ecliptic normal (the M1 ephemeris places all planets in the z=0 plane). For real ephemeris in M6 this assumption breaks (Mars has 1.85° inclination); the dynamic frame is then anchored to instantaneous Sun-Earth geometry per spec §12(c).
- `ω×r` for ω along z is `(-ω·r_y, ω·r_x, 0)` — written out by hand to avoid `np.cross` overhead in the round-trip test.
- Both functions accept 1-D float64 arrays of shape (3,). They return new arrays (no in-place mutation).

---

### 3.3 `search/construct.py`

Patched-conic cycler construction. Given a sequence and encounter times, produces a `Cycler`. No search, no optimisation.

#### Signatures

```python
def construct_cycler(
    sequence: list[str],
    encounter_times_sec: list[float],
    ephem: Ephemeris,
    mu_sun: float = MU_SUN_KM3_S2,
    max_revs_per_leg: list[int] | None = None,
    branch_per_leg: list[str] | None = None,
) -> Cycler:
    """Patched-conic construction of a cycler from a fixed schedule.

    Each leg is solved by Lambert between consecutive encounter positions. V∞ vectors
    at each encounter are derived from the spacecraft heliocentric velocity minus the
    planet's heliocentric velocity at that instant.

    Inputs:
        sequence: ordered body codes, e.g. ["E", "M", "E"] for a 2-synodic E-M-E cycler.
                  len(sequence) == n_encounters; legs[i] connects sequence[i] -> sequence[i+1].
        encounter_times_sec: encounter epochs in seconds from t=0. Strictly increasing.
                             len == len(sequence).
        ephem: M1 Ephemeris (circular-coplanar in M3).
        max_revs_per_leg: optional per-leg multi-rev hint. Default = [0] * (n_legs).
        branch_per_leg: optional per-leg branch selection in {"single","low","high"}.
                        Default = ["single"] * n_legs (0-rev solutions).

    Returns:
        Cycler with encounters and legs fully populated. period == t[-1] - t[0]
        when the sequence is open; for a closed-loop representation (sequence[-1]
        equals sequence[0] same body, by convention), period is still t[-1] - t[0].
    """

def build_aldrin_seed(
    ephem: Ephemeris,
    t_start_sec: float = 0.0,
) -> Cycler:
    """Build the canonical Aldrin cycler as a 2-encounter E -> M slice
    (one synodic period; the "outbound escalator" leg).

    Places Earth at t_start using the supplied ephemeris, places Mars at
    t_start + 146 days, and Lambert-solves the connecting leg with n_revs=0.

    Returns a Cycler with bodies=["E","M"], 2 encounters, 1 leg. This is
    enough to recover Aldrin's (a, e, perihelion, aphelion) — the test
    extracts them from the leg's departure (r, v) and asserts spec §9 values.

    The full closed Aldrin cycle would be E -> M -> E over one synodic period;
    the M3 gate test verifies only the E -> M slice (matching the published
    elements). The closed-loop case is exercised by the 2-synodic E-M-E test
    in test_construct.py, which uses construct_cycler directly.
    """
```

#### Algorithm (construct_cycler)

1. Validate inputs: `len(sequence) == len(encounter_times_sec) >= 2`; times strictly increasing; every body code in `constants.PLANETS`.
2. Default `max_revs_per_leg = [0] * (n - 1)`, `branch_per_leg = ["single"] * (n - 1)`.
3. For each encounter `i`, call `ephem.state(sequence[i], encounter_times_sec[i])` → `(r_i, v_planet_i)`.
4. For each leg `j` (from encounter `j` to encounter `j+1`):
   - Compute `tof = encounter_times_sec[j+1] - encounter_times_sec[j]`.
   - Call `lambert(r_j, r_{j+1}, tof, mu=mu_sun, prograde=True, max_revs=max_revs_per_leg[j])`.
   - Select the solution matching `branch_per_leg[j]` (for `n_revs=0`, only `"single"` exists; for `n_revs ≥ 1`, pick the matching branch). Raise `ValueError` if no matching branch.
   - Store `(v_depart, v_arrive) = (sol.v1, sol.v2)`.
5. Build `Encounter` objects:
   - First encounter: `vinf_in = v_depart_0 - v_planet_0` (placeholder — no preceding leg; for a closed cycler this would be the previous lap's V∞_out, which equals V∞_out_0 in magnitude). Set `vinf_in = vinf_out` for the first encounter in the open-sequence case so `maintenance_dv` doesn't spuriously charge it. Document this convention in the docstring.
   - Intermediate encounters `i`: `vinf_in = v_arrive_{i-1} - v_planet_i`, `vinf_out = v_depart_i - v_planet_i`.
   - Last encounter: `vinf_in = v_arrive_{last-1} - v_planet_last`, `vinf_out = vinf_in` (same convention, opposite end).
6. Build `Leg` objects with `from_body`, `to_body`, `t_depart`, `t_arrive`, `v_depart`, `v_arrive`, `n_revs`, `branch`.
7. Period = `encounter_times_sec[-1] - encounter_times_sec[0]`.
8. Return `Cycler(bodies=sequence, period=period, encounters=encounters, legs=legs)`.

#### Aldrin seed details

The published Aldrin numbers (spec §9) use Earth–Mars *circular* heliocentric orbits with semi-major axes 1.0 AU and 1.524 AU. The M1 ephemeris uses the constants from `PLANETS["E"].sma_au` (1.00000261 AU) and `PLANETS["M"].sma_au` (1.52371034 AU) — close to the canonical literature numbers but not identical. We accept this small mismatch; the test tolerances (§4) absorb it.

For `t_start_sec = 0`, we place Earth at the standard `(a_E, 0, 0)` and Mars at `t_start + 146 days` at its circular-coplanar position (Mars's angular position at that moment depends on the M1 ephemeris's epoch convention — Mars at t=0 sits somewhere on its orbit, propagated forward by 146 days at `mean_motion`). The test does not care about absolute longitudes; it only cares about the orbital elements of the resulting transfer.

For the 2-synodic E-M-E test, encounter times are `[0, 146 d, synodic_period_EM]` (Earth → Mars → Earth, where the second Earth encounter is exactly one synodic period after the first). This gives a single Lambert call for each leg; the second leg returns Earth → Earth across (synodic - 146 d), which is a long transit and may require `n_revs ≥ 1` (TBD in the test — the test first attempts `n_revs=0`, falls back to `n_revs=1` if no solution converges; this fallback path validates the multi-rev branch selection from spec §12(b)).

---

### 3.4 Imports / dependency graph after M3

```
constants.py    (M0)  ← root, no deps
ephemeris.py    (M1)  ← constants
lambert.py      (M1)  ← constants
kepler.py       (M1)  ← constants
flyby.py        (M2)  ← constants
frames.py       (M3)  ← constants                 [new]
tisserand.py    (M2)  ← constants, flyby
resonance.py    (M2)  ← constants
model/cycler.py (M3)  ← constants, frames         [new]
construct.py    (M3)  ← constants, ephemeris,     [new]
                       lambert, model/cycler
```

No cycles; `model/` depends on `core/`; `search/` depends on `core/` and `model/`; nothing in `core/` depends on `model/` or `search/`.

---

## 4. Tests + gate

Three test files, mirroring the three new source modules. Tolerances are **named module-level constants** at the top of each test file so loosening or tightening is a one-line change.

### 4.1 `tests/core/test_frames.py`

Risk-focused. Spec §10 calls out "closure frame correctness — unit-test the rotating-frame transform; a wrong frame silently fakes/breaks closure." This file is that test.

```python
TOL_ROUND_TRIP_KM = 1e-10        # round-trip identity for (r, v) at 1 AU scale
TOL_ROUND_TRIP_KMS = 1e-13       # km/s

def test_round_trip_identity():
    # Random states at 1 AU scale, various t and omega.
    # for each: r2, v2 = from_rotating(*to_rotating(r, v, t, omega), t, omega)
    # assert ‖r2 - r‖ < TOL_ROUND_TRIP_KM and ‖v2 - v‖ < TOL_ROUND_TRIP_KMS.

def test_omega_zero_is_identity():
    # to_rotating(r, v, t, 0.0) == (r, v) for any t.

def test_circular_orbit_is_stationary_in_its_own_frame():
    # A particle on a circular orbit at radius a with v_circ = sqrt(mu/a),
    # transformed into the frame rotating at its own mean motion, has v_rot ≈ 0.

def test_synodic_omega_earth_matches_constants():
    # synodic_omega("E") == PLANETS["E"].mean_motion_deg_day converted to rad/s,
    # within 1e-12 rel tolerance.

def test_synodic_omega_venus_raises():
    # synodic_omega("V") is M8 work; raises NotImplementedError in M3.
```

The round-trip test is the load-bearing one — it directly addresses the spec §10 risk. The "circular orbit stationary" test is a physics sanity that catches the most common sign-error bug (the `ω×r` term having the wrong sign would silently survive the round-trip but blow up here).

### 4.2 `tests/model/test_cycler.py`

Dataclass behaviour + closure-residual semantics. No physics; just structural correctness.

```python
TOL_FROZEN = "dataclasses.FrozenInstanceError"   # symbolic — assert raises this

def test_encounter_is_frozen():
    e = Encounter(...)
    with pytest.raises(FrozenInstanceError):
        e.body = "M"

def test_leg_defaults():
    leg = Leg(from_body="E", to_body="M", ...)   # no n_revs, no branch
    assert leg.n_revs == 0
    assert leg.branch == "single"

def test_cycler_max_vinf_returns_largest_magnitude():
    # Construct a Cycler with hand-crafted Encounters of known vinf magnitudes.
    # Assert max_vinf() returns the maximum.

def test_cycler_maintenance_dv_zero_when_vinf_in_equals_vinf_out():
    # Hand-crafted Cycler where every encounter has vinf_in == vinf_out (ballistic).
    # Assert maintenance_dv() < 1e-12 km/s.

def test_cycler_radial_span_two_body_circular():
    # Hand-crafted Cycler whose single leg is a known Hohmann transfer.
    # Assert radial_span() matches the Hohmann (departure, arrival) radii to 1e-6 AU.

def test_closure_residual_uses_default_earth_omega():
    # closure_residual(omega=None) should match closure_residual(omega=synodic_omega("E")).
```

### 4.3 `tests/model/test_aldrin.py` — the gate

```python
# Tolerances (spec §9 anchors)
TOL_A_AU       = 0.01    # |a - 1.659| < 0.01
TOL_E          = 0.02    # |e - 0.41| < 0.02
TOL_PERI_AU    = 0.02    # |perihelion - 0.98| < 0.02
TOL_APO_AU     = 0.02    # |aphelion - 2.34| < 0.02
TOL_TOF_DAYS   = 1.0     # leg time-of-flight within ±1 d of 146

def test_aldrin_E_to_M_leg_elements():
    ephem = Ephemeris(model="circular")
    cyc = build_aldrin_seed(ephem, t_start_sec=0.0)
    leg = cyc.legs[0]
    r_dep = cyc.encounters[0].r
    v_dep = leg.v_depart
    a_au, e = orbital_elements_from_rv(r_dep, v_dep, mu=MU_SUN_KM3_S2)
    peri_au = a_au * (1 - e)
    apo_au  = a_au * (1 + e)

    assert abs(a_au - 1.659) < TOL_A_AU
    assert abs(e    - 0.41 ) < TOL_E
    assert abs(peri_au - 0.98) < TOL_PERI_AU
    assert abs(apo_au  - 2.34) < TOL_APO_AU

    tof_days = (leg.t_arrive - leg.t_depart) / SECONDS_PER_DAY
    assert abs(tof_days - 146.0) < TOL_TOF_DAYS

def test_aldrin_closure_residual_small():
    ephem = Ephemeris(model="circular")
    cyc = build_aldrin_seed(ephem, t_start_sec=0.0)
    # The 1-synodic Aldrin "slice" (E->M) is not a closed cycler on its own;
    # this test instead exercises the closure path on the 2-synodic E-M-E
    # construction via test_construct.py. Here we just assert that
    # closure_residual() is callable and returns a finite float >= 0.
    r = cyc.closure_residual()
    assert r >= 0.0
    assert math.isfinite(r)
```

`orbital_elements_from_rv` is a 10-line helper inside the test file (not promoted to a module — it's M5 work to add a proper orbital-elements module). It computes `(a, e)` from the vis-viva equation and the eccentricity vector.

### 4.4 `tests/model/test_construct.py` — 2-synodic E–M–E

```python
TOL_CLOSURE_KMS = 0.05      # 50 m/s residual is comfortably tight for circular-coplanar

def test_construct_two_synodic_em_cycler():
    ephem = Ephemeris(model="circular")
    T_syn = synodic_period_seconds("E", "M")     # from M2 resonance
    times = [0.0, 146.0 * SECONDS_PER_DAY, T_syn]
    cyc = construct_cycler(
        sequence=["E", "M", "E"],
        encounter_times_sec=times,
        ephem=ephem,
    )
    assert len(cyc.encounters) == 3
    assert len(cyc.legs) == 2
    assert cyc.period == pytest.approx(T_syn)

    # Closure residual: small but not zero (only Lambert numerical error contributes).
    assert cyc.closure_residual() < TOL_CLOSURE_KMS

def test_construct_validates_input_lengths():
    ephem = Ephemeris(model="circular")
    with pytest.raises(ValueError):
        construct_cycler(["E", "M"], [0.0], ephem)        # length mismatch
    with pytest.raises(ValueError):
        construct_cycler(["E", "M"], [100.0, 0.0], ephem) # non-monotonic

def test_construct_unknown_body_raises():
    ephem = Ephemeris(model="circular")
    with pytest.raises(ValueError):
        construct_cycler(["E", "X"], [0.0, 1e6], ephem)

def test_construct_multi_rev_branch_selection():
    # Long E->E transit where n_revs=0 has no solution; verify n_revs=1, branch="low"
    # returns a valid Leg with the expected branch label.
    ...
```

The `test_construct_two_synodic_em_cycler` test is the *second half* of the spec §8 gate ("reproduce ... a 2-synodic E–M cycler"). The `test_construct_multi_rev_branch_selection` test is the M3 entry point for the multi-rev Lambert branching from spec §12(b) — it does not exercise optimisation over branch choice (M4 work), only that the constructor respects an explicit branch argument.

### 4.5 Gate restated

The M3 gate passes when **all** of the following are true:

- `tests/core/test_frames.py::test_round_trip_identity` passes at `TOL_ROUND_TRIP_KM = 1e-10`.
- `tests/model/test_aldrin.py::test_aldrin_E_to_M_leg_elements` passes at the spec §9 tolerances above.
- `tests/model/test_construct.py::test_construct_two_synodic_em_cycler` shows a closure residual `< 0.05 km/s` (numerically-limited only; the circular-coplanar geometry is exact).
- `uv run pytest`, `uv run ruff check .`, and `uv run mypy src tests` are all clean (carrying forward the M0 standard).

---

## 5. Risks

| # | Risk | Headline | Mitigation in M3 |
|---|------|----------|------------------|
| 1 | **Closure frame correctness** (spec §10) | A wrong sign or wrong-axis rotation matrix silently makes closure "work" or "fail" by the wrong amount. This is the single most likely source of an undetected bug. | Dedicated `test_round_trip_identity` at `1e-10 km` tolerance, plus `test_circular_orbit_is_stationary_in_its_own_frame` as an independent physical anchor (the sign error survives the round-trip but blows up here). |
| 2 | **Mis-identifying the frame's body for ω** | Using Mars's mean motion instead of Earth's silently produces a wrong-but-plausible residual. | `synodic_omega("E")` is the documented default; `Cycler.closure_residual(omega=None)` uses it explicitly. Test `test_closure_residual_uses_default_earth_omega` pins this. |
| 3 | **Aldrin elements depend on choice of semi-major axes** | Literature numbers use 1.0 and 1.524 AU exactly; M1's `PLANETS["E"].sma_au = 1.00000261` and `PLANETS["M"].sma_au = 1.52371034` shift the result by ~0.001 AU. | Test tolerances (`TOL_A_AU = 0.01`) absorb this; M6 will document the shift explicitly. |
| 4 | **V∞ convention at first/last encounter of open sequence** | If `vinf_in` on encounter 0 is mistakenly set to `(0,0,0)` (no preceding leg), `maintenance_dv()` will spuriously charge `‖vinf_out‖` at the first node. | The constructor sets `vinf_in = vinf_out` on the boundary encounters of an open sequence (documented in the `construct_cycler` docstring). `test_cycler_maintenance_dv_zero_when_vinf_in_equals_vinf_out` covers the dataclass; the Aldrin construction implicitly relies on it. |
| 5 | **Multi-rev Lambert branch selection** (spec §10) | The 2-synodic E–E leg may require `n_revs=1`; picking the wrong branch silently gives a different orbit. | M3 only respects an explicit `branch_per_leg` argument; it does not search. `test_construct_multi_rev_branch_selection` verifies the constructor passes the branch through correctly. The optimisation across branches is M4/M5. |
| 6 | **Conflating "uniform frame" with "dynamic frame"** (spec §12(c)) | A reviewer reading `closure_residual` and assuming it works for real ephemeris would be silently wrong. | The `frames.py` module docstring explicitly states the uniform-frame scope. M6 plan will introduce a `Frame` abstraction (uniform vs dynamic) and rewire `closure_residual` to take a `Frame` rather than an `omega` scalar. |
| 7 | **Frozen dataclass + numpy arrays make equality awkward** | Two `Cycler` instances with identical arrays compare unequal under `==` because numpy `array_equal` is not used by the dataclass-generated `__eq__`. | Documented (no `__eq__` override in M3). Tests use field-by-field comparison; identity matching is the catalogue's job in M7 (spec §16.2). |

---

## 6. Dependency additions

**None.** M3 uses only `numpy` (M0) and the M1/M2 modules. No new entries in `pyproject.toml`; no `uv lock` regeneration needed.

---

## 7. Order of work

The todo.md mirrors this with checkboxes.

1. **Frame transforms first.** `core/frames.py` + `tests/core/test_frames.py`. The round-trip test is the M3 risk gate; landing it first means every later module is built on a frame we trust.
2. **Dataclasses next.** `model/cycler.py` + `tests/model/test_cycler.py`. Structural tests only — no Lambert calls. Defines the contract `construct.py` will fulfil.
3. **Constructor and gate.** `search/construct.py` + `tests/model/test_aldrin.py` + `tests/model/test_construct.py`. Build `construct_cycler` first (it's the general routine); `build_aldrin_seed` is a one-screen wrapper over it.
4. **Re-export aggregator.** `model/__init__.py` re-exporting `Cycler`, `Leg`, `Encounter`. Trivial; landed last so the public surface is settled.
5. **Run the gate.** `uv run pytest tests/model/test_aldrin.py -v` — confirm spec §9 values reproduce.
6. **Lint and type-check sweep.** `uv run ruff check .`, `uv run ruff format --check .`, `uv run mypy src tests`. The `NDArray[np.float64]` annotations on dataclass fields will be the most likely source of mypy surprises; resolve in module, not via `# type: ignore`.
7. **Update milestone status.** `docs/overview.md` table: M3 status `planned` → `completed`; M4 row `not yet planned` → `planned`.
8. **Hand-off note.** Append `## Hand-off to M4` to `todo.md` noting anything M4 needs to know (e.g. which Lambert branch the E-M-E test actually used, whether the Aldrin tolerances were tight or comfortable).

The order is "frames → cycler → construct" deliberately: each module imports only what came before it, so any test failure localises immediately to the module just landed.

---

## 8. Exit checklist (the gate, restated)

Before declaring M3 done:

- [ ] `uv run pytest tests/core/test_frames.py` green; `test_round_trip_identity` passes at `1e-10 km`.
- [ ] `uv run pytest tests/model/test_cycler.py` green; frozen-ness, defaults, and metric semantics all verified.
- [ ] `uv run pytest tests/model/test_aldrin.py` green; `a, e, perihelion, aphelion, tof` all within spec §9 tolerances.
- [ ] `uv run pytest tests/model/test_construct.py` green; 2-synodic E–M–E built, closure residual `< 0.05 km/s`.
- [ ] `uv run ruff check .` clean.
- [ ] `uv run ruff format --check .` clean.
- [ ] `uv run mypy src tests` clean (no new `# type: ignore` introduced).
- [ ] CI on the M3 branch is green.
- [ ] `docs/overview.md` updated: M3 status = `completed`; M4 row = `planned`.
- [ ] `## Hand-off to M4` section appended to `phases/m3-model-construct/todo.md` with anything M4 needs to know (branch choices, residual values observed, tolerance margins).

(Writing the M4 plan doc is the first task of M4, not an M3 exit criterion.)
