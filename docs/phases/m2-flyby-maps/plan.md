# M2 — Flyby + maps

**Spec reference:** spec.md §3 (background — synodic periods, V∞, bend formula), §4 (architecture — `core/flyby.py`, `search/tisserand.py`, `search/resonance.py`), §5 (algorithm pipeline steps 1–2), §6 (interfaces sketch — `flyby.py`), §7 (tech stack), §8 (M2 milestone definition), §9 (validation anchors), §13.3 (Tisserand pruning — the role of `linkable`).

**Purpose:** deliver the **energetics layer** that gates everything M3+ will build on. M2 produces the gravity-assist mechanics (`flyby`), the V∞-contour graph (`tisserand`), and the resonance arithmetic (`resonance`) needed to (a) compute per-encounter maximum bend angles, (b) decide whether two bodies are linkable at a given V∞ (the M4 pruning predicate, spec §13.3), and (c) generate the candidate cycler periods (k·synodic and multi-body beats) that seed M3 construction.

**Gate (definition of done):** the four spec §9 anchors that fall in M2 all pass as pytest assertions:

- Earth–Mars synodic period = 2.135 yr (tol 0.001 yr).
- Earth–Venus synodic period = 1.599 yr (tol 0.001 yr).
- VEM beat: `multi_body_beat_days(["V","E","M"], k_max=6)` yields a `(3, 4)` (or equivalent ordering) commensurability ≈ 6.406 yr (tol 0.01 yr).
- Mars max bend at V∞ = 7 km/s, `rp_min = R_Mars + 300 km` ≈ 24° (tol 1°).

Plus the sanity check that Earth/Venus max bend at V∞ = 7 km/s lies in [60°, 63°], and that `flyby_dv` is exactly zero on a ballistic-feasible pair and strictly positive on an over-bent pair.

---

## 1. What this milestone delivers

Three modules and their tests:

1. **`src/cyclerfinder/core/flyby.py`** — gravity-assist mechanics: max bend, ballistic-feasibility predicate, powered-flyby ΔV, and a planet-aware convenience wrapper. Pure functions; no I/O.
2. **`src/cyclerfinder/search/tisserand.py`** — V∞-graph contours in (a, e) space, the pairwise `linkable(body_a, body_b, vinf)` predicate (the **M4 pruning gate**, spec §13.3), and a `linkable_region` interval helper. Includes an optional plotting helper guarded behind the `viz` extra.
3. **`src/cyclerfinder/search/resonance.py`** — synodic periods, k-multiples, and a multi-body-beat finder that returns integer commensurability tuples within a tolerance.

M2 also introduces the `search/` subpackage (its first occupants) and the optional `viz` dependency group (`matplotlib`), per the *what exists = what works* convention from M0.

**Explicitly out of scope for M2** (later milestones — do not implement here):

- `search/sequence.py` — enumeration of body orderings × revolution counts × branches. M4.
- `search/construct.py` — Lambert-built leg construction from a (sequence, times) tuple. M3 (single-sequence form) and M4 (over enumeration output).
- `model/cycler.py` — the `Cycler` / `Leg` / `Encounter` dataclasses and closure residual. M3.
- The optimiser, the scoring module, the ephemeris bridge. M5–M6.
- 3D Tisserand with non-zero inclination. M6+ (see §3.2 below for the explicit coplanar assumption).

---

## 2. File tree after M2

```
cyclers/
├── … (M0/M1 layout preserved)
├── src/
│   └── cyclerfinder/
│       ├── __init__.py
│       ├── core/
│       │   ├── __init__.py
│       │   ├── constants.py            # M0
│       │   ├── ephemeris.py            # M1
│       │   ├── lambert.py              # M1
│       │   ├── kepler.py               # M1
│       │   └── flyby.py                # NEW (M2)
│       └── search/                     # NEW subpackage (M2)
│           ├── __init__.py             # empty
│           ├── tisserand.py            # NEW (M2)
│           └── resonance.py            # NEW (M2)
└── tests/
    ├── … (M0/M1 tests preserved)
    ├── test_flyby.py                   # NEW
    ├── test_tisserand.py               # NEW
    └── test_resonance.py               # NEW
```

Subpackages `model/`, `verify/`, `data/`, `viz/` are still **not** created — they appear when their first module lands. The `viz/` subpackage in particular waits until M8; M2's plotting helper lives inside `tisserand.py` itself, gated by an import-guard on `matplotlib`.

---

## 3. Module designs

### 3.1 `core/flyby.py`

Pure gravity-assist mechanics on heliocentric V∞ vectors. The full spec §6 interface plus a planet-aware convenience wrapper that reads `constants.PLANETS` and `constants.SAFE_PERIHELION_KM` (introduced in M0).

#### Signatures

```python
from numpy.typing import NDArray
import numpy as np

def max_bend(mu_planet: float, rp_min: float, vinf: float) -> float:
    """Max ballistic deflection angle (rad).
    sin(delta/2) = 1 / (1 + rp_min * vinf**2 / mu_planet).
    Returns pi as vinf -> 0; returns ~0 as vinf -> infinity."""

def bend_angle(vin_vec: NDArray[np.float64], vout_vec: NDArray[np.float64]) -> float:
    """Angle (rad) between two V∞ vectors of (ideally) equal magnitude.
    Numerically robust via clip on the acos argument."""

def is_ballistic_feasible(
    vin_vec: NDArray[np.float64],
    vout_vec: NDArray[np.float64],
    mu_planet: float,
    rp_min: float,
    speed_tol: float = 1e-6,
) -> bool:
    """True iff |vin| == |vout| within speed_tol (km/s) AND
    bend_angle(vin,vout) <= max_bend(mu_planet, rp_min, |vin|)."""

def flyby_dv(
    vin_vec: NDArray[np.float64],
    vout_vec: NDArray[np.float64],
    mu_planet: float,
    rp_min: float,
) -> float:
    """km/s. Powered-flyby ΔV required to convert vin into vout.

    Decomposition (matches spec §6 sketch):
      * magnitude component: |vout| - |vin| (may be sign-mixed); the
        scalar cost is the impulse magnitude required at the periapsis
        burn to change the asymptotic speed by that amount.
      * direction component: any bend in excess of max_bend(|v_burn|, ...)
        is paid as a second tangential impulse rotating the asymptote
        until the residual bend is within the achievable cone.
    Returns 0.0 exactly when is_ballistic_feasible is True."""

def flyby_dv_for(
    code: str,
    vin_vec: NDArray[np.float64],
    vout_vec: NDArray[np.float64],
) -> float:
    """Convenience wrapper: looks up mu_planet from PLANETS[code]
    and rp_min from SAFE_PERIHELION_KM[code]."""
```

#### Design notes

- **Single-impulse decomposition.** The full powered-flyby is an open optimisation; we adopt the standard patched-conic surrogate (a single tangential burn at periapsis that adjusts hyperbolic excess magnitude, plus a residual bend-deficit penalty if the remaining geometry still over-bends). This is **good enough for M4 scoring** and matches the spec §6 sketch (`flyby_dv` is zero in the ballistic case, positive otherwise). M5 may refine this if optimiser behaviour demands it; we flag this as a model-fidelity choice rather than a bug.
- **Floating-point near-zero.** `vinf` near zero causes `sin(delta/2) -> 1` (bend approaches π); guard with a small `vinf > 0` check that returns `np.pi` for `vinf == 0`. Equal-magnitude tolerance is in km/s, with default `1e-6` matching machine-precision Lambert outputs from M1.
- **Vector contract.** All `*_vec` arguments are 1-D `NDArray[np.float64]` of length 3. mypy-strict means we declare those via `numpy.typing.NDArray[np.float64]`. Scalar inputs/returns are plain `float`, never `np.float64`, per the M0 convention.
- **No singleton planet lookup in the core helpers** — only `flyby_dv_for` reaches into `constants.PLANETS`. Keeps the math primitive testable without the planet table.

#### Validation hooks (used by tests, spec §9)

| Anchor | Implementation check |
|--------|----------------------|
| Mars max bend @ V∞ = 7 km/s ≈ 24° | `np.degrees(max_bend(PLANETS["M"].mu_km3_s2, SAFE_PERIHELION_KM["M"], 7.0))` ∈ [23, 25] |
| Earth max bend @ V∞ = 7 km/s ≈ 60–63° | range assertion ∈ [60, 63] |
| Venus max bend @ V∞ = 7 km/s ≈ 60–63° | range assertion ∈ [60, 63] |
| Ballistic-feasible pair → ΔV = 0 | hand-construct vin/vout with equal magnitude and bend < max_bend |
| Over-bent pair → ΔV > 0 | rotate vout past max_bend; assert returned value > 0 and finite |

### 3.2 `search/tisserand.py`

V∞-graph contours plus the **`linkable`** predicate that spec §13.3 identifies as the linchpin of the M4 pruning gate.

#### The coplanar reduction (mandatory scope note)

The Tisserand parameter at body p (semi-major axis `a_p`) for a spacecraft heliocentric orbit `(a, e, i)` is

```
T_p(a, e, i) = a_p/a + 2 * cos(i) * sqrt((a/a_p) * (1 - e**2))
```

**M2 implements the coplanar case only**, fixing `i = 0` so `cos(i) = 1` throughout. The 3D Tisserand (`i != 0`, with its extra DOF) is M6+ work, deferred per the spec §12 ephemeris-realisation pipeline (where inclination is absorbed into the b-plane during phase-matching, not the Tisserand cell-enumeration). This restriction must be **explicit in the module docstring**, in the function docstrings, and re-stated in the M2 hand-off note to M3.

Why coplanar is sufficient for M2: cell enumeration (M4) walks (body-set, sequence, k, revs, branches) over the coplanar V∞ graph; the M3 Aldrin reproduction is itself coplanar; the VEM enumeration (M8) uses the coplanar graph as the structure generator and then phase-matches inclined ephemeris in M6.

#### The contour parameterisation

For a given body `p` (so `a_p` is fixed) and a target Tisserand value `T_target` (equivalently, a target V∞: `V∞² = μ_sun/a_p · (3 - T_target)`), the constant-`T_p` curve in (a, e) space is **the** parameterisation choice. Two natural options:

- **Parameterise by `e`** (chosen): for each `e ∈ [0, 1)`, solve the algebraic equation `a_p/a + 2 * sqrt((a/a_p) * (1 - e**2)) = T_target` for `a`. This is a quartic-like equation that reduces to a quadratic in `sqrt(a/a_p)`:

  Let `u = sqrt(a/a_p)`. Then `1/u² + 2u * sqrt(1 - e²) = T_target` ⇒ `1 + 2u³ * sqrt(1 - e²) = T_target * u²`. Cubic in `u`. We solve numerically via `scipy.optimize.brentq` on a bracketed root (the physically meaningful branch with `u > 0` and the resulting `a` such that perihelion `a(1-e) ≤ a_aphelion_other_body` constraints are not enforced here — they're enforced inside `linkable`).

- Parameterise by `a` (rejected for the primary path): the curve doubles back on itself in `e` for some values of `a`, making sample-by-`a` lossy. The `e`-parameterised sweep is single-valued in `a` per branch.

We pick the `e`-parameterised sweep; `n_points` samples of `e` from `0` to `1 − eps` (default `n_points = 200`, `eps = 1e-4`). Sample points where no real `a > 0` solution exists are dropped (filtered out, **not** raised — this is a normal regime, not an error).

#### Signatures

```python
from numpy.typing import NDArray
import numpy as np

def vinf_to_tisserand(body: str, vinf_kms: float) -> float:
    """T_p = 3 - V∞² * a_p / μ_sun. Uses constants.PLANETS[body].sma_au
    converted to km and constants.MU_SUN_KM3_S2."""

def tisserand_to_vinf(body: str, T_p: float) -> float:
    """Inverse. Returns 0.0 if T_p >= 3 (no real V∞)."""

def vinf_contour(
    body: str,
    vinf_kms: float,
    a_range_au: tuple[float, float] = (0.3, 5.0),
    n_points: int = 200,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Constant-V∞ contour at `body` in (a [AU], e) space, parameterised
    by e ∈ [0, 1 - 1e-4) with n_points samples. Returns two equal-length
    arrays; entries where no a in a_range_au solves the Tisserand equation
    are omitted (arrays may be shorter than n_points). Empty arrays if
    no contour exists at this V∞."""

def linkable(
    body_a: str,
    body_b: str,
    vinf_kms: float,
    tol_au: float = 0.01,
    tol_e: float = 0.01,
) -> bool:
    """True iff the constant-V∞ contours of body_a and body_b intersect
    in (a, e) space within (tol_au, tol_e). Intersection means: there
    exists an orbit reaching both bodies at this V∞ — hence a flyby
    sequence between them at this V∞ is energetically possible.

    Implementation: for each body, the contour `e -> a(e)` is a 1-D
    curve; intersection is the zero-set of `a_a(e) - a_b(e)` on shared
    e-support. Use scipy.optimize.brentq on a sign change of the
    difference function.

    This predicate is the M4 Tisserand-pruning gate (spec §13.3)."""

def linkable_region(
    body_a: str,
    body_b: str,
    vinf_cap_kms: float,
    n_vinf: int = 50,
) -> list[float]:
    """Sample V∞ ∈ (0, vinf_cap_kms] and return the V∞ values where
    linkable(body_a, body_b, V∞) is True. Returned as a list (possibly
    empty); M4 will reduce this to interval edges via run-length encoding.

    M2 returns the raw True-set sampled at n_vinf points; the M4 caller
    can post-process. The reason for the raw form here is that M2 has
    no use for the intervals yet — exposing them prematurely would
    couple M2 to an M4-shaped API."""

def plot_tisserand(
    bodies: list[str],
    vinf_levels_kms: list[float],
    ax: "Axes | None" = None,  # matplotlib.axes.Axes; optional dep
) -> "Axes":
    """Diagnostic plot: overlays contours for every (body, V∞) pair.
    Imports matplotlib lazily inside the function; raises ImportError
    with a hint to install the `viz` extra if matplotlib is absent."""
```

#### What `linkable` returning False means (spec §13.3, called out by request)

If `linkable(body_a, body_b, V∞) == False`, then **no orbit at that V∞ reaches both bodies**, so **no flyby sequence between them at that V∞ can exist**, so the M4 enumerator can discard every cell whose consecutive bodies require linking at that V∞. This is the structural pruning that makes the M4 search finite-effort: the vast majority of (body-set, sequence, k, revs, branches) cells die at this gate without ever entering optimisation. The correctness of this predicate is therefore the linchpin of search tractability. Its tests are correspondingly tight.

#### Numerical edge cases

| Case | Handling |
|------|----------|
| `V∞` too small for a real contour at `body` | `vinf_contour` returns empty arrays; `linkable` returns False; never raises. |
| `e -> 1` (parabolic limit) | Sweep stops at `1 - 1e-4`; the contour family approaches it asymptotically. |
| `a -> 0` or `a > a_range_au[1]` | Root-finder rejects the sample (out of bracket); contour just doesn't extend there. |
| Both contours empty | `linkable` returns False (vacuous). |
| Contours touch tangentially | Counted as linkable; `brentq` sign-change check uses tolerance `tol_e`. |

### 3.3 `search/resonance.py`

Synodic periods and multi-body beats. All times in days; year conversions use `constants.DAYS_PER_JULIAN_YEAR` (365.25, per M0).

#### Signatures

```python
def synodic_period_days(body_a: str, body_b: str) -> float:
    """1 / |1/T_a - 1/T_b| where T_x = 360 / constants.PLANETS[x].mean_motion_deg_day.
    Raises ValueError on body_a == body_b."""

def synodic_period_years(body_a: str, body_b: str) -> float:
    """synodic_period_days / DAYS_PER_JULIAN_YEAR. Convenience."""

def k_synodic_periods_days(body_a: str, body_b: str, k_max: int) -> list[float]:
    """[k * synodic for k in range(1, k_max + 1)]. The candidate cycler
    period bank for a single body pair."""

def multi_body_beat_days(
    bodies: list[str],
    k_max: int = 6,
    tol_frac: float = 0.02,
) -> list[tuple[int, ...]]:
    """For a body set of size >= 2, find integer tuples
    (k_1, k_2, ...) with each k_i in 1..k_max where
    k_1 * synodic(body_1, ref) ≈ k_2 * synodic(body_2, ref) ≈ ...
    within fractional tolerance tol_frac.

    The 'ref' body is chosen as the body shared across all pairs
    (Earth is the natural reference for V/E/M sets). Implementation:
      * pick ref = bodies[1] for a 3-body set (E in [V,E,M]);
      * compute T_syn(V,E) and T_syn(E,M);
      * search k_V, k_M in 1..k_max for the closest match.

    For [V, E, M] returns [(3, 4)] (or its equivalent ordering) ⇒
    3 * T_syn(E,M) ≈ 4 * T_syn(E,V) ≈ 6.406 yr. Ranked by mismatch
    ascending; usually a single dominant tuple but multiple are kept
    for transparency."""

def beat_period_days(bodies: list[str], k_tuple: tuple[int, ...]) -> float:
    """Reverse helper: given a k-tuple from multi_body_beat_days,
    return the mean of k_i * T_syn(body_i, ref) — i.e. the beat
    period in days. The 'period' a returned tuple corresponds to."""
```

#### Design notes

- **Reference-body choice for >2 bodies.** For a 3-body set we treat the **middle body** as the reference (so `[V, E, M]` uses E as reference; `T_syn(V,E)` and `T_syn(E,M)` are the two co-resonances). For more than 3 bodies the rule becomes "the body that maximises pairwise coverage"; M2 only validates the 3-body case and documents the rule, deferring N-body to M8.
- **Tolerance.** `tol_frac = 0.02` (2 %) accommodates the small irrational mismatch between 3·E–M and 4·E–V (the actual ratio is ≈ 1.00097, so a 0.1 % tolerance would also work). The looser default lets near-miss tuples surface for inspection; M4 may tighten this when filtering cells.
- **Period bank for M3/M4.** `k_synodic_periods_days` is what M3 calls to get the Aldrin candidate period (k=1 on E-M ≈ 2.135 yr) and M4 calls to enumerate candidate periods per cell.
- **No frequency analysis, no FFT, no continued fractions.** This is a finite integer search over `(k_1, ..., k_N) ∈ [1, k_max]^N`; for `N = 3, k_max = 6` that's 36 tuples — trivially exhaustive. Avoid overengineering.

---

## 4. Tests + gate

All four modules get a dedicated test file. The gate assertions are explicit; the supporting tests cover edge cases and the contracts that M3/M4 will rely on.

### 4.1 `tests/test_flyby.py`

| Test | Assertion | Tolerance |
|------|-----------|-----------|
| `test_mars_max_bend_24deg_at_7kms` | `np.degrees(max_bend(MU_MARS, RP_MARS_SAFE, 7.0)) == pytest.approx(24, abs=1)` | ±1° (spec §9) |
| `test_earth_max_bend_in_range_at_7kms` | bend ∈ [60, 63]° | sanity range (spec §9) |
| `test_venus_max_bend_in_range_at_7kms` | bend ∈ [60, 63]° | sanity range (spec §9) |
| `test_max_bend_zero_at_infinite_vinf` | `max_bend(..., vinf=1e6)` < 1e-4 rad | asymptote |
| `test_max_bend_pi_at_zero_vinf` | `max_bend(..., vinf=0.0) == pytest.approx(np.pi, abs=1e-10)` | guard branch |
| `test_ballistic_feasible_zero_dv` | hand-built equal-mag pair with bend = 0.5 * max_bend → `flyby_dv == 0.0` | exact |
| `test_overbent_pair_positive_dv` | same magnitude, bend = 1.5 * max_bend → `flyby_dv > 0` and finite | strict |
| `test_speed_mismatch_positive_dv` | `|vin| = 5, |vout| = 7` with bend 0 → `flyby_dv >= 2.0 - eps` | lower bound |
| `test_flyby_dv_for_matches_explicit` | `flyby_dv_for("M", vin, vout) == flyby_dv(vin, vout, MU_MARS, RP_MARS_SAFE)` | exact |
| `test_is_ballistic_feasible_consistency` | for 50 random pairs, `(flyby_dv == 0) iff is_ballistic_feasible` | parameterised |

### 4.2 `tests/test_resonance.py`

| Test | Assertion | Tolerance |
|------|-----------|-----------|
| `test_em_synodic_2135yr` | `synodic_period_years("E","M") == pytest.approx(2.135, abs=0.001)` | spec §9 gate |
| `test_ev_synodic_1599yr` | `synodic_period_years("E","V") == pytest.approx(1.599, abs=0.001)` | spec §9 gate |
| `test_synodic_symmetric` | `T_syn(a,b) == T_syn(b,a)` | exact |
| `test_synodic_self_raises` | `synodic_period_days("E","E")` raises `ValueError` | API contract |
| `test_vem_beat_yields_3_4` | `(3, 4)` (or equivalent ordering) ∈ `multi_body_beat_days(["V","E","M"], k_max=6)` | spec §9 gate |
| `test_vem_beat_period_6_406yr` | `beat_period_days(["V","E","M"], (3, 4)) / 365.25 == pytest.approx(6.406, abs=0.01)` | spec §9 gate |
| `test_k_synodic_monotone` | `k_synodic_periods_days("E","M",5)` strictly increasing | sanity |

### 4.3 `tests/test_tisserand.py`

| Test | Assertion | Tolerance |
|------|-----------|-----------|
| `test_vinf_to_tisserand_inverse` | `tisserand_to_vinf("E", vinf_to_tisserand("E", 5.0)) == pytest.approx(5.0, abs=1e-9)` | round-trip |
| `test_contour_returns_arrays_or_empty` | `vinf_contour("E", 0.01)` returns two empty arrays (V∞ too small) | edge case |
| `test_contour_non_empty_at_reasonable_vinf` | `vinf_contour("E", 5.0)` returns ≥ 50 points within default a_range | smoke |
| `test_contour_a_at_e0_matches_homann` | for `e=0` (circular orbit), the contour passes through `a == a_p_body` only at V∞ = 0 | derived |
| `test_linkable_em_at_5_5_kms_true` | `linkable("E","M", 5.5)` is True (Aldrin family lives here, spec §9) | physics |
| `test_linkable_em_at_0_5_kms_false` | `linkable("E","M", 0.5)` is False (too little energy to span) | physics |
| `test_linkable_symmetric` | `linkable(a,b,V∞) == linkable(b,a,V∞)` for several samples | API contract |
| `test_linkable_region_em_non_empty` | `len(linkable_region("E","M", 12.0))` > 0 | smoke |
| `test_linkable_never_raises` | parametrised over impossible V∞ (0.0, 1e-6, 1e6) → False, no exception | robustness |

### 4.4 Plotting test (viz extra, conditional)

`tests/test_tisserand.py::test_plot_runs_when_matplotlib_present` is decorated with `@pytest.mark.skipif(matplotlib is None, ...)`. Asserts only that the call returns an `Axes` without raising — the figure content is a human-inspection artefact, not a CI gate. Marked `viz` so it can be filtered.

### 4.5 What CI runs

The standard set from M0 (`ruff check`, `ruff format --check`, `mypy src tests`, `pytest`). The viz tests run only if `uv sync --all-extras` includes the `viz` group (it does — `--all-extras` pulls everything declared). On a fresh checkout without extras the skip-condition fires.

---

## 5. Risks

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| **Contour root-finder fails to bracket** for an extreme V∞ | medium | medium — would raise instead of returning empty | Use try/except around `scipy.optimize.brentq`; on `ValueError` (no sign change in bracket) skip that sample silently. Document this in the function. |
| `linkable` false-negative (declares two bodies un-linkable when they are) | low | **high — would prune real cyclers in M4** | Add tests at known-linkable V∞ values (E–M @ 5.5 km/s is Aldrin; E–V @ 5.0 km/s; M–V @ 4.0 km/s). Sample contours densely (`n_points=200` default). Use a generous intersection tolerance (`tol_au=0.01`, `tol_e=0.01`). |
| `linkable` false-positive (declares un-linkable bodies linkable) | medium | medium — wastes M4 compute but doesn't lose cyclers | Acceptable failure direction. The M4 enumerator will discover the lack of a real construction during inner search and discard the cell with a `failed` ledger entry. |
| **Single-impulse `flyby_dv` model is too coarse** | medium | low (M2) → medium (M5) | Document the surrogate explicitly; flag as a model-fidelity revisit in M5 if optimiser pathologies appear. The M2 gate only requires "zero on ballistic-feasible, positive otherwise" — both satisfied. |
| Coplanar restriction misunderstood downstream | medium | **high — silently using i=0 on inclined real ephemeris would underestimate ΔV** | Loud docstring on every Tisserand entry point. Hand-off note to M3 explicitly states "M2 Tisserand is coplanar." M6 will introduce 3D Tisserand under a new entry point, not by quietly upgrading these. |
| VEM beat tolerance too tight rejects the real (3, 4) tuple | low | medium — gate fails | The real ratio 3·E–M / 4·E–V ≈ 1.00097; `tol_frac = 0.02` accepts up to 2 % mismatch. Confirmed sufficient by hand calculation; if the gate fails the cause will be a wrong synodic, not a wrong tolerance. |
| `matplotlib` import in CI without viz extra | low | low | Guard the import inside `plot_tisserand`; never import at module top. Test skipped if absent. |
| numpy 2.x vs scipy version mismatch under `mypy --strict` | medium | low | Pin `scipy>=1.13` (numpy-2.x compatible) per M0 §7 forecast. Use `from numpy.typing import NDArray` everywhere; avoid bare `np.ndarray` in signatures. |

---

## 6. Dependency additions

Per M0 plan §7. Additions in **this milestone** only:

| Package | Group | Reason | Version pin |
|---------|-------|--------|-------------|
| `scipy` | runtime | `scipy.optimize.brentq` for contour root-finding and Tisserand intersection. Likely already present from M1 (Lambert universal-variable Newton iterations), but explicit declaration here in case M1 inlined the solver. | `>=1.13` |
| `matplotlib` | optional `viz` extra | `plot_tisserand` diagnostic helper. CI runs viz tests when the extra is installed; not required for the M2 gate. | `>=3.9` |

`pyproject.toml` edits:

```toml
[project.dependencies]
# numpy >=2.0 already declared in M0
# scipy >=1.13 likely added in M1; add here if not present
scipy = ">=1.13"

[project.optional-dependencies]
# dev = [...] from M0 unchanged
viz = ["matplotlib>=3.9"]
```

Run `uv sync --all-extras` after editing to refresh `uv.lock`; commit the lockfile change.

CI (`.github/workflows/ci.yml`) already runs `uv sync --frozen --all-extras` per M0, so the viz extra is picked up automatically.

---

## 7. Order of work

The todo.md mirrors this as checkboxes.

1. **Read predecessor docs.** Re-read M0 `plan.md` constants section; re-read M1 `plan.md`/`todo.md` to confirm `scipy` is already a dep and to inherit any conventions M1 established (e.g. vector type aliases).
2. **Update `pyproject.toml`** with `scipy` (if absent) and the `viz` extra. Run `uv sync --all-extras`; commit `uv.lock`.
3. **`core/flyby.py`:**
   1. Implement `max_bend`, `bend_angle`.
   2. Implement `is_ballistic_feasible`.
   3. Implement `flyby_dv` with the single-impulse decomposition.
   4. Implement `flyby_dv_for` planet-aware wrapper.
   5. Write `tests/test_flyby.py`; assert Mars bend = 24° gate passes locally.
4. **`search/__init__.py`** (empty file to register the subpackage).
5. **`search/resonance.py`:**
   1. Implement `synodic_period_days`, `synodic_period_years`, `k_synodic_periods_days`.
   2. Implement `multi_body_beat_days`, `beat_period_days`.
   3. Write `tests/test_resonance.py`; assert all three resonance gate anchors pass (E–M 2.135 yr, E–V 1.599 yr, VEM 6.406 yr).
6. **`search/tisserand.py`:**
   1. Implement `vinf_to_tisserand`, `tisserand_to_vinf`.
   2. Implement `vinf_contour` (the `e`-parameterised cubic-in-`u` root-find via `brentq`).
   3. Implement `linkable` (sign-change search on `a_a(e) − a_b(e)`).
   4. Implement `linkable_region` (sample V∞ grid).
   5. Implement `plot_tisserand` (lazy matplotlib import).
   6. Write `tests/test_tisserand.py`; assert contour smoke + `linkable` physics asserts pass.
7. **Run the full local quality gate:** `uv run pytest`, `uv run ruff check .`, `uv run ruff format --check .`, `uv run mypy src tests`.
8. **Commit:** `m2: flyby mechanics, Tisserand contours, synodic + multi-body beat (gate: E-M 2.135 yr, E-V 1.599 yr, VEM 6.406 yr, Mars bend 24°)`. Push; confirm CI green.
9. **Update `docs/overview.md`** §4 milestone table: M2 status → completed; M3 status → planned.
10. **Hand-off note** appended to `todo.md` under `## Hand-off to M3`.

---

## 8. Exit checklist

Before declaring M2 done:

- [ ] `uv run pytest` green locally on a fresh checkout (`uv sync --all-extras` from scratch).
- [ ] All four spec §9 gate anchors pass as **explicitly named** test functions (Mars 24° bend, E–M 2.135 yr, E–V 1.599 yr, VEM 6.406 yr).
- [ ] Earth & Venus bend-range sanity tests pass.
- [ ] `flyby_dv == 0` ↔ `is_ballistic_feasible` cross-check parametrised test passes.
- [ ] `linkable("E","M", 5.5)` returns True (Aldrin neighbourhood); `linkable("E","M", 0.5)` returns False; `linkable` symmetric.
- [ ] `uv run ruff check .` clean.
- [ ] `uv run ruff format --check .` clean.
- [ ] `uv run mypy src tests` clean — including new modules and tests (numpy typing for vector arguments).
- [ ] `viz` extra documented in README quick-start.
- [ ] CI green on push for the M2 commit.
- [ ] `docs/overview.md` updated: M2 → completed; M3 row marked `planned`.
- [ ] `## Hand-off to M3` section appended to `phases/m2-flyby-maps/todo.md` covering: the coplanar-only restriction (must be re-stated when M3 starts using `linkable`), any model-fidelity caveats discovered in `flyby_dv`, and the period bank that `k_synodic_periods_days("E","M", 2)` gives M3 for the Aldrin reproduction.

(Writing the M3 plan doc is the first task of M3, not an M2 exit criterion.)
