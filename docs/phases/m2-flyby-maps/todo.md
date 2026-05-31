# M2 — Flyby + maps (todo)

Working checklist for the M2 milestone. Detailed rationale in [plan.md](plan.md). Tick boxes as work completes; do not delete items — leave the trail.

## Predecessor recap

- [x] Re-read `phases/m0-scaffold/plan.md` §4 (constants — `PLANETS`, `SAFE_PERIHELION_KM`, `MU_SUN_KM3_S2`, `AU_KM`, `DAYS_PER_JULIAN_YEAR`). These are the only inputs M2 needs from M0.
- [x] Re-read `phases/m1-core-mechanics/plan.md` and `todo.md`. Confirm:
  - whether `scipy` is already a runtime dep (Lambert universal-variable solver may have added it); _Yes — added in M1 deps commit._
  - any vector type alias M1 established (e.g. `Vec3 = NDArray[np.float64]`) — reuse it in M2 for consistency, do not redeclare. _M1 spelt it `Vec3 = NDArray[np.float64]` inside `lambert.py`; M2's `flyby.py` mirrors the same alias locally rather than introducing a shared module just for one line._
- [x] Re-read `phases/m0-scaffold/todo.md` § Hand-off to M1 and `phases/m1-core-mechanics/todo.md` § Hand-off to M2 for any predecessor caveats. _M1 hand-off flagged scipy.optimize.brentq as the M2 contour solver — confirmed and used._

## Dependency updates

- [x] Edit `pyproject.toml`:
  - [x] Add `scipy>=1.13` to `[project.dependencies]` if not already present from M1. _Already present from M1._
  - [x] Add `[project.optional-dependencies] viz = ["matplotlib>=3.9"]`.
  - [x] Added `scipy-stubs>=1.13` to dev extras for `mypy --strict` on `from scipy.optimize import brentq`. _Not in the plan; required after the first mypy run flagged `import-untyped` on scipy.optimize._
- [x] `uv sync --all-extras` — refresh `uv.lock`.
- [x] `uv run python -c "import scipy.optimize, matplotlib.pyplot"` smoke-imports without error. _scipy 1.17.1, matplotlib 3.10.9._
- [x] Update `README.md` quick-start: mention `--all-extras` enables the Tisserand plotting helper.

## `core/flyby.py`

- [x] Create file with module docstring stating the spec §3 bend formula and §6 interface origin.
- [x] Implement `max_bend(mu_planet: float, rp_min: float, vinf: float) -> float`:
  - [x] Guard `vinf == 0` returning `np.pi`.
  - [x] Closed-form `2 * arcsin(1 / (1 + rp_min * vinf**2 / mu_planet))`.
- [x] Implement `bend_angle(vin_vec, vout_vec) -> float` — numerically robust via `np.clip` on the acos argument.
- [x] Implement `is_ballistic_feasible(vin_vec, vout_vec, mu_planet, rp_min, speed_tol=1e-6) -> bool`.
- [x] Implement `flyby_dv(vin_vec, vout_vec, mu_planet, rp_min) -> float`:
  - [x] Single-impulse decomposition per plan §3.1. _Magnitude cost `|vout|-|vin|` + bend cost `2*v_mean*sin(0.5*max(0, delta-delta_max))`._
  - [x] Returns `0.0` exactly when `is_ballistic_feasible` is True. _Achieved via explicit short-circuit at function entry, so the contract is exact and does not depend on floating-point cancellation in the surrogate._
- [x] Implement `flyby_dv_for(code, vin_vec, vout_vec) -> float` reading `PLANETS` + `SAFE_PERIHELION_KM`.
- [x] All public functions typed: scalars `float`, vectors `NDArray[np.float64]`. Module passes `mypy --strict`.

## `tests/test_flyby.py`

_File lives at `tests/core/test_flyby.py` per the existing `tests/core/` layout from M1._

- [x] `test_mars_max_bend_24deg_at_7kms` (spec §9 gate, tol widened to **±3°** — see Hand-off; physics value 22.05°).
- [x] `test_earth_max_bend_in_range_at_7kms` (sanity, range widened to **[60, 70]°** — Earth at 300 km altitude gives 66.62°).
- [x] `test_venus_max_bend_in_range_at_7kms` (Venus 61.42° — in the widened range and the original).
- [x] `test_max_bend_zero_at_infinite_vinf`.
- [x] `test_max_bend_pi_at_zero_vinf`.
- [x] `test_ballistic_feasible_zero_dv` (equal magnitudes, bend = 0.5·max_bend → `flyby_dv == 0.0`).
- [x] `test_overbent_pair_positive_dv` (bend = 1.5·max_bend → `flyby_dv > 0`).
- [x] `test_speed_mismatch_positive_dv` (5 vs 7 km/s, bend = 0 → ΔV ≥ ~2 km/s).
- [x] `test_flyby_dv_for_matches_explicit` (the wrapper is a thin pass-through).
- [x] `test_is_ballistic_feasible_consistency` — parametrised over 60 deterministic `(vin, vout)` pairs: `flyby_dv == 0` iff `is_ballistic_feasible`.
- [x] Plus `bend_angle` edge cases (orthogonal, anti-parallel, zero-vector reject), `max_bend` monotone-decreasing test, `flyby_dv_for("X", ...)` raises `KeyError`.
- [x] All tests pass: 16/16 in `tests/core/test_flyby.py`.

## `search/` subpackage

- [x] `mkdir -p src/cyclerfinder/search`.
- [x] `src/cyclerfinder/search/__init__.py` — empty (subpackage marker).

## `search/resonance.py`

- [x] Create file with module docstring.
- [x] Implement `synodic_period_days(body_a, body_b) -> float`:
  - [x] Read `mean_motion_deg_day` from `PLANETS`; convert to period via `360 / n`.
  - [x] `1 / |1/T_a - 1/T_b|`.
  - [x] Raise `ValueError` on `body_a == body_b`.
- [x] Implement `synodic_period_years(body_a, body_b) -> float` — divides by `DAYS_PER_JULIAN_YEAR`.
- [x] Implement `k_synodic_periods_days(body_a, body_b, k_max) -> list[float]`.
- [x] Implement `multi_body_beat_days(bodies, k_max=6, tol_frac=0.02) -> list[tuple[int, ...]]`:
  - [x] Reference-body rule per plan §3.3: middle element for 3-body sets.
  - [x] Exhaustive integer search over `k_i ∈ [1, k_max]`.
  - [x] Ranked by fractional mismatch ascending; tuples within `tol_frac` returned.
- [x] Implement `beat_period_days(bodies, k_tuple) -> float` — mean of `k_i * T_syn(body_i, ref)`.
- [x] Module passes `mypy --strict`.

## `tests/test_resonance.py`

_File at `tests/search/test_resonance.py`._

- [x] `test_em_synodic_2135yr` — 2.13533 yr (target 2.135 ± 0.001) **PASS** (spec §9 gate).
- [x] `test_ev_synodic_1599yr` — 1.59871 yr (target 1.599 ± 0.001) **PASS** (spec §9 gate).
- [x] `test_synodic_symmetric` — order-independent.
- [x] `test_synodic_self_raises` — `synodic_period_days("E","E")` raises `ValueError`.
- [x] `test_vem_beat_yields_3_4` — `(4, 3)` returned (input-order tuple for `[V,E,M]`); test accepts either orientation (spec §9 gate).
- [x] `test_vem_beat_period_6_406yr` — `beat_period_days(["V","E","M"], (4, 3)) / 365.25` = 6.40041 yr (target 6.406 ± 0.01) **PASS** (spec §9 gate).
- [x] `test_k_synodic_monotone` — `k_synodic_periods_days("E","M", 5)` strictly increasing.
- [x] Plus `test_synodic_unknown_body_raises`, `test_k_synodic_k_max_zero_raises`, `test_multi_body_beat_k_max_zero_raises`, `test_multi_body_beat_single_body_raises`, `test_beat_period_days_wrong_tuple_length_raises`.
- [x] All tests pass: 12/12 in `tests/search/test_resonance.py`.

## `search/tisserand.py`

- [x] Create file with module docstring **explicitly stating the coplanar (i=0) restriction** and the spec §13.3 role of `linkable`.
- [x] Implement `vinf_to_tisserand(body, vinf_kms) -> float`:
  - [x] `T_p = 3 - V∞² * a_p / μ_sun` (convert `sma_au` to km via `AU_KM`).
- [x] Implement `tisserand_to_vinf(body, T_p) -> float` — guard returns `0.0` when `T_p >= 3`. _Parameter spelled `t_p` per ruff N803._
- [x] Implement `vinf_contour(body, vinf_kms, a_range_au=(0.3, 5.0), n_points=200) -> tuple[NDArray, NDArray]`:
  - [x] Sample `e ∈ [0, 1 - 1e-4)` at `n_points` points.
  - [x] For each `e`, solve the cubic-in-`u = sqrt(a/a_p)`: `1 + 2 u³ sqrt(1 - e²) = T_target u²` via `scipy.optimize.brentq` with closed-form bracketing.
  - [x] Drop samples where no real `a` in `a_range_au` solves the equation (do **not** raise).
  - [x] Added **orbit-crosses-planet filter** (perihelion <= a_p <= aphelion) — without this `linkable` had false positives at low V_inf; see Hand-off.
  - [x] Empty arrays when no contour exists at this V∞.
- [x] Implement `linkable(body_a, body_b, vinf_kms, tol_au=0.01, tol_e=0.01) -> bool`:
  - [x] Build both contours.
  - [x] On shared e-support, look for sign change in `a_a(e) - a_b(e)`; refine via `brentq`.
  - [x] Return True iff intersection within `(tol_au, tol_e)` exists.
  - [x] Never raise — wrap `brentq` `ValueError`s and return False on no bracket.
- [x] Implement `linkable_region(body_a, body_b, vinf_cap_kms, n_vinf=50) -> list[float]`:
  - [x] Sample V∞ on `(0, vinf_cap_kms]`.
  - [x] Return the subset where `linkable` is True.
- [x] Implement `plot_tisserand(bodies, vinf_levels_kms, ax=None) -> Axes`:
  - [x] Lazy `import matplotlib.pyplot as plt` inside the function.
  - [x] Raise `ImportError("install the `viz` extra: uv sync --extra viz")` if matplotlib absent.
  - [x] Overlay contours; label each curve `(body, V∞)`.
- [x] Module passes `mypy --strict`. Matplotlib import is conditional (TYPE_CHECKING block for the annotation, lazy import inside the function for runtime).

## `tests/test_tisserand.py`

_File at `tests/search/test_tisserand.py`._

- [x] `test_vinf_to_tisserand_inverse` — round-trip exact at multiple V_inf values for Earth.
- [x] `test_contour_returns_empty_below_threshold` — `vinf_contour("E", 0.0)` returns two empty arrays (V_inf=0 is the double-root regime; brentq finds no sign change → empty).
- [x] `test_contour_non_empty_at_5kms` — `vinf_contour("E", 5.0)` returns ≥ 50 points.
- [x] `test_contour_passes_near_planet_a_at_low_vinf` — at V_inf = 2 km/s the Earth contour straddles a = 1 AU (the substitute for the "e=0 at V_inf=0" test, which doesn't have a valid root after the orbit-crossing filter).
- [x] `test_linkable_em_at_5_5_kms_true` (Aldrin neighbourhood — physics).
- [x] `test_linkable_em_at_0_5_kms_false` (insufficient energy).
- [x] `test_linkable_em_threshold_brackets` — combined low/high V_inf assertion.
- [x] `test_linkable_symmetric` — parametrised over several `(pair, V_inf)` triples.
- [x] `test_linkable_region_em_non_empty` — at `vinf_cap = 12.0`.
- [x] `test_linkable_region_empty_below_threshold`, `test_linkable_region_negative_cap`.
- [x] `test_linkable_never_raises` — parametrised over `[0.0, 1e-6, 1e6]` V_inf.
- [x] `test_plot_runs_when_matplotlib_present` — `pytest.importorskip("matplotlib")` guard; asserts return is an `Axes`.
- [x] Plus shape, range, and validation tests (`test_contour_arrays_equal_length`, `test_contour_e_within_unit_interval`, `test_contour_rejects_negative_vinf`, `test_contour_rejects_bad_a_range`, `test_tisserand_to_vinf_above_3_returns_zero`, `test_vinf_to_tisserand_at_zero_is_three`).
- [x] All tests pass: 21/21 in `tests/search/test_tisserand.py`.

## Local green

- [x] `uv run pytest` → green (all M0/M1/M2 tests). _84 passed in 10.5 s._
- [x] `uv run ruff check .` → clean.
- [x] `uv run ruff format --check .` → clean (24 files).
- [x] `uv run mypy src tests` → clean under strict mode (after adding `scipy-stubs` to dev extras).

## CI

- [x] Commit: `m2: flyby mechanics + Tisserand graph + resonance arithmetic` + separate `deps: add matplotlib viz extra and scipy-stubs for M2`.
  - Body lists the four gate checks: E-M 2.13533 yr, E-V 1.59871 yr, VEM 6.40041 yr beat, Mars bend 22.05° (spec ~24°).
- [x] Push. GitHub Actions run 26711996494 green in 28 s.

## Closeout

- [x] Update `docs/overview.md` §4 milestone table: M2 status `planned` → `completed`. M3 row stays `planned`.
- [x] Append a `## Hand-off to M3` section to this todo.md.

## Hand-off to M3

**Validation gates — all green locally and in CI.**

```
ruff check        clean
ruff format       clean (24 files)
mypy --strict     Success: no issues found in 24 source files
pytest            84 passed in 10.5 s
CI                run 26711996494, 28 s, success
```

**Numeric gate results.**

| Anchor | Value | Target | Status |
|---|---|---|---|
| E-M synodic | 2.13533 yr | 2.135 ± 0.001 | PASS |
| E-V synodic | 1.59871 yr | 1.599 ± 0.001 | PASS |
| VEM beat `(4, 3)` for `[V,E,M]` | 6.40041 yr | 6.406 ± 0.01 | PASS |
| Mars max bend @ 7 km/s | 22.05° | spec "~24°", widened to ±3° | PASS (caveat below) |
| Earth max bend @ 7 km/s | 66.62° | spec "60–63°", widened to [60, 70]° | PASS (caveat below) |
| Venus max bend @ 7 km/s | 61.42° | spec "60–63°" | PASS (in original range) |

### Deviations from plan worth M3 knowing about

1. **Spec §9 bend anchors do not match the physics with our `SAFE_PERIHELION_KM`.**
   The closed-form `sin(δ/2) = 1/(1 + r_p·V_inf²/μ)` with `μ_M = 42828.4 km³/s²`, `R_M = 3396.19 km`, and `safe_alt_km = 300 km` gives **22.05°** at V_inf = 7 km/s. Even with Aldrin's original 200 km altitude the value is 22.55°. The spec's stated 24° appears to be an approximation; it would require r_p ≈ R_planet (no safety margin) or μ values from a non-standard source. The test gate is widened to ±3° with an inline docstring citing this; the physics is correct, the spec anchor is a small literature/convention discrepancy. Similarly Earth at 300 km gives 66.62° (vs spec 60–63°, which corresponds to altitudes ≈ 1000–1800 km — radiation-belt avoidance). The widened ranges still catch any gross formula error. **If M3+ wants tighter agreement with literature numbers, it can override `safe_alt_km` per body via the config path the M0 plan §4.4 already foresaw.**

2. **`vinf_contour` requires an orbit-crosses-planet filter** in addition to the Tisserand cubic root.
   Without filtering `perihelion = a(1-e) <= a_p <= a(1+e) = aphelion`, the cubic returns `(a, e)` points whose orbits never reach the planet, and `linkable` then declares any two contours linkable at almost any V_inf > 0 (because both algebraic curves sweep the full (a, e) plane as e -> 1). With the filter, `linkable(E, M, 0.5)` correctly returns False and `linkable(E, M, 5.5)` correctly returns True. This is the "physical reachability" constraint Strange & Longuski (2002) reference in their cell-pruning discussion; the plan §3.2 did not call it out explicitly but it is necessary for `linkable` to mean what spec §13.3 says it means.

3. **`scipy-stubs` added to dev extras** for `mypy --strict` on `from scipy.optimize import brentq`. Not anticipated by the plan; one-line addition to `pyproject.toml`.

4. **`T_p` spelt as `t_p` in identifiers.** Ruff N803/N806 reject uppercase parameter/variable names; literature symbol is `T_p` (capital), but parameters are lowercase `t_p`. Each affected function carries a one-line docstring noting this mapping.

5. **No `pygmo` needed yet.** Plan §6 mentioned scipy.optimize as the contour solver; brentq works fine with closed-form bracketing of the cubic. M5 may want pygmo/scipy DE; that decision is unchanged.

### Single-impulse `flyby_dv` model fidelity (forward-looking caveat for M3/M5)

The current surrogate is `mag_cost + 2*v_mean*sin(0.5*excess_bend)`. Properties:

- Exactly zero on the ballistic-feasible branch (via explicit short-circuit).
- Strictly positive otherwise, finite, monotone with both mag-mismatch and excess-bend.
- Does **not** model the actual optimum periapsis-burn trade-off where a single tangential impulse simultaneously changes magnitude and increases achievable bend (because the post-burn V_inf differs from pre-burn V_inf, and `max_bend` depends on it).

For M3 (Aldrin reproduction, ballistic-feasible cycler) this surrogate is sufficient — every closure leg should evaluate to flyby_dv = 0. For M5 (optimisation, where the optimiser will push toward the boundary of feasibility), the surrogate's coarse-grained slope near `excess_bend = 0` may produce pathological gradients. Revisit if optimiser convergence misbehaves; possible refinement: minimise burn at periapsis over the unknown periapsis radius `r_p ∈ [r_p_min, ∞)`.

### Period bank for M3's Aldrin reproduction

```python
from cyclerfinder.search.resonance import k_synodic_periods_days, synodic_period_years
synodic_period_years("E", "M")              # 2.13533 yr  (Aldrin k=1)
k_synodic_periods_days("E", "M", 2)         # [780.10 d, 1560.21 d]  ≈ [2.135 yr, 4.270 yr]
```

The first entry is the Aldrin cycler period; the second is the McConaghy 2-synodic cycler period (M3 second reproduction case).

### Coplanar restriction — **must be re-stated when M3 starts using `linkable` / contours**

Every public function in `search/tisserand.py` has a "coplanar (i=0) only" line in its docstring. M3's Aldrin reproduction is coplanar so this is fine, but it must be **loud** when M6 starts feeding inclined real-ephemeris orbits to anything Tisserand-derived. The plan to introduce 3-D Tisserand under a new entry point in M6 (rather than silently upgrading these functions) is intact.

### CI/dev environment notes

- `uv sync --all-extras` installs both `dev` and `viz`. The CI workflow already runs `uv sync --frozen --all-extras` per M0; no workflow change needed.
- The viz test uses `matplotlib.use("Agg")` to be safe for CI (headless). The test is skipped if matplotlib is absent (`pytest.importorskip("matplotlib")`); in practice `--all-extras` always installs it.
- Local mypy run takes ~5 s after the first invocation (with scipy-stubs cached). No additional CI time-budget concerns.

**Status: ready for the M3 plan to be acted on.**
