# M2 — Flyby + maps (todo)

Working checklist for the M2 milestone. Detailed rationale in [plan.md](plan.md). Tick boxes as work completes; do not delete items — leave the trail.

## Predecessor recap

- [ ] Re-read `phases/m0-scaffold/plan.md` §4 (constants — `PLANETS`, `SAFE_PERIHELION_KM`, `MU_SUN_KM3_S2`, `AU_KM`, `DAYS_PER_JULIAN_YEAR`). These are the only inputs M2 needs from M0.
- [ ] Re-read `phases/m1-core-mechanics/plan.md` and `todo.md`. Confirm:
  - whether `scipy` is already a runtime dep (Lambert universal-variable solver may have added it);
  - any vector type alias M1 established (e.g. `Vec3 = NDArray[np.float64]`) — reuse it in M2 for consistency, do not redeclare.
- [ ] Re-read `phases/m0-scaffold/todo.md` § Hand-off to M1 and `phases/m1-core-mechanics/todo.md` § Hand-off to M2 for any predecessor caveats.

## Dependency updates

- [ ] Edit `pyproject.toml`:
  - [ ] Add `scipy>=1.13` to `[project.dependencies]` if not already present from M1.
  - [ ] Add `[project.optional-dependencies] viz = ["matplotlib>=3.9"]`.
- [ ] `uv sync --all-extras` — refresh `uv.lock`.
- [ ] `uv run python -c "import scipy.optimize, matplotlib.pyplot"` smoke-imports without error.
- [ ] Update `README.md` quick-start: mention `--all-extras` enables the Tisserand plotting helper.

## `core/flyby.py`

- [ ] Create file with module docstring stating the spec §3 bend formula and §6 interface origin.
- [ ] Implement `max_bend(mu_planet: float, rp_min: float, vinf: float) -> float`:
  - [ ] Guard `vinf == 0` returning `np.pi`.
  - [ ] Closed-form `2 * arcsin(1 / (1 + rp_min * vinf**2 / mu_planet))`.
- [ ] Implement `bend_angle(vin_vec, vout_vec) -> float` — numerically robust via `np.clip` on the acos argument.
- [ ] Implement `is_ballistic_feasible(vin_vec, vout_vec, mu_planet, rp_min, speed_tol=1e-6) -> bool`.
- [ ] Implement `flyby_dv(vin_vec, vout_vec, mu_planet, rp_min) -> float`:
  - [ ] Single-impulse decomposition per plan §3.1.
  - [ ] Returns `0.0` exactly when `is_ballistic_feasible` is True.
- [ ] Implement `flyby_dv_for(code, vin_vec, vout_vec) -> float` reading `PLANETS` + `SAFE_PERIHELION_KM`.
- [ ] All public functions typed: scalars `float`, vectors `NDArray[np.float64]`. Module passes `mypy --strict`.

## `tests/test_flyby.py`

- [ ] `test_mars_max_bend_24deg_at_7kms` (spec §9 gate, tol ±1°).
- [ ] `test_earth_max_bend_in_range_at_7kms` (∈ [60, 63]°).
- [ ] `test_venus_max_bend_in_range_at_7kms` (∈ [60, 63]°).
- [ ] `test_max_bend_zero_at_infinite_vinf`.
- [ ] `test_max_bend_pi_at_zero_vinf`.
- [ ] `test_ballistic_feasible_zero_dv` (equal magnitudes, bend = 0.5·max_bend → `flyby_dv == 0.0`).
- [ ] `test_overbent_pair_positive_dv` (bend = 1.5·max_bend → `flyby_dv > 0`).
- [ ] `test_speed_mismatch_positive_dv` (5 vs 7 km/s, bend = 0 → ΔV ≥ ~2 km/s).
- [ ] `test_flyby_dv_for_matches_explicit` (the wrapper is a thin pass-through).
- [ ] `test_is_ballistic_feasible_consistency` — parametrised over ≥50 random `(vin, vout)` pairs: `flyby_dv == 0` iff `is_ballistic_feasible`.
- [ ] All tests pass: `uv run pytest tests/test_flyby.py`.

## `search/` subpackage

- [ ] `mkdir -p src/cyclerfinder/search`.
- [ ] `src/cyclerfinder/search/__init__.py` — empty (subpackage marker).

## `search/resonance.py`

- [ ] Create file with module docstring.
- [ ] Implement `synodic_period_days(body_a, body_b) -> float`:
  - [ ] Read `mean_motion_deg_day` from `PLANETS`; convert to period via `360 / n`.
  - [ ] `1 / |1/T_a - 1/T_b|`.
  - [ ] Raise `ValueError` on `body_a == body_b`.
- [ ] Implement `synodic_period_years(body_a, body_b) -> float` — divides by `DAYS_PER_JULIAN_YEAR`.
- [ ] Implement `k_synodic_periods_days(body_a, body_b, k_max) -> list[float]`.
- [ ] Implement `multi_body_beat_days(bodies, k_max=6, tol_frac=0.02) -> list[tuple[int, ...]]`:
  - [ ] Reference-body rule per plan §3.3: middle element for 3-body sets.
  - [ ] Exhaustive integer search over `k_i ∈ [1, k_max]`.
  - [ ] Ranked by fractional mismatch ascending; tuples within `tol_frac` returned.
- [ ] Implement `beat_period_days(bodies, k_tuple) -> float` — mean of `k_i * T_syn(body_i, ref)`.
- [ ] Module passes `mypy --strict`.

## `tests/test_resonance.py`

- [ ] `test_em_synodic_2135yr` — `synodic_period_years("E","M")` ≈ 2.135 yr ±0.001 (**spec §9 gate**).
- [ ] `test_ev_synodic_1599yr` — `synodic_period_years("E","V")` ≈ 1.599 yr ±0.001 (**spec §9 gate**).
- [ ] `test_synodic_symmetric` — order-independent.
- [ ] `test_synodic_self_raises` — `synodic_period_days("E","E")` raises `ValueError`.
- [ ] `test_vem_beat_yields_3_4` — `(3, 4)` (or equivalent ordering) appears in `multi_body_beat_days(["V","E","M"], k_max=6)` (**spec §9 gate**).
- [ ] `test_vem_beat_period_6_406yr` — `beat_period_days(["V","E","M"], (3, 4)) / 365.25` ≈ 6.406 yr ±0.01 (**spec §9 gate**).
- [ ] `test_k_synodic_monotone` — `k_synodic_periods_days("E","M", 5)` strictly increasing.
- [ ] All tests pass: `uv run pytest tests/test_resonance.py`.

## `search/tisserand.py`

- [ ] Create file with module docstring **explicitly stating the coplanar (i=0) restriction** and the spec §13.3 role of `linkable`.
- [ ] Implement `vinf_to_tisserand(body, vinf_kms) -> float`:
  - [ ] `T_p = 3 - V∞² * a_p / μ_sun` (convert `sma_au` to km via `AU_KM`).
- [ ] Implement `tisserand_to_vinf(body, T_p) -> float` — guard returns `0.0` when `T_p >= 3`.
- [ ] Implement `vinf_contour(body, vinf_kms, a_range_au=(0.3, 5.0), n_points=200) -> tuple[NDArray, NDArray]`:
  - [ ] Sample `e ∈ [0, 1 - 1e-4)` at `n_points` points.
  - [ ] For each `e`, solve the cubic-in-`u = sqrt(a/a_p)`: `1 + 2 u³ sqrt(1 - e²) = T_target u²` via `scipy.optimize.brentq`.
  - [ ] Drop samples where no real `a` in `a_range_au` solves the equation (do **not** raise).
  - [ ] Empty arrays when no contour exists at this V∞.
- [ ] Implement `linkable(body_a, body_b, vinf_kms, tol_au=0.01, tol_e=0.01) -> bool`:
  - [ ] Build both contours.
  - [ ] On shared e-support, look for sign change in `a_a(e) - a_b(e)`; refine via `brentq`.
  - [ ] Return True iff intersection within `(tol_au, tol_e)` exists.
  - [ ] Never raise — wrap `brentq` `ValueError`s and return False on no bracket.
- [ ] Implement `linkable_region(body_a, body_b, vinf_cap_kms, n_vinf=50) -> list[float]`:
  - [ ] Sample V∞ on `(0, vinf_cap_kms]`.
  - [ ] Return the subset where `linkable` is True.
- [ ] Implement `plot_tisserand(bodies, vinf_levels_kms, ax=None) -> Axes`:
  - [ ] Lazy `import matplotlib.pyplot as plt` inside the function.
  - [ ] Raise `ImportError("install the `viz` extra: uv sync --extra viz")` if matplotlib absent.
  - [ ] Overlay contours; label each curve `(body, V∞)`.
- [ ] Module passes `mypy --strict`. Matplotlib import is conditional so mypy must not require the stub on the default install.

## `tests/test_tisserand.py`

- [ ] `test_vinf_to_tisserand_inverse` — round-trip exact at `V∞ = 5.0` for Earth.
- [ ] `test_contour_returns_empty_below_threshold` — `vinf_contour("E", 0.01)` returns two empty arrays.
- [ ] `test_contour_non_empty_at_5kms` — `vinf_contour("E", 5.0)` returns ≥ 50 points.
- [ ] `test_contour_a_e0_at_zero_vinf` — at e=0 and V∞=0, contour passes through `a = a_p`.
- [ ] `test_linkable_em_at_5_5_kms_true` (Aldrin neighbourhood — physics).
- [ ] `test_linkable_em_at_0_5_kms_false` (insufficient energy).
- [ ] `test_linkable_symmetric` — parametrised over several `(pair, V∞)` triples.
- [ ] `test_linkable_region_em_non_empty` — at `vinf_cap = 12.0`.
- [ ] `test_linkable_never_raises` — parametrised over `[0.0, 1e-6, 1e6]` V∞.
- [ ] `test_plot_runs_when_matplotlib_present` — `@pytest.mark.skipif` guard on import; asserts return is an `Axes`.
- [ ] All tests pass: `uv run pytest tests/test_tisserand.py`.

## Local green

- [ ] `uv run pytest` → green (all M0/M1/M2 tests).
- [ ] `uv run ruff check .` → clean.
- [ ] `uv run ruff format --check .` → clean (run `uv run ruff format .` if it isn't).
- [ ] `uv run mypy src tests` → clean (resolve any strict-mode complaints in new modules — typical: numpy generic args, missing return types on test functions).

## CI

- [ ] Commit: `m2: flyby mechanics, Tisserand contours, synodic + multi-body beat`.
  - Body lists the four gate checks: E-M 2.135 yr, E-V 1.599 yr, VEM 6.406 yr beat, Mars bend 24°.
- [ ] Push. Confirm GitHub Actions runs green.

## Closeout

- [ ] Update `docs/overview.md` §4 milestone table: M2 status `planned` → `completed`. M3 row `not yet planned` → `planned`.
- [ ] Append a `## Hand-off to M3` section to this todo.md noting:
  - **Coplanar restriction must be re-stated** when M3 starts using `linkable` / contours — M3's Aldrin reproduction is coplanar so this is fine, but it must be loud so it isn't carried silently into M6.
  - The period bank `k_synodic_periods_days("E","M", 2)` (≈ 2.135 yr and 4.27 yr) is what M3 uses for the Aldrin (k=1) and 2-synodic (k=2) reproductions.
  - Any model-fidelity caveats discovered during M2 in the single-impulse `flyby_dv` decomposition that M3/M5 should revisit.
  - Anything that didn't go as the plan predicted (e.g. `brentq` bracketing surprises in contour sampling, mypy/scipy stub friction).
  - Confirmed: ready to write the M3 plan doc.

## Hand-off to M3

_To be filled in when M2 completes._
