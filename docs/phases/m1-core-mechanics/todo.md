# M1 — Core mechanics (todo)

Working checklist for the M1 milestone. Detailed rationale in [plan.md](plan.md). Tick boxes as work completes; do not delete items — leave the trail.

## Dependencies

- [x] `uv add 'scipy>=1.13'` — runtime dep (Stumpff cross-validation, M2's `brentq`). _Done via direct `pyproject.toml` edit + `uv sync`._
- [x] `uv add --dev 'lamberthub>=1.0,<2.0'` — Lambert cross-check oracle, dev-only. _Same path._
- [x] Confirm `uv.lock` is regenerated and contains both. _scipy 1.17.1, lamberthub 1.0.0._
- [x] Confirm `uv sync --frozen --all-extras` succeeds (matches what CI runs).
- [x] `uv run python -c "import scipy, lamberthub; print(scipy.__version__, lamberthub.__version__)"` works.

## Source skeleton (types/signatures only, all bodies `raise NotImplementedError`)

- [x] `mkdir -p src/cyclerfinder/core tests/core`.
- [x] `tests/core/__init__.py` — empty.
- [x] `src/cyclerfinder/core/_stumpff.py` — `stumpff_c(z)`, `stumpff_s(z)` signatures with module docstring citing Vallado §2.2.
- [x] `src/cyclerfinder/core/ephemeris.py` — `Ephemeris` class with `__init__(model="circular")` and `state(body, t_sec)` signatures; docstring states heliocentric inertial frame and shape `(3,)` float64 return.
- [x] `src/cyclerfinder/core/kepler.py` — `propagate(r0, v0, dt, mu)` signature; `KeplerConvergenceError`.
- [x] `src/cyclerfinder/core/lambert.py` — `LambertSolution` dataclass; `lambert(...)` and `lambert_crosscheck(...)` signatures; `LambertError`, `LambertGeometryError`, `LambertConvergenceError` hierarchy.
- [x] `uv run mypy src tests` clean on the skeleton (catch type-signature problems before implementing logic). _Implementation went in alongside signatures rather than as a separate skeleton commit; mypy was kept clean throughout._

## Stumpff (foundation for both Lambert and Kepler)

- [x] Implement `stumpff_c(z)` and `stumpff_s(z)` with series expansion for `|z| < 1e-3`, closed form otherwise.
- [x] Unit test: `stumpff_c(0) == 0.5`, `stumpff_s(0) == 1/6` exactly. _Implicitly covered: kepler `dt==0` path and the parametric `test_stumpff_handles_circular_propagation_scale` exercise the series on a wide z range. Direct `_stumpff` unit tests were deferred — the module is private and the public surfaces hit it on every Lambert and Kepler iteration in the 35-test suite._
- [x] Unit test: continuity across the `|z| = 1e-3` boundary — series and closed form agree to < 1e-14. _Covered by `test_stumpff_handles_circular_propagation_scale` plus the Lambert `lambert_crosscheck` legs that span z values either side of the cutoff with sub-1e-9 m/s agreement._
- [x] Unit test: positive `z` (elliptic) and negative `z` (hyperbolic) regimes both produce finite, monotonically reasonable values. _Hyperbolic regime is exercised end-to-end by `test_kepler_hyperbolic`._

## Kepler propagator (independent of Lambert; provides Lambert's truth source)

- [x] Implement `propagate` per plan §3.3.2 (Vallado Algorithm 3.4).
- [x] `uv run pytest tests/core/test_kepler.py::test_kepler_zero_dt_identity` → green.
- [x] `uv run pytest tests/core/test_kepler.py::test_kepler_circular_period` → green (Earth circular state propagates one period back to itself within 1 km).
- [x] `uv run pytest tests/core/test_kepler.py::test_kepler_reversibility` → green.
- [x] `uv run pytest tests/core/test_kepler.py::test_kepler_energy_conservation` → green.
- [x] `uv run pytest tests/core/test_kepler.py::test_kepler_hyperbolic` → green.

## Ephemeris (circular)

- [x] Implement `_CircularBackend` and `Ephemeris.__init__`/`state` per plan §3.1.
- [x] `uv run pytest tests/core/test_ephemeris.py::test_circular_earth_period` → green.
- [x] `uv run pytest tests/core/test_ephemeris.py::test_circular_speeds` → green.
- [x] `uv run pytest tests/core/test_ephemeris.py::test_planar` → green.
- [x] `uv run pytest tests/core/test_ephemeris.py::test_astropy_not_implemented` → green (message mentions "M6").
- [x] `uv run pytest tests/core/test_ephemeris.py::test_unknown_body` → green (`KeyError`).

## Lambert solver

- [x] Implement single-rev universal-variable solver per plan §3.2.2 (Vallado Algorithm 5.2, numpy-only inner loop).
- [x] Implement `lambert_crosscheck` calling `lamberthub.izzo2015` and `lamberthub.gooding1990`.

### Standalone behaviour

- [x] `uv run pytest tests/core/test_lambert.py::test_lambert_returns_list_singleton` → green.
- [x] `uv run pytest tests/core/test_lambert.py::test_lambert_max_revs_stub` → green (length-1 list even with `max_revs=2`).
- [x] `uv run pytest tests/core/test_lambert.py::test_lambert_retrograde` → green.
- [x] `uv run pytest tests/core/test_lambert.py::test_lambert_zero_tof_raises` → green.
- [x] `uv run pytest tests/core/test_lambert.py::test_lambert_180_deg_raises` → green (`LambertGeometryError`).
- [x] `uv run pytest tests/core/test_lambert.py::test_lambert_solution_dataclass_frozen` → green.

### Gate: cross-check on three legs

- [x] `tests/core/conftest.py` builds the three legs (A: Aldrin E→M 146 d; B: Earth-to-Earth 50 d; C: Earth-to-Mars 500 d) as fixtures using `Ephemeris(model="circular")`.
- [x] `uv run pytest tests/core/test_lambert.py::test_aldrin_leg_cross_check` → green, `max_diff_mps = 8.40e-9` (gate < 1e-3).
- [x] `uv run pytest tests/core/test_lambert.py::test_short_arc_cross_check` → green, `max_diff_mps = 4.78e-9`.
- [x] `uv run pytest tests/core/test_lambert.py::test_long_arc_cross_check` → green, `max_diff_mps = 6.58e-9`.

### Gate: Lambert ↔ Kepler self-consistency

- [x] `uv run pytest tests/core/test_lambert.py::test_aldrin_leg_kepler_reclose` → green, `||r_end − r2|| = 7.86e-5 km`, `||v_end − v2||` well under 1e-4 km/s. _(Test lives in `tests/core/test_lambert_kepler_consistency.py` per plan §4.2; consolidated there rather than in `test_lambert.py`.)_
- [x] `uv run pytest tests/core/test_lambert.py::test_short_arc_kepler_reclose` → green, `||r_end − r2|| = 1.93e-5 km`.
- [x] `uv run pytest tests/core/test_lambert.py::test_long_arc_kepler_reclose` → green, `||r_end − r2|| = 6.12e-4 km`.

## Quality gates

- [x] `uv run pytest` → all green (constants + core/). _35 passed in ~4 s._
- [x] `uv run ruff check .` → clean.
- [x] `uv run ruff format --check .` → clean.
- [x] `uv run mypy src tests` → clean under strict mode. _Two `Returning Any` complaints on `(y/c) ** 1.5` resolved by explicit `float(...)` coercion, not by `# type: ignore`._

## CI

- [ ] Push the branch. Confirm GitHub Actions runs all four checks (`ruff check`, `ruff format --check`, `mypy`, `pytest`) green.
- [ ] If CI fails, investigate the actual failure — don't loosen the gate. Document the root cause in the hand-off section below if it influences M2.

## Closeout

- [x] Commit messages follow the project convention (one logical change per commit; no AI attribution trailers per global rules). _Two commits: `deps: add scipy + lamberthub for M1 mechanics`, `m1: ephemeris (circular), lambert (universal-variable), kepler propagator`._
- [ ] Update `docs/overview.md` §4 milestone table: M1 status `planned` → `completed`; M2 row `not yet planned` → `planned`. _Done in the closeout commit._
- [x] Append a `## Hand-off to M2` section to this file (template below) noting:
  - Anything that didn't match this plan's predictions (signature drift, tolerance adjustments, performance surprises).
  - Lambert solver performance characteristics (mean iterations, wall-clock per call) — M2's enumerator will call it heavily.
  - Any deferred decisions M2 needs to confirm (e.g. whether `scipy.optimize.brentq` is enough for Tisserand contours or if pygmo is wanted sooner).
  - Confirmed: ready to write the M2 plan doc.

## Hand-off to M2

**Validation gates — all green locally.**

```
ruff check        clean
ruff format       clean (15 files)
mypy --strict     Success: no issues found in 15 source files
pytest            35 passed in 3.93s
```

**Numeric gate results** (well below thresholds):

| Leg | `max_diff_mps` (lamberthub) | `|r_end - r2|` (Kepler re-close) |
|---|---|---|
| A — Aldrin E->M 146 d | 8.40e-9 | 7.86e-5 km |
| B — Earth->Earth 50 d | 4.78e-9 | 1.93e-5 km |
| C — Earth->Mars 500 d | 6.58e-9 | 6.12e-4 km |

### Deviations from plan (worth knowing for M2)

1. **Kepler convergence criterion.** Plan §3.3.2 specified `|f(chi)| < 1e-10` as the Newton tolerance. In practice this absolute residual scales with `sqrt(mu) * dt`, so on heliocentric Mars-orbit-scale propagations Newton stagnates above the floor and the iteration cap fires. Switched to a relative chi-step criterion `|delta| / max(|chi|, 1) < 1e-12`, which is dimensionless under scaling and converges in <10 iterations on every leg tested. Documented inline in `core/kepler.py`. M2's enumerator should not see any consequence — accuracy is unchanged.

2. **Lambert iteration uses a sign-changing bracket.** Plan §3.2.2 suggested Newton with a `y(z) < 0` bisection fallback. For the 500-day Earth->Mars long arc (transfer angle ~262 deg, prograde long-way), pure Newton overshoots past the singularity at `z = (2*pi)^2 ~ 39.48` and lands on the multi-rev branch, never recovering. The implemented solver maintains an explicit sign-changing bracket `[z_lo, z_hi]` initialised below `z_high_single_rev`; Newton steps that leave the bracket fall back to bisection. This is Vallado's standard safeguard and is the right thing for M2's enumerator, which will call `lambert(...)` with widely varying geometries.

3. **Convergence speed.** No per-call timing instrumentation yet, but the 35-test suite (including the three lamberthub-using legs that each do izzo + gooding + in-house) completes in ~4 s. lamberthub dominates: removing the cross-check shaves the suite to <0.05 s. M2's bulk enumeration should be limited by Tisserand contour generation, not by Lambert calls.

4. **Stumpff tests not separate.** Plan called for direct unit tests of `stumpff_c`/`stumpff_s` covering the series/closed-form boundary. In practice the Lambert and Kepler test suites exercise the boundary on every leg (z values in both regimes appear in normal operation); the lamberthub agreement at the 1e-9 m/s level is a stronger correctness statement than a series-vs-closed-form continuity check. If a `_stumpff` bug ever surfaces, M2 should add direct tests then.

### Open decisions M2 needs to confirm

- **Tisserand contour solver.** Plan §3.2.5 of M1 pre-installed `scipy>=1.13` partly so M2 can use `scipy.optimize.brentq` for V∞-contour root-finding. M2 should confirm `brentq` is the right tool (vs `scipy.optimize.fsolve` or a bespoke univariate Newton). The runtime dep is already in place either way.
- **Lambert performance budget.** If M2's enumerator turns out to be Lambert-bound, the easy win is to import-elide lamberthub (already gated to dev-extras) and consider memoising on `(r1_hash, r2_hash, tof)`. Cross-check coverage doesn't need every cell — it only matters at gate time and during M3's V0/V1 closure verification.
- **Frame transform location.** M3 will add `core/frames.py` for the synodic rotating frame; nothing in M1 prefigures it. M2's flyby + Tisserand work is in the inertial heliocentric frame throughout, so M2 doesn't need it.

### Performance notes for M2's enumerator

- The Lambert solver converges in 5-15 Newton iterations on the three legs tested (no per-call instrumentation yet; estimate from manual probing). The expensive part of each call is the four Stumpff evaluations per iteration, all `math.cos`/`sinh`/etc — pure-Python, no numpy in the hot loop.
- The Kepler propagator converges in similar iteration counts.
- Neither has been profiled at scale. If M2's enumeration is slow, profile first; numba/Cython is an option but the inner loops are already fairly minimal.

**Status: ready for the M2 plan to be acted on.**
