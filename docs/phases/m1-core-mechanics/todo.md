# M1 — Core mechanics (todo)

Working checklist for the M1 milestone. Detailed rationale in [plan.md](plan.md). Tick boxes as work completes; do not delete items — leave the trail.

## Dependencies

- [ ] `uv add 'scipy>=1.13'` — runtime dep (Stumpff cross-validation, M2's `brentq`).
- [ ] `uv add --dev 'lamberthub>=1.0,<2.0'` — Lambert cross-check oracle, dev-only.
- [ ] Confirm `uv.lock` is regenerated and contains both.
- [ ] Confirm `uv sync --frozen --all-extras` succeeds (matches what CI runs).
- [ ] `uv run python -c "import scipy, lamberthub; print(scipy.__version__, lamberthub.__version__)"` works.

## Source skeleton (types/signatures only, all bodies `raise NotImplementedError`)

- [ ] `mkdir -p src/cyclerfinder/core tests/core`.
- [ ] `tests/core/__init__.py` — empty.
- [ ] `src/cyclerfinder/core/_stumpff.py` — `stumpff_c(z)`, `stumpff_s(z)` signatures with module docstring citing Vallado §2.2.
- [ ] `src/cyclerfinder/core/ephemeris.py` — `Ephemeris` class with `__init__(model="circular")` and `state(body, t_sec)` signatures; docstring states heliocentric inertial frame and shape `(3,)` float64 return.
- [ ] `src/cyclerfinder/core/kepler.py` — `propagate(r0, v0, dt, mu)` signature; `KeplerConvergenceError`.
- [ ] `src/cyclerfinder/core/lambert.py` — `LambertSolution` dataclass; `lambert(...)` and `lambert_crosscheck(...)` signatures; `LambertError`, `LambertGeometryError`, `LambertConvergenceError` hierarchy.
- [ ] `uv run mypy src tests` clean on the skeleton (catch type-signature problems before implementing logic).

## Stumpff (foundation for both Lambert and Kepler)

- [ ] Implement `stumpff_c(z)` and `stumpff_s(z)` with series expansion for `|z| < 1e-3`, closed form otherwise.
- [ ] Unit test: `stumpff_c(0) == 0.5`, `stumpff_s(0) == 1/6` exactly.
- [ ] Unit test: continuity across the `|z| = 1e-3` boundary — series and closed form agree to < 1e-14.
- [ ] Unit test: positive `z` (elliptic) and negative `z` (hyperbolic) regimes both produce finite, monotonically reasonable values.

## Kepler propagator (independent of Lambert; provides Lambert's truth source)

- [ ] Implement `propagate` per plan §3.3.2 (Vallado Algorithm 3.4).
- [ ] `uv run pytest tests/core/test_kepler.py::test_kepler_zero_dt_identity` → green.
- [ ] `uv run pytest tests/core/test_kepler.py::test_kepler_circular_period` → green (Earth circular state propagates one period back to itself within 1 km).
- [ ] `uv run pytest tests/core/test_kepler.py::test_kepler_reversibility` → green.
- [ ] `uv run pytest tests/core/test_kepler.py::test_kepler_energy_conservation` → green.
- [ ] `uv run pytest tests/core/test_kepler.py::test_kepler_hyperbolic` → green.

## Ephemeris (circular)

- [ ] Implement `_CircularBackend` and `Ephemeris.__init__`/`state` per plan §3.1.
- [ ] `uv run pytest tests/core/test_ephemeris.py::test_circular_earth_period` → green.
- [ ] `uv run pytest tests/core/test_ephemeris.py::test_circular_speeds` → green.
- [ ] `uv run pytest tests/core/test_ephemeris.py::test_planar` → green.
- [ ] `uv run pytest tests/core/test_ephemeris.py::test_astropy_not_implemented` → green (message mentions "M6").
- [ ] `uv run pytest tests/core/test_ephemeris.py::test_unknown_body` → green (`KeyError`).

## Lambert solver

- [ ] Implement single-rev universal-variable solver per plan §3.2.2 (Vallado Algorithm 5.2, numpy-only inner loop).
- [ ] Implement `lambert_crosscheck` calling `lamberthub.izzo2015` and `lamberthub.gooding1990`.

### Standalone behaviour

- [ ] `uv run pytest tests/core/test_lambert.py::test_lambert_returns_list_singleton` → green.
- [ ] `uv run pytest tests/core/test_lambert.py::test_lambert_max_revs_stub` → green (length-1 list even with `max_revs=2`).
- [ ] `uv run pytest tests/core/test_lambert.py::test_lambert_retrograde` → green.
- [ ] `uv run pytest tests/core/test_lambert.py::test_lambert_zero_tof_raises` → green.
- [ ] `uv run pytest tests/core/test_lambert.py::test_lambert_180_deg_raises` → green (`LambertGeometryError`).
- [ ] `uv run pytest tests/core/test_lambert.py::test_lambert_solution_dataclass_frozen` → green.

### Gate: cross-check on three legs

- [ ] `tests/core/conftest.py` builds the three legs (A: Aldrin E→M 146 d; B: Earth-to-Earth 50 d; C: Earth-to-Mars 500 d) as fixtures using `Ephemeris(model="circular")`.
- [ ] `uv run pytest tests/core/test_lambert.py::test_aldrin_leg_cross_check` → green, `max_diff_mps < 1e-3`.
- [ ] `uv run pytest tests/core/test_lambert.py::test_short_arc_cross_check` → green, `max_diff_mps < 1e-3`.
- [ ] `uv run pytest tests/core/test_lambert.py::test_long_arc_cross_check` → green, `max_diff_mps < 1e-3`.

### Gate: Lambert ↔ Kepler self-consistency

- [ ] `uv run pytest tests/core/test_lambert.py::test_aldrin_leg_kepler_reclose` → green, `||r_end − r2|| < 1.0` km, `||v_end − v2|| < 1e-4` km/s.
- [ ] `uv run pytest tests/core/test_lambert.py::test_short_arc_kepler_reclose` → green.
- [ ] `uv run pytest tests/core/test_lambert.py::test_long_arc_kepler_reclose` → green.

## Quality gates

- [ ] `uv run pytest` → all green (constants + core/).
- [ ] `uv run ruff check .` → clean.
- [ ] `uv run ruff format --check .` → clean (run `uv run ruff format .` if not).
- [ ] `uv run mypy src tests` → clean under strict mode. (Document any targeted `# type: ignore` in the file with a justification comment.)

## CI

- [ ] Push the branch. Confirm GitHub Actions runs all four checks (`ruff check`, `ruff format --check`, `mypy`, `pytest`) green.
- [ ] If CI fails, investigate the actual failure — don't loosen the gate. Document the root cause in the hand-off section below if it influences M2.

## Closeout

- [ ] Commit messages follow the project convention (one logical change per commit; no AI attribution trailers per global rules):
  - `m1: add scipy + lamberthub deps`
  - `m1: stumpff helpers` (if split out)
  - `m1: kepler universal-variable propagator`
  - `m1: circular ephemeris backend`
  - `m1: lambert universal-variable solver + cross-check`
  - `m1: gate tests (lamberthub cross-check + kepler re-closure)`
- [ ] Update `docs/overview.md` §4 milestone table: M1 status `planned` → `completed`; M2 row `not yet planned` → `planned`.
- [ ] Append a `## Hand-off to M2` section to this file (template below) noting:
  - Anything that didn't match this plan's predictions (signature drift, tolerance adjustments, performance surprises).
  - Lambert solver performance characteristics (mean iterations, wall-clock per call) — M2's enumerator will call it heavily.
  - Any deferred decisions M2 needs to confirm (e.g. whether `scipy.optimize.brentq` is enough for Tisserand contours or if pygmo is wanted sooner).
  - Confirmed: ready to write the M2 plan doc.

## Hand-off to M2

_To be filled in when M1 completes._
