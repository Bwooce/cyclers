# M0 — Scaffold (todo)

Working checklist for the M0 milestone. Detailed rationale in [plan.md](plan.md). Tick boxes as work completes; do not delete items — leave the trail.

## Repo init

- [x] `git init` in the repo root (if not already a repo).
- [x] Add `nohup.out` to `.gitignore` early so it doesn't get committed.
- [x] Initial commit: `docs: project spec + overview + M0 plan` (just the `docs/` tree). _Realised as `docs: project spec, overview, and phase plans for M0-M3` once the M1/M2/M3 plans were also written in this slice — single planning-baseline commit rather than docs-then-plans._

## Top-level config files

- [x] Create `.python-version` containing `3.11`.
- [x] Create `.gitignore` (standard Python + `out/` + `nohup.out` + `.venv/`).
- [x] Create `README.md` (one-screen intro per plan §3.5).

## Packaging

- [x] Create `pyproject.toml` with:
  - [x] `[project]` block (name, version `0.0.1`, requires-python `>=3.11,<3.12`, description from spec §1).
  - [x] `[project.optional-dependencies] dev = [...]` (pytest, pytest-cov, ruff, mypy, numpy).
  - [x] `[build-system]` using `hatchling`.
  - [x] `[tool.hatch.build.targets.wheel] packages = ["src/cyclerfinder"]`.
  - [x] `[tool.ruff]` and `[tool.ruff.lint]` per plan §3.1.
  - [x] `[tool.mypy]` strict per plan §3.1.
  - [x] `[tool.pytest.ini_options]` per plan §3.1.
- [x] Run `uv sync --all-extras` locally to generate `uv.lock`. Verify it produces a `.venv/`.
- [x] Confirm `uv run python -c "import cyclerfinder"` works (will fail until source tree exists — that's next).

## Source tree

- [x] `mkdir -p src/cyclerfinder/core tests`.
- [x] `src/cyclerfinder/__init__.py` — exports `__version__ = "0.0.1"`.
- [x] `src/cyclerfinder/core/__init__.py` — empty.
- [x] `src/cyclerfinder/core/constants.py` — universal constants + `PlanetData` dataclass + `PLANETS` dict + `SAFE_PERIHELION_KM` per plan §4. Module docstring cites sources (IAU 2015, IAU 2012, JPL DE441).
- [x] `tests/__init__.py` — empty.
- [x] `tests/test_constants.py` — six assertions per plan §5.

## Local green

- [x] `uv run pytest` → green. (6 passed in 0.01s.)
- [x] `uv run ruff check .` → clean.
- [x] `uv run ruff format --check .` → clean (run `uv run ruff format .` if it isn't). (5 files already formatted.)
- [x] `uv run mypy src tests` → clean. (Resolve any strict-mode complaints — typical: missing return types on test functions, untyped dataclass defaults.) (Success: no issues found in 5 source files. No tweaks were needed — the `from __future__ import annotations` + `Final` + `dataclass(frozen=True)` choices in `constants.py` and explicit `-> None` returns in the tests satisfied strict mode on the first pass.)

## CI

- [x] `mkdir -p .github/workflows`.
- [x] `.github/workflows/ci.yml` per plan §3.4.
- [x] Commit: `m0: scaffold package, lint, type-check, test, CI`.
- [ ] Push to GitHub remote (if/when remote is configured). Confirm Actions run green. _Deferred: no remote is configured on this clone yet. User will run `gh repo create` (or equivalent) outside the M0 closeout._

## Closeout

- [x] Update `docs/overview.md` §4 milestone table: M0 status `planned` → `completed`. M1 row `not yet planned` → `planned`.
- [x] Append a `## Hand-off to M1` section to this todo.md noting:
  - Anything that didn't go as the plan predicted (e.g. mypy strict needed a tweak we didn't anticipate).
  - Any decisions deferred during M0 that M1 needs to make.
  - Confirmed: ready to write the M1 plan doc.

## Hand-off to M1

**Validation gates — all green on first try.** No fixes needed.

```
ruff check        clean
ruff format       clean (5 files already formatted)
mypy --strict     Success: no issues found in 5 source files
pytest            6 passed in 0.01s
```

### Deviations from plan (worth knowing for M1)

1. **`mean_motion_deg_day` is derived, not hand-copied.** Plan §4.2 listed mean motions as a column in the per-planet table with copy-paste values (Venus 1.602131, Earth 0.985609, Mars 0.524033). The implementation in `core/constants.py` derives `n` from each planet's `sma_au` and `MU_SUN_KM3_S2` via Kepler's third law in a private `_mean_motion_deg_day()` helper, so the three planets' mean motions are mathematically consistent with the AU and μ_Sun we adopted, rather than independently sourced numbers that could drift apart in the next digit. The plan-documented values still round-trip within the test's tolerance. This is a strict improvement on the plan but worth noting if M1 ever wants to compare against an external `n` it should expect ≤ 1 part-in-10^5 differences from textbook tables. The docstring in `constants.py` calls this out explicitly.

2. **The planning-baseline commit covered M0-M3 plans, not just M0.** The plan's todo §1 imagined committing the docs early as `docs: project spec + overview + M0 plan`, before writing the M1+ plans. In practice all four phase plans were written up-front in the same planning pass, so the single planning baseline commit is `docs: project spec, overview, and phase plans for M0-M3`. The principle (planning is committed separately from code) is preserved.

3. **Constants module also exposes `SAFE_PERIHELION_KM` (per plan §4.3).** No deviation here — flagged for visibility because M2's flyby module will consume it directly.

### Open decisions M0 deferred to M1

- **Exact pin for `lamberthub`.** M0 plan §7 notes lamberthub will be added as a dev-extra in M1 for cross-checking. The version is unpinned in the plan; M1 should pick the latest stable at the time of writing and pin it.
- **`scipy` version floor.** Same situation — runtime dep added in M1, version floor TBD. Suggest `scipy>=1.13` (current stable, plays nicely with numpy 2.x which is what M0 pinned).
- **numpy typing plugin for mypy.** M0 plan §3.1 explicitly left `plugins = []` in `[tool.mypy]` because M0 used no numpy arrays. M1 will pass numpy arrays around (Lambert state vectors, etc.). M1 should decide whether to install `numpy.typing` plugin support or use looser `numpy.typing.NDArray[np.float64]` annotations without a plugin. Recommended: try without the plugin first, since strict mode has handled the `Final[float]` patterns cleanly so far.
- **Universal-variable Lambert vs Battin / Izzo.** The spec calls for "universal variable" Lambert; M1's plan doc will need to commit to a specific formulation before code is written.

### Remote / CI status

- CI workflow file (`.github/workflows/ci.yml`) lands on disk in the M0 scaffold commit. No GitHub remote is configured yet, so Actions has not actually run. First Actions run will happen when the user pushes to a remote — at that point the four checks (ruff check, ruff format --check, mypy, pytest) should all be green since they all pass locally on the committed tree.

**Status: ready for the M1 plan to be acted on.**
