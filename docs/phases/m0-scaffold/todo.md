# M0 — Scaffold (todo)

Working checklist for the M0 milestone. Detailed rationale in [plan.md](plan.md). Tick boxes as work completes; do not delete items — leave the trail.

## Repo init

- [ ] `git init` in `/home/bruce/dev/cyclers` (if not already a repo).
- [ ] Add `nohup.out` to `.gitignore` early so it doesn't get committed.
- [ ] Initial commit: `docs: project spec + overview + M0 plan` (just the `docs/` tree).

## Top-level config files

- [ ] Create `.python-version` containing `3.11`.
- [ ] Create `.gitignore` (standard Python + `out/` + `nohup.out` + `.venv/`).
- [ ] Create `README.md` (one-screen intro per plan §3.5).

## Packaging

- [ ] Create `pyproject.toml` with:
  - [ ] `[project]` block (name, version `0.0.1`, requires-python `>=3.11,<3.12`, description from spec §1).
  - [ ] `[project.optional-dependencies] dev = [...]` (pytest, pytest-cov, ruff, mypy, numpy).
  - [ ] `[build-system]` using `hatchling`.
  - [ ] `[tool.hatch.build.targets.wheel] packages = ["src/cyclerfinder"]`.
  - [ ] `[tool.ruff]` and `[tool.ruff.lint]` per plan §3.1.
  - [ ] `[tool.mypy]` strict per plan §3.1.
  - [ ] `[tool.pytest.ini_options]` per plan §3.1.
- [ ] Run `uv sync --all-extras` locally to generate `uv.lock`. Verify it produces a `.venv/`.
- [ ] Confirm `uv run python -c "import cyclerfinder"` works (will fail until source tree exists — that's next).

## Source tree

- [ ] `mkdir -p src/cyclerfinder/core tests`.
- [ ] `src/cyclerfinder/__init__.py` — exports `__version__ = "0.0.1"`.
- [ ] `src/cyclerfinder/core/__init__.py` — empty.
- [ ] `src/cyclerfinder/core/constants.py` — universal constants + `PlanetData` dataclass + `PLANETS` dict + `SAFE_PERIHELION_KM` per plan §4. Module docstring cites sources (IAU 2015, IAU 2012, JPL DE441).
- [ ] `tests/__init__.py` — empty.
- [ ] `tests/test_constants.py` — six assertions per plan §5.

## Local green

- [ ] `uv run pytest` → green.
- [ ] `uv run ruff check .` → clean.
- [ ] `uv run ruff format --check .` → clean (run `uv run ruff format .` if it isn't).
- [ ] `uv run mypy src tests` → clean. (Resolve any strict-mode complaints — typical: missing return types on test functions, untyped dataclass defaults.)

## CI

- [ ] `mkdir -p .github/workflows`.
- [ ] `.github/workflows/ci.yml` per plan §3.4.
- [ ] Commit: `m0: scaffold package, lint, type-check, test, CI`.
- [ ] Push to GitHub remote (if/when remote is configured). Confirm Actions run green.

## Closeout

- [ ] Update `docs/overview.md` §4 milestone table: M0 status `planned` → `completed`. M1 row `not yet planned` → `planned`.
- [ ] Append a `## Hand-off to M1` section to this todo.md noting:
  - Anything that didn't go as the plan predicted (e.g. mypy strict needed a tweak we didn't anticipate).
  - Any decisions deferred during M0 that M1 needs to make.
  - Confirmed: ready to write the M1 plan doc.

## Hand-off to M1

_To be filled in when M0 completes._
