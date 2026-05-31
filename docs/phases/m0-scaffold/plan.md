# M0 — Scaffold

**Spec reference:** spec.md §4 (architecture), §7 (tech stack), §8 (M0 milestone definition).

**Purpose:** stand up a runnable Python package with packaging, lint, type-check, test, and CI infrastructure, plus the physical constants that the rest of the codebase will import. M0 produces no domain logic beyond constants — it produces the *ground* on which M1+ will build.

**Gate (definition of done):** `git clone` → `uv sync` → `uv run pytest` is green; `uv run ruff check` is clean; `uv run mypy src tests` is clean; GitHub Actions runs all three on push and passes.

---

## 1. What this milestone delivers

A repository that:

1. Installs with `uv sync` on a fresh machine pinned to Python 3.11.
2. Exposes `cyclerfinder` as an importable package using the `src/` layout.
3. Provides `cyclerfinder.core.constants` with every physical constant M1–M3 will need (Sun + Venus + Earth + Mars: GM, equatorial radius, semi-major axis, mean motion, default safe-flyby altitude).
4. Passes a smoke test confirming constants exist and pass sanity checks (e.g. AU within IAU value, planet GM positive, Mars semi-major axis ∈ [1.4, 1.6] AU).
5. Runs `ruff check`, `mypy --strict`, and `pytest` on every push via GitHub Actions.

**Explicitly out of scope for M0:** ephemeris, Lambert solver, Kepler propagator, flyby logic, anything Tisserand-related, any `Cycler` dataclass. Those are M1–M3.

---

## 2. File tree after M0

```
cyclers/
├── .github/
│   └── workflows/
│       └── ci.yml
├── .gitignore
├── .python-version                  # "3.11"
├── README.md                        # one-screen project intro + setup
├── pyproject.toml                   # package + tool config (ruff, mypy, pytest)
├── uv.lock                          # committed lockfile
├── docs/                            # already exists (spec, overview, this plan)
│   └── …
├── src/
│   └── cyclerfinder/
│       ├── __init__.py              # version string
│       └── core/
│           ├── __init__.py
│           └── constants.py
└── tests/
    ├── __init__.py
    └── test_constants.py
```

Directories for `search/`, `model/`, `verify/`, `data/`, `viz/` are **not** created in M0. They appear when the first module of each lands (per *what exists = what works*). The spec §4 architecture is the target, not a requirement to scaffold dead directories.

---

## 3. Tool configuration

### 3.1 `pyproject.toml` (PEP 621)

Sections:

- `[project]`: name `cyclerfinder`, version `0.0.1`, Python `>=3.11,<3.12`, description from spec §1.
- `[project.dependencies]`: M0 needs only `numpy>=2.0` (used by constants typing). Heavier deps (`scipy`, `astropy`, `lamberthub`, `matplotlib`, `pyyaml`) added as their phases arrive — see §7 below.
- `[project.optional-dependencies]`: `dev = ["pytest>=8", "pytest-cov>=5", "ruff>=0.6", "mypy>=1.11", "numpy>=2.0"]`.
- `[build-system]`: `hatchling` (uv-native, simplest).
- `[tool.hatch.build.targets.wheel]`: `packages = ["src/cyclerfinder"]`.
- `[tool.ruff]`: target-version `py311`, line-length `100`. `[tool.ruff.lint]` selects `E,F,W,I,N,UP,B,SIM,RUF` — standard sensible set, errors + flake8-bugbear + pyupgrade + import sort. `[tool.ruff.format]` left at defaults (PEP 8 compatible) — `ruff format` is the formatter and `ruff check` is the linter; both run in CI.
- `[tool.mypy]`: `python_version = "3.11"`, `strict = true`, `warn_unused_ignores = true`, `disallow_any_unimported = true`, `plugins = []` (numpy plugin not strictly required for M0; revisit when scipy/astropy arrive).
- `[tool.pytest.ini_options]`: `testpaths = ["tests"]`, `addopts = "-ra -q --strict-markers"`.

### 3.2 `.python-version`

Single line: `3.11`. Picked up by `uv` for environment creation.

### 3.3 `.gitignore`

Standard Python ignore set (`__pycache__/`, `*.egg-info/`, `.pytest_cache/`, `.mypy_cache/`, `.ruff_cache/`, `dist/`, `build/`, `.venv/`), plus `out/` (per spec §7 — results land there), plus `nohup.out` (pre-existing in the repo root).

### 3.4 `.github/workflows/ci.yml`

Single job, ubuntu-latest, Python 3.11. Steps:

1. `actions/checkout@v4`.
2. `astral-sh/setup-uv@v3` (installs uv, restores cache).
3. `uv sync --frozen --all-extras` — installs from `uv.lock`, fails if drift.
4. `uv run ruff check .`.
5. `uv run ruff format --check .`.
6. `uv run mypy src tests`.
7. `uv run pytest`.

Triggers: `push`, `pull_request` on all branches. No matrix in M0 (single Python version).

### 3.5 `README.md`

One screen:
- Title and one-sentence purpose (from spec §1).
- Status: M0 scaffold; not yet usable.
- "Read [docs/spec.md](docs/spec.md) for the full specification, [docs/overview.md](docs/overview.md) for sequencing."
- Local setup: `uv sync --all-extras`, `uv run pytest`.
- CI badge (optional — fine to skip until repo is on GitHub).
- License section deferred until M8 (the spec mentions dissemination but not a license).

---

## 4. `core/constants.py` content

Module-level constants only — no functions, no classes. All values typed as `Final[float]` (or `Final[dict[str, ...]]` for the planet table). Units stated in attribute docstrings/comments.

### 4.1 Universal

| Name | Value | Source | Notes |
|------|-------|--------|-------|
| `MU_SUN_KM3_S2` | `1.32712440018e11` | IAU 2015 nominal `GM_Sun` | km³/s² |
| `AU_KM` | `1.49597870700e8` | IAU 2012 definition | exact |
| `SECONDS_PER_DAY` | `86400.0` | exact | |
| `DAYS_PER_JULIAN_YEAR` | `365.25` | exact | matches `astropy.units.yr` Julian year |

### 4.2 Planets (Venus, Earth, Mars — others added as needed)

A `Final[dict[str, PlanetData]]` keyed by short code `"V"`, `"E"`, `"M"`. `PlanetData` is a `dataclass(frozen=True)` with fields:

| Field | Unit | Description |
|-------|------|-------------|
| `name` | — | e.g. `"Venus"` |
| `code` | — | one-letter |
| `mu_km3_s2` | km³/s² | gravitational parameter |
| `radius_eq_km` | km | mean equatorial radius |
| `sma_au` | AU | mean heliocentric semi-major axis (J2000) |
| `mean_motion_deg_day` | deg/day | 360 / orbital period in days |
| `safe_alt_km` | km | min flyby altitude (atmosphere top + margin) |

Values:

| Body | μ (km³/s²) | μ source | R (km) | R source | a (AU) | a source | n (°/day) | safe alt (km) |
|------|-----------|----------|--------|----------|--------|----------|-----------|---------------|
| Venus | 3.24858592e5 | JPL DE440 | 6051.8 | IAU 2015 | 0.72333566 | Standish J2000 | 1.602131 | 300 |
| Earth | 3.98600435507e5 | IERS 2010 | 6378.137 | WGS84 | 1.00000261 | Standish J2000 | 0.985609 | 300 |
| Mars | 4.282837521e4 | JPL DE440 | 3396.19 | IAU 2015 | 1.52371034 | Standish J2000 | 0.524033 | 300 |

Each value is cited inline in the constants module docstring with the same attribution. Safe altitude 300 km is a conservative default — Aldrin's original work used 200 km; we add margin and let later phases override via config.

### 4.3 Derived helpers

A single derived dict `SAFE_PERIHELION_KM = {code: data.radius_eq_km + data.safe_alt_km for code, data in PLANETS.items()}` exposed as `Final` for use by the flyby module later. Adding this here (instead of in `flyby.py`) keeps the constants module the single source of truth for physical numbers.

### 4.4 What's deliberately not here

- No Sun radius (not needed until viz).
- No planet inclination/eccentricity (the circular-coplanar model treats both as zero; the ephemeris backend in M6 will pull these from JPL DE).
- No moons, no asteroids.

---

## 5. Tests (`tests/test_constants.py`)

Six small assertions — enough to fail loudly if a number is corrupted or a refactor breaks an import. Not enough to be a substitute for downstream physics tests.

```python
def test_au_matches_iau():
    assert constants.AU_KM == 1.49597870700e8  # exact IAU 2012

def test_sun_gm_in_range():
    # IAU 2015 nominal, allow ±0.1% drift if a later source is adopted
    assert 1.326e11 < constants.MU_SUN_KM3_S2 < 1.328e11

def test_mars_sma_in_au_range():
    assert 1.50 < constants.PLANETS["M"].sma_au < 1.55

def test_planet_codes_unique_and_match_names():
    codes = [p.code for p in constants.PLANETS.values()]
    assert len(codes) == len(set(codes))
    for code, p in constants.PLANETS.items():
        assert code == p.code

def test_all_planet_gms_positive():
    assert all(p.mu_km3_s2 > 0 for p in constants.PLANETS.values())

def test_safe_perihelion_above_planet_radius():
    for code, p in constants.PLANETS.items():
        assert constants.SAFE_PERIHELION_KM[code] > p.radius_eq_km
```

No fixtures, no parameterization, no fancy markers. M0 is sanity, not science.

---

## 6. Risks specific to M0

| Risk | Mitigation |
|------|------------|
| `uv` not installed on dev machine | README says `pipx install uv` first; CI uses `astral-sh/setup-uv` action. |
| `mypy --strict` chokes on numpy 2.0 generics | M0 doesn't use numpy types yet — only the `PlanetData` dataclass. Defer the `numpy.typing` work to M1 when arrays appear. |
| Constants drift between sources (JPL vs IAU) | Each value cites its source in the table above and the module docstring. Tests use ranges, not exact equality, for derived values. |
| Pre-existing `nohup.out` in repo root | Added to `.gitignore` so it's not committed by accident; the file itself stays untracked. |
| `hatchling` vs `setuptools` choice locks in for life | Both are fine; `hatchling` is lighter and uv-native. Switching later is a one-file change in `pyproject.toml`. |

---

## 7. Dependency policy

Add a runtime dep to `pyproject.toml` only when the first module that needs it is being written. Keeps the lockfile honest and CI install fast.

Anticipated additions (do NOT add now):

| Phase | Adds |
|-------|------|
| M1 | `scipy` (root-finding for Lambert), `lamberthub` (cross-check, dev-extra is fine) |
| M2 | `matplotlib` (Tisserand plots — optional viz extra) |
| M6 | `astropy` (JPL ephemeris) |
| M5 (decision deferred) | `pygmo` *or* stick with `scipy.optimize.differential_evolution` |
| M8 | `pyyaml` (config files) |
| Stretch | GMAT bridge runtime (subprocess only, no Python dep) |

---

## 8. Order of work

The todo.md mirrors this with checkboxes. High-level:

1. `git init` (if not already), commit the docs that already exist as the first commit (`docs: project spec + overview + M0 plan`).
2. Create `.python-version`, `.gitignore`, `README.md`.
3. Create `pyproject.toml` with the configuration above. Run `uv sync --all-extras` locally to generate `uv.lock`.
4. Create the `src/cyclerfinder/` and `tests/` layout with `__init__.py` files.
5. Write `core/constants.py`.
6. Write `tests/test_constants.py`; confirm `uv run pytest` is green locally.
7. Write `.github/workflows/ci.yml`.
8. Commit (`m0: scaffold package, lint, type-check, test, CI`).
9. Push to GitHub (when the remote exists) and confirm Actions go green.
10. Update `docs/overview.md` milestone table: M0 → completed; mark M1 as next.

---

## 9. Exit checklist (the gate, restated)

Before declaring M0 done:

- [ ] `uv run pytest` green locally on a fresh checkout.
- [ ] `uv run ruff check .` clean.
- [ ] `uv run mypy src tests` clean.
- [ ] CI workflow runs all four checks (ruff lint, ruff format, mypy, pytest) on push and passes.
- [ ] `docs/overview.md` updated: M0 status = completed.
- [ ] Hand-off note appended to `phases/m0-scaffold/todo.md` under `## Hand-off to M1` — anything unexpected encountered during scaffold that M1 needs to know.

(Writing the M1 plan doc is the first task of M1, not an M0 exit criterion.)
