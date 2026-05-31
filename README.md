# cyclerfinder

Systematically find, rank, and verify planetary cycler trajectories — heliocentric
orbits that repeatedly re-encounter two or more planets on a fixed schedule,
maintained by gravity assists with little or no propellant. Primary targets are
Earth–Mars cyclers (used as validation) and the under-explored Venus–Earth–Mars
(VEM) family.

## Status

**M0 scaffold.** The package installs, lints, type-checks, and runs a smoke test
over the physical constants. No domain logic yet — `ephemeris`, `lambert`,
`kepler`, flyby and cycler construction land in M1–M3.

## Read first

- [docs/spec.md](docs/spec.md) — the canonical project specification.
- [docs/overview.md](docs/overview.md) — decisions made during planning and the
  milestone roadmap.
- [docs/phases/m0-scaffold/plan.md](docs/phases/m0-scaffold/plan.md) — what M0
  delivers and how.

## Local setup

Requires [`uv`](https://github.com/astral-sh/uv). The repo pins Python 3.11
via `.python-version`; `uv` will install it on first sync if it is missing.

```sh
uv sync --all-extras
uv run pytest
uv run ruff check .
uv run ruff format --check .
uv run mypy src tests
```

`--all-extras` includes the optional `viz` group (matplotlib), which enables
the Tisserand contour plotting helper in `search/tisserand.py`. The core
code does not require it; the plotting function imports matplotlib lazily.

CI runs the same four checks on every push (see `.github/workflows/ci.yml`).
