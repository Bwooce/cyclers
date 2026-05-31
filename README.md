# cyclerfinder

Systematically find, rank, and verify planetary cycler trajectories — heliocentric
orbits that repeatedly re-encounter two or more planets on a fixed schedule,
maintained by gravity assists with little or no propellant. Primary targets are
Earth–Mars cyclers (used as validation) and the under-explored Venus–Earth–Mars
(VEM) family.

## Status

| Milestone | Status |
|---|---|
| M0 — scaffold (uv + pyproject + ruff + mypy --strict + pytest + CI) | ✓ done |
| M1 — core mechanics (ephemeris circular, lambert universal-variable, kepler) | ✓ done |
| M2 — flyby + Tisserand + resonance | ✓ done |
| M3 — Cycler/Leg/Encounter model + rotating frame + construct; reproduces Aldrin | ✓ done |
| M4 — cell enumeration + Tisserand pruning + scoring + ranking | ✓ done |
| M5 — optimisation (scipy DE + SLSQP with hard constraints) | in progress |
| M6 (slice) — astropy Ephemeris backend + phase_match.find_real_windows | ✓ done |
| M6a/M6b — multi-lap propagation, ephemeris-mode TCM minimisation | planned |
| M7 — catalogue loader, signature matching, novelty | planned |
| M8 — VEM campaign + CLI + viz | planned |

**Companion catalogue** at [`data/seed_cyclers.yaml`](data/seed_cyclers.yaml) carries the
published-cycler seed library (Aldrin family, Russell-Ocampo Table 3.4, McConaghy SnLm
broad classes, Niehoff VISIT, Jones VEM family, Hollister–Menning, plus lunar and Jovian
family seeds). Every numerical value carries a source quote per [`data/README.md`](data/README.md)
conventions. Real-ephemeris launch windows for each ballistic Earth-touching entry are
auto-published to <https://cyclers.space/launch-windows/> (weekly cron sync).

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
