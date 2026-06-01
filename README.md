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
| M5 — optimisation (scipy DE + SLSQP with hard constraints) | ✓ done |
| M6 (slice) — astropy Ephemeris backend + phase_match.find_real_windows | ✓ done |
| M6a — idealized closure verification (multi-lap rotating-frame drift check) | ✓ done |
| M6b — real-ephemeris closure; powered Aldrin solver landed; drift closure proven physically unreachable for k=1 (retargeted each synodic period w/ maintenance ΔV) | scaffolding shipped |
| M7 — catalogue loader, signature matching, novelty (crosscheck + writeback + discover landed) | partial |
| M8 — VEM campaign + CLI + viz | planned |

**Companion catalogue** at [`data/catalogue.yaml`](data/catalogue.yaml) carries the
published-cycler seed library (Aldrin family, Russell-Ocampo Table 3.4, McConaghy SnLm
broad classes, Niehoff VISIT, Jones VEM family, Hollister–Menning, plus lunar and Jovian
family seeds). Every numerical value carries a source quote per [`data/README.md`](data/README.md)
conventions. Real-ephemeris launch windows for each ballistic Earth-touching entry are
auto-published to <https://cyclers.space/launch-windows/> (weekly cron sync).

## Read first

- [docs/spec.md](docs/spec.md) — the canonical project specification.
- [docs/overview.md](docs/overview.md) — decisions made during planning and the
  milestone roadmap.
- [docs/phases/](docs/phases/) — per-phase plan/todo docs for the active phases
  (m6b, m7); completed-phase docs are retired into the roadmap + spec.

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
