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
| M6b — real-ephemeris closure; powered Aldrin solver landed; `optimise_cell_ephemeris` (real-DE440 cell optimiser) implemented with asymmetric `tof_seed_days`; drift closure proven physically unreachable for k=1 (retargeted each synodic period w/ maintenance ΔV) | done |
| M7 — catalogue loader, signature matching, novelty (crosscheck + writeback + discover landed; `discover(optimiser="ephemeris")` wired) | done |
| M8-Core — VEM 3-body search core (`period_basis` beat dispatch, same-body Tisserand bypass, `CONSTRUCTIBLE_MULTIBODY` admission, sourced 12.8-yr Jones gate) | ✓ done |
| M8-UX — VEM campaign CLI + viz + reporting | planned |
| Low-thrust (v2) — Sims-Flanagan model: leg model, feasibility/NLP constraints, two-phase DE+SLSQP solve, powered-maintenance evaluator (machinery only — no sourced powered rows exist) | ✓ done |

**Companion catalogue** at [`data/catalogue.yaml`](data/catalogue.yaml) carries the
237-entry published-cycler seed library (Aldrin family, Russell-Ocampo Table 3.4, McConaghy SnLm
broad classes, Niehoff VISIT, Jones VEM family, the 15-orbit Hollister–Menning Earth–Venus
family, plus lunar and Jovian family seeds). Every numerical value carries a source quote per [`data/README.md`](data/README.md)
conventions. Real-ephemeris launch windows for each ballistic Earth-touching entry are
auto-published to <https://cyclers.space/launch-windows/> (weekly cron sync).

## Read first

- [docs/spec.md](docs/spec.md) — the canonical project specification.
- [docs/overview.md](docs/overview.md) — decisions made during planning and the
  milestone roadmap.
- [docs/phases/](docs/phases/) — per-phase plan docs; the active phase is m8
  (`m8-multibody-vem/plan.md`). Completed-phase working docs (M0–M7) are retired
  into the roadmap + spec; their `plan.md` files remain as milestone history.

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
