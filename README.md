# cyclerfinder

Systematically find, rank, and verify planetary cycler trajectories — heliocentric
orbits that repeatedly re-encounter two or more planets on a fixed schedule,
maintained by gravity assists with little or no propellant. Primary targets are
Earth–Mars cyclers (used as validation) and the under-explored Venus–Earth–Mars
(VEM) family.

**Catalogue scope (schema v4.7, expanded 2026-06-15):** four orbit classes,
admitting the mission-actionable epoch-locked literature alongside strict cyclers.

| Class | Period? | Epoch-locked? | Returns | Type case |
|---|---|---|---|---|
| `cycler` | strictly periodic | NO | ∞ | Aldrin Earth-Mars; Russell-Ocampo S1L1 |
| `quasi_cycler` | closes-up-to-rotation | YES (10–15 yr) | 3–15 | cyclers-of-opportunity |
| `precursor_mga` | non-repeating | YES (launch window) | 1 (insertion) | one-shot insertion into an extant cycler |
| `mga_tour` | non-repeating | YES (launch window) | 1 (terminal) | Galileo VEEGA; Tito 2018 Mars free-return |

See [`docs/notes/2026-06-16-catalogue-scope-taxonomy.md`](docs/notes/2026-06-16-catalogue-scope-taxonomy.md)
for the full taxonomy, V0–V5 gauntlet extension, and migration record.

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
| M8-UX — `cyclerfinder` CLI (enumerate/solve/discover/report/viz) + campaign reports + optional-extra viz | ✓ done |
| Low-thrust (v2) — Sims-Flanagan model: leg model, feasibility/NLP constraints, two-phase DE+SLSQP solve, powered-maintenance evaluator (machinery only — no sourced powered rows exist) | ✓ done |

**Companion catalogue** at [`data/catalogue.yaml`](data/catalogue.yaml) carries the
282-entry published-cycler seed library (incl. Tito 2018 Mars free-return + Heaton-Longuski 2003 Uranian satellite tour U00-01 as the two `mga_tour` rows under schema v4.7's expanded scope) (Aldrin family, Russell-Ocampo Table 3.4, McConaghy SnLm
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
