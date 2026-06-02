# Cyclerfinder — Project Overview

This is the navigation document for the cyclerfinder project. The canonical specification is **[spec.md](spec.md)** (authored by the project owner, preserved verbatim). This overview captures the decisions made during planning, the milestone roadmap as actually sequenced, and the structure of the per-phase planning documents.

If you are picking up the project cold, read in this order:
1. **spec.md** — what we're building and why.
2. **This overview** — sequencing, decisions, document layout.
3. **phases/m6b-real-ephemeris-closure/plan.md** / **phases/m7-catalogue-novelty-matching/plan.md** (the active phases) — the immediate work.

---

## 1. What we're building

A tool that **systematically finds, ranks, and verifies planetary cycler trajectories** — heliocentric orbits that repeatedly re-encounter two or more planets on a fixed schedule, maintained by gravity assists with little or no propellant.

Primary targets:
- **Earth–Mars cyclers** (well-documented, used as validation).
- **Venus–Earth–Mars (VEM)** triple cyclers (genuinely under-explored).

The structural insight is that blind optimisation does not find ballistic cyclers — you must enumerate **energetically feasible flyby structures first** (Tisserand-pruned), then construct V∞-matched legs, then optimise inside each surviving structure. See spec.md §10 and §13 for the failure modes this design is built to avoid.

---

## 2. Decisions made during planning

These supplement spec.md with the concrete choices needed to start building.

| Decision | Choice | Rationale |
|----------|--------|-----------|
| **First slice** | M0–M3 (scaffold → reproduce Aldrin cycler) | Meaningful end-to-end deliverable. Validates the foundation against a hard published number. No optimiser yet, so the slice is bounded. |
| **Prior prototype reuse** | Rewrite from scratch against spec §9 anchors | Prototype scripts (`lambert.py`, `scan*.py`) not accessible. Clean build against validation anchors. |
| **Python version** | 3.11 (pinned) | Per spec §7. Keeps the pykep/pygmo door open for M5+. |
| **Toolchain** | `uv` + `pyproject.toml` + `src/` layout | Fast, lockfile-backed reproducibility (uv.lock committed). `src/` layout prevents accidental imports of unbuilt package. |
| **Lambert interface** | §12(b) multi-rev list-return from day one | `lambert(...) -> list[LambertSolution]` with `max_revs` parameter. M1 only exercises `max_revs=0`, but the interface is stable through M4. |
| **Lint / format** | `ruff` (`ruff check` + `ruff format`) | One tool for both, fast, no config sprawl. Both run in CI. |
| **Type checking** | `mypy --strict` | Strict from day one; numerical signatures benefit from explicit array shape contracts. Pay the discipline cost early. |
| **Test framework** | `pytest` | Per spec §7. |
| **CI** | GitHub Actions on push | Standard; runs ruff + mypy + pytest. Detailed config in M0 plan. |
| **Visual companion** | Not used for planning docs | Spec was prose; planning artifacts are prose. Re-evaluate when we start designing UI/viz in M8. |

Deferred decisions (revisit when their phase begins):

- **Ledger backend** — SQLite vs JSONL (decided in M4/M7 when ledger is built).
- **Global optimiser** — scipy DE vs pygmo (decided in M5).
- **GMAT bridge transport** — file generation strategy (decided in stretch/V4 work).
- **Public site stack** — static generator choice (decided when M8 completes and dissemination begins).

---

## 3. Scope of the first slice (M0–M3)

What's **in**:
- Project scaffold, packaging, CI, test harness.
- Constants (μ, AU, planet GM/radii/a).
- Circular-coplanar ephemeris.
- Lambert solver (universal-variable, multi-rev list-return) with lamberthub cross-check.
- Kepler propagator (universal-variable).
- Flyby max-bend, feasibility, powered-ΔV.
- Tisserand graph (contours + `linkable(body_a, body_b, vinf)` predicate).
- Resonance: synodic periods, k-multiples, VEM beat.
- `Cycler` / `Leg` / `Encounter` dataclasses; closure residual.
- Uniform rotating frame (for circular-coplanar closure check).
- Patched-conic leg construction; **reproduce Aldrin cycler** as the gate.

What's **out** of this slice (later milestones):
- Sequence enumeration (M4).
- Global + local optimisation (M5).
- Real ephemeris backend, multi-lap verification, phase-matching (M6).
- Catalogue, signatures, novelty (M7).
- VEM campaign, CLI, viz (M8).
- GMAT bridge, low-thrust, public site (stretch).

The cli skeleton from spec §4 is scaffolded in M0 but stays empty until M8 — there's nothing to expose as commands until then.

---

## 4. Milestone roadmap

Authoritative milestone definitions and gates: **spec.md §8**. Planning status below.

| ID | Title | Status | Plan doc |
|----|-------|--------|----------|
| **M0** | Scaffold | completed | — |
| **M1** | Core mechanics: ephemeris, lambert, kepler | completed | — |
| **M2** | Flyby + maps: flyby, tisserand, resonance | completed | — |
| **M3** | Model + construction: cycler, frames, construct; reproduce Aldrin | completed | — |
| **M4** | Enumeration + scoring | completed | — |
| **M5** | Optimisation (rediscover 2-synodic E–M from scratch) | **completed** (commit `e6412c4`); `test_multi_start_grid_distinct` xfail diagnosed as load-bearing collision (task #53); binding gate `test_2syn_em_rediscovers_5_65_kms_earth` is a pre-existing slow-marked regression (task #54) | — |
| **M6 (slice)** | astropy `Ephemeris` backend + `phase_match.find_real_windows` (geometric launch-window dates); ICRS→ecliptic rotation fix `4fe901d` (2026-06-01) | **completed** (commits `9b2611d`, `4fe901d`) | — |
| **M6a** | Idealized closure verification (multi-lap propagation in dynamic rotating frame, bounded closure-drift check; `DRIFT_TOLERANCE_KM = 50_000`) | **completed** (commits `ba01f37` + `1852750`, 2026-06-01); 7/7 binding gates pass, dynamic-frame round-trip 1e-10 rel, circular-Aldrin drift ~3e-7 km | — |
| **M6b** | Real-ephemeris closure verification; Lambert-chain across DE440 with `REAL_DRIFT_TOLERANCE_KM = 200_000` over 2 cycles; Pascarella 2024 template | **done**; powered-cycler solver landed (`solve_powered_periodic_cycler`, commit `83f6272`); `optimise_cell_ephemeris` (real-DE440 cell optimiser) implemented with asymmetric `tof_seed_days` + Aldrin parity test (was a stub). **Proven:** rotating-frame drift closure is physically unreachable for k=1 Aldrin on DE440 (Mars heliocentric radius breathes ≈0.117 AU/cycle; measured drift ≈7.24e7 km ≈ 362× tol) — this is *why* the real Aldrin cycler is retargeted each synodic period with maintenance ΔV; literature never claims zero-maintenance periodicity. Multi-rev Lambert (S1L1) remains stretch. | [phases/m6b-real-ephemeris-closure/plan.md](phases/m6b-real-ephemeris-closure/plan.md) |
| **M7** | Catalogue loader, canonical signature matching, novelty scoring (spec §8 + §14 V1); also absorbs M6a-deferred `verify/crosscheck.py` + the spec §13.6/§13.8 JSONL ledger | **done** — `data/catalog.py`, `ledger.py`, `writeback.py`, `discover.py` (with `optimiser="ephemeris"`), `verify/crosscheck.py` all shipped; 233-entry census frozen as a ratchet | [phases/m7-catalogue-novelty-matching/plan.md](phases/m7-catalogue-novelty-matching/plan.md) |
| **M8** | VEM campaign + CLI + viz | planned (DEFERRED by user decision; plan written + revised) | [phases/m8-multibody-vem/plan.md](phases/m8-multibody-vem/plan.md) |
| Live | `cyclers.space` public site — catalogue browser + planet filter + real-ephemeris launch windows | **shipped** ([cyclers.space](https://cyclers.space)) | — |
| Stretch | GMAT bridge (V4 of validation gauntlet) | planned | — |
| Stretch | Low-thrust v2 scope expansion | refs queued (Yam 2010, Pascarella 2024) | — |

Completed-phase `todo.md` working checklists (M0–M7) have been retired; their durable outcomes live in this table, spec.md, and the code itself. Each milestone's `plan.md` is kept as history. The active phase is M8 (`phases/m8-multibody-vem/plan.md`, currently deferred).

### 4.1 Schema-v2 + representation framework (2026-06-01)

A major catalogue / spec rev shipped 2026-06-01 (commits `debd285` schema v2 + `5145bb1` spec §12.2):

- **Catalogue schema v2** (`data/catalogue.yaml`): six additive optional field categories — `model_assumption`, `flyby_mechanics[]`, `delta_v_kms` + `v_infinity_leveraging_dv_kms`, `periapse_km`/`apoapse_km`, RAAN/ω/ν/epoch, `fleet_size`. All 233 entries mechanically backfilled where derivable; remaining nulls are honest gaps. See `data/README.md` "Schema v2" for full convention. spec.md §16.1 carries the JSON schema; §16.2 explicitly carves the new fields OUT of the canonical signature so existing matches stay valid.
- **Representation framework** (spec §12.2): explicit architectural rule that the catalogue carries multiple legitimate representations per cycler (idealised circular-coplanar, real-ephemeris instances on the site, analytic-ephemeris from Rogers 2012, CR3BP from Arenstorf/Saturnian). Downstream consumers pick the right model via the `model_assumption` field; we deliberately do NOT standardise on one form. "Derived views, not conversion."
- **Future scope references** (`docs/v2-future-references.md`): bibliography for the v2 low-thrust scope expansion (Yam 2010 Sims-Flanagan, Pascarella 2024 medium-fidelity pipeline, Izzo 2015 GTOC methods, Burhani 2023 inclined low-thrust, Hollister–Menning 1969-71 foundational lineage).

---

## 5. Validation anchors (compact reference)

These are the hard numbers any implementation must hit. Authoritative list in spec.md §9.

| Anchor | Value | Used in |
|--------|-------|---------|
| Lambert vs lamberthub (izzo, gooding) | max Δ\|v\| < 1e-3 m/s on test legs | M1 gate |
| Earth–Mars synodic | 2.135 yr | M2 gate |
| Earth–Venus synodic | 1.599 yr | M2 gate |
| VEM beat (3× E–M ≈ 4× E–V) | ≈ 6.40 yr | M2 gate |
| Mars max bend at V∞ = 7 km/s | ≈ 24° | M2 gate |
| Earth/Venus max bend at V∞ = 7 km/s | ≈ 60–63° | M2 sanity |
| Aldrin cycler semi-major axis | a ≈ 1.60 AU | M3 gate |
| Aldrin cycler eccentricity | e ≈ 0.393 | M3 gate |
| Aldrin perihelion / aphelion | ≈ 0.97 / 2.23 AU | M3 gate |
| Aldrin E→M leg time-of-flight | ≈ 146 d | M3 gate |
| Published 2-synodic E–M V∞ | ≈ 5.65 km/s (E), 3.05 km/s (M) | M5 gate |
| Degenerate-solution guard | reject "ballistic" closure if V∞ > ~11 km/s | M5 constraint |

> **Aldrin anchor reconciliation:** the values above (a ≈ 1.60 AU, e ≈ 0.393, perihelion/aphelion ≈ 0.97/2.23 AU) are the literature anchors from Byrnes/Longuski/Aldrin 1993 as tabulated in Rogers et al. 2012, and are what the M3 gate tests assert. Earlier drafts cited a ≈ 1.659 / e ≈ 0.41 — a resonance-construction artifact; see spec.md §9.1 for the full reconciliation.

---

## 6. Document layout

```
docs/
├── spec.md                      # Authoritative project spec (owner-authored, verbatim)
├── overview.md                  # This file: navigation + decisions + roadmap status
├── errata-investigation.md      # Full record of the Aldrin anchor reconciliation (spec §9.1 defers here)
├── v2-future-references.md      # Seed bibliography for the v2 low-thrust scope expansion
├── refs/                        # Cached primary-source PDFs (Russell 2004, Rogers 2012, Hollister-Menning 1970, ...)
└── phases/
    └── <id>-<short-name>/        # Per-phase docs; active phase is m8 (m6b/m7 retain plan.md as history)
        └── plan.md              # Detailed implementation plan for the milestone
```

Each milestone gets a `phases/<id>-<short-name>/` directory whose `plan.md`
carries module-level design decisions, file structure, tests, gates, and risks.
The per-phase `todo.md` working checklists are deleted once the milestone ships
(their outcome lives in the roadmap table, spec.md, and the code).

Once a milestone is complete, its plan/todo working docs are retired — the durable outcomes live in the milestone table above, spec.md, and the code. Only active phases keep their phase docs.

---

## 7. How to use this repo as an implementing agent

1. Read **spec.md** end-to-end. The spec includes the architecture, algorithm pipeline, and validation anchors that constrain every implementation choice.
2. Read **overview.md** (this file) for what's already decided and what's deferred.
3. Open the current phase's **plan.md** + **todo.md**. Work the todo in order; mark items as done.
4. If you hit a design ambiguity not resolved by spec.md or the phase plan, **ask** before guessing — the spec's §10 risks list calls out failure modes that "reasonable guesses" have historically triggered (notably blind DE on a guessed sequence).
5. Every milestone has a **gate**: a hard numeric check from spec.md §9. Don't claim a milestone done until its gate passes in `pytest`.
6. Commit per milestone, not per file. The git history should read as the project's structural progress.
