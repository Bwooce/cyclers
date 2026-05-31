# Cyclerfinder — Project Overview

This is the navigation document for the cyclerfinder project. The canonical specification is **[spec.md](spec.md)** (authored by the project owner, preserved verbatim). This overview captures the decisions made during planning, the milestone roadmap as actually sequenced, and the structure of the per-phase planning documents.

If you are picking up the project cold, read in this order:
1. **spec.md** — what we're building and why.
2. **This overview** — sequencing, decisions, document layout.
3. **phases/m0-scaffold/plan.md** (or whichever phase is current) — the immediate work.

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
| **M0** | Scaffold | completed | [phases/m0-scaffold/plan.md](phases/m0-scaffold/plan.md) |
| **M1** | Core mechanics: ephemeris, lambert, kepler | completed | [phases/m1-core-mechanics/plan.md](phases/m1-core-mechanics/plan.md) |
| **M2** | Flyby + maps: flyby, tisserand, resonance | completed | [phases/m2-flyby-maps/plan.md](phases/m2-flyby-maps/plan.md) |
| **M3** | Model + construction: cycler, frames, construct; reproduce Aldrin | completed | [phases/m3-model-construct/plan.md](phases/m3-model-construct/plan.md) |
| **M4** | Enumeration + scoring | planned | — |
| **M5** | Optimisation (rediscover 2-synodic E–M from scratch) | future slice | — |
| **M6a** | Idealized closure verification | future slice | — |
| **M6b** | Phase-match + ephemeris-mode TCM over 3–5 laps | future slice | — |
| **M7** | Catalogue, signatures, novelty | future slice | — |
| **M8** | VEM campaign + CLI + viz | future slice | — |
| Stretch | GMAT bridge | future slice | — |
| Stretch | Low-thrust | future slice | — |
| Stretch | `cyclers.space` public site | future slice | — |

Per the agreed first slice, **only M0's plan is written now**. M1–M3 plans are written as each previous phase completes, so each plan can incorporate what was actually learned from the predecessor (file boundaries that emerged, types that turned out wrong, etc.). Writing all four up-front risks four plans built on assumptions that get invalidated by M0.

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
| Aldrin cycler semi-major axis | a ≈ 1.659 AU | M3 gate |
| Aldrin cycler eccentricity | e ≈ 0.41 | M3 gate |
| Aldrin perihelion / aphelion | ≈ 0.98 / 2.34 AU | M3 gate |
| Aldrin E→M leg time-of-flight | ≈ 146 d | M3 gate |
| Published 2-synodic E–M V∞ | ≈ 5.65 km/s (E), 3.05 km/s (M) | M5 gate |
| Degenerate-solution guard | reject "ballistic" closure if V∞ > ~11 km/s | M5 constraint |

---

## 6. Document layout

```
docs/
├── spec.md                      # Authoritative project spec (owner-authored, verbatim)
├── overview.md                  # This file: navigation + decisions + roadmap status
└── phases/
    └── m0-scaffold/
        ├── plan.md              # Detailed implementation plan for the milestone
        └── todo.md              # Actionable checklist derived from plan.md
```

Each milestone that gets started gets a `phases/<id>-<short-name>/` directory containing:
- `plan.md` — module-level design decisions for that milestone, file structure, tests, gates, risks.
- `todo.md` — concrete actionable items, checked off as work progresses.

Closed milestones keep their docs in place as a record. The plan documents are append-only after the milestone completes (notes added if subsequent work reveals issues, but the original plan stands as it was).

---

## 7. How to use this repo as an implementing agent

1. Read **spec.md** end-to-end. The spec includes the architecture, algorithm pipeline, and validation anchors that constrain every implementation choice.
2. Read **overview.md** (this file) for what's already decided and what's deferred.
3. Open the current phase's **plan.md** + **todo.md**. Work the todo in order; mark items as done.
4. If you hit a design ambiguity not resolved by spec.md or the phase plan, **ask** before guessing — the spec's §10 risks list calls out failure modes that "reasonable guesses" have historically triggered (notably blind DE on a guessed sequence).
5. Every milestone has a **gate**: a hard numeric check from spec.md §9. Don't claim a milestone done until its gate passes in `pytest`.
6. Commit per milestone, not per file. The git history should read as the project's structural progress.
