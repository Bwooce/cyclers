# VILM-Leveraging Endgame Solver — Design Spec

**Date:** 2026-06-09
**Status:** Design approved (brainstorming), pending implementation plan.
**Sub-project:** 1 of 2 in the "genome upgrade" track. Spec 2 (multi-arc DSM into
the corrector — the catalogue lift) is a separate, later spec → plan → build cycle.

## Goal

Turn the existing V-infinity-leveraging machinery from a **feasibility filter**
into **search degrees of freedom**, so the discovery search can *lower* the
encounter V∞ of a moon-tour cycler into bend-feasible reach — reopening the
Jovian/Saturnian "empty regions" with a phase-full Campagnola–Russell endgame
solver. The product is a **powered** moon-tour cycler family (V3-powered class),
not a ballistic one.

## Why this, why now

The Forge Phase 6 ballistic sweeps (#172, #178) are rigorous **bounded
negatives**: the no-leveraging single-ellipse genome closes only high-V∞ moon
families that are **bend-infeasible** — best achieved V∞ 18.98–27.16 km/s against
a ~6.0 km/s floor (12.98–21.16 km/s gaps). The closest-to-floor closure seen was
the Saturnian Enceladus-Dione-Rhea family at 7.52 km/s — still infeasible but the
smallest gap. The blocker is a **genome** limitation, not infrastructure: the
search cannot actuate the V∞-leveraging that real moon-tour designs use to walk
V∞ down. This is the **no-papers discovery lever** (see memory
`validation-ceiling`): pure engineering, no acquisitions required.

The primitives already exist but are not wired as search DOFs:
- `src/cyclerfinder/search/vilm.py` — the **phase-free** half of the
  Campagnola–Russell "Endgame Problem" Part 1 (Γ function, V̄∞ thresholds, the
  theoretical-minimum ΔV quadrature). Currently used **only** in `moon_prune.py`
  as a prune gate.
- `src/cyclerfinder/search/dsm_leg.py` — Takao one-DSM-per-leg primitive (used by
  `free_return_chain.py`, not by discovery; reserved for Spec 2).

This spec builds the **phase-full** half and wires it into discovery.

## Binding constraints (orbit-closure-discipline, non-negotiable)

These are baked into the acceptance gates, not optional:
1. A **same-model sourced golden** must be reproduced before any discovery is
   trusted (avoids the #164 S1L1 overclaim).
2. **Independent cross-check** mandatory — here, the phase-free Γ-quadrature is a
   genuinely different code path from the phase-full leg evaluator.
3. **No catalogue writeback** until n-body confirms AND a same-model golden
   exists for the family.
4. **No tolerance/budget/cap loosening** to manufacture a survivor. A clean
   negative (no route to the V̄∞ floor) is a SUCCESS, recorded as a
   method-versioned EMPTY-region report.
5. A leg/route cheaper than the Γ-quadrature analytic minimum is a **bug**, never
   a discovery.

## Sourced golden anchors

From `docs/notes/2026-06-05-endgame-tisserand-mining.md` (Campagnola & Russell,
"The Endgame Problem" Parts 1–2). All transcribed verbatim from the papers'
tables/text. Circular-coplanar model — the same model this solver uses at its
first fidelity rung.

- **A1** (Part 1, Table 3): per-moon V̄∞ E/I thresholds — the **leveraging-validity
  lower edge** (the minimum V∞ below which a VILM stops helping), NOT the route
  target. Europa **0.277 / 0.290 km/s**. Used as a leg-evaluator golden and a
  termination sanity bound — distinct from `target_vinf_floor_kms` (below).
- **A2** (Part 1, Table 1): per-transfer ΔV min/max, no gravity assist —
  Ganymede→Europa **1.71** / 2.18, Callisto→Europa 1.94 / 3.75,
  Europa→Io 1.76 / 2.54 km/s.
- **A3** (Part 1, Table 2): same, *with* intermoon gravity assists
  (e.g. Callisto-G-Europa 1.61 / 2.07).
- **A5** (Part 2, worked scalars): Ganymede→Europa quasi-ballistic budget —
  Δv_escape ≅ 0.72, Δv_capture ≅ 0.51, total **Δv_TOT ≈ 1.25 km/s**, 291+82 days;
  VILM theoretical minimum **1.71 km/s** (matches A2).
- **A4** (Part 2, Table 1): orbit-insertion ΔV — contains the **two flagged
  suspect cells; NEVER used as goldens** (asserted in the test module docstring).

Proving ground: **Jupiter** (Europa/Ganymede/Callisto/Io) — where the worked
numbers live and where the empty region sits.

## Architecture

Two new focused modules under `src/cyclerfinder/search/`, plus reuse.

### Component 1 — `leveraging_leg.py` (phase-full single-leg evaluator)

A VILM leg: depart moon M with `V∞_in`, fly an orbit resonant with M
(`n:m_K±` — n spacecraft revs : m moon revs, K the encounter/apsidal index,
± exterior/interior), apply a DSM at the apse (apoapsis `+` / periapsis `−`),
return to M with changed `V∞_out`. The apse DSM **is** the leveraging maneuver —
the leg is powered by construction.

```python
@dataclass(frozen=True)
class LeveragingLegResult:
    dv_dsm_kms: float            # CONSTRAINED — the apse leveraging burn
    vinf_out_kms: float          # EMERGED evidence
    resonance_residual: float    # phasing closure (return to M) — must be ~0
    apse_radius_km: float        # EMERGED
    taxonomy: VilmLeg            # the n:m_K± classification (from vilm.py)
    converged: bool              # by residual magnitude (mirrors dsm_leg)
    gamma_floor_ok: bool         # dv_dsm >= Γ-quadrature theoretical min

def evaluate_leveraging_leg(
    moon: str,
    taxonomy: VilmLeg,
    vinf_in_kms: float,
    vinf_out_target_kms: float,
    epoch_sec: float,
    ephem: Ephemeris,
) -> LeveragingLegResult: ...
```

Given `(n:m_K±, V∞_in, V∞_out_target)`: the resonant orbit period is fixed by the
n:m resonance with M's period (`satellites.py`); the apse radius follows from the
energy yielding `V∞_out_target`; the DSM is the velocity change at the apse
between inbound and outbound resonant arcs. Conventions mirror `dsm_leg.py`:
frozen dataclass separating CONSTRAINED from EMERGED quantities; `converged` by
residual magnitude; full audit trail. `gamma_floor_ok` asserts `dv_dsm_kms` ≥ the
`vilm.py` Γ-quadrature minimum for that V∞ band — the built-in second opinion.

**Golden gate (blocks Component 2):** reproduce A1 thresholds and A2 per-leg ΔV
minima to a few-% tolerance (tight — same circular-coplanar model).

### Component 2 — `endgame_graph.py` (branch-and-bound route search)

```python
@dataclass(frozen=True)
class EndgameRoute:
    legs: tuple[LeveragingLegResult, ...]
    total_dv_kms: float
    vinf_entry_kms: float
    vinf_final_kms: float        # at/below target_vinf_floor_kms
    lower_bound_kms: float       # Γ-quadrature bound for the whole route

def solve_endgame(
    moon_system: str,            # e.g. "Jupiter"
    entry_moon: str,
    vinf_entry_kms: float,
    target_vinf_floor_kms: float,  # the BEND-FEASIBLE CAPTURE target (the goal)
    dv_budget_kms: float,
    ephem: Ephemeris,
) -> EndgameRoute | None: ...
```

`target_vinf_floor_kms` is the **bend-feasible capture target** — the V∞ at which
the moon-tour's flybys become ballistically linkable (the same feasibility the
Phase 6 ballistic sweep tested against, e.g. ~6.0 km/s at the Jovian moons). It is
distinct from A1's V̄∞ leveraging-validity edge. Best-first branch-and-bound.
Node = `(current_moon, current_V∞)`. Expand by candidate leveraging legs + direct
intermoon transfers. **Admissible lower bound** = the Γ-quadrature
`ΔV_min(V∞_current → target_vinf_floor)` (`vilm.py`, Eq. 12-13): no real leg beats
the analytic minimum, so any branch with `cost_so_far + Γ_bound ≥ incumbent` is
pruned soundly. Returns the cheapest route to the target floor, or `None` →
method-versioned EMPTY-region report.

**Golden gate (blocks Component 3):** reproduce the A5 Ganymede→Europa endgame —
total **Δv_TOT ≈ 1.25 km/s**, bounded below by the **1.71** VILM minimum — within
tolerance.

### Component 3 — Discovery integration

`discover_novel_moon` gains an `endgame_genome=True` path: for a moon-tour
topology that closes only at bend-infeasible high V∞, the solver walks V∞ into
bend-feasible reach, yielding a **powered** moon-tour cycler. The Phase 6
`MethodCapability` gains a `"leveraging"` capability tag; since
leveraging ⊐ no-leveraging, `should_sweep` returns `True` for the
`jovian-IEG-vilm-*`, `jovian-perm-vilm-*`, and `saturnian-titan-vilm-*` regions,
which re-open automatically. A survivor → `review_queue.jsonl` (SILVER-ARTIFACT
discipline) → n-body. The existing **ballistic** EMPTY entries are left
untouched; the leveraging run writes its own method-versioned region records
alongside them (no contradiction — a new, more-capable lane).

### Component 4 — Validation harness (binding)

Three independent gates, in order:
1. **Same-model sourced golden** — A1/A2 (leg), A5 (graph). Tight tolerance.
2. **Independent cross-check** — every leg ΔV ≥ Γ-quadrature floor
   (`gamma_floor_ok`); a leg/route beating the analytic minimum is rejected. A
   different code path (closed-form quadrature vs Lambert/propagation).
3. **n-body confirmation** — survivor propagated in the REBOUND moon-system
   harness; only a closing trajectory is trusted. No catalogue writeback until
   n-body passes and a same-model golden exists for the family.

A clean negative (no route to floor within budget) → success → method-versioned
EMPTY-region report.

## Reused unchanged

`vilm.py` (Γ-quadrature, n:m_K± taxonomy `VilmLeg`, V̄∞ thresholds — now both the
B&B lower bound and the cross-check), `core/lambert`, `core/ephemeris`,
`core/satellites` (`PRIMARIES`, `SATELLITES`), the REBOUND n-body harness
(moon-system variant), the Phase 6 sweep + `should_sweep` gate +
`empty_regions`/`review_queue` infrastructure.

## Data flow

topology + high entry-V∞ + target V̄∞ floor (A1)
→ `solve_endgame` (calls `evaluate_leveraging_leg` per edge; Γ-quadrature prunes)
→ route + total ΔV + per-leg breakdown
→ Phase 6 integration re-sweeps Jovian/Saturnian (leveraging capability)
→ survivor → `review_queue.jsonl` → n-body confirm
→ (only then) catalogue writeback as V3-powered.

## Testing (TDD, bottom-up per Approach A)

- **`leveraging_leg`:** resonance geometry (return to M, residual ~0); DSM ≥ Γ
  floor; A1/A2 goldens; infeasible resonance → `converged=False` with recorded
  reason.
- **`endgame_graph`:** route found on a toy 2-moon graph; **B&B pruning soundness**
  (admissible bound never prunes the true optimum — tested against brute force on
  a small graph); A5 golden; no-route → `None`.
- **integration:** a leveraging run emits a `"leveraging"`-tagged region record;
  `should_sweep` flips the three ballistic-EMPTY regions to re-sweepable; survivor
  routes to the review queue, **not** the catalogue.

## Error handling

- Infeasible resonance / no real apse solution → recorded prune reason (mirror
  `discover`'s per-leg reason discipline); never a silent drop.
- No route within `dv_budget` → method-versioned EMPTY-region report.
- A4's two suspect cells never used as goldens (asserted in the test docstring).
- Concurrent-agent git rules + the no-tolerance-loosening firewall carry over
  from the Phase 6 sweep.

## Out of scope (YAGNI / deferred)

- **Real DE440 (phase-full, non-circular) endgame** — a later fidelity rung. This
  spec proves the solver on the circular-coplanar model where the goldens live.
- **Broken-plane / inclined leveraging** — coplanar first.
- **Spec 2 (multi-arc DSM into the corrector)** — the catalogue-lift sub-project,
  its own spec.
- **Heliocentric (Earth-Mars) leveraging cyclers** — the moon-system endgame is
  the proving ground; heliocentric reuse is a later extension.

## References

- `docs/notes/2026-06-05-endgame-tisserand-mining.md` — Campagnola & Russell
  "Endgame Problem" Parts 1–2 transcription (goldens A1–A5).
- `docs/notes/2026-06-08-forge-phase6-discovery-design.md` — Phase 6 pipeline +
  §6b capability-subsumption re-sweep gate.
- `src/cyclerfinder/search/vilm.py`, `dsm_leg.py`, `moon_prune.py`.
- Memory: `validation-ceiling`, `orbit-closure-discipline`,
  `negative-results-registry`.
