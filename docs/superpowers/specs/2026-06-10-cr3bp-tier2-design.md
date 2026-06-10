# CR3BP Tier-2 — Design Spec

**Date:** 2026-06-10
**Status:** Design approved (brainstorming), pending implementation plan.
**Track:** Genome-capability extension #2 (the new dynamical model). Follows
Spec 1 (#179 VILM endgame) and the corrected Spec 2 (#181 ToF-fix catalogue lift).

## Goal

Build the circular-restricted-three-body (CR3BP) capability the catalogue's
`non-keplerian` / Tier-2 rows need: a rotating-frame propagator + Jacobi constant
+ a periodic-orbit differential corrector, validated against a published
Arenstorf orbit, then used to (a) backfill the citation-only Earth-Moon CR3BP rows
and (b) discover new midsize-moon (Saturnian) cyclers the patched-conic genome
cannot represent.

## Why this, why now

Everything the project has built is patched-conic or n-body-inertial. The Tier-2
rows — Earth-Moon **Arenstorf / Genova-Aldrin 3-petal / Wittal** periodic orbits
and the Saturnian midsize-moon (Mimas/Enceladus/Tethys) members — live in the
CR3BP, where the conserved quantity is the **Jacobi constant**, not patched-conic
V∞. They are currently **citation-only** (`period.years=null`, `vinf=null`, no
`orbit_elements.cr3bp{}` populated) because we have no CR3BP machinery. The
moontour Tier-1 note (`docs/notes/2026-06-08-moontour-tier1-complete.md` §"What
Tier-2 holds") explicitly defers these pending "a CR3BP propagator, Jacobi
constant, and the `orbit_elements.cr3bp{}` backfill." This spec builds exactly
that. Pure engineering — no paper acquisition required.

## Binding constraints (orbit-closure-discipline)

1. **Same-model sourced golden:** the propagator/corrector is validated by
   reproducing a PUBLISHED Arenstorf orbit — the canonical planar CR3BP periodic
   orbit (μ = 0.012277471; initial conditions x₀ ≈ 0.994, ẏ₀ ≈ −2.0015851,
   period T ≈ 17.0652166 in nondimensional CR3BP units), reproduced in standard
   numerical-ODE references (Arenstorf 1963; Hairer-Nørsett-Wanner, *Solving ODEs
   I*). The EXPECTED values cite the publication, never our own computation
   (`golden-tests-sourced-only`). The agent MUST confirm the exact sourced numbers
   from a citable reference before using them as the golden — if it cannot, it
   reports NEEDS_CONTEXT rather than using remembered constants.
   **μ caveat:** 0.012277471 is the *Arenstorf test-problem* mass ratio (the
   published orbit is specific to it) — use it directly for the propagator/corrector
   GOLDEN. It is distinct from the *physical* Earth-Moon μ (~0.01215) that
   `cr3bp_system("Earth","Moon")` derives from `satellites.py`, which is what
   backfilling the REAL Earth-Moon rows uses. Do not conflate the two.
2. **Integrator self-check:** the Jacobi constant is conserved to ≤ 1e-10 over an
   orbit period (a model invariant, not table-gated).
3. **Independent cross-check:** a discovered/backfilled CR3BP orbit is
   re-propagated in the INERTIAL n-body harness (REBOUND, Earth+Moon or
   Saturn+moon on rails) and must stay bounded/periodic — a different code path
   and frame from the rotating-frame CR3BP integrator.
4. **No catalogue writeback** until reviewed. A discovered NEW orbit with no
   sourced anchor routes to the review queue (SILVER), never auto-promoted.
5. **No tolerance loosening.** A family that does not yield a periodic orbit is a
   recorded negative.

## Architecture — three phases on a shared foundation

### Phase 1 — CR3BP dynamics core (`src/cyclerfinder/core/cr3bp.py`)

```python
@dataclass(frozen=True)
class CR3BPSystem:
    mu: float                 # mass ratio m2/(m1+m2)
    primary: str              # "Earth" | "Saturn" | ...
    secondary: str            # "Moon" | "Enceladus" | ...
    l_km: float               # characteristic length (secondary SMA) for re-dimensioning
    t_s: float                # characteristic time

def cr3bp_eom(t, state6, mu) -> np.ndarray        # rotating-frame EOM (pos,vel)
def cr3bp_stm_eom(t, state42, mu) -> np.ndarray   # state + 6x6 STM variationals
def jacobi_constant(state6, mu) -> float          # C = 2*Omega - v^2
def propagate(system, state6, t, *, with_stm=False, events=None) -> CR3BPArc
def cr3bp_system(primary, secondary) -> CR3BPSystem  # μ + scales from satellites.py / PRIMARIES
```
Integrated with scipy `solve_ivp` (high-order, tight rtol/atol, dense output,
event detection for plane crossings). μ and the length/time scales come from the
existing `core/satellites.py` (`SATELLITES`, `PRIMARIES`) so Earth-Moon and the
Saturnian pairs are both expressible. *Golden: Jacobi conserved ≤ 1e-10 over a
period for the Arenstorf orbit.*

### Phase 2 — periodic-orbit differential corrector (`src/cyclerfinder/search/cr3bp_periodic.py`)

```python
@dataclass(frozen=True)
class PeriodicOrbit:
    state0: np.ndarray        # corrected periodic initial state (6,)
    period: float             # nondimensional
    jacobi: float
    converged: bool
    closure_residual: float   # |X(T) - X(0)|

def correct_periodic(system, state0_guess, period_guess, *,
                     symmetry="xz", tol=1e-10, max_iter=30) -> PeriodicOrbit
```
STM-based single-shooting: enforce the periodicity/symmetry constraint (e.g. the
classic perpendicular-x-axis-crossing condition for symmetric planar orbits) by
Newton steps using the monodromy/STM. *Golden (blocks Phase 3): from the published
Arenstorf guess it converges to the published `state0`/`period`/Jacobi to tight
tolerance; closure residual ≤ tol.*

### Phase 3a — Earth-Moon backfill (`scripts/cr3bp_backfill.py`)

For the 3 citation-only Earth-Moon rows, populate `orbit_elements.cr3bp{}`
(μ, state0, period, Jacobi) by correcting from each row's sourced guess. Each row
earns a level only where its OWN sourced data anchors it: Arenstorf has the
published IC golden (→ V1 reproduction); Genova-Aldrin 2015 / Wittal 2022 depend
on what those papers published (if no sourced IC/Jacobi, the orbit is reproduced
but stays a documented backfill, not a sourced-validated promotion — flag
honestly). Proposed `_LEVEL_EVIDENCE` + cr3bp{} written to a runlog/note; NO
catalogue writeback.

### Phase 3b — midsize-moon discovery (`scripts/cr3bp_moontour_run.py`)

Seed periodic-orbit families in the Saturnian midsize system (Mimas/Enceladus/
Tethys CR3BP pairs), continue along each family (natural-parameter continuation on
Jacobi), run the gauntlet. A NEW periodic orbit has no sourced anchor → routes to
`review_queue.jsonl` (SILVER) + the inertial-n-body cross-check, never
auto-promoted. EMPTY (no periodic family found in a region) → a method-versioned
record, mirroring Forge Phase 6.

### Phase 4 — validation harness (cross-cutting)

Arenstorf golden (Phase 1+2), Jacobi conservation (Phase 1), inertial-n-body
re-propagation cross-check (Phase 3), review-queue routing for discoveries. No
writeback until reviewed.

## Data flow

`cr3bp_system(primary, secondary)` → μ, scales →
`correct_periodic(guess)` → `PeriodicOrbit` (state0, period, Jacobi) →
[backfill: cr3bp{} for known rows] / [discovery: family continuation → gauntlet] →
inertial-n-body cross-check → runlog + proposed evidence → [review] → writeback.

## Out of scope (YAGNI / deferred)

- **Bicircular / full-ephemeris 4-body** refinement of CR3BP orbits — Tier-3.
- **The T > 3 ballistic-transfer region** (Endgame Part 2) beyond what a periodic
  family naturally covers.
- Low-thrust / broken-plane genomes — the separate #3 sub-project.
- Auto-promotion of discovered orbits — always review-gated.

## Testing (TDD)

- **`cr3bp.py`:** Jacobi conserved over a period (Arenstorf); EOM matches a known
  derivative at a sample state; `cr3bp_system("Earth","Moon").mu` ≈ 0.01215 (from
  satellites.py masses) within tolerance of the canonical Earth-Moon μ.
- **`cr3bp_periodic.py`:** from the published Arenstorf guess, `correct_periodic`
  returns the published state0/period/Jacobi (sourced golden) with residual ≤ tol;
  a non-periodic guess → `converged=False` (honest negative).
- **backfill:** the Arenstorf row's cr3bp{} is populated and its Jacobi matches the
  golden; a row with no sourced IC is flagged backfill-only, not promoted.
- **discovery:** a seeded Saturnian family yields a `PeriodicOrbit` OR a clean
  EMPTY; a discovery routes to the review queue, never the catalogue.
- **cross-check:** a corrected orbit re-propagated inertially (REBOUND) stays
  bounded over several periods.

## References

- `docs/notes/2026-06-08-moontour-tier1-complete.md` — Tier-2 scope + what's deferred.
- `src/cyclerfinder/core/satellites.py` — `SATELLITES`/`PRIMARIES` (μ, SMA scales).
- `src/cyclerfinder/data/catalog.py` — `orbit_elements.cr3bp{}` schema (#88).
- Arenstorf (1963); Hairer-Nørsett-Wanner *Solving ODEs I* (the canonical IC).
- Catalogue rows: `arenstorf-em-figure8-1963`,
  `genova-aldrin-2015-em-3petal-cycler`, `wittal-2022-em-cycler-family`.
- Memory: `orbit-closure-discipline`, `golden-tests-sourced-only`, `gmat-install`.
