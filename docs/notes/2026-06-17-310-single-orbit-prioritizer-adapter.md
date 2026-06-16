# 2026-06-17: #310 Single-orbit prioritizer adapter (Phase 1)

## Scope

Architectural-gap closure for the #282 five-tier accessibility prioritizer
(`src/cyclerfinder/search/five_tier_prioritizer.py`). The gap was surfaced
by the #284 asymmetric-corrector scan but is structural: the five-tier
stack scores either patched-conic Lambert LEGS (Tier 0 only) or
representative orbit PAIRS (all five tiers). A single-orbit discovery
candidate — one CR3BP IC + period + Jacobi, the shape emitted by
asymmetric correctors (#284), 3D corrector spikes (#287), or any sub-
family scoring pass — fits NEITHER input mode.

Phase 1 ships the adapter, the smoke tests, and this disposition. No
catalogue writeback. No novelty claims. The adapter is GLUE; the scoring
physics is unchanged.

## The architectural gap from #284

`five_tier_prioritizer.FiveTierPrioritizer` exposes:

* `score_leg(leg)` — Tier 0 (NN ΔV) on one `PatchedConicLeg` (SI r/v
  endpoints + dt + central-body mu). Tiers 1-5 are documented as
  inapplicable: their inputs are CR3BP representative orbits with
  `state0` + `period` (and Floquet structure for Tier 5), and one
  Lambert leg is not such an anchor.
* `score_pair_full(rep_from, rep_to, …)` — All five tiers on an ordered
  PAIR of representative-shape objects.

A single-orbit candidate is a 6-vector + period + (system, Jacobi). It is
neither a Lambert leg nor a representative pair. The five-tier API has
no entry point that accepts it. That blocks:

* Batch scoring a #284-style asymmetric-corrector scan (one row per
  converged single-orbit candidate; the scan's output went unscored
  because there was no input shape that matched).
* Sub-family member ranking (e.g. would have been wanted by #301).
* Any single-orbit shortlist downstream of the SILVER admission widening
  (#339).

## Two adaptation strategies

The adapter is `src/cyclerfinder/search/single_orbit_prioritizer.py`. It
exposes `score_single_orbit(...)` returning a `SingleOrbitScore` dataclass
with per-tier verdicts + a `notes` audit string. Two strategies run
side-by-side:

### Strategy 1 — surrogate-pair neighbor

Finds the nearest catalogue pair-shape row to the candidate by the tuple
distance `(mu_ratio, jacobi_C, period_nondim)` (L2 norm of the relative
differences). Mu is the system-identity coordinate; Jacobi is the energy-
level coordinate; period is the characteristic-timescale coordinate. The
surrogate is substituted into one slot of a representative pair; the
candidate fills the other. Tiers 0-3 then run on the pair.

The surrogate-pair verdict is HYBRID: it is `candidate → surrogate`, not
`candidate → candidate's eventual partner`. The adapter's docstring is
explicit about that.

Defensible-not-sourced: no source paper publishes a single-orbit-to-pair
surrogate-distance recipe (the question doesn't arise in the paper batch).
The tuple choice is the smallest physically-meaningful descriptor; the
L2-with-relative-differences norm is the standard multi-axis distance
when the axes have different units.

### Strategy 2 — parallel single-orbit pipeline

Tier 4 (Canales-Howell 2023 FTLE corridor strength) and Tier 5 (Hiraiwa-
Bando 2026 lobe-overlap bottleneck flux) do not require a partner orbit.
FTLE samples the candidate's planar position in the field (the scorer's
existing degenerate-same-point geodesic branch reports the local cell
value). The lobe scorer self-graphs the candidate's own manifold
partition. The adapter builds a synthetic `ResonantMember` from the
candidate's `(state, period)` via `resonance_network._planar_floquet` so
the lobe partition can be computed.

Both strategies run in the default policy: surrogate-neighbor lookup
provides tiers 0-3 (when feasible), parallel pipeline always provides
tier 4 (when an FTLEScorer is wired) and attempts tier 5.

## Smoke test verdicts

Test suite: `tests/search/test_single_orbit_prioritizer.py` (4 tests, all
pass in ~10 s on a coarse 8×8 FTLE grid).

Smoke IC: the #287 spike's `planar_baseline_z0eq0` row — the Braik-Ross
C11a-cycler at C=3.1294 (state and period reproduced by the #249
corrector against Braik-Ross 2026 Table 2 to 0.0011%). The state matches
the catalogue row `braik-ross-c11a-cycler-2026`'s
`orbit_elements.cr3bp.state_nd` verbatim — a sourced golden anchor.

| Tier | Verdict on smoke IC | Reason |
|------|--------------------|--------|
| 0 (NN ΔV) | skip | Pure-CR3BP candidate carries no heliocentric SI state; NN returns `fallback_used=no-heliocentric-state`. |
| 1 (Braik-Ross overlap) | skip | No Tier-1/2 scorer attached on the smoke prioritizer (atlas build is expensive; out of smoke budget). |
| 2 (Zhou-Armellin ΔV) | skip | Same — no two-tier prioritizer attached. |
| 3 (Kumar perigee) | skip | No resonance-network scorer attached. |
| 4 (Canales-Howell FTLE) | 0.0 | Candidate position (x=-0.81, y=0.0) at C=3.05 sits in an L4-side region of the smoke 8×8 grid where the chaos class is `escape`; corridor strength = 0 (the smoke grid is too coarse to resolve the (1,1)a-cycler basin, but the scorer ran and returned a finite [0,1] score). |
| 5 (Hiraiwa lobe flux) | 3.13 | Self-graph lobe-overlap on the synthetic ResonantMember; the C11a UNSTABLE orbit has Floquet `|λ|=2.58e4` (Braik-Ross σ=1.0482 TU^-1), so its manifold tubes produce overlapping lobes at the default `r_lobe_threshold=0.002`. |
| combined_rank | 1.0 | Single-candidate rank-product (every tier ranks the lone candidate at rank 1; meaningful composite requires a multi-candidate batch). |
| surrogate | `braik-ross-c11a-cycler-2026` | Tuple-distance 0 (exact match — the spike IC is the verbatim catalogue state). |

The four tests cover:

1. **`test_smoke_287_spike_returns_score`** — full default run; tier 4
   must populate, tier 5 may populate (or skip with a notes audit). PASS.
2. **`test_surrogate_neighbor_matches_c11a`** — surrogate-neighbor lookup
   pins `braik-ross-c11a-cycler-2026` at tuple distance < 1e-4. PASS.
3. **`test_skip_tiers_yields_partial_score`** — `skip_tiers=(0,1,2,3)`
   collapses to parallel pipeline only; surrogate lookup bypassed. PASS.
4. **`test_no_neighbor_fallback_to_parallel_pipeline`** — tiny-mu Sun-
   Earth-like system; no catalogue row within `max_distance=0.1`; tiers
   0-3 stay None; tier 4 still runs on the EM FTLE field. PASS.

## Module signatures

```
@dataclass(frozen=True)
class SingleOrbitScore:
    candidate_id: str
    surrogate_pair_neighbor_id: str | None
    tier_0_score: float | None
    tier_1_score: float | None
    tier_2_score: float | None
    tier_3_score: float | None
    tier_4_score: float | None
    tier_5_score: float | None
    combined_rank: float | None
    notes: str = ""

def score_single_orbit(
    candidate_state0: NDArray[np.float64],
    candidate_period_nondim: float,
    candidate_system: cr3bp.CR3BPSystem,
    *,
    candidate_id: str = "single-orbit-candidate",
    candidate_jacobi: float | None = None,
    prioritizer: ftp.FiveTierPrioritizer | None = None,
    use_surrogate_neighbor: bool = True,
    surrogate_neighbor_max_distance: float | None = None,
    catalogue_path: Path | None = None,
    extra_surrogates: list[_SurrogateCandidate] | None = None,
    skip_tiers: tuple[int, ...] = (),
    notes: str = "",
) -> SingleOrbitScore: ...

def find_surrogate_neighbor(
    candidate_mu: float,
    candidate_jacobi: float,
    candidate_period_nondim: float,
    *,
    max_distance: float | None = None,
    catalogue_path: Path | None = None,
    extra_surrogates: list[_SurrogateCandidate] | None = None,
) -> tuple[_SurrogateCandidate, float] | None: ...
```

## Catalogue read-only seam

`_load_catalogue_surrogates(...)` scans `data/catalogue.yaml` for rows
whose `orbit_elements.cr3bp` block exposes
`(mass_ratio, jacobi_constant, period_nd, state_nd)`. Eight such rows
are presently available (all in the Braik-Ross / Ross-RT Earth-Moon
block):

```
ross-rt-em-cycler-11-2025
ross-rt-em-cycler-21-2025
ross-rt-em-cycler-31-2025
ross-rt-em-cycler-32-2025
ross-rt-em-cycler-33-2025
braik-ross-c11a-cycler-2026
braik-ross-c11b-cycler-2026
braik-ross-c32-cycler-2026
```

Other catalogue rows have `state_nd: null` (not in source; v4.2 backfill
gap) and don't qualify. Adding pair-shape candidates from the negative-
registry is a future seam (`extra_surrogates` parameter).

## Phase 2 recommendation

Re-run #284's asymmetric-corrector scan with the new single-orbit
prioritizer attached. The scan emits one converged-orbit row per cell
to `data/scan_284.jsonl`; piping each row through `score_single_orbit`
(with a real Tier 1-3 scorer set wired to the prioritizer, not just the
smoke Tier 4 + 5 used here) gives the per-row accessibility verdict that
the original #284 scan went without. The surrogate-pair Strategy 1 is
particularly valuable here: every asymmetric scan cell at the Braik-Ross
common Jacobi (C=3.1294) has a sourced catalogue neighbor for the pair
slot.

Open question for Phase 2: the rank-product composite is only meaningful
ACROSS a batch (single-candidate rank-product is trivially 1.0). The
caller, not this adapter, should compose the per-candidate verdicts via
`five_tier_prioritizer.rank_product_score` on whichever tier-score columns
are populated across the batch.

## Files

- `src/cyclerfinder/search/single_orbit_prioritizer.py` (new)
- `tests/search/test_single_orbit_prioritizer.py` (new)
- `docs/notes/2026-06-17-310-single-orbit-prioritizer-adapter.md` (this file)

## Discipline log

- READ-ONLY on `five_tier_prioritizer.py`, `physical_sanity.py`, all
  Phase 1 scorer modules.
- READ-ONLY on `data/catalogue.yaml`.
- No catalogue writeback. No novelty claims.
- Stayed out of #320 scan paths and #340 silver-validation paths.
- `tests/test_catalogue_rediscovery.py` not touched (frozen census).
- `uv run ruff check` + `ruff format --check` both clean.
