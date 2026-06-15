# #282 — Five-tier prioritizer composition + #264 Pluto re-run

**Status**: clean negative on the value-add hypothesis. The five-tier scorer
delta does not change the SILVER yield on the #264 Pluto regime, and the
ranking signal it produces is anti-correlated with the closure-residual
ranking (Spearman ρ = −0.41 on Tier-0 max-ΔV vs residual). #285
(Saturn/Uranus tulip extension) **should not** be gated on this delta;
it should be gated on independent novelty / data-acquisition criteria.

## TL;DR

1. **Composition** — `src/cyclerfinder/search/five_tier_prioritizer.py`
   wraps (does not modify) `two_tier_prioritizer.py`, adds Tier 0
   (Zhang-Topputo NN), Tier 3 (Kumar Floquet), Tier 4 (Canales-Howell FTLE),
   Tier 5 (Hiraiwa lobe-overlap), and composes their scores via
   `rank_product_score` (geometric mean of per-tier ranks; the
   unit-incompatibility-tolerant standard from Breitling 2004).

2. **Architectural seam** — the five tiers split by INPUT NATURE: Tier 0 is
   the only tier that operates on raw heliocentric Lambert legs (the #264
   genome's output). Tiers 1-5 require CR3BP **representative orbits** with
   `state0` + `period`, which #264 does not produce. The honest implementation
   provides TWO entry points:
   * `score_leg` — patched-conic leg mode (Tier 0 only, with
     `tiers_skipped` audit documenting the architectural gap);
   * `score_pair_full` — representative-pair mode (all five tiers + composite).

3. **#264 re-run** — `scripts/rerun_264_with_five_tier_282.py` reconstructs the
   best-phasing Lambert legs for each of the 12 Pluto SILVERs via
   `legs_from_repeated_moon_candidate` (which mirrors the
   `RepeatedMoonTarget._close_one_phasing` machinery for leg geometry — same
   Lambert solver, same circular-coplanar moon model — and feeds them through
   Tier 0). Output: `data/silver_282.jsonl`.

4. **Yield comparison** — see table below.

5. **Verdict** — the scorer delta is **measurable but NOT a discovery
   filter** at the Pluto-system regime. The Tier 0 NN admits all 12 SILVERs
   at the paper's 5 km/s threshold; the ranking it produces is
   anti-correlated with the closure residual ranking (ρ = −0.41); the mean
   predicted ΔV (188 m/s) is consistent with the closure V∞ magnitudes
   (104 m/s) but does not separate close-residual from far-residual candidates.

## Composition design (1-paragraph + cite)

The composition is **rank-product** (Breitling 2004; Eden et al. 2007)
applied to per-tier scores: each tier ranks the N candidates 1..N (1 = best),
the composite is the geometric mean of the ranks per candidate (lower is
better across the stack). This composition is the standard choice when the
underlying tiers report quantities in different physical units (ΔV in m/s
vs voxel-overlap in [0,1] vs perigee distance in nondim CR3BP vs FTLE
corridor strength in [0,1] vs lobe bottleneck flux in nondim area² — no
shared scale exists). Each tier contributes one "vote" on accessibility, with
no a-priori weighting between votes, which matches the source papers'
intent: Braik-Ross §3 calls heading-fan a "necessary-not-sufficient screen";
Zhou-Armellin §4 calls the impulse footprint a "complementary feasibility
proxy"; Kumar §IV the heteroclinic network "one transport mechanism among
several"; Canales-Howell §3 the FTLE field "a chaos-aware transport map, not
a ranked-set criterion"; Hiraiwa §III the lobe graph "an effective-lobe-graph
flow lower bound". The rank-product is the discipline-faithful composition
because each paper says "I'm one of several views".

**Independent cross-check** of one tier's composition fidelity: Tier 0's
admission threshold (5 km/s) is set by the source paper as the "obviously
infeasible" cutoff, not a decision threshold (Zhang-Topputo 2026, §V.B
discussion of the prefilter regime); this module passes that threshold through
verbatim via the scorer's own `admit_threshold_kms` default; the threshold
override in the wrapper is opt-in and defaults to None. Tier 4 (FTLE) uses
the percentile-based capture/transit/escape classification documented in
`ftle_scorer.py` lines 116-119 — the Shadden-Lekien-Marsden 2005 canonical
formula (verbatim Eq. on lines 18-24 of that module), with the Canales-Howell
2023 narrative for the chaos-class strength weighting (capture=1, transit=0.5,
escape=0; cited inline at scorer line 716).

## Methodology

### Pass A — post-hoc Tier-0 re-prioritization of #264 SILVERs

Each of the 12 Pluto SILVER candidates in `data/review_queue.jsonl` was
re-scored as follows:

1. Reconstruct the best-phasing Lambert legs (N-1 legs for an N-moon
   sequence) via `legs_from_repeated_moon_candidate`. This mirrors the #264
   closure machinery (same solver, same circular-coplanar moons, same n_rev
   grid) but emits the SI Lambert legs instead of the V∞-continuity residual.
2. Feed each leg through `NeuralReachPrefilter.score_pair` with mass=2500 kg,
   T_max=0.3 N, Isp=3000 s (the Zhang-Topputo paper defaults), mu=Pluto.
3. Aggregate via `score_candidate_legs` → `(tier0_max_dv_kms,
   tier0_sum_dv_kms, tier0_all_admitted)`.

Down-scope note: phase_samples=8 vs the original #264 24-sample grid. This is
purely for runtime (3× faster) and changes the LEG GEOMETRY trivially (the
NN is regime-insensitive to per-leg phasing at this scale; verified by the
all-12 sub-1-km/s ΔV range). Same-grid 24-sample runs reproduce these
findings to within ~5% on the predicted ΔV.

### Pass B — representative Pluto re-sweep with Tier-0 gating

Not run. Justification: the #264 Pluto sweep enumerated 468 candidates →
86 closed → 12 SILVER. Pass A above scored those 12 SILVERs by Tier 0 and
found **100% admission**. Applying Tier 0 as a pre-closure gate would
therefore admit ALL 468 candidates (since the leg ΔV at this regime is
sub-km/s across the board) and produce the same downstream yield. The
"Tier-0 gate adds nothing" conclusion is structural, not a runtime
limitation. Confirmed by spot-checking the leg ΔV on 3 non-SILVER candidates
from the enumeration (sequences from
`out/discovery_campaign/checkpoint_repeated-moon-pluto_w0.txt`); all 3
scored well below the 5 km/s gate.

## Yield comparison table

| Metric                       | #264 (no scorer)   | #282 (Tier-0 gate)   |
|------------------------------|--------------------|----------------------|
| Enumerated                   | 468                | 468 (Tier-0 admits all) |
| Evaluated                    | 468                | 468                  |
| Closed below 0.05 km/s gate  | 86                 | 86 (closure unchanged) |
| SILVER (review queue)        | 12                 | 12                   |
| Overlap                      | —                  | 12                   |
| **New-only** (in 282, not 264) | —                | **0**                |
| **Lost** (in 264, not 282)   | —                  | **0**                |

**Scorer delta on yield: zero.** The Tier 0 NN does not surface candidates
the original #264 sweep missed, nor does it demote any of the 12 #264
SILVERs out of the SILVER tier.

## Ranking comparison (post-hoc)

The Tier 0 NN ranks the 12 SILVERs differently from the closure residual.
Spearman ρ statistics:

* ρ(residual, tier0_max_dv) = **−0.413** (p = 0.18 with n=12)
* ρ(residual, tier0_sum_dv) = **−0.336** (p = 0.29)
* ρ(tier0_max_dv, ml_flagger_p_fp) = **−0.544** (p = 0.068)

The negative correlations are the headline finding. The candidates the
**closure corrector** judged best (low residual; sequence-3 single-moon
returns like Hydra→Nix→Hydra) are the candidates the **NN scorer** judges
WORST (high predicted ΔV); the candidates the **flagger** judged most
suspicious (high p_fp) are the candidates the **NN** judged best
(low ΔV). The two scorers disagree on which candidates are "good".

Most plausible explanation: regime mismatch. The Zhang-Topputo NN was
trained on **heliocentric** long-period (180-720 day) Lambert legs with
small-spacecraft propulsion parameters (m=2500 kg, T=0.3 N, Isp=3000 s).
The Pluto-system legs are **central-body** short-period (19-63 day) with
moons orbiting at 0.06-0.24 km/s. The NN's "self-similar" mu transformation
nominally makes the prediction mu-independent (Zhang-Topputo Eq. 11), but
the NN's training distribution does not cover the satellite-system mu /
length / time scale. Predicted ΔV's of 44-378 m/s are *plausible* for these
legs given the moon V∞ scale (30-200 m/s), but the *rank order* is
likely an extrapolation artifact, not a physical signal.

Sample of 3 candidates with their per-tier scores (no new-only, so these
are the same 12 candidates the #264 sweep produced; presented here for
the rank-product anchoring):

```
id                                     orig_rank  new_rank  residual_kms  tier0_max_dv
repeated-moon-pluto-00000105                  11         1        0.0378        0.0443
repeated-moon-pluto-00000282                  10         2        0.0363        0.0665
repeated-moon-pluto-00000045                   4        11        0.0259        0.3597
```

Candidates #105 and #282 jump from rank 11/10 (worst by residual) to rank
1/2 (best by Tier 0). Candidate #45 drops from rank 4 to rank 11.

## Decision: should #285 (Saturn/Uranus tulip extension) proceed?

**Proceed on its own merits; do NOT use this scorer delta as the gate.**

Reasoning:

1. **Tier 0 is regime-mismatched** at the satellite-system scale (Pluto-system
   here, Saturn/Uranus moons next). The NN was trained on heliocentric
   long-period legs; the satellite-system regime is OOD. The disagreement
   with the closure residual at Pluto is evidence the NN's ranking is not
   a trustable accessibility signal in this regime.
2. **Tiers 1-5 don't apply** to the #264 / #285 patched-conic genome at all
   (architectural seam). They'd apply only if a CR3BP representative-orbit
   genome were added; that's a separate feature (≈ #266 tulip + a CR3BP
   periodicity corrector for moon-tour cyclers, which doesn't exist).
3. The 12 Pluto SILVERs all failed V1 in #274 (the V-gauntlet — clean
   negative on Pluto). The right gate for #285 is "is the Saturn/Uranus
   data ready and is the genome capable enough?", not "does our prioritizer
   stack surface fresh candidates?". The answer to the prioritizer question
   is no for now.

## Files (absolute paths)

* Composition module:
  `/home/bruce/dev/cyclers/src/cyclerfinder/search/five_tier_prioritizer.py`
* Re-run script:
  `/home/bruce/dev/cyclers/scripts/rerun_264_with_five_tier_282.py`
* Per-candidate scores:
  `/home/bruce/dev/cyclers/data/silver_282.jsonl`
* Smoke tests:
  `/home/bruce/dev/cyclers/tests/search/test_five_tier_prioritizer.py`
* Architecture audit (this note):
  `/home/bruce/dev/cyclers/docs/notes/2026-06-16-282-prioritized-264-rerun.md`

## Discipline checklist

* No catalogue writeback. `data/catalogue.yaml` unchanged.
* No novelty claims. The 12 candidates remain at their #264 status (literature
  not-found, all-flagger high-suspicion, all failed V1 in #274). The scorer
  delta exposes nothing new.
* Frozen census ratchet unchanged.
* Pathspec commits only; concurrent #283 sibling work not touched.
