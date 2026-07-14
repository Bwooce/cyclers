# #590: bounded follow-up on the 6 #588-flagged unmatched Sun-Earth ER3BP clusters

**Date:** 2026-07-14
**Origin:** `docs/notes/2026-07-14-588-final-adjudication.md`. Bounded, non-novelty-claiming
follow-up on clusters 40/42/43 (primary) and 24/25/39 (secondary) from #588's 20-cluster
unmatched pool. Scope was fixed by the task: (1) a live, independently-verified literature
search against the classical quasi-satellite-orbit (QSO) corpus, (2) a cheap IC-interpolation
connectivity heuristic against the nearest known Gurfil-Kasdin family, (3) a longer-horizon
drift re-check for clusters 40 and 42. **No catalogue rows changed, no novelty claimed** --
this narrows the open question for a future adjudication pass.

## Step 1: targeted quasi-satellite-orbit literature search

Searched WebSearch + independently verified every DOI via direct CrossRef API resolution
(`https://api.crossref.org/works/<doi>`), then read full text (or the free arXiv preprint
where the publisher page was paywalled) rather than trusting title/abstract alone, per
[[feedback_ground_citations_against_content]].

**Verified and added as new `CorpusAnchor` entries** in `src/cyclerfinder/search/literature_check.py`
(`provenance="verified-against-source"` on all three -- full text read, not just CrossRef
metadata):

1. **Mikkola, Innanen, Wiegert, Connors & Brasser 2006**, "Stability limits for the
   quasi-satellite orbit," MNRAS 369(1):15-24, DOI `10.1111/j.1365-2966.2006.10306.x`.
   CrossRef-confirmed (title/authors/journal/volume/pages match exactly). Full text read via
   the free copy at `physics.uwo.ca/~pwiegert/papers/2006MNRAS.369.15.pdf`. Key facts:
   - A secular-averaged perturbation theory in osculating heliocentric elements
     (a, e, I, omega, Omega, M), for the SAME Sun-Earth mass ratio the Gurfil-Kasdin
     `MU_SUN_EARTH_GURFIL_KASDIN` constant uses (mu = 3.003e-6).
   - Real Earth quasi-satellite examples (2003 YN107, 2004 GU9) reach geocentric distances
     of order sigma\*a1 ~ 0.1-0.2 AU (1.5e7-3e7 km) -- the SAME order of magnitude as the
     4-13M km #588 cluster pool.
   - **Headline theoretical result, directly relevant to #588/#590:** permanently stable QS
     motion requires SMALL inclination (I < |e - e1|); for large inclination -- the #588
     3D-DRO regime clusters 40/42/43 sit in -- only TEMPORARY capture occurs, escaping to
     horseshoe/unbound motion on a case-dependent timescale (years to ~1000 yr in their real
     examples).
2. **Sidorenko, Neishtadt, Artemyev & Zelenyi 2014**, "Quasi-satellite orbits in the general
   context of dynamics in the 1:1 mean motion resonance: perturbative treatment," CMDA
   120:131-162, DOI `10.1007/s10569-014-9565-4`. CrossRef-confirmed. Abstract independently
   fetched (full text paywalled): a **spatial (3D)** but circular restricted three-body
   problem "Sun-planet-asteroid," generic to any planet, double-averaged treatment of
   small-eccentricity/small-inclination quasi-satellite motion and its transitions to other
   1:1-resonance regimes. Corroborates Mikkola et al.'s inclination-instability picture from
   an independent analytical route.
3. **Pousse, Robutel & Vienne 2017**, "On the co-orbital motion in the planar restricted
   three-body problem: the quasi-satellite motion revisited," CMDA 128:383-407, DOI
   `10.1007/s10569-016-9749-1` (arXiv:1603.06543). CrossRef-confirmed; full text read via the
   free arXiv preprint. **PLANAR ONLY** -- 3D/inclination effects are explicitly deferred to
   "a forthcoming work" in the paper's own conclusion. Confirms Henon's (1969) "family f" --
   already anchored in the corpus from #588 -- is one continuous one-parameter curve from an
   infinitesimal Earth neighbourhood to Sun-collision (X in [0,1] AU), and that for
   small-mass-ratio (terrestrial-planet-like) systems the heliocentric quasi-satellite domain
   dominates nearly the whole family-f curve. This corroborates (does not newly establish)
   the #588 finding that Gurfil-Kasdin's planar families sit on Henon's curve and that
   extending it is expected, not novel -- but being planar-only, it does NOT itself cover the
   3D/inclined regime clusters 40/42/43 sit in.

**Deliberately NOT added** (a genuine, honest negative, not a gap papered over): Lidov &
Vashkov'yak 1993 ("Theory of perturbations and analysis of the evolution of quasi-satellite
orbits in the restricted three-body problem," Kosmicheskie Issledovaniia 31:75-99) and 1994a/b
(Astronomy Letters 20:188-198 and 20:676-690). Both a direct DOI-shaped CrossRef query and a
bibliographic title/author CrossRef search returned no match for any of the three -- these are
apparently pre-DOI-era Russian-journal papers CrossRef has not indexed (the identical honest
gap the existing Henon-1969 anchor already documents). Their existence/titles/journals are
corroborated only by cross-checking Mikkola et al. 2006's and Pousse-Robutel-Vienne's own
reference lists against each other, which is secondary citation, not source-grounding. Left
un-anchored per [[feedback_ground_citations_against_content]] rather than cited on inherited
authority.

**Verdict on Step 1:** none of these papers provide a tabulated orbit-family catalogue
comparable to Gurfil-Kasdin's own Table 4 (no rmin/rmax entries to check clusters 40/42/43
against directly) -- so **this does not "close" the #588 gap** in the sense of confirming or
refuting any specific cluster IC as a published family member. What it DOES establish: the
inclined 3D quasi-satellite/DRO regime these clusters sit in is a well-studied dynamical
class (not itself novel structure), and its own theoretical literature independently predicts
the long-horizon numerical result found in Step 3 below (see there). This is a genuine,
useful, but non-closing finding.

## Step 2: cheap IC-interpolation connectivity heuristic

**Method (and its limits, stated up front):** linearly interpolated the 6D geocentric IC
vector in 15 steps between each cluster's IC and the nearest Gurfil-Kasdin `TABLE34` family
IC (Family J for 40/42/43; Family I for 24; Family C for 25; Family B for 39), holding
`theta0` similarly interpolated, then ran `classify_bounded_drift` (`n_revs=20`, ~20-year
horizon) at each point and tracked the defining ratio ydot0/x0. **This is a weak heuristic
plausibility probe, not a rigorous differential-correction continuation** -- a real
predictor-corrector chain could very plausibly find a genuinely bounded path through a region
this straight-line probe reports as "escape," and conversely a smooth straight-line path is not
proof of a real family relationship. Script: see the step2/step3 scratch scripts referenced
below (not committed; reproducible from the code cited).

**Primary (Family J anchor, `matched=True`/`practically_stable=True` per #581's own
reproduction, i.e. a validated anchor):**

| cluster | ydot0/x0 range (J -> cluster) | path result |
|---|---|---|
| 40 | -1.987 -> -2.076 | bounded at both ends and nearly everywhere; one non-escaping rough patch (growth_ratio up to ~2.06) around t=0.40-0.60, no true escape anywhere on the path |
| 42 | -1.987 -> -2.034 | bounded at both ends, but the path passes through a genuine **escape** region for 6 consecutive interior points (t=0.27-0.67), plus one more escape blip at t=0.80 |
| 43 | -1.987 -> -2.046 | bounded almost everywhere; only a mild non-escaping dip (growth up to ~1.38) near t=0.07-0.20; the cleanest of the three |

The ydot0/x0 ratio varies monotonically along all three paths -- worth noting but not strong
evidence on its own: a Mobius-transform ratio of two linearly-interpolated quantities is
monotonic wherever the denominator (x0) doesn't change sign, which holds trivially here
regardless of whether a real family connects the endpoints.

**Reading:** cluster 43's connectivity path is the cleanest ("same curve, unsampled interior"
is the more natural reading); cluster 42's is the roughest (passes through actual escape, so
"distinct branch" is at least as plausible as "same curve" from this heuristic alone); cluster
40 is intermediate. None of this is decisive -- it is exactly the kind of ambiguous evidence
the task anticipated.

**Secondary (Families B, C anchors -- both `matched=True`/`practically_stable=True`, valid
anchors):** clusters 25 (-> C) and 39 (-> B) both show the SAME pattern as cluster 42:
bounded at both endpoints, escape region in the middle (7-8 interior points for 25, 8 for 39).
Weak, ambiguous evidence, same caveats as above.

**Secondary (Family I anchor for cluster 24) -- invalid anchor, reported as a finding, not
worked around:** Family I is one of the **3 Gurfil-Kasdin families that do NOT reproduce**
under this repo's own #581 pipeline (`data/found/581_niching_ga/analysis_summary.json`:
`"matched": false`, `rmin_ratio=0.216` vs. Table 4, and even the GA's own best-fitness
individual found near Family I is `"practically_stable": false` -- it escapes to 74.8M km
within 5 years). Family I's own literal Table 3 IC (any theta0 phase tried, 0 to 2*pi)
immediately classifies DIVERGENT/escape under `classify_bounded_drift`. **There is no valid
Family-I anchor to run a connectivity check against in this model** -- this is itself a useful
finding: cluster 24's #588 "band-adjacent to Family I" match was already the weakest of the
three secondary candidates, since Family I itself is not a confirmed-stable structure in this
exact numerical model to begin with. No connectivity conclusion is drawn for cluster 24;
this is reported as an honest gap, not silently skipped.

## Step 3: longer-horizon stability re-check (clusters 40 and 42)

**Important correction to the #588 note's framing, made explicit here so it does not
propagate further:** the #590 task text (following #588's note) describes clusters 40/42's
"growing 5yr r-band drift trend" as the motivation for this step. That is a **misnomer** for
what is actually in the data. The `drift` field in
`data/found/583_widened_search/unmatched_20_for_adjudication.json` is the output of
`classify_bounded_drift` run at `n_revs=N_REVS_DEFAULT=50` (confirmed:
`"n_windows_complete": 50` in the JSON, and `scripts/run_583_widened_bounded_drift_search.py`
calls `classify_bounded_drift(..., n_revs=N_REVS_DEFAULT)`), i.e. an **already-50-year**
horizon, not 5 years. The separate `rmin_km_5yr`/`rmax_km_5yr` fields in the same JSON are
from a DIFFERENT, unrelated quick 1yr/5yr `characterize()` snapshot
(`scripts/run_581_gurfil_reproduction.py`). Re-ran `classify_bounded_drift` at 50 (reproducing
the original numbers exactly, confirming no mismatch), 100, 150 and 300 years:

| cluster | n_revs=50 (original) | 100 | 150 | 300 |
|---|---|---|---|---|
| 40 | bounded, growth=1.046, trend=0.066 | bounded, growth=1.050, trend=0.065 | bounded, growth=1.024, trend=0.019 | **escapes** at year ~213 |
| 42 | bounded, growth=1.076, trend=0.107 | bounded, growth=0.971, trend=0.064 | bounded, growth=1.034, trend=0.029 | **escapes** at year ~220 |

Both clusters remain bounded through 150 years (the growth trend actually SHRINKS relative to
the 50-year figure, not a monotonic escalation) but both do eventually escape somewhere between
150 and 300 years.

**Critical control, run to check whether this is specific to clusters 40/42 or a generic
property of this whole orbit regime:** the SAME classifier applied to Family J's own literal,
paper-printed Table 3 IC (the validated anchor itself, `theta0=0`, `mu`/`e` per
`SUN_EARTH_ER3BP`) escapes far EARLIER -- at year ~29.7, confirmed both via
`classify_bounded_drift` (n_revs=50/100/150/300 all report `n_windows_complete=30`, escape) and
via an independent hand-rolled 1/5/20/30-year `solve_ivp` check that also reproduces the
paper's own Table 4 rmin/rmax to 5-6 significant figures at 1 year (6,892,107 km / 8,511,040 km
vs. published 6,892,060 / 8,510,975), ruling out an IC-conversion bug.

**As a bonus (not originally in the task's step-3 scope, but cheap to check given the above
control result), cluster 43** was also re-run at longer horizons: bounded through 150, 300, 500
and even **1000 years** (growth_ratio 0.97-1.16 throughout, no escape).

**Interpretation:** "bounded" is a horizon-relative verdict for essentially every member of
this inclined 3D DRO/quasi-satellite regime, INCLUDING the published, validated anchor
(Family J's own literal IC escapes at ~30 years -- far sooner than clusters 40 or 42, which
last 213-220 years, and much sooner than cluster 43, which is still bounded at 1000 years).
This is not a numerical artifact unique to the unmatched clusters; it is exactly what Step 1's
literature predicts: Mikkola et al. (2006) establish that permanently stable QS/DRO motion
requires small inclination, and that inclined ("large I") quasi-satellite orbits are
generically only TEMPORARILY captured before escaping, on a case-dependent timescale. The
original #588 framing -- treating clusters 40/42's 50-year drift trend as a red flag relative
to an implicitly-assumed-permanent published family -- was subtly miscalibrated: the anchor
itself is not permanently bounded at this precision either. This reframes, rather than
resolves, the open question: clusters 40/42/43 are not obviously LESS stable than the
published anchor (43 is dramatically MORE stable); the escape-timescale variation across all
four (30 / 213 / 220 / 1000+ years) looks more like the expected spread within a single
metastable-orbit continuum than like a bright line separating "real family member" from "GA
artifact."

## Summary for the next adjudication pass

- **Literature (Step 1):** genuinely relevant QSO prior art now anchored in the corpus
  (Mikkola 2006, Sidorenko 2014, Pousse-Robutel-Vienne 2017), but none of it tabulates a
  direct comparison point for clusters 40/42/43 -- the gap is not closed, only better
  characterized. Lidov & Vashkov'yak honestly left un-anchored (CrossRef does not resolve
  them).
- **Connectivity (Step 2):** ambiguous, as expected from a heuristic straight-line probe.
  Cluster 43's path to Family J is the cleanest; cluster 42's is the roughest (genuine escape
  region). Cluster 24's nominal anchor (Family I) is itself not a validated bounded structure
  in this model, so no connectivity conclusion is possible there at all.
- **Long-horizon stability (Step 3):** decisive, and reframes the question rather than
  resolving it. Clusters 40 and 42 both eventually escape (~213, ~220 yr) -- but so does the
  validated Family J anchor itself, at a MUCH shorter ~30 yr. Cluster 43 remains bounded out to
  at least 1000 years, more robust than the anchor. "Escapes eventually" is the literature-
  predicted generic behaviour of inclined QS/DRO motion, not a discriminator between real
  family members and artifacts.
- **Recommendation for the next pass:** none of the three steps yields a clean verdict on
  "same family vs. distinct branch" for any of the 6 clusters. A genuine differential-
  correction continuation (not this task's straight-line heuristic) from Family J to each of
  40/42/43 would be the natural next escalation if this thread is revisited, but that is new
  capability-building work, out of this task's bounded scope.

## #590 status: CLOSED (bounded follow-up complete; no catalogue rows changed, no novelty
claimed). The 6-cluster open question from #588 is narrowed but not resolved -- left for a
future adjudication pass to weigh this evidence.
