# #583 (stage 3b of #581): corpus anchors, drift classifier, widened-domain build + positive control

**Date:** 2026-07-13
**Scope of this note:** the corpus-anchor prerequisite + build + positive-control
deliverable only, per #583's own gate ("NOT the full novel-territory sweep,
which is a longer job the coordinator will launch and own separately"). No
novelty claim is made anywhere in this note.

## Part 1: corpus-anchor prerequisite (mandatory, done first)

Filed 3 new anchors in `src/cyclerfinder/search/literature_check.py`'s
`KNOWN_CORPUS`, all DOIs independently verified via direct CrossRef API
resolution 2026-07-13 (title/authors/journal/volume/pages cross-checked
against the CrossRef record, not inherited from a prior note):

1. **Gurfil & Kasdin 2002** (`10.1016/S0045-7825(02)00481-4`) — the 14
   published families' own r-band footprint. Without this anchor, a widened
   search would rediscover the paper's own results and wrongly clear them
   "novel" (the #577 Io-Callisto false-novel-factory trap).
2. **Sun-Earth co-orbital dynamics** — Wiegert, Innanen & Mikkola 1997
   (*Nature* 387:685-686, `10.1038/42662`, discovery of 3753 Cruithne's
   horseshoe orbit) and Namouni 1999 (*Icarus* 137(2):293-314,
   `10.1006/icar.1998.6032`).
3. **Henon family-f distant retrograde orbits** — citation gap flagged
   honestly rather than fabricated: Henon's original family-f papers
   (1969-1970, *Astronomy & Astrophysics*) predate CrossRef's DOI coverage
   for that era and could not be independently verified via the CrossRef API
   in the time available for this dispatch. The anchor is filed with the
   commonly-cited title/venue but WITHOUT an unverified DOI — a TODO for
   whoever next touches this corpus area to track down and verify a proper
   DOI or archival reference, not a blocking gap for #583 itself (DROs are
   still correctly excluded from novelty claims via the anchor's keyword/
   topology matching even without a DOI).

## Part 2: drift classifier

`src/cyclerfinder/data/validation/er3bp_drift_classifier.py` — bounded-vs-
divergent classifier for quasi-periodic candidates (which have no period to
gate on, unlike stage 2's family reproduction). Propagates 50 years past the
1-year fitness window, classifies bounded (stationary geocentric r-band) vs.
divergent (secular rmax growth or escape) via windowed quartile-ratio +
trend tests. Own positive control (per #583's own mandate): confirmed
stage 2's known-good families classify bounded, and a deliberately-escaping
test IC classifies divergent, before trusting it on anything new — see
`tests/data/test_er3bp_drift_classifier.py`. Also includes
`spot_check_theta0_robustness()`, a cheap diagnostic (not a gate) re-testing
a bounded verdict at 2 other launch phases, per
[[project_388_wall_energy_selective]]'s epoch-fragility lesson.

## Part 3: widened-domain search + positive control

`scripts/run_583_widened_bounded_drift_search.py` — 7 partitions (`P1`-`P7`)
each targeting a subset of Gurfil-Kasdin's 14 families by freeing a
different subset of the 7-gene genome (6 interleaved state slots + free
`theta0`). Reuses `search/niching_ga.py::run_deterministic_crowding` and
`core/er3bp_geocentric.py::gurfil_kasdin_fitness` UNMODIFIED, per #583's own
scope.

### First attempt: FAILED (0/6), root-caused and fixed

The first P1 run (population=200, 400 generations) converged the entire
population to fitness mean=0.999188/max=1.000000/std=0.003418 -- but
`--analyze` showed **0/6 families matched**, and every one of the 6
per-family comparisons showed an IDENTICAL candidate block. An independent
Fable review (dispatched by the coordinator) diagnosed this precisely: the
originally-chosen `LO_WIDE = 0.002` AU (~299,200 km) sits well inside Earth's
Hill sphere and re-admits a physically trivial, strictly fitness-dominant
quasi-circular basin -- `gurfil_kasdin_fitness` (Eq. 15) rewards only
1-year annulus thinness with no periodicity/family content, and deep inside
the Hill sphere essentially ANY near-circular orbit scores ~1 (the trivial
basin's deficit was ~5.6e-9 vs. ~8.7e-7 for the best genuine target family --
a 150x fitness gap). Deterministic crowding cannot protect a lower-fitness
niche against a larger, strictly-higher-fitness basin under
`child >= parent` replacement, so by generation ~62 the entire population
collapsed onto the trivial solution. Gurfil-Kasdin's own `LO_R` (1e6 km) is
what excludes this attractor in the original paper -- the widened choice of
`0.002` AU (chosen to explore "between and around" the 12 published boxes)
was an unvetted overreach into territory the source paper deliberately
excluded, not a valid widening.

The same review caught a second, independent bug: `theta0` (free over
`[0, 2*pi)` in this widened search, vs. stage 2's fixed `{0, pi}`) was
compared via a raw linear normalized difference in
`match_family_in_widened_population`, so a genuine candidate that landed
near `theta0 = 2*pi - epsilon` against a published `theta0 = 0` target was
spuriously rejected as far away.

### Fixes applied

1. `LO_WIDE` reverted to `LO_R` (imported from `run_581_gurfil_reproduction.py`
   unmodified) -- removes the trivial deep-Hill basin. `HI_WIDE` (widened
   above the paper's `HI_R`), full-signed `V_WIDE` per free velocity
   component, and free `theta0` all remain as the genuine widening.
2. `match_family_in_widened_population`'s theta0 comparison now wraps the
   normalized difference into `[-0.5, 0.5]` (circular distance) before
   computing the RMS distance, instead of a raw linear difference.

### Corrected result (P1 partition, same seed, full 400 generations)

```
family A: MISS  ic_dist=0.1330(FAR) type=DRO vs DRO rmin_ratio=0.168 rmax_ratio=0.094
family B: MISS  ic_dist=0.0330(ok)  type=DRO vs DRO rmin_ratio=0.443 rmax_ratio=0.337
family C: MATCH ic_dist=0.0005(ok)  type=DRO vs DRO rmin_ratio=0.980 rmax_ratio=0.977
family D: MISS  ic_dist=0.0230(ok)  type=DRO vs DPO rmin_ratio=3.051 rmax_ratio=1.099
family E: MISS  ic_dist=0.0232(ok)  type=DRO vs DPO rmin_ratio=3.184 rmax_ratio=1.135
family F: MISS  ic_dist=0.0034(ok)  type=DRO vs DRO rmin_ratio=4.416 rmax_ratio=1.000
drift classification: 200/200 high-fitness members bounded at 50yr
reproduction: 1/6
```

Population converged to mean=0.998937/max=0.999999/std=0.003960 by gen 400 --
population index 0,2,3,4,5 all sit within 0.0003 (normalized) of the family-C
target with fitness 0.999999; verified directly from the checkpoint.

**Verdict (independent Fable adjudication, second pass): GO.** A positive
control validates machinery, not completeness. Family C is an excellent
end-to-end demonstration (ic_dist 5e-4, both feature ratios ~0.98) that, once
the trivial basin was removed and the theta0 wrap fixed, the pipeline finds a
genuine bounded family and the matcher recognizes it under the SAME
pre-registered criterion stage 2 used. This passes #583's own gate text
("recovers its neighboring known family/families") -- reported here as
"machinery validated on one family," not a 6/6 reproduction claim.

**Why only 1/6, and is that itself a further bug?** No -- ruled out a
systematic `characterize()` measurement bug (the same function scored C
correctly and was already validated per-family in stage 2's own 12-set
reproduction). The B/F near-miss signature (good IC distance + one matching
feature ratio, one badly off) is exactly what an argmin-by-normalized-IC
matcher produces when the population is NOT actually sitting on those
families: a small *normalized* RMS distance in a widened box is still a
large *physical* offset, and DRO rmin is the more sensitive feature: F's
`rmax_ratio=1.000` is a near-circular orbit that coincidentally shares F's
outer radius -- single-feature agreement, correctly rejected by the
combined criterion (the "it matched!" trap working as designed, not
failing). The deeper cause, per the second Fable pass: `gurfil_kasdin_fitness`
is a *boundedness* measure that saturates near 1.0 across the whole bounded
continuum in this box -- nothing in the objective itself steers the search
toward specific published families once the trivial basin is removed, so
the population settles wherever the boundedness landscape happens to peak
first (family C here), not necessarily distributed across all 6 targets.
Population-size/niche-capacity limits (200 covering 6 families in one box,
vs. stage 2's 200 per 1-3 families in a narrower box) compound this but are
not the sole cause. **Consequence for any future full sweep:** per-family or
radial partitioning (closer to stage 2's own per-set structure) is needed,
not one wide box per family group -- this is scoped OUT of the current
dispatch, which only required proving the corrected machinery recovers at
least one known family.

**theta0 spot-check caveat:** one of 5 spot-checked bounded members
("member 1" in the driver's diagnostic loop, a population-index-order pick,
not necessarily related to any specific family match) flipped
bounded/divergent classification across launch phases. Directly verified
this is population index 1 specifically: fitness 0.982 (vs. 0.999999 for
its neighbors) and normalized distance 0.262 from the family-C target --
i.e. a distinct, weaker, more marginal individual, NOT the family-C match
itself (which sits at indices 0/2/3/4/5, all within 0.0003 of the C target
and confirmed stable). This is a documented caveat about classification-
boundary fragility for marginal candidates, not a defect in the validated
C-match or a blocker for this gate.

## Reproduction

```
uv run python scripts/run_583_widened_bounded_drift_search.py --partition P1 --workers 8
uv run python scripts/run_583_widened_bounded_drift_search.py --analyze --partition P1
uv run pytest tests/data/test_er3bp_drift_classifier.py -v
```

## Not yet dispatched

The full widened-domain novelty sweep across all 7 partitions, with the
family-targeted/radial partitioning redesign this positive control shows is
needed for genuine multi-family recovery -- a longer job for a future,
separately-scoped dispatch, per #583's own gate.
