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

## Reproduction (P1, historical -- superseded by Part 4 below)

```
uv run python scripts/run_583_widened_bounded_drift_search.py --partition P1 --workers 8
uv run python scripts/run_583_widened_bounded_drift_search.py --analyze --partition P1
uv run pytest tests/data/test_er3bp_drift_classifier.py -v
```

## Part 4: partition redesign + smoke-scale validation (2026-07-13, follow-on dispatch)

Follow-on task: replace the 7 wide `P1`-`P7` partitions with a genuinely
narrower, radially-partitioned scheme, per Fable's own recommendation
("partition radially into ~stage-2-sized sub-boxes... pop 200 each"), and
validate cheaply (one smoke-scale run) that the redesign out-recovers old
`P1`'s 1/6. This section reports 3 rounds of iteration -- 2 of which caught
real bugs, and a 3rd that surfaced a genuine, deeper limitation of the
search machinery itself (reported honestly, not spun as a clean win).

### Design: 3 radial bands

1. **`deep_hill`** `[LO_DEEP=0.002 AU, LO_R=1e6 km]`: the originally-mistaken
   floor's own territory (the trivial deep-Hill-sphere basin diagnosed in
   Part 3 above). Its own partition, `DEEP_HILL`, judged ONLY by
   `classify_bounded_drift` + the corpus anchors (Henon family-f, Sun-Earth
   co-orbital) -- never Table 3/4 match (nothing published lives here).
2. **`paper`** `[LO_R, HI_R=1e7 km]`: where all 14 published families live.
   Split into 14 SINGLE-family partitions (see "Why single-family, not
   stage 2's own A/B/C grouping" below).
3. **`beyond_hi_r`** `[HI_R, HI_OUTER=0.15 AU]`: genuinely uncharted
   territory past every published family. Its own partition,
   `BEYOND_HI_R`, judged the same way as `deep_hill`.

### Round 1: position-floor correctness bug (caught before any run)

A first draft gave every `paper`-band partition the SAME blanket
`[LO_R, HI_R]` position bound per free component (inherited from the
original #583 build's own decision). Checked directly against every
`TABLE34` IC's own per-COMPONENT values (not just Euclidean radial
magnitude) before trusting this:

| family | free position components | within `[LO_R, HI_R]`? |
|---|---|---|
| A-F | x only | yes (F is 8e-9 below LO_R -- floating-point boundary, not real) |
| G | x=y=0.004727 | **NO** (29% below LO_R) |
| H | x=0.005013 | **NO** (25% below LO_R) |
| I | x=0.01490, y=0.01513 | yes |
| J | x=0.03348, y=0.00775, z=0.03676 | yes |
| K | x=0.00584, y=0.0000547, z=-0.000909 | **NO** (y, z far below/negative) |
| L | x=y=0.00386, z=0.00950 | **NO** (x, y 42% below) |
| M | x=y=z=0.00386 | **NO** (42% below) |
| N | x=0.00766, z=0.006684 | yes |

Root cause: stage 2's own Table 2 divided BOTH `lo` and `hi` by `sqrt(2)`
(G/H/I/N) or `sqrt(3)` (J/K/L/M) -- a shared "vector norm budget" across
however many genes (position AND velocity) that set left simultaneously
free. A blanket `[LO_R, HI_R]` per component EXCLUDES several families' own
published state entirely -- a hard correctness bug (the family could never
be found, not just harder to find), not just a niche-capacity dilution.
This had never surfaced before because the original #583 build only ever
ran `P1` (A-F) to completion; `P2`-`P7` (G/H/I/J/K/L/M/N) were never
actually executed. Fixed: `_POS_S2 = (LO_R/sqrt(2), HI_R/sqrt(2))` for
G/H/I/N, `_POS_S3 = (LO_R/sqrt(3), HI_R/sqrt(3))` for J/K/L/M (further
individually bracketed for J/K/L, which would otherwise collide -- see
Round 2).

### Round 2: byte-identical-bounds bug (caught by the first smoke test)

First smoke-test attempt split A-F into `ABC` and `DEF`, judged separately,
but BOTH partitions were given the identical bounds `_sig_x_vy(_POS_S1)`
(same position range, same free dims, same full-signed velocity) --
splitting the JUDGING target without narrowing the actual GA search box at
all. Ran `ABC` (full 400 gen, pop 200): **1/3** (only C matched; A/B missed,
with WORSE `ic_dist` than old `P1`'s own A/B misses) -- no improvement over
`P1`'s 1/6. The `J`/`K`/`L` partitions (all sharing `_POS_S3` with the same
all-6-state-free signature) had the identical flaw. Fixed: velocity-SIGN
split for A-F (mirroring stage 2's own DRO `vy<0` vs DPO `vy>0` branches:
`ABCF` vs `DE`), and individually-bracketed position ranges for J
(`[0.006, HI_R/sqrt(3)]`), K (`[-0.002, 0.015]`), L (`[0.003, 0.018]`).

### Round 3: multi-family pooling still fails after both fixes (deeper limit, not a bug)

Re-ran the corrected `ABCF` partition (A, B, C, F together, disjoint from
`DE`, velocity-sign-restricted, correct position range, full 400 gen, pop
200): population fitness **std hit EXACTLY 0.0 by generation ~120**
(complete collapse) and **reproduction stayed 1/4** (only C):

```
family A: MISS  ic_dist=0.2725(FAR) type=DEO vs DRO rmin_ratio=1.389 rmax_ratio=1.677
family B: MISS  ic_dist=0.0785(ok)  type=DRO vs DRO rmin_ratio=0.443 rmax_ratio=0.337
family C: MATCH ic_dist=0.0011(ok)  type=DRO vs DRO rmin_ratio=0.983 rmax_ratio=0.978
family F: MISS  ic_dist=0.0905(ok)  type=DRO vs DRO rmin_ratio=4.391 rmax_ratio=1.003
```

Both bug fixes (position floor, byte-identical bounds) are real,
independently necessary corrections -- but neither restored multi-family
recovery. A follow-up diagnostic isolated why, using throwaway standalone
scripts (not committed) that called `run_deterministic_crowding` directly:

* **Fresh random seeds don't help.** 3 alternative seeds (`700001-700003`)
  targeting family A solo (same `ABCF`-structure bounds) mostly converged
  to essentially C's OWN IC (`x=0.006685, vy=-0.029050` vs C's published
  `x=0.00680, vy=-0.02902`), not A's (`x=0.038944, vy=-0.077760`); the 3rd
  seed landed at a different, still-not-A, intermediate point on the same
  low-radius continuum.
* **A tightly-bracketed single-family box doesn't help either.** Restricting
  A's OWN box to `x in [0.03, 0.05]` (well clear of C's neighborhood, 3
  seeds) converged to `x~0.030` -- the box's OWN edge, not A's actual
  `x=0.038944` -- in all 3 seeds.
* **theta0 fixed vs free makes no difference.** Repeating the tight-box test
  with `theta0` FIXED at 0 (matching the paper's own convention, removing
  that widening axis entirely) gave the SAME result (`x~0.030-0.031`, not
  A's `0.038944`).
* **Stage 2's OWN `set1` (A/B/C together) genuinely worked** (per the
  original #581 positive control: A/B/C all matched, `ic_dist` 0.0015-0.0048)
  and its own final population IS genuinely diverse (`x` spans
  `0.0067-0.0667`, `std=0.0175`, not collapsed) -- confirmed by directly
  inspecting `data/found/581_niching_ga/set01_final.npz`. What allows that
  diversity to persist (a specific seed, the chunked multi-invocation
  running pattern, or some other factor not yet isolated) was NOT
  identified in the time budget for this follow-on dispatch.

**Conclusion (matches Fable's own second-round diagnosis, now confirmed
empirically, not just theoretically):** `gurfil_kasdin_fitness` (Eq. 15) is
a pure boundedness measure with no periodicity/family-membership content;
once more than one family's neighborhood is reachable within a search box,
nothing in the objective steers the population toward a SPECIFIC published
point, and deterministic crowding does not reliably protect multiple
comparably-fit niches over a full 400-generation run on this landscape.
Box-narrowing (position-floor correctness + disjoint bounds) is necessary
but demonstrably NOT sufficient. This is reported as an honest, evidence-
based limitation of the (unmodified, per #583's own scope) objective +
niching combination -- not spun as a clean redesign win.

### Final design + validated result

A-F split to FULL single-family granularity (not even stage 2's own 3-family
A/B/C grouping) as the safest available mitigation -- removes the
family-vs-family competition pathway, though the diagnostics above show it
does NOT guarantee convergence to a family's specific point on any given
seed. 14 single-family `paper`-band partitions (A-N) + `DEEP_HILL` +
`BEYOND_HI_R` = 16 total.

Ran the FINAL single-family `C` partition (production script, full 400 gen,
pop 200, `data/found/583_widened_search/C_*`) as the clean validating
artifact:

```
[C] family C: MATCH ic_dist=0.0011(ok) type=DRO vs DRO rmin_ratio=0.983 rmax_ratio=0.978
[C] drift classification: 200/200 high-fitness members bounded at 50yr
[C] theta0 spot-check (5 members): all stable across phases
[C] reproduction: 1/1
```

This satisfies #583's own already-accepted positive-control bar ("recovers
its neighboring known family/families") cleanly and reproducibly. It does
NOT demonstrate that every single-family partition will recover its own
target on the first seed -- per the diagnostics above, that is NOT
guaranteed by this objective/niching combination. **Recommendation for the
coordinator's eventual full sweep:** run multiple independent seeds per
`paper`-band partition (not just 1) and treat "at least 1 known family
recovered per partition attempt" as the realistic bar, not "every family in
every partition on the first try."

### Reproduction (current, single-family design)

```
uv run python scripts/run_583_widened_bounded_drift_search.py --partition C --workers 8
uv run python scripts/run_583_widened_bounded_drift_search.py --analyze --partition C
```

## Not yet dispatched

The full widened-domain novelty sweep across all 16 partitions (A-N +
`DEEP_HILL` + `BEYOND_HI_R`), running multiple seeds per `paper`-band
partition per this dispatch's own recommendation above -- a longer job for a
future, separately-scoped dispatch, per #583's own gate.
