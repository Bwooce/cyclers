# #588: two Gurfil-Kasdin companion papers, digested and checked against the 20 unmatched clusters

**Date:** 2026-07-14
**Trigger:** live literature-check phase of #588 (dedup + novelty adjudication of the 264-then-45-
candidate pool from #583/#586's widened Sun-Earth ER3BP sweep). Cross-referencing the 45 deduped
clusters against the 14 known Gurfil-Kasdin families' own r-bands left 20 clusters with no
same-type/overlapping-range match. A live search for other Gurfil-Kasdin Sun-Earth work surfaced two
candidate papers not previously in this project's corpus. Both were independently confirmed via
CrossRef (title/authors/journal/volume/page/DOI, not just search-engine synthesis) before acquisition,
per [[feedback_ground_citations_against_content]]'s same-author/same-year collision warning, then
read in full (user-supplied PDFs).

## Paper 1: Gurfil & Kasdin (2003), SPIE Proc. 4854:251-261, DOI 10.1117/12.459820

"Practical Deep-Space Geocentric and Out-of-Ecliptic Orbits in the Sun-Earth Restricted
Three-Body Problem." **Confirmed to be the SAME 14-family census already anchored in this repo**
(`gurfil-kasdin-2002-geocentric-er3bp` in `literature_check.py`, reproduced in
`scripts/run_581_gurfil_reproduction.py`), not independent data — the conference precursor to the
CMAME 191:5683-5706 journal paper, before its title changed. Verified by exact numeric match: this
paper states Family F's minimum Earth approach as "224,900 km" and Family C's orbit as "within
almost a constant distance of 1 million km" — both match `run_581_gurfil_reproduction.py`'s
`TABLE34` dict verbatim (Family F rmin=224,900 km; Family C rmax≈1,126,563 km). Its own reference
[27] ("Practical Geocentric Orbits in the Spatial Elliptic Restricted Three-Body Problem," submitted
to JAS) is almost certainly an earlier working title for the paper that became the CMAME 5683-5706
journal article. **No new corpus gap; closes zero of the 20 unmatched clusters.**

## Paper 2: Gurfil & Kasdin (2002), CMAME 191:2141-2158, DOI 10.1016/S0045-7825(01)00380-2

"Characterization and design of out-of-ecliptic trajectories using deterministic crowding genetic
algorithms." A genuinely different paper from the already-anchored CMAME 191:5683-5706 one — same
journal/year/authors, different page range, different problem formulation. This one:

- Uses the **circular** (not elliptic) Sun-Earth RTBP, Earth-centered rotating frame.
- Objective function is **maximize z² (out-of-ecliptic displacement)**, not periodicity/boundedness
  — a fundamentally different search than the family-census paper's `gurfil_kasdin_fitness`.
- Search box extends out to **1-5 AU from Earth** (`r_max` in Table 2) with initial speeds up to
  15.4 km/s — far beyond the family census's ~5.8e5-1.17e7 km geocentric footprint.
- Reports only a **handful of individual best-fitness trajectories** per constraint set (8
  characterization sets + 2 design points), classified qualitatively into Type I/II/III by relative
  in-plane vs. out-of-plane frequency — not a bounded-family catalogue with reproducible per-family
  r-bands the way the census paper is.
- Its reported trajectories reach z-excursions of 0.223-0.374 AU (33.4M-56.0M km) and Earth
  distances up to several AU — an order of magnitude beyond our largest unmatched cluster (cluster
  36, rmax=40.8M km ≈ 0.27 AU).

**Conclusion: real, distinct, now correctly in corpus — but does not structurally cover our 20
unmatched clusters.** The search formulations don't correspond closely enough to check family-by-
family (this paper doesn't report families, just isolated maximum-displacement points), and the
scale/objective mismatch (unbounded z-maximizing drift over an Earth-departure box out to 5 AU, vs.
our bounded/periodic DRO-DPO-ERO-DEO characterization) means "not found in this paper" is a weak
signal either way — filed as a corpus anchor for future widened-domain searches to check against
(especially any future search that reaches the ~30-60M km z-excursion regime), not treated as a
closure of the current gap.

## Net effect on #588

Neither paper reduces the 20-cluster unmatched set. The live-literature-check phase (mandatory per
[[feedback_literature_novelty_check_baseline]]) surfaced two real companion works, both now
digested, filed in `cyclers_pdf/papers/` and `CORPUS_INDEX.md`, and one (the out-of-ecliptic paper)
added as a new `CorpusAnchor` in `literature_check.py` so future searches into that z-excursion
regime engage it structurally. The 20 clusters proceed to Opus/Fable adjudication as originally
scoped, with this search now honestly exhausted rather than skipped.
