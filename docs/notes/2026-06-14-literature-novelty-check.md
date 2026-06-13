# Literature-novelty check for discovery-daemon SILVER candidates (#261)

**Date:** 2026-06-14
**Status:** built + self-validated (offline deterministic tests + live WebSearch)
**Files:** `src/cyclerfinder/search/literature_check.py`,
`scripts/literature_check_review_queue.py`,
`tests/search/test_literature_check.py`

## Problem

The #253 discovery daemon's "novel" verdict only deduplicates a candidate
against (a) our 277-row **sourced** catalogue and (b) the method-versioned
negative registry. Both are a *subset* of the published cycler literature, so a
daemon SILVER could be a cycler **published elsewhere** — a rediscovery
mislabelled novel. The `review_queue` schema already reserved a
`literature_check` field (`None`). This task builds the check that populates it,
so no candidate can be claimed novel until it has been searched against the
published record.

## Design

### The signature (structural fingerprint, not numbers)

`CandidateSignature` = the *physics identity* of a trajectory:
`primary` (Sun / Jupiter / …), the body-or-moon `sequence` (encounter tour),
`period_k` / `period_years`, `vinf_per_encounter_kms` (the V∞ regime), and the
per-leg `resonances` + `n_rev` (the multi-rev structure). Two trajectories with
the same fingerprint are the same family regardless of print-precision digits.

### The check (`check_literature`)

1. `build_queries` turns the fingerprint into an ordered query trail,
   most-specific → most-generic: exact tour + period; resonance / multi-rev
   structure; V∞ regime; named-corpus author/keyword anchors for overlapping
   families; then a generic arXiv/ADS/NTRS fallback.
2. Each query is run through an **injected** `SearchFn` (the real WebSearch in
   the live path; a deterministic corpus in tests). The search backend is
   injected, never hardcoded, so the module never silently "passes" with no
   search performed.
3. `_result_matches_fingerprint` scores every hit on *structural* overlap with
   the signature — mandatory "cycler" floor, then each tour body named, then a
   curated-corpus author/keyword, then the primary — **not** raw keyword match.
4. Verdict:
   - any hit ≥ `MATCH_THRESHOLD` (0.70) → **published** (+ citation/doi),
   - best hit in `[0.45, 0.70)` → **inconclusive** (cycler-adjacent but
     unconfirmed — a human must adjudicate),
   - cycler-adjacent material absent, with searches that DID return results →
     **not-found**,
   - searches that returned *nothing at all* → **inconclusive** (we never
     actually searched, so a "not-found" would be a lie).

`KNOWN_CORPUS` is a hand-curated set of `CorpusAnchor`s (Aldrin; Russell-Ocampo /
McConaghy-Landau SnLm; Liang CGE; Hernandez/Jones/Jesick IEG; Strange/Campagnola/
Russell moon-tour & V∞-leveraging; Jones VEM). Each carries the **publication**
facts (authors / venue / DOI) — these are citation facts, never values our own
code computed.

### The SILVER gate

`is_novelty_claimable(literature_check) -> bool` is the gate the daemon / review
flow MUST consult:

> novelty is claimable **iff** the block is populated (`checked is True`) AND its
> status is not `published` AND not `inconclusive`.

An unpopulated (`None`) block, a `published` status, or an `inconclusive` status
all return `False`. `to_review_block()` maps `not-found → result:"no-match"` so
the pre-existing `review_queue.is_promotion_eligible` gate (which keys on
`result == "no-match"`) lines up unchanged.

### The driver

`scripts/literature_check_review_queue.py` walks `data/review_queue.jsonl`,
builds each SILVER entry's signature (`signature_from_review_entry`), runs the
check, and rewrites the entry with `literature_check` populated. Backends:

- `--backend offline` (default): conservative corpus-only filter — flags
  candidates inside a known published family as `published`, everything else
  `inconclusive` (an offline run never emits a false `not-found`).
- `--backend results-json PATH`: replays pre-collected **live** WebSearch hits
  (`{query: [{title,url,snippet}]}`) — the path that produces a *trustworthy*
  `not-found` (a real search returned nothing).

No catalogue writeback; only the review queue is touched.

## SELF-VALIDATION (reproduce-before-trust)

### Deterministic offline tests

`tests/search/test_literature_check.py` (14 tests, all pass) feeds the checker
a fixed corpus mirroring the real published hits and asserts:

| signature | expected | got |
|---|---|---|
| Aldrin E–M (k=1, V∞ 6.5/9.7) | published | **published**, DOI 10.2514/6.2002-4420 |
| Russell/McConaghy E-E-M-M (k=2, V∞ 4.1/2.0) | published | **published**, Russell&Ocampo 2004 / McConaghy 2006 |
| Liang CGE (Callisto-Ganymede-Callisto-Europa-Callisto) | published | **published**, DOI 10.2514/1.G008387 |
| fabricated Neptune Triton-Nereid-Proteus | not-found | **not-found** |
| empty-search-everywhere | inconclusive | **inconclusive** |

Plus gate tests (unpopulated / published / inconclusive → not claimable;
clean not-found → claimable).

### LIVE WebSearch self-validation (real queries + verdicts)

Run 2026-06-14 with the harness WebSearch tool:

1. **Aldrin** — query `Aldrin Earth-Mars cycler trajectory gravity assist`
   → top hits: Wikipedia "Mars cycler" (names the Aldrin cycler, 1L1), AIAA
   `arc.aiaa.org/doi/10.2514/3.25519` "Cycler orbit between Earth and Mars",
   Purdue "Powered Earth-Mars Cycler". **Verdict: published** (structural match —
   "Earth-Mars cycler" + Aldrin named). ✔

2. **Liang CGE** — query
   `Callisto-Ganymede-Europa triple cycler Jovian moon trajectory`
   → top hit: "Callisto-Ganymede-Europa Triple Cyclers" (the Liang/Yang et al.
   paper, our DOI 10.2514/1.G008387), describing the exact
   Europa→Ganymede→Callisto→return-to-Callisto multi-rev Lambert construction.
   **Verdict: published.** ✔

3. **Russell/McConaghy SnLm** — query
   `ballistic Earth-Mars cycler systematic catalog Russell Ocampo two synodic`
   → top hits: "Systematic Method for Constructing Earth-Mars Cyclers Using
   Free-Return Trajectories", "Notable Two-Synodic-Period Earth-Mars Cycler",
   "Analysis of Various Two Synodic Period Earth-Mars Cycler Trajectories"
   (AIAA 2002-4423). Search text even quotes "24 purely ballistic cyclers with
   periods of two to four synodic periods". **Verdict: published.** ✔

4. **Fabricated Neptune** — query
   `Triton Nereid Proteus Neptune moon cycler trajectory repeating`
   → all hits are moon-*discovery* / orbital-origin pages (Wikipedia "Moons of
   Neptune", NASA Proteus, A&A capture dynamics); **zero** carry "cycler
   trajectory" content. The search engine itself flagged: "results focus on
   orbital characteristics … rather than specific cycler trajectory
   information." No hit clears the "cycler" structural floor. **Verdict:
   not-found.** ✔

**Conclusion: the check correctly flags all three known-published cyclers as
`published` with plausible real citations, and correctly returns `not-found` for
the fabricated signature.** The live behaviour matches the deterministic test
corpus, so the offline tests are a faithful proxy. The check is trustworthy
enough to USE as a rediscovery filter.

## DISCIPLINE — read before trusting any `not-found`

**`not-found` is NECESSARY-NOT-SUFFICIENT for novelty.** Absence of a search hit
is *not* evidence of absence in the literature: web search is partial,
conference proceedings and paywalled venues are under-indexed, and structural
fingerprints do not always survive into searchable text. This module is a
**FILTER against obvious rediscoveries, not a novelty proof.** A `not-found`
only *clears* a candidate to continue to the human + V0–V5 gauntlet; it never
*certifies* novelty. The human and the gauntlet still govern. Citations are by
publication (DOI / arXiv / public URL) only — never a private repo path.

## Decisions under ambiguity

- **Offline default emits `inconclusive`, not `not-found`, for unknown
  candidates.** An offline run has not searched the web; emitting a clean
  `not-found` there would be a false bill of novelty. Trustworthy `not-found`
  requires the `results-json` backend fed by real WebSearch hits.
- **`inconclusive` is NOT novelty-claimable.** A weak/ambiguous search result is
  treated like an unfinished check, forcing a re-run or human look rather than
  defaulting to "novel".
- **WebSearch is a harness tool, not a Python import**, so the checker takes an
  injected `SearchFn` and the live search is run deliberately (recorded to
  `results-json`) rather than fired inline from the quota-bound daemon.
- The review queue was empty at build time (the #259 purge removed the four
  spurious all-zero-rev Jovian SILVERs), so the driver was validated by dry-run
  + unit tests rather than on live rows.
