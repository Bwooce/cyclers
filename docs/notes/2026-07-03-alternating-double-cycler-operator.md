# Build note: alternating-double-cycler genome operator (#526 Part 2)

**Task:** #526 asked to (1) ground the citation Liang et al., JGCD 2024,
DOI 10.2514/1.G008387, and (2) encode its "alternating-double-cycler"
construction as a reusable genome operator in `src/cyclerfinder/genome/`.

## 1. Citation grounding (independent re-verification, 2026-07-03)

The project was previously burned by a same-author/same-year citation
collision (`feedback_ground_citations_against_content`), so #526 explicitly
asked for independent verification before trusting the citation. Checked
against THREE independent sources, none of which is this project's own
prior work:

* **AIAA/ARC** (`https://arc.aiaa.org/doi/10.2514/1.G008387`, fetched via
  browser 2026-07-03): title "Callisto-Ganymede-Europa Triple Cyclers",
  authors Guoliang Liang, Hongwei Yang, Shuang Li, Xiaoli Bai, Limin Qin;
  *Journal of Guidance, Control, and Dynamics*, Vol. 48 No. 1 (Jan 2025),
  published online 12 Sep 2024; a Technical Note; a portion presented as
  Paper IAC-23,C1,9,9,x76777, 74th IAC, Baku, 2-6 Oct 2023. The publisher's
  free first-page preview (page 146 of the printed volume) additionally
  confirms the Sec. I/II content matches this project's existing mining
  note word-for-word in substance (near-resonance of Callisto-Ganymede vs.
  Ganymede-Europa synodic periods, Eqs. 1-2).
* **Semantic Scholar** (`api.semanticscholar.org` DOI lookup): same title,
  same 5-author list (Guoliang Liang, Hongwei Yang, Shuanglin Li [sic,
  Semantic Scholar's OCR of "Shuang Li"], Xiaoli Bai, Limin Qin), venue
  "Journal of Guidance Control and Dynamics", year 2024.
* **Unpaywall** (`api.unpaywall.org` DOI lookup): confirms closed access
  (`is_oa: false`, `has_repository_copy: false`) with corresponding
  author Guoliang Liang, Nanjing University of Aeronautics and
  Astronautics -- i.e. there is exactly ONE paper at this DOI and no
  legitimate free full copy exists anywhere (ruling out a mix-up with a
  different, possibly open, paper).

**Verdict: the citation is real and says what the task assumed.** No
same-author/same-year collision (cf. the Hernandez-Jones-Jesick 2017
VEM-vs-Galilean-triple-cycler precedent this project already knows about).

## 2. The paper was already filed, digested and mined under #216

Re-checking `docs/notes/CORPUS_INDEX.md` before doing anything else (the
corpus-document-policy discipline) found the paper is **already**:

* Filed: `papers/liang-2024-callisto-ganymede-europa-triple-cyclers-JGCD.pdf`
  (text-layer).
* Deeply digested: `docs/notes/2026-06-11-liang-2024-cge-triple-cyclers-mining.md`
  (Eqs. 1-20, Tables 1-7, character-by-character rescanned 2026-06-12) and
  `docs/notes/2026-06-13-liang-abc-reproduction.md` (numeric reproduction
  results).
* Reproduced exactly: `src/cyclerfinder/search/cge_scaffold.py` (Members
  A/B/C from their printed Tables 2-7) and consumed by
  `src/cyclerfinder/search/moon_cycler_genome.py` (the repeated-moon
  multi-revolution SEARCH genome, `MoonSystem` / `MoonCyclerGenome` /
  `EncounterGene` / `LegGene`, plus the reproduce-before-search gate).
* Catalogue rows: `liang-2024-cgcec-{111-highperijove,110-highperijove,
  111-lowperijove,ephemeris-2033}` (V1/V0).

None of that is duplicated here. What #526 actually asks for -- and what
was missing -- is a REUSABLE operator generalizing the paper's
construction *strategy* beyond the one hard-coded Callisto-Ganymede-Europa
case.

## 3. What was built

`src/cyclerfinder/genome/alternating_double_cycler.py` (+
`tests/genome/test_alternating_double_cycler.py`, 13 tests):

* `synodic_period_days` / `mean_motion_rad_day` -- primitives.
* `analyze_near_resonance(body_a, body_b, body_c, mean_motions, ...)` --
  generalizes Liang's Eq. 1-2 near-commensurability analysis to ANY
  3-body chain, not hard-coded to Jupiter's moons. Uses a
  lowest-order-continued-fraction-convergent search (see design note
  below), not a naive "closest fraction with denominator <= N" search.
* `build_alternating_double_cycler_seed(system, half_cycle_1, half_cycle_2,
  n_cycles=1, legs=None)` -- the "key trick" (Sec. III.B) as a composable
  operator: stitches two hub-sharing double-cycler encounter half-sequences
  into one repeating `MoonCyclerGenome`, dropping the duplicated junction
  hub encounter at each seam.
* `liang_cge_near_resonance()` / `liang_cge_alternating_seed()` -- the
  sourced positive-control anchors (below).
* Added `SATURN_MOONS` / `saturn_system()` to
  `search/moon_cycler_genome.py` (small additive registry helper) to
  support the reusability demonstration on a second moon system.

### Design pitfall caught and fixed: naive closest-fraction search is WRONG

The first implementation used
`Fraction(ratio).limit_denominator(max_denominator)` directly. Run on
Liang's own Table 1 numbers this returns **16/9** (mismatch fraction
0.09%), not the paper's published **7/4** (mismatch fraction ~1.5%) --
16/9 IS the numerically closest fraction with denominator <= 20, but it is
NOT what the paper reports. This would have been a silent, wrong
generalization (exactly the fabrication risk flagged in the module's own
"what this module does NOT do" section). Root cause: Liang et al. pick the
LOWEST-ORDER commensurability that is already tight (Eq. 3, p. 3, ties the
quasi-period to small-integer multiples of ALL THREE orbital periods, not
just a blind best-rational-approximation of the two-period ratio) --
standard practice in resonance identification (real/physically-relevant
resonances are low-order; a higher-order match is usually numerical
coincidence). Fixed by scanning continued-fraction convergents in
increasing-denominator order and stopping at the first one whose
scale-free mismatch clears a tolerance (default 5%) -- this reproduces
7/4 for Liang's own numbers without hand-coding "4" anywhere. Caught by
running the positive control against the paper's own numbers BEFORE
trusting the implementation (the project's `check-don't-guess` /
`verify-a-gauntlet-with-a-positive-control` discipline in miniature).

### Positive controls (sourced-only, never circular)

1. `liang_cge_near_resonance()` on Liang's own Table 1 mean motions
   recovers `(p, q) == (7, 4)` and `mismatch_days ~= 0.7365` (Eq. 2,
   within ~3e-3 d of the paper's printed value -- the residual difference
   traces to using the exact `2*pi/dn` formula on the 4-decimal-quantized
   Table 1 inputs rather than the paper's own possibly-higher-precision
   intermediate values).
2. Independent (non-golden) cross-check: the SAME analysis run on this
   project's own JPL-SSD-sourced registry mean motions (NOT copied from
   the paper) recovers the same 7:4 resonance with a mismatch within 2.4%
   and a quasi-period within 0.004% of the paper's own numbers -- two
   independent data sources agreeing.
3. `liang_cge_alternating_seed()` stitched from the two published CGE
   halves (`Callisto-Ganymede-Callisto`, `Callisto-Europa-Callisto`)
   reproduces `moon_cycler_genome.CGE_SEQUENCE`
   (`Callisto-Ganymede-Callisto-Europa-Callisto`) exactly.
4. Reusability (explicitly NOT a novelty claim): `analyze_near_resonance`
   runs body-agnostically on a Saturnian chain (Enceladus-Dione-Rhea) via
   the registry and returns a well-formed result -- no resonance is
   asserted or claimed for that chain, this only demonstrates the
   function is not Jupiter-specific.

### What this deliberately does NOT do

* Does not re-implement Liang et al.'s Eq. 16 specific epoch-multiplier
  law or their Sec. III.C sequential Lambert/differential-evolution
  optimizer -- those stay exactly as already reproduced, unmodified, in
  `cge_scaffold.py`. Re-deriving a "generalized Eq. 16" without the full
  derivation in hand would risk fabricating a formula the paper never
  published.
* Does not construct a physical double cycler from scratch -- it composes
  two caller-supplied encounter half-sequences.
* No novelty claim, no catalogue writeback -- pure scaffolding.

## 4. Validation status

* `uv run pytest tests/genome/test_alternating_double_cycler.py -q` --
  13/13 pass.
* `uv run pytest tests/data tests/search tests/genome -q` -- full ratchet,
  run alongside this change (see #526 OUTSTANDING entry for the result).
* `uv run ruff check` / `ruff format --check` / `uv run mypy src tests` --
  clean.
