# #972 literature_check gate v2 fixes — PRE-REGISTRATION (2026-10-07), committed before any code

Follows `2026-10-07-942-943-literature-gate-scope-preregistration.md` (v1, commit 2e16b56b) and the
Fable adversarial review (findings 2-7). Lead approval 2026-10-07: option (b), fixes 2-7. Whatever the
gate returns after these changes is the result. A waiver of spec sec. 16.5 for an "inconclusive"
candidate is the OWNER's decision only.

## 1. The fixes (one rule each; the probe that exposed it becomes a pinned test)

| # | Rule | Source / basis | Probe (pinned test) |
|---|---|---|---|
| F2 | `_architecture_anchors` returns nothing when ANY anchor with the same primary and the same `working_bodies_scope` has a body set containing the candidate's bodies, whatever its other scopes. A system the source treated is never "never-treated". | the meaning of #875 (ii), spec sec. 16.4 | Jupiter, (Ganymede, Callisto, Ganymede), working "one", topology {"resonant-hopping"}: today "known-architecture-new-system" citing R-S G-Io/G-E; must NOT return it |
| F3 | The alternating test first drops a final encounter equal to the first (a closing repeat), then reads the sequence cyclically. | the catalogue convention that repeats the closing body (e.g. "Callisto-Ganymede-Callisto-Europa-Callisto") | Jupiter, (G, C, G, C, G), working "two": today not-found; must return "published" (Campagnola GCGC), like (G, C, G, C) |
| F4 | A search hit synthesised from a corpus anchor (`anchor_name` set) is scored only if that anchor passes `_candidate_anchors` for the signature (same primary, the candidate's bodies within the anchor's body set, no declared-scope exclusion). Real web hits (no anchor identity) are unchanged. | the matcher's own structural-fingerprint rule (`_candidate_anchors` docstring: "A heliocentric Earth-Mars candidate cannot collide with a Jovian moon anchor and vice versa") | Sun, (E, J), unlabelled: today "inconclusive" 0.575 via Koon-Lo-Marsden (body set {Jupiter}); must not be scored from that anchor |
| F5 | Drop the Hughes anchor's `n_bodies_scope=3` tag (unconditional, lead). | the tag was not quoted from its source | — |
| F6 | Pinned tests under `tests/search/` (each < 2 s): each tag, positive and negative (the n-bodies, working-bodies, return-types and alternating exclusions, with a mutated-tag negative check), the new status, F2, F3, F4 and F7 probes, and the controls below. | lead/reviewer | the probes above |
| F7 | If an anchor is excluded ONLY by "working-bodies" and its body set EQUALS the candidate's body set, a literal result of "not-found" or "known-architecture-new-system" becomes "inconclusive", naming the anchor. A different-architecture object at a system the source treated goes to a human. | reviewer finding 7: the gc-2/GanCal#5 case | gc-2 |

## 2. A further finding, for the lead's decision (F8, NOT applied unless approved)

The Hughes anchor ("Hughes-Edelman-Longuski VEM cycler extensions (2014)") cites "AAS 14-822 /
'Venus-Earth-Mars Cyclers' extension paper (2014)", with provenance "inherited-unverified".
- A web search finds no AAS 14-822 by these authors.
- Their held 2014 paper is AIAA 2014-4109, "Fast Mars Free-Returns via Venus Gravity Assist". It is
  one-shot free returns, not cyclers (CORPUS_INDEX rows 152 and 398; digest of the JSR 2015 version).
- With F5 the anchor is untagged and body-set {V, E, M}, so it would claim any E-V candidate as
  "published" (confidence 0.85) on an unverified citation.
- PROPOSED F8: ground the anchor to AIAA 2014-4109 (citation, doi 10.2514/6.2014-4109, provenance
  verified-against-source after a page check), with topology {"mga-tour"} (one-shot free returns).
- Expected outcomes are given below both WITH and WITHOUT F8.

## 3. Positive controls (each must return "published")

- The 9 recovered members: GanCal#5, GanCal#1, GanEur#43, EurGan#131, GanEur#316, VenMar#45, Hollister
  1H, Hollister 2H, R-O 2.5.1.+0 (labels derived as in v1 + A1).
- The 15 catalogued H&M rows, unlabelled.
- NEW triple-cycler control: catalogue row `jones-2017-vem-emevve-outbound` (Jones, Hernandez & Jesick
  2017, AAS 17-577; held, mined in `2026-06-05-jones-aas17-577-vem-mining.md`).
  - Signature: primary Sun, sequence E, M, E, V, V, E (3 bodies), k = 2, V_inf from the row, unlabelled.
  - Expected: published via the Jones VEM anchor (n_bodies 3).
  - Negative: with the Jones tag mutated to n_bodies 2, it must NOT match Jones.
- Also: `liang-2024-cgcec-111-highperijove` (sequence C, G, C, E, C): expected published via Liang
  (n_bodies 3).

## 4. Expected candidate results (written before running)

| Candidate | Labels | Without F8 | With F8 | Why |
|---|---|---|---|---|
| gc-1 | two; {FR, SY} | INCONCLUSIVE | INCONCLUSIVE | F7: R-S Ganymede-Callisto (one) is excluded only by working-bodies and has the same body set {G, C} |
| gc-2 | two; {SY} | INCONCLUSIVE | INCONCLUSIVE | F7, the same anchor (the lead's stated expectation) |
| ev-A | two; {FR, GEN} | PUBLISHED via the Hughes anchor (0.85) | NOT-FOUND | F4 removes the Aldrin floor; without F8 the untagged Hughes {V, E, M} anchor matches; with F8 it is excluded by topology. H&M is excluded by return-types (not working-bodies, so F7 does not apply); F2 blocks the architecture match because H&M (two) contains {E, V} |
| ev-B | two; {FR, FR-3:2, GEN} | PUBLISHED via Hughes | NOT-FOUND | as ev-A |
| ev-C | one; {GEN} | PUBLISHED via Hughes | INCONCLUSIVE | with F8: F7 applies, because H&M (two) is excluded only by working-bodies and has body set {E, V}. No "known-architecture-new-system": the gate escalates the same-system, different-architecture case to a human |

So after these fixes, ev-A and ev-B (with F8) would pass the gate. gc-1, gc-2 and ev-C would be
"inconclusive": rows for gc-1 and ev-C then need the owner's decision. gc-2 is known-class-member by
the owner's ruling, so the gate's escalation is consistent with it.

Note for the owner (not applied): R-S AAS 07-118 p.18 names the massive target as future work, so
gc-1, gc-2 (at R-S's Ganymede-Callisto) and ev-C are case (i) or (ii) of #875 by the human reading. A
sourced "future-work architecture" field on the R-S anchors could let the gate express that. It is not
part of this pre-registration.

## 5. Process

- Implement F2-F7 (and F8 only if the lead approves).
- Run the controls and candidates (`scripts/litcheck_942_943_scope.py`, extended with the triple-cycler
  controls), the pinned tests, one full suite (tee'd), ruff and full mypy.
- Re-review (Fable). Its verdict is recorded here before any row is drafted.

## 6. Amendment (lead addendum 2026-10-07, after the reviewer's full report), committed before any code

| # | Rule | Basis | Probe (pinned test) |
|---|---|---|---|
| F13 | Seven callers compare the status literally (`low_thrust_cycler_search.py:319`, `cislunar_bct_search.py:216`, `precursor_matcher.py:135`, `scripts/campaign_468_multirev_tour.py:339`, `scripts/verify_327_umbriel_silver.py:492`, `scripts/run_299_lit_check_3d_family.py:215`, `scripts/run_435_high_e_er3bp.py:253`). They go through ONE helper, `literature_check.is_literature_fresh(status)` = status in {"not-found", "known-architecture-new-system"}. run_435's fixed-key counter gets the new key. | reviewer finding 13 | a test that the helper accepts both fresh statuses and rejects "published" and "inconclusive"; and an AST/grep ratchet over src/ and scripts/ that no code compares a literature status to the literal "not-found" outside the helper |
| F14 | An anchor that declares an architecture scope (`working_bodies_scope` or `return_types_scope`) is matched only by a signature that DECLARES the corresponding label. An unlabelled signature cannot be checked against that scope, so it does not match that anchor. | reviewer finding 14: an unlabelled Sun {V, E} signature (e.g. a VEM row with a partial derived sequence) must not be cited to Hollister/H&M at 0.85 | Sun, (V, E), unlabelled: today "published" 0.85 via Hollister/H&M; must not match H&M |

Consequences of F14 for the controls (rule fixed now):
- The 15 catalogued H&M rows are then labelled MECHANICALLY from the sourced per-encounter table
  `data/sources/hollister-menning-1970-table3.yaml` (Table 3 transcription):
  - working_bodies "two" if both planets have a turn theta >= 0.05 deg somewhere;
  - return types from each consecutive same-planet encounter interval divided by that planet's period
    (Earth 365.25 d, Venus 224.70 d): within 5 % of 1 -> "FR"; in (1, 2) -> "SY"; else "GEN".
  - Expected: all 15 labelled "two" with return types within {FR, SY}, and "published" via H&M.
- The R-S anchors (working_bodies_scope "one"): every recovered R-S control is labelled, so unaffected.
- F14 also applies to the GCGC anchor (working "two"). No control is unlabelled there.

Hygiene (reviewer 11, 12; lead):
- H11: v1 note sec. 2 (the 1e-6-deg threshold) is marked SUPERSEDED by its sec. 8 (A1). A1 states that
  ev-C's label flips ("two" to "one") and its outcome does not (inconclusive 0.475 under either label;
  reviewer probe). Done in this commit.
- H12: the R-S 2009 anchor comments say the AAS 07-118 statement is quoted for the 2009 rows (the same
  rows and nomenclature, 2007 digest sec. 0). The R-S 2009 Titan-Enceladus anchor gets the same
  `working_bodies_scope="one"` (the same source statement covers it).

Expected candidate results: unchanged from sec. 4 (every candidate is labelled, so F13/F14 do not
apply to them).

## 7. Lead rulings 2026-10-07 (F7 reading, F8), committed before the F7/F8 code

**F7: correction of the sec. 4 expectation, and the ruled reading.**
- Discrepancy found before coding: sec. 4 says H&M excludes ev-C "only by working-bodies". That is
  wrong. ev-C is labelled working "one", return types {GEN}. H&M declares working "two" AND return types
  {FR, SY}, so it excludes ev-C by working-bodies AND by return-types.
- Under the literal F7 text ("ONLY by working-bodies"), F7 would not fire for ev-C. ev-C would then
  return known-architecture-new-system via the R-S VenMar/VenMer "one" anchors and pass the gate.
- Lead ruling (b), the intent reading. RULE TEXT: "if working-bodies is AMONG the exclusion reasons
  for an anchor whose body set equals the candidate's, the result is inconclusive, listing that anchor."
  Basis: a different-architecture object at a system the source treated goes to a human. H&M treated
  Earth-Venus cyclers. Choosing the reading that passes our own candidate after seeing the discrepancy
  would be tailoring.
- Expected for ev-C: INCONCLUSIVE (H&M listed), then the owner's recorded decision. gc-1, gc-2, ev-A and
  ev-B are unchanged under either reading.

**F8 approved.** Page check of the held PDF
(`hughes-edelman-longuski-2014-fast-mars-free-returns-venus-ga-AIAA-2014-4109.pdf`, p. 1):
- title "Fast Mars Free-Returns via Venus Gravity Assist";
- AIAA 2014-4109, AIAA/AAS Astrodynamics Specialist Conference, 4-7 Aug 2014, San Diego;
- DOI 10.2514/6.2014-4109 (printed on the page);
- authors Hughes, Edelman, Longuski, Loucks, Carrico, Titok.
- The abstract describes one-shot free returns on E-V-M-E (and Venus-after-Mars) paths, 2015-2060
  launches. It does not describe cyclers.

The anchor is grounded to this paper:
- body set {V, E, M} (kept); topology {"mga-tour"}; no n-bodies tag (F5);
- citation and keywords rewritten to the paper; provenance verified-against-source.

The old anchor comment "#350: extends Jones-Hernandez-Jesick AAS 17-577 VEM cycler family" was
anachronistic: a 2014 citation cannot extend a 2017 paper. The lead checked the corpus: no held Hughes
work treats VEM cyclers. The 2016 thesis preview's cycler chapters (IEG triple cyclers; Earth-Mars
establishment) are covered by other anchors. The JSR 2015 version credits Hollister/H&M for E-V cyclers.

**Parked, NOT implemented:** a sourced "future-work architecture" field on the R-S anchors. R-S AAS
07-118 p. 18 names the massive-target architecture as future work. Such a field would let the gate
express #875 (i)/(ii) for gc-1/gc-2/ev-C mechanically. Recorded as an option only.

**After the v2 run (lead):** gc-1 and ev-C "inconclusive" is the spec's human-review path. The lead
takes the gate diagnosis to the owner for a recorded decision. No row is written before that.

## 8. Results and after-the-fact amendment (2026-10-08), recorded after the run

**Run** (`scripts/litcheck_942_943_scope.py`, output in `data/942_943_litcheck_scope.json`):
- All 26 controls return "published":
  - the 9 recovered members;
  - the 15 H&M rows, all labelled "two" {FR, SY};
  - Jones VEM and Liang CGE.
- The Jones mutated-tag negative (n_bodies 2) does not match Jones.
- The candidates match the pre-registered expectations (sec. 4 with F8, sec. 7 for ev-C):

| Candidate | Result | Anchor named |
|---|---|---|
| gc-1 | inconclusive | R-S 2009 Ganymede-Callisto (F7) |
| gc-2 | inconclusive | R-S 2009 Ganymede-Callisto (F7) |
| ev-A | not-found | |
| ev-B | not-found | |
| ev-C | inconclusive | H&M (F7, ruling (b)) |

So ev-A and ev-B pass the gate. gc-1, gc-2 and ev-C go on the owner path.

**Deviations (stated, accepted by the lead 2026-10-08):**
1. F13: there were 12 literal comparison sites, not 7. All go through `is_literature_fresh`. The
   five extra sites are `gauntlet_run_274`, `branch_c32_b0_p389_2_gates`,
   `run_301_subfamily_validation`, `literature_check_review_queue` and `run_432_er3bp_discovery`.
   The AST ratchet covers all of src/ and scripts/. Positive control: it flags the old code at
   `branch_c32...:409` and `precursor_matcher.py:135`.
2. F14 H&M labels: computed from `hollister_menning_1970.load_table3`, which is the Table 3 YAML with
   that module's documented print-error date fixes (they predate #972). The lead ruled this the right
   source. With the raw printed dates, orbits 2, 6 and 8 give a GEN interval, each at a documented
   print error (e.g. orbit 8 "4870" for 5870). With the fixes, all 15 are {FR, SY}, which agrees
   with the module's `block_types`.
3. F14 + F7: an unlabelled signature against a working-scope anchor gives the reason
   "working-bodies-unlabelled". When the body set is equal, F7 makes the result "inconclusive", not
   "not-found". This is the conservative combination.

**Consequence NOT predicted by this pre-registration:** under F14 (with H12's Titan-Enceladus tag),
every old pipeline that feeds UNLABELLED signatures loses the R-S citation on R-S rediscoveries.
- Same body set: the result is "inconclusive" via F7. Otherwise, the R-S author bonus is lost on web
  hits.
- The direction is safe: the result is never a false "not-found".
- The fix is to label signatures at the source. Registered by the lead as #1035: label the
  signatures in the old callers (`low_thrust_cycler_search`, `cislunar_bct_search`,
  `precursor_matcher`, the `run_*` scripts), derived mechanically from turns as in the scope script.

Two existing tests in `tests/search/test_literature_check.py` were retrofitted (lead ruling (A); not
exempted). Every other assertion is kept:
- `test_russell_strange_2009_double_cyclers_flagged_published` now uses labelled "one" signatures
  and still expects "published" via R-S.
- The new `test_russell_strange_2009_unlabelled_signature_not_cleared`: unlabelled Ganymede-Io
  returns "inconclusive" and names the R-S paper (title and DOI URL).
  - The unlabelled Titan-Enceladus case on that test's web backend still reads "published" (0.75).
  - Reason: that hit has no anchor identity and scores on its text alone. Real web hits are
    unchanged by design (F4).
- `test_880_same_scope_anchor_still_flags_published` now uses a labelled "one" Titan-Enceladus
  signature and expects "published".
- The new `test_880_unlabelled_same_system_anchor_goes_to_a_human`: unlabelled returns
  "inconclusive", naming the R-S Titan-Enceladus anchor.

**Full suite** (tests/data tests/search tests/scripts, tee'd, 2026-10-07 21:07 to 2026-10-08 00:47
AEDT, load average up to 48):
- 14 failures: the 2 above, plus 12 pytest timeouts in heavy CR3BP tests. Those tests do not import
  `literature_check`; they mention it only in docstrings.
- The 12 were re-run outside the suite (call times; setup in brackets where it was large):
  - Re-run 1, 2026-10-08 01:10-01:29 AEDT, 4 workers, load average 34-70: 10 passed.
    - Neptune-Triton families: `test_continue_32_esm4_hc1_family_777_reaches_jacobi_bound` 920.6 s
      (its own 1200 s limit), `test_continue_23_family_reaches_jacobi_bound_and_matches_printed_members`
      414.9 s, `test_continue_47_esm4_pair_777_reaches_jacobi_bound` 212.9 s.
    - Neptune-Triton connections: `test_find_homoclinic_returns_known_primary_combo` 255.3 s.
    - Earth-Moon class 1: `test_known_close_pair_73b_plateau_independent_of_fd_step` 379.5 s (+172.8 s
      setup), `test_find_homoclinic_narrow_diagonal_window_is_honest_empty` 149.9 s (+77.0 s).
    - #549 binary sweep: `test_549_orcus_vanth_32_clean_negative_aware` 278.2 s,
      `test_549_didymos_dimorphos_11_clean_negative_aware` 121.7 s.
    - Jovian: `test_find_homoclinic_default_target_unchanged` 221.8 s,
      `test_5_6_li_continues_smoothly_toward_c_flyby` 215.6 s.
    - Two timed out again at 600 s under that load: `test_known_close_pair_73c_plateaus_just_outside_guard`
      and `test_find_homoclinic_default_k_range_is_too_narrow_for_this_orbit`.
  - Re-run 2, 2026-10-08 08:01 AEDT, serial (-n 0), load average about 4: both passed, 133.4 s and
    111.4 s (+16 s setup each).
  - All 12 pass. Their timeouts in the suite were CPU contention.
