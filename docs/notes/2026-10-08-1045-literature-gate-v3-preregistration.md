# #1045 literature_check gate v3 — PRE-REGISTRATION (2026-10-08), committed before any code

Follows #972 v2 (`2026-10-07-972-literature-gate-v2-preregistration.md`, code 5fc5431d) and its Fable
re-review (v2 note sec. 9 and 9.1, findings 1-7, verdict PASS-WITH-NITS).

Lead ruling 2026-10-08:
- v3 fixes the two MAJOR findings, the three MINOR findings and the two NITs.
- Pinned tests; the 26 controls and the five candidates re-run; ruff, mypy; one more read-only review.
- The expected candidate statuses are IDENTICAL to v2. If any candidate changes, STOP and report.
- v2 is not used on any other candidate set (#1025, #1036) until v3 lands.

## 1. The fixes (one rule each; the probe becomes a pinned test)

| # | Finding | Rule | Probe (pinned test) |
|---|---|---|---|
| G1 | 1 (MAJOR) | A TOUR-ONLY anchor (declared `topology_label` non-empty and contained in {"mga-tour", "pump-tour", "ephemeris"}) cannot make a signature with NO topology label "published". Such a hit is scored, but its confidence is capped below `MATCH_THRESHOLD`, so at most it gives "inconclusive", naming the anchor. Signatures that declare a topology are unchanged (the #349 topology filter already applies). Chosen over a "topology-unlabelled" EXCLUSION: an exclusion would turn an untopologied TOUR candidate that matches a tour anchor into a false not-found. | Sun (E, V), no labels, no topology: not "published", and not citing Hughes as published. Sun (E, V, E), the same. A topology-declared mga-tour signature still matches a tour anchor as before. |
| G2 | 2 (MAJOR) | F7 fires on body-set INCLUSION: an anchor of the same primary whose body set CONTAINS the candidate's bodies (`seq_set <= a.body_set`, the `_candidate_anchors` structural test), with working-bodies (declared-different or unlabelled) among its exclusion reasons, makes a not-found or known-architecture-new-system result "inconclusive", naming the anchor. This replaces body-set EQUALITY. | Saturn ("Titan",) repeated-moon, unlabelled: inconclusive naming R-S Titan-Enceladus (v2: not-found). Jupiter ("Ganymede", "Ganymede") repeated-moon, unlabelled: inconclusive (v2: not-found). |
| G3 | 3 (MINOR) | `signature_from_review_entry`: a missing `verdict_audit['primary']` gives primary "" (never "Sun"). `check_literature` returns "inconclusive" for an empty primary ("primary unknown: rerun with the primary stamped"), before any search. | A review entry with no audit primary: inconclusive; a stamped one is unchanged. |
| G4 | 4 (MINOR) | The H&M anchor's `return_types_scope` gains "HR". Source: Menning 1968 ch. 2, pp.7-10 (digest `2026-10-05-digest-menning-1968-mit-thesis-earth-venus-periodic-orbits.md` lines 41-43): "The half-revolution return is a special case" of the full-revolution return. | Sun (E, V), two, {FR, HR}: published via H&M (v2: excluded by return-types). |
| G5 | 5 (MINOR) | The F13 ratchet also flags each of these: a `match` statement with a `"not-found"` case; a `.startswith("not-found")` / `.endswith(...)` call; a comparison of a name ending in `status` (or a `["status"]` / `.status` / `.get("status")` access) with a NAME rather than a literal, unless that name is `FRESH_STATUSES`. | Synthetic source snippets for each form: flagged. The current tree: clean. |
| G6 | 6 (NIT) | `_declared_scope_exclusion(s)` docstrings corrected (the list form returns [] when nothing excludes). | — |
| G7 | 7 (NIT) | test_f4 asserts status "not-found" and the footprint note. | — |

`phase_match.py:222` (also cited under finding 3) is NOT changed. It builds a phase signature from a
catalogue row, where an absent `primary` is the catalogue convention for heliocentric rows. It is not a
literature status.

## 2. Expected results (written before running)

- All 26 controls: "published", unchanged. Each control declares topology {"repeated-moon"} (G1 does
  not apply). Its body set equals an anchor's (G2 changes nothing for a published result). None has an
  HR return (G4). The Jones mutated-tag negative: unchanged.
- The five candidates: identical to v2.

| Candidate | v2 | v3 expected | Why unchanged |
|---|---|---|---|
| gc-1 | inconclusive (R-S G-C) | inconclusive (R-S G-C) | G2: equality already held |
| gc-2 | inconclusive (R-S G-C) | inconclusive (R-S G-C) | as gc-1 |
| ev-A | not-found | not-found | H&M ({E, V}, contains ev-A's bodies) excludes by return-types only (the same working label, "two"), so G2 does not fire. G4: ev-A has GEN, not HR. Hughes ({V, E, M}, contains them) declares no working scope. Topology declared (G1 n/a). |
| ev-B | not-found | not-found | as ev-A ({FR, FR-3:2, GEN}) |
| ev-C | inconclusive (H&M) | inconclusive (H&M) | F7 already fired at equality |

- The 393-row catalogue probe (the reviewer's method, unlabelled, v2 against v3), run as a check: no
  row may move from "published" or "inconclusive" to "not-found" or "known-architecture-new-system".
  The expected direction of every change is toward "inconclusive".

## 3. Process

- Implement G1-G7; pinned tests under tests/search/ (< 2 s each).
- Run: the controls and candidates (`scripts/litcheck_942_943_scope.py`); the catalogue probe; the
  literature test files.
- Then one full suite (tee'd, one at a time), ruff and full mypy; commit; one more read-only review.

## 4. Results (2026-10-08), recorded after the run

- Controls and candidates (`scripts/litcheck_942_943_scope.py`): all 26 controls return "published", and
  the Jones mutated-tag negative holds. The five candidates are IDENTICAL to v2 (gc-1, gc-2 and ev-C
  inconclusive; ev-A and ev-B not-found). `data/942_943_litcheck_scope.json` is byte-unchanged.
- Catalogue probe, v2 (5fc5431d) against v3:
  - 790 signatures: every row with a sequence, unlabelled, once with topology {repeated-moon} and once
    with none (`scratchpad/twobody-gen2-opus/probe_v2_v3.py`, `.out`).
  - 17 change, all published -> inconclusive, all in the no-topology half: the 15 H&M rows, ev-A and
    ev-C, which v2 cited to the Hughes one-shot free-return paper. That is finding 1, fixed by G1.
  - 0 rows move toward "not-found" or "known-architecture-new-system".
- Pinned tests: `tests/search/test_1045_literature_gate_v3.py` (G1-G5, 7 tests), test_f4 tightened
  (G7); the F13 ratchet over src/ and scripts/ stays clean under the wider G5 rule (0.6 s).
- Full set (tee'd, 2026-10-08, load average up to 22):
  - tests/data, tests/scripts and test_catalogue_rediscovery: EXIT 0 (14:15-14:21 AEDT).
  - tests/search: EXIT 1 (14:21-15:14). The only failures were the two load-sensitive Earth-Moon class 1
    tests, both 600-s timeouts: `test_known_close_pair_73c_plateaus_just_outside_guard` and
    `test_find_homoclinic_default_k_range_is_too_narrow_for_this_orbit`.
  - The lead re-ran both alone (15:15-15:19 AEDT, load 4.5): EXIT 0. Their serial call times on the
    #972 re-run were 133.4 s and 111.4 s.
- ruff clean; mypy src tests clean (939 files).
- Next: a read-only review of the v3 diff.
