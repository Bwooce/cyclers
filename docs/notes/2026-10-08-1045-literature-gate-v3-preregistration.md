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

## 5. Read-only review verdict (fresh Fable reviewer, 2026-10-08)

**VERDICT: FAIL (narrow).** One MAJOR remains: G1 does not close finding 1 on the Jupiter pipeline.
- The 26 controls, the negative and the five candidates are unchanged (independently re-run;
  the JSON is byte-identical), so the row results stand.
- G1-G5 and G7 match the pre-registered rule text; the G1 problem is in the RULE, not the code.
- Probes, v2 against v3:
  - unlabelled: 790 signatures, reproducing the 17 changes;
  - labelled: 1364 signatures (Sun, Jupiter and Saturn rows, both topologies, working bodies "one"
    and "two"): 34 changes, all published -> inconclusive via G1, 0 toward fresh.

Findings (summarised; the full text is the reviewer's hand-back of 2026-10-08):
1. MAJOR. The Niehoff anchor is tagged {"mga-tour", "pump-tour", "resonant"}, so it is not
   "tour-only" and the cap never applies.
   - Jupiter (Ganymede, Callisto, Ganymede) untopologied, period 3 (the review-queue gancal shape)
     gives "published" 0.95 citing Niehoff, "Touring the Galilean Satellites". The same holds for
     (Europa, Ganymede), (Io, Europa) and the 4-moon sets.
   - Fix: cap anchors with any tour label and NO cycler-class label, or drop "resonant" from Niehoff.
     Pin a Jupiter probe.
   - Known limit, outside G1: Uranus (Umbriel, Titania) untopologied gives published via Pergola et
     al. (halo-tagged), and Pluto (Charon, Nix) via Howett (no topology label).
2. MINOR. The cap value shadows F7. Sun (E, V) untopologied gives inconclusive 0.69 citing the
   Hughes hit; H&M appears only in the excluded-anchor list. Fix: when the best hit is capped, also
   name the F7 anchors.
3. MINOR. G6 is incomplete: `_declared_scope_exclusions`'s first docstring line still says "or None".
4. MINOR. G5 evasions remain:
   - an alias (`st = r.status; st in FRESH`), the realistic one;
   - an Attribute set or constant (`mod.OTHER_SET`, `consts.NOT_FOUND`);
   - a dict, a walrus, `.startswith('not-')`.
   Extended to tests/, the rule would flag 11 legitimate literal asserts in 7 files (report only).
5. NIT. G4 (HR) is faithful and has a positive control: the HV form of Hollister 1H
   (k2|RE/1:1|LE>V/0s|HV/1,0,a|HV/3,1,a|LV>E/0s) reads not-found under v2 and published via H&M
   under v3. So v2 gave a false not-found there. 3 of the 6 HR-bearing gauntlet candidates flip
   not-found -> published; none is a writeback candidate. Fix: the anchor comment should say Menning
   p.42 names HR variations but Table 3 computes none.
6. NIT. G2: `seq_set <= a.body_set` is vacuously true for an empty sequence. Fix: add `seq_set and`.
7. G3: no caller relies on the old "Sun" default.
8. Candidates confirmed: gc-1, gc-2, ev-C inconclusive; ev-A, ev-B not-found; controls 26/26.

Disposition: lead's ruling.

## 6. Amendment A1 (lead approval 2026-10-08), committed before the A1 code

| # | Finding | Rule | Probe (pinned test) |
|---|---|---|---|
| G1' | 1 | Replaces G1's test. An anchor is capped (for a signature with no topology label) when its declared `topology_label` contains a TOUR label ({"mga-tour", "pump-tour", "ephemeris"}) and NO cycler-class label. CYCLER_CLASSES = {"repeated-moon", "halo", "nrho", "tulip", "binary-coorbital", "quasi-satellite", "retrograde-satellite", "axisymmetric", "planar"}. "resonant" alone is not a cycler class. | Jupiter (Ganymede, Callisto, Ganymede), untopologied, period 3: not "published", and the citation is not Niehoff. The Sun E-V probes of G1 still hold. |
| A1-2 | 2 | When the best hit came from a capped anchor and the result is "inconclusive", the notes also name the F7 anchors (`_different_architecture_same_system`), when there are any. | Sun (E, V), untopologied and unlabelled: the inconclusive notes name the H&M anchor. |
| A1-3 | 3 | The `_declared_scope_exclusions` docstring's first line is corrected. | — |
| A1-4 | 4 | The F13 ratchet (in the test) also flags: (a) a comparison of a name assigned, in the same function or module, from a status access; (b) a comparison of a status access with an Attribute (e.g. `mod.OTHER_SET`, `consts.NOT_FOUND`). | Synthetic snippets: `st = r.status` then `st in FRESH`; `r.status in mod.OTHER_SET`; `r.status == consts.NOT_FOUND`: all flagged. The current src/ and scripts/ tree: clean. |
| A1-5 | 5 | The H&M anchor comment says Menning p.42 names HR variations but Table 3 computes none ("anticipated, not tabulated"). | — |
| A1-6 | 6 | F7's inclusion test requires a non-empty sequence (`seq_set and seq_set <= a.body_set`). | Jupiter, (), working "two": not inconclusive via F7. |

Expected: the 26 controls and the five candidates unchanged. The catalogue probe (v2 against v3-A1)
moves rows only toward "inconclusive". If any candidate or control changes: STOP and report.

The untopologied Pergola (Uranus) and Howett (Pluto) limit is recorded, not fixed here; the lead
registered it as #1052.

## 7. A1 results (2026-10-08), recorded after the run

- Controls and candidates UNCHANGED: 26/26 published; the Jones negative holds; gc-1, gc-2 and ev-C
  inconclusive; ev-A and ev-B not-found. `data/942_943_litcheck_scope.json` is byte-unchanged.
- Catalogue probe, v2 against v3-A1 (790 signatures, `probe_v2_v3a1.out`): 29 change, all
  published -> inconclusive, all in the no-topology half. 0 move toward fresh.
  - The 17 of sec. 4.
  - The 10 R-S Jovian rows (ganio-53/185/403, ganeur-5/43/316, eurgan-131/159, gancal-1/5) and
    gc-1/gc-2. Without a topology label, v2 cited them to tour papers (Niehoff etc.). This is the
    G1' fix.
- Jupiter (Ganymede, Callisto, Ganymede) untopologied: inconclusive 0.69 (the capped
  Strange/Campagnola/Russell tour anchor), not Niehoff. Sun (E, V) untopologied: the inconclusive
  notes name the H&M anchor (A1-2).
- Pinned tests: test_1045 (+3 A1 tests, +3 ratchet forms), test_972 and test_literature_check pass.
  The src/ and scripts/ ratchet is clean.
- Full set (tee'd, 2026-10-08):
  - tests/data, tests/scripts and test_catalogue_rediscovery: EXIT 0 (15:57-16:04 AEDT).
  - tests/search: EXIT 1 (16:04-17:03, load up to 23 with the #1025 batch and CI on the box). The
    only failure was the known load-sensitive Earth-Moon test
    `test_known_close_pair_73c_plateaus_just_outside_guard` (600-s timeout).
  - The lead re-ran it alone at 17:04-17:08 AEDT, load 5: EXIT 0.
- ruff clean; mypy src tests clean (939 files).
- Next: the read-only review of the A1 diff.
