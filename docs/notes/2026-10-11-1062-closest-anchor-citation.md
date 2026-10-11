# #1062 Literature gate: cite the closest anchor, not the first in corpus order — PRE-REGISTRATION (2026-10-11), committed before any code

Registered by the lead (ruling on #1047). Diagnosis: `docs/notes/2026-10-11-1047-earth-moon-families-literal-table-and-v3-gate.md`
secs. 7-8.
- `build_queries` adds the keyword queries of the in-footprint anchors in `KNOWN_CORPUS` order.
- `check_literature` stops at the first hit at or above `MATCH_THRESHOLD`, and cites the first anchor
  in corpus order that the hit corroborates.
- For every Earth-Moon resonant signature that is Braik & Ross 2026, even for a literal Casoliva 2010
  member (7-3b).
- The status is not affected; the citation is.

## 1. The rule (one change, ranking only)

1. A new optional `CorpusAnchor` field, `resonances_scope: frozenset[str] | None`: the p:q resonances
   the SOURCE computes, in the signature's "p:q" convention (q T_moon = p T_spacecraft). Declared only
   where a held digest prints them:
   - Casoliva et al. 2010 (JGCD 33(5)), Table 3: {"1:2", "2:1", "3:2", "7:3"}. Digest
     `docs/notes/2026-07-27-725-casoliva-earth-moon-cycler-families-digest.md` lines 78-88: 1-2c/d/e,
     2-1a/b, 3-2c, 7-3a/b/c.
   - Braik & Ross 2026, Table 2: {"2:1", "3:1", "5:2"}. Digest
     `docs/notes/2026-06-13-braik-ross-2026-orbital-networks-mining.md` line 26: R21-S/U, R31-S/U,
     R52-S/U. Its cyclers C11a/C11b/C21/C32 are (k1, k2) winding labels, not p:q resonances, and are
     not included.
   - Every other anchor: undeclared (None).
2. CLOSENESS RANK of an in-footprint anchor for a signature that carries `resonances`:
   - 0: the anchor declares `resonances_scope` and it intersects the signature's resonances;
   - 1: the anchor does not declare it;
   - 2: the anchor declares it and it is disjoint.

   Ties keep `KNOWN_CORPUS` order. A signature without `resonances` ranks every anchor 1, so the
   order is unchanged.
3. The rank orders (a) `build_queries` step 4, the anchor keyword queries, and (b) the anchors tried
   for the published citation in `check_literature`.
4. `resonances_scope` is NOT an exclusion. It never removes an anchor and never changes a confidence
   score, so the status cannot change through it.

## 2. Probes (pinned tests under tests/search/, each < 2 s)

| Signature (Earth, ("Moon",), topology {"resonant"}, working "one") | Expected citation |
|---|---|
| 7:3, period_k 3 (the #1047 Casoliva 7-3b control) | Casoliva 2010 (rank 0; Braik rank 2) |
| 3:1 | Braik & Ross (the only rank 0) |
| 2:1 | Braik & Ross (both rank 0; corpus order) |
| 4:3 (no anchor declares it) | the first UNDECLARED in-footprint anchor in corpus order (rank 1): Kumar-Rawat-Rosengren-Ross cislunar resonant transport (2026) |
| no resonances | unchanged from v3-A1 (Braik & Ross) |

The status is "published" in every row, unchanged.

## 3. Controls and expected outcomes (written before running)

- `scripts/litcheck_942_943_scope.py`: the 26 controls, the negative and the five candidates.
  `data/942_943_litcheck_scope.json` stays byte-identical, because none of these signatures carries
  `resonances`.
- `scripts/run_1047_gate.py`: all 18 Earth-Moon entries keep status "published". Citations become:
  - 7:3 entries (the three 7:3 ret families, D, 7:3 pro symmetric) -> Casoliva 2010;
  - 2:1, 3:1 and 5:2 entries -> Braik & Ross;
  - 2:3, 4:3, 5:3, 7:2, 8:3 and the one 6:3 orbit -> Kumar-Rawat-Rosengren-Ross (2026), the first
    undeclared anchor.
  - The positive control `positive_cites_casoliva` becomes true. The negative controls are unchanged.
- The literature test files (`tests/search/test_literature_check.py`, `test_972`, `test_1045`) stay
  green.
- If any STATUS changes anywhere: STOP and report.

## 4. Process

- Implement; add the pinned tests; run the two scripts and the literature tests; then ruff and full
  mypy; commit.
- The full suite runs only in a lead window, if the lead wants one. This is a code-only change, with
  no catalogue or census effect.

## 5. Results (2026-10-11 14:12 AEDT)

- `scripts/litcheck_942_943_scope.py`: `data/942_943_litcheck_scope.json` is BYTE-IDENTICAL (cmp).
- `scripts/run_1047_gate.py` (`data/1047_literature/gate_results.json` updated):
  - `positive_cites_casoliva` false -> TRUE; `controls_pass` true (unchanged);
  - both negatives unchanged;
  - all 18 entries keep status "published" under both topologies (asserted).
  - Citations under {"resonant"}, as pre-registered:
    - 7:3 entries -> Casoliva 2010;
    - 2:1, 3:1, 5:2 -> Braik & Ross;
    - 2:3, 4:3, 5:3, 7:2, 8:3 and the 6:3 orbit -> Kumar-Rawat-Rosengren-Ross (2026).
  - Under {"repeated-moon"} the undeclared group cites Roberts-Tsoukkas & Ross instead. That anchor is
    in the footprint only for that topology, and it precedes Kumar there in corpus order: the same
    rule (first undeclared anchor). Sec. 3 did not state the repeated-moon variant.
- Pinned tests: `tests/search/test_1062_closest_anchor.py`, 5 cases, pass.
- The literature tests (`test_literature_check`, `test_972`, `test_1045`, `test_1062`) and
  `tests/scripts/test_1038_scripts_import_without_literature_check.py`: EXIT 0.
- ruff clean; mypy src tests clean (946 files).
- No status changed anywhere.
