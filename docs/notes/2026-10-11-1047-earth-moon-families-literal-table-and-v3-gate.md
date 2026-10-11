# #1047: Earth-Moon families: literal collision table and the v3 literature gate (2026-10-11)

Owner of this note: earthmoon-opus. No catalogue edits. Nothing in this note is called novel.

## 1. Scope

The cycler-class families of #1000 (data/1000_complement/shooting/grid.jsonl) and #1048
(gaps.jsonl): families A (5:2 pro asym), B (5:3 ret asym) and D (7:3 pro asym), the symmetric
families, the 3:1 retrograde family and the 7:2 prograde family. A and B and D had a literal table
before (f0c24e21, #1000 sec. 5.1b). The others had not been started.

## 2. Literal collision table (done)

Script: `scripts/run_1047_literal.py`. Output: `data/1047_literature/literal_table.{json,md}`.

Every published Earth-Moon p:q orbit that is held as numbers (the rows of
data/1000_complement/literal_check.json: Casoliva 2010 Table 3, Vaquero 2013, Newton 1959, Casoliva
2008 seeds, Liang 2017, Liang 2020, Hoelker-Winston 1968, Genova-Aldrin 2015) is compared with every
family, row by row. Each family's p:q and sense are computed again from its converged orbits
(N = w_E + T/2pi). The match rule is the #1000 sec. 5.1b rule: same p:q and sense, published C inside
the family's C range +- 0.01, and at that C |dT|/T < 1e-3 with perigee and periselene altitudes
within 10%.

Results:

| family | literal members | nearest published rows |
|---|---|---|
| 7:3 ret wE-10 wM-1 asym | Casoliva 2010 7-3b, 7-3c | (positive check: these are the #1053 rows) |
| 7:3 ret wE-10 wM-1 sym | Casoliva 2010 7-3a | 7-3b same p:q, dT/T 1.2e-3 |
| 2:1 pro wE1 wM-1 sym | none | Casoliva 2-1a/b opposite sense; Vaquero and Newton F3 rows outside the C range; known class through #997 F3 and Ross & Roberts-Tsoukkas |
| 3:1 pro wE2 wM0 sym (#1048) | none | Vaquero c254, c313: C outside the computed range; known class through Ross & Roberts-Tsoukkas (#1048) |
| 3:1 ret wE-4 wM0 sym (#1048) | none | Vaquero rows have the opposite sense |
| 5:2 pro sym, A (5:2 pro asym), 5:2 ret (two families) | none | Casoliva 2008 52a (mu = 1e-6, class level only); Liang 2017 5:2 PLPO (no literal match, #1050) |
| 2:3 ret | none | Casoliva 2008 23a/b (mu = 1e-6, class level only) |
| 7:3 pro sym, D (7:3 pro asym) | none | Casoliva 2010 7:3 rows have the opposite sense |
| 4:3 ret, B (5:3 ret), 7:2 pro (#1048), 8:3 ret | none | no published row with the same p:q |
| 7:3 ret wE-10 wM0 sym (1 orbit) | none | Casoliva 7-3a/b/c C outside the range |
| "7:3 ret wE-9 wM0 sym" (1 orbit) | none | its computed p:q is 6:3, not its 7:3 skeleton key (recorded, not a family) |

## 3. v3 literature gate: the negative control FAILED, run stopped

Script: `scripts/run_1047_gate.py`. Output: `data/1047_literature/gate_results.json`.

Signature: primary "Earth", sequence ("Moon",), resonances ("p:q",), period_k = q, period band from
the computed members, topology {"resonant"} (sensitivity: {"repeated-moon"}), working_bodies "one".

| control | required | result |
|---|---|---|
| positive: Casoliva 7-3b (7:3, T 18.8496), {"resonant"} | published | published, Braik & Ross 2026, 0.95 (pass, but it cites Braik & Ross, not Casoliva) |
| positive, {"repeated-moon"} | published | published, Braik & Ross 2026, 0.95 |
| negative: same signature, primary "Jupiter" | not published | known-architecture-new-system (pass) |
| negative: same signature, topology {"pump-tour"} | not published | **published, Braik & Ross 2026, 0.95 (FAIL)** |

Cause: nine Earth-Moon anchors in the candidate footprint carry no topology label: Braik-Ross 2026,
the four Kumar/Rawat/Rosengren/Ross papers, Casoliva et al. 2010, Rosengren et al. xGEO 2024,
Hiraiwa et al. 2026 and Sanaga-Park-Howell 2026. The topology filter
(`_declared_scope_exclusion`) excludes an anchor only when both sides declare a label, so an
unlabelled anchor matches any topology. For an Earth primary with the Moon as the only body, the
candidate's topology label has no effect, and Braik & Ross (first confident hit) is returned for every
signature.

Procedural record: the first run used a looser negative test (the tour-label control must not cite
Casoliva). That test passed, so the first run also scored all 19 families: every one returned
"published", Braik & Ross, 0.95. Those family results are not interpretable and are not kept. The
script now uses the lead's criterion (a mislabelled signature must not return "published") and stops
after the controls.

## 4. Proposed #875 classification (provisional; the gate step is not done)

From the literal table only:

- known-reproduction: 7:3 ret wE-10 wM-1 asym (contains Casoliva 7-3b, 7-3c) and 7:3 ret wE-10 wM-1
  sym (contains Casoliva 7-3a).
- known-class-member: 2:1 pro sym (#997 F3, Ross & Roberts-Tsoukkas), 3:1 pro sym (Ross &
  Roberts-Tsoukkas, #1048).
- undecided, waiting for a gate that discriminates: all other families, including A, B and D. The
  Earth-Moon resonant-orbit literature (Casoliva 2008/2010, Kumar/Rawat/Rosengren/Ross, Antoniadou &
  Libert) covers p:q resonant periodic orbits as a class, so #875 (iii) "known-class-member" is the
  likely outcome for most of them; a family is candidate-novel only if a labelled anchor set shows that
  no published work covers its p:q and sense.

## 5. What is needed to finish (ruling requested)

(a) Give topology labels to the nine unlabelled Earth-Moon anchors, from each paper's actual scope
(a #1045 gate change in literature_check.py, outside my ownership), then rerun the controls and the
families; or (b) accept the literal table plus a manual class-level reading of the anchors as the
literature step for these families.

## 6. Ruling (a) and pre-registration of the anchor labels (written before any rerun)

Lead ruling 2026-10-11: option (a). I label the nine Earth-Moon anchors (the Earth-Moon part of
#1052); I change the anchor entries in `literature_check.py` only.

Rule for a label: the trajectory types the paper itself computes or designs, in the gate vocabulary.
"resonant" = p:q Earth-resonant periodic orbits; "repeated-moon" = periodic orbits with a lunar
encounter each period (the cycler class used for the Ross & Roberts-Tsoukkas anchors); "transfer" =
a new anchor-only label for transfers or connections between orbits (no periodic family); a paper
that also computes resonant orbits keeps "resonant" next to "transfer". "model-transition" = a new
anchor-only label for a paper about moving trajectories between dynamical models.

| anchor | label | grounding |
|---|---|---|
| Braik & Ross 2026, orbital networks | {"transfer", "resonant", "repeated-moon"} | mining note 2026-06-13-braik-ross-2026-orbital-networks-mining.md: Secs. 2-7 build a transfer-accessibility network; Table 2 holds 13 representative planar orbits at C_J 3.1294: Lyapunov, cyclers C11a/C11b/C21/C32, resonant R21/R31/R52 (S and U), DPO. Q1: no new cycler tuples |
| Kumar, Rawat, Rosengren & Ross 2026 (ASR), resonant transport | {"resonant", "transfer"} | digest 2026-06-20-digest-kumar-2025.md: 4:1, 3:1, 2:1 unstable resonant orbits (Table 6, Figs. 7-12); 3:1-2:1 heteroclinics (Table 1) |
| Kumar et al. 2024 (IAC), interior MMR | {"resonant", "transfer"} | digest 2026-07-15-kumar-2024-interior-mmr-heteroclinic-digest.md: 4:1, 3:1, 2:1 interior resonant families; heteroclinic connections |
| Rawat et al. 2025 (AAS 25-569), exterior MMR | {"resonant", "transfer"} | digest 2026-07-15-rawat-2025-exterior-mmr-heteroclinic-digest.md: 1:3, 1:4 exterior resonant orbits; 1:3 to 2:1/3:1 heteroclinics through L1-L2 tubes |
| Rawat et al. 2024 (AAS 24-368), resonance widths | {"resonant", "transfer"} | digest 2026-07-15-rawat-2024-resonance-widths-chaotic-zones-digest.md: 2:1, 3:1 stable and unstable resonant orbits and widths; 3:1-2:1 heteroclinic transfer times |
| Casoliva et al. 2010 (JGCD), two classes | {"resonant", "repeated-moon"} | digest 2026-07-27-725-casoliva-earth-moon-cycler-families-digest.md: Class 1 p-q resonant cyclers (sec. IV), Class 2 L1-Lyapunov homoclinic cyclers (sec. V); both periodic with a lunar encounter each period |
| Rosengren et al. 2024 (AMOS), xGEO | {"resonant"} | digest 2026-07-15-rosengren-2024-xgeo-resonant-structure-digest.md: MMR widths (sec. 2.1.4), Kozai-Lidov (2.1.3), TLE case studies (sec. 4); no transfer design |
| Hiraiwa et al. 2026, lobe dynamics | {"transfer", "resonant"} | mining note 2026-06-07-hiraiwa-lobe-dynamics-method-mining.md: low-energy transfers by lobe sequences (Ch. 3-4, Table 1); the unstable 7:2 and 3:1 resonant orbits are computed as manifold sources (Ch. 3.2) |
| Sanaga, Park & Howell 2026 (JAS) | {"model-transition"} | no digest and no copy in the corpus; Crossref record of DOI 10.1007/s40295-026-00571-5 gives the title "Systematic Numerical Transitions of Cislunar Trajectories Across Dynamical Models" (no abstract). Title-level grounding only; flagged for a digest |

Working-body scope: left unset (None) for all nine. Each paper has one secondary (the Moon), so the
one/two working-body distinction (#942/#943, defined for two flyby bodies) does not apply; setting
"one" would only add a "working-bodies-unlabelled" exclusion for unlabelled Earth-Moon signatures.

Expected outcomes (written before the run):
- Casoliva 7-3b positive ({"resonant"}): "published", citing Casoliva 2010 (not Braik & Ross). Risk
  recorded: Braik & Ross now also carries "resonant" (its Table 2 orbits), and the citation is the
  first corroborated anchor in corpus order, so it may still cite Braik & Ross; if it does, that is a
  failed expectation and I report it, not change a label to fix it.
- Jupiter negative: "known-architecture-new-system" (unchanged).
- pump-tour negative: not "published" (no Earth-Moon anchor carries "pump-tour" or an empty label).
- The gate's 26 controls (scripts/litcheck_942_943_scope.py; none is Earth-Moon): unchanged, all
  "published"; the negative and the five candidates unchanged. Before-run saved.
- tests/search -k literature and tests/scripts: green.
- Families: the gate is class-level (it does not compare p:q), so every family with topology
  {"resonant"} is expected to return "published" (known-class) from a "resonant" Earth-Moon anchor.

## 7. Rerun after the labels (3b434aec): negative controls pass, positive citation expectation FAILS

| check | expected (sec. 6) | result |
|---|---|---|
| Casoliva 7-3b positive, {"resonant"} | published, citing Casoliva | published 0.95, citing **Braik & Ross** (FAIL on the citation) |
| same, {"repeated-moon"} | published | published 0.95, Braik & Ross |
| Jupiter negative | known-architecture-new-system | known-architecture-new-system (pass) |
| pump-tour negative | not published | not-found (pass) |
| 26 controls + negative + 5 candidates (litcheck_942_943_scope.py) | unchanged | output identical to the before-run (diff empty) |
| tests/search -k literature | green | EXIT 0 |
| tests/scripts | green | not completed in my session: two attempts (full, then a 10-file chunk) passed the 7.5-min limit under load (load average ~38, CI runner active). The one literature-related file, test_1038_scripts_import_without_literature_check.py, passes alone. Full run handed to the lead. |

Cause of the citation failure (diagnosed, not fixed): `build_queries` adds a keyword query from the
first in-footprint anchor, and the citation is the first anchor in corpus order that the hit
corroborates. In-footprint order for this signature: Braik-Ross, the four Kumar/Rawat papers,
Casoliva 2010, xGEO, Hiraiwa. Braik & Ross carries "resonant" honestly (its Table 2 holds the 2:1,
3:1 and 5:2 resonant representatives), so it is the first match and it is cited for every Earth-Moon
resonant signature. The citation reflects corpus order, not which paper is closest. The status
("published") does not depend on it. I did not change a label or the order to make the control pass.

Per the pre-registration this is a failed control, so the run stopped before the families
(`data/1047_literature/gate_results.json`: controls_pass false, positive_cites_casoliva false). The
families are not scored.

## 8. Ruling (i), family scores and the proposed #875 classification

Lead ruling 2026-10-11: (i). The positive control **passed at status level** ("published") and
**failed the pre-registered citation** (it cites Braik & Ross, not Casoliva). Cause as in sec. 7:
the gate cites the first in-footprint anchor in corpus order. The citation choice is registered as
#1062 (c93b4f65, backlog for twobody-gen2; test case Casoliva 7-3b). The script now records the
citation (`positive_cites_casoliva: false`) and gates on status only.

Count correction: the table holds **18** entries, not the "19 families" in secs. 3-5 and my
messages. 17 are families. The 18th is the one-orbit "7:3 ret wE-9 wM0" entry, whose computed p:q is
6:3.

All 18 entries return "published" 0.95, citing Braik & Ross, for both {"resonant"} and
{"repeated-moon"} (`data/1047_literature/gate_results.json`). **The gate does not compare p:q,
sense or symmetry.** It matches on primary, body set and topology class only, so "published" here
means: the Earth-Moon resonant / repeated-lunar-encounter class is published (Casoliva 2010, Kumar,
Rawat, Rosengren & Ross, Braik & Ross). It does not mean that this family's own p:q, or a member of it, is
published. That information is in the literal table (sec. 2).

| family (p:q, sense, winding, symmetry) | candidates | literal-table outcome | gate | proposed #875 |
|---|---|---|---|---|
| 7:3 ret wE-10 wM-1, asymmetric | 13 | literal members Casoliva 2010 7-3b, 7-3c | published (class) | known-reproduction |
| 7:3 ret wE-10 wM-1, symmetric | 70 | literal member Casoliva 2010 7-3a | published (class) | known-reproduction |
| 2:1 pro wE1 wM-1, symmetric | 4 | no member; Casoliva 2-1a/b opposite sense; Vaquero/Newton F3 rows outside the C range; class through #997 F3 and Ross & Roberts-Tsoukkas | published (class) | known-class-member |
| 3:1 pro wE2 wM0, symmetric (#1048) | 18 | no member; Vaquero c254/c313 C outside the computed range | published (class) | known-class-member |
| 3:1 ret wE-4 wM0, symmetric (#1048) | 17 | no member; Vaquero rows opposite sense | published (class) | known-class-member |
| 2:3 ret wE-5 wM-2, symmetric | 10 | no member; Casoliva 2008 23a/b at mu 1e-6 (class only) | published (class) | known-class-member |
| 4:3 ret wE-7 wM-2, symmetric | 4 | no member; no published row with this p:q | published (class) | known-class-member |
| 5:2 pro wE3 wM-1, symmetric | 11 | no member; Casoliva 2008 52a (mu 1e-6), Liang 2017 5:2 (no literal match, #1050) | published (class) | known-class-member |
| A: 5:2 pro wE3 wM0, asymmetric | 3 | no member; as above | published (class) | known-class-member |
| 5:2 ret wE-7 wM-1, symmetric | 12 | no member; as above | published (class) | known-class-member |
| 5:2 ret wE-7 wM0, symmetric | 26 | no member; as above | published (class) | known-class-member |
| B: 5:3 ret wE-8 wM-2, asymmetric | 4 | no member; no published row with this p:q | published (class) | known-class-member |
| 7:2 pro wE5 wM0, symmetric (#1048) | 15 | no member; no published row with this p:q | published (class) | known-class-member |
| D: 7:3 pro wE4 wM-1, asymmetric | 2 | no member; Casoliva 7:3 rows opposite sense | published (class) | known-class-member |
| 7:3 pro wE4 wM-1, symmetric | 50 | no member; Casoliva 7:3 rows opposite sense | published (class) | known-class-member |
| 7:3 ret wE-10 wM0, symmetric | 1 | no member; Casoliva 7-3a/b/c C outside the range | published (class) | known-class-member |
| 8:3 ret wE-11 wM0, symmetric | 2 | no member; no published row with this p:q | published (class) | known-class-member |
| (one orbit, computed 6:3 ret wE-9 wM0) | 1 | not a family; no member | published (class) | not classified (one orbit, p:q inconsistent with its skeleton) |

Reading of the known-class-member rows: under #875 the literal check comes first, and none of these
has a literal member. The class (planar Earth-Moon p:q resonant periodic orbits with a lunar pass
each period) is published, so they are known-class-member under #875 (iii). For 4:3 ret, B (5:3
ret), 7:2 pro and 8:3 ret no published orbit with the same p:q is held as numbers. The class covers
them, and the specific p:q is not documented in the held sources. Nothing here is candidate-novel.
No catalogue edits.

Count: 2 known-reproduction + 15 known-class-member + 1 unclassified orbit = 18. (The commit message
of bba3eeb8 says "16 families known-class-member"; the correct number is 15.)

## 9. Owner decision 2026-10-11 and the window patch for rows A, B, D (pre-registration)

Owner decision (verbatim label): "Rows for A, B, D only". Three rows, our_status
known-class-member (#875 (iii)), validation V1. The other 12 symmetric families and the two
7:3 retrograde wM-1 families stay here as census.

Rows (ids):
- A: earth-moon-5-2-prograde-asymmetric-cycler-1000a-2026
- B: earth-moon-5-3-retrograde-asymmetric-cycler-1000b-2026
- D: earth-moon-7-3-prograde-asymmetric-cycler-1000d-2026

Representative: the #1000 grid member that seeded each family's continuation, re-closed with
`correct_periodic`. V1 evidence: `scripts/writeback_1047_rows.py`, `data/1047_rows/writeback_v1.json`,
`tests/search/test_1047_rows_writeback_v1.py` (305adc7c).

| row | C | T | k (Casoliva Eq. 8) | perigee alt (km) | periselene alt (km) | closure | Radau dJ | Radau arcs |
|---|---|---|---|---|---|---|---|---|
| A | 2.4818625406 | 12.7126038450 | -110.10 | 3,332 | 13,243 | 3.7e-11 | 1.9e-12 | 4.6e-11 |
| B | 0.8604292804 | 18.8866755266 | 115.38 | 3,569 | 21,122 | 6.7e-11 | 2.3e-12 | 4.6e-10 |
| D | 2.5472667556 | 18.6647180953 | 717.95 | 5,754 | 13,774 | 5.3e-12 | 2.3e-12 | 1.8e-11 |

One deviation from the Casoliva rows' V1 recipe, stated in each row: `crosscheck_periodic`'s whole-period
Radau closure (tolerance 1e-8) is 1.6e-9 for A but 6.4e-7 for B and 4.2e-8 for D. A plain DOP853
run (without the STM) gives 9.9e-8, 7.1e-6 and 8.4e-7. These orbits have |lambda| of about 110-720,
so a whole-period closure depends on the integrator's step sequence. The independent-integrator check is
therefore made arc by arc (4p arcs, each Radau arc from the DOP853 node; <= 4.6e-10). This is the same
reason #1000 used multiple shooting.

Window patch (apply script in my scratchpad, one match asserted per edit):
1. `data/catalogue.yaml`: the three rows inserted after casoliva-7-3c-em-cycler-2010 (rows generated
   by `scripts/writeback_1047_rows.py --rows-only --rows-out`).
2. `src/cyclerfinder/data/validate.py`: a `_V1_1047` evidence string and three `_LEVEL_EVIDENCE`
   entries (V1).
3. Ratchet deltas, pre-registered:
   - `tests/data/test_cycler_class_census.py`: non-keplerian 56 -> 59 (expected dict and docstring).
   - `tests/data/test_validation_tier_census.py`: unvalidated 144 -> 147 (this-project rows, no
     provenance-registry tags).
   - `tests/data/test_schema_v45_fields.py`: live V1 census + the three ids at V1.
   - `tests/test_catalogue_rediscovery.py`: NON_HELIOCENTRIC 111 -> 114 (Earth primary).
   - No other ratchet is expected to move. Catalogue 423 -> 426 rows; candidate-novel unchanged
     (29); known-class-member +3.
4. `tests/search/test_1047_rows_writeback_v1.py`: its catalogue check stops skipping and must pass.

### 9.1 Dry run (2026-10-11 15:08-16:04 AEDT; scratch worktree at 7feb2bf5, removed)

- Apply script: every anchor matched once (8 edits after the fix below). Patch: 6 files, +487 / -4.
- check-jsonschema on the patched catalogue: ok. ruff check and format on the touched files: clean.
  mypy (the hook's file set, src + tests, worktree src first on the path): no issues in 944 files.
- `tests/search/test_1047_rows_writeback_v1.py` plus `tests/test_catalogue_rediscovery.py`: pass. The
  writeback test's catalogue check runs, with no skip.
- `tests/data` (82 files, in chunks under 8 minutes): first chunk FAILED
  `test_non_keplerian_ids_match_ratchet`. **Pre-registration miss:** the id-level ratchet
  `NON_KEPLERIAN_IDS` (56 -> 59, its own size assertion) also moves with the non-keplerian count; I
  had listed only the count. Added to the apply script (the three ids, size 59); the worktree was
  rebuilt from HEAD and re-patched. Then all 82 files passed (EXIT 0 in every chunk).
- Final ratchet deltas: non-keplerian count 56 -> 59 and NON_KEPLERIAN_IDS +3 (59); unvalidated
  144 -> 147; V1 census +3; NON_HELIOCENTRIC 111 -> 114. Nothing else moved.

## 10. What the three rows mean (for the owner summary)

**What they are.** Rows A, B and D are three asymmetric periodic orbits of the planar Earth-Moon
CR3BP. Each passes the Moon once per period, inside the Hill radius. They are p:q resonant: A
makes 5 inertial revolutions about the Earth while the Moon makes 2, prograde. B is 5:3
retrograde. D is 7:3 prograde. Each row is one member of a family that #1000 and #1048 continued
over a range of Jacobi constant. All three are strongly unstable (|k| about 110, 115 and 720).

**How they relate to Casoliva's class.** Casoliva et al. 2010 describe two classes of Earth-Moon
cyclers. Their Class 1 is the near-Keplerian p-q resonant cycler: a spacecraft on an Earth orbit
in p:q resonance with the Moon, corrected into an exact periodic orbit of the CR3BP. That class
includes asymmetric orbits: their 7-3b and 7-3c are one asymmetric orbit and its mirror (#1053),
and our family C contains them. A, B and D are the same kind of object. Casoliva's Table 3
tabulates only p-q = 1-2, 2-1, 3-2 and 7-3, and her 7-3 rows are retrograde. So A (5:2), B (5:3)
and D (7:3 prograde) are at p:q values or a sense that the source did not tabulate.

**Why known-class-member and not candidate-novel.** Under the #875 policy the literal check comes
first, then the class:
- Literal: no published orbit held as numbers is a member of any of the three families (#1047 literal
  table: Casoliva 2008 and 2010, Vaquero 2013, Newton 1959, Liang 2017 and 2020, Hoelker-Winston
  1968, Genova-Aldrin 2015). So none of them is known-reproduction.
- Class: the object class is published (Casoliva Class 1, including its asymmetric members). A
  p:q resonant Earth-Moon periodic orbit is also the subject of the Kumar, Rawat, Rosengren & Ross
  papers, and of Liang 2017 for the same 5:2 and 7:3 ratios, with different orbits. A new p:q value
  or sense inside a published class is "a computed member of a published class" (spec sec. 16.4;
  #875 (iii)). It is not an author-excluded class (#875 (i)) and not a new system (#875 (ii)).
- The two databases that would hold such orbits leave them out by their own construction rules.
  Restrepo-Russell keeps axisymmetric orbits only. Franz-Russell drops orbits that go further than
  about 350,000 km from the Moon, and these reach 795,000-957,000 km. Their absence from those
  databases is therefore not evidence of novelty.
- The v3 literature gate returns "published" at class level. It does not compare p:q, so it
  supports the class reading and says nothing about the specific families.

**What would change the classification.**
- To known-reproduction: a published orbit that falls on one of these families (a member of A, B or
  D), found in a source not yet held as numbers. Examples are a Liang-type averaged orbit continued
  through A's suspected fold near C 3.097 (#1050), or another author's p:q table.
- To candidate-novel: only if the class reading is wrong. That would need evidence that these
  orbits are not of Casoliva's Class 1 kind, for example a topology the near-Keplerian
  resonant-cycler construction cannot produce. Nothing found so far points that way. Their windings,
  resonance and continuation from p:q Kepler skeletons are exactly that construction.
- Withdrawal: if the V1 evidence fails on a re-check. An example is a closure that does not hold arc
  by arc with an independent integrator, or a periselene outside the Hill radius at the row's
  state. The writeback test re-checks both on every run.
