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
