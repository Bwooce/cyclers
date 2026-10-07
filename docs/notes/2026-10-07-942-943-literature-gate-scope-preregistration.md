# #942 / #943 literature_check gate: declared-scope tags and a Hollister/H&M anchor — PRE-REGISTRATION (2026-10-07)

Committed BEFORE any change to `src/cyclerfinder/search/literature_check.py`. Lead ruling 2026-10-07:
option (a), with the conditions below. Whatever the gate returns after the change is the result. The tags
are not adjusted after the candidate outputs are seen. If a positive control fails, the failing tag is fixed
against its SOURCE, and ALL controls are re-run.

## 0. The problem (diagnosed 2026-10-07, before this note)

`check_literature` with the offline backend returns "published" for every #942/#943 candidate, and not
because of prior art:
- gc-1 and gc-2 (Ganymede-Callisto): confidence 0.95 against Liang et al. 2024 Callisto-Ganymede-Europa
  TRIPLE cyclers. The candidate's two-moon set is a subset of the anchor's three moons.
- ev-A, ev-B and ev-C (Earth-Venus): confidence 0.85 against Jones et al. 2017 Venus-Earth-Mars TRIPLE
  cyclers (and Hughes et al. 2014), again by body subset.
- KNOWN_CORPUS has no Hollister 1969 / Hollister & Menning 1970 Earth-Venus anchor at all. The actual
  published E-V cycler class is missing.

## 1. New declared-scope fields (#880 mechanism: a filter applies only if BOTH sides declare)

On `CorpusAnchor` (default None = undeclared = no filter):
- `n_bodies_scope: int | None`: the number of distinct encounter bodies in the anchor's family.
- `working_bodies_scope: "one" | "two" | None`: whether the anchor's ideal model has one working (flyby)
  body with a massless target ("one"), or every encounter body bends ("two").
- `return_types_scope: frozenset[str] | None`: the same-body return types the anchor's authors COMPUTED.
- `alternating_scope: bool | None`: True if the anchor's family visits the bodies strictly alternately
  (no two consecutive encounters at the same body).

On `CandidateSignature` (default None = undeclared), plus one quantity derived from existing data:
- n_bodies: derived, `len(set(sig.sequence))`. Always declared.
- `working_bodies: "one" | "two" | None`.
- `return_types: frozenset[str] | None`.
- has a consecutive same-body encounter: derived from `sig.sequence` read cyclically. Always declared
  when the sequence lists encounters in order.

Exclusion rules, added to `_declared_scope_exclusion`, each returning a reason string:
- "n-bodies": anchor declares `n_bodies_scope` and it differs from `len(set(sig.sequence))`.
- "working-bodies": both declare, and the values differ.
- "return-types": both declare, and `sig.return_types` is not a subset of `anchor.return_types_scope`.
- "alternating": anchor declares `alternating_scope=True`, and `sig.sequence` (cyclic) has two
  consecutive equal bodies.

## 2. Mechanical derivation of a candidate's labels (no per-row hand setting)

From the cycle key and the ideal-model flyby table (`scripts/gauntlet_942.py`, the same code for every
candidate):
- `working_bodies` = "two" if every body in the key has a demanded turn > 1e-6 deg at some encounter,
  else "one".
- `return_types`: one entry per same-body leg in the key:
  - `R<b>/1:1` -> "FR" (Menning's full-revolution return: one planet period);
  - `R<b>/n:m` other -> "FR-n:m";
  - `H<b>/...` -> "HR" (half-revolution);
  - `L<b>><b>/1x` with flight time / body period in (1, 2) -> "SY" (Menning's symmetric return:
    "time of flight between one and two planet periods", one extra revolution; Menning 1968 ch. 2-3,
    digest `2026-10-05-digest-menning-1968-mit-thesis-earth-venus-periodic-orbits.md` sec. 1);
  - any other same-body Lambert leg -> "GEN".
- `sequence` = the bodies of the cycle's flybys in time order, as the gauntlet already passes.

## 3. Tags on anchors, each with its source

| Anchor | Tag | Source (quoted or read) |
|---|---|---|
| Liang et al. 2024 CGE triple cyclers | n_bodies_scope = 3 | title "Callisto-Ganymede-Europa Triple Cyclers" |
| Hernandez/Jones/Jesick IEG triple cyclers | n_bodies_scope = 3 | title "One Class of Io-Europa-Ganymede Triple Cyclers" |
| Jones et al. 2017 VEM triple cyclers | n_bodies_scope = 3 | title "Low Excess Speed Triple Cyclers of Venus, Earth, and Mars" |
| Hughes-Edelman-Longuski VEM extensions (2014) | n_bodies_scope = 3 | the anchor's own body set {V, E, M}; the paper's E-V-M sequences (digest `2026-10-06-digest-hughes-et-al-2015-fast-free-returns-mars-venus.md`) |
| Russell-Strange 2009 Ganymede-Io, Ganymede-Europa, Ganymede-Callisto; R-S 2007 Venus-Mars, Venus-Mercury | working_bodies_scope = "one" | AAS 07-118 p.2: "Free-return cyclers ... where one of the two orbiting celestial bodies in the ideal model is considered massless"; p.18 names removing it as future work (digest `2026-10-05-digest-russell-strange-2007-aas-07-118-planetary-moon-cyclers.md`) |
| Campagnola et al. 2019 GCGC | working_bodies_scope = "two"; alternating_scope = True | Fig. 11: flybys "G1, C2, G3, C4, G5"; "repeating two GCGC cycles (G1-G5 and G5-G9)"; both moons host flybys (digest `2026-10-05-digest-campagnola-2019-europa-clipper-tour-design-techniques.md` sec. 2) |
| NEW: Hollister 1969 / Hollister & Menning 1970 (Menning 1968) Earth-Venus periodic swing-by orbits | primary Sun, body_set {E, V}, topology_label {"repeated-moon"}, n_bodies_scope = 2, working_bodies_scope = "two", return_types_scope = {"FR", "SY"} | Menning 1968 App. A-1: every computed orbit is built from FR, SY, FRSY and TFR returns; flybys at both planets (p.3-6); Hollister 1969 orbits I-III; H&M 1970 Table 3 |

The Hollister/H&M anchor is added regardless of the outcome (lead ruling: a real gap). Its provenance is
"verified-against-source" (digests of Hollister 1969, Menning 1968 and the H&M 1970 Table 3 recheck).

## 4. One new return value (lead ruling: the gate must express "known architecture at a new system")

- New status: `"known-architecture-new-system"`.
- It is returned when:
  - the literal check would return "not-found"; AND
  - the signature declares `working_bodies`; AND
  - some anchor with the SAME primary and the same declared `working_bodies_scope` has a body set that
    does not contain the signature's bodies (a different system).
- The citation is that anchor's (all such anchors are listed in `notes`).
- Precedence: "published" (literal) > "inconclusive" > "known-architecture-new-system" > "not-found".
  An inconclusive literal result is left for a human.
- `is_novelty_claimable`: "known-architecture-new-system" is claimable. That is the owner's `#875` (ii);
  the attribution obligation is in the row's notes.

## 5. Positive controls (judged by the same gate; each must return "published" after the change)

Signatures built exactly as the gauntlet builds them (ideal-model zeros from the cell runs), with the
derived labels:
- GanCal#5 and GanCal#1 (gc cell, Callisto turn 0) -> R-S Ganymede-Callisto.
- GanEur#43 and EurGan#131 (ge cell) -> R-S Ganymede-Europa.
- GanEur#316 (ge cell, `data/943_ganeur316_recall.json`) -> R-S Ganymede-Europa.
- VenMar#45 (vm2 cell, Mars turn 0) -> R-S Venus-Mars.
- Hollister 1H (k2|RE/1:1|LE>V/0s|RV/1:1|RV/1:1|LV>E/0s) and 2H (k2|RE/1:1|LE>V/0s|LV>V/1l|RV/1:1|LV>E/0s;
  its V-V leg is 1.46-1.47 Venus periods = SY) (ev cell) -> the new Hollister/H&M anchor.
- The catalogued H&M orbits (hollister-menning-1970-ev-orbit-01..15), signatures from catalogue rows
  (labels undeclared) -> Hollister/H&M.
- R-O 2.5.1.+0 (em recall, `data/942_em_ro_recall.json`) -> R-O / McConaghy Earth-Mars.
- The existing test suites: `tests/data/test_citation_integrity.py` and every literature_check test in
  tests/search and tests/data, unchanged.

## 6. Expected results for the candidates (written before running)

| Candidate | Labels (derived) | Expected status | Why |
|---|---|---|---|
| gc-1 | 2 bodies, two, {SY, FR}, consecutive same-body | not-found | Liang excluded (n-bodies); R-S GanCal excluded (working-bodies); GCGC excluded (alternating) |
| gc-2 | 2 bodies, two, {SY} (G-G 10.85 d = 1.52 P_G; C-C 21.51 d = 1.29 P_C), consecutive same-body | not-found | as gc-1. The owner's known-class ruling stands on the human analysis (6.20-6.21) |
| ev-A | 2 bodies, two, {FR, GEN}, consecutive same-body | not-found | Jones/Hughes excluded (n-bodies); Hollister/H&M excluded (return-types: GEN) |
| ev-B | 2 bodies, two, {FR, GEN, FR-3:2}, consecutive same-body | not-found | as ev-A |
| ev-C | 2 bodies, one, {GEN} | known-architecture-new-system (R-S 2007 Venus-Mars / Venus-Mercury) | Hollister/H&M excluded (working-bodies); triples excluded; R-S heliocentric one-working-body anchors are at other pairs |

(gc-1's G-G leg is 10.58 d = 1.48 Ganymede periods, so SY. The expected statuses do not depend on SY
against GEN for gc-1 and gc-2; the labels are whatever the derivation gives.)

## 7. Process

- Implement, run the controls and the candidates (`scripts/litcheck_942_943_scope.py`, output to
  `data/942_943_litcheck_scope.json`), then all ratchets (one full suite, tee'd), ruff and full mypy.
- A fresh adversarial reviewer (Fable; Opus if unavailable) reviews the diff and this note before the
  catalogue writeback.
- The writeback rows record the gate result as returned.

## 8. Amendment A1 (2026-10-07, after the first control run; the controls failed, so this is allowed by sec. 0)

- First run: every control returned "published" EXCEPT EurGan#131, which returned
  "known-architecture-new-system". The cause is in the label derivation, not a tag:
  - its Ganymede turn in the ge cell is a round-off 1.2e-6 deg, above the 1e-6-deg threshold of sec. 2;
  - ev-C's Earth turn (1.5e-6 deg) is mislabelled the same way.
- Fix: a body "works" if its demanded turn is >= 0.05 deg at some encounter. That is half the 0.1-deg
  turn resolution with which the results note's sec. 6.1 dedupe rule (and every "turn 0" classification
  in secs. 6.3-6.37) defines a zero turn.
- No tag was changed. ALL controls and the candidates are re-run with this derivation.
- The first-run candidate outputs (recorded in the scratch log for transparency):
  - gc-1 and gc-2: not-found;
  - ev-A, ev-B, ev-C: "inconclusive" at 0.475, the best hit being the Aldrin Earth-Mars anchor's
    synthetic title ("cycler" plus "Earth" named = 0.30 + 0.175, above the 0.45 inconclusive floor).
  - Nothing is adjusted in response to them.
