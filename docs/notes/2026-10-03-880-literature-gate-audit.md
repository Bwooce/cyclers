# #880 — offline literature gate: defect, fix, and audit of past verdicts

**Date:** 2026-10-03. **Fix commit:** `d8378af8` (CI green). **Found by:** `#865`, while populating
the spec §16.5 `literature_check` block on the six Uranian quasi_cycler rows.

## The defect

Every offline backend synthesised a search hit from any corpus anchor sharing two body names with
the query and titled it `"<anchor name> (<bodies> cycler)"`. `check_literature` then scored that
text with no reference to the anchor's declared scope:

| term | score |
|---|---|
| the word "cycler" (appended by the backend to every anchor) | 0.30 |
| both moons named | 0.35 |
| primary named | 0.10 |
| **total** | **0.75**, against `MATCH_THRESHOLD = 0.70` |

The `#349` topology-label filter existed only in `_candidate_anchors`. An anchor it correctly
excluded was therefore resurrected through its own synthetic hit. Measured before the fix, with
`topology_label={"repeated-moon"}`: all six catalogued Uranian rows and probe candidates at Uranus,
Saturn, Jupiter, Neptune and Pluto read `published`. At Uranus the citing anchor was
Heaton-Longuski 2003, which the corpus tags `mga-tour` and describes as "NOT a periodic cycler",
with zero structurally overlapping anchors. The offline gate could not return `not-found` for any
moon-system candidate.

## The fix

- `SearchResult.anchor_name`: a corpus-derived hit names the anchor it was synthesised from.
- `_declared_scope_exclusion(sig, anchor)`: the three filters under which BOTH sides declare a
  scope and the two are incompatible (period band `#301`, topology label `#349`, 3D topology
  `#434`), factored out of `_candidate_anchors` with no change of behaviour there.
- `check_literature` does not score a corpus-derived hit whose anchor is excluded that way. Such
  hits still count as "the search returned results" and are listed in the result's `notes`.
- One canonical `offline_corpus_search` in `literature_check.py`; the campaign module and the
  review-queue script import it.

Deliberately unchanged: web hits (no anchor identity), signatures that declare no scope, anchors
that declare no scope, and partial body overlap. The gate stays conservative wherever nothing is
declared.

Measured after the fix (repeated-moon signatures, canonical backend):

| candidate | before | after |
|---|---|---|
| six Uranian catalogue rows | published (Heaton-Longuski tour) | not-found |
| Uranus Miranda-Ariel | published | not-found |
| Saturn Titan-Rhea, Dione-Rhea | published (Cassini tour) | not-found |
| Saturn Titan-Enceladus | published | published (Russell-Strange 2009) |
| Jupiter Europa-Callisto, Ganymede-Io | published | published |
| Neptune Triton-Proteus | published | published (unlabelled anchor, see below) |
| Pluto Charon-Nix | published | published (unlabelled anchors, see below) |
| Uranian row, NO declared topology | published | published (conservative behaviour kept) |

## Audit: which past verdicts consumed the artefact

Method: every tracked data file with a recorded literature `published`/`match` verdict, plus every
call site of `check_literature`. The fix itself only changes outcomes for signatures that declare a
scope; the audit also covers unlabelled signatures, where the same mechanism operates and is left
in place on purpose.

**1. Campaign scorer (`saturn_uranus_campaign.score_candidate`, demotes SILVER to BRONZE on
`published`).** Only the scans that call `run_prioritized_scan` use it: `scan_285` (Saturn, Uranus),
`scan_311` (Saturn pairs), `scan_312` (Uranus pairs), `scan_492` (Neptune, Pluto). `scan_558`,
`enumerate_600` and `scan_816` import only the DOP853 cross-check and are NOT affected (the `#880`
ledger bullet first listed them; corrected).

- Gate-passing records across all of those outputs: one SILVER (`#312` Umbriel-Oberon, `not-found`,
  run before the tour anchors existed) and nine BRONZE in `data/scan_492_pluto.jsonl`.
- **`#492` Pluto: all nine closures were demoted with the false hit as their ONLY recorded reason**
  (`verdict_reasons` = "literature_check published: Howard, Stern et al., Persephone ..."; the ML
  flagger score was 0.58, under its 0.75 limit). Persephone is a mission concept for a Pluto
  orbiter, not a cycler.
- **No verdict flips.** The `#492` note disqualifies the nine on three grounds, and the first two
  are independent of literature: the small-moon-flyby closures are physically infeasible (Nix and
  Hydra cannot supply a gravity assist), and the Charon-flyby closures are model-invalid (Pluto-
  Charon is a binary; the point-primary assumption fails). Ground 3 ("the Pluto-system tour regime
  is published, not virgin ground") is the artefact and is withdrawn; an addendum is on that note.
- Neptune had zero gate-passing closures, so literature was never consulted.

**2. `#641` Sun-Jupiter seed census.** All five clusters read `published` against "Strange-Russell
Tisserand pump-tour graph (2007)" through that script's own backend, which appends "(cycler
trajectory)" to its hits. Same mechanism, mostly on unlabelled (planar) signatures. A pump-tour
graph paper says nothing about Lyapunov or DRO-type periodic orbits. `#641`'s working conclusion
("nothing genuinely novel", classical Sun-Jupiter CR3BP families, already weakened once by `#647`)
does not rest on that citation, but the citation must not be quoted as support for it.

**3. Scope-consistent, not artefacts.** `#468` powered moon tours (signature labelled
`pump-tour`; tour anchors are the right match). `#299`/`#301` Earth-Moon 3D families (265 + 149
records matched to the Antoniadou and Braik-Ross family anchors, separately adjudicated by `#579`).
`#302`/`#307`/`#430` precursor lanes (heliocentric precursors of Aldrin/S1L1, matched to the cycler
they insert into, as intended). These were not re-run.

**4. No recorded verdicts found.** The Titan corridor scripts (`#627`/`#629`/`#633`) call the gate
with a one-body signature, which the two-body route cannot trigger; no `published` status appears
in their tracked outputs.

## Left open

1. **Unlabelled anchors keep Neptune and Pluto at `published`.** "Voyager 2 Triton encounter +
   Trident / Triton-Hopper concept" and the Pluto mission anchors (Persephone and two others) carry
   no topology label. They are flyby and mission references. Until they are labelled, `#868`'s
   literature gate cannot clear.
2. **The campaign scorer still demotes**, because it builds an unlabelled signature and that
   behaviour is deliberately unchanged. `#870`'s charter must declare the genome's topology label.
3. **Historical script copies** of the backend (`run_627`, `run_629`, `run_633`) and the bespoke
   backends in `run_641`, `run_436`, `campaign_468`, `run_299`, `run_301` are untouched run records.
   None sets `anchor_name`, so a re-run would reproduce the old behaviour.
4. 19 of the 69 corpus anchors declare no topology label at all. A labelling pass over them is the
   real completion of `#349`; it is corpus curation, one reviewed decision per anchor.

## Addendum, same day (`#881`)

- Item 1 above is done for Neptune: the flyby anchor is scoped `mga-tour` and the families that
  are published at Neptune-Triton (Miceli & Bosanac 2026; Spear 2021) are now anchors, so a
  `resonant` candidate there reads `published` and a `repeated-moon` one reads `not-found`.
- The "Howard, Stern et al." Persephone citation quoted above is the string recorded in the run
  data and is wrong. CrossRef: Howett, Robbins, Holler et al., *Planetary Science Journal*
  2(2):75 (2021), DOI 10.3847/PSJ/abe6aa. The DOI the anchor carried (10.3847/PSJ/abf837) belongs
  to an unrelated paper. The anchor is corrected; its scope stays undeclared because the paper is
  not in the corpus and has not been read, so Pluto-system candidates still read `published`.
- 16 anchors remain without a topology label (19 before).
