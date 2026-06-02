# Outstanding questions — cyclerfinder catalogue

Long-form log of research questions, source-access gaps, parameter
contradictions, and out-of-paradigm flags encountered while compiling
`catalogue.yaml`. The YAML's per-entry `notes:` field carries
short-form caveats; this file carries the discussion threads that
don't fit there.

**Resolution policy:** when a question is resolved, prefix its heading
with `✓ Resolved (YYYY-MM-DD)` and add a one-line pointer to the
resolution (commit SHA, spec section, errata-investigation note, etc.).
Do not delete the original question text — the audit trail matters.

---

# Project state at a glance (updated 2026-06-03)

This top section is the orientation map for a contributor returning to
the project: what's done, what's in progress (with plan-file pointers),
what's blocked and why, and the prioritised human-actionable items. The
lettered Q&A log (A–H) below it is the per-entry catalogue-sourcing
audit trail and is unchanged in spirit.

## Done

- **Aldrin cycler replicated on the real ephemeris** — M6b binding gate;
  patched-conic constructor reproduces the literature ellipse to ±0.002
  (see Q&A item A, commit `ba12554`).
- **Canonical 2-synodic Earth–Mars anchor replicated idealised**
  (circular-coplanar model).
- **Multi-revolution Lambert solver landed** — universal-variable
  multi-rev solver with sourced cross-check, single-rev Newton restored,
  multi-rev threaded through `verify` (commits `7d620e8`, `9f618f3`).
- **S1L1 topology corrected** — the catalogue no longer carries a
  spurious direct Mars→Earth "return leg" (commit `f4074d0`). See
  Key findings below and the S1L1-nomenclature memory note.
- **S1L1 not hostable in the idealised model — characterised and
  documented** — `scripts/characterise_s1l1.py` shows no circular-coplanar
  topology reproduces the published 5.65 / 3.05 km/s V∞ anchors within
  ±0.3 km/s. The E-M-E-E topology was additionally tested
  (`scripts/characterise_s1l1_emee.py`) and does NOT close in
  circular-coplanar (V∞_E ≈ 25-39 km/s): the blocker is the MODEL (the
  ~154-d outbound leg is near-hyperbolic, exactly like the Aldrin case),
  NOT the topology. Real-ephemeris closure additionally needs multi-rev
  support in the maintenance engine. S1L1 gates remain xfail. This is a
  model limitation, not a solver bug.
- **Hollister & Menning Earth-Venus family individuated** — the single
  E-V placeholder was expanded into the 15-orbit
  `hollister-menning-1970-ev-orbit-01..15` family from the PRIMARY
  Hollister & Menning 1970 paper (now in `docs/refs/`), with V∞ from
  Table 3 (Vr × 29.785 EMOS) and a shared period 16 yr / k=10 (corrected
  from a wrong secondary "3.2 yr", which was the coplanar sub-orbit).
  New V0 data-integrity tests: `tests/data/test_hollister_family.py`,
  `tests/data/test_aldrin_establishment.py`.
- **Real-ephemeris cell optimiser IMPLEMENTED** — `optimise_cell_ephemeris`
  (was a `NotImplementedError` stub) now runs against DE440 with an
  asymmetric `tof_seed_days` option and an Aldrin parity test; `discover()`
  supports `optimiser="ephemeris"` (Part 2 wired).
- **Catalogue census frozen as a ratchet** — `tests/test_catalogue_rediscovery.py`
  (`EXPECTED_COVERAGE`) pins the 233-entry distribution; any catalogue
  change must update it in the same commit.

## In progress (with plan references)

- **SnLm multi-rev rediscovery** —
  `docs/superpowers/plans/2026-06-02-snlm-multirev-rediscovery.md`.
  Phase 1 Task 1 EXECUTED (S1L1 characterisation above). Phase 2 —
  plumbing multi-rev enumeration through `data/discover.py` and the
  sweep gate so the 202-entry `MULTI_ENCOUNTER_SEQUENCE` bucket becomes
  reachable by `discover` — IN PROGRESS.
- **Catalogue data-gap sourcing workflow** — a workflow is extracting
  Russell 2004 + Rogers 2012 tables and filling sourced catalogue gaps
  (goal: more usable golden orbits). Tracked in `data/MISSING_DATA.md`.

Two further plan files exist alongside the above and feed the same
effort: `2026-06-02-ml-multirev-lambert.md` and
`2026-06-02-multirev-3d-vem-ephemeris-roadmap.md`.

## Deferred (by explicit user decision — not being implemented now)

- **M8-Core (VEM 3-body)** — `docs/phases/m8-multibody-vem/plan.md`.
  Plan is written and revised (`period_basis` on `Cell`, same-body
  Tisserand bypass, N-agnostic loader, M8-Core / M8-UX split) but is
  DEFERRED. The Jones 2017 VEM member-list sourcing gap (Q&A item D)
  remains the data blocker for this phase.

## Blocked — needs HUMAN institutional / paywalled access (prioritised)

These are the highest-value human-actionable items: each unblocks
catalogue gaps that no amount of computation or open-web fetching can
fill. Listed in rough priority order.

1. **Russell 2004 dissertation tables (dominant gap source).** The PDF
   itself is open-access (`docs/refs/russell-2004-dissertation.pdf`;
   UT Austin `http://hdl.handle.net/2152/1253`) but binary-compressed —
   the in-progress extraction workflow is handling it. 772 of 790 gaps
   trace here. NOT human-blocked, but flagged as the single
   highest-leverage target; see `data/MISSING_DATA.md` §3.1.
2. **Friedlander/Niehoff/Byrnes/Longuski 1986, AIAA 86-2009-CP** —
   VISIT-1 / VISIT-2 V∞ at Earth and Mars. Paywalled
   (DOI `10.2514/6.1986-2009`). 4 "unknown" V∞ gaps; never digitised
   online (Q&A item C still partly open).
3. **McConaghy 2005 Purdue dissertation (AAI3166673) / AIAA 2002-4420
   full text** — U0L1 and Case 1 steady-state V∞ (and U0L1 return ToF).
   Both ResearchGate / ProQuest restricted. Note: the McConaghy 2006
   S1L1 orbital-element gap this would have closed is now resolved via
   Russell Table 4.9 (Q&A item B), so this is needed only for the
   remaining U0L1 / Case 1 / SnLm steady-state V∞ gaps. **NB S1L1,
   U0L1, and "Case 1" are THREE DISTINCT orbits** (Rogers 2012 Table 1:
   S1L1 a=1.30/e=0.257, Case 1 a=1.22/e=0.238, U0L1 a=2.05/e=0.563),
   not the same trajectory.
4. **Jones / Hernandez / Jesick 2017 (AAS 17-577)** — VEM triple-cycler
   member list (Q&A item D). Needed to lift M8 out of placeholder
   territory; lower priority while M8 is deferred.

## Catalogue census (frozen ratchet — 233 entries)

Source of truth: `tests/test_catalogue_rediscovery.py` `EXPECTED_COVERAGE`.

| Exclusion reason | Count |
|---|---|
| `multi_encounter_sequence` | 202 |
| `missing_leg_tofs` | 15 |
| `non_heliocentric` | 6 |
| `missing_vinf` | 5 |
| `constructible` | 2 |
| `not_two_body` | 2 |
| `missing_period` | 1 |
| **Total** | **233** |

The new `missing_leg_tofs` bucket (15) holds the individuated Hollister &
Menning E-V family (orbits 1-15): each has matched V∞ at Earth and Venus
but no per-leg ToFs encoded in `legs`. Expanding the old single E-V
placeholder moved one row out of `missing_vinf` (6 → 5) and one out of
`missing_period` (2 → 1), and added the 15-row bucket.

Through the idealised rediscovery gauntlet ~0/233 currently pass green:
the 2 `constructible` (Aldrin) entries are model-limited skips
(`EXPECTED_SKIPS`). Aldrin IS replicated on the real ephemeris (M6b)
and the canonical 2-synodic E–M anchor is replicated idealised. The
dominant blocker is the 202-entry (87%) `multi_encounter_sequence` /
SnLm multi-rev modelling bucket — exactly what the in-progress SnLm
multi-rev plan targets.

## Data-gap accounting (see `data/MISSING_DATA.md`)

788 data gaps across 207 entries:

- **403 "unknown"** — need sourcing (the bulk from Russell 2004, the
  rest from Rogers 2012 / McConaghy 2005 / VISIT).
- **201 "derive"** — COMPUTABLE by our Lambert solver, no sourcing
  needed.
- **184 "uncertain"** — topology-provisional; resolve as a by-product
  of extracting Russell per-leg elements.

The bulk of the gaps trace to the Russell 2004 dissertation. Source PDFs live in
`docs/refs/`: `russell-2004-dissertation.pdf`,
`rogers-2012-vinf-leveraging-cyclers-AIAA-2012-4746.pdf`,
`genova-aldrin-2015-earth-moon-cycler-AAS-15.pdf`,
`AAS-22-015-pascarella-pony-express.pdf`.

## Key findings

- **S1L1 nomenclature.** In McConaghy / Longuski / Byrnes, `S` and `L`
  denote consecutive Earth-to-Earth RESONANT INTERVALS, not Earth↔Mars
  legs. S1L1 therefore has NO direct Mars→Earth return leg — each
  vehicle is one-way; the crewed return needs a mirrored conjugate
  L1S1 cycler. The catalogue's old direct M→E "return leg" was a
  mismodelling, corrected in commit `f4074d0`.
- **S1L1, U0L1, and "Case 1" are three distinct orbits** — per Rogers
  2012 Table 1: S1L1 a=1.30/e=0.257, Case 1 a=1.22/e=0.238, U0L1
  a=2.05/e=0.563. No doc or entry should treat them as the same
  trajectory.
- **Idealised-model limitations are real, not bugs.** Both the Aldrin
  cycler and S1L1 have published anchors that the circular-coplanar
  idealised model cannot host. These are documented xfails / skips, and
  the now-implemented real-ephemeris cell optimiser is the intended way
  to close them.
- **Residual true data gaps** (need human / restricted access):
  VISIT-1/2 V∞ (Friedlander/Niehoff 1986 AIAA 86-2009-CP); U0L1 &
  Case 1 V∞ (McConaghy 2005 dissertation / AIAA 2002-4420 full text,
  ResearchGate/ProQuest restricted). The Aldrin establishment Mars-V∞
  is N/A (Earth-side V∞ leveraging). The 201 catalogue `derive`-kind
  data_gaps are computable (no sourcing). New reference PDFs added to
  `docs/refs/` this session: russell-2004-dissertation, rogers-2012,
  hollister-menning-1970 primary scan, hollister-rall-1970 NASA-CR,
  landau-longuski-2006-pt1, vasile-2005, genova-2016-phobos-deimos.

## Infra note

Background subagents in this environment appear to lack the `Bash` tool
(under test). Implementation and verification (pytest, git) are run in
the foreground by the orchestrator, not by dispatched subagents.

---

## ✓ Resolved (2026-05-31) — A. Aldrin orbital-element discrepancy

**Resolution:** commit `32c5eab` (errata: reconcile spec §3/§9/§9.1/§16.4)
plus the M3 implementation in `ba12554` which reproduces the literature
ellipse to ±0.002. See `docs/errata-investigation.md` §1 for the full
analysis: the spec's `a=1.659, e=0.41` is a resonance-construction
choice that is internally inconsistent with the same spec's 146-d
Earth-Mars leg (those elements imply 138.9 d). The literature value
`a=1.60, e=0.393` is internally consistent at 146 d and is what the
M3 patched-conic constructor produces. spec.md §9 now carries the
literature values; §9.1 documents the reconciliation.

**Original question (preserved as audit trail):**

`spec.md` §9 anchored the M3 gate to:

> Aldrin cycler: a ≈ 1.659 AU, e ≈ 0.41, perihelion ≈ 0.98 AU,
> aphelion ≈ 2.34 AU, E→M leg ≈ 146 d.

But the literature consistently reports:

> a = 1.60 AU, e = 0.393, perihelion = 0.97 AU, aphelion = 2.23 AU,
> E→M = 146 d
> (Rogers et al. 2012 Table 1; Russell 2004 Table 3.4 via Aphelion
> Ratio 1.47; Wikipedia citing Byrnes/Longuski/Aldrin 1993)

The gap on *a* is 0.06 AU — six times the M3 gate's `TOL_A_AU = 0.01`.
Either the spec's value set is wrong, the literature set is wrong, or
both refer to different "Aldrin cyclers".

---

## ✓ Resolved (2026-06-01) — B. McConaghy 2006 orbital elements (medium priority)

**Fully resolved by Russell 2004 Chapter 4 tables ingest** (commit
pending, 2026-06-01). Russell 2004 dissertation Table 4.9 (page 127)
row 1 carries the orbital data for the McConaghy 2006 "Notable" (S1L1)
cycler under Russell's own nomenclature `4.991gG2`:

- aphelion = 1.64 AU (row 1 of Table 4.9)
- V_inf at Earth = 4.99 km/s (Russell), vs McConaghy 2006 abstract = 4.7 km/s
- V_inf at Mars = 5.10 km/s (Russell), vs McConaghy 2006 abstract = 5.0 km/s
- E-M ToF = 150 days (Russell), vs McConaghy 2006 abstract = 153 days

Russell explicitly cross-references the two: dissertation line 7416
states "cycler 4.991gG2(#83) ... Also known as the 'S1L1' cycler",
line 5476 says "notable 'S1L1' cycler ... discovered first by
McConaghy et al. in Ref. 15", and line 8008 lists "the S1L1 cycler
(4.991Gg2), 8.049gGf2, and the Aldrin cycler" as the three most
promising designs.

**Entry 2 (`mcconaghy-2006-em-k2`) updated:**

- `orbit_elements.aphelion_au` backfilled with 1.64 AU.
- `orbit_elements.note` updated to cite Russell Table 4.9 row 1 as
  the source.
- `orbit_elements.a_au`, `e`, `perihelion_au` retain `null` because
  the cycler is a piecewise sequence of two generic-return arcs
  (g(1.4612,526.02,Ll) + G(2.8096,651.46,U)), not a single Keplerian
  ellipse. Each leg has its own (a, e), so the whole-cycler (a, e)
  are not well-defined; only the maximum aphelion is.
- `source_quotes.orbit_elements.aphelion_au` added with Russell Table
  4.9 citation.
- `notes:` block updated with the V_inf discrepancy analysis
  (McConaghy 4.7/5.0 vs Russell 4.99/5.10) and points to the
  new sibling entry `russell-ch4-4.991gG2` carrying Russell's
  circular-coplanar reference values.

The new entry `russell-ch4-4.991gG2` catalogues the same cycler
under Russell's framing — both entries are preserved per the
intended M7 collapse-via-canonical-signature semantics.

**Discrepancy noted, NOT silently resolved:** the 0.29 km/s V_inf E
difference and 3-day ToF difference between McConaghy 2006 and
Russell 2004 are larger than rounding alone could explain. The most
plausible reading per Russell's own text (line 7418 "essentially
ballistic for all launch dates ... consistent with the findings in
Ref. 15") is that McConaghy reports ephemeris-optimised values for
a realistic launch while Russell reports circular-coplanar simple-
model reference values. Both characterise the same trajectory.
Captured verbatim in both entries' `notes:` blocks for audit.

**Decision: Purdue dissertation acquisition is no longer required.**
Russell Table 4.9 provides the McConaghy 2006 orbital data with
acceptable precision for M7 matching.

**Audit trail (preserved):**

- 2026-06-01 morning: McConaghy 2005 SnLm broad-class family ingest
  landed (entries 45-47 = `mcconaghy-2005-em-case1`,
  `mcconaghy-2005-em-u0l1`, `mcconaghy-2005-em-snlm-broadclass-family`).
  Partial closure: documented SnLm class membership but did not
  backfill entry 2's null orbital elements (Purdue dissertation
  full text not accessible).
- 2026-06-01 evening: Russell 2004 Chapter 4 tables ingest landed
  (this work). Closed the gap by transcribing aphelion = 1.64 AU
  from Russell Table 4.9 row 1 and adding the cross-reference
  `russell-ch4-4.991gG2` entry. The Purdue dissertation acquisition
  is no longer the gating work; if acquired later it would only
  refine the aphelion value to higher precision and potentially
  add a true (a, e) per-leg breakdown.

**Original question (preserved as audit trail):**

The McConaghy 2006 abstract gives V∞ at Earth (4.7 km/s), V∞ at Mars
(5.0 km/s), and Earth–Mars ToF (153 d), but no orbital elements (a, e,
peri, apo). The full paper is paywalled at AIAA. **Without access to the
paper, we cannot fully specify the canonical signature for M7
matching** — finders that hit this cycler will get `null` matches on
the leg_elements field.

**Recommendation (historical):** the McConaghy 2005 Purdue PhD dissertation
(e-Pubs AAI3166673) is the open-access alternative containing the
broader SnLm taxonomy — queued for future ingest (task #34). When
ingested, the McConaghy 2006 "Notable" cycler should be cross-derived
from its dissertation analog (SnLm sibling family). _Now superseded by
the Russell Table 4.9 ingest._

---

## ◐ Partly resolved (2026-05-31) — C. VISIT-1 / VISIT-2 parameter inconsistency

**Resolution:** commit `b388b8d` — both VISIT entries now carry
arithmetic verification in `period.note` showing that the Rogers 2012
elements (a, e) are internally consistent with the 7-synodic / 14.95 yr
repeat period. The "Wikipedia swap" appears to be a different VISIT
variant (multiple "VISIT-1"-named cyclers exist in the literature).
The Rogers 2012 numbers are taken as authoritative and the YAML now
documents this choice. **Still open:** the original Niehoff 1985 / 1986
sources have not been consulted; if a future ingest of those originals
yields different elements, the entries should be re-evaluated.

**Original question (preserved):**

Wikipedia (citing McConaghy/Longuski/Byrnes 2002 p. 6) and Rogers et al.
2012 Table 1 give contradictory aphelion radii for VISIT-1 and VISIT-2:

| Source | VISIT-1 aphelion | VISIT-2 aphelion |
|---|---|---|
| Wikipedia (citing McConaghy 2002) | 1.89 AU | 1.45 AU |
| Rogers 2012 Table 1 | 1.40 AU | 1.67 AU |

The values appear to be *swapped* — i.e. Wikipedia's "VISIT-1" is
Rogers's "VISIT-2" and vice versa, OR they refer to different
"VISIT"-named cyclers in different papers (Niehoff published several
slightly different variants over 1985-91). Without the original
Niehoff documents (none online), this cannot be resolved.

---

## D. Jones 2017 VEM triple cyclers — full member list (high priority for M8)

The Jones/Hernandez/Jesick 2017 paper reports "thousands" of VEM triple
cyclers but the abstract gives only a family-level summary (average
transit V∞ < 5 km/s). The NTRS record's "downloads" section returned
HTTP 404. The ResearchGate PDF returned HTTP 403. **The catalogue
currently has only a family-seed entry; M7 matching will tag any VEM
finder hit as "probable-novel, flag for human review against Jones
2017" until member-level data can be ingested.**

The 2026-05-31 lunar+Jovian expansion added the EMEEVE 3-synodic
archetype as a sibling entry capturing the 6.4-yr branch geometry; the
2-synodic family-seed entry remains as-is. Both seeds are
deliberate-placeholders.

**Recommendation:** Obtain the Jones 2017 paper via NASA STI, JPL Open
Repository (handle hdl:2014/46418), or AIAA institutional access. Once
the member list is available, add each member as its own YAML entry
with full (sequence, period, V∞ multiset, leg elements).

---

## ✓ Resolved (2026-05-31) — E. spec.md §16.4 attribution correction

**Resolution:** commit `32c5eab`. spec.md §16.4 now cites "Jones,
Hernandez, Jesick (AAS 17-577, 2017)". The Longuski mis-attribution is
documented in `docs/errata-investigation.md` §3.

**Original question (preserved):**

spec.md §16.4 attributed the 2017 triple cyclers paper to "Longuski et
al." Per the NTRS record and the paper's title page, the authors are
Drew R. Jones, Sonia Hernandez, and Mark Jesick (all JPL). Longuski is
not an author.

---

## ✓ Resolved (2026-05-31) — F. spec.md §3 VEM beat period vs. Jones 2017 findings

**Resolution:** commit `32c5eab`. spec.md §3 was updated to clarify that
6.4-yr is the *lowest* natural beat with longer commensurabilities
(12.8 yr, 32 yr) also supporting closure per Jones 2017. The catalogue
preserves both the 2-synodic Jones family-seed entry and the
3-synodic EMEEVE archetype as separate records, accommodating both
readings of the abstract. See `docs/errata-investigation.md` §4 and the
EMEEVE entry's `period.note` for the detailed reconciliation.

**Original question (preserved):**

spec.md §3 says "the natural beat is ≈ 6.4 yr (3 × E–M ≈ 4 × E–V)." But
Jones et al. 2017 found 2-synodic-E-M (4.27 yr) VEM triple cyclers,
which is NOT the 6.4-yr beat. The beat period is sufficient for closure
in the simplified circular-coplanar model with strict commensurability,
but real eccentricities/inclinations + the b-plane DOF open up shorter
periods. M8's enumerator should NOT hard-code the 6.4-yr beat as the
only feasible VEM period.

---

## ✓ Resolved (2026-05-31) — G. Long Mars→Earth return leg of the Aldrin cycler

**Resolution:** commit `b388b8d`. The Aldrin outbound entry's
`legs[1].tof_days` was corrected from 519 d to 634 d with a derivation
source quote: T_cycler (779.8 d) − tof_outbound (146 d) = 633.8 d ≈
634 d. The simplified circular-coplanar model treats the return as a
single aggregate leg without sub-segment breakdown; a future ingest
of the Byrnes/Longuski/Aldrin 1993 paper could split it further if
needed.

**Original question (preserved):**

The 146-day Earth→Mars leg of the Aldrin cycler is well-cited. The
complementary Mars→Earth return — qualitatively described as "16 months
beyond Mars" by Wikipedia — was not cleanly tabulated in any single
primary source accessed during initial compilation. The YAML originally
recorded the return as `tof_days: 519` with an explicit "UNVERIFIED"
note.

---

## H. Out-of-paradigm work flagged 2026-05-31 (NOT in the catalogue)

When the catalogue was extended on 2026-05-31 to carry non-heliocentric
(lunar + Jovian + family-seed Saturnian) cyclers, additional bodies of
work were identified as **adjacent but out of the current cyclerfinder
paradigm**. They are recorded here as awareness for the user / future
implementers, but they are **deliberately NOT added to
`catalogue.yaml`** because:

1. cyclerfinder v1 models cyclers as patched-conic gravity-assist
   sequences with a V∞ + bend-angle abstraction at each flyby.
2. The papers below use fundamentally different mathematical paradigms
   (CR3BP invariant manifolds; low-thrust / solar sail) for which the
   V∞ + bend-angle signature is undefined.
3. Including them as YAML entries would make M7 novelty matching
   meaningless against them: any heliocentric or planet-centric finder
   hit would either falsely match them (signature comparison undefined)
   or never match them (`null` signatures everywhere). Better to flag
   them here and re-evaluate when / if the project adopts those
   modelling paradigms (cf. spec §2 stretch goals).

### H.1 Fantino, Alessi, Peláez Álvarez 2019 — Saturnian CR3BP manifold connections

| Field | Value |
|---|---|
| Title | "Connecting low-energy orbits in the Saturn system" |
| Authors | Elena Fantino, Elisa Maria Alessi, Jesús Peláez Álvarez |
| Venue | 18th Australian International Aerospace Congress (ISSFD-AIAC18), Melbourne, Australia, 24-26 February 2019, paper AIAC18 |
| URL (open) | <https://issfd.org/ISSFD_2019/ISSFD_2019_AIAC18_Fantino-Elena.pdf> |
| Mirror | <https://oa.upm.es/56463/> (Universidad Politécnica de Madrid open repository) |
| Methodology | CR3BP planar Lyapunov orbits + hyperbolic invariant manifolds + low-thrust patches; demonstrates a Tethys→Dione connection of 50 d using 9 kg propellant at 25 mN continuous thrust |
| Why excluded | The patched-conic + V∞ abstraction does not apply to manifold-based low-energy trajectories; the conserved quantity is the Jacobi constant, not V∞ |
| Re-evaluate when | The project adopts CR3BP modelling (would be the natural entry point for CR3BP catalogue ingestion) |

### H.2 Vergaaij & Heiligers 2018 — TU Delft solar-sail Earth-Mars cycler

| Field | Value |
|---|---|
| Title | "Time-optimal solar sail heteroclinic-like connections for an Earth-Mars cycler" |
| Authors | Merel Vergaaij, Jeannette Heiligers |
| Venue | *Acta Astronautica*, 2018, DOI <https://doi.org/10.1016/j.actaastro.2018.06.011> (ScienceDirect S0094576518303734) |
| URL | <https://research.tudelft.nl/en/publications/time-optimal-solar-sail-heteroclinic-like-connections-for-an-eart/> |
| Methodology | Direct pseudospectral optimisation + dynamical-systems heteroclinic connections between Earth-Moon L2 and Sun-Mars L1 libration-point orbits; requires a solar sail to close the connection (no ballistic solution exists). Time-optimal cyclers span ~3 synodic Earth-Mars periods. |
| Why excluded | Low-thrust / solar-sail propulsion is a spec §2 stretch goal, out of v1 scope. The trajectories are not ballistic and have no patched-conic V∞ signature. |
| Re-evaluate when | The project adopts low-thrust modelling. NB the original task brief referred to an "Earth-asteroid" version of this paper; the closest TU Delft paper actually found is this Earth-Mars version. If a distinct Earth-asteroid TU Delft paper exists, it would have the same out-of-paradigm classification. |

These flags are deliberately separate from the within-paradigm questions
A–G because those are gap-filling (missing numerics, inconsistent
secondary sources, attribution corrections) while H.1 and H.2 are
paradigm mismatches.
