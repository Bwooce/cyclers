# Outstanding questions — cyclerfinder catalogue

Long-form log of research questions, source-access gaps, parameter
contradictions, and out-of-paradigm flags encountered while compiling
`seed_cyclers.yaml`. The YAML's per-entry `notes:` field carries
short-form caveats; this file carries the discussion threads that
don't fit there.

**Resolution policy:** when a question is resolved, prefix its heading
with `✓ Resolved (YYYY-MM-DD)` and add a one-line pointer to the
resolution (commit SHA, spec section, errata-investigation note, etc.).
Do not delete the original question text — the audit trail matters.

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

## ◐ Partly resolved (2026-06-01) — B. McConaghy 2006 orbital elements (medium priority)

**Partial resolution:** SnLm broad-class family ingest landed on
2026-06-01 (commit pending) as catalogue entries 45-47:

- entry 45 `mcconaghy-2005-em-case1` — S2L1 family member Case 1, with
  a=1.22, e=0.238, peri=0.93, apo=1.51 from Rogers 2012 Table 1.
- entry 46 `mcconaghy-2005-em-u0l1` — U0L1 high-energy zero-loop
  one-synodic cycler, with a=2.05, e=0.563, peri=0.90, apo=3.20.
- entry 47 `mcconaghy-2005-em-snlm-broadclass-family` — citation-only
  family-seed for the McConaghy 2005 Purdue PhD dissertation
  (AAI3166673) capturing the SnLm taxonomy as the comprehensive
  open-access source.

**Still open:** the McConaghy 2006 'Notable' cycler's specific
orbital elements (a, e, peri, apo) are STILL null on entry 2
`mcconaghy-2006-em-k2`. The McConaghy 2005 dissertation full text
was not directly accessible at the 2026-06-01 ingest (Purdue
docs.lib.purdue.edu blocked to WebFetch; ProQuest dissertation full
text behind subscription; AIAA arc.aiaa.org 403 across all relevant
DOIs). The dissertation is expected to contain the per-cycler
orbital elements but those values were not transcribable from any
accessible secondary source.

The SnLm family-seed entry 47 documents which McConaghy SnLm class
the 'Notable' cycler belongs to (it is one of the S2L1 family
alongside Cases 1, 2, 3 with V_inf 4.7/5.0 km/s as a distinguishing
fingerprint), so M7 matching is now improved: a finder hit matching
the entry 2 V_inf multiset can be cross-referenced against the four
catalogued S2L1 members (entries 2, 4, 5, 45) by orbital-element
proximity, even though entry 2's own (a, e, peri, apo) remain null.

**Recommendation (still open):** acquire the McConaghy 2005 dissertation
full text (Purdue library print copy; ProQuest institutional
subscription; or direct request to the Purdue Advanced Astrodynamics
Concepts lab where the dissertation originated) and backfill entry 2's
orbital elements. Same source would unlock per-member numerics for
the S3L1 ballistic and S1L2 cycler families (currently not catalogued
as per-member entries; family seed entry 47 carries the citation).

**Original question (preserved as audit trail):**

The McConaghy 2006 abstract gives V∞ at Earth (4.7 km/s), V∞ at Mars
(5.0 km/s), and Earth–Mars ToF (153 d), but no orbital elements (a, e,
peri, apo). The full paper is paywalled at AIAA. **Without access to the
paper, we cannot fully specify the canonical signature for M7
matching** — finders that hit this cycler will get `null` matches on
the leg_elements field.

**Recommendation:** the McConaghy 2005 Purdue PhD dissertation
(e-Pubs AAI3166673) is the open-access alternative containing the
broader SnLm taxonomy — queued for future ingest (task #34). When
ingested, the McConaghy 2006 "Notable" cycler should be cross-derived
from its dissertation analog (SnLm sibling family).

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
`seed_cyclers.yaml`** because:

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
