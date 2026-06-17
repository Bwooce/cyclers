# Russell 2012 — Survey of Spacecraft Trajectory Design in Strongly Perturbed Environments

**Recovery-digest verdict note, 2026-06-17 AET.** Supersedes the file-only
entry in `2026-06-17-mars-cycler-wave-digest.md`. **JOURNAL-SURVEY
SCOPE:** every-page read (38-page JGCD article, fully readable in tool-
budget); no chapter-summary caveat needed.

## Header

- **Title:** *Survey of Spacecraft Trajectory Design in Strongly
  Perturbed Environments*
- **Author:** Ryan P. Russell
- **Affiliation (at publication):** University of Texas at Austin,
  Department of Aerospace Engineering and Engineering Mechanics
- **Author bio:** Russell was previously a member of the JPL Guidance,
  Navigation, and Control Section; Assistant Professor at Georgia Tech
  2007-2011 and then UT Austin from 2011. Senior Member AIAA. The same
  Russell whose 2004 UT Austin dissertation (hdl:2152/1253) anchors a
  large fraction of our catalogue's Earth-Mars cycler rows.
- **Venue:** *Journal of Guidance, Control, and Dynamics*, Vol. 35,
  No. 3, 2012, pp. 705-720
- **DOI:** 10.2514/1.56813
- **Manuscript version note:** the published PDF in
  `cyclers_pdf/papers/` is the **approved manuscript with errata for
  Table 1 and Figure 1** corrected (footnotes on pp. 4 and 5 of the
  manuscript flag the JGCD-printed typos). The manuscript paginates
  1/38-38/38 vs. the journal pages 705-720.
- **References:** 109 numbered.

## What the paper actually is

A career-survey article — Russell consolidates the methodology and
mission-design context for trajectories in **strongly perturbed
dynamical environments**: lunar / planetary-satellite systems, comets,
asteroids, and grand-tour transfers. The organisation is by
**dominant perturbing force**:

- §2.2-2.4 third-body gravity (TB) including patched-conics, grand
  tours, Tisserand maps, V∞ leveraging.
- §2.5-2.7 non-spherical gravity (NSG): zonal harmonics for
  large bodies, polyhedral / mascon models for small irregular
  bodies, frozen-orbit theory.
- §2.8 the Restricted Three Body Problem (RTBP) / Hill model /
  Elliptical RTBP, with stable/unstable manifold methodology and the
  endgame / leveraging tour design technique.
- §2.9 Earth-Moon system as a special case (large mass ratio,
  invalid zero-radius patched conics, Belbruno WSB rescue of Hiten).
- §2.10 solar radiation pressure (SRP), including small-body and
  solar-sail extremes.

The article's principal **contribution** (and the part most often cited)
is **Table 1 with its companion Figure 1** (pp. 4-5 of the manuscript):
a one-glance comparison of physical and dynamic parameters for ~45
representative target bodies (planets, Galilean / Saturnian / Uranian
/ Neptunian / Plutonian moons, asteroids, comets, and 2 NEAs) with
**quantitative α_TB and α_SRP perturbation ratios** evaluated at the
body radius and the Hill radius. The table is a **mission-design
calibration reference** — for a given target body, one can read off the
relative importance of TB vs. NSG vs. SRP and choose the right
modelling fidelity (averaged, periodic-orbit, full-ephemeris) before
investing in expensive numerical machinery.

The remainder of the article is a methodology survey, anchored in 109
references, **almost all of which are existing JPL / Russell-school /
Strange-Campagnola-Russell / Scheeres-Hu / Lara works** — the article
does not introduce new techniques so much as it situates Russell's
career body of work in a unified perturbations-dominated framing.

## Catalogue / KNOWN_CORPUS relevance

### Cycler-specific content (extremely light)

The word "cycler" appears **once** in the article body (Figure 4
caption, manuscript p. 10): the figure presents four example
trajectories side-by-side, the second of which is labelled "Titan and
Enceladus cycler tour at Saturn". No body-text section is devoted to
cyclers. The associated discussion (§2.3, manuscript pp. 11-13) treats
this as a patched-conics grand-tour example, **not** a continuous-
trajectory cycler in the Aldrin / SnLm sense.

### Russell & Strange "Planetary Moon Cycler Trajectories" — reference [21]

The single load-bearing cycler citation in Russell 2012 is **[21]:
Russell, R. P., Strange N.J., "Planetary Moon Cycler Trajectories,"
*Journal of Guidance, Control, and Dynamics*, Vol. 32, No. 1, 2009,
pp. 143-157**. This is invoked in §2.3 (manuscript p. 12) only as a
reference for the **V∞ globe visualisation tool** — *"The v∞ globe is
useful for both the design and post-process analyses of interplanetary
and inter-satellite gravity-assisted tours [18, 19, 20, 21]"*. No new
analytical content from the cycler paper is presented in Russell 2012.

**Reference [21] is the JGCD version of the 2008 AAS Honolulu paper**
already tracked in our `data/MISSING_DATA.md` §3.7 (Russell & Strange
2009 Jovian/Saturnian moon cyclers, classified `non_heliocentric`).
The two `russell-strange-2009-*-multimoon-family` rows in
`catalogue.yaml` (lines 8313 and 8453) are seed-only family entries
with all numerical fields null pending the AIAA-paywalled per-member
extraction.

### Strange-Campagnola-Russell KNOWN_CORPUS anchor

The Strange-Campagnola-Russell anchor in
`src/cyclerfinder/search/literature_check.py` (lines 308-325) is
**already** carrying:

- Strange-Russell-Buffington "Mapping the V-infinity globe" AAS 07-277
- Campagnola-Russell "Endgame Problem Part 1 + Part 2" JGCD 2010
- DOI `10.2514/1.45645`

Russell 2012 **adds no new references** to that anchor that aren't
already present in either the anchor's citation list or the broader
catalogue. The only relevant new reference would be the **survey
article itself** as a "tour through Russell's framework" citation. The
KNOWN_CORPUS anchor is structural / topology-bearing — it asserts
*"this body set + topology + method is published"* — so a survey
article that cites these methods doesn't strengthen the anchor; it
just describes it. **No KNOWN_CORPUS edit warranted by Russell 2012.**

### Direct catalogue impact

Zero. The article presents no new sourced cycler families. The
*"Titan and Enceladus cycler tour at Saturn"* figure caption refers
to existing Russell & Strange 2009 work which is already
catalogue-anchored (seed rows). No catalogue VO certifiable from
Russell 2012 directly.

### Possible reference resource

- **Table 1 + Figure 1** (manuscript pp. 4-5) are an excellent
  one-stop mission-design calibration reference for any new
  body-system extension to the catalogue. If a future task (e.g. an
  asteroid-leveraging or NEO cycler-precursor scan) needs to compute
  TB and SRP perturbation magnitudes for a target body, Table 1
  provides 45 rows of pre-computed values.
- **§2.9 (manuscript p. 26-27)** on the Earth-Moon system is a
  textbook-grade short summary of why the Earth-Moon system is
  uniquely difficult: 1/81 mass ratio (largest of all planetary
  satellites except Charon), 18o-28o inclination cycling, eccentric
  lunar orbit. Useful as a citable summary for any cislunar
  documentation.
- **§2.4 (manuscript pp. 14-15)** is the V∞-leveraging / Tisserand-
  graph summary, with the endgame methodology paper citations
  ([27, 30, 35]). Already known to our corpus.

### Verdict

**No KNOWN_CORPUS edits.** **No catalogue row changes.** The article is
a methodology survey, not a corpus-extension source. It cites work
already anchored in `literature_check.py`.

It is, however, an extremely useful **single-source reference card**
for the Russell-school perturbations-dominated mission-design
framework. Worth filing as the canonical citation when explaining the
Russell-Strange-Campagnola moon-tour endgame methodology to readers
new to the field.

## Errata vs the pre-fire survey

The pre-fire Mars-cycler wave digest filed Russell 2012 as: *"JGCD
survey paper on trajectory design in strongly perturbed environments
(Russell's career-survey article). Could yield methodology anchors for
KNOWN_CORPUS Strange/Campagnola/Russell entry. Light-read tagged for
future KNOWN_CORPUS extension."*

The deep-read **partially refines** that verdict:

1. **Yes**, the article is a Russell career-survey article, as the
   pre-read predicted.
2. **No**, it does **not** yield new methodology anchors for the
   Strange/Campagnola/Russell KNOWN_CORPUS entry. The corpus anchor
   already carries the Strange-Russell-Buffington 2007 V∞-globe paper
   and the Campagnola-Russell 2010 endgame papers (lines 309-325 of
   `literature_check.py`). Russell 2012 references the same papers but
   doesn't add new topology / method ground.
3. The article **does** introduce one citation worth flagging for any
   future cislunar / Hiten lit-pass: **[105] Wilson, R.S.,
   "Trajectory Design in the Sun-Earth-Moon Four Body Problem," PhD
   Thesis, Purdue University, Dec. 1998** — a Sun-Earth-Moon
   four-body PhD thesis predating Parker 2007's BLT thesis. This is
   already known territory (Hiten lineage) but is a useful pre-Parker
   cislunar-BCT reference.

The "could yield methodology anchors" framing was overly optimistic —
the article cites the methodology, but the anchor citations are
already in our corpus.

## Action items for parent

1. **No KNOWN_CORPUS additions** from Russell 2012 itself. The
   methodology citations the article relies on are already anchored
   in `src/cyclerfinder/search/literature_check.py` lines 309-325.
2. **Light citation enhancement (optional):** if the
   Strange/Campagnola/Russell KNOWN_CORPUS anchor wanted a **canonical
   survey reference** to point readers at, Russell 2012 JGCD
   (10.2514/1.56813) is the obvious choice — but adding it to the
   anchor's `citation` field is purely cosmetic; the topology /
   body-set / authors-tuple is already correctly captured.
3. **No catalogue row changes.** The two seed rows
   `russell-strange-2009-jovian-multimoon-family` and
   `russell-strange-2009-saturnian-multimoon-family` remain
   null-numeric pending the AIAA-paywalled per-member extraction.
   Russell 2012 does **not** provide that extraction (it just cites
   the underlying paper).
4. **Acquisitions wishlist addition (low priority):** *Wilson, R.S.,
   "Trajectory Design in the Sun-Earth-Moon Four Body Problem," PhD
   Thesis, Purdue University, Dec. 1998* — Russell 2012 ref [105]. A
   pre-Parker cislunar four-body methodology thesis. Useful if #316
   cislunar BCT lit-pass needs a pre-Parker-2007 methodology
   reference. Probably already known to the Parker-Anderson lineage.
5. **Reference resource only.** Table 1 and Figure 1 (manuscript
   pp. 4-5) are valuable as a body-system perturbation-strength
   reference for any future catalogue extension to new primaries.
6. **No tracked task created.** Pure reference resource.

## Source-page references

- Title page, DOI, errata note: manuscript p. 1.
- Table 1 (45-body perturbation reference): manuscript p. 4.
- Figure 1 (SRP vs. TB plot): manuscript p. 5.
- §2.3 patched conics / grand tours / Figure 4 cycler caption:
  manuscript pp. 10-13.
- §2.4 Tisserand maps + V∞ leveraging: manuscript pp. 14-15.
- §2.8 Restricted Three-Body Problem + manifolds: manuscript
  pp. 20-23.
- §2.8.2 Hill model + Hill sphere + planetary-satellite orbiters:
  manuscript pp. 23-25.
- §2.9 Earth-Moon system: manuscript pp. 26-27.
- Reference [21] Russell-Strange "Planetary Moon Cycler Trajectories"
  JGCD 2009: manuscript p. 31.
- Reference [27] Strange-Campagnola-Russell "Leveraging Flybys of Low
  Mass Moons to Enable an Enceladus Orbiter" AAS 09-435 2009:
  manuscript p. 32.
- Reference [105] Wilson 1998 Sun-Earth-Moon four-body PhD thesis
  (Purdue): manuscript p. 38.
- Reference [107] Belbruno 2004 *Capture Dynamics* textbook:
  manuscript p. 38.
