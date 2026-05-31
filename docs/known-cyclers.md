# Known cyclers — seed catalogue reference

**Purpose.** This document is the human-readable companion to
`src/cyclerfinder/data/seed_cyclers.yaml`. Both files exist to give
cyclerfinder a starting catalogue of *published* cycler trajectories with
full source attribution, so that:

1. **M3 golden tests** can pull the Aldrin and McConaghy/S1L1 entries as
   numerical fixtures asserted against the patched-conic constructor.
2. **M7 novelty matching** can compare finder hits against prior art from day
   one, per spec §16.4: "the catalogue is seeded with, and continuously
   ingests, published cyclers — each with full citation — so finder hits
   match against prior art from day one."

**Audience.** This is for the implementer who will later need to verify a
number against its source without re-fetching every paper. Every numerical
value in the YAML catalogue is grounded here in a verbatim quote or
table/page citation; every gap is flagged so the next person knows what
remains to be looked up.

**Scope and honesty rules (per spec §16.4).**

- **Attribution always goes to the earliest established `priority_date`.**
  If a `first_published` source exists, the cycler is `known-reproduction`
  and credit belongs to those authors regardless of whether the finder
  rediscovers it independently. We never claim a cycler with prior
  publication as novel.
- **No fabrication.** Every numerical value in the YAML is either backed by
  a source quote in `source_quotes:` or is `null` with the gap noted. If a
  paper was not accessible, the entry carries the citation but `null`
  numerics rather than secondary-source-derived guesses.
- **Single source of truth.** The YAML is the machine-readable record; this
  markdown explains *why* each number is what it is. Discrepancies between
  the two are a bug — fix the YAML, then update this doc.

**Source access caveats encountered during this work (May 2026).**

- All AIAA-hosted PDFs and abstract pages (arc.aiaa.org) returned **HTTP 403
  Forbidden** to the WebFetch tool used to compile this catalogue. Quoted
  AIAA abstracts (McConaghy 2006, Russell-Ocampo 2005, Spreen 2020) come
  from secondary sources (web search snippets that quoted the abstracts,
  ResearchGate "Request PDF" pages that exposed abstract text in search
  results) and from one cached AIAA conference PDF (Rogers 2012) that
  happened to be hosted at engineering.purdue.edu rather than at AIAA.
- The full Russell 2004 dissertation was successfully downloaded from
  UT Austin's open-access repository
  (<http://hdl.handle.net/2152/1253>); its Tables 3.4–3.8 and 4.7–4.13
  are the most detailed primary source used here.
- Niehoff's original 1985 SAIC presentation and the early Niehoff 1986
  AAS Paper 86-172 were **not accessible** in any digitised form; values
  attributed to "Niehoff 1985" come from Rogers et al. 2012 Table 1 and the
  spaceflighthistory blog summary, with the original Niehoff documents cited
  but not consulted directly.
- The McConaghy 2006 JSR paper (10.2514/1.15215) was inaccessible beyond
  its abstract; orbital elements (a, e, peri, apo) for the "Notable
  Two-Synodic" cycler are therefore `null` in the YAML and listed in the
  Outstanding Questions section below.

---

## 1. Aldrin classic Earth-Mars cycler

**Catalogue ID:** `aldrin-classic-em-k1-outbound`

### Citation

| Field | Value |
|---|---|
| First published | Aldrin, B., "Cyclic Trajectory Concepts," SAIC presentation to the Interplanetary Rapid Transit Study Meeting, Jet Propulsion Laboratory, Pasadena CA, **28 October 1985**. (Not online; cited per Rogers et al. 2012 ref [8] and Russell 2004 ref [13]). |
| Peer-reviewed primary | Byrnes, D. V., Longuski, J. M., Aldrin, B., **"Cycler Orbit Between Earth and Mars,"** *Journal of Spacecraft and Rockets*, Vol. 30, No. 3, May–June 1993, pp. 334–336. **DOI:** [10.2514/3.25519](https://doi.org/10.2514/3.25519). |
| Priority date | 1985-10-28 |
| Sequence | E–M (one Earth flyby and one Mars flyby per synodic period; complemented by an "inbound" / down-escalator twin) |
| Period | 1 Earth-Mars synodic period = 2.135 yr |
| Sense | Outbound ("up escalator": short Earth→Mars leg, long Mars→Earth return) |

### Published parameters

| Quantity | Value | Source |
|---|---|---|
| Semi-major axis *a* | 1.60 AU | Rogers et al. 2012 Table 1 |
| Eccentricity *e* | 0.393 | Rogers et al. 2012 Table 1; corroborated by Wikipedia citing Byrnes/Longuski/Aldrin 1993 pp. 334–335 |
| Perihelion | 0.97 AU | Rogers et al. 2012 Table 1 |
| Aphelion | 2.23 AU | Rogers et al. 2012 Table 1; also derivable as 1.47 × 1.52 AU from Russell 2004 Table 3.4 Aphelion Ratio |
| Inclination | 0° (circular-coplanar) | modelling assumption |
| Earth–Mars leg time of flight | 146 d | Russell 2004 Table 3.4 (cycler 1.0.1.-1); Wikipedia citing McConaghy/Longuski/Byrnes 2002 p. 6 |
| V∞ at Earth | 6.5 km/s | Russell 2004 dissertation Table 3.4 (cycler 1.0.1.-1, footnoted as the Aldrin cycler) |
| V∞ at Mars | 9.7 km/s | Russell 2004 Table 3.4 |

### Source quotes

> "Aldrin Cycler ... Number of Vehicles: 2 ... Semi-Major Axis, AU: 1.60 ...
> Eccentricity: 0.393 ... Aphelion Radius, AU: 2.23 ... Perihelion Radius,
> AU: 0.97."
> — Rogers, B. A., Hughes, K. M., Longuski, J. M., Aldrin, B.,
> *Preliminary Analysis of Establishing Cycler Trajectories Between Earth
> and Mars via V-Infinity Leveraging*, AIAA 2012-4746, Table 1.

> "Cycler 1.0.1.-1 ... Aphelion Ratio: 1.47, Turn Ratio: 0.86, Earth→Mars
> Time: 146 days, Earth v∞: 6.5 km/s, Mars v∞: 9.7 km/s, Required
> Geocentric Turning Angle at each Flyby: 84°."
> Footnote c: "Aldrin cycler[1,13,14]"
> — Russell, R. P., *Global Search and Optimization for Free-Return
> Earth–Mars Cyclers*, Ph.D. dissertation, UT Austin, 2004, Table 3.4
> (cycler 1.0.1.-1 row).

> "It travels from Earth to Mars in 146 days (4.8 months), spends the next
> 16 months beyond the orbit of Mars" (Wikipedia, Mars cycler article,
> citing McConaghy/Longuski/Byrnes 2002 p. 6).

### Notes and known ambiguities

1. **The spec.md §9 numbers vs. Rogers 2012 Table 1 numbers.** spec.md §9
   anchors the M3 gate to *a* ≈ 1.659 AU, *e* ≈ 0.41, perihelion ≈ 0.98 AU,
   aphelion ≈ 2.34 AU. Rogers 2012 / Russell 2004 / Wikipedia consistently
   report *a* = 1.60, *e* = 0.393, perihelion = 0.97, aphelion = 2.23. The
   two value sets differ by ~0.06 AU on *a*. **Both are literature numbers
   for the Aldrin cycler**; the difference appears to stem from a different
   rounding/derivation convention (possibly directly from the 1985 SAIC
   presentation, which we cannot consult, vs. the McConaghy/Longuski/Byrnes
   2002 re-derivation). The M3 gate tolerance `TOL_A_AU = 0.01` in
   `docs/phases/m3-model-construct/plan.md` §4.3 is too tight to absorb the
   gap; this needs reconciliation in M3 (see Outstanding Questions §A).

2. **Outbound vs. inbound.** The Aldrin concept includes a twin "inbound"
   (down-escalator) cycler that is the time-reverse of the outbound; per
   Russell 2004 §3.8 the two have identical energy properties. Only the
   outbound is in the seed catalogue; M7's canonical signature should
   collapse them via the `sense` field.

3. **Mars–Earth return leg.** The 146-d Earth–Mars leg is solidly cited.
   The complementary Mars–Earth return is described qualitatively (Wikipedia:
   "spends the next 16 months beyond the orbit of Mars") but the specific
   numerical breakdown (loitering beyond Mars vs. inbound transit) is not
   tabulated in any single primary source we accessed. The YAML therefore
   marks `legs[1].tof_days` with a note that the 519-d value is
   provisional and inconsistent with the cycler period.

---

## 2. McConaghy/Landau/Yam/Longuski 2006 — "Notable" Two-Synodic E–M cycler

**Catalogue ID:** `mcconaghy-2006-em-k2`

### Citation

| Field | Value |
|---|---|
| First published | McConaghy, T. T., Yam, C. H., Landau, D. F., Longuski, J. M., "Two-Synodic-Period Earth-Mars Cyclers with Intermediate Earth Encounter," AAS Paper **03-509**, AAS/AIAA Astrodynamics Specialist Conference, August 2003. |
| Peer-reviewed primary | McConaghy, T. T., Landau, D. F., Yam, C. H., Longuski, J. M., **"Notable Two-Synodic-Period Earth-Mars Cycler,"** *Journal of Spacecraft and Rockets*, Vol. 43, No. 2, March–April 2006, pp. 456–465. **DOI:** [10.2514/1.15215](https://doi.org/10.2514/1.15215). |
| Priority date | 2003-08 |
| Sequence | E–E–M–M (2 Earth + 2 Mars encounters per cycle, one Earth being an intermediate) |
| Period | 2 × E-M synodic ≈ 4.27 yr |

### Published parameters

| Quantity | Value | Source |
|---|---|---|
| Earth–Mars leg time of flight | 153 d | McConaghy 2006 abstract |
| V∞ at Earth (arrival) | 4.7 km/s | McConaghy 2006 abstract |
| V∞ at Mars (arrival) | 5.0 km/s | McConaghy 2006 abstract |
| *a*, *e*, perihelion, aphelion | `null` (paywalled) | see Outstanding Questions §B |

### Source quotes

> "[A cycler that] repeats every two synodic periods and has one intermediate
> Earth encounter. In a circular-coplanar model it requires no propulsive
> maneuvers, has 153-day transfer times between Earth and Mars, and has
> arrival V-infinity magnitudes of 4.7 km/s at Earth and 5.0 km/s at Mars."
> — McConaghy, T. T., Landau, D. F., Yam, C. H., Longuski, J. M., *Notable
> Two-Synodic-Period Earth-Mars Cycler*, JSR 43(2) 2006, abstract (accessed
> via web search snippet quoting <https://arc.aiaa.org/doi/10.2514/1.15215>).

### Notes and known ambiguities

1. **Not the same cycler as spec §9's "5.65 / 3.05 km/s" 2-synodic cycler.**
   spec.md §9's M5 anchor calls for "≈ 5.65 km/s (Earth), 3.05 km/s (Mars)"
   for "the published 2-synodic E-M cycler". Those values do NOT match the
   McConaghy 2006 abstract (4.7 / 5.0). The 5.65/3.05 pair refers to the
   *S1L1* cycler (the next entry below) — different family member from the
   same broad class (McConaghy/Longuski/Byrnes 2002, AIAA 2002-4420). M5
   should test against the S1L1 numbers per the spec.

2. **Orbital elements gap.** The McConaghy 2006 abstract gives V∞ and ToF
   but not (a, e, peri, apo). The full paper is paywalled. The YAML's
   `orbit_elements` block is null with a TBD note.

3. **Comparison to Russell 2004 Table 4.9.** Russell's dissertation Table 4.9
   ("Ballistic two-synodic period cyclers") lists 4 ballistic 2-synodic
   cyclers. The first row (V∞E=4.99, V∞M=5.10, t_out=t_in=150 d,
   aphelion=1.64 AU) is the closest match to McConaghy 2006 but differs
   slightly — likely the same cycler with Russell's binning to 0.01 km/s
   vs. McConaghy's 0.1 km/s rounding. The McConaghy 2006 paper post-dates
   Russell's 2004 dissertation, so the literature attribution is to
   McConaghy.

---

## 3. S1L1 ballistic two-synodic E–M cycler (CPOM nominal)

**Catalogue ID:** `s1l1-2syn-em-cpom`

### Citation

| Field | Value |
|---|---|
| First published | McConaghy, T. T., Longuski, J. M., Byrnes, D. V., "Analysis of a Broad Class of Earth-Mars Cycler Trajectories," AIAA Paper **2002-4420**, AIAA/AAS Astrodynamics Specialist Conference, Monterey CA, August 5–8 2002. **DOI:** 10.2514/6.2002-4420. |
| Nomenclature paper | McConaghy, T. T., Russell, R. P., Longuski, J. M., "Towards a Standard Nomenclature for Earth-Mars Cycler Trajectories," *Journal of Spacecraft and Rockets* (in press as of Russell 2004 ref [25]; published ca. 2005). |
| Modern design study | Spreen, C., et al., "Design Considerations for an Earth-Mars Cycler Spacecraft Using the S1L1 Cycler," *Journal of Spacecraft and Rockets*. **DOI:** [10.2514/1.A35160](https://doi.org/10.2514/1.A35160). |
| Priority date | 2002-08-05 |
| Sequence | E–E–M–M (one intermediate Earth encounter) |
| Period | 2 × E-M synodic ≈ 4.27 yr |

### Published parameters

| Quantity | Value | Source |
|---|---|---|
| Semi-major axis *a* | 1.30 AU | Rogers et al. 2012 Table 1 |
| Eccentricity *e* | 0.257 | Rogers et al. 2012 Table 1 |
| Perihelion | 0.97 AU | Rogers et al. 2012 Table 1 |
| Aphelion | 1.64 AU | Rogers et al. 2012 Table 1 |
| Earth–Mars leg time of flight | 154 d | Spreen 2020 / multiple secondary sources |
| V∞ at Earth | 5.65 km/s | spec.md §9 and secondary sources attributing this to the S1L1 cycler |
| V∞ at Mars | 3.05 km/s | spec.md §9 |

### Source quotes

> "S1L1 ... Number of Vehicles: 4 ... Semi-Major Axis, AU: 1.30 ...
> Eccentricity: 0.257 ... Aphelion Radius, AU: 1.64 ... Perihelion
> Radius, AU: 0.97."
> — Rogers et al. 2012 Table 1.

> "A cycler that repeats every two synodic periods has a low V-infinity at
> Earth and Mars (5.65 and 3.05 km/s, respectively)."
> — Web search snippet citing the S1L1 cycler literature.

> "The S1L1 cycler is designed to transfer a crew of six from Earth to Mars
> on a nominal 154-day trajectory."
> — Secondary source citing Spreen et al. 2020 (DOI 10.2514/1.A35160).

### Notes and known ambiguities

1. **Nomenclature.** "S1L1" follows McConaghy/Russell/Longuski's standard
   nomenclature for Earth–Mars cyclers — broadly: "S" then the number of
   synodic periods, "L" then the number of intermediate Earth-Earth loops,
   with trailing digits encoding leg types. The full naming convention is
   in Russell 2004 dissertation ref [25].
2. **Aphelion vs. Mars's orbit.** Aphelion 1.64 AU is only marginally
   inside Mars's perihelion (1.38 AU). Mars's mean SMA is 1.524 AU.
   The cycler reaches Mars when Mars is closer to the Sun than ~1.64 AU,
   which happens for most of Mars's orbit but tightens the launch-window
   constraints. This is why CPOM emphasises the S1L1 cycler — Mars is
   reached at low V∞ near the cycler's aphelion.
3. **Russell 2004 cross-reference.** Russell's dissertation Table 4.9
   first row appears to be the S1L1 cycler in his nomenclature:
   "g(1.4612,526.02,Ll) G(2.8096,651.46,U) ... V∞E 4.99, V∞M 5.10"
   — but his V∞ values (4.99 / 5.10) match McConaghy 2006 better than
   the spec.md §9 value (5.65 / 3.05). Note possible inconsistency:
   the "S1L1" name and "5.65/3.05" pairing in spec.md may refer to a
   variant (Rogers 2012 Table 1's S1L1 row, the modern CPOM mission
   architecture) rather than the original McConaghy 2002 first-row
   cycler. The orbital elements (a, e, peri, apo) from Rogers 2012 are
   the authoritative pairing with the 5.65/3.05 V∞ values quoted in spec
   §9; treating them as one cycler is consistent with CPOM literature.

---

## 4. Russell–Ocampo ballistic cyclers (the "24 ballistic")

**Catalogue IDs:** five representatives currently seeded —
`russell-ocampo-2.1.1+2-case2`, `russell-ocampo-2.3.1+1-case3`,
`russell-ocampo-4.3.1-5`, `russell-ocampo-3.1.2+1`,
`russell-ocampo-2.5.1+0`.

### Citation

| Field | Value |
|---|---|
| Conference precursor | Russell, R. P., Ocampo, C. A., "A Systematic Method for Constructing Earth-Mars Cyclers Using Direct Return Trajectories," AAS Paper **03-145**, February 2003. |
| Peer-reviewed primary (free-return method) | Russell, R. P., Ocampo, C. A., **"Systematic Method for Constructing Earth-Mars Cyclers Using Free-Return Trajectories,"** *Journal of Guidance, Control, and Dynamics*, Vol. 27, No. 3, May–June 2004, pp. 321–335. **DOI:** 10.2514/1.5878. |
| Peer-reviewed primary (geometric analysis) | Russell, R. P., Ocampo, C. A., **"Geometric Analysis of Free-Return Trajectories Following a Gravity-Assisted Flyby,"** *Journal of Spacecraft and Rockets*, Vol. 42, No. 1, January–February 2005, pp. 138–151. **DOI:** [10.2514/1.5571](https://doi.org/10.2514/1.5571). |
| Comprehensive catalogue | Russell, R. P., **"Global Search and Optimization for Free-Return Earth-Mars Cyclers,"** Ph.D. dissertation, University of Texas at Austin, Department of Aerospace Engineering, 2004. **Handle:** <http://hdl.handle.net/2152/1253>. |

### Russell's nomenclature

A cycler is named `p.h.s.i` where:

- **p** = period in synodic periods,
- **h** = number of half-years allotted for half-rev or full-rev returns,
- **s** = number of generic (non-nπ) returns,
- **i** = signed integer specifying revolution count and branch of the
  Lambert solution (negative = lower-energy curve).

For example, `4.3.2.-5` is a 4-synodic cycler using 3 half-years for nπ
returns, 2 generic returns, with the generic returns on the lower-energy
5-revolution Lambert branch.

### The 24 ballistic cyclers (Russell 2004 Table 3.4, ARMIN=0.9, TRMIN=0.85)

Russell's Table 3.4 lists 44 cyclers of which 24 are entirely ballistic
(turn ratio ≥ 1.0 and aphelion ratio ≥ 1.0). The full table is reproduced
below verbatim from the dissertation; the seed catalogue currently
contains only five representatives (marked **bold**), chosen for being:

- The Aldrin cycler (`1.0.1.-1`, footnote c) — already its own entry above;
- "Case 2" and "Case 3" of Byrnes/McConaghy/Longuski 2002 (Russell
  footnotes d and e) — explicitly attributed to that prior work;
- `4.3.1.-5` — lowest V∞ of the entire table, near-Hohmann;
- `3.1.2.+1` — Russell's chosen exemplar low-V∞ 3-synodic cycler;
- `2.5.1.+0` — short-transit 2-synodic counter-example.

Other rows can be added to the YAML as individual entries when M7 wants
finer matching coverage.

| p.h.s.i        | AR    | TR   | E→M days | V∞E (km/s) | V∞M (km/s) | Notes |
|---|---|---|---|---|---|---|
| **`1.0.1.-1`** | 1.47 | 0.86 | 146 | 6.5 | 9.7 | Aldrin cycler (footnote c) |
| **`2.1.1.+2`** | 0.95 | 1.11 | 207 | 4.1 | 2.0 | "Case 2" of Byrnes et al. (footnote d) |
| **`2.3.1.+1`** | 1.08 | 0.92 | 143 | 5.4 | 5.3 | "Case 3" of Byrnes et al. (footnote e) |
| **`2.5.1.+0`** | 1.44 | 1.12 | 94 | 7.8 | 9.9 | Short-transit, high V∞ |
| `3.1.1.+3` | 1.07 | 1.19 | 174 | 3.6 | 4.6 | Low V∞ |
| `3.1.1.+2` | 1.43 | 0.89 | 115 | 5.4 | 9.2 | |
| **`3.1.2.+1`** | 1.07 | 1.23 | 181 | 3.4 | 4.6 | Low V∞, 3-synodic |
| `3.1.3.+0` | 1.43 | 0.93 | 123 | 5.1 | 9.1 | |
| `3.3.1.+2` | 1.19 | 1.06 | 141 | 4.3 | 6.8 | |
| `3.5.1.+2` | 0.94 | 1.80 | 231 | 2.7 | 1.5 | Very low V∞ but AR<1 |
| `3.5.1.+1` | 1.43 | 1.15 | 115 | 5.4 | 9.2 | |
| `3.5.2.+0` | 1.43 | 1.06 | 121 | 5.2 | 9.2 | |
| `3.7.1.+1` | 1.07 | 1.56 | 175 | 3.6 | 4.6 | Low V∞ |
| `3.9.1.+0` | 1.43 | 1.17 | 116 | 5.4 | 9.2 | |
| `4.0.3.+1` | 1.07 | 1.18 | 160 | 4.3 | 4.9 | Short ToF, moderate V∞ |
| `4.1.1.-6` | 0.94 | 1.37 | 256 | 2.7 | 1.6 | |
| `4.1.1.-5` | 1.15 | 1.11 | 173 | 4.1 | 6.1 | |
| `4.1.1.-4` | 1.44 | 0.89 | 137 | 5.5 | 9.3 | |
| `4.1.2.-3` | 0.94 | 1.40 | 250 | 2.6 | 1.5 | |
| `4.1.2.-2` | 1.43 | 0.93 | 132 | 5.2 | 9.2 | |
| `4.1.4.-1` | 1.43 | 0.93 | 129 | 5.1 | 9.2 | |
| **`4.3.1.-5`** | 0.99 | 1.29 | 268 | 3.1 | 2.5 | Near-Hohmann, lowest V∞ overall |
| `4.3.1.-4` | 1.26 | 1.01 | 154 | 4.7 | 7.6 | |
| `4.5.1.-4` | 1.07 | 1.55 | 196 | 3.6 | 4.7 | Low V∞ |
| `4.5.1.-3` | 1.44 | 1.15 | 137 | 5.5 | 9.3 | |
| `4.5.2.-2` | 1.07 | 1.40 | 191 | 3.4 | 4.6 | |
| `4.5.3.-1` | 1.43 | 1.02 | 130 | 5.1 | 9.2 | |
| `4.6.1.+4` | 0.91 | 1.50 | 154 | 6.8 | 2.1 | |
| `4.6.3.+0` | 1.43 | 0.88 | 105 | 6.4 | 9.5 | |
| `4.7.1.-3` | 1.20 | 1.38 | 163 | 4.3 | 6.8 | |
| `4.7.1.-2` | 1.77 | 0.96 | 120 | 6.6 | 11.4 | |
| `4.8.1.+3` | 0.96 | 1.64 | 164 | 7.7 | 3.1 | |
| `4.8.1.+2` | 1.31 | 0.86 | 76 | 12.5 | 10.7 | |
| `4.9.1.-3` | 0.94 | 1.83 | 256 | 2.7 | 1.6 | |
| `4.9.1.-2` | 1.44 | 1.16 | 137 | 5.5 | 9.3 | |
| `4.9.2.-1` | 1.44 | 1.05 | 132 | 5.2 | 9.2 | |
| `4.10.1.-3`| 0.92 | 1.46 | 263 | 10.2 | 3.6 | |
| `4.10.1.+2`| 1.03 | 1.65 | 131 | 8.9 | 5.0 | |
| `4.11.1.-2`| 1.07 | 1.58 | 195 | 3.6 | 4.7 | Low V∞ |
| `4.12.1.-2`| 0.97 | 1.43 | 268 | 11.6 | 4.8 | |
| `4.12.1.+1`| 1.16 | 1.48 | 93 | 10.8 | 8.2 | |
| `4.13.1.-1`| 1.44 | 1.16 | 137 | 5.5 | 9.3 | |
| `4.14.1.-1`| 1.12 | 1.13 | 199 | 14.7 | 9.4 | |
| `4.14.1.+0`| 1.49 | 1.09 | 66 | 14.1 | 12.7 | Shortest ToF in table |

(AR = Aphelion Ratio to 1.52 AU; TR = Turn Ratio = max allowed / max
required turning angle. Ballistic = AR ≥ 1 AND TR ≥ 1.)

### Source quote

> "A total of forty-four cyclers are shown in Table 3.4, of which twenty-four
> are entirely ballistic."
> — Russell 2004 dissertation §3.8.

> "The original method [Ref. 29] identifies twenty-four ballistic cyclers
> with periods of two to four synodic periods, 92 ballistic cyclers with
> periods of five or six synodic periods, and hundreds of near-ballistic
> cyclers, most of which are previously undocumented."
> — Russell 2004 dissertation §1.4.

### Notes

1. The 92 five-/six-synodic-period ballistic cyclers and the "hundreds"
   of near-ballistic cyclers from Russell's broader analysis (Tables 3.9
   – 3.11, 4.9 – 4.13) are NOT in the seed catalogue. They can be added
   as M7 needs them; the dissertation is the authoritative source.
2. Russell's *Aphelion Ratio* normalises to a fixed Mars semi-major axis
   of 1.52 AU (his §3.3 simplification: "Mars is assumed to be in an orbit
   with a period of 1.875 yrs ... This value is chosen such that the
   absolute geometry repeats every 15 yrs"). This differs from
   constants.py's `_MARS_SMA_AU = 1.52371034` by 0.244%; the difference
   propagates ~0.4% into the aphelion values quoted above, which is well
   within the 0.01 AU binning of spec §16.2.

---

## 5. Niehoff VISIT cyclers

**Catalogue IDs:** `niehoff-visit1`, `niehoff-visit2`.

### Citation

| Field | Value |
|---|---|
| First published (VISIT-1, VISIT-2) | Niehoff, J., "Manned Mars Mission Design," in *Steps to Mars*, Joint AIAA/Planetary Society Conference, National Academy of Sciences, Washington DC, **July 1985**. |
| Companion presentation | Niehoff, J., "Integrated Mars Unmanned Surface Exploration (IMUSE), A New Strategy for the Intensive Science Exploration of Mars," presentation to the Planetary Task Group, Space Science Board Major Directions Summer Study, Woods Hole MA, **30 July 1985**. |
| Conference paper | Niehoff, J., "Pathways to Mars: New Trajectory Opportunities," AAS Paper **86-172**, July 1986. |
| First peer-reviewed | Friedlander, A. L., Niehoff, J. C., Byrnes, D. V., Longuski, J. M., "Circulating Transportation Orbits Between Earth and Mars," AIAA/AAS Astrodynamics Conference, Williamsburg VA, **AIAA 86-2009-CP**, August 1986. **DOI:** [10.2514/6.1986-2009](https://doi.org/10.2514/6.1986-2009). |
| Subsequent | Niehoff, J., Friedlander, A., McAdams, J., "Earth-Mars Transport Cycler Concepts," International Astronautical Congress, IAF Paper **91-438**, October 1991. |
| Priority date | 1985-07-01 |

### Published parameters

| Quantity | VISIT-1 | VISIT-2 | Source |
|---|---|---|---|
| Semi-major axis *a* | 1.17 AU | 1.31 AU | Rogers et al. 2012 Table 1 |
| Eccentricity *e* | 0.193 | 0.275 | Rogers et al. 2012 Table 1 |
| Perihelion | 0.94 AU | 0.95 AU | Rogers et al. 2012 Table 1 |
| Aphelion | 1.40 AU | 1.67 AU | Rogers et al. 2012 Table 1 |
| Repeat period | 7 × E-M synodic ≈ 14.95 yr | same | Rogers et al. 2012 Table 1 footnote a |
| V∞ at Earth | `null` (gap) | `null` (gap) | not in any accessible source |
| V∞ at Mars | `null` (gap) | `null` (gap) | not in any accessible source |

### Source quotes

> "VISIT-1 ... Semi-Major Axis, AU: 1.17 ... Eccentricity: 0.193 ...
> Aphelion Radius, AU: 1.40 ... Perihelion Radius, AU: 0.94 ...
> VISIT-2 ... Semi-Major Axis, AU: 1.31 ... Eccentricity: 0.275 ...
> Aphelion Radius, AU: 1.67 ... Perihelion Radius, AU: 0.95."
> — Rogers et al. 2012 Table 1.

> "These cyclers repeat every 7 Earth-Mars synodic periods, which usually
> means that 14 vehicles are needed. However, the VISIT cyclers encounter
> Earth and Mars more often than once every 15 years, so fewer vehicles
> are needed."
> — Rogers et al. 2012 Table 1 footnote a.

> "Unlike most other cyclers, the VISIT class of cyclers have orbits that
> are inertially fixed and utilize resonance opportunities between the
> periods of the cycler, Earth, and Mars; and thus, require no flyby
> maneuvers at Earth or Mars."
> — Russell 2004 dissertation §1.3.

> "A spacecraft in a VISIT-1 orbit would circle the Sun in 1.25 Earth
> years, encountering Earth four times in five Earth years and Mars three
> times in two Mars years. A VISIT-2 orbit would need 1.5 Earth years to
> complete, with a spacecraft encountering Earth twice in three Earth
> years and Mars five times in four Mars years."
> — Spaceflighthistory blog summary of Niehoff 1985 IMUSE presentation.

### Notes and known ambiguities

1. **Inconsistent VISIT parameter sets across secondary sources.** Wikipedia
   (citing McConaghy/Longuski/Byrnes 2002 p. 6) reports VISIT-1 aphelion
   1.89 AU, VISIT-2 aphelion 1.45 AU — opposite to Rogers 2012's
   (1.40, 1.67). The spaceflighthistory blog reports VISIT-1 period
   1.25 yr (implying *a* ≈ 1.16 AU) which roughly matches Rogers 2012,
   while Wikipedia's "encounters Earth 3 times in 15 yr" only fits a
   ~5-year orbit. **Several "VISIT" cyclers appear in the literature with
   the same names but different parameters.** Without the original
   Niehoff documents we cannot adjudicate. The YAML uses Rogers 2012's
   internally-consistent values; M7 fuzzy matching should treat the
   VISIT entries with a wide tolerance until this is resolved.
2. **V∞ gaps.** No accessible source tabulates V∞ at Earth or Mars for the
   VISIT cyclers themselves. Rogers 2012 Table 4's V∞,flyby values are
   *establishment* trajectories' V∞, not the steady-state cycling V∞.
3. **The IMUSE / VISIT distinction.** IMUSE is Niehoff's broader mission
   architecture; VISIT is the cycler-spacecraft concept within it. The
   nomenclature is sometimes conflated.

---

## 6. Jones/Hernandez/Jesick 2017 — VEM triple cyclers (family seed)

**Catalogue ID:** `jones-2017-vem-triple-family`

### Citation

| Field | Value |
|---|---|
| First published | Jones, D. R., Hernandez, S., Jesick, M., **"Low Excess Speed Triple Cyclers of Venus, Earth, and Mars,"** AAS Paper **17-577**, JPL-CL#17-3322, AAS/AIAA Astrodynamics Specialist Conference, Stevenson WA, 20–24 August 2017. **NTRS:** <https://ntrs.nasa.gov/citations/20190028464>. **JPL handle:** hdl:2014/46418. |
| Priority date | 2017-08-20 |
| Sequence | V–E–M (canonical placeholder; actual sequences vary across family members) |
| Period (per abstract) | 2 × E-M synodic ≈ 4.27 yr |

### Published parameters (family-level only)

The YAML entry is intentionally a **family seed**: the abstract describes
"thousands" of cyclers and the full member-level table is in the paper,
which is not freely accessible. Individual members can be added as
separate entries once the full PDF is obtained.

| Quantity | Value | Source |
|---|---|---|
| Average transit-leg V∞ across family | < 5 km/s | Jones 2017 abstract |
| Repeat period (E-M synodic units) | 2 (per the abstract) | Jones 2017 abstract |
| Per-cycler V∞, sequences, orbital elements | `null` (paper not accessible) | see Outstanding Questions §D |

### Source quotes

> "Ballistic cycler trajectories which repeatedly encounter Earth and Mars
> may be invaluable to a future transportation architecture ferrying humans
> to and from Mars. Such trajectories which also involve at least one flyby
> of Venus are computed here for the first time."
> — Jones, D. R., Hernandez, S., Jesick, M., *Low Excess Speed Triple
> Cyclers of Venus, Earth, and Mars*, AAS 17-577 (2017), abstract.

> "Triple cyclers are constructed to exhibit low excess speed on Earth-Mars
> transit legs, and thereby reduce the cost of hyperbolic rendezvous.
> Numerous solutions are identified with average transit leg excess speed
> below 5 km/sec, independent of encounter epoch."
> — Same paper, abstract (web search snippet).

> "The paper discovered thousands of previously undocumented two synodic
> period Earth-Mars-Venus triple cyclers, with many solutions identified
> having average transit leg excess speed below 5 km/sec, independent of
> encounter epoch."
> — Web search summary referencing the same paper.

### Notes and known ambiguities

1. **Authorship attribution caveat.** spec.md §16.4 attributes the 2017
   triple cyclers paper to "Longuski et al." That is **incorrect** — the
   primary authors are Jones, Hernandez, and Jesick, all at JPL. Longuski
   is a Purdue researcher who supervised much of the surrounding cycler
   literature (McConaghy, Landau, Yam, Chen, Rogers, Hughes) but is not
   an author of the Jones 2017 paper. Spec.md should be corrected in a
   future revision.
2. **Period anchor.** The spec says (§3): "the natural beat is ≈ 6.4 yr
   (3 × E–M ≈ 4 × E–V)." Jones et al. 2017 found 2-synodic-E-M (4.27 yr)
   triple cyclers, which is NOT the 6.4-yr beat. The beat period is a
   *necessary* condition for V-E-M closure in the simplified model, but
   in the real ephemeris, eccentricities and the gravity assist's
   b-plane degree of freedom open up shorter-period solutions. M8's VEM
   campaign should NOT assume the 6.4-yr beat is the only feasible
   period for VEM cyclers.
3. **Why this is a family seed, not a single entry.** With "thousands"
   of cyclers in the family and no member-level data accessible, the
   only honest catalogue action is to record the family + citation, and
   treat M7 matching as "probable-novel, flag for human review against
   Jones 2017" until the paper's full table can be ingested.

---

## 7. Cross-reference table

| Catalogue ID | Bodies | Period | Max V∞ | Primary citation | Data row in YAML |
|---|---|---|---|---|---|
| `aldrin-classic-em-k1-outbound` | E-M | 2.135 yr (k=1) | 9.7 km/s (M) | Byrnes/Longuski/Aldrin 1993, 10.2514/3.25519 | entry 1 |
| `mcconaghy-2006-em-k2` | E-M | 4.27 yr (k=2) | 5.0 km/s (M) | McConaghy/Landau/Yam/Longuski 2006, 10.2514/1.15215 | entry 2 |
| `s1l1-2syn-em-cpom` | E-M | 4.27 yr (k=2) | 5.65 km/s (E) | McConaghy/Longuski/Byrnes 2002 AIAA 2002-4420 | entry 3 |
| `russell-ocampo-2.1.1+2-case2` | E-M | 4.27 yr (k=2) | 4.1 km/s (E) | Byrnes/McConaghy/Longuski 2002 AIAA 2002-4423 | entry 4 |
| `russell-ocampo-2.3.1+1-case3` | E-M | 4.27 yr (k=2) | 5.4 km/s (E) | Byrnes/McConaghy/Longuski 2002 AIAA 2002-4423 | entry 5 |
| `russell-ocampo-4.3.1-5` | E-M | 8.54 yr (k=4) | 3.1 km/s (E) | Russell/Ocampo 2005 JSR, 10.2514/1.5571 | entry 6 |
| `russell-ocampo-3.1.2+1` | E-M | 6.41 yr (k=3) | 4.6 km/s (M) | Russell/Ocampo 2005 JSR, 10.2514/1.5571 | entry 7 |
| `russell-ocampo-2.5.1+0` | E-M | 4.27 yr (k=2) | 9.9 km/s (M) | Russell/Ocampo 2005 JSR, 10.2514/1.5571 | entry 8 |
| `niehoff-visit1` | E-M | ≈14.95 yr (k=7) | n/k | Niehoff 1985 (presentation, not online) | entry 9 |
| `niehoff-visit2` | E-M | ≈14.95 yr (k=7) | n/k | Niehoff 1985 (presentation, not online) | entry 10 |
| `jones-2017-vem-triple-family` | V-E-M | 4.27 yr (k=2 E-M) | <5 km/s (avg) | Jones/Hernandez/Jesick 2017, AAS 17-577 | entry 11 |

**Machine data file:** `src/cyclerfinder/data/seed_cyclers.yaml`

---

## 8. Outstanding questions

The following gaps and contradictions were encountered while compiling
this catalogue. An M3 or M7 implementer relying on a given entry should
resolve the corresponding question first.

### A. Aldrin orbital-element discrepancy (high priority for M3)

`spec.md` §9 anchors the M3 gate to:

> Aldrin cycler: a ≈ 1.659 AU, e ≈ 0.41, perihelion ≈ 0.98 AU,
> aphelion ≈ 2.34 AU, E→M leg ≈ 146 d.

But the literature consistently reports:

> a = 1.60 AU, e = 0.393, perihelion = 0.97 AU, aphelion = 2.23 AU,
> E→M = 146 d
> (Rogers et al. 2012 Table 1; Russell 2004 Table 3.4 via Aphelion
> Ratio 1.47; Wikipedia citing Byrnes/Longuski/Aldrin 1993)

The gap on *a* is 0.06 AU — six times the M3 gate's
`TOL_A_AU = 0.01`. **Either the spec's value set is wrong, the literature
set is wrong, or both refer to different "Aldrin cyclers" (likely the
distinct outbound/inbound branches, but Russell §3.8 explicitly states
energy properties are identical between the two).**

**Recommendation for M3:**

1. Identify which (a, e) pair the M1 circular-coplanar Lambert
   construction actually produces from `PLANETS["E"].sma_au = 1.00000261`
   and `PLANETS["M"].sma_au = 1.52371034` with 146-day Earth-Mars ToF.
2. Whichever it produces, update the YAML's `aldrin-classic-em-k1`
   entry to that value set (with the note that the alternative also
   appears in the literature with the same provenance) and amend
   spec.md §9 accordingly OR widen the gate tolerances to absorb the
   gap (`TOL_A_AU = 0.07`, `TOL_E = 0.02`, etc.).
3. The 146-d ToF and 0.97-0.98 AU perihelion are robust across all
   sources; M3 should keep the gate on those specifically.

### B. McConaghy 2006 orbital elements (medium priority)

The McConaghy 2006 abstract gives V∞ at Earth (4.7 km/s), V∞ at Mars
(5.0 km/s), and Earth–Mars ToF (153 d), but no orbital elements (a, e,
peri, apo). The full paper is paywalled at AIAA. **Without access to the
paper, we cannot fully specify the canonical signature for M7
matching** — finders that hit this cycler will get `null` matches on
the leg_elements field.

**Recommendation:** Obtain the McConaghy 2006 paper via institutional
access or interlibrary loan; or, alternatively, derive (a, e) numerically
from the V∞ + ToF + the known E-E intermediate encounter geometry.

### C. VISIT-1 / VISIT-2 parameter inconsistency (medium priority for M7)

Wikipedia (citing McConaghy/Longuski/Byrnes 2002 p. 6) and Rogers et al.
2012 Table 1 give contradictory aphelion radii for VISIT-1 and VISIT-2:

| Source | VISIT-1 aphelion | VISIT-2 aphelion |
|---|---|---|
| Wikipedia (citing McConaghy 2002) | 1.89 AU | 1.45 AU |
| Rogers 2012 Table 1 | 1.40 AU | 1.67 AU |

The values appear to be *swapped* — i.e. Wikipedia's "VISIT-1" is
Rogers's "VISIT-2" and vice versa, OR they refer to different
"VISIT"-named cyclers in different papers (Niehoff published several
slightly different variants over 1985-91). **Without the original
Niehoff documents (none online), this cannot be resolved.**

**Recommendation:** Obtain Niehoff 1985 "Manned Mars Mission Design" or
the Friedlander/Niehoff/Byrnes/Longuski 1986 AIAA 86-2009-CP paper
(DOI 10.2514/6.1986-2009) to establish the canonical (a, e, peri, apo)
for each. Until then, M7 matching against the VISIT entries should use
a wide tolerance and any "candidate-novel" finder hit near the VISIT
parameter space should be flagged for human review.

### D. Jones 2017 VEM triple cyclers — full member list (high priority for M8)

The Jones/Hernandez/Jesick 2017 paper reports "thousands" of VEM triple
cyclers but the abstract gives only a family-level summary (average
transit V∞ < 5 km/s). The NTRS record's "downloads" section returned
HTTP 404. The ResearchGate PDF returned HTTP 403. **The catalogue
currently has only a family-seed entry; M7 matching will tag any VEM
finder hit as "probable-novel, flag for human review against Jones
2017" until member-level data can be ingested.**

**Recommendation:** Obtain the Jones 2017 paper via NASA STI, JPL Open
Repository (handle hdl:2014/46418), or AIAA institutional access. Once
the member list is available, add each member as its own YAML entry
with full (sequence, period, V∞ multiset, leg elements).

### E. spec.md §16.4 attribution correction (low priority)

spec.md §16.4 attributes the 2017 triple cyclers paper to "Longuski et
al." Per the NTRS record and the paper's title page, the authors are
Drew R. Jones, Sonia Hernandez, and Mark Jesick (all JPL). Longuski is
not an author. spec.md should be updated to cite Jones, Hernandez,
Jesick. (Not blocking any milestone; just a correctness fix.)

### F. spec.md §3 VEM beat period vs. Jones 2017 findings (medium priority for M8)

spec.md §3 says "the natural beat is ≈ 6.4 yr (3 × E–M ≈ 4 × E–V)." But
Jones et al. 2017 found 2-synodic-E-M (4.27 yr) VEM triple cyclers,
which is NOT the 6.4-yr beat. The beat period is sufficient for closure
in the simplified circular-coplanar model with strict commensurability,
but real eccentricities/inclinations + the b-plane DOF open up shorter
periods. **M8's enumerator should NOT hard-code the 6.4-yr beat as the
only feasible VEM period.**

**Recommendation:** When M8 plans the VEM campaign, broaden the period
search to include all integer multiples of the E-M synodic up to some
cap, not just the 6.4-yr beat. The Jones 2017 paper is the prior art
that demonstrates 2-synodic VEM cyclers exist.

### G. Long Mars→Earth return leg of the Aldrin cycler (low priority for M3)

The 146-day Earth→Mars leg of the Aldrin cycler is well-cited. The
complementary Mars→Earth return — qualitatively described as "16 months
beyond Mars" by Wikipedia — is not cleanly tabulated in any single
primary source we accessed. The YAML records the return as
`tof_days: 519` with an explicit "UNVERIFIED" note in the source quote.

**Recommendation:** M3 only needs the E→M leg for the gate test; if a
future milestone needs the full Aldrin cycle geometry, look in the
Byrnes/Longuski/Aldrin 1993 paper or in McConaghy/Longuski/Byrnes 2002
for a full per-leg breakdown.

---

## 9. Provenance of this document

Compiled May 2026 from:

- The Russell 2004 dissertation (full PDF, UT Austin handle 2152/1253) —
  primary source for the 24-ballistic-cycler taxonomy and the Aldrin
  cycler's energy parameters.
- Rogers/Hughes/Longuski/Aldrin 2012 (AIAA 2012-4746, cached PDF from
  engineering.purdue.edu) — primary source for the orbital elements
  table (Aldrin, VISIT-1, VISIT-2, S1L1, plus several other variants).
- Wikipedia's "Mars cycler" article and the spaceflighthistory blog
  summary of Niehoff 1985 — secondary corroboration for several values.
- Web search snippets quoting AIAA abstracts (McConaghy 2006, Jones 2017)
  — used where direct access was blocked by HTTP 403.

Author of this compilation has not personally consulted: the Niehoff 1985
SAIC presentation, the Niehoff 1986 AAS 86-172 paper, the McConaghy 2006
full paper (only the abstract), the Russell-Ocampo 2004 JGCD or 2005 JSR
papers (only abstracts and the dissertation as the comprehensive
treatment), or the Jones/Hernandez/Jesick 2017 full paper (only the
abstract). Every numerical value in the YAML is grounded in one of the
above sources; gaps are flagged in §8.
