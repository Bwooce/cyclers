# Missing Data — Cycler Catalogue Sourcing Report

*Updated 2026-06-03. Based on analysis of 233 catalogue entries; 788 data_gaps across 207 entries.*

---

## 1. Summary

The catalogue carries **788 data gaps** distributed across 207 entries. The bulk originate from Russell-2004-dissertation entries whose per-leg orbital elements (`a_au`, `e`) and per-leg time-of-flight (`tof_days`) were never fully extracted from the dissertation tables. This makes the **Russell 2004 PhD dissertation** (UT Austin, open-access PDF at `http://hdl.handle.net/2152/1253`; now cached at `docs/refs/russell-2004-dissertation.pdf`) the single highest-leverage sourcing target: a systematic extraction of Tables 3.4, 3.9–3.11, 4.9–4.13 would close the bulk of all gaps. The remaining gaps are:

- **201 "derive" gaps** — fully computable by the project's own multi-rev Lambert solver; no sourcing needed.
- **184 "uncertain/topology-provisional" gaps** — resolved once per-leg elements are confirmed from Russell.
- **U0L1 and "Case 1" steady-state V∞** (McConaghy 2005 dissertation / AIAA 2002-4420 full text) and return-leg ToF for the establishment variants (Rogers 2012 Table 4). NB **S1L1, U0L1, and "Case 1" are three distinct orbits** (Rogers 2012 Table 1: S1L1 a=1.30/e=0.257, Case 1 a=1.22/e=0.238, U0L1 a=2.05/e=0.563).
- **4 "unknown" V∞ gaps** for VISIT-1 and VISIT-2 (Friedlander/Niehoff 1986 AIAA 86-2009-CP; Niehoff 1985/1986 originals never digitised online).

**Resolved this session:** the Hollister & Menning 1970 Earth-Venus periodic-orbit gap (period `k` / years / elements / V∞) is now closed — the single placeholder was individuated into the 15-orbit `hollister-menning-1970-ev-orbit-01..15` family from the PRIMARY paper (now in `docs/refs/`), with V∞ from Table 3 (Vr × 29.785 EMOS) and a shared period 16 yr / k=10 (corrected from a wrong secondary "3.2 yr").

The Russell dissertation PDF (3.45 MB, `russellr74662.pdf`) is **freely downloadable from UT Austin repositories** but is binary-compressed and requires a PDF reader or extraction tool to read the tables — web-fetch cannot decode it. All other primary sources listed below are either paywalled (AIAA, JSR) or inaccessible conference-paper archives; the Rogers 2012 preprint on Purdue's server is the only other freely accessible key source.

---

## 2. Gap Table by Category

| Category | Count | Sourcing status |
|---|---|---|
| **unknown** — needs sourcing | 403 | bulk from Russell 2004 dissertation tables; the rest from Rogers 2012 / McConaghy 2005; VISIT V∞ |
| **derive** — computable by Lambert solver | 201 | No sourcing needed; internal computation |
| **uncertain** — topology-provisional | 184 | Resolved as a by-product of extracting Russell per-leg `a_au` |

| Source family | Gap count | Primary tables needed |
|---|---|---|
| Russell 2004 dissertation | 772 | Tables 3.4, 3.9–3.11, 4.9–4.13 |
| McConaghy 2005 / AIAA 2002-4420 | 7 | Steady-state V∞ for U0L1 and SnLm members |
| Rogers 2012 AIAA 2012-4746 | 6 | Table 4 return-leg ToF for 4:3(2)− and 3:2(1)− variants |
| task#54 backfill items | 3 | Rogers 2012 Tables 3/4 |
| Niehoff 1985/1986 (VISIT V∞) | 4 | Original conference proceedings |
| Hollister–Menning 1970 (period k) | 0 | RESOLVED — primary paper cached in `docs/refs/`, family individuated |

| Structural gap path | Count |
|---|---|
| `trajectory.segments[].a_au` | 216 |
| `trajectory.segments[]` (whole segment) | 198 |
| `trajectory.segments[].tof_days` | 184 |
| `trajectory.segments` (segment list) | 183 |
| `trajectory.maneuvers[].dv_kms` | 4 |
| `trajectory.segments[].n_revs` | 4 |
| `vinf_kms_at_encounters` | 1 |

---

## 3. Per-Source Sourcing Sections

### 3.1 Russell 2004 PhD Dissertation (772 gaps — dominant)

**Confirmed open-access record:**
> Russell, Ryan Paul (2004). *Global Search and Optimization for Free-Return Earth-Mars Cyclers*. Ph.D. dissertation, Department of Aerospace Engineering, University of Texas at Austin.
> Handle: `http://hdl.handle.net/2152/1253`
> Repository item: `https://repositories.lib.utexas.edu/items/1b92593d-d1c5-4183-bdff-319f72fd2f62`
> Direct PDF: `https://repositories.lib.utexas.edu/bitstreams/6920bcd8-7ec8-47b1-9f35-eb9e7eef60c8/download`
> File: `russellr74662.pdf` (3.45 MB)
> **Status: FREELY DOWNLOADABLE — no authentication required**

The dissertation finds 24 strictly-ballistic cyclers with 2–4 synodic periods and 92 ballistic cyclers with 5–6 synodic periods, plus hundreds of near-ballistic cyclers. It uses an E-E free-return construction patched with gravity-assist flybys and employs an "Aphelion Ratio" (AR) and "Turn Ratio" (TR) to classify ballistic feasibility (AR ≥ 1 means the cycler reaches Mars; TR ≥ 1 means all flybys have altitudes ≥ 200 km).

**Tables to extract (confirmed from OUTSTANDING.md research log):**
- **Table 3.4** — the 44-cycler summary for the 2-synodic period family: columns include Russell's `K.J.M+/-N` designator, AR, TR, V∞ at Earth, and per-cycler aphelion. This is the master index for the `russell-ocampo-*` catalogue entries.
- **Tables 3.9–3.11** — per-leg breakdown for the 3-synodic and 4-synodic ballistic families (semi-major axis and tof per segment expected, paralleling Table 3.4 structure).
- **Tables 4.9–4.13** — the Chapter 4 extension covering near-ballistic cyclers and the Chapter 4 family including `4.991gG2` (= S1L1 / McConaghy "Notable"). Table 4.9 row 1 was previously verified to carry aphelion = 1.64 AU and V∞ = 4.99/5.10 km/s for the S1L1 cycler (per OUTSTANDING.md §B resolution).

**What is still needed from Russell Tables 3.4, 3.9–3.11, 4.9–4.13:**
- Per-leg semi-major axis (`a_au`) and eccentricity for every `trajectory.segments[]` entry currently `null`. These are the 216 `trajectory.segments[].a_au` gaps.
- Per-leg time-of-flight (`tof_days`) for 184 segments currently `null`.
- V∞ at Mars for Russell-family entries (currently `null` for most entries — one entry in the VISIT family is noted as a specific gap).

**Note on the published companion journal paper:**
> Russell, R. P., and Ocampo, C. A. (2004). "Systematic Method for Constructing Earth-Mars Cyclers Using Free-Return Trajectories." *Journal of Guidance, Control, and Dynamics*, Vol. 27, No. 3, pp. 321–335.
> DOI: `10.2514/1.1011`
> **Status: PAYWALLED (AIAA). ResearchGate link exists but 403 on direct fetch.**

This companion JGCD paper is a condensed version of the dissertation. For the per-leg per-segment data, the dissertation is the authoritative source (the journal paper condenses the results into summary tables). The dissertation PDF is the more productive target.

---

### 3.2 McConaghy / Longuski / Byrnes — AIAA 2002-4420 (7 gaps)

**AIAA 2002-4420 (conference version, 2002):**
> McConaghy, T. T., Longuski, J. M., and Byrnes, D. V. (2002). "Analysis of a Broad Class of Earth-Mars Cycler Trajectories." *AIAA/AAS Astrodynamics Specialist Conference*, Monterey CA, 5–8 August 2002. AIAA 2002-4420.
> DOI: `10.2514/6.2002-4420`
> NTRS record: `https://ntrs.nasa.gov/citations/20060029711` — **no PDF available ("No Preview Available")**
> JPL repository attempt: `hdl:2014/8886` — returned HTTP 500 at prior ingest and again here.
> **Status: NOT FREELY ACCESSIBLE. AIAA paywalled; NTRS metadata-only; JPL handle broken.**

**JSR 2004 journal version:**
> McConaghy, T. T., Longuski, J. M., and Byrnes, D. V. (2004). "Analysis of a Class of Earth-Mars Cycler Trajectories." *Journal of Spacecraft and Rockets*, Vol. 41, No. 4 (confirmed from search; the OUTSTANDING.md note citing Vol. 41, No. 6 may be a minor error — see note below), pp. 622–628.
> DOI: `10.2514/1.11939`
> **Status: PAYWALLED (AIAA arc.aiaa.org; ResearchGate returns 403).**

*Note on volume/issue discrepancy: search engine snippets and OUTSTANDING.md disagree on Vol. 41, No. 4 vs. No. 6. This is a copy-error risk — verify the issue number when you access the paper.*

**McConaghy 2005 Purdue PhD dissertation:**
> McConaghy, T. T. (2005). *Design and Optimization of Interplanetary Spacecraft Trajectories*. Ph.D. dissertation, School of Aeronautics and Astronautics, Purdue University, West Lafayette, Indiana.
> ProQuest ID 3166673; Purdue e-Pubs identifier AAI3166673.
> URL: `https://docs.lib.purdue.edu/dissertations/AAI3166673/`
> **Status: RESTRICTED — requires Purdue institutional proxy login. No open-access PDF available.**

The dissertation abstract confirms it covers Earth-Mars cycler taxonomy in its second half and introduces the "ballistic S1L1 cycler" explicitly. But no numerical data (V∞, a, e, tof per leg) is available in the abstract.

**Specific gap for `mcconaghy-2005-em-u0l1`:**
- `trajectory.segments[ret-me].tof_days` (return M→E leg ToF) — not in Rogers 2012 or any accessible source.
- `vinf_kms_at_encounters[E].vinf_kms` and `vinf_kms_at_encounters[M].vinf_kms` — steady-state cycling V∞ (not the establishment V∞).
- Source: McConaghy 2005 dissertation or AIAA 2002-4420 — **both inaccessible online**.

**Companion nomenclature paper:**
> McConaghy, T. T., Russell, R. P., and Longuski, J. M. (2005). "Towards a Standard Nomenclature for Earth-Mars Cycler Trajectories." *Journal of Spacecraft and Rockets*, Vol. 42, No. 4, pp. 694–698.
> DOI: `10.2514/1.8123`
> **Status: PAYWALLED (AIAA).**

---

### 3.3 Rogers / Hughes / Longuski / Aldrin 2012 AIAA 2012-4746 (6 gaps)

**Conference preprint (open access):**
> Rogers, B. A., Hughes, K. M., Longuski, J. M., and Aldrin, B. (2012). "Preliminary Analysis of Establishing Cycler Trajectories Between Earth and Mars via V-Infinity Leveraging." *AIAA/AAS Astrodynamics Specialist Conference*, Minneapolis MN, 13–16 August 2012. AIAA 2012-4746.
> DOI: `10.2514/6.2012-4746`
> **Open-access preprint PDF: `https://engineering.purdue.edu/AAC/wp-content/uploads/2012/09/EstablishingCyclerTrajectoriesBetweenEarthAndMarsViaV-InfinityLeveraging-AIAA-2012-47461.pdf`**
> **Status: FREELY ACCESSIBLE (Purdue AAC server). PDF fetched successfully in binary form — needs a PDF reader.**

**Journal version (2015 Acta Astronautica):**
> Rogers, B. A., Hughes, K. M., and Longuski, J. M. (2015). "Establishing Cycler Trajectories Between Earth and Mars." *Acta Astronautica*, Vol. 112, pp. 114–125. (July 2015)
> DOI: `10.1016/j.actaastro.2015.01.023` (ScienceDirect pii S0094576515000831)
> ADS record: `https://ui.adsabs.harvard.edu/abs/2015AcAau.112..114R/abstract`
> **Status: PAYWALLED (ScienceDirect).**

**Specific gaps for `aldrin-4-3-2-establishment` and `aldrin-3-2-1-establishment`:**
Both entries track three data_gaps each (todo_ref `#54-backfill`):
1. `trajectory.segments[ret-me].tof_days` — return M→E leg ToF for the 4:3(2)− and 3:2(1)− establishment timelines. **This is the most actionable gap in the Rogers family.** The open-access preprint PDF at Purdue contains Table 4; a human reading that PDF will find this value directly. A search snippet found the 4:3(2)− outbound ToF = 161 d and confirms V∞,flyby in the 5.51–6.55 km/s range, but the return leg ToF is not visible in snippets.
2. `trajectory.maneuvers[E].dv_kms` — the V∞-leveraging establishment ΔV; in Rogers Tables 3/4.
3. `trajectory.segments[out-em].a_au` — per-segment establishment-epoch (a, e); not tabulated in accessible sources.

**Rogers 2013 PhD dissertation (Blake Rogers):**
> Rogers, B. A. (2013). *Design of Cycler Trajectories and Analysis of Solar Influences on Radioactive Decay Rates During Space Missions*. Ph.D. dissertation, School of Aeronautics and Astronautics, Purdue University.
> Purdue e-Pubs AAI3636504.
> URL: `https://docs.lib.purdue.edu/dissertations/AAI3636504/`
> **Status: RESTRICTED — requires Purdue institutional proxy. No open-access PDF.**

The Rogers 2013 dissertation is a superset of AIAA 2012-4746 and Acta Astronautica 2015, and likely contains the complete establishment-trajectory tables including the return-leg ToF values.

---

### 3.4 Niehoff 1985 / 1986 — VISIT-1 and VISIT-2 V∞ (4 gaps)

The VISIT entries (`niehoff-visit1`, `niehoff-visit2`) have `vinf_kms_at_encounters` null at both Earth and Mars. These are the **steady-state cycling V∞** values, not the establishment V∞.

**Primary sources:**

1. Niehoff, J. (1985). "Manned Mars Mission Design." *Joint AIAA/Planetary Society Conference 'Steps to Mars'*, National Academy of Sciences, Washington DC, July 1985.
   - **Status: NOT ONLINE. No digital archive found. This is an informal conference briefing, not a peer-reviewed paper. The NTRS does not index it.**
   - The spaceflighthistory.blogspot.com post (`http://spaceflighthistory.blogspot.com/2016/01/bridging-gap-between-space-station-and.html`) provides an informal secondary summary confirming VISIT-1 orbital period ~1.25 yr, encounter frequencies, but **no V∞ values**.

2. Niehoff, J. (1986). "Pathways to Mars: New Trajectory Opportunities." *American Astronautical Society*, AAS Paper 86-172.
   - Semantic Scholar record found: `https://www.semanticscholar.org/paper/Pathways-to-Mars:-new-trajectory-opportunities.-Niehoff/b4820fba431adfde41bf54708e5faee4f4daadd6` — **page returned empty content on fetch**. No PDF available.
   - **Status: NOT FREELY ACCESSIBLE. AAS conference proceedings archive not online for 1986.**

3. Friedlander, A. L., Niehoff, J. C., Byrnes, D. V., and Longuski, J. M. (1986). "Circulating Transportation Orbits Between Earth and Mars." *AIAA/AAS Astrodynamics Conference*, Williamsburg VA, August 1986. AIAA 86-2009-CP.
   - DOI: `10.2514/6.1986-2009`
   - AIAA record: `https://arc.aiaa.org/doi/10.2514/6.1986-2009`
   - **Status: PAYWALLED (AIAA; returned HTTP 403 on fetch).**
   - This is the **most promising** source for VISIT V∞ values. The abstract search snippet confirms: "The paper presents specific examples of real trajectory data including flight times, encounter frequency, **hyperbolic velocities**, closest approach distances, and Delta V maneuver requirements." The word "hyperbolic velocities" strongly implies that per-encounter V∞ values for VISIT-1 and VISIT-2 are in this paper.

4. Niehoff, J., Friedlander, A., and McAdams, J. (1991). "Earth-Mars Transport Cycler Concepts." *International Astronautical Congress*, IAF Paper 91-438.
   - **Status: NOT FREELY ACCESSIBLE. IAF/IAC 1991 proceedings not online.**

**CANDIDATE sourcing note:** The abstract snippet for AIAA 86-2009-CP explicitly mentions "hyperbolic velocities" and "encounter frequency" for VISIT-1 and VISIT-2. Any V∞ values extracted from that paper should be presented as CANDIDATE (verify before use) with citation:
> Friedlander, A. L., Niehoff, J. C., Byrnes, D. V., and Longuski, J. M. (1986). "Circulating Transportation Orbits Between Earth and Mars." AIAA 86-2009-CP. DOI: `10.2514/6.1986-2009`.

Access route: AIAA institutional subscription or interlibrary loan.

---

### 3.5 Hollister & Menning 1970 — Earth-Venus periodic orbits (RESOLVED 2026-06-03)

**RESOLVED.** The primary paper was obtained and is cached at
`docs/refs/hollister-menning-1970-periodic-swingby-earth-venus-JSR-7-10.pdf`.
The single placeholder entry was individuated into the 15-orbit family
`hollister-menning-1970-ev-orbit-01..15`, with V∞ taken from the paper's
Table 3 (Vr × 29.785 EMOS) and a shared period of **16 yr / k=10**
(corrected from a wrong secondary "3.2 yr", which was the coplanar
sub-orbit). New V0 data-integrity test: `tests/data/test_hollister_family.py`.

> Hollister, W. M., and Menning, M. D. (1970). "Periodic Swing-By Orbits between Earth and Venus." *Journal of Spacecraft and Rockets*, Vol. 7, No. 10, pp. 1193–1199.
> DOI: `10.2514/3.30134`

The precursor paper:
> Hollister, W. M. (1969). "Periodic Orbits for Interplanetary Flight." *Journal of Spacecraft and Rockets*, Vol. 6, No. 4, pp. 366–369.
> DOI: `10.2514/3.29664`
> **Status: PAYWALLED (AIAA).** No open-access version found.

Access route: AIAA institutional subscription or university library JSTOR/AIAA access. Both papers are 1960s–1970s AIAA Journals; many university libraries have physical or digitised backfile access.

---

### 3.6 Jones / Hernandez / Jesick 2017 — VEM Triple Cyclers (OUTSTANDING.md §D)

This is a `not_two_body` entry (`jones-2017-vem-triple-family`), not in the `missing_vinf` bucket, but recorded here because it is the largest qualitative gap in the catalogue.

> Jones, D. R., Hernandez, S., and Jesick, M. (2017). "Low Excess Speed Triple Cyclers of Venus, Earth, and Mars." *AAS/AIAA Astrodynamics Specialist Conference*, Stevenson WA, 20–24 August 2017. AAS 17-577.
> NTRS: `https://ntrs.nasa.gov/citations/20190028464` — **no PDF download available ("There are no available downloads for this record")**
> JPL handle: `hdl:2014/46418` — **HTTP 500 on fetch (server error)**
> ResearchGate direct PDF link: `https://www.researchgate.net/profile/Drew-Jones-4/publication/321366190.../Low-Excess-Speed-Triple-Cyclers-of-Venus-Earth-and-Mars.pdf` — **HTTP 403 on fetch (authenticated)**

The NTRS abstract confirms: "Ballistic cycler trajectories which repeatedly encounter Earth and Mars ... also involve Venus flybys ... constructed to exhibit low excess speed on Earth-Mars transit legs ... solutions showing average transit leg excess speed below 5 km/s." The paper reportedly tabulates "thousands" of VEM triple cyclers. No individual cycler data can be extracted without the full text.

**Access route:** Contact Drew Jones (JPL), Sonia Hernandez (JPL), or Mark Jesick (Caltech/JPL) directly for a preprint. Alternatively, request via JPL Technical Reports Server or AIAA institutional access.

---

### 3.7 Russell & Strange 2009 — Jovian/Saturnian Moon Cyclers (OUTSTANDING.md §H-adjacent)

> Russell, R. P., and Strange, N. J. (2009). "Cycler Trajectories in Planetary Moon Systems." *Journal of Guidance, Control, and Dynamics*, Vol. 32, No. 1, pp. 143–157.
> DOI: `10.2514/1.36610`
> JPL Open Repository (pre-print): `https://trs.jpl.nasa.gov/handle/2014/40318` → redirects to `hdl:2014/40318` → **HTTP 500 on fetch**
> ResearchGate: `https://www.researchgate.net/publication/245432979_Cycler_Trajectories_in_Planetary_Moon_Systems` — **HTTP 403 on fetch**
> **Status: NOT FREELY ACCESSIBLE via automated fetch; JPL handle currently broken.**

The Russell-Strange 2009 paper covers the Jovian moon (Io, Europa, Ganymede) and Saturnian moon (Titan, Enceladus) cycler families. These entries (`russell-strange-2009-jovian-multimoon-family`, `russell-strange-2009-saturnian-multimoon-family`) are classified `non_heliocentric` and are NOT in the `missing_vinf` or `missing_period` priority buckets. They are noted here only for completeness.

Companion conference paper (likely more accessible):
> Russell, R. P., and Strange, N. J. (2008). "Planetary Moon Cycler Trajectories." *AIAA/AAS Astrodynamics Specialist Conference*, Honolulu HI. AAS 08-222.
> JPL Open Repository: `https://trs.jpl.nasa.gov/handle/2014/40318` (same handle as above, may refer to conference preprint)

---

## 4. Prioritised Human Checklist

The following items are ordered by number of catalogue gaps they would close. Items marked **FREELY ACCESSIBLE** can be done immediately with a PDF reader; others require library access.

---

### Priority 1 — Russell 2004 dissertation, Tables 3.4, 3.9–3.11, 4.9–4.13
**Gaps closed: ~772 (all Russell-family `a_au` and `tof_days` gaps)**
**Status: FREELY DOWNLOADABLE**

1. Download `russellr74662.pdf` (3.45 MB) from:
   `https://repositories.lib.utexas.edu/bitstreams/6920bcd8-7ec8-47b1-9f35-eb9e7eef60c8/download`
2. Open in a PDF reader (not a web fetcher — the binary is compressed).
3. Navigate to Table 3.4 (Chapter 3, 2-synodic ballistic family): extract all rows — columns are Russell designator `K.J.M+/-N`, AR, TR, V∞ at Earth, V∞ at Mars, and per-leg semi-major axis and ToF.
4. Repeat for Tables 3.9, 3.10, 3.11 (3-synodic and 4-synodic families).
5. Repeat for Tables 4.9–4.13 (near-ballistic Chapter 4 families; Table 4.9 row 1 = S1L1 = `4.991gG2`).
6. Cross-reference each row to the `russell-ocampo-*` catalogue entries using the `K.J.M+/-N` designator already in the entry `id` field.
7. Backfill `trajectory.segments[].a_au`, `trajectory.segments[].e`, `trajectory.segments[].tof_days` for each entry. These are sourced values — cite as: "Russell, R. P. (2004). *Global Search and Optimization for Free-Return Earth-Mars Cyclers*. Ph.D. dissertation, UT Austin. Table 3.4 (or 4.9 etc.), row [N]."

**Do NOT** treat "derive"-kind gaps (marked `kind: derive` in data_gaps) as needing sourcing — those are computed from the now-sourced elements.

---

### Priority 2 — Rogers 2012 AIAA preprint, Tables 3–6 (return leg ToF)
**Gaps closed: 6 (aldrin-3-2-1-establishment, aldrin-4-3-2-establishment)**
**Status: FREELY ACCESSIBLE (Purdue server)**

1. Open the Rogers 2012 preprint PDF (already downloaded to the project session as binary — use a PDF reader):
   `https://engineering.purdue.edu/AAC/wp-content/uploads/2012/09/EstablishingCyclerTrajectoriesBetweenEarthAndMarsViaV-InfinityLeveraging-AIAA-2012-47461.pdf`
2. Locate Table 4 (STOUR analytic-ephemeris results).
3. For each row labelled "Aldrin 4:3(2)−" and "Aldrin 3:2(1)−": extract the **total establishment timeline** and derive the return M→E leg ToF as: `T_establishment − ToF_outbound`. Alternatively, if Rogers tabulates return-leg ToF directly, extract that value.
4. Also extract from Tables 3/4: the V∞-leveraging establishment ΔV (`maneuvers[E].dv_kms`) for each variant.
5. Add as CANDIDATE values with citation: "Rogers, B. A., Hughes, K. M., Longuski, J. M., and Aldrin, B. (2012). 'Preliminary Analysis of Establishing Cycler Trajectories Between Earth and Mars.' AIAA 2012-4746. Table 4, row [Aldrin 4:3(2)−]."

---

### Priority 3 — Friedlander et al. 1986 AIAA 86-2009-CP (VISIT V∞)
**Gaps closed: 4 (niehoff-visit1 and niehoff-visit2 V∞ at Earth and Mars)**
**Status: PAYWALLED — requires AIAA institutional access or interlibrary loan**

1. Request AIAA 86-2009-CP via your institution's library AIAA subscription or interlibrary loan.
   DOI: `10.2514/6.1986-2009`
2. Locate the tables of "hyperbolic velocities" (confirmed in abstract snippet) for VISIT-1 and VISIT-2.
3. Extract V∞ at Earth and V∞ at Mars for both cyclers.
4. Backfill `vinf_kms_at_encounters[E].vinf_kms` and `vinf_kms_at_encounters[M].vinf_kms` in both VISIT entries.
5. Add as CANDIDATE values with citation: "Friedlander, A. L., Niehoff, J. C., Byrnes, D. V., and Longuski, J. M. (1986). 'Circulating Transportation Orbits Between Earth and Mars.' AIAA 86-2009-CP. DOI: 10.2514/6.1986-2009. Table [N]."

**Alternative access route:** The Niehoff 1991 IAC paper (IAF 91-438) is a later treatment of the same material and may also contain V∞ tables. No online access found for IAC 1991 proceedings either.

---

### Priority 4 — McConaghy 2005 Purdue dissertation (U0L1 V∞ and SnLm per-member data)
**Gaps closed: 3 (mcconaghy-2005-em-u0l1 return ToF, steady-state V∞ at E and M)**
**Status: RESTRICTED — Purdue proxy required**

1. Access via Purdue institutional proxy: `https://docs.lib.purdue.edu/dissertations/AAI3166673/`
2. Locate the Earth-Mars cycler chapters (second half of dissertation per abstract).
3. Extract per-SnLm-member steady-state cycling V∞ at Earth and Mars (not the establishment V∞) for U0L1, S1L1 (the ballistic version), Cases 1–3.
4. Extract the U0L1 return M→E leg ToF.
5. If S3L1-ballistic and S1L2 members are tabulated, create new catalogue entries following the `mcconaghy-2005-em-<snlm>` naming convention (per the member map in `mcconaghy-2005-em-snlm-broadclass-family`).
6. Add all extracted values as CANDIDATE with citation: "McConaghy, T. T. (2005). *Design and Optimization of Interplanetary Spacecraft Trajectories*. Ph.D. dissertation, Purdue University. ProQuest AAI3166673. Table [N]."

**Alternative if Purdue proxy is unavailable:** Contact Thomas McConaghy directly (LinkedIn: former JPL, Google X) for a copy of the cycler chapters.

---

### Priority 5 — Hollister & Menning 1970 (Earth-Venus periodic orbits) — ✓ RESOLVED 2026-06-03
**Gaps closed: family individuated (period k=10 / 16 yr, V∞ from Table 3, 15 orbits)**
**Status: DONE — primary paper cached in `docs/refs/`**

The Hollister & Menning 1970 JSR paper and the Hollister-Rall 1970 NASA-CR
companion were obtained and cached (`docs/refs/`). The single placeholder was
expanded into `hollister-menning-1970-ev-orbit-01..15` with V∞ from Table 3
(Vr × 29.785 EMOS) and shared period 16 yr / k=10. No further sourcing needed.

---

### Priority 6 — Jones / Hernandez / Jesick 2017 AAS 17-577 (VEM triple cycler member data)
**Gaps closed: 0 current data_gaps (family-seed has no numeric fields), but unlocks M7 matching**
**Status: NOT FREELY ACCESSIBLE — JPL handle broken, NTRS metadata-only**

1. Try the JPL direct URL (check if server is restored):
   `https://hdl.handle.net/2014/46418`
2. If still 500 error, contact Drew Jones (JPL) or Mark Jesick (JPL/Caltech) for a preprint.
3. When the paper is obtained, extract per-member cycler data: sequence (V-E-M or E-M-V-E etc.), period, V∞ at each planet, leg ToF.
4. Add each member as a new YAML entry with id `jones-2017-vem-<sequence>-<n>` following the catalogue naming conventions.

---

## 5. Notes on Sources NOT Requiring Human Sourcing

The following gap types are **computable by the project's multi-rev Lambert solver** and should NOT be sourced:

- `trajectory.segments[].tof_days` where `kind: derive` — these 201 entries are calculated from the `a_au` and encounter geometry once `a_au` is backfilled from Russell.
- `trajectory.segments[].n_revs` where `kind: derive` — same, resolved by the Lambert solver.

These are internal engineering computations; treating them as "sourcing targets" would violate the golden-test discipline (circular: a value cannot be cited to itself). Once Priority 1 provides the sourced `a_au`, run the Lambert solver to fill the derive gaps programmatically.

---

## 6. Source Access Summary

| Source | Open access? | URL |
|---|---|---|
| Russell 2004 dissertation PDF | **YES — free** | `https://repositories.lib.utexas.edu/bitstreams/6920bcd8-7ec8-47b1-9f35-eb9e7eef60c8/download` |
| Rogers 2012 AIAA preprint PDF | **YES — free** | `https://engineering.purdue.edu/AAC/wp-content/uploads/2012/09/EstablishingCyclerTrajectoriesBetweenEarthAndMarsViaV-InfinityLeveraging-AIAA-2012-47461.pdf` |
| Russell & Ocampo 2004 JGCD | NO — AIAA paywall | DOI `10.2514/1.1011` |
| McConaghy et al. 2002 AIAA 2002-4420 | NO — NTRS metadata-only | DOI `10.2514/6.2002-4420` |
| McConaghy et al. 2004 JSR | NO — AIAA paywall | DOI `10.2514/1.11939` |
| McConaghy 2005 Purdue dissertation | NO — Purdue proxy | `https://docs.lib.purdue.edu/dissertations/AAI3166673/` |
| Friedlander et al. 1986 AIAA 86-2009 | NO — AIAA paywall | DOI `10.2514/6.1986-2009` |
| Niehoff 1985 "Steps to Mars" | NO — no digital archive | N/A |
| Niehoff 1986 AAS 86-172 | NO — no digital archive | N/A |
| Hollister & Menning 1970 JSR | NO — AIAA paywall | DOI `10.2514/3.30134` |
| Rogers 2013 Purdue dissertation | NO — Purdue proxy | `https://docs.lib.purdue.edu/dissertations/AAI3636504/` |
| Jones et al. 2017 AAS 17-577 | NO — JPL handle broken | `hdl:2014/46418` |
| Russell & Strange 2009 JGCD | NO — AIAA paywall + JPL broken | DOI `10.2514/1.36610` |
| McConaghy/Russell/Longuski 2005 JSR | NO — AIAA paywall | DOI `10.2514/1.8123` |
| Rogers et al. 2015 Acta Astronautica | NO — ScienceDirect paywall | DOI `10.1016/j.actaastro.2015.01.023` |
