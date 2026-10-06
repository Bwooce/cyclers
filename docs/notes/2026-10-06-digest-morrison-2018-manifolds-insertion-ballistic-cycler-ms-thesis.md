# Digest: Morrison 2018, "Use of Manifolds in the Insertion of Ballistic Cycler Trajectories"

- **Citation:** Morrison, O. (2018). *Use of Manifolds in the Insertion of Ballistic Cycler Trajectories*. MS thesis, Aerospace Engineering, California Polytechnic State University, San Luis Obispo. May 2018. doi 10.15368/theses.2018.80. Committee chair: K. Abercromby.
- **Pages:** 98 PDF pages (i–xi front matter, 87 numbered pages). PDF page = printed page + 11.
- **md5:** `c5f2c1e8db1f40753b52c50ca482fbc0`
- **File:** `cyclers_pdf/papers/morrison-2018-use-of-manifolds-insertion-ballistic-cycler-trajectories-ms-thesis-cal-poly-doi-10.15368-theses.2018.80.pdf`
- **What I read:**
  - The full pdftotext output of all 98 pages.
  - Page images at 150 dpi for every table quoted here: Tables 2.1 (p.17), 4.1 (p.36), 4.2 (p.37), 4.3 (p.40), 4.4 (p.42), 5.1 (p.50), 5.2 and Eq. 5.1 (p.60), 6.1 (p.68), A.1 (p.81) and A.2 (p.82).
  - For cross-checks I also read page images of McConaghy, Landau, Yam & Longuski 2006, JSR 43(2), Tables 2–6 (pp.460–461), from the corpus PDF.
  - Catalogue rows `mcconaghy-2006-em-k2` (L210), `s1l1-2syn-em-cpom` (L484) and `russell-ch4-4.991gG2` (L40799).
  - Notes: `2026-06-17-digest-mcconaghy-2006.md` and `2026-06-19-digest-byrnes-mcconaghy-longuski-2002-two-synodic-cyclers.md`.
- **Helper scripts (scratch, not kept):**
  - `tofcheck.py` checks that each itinerary's dates agree with its printed leg TOFs. All 5 tables agree to within 1 day.
  - `diffmc.py` and `diffmc.out` match the thesis itineraries row by row against McConaghy 2006.
- **Scope:** The thesis has no catalogue-grade new cycler. It is relevant to S1L1 insertion and to the provenance of S1L1 ephemeris itineraries.

## 0. Verdict

- **Cycler:** Only the ballistic **S1L1 cycler** ("S1L1-B") of McConaghy, Longuski & Byrnes 2002. Outbound and inbound vehicles are both used. The Aldrin cycler appears only as a warm-up example (circular model, App. A.1). No other Earth–Mars cycler is used.
- **Cycler model:** There are two models.
  1. A circular-coplanar check (App. A.2). It uses a synodic period of exactly 2 1/7 yr and multi-rev Izzo–Gooding Lambert arcs. It finds τ = 2.8277 yr.
  2. A **DE405 ephemeris, patched-conic** model. It uses Lambert legs, instantaneous V∞ rotation at flybys and V∞ matching ("V∞ surface plots" plus a "timeline" of contour intersections). This is a brute-force search with whole-day steps. It uses no optimiser and no integration.
- **Manifold model:** Sun–Earth **CR3BP**.
  - The hub is a **southern-class halo orbit about Sun–Earth L2**. It is **not L1** and it is **not a Lyapunov orbit**. L1 appears only in the illustration Figs 3.3–3.4 and in Future Work (§8.2).
  - The halo code is adapted from M. Rund (2018, Cal Poly).
  - The leg to the halo is an **invariant stable manifold** branch. A brute-force search finds the branch that comes closest to Earth, and then a parking orbit is chosen.
  - The leg away from the halo is an **unstable manifold** branch, propagated for 6 TU (≈ 334 d).
  - The unstable branch is treated as planar. Earth's orbit is taken as circular until the cycler burn. The branch is rotated about z (angle θz) so that the departure date matches the first leg of the cycler. The search minimises score = 0.7·dc/10³ + 0.3·ΔV (Eq. 5.1, p.60), using MATLAB fminbnd.
  - The halo Az is chosen with fminsearch (Table 5.1).
- **Insertion ΔV results (Table 6.1, p.68, page image):**
  - Parking orbit to halo is 2.75 km/s. This uses the chosen halo Az = 4.197894e5 km, parking-orbit altitude hp = 2413 km and TOF = 31 d.
  - Manifold to cycler is 2.97–3.15 km/s.
  - Each total is the sum of these two: 5.72–5.90 km/s.
  - The thesis compares these with the first Earth V∞ of each cycler, 3.33–4.52 km/s. **The thesis's conclusion:** the manifold route costs more than a direct burn from Earth orbit. It is better only if the cycler vehicle starts in (or is built in) the halo orbit.
- **Caveats (mine):**
  - (a) The "direct-burn" baseline is the bare V∞ (Ch 4.3: "the ΔV value would be 4.01 km/s"). It is not a hyperbolic injection from the parking orbit. With the Oberth effect from the thesis's own 2413 km parking orbit (my computation: μE = 398600 km³/s², RE = 6378 km), direct injection is about 3.60 km/s for V∞ 4.01, 3.36 for 3.33, 3.81 for 4.52 and 3.71 for 4.29. So the thesis's negative verdict on the full route is **stronger** than stated.
  - (b) The thesis itself says the manifold-to-cycler ΔV comes from a velocity-vector difference only, with a position gap dc of 132–841 km. So "deceptive" is the thesis's own word for it.
  - (c) The manifold phase is planar CR3BP with a circular Earth. The cycler is 3-D ephemeris. The halo has a 4.2e5 km out-of-plane amplitude, but the departure is analysed in the x–y plane.
  - (d) Outbound-1 θz is 0.741596 rad in the text (p.61) and 0.7632 rad in Table 6.1. The full-fidelity rerun probably moved it. This is a minor point.
- **Against the catalogue's S1L1 insertion rows:** The thesis says the only existing insertion method is a direct burn (abstract and §1.1). The catalogue already holds sourced S1L1 establishment trajectories that use V∞ leveraging, all with `inserts_into: s1l1-2syn-em-cpom`. They come from Rogers 2015 Tables 3–4, and Rogers et al. 2012 (AIAA 2012-4746) is earlier than this thesis.
  - `s1l1-5-4-3-establishment-cc` (L4241): launch V∞ 2.042, ΔV_DSM 0.258, ΔV_total 3.647 km/s.
  - `s1l1-3-2-2-establishment-cc` (L4321): launch V∞ 3.340, ΔV_DSM 0.160, ΔV_total 3.859 km/s.
  - `s1l1-4-3-2-establishment` (L3505, ephemeris): launch V∞ 2.492 on 2022-12-20, ΔV_DSM 0.182 km/s.
  - I did not check here how Rogers defines ΔV_total. Even so, the leveraging launch V∞ of 2.0–3.3 km/s plus a DSM of 0.16–0.26 km/s is far below the manifold route's 5.72–5.90 km/s from parking orbit. It is also below the 2.97–3.15 km/s halo-to-cycler burn alone. The thesis does not cite Rogers.
- **Most useful content for us:**
  - The thesis prints two **30-year DE405 S1L1 itineraries** that it says came from T. T. McConaghy by personal communication (Table 2.1 outbound, Table A.1 inbound). It also prints three "generated" ones (Tables 4.3, 4.4 and A.2).
  - **All five match McConaghy et al. 2006 JSR Tables 2–5 row for row.** The usual match is 0–1 day, ≤0.04 km/s and ≤200 km. See §1.
  - So they are **not independent** of McConaghy 2006. Do not use them as a second source for any catalogue value, because that would be circular.
  - The "personal communication" tables add **three earlier encounters (2005–2008)** that McConaghy 2006 does not print. They also have a different tail after 2031. These rows are probably from the AAS 03-509 conference data that McConaghy 2006 calls "omitted for brevity". This is a provenance lead only.
- **Two project records corrected in batch 29 (row removed; Ross digest fixed); before that:** `docs/notes/2026-10-05-960-wanted-papers.md:117` (row 59) and `docs/notes/2026-10-06-digest-ross-2021-cycler-quartet-venus-sol-terra-l1.md` (~L73–75). Both describe this thesis as "S1L1 insertion from STL1" or "Sun-Earth L1". The thesis uses a **Sun–Earth L2 halo**.

## 1. Cycler states against the catalogue

### 1.1 What the thesis gives, and what it does not

- **Circular-coplanar S1L1 (App. A.2, p.78–79; also §2.6.1, p.13–14).**
  - Units: AU, with a_e = 149 598 023 km. Times are in years, and S = 2 1/7 yr exactly.
  - The positions are Earth positions at the cycler's Earth encounters, and the Mars encounter position:
    - R0 = a_e[1, 0, 0]
    - R(τ = 2.8277 yr) = a_e[0.4688, −0.8833, 0]. In §2.6.1 the values are τ = 2.8276 and [0.4690, −0.8832].
    - R(4 2/7 yr) = a_e[−0.2225, 0.9749, 0]
    - Rm = a_e[−0.5192, 1.4292, 0]. |Rm| = 1.5206 AU, which is the circular Mars radius. It is **not an aphelion**.
  - The leg structure is in the thesis's own words (Eqs 2.6–2.8, p.13). Leg 1 runs from 0 to τ: E → M (no flyby) → E, 2.8277 yr. It uses the "left branch, long way" multi-rev Lambert solution. Leg 2 runs from τ to 4 2/7 yr: E → E, 1.4580 yr, "right branch, short way".
  - **The thesis prints no circular-model V∞, velocity vector, turn angle, flyby altitude or aphelion.** It also prints no circular E→M transfer time.
- **Ephemeris S1L1 (DE405).**
  - The thesis gives only itineraries: epoch (calendar date), body, V∞ (km/s), "closest approach distance" (km) and leg TOF (days).
  - The "distance" column holds **altitudes**. McConaghy 2006 gives the same numbers under the heading "closest approach altitude, km".
  - **The thesis prints no heliocentric position or velocity vectors, no turn angles and no aphelia.**
  - The whole itinerary spans about 33 yr, which the thesis calls the "Inertial Period" (my date spans: 32.6–32.9 yr).
- **Summary statistics from the itineraries.**
  - Table 4.1 (p.36, image) gives transit times as mean | (min, max) in days:

    | | Outbound | Inbound |
    |---|---|---|
    | E→M | 161 (115, 231) | 868 (788, 938) |
    | M→E | 862 (788, 930) | 159 (108, 225) |
    | E→E | 536 (529, 542) | 536 (529, 542) |

  - Table 4.2 (p.37, image) gives the Mars-minus-Earth heliocentric angle just before the cycler's E→M transfer:
    - Outbound: 4.5392 rad (4.3309, 4.7507) = 260.08° (248.14, 272.20)
    - Inbound: 3.5866 rad (3.4007, 3.6737) = 205.50° (194.84, 210.49)

### 1.2 Itineraries (page images): first rows and ranges

Each table's dates and printed leg TOFs agree to within 1 day (`tofcheck.py`). The full transcriptions are given below in §1.2a. I re-read the 2005–2008 leading rows of Tables 2.1 and A.1 from 300 dpi crops, because no published table covers them.

| Thesis table (page) | Name in Ch 6 | Source as stated | Start (Earth V∞) | Earth V∞ range | Mars V∞ range | Min Earth alt | Min Mars alt |
|---|---|---|---|---|---|---|---|
| 2.1 (p.17) | Outbound-1 | McConaghy, personal comm. | 2005-08-13 (4.01) | 4.01–7.23 | 2.77–7.85 | 2 756 | 1 770 |
| A.1 (p.81) | Inbound-1 | McConaghy, personal comm. | 2005-04-01 (3.33) | 3.33–9.12 | 2.75–7.71 | 10 950 | 678 |
| 4.3 (p.40) | Outbound-2 ("Vehicle 1") | "generated" by author | 2007-09-30 (4.52) | 3.98–11.05 | 3.00–7.87 | 600 | 1 450 |
| 4.4 (p.42) | Inbound-2 ("Vehicle 2") | "generated" by author | 2007-05-29 (4.29) | 3.40–7.16 | 2.76–7.72 | 8 150 | 2 100 |
| A.2 (p.82) | ("Vehicle 3") | "generated" by author | 2005-09-09 (5.32) | 3.76–7.09 | 2.77–7.70 | 9 630 | 7 590 |

- Table 6.1's V∞ column (4.01, 3.33, 4.52, 4.29) equals the first-row Earth V∞ of Tables 2.1, A.1, 4.3 and 4.4. The tables are consistent.
- Table 4.3, encounter 18 (Earth 2031-11-14), prints V∞ "446". This is a typo for 4.46. McConaghy 2006 Table 2 gives 4.43.
- Each thesis table has 24 encounters. McConaghy 2006 has 22 per vehicle. For Table 2.1 this means three extra leading rows plus 21 rows that overlap McC Table 3. McC Table 3 also has one trailing row (Earth 2039-10-26) that the thesis lacks. In A.1 the last leg is M→E of 96 d (2037-09-08 to 2037-12-13) with V∞ 9.12. That is outside the thesis's own inbound range of 108–225 d.

### 1.2a Full transcriptions (from page images)

**Table 2.1** (encounter numbers as printed; alt = closest-approach altitude, km; TOF = leg TOF printed on the arrival row, days)

| # | Body | Date | V∞ km/s | Alt km | Leg TOF d |
|---|---|---|---|---|---|
| 0 | Earth | 2005-08-13 | 4.01 | -- | -- |
| 1 | Mars | 2006-02-27 | 3.02 | 4816 | 198 |
| 2 | Earth | 2008-06-09 | 6.89 | 20130 | 833 |
| 3 | Earth | 2009-12-03 | 6.9 | 31110 | 541 |
| 4 | Mars | 2010-06-06 | 4.31 | 17710 | 186 |
| 5 | Earth | 2012-08-24 | 6.42 | 26490 | 809 |
| 6 | Earth | 2014-02-14 | 6.43 | 41520 | 539 |
| 7 | Mars | 2014-07-03 | 7.14 | 12190 | 138 |
| 8 | Earth | 2016-12-09 | 4.01 | 27730 | 890 |
| 9 | Earth | 2018-05-22 | 4.03 | 19920 | 530 |
| 10 | Mars | 2018-09-15 | 6.47 | 11580 | 115 |
| 11 | Earth | 2021-04-06 | 4.61 | 22990 | 934 |
| 12 | Earth | 2022-09-20 | 4.59 | 14780 | 532 |
| 13 | Mars | 2023-05-01 | 2.77 | 7601 | 223 |
| 14 | Earth | 2025-07-02 | 7.08 | 23860 | 793 |
| 15 | Earth | 2026-12-26 | 7.09 | 35120 | 542 |
| 16 | Mars | 2027-06-14 | 5.27 | 13840 | 170 |
| 17 | Earth | 2029-09-21 | 5.8 | 26850 | 830 |
| 18 | Earth | 2031-03-12 | 5.8 | 37520 | 537 |
| 19 | Mars | 2031-07-15 | 7.85 | 8802 | 125 |
| 20 | Earth | 2034-01-15 | 4.21 | 24870 | 915 |
| 21 | Earth | 2035-06-28 | 4.2 | 2756 | 529 |
| 22 | Mars | 2035-11-12 | 5.87 | 1770 | 137 |
| 23 | Earth | 2038-05-05 | 7.23 | -- | 906 |

**Table A.1** (encounter numbers as printed; alt = closest-approach altitude, km; TOF = leg TOF printed on the arrival row, days)

| # | Body | Date | V∞ km/s | Alt km | Leg TOF d |
|---|---|---|---|---|---|
| 0 | Earth | 2005-04-01 | 3.33 | -- | -- |
| 1 | Mars | 2007-10-05 | 7.25 | 12140 | 918 |
| 2 | Earth | 2008-02-15 | 6.22 | 44630 | 133 |
| 3 | Earth | 2009-08-07 | 6.23 | 19040 | 539 |
| 4 | Mars | 2011-11-05 | 4.78 | 6710 | 820 |
| 5 | Earth | 2012-04-30 | 7.05 | 29830 | 177 |
| 6 | Earth | 2013-10-24 | 7.05 | 24830 | 542 |
| 7 | Mars | 2016-01-08 | 2.75 | 9870 | 805 |
| 8 | Earth | 2016-08-13 | 4.15 | 13410 | 218 |
| 9 | Earth | 2018-01-25 | 4.17 | 22550 | 530 |
| 10 | Mars | 2020-08-20 | 7.19 | 10440 | 938 |
| 11 | Earth | 2020-12-06 | 4.69 | 31100 | 108 |
| 12 | Earth | 2022-05-22 | 4.67 | 16020 | 532 |
| 13 | Mars | 2024-10-16 | 6.64 | 3854 | 877 |
| 14 | Earth | 2025-03-11 | 6.72 | 40750 | 146 |
| 15 | Earth | 2026-09-03 | 6.71 | 22390 | 541 |
| 16 | Mars | 2028-11-11 | 3.79 | 10580 | 800 |
| 17 | Earth | 2029-05-29 | 6.19 | 25610 | 199 |
| 18 | Earth | 2030-11-18 | 6.18 | 30640 | 539 |
| 19 | Mars | 2033-04-07 | 3.77 | 15350 | 870 |
| 20 | Earth | 2033-09-16 | 3.87 | 10950 | 162 |
| 21 | Earth | 2035-02-27 | 3.84 | 24270 | 529 |
| 22 | Mars | 2037-09-08 | 7.71 | 678 | 924 |
| 23 | Earth | 2037-12-13 | 9.12 | -- | 96 |

**Table 4.3** (encounter numbers as printed; alt = closest-approach altitude, km; TOF = leg TOF printed on the arrival row, days)

| # | Body | Date | V∞ km/s | Alt km | Leg TOF d |
|---|---|---|---|---|---|
| 1 | Earth | 2007-09-30 | 4.52 | -- | -- |
| 2 | Mars | 2008-05-19 | 3 | 6600 | 231 |
| 3 | Earth | 2010-07-16 | 7.05 | 25000 | 788 |
| 4 | Earth | 2012-01-09 | 7.06 | 37500 | 542 |
| 5 | Mars | 2012-06-17 | 5.89 | 9800 | 160 |
| 6 | Earth | 2014-10-11 | 5.33 | 25500 | 846 |
| 7 | Earth | 2016-03-29 | 5.31 | 35200 | 535 |
| 8 | Mars | 2016-07-25 | 7.87 | 9600 | 118 |
| 9 | Earth | 2019-02-04 | 3.99 | 22800 | 924 |
| 10 | Earth | 2020-07-17 | 3.98 | 4270 | 530 |
| 11 | Mars | 2020-12-18 | 4.36 | 5150 | 154 |
| 12 | Earth | 2023-05-25 | 6.08 | 20050 | 887 |
| 13 | Earth | 2024-11-13 | 6.1 | 27300 | 538 |
| 14 | Mars | 2025-06-01 | 3.71 | 16100 | 200 |
| 15 | Earth | 2027-08-08 | 6.72 | 25900 | 799 |
| 16 | Earth | 2029-01-30 | 6.73 | 41300 | 541 |
| 17 | Mars | 2029-06-27 | 6.62 | 16700 | 148 |
| 18 | Earth | 2031-11-14 | 446 (sic; = 4.46) | 32000 | 871 |
| 19 | Earth | 2033-04-29 | 4.47 | 24900 | 532 |
| 20 | Mars | 2033-08-18 | 7.58 | 7100 | 111 |
| 21 | Earth | 2036-03-06 | 4.72 | 22300 | 930 |
| 22 | Earth | 2037-08-20 | 4.74 | 600 | 532 |
| 23 | Mars | 2038-01-22 | 5.66 | 1450 | 155 |
| 24 | Earth | 2040-05-11 | 11.05 | -- | 840 |

**Table 4.4** (encounter numbers as printed; alt = closest-approach altitude, km; TOF = leg TOF printed on the arrival row, days)

| # | Body | Date | V∞ km/s | Alt km | Leg TOF d |
|---|---|---|---|---|---|
| 1 | Earth | 2007-05-29 | 4.29 | -- | -- |
| 2 | Mars | 2009-10-17 | 5.83 | 4500 | 872 |
| 3 | Earth | 2010-03-25 | 6.94 | 40600 | 159 |
| 4 | Earth | 2011-09-18 | 6.95 | 23400 | 542 |
| 5 | Mars | 2013-11-17 | 3.25 | 12000 | 791 |
| 6 | Earth | 2014-06-15 | 5.81 | 23700 | 210 |
| 7 | Earth | 2015-12-04 | 5.79 | 18300 | 537 |
| 8 | Mars | 2018-05-27 | 4.6 | 6960 | 905 |
| 9 | Earth | 2018-10-15 | 3.78 | 8150 | 141 |
| 10 | Earth | 2020-03-26 | 3.79 | 26100 | 529 |
| 11 | Mars | 2022-09-28 | 7.72 | 12500 | 916 |
| 12 | Earth | 2023-01-27 | 6.05 | 43300 | 121 |
| 13 | Earth | 2024-07-18 | 6.07 | 11200 | 538 |
| 14 | Mars | 2026-11-03 | 5.57 | 2100 | 838 |
| 15 | Earth | 2027-04-17 | 7.16 | 31000 | 165 |
| 16 | Earth | 2028-10-10 | 7.15 | 24400 | 542 |
| 17 | Mars | 2030-12-08 | 2.76 | 9600 | 788 |
| 18 | Earth | 2031-07-21 | 4.66 | 17300 | 225 |
| 19 | Earth | 2033-01-03 | 4.66 | 21800 | 532 |
| 20 | Mars | 2035-07-26 | 6.34 | 9400 | 934 |
| 21 | Earth | 2035-11-19 | 3.96 | 19000 | 116 |
| 22 | Earth | 2037-05-02 | 3.95 | 26200 | 530 |
| 23 | Mars | 2039-10-14 | 7.08 | 41800 | 895 |
| 24 | Earth | 2040-04-07 | 3.4 | -- | 176 |

**Table A.2** (encounter numbers as printed; alt = closest-approach altitude, km; TOF = leg TOF printed on the arrival row, days)

| # | Body | Date | V∞ km/s | Alt km | Leg TOF d |
|---|---|---|---|---|---|
| 1 | Earth | 2005-09-09 | 5.32 | -- | -- |
| 2 | Mars | 2006-03-03 | 3 | 9960 | 175 |
| 3 | Earth | 2008-06-09 | 6.92 | 20500 | 830 |
| 4 | Earth | 2009-12-03 | 6.93 | 31200 | 542 |
| 5 | Mars | 2010-06-06 | 4.31 | 17800 | 185 |
| 6 | Earth | 2012-08-24 | 6.42 | 26500 | 809 |
| 7 | Earth | 2014-02-14 | 6.43 | 41500 | 540 |
| 8 | Mars | 2014-07-03 | 7.14 | 12200 | 138 |
| 9 | Earth | 2016-12-09 | 4.01 | 27700 | 890 |
| 10 | Earth | 2018-05-22 | 4.03 | 19900 | 530 |
| 11 | Mars | 2018-09-15 | 6.47 | 11600 | 115 |
| 12 | Earth | 2021-04-06 | 4.61 | 23000 | 934 |
| 13 | Earth | 2022-09-20 | 4.59 | 14800 | 532 |
| 14 | Mars | 2023-05-01 | 2.77 | 7590 | 223 |
| 15 | Earth | 2025-07-02 | 7.08 | 23900 | 793 |
| 16 | Earth | 2026-12-26 | 7.09 | 35200 | 542 |
| 17 | Mars | 2027-06-14 | 5.26 | 13800 | 170 |
| 18 | Earth | 2029-09-21 | 5.78 | 26800 | 830 |
| 19 | Earth | 2031-03-12 | 5.78 | 39000 | 537 |
| 20 | Mars | 2031-07-15 | 7.7 | 10600 | 125 |
| 21 | Earth | 2034-01-15 | 3.78 | 23000 | 915 |
| 22 | Earth | 2035-06-28 | 3.76 | 9630 | 529 |
| 23 | Mars | 2035-11-13 | 4.68 | 15700 | 138 |
| 24 | Earth | 2038-05-08 | 5.54 | -- | 907 |

### 1.3 Provenance: row-for-row match with McConaghy et al. 2006 JSR Tables 2–5

From `diffmc.out`. The McConaghy values are from page images of JSR 43(2) pp.460–461.

| Thesis table | Matches McC 2006 | Rows matched closely | Typical difference | Where it differs |
|---|---|---|---|---|
| 2.1 (pers. comm., outbound) | Table 3 (vehicle 2, outbound) | McC Earth-1 2009-12-03 = thesis row 3. Matches to Earth 2031-03-12 (16 rows) | 0–1 d, ≤0.03 km/s, ≤140 km | Thesis rows 0–2 (Earth 2005-08-13, Mars 2006-02-27, Earth 2008-06-09) are not in McC 2006. From Mars 2031-07-15 the tail differs: ΔV∞ +0.16 to +1.70 km/s and Mars altitude 1 770 vs 16 200 km. |
| A.2 ("generated" Vehicle 3) | Table 3 | All 21 overlapping rows, tail included | 0–1 d, ≤0.06 km/s, ≤500 km | Only rows 0–1 differ (launch 2005-09-09, V∞ 5.32). Row 2 = Earth 2008-06-09, as in Table 2.1. The tail of A.2 matches the published McC Table 3, but the tail of its seed (Table 2.1) does not. So Table 2.1 is probably an earlier or different solution from the one McConaghy later published. |
| 4.3 ("generated" Outbound-2) | Table 2 (vehicle 1, outbound) | Rows 2–17 | 0–1 d, ≤0.04 km/s, ≤200 km | Launch 2007-09-30 at 4.52 km/s vs McC 2007-10-14 at 5.19. First Mars flyby altitude 6 600 vs 11 200 km. Tail drifts from 2033: up to −7 d and +0.49 km/s, and Earth altitude 600 vs 8 900 km. |
| 4.4 ("generated" Inbound-2) | Table 4 (vehicle 3, inbound) | Rows 2–21 | 0–1 d, ≤0.04 km/s | Launch 2007-05-29 at 4.29 km/s vs McC 2007-06-07 at 5.80 **plus a 0.17 km/s DSM on 2006-04-04**. First Mars flyby 2009-10-17 at 4 500 km vs 2009-10-23 at 300 km. |
| A.1 (pers. comm., inbound) | Table 5 (vehicle 4, inbound) | McC Earth-1 2009-08-07 = thesis row 3. Matches to Mars 2028-11-11 (14 rows) | 0–1 d, ≤0.03 km/s | Rows 0–2 (Earth 2005-04-01, Mars 2007-10-05, Earth 2008-02-15) are not in McC 2006. From 2029 the tail differs: ΔV∞ up to −0.56, dates up to −27 d, and the last Mars altitude is 678 km. |

What the comparison shows:

- The author says his "generated" itineraries come from his own search (§4.2). But he bounded that search with transit-time and phase-angle windows (Tables 4.1 and 4.2) taken from the McConaghy itineraries. So the "generated" itineraries converge back onto McConaghy's published solutions. They are a reproduction with day-step resolution, not an independent solution family.
- The thesis does **not** cite McConaghy 2006 JSR.
- Bibliography item [14] gives the authors of "Analysis of Various Two Synodic Period Earth-Mars Cycler Trajectories" (2002) as "McConaghy and Longuski". The project holds this paper as **Byrnes**, McConaghy & Longuski 2002. This is a small attribution slip.
- §2.6.2 credits the SNOPT optimisation to [14]. In the corpus, SNOPT appears in McConaghy 2006 (digest, p.458–463).
- The personal-communication rows for 2005–2008 come before McC 2006's Earth-1 dates. The 2008-06-09 Earth encounter in Tables 2.1 and A.2 is the start date of the McConaghy 2004 JSR Table 6 S1L1 (see the McC 2006 digest). This suggests the rows come from the AAS 03-509 data set (McC 2006 ref 19).

### 1.4 Field-by-field comparison with the catalogue

**Row `mcconaghy-2006-em-k2` (L210).** This is the correct row for the ephemeris itineraries.

| Field (line) | Catalogue | Thesis | Agreement |
|---|---|---|---|
| model_assumption (L214) | circular-coplanar | DE405 ephemeris itineraries, plus a circular check in App. A.2 | Different model. The thesis gives ephemeris data only. |
| vinf E (L275) | 4.7 km/s (McC 2006 abstract, circular) | Ephemeris Earth V∞ 3.33–9.12. The McC-matched rows give 3.75–7.16. No circular value is printed. | The circular value is not tested. 4.7 lies inside the ephemeris range. |
| vinf M (L278) | 5.0 km/s | Ephemeris Mars V∞ 2.75–7.87. No circular value. | Inside the range only. |
| transit_times_days (L226) | [153, 153] | Outbound E→M mean 161 d (115–231). Inbound M→E mean 159 d (108–225). | The ephemeris means are 6–8 d longer. McC 2006 p.461 gives an E→M range of 113–223 d. |
| period.years (L271) | 4.27 (2 × 2.135) | 4 2/7 = 4.2857 yr (circular, S = 2 1/7 yr). Mean ephemeris cycle (E→M + M→E + E→E): outbound 161 + 862 + 536 = 1559 d = 4.268 yr. | Ephemeris mean agrees to 0.002 yr. The circular check uses a synodic period 0.0079 yr longer. |
| loop-ee tof_days (L351) | 533.7 d (Russell g 1.4612 yr) | Circular E→E leg 4.2857 − 2.8277 = 1.4580 yr = 532.5 d. Ephemeris E→E mean 536 d (529–542). | Circular: −1.2 d. Ephemeris mean: +2.3 d. Range of 13 d. |
| free_return_arcs G tof (L235) | 2.8096 yr = 1026.2 d | Circular leg 1 τ = 2.8277 yr = 1032.8 d | +6.6 d. Most of the gap comes from the different synodic period (4.2857 vs 4.2708 yr). As a fraction of the period: 0.6598 vs 0.6579. |
| aphelion_au (L284) | 1.64 | Not given | n/a |
| dv_band (L238) | essentially_ballistic | The thesis calls the S1L1 "ballistic" (§2.6). It does not report maintenance ΔV. The inbound DSM in McC Table 4 is missing from its own Table 4.4. | Neither confirms nor contradicts. |
| flyby_altitudes_km (L240–260, computed-m7) | M7-computed, not sourced | Thesis/McC altitudes, e.g. Earth 31 110, Mars 17 710 km (2009–2010) | Different realisations. Do not compare node by node. |

**Row `russell-ch4-4.991gG2` (L40799).**
- The leg times L40820 and L40824 are 1.4612 and 2.8096 yr. The thesis's circular values differ as shown above: −1.2 d and +6.6 d.
- V∞ at L40840 and L40843 is 4.99 and 5.10. The thesis has no circular V∞ to compare.
- The transit times at L40815 are [150, 150]. The thesis has no circular transit time. Its ephemeris means are 161 and 159 d.

**Rows `s1l1-*-establishment*` (L3505, L4241, L4321).** These are S1L1 insertion rows. See §0 for the ΔV comparison. The thesis gives no leveraging data, so there is no field to compare.

**Row `s1l1-2syn-em-cpom` (L484).**
- V∞ 5.65 and 3.05 (L510, L513) fall **inside** the thesis's ephemeris envelopes (E 3.33–9.12, M 2.75–7.87). This is only an envelope, the same kind of support as Byrnes 2002 Case 3. It gives **no literal support** for that row.
- E→M 154 d (L542) is close to the thesis's 161 d ephemeris mean. The thesis cites neither 154 d nor 153 d.

**Note on sources.** The thesis takes its S1L1 from:
- McConaghy, Longuski & Byrnes 2002 [13]
- "McConaghy & Longuski 2002" [14] (in fact Byrnes-led)
- McConaghy, personal communication [15]
- Russell & Ocampo [18], for background only

It does **not** use Rogers 2012 or McConaghy 2006 JSR.

## 2. Content by chapter

- **Ch 1 Introduction (p.1–2).**
  - Cyclers need large vehicles. Insertion has only been studied as a direct (Hohmann-like) burn.
  - Proposed route: parking orbit → Sun–Earth L2 halo → unstable manifold → burn onto the first leg of the cycler.
  - Scope: S1L1-B only, and Sun–Earth L2 only.
- **Ch 2 Cycler trajectory examination (p.3–18).**
  - Background on Aldrin, Russell–Ocampo and the Purdue work.
  - Flyby turn angle: Eq 2.1 δ = 2 asin(μ/(μ + rp v∞²)) and Eq 2.2, the atan2 form. Ballistic condition, Eq 2.3: |v∞,out| − |v∞,in| = 0.
  - Lambert solver: Izzo–Gooding with Lancaster–Blanchard elements [11]. It handles multi-rev, short/long way and left/right branch (Fig 2.2).
  - Nomenclature: nPr, P1r1P2r2, outbound/inbound, and ballistic/powered.
  - Aldrin cycler in the circular model (Eqs 2.4–2.5, Fig 2.3): 147 d to Mars, 633 d back, 779 d in total.
  - S1L1 circular solution (Eqs 2.6–2.8, τ = 2.8276 yr, Fig 2.4).
  - Ephemeris assumptions: JPL DE405, flybys at Earth and Mars, instantaneous rotation.
  - Table 2.1 is McConaghy's outbound itinerary. Fig 2.5 shows the first two legs.
- **Ch 3 Manifold dynamics (p.19–28).**
  - CR3BP in canonical units: Eqs 3.1–3.5 for μ*, DU, MU, TU and VU, and Eqs 3.6–3.10 for the equations of motion. Libration points: Eqs 3.11–3.12.
  - Halo orbits: northern/southern classes and Az. The method is from Rund [17]: analytic first guess plus a differential correction.
  - Stable and unstable invariant manifolds (Fig 3.4, L1 example).
  - No new equations.
- **Ch 4 Cyclers (p.29–44).**
  - Verification method: pork-chop plots (Figs 4.1–4.2) are merged into "V∞ surface plots" (Fig 4.3). The intersection contour gives the dates on which a ballistic flyby works. All contours are drawn on a "timeline" with x = y (Fig 4.4, App. A.4). The closest point is used when a contour misses the line by less than a day.
  - The verification reproduces Table 2.1 "with some admitted variation". No numbers are given.
  - Generation: transit-time windows (Table 4.1) and the E–M phase angle φp = ω + θ − 2πm (Eq 4.1, Table 4.2) bound the search. An 8-step brute-force procedure (p.38) gives Tables 4.3 and 4.4 and Figs 4.5–4.6.
  - Insertion baseline: ΔV is taken as the first V∞, 4.01 km/s (§4.3).
- **Ch 5 Manifolds (p.45–65).**
  - §5.1, halo insertion. A 7-step search finds the stable-manifold branch closest to Earth and its parking orbit. Example: Az = 476 083 km gives hp = 330.5 km, TOF = 197.0 d and ΔV = 3.1629 km/s (p.46).
  - §5.1.1, fminsearch over Az. Table 5.1 (p.50, image):

    | Az0 (1e5 km) | Az (1e5 km) | ΔV (km/s) | hp (km) | TOF (d) |
    |---|---|---|---|---|
    | 3 | 2.999707 | 3.18390 | 235 | 197 |
    | 4 | 4.197894 | 2.75254 | 2413 | 31 |
    | 5 | 5.183945 | 3.08510 | 353 | 169 |
    | 6 | 6.300286 | 2.31501 | 5799 | 738 |
    | 7 | 7.746880 | 1.58907 | 11042 | 596 |

  - Constraints: TOF < 1 yr, and the parking orbit no higher than low MEO. The Az ≈ 4.2e5 km case is chosen. A refined run gives Az = 4.197279e5, ΔV = 2.7513, TOF = 172 d and hp = 2430. It is rejected because the TOF grows.
  - The direct burn from the same parking orbit is given as 4.01–4.52 km/s, mean 4.26 km/s (p.52). These numbers equal the first-row V∞ of Outbound-1 (4.01) and Outbound-2 (4.52), and 4.26 = (4.01 + 4.52)/2. So the burn seems to be taken as the bare V∞, with no Oberth correction. Inbound-1's lower V∞ of 3.33 is not in this range, so only the two outbound vehicles appear to be used.
  - §5.2, departure. Earth's orbit is "circular enough" (Fig 5.7). The unstable manifold is propagated for 6 TU (334 d, inside the 1-yr limit after the 31-d leg to the halo). The shape is fixed, so the departure date becomes a z-rotation θz.
  - Score (Eq 5.1): score = w1·dc/10³ + w2·ΔV, with w1 = 0.7 and w2 = 0.3.
  - Table 5.2 (p.60, image) is a θz scan:

    | θz | 0 | π/3 | 2π/3 | π | 4π/3 | 5π/3 |
    |---|---|---|---|---|---|---|
    | ΔV (km/s) | 3.38 | 3.63 | 3.64 | 3.57 | 3.62 | 3.62 |
    | dc (1e4 km) | 0.57 | 3.29 | 4.01 | 4.05 | 3.55 | 3.75 |
    | score | 1.41 | 3.39 | 3.90 | 3.97 | 3.57 | 3.71 |

  - fminbnd over (0, 2π) and (−π, π) gives θz = 0.741596 rad, score 1.0639, dc = 1500 km and ΔV = 3.51 km/s. A full-resolution rerun gives score 0.91, dc = 132 km and ΔV = 3.02 km/s.
  - Total 5.77 km/s against the 4.01 baseline (p.63). Figs 5.13–5.14.
- **Ch 6 Test cases (p.66–70).**
  - Table 6.1 (p.68, image):

    | Vehicle | θz (rad) | ΔV (km/s) | ΔV total (km/s) | V∞ (km/s) | dc (km) | TOF (d) |
    |---|---|---|---|---|---|---|
    | Outbound-1 | 0.7632 | 3.02 | 5.77 | 4.01 | 132 | 298 |
    | Inbound-1 | −0.4909 | 3.15 | 5.90 | 3.33 | 286 | 248 |
    | Outbound-2 | 3.0434 | 2.97 | 5.72 | 4.52 | 841 | 265 |
    | Inbound-2 | 5.4978 | 3.09 | 5.84 | 4.29 | 137 | 332 |

    Each ΔV total equals 2.75 + ΔV. I checked this for all four rows.
  - Results for inbound and outbound vehicles are very similar. Only a start in the halo makes the route worth using.
- **Ch 7 Conclusions (p.71–72).** The route is feasible. It does not save ΔV from Earth orbit, but it does save ΔV if the vehicle is assembled in the halo.
- **Ch 8 Future work (p.73–74).**
  - Choose the best Earth-departure leg.
  - Find a faster way to generate the cycler.
  - Try other cycler types.
  - Try Sun–Earth L1 and the Earth–Moon libration points. Try halo-to-halo transfers from Earth to Mars.
  - Optimise the stable and unstable manifolds jointly.
  - Use a higher-fidelity model.
- **App A (p.77–84):**
  - Aldrin computation. Results: R(2 1/7) = a_e[0.6235, 0.7818, 0] and Rm = a_e[−1.0163, 0.2505, 0].
  - S1L1 circular computation (see §1.1). Earth velocity is taken from the vis-viva equation. The ballistic τ is found by minimising |v∞1 − v∞2|.
  - Tables A.1 and A.2.
  - Timeline examples (Figs A.3–A.4).
- **App B (p.85–87):** Plots of the test cases in the synodic frame.

## 3. Citation mining

The bibliography has 20 items (p.75–76). Marks:
- ★★ = directly relevant (cycler or manifold insertion)
- ★ = background
- (held) = in CORPUS_INDEX or `cyclers_pdf`, checked by grep

| # | Reference | Relevance |
|---|---|---|
| [1] | Aldrin Mars Cycler, buzzaldrin.com web page (accessed 2018-03-20) | ★ history |
| [2] | Cassini Legacy: gravity assists, JPL web page | — |
| [3] | "Circular Restricted Three-Body Problem", *Interplanetary Mission Design* (course notes; no details given) | — |
| [4] | "Grand Theft Pluto", NASA science news 2007 | — |
| [5] | "Gravity assist", Planetary Society blog 2013 | — |
| [6] | ISEE-3, NSSDC catalogue page | — |
| [7] | JPL Planetary and Lunar Ephemerides (DE405) | ★ model |
| [8] | SNOPT product page | — |
| [9] | B. Aldrin, "Cyclic Trajectory Concepts", Interplanetary Rapid Transit Study, Aerospace Systems Group (SAIC), 1985 | ★★ cycler origin. Check whether held: the grep found no "Cyclic Trajectory" entry in CORPUS_INDEX. |
| [10] | Gómez, Masdemont & Lo (eds), *Libration Point Orbits and Applications*, World Scientific 2003 | ★ manifolds. "Gomez" appears in CORPUS_INDEX; the exact volume is not checked. |
| [11] | D. Izzo, "Revisiting Lambert's Problem", ESA tech. report 2014 (CMDA 2015) | ★ tool. Not in CORPUS_INDEX. |
| [12] | printed as "B. V. Jonathan Brown, Jeremy Peterson and W. Yu", "Seasonal Variations of the James Webb Space Telescope Orbital Dynamics", NASA GSFC 2015 | ★ Sun–Earth L2 transfers |
| [13] | McConaghy, Longuski & Byrnes, "Analysis of a Broad Class of Earth-Mars Cycler Trajectories", AIAA 2002-4420 | ★★ (held) S1L1 primary source |
| [14] | Cited as "McConaghy & Longuski, Analysis of Various Two Synodic Period Earth-Mars Cycler Trajectories", 2002. Lead author is in fact Byrnes. | ★★ (held, as byrnes-mcconaghy-longuski-2002) |
| [15] | T. T. McConaghy, personal communication (itineraries) | ★★ Not obtainable. The data match McC 2006 JSR Tables 3 and 5 plus extra 2005–2008 rows (§1.3). |
| [16] | M. R. Roussel, "Invariant Manifolds", technical report 2005 (as printed) | — |
| [17] | M. Rund, "Numerical Halo Orbit Computation", Cal Poly, March 2018 (senior project or report) | ★ halo/manifold code that this thesis uses. Grey literature, probably on Cal Poly DigitalCommons. |
| [18] | Russell & Ocampo, "A Systematic Method for Constructing Earth-Mars Cyclers Using Direct Return Trajectories", 2004 (as printed: "Technical report, UT Austin"; the journal version is J. Guidance, Control, and Dynamics 27(3), 2004; the corpus holds the AAS 03-145 conference version) | ★★ (held) |
| [19] | N. Shupe, "Earth-Mars Cycler Trajectories", ASEN 5050 research paper (CU Boulder course paper) | ★ grey literature, low value |
| [20] | V. Szebehely, *Theory of Orbits*, Academic Press 1967 | ★ (held, many entries) |

**Leads that this thesis does not cite but that the comparison points to:**
- McConaghy, Yam, Landau & Longuski, AAS 03-509 (2003), "Two-Synodic Period Earth-Mars Cyclers with Intermediate Earth Encounter". This is the likely source of the 2005–2008 itinerary rows (McC 2006 ref 19).
- McConaghy et al. 2006 JSR 43(2) (held). The thesis's itineraries reproduce it.

**Forward-citation note:** Ross (c. 2021) cites this thesis. See `2026-10-06-digest-ross-2021-…`, where the batch-16 entry said "STL1". That entry is now corrected to Sun–Earth L2.

**Wanted-list action:** none. No reference here is new and high-value except AAS 03-509, which is already known through the McC 2006 digest.
