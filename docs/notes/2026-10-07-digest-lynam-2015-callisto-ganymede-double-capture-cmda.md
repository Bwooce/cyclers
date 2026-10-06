# Digest: Lynam 2015/2016, "Broad search for trajectories from Earth to Callisto-Ganymede-JOI double-satellite-aided capture at Jupiter from 2020 to 2060" (CMDA) (#960, #943)

A. E. Lynam (West Virginia University), "Broad search for trajectories from Earth to Callisto-Ganymede-JOI
double-satellite-aided capture at Jupiter from 2020 to 2060", Celestial Mechanics and Dynamical Astronomy 124(1):33-50.
Online 9 September 2015; issue dated January 2016 (Crossref published-print 2016-01). doi 10.1007/s10569-015-9649-9
(printed on p.1; Crossref title, author, volume, issue and pages agree). Received 30 January 2015, revised 2 June 2015,
accepted 28 August 2015.
- File given: `5de2dd22-lynam2015.pdf`, 18 pages, text layer (online-first PDF, no printed page numbers).
  md5 3b51e60b7f2d2d5343d1186a2814594c. I cite PDF pages; journal page = 32 + PDF page.
- **Proposed corpus filename** (year = online year, which is also the year in the DOI and in the file name given):
  `cyclers_pdf/papers/lynam-2015-broad-search-earth-callisto-ganymede-joi-double-satellite-aided-capture-2020-2060-cmda-124-33-doi-10.1007-s10569-015-9649-9.pdf`
  (If the filer prefers the issue year, use `lynam-2016-...`; cite as "Lynam (2016), CMDA 124(1):33-50, online 2015".)
- **Wanted list:** not listed.
- **How I read it.**
  - I read the whole text layer. I read p.5 (phase-angle ranges), p.13 (Fig. 5, the 50-day statement), p.14 (Table 2,
    also at 300 dpi; Table 3), p.15 (Table 4, discussion) and p.2 (Table 1) on renders. Every number quoted from those
    pages was read on the image. I counted the Fig. 5 points on a 300-dpi crop.
  - Arithmetic checks: `checks_lynam.py` secs. 1-3, 5-9, output in `checks_lynam.out`.
  - The #943 collision check is in `collision.md` (written first). The verdict is **no collision**.

## 0. Verdict

**A capture-trajectory paper: Earth launch to one Callisto flyby, one Ganymede flyby and a JOI burn.**
- **#943 collision: no collision** with gc-1, gc-2, ge-1, ge-2 or ge-3. The Callisto-to-Ganymede leg is a single 19-21 h
  transfer on the inbound hyperbola. No moon v_inf is printed. My patched-conic estimates (INFERRED): Ganymede v_inf about
  11.4-13.0 km/s at G1, 7.6-11.7 km/s on the ~205-day orbits and 6.8-9.1 km/s on the post-G2 orbits; Callisto
  7.1-10.1 km/s on the captured orbits. The candidates have 1.4-3.9 km/s.
- **Related (not a cycler):** after JOI and a 200 m/s perijove raise, each trajectory is targeted to a second Ganymede
  flyby (G2) 207.97-208.05 d after the first, which is 29.07 Ganymede periods (my check). This is a targeted one-shot
  re-encounter. It is not periodic.
- **Phasing for X1 (the main value for the project):** consecutive CGJ captures are 50.12-50.13 d apart (Table 2). That is
  4 Ganymede-Callisto synodic periods (50.09 d), which is also 3 Callisto periods and 7 Ganymede periods. This is the
  inertial repeat of the C-G geometry. gc-1 and gc-2 use 3 synodic periods (37.57 d), which repeats only in the rotating
  frame. See `collision.md` sec. 3.
- **Catalogue implication (PROPOSAL only):** none.
- **Lynam 2012 PhD (wanted row 14; the brief says 15, the current list at 52d51c6f has 14):**
  - Post-thesis work (received January 2015, from WVU).
  - It extends thesis secs. 2.6 (double-satellite-aided capture) and 3.4 (targeting double sequences, incl. 3.4.1
    interplanetary windows and 3.4.2 nested backward targeting), and it is the "broad-search algorithm for ballistic
    interplanetary trajectories" of future-work item 6.1 (preview TOC p.vi), for the double case.
  - It does not cover Ch. 4 (navigation; only a qualitative sec. 4.1 here) or Ch. 5 (cyclers). It is not among thesis
    refs [10]-[16].

## 1. Method (READ, pp.3-10)

- **Why CGJ (p.2-3, Table 1, image-read; identical to Table 7 of Lynam 2015 doi 10.1007/s10569-015-9602-y):**
  CGJ needs 568 m/s JOI and 661 m/s JOI + PJR to 14 R_J, with a first perijove of 13 R_J, above most of the radiation.
  That is about 100 m/s more than the triples (CGIJ 563, CGEJ 569) and 250-300 m/s less than single flybys.
- **Phase angles (eqs. 1-5).** Delta-lambda(Ca,Ga) and Delta-lambda(Ca,Sun) (the Sun treated as orbiting Jupiter).
- **Search-space reduction (p.4-5).** Ephemeris-free circular coplanar patched conics. Incoming asymptote rp 6.4-15 R_J,
  V_inf 5.1-6.3 km/s (Hohmann 5.6 km/s). Maximally energy-reducing Callisto flyby at 100 km. Results:
  Delta-lambda(Ca,Ga) must be 8.3-20.3 deg (96.67% of time excluded); Sun-asymptote angle 80-100 deg gives
  Delta-lambda(Ca,Sun) 291.2-352.6 deg (82.94% excluded); combined 99.43% excluded.
- **Windows (eqs. 6-7).** Window start and end from the linear drift of Delta-lambda(Ca,Ga) at n_Ca - n_Ga (negative).
- **Lambert (p.6).** Oldenhuis's MATLAB solver (Izzo's method plus Lancaster-Blanchard-Gooding), compiled; mean
  2.36e-5 s per solve.
- **Backward propagation (pp.6-7).** Callisto B-plane angle 170-190 deg (backward flyby), pre-flyby V_inf by MacDonald &
  McInnes 2005 (eq. 9), Jupiter-centred then heliocentric f and g functions back to perihelion.
- **Perigee finder (pp.7-9, eqs. 19-32).** Newton iteration on heliocentric eccentric anomaly to make the geocentric
  position and velocity perpendicular; 15 fixed iterations; convergence "assumed but not proven"; elliptic only.
- **Grid search (pp.9-11, Figs. 2-4).** Over (rp, V_inf), find the quadrilateral containing the target Earth B-plane
  (B.T = 20,000 km, B.R = 0) and refine up to 10 times to within 1000 km.
- **GMAT (p.12-13).** Five patched-conic solutions re-converged in GMAT. Asymptote epoch, inclination, RAAN and argument of
  perijove target the Callisto and Ganymede B-planes; semi-major axis and eccentricity are hand-tuned to hit Earth.

## 2. Results (READ; image-checked)

- **Patched-conic census (Fig. 5, p.13; conclusions p.17):** 29 trajectories in 8 launch windows, 2017-2057, C3 below
  120 km^2/s^2. My count on a 300-dpi crop: about 2017 (1), 2023 (6), 2030 (2), 2036-37 (6), 2041 (4), 2047-48 (5),
  2054 (3), 2057 (2) = 29 in 8. A January 2017 launch has C3 81 km^2/s^2 (p.12).
- **Table 2 (p.14), five GMAT trajectories, image-read at 300 dpi:**

| | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| launch (2023) | 11 Jul | 13 Jul | 15 Jul | 18 Jul | 20 Jul |
| C3, km^2/s^2 | 85.860 | 85.605 | 85.653 | 85.893 | 86.326 |
| Dec / RA, deg | 1.326 / 33.913 | 1.968 / 35.518 | 2.567 / 37.346 | 3.140 / 39.296 | 3.713 / 41.386 |
| C1 | 10 Feb 2026 14:53 | 1 Apr 2026 18:00 | 21 May 2026 21:01 | 10 Jul 2026 23:58 | 30 Aug 2026 02:47 |
| C1 alt., km | 94.197 | 101.113 | 98.457 | 96.508 | 95.270 |
| G1 | 11 Feb 2026 10:14 | 2 Apr 2026 13:57 | 22 May 2026 17:33 | 11 Jul 2026 20:58 | 31 Aug 2026 00:09 |
| G1 alt., km | 237.303 | 230.156 | 223.079 | 216.624 | 211.409 |
| PJ1 time | 12 Feb 2026 03:42 | 3 Apr 2026 07:41 | 23 May 2026 11:20 | 12 Jul 2026 14:40 | 31 Aug 2026 17:42 |
| PJ1, R_J | 5.730 | 6.618 | 7.390 | 8.029 | 8.519 |
| JOI, m/s | 453 | 476 | 498 | 533 | 575 |
| AJ time | 26 May 2026 01:30 | 15 Jul 2026 06:26 | 3 Sep 2026 08:15 | 23 Oct 2026 11:07 | 12 Dec 2026 14:43 |
| AJ, R_J | 276.187 | 274.987 | 273.835 | 272.823 | 271.998 |
| Period 1, d | 205.829 | 205.310 | 204.789 | 204.297 | 203.975 |
| PJR, m/s | 200 | 200 | 200 | 200 | 200 |
| G2 | 7 Sep 2026 09:45 | 27 Oct 2026 13:17 | 16 Dec 2026 17:11 | 4 Feb 2027 21:18 | 27 Mar 2027 01:27 |
| G2 alt., km | 147.414 | 137.711 | 95.922 | 136.524 | 91.046 |
| PJ2 time | 8 Sep 2026 05:26 | 28 Oct 2026 08:48 | 17 Dec 2026 12:24 | 5 Feb 2027 15:57 | 27 Mar 2027 19:39 |
| PJ2, R_J | 8.538 | 9.276 | 10.001 | 10.722 | 11.341 |
| Period 2, d | 59.413 | 56.370 | 53.025 | 51.301 | 48.382 |

  (The G1 altitude of trajectory 5 is printed without "km".)
- **Table 4 (p.15), against the nominal Europa Clipper GJ capture (Buffington 2014: 500 km Ganymede flyby, 839 m/s JOI,
  122 m/s PJR):** total Delta-V GJ 961 m/s against CGJ 653, 676, 698, 733, 775 m/s; perijove 1: 12.1 against 5.730-8.519 R_J;
  perijove 2: 13.0 against 8.538-11.341 R_J. CGJ saves 186-308 m/s but has lower perijoves.
- **Run time (Table 3, p.14):** 1.6 h in MATLAB; about half in SPICE reads (cspice_spkezr, mice.mex).
- **Navigation (p.16, qualitative):** a +/-20 km Callisto altitude error becomes +/-200 to +/-280 km at Ganymede
  (geometric error growth). Either target Callisto to +/-5 km (3 sigma) with no TCM, or to +/-15 km with a TCM in the
  19-21 h between flybys.
- **Note on the CGEJ alternative (p.15):** "the only potentially plausible CGEJ captures with perijoves inside that range
  occur at 2022, 2050, and 2077" (from Lynam 2015, the triple/quadruple paper).

## 3. Checks (`checks_lynam.out`)

- Search-space percentages: 12/360 gives 96.67%; (352.6 - 291.2)/360 gives 82.94%; combined 99.43%. All agree.
- Table 4 totals = Table 2 JOI + 200 m/s PJR for all five (653, 676, 698, 733, 775). Nominal GJ 839 + 122 = 961. All agree.
- C1 to G1: 19.35, 19.95, 20.53, 21.00, 21.37 h (text "19-21 h", p.16). G1 to PJ1: 17.47-17.78 h (text "15-17 h" for
  the Ganymede-to-JOI interval, p.16). **Small disagreement:** the table gives 17.5-17.8 h, slightly above the text's
  15-17 h. I re-read the G1 and PJ1 rows at 300 dpi; the table values stand. Probably loose wording in the text.
- C1 to C1 steps: 50.130, 50.126, 50.123, 50.117 d. 4 G-C synodic periods = 50.093 d; 3 P_C = 50.067 d; 7 P_G = 50.082 d.
- PJ1 to AJ, doubled: 205.817, 205.896, 205.743, 205.704, 205.751 d. "Period 1": 205.829, 205.310, 204.789, 204.297,
  203.975 d. Trajectory 1 agrees to 0.01 d; trajectories 2-5 differ by 0.6-1.8 d, and the gap grows steadily.
  A two-body period from rp and ra (R_J = 71,492 km, GM 126,686,534 km^3/s^2) gives 205.2-206.7 d.
  I re-read these rows at 300 dpi; the values stand. **Not a misprint, in my view:** "Period 1" is most likely the
  osculating period at one epoch of a GMAT run with solar perturbation, while the half-period PJ-to-AJ time is
  integrated (INFERRED). The paper does not define "Period 1".
- Sec. 4.1 (p.16, image) says "The Callisto flybys in Table 3 are low (95-100 km altitude)". The flyby altitudes are in
  Table 2; Table 3 is the run-time profile. A table-number slip.
- G1 to G2: 207.97-208.05 d = 29.07 Ganymede periods = 16.61 G-C synodic periods. PJ1 to PJ2: 208.04-208.08 d.

## 4. Citation mining (28 references, pp.17-18)

Checked with `ls cyclers_pdf/papers | grep -i` and `grep -i docs/notes/CORPUS_INDEX.md`.

| work | status |
|---|---|
| Lynam 2014 Part I (Acta 94:246-252, doi 10.1016/j.actaastro.2013.07.018) | not held; not on the wanted list. **New candidate (low priority).** |
| Lynam 2014 Part II; Lynam 2015 (CMDA 121:347-363) | in batch 35 (this agent) |
| Lynam, Kloster & Longuski 2011 (CMDA 109); Lynam & Longuski 2011 (JGCD 34); Didion & Lynam 2014 (AIAA 2014-4106) | in batch 35 with a sibling agent |
| Lynam & Longuski 2012 (Acta 70:33-43, navigation) | not held; not on the wanted list. Low priority. |
| Patrick & Lynam 2014 (AIAA 2014-4218) | not held. Low priority. |
| Buffington 2014 (AIAA 2014-4105) | HELD |
| Acton 1996 (NAIF ancillary data services) | HELD (`acton-1996-ancillary-data-services-NASA-NAIF-SPICE-planet-space-sci-doi-10.1016-0032-0633(95)00107-7.pdf`) |
| Russell & Arora 2008, "FIRE: a fast, accurate, and smooth planetary body ephemeris interpolation system", AIAA 2008-6278 | not held (no hit). Tooling; not needed. |
| Jacobson, Haw, McElrath & Antreasian 1999 (Galileo orbit reconstruction); Folkner 2010 (ephemeris uncertainties) | not held (the held Jacobson papers are Uranian ephemerides). Navigation context only. |
| MacDonald & McInnes 2005 (JGCD 28:365); Lancaster & Blanchard 1969 (NASA TN D-5368); Gooding 1990 (CMDA 48:145) | not held. Method sources; not needed. |
| Izzo et al. 2013 (GECCO); Schadegg, Russell & Lantoine (JSR, doi 10.2514/1.A32962); Strange et al. 2012; Landau, Strange & Lam 2010; Longman 1968; Nock & Uphoff 1979; Johannesen & D'Amario 1999; Garrett et al. 2012; Riedel et al. 2006; Hughes 2008; Vallado 2008; Wilson et al. 1997 | not held. None needed for #943. |

- The reference list gives Jacobson et al. 1999 the same paper number (AAS 99-330) as Johannesen & D'Amario 1999. One of the
  two is probably wrong (not checked further; neither matters here).
- New candidates: Part I (low). Nothing for #943.

*Wanted-list row numbers in this digest are those of the list at commit 52d51c6f.*

*Filed as `cyclers_pdf/papers/lynam-2015-earth-callisto-ganymede-joi-double-satellite-aided-capture-2020-2060-cmda-124-33-doi-10.1007-s10569-015-9649-9.pdf`.*

*Wanted-list row numbers in this digest are the batch-34 numbering; the list was renumbered in batch 35.*
