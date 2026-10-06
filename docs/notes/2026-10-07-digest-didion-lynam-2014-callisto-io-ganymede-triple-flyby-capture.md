# Digest: Didion & Lynam 2014, "Impulsive Trajectories from Earth to Callisto-Io-Ganymede Triple Flyby Capture at Jupiter" (AIAA 2014-4106) (#960)

A. M. Didion & A. E. Lynam (West Virginia University), "Impulsive Trajectories from Earth to Callisto-Io-Ganymede
Triple Flyby Capture at Jupiter", AIAA/AAS Astrodynamics Specialist Conference (part of AIAA SPACE 2014, San Diego,
Aug. 2014), AIAA Paper 2014-4106, doi 10.2514/6.2014-4106.
- **Identifier check.** The PDF carries no paper number, no DOI and no conference header (Word metadata title
  "Preparation of Papers for AIAA Technical Conferences"; an author copy). The acknowledgement (p.9) names the
  "Space and Astronautics Forum and Exposition 2014 conference in San Diego". Crossref: doi 10.2514/6.2014-4106,
  same title, "AIAA/AAS Astrodynamics Specialist Conference", published 1 Aug 2014, **author order Lynam, Didion**.
  The PDF lists **Didion first**. I give the PDF order and note the Crossref order.
- **Journal form (not given):** Didion & Lynam (2015), "Impulsive Trajectories from Earth to Callisto-Io-Ganymede
  Triple Flyby Jovian Capture", JSR 52(3):746-753, doi 10.2514/1.A33159 (Crossref: Didion first). Also Didion's WVU
  thesis "Mission Design, Guidance, and Navigation of a Callisto-Io-Ganymede Triple Flyby Jovian Capture",
  doi 10.33915/etd.5495.
- File given: `26289618-ImpulsiveTrajectories2014.pdf`, 10 pages, text layer (MS Word 2007). md5 5500d8e1392f62eb9bfe1e92c18ec7c8.
- **Proposed corpus filename:** `cyclers_pdf/papers/didion-lynam-2014-impulsive-trajectories-earth-callisto-io-ganymede-triple-flyby-capture-jupiter-aiaa-2014-4106-doi-10.2514-6.2014-4106.pdf`
- **How I read it.**
  - Full text layer read. Pages 1-4, 7, 8, 9 rendered at 200 dpi; Tables 1-2 (p.3), 5-8 (p.7), 11-12 (p.9) read on the
    image. Tables 3-4, 9-10 (p.4, p.8) from the text layer, checked against the renders of p.4 and p.8.
  - Arithmetic: `checks_lynam.py` / `checks_lynam.out` ("Didion & Lynam 2014" block).
  - #943 collision check: `collision.md`. Verdict: **no collision** for any candidate.

## 0. Verdict

**One worked Earth-to-Jupiter mission: launch 25 Jun 2022, one broken-plane manoeuvre, then Callisto, Io, JOI and
Ganymede flybys on 6-8 Feb 2025, into a 198.9-d orbit. Total deterministic dV 275.8 m/s. Built in GMAT with an
fmincon optimiser.**
- No cycler; each moon is met once. **#943: no collision** (moon-relative speeds on this arrival are about 12.6 km/s at
  Callisto and 14.5 km/s at Ganymede by our patched-conic estimate; the candidates are 1.8-3.9 km/s).
- **What it gives the project:** nothing for the catalogue. It is a citation leaf in the Lynam capture line. It does
  show that GMAT with Galilean moons added from SPICE handles this kind of multi-flyby targeting.
- **"Novel" claim.** The paper calls the Callisto-Io-perijove-Ganymede sequence "novel" and "one of the other six
  permutations" not found by Lynam [25-26]. The CMDA paper (Table 5, p.19) already lists **CIJG** (Callisto-Io-JOI-
  Ganymede) in patched conics: 211/236/248/236/191 m/s at R_p 5/4/3/2/1.01 R_J. So the novelty is the integrated
  Earth-to-Jupiter trajectory, not the sequence. The achieved JOI, 264 m/s at 3.237 R_J (with 55-279 km flybys and
  a plane change at Io), is close to the patched-conic 236-248 m/s.
- **Catalogue implication (PROPOSAL only):** none.
- **Lynam 2012 PhD (wanted row 14 in the current file; the brief says 15):** none of the thesis. This paper comes after
  the thesis. It carries on the thesis's Ch. 6.1 future work ("Broad-search algorithm for ballistic interplanetary
  trajectories for multiple-satellite-aided capture") through Lynam's 2014 Acta papers, and the Ch. 3 problem
  (Earth-to-Jupiter legs for triple captures), now with an impulsive broken-plane manoeuvre and GMAT instead of STK.

## 1. Content (READ)

- **Setup (sec. II).** GMAT with the four Galilean moons added (JPL SSD gravitational parameters and radii, SPICE
  ephemerides); Prince-Dormand 7(8) integrator; a "JupiterOnly" fallback when a trial passes below a moon's surface.
  Point masses: Sun, planets and Pluto outside Jupiter's SOI; Sun, Jupiter and the four moons inside.
- **Candidate selection.** About 20,000 "theoretically feasible capture orbits" for 2024-2028 (made "using a similar
  method to that of Lynam [25-26]"), cut to about 600 by Callisto B-plane angle, interplanetary dV and R_p; the
  highest-R_p one was chosen.
- **Targets (Tables 3-4).** Callisto 100 km, theta 0; Io 300 km, theta 0; JOI 0.230 km/s; Ganymede 100 km, theta 180;
  final period about 200 d, inclination about 0, JOI perijove above 3 R_J.
- **Solver.** Two scripts: forward (triple flyby) and backward (to the BPM and Earth), meeting at the Jupiter-approach
  state. Inner differential corrector: epoch, i, RAAN, AOP to hit Callisto and Io B.T/B.R. Outer: e and the Io B-plane
  angle. fmincon minimises JOI dV over flyby altitudes, and Earth periapsis radius over the Jupiter-hyperbola a.
- **Results (image-read):**
  - Initial guess (Tables 1-2, p.3): 3 Feb 2025 02:01:24 UTC (MJD 30709.584); a = -3,774,011.246 km, e = 1.07364,
    i = 1.950 deg, R_p = 3.887 R_J; printed V_inf 5.622 km/s.
  - Targeted state (Tables 5-6, p.7): 2 Feb 2025 21:19:19 (MJD 30709.388); a = -3,894,011.246 km, e = 1.07342,
    i = 1.956 deg, V_inf 5.704 km/s, R_p 3.999 (unit printed "[km]").
  - Achieved (Table 7): Callisto 54.8 km, theta 0.0; Io 279.0 km, theta -26.3 deg; JOI 0.264 km/s; Ganymede 97.3 km,
    theta 168.6 deg. Final orbit (Table 8): period 198.913 d, i 2.329 deg, JOI perijove 3.237 R_J.
  - BPM (Table 9): "TOF 610.0 d" (but see sec. 2, item 4), dV (7.50, -8.36, 3.50) m/s = 11.76 m/s. Earth escape (Table 10):
    from a 20,000-km circular orbit, V_inf 9.253 km/s ("requiring a NASA SLS launch").
  - Timeline (Table 12, p.9): Earth escape 25 Jun 2022 01:53:32; BPM 8 Nov 2023; Jupiter SOI 19 Nov 2024; Callisto
    6 Feb 2025 02:05:20; Io 7 Feb 07:02:38; JOI 7 Feb 11:29:08; Ganymede 8 Feb 03:54:36; first apojove 18 May 2025
    03:54:06 (Table 11: a = 9,823,607.4 km, e = 0.98388, R_p 2.215 R_J, true anomaly 180 deg).
- **Conclusion.** Total 275.8 m/s against 330 m/s for a Ganymede-Io-JOI double capture; 2.6 yr Earth escape to
  capture. The scripts can be rerun for other windows and for Callisto-Ganymede-Io orderings.

## 2. Checks (`checks_lynam.out`)

- Total dV: 264 + 11.76 = 275.76 m/s (printed 275.8). BPM magnitude sqrt(7.50^2 + 8.36^2 + 3.50^2) = 11.76.
- Table 12 MET against UTC: all seven events agree to 1 s (rounding).
- Durations: escape to first apojove 1058.08 d = 2.90 yr ("2.9 years"); escape to Ganymede 959.08 d = 2.63 yr ("2.6").
- GMAT ModJulian (epoch 5 Jan 1941 12:00 UTC): MJD 30709.584 = 3 Feb 2025 02:01 (Table 1: 02:01:24); 30709.388 =
  2 Feb 2025 21:19 (Table 5: 21:19:19); 29755.579 = 25 Jun 2022 01:54 (Table 12 escape 01:53:32); 30813.663 =
  18 May 2025 03:55 (Table 12 apojove 03:54:06). All agree to within the 3-decimal rounding of the MJD.
- Table 2 and Table 6: R_p = a(e - 1) gives 3.887 and 3.999 R_J. Agree. Table 6's unit "[km]" is a misprint for R_J.
- Table 11: R_p = a(1 - e) = 2.215 R_J; period 2 pi sqrt(a^3/mu) = 198.93 d (Table 8: 198.913). Agree.
- Table 10: V_inf = sqrt(mu_E / |a|) = 9.253 km/s; periapsis a(e - 1) = 20,099 km (a 20,000-km orbit). Agree.
- **Disagreements (re-read on the 200-dpi image):**
  1. **Table 1 V_inf 5.622 km/s** does not fit its own data. The Cartesian state gives 5.795 km/s and the Table 2 a
     gives 5.794 km/s. Table 5 (the targeted state) is self-consistent: state 5.706, a 5.704, printed 5.704.
     So 5.622 is a misprint or a different quantity; I cannot tell which.
  2. **Table 11 caption** says "Keplerian elements at first perijove"; the text (p.9) says "at first apojove", the true
     anomaly is 180 deg, and the epoch matches the first apojove in Table 12. The caption is wrong.
  3. Tables 2 and 6 differ in a by exactly 120,000.000 km (-3,774,011.246 vs -3,894,011.246 km, same decimals).
     Possible, since a is an optimiser variable, but notable.
  4. **Table 9 BPM TOF 610.0 d** (p.8, image-read) does not match Table 12. From the BPM (8 Nov 2023 11:02:34) to the
     targeted initial state (2 Feb 2025 21:19:19) is 452.4 d; to Jupiter SOI entry 376.8 d; Earth escape to BPM is
     501.4 d. None is 610 d. Probably a value from an earlier iteration of the backward script; I cannot tell.
  5. Reference slips: [23] gives the JGCD paper as "Vol. 34, 2001" (it is 2011); [24] gives Lynam & Longuski,
     "Preliminary Analysis for the Navigation of Multiple-satellite-aided Capture ...", as "Acta Astronautica Vol. 70,
     2012, pp. 33-43"; Crossref has Acta Astronaut. **79**:33-43 (2012), doi 10.1016/j.actaastro.2012.04.012.

## 3. Citation mining (30 references)

Held status checked by title and author (`ls cyclers_pdf/papers | grep -i`; `grep -i CORPUS_INDEX.md`).

| ref | work | held? | note |
|---|---|---|---|
| 21 | Lynam, Kloster & Longuski (2009), AAS 09-424 | not held | see CMDA digest |
| 22 | Lynam, Kloster & Longuski (2011), CMDA 109:59-84 | **filed in this batch** (author manuscript) | |
| 23 | Lynam & Longuski (2011), JGCD 34(5):1485-1494 | **filed in this batch** | |
| 24 | Lynam & Longuski (2012), Acta Astronaut. 79:33-43, doi 10.1016/j.actaastro.2012.04.012 | not held | = thesis Ch. 4 (navigation). Low. |
| 25 | Lynam (2014), "Broad-search algorithms ... Callisto-Ganymede-Io triple flyby sequences from 2024 to 2040, Part I", Acta Astronaut. 94:246-252, doi 10.1016/j.actaastro.2013.07.018 | not held | **new candidate, low-medium**: the broad C-G-I search over 16 yr; capture only, but it is the largest printed set of Callisto-Ganymede flyby timings in this line. |
| 26 | Lynam (2014), Part II, Acta Astronaut. 94:253-261, doi 10.1016/j.actaastro.2013.07.020 | not held | same |
| 19 | Landau, Strange & Lam (2010), AAS 10-169 | not held | low |
| 20 | Strange, Landau, Hofer et al. (2012), "Solar Electric Propulsion Gravity-Assist Tours for Jupiter Missions", AIAA 2012-4518, doi 10.2514/6.2012-4518 | not held | Jovian SEP tours; low |
| 1-8 | Galileo and Cassini operations papers (Potts & Wilson 1993; Barber et al. 1993; Wilson et al. 1997; Haw et al. 2000; Goodson et al. 2000; Roth et al. 2005; Ballard et al. 2010; Williams et al. 2009) | not held (the held Galileo and Cassini files are other works) | low |
| 9-12, 27-30 | Anderson et al. gravity-field papers; JPL SSD satellite parameters; Thomas et al. 1998 (Io shape) | not held | constants |
| 13-18 | Longman 1968; Longman & Schneider 1970; Cline 1979; Nock & Uphoff 1979; Johannesen & D'Amario 1999; Yam 2008 | not held | low |

- Related, not cited here (found by Crossref while checking the identifier): **Lynam (2015), "Broad search for
  trajectories from Earth to Callisto-Ganymede-JOI double-satellite-aided capture at Jupiter from 2020 to 2060",
  CMDA 124(1):33-50, doi 10.1007/s10569-015-9649-9** (not held). It is the one Lynam paper whose subject is the
  Callisto-Ganymede pair itself. Capture only, but worth a quick look for the gc prior-art file. **New candidate, medium.**
  Also Patrick & Lynam (2014), "Optimal SEP Trajectories from Earth to Jupiter with Triple Flyby Capture",
  doi 10.2514/6.2014-4218 (not held; low).
- Proposals: add Lynam 2015 CMDA 124 (medium, gc completeness) and Lynam 2014 Acta 94 Parts I-II (low-medium) as
  wanted rows; note the JSR 2015 journal form of this paper (doi 10.2514/1.A33159) for attribution only.

*Filed as `cyclers_pdf/papers/didion-lynam-2014-impulsive-trajectories-earth-callisto-io-ganymede-triple-flyby-capture-jupiter-aiaa-2014-4106-doi-10.2514-6.2014-4106.pdf`.*

*Wanted-list row numbers in this digest are the batch-34 numbering; the list was renumbered in batch 35.*
