# Digest: Lynam & Longuski 2011, "Interplanetary Trajectories for Multiple Satellite-Aided Capture at Jupiter" (JGCD) (#960)

A. E. Lynam & J. M. Longuski (Purdue), "Interplanetary Trajectories for Multiple Satellite-Aided Capture at Jupiter",
Journal of Guidance, Control, and Dynamics 34(5):1485-1494, September-October 2011, doi 10.2514/1.53251 (printed on
p.1485 and every page margin; Crossref checked: title, authors, volume 34, issue 5, pages 1485-1494 match).
Presented as AIAA Paper 2010-8254, AIAA/AAS Astrodynamics Specialist Conference, Toronto, 2-5 Aug. 2010; received
25 Nov. 2010; revised 21 Apr. 2011; accepted 26 Apr. 2011.
- File given: `05b40581-lynam2011_1.pdf`, 10 pages (pp. 1485-1494), publisher text layer (AIAA ARC download,
  "Downloaded by UNIVERSITY OF PITTSBURGH on January 29, 2015"). md5 0f29f48419f03cecf7eea8ee3e27494e.
- **Proposed corpus filename:** `cyclers_pdf/papers/lynam-longuski-2011-interplanetary-trajectories-multiple-satellite-aided-capture-jupiter-jgcd-34-5-1485-doi-10.2514-1.53251.pdf`
- **Wanted list:** not listed by name. (Do not confuse with the held Lynam & Longuski 2011 Acta Astronautica 69
  triple-cycler paper, doi 10.1016/j.actaastro.2011.03.011; that is a different paper.)
- **How I read it.**
  - Full text layer read. Pages 1486, 1490, 1491, 1492, 1493 rendered at 200 dpi and read on the image: Table 1,
    Tables 2-8, eq. (3), Fig. 6 date labels, Fig. 7 caption.
  - Arithmetic: `checks_lynam.py` / `checks_lynam.out` (JGCD block).
  - #943 collision check: `collision.md`. Verdict: **no collision** for any candidate.

## 0. Verdict

**The interplanetary half of Lynam's capture work: it connects the Jupiter capture sequences of the CMDA paper to
real Earth departures by backward propagation in STK, and compares total mission dV.**
- No cyclers and no repeating moon tour. **#943: no collision** (one-shot captures; moon-relative speeds 10-24 km/s).
- Its Fig. 6 is the cleanest printed evidence of the 37.6-d Callisto-Ganymede capture clock: eight
  Callisto-Ganymede-Io-JOI opportunities 25 Aug 2023 - 14 May 2024, spaced 37-38 d (mean 37.571 d = 3 S_Ca,Ga).
  These are separate captures, so this is a period coincidence with gc-1/gc-2 only (`collision.md`).
- **What it gives the project:** context only. Mission-level costs: GIJ double capture 347 m/s total with a ballistic
  Hohmann from Earth; GIJE triple 367 m/s with a 39-m/s DSM; quadruple captures estimated once per 7.4 yr.
- **Catalogue implication (PROPOSAL only):** none.
- **Lynam 2012 PhD (wanted row 14 in the current file; the brief says 15):** this paper is **Chapter 3**
  ("Interplanetary trajectory design for multiple-satellite-aided capture"), apparently whole. JGCD Tables 1-8 match
  thesis Tables 3.1-3.8 one to one by title ("dV for Best Jupiter Capture", "Integrated GIJ Flyby Sequence: Flyby
  Parameters", "... JOI maneuver", "Backward-propagated GIJE Flyby Sequences", "2:1 and 3:1 dV-EGA", "VEEGA",
  "Interplanetary Trajectories for Triple-Satellite-Aided Capture"). Section structure also matches (thesis 3.2-3.8 =
  JGCD II-VIII). Its Table 1 also appears as thesis Table 2.13.

## 1. Content (READ)

- **Table 1 (p.1486, image-read): best JOI dV, m/s, at R_p = 5, 4, 3, 2, 1.01 R_J** (C3 31.36 km2/s2 = 5.6 km/s,
  300-km flybys, 200-d orbit, R_J = 71,492 km):
  unaided 825/735/641/524/371; best single 556/526/483/416/308; best double 330/340/333/299/228;
  best triple 202/232/245/234/190; best quadruple -/-/-/175/160.
  Best triple at 5 R_J (202 m/s) is Ganymede-Io-Callisto; with 10-km flybys it could capture into a 560-d orbit with no
  deterministic dV (p.1486).
- **Method (secs. III-IV).** STK incoming-asymptote initial state (epoch, RA, Dec, velocity azimuth at perijove, true
  anomaly, R_p, C3). Inner targeter: four B-plane targets with epoch, Dec, RA, azimuth. Outer loop: C3 and R_p
  back-target the Earth B-plane ("nested backward targeting", Fig. 5). "Normalized RA" petal plots against the Hohmann
  RA give the windows.
- **Double captures (sec. IV).** Ganymede-Io-JOI: three windows per week, 120 deg apart, drifting about 6 deg/week
  (Fig. 4; Jul 2020 example). One of the four G-I orderings is available "at virtually any Jupiter arrival date";
  GIJ + GJI cover about 75 % of RA (Fig. 7). Converged GIJ (Tables 2-3, p.1490, image-read): Earth 15 Jun 2022 4:29;
  Ganymede 29 Aug 2024 19:19 (h_p 301.0 km); Io 30 Aug 2024 5:48 (331.5 km); JOI 30 Aug 2024 10:05 at 2.0439 R_J,
  347.22 m/s, 200-d orbit. No DSM.
- **Triple captures (sec. V).** Callisto-Ganymede-Io-JOI: 8 opportunities Aug 2023 - May 2024 (Fig. 6), none near the
  Hohmann RA; none for 1.5 yr either side. Forcing one costs 607 m/s total. Laplacian GIJE/EJIG have R_p ~2.1 R_J,
  GEJI/IJEG ~1.15 R_J; the resonance geometry drifts 6 deg per 7.05 d, so each aligns with a Hohmann about once every
  60 weeks. Table 4 (p.1491): four GIJE perijoves 16 Aug - 6 Sep 2024 at C3 33.0625; the 30 Aug case is retargeted
  (C3 34.4125) with a DSM of 39 m/s 0.2 yr before arrival; JOI 327 m/s; total 366.58 m/s.
- **Quadruple frequency (eq. 3, p.1491):** (4 trajectories/window) x (4 triple types / 430 d) x (120,000 km /
  2 pi x 1,882,700 km) = one per 7.4 yr.
- **Extended trajectories (sec. VI), Tables 5-8 (p.1492-1493, image-read):**

| type | launch | TOF, yr | launch C3, km2/s2 | total dV incl. JOI, m/s |
|---|---|---|---|---|
| Hohmann | (16 Jun 2022) | 2.2 | 84.0 | 367 |
| 2:1 dV-EGA | 13 Apr 2020 | 4.4 | 33.8 | 967 |
| 3:1 dV-EGA | 23 May 2019 | 5.3 | 51.4 | 567 |
| VEEGA | 25 Oct 2016 | 7.8 | 26.3 | 469 |
| Galileo (ref.) | - | 6.1 | 9.6 | 801 |

  All arrive 30 Aug 2024 10:10 on the same GIJE capture (JOI 327 m/s, DSM 2 of 40 m/s on 14 Jun 2024).
- **Conclusion.** Double capture is almost always available and usually the best total; triple or quadruple capture
  only sometimes beats it because of plane changes and DSMs.

## 2. Checks (`checks_lynam.out`)

- C3 31.36 gives 5.600 km/s; 33.0625 gives 5.750; 34.4125 gives 5.866.
- Table 1 agrees cell for cell with the CMDA manuscript: unaided = Table 2 row; best single = min(Io, Ganymede) per
  column (308 at 1.01 R_J is Ganymede); best double = JIG row; best triple = column minimum of CMDA Table 5
  (202, 232, 245, 234, 190); best quadruple = CMDA Table 7 minima at 2.16 and 1.15 R_J (175, 160).
- "Saves about 75 %": 1 - 202/825 = 75.5 %.
- Table 4 perijove spacing 7.0514, 7.0507, 7.0507 d = the 7.0509-d Laplace period.
- Table 8 totals: 39 + 327 = 366 (366.58; 367); 450 + 150 + 40 + 327 = 967; 200 + 0 + 40 + 327 = 567;
  2 + 100 + 40 + 327 = 469. All agree.
- Table 8 TOFs from the table dates to 30 Aug 2024: 2.21 (Hohmann, from 15 or 16 Jun 2022), 4.38, 5.27, 7.85 yr;
  printed 2.2, 4.4, 5.3, 7.8. Agree.
- Eq. (3): 430/16 = 26.9 d per triple; 120,000 / (2 pi x 1,882,700) = 1.014 %; gives 2649 d = 7.25 yr. Printed
  2700 d, 7.4 yr (uses 1 % and 27 d). Fine as an estimate.
- Fig. 6 dates: gaps 37, 38, 37, 38, 38, 37, 38 d; mean 37.571 d.
- **Possible copy error (re-read on the image):** Table 4 rows 3 (30 Aug, R_p 2.1198, C3 33.0625) and 6 (30 Aug,
  R_p 2.1180, C3 34.2125) print the same Earth B-plane values, B.R = -3.29e6 km and B.T = 10.5e6 km. Different C3
  should give different B-plane values. One row was probably copied.
- **Pointer:** p.1490 cites "Table 8 of Lynam et al. [16,17]" for the Laplacian perijoves. In the CMDA manuscript the
  matching table is Table 4 (see the CMDA digest).
- Ref. [23] Laplace is translated "by J. Pond" (ok). The page text calls the resonance "an exact 1:2:4 orbital resonance";
  the mean motions are exactly Laplace-locked, but the pairwise 2:1s are not exact (perijove precession -0.74 deg/d).
  Loose wording only.

## 3. Citation mining (28 references)

Held status checked by title and author (`ls cyclers_pdf/papers | grep -i`; `grep -i CORPUS_INDEX.md`).

| ref | work | held? | note |
|---|---|---|---|
| 16 | Lynam, Kloster & Longuski (2009), AAS 09-424 | not held | see CMDA digest; low-medium candidate |
| 17 | Lynam, Kloster & Longuski (2011), CMDA 109(1):59-84, doi 10.1007/s10569-010-9307-1 | **filed in this batch** (author manuscript) | |
| 15 | Kloster, Petropoulos & Longuski (2011), "Europa Orbiter Tour Design with Io Gravity Assists", Acta Astronaut. 68(7-8):931-946, doi 10.1016/j.actaastro.2010.08.041 | not held (no "kloster" file) | Jovian tour; low |
| 26 | Petropoulos, Longuski & Bonfiglio (2000), JSR 37(6):776-783, doi 10.2514/2.3650 | not held | **already wanted, row 61** |
| 27 | Sims, Longuski & Staugler (1997), JGCD 20(3):409, doi 10.2514/2.4064 | **HELD** (`sims-longuski-staugler-1997-...-doi-10.2514-2.4064.pdf`) | |
| 8 | Landau, Strange & Lam (2010), "Solar Electric Propulsion with Satellite Flyby for Jovian Capture", AAS 10-169 | not held (held Landau files are other works) | low |
| 9 | Gao (2007), JGCD 30(6):1814, doi 10.2514/1.26427 | not held | low |
| 1-7, 10-14 | Longman 1968 (RAND RM-5479-PR); Longman & Schneider 1970 (doi 10.2514/3.29992); Cline 1979 (doi 10.1007/BF01231017); Nock & Uphoff 1979; Macdonald & McInnes 2005 (doi 10.2514/1.11866); Yam 2008; Okutsu et al. 2007; Wilson et al. 1997; Sweetser et al. 1997; Johannesen & D'Amario 1999; Heaton et al. 2002 (doi 10.2514/2.3801); Whiffen & Lam 2006 | not held (checked each by title words; same-surname hits are other works) | satellite-aided capture and mission background; low |
| 18 | Carrico & Fletcher (2002), STK Astrogator | not held | low |
| 19-22 | Anderson et al. gravity-field papers | not held | constants |
| 23-25 | Laplace 1809; Showman & Malhotra 1997 (doi 10.1006/icar.1996.5669); Musotto et al. 2002 (doi 10.1006/icar.2002.6939) | not held | Laplace-resonance dynamics; low |
| 28 | Potts & Wilson (1993), AAS 93-566, Galileo VEEGA | not held | low |

- No new candidate of note. Proposal: none beyond the AAS 09-424 row proposed in the CMDA digest.

*Filed as `cyclers_pdf/papers/lynam-longuski-2011-interplanetary-trajectories-multiple-satellite-aided-capture-jupiter-jgcd-34-5-1485-doi-10.2514-1.53251.pdf`.*

*Wanted-list row numbers in this digest are the batch-34 numbering; the list was renumbered in batch 35.*
