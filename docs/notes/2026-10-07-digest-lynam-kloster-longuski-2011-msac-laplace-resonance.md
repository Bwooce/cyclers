# Digest: Lynam, Kloster & Longuski, "Multiple-Satellite-Aided-Capture Trajectories at Jupiter Using the Laplace Resonance" (CMDA author manuscript) (#960)

A. E. Lynam, K. W. Kloster & J. M. Longuski (Purdue), "Multiple-Satellite-Aided-Capture Trajectories at Jupiter Using
the Laplace Resonance". Journal form: Celestial Mechanics and Dynamical Astronomy 109(1):59-84, January 2011 (online
23 Oct 2010), doi 10.1007/s10569-010-9307-1 (Crossref checked: title, authors Lynam/Kloster/Longuski, volume, issue,
pages all match). The conference form is AAS 09-424 (AAS/AIAA Astrodynamics Conference, Pittsburgh, Aug. 2009),
which this text cites as "Lynam et al. (2009)" for many tables.
- File given: `5d107f14-MSAC_Laplace_Resonance.pdf`, 26 pages, LaTeX text layer (Springer `svjour` manuscript
  layout: "Celest Mech Dyn Astr manuscript No. (will be inserted by the editor)", "Received: date / Accepted: date").
  md5 b17ce60318f8c7755196c8cbb3bb61fc.
- **Manuscript vs journal form.** This is an author manuscript, not the typeset article. What I can tell:
  - Crossref lists **26 references** for the journal article; the manuscript has **25**. The journal adds Carrico &
    Fletcher (2002), "Software architecture and use of Satellite Tool Kit's Astrogator Module ...". So the manuscript
    is an earlier state than the printed text (at least one change).
  - The other 25 references match the Crossref list one for one, including the misnamed "Malcolm, M., McInnes, C.R."
    (the JGCD paper correctly gives Macdonald & McInnes, JGCD 28:365).
  - Page count is 26 in both (journal pp. 59-84), so the body is probably close.
  - The JGCD paper (p.1490) points to "Table 8 of Lynam et al. [16,17]" for the Laplacian triples. In this manuscript
    that is Table 4 (Table 8 is the GIJEC initial state). Either the journal renumbered, or the pointer means AAS 09-424.
  - I cannot fully tell beyond this. Crossref has no abstract.
- **Proposed corpus filename:** `cyclers_pdf/papers/lynam-kloster-longuski-2011-multiple-satellite-aided-capture-jupiter-laplace-resonance-cmda-109-59-doi-10.1007-s10569-010-9307-1-author-manuscript-26pp.pdf`
  (suffix as in the held `marchand-howell-wilson-2007-...-author-manuscript-33pp.pdf`).
- **How I read it.**
  - Full text layer read (pp. 1-26). Pages 8, 10, 13, 16-23 rendered at 200 dpi (pp. 10, 20, 21 at 110 dpi) and
    read on the image: Table 1 (all cells), eqs. (19)-(26), (42)-(45), Table 4, Fig. 4 and 5 captions, Tables 5-10,
    Fig. 2, 6, 7 captions.
  - Arithmetic: `checks_lynam.py`, output `checks_lynam.out` (shared with the other two digests).
  - #943 collision check: `collision.md` (written first). Verdict: **no collision** for gc-1, gc-2, ge-1, ge-2, ge-3.

## 0. Verdict

**A systematic survey of Jupiter capture using two, three or four Galilean-moon flybys plus one JOI burn. It is the
source of the Laplace-resonance phase-angle relations and of the Callisto-Ganymede-Io near-resonance numbers that
Lynam's later papers use.**
- No cyclers. Every sequence is a one-shot capture from a 5.6 km/s Jupiter arrival into a ~200-d orbit.
- **#943: no collision** with any gc or ge candidate. Capture flybys are at 10-24 km/s relative to the moon (our
  patched-conic computation); the candidates are at 1.4-8.4 km/s. One note for the gc file: the paper's 37.6-d
  near-resonance clock (16 S_Ga,Io ~ 3 S_Ca,Ga) has the same length as the gc period, but it runs on Io. It is a
  period coincidence only (`collision.md`).
- **What it gives the project:**
  - A clean, checked table of Galilean mean motions and synodic periods (Table 1; all six reproduce).
  - The Laplace phase-angle algebra (eqs. 13-41): given one pair's phase, the third moon has 2 or 3 allowed
    positions. This is a ready test for any I-E-G sequence (for example the held Lynam-Longuski 2011 triple cyclers).
  - Prior art for "repeat a moon-flyby geometry every synodic period" as a capture-window clock.
- **Catalogue implication (PROPOSAL only):** none. No periodic orbit. If the gc note gains a prior-art line, cite
  sec. 7.4 here for the 37.6-d Io-Ganymede-Callisto clock.
- **Lynam 2012 PhD (wanted row 14 in the current file; the task brief says row 15, so the list has been renumbered):**
  this paper is **Chapter 2** ("Mission design of multiple-satellite-aided captures"), nearly complete. The thesis
  List of Tables matches one to one: manuscript Tables 1-11 = thesis Tables 2.2-2.12 (same titles, e.g. 2.2 "Synodic
  Periods(S) between Galilean Moons", 2.5 "Laplace Resonance Capture Sequences", 2.9-2.11 "Integrated GIJEC ...").
  The thesis adds Table 2.1 (STK gravity fields; here only by reference to AAS 09-424) and Table 2.13 (best capture
  with zero, one or several moons; this is JGCD 2011 Table 1). It covers thesis secs. 2.2-2.10, including 2.7.1-2.7.2
  (Laplace derivation and phase angles) and 2.7.4 (Callistan near-resonance). Chapter 5 (triple cyclers) is not here.

## 1. Content (READ)

- **Model (sec. 2-4).** Circular coplanar moons, patched conics, ephemeris-free search; then STK/Astrogator 8.1.3
  integration with full gravity fields of Sun, Jupiter and the four moons, SRP and relativity. JOI is a finite burn
  (Isp 323 s, 890 N, 4317.3 kg). V_inf at Jupiter 5.6 km/s, flyby altitude 300 km, capture into a 200-d orbit unless noted.
- **Synodic periods, Table 1 (p.8), image-read:** Io-Eu 3.5255, Io-Ga 2.3503, Io-Ca 1.9789, Eu-Ga 7.0509, Eu-Ca 4.5111,
  Ga-Ca 12.5232 d; n = 203.4890, 101.3747, 50.3176, 21.5711 deg/d.
- **Double captures (sec. 6).** Io-Ganymede best (Table 2, p.9): JOI dV for GIJ/GJI/IJG/JIG at R_p = 5, 4, 3, 2, 1.01 R_J,
  e.g. JIG 330/340/333/299/228 m/s against unaided 825/735/641/524/371. Ganymede-Callisto best above 8 R_J (numbers
  only in AAS 09-424). Petal plots of approach right ascension: Ganymede-Io-JOI petals drift about 6 deg/week
  (Fig. 2); Callisto-Ganymede-JOI petals about 12.5 d apart, near-repeat every 50 d (the 3:4:7 near-resonance), next
  usable window up to two years away (Fig. 3, p.11).
- **Navigation proxy (sec. 5).** "Delta-Delta-v" (Buffington et al. 2005): dV_eq = 2 mu V_inf / (mu + V_inf^2 (r + h_p))
  (eq. 8); B-plane error = delta(dV) x time of flight (eq. 10). Table 3: a 10-km first-flyby error gives 47-370 km at the
  second flyby.
- **Laplace triple captures (sec. 7.1-7.3).** Relations eqs. 13-30 (see `collision.md` for the list). Of the 8 orderings
  of G, E, I and JOI, 4 exist and each is unique (Table 4, p.16, image-read): GEJI 213 m/s at 1.15 R_J; GIJE 254 m/s at
  2.16 R_J; EJIG 244 m/s at 2.11 R_J; IJEG 245 m/s at 1.15 R_J. GEIJ, GJIE, EIJG, JIEG are impossible (Fig. 5: Io is
  on the wrong side). They repeat every 7.0509 d.
- **Callistan triples (sec. 7.4-7.5).** Near-resonance S_Ca,Ga/S_Ga,Io = 5.3283 ~ 16/3; mismatch 0.0352 d; 24,900 km;
  window repeat about 1110 d = 3.05 yr (eqs. 42-45). Table 5 (p.19, image-read): eight C-G-I orderings, JOI dV 190-248
  m/s; best GJIC 202 m/s at 5 R_J. Integrated windows show "staircase" perijove trends (Figs. 6-7).
- **Quadruple captures (sec. 8).** Eight sequences, JOI dV 159-177 m/s (Table 7, p.22, image-read). Integrated GIJEC
  (Tables 8-10, p.23, image-read): start 21 Sep 2022 5:14, target R_p 2.5045 R_J, V_inf 5.600 km/s; flybys Ganymede
  24 Sep 17:13 (300.3 km), Io 25 Sep 3:46 (255.8 km), Europa 25 Sep 16:15 (245.7 km, B-plane angle 134.7 deg), Callisto
  26 Sep 18:06 (276.3 km); JOI 25 Sep 8:05 at 2.1439 R_J, 220.3 m/s, orbit 167.6 d. Only one quadruple found in
  Sep 2020 - Sep 2022.
- **Conclusions.** Savings over unaided capture from 210 m/s (quadruple, ~1 R_J) to 847 m/s (Callisto-Ganymede, 14 R_J).

## 2. Checks (`checks_lynam.out`)

- **Table 1:** all six S values reproduce from 360/(n1 - n2) (Ga-Ca 12.5233 vs 12.5232: last-digit rounding).
  The Eu-Ga row prints n_Ga = 50.3171 where the other rows print 50.3176; the period 7.0509 d is the same with
  either. A harmless typo.
- **Laplace rates:** 2n_Ga - 3n_Eu + n_Io = 0.0001 deg/d; 2n_Eu - n_Io = -0.7396, 2n_Ga - n_Eu = -0.7395 (printed -0.7395).
- **P_Ga = 7.1546 d and 3 S_Io,Ga = 7.0509 d** both reproduce (p.10).
- **Eqs. 42-45** reproduce: 5.3283; 0.0352 d; 24,950 km (printed 24,900); 1114 d = 3.05 yr (printed 1110 d, 3.05 yr).
- **3:4:7:** 3 P_Ca = 50.07, 4 S_Ca,Ga = 50.09, 7 P_Ga = 50.08 d.
- **Petal drift:** Ganymede falls 5.22 deg behind per 7.0509 d; with Jupiter's 0.59 deg this is 5.8 deg in the
  Sun-Jupiter frame. Agrees with "5 degrees every third synodic period" (p.10) and "about 6 degrees every week" (Fig. 2).
- **Savings:** 371 - 160 = 211 m/s (printed 210).
- **Disagreements (all re-read on the 200-dpi page image):**
  1. **Fig. 4 caption (p.17)** gives the GIJE sequence "a dV_JOI of 213 m/s". Table 4 (p.16) gives GIJE 254 m/s and
     GEJI 213 m/s. The caption's 2.16 R_J matches GIJE. So the caption dV is probably copied from the GEJI row.
  2. **Table 6 (p.22)** perijoves EJIG 2.15 and IJEG 1.11 R_J; Table 4 gives 2.11 and 1.15. Table 7 also uses 2.11 and
     1.15. Possibly Table 6 uses integrated values; the text does not say. I report it as an inconsistency.
  3. **GIJEC cost:** the text (p.23) says the integrated trajectory "was 30 m/s more costly than an equivalent
     patched-conic trajectory". Table 10 gives 220.3 m/s; Table 7's patched-conic GIJEC is 175 m/s, a 45 m/s gap.
     But the integrated orbit is 167.6 d, not 200 d, so the "equivalent" patched conic is not the Table 7 row. Not a
     clear misprint.
  4. **"7.055 days" (p.17)** against 7.0509 d everywhere else. A slip.
  5. **Fig. 6 caption (p.20)** says "a two-year gap"; the text (p.20) says "1.5 year"; the plot shows no trajectories
     from 09/2022 to 03/2024, i.e. 1.5 yr. The caption is loose.
  6. Eq. (29)'s note says "three solutions for Delta lambda_Ga,Eu for a given value of Delta lambda_Ga,Io"; the equation
     solves for Delta lambda_Eu,Io. A wording slip.

## 3. Citation mining (25 references in the manuscript; 26 in the journal)

Held status checked by title and author with `ls cyclers_pdf/papers | grep -i` and `grep -i CORPUS_INDEX.md`.

| work | held? | note |
|---|---|---|
| Lynam, Kloster & Longuski (2009), AAS 09-424, "An Assessment of Multiple Satellite-Aided Capture at Jupiter" | not held (no "09-424" or "satellite-aided" hit beyond the Lynam preview) | **new candidate, low-medium.** Holds the Ganymede-Callisto double-capture tables (Tables 4, 7) and C-G-E table (10). Capture only; useful only for completeness of the gc prior-art file. |
| Cline (1979), Celest. Mech. 19:405 | not held | satellite-aided capture background; low |
| Nock & Uphoff (1979), AAS 79-165 | not held | double-capture phasing origin; low |
| Longman (1968) RAND; Longman & Schneider (1970) JSR 7(5):570 | not held | first satellite-aided capture; low |
| Malcolm [= Macdonald] & McInnes (2005), JGCD 28:365 | not held | low |
| Yam (2008) Purdue PhD; Okutsu, Yam & Longuski (2007) AAS 07-258 | not held (the held Yam and Okutsu files are other works) | low |
| Broucke (1988), AIAA 88-4220 | not held (held Broucke files are 1968-1969 and Broucke-Prado 1993) | gravity-assist mechanics; low |
| Buffington, Strange & Ionasescu (2005), AAS 05-270 | not held (held Buffington files are Europa Clipper) | delta-Delta-v method; low |
| Heaton, Strange, Longuski & Bonfiglio (2002), JSR 39(1):17 | not held (held Heaton file is Heaton-Longuski 2003) | Europa orbiter tour; low |
| Sweetser et al. (1997) AAS 97-174; Johannesen & D'Amario (1999) AAS 99-330; Whiffen & Lam (2006) AAS 06-186; Clark et al. (2009) JEO study | not held | mission context; low |
| Wilson et al. (1997), Galileo approach | not held (held Galileo files are other works) | low |
| Laplace (1809); Sinclair (1975) Celest. Mech. 12:89; Showman & Malhotra (1997) Icarus 127:93; Musotto et al. (2002) Icarus 159:500 | not held | Laplace-resonance dynamics; Sinclair and Showman-Malhotra are the sources of eqs. 14-18. Low unless the project models the resonance libration. |
| Stastny & Geller (2008), JSR 45(2):290 | not held | navigation; low |
| Anderson et al. 1996, 1998, 2001a, 2001b (gravity fields) | not held | constants only |
| (journal only) Carrico & Fletcher (2002), Astrogator | not held (held Carrico file is Tito et al. 2013) | low |

- None of these is on the current wanted list by name. Petropoulos et al. 2000 (row 61) is cited by the JGCD paper, not this one.
- **Proposal:** add AAS 09-424 as a low-priority row ("G-C double-capture tables; capture only; completeness of the
  gc prior-art file").

*Filed as `cyclers_pdf/papers/lynam-kloster-longuski-2011-multiple-satellite-aided-capture-jupiter-laplace-resonance-cmda-109-59-doi-10.1007-s10569-010-9307-1-author-manuscript.pdf`. Check scripts, outputs and other files named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`.*

*Wanted-list row numbers in this digest are the batch-34 numbering; the list was renumbered in batch 35.*
