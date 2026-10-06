# Digest: Niehoff 1971, "Touring the Galilean Satellites" (J. Spacecraft Rockets) (#960)

J. C. Niehoff (IIT Research Institute, Chicago), "Touring the Galilean Satellites", Journal of Spacecraft and Rockets
8(10):1021-1027, October 1971. doi 10.2514/3.59764 (printed on every page margin). Presented as AIAA Paper 70-1070,
AAS/AIAA Astrodynamics Conference, Santa Barbara, 19-21 August 1970; submitted 26 October 1970; revised June 1971.
- File given: `60d87d12-niehoff1971.pdf`, 7 pages, text layer (AIAA reprint, downloaded 2015). md5 aeb93eab95d868644038f923b1ad2eeb.
  Page 7 also carries the first page of the next article (Powers & McDanell, "Switching Conditions and a Synthesis
  Technique for the Singular Saturn Guidance Problem"). That article is not part of this digest.
- **Proposed corpus filename:** `cyclers_pdf/papers/niehoff-1971-touring-galilean-satellites-jsr-8-10-1021-doi-10.2514-3.59764.pdf`
  (beside the held `niehoff-1970-touring-galilean-satellites-AIAA-paper-70-1070.pdf`).
- **Wanted list:** row 66 ("Content held in another form (AIAA 70-1070); acquire only for attribution"). This file closes it.
- **How I read it.**
  - I read all 7 pages on 110-dpi renders. Then I read eq. (1), Table 1, Table 2 and Fig. 4 (labels and T_p values) on 300-dpi crops.
  - I read all 9 sheets of the held conference paper (AIAA 70-1070, two printed pages per sheet, printed pp. 1-14) at 100 dpi.
    I read its Table 1 (physical data), Table 2 (orbits), eq. (18) and the conclusions on 300-400 dpi crops.
  - Arithmetic checks: `checks_niehoff1971.py`, output in `checks_niehoff1971.out`.
  - The #943 collision check is in `collision.md` (written first). The verdict is **no collision**.

## 0. Verdict

**This is the journal form of the held AIAA 70-1070. The method, orbits and conclusions are the same, in a shorter text.**
- The idea: a Jupiter orbit whose period is a whole number of Io periods. If it meets Io, Europa and Ganymede once, it meets
  them on every revolution, because their periods are close to 1:2:4. This is the earliest held statement of
  "repeat moon encounters by orbit-period commensurability".
- The orbits are capture ellipses (periapse 1.2-2.6 R_J, apoapse 16-61 R_J). They are not moon-to-moon cyclers.
  There are no gravity assists. Moon masses are a perturbation to cancel (about 60 m/s over 12 revolutions).
- **#943 collision: no collision.** The relative speeds at Ganymede and Europa are 6-14 and 11-18 km/s (our computation).
  The candidates have 1.4-3.9 km/s at Ganymede. Callisto is not in the resonance (see `collision.md`).
- **What it gives the project:** attribution only (the wanted row says so). It also gives a clean two-body check case,
  because Table 1 is self-consistent to 3 decimals (sec. 3).
- **Catalogue implication (PROPOSAL only):** none. If a Galilean-tour row ever cites Niehoff, cite the journal form for the
  method and the conference form for the figures the journal dropped (sec. 2).
- **The held 2026-06-17 digest of the conference paper has errors.** I list them in sec. 4 as proposed corrections.

## 1. Content (READ)

- **Satellite model (p.1021).** Circular, coplanar orbits at mean rates. Radii from JPL TR 32-1306 in units of R_J = 71,372 km.
  Mean longitudes from Sampson 1920, eq. (1), image-read:
  - L_I = 142.59987 + 203.48895 (t - 2415020.0)
  - L_E = 99.55081 + 101.37472 (t - 2415020.0)
  - L_G = 168.02628 + 50.31761 (t - 2415020.0)
  - L_C = 234.40790 + **21.51707** (t - 2415020.0). This rate is a misprint (sec. 3).
- **Laplace relation (eq. 2) and syzygy (eq. 3).** The syzygy line turns retrograde with a 437-d period relative to the
  Jupiter-Sun line (p.1021) and 486 d relative to the stars (p.1024).
- **Intercept equations (eqs. 4-16, p.1022-1023).** Intercept at a crossing true anomaly; phase each moon by Kepler's equation;
  period P = k x (Io period). The unknowns lambda_I, lambda_E, lambda_G and R_p come from three eqs. (12) plus the
  Laplace equality (15). Iterate on R_p until F = 0.
  - Eq. (16) is printed "F = lambda_E - 3 lambda_E + 2 lambda_G - pi". The first term should be lambda_I (eq. 15). This is a typo.
- **Eight intercept orders ("Modes", Fig. 2 table):** 1 A-B-C, 2 C-E-F, 3 A-C-E, 4 A-D-E; 5-8 mirror 1-4 about the apse line.
  Modes 1/5 and 2/6 need a periapse below 1 R_J, so they are not usable.
- **Table 1 (p.1023), image-read.** Periods 3.556, 7.111, 14.222, 21.333 d, headed k = 2, 4, 8, 16.

| | 3.556 | 7.111 | 14.222 | 21.333 |
|---|---|---|---|---|
| Mode 3/7 periapse, R_J | 2.556 | 2.391 | 2.290 | 2.255 |
| apoapse, R_J | 16.260 | 27.483 | 45.131 | 59.884 |
| capture impulse, km/s (v_inf 7.25) | 3.324 | 2.250 | 1.624 | 1.384 |
| radiation lifetime, revs / days | 182 / 648 | 135 / 957 | 109 / 1548 | 101 / 2149 |
| lambda_I(t0), lambda_E(t2), lambda_G(t1), deg | 110.418, 223.277, 165.218 | 107.877, 230.658, 147.132 | 107.035, 233.548, 141.618 | 106.717, 234.535, 139.863 |
| Mode 4/8 periapse, R_J | 1.368 | 1.250 | 1.213 | 1.201 |
| apoapse, R_J | 17.451 | 28.624 | 46.208 | 60.938 |
| capture impulse, km/s | 2.400 | 1.618 | 1.180 | 1.009 |
| radiation lifetime, revs / days | 4.9 / 17.5 | 3.1 / 22.2 | 2.7 / 38.6 | 2.6 / 55.3 |
| lambda_I, lambda_E, lambda_G, deg | 132.012, 211.283, 193.745 | 130.566, 215.560, 203.607 | 129.317, 217.894, 207.486 | 128.818, 218.757, 208.817 |

- **Choice (p.1023):** the 14.222-d Mode 3 orbit is "the best compromise between capture impulse, radiation lifetime, and
  potential number of satellite encounters". The 3.556-d orbits never reach Callisto. The Mode 4 lifetimes are too short.
- **Table 2 (p.1024):** two opportunities with 760-d Earth-Jupiter transfers. Arrival 1 Feb 1984 ("1981-82") and 27 Feb 1985 ("1983").
  Approach speed 7.26 and 7.12 km/s. Approach declination is set to zero. The real plane change costs up to 350 m/s
  in a bad year (1984) and under 25 m/s in a good one (1980-81).
- **Encounter histories (Fig. 4, p.1024).** These are 170-d sequences starting 20 d after arrival.
  - **The 7-day orbit gives 73 encounters** closer than the "lunar disc size" distance, one every 2.3 days.
  - The 14-day orbit gives half as many, and the 21-day orbit one third as many.
  - The periapse time was tuned to keep minimum distances outside each moon's sphere of influence where possible.
- **Opportunities (Fig. 5, p.1025):** Modes 3, 4, 8, 7 begin 20, 108, 256 and 366 d after arrival. The second Mode 3 window
  begins 79 d after the 1983-opportunity arrival.
- **Geometry (Fig. 6, p.1025):** the spacecraft approaches Ganymede from the Sun side. Closest approach is near the terminator,
  with more than 6 h for imaging, and occultations on early passes. Callisto encounters are "more or less random".
- **Control (pp.1025-1027).**
  - A 10 m/s capture error changes the period by 5.2 h and destroys the sequence (Fig. 7).
  - Two-dimensional integration with all four moon masses shows that the perturbations alone upset the sequence (Fig. 8a).
  - Fix: at each periapse, change the period by the opposite of the predicted drift, using eq. (17):
    dV ~ (R_p V_p / 3 R_a P) dP. Total "just under 60 m/s" over 12 revolutions, under 4% of the 1.624 km/s capture.
    The two largest burns, about 15 m/s each, follow the close Europa passes on revolutions 4 and 5.
  - With control, 13 encounters are within 100,000 km. Without control, only one is (Ganymede, revolution 5).
- **Conclusions (p.1027):** 14.222-d orbit, Mode 3 or 7, "typically provides 35 close encounters", "about one every 5 days";
  "the approach direction is always from the sun".

## 2. What differs from the conference paper (AIAA 70-1070, held)

Read on both sets of images.
- **The same in both:** the method and the equations. The numbering differs (300-dpi crops of conference pp.3-4):
  - conference eq. (1) is the general L_s = L_0 + n_s(t - 2415020.0), with L_0 and n_s in its Table 1. The journal writes it out as (1a)-(1d).
  - conference eqs. (10) and (11) give lambda_E(t0) and lambda_G(t0) explicitly. The journal gives (10) for Europa and states the
    Ganymede case in words. So conference eqs. (12)-(18) are journal eqs. (11)-(17).
  - Table 2 = journal Table 1, Table 3 = journal Table 2, the 73-encounter 7-day sequence,
  the 486-d syzygistic period, the control law (18) = (17), the 60 m/s total.
- **Cut from the journal:**
  - conference Table 1 (masses, radii, a, sidereal periods, eccentricities, L_0, n_s, with four sources). The journal
    folds L_0 and n_s into eq. (1) and the radii into Fig. 1 and the Fig. 4 key.
  - a paragraph on the libration of the Laplace argument (6-yr period, "certain that ... the amplitude ... is small").
  - Figs. 6-9: 14-day encounter histories for Modes 4, 8, 7 and the 1983 Mode 3 window. The conference text notes a Callisto
    pass inside the satellite's radius on orbit 9 of Fig. 9, and says this is "doubtful" because the model is only good to
    about seven Callisto diameters.
  - Figs. 10-13: flyby paths at all four moons (the journal keeps only Ganymede, as Fig. 6).
  - Figs. 14-17: the four injection-error cases (1000 km periapse radius, 1 deg apse line, 1 h periapse time, 5.2 h period).
    The journal keeps only the 5.2-h case (Fig. 7) and refers to "Ref. 4" (the conference paper) for the rest.
  - the 5-step control algorithm (integrate, compare periapse time, adjust, re-integrate, compare).
  - Figs. 19-21: controlled-encounter comparisons at Io, Europa and Ganymede. The journal keeps only Ganymede (Fig. 9).
    The conference also gives the closest controlled Callisto encounter: about 50,000 km against a predicted 36,000 km.
  - conclusions on timing: the Mode 3 sequence repeats 451 d later, then about every 450 d; Mode 7 opportunities come about
    120 d earlier.
- **Changed numbers or words:**
  - Encounter count: the conference conclusion says "an average of 36". The journal abstract and conclusion say 35.
    The conference abstract already says 35.
  - Approach direction: the conference conclusion says "always away from the sun". This contradicts its own body text.
    The journal says "always from the sun", which agrees with Fig. 6.
  - **Callisto mean motion:** conference Table 1 prints 21.57107 deg/d. Journal eq. (1d) prints 21.51707 deg/d. The journal is wrong (sec. 3).
  - **Mode 3, 3.556-d periapse:** the conference prints 2.559 and the journal prints 2.556. The conference value fits the other
    columns better (sec. 3).
  - Mode 4, 7.111-d apoapse: the journal prints 28.624. The conference scan is not clear (28.62?). Mode 3, 7.111-d lambda_E:
    the journal prints 230.658. The conference scan could be 230.558 or 230.658. I did not resolve these two conference cells.
  - Eq. (16) typo: the journal prints "F = lambda_E - 3 lambda_E + ...". Conference eq. (17) prints the correct "F = lambda_I - 3 lambda_E + 2 lambda_G - pi". So the typo is in the journal only.
  - The acknowledgement to J. Williams (JPL) and A. Friedlander (IITRI) moved from the end of the paper to the title-page footnote.

## 3. Checks (`checks_niehoff1971.py` / `.out`; Jupiter GM 1.26687e8 km^3/s^2 is our assumption)

- **Eq. (1):** L_I - 3 L_E + 2 L_G = 180.00000 deg at the epoch, and n_I - 3 n_E + 2 n_G = 1e-5 deg/d. Both agree with eqs. (2)-(3).
- **Callisto rate misprint.** 360/21.51707 = 16.731 d. 360/21.57107 = 16.689 d, which is the conference Table 1 period
  16.689018 d and the "16.7d" of journal Fig. 1. **The journal's 21.51707 has two digits transposed.** Keep the printed value
  and use 21.57107.
- **Syzygy periods.** 360/|n_I - 2 n_E| = 486.8 d ("486 days with respect to the fixed stars"). Relative to the Sun line
  (Jupiter period 11.862 yr), 437.6 d ("437 days"). Both agree.
- **Table 1 periods and radii.** For every column, a from the period agrees with (R_p + R_a)/2 to 0.003 R_J.
- **Table 1 capture impulses.** The burn from v_inf 7.25 km/s at periapse into each ellipse reproduces all eight printed
  impulses to 0.002 km/s.
- **"k = 16" is a misprint in both versions.** The base period is the average of P_I, P_E/2 and P_G/4 = 1.77779 d (from the
  eq. (1) rates). The four periods are 2.000, 4.000, 8.000 and 12.000 base periods. The last column is k = 12. Eq. (14) also lists 12 among the allowed k.
- **Mode 3, 3.556-d periapse.** With R_a = 16.260: R_p = 2.556 gives a period of 3.5550 d and a capture dV of 3.322 km/s.
  R_p = 2.559 gives 3.5558 d (target 3.5556) and 3.324 km/s (printed 3.324). **The conference's 2.559 fits. The journal's 2.556
  is probably a misprint.** This depends a little on the GM value, so I call it "probably".
- **Lifetimes.** revs x P against printed days: 647.2/648, 960.0/957, 1550.2/1548, 2154.6/2149; Mode 4: 17.4/17.5,
  22.0/22.2, 38.4/38.6, 55.5/55.3. Every difference (up to 5.6 d) is inside the rounding of the revolution counts.
- **Eq. (17)** is the exact vis-viva result, not only an approximation: with mu = R_p V_p^2 a / R_a, mu/(3 a V_p P) =
  R_p V_p/(3 R_a P). For the 14-day orbit, dP = 5.2 h needs 9.9 m/s ("10 m/sec").
- **Dates.** Fig. 4 T_p(1) = JD 2445751.852-.918, which is 20.35-20.42 d after 1984 Feb 1.0 ("20 days after arrival").
  A 760-d transfer arriving on 1 Feb 1984 leaves on 2 Jan 1982. So the Table 2 label "1981-82" fits. The text's "1980-81
  opportunity" (pp.1023-1025) does not fit the 1984 arrival. It is a label slip in both versions. A 760-d transfer arriving
  27 Feb 1985 leaves 29 Jan 1983, which matches "1983".
- **Small internal point (INFERRED).** 3.556 d is 0.497 Ganymede periods. So a 3.556-d orbit can meet Ganymede at the same
  crossing at most every second revolution, not "on every orbit" as the general principle (p.1022) states. This does not
  affect the 7-day and longer orbits that the paper uses.

## 4. Proposed corrections to the held digest `2026-06-17-digest-niehoff-1970.md` (each checked on the conference image)

- Sec. 3.1 reference longitudes: it prints 182.59587 / 188.02628 / 294.40790 for Io / Ganymede / Callisto. The image shows
  **142.59987 / 168.02628 / 234.40790** (Europa 99.55081 is right).
- Sec. 3.1 Callisto sidereal period: it prints 16.659018. The image shows **16.689018**.
- Sec. 2 "Method": it gives "Mode 3/Mode 7 (long period): k=4 (Po = 3.556 d) and k=16 (Po = 14.222 d)" and "Mode 4/Mode 8: k=8
  (7.111 d) and k=18 (21.333 d)". The table shows that **both mode pairs have all four periods**, with k = 2, 4, 8, 16 printed
  (3.556, 7.111, 14.222, 21.333 d). The last is really k = 12 (sec. 3).
- Sec. 2 "The 14-day orbit recommendation": it says the 14-day orbit gives "73 encounters ... one encounter every 2.3 days".
  **The 73 encounters belong to the 7-day orbit** (conference p.6; journal p.1024). The 14-day orbit gives 35 (abstract) or 36
  (conference conclusion), about one every 5 days.
- Header: it says 14 pages. The scan has 9 sheets carrying printed pages 1-14 (two per sheet) plus a cover.
- Sec. 4 "n_returns: 12 ... V_inf 7.26 km/s": 7.26 km/s is the Jupiter approach speed, not a moon v_inf. No moon v_inf is printed.

## 5. Citation mining (4 references, p.1027)

Checked with `ls cyclers_pdf/papers | grep -i` and `grep -i docs/notes/CORPUS_INDEX.md`.

| ref | work | status |
|---|---|---|
| 1 | Niehoff, J. C. et al., "First Generation Jupiter Orbiter Missions", Astro Sciences Rept. M-20, IITRI, "to be published" | not held (no "first generation jupiter" hit; the "M-20" hits are other files). Not on the wanted list. The source of the radiation-lifetime model. Low priority. |
| 2 | Melbourne, W. G. et al. (1968), "Constants and Related Information for Astrodynamics Calculations, 1968", JPL TR 32-1306 | not held. Constants only. |
| 3 | Sampson, R. A. (1920), "Theory of the Four Great Satellites of Jupiter", Mem. RAS 63 | not held. Ephemeris only. |
| 4 | Niehoff, J. C. (1970), AIAA Paper 70-1070 | HELD (`niehoff-1970-touring-galilean-satellites-AIAA-paper-70-1070.pdf`) |

- No new candidates. The conference version also cites Brouwer & Clemence 1961 and Price (NASA TR 1970, radii). Neither is held (no file hit for brouwer, clemence or price; the one CORPUS_INDEX "brouwer" hit is inside the Danby 1965 row). Both are data sources only.
- **Proposal for row 66:** mark it received; remove it.

*Filed as `cyclers_pdf/papers/niehoff-1971-touring-galilean-satellites-jsr-8-10-1021-doi-10.2514-3.59764.pdf`. Check scripts, outputs and other files named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`.*

*Wanted-list row numbers in this digest are the batch-30 numbering; the list was renumbered after batch 34.*
