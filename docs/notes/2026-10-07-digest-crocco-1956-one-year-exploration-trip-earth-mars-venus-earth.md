# Digest: Crocco 1956, "One-Year Exploration-Trip Earth-Mars-Venus-Earth" (Proc. VII IAC, Rome) (#960 batch 35)

G. A. Crocco, "One - Year Exploration - Trip Earth - Mars - Venus - Earth", in *Rendiconti del VII Congresso
Internazionale Astronautico / Proceedings of the VII International Astronautical Congress*, Roma, 17-22 settembre
1956, pp. 227-252. Published by the Associazione Italiana Razzi, Roma, 1956. No DOI, no report number.
- **Volume metadata from the scan.** Scan p. 30 is the volume title page. It gives the title in Italian, English,
  German and French, "Associazione Italiana Razzi" at the head, and "Roma 1956" at the foot. **No editor is named
  on the scan.** I did not supply one.
- **Page range.** The article runs pp. 227-252. Scan p. 1 is p. 227 (title page, no printed number). Scans 2-20 =
  pp. 228-246. Scans 21-23 are unnumbered fold-out plates (Figs. 8, 9 and 10). Scans 24-29 = pp. 247-252.
  Scan 30 = the volume title page.
- **Language.** The text is English (translated by Capt. Glauco Partel, p. 252). The figure captions and figure labels
  are Italian ("giorni" = days, "Orbita dell'Astronave" = orbit of the spaceship, "Ritardo" = delay).
- **File given:** `80cd3bba-Ref._2-34.pdf`, 30 pages, image-only photocopy with dark margins (HP scan, 2004).
  md5 `e28c7f80405a9a6cf599fb3c4ee92e7b`.
- **OCR copy (the file to be filed):** `crocco-ocr.pdf` in my folder, made with
  `ocrmypdf --force-ocr -l eng+ita` (Italian for the captions). md5 `236e5c712bf92e2576e54dc16c1df166`, 30 pages.
  I rendered scan pp. 1, 15 and 22 of the OCR copy at 60 dpi and compared them with the same pages of the original.
  The images survived (mean pixel difference 1.9-4.5 of 255, from re-compression only; I also viewed p. 15).
- **Proposed corpus filename:**
  `cyclers_pdf/papers/crocco-1956-one-year-exploration-trip-earth-mars-venus-earth-proc-vii-iac-rome-227-252.pdf`
- **Wanted list:** this paper is **row 10** in the current file (`2026-10-05-960-wanted-papers.md`, line 53, Tier A).
  The task brief and the TM X-53049 digest call it row 11. That is the older numbering; row 11 is now Strange & Sims 2001.
  This file closes row 10.
- **How I read it.**
  - I read all 30 scan pages on 100-dpi renders.
  - I read these on 300-dpi crops: eq. (1); eqs. (7)-(10); the p. 228 impulse figures (1760, 580, 1140);
    the p. 235 values (W = 11.7, V_E = 11.2, V_L = 16.2, 1650 s); p. 237 text heights and the Fig. 5 labels;
    Fig. 4 angles; p. 242 (407 and 320 radii); p. 244 (b/R = 2.26 and 1.59); p. 246 (Phi = 13 deg 30', 7 deg 30');
    p. 248 (MS = 206.9, a_M = 228, a = 149.5, b = 137.8, t = 125.42); p. 249 (V_M, V_A, u, 2 psi, a*, b*, t*, t_T);
    p. 250 (V_V, V_A, u, V'_A, a'_1, T_1, Phi); p. 251 (t_1, 296.4, 426, 442.9, 16.9, d = 1.4, 4, 71);
    the trial labels of Figs. 9 and 10.
  - Arithmetic checks: `checks_crocco1956.py`, output in `checks_crocco1956.out`. Sun GM and planet constants in the
    script are modern values (our assumption).

## 0. Verdict

**What it is.** A 1956 graphical and analytic study of a one-shot crewed reconnaissance flight from Earth past Mars and
Venus and back to Earth in about one year, with no stop at either planet.
- The idea: give the ship an orbit with the **same semi-major axis as Earth's orbit**, so its period is one year. Put
  its aphelion at Mars's closest point to Earth's orbit. The same ellipse crosses Venus's orbit on the way in.
  Without perturbations the ship meets Mars after 113 d, Venus after 154 d more, and Earth after 98 d more (365 d).
- Part II adds the flybys. A Mars pass changes the orbit and makes the ship late at Earth (14.7 d for a grazing pass).
  A suitably aimed Venus pass then cancels this. The result is a ballistic round trip of about 12.5-13.5 months.
- The abstract says the Brera Observatory computed "a favourable chance" in **June 1971**.

**Does the trip repeat? No.** It is a one-shot round trip. The paper never proposes to fly it again.
- The flybys move the arrival point: "Therefore, the return point on the Earth will be however different from the
  departure point B." (p. 236)
- "Only for d=∞, that is without perturbations neither from Mars nor from Venus, the contact with the departure planet
  will be allowed to happen." (p. 249, sic)
- The only "recurrence" in the paper is of the planetary alignment, not of the trajectory: "First of all, the
  existence and the recurrence of a favourable astronomical chance of the three planets, as required by the proposed
  astronautical contacts, are to be verified." (p. 252)
- The unperturbed ellipse does return to B every year (a = 1 AU), but the paper does not use that as a repeat. It says
  the one-year figure is kept "only for the title, as an attraction to read our paper; (title appeal)" (p. 236).

**R1(c) history for ev-A, ev-B, ev-C: none (no collision; not related prior art).**
- Crocco's trip is one-shot, has a Mars node, has one Venus pass used only as a timing correction, has no Venus-Venus
  return and no structure of 2-3 Earth-Venus synodic periods. ev-A/B/C are periodic Earth-Venus cyclers with
  Venus-Venus legs and periods of 1167.9 d and 1751.8 d.
- A one-year orbit is an Earth 1:1 resonance by definition. I do not count that as a link to ev-B's `RE/1:1` leg or
  ev-A's `RV/1:1` leg. It is the generic "same period as Earth" idea, as for Minovitch 1963's one-year E-V-E free
  returns (that digest also found nothing beyond the generic idea).
- How the Earth-Venus literature itself uses Crocco: Hollister & Menning 1970 (held) cite him, with Hohmann 1925, only
  as an early proposal of "interplanetary fly-by missions that would take a vehicle past both Mars and Venus before
  returning to Earth" (p. 1, ref. 2). Menning 1968 (held) says "Although Crocco suggested a trip visiting Mars before
  Venus, relaxed velocity requirements can be gained by first exploiting the stronger gravitational field of Venus".
  Neither treats Crocco as a periodic-orbit precursor.

**What it gives the project.**
- History and attribution only: the earliest held statement (1956) of a multi-planet flyby round trip in which a second
  planet's flyby is used to correct the error caused by the first. It is cited for that by Hollister & Menning 1970,
  Menning 1968, Hollister & Prussing 1965 and VanderVeen 1969.
- A small two-body check case: the ideal leg times reproduce to 0.5 d (sec. 3).

**Catalogue implication (PROPOSAL only): none.** Nothing is periodic. If the project ever keeps dated historical
`mga_tour` rows, this could be a V0 row "E-M-V-E, ~1 yr, June 1971", but the paper gives no launch date, no v_inf at
Earth departure for the executive scheme, and only graphical solutions. I do not recommend a row.

## 1. Content (READ)

### 1.1 Motivation (pp. 227-229)
- Hohmann-type Mars trip "according to Clarke": 259 d out, 455 d waiting at Mars, 259 d back, "almost three years".
  "Stuhlinger's nuclear space ship requires a total time of two years."
- Lawden's characteristic velocity for the round trip: 17.3 km/s, or a "necessary impulse" (velocity / g0) of 1760 s.
  "About 580 seconds" of this are for the Mars orbit insertion and approach. The outbound transfer alone is 1140 s,
  "aggravated by over a half" (p. 228).
- Crocco's earlier proposal: a reconnaissance flyby with no stop (Il Tempo, 4 April 1954, and a "recent paper of mine"
  with the transfer major axis along the Earth-Mars line of nodes). He says the flyby alone does not help the return:
  "three years are not sufficient to find again the coincidence with the Earth" (p. 228). Hence the one-year scheme.

### 1.2 The ideal one-year trip without perturbations (Part I, pp. 229-235, Figs. 1-4)
- Eq. (1), Kepler's third law, T in days and a in Mkm, "K ... equals about 4.98". **Printed as "T = K a^(2/3)"**
  (300-dpi crop). Only K T = a^(3/2) works with K about 5 (sec. 3). Eq. (11) later prints the right form, "5 T = a_A^(3/2)".
- For a one-year period the ship's semi-major axis equals Earth's orbit radius: a = s_0 (Earth's orbit taken circular).
- Mars is taken coplanar. A (the ship's aphelion) is the point of Mars's orbit closest to Earth's orbit, AD = d.
  So the centre C is at AC = s_0, CS = d, and e = d/a = sin(sigma) (eq. 2).
- The ship leaves Earth at B, the end of the minor axis. There r = a, so the speed equals Earth's (eq. 3), and only the
  direction changes: Earth's velocity is turned by sigma. The impulse is W = 2 V_T sin(sigma/2) (eq. 4),
  "about sqrt(1+e) - sqrt(1-e)" times V_T.
- Escape from a parking orbit: V_L^2 = W^2 + V_E^2 (eq. 5). With W = 11.7 and V_E = 11.2, V_L about 16.2 km/s,
  "a necessary impulse of 1650 seconds", against Lawden's 1760 s for the minimum-energy round trip (p. 235).
  A Venus-only one-year scheme would need 1435 s, against Lawden's 2260 s for the Venus cotangential case.
- Departure hyperbola: semi-axis a = g_0 R^2 / W^2, asymptote angle psi = arcsin(a/(r_0 + a)) (p. 232).
- **Fig. 3 ("I tre tempi del viaggio senza perturbazioni"):** B to A (Mars) 113 d, A to J_1 (Venus) 154 d, J_1 to B
  (Earth) 98 d. The text: "113 terrestrial days ... other 154 days until grazing Venus, and at last other 98 days until
  finding again the Earth, after one year, at point B" (p. 234). The ellipse meets Venus's orbit twice (J_1 and J_2);
  J_1 is the one used.
- **Fig. 4** (positions at departure, angles measured from the A direction): Mars 301 deg, Earth 292 deg, Venus 51 deg.
- Payload argument (p. 235): one-year trip 2 t equipment (3 observers) + 2 t supplies = 4 t; the three-year trip
  2 + 6 = 8 t. So the propellant for one year is "less than half".

### 1.3 Plane changes (pp. 235-239, Fig. 5)
- The real points of contact move: the return point differs from B, Mars is met near but not at its minimum distance,
  and the Venus point moves too (p. 236).
- Departure at 407 Earth radii, where Earth's pull is 1/100 of the Sun's (from Crocco's Lincei 1955 paper).
- **Fig. 5:** the ship's plane changes three times: ecliptic to plane B-S-M at departure, then M-S-V at Mars, then
  V-S-T at Venus. The text gives Mars below the ecliptic at Z = -6.40 Mkm and Venus above at +1.96 Mkm (p. 237).
  **Fig. 5 prints z = -6.23 and z = +3.04.** The node lines are drawn at 78 deg 54' (Venus) and 74 deg 21' (Mars, first
  digit not sure) from the Mars perihelion axis.
- The third (Venus) plane change "could be of a remarkable nature" unless the contact points are placed well relative
  to the lines of nodes (p. 238).
- Passing beyond "four hundred radii" would avoid perturbations but would defeat the purpose (telescope imaging; looking
  for signs of "intelligent beings"). So the perturbations are used: "it is logical to think of the second one in order
  to compensate the first one for the best" (p. 239). The Venus contact, optional in Part I, "becomes hence an integral
  part of the executive scheme" (p. 239).

### 1.4 Flyby model (Part II, pp. 239-245, Figs. 6-7)
- At Mars the planet is faster than the ship and overtakes it; at Venus the ship is faster (p. 240).
- Sphere of influence: 407 Earth radii; "almost equal for Venus, and remarkably lower for Mars, that is about 320
  Martian radii" (p. 242).
- Eqs. (6): vis-viva for Mars and the ship. Eqs. (7)-(9): hyperbola with tan psi = a/b, sin psi = a/(R+d+a),
  a/R = V_E^2/(2 u^2) (eq. 8), so **sin psi = V_E^2 / (V_E^2 + 2(1 + d/R) u^2)** (eq. 9), 2 psi = turn angle. This is
  the modern turn-angle formula.
- Eq. (10): b/R = sqrt((1 + d/R)(1 + d/R + V_E^2/u^2)), the aiming offset ("yaw of the inlet asymptote").
  For Mars, V_E = 5 km/s, u = 6.78 km/s: b/R = 2.26 for d/R = 1 and "1.59" for d/R = 0 (p. 244).
- Two Mars passes I_1 and I_2 (Fig. 6). Crocco excludes I_1 "since the increase in velocity causes an increment of the
  major axis" and a longer period (eq. 11), and chooses I_2, "where the periodic time is decreased" (pp. 244-245).
- Fig. 7 shows two Venus passes, J_1 (overtaken on the right) and J_2 (on the left). At Venus the ship flies by "on an
  asymptote almost directed towards the Sun" (p. 247). Only J_1 is computed.

### 1.5 The numerical example (pp. 246-251, Figs. 8-10)
Every value below was read on 300-dpi crops.
- **Mars pass (p. 248-249, Fig. 8).** M is 20 Mkm past Mars's perihelion, MS = 206.9 Mkm, a_Mars = 228, a_Earth = 149.5.
  The ship's ellipse: a = 149.5, b = 137.8. Time B to M: t = 125.42 d for T = 365 (planimeter on the drawing, eq. 12).
  - V_M = 26.46, V_A = 19.87, u = 6.78 km/s; V_E(Mars) = 5 km/s; eq. (9) gives 2 psi = 24 deg 42'.
  - Grazing pass (d = 0): a* = 149.8, b* = 138, t* = 268 ("in round figures"). t + t* = 393.4 d; Earth reaches the
    crossing point B'' at t_T = 378.7 d. **Delay 14.7 d.** The apse line turns by Phi = 13 deg 30'.
  - d = R: t* = 254.1, t + t* = 379.5, t_T = 374.4. **Delay 5.1 d.** Phi = 7 deg 30' (not plotted).
- **Venus pass (pp. 249-251, Fig. 9, Mars d = 0).** V_V = 35.04, V_A = 39.63, u = 11.8 km/s, "remarkably oriented into
  the direction of the Sun". V_E(Venus) = 10.4 km/s.
  - Trial 1, (2 psi)_V = 16 deg: V'_A = 36.8, a'_1 = 120.4, Phi = 41 deg 30', T_1 = 264.2 d, t_1 = 129.6 d.
    The text prints this velocity as "V_A=36.8" and says T_1 comes from "equation (12)"; it is V'_A and eq. (11).
    296.4 d (Earth-Mars-Venus) + 129.6 = 426 d; Earth reaches point 1 at 442.9 d. **Advance 16.9 d.** The 16 deg turn
    needs a pass at d = 1.4 R_V.
  - Fig. 9 trial labels (days; minus = early): 16 deg -16.8; 12 deg -7.2; 8 deg +1.3; zero turn +14.7 (= the Mars-only
    delay). X ("zero") lies just past trial 3. The text: X needs d_V = 4 R_V, and the trip lasts "about thirteen months
    and a half" (p. 251).
- **Venus pass with Mars d = R (Fig. 10).** Labels: 16 deg -30.8; 8 deg -11.5; 4 deg -2.8; zero turn +5.1. "Only the
  result is indicated": d_V = **71 R_V** and "about twelve months and a half" (p. 251).
- Conclusion (p. 252): "the research is to be considered like a laying out of the problem and not like a solution."
  Open items: the planetary alignment and its recurrence, plane-change costs, guidance instruments, and the cases set
  aside. He also says the a = 1 AU ellipse idea "should be dived deep in our search for the best".
- Acknowledgements: chief draughtsman Spartaco Migani did "all the graphical and analytical calculations"; Capt.
  Glauco Partel translated.

## 2. Comparison with held digests

- **Minovitch 1963 (JPL TR 32-464, held, digest 2026-10-07).** Minovitch computes E-V-M-E round trips with dated
  patched conics and the flyby as a free "gravity thrust" (Tables 16-17, 23). Crocco has the opposite order (E-M-V-E),
  no dates, and graphical solutions. Crocco treats the planet passes as unavoidable perturbations to be cancelled, not
  as a resource to plan with, though he notes they give "exceptional chances of free manoeuvres" after Lawden (p. 239).
  The Minovitch digest does not mention Crocco (grep), so I have no sign that Minovitch cites him.
  Both are one-shot. Minovitch's Table 14 one-year E-V-E free returns and Crocco's one-year E-M-V-E trip share only the
  generic "period of one year" idea.
- **Sohn 1964 and Gillespie & Ross 1967 (held, digests 2026-10-05), Breakwell-Gillespie-Ross 1961 (held).** These are
  E-V-M-E or E-M-V-E Mars round trips with ONE Venus swingby, chosen to cut the Mars mission's energy. Crocco is the
  earlier, cruder form of the same mission class: a Mars flyby (no landing, no stop) plus a Venus flyby on the way home.
  None is periodic; none collides with ev-A/B/C or vm2-1. The Sohn, Gillespie-Ross and Breakwell-Gillespie-Ross digests
  do not mention Crocco (grep).
- **Hollister & Menning 1970 and Menning 1968 (held).** They cite Crocco as history only (sec. 0).

## 3. Checks (`checks_crocco1956.py` / `.out`)

Agreements:
- **Leg times of the ideal scheme.** With e = d/a = 0.3814 (d from the modern Mars perihelion, 57.1 Mkm) and T =
  365.25 d: B to A 113.5 d (printed 113), A to J_1 153.8 d (154), J_1 to B 98.0 d (98). All agree.
- **Fig. 4 angles.** Angle A-S-B = acos(e) = 67.6 deg; Fig. 4 has Earth at 292 = 68 deg before A. Mars at mean motion
  moves 59.2 deg in 113 d; Fig. 4 has Mars at 301 = 59 deg before A. Venus moves 427.8 deg in 267 d and so arrives
  118.8 deg past A; J_1 is at 118.4 deg past A. All agree.
- **June 1971 (our computation, low-precision Standish 1800-2050 mean elements).** The best day in 1956-2000 for the
  Fig. 4 layout is **12 June 1971** (rms 8.9 deg: Earth -75, Mars -55, Venus +54 deg from A). The next-best year is 1969
  (rms 38.5 deg). This agrees with the abstract's "June 1971". The paper's own June 1971 date appears only in the
  abstract; the body gives no date.
- W = 2 V_T sin(sigma/2) = 11.58 km/s with e = 0.3814 (printed 11.7; with the example's e = 0.388 it is 11.79).
  Within the precision of d. V_L = sqrt(11.7^2 + 11.2^2) = 16.20 (16.2). 16.2 km/s / g0 = 1652 s (1650).
  1760 s x g0 = 17.26 km/s (17.3).
- 407 Earth radii: the 1/100 distance at 1 AU is 2.593 Mkm = 407 R_E. Agrees.
- Mars pass: V_M = 26.47 (26.46), V_A = 19.88 (19.87). With flight-path angles -3.39 deg (ship) and +0.47 deg (Mars),
  u = 6.77 (6.78). Eq. (9): 2 psi = 24 deg 41' (24 deg 42'). b/R at d/R = 1: 2.256 (2.26).
- t(B to M) by Kepler's equation on a = 149.5, b = 137.8, T = 365, with M just past the aphelion: 125.16 d
  (printed 125.42 from the planimeter). Agrees to 0.26 d. This also confirms M is past A, as Fig. 8 draws it.
- Sums: 125.42 + 268 = 393.4; 393.4 - 378.7 = 14.7; 125.42 + 254.1 = 379.5; 379.5 - 374.4 = 5.1; 296.4 + 129.6 = 426;
  442.9 - 426 = 16.9. All agree.
- Venus pass: V_V = 35.04 needs r = 108.14 Mkm (vis-viva, a_V = 108.2). V'_A = 36.8 at the r where V_A = 39.63 gives
  a'_1 = 120.40 (120.4). T_1 = 120.4^1.5 / 5 = 264.2 d (264.2): Crocco used eq. (11) with K = 5 (two-body gives 263.7).
- Venus d for 2 psi = 16 deg: d/R = 1.40 (1.4). For 8 deg: 4.18; the Fig. 9 zero by linear interpolation is at 8.6 deg,
  d/R = 3.8 (printed 4). Agree.

Disagreements (each re-read at 300 dpi; none affects the verdict):
1. **Eq. (1) is misprinted.** It prints T = K a^(2/3). With K = 4.98 that gives 140 d for Earth. K T = a^(3/2) gives 367 d.
   Eq. (11) prints the right form, 5 T = a^(3/2). K = a^1.5/T = 5.01 for Earth.
2. **b/R = 1.59 at d/R = 0 (p. 244) does not follow from eq. (10).** Eq. (10), read at 300 dpi, gives 1.24. The same
   equation gives the printed 2.26 at d/R = 1. Unresolved. (1 + V_E^2/u^2 without the square root is 1.54, so a
   dropped root is one possible cause; this is a guess.)
3. **Out-of-plane heights.** The text gives Z = -6.40 and +1.96 Mkm; Fig. 5 prints -6.23 and +3.04. Our Mars height is
   -6.40 at the perihelion and -6.19 at 5.5 deg past it, so both Mars values fit slightly different positions of M.
   Venus's maximum height is 6.40 Mkm, so both Venus values are possible. I cannot say which is right.
4. **Advance 16.9 d (text) against -16.8 (Fig. 9 label).** A 0.1-d rounding difference.
5. **p. 251 prints "d_M = 1.4 R_V, calculated by means of equation (4)".** The pass is at Venus (subscript M is a slip)
   and the formula used is eq. (9) (eq. (4) is the departure impulse). The value 1.4 agrees with eq. (9).
6. **d_V = 71 R_V for Fig. 10.** The Fig. 10 labels (-2.8 at 4 deg, +5.1 at 0 deg) put the zero near 2.6 deg by linear
   interpolation. With u = 11.8 and V_E = 10.4 that needs d/R of about 16; d/R = 71 needs 2 psi = 0.61 deg. My u comes
   from the d = 0 Mars case; the d = R case meets Venus on a different orbit, so its u is different, and I cannot check
   it. "17" (2 psi = 2.42 deg) would fit, so a digit transposition is possible. **Probable misprint, depends on u.**
7. **Mars sphere "about 320 Martian radii".** The 1/100 rule gives 347 radii at 206.9 Mkm and 382 at 227.9 Mkm. Venus
   "almost equal" to Earth: the rule gives 280 Venus radii (1.69 Mkm) against Earth's 407 radii (2.59 Mkm). These are
   loose round figures.
8. **V_A = 39.63 at Venus.** At the r that gives V_V = 35.04, vis-viva with a* = 149.8 gives 39.60. A 0.03 km/s
   difference, inside the graphical precision.
9. **The Mars-pass energy change.** The text says case I_2 is chosen because V_A decreases and the period gets shorter
   (pp. 244-245). But the worked example gives a* = 149.8 > a = 149.5 at the same r, so V'_A > V_A and the period is
   longer. Our check: u is only about 12 deg from anti-parallel to V_M, so a 24.7-deg turn cannot lower V_A much.
   Turning u one way gives V'_A = 19.93 km/s and a' = 149.87 (this matches a* = 149.8); the other way gives 21.37 and
   160.6. So the example's numbers agree with a small increase in period, not with the stated decrease. The delay of
   14.7 d comes from this and from the turn of the apse line (Phi). I record both facts. I do not state "the Mars pass
   shortens the period" as the paper's result.
10. Small: 1140 + "about 580" = 1720, not 1760 (p. 228). The 580 is "about".

## 4. Citation mining (footnotes; there is no reference list)

Checked with `ls cyclers_pdf/papers | grep -i` and `grep -i docs/notes/CORPUS_INDEX.md` for crocco, lawden, clarke,
stuhlinger, lincei and "il tempo". **No hits for any of them.** None is on the wanted list.

| where | work | status |
|---|---|---|
| p. 227 | Clarke (the 259 + 455 + 259 d Mars mission); no title given | not held. Popular source. Not needed. |
| p. 227 | Stuhlinger's nuclear (ion) space ship, two years; no title | not held. History only. |
| p. 228 | D. F. Lawden, characteristic velocity 17.3 km/s (1760 s) for the Mars round trip; no title | not held. Probably Lawden's 1950s JBIS papers. Not needed. |
| p. 228 | Il Tempo (daily paper), 4 April 1954 | Crocco's own popular piece. Not held. Not needed. |
| p. 228 | "(**) See note on page 11" | refers to a page of an earlier Crocco text, unidentified. |
| p. 236 | G. A. Crocco, "Formulazioni di Meccanica Astronautica", Proceedings Accademia Lincei, July-August 1955 | not held. The 407-radii rule and the earlier line-of-nodes Mars flyby. **New candidate, low priority** (history only). |
| p. 239 | F. Lawden, "Perturbation manoeuvres", BIS (J. Brit. Interplanet. Soc.), Nov. 1954 | not held. **New candidate, low priority.** An early gravity-assist ("free manoeuvres") paper; history for the flyby idea, not for cyclers. |
| p. 227 abstract | Brera Observatory calculation of the June 1971 opportunity | no publication cited. |

- Works that cite Crocco 1956 and are held: Hollister & Menning 1970 (ref. 2), Menning 1968, Hollister & Prussing 1965
  (ref. 1), VanderVeen 1969 (ref. list).
- **Proposal for row 10:** mark it received and remove it.

*Proposed filing: `cyclers_pdf/papers/crocco-1956-one-year-exploration-trip-earth-mars-venus-earth-proc-vii-iac-rome-227-252.pdf`
(the OCR copy, md5 236e5c712bf92e2576e54dc16c1df166), with `checks_crocco1956.py` and `checks_crocco1956.out` filed
beside it as `<pdf stem>-<file name>`.*

*Filed as `cyclers_pdf/papers/crocco-1956-one-year-exploration-trip-earth-mars-venus-earth-proc-vii-iac-rome-227-252.pdf`. Check scripts, outputs and other files named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`.*

*Wanted-list row numbers in this digest are the batch-34 numbering; the list was renumbered in batch 35.*
