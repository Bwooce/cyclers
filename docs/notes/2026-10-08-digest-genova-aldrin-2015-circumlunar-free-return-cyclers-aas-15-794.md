# Digest: Genova & Aldrin 2015, "Circumlunar Free-Return Cycler Orbits for a Manned Earth-Moon Space Station" (AAS 15-794) (#960, #948)

A. L. Genova (Mission Design Division, NASA Ames Research Center) and B. Aldrin (Buzz Aldrin Enterprises), "Circumlunar
Free-Return Cycler Orbits for a Manned Earth-Moon Space Station", **AAS 15-794** (printed top right of p.1),
AAS/AIAA Astrodynamics Specialist Conference, Vail CO, 9-13 Aug 2015. NTRS 20160004674. The NTRS record lists report
numbers ARC-E-DAA-TN31042 and AAS-15-794 and the publication date 2015-08-09. No DOI: Crossref has no record of it, and
the AAS proceedings volume is not DOI-registered.
- File: upload `ea5d7dce-20160004674_1.pdf`, 20 pp., md5 ee35eda13c0aa2dcc105b4904faba63a. It is a Word 2007 digital
  PDF with a publisher text layer, so no OCR is needed.
- Proposed filename: `cyclers_pdf/papers/genova-aldrin-2015-circumlunar-free-return-cycler-orbits-manned-earth-moon-space-station-aas-15-794-ntrs-20160004674.pdf`.
- How I read it: the full text layer, then pp. 1-5, 8, 10, 11 and 13-18 on 100-110-dpi images, with 250-400-dpi zooms
  of Fig. 2 (centre), Fig. 14 and Table 1. Every number quoted from AAS 15-794 is on one of those pages. Pp. 6-7, 9 and 12
  hold only figures and captions, and pp. 19-20 are the reference list; for those I used the digital text layer.
  Held preprint: pp. 2, 3, 4, 6 and 9 read on 100-dpi images (the 50x50 fields, the 13/14/8 m/s burns, "20 to 62 m/s per
  cycle", the 19 m/s per 2 months and the 947 m/s total); other held values are from its text layer. Table 1 was read on the image and matched to the text layer cell by cell (`witness-comparison.tsv`, 30/30
  and the footnote). The arithmetic checks are in `scripts/check_genova_arith.py` and `.out`.
- Tables transcribed: `genova-aldrin-2015-aas-15-794-tables.yaml` (Table 1 and the per-cycler text values).
- **This is not the held Genova-Aldrin paper.** The held file is `genova-aldrin-2015-earth-moon-cycler-AAS-15.pdf`:
  "A Free-Return Earth-Moon Cycler Orbit for an Interplanetary Cruise Ship", AAS 15-___, 11 pp., NTRS 20150018049,
  md5 3efc6cb337fb442a90bf15e52730fdcc. The relation between the two is in sec. 2.

## 0. Verdict

**A full-ephemeris design paper with deterministic maintenance burns. It gives no initial state and no CR3BP
quantities. It cannot be integrated, and none of its cyclers can be the "same orbit" as a CR3BP catalogue row.**
- **Model (pp.2-3).** STK/Astrogator with Earth and Moon 30x30 gravity, Sun 4x0, Jacchia-Roberts drag, solar and thermal
  radiation pressure, and "the best-known ephemeris data", propagated with RK 8/9 over 2018-2022. It is neither CR3BP
  nor bicircular.
- **Data.** No state vector, Jacobi constant, epoch state or mass ratio is printed for any cycler. mu = 0.01215 appears
  only in the history on p.1. Table 1 and the text give periods, perigee gaps, target altitudes and ΔV ranges.
- **Five cyclers (Table 1).** All are flown with impulsive maintenance burns; Arenstorf's two orbits are flown as
  comparisons.

| Cycler | Resonance / shape | Cycle | Perigee passes per cycle | Perigee alt. | Perilune alt. | Maintenance ΔV |
|---|---|---|---|---|---|---|
| Shamrock (3-leaf clover; Aldrin's C-2-R) | 3:1, three lobes in the rotating frame | 26.3 d (p.5) | 3 (two 9.5-d holding orbits plus a ~7-d figure-8 leg) | 3,000 km target (p.4, p.17) | 3,000 km, lunar far side (p.4) | 19-60 m/s per cycle (p.5); Table 1 ~20/~62 m/s per sidereal month |
| Mushroom (Aldrin's C-1-R) | 2:1 (INFERRED: p.3 lists "2:1 and 3:1 ... C-1-R and C-2-R" in that order, and the shamrock is C-2-R and 3:1) | 25.7 d (p.5) | 2 (one phasing orbit, max gap 18.5 d, plus the figure-8 leg) | 5,000 km "for the other cyclers" (p.18) | 3,000 km implied (shared figure-8 segment) | 29-70 m/s per cycle (p.5); Table 1 ~32/~74 |
| Hybrid (mushroom and shamrock) | switches between the two | ~1 sidereal month | - | - | - | max 45 (p.5) or 47 (p.17, Table 1); average 31.5 m/s per month (p.17); Table 1 min ~18 |
| 4-leaf clover (Arenstorf 1963, ref. 18) | far-side figure-8 and front-side reverse figure-8 in turn | 55 d (p.8); Table 1 ~2 months | 4 per 55 d (gaps 13 and 22 d; ours from p.8) | 5,000 km (p.18) | not stated | "discovered without deterministic ΔV"; up to 55 m/s per month near lunar perigee (p.8); Table 1 ~4/~55 |
| Hybrid (3- and 4-leaf) | switches | ~1 or ~2 months | - | - | - | max ~37, average ~18 m/s per month (p.8); Table 1 min ~4 |
| (comparison) Arenstorf monthly "Big Loop" (refs. 15-17) | "2:1 resonance cycler" (p.17); Fig. 14 | "monthly" | - | "≈40,000 km" (p.10) | not stated | not given |

- **The Arenstorf 3:1 reference (`#948` R4).**
  - The source is ref. 18: Arenstorf, "Periodic Trajectories Passing near both Masses of the Restricted Three-Body
    Problem", Proc. XIV IAC, Paris 1963, Vol. IV, p.85. The held preprint cites the same work as its ref. 14.
  - It is wanted-list row 23, still not held.
  - Claim (p.3): "Arenstorf presented a 3:1 resonance orbit but without the required close Earth passes (Fig. 2,
    center)". The Fig. 2 caption reads "Arenstorf shows no close-Earth approaches are possible for a 3:1 resonance cycler
    orbit in the restricted three body problem".
  - Fig. 2 centre is a photograph of a printed figure: a three-lobed rotating-frame orbit with a small loop around M and
    no pass near E. The figure carries no numbers.
  - The same paper (ref. 18) is cited as the source of the 4-leaf clover (p.8, Fig. 10).
  - The paper reports no search. The "not shown to exist" sentence (p.3) is an assertion.
  - The sketch does not fit an interior 3:1 orbit. It loops round M and, measured on the image, reaches about 1.41 L
    from E (`em-literal-check.md` sec. 4). An interior 3:1 Kepler ellipse (a = 0.48, apogee <= 0.96 L at mu = 0) can do
    neither, so the a > 1/2 bound in the held Arenstorf AIAA J digest does not apply to this drawing. What Arenstorf
    meant by "3:1" here stays open until the Proc. XIV IAC paper (wanted-list row 23) is in hand. R4's framing is
    unchanged.
- **What it gives the project:**
  - The final, full version of the conference paper that the catalogue's `genova-aldrin-2015-em-3petal-cycler` row is
    built on (sec. 2). Two catalogue items rest on the preprint only and should be revisited (sec. 4).
  - A named, figure-level identification of the Arenstorf "monthly"/"Big Loop" 2:1 cycler with the `#997` F3 family.
    The Fig. 14 top-right drawing is Arenstorf 1963 AIAA J Fig. 1 ("Closed path of P in rotating co-system with m = 1,
    k = 2"; I compared the held AIAA J p.239 image: same shape and labels). It matches F3 members at C = 2.59 to
    0.008 L RMS, with the same sense of motion as the drawn arrows. See `em-literal-check.md`.
  - No positive control for a corrector, since no states are printed.

## 1. Content (by section)

- **Introduction (pp.1-3).**
  - History of Earth-Moon periodic orbits: Darwin, Moulton, Stromgren, Egorov 1958, Message, Newton 1959, Thuring and
    others.
  - The "figure-8" (Apollo) free return is traced to Egorov's 1953-55 dissertation (published 1958) and to Lisovskaya
    1957 (6,100 km perilune, Fig. 1).
  - Chebotarev 1957 computed a circumlunar orbit with a 29,860 km perilune, which is "not eight-shaped" (STK
    re-creation).
  - Lieske (RAND P-1293, 1958) drew a figure-8 that "purposely avoided an Earth reentry", relevant to the cyclers here.
  - Requirements for the station (p.2): the cycler lies in the lunar orbit plane, passes "near (≈5,000 km altitude)
    both the Earth and Moon", and meets the Moon "at least once per month".
- **Shamrock and mushroom (pp.3-7).**
  - Aldrin 1985 (ref. 27, SAIC presentation) proposed 2:1 and 3:1 free-return cyclers, C-1-R and C-2-R. Uphoff lists
    the same resonances (ref. 28, Uphoff & Crouch 1993).
  - Shamrock maintenance, three burns per cycle (p.4):
    - an anti-velocity burn at perigee (Fig. 3, C), which lowers the energy for lunar phasing;
    - a normal burn at the apogee before the figure-8, which keeps the plane;
    - a velocity-direction burn at perigee (Fig. 3, E) after two 9.5-d holding orbits, which puts the perilune at
      3,000 km on the far side.
    - No burn magnitudes are given.
  - The mushroom reverses the in-plane burns (p.4).
  - "The Moon's varying Earth-range (due to the eccentricity of the lunar orbit) was observed as the primary cause of
    variance in the ΔV requirements" (p.4).
  - The analysis period is 553 d, "the time needed for the shamrock cycler to repeat itself in an inertial frame" (p.5).
  - The shamrock has its ΔV minimum near lunar perigee and the mushroom near lunar apogee. Switching between them gives
    the hybrid.
- **4-leaf clover (pp.8-10).** "In 1963, Arenstorf presented a cycler orbit that resembles a four-leaf clover" (p.8).
  - It has a reverse figure-8 segment (front-side lunar pass, counterclockwise viewed from north). Egorov 1958 had the
    same segment in a periodic circumlunar orbit (Fig. 11).
  - Earth gaps are 22 d and 13 d. Lunar encounters are every ≈27.5 d. It repeats every 55 d.
  - A hybrid with the shamrock lowers the maximum to ≈37 m/s per month.
- **Taxi rendezvous (pp.10-12).**
  - C3 is "≈ -2.3 km2/s2" for Arenstorf's high-perigee cycler and "≈ -1.75 km2/s2" for a low-perigee one (p.10).
  - Coplanar taxi rendezvous ΔV, 6-48 h after launch, is plotted (Fig. 16): Arenstorf Big Loop about 1,600 m/s at
    6 h down to about 380 m/s at 48 h; 3-leaf clover about 100 m/s or less. These are figure readings.
  - Fig. 17 shows lunar-plane inclination between ≈18 and ≈28.5 deg over 18.6 yr. Both cyclers were "solved with a
    figure-8 lunar encounter occurring in March of 2034".
- **Lunar stations (pp.13-15).**
  - Equatorial LLO from the cycler: 326 + 88 + 492 m/s, plus 33 m/s, for a total of 939 m/s.
  - Equatorial LLO to the cycler: 51 + 281 + 533 + 46 m/s, printed total 860 m/s. The parts sum to 911; see sec. 3.
  - Polar LLO to the cycler: 801 + 100 = 901 m/s.
  - WSB transfer from polar to equatorial: 573, 87 and 342 m/s, printed total 2,120 m/s. The remainder is not itemised.
- **EM-L2 halo (pp.15-16).** 11 + 184 + 54 = 249 m/s, 17 d (2018 Jun 17 to Jul 4). The halo state came from D. W. Dunham:
  "northern Class I quasi-periodic halo", Z amplitude ≈7,000 km, Y amplitude ≈33,000 km.
- **Mars departure (pp.16-17).** Injection of 736 m/s at perigee (against ≈4,100 m/s from LEO). There is a 265 m/s plane
  change at a ≈665,000 km apogee, a 28-d transfer and a 2020 Aug 3 injection. Switching to "Arenstorf's 2:1 resonance
  cycler (Fig. 14)" uses its "relatively high argument of perigee rotation rate".
- **Conclusions and Table 1 (pp.17-18).**
  - "The perigee altitudes directly following the figure-8 segments were targeted to 3,000 km for the shamrock cycler
    orbit (and 5,000 km for the other cyclers) to avoid perigee maintenance maneuvers."
  - "the orbit types presented are quite sensitive to errors in velocity [8, 48] but this sensitivity has not been
    quantified in this paper."

## 2. Relation to the held preprint (AAS 15-___, NTRS 20150018049)

Both NTRS records give the same conference (Vail, 9-13 Aug 2015), the same publication date (2015-08-09) and the same
authors. The held 11-pp file ends its introduction with "More destinations, such as distant retrograde orbits (DROs),
near-Earth asteroids (NEAs), and Mars will be explored in the full paper, if accepted" (held p.2). AAS 15-794 contains
the Mars, LLO and WSB material. **So AAS 15-794 is the full conference paper, and the held file is the earlier
submission (INFERRED from these facts; neither file says so).** They are not content-identical. The differences:

| Item | Held preprint (AAS 15-___) | AAS 15-794 |
|---|---|---|
| Title | "A Free-Return Earth-Moon Cycler Orbit for an Interplanetary Cruise Ship" | "Circumlunar Free-Return Cycler Orbits for a Manned Earth-Moon Space Station" |
| Gravity model | Earth 50x50, Moon 50x50, Sun 4x0 (held p.2) | Earth 30x30, Moon 30x30, Sun 4x0, plus Jacchia-Roberts, SRP, TRP (p.2) |
| Why the 3:1 cycler exists | "the addition of solar gravity to the astrodynamics model and a modest ∆V maneuver" (held p.3) | "the addition of modest ∆V maneuvers" (p.3); solar gravity is not named |
| Cause of the ΔV variance | "Solar gravity perturbations cause variance in the ∆V requirements" (held p.4) | "The Moon's varying Earth-range (due to the eccentricity of the lunar orbit) was observed as the primary cause" (p.4) |
| 3:1 maintenance burns | 13 m/s at perigee C, 14 m/s at perigee E, 8 m/s out-of-plane 0.5 d later (held p.3) | perigee C anti-velocity; apogee normal burn before the figure-8; perigee E velocity burn; no magnitudes (p.4) |
| 3:1 ΔV range | "from 20 to 62 m/s per cycle ... or 26 days" (held p.4); abstract "average of 39 m/s per month" | "19 to 60 m/s per cycle (26.3 days)" (p.5); Table 1 ~20/~62 m/s per lunar sidereal month |
| Holding-orbit Earth gap | "every 7 or 10 days" (abstract); "9.5 days" (held p.3); "just under 10 days" (Fig. 4) | "≈7 or ≈9.5 days" (p.3); 9.7 d (Table 1) |
| Name of the 3:1 cycler | "3-petal" | "shamrock" / "3-leaf clover" |
| 4-leaf / 4-petal | "∆V cost of 19 m/s per 2 months"; 3-to-4 petal transition 68 m/s, 58 d (held p.6) | "as much as 55 m/s per lunar sidereal month" near lunar perigee (p.8); Table 1 min ~4, max ~55 |
| Other cyclers | Arenstorf's 5-petal (ref. 13) and 4-petal (ref. 14); DLS (Farquhar-Dunham); Uphoff | adds the mushroom, two hybrids and the Arenstorf "Big Loop" (refs. 15-17); drops the 5-petal |
| Lunar station rendezvous | 325 + 88 + 492 m/s, plus 32 m/s; printed total 947 m/s (held p.9). The held Fig. 12 caption says "82 m/s needed above to change inclination about 7 degrees", against 88 m/s in its text | 326 + 88 + 492 m/s, plus 33 m/s; printed total 939 m/s (p.13) |
| Inclination | "the 24 degrees contained in the cycler"; "about 174 degrees in the Moon TOD frame" (held p.8) | lunar-plane inclination ≈18-28.5 deg over 18.6 yr (p.10, Fig. 17); "174 degrees" for the 12-h lunar orbit (p.13) |
| Extra sections | crew launch rendezvous (84 + 43 m/s) | LLO-to-cycler, polar LLO, WSB, Mars departure, taxi rendezvous vs inclination |

The ΔV units do not agree. The preprint's "20 to 62 m/s per cycle" are the same two numbers that AAS 15-794 prints per
lunar sidereal month (Table 1). AAS 15-794's per-cycle range 19-60 m/s, scaled by 27.32/26.3, gives 19.7-62.3 m/s per
month (ours). So the preprint's "per cycle" label is probably a per-month figure, or comes from a different run.
`identity-genova-aldrin-2015.md` states the identity verdict.

## 3. Checks (ours; `scripts/check_genova_arith.out`)

- ΔV sums:
  - 326 + 88 + 492 + 33 = 939, as printed.
  - 801 + 100 = 901, as printed.
  - 11 + 184 + 54 = 249, as printed.
  - **LLO to cycler: 51 + 281 + 533 + 46 = 911, but the paper prints 860 (p.14).** 860 = 281 + 533 + 46, so the 51 m/s
    plane change is left out of the total.
  - The held preprint's 325 + 88 + 492 + 32 = 937 is printed as 947.
  - The WSB total of 2,120 m/s holds 1,002 m/s of itemised burns. The other 1,118 m/s (presumably the circular-orbit
    changes) is not itemised, so this is not an error.
- Per-cycle to per-month conversion (x 27.32/26.3 and x 27.32/25.7):
  - Shamrock: 19.7-62.3 m/s, against ~20/~62 in Table 1. Consistent.
  - Mushroom: 30.8-74.4 m/s, against ~32/~74. The minimum is 1.2 m/s higher than the conversion gives.
- Cadences:
  - Shamrock: 2 x 9.5 + 7 = 26.0 d; with 9.7, 26.4 d; the paper says 26.3 d. 553/21 = 26.33 d.
  - Mushroom: 25.7 - 18.5 = 7.2 d for the figure-8 leg.
  - 4-leaf: 2 x 27.5 = 55 d. The decomposition 7 + 13 + 22 + 13 = 55 d is ours, built from the p.8 gaps.
- The hybrid maximum is "45 m/s" on p.5 and "47 m/s" on p.17 and in Table 1.

## 4. Catalogue implications (PROPOSALS only; catalogue not edited)

1. **`genova-aldrin-2015-em-3petal-cycler`, `model_assumption: bicircular`.** The 2026-06-10 `#184` reclassification
   rests on the preprint's "addition of solar gravity" sentence.
   - The full paper drops that sentence and gives a full-ephemeris model (30x30 fields, drag, SRP, TRP). The preprint
     also used a full ephemeris (50x50 fields).
   - Propose `analytic-ephemeris`, following the Wittal precedent (`#211`/`#216`).
   - Change the notes so they no longer say the orbit "exists only when solar gravity is added". AAS 15-794 credits
     the burns, and calls the lunar eccentricity the main driver of the ΔV variance.
2. **The same row, `corroborating_sources`.** It now reads `authors: ["Genova, A. L."]`, year 2016, venue
   "NTRS 20160004674", note "Companion / follow-up report ... covering the same orbit family". Propose:
   - authors Genova, A. L. and Aldrin, B.; year 2015; AAS 15-794, AAS/AIAA Astrodynamics Specialist Conference, Vail CO,
     Aug 2015; NTRS 20160004674; ARC-E-DAA-TN31042.
   - note: "full conference paper; the held NTRS 20150018049 file is the earlier submission".
   - Consider making AAS 15-794 the `first_published` record and the preprint the corroborating one. The row's
     `first_published.note` says the AAS number "is not exposed"; it is 15-794.
3. **The same row, ΔV fields.** "20 to 62 m/s per cycle" (preprint) conflicts with AAS 15-794 Table 1 (~20/~62 per
   lunar sidereal month) and its text (19-60 per 26.3-d cycle). The per-cycle burns 13/14/8 m/s are preprint-only.
   Propose citing AAS 15-794 for the range and flagging the unit.
4. **`arenstorf-em-figure8-1963`.** The row is vague: no state, a "figure-8" family label, and a period "on the order of
   one Earth-Moon synodic period". AAS 15-794 distinguishes three Arenstorf orbits:
   - the monthly 2:1 "Big Loop" (refs. 15-17, ≈40,000 km perigee);
   - the 4-leaf clover (ref. 18);
   - the 3:1 sketch (ref. 18).
   The Big Loop sketch matches `#997` F3 (see `em-literal-check.md`). PROPOSAL for the lead or `#972`: tie the row to
   F3 at family level, or split it. No new row.
5. **No new rows from this paper.** Every cycler here needs maintenance burns in a full ephemeris, and no state is
   printed.

## 5. Citation mining (48 references)

Held or not held, checked with `ls cyclers_pdf/papers | grep -i` and CORPUS_INDEX; wanted-list rows where listed.
Only the references that bear on cyclers or Earth-Moon periodic orbits are listed. Refs. 1-7, 21-26, 29-33, 36-46 are
history, station or habitat context, or lunar science, and none is a cycler source.

| Ref | Item | Status |
|---|---|---|
| 8 | Egorov, V. A. (1958), "Certain Problems on Moon Flight Dynamics", Russian Literature of Satellites Part I, pp.115-175 | not held; figure-8 and reverse-figure-8 periodic circumlunar orbits (Fig. 11); not on the wanted list. Candidate, low priority (historical, no states expected) |
| 9 | Message 1958, AJ 63:443 | not held (asymmetric periodic orbits; not Earth-Moon cyclers) |
| 10 | Newton 1959, Smithson. Contrib. Astrophys. 3:69 | HELD (`newton-1959-...`); the `#997` F2/F3 seed source |
| 11 | Thuring 1959, Astronautica Acta 5:241 ("Mondeinfangbahnen") | not held; lunar capture orbits; low priority |
| 12 | Huang 1962, AJ 67:304 | not held; Moon-probe orbits; low priority |
| 13 | Ehricke 1962, Space Flight Vol. 2 | not held; textbook |
| 14 | Broucke, "Recherches d'orbites periodiques ... (systeme Terre-Lune)", Louvain | not held as such; Broucke 1968 JPL TR 32-1168 (Earth-Moon periodic orbits) IS held |
| 15 | Arenstorf 1962, ARS Paper 2604-62 | not held; conference form of ref. 16 |
| 16 | Arenstorf 1963, AIAA J 1(1):238 | HELD (`arenstorf-1963-existence-...aiaa-j-1-238...`) |
| 17 | Arenstorf 1963, "Periodic Solutions ... Analytic Continuations of Keplerian Elliptic Motions", NASA TN (May 1963) | NASA TN D-1859 form not held; the Amer. J. Math. 85:27 version IS held (`arenstorf-1963-periodic-solutions-...amer-j-math-85-27...`) |
| 18 | Arenstorf 1963, "Periodic Trajectories Passing near both Masses ...", Proc. XIV IAC Paris, Vol. IV p.85 | **not held; wanted-list row 23**. Holds the 3:1 sketch (Fig. 2 centre) and the 4-leaf clover (Fig. 10) |
| 19 | Davidson 1964, "Numerical Examples of Transition Orbits in the Restricted Three-Body Problem", Astronaut. Acta 10:308 | not held; not on the wanted list. The held preprint credits Arenstorf "with the help of Davidson" for the 4-petal. Candidate, medium priority for `#948`/`#972`: likely numerical Arenstorf-type orbits |
| 20 | Deprit & Henrard 1965, AJ 70:271 (symmetric double-asymptotic orbits) | not held; low priority |
| 27 | Aldrin 1985, "Cyclic Trajectory Concepts", SAIC presentation, JPL, 28 Oct 1985 | not held; C-1-R and C-2-R sketches (Figs. 2, 5); grey literature |
| 28 | Uphoff & Crouch 1993, "Lunar Cycler Orbits with Alternating Semi-Monthly Transfer Windows", JAS 41(2):189-205 (AAS 91-105) | not held (`ls | grep -i uphoff`: only the Uphoff-Roberts-Friedman 1976 tour paper, on the wanted list as row 40); no Crossref DOI found. **A lunar cycler paper. Candidate, medium-high priority**: inclined "back-flip" lunar cyclers, not in the catalogue |
| 34 | Casoliva et al. 2010, JGCD 33(5) | HELD |
| 35 | Carrico et al. 2011, AAS 11-454 (IBEX 3:1 lunar-resonant orbit) | not held; real-mission 3:1 resonant orbit kept away from the Moon; low priority |
| 47 | Farquhar & Dunham 1981, J. Guid. Control 4(2):192 (double lunar swingby) | not held; DLS background; low priority |
| 48 | Lieske 1957, "Accuracy Requirements for Trajectories in the Earth-Moon System", Convair-OSR symposium | not held; with ref. 32 (Lieske 1958, RAND P-1293, figure-8 avoiding re-entry) |

Also found: Thangavelu 2022, "BuzzCraft: Evolution of A Cislunar Cycler Architecture for Permanent Lunar Settlement
Logistics", ASCEND 2022, doi 10.2514/6.2022-4345 (Crossref). It is not held, not on the wanted list and not cited by this
paper. It is a later cislunar-cycler architecture paper in the Aldrin line. Candidate, low-medium priority.

*Filed as `cyclers_pdf/papers/genova-aldrin-2015-circumlunar-free-return-cycler-orbits-manned-earth-moon-space-station-aas-15-794-ntrs-20160004674.pdf`. Check scripts, outputs and notes named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`. Table transcription: `data/sources/genova-aldrin-2015-aas-15-794-tables.yaml`.*
