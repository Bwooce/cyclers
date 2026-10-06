# Digest: Hollister 1963, "The Mission for a Manned Expedition to Mars", MIT Sc.D. thesis (#960)

W. M. Hollister, "The Mission for a Manned Expedition to Mars", Sc.D. thesis (Instrumentation), Department of
Aeronautics and Astronautics, MIT. Submitted 10 May 1963. Supervisors: W. Wrigley, R. H. Miller, E. J. Frey.
MIT DSpace hdl 1721.1/15698. The handle resolves (hdl.handle.net API, 2026-10-07) to
dspace.mit.edu/handle/1721.1/15698. The DSpace page is captcha-gated, so I could not check the record's title
there. No DOI (thesis). Computer work: IBM 7090, MIT Computation Center, Problem M 2600. Funding: NASA NsG 254-62
(DSR 9406).
- File given: `10973410-MIT.pdf` (scratch), 158 pages, image-only typewritten scan, md5 da8075f7d43ccee2f039b992429b0e32
  (not modified).
- **OCR copy (the file to file):** `hollister1963-ocr.pdf` in this folder, md5 **dfea26d59e86c017b9326caab55475b8**.
  I made it with `ocrmypdf --force-ocr -l eng` (4 min 21 s). I rendered pages 1, 60 and 150 of both files at 50 dpi:
  the sizes are the same and the mean pixel difference is 3-6 grey levels in 255. This is the image optimisation.
  I looked at page 60 of the copy, and it is clean.
- **Proposed corpus filename:**
  `cyclers_pdf/papers/hollister-1963-mission-manned-expedition-mars-mit-scd-thesis-dspace-1721.1-15698.pdf`
- **Wanted list:** row 42 (the numbering in the file as read on 2026-10-07, after batch 34). This file closes it.
- **Page convention:** the thesis's printed page, with the PDF page in brackets. PDF = printed + 8 for the body.
- **How I read it.**
  - I read the whole OCR text (it is fair, but it often turns 4 into 1, for example "1.00 days" for 400).
  - On page images (90-200 dpi; 200-dpi crops for tables):
    - the abstract;
    - pp.6-7 [14-15] (literature and priority claim);
    - p.12 [20] (synodic periods);
    - pp.50-51 [58-59] (flyby equations, Vc values);
    - p.65 [73] (Table 9.1);
    - pp.66-67 [74-75] (eq. 9.1, braking limits);
    - pp.75-82 [83-90] (all of Ch. 11 and Table 11.1);
    - pp.83, 85, 86, 88, 89 [91, 93, 94, 96, 97] (Ch. 12, Tables 12.1 and 12.2);
    - p.96 [104];
    - pp.110-112 [118-120] (App. D and Table D.1);
    - Figs. 2, 17 and 18 [129, 145, 146];
    - references pp.145-149 [153-157].
  - Chapters 3-8 and 10, 13 and Apps. A-C and E I read in the OCR text only. I quote no numbers from them, except the
    flyby equations, which I read on the image.
  - Arithmetic checks: `checks.py`, output `checks.out`.
  - I wrote the #942 check first: `evm-repeat.md`. The verdict there is **no repeating orbit and no collision**.

## 0. Verdict

**A one-shot Mars mission study, and the earliest Hollister work on Venus flybys. It has no periodic orbit and
no cycler idea.**
- **What it is:** a 1963 mission-design thesis. Its main proposal is to make one leg of a 400-day Mars round trip
  "bi-elliptical". In the best case, the mid-course burn of that leg is moved to a close Venus pass. The
  recommended programme is: launch 1 Sep 1970, powered Venus flyby 31 Dec 1970, Mars 19 May 1971, and a return
  via a second powered Venus flyby (3 Apr 1972) to Earth on 22 Jul 1972.
- **Priority claim (p.7 [15]):** "the author has found no mention in the literature of ... trips to Mars via
  bi-elliptical transfer or via a Venus encounter that includes a significant velocity change near Venus." This
  predates Sohn 1964 (held), whose Venus swingby is unpowered.
- **#942 R1(c):** no collision with ev-A, ev-B, ev-C or vm2-1 (`evm-repeat.md`). Every trajectory is one-shot.
  What repeats is the launch OPPORTUNITY: about every six years, "prior to every third opposition of Mars"
  (p.75 [83]).
- **What it gives the project:**
  - the origin of the Hollister group's Venus work and of its units and flyby conventions (EMOS; distances in
    planet radii; speeds scaled by the circular surface speed);
  - a clean, self-consistent patched-conic itinerary with full v_inf vectors at both Venus flybys. All the
    parameters I could check reproduce to 0.001-0.002 EMOS (sec. 3);
  - the earliest held statement of the "every third opposition / about six years" E-V-M opportunity cycle (sec. 4).
- **Catalogue implication (PROPOSAL only):** none. It is not a cycler, a quasi-cycler or a precursor that
  inserts into one. It could be cited as history (with Crocco 1956 and Minovitch 1963) in any note on the
  origins of Venus-assisted Mars trips.

## 1. Content, chapter by chapter (READ)

- **Ch. 1 Introduction (pp.1-9 [9-17]):**
  - The reasons to go: Mars as the only likely life-bearing planet besides Earth.
  - Prior mission concepts:
    - von Braun 1953 and Ley-von Braun 1956: Hohmann, 260 d each way, 450-d wait;
    - Himmel et al. 1961: a 400-d nuclear mission, with 70 t of passive shielding;
    - the NASA-Marshall EMPIRE contractor studies: Aeronutronic [1] studied dual-planet flybys "of the Crocco(5)
      and symmetric types(1)" (p.6 [14]); Lockheed [4]; General Dynamics [6].
  - The priority claim above.
  - The approach: build simple models first, draw conclusions, then check them on the computer.
- **Ch. 2 Orbital geometry (pp.10-12 [18-20]):**
  - Units: a.u., planet radius, Julian date, and EMOS (Earth mean orbital speed, "about 100,000 feet per
    second").
  - Periods: Mars 687 d, Earth 365 d, Venus 224 d. Synodic periods: "about 780" d (E-M), "about 580" d (E-V) and
    "about 340" d (V-M alignments) (p.12 [20], image).
  - Fig. 2 [129] plots the E-M, E-V and V-M alignment dates for 1964-1980, with tick marks at 1965, 1971 and
    1978.
- **Ch. 3 Orbital parameters of all possible transfers (pp.13-20 [21-28]):** a "V_x-V_y plane" chart. Each point
  is a departure v_inf in the planet frame, and loci of constant period, semi-latus rectum, eccentricity,
  flight time, trip angle and launch/arrival configuration are plotted on it. It follows Vertregt 1958 and
  Dugan 1960. "All combinations of flyby and stopover missions may be analyzed by pairing inbound and outbound
  velocity vectors at the destination planet" (p.20 [28], OCR).
- **Ch. 4 Simplified approximate model (pp.21-28 [29-36]; App. A):**
  - Linearised motion about a circular-orbit Earth in a rotating frame (Hill-Clohessy-Wiltshire type, in a.u.,
    years and EMOS): q = 4 sin 2pi t - 6 pi t, r = 2 - 2 cos 2pi t, s = sin 2pi t.
  - Superposition of impulses. This comes from his 1959 S.M. thesis on rendezvous [13].
- **Ch. 5 Round trip (pp.29-36 [37-44]):**
  - With Dugan's velocity-versus-trip-time curve and a payload-weight model, the rocket equation with a nuclear
    c = 30,000 ft/s gives a minimum-weight trip of about 400 days.
  - Conclusion: "about 400 days' duration or less".
- **Ch. 6 Paths in space (pp.37-41 [45-49]):** Earth- and Mars-centred frames (Figs. 8-10).
  - Fast (under 400 d) trips arrive near opposition.
  - The best 400-d pair has one leg that launches almost straight sunward, with a trip angle of about 270 deg and
    a perihelion near Venus's orbit. It is "uneconomical", but needed so that the two legs balance the drift
    ahead of and behind Earth.
- **Ch. 7 Bi-elliptical transfer (pp.42-47 [50-55]):**
  - Replace the uneconomical leg by two ellipses with a common perihelion m and a burn there (Barrar 1963:
    two impulses beat one).
  - The angle mSM (120, 135 or 150 deg in the program) sets the Mars arrival speed. This can be matched to
    the drag-brake limit.
  - The burn would be better made deep in a planet's well: "the bi-elliptical transfer would be even more
    effective if the second orbital transfer could be accomplished during an encounter with Venus" (p.46 [54]).
    The geometry recurs with "the synodic period of their relative positions being about six years".
  - The Venus mission "was conceived as an extension of the bi-elliptical transfer" (p.47 [55]).
- **Ch. 8 Planetary encounter (pp.48-56 [56-64]):**
  - Patched conics.
  - Vc (circular speed at the surface) = 0.266 (Earth), 0.243 (Venus), 0.121 (Mars) EMOS (image).
  - Equations, read on the image:
    - r_pi = (Vc/Vs)^2 (8.1);
    - V_E^2 = 2 Vs^2 (8.2);
    - V_pi^2 = V_E^2 + V_H^2 (8.3);
    - sin delta = 1/(1 + (V_H/Vs)^2) (8.4);
    - e = 1 + (V_H/Vs)^2 (8.5).
  - Parking-orbit burn dV = V_pi - Vs (8.8), at parking radii of 1.1 (Earth), 1.1 (Venus) and 1.3 (Mars).
  - Powered flyby: inbound and outbound hyperbolas share a periapsis, with an impulse there (8.9-8.13).
  - Rule of thumb: the turn is at most 60 deg when V_H > Vc.
  - "Most agencies studying flyby missions" treat only the unpowered case.
- **Ch. 9 Atmospheric braking (pp.57-67 [65-75]):**
  - Chapman corridor theory.
  - Table 9.1 (image):

    | Vehicle | L/D | Altitude | g | Heat protection | W/(C_D A) |
    |---|---|---|---|---|---|
    | drag brake | 0.2 | 50 n mi | 6 | reradiation | 5 psf |
    | Apollo type | 0.5 | 40 n mi | 8 | ablative shield | 50 psf |
    | hypersonic glider | 2 | 40 n mi | 4 | ablative nose | 500 psf |

  - Mars capture uses a Lockheed rocket-augmented drag brake. Mass ratio M_i/M_f = 1.9 (V_pi - 0.04)^(1/4)
    (9.1), from Cottrell & Olson [30].
  - Limits (p.67 [75]):
    - Mars entry V_pi up to 0.34 EMOS ("about three times" Vc), that is V_H 0.3.
    - Earth entry V_pi 0.45 EMOS ("about 1.8 times" Vc), that is V_H 0.25. Above that, a rocket brakes the
      vehicle to 0.45.
- **Ch. 10 Use of the computer (pp.68-73 [76-81]):**
  - FORTRAN "General Trajectory Program": a Lambert solver following Battin, with mean elements for the epoch
    1971 Feb 18.0 (JD 2441000.5) taken from the 1964 American Ephemeris and the 1963 Annuaire.
  - A bi-elliptical program (App. E listing).
  - "All practical trips between planets have been computed at ten-day intervals through 1980". The print-out
    "read[s] very much like a train schedule".
  - Multi-leg trips are built by pairing legs at a planet with the Ch. 8 encounter analysis.
- **Ch. 11 Transfer via Venus (pp.74-82 [82-90]).** See sec. 2.
- **Ch. 12 Recommended mission (pp.83-91 [91-99]).** See sec. 2.
  - Alternatives (p.89 [97]): if 1971 is missed, transfer via Venus is possible again in six years, but not as
    attractively until after 2000.
  - In other oppositions, use one bi-elliptical leg ("about 0.03 EMOS or 3000 feet per second").
  - The 1973 or 1975 windows "could double or triple" the launch weight.
  - Manned flybys are judged poor value. App. D gives a single Mars flyby example.
- **Ch. 13 Summary and conclusions (pp.92-97 [100-105]):**
  - The most economical round trips make a trajectory change at a close Venus approach. This is possible
    "every six years or during every third opposition period of Mars" (p.94 [102], OCR).
  - Navigation must be extended to handle powered flybys (Stern [19], Battin [34]).
  - Ross [37] suggested not burning exactly at closest approach, so as not to disturb observations
    (p.96 [104], image).
  - EMPIRE should be accelerated.
- **App. B:** accuracy of the models. Lawden [40]: neglecting Jupiter costs about 5 ft/s.
- **App. C:** a calculus-of-variations optimisation in the simple model. It shows no gain from a separate
  90-deg plane change.
- **App. D (pp.110-112 [118-120]):**
  - A 1.5-yr Earth-Mars-Earth free return in 1970-71. In the simple model it is the root t = 1.41 yr of
    4(1 - cos 2pi t) = 3 pi t sin 2pi t (D.2).
  - Table D.1, image:

    | Parameter | Simple model | Accurate two-body |
    |---|---|---|
    | launch JD | 2440920 | 2440930 |
    | launch speed | 0.304 | 0.238 |
    | V_x, V_y, V_z | -0.298, 0.045, 0.047 | -0.234, 0.021, 0.037 |
    | Mars JD | 2441175 | 2441180 |
    | Earth return JD | 2441430 | 2441440 |
    | return speed | 0.304 | 0.295 |
    | V_x, V_y, V_z | 0.298, 0.045, -0.047 | 0.294, 0.005, -0.033 |

  - He notes that Lockheed [4] calls these "symmetric" round trips and needs a computer to find them. The simple
    model finds them by hand.
- **App. E:** the FORTRAN listing and sample output of the bi-elliptical program. I did not transcribe it.

## 2. The Venus-swingby mode (Chs. 11-12)

- **Availability rule (p.75 [83], image):**
  - A low-energy transfer has the launch-arrival alignment about midway through the transfer. So E-V-M needs the
    E-V, V-M and M-E alignments "in that order with a spacing of about 140 days between them".
  - This happens "about every six years or prior to every third opposition of Mars". When it is right for E-V-M
    before an opposition, it is right for M-V-E after it.
  - Examples: 1965, 1971 and 1978. 1971 is best, because it is near Mars perihelion. "An orientation which
    allows an Earth-Venus-Mars transfer will not occur in so favorable a position again until after the turn of
    the century."
- **Selection (p.76 [84], image):**
  - Drop legs with Earth V_H > 0.2, Mars V_H > 0.3, or perihelion < 0.7 a.u.
  - Keep Venus dates that have both legs, then pair them on the Venus V_x-V_y plane.
  - In 1971 the window is 120 d at Earth, 50 d at Venus (OCR misreads this as 60) and 120 d at Mars.
- **Unpowered case (p.77 [85]):** equal |v_inf| is needed. That is "the type of mission reported by the agencies
  studying dual-planet 'flyby' missions" [1][4], and its window "extends for only a day or two".
- **Typical trips, Table 11.1 (p.80 [88]) and text pp.78-79 [86-87], all read on the image.** Units are EMOS.
  - Outbound:
    - Earth 1 Sep 1970 (JD 2440830), longitude 338 deg, V_H 0.126, burn 0.126.
    - Venus 31 Dec 1970 (2440950), longitude 127 deg, at 1.74 Rv on the southern sunlit side. In
      (-0.027, 0.084, -0.146) |0.171|; impulse +0.024; out (0.027, 0.208, -0.020) |0.211|. The vector is
      "rotated 57 degrees".
    - Mars 19 May 1971 (2441090), longitude 267 deg, V_H 0.284. Drag brake, mass ratio 1.42.
  - Return:
    - Mars 6 Oct 1971 (2441230), longitude 353 deg, V_H 0.196, burn 0.141.
    - Venus 3 Apr 1972 (2441410), longitude 144 deg, at 1.2 Rv on the southern sunlit side. In
      (-0.196, 0.154, -0.089) |0.265|; impulse 0.045 SUBTRACTED; out (-0.095, 0.106, 0.121) |0.186|.
    - Earth 22 Jul 1972 (2441520), longitude 299 deg, V_H 0.226. Direct entry.
  - The Mars arrival is "84 days prior to" and the departure "56 days after" the 1971 opposition (p.77 [85]).
  - A drag-brake pass at Venus is "conceivable (but not recommended)".
- **Advantages (p.79-82 [87-90]):**
  - The single-ellipse trip arriving on the same day needs Earth V_H 0.366, which is a burn of 0.26, "almost
    twice".
  - The typical saving is "about 0.07 EMOS" for equal duration: 25 % in initial weight (nuclear), and half
    (chemical).
  - The return via Venus is "not as spectacular".
  - "Three and a half months to Venus and four and a half months to Mars instead of eight months."
- **Programme (Ch. 12):**
  - The optimistic option is 140 days at Mars and a return via Venus: 690 d in total, with "three close approaches
    to Venus" across two vehicles.
  - Early-return trips, Table 12.1 (p.85 [93], image). Mars departure, Earth arrival, (V_H; burn):

    | Stay | Mars departure (V_H; burn) | Earth arrival (V_H; burn) | Total |
    |---|---|---|---|
    | 0 d | 19 May 1971 (0.170; 0.121) | 27 Aug 1971 (0.235) | 360 d |
    | 10 d | 29 May 1971 (0.171; 0.122) | 6 Sep 1971 (0.198) | 370 d |
    | 40 d | 28 Jun 1971 (0.180; 0.128) | 25 Nov 1971 (0.177) | 450 d |
    | 90 d | 17 Aug 1971 (0.185; 0.132) | 24 Mar 1972 (0.276; 0.012) | 570 d |

  - Second launch, Table 12.2 (p.86 [94], image):
    - Earth 18 Jun 1971 (V_H 0.167, burn 0.142);
    - Mars 26 Sep 1971 (0.236);
    - 10 days at Mars;
    - then the same return as the first vehicle. "Total mission = 400 days".
  - The 400-d trip of the abstract is this one, shown in Fig. 18. Fig. 17 is the 370-d trip: out via Venus, back
    direct.
  - The single-ellipse alternative to the itinerary would need +45 % (nuclear) or x3 (chemical) initial weight.
    The best 350-d single-ellipse trip needs +25 % or about x2 (p.88 [96], image; the OCR reads "15 per cent").

## 3. Checks (`checks.py` / `checks.out`)

Inputs are image-read. 1 EMOS = 29.785 km/s = 97,720 ft/s is my value.

- **Julian dates:** 14 of 15 printed dates agree with the calendar dates.
  - **One disagreement: JD 2440950 is 30 Dec 1970 (0h on JD 2440950.5), but it is printed "31 Dec. 1970"** in
    Table 11.1 and the text. All the other dates use the JD n.5 = date convention. This is probably a one-day
    calendar slip. The legs use the JD grid (10-day steps), so no other number is affected.
- **Leg days:** E-V 120, V-M 140, Mars stay 140, M-V 180, V-E 110 (690 in total).
  - Tables 12.1 and 12.2: totals 360, 370, 450, 570 and 400, and stays 0, 10, 40, 90 and 10. All agree.
  - p.86: 290 d and 170 d before the second launch, which is 30 d after the first arrival. All agree.
  - Opposition: 84 d before the Mars arrival plus 56 d after the departure gives JD 2441174. App. D uses
    2441175. The two agree within a day.
- **Vector magnitudes:** 0.1706, 0.2107, 0.2647 and 0.1868, against the printed 0.171, 0.211, 0.265 and 0.186.
  - The outbound angle is 56.8 deg (printed 57). The return angle is 60.7 deg (not printed).
- **Powered flyby with eqs. 8.1-8.5 (Vc = 0.243):**
  - Outbound at 1.74 Rv: dV 0.0236 (printed 0.024), turn 58.1 deg (vectors 56.8).
  - Return at 1.2 Rv: dV 0.0460 (printed 0.045), turn 60.3 deg (vectors 60.7).
  - The 1-1.3 deg turn mismatches are within the 3-decimal rounding of the vectors and of the radius. Not a
    disagreement.
- **Parking-orbit burns (eq. 8.8)** reproduce 0.126, 0.142, 0.141, 0.121, 0.122, 0.128, 0.132, and 0.26 for the
  0.366 single-ellipse case, to 0.001.
- **Rocket equation:** +45 % nuclear means 0.114 EMOS, which gives x3.05 chemical. +25 % means 0.069 EMOS
  ("about 0.07"), which gives x1.95 chemical ("half", "about twice"). The weight claims agree.
- **Rules:** the 60-deg rule gives 2 x 30 deg exactly. Mars 0.34 EMOS gives V_H 0.294 ("0.3"), and
  0.34/0.121 = 2.81 ("about three"). Synodic periods from his own periods: 779 d (780) and 580 d (580). Three E-M
  synodic periods are 2340 d (6.41 yr, "about every six years"). The 1965, 1971 and 1978 oppositions are each
  three oppositions apart.
- **Table D.1:** the launch and return magnitudes from the components agree (0.3050/0.304, 0.2378/0.238,
  0.2959/0.295).
  - D.2 root t = 1.4067 yr, which gives q = -24.30, r = 3.67, s = 0.553 and V_x/V_y = -6.63. Printed: 1.41,
    -24.35, 3.67, 0.55, -6.67. These agree to his rounding.
  - The day spans (255/510 d simple, 250/510 accurate) against 0.705 and 1.41 yr (257.5 and 515 d). These agree
    roughly.
- **Disagreements, all UNRESOLVED.** I make no misprint claim. In each case my own model or reading may be the
  cause.
  1. **"31 Dec. 1970" for JD 2440950** (see above).
  2. **Mars mass ratio:** 1.42 is printed and clear on the image. Eq. 9.1 with V_pi from V_H 0.284 gives
     1.39-1.40 for an entry radius of 1.0-1.1 Rm. He may have read it off Fig. 16 rather than the equation.
  3. **Table 12.1, 90-d return, Earth burn 0.012:** if the burn only slows the vehicle to the 0.45 EMOS entry
     limit, V_H 0.276 needs 0.017 at r = 1. An entry radius slightly above 1 narrows the gap. His exact rule for
     this burn is not stated.
  4. **V-M synodic period "about 340" d:** his own periods give 332 d, and the true value is about 334 d. This is
     probably loose rounding.
  5. **Earth limit "about 1.8 times" Vc:** 0.45/0.266 = 1.69.
  6. **App. D, V_z = 0.047 (text and table):** z = s V_z at t = 0.705 with z = -0.04 gives 0.042. The printed
     |V| 0.304 fits 0.042 (0.3043) a little better than 0.047 (0.3050). I did not find his exact Mars z or t.
  7. **"Three and a half months to Venus":** the E-V leg is 120 d, which is 3.9 months. (The V-M leg is 140 d,
     4.6 months, against "four and a half".) This is a loose phrase.

## 4. Lineage: how this thesis leads to the Hollister group's periodic orbits

What is READ in held files is marked; what is my interpretation is marked INFERRED.
- **Hollister 1963 (this file):**
  - powered Venus flybys inside E-M round trips;
  - the EMOS unit and Vc-scaled flyby equations;
  - leg pairing at a planet by v_inf vectors;
  - the "every third opposition" E-V-M rule;
  - a 10-day-grid Lambert "train schedule" through 1980.
- **Hollister 1964, AIAA 64-647, "Predicting launch dates for Mars transfer via Venus":** NOT held (wanted row
  43). By its title, it is the follow-on to the Ch. 11 availability rule (INFERRED).
- **Hollister & Prussing 1965, AIAA 65-700 (held; digest 2026-10-05):**
  - It uses the same EMOS and Vc-normalised flyby units, and cites this thesis as [5].
  - It optimises the impulse that connects two given hyperbolas. The common-periapsis impulse used here is
    "within 3 % of the optimum".
  - It widens the rule: Venus "is available on one leg during every Mars opposition period", and on both legs in
    "every third".
  - It finds that where a pure flyby is possible, thrust at Venus saves only "a few hundred feet per second".
  - This does not contradict the thesis's 0.07 EMOS. Here the leg dates are fixed by the round trip, so the
    inbound and outbound v_inf differ (0.171/0.211 and 0.265/0.186). The impulse is forced by the date pairing.
    The saving is measured against a single-ellipse trip, not against a nearby unpowered flyby (INFERRED).
- **Sohn 1964 and Gillespie & Ross 1967 (both held):** one-shot unpowered Venus swingbys. G-R 1967 cites this
  thesis. The 2026-10-05 G-R digest (sec. 2) calls G-R "the first published statement (to my knowledge, INFERRED)
  of the 2338-d / 32-yr E-V-M commensurability".
  - **Proposed qualification:** this thesis states the about-6-yr, every-third-opposition E-V-M opportunity cycle
    four years earlier (p.46 [54], p.75 [83], p.94 [102]). It does not give the 2338-d value or the 32-yr cycle.
    So G-R stays the first held source for those two, and this thesis is the first held source for the 6-yr rule.
- **Hollister 1969, JSR 6(4) (held):** the first published Earth-Venus periodic orbits. Its Orbit I has
  v_inf 0.126 EMOS, the same number as this thesis's Earth departure. The legs differ (Orbit I: 0.485 yr; here:
  120 d), so this is a numerical coincidence. I build no lineage on it.
- **Menning 1968 S.M. thesis and Hollister & Menning 1970 (held):** the 15 E-V periodic orbits, built by pairing
  legs at each flyby with equal |v_inf|.
  - In this thesis the equal-|v_inf| unpowered flyby is the "restrictive" special case (p.77 [85]).
  - The later periodic-orbit work makes that special case the whole construction, and repeats it (INFERRED).
- **Rall 1969 Sc.D. thesis (held; Hollister chairman) and Rall & Hollister 1971 (held):** the Earth-Mars periodic
  orbits. Rall's sec. 4.2/4.4 records failed E-V-M and V-M periodic attempts (2026-06-07 mining note). The
  E-V-M chain that this thesis flew once was later tried as a repeating orbit and did not close (READ in that
  note).
- **Ross 1963 [37]:** cited here only on the timing of flyby burns. It is not cited for symmetric round trips,
  which the thesis credits to Lockheed [4] (App. D, p.111 [119]).

## 5. Proposed corrections to held notes

- `2026-10-05-digest-hollister-prussing-1965-optimum-transfer-mars-via-venus.md` sec. 4 item 1, and
  `2026-10-05-digest-gillespie-ross-1967-venus-swingby-mission-mode.md` sec. 4 (the "Hollister 1963 is not held"
  line and item 1): mark Hollister 1963 as received (this filing).
- The G-R digest sec. 2 priority sentence: qualify it as in sec. 4.
- The H-P digest sec. 4 item 2 gives Hollister 1964's title as "Mars Transfer via Venus". The wanted list (row 43)
  already corrects it to "Predicting launch dates for Mars transfer via Venus" (refs-sonnet).
- Wanted list row 42: mark it received and remove it.

## 6. Citation mining (references 1-40, pp.145-149 [153-157], all read on the image)

Held status: checked with `ls cyclers_pdf/papers | grep -i <author>` and `grep -i <author> CORPUS_INDEX.md`.
Wanted-list rows use the numbering of the file as read on 2026-10-07. Watch for these false hits:
- `luidens-1964-mars-nonstop-round-trip-trajectories` is NOT refs 26-27 (Luidens 1961 TNs);
- `stern-tapley-...-2020` is NOT R. G. Stern 1963;
- `putnam-braun-2005` is NOT von Braun.

| ref | work | status |
|---|---|---|
| 1 | Aeronutronic (Ford), "EMPIRE - A Study of Early Manned Interplanetary Missions", Publ. U-1951, 21 Dec 1962 | not held; not on the list. **New candidate (R1(c) history, low):** it treats the Crocco and symmetric E-V-M-E flyby types |
| 2 | de Vaucouleurs 1954, Physics of the Planet Mars | not held; not relevant |
| 3 | Taylor & Blockley 1959, "Crew Performance in a Space Vehicle" | not held; not relevant |
| 4 | Lockheed, Final Report "A Study of Interplanetary Transportation Systems", 3-17-62-1, 2 Jun 1962 | not held; not on the list. Source of the "symmetric" round-trip name (App. D). Low |
| 5 | Crocco, G. A. (1956), "One Year Exploration Trip", Proc. VII IAC, Rome | not held; **wanted row 10** |
| 6 | General Dynamics/Astronautics, "A Study of Early Manned Interplanetary Missions", AOK 63-0001, 31 Jan 1963 | not held; not on the list. Low |
| 7 | von Braun 1953, The Mars Project | not held; history |
| 8 | Ley & von Braun 1956, The Exploration of Mars | not held; history |
| 9 | Himmel, Dugan, Luidens & Weber 1961, Aerospace Engineering 20(7) | not held; not relevant |
| 10 | Finney 1962, New York Times | not relevant |
| 11 | Vertregt 1958, "Interplanetary Orbits", JBIS 16(6) | not held; method background. Low |
| 12 | Dugan 1960, NASA TN D-281 | not held; round-trip parameter charts. Low |
| 13 | Hollister 1959, MIT S.M. thesis (rendezvous control) | not held; source of the simple model. Low |
| 14 | Lawden 1959, "Interplanetary Rocket Trajectories", Adv. Space Sci. 1 | not held; optimisation background. Low |
| 15 | Long 1960, Astronautica Acta 6(2-3) | not held. Low |
| 16 | Johnson & Smith 1960, "Round-Trip Trajectories for Mars Observation", Adv. Astronaut. Sci. 5 | not held. Low |
| 17-18 | Wallner & Kaufman 1961 (TN D-681); Kash & Tooper 1962 | not held; shielding, not relevant |
| 19 | Stern, R. G. (1963), MIT Sc.D. thesis, midcourse guidance | not held; not relevant |
| 20 | Hayes, Hoffman & Stuhlinger 1961 | not relevant |
| 21 | Barrar 1963, AIAA J 1(1), two-impulse vs one-impulse | not held. Low |
| 22 | Zee 1963, AIAA J 1(1), finite thrust | not held; not relevant |
| 23-28 | Eggers et al. 1957; Chapman 1959 x2; Luidens 1961 x2; Hayes & Vander Velde 1962 | not held; entry, not relevant |
| 29-30 | Ragsac & Titus 1962 (ARS); Cottrell & Olson 1963 (Lockheed) | not held; mission optimisation and drag-brake mass. Low |
| 31-35 | 1964 American Ephemeris; 1963 Annuaire; Moulton 1914; Battin 1963; IBM FORTRAN | data and texts |
| 36 | Barrett & Lilley 1963, Sky and Telescope | not relevant |
| 37 | Ross, S. (1963), "A Systematic Approach to the Study of Nonstop Interplanetary Round Trips", AAS 9th Annual Meeting, Jan 1963 | not held; **wanted row 4** |
| 38 | Hildebrand 1952 | text |
| 39 | Lockheed, "Early Manned Interplanetary Mission Study [EMPIRE]", 8-32-62-1, 4 Dec 1962 | not held; not on the list. Low |
| 40 | Lawden 1959, "Manned Navigation and Guidance in the Solar System", JBIS 17(5) | not held; not relevant |

- No cited work is held.
- The only new candidate worth listing is [1], the Aeronutronic EMPIRE report, as R1(c) history: it is a contractor
  study of Crocco-type E-V-M-E flybys, and it is not periodic. [4] and [39] (Lockheed) are lower.
- The Hollister group's own follow-ons (Hollister 1964, Hollister & Prussing 1965/1966) are not cited, because
  they came later.

*Proposed filing: the OCR copy as `cyclers_pdf/papers/hollister-1963-mission-manned-expedition-mars-mit-scd-thesis-dspace-1721.1-15698.pdf`,
with `checks.py`, `checks.out` and `evm-repeat.md` beside it as `<pdf stem>-<file name>`.*

*Filed as `cyclers_pdf/papers/hollister-1963-mission-manned-expedition-mars-mit-scd-thesis-dspace-1721.1-15698.pdf`. Check scripts, outputs and other files named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`.*

*Wanted-list row numbers in this digest are the batch-34 numbering; the list was renumbered in batch 35.*
