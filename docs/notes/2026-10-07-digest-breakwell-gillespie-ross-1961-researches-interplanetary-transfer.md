# Digest: Breakwell, Gillespie & Ross 1961, "Researches in Interplanetary Transfer" (ARS Journal) (#960)

J. V. Breakwell, R. W. Gillespie and S. Ross (Lockheed Missiles and Space Division, Sunnyvale), "Researches in
Interplanetary Transfer", ARS Journal 31(2):201-208, February 1961. doi 10.2514/8.5428 (printed on every page margin).
"Presented at the ARS 14th Annual Meeting, Washington, D.C., Nov. 16-20, 1959" (p.201 footnote).
- File given: `5f61a330-breakwell1961.pdf`, 8 pages, text layer (AIAA reprint, downloaded 2015). md5 f97b4fd8642aa5c2be917efb8f7625e1.
  The paper ends on p.208 (printed p.207 plus the start of p.208). The rest of p.208 is the first page of R. V. Warden,
  "Ballistic Re-Entries With a Varying W/CDA". That article is not part of this digest.
- **Proposed corpus filename:** `cyclers_pdf/papers/breakwell-gillespie-ross-1961-researches-interplanetary-transfer-ars-j-31-2-201-doi-10.2514-8.5428.pdf`
- **Wanted list:** row 48, "Breakwell, Gillespie & Ross (1959), ARS paper 954-959 (journal form: ARS J. 31(2):201-208, 1961)".
  The journal prints **no preprint number**, only the meeting. Ref. 4 is cited in the ARS style "ARS preprint 870-59", so the row's
  "954-959" probably means "954-59". This file cannot confirm the number. The journal form is in hand, so the row can close.
- **How I read it.** I viewed all 8 pages on 110-dpi renders (pp.201-206 and 208 in full; p.207 as a 250-dpi crop of the
  E(1) = E(2) passage). I read the prose from the text layer and checked the header, the reference list and the quoted numbers
  against the images. I read the labelled minimum points of Figs. 3, 5 and 7 on 250-dpi crops. Arithmetic checks:
  `checks_bgr1961.py`, output in `checks_bgr1961.out`.

## 0. Verdict

**A method-and-charts paper: fast patched-conic (zero-SOI) Lambert solutions and early contour ("pork-chop") maps of
Earth-Mars and Mars-Earth hyperbolic excess speed against departure date and flight time, for 1960-61.** They are not the first:
"Fig. 3 is similar in nature to Fig. 15 of (3)" (Gunkel et al. 1959; p.202).
- **Focus questions, answered plainly:**
  - **Mars or Venus round trips:** none worked. Round trips are deferred: "These trips will be investigated in a paper now in
    preparation" (p.204). That points to Ross 1963 (wanted row 5) and the later Gillespie & Ross 1967.
  - **Venus:** not in the paper at all.
  - **Free returns:** none.
  - **Symmetric or repeating trajectories:** none as missions. Two structural results come close:
    1. The 180-deg ridge lines repeat once per synodic period, with slope d(dt)/dt1 = n1/n2 - 1 (p.202).
    2. Appendix (p.207): at the crossing of two Lambert branches, the two transfer arcs together "form one complete ellipse in
       space". I show below that this ellipse has period exactly 2 dt. Planet motion is not included, so this is geometry, not a round trip.
- **What it gives the project:** history and the method origin of the Lockheed group's charts.
  - The speeds are in units of Earth's mean orbital speed. This is the "EMOS" unit of Gillespie & Ross 1967.
  - It counts multi-revolution Lambert solutions: exactly two for less than one revolution, and up to four more for each extra revolution (pp.201, 206).
  - The low-energy Earth-Mars minima for 1960 are in sec. 1.
- **Catalogue implication (PROPOSAL only):** none. No periodic or repeating trajectory. No #943 or #942 bearing.

## 1. Content (READ)

- **Method (pp.201, 205-208).**
  - Three one-centre conics. The heliocentric arc is extended to the planet centres, and each planet hyperbola meets it at its asymptote (zero-radius sphere of influence).
  - Planets move on fixed, mutually inclined ellipses with 1960 elements.
  - A "modified form of Lambert's theorem" gives every heliocentric arc between the two dates. The paper gives the formulas
    for cases 1A, 1B, 2A, 2B and 1H, 2H, with m whole revolutions [A-2]-[A-8]. It proves, from the monotonicity of T(E),
    that there are exactly two solutions for m = 0 (pp.205-206).
  - The excess-velocity vector is resolved into right ascension and declination in the planet's equatorial frame [A-14].
  - About 0.1 s per orbit on an IBM 7090 (Fortran).
- **Charts (Figs. 1-8; contour step 0.1 of Earth's mean orbital speed):**
  - Fig. 1, circular coplanar orbits, Earth to Mars, departure March 1960-April 1961.
    - The "Hohmann" point is 1 Oct 1960, 259 d, 0.099 (image).
    - The 180-deg lines are thin ridges. Their slope is "about 0.88" for Earth to Mars and -0.468 for Mars to Earth.
    - There are also ridges at 0/360 deg.
  - Fig. 2, eccentric coplanar orbits: the minimum moves to 26 Sep 1960, 360 d, 0.117, with a transfer angle of about 210 deg.
    This is "an increase of 18 per cent" over the Hohmann value.
  - Fig. 3, eccentric and inclined orbits (image).
    - The Hohmann point splits into two local minima at about +-5 d from 1 Oct 1960: **0.118 on 26 Sep 1960** (transfer angle about 220 deg, almost in the ecliptic)
      and **0.147 on 5 Oct 1960**, which the text gives as a 209-d trip.
    - The nodal-transfer points are marked at 26 Sep 1960 and 13 Oct 1961.
    - The text also says the 1966-68 minimum is 0.098, "smaller ... than the Hohmann case", because Fig. 1 used Mars' mean radius.
  - Fig. 4: long transfer times, 600-1000 d, where the contours become nearly parallel to the ridges.
  - Fig. 5, Earth to Mars arrival speeds (image): minima **0.080 on 24 Aug 1960** and **0.079 on 14 Oct 1960**.
  - Fig. 6: Mars to Earth departure speeds. Nodal transfer on 11 May 1961; no labelled minimum.
  - Fig. 7, Mars to Earth arrival speeds (image): minima **0.098 on 5 Aug 1960** and **0.096 on 16 Sep 1960**.
  - Fig. 8: Earth to Mars retrograde launches (circular, coplanar), 1.50-2.40.
- **The round-trip remark (pp.203-204).** Short outbound trips arrive at Mars slowly, but a matching low departure speed on the same day exists
  only for long return trips. "Thus, the choice of either a short trip going or a short return trip, but not both, is possible in such cases."
  This is about a zero-stay round trip, not a free return.

## 2. Checks (`checks_bgr1961.py` / `.out`)

- Ridge slopes: n_E/n_M - 1 = 0.8808 ("about 0.88"). n_M/n_E - 1 = -0.4683 (-0.468). Both agree.
- Hohmann transfer, a_M = 1.5237 AU: 258.9 d (259 d). The departure excess speed is 0.0989 of Earth's mean speed (0.099). Both agree.
- 0.117/0.099 = 1.182 ("18 per cent"). (0.118 - 0.117)/0.117 = 0.85% ("only 1 per cent"). Both agree.
- **Appendix result (p.207, image):** at the crossing of the 1B and 2A branches, E(1) = E(2) = -(pi/T)^(2/3), independent of m.
  - With E = -s/2a and T = n dt/(s/2)^(3/2) (p.205), this gives a = (n dt/pi)^(2/3) in AU.
  - So the period is 2 pi a^(3/2)/n = 2 dt. The script confirms 2.000000 dt for three cases.
  - So "one complete ellipse" means one closed orbit through both points. The spacecraft flies one arc out and the other arc back, each in dt.
    The paper does not draw this conclusion.
- Not checked: the contour charts themselves (no table), and the right-ascension/declination transform.

## 3. Relation to held work

- **Gillespie & Ross 1967** (JSR 4(2):170, HELD, `gillespie-ross-1967-venus-swingby-...`; digest `2026-10-05-digest-gillespie-ross-1967-venus-swingby-mission-mode.md`).
  - It shares two authors, the unit (EMOS, "Earth mean orbital speed") and the 0.1-unit contour charts.
  - It adds what this paper lacks: Venus swingbys, the 2338-d syzygistic period, the symmetry principle from Ross 1963, and round trips.
  - This 1961 paper is the chart method that the 1967 paper's Figs. 3-12 build on (INFERRED; the 1967 paper does not cite it).
- **Ross 1963** (wanted row 5, not held): the most likely form of the "paper now in preparation" on round trips (INFERRED from the
  author list and the topic; the 1961 text gives no title).

## 4. Citation mining (9 references, p.208)

Checked with `ls cyclers_pdf/papers | grep -i` and `grep -i docs/notes/CORPUS_INDEX.md`, and against the wanted list.

| ref | work | status |
|---|---|---|
| 1 | Plummer 1918, An Introductory Treatise on Dynamical Astronomy, ch. V | not held; textbook (Lambert's theorem) |
| 2 | Battin, R. H. (1959), "The Determination of Round-Trip Planetary Reconnaissance Trajectories", J. Aero/Space Sci. 26(9):545-567 | not held. **On wanted row 42 ("Battin 1959/1999").** By its title (not read), an early round-trip reconnaissance paper. It is the one reference here that may bear on free returns. Suggest the row names it in full. |
| 3 | Gunkel, Lascody & Merrilees (1959), "Impulsive Midcourse Correction of an Interplanetary Transfer", Proc. 10th IAC, London | not held; not on the wanted list. Its Fig. 15 is the earlier version of Fig. 3 here. Low priority (history of pork-chop charts). |
| 4 | Karrenberg & Arthur (1959), "Interplanetary Ballistic Orbits", ARS preprint 870-59 | not held; not on the wanted list. Low priority. |
| 5 | Ehricke (1960), "Zur Auswahl von Flugbahnen für bemannte Raumfahrzeuge zu den Planeten Mars und Venus", Raketentechnik und Raumfahrtforschung 4(1):16-22 | not held. Not on the wanted list (rows 47 and 50 are other Ehricke items). Low priority; manned Mars/Venus trajectory choice. |
| 6 | Vertregt (1958), "Interplanetary Orbits", JBIS 16(6):326-354 | not held; not on the wanted list. Low priority. |
| 7-9 | American Ephemeris 1960; Planetary Coordinates 1960-1980 (HMSO); Astron. Ezheg. SSSR 1960 | ephemerides; not needed |

- No new high-priority candidates. Battin 1959 is already covered by row 42.
- **Proposal for row 48:** mark it received (journal form; the preprint number is not printed); remove it.

*Filed as `cyclers_pdf/papers/breakwell-gillespie-ross-1961-researches-interplanetary-transfer-ars-j-31-2-201-doi-10.2514-8.5428.pdf`. Check scripts, outputs and other files named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`.*

*Wanted-list row numbers in this digest are the batch-30 numbering; the list was renumbered after batch 34.*
