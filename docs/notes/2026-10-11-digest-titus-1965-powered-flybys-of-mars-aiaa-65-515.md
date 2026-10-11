# Digest: Titus 1965, "Powered Flybys of Mars", AIAA Paper 65-515 (#960)

Richard R. Titus, Research Laboratories, United Aircraft Corporation, East Hartford, Connecticut. AIAA Paper No. 65-515,
AIAA Second Annual Meeting, San Francisco, 26-29 July 1965. DOI 10.2514/6.1965-515. 18 printed pages of text and 23 figures.
Wanted-list row 56 (Titus is the first item; Wood et al. 1964, NASA TM X-53055, stays wanted).

- **Filed as** `cyclers_pdf/papers/titus-1965-powered-flybys-of-mars-aiaa-65-515-doi-10.2514-6.1965-515.pdf`. Table
  transcription: `data/sources/titus-1965-tables.yaml`.
  - The filed file is an `ocrmypdf --redo-ocr -l eng --rotate-pages -O1` copy of the received scan, md5
    39880f88938228c7ec5261169bf65e94 (554,745 bytes), 16 pdf pages (two printed pages per pdf page, 783 x 610 pt).
    The received scan has md5 7aab3cf16abeb70160c67e210ddede3e (529,072 bytes) and almost no text layer.
  - A first try with `--force-ocr --deskew` made a 24.9 MB file and was discarded. The page images of the filed copy
    were not changed (size ratio 1.05).
  - PDF page 1 = cover, 2-10 = printed pp.1-18 in pairs (text, equations, list of symbols, two references), 11-16 =
    Figs. 1-23.
- **How it was read:** every page was read on a 100 dpi render. The three tables, the constant-mass list and every
  figure annotation (Figs. 5-9 and 16-20) were read again on 250 dpi crops. The text layer was the third witness for the
  text and the tables. No table digit differs between the three. The OCR errors were '1919' for 1975, '1960' for 1980,
  '1,470,000 ld' and a stray space in '12, 200'.
- **Independent arithmetic checks, all ours:**
  - The constant masses sum to 60,000 lb. Adding the 14,500 lb life support of the 440-day example gives the printed
    74,500 lb. The life support formula (3500 lb plus 25 lb/day) gives 14,500 lb at 440 days.
  - The dV of Figs. 5-9, in EMOS (Earth's mean orbital speed, taken as 29.78 km/s = 97,700 ft/s, our value), gives
    5,471, 8,891, 12,310, 12,213 and 11,236 ft/s. The printed table has 5,470, 8,900, 12,300, 12,200 and 11,200.
    All five agree to 0.1 percent. The table is a second witness for the figure labels and the other way round.
  - The Julian dates of Figs. 16-20 give total trips of 430, 450, 450, 460 and 470 days. These match the powered points of
    Fig. 12 (read off the plot). 2443290 + 625 = 2443915, as printed for the 1978 Venus case.
  - Fig. 17 has two labels cut by a trajectory line ('244 12340' and '244 2247'). They read as 2442340 and 2442247
    from the clean partner labels and the 450-day total.

## 0. Verdict

**A one-shot Earth-Mars-Earth powered-flyby study. It contains no periodic or repeating trajectory and no cycler.**
- It takes the ballistic Mars flyby (Earth escape, Mars flyby, Earth entry) of the opposition years 1971, 1973, 1975, 1978
  and 1980 and adds an impulse inside the Mars sphere of influence. The ballistic flyby fixes the passage distance
  once the dates are fixed. The impulse removes that restriction.
- Result (text, p.1 and p.16): in the unfavorable years 1978 and 1980 the powered flyby needs less than half the mass in
  Earth orbit of the ballistic one, and the trip is about three months shorter. In 1971 and 1973 the mass saving is small
  and the time saving is the same. The price is a high Earth entry speed.
- A string search of the OCR text for "period", "repeat", "recur", "cycl" and "synodic" finds only "most favorable period" for
  surface observation. The paper makes no statement that a trajectory repeats.

## 1. What the paper does

- **Method (pp.2-6, Eqs. 1-9, Figs. 1-4).** A planetocentric two-body, two-dimensional problem inside the Mars sphere
  of influence. Two hyperbolas with the same point of transfer at radius R are matched. Independent variables are the
  two excess speeds V_inf1 and V_inf2 (in EMOS), the bend angle K, the passage distance r_p and the true anomaly of the
  transfer point. Two modes: Mode I (Fig. 1, covered by Gobetz, Ref. 1) and Mode II (Fig. 2, new here).
  - Passage distance fixed at r_p = 1.1 Mars radii, just above an assumed atmosphere edge at 1.09.
  - The transfer dV is the vector difference of the two hyperbola speeds at the transfer point (Eq. 8). The program
    runs backwards (indices swapped) to get the cases where the governing passage is on the second hyperbola.
  - Fig. 3 shows the dV curve against the true anomaly of the transfer point for both modes. Fig. 4 shows the true
    anomaly pairs for fixed r_p = 1.1.
- **Transportation system (p.7).** Constant masses: entry vehicle 15,000 lb, command module 10,000, power 1,000, radiation
  shelter 20,000, service module 5,000, added structure 4,000, corrective propulsion 5,000. Life support 3500 lb plus 25
  lb/day for a crew of three. LH2/LOX at 430 s. The mass in Earth orbit is the optimisation measure. Mass is built
  backwards from the Earth return.
- **Results (pp.8-10).** The hyperbolic transfer dV and Earth escape dV by opposition year, in ft/s (table, p.9-10):

  | Year | Earth escape dV | Hyperbolic transfer dV |
  |---|---|---|
  | 1971 | 12,000 | 5,470 |
  | 1973 | 12,800 | 8,900 |
  | 1975 | 13,200 | 12,300 |
  | 1978 | 13,100 | 12,200 |
  | 1980 | 12,800 | 11,200 |

  The two dV are nearly equal from 1975, so the hyperbolic stage could serve as an abort stage (Conclusion 5). The mass in
  Earth orbit is nearly constant from 1975 to 1980 (Conclusion 6). 1978 optimum powered: 750,000 lb; 1978 ballistic:
  1.8 million lb.
- **The Earth-entry-speed problem and four cures (pp.10-15).** Powered flybys have entry speeds above 65,000 ft/s
  (Fig. 13, powered about 52, 65, 72, 73 and 73 thousand ft/s, plot reading).
  1. **Earth-entry retro:** up to 9000 ft/s, applied to the 15,000 lb entry vehicle.
  2. **Perihelion retro** on the Mars-Earth leg (Figs. 14-15).
  3. **Nonoptimum powered flybys** (table p.14): 1975, 1978, 1980 with totals of 510, 530 and 520 days, Earth departure
     19,100, 21,000 and 17,150 ft/s, dV at Mars 6550, 6850 and 6650 ft/s, entry 64,500, 64,900 and 64,800 ft/s.
  4. **A Venus swingby on the Mars-Earth leg (pp.11-13, Figs. 16-21).** Venus is treated as a coplanar perturber met by
     chance on the way back, not as a leg of the design.
     - 1971: spacecraft and Venus are close at the crossing; Earth entry is already below the limit, so a Venus flyby
       at five or more radii is suggested only for a close look.
     - 1975: a unique Venus passage at 1.8 Venus radii gives an Earth rendezvous, 23 days longer, entry 57,200 ft/s.
     - 1978: Venus is 15 days behind at the second crossing. Revising the Mars-Earth leg needs a Venus passage below
       1 radius for most cases, so 1.2 radii is imposed. All combinations with entry below 65,000 ft/s need 0.5 radii or
       less. The case computed (depart JD 2443290, Mars +270 d, Venus +200 d at 1.4 radii, Earth JD 2443915, 625 d)
       needs 2,935,000 lb in Earth orbit, against 750,000 lb optimum powered and 1.8 million lb ballistic. Judged
       not justified.
     - 1973 and 1980: Venus is too far away.
     - A class of Venus-Earth transfers of about 400 days with entry below the limit exists but the totals exceed 700 days.
  - **Combination:** a perihelion retro of 6250 ft/s plus an Earth retro of 8000 ft/s is the best mix for the 1980 case with a
    50,000 ft/s limit (1,470,000 lb, against 2,511,000 lb for Earth retro alone and 1,561,000 lb for perihelion retro
    alone, p.15).
- **Conclusions 1-12 (pp.16-17).** Besides the ones above: the required maneuver happens after pericenter and after
  the best surface lighting (3-4); the combination retro is the usual optimum (10); and the Venus swingby should be used
  "only if Venus is normally near the intersection of the spacecraft with Venus' orbit" (12).
- **References:** only two. Gobetz, UAC Research Laboratories Report A-110058-5 (November 1962), and Willis, NASA Lewis,
  AIAA Paper 65-66. Neither is in the wanted list as a row (see sec. 4).

## 2. Figures (described, not digitised except the printed labels)

- **Figs. 1-2 (PDF p.11):** the two hyperbola geometries, mode I and mode II. **Fig. 3:** the dV curve against the true anomaly of the transfer point, with the omitted region near the ordinate. **Fig. 4:** a family of parabola-like curves of psi_2 against psi_1 with the psi_1 minimum limit.
- **Figs. 5-9 (PDF pp.12-13):** the Mars hyperbola of each optimum powered flyby, with the Mars disc shaded for the
  illuminated part. Printed labels are in the tables file (V_inf1 0.150 to 0.190, V_inf2 0.212 to 0.307, dV 0.056 to 0.126
  EMOS). The transfer point moves further from Mars as the year goes from 1971 to 1980 (text p.8, Fig. 10).
- **Fig. 10:** dV against true anomaly of the transfer point for five years. **Fig. 11:** mass in Earth orbit as bars for the
  optimum, Earth-retro, perihelion-retro and nonoptimum modes, with the ballistic bar for each year. **Fig. 12:** mass
  against total trip time. **Fig. 13:** entry speed bars. **Figs. 14-15, 22-23:** retro trade curves.
- **Figs. 16-20 (PDF pp.14-15):** plan views with the orbits of Mercury, Venus, Earth and Mars. The heavy line is the
  powered flyby trajectory. The dashed line is the ballistic one. Stars mark crossings of the Venus or Mercury orbit,
  with the planet position on the same date beside them. The Julian dates printed are in the tables file. **Fig. 21:**
  Venus longitude against date as a sawtooth, with the five spacecraft intercepts.

## 3. Relevance to this project

- **Not a cycler, not periodic, no collision with any catalogue row.** Every trajectory is a single outbound-and-return
  journey with one Mars encounter. Dates are tied to one opposition each. No period, no repetition, no resonance.
- **R1(b) (the Hollister lineage):** Titus is a precursor of the Venus-swingby round trips, not a member. His Venus swingby
  is a post-hoc perturbation of the Mars-Earth leg of the optimum powered flyby, tested to reduce entry speed, and it
  fails the mass test in 1978. Hollister & Prussing 1965 (held) and Sohn 1964 (held) design the Venus leg in. It is the
  closest published statement from this period that a Mars flyby leg can be perturbed by an incidental Venus encounter.
- **For Young 1966 (held):** Young lists Titus as his Ref. 2 and gives no number that I could find from him. Young's
  reference list is the only place. Titus gives no value that Young reuses, and the digest of Young (2026-10-06) records
  that Titus was only cited as lineage. The numbers in this file are therefore new to the corpus.
- **Numbers useful as a sanity set (V0 level, ballistic comparisons):**
  - The powered-flyby total times of 430 to 470 days with Mars at 150 to 210 days and the dates in Figs. 16-20 are a
    published set of 1971-1980 Earth-Mars-Earth flybys for any later check of our Mars-flyby legs. They are
    plot-and-label values; the paper gives no state vectors.
  - The dates also show the near-repeat of Earth-Mars geometry: Earth departures at JD 2441090, 2441890, 2442670, 2443430
    and 2444190 are 800, 780, 760 and 760 days apart (ours). This is about the 2.1 to 2.2 year Earth-Mars synodic cadence,
    not a statement in the paper.
- **Risk of misuse:** the unit "EMOS" is Earth's mean orbital speed, which the paper never states. The Mars dV
  (0.056 to 0.126 EMOS) and the table in ft/s are consistent only if EMOS = 29.78 km/s. The paper uses imperial
  units; the project values in km/s are ours (see the tables file).

## 4. Citation mining (checked against CORPUS_INDEX.md and the wanted list, 2026-10-11)

- **Gobetz, F. W., UAC Research Laboratories Report A-110058-5, November 1962 (Ref. 1):** not in CORPUS_INDEX and not on
  the wanted list. Gobetz AIAA 63-422 is cited as Young 1966 Ref. 5 and is also not a row. This is a United Aircraft
  company report, so a public copy is unlikely. Low priority; **not added as a row** because Titus itself says only that
  it covers mode I of the hyperbolic transfer geometry.
- **Willis, E. A., Jr., AIAA Paper 65-66, "New Class of Optimal Interplanetary Trajectories with Specified Trip Time"
  (Ref. 2):** not in CORPUS_INDEX or the wanted list. Optimal trip-time-constrained trajectories; no cycler content
  expected. Low priority; **not added**.
- Nothing else: the paper has two references.
- Held neighbours (CORPUS_INDEX): Hollister & Prussing 1965, Sohn 1964, Hollister 1963 thesis, Hollister 1966
  (JSR), Young 1966. Deerwester AIAA 65-89 stays wanted.

*Filed as `cyclers_pdf/papers/titus-1965-powered-flybys-of-mars-aiaa-65-515-doi-10.2514-6.1965-515.pdf`. Table transcription: `data/sources/titus-1965-tables.yaml`.*
