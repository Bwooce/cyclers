# Digest: Causey 1964, "Characteristic Features of Some Periodic Orbits in the Restricted Three Body Problem" (#960; consumers #997, #1000, #948)

W. E. Causey (Aero-Astrodynamics Laboratory, NASA George C. Marshall Space Flight Center, Huntsville AL),
NASA TM X-53108, 17 August 1964. NTRS 19640019815, accession N64-29729. No DOI: a Crossref search (2026-10-08)
finds no match. Approved by R. F. Hoelker (Chief, Astrodynamics & Guidance Theory Division) and E. D. Geissler
(p.59). Distribution includes Hoelker, Schwaniger, Arenstorf and Davidson (p.60). (pp. 59-60 read in the
text layer and on a 50 dpi sheet only.)
- Supplied file `a4ac84c7-19640019815.pdf`, 65 PDF pages (covers, front matter, report pp. 1-60), md5
  8e5b0c2ad86328eee23c740c29ed70d4. NTRS scan, 1-bit CCITT at 300 dpi, Acrobat Capture text layer. The old text
  layer is letter-spaced, garbles fractions ("ratio $", "Q", "F" for 1/2, 3/5, 2/5) and some numbers
  ("l3 2 , O O O" for 132,000), and is useless on the sideways table and the graphs.
- OCR copy (the file to be filed), made with `ocrmypdf --force-ocr -l eng`:
  `cyclers_pdf/papers/causey-1964-characteristic-features-periodic-orbits-restricted-three-body-problem-nasa-tm-x-53108-ntrs-19640019815.pdf`,
  65 pp., md5 28fc3a2cc57a9da4833b648ae3757037, PDF/A-2b, 10.2 MB, made with `-O 3` added (section 7).
  Render check at 100 dpi on PDF pp. 1, 30 and 63: same pixel size as the original, mean grey difference 7.0,
  2.1 and 8.5 of 255, ink fraction equal within 0.001. The page images survive.
- Proposed filename: `cyclers_pdf/papers/causey-1964-characteristic-features-periodic-orbits-restricted-three-body-problem-nasa-tm-x-53108-ntrs-19640019815.pdf`.
- How I read it:
  - Text (report pp. 1-7, 59-60): old text layer for structure; every number and statement quoted below was read
    on 170 dpi renders of PDF pp. 7-12.
  - Fig. 51 (p.58, the only table): 400 dpi render, rotated, two OCR witnesses (`witness-comparison.tsv`).
  - Graphs: Figs. 1-50 seen on 50 dpi contact sheets for type and axes; Fig. 10 (p.17), Fig. 44 (p.51) and
    Fig. 49 (p.56) at 150 dpi for the unit legends and axis calibration. No graph value is transcribed.
- Not on the wanted list (grep for "causey": no hit). Not in the corpus. The held Arenstorf 1963 (AIAA J) and
  Hitzl 1977 digests already mention "Davidson and Causey" computations (Hitzl 1977 Figs. 1-4 redraw two
  Arenstorf-Davidson-Causey Earth-Moon orbits at mu = 1/82.30).

## 0. Verdict

**An early (1964) atlas of Earth-Moon symmetric periodic orbit families at Arenstorf's mass ratio mu = 1/81.45,
with a table of orbits that pass both a 100 n.mi. injection perigee and the Moon. (Newton 1959 and Schwaniger
1963 are earlier sources of Earth-Moon periodic orbits; Newton's perigees are 13,400 km and up, and Schwaniger has
one retrograde orbit.) The only numeric table, Fig. 51, has 7
orbits. All 7 reproduce in the planar CR3BP: period and time to the Earth to the printed digits, second lunar
approach within 1-76 km, perigee within 22 km, perigee speed consistent with the perigee. Five of the 7 pass the
lunar test of the Hoelker-Winston digest (a pass inside the lunar Hill radius every period); their perigees are
161-230 km altitude, at or just under that test's 180 km (6,558 km) floor. Row 1 is in the same 2:1 family as the catalogue's `vaquero-21-*` rows and as Newton
1959's "direct" type-1/2 orbits: continuations in mu and x0 join Newton's last direct row to Causey row 1, and
Causey row 1 to the catalogue's `vaquero-21-c198` state, to 10 digits (section 4).**

- **Model (pp. 2-3).** Planar circular restricted problem. Earth/Moon mass ratio 80.45 (so mu = 1/81.45 =
  0.01227747, the "Arenstorf" mu); 385,000 km for physical units; primaries' period 2 pi. Figure legends:
  1 T.U. of time = 105.13 hr, 1 T.U. of velocity = 1017.42 m/sec. These two units disagree by 0.02 per cent
  (385,000 km / 105.13 hr = 1017.26 m/s).
- **Method (pp. 2-3).** Start perpendicular to the Earth-Moon line on the far side of the Moon; vary the
  starting speed until a second perpendicular crossing occurs; then vary the start radius to trace a family.
  Orbits retrograde at the Earth are excluded (p.3). Every orbit is symmetric about the x axis.
- **Classification (pp. 3-4)** follows Arenstorf (ref. 3, the 1963 IAC paper) and M. C. Davidson (private
  communication): **ratio m/k** (Moon makes m revolutions while the probe makes k), **order n** (ratio nm/nk is
  "ratio m/k order n"), **class A/B** (second perpendicular crossing behind / in front of the Moon). Order 1 =
  Arenstorf's periodic orbits of the second kind; higher orders are "orbits of the restricted three body problem
  proper" (p.4). Families studied: ratios 1/2 (orders 1, 2A, 2B, 3, 4A), 2/3 (orders 1, 2A, 2B), 2/5 (order 1)
  and 3/5 (orders 1, 1*, 2B).
- **Data.** Figs. 5-43: closest approach to Earth, second lunar approach, period, per cent time on the inner leg
  and the perigee spread Delta R, each against the start radius. Figs. 44-50: rotating-frame starting speed (TU)
  against start position. Fig. 51: 7 orbits with about 100 n.mi. injection perigees. Only Fig. 51 is numeric.
- **Reproduction (section 3).** For each Fig. 51 row I solve the starting speed at the printed start radius.
  The root is selected by the printed period (which crossing) and the printed perigee (which root). The checks
  that remain independent all pass: the period to 0.1 d, the second lunar approach, the time to the Earth, and
  the perigee speed (to the 10-22 km perigee offsets of rows 2, 3 and 7). For row 1 the rotating-frame starting
  speed also agrees with the Fig. 49 curve; the other rows were not compared with the graphs.
- **A printing error: the "Velocity (Space-Fixed) (m/sec)" column is the geocentric inertial perigee speed
  divided by 10.** Our integration gives 10,888.7 m/s for printed 1088.9 (row 1), 10,892.7 for 1089.3 (row 4),
  10,881.7 for 1088.2 (row 5), 10,836.3 for 1083.6 (row 6). It is not the starting speed (our rotating-frame
  ydot0 runs from 1.06 to 2.46 TU over the rows). The "Time On Inner Leg (%)" column was not checked.
- **Cycler verdict (section 5):** rows 2, 3, 4, 6 and 7 combine a 100 n.mi. perigee with a lunar pass inside the
  lunar Hill radius in every period. Rows 1 and 5 have no lunar pass closer than 89,080 and 182,551 km.
  All perigees (6,539-6,599 km) sit at the low edge of the Vaquero LEO-GEO band (6,558-42,164 km); rows 1 and 4
  (6,552 and 6,539 km) are just below 6,558 km.
- **Catalogue (PROPOSALS only):**
  1. For the 2:1 interior family held as `vaquero-21-c198/c246/c247/c266` (catalogue lines 56588-56900): cite
     Newton 1959 (direct type 1/2, Table 1) as the earliest printed members found so far, and Causey 1964 Fig. 51
     row 1 (ratio 1/2 order 1) as the first printed member with a LEO perigee. Causey gives a sourced start radius
     (89,080 km), period (27.2 d), perigee (6,552 km) and Earth-to-Moon time (155 h) at mu = 1/81.45.
     The starting speed is ours, not printed, so no `state_nd` can be called SOURCED. The continuations in
     section 4 are the evidence that these are one family; it also reproduces the c198 DERIVED state to 10 digits
     from an independent starting point.
  2. The Arenstorf row (`arenstorf-em-figure8-1963`, lines 9077-9230) has null period, Jacobi and state. Causey
     does not fill them: the row does not say which member it means, and Causey's Fig. 1 orbit (ratio 1/2
     order 1, p.8) is a chain of three loops along the x axis in the rotating frame (a small loop at the Moon, a
     loop round the Earth, a large loop beyond the Earth; two self-crossings on the axis), with a 27-d period. It is
     not the textbook 17.07-TU "Arenstorf orbit". That textbook orbit (the standard ODE test case, the orbit behind
     the catalogue's Cook 2020 source; IC typed from memory, confirmed only by its closure to 8e-11, citation not
     checked) has its closest Earth approach at 178,000 km at the same mu (`arenstorf_hairer_compare_output.txt`); it is not in Causey's Fig. 51. Optionally,
     add Causey 1964 to that row's notes as the MSFC atlas of Arenstorf-classified families at the Arenstorf mu.
  3. Rows 2-4, 6 and 7 are Earth-Moon orbits with a 100 n.mi. perigee and a lunar pass of 1,900-19,500 km every
     period that are not in the catalogue. Before any row: a literature-novelty check, a continuation test against
     the Casoliva/Vaquero/Ross rows, and a stability check (none printed). See section 5.

## 1. Content

- Sec. I (p.1): orbits "offer repeated approaches to both earth and moon and could be used in instrumented
  exploration of earth-moon space for meteoroid concentration and radiation belts". Periods 1 to 4 months.
- Sec. II.A (p.2): Poincare's first sort (near one mass); Arenstorf's second kind (ref. 2, TN D-1859, near
  rotating Kepler ellipses); and orbits that exist only for mu > 0 but degenerate into second-kind orbits as the
  lunar disturbance goes to zero. The report presents the last two kinds.
- Sec. II.C (pp. 3-4): Kepler's third law gives the smallest m/k for an orbit that encloses both bodies,
  (1/2)^(3/2) = .354 (with a >= 1/2). (Our check: 0.3536.)
- Sec. II.D (p.4): a lifetime of a year or more is wanted; the real Moon orbit is eccentric, so a station-keeping
  budget is needed, except for orbits whose period is a whole number of lunar months, which close in both the
  rotating and the inertial frame. p.5: close Earth approaches come in pairs; Fig. 3 has four, with a 1,240 km
  spread in altitude.
- Sec. III (pp. 5-7), family ranges read on the page:
  - 1/2 order 1: start radius from the lunar surface to 89,400 km; above that, Earth impact, then retrograde orbits.
  - 1/2 order 2 class B: 3,075-57,000 km (Earth collision outside); contains an orbit of period 4 pi (2 lunar
    months), so it persists for an eccentric lunar orbit.
  - 1/2 order 2 class A: Earth collision at 1,928 km start, where the second crossing is behind the Moon at
    132,000 km; the two crossings coincide at about 15,000 km; beyond that, duplicates.
  - 1/2 order 3: Earth collision below 2,310 km; the second crossing is behind the Earth; a sister family with the
    second lunar approach on the ascending leg exists (Fig. 16) but has no data.
  - 1/2 order 4 class A: closest Earth approach 11,895 km at zero lunar altitude.
  - 2/3 order 1: from the lunar surface to 183,000 km.
  - 2/3 order 2 class A: closest Earth approach 93,610 km at the smallest start radius (1,738 km); duplicates
    beyond 8,900 km.
  - 2/3 order 2 class B: an orbit of period 8 pi at 51,000 km; dual solutions near 1,994 km and 110,000 km.
  - 2/5 order 1: 7,800-19,800 km; Earth collision outside.
  - 3/5 order 1: 3,187-74,026 km, two solutions per radius (the faster marked *), differing by 6.4-38.4 m/s.
  - 3/5 order 2 class B: closest Earth approach at a start radius of 2,802 km.
- Sec. IV (p.7): higher orders and other ratios "will be described in a later paper". Fig. 51 summarises orbits
  with injection altitudes of about 100 n.mi.

## 2. Fig. 51 transcription

`causey-1964-tables.yaml` (7 rows, 10 columns, numbers as printed, conventions block). Counted by
`checks/witness_counts.py` over the 49 value cells: 30 identical in all three readings; 7 differ only by
punctuation; 12 have a digit error in one OCR witness (5/9, 7/1, dropped, doubled or run-on digits). In all 12 the
other OCR witness agrees with my image reading, and, where they apply, the family ranges of pp. 5-6 rule the OCR
reading out. Ratio, order and class were read on the image only. Row 7's order digit "3"
looks hand-written over the type; order 3 fits the text (2,310 km limit) and the Fig. 14 period axis.

## 3. Integration check (`checks/causey_fig51_check.py`, output `checks/causey_fig51_check_output.txt`)

Planar CR3BP, mu = 1/81.45, DOP853 rtol = atol = 1e-12 (scan rtol 1e-10). Start x0 = (1 - mu) + r/385,000,
y = 0, xdot = 0. ydot0 is scanned on [-2.6, -0.3] (step 0.003) and [+0.3, +2.6] (step 0.01) for rows 1, 4, 5
and 6; rows 2, 3 and 7 use narrow negative windows (script header), because their near-Moon starts made the full
scan too slow and a coarse step lost the row-2 bracket to an Earth-impact neighbour. I keep the roots
with xdot = 0 at the y = 0 crossing nearest the printed half period, then the root nearest the printed perigee.
The Jacobi constant below is C = x^2 + y^2 + 2(1-mu)/r1 + 2 mu/r2 - v^2. Closure is the largest state error
after one full period.

| row | ratio, order, class | start r (km) | our ydot0 (TU) | P ours / printed (d) | perigee ours / printed (km) | 2nd lunar ours / printed (km) | first Earth pass ours / printed (h) | perigee speed ours (m/s) / printed x 10 | C | closure |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1/2, 1 | 89,080 | -1.1244880 | 27.236 / 27.2 | 6,552 / 6,552 | none / - | 155 / 155 | 10,888.7 / 10,889 | 1.93211 | 3e-11 |
| 2 | 1/2, 2, A | 1,931 | -2.4564039 | 53.260 / 53.3 | 6,608 / 6,586 | 130,319 / 130,395 | 70 / 70 | 10,847.5 / 10,867 | 1.81292 | 1e-8 |
| 3 | 1/2, 2, B | 3,081 | -2.0111348 | 54.343 / 54.3 | 6,546 / 6,560 | 3,340 / 3,341 | 76 / 76 | 10,891.6 / 10,879 | 1.97494 | 7e-9 |
| 4 | 1/2, 2, B | 56,713 | -1.0601479 | 55.160 / 55.2 | 6,540 / 6,539 | 4,352 / 4,352 | 136 / 136 | 10,892.7 / 10,893 | 2.05288 | 7e-11 |
| 5 | 2/3, 1 | 182,551 | -1.3574971 | 55.028 / 55.0 | 6,574 / 6,574 | none / - | 208 / 208 | 10,881.7 / 10,882 | 1.68613 | 1e-10 |
| 6 | 2/5, 1 | 19,484 | -1.1194322 | 54.411 / 54.4 | 6,600 / 6,599 | none / - | 107 / 107 | 10,836.3 / 10,836 | 2.19049 | 2e-10 |
| 7 | 1/2, 3, A | 2,317 | -2.2688610 | 80.446 / 80.4 | 6,587 / 6,577 | 91,485 / 91,507 | 72 / 72 | 10,861.4 / 10,869 | 1.88354 | 2e-8 |

- Ours are DERIVED values (the starting speed is not printed). Periods use 1 TU = 105.13 hr. "First Earth pass" is
  the time from the start point to the first local minimum of the Earth distance; "none" means no lunar minimum
  inside 100,000 km other than the start point.
- **Every printed period and first-pass time is reproduced to the printed digits.** Second lunar approaches agree
  to 0-76 km (row 2: 130,319 vs 130,395; row 7: 91,485 vs 91,507; rows 3-4 to 1 km). Perigees agree to 0-22 km.
- Rows 2, 3 and 7 start 1,931-3,081 km from the Moon's centre and graze the Earth; the perigee there changes by
  about 200 km per 0.001 TU of starting speed. Their perigees are 10-22 km off. In each of these rows the printed
  speed differs from ours by the amount the perigee offset explains (v roughly proportional to r^-1/2:
  22 km gives about 18 m/s, against 19.5 m/s seen for row 2). So the printed perigee and speed are consistent with
  each other; the 1964 orbit is a few kilometres from our root. This is an integration-accuracy difference in
  1964 or a rounding of the start radius, not a misprint.
- Rows 2, 3 and 7 were run with a targeted scan (see the script header). The output file holds, in order: the
  full scans of rows 1, 4, 5, 6, the targeted runs (row 2 coarse, which missed, then rows 3 and 7), and
  the fine row-2 run. The full scans of rows 1, 4, 5 and 6 found
  5-26 symmetric roots near the printed half period. In each, the selected root is the only one with a perigee
  within 1,000 km of the printed value; the next roots have perigees of 2,000-47,000 km.
- **Time column.** "Time From Earth to Moon (hr)" equals the time from the FIRST Earth passage to the start point
  behind the Moon (by symmetry, Moon to Earth). Where an orbit has two different perigees, the first passage is
  not always the lowest one. Rows 4 and 6 show this: their first passages are at 7,037 km (136 h) and 7,629 km
  (107 h), while the printed close approach is the lower, later perigee.
- **Fig. 49 cross-check (row 1).** Our ydot0 = -1.1245 TU at x0 = 1.2191. Fig. 49's "Ratio 1/2, Order 1" curve
  reads about 1.13 at 1.22 TU, so the graph plots |ydot0| in the rotating frame, and the sign is negative (motion
  toward -y behind the Moon).
- The Jacobi constant is conserved to 1e-10 or better in every row.

## 4. Relation to held Earth-Moon digests and catalogue rows

- **Newton 1959** (held digest): "type 1/2" orbits at mu = 1/82.45, sidereal period pi about E, one pass at M per
  period, perigees 13,400 km and up. **Same family, by continuation** (`checks/newton_to_causey_row1.py`, output
  `_output.txt`): from Newton's last "direct" row (rho_M0 = 0.212129, ydot0 = -1.049, T = 6.192; re-solved to
  -1.04935900, P = 6.192170) I stepped mu to 1/81.45 and then rho to Causey's 89,080 km offset; it lands on
  ydot0 = -1.1244879632, P = 6.21775871, C = 1.93211238, which is our Causey row-1 root to all printed digits.
  Newton's "retrograde about E" rows are retrograde at the Earth, which Causey excludes (p.3); not continued. So Newton 1959 holds the earliest printed members of
  this family, and Causey extends it to a 100 n.mi. perigee (6,552 km at an 89,080 km start).
- **Hoelker & Winston 1968** (held digest): mu = 1/80, every orbit starts 1.5 units from Earth, labelled by the
  Kepler mean motion n*. No HW orbit has both a close Earth pass and a close lunar pass. Causey (same MSFC group:
  Hoelker approved this report) is the companion that does. By the definitions, HW's n* = 2 should match Causey's
  ratio 1/2 (two probe revolutions per lunar revolution); I did not compare orbits.
- **Schwaniger 1963** (held digest): one cislunar periodic orbit, retrograde (counter-rotation) at the Earth,
  177 km perigee altitude, lunar pass about 2,200 km, period 26 d. Causey excludes orbits retrograde at the Earth
  (p.3), so Schwaniger's orbit is in a family Causey does not treat.
- **Catalogue, Earth-Moon rows (read-only):**
  - `vaquero-21-c198-em-resonant-po-2013` (line 56588): DERIVED state x0 = 1.21394, ydot0 = -1.10114,
    period 6.2111 TU, C = 1.98, perigee 8,409 km, at mu = 0.0121506. Causey row 1 at mu = 1/81.45: x0 = 1.21910,
    ydot0 = -1.12449, period 6.218 TU, C = 1.932, perigee 6,552 km. Same shape and same period.
    **Continuation (`checks/causey_row1_to_vaquero.py`, output `_output.txt`):** I stepped mu from 1/81.45 to
    0.01215058 in 10 steps at a fixed start offset behind the Moon. Then I stepped x0 to the catalogue x0 in 20
    steps, re-solving ydot0 at each step. It lands on ydot0 = -1.1011440751, P = 6.2111417762, C = 1.9800000001.
    The catalogue gives -1.1011440751, 6.2111417763 and 1.98. **So Causey row 1 and `vaquero-21-c198` are the
    same family (one continuation; every step converged with a small ydot0 change; I did not test for folds
    between steps).** Causey 1964 is a
    printed member of this family 49 years before Vaquero 2013, with an Earth-to-Moon time and a 100 n.mi.
    perigee. This also confirms the catalogue's DERIVED state for c198 independently: our corrector, starting
    from a different printed source, gives the same digits.
  - `casoliva-*` (55643-56483): (1,2), (2,1), (3,2) members have no lunar encounter; the (7,3) members pass the
    Moon at 13,000-27,000 km. No (7,3) ratio is in Causey. Causey 1/2 order 2 class B contains a period-4 pi orbit,
    the same resonance as Casoliva's (1,2) label, but Casoliva's (1,2) rows never come near the Moon, so they are
    not Causey's class-B orbits.
  - `ross-rt-*`, `braik-ross-*` (48219-49164): stable prograde (k1, k2) cyclers. None of these rows states a
    LEO perigee (`periapse_km` is null; the (3,3) row states perigee altitudes of 112,400-113,500 km). I did not
    integrate them, so a collision with Fig. 51 rows is not excluded by computation, only by the stated numbers.
  - `arenstorf-em-figure8-1963`: see proposal 2. `genova-aldrin` (bicircular, 3:1) and `wittal-2022` (petal
    family, period 163 d): different models or ratios.

## 5. Cycler verdict per row (test from the Hoelker-Winston digest: perigee in 6,558-42,164 km and a lunar pass
inside the lunar Hill radius)

Lunar Hill radius (mu/3)^(1/3) x 385,000 km = 61,600 km (ours). Lunar passes counted per period, from our
integration (start point included).

| row | lunar passes inside 61,600 km per period | Earth passes per period (km) | verdict |
|---|---|---|---|
| 1 | none (start 89,080 km is the closest) | 2 x 6,552 | not a cycler: no lunar encounter |
| 2 | 1 (1,931 km) | 7,387; 6,608; 6,608; 7,387 | candidate |
| 3 | 2 (3,081 and 3,340 km) | 7,551; 6,546; 6,546; 7,551 | candidate |
| 4 | 1 (4,352 km; the start at 56,713 km is the far side) | 7,037; 6,540; 6,540; 7,037 | candidate |
| 5 | none (start 182,551 km) | 6,574; 7,498; 6,574 | not a cycler: no lunar encounter |
| 6 | 1 (19,484 km) | 7,629; 6,763; 6,600; 6,763; 7,629 | candidate |
| 7 | 1 (2,317 km) | 7,506; 6,587; 9,205; 9,205; 6,587; 7,506 | candidate |

- "candidate" = passes the lunar test. All Earth passes are in the 6,539-9,205 km range, at the very bottom of the
  6,558-42,164 km band: the lowest perigees are 161-230 km altitude (with 6,378 km Earth radius), so the 180 km
  floor (6,558 km) is met by rows 2, 3, 5, 6 and 7 as printed (row 3 only just, at 6,560 km; our row-3 root is
  6,546 km) and missed by rows 1 and 4 (printed 6,552 and 6,539 km). These are 100 n.mi. injection orbits by design (caption), so every
  Earth pass is a low pass.
- Periods are 27-80 d (1-3 lunar months). Row 4's family contains a period-4 pi member (p.5) and the 2/3 order
  2B family a period-8 pi member (p.6); those would also close in the inertial frame.
- No stability is printed. Hitzl 1977 (held digest) says close-encounter periodic orbits "in general ... are
  highly unstable"; expect the same here. Not computed.
- The point-mass model is used throughout; rows 2 and 7 start 193 and 579 km above the lunar surface
  (radius 1,738 km).

## 6. Citation mining (p.7; `ls cyclers_pdf/papers | grep -i` and CORPUS_INDEX grep)

| Ref. | Citation | Held? |
|---|---|---|
| 1 | Poincare, Les Methodes Nouvelles de la Mecanique Celeste, Vol. I, p.97 | no (classic; not needed) |
| 2 | Arenstorf, "Periodic Solutions of the Restricted Three Body Problem Representing Analytic Continuation of Keplerian Elliptic Motion", NASA TN D-1859, May 1963 | yes, as the Amer. J. Math. 85:27 version (doi 10.2307/2373181; index line 245) |
| 3 | Arenstorf, "Periodic Trajectories Passing Near Both Masses of the Restricted Problem of Three Bodies", XIVth IAC, Paris, 26 Sep - 1 Oct 1963 | no; wanted-list row 23 (already listed). Causey says his classification follows it. |
| footnote | M. C. Davidson, MSFC Computation Laboratory, private communication | n/a |

Related, found while checking: Arenstorf 1963, "Existence of periodic solutions passing near both masses of the
restricted three-body problem", AIAA J 1:238 (doi 10.2514/3.1516) is held (index line 244); it is probably the
journal note of ref. 3's result, not the IAC paper itself. Crossref also lists Arenstorf, "Trajectories Around Both
Masses in the Restricted Three Body Problem", in "Trajectories of Artificial Celestial Bodies" (Springer 1966),
p.261, doi 10.1007/978-3-642-49326-3_26; not held; it may be the IAC paper's published form. Suggest adding the
DOI to wanted-list row 23 as a candidate.
No "later paper" by Causey on higher orders was found in this search; NTRS was not searched.

## 7. Filing note

A plain `--force-ocr` copy was 41.4 MB (25 times the 1.6 MB scan: ocrmypdf rasterises the 1-bit pages to
400 dpi colour). The filed copy adds `-O 3` (still `--force-ocr`, not `--redo-ocr`) and is 10.2 MB; its render
check is above. Precedent OCR copies in the corpus are 0.9-2.7 MB, so it is still large. The plain copy was
deleted. Logs: `ocr-causey-O3.log`.

*Filed as `cyclers_pdf/papers/causey-1964-characteristic-features-periodic-orbits-restricted-three-body-problem-nasa-tm-x-53108-ntrs-19640019815.pdf`. Check scripts, outputs and notes named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`. Table transcription: `data/sources/causey-1964-tables.yaml`.*
