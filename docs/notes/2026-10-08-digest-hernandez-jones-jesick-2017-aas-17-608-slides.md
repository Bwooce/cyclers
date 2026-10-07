# Digest: Hernandez, Jones & Jesick 2017, AAS 17-608 slide deck (supplement to the held paper)

S. Hernandez, D. R. Jones and M. Jesick, "Families of Io-Europa-Ganymede Triple Cyclers", slides for
paper **AAS 17-608** (printed bottom left of every slide). Jet Propulsion Laboratory, California
Institute of Technology. "(c) 2017 California Institute of Technology. Government sponsorship
acknowledged." The slides do not print the meeting, so the venue comes from the paper's own reference
list (AAS/AIAA Astrodynamics Specialist Conference, Columbia River Gorge, WA, August 2017); I did not
check that on Crossref. A Crossref search for the deck found no DOI. The paper is held (see
`identity-hernandez2017-aas17-608.md`).

- Source file: upload `258c49ec-CL17-4128.pdf`, 15 pp, md5 8301ae0708a9bc67b4e822dbc5ec0d8a, 720 x 540 pt slides
  (4:3), JPL template. The pages are 11 numbered slides plus animation builds.
- Proposed corpus filename: `hernandez-jones-jesick-2017-one-class-io-europa-ganymede-triple-cyclers-AAS-17-608-slides.pdf`.
- How I read it. The text layer is usable but carries ligature damage ("IntroducBon", "ﬂight ;me") and
  the slide builds overprint each other in the text. I did not OCR it, because the page images are clean
  and a copy adds nothing. I read pages 2 and 4 to 15 on 100-dpi page images; pages 1 and 3 (title slide
  and an animation build of slide 2) I read from the text layer. Every number below was read on an
  image and matched against the held paper's text layer (grep of the same digits).

## 0. Verdict

The deck adds **almost nothing** to the held paper. Every table and number in it appears in the paper.
Two small items are new: the itinerary examples by synodic count (slide 2 build) and the explicit
"3 parameters to vary / n parameters for an n-sequence cycle" statement (slide 7). Keep it as a
supplement. It needs no catalogue or code change.

A useful side fact: the deck follows the **final** paper, not the preprint. It carries the
sequence-count formula s(n) = 3^(n-1) - 2^n + 1 and the 3/4/5/6-encounter table (2, 12, 50, 180),
which the preprint (item b) does not have. Its title, however, is the preprint title ("Families of ...").

## 1. Slide map and what is new

| Slide (PDF page) | Content | New vs paper? |
|---|---|---|
| 1 (p2) Introduction | Cycler definition, crew and cargo Earth-Mars use, inter-moon robotic tours, "triple cyclers" | No (paper intro) |
| 2 (p3-p5) IEG system | 1:2:4 resonance, T_syn = 7.05 days, inertial shift 5.2 degrees per synodic period; cycle and cycler definitions | Numbers are in the paper (7.05 d, 5.2 deg). **Itinerary examples are new**, below |
| 3 (p6) Tisserand graph | Io, Europa, Ganymede contours; search box r_p 1 to 6 R_J, r_a 15 to 25 R_J; EGIGE example with v_inf 10 (Europa), 6 (Ganymede), 9 (Io), 7 km/s | No (paper Fig. 1 and text) |
| 4 (p7) Strategy | s(n) formula, table 2/12/50/180 for 3/4/5/6 encounters, four-step strategy (phase-space reduction, sequence plus time of flight, patched conic plus Lambert, Monte Carlo) | No (paper Eq. 2 and Table 2) |
| 5 (p8) Conic search | n_rev T_i = n_syn T_syn; r_pmin = R_Jup, r_pmax = a_Io, r_amin = a_Gan, r_amax = infinity; two argument-of-periapsis options (omega1 outbound, omega2 inbound); revs table | No (paper Table 1 and p.4 equations) |
| 6 (p9) Conic search | "3 independent variables looped"; sequence and approximate time of flight per cycle | No |
| 7 (p10-p11) Initial guess to Lambert | Zero-sphere-of-influence patched conic; input n_syn, sequence, tof; output delta-V, r_p; EIGE chain t_dep, tof(E-I), tof(I-G), tof(G-E); fixed T_cycle = sum of tofs | **"3 parameters to vary" and "n parameters to optimize for n sequence cycle" are slide wording** (paper says the same in prose) |
| 8 (p12) One synodic period EGIEIE | Same four panels and table as paper Fig. 3 and Table 3 | No |
| 9 (p13) Four synodic period EGGIE | Same as paper Fig. 4 and Table 4 | No |
| 10 (p14) High fidelity | Ideal delta-V = 0 m/s; real ephemeris 10 repeat cycles delta-V = 30 m/s; high fidelity ballistic for 2 cycles | No (paper p.10) |
| 11 (p15) Conclusion | Same as paper conclusions, plus future work | No |

## 2. Numbers read on the slides (all present in the paper)

- T_syn = 7.05 days; shift 5.2 degrees; 1:2:4 resonance (slide 2).
- Table 1 revs (slide 5): 1 syn: 1:1, 1:2; 2 syn: 2:1, 2:3, 2:5; 3 syn: 3:1, 3:2, 3:4, 3:5, 3:7; 4 syn: 4:1, 4:3, 4:5, 4:7, 4:9.
- Slide 8 EGIEIE, departing 03-Oct-2020: Europa 12.27 km/s, 9,260 km, 1.8 m/s; Ganymede 2.03 d, 7.29, 69,990, 0.1; Io 0.85, 15.77, 3,176, 4.4; Europa 0.66, 12.15, 2,259, 0.6; Io 2.87, 15.86, 494, 0.0; Europa 0.65; total 7.06 d and 6.90 m/s. Identical to paper Table 3.
- Slide 9 EGGIE: Europa 9.12 km/s, 1,444 km, 0.00; Ganymede 1.59 d, 7.07, 2,155, 0.60; Ganymede 8.60, 7.07, 6,263, 0.00; Io 7.34, 8.38, 653, 0.10; Europa 10.69; total 28.22 d and 0.70 m/s. Identical to paper Table 4. (Slide title "Four Synodic Period: EGGIE Cycle".)
- Slide 10: EIGE one synodic period, ideal delta-V 0 m/s; real ephemeris 30 m/s after 10 repeat cycles. Same as paper p.10.

## 3. The only new content

Slide 2 (page 5 build) lists example itineraries by synodic period count:
- 1 synodic: EGIE, EIGIE
- 2 synodic: EGIE, EIGIE, GEIIEIG, IEIGIEI
- 3 synodic: EGIE, EIGIE, GEIIEIG, IEIGIEI, EIGIGEIGE

These strings do not appear in the paper (grep of GEIIEIG, IEIGIEI, EIGIGEIGE found none). Treat them as illustrative
only. The slides give no trajectory for them.

## 4. Citation mining

The deck has no reference slide. Its content cites nothing beyond the held paper.
Held check: `ls papers | grep -i hernandez` finds only the held AAS 17-608 paper; the index row for it is in
`CORPUS_INDEX.md` (hernandez-jones-jesick-2017-one-class-...-AAS-17-608.pdf, digest
2026-06-26-digest-hernandez-2017-ieg-triple-cyclers-aas-17-608.md). The deck itself is not held and not
listed (no filename with "slides" for this paper). It is not in the wanted list
`2026-10-05-960-wanted-papers.md` (grep for 17-608 and Hernandez found nothing).

## 5. Proposals (not applied)

- File the deck as a supplement under the proposed filename and add one index row pointing to the 2026-06-26 digest. No new digest is needed beyond this note.
- No catalogue or code change.

*Filed as `cyclers_pdf/papers/hernandez-jones-jesick-2017-one-class-io-europa-ganymede-triple-cyclers-AAS-17-608-slides-jpl-cl-17-4128.pdf`.*
