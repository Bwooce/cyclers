# Digest: Miele and Wang 1999, "Optimal trajectories and mirror properties for round-trip Mars missions" (Acta Astronautica 45:655-668) (#960, owner-supplied)

A. Miele and T. Wang (Aero-Astronautics Group, Rice University), Acta Astronautica 45(11):655-668 (December 1999; copyright 2000),
PII S0094-5765(99)00117-4, DOI 10.1016/S0094-5765(99)00117-4. Received 5 September 1998; presented as paper IAF-98-A.1.07, 49th International
Astronautical Congress, Melbourne, 1998. Title, authors, volume, issue, pages and date agree with Crossref (checked 2026-10-11). Not a wanted-list
row: the owner supplied it, so there is no row to remove. No other Miele paper is held (CORPUS_INDEX searched; Miele appears only as a citation).

- **Filed as** `cyclers_pdf/papers/miele-wang-1999-optimal-trajectories-mirror-properties-round-trip-mars-acta-astronautica-45-655-doi-10.1016-S0094-5765-99-00117-4.pdf`.
  Table transcription: `data/sources/miele-wang-1999-tables.yaml`.
  - Publisher PDF as received, md5 b8838006eb5aca17fe6622d9d74ef4a0, 14 pages (Letter), text layer present, no OCR. PDF page n = journal page 654 + n.
    The text layer garbles minus signs and some symbols.
- **How it was read:** the whole text layer; every table and Fig. 5 on renders (250 dpi crops of the four tables, 70 dpi pages). Three witnesses for Tables 1-4:
  the text layer, the 250 dpi crops read by eye, and tesseract on the crops. All cells agree.

## 0. Verdict

**A one-shot Earth-LEO to Mars-LMO and back mission, optimised for total characteristic velocity in a restricted four-body model with circular coplanar
planets. The optimum is a Hohmann-class transfer each way with a 454-day wait: 11.30 km/s and 970 days. Nothing repeats. The "mirror property" is a
time-reversal and reflection symmetry that the paper shows graphically and does not prove.**
- Result: each leg 5.65 km/s and about 258 days; stay in low Mars orbit 454.3 days; total 11.30 km/s, 970.0 days (2.66 years).
- The mirror property (Sec. 8): the return trajectory is the outgoing trajectory after a rotation by 238 deg about the Sun, a reflection about the
  Earth-Sun-Mars axis, and an inversion of the sense of motion. It holds in interplanetary space and in the near-Earth and near-Mars frames.
- The authors call it "a qualitative property", verified only as a graphic superposition at the scale of Figs. 5-7 (p.668).

## 1. Content

- **Model (Sec. 2, Eqs. 1-9):** Sun fixed; Earth and Mars on circular coplanar orbits (Eqs. 2, 4); the spacecraft attracted by Sun, Earth and Mars
  along the whole path, in Sun-centred polar variables (r, phi, V, gamma). No patched conics. Two-impulse legs, impulses tangential at the circular
  terminal orbits (hypotheses A1-A7).
- **Boundary conditions and transformations (Secs. 3-4, Eqs. 10-30):** four epochs t0 (leave LEO), t1 (arrive LMO), t2 (leave LMO), t3 (arrive LEO), with
  the Sun-to-Earth and Sun-to-Mars frame transformations.
- **Optimisation (Sec. 5):** minimise the sum of the four impulses. Nine parameters (four impulses, three times, the Mars-Earth phase difference at
  departure, the spacecraft phase angle at departure), six constraints, three degrees of freedom. Sequential gradient-restoration algorithm (SGRA)
  with a variable-stepsize integrator that shrinks the step with the total gravitational acceleration.
- **Data (Sec. 6):** the constants of the tables file (mu_Sun 1.327E11, mu_Earth 3.986E05, mu_Mars 4.283E04 km^3/s^2; r_E 1.496E08, r_M 2.279E08 km;
  LEO radius 6841 km, LMO radius 3597 km).
- **Results (Sec. 7, Tables 1-4):** outgoing 3.55 + 2.10 = 5.65 km/s in 257.9 days (Table 1); stay of 454.3 days with 0.00 km/s (Table 2); return
  2.10 + 3.55 = 5.65 km/s in 257.8 days (Table 3); summary (Table 4) with mass ratios 4.5, 1.0, 4.5 and 20.2 (Isp 450 s, structural factor 0.10).
  The Mars-Earth phase difference is +43.9 deg at departure, -75.1 deg at arrival, +75.0 deg at departure from Mars and -43.9 deg at return.
- **Mirror properties (Sec. 8):** (i) equal impulses, (ii) equal flight times, (iii) phase differences equal and opposite, (iv) the wait follows from (iii):
  the phase difference must move from -75.1 to +75.0 deg, which cannot be done backwards, so it changes by 360 - 150 = 210 deg at 0.462 deg/day, or 454 days.
  Mirror construction in Figs. 5-7, four steps each. The three-body precedent is Miele 1960 (Ref. 10) for the Earth-Moon problem.

## 2. Table values (all in the tables file)

| Leg | Time (days) | Impulses (km/s) | Mars-Earth phase difference (deg) |
|---|---|---|---|
| Outgoing LEO to LMO | 257.9 | 3.55 + 2.10 = 5.65 | +43.9 to -75.1 |
| Stay in LMO | 454.3 | 0.00 | -75.1 to +75.0 |
| Return LMO to LEO | 257.8 | 2.10 + 3.55 = 5.65 | +75.0 to -43.9 |
| Total | 970.0 | 11.30 | |

## 3. Our checks (all computed by us; nothing taken from the paper's code)

- **The optimum is a Hohmann transfer.** With the printed mu and radii, a circular Hohmann transfer takes 258.8 days and needs v_inf of 2.943 km/s at
  Earth and 2.648 km/s at Mars. Impulses from the printed LEO and LMO radii are 3.556 and 2.101 km/s, total 5.657 km/s. The paper has 3.55, 2.10, 5.65
  and 257.9 days. The difference is 0.4 percent in time and 0.1 percent in impulse. The paper never uses the words Hohmann, conjunction or opposition.
  Its swept angles (179.0 and 179.1 deg) are the Hohmann 180 deg.
- **Phase bookkeeping closes.** The phase difference falls at the Earth-Mars rate 0.4614 deg/day: +43.9 - 0.4614 x 257.9 = -75.1 deg; -75.1 - 0.4614 x 454.3 =
  -284.7 deg, equal to +75.3 deg modulo 360 (printed +75.0, rate rounding); 75.0 - 0.4614 x 257.8 = -44.0 deg (printed -43.9).
- **Mid-mission geometry.** At 485 days (half of 970) the phase difference is -179.9 deg. So Earth and Mars are on opposite sides of the Sun at the middle of
  the stay: Sun conjunction. The three-step construction (rotate by 238 deg, reflect, invert) is then a single reflection about the Sun-Earth-Mars line
  through that instant, combined with time reversal about it. This single-axis reading is ours; the paper states the three steps. It also explains
  why the paper's own remark "Earth and Mars aligned with the Sun, within one degree" is needed.
- **Sums and units:** 257.9 + 454.3 + 257.8 = 970.0; 5.65 x 2 = 11.30; 4.5 x 4.5 = 20.25 (printed 20.2); the synodic period from the printed rates is
  780.2 days, so the mission is 1.24 synodic periods. The mass ratios are not reproduced (the paper's mass-ratio formula is not stated).

## 4. Relevance to this project

- **Mirror property, free-return and cycler symmetry.**
  - It is the Sun-Earth-Mars version of the image theorem that Miele gave for the Earth-Moon problem in 1960 and that Schwaniger 1963 uses (digest
    2026-10-07, "Two kinds of free return"). In the planar circular restricted problem the symmetry is a reflection about the x axis with time reversal.
    Here it is the same pair of operations about the instant at which the two planets sit on one line through the Sun.
  - **Why it holds (ours):** circular coplanar planetary motion with a fixed Sun is invariant under reflection about any Sun-planet line combined with time
    reversal about the instant the planets lie on that line. A minimum-impulse trajectory is then expected to be symmetric about that instant, or to come in
    a mirror pair. The paper's observation (equal impulses, equal times, opposite phases) is what this predicts. It is not a proof; the paper offers none.
  - **Link to #1047 sec. 9:** that section's rows A, B and D are asymmetric Earth-Moon periodic orbits and the 7-3b and 7-3c orbits are one asymmetric orbit and its
    mirror. The symmetry there maps an orbit to a distinct mirror orbit, and a symmetric orbit maps to itself. Miele and Wang's round trip is the second case
    (the trajectory maps to itself with legs swapped), but it is not periodic, so it is no closure test. They are analogous, not the same object.
  - **Hollister lineage:** a periodic Earth-Mars orbit that is symmetric about an Earth-Mars-Sun alignment is the repeating cousin of this single round trip. This
    paper has no flyby, no Venus and no resonance, so it adds nothing for the swingby cyclers (R1(b)). I did not re-read the Hollister digests for this note.
- **Baseline for the Earth-Mars rows:** the 11.30 km/s, 970 days and 454-day wait are a published circular-coplanar conjunction-class value that any
  Earth-Mars row's round-trip cost can be set against. The impulses include LEO escape (3.55 km/s) and LMO capture (2.10 km/s), so a fair comparison with a
  cycler row is in v_inf (our Hohmann values 2.94 and 2.65 km/s) or with the same terminal orbits. It is not an ephemeris result and gives no calendar dates.
- **Does any trajectory repeat?** No. One round trip, no periodicity, no resonance. The 454-day wait is the part a cycler avoids by letting the vehicle
  keep orbiting. The geometry repeats at the synodic period (780 days by our calculation), which the paper does not mention.
- **Catalogue implication:** none. No row is a cycler. At most a baseline citation for the conjunction-class round trip.

## 5. Citation mining (against CORPUS_INDEX.md and the wanted list, 2026-10-11; proposals only, no row edited)

- **[10] Miele, A. (1960), "Theorem of image trajectories in the Earth-Moon space", Astronautica Acta 4(5):225-232.** Neither held nor wanted. The
  Schwaniger digest cites it as a Boeing report (Flight Sciences Laboratory Report No. 21, January 1960); this is probably the journal form of the same theorem (not confirmed). **Proposed row (medium):**
  it is the symmetry theorem behind the project's perpendicular-crossing correctors and the Earth-Moon free-return families.
- **[9] Miele and Wang (1997), "Optimal transfers from an Earth orbit to a Mars orbit", IAF-97-A.6.05, 48th IAC Turin (Acta Astronautica to appear).**
  Neither held nor wanted. The outgoing-only precursor of this paper. Low priority; propose a row only if the baseline becomes load-bearing.
- **[2] Niehoff (1988), "Pathways to Mars: new trajectory opportunities", in The NASA Mars Conference, Univelt, pp. 381-401:** already wanted (row 9, which
  records the 1988 Univelt volume as the library form).
- **[8] Casalino, Colasurdo and Pastrone (1998), "Optimization procedure for preliminary design of opposition-class Mars missions", JGCD 21(1):134-140.**
  Neither held nor wanted. Opposition-class round trips; the usual counterpart of the conjunction-class case here. Low priority.
- **[6] Wagner and Munk 1996, [5] Lee and Sidney 1996 (AAS Spaceflight Mechanics volume), [4] Striepe and Desai 1996 (JAS 44:207), [7] Wercinski 1996 (JSR 33:381),
  [3] Braun et al. 1995:** neither held nor wanted; mission designs and entry analysis, not cycler or symmetry content. Not proposed.
- **[16] Miele, Huang and Heideman 1969 (JOTA 4:213) and [17] Miele, Wang and Basapur 1986 (Acta Astronautica 13:491):** the SGRA method papers; held-neither,
  not proposed (method only).
- **[11] Roy (1988), Orbital Motion; [12] Miele (1962), Flight Mechanics; [13] Gurzadyan (1996):** books, not proposed. **[14] and [15]:** Miele and Wang aeroassisted
  transfer papers, not relevant.
- Held neighbours (CORPUS_INDEX): Schwaniger 1963, Rall and Hollister 1971, Hollister and Prussing 1965, Sohn 1964, Young 1966 (Venus-swingby Mars round trips). I did not check whether any of them cites this paper.

*Filed as `cyclers_pdf/papers/miele-wang-1999-optimal-trajectories-mirror-properties-round-trip-mars-acta-astronautica-45-655-doi-10.1016-S0094-5765-99-00117-4.pdf`. Table transcription: `data/sources/miele-wang-1999-tables.yaml`.*
