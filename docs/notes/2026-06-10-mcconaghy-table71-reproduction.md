# McConaghy 2004 Table 7.1 reproduction on DE440 — second-source check on S1L1 (#94)

**Date:** 2026-06-10
**Task:** #94 follow-on — independently validate S1L1 against a SECOND published
source: reproduce McConaghy 2004 (Purdue PhD, UMI 3166673) Table 7.1 — the
24-encounter, 33-year DE405 itinerary of S1L1 outbound cycler vehicle 1 — per-leg on
DE440, with REBOUND/IAS15 n-body checks. Data source: the mining note
`docs/notes/2026-06-10-mcconaghy-2004-dissertation-mining.md` §2.2 (the dissertation
PDF itself was NOT available to this task). Tool:
`scripts/reproduce_mcconaghy_table71.py`.

**Headline verdict: PARTIAL-CONFIRMED, with the discrepant tail ADJUDICATED.**

- **Encounters 1–19 (2005-08-13 → 2031-03-12, six full E→M→E→E cycles) REPRODUCE**
  on DE440: |V∞| residual ≤ 0.195 km/s (mean 0.047 over 36 leg-endpoints), ballistic
  flyby |V∞| continuity ≤ 0.19 km/s, implied periapsis altitude within ≤ 1,800 km of
  the printed closest approach at every interior encounter, all Mars encounters
  inside the 3-SOI band on the independent integrator, every E→E g leg sub-Mars.
- **Encounters 20–24 as transcribed DO NOT reproduce** (residual grows 0.15 → 1.56
  km/s; implied closest approach inconsistent at E22/M23). This is a clean negative
  on the transcribed values — and a third source (Russell App-C #83, independently
  DE440-confirmed in #167) **overlaps these exact dates and sides with OUR emerged
  values**, not the transcribed ones (§4). The tail rows need re-verification against
  the dissertation PDF; nothing was tuned to hide this.

**NO catalogue writeback.** Level implications in §6.

---

## 1. Method (the #167 recipe, adapted to McConaghy's data shape)

McConaghy prints per encounter: calendar date (day precision), |V∞| km/s,
closest-approach **altitude** km (Mars-23's 1,770 km < Mars's 3,396 km radius pins
the altitude convention), and arriving-leg TOF days. Unlike Russell App-C there are
no v∞ *vectors*, so per-leg seeding is Lambert between the published encounter dates
on DE440 (epochs at 12:00 TDB on the printed date; McConaghy used DE405 — the Δ is
reported, not tuned away):

1. **Leg construction:** Lambert (multi-rev to N=2; the S1L1 legs sweep > 360°,
   McConaghy Table 6.2: θ = 657.97°/522.29°) between each consecutive encounter
   pair. The discrete (n_revs, branch) choice is selected by best match to the
   published endpoint |V∞| — a topology selection among well-separated conics
   (runner-up gaps 1.2–16.3 km/s, reported per leg). Within the chosen branch
   nothing is fit: the emerged V∞ at both ends is evidence (published = golden).
   The recovered structure is perfectly regular: E→M = 0-rev, M→E = 1-rev/high,
   E→E = 1-rev/low, identically in all 7–8 cycles.
2. **Binding constraints, all in the residual:** encounter dates are inputs
   (imposed), so the date constraint is checked via the integer-day TOF
   cross-check (dates vs printed leg TOF agree to ±1 d on all 23 legs — pure
   day-rounding of the printed dates). What emerges and is checked: |V∞| at both
   leg ends (vs published), flyby |V∞_in| vs |V∞_out| continuity (ballistic
   requirement), implied periapsis altitude vs printed closest approach (from the
   bend angle, sin(δ/2) = 1/(1 + r_p v²/μ) — never used in selection), and n-body
   arrival miss. "Closed on a subset" is the danger signal; all of these are
   tabulated below for every encounter.
3. **Independent n-body:** every leg's departure state propagated by REBOUND/IAS15
   over DE440 — Sun-only (the patched-conic cruise model) for all 23 legs, plus
   real DE440 Mars as a continuous perturber for the 8 E→M legs. Mars band =
   3 × r_SOI(Mars) ≈ 0.0116 AU ≈ 1.73e6 km, the SAME band as #165/#167, never
   loosened.

Caveat on circularity, stated plainly: a Lambert leg lands on the target planet by
construction, so the Sun-only IAS15 miss (~0 km) verifies integrator/solver
consistency, not geometry. The non-circular content is (a) the emerged V∞ vs
published at 46 leg-endpoints, (b) flyby continuity, (c) implied vs printed closest
approach, (d) the Mars-perturbed misses (real flyby gravity acting continuously),
(e) the g-leg sub-Mars structure.

## 2. Per-leg results (emerged vs published V∞)

```text
leg    kind  TOF_d  Nrev/branch  vinf_dep  pub_dep  d_dep | vinf_arr  pub_arr  d_arr
 1->2   E->M    198  0/single       4.002    4.01  -0.008 |   3.022    3.02  +0.002
 2->3   M->E    833  1/high         3.039    3.02  +0.019 |   6.900    6.89  +0.010
 3->4   E->E    542  1/low          7.029    6.89  +0.139 |   7.038    6.90  +0.138
 4->5   E->M    185  0/single       6.929    6.90  +0.029 |   4.322    4.31  +0.012
 5->6   M->E    810  1/high         4.269    4.31  -0.041 |   6.350    6.42  -0.070
 6->7   E->E    539  1/low          6.295    6.42  -0.125 |   6.309    6.43  -0.121
 7->8   E->M    139  0/single       6.411    6.43  -0.019 |   7.097    7.14  -0.043
 8->9   M->E    890  1/high         7.139    7.14  -0.001 |   4.001    4.01  -0.009
 9->10  E->E    529  1/low          3.815    4.01  -0.195 |   3.837    4.03  -0.193
10->11  E->M    116  0/single       3.998    4.03  -0.032 |   6.434    6.47  -0.036
11->12  M->E    934  1/high         6.462    6.47  -0.008 |   4.607    4.61  -0.003
12->13  E->E    532  1/low          4.573    4.61  -0.037 |   4.552    4.59  -0.038
13->14  E->M    223  0/single       4.598    4.59  +0.008 |   2.762    2.77  -0.008
14->15  M->E    793  1/high         2.735    2.77  -0.035 |   6.996    7.08  -0.084
15->16  E->E    542  1/low          7.036    7.08  -0.044 |   7.041    7.09  -0.049
16->17  E->M    170  0/single       7.126    7.09  +0.036 |   5.248    5.27  -0.022
17->18  M->E    830  1/high         5.263    5.27  -0.007 |   5.774    5.80  -0.026
18->19  E->E    537  1/low          5.817    5.80  +0.017 |   5.812    5.80  +0.012
19->20  E->M    125  0/single       5.820    5.80  +0.020 |   7.705    7.85  -0.145   <- tail
20->21  M->E    915  1/high         7.693    7.85  -0.157 |   3.771    4.21  -0.439   <- tail
21->22  E->E    529  1/low          3.824    4.21  -0.386 |   3.805    4.20  -0.395   <- tail
22->23  E->M    137  0/single       3.784    4.20  -0.416 |   4.746    5.87  -1.124   <- tail
23->24  M->E    906  1/high         4.749    5.87  -1.121 |   5.670    7.23  -1.560   <- tail
```

Encounters 1–19 (36 endpoints): max |Δ| 0.195, mean 0.047 km/s — consistent with the
±0.5 d date quantisation (measured sensitivity ~0.03–0.07 km/s per day, §4). The
slightly larger residuals concentrate on the 1-rev E→E legs, whose multi-rev
geometry is the most date-sensitive. Encounters 20–24: max |Δ| 1.56 km/s — far
beyond any day-rounding explanation.

## 3. Per-encounter flyby checks (continuity + implied closest approach)

```text
enc  body  date        pub_vinf  vinf_in  vinf_out  |in-out|  bend_deg  implied_alt_km  pub_CA_km
  2  M     2006-02-27    3.02     3.022    3.039    0.018     39.62           5700       4816
  3  E     2008-06-09    6.89     6.900    7.029    0.130     27.96          19415      20130
  4  E     2009-12-03    6.90     7.038    6.929    0.109     20.99          30321      31110
  5  M     2010-06-06    4.31     4.322    4.269    0.053     11.61          17236      17710
  6  E     2012-08-24    6.42     6.350    6.295    0.055     26.80          26679      26490
  7  E     2014-02-14    6.43     6.309    6.411    0.102     19.17          42943      41520
  8  M     2014-07-03    7.14     7.097    7.139    0.042      6.02          11862      12190
  9  E     2016-12-09    4.01     4.001    3.815    0.186     49.78          29525      27730
 10  E     2018-05-22    4.03     3.837    3.998    0.161     58.11          21130      19920
 11  M     2018-09-15    6.47     6.434    6.462    0.028      7.22          11928      11580
 12  E     2021-04-06    4.61     4.607    4.573    0.033     45.87          23251      22990
 13  E     2022-09-20    4.59     4.552    4.598    0.045     56.14          15049      14780
 14  M     2023-05-01    2.77     2.762    2.735    0.027     40.21           7426       7601
 15  E     2025-07-02    7.08     6.996    7.036    0.040     24.59          23558      23860
 16  E     2026-12-26    7.09     7.041    7.126    0.084     18.11          36144      35120
 17  M     2027-06-14    5.27     5.248    5.263    0.015      9.55          13689      13840
 18  E     2029-09-21    5.80     5.774    5.817    0.043     30.73          26548      26850
 19  E     2031-03-12    5.80     5.812    5.820    0.008     23.70          39234      37520
 20  M     2031-07-15    7.85     7.705    7.693    0.012      5.68          10469       8802   <- tail
 21  E     2034-01-15    4.21     3.771    3.824    0.053     58.61          22459      24870   <- tail
 22  E     2035-06-28    4.20     3.805    3.784    0.021     79.05           9442       2756   <- tail
 23  M     2035-11-12    5.87     4.746    4.749    0.003     11.43          13782       1770   <- tail
```

Every interior flyby is ballistically consistent (max |in−out| 0.186, mean 0.058
km/s — the date-noise level), **including the tail**: the date-anchored DE440
geometry forms a perfectly good ballistic cycler throughout; in the tail it is
simply not the one the transcribed V∞/CA values describe. In the head region the
implied periapsis altitude tracks the printed closest approach to ≤ 1,800 km
(shallow-bend flybys make this quantity very sensitive to the ~0.1 km/s V∞ noise);
at E22/M23 the printed deep flybys (2,756 / 1,770 km) are irreconcilable with the
printed dates (implied 9,442 / 13,782 km).

**E→E (g) free-return legs — corrected-topology signature (#166/#167):** all seven
have aphelion 1.175–1.327 AU (< Mars's 1.52) and closest approach to real DE440 Mars
0.66–0.96 AU. Same structure as the #167 Russell reproduction (aphelion < 1.4,
clearance 0.67–1.05 AU).

## 4. Adjudicating the tail: Russell App-C #83 overlaps these dates

Diagnostics first (script + `/tmp` probe, both reproducible):

- **No alternative Lambert branch matches the tail.** All (N ≤ 2) × branch
  candidates enumerated; the non-selected ones are 6.7–16.3 km/s worse. The printed
  tail values correspond to NO conic through the printed dates.
- **Date sensitivity is ~0.03–0.07 km/s per day** (measured on leg 22→23 over ±4 d).
  The 1.12 km/s gap at M23 would need ~15-day date errors, yet the printed dates and
  printed TOFs are mutually consistent to ±1 d across the whole table. The tail rows
  are therefore **internally inconsistent as transcribed** (dates+TOF vs V∞+CA).
- **Third source.** Russell 2004 App-C #83 — the per-leg real-ephemeris S1L1 block
  #167 independently confirmed on DE440 to 4-decimal v∞ — turns out to overlap
  Table 7.1's calendar dates from 2027 on (the two publications describe near-identical
  members of the same outbound-vehicle-1 trajectory; Russell's itinerary starts
  2025-06-26, 5.7 d before McConaghy's E15):

| enc | body | date | Russell node Δt (d) | Russell v∞ | McC printed | our emerged |
|---|---|---|---|---|---|---|
| 15 | E | 2025-07-02 | −5.7 | 5.802 | 7.08 | 7.036 |
| 17 | M | 2027-06-14 | −0.9 | 5.248 | 5.27 | 5.263 |
| 18 | E | 2029-09-21 | −0.3 | 5.767 | 5.80 | 5.817 |
| 19 | E | 2031-03-12 | −0.5 | 5.764 | 5.80 | 5.820 |
| 20 | M | 2031-07-15 | +0.0 | 7.693 | **7.85** | 7.693 |
| 21 | E | 2034-01-15 | +0.0 | 3.771 | **4.21** | 3.824 |
| 22 | E | 2035-06-28 | −0.2 | 3.752 | **4.20** | 3.784 |
| 23 | M | 2035-11-12 | +1.2 | 4.657 | **5.87** | 4.749 |
| 24 | E | 2038-05-06 | +1.9 | 5.512 | **7.23** | (terminal) |

At every coincident node from 2027 on, Russell's printed (and #167-confirmed) v∞
agrees with OUR emerged values to ≤ 0.16 km/s, and from encounter 20 the transcribed
McConaghy values are the outlier (up to 1.7 km/s off at E24). Encounter 15 is the one
overlap point where Russell differs from both (his itinerary's start node, 5.7 d
earlier — the family's launch-window DOF; McConaghy's printed 7.08 matches OUR 7.04
there, and 7.09 at E16 matches the dissertation's own §7.3 family bound, smallest max
taxi-Earth V∞ = 7.0864).

**Most likely explanation:** transcription drift in the mining note's last table
rows (rendered-PDF table reading; possibly bleed from the adjacent Table 7.5 — whose
boxed Mars maximum 7.70 km/s at Mars-20 matches our emerged 7.693–7.705 exactly), or
less likely a misprint in the dissertation itself. **Unresolvable from the mining
note alone — the PDF is not available to this task.** Follow-up: re-verify Table 7.1
rows 20–24 (V∞ and closest-approach columns) against the PDF
(`cyclers_pdf` private repo) in a task with source access.

## 5. Independent n-body (REBOUND/IAS15 over DE440)

All 23 legs converge (energy drift 0.0). Sun-only misses are ~0 km (Lambert BVP +
GOLDEN-GATE-1-pinned integrator — consistency, not geometry; see §1 caveat).
Mars-perturbed (real DE440 Mars as continuous perturber, the flyby gravity acting
along the whole approach), for the 8 E→M legs:

| Mars arrival | M2 | M5 | M8 | M11 | M14 | M17 | M20 | M23 |
|---|---|---|---|---|---|---|---|---|
| miss (km) | 44,942 | 24,041 | 9,442 | 11,170 | 52,958 | 16,825 | 8,006 | 19,749 |

All 8 inside the 3-SOI band (1.73e6 km) by ≥ 33× — the same order as #167's
Mars-perturbed run (6,600–40,900 km). The band was NOT loosened.

## 6. What this supports (level implications — held for main session)

- **This is a second-source check on the #167 V3** for S1L1 /
  `russell-ch4-4.991gG2` (#94). It now stands confirmed on DE440 from TWO published
  per-member itineraries: Russell App-C #83 (2025–2055, v∞-vector seeding, #167) and
  McConaghy Table 7.1 encounters 1–19 (2005–2031, date-anchored Lambert, this work).
  The McConaghy reproduction covers 2005–2025 — a window Russell's block does NOT —
  so the validated span of the S1L1 vehicle now extends back to the 2005-08-13
  launch.
- **Independence caveat, stated plainly:** Russell and McConaghy are collaborators
  and their published itineraries coincide (≤ 1.2 d, ≤ 0.1 km/s) from 2027 onward —
  over the overlap window these are two optimisation runs of essentially the same
  physical trajectory, not fully independent solutions. The genuinely additional
  validation content here is (a) the 2005–2025 segment, (b) the closest-approach
  cross-check (Russell App-C prints no CA), (c) reproduction from a data shape with
  no v∞ vectors (dates-only Lambert), which exercises a different reconstruction
  path.
- **No V4+ claim.** Nothing here is a continuous one-seed propagation or an
  in-environment (GMAT-class) validation; per-leg re-anchoring at published nodes is
  exactly the V3-class evidence #167 established. The #169 continuous-chain measure
  is unaffected.
- **Proposed (NOT performed) writeback:** add McConaghy 2004 Table 7.1 (encounters
  1–19 reproduced on DE440, this note) as a corroborating source on the S1L1 row,
  alongside the #167 evidence chain; record the rows-20–24 transcription flag in the
  mining note's ledger after PDF re-verification.

## 7. Honesty ledger

- EXPECTED side = the transcribed Table 7.1 values only (golden, from the mining
  note, which traces every value to the printed page); nothing our code computed was
  used as a target. The emerged V∞ / continuity / implied CA / n-body misses are
  EVIDENCE.
- The discrete Lambert (n_revs, branch) selection per leg uses the published V∞ —
  it picks among 1–5 conics separated by 1.2–16.3 km/s, then everything continuous
  emerges. Selection transparency: the pattern is uniform across cycles
  (E→M 0-rev, M→E 1-rev/high, E→E 1-rev/low) and per-leg runner-up gaps are printed.
- The tail failure is REPORTED, adjudicated with a third sourced reference, and
  flagged for PDF re-verification — not tuned away, not absorbed by loosening any
  band, and the clean-negative finding (printed tail values inconsistent with
  printed tail dates) stands on the internal evidence alone (Lambert exhaustion +
  measured date sensitivity) even without the Russell overlap.
- 12:00-TDB day-precision epochs and DE405→DE440 are acknowledged conventions; their
  measured effect (~0.05 km/s mean) bounds the head-region residuals.
- NO catalogue writeback; no test was added to the gated suites (this is a
  diagnostic/reproduction script, the #167 gate tests are untouched).

## 8. Reproduce

```sh
uv run python scripts/reproduce_mcconaghy_table71.py            # full (needs rebound)
uv run python scripts/reproduce_mcconaghy_table71.py --no-nbody # Lambert stage only
```
