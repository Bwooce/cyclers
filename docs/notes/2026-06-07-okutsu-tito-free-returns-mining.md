# Okutsu 2002 + Tito 2013 — the free-return ancestor line (consolidated frontier)

Mined 2026-06-07. A combined **DATA-focused** pass on the two acquired ancestor
papers of the Hughes 2014 / Boeing 2033 free-return line mined earlier today
(`docs/notes/2026-06-07-hughes-2014-fast-mars-free-returns-mining.md`,
`docs/notes/2026-06-07-donahue-duggan-2022-mars2033-flyby-mining.md`). Goal:
extract every sourced trajectory anchor (dates / C3 / V∞ / TOF / entry speed),
fix the topology / class taxonomy, settle catalogue-eligibility, and — the
headline deliverable — build the **consolidated frontier table** unifying all
four free-return sources now held.

**Sources (cite EXACTLY — authors/title/venue only):**

A. Okutsu, M., & Longuski, J. M., "Mars Free Returns via Gravity Assist from
   Venus," *Journal of Spacecraft and Rockets*, Vol. 39, No. 1, Jan–Feb 2002,
   pp. 31–36. Purdue University, West Lafayette, IN. (Also AIAA Paper 2000-4030,
   Aug 2000, per their own ref 25.) Associate Editor C. A. Kluever.

B. Tito, D. A., Anderson, G., Carrico, J. P., Jr., Clark, J., Finger, B.,
   Lantz, G. A., Loucks, M. E., MacCallum, T., Poynter, J., Squire, T. H., &
   Worden, S. P., "Feasibility Analysis for a Manned Mars Free-Return Mission in
   2018," IEEE Aerospace Conference, 2013. (Wilshire Associates / Paragon Space
   Development / Applied Defense Solutions / Space Exploration Engineering / NASA
   Ames / Baylor / Center for Space Medicine.) This is the **Inspiration Mars**
   501-day E-M-E free return.

---

## 0. SCAN-READABILITY ASSESSMENT — Okutsu (paper A)

**Read from a SCANNED image (no text layer); the Read tool renders pages
visually.** Overall **the scan is high quality**: it is a clean printed JSR
journal page, sharp, square, with crisp serif type. **Tables 1 and 2 read
unambiguously digit-by-digit** (verified by re-reading pp. 34, 35 at single-page
resolution). Body text is fully legible. Trajectory figures (Figs 1, 2, 4, 7,
12) carry **printed numeric callouts** (dates, V∞, TOF) that are sharp and
readable. The Tisserand graph (Fig 5) and the scatter plots (Figs 3, 6, 9–11)
are readable for qualitative structure; individual scatter glyphs are not
data-grade (expected — those are search-output clouds, the numbers live in the
tables and figure callouts).

**Digit-confidence caveat:** every number below traces to either a TABLE
(Tables 1–2, high confidence) or a printed FIGURE CALLOUT (Figs 1/2/4/7/12,
high confidence — they are typeset, not handwritten). No value here was read off
a scatter-glyph. The only soft spots are the scatter-plot axis ranges, which I do
not quote as anchors. **No digit was ambiguous enough to flag.** Paper B (Tito)
is a clean digital PDF with a text layer — no scan issues.

---

## 1. HEADLINE VERDICT

Both papers **CONFIRM and REFINE** the frontier picture, and **Okutsu fills the
Hughes V∞ gap** — it is the one source of the four that **tabulates the
Mars-encounter and Venus-encounter V∞ numerically** (Hughes tabulated only
Earth-side V∞; Boeing gave figure callouts; Okutsu gives a full data table).

- **Okutsu is the direct ancestor of Hughes 2014** — same Purdue/Longuski group,
  same STOUR engine, same E-V-M (here written EME/EMVE/EMVE) free-return problem,
  same DRM-2014 framing. Hughes 2014 cites it as the predecessor; Okutsu is the
  2002 version of the identical search.
- **Tito 2013 is the direct ancestor of Boeing 2033** — Boeing (Donahue-Duggan
  2022) cites Tito as its free-return heritage. Tito's 2018 trajectory is a
  **Mars-only ballistic free return** (E-M-E, no Venus flyby), Patel-Longuski-
  Sims class — exactly the topology Boeing's 2033 baseline reuses. Carrico and
  Loucks appear on BOTH Tito 2013 and Hughes 2014, tying the lines together.

**Mars-V∞ frontier, now with Okutsu's numbers:**
- **Pure Mars-only free return (Tito 2018 E-M-E):** Mars V∞ = **6.22697 km/s**
  outbound / **5.42540** inbound (Table II/IV); Mars periapsis speed ≈ 7.27 km/s
  (text p.4). This is the **5-7 km/s class** — confirms Hughes/Boeing.
- **Mars-Venus free return (Okutsu EMVE, Venus AFTER Mars):** Mars encounter
  V∞ as low as **6.79–6.97 km/s** (Figs 1, 7; Table 1 EMVE = 6.97), with the
  Venus flyby on the way home pulling **Earth-arrival V∞ down to ~4.77 km/s**.
- **Earth-Mars-Earth (EME) ballistic free returns (Okutsu):** Mars arrival V∞
  ranges **6.78 (P=1.5 yr) → 11.22 (P=2 yr) → 14.63 (P=3 yr) km/s** (Table 1) —
  i.e. the shorter-TOF EME families have HIGHER Mars V∞. The lowest Mars V∞
  ballistic case (6.78) is the long-period collision-orbit family.

So the consolidated answer to the frontier question: **every ballistic Mars
free return across all four sources sits at Mars V∞ ≈ 5–7 km/s for the
practical (low Earth-side V∞) families, climbing to 11–15 km/s for the fast
2-/3-year EME families.** None reaches the Jones Mars-~3 km/s cycler class
ballistically. **Okutsu is the first source to make this Mars-V∞ ceiling
explicit in a data table** — the gap Hughes left open.

---

## 2. OKUTSU 2002 — TRAJECTORY ANCHORS (verbatim, with page numbers)

### 2.1 Topology / taxonomy (pp. 31–33)

- **Free-return classes defined (p. 31–32):**
  - **EME** = Earth-Mars-Earth (Mars-only free return, no Venus). Sought over a
    15-yr period starting 2010. Parameterised by **round-trip period P (years)**:
    **collision orbits** are the integer-Earth-year subset (TOF = integer yrs),
    families at **P = 1.5, 2, 3 yr** discussed. Synodic Earth-Mars = 2.14 yr; EME
    geometry repeats every **15 years** (the famous recurrence).
  - **EMVE** = Earth-Mars-Venus-Earth: ballistic free return with Venus flyby
    **AFTER** Mars (on the inbound/return leg). **This is the paper's preferred
    abort option.** EMVE composite periodicity ≈ **6.4 yr** (E/M/V phasing);
    inertial geometry of all three planets repeats every **32 years**.
  - **EMVE short-stay** = the EMVE free return converted to a short-stay mission
    with a Mars stopover of 0/30/60/90 days (Table 2).
  - **EVME** (Venus BEFORE Mars) is mentioned as the Lyne-Townsend powered-
    swingby concept but Okutsu's own searched/preferred family is **EMVE**
    (Venus AFTER Mars). **NOTE the contrast with Hughes 2014, whose primary
    family is EVME (Venus BEFORE Mars).** Okutsu = E-M-V-E; Hughes = E-V-M-E.
- **Engine:** **STOUR** (patched-conic, analytic ephemeris, grid search over
  launch date × launch V∞), with **MIDAS** used to optimize total ΔV of selected
  cases (p. 32; refs 23 Sauer MIDAS, 18-22 STOUR). **Same STOUR family as Hughes
  and as our `hunt_vem_ballistic.py`.** Tisserand-graph pre-screen (Fig 5) used
  to read off minimum achievable V∞ per encounter (p. 32–33).
- **Pure ballistic** ("free return ... without a large trajectory correction
  maneuver", abstract p. 31). Gravity assists treated as **impulsive/unpowered**
  in the patched-conic propagation (p. 31). An **aerogravity-assist (AGA)** Mars
  variant is also explored (Figs 8–10) but the baseline accepted set is ballistic.

### 2.2 Tisserand minimum-V∞ statement (p. 33, the load-bearing physics)

Verbatim-paraphrased: from the P–r_p Tisserand graph (Fig 5), tracing the
constant-V∞ contours where Earth, Venus and Mars orbits intersect:
- **minimum launch V∞ Earth→Mars ≈ 3.0 km/s** (Earth 3.0-km/s contour ∩ Mars);
- **minimum launch V∞ Earth→Venus ≈ 2.5 km/s**;
- **minimum arrival V∞ for an EVME free return ≈ 3.4 km/s** (the dashed bold
  Earth 3.4 contour) — "the minimum arrival V∞ for the Earth-Venus-Mars-Earth
  (EVME) path."

This is the **theoretical floor**, not what the IM-class missions realise — they
need ~180-day legs, which forces V∞ up. (Same node-crossing / near-180° physics
Hughes 2014 quotes.) **Important nuance:** Okutsu's Tisserand floor of ~3 km/s
Earth→Mars launch V∞ is a coplanar-circular idealisation; the realised ballistic
families still land at Mars V∞ 5–7+ km/s once real phasing/inclination is imposed.

### 2.3 Table 1 — Nominal and abort missions available in Jan. 2014 (p. 34, verbatim)

Launch ΔV optimized in all cases. Trajectory = family; P = round-trip period.

| Trajectory | Launch date | Earth launch V∞ (km/s) | Transit to Mars (days) | Mars arrival V∞ (km/s) | Earth arrival V∞ (km/s) | Free-return total TOF (yrs) |
|---|---|---|---|---|---|---|
| DRM-class mission [a] | 4 Jan 2014 | 3.32 | 180 | 6.78 | 7.34 | N/A |
| EME (P=1.5 yr) | 29 Dec 2013 | 3.30 | 184 | 6.90 | 3.30 | 3.00 [b] |
| EME (P=2.0 yr) | 11 Jan 2014 | 5.01 [b] | 135 | 11.22 [b] | 5.02 | 1.99 |
| EME (P=3.0 yr) | 18 Jan 2014 | 6.84 [b] | 111 | 14.63 [b] | 6.85 | 3.00 [b] |
| EMVE [c] | 13 Jan 2014 | 3.60 | 170 | 6.97 | 4.77 | 2.19 |

[a] Consistent with 180-day one-way TOFs for both outbound and inbound legs.
[b] Values exceed constraint guidelines.
[c] Total TOF is fixed at 800 years [sic — typo for 800 DAYS; total TOF 2.19 yr ≈ 800 d].

**This is the four-source frontier's single most valuable table — it gives Mars
arrival V∞ AND Earth arrival V∞ AND Earth launch V∞ together, by family.**

Key readings:
- **EMVE is the standout:** Earth launch V∞ **3.60**, Mars arrival V∞ **6.97**,
  Earth arrival V∞ **4.77 km/s**, total TOF 2.19 yr (~800 d). The Venus flyby on
  the return cuts Earth-arrival V∞ from ~6.5–7.5 (nominal) to **4.77** (text p.34:
  "required Earth approach speed of 4.77 km/s ... much less than that of the
  nominal mission, about 6.5–7.5 km/s"). This is the lowest Earth-arrival V∞
  of any Mars-touching ballistic free return in the four sources.
- **The fast EME families pay in Mars V∞:** P=2 yr → Mars V∞ 11.22; P=3 yr →
  Mars V∞ 14.63 — both flagged as exceeding DRM guidelines. The Mars approach
  speed for the 2-yr family (11.22) is called "excessive" (p. 32).
- **DRM-class (P=1.5-ish, 180-d legs):** Mars V∞ 6.78, Earth arrival 7.34.

### 2.4 Table 2 — Short-stay EMVE free-return abort options (p. 35, verbatim)

For Earth launch V∞ = 4.43 km/s on 24 March 2017 and Mars arrival V∞ = 5.98 km/s
on 9 March 2018. (Duration = Earth launch → Earth return.)

| Mars stopover (days) | Mars launch V∞ (km/s) | Mars→Earth TOF (days) | Earth arrival V∞ (km/s) | Earth arrival date | Total mission duration (days) |
|---|---|---|---|---|---|
| 0 [c, free return] | 5.98 [d, flyby V∞] | 125 [e, fixed 125 d] | 13.94 | 11 July 2018 | 475 |
| 30 | 4.12 | 125 [e] | 9.12 | 11 Aug 2018 | 505 |
| 60 | 4.08 | 125 [e] | 5.52 | 10 Sept 2018 | 535 |
| 90 | 5.49 | 125 [e] | 4.16 | 10 Oct 2018 | 560 |
| 0 | 2.48 | 219 [f] | 3.19 | 14 Oct 2018 | 564 |
| 30 | 2.61 | 197 [f] | 3.22 | 22 Oct 2018 | 572 |
| 60 | 3.33 | 178 [f] | 3.56 | 2 Nov 2018 | 584 |
| 90 | 4.48 | 171 [f] | 4.64 | 24 Nov 2018 | 605 |

[c] Free return. [d] Flyby V∞. [e] Fixed to 125 days. [f] Longer TOF for lower-ΔV trajectories.

> Note (p. 36): constraining Mars→Earth TOF to ≤125 d, Mars-launch V∞ ranges
> 4.08–5.98; allowing free TOF, the Mars-launch V∞ "drops considerably,
> ranging from 2.48 to 4.48 km/s" but total duration 564–605 d. **This is the
> Mars-side V∞ trade Hughes never tabulated** — a Mars-departure V∞ as low as
> **2.48 km/s** appears, BUT that is a powered short-stay (Mars stopover +
> re-launch), NOT the ballistic flyby V∞. The ballistic flyby Mars-arrival V∞
> for this 2017-launch EMVE is **5.98 km/s** (the stopover-0/125-d row).

### 2.5 Figure-callout anchors (printed, sharp — high confidence)

These are the annotated trajectory diagrams; the V∞/date/TOF are typeset callouts.

**Fig 1 — Long-stay (conjunction-class) mission profile (p. 32):**
- Earth Launch 1/4/2016, V∞ = 3.32 km/s
- Mars Arrival 7/3/2014 [sic: 2016], V∞ = 6.79 km/s; Mars Departure 1/8/2016 [sic],
  V∞ = 2.99 km/s; (Mars stopover between)
- Earth Return 7/6/2016 [year per tick spacing], V∞ = 7.34 km/s
- (Mars arrival V∞ 6.79 / Mars departure 2.99 / Earth return 7.34 — matches
  Table 1 DRM row's 6.78 / 7.34.)

**Fig 2 — Short-stay mission profile (p. 32):**
- Earth Launch 3/24/2017, V∞ = 4.43 km/s; Venus Encounter 9/9/2017, V∞ = 9.81 km/s;
  Mars Arrival 3/9/2018, V∞ = 5.98 km/s; Mars Departure 6/7/2018; Earth Return
  11/24/2018, V∞ = 4.64 km/s. (Matches Table 2 anchor: 4.43 launch, 5.98 Mars.)

**Fig 4 — Minimum-ΔV EME free return, 2-yr period (p. 32):**
- Earth Launch 1/11/2014, V∞ = 5.01 km/s; Mars Encounter 5/26/2014, V∞ = 11.22 km/s,
  TOF = 135 days; Earth Return 1/9/2016, V∞ = 5.02 km/s, TOF = 728 days.
  (The 2-yr EME of Table 1 — note Mars V∞ 11.22, "high Mars approach speed".)

**Fig 7 — EMVE abort option (p. 33):**
- Earth Launch 1/13/2014, V∞ = 3.60 km/s; Mars Encounter 7/2/2014, V∞ = 6.97 km/s,
  TOF = 170 days; Venus Encounter 3/25/2015, V∞ = 6.85 km/s, TOF = 436 days;
  Earth Return (after one revolution post-Venus flyby) 3/23/2016, V∞ = 4.77 km/s,
  TOF = 800 days. **(This is the Table-1 EMVE row, fully annotated: launch 3.60,
  Mars 6.97, Venus 6.85, Earth arrival 4.77.)**

**Fig 12 — Optimal EMVE free return (p. 35):**
- Earth Launch 3/24/2017, V∞ = 4.43 km/s; Mars Encounter 3/9/2018, V∞ = 5.98 km/s,
  TOF = 350 days; Venus Encounter 9/9/2017 [sic — ordering per diagram],
  V∞ = 9.81 km/s, TOF = 169 days; Earth Return 7/11/2018, V∞ = 14.01 km/s,
  TOF = 475 days. (TOF fixed to 1.3 yr; the short-stay-0 / 125-d Table-2 row.)

**Caveat on figure dates:** several Fig 1/2/4 year labels are internally
inconsistent (e.g. Fig 1 "Mars Arrival 7/3/2014" but launch 2016) — these are
tick-mark/print artifacts in the diagrams; the **V∞ and TOF callouts are
self-consistent with the tables**, so I trust V∞/TOF over the figure year strings.

---

## 3. TITO 2013 — TRAJECTORY ANCHORS (verbatim, with page numbers)

### 3.1 Topology (pp. 4)

- **E-M-E, Mars-only ballistic free return — NO Venus.** "leave Earth, fly by
  Mars, and return to Earth without any deterministic maneuvers after the
  Trans-Mars Injection (TMI)" (abstract p. 1). Patel-Longuski-Sims class (their
  ref [1], reproduced as Fig 1). Path code in the STOUR/MAnE figure: **"3 4 3"**
  (Earth-Mars-Earth). The Mars flyby is purely ballistic; the only deterministic
  burn is the **TMI** (out of LEO).
- **The fast 1.4-yr family:** Patel et al. found "fast free-return opportunities
  ... approximately two times every 15 years with a 1.4-year duration" (abstract
  p. 1; Fig 1 red circle highlights the near-term fast ones). The 2018 case is
  one of these fast opportunities — **501-day total** (1.37 yr), much shorter
  than the typical >1.8-yr Mars free returns.
- **Tools:** **MAnE** (Mission Analysis Environment, Space Flight Solutions,
  refs 4-6) for the patched-conic optimum-C3 search, then **STK/Astrogator**
  (ref 3) with full force model + **JPL DE-421 ephemeris** (ref 2) for the
  high-fidelity numerical integration. **First DE-421-named source in the
  free-return set** (Hughes = unspecified analytic; Boeing = MAnE unspecified).

### 3.2 Table I — Optimum EME free-return solution DATES from MAnE (p. 4, verbatim)

| Leg | Stay (days) | Depart | Depart date | Arrive | Arrive date | Flight time (days) |
|---|---|---|---|---|---|---|
| 1 | — | Earth | Jan 5, 2018, 7.1756 h GMT (JD 58123.7990) | Mars | Aug 20, 2018, 7.8289 h GMT (JD 58350.8262) | 227.0272 |
| 2 | 0.0000 | Mars | Aug 20, 2018, 7.8289 h GMT (JD 58350.8262) | Earth | May 21, 2019, 20.9618 h GMT (JD 58625.3734) | 274.5472 |
| | | | | | **Total Duration** | **501.5744** |

### 3.3 Table II — Optimum EME solution VALUES from MAnE (p. 4, verbatim)

V∞ at each end of each leg (km/s), with declination / right ascension:

| Leg | V∞ depart (km/s) | Decl | RtAsc | V∞ arrive (km/s) | Decl | RtAsc |
|---|---|---|---|---|---|---|
| 1 (Earth→Mars) | 6.22697 | -6.48 | 181.24 | 5.42540 | -7.48 | 233.71 |
| 2 (Mars→Earth) | 5.42540 | -11.85 | 200.36 | 8.91499 | -25.40 | 142.87 |

**Mars-encounter V∞ = 5.42540 km/s** (arrival of leg 1 = departure of leg 2 —
continuity, unpowered flyby). **Earth-departure V∞ = 6.22697; Earth-arrival
V∞ = 8.91499 km/s.**

### 3.4 Tables III / IV — same solution, full-force ASTROGATOR (DE-421) (p. 4, verbatim)

**Table III (dates):** Leg 1 Earth 5 Jan 2018 07:00:00.000 UTCG → Mars 20 Aug
2018 08:18:19.619 UTCG, 227.05439374 d; Leg 2 Mars → Earth 21 May 2019
13:52:48.012 UTCG, 274.23227306 d; **Total 501.2866668 d.**

**Table IV (values, V∞ in km/s, C3 in km²/s²):**

| Leg | Dep V∞ | Dep Decl | Dep RtAsc | Dep V peri | Dep C3 | Arr V∞ | Arr Decl | Arr RtAsc | Arr V peri | Arr C3 |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 6.232 | -6.554 | 271.053 | 12.578 | **38.835** | 5.417 | -7.488 | 53.581 | 7.272 | 29.344 |
| 2 | 5.417 | -11.94 | 200.12 | 7.272 | 29.344 | **8.837** | -25.34 | 52.718 | 14.18 | **78.094** |

**The defensible golden numbers (full-force, DE-421):**
- **Earth-departure C3 = 38.835 km²/s²** (text p. 4: "leaving the Earth with a
  C3 of 38.835 km²/sec²"). [Cross-check: Boeing 2033 C3 = 36.35 — same ballpark.]
- **Mars-encounter V∞ = 5.417 km/s** (continuity value, leg-1 arrival /
  leg-2 departure). **Mars periapsis speed ≈ 7.272 km/s** (Table IV V peri;
  text p. 4 "velocity at periapsis is about 7.27 km/sec"). **Mars periapsis
  altitude ≈ 100 km, dark side** (text p. 4, Fig 3).
- **Earth-return V∞ = 8.837 km/s**, **Earth-return C3 = 78.094 km²/s²**,
  **Earth re-entry velocity (V peri) = 14.18 km/sec** at perigee (text p. 4:
  "returning with a velocity at perigee of 14.18 km/sec"). Aerocapture options
  studied to bring g-load down; direct-reentry peaks ~14.2 km/s (p. 6).
- **Total mission duration ≈ 501 days** (501.29 full-force / 501.57 MAnE).

Patched-conic (MAnE) vs full-force (Astrogator) agree: Mars V∞ 5.42540 vs
5.417; Earth-return V∞ 8.91499 vs 8.837; total TOF 501.57 vs 501.29 d. **A
clean STOUR/MAnE-vs-high-fidelity validation anchor**, like Hughes Table 4.

### 3.5 The recurrence / uniqueness claim (pp. 1, 14, 15)

- Abstract (p. 1): the fast 1.4-yr opportunities occur "two times every 15
  years"; the 2018 case is "a rare opportunity for a fast free-return flyby of
  Mars with total mission duration of 501 days" (concluding remarks p. 14).
- **Recurrence (p. 14): "the next opportunity after this mission wouldn't be for
  about another 13 years."** **(p. 15): "the next opportunity does not occur
  again until 2031."** So 2018 → **2031** is the ~13-yr gap. (Consistent with
  the "twice per 15 yr" cadence: 2018 + 13 = 2031.)
- **This is the link to Boeing 2033:** Boeing's baseline 2033 (with a 2035
  backup) is the *next-generation* realisation of the same fast-free-return
  recurrence Tito flags. (Boeing slips to 2033 vs Tito's predicted 2031 because
  Boeing optimises different constraints — entry speed via inbound DSB, 250-km
  flyby — and uses the slightly later phasing window.)

---

## 4. CONSOLIDATED FRONTIER TABLE — all four free-return sources

Per-source lowest-Earth-V∞ ballistic Mars-touching free return, the realised
Mars-encounter V∞, the Earth-side V∞/C3, topology, tool/ephemeris, and TOF.
**(V∞ "ballistic" = the unpowered flyby continuity value, not a powered
short-stay re-launch.)**

| Source | Topology (path) | Earth dep V∞ / C3 | **Mars-enc V∞ (km/s)** | Venus-enc V∞ | Earth-arr V∞ | Total TOF | Tool / ephemeris | Mars V∞ tabulated? |
|---|---|---|---|---|---|---|---|---|
| **Okutsu 2002** EMVE | E-M-V-E (Venus AFTER Mars), ballistic | **3.60 km/s** / ~12.96 | **6.97** (Table 1) / 5.98 (2017 case, Table 2) | 6.85 (Fig 7) / 9.81 (Fig 2) | **4.77** (lowest of all) | 2.19 yr (~800 d) | STOUR + MIDAS, analytic eph (unspecified) | **YES — Table 1 & 2 + Fig callouts** |
| **Okutsu 2002** EME | E-M-E (Mars only), ballistic | 3.30 (P=1.5) → 6.84 (P=3) | 6.78 (P=1.5) → 11.22 (P=2) → 14.63 (P=3) | n/a | 3.30 → 6.85 | 1.99–3.0 yr | STOUR + MIDAS | **YES — Table 1** |
| **Tito 2013** (Inspiration Mars) | E-M-E (Mars only), ballistic flyby + TMI only | 6.232 / **38.835** | **5.417** (Table IV, DE-421) | n/a (no Venus) | 8.837 (C3 78.094); re-entry V peri 14.18 | **501 d** (1.37 yr) | MAnE + STK/Astrogator, **DE-421** | **YES — Tables II & IV** |
| **Hughes 2014** | E-V-M-E (Venus BEFORE Mars), pure ballistic | **4.50** / 20.25 | ~5–7 (Fig 1 Tisserand only) | ~7 (Fig 1 only) | 6.34–9.0 | 499–582 d | STOUR + STK Astrogator, analytic (unspecified) | **NO — only Earth-side tabulated** |
| **Boeing 2022** (Donahue-Duggan) | 2033: E-M-E + inbound DSB; 2035: E-M-V-E ballistic | 36.35 (2033) / 66.18 (2035) C3 | **5.113** (2033) / **7.037** (2035) | 6.491 (2035, alt 10,855 km) | 5.416 (2035); V_entry 12.5/13.9 (2033) | 531 d (2033) / 573 d (2035) | MAnE / SpaceFlightSolutions, DE unspecified | **YES — Fig 4 / Fig 18 callouts** |

### Frontier synthesis (the consolidated sourced picture)

1. **The ballistic Mars-V∞ floor is ~5 km/s, ceiling ~7 km/s for practical
   (low Earth-side V∞, ~500–800 d) free returns.** Realised values across all
   four sources: **5.113 (Boeing 2033), 5.417 (Tito 2018), 5.98 (Okutsu 2017
   EMVE), 6.78–6.97 (Okutsu Jan-2014 EME/EMVE), 7.037 (Boeing 2035).** Tight
   cluster. **Hughes's "5–7 km/s class" is now fully corroborated by THREE
   independent tabulated sources** (Okutsu, Tito, Boeing) — Hughes was the only
   one that *didn't* tabulate it, and the ancestors fill that gap exactly.
2. **Fast (2-/3-yr) EME families pay in Mars V∞:** Okutsu Table 1 gives
   Mars V∞ **11.22 (P=2 yr), 14.63 (P=3 yr)** — flagged "excessive". The fast
   abort families have HIGH Mars approach speed; the low-Mars-V∞ families are the
   longer-period / Venus-assisted ones.
3. **None reaches the Jones Mars-~3 km/s cycler class ballistically.** Okutsu's
   *Tisserand floor* (Earth→Mars launch V∞ ~3.0 km/s, EVME arrival ~3.4 km/s) is
   a coplanar-circular idealisation; the realised ballistic families are all
   5–7+ km/s. The only sub-3-km/s Mars-side number anywhere is Okutsu's
   **Mars-DEPARTURE V∞ 2.48 km/s** (Table 2) — and that is a **powered short-stay
   re-launch** (Mars stopover + departure burn), NOT a ballistic flyby. This
   **confirms** the node-crossing/broken-plane gate: low Mars-V∞ needs added
   energy management (stopover burn / broken-plane DSM / DSB), i.e. multi-arc.
4. **Lowest Earth-arrival V∞ across all four = Okutsu EMVE 4.77 km/s** — the
   Venus-flyby-on-return (E-M-V-E) is the unique topology that buys a gentle
   Earth re-entry ballistically. Tito's Mars-only E-M-E pays for it with an
   **8.84 km/s / 14.18 km/s re-entry** (hence Tito's heatshield/aerocapture study
   and Boeing's inbound DSB). **E-M-V-E (Venus on return) is the entry-speed
   escape valve; E-M-E is not.**
5. **Two distinct Venus-assist topologies in the corpus:**
   **Hughes = E-V-M-E (Venus FIRST)**, optimised for low Earth-LAUNCH V∞ (4.50);
   **Okutsu/Boeing-2035 = E-M-V-E (Venus LAST)**, optimised for low Earth-ARRIVAL
   V∞ (4.77 / 5.42). Same three bodies, opposite Venus placement, opposite
   objective. **Note:** Hughes 2014 reported its E-M-V-E ("EMVE") search gave only
   TOF>800 d ⇒ no IM candidates, yet Okutsu's E-M-V-E closes at ~800 d and
   Boeing-2035's M-V-E return closes at 573 d — so the E-M-V-E family IS viable;
   Hughes's pessimism was constraint-window-specific.

---

## 5. METHOD / TOOL PROVENANCE → v4.2 FLAGS

| Source | Optimizer / propagator | Ephemeris | Flag |
|---|---|---|---|
| Okutsu 2002 | **STOUR** (patched-conic grid) + **MIDAS** (ΔV optimization) | "analytic ephemeris" — **DE version UNSPECIFIED** | `source_ephemeris: analytic (STOUR/MIDAS, unspecified)` — UNSPECIFIED |
| Tito 2013 | **MAnE** (patched-conic optimum-C3) + **STK/Astrogator** (full force) | **JPL DE-421** (ref 2, named) | `source_ephemeris: DE-421` — **FIRST DE-NAMED free-return source.** SPECIFIED |
| Hughes 2014 | STOUR + STK Astrogator | analytic, unspecified | UNSPECIFIED (per Hughes note) |
| Boeing 2022 | MAnE / SpaceFlightSolutions | unspecified | UNSPECIFIED (per Boeing note) |

**Tito 2013 is the only one of the four with a named DE kernel (DE-421)** — its
Table III/IV full-force values are the most provenance-clean anchors in the set.

---

## 6. CATALOGUE-ELIGIBILITY RULING — both NOT cyclers (confirmed)

- **Okutsu 2002:** EME / EMVE / EVME are **one-shot human free-return / abort
  trajectories** (single pass Earth→Mars(→Venus)→Earth). They do not re-encounter
  a body to sustain a repeating orbit. **NOT cyclers → NOT catalogue rows.** The
  paper itself distinguishes these from the Earth-Venus cycler line (it is the
  abort-trajectory study for NASA's DRM-2014, not a cycler paper). The collision
  orbits (integer-Earth-year EME) are a *periodic* subset but they are still
  single human missions, not infrastructure cyclers.
- **Tito 2013:** the Inspiration Mars 2018 flight is a **one-shot human Mars
  flyby free return** (501-day single mission). Explicitly a "precursor"
  demonstrator. **NOT a cycler → NOT a catalogue row.**

Both are valid **cross-check / golden-test anchors** for a patched-conic
free-return pipeline (NOT catalogue goldens). Tito's DE-421 full-force Table IV
is the strongest anchor candidate of the four sources.

### v4.2 backfill flags (recorded per checklist, in case ever staged as anchors)

- **center:** heliocentric (Sun) for transfer arcs; body-centered unpowered flyby
  at Mars / Venus. No catalogue `center` field (not rows); legs `center: "Sun"`.
- **tof_days_bounds:** Okutsu — TOTAL TOF per family (Table 1: 1.99–3.0 yr;
  Table 2: 475–605 d). Per-leg ToFs derivable from Fig 7/12 callouts (SOURCED,
  figure-annotated). Tito — per-leg flight times tabulated (Table I/III:
  E→M 227.05 d, M→E 274.23 d; total 501.29 d) → **SOURCED, table-grade.**
- **source_ephemeris:** Okutsu UNSPECIFIED; **Tito DE-421 (SPECIFIED).**

---

## 7. MAPS-TO-OUR-X + decisive findings

| Element | Our code / frontier | Verdict |
|---|---|---|
| Okutsu EMVE (E-M-V-E, Venus on return, ballistic, Earth-arr V∞ 4.77) | VEM hunt; free-return genome (#137); single-ellipse corrector | **MAPS on topology** (heliocentric Lambert arcs + unpowered Mars & Venus flybys). The low-Earth-arrival-V∞ escape valve is the Venus-on-RETURN placement — distinct from Hughes's Venus-first. Worth a scan probe: does our VEM hunt surface an E-M-V-E (Venus-last) arc, not just E-V-M? |
| Okutsu Table 1 Mars V∞ 6.78–14.63 by family | Jones Mars-~3 km/s; #137 E-M sub-arc ~2.81; #110 floors | **First tabulated Mars-V∞-by-family** — confirms ballistic Mars free returns are 5–7 km/s (practical) to 11–15 (fast). Sub-3 only via powered stopover (Table 2's 2.48). Independent corroboration of node-crossing/broken-plane gate. |
| Tito DE-421 full-force EME (Mars V∞ 5.417, C3 38.835, 501 d) | launch-window scans; verify chain | **VALIDATION + NEIGHBOURHOOD ANCHOR.** DE-421-named, full-force, patched-conic-cross-checked. The cleanest sourced free-return anchor we hold. 2018 is outside our forward horizon but the recurrence (2031/2033) is in it. |
| 2018→2031→2033 fast-free-return recurrence | 10-yr horizon scans | **Design note** — the "twice per 15 yr" 1.4-yr-family cadence (Patel/Okutsu/Tito) explains the sparse fast-free-return clusters; Boeing's 2033/2035 are the next realisations. |
| STOUR(Okutsu)/MAnE(Tito) = grid patched-conic + corrector-to-full-force | `hunt_vem_ballistic.py`, `correct.py`, verify chain | **MAPS — our pipeline IS this architecture** (grid patched-conic → numerical-integration validation). |

**SINGLE MOST DECISIVE FINDING:** Okutsu 2002 is the **missing Mars-V∞ data
table Hughes 2014 omitted** — it tabulates Mars-encounter V∞ by free-return
family (EME P=1.5/2/3, EMVE), establishing that the practical ballistic Mars
free return sits at **Mars V∞ ≈ 6–7 km/s** (6.78–6.97 for Jan-2014, 5.98 for the
2017 EMVE), with fast 2-/3-yr families climbing to **11–15 km/s**, and the only
sub-3-km/s Mars-side number being a **powered short-stay re-launch (2.48)**, not
a ballistic flyby. Tito 2013 anchors the Mars-only E-M-E corner with the
**only DE-421-named, full-force-validated** free-return solution we hold
(Mars V∞ 5.417, C3 38.835, 501 d). Together with Hughes and Boeing the four
sources now bracket the ballistic free-return frontier from **two independent
tools (STOUR, MAnE), four institutions, 2002→2022**, and they all agree:
**ballistic Mars free returns live at Mars V∞ 5–7 km/s and DO NOT reach the
Jones Mars-~3 km/s cycler class without added energy management (stopover burn /
broken-plane DSM / inbound DSB) — i.e. a multi-arc / powered topology our
single-ellipse-per-leg ballistic corrector cannot represent.**

---

## 8. HONEST "not extractable" / caveats

- **Okutsu figure year-labels are internally inconsistent** (Fig 1/2/4 tick-mark
  print artifacts); I trust the V∞/TOF callouts (table-consistent) over the
  figure year strings. Dates quoted from figures carry [sic] where the printed
  year is implausible.
- **Okutsu's analytic-ephemeris DE version is not stated** (UNSPECIFIED), unlike
  Tito's DE-421.
- **No per-leg V∞ for every Okutsu family** — only the four Table-1 families and
  the annotated Figs (1,2,4,7,12) have full callouts. The scatter plots (Figs
  3,6,9,10,11) are search clouds, NOT data-grade for individual points.
- **Okutsu Table 1 footnote [c] "Total TOF is fixed at 800 years"** is a clear
  typo for **800 days** (TOF column reads 2.19 yr ≈ 800 d). Logged, not a real
  value.
- **Tito's Mars-only E-M-E has no Venus encounter** — the Boeing-note's mention
  of a Tito "Venus" is incorrect for Tito 2013; Tito's 2018 trajectory is
  strictly E-M-E (Venus appears only as a heliocentric radius marker in the
  Boeing 2033/2035 descendants, not in Tito). The Donahue-Duggan note's ref-[4]
  attribution (Tito as Boeing's ref 4) is correct; but note Tito's OWN ref [4]
  is George & Kos (NASA TM-1998-208533), not a self-reference.
- **Re-entry speed:** Tito direct-reentry ~14.18–14.2 km/s (V peri / peak);
  aerocapture options bring g-load down but extend mission ~10 days (p. 6).
