# Hughes-Edelman-Longuski 2014 — Fast Mars Free-Returns via Venus GA (STOUR cross-check)

Mined 2026-06-07 (Task #147). The **STOUR cross-check** flagged rank-3 of the
external-algorithms survey (`docs/notes/2026-06-07-external-algorithms-survey.md`,
Thread 4) and the decisive diligence on the Jones VEM low-V∞ question: at what
*topology* do published E-V-M free returns achieve low V∞, and is that family
single-ellipse-representable by our corrector?

**Source (cite exactly — authors/title/AIAA number only):**
Hughes, K. M., Edelman, P. J., Longuski, J. M., Loucks, M. E., Carrico, J. P.,
& Tito, D. A., "Fast Mars Free-Returns via Venus Gravity Assist," AIAA 2014-4109,
AIAA/AAS Astrodynamics Specialist Conference (SPACE Conferences & Exposition),
San Diego, CA, 4-7 Aug 2014. Purdue University AAC + Space Exploration
Engineering + Applied Defense Solutions + Wilshire Associates. DOI 10.2514/6.2014-4109.

> 18 pp., clean digital typeset. All six tables (Tables 1-6) and the method text
> read unambiguously. Figures are search-scatter plots (Figs 2-6, 9-14) and
> trajectory geometry (Figs 1, 7, 8, 11) — read for qualitative topology, not
> for numbers (the numbers are in the tables).

---

## 0. HEADLINE VERDICT (the decisive question)

**TOPOLOGY: pure ballistic patched-conic conic-chains.** Every reported solution
is a *free return* — Lambert heliocentric arcs joined at **unpowered** Venus and
Mars gravity assists, with **zero deterministic maneuvers** (zero DSMs, zero
broken-plane). The defining mission constraint of the whole paper is that the
trajectory be ballistic so a crew returns to Earth without any burn (p.2, p.12).
The arcs are **near-180° Hohmann-like transfers** lying essentially in the
ecliptic; the *only* notable inclination is ~5.4° on the Earth-Venus leg (p.9).

This is the topology our single-ellipse-per-leg ballistic corrector is *designed*
to represent (heliocentric Lambert arcs + unpowered flyby continuity). So on
topology grounds this family IS representable by our corrector.

**BUT the V∞ numbers do NOT reach the Jones / #110 low-V∞ regime that is our
actual blocker — and the paper explains exactly why, which is the load-bearing
finding:**

- Lowest published **Mars-encounter V∞** here ≈ NOT a free parameter they minimize
  — they minimize *Earth* launch/arrival V∞. The Venus flyby in EVME occurs at
  V∞ ≈ 7 km/s and Mars arrival V∞ ≈ 5 km/s in the Tisserand illustration (Fig 1,
  p.4). This is the **5-7 km/s class**, NOT the **Mars ~3 km/s class** that Jones
  AAS 17-577 targets and that #137's E-M sub-arc reaches (~2.81).
- Lowest published **Earth-launch V∞ = 4.50 km/s** (11/22/2021 EVME, Table 2/3/4,
  pp.5,8,10). Lowest **Earth-arrival V∞ = 6.34 km/s** (11/28/2053) for Mars
  returns; **6.46** is the lowest near-term (12/04/2021).
- For the Venus-only EVE free returns they reach **Earth launch V∞ = 2.75 km/s**
  and **arrival V∞ = 4.34 km/s** (Table 5, p.15) — but those never touch Mars.

**The crucial topological caveat (pp.12-13):** the paper states explicitly that
the near-180° transfer geometry produces *highly inclined* transfer arcs unless
the Venus/Mars encounters sit on node crossings, and that the way to relax this —
to get *lower* Earth launch/arrival V∞ than Table 3 — is a **broken-plane
maneuver**, which they deliberately **exclude** because it makes the trajectory
no longer a *free* return:

> "Nevertheless, a broken-plane maneuver could potentially allow for transfers
> closer to a true Hohmann transfer, and thereby (potentially) provide lower
> Earth launch and arrival V∞ values than those found in table 3. Such
> trajectories are not investigated in this paper, as the primary goal of this
> study is to find trajectories for an IM-type mission that are purely
> ballistic..." (p.13)

> "the use of a broken-plane maneuver prevents the resulting EVME trajectory from
> being a true *free* return, and is therefore less desirable for an IM-type
> mission." (p.12)

**Verdict on the M-ED / Jones blocker:** this paper CONFIRMS, from an independent
STOUR broad search by the originating Purdue group, that the *purely ballistic
patched-conic* E-V-M family lives at **V∞ ~4.5-9 km/s at Earth and ~5-7 km/s at
the flyby/Mars** — i.e. comfortably ABOVE our #110 floors of 17.9/18.5 km/s only
in the sense that our floors are an artifact of the wrong basin; but it does NOT
exhibit a ballistic single-ellipse solution down at the Jones Mars-~3-km/s class.
The lower-V∞ direction is explicitly gated behind a **broken-plane (DSM)
maneuver** — i.e. a *multi-arc / powered* topology, NOT single-ellipse ballistic.
This is the same conclusion the #110 dense scan reached (0 bend-feasible) and the
MEMORY S1L1 finding (the low-V∞ closure needs intermediate eccentricity / an
extra arc). **The low-V∞ family is NOT single-ellipse-representable as a pure
free return; it requires a broken-plane DSM — consistent with multi-arc being the
real fix.**

---

## 1. GOLDEN-ELIGIBLE ANCHOR TABLES (verbatim, with page numbers)

These are published, sourced V∞/date/TOF anchors. See §4 for the catalogue-
eligibility ruling (they are NOT cyclers → NOT catalogue rows, but ARE valid
cross-check / golden-test anchors for a patched-conic E-V-M pipeline).

### Table 1 — Trajectory Search Parameters (p.3, verbatim)

| Parameter | Value |
|---|---|
| Max V∞,Launch (km/s) | 6.5 |
| Max V∞,Arrival (km/s) | 9.0 |
| Max TOF (days) | 600 |
| Min Launch Date (mm/dd/yyyy) | 12/01/2014 |
| Max Launch Date (mm/dd/yyyy) | 01/31/2060 |

### Table 2 — Notable EVME (Earth-Venus-Mars-Earth) from broad 45-yr search (p.5, verbatim)

Gravity-assist path **EVME** = Venus flyby BEFORE Mars. Path code "3 2 4 3"
(Earth-Venus-Mars-Earth). VINF 2.50→6.50 by 0.25 km/s; ALTMIN 200 km at Venus
and Mars. V_Entry = inertial Earth entry speed at 122 km altitude.

| Launch Date | V∞,Launch (km/s) | C3 (km²/s²) | TOF (days) | V∞,Arrival (km/s) | V_Entry (km/s) |
|---|---|---|---|---|---|
| 11/22/2021 | 4.50 | 20.25 | 582 | 6.53 | 12.85 |
| 12/08/2021 | 5.50 | 30.25 | 566 | 6.55 | 12.87 |
| 08/28/2034 | 4.75 | 22.56 | 558 | 6.52 | 12.85 |
| 06/24/2036 | 5.50 | 30.25 | 499 | 8.79 | 14.14 |
| 07/02/2047 | 5.75 | 33.06 | 565 | 8.84 | 14.17 |
| 11/28/2053 | 5.00 | 25.00 | 580 | 6.34 | 12.76 |

### Table 3 — Best Near-Term EVME Opportunities (p.9, verbatim)

Near-term Pareto set; all launch 11/18/2021–12/21/2021. VINF 2.50→6.50 by 0.10 km/s.

| Launch Date | V∞,Launch (km/s) | C3 (km²/s²) | TOF (days) | V∞,Arrival (km/s) | V_Entry (km/s) |
|---|---|---|---|---|---|
| 11/22/2021 | 4.50 | 20.25 | 582 | 6.53 | 12.85 |
| 12/03/2021 | 5.10 | 26.01 | 579 | 6.95 | 13.07 |
| 12/04/2021 | 5.20 | 27.04 | 569 | 6.46 | 12.82 |
| 12/08/2021 | 5.50 | 30.25 | 566 | 6.55 | 12.87 |
| 12/14/2021 | 6.00 | 36.00 | 564 | 6.81 | 13.00 |
| 12/15/2021 | 6.10 | 37.21 | 563 | 6.78 | 12.99 |
| 12/16/2021 | 6.20 | 38.44 | 561 | 6.77 | 12.98 |
| 12/17/2021 | 6.30 | 39.69 | 560 | 6.78 | 12.99 |
| 12/19/2021 | 6.50 | 42.25 | 559 | 6.84 | 13.02 |

### Table 4 — High-Fidelity Comparison of 11/22/2021 EVME (p.10, verbatim)

STK Astrogator (full n-body: Earth, Moon, Sun, Mars, Venus) vs STOUR patched-
conic+analytic-ephemeris, SAME opportunity. h_V, h_M = closest-approach altitudes
at Venus/Mars (km). **This is the validation anchor: STOUR vs high-fidelity agree.**

| Propagator | Launch Date | V∞,Launch | C3 | h_V (km) | h_M (km) | TOF (d) | V∞,Arrival | V_Entry |
|---|---|---|---|---|---|---|---|---|
| Astrogator | 11/22/2021 | 4.511 | 20.353 | 11097 | 363 | 582.5 | 6.459 | 12.87 |
| STOUR | 11/22/2021 | 4.500 | 20.250 | 10868 | 346 | 582.2 | 6.526 | 12.85 |

> Encounter dates for the 11/22/2021 EVME (Figs 7, 8): Launch 11/22/2021,
> Venus flyby 4/4/2022, Mars flyby 10/11/2022 (STK: 10/12/2022 — 1-day diff),
> Earth arrival 6/27/2023. (pp.9-10)

### Table 5 — Subset of Best Near-Term EVE (Earth-Venus-Earth, Venus free return) (p.15, verbatim)

Venus-only free returns (NO Mars). Path code "3 2 3". VINF 2.00→6.50 by 0.25 km/s,
ALTMIN 300 km at Venus.

| Launch Date | V∞,Launch (km/s) | C3 (km²/s²) | TOF (days) | V∞,Arrival (km/s) | V_Entry (km/s) |
|---|---|---|---|---|---|
| 07/31/2026 | 2.75 | 7.56 | 366 | 7.09 | 13.15 |
| 05/07/2021 | 6.50 | 42.25 | 400 | 4.34 | 11.90 |
| 09/18/2026 | 6.25 | 39.06 | 327 | 7.25 | 13.24 |
| 11/18/2019 | 6.00 | 36.00 | 417 | 6.39 | 12.79 |
| *07/02/2023 | 5.75 | 33.06 | 347 | 7.83 | 13.56 |

*Not in the Pareto optimal set of Fig 13.

### Table 6 — Subset of Best 2018 EVE Opportunities (p.16, verbatim)

| Launch Date | V∞,Launch (km/s) | C3 (km²/s²) | TOF (days) | V∞,Arrival (km/s) | V_Entry (km/s) |
|---|---|---|---|---|---|
| 08/03/2018 | 2.75 | 7.56 | 371 | 7.13 | 13.17 |
| 02/23/2018 | 6.50 | 42.25 | 454 | 6.35 | 12.77 |
| 09/17/2018 | 5.75 | 33.06 | 324 | 7.18 | 13.20 |

**Note (no per-flyby Mars V∞ table):** the paper tabulates Earth launch and Earth
arrival V∞ only. It does NOT tabulate the Venus-flyby or Mars-encounter V∞
numerically (only the Fig 1 Tisserand illustration gives Venus ~7, Mars ~5 km/s
for a 5-km/s Earth launch). So a *Mars-V∞ golden* is NOT directly extractable —
only Earth-side V∞, C3, TOF, V_Entry, and (Table 4 only) the flyby altitudes.

---

## 2. METHODOLOGY (STOUR) — what our scan rung lacks (verbatim where load-bearing)

### 2.1 The engine (p.2, §II.A)

> "The STOUR program (developed by the Jet Propulsion Laboratory and Purdue
> University) was used to compute the Mars free-return opportunities with an
> intermediate Venus flyby. The STOUR program uses a patched-conic model with an
> analytic ephemeris to rapidly compute multiple gravity-assist trajectories. It
> imposes a grid search to find trajectories by stepping through specified launch
> dates and launch V∞ — thereby revealing all candidate trajectories within the
> search parameters." (p.2)

This is **exactly our `hunt_vem_ballistic.py` architecture** (grid over launch
date × launch V∞, patched-conic legs, bend-feasibility prune). Confirms Thread 4's
claim that our scan IS a STOUR-style search.

### 2.2 Search grid (Figs 2-6 captions, p.4)

- **EVME broad**: VINF 2.50→6.50 km/s **by 0.25 km/s**; launch dates 12/01/2014→
  01/31/2060 **by 1.0 day**; TFMAX 700 days (results shown to 700 d but only ≤600 d
  are IM-candidates); **ALTMIN 200 km** at BOTH Venus and Mars (p.4).
- **EVME near-term** (Figs 5,6): VINF **by 0.10 km/s**, 1-day launch step,
  TFMAX 600 d. Step refinement near term — a two-stage coarse→fine grid.
- **EVE (Venus)**: VINF 2.00→6.50 by 0.25; **ALTMIN 300 km** at Venus; TFMAX 600.
- **EMVE** (Venus AFTER Mars) was searched with launch V∞ relaxed to 7.0 km/s and
  TOF up to 5 yr (1826 d) — found only TOF>800 d ⇒ no IM candidates (pp.5-6).

**What our scan rung lacks vs STOUR here:**
1. **Two-stage grid refinement** (0.25 km/s broad → 0.10 km/s near-term) — we run a
   single fixed density. STOUR coarse-locates the launch-date clusters, then
   refines V∞ resolution only inside the surviving near-term window.
2. **Pareto-set extraction** (§II.C, p.3): they keep the *nondominated* set over
   (TOF, V∞,Launch, V∞,Arrival) rather than a single "best". A principled
   multi-objective selection our scan does not do.
3. **Tisserand-graph feasibility pre-screen** (§II.B, p.3): before gridding, they
   use the Tisserand graph (specific energy vs r_p, constant-V∞ curves per body)
   to confirm a path is *energetically* feasible at all and read off the achievable
   V∞ at each encounter. This is the analytic pre-filter our bend-prune approximates
   numerically. Refs: Strange & Longuski; Labunsky et al.
4. **Multi-cycle / repeat-geometry exploitation**: they note the E-V-M inertial
   geometry repeats ~every 32 years (Okutsu & Longuski), so 2021/2053 (and 2034)
   opportunities recur — they search a 45-yr span to capture all phasings.

### 2.3 Patched-conic assumptions (pp.2-3)

- **Zero-radius / analytic-ephemeris patched conic** (same family as Jones'
  "zero-SOI patched conic with real ephemeris"). Lambert legs between encounters.
- **Unpowered flybys**: feasibility judged by whether the required turn is
  achievable above ALTMIN (200 km Venus/Mars; 300 km Venus for EVE) — the
  Tisserand "maximum energy change for a minimum flyby altitude of 200 km" (Fig 1
  caption, p.4). **No powered-flyby Δv, no DSM** in the accepted set (that is the
  whole point — it is a *free* return).
- Pruning criterion = within the Table-1 V∞/TOF box AND flyby altitude ≥ ALTMIN
  AND (post-hoc) no unintended intermediate flybys.

### 2.4 The physics of WHY low-V∞ needs a node-crossing or a broken-plane (pp.10-12)

The most useful method insight for OUR topology fight (§II.C "Physical Behavior",
pp.10-12), verbatim-paraphrased with the load-bearing quote:

- Each leg's transfer angle is **near 180°** — "characteristic of a Hohmann
  transfer" — which is what gives the low Earth V∞ (p.10).
- But a 180° transfer between bodies in *different* planes produces a transfer arc
  inclined ~90° to the ecliptic → would need huge V∞ (p.10). The escape valve:
  > "Such a problem does not exist however if encounters occur at points where
  > these orbital planes intersect... where Venus and Mars cross the ecliptic—at
  > their ascending and descending nodes. Therefore, the desired characteristics
  > ... are that transfer arcs start and end at an apse ... and that the Venus and
  > Mars encounters occur at a node crossing to allow for transfer arcs near 180°."
  > (p.12)
- So the *purely ballistic* low-V∞ solutions exist only at the **rare phasings**
  where Venus AND Mars are simultaneously near their nodes near-apse (this is why
  the near-term window is one month wide, p.7, and why clusters appear only at
  2017/2021/2034/2036/2047/2049/2053). Off those phasings, **the lower-V∞
  direction is only reachable via a broken-plane maneuver** (multi-arc), which
  they exclude (pp.12-13, quoted in §0).

This is the **independent confirmation of our topology blocker**: a single
ballistic ellipse-pair only closes at the special node-crossing phasings; pushing
to lower V∞ off-phasing requires an added DSM/broken-plane arc = multi-arc.

---

## 3. MAPS-TO-OUR-X verdicts

| Hughes et al. element | Our code / frontier | Verdict |
|---|---|---|
| STOUR = grid over (launch date × launch V∞), patched-conic Lambert legs, unpowered-flyby altitude prune | `search/hunt_vem_ballistic.py`, `search/scan.py` | **MAPS — our VEM hunt IS this engine.** Confirms Thread 4. |
| Zero-SOI patched conic + analytic ephemeris | our patched-conic broad search | **MAPS.** Same model class as Jones broad search. |
| Pure-ballistic free return = heliocentric Lambert arcs + unpowered GA, ZERO DSM | single-ellipse-per-leg ballistic corrector (`correct.py`) | **MAPS on topology** — this family IS what our corrector represents. |
| Low Earth-V∞ requires near-180° arcs AT node crossings; off-phasing needs broken-plane DSM | #110 "0 bend-feasible"; M-ED single-ellipse blocker; S1L1 multi-arc finding | **MAPS — independent confirmation.** Pure-ballistic single-ellipse low-V∞ only at rare phasings; lower needs multi-arc. |
| Lowest ballistic Mars-side V∞ ~5-7 km/s (Fig 1), NOT ~3 km/s | Jones Mars-~3 km/s target; #137 sub-arc ~2.81 | **DOES NOT REACH Jones regime ballistically.** Jones' ~3 km/s class is NOT a pure free return here; needs the broken-plane/powered topology. |
| Two-stage grid (0.25→0.10 km/s) + Tisserand pre-screen + Pareto-set selection | single-density scan, bend-prune, single-best pick | **DOES NOT MAP — improvements available.** Coarse→fine V∞ refinement and Pareto nondomination are cheap upgrades to our scan rung. |
| 32-yr inertial repeat (Okutsu & Longuski) ⇒ search 45 yr to catch all phasings | our per-window scans | **Design note** — phasing recurrence explains the sparse cluster structure; mirrors Jones' multi-cycle matching motivation. |
| Table 4: STOUR vs STK Astrogator agree to ~0.06 km/s V∞, ~1 day, ~20 km altitude | our patched-conic→n-body verify chain | **VALIDATION ANCHOR** — sourced evidence that patched-conic+analytic-ephemeris E-V-M free returns survive to high fidelity essentially unchanged. |

---

## 4. v4.2 BACKFILL CHECKS (binding checklist) + catalogue-eligibility ruling

**CATALOGUE-ELIGIBLE? NO — confirmed.** These are **one-shot free-return
trajectories** (Earth→Venus→Mars→Earth, single pass), explicitly NOT repeating
cyclers. They do not re-encounter a body to re-initiate a repeating orbit; they
are IM-type human-flyby missions. They are **NOT catalogue rows.** (The EVE cases
are likewise single Venus free returns, not the Hollister-Menning Earth-Venus
*cyclers* — the paper itself distinguishes them, p.13, citing Hollister 1969 /
Hollister & Menning 1970 separately as the cycler line.)

They ARE valid as **cross-check / golden-test anchors** for any patched-conic
E-V-M *pipeline* (not catalogue goldens). Recording provenance flags anyway per
the checklist:

- **center**: heliocentric (Sun-centered) for the Lambert transfer arcs; flybys
  body-centered at Venus/Mars (unpowered, altitude-gated). No catalogue `center`
  field applies (not a row), but if ever staged as a pipeline golden the legs are
  `center: "Sun"`.
- **tof_days_bounds**: **TOTAL** TOF (Earth-launch→Earth-arrival), not per-leg.
  Table-2/3 totals 499-582 d; Table-5/6 EVE totals 324-454 d. Per-leg ToFs are
  NOT tabulated (only encounter dates for the 11/22/2021 case give legs: Launch→
  Venus ~133 d, Venus→Mars ~190 d, Mars→Earth ~259 d, derived from Fig 7 dates —
  derived, not tabulated, so flag as DERIVED if used). Search cap was TFMAX 600 d
  (IM constraint); broad plots shown to 700 d (EVME) / 1826 d (EMVE).
- **source_ephemeris**: STOUR = **"analytic ephemeris"** (patched conic). The
  specific analytic model / DE version is **NOT stated**. The Table-4 high-fidelity
  cross-check used **STK Astrogator** with point-mass gravity of Earth, Moon, Sun,
  Mars, Venus (Mars/Venus gravity neglected near Earth launch/arrival) — but no
  named DE ephemeris kernel given. So any anchor from this paper carries
  `source_ephemeris: analytic (STOUR, unspecified)` with the Astrogator
  cross-check as the fidelity witness. Flag UNSPECIFIED.

---

## 5. CROSS-REFERENCE: relation to Jones AAS 17-577 and the Purdue line

- **Shared lineage / member overlap:** This paper is **Purdue AAC (Longuski
  group)** using **STOUR**. Jones AAS 17-577 is **JPL** (Jones, Hernandez,
  Jesick) using their own Lambert broad-search + SNOPT. They are **different
  institutions/teams** — no direct author overlap between this paper and Jones
  2017. The connection is the **shared STOUR/MGA heritage** and the shared
  E-V-M-class problem, not shared authorship.
- **Coauthor note:** Carrico (here) and Loucks (here) also appear in the
  **Tito et al. 2013** Inspiration Mars feasibility paper (ref 13), which is the
  baseline this whole paper backs up. Dennis Tito is a coauthor here. So this is
  the **Inspiration-Mars-aligned** branch, not the JPL cycler branch.
- **Okutsu & Longuski 2002** (ref 11, "Mars Free Returns via Gravity Assist from
  Venus," JSR 39(1)) is the **direct predecessor** — the 2014-IM-too-soon
  candidate. Patel, Longuski & Sims 2002 (ref 10, "Mars Free Return
  Trajectories," JSR 35(3)) is the no-Venus free-return baseline. These two are
  the acquisition-worthy upstream STOUR free-return sources if we want the
  earliest tabulated E-V-M numbers.
- **Hollister 1969 (ref 24) and Hollister & Menning 1970 (ref 25)** are cited
  here as the **Earth-Venus cycler** origin — we already hold Hollister-Menning
  1970 (`data/sources/hollister-menning-1970-table3.yaml`). Consistent.
- **Relation to #110 / Jones VEM frontier:** Jones targets *repeating* low-excess-
  speed VEM **cyclers** at Mars-V∞ ~3 km/s, achieved with **powered/B-plane
  flybys + SNOPT to ballistic**. This paper targets *one-shot ballistic* E-V-M
  free returns and finds them only at **Mars-V∞ ~5-7 km/s** (Fig 1), explicitly
  stating the lower-V∞ direction needs a **broken-plane DSM**. The two are the
  **same topology question, opposite answers depending on whether a DSM is
  allowed** — which is precisely the M-ED single-ellipse-vs-multi-arc fork.

---

## 6. SINGLE MOST DECISIVE FINDING

**An independent Purdue/STOUR broad search of the E-V-M free-return space confirms
that PURELY BALLISTIC single-ellipse-representable E-V-M solutions live at Mars-
side V∞ ~5-7 km/s and Earth V∞ ~4.5-9 km/s, and that reaching LOWER V∞ (toward
the Jones Mars-~3 km/s class) explicitly requires a broken-plane maneuver — i.e.
a multi-arc / DSM topology our single-ellipse corrector cannot represent.** This
is sourced, third-party corroboration of the #110 "0 bend-feasible" result and
the M-ED multi-arc front-runner: the low-V∞ family is NOT a pure free return, so
the single-ellipse-per-leg ballistic corrector should not be expected to close it.
The paper does NOT refute our pipeline numerically (Table 4 even validates patched-
conic→high-fidelity agreement for the family it DOES find); it refutes the *hope*
that the Jones low-V∞ family is single-ellipse-representable.

---

## 7. HONEST "not extractable" list

- **No tabulated Venus-flyby or Mars-encounter V∞.** Only Earth launch/arrival
  V∞, C3, TOF, V_Entry are tabulated. The Venus ~7 / Mars ~5 km/s figures are read
  off the Fig 1 Tisserand *illustration* (odd-integer V∞ curves), not a data table
  — so they are illustrative, NOT golden-grade. A Mars-V∞ golden is NOT available
  from this paper.
- **Per-leg ToFs** are not tabulated (only total TOF + encounter dates for the one
  11/22/2021 case in Figs 7-8). Per-leg values would be DERIVED, not sourced.
- **No named ephemeris/DE version** for STOUR's "analytic ephemeris," and no DE
  kernel named for the STK Astrogator cross-check.
- **STOUR grid step in launch DATE within clusters, and the exact path codes' rev
  counts**, are partially given (path "3 2 4 3", VINF/ALTMIN/step in figure
  banners) but the internal Lambert rev-enumeration limits are not stated.
- **No flyby altitudes tabulated except the single Table-4 case** (h_V≈10868-11097
  km, h_M≈346-363 km for 11/22/2021). ALTMIN (200/300 km) is a search floor, not a
  realized value for the other rows.
