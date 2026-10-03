# `#885` — Bridge from a published Uranus tour end-state into a catalogued quasi-cycler

Task `#885` (ledger entry in `data/OUTSTANDING.md`). Results note. Sections 1 and 2 were written and
committed BEFORE any bridge or control was computed; the commit that adds them is the
pre-registration record. Later sections are appended below them, and Sections 1 and 2 are not
edited after results exist (a later opinion on them goes in a separate section).

## 1. Pre-registered threshold (written before computing any bridge)

**What had been computed when this was written.** Only the circular-coplanar Tisserand points
(semi-major axis, eccentricity, periapsis, apoapsis) implied by the published V-infinity pairs and by
the catalogued rows, plus per-flyby maximum bend at 50 km. Those show the size of the V-infinity gap
but not what a bridge costs. No leveraging leg, Lambert leg, control or search had been run.

**The bridge is "worth a row" (a `precursor_mga` row with `inserts_into` the quasi-cycler) if BOTH:**

- **total deterministic delta-v <= 150 m/s**, counted from the last published flyby to the first
  encounter of the quasi-cycler, including any entry correction needed beyond the bend a flyby can
  supply; and
- **duration <= 365 days** over the same span.

**A "strong" bridge** (worth stating in the row as an easy add-on) is one with total delta-v <= 60 m/s
and duration <= 180 days.

**Per-manoeuvre limit (a constraint, not part of the verdict):** no single manoeuvre above 60 m/s
(Landau, Davis & Karimi 2023, Table 2, "max. single delta-V 60 m/s"). If the best bridge needs a
larger single burn, it is reported as violating that published limit, and the verdict is still
judged on the two totals above.

**Derivation.** The bridge is an extra phase appended to an existing tour, so it is measured
against what the published Uranus Orbiter and Probe tours pay for comparable things:

- McAdams et al. 2011 (AAS 11-188), Table 5: targeting delta-v to add each new moon to the tour is
  89.0, 102.8, 103.0, 140.9 and 183.3 m/s (619 m/s for the whole 424-day tour). Adding one new
  target to a published tour costs 89-183 m/s; 150 m/s is inside that range, about one such leg.
- Landau et al. 2025 (IEEE Aerospace), Table 7: moon-tour delta-v2 about 0.73-0.78 km/s including
  statistical, margin and attitude control. 150 m/s is about a fifth of a whole moon tour's budget,
  and about 5 percent of the ~2.8 km/s total mission delta-v of the chosen design.
- Landau, Davis & Karimi 2023 (AAS 23-460), Figure 6: the in-tour manoeuvres between the ten moon
  flybys sum to 249 m/s by my addition (23, 59, 15, 19, 15, 8, 5, 45, 60). 150 m/s is about 60 percent
  of that tour's inter-flyby spending; the "strong" level, 60 m/s, is the published single-manoeuvre
  cap, i.e. a bridge that could be flown as one maximum-size published burn in total.
- Ellison et al. 2025 (AAS 25-668): 13.37 m/s deterministic for the 1.69-year science tour, 18.32 m/s
  for the whole 3.17-year moon phase. This tour shows the published designs CAN be nearly ballistic;
  a 150 m/s bridge would be about eight times that tour's deterministic spending, which is why
  150 m/s and not more is the upper edge.
- Duration: Landau 2023 caps the tour at 500 days since orbit insertion; Ellison's science tour is
  1.69 years; Landau 2025's example tour spans 703.5 days of flybys and its magneto-tour adds about
  2 years. A bridge longer than a year is comparable to a whole published tour phase, so 365 days is
  the upper edge, and 180 days (about a quarter of the Landau 2025 tour span) the "strong" level.

## 2. Pre-registered positive control (written before computing it)

The leveraging primitive used for the bridge is Landau, Davis & Karimi's own tour step (their
seven-step algorithm): ballistic flyby at the minimum altitude, Kepler propagation to the next
apoapsis, an impulsive manoeuvre there, and a Lambert arc of at most one revolution to the next moon.
The control applies that same code to three equatorial legs of their Figure 6 (after the
inclination-reduction phase), with the 100 km altitude of that source:

| Leg | Published | delta-v | Flight time |
|---|---|---|---|
| C1 | Ariel 3.990 -> Ariel 3.897 km/s | 15 m/s | 26 d (03/04 -> 03/30/2050) |
| C2 | Umbriel 4.095 -> Umbriel 4.070 km/s | 5 m/s | 23 d (04/23 -> 05/16/2050) |
| C3 | Ariel 3.897 -> Umbriel 4.095 km/s | 8 m/s | 24 d (03/30 -> 04/23/2050) |

The incoming V-infinity direction is not published, so the control minimises over it (and over the
inbound/outbound crossing and flight time +-1 day, because the dates are whole days). **PASS** if,
for all three legs, our minimum delta-v <= published + 3 m/s, AND for at least two of the three it
is >= 0.5 x published - 3 m/s (so a primitive that finds everything free does not pass). Anything
else is FAIL and is reported as such.

**Known limitation of the existing machinery, stated before use.** The Campagnola-Russell phase-free
VILM functions in `cyclerfinder.search.vilm` (Gamma, the Eq. 13 quadrature, `vilm_dv_min`, the
Gamma-floor flag used by `leveraging_leg` and `leveraging_chain`) are not used here. The exterior
Gamma of Eq. 25 has a pole near adimensional V-infinity 0.41 and a zero near 1.38, and every
V-infinity in this problem sits in or across that band (Ariel 4.325/5.509 = 0.79, Oberon
3.448/3.151 = 1.09, cycler Ariel 0.28, cycler Oberon 0.58). Those formulas were derived and
validated for small adimensional V-infinity (Europa endgame, 0.06-0.13). Only plain Keplerian
arithmetic and Lambert arcs are used.

---

*Sections 3 onward were written after the computations. Sections 1 and 2 above are unchanged
from the pre-registration commit `a6fbf9f2`.*

## 3. What was run

Code: `src/cyclerfinder/search/uranus_bridge_885.py` (planar Kepler, the Landau-Davis-Karimi tour
step `landau_leg`, flyby rotation against `core.flyby.max_bend`, Tisserand points, the phase-free
flyby ladder, the two-moon cycler entry geometry), `scripts/screen_885_uranus_tour_to_cycler_bridge.py`
(stages `turns`, `control`, `ideal`), `tests/search/test_uranus_bridge_885.py` (7 tests, under 1 s).
Outputs in `data/found/885_uranus_bridge/`: `turn_check_six_rows.jsonl`,
`control_landau2023_fig6.jsonl`, `ideal_model.jsonl`. Model throughout: circular coplanar moons about
Uranus with the registry constants of `core/satellites.py`, zero-sphere-of-influence patched conics.
For the control the moon longitudes were anchored to URA111 at the published dates, projected onto
the Uranian equatorial plane.

## 4. Headline: the catalogued quasi-cyclers' flybys cannot be flown ballistically

While setting up the cycler entry state I computed the turn each catalogued row demands at its
encounters. **In all six Uranian (1,1) rows the demanded turn exceeds the maximum ballistic bend
at 50 km by a factor of 1.9 to 28.** For the bridge target, `ariel-oberon-1-1-uranian-quasi-cycler-2026`,
the spacecraft arrives at Oberon on one leg and must leave on the next with its V-infinity turned
by 122.8 deg; Oberon can turn a 1.83 km/s V-infinity by 8.07 deg at 50 km (8.56 deg at the
surface). At Ariel the demand is 36.1 deg against 6.22 deg.

| Row | Middle moon: V-inf, demanded / max (50 km) | Closing moon: V-inf, demanded / max | Impulse beyond bend per lap (km/s) |
|---|---|---|---|
| umbriel-oberon (#312) | Oberon 0.960, 97.0 / 24.9 deg | Umbriel 0.895, 45.3 / 16.5 deg | 1.13 + 0.46 |
| titania-oberon | Oberon 1.968, 152.6 / 7.0 | Titania 2.162, 177.2 / 6.3 | 3.76 + 4.31 |
| ariel-umbriel | Umbriel 1.300, 155.7 / 8.4 | Ariel 0.979, 102.7 / 14.0 | 2.50 + 1.37 |
| ariel-titania | Titania 1.719, 131.8 / 9.6 | Ariel 1.231, 17.9 / 9.3 | 3.01 + 0.19 |
| ariel-oberon | Oberon 1.829, 122.8 / 8.1 | Ariel 1.521, 36.1 / 6.2 | 3.08 + 0.79 |
| umbriel-titania | Titania 1.006, 116.8 / 24.4 | Umbriel 1.230, 137.5 / 9.4 | 1.45 + 2.21 |

The impulse column is an estimate: the smallest impulse at the sphere of influence that completes
the turn after the maximum ballistic bend (`entry_correction_kms`). A burn at flyby periapsis would
be somewhat cheaper, but every entry is km/s-class per lap, and capture-and-re-escape at the moon
is also km/s-class, so no patched-conic repair is cheap.

**Why this is believed (checks run):**

- The rows were rebuilt in their own construction (`v2_moontour` convention: sorted moons,
  `phase0` = 30 deg, `rel_offset`, each leg the row's Lambert at its `n_rev`) and every encounter
  V-infinity reproduces the evidence files to four decimals, including #312's 0.9199 / 0.9604 /
  0.8947.
- Two independent routes give the same turn. From the Lambert vectors directly, and from the
  geometry: a symmetric closure arrives on one crossing of the middle moon's orbit and leaves on
  the mirror crossing, so the demanded turn is 2 min(pump, 180 deg - pump) for the conic's pump
  angle there (Oberon 118.6 deg gives 122.8 deg; Ariel 18.1 deg gives 36.1 deg). The identity is
  a test (`test_mirror_closure_turn_identity`).
- A flyable control: Lambert legs on either side of a point of a single conic (the Ellison
  Oberon/Ariel orbit continued through Ariel) give a demanded turn of 4e-13 deg, as they must.
  No published flyable moon cycler was reconstructed as a further control.
- Why the gauntlet did not catch it: V2 (`data/validation/v2_moontour.py::_cycle_residual`) scores
  only the V-infinity MAGNITUDE continuity at each flyby. V4 (`data/validation/v4_uranus.py`, the
  leg loop at line 555: "IC for the V4 leg: ... using Lambert's v-out from moon-A. Identical to
  V3's choice") and V4-strict (`data/validation/v4_uranus_strict.py`, the leg loop at line 592
  calling `_select_leg_transfer`) restart every leg from that leg's own Lambert departure velocity,
  so the velocity is never carried across a flyby. The `#324` physical gate
  (`data/rerun_324_physical_gate.jsonl`, `min_useful_bend_deg` 5) tests that the moon CAN bend
  at least 5 deg, which is capacity, not the demanded turn.
- Difference from the cited architecture (an inference from `scripts/compare_576_russell_strange_galilean.py`,
  not re-read in the paper here): in Russell & Strange's free-return cyclers the target body is
  massless and is passed on the same conic, and only the flyby body bends. The project's 2-leg
  symmetric construction makes the target body turn the trajectory too, which is where the demand
  comes from.

Consequence for the ledger (not edited here): `#879` is registered as "small" with "no
validation-tier or novelty consequence". That understates it: if the numbers above stand, the six
rows are not ballistic trajectories in the model they are validated in.

## 5. Ideal-model picture (circular coplanar, phase-free)

**Assumptions.** The published tables give V-infinity magnitudes, dates and altitudes, not vectors.
Each published stretch is taken as a single coplanar conic through its last two flybys (no
manoeuvre is listed between them), which fixes (a, e) and the pump angle up to the inbound/outbound
sign. Both tours are equatorial by then: AAS 25-668 equatorialises with the first eleven Titania
flybys before the science tour; Landau et al. 2023 split their tours into inclination reduction
then "equatorial flybys". Miranda (inclined 4.2 deg) is left out. The check that the coplanar
reading is sound: the Ellison Oberon 3.448 / Ariel 4.325 conic predicts Umbriel 4.388 and Titania
3.884 km/s, against that tour's own Umbriel 4.392-4.396 and Titania 3.880 / 3.922 (Miranda is off,
3.015 against 3.31-3.37, as expected for an inclined moon). That is evidence for the assumption, not
a control of the leveraging calculation.

| Point | a (km) | e | rp (km) | ra (km) | Period (d) | energy (km2/s2) | V-inf Ariel / Umbriel / Titania / Oberon (km/s) |
|---|---|---|---|---|---|---|---|
| AAS 25-668 Oberon 3.448 -> Ariel 4.325 | 1,079,637 | 0.886 | 123,265 | 2,036,010 | 33.9 | -2.68 | 4.325 / 4.388 / 3.884 / 3.448 |
| ariel-oberon cycler | 478,690 | 0.604 | 189,768 | 767,613 | 10.0 | -6.05 | 1.521 / 2.578 / 2.357 / 1.829 |
| Landau 2025 Umbriel 4.1 -> Oberon 2.6 | 495,761 | 0.785 | 106,638 | 884,883 | 10.6 | -5.84 | 4.331 / 4.100 / 3.257 / 2.600 |
| umbriel-oberon cycler | 444,920 | 0.407 | 264,064 | 625,777 | 9.0 | -6.51 | - / 0.920 / 1.520 / 0.960 |

**Gap, Ariel-Oberon.** Energy -2.68 to -6.05 km2/s2, angular momentum 1.161e6 to 1.328e6 km2/s;
pump angle at Ariel 82.4 to 18.1 deg, at Oberon 109.7 to 118.6 deg. A flyby moves along its own
moon's V-infinity contour; moving between contours of one moon needs another moon's flyby or a
manoeuvre. At the published end-state the moons bend very little (50 km: Ariel 0.81 deg at
4.325 km/s, Umbriel 0.91, Titania 2.0, Oberon 2.39 deg), so contour-walking is slow.

**Lower bound on delta-v: zero, given enough flybys.** Alternating flybys of Ariel and Oberon walk
the Tisserand graph down to the cycler point without any manoeuvre, with periapsis held above the
103,000 km floor (Landau 2023 Table 2) throughout: a breadth-first search over arrival states
(moon, V-infinity, pump angle), with flyby bends sampled at 0, +-1/2 and +-1 times the maximum,
reaches a state from which one Oberon flyby turns onto the cycler's departure in **16 flybys** at
50 km (17 at 100 km; Umbriel and Titania do not shorten it). The route steps V-infinity at Oberon
3.45 -> 3.29 -> 3.20 -> 3.09 and at Ariel 4.33 -> 3.85 -> 3.29 -> 2.58, then eight Ariel flybys at
2.58 km/s rotate the pump angle by about 18 deg before the final Oberon arrival at 1.83 km/s. This
is phase-free: it assumes every encounter can be phased. If each leg takes one revolution of the
orbit it is flown on, the 15 legs between the 16 flybys sum to about 323 days (355 at 100 km). So the ideal model offers
a ballistic route at about the pre-registered duration limit, before any phasing cost. For scale,
going straight from the published conic to the cycler conic with two apse-to-apse burns and no
flybys costs 603 m/s.

**Umbriel-Oberon (step 5).** From Landau et al. 2025 Table 3's last flyby (Oberon 2.6 km/s, on
the Umbriel 4.1 / Oberon 2.6 conic) the ladder reaches the umbriel-oberon cycler point in **12
flybys** (Oberon x2, Ariel x2, Oberon, Ariel x4, Oberon, Umbriel, Oberon), about 117 days at one
revolution per leg; the two-burn no-flyby reference is 806 m/s. It looks CHEAPER in time and flyby
count than the Ariel-Oberon bridge, because the Landau end-state is already a 10.6-day orbit and
bends grow as V-infinity falls. (The same turn defect applies to its target: section 4.)

## 6. Positive control of the leveraging step (pre-registered in section 2)

Stage `control`, 100 km altitude, rp floor 103,000 km, incoming pump angle on a 0.5 deg grid with
both crossing signs, 9 bend samples, flight time +-1 day at 0.02 d; departure at 12:00 UTC on the
published date. The hour is not published; for C3 (two different moons, so the relative moon phase
matters) the hour was scanned over 0, 6, 12 and 18 UTC. That scan was added after a coarse
exploratory run, so the 12:00 value is reported as well.

| Leg | Published | Ours, grid (12:00) | Ours, grid min over hours | Exact re-solve of the grid optimum |
|---|---|---|---|---|
| C1 Ariel 3.990 -> 3.897, 26 d | 15 m/s | 8.6 | 8.6 | 7.5 m/s (tof 26.75 d) |
| C2 Umbriel 4.095 -> 4.070, 23 d | 5 m/s | 6.5 | 6.5 | 3.9 m/s (tof 22.67 d) |
| C3 Ariel 3.897 -> Umbriel 4.095, 24 d | 8 m/s | 9.5 | 9.5 (hours 0/6/18: 12.6, 10.8, 11.1) | 9.4 m/s (tof 24.22 d) |

**Verdict: PASS** by the pre-registered criterion (all three <= published + 3 m/s; all three >= 0.5 x
published - 3 m/s). The burns sit at apoapsis, 1.57-1.73 million km, as in the published method.
Disclosed: before the stage run, two coarse exploratory grids (2 deg pump step, 0.05-0.1 d) gave
C1 13.9-14.8, C2 47.9-56.0 and C3 118 m/s. The minima are narrow in pump angle and flight time, and
the coarse grids missed them. The control is one-sided by nature: our minimum is over an incoming
direction the paper does not give, so ours below the published value (C1) is expected.

## 7. Patched-conic realisation (plan step 3): not done

Not built, deliberately: section 4 means the entry state would put the spacecraft on a trajectory
that cannot continue past its next encounter without a km/s-class burn, so pricing the entry is
moot. The pieces exist and are tested (`landau_leg`, `CyclerGeometry.entry_targets` for the
cycler-compatible entry epochs at either moon, `entry_correction_kms` for the entry mismatch).
What remains is a beam search over (bend, next moon, arrival time) with the ladder as the
heuristic: about half a day to write and an estimated 30-60 minutes on 4 workers, staged by depth.
It is worth doing only for a target whose flybys are flyable.

## 8. Verdict against the threshold

**Not worth a row.** The decisive reason is the target, not the bridge. The ideal model says a
ballistic bridge to the ariel-oberon cycler point exists in about 16 flybys and roughly a year of
phase-free flight, which is at the pre-registered duration limit (365 days) with zero manoeuvre
delta-v; no realised trajectory was computed, so the threshold was not tested on a real bridge.
But the quasi-cycler it would insert into demands turns 15 times (Oberon) and 6 times (Ariel) what
the moons can supply, about 3.9 km/s per lap by the estimate in section 4. A `precursor_mga` row
that `inserts_into` that row would be inserting into something that is not a ballistic orbit.

Separately, AAS 25-668 Table 7 requires at least 30 days between encounters (downlink). The
ariel-oberon row's legs are 7.75 days and umbriel-oberon's 14.94 days, so the cyclers break an
operations constraint of the very tour the bridge starts from. Landau et al. 2025 impose no such
rule (their tour has a 9.1-day interval).

## 9. Other things found on the way (for the coordinator; nothing edited)

- `cyclerfinder.search.vilm` (Gamma, the Eq. 13 quadrature, `vilm_dv_min`) and everything built
  on its Gamma floor (`leveraging_leg.gamma_floor_ok`, `leveraging_chain.walk_vinf_down`,
  `endgame_graph`) are not valid at Uranian V-infinity: the exterior Gamma has a pole between
  adimensional V-infinity 0.40 and 0.42 (48.97 and -122.9) and is negative from there to about
  1.38; the interior Gamma has a pole between 0.79 and 1.09.
- `data/empty_regions.jsonl` entry `uranus-neptune-regular-moon-endgame-vilm-2026-06-23` reads
  "contours disjoint at every probed vinf 4-15 km/s" for Uranian moon pairs. Both published tours
  link Ariel, Umbriel, Titania and Oberon at 2.6-4.5 km/s, and the ladder above links Ariel and
  Oberon ballistically. The entry's interpretation looks wrong; its `linkable` test was not
  diagnosed here.
