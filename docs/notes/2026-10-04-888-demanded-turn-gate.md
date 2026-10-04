# `#888` — A demanded-turn gate, its controls, and a re-screen of the two-moon enumerations

Date: 2026-10-04. Code: `src/cyclerfinder/verify/turn_gate.py`, `src/cyclerfinder/verify/turn_gate_closures.py`,
tests `tests/verify/test_turn_gate.py` (51 tests, under 2 s with `-n 4`),
script `scripts/screen_888_turn_gate_rescreen.py`, outputs `data/found/888_turn_gate/`.
Background: `#885` section 4 and its section 10 (`docs/notes/2026-10-04-885-uranus-tour-to-cycler-bridge.md`).

## 1. What the gate checks

A patched-conic chain can only be flown without propellant if, at every flyby, the body's
gravity can turn the incoming V-infinity vector onto the outgoing one. For each encounter the gate
(`demanded_turn_gate`) takes the body's GM, radius and altitude floor and the two V-infinity
VECTORS (any common frame, 2-D or 3-D) and reports:

| quantity | definition |
|---|---|
| demanded turn | angle between `v_in` and `v_out` (`core.flyby.bend_angle`) |
| available bend | `2 asin(1/(1 + r_p v^2/GM))`, `r_p = R + floor`, `v = min(|v_in|, |v_out|)` (`core.flyby.max_bend`) |
| magnitude mismatch | `abs(|v_in| - |v_out|)`, reported separately |
| ratio | demanded / available (> 1: cannot be flown unpowered) |
| required altitude | `GM/v^2 (1/sin(turn/2) - 1) - R`: the periapsis altitude that would give exactly the demanded turn; below the floor (often negative, inside the body) when infeasible |
| impulse beyond bend | one impulse outside the sphere of influence after a maximal unpowered bend towards `v_out` in the plane of the two vectors: `|v_out - R(min(turn, bend)) v_in|`; equals `2 v sin((turn - bend)/2)` at equal magnitudes (Strange and Longuski 2002, eq. 9) and also pays any magnitude mismatch. A generalisation of `uranus_bridge_885.entry_correction_kms` to 3-D |
| periapsis impulse | a cheaper alternative, the Oberth-credited periapsis pair of `core.flyby.dv_powered_flyby_periapsis` (equal magnitudes only) |

Verdict: `turn_feasible` iff every demanded turn is no larger than its available bend;
`ballistic` iff, in addition, every magnitude mismatch is within 1e-3 km/s. Floors come from the
registries (`PLANETS[...].safe_alt_km`, `SATELLITES[...].safe_alt_km`: Earth 200 km, Uranian tour
moons 50 km, Galilean 100 km except Callisto 200 km, Saturnian midsize 100 km, Titan 1500 km).

Existing helpers, and why none of them was the gate:

* `core.flyby.max_bend`, `bend_angle`, `dv_from_turn_deficit`, `dv_powered_flyby_periapsis`: correct
  primitives, reused unchanged.
* `core.flyby.is_ballistic_feasible`: the right test but a bare boolean with a 1e-6 km/s magnitude
  tolerance, no diagnostics; it is not called by the moon-tour lane.
* `search.physical_sanity.candidate_passes_physical_gate` (the `#324` gate): CAPACITY only, i.e.
  that `max_bend` at the floor exceeds 5 degrees. It never sees a direction. It is the gate the
  `#558`/`#563`/`#571`/`#575`/`#576`/`#655` enumerations used.
* `search.correct._bend_feasible` and `search.turn_ratio_check.measure_turn_ratio` (`#833`): real
  demanded-turn checks for the heliocentric corrector, but on intermediate encounters only; the
  periodicity wrap is checked only when `closure_turn_ratio(include_wrap=True)` is asked for.
  `turn_gate.encounters_from_vinf_nodes` adapts their node dictionaries (with an optional wrap).
* `search.uranus_bridge_885`: planar, Uranus-only primitives; left unchanged.

## 2. Positive control 1: McConaghy, Longuski and Byrnes (2002)

Source: "Analysis of a Broad Class of Earth-Mars Cycler Trajectories", AIAA 2002-4420, filed in the
private paper corpus as
`mcconaghy-longuski-byrnes-2002-analysis-broad-class-earth-mars-cycler-trajectories-AIAA-2002-4420.pdf`.
Same ideal model as the moon enumerations (circular coplanar, Lambert arcs between successive
encounters, pp.2-5), same criterion: the Earth flyby must rotate V-infinity through the line-of-apsides
rotation `Delta Psi = (n/7) 360 deg` (eq. 4, p.5), compared with the turn of a 200 km flyby (p.5).
Table 4 (p.6) gives, for 19 nPr cyclers, the required and the maximum turn angle; footnote e marks
three as ballistic (6S7, 6S8, 6S9). The rebuild (`mcconaghy_npr_cycler`) identifies each arc by the
columns that do not involve the turn (aphelion and Earth V-infinity, all within 0.01) and then gates it.

| cycler | turn, ours / published (deg) | max at 200 km, ours / published | verdict, ours / published |
|---|---|---|---|
| 1L1 (Aldrin) | 83.69 / 84 | 71.81 / 72 | powered / powered |
| 2L2 | 134.11 / 134 | 43.99 / 44 | powered / powered |
| 2L3 | 135.15 / 135 | 81.87 / 82 | powered / powered |
| 3L4, 3L5, 3S5 | 166.69, 166.79, 167.28 / 167 | 35.39, 61.54, 33.37 / 35, 62, 33 | powered / powered |
| 4S5, 4S6 | 166.70, 166.77 / 167 | 37.89, 54.19 / 38, 54 | powered / powered |
| 5S4 ... 5S8 | 133.98, 134.34, 134.71, 135.10, 135.52 / 134, 134, 135, 135, 136 | 40.89, 50.09, 62.41, 79.31, 102.92 / 41, 50, 62, 79, 103 | powered / powered |
| 6S4, 6S5, 6S6 | 82.89, 83.45, 84.03 / 83, 84, 84 | 58.74, 67.73, 78.17 / 59, 68, 78 | powered / powered |
| **6S7** | 84.61 / 85 | 90.29 / 90 | **ballistic / ballistic** |
| **6S8** | 85.21 / 85 | 104.31 / 104 | **ballistic / ballistic** |
| **6S9** | 85.82 / 86 | 120.38 / 120 | **ballistic / ballistic** |

All 38 angles within 0.6 deg of the integer-degree table (test tolerance 1 deg, set before running),
19/19 verdicts agree. The Aldrin deficiency is also reproduced as a NEGATIVE altitude: McConaghy,
Longuski and Byrnes, JSR 41(4) 2004, p.627 (filed as
`mcconaghy-longuski-byrnes-2004-analysis-class-earth-mars-cycler-trajectories-jsr-doi-10.2514-1.11939.pdf`):
"the required flyby altitude is -1731 km"; the gate gives -1723 km (test tolerance 100 km, from the
~135 km/deg sensitivity times the 0.5 deg table rounding). This control is stronger than a bare pass:
the gate passes the three published-ballistic cyclers and fails the published-powered ones by the
published margins. The impulse estimates for 1L1 (1.35 km/s outside the SOI, 1.00 km/s with the
periapsis model) are our numbers and are NOT checked against anything: McConaghy defers the
Aldrin delta-v and the project already bans using it as a golden.

## 3. Positive control 2: Russell and Strange (2009) moon cyclers

Source: "Cycler Trajectories in Planetary Moon Systems", JGCD 32(1) 2009, filed as
`russell-strange-2009-cycler-trajectories-planetary-moon-systems-JGCD-32-doi-10.2514-1.36610.pdf`.
Every Table 3/4 cycler built only from generic free-return legs (no resonant `f`/`h` legs, whose crank
is a free parameter) was rebuilt in the paper's own ideal model: Table 2 constants (p.148), legs from
the Tables 5-6 nomenclature (p.151) `g<N>;<theta>;<label>` read as a Lambert arc of `N` flyby-body
periods with spacecraft transfer angle `theta` (revolutions `floor(theta/360)`). Turns are measured
only at the flyby body (the target body is massless in their ideal model), between consecutive legs and
across the wrap, in the body's local frame (their ground tracks repeat each cycle, p.149). The Lambert
branch was chosen by the PUBLISHED V-infinity, and confirmed against the published period and the
published minimum and maximum distance to the primary (all within 0.1 %); the published altitude was
never used to choose. Golden: Table 3/4 "Min flyby alt. at body A".

| cycler | V-inf published / rebuilt (km/s) | min flyby altitude published / gate's required altitude (km) | demanded / available at source floor (deg) |
|---|---|---|---|
| GanCal#5 | 3.24 / 3.238 | 328 / 328 | 27.95 / 30.57 |
| GanEur#5 | 1.66 / 1.658 | 1819 / 1819 | 53.09 / 70.53 |
| GanEur#43 | 1.87 / 1.872 | 8861 / 8860 | 22.74 / 62.29 |
| GanIo#53 | 3.90 / 3.900 | 518 / 518 | 19.69 / 22.83 |
| GanIo#403 | 4.29 / 4.281, 4.292 | 540 / 557 | 16.63 / 19.57 |
| TitEnc#37 | 2.50 / 2.505 | 1377 / 1377 | 30.83 / 33.22 |
| TitEnc#314 | 3.59 / 3.589 | 1218 / 1219 | 17.86 / 18.78 |
| TitEnc#510 | 5.03 / 5.026 | 1852 / 1852 | 8.52 / 10.38 |
| TitEnc#552 | 5.16 / 5.161 | 3784 / 3784 | 5.77 / 9.89 |
| TitEnc#586 | 5.43 / 5.429 | 3874 / 3875 | 5.17 / 9.01 |

Nine of ten reproduce the published minimum altitude to within 1 km. GanIo#403 is the exception
(557 against 540 km): its first leg is labelled `Ll`, both of its branches are within 0.06 km/s of the
published V-infinity, and the two legs' V-infinities differ by 0.011 km/s; the test holds it to 25 km
and says so. All ten pass the gate at the source floor (1000 km at Titan, p.144; the surface at the
Galilean moons, where every published altitude is positive).

**Floor conflict (flagged, nothing changed):** at the project's Titan floor of 1500 km, TitEnc#37
(1377 km) and TitEnc#314 (1218 km) would FAIL. Russell and Strange used 1000 km. The project floor is
a stricter convention, not a property of the trajectory; anyone using the gate on Titan cyclers must
choose the floor deliberately.

## 4. Regression: the six withdrawn Uranian rows

Rows read from `data/withdrawn/*-1-1-uranian-quasi-cycler-2026.yaml`; the arcs rebuilt by
`symmetric_closure` with rel_offset and n_rev chosen by reproducing the row's own four-decimal
V-infinities (reproduced to 5e-5 km/s). Two data defects in the withdrawn rows were found on the way:
the rows carry no rel_offset field (only prose notes), and the `#312` row's `legs` block lists
`n_revs: 0` while the trajectory (and its name, and `#885`) is (1, 1).

| row | first flyby: demanded / available (ratio), required altitude | closing flyby | impulse beyond bend per lap (km/s) |
|---|---|---|---|
| ariel-oberon | Oberon 122.8 / 8.1 (15.2x), -753 km | Ariel 36.1 / 6.2 (5.8x), -499 km | 3.08 + 0.78 |
| ariel-titania | Titania 131.8 / 9.6 (13.7x), -782 km | Ariel 17.9 / 9.2 (1.9x), -280 km | 3.01 + 0.19 |
| ariel-umbriel | Umbriel 155.7 / 8.4 (18.5x), -584 km | Ariel 102.7 / 14.0 (7.4x), -554 km | 2.50 + 1.37 |
| titania-oberon | Oberon 152.6 / 7.0 (21.7x), -760 km | Titania 177.2 / 6.3 (28.3x), -789 km | 3.76 + 4.31 |
| umbriel-oberon (#312) | Oberon 97.0 / 24.9 (3.9x), -687 km | Umbriel 45.3 / 16.5 (2.7x), -415 km | 1.13 + 0.46 |
| umbriel-titania | Titania 116.8 / 24.4 (4.8x), -750 km | Umbriel 137.5 / 9.3 (14.7x), -581 km | 1.45 + 2.21 |

Every encounter fails, at the 50 km floor, with ratios 1.9 to 28; every required periapsis is inside
the moon. This agrees with `#885` section 4 to the printed digits. Section 10's 39.5 / 15.6 deg at
Umbriel for #312 came from assuming the closure repeats exactly in the rotating frame; its leg time is
not exactly a half-integer number of synodic periods, and the directly solved closing leg gives
45.3 / 16.5 (the gate's local-frame cross-check gives 40.1). The test asserts only the verdict, that
every ratio exceeds 1.5 and every required altitude is negative.

## 5. Re-screen of the stored enumerations

Every stored gate-passing record was rebuilt (`symmetric_closure`, the `residual_at_point`
construction: circular-coplanar moons from `discovery_campaign._moon_state`, minimum `|v1 - v_moon|`
branch, `max_revs = max(n_rev, 1)`; the third leg solved directly at `2 tof`) and its magnitudes
checked against the stored `vinf_per_encounter_kms` (worst error 1.0e-13 km/s across all sets) before
the directions were used. A record and its mirror (`A-B-A` and `B-A-B` at the same leg time) are one
periodic chain; "distinct" counts chains.

| set | records (distinct chains) | pass at project floor | pass at 50 km | smallest worst-encounter ratio |
|---|---|---|---|---|
| `#563` Uranus (Ariel, Umbriel, Titania, Oberon) | 60 (30) | 0 | 0 | 3.46 |
| `#576` Jupiter (Galilean) | 36 (18) | 0 | 0 | 1.57 |
| `#575` Saturn Titan-Iapetus | 18 (9) | 0 | 0 | 10.2 |
| `#655` Saturn Rhea-Titan | 6 (3) | 0 | 0 | 7.07 |
| `#655` Dione-Rhea, Enceladus-Tethys, Tethys-Dione; `#599` Triton-Proteus; `#609` Phobos-Deimos | 0 stored | - | - | - |

The Saturnian symmetric sets are `#575` and `#655`, not `#571`: `#571` is a free grid scan
(`data/scan_571_*.jsonl`, non-symmetric offsets such as 228 deg and non-commensurate leg times), whose
closing encounter is not a periodic flyby; it was not re-screened. All 66 stored symmetric closures fail.

### 5.1 Without the capacity gate, same ranges (`ungated.jsonl`)

The capacity gate pruned exactly the class most likely to pass (a near-tangent encounter demands
almost no turn whatever the moon's mass), so each enumeration was re-run over its own `_meta` ranges
(n_rev 0-3 per leg, rel_offset 0/180 deg, leg time `n T_syn/2` up to the file's `tof_scale_max`:
3.0, except 2.0 at Jupiter) with only the 0.05 km/s closure gate and the turn gate; Miranda added at
Uranus. Result: 852 residual-gate closures (422 distinct), **0 pass** at either floor. Smallest
worst-encounter ratios: Uranus 3.46, Jupiter 1.57, Titan-Iapetus 10.2, Rhea-Titan 7.07, Dione-Rhea
42.7, Tethys-Dione 45.8, Enceladus-Tethys 197, Phobos-Deimos about 26,000.

### 5.2 Extended ranges (`extended_<primary>.jsonl`)

Same ideal model and construction, widened: n_rev 0-6 per leg, BOTH multi-revolution Lambert
branches per leg (the original rule took only the minimum `|v1 - v_moon|` one), leg time up to
`tof_scale` 6 (i.e. up to 6 sqrt(P_A P_B) per leg), rel_offset 0/180 deg. Neptune is excluded (see
section 7).

| system | residual-gate closures (distinct) | pass at floor | pass at 50 km | smallest worst ratio |
|---|---|---|---|---|
| Uranus (Miranda, Ariel, Umbriel, Titania, Oberon) | 3284 (1641) | **1 chain** | **1 chain** | **0.572** |
| Jupiter (Galilean) | 2084 (1042) | 0 | 0 | 1.28 |
| Saturn Titan-Iapetus | 468 (234) | 0 | 0 | 1.83 |
| Saturn Rhea-Titan / Dione-Rhea / Enceladus-Tethys / Tethys-Dione | 392 / 180 / 84 / 112 | 0 | 0 | 5.08 / 7.46 / 96.4 / 45.8 |
| Mars Phobos-Deimos | 384 (192) | 0 | 0 | 8777 |

**The one passing chain** (`candidates.jsonl`): Titania-Oberon-Titania, equal legs of
61.580 d (= 5 T_syn/2; cycle 123.16 d = 5 synodic periods), 5 revolutions per leg on the
higher-energy ("high") multi-revolution branch, rel_offset 0 deg (Oberon-Titania-Oberon at 180 deg is
the same chain). Both legs lie on congruent ellipses with a = 511,105 km, e = 0.1485, periapsis
435,210 km (just inside Titania's orbit), apoapsis 587,001 km (just outside Oberon's): a near-Hohmann
orbit with period 11.04 d, tangent-ish to both moon orbits, which is why the demanded turns are small.

| encounter | V-inf in / out (km/s) | demanded / available at 50 km (deg) | ratio | required altitude |
|---|---|---|---|---|
| Oberon (middle) | 0.2689 / 0.2689 | 58.44 / 102.11 | 0.57 | 2215 km |
| Titania (closing) | 0.2747 / 0.2747 | 42.16 / 102.86 | 0.41 | 4562 km |

Closure residual 7e-16 km/s; magnitude mismatches below 1e-14 km/s. Inside each leg the conic's
closest re-approach to Titania is 129,600 km (12.6 Hill radii) and to Oberon 337,800 km; the orbit
crosses no other regular moon's orbit (Umbriel is at 266,000 km).

What this is and is not (inference, not computed): it passes the gate in the ideal circular-coplanar
patched-conic model, the model the withdrawn rows also lived in. It is NOT yet a trajectory. Its
V-infinity is only 0.075 (Titania) and 0.085 (Oberon) of the moon's orbital speed, lower than any
published control above (Russell-Strange's lowest is GanEur#5 at 0.15), and its flyby periapses sit at
about half a Hill radius (Titania r_p 5,350 km, Hill radius 10,270 km). In that regime the
instantaneous-flyby patched conic is a poor model; the next check is a full Uranus-Titania-Oberon
point-mass (or CR3BP-patched) propagation of one cycle from these states, then the literature check
(near-Hohmann Titania-Oberon transfers are the obvious low-energy tour building block, so novelty is
doubtful). A concrete first step: `src/cyclerfinder/verify/flyby_integrate.py` can integrate one flyby
in Uranus plus moon from several Hill radii in to several out and compare the exit direction with the
patched-conic prediction. No catalogue row; nothing written to `empty_regions.jsonl`. The same
pipeline returning a pass here also shows that the zero-pass systems are not a stuck filter.

The ranges above are the whole claim: within symmetric (rel_offset 0/180, equal-leg, `n T_syn/2`)
two-moon closures, n_rev 0-6, both branches, legs up to 6 sqrt(P_A P_B), no other turn-feasible closure
exists in the six systems at their project floors or at 50 km. Asymmetric closures, unequal legs,
three-moon sequences (`#600`) and the small-body set (`#607`) were not re-screened.

## 6. Where the gate belongs in the validation lane (not wired in this task)

1. **V2 moon tour**, `src/cyclerfinder/data/validation/v2_moontour.py::_cycle_residual` (line 240).
   The function already holds `best.v1` and `best.v2` of every leg; keep the vectors (`best.v1 - v_a`,
   `best.v2 - v_b_moon`) instead of only their norms, build an `Encounter` at every intermediate
   encounter (arrival of leg k against departure of leg k+1) and at the anchor wrap (arrival at the last
   encounter against the next cycle's first departure: either solve cycle k+1's first leg, which the
   multi-cycle loop in `run_v2_moontour` at line 464 already does, or compare in `to_body_local`
   coordinates when the cycle repeats exactly), call `demanded_turn_gate`, and make `run_v2_moontour`'s
   pass require `turn_feasible` in addition to the magnitude residual. This is the single change that
   would have stopped the six rows at V2.
2. **V4 / V4-strict**, `v4_uranus.py::_cycle_v4` (leg loop near line 555, "IC for the V4 leg ... using
   Lambert's v-out") and `v4_uranus_strict.py::_cycle_v4_strict` (line 549, via `_select_leg_transfer`
   at line 437), and the Saturn twins `v4_saturn.py::_cycle_v4` (line 235) and
   `v4_saturn_strict.py::_cycle_v4_strict` (line 346): every leg restarts from its own Lambert
   departure velocity, so a flyby is never flown. Either propagate across the flyby (rotate the arrival
   V-infinity by at most the available bend and carry the result) or, at minimum, apply the gate to
   (arrival velocity of leg k, departure velocity of leg k+1) at each encounter and fail the tier when
   it fails.
3. **Discovery**, `scripts/scan_558_uranus_all_pairs_offset_sweep.py::gate_candidate` and its
   generic callers (`#563`, `#575`, `#576`, `#599`, `#609`, `#655`): replace the `#324` capacity call
   (`physical_sanity.candidate_passes_physical_gate`) with `symmetric_closure` +
   `demanded_turn_gate`, as `scripts/screen_888_turn_gate_rescreen.py` does.
4. **Heliocentric corrector**, `search/correct.py::_bend_feasible` (line 349): checks intermediate
   encounters only; for any multi-lap claim add the wrap (`turn_gate.encounters_from_vinf_nodes` with
   the rotation of `turn_ratio_check.wrap_node_turn`).

## 7. Other findings and open items

* `#599` (Neptune): `discovery_campaign._moon_state` returns the prograde circular velocity for every
  moon, while the `#599` construction advances retrograde Triton with a negative mean motion. Triton's
  position therefore moves clockwise while its velocity points counter-clockwise, so its V-infinity
  vectors (and magnitudes) are not physical. Checked by a run: a one-minute finite difference of
  Triton's position is antiparallel to the returned velocity (cosine -0.99999993). `#599` stored no
  passes; the extended run here skipped Neptune rather than reproduce the defect.
* The frozen-gate tests `tests/verify/test_566_five_representatives_v4.py` and `test_silver_327_*`
  test the insufficient gauntlet (item 4 of the `#888` ledger) and are untouched.
* `data/found/650_transfer_network/` still carries edges between the withdrawn rows (untouched).
* Wiring must gate on `turn_feasible` plus the lane's own magnitude residual, not on `ballistic` with
  the 1e-3 km/s default: GanIo#403 (published ballistic) has a 0.011 km/s leg-to-leg V-infinity spread
  from the rounding of its printed leg parameters, so `ballistic` is False there while
  `turn_feasible` is True. `mag_tol_kms` must match the precision of the source.
* Earlier negatives now rest on a different ground. The symmetric-closure empty regions stamped under
  the capacity gate (`#575`, `#599`, `#609`, `#655`) still hold within the ranges of section 5, but
  because of the demanded turn (smallest worst ratios in 5.1 and 5.2), not capacity. The `#576` set
  (36 records, 18 chains) is not flyable at all; `#577` read it as Russell-Strange class members, but in
  Russell-Strange the target body is massless and passed on the same conic, whereas in these closures
  both moons must turn, so that framing does not hold.
* Plan item 2(b) (gating the catalogue's heliocentric rows) was not run: only the adapter
  `encounters_from_vinf_nodes` was built and tested on synthetic nodes. The `#833`
  `turn_ratio_check.closure_turn_ratio` already measures demanded turns on that lane against
  published turn ratios; extending it with the wrap is item 4 of section 6.
* `uranus_bridge_885.py` was left unchanged; its `entry_correction_kms` is the planar special case of
  `turn_gate.impulse_beyond_bend_kms`.
