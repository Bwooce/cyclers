# #1022d: pre-registration, store the Table 11 mean motion beside a for the four Pluto small moons (2026-10-11)

Status: DRAFT, written before any code or data change. Owner decision on #1022 (2026-10-11, verbatim label):
"Store the fitted period beside a". Nothing is implemented or committed until the lead reports window 3 green.

## 1. What is stored, and where it comes from

Source: Brozovic and Jacobson (2024), AJ 167:256, Table 11 "Mean Precessing Ellipse Model Orbital Elements for the
Satellites of Pluto" (journal p.10); digest
`docs/notes/2026-10-11-digest-brozovic-jacobson-2024-post-new-horizons-orbits-masses-pluto-satellites-aj.md`
(section "Mean elements", Table 11 row); transcription `data/sources/brozovic-jacobson-2024-tables.yaml`
(table 11). The stored quantity is the printed mean motion n (deg/day), because it carries more digits than the
printed period:

| Moon | n (deg/day), Table 11 | 1-sigma (deg/day) | P = 360/n (d) | Table 11 printed P (d) | #998 plu060 fit (d) |
|---|---|---|---|---|---|
| Styx | 17.855 | 2e-5 | 20.1624 | 20.16 | 20.16195 |
| Nix | 14.484 | 4e-6 | 24.8548 | 24.85 | 24.85472 |
| Kerberos | 11.191 | 7e-6 | 32.1687 | 32.17 | 32.16798 |
| Hydra | 9.424 | 2e-6 | 38.2003 | 38.20 | 38.20202 |

(The tabulated n has 3-5 significant digits; the relative difference between 360/n and the #998 fit is 2e-5 to 5e-5.)
Charon is not changed: its registry Kepler period is 0.03 percent from Table 11's (n 56.362530), a difference
below any sweep's resolution.

Added beside the existing printed a (Table 10 mean osculating a, kept as is):
- `SatelliteData.mean_motion_source_deg_day: float | None = None` (the sourced n) with a source string field;
- `mean_motion_deg_day` follows the sourced n when present, and is derived from a (Kepler III) otherwise, exactly as
  today;
- `kepler_sma_km`: the Kepler radius of that mean motion about the registry primary GM,
  (mu_primary / n^2)^(1/3), and `sma_km` itself when no sourced n is stored (bit-identical for all other moons).
  For the four moons: Styx 42168, Nix 48481, Kerberos 57577, Hydra 64569 km (0.6, 0.4, 0.3 and 0.2 percent below the Table 11
  a: 42410, 48689, 57749, 64719 km, which the paper's own consistency shows is not the Kepler radius of n).
- one accessor, `orbit_circle(moon) -> (radius_km, mean_motion_rad_day)`: `(kepler_sma_km, n)` where n is stored,
  and exactly `(sma_km, sqrt(mu / sma^3) * 86400)` otherwise (the expression the consumers compute today).

## 2. Consumers (audited by grep of `SATELLITES`, `sma_km`, `mean_motion`)

Change to use `orbit_circle` (they build a moon circle from a and a rate):
1. `core/satellites.py` `_sat` and the dataclass (the field and the accessor).
2. `core/ephemeris.py` moon backend (line 215-216: radius `sma_km`, rate `mean_motion_deg_day`): radius becomes
   `kepler_sma_km`, so the circle keeps `v = n a` consistent with Kepler.
3. `search/discovery_campaign.py` `_moon_consts` (and `_mean_motion_rad_day` stays the Kepler helper for unsourced
   moons); `search/saturn_uranus_campaign.py` takes its constants through `discovery_campaign`.
4. `search/five_tier_prioritizer.py` (about line 445, mean motion from `sqrt(mu / sma^3)`).
5. `scripts/scan_320_epoch_aware_moon_systems.py` (`consts[m] = (sat.sma_km, _mean_motion_rad_day(mu, sat.sma_km))`).
6. `search/two_working_body.py` `moon_circular` (radius from a, Kepler period).
7. `verify/turn_gate.py` `tidal_speed_kms` reads `s.mean_motion_deg_day` and `s.sma_km`: it follows the field
   automatically; it is only called for flyby bodies (Charon, unchanged), so no value changes.

Do not change (they use a as an orbital element, for which the osculating mean a is the right quantity, or as a
length scale): `search/tisserand.py`, `search/moon_prune.py`, `search/vilm.py`, `search/releg_moontour.py`
(checked: uses a for Tisserand-type screens), the CR3BP length scales (`core/cr3bp.py`, `ccr4bp*.py`, none is at
Pluto's small moons), `search/two_working_body.py:1243` (SOI fallback; the system's own constants come first since
#1060), validation modules for Saturn and Uranus.

#998 driver (`scripts/run_998_enumerate.py`): it keeps its plu060-fitted `PERIODS_D` (lead ruling 2026-10-11). The two
values answer different questions and both are kept:
- the registry's stored n is the SOURCED value (Brozovic and Jacobson 2024 Table 11, 3-5 significant digits, citeable
  in a catalogue or a note);
- the kernel fit is the #998 model's OWN input: more precise (it is a fit of the plu060 states over 2030-2033, the
  same kernel the real-ephemeris step uses), and changing it would move every #998 result by 2e-5 to 5e-5 in
  relative terms for no gain.
A test pins that the driver's periods agree with the registry's stored n to the printed digits (relative 5e-5), so
neither can drift from the other unseen.

## 3. Tests that pin the old values, and what each becomes

Expected values come from the source (the Table 11 n above, transcribed in the test), never from the code under test.
1. `tests/core/test_satellites_registry.py` line ~70: for every moon, `mean_motion_deg_day` equals Kepler III from
   `sma_km` to 1e-12. Becomes: the same for every moon except the four; for the four, equals the Table 11 n literal.
2. `tests/core/test_satellites_pluto_small_moons.py`:
   - `test_registry_period_within_documented_offset_of_printed_period` (Kepler-from-a within 1-4.5 percent of the
     printed period): becomes two tests, (i) the registry mean-motion period equals 360/n of Table 11 and agrees
     with the page's printed period to the printed digits, (ii) the Kepler period from the printed a stays 1-4.5
     percent above it (this documents osculating against mean, and is kept as a regression on `sma_km`);
   - `test_registry_matches_jpl_pages` (GM, radius, a): unchanged;
   - `test_task_998_driver_uses_registry_values`: unchanged, plus the new check that `PERIODS_D` agree with the stored
     n to 5e-5 relative.
   - New: `kepler_sma_km` equals (mu / n^2)^(1/3) of the literal n; `orbit_circle` is bit-identical to the old
     expression for every other moon in the registry (a loop over all 37 minus 4).
3. `tests/core/test_satellites_type.py`: unchanged (builds its own SatelliteData with no sourced n).
4. `tests/search/test_five_tier_prioritizer.py` (Pluto Hydra-Nix-Hydra candidate, lines 141-162): reconstructs
   legs from a stored candidate; it asserts leg labels, not rates. Expected unchanged; confirm by running it.
5. `tests/search/test_discovery_campaign.py`: no Pluto small-moon values found by grep; expected unchanged; run.
6. `tests/search/test_1060_soi_system_constants.py`, `test_1024_barycentric_binary.py`,
   `tests/scripts/test_1064_*`: build their own systems from literal numbers; unchanged.
7. `tests/core/test_satellites_uranian_sources.py`: Uranian, unchanged.
Regression gate for the refactor: a bitwise comparison of `orbit_circle` and `mean_motion_deg_day` for every
registered moon other than the four against the values before the change (recorded before editing).

## 4. Re-run plan for the past Pluto results (old files stay as records)

All three used registry-Kepler rates for Nix and Hydra. After the change the registry gives Nix and Hydra mean
motions 2.5 and 1.5 percent faster, and Charon is unchanged (Charon is not in the four).
- Nix-Hydra synodic period: 74.42 d to 71.14 d, i.e. 4.4 percent shorter (the old value is 4.6 percent longer
  than the new one). Charon-Nix 8.53 to 8.60 d (+0.8 percent), Charon-Hydra 7.65 to 7.67 (+0.3 percent). Styx and Kerberos
  enter any registry-derived Pluto campaign from now on (since #1037) with Charon-Styx 9.20 to 9.35 d.
- **#320 Vector B Pluto** (`scripts/scan_320_epoch_aware_moon_systems.py`, Pluto target only, output to a new
  file `data/1022d_rerun/scan_320_epoch_aware_pluto.jsonl`; the committed `data/scan_320_epoch_aware_pluto.jsonl`
  stays): 90 cells; runtime 27 s measured in the committed file's own summary (about 50 s under load). Expected
  direction: the leg times tof_scale*sqrt(Pa*Pb) for the Hydra-Nix pair shorten by about 2 percent (47.15 d to
  about 46.2 d for tof_scale 1.5), the two Nix-Hydra SILVERs (0.00137 and 0.00070 km/s) sit at a near-resonant
  phasing that a 4.4 percent change in the synodic period moves by about 16.6 degrees per synodic period; they are
  expected NOT to survive as SILVERs at the same cell, and the best Nix-Hydra residual is expected to rise, with a
  different cell possibly closing. The Charon rows move less (a 0.3 to 0.8 percent change in the pair synodic
  period). No outcome is pre-registered beyond that direction and the control below.
- **Control for the refactor (run first):** the same scan for Saturn and Neptune (not touched by the change) must
  be bitwise identical to `data/scan_320_epoch_aware_saturn.jsonl` and `..._neptune.jsonl` apart from `git_sha` and
  `elapsed_s`. Runtime: Neptune seconds; Saturn about 3 to 4 minutes (about 3.3 times the Pluto run); both within
  10 minutes per call, run as two calls.
- **#492 Pluto** (`scripts/scan_492_neptune_pluto.py`, Pluto target only, new output file): 10 rows committed (9
  with a small moon). Runtime: not recorded in the committed file; a 5-structure pilot is run first and the
  full run is handed to the lead if the estimate exceeds 10 minutes. Expected direction as for #320.
- **#264/#282 Pluto SILVERs** (`scripts/rerun_264_with_five_tier_282.py`, `data/silver_282.jsonl`, the 12 rows of
  `data/review_queue.jsonl`): Pass A (re-scoring the 12 stored SILVERs) is cheap; Pass B is one shard of the #264
  Pluto enumeration; runtime unknown, pilot first, hand any run over 10 minutes to the lead. Expected direction: the
  12 stored SILVERs were produced with the old rates, so under the new rates some no longer close at the stored phasing;
  Pass A scores the stored phasing as given and will show the per-leg changes, it does not re-close.
- A result file is written for each rerun with the old file's name as `rerun_of`, and a short comparison table (rows
  that survive, rows that change by more than the gate, new rows) in the note that reports it. Nothing is promoted
  and no catalogue row changes; the verdicts of these rows (V0-known, not novel) do not depend on the rates.

## 5. Order of work (after "green")

1. Record before-values for every moon (bitwise reference file in the tests).
2. Registry change + accessor + tests (section 3), ruff, mypy src tests, tests/core and tests/search files above.
3. Consumers 2-6 one by one, each with its test.
4. Reruns (section 4), control first.
5. Update `docs/notes/2026-10-11-1022-*` notes and the OUTSTANDING #1022 line.
Each step is its own pathspec commit with ruff, mypy and the relevant tests green first.

## 6. Rulings and what stays out (lead, 2026-10-11)

- The #998 driver keeps `PERIODS_D` plus the agreement test (sec. 2).
- The Styx and Kerberos registry GM (Table 8 upper bounds, 7.5 and 3.3 times the values held in PLU060) stay OUT of
  this change. Registered separately as #1067 (BACKLOG, small): candidate fixes (a) a "held value used by PLU060"
  field beside the bound, (b) keep the bound with a comment at the entry; consequence: bending overstated if either
  moon is ever a flyby body; no action now.
- Order after "green": this pre-registration first, then the code and tests, then the control and the reruns, in that
  order; anything over 10 minutes is handed to the lead.
