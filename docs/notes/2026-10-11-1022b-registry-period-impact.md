# #1022b: what the registry's slow small-moon mean motions touch, and the two options (2026-10-11)

Read-only. No code, registry or catalogue change. Context: `docs/notes/2026-10-11-1022-pluto-small-moon-registry-verification.md`
sec. 3-4. The registry stores the JPL page's printed a for Styx, Nix, Kerberos and Hydra and derives
`mean_motion_deg_day` from it by Kepler III about the system GM 975.5 km^3/s^2. The kernel-fitted
(plu060, 2030-2033) periods are shorter, so the derived mean motions are slow.

## 1. Size of the error

| Moon | Registry-derived period (d) | Fitted period (d) | Registry mean motion against fitted |
|---|---|---|---|
| Charon | 6.3891 | 6.38722 | -0.03 percent |
| Styx | 20.9063 | 20.16195 | -3.56 percent |
| Nix | 25.4872 | 24.85472 | -2.48 percent |
| Kerberos | 32.7760 | 32.16798 | -1.85 percent |
| Hydra | 38.7636 | 38.20202 | -1.45 percent |

(The lead's range 1.5-3.7 percent is these four, rounded.) Phase error a moon accumulates in its own
orbit when advanced at the registry rate, by horizon:

| Moon | 25 d | 50 d | 100 d |
|---|---|---|---|
| Charon | 0.4 deg | 0.8 | 1.6 |
| Styx | 15.9 deg | 31.8 | 63.6 |
| Nix | 9.0 deg | 18.0 | 35.9 |
| Kerberos | 5.2 deg | 10.4 | 20.8 |
| Hydra | 3.4 deg | 6.8 | 13.7 |

What matters to a cycler search is the pair's synodic period, because a cycle is closed over whole
synodic periods:

| Pair | Registry synodic (d) | Fitted synodic (d) | Registry error | Extra relative phase per registry synodic period |
|---|---|---|---|---|
| Charon-Styx | 9.20 | 9.35 | -1.6 percent | 5.8 deg |
| Charon-Nix | 8.53 | 8.60 | -0.8 percent | 2.9 deg |
| Charon-Kerberos | 7.94 | 7.97 | -0.4 percent | 1.4 deg |
| Charon-Hydra | 7.65 | 7.67 | -0.3 percent | 1.1 deg |
| Nix-Hydra | 74.42 | 71.14 | +4.6 percent | 16.6 deg |
| Styx-Nix | 116.32 | 106.79 | +8.9 percent | 32 deg |
| Kerberos-Hydra | 212.19 | 203.66 | +4.2 percent | 15 deg |

Charon pairs are mild (the error is the moon's, Charon is right); pairs of two small moons carry 4-9
percent.

## 2. Every path and past result that derives a small-moon rate or period from the registry a

Method: grep of src, scripts, tests and docs/notes for the accessors (`mean_motion_deg_day`,
`SATELLITES[...]`, `sma_km`, `_mean_motion_rad_day`, `moon_circular`) and for Styx, Nix, Kerberos,
Hydra; then a check of which data files carry a Pluto small-moon row.

Code that derives the rate or period from the registry a (the moon is any registry moon):
- `core/satellites.py` `mean_motion_deg_day_about` / `_sat` (the source of the derived value).
- `search/discovery_campaign.py` `_mean_motion_rad_day` and the campaign constants (`{moon: (sma_km,
  mean_motion)}` from the registry); `_registry_moons_for("Pluto")` enumerates every registered Pluto moon,
  so since #1037 a campaign would also include Styx and Kerberos.
- `search/saturn_uranus_campaign.py` `run_prioritized_scan` (used by `scripts/scan_492_neptune_pluto.py`).
- `search/five_tier_prioritizer.py` (about line 445): Kepler mean motion from the registry sma for the
  feasibility re-check of a candidate.
- `scripts/scan_320_epoch_aware_moon_systems.py` (imports `_mean_motion_rad_day` and the registry sma).
- `search/two_working_body.py` `moon_circular`: radius from the registry sma, period from Kepler III. It is
  not used for any Pluto cell: the #998 driver builds its own circles from the fitted periods.

Not affected: the #998 enumeration, gauntlet and real-ephemeris columns (circles from fitted periods; the
SOI now uses the system's own circle radius, #1060; Charon constants agree with the registry to 0.02
percent); the Pluto-Charon CR3BP/ER3BP work (`pluto_charon_kk_sweep`, `real_binary_kk_sweep`,
`pluto_charon_realeph`, the catalogue's two Pluto-Charon rows) which uses Charon only; `turn_gate`
`tidal_speed_kms` (a massless target is not gated; Charon is correct). `data/er3bp_direct_436.jsonl`,
`scan_334`, `tulip_sweep_281` have Pluto rows with no small moon.

Past results that carry the registry rates:
- #320 Vector B, `data/scan_320_epoch_aware_pluto.jsonl`: 53 rows, all with a small moon (Charon, Hydra,
  Nix, 6 cycles, 90 cells). The two SILVERs (Hydra-Nix-Hydra (1,1) 0.00137 km/s, Nix-Hydra-Nix (1,1)
  0.00070 km/s) are two-leg cycles of 47.15 d per leg: over their 94.3 d the registry Nix is 33.9 deg
  and Hydra 12.9 deg behind, a relative error of about 21 deg; their Nix-Hydra synodic period is
  74.4 d against the true 71.1 d. Verdict already V0-known and "degenerate co-orbital closure, not
  operational" (#320 notes), so no conclusion rests on the phasing. The other near-miss rows (51 in all) use
  legs of 19-63 d; a Nix or Hydra phase carries 3-18 deg of the error above over one leg-pair of 25-50 d.
- #492 Pluto sweep, `data/scan_492_pluto.jsonl` (task label "#285 prioritized repeated-moon scan"): 9
  small-moon rows, legs 19-63 d (Hydra-Nix-Hydra 47.1 and 62.9 d). Same size of error as #320.
- #264/#282 campaign Pluto SILVERs, `data/silver_282.jsonl`, repeated in `data/gauntlet_queue.jsonl` and
  `data/review_queue.jsonl`: 12 rows, all with Nix and Hydra (two- and three-moon sequences, Charon-Hydra-
  Nix-Charon and similar). Their registry-model phasing is wrong by the pair errors above; no row was
  promoted on it (the Pluto verdicts are V0-known or not-novel in the #320/#492 notes).
- #998: one place only. `scripts/gauntlet_998.py` `collide_320` compares a #998 candidate's period (fitted
  rates) with the total flight time of a #320 row (registry rates); the 454 NEAR tags on Nix come from the
  51.04 d (registry) against 51.58 d (fitted) k = 6 period. The tag was already stated to carry no weight.
- No catalogue row depends on a Pluto small-moon rate.

## 3. The two options

Because the registry also has to feed Kepler-consistent code (a circle of radius a moved at the speed
sqrt(mu/a)), storing a fitted period beside the printed a cannot simply replace the derived mean motion:
position, speed and period of a circle would disagree. The workable form of (a) stores the fitted period and
derives a Kepler radius from it (what #998 does), keeping the printed JPL a as the sourced value.

(a) Store fitted periods beside a, and let period-based code prefer them.
- Add a `period_days` (fitted, source: plu060 fit, span and method stated) to the four small moons and a
  derived `kepler_sma_km`; `mean_motion_deg_day` for those moons follows the period; period-based code
  (campaign constants, five-tier prioritizer, scan_320, `moon_circular`) uses the fitted rate and the
  Kepler radius.
- Consequence: a rerun of #320, #492 or the #264/#282 campaign would give different cells (the Nix-Hydra
  synodic period changes by 4.6 percent, so resonance grids and `tof_resonance_grid` multiples move); the
  committed JSON stays as the record of the registry-Kepler model and is no longer bit-reproducible from
  the code. Tests that pin Pluto small-moon values (`tests/core/test_satellites_registry.py`,
  `test_satellites_pluto_small_moons.py`, `tests/search/test_five_tier_prioritizer.py`, anything frozen on
  `scan_320`) need review. A source-pinned test (fitted period against the page's printed period to the
  printed digits) is easy. The same offset may exist in other registry systems (Jovian, Saturnian, Uranian
  moons); not checked here, and an option (a) for the Pluto moons only leaves them inconsistent.
- Cost: a registry schema field, a migration of four consumers, a test review. Benefit: any future Pluto
  small-moon sweep is phased correctly without each driver rebuilding its own circles.

(b) Leave as is and document.
- Add a comment at the four entries (and in `data/README.md` if it lists the registry) saying the derived
  mean motion is 1.5-3.6 percent slow and that period-sensitive work must use the fitted periods (as #998 does).
- Consequence: the three past sweeps stay reproducible; their registry-model phasing stays as described in
  sec. 2 and is documented, not changed. Risk: the next registry-based Pluto sweep inherits the offset
  silently unless the author reads the comment; #1037's two new moons make the registry-based campaign larger,
  with the largest error (Styx, 3.6 percent; Styx-Nix 8.9 percent in synodic period).

Not decided here; owner decision, registered under #1022 in OUTSTANDING.

## Addendum 2026-10-11: the cause of the Pluto offset is sourced

Brozovic and Jacobson (2024), Tables 10 and 11 (digest `docs/notes/2026-10-11-digest-brozovic-jacobson-2024-post-new-horizons-orbits-masses-pluto-satellites-aj.md`): the registry's a is the MEAN OSCULATING a, and its
Kepler period reproduces the paper's osculating period (Table 10) to 0.05-0.1 percent; the printed periods are the
precessing-ellipse (true mean-motion) periods of Table 11. The offsets in this note are the osculating-versus-mean
difference, not a registry error. Where this note says the cause is "the binary quadrupole, remainder unexplained",
read: the binary forces the osculating elements away from the mean ones, and the paper's tables give both. See
`docs/notes/2026-10-11-1022-pluto-small-moon-registry-verification.md` sec. 5. The numbers and the options are
unchanged; for option (a) the sourced period is Table 11's.
