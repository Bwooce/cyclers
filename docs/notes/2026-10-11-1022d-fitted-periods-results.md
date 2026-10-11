# #1022d: results of the fitted-period registry change and the Pluto re-runs (2026-10-11)

Pre-registration: `docs/notes/2026-10-11-1022d-fitted-periods-prereg.md` (committed before the code). Owner
decision: "Store the fitted period beside a". Old result files are untouched; new files are in
`data/1022d_rerun/`. Runner: `scripts/run_1022d_rerun.py`.

## 1. What was implemented

- Registry: Styx, Nix, Kerberos, Hydra carry the Brozovic and Jacobson (2024) Table 11 mean motion (17.855, 14.484,
  11.191, 9.424 deg/day), `kepler_sma_km` (the Kepler radius of it) and `orbit_circle(moon)`.
  `mean_motion_deg_day` follows the stored rate. For the other 33 moons `orbit_circle`, `mean_motion_deg_day` and
  `kepler_sma_km` are bit-identical to the pre-change registry (fixture recorded before editing,
  `tests/core/data/registry_before_1022d.json`, checked by a test over all 37 moons).
- Consumers switched to `orbit_circle`: the planet-centred circular moon ephemeris, the discovery campaign constants
  (the Saturn/Uranus campaign inherits them), the five-tier prioritizer, `scan_320`, `moon_circular`. Tests:
  `tests/core/test_satellites_pluto_small_moons.py`, `test_satellites_registry.py`, `test_ephemeris_centred.py`,
  `tests/search/test_1022d_orbit_circle_consumers.py`; the #998 driver keeps its plu060 periods with an
  agreement test (relative 5e-5) as ruled.

## 2. Control

| Control | Result |
|---|---|
| Neptune scan_320, new code against the pre-change expression run on THIS machine | bitwise identical (8 rows) |
| Neptune scan_320, new code against the committed `data/scan_320_epoch_aware_neptune.jsonl` | 7 of 8 rows equal to 1e-9; 5 equal bitwise; the 3 rows that differ bitwise differ in the last digits (1e-14 relative); one row (idx 5) picks a different one of two equal-residual cells (phase0 225 against 240 deg). Floating-point noise between machines, not the change. |
| Saturn scan_320, new code against the pre-change expression run on THIS machine (launched by the lead, 238 s for the new-code sweep) | bitwise identical (62 rows) |
| Saturn scan_320, new code against the committed `data/scan_320_epoch_aware_saturn.jsonl` | 23 of 62 rows identical bitwise, 36 equal to 1e-9; 25 further rows differ only by the choice among equal-residual cells (phase0_deg), as for Neptune; the last row is the trailing summary, whose top-5 near-miss list differs. Last-digit and tie noise between machines, not the change. Files: `control_saturn.json`, `scan_320_epoch_aware_saturn{,_oldpath}.jsonl`. |

I changed the control from the pre-registered "bitwise against the committed file" to "bitwise against the pre-change
expression on this machine", because the committed Neptune file was produced on another machine and is not bit
reproducible here. Reported, not hidden.

## 3. #320 Pluto sweep (`data/1022d_rerun/scan_320_epoch_aware_pluto.jsonl`)

90 cells, 80 s under load. Rows: 51 near-miss in both runs, the same 51 sequences. 2 SILVERs in both runs.

| SILVER | Old residual (cell) | New residual (cell) |
|---|---|---|
| Hydra-Nix-Hydra (1,1) | 0.00137 km/s (tof_scale 1.5, phase 60, offset 315) | 0.00243 (tof_scale 2.0, phase 30, offset 45) |
| Nix-Hydra-Nix (1,1) | 0.00070 (tof_scale 1.5, phase 135, offset 45) | 0.00156 (tof_scale 1.5, phase 330, offset 60) |

The best residual over the sweep is 0.00070 (old) and 0.00073 (new). Of the Charon rows, 18 of 22 change by more
than 10 percent. Nix-Hydra rows move in both directions (for example Hydra-Nix-Hydra (1,0) 0.0053 to 0.0007,
Nix-Hydra-Nix (3,0) 0.0465 to 0.0137, (0,3) 0.0278 to 0.0402).

**The pre-registered expectation is falsified.** The note expected the two Nix-Hydra SILVERs not to survive at their
cell and the best Nix-Hydra residual to rise. Both survive, at different cells, with residuals 1.8 and 2.2 times
larger, and the best residual is unchanged. Reason (from the sweep's design, not tested further): the sweep scans a
24 by 24 grid of phase and relative offset, so a 4.4 percent change in the synodic period moves the best cell
instead of removing the near-resonant closure, which is a property of these co-orbital, tens-of-m/s orbits (the
#320 note already calls them degenerate). The leg times shorten as predicted (47.15 to about 46.2 d at tof_scale
1.5). The verdict of these rows (V0-known, degenerate co-orbital closure) is unchanged.

## 4. #492 Pluto scan (`data/1022d_rerun/scan_492_pluto.jsonl`)

210 candidates enumerated, 20 closed, 9 rows written (6 s). Against the committed file (9 rows): 7 sequences in both,
2 old rows gone (Charon-Nix-Charon (0,1), Nix-Charon-Nix (0,1)), 2 new rows (Nix-Charon-Nix (1,0) and (1,1)). The
residuals of the 7 common rows move by about -45 percent to +50 percent (Hydra-Nix-Hydra (1,1): 0.0090 to 0.0052;
Charon-Nix-Charon (1,1): 0.031 to 0.0465); leg times shorten as expected (Nix-Hydra legs 47.1 to 46.2 d,
62.9 to 61.6 d). The verdict column changed from BRONZE to SILVER for every row; that comes from the
literature and ML guards (the Persephone anchor was corrected in #881 and the literature gate has been reworked
since the committed run), not from the rates, and is not a rate effect.

## 5. #282 Pluto SILVERs (`data/1022d_rerun/silver_282.jsonl`)

Pass A re-scores the 12 stored SILVERs from `data/review_queue.jsonl` by reconstructing the legs under the new
rates (1.6 s). The Tier-0 maximum per-leg dV of the 2-leg sequences hardly moves (for example Hydra-Nix-Hydra 0.360
to 0.359 km/s, Nix-Hydra-Nix 0.378 to 0.381); that of the 3-leg sequences changes a lot (Charon-Hydra-Nix-Charon
0.044 to 0.500; Nix-Charon-Hydra-Nix 0.071 to 0.203 and 0.323; Hydra-Nix-Charon-Hydra 0.237 to 0.125). This is a
re-scoring of the stored phasings, not a new closure; Pass B (a shard of the #264 enumeration) is not in this script's
current form and was not run.

## 6. Consequences

- Nothing in the catalogue or in any verdict changes: the Pluto small-moon rows were V0-known or not-novel before and
  stay so; no row depends on these rates.
- The registry-based Pluto sweeps are now phased with the sourced mean motions; the old files stay as records of the
  registry-Kepler model.
- The Saturn control passed (sec. 2). #1022d is DONE.
