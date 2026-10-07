# #998: Pluto-Charon one-working-node cyclers with passive small moons (2026-10-07)

Status: IN PROGRESS. Sections 1-2 written before the enumeration; results are added below as each
moon finishes. The literature step is DEFERRED (the gate is being reworked under #972). Nothing here
is called novel.

Driver: `scripts/run_998_enumerate.py`. Data: `data/998_pluto_smallmoons/`.

## 1. Model

Circular, coplanar, ideal model. Primary mass: Pluto+Charon system GM 975.5 km^3/s^2 (registry
`PRIMARIES["Pluto"]`). Each body rides a circle whose radius follows from its period by Kepler III
(the Russell-Strange convention). Charon is the only massive flyby body (GM 106.1, R 606 km, floor
100 km, registry). The small moon is a massless target.

Periods are fitted to `plu060.bsp` (2030-01-01 plus 3 yr, 3000 samples, J2000, relative to the
Pluto-system barycentre NAIF 9), not taken from the registry:

| Body | Fitted period (d) | Ratio to Charon (fit) | Kepler a (km) | Registry a (km) | Registry ratio |
|---|---|---|---|---|---|
| Charon | 6.38722 | 1 | 19596 | 19600 | 1.000 |
| Styx | 20.16195 | 3.1566 | 42168 | not in registry | - |
| Nix | 24.85472 | 3.8913 | 48481 | 49300 | 3.989 |
| Kerberos | 32.16798 | 5.0363 | 57577 | not in registry | - |
| Hydra | 38.20202 | 5.9810 | 64569 | 65200 | 6.067 |

The brief's ratios (3.16, 3.89, 5.04, 5.98) agree with the fit. The registry Nix and Hydra values
are 1.7 and 1.0 percent off in a (papercut filed). The registry is not edited.

The four small moons are coplanar with Charon's orbit: the angle between each moon's mean orbit
normal and Charon's is 0.006 (Styx), 0.010 (Nix), 0.348 (Kerberos), 0.260 (Hydra) degrees (same
fit). The coplanar model is therefore a good fit to the geometry.

Known model offset: Charon is at 19596 km in the ideal model (Kepler a of the relative orbit) but
at a mean 17464 km from the barycentre in the kernel. Its speed in the ideal model is therefore
about 12 percent above its barycentric speed. Wobble of Pluto is ignored. The real-ephemeris step
(sec. 6) is where this shows.

## 2. Controls

Two controls. Neither is a Charon-flyby control; the #320 result cannot be recalled by the #998
generator (reasons below). The lead accepted this on 2026-10-07.

1. #320 self-regression. `run_998_enumerate.py control320` re-runs #320's own Pluto sweep through
   its own code (`scan_320_epoch_aware_moon_systems._per_system_sweep`; the module-level import of
   `literature_check` is satisfied by an inert stand-in, so the anchor-overlap field is not
   reproduced). 90 cells evaluated. All 51 committed near-miss rows are reproduced: residual and
   per-encounter V_inf agree to 4.3e-14 km/s, the physical-gate flag and SILVER flag agree on every
   row, no extra rows. The two SILVERs recur: Hydra-Nix-Hydra (1,1) residual 0.00137 km/s and
   Nix-Hydra-Nix (1,1) 0.00070 km/s. Data: `control_320_pluto_rerun.jsonl`,
   `control_320_pluto_rerun.compare.json`.
   Why this is a self-regression only and not a recall by the #998 pipeline:
   - neither SILVER contains Charon; Hydra and Nix are both flyby bodies with GM about 0.002 at
     V_inf of 15-30 m/s;
   - #320 closes V_inf magnitudes plus an anchor wrap, while the generator needs an exact periodic
     closure over k synodic periods;
   - the two legs sum to 94.3 d, and the registry Nix-Hydra synodic period is 74.4 d (ratio 1.267,
     not an integer).
2. Generator code path. `run_998_enumerate.py controlgen` solves VenMar#45 (R-S 2007, Venus flyby,
   Mars massless, key `k2|LV>M/0s|LM>V/0s`) through the same `solve_structure` and `assess`. Exact
   zero (residual 7.1e-15 km/s), gate pass, V_inf Venus 8.2196 and Mars 12.9648 km/s; published
   8.22 and 12.96 (R-S 2007 Table 3), stored #942 values identical to 1e-9. Data:
   `control_generator_venmar45.json`.
