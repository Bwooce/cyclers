# JPL Three-Body Periodic Orbit Catalog — independent-oracle cross-check (#116)

**Date:** 2026-06-13
**Client:** `src/cyclerfinder/verify/jpl_periodic_orbits.py` (live API, urllib;
parser pinned by `tests/verify/test_jpl_periodic_orbits.py` against a saved
4-row fixture).
**Driver:** `scripts/jpl_oracle_crosscheck.py` (live, network).
**Source:** JPL SSD periodic_orbits API
(https://ssd-api.jpl.nasa.gov/doc/periodic_orbits.html) — machine-readable,
JPL-published CR3BP ICs (x,y,z,vx,vy,vz, Jacobi, period, stability) per family.

## What this is (and is NOT)

The JPL catalog is an **independent reproduction by a separate group** of CR3BP
periodic orbits — exactly the per-interface external anchor the false-consensus
discipline wants. Used here as a cross-check ORACLE for our CR3BP propagator +
Jacobi, NOT as a catalogue import: halo/DRO/Lyapunov orbits are not cyclers.
Only a genuine *cycler* family (none in the standard family list) would feed the
catalogue, and then as V0-sourced seeds.

## Convention reconciliation (mandatory, done first)

JPL Earth-Moon mass ratio `mu = 0.012150585609624`; ours
`0.0121505843946952` (Ross-sourced, `core.cr3bp.cr3bp_system`). **abs diff
1.215e-9, rel diff 1.0e-7.** An IC tuned to one mu and propagated under the other
cannot re-close better than this offset allows — so the cross-check carries BOTH
mu values explicitly (`reconcile_mu`) and reports closure under each.

## Result (Earth-Moon L2 halo, N branch, 1535 orbits; sample across the family)

| idx | JPL C | our C @ our mu | dC | closure @ our mu | closure @ JPL mu |
|---|---|---|---|---|---|
| 0    | 3.015178 | 3.015178 | -6.5e-9  | 1.42e-7 | 1.50e-11 |
| 383  | 3.027259 | 3.027259 | -7.4e-9  | 6.87e-8 | 6.47e-12 |
| 767  | 3.061024 | 3.061024 | -6.7e-9  | 1.11e-6 | 2.43e-11 |
| 1151 | 3.114546 | 3.114546 | -7.4e-9  | 4.09e-6 | 1.21e-10 |
| 1534 | 3.158445 | 3.158445 | -1.7e-8  | 3.67e-8 | 2.76e-10 |

**Verdict — our CR3BP stack reproduces JPL's catalog independently.**
- **Closure @ JPL mu (mu matched): worst 2.76e-10 nd** over the sampled family.
  Our DOP853 propagator re-closes JPL's published ICs to sub-nanometre-scale
  nondimensional residual — an independent confirmation of the propagator.
- **Closure @ our mu degrades to ~1e-6** — purely the 1e-7 mu offset acting over
  ~one period (it does NOT indicate integrator disagreement; the mu-matched
  column is the controlled comparison).
- **Jacobi agreement: dC ~ 6e-9..2e-8** across C ∈ [3.015, 3.158] — our
  `jacobi_constant` matches JPL's published C to ~1e-8 (the residual tracks the
  mu offset, as expected since C depends on mu).

## Use going forward

- **Cross-check oracle** for the CR3BP continuation campaign (#218/#219) and the
  Ross/Cuevas rows: re-derive a row, then confirm it against the nearest JPL
  family member at matched mu (the V5 "independent reproduction by a separate
  group" angle). Always reconcile mu first.
- **Targeted sourced ICs**: if a cycler family ever appears in the JPL family
  list, pull those specific ICs as V0-sourced seeds (not a bulk import).
- The client is offline-safe for CI (tests use the fixture); only the driver
  script hits the network.

Closes the #116 "JPL 3-Body catalog" line as a *reachable capability* (it is a
live API, not a document to acquire).
