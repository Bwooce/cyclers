# #430 Global MGA-DSM precursor engine — Aldrin / S1L1 deliverable verdict

**Date:** 2026-06-24
**Status:** COMPLETE — capability advance over #307; honest partial result (no novel, no clean ballistic insertion). Both rows stay V0.

## What was run

The unified global MGA-DSM engine (`global_precursor_engine.py`, spec
`docs/superpowers/specs/2026-06-23-global-precursor-engine-design.md`) re-ran the
Aldrin (`aldrin-classic-em-k1-outbound`) and S1L1 (`s1l1-2syn-em-cpom`) precursor
scans via `find_cycler_precursors(..., use_global_engine=True)` over the
2030–2034 launch window. The engine: enumerates candidate sequences, then for
each runs a global `differential_evolution` over
`(epoch, per-leg TOFs, per-leg DSM)` against the real-DE440 `close_epoch_locked`
oracle, ranking survivors by total ΔV / `dv_band`. Offline KNOWN_CORPUS
literature check on each survivor. 14 survivors each; outputs
`data/precursor_430_{aldrin,s1l1}_global.jsonl`.

## Result vs the #307 baseline

| metric | #307 single/multi-rev (local) | #430 global engine |
|---|---|---|
| min flyby-continuity (Aldrin) | 4.15 km/s | **~0.00 km/s** |
| min flyby-continuity (S1L1) | 3.95 km/s | **~0.00 km/s** |
| precursor type found | none sub-gate | low_maintenance / powered_dsm |
| best Earth-launched vinf-match (Aldrin) | n/a (wall) | **0.21 km/s** (E-V-V-E, 75 m/s DSM) |
| best Earth-launched vinf-match (S1L1) | n/a (wall) | **0.44 km/s** (E-V-E, 259 m/s DSM) |
| literature-fresh | 0 / 394 | **0 / 14** |

**The #307 flyby-continuity wall fell.** The global search (the load-bearing
lever — *not* the eccentric seeds; see Caveat) plus per-leg DSMs drives
flyby-continuity to ~0 on every survivor and finds genuine **low_maintenance /
powered** precursors into both classic E-M cyclers — geometries the #307
local-optimiser-from-circular-seed method could not reach.

## Reading (with discipline scrutiny)

The run summary's headline ("closure 0.0022 km/s, strictly_ballistic") is
**misleading** and was checked before any claim (per `feedback_orbit_closure_discipline`
— "it closed!" is the danger signal):

- The engine's internal `closure_residual` is measured against each candidate's
  *own enumerated seed V∞ bins*, not the cycler's published seed. The honest
  insertion metric is `vinf_match_residual_kms` (terminal V∞ vs the cycler's
  6.5 / 5.65 km/s seed).
- The `closure≈0.0022` "best" is a `('V','E')` chain that is **not
  Earth-launched** (it starts at Venus) and whose `vinf_match` is 0.50 km/s.
  Half the 14 survivors per target are non-Earth-launched enumerator artifacts
  (the engine's sequence enumeration does not enforce an Earth start — a known
  scope gap inherited from `find_mga_chains`, logged as follow-on).
- Among **valid Earth-launched** chains, the best reaching the cycler's seed
  V∞ are **powered**: Aldrin `E-V-V-E` (vinf_match 0.21 km/s, 75 m/s DSM,
  low_maintenance); S1L1 `E-V-E` (0.44 km/s, 259 m/s DSM) / `E-E` (0.51 km/s,
  82 m/s DSM). The `strictly_ballistic`-banded survivors either miss the seed
  V∞ badly (`E-E-E` resonant loops, vinf_match ~6.5) or are the non-Earth-
  launched artifacts.

So: **no clean sub-gate *ballistic* insertion at the seed V∞**, and **0
literature-fresh** (every survivor resembles a published anchor — expected for
the well-studied E-M precursor space). The powered survivors are real low-ΔV
insertions but are (a) not novel and (b) ~0.2–0.5 km/s off the seed V∞, so none
warrants V0-V5 gauntlet promotion (the gauntlet is for *novel* candidates;
closure discipline: no novelty → no catalogue claim).

## Caveat — the "eccentric Tisserand" half is scaffolded, not wired

The build's headline pairs "eccentric-body Tisserand seeds + global search". The
**global search is fully delivered and is what produced the gain.** The
**eccentric-radius re-screen is NOT yet wired** into the seeder
(`eccentric_tp_linkable_radius_au` exists + is golden-tested, but
`eccentric_tp_seeds` still enumerates on the circular/mean-`a` Tisserand graph).
The result therefore reflects *global-search-over-circular-seeds + DSM*, which
already cleared the #307 wall — empirically confirming the prior analysis that
the global optimiser, not the seed eccentricity, was the binding lever. Wiring
the real-radius re-screen is a documented follow-on; its marginal value is now
testable (and likely small, given the global search already escapes the seed
basin).

## Disposition

- Both `aldrin-classic-em-k1-outbound` and `s1l1-2syn-em-cpom` precursor scans
  stay **V0**. No promotion, no novelty claim (0 fresh).
- Registered the method-versioned result in `data/empty_regions.jsonl` (the
  global engine capability-subsumes the #307 multi-rev negative): the
  capability advance is real, but the *novel ballistic precursor* remains empty.
- Follow-ons (logged, not blocking): (1) enforce Earth-launch in the engine's
  sequence enumeration (drop the non-Earth-start artifacts); (2) wire the
  eccentric-radius re-screen into `eccentric_tp_seeds`; (3) the powered
  survivors (E-V-V-E etc.) are candidates for a human design-reference check,
  not catalogue admission.
