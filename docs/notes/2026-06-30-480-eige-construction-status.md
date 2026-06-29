# #480 — EIGE positive-control construction: scoped + first-attempt obstacle (resume fresh)

**Date:** 2026-06-30. Status: the EIGE positive control (to validate the maintenance-ΔV
method before trusting an EGGIE number) requires a full EIGE resonant-conic CONSTRUCTION.
Energy is solved; the encounter-topology is the open work. Pausing to resume fresh (a clean
sub-project) rather than grind it error-prone at the tail of a long session.

## What's solved (constructed, not guessed)
- EIGE is 1-synodic. Resonance 1:1 (Table 1) → `resonant_sma(1,1,T_syn)` = **a_sc = ideal
  a_Ganymede = 1,055,289 km**, T_sc = 7.004 d. Eccentricity band [0.60, 0.93].
- Analytic conic V∞ over e: the EIGE/1-synodic regime is **HIGH excess speed ~12-16 km/s**
  (the paper's "hard to fly" 1-synodic class — Table 3 EGIEIE is 12-16), NOT the 7-9 km/s of
  the 4-synodic EGGIE. (Corrects a prior spec error that asserted 7-9 for EIGE.)
- Sourced gates for the eventual control (cite, do not assert): Fig-5 alts Io 2817 /
  Ganymede 13180 / Europa 470 km; ballistic first cycle; maintenance ΔV → ~30 m/s over 10
  cycles (digest pp.10-11).

## The open obstacle (first attempt failed, informatively)
A hand-derived static topology (E inbound, I inbound, G outbound) gives non-monotonic
encounter times: over ONE rev the conic radius crossings occur I-out, E-out, G-out, G-in,
E-in, I-in, and **E-before-I-before-G does not appear in a single rev's static order**. So:
- The EIGE encounter order is set by the moons' MOTION (resonance phasing), not a static
  crossing order — it needs the full phasing enumeration (the EGGIE-class hard part), with
  the self-consistency guard (`search/tour_self_consistency.py`) on every encounter incl. the
  repeated Europa.
- Re-check the rev count: the sibling EGIEIE (Table 3) is 1-synodic / **2** revs (paper p.7),
  so EIGE's "1 rev" may be 1:2, not 1:1 — verify against the PDF before fixing the topology.

## Resume plan (fresh sub-project, no guessing)
1. Confirm EIGE rev count from the PDF (1:1 vs 1:2).
2. Enumerate conic crossings over the rev(s); solve the (e, moon-phase) configuration where
   the E-I-G-E encounters all coincide with moons within SOI (the resonance-phasing search),
   validated by `assert_encounters_self_consistent`. Derive V∞ (output); gate altitudes vs the
   Fig-5 sourced values.
3. Run the (CGCEC-validated, #223) chained-maintenance lane on the constructed EIGE → must
   reproduce ballistic-first-cycle + ~30 m/s/10 (the positive control). Only then run EGGIE.

## Context
The maintenance MACHINERY is already validated (chain_cycles on Liang CGCEC, #223;
generalized to arbitrary sequences in `cf1f72a`). The EGGIE seed is constructed (not guessed)
and ~0.5 km/s below Table-4. Core #480 is settled (bug fixed + guarded, class reproduced,
exact member not) — this maintenance-ΔV work is refinement with no catalogue impact.
