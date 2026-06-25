# Next-frontier prioritization — where new orbits actually live (2026-06-26)

Successor to `2026-06-25-discovery-strategy-prioritization-design-draft.md`. Written
after the 2026-06-25→26 session mined the canonical CR3BP/ER3BP/3D periodic-orbit
frontier to its floor. Purpose: aim future discovery at the genuinely-fertile
territory instead of re-mining the dry one.

## The load-bearing finding
The **canonical CR3BP / ER3BP / 3D periodic-orbit frontier** (Earth-Moon, Sun-planet)
is largely **mined out**. The session's three discovery regions all returned
negative-for-rows:
- Region A (isolated elliptic resonant, #442/#457): now *reachable* (a real capability
  win) but the members are published reproductions — **not novel, not cyclers**.
- Region B (stable 3D C21/C32, #447): published-class members, **no transport utility**.
- Region C (high-e Sun-planet ER3BP, #448): families persist, **no V∞ benefit**.

The DA/HOTM lane (#450) likewise proved a capability (re-opens EM C≈3.0 from a coarse
grid) but finds **periodic orbits, not cyclers**. Net: this model-class yields
capability + honest negatives, **not new catalogue rows**. More sweeps here have low
expected yield.

## The proof-point — where the ONE novel row came from
The only genuinely-novel catalogue row produced across this entire arc is
**`umbriel-oberon-1-1-uranian-quasi-cycler-2026` (#339)**: `source: discovered`,
`first_published: cyclerfinder 2026`, `corroborating_sources: []`, cleared the
41-anchor lit-novelty check, **V4-validated** on real URA111 ephemeris. It did NOT
come from the strict-cycler CR3BP frontier — it came from the **quasi-cycler /
epoch-locked (QP-tori) frontier**. That is the single strongest signal we have about
where the live territory is.

## Fertile directions, ranked by expected new-row yield

1. **Moon-tours via the releg genomes — HIGHEST.** (#449 powered/DSM releg, building
   now; #451 3D-inclined, designed.) A *different object class* from periodic orbits:
   multi-moon resonant tours in the Jovian/Saturnian/Uranian/Neptunian satellite
   systems — the largest untapped combinatorial space we have, with sourced VILM
   goldens. #449 may flag a catalogue candidate directly. **This is the right primary bet.**
2. **Quasi-cyclers / epoch-locked (QP-tori) — HIGH (proven fertile).** The #339
   frontier. #333 (QP-GMOS continuation) extends exactly the lane that produced our one
   novel row; #320-class sweeps over more satellite systems are the natural follow-on.
3. **Cislunar BCR4BP / Belbruno BCT — MEDIUM (under-mined, harder).** (#378, design in
   progress.) The Sun-Earth-Moon coupled cislunar region via weak-stability-boundary /
   ballistic-capture transfers — methodology digested, lane unbuilt. Honest caveat: less
   directly cycler-shaped; may be capability/transfer-construction more than catalogue rows.
4. **Low-thrust / indirect lanes — MEDIUM.** (#309 built; #359 indirect TPBVP companion
   pending.) Opens powered transport objects the ballistic lanes can't host.

## Deprioritize (the dry well)
More CR3BP / ER3BP / 3D periodic-orbit enumeration in Earth-Moon / Sun-planet systems.
Mined to the floor; rigorously stamped in `empty_regions.jsonl`. Re-open ONLY when a
strictly-more-capable method ships (the registry already encodes those re-open keys).

## Two structural truths to plan around
- **Much of what "remains" is infrastructure, not discovery.** The #305/#306 gauntlets
  *promote* already-found candidates (BCR4BP families from #303/#304; moontour/3D
  candidates) to validated catalogue rows. Real value — but it grows the catalogue from
  existing finds, not new search. Build these to *unlock promotions*, not expecting them
  to find orbits.
- **Growth is increasingly new-input-gated.** Past this point, new rows need new *inputs*
  — new systems, higher model fidelity, new methods/papers — not more iteration in the
  mined space (the validation-ceiling lesson, applied to discovery). The literature-
  monitoring lane (#265) and targeted acquisitions matter more now than another sweep.
  This is a maturing-catalogue signal, not a dead end.

## Recommended arc
`#449` (in flight, primary bet) → `#333` (plan → implement; proven-fertile quasi-cycler
lane) → `#378` (design → implement; under-mined cislunar) → the `#305/#306` gauntlets to
promote whatever the above + prior campaigns surfaced. Treat new-input acquisition (#265
monitoring, fresh systems/fidelity) as a standing parallel lever, not an afterthought.
