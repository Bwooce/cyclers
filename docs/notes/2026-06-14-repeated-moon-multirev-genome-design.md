# Repeated-moon multi-rev cycler genome — design (Track A)

**Goal:** a SEARCH genome that can represent and close **repeated-moon
multi-revolution cyclers** — the Liang Callisto-Ganymede-Europa (CGE) class:
trajectories that orbit a planet doing several resonant revolutions between
*repeated* encounters with the *same* moons in a fixed sequence. Our current
genome (zero-rev / single-encounter / no repeated body) structurally cannot
express these, which is exactly why the Jovian region read EMPTY in the sweeps.

**Why this one first (evidence × tractability):** Liang et al. 2024 give **4
concrete ballistic CGE members** (per-flyby V∞/ToF/phases, ≥10 cycles) — a
*proven non-empty target*. We already adopted 3 at V1 (#222) + Member D n-body
(#223), but only as ROWS; we cannot yet SEARCH for new CGE-class members. This
genome turns the held machinery (moon-tour support #76/#117, VILM endgame #179,
Tisserand/V∞ graph, body-agnostic compute #75) into a discovery genome, and it
triggers the empty-region §6b capability-subsumption re-sweep of the Jovian bucket.

## Representation
A repeated-moon cycler genome is:
- **Moon sequence** `[m_1, m_2, …, m_k]` over a planet's moons (Jupiter: Io,
  Europa, Ganymede, Callisto), repeating each cycle.
- **Per inter-encounter leg:** a planet-centric resonance `p:q` + revolution count
  `n_rev`, realised as a resonant arc connected by a **V∞-leveraging transfer**
  (VILM) — i.e. the leg is a resonant orbit segment in the planet-centric
  patched-conic / CR3BP moon model.
- **Per encounter:** a flyby of moon `m_i` that bends the trajectory (B-plane /
  turn angle) to set up the next resonant leg, with **V∞-magnitude continuity**
  (ballistic) or a small bounded maneuver.
- **Periodicity (the cycler condition):** after the full k-encounter sequence the
  trajectory returns to the same relative geometry — the cycle repeats.

Decision vector: per-leg `(p:q discrete, n_rev discrete)`, per-encounter
`(moon id discrete, flyby B-plane continuous)`, plus overall epoch/phasing and the
repeat period. Discrete combinatorics (moon-seq × resonance) enumerated up front;
continuous DOF (flyby, timing) closed by the corrector.

## Closure / solver
- **Enumerate feasible sequences** with the Tisserand / V∞ graph (which moons link
  at which V∞ — body-agnostic Tisserand we already have). Bounded combinatorics.
- **Close each candidate** to periodicity (the repeated sequence returns to itself)
  via the VILM endgame solver (#179) + the planet-centric multiple-shooting /
  match-point corrector, with **flyby continuity** — the #226 FBS flyby-coupling
  block gives analytic flyby gradients here, a natural reuse.
- **One canonical residual** (periodicity match-point, km/s) + a gate, mirroring the
  #248 harness discipline (epoch-safe eval, no ad-hoc metrics).

## VALIDATION GATE (reproduce-before-search — mandatory)
Before searching for new members, the genome MUST reproduce **Liang's 4 known CGE
members** (the V1 rows from #222/#223) from sourced (V∞, ToF, phase) data — same
discipline as the Braik-Ross #236 gate. The published members are the goldens
(sourced; not our own output). Only after the genome recovers them do we trust it
to find new ones. A failure to reproduce = stop and fix, don't search.

## Build sequence (bite-sized, commit each)
1. **Moon-system registry + Tisserand/V∞ graph** for a planet's moons (reuse the
   body-agnostic Tisserand machinery). Test: known CGE V∞ links.
2. **Genome representation + decision vector** (moon-seq, per-leg p:q/n_rev,
   flyby). Test: round-trip + the Liang members encode validly.
3. **Periodicity corrector** for the repeated sequence (VILM + flyby continuity,
   canonical residual). Test: a single known leg closes.
4. **VALIDATION GATE:** reproduce Liang's 4 CGE members from sourced data. MUST
   pass before step 5.
5. **Search:** enumerate moon-seq × resonance, close, SILVER survivors → the
   unchanged V0–V5 gauntlet.
6. **Daemon-host** the search (Track-C) — it's long combinatorial compute, the
   daemon's natural job (NOT one-shot agents; see #248's three hangs).

## Scope & discipline
- Start with **ONE planet system (Jupiter)** — Liang gives the validation gate.
  Generalize to Saturn (Enceladus/Titan tours) afterward.
- Reproduce-before-trust (Liang gate); same-model / sourced goldens only.
- **NO catalogue writeback** until a found member passes the gauntlet — a closed
  search candidate is SILVER, not a validated cycler.
- This is the Track-A genome effort the discovery-program spec calls the highest
  evidenced-payoff lever (the only genome with a *proven* non-empty target).

## Success criterion
The genome reproduces Liang's CGE members (gate), THEN the bounded search returns
≥1 NEW repeated-moon cycler that survives to the gauntlet — a genuine novel cycler
in a topology the prior genome could not represent. A clean negative (the gate
passes but the search finds nothing new in the swept moon-seq/resonance box) is
also valuable: it bounds the Jovian region under the *new* genome (a stronger
empty-region claim than the zero-rev sweep).
