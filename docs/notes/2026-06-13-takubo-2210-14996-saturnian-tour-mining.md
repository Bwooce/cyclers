# Takubo, Landau & Anderson 2023/2024 — "Automated Tour Design in the Saturnian System" (discovery-capability digest)

**Date:** 2026-06-13
**Task:** #250 discovery-capability sweep (background candidate 5; cross-check
vs our VILM endgame, #179).
**Source (free):** Y. Takubo, D. Landau, B. Anderson, arXiv:2210.14996 (26 Oct
2022, math.OC); AAS Space Flight Mechanics 2023 (Austin); published *Celestial
Mechanics and Dynamical Astronomy* 136:8 (2024), DOI 10.1007/s10569-023-10179-8.
**Access:** abstract + ADS/Springer metadata via WebSearch; `WebFetch` denied.
**Writeback: NONE** (digest only).

---

## Verdict: BACKGROUND (moon-tour lane backfill; method we substantially have)

## What it does

Multi-objective **moon-tour** design from **Titan to Enceladus**, using
"resonance family hopping": progressively transfer between moons by exploiting
V∞ on resonant orbits. The search space (combinations of resonances × flyby
speeds) is large, so the method:

1. precomputes a **database of V∞-leveraging (VILM) transfers**, then
2. runs **grid-based dynamic programming** over that database to optimize the
   whole Titan→Enceladus tour (multi-objective: ΔV, TOF), yielding a complete
   trade space at global optimality in tractable time, and
3. **validates the result in a full-ephemeris model.**

## Cross-check vs our VILM endgame (#179)

This is the *application* layer on top of the VILM/Leveraging-Graph theory we
already mined (Strange/Campagnola, `2026-06-05-endgame-tisserand-mining.md`):

- We hold the **VILM closed-form maneuver formulae, the Leveraging Graph, the
  n:m±(K) VILM taxonomy, and the theoretical-minimum-ΔV multi-VILM result.**
- Takubo et al. add a **DP search wrapper over a VILM database** to assemble VILM
  legs into a full optimized moon tour. The endgame *primitive* is the same; the
  contribution is the global-optimization tour assembly + the explicit
  CR3BP→ephemeris transition validation.

So the *capability* (VILM-based moon-tour construction) is already in our
toolkit at the primitive level; this paper is a known-good reference for how to
*assemble* those primitives into a globally-optimal tour, plus a Saturnian
worked example.

## Discovery vs anchor

**Neither, for cyclers.** A Titan→Enceladus tour is a one-way transfer tour, not
a repeating cycler. It produces no cycler catalogue rows. Relevant only as
moon-tour-lane context and as a VILM-assembly method reference.

## Reproducible data

The VILM-transfer database structure + the DP formulation are the reproducible
assets. The Saturnian worked example (Titan/Enceladus resonances, flyby V∞
values) is a potential cross-check set for our VILM leg code if we ever exercise
the Saturnian moon-tour lane, but it is not cycler data. Full numeric tables not
transcribable here (WebFetch denied) — note for optional user fetch if the
Saturnian tour lane is ever revived.

## Proposed follow-on

**Background only.** Cross-reference from `2026-06-05-endgame-tisserand-mining.md`
as the tour-assembly application of the VILM theory. No adoption: the DP
tour-assembly is downstream of moon-tour-lane work that is not on the current
discovery critical path (which is Earth-Moon ballistic cyclers + the
prioritizer + the daemon). Revisit only if a Saturnian/Jovian moon-tour
*cycler* discovery campaign is opened.
