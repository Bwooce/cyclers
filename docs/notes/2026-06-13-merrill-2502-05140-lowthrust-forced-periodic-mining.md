# Merrill, Kulik, Bryan & Savransky 2025 — "Mass-Optimal Low-Thrust Forced Periodic Trajectories in the Earth-Moon CR3BP" (discovery-capability digest)

**Date:** 2026-06-13
**Task:** #250 discovery-capability sweep (head-start candidate 1, Track A).
**Source (free):** C. C. Merrill, J. Kulik, M. J. Bryan, D. Savransky,
arXiv:2502.05140 (7 Feb 2025), Cornell. Companion: Merrill & Kulik, "Generation
of Energy-Optimal Low-Thrust Forced Periodic [Trajectories]" (Semantic Scholar
2114c448…) — the energy-optimal predecessor this mass-optimal paper supersets.
**Access:** abstract + literature-review mirror (themoonlight.io) via WebSearch;
`WebFetch` denied, so the optimal-control formulation detail below is from the
abstract/review, not the full derivation.
**Writeback: NONE** (digest only).

---

## Verdict: BACKGROUND (Track-A-adjacent, but NOT a cycler genome)

This is a **forced-periodic low-thrust** paper, and the distinction is
load-bearing for us: a "forced periodic trajectory" is one whose periodic
structure is *maintained by continuous thrust* — it is not a ballistic
(natural) periodic orbit and not a ballistic cycler. The project deliverable is
**novel ballistic cyclers** (repeating, propulsion-free or impulsive-flyby
geometry). Powered station-kept "cyclers" are explicitly out of catalogue scope
per the cyclers-only rule (cf. the Genova-Aldrin reclassification and the
withheld powered concepts). So the *output* class is not catalogue-grade for us.

The *method*, however, touches two existing tracks:

- **Track A (low-thrust genome).** The spec lists a low-thrust genome as a
  Track-A priority ("Sims-Flanagan machinery exists; no sourced powered rows —
  but discovery need not be sourced"). This paper is a mass-optimal
  continuous-thrust formulation that *enforces periodicity as a constraint* in
  the Earth-Moon CR3BP. If we ever build the low-thrust discovery genome, this
  is a relevant reference formulation — but it is one of several (Sullivan 2022
  PhD, Das-Sharma 2022 AeMiS, our existing Sims-Flanagan machinery already
  cover the construction primitives). It adds no capability our FBS/Sims-Flanagan
  stack lacks for *finding* powered periodic structure.
- **Track B (reachable sets).** The paper's headline result is reachable-set
  geometry: "the thrust-limited mass-optimal reachable set is a SUPERSET of the
  energy-limited energy-optimal reachable set in the xy-plane." This is a
  reachable-set ordering result. We already hold two reachable-set scorers —
  Braik-Ross energy-preserving (#236, ADOPTED) and Zhou-Armellin impulsive
  (#239). This adds a *mass-optimal low-thrust* reachable set as a third
  flavour. Marginally interesting as a future reachable-set complement, but the
  impulsive/energy-preserving pair already covers our prioritizer needs, and the
  mass-optimal variant is the most expensive to evaluate (full low-thrust
  optimal-control inner loop per seed) — wrong cost profile for a screen.

## Discovery vs anchor

**Neither, for our deliverable.** It does not find ballistic cyclers (its
periodic trajectories are thrust-forced, out of scope), and it provides no
sourced catalogue rows. It is a method reference for a track (low-thrust genome)
that is itself deferred behind Tracks A-ballistic and B in the discovery spec.

## Reproducible data

The paper studies reachable sets *around a reference trajectory*; it is a
methodology paper, not a family catalogue. No (μ, C, T) cycler tuples, no
ballistic ICs. The optimal-control formulation (mass-optimal indirect/Pontryagin
with periodicity constraints) could be reproduced, but only matters if/when the
low-thrust discovery genome is built — which is downstream of Tracks A+B per the
spec.

## Proposed follow-on

**Background only.** Park as a low-thrust-genome reference alongside Sullivan
2022 and the Das-Sharma low-thrust periodic work. Do NOT prioritize: the
low-thrust genome sits behind ballistic Track-A (repeated-moon multi-rev) and
Track-B in the discovery program, and this paper's forced-periodic output is
out of catalogue scope. Note the companion energy-optimal paper (Merrill &
Kulik) as the same lineage. Adjacent newer item from the same search cluster:
arXiv:2605.23770 "Reachability for Low-Thrust Trajectories via Maximum Initial
Mass" (2026) — same reachability-of-low-thrust theme, same background tier.
