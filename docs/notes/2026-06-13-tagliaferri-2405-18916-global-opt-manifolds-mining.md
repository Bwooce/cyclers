# Tagliaferri, Blazquez, Acciarini & Izzo 2024 — "Global Optimization for Trajectory Design via Invariant Manifolds in the Earth-Moon CR3BP" (discovery-capability digest)

**Date:** 2026-06-13
**Task:** #250 discovery-capability sweep (head-start candidate 2, Track B).
**Source (free):** F. Tagliaferri, E. Blazquez, G. Acciarini, D. Izzo, ESA
Advanced Concepts Team. arXiv:2405.18916 (29 May 2024); also ISSFD 2024 paper
18-6 (issfd.org/ISSFD_2024/ISSFD2024_18-6.pdf).
**Access:** abstract + ISSFD/ResearchGate metadata via WebSearch; `WebFetch`
denied (no table transcription).
**Writeback: NONE** (digest only).

---

## Verdict: BACKGROUND (capability we already have)

The head-start brief framed this as "multiparameter continuation of transfer
FAMILIES → family-continuation discovery method." The actual content is
narrower than that framing:

- **Problem:** impulsive minimum-propellant TRANSFERS (not cyclers, not new
  periodic families) between **halo** and **vertical-Lyapunov** orbits about
  Earth-Moon L1/L2, at the **same Jacobi constant**, via their invariant
  manifolds.
- **Method:** Monotonic Basin Hopping (MBH) + an SQP local solver in a parallel
  optimization framework, swept over an interval of Jacobi constants and a max
  TOF. The "families" are families of *transfers* (the low-ΔV solutions vary
  continuously with C_J and TOF), not new families of orbits/cyclers.
- **Result:** a systematic map of very-low-ΔV manifold-mediated transfers across
  the chosen C_J / TOF window.

Every ingredient is already in our earned toolkit:

- **MBH** — we have it (`docs/notes/2026-06-07-mbh-wrapper.md`,
  `2026-06-07-englander-2014-mbh-tuning-mining.md`).
- **Invariant-manifold seeding** — our CR3BP lane + the Ross tube-overlap
  seeding (ADOPTED) already do manifold-based seeding.
- **Same-Jacobi family-to-family accessibility** — this is exactly the
  Braik-Ross orbital-network question (#230, ADOPTED reachable-set scorer), at a
  finer resolution but for transfers rather than the network abstraction.

So this paper restates, with an ESA-ACT MBH+SQP wrapper, a transfer-design
capability we already possess. It does NOT introduce a new genome, a new cycler
topology, or a discovery method that reaches a region our current stack cannot.

## Discovery vs anchor

**Neither.** Output is halo↔vertical-Lyapunov transfers — not cyclers, not new
orbit families, no catalogue rows. It is a transfer-optimization study.

## Reproducible data

The reproducible artifact is the MBH+SQP pipeline and the low-ΔV transfer maps
over (C_J, TOF). Of marginal interest: the ESA-ACT lineage means the code style
overlaps Izzo's `pykep`/`pygmo` ecosystem, which we could lift MBH-tuning ideas
from — but we already tuned MBH (Englander 2014). No cycler ICs/tuples.

## Method capability

None new for us. The one transferable nuance is the parallel MBH+SQP framing for
manifold-seeded transfers at fixed C_J — a possible micro-optimization of our
transfer solver, not a discovery capability. The genuinely discovery-relevant
sibling in this lineage is **Li, Gondelach, Izzo et al. diffusion-model global
trajectory search** (arXiv:2410.02976, surfaced in the same search cluster) —
that is the ML-generative-prior direction, downstream of Track B per the spec
(#240-adjacent), and a separate background item.

## Proposed follow-on

**Background only.** Cross-reference from the MBH wrapper note as an
ESA-ACT manifold-transfer reference. Do not adopt — no capability gap is filled.
Re-derivation effort would duplicate the Braik-Ross accessibility work and our
existing MBH stack.
