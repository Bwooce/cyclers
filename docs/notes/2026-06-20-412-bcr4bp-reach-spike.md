# #412 — BCR4BP Phase-B reach spike (EM-libration family)

Date: 2026-06-20. Scoping gate for the full Phase-B build (a both-region-spanning,
synodic-resonant cross-system cycle in the coherent BCR4BP). Necessary-condition
measurement; asserts NO orbit. Reproduce: `uv run python scripts/spike_412_bcr4bp_reach.py`.

## Question
The #405/#411 cross-system cycle needs an object whose support spans BOTH the Earth-Moon
libration region AND the Sun-Earth-L region (~3.9 lunar distances, LD, from Earth = the
1.5e6 km patch plane). Does Sun-perturbing an **Earth-Moon-libration** orbit (EM-L1
Lyapunov) stretch its Earth-relative reach toward the SE-L region as μ_sun turns on? If yes,
the EM-libration family is a viable Phase-B vehicle; if it stays Earth-Moon-Hill-bounded
(~1 LD), that family is the wrong seed and a full Phase-B build off it is not justified.

## Method
For four EM-L1 Lyapunov amplitudes (C = 3.18, 3.1294, 3.05, 2.95): build the CR3BP orbit
(`correct_symmetric_fixed_jacobi`), cast it as a BCR4BP-at-μ_sun=0 seed
(`correct_bcr4bp_periodic`), geometrically continue μ_sun toward the full Sun value
(`continue_bcr4bp_family_in_musun`, 40 steps), and measure max Earth-relative distance over
one period (nondim = LD) for the seed and the last converged member. Reuse-only; every
corrector/continuation routine is already unit-tested.

## Result — STAYS EM-BOUNDED
| C | CR3BP reach | μ_sun reached (members) | BCR4BP reach |
|---|---|---|---|
| 3.18    | — (BCR4BP@0 seed failed to close; skip) | — | — |
| 3.1294  | 0.90 LD | 2.30e1 (10) | 0.90 LD |
| 3.05    | 0.96 LD | 4.42e1 (12) | 0.96 LD |
| 2.95    | 1.16 LD | 1.64e1 (9)  | 1.16 LD |

**Best reach 1.16 LD vs the 3.90 LD SE-L target.** Two independent signals, both pointing
the same way:

1. **Reach is flat under the Sun term.** Over the converged μ_sun range the Earth-relative
   reach does not trend upward toward 3.9 LD — the orbit stays an Earth-Moon-Hill-scale
   object (~1 LD). Turning on the Sun does not stretch an EM-libration orbit out to the
   Sun-Earth region.
2. **The family does not survive to full Sun strength.** Each continuation broke after
   9–12 of 40 steps (geometric stepping off a μ_sun≈1 floor reached only μ_sun≈16–44, far
   short of the full 328900). A divergent step ends the walk — the EM-libration family
   breaks down well before full solar perturbation.

## Conclusion → next step
The **EM-libration family is the wrong vehicle** for a both-region-spanning cross-system
object: it stays Earth-Moon-bounded (~1 LD) and breaks under μ_sun continuation before full
Sun strength. A full Phase-B build seeded from EM-libration orbits is **NOT justified**.

This does not kill BCR4BP Phase B universally — a Sun-Earth-scale seed (SE-libration orbit
or a transfer manifold) is the geometrically correct vehicle, but that is a from-scratch
build, not a continuation of the EM family. Given that, the productive near-term path
remains the **#411 full multi-revolution θ-closure corrector** in the patched model, whose
feasibility was already confirmed (n_em=41, n_se=19 → ~11 yr; see
2026-06-20-411-theta-closure-feasibility.md). #412 stays open as the longer-horizon
coherent-model venue, re-scoped to a SE-scale seed.

## Honesty caveat
Necessary-condition reach measurement only; no orbit asserted. The C=3.18 BCR4BP@μ_sun=0
seed did not close (skipped) — not investigated, as the three converged amplitudes already
answer the scoping question consistently.
