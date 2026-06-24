# #438 broadened 3D broken-plane sweep — verdict

**Date:** 2026-06-25. The #434 follow-on: lift ALL FOUR Braik-Ross Earth-Moon
cyclers (C11a/C11b/C21/C32) into z≠0, not just C11a. Result: **the broadening
worked — C11b/C21/C32 yield NEW 3D topologies and (notably) STABLE out-of-plane
families.** This RETIRES the #434 "seed-limited to one planar root" caveat and
surfaces the strongest 3D-cycler novelty candidates to date. Report-only; novelty
claim gated on the 3D V0-V5 gauntlet (#306) + a real lit-check.

## What ran
`scripts/scan_434_3d_broken_plane_em.py` (broadened by #438 + parallelised),
`braik_ross_system()` (all 4 cyclers confirm-recovered), 22 units across 12
workers, wall 12067 s. **2260 converged members**, every one passing the compound
gate (corrector + independent Radau closure); closure median 2.75e-11, max 2.06e-8.
Output `data/scan_434_3d_broken_plane_em.jsonl`.

## Result — per cycler (the broadening is the story)

| seed | members | (k1,k2,k_z) | z0 range | C range | stability |
|---|---|---|---|---|---|
| C11a (1,1) ×5 z0 | 912 | (1,0,6), (1,1,8), (1,1,10), (1,2,10) | −0.24..0.24 | 2.99..3.12 | hyperbolic (the known #287 spike + branches) |
| C11b (1,1) | 201 | (2,0,10) | −0.69..−0.24 | 2.03..3.02 | hyperbolic |
| **C21 (2,1) ×2** | **402** | **(2,0,10)** | −0.65..−0.21 | 2.15..3.03 | **107/201 STABLE** + 69 hyp + 25 unstable |
| **C32 (3,2) z0_0.24** | **201** | **(8,0,26)** | 0.14..0.35 | 2.68..3.06 | **164/201 STABLE** + 37 hyp |
| C32 (3,2) z0_0.15 | 201 | (0,0,6) | −0.13..0.13 | 1.73..2.42 | hyperbolic |
| C32 (3,2) z0_0.10 | 142 | (66,0,136) ⚠ | −0.19..−0.12 | **1.79..10.75 ⚠** | 45 stable + 97 hyp |
| L1 vertical-Lyapunov | 201 | (0,0,2), (0,1,2) | −0.19..0.23 | 3.00..3.17 | hyperbolic + 14 stable |

## Reading

**The headline:** broadening from C11a to all four cyclers surfaced **stable
out-of-plane families off the differently-wound C21 (2,1) and C32 (3,2) planar
cyclers** — new (k1,k2,k_z) topologies the C11a-only #434 run could never reach.
Stable 3D cycler-class families are the genuine prize (most 3D extensions are
hyperbolic). These are the strongest 3D-cycler novelty candidates the project has
produced.

**Honest caveats (the candidates are NOT yet a discovery claim):**
1. **They are 3D extensions of KNOWN planar Aldrin/Braik-Ross cyclers** (C21, C32),
   so "novel" means *novel out-of-plane structure on a known planar root*, not a
   from-scratch novel cycler.
2. **Earth-Moon, not coordinate-matchable** to the 3-anchor spatial-CR3BP corpus
   (Howell/Folta/Antoniadou are µ≠EM or taxonomy-only) — a lit-check would be the
   same sparse-corpus false-negative as #434. Real novelty adjudication needs an
   Earth-Moon 3D catalogue (Breakwell & Brown 1979; the JPL 3-body DB).
3. **The (66,0,136) C32 z0_0.10 family is SUSPECT:** its Jacobi spans an order of
   magnitude (1.79→10.75) within "one" continuation, and the winding numbers are
   extreme — almost certainly the continuation wandering across families / a
   winding-topology miscount on a near-degenerate arc, NOT a single coherent
   family. Excluded from the candidate set pending scrutiny.
4. Per the #436 closure-discipline lesson the closure gate is already satisfied
   per-member (genuine POs), but a *novelty* claim still needs independent
   re-verification of the stable C21/C32 families + the V0-V5 3D gauntlet.

## Disposition
- **No catalogue writeback** (report-only; matches discipline).
- **Retires the #434 "seed-limited" caveat** — broadening to all four planar
  roots DID surface new structure (the caveat's prediction confirmed).
- **The stable C21 (2,0,10) and C32 (8,0,26) out-of-plane families are the
  strongest 3D-cycler candidates** → run them through the 3D V0-V5 gauntlet
  (**#306**, pending) + an Earth-Moon-3D lit-check before any novelty/catalogue
  claim. That gauntlet is now the gating next step for a 3D discovery.
- Method-versioned registry entry to follow (supersedes the C11a-only #434 entry).
- Follow-up scrutiny: confirm the (66,0,136) family is a continuation artifact;
  independently re-verify a sample of the stable C21/C32 members.
