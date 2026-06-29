# #480 correction — the resonant-conic seed is structurally NON-ballistic (and the 3-D B-plane hypothesis is withdrawn)

**Date:** 2026-06-29 (later same day; supersedes the "next lever" of
`docs/notes/2026-06-29-480-eggie-zerosoi-verdict.md` and walks back the
"basin solved" framing of the Stage-1/2/3 notes).
**Method:** forward-verification (not search) — `scripts/_eggie_fwd_verify_480.py`.

## What I checked

For each e, build the rigid resonant-conic guess `eggie_initial_guess(e)`
(`theta_dep_europa` is a pure gauge — V∞/ToFs invariant — so e is the only DOF),
and report the EXACT paper Eq 3-5 flyby ΔV + periodicity wrap + per-encounter V∞ +
flyby altitudes against Table 4. No optimization, no corrector — does the
construction itself realize a ballistic EGGIE?

## Result — NO, and it exposes the real gap

At e≈0.62 (the spike-4 "all V∞ on target" optimum):

| node | guess V∞ | Table 4 | guess alt | Table 4 alt |
|---|---|---|---|---|
| Europa depart | 9.08 | 9.12 | — | 1444 |
| Ganymede #1 | 6.78 | 7.07 | 13109 km | 2155 |
| **Ganymede #2** | **4.17** | **7.07** | 15661 km | 6263 |
| Io | 6.23 | 8.38 | **−535 km (sub-surface!)** | 653 |
| Europa arrive | 7.51 | 9.12 | — | — |

interior Eq 3-5 flyby ΔV ≈ 2493 / 1913 / 1243 m/s, periodicity wrap ≈ 2596 m/s,
total ≈ 8.2 km/s; `all_flybys_feasible=False`. Across e∈[0.55,0.75] the two Ganymede
V∞ are NEVER equal and never both 7.07 (Ganymede #2 swings 4.0↔9.8 km/s as the
Lambert branch flips at e≈0.64), and the Io flyby is sub-surface over much of the
range.

## The real diagnosis (corrects three earlier claims)

1. **"Basin solved" was over-stated.** Spike-4's "single resonant conic puts all
   three V∞ on the Table-4 targets" is true of the **analytic conic magnitude**
   |V∞|(a,e,r) — but that is **necessary, not sufficient**. The realized tour is a
   chain of Lambert legs whose V∞ at each node is set by the leg, and the
   construction only pins the *departure* V∞ near the conic (branch chosen closest
   to the conic velocity); it does **not** enforce the **equal in/out V∞**
   ballistic-flyby constraint. So Ganymede #2 (arrival of the G→G leg) collapses to
   4.17, the equal-7.07-km/s Ganymede property (the paper's defining ballistic
   signature) is broken, and the Io flyby is infeasible.
2. **The 3-D B-plane DOF hypothesis is WITHDRAWN.** The paper's ideal model is
   circular+coplanar, and Lambert legs between coplanar moons are coplanar (transfer
   plane = the moon plane), so θ_B has no out-of-plane freedom there — yet the paper
   gets 0.70 m/s. 3-D is not the differentiator; the seed construction is.
3. **Why no corrector closed it:** all four correctors (FD/STM/epoch/sub-arc) and the
   zero-SOI search were started from a seed that is structurally non-ballistic (G#2
   and Io off by km/s, Io sub-surface). Local correction from a structurally-wrong
   seed lands in the documented ~0.1 km/s (n-body) / ~1 km/s (zero-SOI) basins — it
   was never going to reach 0.70 m/s. The wall was the SEED, not the solver.

## The actual next lever (fresh build)

Construct the EGGIE as a genuine **ballistic patched-conic cycler**: a chain where
each flyby is enforced to preserve |V∞| (equal in/out) AND deliver the required bend
within the 25–70,000 km altitude window, AND the cycle closes — solved as a coupled
system over (e, the resonant phasing, per-leg ToFs, Lambert branch/rev per leg),
seeded by the conic guess. This is the paper's full method (Tisserand/conic guess →
Lambert + Eq 3-7 flyby → Monte-Carlo closure) done with the ballistic constraint
IN the system, not just V∞-magnitude matching. The reusable pieces all exist
(`resonant_conic`, `flyby_maneuver_dv`, `jovian_stm`); what's missing is the
ballistic-cycler **constraint formulation** (equal-V∞ + bend-feasible + periodic),
which the current `eggie_initial_guess`/`refine_eggie` do not impose.

Best treated as a fresh, well-scoped build (the equal-V∞ G→G resonant-return leg is
the crux — it must be a true resonant transfer, not a min-distance Lambert pick).
Scratch: `scripts/_eggie_fwd_verify_480.py` (removed).
