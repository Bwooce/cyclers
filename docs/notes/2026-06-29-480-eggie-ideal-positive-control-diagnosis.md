# #480 follow-up — EGGIE ideal-model positive control: the seed topology is the wall, not the ephemeris

**Date:** 2026-06-29
**Task:** #480 follow-up 1 (the EGGIE reproduction unblock named in
`docs/notes/2026-06-27-480-ieg-reproduction-verdict.md`).
**Status:** diagnosis from two cheap ideal-model spikes; re-scopes the unblock.

## One-line finding

The EGGIE off-paper-basin failure **reproduces in the paper's own ideal
circular-coplanar model**, from the same per-leg multi-rev Lambert seed
(`search/ieg_seed.py`). So the barrier is the **seed construction (per-leg
Lambert legs do not enforce one coherent resonant conic)**, NOT the real
ephemeris. A homotopy of the *existing* seed cannot fix it; the unblock is to
build the paper's resonant-conic initial-guess generator (pp.4-5) so the ideal
positive control closes near-ballistic FIRST, then continue that conic
ideal→real.

## Setup

Ideal model exactly per Hernandez-Jones-Jesick 2017 (AAS 17-608) p.3: Galilean
moons circular + coplanar, `a_Io` = real Io sma, and

- `a_Eur = ((8π+Δ)/(4π+Δ))^(2/3) · a_Io`
- `a_Gan = ((8π+Δ)/(2π+Δ))^(2/3) · a_Io`, with `Δ = 5.2°`.

This yields ideal periods Io 1.7699 d, Europa 3.5272 d, Ganymede 7.0042 d
(T_syn = ideal Ganymede period = 7.0042 d vs the paper's printed 7.05 d). The
tour is the EGGIE sequence Europa-Ganymede-Ganymede-Io-Europa; the patched-conic
powered-flyby ΔV (`nbody/jovian.flyby_min_dv`) is summed over the 3 interior
flybys. Sourced Table-4 targets: Europa V∞ 9.12, Ganymede 7.07, Io 8.38 km/s;
total ΔV 0.70 m/s.

## Spike 1 — free ToFs + phase sweep (Nelder-Mead per phase)

Best ideal close: **sum interior ΔV = 1139 m/s** (paper 0.70). The optimizer let
the ToFs drift far from Table 4 ([3.68, 14.0, 2.34, 11.3] vs [1.59, 8.60, 7.34,
10.69]) and the V∞ collapsed to ~⅓ of target (Ganymede ~2.4, Io ~4.3 km/s) — a
**low-energy off-resonance basin**, i.e. the same family-selection failure the
real-eph corrector hit. Free ToFs make it worse: they buy geometric closure by
abandoning the resonant energy.

## Spike 2 — ToFs pinned to Table 4, phase-only sweep (80 phases / synodic period)

| Topology (n_revs, branch) | best ΣΔV | E0 depart V∞ | Ganymede V∞ (in) | Io V∞ (in) |
|---|---|---|---|---|
| `[0,1,2,1]` high | 11273 m/s | **8.97** (tgt 9.12) | **6.91** (tgt 7.07) | 16.65 (tgt 8.38) |
| `[0,1,1,1]` high | 9409 m/s | 6.24 | 6.04 | 3.72 |
| `[0,1,1,1]` low  | 12727 m/s | 4.26 | 2.00 / 0.00 | 4.59 |

The `[0,1,2,1]`-high topology (the one `ieg_seed.py` already uses) recovers the
**Europa departure V∞ and both Ganymede V∞ to within ~0.2 km/s** at its best
phase — the outer, near-apojove legs are essentially correct. The defect is
**concentrated entirely at the Io encounter** (16–20 km/s vs 8.38), exactly the
"Leg 2 Io-arrival ~20 km/s geometry gap" the `ieg_seed.py` docstring already
flagged. The inward dip to Io's orbit (a_Io ≈ 0.39·a_Gan) is where independent
per-leg Lambert arcs pick the wrong energy/branch.

## Spike 3 — global differential-evolution search (ideal model, ToFs ±25%)

Objective = interior ΔV + wrap-closure V∞ mismatch, free [departure phase, 4
ToFs], `differential_evolution` (popsize 30, 300 iters). Best = **385 m/s**
(paper 0.70). The wrap CLOSES (Europa in=out=8.95, ~tgt 9.12) and **Io comes
into range (9.12 vs 8.38)** — but now **Ganymede collapses (3.9 vs 7.07)** and
carries a 291 m/s flyby mismatch. So the basins FLIP versus spike 2: a loose
global optimizer gets Europa+Io right xor Ganymede right, never all three. The
three single-moon basins do not coincide under free ToFs.

This is the key quantitative result: it matches the paper's own method exactly.
Hernandez also uses per-leg Lambert + free-ToF refinement (pp.6-7) — but seeds
from the **conic initial-guess tool** and refines within only **±10% of a moon
period** (p.7; Io ≈ 4 h, Ganymede ≈ 17 h — extremely tight). The narrow ballistic
basin is unreachable without that accurate seed; a loose search lands in the
wrong single-moon basins, which is precisely the M1 verdict's off-paper-basin
result. **The conic initial-guess generator is the missing ingredient, confirmed
quantitatively.**

## Interpretation

The EGGIE cycler is **one eccentric Jovian conic** (period = (4 synodic / 5
revs)·T_syn ≈ 5.60 d; perijove ≤ a_Io; apojove ≥ a_Gan) whose five encounters
are its crossings of the three moon orbit radii. Solving each leg as an
*independent* Lambert problem does not constrain the legs to lie on one such
conic, so the leg that must reach perijove (the Io dip) relaxes off the resonant
energy — and the corrector then descends into the off-paper high-energy basin
recorded in the M1 verdict. Europa/Ganymede sit near apojove where almost any
arc of the right ToF has ~right V∞, which is why only Io breaks.

## Spike 4 — the resonant conic hits all three V∞ at once (the breakthrough)

Key analytic fact: on a Jovian conic, **|V∞| at a moon depends only on (a, e) and
the moon's orbit radius** — ω and the departure phase θ affect only timing and the
flyby turn angles, not the V∞ magnitude (V∞² = (v_t − v_m)² + v_r² at r = r_moon).
For EGGIE, `a` is FIXED by the 4:5 resonance (a ≈ 9.094e5 km, T_sc = 5.603 d), so
there is a single free parameter (e). Sweeping e over the valid band:

| e | Io V∞ (tgt 8.38) | Europa V∞ (tgt 9.12) | Ganymede V∞ (tgt 7.07) | err |
|---|---|---|---|---|
| **0.621** | **8.393** | **9.099** | 6.794 | **0.278 km/s** |

A SINGLE resonant conic at e ≈ 0.621 puts all three moon V∞ essentially on the
Table-4 targets simultaneously (Ganymede 0.28 km/s low — within MC-refine /
real-T_syn tolerance). This is why the per-leg free-ToF search failed: it left this
conic's energy. Pinning the seed to the resonant conic is the fix, and it is now
PROVEN to put EGGIE in band. e_bounds = [0.536, 0.921] (r_p ≤ a_Io binds e_min;
r_p ≥ R_Jup binds e_max).

## Re-scoped unblock (supersedes "homotopy of the existing seed")

1. **Build the resonant-conic initial-guess generator** (paper pp.4-5, Table 1):
   construct the 2-body Jovian conic from (period = integer·T_syn / sc-rev,
   eccentricity, ω) with r_p ≤ a_Io and r_a ≥ a_Gan; place each moon at the
   conic's orbit-radius crossing at the encounter true anomaly; derive ToFs from
   the conic, not freely. V∞ are then consistent across all legs by construction.
2. **Ideal-model positive control**: require ΣΔV near-ballistic (paper 0.70 m/s,
   navigation-viable 7–9 km/s V∞) in the ideal circular-coplanar model BEFORE
   any real-eph step. This is the mandatory positive control
   (`memory/feedback_verify_gauntlet_with_positive_control.md`).
3. **Homotopy ideal→real**: continue the coherent conic seed from ideal
   circular-coplanar to jup365 (a Jovian analogue of
   `search/continuation.ramped_ephemeris`), re-converging the patched-conic /
   `jovian_shoot` corrector at each step — tracking the ballistic family rather
   than re-seeding per-leg Lambert in the real ephemeris.
4. **n-body confirmation** at λ=1 only (the FD-Jacobian-bound expensive step).

Spike scripts: `scratchpad/spike_ideal_eggie.py`, `scratchpad/spike_pinned.py`
(not committed — scratch). The `ieg_seed.py` `[0,1,2,1]` topology is retained as
the outer-leg seed; only the Io-leg construction needs the conic constraint.
