# #411 — Cross-system θ-closure is feasible at SINGLE revolution (amplitude knob)

Date: 2026-06-20. Supersedes the fixed-amplitude (41,19) feasibility map
(2026-06-20-411-theta-closure-feasibility.md). Necessary-condition analysis;
asserts NO closed orbit. Reproduce:
`uv run python scripts/analyze_411_amplitude_theta_closure.py`.

## The (41,19) result was an artifact of freezing the amplitudes
The earlier map held both Lyapunov amplitudes fixed and forced integer revolution
counts to absorb the entire ~2.2 rad single-rev θ-time-consistency gap, giving
(n_em=41, n_se=19) → ~11 yr. That is the **wrong model** on two counts.

**Infeasible.** A manifold seeded at ε shadows its orbit for only
≈ ln(ε_range)/ln(|λ_u|) revolutions before exponential departure. Measured
Floquet multipliers:
- EM-L2 Lyapunov: |λ_u| ≈ 1.2e3 (ln ≈ 7.1/rev) → shadow budget **~2–4 rev**
- SE-L2 Lyapunov: nearly neutral off the Canalias bifurcation (|λ_u| ≈ 1.000x)
  → budget effectively **unbounded** (10⁴–10⁶ rev), though it stiffens to ~3 rev
  approaching the bifurcation C.

So n_em=41 is **>10× over the EM shadow budget** — you cannot loiter 41 revs on an
orbit that unstable; the trajectory left long ago. (41,19) is physically void.

## Unnecessary — amplitude is a continuous phase knob
Over each Lyapunov family the per-revolution relative-phase advance Δθ(C) mod 2π
sweeps essentially the whole circle (the SE term is ω_rel·T_se ≈ 38 rad before the
mod, so a *tiny* SE-period change sweeps all of [0,2π)). Measured over the families:

| family | members | Δθ mod 2π span | shadow budget |
|---|---|---|---|
| EM-L2 | 11 | [0.26, 5.88] rad = **89% of 2π** | 2.1–3.9 rev |
| SE-L2 | 12 | [1.52, 5.88] rad = 69% of 2π | 3.2 rev → ∞ |

The θ-closure condition

    gap(c_em, c_se) + n_em·Δθ_em(c_em) + n_se·Δθ_se(c_se) ≡ 0   (mod 2π)

with **n_em = n_se = 1** is ONE equation in TWO continuous amplitude knobs. The
EM family alone sweeps ~full circle *within its 2–4 rev shadow budget*, so for any
fixed c_se, sweeping c_em crosses 0 — a **1-D solution curve exists at single
revolution**. One covering knob within budget is sufficient; the EM side supplies it.

## Verdict — FEASIBLE at single revolution
The θ-closure precondition is met without any multi-revolution loitering. The
productive #411 build is a **coupled single-rev corrector**: free variables
(c_em, c_se) + the two connection legs' phases + θ at n_em=n_se=1; residual = both
inertial position gaps PLUS the θ-time-consistency closure. Solve along the closure
curve and read off the patch ΔV.

## Honesty caveat — necessary, not sufficient
This establishes only that a θ-consistent *curve* exists at low rev. **Sufficiency
is open**: whether both connection legs still converge at *low patch ΔV* on that
curve is exactly what the corrector must determine — the amplitude that closes θ may
not be the amplitude where the manifolds co-reach the patch cheaply. A clean
negative there (θ closes but ΔV stays high everywhere on the curve) is still a real
result. Build the corrector; verify with an independent Radau pass before any
closed-cycle claim.
