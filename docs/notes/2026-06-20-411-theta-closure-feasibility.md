# #411 — Cross-system θ-closure feasibility (patched-CR3BP)

Date: 2026-06-20. Gates the full multi-revolution closure corrector for the #405
cross-system (Sun-Earth ↔ Earth-Moon) heteroclinic cycle. Necessary-condition analysis;
asserts NO orbit. Reproduce: `uv run python scripts/analyze_405_theta_closure.py`.

## Context
After the #411 direction-bug fix, BOTH cross-system connection legs converge at low ΔV
(forward EM-L2→SE-L2 ≈0.36 km/s; return SE-L2→EM-L2 ≈0.15 km/s). The open gate for a
CLOSED cycle is θ-time-consistency: the legs converge at relative phases that, accounting
for the ~52-day forward transit's phase advance, leave a **~2.2 rad single-revolution
gap**. A closed cycle requires absorbing that gap via multi-revolution phasing (extra
orbit revolutions each advance the SE-EM relative phase) — the patched-model analog of the
~19 yr Metonic (235:19) commensurate return predicted in note 2026-06-17-316.

## Method
ω_rel = 1/t_s(EM) − 1/t_s(SE) = 2.466e-6 rad/s (relative SE-EM line rate). Each extra
revolution on an orbit advances θ by ω_rel·T_orbit:
- EM-L2 period 14.85 d → **Δθ_em = 3.165 rad/rev**
- SE-L2 period 0.489 yr → **Δθ_se = 0.330 rad/rev**

Grid-search integers 0 ≤ n_em, n_se ≤ 60 for the minimal mod-2π residual of
`gap + n_em·Δθ_em + n_se·Δθ_se` (pure arithmetic; helper `theta_commensurability`,
unit-tested).

## Result — FEASIBLE
Best lattice point **(n_em = 41, n_se = 19)** nulls the 2.2 rad gap to a residual of
**2e-4 rad** (well inside any achievable phasing tolerance). Resulting cycle duration
**≈ 10.95 years** — the same order as the #316 Metonic prediction (a near-commensurate
SE-EM phasing). So a θ-commensurate closed cross-system cycle is NOT excluded by the
phase arithmetic; it is admissible at ~11 yr with 41 EM-orbit + 19 SE-orbit revolutions.

## Honesty caveat (necessary, not sufficient)
This is a NECESSARY-condition feasibility check on the relative-phase arithmetic only. It
does NOT assert a verified closed orbit. Sufficiency is unproven: the connection legs were
converged at SINGLE revolution; whether they re-converge with 41/19 orbit loops before
departure (a multi-revolution manifold-departure corrector) is the open question. 41/19
revolutions is also a lot — the manifold would shadow each orbit many times, where
accumulated instability may degrade the connection.

## Conclusion → next step
The phasing permits closure at ~11 yr, so the FULL multi-revolution θ-closure corrector
(#411) is justified to build: add n-revolution manifold departure, jointly solve both legs
+ θ-time-consistency at (n_em, n_se) near (41, 19), and rigorously verify (independent
Radau; the cycle actually returns) before any closed-cycle claim. The alternative venue is
#412 BCR4BP Phase B, where closure is well-posed by construction.
