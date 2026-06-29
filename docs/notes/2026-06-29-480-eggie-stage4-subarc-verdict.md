# #480 Stage 4 verdict — sub-arc multiple shooting does NOT break the wall; the ~0.1 km/s floor is robust

**Date:** 2026-06-29
**Task:** #480 follow-up, Stage 4 (sub-arc multiple shooting) — the last gradient/
discretisation lever of the corrector program.
**Plan:** `docs/superpowers/plans/2026-06-29-480-eggie-analytic-stm-corrector-plan.md`
**Builds on:** Stage 1 conic (`535d2fb`), STM co-integrator (`d619c44`), Jacobian
(`9f98bb1`), `jacobian="stm"` wiring (`193569c`), Stage-4 sub-arc residual+Jacobian
(`428ca31`). **Reads with** the Stage-2/3 verdict notes.

## One-line verdict — NULL (sub-arc does not close EGGIE); program-level characterized negative

Adding K interior shooting nodes per leg (K=2,4,6; M=9,17,25 nodes) lowers the total
residual (4.35e2 → 1.40e2) and the max leg *velocity* continuity (0.41 → 0.058 km/s)
**only slowly**, and does so by letting the corrector **drift the encounter V∞ off
the Table-4 targets and keep ΔV at 3000–6000 m/s** — i.e. it buys continuity by going
non-physical, not by finding the ballistic cycler. The ~0.05–0.1 km/s
velocity-continuity floor is **robust across all four correctors built this round**
(FD → analytic STM → epoch/ToF-free → sub-arc). EGGIE is **not** reproduced as a
ballistic n-body cycler; the barrier is a genuine basin/model feature localized to
the high-speed Io perijove, not a solver-gradient or discretisation artifact.

## Stage-4 numbers (ideal model, `scripts/_stm_subarc_480.py 2,4,6 400`, jac="stm")

| K (sub-arcs/leg) | nodes M | final \|r\| | leg pos defect | leg vel defect | wrap vel | corr ΔV | V∞ Europa (tgt 9.12) |
|---|---|---|---|---|---|---|---|
| 1 (Stage 3) | 5 | 1.99e2 | 0.5–4.1 km | 0.057–0.122 km/s | 0.057 | 1450 m/s | 9.26 |
| 2 | 9 | 4.35e2 | 0.27–6.3 km | 0.0096–**0.41** km/s | 0.080 | 3167 m/s | 9.06 |
| 4 | 17 | 1.83e2 | 0.05–4.3 km | 0.0039–**0.0998** km/s | 0.058 | 6136 m/s | 7.50 (drift) |
| 6 | 25 | 1.40e2 | 0.05–2.5 km | 0.0040–**0.058** km/s | 0.039 | 4908 m/s | 8.08 (drift) |

(Jacobian-vs-FD parity for the sub-arc Jacobian passed before these runs — `428ca31`.)

Observations:
- The **max** leg velocity defect falls with K (0.41 → 0.0998 → 0.058 km/s) — sub-arcs
  do spread the perijove sensitivity — but it does NOT collapse toward the 1e-6 km/s
  floor; it crawls. At K=6 it is still ~0.058 km/s ≈ 58 m/s, ~80× the paper's 0.70 m/s.
- As K grows the corrector **relaxes the encounter V∞ off Table 4** (Europa 9.06 → 7.50
  → 8.08) and keeps ΔV high (3–6 km/s). The extra interior freedom is spent reducing
  the total-norm by deforming the trajectory away from the physical resonant family,
  not by closing it ballistically — the signature of a degenerate / non-unique least-
  squares minimum near the seed, not convergence to the cycler.
- At K=2 the V∞ stay near target (9.06/7.09/7.14/8.68/9.01) but one sub-arc holds
  0.41 km/s — the stubborn region is the **Io perijove** (the fastest, deepest pass,
  where the zero-SOI patched conic and continuous-gravity n-body diverge most).

## Program-level conclusion (the four levers)

| lever | result | leg-vel floor |
|---|---|---|
| FD-Jacobian (Stage 2) | plateau | 0.10–0.19 km/s |
| epoch/ToF freedom (Stage 2b) | no help | ~0.1–0.3 km/s |
| analytic STM (Stage 3) | broke FD plateau; V∞ in-band, ΔV ⅓ | 0.057–0.122 km/s |
| sub-arc nodes (Stage 4) | total↓, but V∞ drift + high ΔV | 0.058 km/s (K=6) |

The exact analytic Jacobian was a real win (FD noise WAS a wall — position 10×, V∞
in-band, ΔV to ⅓). But beneath it sits a ~0.06–0.1 km/s velocity-continuity floor
that NONE of gradient-quality, epoch freedom, or discretisation refinement breaks.
That robustness is the evidence: this is a **basin/model property**, not a numerics
problem. The local multiple-shooting family around the resonant-conic seed does not
contain a ballistic continuous-gravity EGGIE; either the Io-perijove patched-conic↔
n-body gap genuinely carries ~tens of m/s in this model, or the true ballistic cycler
sits in a different sub-basin not reachable by local correction from this seed.

## Status of #480

- **Stage 1 (resonant conic): SUCCESS** — all three V∞ on Table 4 (seed/basin solved).
- **Corrector program (Stages 2–4): infrastructure SHIPPED** — analytic state+STM
  Jovian co-integrator (parity-gated at a real flyby), block-bidiagonal Jacobian
  (~40× faster than FD, FD-parity 1e-6), opt-in `jacobian="stm"`, and sub-arc
  multiple shooting. All reusable for any future Jovian moon-tour.
- **EGGIE ballistic reproduction: characterized NEGATIVE (advanced).** From M1's
  off-basin 5.9 km/s to a right-basin solution with V∞ within ~0.1–0.3 km/s of Table 4
  and a robust ~0.06–0.1 km/s velocity-continuity wall. NOT a ballistic cycler. No
  catalogue change; golden (`tests/verify/test_ieg_reproduction_golden.py`) stays
  skipped.

## Remaining lever (one, distinct in kind) + honest prior

The only untried lever is **gravity homotopy** (ramp moon GM 0→1, tracking the family
from the trivially-closed 2-body case). Honest prior: it is *continuation*, which
helps the corrector TRACK a family and avoid local minima — but the GM=1 endpoint is
the same model whose inherent ~0.1 km/s Io-perijove residual the four direct correctors
already hit, so homotopy changes the *path*, not the GM=1 model's floor. Lower expected
leverage than the STM was. Alternatives if pursued: the paper's exact zero-SOI
instantaneous-flyby patched-conic model (matches their 0.70 m/s definition; our n-body
integrates continuous gravity, a genuinely different — and harder — model), or a
global/family-targeted method. Scratch: `scripts/_stm_*_480.py` (removed).
