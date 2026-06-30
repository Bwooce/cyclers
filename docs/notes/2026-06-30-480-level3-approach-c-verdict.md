# #480 — Approach C (level-3 high-fidelity maintenance): ATTEMPTED, blocked at the conversion

**Date:** 2026-06-30. Status: Approach C (the paper's full level-3 high-fidelity
maintenance — the only remaining #480 lever) was attempted with fresh diagnostics and is a
**characterized NEGATIVE**: the blocker is the patched-conic→continuous-gravity n-body
CONVERSION, re-confirmed three independent ways. A faithful level-3 maintenance number
needs the paper's B-plane SNOPT NLP with a continuation/good-seed strategy (scope Approach
C, "weeks, last resort") — not achievable by our corrector. The #480 maintenance result
stands at level-2.

## What was tried (the maintenance reframe at level-3)
The scope doc's crux: a maintenance ΔV is forward-propagate + per-cycle re-target, NOT a
periodic-orbit closure. So I tested the forward-propagate framing directly at level-3 — for
BOTH cyclers — using `JovianRestrictedNBody` (REBOUND/IAS15, moons-on-rails, flybys
integrated) seeded from the feasible patched-conic members via `periapsis_node`.

## Finding 1 — the patched-conic seed is NOT an n-body-valid trajectory
Forward-propagating the patched-conic nodes through the real n-body drifts enormously:

| cycler | leg | tof (d) | |Δr| (km) | |Δv| (m/s) |
|---|---|---|---|---|
| EGGIE | E→G | 1.63 | 18,554 | 96 |
| EGGIE | G→G | 8.68 | 189,926 | 1,801 |
| EGGIE | G→I | 7.16 | 97,396 | 3,337 |
| EGGIE | I→E | 10.98 | 136,549 | 2,749 |
| EIGE | E→I | 0.46 | 3,439 | **11,608** |
| EIGE | I→G | 1.08 | 1,129,072 | 15,172 |
| EIGE | G→E | 5.54 | 168,688 | 3,488 |

The km/s velocity drifts (and the EIGE E→I leg's 11.6 km/s over 0.45 d) are **uncontrolled
close encounters** — the `periapsis_node` placement (a near-zero-turn node clamped to
~0.6·SOI offset; a weak-gravity Io node) sends the integrated trajectory diving into a moon.
The zero-SOI patched-conic (legs between moon CENTRES, no moon gravity along the leg) is
qualitatively far from any continuous-gravity trajectory. EIGE's shorter legs do NOT help —
the weak-Io / near-zero-turn-Europa seed nodes are the dominant error, not leg length.

## Finding 2 — the corrector that would convert it does not converge
- The Jovian multiple-shooting corrector (`_cycle_residual` / `shoot_cycle`: 31 free vars =
  node velocities + interior node states + epoch offsets, with leg continuity + encounter
  hinges) is the step that would pull the patched-conic seed onto a continuous n-body cycle.
- **FD is compute-infeasible at this scale:** on the clean feasible 2020-09-22 EGGIE member,
  the FD corrector (31-var Jacobian = 32 × 4 n-body propagations per build) did not complete
  60 nfev in 10 minutes (killed). "FD Jacobian dominates" ([[project_dsm_closure_modeljump_blocker]]).
- **The analytic STM does not cure the plateau:** previously documented
  (`2026-06-30-480-eggie-level3-nbody.md` correction, `cc4f241`) — STM build 0.5 s vs FD
  16.6 s, parity-verified, but the residual PLATEAUS at ~0.1-0.4 km/s continuity (2nd chunk
  moved 3.4%). Not an iteration-count problem; the corrector does not converge to a
  continuous n-body periodic orbit.

## Finding 3 — EIGE's level-2 "sub-surface" is intrinsic, not just patched-conic
The real-eph EIGE Io flyby needs a ~98° turn at V∞≈5.3 km/s; Io's weak gravity
(μ≈5960 km³/s²) puts the patched-conic periapsis at ~70 km from Io's CENTRE (sub-surface)
even for that sub-180° turn. The B-plane freedom at level-3 chooses the turn PLANE, but the
turn-angle ceiling at min altitude is fixed by |V∞| and μ — so the near-anti-parallel Io
in/out geometry is hard to realise feasibly regardless of model. (Level-3 could still find a
DIFFERENT feasible EIGE by re-routing through the flyby cone, but only via the full NLP.)

## Conclusion / #480 final standing
- **Approach C is blocked at the patched-conic→n-body conversion** — confirmed by forward-prop
  drift (both cyclers), the FD-corrector compute wall, and the documented analytic-STM plateau.
  A faithful level-3 maintenance number requires the paper's full B-plane SNOPT NLP with a
  continuation/good-seed strategy (genuinely different + heavier method; weeks). NOT obtainable
  by re-running our multiple-shooting corrector.
- **The #480 maintenance gap is closed at LEVEL-2** (the validated lane): method validated on
  Liang Member D (#223); EGGIE maintenance characterized (ballistic for exactly 2 cycles, then
  large impulses, robust — `2026-06-30-480-eggie-maintenance-verdict.md`); EIGE construction +
  the 1-rev B-plane wall characterized. No exact-paper-number (EIGE ~30 m/s/10) reproduction;
  no catalogue impact.
- **What a faithful level-3 completion needs** (if ever pursued): (1) a continuation/good-seed
  strategy to land an n-body-valid one-cycle trajectory (the periapsis_node seed is too far);
  (2) the analytic-STM corrector generalised to arbitrary sequence with B-plane + altitude
  constraints; (3) the full ΔV-minimising NLP (FBS gradients + SLSQP, the paper used SNOPT).
  This is the scope's Approach C, correctly rated "weeks, last resort".
