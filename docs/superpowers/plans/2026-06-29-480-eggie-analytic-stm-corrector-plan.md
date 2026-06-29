# #480 — EGGIE analytic-STM Jovian corrector (plan to break the FD closure wall)

**Goal:** close the EGGIE ideal-model triple cycler to the continuity floor (and
then real ephemeris), breaking the FD-Jacobian noise plateau characterized in
`docs/notes/2026-06-29-480-eggie-stage2-nbody-verdict.md`.

**Established:** Stage 1 (resonant conic) puts all 3 V∞ on Table-4 (`535d2fb`).
Stage 2 corrector plateaus at ~0.1–0.3 km/s velocity continuity; epoch/ToF freedom
ruled out as the cause; the FD-Jacobian noise floor is the wall (same as
`memory/project_dsm_closure_modeljump_blocker`). The heliocentric shooter already
has an analytic block-bidiagonal STM (`nbody/shooter.py::_stm_jacobian`, opt-in
`jacobian="stm"`) — but it is hardwired to RestrictedNBody (MU_SUN/PLANETS) and
uses REBOUND. The Jovian propagator (`JovianRestrictedNBody`) is REBOUND with a
Python `additional_forces` moon callback, which REBOUND's variational particles do
NOT differentiate (`memory/reference_rebound_variation_custom_force_gotcha`) — so
we CANNOT reuse REBOUND variational STM for the moons.

## Approach: self-contained state+STM co-integrator (avoids the REBOUND gotcha)

Integrate the spacecraft state AND the 6×6 STM together with an analytic
gravity-gradient — no REBOUND variational particles.

**Dynamics** (Jupiter-central, moons as point masses on ephemeris rails r_m(t)):
- a(r,t) = −μ_J r/|r|³ + Σ_m μ_m [ (r_m−r)/|r_m−r|³ − r_m/|r_m|³ ]  (indirect term).
- Gravity gradient G(r,t) = ∂a/∂r = −μ_J/|r|³ (I − 3 r̂r̂ᵀ)
  + Σ_m −μ_m/|d_m|³ (I − 3 d̂_m d̂_mᵀ),  d_m = r_m − r.
- STM ODE: dΦ/dt = [[0, I],[G(r,t), 0]] Φ,  Φ(t0)=I₆. Co-integrate [r,v,Φ] with
  scipy `solve_ivp` (DOP853 or Radau, rtol≈1e-11) or a fixed high-order RK. r_m(t)
  from the injected ephemeris (ideal `IdealJovianEphemeris` or jup365), same model
  the residual uses (so the STM matches the propagated trajectory exactly).
- Surface softening: keep the moon-surface clamp for below-surface dives, but real
  ≥25 km flybys integrate exactly (STM valid there).

## Build stages (TDD, gated)

1. **`nbody/jovian_stm.py`: state+STM co-integrator.** `propagate_with_stm(r0,v0,
   t0,t1,ephem,moons) -> (rf, vf, Phi)`.
   **GATE (mandatory, the gotcha's perturbed-parity test):** for random δ (incl. a
   leg that passes WITHIN a moon SOI / near a flyby), |propagate(x0+δ) −
   propagate(x0) − Φ·δ| / |δ|² is O(δ) (second-order) — i.e. Φ matches the FD
   Jacobian of the REBOUND propagator to ~1e-6 rel AT a flyby, not just in deep
   space. A Sun/Jupiter-only parity test PASSES even with the moon bug, so the gate
   MUST include a perturber-close leg (memory rule).
2. **Analytic block-bidiagonal Jacobian** for the EGGIE multiple-shooting residual:
   ∂c_i/∂node_i = Φ_i, ∂c_i/∂node_{i+1} = −I; flyby-hinge + periodicity-wrap blocks
   (analytic or FD on those few rows). Mirror `_stm_jacobian` layout. Validate the
   assembled Jacobian vs `_fd_jacobian` on the EGGIE seed (parity oracle).
3. **Wire `jacobian="stm"` into `ideal_eggie_shoot`** (pass `jac=` to
   `least_squares`, method `trf`/`lm`). Re-run the EGGIE close.
   **GATE:** does the defect break below the ~3e2 FD plateau toward the continuity
   floor (leg velocity ≪ 0.1 km/s; correction ΔV → small)? Honest numbers; do not
   loosen. If it converges near-ballistic → the ideal EGGIE is reproduced.
4. **If ideal closes:** Stage 3 homotopy ideal→real (ramped Jovian ephemeris) +
   real-eph confirm with the STM corrector; then Stage 4 golden/verdict (un-skip
   only on a real converged ballistic close, tolerances unchanged; no catalogue
   self-admission — human-admitted V4-ceiling reproduction).
   **If it still plateaus:** the wall is deeper than FD noise (genuine model /
   local-minimum); record honestly and reassess (sub-arc nodes, gravity homotopy,
   or a global/family method).

## Discipline
- Perturbed-parity gate is NON-NEGOTIABLE before trusting the STM (the gotcha makes
  a wrong STM look right in deep space). Keep FD as the parity oracle.
- Goldens sourced to Table 4; work on main; incremental pathspec commits; ruff+mypy
  pre-commit; `uv run pytest tests/nbody tests/search -q -m "not slow"` green.
- Long runs: instrument (runlog append+flush), bound trf max_nfev, cap wall.
