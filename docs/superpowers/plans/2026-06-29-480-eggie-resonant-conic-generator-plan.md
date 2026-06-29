# #480 follow-up 1 — EGGIE resonant-conic generator + ideal→real homotopy (plan)

**Goal:** reproduce the Hernandez-Jones-Jesick 2017 (AAS 17-608) EGGIE Io-Europa-
Ganymede triple cycler — first near-ballistic in the paper's ideal model (the
mandatory positive control), then continued to real jup365 ephemeris and confirmed
in the Jupiter-central n-body model.

**Why this plan exists:** the M1 verdict
(`docs/notes/2026-06-27-480-ieg-reproduction-verdict.md`) and the diagnosis
(`docs/notes/2026-06-29-480-eggie-ideal-positive-control-diagnosis.md`) showed the
existing per-leg Lambert seed lands in an off-paper basin — even in the ideal model
(spikes: best 385 m/s vs paper 0.70; basins for the three moons don't coincide under
free ToFs). The missing ingredient is the paper's **conic initial-guess tool** that
seeds the optimizer inside the narrow ballistic basin.

## Sourced algorithm (paper pp.4-7, Eqs 1-7; Table 1)

**Conic initial-guess tool**
1. Semi-major axis: `n_rev · T_sc = n_syn · T_syn` (Eq.1), `T_sc = 2π√(a³/μ)`.
   EGGIE: `n_syn=4, n_rev=5` (Table 1 "4:5") → `T_sc = 0.8·T_syn`, a ≈ 9.2e5 km.
2. Eccentricity bounds (orbit must intersect all three moon orbits):
   - `r_pmin = R_Jup`, `r_pmax = a_Io`, `r_amin = a_Gan`, `r_amax = ∞`
   - `e_min = max(1 − a_Io/a, a_Gan/a − 1)`
   - `e_max = min(1 − R_Jup/a, ∞) = 1 − R_Jup/a`
3. Argument of periapsis ω: for a Europa departure phase θ ∈ [0,2π] and chosen
   (a,e), the departure true anomaly ν_dep solves `a_Eur = a(1−e²)/(1+e cos ν)`;
   two prograde branches (outbound `+ν_dep`, inbound `−ν_dep`); `ω = θ − (±ν_dep)`.
4. Loop the conic over `n_rev` revolutions; at each crossing of a moon-orbit radius
   record (moon, time, position). The crossing sequence predetermines the flyby
   sequence + approximate ToFs (the "beauty": no combinatorial sequence search).
   Filter to the target sequence E-G-G-I-E.

**Lambert search:** connect adjacent encounters via zero-SOI Lambert; rev count
fixed by ToF; fast/slow arc chosen to minimise flyby ΔV.

**Flyby evaluation (Eqs 3-7):** tangential periapsis maneuver for the v∞ mismatch.
`δ = angle(v∞⁻, v∞⁺)` (Eq.3); r_p iterative (Eq.4); v_p± (Eq.5); B-plane (Eq.6);
θ_B (Eq.7). `nbody/jovian.flyby_min_dv` already implements the turn-bounded
defect equivalent — reuse it; only add the explicit r_p(δ) solve if needed.

**Monte-Carlo refine:** optimise departure phase + (N−1) ToFs, total pinned to
`n_syn·T_syn`, ToF perturbation ≈ 10% of the flyby body period (Io ~4 h, Ganymede
~17 h — TIGHT), flyby altitudes 25 km–70,000 km. Thousands of samples.

**Ideal model (p.3):** moons circular+coplanar; `a_Io` = real Io sma;
`a_Eur=((8π+Δ)/(4π+Δ))^(2/3)·a_Io`, `a_Gan=((8π+Δ)/(2π+Δ))^(2/3)·a_Io`, Δ=5.2°.
In the ideal model the moon PHASES are free ICs set so each moon sits at its
conic crossing at the encounter time (repeated-moon encounters are the resonance
self-consistency constraint).

## Build stages (TDD; each gated)

**Stage 1 — resonant-conic generator (`search/resonant_conic.py`).**
- `resonant_sma(n_syn, n_rev, t_syn, mu)`; `ecc_bounds(a, r_pmin, r_pmax, r_amin)`.
- conic state `(a,e,ω,ν) → (r,v)`; crossing finder `ν: r(ν)=r_moon`.
- `eggie_initial_guess(...) → (sequence, epochs, moon_phases, sc_node_states)`.
- **GATE (positive control):** feeding the guess through per-leg Lambert + tight
  (±10% body-period) MC refine in the IDEAL model yields ΣΔV near-ballistic AND
  all three V∞ in band (Europa 9.12, Ganymede 7.07, Io 8.38 km/s ± ~0.5). Target:
  ΣΔV ≲ a few m/s (paper 0.70). Golden values SOURCED to Table 4. If this gate
  does NOT pass, STOP and report — the construction is wrong, not the ephemeris.

**Stage 2 — ideal→real homotopy (`search/jovian_continuation.py`).**
- Jovian analogue of `search/continuation.ramped_ephemeris`: ramp moon (a,e,i,phase)
  from ideal circular-coplanar (λ=0) to jup365 (λ=1); re-run Lambert+MC refine at
  each λ, seeded from the previous λ. **GATE:** λ=1 patched-conic ΔV stays in-band
  (paper EIGE real-eph: ballistic first cycle → ~30 m/s over 10 cycles).

**Stage 3 — n-body confirmation (reuse `nbody/jovian.jovian_shoot`).**
- Feed the λ=1 family-correct seed (NOT the off-basin Lambert-real seed) to
  `jovian_shoot`. **GATE:** converges in-family (V∞ near Table 4, correction ΔV
  ≪ the 5.9 km/s off-basin result). Expensive (FD-Jacobian ~20 s/eval) → detached,
  checkpointed runlog, instrumented (`memory/feedback_incremental_progress_reports`).

**Stage 4 — golden + verdict.**
- Un-skip `tests/verify/test_ieg_reproduction_golden.py::
  test_eggie_reproduction_matches_published_invariants` WITHOUT loosening tolerances
  (only if Stage 3 converges in-family). Else: upgrade the characterized negative
  with the sharper "narrow-basin / conic-seed-required" finding.
- NO catalogue self-admission (a reproduced published tour is V4-ceiling, human-
  admitted only). Update OUTSTANDING.md.

## Discipline
- Positive control before trusting any negative (`feedback_verify_gauntlet_with_positive_control`).
- Goldens sourced to Table 4, never to our own code (`feedback_golden_tests_sourced_only`).
- `uv run pytest tests/data tests/search -q` + ruff before each commit; work on main
  (`feedback_never_branch`); incremental pathspec commits.
