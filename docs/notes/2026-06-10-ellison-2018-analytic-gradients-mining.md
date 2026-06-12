# Ellison et al. 2018 — analytic match-point gradients for bounded-impulse shooting (mining)

**Date:** 2026-06-10
**Source:** D. H. Ellison, B. A. Conway, J. A. Englander, M. T. Ozimek, "Analytic
Gradient Computation for Bounded-Impulse Trajectory Models Using Two-Sided
Shooting," *Journal of Guidance, Control, and Dynamics*, Vol. 41, No. 7, 2018,
pp. 1449–1462. doi:10.2514/1.G003077. Read from the NASA Public Access author
manuscript (41 pp.); **page numbers below are author-manuscript pages, equation/
table/figure numbers are publication-stable** — cite the eq. numbers, not the pages.
**Writeback: NONE.** Extraction + adoption assessment only; no code or catalogue
edits in this pass.

**TRANSCRIPTION RESCAN (2026-06-12, spot-check):** the cited headline numbers
re-read against the PDF — **MATCH**. Table 4 (p. 41) verbatim: best mass delivered
**2104 kg**, mean feasible solutions **55**, success rate **94%**, mean time to
best **16.6 min**, finite-differenced column all dashes (zero feasible in 100
runs); corroborated independently by the Sec. V prose (p. 20: "None of the runs
using finite differencing identified any feasible solutions"; "94 identified a
feasible solution"; "55 feasible solutions in its 2 h search"; "best solution
after 16.6 min"; FD seeded with the analytic optimum declared infeasible and
diverged). Fig. 8 (p. 37) read: FD/analytic ≈160 at 10 segments rising to ≈1650
at 200 segments (AD/analytic ≈25→≈730) — the §6 figure-read caveat stands. Fig. 6
(p. 35) event labels and Table 3 (p. 40) inputs also match §6 as transcribed.

**Verdict: USEFUL.** This is the published algorithm specification for exactly the
thing flagged as our known performance bottleneck: replacing finite-difference
Jacobians in STM-based multiple shooting with analytic STM/MTM chain products. It
is the gradient layer of the EMTG lineage (Englander & Conway [46]) whose MBH we
already implement (`search/mbh.py`, Englander 2014 tuning = their ref [64]). The
flyby-continuity gradients (Appendix, Eqs. A1–A6) and the impulsive-maneuver
transition matrix (Eqs. 29–30, 41–42) are directly implementable. The published
test case (67P sample return, Tables 3–4, Fig. 6) is an end-to-end optimizer
benchmark, **not** a clean unit-level golden — see §6.

---

## 1. Two-sided (forward–backward) shooting structure (Sec. II.A–B, pp. 3–4)

- A mission is `N_p` phases; phase boundaries are **control points** (massive
  bodies — planets, moons, asteroids — or free points in space) (Sec. II.A).
- Decision variables at each control point: spacecraft **mass** and **v∞ vector
  relative to the body** on either side of the phase (Table 1: t_launch, v∞
  magnitude, RA/DEC of launch asymptote, Δt_p per phase, m_f per phase,
  v∞_p(t0) 3(N_p−1), v∞_p(tf) 3N_p). Combined with ephemeris data these form the
  full state at both phase ends, which is propagated **forward from the left
  boundary and backward from the right boundary** (Fig. 1, p. 30).
- State vector (Eq. 1): `X = [r; v; m]`.
- The hallmark constraint is the **match-point defect** (Eq. 2):

  ```
  c_mp = X^B_mp − X^F_mp = [r^B−r^F; v^B−v^F; m^B−m^F] = 0
  ```

  The match point sits at mid-phase (between segments N/2 and N/2+1 in the MGALT
  discretization, Fig. 2 p. 31).
- **Flyby continuity** at a phase boundary ending in a flyby — two constraints:
  - v∞ magnitude continuity (Eq. 3): `c_v∞ = v∞⁺ − v∞⁻ = 0` (patched-conic flyby
    modeled as a velocity-vector discontinuity; the *direction* change is free,
    the magnitude is matched).
  - Minimum-altitude / realizable-turn (Eq. 4):

    ```
    c_flyby-alt = r_periapse − (r_body + h_safe)
                = (μ_body/v∞⁺²)·[1/sin(δ/2) − 1] − (r_body + h_safe) ≥ 0
    ```

    with turn angle (Eq. 5): `δ = acos( v∞⁻·v∞⁺ / (v∞⁻ v∞⁺) )`.
- Total-TOF box constraint (Eq. 6): `Δt_min ≤ Σ_i Δt_p_i ≤ Δt_max`.
- **What it buys numerically** (Sec. II.A + III): (a) each STM/MTM chain runs only
  half a phase, halving the sensitivity growth of pure forward shooting; (b) both
  phase endpoints are pinned to bodies through v∞ decision variables, so flyby
  constraints (Eqs. 3–5) are *algebraic in the decision variables* — no
  propagation enters their gradients (Appendix); (c) **no Lambert solver
  anywhere** — interior-maneuver Δv's are decision variables and continuity is
  enforced by the defect, so the transcription never needs Lambert partials.
  Sec. V (p. 19) uses exactly our stack on top: MBH (Yam et al.; Englander &
  Englander tuning) + SNOPT as local NLP.

## 2. Bounded-impulse transcriptions (Sec. II.C–D, pp. 4–6)

Two specializations of the FBS phase:

### MGALT (Sims–Flanagan low-thrust; Sec. II.C)
Phase split into `N` equal-time segments; continuous thrust per segment
approximated by one bounded impulse at segment center, with bound (Eq. 7):
`Δv_max_k = N_active · D · T_max_k · (t_f−t_0) / (m_k N)` (D = duty cycle).
Impulse insertion (Eqs. 8–9): `v_k⁺ = v_k⁻ + Δv_max_k u_k` (forward) /
`v_k⁻ = v_k⁺ + Δv_max_k u_k` (backward), with up-to-unit control `‖u_k‖ ≤ 1`
(Eq. 10). Mass recursion (Eq. 11): `m_k = m_{k−1} − ‖u_{k−1}‖ Δm_max_{k−1}`
(forward), `Δm_max = D Δt_k ṁ_max_k`. Thrust/flow-rate depend on the solar-power
model `T_max(P(r_s/⊙, t))` (Eq. 12); power model Eq. (73) p. 20.

### MGAnDSMs (n deep-space maneuvers, chemical; Sec. II.D) — **our DSM lane's cousin**
- Each phase carries up to `n` impulsive midcourse maneuvers; the Δv at each is
  **directly a decision variable** (Eq. 13): `v_k⁺ = v_k⁻ + Δv_k`.
- Mass by rocket equation (Eq. 14): `m_{k+1} = m_k e^(−Δv_k/v_e)`.
- Maneuver placement: inter-maneuver times Δt_k with `Σ_{k=1}^{n+1} Δt_k − Δt_p = 0`
  (Eq. 15); in practice the decision variables are **burn-index fractions
  α_k ∈ [0,1] with Δt_k = α_k Δt_p** (p. 6). Table 2 (p. 39): per phase, decision
  variables are Δt_1 (time to first maneuver, 1), Δt_2…Δt_n (n−1), Δt_{n+1} (1),
  and Δv_1…Δv_n (3n).
- **Mapping to our Takao η-coordinate** (`search/dsm_leg.py`): our η is exactly
  their α_1 for n=1; the difference is that we Lambert-solve the post-DSM arc
  (one-sided shooting, Δv emerges from the Lambert mismatch) whereas MGAnDSMs
  propagates both ends and lets the NLP close the mid-phase defect with Δv as an
  explicit variable. Their form has 6–7 extra constraints per phase but fully
  analytic derivatives; ours has fewer variables but a Lambert solve whose
  partials this paper does **not** provide.

## 3. Match-point partial derivatives — the core formulas (Sec. III, pp. 7–18)

Master expression (Eq. 16): `∂c_mp/∂x = ∂X^B_mp/∂x − ∂X^F_mp/∂x`.
Derivatives are propagated to the match point by alternating two matrix types
(Fig. 4, p. 33): **STMs across Keplerian coast arcs** and **maneuver transition
matrices (MTMs) across the velocity discontinuities** (r ∈ C⁰, v ∈ C⁻¹ at each
impulse).

1. **Two-body STM** (Eq. 17): `Φ(t_k, t_{k−1}) = [R̃ R; Ṽ V]` in Battin's
   notation — the 3×3 quadrants computed analytically (their refs [51–56]:
   Battin; Shepperd 1985 universal-variables STM; Goodyear; Sperling; Herrick;
   Lemmon & Brooks).
2. **Augmented state** (Eq. 18): `X = [r; v; m; Δt_p; Δt_previous]` (9×1) — two
   time slots so TOF sensitivities ride the same chain. Augmented STM (Eq. 19,
   9×9): top-left the 6×6 two-body STM; a 6×1 zero column for mass; a 6×1 column
   `[∂r_k⁻/∂Δt_p; ∂v_k⁻/∂Δt_p]`; a 6×1 zero column for Δt_previous; bottom rows
   `[0_{3×6} | I_{3×3}]`. Backward half-phase convention (Eq. 20):
   `Φ_k = ∂X_k⁺/∂X_{k+1}⁻`. Note (p. 9): the augmentation is bookkeeping only —
   c_mp stays a defect in r, v, m, and the Δt_previous slot's computations are
   valid for *all* previous phase TOFs and the launch epoch.
3. **MGALT MTM** (Eq. 21, 9×9): identity in r; `∂v_k⁺/∂r_k⁻` (Eq. 22, the
   power-model coupling `u_k (D Δt_k/m_k)(∂T_max/∂P)(∂P/∂r_s/⊙)(∂r_s/⊙/∂r)`,
   with `r_s/⊙ = r − r_⊙` Eq. 23 making it central-body agnostic; backward
   variant Eq. 24 adds the ṁ_max term); `∂v_k±/∂m = ∓u_k Δv_max_k/m_k`
   (Eqs. 25–26); mass–position couplings Eqs. 27–28; TOF couplings M_k24/M_k34
   (Eqs. 47–48, 52) and epoch couplings M_k25/M_k35 (Eqs. 62–64).
4. **MGAnDSMs MTM** (Eq. 29) — *the impulsive case, trivially simple*: identity
   everywhere except `M_k33 = ∂m_{k+1}/∂m_k = e^(−Δv_k/v_e)` (Eq. 30). The Δv
   decision-variable partials: `∂m_{k+1}/∂Δv_k = −(m_k/v_e)(Δv_kᵀ/Δv_k)
   e^(−Δv_k/v_e)` (Eq. 41) and connection matrix `Ξ_k = ∂X_k⁺/∂Δv_k =
   [0_{3×3}; I_{3×3}; ∂m_{k+1}/∂u_k; 0_{2×3}]` (Eq. 42).
5. **STM–MTM chain** (Eqs. 31–32): sensitivity of the forward match-point state
   to the state just after impulse k:

   ```
   ∂X^F_mp/∂X_k⁺ = Φ_{N/2+1} M_{N/2} Φ_{N/2} · … · M_{k+1} Φ_{k+1}      (31)
   ```

   and a match-point Jacobian column for decision variable x_i is that chain
   right-multiplied by the **derivative connection matrix** Ξ_k = ∂X_k/∂x_i
   (Eq. 32). MGALT control-variable cases: Eqs. 33–36 (forward: Ξ_k from
   `∂v_k⁺/∂u_k = Δv_max_k I_{3×3}` Eq. 33 and `∂m_{k+1}/∂u_k =
   −(u_kᵀ/‖u_k‖) D Δt_k ṁ_max_k` Eq. 34); backward: Eqs. 37–40 (note Eq. 37's
   `u_k u_kᵀ` rank-1 correction).
6. **Phase-TOF partials** (Sec. III.C): the per-arc time derivative is Pitkin's
   Lagrange-coefficient form (their refs [57] Pitkin 1966, [58] Lantoine &
   Russell):

   ```
   ∂r_k⁻/∂Δt_p = ±[Ḟ r_{k−1}⁺ + Ġ v_{k−1}⁺]·(∂Δt_k/∂Δt_p)    (43)
   ∂v_k⁻/∂Δt_p = ±[F̈ r_{k−1}⁺ + G̈ v_{k−1}⁺]·(∂Δt_k/∂Δt_p)    (44)
   ```

   with `∂Δt_k/∂Δt_p = 1/N` (full segments) or `1/(2N)` (the four half-segments)
   for MGALT (Eq. 45), and `∂Δt_k/∂Δt_p = α_k` for MGAnDSMs (Eq. 58). Segment
   epoch sensitivity `∂t_k/∂Δt_p = (k−0.5)/N` forward (Eqs. 49–50), `(0.5−k)/N`
   backward (Eq. 51). Chain assembly Eqs. 53–57 (Ξ_1 seeds the chain with the
   Δt_p column and a 1 in the Δt_p slot; φ^B Eq. 57 handles the moving right
   boundary through `∂X_f/∂Δt_p` = ephemeris velocity of the boundary body).
7. **Previous-phase TOFs / launch epoch** (Sec. III.D, Eqs. 59–61): these shift
   the boundary body along its ephemeris (`∂X_0/∂Δt_previous` in Eq. 60) but do
   **not** enter the Kepler propagation times; `∂t_k/∂Δt_previous = 1` (Eq. 61).
8. **MGAnDSMs burn-index partials** (Sec. III.E, Eqs. 65–68): each α_k stretches
   exactly one coast arc in a half-phase; same chain with Ξ_k built from
   Eqs. 43–44 applied to that arc.

## 4. Flyby-continuity gradients (Appendix, Eqs. A1–A6, pp. 22–24)

Closed-form gradients of the altitude constraint (Eq. 4) w.r.t. the Cartesian
components of v∞⁻ (Eqs. A1–A3, with scalars `α = (v∞⁻·v∞⁺)/(β^{1/2}γ^{1/2})`,
`β = ‖v∞⁺‖²`, `γ = ‖v∞⁻‖²`) and w.r.t. v∞⁺ (Eqs. A4–A6, with `ε = φ/(ψ^{1/2}ξ^{1/2})`,
`ξ = ‖v∞⁺‖²`, `φ = v∞⁻·v∞⁺`, `ψ = ‖v∞⁻‖²`; the v∞⁺ case has the two extra
`2v∞⁺μ/(r_p ξ²)` terms because r_periapse itself depends on ‖v∞⁺‖). Gradient of
Eq. (3) is trivial (unit vectors). These are purely algebraic in the v∞ decision
variables — no STM chain involvement. Directly portable to any of our correctors
that carry flyby turn/altitude feasibility (e.g. the powered-flyby feasibility
checks in the free-return/multi-arc lanes).

## 5. Path (distance) constraint gradients (Sec. IV, Eqs. 69–72, pp. 18–19)

Min/max distance to any body, `d_LB ≤ r_s/c−body ≤ d_UB` (Eq. 69), enforced at
every maneuver point; needs *intermediate* STM–MTM chains `∂X_j⁻/∂X_k⁺`
(Eq. 70 forward, Eq. 71 backward), combined with Ξ_k as
`∂d_j/∂u_k = (∂d_j/∂X_j⁻)(∂X_j⁻/∂X_k⁺) Ξ_k` (Eq. 72). Operational warning
(p. 18): in large problems these chains can be **>50% of total optimization
runtime**; the same machinery generalizes to any path constraint. Relevant if we
ever add perihelion-floor or Earth-comm constraints to cycler legs.

## 6. Published numeric test case (Sec. V–VI, Tables 3–4, Figs. 6–8)

67P/Churyumov–Gerasimenko round-trip sample return (New Frontiers 4 CONDOR-like),
MGALT transcription, MBH + SNOPT.

- **Inputs (Table 3, p. 40):** max initial mass 10,000 kg; earliest launch
  16 June 2024, 1-day window; 67P arrival 1–30 Apr 2029; fixed 67P departure
  8 Dec 2033; Earth return 8 Nov 2034–9 Nov 2036; max flight 12.5 yr; Atlas
  V(431); max C3 25 km²/s²; DLA ∈ [−28.5°, 28.5°]; low-thrust rendezvous; 3×
  BPT-4000 extra-high-Isp, duty cycle 0.9, min-thruster throttle logic; solar
  power γ0=1, γ1..4=0 (1/r² model, Eq. 73), 30 kW BOL at 1 AU; sequence E–67P–E;
  60-day post-launch forced coast; **80 segments/phase**; SNOPT feas. tol 1e-5;
  objective max returned mass; MBH 2 h, SNOPT 120 s / 8000 major iterations.
- **Outputs:** best solution (Fig. 6, p. 35): launch 17 June 2024, C3 = 6.101
  km²/s², DLA = 11.1°, m0 = 4653 kg; 67P rendezvous 28 Apr 2029, m = 2610 kg;
  depart 67P 8 Dec 2033 (m = 2610 kg); Earth intercept 5 Nov 2036, v∞ = 5.620
  km/s, DEC = −23.0°, m = 2101 kg. Statistics over 100×2 h runs (Table 4,
  p. 41): analytic — best mass 2104 kg, 55 mean feasible solutions/run, 94%
  success, 16.6 min mean time-to-best; **finite-differenced — zero feasible
  solutions in all 100 runs** (Table 4 dashes; p. 20: even seeding FD with the
  analytic optimum, SNOPT declared it infeasible and diverged immediately).
- **Jacobian timing (Fig. 8, p. 37, read off the plot, approximate):** FD/analytic
  time ratio grows ≈160 (10 segments) → ≈1650 (200 segments); AD/analytic ≈25 →
  ≈730. CD relative accuracy typically no better than 1e-8 (p. 21, citing
  Pellegrini & Russell [69]).
- **Golden assessment (golden-rule discipline):** Table 3 + Fig. 6 is a fully
  sourced end-to-end input/output pair, but reproducing it requires the whole
  MGALT + SEP-hardware + Atlas-V + MBH stack and is stochastic (MBH) — usable
  only as a *stack-level* anchor if we ever build MGALT, not as a unit golden.
  **No leg-level numeric gradient values are published** (no worked STM/MTM
  example with numbers). For verifying an analytic-Jacobian implementation the
  paper's own recommended pattern (Sec. VI, p. 22) is cross-checking against
  AD/complex-step/central differences — an internal-consistency check, which is
  the appropriate tool for derivative code (distinct from catalogue golden
  values, which must stay externally sourced).

## 7. Adoption assessment for cyclerfinder

Current state: `search/dsm_leg.py` (Takao one-DSM legs, forward propagate +
Lambert back-fill), `search/correct.py` and the free-return chains all use scipy
`least_squares` with finite-difference Jacobians (the known bottleneck; we
parallelise FD columns). Core has CR3BP integrated STMs (`core/cr3bp.py`) but
**no analytic heliocentric two-body STM**.

- **Path A (minimal, keep our transcription):** analytic Jacobians for the
  existing Lambert-based DSM legs need *Lambert partials*, which this paper
  deliberately avoids and does not supply. Not served by this source (Lantoine &
  Russell [58] / Pellegrini & Russell [69] in its reference list are the leads).
- **Path B (paper-faithful re-transcription, recommended for the DSM lane):**
  recast each DSM leg/chain as an MGAnDSMs FBS phase — decision variables
  (v∞-out, v∞-in, Δv_DSM, α_k, Δt_p), mid-phase 6-element defect (Eq. 2, mass
  row dropped since our cycler legs don't track mass — MTM Eq. 29 then collapses
  to identity + the Eq. 42 velocity-slot connection), flyby continuity Eqs. 3–4
  with Appendix gradients. Ingredients we'd need to add to `core/`:
  1. analytic two-body STM (Battin Eq.-17 partition; Shepperd 1985 universal
     variables is the robust choice),
  2. Lagrange-coefficient time derivatives (Pitkin form, Eqs. 43–44),
  3. the chain assembler (Eqs. 31–32) — pure matrix products, trivially testable
     against our existing FD machinery.
  This removes Lambert from the DSM lane entirely and gives exact sparsity to
  the corrector. The MGA*n*DSMs piece is small: for massless legs the only
  nontrivial matrices are the STMs.
- **Immediately portable regardless of path:** the Appendix flyby-altitude
  gradients (Eqs. A1–A6) — algebraic, self-contained, useful anywhere we
  constrain flyby turn feasibility.
- **Strategic corroboration:** Table 4 (FD: 0/100 feasible; analytic: 94/100) is
  the strongest published evidence that our FD-Jacobian MBH lane is leaving
  robustness (not just speed) on the table — the same MBH+NLP architecture we
  run.

## 8. Reference leads worth pulling later

- [50] Ellison, Conway, Englander, Ozimek, "Application and Analysis of
  Bounded-Impulse Trajectory Models with Analytic Gradients," JGCD (companion).
- [52] Shepperd, "Universal Keplerian State Transition Matrix," Celestial
  Mechanics 35(2), 1985, 129–144 — the STM we'd implement.
- [57] Pitkin, "Second Transition Partial Derivatives via Universal Variables,"
  J. Astronautical Sciences 13(5), 1966 — the Eq. 43–44 time partials.
- [69] Pellegrini & Russell, "On the Computation and Accuracy of Trajectory
  State Transition Matrices," JGCD 39(11), 2016 — FD/CSD/analytic STM accuracy
  study (relevant to Path A and to validating Path B).
- [44] Vavrina, Englander, Ellison, "Global Optimization of N-Maneuver,
  High-Thrust Trajectories Using Direct Multiple Shooting," AAS 16-272 — the
  MGAnDSMs transcription origin.
