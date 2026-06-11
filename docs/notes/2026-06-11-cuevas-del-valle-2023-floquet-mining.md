# Mining Note: Cuevas del Valle, Urrutxua & Solano-López (2023), Aerospace

**Source:** Sergio Cuevas del Valle, Hodei Urrutxua, Pablo Solano-López, "Optimal
Floquet Stationkeeping under the Relative Dynamics of the Three-Body Problem,"
*Aerospace* 2023, 10(5), 393, doi 10.3390/aerospace10050393. Open access (CC-BY). 26 pp.
Received 30 Mar 2023, published 25 Apr 2023.

**Mined:** 2026-06-11 from the full PDF. This is ref [47] of the same group's EuroGNC
2026 MPC paper — i.e., acquisition-lead item 3 of
`docs/notes/2026-06-10-cuevas-del-valle-2026-cr3bp-mpc-mining.md` is now **CLOSED**.

**Verdict: USEFUL.** Yields: (a) a complete, implementable Floquet station-keeping
formulation whose key trick — augmenting the unstable Floquet mode with a
**Jacobi-constant error state** — independently confirms the Jacobi constant as the
right energy coordinate for the task #190 maintenance-ΔV lane; (b) one published
station-keeping ΔV table (Table 1, p. 19: 18.5–22.7 m/s to absorb a ~330 km insertion
error on an Earth-Moon L2 halo within one period, low-thrust) plus a Monte-Carlo cost
range — our first sourced CR3BP maintenance-cost anchor; (c) two Earth-Moon halo ICs
(5 s.f., **seeds only** — μ and Jacobi never printed); (d) a cross-paper finding: the
2026 paper's "L2 southern halo" IC is almost certainly an **L1** halo (Sec. 3.3 below).

---

## 1. Floquet station-keeping formulation (mining lens a)

All in the nondimensional CR3BP synodic frame (Sec. 2.1, pp. 4–5; EoM Eq. (2), p. 5).
The station-keeping problem is recast as regulation of the *relative* state s → 0
against a virtual target on the reference orbit (Sec. 3, p. 15).

**Monodromy eigenstructure use (p. 15, Sec. 3):**

- Floquet decomposition of the monodromy matrix: `Φ(T,0) = F(0) e^{JT} F^{-1}(0)`,
  with F the (T-periodic) Floquet modal matrix and J the constant Jordan form of
  Floquet exponents.
- Floquet basis vectors from the STM eigenpairs (e_i, λ_i):
  `F_i(t) = e^{-ω_i t} e_i`, `ω_i = (1/T) log(λ_i)` (p. 15). The Floquet
  representation is preferred over raw monodromy "due [to] its stable numerical
  behavior" (p. 15).
- STM transport along the orbit: `Φ(t,t0) = F(t) e^{J(t-t0)} F^{-1}(t0)` (p. 15).
- Coordinate change `Fα = s` turns the linearized relative dynamics (RLM, Sec. 2.4)
  into `α̇ = Jα + F†Bu` (p. 15); basis evolution `Ḟ = (∂f/∂s)F − FJ` (Eq. (16), p. 16).

**Control formulation (pp. 16–17):**

- Naive approach — SDRE regulation of the unstable Floquet component α₁ alone
  (`α̇₁ = (J₁₁ − F₁(t)ᵀΠ(t))α₁`) — **fails**: with the center manifold unactuated, the
  spacecraft drifts along the orbit family's continuation direction even though α₁ → 0
  (Fig. 4, p. 16). A documented failure mode worth remembering for any "kill the
  unstable mode" maintenance scheme.
- **Fix (the paper's core idea, p. 17):** LPOs form one-parameter families with the
  Jacobi constant C (or period T) as continuation parameter; tie the controller to the
  reference orbit by augmenting the state with the energy error `e_C = C − C_ref`
  (only C is used, since T(state) has no analytic expression). The SDRE plant becomes

  ```text
  d/dt [α₁; e_C] = [[J₁₁, 0]; [C₁₁, 0]] [α₁; e_C] + [F₁; ∇C·B] u
  ```

  with C = ∇C·F·J, ∇C the Jacobi-constant gradient w.r.t. the absolute chaser state
  (p. 17). Hybrid SDRE + energy-shaping/Koopman control; ensures the controlled orbit
  stays energetically tied to C_ref — "a necessary condition for transport phenomena in
  the CR3BP" (p. 17). A 2-state, 3-input reduced-order controller.
- Practical settings (pp. 18–19): penalties Q = I₂, R = 5·10⁻³ I₃; maneuvers executed
  only over 30% of the orbital period, at perilune, after insertion; max acceleration
  constrained ≤ 0.5 mm/s². The LTI/LQR variant (linearized around the monodromy
  eigenbasis at insertion) is valid only over ~10% of T (p. 19); LQR comparison weights
  Q_LQR = 2·I₉, R_LQR = 0.02·I₃ (p. 20).

**Implementability for us:** we already compute the STM
(`src/cyclerfinder/core/cr3bp.py`); the monodromy eigendecomposition gives F(0), and
∇C is one differentiation of our existing `jacobi_constant`. The only genuinely new
machinery is the Floquet transport (Eq. (16)) — cheap. For task #190 this paper is the
maintenance-lane half: the 2026 EuroGNC paper uses the Jacobi gap as a transfer-ΔV
*lower bound*; this one uses the Jacobi error as the *controlled state* for maintenance.
Same coordinate, two complementary sourced uses.

## 2. Published ΔV budgets (mining lens b) — the maintenance-cost anchor

**Scenario (Sec. 3 "Stationkeeping Examples", pp. 17–18):** Earth-Moon L2 *southern*
halo, out-of-plane amplitude A_z = 30,000 km, orbit "as given in [26]" (Nazari et al.
2014, AIAA 2014-4140). Insertion-error correction over **one orbital period** (mission
duration 3.76 nd units = one reference-halo period, p. 18). Insertion dispersions
1σ = 190 km position / 4 m/s velocity (p. 17, "similar values to those used in [26]");
the deterministic case uses the specific relative IC of Sec. 3.2 below (initial range
332.55 km).

**Table 1, p. 19 — performance of the FSK schemes:**

| Controller | TOF | ISE | IAE | ΔV_T [m/s] | ‖u‖_min [mm/s²] | ‖u‖_max [mm/s²] | Comp. time [s] |
|---|---|---|---|---|---|---|---|
| SDRE (RFSK) | 0.3 T | 0.00025 | 0.02878 | 22.73 | 0.03 | 0.5 | 28.93 |
| LQR | 0.1 T | 0.00011 | 0.01817 | 18.47 | 0.04 | 0.5 | 1.10 |
| None | T | 0.68393 | 1.02245 | 0.00 | 0.00 | 0.00 | 0.07 |

ΔV_T = ∫‖u‖dt over the flight (Appendix B, p. 23). Uncontrolled, the spacecraft
"naturally diverged to the Earth-Moon realm" within the period (p. 18) — i.e., for this
halo the divergence timescale under ~200 km-class insertion errors is **one period**.

**Monte Carlo (p. 19, Figs. 8–10, pp. 20–21):** 1000 draws, same orbit/parameters,
dispersions 1σ = 100 km / 0.5 m/s. Costs are scatter-plotted only (no table): ΔV ranges
roughly 0.5–6 m/s for both RFSK and LQR, growing with initial position error
(Fig. 8, p. 20); RFSK clearly better in final position/velocity error at equal-or-lower
control cost (Fig. 9, p. 20; conclusion, p. 21).

**Caveats before using as a maintenance anchor:**

1. This is **insertion-dispersion correction over one period**, not steady-state annual
   station-keeping — it cannot be quoted as an m/s-per-year rate.
2. The deterministic 18.5–22.7 m/s row corresponds to a large (332.55 km) initial
   offset; the 100 km-class Monte-Carlo errors cost an order of magnitude less.
3. Low-thrust (≤ 0.5 mm/s²), maneuvering window restricted to 30% of T — costs are
   policy-dependent, not a dynamical floor (no Jacobi-gap bound is computed here; that
   technique is in the 2026 paper).

Usable anchor statement: *correcting a 100–300 km-class insertion error on an
Earth-Moon L2 halo costs O(1–25) m/s within one orbit, low-thrust (Table 1 +
Figs. 8–10).* Order-of-magnitude, scenario-tied.

## 3. CR3BP periodic-orbit ICs (mining lens c)

### 3.1 Earth-Moon L1 standard halo, A_z = 20,000 km (Sec. 2.5 example, p. 12)

Target absolute state (nd synodic, plane-crossing form) and relative state:

```text
S_t(t0) = [0.82413, 0, 0.05680, 0, 0.16725, 0]
s(t0)   = [0.00058, 0, 0.1,     0, 0.01139, 0]
```

Used for the model-accuracy comparison of Fig. 2 (p. 12); integrated "for an orbital
period of the chaser periodic orbit" — period value not printed.

### 3.2 Earth-Moon L2 southern halo, A_z = 30,000 km (Sec. 3 examples, p. 17)

"As given in [26]" (Nazari, Anthony & Butcher 2014):

```text
S_t(t0)  = [1.08238, 0, 0.06460, 0, 0.28198, 0]
s_c(t0)  = [0.00050, -0.00050, 0.00050, 0.01268, 0.00049, 0.00048]   (initial range 332.55 km)
```

Period: stated only indirectly — mission duration "3.76 non-dimensional units, which
corresponded to an orbital period of the reference halo orbit" (p. 18). Three s.f.,
and unusually long for a small-amplitude EM L2 halo (typical ≈ 3.4 nd) — treat with
suspicion until checked against the upstream source [26].

### Why these are SEEDS, not goldens (per golden-tests-sourced-only)

- **μ is never printed** (defined symbolically only, p. 5).
- **No Jacobi constant value is printed anywhere** — C appears as a controlled quantity
  but its numeric values are never tabulated.
- States are 5 s.f.; the only period is 3 s.f., indirect, and physically suspect.
- The 3.2 orbit has an upstream source ([26]) — precision should be pulled from there,
  not from this transcription.

Status: both usable as `correct_periodic` seeds with standard EM μ; converged orbits
carry sourced-seed provenance only.

### 3.3 Cross-paper IC finding: the 2026 paper's "L2" label is likely wrong

The 2026 EuroGNC note records its Scenario II target as an "L2 southern halo" at
x₀ = 0.824024728136525. But x ≈ 0.824 lies Earth-side of the Moon near L1
(x_L1 ≈ 0.8369, Moon at 1−μ ≈ 0.98785) — an L2 halo has no y = 0 crossing there. This
2023 paper explicitly labels its nearly identical state (x = 0.82413, Sec. 3.1 above,
same family neighborhood, opposite z sign) an **L1 standard halo** (p. 12). Conclusion:
the 2026 paper's "L2" scenario label is most plausibly a misnomer; the seed run proposed
in the 2026 note should treat that state as an Earth-Moon **L1 southern** halo and
expect the corrector to converge in the L1 family.

## 4. Relative-dynamics model stack (canonical-source upgrade)

The model tiers summarized from the 2026 paper in the prior note are *derived* here —
this 2023 paper is the canonical citation for the stack:

- Exact nonlinear co-orbital ODEs: Eq. (3), p. 6 (= the 2026 paper's Eq. 2); compact
  control-affine form ṡ = f(μ,s,r_t) + u (p. 6); variational/STM Eq. (4), p. 6; the STM
  is symplectic (det = 1, reciprocal eigenvalue pairs, p. 6).
- Hamiltonian formulation: relative Hamiltonian Eq. (7), p. 7; Legendre-series
  potential with relative coefficients c̄_{n,i} = μ_i/‖R_i − r_t‖^{n+1} (Eqs. (8)–(10),
  p. 8). Richardson's classical libration-point model is recovered as the special case
  of a target fixed at the libration point (p. 8).
- **The relative problem has NO Jacobi constant** (time-dependent Hamiltonian, p. 8) —
  only an energy-like integral C = −2Ũ − V² − ∫(∂Ũ/∂t)dt (pp. 8–9), and zero-velocity
  surfaces are "much more complex or even unavailable" (p. 9). Relevant caution for any
  relative-dynamics costing work.
- RLM (Rendezvous Linear Model): p. 9, Σ = −(κ₁+κ₂)I + 3κ₁(e₁⊗e₁) + 3κ₂(e₂⊗e₂),
  κ_i = μ_i/‖r_t − R_i‖³ — time-varying along the target orbit.
- RLLM: Eq. (12), p. 10, Σ = diag(1+2c₂, 1−c₂, −c₂);
  c₂ = Σ_i μ_i/‖R_i − R₂ ± γ‖³ recovered analytically (p. 13). (The 2026 paper's
  numeric c₂ = 3.190425213622208 for EM L2 remains the only printed value across the
  two papers.)
- Higher-order co-orbital models: second-order Eq. (13) and third-order Eq. (14),
  p. 11; N-th order Eq. (11), p. 10. Accuracy comparison Fig. 2, p. 12: linear best
  short-term, higher orders better long-term.
- Normal form of the second-order relative Hamiltonian (p. 13):
  H₂ = λq₁p₁ + (ω_p/2)(q₂²+p₂²) + (ω_v/2)(q₃²+p₃²) — relative motion near the collinear
  points analytically inherits the saddle × center × center structure of the absolute
  problem; (λ, ω_p, ω_v) are the eigenvalues of the relative STM at the rendezvous
  condition. First-order neutrally-stable relative motion is a Lissajous curve
  (Eq. (15), p. 14).
- **Encke formulation of the relative CR3BP** (Appendix A, p. 22): full equations for
  numerically well-behaved integration when ρ ≪ ‖r_t − R_i‖, including the analytic
  Jacobian for the STM. Directly reusable if we ever integrate relative states.
  Integrator used: variable-step variable-order Adams-Bashforth-Moulton, orders 1–13,
  Matlab 2021b (pp. 22–23); hardware i7-1165G7 @ 2.80 GHz (p. 23).

## 5. Overlap / diff vs the 2026 EuroGNC paper (mining lens d)

- **Lineage:** this paper = ref [47] of the 2026 paper; the prior note's acquisition
  lead 3 is closed by this mining.
- **Overlap:** the relative-dynamics model stack (exact ODEs, RLM, RLLM, c₂) — the 2026
  paper cites it; the derivations live here. Cite the 2023 paper for the models.
- **Complementary, not duplicate:** 2023 = continuous low-thrust *station-keeping*
  (Floquet + SDRE + Jacobi energy-shaping) with a published SK ΔV table; 2026 =
  impulsive fuel-optimal *rendezvous* (Neustadt/ADMM/MPC) with rendezvous ΔV tables and
  the Jacobi-gap ΔV lower bound. For task #190: maintenance-lane anchor ← 2023 (this
  note, Sec. 2); ΔV lower bound ← 2026 (prior note, Sec. 2).
- **No numeric overlap:** different orbits and scenarios; no value here cross-checks a
  value there. The only shared numeric thread is the x ≈ 0.824 halo IC pair, which
  exposes the 2026 L1/L2 labeling problem (Sec. 3.3).

## 6. Actions / leads

1. **Task #190 maintenance lane:** adopt Table 1 (p. 19) + the Monte-Carlo range
   (Figs. 8–10) as the sourced order-of-magnitude SK-cost anchor for Earth-Moon
   halo-class orbits, with the Sec. 2 caveats attached; pair with the 2026 paper's
   Jacobi-gap lower bound for a floor + empirical-cost bracket.
2. **Amend the prior note's seed-run action:** treat the 2026 Scenario II target state
   as an Earth-Moon **L1** (not L2) southern halo (Sec. 3.3); same for the derived
   northern chaser IC.
3. **Optional seed runs:** Secs. 3.1–3.2 states through `correct_periodic` with
   standard EM μ — sourced-seed provenance only. For 3.2, the printed 3.76 nd period is
   suspect; let the corrector determine the period and compare.
4. **Acquisition leads:**
   - Ref [26]: Nazari, Anthony & Butcher 2014, "Continuous Thrust Stationkeeping in
     Earth-Moon L1 Halo Orbits Based on LQR control and Floquet Theory," AIAA 2014-4140
     — the upstream source of the Sec. 3.2 halo and the dispersion model; note its
     title says **L1** while this paper's usage says L2 — pulling it would resolve both
     the labeling and the precision questions. **Priority lead.**
   - Ref [74]: Folta, Pavlak, Haapala, Howell & Woodard 2014, "Earth-Moon libration
     point orbit stationkeeping: Theory, modeling, and operations," Acta Astronautica
     94, 421–433 — operational (ARTEMIS) EM LPO SK costs; the best candidate for a
     real-mission maintenance-ΔV anchor.
   - Ref [30]: Shirobokov, Trofimov & Ovchinnikov 2017, "Survey of Station-Keeping
     Techniques for Libration Point Orbits," JGCD 40, 1085–1105 — survey with cost
     comparisons across methods; candidate compendium source for the maintenance lane.
5. **No catalogue writeback.** No complete (μ, state, period, Jacobi) tuple is printed;
   nothing meets v4.2 backfill standards. Negatives: no μ value, no Jacobi values, no
   stability indices/eigenvalue tables, no orbit-family tables anywhere in the paper.
