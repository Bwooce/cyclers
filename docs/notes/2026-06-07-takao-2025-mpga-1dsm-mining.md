# Takao 2025 — MPGA-1DSM transcription mining (multi-arc / one-DSM-per-leg)

Mined 2026-06-07. **ALGORITHM-focused** mine: extracts the MPGA-1DSM decision
vector, the powered-flyby / Oberth bend-feasibility model, the PSO-SQP-MBH global
search, and the SPICE usage pattern — for the multi-arc topology frontier (the
front-runner explanation for the S1L1 / M-ED closure blocker). Context read
first: `docs/notes/2026-06-07-marginal-papers-triage.md` (entry #4),
`docs/notes/2026-06-07-external-algorithms-survey.md`,
`docs/notes/2026-06-07-jones-aas17-577-method-deepdive.md` (our B-plane /
powered-flyby kernel).

**Source (cite exactly, author/title/arXiv only — no file path):**
Y. Takao, "Mission Analysis for the First-Ever Saturn Trojan 2019 UO14,"
arXiv:2501.06586 (astro-ph.EP), 11 Jan 2025. Dept. of Aeronautics and
Astronautics, Kyushu University.

> 21 pp., clean digital typeset. The whole transcription is on **pp.5-8**
> (Eqs.1-17 + Figs.4-6); optimizer settings **p.10 (Table 4)**; results tables
> **pp.11-13 (Tables 5-13)**; low-thrust on pp.14-15; bound-generation recipe in
> the **Appendix (Eqs. A.1-A.3, p.19)**. Every equation transcribed verbatim.

---

## 0. HEADLINE: the transcription in 5 lines

1. **Genome** (Eq.1-2, p.5): per leg `Y_i = [V∞, α, β, τ, η]_i`; full vector
   `X = [t0, Y1, …, Y_N]`. `V∞,α,β` = outgoing hyperbolic excess speed +
   azimuth + elevation at body i; `τ` = leg time-of-flight; `η ∈ [0,1]` = DSM
   timing fraction along the leg. One DSM per leg, one decision-vector block per
   leg.
2. **Epoch chaining** (Eq.3): `t_j = t0 + Σ τ_i`; body states `r_p,j, v_p,j`
   pulled from **SPICE** at those epochs. Outgoing state at body i is
   `v = v_p,i + V∞,out,i` with `V∞,out,i = v∞·[cosα cosβ, sinα cosβ, sinβ]`
   (Eqs.4-5).
3. **DSM construction** (Eqs.6-7): propagate the outgoing state ballistically
   (2-body heliocentric) for `η·τ`; at the DSM point solve **Lambert** to the
   next body over the remaining `(1−η)·τ`; `ΔV_DSM = ‖v21 − v12‖` (velocity
   mismatch at the DSM point); incoming hyperbolic `V∞,in,i+1 = v22 − v_p,i+1`.
4. **Powered flyby** (Eqs.8-13): deflection feasibility from `δ = sin⁻¹(1/e)`,
   `e = 1 + rp·‖V∞‖²/μ`; if required bend ≤ `δ̄ = δ_in + δ_out` an **Oberth
   periapsis P-FB ΔV** (Eq.11) is charged; else a deflection P-FB ΔV (Eq.13).
   rp is *not* a free variable — it is solved implicitly by minimizing total ΔV.
5. **Objective + optimizer** (Eqs.15-17): minimise
   `f = ΣΔV_DSM + ΣΔV_P-FB + ΔV_arr` (Eq.15) plus a TOF penalty `g` (Eq.16),
   box-bounded, solved by a **parallelised PSO → SQP → MBH** hybrid with
   evolutionary-branching search-space pruning (Fig.6); body states from SPICE,
   patched-conic / zero-SOI throughout.

**MAPS-TO-OUR-GENOME VERDICT: a one-DSM-per-leg variant of our free-return /
multi-arc genome is a SMALL, well-specified extension.** Our free-return
corrector already carries (t0, per-leg V∞ direction, per-leg ToF); the only
genuinely new genome coordinate is **η (DSM timing fraction)** plus the
"propagate-η·τ-then-Lambert" leg evaluator. The flyby physics is the **same
family** as our `core/flyby.py` `max_bend` (their Eq.8-9 ≡ our `max_bend`), but
their **Oberth periapsis P-FB ΔV (Eq.11) is a strict improvement** over our
asymptote-rotation deficit Δv. Scope: **MED** (genome+evaluator LOW; full
PSO-SQP-MBH+pruning optimiser MED-HIGH; we already have an MBH wrapper).

**CATALOGUE-ELIGIBLE? NO** — one-shot probe to a single Saturn Trojan; no
re-encounter / repeating orbit. Tables 5-13 are **cross-check fixtures** for a
future MPGA-1DSM corrector variant, NOT catalogue rows. Confirmed §6.

---

## 1. THE GENOME, VERBATIM (pp.5-6)

### 1.1 Decision vector (Eqs.1-2, p.5)
> "For each trajectory leg between a pair of celestial bodies … the following
> decision variables are defined:
> `Y_i = [V∞, α, β, τ, η]_i`  (1)
> where V∞, α, and β represent, respectively, the magnitude, azimuth, and
> elevation angle of the outgoing V∞ vector at the ith body; τ ∈ [0,1] indicates
> the timing of DSM; … i = 1,…,N is the trajectory leg index. With the departure
> epoch denoted as t0, the complete decision vector for the MPGA-1DSM problem is
> defined as
> `X = [t0, Y1, …, Y_N]`  (2)"

> NOTE: the paper's prose has a typo — it says "τ ∈ [0,1] indicates the timing
> of DSM" but Eqs.4-7 and the Appendix make clear **τ is the leg ToF** and
> **η ∈ [0,1] is the DSM timing fraction**. Fig.4 confirms: `TOF τ` split into
> `(1−η)τ` (Lambert part) and `ητ` (propagate part). The Appendix (A.1) bounds
> `η ∈ [0,1]` explicitly. Use η for the fraction, τ for the ToF.

### 1.2 Epoch chaining + SPICE body states (Eq.3, p.6)
> "The epochs when the spacecraft reaches the jth body is given by
> `t_1 = t0 (j=1);  t_j = t0 + Σ_{i=1}^{j-1} τ_i  (j≥2)`  (3)
> … The position and the velocity of the jth body, r_p,j and v_p,j, can be
> computed using the SPICE Toolkit, corresponding ephemeris files, and the
> epochs t_j."

### 1.3 Outgoing-velocity parameterisation (Eqs.4-5, p.6)
> `v_{11,i} = v_p,i + V∞,out,i`  (4)
> `V∞,out,i = v∞,i · [cosα_i cosβ_i, sinα_i cosβ_i, sinβ_i]ᵀ`  (5)

### 1.4 DSM leg evaluator (Eqs.6-7, p.6)
> "Starting from the ith body with the velocity v_{11,i}, the spacecraft's
> trajectory is propagated for **ητ** under the assumption of a heliocentric
> two-body problem. This yields the position vector r_{12,i} and the velocity
> vector v_{12,i} at the point where the DSM occurs. Next, Lambert's problem is
> solved between the DSM position r_{12,i}=r_{21,i} and the (i+1)th body with a
> transfer time **(1−η)τ**, which gives the starting velocity v_{21,i} … and the
> arrival velocity v_{22,i} … The difference in velocity vectors at the DSM
> position corresponds to the DSM ΔV.
> `ΔV_DSM,i = ‖v_{21,i} − v_{12,i}‖`  (6)
> The incoming hyperbolic velocity at the (i+1)th body is
> `V∞,in,i+1 = v_{22,i} − v_p,i+1`  (7)"

**This is the clean recipe for inserting an interior impulse into a leg:**
propagate the front fraction ballistically, Lambert the remainder, charge the
mismatch. Exactly the "interior impulse on a return arc" primitive our multi-arc
return modelling needs.

---

## 2. THE POWERED-FLYBY / BEND-FEASIBILITY MODEL (Eqs.8-13, Fig.5, pp.6-7)

### 2.1 Deflection feasibility (Eqs.8-9)
> `δ_in,i+1 = sin⁻¹(1/e_in,i+1);  δ_out,i+1 = sin⁻¹(1/e_out,i+1)`  (8)
> `e_in,i+1 = 1 + (rp,i+1/μ_i+1)·‖V∞,in,i+1‖²;  e_out = 1 + (rp/μ)·‖V∞,out‖²`  (9)

### 2.2 Max-deflection test + Oberth P-FB (Eqs.10-11)
> "the maximum deflection angle `δ̄_i+1 = δ_in,i+1 + δ_out,i+1` can be determined
> … If the actual deflection angle
> `δ_i+1 = cos⁻¹( V∞,in·V∞,out / (‖V∞,in‖‖V∞,out‖) )`  (10)
> is smaller than its upper limit δ̄_i+1, the **Oberth effect can be utilized by
> conducting the ΔV maneuver at the periapsis passage. The P-FB ΔV in this case
> is given by
> `ΔV_P-FB = | √(‖V∞,in‖² + 2μ/rp) − √(‖V∞,out‖² + 2μ/rp) |`  (11)"

### 2.3 Deflection-deficit cases (Eqs.12-13, Fig.5 a/b/c)
If `δ` exceeds the bend cone, a deflecting ΔV is charged. When
`δ < 2·max(δ̄_in, δ̄_out)` (cases b1/b2) a simple speed-change suffices:
> `ΔV_P-FB,i+1 = | ‖V∞,out‖ − ‖V∞,in‖ |`  (12)
Otherwise (cases c1/c2) a vector-deflection ΔV (law of cosines):
> `ΔV_P-FB = √( ‖V∞,in‖² + ‖V∞,out‖² − 2‖V∞,in‖‖V∞,out‖·cos{δ − 2max(δ̄_in,δ̄_out)} )`  (13)

### 2.4 KEY: rp is implicit, not free
> "In the present study, a ΔV is permitted at each intermediate flyby to let the
> flyby be feasible; thus, the periapsis distance constraint is **implicitly
> incorporated by minimizing the total ΔV**. … The actual periapsis distance
> rp,i+1 in Eq.(11) can be computed by numerically solving the following
> equation: `δ_i+1 = δ_in,i+1 + δ_out,i+1`." (p.7)

So unlike Jones (rp solved from the required turn) or our `rp_min` clamp, Takao
lets the optimiser trade rp against ΔV; the periapsis-altitude floor is enforced
only as a box/penalty.

### 2.5 Arrival + objective (Eqs.14-17)
> `ΔV_arr = ‖V∞,in,N+1‖` (rendezvous) or `0` (flyby)  (14)
> `f(X) = Σ_{i=1}^N ΔV_DSM,i + Σ_{i=1}^{N-1} ΔV_P-FB,i+1 + ΔV_arr`  (15)
> `g(X) = 0 if Στ ≤ τ_max ; w(Στ − τ_max) otherwise`  (16), weight `w`.
> `min_X f(X)+g(X)  s.t.  X_lb ≤ X ≤ X_ub`  (17)
P-FB low-velocity-flyby penalty deliberately **dropped** in favour of a
"sufficiently high lower bound on V∞" (2 km/s for flybys, see §4) to keep the
hyperbolic assumption valid — a cleaner alternative to Jones's penalty.

---

## 3. GLOBAL SEARCH (pp.7-8, Fig.6, Table 4 p.10)

- **Optimiser: hybrid PSO → SQP → MBH** per subproblem. "Initially, Eq.(17) is
  solved using conventional PSO … The PSO-derived solution is then refined using
  a gradient-based local optimization method … SQP … further refined through
  **Monotonic Basin Hopping (MBH)**." (pp.8-9)
- **Outer loop = evolutionary-branching search-space pruning** (Vasile & De
  Pascale [15]): solve Eq.17 on user bounds → X0; contract domain by scale
  factor **λ=0.6** around the best → solve again → X1 → X2 …, fixed
  **5 iterations** per execution; **30 independent executions** per sequence,
  parallelised across threads (Fig.6). "different optima … obtained in parallel."
- **Settings (Table 4, p.10):** PSO individuals **1000**, generations **2000**;
  MBH perturbation **0.2**, max iterations **100**; contraction **λ=0.6**, TOF
  penalty weight **w=10**. Software: **PAGMO [22]** (Biscani & Izzo), "default
  values were employed."
- **Performance:** 5-flyby case, one execution (5 iters) ≈ 10 min on a laptop;
  30 executions in 10 min on a 32-thread cluster of 3 commodity machines; whole
  sequence-enumeration (≤200 sequences, 30×5 trials each) in 1-2 days "without
  the need for outer-loop optimization as in [14]."
- **Sequence enumeration:** integer body codes (Table 2); per-position code sets
  (Tables 3, 9); duplicates / >3-consecutive-same-planet / inner-after-outer
  sequences pruned → **152 unique flyby sequences**, **190 / rendezvous**.

**vs our scan/MBH stack:** Our MBH wrapper (`docs/notes/2026-06-07-mbh-wrapper.md`)
is the inner-loop analogue. What we **lack**: (a) the **PSO global pre-seed**
before SQP/MBH; (b) the **evolutionary-branching domain-contraction** outer loop
(λ=0.6, 5 iters) as an alternative to grid; (c) the **embarrassingly-parallel
30-execution** pattern that buys global-optimum probability for free; (d) the
**automatic bound generation** from the body sequence (Appendix, §5).

---

## 4. SPICE / FIDELITY PATTERN (pp.4, 6) — vs our DE440/rails split

- **Single fidelity throughout search:** body states `r_p,j, v_p,j` from the
  **SPICE Toolkit + ephemeris files** at chained epochs (Eq.3); **zero-SOI
  patched conic**, heliocentric 2-body propagation for the DSM front-arc; Lambert
  for the back-arc. Orbit-analysis §2 used "the **SPICE kernel provided by the
  JPL Horizons System**" (p.4) — DE version **not named**.
- **No separate validation ephemeris.** Unlike Jones (zero-SOI search → full
  n-body SNOPT) or our DE440-rails split, Takao does **not** re-optimise in
  higher fidelity; the patched-conic optimum *is* the result. The low-thrust §4
  adds direct collocation but still in 2-body heliocentric dynamics.
- **Comparison:** our search/validation split is *more* rigorous; Takao's value
  is the **single-fidelity ephemeris-patched-conic genome**, not the validation
  story. For a future MPGA-1DSM variant we would keep our DE440-rails validation
  and borrow only the genome/evaluator.

---

## 5. AUTOMATIC BOUND GENERATION (Appendix Eqs. A.1-A.3, p.19) — directly adoptable

> `α ∈ [−π, π];  β ∈ [−π/2, π/2];  η ∈ [0, 1]`  (A.1)
> Inner-pair leg ToF: `τ ∈ [30 days, P_s + P_H]`  (A.2)  (P_s synodic, P_H Hohmann)
> Outer / target leg ToF: `τ ∈ [0.3 P_H, 1.3 P_H]`  (A.3)
- Departure V∞ at Earth: **[1, 5.1] km/s** (5.1 = 1:2 Earth-resonant cap; 1 =
  prevent low-velocity flybys); 5.1 km/s ↔ C3=26 km²/s². Flyby V∞ lower bound
  **2 km/s** uniformly (prevents hyperbolic-collapse). V∞ upper bound
  **1.1× escape velocity** for inner planets, **9 km/s** for Jupiter/Saturn.
- "users need to specify only the epoch and the V∞ magnitude at launch. The
  bounds on other parameters are automatically determined based on the sequence
  of bodies, allowing the fully automated design of MPGA-1DSM trajectories."

**MAPS:** these are concrete, source-traced bound recipes for our scan — the
Hohmann/synodic-keyed ToF windows (A.2/A.3) and the 2-km/s flyby-V∞ floor are
directly adoptable, complementing Jones's "±50 d from Hohmann, V∞<5 km/s" seed
window.

---

## 6. GOLDEN-ELIGIBLE NUMBERS (cross-check fixtures, NOT catalogue rows)

Two fully-tabulated itineraries (verbatim). P-FBs/DSMs reported only if total ΔV
> 10 m/s. **CATALOGUE-INELIGIBLE** (one-shot, no re-encounter).

### Table 7 (p.11) — Flyby mission, EEMEJA, launch 2033-Nov-27
| Event | Epoch | V∞ [km/s] | Flyby alt [km] | ΔV [m/s] |
|---|---|---|---|---|
| Earth | 2033-Nov-27 | 5.10 | — | 0 |
| DSM-0 | 2034-Oct-25 | — | — | 61.87 |
| Earth | 2035-Nov-03 | 5.66 | 20293.61 | 0 |
| Mars  | 2038-Jun-10 | 9.56 | 1135.62 | 0 |
| Earth | 2039-Dec-18 | 11.31 | 638.17 | 0 |
| Jupiter | 2041-Oct-18 | 8.45 | 650639.49 | 29.36 |
| 2019 UO14 | 2044-Nov-17 | 10.02 | — | 0 |
Total ΔV (Table 5) = **91.82 m/s** (best flyby, 1st opportunity).

### Table 8 (p.12) — Flyby mission, EEEA, launch 2028-Dec-10
| Event | Epoch | V∞ [km/s] | Flyby alt [km] | ΔV [m/s] |
|---|---|---|---|---|
| Earth | 2028-Dec-10 | 5.10 | — | 0 |
| DSM-0 | 2029-Dec-15 | — | — | 222.78 |
| Earth | 2031-Jan-09 | 6.95 | 9588.04 | 0 |
| DSM-1 | 2032-Jun-08 | — | — | 542.97 |
| Earth | 2033-Nov-16 | 11.77 | 637.82 | 0 |
| 2019 UO14 | 2039-Aug-21 | 6.97 | — | 0 |
Total ΔV (Table 5 No.6) = **765.75 m/s**.

### Table 12 (p.13) — Rendezvous, EMEEJSA, launch 2033-Aug-03 (via Saturn)
| Event | Epoch | V∞ [km/s] | Flyby alt [km] | ΔV [m/s] |
|---|---|---|---|---|
| Earth | 2033-Aug-03 | 4.50 | — | 0 |
| Mars | 2034-Nov-02 | 5.18 | 1105.85 | 0 |
| Earth | 2036-Nov-19 | 11.89 | 13515.88 | 0 |
| Earth | 2038-Nov-20 | 11.88 | 749.71 | 0 |
| Jupiter | 2040-Aug-14 | 8.98 | 4.22×10⁶ | 0 |
| Saturn | 2043-Aug-14 | 5.69 | 254014.43 | 154.33 |
| 2019 UO14 | 2058-Mar-23 | 1.82 | — | 1822.46 |
Total ΔV (Table 10) = **1977.06 m/s**, TOF 24.6 yr (the SGA route).

### Table 13 (p.13) — Rendezvous, EMEEJA, launch 2035-Sep-18 (via Jupiter)
| Event | Epoch | V∞ [km/s] | Flyby alt [km] | ΔV [m/s] |
|---|---|---|---|---|
| Earth | 2035-Sep-18 | 4.40 | — | 0 |
| Mars | 2037-Feb-04 | 5.36 | 767.77 | 0 |
| Earth | 2039-Jan-10 | 10.22 | 13288.37 | 0 |
| Earth | 2041-Jan-09 | 10.23 | 1631.59 | 0 |
| Jupiter | 2042-Dec-08 | 8.05 | 1.73×10⁶ | 0 |
| 2019 UO14 | 2049-Feb-24 | 2.72 | — | 2716.44 |
Total ΔV (Table 11) = **2716.56 m/s**, TOF 13.4 yr (the JGA route).

Top-10 totals also available: Tables 5/6 (flyby 1st/2nd opp), 10/11 (rendezvous
SGA/JGA). Notable Mars-leg V∞ here: **5.18 / 5.36 km/s** (Tables 12/13) — same
5-7 km/s ballistic-Mars class as Hughes 2014 / Boeing 2033 (corroboration, not a
cycler).

**Bodies/constants (Table 2, p.9, source-traced GM/radius/min-alt):** Earth
GM 3.986004418e14 m³/s² R 6378 km alt-floor 637.8 km; Venus 3.24859e14 / 6052 /
605.2; Mars 4.2828e13 / 3397 / 339.7; Jupiter 1.26686534e17 / 71492 / 571936;
Saturn 3.7931187e16 / 60330 / 6033. (min-alt = 0.1× radius.) These are usable as
a sourced constants fixture; UO14 elements in Table 1 (p.3, ex Hui et al. 2024).

---

## 7. MAPS-TO-OUR-X verdicts

| Takao element | Our code / frontier | Verdict |
|---|---|---|
| Genome `[V∞,α,β,τ,η]` per leg + `t0` (Eq.1-2) | `search/free_return.py`, `correct.py` genome | **MAPS — small extension.** We already carry t0, per-leg V∞ dir + ToF. The ONE new coordinate is **η (DSM fraction)**. A "one-DSM-per-leg" free-return variant = add η + the Eq.6-7 evaluator. |
| DSM evaluator: propagate η·τ (2-body) → Lambert (1−η)·τ → ΔV=‖v21−v12‖ (Eq.6-7) | `core/lambert.py`, propagate; (no interior-impulse leg) | **PARTIAL → buildable.** We have Lambert + 2-body propagate; we lack the propagate-then-Lambert interior-impulse leg primitive. Clean recipe; LOW effort. This is the multi-arc-return primitive the S1L1 multi-arc finding needs. |
| Bend feasibility `δ=sin⁻¹(1/e)`, `e=1+rp·V∞²/μ` (Eq.8-9) | `core/flyby.py::max_bend` (`sin(δmax/2)=1/(1+rp·v∞²/μ)`) | **MAPS — same physics.** Their per-asymptote `sin⁻¹(1/e)` summed (δ_in+δ_out) equals our half-angle `max_bend` for the symmetric-rp case. Confirm the doubling convention when porting. |
| **Oberth periapsis P-FB ΔV** `|√(V∞in²+2μ/rp) − √(V∞out²+2μ/rp)|` (Eq.11) | `core/flyby.py::dv_from_turn_deficit` = `2·V∞·sin(deficit/2)` (asymptote rotation) | **DOES NOT MAP — and Takao is the IMPROVEMENT.** Our deficit Δv rotates the asymptote at infinity (over-estimates, no Oberth credit). Eq.11 charges the burn **at periapsis** (Oberth-aware, deep in the well) → strictly cheaper for the speed-change case. **Adopt Eq.11/12/13 as an Oberth-aware powered-flyby cost.** |
| rp implicit (solved from δ=δ_in+δ_out, traded vs ΔV) not clamped | our `rp_min` floor | **DESIGN DIFFERENCE.** Takao optimises rp; we clamp. Worth a flag — implicit-rp gives the optimiser a lever we currently deny it. |
| Objective `ΣΔV_DSM+ΣΔV_PFB+ΔV_arr` + TOF penalty (Eq.15-16) | `model/score.py`, `correct.py` residual | **MAPS — additive ΔV objective + box-bounded TOF penalty** is our shape; w=10, τ_max box adoptable. |
| PSO→SQP→MBH hybrid + evolutionary-branching pruning (λ=0.6, 5 iters, 30 exec) | MBH wrapper (Thread 1), grid scan | **PARTIAL.** We have MBH; we lack PSO pre-seed + domain-contraction outer loop + 30-way parallel-execution pattern. MED effort; PAGMO ([22]) is the off-the-shelf engine. |
| Automatic bound generation from body sequence (A.1-A.3) | `search/scan.py` bounds (hand-set) | **MAPS — adopt.** Synodic/Hohmann-keyed ToF windows + 2-km/s flyby-V∞ floor + 1.1×v_esc caps are source-traced, sequence-keyed bounds. |
| Single-fidelity SPICE patched conic (no separate validation) | DE440 search + rails validation | **OUR SPLIT IS STRONGER.** Borrow genome/evaluator only; keep our DE440-rails validation. |
| Low-thrust direct-collocation final arc (Eqs.18-28) | `core/sims_flanagan.py`, Pony-Express lane | **OUT-OF-SCOPE for the multi-arc mine** — ballistic frontier. Noted only: they convert the final ballistic arc to LT via DSM/P-FB-as-equality-constraints (Eq.19), normalized a0/Isp. |

---

## 8. SCOPE ESTIMATE — adopting the transcription

| Component | Effort | Note |
|---|---|---|
| η coordinate + propagate-then-Lambert leg evaluator (Eq.6-7) | **LOW** | One new genome coord + one leg primitive; Lambert + propagate already present. The minimal "one-DSM-per-leg" free-return variant. |
| Oberth-aware P-FB ΔV (Eq.11/12/13) replacing/augmenting `dv_from_turn_deficit` | **LOW-MED** | Drop-in alternative cost in `core/flyby.py`; needs the case logic (Fig.5 a/b1/b2/c1/c2) + implicit-rp solve `δ=δin+δout`. |
| Automatic sequence-keyed bounds (A.1-A.3) | **LOW** | Bound-generator helper for `search/scan.py`. |
| Full PSO→SQP→MBH + evolutionary-branching pruning, 30-way parallel | **MED-HIGH** | We have MBH; PSO + outer-loop pruning + parallel-execution harness is new. PAGMO available off-the-shelf but is a dependency decision. |
| Single-impulse-per-leg → full MPGA-1DSM multi-flyby integration end-to-end | **MED** | Composes the above; the genome already strings legs (Eq.2-3). |

**Overall: MED.** The genome + DSM evaluator + Oberth P-FB cost (the parts the
multi-arc frontier actually needs) are **LOW-MED**; only the full bespoke
optimiser is MED-HIGH and is largely redundant with our existing MBH.

---

## 9. v4.2 BACKFILL CHECKS + catalogue ruling

**CATALOGUE-ELIGIBLE? NO — confirmed.** One-shot flyby/rendezvous probe to a
single Saturn Trojan; no body re-encounter sustaining a repeating orbit. Tables
5-13 are **cross-check fixtures** for a future MPGA-1DSM corrector variant, never
catalogue rows. Provenance flags recorded per the v4.2 checklist in case any
itinerary is ever staged as a pipeline fixture:

- **center**: heliocentric (Sun-centered) 2-body for all transfer/DSM arcs;
  body-centered hyperbolic for flybys (rp, V∞ body-relative). If staged, legs
  `center: "Sun"`; flyby states body-centered.
- **tof_days_bounds**: per-event epochs ARE tabulated (Tables 7/8/12/13) → per-leg
  ToFs derivable and **SOURCED (table-listed)**, not derived. Search bounds keyed
  by Eqs. A.2 (`[30 d, P_s+P_H]` inner) / A.3 (`[0.3 P_H, 1.3 P_H]` outer-or-UO14);
  total-TOF caps: flyby τ_max **4000 d**, rendezvous **9000 d (SGA) / 5500 d (JGA)**.
- **source_ephemeris**: **SPICE Toolkit + JPL Horizons kernel, DE version
  UNSPECIFIED** (paper names neither DE440 nor a number). Flag **UNSPECIFIED**.
  Patched-conic / zero-SOI gravity model; `model_assumption: analytic-ephemeris
  patched-conic`. (Distinct from our DE440 rails — do NOT assume DE440.)

---

## 10. HONEST "not extractable" / caveats

- **τ vs η typo** (§1.1): the prose mislabels τ as the DSM timing; Eqs.4-7 +
  Fig.4 + Appendix A.1 make η the fraction and τ the ToF. Resolved above; flag
  when porting so the genome layout is not mis-read.
- **No worked numeric example of Eqs.8-13** (no tabulated δ/rp/ΔV_P-FB for a
  single flyby) → those equations are self-consistency checks when ported, not
  goldens, exactly as for Jones Eqs.1-5.
- **DE/ephemeris version never stated** → the itinerary epochs are reproducible
  only up to the unknown kernel; treat Table 7/8/12/13 dates as **good sourced
  anchors, second-tier to a kernel-pinned listing.**
- **PSO/MBH "default values"** beyond Table 4 not enumerated (PAGMO defaults).
- **Implicit-rp solve** (`δ=δin+δout`) is stated but the numerical method is not
  detailed ("numerically solving the following equation").
- Single most decisive finding → §0 headline + the **Eq.11 Oberth P-FB** as the
  one concrete physics improvement over our current `dv_from_turn_deficit`.
