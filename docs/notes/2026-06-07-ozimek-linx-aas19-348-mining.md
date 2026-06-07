# Ozimek / Riley / Arrieta 2019 — LInX low-thrust MGA optimization (algorithm mining)

Mined 2026-06-07 (Task #142). **NEVER MINED before.** Method mining for our
low-thrust v2 (Sims-Flanagan machinery in `src/cyclerfinder/core/sims_flanagan.py`)
and the Forge.

**Source (cite exactly, no file path):**
Ozimek, M. T., Riley, J. F., & Arrieta, J., "The Low-Thrust Interplanetary
Explorer: A Medium-Fidelity Algorithm for Multi-Gravity Assist Low-Thrust
Trajectory Optimization," (Preprint) AAS 19-348, AAS/AIAA Space Flight Mechanics
Meeting, 2019. JHU/APL + Nabla Zero Labs. 20 pages.

> Quality: PDF renders cleanly; equations and the Table 1 benchmark numbers
> read unambiguously. Vision read of all 20 pages.

---

## 1. The algorithm in 3 lines

Trajectory = a **directed graph of nodes (states) + edges (legs)**; each edge is
propagated by **Generalized Sims-Flanagan Transcription (GSFT)** — N segments,
one impulse at each segment midpoint, **forward+backward propagation to a match
point**, defect = match-point residual. The whole graph becomes one **sparse NLP**
(min f(x) s.t. node residuals + edge residuals + bounds) solved by **SNOPT** (SQP,
BFGS quasi-Newton Hessian) with **fully-analytic Jacobian via chained STM/MTM**.
Mapping verdict: **MAPS DIRECTLY and is the canonical reference architecture for our
v2** — our `sims_flanagan.py` is a faithful subset (single-leg, midpoint impulse,
fwd/bwd match point, 7-defect, `(8+3N)·M`/`7M` NLP sizing); LInX is the full
graph+SNOPT+analytic-partials superstructure we have not yet built.

---

## 2. The algorithm stated precisely

### 2.1 Trajectory representation (pp.2-3, Fig.2)
- A path is **nodes + edges in a graph**. Nodes = "locations in time and space",
  typed **initial** (launch/departure), **intermediate** (flyby), **final**
  (intercept/rendezvous) (p.2). Edges "adjoin two nodes and represent the
  spacecraft trajectory leg between the outbound node and the inbound node" (p.2).
- "Feasible paths require vectors of constraints, collectively known as
  **residuals**, to be satisfied on each node and edge" (p.2).
- Node state representation (Fig.2b): a node is `(t, v_inf, m)` — time, excess
  velocity, mass — "a common representation in **LInX**".

### 2.2 GSFT edge propagation (pp.3-5)
The four GSFT assumptions, verbatim (p.3):
1. "Each edge ... is subdivided into a prescribed number of *segments*, N."
2. "The effect of engine thrust on a particular segment is modeled by an impulse
   maneuver at the midpoint of the segment. The result is a
   coast-impulse-coast arc emulating a fixed-direction finite-burn arc."
3. "Coast arcs are represented by the designer's choice of trajectory
   propagation scheme."
4. "The spacecraft mass changes instantaneously at the midpoints of the segments
   in accordance with Tsiolkovsky's Rocket Equation."

State vector (Eq.1): `x = [r^T, v^T, m]^T` (7 variables — **matches our 7-defect**).

Time discretization (Eqs.2-3):
- `Δt = Δt_k = (t_m − t_o)/N` (p.3) — uniform in practice.
- `t_k = t_o + (k − 1/2)Δt` — impulse at segment midpoint (p.4).

**Forward+backward to match point** (p.3): "the spacecraft trajectory edge is
propagated forwards and backwards from the source and target node, respectively,
up to a state residual known as the match point. This ... places the match point
into a region of space and time that is distant from potential dynamical
bifurcation at the nodes, **reducing overall problem sensitivity**." — this is the
explicit *rationale* for our `match_index` (we default to temporal midpoint).

Coast propagation (Eqs.4-5): propagator `P` is the "designer's choice", with
"Common choices ... Keplerian universal variable propagation for speed or a
numerically integrated method" (p.4). **Our `core/kepler.propagate` = the
universal-variable choice.**

Impulse / throttle model (Eqs.6-9, p.4):
- `[Δv_k, Δm_k]^T = I(T_k, ṁ_k, m_k, Δt)` at each midpoint.
- `T_k = T_max,k · T_n,k` — normalized throttle vector `T_n,k`.
- `ṁ_k = ṁ_max,k · sqrt(T_n,k · T_n,k + λ²)` — `λ` is a "fictitious **mass leak**
  ... a small positive number ... crucial for gradient-based optimization"
  (prevents undefined mass-flow partials at zero throttle) (p.4).
- Throttle residual `δT_n,k = T_n,k · T_n,k` with `0 ≤ δT_n,k ≤ 1` (Eq.9) —
  enforces `|throttle| ≤ 1` as an NLP constraint.

Match-point residual (Eq.11, p.5): `δx_m = [(r+−r−)^T, (v+−v−)^T, (m+−m−)]^T`;
"For full state continuity, it is required that `δx_m = 0`." **This is our
`match_point_defect` 7-vector exactly.**

### 2.3 Analytic partials — STM/MTM chaining (Eqs.12-14, p.5)
Augmented state `X = [r^T, v^T, m, t_o, Δt, t_L]^T` (Eq.12). A **10×10 State
Transition Matrix Φ** (coast arcs) and **Maneuver Transition Matrix M** (impulse
maneuvers) (Eq.13) are **chained** to get match-point partials analytically:
`∂x_m+/∂v_inf,o = Φ_3 M_3 Φ_2 M_2 Φ_1 M_1 Φ_0 ∂X_0/∂v_inf,o` (Eq.14).
"a fully-analytic Jacobian in LInX results in a tremendous improvement in speed
and robustness of the NLP for LTTO problems" (p.8). **We do NOT yet have analytic
partials** — this is the single most valuable extension target (see §5).

### 2.4 Node residuals (Eqs.15-19, pp.5-6)
General node residual `δn = N(t_i, m_i, v_inf,i | t_o, m_o, v_inf,o)` (Eq.15) —
inbound-edge vs outbound-edge quantities.
- **Massive flyby (planet) node, patched-conic** (Eqs.16-17, p.5):
  - `||v_inf,i|| = ||v_inf,o|| = v_inf` (Eq.16) — **excess-speed-magnitude
    continuity. Identical to our N-arc corrector's magnitude-continuity residual
    (`search/correct.py`).**
  - turn angle: `||v_inf,o − v_inf,i|| = 2 v_inf sin(β_fb/2) = 2v_inf / (1 + (r_bod
    + h_fb) v_inf²/μ_bod)` (Eq.17) — **this is OUR `core/flyby.max_bend` /
    `sin(δ_max/2)=1/(1+rp·V_inf²/μ)` formula, identical.**
  - "used to impose constraints on flyby altitude (h_fb), turn angle (β_fb), and
    conservation of angular momentum" (p.5). For small bodies the turn angle is
    "constrained to be zero" (p.6).
- **Rendezvous node** (Eqs.18-19): match incoming/outgoing position+velocity;
  `m_o = m_i + δm` (user mass loss/gain), `t_o = t_i + δt` (δt=0 ⇒ instantaneous
  GA; δt≥0 ⇒ stay time).

### 2.5 The NLP (Eqs.30-32, p.8)
```
min_x  f(x)
s.t.   N_i^L ≤ N_i(x) ≤ N_i^U    i = 1..N_nodes     (node residuals)
       E_j^L ≤ E_j(x) ≤ E_j^U    j = 1..N_edges     (edge residuals)
       x^L ≤ x ≤ x^U
```
- "the majority of the NLP variables correspond to normalized thrust vector
  components within each GSFT segment" (p.8).
- "this structure leads to a **sparse NLP** of modest size that can be solved with
  conventional sparse solvers" (p.8). Fig.6: Jacobian sparsity — match-point
  partials (red) "comprise the majority of the Jacobian evaluation effort".
- Solver: **SNOPT** (Gill/Murray/Saunders) — "a sequential quadratic programming
  method that employs a BFGS quasi-Newton update to the Hessian"; LInX does NOT
  compute an analytic Hessian (p.8). `NLPSolver` interface also targets Ipopt /
  Knitro (p.9).

### 2.6 Mesh-refinement convergence loop (p.8) — IMPLEMENTABLE PROCEDURE
> "obtain a solution starting from a reasonably coarse impulse grid, progressively
> increase the number of impulses and re-solve the previous optimal solution, and
> iterate until the objective value no longer shows improvement." (p.8)
Caveat (footnote p.8): "there can be pathological cases where adding impulses can
lead to ... worse performance." Trade-off: few impulses ⇒ fewer vars but higher
sensitivity (each propagation spans longer); more impulses ⇒ larger but more
robust NLP, equal-or-better objective (p.8).

### 2.7 Hardware models (Eqs.20-28, pp.6-7) — for the Forge's fidelity layer
- `m_L,max = f(C_3)`, `T_max = f_dc·f(P_av)`, `ṁ_max = f_dc·f(P_av)` (Eqs.20-22);
  `C_3 = v_inf,o · v_inf,o` (Eq.23).
- Polynomial performance models `P_n(x)` (Eq.24): "Third-degree polynomials (n=3)
  ... for 5 kg of launch mass resolution and for engine throttling, n=6 is
  typically representative" — coefficients by **least-squares fit on performance
  data** (p.6).
- Solar array power model (Eqs.25-28): `P_av` with bus reservation `P_sc` and
  margin `E_mg`; radial degradation `P_rel = α_p κ(r_s)/r_s²` (α_p = 1 AU², Eq.28);
  time degradation `ε = β_1 + β_2 e^{β_3 t_L} + β_4 t_L` (Eq.27).

---

## 3. Maps to our X / does not map

| LInX construct | Our code | Verdict |
|---|---|---|
| GSFT edge (N seg, midpoint impulse, fwd/bwd match point) | `core/sims_flanagan.py` `SimsFlanaganLeg`, `propagate_forward/backward` | **MAPS — same transcription.** Their match-point rationale ("distant from dynamical bifurcation") justifies our `match_index` midpoint default. |
| `x = [r,v,m]^T`, `δx_m = [Δr,Δv,Δm]` = 0 | `match_point_defect` 7-vector | **MAPS exactly** (their Eq.11 = our defect). |
| NLP `(node residuals, edge residuals)` sparse, SNOPT | `nlp_dimensions` (`(8+3N)·M` vars, `7M` constraints), `chain_defect` | **MAPS — our sizing/assembly matches; we have no SNOPT driver yet** (we use SLSQP/least_squares in `search/`). |
| Patched-conic flyby node: `‖v_inf,i‖=‖v_inf,o‖`, turn `2v_inf/(1+(r+h)v_inf²/μ)` (Eqs.16-17) | `core/flyby.max_bend`, `flyby_bend_slacks`, `search/correct.py` magnitude continuity | **MAPS exactly** — their Eq.17 is our bend formula; their Eq.16 is our magnitude-continuity residual. |
| Mass-leak `λ` in `ṁ` (Eq.8) | (absent) | **DOES NOT MAP — gap.** We have no throttle-vector mass-flow model (our leg takes Δv schedule directly, no `T_n,k`/`λ`). Needed only if we add a thruster/power model. |
| Analytic STM/MTM Jacobian (Eqs.12-14) | (absent — we rely on numerical least_squares Jacobian) | **DOES NOT MAP — highest-value gap (see §5).** |
| Hardware models (launch C_3 poly, SEP power) | (absent) | **DOES NOT MAP** — Forge currently ballistic; this is the SEP-fidelity recipe if/when v2 grows a propulsion model. |
| Graph-of-nodes/edges path object | our chains are ordered tuples/sequences | **PARTIAL** — we model chains linearly (a list), LInX uses a full directed graph holding alternative paths (Fig.2a). Their graph buys "holding alternatives"; not needed for our single-family corrector. |
| Mesh refinement (coarse→fine impulse grid) | (absent) | **DOES NOT MAP — cheap, implementable win** for our SF leg solves. |

---

## 4. Candidate test anchor (tabulated → golden-eligible)

**Table 1 (p.17) — Dawn-like Earth→Mars→Vesta→Ceres benchmark.** This is the
ONLY tabulated numeric solution, given for three independent tools
(LInX / MALTO / HILTOP). Because the EXPECTED values trace to a published source
(and to two *independent* tools — MALTO direct, HILTOP indirect — cross-checking
each other), this is **golden-discipline-clean** as a low-thrust regression anchor
*if we ever implement the full SF+NLP+hardware stack*. Verbatim numbers:

Setup (Table 1a): 1× NEXT thruster; P_0 = 10.3 kW, P_sc = 0.25 kW; launch
vehicle Delta II 2925-9.5 (no contingency); launch declination limits ±28.5°;
thruster startup delay 45 days; launch date 9/27/2007; Mars flyby date
2/17/2009; min flyby altitude 300 km; arrive Vesta Aug 2011, stay 270 d; arrive
Ceres Feb 2015; planetary ephemerides Horizons SPK; mass scale 1000 kg;
λ = 1.0×10⁻⁴; NLP solver SNOPT; major feasibility tol 1.0×10⁻⁸, optimality tol
1.0×10⁻⁶; major iterations required 915; max C_3 allowed 5.1529 km²/s²;
launch mass 1114.4 kg; segments per leg (34, 60, 70).

Solution comparison (Table 1b), **LInX / MALTO / HILTOP**:
- Leg 1 Earth-Mars: launch date 9/27/2007 (all); launch C_3 5.1529 / 5.1529 /
  5.2285 km²/s²; launch declination 28.5 (all) deg; launch mass 1114.4 / 1114.4 /
  1105.2 kg; flight time 510 / 510 / 510 d; arrival mass 1039.0 / 1039.8 /
  1032.7 kg; propellant 75.4 / 74.6 / 72.4 kg.
- Leg 2 Mars-Vesta: flyby date 2/18/2009 (all); flyby v_inf 4.11 / 4.10 / 4.11
  km/s; passage altitude 300 (all) km; flight time 894 / 894 / 827.7 d; arrival
  date 8/1/2011 / 8/1/2011 / 5/26/2011; arrival mass 906.6 / 907.3 / 901.4 kg;
  propellant 132.4 / 132.5 / 131.3 kg; stay 270 (all) d.
- Leg 3 Vesta-Ceres: departure date 4/27/2012 (all); flight time 1,037 / 1,038 /
  1,038 d; arrival 2/28/2015 (all); arrival mass 806.7 / 807.2 / 802.3 kg;
  propellant 99.9 / 100.1 / 99.1 kg.
- **Total propellant 307.7 / 307.2 / 302.8 kg; mission duration 2,711 d (all).**

Headline agreement claim (p.17): "the objective function of 806.7 kg found in
**LInX** is within 0.5 kg or **0.06%** of the MALTO solution."

> Golden caveat: this is an Earth-Mars-**Vesta-Ceres** asteroid mission, NOT an
> Earth-Mars cycler. It is a low-thrust SF-transcription regression anchor for the
> v2 machinery, not a cycler-catalogue row. Use only if/when we build the
> SF+SNOPT+hardware stack; it cannot anchor anything in the current ballistic
> catalogue.

No other tabulated trajectory exists in the paper (Fig.10 is a plot only).

---

## 5. Single most implementable finding (this paper)

**Analytic match-point partials via chained STM/MTM (Eqs.12-14), plus the
coarse→fine mesh-refinement loop (§2.6).** Our `sims_flanagan.py` already has the
exact transcription LInX uses; the missing piece that LInX calls out as "a
tremendous improvement in speed and robustness" is the **analytic Jacobian** of
the 7-defect w.r.t. the decision variables (Δv schedule, t_o, Δt). The detailed
derivation is in their ref [18] (Ellison/Englander/Ozimek/Conway, "Analytical
Partial Derivative Calculation of the Sims-Flanagan Transcription Match point
constraints," AAS 14-310). Implementing even a forward-mode STM/MTM chain for our
single-leg `match_point_defect` would replace the numerical least_squares Jacobian
and is a self-contained, testable addition (regression: analytic vs finite-diff
partials agree to tol). The mesh-refinement wrapper (re-solve on a finer impulse
grid, accept while objective improves) is an even cheaper near-term win.

---

## 6. v4.2 backfill checks

This paper carries **no catalogue rows** (it is a method/software paper; its one
numeric example is a non-cycler asteroid mission). So no center / tof_days_bounds
/ source_ephemeris backfill applies to catalogue.yaml from this source.

For completeness, the source_ephemeris LInX *uses* is NAIF SPICE / Horizons SPK
with DE430/DE431 (refs [12][29], Table 1a "planetary ephemerides: Horizons SPK") —
relevant only as provenance if the Dawn benchmark is ever added as a v2 anchor
(it would be `source_ephemeris: DE430/DE431 via Horizons SPK`).

---

## 7. Honest "not extractable" list

- No cycler trajectories, no orbital elements (a,e,i,…) for any solution.
- Throttle/thrust time histories are plotted (Fig.10b) but not tabulated.
- The STM/MTM matrix entries (Eq.13) are structural (block layout) only; the
  scalar partial expressions are deferred to ref [18] (not in this PDF).
- The hardware-model polynomial coefficients (β's, α_p aside) are not given
  numerically — they are fit per-mission.
