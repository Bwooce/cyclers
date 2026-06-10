# Mining Note: Cuevas del Valle, Urrutxua & Solano-López (2026), EuroGNC

**Source:** Sergio Cuevas del Valle, Hodei Urrutxua, Pablo Solano-López, "Fuel-optimal
Rendezvous in the CR3BP via MPC and Proximal Operators," CEAS EuroGNC 2026 conference,
Madrid, Spain, May 5–7 2026. Paper CEAS-GNC-2026-012, doi: 10.82124/CEAS-GNC-2026-012.
CC-BY 4.0. 25 pp. (21 pp. body + references).

**Mined:** 2026-06-10 from full PDF.

**Verdict: USEFUL (narrow).** This is a GNC/control paper (ADMM + primer vector theory +
MPC), not an orbit-catalogue paper — the solver itself is out of scope for a
search/validation codebase. But it yields: (a) one full-precision Earth-Moon L2 *southern*
halo initial condition (+ a derivable northern chaser IC) usable as a **corrector seed**,
(b) a sourced Richardson `c2` value for Earth-Moon L2, (c) a cheap, directly applicable
**Jacobi-gap minimum-ΔV lower bound** technique for the maintenance-ΔV lane, and (d) an
acquisition lead (their 2023 Floquet station-keeping paper, ref [47]).

---

## 1. Published CR3BP initial conditions (mining lens 1)

### 1.1 Earth-Moon L2 southern halo — target spacecraft (Sec. 7.3, p. 18)

The Scenario II target is on an L2 **southern** halo orbit (text: "a chaser and a target
spacecraft in an L2 northern and southern halo orbit, respectively"). Nondimensional
synodic state, printed at full precision (p. 18):

```text
[r_t, rdot_t](t0) = [0.824024728136525, 0, -0.054501847320725,
                     0, 0.164671964079122, 0]
```

This is a plane-crossing state (y = 0, vx = vz = 0, z < 0 ⇒ southern member) — exactly
the symmetric form our STM corrector in `src/cyclerfinder/search/cr3bp_periodic.py`
expects as a seed.

### 1.2 Derived chaser IC (northern halo)

The co-orbital (relative) state is also printed at full precision (p. 18):

```text
s(t0) = [-0.000385289210804, 0, 0.097783417040814,
         0, -0.012103983186502, 0]
```

Since s = chaser − target (Sec. 2, p. 4: ρ = r_c − r_t), the chaser (northern halo) IC is
target + s:

```text
[r_c, rdot_c](t0) = [0.823639438925721, 0, 0.043281569720089,
                     0, 0.152567980892620, 0]   (derived, not printed)
```

The chaser's orbital period is stated only as "approximately t_f − t_0 = 2.7549
nondimensional time units" (p. 18). Initial miss-distance: 37,594 km (p. 18).

### 1.3 Why these are SEEDS, not goldens

Per the golden-tests-sourced-only rule, these cannot serve as golden fixtures:

- **μ is never printed.** The Earth-Moon mass parameter is implicit throughout; the
  paper never tabulates it (Sec. 2 defines μ symbolically only, p. 4).
- **No Jacobi constant is tabulated** for either halo (the Jacobi constant is *used* in
  Sec. 7.2 to bound ΔV but its values are not printed).
- **Period is "approximately" 2.7549 nd**, for the chaser only, four significant figures.
- The chaser IC above is *derived* by us (target + s), not printed by the authors.

Status: usable as corrector seeds (run our `correct_periodic` with standard EM
μ ≈ 0.012150585; if it converges with a small correction, the converged orbit is *ours*,
sourced-seed provenance only).

### 1.4 Richardson c2 coefficient, Earth-Moon L2 (Sec. 7.2, p. 15)

```text
c2 = 3.190425213622208
```

This is the fundamental-frequency coefficient of the Relative Libration Linear Model
(RLLM, Eq. 4, p. 5), per Richardson 1980 ([48], Celestial Mechanics 22:231–236,
doi 10.1007/BF01229509). It depends only on the libration point and μ, so it is a sourced
check value if we ever implement linearized collinear-point dynamics — and it implicitly
pins the authors' Earth-Moon μ.

### 1.5 Scenario I relative states (Sec. 7.2, p. 15)

Dimensional SI co-orbital boundary states near Earth-Moon L2 (relative states only — not
absolute orbit ICs, so not catalogue-relevant; recorded for completeness):

```text
s(t0) = [6449.40, 65117.03, 22814.91, -0.0312, 0.0392,  0.2114]
s(tf) = [59066.09, 67728.64, 84015.47, -0.1087, 0.1616, -0.1730]
```

(positions m? — the paper says "dimensional S.I. units"; the plotted trajectories in
Fig. 4, p. 18, span ~100 km, consistent with metres.)

---

## 2. Relative-motion formulation (mining lens 2)

Three model tiers, all in the synodic frame (Sec. 2, pp. 3–5), citable for a future
maintenance-ΔV lane on CR3BP cyclers:

1. **Exact nonlinear co-orbital ODEs** (Eq. 2, p. 4): Hamiltonian 2nd-order ODEs for the
   relative state s = [ρ, ρ̇] between target and chaser, both under full CR3BP gravity.
   Attributed to Luquette 2006 [46] and the authors' own Aerospace 2023 paper [47].
2. **Rendezvous Linear Model (RLM)** (Eq. 3, p. 5): linearization about the target's
   position; Hessian Σ = −(κ1+κ2)I + 3κ1(e1⊗e1) + 3κ2(e2⊗e2) with
   κ_i = μ_i/‖r_t − R_i‖³ — i.e., the CR3BP analogue of Clohessy-Wiltshire, valid along
   any target trajectory (time-varying A(t), STM by numerical integration).
3. **Relative Libration Linear Model (RLLM)** (Eq. 4, p. 5): Σ constant near a collinear
   point, diag(1+2c2, 1−c2, −c2) — an LTI model decoupled from the exact trajectories.

Relevance to us: if we ever cost station-keeping for our 14 SILVER Saturnian Lyapunov
members, the RLM (tier 2) is the published, citable linearization to use along the
periodic orbit — and we already compute the STM it needs in
`src/cyclerfinder/core/cr3bp.py`.

**Directly useful technique (Sec. 7.2, pp. 15–16):** the Jacobi constant gives a
**minimum l2 ΔV lower bound** to bridge two CR3BP states: "The existence of an
energy-integral in the CR3BP —known as the Jacobi constant— allow us to establish the
minimum l2 ΔV-requirement to accomplish the transfer/bridge the energy gap between the
two initial conditions in an l2-sense." For their Scenario I gap this gives
ΔV_min = 0.2734 m/s (p. 16). This is a one-line computation on top of our existing
`jacobi_constant` and would give every transfer/maintenance estimate a rigorous,
sourced floor. **Recommended adoption.**

---

## 3. Solver structure — implement-for-us assessment (mining lens 3)

Architecture (Secs. 3–6): fuel proxy = L_{1,p} norm of control; Neustadt's moment-problem
relaxation (Problem IV SICP, Eq. 8, p. 7) solved by ADMM in closed form ("Neustadt
solver", Algorithm 1, p. 12); a primal ADMM solver for control-bounded problems
(Eq. 15, p. 13); both embedded in a shrinking-horizon nonlinear MPC loop (Sec. 6, p. 14).
Implementation: Matlab 2021b, ode113, modified OSQP (Sec. 7.1, p. 15).

**Verdict: do not implement.** This is a real-time embedded-guidance engine; our use case
is offline search/validation. Two citable theory facts worth keeping:

- **Impulse-count bound** (Theorem 1, p. 8, after Neustadt 1964): the fuel-optimal
  impulsive solution has at most n impulses, n = dimension of the fixed final condition
  (= 6 for full-state rendezvous). Useful sanity bound if we ever cost impulsive
  orbit-to-orbit transfers.
- The primer vector p(t) = Yᵀ(t)λ obeys the same linear dynamics/STM as the state
  (p. 8) — so our existing STM machinery is sufficient to evaluate primer-vector
  optimality of any impulse sequence at near-zero extra cost.

---

## 4. Numeric results tables (mining lens 4)

### Scenario I — L2 vicinity rendezvous, RLLM dynamics (Sec. 7.2)

| Quantity | Neustadt | Direct (primal) | Serra et al. [35] | Source |
|---|---|---|---|---|
| Sequence cost, p=1 [m/s] | 1.8714 | 1.8397 | 1.6384 | Table 1, p. 16 |
| Sequence cost, p=2 [m/s] | 1.3050 | 1.3050 | 1.2251 | Table 2, p. 17 |
| N* impulses (p=1 / p=2) | 4 / 2 | 4 / 3 | 4 / 2 | Tables 1–2 |
| Jacobi-gap ΔV_min [m/s] | 0.2734 | — | — | p. 16 |

**Provenance caveat (p. 16):** the authors state their results "cannot be reproduced
exactly in quantitative terms" against Serra et al. [35] because [35]'s unit system was
"not included in the original manuscript" — i.e., the Serra benchmark itself has a
known unit-provenance gap. Do not treat the Serra column as an independent golden.

### Scenario II — halo-to-halo transfer, MPC closed loop (Sec. 7.3)

| Quantity | Value | Source |
|---|---|---|
| Transfer time | 2.7549 nd (≈ 1 chaser period) | p. 18 |
| MPC sampling time T_s | 0.028 nd (≈ 2 h 53 min) | p. 18 |
| Total cost, MPC-Primal (RLLM), p=2 | 304.5375 m/s | p. 19; Table 3, p. 20 |
| Total cost, MPC-Primal (RLM), p=2 | 429.9452 m/s | Table 4, p. 20 |
| Chebyshev PS benchmark (full nonlinear) | 211.2 m/s | Tables 3–4, p. 20 |
| Final rendezvous error | 439.92 m, 2 mm/s | p. 19 |
| MPC suboptimality vs PS | up to ~30% (RLLM) | p. 20 |

The 211.2 m/s Chebyshev pseudospectral figure (Gong/Ross/Fahroo method, ref [71]) is the
closest thing to a "true" fuel-optimal halo-to-halo cost here, but it is for *this
specific* (untabulated-μ) IC pair — citable as context, not as a golden.

---

## 5. Actions / leads

1. **Seed run (cheap, optional):** feed Sec. 1.1's southern-halo state to
   `correct_periodic` with standard EM μ; if it converges, we gain an Earth-Moon L2
   southern halo member with sourced-seed provenance (V-level per our ladder — not a
   golden; no published Jacobi/period/μ to validate against).
2. **Adopt the Jacobi-gap ΔV lower bound** (Sec. 2 above) in the maintenance-ΔV /
   transfer-costing lane — one line on top of existing `jacobi_constant`.
3. **Acquisition lead:** ref [47] — S. Cuevas del Valle, H. Urrutxua, P. Solano-López,
   "Optimal Floquet stationkeeping under the relative dynamics of the three-body
   problem," Aerospace 10(5), 2023, doi 10.3390/aerospace10050393 (open access). That is
   the station-keeping-cost paper this one builds on, and a better match for our
   "maintenance-ΔV costing for CR3BP periodic orbits" question than this GNC paper.
   Secondary lead: ref [35], Serra, Arzelier, Bréhard & Joldes 2018, "Fuel-optimal
   impulsive fixed-time trajectories in the linearized circular restricted 3-body
   problem" (the SICP benchmark; note the unit-system caveat above).
4. **No catalogue writeback from this paper.** No cycler orbits, no complete (μ, state,
   period, Jacobi) tuples; nothing meets v4.2 backfill standards.
