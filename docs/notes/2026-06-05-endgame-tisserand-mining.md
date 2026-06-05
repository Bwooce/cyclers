# Data Mining: Campagnola & Russell "Endgame Problem" Parts 1 & 2 (2010)

Extracted 2026-06-05. Method: Read tool with `pages=` on the AAS 2009 preprint PDFs
(text-layer PDFs; both 20 pp). All quotes below are verbatim transcriptions from the page
images, with page numbers as printed on each page.

---

## Full Citations

**[CR-PART1]**
S. Campagnola and R. P. Russell, "The Endgame Problem Part 1: V∞-Leveraging Technique
and the Leveraging Graph," *Journal of Guidance, Control, and Dynamics*, Vol. 33, No. 2,
2010, pp. 463-475, DOI 10.2514/1.44258. (Preprint AAS 09-224, "The Endgame Problem
Part A: V-Infinity Leveraging Technique and the Leveraging Graph," 2009 AAS/AIAA Space
Flight Mechanics Meeting.)

**[CR-PART2]**
S. Campagnola and R. P. Russell, "The Endgame Problem Part 2: Multibody Technique and
the Tisserand–Poincaré Graph," *Journal of Guidance, Control, and Dynamics*, Vol. 33,
No. 2, 2010, pp. 476-486, DOI 10.2514/1.44290. (Preprint AAS 09-227, "The Endgame
Problem Part B: The Multi-Body Technique and the T-P Graph," 2009 AAS/AIAA Space Flight
Mechanics Meeting.)

Both already cited (Part 2) in `src/cyclerfinder/search/tisserand.py` docstring.

---

## Paper 1 (Part A / Part 1): V∞-Leveraging Technique and the Leveraging Graph

### What the paper is

Derives new closed-form formulae for the V∞-leveraging maneuver (VILM), builds the
**Leveraging Graph** as a reference tool for endgame/begingame design, proves that VILM
sequence cost decreases when favouring high-altitude flybys, and gives a **quadrature
formula for the theoretical-minimum ΔV** of a multi-VILM transfer between moons. Worked
on the Jupiter (Galilean) and Saturn moon systems.

Abstract claim (p.1): "we find a simple quadrature formula to compute the minimum DV
transfer between moons using VILMs, which is the main result of the paper. The Leveraging
Graphs and associated formulae are derived in canonical units and therefore apply to any
celestial system with a smaller body in a circular orbit around a primary."

### Model assumptions (applicability limits)

p.2: "VILMs are typically modeled in the *linked-conics model* (or zero radius sphere of
influence, patched-conics model) where the minor body is considered massless and is on a
circular orbit around the major body. The spacecraft trajectory is coplanar and starts and
ends at the minor body. The gravity assist is modeled as an instantaneous change in the
direction of the relative velocity V∞ by the deviation angle δ."

p.3 (assumptions baked into the VILM model): "We assume that the impulsive maneuver is
tangential and is performed exactly at the apses. This assumption is typically included
when studying VILMs because the Jacobi constant in the rotating frame is maximally changed
by performing the maneuver when the rotating velocity is the greatest - this occurs at
apses. We also assume that the spacecraft departs/arrives at point L tangent to the minor
body. This condition guarantees the lowest V∞L and greatly simplifies the tour problem
because we can decouple each VILM as opposed to having to optimize a large sequence of
VILMs altogether."

So: **coplanar, circular-orbit minor body, patched-/linked-conics, tangential impulse at
apses.** This is the same regime as our coplanar `linkable` (i=0). Note the leveraging-graph
ΔV formulae are linked-conic; the CR3BP refinement is Part 2.

### Adimensional variables (p.3)

Length scale `l_scale = ã_M` (minor-body SMA); time scale `t_scale = sqrt(ã_M^3 / μ̃_P)`;
velocity scale = minor-body velocity Ṽ_M. In adimensional units the minor body's SMA,
velocity, and the major-body GM are unity.

"V_c as the adimensional velocity of the circular orbit of radius r̃_π = r̃_M + h̃_π around
the minor body: V_c = sqrt(μ̃_M / r̃_π) / Ṽ_M" — V_c "groups the problem dependency on the
minor body gravity constant, minor body radius and altitude of the final/initial orbit
insertion/escape." (p.3)

### Key VILM equations (quoted by number)

- Gravity-assist deflection (adimensional), p.3:
  "δ = 2 arcsin(1 / (1 + (V∞/V_c)^2))"
- Pericenter speed at the minor body, **Eq. (1)**, p.3:
  "V_π(V∞, h_π) = sqrt(V∞^2 + 2 V_c^2)"
- VILM boundary speeds, **Eqs. (2)–(3)**, p.4:
  "V_L^(E,I) = 1 ± V∞L"  (2)
  "V_A^(E,I) = V_B ± ΔV_AB"  (3)
  "where the upper sign refers to the Exterior VILM and the lower sign refers to the Interior
  VILM." Boundary values: "0 < V∞L < sqrt(2) − 1 for the Exterior VILM for r_A to be
  bounded, and 0 < V∞L < 1 for the Interior VILM for V_L to be positive."
- The Γ function (phase-free), p.5, defined fully at **Eq. (25)** (Appendix A), p.19:
  "Γ^(E,I)(V∞L) ≡ ±(r_A − V_A) = V∞L (V∞L^3 ± 3 V∞L^2 − V∞L ∓ 7) / (V∞L^3 ± 3 V∞L^2 +
  V∞L ∓ 1)"  — "Γ^(E) is computed for the Exterior VILM, and Γ^(I) for the interior VILM…
  Γ is a positive strictly monotonic function of V∞."
- **The two central phase-free formulae**, p.5:
  - **Eq. (4)**: "V∞H(V∞L, ΔV_AB) = sqrt((V∞L)^2 + (ΔV_AB)^2 + 2 ΔV_AB Γ)"
  - **Eq. (5)**: "ΔV_AB(V∞L, V∞H) = −Γ + sqrt(Γ^2 + (V∞H^2 − V∞L^2))"
- Vis-viva / relative-velocity identities, **Eq. (10)**, p.7:
  "V∞H^2 = V_πH^2 − 2 V_c^2 ,  V∞L^2 = V_πL^2 − 2 V_c^2"

### VILM taxonomy (the n:m±(K) classification)

This is the n:m taxonomy the prompt asked for; quoted from pp.3-4.

Four types by two binary features (p.4):
- "*Forward (Backward)* if the ΔV_AB is in the same (opposite) direction of the spacecraft
  velocity."
- "*Exterior (Interior)* if the ΔV_AB occurs at apocenter (pericenter), thus if r_A > a_M
  (r_A < a_M)."
- "Forward-Exterior V∞−leveraging and Backward-Interior V∞−leveraging *decrease* the V∞,
  while the Forward-Interior V∞−leveraging and Backward-Exterior V∞−leveraging *increase*
  the V∞." (p.4)

Per-VILM labels (p.4):
- "the resonant ratio: n : m, where n (m) is the approximate number of the minor body
  (spacecraft) revolutions during the VILM."
- "K number of full revolution in the arc H − B."
- "the point H− or H+ where the spacecraft encounters the minor body, resulting in a
  *long-transfer* VILM or *short-transfer* VILM respectively. Exterior, long-transfer VILMs
  and Interior, short-transfer VILMs are linked by prograde gravity assists. Exterior,
  short-transfer VILMs and Interior, long-transfer VILMs are linked by retrograde gravity
  assists."

Notation convention adopted in the paper (p.4): "we refer to 'Backward/Forward,
Interior/Exterior n : m_K±' VILMs. For example the Europa endgame when approached from
Ganymede is a sequence of Forward Exterior VILMs."
Footnote, p.4: "In literature we can find a different choice of letters: K : L(M)± where
K ≡ n, L ≡ m, and M ≡ K for Exterior VILM."

Example tokens shown in Fig. 3 caption (p.4): "5 : 4_2^+" (ΔV_AB after two full revolutions
on leg H+−B; transfer lasts a bit more than 5 minor-body revolutions) and "5 : 4_1^−"
(ΔV_AB after one full revolution on leg H−−B; lasts a bit less than 5 minor-body
revolutions).

> Cross-reference to our schema: our `free_return_arcs[].resonance` field stores M:N tokens
> like `F(3:2,82.487,180.000)` (Russell descriptor, per catalogue.schema.json). That n:m is
> the *same kind of object* as the VILM resonant ratio here (body-revs : spacecraft-revs),
> but the source families differ — our descriptor is from the McConaghy/Russell cycler work,
> not from this endgame paper. See "what this does NOT give us" below.

### Efficiency definitions

- Phase-free efficiencies, **Eqs. (6)–(8)**, p.5 (ε_BI−FE and ε_BE−FI).
- Phase-fixed efficiency, p.7: "E = (V∞H − V∞L) / ΔV_AB".
- **Theorem (p.11)**: "The total ΔV of a sequence of VILMs decreases when favoring VILMs
  with high altitude gravity assists." (Proven via signs of the efficiency derivatives.)

### Minimum-V∞L proposition (when is a VILM worth doing)

p.7 **Proposition**: "The VILM strategy is efficient iff V∞L > V̄∞ where V̄∞ =
sqrt(V̄_π^2 − 2 V_c^2) and V̄_π(V_c) is the root of the function: f(V_π) = Γ ∘ V∞(V_π;V_c)
− V_π" (**Eq. (9)**).
p.8 cubic-spline approximations of V̄∞(V_c) for Exterior and Interior cases are given (two
cubics in V_c).

### Theoretical-minimum-ΔV quadrature (the headline result)

p.13 **Eq. (12)**: "ΔV_min^(E,I)(V∞L, V∞H) = ∫_{V∞L}^{V∞H} V∞ / Γ^(E,I)(V∞) dV∞"
p.13 **Eq. (13)** (same integral with Γ expanded):
"ΔV_min^(E,I) = ∫ (V∞^3 ± 3V∞^2 + V∞ ∓ 1) / (V∞^3 ± 3V∞^2 − V∞ ∓ 7) dt"
(integrated V∞L→V∞H; "the integral can be solved numerically with quadrature or with
partial fractions"). Domain: "0 ≤ V∞ ≤ sqrt(2) − 1 for the Exterior VILM, and
0 ≤ V∞ ≤ 1 for the interior VILM."
p.13 **Eq. (14)** (the maximum, = single VILM): "ΔV_max(V∞L, V∞H) = −Γ +
sqrt(Γ^2 + (V∞H^2 − V∞L^2))".

### Algorithmic content (Forge-relevant) — Part 1

- **Phase-fixed solver** (p.5, "the algorithm to compute the numerical solution to the
  constrained problem is a variation of the one described for instance in [3], so we skip
  the details"): guess flight-path angle γ at H → orbital params of leg H−B → orbital params
  of leg L−A → transfer time and t_L → "differentially correct the path angle γ until the
  distance vanishes." A shooting/differential-correction loop, not graph search.
- **Branch & bound over the Leveraging Graph** (pp.9-10): "This recursive strategy is
  well-suited for a *branch&bound* search, because starting from a fixed Ṽ∞INITIAL the
  algorithm recursively applies Forward-Exterior VILMs and stores the ToF and total cost of
  the Endgame." (p.10) "The branch and bound solutions from Figure 9 on the right agree
  qualitatively with those from [23] that are found using an enumerative method based on
  dynamic programming principles." (p.10) — explicit pointer that the Pareto front (ToF vs
  total ΔV) can also be produced by **dynamic programming**.
- **Minimum-ΔV via quadrature** (pp.12-14): instead of searching, integrate Eq. (13).
  "Thus the cheapest way to move from an initial to a final V∞ is by zigzagging 'low' on the
  x−axis." (p.11) For moon-to-moon (p.14): "the logical strategy for the minimum ΔV transfer
  consists of a sequence of Interior VILMs at M2, followed by the Hohmann transfer, and
  finally a sequence of Exterior VILM at M1." Use Eq. (13) twice + Eq. (15) for the Hohmann
  V∞.

---

## Paper 2 (Part B / Part 2): Multibody Technique and the Tisserand–Poincaré (T-P) Graph

### What the paper is

Introduces the **Tisserand–Poincaré (T-P) graph**, a CR3BP extension of the classic
Tisserand graph, to design *ballistic / quasi-ballistic* endgames and intermoon transfers.
Core results: (a) ballistic intermoon transfers are energetically possible at low energy
even though linked-conics says they are not; (b) the T-P graph supplies a target patch
point for begingame+endgame; (c) two worked Ganymede–Europa transfers cheaper than the
VILM minimum.

### How the T-P graph differs from the classic Tisserand graph

p.7: classic Tisserand graph plots pericenter r_p and period T (or, as used here, an
apocenter–pericenter r_a–r_p plane) of *Keplerian coplanar* orbits — it is a **linked-conic**
construct. The T-P graph is built from **Poincaré sections of CR3BP** trajectories:

p.2: "we introduce a Poincaré section in the negative x-axis of the rotating reference frame
of the CR3BP. Far from the minor body the spacecraft trajectory is very similar to a
Keplerian orbit; thus we can compute the osculating orbital elements of the spacecraft as it
crosses the section, and plot them in a pericenter vs apocenter graph. On the same graph we
plot Tisserand parameter level sets… The result is the Tisserand-Poincaré (T-P) graph, which
is a natural extension of the v∞ level sets of the Tisserand graph, as the v∞ level sets are
synonymous to Tisserand parameter level sets noting that T = 3 − v∞^2 (see appendix). Yet the
Tisserand level sets extend beyond the v∞ curves well into the regions where v∞ is not
feasible in the linked-conic models (v∞^2 < 0 if T > 3). Therefore… the T-P graph
demonstrates that ballistic transfers between moons are energetically possible despite the
contrary conclusion derived from linked-conics theory. This is the first important result of
the T-P graph."

**Axes**: apocenter (x) vs pericenter (y), i.e. r_a–r_p, *not* energy-vs-rp. (Figs. 3, 4, 8,
11, 13; p.7-8.) Period level sets are diagonal lines; resonance lines have slope −1.

### Tisserand-parameter equations (T-P graph)

- **Eq. (10)**, p.7 (general): "T(a,e,i) = 1/a + 2 sqrt(a(1−e^2)) cos i"  (normalized
  a = ã/ã_M).
- **Eq. (11)**, p.7 (planar, in r_a, r_p): "T(r_a, r_p) = 2/(r_a + r_p) +
  2 sqrt(2 r_a r_p / (r_a + r_p))".
- **Eq. (12)**, p.8 (w.r.t. minor body M, dimensional inputs): "T_M = 2ã_M/(r̃_a + r̃_p) +
  2 sqrt(2 r̃_a r̃_p / ((r̃_a + r̃_p) ã_M))".
- Resonance lines, p.8: "a = (n/m)^(2/3) → r_p = −r_a + 2(n/m)^(2/3)" (n body revs,
  m spacecraft revs; lines of slope −1).
- **Tisserand ↔ V∞ ↔ Jacobi identity**, p.8: "T = 3 − v∞^2 ≈ J"; derived in full at
  **Eq. (42)**, p.20: "J ≈ T = 3 − v∞^2". This is *exactly* the relation our
  `vinf_to_tisserand` implements (`T_p = 3 − V∞^2 a_p/μ_sun` in dimensional form).
- **3-D T-P**, **Eq. (14)**, p.9: "T_M = 2ã_M/(r̃_a + r̃_p) + 2 sqrt(2 r̃_a r̃_p /
  (r̃_a + r̃_p)) cos i". (Used to visualise asteroid families / Solar Orbiter; Fig. 6.)

### Regions of motion (applicability structure)

p.8-9: Tisserand level sets T = J_Li (i=1..4; note J_L4 = J_L5 = 3) divide the r_a–r_p plane
into regions. "Transfers to the minor body are possible only when the spacecraft is in the
regions II^i, II^e, III. In particular, we expect low-energy transfer and capture
trajectories to occur in the region II^i (if coming from the inner moons) or II^e (if coming
from the outer moons)." Patch point between two moons = intersection of their two T_M level
sets, **Eq. (13)**, p.9 (solve the 2×2 system for r_a, r_p).

### The multibody / ballistic-endgame mechanism (the "ballistic endgame paradox")

p.11 **Paradox statement** (verbatim, displayed): "Given a fixed Tisserand energy and arrival
circular orbit altitude, the insertion maneuver costs remain essentially fixed for all
possible arrival geometries."

p.12 resolution math: orbit-insertion ΔV at angle θ, **Eq. (19)**: "Δv_π = v_π − v_c =
V + σ r_π − v_c". "(V^2)_MAX − (V^2)_MIN = r_π^2 (1−μ) (3+r_π)/(1+r_π), and because r_π^2 is
small compared to other terms… V_MAX ≈ V_MIN, i.e. the velocity and thus the orbit insertion
maneuver doesn't depend significantly on the angle θ." Table 1 (see anchors below) shows
"the difference in the orbit insertion maneuver is just a few meters per seconds or less."

p.13 resolution (why resonant orbits are still needed): "if the spacecraft has the right
phasing, it can use Europa perturbing force to slightly lower its apocenter AND pericenter,
thus moving to the left in the T-P graph, along the level set. Such maneuver is in fact a
high altitude flyby performed close to the pericenter of the spacecraft orbit. When several
high altitude flybys are linked together by free-return orbits, the pericenter can be lowered
to the point where a 100 km approach at Europa is possible. Thus the high altitude flybys are
necessary *to reduce the pericenter*… while the resonant orbits simply provide a mechanism to
achieve multiple flybys."

> **Direct resonance/free-return overlap with our schema.** p.13: "the low-cost endgame
> orbit must be approximately *resonant*; whereas a non-resonant returns would necessarily
> have two intersection points between the spacecraft and minor body orbits." This is the
> same physical content as our `free_return_arcs[]` M:N resonant-return fields: a resonant
> (single-tangency) return vs a non-resonant (two-crossing) return. The "free-return orbits"
> linking high-altitude flybys are literally the begingame/endgame analogue of our
> Earth–Earth resonant arcs.

### CR3BP / patched-conic model assumptions (applicability limits)

- p.3: planar CR3BP, EOM **Eq. (2)**, Jacobi constant **Eq. (4)**:
  "J = 2Ω − V^2 = (X^2 + Y^2) + 2(1−μ)/R_1 + 2μ/R_2 + (1−μ)μ − V^2".
- p.4: "unless specified, we use the *planar* CR3BP, noting that the spacecraft trajectory
  is very close to the minor body orbital plane for the applications of interest here."
- p.4: "*patched CR3BP model*, i.e. we split the trajectory in phases where only one minor
  body at a time affects the motion… the boundary points of contiguous phases are patched
  together, sometimes using impulsive maneuvers."
- Tisserand ≈ Jacobi approximation validity, p.7: "increasingly accurate for smaller mass
  parameters μ and when the spacecraft is far from the minor body."
- VILM-vs-CR3BP discrepancy, p.6: "their total Δv might differ of as much as 10% when
  computed in a more accurate model"; p.6: "the cost of the VILM endgames can be off up to
  ±5% when compared to the more accurate CR3BP solutions" (footnote: "consistent with the
  ±10% difference observed during the design of the Cassini tour").

> Relevant to our catalogue's CR3BP-vs-patched-conic flag on moon-tour rows: this paper is
> the canonical statement that the *same* endgame can read 5–10% different ΔV between the two
> models, and that ballistic (Δv=0) intermoon transfers exist in CR3BP that linked-conics
> declares impossible. Any moon-tour row we tag "patched-conic" carries that 5–10% optimism.

### Algorithmic content (Forge-relevant) — Part 2

- **Flyby reproduction in CR3BP** (pp.4-5): closed-form map from linked-conic flyby params
  (r_π, σ, v∞, γ) to the CR3BP state at closest approach, **Eqs. (5)–(9)**. Lets a
  linked-conic VILM seed a CR3BP trajectory.
- **Endgame optimization** (pp.5-6): VILM solution → first guess → multiple-shooting in
  CR3BP, "nonlinear parameter optimization problem… control variables are the times,
  altitudes r_π−r̃_M, speeds V_π and angles θ of all the closest approaches and the times of
  the mid-course maneuvers", solved with Matlab `fmincon`. Caveats (p.6): "quasi-ballistic
  endgames cannot be found by simply designing a VILM endgame and optimizing it in the
  CR3BP… the multi-resonant transfers are chaotic in nature where the design space is plagued
  by multiple local minima that can easily trap gradient based optimizers… Instead we should
  seek solutions that start in the correct basin!"
- **T-P-graph transfer search (the closest thing to a Forge pathfinder)**, pp.13-14: "we
  implement a simple search to find such trajectories". Procedure (p.14): fix escape energy
  (region II^e) and capture energy (region II^i); from Jacobi constants get pericenter
  velocities via Eq. (18) and insertion/escape costs via Eq. (19); "We then scan the angles
  θ_Ga, propagate the initial conditions and store the transfers that decrease the pericenter
  the most in the shortest time. We also scan the angles θ_Eu, propagate *backwards* the
  initial conditions and store the transfers that increase the apocenter the most in the
  shortest time. In both the forward and backward propagations we have a precalculated target
  value for r_a and r_p respectively - from the intersection point in the TP graph - found
  from the solution to Eq. (13)." Then keep the Pareto front (p.14): "we only plot the set of
  Pareto-optimum points (shortest time, highest apocenter) reached by all the solutions."
- **Explicitly NOT manifold-based** (p.2 abstract; p.9): "the patch point calculation is
  straight forward and does not require the often tedious generation of manifolds and their
  associated intersections." This is a forward/backward greedy scan toward a precomputed
  T-P-graph target point, keeping a Pareto front — not Dijkstra, not full enumeration, not
  invariant-manifold intersection. Future work (p.16): "methods that allow specification of
  the resonant paths to reduce the computation requirements and provide more systematic
  searches."

---

## Candidate Golden Anchors

Only values transcribed directly from the papers' tables/text appear here. None were
computed by me. These are linked-conic / CR3BP *moon-system* numbers — appropriate for a
future VILM/T-P implementation, NOT for the Earth–Mars cycler catalogue.

### A1 — Part 1, Table 3 (p.17): Moon physical data used for the ΔV computations

Verbatim (columns: Moon, μ̃_M [km³/s²], ã_M [km], Ṽ_M [km/s], r̃_π [km], V̄∞ E/I [km/s]):

| Moon | μ̃_M | ã_M (km) | Ṽ_M (km/s) | r̃_π (km) | V̄∞ E/I (km/s) |
|---|---|---|---|---|---|
| Io | 5960 | 421800 | 17.330 | 1922 | 0.351 / 0.368 |
| Europa | 3203 | 671100 | 13.739 | 1661 | 0.277 / 0.290 |
| Ganymede | 9888 | 1070400 | 10.879 | 2731 | 0.372 / 0.404 |
| Callisto | 7179 | 1882700 | 8.203 | 2510 | 0.328 / 0.361 |
| Enceladus | 7 | 238040 | 12.624 | 352 | 0.029 / 0.029 |
| Tethys | 41 | 294670 | 11.346 | 633 | 0.052 / 0.052 |
| Dione | 73 | 377420 | 10.025 | 662 | 0.067 / 0.068 |
| Rhea | 154 | 527070 | 8.484 | 864 | 0.085 / 0.087 |
| Titan | 8978 | 1221870 | 5.572 | 4076 | 0.283 / 0.321 |

(Source for the physical data, per the table footnote p.17: http://ssd.jpl.nasa.gov/.
Circular orbits at 100 km altitude except Titan at 1500 km.)

> Best golden-anchor candidates here: V̄∞ E/I — the minimum V∞ at which a VILM becomes
> efficient — is a *derived* quantity (root of Eq. (9) for the listed V_c), so it doubles as
> a check on a future `min_vinf_for_vilm(moon)` implementation. ã_M and Ṽ_M are inputs.

### A2 — Part 1, Table 1 (p.16): Min/Max ΔV for VILM intermoon transfers (no gravity assists)

Min ΔV assumes infinite transfer time (sum of escape+begingame+endgame+capture);
Max ΔV = the Hohmann transfer without VILMs. Initial/final circular orbits at 100 km
(Titan at 1500 km). Verbatim (km/s):

| Transfer | ΔV_min | ΔV_max | ΔV_escape | ΔV_begingame | ΔV_endgame | ΔV_capture |
|---|---|---|---|---|---|---|
| Callisto-Ganymede | 1.81 | 2.13 | 0.73 | 0.13 | 0.13 | 0.81 |
| Callisto-Europa | 1.94 | 3.75 | 0.73 | 0.3 | 0.31 | 0.59 |
| Callisto-Io | 2.43 | 6.00 | 0.73 | 0.46 | 0.48 | 0.75 |
| Ganymede-Europa | 1.71 | 2.18 | 0.82 | 0.14 | 0.16 | 0.59 |
| Ganymede-Io | 2.3 | 4.38 | 0.82 | 0.36 | 0.37 | 0.75 |
| Europa-Io | 1.76 | 2.54 | 0.6 | 0.21 | 0.2 | 0.75 |
| Titan-Rhea | 1.15 | 2.19 | 0.64 | 0.15 | 0.18 | 0.18 |
| Titan-Dione | 1.28 | 3.33 | 0.64 | 0.23 | 0.27 | 0.14 |
| Titan-Tethys | 1.37 | 4.31 | 0.64 | 0.29 | 0.33 | 0.11 |
| Titan-Enceladus | 1.43 | 5.27 | 0.64 | 0.33 | 0.4 | 0.06 |
| Rhea-Dione | 0.52 | 1.12 | 0.18 | 0.10 | 0.10 | 0.14 |
| Rhea-Tethys | 0.66 | 2.3 | 0.18 | 0.19 | 0.19 | 0.11 |
| Rhea-Enceladus | 0.78 | 3.53 | 0.18 | 0.27 | 0.27 | 0.06 |
| Dione-Tethys | 0.42 | 0.97 | 0.14 | 0.08 | 0.09 | 0.11 |
| Dione-Enceladus | 0.55 | 2.19 | 0.14 | 0.17 | 0.18 | 0.06 |
| Tethys-Enceladus | 0.34 | 1.00 | 0.11 | 0.08 | 0.09 | 0.06 |

### A3 — Part 1, Table 2 (p.16): Min/Max ΔV for VILM transfers WITH intermoon gravity assists

Same ΔV_min/ΔV_max convention; transfer routed through intermediate moon(s). Verbatim (km/s):

| Transfer | ΔV_min | ΔV_max | ΔV_escape | ΔV_begingame | ΔV_endgame | ΔV_capture |
|---|---|---|---|---|---|---|
| Callisto-G-Europa | 1.61 | 2.07 | 0.73 | 0.13 | 0.16 | 0.59 |
| Callisto-G-E-Io | 1.81 | 2.35 | 0.73 | 0.13 | 0.2 | 0.75 |
| Ganymede-E-Io | 1.91 | 2.45 | 0.82 | 0.14 | 0.2 | 0.75 |
| Titan-R-Dione | 1.03 | 1.55 | 0.64 | 0.15 | 0.099 | 0.14 |
| Titan-R-D-Tethys | 0.98 | 1.47 | 0.64 | 0.15 | 0.086 | 0.11 |
| Titan-R-D-T-Enceladus | 0.93 | 1.5 | 0.64 | 0.15 | 0.086 | 0.061 |
| Rhea-D-Tethys | 0.47 | 1.04 | 0.18 | 0.097 | 0.086 | 0.11 |
| Rhea-D-T-Enceladus | 0.43 | 1.07 | 0.18 | 0.097 | 0.086 | 0.061 |
| Dione-T-Enceladus | 0.37 | 1 | 0.14 | 0.084 | 0.086 | 0.061 |

### A4 — Part 2, Table 1 (p.13): Max/Min orbit-insertion ΔV (the ballistic-endgame paradox)

"The maximum and minimum orbit insertion maneuver (m/s) for given altitudes and Jacobi
constant at Europa and at Titan." Cells are ΔvMAX/ΔvMIN (km/s) verbatim:

| Moon | J_L1@100km | J_L4@100km | J_L1@1000km | J_L4@1000km |
|---|---|---|---|---|
| Europa | 421.1/420.1 | 606.5/605.5 | 276.7/273.7 | 513.7/511.1 |
| Titan | 668.6/668.5 | 766.5/776.4 | 553.7/553.5 | 667.6/667.5 |

> Anchor value: the paradox claim is the tiny max−min spread (≈1 m/s at J_L1@100km for both
> moons). Note the Titan J_L4@100km cell reads "766.5/776.4" — MAX < MIN, which is almost
> certainly a transcription/typesetting inversion in the source (every other cell has
> MAX ≥ MIN). Flag, do not "correct".

### A5 — Part 2, worked-transfer scalar anchors (text, pp.13-16)

- Ganymede→Europa, 100 km circular both ends (p.13): "A direct Hohmann transfer would
  require only a few days, but would cost 2.18 km/s. A VILM strategy can reduce this Δv to a
  theoretical minimum of 1.71 km/s." (Note: 1.71 matches Ganymede-Europa ΔV_min in A2.)
- Same, T-P graph quasi-ballistic design (pp.13-14): escape "Δv_Escape ≅ 0.72 km/s",
  capture "Δv_Capture ≅ 0.51 km/s"; patch "a Δv of some 70 m/s to patch the 82 day begingame
  with a 252 day endgame"; "(291 + 82) days to transfer"; total "Δv_TOT ≈ 1.25 km/s, almost
  500 m/s less than the VILM theoretical minimum Δv, and almost 1 km/s less (but one year
  more!) than the direct Hohmann transfer." (p.15)
- Halo-to-halo Ganymede→Europa (pp.15-16): Jacobi constants "Jupiter-Ganymede CR3BP of
  J=3.0052" and "Jupiter-Europa CR3BP of J=3.0023"; target intersection
  "(r_a,r_p) = (694641 km, 1021834 km)" from Eq. (13); "the whole transfer takes less than
  300 days"; "3 high altitude flybys of Ganymede and 6 high altitude flybys of Europa";
  Ganymede begingame "takes 191 days", Europa endgame "lasts 110 days".

### A6 — Part 1, worked Europa endgame scalar anchors (text, pp.6, 10, 11-12)

- 3:2 Exterior VILM example (p.6): "the 3 : 2 Exterior VILM reduces the V∞ from V∞H = 0.131
  to V∞L = 0.1135… we estimate the ΔV_AB ≈ 0.0022. For a VILM at Europa, we multiply these
  values by the average velocity of Europa of approx 13.7 km/s to find that we decrease
  V̄∞H = 1.8 km/s to V̄∞L = 1.56 km/s using approx 30 m/s."
- Europa endgame, 3-VILM design (pp.10-11): total transfer time 46 days, total ΔV 154 m/s,
  reduces V̄∞ from 1.8 km/s to 0.77 km/s; alternative 14-VILM design "each using 10m/s for a
  total of 140 m/s" (cheaper, much longer).
- (Part-2 cross-check, p.6:) the same two endgames re-optimized in CR3BP: long-transfer
  147 m/s (VILM 154), short-transfer 165 m/s (VILM 155).

---

## Forge-relevant algorithms (consolidated)

What these two papers actually offer a graph-search trajectory finder:

1. **Two graphs, one identity.** Both reduce to the V∞↔Tisserand relation we already
   implement (`T = 3 − v∞^2`, Part 2 Eq. (42); our `vinf_to_tisserand`). The Leveraging
   Graph adds ΔV as a third axis on top of the V∞ level sets; the T-P graph adds CR3BP
   reachability (Tisserand level sets that extend past the linked-conic feasible region).
2. **Branch & bound on the Leveraging Graph** (Part 1, p.10) with an explicit note that a
   **dynamic-programming enumeration** ([23]) gives the same Pareto front (ToF vs ΔV). This
   is the cleanest "graph search over Tisserand space" statement in either paper and the
   most directly Forge-shaped.
3. **Quadrature theoretical-minimum ΔV** (Part 1, Eq. (13)) — a closed-form floor for any
   VILM tour leg, computable without search. Good as a *bound/heuristic* inside a Forge
   search (admissible lower bound for A*-style pruning), and as a fast feasibility screen.
4. **T-P-graph forward/backward greedy scan to a precomputed patch point** (Part 2, p.14):
   not Dijkstra/A*, not full enumeration, not manifold intersection — a Pareto-front-keeping
   shooting search seeded by the T_M1=T_M2 intersection (Eq. (13)). Explicitly avoids
   manifolds (p.9). Future-work admits it is not yet "systematic" (p.16).
5. **The linkable-pair idea is theirs too.** Part 2 Eq. (13) (intersect two moons' Tisserand
   level sets) is the multi-body generalization of our coplanar `linkable` contour
   intersection — same predicate, different feasibility region (T can exceed 3 in CR3BP).

---

## What this does NOT give us

- **No cyclers.** Neither paper mentions cyclers, Earth–Mars transfers, or any heliocentric
  cycler. The shared DNA is the resonant-return / V∞-preserving machinery, not cycler orbits.
  Do not add any catalogue rows from these papers.
- **No Earth–Mars / heliocentric numbers.** Every tabulated number is for Jupiter or Saturn
  moons (plus an asteroid-family illustration, Fig. 6). Our V∞-matching for the Earth–Mars
  catalogue gains nothing numerically here.
- **The n:m here is body-revs:spacecraft-revs for a single VILM**, *within one moon's
  endgame*. Our `free_return_arcs[].resonance` n:m is for Earth–Earth resonant return arcs
  of a heliocentric cycler. Same algebraic object (resonance ratio), different physical
  layer — do not conflate the tokens when populating the schema.
- **No verification of S1L1.** Nothing here touches McConaghy 2006 or S1L1 elements.
- **Tables 1–3 are golden-anchor candidates only for a *future* VILM/T-P module**, which we
  do not yet have. They cannot validate the current `tisserand.py` directly because that
  module computes contour *intersection* (a boolean), not ΔV. To use A1's V̄∞ E/I or A2/A3's
  ΔV_min as goldens we would first need to implement Eq. (9) (V̄∞ root) and Eq. (13)
  (quadrature). Logged as candidates, not wired in.
- **Model caveat for any future use.** All Part-1 ΔV are linked-conic; Part-2 shows the same
  endgame can differ 5–10% in CR3BP and that ballistic (Δv=0) intermoon transfers exist that
  linked-conics calls impossible. Tables 1–3 are therefore *upper-bound-ish* linked-conic
  references, consistent with our moon-tour rows' CR3BP-vs-patched-conic flag.
- **Two suspected source typos preserved verbatim** (do not silently fix): Part-2 Table 1
  Titan J_L4@100km cell "766.5/776.4" (MAX < MIN); and the Part-1 down-stream note is not
  affected. Both flagged in-place above.

---

## Summary

| Paper | Biggest finding | Cycler data? | Golden anchors | Forge algorithm |
|---|---|---|---|---|
| Part 1 (AAS 09-224) | Closed-form VILM ΔV (Eqs. 4-5) + quadrature theoretical-min ΔV (Eq. 13); n:m_K± taxonomy; branch&bound on Leveraging Graph | No | Tables 1, 2, 3 (moon ΔV + V̄∞) | Branch&bound / DP over Leveraging Graph; quadrature min-ΔV bound |
| Part 2 (AAS 09-227) | T-P graph (CR3BP extension of Tisserand graph, r_a-r_p axes); T=3−v∞²≈J; ballistic-endgame paradox; ballistic intermoon transfers feasible | No | Part-2 Table 1 (insertion ΔV); A5/A6 scalar transfer costs | Forward/backward greedy scan to T_M1=T_M2 patch point (Eq. 13), Pareto front; explicitly no manifolds |
