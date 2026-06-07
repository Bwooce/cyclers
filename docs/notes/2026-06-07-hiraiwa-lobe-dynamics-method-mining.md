# Hiraiwa et al. 2026 — lobe-dynamics directed-graph sequencing (method mining)

Mined 2026-06-07 (Task #142). The lobe-dynamics theory (pp.3-16) + the
weighted-directed-graph sequencing (pp.19-25), deliberately skipped before.
Maps to #76 Tier 2 / low-energy transfer sequencing.

**Source (cite exactly, no file path):**
Hiraiwa, N., Bando, M., Sato, Y., & Hokamoto, S., "Design of low-energy transfers
in cislunar space using sequences of lobe dynamics," preprint submitted to Acta
Astronautica (arXiv:2602.17444v3, 13 May 2026). Kyushu Univ. / Hokkaido Univ.
54 pp.

> Clean digital typeset; equations and Table 1 read unambiguously. Vision read of
> pp.1-25 (the theory + graph method + test-problem results).

---

## 1. The algorithm in 3 lines

In the planar CR3BP, the stable/unstable manifolds of **unstable resonant orbits**
form **lobes** on a periapsis Poincaré map (Delaunay `g_d`-`G_d` plane); a
**lobe sequence** is the forward/backward `F`-images of one lobe; keep only
**effective lobes** (radius `r_{L,i} > r_L*`) for robustness. Build a **weighted
directed graph**: nodes = {start, effective-lobe centroids, goal}, edges =
allowable transfers (within a sequence = near-free natural dynamics; between
sequences/from start = a small **targeting Δv**), weights = Δv (or ToF). Find the
fuel-optimal path by **exhaustive search** minimizing `J = Σ w_k` s.t. each edge
`w_k < w*` and "use lobe sequences properly" (≥2 adjacent lobes per sequence
visited). Mapping verdict: **the directed-graph-of-natural-arcs + targeting-Δv-
edges + threshold-pruned exhaustive search MAPS to our sequence enumeration /
Tisserand-graph (#76 Tier 2) as a low-energy multi-body analogue; the CR3BP/lobe
dynamics itself does NOT map to our heliocentric patched-conic cycler model.**

---

## 2. The method stated precisely

### 2.1 Setup — CR3BP, resonant orbits, periapsis Poincaré map (Ch.2)
- CR3BP rotating frame, pseudo-potential `U` (Eqs.1-4), Jacobi constant
  `C_J = 2U − (ẋ²+ẏ²+ż²)` (Eq.5); zero-velocity curve at `2U−C_J=0` divides
  interior/moon/exterior realms (Eq.6).
- `p:q` **resonant orbit**: `pT_3^I ≃ qT_2`, `a ≃ (q/p)^{2/3}(1−μ)^{1/3}` (Eq.7).
  Earth-Moon μ = 1.21509e-2.
- **Periapsis Poincaré map** (Ch.2.3): surface of section at periapsis passage
  (`l_d = 0`, `XẎ−YẊ>0`), plotted in canonical **Delaunay elements**
  `(l_d, g_d, L_d, G_d)` (Eqs.12-16): `l_d=M`, `g_d=ω`, `L_d=√((1−μ)a)`,
  `G_d=√((1−μ)a(1−e²))` (= angular momentum). Map drawn in the `g_d`-`G_d` plane
  at fixed `C_J` (=3.16) → divides into **resonance regions** (quasi-periodic/
  periodic) and **chaotic zones** (Fig.5).
- BCR4BP (Ch.2.4, Eqs.17-21) used only for the final Sun-perturbed comparison.

### 2.2 Lobe dynamics theory (Ch.3, Defs 1-7)
On a 2-D area/orientation-preserving map `F` (here the Poincaré map), with
unstable manifold `W^u_{p1}` and stable manifold `W^s_{p2}` forming a hetero/
homoclinic tangle:
- **Def 1 heteroclinic/homoclinic point**: `q_i ∈ W^u_{p1} ∩ W^s_{p2}`.
- **Def 2 primary intersection point (pip)**: `q_i` such that the manifold
  segments `U[p1,q_i]` and `S[p2,q_i]` intersect only at `q_i`.
- **Def 3 lobe**: region bounded by `U[q_0,q_1]` and `S[q_0,q_1]` for adjacent
  pips q_0, q_1.
- Lobe transport governed by `F`: lobe `(U[q0,q1],S[q0,q1]) → (U[Fq0,Fq1],
  S[Fq0,Fq1])`. Lobes crossing resonance-region boundaries are **turnstiles**
  governing transport between resonance regions and chaotic zones.
- **Def 4 lobe sequence** `Λ(F,L) = {…, F^{-1}(L), L, F(L), …}`.
- **Def 5 radius** `r_{L,i}` = max ε-ball contained in lobe (≈ min distance
  centroid→boundary).
- **Def 6 effective lobe**: `r_{L,i} > r_L*` (threshold). **Def 7 effective lobe
  sequence**: finite subsequence of effective lobes. Rationale: transfers between
  lobe centroids tolerate control/OD errors smaller than `r_L*` and stay inside
  the lobe → **robust** point-to-point transfer. (`r_L* = 0.002` representative.)
- Construction (Ch.3.2): pick large-resonance unstable orbits (7:2, 3:1); grid-
  search Poincaré map → differential-correct exact ICs; propagate manifolds →
  heteroclinic tangle → parametrize manifolds → find pips → lobes → propagate
  one lobe ± in time to reveal the whole sequence.

### 2.3 The trajectory-design graph (Ch.4, Eqs.22-24) — the namesake sequencing
"a **weighted, directed graph** with three components: nodes (discrete state
points), directed edges (transfer paths), and weights (transfer costs)."
- **Nodes** = potential start/goal points + **centroids of effective lobes**.
- **Edges** = allowable transfer paths, predesigned from dynamical geometry;
  edges designed so the spacecraft **monotonically increases `G_d`** (climbs
  toward escape).
- **Weights** = transfer cost (Δv or ToF).
- **Combinatorial optimization** (Eqs.22-24):
  `minimize J = Σ_k w_k` s.t. `w_k < w*` (Eq.23, threshold-prune infeasible
  edges) and **"use lobe sequences properly"** (Eq.24 — if a lobe sequence is
  entered, ≥2 adjacent lobes of it must be in the path, forcing full exploitation
  of the natural lobe dynamics). Solved by **exhaustive search** over all
  start→goal paths. `w* = 100 m/s` (paths with any edge >100 m/s ignored).
- **Edge construction by targeting** (Ch.4.2, Fig.17): for a controlled transfer
  `x_1 → x_2`, propagate `x_1` forward over `[0, 2π]` and `x_2` backward over
  `[0, −2π]`; find a **crossing point in position space**; apply an impulsive
  **ΔV at the crossing** to adjust velocity direction while keeping `C_J`; edge
  weight = |ΔV|. (This is the classic "targeting" of two manifold legs to a
  position-space crossing.)
  - **Within a lobe sequence** edges use natural dynamics: Case 1 (no targeting)
    propagate `x_1` until it reaches the lobe of `x_2` (weight 0, `x_2` updated to
    the propagated endpoint); Case 2 (with targeting) use Fig.17 (small weights,
    since arcs follow lobe dynamics).

### 2.4 Test problem + results (Ch.4.2-4.3, Table 1)
- Departure = **7:2 stable resonant orbit** (periapses ~27,279-37,746 km alt ≈
  GEO); destination = left half of the **L1 Lyapunov stable manifold** (exit to
  moon realm / deep space). Goal node = centroid of the 8th effective lobe in
  sequence L2.6. Earth-Moon CR3BP, `C_J=3.16`.
- 8 effective lobe sequences (L1.1, L1.2, L2.1-L2.6) selected (Fig.10), `r_L*=0.002`.
- **Table 1 — number of transfers found**: Estimated / Without-targeting /
  With-targeting: without Eq.24 = **544 / 528 / 544**; with Eq.24 = **388 / 376 /
  388**. ("Estimated" assumes intra-sequence edge weights are exactly zero.)
- Targeting (Case 2) reconstructs ALL graph paths into real trajectories; Case 1
  fails for some. Transfer costs cluster ~150-420 m/s (Fig.18 detailed panels);
  Eq.24 prunes count but keeps low-cost solutions.

---

## 3. Maps to our X / does not map

| Hiraiwa element | Our code / concept | Verdict |
|---|---|---|
| Weighted directed graph of natural arcs; nodes=lobe centroids+start/goal; min Σw_k, threshold-prune w<w* (Eqs.22-24) | `search/sequence.py` / `search/scan.py` sequence enumeration; #76 Tier 2 low-energy sequencing | **MAPS (as analogue).** The directed-graph-of-dynamical-arcs + threshold-pruned exhaustive path search is exactly the structural pattern for low-energy multi-body sequencing. Our cycler sequence enumeration is the heliocentric flyby-chain version; this is the cislunar manifold version. |
| Targeting edge: propagate fwd [0,2π] + bwd [0,−2π], crossing in position space, impulsive ΔV at crossing keeping C_J (Fig.17) | (no manifold-targeting; our corrector matches V∞ at flybys) | **DOES NOT MAP — different domain.** Targeting two manifold legs to a position crossing is a CR3BP technique; our Lambert legs connect planet positions directly. Conceptually parallel to a deep-space-maneuver patch (Vasile BS3 DSM), though. |
| "Use lobe sequences properly" (≥2 adjacent lobes per sequence, Eq.24) | (no analogue) | **DOES NOT MAP** — specific to exploiting lobe turnstiles; no cycler equivalent. |
| Effective lobe / threshold radius r_L* for robustness | our verify/adversarial robustness gates | **PARTIAL (philosophy only).** The "stay-within-error-ball ⇒ robust transfer" idea parallels our adversarial/robustness checks, but the lobe-radius construction is CR3BP-specific. |
| Periapsis Poincaré map in Delaunay (g_d,G_d), resonance regions vs chaotic zones | `search/tisserand.py`, `search/resonance.py` | **PARTIAL — same conceptual space (resonance structure), different map.** Tisserand graph (v_rel vs period/pericenter) is the patched-conic analogue; the periapsis-Poincaré-Delaunay map is the CR3BP analogue. Both organize resonant transfers. |
| `p:q` resonant orbit, `a≃(q/p)^{2/3}` | `search/resonance.py` resonance ratios | **MAPS conceptually** — same resonance bookkeeping (cf. Vasile SOT ρ=n/m, Campagnola endgame). |
| CR3BP / BCR4BP dynamics, Jacobi constant, manifolds, lobes | our heliocentric patched-conic / Lambert model | **DOES NOT MAP** — fundamentally different dynamical model (rotating-frame 3-body vs Sun-centred 2-body patched conic). This paper is a *cislunar* method; our cyclers are *interplanetary*. The transferable content is the GRAPH/SEQUENCING layer, not the dynamics. |

---

## 4. Candidate test anchors

The only tabulation is **Table 1 (p.24): transfer counts** (544/528/544 without
Eq.24; 388/376/388 with Eq.24). These are combinatorial-search counts for a
specific Earth-Moon test problem (8 lobe sequences, C_J=3.16) — usable only as a
reproduction check IF we ever implement the exact lobe-graph (not applicable to
the cycler catalogue). No orbital-element or cycler tables.

Resonant-orbit anchors (closed-form): `a ≃ (q/p)^{2/3}(1−μ)^{1/3}` (Eq.7);
Earth-Moon μ=1.21509e-2; test orbits 7:2 and 3:1 at C_J=3.16. These are
CR3BP-internal, not cycler-catalogue rows.

> Nothing here is golden-eligible for our heliocentric cycler catalogue. The
> paper is a cislunar low-energy-transfer method; its value is methodological
> (the graph-sequencing pattern), not as data anchors.

---

## 5. Single most implementable finding (this paper)

**The weighted-directed-graph-of-natural-arcs with threshold-pruned exhaustive
path search (Eqs.22-24) as the sequencing layer for #76 Tier 2 low-energy
cyclers.** The reusable pattern: (1) discretize the solution space into nodes =
centroids of robust dynamical structures (for us: feasible flyby/resonance states
on the Tisserand graph rather than lobe centroids); (2) edges = predesigned
allowable transfers with a Δv (or ToF) weight; (3) prune any edge above a
threshold `w*` (their 100 m/s) to keep the combinatorics tractable; (4)
exhaustive (or Dijkstra) shortest-path for the fuel-optimal sequence. This is a
clean, domain-agnostic recipe that our `search/sequence.py` /`tisserand.py`
sequencing could adopt directly — replacing ad-hoc enumeration with a weighted
graph + threshold prune + shortest-path. (Same graph-of-arcs idea appears in
Ozimek LInX and Russell's enumeration; Hiraiwa gives the cleanest explicit
formulation with the `w<w*` prune.)

---

## 6. v4.2 backfill checks

- **center**: CR3BP rotating frame (Earth-Moon barycenter); periapses computed
  w.r.t. Earth (`X-Y` Earth-centered inertial frame, Eq.8). BCR4BP for the Sun-
  perturbed step. **Not heliocentric — no overlap with our cycler catalogue
  center conventions.**
- **tof_days_bounds**: transfer times in non-dimensional CR3BP units / lobe-map
  iterations; manifold propagation over `[0, ±2π]` (one synodic-ish period). Not
  convertible to cycler tof_days_bounds without the Earth-Moon timescale, and not
  relevant to any catalogue row.
- **source_ephemeris**: **none** — analytic CR3BP/BCR4BP (μ=1.21509e-2, Sun mass
  3.28914e5, Sun distance 3.89173e2 ND). No DE ephemeris. Any (hypothetical)
  catalogue note must record `model_assumption: CR3BP/BCR4BP (analytic, no
  ephemeris)`.

No catalogue rows trace to this paper (cislunar, not Earth-Mars cycler).

---

## 7. Honest "not extractable" list

- No Earth-Mars cycler content; entirely Earth-Moon cislunar.
- The BCR4BP optimal-transfer comparison (Ch.5-6, pp.26+) and the multiple-
  shooting CR3BP→BCR4BP conversion were NOT read this pass (the task scoped
  pp.3-32; the lobe theory + graph method through p.25 is the core). The §5-6
  BCR4BP results / literature comparison remain unmined if a low-energy follow-up
  wants them.
- Lobe radii are given as histograms (Fig.9) not tabulated values.
- Effective-lobe-sequence membership (which lobes, what G_d) is shown in Fig.10
  plots only, not tabulated.
