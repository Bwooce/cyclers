# Braik & Ross 2026 — "Orbital Networks in the Three-Body Problem" — mining note (#230)

**Source:** Abdullah Braik & Shane D. Ross, *Orbital Networks in the Three-Body
Problem*, arXiv:2605.31543 (v1, 29 May 2026; dated June 1, 2026). Virginia Tech.
**Lineage:** First forward citation of Ross & Roberts-Tsoukkas AAS 25-621 (same
lab). Cites it (and JPL periodic-orbit DB, Rawat et al. 2026) only as the *source
of the cycler representative orbits*.

Read in full (58 pp): abstract, Secs. 2–8, Tables 1–8, Figs. 1–18, Appendices A–C.

---

## What the paper does

Introduces a **reachable-set-based accessibility-network** framework for the
planar Earth–Moon CR3BP. Rather than designing point-to-point transfers, it asks:
across many coexisting periodic-orbit families *at one fixed Jacobi energy*, which
families are hubs / gateways / relays, and which are persistently hard to reach?

Pipeline (Secs. 2–5):

1. **Common-energy representative set.** Thirteen representative planar periodic
   orbits, one per family, all re-selected at a **single common Jacobi constant
   C_J = 3.1294** (so accessibility differences reflect geometry/transport, not
   energy mismatch). Families: LL1, LL2 (Lyapunov); C11a, C11b, C21, C32
   (cyclers); R21-S/U, R31-S/U, R52-S/U (2:1, 3:1, 5:2 resonant, stable+unstable);
   DPO (distant prograde). Table 2 lists period (days) and Floquet instability
   rate σ = ln(λ_max)/T_PO.
2. **Reduced (x,y,θ) model.** At fixed C_J, speed v = √(−2Ū−C_J) is determined by
   position, so the state collapses to (x, y, θ) with θ the velocity heading
   (Eq. 8). Time-reversal symmetry R(x,y,θ)=(x,−y,π−θ) (Eq. 14).
3. **Local reachable-set atlases (Sec. 3).** From arc-length-spaced seeds on each
   representative orbit, apply an **energy-preserving instantaneous heading-change
   maneuver** (a pure rotation of the rotating-frame velocity at fixed speed →
   stays on the same C_J manifold; ΔV_turn = 2v sin(|δ|/2), Eq. 26), propagate
   one-sided to horizon T_a, and log every crossed voxel on a (x,y,θ) grid.
   Backward reachable set obtained for free via the time-reversal symmetry.
4. **Pairwise overlap → proxy cost (Sec. 4).** Family A→B is "directly accessible"
   iff forward-reachable(A) ∩ backward-reachable(B) ≠ ∅ on the shared voxel grid.
   Each overlap voxel gets a **proxy ΔV** (sum of min source-side + min
   target-side turning costs + a voxel-scale heading-mismatch patch term) and a
   **proxy time**. The pair value is the minimum-proxy-ΔV voxel. This is an
   explicitly *necessary-not-sufficient*, *screening* metric — not a transfer
   solution.
5. **Cislunar orbital network (Sec. 5).** The 13×13 symmetric proxy-ΔV matrix is a
   weighted undirected graph. Cost-aware centralities: **strength** (hub,
   reciprocal-cost), **harmonic closeness** (gateway/staging), **betweenness**
   (relay). At the max-budget reference case (ΔV_cap = 409.3 m/s, T_cap = 27.32 d)
   the graph is nearly complete (75 of 78 pairs present).
6. **Budget-dependent regimes (Sec. 6).** Across the (ΔV_cap, T_cap) plane, direct
   accessibility, graph connectedness, and budget-feasible multileg closure
   separate into three distinct operational regimes (time-critical retasking /
   sustained coverage / time-flexible redistribution).
7. **Validation (Sec. 7).** Selected proxy-supported connections are refined into
   concrete patched trajectories by local differential correction (Appendix C);
   **corrected ΔV < proxy ΔV in every tested case** → proxy is conservative.

**Headline structural results:**
- **(3,2)-cycler (C32) is the dominant family** at max budget: rank 1 in strength,
  harmonic closeness, AND betweenness (Table 4; Fig. 10). Hub + gateway + relay.
- **(1,1)a-cycler (C11a) dominates the low-time-of-flight regime** and is the
  primary *relay* bridge in the two-leg examples.
- **2:1 stable resonant (R21-S) is the persistent hard-access family** (last in
  strength & closeness, zero betweenness, the only 3 missing direct edges all
  involve R21-S: LL1–R21-S, LL2–R21-S, DPO–R21-S).
- Dynamical reading: the 10 *unstable* representatives share a single connected
  chaotic region of the manifold (low-cost mutual access); the 3 *stable* resonant
  representatives sit inside resonant tori that resist low-energy invasion from
  the chaotic sea → the hard-access edge.

---

## Digest answers

### Q1 — Does it restate or EXTEND the AAS 25-621 cycler tuples? **NO new sourced tuples.**

**Verdict: N.** The paper does **not** reproduce or extend any C/T/IC tuples from
Ross & Roberts-Tsoukkas 2025. It *uses* the cycler families as named
representatives but:

- It re-selects each family's representative orbit at a **single common
  C_J = 3.1294**, which is a *different energy* from the per-family stable Jacobi
  values published in AAS 25-621 (and carried in our catalogue rows). Example: our
  `ross-rt-em-cycler-11-2025` is at C_J = 3.1511759 with T = 10.2920692 TU =
  44.7538 d (the family's nu≈0 stable midpoint). Braik-Ross C11a is at
  C_J = 3.1294 with T = 42.140 d (= 9.6909 TU). **These are different members of
  the same family on different energy manifolds — not the same orbit.**
- No initial-condition state vectors are printed for any family (Table 2 carries
  only period-in-days and σ). The reduced (x,y,θ) model and heading-fan atlas do
  not require, and do not report, x-axis-crossing ICs.
- The labeling differs from our catalogue. Braik-Ross uses **C11a, C11b, C21,
  C32** (four cycler representatives). Our catalogue has **five** rows
  (-11,-21,-31,-32,-33). C11a/C11b are evidently two distinct (1,1) members; there
  is **no (3,3) or split 3-x representative** in Braik-Ross. Mapping is therefore
  *not 1:1* and would have to be reconstructed carefully from the source AAS
  25-621 (k1,k2) definitions before any cross-use.

So nothing here becomes a catalogue row or a golden cross-check on the IC level:
the only quantitative cycler data (period-in-days, σ) are at a non-catalogue
energy and carry no state vector. **Do not wire.**

#### Transcribed candidate data (REPRODUCE-BEFORE-USE — all at C_J = 3.1294, NOT our catalogue energy)

Table 2 (cycler rows only; period in days, σ in TU⁻¹ and day⁻¹):

| Family | Abbrev | T_PO [d] | σ [TU⁻¹] | σ [day⁻¹] |
|---|---|---|---|---|
| (1,1)a-cycler | C11a | 42.140 | 1.0482 | 0.2410 |
| (1,1)b-cycler | C11b | 55.995 | 0.9255 | 0.2128 |
| (2,1)-cycler  | C21  | 84.533 | 0.1358 | 0.0312 |
| (3,2)-cycler  | C32  | 78.613 | 0.6886 | 0.1583 |

These are **method inputs at a single shared energy**, not sourced family
representatives in our catalogue's per-family-stable-C convention. Treat as
context only. No state vector → not independently reproducible as ICs; would
require re-running a fixed-C_J=3.1294 family continuation to recover the member.

### Q2 — The reachable-set accessibility-network method, and its relevance to our family-selection / basin / continuation problem.

**Method, precisely:** at one fixed C_J, reduce to (x,y,θ); from arc-length seeds
on each representative orbit, fan out bounded **energy-preserving heading-change**
maneuvers (|δ| ≤ δ_max, ΔV_turn = 2v sin(|δ|/2)); finite-time-propagate and
voxel-log to get a forward reachable set; mirror via time-reversal symmetry for
the backward set; declare A→B accessible iff F(A)∩B(B)≠∅ on the shared grid;
assign each overlap voxel a proxy (ΔV, T); reduce each family pair to its
min-proxy-ΔV voxel; assemble the symmetric proxy-ΔV matrix into a weighted graph;
read hub/gateway/relay from strength/closeness/betweenness; sweep the
(ΔV_cap, T_cap) budget plane for regime transitions.

**Relevance verdict: METHOD-ONLY, and directly relevant — *positive*.** This is
exactly the family-to-family accessibility-network method the task hoped for, and
it speaks to two of our open lanes:

- **Family selection / which families to continue (#162, #218/#219).** Our
  continuation campaign found 0 distinct new families from sourced seeds; the open
  problem is *which* basins/families are worth seeding from. Braik-Ross gives a
  cheap, energy-fixed *screening layer* that ranks families by how reachable they
  are from the rest of the cislunar set and identifies the dominant
  hub/gateway/relay (C32, C11a) and the hard-access family (R21-S). Implementing a
  reduced (x,y,θ) heading-fan reachable-set overlap at a chosen C_J would let us
  *score candidate seed families before* spending corrector/continuation effort —
  i.e. a prioritizer for the continuation lane rather than a new family generator.
- **Basin-of-attraction framing (#162).** The dynamical reading (unstable
  representatives share one chaotic sea = low-cost mutual access; stable resonant
  tori resist invasion = hard access) is a *qualitative* basin story we can test
  against our own stack: it predicts our stable cyclers should be relatively
  hard-access nodes and our unstable resonant/Lyapunov members easy-access.
- **Connecting the Ross EM cyclers into the wider cislunar family set.** Yes — this
  is precisely what the network does (cyclers + Lyapunov + DPO occupy the low-cost
  block; resonants the expensive block). It would let us situate our 5 catalogue
  cycler rows relative to Lyapunov/resonant/DPO families on a common manifold.

**Caveats / cost of adoption:** (a) it is *planar only*, single common C_J,
representative-orbit (not full continued family); (b) the maneuver model is a
strong restriction (pure heading rotation at fixed rotating-frame speed — *not*
general impulsive ΔV; authors note this is often a large normal-component burn,
not the tangential burn 2-body intuition expects); (c) proxy costs are *screening*
quantities, not optimized transfers (validated only as a conservative upper bound,
corrected < proxy); (d) builds a 3-D voxel grid (Δx=Δy=0.001 LU≈384 km, Δθ=1°)
over the admissible domain — implementable but non-trivial infra. It is a
*screening/initialization layer*, by the authors' own framing, not a transfer
designer.

### Q3 — Cross-checkable numbers against our CR3BP stack / JPL oracle.

Several, all at **C_J = 3.1294** (a *new* energy for us; our catalogue cyclers are
at their per-family stable C, so this is a fresh independent verification target,
not a same-orbit check):

- **μ = 0.012150584270572** (Table 1). Matches our Earth-Moon μ
  (0.01215058439469525) to 9 sig figs; the 10th-digit difference is a DE source/
  rounding difference — worth noting but consistent. (Our catalogue rows use
  1.2150584270572e-2 — *identical* to Braik-Ross Table 1.)
- **LU = 384400 km, T_EM = 27.321661 d, TU = 4.34837740 d, R_E = 0.0165921,
  R_M = 0.00451873** (Table 1) — all standard; R_E/R_M check against
  6378/384400 and 1737/384400.
- **Equations of motion / Jacobi convention** (Eqs. 1–5): C_J = −2Ū − (ẋ²+ẏ²)
  with Ū = −½(x²+y²) − (1−μ)/r1 − μ/r2. **Sign convention check:** this is the
  *negative-Ū* form, i.e. C_J = (x²+y²) + 2(1−μ)/r1 + 2μ/r2 − v², which **matches
  our catalogue jacobi_constant convention** (already noted in row -11). Good
  consistency anchor.
- **Periods (days) at C_J=3.1294, Table 2** — could be cross-checked by running our
  fixed-C_J continuation/corrector at C_J=3.1294 for, e.g., LL1 (12.811 d), LL2
  (15.117 d), DPO (11.184 d), R21-S (26.500 d). These are *clean independent
  targets at a single sourced energy with sourced periods* — a candidate golden
  set IF we reproduce the member at C_J=3.1294 first. The instability rates σ
  (Floquet exponents, Table 2) are also checkable against our Barden/monodromy
  routine at that energy.
- **JPL oracle relevance:** Table 2 representatives (Lyapunov, resonant, DPO) are
  standard JPL 3-Body Periodic-Orbit families; the JPL oracle
  (`src/cyclerfinder/verify/jpl_periodic_orbits.py`, #116) can pull L1/L2 Lyapunov
  and resonant ICs at C_J=3.1294 to cross-check Braik-Ross periods. Cyclers are
  NOT in the JPL DB, so only the Lyapunov/resonant/DPO subset is JPL-checkable.

No state vectors are published, so any reproduction must first *recover the member*
at C_J=3.1294 (corrector to the sourced period), then compare. The check would
verify period+σ, not ICs.

---

## PROPOSED FOLLOW-ON

This paper carries **no new sourced family tuples** (Q1 = N), so the follow-on is
**method-informed**, not a row ingest. Two proposed tasks, in priority order:

1. **[Method, family-selection lane — recommended] Prototype a reduced (x,y,θ)
   heading-fan reachable-set overlap *scorer* as a continuation prioritizer.**
   Scope: at a chosen C_J, build forward/backward reachable sets (heading fan +
   time-reversal mirror) for our existing sourced representatives + a small set of
   candidate families, and compute the pairwise overlap proxy-ΔV matrix. Use the
   resulting strength/closeness/betweenness ranking to *pre-screen which families
   to seed* in the next continuation pass (#218/#219), instead of seeding blindly.
   Deliverable is a screening tool + a ranked candidate list, explicitly framed as
   prioritization, NOT a new-family generator. This directly attacks the
   family-selection half of #162 with a published, validated method. Note the
   planar/single-C_J/heading-only limitations up front. Reproduce-before-trust:
   first replicate Braik-Ross's C32-dominant / R21-S-hard-access result at
   C_J=3.1294 as a method-validation gate before trusting it on our families.

2. **[Cross-check, low cost] Independent period + σ verification at C_J=3.1294.**
   Run our fixed-C_J corrector/continuation at C_J=3.1294 for the JPL-available
   subset (LL1, LL2, DPO, R21-S, R31-S/U, R52-S/U) and the four cyclers, recover
   each member to the sourced period (Table 2), and compare period + Floquet σ.
   Pull the Lyapunov/resonant/DPO ICs from the JPL oracle (#116) as the recovery
   seed. This is a *new-energy* cross-check (our catalogue cyclers are at their own
   stable C, so no same-orbit overlap) — a fresh independent anchor on our stack,
   and a prerequisite gate for task 1. Cyclers have no JPL entry, so recover them
   from our own corrector at C_J=3.1294.

**Not proposed:** ingesting C11a/C11b/C21/C32 Table-2 data as catalogue rows —
they are at a non-catalogue energy, carry no state vector, and the labeling does
not map 1:1 to our -11/-21/-31/-32/-33 rows. Adopting them would create
energy-mismatched duplicate cycler entries. Hold.

---

## One-line summary

Method-only paper (reachable-set family-accessibility network); **no new sourced
cycler tuples** — uses the AAS 25-621 cyclers as named representatives re-selected
at a single common C_J=3.1294 with periods/σ only, no ICs, labeling not 1:1 to our
catalogue. Directly relevant to the family-selection/continuation lane (#162,
#218/#219) as a screening prioritizer; periods+σ at C_J=3.1294 are a fresh
independent cross-check target (reproduce-the-member-first).
