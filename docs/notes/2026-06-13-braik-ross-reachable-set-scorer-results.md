# Braik-Ross reachable-set accessibility-network scorer — build + validation-gate results (#236)

**Date:** 2026-06-13
**Method source:** Abdullah Braik & Shane D. Ross, *Orbital Networks in the
Three-Body Problem*, arXiv:2605.31543 (2026).
**Spec:** PROPOSED FOLLOW-ON task 1 of
`docs/notes/2026-06-13-braik-ross-2026-orbital-networks-mining.md`.
**Framing:** a family-selection / continuation *prioritizer* (a screening tool),
NOT a transfer designer and NOT a new-family generator.

---

## TL;DR

- The scorer is **built and unit-tested** (reduced (x,y,θ) model, Eq. 26 turn
  cost, Eq. 14 time-reversal mirror, voxel overlap → proxy ΔV/time, N×N matrix,
  cost-aware strength / harmonic-closeness / betweenness).
- **Member recovery at C_J = 3.1294 independently CONFIRMS Braik-Ross Table 2.**
  LL1, LL2, DPO, R21-S and the C11b cycler recover to the sourced periods to
  ≤ 0.005 d (essentially exact); C32 recovers to ~1.1%.
- **Validation gate: PARTIAL / faithful negative.**
  - REPRODUCED: **R21-S is a hard-access (bottom) node** (last/second-last in
    strength and harmonic closeness) — the "stable resonant tori resist
    invasion" reading holds.
  - NOT REPRODUCED: **C32 is NOT the dominant family** on the recoverable subset
    (it ranks 5th of 6 in strength/closeness; betweenness is zero for all nodes
    on the small complete graph).
- **Root cause of the negative:** C32's hub/gateway/relay dominance is a property
  of the *full 13-node chaotic-sea network*. Three of the four cyclers (C11a,
  C21) and the R21-U / R31-S/U / R52-S/U resonants do **not** recover via the
  available 1-DOF perpendicular-x-crossing symmetric corrector at this off-stable
  common energy, so the scored network is only **6 nodes**. C32's dominance does
  not survive that truncation.
- **Per the discipline, parameters were NOT tuned to force the published answer,
  and the scorer is left GATED**: it was **not** applied to rank our families,
  because the headline gate (C32 dominant) did not pass.

---

## What was built

`src/cyclerfinder/search/reachable_network.py` — the scorer:
- `reduced_speed` (Eq. 8 `v = √(−2Ū − C_J)`), `ubar`, `is_admissible` (Hill
  region), `dv_turn` (Eq. 26 `2v·sin(|δ|/2)`), `time_reversal` (Eq. 14
  `R(x,y,θ)=(x,−y,π−θ)`), heading/angle helpers.
- `VoxelGrid` — uniform (x,y,θ) grid (modular θ axis), index ↔ center.
- `ReachableSet` + `build_reachable_set` — arc-length-spaced seeds, heading fan
  in `[−δ_max, δ_max]`, one-sided propagation to a common horizon T_a, voxel-log
  of min source-side turn cost + crossing time.
- `mirror_reachable_set` — backward set = time-reversal mirror (free).
- `pair_proxy` / `proxy_matrix` — overlap → min-proxy-ΔV (source turn + target
  turn + half-voxel heading-mismatch patch), symmetric N×N matrix.
- `centralities` — strength (Σ 1/w), harmonic closeness (Σ 1/shortest-path),
  Brandes betweenness; `_floyd_warshall` all-pairs.

`src/cyclerfinder/search/reachable_representatives.py` — member recovery at
C_J = 3.1294: JPL-oracle seeds (LL1/LL2/DPO/R21-S) re-corrected under our μ, plus
cycler recovery from the Ross-RT seed regions, each confirmed against the sourced
Table-2 period.

Tests:
- `tests/search/test_reachable_network.py` — 19 hand-checkable mechanics unit
  tests (Eq. 8 / 14 / 26, voxel indexing, overlap/proxy, centralities on tiny
  star/path graphs).
- `tests/search/test_reachable_network_gate.py` — the C_J=3.1294 gate (marked
  `slow`): period-confirmation, R21-S hard-access (passes), C32-dominant
  (`xfail`, documented).

---

## Member recovery at C_J = 3.1294 (independent cross-check of Table 2)

No state vectors are published; we recovered each member and confirmed against the
sourced period (days). Recovered to:

| Family | Sourced T [d] | Recovered T [d] | Source of seed | Status |
|---|---|---|---|---|
| LL1   | 12.811 | 12.811 | JPL L1 Lyapunov | exact |
| LL2   | 15.117 | 15.117 | JPL L2 Lyapunov | exact |
| DPO   | 11.184 | 11.184 | JPL DPO         | exact |
| R21-S | 26.500 | 26.500 | JPL resonant br=21 (stable) | exact |
| C11b  | 55.995 | 55.995 | Ross-RT (1,1) seed | exact |
| C32   | 78.613 | ~79.50 | Ross-RT (3,2) seed | ~1.1% (off-stable branch) |
| C11a  | 42.140 | — | Ross-RT (1,1) seed | **not recovered** |
| C21   | 84.533 | — | Ross-RT (2,1) seed | **not recovered** |

The LL1/LL2/DPO/R21-S exact matches are a clean, *independent* confirmation of
Braik-Ross Table 2 at a fresh energy (our catalogue cyclers live at their own
per-family stable C, so this is not a same-orbit check) — proposed task 2 of the
mining note, delivered as a side effect.

C11a, C21 (and C32 to ~1.1%) are a genuine **recovery gap**: the 1-DOF
perpendicular-x-axis-crossing symmetric corrector that recovers the *stable*
members at their own energy does not cleanly recover these cyclers at the
off-stable common energy C_J=3.1294. A multi-segment / 2-D shooting corrector
(the kind already flagged for the multi-arc S1L1 problem) would be needed. They
were **excluded** from the scored network rather than faked with an
unconfirmed member.

---

## Validation-gate result (the reproduce-before-trust gate)

Scored network = the 6 source-confirmable nodes
{LL1, LL2, DPO, R21-S, C11b, C32}, common accessibility horizon
T_a = 2π TU ≈ 27.32 d (paper reference T_cap), grid Δx=Δy=0.02 LU, Δθ=10°,
n_seeds=10, n_fan=9, δ_max=30°.

Cost-aware centralities (high → low):

```
strength   : C11b=38.16, DPO=36.84, LL1=35.95, LL2=32.20, C32=24.02, R21-S=23.86
harmonic   : C11b=38.16, DPO=36.84, LL1=35.95, LL2=32.20, C32=24.02, R21-S=23.86
betweenness: all zero (complete 6-node graph — every shortest path is a direct edge)
```

- **R21-S hard-access: REPRODUCED.** R21-S is last in both strength and harmonic
  closeness, consistent with the published "stable resonant tori resist
  low-energy invasion" reading. (`test_r21s_is_hard_access_on_subset` passes.)
- **C32 dominant: NOT REPRODUCED.** C32 ranks 5th of 6 (just above R21-S), not
  rank 1. (`test_validation_gate_c32_dominant` is recorded as `xfail`.)
- Betweenness carries no signal here: with 6 nodes the proxy graph is complete,
  so no node sits on another pair's shortest path. The paper's nonzero
  betweenness arises at *budget caps* on the *13-node* network where some direct
  edges (notably the 3 R21-S edges) drop out and force relay routing — neither
  condition is present on this small, fully-connected subset.

### Why the negative is the subset, not the method

C32's headline dominance is intrinsically a full-network statement: among the
*ten unstable* representatives sharing one connected chaotic region, C32 is the
best-connected. On a 6-node subset that omits C11a, C21 and the four extra
resonants, that connectivity advantage is gone — the remaining low-cost block is
dominated by the compact, dense-coverage orbits (C11b, DPO, the Lyapunovs), and
the long-period unstable C32 (whose heading-fanned arcs diverge and partly escape
the grid over the horizon) scores low. The result is stable under switching from
per-orbit-period to a common horizon, so it is not a horizon artifact.

This is a **faithful negative on the headline claim** with the **secondary claim
reproduced**. No parameters were tuned to recover C32's rank.

---

## Consequence: the scorer stays gated

Because the headline gate (C32 dominant) did not reproduce, the scorer was **not**
applied to rank our sourced families. Doing so would be unjustified — the gate
exists precisely to license that step. To ungate it, the missing nodes must first
be recovered:

1. Recover C11a, C21 (and a tightened C32) at C_J=3.1294 with a multi-segment /
   2-D shooting corrector (the perpendicular-x-crossing 1-DOF corrector is
   insufficient for these cyclers at the off-stable energy).
2. Recover R21-U, R31-S/U, R52-S/U. R21/R31 are in the JPL DB (branches `21`,
   `31`) with stable+unstable members at C_J=3.1294; **5:2 (R52) is not exposed
   by the JPL periodic_orbits API** under the tried branch codes and would need
   another source.
3. Re-run the gate on the full 13-node network. Only if C32 then ranks 1 across
   strength/closeness/betweenness should the scorer be used to prioritise our
   continuation seeds.

---

## Limitations (carried up front, per the paper's own framing)

- **Planar only**, single common **C_J = 3.1294** (a non-catalogue energy).
- **Heading-only maneuver**: a pure rotation of the rotating-frame velocity at
  fixed speed (Eq. 26) — *not* general impulsive ΔV (often a large
  normal-component burn, not the tangential burn 2-body intuition expects).
- **Screening, not transfer**: the proxy ΔV is necessary-not-sufficient and a
  conservative upper bound (the paper's corrected ΔV < proxy in every tested
  case). It pre-screens which families to seed; it does not design transfers.
- **No catalogue writeback**: Braik-Ross carries no new sourced tuples; the
  labeling is not 1:1 with our rows and the energy differs. This is a tool, not
  data. Nothing was written to `data/catalogue.yaml`.

---

## How to use it (once ungated)

```python
from cyclerfinder.search import reachable_network as rn, reachable_representatives as rr

sysm = rr.braik_ross_system()
reps = [...]  # recovered Representative objects (state0, period) at a chosen C_J
grid = rn.VoxelGrid(dx=0.02, dy=0.02, dtheta=math.radians(10.0))
fwd  = [rn.build_reachable_set(sysm, r.state0, r.period, grid, C_J, horizon=2*math.pi) for r in reps]
back = [rn.mirror_reachable_set(f, grid) for f in fwd]
cent = rn.centralities(rn.proxy_matrix(fwd, back, grid))
# rank families by cent.strength / cent.harmonic_closeness / cent.betweenness
```
