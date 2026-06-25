# Multi-rev resonant-leveraging releg — the in-band moon-tour capability (DESIGN DRAFT, #465)

**Date:** 2026-06-26
**Status:** DESIGN-DRAFT — for user review. No production code written by this doc.
**Issue:** #465. Direct follow-on to **#449** (single-DSM releg) and **#464**
(single low-thrust releg), which both **closed the Galilean moon-tour
geometrically but OUT-OF-BAND** — 13.18 km/s/cycle (DSM) and 12.03 km/s/cycle
(SF), vs the 3.5 km/s/cycle powered ceiling.

Sibling note: #466 (energy-moving QP walk) edits
`src/cyclerfinder/genome/qp_tori_arclength.py` — a different subsystem. This draft
touches only the `search/` releg area.

---

## 0. TL;DR

- **The #449/#464 gate, precisely.** Both prior relegs replace ONE ballistic
  Lambert leg with ONE powered retarget arc that sheds the leg's whole V∞ defect
  in a single impulse (DSM) or one distributed train (SF). The #464 verdict
  (`2026-06-26-464-lowthrust-releg-verdict.md` §4) is explicit: the 12-13 km/s
  cost is **the irreducible magnitude of the V∞ retarget itself**, not
  impulse-vs-distributed inefficiency. At a common flyby target of V∞≈4 km/s the
  individual Galilean legs naturally want arrival V∞ of 5-9 km/s, so each single
  leg must shed 2-5 km/s of V∞ *at once*. A single retarget pays that shed at the
  expensive *single-VILM maximum* (`vilm.py` Eq.14, `ΔV_max`), not the cheap
  *multi-VILM minimum* (Eq.13 quadrature).
- **What the literature does instead (the cheap closure).** Campagnola-Russell
  "Endgame Problem" Part 1 (digested `2026-06-05-endgame-tisserand-mining.md`,
  golden in `search/vilm.py`) walks V∞ down across MANY revolutions with a CHAIN
  of resonant gravity-assist legs (a VILM endgame), each shedding a few tens to a
  few hundred m/s. The total is the Eq.(13) quadrature — the
  **theoretical-minimum** ΔV — which is roughly **an order of magnitude cheaper**
  than the single-impulse shed.
- **#465 = the chained releg.** A new `MultiRevLeveragingReleg` backend behind the
  EXISTING `Releg` protocol (`search/releg_solver.py`). Where `DsmReleg` does one
  retarget impulse, this backend internally **chains N resonant-hop legs**
  (`search/leveraging_leg.py::evaluate_leveraging_leg`, the #179 phase-full apse
  burn) to walk V∞ from the leg's natural arrival value down to the common flyby
  target, summing the per-hop ΔV. It returns total delivered ΔV in the SAME
  `RelegResult` so `releg_moontour.close_powered_cycle` scores it against the
  3.5 km/s band exactly like the other backends — **no driver rewrite**.
- **In-band evidence (computed from the sourced quadrature, this draft §6).** For
  a flyby cycler the escape+capture insertions are NOT paid (it never inserts into
  orbit at a moon — it flies by). The cycler's per-transfer cost is the
  **leveraging-only** (begingame + endgame) quadrature. Computed via `vilm.py`:
  Ganymede↔Europa **0.31 km/s**, Europa↔Io **0.41 km/s**, Ganymede↔Io
  **0.74 km/s**, Callisto↔Ganymede **0.27 km/s**. A full Io-Europa-Ganymede-Io
  cycle leverages at **≈1.5 km/s/cycle — comfortably IN-BAND** (vs the 13.18 km/s
  single-leg DSM, vs the 3.5 km/s ceiling). This is the ~25× gap the prior relegs
  could not capture because they shed V∞ in one impulse.
- **HONEST confidence (this is the gate — see §8).** **HIGH** that multi-rev
  leveraging brings the **Galilean and Saturnian** tours in-band: the quadrature
  numbers above are sourced and well under the ceiling. **STRUCTURALLY NEGATIVE
  and unchanged** for **Uranus/Neptune**: the VILM-floor prefilter still forbids
  the bridge (disjoint Tisserand contours at every V∞), and chaining cannot bridge
  contours that don't intersect — leveraging walks V∞ *within* a contour, it can't
  jump between disjoint ones. So #465 re-opens Jovian/Saturnian in-band and leaves
  Uranus/Neptune a (now stronger) powered-empty.
- **The honest residual risk:** the quadrature is the *infinite-VILM continuous
  floor*; a finite chain of integer-resonance hops costs MORE (the published
  Europa endgame is 154 m/s discrete vs the 128 m/s continuous floor — a +20%
  finite-chain penalty). The Galilean cycle has comfortable margin (1.5 km/s vs
  3.5), so the finite penalty does not breach the band. The risk is real only if
  the integer-resonance hop sequence cannot reach the required V∞ in a feasible
  number of revs (phasing infeasibility), which is the kill-criterion (§7).

---

## 1. The mechanism — precisely how the chain replaces the single retarget

### 1.1 What the single-leg releg does (and why it is expensive)

`DsmReleg`/`LowThrustReleg` (today): given a leg `M_k → M_{k+1}`, the driver pins
departure V∞ to a common target `T` and asks the backend to **retarget the arrival
V∞ to `T` in one arc**. The arc spends `|ΔV|` to change the V∞ magnitude by
`|V∞_arrival_natural − T|` (2-5 km/s in the Galilean basin). That single change is
the *single-VILM maximum* cost: `vilm.py` Eq.(14), `ΔV_max = −Γ + sqrt(Γ² +
(V∞H² − V∞L²))`. One big impulse.

### 1.2 What the multi-rev chain does (and why it is cheap)

Replace the single arc with a SEQUENCE of resonant-hop legs that each shed a small
slice of V∞ over one resonance (one to a few revolutions), so the V∞ walks down
the leveraging graph step by step:

```
  arrival V∞ = V∞_H  (natural, high)
        │  hop 1: n1:m1 resonant leg at the moon, apse burn δ1  → V∞ drops a little
        ▼
      V∞_1
        │  hop 2: n2:m2 resonant leg, apse burn δ2             → V∞ drops a little
        ▼
      V∞_2
        ⋮   (K hops, each a phase-full VILM at the shared flyby moon)
        ▼
  departure V∞ = V∞_L = T  (common flyby target, low)
```

Each hop is exactly one `evaluate_leveraging_leg(moon, n, m, vinf_in, vinf_out_target,
exterior)` call (#179, `search/leveraging_leg.py`): a resonant orbit `a=(n/m)^(2/3)`,
a tangential apse burn solving the Tisserand quadratic for the target V∞_out, the
burn magnitude `δ = |v' − v|`, with the analytic Γ-floor cross-check
(`gamma_floor_kms`). The total chain ΔV is `Σ δ_i`, which by the Eq.(13)
construction approaches the quadrature minimum `∫_{V∞L}^{V∞H} V∞/Γ dV∞` — the
**theoretical minimum**, ~an order of magnitude under the single-impulse `ΔV_max`.

**This is the V∞-leveraging endgame the prior relegs could not express:** the
single-leg `Releg` contract has no notion of intermediate revolutions, so it
could only pay `ΔV_max`. The chain spends `ΔV_min`.

### 1.3 The cycler subtlety — no escape/capture (the in-band lever)

Campagnola-Russell Table 1 (A2) tabulates the FULL moon-to-moon ΔV_min =
escape + begingame + endgame + capture, because their use case INSERTS into orbit
at each end (e.g. Europa orbiter). A **cycler flies by** every moon — it never
escapes a parking orbit or captures into one. The cycler pays ONLY the leveraging
quadratures (begingame at the outer moon + endgame at the inner moon). That is the
`_leverage_dv_kms` half of `vilm_dv_min`, and it is what §6 evaluates: 0.27-0.74
km/s per transfer, vs the 1.7-2.4 km/s full Table-1 value. Dropping escape+capture
is the single biggest reason the cycler closes in-band where an orbiter mission
would not.

---

## 2. The `Releg`-seam plug-in (reuse, no new protocol)

The backend is a NEW `Releg` implementation behind the SAME protocol
(`search/releg_solver.py::Releg`), so `releg_moontour.close_powered_cycle`
consumes it unchanged (the driver already selects "powered iff not
`BallisticReleg`", `releg_moontour.py:263`, and pins/retargets via
`vinf_depart_mag`/`vinf_target_in`).

```
class MultiRevLeveragingReleg:                 # NEW backend, search/releg_solver.py
    def solve(r_a, v_a, r_b, v_b, tof_s, mu, *, n_rev=0,
              vinf_target_in=None, vinf_depart_mag=None) -> RelegResult:
        # 1. Seed: the ballistic Lambert leg fixes the geometry + the NATURAL
        #    arrival V∞ magnitude V∞_H (lowest-energy branch — same seed the other
        #    backends use, releg_solver.py:165-174).
        # 2. If vinf_target_in is None: ballistic-equivalent (no walk needed).
        # 3. Else: CHAIN of leveraging hops at the arrival moon B from V∞_H down to
        #    vinf_target_in, each evaluate_leveraging_leg(...) (#179); sum δ_i.
        #    Pick the hop resonances (n:m) by a small descent that keeps each hop
        #    feasible (real near-root, bound post-burn orbit, moon-crossing) — the
        #    leveraging-graph "zigzag low on the x-axis" (mining note p.11).
        # 4. dv_kms = Σ δ_i ; vinf_in = vinf_target_in ; vinf_out = vinf_depart_mag.
        return RelegResult(vinf_out, vinf_in, dv_kms=sum_delta, feasible=...)
```

Key contract points (preserve the driver's by-construction continuity):

- `vinf_depart_mag` PINS departure V∞ to `T`; `vinf_target_in` RETARGETS arrival to
  `T`. The chain walks the ARRIVAL side from `V∞_H` (natural) down to `T`. (The
  departure-side pin is by construction — the cycle is continuous because every
  leg departs and arrives at the same `T`, the driver's existing model,
  `releg_moontour.py:24-32`.)
- The chain's ToF is the sum of the hop periods (`Σ n_i·T_M`), which the driver can
  carry but does NOT gate on for the capability proof (the published endgame ToFs
  are 46-291 days; the cycler's per-cycle period is a downstream V2 concern, not a
  closure gate). The released ΔV is what the band gates.
- Zero-retarget limit: `vinf_target_in is None` (or == natural V∞_H) ⇒ no hops ⇒
  `dv_kms = 0` ⇒ reproduces `BallisticReleg` (the regression limit every backend
  must honour, `feedback_orbit_closure_discipline`).

**No new optimiser, no new cost model.** The hop primitive (`leveraging_leg.py`,
#179) and the floor (`vilm.py`, golden) both exist. The new work is the **descent
orchestration**: choosing the hop resonances and summing δ. The Γ-floor
(`gamma_floor_kms`) is the per-hop sanity lower bound (a hop ΔV below its Γ floor
is non-physical — already enforced in `leveraging_leg.py:202`).

---

## 3. Which dead regions this re-opens IN-BAND

| Region (registry id) | #449/#464 single-leg | #465 multi-rev chain (this draft) | In-band? |
|---|---|---|---|
| Jovian Galilean (`jovian-IEG-vilm-2026-06-09`, `jovian-perm-vilm-2026-06-09`, `jupiter-galilean-amalthea-…`) | closed at 13.18 / 12.03 km/s/cycle — OUT | leveraging-only ≈1.5 km/s/cycle (§6) | **YES (re-open)** |
| Saturnian (`saturnian-titan-vilm-2026-06-09`, `saturnian-titan-endgame-vilm-2026-06-10`) | high-V∞ basin, gap ~7-13 km/s | Titan-Rhea-Dione lev-only ≈0.5-0.7 km/s/cycle (§6) | **YES (re-open)** |
| `repeated-moon-{jupiter,saturn}-sweep` | empty (ballistic reach) | powered chain re-test, in-band where contours link | **YES (re-test)** |
| Uranian / Neptunian (`uranus-neptune-regular-moon-endgame-vilm-2026-06-23`, `repeated-moon-{uranus,neptune}-sweep`) | disjoint Tisserand contours (STRUCTURAL) | **STILL EMPTY** — chaining walks V∞ *within* a contour, cannot jump disjoint contours; prefilter (`moon_prune.moon_leg_admissible`) skips before any chain solve | **NO (stronger powered-empty)** |
| `repeated-moon-mars-sweep` | empty | Mars has one moon-pair of interest; leveraging needs ≥1 shared-primary resonance — likely stays empty | (re-test, expect empty) |

The re-open is **honest and selective**: Jovian + Saturnian come in-band; Uranus +
Neptune are re-stamped as a *stronger* powered-empty (a powered chain that ALSO
can't bridge is a more final negative than the ballistic one), via the EXISTING
`build_powered_empty_restamp` (`releg_moontour.py:354`) with the new backend's
capability tags.

---

## 4. Data flow

1. **Input:** a tour skeleton `(primary, sequence, per-leg tof, n_rev, phasing)` —
   the SAME existing `repeated-moon-*-sweep` / VILM-negative skeletons #449 relegs.
   #465 does not enumerate skeletons; it re-legs given ones with a cheaper backend.
2. **VILM/linkability prefilter** (`moon_prune.moon_leg_admissible`, unchanged):
   skip a skeleton whose any leg is unbridgeable at every probed V∞ (the
   Uranus/Neptune disjoint case). This is the structural negative, kept exactly.
3. **Per-leg multi-rev chain solve:** for surviving legs, `MultiRevLeveragingReleg`
   walks the arrival V∞ from natural down to the common `T` via N resonant hops,
   returning `Σ δ` in `RelegResult.dv_kms`.
4. **Cycle close + wrap:** `close_powered_cycle` enforces post-retarget
   V∞-continuity (by construction at `T`) + the closed-cycle wrap, sums per-cycle
   ΔV — unchanged driver logic.
5. **Powered dv-band gate** (`verify/dv_band_acceptance.classify_dv_band`): the
   per-cycle ΔV must land in `[300 m/s, 3.5 km/s]`.
6. **Novelty pipeline** (unchanged): signature dedup → `literature_check` (V0,
   necessary-not-sufficient, mandatory `feedback_literature_novelty_check_baseline`)
   → V2 moontour (releg-aware) → SILVER holding.
7. **Output:** classified ledger; an in-band Jovian/Saturnian close is a
   `quasi_cycler`/`mga_tour` SILVER candidate (V0-known if it matches a published
   tour, e.g. Callisto-Ganymede-Callisto; SILVER-novel otherwise). A Uranus/Neptune
   empty re-stamps the registry with the subsuming method+version.

---

## 5. Reuse map (what already exists vs the new ~1 file)

| Need | Already built | New work |
|---|---|---|
| Per-hop resonant leverage leg (apse burn, Tisserand quadratic, Γ-floor check) | `search/leveraging_leg.py::evaluate_leveraging_leg` (#179, validated) | — |
| Quadrature ΔV floor / golden | `search/vilm.py::vilm_dv_min`, `_leverage_dv_kms`, `europa_endgame_dv` (golden vs Tables 1/2 + 154/147 m/s) | — |
| Releg protocol + driver + prefilter + dv-band gate + empty re-stamp | `search/releg_solver.py`, `search/releg_moontour.py`, `search/moon_prune.py`, `verify/dv_band_acceptance.py` | — |
| **Chain orchestration** (choose hop resonances, walk V∞_H→T, sum δ) | — | **`MultiRevLeveragingReleg.solve` in `releg_solver.py` (one new class)** |
| Capability tags for the subsumption record | `data/method_capability` edges (`leveraging` ⊐ `single-arc`) | one new tag set `multi-rev-leveraging` ⊐ `one-dsm-per-leg`, `leveraging` |

The single new object is the chain orchestrator. Everything else is import-and-call.

---

## 6. Validation / golden (sourced, non-circular)

**Primary golden — the leveraging quadrature (already golden-validated in
`vilm.py`).** The chain's summed ΔV for a moon-A→moon-B leveraging walk must:

- be **≥** the Eq.(13) continuous quadrature floor `vilm._leverage_dv_kms`
  (begingame + endgame) — a finite chain of integer-resonance hops cannot beat the
  infinite-VILM continuous minimum (`feedback_golden_tests_sourced_only`: the
  EXPECTED lower bound traces to the paper's Eq.13, never a self-computed number);
- be within the **+20% finite-chain band** the published Europa endgame sets: the
  3-VILM discrete design is **154 m/s** vs the **128 m/s** continuous floor
  (`vilm.europa_endgame_dv`, sourced A6, mining note 436-438). So the golden is
  `floor ≤ chain ≤ 1.2·floor + slack`, with the EXPECTED edges both sourced
  (floor = Eq.13; ceiling = the published discrete-VILM penalty).

**The sourced in-band number (the gate evidence).** Computed from `vilm.py` for the
flyby-cycler leveraging-only cost (NO escape/capture):

| Inter-moon transfer | begingame | endgame | **lev-only ΔV (km/s)** | full Table-1 (orbiter) |
|---|---|---|---|---|
| Ganymede↔Europa | 0.122 | 0.183 | **0.305** | 1.715 |
| Europa↔Io | 0.180 | 0.229 | **0.409** | 1.760 |
| Ganymede↔Io | 0.258 | 0.477 | **0.735** | 2.305 |
| Callisto↔Ganymede | 0.111 | 0.154 | **0.265** | 1.793 |
| Titan↔Rhea | — | — | **0.334** | 1.149 |
| Rhea↔Dione | — | — | **0.196** | 0.519 |

A Io-Europa-Ganymede-Io cycle leverages at ≈ 0.41 + 0.31 + 0.74 ≈ **1.45
km/s/cycle**. A Saturnian Titan-Rhea-Dione-Titan cycle at ≈ 0.33 + 0.20 + (Titan-Dione)
≈ **0.8-1.0 km/s/cycle**. Both **inside** the 3.5 km/s/cycle ceiling — the in-band
closure the single-leg relegs (13.18 / 12.03 km/s/cycle) missed.

**Capability golden (positive control).** Run `MultiRevLeveragingReleg` on the
Galilean positive control through `close_powered_cycle`; assert the cycle closes
(continuity below the 0.05 km/s gate, by construction) AND total ΔV is inside the
powered band — proving the chain does what the single leg could not.

**Structural-emptiness golden (the honesty test, unchanged).** The driver on a
Uranian Ariel→Umbriel leg must report unbridgeable via the prefilter (disjoint
contours), reproducing `uranus-neptune-regular-moon-endgame-vilm-2026-06-23` as a
powered re-test that STILL finds it empty. The chain must NOT fabricate a bridge.

**Why golden-clean.** Every EXPECTED value traces to Campagnola-Russell's printed
tables / the digested mining note (the Eq.13 quadrature floor and the 154/147 m/s
discrete endgame), never to a number this backend computed.

---

## 7. Risks + kill-criteria

| Risk | Mitigation / kill-criterion |
|---|---|
| **Finite integer-resonance chain cannot reach the target V∞** (no feasible n:m hop sequence walks V∞_H→T in a sane number of revs — phasing infeasibility). | Each hop is gated by `evaluate_leveraging_leg` feasibility (real near-root, bound + moon-crossing post-burn orbit). If the descent stalls before reaching `T`, the leg is INFEASIBLE (return `feasible=False`) — never a fabricated close. KILL the in-band claim for a region if no chain reaches `T` under the band; re-stamp powered-empty. |
| **The chain's ΔV beats its own quadrature floor** (a bug — a finite chain CANNOT beat the continuous minimum). | Golden asserts `chain ≥ floor − tol`. A sub-floor result is a correctness failure, not a discovery. |
| **ToF blows up** (a deep V∞ walk needs hundreds of revs / many years). | ToF is REPORTED and carried to V2 (the published endgames are 46-291 days; a multi-year cycler period is a V2/drift concern, not a closure gate). A pathologically long chain (e.g. >N_max revs) is a feasibility kill for that leg. |
| **"In-band" is read as "novel catalogue cycler".** | Hard boundary: the lit-novelty gate (`literature_check`, mandatory) + the V0-known vs SILVER-novel split stay. Callisto-Ganymede-Callisto is a V0-known admission; only a literature-fresh in-band tour is a SILVER candidate. |
| **Uranus/Neptune "rescued" by chaining** (it cannot be). | The prefilter (`moon_prune`) skips disjoint-contour legs BEFORE any chain solve; the structural negative is preserved and re-stamped, not overturned. Chaining moves V∞ within a contour only. |
| **Continuous-floor optimism** (Eq.13 is linked-conic; CR3BP differs ±5-10%, mining note p.6). | The Galilean cycle has ~2.3× band margin (1.5 vs 3.5 km/s); a +10% CR3BP correction does not breach it. Saturnian has even more. Flag the model band on any resulting row (the existing CR3BP-vs-patched-conic catalogue flag). |

---

## 8. HONEST confidence assessment (the gate)

**Does multi-rev leveraging actually bring the Galilean tour in-band?**

**YES — HIGH confidence — for Jovian and Saturnian.** The evidence is the sourced
Eq.(13) quadrature, computed in §6: the flyby-cycler leveraging-only cost is
0.27-0.74 km/s per inter-moon transfer, ≈1.5 km/s for a full Galilean cycle —
roughly **2.3× under** the 3.5 km/s ceiling and **~9× under** the 13.18 km/s
single-DSM closure. The gap between #449/#464 and #465 is exactly the gap between
the single-VILM *maximum* (Eq.14, one impulse) and the multi-VILM *minimum*
(Eq.13, the chain). The escape/capture drop (§1.3) adds further margin. The
finite-chain +20% penalty (154 vs 128 m/s, sourced) does not breach the band given
the 2.3× margin.

**NO — and honestly so — for Uranus/Neptune.** This is NOT brought in-band by
chaining and the draft does not pretend otherwise. The block is *disjoint Tisserand
contours*: there is no V∞ at which consecutive Uranian/Neptunian moons share a
resonance contour, and leveraging only walks V∞ *along* a contour — it cannot jump
between disjoint ones. The prefilter forbids it before any chain runs. #465
re-stamps these as a stronger powered-empty, not a rescue.

**The one genuine open risk** (not band magnitude, but reachability): whether a
*finite* sequence of *integer-resonance* hops can actually walk V∞ from the natural
high value to the common target in a feasible number of revolutions for each
specific Galilean leg. The continuous quadrature says the *cost* is in-band; the
finite chain must realise it with real n:m hops. The Europa endgame (sourced 3-VILM
/ 14-VILM designs, A6) is the existence proof that such finite chains exist for at
least one Galilean endgame — which is strong but not a blanket guarantee for every
leg. Confidence the chain *realises* the in-band cost on the Galilean positive
control: **MEDIUM-HIGH** (the published Europa endgame realises it; the others are
expected to by the same construction, to be confirmed by the positive-control
test).

**Net:** #465 is the genuine path to an in-band Jovian/Saturnian moon-tour cycler.
The capability-and-cost case is sourced and strong; the remaining empirical
question is finite-chain reachability per leg, which the positive-control golden
answers. Uranus/Neptune stay (more firmly) empty.

---

## 9. Prerequisite corpus / environment acquisitions

**NONE REQUIRED.** Every primitive is built, tested, and golden-anchored to already-
digested material:

- Per-hop leverage leg: `search/leveraging_leg.py` (#179, phase-full apse VILM,
  validated against the Γ-floor and brute force).
- Quadrature floor + endgame golden: `search/vilm.py` (Campagnola-Russell Endgame
  Part-1, digest `2026-06-05-endgame-tisserand-mining.md`, golden to <10% band).
- Releg protocol, driver, prefilter, dv-band gate, empty re-stamp: all present
  (`search/releg_solver.py`, `releg_moontour.py`, `moon_prune.py`,
  `verify/dv_band_acceptance.py`).

No new paper, no install gate. Pure-Python on the existing scipy stack. (The
finite-chain *phasing* realisation — the integer n:m sequence — is computed from
the existing `leveraging_leg` feasibility; if a future tighter golden on the
DISCRETE chain were wanted, the published Europa 3-VILM/14-VILM ToF+ΔV breakdown is
already digested in A6 and suffices as a bracket.)

---

## 10. Open questions for the user

1. **Hop-selection policy:** greedy "zigzag low" descent (mining note p.11) picking
   the cheapest feasible n:m at each step, vs a small branch&bound over the
   leveraging graph (Part-1 pp.9-10)? Recommendation: **greedy first** (cheapest,
   matches the quadrature in the limit); branch&bound only if greedy stalls.
2. **ToF gating:** carry chain ToF to V2 only, or add a per-leg max-revs feasibility
   cap in the backend? Recommendation: **max-revs cap in the backend** (kills
   pathological deep walks cleanly) + carry ToF to V2.
3. **Scope:** ship Jovian+Saturnian in-band capability + the Uranus/Neptune
   stronger-empty re-stamp; defer the discovery CAMPAIGN (relegging the
   `repeated-moon-*-sweep` skeletons at scale) to a follow-on issue, as #449 did.
   Confirm.
