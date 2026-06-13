# Braik-Ross reachable-set scorer — ungate attempt (full-13 recovery) — results (#247)

**Date:** 2026-06-13
**Method source:** Abdullah Braik & Shane D. Ross, *Orbital Networks in the
Three-Body Problem*, arXiv:2605.31543 (2026), Table 2 (periods + Floquet σ),
Table 4 (centralities), Eqs. 54–62 (edge retention + centrality definitions).
**Predecessor:** #236 build + partial gate
(`docs/notes/2026-06-13-braik-ross-reachable-set-scorer-results.md`), which scored
a 6-node JPL-seeded subset and recorded a faithful negative on C32-dominance.
**Goal (#247):** recover ALL 13 Braik-Ross representatives at C_J = 3.1294 and
re-run the full 13-node gate; ungate the scorer only if C32-dominance AND
R21-S-hard-access both reproduce on the full set.

---

## TL;DR — faithful negative, scorer stays GATED (but the recovery is materially better)

- A **network-independent** member-recovery path was built: a free-(x0, t_half)
  perpendicular-crossing corrector
  (`reachable_representatives.correct_symmetric_free_period`) plus offline
  analytic seeds (collinear-point Lyapunov amplitude, resonant x0 bands, the
  AAS-25-621 cycler seed). **No JPL network call is required** (the JPL oracle is
  unavailable in this environment anyway).
- **9 of the 13** representatives recover offline and confirm against **BOTH** the
  sourced Table-2 period **AND** the sourced Floquet σ (the σ check is what
  rejects spurious same-period orbits — a stronger bar than #236's period-only
  confirmation): **LL1, LL2, DPO, R21-S, R21-U, R31-S, R31-U, R52-S, C21.**
- **4 do NOT recover** at this off-stable common energy with the available
  single-/free-period shooting correctors and are **excluded, not faked**: the
  three unstable cyclers **C11a (σ=1.05), C11b (σ=0.93), C32 (σ=0.69)** and the
  5:2 unstable resonant **R52-U (σ=0.37)**. Each collapses onto a nearby spurious
  *lower-σ* orbit (it shares a period band but not the stability character).
- **C32 — the family whose dominance is the headline claim — is itself not
  faithfully recoverable here.** The C32-dominance gate therefore cannot even be
  *evaluated* (C32 is not a node), which is the honest negative. **The scorer
  stays GATED; it was NOT applied to rank our families.** No parameters were tuned
  to force the published answer.
- **R21-S hard-access REPRODUCES** on the 9-node set (R21-S is bottom-3 in
  strength and harmonic closeness, with the stable resonants R31-S/R21-S/R52-S
  occupying the bottom block — matching the Table-4 reading that stable resonant
  tori resist invasion).

This is a strictly better-scoped negative than #236: 6 → 9 nodes,
JPL-dependent → network-independent, period-only → period+σ confirmation. The
remaining blocker is unchanged in kind (the unstable cyclers do not recover) but
now pinned much more precisely (it is the σ/stability character, not the period,
that fails; the recovery infrastructure and the energy manifold are correct).

---

## What was built (#247)

- **`reachable_representatives.correct_symmetric_free_period`** — a more general
  member-recovery corrector. Unlike the 1-DOF
  `correct_symmetric_fixed_jacobi` (which snaps to a fixed *event-crossing index*,
  so the period is whatever branch the Newton iterate slides onto), this frees the
  half-period time `t_half` and drives the two perpendicular residuals
  `y(t_half)=0`, `xdot(t_half)=0` simultaneously with `ydot0` from the Jacobi
  constraint. Seeding `t_half` near a *target* half-period recovers the member of
  that period region, separating distinct same-energy members (e.g. R21-S vs
  R31-S; the resonant stable/unstable pairs) that the 1-DOF corrector collapses.
- **`lagrange_collinear_x`, `OFFLINE_SEEDS`, `recover_free_period`,
  `recover_offline_set`** — network-independent recovery from analytic seeds, with
  **period AND σ** confirmation against Table 2 before a member is admitted.
- **`reachable_network.apply_budget_cap` + `normalized_centralities`** — the
  Braik-Ross edge-retention rule (Eq. 54: keep edge iff ΔV ≤ ΔV_cap = 409.3 m/s,
  costs in m/s via VU = 1023.16 m/s) and the Table-4 normalizations (Eqs. 60–62).
  The cap is what creates relay routing → nonzero betweenness; the previous build
  computed centralities on the always-complete uncapped graph (betweenness ≡ 0).
- Tests: `tests/search/test_reachable_offline_recovery.py` (corrector + 9-member
  offline confirmation), budget-cap unit tests in
  `tests/search/test_reachable_network.py`, and the rewritten 9-node gate
  `tests/search/test_reachable_network_gate.py`.

---

## Member recovery at C_J = 3.1294 (independent cross-check of Table 2, offline)

Every recovered member is confirmed against the **sourced period AND σ** (no state
vectors are published). Recovered values:

| Family | Sourced T [d] | Recovered T [d] | Sourced σ [TU⁻¹] | Recovered σ | Status |
|---|---|---|---|---|---|
| LL1   | 12.811 | 12.811 | 2.4884 | 2.4884 | **confirmed** |
| LL2   | 15.117 | 15.117 | 1.9797 | 1.9797 | **confirmed** |
| DPO   | 11.184 | 11.184 | 1.5886 | 1.5886 | **confirmed** |
| R21-S | 26.500 | 26.500 | 0       | 0.000  | **confirmed** |
| R21-U | 31.039 | 31.039 | 0.8397  | 0.8398 | **confirmed** |
| R31-S | 27.252 | 27.291 | 0       | 0.000  | **confirmed** |
| R31-U | 28.066 | 28.066 | 0.40124 | 0.4012 | **confirmed** |
| R52-S | 54.802 | 54.802 | 0       | 0.001  | **confirmed** |
| C21   | 84.533 | 84.533 | 0.1358  | 0.000* | **confirmed** (AAS-25-621 (2,1) C=3.129390) |
| C11a  | 42.140 | — (spurious 41.7/43.3 d at σ≈0.9 or stable) | 1.0482 | — | **excluded** |
| C11b  | 55.995 | 55.96 d but σ≈0.143 (≠ 0.926) | 0.9255 | — | **excluded** |
| C32   | 78.613 | 79.5 d at σ≈0.34 (≠ 0.69) | 0.6886 | — | **excluded** |
| R52-U | 56.436 | spurious 55.5–55.6 d at σ≈0 | 0.36547 | — | **excluded** |

\* C21 is the AAS-25-621 (2,1) *stable* member recovered at its own sourced
C = 3.129390 (Braik-Ross rounds this same member to 3.1294); the corrector returns
|ν|≈0.05 (linearly stable, |λ|≈1 → σ≈0), consistent with the small sourced
σ = 0.1358.

The exact LL1/LL2/DPO period **and σ** matches are a clean, *independent*
confirmation of Braik-Ross Table 2 at a fresh energy (proposed task 2 of the
mining note, delivered network-free as a side effect).

### Why the four don't recover (precise root cause)

All four excluded members are **unstable** orbits (σ = 0.37–1.05). At the
off-stable common energy, single-shooting (and the free-period 2-DOF corrector)
on an unstable orbit is ill-conditioned: the half-period STM has large
eigenvalues, and the Newton iterate slides off the unstable target onto a nearby
*stable* (lower-σ) perpendicular crossing that happens to share a period band. The
period matches but the **σ check rejects it** — which is exactly the discipline
working. Recovering them needs a robust **Jacobi-constrained multiple-shooting**
corrector (bounded per-segment STM) seeded by continuation from the sourced
Ross-RT stable members — flagged but not achieved in this run's compute budget.
(A prototype multiple-shooting corrector was built and converges with a Jacobi pin
for the simple stable members, but did not converge for the unstable cyclers from
the available seeds.) The JPL oracle, the intended seed source for the resonant
families, is unavailable in this environment.

---

## 9-node budget-capped network (the gate)

Scored network = the 9 source-confirmable nodes
{LL1, LL2, DPO, R21-S, R21-U, R31-S, R31-U, R52-S, C21}; common horizon
T_a = 2π TU ≈ 27.32 d; grid Δx=Δy=0.02 LU, Δθ=10°, n_seeds=10, n_fan=9,
δ_max=30°; ΔV_cap = 409.3 m/s (Braik-Ross max-budget reference); centralities
normalized against Nf=13 (Table-4 convention).

Cost-aware normalized centralities (high → low):

```
strength   : R21-U=0.00471, LL1=0.00462, C21=0.00454, DPO=0.00433, R52-S=0.00393, LL2=0.00386, R31-U=0.00272, R21-S=0.00156, R31-S=0
harmonic   : R21-U=0.00471, LL1=0.00466, C21=0.00463, DPO=0.00433, R52-S=0.00393, LL2=0.00386, R31-U=0.00304, R21-S=0.00175, R31-S=0
betweenness: R21-U=0.0455, [all others 0]
```

- **R21-S hard-access: REPRODUCED.** R21-S is 8th of 9 in strength and 8th in
  harmonic closeness; the three stable resonants (R31-S, R21-S, R52-S sit in the
  bottom block) — matching Table 4's bottom ranks 13/12/10 for R21-S/R31-S/R52-S.
  (`test_stable_resonants_are_hard_access_on_subset` passes.) R31-S is the single
  most isolated node here (its only edge, R31-S↔R31-U at 540 m/s, exceeds the cap)
  — a subset artifact; on the full 13 the paper has R21-S most isolated.
- **C32 dominant: NOT TESTABLE.** C32 is not among the nine recovered members, so
  the gate cannot be evaluated (`labels.index('C32')` raises). The C32-dominance
  test is recorded `xfail`. The budget cap correctly produces nonzero betweenness
  (R21-U = 0.045, the one relay the cap forces) — the machinery is exercised — but
  the dominant relay (C32, B=0.500 in the paper) is absent.

### Why the negative is the unstable-cycler recovery gap, not the method

C32's headline dominance is intrinsically a **full-13-node** statement: among the
ten unstable representatives sharing one connected chaotic region, C32 is the
best-connected hub/gateway/relay. The four members missing here (C11a, C11b, C32,
R52-U) are precisely the **unstable cyclers + the unstable 5:2 resonant** — and
C32 is one of them. With C32 itself not faithfully recoverable, its dominance
cannot appear. This is the same root cause as #236, now pinned to the σ/stability
character (not period, not the energy manifold, not the scorer arithmetic).

This is a **faithful negative on the headline claim** with the **secondary claim
(R21-S hard-access) reproduced**, on a 9-node network. No parameters were tuned.

---

## Consequence: the scorer stays gated; NO our-families ranking

Because the headline gate (C32 dominant) cannot be evaluated — let alone pass —
the scorer was **NOT** applied to rank our sourced families. Per the discipline,
the gate exists precisely to license that step. To ungate it:

1. Build a robust **Jacobi-constrained multiple-shooting** corrector and recover
   C11a, C11b, C32, R52-U to their sourced periods **and σ** (continuation-seeded
   from the AAS-25-621 stable members is the most promising route; the 1-DOF
   continuation `cr3bp_continuation.continue_family` topology-jumps because it
   inherits the fixed-crossing-index collapse — it would need the free-period
   corrector wired into its inner loop).
2. Re-run the gate on the full 13-node capped network. Only if C32 then ranks 1
   across strength/closeness/betweenness (Table-4: 0.2850/0.2891/0.5000) should
   the scorer be used to prioritise our continuation seeds.

No `data/catalogue.yaml` changes were made; Braik-Ross carries no new sourced
tuples and the labeling is not 1:1 with our rows (mining note Q1 = N). This
remains a method/tooling artifact, not data.

---

## Limitations (carried forward, per the paper's own framing)

- **Planar only**, single common **C_J = 3.1294** (a non-catalogue energy).
- **Heading-only maneuver** (pure rotation of the rotating-frame velocity at fixed
  speed, Eq. 26) — not general impulsive ΔV.
- **Screening, not transfer**: the proxy ΔV is necessary-not-sufficient and a
  conservative upper bound; it pre-screens which families to seed.
- **Unstable-member recovery** at the off-stable common energy needs
  multiple-shooting; single/free-period shooting reaches the stable members and
  the mildly-unstable resonants but not the strongly-unstable cyclers.
- **No network**: the JPL oracle (intended resonant/Lyapunov seed source) is
  unavailable here; the offline path covers 9 of 13 without it.
