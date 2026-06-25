# Multi-rev resonant-leveraging releg — VERDICT (#465)

**Date:** 2026-06-26
**Design draft:** `docs/superpowers/specs/2026-06-26-465-multirev-leveraging-releg-design-draft.md`
**Plan:** `docs/superpowers/plans/2026-06-26-465-multirev-leveraging-releg-plan.md`
**Status:** SHIPPED — capability + golden + positive controls; honest reachability
characterisation. No catalogue row self-admitted (a candidate is FLAGGED for human
gauntlet review below, not promoted).

---

## 1. The key result (the #465 gate)

**The Galilean Io-Europa-Ganymede-Io cycle CLOSES IN-BAND with the multi-rev
leveraging chain.** Measured through the real `close_powered_cycle` driver with the
new `MultiRevLeveragingReleg` backend (NOT a cost estimate — an actual closed,
V∞-continuous cycle):

| Skeleton | Backend | Closes? | Continuity residual | Total ΔV/cycle | In band (<3.5 km/s)? |
|---|---|---|---|---|---|
| Io-Europa-Ganymede-Io | `MultiRevLeveragingReleg` | YES | 0.0 km/s (by construction) | **≈0.71 km/s** (driver picks cheapest T) | **YES** |
| Io-Europa-Ganymede-Io | `DsmReleg` (#449) | YES | exact | ~13 km/s (single-VILM max) | NO |
| Titan-Rhea-Dione-Titan (Saturnian) | `MultiRevLeveragingReleg` | YES (scale 1.2) | 0.0 km/s | **≈0.65 km/s** | **YES** |
| Ariel-Umbriel-Ariel (Uranian) | `MultiRevLeveragingReleg` | NO — prefilter EMPTY | — | inf | — (stronger powered-empty) |

The chain's per-cycle ΔV is **~18× cheaper** than the single-DSM close (0.71 vs
~13 km/s) and **~5× under** the 3.5 km/s powered ceiling. This is the in-band
closure #449 (13.18 km/s) and #464 (12.03 km/s) structurally could not reach,
because they shed the whole V∞ defect in ONE impulse (the single-VILM *maximum*,
Eq.14) where the chain walks it down across many resonant hops (the multi-VILM
*minimum*, Eq.13).

**Outcome 1 of the plan's three honest outcomes is realised for Jovian + Saturnian.**

---

## 2. Validation golden status (sourced, non-circular)

- **Per-hop floor:** every realised chain hop respects its own Γ-floor
  (`leveraging_leg.gamma_floor_ok`), so the summed chain ΔV is ≥ the sum of the
  per-hop Eq.(13) quadrature floors — a finite integer-resonance chain cannot beat
  the continuous minimum (`test_multirev_releg_dv_geq_leverage_floor`).
- **Europa endgame bracket:** the chain's Europa V∞ walk 1.8→0.77 km/s realises
  **≈137 m/s**, inside the sourced bracket `[128 m/s continuous floor (Eq.13),
  154 m/s published 3-VILM discrete (A6)]` — both EXPECTED edges trace to
  Campagnola-Russell, never a chain-computed number
  (`test_chain_walks_vinf_down_at_floor`). Step sensitivity: 0.1→139.85, 0.05→136.87,
  0.02→135.43, 0.01→134.74 m/s — monotone toward the floor as the hop step shrinks,
  always above it.
- **Decomposition golden:** `leverage-only + escape + capture == published Table-1
  ΔV_min` reproduces the printed Part-1 Table 1 values (Ganymede-Europa 1.715,
  Europa-Io 1.760, Rhea-Dione 0.519 km/s) to the paper's 2-sig-fig precision; the
  cycler's leverage-only per-transfer cost is the published total minus the sourced
  escape/capture, reproducing the design-draft §6 numbers exactly (Ganymede-Europa
  0.305, Europa-Io 0.409, Ganymede-Io 0.735, Rhea-Dione 0.196 km/s)
  (`test_golden_decomposition_reconstructs_published_table1`,
  `test_leverage_only_costs_match_design_section6`).

All goldens GREEN.

---

## 3. Finite-chain reachability verdict (the genuine open risk, design MEDIUM-HIGH)

The design flagged finite-chain reachability — *can integer-resonance hops realise
the in-band cost in a feasible number of revs?* — as the real open risk. The honest
finding:

- **It is real and it bites at HIGH arrival V∞ relative to the moon's V_M.** The
  greedy descent stalls immediately (zero feasible advancing hops) when the natural
  arrival V∞ exceeds roughly ~2·V_M at the flyby body: the high-V∞ resonant orbits
  no longer cross the moon for an exterior apse burn. Concretely, a Dione→Titan leg
  whose phasing forces a 12.2 km/s arrival at the heavy Titan (V_M≈5.6 km/s,
  V∞/V_M≈2.2) CANNOT be walked down — the chain returns infeasible (honestly, never
  a fabricated bridge).
- **It does NOT bite the Galilean cycle, and is a PHASING artifact for the
  Saturnian one.** At the Galilean moons the per-leg arrival V∞ (≈10-12 km/s at
  Europa, V_M≈13.7) is below the reach ceiling, so the chain walks every leg down
  (the cycle closes). For the Saturnian skeleton, a tof-scale of 1.0 produced a
  12.2 km/s arrival at Titan that stalled; scale 1.2 produced a ≈4-6 km/s arrival
  that the chain reaches cleanly. The high-V∞ stall was a phasing/ToF mismatch (the
  moon not where the transfer arrives), not a fundamental block — a discovery
  campaign that sweeps ToF/phasing finds the reachable skeletons.
- **Net:** the in-band cost is realised by a real finite chain for the Galilean and
  (phasing-tuned) Saturnian positive controls. The reachability ceiling
  (V∞ ≲ ~2·V_M for an exterior leveraging hop) is a genuine constraint the at-scale
  campaign must respect — skeletons whose phasing drives arrival V∞ above it are
  infeasible for the coplanar chain and stay empty.

---

## 4. Uranus / Neptune — stronger powered-empty (Outcome 3, unchanged)

The Uranian/Neptunian disjoint-Tisserand-contour systems are NOT rescued and the
backend does not pretend otherwise: the `moon_prune` prefilter skips the
unbridgeable legs BEFORE any chain solve (chaining walks V∞ *within* a contour, it
cannot jump disjoint ones). The `multirev_leveraging_method_capability` STRICTLY
SUBSUMES the single-DSM powered-empty (`multi-rev-leveraging` ⊐ `one-dsm-per-leg`),
so `build_powered_empty_restamp` records these as a *stronger* powered-empty
negative (`test_multirev_uranus_powered_empty_restamp`).

---

## 5. Catalogue candidate (FLAGGED for human gauntlet review — NOT self-admitted)

The in-band Galilean Io-Europa-Ganymede-Io leveraging-cycler is a genuine candidate.
Per discipline it is NOT promoted here; it is flagged for the human gauntlet:

- **Proposed class:** `mga_tour` / `quasi_cycler` (a flyby resonant tour; a cycler
  in the loose sense — it returns to Io and repeats — but powered, so not a free
  ballistic cycler).
- **Measured ΔV:** ≈0.71 km/s/cycle (driver-minimised over the common flyby target),
  `powered_dsm` band.
- **Lit-novelty check (the mandatory baseline,
  `feedback_literature_novelty_check_baseline`):** a Galilean resonant leveraging
  tour is very likely a **REPRODUCTION** of the published Campagnola/Strange Jovian
  endgame tours (Campagnola-Russell "The Endgame Problem" Part-1/2 explicitly design
  Galilean Io-Europa-Ganymede VILM endgames; the Strange-Russell Tisserand-graph
  tours cover the same moons). The decomposition golden in fact ANCHORS to those
  papers' Table 1. So the honest expectation is **V0-known (a reproduction of a
  published tour), not SILVER-novel.** `search/literature_check.py` must clear it
  before any admission; the strong prior is that it matches the published record.

**No catalogue row is written by #465.** The at-scale discovery campaign +
lit-check + V2-moontour gauntlet is the explicit follow-on (as #449 did).

---

## 6. What shipped

| Task | Commit | Content |
|---|---|---|
| 1 | `53e2fb1` | `search/leveraging_chain.py` — resonant-hop V∞ descent (`walk_vinf_down`) + tests |
| 2 | `fe2c137` | `MultiRevLeveragingReleg` backend behind the `Releg` protocol + `arrival_moon` protocol add + `chain_hops` on `RelegResult` |
| 3 | `ddd6574` | driver wires `arrival_moon=sequence[k+1]` into `close_powered_cycle`; Galilean in-band + Uranus-empty regression |
| 4 | `e14b074` | sourced in-band leveraging golden (decomposition) + Saturnian positive control |
| 5 | `9e3b44d` | `multi-rev-leveraging` capability edges + `multirev_leveraging_method_capability` + stronger powered-empty re-stamp |
| 6 | (this) | verdict doc + OUTSTANDING capability note |

Reuse honoured: the per-hop primitive (`leveraging_leg.py`, #179), the VILM cost
model / floor (`vilm.py`, golden), the `Releg` protocol + driver + prefilter +
dv-band gate + empty re-stamp — all imported, not re-derived. The new code is the
chain ORCHESTRATOR only.

---

## 7. Out of scope (follow-on issues)

- The at-scale discovery CAMPAIGN (relegging the `repeated-moon-*-sweep` skeletons,
  sweeping ToF/phasing to find the reachable in-band skeletons, re-stamping the
  registry) — this plan ships the capability + golden; the campaign spends it.
- 3D / inclined leveraging (the coplanar chain is the `leveraging_leg`/`vilm`
  regime).
- Discrete-chain ToF Pareto optimisation (the greedy descent suffices for the
  in-band proof).
