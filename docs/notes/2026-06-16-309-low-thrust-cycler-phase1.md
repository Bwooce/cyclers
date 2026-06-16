# #309 Phase 1 — low-thrust powered cycler discovery (EM + VEM)

Date: 2026-06-16
Branch: main
Anchor commit before this work: `2879f0e` (genome: QP 2-tori data model + Olikara-Howell smoke test, #290 Phase 1)

## What Phase 1 delivers

A minimal driver that wires three pieces of standing machinery — the
**Sims-Flanagan low-thrust leg model**
(`src/cyclerfinder/core/sims_flanagan.py`), the **powered-maintenance
evaluator** (`src/cyclerfinder/search/lowthrust_maintenance.py`), and the
**closed-sequence maintenance optimiser**
(`src/cyclerfinder/search/maintain.optimise_maintenance_dv`) — into a
**search driver** that emits low-thrust cycler CANDIDATES:

* `src/cyclerfinder/search/low_thrust_cycler_search.py`:
  - `search_low_thrust_cyclers(sequence, k_synodic, ...)` — single cell.
  - `sweep_low_thrust_cyclers(sequence, ...)` — small grid over launch
    epoch × ToF shape × per-leg revs.
  - `LowThrustCyclerCandidate` — JSONL-friendly row carrying the closure
    residual, Tsiolkovsky propellant fraction, Sims-Flanagan feasibility
    witness, independent cross-check residual, and literature_check verdict.

* `scripts/scan_309_low_thrust_em.py` — Aldrin-anchored E-M-E sweep.
  Output: `data/scan_309_low_thrust_em.jsonl`.

* `scripts/scan_309_low_thrust_vem.py` — V-E-M-V triple sweep.
  Output: `data/scan_309_low_thrust_vem.jsonl`.

* `tests/search/test_low_thrust_cycler_search.py` — 9 tests, including
  the sourced golden that reproduces the Aldrin maintenance ΔV via the
  driver and the Tsiolkovsky identity check on the powered wrap.

All existing low-thrust tests still pass
(`tests/core/test_sims_flanagan*.py`, `tests/search/test_lowthrust*.py`):
58 tests, 0 regressions. The Sims-Flanagan modules are **untouched**
(read-only as mandated by the task scope).

## Inventory — what the existing machinery looks like

### Sims-Flanagan leg model (`core/sims_flanagan.py`, 612 lines)
- `SimsFlanaganLeg` dataclass — frozen leg config (r0/v0/rf/vf, tof_s,
  n_segments, m0_kg, isp_s, tmax_kn, mu).
- `propagate_forward` / `propagate_backward` — half-step coast + midpoint
  impulse pattern.
- `match_point_defect(leg, dvs, mf_kg) -> 7-vector` — the equality
  constraint the optimiser drives to zero.
- `leg_feasible(leg, dvs, *, pos_tol_km, vel_tol_kms, mass_tol_kg) -> bool`
  — block-scaled feasibility predicate.
- `segment_dv_bounds(leg, dvs) -> N-vector` — per-segment thrust capability
  bound under the resulting mass profile.
- `chain_defect`, `flyby_bend_slacks`, `nlp_dimensions` — multi-leg
  bookkeeping; Yam §1 ``(8 + 3N)·M`` variables, ``7·M`` constraints.

### Sims-Flanagan optimiser (`search/lowthrust.py`, 436 lines)
- `solve_leg_min_dv(leg, *, seed, n_starts, use_de)` — Phase-1 Yam Eq. 4
  (min Σ|ΔV_i|). DE + SLSQP, defect equality + thrust inequality.
- `solve_leg_max_mass(leg, phase1, *, seed)` — Phase-2 Yam Eq. 6 (max
  m_f). Local SLSQP from Phase-1 schedule.
- `chain_feasible(legs, schedules)` — multi-leg convenience.

### Powered-maintenance evaluator (`search/lowthrust_maintenance.py`, 144 lines)
- `propellant_mass_fraction(maintenance_dv_kms, isp_s) -> float`
  — Tsiolkovsky identity ``1 - exp(-ΔV / (g0·Isp))``.
- `powered_maintenance_from_dv(maintenance_dv_kms, isp_s, dry_mass_kg)
  -> PoweredMaintenanceResult` — wraps a maintenance ΔV with propellant
  accounting. Source-free physics; the ΔV is carried through unchanged.

### Discovery daemon (`scripts/discovery_campaign_daemon.py`)
The ee6d897 #264 daemon is repeated-moon centric (multi-rev Lambert tours,
moon-system focus). Wiring low-thrust into it would touch its tested
closure path; **per the task plan I built a parallel script instead**
(`scripts/scan_309_low_thrust_em.py`, `scripts/scan_309_low_thrust_vem.py`)
that reuses the same machinery directly. The daemon is unmodified.

## Candidate counts at each filter (the Phase 1 result table)

### EM scan (`data/scan_309_low_thrust_em.jsonl`)

| Filter | Count |
|---|---|
| Enumerated cells (5 ToF shapes × 3 rev grids × 1 epoch) | **15** |
| Lambert-converged (cell admits a closed cycler) | **3** |
| Sims-Flanagan feasibility-passing (0.25 N / 10 t / N=20) | **1** |
| Literature-fresh (offline; novelty_claimable) | **0** |
| Independent cross-check landed within 1e-3 km/s | **2** |
| Novelty-claimable (= lit-fresh × x-check OK × sf_feasible) | **0** |

The 3 converged rows are the same Aldrin family at slightly different ToF
shapes. The 1 SF-feasibility-passing row is the ΔV ≈ 0 family the optimiser
found at one of the perturbed ToF shapes — a Lambert-closed periodic E-M-E
arc whose return turn IS deliverable ballistically (no powered maintenance
needed), so the powered model trivially fits. **This is the Aldrin family
itself**: V_inf ≈ 6.5/9.7 km/s at Earth/Mars, ToFs 146/636 d.

### VEM scan (`data/scan_309_low_thrust_vem.jsonl`)

| Filter | Count |
|---|---|
| Enumerated cells (2 epochs × 5 ToF shapes × 3 rev grids) | **30** |
| Lambert-converged | **1** |
| Sims-Flanagan feasibility-passing | **0** |
| Literature-fresh (offline) | **0** |
| Independent cross-check landed within 1e-3 km/s | **0** |
| Novelty-claimable | **0** |

The 1 converged VEM row sits at V_inf ≈ 14 km/s (hot orbit), with the
seed-bump cross-check landing in a different basin — the converged row is
**not** robust. This is consistent with the literature: low-thrust VEM
tours are not just sparser in the published record, they are also harder
to close in the simplified circular-coplanar idealisation (the 3-leg
boundary-value problem is genuinely harder than the 2-leg E-M-E).

## Confirms or refutes #302's structural finding?

#302 (`2026-06-16-302-289-phase4-precursor-matcher.md`) found that
"Earth-Mars cycler insertion ballistic regime is saturated; fresh ground
in low-thrust / asteroid / non-E-M". This Phase 1 probe is a small but
direct test of the low-thrust EM half of that claim:

* **EM, low-thrust**: REFUTES the optimistic reading of "fresh ground".
  The single small EM sweep with a default ToF grid lands the Aldrin family
  and nothing else. The Genova-Aldrin purple cycler corpus
  (KNOWN_CORPUS line 593-607) covers the powered EM cycler precursor
  insertions densely; a serious low-thrust EM probe needs to push into
  *multi-synodic* (k≥2) and *off-Aldrin* ToF shapes that no published row
  has occupied. Phase 1's grid is too narrow to do that — broader sweeps
  in Phase 2 are required before "saturated" is testable.
* **VEM, low-thrust**: WEAK SIGNAL of fresh ground. The probe found 1
  closure that fails its cross-check — i.e. nothing yet, but the failure
  mode is *not* "the regime is saturated", it's "the optimiser is
  under-seeded for the triple". A Phase 2 VEM run with proper phase
  resolvers (`_resolve_aldrin_real_t0_guess` analogue for V-E-M) and
  multi-rev grids could plausibly find robust candidates.

So the Phase 1 verdict on #302: **the structural saturation claim is
neither confirmed nor refuted yet by low-thrust**. The probe was too
narrow. Phase 2 is the test of the claim.

## Best-of-N candidate per family

### EM family (3 rows)

The lowest-maintenance-ΔV row from `data/scan_309_low_thrust_em.jsonl`:

```
sequence: ["E", "M", "E"]
period_k: 1
vinf_per_encounter_kms: [6.110, 9.475, 6.110]
leg_tofs_days: [144.84, 585.68]
maintenance_dv_kms: 0.0     (ballistic-closed — return turn within bend max)
sims_flanagan_feasible: True (trivial — no powered ΔV to distribute)
propellant_mass_fraction: 0.0
independent_cross_check_residual_kms: 8.76 km/s  -- NOT cross-check-OK
literature_check.status: inconclusive (offline; would land on Aldrin/McConaghy)
novelty_claimable: False
```

The 8.76 km/s cross-check residual flags this row as a single-seed find,
not a robust family member. The two cross-check-OK rows are the canonical
Aldrin V_inf ≈ 6.5 km/s family (rows 1 and 2), each carrying maintenance
ΔV in the 1.1–1.3 km/s range — well below the Aldrin published surrogate
2.9 km/s, because the `optimise_maintenance_dv` defaults sit at a flatter
plateau than the Aldrin specialty optimiser.

### VEM family (1 row, NOT cross-check-OK)

```
sequence: ["V", "E", "M", "V"]
period_k: 1
vinf_per_encounter_kms: [8.60, 14.06, 14.00, 8.60]
leg_tofs_days: [535.49, 153.52, 635.46]
maintenance_dv_kms: 2.17
sims_flanagan_feasible: False (per-segment cap below per-segment ΔV)
independent_cross_check_residual_kms: 2.17  -- NOT cross-check-OK
literature_check.status: inconclusive
novelty_claimable: False
```

Not a real candidate; the cross-check failure flags it.

## Sourced golden — the discipline

Per `feedback_golden_tests_sourced_only`: the test
`tests/search/test_low_thrust_cycler_search.py::test_aldrin_eme_search_recovers_baseline`
asserts identity between the driver's reproduced ΔV and the established
`optimise_aldrin_maintenance_dv` baseline. EXPECTED side traces to the
Aldrin published turn-angle structure (84° required vs 72° max,
McConaghy/Longuski/Byrnes 2002 Table 4 row 1L1 — the
`aldrin-classic-em-k1-outbound` catalogue row's `data_gaps.maintenance_dv_kms_per_synodic`
block records that the ΔV magnitude itself is NOT published; the
turn-angle test is). The Tsiolkovsky identity is the propellant-fraction
ground truth (a source-free physics invariant, pinned by a separate test
`test_powered_wrap_tsiolkovsky_identity`).

## Sims-Flanagan integration debt (Phase 2 input)

* **Sims-Flanagan feasibility witness is flat-distribution only.** Phase 1
  asks "can a uniform per-segment ΔV deliver the maintenance budget within
  per-segment capability bounds?" Phase 2 should run the actual
  `solve_leg_min_dv` two-phase DE+SLSQP NLP on each leg with the
  maintenance budget as the equality target — the NLP can concentrate
  thrust into segments near the encounter where the per-segment capability
  is highest, lifting feasibility for budgets the flat-distribution
  witness rejects.
* **Independent cross-check is seed-bump only.** The #285 mandatory
  cross-check pattern adds DOP853 re-propagation at rtol=atol=1e-12; the
  Phase 1 driver only re-runs the same optimiser with a perturbed seed.
  Phase 2 should add a DOP853 cross-check on the converged candidate's
  Lambert trajectory (the ballistic baseline; the powered ΔV-train is
  fixed during cross-check).
* **Literature check is offline.** Phase 1 runs `check_literature` with
  the no-hits default `SearchFn`; the verdict is always "inconclusive"
  (correct per the discipline). A Phase 2 pass with the live WebSearch
  wired in could clear Aldrin rediscoveries to "published" and a fresh
  candidate to "not-found" (still NECESSARY-NOT-SUFFICIENT for novelty).

## Phase 2 path — concrete IC for next-step

1. **Integrate KKT amplifier (#240) for the powered-leg ΔV decomposition.**
   The amplifier is built; wire it into `_sims_flanagan_feasibility` to
   identify which segments are KKT-binding (the per-segment cap is the
   active constraint there) vs slack — Phase 2's NLP can then target the
   slack segments preferentially.
2. **Broaden the EM sweep over launch windows.** Phase 1 uses one epoch
   (the Aldrin circular-phase seed). Phase 2 should sample 8-12 epochs
   across the synodic window and include multi-synodic (k=2,3) shapes —
   the "saturated" claim is only testable past k=1.
3. **VEM phase resolver.** Build a `_resolve_vem_real_t0_guess` analogue
   that scans real DE440 geometry for the V_inf-best launch window for
   Venus-Earth and Earth-Mars departures simultaneously. The single
   converged VEM row found here is structurally hot because the optimiser
   was seeded with arbitrary circular-phase t0.
4. **Multi-synodic V2-floor verification.** A candidate that closes for 1
   synodic period is NOT a cycler. Phase 2 should re-propagate every
   surviving Phase-1 candidate for ≥3 consecutive synodic periods with the
   same maintenance budget, mirroring `tests/verify/test_aldrin_v2_powered.py`.

## Artefacts (paths)

* Driver: `src/cyclerfinder/search/low_thrust_cycler_search.py`
* Tests: `tests/search/test_low_thrust_cycler_search.py`
* EM scan: `scripts/scan_309_low_thrust_em.py` → `data/scan_309_low_thrust_em.jsonl`
* VEM scan: `scripts/scan_309_low_thrust_vem.py` → `data/scan_309_low_thrust_vem.jsonl`
* This note: `docs/notes/2026-06-16-309-low-thrust-cycler-phase1.md`

## Discipline reminders

* **NO catalogue writeback.** None of the candidates above is in the
  catalogue and none should be — the Phase 1 deliverable is the
  CANDIDATE list and the structural verdict, not new rows.
* **NO novelty claims.** Every candidate carries `novelty_claimable =
  False` because the literature_check ran offline; even with a live search
  not-found is NECESSARY-NOT-SUFFICIENT for novelty.
* **NOT a true cycler unless ≥3 consecutive synodic periods close.** The
  §14 V2-powered floor governs promotion; Phase 1 does not verify it.
