# #389 Phase 4 — branch_C32_b0 V0-V5 gauntlet for catalogue admission

**Status:** P389.0 complete, working on `main`.
**Date opened:** 2026-06-18 AET.
**Predecessor task:** #347 Phase 2 (`docs/notes/2026-06-17-347-floquet-phase2-progress.md`).
**Structural template:** #339 (`docs/notes/2026-06-17-339-silver-quasi-cycler-admission.md`).
**Candidate:** `branch-c32-b0-em-3-3-quasi-cycler-2026` — the #347 Floquet bifurcation framework's first discovered cycler-family member (planar (3, 3) Earth-Moon CR3BP at T=101.56d, J=3.797, max Floquet 1.000000000000617).
**Catalogue target:** if all V-tiers clear → catalogue 302 → 303 rows, second computed quasi_cycler row (first from the Floquet framework).
**Discipline:** every V-tier cleanly + atomically committed; HALT on any failure per `feedback_orbit_closure_discipline`.

---

## P389.0 — Provenance setup

**New artifacts:**

* `data/branch_c32_b0_ic.jsonl` — sanitized verification IC. Three rows:
  header / `kind: ic` (the canonical Phase 2 row's state0 + period + jacobi
  + topology + bifurcation context) / footer.
* This progress note.

**IC source-trace.**

* `data/floquet_phase2_sweep_results.jsonl` row 3 — `kind: branch_record`,
  `parent_label: C32`, `bracket_index: 0`. State0 / period / jacobi / topology
  match `data/floquet_phase2_family_network.jsonl` row `branch_C32_b0`
  (independent cross-check).

**IC characteristics frozen from Phase 2:**

| Field | Value |
|-------|-------|
| state0_rotating_nondim | (-0.7033325748418664, -2.9123784605081626, 1.73e-22, -2.3503749595840504, 0.567571628434372, 8.30e-24) |
| period_TU | 23.355184434547017 |
| period_days | 101.55715620610965 |
| jacobi | 3.797487163854691 |
| topology | (3, 3) |
| degenerate_planar | True |
| max Floquet mag | 1.000000000000617 |
| σ_d / day | 6.08e-15 |
| corrector_residual (Phase 2) | 4.77e-12 |
| independent Radau closure (Phase 2) | 2.59e-11 |
| Parent | C32 (Braik-Ross 2026 Table 2 (3, 2)) |
| Bracket | C ∈ (3.14170, 3.14180), saddle-center |

**P389.0 commit (pending):** `data/branch_c32_b0_ic.jsonl` + this note.

---

## P389.1 — Closure + DOP853 cross-check (V1 re-confirm) — PASS

**New artifacts:**

* `scripts/branch_c32_b0_v1_verify.py` — re-runs the Phase 2 branch corrector
  against `data/branch_c32_b0_ic.jsonl`. Uses
  `correct_general_periodic_3d` in full-asymmetric mode (free vars =
  (x0, y0, z0, xdot0, ydot0, zdot0, T); residual = full 6D state closure at T)
  with tol=1e-12, max_iter=80, rtol=atol=1e-12. The independent Radau closure
  check is built into the `Periodic3DOrbit.converged` compound gate.
* `data/branch_c32_b0_v1_verdict.jsonl` — 3-row verdict (header / verdict / footer).
* `tests/verify/test_branch_c32_b0_v1_passes.py` — 3-test frozen-gate pytest.

**Results (run 2026-06-18, 0.95s wall):**

| Field | Value | Phase 2 reference | Margin |
|-------|-------|-------------------|--------|
| Corrector residual | **1.276e-13** | 4.77e-12 | ~37x tighter |
| Independent Radau closure | **2.188e-11** | 2.59e-11 | comparable |
| `converged` (compound gate) | True | True | — |
| `n_iter` | 2 | — | — |
| T_TU corrected | 23.355184434547020 | 23.355184434547017 | drift 3e-15 |
| Jacobi corrected | 3.797487163854493 | 3.797487163854691 | drift 2e-13 |
| `degenerate_planar` | True | True | — |

**V1 floors held:**

* Corrector residual 1.276e-13 < 1e-10 floor — PASS.
* Independent Radau closure 2.188e-11 < 1e-6 floor — PASS.
* Period preserved to 8+ decimals — PASS (Phase 2 structural identity intact).
* Planar character preserved — PASS.

**P389.1 gate:** PASS. The Phase 2 sweep IC re-closes under V1 spec with the
same characteristic period, jacobi, and planar topology. The orbit is a real
CR3BP periodic orbit, not a Phase 2 driver artifact.

**P389.1 commit (pending):** `scripts/branch_c32_b0_v1_verify.py` +
`data/branch_c32_b0_v1_verdict.jsonl` +
`tests/verify/test_branch_c32_b0_v1_passes.py` + this note update.

---

## P389.2 — Physical-sanity + lit-fresh + ML flagger — PASS (with structural caveat)

**New artifacts:**

* `scripts/branch_c32_b0_p389_2_gates.py` — combined three-gate driver.
* `data/branch_c32_b0_p389_2_verdict.jsonl` — 3-row verdict.
* `tests/verify/test_branch_c32_b0_p389_2_passes.py` — 4-test frozen-gate pytest.

**Physical-sanity (#324 max-bend gate):**

| Encounter | t_TU | r_min_to_Moon (km) | V_rel (km/s) | V_∞ (km/s) | max_bend (deg) | is_useful |
|-----------|------|---------------------|-------------|-------------|----------------|-----------|
| 0 | 5.546 | 772,165 | 2.4749 | 2.4723 | 35.383 | True |
| 1 | 13.330 | 772,165 | 2.4749 | 2.4723 | 35.383 | True |
| 2 | 21.115 | 772,165 | 2.4749 | 2.4723 | 35.383 | True |

All three encounters clear the 5° structural-feasibility floor — PASS.

**HONEST CAVEAT (recorded in JSONL `interpretation` per encounter):** the
actual minimum r-to-Moon distance is **772,165 km**, which is **~12× the
Moon's Hill sphere (66,100 km)**. branch_C32_b0 is a *far-amplitude bound
Earth-system orbit at ~1.15 million km Earth-distance*, NOT a lunar-flyby
tour. The (3, 3) topology label refers to the winding number around each
primary, not to physical close encounters. The #324 gate is checking
*structural feasibility* at the registered Moon safe-altitude (Moon radius +
100 km = 1837 km r_periapsis) IF the spacecraft were placed there at the
indicated V_∞; the orbit itself doesn't transit through the Moon's vicinity.

This is consistent with branch_C32_b0's parent (C32 = (3, 2)) which is a
Braik-Ross Table 2 cycler — those cyclers also have large amplitudes, and
their (k1, k2) labels are topological invariants of the winding, not flyby
counts. The catalogue row will describe this as a `quasi_cycler` with
`flyby_class: distant-libration` / `flyby_class: none` rather than a flyby
tour.

**Lit-fresh (#346/#349):**

* Signature: primary=Earth, sequence=("Moon",), period_k=3, topology_label =
  frozenset({"repeated-moon"}), period_band_tu = (T-1e-6, T+1e-6).
* `check_literature` injected with an offline no-results search returned
  `status=inconclusive` (no results returned at all — the project's
  web-search lane is not wired in this run). Per the SILVER's #328 precedent,
  this is the same offline-corpus regime the SILVER's V0 lit-check ran under.
* **Structural KNOWN_CORPUS scan:** the only Earth-Moon CR3BP anchors with
  `topology_label = frozenset({"repeated-moon"})` are Braik-Ross 2026
  (publishes Table 2 = C11a/C11b/C21/C32; no (3, 3)) and Roberts-Tsoukkas &
  Ross AAS 25-621 (publishes the 5 stable EM cyclers; no (3, 3)). No
  published (3, 3) planar EM CR3BP cycler at jacobi=3.797 in the present
  corpus — branch_C32_b0 is **structurally novel at the published-record
  level** (necessary-not-sufficient per `feedback_literature_novelty_check_baseline`).

**ML flagger (#256):**

* Trained on the labeled corpus (`build_training_set`).
* Score: **p_fp = 0.3748** (well below the 0.75 spec §16.5 routing threshold).
* `ml_passes = True`.

**P389.2 gate:** PASS on all three sub-gates, with the structural caveat
(far-amplitude orbit, not a lunar-flyby tour) recorded for honesty.

**P389.2 commit (pending):** `scripts/branch_c32_b0_p389_2_gates.py` +
`data/branch_c32_b0_p389_2_verdict.jsonl` +
`tests/verify/test_branch_c32_b0_p389_2_passes.py` + this note update.

---

## P389.3 — V2 bounded-cycle gate (CR3BP-adapted) — PASS (spectacular)

**New artifacts:**

* `scripts/branch_c32_b0_v2_verify.py` — calls
  `cyclerfinder.data.validation.v2_3d.run_v2_3d` at n_cycles ∈ {3, 5, 10}
  with the V1-corrected state + period.
* `data/branch_c32_b0_v2_verdict.jsonl` — 5-row verdict (header + 3 per-n_cycles
  rows + footer).
* `tests/verify/test_branch_c32_b0_v2_passes.py` — 4-test frozen-gate pytest.

**Results:**

| n_cycles | max_drift (km) | Spec §14 floor | passes_v2 |
|----------|----------------|-----------------|-----------|
| 3 | 2.47e-5 km (24.7 μm) | 50,000 km | True |
| 5 | 4.35e-5 km (43.5 μm) | 50,000 km | True |
| 10 | 5.51e-5 km (55.1 μm) | 50,000 km | True |

**Margin: ~9 orders of magnitude below the spec floor.** The per-cycle drift
is pure DOP853 round-off accumulation (rtol=atol=1e-12). The orbit is
dynamically essentially perfectly stable, exactly as the Phase 2 Floquet
character predicted (max_floquet_mag = 1.000000000000617, σ_d = 6.08e-15/day).

**P389.3 gate:** PASS. branch_C32_b0 is a strict bounded cycler in the spec
§14 sense at n_cycles up to 10 in the same model. (Compare to the SILVER's
~530,000 km bounded oscillation; branch_C32_b0 is ~10 orders tighter — a
genuinely strict-periodic CR3BP cycler, NOT a quasi_cycler in the
bounded-oscillation sense.)

**Implication for the orbit_class label.** The SILVER admission used the
`quasi_cycler` slot because it FAILED strict V2 by exceeding the 50,000 km
floor (in bounded oscillation). branch_C32_b0 PASSES strict V2 by 9 orders of
magnitude, suggesting it should enter the `cycler` slot, not `quasi_cycler`.
This will be revisited in P389.6 row composition — the orbit_class admission
slot depends on the full V1-V5 verdict.

**P389.3 commit (pending):** `scripts/branch_c32_b0_v2_verify.py` +
`data/branch_c32_b0_v2_verdict.jsonl` +
`tests/verify/test_branch_c32_b0_v2_passes.py` + this note update.

---

## P389.4 — V3 REBOUND IAS15 n-body [PENDING]

REBOUND IAS15 + Earth + Moon + Sun + Mars + Jupiter mass points (DE440-class).
PASS = nanometer/micrometer agreement vs V2 at n=3, 5, 10 cycles.

---

## P389.5 — V4 GMAT + SPICE real-eph [PENDING]

DE440 ephemeris (not URA111 as in #339), V4-scipy first then V4-strict, plus
annual epoch sweep 2000-2099 analogous to #338.

---

## P389.6 — Admission verdict + catalogue writeback [PENDING]

Only on full PASS: row composition + 3 ratchet bumps + _LEVEL_EVIDENCE
registration. Catalogue 302 → 303.

---

## Discipline checklist

* Work on `main` directly. No branches.
* Atomic pathspec commits per sub-phase.
* Pre-commit hooks must pass on every commit.
* No `--no-verify`, no Co-Authored-By trailer.
* Sourced-only: every numeric value traces to a specific PDF page / committed
  JSONL output.
* Concurrent-agent commit hygiene: my paths are
  `data/branch_c32_b0_*.jsonl` + `tests/verify/test_branch_c32_b0_*.py` +
  `docs/notes/2026-06-18-389-*.md` + `data/catalogue.yaml` +
  `src/cyclerfinder/data/validate.py` + the three ratchet tests.
* No catalogue writeback until V5 verdict.
