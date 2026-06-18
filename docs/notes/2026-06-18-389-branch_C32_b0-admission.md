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

## P389.2 — Physical-sanity + lit-fresh + ML flagger [PENDING]

Per-encounter V_∞ at the Moon for the 3 lunar encounters in the (3,3) cycle;
#324 max-bend gate (5° threshold); offline `check_literature` against the
KNOWN_CORPUS with `topology_label = frozenset({"repeated-moon"})`; #256 ML
flagger score.

---

## P389.3 — V2 bounded-cycle gate (CR3BP-adapted) [PENDING]

`run_v2_3d` over n_cycles ∈ {3, 5, 10} — DOES the orbit remain bounded under
successive DOP853 propagations WITHOUT re-correction? Spec §14 V2 floor =
50,000 km (same-model).

For an *essentially stable* (max Floquet ~= 1 + 6e-13) orbit the answer should
be a clean PASS — the per-cycle drift should grow only as the cumulative round-
off error × Floquet amplification factor, which for σ_d~6e-15/day is
effectively zero.

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
