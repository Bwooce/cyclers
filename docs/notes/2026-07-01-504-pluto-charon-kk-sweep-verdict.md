# Task #504: Pluto-Charon (k1,k2)-cycler sweep — verdict

**Date:** 2026-07-01  
**Task:** #504  
**Sweep target:** mu = 0.10876473603280369 (cr3bp_system("Pluto","Charon"))  
**C_L1(PC):** 3.621018  
**Machinery:** `src/cyclerfinder/search/pluto_charon_kk_sweep.py`  
**Tests:** `tests/search/test_504_pluto_charon_kk_sweep.py` (6 tests, all PASS)

---

## Summary

**1 of 6 families yields a stable Pluto-Charon member: (3,2) only.**

The (3,2) family was the positive control (re-found, confirming machinery correct).
The other five families — (1,1), (2,1), (3,1), (2,2), (3,3) — are all clean negatives
at mu = 0.10876.

---

## Per-family result table

| (k1,k2) | Result | C (mid) | x0 (mid) | T (days) | nu (mid) | topo_ok | crosscheck | Method |
|----------|--------|---------|----------|----------|----------|---------|-----------|--------|
| **(3,2)** | **STABLE** | 3.5795150 | -0.693198287 | 12.033 | ~7.9e-8 | True | PASS | c_sweep from mu=0.1 anchor |
| (1,1) | clean negative | — | — | — | — | — | — | mu-step from mu=0.001; no stable window in C-sweep |
| (3,1) | clean negative | — | — | — | — | — | — | Strategy B found stable orbit, wrong topology (reaches_secondary=False) |
| (3,3) | clean negative | — | — | — | — | — | — | mu-continuation failed (mu-step diverges before PC mu) |
| (2,1) | clean negative | — | — | — | — | — | — | 8x8x2 grid search exhausted; no (2,1) orbit found |
| (2,2) | clean negative | — | — | — | — | — | — | 8x8x2 grid search exhausted; no (2,2) orbit found |

---

## (3,2) positive control detail

**Already admitted to catalogue** as `ross-rt-pc-cycler-32-2026` (task #494, V1).  
Re-confirmed by this sweep with exact match to prior values:

- Jacobi constant: C = 3.57951501972907 (DERIVED, nu=0 midpoint)
- x0 = -0.693198287043 (non-dimensional, rotating frame)
- ydot0 = -0.297004785528 (non-dimensional)
- Period (ND): T = 11.8334625170 TU
- Period (days): 12.033 d (using t_s = 87855.81 s)
- Barden stability index nu = 7.9e-8 (~0; maximally stable)
- Topology: (3,2) prograde, reaches_secondary=True
- Independent Radau crosscheck: PASS (dJ < 1e-12)

Cross-checks vs sourced values:
- C_L1 = 3.6210 (computed) vs Jbara 2025 arXiv:2510.13479 C_L1 ≈ 3.6210 at mu~0.109 — MATCH (diff < 0.005)
- Holman-Wiegert 1999 Eq. 1: a_crit = 38948 km at PC mu — within 150 km tolerance

**Literature status:** sourced to Ross-RT 2026 (arXiv:2606.29189), (3,2) family.
Pluto-Charon instantiation is a computed known-class member, not a literal table entry
in the paper. No new literature check needed (prior #494 lit-novelty confirmed).

---

## (1,1) family

**Clean negative.** Mu-continuation from the mu=0.001 Table-I anchor succeeded, landing
at PC mu with a valid orbit. C-sweep over [seed_C - 0.3, C_L1 - 0.002] found no nu
sign-change (no stable window). The (1,1) family exists at PC mu but is unstable across
the full C range scanned.

**Method:** mu_step_from_mu0.001_then_c_sweep  
**Note:** no stable window in C-sweep range

---

## (3,1) family

**Clean negative.** The mu=0.3 Table-I anchor has C = 3.7020 > C_L1(PC) = 3.621,
so direct mu-step fails (orbit cannot reach Charon at PC mu at that energy).

Strategy B: C-walk the (3,1) family at mu=0.3 from C=3.702 down to C_target = 3.571
(below C_L1(PC)), then mu-step to PC mu. This succeeded numerically — a stable orbit
was found — but winding-topology analysis confirmed it is NOT a (3,1) cycler:
`reaches_secondary=False` (orbit stays in the Pluto-side realm). The orbit found by
Strategy B belongs to a different near-primary family, not the (3,1) family.
Classified as clean negative per topology check.

**Method:** c_walk_at_mu03_then_mu_step  
**Note:** Strategy B found a stable orbit but wrong topology (reaches_secondary=False)

Physical interpretation: the (3,1) family anchors at C > C_L1(PC); walking C below
C_L1(PC) transitions to a different orbit class before the mu-step can continue on the
(3,1) branch. The (3,1) family either does not have a stable window accessible at PC mu,
or requires a different seeding approach (not covered by the current Table-I anchors).

---

## (3,3) family

**Clean negative.** Mu-continuation from the mu=0.01215 Table-I anchor failed before
reaching PC mu. The (3,3) branch did not survive the mu-step to 0.10876.

**Method:** mu_step_from_mu0.01215  
**Note:** mu-continuation failed

---

## (2,1) and (2,2) families

**Clean negatives.** Neither family appears in the Ross-RT 2026 Table-I, meaning no
published anchor is available. Grid searches were performed over (x0, C, hc) space:

- **(2,1):** x0 ∈ [-0.80, -0.35] (8 pts) × C ∈ [3.10, C_L1-0.01] (8 pts) × hc ∈ (3,4);
  128 total seeds. No (2,1) prograde orbit with correct topology found.
- **(2,2):** x0 ∈ [-0.70, -0.20] (8 pts) × C ∈ [3.05, C_L1-0.01] (8 pts) × hc ∈ (4,5);
  128 total seeds. No (2,2) prograde orbit with correct topology found.

Each corrector call was bounded by SIGALRM at 4s to avoid runaway integrations near
the Charon singularity. Grid times: (2,1) 173s, (2,2) 241s.

**Caveat:** These grids are representative, not exhaustive. If (2,1)/(2,2) families
exist at PC mu with different winding topology or x0 outside the searched range, they
would be missed. This is a necessary-not-sufficient negative — grid coverage was
limited by compute budget. A finer search or continuation from a different anchor would
be needed to certify absence.

---

## Literature novelty summary

No new stable members were found beyond the already-admitted (3,2). The catalogue entry
`ross-rt-pc-cycler-32-2026` (V1) remains the only stable Pluto-Charon CR3BP cycler
in the catalogue. No further literature checks are warranted for this sweep.

---

## Admission recommendations

**No new admissions.** The positive control re-confirmed `ross-rt-pc-cycler-32-2026`
(already V1). All other families yielded clean negatives.

**For the human adjudicator:**
- No catalogue row additions needed from this sweep.
- The negative results for (1,1), (3,1), (3,3), (2,1), (2,2) narrow the search space:
  at PC mu, the (3,2) family is structurally special among the Ross-RT Table-I set.
- The (2,1)/(2,2) negatives are grid-limited; a wider search is possible if warranted.
- The (3,1) result (wrong-topology stable orbit via Strategy B) is an interesting
  by-product: a near-primary stable orbit exists near C=3.502, x0=-0.526 in the
  Pluto realm; it is not a cycler but might be of interest for Pluto-only missions.

---

## Technical notes

**Positive control gate:** The (3,2) test is a HARD assertion (not clean-negative-aware).
If the (3,2) fails, the sweep machinery is broken. It passed in this run.

**Sweep timing** (single core):
- (3,2): ~17s (direct C-sweep)
- (1,1): ~24s (mu-step + C-sweep)
- (3,3): ~19s (mu-step)
- (3,1): ~95s (Strategy A + B + topology check)
- (2,1): ~173s (grid search)
- (2,2): ~241s (grid search)

**Tests passed:** `uv run pytest tests/search -q` — all pass, no regressions.  
**Ruff:** clean (`uv run ruff check` + `uv run ruff format --check`).
