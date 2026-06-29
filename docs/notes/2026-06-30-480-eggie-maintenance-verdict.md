# #480 — EGGIE maintenance-ΔV curve (real ephemeris, Approach A): RESULT

**Date:** 2026-06-30. Status: **DONE.** The scope doc's Approach A (quantify the per-cycle
maintenance ΔV of a feasible real-ephemeris EGGIE by chaining it forward) is executed. The
result **quantitatively reproduces the paper's qualitative maintenance claim**: a feasible
real-eph EGGIE is ballistic for exactly **2 cycles**, then requires large maintenance
impulses. Method-validated by the Liang Member D positive control (#223). No catalogue impact.

## Method (the validated lane, no steering)
- **Positive control FIRST** ([[feedback_verify_gauntlet_with_positive_control]]):
  `scripts/liang_member_d_run.py --n-cycles 10` reproduces the Liang CGCEC cycler at
  sub-nm/s defect for 9/10 cycles (one single-cycle optimizer hiccup) — the maintenance
  method (`nbody.jovian.chain_cycles`) is sound.
- **EGGIE lane** `scripts/eggie_maintenance_480.py`: Stage 1 scans 6 branch plans × 31
  epochs near the paper era (2020-09-15 +30 d) for a feasible ballistic EGGIE cycle-1 (V∞
  an OUTPUT, no steering); Stage 2 chains 10 cycles at the best seed. Seed:
  sequence `(Europa, Ganymede, Ganymede, Io, Europa)`, plan
  `((0,single),(1,high),(1,low),(1,high))`, departure 2020-09-22 (+7.00 d).

## The discovered feasible EGGIE member (cycle 1)
- Ballistic to **9.2e-10 m/s**, all flyby altitudes feasible (min 5,994 km).
- V∞ (km/s): Europa 9.38 → Ganymede **6.66** → Ganymede **6.66** → Io 7.35 → Europa 8.87.
  The two Ganymede flybys are at **EQUAL V∞** (the paper's defining ballistic-cycler
  property, Table-4 both 7.07). This member sits ~0.4 km/s BELOW Table-4 Ganymede (6.66 vs
  7.07) — consistent with the prior unguided discovery (the ballistic family is ~0.5 km/s
  below Table 4; `2026-06-30-480-eggie-realeph-unguided-discovery.md`). Not the exact
  Table-4 member (no such claim).

## The maintenance curve (the result)
| cycle | per-cycle ΔV (m/s) | cumulative (m/s) | min alt (km) |
|---|---|---|---|
| 1 | 0.0 | 0.0 | 5994 |
| 2 | 0.0 | 0.0 | 4895 |
| 3 | 167.7 | 167.7 | 50 (floor) |
| 4 | 759.0 | 926.8 | 50 |
| 5 | 656.0 | 1583 | 50 |
| 6 | 548.2 | 2131 | 50 |
| 7 | 359.9 | 2491 | 50 |
| 8 | 259.9 | 2751 | 50 |
| 9 | 286.2 | 3037 | 50 |
| 10 | 258.7 | 3296 | 50 |

**Ballistic for exactly 2 cycles, then large per-cycle impulses (≈170-760 m/s).** Cumulative
~3.3 km/s over 10 cycles.

**Robustness:** identical curve at retarget budgets of 2, 4, and 8 days (the breakdown is
invariant to the maintenance budget → a genuine geometric horizon, NOT an under-resourced
optimizer). All 6 branch plans converge to the same +7.00 d feasible seed.

## Comparison to the paper (honest)
- The paper prints **no EGGIE maintenance number** — only the general statement that the
  real-ephemeris solution's "ballistic repeatability... will only last for a few cycles,"
  and, for the high-fidelity EIGE specifically, "remains ballistic for two cycles, after
  which large impulses are required." Our patched-conic (level-2) EGGIE shows **exactly the
  same 2-cycle ballistic horizon then large impulses** — a qualitative + cross-cycler
  confirmation of the paper's central maintenance finding.
- This is a **NOVEL maintenance curve** (the EXPECTED values are our own computation — not a
  golden), not a reproduction of a printed number. The paper's one quantified maintenance
  figure (EIGE ~30 m/s/10) is for the 1-synodic EIGE, which the patched-conic lane cannot
  reproduce feasibly (the 1-rev B-plane wall, `2026-06-30-480-eige-realeph-maintenance-verdict.md`).
- Magnitude caveats: level-2 patched-conic (not the paper's level-3 high-fidelity); the
  per-cycle impulses are the cost the chained re-targeter pays once the geometry drifts off
  the ballistic resonance. The robust, model-independent finding is the SHAPE (2 ballistic
  cycles → breakdown), not the exact m/s.

## Standing
Approach A complete. The #480 maintenance-ΔV gap is now closed at level-2: the maintenance
method is validated (Liang #223) and the EGGIE maintenance behavior is characterized
(2-cycle ballistic horizon then large impulses, robust). No catalogue impact; no
exact-Table-4 / exact-paper-number claim. Driver `scripts/eggie_maintenance_480.py` is
reusable. Closing level-3 high-fidelity maintenance (the paper's exact regime) would need
the full B-plane NLP — out of scope (Approach C, "weeks, last resort").
