# #464 — Sims-Flanagan low-thrust releg branch: VERDICT

**Date:** 2026-06-26
**Task:** #464 (the deferred Task 7 of the #449 plan — the SF low-thrust releg backend).
**Plan:** `docs/superpowers/plans/2026-06-25-449-lowthrust-dsm-releg-genome-plan.md` (Task 7)
**Design:** `docs/superpowers/specs/2026-06-25-449-lowthrust-dsm-releg-genome-design-draft.md` (§2 SF backend, §6 golden, §7 risks)

---

## 1. One-line verdict

**Low-thrust is real and slightly cheaper, but it does NOT bring the Galilean
moon-tour closure in-band.** The Sims-Flanagan low-thrust releg narrows the
Galilean positive-control cycle from the DSM branch's **13.18 km/s/cycle** to
**12.03 km/s/cycle** (~9% cheaper) — still **~3.4× above** the 3.5 km/s/cycle
powered dv-band ceiling. The moon-tour closure is **powered-EMPTY even with
low-thrust**. This is a characterized, honest NEGATIVE, not an in-band hit.

No tour reached in-band; nothing is flagged for human gauntlet review.

---

## 2. What was built (Task 7)

`LowThrustReleg`, a third `Releg` backend in
`src/cyclerfinder/search/releg_solver.py`, behind the same protocol as
`DsmReleg` — swappable into the `releg_moontour` driver with no driver rewrite
(the driver's powered-backend selection now keys off "not `BallisticReleg`", so
SF swaps in for DSM automatically).

* It wraps the existing #309 Sims-Flanagan leg solver
  (`core.sims_flanagan.SimsFlanaganLeg` + `search.lowthrust.solve_leg_min_dv`) —
  no new optimiser was written.
* Boundary states pin the departure V∞ magnitude to `vinf_depart_mag` and
  retarget the arrival V∞ magnitude to `vinf_target_in` along the
  ballistic-Lambert seed directions (the driver's continuity-by-construction
  contract — identical to `DsmReleg`).
* The deliverable ΔV is distributed across the SF thrust train (the source of
  its efficiency advantage over the single-impulse DSM).
* Thrust capability is auto-scaled from the V∞ delta so the leg has the reach to
  deliver the retarget (a too-weak train cannot close the boundary defect).

Committed in git history (see §6 on the concurrent-agent commit-attribution
note).

---

## 3. Validation status

| Test | Status | What it asserts |
|---|---|---|
| `test_lowthrust_releg_zero_retarget_matches_ballistic` | PASS | Coplanar/zero-retarget regression: the SF leg's all-zero thrust limit closes at ΔV ≈ 0 and its V∞ chain matches the ballistic lowest-energy branch (the SF zero-thrust limit IS the ballistic leg, like `BallisticReleg`). |
| `test_lowthrust_releg_brackets_dsm` | PASS | Bracket golden (sourced): on a Ganymede→Europa retarget (T=4 km/s) the SF ΔV (1.97 km/s) brackets `[vilm_dv_min floor (1.71 km/s), DsmReleg cost (4.20 km/s)]`. |

**On the golden's reach (honest).** The SF leg model has **no clean state-level
literature anchor** — its own `lowthrust.py` docstring says so, and the
Vasile-Campagnola 2009 DFET transcription (the only digested low-thrust
outer-moon tour) "DOES NOT MAP" to our Sims-Flanagan leg (digest
`2026-06-07-vasile-campagnola-dfet-method-mining.md` §2.6: a different
transcription — polynomial finite elements vs midpoint-impulse). So the SF
golden is **bracket-only**, and the bracket's lower edge is the SAME sourced
Campagnola-Russell Part-1 Table 1 VILM floor (1.71 km/s, Ganymede-Europa ΔV_min)
the tight DSM golden uses — never a number SF itself computed
(`feedback_golden_tests_sourced_only`). The upper edge is the single-impulse DSM
cost (distributing ΔV across a thrust train is strictly more efficient than one
impulse). This is the honest, falsifiable, non-circular bar; we assert the
bracket and say plainly it is a bracket, not an equality.

---

## 4. THE KEY EXPERIMENT (the decisive result)

`LowThrustReleg` run on the SAME Galilean positive control the DSM branch closed,
through the same `releg_moontour.close_powered_cycle` driver (common flyby V∞
target swept over the probe band, continuity-by-construction, total ΔV summed):

**Skeleton:** Io → Europa → Ganymede → Io, per-leg ToF = 1.5× geometric-mean of
the orbital periods (3.76 / 7.56 / 5.34 days). Both legs link at V∞ = 4 km/s.

| Backend | Continuity residual | Per-leg ΔV (km/s) | **Total ΔV/cycle** | In-band (≤ 3.5)? |
|---|---|---|---|---|
| DSM (single impulse) | 0.000 km/s | 4.99 / 2.58 / 5.61 | **13.18 km/s** | **NO** |
| **Low-thrust (SF)** | 0.000 km/s | 4.64 / 2.21 / 5.19 | **12.03 km/s** | **NO** |

**Powered dv-band:** `[300 m/s, 3.5 km/s/cycle]` (`verify/dv_band_acceptance.py`).

**Finding.** Low-thrust narrows the closure by ~1.15 km/s (~9%), from 13.18 to
**12.03 km/s/cycle** — but it stays **~3.4× above** the 3.5 km/s/cycle ceiling.
The Galilean moon-tour is powered-EMPTY even with the more efficient distributed
thrust.

**Why low-thrust does not rescue it.** The cost is dominated not by
impulse-vs-distributed inefficiency but by the *irreducible magnitude of the V∞
retarget itself*. At a common flyby target of V∞ = 4 km/s, the individual legs
naturally want arrival V∞ of 5-9 km/s (the high-V∞ Galilean basin), so each leg
must shed 2-5 km/s of V∞ to enforce continuity. Distributing that shed across a
thrust train saves ~0.3-0.4 km/s per leg — a real efficiency gain, but an order
of magnitude smaller than the retarget magnitude it cannot avoid. The wall is the
V∞-basin geometry, which neither transcription escapes.

This corroborates the #449 verdict (commit `071064e`, the powered-empty finding)
from a second, independent, more-efficient transcription: the Galilean
moon-tour closure is a *structural* powered-empty, robust across the DSM and the
SF low-thrust branches.

---

## 5. Honesty / kill-criterion (design §7)

Per design §7's kill-criterion ("KILL the 'powered helps' claim for a region if
the only closures need ΔV above the 3.5 km/s/cycle powered ceiling — then it is a
*powered* empty"): the Galilean control closes only at 12.03 km/s/cycle, well
above the ceiling, **under BOTH** the DSM and the SF low-thrust backends. The
region is therefore a *powered empty* (now confirmed by low-thrust), suitable for
a capability-subsumption re-stamp of the registry — but the actual re-stamp
writeback is a campaign-issue action, not a build action (the #449 plan ships the
capability; the campaign spends it).

No fudging of the band or the golden: the measured ΔV is reported as-is, above
the band, and the golden is asserted only as the bracket the source supports.

---

## 6. Concurrent-agent commit-attribution note (transparency)

During the Task 7 commit a sibling agent committed concurrently and its
broadly-staged commit `a263304` ("genome: #333 augmented QP-tori ...") **swept in
the three Task 7 files** (`releg_solver.py`, `releg_moontour.py`,
`test_releg_solver.py`) alongside its own — a shared-tree staging race (the
hazard recorded in `feedback_concurrent_agent_git_rules`). The Task 7 CODE is
correctly committed and all 7 releg-solver tests pass at HEAD; only the
commit-message attribution is tangled. History was NOT rewritten (the memory rule
forbids resetting past a sibling's commit on the shared `main`), so the tangle is
left in place rather than risk the sibling's #333 work. This verdict doc is the
clean #464-attributed record of the Task 7 deliverable.

---

## 7. Definition-of-done check (Task 7)

- [x] `LowThrustReleg` backend behind the `Releg` protocol, wrapping the #309 SF
      solver (no new optimiser).
- [x] Wired as a swappable backend into `releg_moontour` (selector keys off
      "not BallisticReleg").
- [x] Coplanar-limit / regression test (zero-retarget ≈ ballistic, ΔV ≈ 0).
- [x] Bracket golden (SF ΔV brackets the VILM floor and the DSM cost; sourced
      lower edge, honestly a bracket because the SF leg model has no clean
      state-level anchor).
- [x] KEY EXPERIMENT run: Galilean control low-thrust ΔV/cycle = **12.03 km/s**,
      vs 13.18 km/s DSM, vs the 3.5 km/s band — **NOT in-band** (honest negative).
- [x] No tour reached in-band; nothing flagged for gauntlet review.
- [x] Verdict written (this doc).
