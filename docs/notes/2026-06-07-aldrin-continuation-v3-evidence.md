# Aldrin continuation V3 evidence — horizon-TCM build (task #161)

**Type:** evidence build + results. Builds the V3 evidence path the #160
assessment defined as missing, runs it, records the numbers, and states whether
the V3 evidence chain is complete. **No catalogue writeback is performed** (that
is the user's call — see "Conclusion").

**What #160 said was missing.** The #158 continuation driver walked the Aldrin
`6.399G1` free-return from circular-coplanar to the true ephemeris (DE440),
closing **ballistically** to a 0.00158 km/s residual at `t0` 26 d from Russell's
sourced Aug-6-2003 window. The #160 promotion assessment
(`docs/notes/2026-06-07-aldrin-promotion-assessment.md` §3) held at V2-powered
because the continuation produces a *single closing arc*, while spec §14 **V3**
requires, verbatim:

> "phase-matched to a real launch window; ephemeris-mode horizon TCM over 3-5
> laps (~20-30 yr) bounded and within ΔV budget" — spec §14

The phase-match half was already present (continued `t0` 26 d from the window).
This task built the **missing half**: the ephemeris-mode horizon TCM over 3-5
laps, chained per-cycle in-family, anchored at the continued ballistic solution.

**Build:** `tests/verify/test_aldrin_continuation_v3.py` (new, `@pytest.mark.slow`).

---

## 1. The build

Two slow tests:

1. **`test_continued_aldrin_start_point_is_ballistic_in_window`** — reproduces the
   continuation driver's winning `nstep=3` rung for the Aldrin `6.399G1` golden
   directly (`ladder=(3,)`, rng-free and deterministic, ~0.3 s — far cheaper than
   the full ladder, *same converged solution*; `243`/other rungs recorded as
   `skipped`, no silent cap). This is the V3 **start point** and confirms it
   reproduces the #158 documented numbers: residual **0.00158 km/s**, converged
   `(a, e) = (1.5249 AU, 0.3616)`, launched **2003-07-11** (26 d from Aug-6-2003,
   inside the 200-day gate), emerged V∞ E=6.08 / M=8.88 km/s.

2. **`test_continued_aldrin_v3_horizon_tcm_bounded_and_within_bar`** — anchors the
   #134 horizon-TCM machinery (`optimise_aldrin_maintenance_dv`) at the continued
   solution's phase-matched window and chains it over 5 laps, re-phasing the
   priority date one cycler period (2.135 yr) per lap, re-solving the in-family
   maintenance ΔV on DE440 each lap. Records the per-lap chained horizon TCM.

**Model fidelity (honesty).** This is **ephemeris-positions two-body + DE440
propagation = V3-class fidelity**, NOT n-body. Russell's Chapter-5 force model
(and the continuation endpoint) is patched-conic two-body legs between
true-ephemeris planet positions; the maintenance chain re-solves Lambert /
free-return legs on DE440 planet states. This is exactly spec §14 V3's "astropy
backend / ephemeris realisation" class — the precursor (not the substitute) for
the planned REBOUND n-body harness (V4-class).

---

## 2. The numbers (verbatim from the run)

Continued start: **2003-07-11** (the phase-matched window anchor).
Cycler repeat period for re-phasing: **2.135 yr** (= 1 E-M synodic, the #134 value).

| lap | priority date | converged | a (AU) | e | maint. ΔV (km/s) | turn req° | turn max° | feasible |
|----:|---------------|:---------:|-------:|------:|-----------------:|----------:|----------:|:--------:|
| 0 | 2003-07-11 | yes | 1.5878 | 0.3932 | **2.9138** | 93.0 | 68.5 | **no** (powered) |
| 1 | 2005-08-28 | yes | 1.5877 | 0.3777 | **0.0000** | 59.1 | 81.1 | yes (ballistic) |
| 2 | 2007-10-17 | yes | 1.5877 | 0.3775 | **0.0000** | 58.5 | 81.3 | yes (ballistic) |
| 3 | 2009-12-05 | yes | 1.5876 | 0.3773 | **0.0000** | 57.8 | 81.5 | yes (ballistic) |
| 4 | 2012-01-24 | yes | 1.5876 | 0.3771 | **0.0000** | 57.2 | 81.7 | yes (ballistic) |

**Chained horizon TCM (continued-start):**

| laps | horizon TCM (km/s) | #134 powered-Aldrin comparison |
|-----:|-------------------:|--------------------------------|
| 3 | **2.9138** | 8.51 |
| 4 | **2.9138** | 11.19 |
| 5 | **2.9138** | 13.79 |

Only **lap 0** contributes; laps 1-4 are ballistically feasible (dv = 0).

---

## 3. The gate verdicts

The gate is qualitative, sourced-floor discipline (NOT an invented absolute
golden), per the #160 definition:

- **(a) BOUNDED — PASS.** Every lap's per-cycle re-solve **converges in-family**
  (a ≈ 1.588 AU, e ≈ 0.377-0.393 — within the sourced Aldrin anchors). No lap
  slides to the degenerate high-energy basin (e→0.95, V∞~38 km/s) that the
  codebase's #114 finding warns about. The horizon TCM is finite and does not
  grow lap-over-lap: it is **flat at 2.9138 km/s** (it does not even accumulate,
  because laps 1-4 add zero). No divergent outlier; horizon TCM < n_laps × 3.0
  for all of 3/4/5 laps.

- **(b) WITHIN BAR — PASS.** Every per-cycle ΔV is **under the engineering
  plausibility bar** `MAINTENANCE_DV_CONVENTION_KMS = 3.0 km/s/cycle`
  (`verify/plausibility.py` — a CONVENTION bar, not a sourced Aldrin budget;
  no Aldrin maintenance-ΔV magnitude is published, McConaghy 2002 defers it,
  catalogue `data_gaps`). Max per-cycle = 2.9138 km/s < 3.0. The bound holds at
  the per-cycle granularity, not just the sum.

---

## 4. The interesting finding (why the continued-start horizon is so much lower)

The task flagged the open question: *the continued solution starts
ballistic-at-epoch, so how does the budget compare to #134's powered-Aldrin
horizon (8.51 / 11.19 / 13.79 km/s over 3/4/5)?* The answer is striking and
honest:

**The continued-start horizon TCM is 2.9138 km/s flat over 3/4/5 laps — ~3-5×
lower than #134's powered chain, because only the first lap is powered.**

- **Lap 0** (the exact continued window, e = 0.3932 = the sourced anchor) lands on
  the **powered** Aldrin member: return-flyby turn required 93.0° > 68.5°
  achievable → not ballistically feasible → maintenance ΔV = 2.9138 km/s (the
  identical #134 outbound number — this lap is the #134 per-cycle solve).
- **Laps 1-4** (re-phased one synodic period forward) land on a **slightly
  lower-eccentricity in-family member** (e ≈ 0.377): its return-flyby turn
  required (~57-59°) is **below** the achievable max (~81°), so the orbit repeats
  **ballistically** and the in-family maintenance solve returns **dv = 0.0**.

This is the same "ballistic dv≈0 in-family neighbour" phenomenon the codebase
documents for the inbound twin (`test_aldrin_v2_v3_campaign.py`,
`validate.py` inbound V1 evidence) — here it appears *favourably* across the later
laps. Two honest readings:

1. **Optimistic:** the maintenance optimiser, re-phased forward from the continued
   window, finds that the Aldrin family is *ballistically maintainable at most
   laps* on the real ephemeris — the horizon TCM is dominated by a single powered
   lap, far inside any plausible budget.
2. **Cautious (the one I flag):** laps 1-4 converge on a **different in-family
   member** (e ≈ 0.377) than lap 0's exact-anchor (e = 0.393). The per-lap solve
   is free to choose the lowest-ΔV in-family solution at each window, so the chain
   is **not** a single fixed cycler re-flown — it is a per-lap *re-optimisation*
   that happens to find ballistic members. This is the **honest caveat**: the
   2.9138 km/s flat horizon is a *lower bound on the maintainable Aldrin family
   across these windows*, not proof that one fixed continued cycler flies 5 laps
   for 2.9138 km/s. (Both #134 and this build share this per-lap-reoptimisation
   structure; #134's chain re-solved to the powered member each lap, this one
   finds ballistic members at laps 1-4.)

Either reading satisfies the §14 V3 words: the horizon TCM **is** bounded (flat,
finite, in-family, no divergence) and **within budget** (< the 3.0 km/s/cycle
engineering bar, with no published Aldrin budget to better it).

---

## 5. The V3 phase-match half (already present — cited, not re-derived)

Spec §14 V3's "phase-matched to a real launch window" is satisfied by the
continued solution itself: `t0 = 2003-07-11` is **26 d** from Russell's sourced
Aug-6-2003 window (well inside the test's 200-day gate). This is the
`test_continued_aldrin_start_point_is_ballistic_in_window` evidence and was
already established by #158; it is reproduced here as the horizon chain's anchor.

---

## 6. Conclusion — is the V3 evidence chain complete?

**The V3 evidence chain is COMPLETE pending user-approved writeback (no writeback
performed here).** Both V3 halves are now on disk as passing mechanical evidence:

- **phase-match** to a real launch window — continued `t0` 26 d from the sourced
  Aug-6-2003 window (`test_continued_aldrin_start_point_is_ballistic_in_window`);
- **ephemeris-mode horizon TCM over 3-5 laps, bounded and within budget** —
  2.9138 km/s flat over 3/4/5 laps, every lap converging in-family under the
  3.0 km/s/cycle engineering bar
  (`test_continued_aldrin_v3_horizon_tcm_bounded_and_within_bar`).

This clears spec §14 V3's words **in the V3-class (ephemeris-positions two-body +
DE440) fidelity** — explicitly NOT n-body (that is V4 / the REBOUND harness).

**Caveats the user must weigh before any writeback (I do NOT promote):**

1. **Per-lap re-optimisation, not one fixed cycler.** Laps 1-4 converge on a
   slightly different in-family member (e ≈ 0.377) than the continued start's
   exact anchor (e = 0.393). The 2.9138 km/s flat horizon is a *lower bound on the
   maintainable Aldrin family across these windows*, not a single continued cycler
   re-flown for 5 laps. This is the same per-lap structure #134 used and is the
   honest framing for the number. A stricter "one fixed cycler, retargeted only"
   horizon would re-solve the *powered* member each lap and recover ~#134's
   8.51/11.19/13.79 — which is *also* bounded and *also* < n × 3.0/cycle, so the
   V3 verdict is robust to which reading is taken, but the magnitude is not.

2. **The 3.0 km/s bar is a CONVENTION, not a sourced Aldrin budget.** No Aldrin
   maintenance-ΔV magnitude is published (McConaghy 2002 defers it; catalogue
   `data_gaps`). "Within ΔV budget" is satisfied against the project's engineering
   plausibility bar, exactly as the #134 note frames the absence of a sourced
   budget. A writeback's `_LEVEL_EVIDENCE` entry must record this honestly (the
   budget is the engineering convention, not a source-attested figure).

3. **Model scope.** The evidence is **true-ephemeris (patched-conic two-body legs,
   DE440)**, the V3-class realisation of the Aldrin family — explicitly distinct
   from the catalogue row's circular-coplanar *defining* model (#160 §1) and
   explicitly **not** n-body. Any writeback must record the model scope as
   true-ephemeris V3-class, same convention as the V1 scoping.

**If/when the user approves a writeback**, the row is `aldrin-classic-em-k1-
outbound` (the only Aldrin catalogue home, #160 §1), the level is **V3**, and the
`_LEVEL_EVIDENCE[("aldrin-classic-em-k1-outbound", "V3")]` entry would cite:
- sourced EXPECTED: Russell 2004 Table 5.5 (p.178) `6.399G1 (#1)` Aldrin, 0 m/s,
  Aug-6-2003 window;
- in-repo EVIDENCE: `tests/verify/test_aldrin_continuation_v3.py` (both passing
  slow tests — phase-match + bounded/within-bar horizon TCM);
- model scope: **true-ephemeris (patched-conic two-body legs, DE440) = V3-class,
  NOT n-body**; budget bar = engineering convention (no sourced Aldrin ΔV).

Until the user approves, the row stays at its earned **V2-powered** (golden
discipline: when in doubt, hold). The mechanical evidence for V3 now exists on
disk — the decision to cite it is the user's.

---

## 7. Run record

- Both slow tests **PASS** (~53 s wall total, well inside the ≤30-min cap; no lap
  cap applied — all 5 laps ran).
- ruff check + ruff format: clean. mypy whole-tree (`src tests`): clean for this
  file (a concurrent sibling's WIP `tests/data/test_moontour_retag.py` has
  unrelated mypy errors; this file commits `--no-verify` with the per-file checks
  verified individually, per the concurrent-agent rule).
- No `search/`, `core/`, `nbody/`, `data/catalogue.yaml`, or `spec.md` touched.
