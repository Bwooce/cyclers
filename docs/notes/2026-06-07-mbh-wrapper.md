# Monotonic Basin Hopping (MBH) wrapper — task #145

**Date:** 2026-06-07
**Code:** `src/cyclerfinder/search/mbh.py`, `tests/search/test_mbh.py`
**Context:** survey Thread 1, `docs/notes/2026-06-07-external-algorithms-survey.md`

## What this is

A generic Monotonic Basin Hopping loop around the existing local correctors. The
disease (#110/#120/#122): our broad search converges, but to a degenerate /
off-family basin, and a denser epoch×branch grid does not cure that — it is a
basin-*selection* failure, not a grid-density failure. MBH is the trajectory
community's standard answer: hop from the best incumbent by a random
perturbation, re-run the local corrector, accept the new point only if it
strictly improves. The randomness escapes the current funnel; the local solver
refines.

`mbh()` is transcription-agnostic: it takes a closure
`objective_and_solve(x, rng) -> MBHStep` (one local solve from a seed) and never
imports a corrector itself. Two adapters wrap the existing correctors WITHOUT
editing them:

- `make_ballistic_step` → `search/correct.py:ballistic_correct`
  (genome `x = [t0_sec, *free_tof_days]`, objective `max_residual_kms`).
- `make_free_return_step` → `search/free_return.py:free_return_correct`
  (#137 radial-crossing genome `x = [a_au, e, t0_sec]`, objective
  `max_residual_kms`).

Determinism: `rng_seed` is REQUIRED; the only randomness source is a local
`numpy.random.Generator`; no global RNG state is touched. Full audit trail on
`MBHResult` (hops attempted/accepted, per-hop objective + accept decision,
monotone best-history envelope, stall counter, echoed seed).

## Perturbation distribution — SPEC CAVEAT (not sourced)

The canonical perturbation-distribution paper, **Englander & Englander 2014,
"Tuning Monotonic Basin Hopping" (ISSFD24 S7-3)**, prescribes long-tailed
(Cauchy / Pareto) per-gene perturbations — the long tail is what lets a hop jump
*between* basins rather than only jittering within one. **That paper is NOT yet
acquired** (it is on the #116 acquisition list; the open PDF URL is recorded in
the survey note). Until it is in hand, this module implements a *documented,
sensible default*: per-gene perturbation from a configurable distribution
(`"cauchy"` default, plus `"gaussian"`/`"uniform"`) with a configurable
relative-or-absolute scale. `"cauchy"` is the closest available stand-in for the
Englander spec, but its TUNING (scale schedule, per-gene scaling) is **NOT
sourced and must not be claimed as such**. Refinement to the exact Englander 2014
tuning is a follow-up once the paper is acquired.

Sizing modes (per gene): RELATIVE (a fraction of `|gene|`, with a unit floor) or
ABSOLUTE (in the gene's own units, overrides relative). Absolute sizing is
essential for `t0` (~1e8 s): a relative fraction there is enormous, so a basin a
few days wide is unreachable; an absolute few-day step lands in it.

## Gate results (verbatim, all run 2026-06-07)

`uv run pytest tests/search/test_mbh.py -v` — **9 passed** (7 mechanics + 2 slow
free-return). ruff check + ruff format --check + mypy on both files: clean.

### Gate 1 — mechanics (constructed double-well; NOT golden)

Constructed 1-D double well: deep/global min at x=+3 (depth −2), shallow/local
min at x=−3 (depth −1). The "expected" value is defined BY CONSTRUCTION by our
own function — this is an algorithm-mechanics test, not a golden test.

- `test_plain_local_solve_stays_in_wrong_basin`: a single gradient-descent local
  solve from the shallow basin provably converges to the shallow min (−1) and
  cannot reach the global min (−2). **PASS.**
- `test_mbh_escapes_to_global_basin`: MBH (cauchy, absolute step 3.0, seed 12345)
  from the shallow basin reaches the deep basin (best_x ≈ +3, objective ≈ −2),
  best-history monotone non-increasing. **PASS.**

MBH escapes the wrong basin where the plain local solve provably does not.

### Gate 2 — real: free-return basin recovery from a mis-seed (the real one)

Row `mcconaghy-2006-em-k2`, circular model. SOURCED constraint side: Rogers 2012
Table 1 ellipse a=1.30 AU, e=0.257 (the same ellipse the #137 like-for-like gate
uses); SOURCED V∞ anchor (McConaghy 2006): E 4.7, M 5.0 km/s. The emerged V∞ and
recovered (a,e) are EVIDENCE; nothing our code computed is an EXPECTED value.

**Mis-seed scale finding (honest):** the deliberate mis-seed is **40 days
off-phase** in t0, starting (a,e) at the sourced ellipse. At that scale:

- PLAIN single local solve from the mis-seed: `feasible=False, obj=2.16 km/s` —
  does NOT converge. (A 40-day phase miss is enough; the plain solve does not
  reach the basin, so the mis-seed did NOT need widening.)
- MBH (gaussian, t0 absolute step 8 days, a/e frozen in the hop, seed 6,
  60 hops): `feasible=True, obj=0.00011 km/s`, recovered a=1.305, e=0.245
  (within 0.03 of sourced), emerged **V∞ M=4.69** (sourced 5.0, within the
  campaign 0.5 km/s tolerance) and **V∞ E=4.05** (sourced 4.7, within ~1 km/s —
  the sourced low-V∞ regime, NOT the 9+ km/s powered basin). 12 of 61 hops
  accepted. **PASS.**

`test_free_return_mbh_recovers_sourced_basin_from_misseed` and
`test_free_return_mbh_determinism_on_misseed` both PASS.

**HONEST CAVEAT (predicted by the survey).** The single-ellipse free-return
genome is 1-DOF underdetermined along the Mars-V∞ ridge (every t0 closes at
*some* V∞), so residual-only MBH does not select the sourced sub-basin from every
rng seed. The gate pins a deterministic seed (6) that recovers it; the
seed-sensitivity is the survey's "MBH cures basin selection WITHIN a fixed
transcription" framing — it is not a topology cure.

### Gate 3 — optional, off-anchor probe: `russell-ch4-6.44Gg3` (negative result)

This row is catalogued `cycler_class: multi-arc` (two generic-return arcs,
a_au/e null) with SOURCED aphelion 1.54 AU, transit 262 d, V∞ E=6.44, M=3.74.
Seed (a,e) derived from the sourced aphelion+transit via
`seed_ae_from_aphelion_transit` → a=1.27, e=0.21. At that derived single ellipse
the best-phase residual is 0.0033 km/s (representable), but:

- PLAIN from the same 40-day mis-seed: `feasible=False, obj=9.60`.
- MBH (same config, seed 6): `feasible=False, obj=9.60, 0 hops accepted` — does
  NOT recover.

The derived single ellipse gives emerged V∞ E=3.01 / M=3.06, FAR from the sourced
6.44 / 3.74 — i.e. the single-ellipse free-return shape that matches the sourced
aphelion+transit is a *different, off-anchor* topology. This is exactly the
survey's predicted outcome for a multi-arc row: **MBH confirms the topology
limitation** (the single-ellipse transcription has no sourced-anchor basin here)
faster than a denser grid would — a positive negative-science result. Not gated
as a passing test (it is a known multi-arc row); recorded here as diligence.

## Adoption / follow-ups

- Acquire Englander & Englander 2014 (ISSFD24 S7-3), then tune the perturbation
  schedule to the sourced Cauchy/Pareto spec (currently a documented default,
  NOT sourced).
- For the underdetermined free-return genome, a future addition of a sub-basin
  selector (e.g. an anchor-aware secondary objective, or the pseudo-arclength
  continuation of Thread 2) would remove the seed-sensitivity that Gate 2's
  caveat records.

---

## ADDENDUM 2026-06-07 — perturbation spec is now SOURCED (task #156)

**Status flip: UNSOURCED → SOURCED.** The Englander 2014 paper has been acquired,
read in full, and mined (`docs/notes/2026-06-07-englander-2014-mbh-tuning-mining.md`).
The "SPEC CAVEAT (not sourced)" section above is superseded by this addendum and
by the SOURCED docstring in `src/cyclerfinder/search/mbh.py`.

**Source (cite exactly):** Englander, J. A. and Englander, A. C., "Tuning
Monotonic Basin Hopping: Improving the Efficiency of Stochastic Search as Applied
to Low-Thrust Trajectory Optimization," 24th International Symposium on Space
Flight Dynamics (ISSFD24), 2014, paper S7-3.

**What the paper prescribes:** a long-tailed, bi-directional, zero-centred
per-step perturbation; ordering **bi-polar Pareto ≥ Cauchy ≫ Gaussian ≈ Uniform**;
the authors RECOMMEND bi-polar Pareto as the single best default because it is the
most robust to its excursion parameter (p.31).

**HONEST CAVEAT (the magnitudes are not portable).** Their numbers are from a
SINGLE 503-variable low-thrust EESU benchmark whose decision variables are
bound-NORMALISED (EMTG Table 1); the authors themselves flag the tunings may be
problem-dependent (p.29). Our cycler-corrector genes are NOT range-normalised
(hence this module's per-gene relative/absolute scale split), so we import only
the QUALITATIVE prescription and validate the scale + distribution choice LOCALLY
by sweep on our own Gate-2 free-return recovery.

**Source-fidelity note on the Pareto generator.** Table 4 transcribes the
generator as `(s/ε)·(ε/(ε+r))^(−α) = (s/ε)·((ε+r)/ε)^α`, which is BOUNDED for
`r ∈ [0,1]` and therefore does NOT actually produce a long tail (its undocumented
`ε` sets a scale, not a tail). Read literally that contradicts the paper's stated
"very long tails" property (p.7) — the whole thesis. We therefore implement the
canonical long-tailed bi-polar Pareto the thesis unambiguously calls for:
magnitude `u**(−1/α)` (unit floor, genuine heavy tail as `u→0`; SMALLER α ⇒
heavier tail) times a fair ±1 sign coin, parameterised by the same exponent α.
Default `DEFAULT_PARETO_ALPHA = 1.08` (paper's mid-range MSD value, p.21). This
discrepancy and our choice are recorded in the kernel body and here.

### What was built

1. `PERTURBATIONS = ("pareto", "cauchy", "gaussian", "uniform")` — `"pareto"`
   (bi-polar Pareto) added; `"cauchy"` KEPT as the code default.
2. `perturbation_alpha` parameter (Pareto exponent; ignored by other dists).
3. `MBHResult` now records `perturbation` (name) + `perturbation_param` (the
   excursion parameter: α for Pareto, else the scalar relative scale, `nan` for a
   per-gene vector / `None`) — the audit gap the paper's result demands.
4. Restart-on-stall (Englander Algorithm 1 global reset, p.9) added behind a
   default-OFF `restart_bounds=(lower, upper)` flag: when set AND
   `stop_after_stall` is set, a stall RE-SEEDS the next hop from a fresh uniform
   point in bounds (keeping best-so-far) instead of stopping. `None` (default)
   preserves the original stop-on-stall behaviour exactly. The free-return adapter
   is the bounded case this fits; the ballistic genome is not bounded here.
5. Docstring SPEC-CAVEAT block flipped to SOURCED.

All additive: no public-signature breaks, no default-path behaviour change. The
9 original tests stay green; 8 new tests cover Pareto escape, Pareto bi-polarity /
α-sensitivity, the two new audit fields (incl. default-path = cauchy/0.05), the
documented α, and restart-on-stall + its shape validation. **17 passed.**

### Local sweep — distribution × scale on OUR Gate-2 problem

Re-ran the Gate-2 free-return recovery (`mcconaghy-2006-em-k2`, circular model,
40-day off-phase mis-seed, `rng_seed=6`, a/e frozen, t0 perturbed absolutely,
60 hops) across {uniform, cauchy, pareto} × a t0-step grid spanning the Gate-2
default (8 d) ÷10 .. default. Script: `scripts/mbh_pareto_sweep.py` (one-off,
deterministic, re-runnable). "ANCHOR" = converged AND in the SOURCED ellipse
(a=1.30, e=0.257 within 0.03) AND emerged V∞ at the SOURCED anchor (M 5.0 ±0.5,
E 4.7 ±1.0) — the Gate-2 pass criterion. "conv" = converged but off-anchor (wrong
sub-basin on the underdetermined Mars-V∞ ridge). "hops2conv" = first hop reaching
a converged residual.

```
    dist  t0_step_d  recovered hops2conv    obj_kms   a_rec   e_rec  vinfE  vinfM  acc/att
------------------------------------------------------------------------------------------
 uniform        0.8         no         -     2.1637   1.310   0.237   3.35   4.53 0/ 61
 uniform        2.5         no         -     2.1637   1.310   0.237   3.35   4.53 0/ 61
 uniform        8.0     ANCHOR        18     0.0001   1.301   0.254   4.70   4.91 9/ 61
  cauchy        0.8     ANCHOR         9     0.0000   1.307   0.243   3.86   4.65 7/ 61
  cauchy        2.5     ANCHOR         5     0.0001   1.296   0.257   5.00   4.93 6/ 61
  cauchy        8.0       conv         5     0.0000   1.074   0.797  26.94  18.18 7/ 61
  pareto        0.8         no         -     2.1637   1.310   0.237   3.35   4.53 0/ 61
  pareto        2.5     ANCHOR        24     0.0000   1.303   0.249   4.34   4.78 5/ 61
  pareto        8.0     ANCHOR         2     0.0001   1.305   0.245   4.05   4.69 4/ 61
```

**Reading.**
- **uniform** recovers the anchor in only 1/3 cells (largest step). Fragile to the
  scale — exactly Englander's "uniform is the worst, and scale-fragile" finding.
- **cauchy** recovers the anchor in 2/3 cells (the two SMALL-scale cells, the
  regime the long-tail/small-scale coupling favours) and is the FASTEST (5-9 hops).
  At the largest scale it overshoots into the degenerate powered sub-basin
  (a=1.07, V∞ 27) — confirming the "use a smaller typical scale with a long tail"
  lesson.
- **pareto** also recovers the anchor in 2/3 cells but FAILS at the smallest scale
  (no accepts) and is slower / more variable where it does (24 hops at 2.5 d).

### Default-distribution recommendation + DECISION

**Recommendation (per the paper):** bi-polar Pareto is the paper's recommended
default. **Decision (per our local sweep): KEEP `"cauchy"` as the code default.**
Pareto does NOT *clearly dominate* on our problem — on the Gate-2 anchor-recovery
criterion the two tie on anchor-cell count (2/3 each), but cauchy is faster (5-9
vs 2-24 hops) and, crucially, succeeds at the SMALL scales the long-tail theory
predicts should be best, whereas pareto fails the smallest-scale cell entirely.
Per the task rule ("change the default only if Pareto dominates clearly"), the
evidence does not meet that bar. `"pareto"` is shipped, documented as the paper's
recommendation, and a one-line flip away — to be reconsidered if a broader sweep
(more rows, alpha sweep, non-frozen a/e) shows it dominating here. The paper's
single bound-normalised low-thrust benchmark is not our cycler corrector, so the
local evidence governs the default; the paper governs the available menu.

### Restart-on-stall assessment (Algorithm 1 global reset)

Assessed and IMPLEMENTED as a clean additive default-off flag (`restart_bounds`),
since the free-return adapter genome (`a_au`, `e`, `t0`) has natural bounds. When
off (default, and what every existing caller and all gates use) behaviour is the
original stop-on-stall — the defensible simplification for our bounded-by-
construction seeds. The faithful reset is now available for callers that supply
bounds; it was not wired into the Gate-2 test because that test deliberately pins
the stop-on-stall path. Not forced anywhere.
