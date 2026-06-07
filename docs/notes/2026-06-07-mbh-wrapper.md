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
