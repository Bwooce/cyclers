# Phase C verdict — Jones n-body flyby-propagation shooter (#133)

Date: 2026-06-07. Closes the build phase of #133 (Jones Phase 3) and the
shooter rung of the n-body harness (#136). Written by the main session from
the bounded gate's captured output after three agent sessions on this lane
(two killed by session limits; incident record below).

## What was built

Commits `259a0d3`, `b476d55`, `41da3bd`, `a941a46`:

- **Multiple-shooting differential corrector** (`nbody/shooter.py::shoot`):
  node full-state defects in restricted n-body (planets on DE440 rails),
  `scipy.optimize.least_squares(method="lm")`, finite-difference Jacobian
  over REBOUND — documented as the Q1 baseline. Node epochs fixed (bounds
  the FD cost). Divergence sentinels: a non-converging arc contributes a
  large finite defect, never a NaN/crash (mirrors `correct.py`).
- **B-plane targeting kernel** per AAS 17-577 Eqs. 1–5 (the #142 deep-dive),
  with the published tolerances (v∞-mismatch 100–200 m/s, altitude bounds);
  Russell 2004 Eq. 5.5 independently corroborates the powered-SOI cost.
- **Near-miss survey + conic seeding** (`near_miss_survey`,
  `shooting_seed_from_near_miss`): per the #135 verdict, the shooter is
  seeded from v∞-mismatch near-misses, never a blind scan.
- **RailsEphemerisCache reuse** across all legs of a Jacobian build —
  turned an 8-minute-timeout shoot into a tractable one without touching
  the golden-gate path.
- **Bounded gate** (`a941a46`): per-leg wall budget + single-member cap +
  per-test `--timeout` — after the unbounded variant proved capable of
  running 60+ CPU-minutes.

Golden gates (two-body limit, energy drift, DE440 anchor, convergence) and
the full non-slow nbody suite are green on the committed state
(31 passed, 2026-06-07 19:50).

## The bounded gate's verdict (verbatim outcome)

Run 2026-06-07 19:34:52 → 19:50:24, `timeout 1200`, exit code 0:

```
XFAIL tests/nbody/test_shooter_jones_gate.py::test_jones_vem_nbody_rediscovers_sourced_multiset[jones-2017-vem-meevem-inbound]
  — n-body Jones shooter: rediscovery to the sourced AAS 17-577 multiset
  within 0.5 km/s. Open per the honesty boundary + #135 verdict
  (seeding/basin problem): for these 6-node VEM rows the near-miss survey
  surfaces only the high-V_inf Lambert-chain basin (max-V_inf ~29-33 km/s,
  0 bend-feasible vs the sourced Jones 2.5-5.2 km/s family), and the
  bounded n-body shoot from those seeds does not reach it. Flip ONLY when
  the member converges within VEM_VINF_TOL_KMS.
1 xfailed in 913.18s (0:15:13)
```

**The headline question — does a ballistic n-body Jones VEM member exist? —
remains OPEN.** The bounded evidence:

- Conic near-miss seeds exist only in the high-V∞ basin (~29–33 km/s);
  the n-body shoot from them does not cross to the sourced 2.5–5.2 km/s
  family within budget.
- This is consistent, not anomalous: the consolidated free-return frontier
  (Okutsu 2002 / Tito 2013 / Hughes 2014 / Donahue-Duggan 2022 — four
  sources, two independent tools, twenty years) floors ballistic Mars
  free returns at Mars V∞ ≈ 5–7 km/s; lower is explicitly gated behind
  added energy management (broken-plane/DSM/powered flybys).
- Jones's own pipeline (AAS 17-577) accepted ≤200 m/s v∞ mismatch at
  stage 1 and SNOPT-corrected in n-body — their members are not claimed
  as exact ballistic conic closures either.
- #157 (multi-rev DSM legs): making a multi-arc topology *representable*
  did not make its basin *reachable* by hopping — representability ≠
  reachability.
- #158 (continuation driver): walking model fidelity up from a closing
  circular-coplanar seed DOES reach sourced ephemeris basins (Aldrin
  6.399G1 reproduced ballistic to 0.00158 km/s). Basins are reached by
  walking, not jumping.

## Incident record (process honesty)

Two prior gate runs produced no surviving output: (1) the unbounded
variant, killed manually at 62 CPU-min; (2) the first bounded run, whose
buffered stdout was lost when a session limit killed its supervising
agent. The third run wrote through `tee` with `-s` and an exit-code
trailer; that is the run recorded above and the pattern to keep.

## Next rungs (in value order)

1. **Shooter performance pass**: parallel FD Jacobian columns
   (embarrassingly parallel, 16 cores, pool pattern exists) and/or REBOUND
   variational equations (one augmented propagation replaces all ~19 FD
   sweeps). Cost model: `n_LM_steps × (3·n_nodes+1) × n_segments ×
   t_segment` ≈ 0.5–5 CPU-h/member today; ~10–100× available. This changes
   which hunts are affordable.
2. **Powered-tolerance acceptance** per Jones: accept members within
   ≤200 m/s applied Δv (their own stage-1 criterion) and report the Δv —
   converts "never converges ballistically" into a measurable, publishable
   bound either way.
3. **Continuation into the VEM family**: #158 validated the homotopy
   driver; a VEM variant (walk e/i/phasing from a circular-coplanar VEM
   construction) is the literature-consistent route to the low-V∞ basin —
   contingent on a circular-coplanar VEM closure existing (#137 Part 2
   showed the radial-crossing primitive does not carry over; this needs
   its own genome work first).
4. **GMAT V4 rung** per the Beeson 2015 handoff recipe (auto-generated
   script carrying formulation + seed; feasibility-first + homotopy; the
   correction *cost* is the artifact signal, not raw closure error).
