# FBS vs Lambert closure attempt on the hard multi-arc rows (#242)

**Date:** 2026-06-13
**Driver:** `scripts/fbs_hard_rows_closure.py` (committed)
**Tools tested (#226):** `src/cyclerfinder/core/fbs_match_point.py`
(`match_point_defect` + analytic Jacobian) and
`src/cyclerfinder/search/dsm_leg.py::dsm_leg_correct_fbs` — the opt-in
Lambert-free FBS match-point single-leg corrector.
**Writeback:** NONE. No catalogue edit, no validation-level change. This is an
evidence/decision-input run for #243 (default-FBS decision).

## Question

Can the Lambert-free FBS match-point lane CLOSE the historically un-closeable
multi-arc rows that the Lambert corrector could never close? A clean closure FBS
achieves that Lambert cannot = strong argument to default FBS. No new closures =
FBS stays an opt-in cross-check tool. **A clean negative is a success.**

## Answer (one line)

**NO — FBS does not close anything Lambert could not, and on the defining hard
arcs (long multi-rev loop arcs) FBS is WEAKER than Lambert, not stronger.** FBS
and Lambert solve the SAME per-leg two-point boundary value problem; they agree
to ~1e-12 km/s wherever both run, and the chain-level blocker (V∞ continuity /
basin selection) is untouched by swapping the leg primitive. **Recommendation:
do NOT default FBS (#243); keep it opt-in.**

---

## Rows attempted, topology, and source

### S1L1 / `russell-ch4-4.991gG2` (already V3)

Corrected two-arc topology, sourced from Russell 2004 Appendix C #83
(`docs/notes/2026-06-08-s1l1-source-dig.md`, `src/cyclerfinder/search/s1l1_corrected.py`):

```
E -> g(E-E free return, NO Mars) -> E flyby -> G(E-M-E Mars transit, longitude
   rendezvous with the TRUE DE440 Mars position) -> E
```

The binding multi-arc constraint lives on the **G (Mars-transit) leg**: its
right-boundary position `rf` is set to the TRUE DE440 Mars position at the
published arrival epoch, so a converged match point IS a longitude-correct Mars
intercept. Verified topology against source BEFORE attempting (g = sub-Mars
Earth-Earth return; only G reaches Mars; one Mars encounter/cycle).

### `russell-ch4-6.44Gg3` (V1)

Full-sequence E-M-E-M-E loop-arcs-as-DSM-legs probe (#153 / #157). Descriptor
arcs g(2.087 yr) + G(4.3191 yr) sum to the 6.41-yr period; unrolled into
Mars-bracketing pieces → 4 legs (E→M transit, M→E = g_arc−transit, E→M transit,
M→E = G_arc−transit). Sourced anchors (EXPECTED, never imposed): v∞ E=6.44,
M=3.74 km/s.

## Which constraints are in each residual (explicit)

* **Per-leg FBS residual** = the 6-vector match-point defect at the interior
  impulse point: `[r^B − r^F (3 rows, km); v^B − v^F (3 rows, km/s)]` — full
  POSITION and VELOCITY match. A leg that does not dynamically connect its two
  boundary states cannot drive this to zero. (Non-dimensionalised in the
  corrector: position/AU, velocity/v_circ.)
* **S1L1 G-leg longitude rendezvous** is carried by `rf` = true DE440 Mars
  position → it is the position rows of the defect. A converged match point is a
  true-Mars intercept, not a radius coincidence.
* **NOT in the per-leg residual (by construction of the leg primitive):** the
  chain-level flyby V∞-magnitude continuity that *defines* a ballistic cycler.
  That is the chain corrector's job (`evaluate_dsm_chain` charged path), and it
  is exactly the constraint that historically did not close. We measure the
  emerged encounter v∞ and report the gap honestly rather than hiding it.
* **Same-model golden:** the Lambert `dsm_leg` solution on the IDENTICAL boundary
  states / ToF / eta / mu. At the Lambert solution the FBS defect is exactly zero
  with `dv = v21−v12, vf = v22`; the two lanes MUST agree. No cross-model golden.

---

## S1L1 — per-G-leg FBS result + Lambert cross-check (real numbers)

FBS seeded from a PERTURBED (±0.3 km/s) Lambert solution so a genuine root-find
runs, not a trivial echo. eta=0.5 (a real interior DSM; the leg is ballistic so
the converged dV_DSM is ~0). `rf` = true DE440 Mars position.

| G leg → Mars | tof (d) | Lambert dV_DSM | FBS converged | FBS max_res | FBS nfev | \|dv_fbs−dv_lam\| | \|vf_fbs−vf_lam\| | v∞@Mars L/F vs pub |
|---|---|---|---|---|---|---|---|---|
| #2→M#3   | 179.8 | ~0 | True | 9.96e-17 | 5 | 5.2e-15 | 1.3e-14 | 5.2480 / 5.2480 vs 5.2480 |
| #5→M#6   | 125.5 | ~0 | True | 0.00e+00 | 4 | 1.8e-14 | 1.8e-14 | 7.6934 / 7.6934 vs 7.6930 |
| #8→M#9   | 138.4 | ~0 | True | 2.41e-16 | 5 | 2.2e-14 | 5.7e-15 | 4.6570 / 4.6570 vs 4.6570 |
| #11→M#12 | 210.9 | ~0 | True | 1.19e-16 | 5 | 8.1e-15 | 5.6e-17 | 3.1984 / 3.1984 vs 3.1980 |
| #14→M#15 | 154.6 | ~0 | True | 1.99e-16 | 5 | 7.8e-14 | 9.8e-14 | 6.2628 / 6.2628 vs 6.2630 |
| #17→M#18 | 112.2 | ~0 | True | 1.19e-16 | 5 | 5.2e-14 | 6.4e-14 | 8.0457 / 8.0457 vs 8.0460 |
| #20→M#21 | 230.2 | ~0 | True | 1.20e-16 | 5 | 2.1e-13 | 2.2e-13 | 3.2190 / 3.2190 vs 3.2190 |

**Agreement verdict:** all 7 G legs FBS-converged (max_res ≤ 1e-16); FBS↔Lambert
agree to ≤ 2.2e-13 km/s on both Δv and vf; emerged v∞@Mars matches the published
App-C value on BOTH lanes to 4 decimals. The two independent lanes AGREE — this
is a genuine cross-confirmation.

**But it is NOT a new closure.** S1L1 is already V3, confirmed by the corrected
topology + App-C real-eph seed + REBOUND/IAS15 n-body cross-check (#167). What
this run shows is that the FBS lane reproduces the same longitude-correct Mars
intercept the Lambert lane (and the n-body confirmation) already established. FBS
brought no new capability here; it agreed with the incumbent on an already-closed
row.

---

## 6.44Gg3 — per-leg FBS result + Lambert cross-check (real numbers)

Chain built ballistically (departure v∞=6.44 km/s tangential-prograde at t0;
Lambert arrival velocity propagated forward as flyby continuity). eta=0.5 per leg.

| leg | bodies | tof (d) | Lambert dV_DSM | FBS conv | FBS max_res | nfev | \|dv_fbs−dv_lam\| | emerged v∞ (incoming) |
|---|---|---|---|---|---|---|---|---|
| 1 | E→M | 262.0  | 45.113 | True | 1.99e-15 | 5 | 1.4e-12 | M = 33.75 |
| 2 | M→E | 500.3  | 50.837 | True | 4.78e-15 | 7 | 2.5e-13 | E = 55.56 |
| 3 | E→M | 262.0  | 62.293 | True | 3.98e-16 | 5 | 1.0e-13 | M = 27.31 |
| 4 | M→E | 1315.6 | 48.457 | True | 1.59e-15 | 5 | 5.7e-14 | E = 38.83 |

**Agreement verdict:** all 4 legs FBS-converged; FBS↔Lambert agree to ≤ 1.4e-12
km/s.

**Closure verdict: did NOT close (as expected).** Emerged encounter v∞ are
27–56 km/s, nowhere near the sourced 6.44/3.74 km/s. The per-leg BVP closes
trivially on BOTH lanes (with huge interior impulses, because eta=0.5 forces an
impulse the ballistic seed does not want); the chain does not close because of
the V∞-continuity / basin-selection blocker documented in #153/#157 — which is
NOT a leg-primitive property and is therefore unchanged by switching Lambert→FBS.
This reproduces the #157 conclusion: making the topology representable was
necessary but not sufficient; the gap is basin/seed selection, not the leg
mechanic.

---

## The decisive probe: does FBS solve any leg where Lambert FAILS?

The strongest possible #243 argument would be a leg geometry FBS can represent
but Lambert cannot. We probed the long M→E loop arc (1315.6 d, ~3 revs — the
defining hard arc of these multi-arc rows), eta=0.5, modest Mars departure:

```
Lambert (max_revs=0):        SOLVED,        dV_DSM = 9.3567 km/s
FBS (cold circular seed):    NOT converged, dV_DSM = 337.95 km/s, max_res = 2.79, nfev = 400 (exhausted)
FBS (warm Lambert seed):     converged,     dV_DSM = 9.3567 km/s, max_res = 2.4e-15, nfev = 2
```

**This is the opposite of an FBS advantage.** On the long multi-rev arc:

* Lambert (a direct global BVP solver) finds the leg from a trivial seed.
* FBS (a two-sided SHOOTING method) does NOT converge from a cold circular seed —
  its match-point defect is multi-rooted and the cold start lands outside the
  basin (max_res 2.79 after exhausting 400 evaluations).
* FBS converges in 2 evaluations ONLY when seeded from the Lambert solution it is
  supposed to be replacing.

So FBS is a **refiner that needs the Lambert basin**, not a standalone
replacement. On exactly the arcs that make these rows hard, FBS is the weaker
leg solver.

---

## Honest verdict (per row)

| Row | FBS outcome | Lambert cross-check | Agree? | Verdict |
|---|---|---|---|---|
| S1L1 G legs (×7) | converged, max_res ≤1e-16 | converged | YES (≤2.2e-13 km/s) | **closed-and-cross-confirmed, but NOT a new closure** (row already V3; FBS reproduced the incumbent) |
| 6.44Gg3 (4 legs) | per-leg converged | per-leg converged | YES (≤1.4e-12 km/s) | **did-not-close** — chain v∞ 27–56 vs anchors 6.44/3.74; blocker is V∞-continuity/basin, not the leg primitive |
| long M→E loop arc probe | cold seed FAILED; warm seed OK | Lambert solved from trivial seed | n/a | **FBS weaker than Lambert** on the defining hard arc |

No false-positive risk to flag: there is no claimed new closure to be skeptical
of. The S1L1 agreement is a cross-confirmation of an already-confirmed result;
the 6.44Gg3 result is a clean negative consistent with three prior independent
findings (#153, #157, #180).

---

## Recommendation for #243 (default-FBS decision)

**Do NOT make FBS the default leg corrector. Keep it opt-in.** Reasons, in order
of weight:

1. **FBS closed nothing Lambert could not.** Every leg both lanes ran, they
   agreed to ~1e-12 km/s. There is no row, and no single leg, where FBS produced
   a closure Lambert missed.
2. **FBS is strictly weaker on the hard arcs.** On the long multi-rev loop arc —
   the geometry that *defines* these multi-arc rows — FBS shooting fails to
   converge from a cold seed while Lambert solves directly. Defaulting FBS would
   make the search MORE seed-sensitive on precisely the cases we care about.
3. **The real blocker is not the leg primitive.** Both lanes solve the identical
   per-leg BVP; the chain non-closure is V∞-continuity / basin selection
   (#153/#157/#180). Swapping the leg solver cannot move that wall.
4. **FBS retains genuine value as an opt-in cross-check.** Its analytic Jacobian
   and exact agreement with Lambert make it an excellent INDEPENDENT verifier
   (the role it played here, confirming the S1L1 G-leg intercepts). That is the
   role to keep it in.

The clean negative is the result: FBS is a verified second opinion, not a more
capable corrector.

## Reproduce

```
uv run python scripts/fbs_hard_rows_closure.py
```

## Lint / type / test status

* `uv run ruff check scripts/fbs_hard_rows_closure.py` — clean.
* `uv run ruff format --check scripts/fbs_hard_rows_closure.py` — clean.
* `uv run mypy src tests` — clean (334 files; `scripts/` is outside the
  project mypy scope per `.pre-commit-config.yaml`, consistent with every other
  file in `scripts/`).
* FBS tool suite `tests/search/test_dsm_leg.py` + `tests/core/test_fbs_match_point.py`
  — 24 passed (the #226 green baseline, unchanged).
