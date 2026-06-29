# #480 — real-ephemeris UNGUIDED discovery: the IEG ballistic-cycler CLASS is reproduced; the exact Table-4 member is not (yet) matched

**Date:** 2026-06-30
**Method:** unguided real-eph (jup365) search — minimize total ballistic flyby ΔV +
feasibility ONLY (no V∞-to-Table-4 pull); V∞ is an OUTPUT compared to Table 4 after.
Driver (scratch): inline foreground (`scripts/_eggie_realeph_480.py` primitives).
**Context:** the corrected ballistic construction (`search/eggie_ballistic.py`, equal-in/out
|V∞| enforced) + the repeated-encounter root-cause fix
(`docs/notes/2026-06-29-480-eggie-ballistic-construction-verdict.md`).

## Why this run (methodology correction)

The earlier real-eph search had a `+4000·verr` term pulling V∞ toward Table 4 — i.e. it
STEERED toward the paper's answer (the subtle cousin of hardcoding; flagged by the user).
It hadn't manufactured a false success (still infeasible), but the method was impure. This
run removes the V∞ target entirely: discover whatever feasible ballistic cycler exists,
then CHECK if Table 4 appears on its own.

## Result (no V∞ steering)

9 feasible, ballistic (ΔV ≲ 0.2 m/s), equal-Ganymede IEG triple cyclers found in real
ephemeris near the paper's epoch, all flyby altitudes above-surface. Distinct members:

| total ΔV | dep off | V∞ E / G(=G) / Io | flyby alts (km) |
|---|---|---|---|
| 0.00 m/s | +10.84 d | 8.88 / 6.45 / 7.28 | 52904 / 21893 / 9089 / 2277 |
| 0.00 m/s | +3.78 d | 9.00 / 6.59 / 7.58 | 20216 / 12489 / 6678 / 112 |
| 0.16 m/s | +10.67 d | 10.03 / 7.57 / 9.54 | 807 / 594 / 9899 / 246 |

**Table-4 member (9.12 / 7.07 / 8.38) discovered unguided? NO** — 0 of 9 within 0.3 km/s.

## Honest reading

- **POSITIVE — the class is reproduced.** Our pipeline INDEPENDENTLY discovers feasible,
  ballistic, equal-Ganymede Io-Europa-Ganymede triple cyclers in real Galilean ephemeris,
  with no steering. That confirms Hernandez-Jones-Jesick 2017's central claim (such
  ballistic IEG triple cyclers exist) and the equal-same-body-V∞ property holds in every
  member. This is the payoff of fixing the repeated-encounter construction bug — the
  earlier "~0.1 km/s wall" / "no feasible closure" verdicts were artifacts of the broken
  (structurally non-ballistic) seed and are now SUPERSEDED
  ([[feedback_bugfix_invalidates_past_searches]]).
- **NOT a reproduction of the exact Table-4 member.** The discovered members BRACKET Table 4
  (E 8.88 ↔ 10.03 around 9.12; G 6.45 ↔ 7.57 around 7.07), so the specific member is almost
  certainly in the family between them — but the coarse scan (45 epochs/plan, Nelder-Mead,
  2 branch plans) did not land it. No claim of Table-4 reproduction is made.

## Discipline note (the steering trap)

Matching a published number by putting it in the objective is the subtle form of hardcoding
("it matched!" = danger, like "it closed!"). The legitimate test is unguided discovery +
post-hoc comparison, recorded here. Captured in
[[feedback_constructed_tour_per_encounter_self_consistency]].

## Next (no steering)

1. **Finer unguided scan** for the exact member: denser epoch grid + finer ToF/branch
   enumeration BETWEEN the bracketing members, V∞ still an output. If a feasible member
   lands at 9.12/7.07/8.38 → exact-member reproduction (then level-3 n-body confirm).
2. **No catalogue self-admission** — a reproduced published tour is a human-admitted
   V4-ceiling result, not self-admitted.
3. Golden (`tests/verify/test_ieg_reproduction_golden.py`) stays skipped until an
   exact-member feasible close (then un-skip WITHOUT loosening tolerances).
