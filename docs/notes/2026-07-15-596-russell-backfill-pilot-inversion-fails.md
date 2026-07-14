# #596 follow-up: Russell Table 3.4 -> catalogue backfill pilot — naive inversion FAILS

**Date:** 2026-07-15
**Context:** #596 identified the Russell dissertation per-leg `a_au`/`tof_days` backfill (~216+38
rows) as "ready to execute, source in hand" after correcting `MISSING_DATA.md`'s staleness. Before
attempting the full backfill, ran a pilot inversion against the Aldrin cycler (row `1.0.1.-1` in
Table 3.4), which is already correctly catalogued (`aldrin-classic-em-k1-outbound`, sourced
`a_au=1.60`, `e=0.393` from multiple independent references) — the obvious positive control.

## Method tried

Table 3.4 gives, per row: AR (Aphelion Ratio), TR (Turn Ratio), "Earth→Mars (or aphelion) Time
(days)", V∞ at Earth, V∞ at Mars, geocentric turn angles. None of these is `(a_au, e)` directly.
Attempted to invert using this project's own `cyclerfinder.search.free_return.free_return_geometry(a_au,
e)` forward map (the same radial-crossing machinery the #593/#544 QBCP work and the free-return genome
already use) via `scipy.optimize.least_squares`, fitting `(a_au, e)` to match Russell's `tof_em_days`
and `V∞_Earth` (with AR/`V∞_Mars` held out as independent cross-checks — per this project's
independent-cross-check-mandatory discipline).

## Result: clean negative, NOT a quick win

- Local fit near the known answer (seed `a=1.6, e=0.4`): converges to `a=1.618, e=0.382` — close on
  AR (1.467 vs. target 1.47) but **badly misses ToF (124.4d vs. 146d target, 15% off) and V∞_Earth
  (5.23 vs. 6.5 km/s target, 20% off)**.
- Global multi-start search (29 starting points spanning `a∈[1.1,2.9]`, `e∈[0.05,0.85]`): found a
  DIFFERENT, unrelated `(a,e)` = `(1.445, 0.308)` that matches ToF almost exactly (145.9 vs. 146d)
  but misses AR badly (1.24 vs. 1.47) and V∞ badly (4.28/7.35 vs. 6.5/9.7 km/s).

**No `(a,e)` pair simultaneously satisfies all three of Russell's Table 3.4 columns (AR, ToF, V∞)
under this project's `free_return_geometry` convention.** This is not a local-minimum/bad-seed
problem (the global search confirms no consistent solution exists in the searched space) — it's a
genuine modeling/convention mismatch.

## Leading hypothesis (NOT verified)

Table 3.4's own column header is a hint: **"Earth→Mars (or aphelion) Time (days)"** — the parenthetical
suggests Russell's ToF column measures time-to-APHELION for some rows and time-to-Mars-crossing for
others (presumably rows where the ellipse's aphelion is near Mars's orbit, AR≈1, vs. rows where it
overshoots). `free_return_geometry`'s `tof_em_days` always measures time-to-Mars-radial-crossing,
never time-to-aphelion. This alone could explain the ToF mismatch. A second candidate cause: Russell's
V∞ definition (Ch.2/Ch.3, per the 2026-06-07 method-mining note) may reference a different frame or
flyby-turn convention than this project's raw heliocentric-relative-velocity-at-crossing.

**Neither hypothesis has been checked against Russell's actual Ch.3 text** (the method-mining note
covers Ch.2.7/Ch.3.5-3.8/Ch.5.2-5.4 but wasn't read with this specific ToF/V∞-definition question in
mind) — this is the next actual step, not a rewrite of the inversion code.

## Disposition

**Do NOT attempt the full 216+38-row backfill until this modeling question is resolved.** Per
`feedback_orbit_closure_discipline` (independent cross-check mandatory, verify topology vs. source
FIRST) and `feedback_golden_tests_sourced_only`, writing derived `(a_au, e)` values into ~250
catalogue rows without a working, positive-control-verified derivation would risk seeding the
catalogue with confidently-wrong numbers across a quarter of all `russell-ocampo-*` entries — a much
worse outcome than the current honest `data_gaps` state.

**Next step if resumed:** re-read Russell 2004 Ch.3.5-3.8 (pp. 71-84, already partially covered by the
method-mining note) specifically for the ToF-column and V∞-column exact definitions, before writing
any new inversion code. This is real domain-reading work, not a mechanical task — scope accordingly
(likely a multi-hour dedicated investigation, not a quick follow-on).
