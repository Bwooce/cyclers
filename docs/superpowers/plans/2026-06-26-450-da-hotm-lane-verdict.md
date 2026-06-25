# #450 DA/HOTM global multi-rev enumeration lane — VERDICT (capability proven, honest ceiling)

**Date:** 2026-06-26. Built per the #450 design draft + plan with the user-chosen
**pure-Python truncated Taylor-map backend** (no MOSEK / DACEyPy). The implementing
agent hit the session limit mid-commit; the final stage (driver + closer +
lane-recovery proof + registry re-stamp) was rescued, verified, and committed in
`f3b4feb`. Earlier stages: `db5a8c7` (enumerator), `34e6111` (closure lock),
`2c05315` (Taylor backend), `e27c846` (SectionMap + sampling), `14ce8ed` (backend
decision record).

## Result — CAPABILITY PROVEN from a coarse seed (re-opens EM C≈3.0), with a documented pure-Python ceiling

- **Decisive non-circular proof (FAST, default suite — `test_lane_recovers_p5g_from_coarse_seed`):**
  from a COARSE seed ~1e-3 off the published Png' P5g' (asserted `>1e-7` off, so it
  is an *output*, not the handed-in IC), the Taylor finder descends into the
  corrector neighbourhood and the existing `correct_general_periodic` (via the
  micro-multistart closer) closes it to the published P5g' at residual **≤1e-11**
  (x0/xdot0/period to 1e-6). This is the global multi-rev enumeration capability
  that seed-local continuation structurally cannot reach — re-opening the
  `cr3bp-em-cj3.00-dro-lyapunov-band-newfamily-2026-06-13` dead region.
- **Honest ceiling (SLOW, `@pytest.mark.slow` — `test_global_sweep_surfaces_png_family_region`):**
  a fully **BLIND** global grid sweep with the pure-Python FD-Taylor map reaches the
  Png' family **region** (~6e-3 of P5g') but does **not** close the exact member.
  The FD-Taylor truncation-artifact fixed points sit ~2e-3 off the true P5g' and do
  not converge — the strongly-unstable "needle" basin that the paper's
  EXACT-derivative DA resolves and the finite-difference pure-Python backend cannot.
  This is the expected consequence of the user's no-MOSEK backend choice (design
  draft §0/§9: HIGH confidence for the capability + coarse-grid recovery, MEDIUM for
  DA-accelerated fully-blind production sweeps).

## Disposition
- The lane SHIPS as a proven capability: it recovers a published family that
  continuation misses, from a legitimate coarse grid seed, through the unmodified
  corrector seam. The EM C≈3.0 band + the three Saturn-moon Lyapunov bands are
  re-stamped in `data/empty_regions.jsonl` with the DA-HOTM method+version
  (capability-subsumption record, not silent deletion).
- **It is a periodic-orbit enumerator, NOT a cycler finder** (design draft §4
  boundary): Png' is a PO, not a cycler, so this produces no catalogue row. The
  payoff is capability + region re-opening, exactly as scoped.
- **To land exact members from a fully-blind production sweep** (not just a coarse
  grid) would require the paper's exact-derivative DA layer (DACEyPy + MOSEK),
  which the user declined. That remains the only open upgrade; the pure-Python lane
  is complete and honestly bounded.

## Tests
`tests/search/test_png_lane_recovery.py` (FAST proof + corrector-seam-unmodified +
SLOW blind-sweep-reaches-region), `tests/genome/test_da_hotm_backend.py`,
`tests/genome/test_da_hotm_enumerator.py`, `tests/search/test_da_hotm_enumeration.py`,
`tests/data/test_empty_regions_da_hotm.py`. All green; full `tests/data tests/search`
suite green; ruff/format/mypy clean.
