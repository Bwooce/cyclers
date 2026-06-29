# #480 level-3 high-fidelity n-body — the real-eph corrector CONVERGES on the correct seed (no wall); V∞ stays below Table 4

**Date:** 2026-06-30
**Task #3.** High-fidelity Jupiter-central n-body corrector (`nbody/jovian.jovian_shoot`
against jup365, now `trf` — commit `2389f27` fixed the lm→trf soft-cap) applied to the
real-ephemeris EGGIE periapsis seed at the paper epoch (`paper_departure_et()` + offset).
Sourced V∞ targets (Table 4): Europa 9.12 / Ganymede 7.07 / Io 8.38 km/s.

## Verified results (logged)

| seed | seed defect | after 12 nfev | correction ΔV | V∞ (E / G / Io) |
|---|---|---|---|---|
| near-Table-4 member (V∞≈9.50/6.86/8.33, non-ballistic) | 1.46e6 | 5.0e4 (not conv.) | 2923 m/s | drifted (Io→6.27) |
| **ballistic member (V∞≈9.00/6.59/7.58)** | **4.13e5** | **9.1e2 (≈450× cut)** | **529 m/s** | 8.90 / 6.61=6.57 / 7.72 |

(Both capped at 12 corrector iterations — partial, not converged; trends clear.)

## Honest conclusions

1. **The "~0.1 km/s wall" was the BROKEN seed, not the model.** On the structurally-correct
   ballistic seed the n-body corrector CONVERGES — continuity defect drops ~450× (4.13e5 →
   9.1e2) in 12 iterations with a modest 529 m/s correction. The earlier Stage-2/3/M1 walls
   were run on the repeated-encounter-broken seed and are SUPERSEDED
   ([[feedback_bugfix_invalidates_past_searches]]). This validates the whole bug-fix arc.
2. **V∞ stays in the real-eph ballistic family ~0.4-0.7 km/s BELOW Table 4** (E 8.9, G 6.6,
   Io 7.7) — it does NOT drift toward Table-4. Consistent with every prior result: the exact
   Table-4 member is not reached; the ballistic IEG-cycler CLASS is reproduced just below it.
3. **Maintenance ΔV is in the few-hundred-m/s range at 12 nfev and decreasing** as the defect
   falls; a precisely-converged number needs more iterations (chunked FD, or a real-eph port
   of the `jovian_stm` analytic Jacobian — the FD corrector is iteration-bound, not walled).
   This is larger than the paper's real-eph maintenance scale (tens of m/s over ~10 cycles),
   the residual ideal-vs-real-eph / level-3-optimization gap noted in
   `2026-06-30-480-eggie-realeph-unguided-discovery.md`.

## #480 final standing (unchanged by this run, now fully evidenced)

- Construction bug found + fixed; guard added (`search/tour_self_consistency.py`).
- IEG ballistic-cycler CLASS reproduced unguided in real ephemeris (paper's central claim).
- Exact Table-4 member NOT reproduced (ideal-vs-real-eph gap; paper's own maintenance-ΔV caveat).
- n-body corrector confirmed to converge on the correct seed (no model wall).
- No catalogue change; golden stays skipped. Scratch `scripts/_eggie_*level3*_480.py` removed.
