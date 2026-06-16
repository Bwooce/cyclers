# #322 -- tulip petal_count topology gate: z0-collapse bug fix

Surfaced by #313's structural negative.  `find_tulip_at_system` was certifying
"3D tulip" matches at very small mu (Mars-Phobos, Mars-Deimos) on the basis of
`petal_count == np_target` alone.  At those systems the symmetric corrector
drives the seed's `z0` to ~0 and the orbit collapses to a PLANAR Np-petal
orbit -- which shares the in-plane petal count with a genuine 3D tulip but is
NOT a tulip.

## The bug

Pre-fix Tier A gate (`src/cyclerfinder/genome/tulip.py`, line 945 pre-fix):

```python
if n_direct == np_target:
    # Tier A hit: the seed IS already the target Np tulip at this mu.
    return FindTulipResult(..., success=True, reason="direct_seed_match")
```

`petal_count` counts local minima of the spacecraft-to-secondary distance in 3D
-- but a planar orbit with z(t) identically zero still produces exactly Np
in-plane perilune minima.  There was NO complementary check on z0 magnitude or
out-of-plane excursion.  Result: spurious `direct_seed_match` claims at every
system where the corrector drove z0 to ~0.

Evidence rows (from `data/tulip_sweep_281.jsonl`, pre-fix):

| System | np_target | z0 (corrected) | petal_count | Old verdict |
| --- | --- | --- | --- | --- |
| Earth-Moon | 2 | 0.173 | 2 | direct_seed_match (CORRECT) |
| Jupiter-Ganymede | 2 | 4.5e-19 | 2 | direct_seed_match (FALSE POS) |
| Saturn-Titan | 2 | 8.5e-23 | 2 | direct_seed_match (FALSE POS) |
| Neptune-Triton | 2 | 3.6e-14 | 2 | direct_seed_match (FALSE POS) |
| Pluto-Charon | 2 | 1.5e-18 | 2 | direct_seed_match (FALSE POS) |

The Earth-Moon row's z0 (0.173) is the only one above the published Koblick
Np=15 family floor (0.046).

## The fix

Added to `src/cyclerfinder/genome/tulip.py`:

1. **Constant** `TULIP_Z_AMPLITUDE_FLOOR_NONDIM = 5e-3` -- the floor for the
   out-of-plane amplitude gate.  Justification: the smallest z0 in the
   Koblick 2023 AMOSTECH Table 4 paper-row table is Np=15 at z0=0.045796
   nondim; 5e-3 sits ~9x below that floor, so it is a CONSERVATIVE cliff
   (anything below 5e-3 is unambiguously planar-collapse rather than a
   sourced family member).  This is a Phase 1 floor; if Koblick's Fig 7
   out-of-plane envelopes are ever digitised we may tighten this against
   the family-specific minimum-z amplitude.  At Earth-Moon LU=389703 km
   this is ~1948 km out-of-plane; at Mars-Phobos LU=9375 km it is 46.9 km.

2. **Function** `out_of_plane_amplitude(state0, period, system)` -- propagates
   the orbit one period in physical time (DOP853, rtol/atol=1e-11) on a
   dense 401-point grid and returns `max|z(t)|`.

3. **Function** `is_three_dimensional(state0, period, system, *, z_floor)` ->
   `(is_3d, max_abs_z)` -- the topology gate.  An orbit is 3D iff BOTH
   `|z0| >= z_floor` AND `max|z(t)| >= z_floor`.  The `|z0|` short-circuits
   cheaply on the planar-collapse case (no integration needed).

4. **`FindTulipResult`** gained two backward-compatible optional fields:
   `topology_verdict: Literal["3D tulip", "planar Np-petal collapse",
   "petal count mismatch", "seed no converge", "unknown"]` and
   `max_abs_z: float | None`.  All five existing callers (`tulip.py` and the
   `tulip_multimu_sweep.py`/scan_313 scripts) read only `.switched`,
   `.success`, `.reason` -- unchanged.

5. **`find_tulip_at_system`** wired the gate: after `n_direct == np_target`,
   call `is_three_dimensional`; only declare `direct_seed_match` if `is_3d`.
   On planar collapse the Tier A path falls through to Tier B (continuation
   + family-switch).  Tier B's switched member is ALSO classified via the
   gate; if Tier B also produces a planar collapse the result is
   `success=False` with `topology_verdict="planar Np-petal collapse"` and
   `reason="tier_b_planar_collapse"`.

## Regression tests

Added to `tests/genome/test_tulip.py`:

* `test_322_marsphobos_paperseed_np4_does_not_claim_tulip` -- the canonical
  #313 false-positive case.  Asserts `success=False` AND
  `topology_verdict == "planar Np-petal collapse"`.
* `test_322_earth_moon_genuine_3d_tulip_still_admitted` -- the canonical
  Earth-Moon Np=2 via the Tier B continuation+family-switch path.  Asserts
  `success=True`, `topology_verdict == "3D tulip"`, `max_abs_z >= floor`.
* `test_322_threshold_sanity_at_half_and_double` -- numeric sanity:
  `is_three_dimensional` rejects z0=0.5*floor and accepts z0=2*floor.

All 10 tulip tests pass (7 prior + 3 new).  All 55 genome tests pass.

## Re-verification: past JSONLs under the new gate

`scripts/reverify_tulip_topology.py` re-evaluates every JSONL row that
claimed `converged=True AND reason in {ok, direct_seed_match}` against the
new gate.  Output: `data/tulip_topology_reverify_322.jsonl`.

| Source JSONL | Total rows | Claimed tulip | Still 3D | Downgraded |
| --- | --- | --- | --- | --- |
| `tulip_sweep_281.jsonl` (#281 multi-mu Np=2) | 14 | 5 | 1 | 4 |
| `tulip_higher_np_283.jsonl` (#283 higher-Np) | 8 | 5 | 4 | 1 |
| `scan_313_mars_phobos.jsonl` (#313) | 6 | 1 | 0 | 1 |
| `scan_313_mars_deimos.jsonl` (#313) | 6 | 1 | 0 | 1 |
| `scan_313_sun_jupiter_io.jsonl` (#313) | 54 | 0 | 0 | 0 |
| `scan_313_sun_jupiter_europa.jsonl` (#313) | 54 | 0 | 0 | 0 |
| TOTAL | 142 | 12 | 5 | 7 |

Downgrade detail (all are planar-collapse artifacts; corrected z0 below 1e-13):

| Source row | System | Np | z0 corrected | max\|z\| |
| --- | --- | --- | --- | --- |
| 281#4 | jupiter-ganymede | 2 | 4.5e-19 | 4.5e-19 |
| 281#8 | saturn-titan | 2 | 8.5e-23 | 8.5e-23 |
| 281#12 | neptune-triton | 2 | 3.6e-14 | 3.6e-14 |
| 281#13 | pluto-charon | 2 | 1.5e-18 | 1.5e-18 |
| 283#5 | saturn-titan | 3 | -2.6e-26 | 2.6e-26 |
| 313#5 (mars-phobos) | -- | (seed Np=4, target Np=2) | -1.4e-14 | 1.4e-14 |
| 313#5 (mars-deimos) | -- | (seed Np=4, target Np=2) | -5.6e-17 | 5.6e-17 |

Confirmed-3D survivors (real tulips):

* `tulip_sweep_281` row 0: Earth-Moon Np=2 (z0=0.173, max\|z\|=0.173).
* `tulip_higher_np_283` rows 0-3: Earth-Moon Np=3,4,5,6 (z0=0.109 to 0.159).

## Honest assessment: catalogue impact

**The bug propagated NO false novelty claims.**  Reasoning:

1. The 4 downgraded #281 rows (Jupiter-Ganymede, Saturn-Titan, Neptune-Triton,
   Pluto-Charon) were never written back to `data/catalogue.yaml` -- they
   were a JSONL discovery sweep awaiting topology confirmation before any
   catalogue action.  Per `feedback_orbit_closure_discipline`: hold writeback
   till confirmed.  These rows confirmed nothing reproducible.
2. The 1 downgraded #283 row (Saturn-Titan Np=3) likewise sat in JSONL only.
3. The 2 downgraded #313 rows (Mars moons) were the structural-negative
   evidence that SURFACED this bug -- never claimed as positive in `#313`'s
   verdict document `docs/notes/2026-06-16-313-multi-system-scouts-verdict.md`.

The 5 surviving Earth-Moon rows (Np=2,3,4,5,6) all reproduce sourced Koblick
2023 / pumpkyn anchors -- they were already cross-checked against the
upstream IC table, so the bug fix does not retroactively change their
provenance.  No catalogue retraction needed.

## Phase 2 (deferred)

The Phase 1 fix handles the obvious z0=0 collapse mode.  Phase 2 candidates:

* If Koblick's Fig 7 out-of-plane envelopes are ever digitised, tighten the
  floor to per-Np values rather than the conservative 5e-3 cliff.
* If a torus-near-tulip topology surfaces in the discovery sweep, the gate
  may need a winding-number cross-check on top of `petal_count` + max\|z\|.
  No evidence such a topology exists in the searched (mu, x0) basin yet.
* Reverify any FUTURE tulip JSONL emitted by Tier-A-only code paths against
  this gate.  The fix lives at `find_tulip_at_system` -- callers that build
  custom Tier A logic must invoke `is_three_dimensional` themselves.
