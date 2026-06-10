# DSM multi-arc closure lane — RESULTS (task #180, Spec 2)

**Date:** 2026-06-10
**Code:** `src/cyclerfinder/search/dsm_descriptor_seed.py` (seed adapter +
`close_row_dsm`), `tests/search/test_dsm_descriptor_seed.py`. Commits
`57eb8c6` / `f3e52ce` (Task 1, charged layout) and `6f9e7b1` (Task 2).
**Writeback:** NONE. The lane promotes 0 rows; the catalogue is unchanged.

## Question

Can the Takao η-DSM multi-arc genome (`dsm_chain_correct`) close the **8
descriptor-bearing off-family `russell-ch4` rows** on-family — i.e. produce a
real-ephemeris two-arc trajectory whose EMERGED V∞ matches the row's sourced
Russell-table anchor — and so lift them off V0? (The 204 `russell-ocampo` rows
have no per-arc descriptor and were never in scope — descriptor-gated publication
gap.)

## Answer: NO — triple-confirmed off-family at the sourced anchors

Three independent methods agree these rows close GEOMETRICALLY but NOT at their
published V∞ anchors:

1. **Self-seeding (#177)** — runlog `data/runs/self-seeding-reachable-20260609T0210Z.jsonl`:
   all 6 REACHABLE rows came back `OFF-FAMILY-AT-ANCHOR-VINF` (`vinf_m_ok=False`).
2. **Local DSM solve** (`close_row_dsm`, charged `dsm_chain_correct` on DE440):
   no row converged. Residuals far above the 0.1 km/s criterion:

   | row | DSM residual (km/s) | sourced M-V∞ | emerged (off) |
   |---|---|---|---|
   | russell-ch4-9.353Gg2 | 131.6 | 10.52 | 21.0 |
   | russell-ch4-9.94Gg3 | 151.2 | 10.76 | 17.6 |
   | russell-ch4-3.78Gg3 | 64.0 | 4.63 | 7.3 |
   | russell-ch4-5.30ggF3 | inf | 5.44 | 11.9 |
   | russell-ch4-5.30gGf3 | — | — | "orbit does not reach body" (descriptor misses Mars) |

3. **MBH global basin-hopping** (`mbh` + `make_dsm_chain_step`, 40 hops,
   restart-bounds): `best_feasible=False`, 0 hops accepted, for both probed rows
   (`9.94Gg3` best_obj 151.2; `6.44Gg3` best_obj 74.8). Global search found NO
   feasible DSM closure anywhere. The departure E-V∞ matches its own free-var seed
   (9.94, 6.44) but downstream encounters diverge to 26–35 km/s — the chain cannot
   close at the sourced anchor.

## Interpretation

The descriptor gives the arc SHAPE but not the real-eph state. Unlike the #170
App-C rows (which carried published per-leg DE440 states to seed from), these rows
have only a coplanar descriptor, and the geometry that closes does not reproduce
the tabulated V∞. The multi-arc DSM degree of freedom — the last plausible
no-papers lever for the descriptor-bearing tail — does not bridge the gap, local
OR global. This is consistent with the validation ceiling: the off-family rows are
off-family because the data needed to land them on-family (per-member real-eph
state) was never published, not because our genome lacked a DOF.

## Status

- The DSM closure lane (seed adapter + `close_row_dsm`) is built and tested; it
  correctly records these rows as negatives (`converged=False`).
- Plan Tasks 3–7 (gate / batch / n-body / writeback) were NOT built: with a
  triple-confirmed 0-row result they would only formalise a 0-promotion outcome.
- Catalogue unchanged. The 204 ocampo rows remain descriptor-gated. Past the
  ceiling stays new-input-gated (#116 paper acquisitions) — see memory
  `validation-ceiling`.

## Reproduce the MBH probe

```python
from cyclerfinder.search import dsm_descriptor_seed as dds, dsm_leg, mbh
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.data.catalog import load_catalog
cat = load_catalog(); eph = Ephemeris("astropy")
seed = dds.seed_dsm_chain_from_descriptor(cat.by_id["russell-ch4-9.94Gg3"].raw)
step = dsm_leg.make_dsm_chain_step(sequence=seed.sequence, ephem=eph,
    bounds=seed.bounds, tol_kms=0.1, charge_flyby_continuity=True)
res = mbh.mbh(step, seed.x0, n_hops=40, rng_seed=42,
    restart_bounds=(list(seed.bounds.lower), list(seed.bounds.upper)))
assert res.best_feasible is False  # no feasible DSM closure exists
```
