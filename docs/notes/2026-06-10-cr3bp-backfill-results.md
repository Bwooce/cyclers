# CR3BP Earth-Moon Backfill Results

**Run timestamp:** 2026-06-10T03:23:36Z
**Script:** `scripts/cr3bp_backfill.py`
**Elapsed:** 0.5 s

**Status:** NO catalogue writeback. Results are PROPOSED for human review only.

---

## Row: `arenstorf-em-figure8-1963`

**Backfill outcome: CONVERGED**

- IC source: Hairer, Nørsett & Wanner 'Solving ODEs I' p.129 (B5); Arenstorf 1963
- μ used: `0.012277471`
- Initial state (nd): `[0.994, 0.0, 0.0, 0.0, -2.0015851063790824, 0.0]`
- Period guess (nd): `17.065216560157964`

**Corrector output:**

| Field | Value |
|---|---|
| converged | `True` |
| closure_residual | `7.749e-11` |
| period_nd (corrected) | `17.0652165601594` |
| jacobi_constant | `2.85641252` |
| state0_nd (corrected) | `[0.9940000000000017, -2.034961361182452e-15, 7.573064689951791e-29, 1.2909925227679276e-17, -2.0015851063788865, 8.213520575414745e-26]` |
| lunit_km | `384400.0` |
| tunit_s | `372931` |

**PROPOSED `orbit_elements.cr3bp` fields (review-gated, NO writeback):**

```yaml
    cr3bp:
      mass_ratio: 0.012277471
      jacobi_constant: 2.85641252
      period_nd: 17.0652165601594
      state_nd: [0.9940000000000017, -2.034961361182452e-15, 7.573064689951791e-29, 1.2909925227679276e-17, -2.0015851063788865, 8.213520575414745e-26]
      lunit_km: 384400.0
      tunit_s: 372931
```

**PROPOSED `_LEVEL_EVIDENCE` line (Arenstorf row only):**

The Arenstorf IC is sourced from Hairer et al. (1993), a citable published reference, and the corrector converges to closure < 1e-10. This meets the criteria for promotion from V0 (citation-only) to a higher validation level once the Jacobi constant is cross-checked against an independent source (e.g. the JPL three-body periodic-orbit catalog).

_Proposed level: **V1** (computed from sourced IC; pending independent Jacobi cross-check for V2/V3)._

---

## Row: `genova-aldrin-2015-em-3petal-cycler`

**Backfill outcome: NO_SOURCED_IC — skipped**

> No published initial conditions available from the accessible NTRS abstract (NTRS 20150018049). Full AAS-15 PDF was inaccessible at ingest; state_nd is null in the catalogue row. Per honesty rules, NO IC is fabricated. Backfill requires the full Genova & Aldrin 2015 paper PDF.

No CR3BP fields can be proposed for this row.

---

## Row: `wittal-2022-em-cycler-family`

**Backfill outcome: NO_SOURCED_IC — skipped**

> No published initial conditions available from the accessible NTRS abstract (NTRS 20220013595 / IAC-22-C1.6.6). The catalogue row is a family seed (state_nd: null). Per honesty rules, NO IC is fabricated. Backfill requires the full Wittal, Miaule & Asher 2022 paper PDF.

No CR3BP fields can be proposed for this row.

---

## Summary

- **Converged (sourced IC + periodic orbit found):** 1
  - `arenstorf-em-figure8-1963`
- **No sourced IC (skipped, not fabricated):** 2
  - `genova-aldrin-2015-em-3petal-cycler`
  - `wittal-2022-em-cycler-family`

**NO writeback to `data/catalogue.yaml` or `validate.py`.**  Promotion is review-gated (separate step).
