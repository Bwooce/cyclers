# M7 Phase 2 coverage scan — catalogue-wide maintenance band is data-gated (#423)

**2026-06-23.** Ran `scripts/m7_phase2_coverage.py` (read-only, no writeback) over
every cycler-type catalogue row: `verify_real_closure(compute_tcm=True, n_cycles=2)`,
recording whether the real-eph cycler CONSTRUCTS, CLOSES, and what horizon TCM M7
returns. This is the diagnostic for "which rows can M7 actually measure a real-eph
maintenance band for".

## Result (318 cycler-type rows)

```
measured (finite TCM)  =   8
diverged (inf)         =   5
skip (construct/window/error) = 305
```

- **305 / 318 skip** — they never construct. The dominant ERROR (196) is the
  `russell-ocampo-*` census family ("catalogue entry … " — descriptor-only rows with
  no per-arc free-return geometry), plus `hollister-menning` (15), `russell-ch*` (14),
  and rows with `0 legs` / malformed leg dicts. This is the **#388 family-selection /
  publication-gap wall**: the rows carry summary invariants, not the rev-correct
  departure velocities and constructed legs M7 needs to seed `target_leg`.
- **The 8 "finite" measurements are not usable maintenance bands.** Two are degenerate
  `0.0` (`russell-ocampo-2.5.1+0`, `jones-2017-vem-emevve-outbound`); the rest are the
  single-rev-Lambert-seed artifact — `aldrin-classic-em-k1-outbound` 53,962 m/s,
  `-inbound` 102,761 m/s, `russell-ocampo-4.3.1-5` 204,547 m/s, `rall-1970-m4-1`
  126,328 m/s — i.e. the same ~4-orders-of-magnitude inflation the S1L1 work
  diagnosed (a single-rev seed for a multi-rev cycler).
- Note `mcconaghy-2006-em-k2` (S1L1) itself shows ERROR in this **generic** batch:
  its strictly-ballistic reproduction required the **hand-supplied rev-correct per-leg
  seeds** (`leg_v_guess`), which the batch does not provide.

## Conclusion

**M7 Phase 1 is solid** — the capability is built and proven on S1L1 (strictly
ballistic at the sourced 200 km floor with rev-correct seeds). **Phase 2
(catalogue-wide measured maintenance bands) is blocked on the same data gap as #388**,
not on M7: the rows do not carry the rev-correct departure velocities / constructed
legs needed to seed the maintenance chain. Supplying a single-rev Lambert seed instead
inflates the TCM by ~4 orders of magnitude, so a naive batch run produces no
trustworthy bands.

This means the dv-band ↔ validation **coupling gate**'s `#423` dependency is landed as
far as it can be without new input: the maintenance-ΔV capability exists and is
verified; broad application is **per-row-seed-data-gated** (the #388 wall — descriptor
acquisition / a family-targeted constructor), exactly like census V0→V3 reproduction.
No amount of re-running the batch changes this; it needs the rev-correct seed data.

Coverage log: `/tmp/m7_phase2_coverage.log` (run 2026-06-23T09:08–09:11 UTC, HEAD
4daa977-era). Read-only; nothing written back to the catalogue.
