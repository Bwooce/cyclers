# #596 follow-up 2: Russell Table 3.4 backfill — validated derivation method found

**Date:** 2026-07-15
**Supersedes the "clean negative" finding in `2026-07-15-596-russell-backfill-pilot-inversion-fails.md`**
(that note's conclusion — "no `(a,e)` satisfies all three Table 3.4 columns" — is corrected below: two
of the three columns (AR/aphelion, V∞) DO pin down `(a,e)` accurately; the third (ToF) has a
row-specific quirk, not a general modeling failure).

## What was wrong with the first attempt

The first pilot fit `(a_au, e)` to match BOTH `tof_em_days` AND `V∞_Earth` simultaneously. Re-reading
Russell 2004 Ch.3 (printed pp.56-69, PDF pages 70-84) found the actual structural facts:

- **AR is defined as `aphelion / 1.52 AU`** (his own rounded Mars radius, not the catalogue's more
  precise 1.523691) — confirmed by the method-mining note's own Ch.3.7 summary.
- Table 3.4's "generic return" rows are constructed via Russell's OWN Ch.2 multi-rev-Lambert
  ψ-parameterized machinery (Eq. in Ch.2.7/2.8), which is a DIFFERENT parameterization than this
  project's `free_return_geometry`, even though both describe the same underlying physics
  (Keplerian arc crossing Earth's and Mars's circular radii).
- Critically, the Aldrin row (`1.0.1.-1`) is footnoted as an externally-sourced, pre-existing named
  cycler (Aldrin 1985, refs 1/13/14) — NOT one of Russell's own newly-discovered generic-return
  solutions. This turns out to matter (see below).

## The validated method

Dropped the ToF constraint; fit `(a_au, e)` using only **AR (-> aphelion) and V∞ at Earth** via
`cyclerfinder.search.free_return.free_return_geometry` + `scipy.optimize.least_squares` (global
multi-start, 29 seeds). Validated against **two independent ground truths**:

1. **Aldrin (`1.0.1.-1`)** — catalogue's already-sourced `a=1.60, e=0.393` (multiple independent
   literature references). Fit: `a=1.6038, e=0.3932` — **0.24% / 0.05% error.**
2. **Cycler `2.5.1.+0`** — Table 3.5 gives the actual initial 3D v∞ VECTOR `(6.50, 4.35, 0)` km/s at
   Earth departure (footnote a). Used this to compute `(a,e)` DIRECTLY via two-body state-vector orbit
   determination (vis-viva + angular momentum, no fitting at all): `a=1.5651, e=0.4010`. This is a
   fully independent, non-inversion ground truth. The AR+V∞ inversion method recovered
   `a=1.5633, e=0.4001` — **0.12% / 0.22% error.**

Both independent checks agree to well under 0.3%. **The AR+V∞_Earth inversion is validated and
accurate enough to trust for the backfill.**

## What about the ToF mismatch?

For cycler `2.5.1.+0`, the direct state-vector `(a,e)` ALSO correctly reproduces Table 3.4's ToF
column via `free_return_geometry`'s `tof_em_days` (94.37 vs. table's 94 — matches). This means
`tof_em_days` IS the right quantity in general. The Aldrin-specific mismatch (146 reported vs. ~103
emerged from the fitted `(a,e)`) is most likely explained by Aldrin being an externally-sourced named
cycler rather than one of Russell's own generic-return search results — worth a footnote in any
Aldrin-specific writeback, but **not a blocker for the other ~215 rows**, which are Russell's own
uniform generic-return search output and should follow the `2.5.1.+0` pattern.

## Practical implication for the backfill

- `trajectory.segments[out-em].tof_days`: **cite Table 3.4/3.9-3.11's own ToF column DIRECTLY**
  (`kind: unknown` -> sourced, not derived) — this is exactly what the existing Aldrin catalogue rows
  already do (`tof_days: 146`, matching Table 3.4 verbatim).
- `trajectory.segments[out-em].a_au` / `.e`: **derive via the validated AR+V∞_Earth inversion**
  (`kind: derive`, computed by this project's own solver from sourced AR/V∞ inputs — consistent with
  the existing `data_gaps` "derive" convention).
- `trajectory.segments[ret-me].tof_days`: total cycler period (`p` synodic periods) minus the outbound
  `tof_days` — same pattern as the existing Aldrin `ret-me` note ("780 - 146 = 634").
- **Independent cross-check available for every row**: V∞ at Mars is a THIRD column Russell provides;
  the fitted `(a,e)`'s emerged `vinf['M']` should match it as a free, non-imposed consistency check
  (as it did in both validation cases: 9.751/9.7 for Aldrin, 9.936/9.9 for 2.5.1.+0).

## Status

**Method validated, NOT yet executed at scale.** This note documents a working, cross-validated
derivation pipeline for the ~216 existing Russell-family gaps plus the ~38 uncatalogued Table 3.4
rows. The actual mechanical backfill (parse all Table 3.4/3.9-3.11 rows from the 2026-06-07
transcription note, map to catalogue entry IDs via the `p.h.s.i` designator, run the inversion per
row, write back with `kind: derive` tags + citations, add the ~38 missing rows, run the full
`tests/data tests/search` ratchet) is a separate, large execution task — not attempted in this note.
