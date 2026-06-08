# Closer sweep — V1 candidates from the 2026-06-08 ingest (Russell T3.4 + Rall)

Date: 2026-06-08
Code version: `3b40d23`
Runlog: `data/runs/closer-sweep-2026-06-08.jsonl` (31 records)
Test gate: `tests/search/test_closer_sweep_v1.py`

## Purpose

Determine which catalogue rows close in their **defining circular-coplanar model**
(the V1-candidate set), to feed the continuation campaign. This sweep covers the
rows ingested on 2026-06-08 that had **not** been run through the #137 free-return
(radial-crossing) closer:

- 16 Russell 2004 Table 3.4 cyclers (`russell-ocampo-*`, commit `5dd4742`)
- 15 Rall 1970 (MIT TE-34) free-fall periodic orbits (`rall-1970-*`, commit `3ebbb7a`)

The original 12 descriptor (`free_return_arcs[]`-bearing) rows are already done
(8 CLOSE-AND-MATCH, 4 already V1; `test_russell12_likeforlike_probe.py`).

## Method (like-for-like, golden-firewalled)

For each row carrying a SOURCED `orbit_elements.aphelion_au` **and**
`invariants.transit_times_days`:

1. Seed the transfer ellipse `(a, e)` from the SOURCED aphelion + outbound transit
   (`seed_ae_from_aphelion_transit`) — these are CONSTRAINTS (imposed).
2. Scan `t0` over one period (4096-point dense floor) for the phase minimising the
   free-return residual at the seeded `(a, e)`.
3. Close the single radial-crossing ellipse (`free_return_correct`,
   `Ephemeris("circular")`, `tol_kms=0.1`).
4. Compare the EMERGED per-body V∞ against the row's SOURCED V∞ anchor
   (EVIDENCE — never imposed; tol 0.5 km/s, not loosened per row).
5. For closed+matched rows, run the §14 V1 mechanics
   (`free_return_v1_mechanics`): lamberthub re-solve (path a) + Kepler
   re-propagation (path c) + Mars-flyby V∞-continuity honesty gate.

EXPECTED side of every match = the row's SOURCED V∞ anchor (Russell 2004 Table
3.4). The closer's emerged V∞ is evidence, compared non-circularly.

## Reachability

| Set | Reachable? | Why |
|---|---|---|
| 16 Russell T3.4 `russell-ocampo-*` | YES | carry SOURCED `aphelion_au` (from AR×1.52 AU) + single `transit_times_days` + 2-body (E,M) V∞ anchors |
| 15 Rall `rall-1970-*` | **NO** | Rall Model I.B does not tabulate per-arc `aphelion_au` (`aphelion_au = None`); closer has nothing to seed `(a, e)` from |

## Per-row outcomes

### Russell 2004 Table 3.4 (16 rows)

| id | aphelion_au | transit_d | seed (a, e) | outcome | res km/s | emerged V∞ E/M | sourced E/M | V1 |
|---|---|---|---|---|---|---|---|---|
| russell-ocampo-3.1.1+2 | 2.174 | 115 | 1.583, 0.373 | CLOSE-AND-MATCH | 0.0089 | 5.49 / 9.23 | 5.4 / 9.2 | **PASS** |
| russell-ocampo-3.1.3+0 | 2.174 | 123 | 1.587, 0.370 | CLOSE-AND-MATCH | 0.0025 | 5.13 / 9.15 | 5.1 / 9.1 | **PASS** |
| russell-ocampo-4.1.1-4 | 2.189 | 137 | 1.595, 0.373 | CLOSE-AND-MATCH | 0.0048 | 5.11 / 9.22 | 5.5 / 9.3 | **PASS** |
| russell-ocampo-4.1.2-2 | 2.174 | 132 | 1.587, 0.370 | CLOSE-AND-MATCH | 0.0106 | 5.08 / 9.14 | 5.2 / 9.2 | **PASS** |
| russell-ocampo-4.1.4-1 | 2.174 | 129 | 1.587, 0.370 | CLOSE-AND-MATCH | 0.0106 | 5.08 / 9.14 | 5.1 / 9.2 | **PASS** |
| russell-ocampo-4.6.3+0 | 2.174 | 105 | 1.573, 0.382 | CLOSE-AND-MATCH | 0.0030 | 6.42 / 9.45 | 6.4 / 9.5 | **PASS** |
| russell-ocampo-4.7.1-2 | 2.69  | 120 | 1.845, 0.458 | CLOSE-MATCH-NO-V1 | 0.0347 | 6.18 / 11.26 | 6.6 / 11.4 | fail (13.9 km/s Mars V∞ break) |
| russell-ocampo-4.8.1+2 | 1.991 | 76  | 1.371, 0.452 | CLOSE-MATCH-NO-V1 | 0.0011 | 12.72 / 10.71 | 12.5 / 10.7 | fail (11.9 km/s Mars V∞ break) |
| russell-ocampo-3.5.1+2 | 1.429 | 231 | seed→e-floor | NO-CLOSE | — | — | 2.7 / 1.5 | — |
| russell-ocampo-4.1.1-6 | 1.429 | 256 | seed→e-floor | NO-CLOSE | — | — | 2.7 / 1.6 | — |
| russell-ocampo-4.1.2-3 | 1.429 | 250 | seed→e-floor | NO-CLOSE | — | — | 2.6 / 1.5 | — |
| russell-ocampo-4.6.1-4 | 1.383 | 154 | seed→e-floor | NO-CLOSE | — | — | 6.8 / 2.1 | — |
| russell-ocampo-4.8.1+3 | 1.459 | 164 | seed→e-floor | NO-CLOSE | — | — | 7.7 / 3.1 | — |
| russell-ocampo-4.9.1-3 | 1.429 | 256 | seed→e-floor | NO-CLOSE | — | — | 2.7 / 1.6 | — |
| russell-ocampo-4.10.1-3 | 1.398 | 263 | seed→e-floor | NO-CLOSE | — | — | 10.2 / 3.6 | — |
| russell-ocampo-4.12.1-2 | 1.474 | 268 | seed→e-floor | NO-CLOSE | — | — | 11.6 / 4.8 | — |

The 8 NO-CLOSE rows all have a low aphelion ratio (aphelion ≈ 1.38–1.47 AU) paired
with a long transit (154–268 d). The aphelion+transit seed bisection drives `e`
to its 0.05 floor where the ellipse no longer reaches Mars — a single
radial-crossing ellipse cannot represent these geometries (they are multi-arc /
different-eccentricity solutions). This is a GEOMETRY fact, not a seeding gap.

The 2 CLOSE-MATCH-NO-V1 rows close and the emerged V∞ matches the sourced anchor,
but the reconstructed single ellipse does **not** close to Earth: the §14 V1
V∞-continuity gate finds a ~12–14 km/s Mars-flyby discontinuity (the return needs
intermediate phasing loops — multi-arc). The honesty gate correctly withholds V1.

### Rall 1970 (15 rows) — all NOT-REACHABLE

`rall-1970-{m4-1, m4-1a, m6-1, m6-2, m6-3, m5-1a..e, m5-2a..e}`. Every row is the
circular-coplanar Model I.B with `aphelion_au = None` (per-arc `(a, e)` not
tabulated; Rall App E gives only encounter dates, speeds, and passing distances).
They are also multi-arc E-M-E with 3 distinct Earth V∞ values (asymmetric), so a
single symmetric free-return ellipse could not represent them even if seedable.

## Counts

- **NEW V1-promotable rows: 6** — `russell-ocampo-3.1.1+2`, `russell-ocampo-3.1.3+0`,
  `russell-ocampo-4.1.1-4`, `russell-ocampo-4.1.2-2`, `russell-ocampo-4.1.4-1`,
  `russell-ocampo-4.6.3+0`.
- CLOSE-AND-MATCH but not V1 (multi-arc continuity break): 2 —
  `russell-ocampo-4.7.1-2`, `russell-ocampo-4.8.1+2`.
- NO-CLOSE (single ellipse cannot reach Mars): 8 (Russell T3.4, listed above).
- **NOT-REACHABLE: 15** — all Rall rows (no per-arc aphelion to seed from).

## Recommended validation_level writebacks (for the main session to apply)

Promote to **V1** (closed, V∞-continuous single free-return arc; §14 V1 mechanics
pass; emerged V∞ matches the sourced anchor — pinned by
`tests/search/test_closer_sweep_v1.py::test_new_row_closes_matches_and_passes_v1`):

```
russell-ocampo-3.1.1+2   -> V1
russell-ocampo-3.1.3+0   -> V1
russell-ocampo-4.1.1-4   -> V1
russell-ocampo-4.1.2-2   -> V1
russell-ocampo-4.1.4-1   -> V1
russell-ocampo-4.6.3+0   -> V1
```

No other writebacks recommended. The 2 CLOSE-MATCH-NO-V1 and 8 NO-CLOSE rows stay
at their current level (V0/None); the V1 honesty gate did not pass.

## Continuation-campaign implications

- The 6 V1-promotable rows are well-conditioned single-ellipse seeds — directly
  usable to start a continuation campaign.
- The 2 CLOSE-MATCH-NO-V1 and 8 NO-CLOSE rows are genuine multi-arc geometries; a
  single-ellipse seed cannot start them. They need a multi-arc (phasing-loop)
  closer to continue.
- The 15 NOT-REACHABLE Rall rows cannot be seeded by the closer at all until their
  per-arc `(a, e)` / aphelion state is recovered — these are the rows that need
  **Appendix C** state before the continuation campaign can touch them.
