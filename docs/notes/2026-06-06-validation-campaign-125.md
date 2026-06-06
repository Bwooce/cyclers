# Validation campaign (#125) — Parts 1-3 outcomes

Date: 2026-06-06. Three-part campaign to move the catalogue beyond "1 of 237
validated". Honest outcomes throughout; golden discipline absolute (the EXPECTED
side of every match check is the row's SOURCED anchor — nothing our own code
computed, nothing loosened to manufacture a match).

## Part 1 — Aldrin inbound V1: PROMOTED

The real-DE440 Aldrin INBOUND cycler (`aldrin-classic-em-k1-inbound`) was built
like-for-like with the outbound twin (same `load_m6b_entries` loader, same
`phase_signature_from_catalogue_entry` + `_resolve_real_t_start` resolution, same
`construct_real_ephemeris_cycler`) and run through
`crosscheck_code_paths`. Both spec §14 V1 halves pass:

- **Path (a) lamberthub** izzo2015 + gooding1990 per-leg agreement:
  max diff `1.23e-08 m/s` (`< V1_TOLERANCE_MPS = 1e-3`).
- **Path (c) Kepler** forward re-propagation residual: pass.

Honest difference from outbound: path (b) (the circular-coplanar resonance
construction) is **unavailable** for inbound (the short ~146 d M->E first leg
yields a heliocentric `(a, e)` that does not map onto the analytic single-ellipse
crossing), not available-but-failing as it is on outbound. The §14 V1 verdict is
unaffected — V1 is defined by the lamberthub-agreement and Kepler-re-propagation
halves, both of which pass; `n_paths_available == 2`, both pass, `agreed == True`.

Evidence test: `tests/verify/test_agreement_lamberthub.py::test_inbound_real_eph_lamberthub_and_kepler_paths_pass`.
Registry: `_LEVEL_EVIDENCE[("aldrin-classic-em-k1-inbound", "V1")]`.
Writeback: catalogue row `validation_level: V0 -> V1`; ratchet
`tests/data/test_schema_v45_fields.py` updated same-commit (V1 set is now the
Aldrin pair).

## Part 2 — descriptor-seeded Russell-12 closure: 0 CLOSE-AND-MATCH

Driver: `scripts/campaign_russell12.py`. For each of the 12 `free_return_arcs[]`
rows: descriptor + trajectory segments -> corrector genome (sequence, per-leg
`(n_revs, branch)`, ToF seeds) -> `ballistic_correct` over a 32-epoch grid
centred on the row's priority date, both residual modes (magnitude AND vector),
16 workers -> compare closed solutions to SOURCED anchors.

Proposed tolerances (documented, applied uniformly, never loosened): V∞ 0.5 km/s,
AR/TR 0.05, transit_times 5 d, period 0.05 yr, corrector convergence floor
0.1 km/s.

Topology mapping: all 12 rows are E-M-E-(E...) multi-arc chains. The descriptor's
first generic arc is the Mars free-return arc, split into the two sourced
transfer legs E->M and M->E (segment ToFs). Subsequent arcs are Earth-Earth
phasing loops (generic -> multi-rev direct E->E; full-rev `M:N` -> M-rev resonant
loop). Longest E->E loop eliminated as the period slack leg (spec §2.1(a)).

### Per-row table (verbatim)

| id | outcome | closing mode | residual km/s | closed mag/vec | seq | E∞ src/ach | M∞ src/ach |
|---|---|---|---|---|---|---|---|
| mcconaghy-2006-em-k2  | CLOSE-OFF-ANCHOR | magnitude | 0.0 | 23/0 | E-M-E-E   | 4.70 / 11.8 | 5.00 / 8.3 |
| russell-ch4-4.991gG2  | CLOSE-OFF-ANCHOR | magnitude | 0.0 | 25/0 | E-M-E-E   | 4.99 / 11.78 | 5.10 / 8.33 |
| russell-ch4-8.049gGf2 | CLOSE-OFF-ANCHOR | magnitude | 0.0 | 17/0 | E-M-E-E   | 8.05 / — | 10.02 / — |
| russell-ch4-9.353Gg2  | CLOSE-OFF-ANCHOR | magnitude | 0.0 | 23/0 | E-M-E-E   | 9.35 / — | 10.52 / — |
| russell-ch4-3.64gGg3  | CLOSE-OFF-ANCHOR | magnitude | 0.0 | 1/0  | E-M-E-E-E | 3.64 / 31.87 | 4.59 / 18.51 |
| russell-ch4-3.78Gg3   | CLOSE-OFF-ANCHOR | magnitude | 0.0 | 15/0 | E-M-E-E-E | 3.78 / 24.07 | 4.63 / 10.98 |
| russell-ch4-5.30gGf3  | CLOSE-OFF-ANCHOR | magnitude | 0.0 | 21/0 | E-M-E-E-E | 5.30 / 23.96 | 9.17 / 10.87 |
| russell-ch4-9.94Gg3   | CLOSE-OFF-ANCHOR | magnitude | 0.0 | 1/0  | E-M-E-E-E | 9.94 / 28.10 | 10.76 / 16.46 |
| russell-ch4-3.66gfF3  | NO-CLOSE         | —         | —   | 0/0  | E-M-E-E-E | 3.66 / — | 4.66 / — |
| russell-ch4-5.30ggF3  | CLOSE-OFF-ANCHOR | magnitude | 0.0 | 28/0 | E-M-E-E-E | 5.30 / 23.96 | 5.44 / 10.87 |
| russell-ch4-5.75ggF3  | CLOSE-OFF-ANCHOR | magnitude | 0.0 | 2/0  | E-M-E-E-E | 5.75 / 15.03 | 9.36 / 6.98 |
| russell-ch4-6.44Gg3   | CLOSE-OFF-ANCHOR | magnitude | 0.0 | 3/0  | E-M-E-E-E | 6.44 / 7.31 | 3.74 / 4.88 |

(Some k=2 rows' achieved values not reprinted in the summary header; full data in
the campaign JSON. The "ach" column is OUR computed V∞ at the closed geometry —
shown only as the off-anchor evidence, never asserted as an EXPECTED.)

### Counts

- **CLOSE-AND-MATCH: 0**
- **CLOSE-OFF-ANCHOR: 11**
- **NO-CLOSE: 1** (`russell-ch4-3.66gfF3`)

### Interpretation (honest)

The corrector reaches the V∞-continuity residual floor (magnitude mode, residual
0.0 km/s) on 11 of 12 rows — the multi-arc *structure* is right and the chain
closes ballistically. But every closed solution lands in the **degenerate
high-V∞ basin**, far off the sourced low-V∞ anchors (E∞ achieved 7-32 km/s vs
sourced 3.6-9.9; transits drift off the sourced symmetric ToFs). This is exactly
the documented S1L1 closure blocker (project memory
`project_s1l1_realeph_closure_blocker.md`): a **family-selection** problem, not
an infrastructure gap. The single-start-per-epoch corrector finds the same
degenerate basin S1L1 floors at; reaching the published family needs the lower
basin the seeding ladder cannot yet seat.

Vector mode (bend-feasibility inside the solve) closed **0** rows: the degenerate
high-V∞ basin is bend-infeasible, so the vector residual never reaches the floor
there — a useful, honest cross-check that the magnitude closures are the
degenerate ones.

The closest approach is `russell-ch4-6.44Gg3` (E∞ off by 0.87, M∞ off by 1.14
km/s) — still outside the 0.5 km/s tolerance. **No promotion**: no row earns
V1-grade multi-arc closure evidence in this pass. The campaign is genuine
first-ever multi-arc validation *evidence* (recorded, auditable), with the honest
verdict that the descriptor-seeded single-start corrector does not select the
sourced family on these E-M multi-arc rows.

## Part 3 — gauntlet tier sweep over all 237 rows

Driver: `scripts/sweep_gauntlet_ledger.py`. For each catalogue row, `run_gauntlet`
is fed the axes its data supports:

- **Axis C (provenance)** — `classify_validation(orbit_source, vinf_source,
  same_fidelity)` + `Corroboration` from the row's tags (every row has these).
- **Axis A (agreement)** — supplied ONLY for the two rows with recorded
  real-closure agreement evidence (the Aldrin pair, `agreed == True` from Parts 1
  and the existing outbound test). All other rows: Axis A unavailable.
- **Axis B / D** — not run / not falsified for any row (no adversarial pass).

Decision consequence (gauntlet rules, `verify/gauntlet.py`): GOLD/SILVER require
machine-confirmation (Axis A available + agreed). Only the Aldrin pair has that;
both are cross_validated (Axis C) -> **GOLD**. Every other row has Axis A
unavailable and no failing axis -> **BRONZE**.

Verdicts are written to `data/gauntlet_ledger.jsonl` (one `LedgerEntry` per row,
`verdict_tier` + `verdict_audit`). Census frozen in
`tests/verify/test_verdict_tier_census.py` (live-ledger census, alongside the
existing decision-matrix census). NO catalogue field changes from this part —
verdicts live in the ledger; `validation_level` changes only via Parts 1-2.

### Verdict census (live, 237 rows)

- **gold: 2** (aldrin outbound + inbound)
- **silver: 0**
- **bronze: 235**
- **rejected: 0**
