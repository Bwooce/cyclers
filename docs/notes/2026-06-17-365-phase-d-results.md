# #365 Phase D — V0→V1 promotion wave results

**Read 2026-06-17 AET.** Phase D of the catalogue admission recovery agent
(deferred from earlier per-candidate same-model corrector runs against
published anchors per §14 V1 like-for-like).

## Scope reduction vs the brief

The brief projected 15 V0→V1 promotion candidates from four digest verdict
notes (Russell-Ocampo 2003 / McConaghy 2004 / McConaghy 2005 / McConaghy 2006).
Closer inspection during Phase D execution found:

* **Russell-Ocampo 2003** (Tables 5-8): 4 candidates, all unique.
* **McConaghy 2004 JSR** (Table 6, S1L1 22-encounter DE405): the S1L1 catalogue
  sink is `mcconaghy-2006-em-k2` (V0) — same physical cycler as
  `russell-ch4-4.991gG2` (already V3 via #167/#94). 1 candidate row.
* **McConaghy 2005 JSR** (Table 2 formal labels): the 2 candidates
  (`russell-ocampo-2.5.1+0`, `russell-ocampo-4.3.1-5`) are DUPLICATES of
  Russell-Ocampo 2003 Tables 5/7. 0 new candidates.
* **McConaghy 2006 JSR** (Tables 2-9 + abstract): 8 vehicles all map to the
  same row `mcconaghy-2006-em-k2`. 0 new candidate rows beyond McConaghy 2004.

**Net unique rows under V1 evaluation in Phase D: 5** —
`russell-ocampo-{2.5.1+0, 3.1.2+1, 4.3.1-5, 4.5.2-2}` + `mcconaghy-2006-em-k2`.

The McConaghy 2006 V3→V4 bonus path is misframed in the brief: the V3 row is
`russell-ch4-4.991gG2`, not `mcconaghy-2006-em-k2` (V0). The actual V3→V4
path requires our DE440 stack to reproduce ONE of the McConaghy 2006 Tables
2-9 DE405 itineraries with documented tolerance — a substantial real-eph
reproduction effort, not in Phase D scope this pass.

## Method (§14 V1 like-for-like)

For each row, run the existing `closer_sweep_v1` substrate:

1. Seed the single circular-coplanar free-return ellipse from the row's
   SOURCED `(aphelion_au, transit_em_days)` (Russell 2004 Table 3.4
   / McConaghy 2006 abstract).
2. Dense phase scan (`PHASE_EPOCHS=4096`) for the best-residual `t0`.
3. Run the corrector. Pass = converges + EMERGED V_inf within 0.5 km/s of
   the sourced E/M anchor + §14 V1 mechanics gate (closed,
   V_inf-continuous single ellipse).

Driver: `scripts/run_365_russell_ocampo_v1.py` (4 candidates) +
`scripts/run_365_mcconaghy_2006_em_k2_v1.py` (1 candidate).
Verdicts: `data/<row-id>_v1_verdict.jsonl` per candidate.

## Per-candidate verdict

| Row | Paper cycler | Sourced E / M | Derived E / M | |ΔE|, |ΔM| | V1 gate | Verdict |
|---|---|---|---|---|---|---|
| `russell-ocampo-2.5.1+0` | Cycler-2-5-1-3 (RO Table 5) | 7.8 / 9.9 | 7.895 / 9.942 | 0.10, 0.04 | PASS | **V0→V1** |
| `russell-ocampo-3.1.2+1` | Cycler-3-1-2-11 (RO Table 6) | 3.4 / 4.6 | 3.416 / 4.616 | 0.02, 0.02 | FAIL vinf-continuity | V0 (CLOSE-NOT-V1) |
| `russell-ocampo-4.3.1-5` | Cycler-4-3-1-20 (RO Table 7) | 3.1 / 2.5 | (no close) | — | FAIL no-close | V0 (NO-CLOSE) |
| `russell-ocampo-4.5.2-2` | Cycler-4-5-2-12 (RO Table 8) | 3.4 / 4.6 | 3.361 / 4.608 | 0.04, 0.01 | FAIL vinf-continuity | V0 (CLOSE-NOT-V1) |
| `mcconaghy-2006-em-k2` | S1L1 ballistic | 4.7 / 5.0 | 4.771 / 5.036 | 0.07, 0.04 | FAIL vinf-continuity | V0 (CLOSE-NOT-V1) |

**1 PASS (V0→V1), 4 honest negatives (stay V0).**

## Why so many honest negatives — topology, not infrastructure

Four of the five rows fail the §14 V1 like-for-like at the SAME gate: the
EMERGED V_inf magnitudes are within 0.5 km/s of the sourced anchor, but the
single circular-coplanar free-return ellipse does NOT close to Earth — the
reconstructed return leg breaks Mars V_inf continuity by tens of km/s.

This is the multi-arc topology signature that the
existing `closer_sweep_v1` `CLOSE_NOT_V1` set already pins
(`russell-ocampo-4.7.1-2`, `russell-ocampo-4.8.1+2`). The
McConaghy-Russell-Longuski 2005 Table 2 formal labels make the topology
explicit:

* S1L1 = `2g(2.8277, 657.97°, U) g(1.4508, 522.29°, L)` — TWO generic
  Earth-Earth arcs joined at the Mars flyby.
* Cycler-3-1-2-11 / Cycler-4-5-2-12: multi-arc by Russell-Ocampo's flyby-count
  topology (3+ Earth flybys per cycle ⇒ 2+ intermediate Earth-Earth loops).

A single radial-crossing ellipse FUNDAMENTALLY cannot represent these. The
honest negatives at the V1 tier are CORRECT — these rows' real elevated-tier
evidence lives elsewhere:

* `mcconaghy-2006-em-k2`'s sibling row `russell-ch4-4.991gG2` is V3 via
  #167/#94's REBOUND/IAS15 real-eph closure of the corrected two-arc topology.
* The three Russell-Ocampo multi-arc rows would require the multi-arc
  closure infrastructure to reach V1 — a documented out-of-scope follow-on.

The PASS row `russell-ocampo-2.5.1+0` clears the gate because its 94-day
Earth-Mars transit + aphelion ratio 1.44 admit a single radial-crossing
ellipse that closes to Earth (high-V_inf short-transit, single-ellipse-
reachable family).

## Discipline note

Per `feedback_orbit_closure_discipline`: the math gives the verdict, not a
tuned pass. No corrector tolerance was loosened to force any candidate
through. The 4 honest negatives are documented FAIL evidence with explicit
fail modes (`CLOSE-NOT-V1` vs `NO-CLOSE`), pinned by frozen-gate pytests
(`tests/verify/test_365_russell_ocampo_v1_promotion.py` +
`tests/verify/test_365_mcconaghy_2006_em_k2_v1.py`) so any future "tuned
pass" would itself break the test (the FAIL is the registry evidence at
this row, at this tier).

## Catalogue census

* V4: 1 (unchanged — SILVER `umbriel-oberon-1-1-uranian-quasi-cycler-2026`)
* V3: 2 (unchanged — `russell-ch4-4.991gG2` S1L1, `russell-ch4-8.049gGf2`)
* V2: 6 (unchanged — Aldrin outbound + 5 Ross-RT)
* V1: 21 → **22** (+1: `russell-ocampo-2.5.1+0`)
* V0: 248 → 247 (−1)

## Pointers

* Driver scripts: `scripts/run_365_russell_ocampo_v1.py`,
  `scripts/run_365_mcconaghy_2006_em_k2_v1.py`.
* Verdict JSONLs: `data/russell-ocampo-*_v1_verdict.jsonl` (×4),
  `data/mcconaghy-2006-em-k2_v1_verdict.jsonl`.
* Frozen-gate tests: `tests/verify/test_365_russell_ocampo_v1_promotion.py`,
  `tests/verify/test_365_mcconaghy_2006_em_k2_v1.py`.
* Registry: `src/cyclerfinder/data/validate.py::_LEVEL_EVIDENCE` (new
  `russell-ocampo-2.5.1+0 → V1` entry).
* Catalogue ratchet: `tests/data/test_schema_v45_fields.py`
  (`russell-ocampo-2.5.1+0 → V1`).
* Source digests: `docs/notes/2026-06-17-digest-russell-ocampo-2003.md`,
  `docs/notes/2026-06-17-digest-mcconaghy-2004.md`,
  `docs/notes/2026-06-17-digest-mcconaghy-2005.md`,
  `docs/notes/2026-06-17-digest-mcconaghy-2006.md`,
  `docs/notes/2026-06-17-366-s1l1-cycler-confusion-correction.md`.
