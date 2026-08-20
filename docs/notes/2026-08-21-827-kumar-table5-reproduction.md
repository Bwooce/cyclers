# `#827`: Kumar-Rawat-Rosengren-Ross (2026) Table-5 3:1 <-> 2:1 heteroclinic reproduction — 7/7 CLEAN

**Task:** `#827`, registered 2026-08-11 (found during `#822`, dispatched 2026-08-12; this run
2026-08-20/21 by a fresh agent after a prior dispatch was killed mid-run by an unrelated process
restart, resuming its already-working, uncommitted machinery). Digit-grade reproduction of

    B. Kumar, A. Rawat, A. J. Rosengren, S. D. Ross (2026). "Cislunar Resonant Transport and
    Heteroclinic Pathways: From 3:1 to 2:1 to L1," *Advances in Space Research* 77(3):3815-3845,
    DOI 10.1016/j.asr.2025.12.005 (= arXiv:2509.12675v2; corpus
    `kumar-2025-arxiv-2509.12675.pdf`; digest `docs/notes/2026-06-20-digest-kumar-2025.md`)

Table 5's (Appendix 8.2) seven printed `Wu(3:1) ∩ Ws(2:1)` manifold-intersection states, at
Earth-Moon `mu = 1.2150584270572e-2` (their own printed mass ratio, Section 2.1), built on
`#822`'s Poincare-section Newton connection machinery
(`cyclerfinder.search.vaquero_em_cycler_connections`, reused unchanged).

**Module:** `src/cyclerfinder/search/kumar_em_resonant_heteroclinics.py` (see its own docstring
for the full method account — targeted seeding from the printed state, the perigee-section match
gate, and the `KUMAR_EPSILON` re-calibration for the 2:1's extreme saddle).
**Driver:** `scripts/screen_827_kumar_table5_reproduction.py`.
**Data:** `data/found/827_kumar_table5_reproduction/results.json` (checkpointed per Jacobi
constant; every row carries the full converged connection, verification battery, and
component-wise diff against the print).
**Evidence tests:** `tests/search/test_kumar_em_resonant_heteroclinics.py` (10 tests, NOT slow —
table/constant structure, Table-6 node re-derivation against the print, and two from-scratch
connection reconstructions seeded only by the recorded phase indices, ~29s total).

---

## Result: all seven printed rows reproduce, digit-grade

| C (Jacobi) | transfer type | matched | match_distance | runner_up_distance | elapsed_s |
|---|---|---|---|---|---|
| 2.54 | 1 (short) | **true** | 7.611e-07 | 1.307e+00 | 123.1 |
| 2.70 | 1 (short) | **true** | 6.819e-05 | 3.932e-01 | 136.7 |
| 2.86 | 1 (short) | **true** | 1.283e-06 | 3.589e-01 | 77.3 |
| 3.00 | 1 (short) | **true** | 1.193e-07 | 3.470e-01 | 301.3 |
| 3.05 | 1 (short) | **true** | 1.352e-06 | 4.757e-02 | 63.9 |
| 3.10 | 2 (long, via 5:2) | **true** | 1.745e-06 | 1.635e+00 | 87.9 |
| 3.15 | 2 (long, via 5:2) | **true** | 4.485e-07 | 1.587e+00 | 62.5 |

7/7 — no honest negatives to report. C=3.00 was the prior agent's single completed row (301.3 s,
before the restart killed it); the other six were run this session, one Jacobi constant per
foreground `uv run` invocation (63-137 s each, ~9.5 min total wall for the six — faster in
practice than the ~5 min/point the module docstring's `#822`-lineage estimate implied).

**Honest scope note on the C=2.70 row:** its match distance (6.819e-05) is the tightest of the
seven — still comfortably under `KUMAR_MATCH_TOL = 1e-4` (a ~1.5x margin, vs 50x-800x for the
other six), but the smallest safety factor in the set. The runner-up separation (3.932e-01) is
still four orders of magnitude larger than the match distance, so this is unambiguously a
specific-point identification, not a proximity coincidence — but it is the one row where a
reader auditing this table should look twice. No gate was loosened to pass it; the achieved
distance is recorded exactly as computed.

Every row's `matched_perigee_state` and `component_diffs` (x, y, xdot, ydot) are in
`results.json` alongside the full converged `HeteroclinicConnection` and `ConnectionEvidence`
records (Newton residual, ghost-guard margins, independent-Radau crosscheck, forward/backward
re-approach, Jacobi drift) — the complete `#822` verification battery passed at all seven points,
not just the perigee-match gate.

## What changed this session

- Fixed the two pre-existing `ruff check`/`ruff format` violations (E501 line-too-long in both
  files) left by the prior killed run — no logic changes.
- Ran the six not-yet-recorded Jacobi constants (`{2.54, 2.70, 2.86, 3.05, 3.10, 3.15}`); merged
  into the same checkpointed `results.json` the prior C=3.00 row already lived in (the driver's
  own per-C read-merge-write convention, unchanged).
- Added `tests/search/test_kumar_em_resonant_heteroclinics.py` (the module had zero tests before
  this session).
- **Record-integrity fixes caught on advisor review**: `KUMAR_REPRODUCTION_CS` still listed only
  the four Table-1 values after all seven had matched — anyone re-running the driver against an
  empty `results.json` would have reproduced 4 rows, not 7. Widened it to all seven and tightened
  the corresponding test from a subset check to an equality check. Also corrected two now-
  falsified docstring claims: Table 6 DOES print both parent ICs at 2.70/2.86 (this session's own
  successful runs at those C values are direct proof), and C=2.54 is NOT outside Vaquero's
  `[2.54, 2.66]` band — it is exactly that band's lower edge (see the `#854` section below).

## Literature-novelty gate

**Not re-run live this session** — per the `#827` dispatch's own guidance and the
compute/adjudicate split this project uses for reproduction tasks: this module carries no
catalogue writeback (no row is added, edited, or promoted by this work), and the object being
reproduced is Kumar et al.'s own printed table from their own paper — the least "novel" claim
this project can make. `#822` (the sibling task building the same connection machinery, dispatched
2026-08-11) already ran `search/literature_check.py`'s live-WebSearch mandatory floor against this
exact paper's heteroclinic-connection concept and returned verdict **`published`**, confidence
0.95, anchored on the Casoliva-2010 lineage (`10.2514/1.46856`). Nothing in this task's scope
introduces a new claim beyond that anchor — every number reported here is a digit-grade match to
a state Kumar et al. themselves printed. **Nothing here is novel and nothing here is claimed
novel.**

## `#839` gate status — C=3.13 is NOT covered by this run

`#839` is gated on `#827` because `vaquero-31-c313-em-resonant-po-2013` sits at **C=3.13**, which
`#839`'s own registration describes as "inside" Kumar et al.'s published `C_J ∈ [3.00, 3.15]`
heteroclinic band. **None of the seven Table-5 rows reproduced here is at C=3.13** — the printed
table only tabulates `{2.54, 2.70, 2.86, 3.00, 3.05, 3.10, 3.15}`, and
`keh.kumar_table5_state6(3.13)` correctly raises `KeyError` (asserted directly by
`test_kumar_table5_state6_rejects_a_c_not_in_the_printed_table`) — there is no printed
intersection state at that Jacobi constant to reproduce. **`#839` still needs its own targeted
run** — not a table lookup but a fresh search for a `Wu(3:1) ∩ Ws(2:1)` (or whichever direction
is relevant) connection AT C=3.13 specifically, using this module's machinery
(`build_kumar_node` + `#822`'s connection search) but with no printed target to seed against or
match digit-grade. That is explicitly out of `#827`'s scope (a reproduction of PRINTED digits,
not a new search at an untabulated C) and is left for `#839`'s own dispatch.

## C=2.54 touches a catalogued row — new follow-up `#854` registered

`KUMAR_TABLE5_31_TO_21` covers three C values (2.54, 2.70, 2.86) that sit inside or at the edge
of `#822`'s own Vaquero overlap band `[2.54, 2.66]` — the module docstring originally (wrongly)
claimed the whole target set was "OUTSIDE" that band; corrected in this session's edit.
**C=2.54 is exactly Vaquero's own lower band edge, and `data/catalogue.yaml` carries a 3:1 row
there: `vaquero-31-c254-em-cycler-2013`.** `#828`'s adjudication (`docs/notes/2026-08-12-828-
vaquero-connection-tier-adjudication.md`, Sec. 1) already identifies this row as one endpoint of
`#822`'s own C=2.54 connection evidence (self-consistency-only, at this project's mu) — but this
`#827` C=2.54 result is a **different kind of evidence on the same node**: Kumar's own printed
mu, Kumar's own Wu(3:1)->Ws(2:1) direction, and a **digit-grade match to a state Kumar et al.
themselves published** (match_distance 7.6e-07), not merely a self-consistent Newton closure.
`#828`'s annotation does not cover this (it predates `#827`'s C=2.54 run and is scoped to `#822`'s
own connection only). Per `#828` Sec. 2's own ruling this is still **not a tier-promotion
question** — but it is new, specific evidence a future adjudicator would want, and it directly
overlaps `#840`'s registered scope (round-tripping `#822`'s two in-band annotations at the
C=2.54/2.66 band edges) for the C=2.54 edge specifically. **Registered as `#854`**: adjudicate
whether this `#827` C=2.54 digit-grade connection warrants its own comment-only `ADDED EVIDENCE`
annotation on `vaquero-31-c254-em-cycler-2013` (parallel to `#828`'s two), and whether it narrows
`#840`'s remaining scope to just the C=2.66 edge. No catalogue writeback performed here — that
adjudication is `#854`'s, not `#827`'s.

## Verification

`tests/search/test_kumar_em_resonant_heteroclinics.py` (10 tests, ~29s), `ruff check .` / `ruff
format --check .` (repo-wide, clean), full `uv run mypy src tests` (clean, 853 source files). No
`data/catalogue.yaml` change in this task, so the full `tests/data tests/search -q` ratchet was
not required by the dispatch; `uv run pytest tests/search -q` (the broader sanity pass) was run
separately — see that invocation's own result for pass/fail detail. See the `#827` commits in
`git log`.
