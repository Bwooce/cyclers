# #859 — Resonant Atlas pilot, Stage A: harness built + smoke-tested

**Date**: 2026-08-21
**Scope**: build and smoke-test only, per this task's own explicit dispatch note. The
full ~50-100 CPU-hour / ~2-3 day Stage A sweep was NOT run — it is far beyond a single
foreground dispatched-agent session, and this project's own standing lesson
(`feedback_subagent_background_is_fatal`) forbids a subagent self-backgrounding a
long-running process. The harness is handed to the coordinating session, which can
launch and monitor a genuine multi-day background run.

Sources read in full before building anything: `#859`'s own registration
(`data/OUTSTANDING.md`), the full `#858` review
(`docs/notes/2026-08-21-858-campaigns-789-790-791-review.md`, esp. Sec. 3.4/6/7),
`src/cyclerfinder/search/campaign_runner.py` (`#788`), `neptune_triton_resonant_families.py`
(`#776`/`#777`, the proven family-continuation reference), `jovian_resonant_families.py`
(the system-agnostic seed/corrector/classifier machinery every family module reuses),
`cr3bp_continuation.py` (`continue_family`'s own gauntlet), and `core/satellites.py` (the
moon registry).

## 1. Final system list + literature/corpus pre-check

`#859`'s registered list (Uranus-Oberon positive control, Jupiter-Ganymede,
Uranus-Titania, Saturn-Rhea band-edge probe) **survived the pre-check unchanged**. All
four are instantiable with zero new code via `cr3bp.cr3bp_system(primary, secondary)`
(registry mu, all four pairs present in `core/satellites.py`'s `SATELLITES`/`PRIMARIES`).

The mandatory check was NOT run through `search/literature_check.py`'s
`check_literature()` — that function needs a concrete `CandidateSignature` (a specific
converged connection), which Stage A does not produce (it is a pre-discovery triage
stage). Instead, per `#858` Sec. 3.4's own framing ("mandatory per-system
`literature_check.py` + corpus grep"), this was a **system-level** clearance: live
WebSearch queries per system's own resonant-orbit-family/homoclinic-connection concept,
plus a `docs/notes/CORPUS_INDEX.md` grep. Results:

- **Uranus-Oberon** (positive control, kept regardless of novelty by design): CONFIRMED
  published. Anderson & Kumar 2024, AAS 24-288 (already corpused, `#728`'s digest)
  computes exactly this system's planar CR3BP unstable resonant-orbit/manifold/
  heteroclinic survey (3:4/4:5/5:6 exterior, 4:3/5:4/6:5 interior). This is the intended
  known-answer check.
- **Jupiter-Ganymede**: the CCR4BP (4-body, Europa-perturbed) resonant-orbit/secondary-
  resonance literature is dense (Kumar et al. 2021/2023, already corpused via `#688`/
  `#727`/`#728`) — but that is a DIFFERENT model (Jupiter-Ganymede-Europa N=4), not the
  plain 2-body Jupiter-Ganymede PCR3BP this module builds. No direct hit for a pure
  Jupiter-Ganymede CR3BP p:q-resonance/homoclinic-connection census. Defensible.
- **Uranus-Titania**: the corpused Anderson-Kumar 2024 paper studies Titania only as
  Oberon's CCR4BP PERTURBER (secondary resonances on the Oberon family), and explicitly
  proposes repeating its OWN Oberon-centered study for Umbriel/Ariel — not Titania as the
  primary secondary body. A separate hit ("Mapping and maneuvering long-term natural
  orbits around Titania", arXiv:2307.06570) is frozen/natural-orbit design, a different
  orbit class (not resonant p:q families). No direct hit for a Uranus-Titania PCR3BP
  resonant-family/homoclinic census. Defensible.
- **Saturn-Rhea**: a 2026 tour-design paper (Pozzi, Pontani, Beolchi, Susanto & Fantino,
  "Low-Energy and Low-Thrust Exploration Tour of Saturnian Moons", arXiv:2603.07085) uses
  a J2-perturbed CR3BP with HALO orbits (3D, L1/L2) as staging points for heteroclinic/
  homoclinic loops across Rhea/Dione/Tethys/Enceladus/Mimas — a different orbit class
  (spatial halo, not planar p:q resonant families) and a different question (multi-moon
  tour design, not a single-pair resonant-family census). Worth citing as related
  territory but does not preempt Stage A's own planar resonant-family scope. Defensible,
  and matches its intended role as the deliberate low-mass-ratio band-edge probe (no
  resonant-family-specific literature found at all).

No system was swapped. This is a system-level web/corpus clearance, not the full
candidate-signature novelty gate `check_literature()` runs on an actual discovered
connection — Stage B (if dispatched) must still run the real gate on any specific
in-band candidate it converges, per this project's standing
`feedback_literature_novelty_check_baseline` discipline.

## 2. Harness design

New module `src/cyclerfinder/search/resonant_atlas_stage_a.py`:

- `coprime_pairs(max_pq)`: every `(p, q)` in `[1, max_pq]` with `gcd(p, q) == 1`.
  43 pairs at `max_pq=8`.
- `build_stage_a_cells(systems, max_pq, n_c_steps, d_jacobi, x0_sign)`: one JSON-safe dict
  cell per `(system, p, q)`, `x0_sign=-1` only (the `+1` seed sits exactly at the
  secondary's own orbital radius — the documented DOP853 step-collapse hazard already
  characterized in `jovian_resonant_families.survey_candidates`'s own docstring and
  confirmed to hang in `#776`'s own scratchpad check; never surveyed by this harness).
- `stage_a_worker(cell)`: for one `(system, p, q)` — build the two-body seed
  (`jovian_resonant_families.two_body_resonant_seed`), converge it at its own natural
  Jacobi constant (`cr3bp_periodic.correct_symmetric_fixed_jacobi`, `half_crossings=None`
  auto-detected), recover the crossing index the corrector used (needed because
  `continue_family` requires an explicit index to hold the branch fixed), then continue
  the family across `n_c_steps` Jacobi-constant steps (`cr3bp_continuation.continue_family`
  — the SAME gauntlet-validated machinery `#776`/`#777`/`#781` used: convergence, period
  bounds, equilibrium gate, Jacobi conservation, independent-Radau cross-check, fold/
  topology-jump detection). Classifies every gauntlet-passing member's Barden `|lambda|`
  against the `[50, 2500]` in-band window. `"hit"` iff >=1 member lands in-band, `"miss"`
  otherwise (including two clean, EXPECTED negative sub-cases — see below — never an
  `"error"`), `"error"` reserved for genuinely unanticipated exceptions.

Driver `scripts/run_859_resonant_atlas_stage_a.py`: argparse (`--systems`, `--max-pq`,
`--n-c-steps`, `--d-jacobi`, `--n-workers`, `--checkpoint-batch-size`, `--max-batches`,
`--pause-seconds-per-batch`, `--thermal-backoff-seconds`, `--out-dir`, `--report`), wired
straight through `campaign_runner.run_grid_campaign` (checkpointed `results.jsonl`,
kill-safe resume via re-invocation with the same arguments). `--report` writes
`census_report.json`: per-system evaluated/in-band/domain-invalid/no-converge/error
counts plus the actual in-band `(p, q)` cells — the Stage A deliverable `#858`/`#859`
both specify.

17 evidence tests in `tests/search/test_resonant_atlas_stage_a.py` (not slow, ~6s total):
coprime enumeration correctness/determinism, in-band boundary behavior (sourced against
`#781`'s own real eigenvalues: 105.05 in-band, 14600 out), cell-grid structure and JSON
round-tripping, both deterministic worker fast-paths, one real cheap end-to-end cell
(Uranus-Oberon 3:2, `n_c_steps=1`), and a checkpoint-resume test through the REAL worker
(not a synthetic stub) over a tiny grid.

## 3. Smoke test — measured results

Ran small slices across all four systems in the foreground (no backgrounding), well
under this session's time budget. Total smoke-test compute: well under 10 CPU-minutes.

**A DETERMINISTIC construction-domain finding, not a bug**: `two_body_resonant_seed`
fixes the seed's periapsis radius at `r=1` (the secondary's own orbital radius); the
vis-viva speed there is only real for `a=(q/p)**(2/3) >= 0.5`. At `p,q<=8` this excludes
7 of the 43 coprime pairs per system — `(3,1), (4,1), (5,1), (6,1), (7,1), (7,2), (8,1)`
— which now fail instantly (~0 s) as a clean `"miss"` (`reason=seed_domain_invalid`)
rather than raising. Originally these DID surface as generic `"error"` cells in the first
smoke-test pass; fixed by catching this specific, deterministic `ValueError` and
recording it as an expected negative (see the module's own docstring for the derivation).

**Positive control (Uranus-Oberon), the known-answer check**: seeded and continued 12
distinct `(p, q)` ratios, including the exact 5 resonances Anderson & Kumar 2024's own
survey studies as unstable saddle families (4:5, 5:6 exterior; 4:3, 5:4, 6:5 interior).
**Every one of them, under this harness's naive two-body seed, converged to a
NEAR-UNIT-CIRCLE member instead** (`|lambda|` in `[0.9999999999999887, 1.2733975...]`,
none in-band) — the SAME topology/branch-identity failure mode `#776`'s own "two-body
seed lineage" finding already documented for Neptune-Triton (naive seeding frequently
lands on the WRONG branch, not the paper's own labeled unstable family), now independently
reproduced on a 4th system. This is an honest, important smoke-test result, not a harness
defect: the corrector/continuation/classification WIRING is demonstrably correct (see
below), but **the naive two-body seed alone is not, by itself, a reliable way to locate
the specific unstable saddle family a `p:q` label nominally refers to** — exactly the
caveat the module's own docstring states up front, now empirically confirmed for
Uranus-Oberon specifically rather than inferred by analogy.

**Wiring correctness, independently confirmed**: with a larger `d_jacobi` (5e-3 to 1e-2,
vs. the harness's own coarse-grid default of 5e-4), several Uranus-Oberon cells hit a
genuine `FOLD_RADICAND` stop within 1-5 continuation steps (a real family boundary, the
Jacobi radicand going negative exactly as `cr3bp_continuation`'s own contract describes) —
confirming the corrector, half-crossing-index recovery, and continuation loop are all
functioning correctly end to end, not merely "not crashing." No genuinely in-band cell
was found in this narrow smoke slice at any step size tried; this is NOT evidence the
in-band band is empty for these systems (12-20 cells is nowhere close to the full
43-cell-per-system grid, and per the caveat above, the naive seed likely under-samples
genuinely unstable branches specifically).

**All 4 systems ran cleanly** (Ganymede, Titania, Rhea, in addition to Oberon) — no
crashes, no unexpected exceptions, confirming `#764`'s own "system-agnostic" finding
extends to this pilot's remaining 3 systems.

## 4. Refined Stage A cost estimate

Measured per-cell wall-clock at the harness's own coarse-grid default
(`n_c_steps=9`, `d_jacobi=5e-4`, single-threaded): **~2-25 s for a converging cell**
(mean roughly 10-12 s across ~15 timed samples spanning all 4 systems), **~0 s** for the
7/43 per-system domain-invalid cells, and **~0.01 s** for the (rare, observed once at
1:1) seed-did-not-converge cells.

Full grid: 4 systems x 43 `(p, q)` pairs = **172 cells**. At ~10-12 s/cell average
(weighted by the ~36/43 real-compute cells per system) this is **roughly 0.5-0.7
CPU-hours total, single-threaded** — **two to three orders of magnitude below** `#858`'s
own inferred 50-100 CPU-hour budget.

This gap needs to be read carefully, not taken as "Stage A is nearly free":

- `#858`'s own 50-100 CPU-hour figure (Sec. 7, explicitly flagged as uncertain, "could be
  off by 2-3x") appears to price a MUCH DEEPER continuation than this harness's own
  coarse default — its own cited reference point, `#777`'s heaviest single
  family-continuation test, ran to **376 s and up to 482 gauntlet-passing members** on
  ONE branch, i.e. a walk covering a family's FULL existence range, not a fixed ~10-step,
  `d_jacobi=5e-4` local neighborhood. At the harness's own default step size, 9 steps
  covers a Jacobi-constant span of only `9 * 5e-4 = 4.5e-3` — nowhere near a typical
  family's own full extent (e.g. the vendored Neptune-Triton DPO family spans `~0.013` in
  C over 8 members; the 3:2 family spans considerably more before its own fold). This
  harness, AS DEFAULTED, samples a narrow local neighborhood of the two-body seed's own
  natural C, not "the family's own existence range" `#859`'s registration literally asks
  for.
- Confirmed directly this task: raising `d_jacobi` by 10-20x (to 5e-3 - 1e-2) makes the
  walk terminate at a genuine fold within just 1-5 steps instead of running the full
  budget with `stop_reason=max_steps` — i.e. the CURRENT small-step default is not even
  reaching a fold, so it is measuring a local eigenvalue value near the seed, not a
  genuine "spans the existence range" survey.

**Refined recommendation for the coordinating session, not resolved by this smoke test**:
before dispatching the real run, decide (and re-time) ONE of:

1. Keep the harness's own small `d_jacobi` (fast, ~0.5-0.7 CPU-hr total) but treat it
   explicitly as "eigenvalue near the seed's own natural C only", not a full existence-
   range survey — cheapest, but a materially weaker triage signal than `#858`/`#859`'s own
   framing implies.
2. Widen `d_jacobi` (e.g. 5e-3 to 1e-2, or an adaptive step) and set `n_c_steps` high
   enough (e.g. 50-100) to let `FOLD_RADICAND`/`FOLD_REVERSAL` terminate the walk
   naturally at the family's own real boundary — closer to `#777`'s own reference cost
   basis (up to 376 s/branch observed there), which would land Stage A's real total
   somewhere between this smoke test's ~0.7 CPU-hr floor and `#858`'s own 50-100 CPU-hr
   estimate; the true number needs its own short calibration run (a handful of cells at
   the wider step, timed directly) before committing to the full grid.
3. Address the two-body-seed topology-identification gap directly (Sec. 3 finding above)
   before spending real budget on Stage A at all — e.g. vendor a small number of
   published anchor rows per system the way `#776`/`#777`/`#781` did for Neptune-Triton,
   at least for the positive control, so the "does this harness recover known unstable
   families" question has a real answer rather than the honest miss reported here.

This module's own `DEFAULT_D_JACOBI`/`DEFAULT_N_C_STEPS` are left at the harness's
original coarse-grid design (`5e-4`/`9`) — a deliberate choice to keep the smoke test's
own footprint small and the module's defaults conservative, NOT a recommendation that
the coordinating session run the real Stage A sweep unchanged. Both are `--d-jacobi`/
`--n-c-steps` CLI flags on the driver script for exactly this recalibration.

## 5. Status

Harness built, wired through `campaign_runner`, smoke-tested across all 4 systems
(foreground, small slices, well under a session's time budget). 17 fast evidence tests
green. `ruff check`/`ruff format --check`/full `mypy src tests` clean. **The real Stage A
sweep has NOT been run.** Handed to the coordinating session per this task's own explicit
scope boundary, along with the cost-model recalibration question in Sec. 4 (not resolved
here — a real decision point, not a rubber stamp).

Commits: `4b80b09c` (harness + driver + tests).
