# Outstanding questions — cyclerfinder catalogue

## CURRENT STATE

**This section is authoritative.** It is a dashboard of every currently-OPEN or actively-relevant
task, reconciled against the full `TASK ALLOCATIONS` ledger (below) and each task's own bullet
entry as of 2026-07-15. Update it IN PLACE whenever a task's status changes — do not just append a
newer note further down the file and leave this section stale; that exact failure mode (a correct
resolution sitting under an old section header while a stale duplicate elsewhere looks more recent
and misleads a reader/planner) is what happened to `#500` and is why this section exists. See the
`#500`/`#594`/`#595`-`#599` history in the ledger below for the concrete incident. Scope: this audit
covers the `TASK ALLOCATIONS` paragraph's `#512`-`#601` range (~90 tasks); of these, ~80 are
CLOSED/DONE/RESOLVED/SUPERSEDED and are not repeated here — see the ledger paragraph itself
(search `TASK ALLOCATIONS`) or each task's own bullet entry for history.

### Ready to dispatch — no blocker
- None currently. (`#600` — corrected 2026-07-15: this line was never updated after `#600` was
  dispatched and closed the same day; see its own `✓ DONE` bullet entry — clean negative, 806,400
  candidates, near-miss residual 0.0531 km/s just outside the gate. Removed from this list.)

### Open but blocked / parked
- `#538` — QBCP SE<->EM cislunar cycler correction — blocked on a not-yet-built **multiple-shooting
  GMOS torus corrector**: single-period GMOS structurally cannot converge the violently-unstable
  EM-L1/L2 region (root-caused by `#544`). Portfolio-**PARKED** 2026-07-10: near-zero novelty payoff
  even on success (most-published corridor in the field, 4+ sessions already spent, zero catalogue
  rows). Not cancelled, just deprioritized behind breadth-on-unscreened-systems work.
- `#539` — generalize `#538`'s corrector for reuse (Jovian/Uranian re-screens) — blocked on `#538`
  actually landing a working corrector first; parked alongside it. (`#540`/`#541`, the two other
  tasks in this same sub-thread, are each independently superseded/dead — see the ledger.)
- `#542` — learned-seed generative GA/diffusion warm-start for corrector basins — deferred,
  speculative; explicitly waiting on more corrector-runlog diversity before it's worth building.
- `#556` — large-rotation-number quasi-halo torus corrector — parking lot, **not auto-fired**, needs
  a user greenlight. Would let the `#522` linking-number screen attempt the one remaining
  frequency-matched Owen & Baresi L1<->L2 reproduction shot (`#555`'s terminal qualified-negative
  left exactly this one blocker: the L1 quasi-halo can't be built at the needed amplitude).
- *(sub-task, not its own `#NNN`)* `#503`'s "expand `#267` goldens" mining step remains open (its
  acquisition sub-task itself is closed, per `#595`).
- *(sub-task, not its own `#NNN`)* `#596`'s Russell Table 3.4/3.9-3.11 backfill: 161/197 candidate
  rows were successfully backfilled 2026-07-15; the remaining 36 (`AR < 1.0`, a genuine model-
  boundary case) and 3 genuinely-uncatalogued rows are open for a future pass.
- ~~`#557`~~ **REMOVED from this list 2026-07-15 — was ALREADY CLOSED, not open.** This dashboard
  entry was wrong: `#557` was fully run and registered as a clean negative on 2026-07-11 (commit
  `069d8d2`); this list's "awaiting a user scoping decision" text was itself a stale-header artifact
  (see `#557`'s own bullet entry, corrected in this same pass) that this dashboard's own audit
  reproduced instead of catching. See `[[feedback_outstanding_current_state_maintenance]]`.
- ~~`#516`/`#517`~~ **REMOVED from this list 2026-07-15 — both ALREADY RAN, both EMPTY, not open.**
  This dashboard entry was wrong (and was propagated into a same-day discovery-strategy planning
  pass's recommendation before being caught): both ran to completion on 2026-07-01, both closed
  with a registered negative in `data/empty_regions.jsonl`, and both scripts were committed the next
  day (`42ca41e`, whose own commit message says plainly "All three ran to completion; none closed").
  The "never run"/"uncommitted working-tree" framing was true only briefly, before that commit —
  see `#515`/`#516`/`#517`'s own corrected bullet entries. (`#518` remains closed, superseded by
  `#522`, as this list previously noted.)
- `#520` — Comprehensive 3D cross-system closure sweep (8,640-point grid,
  `scripts/run_520_comprehensive_3d_search.py`) — **ABORTED** 2026-07-02 after 12+ hours with zero
  output (a scoping failure: too coarse to find anything AND too expensive to finish), explicitly
  **not a negative result**. No re-attempt has been made; would need a timed pilot + incremental
  checkpointing first per the DELTA-2026-07-02 review's own diagnosis.

### In progress
- None identified as actively running as of 2026-07-15. The most recent worked chain (`#594`-`#598`,
  same day) is fully closed out.

### Ambiguous, not resolved here
- None found unresolved in the `#512`-`#601` range as of this audit (2026-07-15). The known
  "stale-duplicate" pattern instances in this range (`#521`'s "Phase 2 not built" line, `#545`'s
  "ready to dispatch" mischaracterization, and the `#500` duplicate itself) were already caught and
  corrected in place by `#594`/`#599` before this audit ran — verified directly against each task's
  own bullet entry, not just the ledger paragraph's summary. If a similar contradiction resurfaces
  (one entry implying a task is done, another implying it's still open), tie-break via
  `docs/notes/` for a dedicated verdict note before trusting either copy, per the process this
  section itself was created to enforce.
- **Two more instances of the same pattern, found and corrected during this pass:** (1) `#589`'s
  `TASK ALLOCATIONS` one-liner still reads "flagged for user review before dispatch, not auto-fired,"
  but `#589`'s own bullet entry ends with a 2026-07-14 Fable design-review **NO-GO** and "**Not
  dispatched; do not dispatch as scoped**" — the proposal is killed, not pending review (its
  salvageable reframe shipped separately as `#591`, already closed). Removed from the blocked-list
  above for this reason — do not re-add it as "awaiting review." (2) `#582`'s `TASK ALLOCATIONS`
  one-liner still reads "full novelty sweep not yet dispatched," but `#582`'s own bullet entry records
  the 5-MMR novelty sweep as completed 2026-07-14 (commit `b0225f4`, clean 0/104 negative, stamped to
  `empty_regions.jsonl`) — `#582` is CLOSED, not open, and is correctly omitted from the dispatch-ready
  list above. Neither the master ledger paragraph's inline text nor a prior audit (`#594`/`#599`) caught
  either of these; both were only found by reading the task's own dedicated bullet entry directly. The
  master ledger paragraph itself has NOT been hand-edited to fix these two lines (out of scope for a
  pure-insertion dashboard edit) — a future editing pass through that paragraph should correct both.
  (3) Same pattern again on `#583`: its `TASK ALLOCATIONS` one-liner says "full novelty sweep not yet
  dispatched, needs radial/per-family partitioning first," but that partitioning redesign, the full
  16-partition x 3-seed sweep, and its downstream adjudication all actually completed via the
  `#583`->`#586`->`#588`->`#590`->`#591` chain (all CLOSED) — `#583` is correctly omitted above too.
  A fourth ledger line (`#585`'s "the 0/104 result is very likely a narrow-search-box artifact, not a
  physical negative") is likewise stale: `#585` itself ran the widened-bounds sweep to completion
  2026-07-14 (commit `83607a1`, 0/78 asymmetric, a materially stronger negative) and is CLOSED, also
  correctly omitted. The master ledger paragraph is accumulating this staleness faster than anyone is
  fixing it — treat its inline one-liners as a rough index only, never as a status source, until
  someone does a dedicated pass to hand-correct all of `#582`/`#583`/`#585`/`#589`'s (and possibly
  more) one-liners in place.

### Pre-#512 items (NOT covered by the #512-#601 audit above — spot-checked only, not exhaustive)
This dashboard's main audit only reconciled the `TASK ALLOCATIONS` paragraph's `#512`-`#601` range.
Tasks below `#512` have their OWN, older ledger paragraph(s) further down this file and have NOT been
comprehensively re-audited here. The three below are known, spot-checked opens from this same
session — do not treat their absence from the sections above as "closed," and do not treat this
short list as a complete pre-#512 audit either.
- `#315` — Circumbinary/binary-star μ-gap sweep. Its own bullet entry (search `**#315**`) still says
  "OPEN" — **that line is STALE**, same pattern as `#500`: `#494`'s entry states "**Closes #315/#252/
  #255** positively" (search `Closes #315`). Treat `#315` as CLOSED; the "OPEN" bullet itself needs
  correcting next time someone is in that part of the file.
- `#316` — Cross-system cycler framework (Sun-Earth <-> Earth-Moon manifolds) — genuinely ambiguous,
  not stale: possibly redundant with the `#405`/`#411` SE<->EM heteroclinic-cycle work, but no
  explicit link exists either way. Needs a human merge/supersede-or-keep-distinct decision.
- `#317` — PINN-based pre-filter for sweep-impossible regions — genuinely OPEN, no resolution found
  anywhere in the file as of last check.

Long-form log of research questions, source-access gaps, parameter
contradictions, and out-of-paradigm flags encountered while compiling
`catalogue.yaml`. The YAML's per-entry `notes:` field carries
short-form caveats; this file carries the discussion threads that
don't fit there.

**Resolution policy:** when a question is resolved, prefix its heading
with `✓ Resolved (YYYY-MM-DD)` and add a one-line pointer to the
resolution (commit SHA, spec section, errata-investigation note, etc.).
Do not delete the original question text — the audit trail matters.

---

# Project state at a glance (updated 2026-07-09)

## DELTA 2026-07-02 (independent Fable-model review — sustainability audit + novel-orbit discovery proposal) — read this first

Dispatched a fresh, isolated-context review agent (read-only; edited nothing) to audit the codebase's
sustainability/correctness and propose a genuinely new discovery method, seeded twice with live evidence from
this session's #513-numbering-collision and #520-abort incidents (see DELTA below). Full transcript not saved;
findings + task allocations recorded here.

**Audit findings (most-severe first):**
- **F1 HIGH — the anti-redundancy safeguard is split-brain, and the half that works isn't used.** There are TWO
  negative registries. `data/empty_regions.jsonl` (50 reports) has a real mechanical gate (`should_sweep()` in
  `src/cyclerfinder/data/method_capability.py:125`) wired into `DiscoveryCampaign`. `data/negative_results.yaml`
  (13 entries, the one everyone actually writes to) has **zero programmatic readers anywhere in `src/` or
  `tests/`** — grep-verified. It is write-only prose; five scripts append to it after finishing, nothing reads it
  before starting. 192 one-off scripts (44 `run_NNN_*.py`, no shared harness) each bypass the working gate
  entirely. This is exactly how #515-517 re-ran the approach #411's own entry says NOT to ("NOT more
  revolutions / finer grids") and how a task number collided twice in 24h with nothing catching it mechanically.
- **F2 HIGH** — independently confirms the #520 12-hour-abort root cause below: a scoping failure (no timing
  pilot before an 8,640-point grid), compounded by the same no-incremental-checkpoint / silent-exception-swallow
  anti-pattern.
- **F3 MEDIUM** — CI never runs `-m slow` (no scheduled job exists anywhere in `.github/workflows/`); V-tier
  evidence tests that happen to carry the slow mark could silently rot with nothing mechanically catching it.
- **F4 MEDIUM — validation-census correction.** Computed directly from `data/catalogue.yaml` (not cached):
  **V0:79 / V1:26 / V2:8 / V3:2 / V4:1 / V5:0 / null:240** (of 356). This does NOT reconcile with prior
  DELTA-cited tallies in this file (e.g. "V0:319/V1:27/V2:7/V3:2/V4:1" below) — prior figures likely conflated
  explicit `V0` tags with UNTAGGED (`null`) rows. 240 rows carry no `validation_level` claim at all and nothing
  ratchets that count, so tier coverage can silently erode undetected. `scripts/backfill_validation_level.py`
  already exists to close this; a null-count ratchet does not yet exist.
- **F5 MEDIUM** — the 192-script sprawl is the structural bypass vector for F1; a shared `scripts/_harness.py`
  (preflight + runlog + checkpointing + registry append) would close it.
- **F6 LOW** — test health is genuinely good: independent full-suite + ratchet reruns both green, matching this
  session's verified **3,061-test** count (correcting the earlier unverified "690 tests passing" claim from the
  #515-520 cleanup pass — see DELTA below).
- Credited strengths (agent's words): the digest discipline, DELTA-structured logs, sourced-only goldens, and the
  #425 negative-registry re-audit are "unusually good research hygiene" — the failures are concentrated
  specifically in ad-hoc concurrent search dispatch, not the core codebase.

**Proposal — ranked, falsifiable, new task numbers (see TASK ALLOCATIONS):**
- **#521 (precondition)** — Pre-flight search gate: `preflight_search()` checked against BOTH negative
  registries + task-number collision + a required timing-pilot/checkpoint declaration for any grid above N
  points, enforced by an AST-based ratchet requiring every `scripts/run_*.py` to call it. Nothing below should be
  dispatched via one-off script until this exists.
- **#522 (P1, recommended first)** — Coherent-model whiskered-torus connection search: reformulate the SE<->EM
  cross-system question as heteroclinic connections between whiskered invariant TORI of ONE time-periodic
  BCR4BP/QBCP model, instead of gluing two autonomous CR3BPs — the phase that has never closed becomes an
  intrinsic torus coordinate BY CONSTRUCTION, dissolving the 1-DOF obstruction (#411/#496/#512/#515-517) rather
  than re-sampling it at finer resolution. Grounded in Kumar/Anderson/de la Llave (SIAM J. Appl. Dyn. Sys.
  24(1):219-258, 2025) and Owen & Baresi (Astrodynamics, 2024). Re-scopes/subsumes #518. Positive control: the
  #405 forward EM-L2->SE-L2 leg must persist as a torus connection under mu_sun-continuation, plus one published
  BCR4BP/QBCP transfer reproduction.
- **#523 (P2)** — Co-orbital-exchange cyclers: search horseshoe/tadpole/quasi-satellite topologies (Saturn
  Janus-Epimetheus swap pair; Earth quasi-satellites; Mars Trojans) where repeated encounters come from slow
  co-orbital exchange, NOT hyperbolic flyby bend — sidesteps the mass-deficit no-go theorem (#489/#308) that
  killed every prior small-body lane. Zero existing search code for this topology; ~1-2 weeks. Settle the "does
  a slow co-orbital pass count as an encounter" criterion in writing BEFORE sweeping (a #339-style criterion
  trap otherwise).
- **#524 (precondition alongside #521)** — Continuation + deflated Newton as the mandatory search primitive:
  replace fixed-grid + independent-Newton-per-cell with pseudo-arclength continuation (follows solution curves
  instead of sampling them — directly answers why the #520 grid at ~120 uu spacing stepped clean over the #496
  6.2 uu convergent strip) and deflated Newton (enumerates basins instead of hoping a seed lands in the right
  one — answers the recurring family-selection pathology, e.g. #343/VEM wrong-family closures). Gate: must
  re-find >=2 known distinct basins (PC (3,2) + the discarded (3,1) near-primary root from #504) on a control
  problem before any "only N families exist" claim from it is trusted.
  **BUILT 2026-07-02/03 (this session) — see the full writeup in TASK ALLOCATIONS below. Gate note: the
  literal ">=2 known distinct basins incl. (3,1)" criterion above is only PARTIALLY met** — deflated Newton
  independently reproduces the (3,2) anchor from sourced #504 data, and is validated on two closed-form
  multi-root positive controls (cubic, circle/line), but re-finding the SPECIFIC (3,1) near-primary root via
  deflation from a shared (3,2) seed was not attempted: (3,1) was originally found via a different corrector
  path (mu-continuation + C-walk, not a shared-seed scalar residual at the same `(C, half_crossings)` as
  (3,2)), and asserting it as a same-basin deflation partner without checking would itself be the kind of
  unverified claim these gates exist to prevent. Whether (3,1) is actually reachable by deflating (3,2) at a
  shared C is an open question for whoever first wires this primitive into a live #522/#532 search, not
  something this build should force-fit.
- **#525 (sequence after #522-524)** — Learned seed generation from the project's own corrector archive:
  diffusion/generative warm-start (cf. Graebner & Beeson, arXiv 2501.07005) trained on accumulated
  runlogs/checkpoints to propose seeds in unknown basins; revisits the 2026-06-11 Ozaki (arXiv 2111.11858)
  "below breakeven" triage, which the reviewer flags as now stale.
- **#526 (cheap)** — GTOC 13 + Liang operator ingestion: digest GTOC 13 (Oct-Nov 2025 ballistic Jovian
  gravity-assist tour competition — this project's exact problem class run as a world competition; methods
  papers due mid-2026) and encode Liang et al.'s alternating-double-cycler construction (JGCD 2024, DOI
  10.2514/1.G008387) as a reusable genome operator for near-commensurate phase-mismatch walls.
- **Recommendation**: #522 first, with #521 + #524 landed as preconditions; #523 in parallel (cheap, disjoint
  code path). Reviewer's closing line: "the last 48 hours demonstrated, twice and expensively, that this
  project's discipline currently lives in prose, and prose does not gate concurrent agents."

**#521 PHASE 1 DONE (2026-07-02):** the negative-registry migration described in F1 above is complete — see the
`#521` entry in TASK ALLOCATIONS below for the full report (overlap analysis, migration detail, an honest
provenance-gap finding on `branch_C32_C_3.1774`, and a permanent regression test proving `should_sweep()` now
mechanically skips an equal-capability re-run of the #405/#411 SE<->EM search — the exact pattern that produced
#515-517). Ratchet suite (`tests/data tests/search`) verified green: 1728 passed, 0 failed.

**#521 PHASE 2 DONE (2026-07-02):** the actual gate is built — `src/cyclerfinder/data/preflight.py`'s
`preflight_search()`, three checks (task-number hygiene against `data/OUTSTANDING.md`'s TASK ALLOCATIONS +
filename self-consistency; registry subsumption via `should_sweep()`; a timing-pilot requirement above
`LARGE_GRID_THRESHOLD=500` points), with an explicit, logged `override_reason` escape hatch and every invocation
appended to `data/runlogs/preflight_runlog.jsonl` (gitignored). Enforced by an AST ratchet
(`tests/scripts/test_scripts_call_preflight.py`) requiring every `scripts/run_*.py` to call it; the 40
pre-existing scripts are in a frozen `_LEGACY_EXEMPT` list (retrofitting them is out of scope here), but the four
scripts whose incidents motivated this gate (#515-517, #520) were retrofitted as the reference pattern —
verified LIVE against the real on-disk registry: running any of #515/516/517 again now correctly raises
`PreflightBlockedError` (region already covered — the exact re-run this gate exists to stop), and running #520
again correctly raises on the missing timing pilot (8,640 points with none measured — the exact scoping failure
that cost 12+ hours). 14 unit tests (`tests/data/test_preflight.py`) + the AST ratchet all pass; full
`tests/data tests/search tests/scripts` ratchet green. Not yet built: a `scripts/_harness.py` shared skeleton
(F5) so future scripts get the preflight call, checkpointing, and registry append for free instead of
copy-pasting the boilerplate — worth doing before the next batch of #522-529 scripts are written.

**Second independent pass (2026-07-02, same day, fresh Fable-model context, explicitly told not to repeat
#521-526):**
- **#527 (its top pick)** — Point an already-built, zero-usage systematic enumerator at an untouched target.
  `src/cyclerfinder/search/da_hotm_enumeration.py` (#450, DA/HOTM: a CR3BP fixed-point/periodic-orbit enumerator
  over a Jacobi-constant band, μ-agnostic, fully tested) has **never been called by any campaign script** —
  grep-verified. Propose pointing it at the Sun-Jupiter Hilda/quasi-Hilda 3:2-exterior-MMR band — a region in
  neither negative registry nor the catalogue. Sidesteps the #489 flyby mass-deficit no-go a different way than
  #523: here Jupiter itself (the massive primary) does the gravitational structuring at each pass, not a small
  body acting as a flyby node. Sourced positive control (independently re-fetched and confirmed real): Guido &
  Efthymiopoulos, arXiv:2604.00679 (2026-04-01), computes manifold/heteroclinic structure for exactly this MMR
  band — the enumerator must recover manifold-consistent unstable POs before any "no Hilda cycler" negative is
  trusted. Cheapest of the three: no new build, a config/target change plus a one-paper digest.
- ✓ Resolved (2026-07-06) **#528** — A low-thrust-native validation-gate gap (not a search gap). Resolved in commit cfac8e0: implemented verify_low_thrust_feasibility and a 500 m/s continuous-thrust ceiling, validated on sanchez-net-2022-eem-cycler1. `dv_band: low_thrust_sep` and
  `trajectory_regime: low-thrust` exist as schema fields and `core/sims_flanagan.py` exists as code, but **zero
  catalogue rows use either** (grep-verified) and the low-thrust branch of the acceptance gate
  (`verify/dv_band_acceptance.py`) literally `return`s `None` — no acceptance window at all. The V0/V2 gates
  hardcode a pure-ballistic "|V∞| preserved across each flyby" assumption. Consequence: if the already-queued
  #519 (VEM low-thrust search) finds a genuine candidate, there is currently no gate that can grade it. Land
  before or alongside #519 actually running. Positive control: reproduce Pascarella et al. 2022 (AAS 22-015,
  already in this project's corpus) under the new gate to a stated tolerance.
- **#529 (defer — real, but scope-creep risk)** — Inter-cycler networks / taxi-transfer as a catalogue concept
  the current 5-class taxonomy cannot express. This project's own corpus (Sanchez Net et al. 2022, already
  digested — `docs/notes/s1l1-target-topology-mining.md:156-263`) describes published fleet-of-cyclers networks
  (3/5/6-cycler, shared downlink schedules) that a prior note explicitly says "could not be catalogued" — the
  schema's only inter-row relation is `precursor_mga.inserts_into` (a single pointer), with no concept of a
  cycler *set* or a taxi transfer between two catalogued cyclers. Reframes as combinatorial phasing-compatibility
  search over existing/new catalogue rows (reuses `model/score.py`'s `taxi_cost_kms`), not a new integrator.
  Ranked last deliberately — it changes what a catalogue row/relation can *be*, so it wants a deliberate scoping
  discussion before a sprint slot, not because the idea is weak.
- **Additive-not-new** (fold into #524/#525's scope, not standalone): a possibly-stronger multi-rev Lambert
  solver claim (MDPI *Dynamics* 6(1):3, 2026-01-05 — lower-tier venue, treat with skepticism) that would close a
  silent false-negative vector in the `run_516_multirev*`-style leg enumeration; a Conley-Zehnder topological
  invariant (Aydin, arXiv:2602.16354) to strengthen #524's "must re-find known basins" gate cheaply from STM
  output the codebase already produces; a parity-based seed-initialization scheme (Park & Howell,
  arXiv:2606.08485) as a non-learned alternative/complement to #525's generative model; and a note that #450's
  DA/HOTM enumerator already embodies much of #524's "enumerate, don't grid-guess" philosophy — worth checking
  before building #524's continuation/deflation machinery from scratch.
- **Recommendation**: #527 first (cheapest — no new build), #528 before #519 actually runs, #529 deferred for a
  scoping conversation. All three still sit behind #521.

## DELTA 2026-07-01 (session C — 3D Lift capability landed + 16-core parallel sweep search) — read this first

- **#515 (working-tree)**: Shipped the 3D lift framework for cross-system cycles (uncommitted).
  - Implemented `correct_cross_cycle_3d` over fixed out-of-plane amplitudes ($z_{em}$, $z_{se}$).
  - Resolved the out-of-plane position gap matching obstruction by physically scaling $z_{se}$ relative to $z_{em}$ by the ratio of system length scales ($384,400 / 149,600,000$).
  - Aligned 3D Floquet eigenvector signs with planar counterparts to enforce consistent manifold directions.
  - Refactored search script `scripts/run_515_cross_system_3d_search.py` using `cyclerfinder.parallel.parallel_sweep` (Loky Backend) to execute on all 16 concurrent threads of the amdnuc machine. All unit tests pass.
- **#516 (New)**: Multi-Revolution 3D Patched Search (n_em, n_se > 1) to bypass the single-revolution phase-closure wall.
- **#517 (New)**: Asymmetric/Mixed Libration Pairs in 3D (scanning EM-L1 <-> SE-L2 and EM-L2 <-> SE-L1).
- **#518 (New)**: 3D BCR4BP Continuation using 3D patched orbits as coherent 4-body seeds.
- **#520 (New, working-tree)**: Comprehensive 3D cross-system closure sweep (`scripts/run_520_comprehensive_3d_search.py`,
  8,640-point grid: 4 libration pairs x 6x6 (c_em, c_se) x 6 z-amplitudes x 10 rev-count pairs). Dispatched by a
  separate concurrent agent; was originally mis-numbered `#516` (collision with the Multi-Revolution task above) —
  renumbered here per [[project_task_numbering_convention]]. CAUTION before trusting a 0/8640 negative from this
  grid: its `C_SE_GRID` (6 uniform points over [3.0003, 3.0009], ~120 uu spacing) steps clean over the #496
  convergent strip [3.000854, 3.00086] (6.2 uu wide) -- no n=1 cell in this grid can even reproduce #496's own
  near-closure. Needs a positive control (recover the #496 n=1 near-closure) and sub-uu c_se resolution near the
  strip before an all-negative result is trustworthy; otherwise it is largely re-confirming #512's clean negative
  at high cost with a real chance of stepping over a solution near the strip.
  **UPDATE 2026-07-02: ABORTED, NOT a negative result.** The run pegged all 16 cores for 12+ hours and never
  finished -- no entry was ever written to `data/negative_results.yaml` (checked: none present), so it produced
  ZERO usable output; do not treat its silence as evidence either way. No OOM kill in kernel logs and the machine
  has not rebooted since before the run started, so this was not a crash -- it was killed (or given up on) after
  running out of patience, and everything computed up to that point is lost because the script has no incremental
  checkpointing (writes to `negative_results.yaml` exactly once, only on a fully-completed run) and suppresses
  per-point output for infeasible/failed points ("silently fail ... to prevent log clutter" -- evaluate_point's
  except branch), so there was no visible progress or ETA during the whole run. Root-caused via direct timing
  probe (2026-07-02, this session): even the CHEAPEST cell in the grid (n_em=1, n_se=1, a config #515 already
  computed) took >16s under load and a single n_em=3/n_se=4 infeasible cell took ~20s just to raise-and-discard
  its exception. 8,640 points / 16 workers = 540 points/worker; at a plausible 60-90s/point average (higher
  n_em/n_se pairs up to (4,4) integrate multi-year multi-revolution arcs with an 8-iteration FD-Jacobian
  bounded_ls solve) that lands at 9-13.5 hours wall-clock -- consistent with the observed ">12h, still running"
  report. This is a SCOPING failure (grid sized without a timing pilot), not a solver hang. Before re-attempting
  #520: (1) time a small (~32-64 point) pilot spanning the n_em/n_se and c_em/c_se ranges to get a real per-point
  cost distribution; (2) append each result to negative_results.yaml (or a runlog) incrementally as it completes,
  not only at the end; (3) then size the grid (or add a wall-clock budget + checkpoint-and-resume) to what that
  pilot says is actually affordable. This is IN ADDITION TO the coarse-grid/no-positive-control problem noted
  above -- #520 as originally scoped was both too coarse to find anything AND too expensive to finish.

## DELTA 2026-07-01 (session B cont. — admission + capability fixes + Ross-corpus mined) — read this first

The #286 CAPABILITY FRONTIER IS NOW ESSENTIALLY COMPLETE (3D / BCR4BP / QP-tori / epoch-MGA all
built; #293 ER3BP — the last un-started axis — building now). The program pivots from BUILDING
capabilities to RUNNING discovery on them. Landed this block:
- **#494 ADMITTED** (e7bca1b): **5 new V1 rows → catalogue 356** (V0:319/V1:27/V2:7/V3:2/V4:1).
  4 binary (k₁,k₂)-cycler μ-family reps (μ=0.001 SJ / 0.1 / 0.3 / 0.5 equal-mass, P1/P2 generic
  primaries) + the **first-ever Pluto-Charon (3,2) cycler** (`ross-rt-pc-cycler-32-2026`, μ=0.10876,
  C=3.5795, ν≈0; fresh real-system instantiation, lit not-found). **Closes #315/#252/#255** positively.
  The 2 EM reps were SKIPPED (already V2 as ross-rt-em-cycler-*). Verified green (tests/data+search).
- **#496** feasibility-first corrector shipped (627228d/1edfe3c): scan-resolution wall BROKEN (both
  legs DO converge at c_em=3.150, c_se=3.00086 with return_scan_n=8; prior 217 km gap was scan
  coarseness, not physics). Phase-closure wall CONFIRMED for n_em=1/n_se=1: convergent strip is only
  6.2 μu wide [3.000854, 3.00086]; best |R|=0.52 rad; Newton step requires c_em≈3.158 (outside EM-L2
  family) AND c_se≈3.000844 (dead zone). **#411 NOT CLOSED.** Next: (n_em,n_se) sweep; 3D z-slicing.
- **#497 + #497b** (f56719b / 349229b / 0d6fba1): cap recalibration DISPROVEN; scorer VALIDATED on
  Braik data (Table 4 exact; dc_refined → C32 dominant); proxy was 30-60× off-scale. **be64207 FIXED
  the pedestal**: heading-mismatch patch was a constant ~89 m/s term (dv_turn used unit speed, not
  per-family min-turn cost). After fix: median ratio ~1.9×, Spearman ρ=0.84 vs Braik. C32 xfail
  re-diagnosed as NODE-SET effect: R52-U (excluded from our 12-node set) carries C32's strongest edge
  (0.62 m/s); removing R52-U from Braik's own matrix reproduces our ranking. Flipping C32 gate needs
  R52-U recovery. Verdict: docs/notes/2026-07-01-497-proxy-rebuild-verdict.md.
- **#498/#499/#503 DONE** (4aef28f / a4cc3c1 / 1b332b0): **16 Ross-group papers fetched from
  ross.aoe.vt.edu (no denials), filed, mined, indexed.** Key reuse unlocked: Ross-Scheeres 2007
  Keplerian map + Grover-Ross 2009 → #500; Gómez 2004 z-slicing (4D→1-param 2D intersections) →
  #496 3D lift + #291/#306; **Fitzgerald 2022 ER3BP L1 IC+monodromy → #293 (unblocked the un-started
  axis)**; Onozaki 2017 full BCR4BP constants → #292; Koon 2002 Petit Grand Tour (1208 m/s, 57%
  Hohmann) → #318; Naik-Lekien 2017 + Lober → transport scoring; KLMR2001 tube↔resonance → #267.
- **#500 DONE** (f9adc1a/86c1bc5): Keplerian map built + **17/17 sourced positive controls PASS**
  (RS07 1:2 fixed point a_res=2^(2/3)=1.587; chaotic migration; GR09 controlled ΔV brackets 160 m/s).
  Genome verdict = **clean negative for catalogue**: it's a Jovian moon-to-moon sub-leg planner, NOT
  an interplanetary cycler genome — banked as a tool for a future #318 Jovian sub-leg pipeline.
- **#293 DONE** (b9a7d27/e3a8896): ER3BP corrector + e-continuation, positive control PASS (Fitzgerald
  2022 ER3BP L1 via the paper's own non-pulsating frame transform, x to 0.078%; monodromy structure
  matched), e-continuation 21/21 smooth. **#286 axis-closure verdict: NO new species** (e>0 refines
  CR3BP families, no novel cycler topology — confirms the prediction). This **FORMALLY COMPLETES the
  #286 capability frontier** (all 5 Track-A axes built).

**TIER-1 DISCOVERY RESULTS (both clean negatives — the frontier is mined):**
- **#501** broadened real-eph Galilean joint search (e5b2339/d0c8def): 6 sequences × 512 Sobol = 3072
  cells / 213 feasible / 26 shot / **0 closed**; positive control (Liang Member D) PASSED → trustworthy.
  Clean empty-region map registered. The strongest novel lane yields no novel Galilean cycler at breadth.
- **#504** deeper Pluto-Charon (k₁,k₂) sweep (2b2f74a/d289fcd): only (3,2) yields a stable member (the
  already-admitted ross-rt-pc-cycler-32-2026); (1,1)/(2,1)/(2,2)/(3,1)/(3,3) all clean negatives (the
  (3,1) wrong-topology near-primary correctly discarded). Pluto-Charon admits exactly ONE cycler family.
  NO new admissions.
- **#505 DONE** (4905f11/3e15444): PC (3,2) promoted V1→V2-ballistic via long-span IAS15 evidence
  (10 synodic periods, pos error <3 km drift at T=600 d). SPICE lane not needed for V2-ballistic.
- **#506 DONE** (5004822): PC (3,2) V3-scope assessment — **stays V2-ballistic**; JPL Horizons not
  available for Charon, requires NAIF SPK kernel (SAT441l) not present on this host. V3 upgrade gated
  on acquiring the kernel; V2-ballistic is the correct current tier.
  **Superseded (2026-07-01, same session): #510 fetched the kernel (`plu060.bsp`) and #511 built +
  ran the actual differential-correction V3 attempt — CONFIRMED STAYS V2-ballistic for a structural
  reason (period incommensurate with Charon's own orbit), not a missing-kernel gap. See #511's
  verdict, `docs/notes/2026-07-01-511-pluto-charon-realeph-verdict.md`. #550 (allocated 2026-07-10)
  re-proposed this as if unresolved — corrected there.
  **⛔ CLOSED — DO NOT REOPEN (bolded 2026-07-15, full-file audit): this has already been
  re-litigated once by mistake (#550) despite the resolution above being clear. PC(3,2) stays
  V2-ballistic for a STRUCTURAL reason (period-incommensurability with Charon), independently
  confirmed twice (#511's differential-correction attempt, then #550's own re-check) — this is
  not a missing-kernel or missing-data gap that a future kernel fetch would fix.**
Net Tier-1: discovery on the built substrate confirms #492 exhaustion — no novel cyclers; the deliverable
is the honest empty-region maps + the confirmation that the admitted rows are the frontier.
- Process note: subagents hit a monthly spend limit mid-session (the #494 admission agent died AFTER
  committing — verified via git, not its message; [[feedback_long_agents_commit_incrementally]]).
  Also a Kumar 2509.12675 DUPLICATE was caught + reverted (cfcd7e5; [[feedback_corpus_check_index_not_filenames]]
  — checked md5 but not the arXiv id; fixed).

**FORWARD PLAN (the post-capability-frontier phase):**
- **Tier 1 (the deliverable — RUN discovery on the built substrate):** #501 broadened real-eph #318
  joint search (not CGCEC-densify); circumbinary Pluto-Charon deeper + other-binary (k₁,k₂) sweep;
  #500 Keplerian-map tour search if its genome is cycler-viable; #519 VEM Multi-Synodic Low-Thrust Signal Search. HONEST: novel hits are rare (#492;
  only #312 Uranus across the whole arc) — expect mostly V0-known + occasional fresh instantiations.
- **Tier 2 (deep fixes):** #496 (n_em,n_se) sweep → 3D z-slicing (Gómez 2004) for cross-system closure;
  #497 R52-U recovery to flip C32 gate (node-set, not proxy).
- **Tier 3 (validation + consolidation):** PC (3,2) V3 upgrade RESOLVED (#510/#511, 2026-07-01) —
  confirmed stays V2-ballistic (structural, not kernel-gated);
  #487 V4_qp gauntlet (de-prioritised).
- **Tier 5 (3D Dynamics & Multi-Rev Search):** #515 3D lift framework for cross-system cycles (uncommitted,
  status unverified — check before reviving); #516 Multi-Revolution 3D Patched Search (n_em, n_se > 1,
  uncommitted); #517 Asymmetric/Mixed Libration Pairs in 3D (uncommitted); #518 3D BCR4BP Continuation
  (**superseded by #522 just below, per #522's own text — not independent forward work**); #520
  Comprehensive 3D sweep (needs a positive control before its negative is trustworthy — see DELTA above).
- **Tier 6 (Novel-Orbit Discovery Reboot, per the 2026-07-02 independent review — see DELTA above):** #521
  pre-flight search gate (precondition — DONE, both phases; new scripts now blocked without it); #522
  coherent-model whiskered-torus connection search (recommended first; re-scopes #518); #523 co-orbital-exchange
  cyclers (RAN, strong-but-incomplete evidence — see below); #524 continuation + deflated-Newton search primitive
  (precondition alongside #521); #525 learned seed generation from the corrector archive (sequence after
  #522-524); #526 GTOC 13 + Liang operator ingestion (cheap); #527 point the idle #450 DA/HOTM enumerator at the
  Sun-Jupiter Hilda band (RAN, clean understood negative — see DELTA above); #528 low-thrust
  validation-gate build (land before/alongside #519 running); #529 inter-cycler network taxonomy extension
  (deferred — needs a scoping discussion); #530 unstable-manifold propagation for co-orbital/resonance cyclers
  (RAN, complete local-manifold negative — see below); #531 heteroclinic/homoclinic connection test reusing the
  existing #314 framework (RAN — a genuine certified/cross-checked homoclinic connection found, still no
  Hill-sphere encounter; see below); #532 multi-orbit resonance-hopping search (the actual Guido-Efthymiopoulos
  hypothesis — cross-family connections, not one orbit's self-connection; medium/uncertain feasibility, needs a
  scoping pilot first).

**#530 (new, allocated 2026-07-02 — motivated by this session's #523/#527 results):** Unstable-manifold
propagation for co-orbital/resonance cyclers. Both #523 (Earth co-orbital horseshoe) and #527 (Sun-Jupiter Hilda
3:2 MMR) independently found the SAME pattern: every CERTIFIED, truly periodic orbit in these families
structurally avoids close (sub-Hill-sphere) approaches to the companion body — that is precisely what makes a
periodic orbit dynamically stable. The genuine close approaches only show up in (a) the TRANSIENT/non-periodic
drift phase (2006 RH120's real minimoon episode, reproduced by direct integration in #523), or (b) per Guido &
Efthymiopoulos (arXiv:2604.00679, #527's positive-control paper), chaotic orbits shadowing HETEROCLINIC
connections between the UNSTABLE manifolds of periodic orbits near a resonance separatrix — a genuinely
different dynamical object that periodic-orbit enumeration (#450 DA/HOTM, used for both #523 and #527) cannot
surface by construction, no matter how much the grid is refined.

Proposal: (1) compute Floquet multipliers (`search/bifurcation_detector.py`'s existing `monodromy` /
`floquet_multipliers`, already used for the #347 Phase 1 work) for the orbits already certified by #523 and
#527, to identify which family members are genuinely UNSTABLE (|λ|>1 for at least one multiplier) rather than
assuming it; (2) for the unstable members, propagate their stable/unstable manifolds (eigenvector-seeded, small
perturbation off the periodic orbit, standard invariant-manifold technique) and test whether trajectories along
them achieve genuine Hill-sphere encounters with the companion body — the actual mechanism both source papers
describe. This reuses existing capability (Floquet/monodromy machinery, the two already-built target systems)
rather than starting from scratch. Distinct from #524 (a search PRIMITIVE for finding more periodic orbits via
continuation) — #530 is about a different dynamical OBJECT (manifolds of already-found orbits), not a faster way
to find more periodic orbits; the two are complementary, not overlapping. Positive control: recover a manifold
consistent with Guido & Efthymiopoulos's reported heteroclinic channel structure (qualitative geometry, not
exact numeric match — the paper does not tabulate precise ICs, per this session's independent check) before any
"no manifold-mediated encounter" negative is trusted. Not yet built.

**TASK ALLOCATIONS (next-unused per [[project_task_numbering_convention]]; #512-#514 committed; #515-#518 for session C working-tree; #519 for low-thrust proposal; #520 for the comprehensive sweep; #521-#526 for the 2026-07-02 review's gate + novel-orbit proposals; #527-#529 for the same-day second-pass review; #530 for the #523/#527-motivated unstable-manifold follow-up; #531 for the #314-reuse heteroclinic-connection follow-up; #532 for the multi-orbit resonance-hopping follow-up; #533 for the genuine QBCP model build; #534 for the #522-split single-system torus connection search; #535 for the transient-drift-phase quasi_cycler search; #536 for the Fable-review-motivated Jovian-moon-tori heteroclinic screening follow-up; #537 for the QBCP cross-system connection search; #538 for the QBCP cross-system periodic orbit correction; #539 for generalizing the #538 corrector + a broadened Jovian-moon re-screen; #540 for a hardened-pipeline Uranian-system re-screen; #541 for a first-pass Saturnian resonant-moon-pair screen; #542 for the #525 learned-seed generative warm-start; #543 for #529's inter-cycler-network scoping discussion; #544 for fixing the EM-L2 QBCP torus mu_sun-continuation convergence blocker found while running #538; #545 for the decoupled Jupiter-Europa/Ganymede CR3BP-level band screen (un-gated from #538/#544); #546 for the decoupled Uranian-system CR3BP-level band screen (un-gated from #538/#544); #547 for resolving #534's transit-vs-non-transit manifold classification blocker (establishing the first genuine positive control for the qp_tori/qp_torus_heteroclinic linking-number method); #548 for a reframed, time-boxed Owen & Baresi positive-control gate on the qp_tori/qp_torus_heteroclinic linking-number pipeline; #549 for a real-binary (k1,k2) genome sweep (Patroclus-Menoetius, Didymos-Dimorphos, Orcus-Vanth, Eris-Dysnomia); #550 for the PC(3,2) V2->V3 NAIF SAT441l kernel fetch; #551 for a GTOC 13 methods-paper corpus mining pass; #552 for a 3D/inclined-releg moontour genome extension (flagged for user review before dispatch, not auto-fired); #553 for a Fable ratification pass on #548's linking-number-pipeline shelve verdict; #554 for the #552-scoping-motivated Neptune/Amalthea empty-region retrograde-correction stamp; #555 for the #553-authorized final, correctly-frequency-matched Owen & Baresi shot at the genuine C=3.15 (final, no-appeal kill criterion); #556 for a large-rotation-number quasi-halo torus corrector (the #555-localized blocker, parking lot -- not auto-fired); #557 for a #535-machinery Sun-Jupiter quasi-Hilda transient-capture quasi_cycler screen (planning phase first); #558 for a #312 family census via the original discovery genome (all Uranian moon pairs, densified rel-offset/phase/tof/n_rev sweep); #559 for the never-dispatched #338 Phase 2 DOY sensitivity scan; #560 for V4-strict Lambert-branch-continuity + degenerate-arc-exclusion robustness fixes (parking lot, not auto-fired); #561 for the #558-mandated Opus/Fable family-vs-coincidence adjudication (>1000 basins need a selection criterion, and #312's own novelty framing is now in question); #562 for the Fable-corrected per-basin continuous-tof commensurability refinement (fixes the grid-aliasing flaw that excluded the two largest/highest-literature-risk directions); #563 for a direct symmetric-closure enumeration replacing grid-search-plus-refinement (the #562 20-count is a lower bound -- narrow basins can fall entirely between #558's original grid points and never surface at all); #564 for the Opus+Fable adjudication of the full #563 30-closure survivor list (Stage-3/writeback gate); #565 for the Fable-corrected re-adjudication of #564 (fixes a semantic inversion of the max_bend_deg triage axis); #566 for the Stage-3 V1-V4-strict gauntlet run on the #565-recommended 5-representative list; #567 for the #566-adjudication-mandated trustworthy epoch-robustness gate (apply #560 fixes + audit-field fix + per-representative epoch scan) as the catalogue-writeback prerequisite; #568 for the #567-diagnosis-corrected step-(4) writeback-readiness characterization (duty-cycle + synodic boundary period framing, not raw pass%/flip%); #569 for the #568-cleared actual catalogue writeback of the #312 Uranian symmetric-closure family (5 new rows + #312 row update + frozen-gate provenance; needs a Fable pass before execution); #570 for the #543-scoped cycler-network schema infrastructure (general relation only, no populated network rows); #571 for a Titan-centered Saturnian moon-pair sweep using the #558-corrected discovery genome (plan under Fable review before dispatch); #572 for the #571-motivated Titan-Iapetus 3D-closure probe (the load-bearing cheap gate deciding whether any #552 build is worth scoping); #573 for the #572-adjudication-mandated widened Titan-Iapetus 3D-closure population probe (measure the real gate-passing closure rate across all Titan-anchored candidates before committing to any #552 corrector build; plan under Fable review before dispatch); #574 for the #573-cleared narrow Titan-Iapetus 3D corrector build (staged: a cheap eccentric-Keplerian closure kill-gate first, then productization + real-ephemeris SPICE SAT441 validation only if it passes — NOT the general n-body/arbitrary-inclination build that was killed; needs a Fable second-opinion on the spec before dispatch; Stage B gauntlet result: 0/15 PASS, root cause a #571 search-methodology gap — free single-closure search with no periodicity constraint, not a genuine symmetric/repeating-cycle search — not a bug in the new validation code); #575 for a #563-method Titan-Iapetus direct symmetric-closure re-search (plan under Fable review before dispatch; FINAL: 9 genuine coplanar symmetric closures found, 0 survive inclination-extension as repeating cycles — clean method-conditional negative, thread closed, `empty_regions.jsonl` stamped, #552 general-capability revival NOT warranted); #576 for a #563-method Galilean-moon direct symmetric-closure search (Jupiter — genericized #563 machinery now proven twice, never applied to this system; plan under Fable review before dispatch; RESULT: 36 gate-passing symmetric closures across all 6 pairs, repeat-instrumentation confirmed genuine, Russell-Strange 2009 comparison found no reproduction among the 2 architecturally-comparable published members); #577 for the Opus+Fable adjudication of the #576 36-closure Galilean survivor list (literature-clearance + pipeline-readiness gate; RESULT: 0/36 novelty-clear, FULL STOP on the pipeline — every pair is published R-S 2009 double-cycler / Liang-Hernandez triple-cycler / pump-tour territory, and the lone Io-Callisto structural clear is a corpus gap not novelty; the 8-mismatched-R-S-members "different topology" framing is incorrect — same ideal model + method, known-class-members); #578 for closing the #577-diagnosed R-S-2009 literature-corpus gap + stamping the Galilean symmetric-closure territory (corpus-accuracy + registry, NOT a pipeline task; plan under Fable review before dispatch; RESOLVED 2026-07-12: 4 new per-pair KNOWN_CORPUS anchors, DOI-coverage ratchet, empty_regions.jsonl stamps, commit 575bd44); #579 for the Fable-diagnosed `literature_check.py` Antoniadou-Voyatzis/Libert anchor mislabeling fix + #287/#301 re-audit (corpus-accuracy, do first per Fable's priority order); #580 for a Richardson-1980 analytic third-order halo-seed generator (fills the halo-branch seed-generator gap Fable independently confirmed); #581 for a Gurfil-Kasdin-2002-style niching-GA search layer, gated on a positive-control reproduction (flagged for user review before dispatch, not auto-fired); #582 for #581 stage 3a, an asymmetric/spatial-isolated 3D CR3BP resonant-family niching-GA search (mu=0.001, Fable-corrected from an initial wrong Earth-Moon scoping; build + positive control DONE 2026-07-13, commit `ba60092`, full novelty sweep not yet dispatched); #583 for #581 stage 3b, a Sun-Earth ER3BP bounded-drift quasi-cycler search widening the validated niching layer beyond Gurfil-Kasdin's own 12 sets (build + positive control DONE 2026-07-13, commit `aafa244`/merge `378b1f1`; a real deep-Hill-sphere basin bug + a theta0-wraparound bug found and fixed via 2 independent Fable reviews; full novelty sweep not yet dispatched, needs radial/per-family partitioning first); #584 for 2 unexplained pre-existing test failures discovered on the new Mac (M3/Accelerate) environment during #582 verification (CLOSED 2026-07-14: CI run 29323181595 confirms 3310 passed/0 failed on the identical tree, mechanism checked and plausible for both — a tolerance-edge V∞ match and an integer winding-number topology flip, both exactly the kind of borderline check a different BLAS backend's rounding can tip; local-Mac-only, no code change); #585 for a #582 follow-up widening the symmetry-breaking z0/xdot0/zdot0 bounds and re-running the 5-MMR sweep (the 0/104 result is very likely a narrow-search-box artifact, not a physical negative); #586 for a #583 follow-up addressing the fitness-landscape/niching discrimination limitation the partitioning redesign exposed (3 candidate approaches scoped, needs a Fable/Opus design-tradeoff read before building); #587 for populating real CR3BP identity fields (Jacobi/period/stability) on the 6 #569 Uranian quasi-cycler rows, found via a cyclers.space website visualization bug (CLOSED 2026-07-14, scope corrected: cr3bp schema field doesn't apply to cycler_class=multi-arc rows — populated invariants/legs tof_days instead, fixed the website placeholder's false tabulation claim); #588 for dedup + live-literature-check + Opus/Fable adjudication of #586's 264-candidate unmatched-bounded pool from the full 16-partition sweep (CLOSED 2026-07-14: dedup 264->45->20 unmatched; live literature-check found 2 Gurfil-Kasdin companion papers, neither closes the gap; Opus+Fable adjudication found 14/20 confidently not novel, 6/20 warrant the #590 bounded follow-up; zero catalogue rows changed); #589 for reproducing/extending the Gurfil-Kasdin 2002 out-of-ecliptic z²-maximization approach as a new search regime distinct from the bounded-family census (flagged for user review before dispatch, not auto-fired); #590 for a bounded follow-up settling whether #588 clusters 40/42/43 (+24/25/39) are genuine Family-J-curve members or single-hit GA artifacts (CLOSED 2026-07-14: 3 new verified QSO CorpusAnchors added, connectivity heuristic ambiguous, long-horizon check found even the published Family J anchor itself escapes at ~30yr — reframes rather than resolves the question, left for a future differential-correction continuation); #591 for Fable's salvageable #589 reframe — ranking already-validated bounded 3D structures (5 known G-K 3D families + #588's cluster 43) by out-of-ecliptic excursion (max |z|), analysis-only, no new search (CLOSED 2026-07-14: cluster 43 beats published Family J by ~16% on excursion while being far more durably bounded per #590 — not novel, but a genuine citable figure-of-merit result); #592 for a QBCP alpha_6 EOM scaling bugfix recovered from an abandoned git stash during a repo-hygiene cleanup pass, verified + applied (CLOSED 2026-07-14); #593 for scoping whether #592's fix changes any past QBCP-based search conclusion (CLOSED 2026-07-14: empirical before/after re-run of the SE-L2 positive control found a genuine, non-trivial impact — same seed converges to a different mean-state solution under fixed vs buggy code — so #533/#537/#544's specific quoted numbers are stale, recommendation is re-derive fresh rather than blind full re-run); #594 for a full-file comprehensive audit of every numbered task (~180 entries, 6 parallel verification passes + a Fable adversarial verification pass same day), correcting #545 (was DEAD+SUPERSEDED, not "ready" as first read), #521/#536/#537/#405/#318/#501/#248/#506/#518/#144/#167 status inconsistencies, the stale #307/#310/#320/#321/#322 duplicate queue, and one self-correction (#344's "genuine gap" claim was itself wrong per Fable's check against literature_check.py's own #349 code comments — resolved not-novel, not open) (CLOSED 2026-07-15, 15 total edits across 2 passes); #595 for discovering #498/#499/#503's "(all MISSING)"/"(MISSING)" framing was stale — triggered when the user re-uploaded 5 PDFs already byte-identical to files acquired+digested 2026-06-30 (`6dfedab`, 16-paper ross.aoe.vt.edu batch); all 12 #498/#499 papers and all 4 #503 papers were already acquired, digested, and CORPUS_INDEX-registered that day but the tracker was never updated (CLOSED 2026-07-15, #498/#499 marked CLOSED, #503's acquisition sub-task closed though its "expand #267 goldens" mining step remains open, #500's #498-gate cleared); #596 for a `data/MISSING_DATA.md` staleness correction (same audit method as #595 -- every PAYWALLED/RESTRICTED source status was stale (all acquired since), but the underlying ~216+38-row Russell backfill gap and the Aldrin/U0L1/VISIT precision dead-ends are still real; follow-on backfill task identified, not executed); #597 for 4 more Ross-group cislunar MMR/heteroclinic papers found via a ross.aoe.vt.edu review, filed + fully deep-mined (CLOSED 2026-07-15: 4 CorpusAnchors registered, digest notes rewritten with full findings); #598 for the resonance_network.py reproduce-before-trust data-gap follow-up found during #597's mining pass (CLOSED 2026-07-15: real seed + wrong-Jacobi-constant bug found and fixed, not just a doc gap); #599 for the Neptune Triton-Proteus retrograde-orbital-sense capability gap found while scoping a Fable discovery-strategy pass (CLOSED 2026-07-15: both same-sense bugs fixed, clean negative -- Proteus GM too small for a useful bend); #600 for the Uranian 3-moon-sequence extension of #563 (CLOSED 2026-07-15: clean negative, near-miss residual 0.0531 km/s just outside the gate); #601 for the #582 niching-GA sweep at Earth-Moon mu (CLOSED 2026-07-15: clean negative, 0/81 asymmetric across the 4 trusted MMRs); #602 for the #590-deferred rigorous differential-correction continuation between cluster 43 and Family J (CLOSED 2026-07-15: genuinely inconclusive on a hard same/distinct verdict -- neither endpoint is close to an exact periodic orbit, a model-mismatch finding, not a tooling gap; real evidence favors same-family over distinct-branch but this is preponderance-of-evidence, not proof); #603 for a Sun-Neptune transient-capture quasi_cycler screen extending #535/#557 (DECLINED 2026-07-15 by user decision after the plan found both a ~1649-2473yr window departure and no real anchor object exists -- see the plan doc; no build, no run, no registry entry); #604 for a V∞/Tisserand compatibility gate check on the #500->#318 Keplerian-map-chaining idea before committing to the multi-day build (CLOSED 2026-07-15: INCOMPATIBLE, real Jupiter-arrival V∞ is ~an order of magnitude above the map's valid regime, do not build the chain); #605 for a broad, genuinely creative discovery-strategy pass -- new methods, new objects/systems, innovative search strategies, distinct from the tactical unswept-body passes (CLOSED 2026-07-15: ranked 5-item shortlist produced, diagnosed the family-selection/basin-wall as the real bottleneck); #606 for #605 shortlist item 1, variational/least-action seedless periodic-orbit discovery (CLOSED 2026-07-16: built, positive-controlled, and crossed a real documented basin wall -- #556's L1 quasi-halo -- to near-machine precision; re-verified live by the coordinating session); #607 for #605 shortlist item 2, triple/quadruple small-body multi-moon systems (CLOSED 2026-07-16: clean negative, 0/97,664 across Sylvia/Eugenia/Kleopatra/Elektra -- mass-limited, same failure mode as Amalthea/Triton; Lempo-Paha-Hiisi excluded structurally); #608 for #605 shortlist item 3, a generative ML seed model trained on this project's own runlogs (CLOSED 2026-07-16: partially/meaningfully viable POC, ~12.25x physically-sane convergence lift over a uniform baseline, reproduced bit-for-bit); #609 for #605 shortlist item 4, hierarchical cycler-of-cyclers piloted via Mars Phobos-Deimos (CLOSED 2026-07-16: clean negative, step 1 closed -- both moons individually too low-mass for a useful bend, step 2 hierarchical pilot cannot proceed); #610 for #605 shortlist item 5, certified non-existence via interval arithmetic (CLOSED 2026-07-16: POC certified -- upgraded the #599 Proteus entry's physical-bend sub-gate to a continuum-strength interval-arithmetic non-existence certificate over two boxes, numbers reproduced live); #611 for a #606 follow-up pointing the new seedless corrector at #538/#544's QBCP wall (CLOSED 2026-07-16: positive control passed, crosses the wall to machine-precision agreement with an independent multi-shooting corrector; does NOT solve the actual named torus target, an explicit scope boundary; picked up from a prior agent killed after a 3.5h silent stall caused by a real `scipy` `max_nfev`/`method="lm"` bug, now fixed); #612 for a user-approved `#611` follow-up extending the seedless spectral method to `#556`'s parking-lot EM L1 quasi-halo quasi-periodic torus (CLOSED 2026-07-16: built a 2D pseudospectral torus corrector with analytic Jacobian and no forward integration in the search -- confirmed the amp>0.01 GMOS wall is a stroboscopic-shooting fragility, parent monodromy spectral radius ~1540; L2 positive control reproduces the GMOS torus, L1 continuation crosses the wall cleanly to GMOS-amp-equivalent ~0.02 with independent closure ~1e-8 where GMOS times out; explicit scope boundary -- crosses the CONVERGENCE wall but does NOT reach O&B's 0.2739, which is energy-pinned-unreachable at C=3.15 per #555, a physics fact not a corrector limit; `qp_tori.py` unmodified, no catalogue writeback); #613 for a `#612` follow-up mapping L1 rotation number vs. Jacobi constant to find whether O&B's 0.2739 is reachable at any energy, and if so wiring the result into the `#522` linking-number screen for the final reproduction attempt (dispatched 2026-07-16); #614 next-unused):**
- **#512** — (n_em, n_se) Resonance Sweep: Run sweep driver and build analytic wrap table for #411 cross-system cycle. (Resolved)
- **#513** — R52-U Recovery: Recover R52-U from sourced Braik-Ross initial conditions to partially flip the C32-dominance gate. (Resolved)
- **#514** — NAIF Kernel-Freshness Checker: Build monthly workflow and document NAIF kernel freshness. (Resolved)
- **#515 ✓ RAN, EMPTY (2026-07-01) — header corrected 2026-07-15**: 3D Lift Framework: cross-system
  cycle corrector over fixed out-of-plane amplitudes, physically scaling system lengths and aligning
  Floquet signs. "(uncommitted working-tree)" was stale — true only briefly on 2026-07-01/02; the
  script (`scripts/run_515_cross_system_3d_search.py`) and its result are both committed (`42ca41e`);
  see `data/empty_regions.jsonl` for the registered negative.
- **#516 ✓ RAN, EMPTY (2026-07-01) — header corrected 2026-07-15**: Multi-Revolution 3D Patched
  Search: n_em/n_se > 1 to bypass the single-revolution phase-closure wall. "(uncommitted
  working-tree)" was stale, same as `#515` — `scripts/run_516_multirev_3d_search.py` is committed
  (`42ca41e`), and `data/empty_regions.jsonl`'s `cross-system-se-em-3d-multirev-patched-cr3bp-
  2026-07-01` entry records the actual result: 24 points evaluated, 0 closed (the added
  revolution knobs shift the residual floor but don't close it, consistent with the planar `#411`
  1-DOF-obstruction finding). This was propagated as "never run" by a same-day discovery-strategy
  planning pass (2026-07-15) that trusted this bullet's stale caption instead of checking
  `empty_regions.jsonl` directly — a reminder that even careful cross-referencing can still miss a
  stale CAPTION on an otherwise-accurate-looking bullet; the commit message itself
  (`42ca41e`) states plainly "All three ran to completion; none closed."
- **#517 ✓ RAN, EMPTY (2026-07-01) — header corrected 2026-07-15**: Asymmetric/Mixed Libration Pairs
  in 3D: EM-L1<->SE-L2 and EM-L2<->SE-L1 crossings. Same stale-caption pattern as `#515`/`#516` —
  `scripts/run_517_asymmetric_3d_search.py` is committed (`42ca41e`); `data/empty_regions.jsonl`'s
  `cross-system-se-em-3d-asymmetric-patched-cr3bp-2026-07-01` entry records the actual result: 48
  points evaluated, 0 closed (legs solve spatially but the time-consistent phase residuals don't
  close simultaneously in the grid scanned).
- **#518** — 3D BCR4BP Continuation: Generate coherent 4-body seeds using 3D patched orbits (uncommitted working-tree).
- ✓ Resolved (2026-07-07) **#519** — VEM Multi-Synodic Low-Thrust Signal Search (commit 5ff3105): Ran scan_309_low_thrust_vem.py. The sweep returned 2 converged candidates, both of which are Sims-Flanagan infeasible due to thrust limit constraints.
- **#520** — Comprehensive 3D Sweep: 8,640-point grid search over `scripts/run_520_comprehensive_3d_search.py` (uncommitted working-tree; dispatched by a separate agent, was mis-numbered #516 on arrival — renumbered). ABORTED 2026-07-02 after 12+ hours with zero output — see DELTA above; not a negative result.
- **#521** — Pre-flight Search Gate: mandatory `preflight_search()` checked against both negative registries
  (`data/negative_results.yaml` + `data/empty_regions.jsonl`) plus task-number collision plus a required
  timing-pilot/checkpoint declaration for any grid above N points, enforced by an AST-based ratchet requiring
  every `scripts/run_*.py` to call it (proposed 2026-07-02 independent review; not yet built).
  **#521 PHASE 1 DONE (2026-07-02, this session):** migrated all 13 `data/negative_results.yaml` entries into
  `data/empty_regions.jsonl` (63 rows now, was 50), each as a validated `EmptyRegionReport` (genome/corrector/
  capability_tags/git_sha + bounded `search_extent` + `prune_gates`), so `should_sweep()` can actually gate them.
  Every physical_reason/resweep_condition/ad-hoc field (`audit_425`, `result_388`, `determination_388`, etc.)
  from the original prose was preserved, either folded into `interpretation` verbatim or kept as a structured
  `reverification` entry. No genuine duplicates found between the two registries or within the 13 (all are
  distinct region+method combinations). One provenance gap was FOUND, not fabricated over: for
  `branch_C32_C_3.1774` (#392), the committed JSONL verdict files at the script's default output path
  (`data/branch_392_v1..v4_verdict.jsonl`) actually hold a *different* candidate's numbers
  (`branch_C11a_C_3.1107`, identical to the #393 `branch_C11a_b0` row) — almost certainly overwritten by a
  later `--target C11a` re-run reusing the same default paths. The V2/V4 gauntlet numbers for
  `branch_C32_C_3.1774` therefore trace only to the prose in
  `docs/notes/2026-06-19-392-floquet-low-amplitude.md`, not to a surviving machine artefact; this is flagged
  in that entry's `interpretation` field rather than silently presented as equally well-sourced. `data/
  negative_results.yaml` is kept (not deleted) as a frozen human-readable index: every entry now carries a
  `migrated_to_empty_regions_jsonl: <region_id>` pointer, and a header comment marks it closed to new entries
  (append to `empty_regions.jsonl` via `append_empty_region()` going forward). Verified query-ability with a
  real `should_sweep()` call against the migrated #405/#411 SE<->EM region (added as
  `tests/data/test_method_capability.py::test_real_registry_skips_equal_capability_rerun_of_cross_system_se_em_search`)
  — an equal-or-weaker re-run of that search, the exact pattern that produced #515-517 this week, is now
  mechanically skipped. `uv run pytest tests/data tests/search -q` stays green.
  **STALE TEXT CORRECTED 2026-07-15 (full-file audit): Phase 2 is NOT open — it was built the
  same day.** See the "#521 PHASE 2 DONE (2026-07-02)" note earlier in this file (the AST ratchet
  `tests/scripts/test_scripts_call_preflight.py` is real, live, and independently re-confirmed
  passing multiple times this session as recently as 2026-07-14/15). This line was simply never
  updated after Phase 2 landed later the same day — left as a historical artifact, not a current
  gap.
- **#522** — Coherent-Model Torus-Connection Search: reformulate SE<->EM cross-system closure as BCR4BP/QBCP
  whiskered-torus connections instead of patched-CR3BP grids, dissolving the 1-DOF phase-closure obstruction by
  construction; re-scopes #518 (proposed 2026-07-02 independent review; not yet built).
  **SCOPING PASS 2026-07-03 (this session): citations grounded (both real), a genuine prerequisite gap found
  before any torus-connection code should be written.** Two independent audits (fresh forks, "check before
  build" discipline) plus direct PDF reads of both grounding papers:
  - **Citations verified real and on-topic** — the first audit's "not found in corpus" was correctly read as
    unacquired, not fabricated, but per [[feedback_ground_citations_against_content]] this needed independent
    confirmation before trusting them, given the project's prior same-author citation-collision incident
    (Hernandez-Jones-Jesick 2017). WebSearch + direct PDF fetch confirm: Kumar/Anderson/de la Llave, SIAM J.
    Appl. Dyn. Syst. 24(1):219-258 (2025) = arXiv:2109.14814v2 (title matches exactly); Owen & Baresi,
    Astrodynamics 8:577-595 (2024), DOI 10.1007/s42064-024-0201-0. Both PDFs acquired, filed in the private
    corpus, and digested — see [[2026-07-03-digest-kumar-anderson-delallave-2025-whiskered-tori-connections]]
    and [[2026-07-03-digest-owen-baresi-2024-knot-theory-heteroclinic]].
  - **The 1-DOF obstruction, precisely located**: `genome/cross_system_cycle.py:48-104`'s `FrameBridge` already
    does an EXACT SE<->EM coordinate transform (real rotation+translation+dimensionalization, gated by
    round-trip + Moon-position tests) — the frame bridge itself is not the problem. The problem is that the
    relative SE-vs-EM phase angle `theta` is a *searched Newton unknown* in `correct_cross_cycle`
    (`cross_system_cycle.py:1346+`), not a state variable with its own consistent dynamics — the SE and EM legs
    are still propagated under two SEPARATE autonomous CR3BPs with no shared clock tying them together, which
    is exactly why the closure residual has an irreducible angular component (the same wall #405/#411/#496/
    #515-517 all independently hit).
  - **Genuine QBCP does not exist and is a real prerequisite, not a coding task**: `core/bcr4bp.py`'s own
    docstring (lines 14-30) states outright that Andreu's coherent QBCP needs 8 Fourier tables `alpha_i(theta_S)`
    that are "NOT in the in-repo digest" and that harvesting them is "a documented *future* step" — this is a
    literature-data-acquisition gap, not a missing function. The existing `core/bcr4bp.py` is EM-centric only
    (`l=1` = Earth-Moon distance, Sun at `a_S=388.8` EM-radii, confirmed by direct read of the normalization
    constants) — SE-scale structures are numerically degenerate in this frame, so it cannot host a genuine
    SE<->EM connection as-is either.
  - **Recommended path, NOT yet executed**: do not attempt the full whiskered-torus/knot-theory machinery on top
    of an unresolved cross-system clock — that would repeat the exact "sophisticated method on an underspecified
    glue" failure pattern #405/#411/#496 already hit. Two sequenced options, needs a decision before code is
    written: (a) the full-cost path, now allocated as **#533** — acquire + digest Andreu's 8 `alpha_i` Fourier
    tables and build a genuine QBCP module (real multi-week acquisition+build, matches the "frontier work is
    multi-week" pattern already established for Track-A axes); or (b) a cheaper proxy — anchor `theta` to the Sun's synodic angle already
    present as a state in the existing EM-frame `core/bcr4bp.py` (i.e. reframe the Sun's phase, which the BCR4BP
    already tracks consistently, as the shared SE/EM clock, avoiding a second base model) — unvalidated,
    proposed as a positive-control-gated spike rather than assumed to work. EITHER WAY, the actual torus/manifold
    connection-search code (Owen & Baresi's linking-number method is the recommended first build — see its
    digest — reusing the existing `genome/qp_tori.py` GMOS torus generator) should first be validated as a
    SINGLE-SYSTEM positive control (reproducing Owen & Baresi's own Earth-Moon-internal example) independent of
    the cross-system clock question, before attempting genuine SE<->EM.
  **PHASE 1 BUILT 2026-07-03 (same session): the single-system linking-number screening machinery is complete
  and lint/mypy/test clean, per the recommended path above.** Four new modules, each independently tested before
  wiring together (the same "generic primitive + closed-form positive control, then CR3BP wiring" discipline
  that worked for #524):
  - `search/linking_number.py` — Owen & Baresi's fan-triangulation + signed segment/triangle crossing count,
    using the standard Moller-Trumbore ray-triangle test (not a byte-for-byte transcription of the paper's own
    scanned-PDF inside-triangle formula, which does not match the standard barycentric same-side test closely
    enough to trust). Validated against a textbook Hopf-link pair of unit circles (linking number exactly +-1,
    sign flips under orientation reversal) and unlinked/nested/side-by-side coplanar pairs (linking number 0).
    Found and documented a real degenerate-geometry edge case: a curve piercing exactly through the other's fan
    centroid over-counts against multiple shared-vertex triangles.
  - `search/torus_map_contours.py` — periodic-grid marching squares (torus angles wrap at 2*pi, which a plain
    image-processing contour tracer does not support) with the standard ambiguous-saddle tie-break. Validated
    against a closed-form circular level set (arc length within 0.05 of the exact 2*pi) and a periodic-wraparound
    case. Found and fixed a real degenerate-vertex bug (grid value exactly equal to the level emits spurious
    zero-length segments).
  - `genome/qp_torus_manifold.py` — the missing per-point "Floquet matrix" step `qp_tori.py` (#290) never
    computed: STM of the flow over one stroboscopic period at any torus point, stable/unstable eigenpair
    extraction with continuity-based sign-fixing, and a manifold-endpoint grid generator (perturb + propagate to
    a surface of section). Validated against the CR3BP monodromy's reciprocal-pair eigenvalue property (a
    structural fact of Hamiltonian/symplectic STMs already documented in `search/bifurcation_detector.py`, not a
    value this module computed), reusing the EXISTING sourced #299 Earth-Moon Neimark-Sacker bracket fixture
    (not deriving a new seed).
  - `genome/qp_torus_heteroclinic.py` — wires the three pieces above into `scan_linking_number()` (sweep a
    scanning variable, track linking-number sign changes) and `sign_change_locations()` (Owen & Baresi's own
    initial-guess extraction rule: average the D-values before/after a sign change). This is a SCREEN only —
    flags candidate connection locations, does not run the differential-correction refinement Owen & Baresi's
    method also needs (out of scope for Phase 1).
  - 18 tests total, all passing; validated as a mechanical wiring smoke test (self-consistency: stable vs.
    unstable manifold of the SAME sourced torus) — **NOT yet the real positive control.** Reproducing Owen &
    Baresi's actual Earth-Moon quasi-halo<->quasi-halo result (`mu=0.012153643`, `C=3.15`, 4 connections, their
    Sec 4.1.1 — see [[2026-07-03-digest-owen-baresi-2024-knot-theory-heteroclinic]] for the full sourced numbers)
    needs sourcing TWO DISTINCT quasi-halo tori at their specific published L1/L2 latitudinal frequencies
    (0.2739 / 0.02163) — a real, not-yet-attempted sub-task (likely needs a targeted continuation/search over the
    existing halo-family + Neimark-Sacker-bracket infrastructure to land on those exact frequencies, not just
    "any" Neimark-Sacker bracket like the #299 fixture used for the mechanical tests above). This is the next
    concrete #522 step if pursued further; the cross-system SE<->EM question (path (a) vs (b) above) remains a
    separate, deferred decision.
  **SPLIT 2026-07-03 (user request): #522 is now considered DONE at Phase 1 scope** (the tested, committed
  single-system linking-number screening machinery). The actual single-system connection-detection attempt (real
  L1/L2 tori, manifold globalization, transit-branch problem) is a DIFFERENT kind of work — open numerical
  research on a specific hard sub-problem, not primitive-building — and is split off as **#534** below. The
  cross-system SE<->EM question (path (a) `#533` vs path (b) synodic-angle proxy) remains #522's own open
  decision point if/when the cross-system goal is resumed; neither is required to call #522 Phase 1 complete.
- **#534** — Single-System Torus Heteroclinic Connection Search (Earth-Moon L1<->L2 quasi-halo transit branch)
  (allocated 2026-07-03, split off #522's Phase 2 exploration at the user's request). Goal: actually run #522's
  built-and-tested linking-number pipeline on a REAL, genuinely isoenergetic Earth-Moon L1/L2 quasi-halo pair and
  report a real connection count — not a mechanical self-test.
  ✓ Resolved (2026-07-07): Implemented parallel grid computation and start-of-row sign alignment to resolve the sign flip blocker, and successfully executed standard and advanced searches via scripts/run_534_torus_connection.py (commit 98a1c57).
  **⚠ STAMP (2026-07-10, #548): this 0-connection result is METHOD-INVALID — DO-NOT-CERTIFY, not an
  empty region.** Two reasons established by #548: (a) #534's committed NRHO seeds correct to
  **C=3.045**, NOT Owen & Baresi's C=3.15 (verified `cr3bp.jacobi_constant`), so this run never
  tested the paper's positive-control configuration; (b) the qp_tori/qp_torus_heteroclinic
  linking-number pipeline was SHELVED by #548 after a reframed positive-control sweep at the
  achievable common energy band found zero sign changes and never produced a nonzero torus-level
  linking number. Do not cite #534's 0-connection as evidence of anything about Earth-Moon
  connections. See `docs/notes/2026-07-10-postmortem-548-linking-number-pipeline.md`.
  **PHASE 2 EXPLORED 2026-07-03 (same session as #522's build, exploratory scripts not yet committed): real
  L1/L2 tori built at matched Jacobi, but exact frequency reproduction abandoned as impractical, and manifold-
  transit search hit a genuine open sub-problem.**
  - Found TWO independently-sourced, ALREADY-VALIDATED Earth-Moon halo seeds elsewhere in this codebase (not
    #299's family, which tops out at `C=3.148` and is not confirmed to be a halo family at all): the L1 southern
    halo seed task #304 uses (`x0=0.824024728136525, z0=-0.054501847320725, ydot0=0.164671964079122`, chained
    from Howell 1984 / Breakwell-Brown 1979), and the L2 seed `search/nrho_continuation.py` reproduces from
    Koblick 2023 AMOSTECH Table 4 Np=1 (`x0=1.023731, z0=0.183250, ydot0=-0.106950`). Re-corrected the L1 seed at
    Owen & Baresi's EXACT `mu=0.012153643` (not this codebase's own slightly different `EM_MU` constant) and
    found it lands at `C=3.150155` — essentially exactly the target energy already. Continued the L2 seed in
    `x0` (natural-parameter stepping via `search/nrho_continuation.correct_symmetric_nrho`, through a real fold
    in `C(x0)` that a naive fixed-step continuation initially missed) to `C=3.150098` — a genuinely ISOENERGETIC
    pair (`delta C ~ 6e-5`), the physically essential precondition Owen & Baresi themselves state ("heteroclinic
    connections can only exist between orbits of the same Jacobi integral").
  - Both halos have a genuine complex-conjugate Floquet pair sitting essentially exactly on the unit circle
    (`|lambda|=1.000000` to printed precision — a normal feature of STABLE halo family members, not a rare
    bifurcation point) and both correct into genuine QP tori via `correct_qp_torus` with excellent residuals
    (invariance `7.9e-9` / `6.0e-8`, independent closure `1.4e-8` / `1.1e-7` — both far inside the existing
    `1e-5`/`1e-3` published-practice gates). Resulting `omega_trans`: L1 `-0.169934`, L2 `-0.035640` — same order
    of magnitude as Owen & Baresi's published `0.2739` / `0.02163` but NOT an exact match.
  - **Decision: abandoned chasing the exact published frequencies as impractical, not as a failure.** The paper
    does not give enough explicit formula/normalization detail (from the pages read) to unambiguously
    reverse-engineer their specific "latitudinal frequency" convention, and their exact seed orbit/amplitude is
    unpublished — matching their literal numbers would require guessing unpublished implementation details, which
    is a worse use of effort than validating the pipeline on genuinely correct, independently-sourced EM 3D data
    at the one number that IS fully specified and reproducible (the shared Jacobi constant).
  - Ran `build_manifold_grids` + `scan_linking_number` end-to-end on the real L1/L2 pair (surface `x=1-mu`,
    `t_max=8` nondim TU): the L2 stable manifold crossed the surface cleanly (50/100 sample points), but the L1
    UNSTABLE manifold along the checked eigenvector direction (BOTH signs tried) never crossed `x=1-mu` at all —
    direct propagation shows it oscillates near the L1 region then escapes toward NEGATIVE x (past Earth), not
    through the Moon's realm toward L2. **This is a real, expected CR3BP phenomenon, not a pipeline bug**: not
    every unstable manifold branch of an L1 orbit is a "transit" trajectory that threads through the Moon's Hill
    sphere toward the L2 realm — distinguishing transit from non-transit branches is itself a known, nontrivial
    classification problem in CR3BP research (Conley-McGehee), not something to force through with a quick sign
    flip. Finding the genuine transit branch (if the correct linking eigen-direction even exists at this specific
    C/orbit pair) is the concrete remaining blocker before a genuine connection count can be reported.
  - **Not yet committed**: this Phase 2 exploration lived in scratchpad scripts only, to avoid committing
    half-verified numerical exploration code. If resumed, the concrete next steps are: (1) sample BOTH real
    eigen-directions' sign combinations and a spread of `t_max`/`eps` values systematically rather than 2 ad hoc
    checks, (2) consider that the genuine L1-to-L2 transit branch may require an intermediate low-energy transit
    through the interior/exterior realms (the classical Koon-Lo-Marsden-Ross patched-manifold picture) rather
    than a single direct unstable-to-stable connection at this specific `C`, and (3) once ANY genuine crossing
    grid is obtained for both branches, `scan_linking_number` itself is already built and tested and needs no
    further work to report a real linking-number sequence. Not yet resumed.
- **#523** — Co-Orbital-Exchange Cyclers: search horseshoe/tadpole/quasi-satellite topologies (Janus-Epimetheus,
  Earth quasi-satellites, Mars Trojans) whose repeated encounters need no flyby bend, sidestepping the
  mass-deficit no-go theorem (proposed 2026-07-02 independent review; not yet built).
  **ATTEMPTED 2026-07-02 (this session): STRONG BUT INCOMPLETE evidence -- NOT registered as a clean negative.**
  First settled the encounter criterion in writing (required per the review's #339-style-criterion-trap warning):
  a certified periodic orbit counts as a genuine co-orbital "encounter" with Earth if it passes within Earth's
  Hill-sphere radius during its cycle (same standard used for #527's Jupiter analysis). Checked the Jesick 2019
  Mars-Trojan paper first (already digested, `docs/notes/2026-06-17-digest-jesick-2019.md`) and correctly
  excluded it from this task: its L4/L5 relay orbit "never encounters Earth or Mars at close range" per the
  digest's own text -- it's a stationary loiter slot, not a repeated-encounter transport structure, so it's a
  separate, already-identified `quasi_cycler` admission opportunity (V0-by-publication, item 6 of that digest's
  action list), not #523 work. Built `scripts/run_523_earth_coorbital_search.py` targeting Earth quasi-satellite/
  horseshoe orbits instead (a real class: de la Fuente Marcos & de la Fuente Marcos 2018, Astron. Nachr., already
  digested, gives sourced Keplerian elements for 13 real Earth co-orbital "Arjuna-class" objects). Seeded from
  2006 RH120 (an object that WAS an actual observationally-confirmed transient Earth minimoon, Jul 2006-Jul 2007)
  -- direct integration of this seed, before any search code was written, independently reproduced its known
  qualitative behaviour (an initial close quasi-satellite episode reaching HALF the Hill radius, then a
  transition to a wider horseshoe libration), unplanned confirmation the model is physically correct.
  Positive control passed. **Certified 3 distinct real periodic horseshoe orbits** (residuals ~1e-12, periods
  31.4/44.0/81.7 nondim ~= 5.0/7.0/13.0 years) near this seed. **All three show the same #527 pattern**: minimum
  Earth approach is 9x/16x/22x Earth's Hill radius -- the certified, truly periodic orbit family stays well
  outside gravitational relevance, while only the TRANSIENT (non-periodic) drift phase gets genuinely close
  (matching the real 2006 RH120 case). **Coverage is incomplete**: only ~2 of the intended 6 Jacobi values were
  actually swept (C=2.9990 thoroughly + one point at C=2.9999) before the run was stopped -- the certification
  corrector (`correct_general_periodic` via a chained-refinement loop) costs 60-100+s per candidate here (much
  more than #527's Sun-Jupiter case), making the full 1,890-point band impractical at the current per-candidate
  cost; two design inefficiencies were found and fixed mid-session (near-duplicate coarse candidates converging
  to the same certified orbit without a post-certification dedup; non-converging candidates costing MORE time
  than converging ones under the original 5-pass/max_iter=100 refinement loop) but the remaining per-candidate
  cost is still too high for a full sweep in reasonable time. Do NOT trust a "no encounter" claim for the full
  band from this alone -- the sample is real and consistent but small. Resweep condition: NOT more brute-force
  grid points at this per-candidate cost; needs either a genuinely faster corrector (e.g. warm-starting
  certification directly from the coarse candidate's already-decent residual with a tighter trust region, or
  #524's planned continuation/deflated-Newton primitive) or acceptance that the 3-orbit sample is sufficient
  preliminary evidence to deprioritise this specific (stable-periodic-orbit) angle in favour of the
  transient/quasi-periodic drift phase, which is where the real close approaches actually live (same conclusion
  #527 reached for the Hilda case) -- but that phase is not "periodic" in the catalogue's `cycler` sense and
  would need the `quasi_cycler` epoch-locked framing instead. Janus-Epimetheus and the broader Mars-Trojan-as-
  transport question remain untouched by this session's work.
  **REWORK PLAN 2026-07-03 (not yet executed) — use #524's now-built continuation primitive instead of a
  brute-force resweep.** `scripts/run_523_earth_coorbital_search.py`'s actual bottleneck is per-candidate
  CERTIFICATION (`_certify_with_refinement`, chained `correct_general_periodic` passes from a coarse DA/HOTM
  grid guess, 60-100s each), not the coarse enumerator pass itself (~0.05s/point). The coarse enumerator's job
  — locating an approximate fixed point at a GIVEN fixed `c_target` — is a poor fit for expensive independent
  re-certification at every new `(C, x0, xdot0)` grid cell when 3 orbits are ALREADY certified nearby in family
  space. Plan: treat the family as a co-dimension-1 curve in the 3-vector `(x0, xdot0, C)` — exactly
  `search/pseudo_arclength.py`'s target case (`M=N-1`: the 2-component Poincare-return residual against 3 free
  unknowns with `C` now free instead of fixed) — and `continue_curve()` OUTWARD from each of the 3 already-
  certified orbits in both `C` directions, replacing "independently re-certify every coarse grid cell from
  scratch" with "cheap predictor-corrector steps from an already-converged point," which should need far fewer
  Newton iterations per new family member than a fresh multi-pass certification from a coarse guess. Concrete
  steps: (1) wire a `residual_fn(z)` wrapping the same Poincare-return condition `_certify_with_refinement`
  drives to zero, parameterized by the free `(x0, xdot0, C)` vector; (2) seed `continue_curve()` from each of
  the 3 certified orbits with a numerically-estimated `jacobian_fn` (or let it fall back to the built-in finite-
  difference Jacobian); (3) walk both directions in `C` across the original `C_BAND` intended range; (4) at each
  new certified point, run the SAME Hill-sphere-encounter classification the original script already has, no
  changes needed there. Expected win: if a single continuation step costs closer to a handful of Newton
  iterations (~1-5s) rather than the original chained-certification's 60-100s, the SAME 6-Jacobi-value coverage
  the original scoping wanted becomes a small-minutes job instead of the ~32-52 CPU-hour brute-force estimate
  for the untried remainder of the 1,890-point grid. Not yet built or timed — this is a plan, not a result;
  the actual per-step continuation cost on this specific residual is unmeasured and should get its own timing
  pilot (per the #521 preflight discipline) before committing to a full re-run.
  **REWORK BUILT + RUN 2026-07-03 (same session): the plan above worked, and the finding now covers the full
  intended band densely.** Built `search/cr3bp_general_periodic_free_c.py`: extends
  `cr3bp_general_periodic.py`'s existing analytic STM Jacobian `dR/d(x0,xdot0)` with a third, EQUALLY analytic
  `dR/dC` column (free — reuses the same STM already computed, no extra propagation), turning the 2x2
  fixed-Jacobi residual into the co-dimension-1 curve `search/pseudo_arclength.py` targets. Validated the new
  `dR/dC` formula against finite differences at the existing positive-control seed before trusting it (an
  internal-consistency check, not a value invented and self-checked) — also cross-checked the `(x0,xdot0)`
  columns independently as a sanity check that the shared STM-projection machinery is used correctly. 4 new
  tests, all passing.
  Then rewrote the sweep as `scripts/run_523_earth_coorbital_continuation.py`: certify ONE orbit from the SAME
  2006-RH120-derived positive control as the original script, then walk `continue_curve()` in both directions
  across `C in [2.9980, 3.0020]` (covers and exceeds the original `C_BAND=[2.9990,3.0010]` with margin).
  **Measured, not estimated: ~1.4s/continuation-step** (dominated by two STM propagations over the ~80-nondim-TU
  horizon) vs. the original 60-100s/candidate — confirms the plan's premise directly. The actual run certified
  **120 points spanning C=2.9980 to 3.0019940 in 224.5s wall-clock** — dense, essentially complete coverage of
  the intended band (vs. the original attempt's 3 orbits across ~2 of 6 intended values), a ~50-70x per-point
  speedup that also avoids continuation's main risk (re-discovering already-covered family members, the way
  independent grid re-sampling does). **Result: EVERY certified point's minimum distance to Earth (0.0488-0.1032
  nondim) is 4.9-10.3x Earth's Hill radius (0.01000 nondim) — zero encounters**, extending the original
  3-orbit finding to essentially complete coverage of the intended Jacobi band. Registered as a genuine clean
  negative in `data/empty_regions.jsonl` (`sun-earth-coorbital-horseshoe-qsat-continuation`, distinct region_id
  from the original incomplete attempt since the METHOD capability differs materially) — round-trip and
  `should_sweep()` gating verified live (a same-capability re-sweep correctly returns `False`). **Same caveat as
  the original attempt, unchanged by this rework**: this negates the STABLE PERIODIC ORBIT mechanism
  specifically; it does NOT test the transient/quasi-periodic drift phase (where the real close approaches
  live, per 2006 RH120's own qualitative behavior and the #527 Hilda finding) — that remains a genuinely
  different search needing the `quasi_cycler` epoch-locked framing. Full ratchet
  (`tests/data tests/search tests/genome tests/scripts`) verified green.
- **#524** — Continuation + Deflated-Newton Search Primitive: replace fixed-grid + independent-Newton-per-cell
  sweeps with pseudo-arclength continuation and deflated Newton as the default search method (proposed
  2026-07-02 independent review).
  **BUILT 2026-07-02/03 (this session), scope corrected mid-build — see below.** Two independent
  half-tasks, not one build:
  - **Continuation half: audit found it already exists, extracted rather than built new.** Before writing
    code, checked for existing pseudo-arclength continuation and found FOUR specialized implementations
    already in the tree (`cr3bp_jacobi_arclength.py` for `(x0,C)`, `mu_continuation.py` for `(x0,C,mu)`,
    plus `qp_tori_arclength.py`, `narc_continuation.py`) — each reimplementing the same
    predict-correct-tangent machinery for its own hardcoded residual. Built
    `src/cyclerfinder/search/pseudo_arclength.py` as a generic, CR3BP-agnostic extraction
    (`continue_curve()`, co-dimension-1 curves only — `M=N-1` constraints, the well-posed unique-tangent
    case; deliberately NOT generalizing `mu_continuation.py`'s ambiguous 2-surface tangent-selection case).
    Positive control: a closed-form fold (unit circle `x^2+y^2-1=0`) — the walk sails through the `x=1`
    fold (where a natural-parameter fix-x-solve-y continuation cannot turn) while staying on the circle to
    1e-9, and a full `2*pi`-arclength loop returns to its start. 6/6 tests, lint/mypy clean
    (`tests/search/test_pseudo_arclength.py`, commit `6f22bd7`).
  - **Deflation half: genuinely absent, built new.** Searched the codebase for existing deflated-Newton /
    basin-repulsion root enumeration; the only two "deflat" string hits (`core/lambert.py`'s deflation-ANGLE
    handling, `search/bifurcation_detector.py`'s `_deflated_determinant`, a Doedel-1991 bifurcation TEST
    FUNCTION for detecting folds along a known branch) are unrelated concepts — this is a real capability
    gap. Built `src/cyclerfinder/search/deflated_newton.py` implementing Farrell, Birkisson & Funke (2015,
    SIAM J. Sci. Comput. 37(4):A2026-A2045): deflated system `G(u)=M(u;{u_i})*F(u)`,
    `M=prod_i(1/||u-u_i||^p+shift)`, repels a Newton iterate from already-found roots without moving other
    roots, so `enumerate_roots()` walking a seed set finds distinct co-existing roots one at a time.
    Positive controls: Farrell et al.'s own two worked examples (cubic polynomial roots `{1,2,3}`; unit
    circle/line-`y=x` intersection, two roots) — both closed-form, EXPECTED values sourced from the
    published method, not self-computed. Then a REAL CR3BP tie-in: plain (undeflated) Newton on the exact
    scalar fixed-Jacobi perpendicular-crossing residual `correct_symmetric_fixed_jacobi` already solves,
    seeded at the sourced #504 Ross-RT 2026 Pluto-Charon (3,2) anchor (`x0=-0.694376003123377`,
    `C=3.573367616904619`, `half_crossings=6`), independently reproduces that corrector's converged `x0` to
    1e-6 — an agreement check between two independent Newton implementations on the same physical residual
    and sourced data, not a fabricated golden. 4/4 tests, lint/mypy clean
    (`tests/search/test_deflated_newton.py`, commit `4f4a1ab`).
  - **Net effect:** #524 delivers one genuinely novel capability (deflated Newton) and one
    already-existing-but-now-reusable one (generic continuation); the original proposal's framing ("replace
    fixed-grid... sweeps") describes future consumers (#522 and any `scripts/run_*.py`), not something this
    task itself rewired into existing campaigns. Neither primitive has been wired into a live discovery
    search yet — that is the #522 (and #532, if pursued) follow-on work.
- **#525** — Learned Seed Generation: diffusion/generative corrector-seed model trained on the project's own
  runlog/checkpoint archive; revisits the 2026-06-11 Ozaki (arXiv 2111.11858) "below breakeven" triage as stale
  (proposed 2026-07-02 independent review; not yet built).
- **#526** — GTOC 13 + Liang Operator Ingestion: digest GTOC 13 (Oct-Nov 2025 ballistic Jovian gravity-assist
  tour competition) methods papers and encode Liang et al.'s alternating-double-cycler construction (JGCD 2024,
  DOI 10.2514/1.G008387) as a reusable genome operator (proposed 2026-07-02 independent review; not yet built).
  **DONE 2026-07-03 (this session), with one correction to the proposal text.** Part 1 (GTOC 13): the
  proposal's "ballistic Jovian gravity-assist tour competition" framing was WRONG — GTOC 13 (JPL, released
  20 Oct 2025, ran 24 Oct–17 Nov 2025) is a ballistic-flyby + optional ideal-solar-sail tour of a fictional
  single-star, 10-planet exoplanetary system ("Altaira"), not a Jovian-moon problem at all; there is no moon
  system in it. Acquired the problem statement PDF (`papers/jpl-gtoc13-team-2025-problem-description-altaira-
  system-gtoc-jpl-net.pdf`, filed/digested/indexed per the corpus policy:
  `docs/notes/2026-07-03-digest-gtoc13-problem-statement.md`, `CORPUS_INDEX.md`). Winner: THU-LAD (Tsinghua);
  101 teams registered (Team 11 "NUAA & Friends" includes Liang and Yang, the Part-2 authors). **Honest
  negative on methods papers: none exist yet** — searched WebSearch/arXiv/the GTOC13 site/the ESA GTOC-portal
  mirror; the portal lists the post-competition workshop as still "TBA", consistent with the historical
  6–18-month workshop→journal-special-issue pipeline for a competition that closed only 8 months ago. No
  cycler-relevant content in the problem itself (single one-shot 200-yr tour maximizing a score, not a
  repeating structure) — a problem-statement-only digest, as the proposal itself anticipated as a live
  possibility. Part 2 (Liang operator): independently re-verified the citation before trusting it (AIAA/ARC +
  Semantic Scholar + Unpaywall DOI lookups all agree on title/authors/venue; Unpaywall confirms closed-access,
  no OA copy exists anywhere) — real paper, no same-author collision. Found the paper was **already** filed,
  deeply mined (`2026-06-11-liang-2024-cge-triple-cyclers-mining.md`) and reproduced
  (`search/cge_scaffold.py`, `search/moon_cycler_genome.py`) under #216 — nothing duplicated. Built the
  missing reusable piece: `src/cyclerfinder/genome/alternating_double_cycler.py`
  (`analyze_near_resonance` — generalizes Eqs. 1–2's near-commensurability analysis to any 3-body chain;
  `build_alternating_double_cycler_seed` — the "switched-double-cycler" stitching operator itself, composing
  two hub-sharing double-cycler half-sequences into one repeating `MoonCyclerGenome`). Design pitfall caught
  by the positive control before trusting the code: a naive closest-fraction search returns 16/9 for Liang's
  own numbers, not their published 7/4 — fixed via a lowest-order-convergent-above-tolerance search (see
  `docs/notes/2026-07-03-alternating-double-cycler-operator.md` for the full account). Positive controls:
  exact reproduction of Liang's Eq. 1–2 (7/4, mismatch 0.7365 d) from their own Table 1; an independent
  (non-golden) cross-check against this project's JPL-SSD registry mean motions recovers the same 7:4
  resonance within 2%; exact reconstruction of the published CGCEC sequence from its two halves; a
  reusability-only (no novelty claim) demo on a Saturn chain (Enceladus-Dione-Rhea) via a newly added
  `saturn_system()` registry helper. 13/13 new tests pass; `uv run pytest tests/data tests/search tests/genome
  -q` ratchet, `ruff check`/`ruff format --check`, and `mypy src tests` all clean. No catalogue writeback; no
  novelty claim.
- **#527** — Sun-Jupiter Hilda-Band DA/HOTM Enumeration: point the already-built, zero-usage #450 DA/HOTM
  enumerator (`src/cyclerfinder/search/da_hotm_enumeration.py`, μ-agnostic, never invoked by any campaign
  script) at the Sun-Jupiter Hilda/quasi-Hilda 3:2-exterior-MMR band — a region in neither negative registry nor
  the catalogue. Positive control: recover manifold-consistent unstable POs per Guido & Efthymiopoulos,
  arXiv:2604.00679 (2026-04-01), before any negative is trusted. Cheapest of the second-pass proposals — a
  config/target change, no new build (proposed 2026-07-02 second-pass independent review; not yet built).
  **RUN 2026-07-02 (this session): CLEAN, HONEST NEGATIVE — registered in `data/empty_regions.jsonl`
  (`sun-jupiter-hilda-32-mmr-dahotm`).** `scripts/run_527_hilda_dahotm_search.py` built + run end-to-end.
  Independently re-verified the positive-control paper is real (WebFetch confirmed title/authors/abstract match
  exactly). Derived the Hilda 3:2 interior-MMR seed from first principles (Kepler III on Jupiter's sourced SMA
  -> a=3.9705 AU, matching the real, independently-documented Hilda asteroid group; back-derived GM_sun matched
  the IAU-standard value exactly -- both are hard independent cross-checks, not invented numbers). Found and
  fixed a real bug before it could cause a false negative: the #450 driver's defaults (`t_max=8.0`,
  `period_guess=2n*2.5`) are tuned for Earth-Moon-scale periods and would have silently found nothing here,
  where the natural section-return period is ~12.6 nondim units (verified by direct integration) -- the script
  builds its own correctly-scaled backend instead of reusing the driver unmodified. Positive control PASSED
  (residual 4.25e-3) before the full scan ran. Full scan (8 Jacobi values x 2 rev counts x 651-point grid =
  10,416 points, timed pilot measured 0.007-0.043 s/point beforehand per #521) certified 33 real periodic-orbit
  family members (residuals ~1e-14, well-converged) across the Hilda-resonance Jacobi band. **Verdict: the
  #527 premise does not hold for this method.** The closest Jupiter approach across all 33 certified orbits is
  1.85 Jupiter Hill radii (98.5M km vs Hill radius 53.1M km = 0.355 AU, itself matching the known literature
  value) -- safely OUTSIDE Jupiter's sphere of gravitational dominance. This is physically expected, not a
  search failure: stable 3:2-resonance libration is WHY the Hilda group avoids close Jupiter encounters, so a
  periodic-orbit enumeration of the STABLE resonant backbone structurally cannot find the close-approach
  "cycler" mechanism #527 hypothesized. The Guido-Efthymiopoulos heteroclinic/chaotic-transport structure this
  proposal was grounded in lives on the UNSTABLE manifolds near the resonance separatrix, a different dynamical
  object DA/HOTM's periodic-orbit enumeration cannot surface. Resweep condition: NOT more grid/Jacobi resolution
  on this stable-orbit search -- a genuinely different method (invariant-manifold propagation from the unstable
  separatrix orbits) would be needed to test the actual chaotic-transport cycler hypothesis; scope as a fresh
  task if pursued.
- ✓ Resolved (2026-07-06) **#528** — Low-Thrust-Native Validation Gate: build a continuous-thrust acceptance criterion for (commit cfac8e0).
  `dv_band: low_thrust_sep` (currently a schema field with zero catalogue rows and a gate that `return`s `None`
  — no acceptance window exists), reusing `core/sims_flanagan.py`. Positive control: reproduce Pascarella et al.
  2022 (AAS 22-015, already in this project's corpus) to a stated tolerance. Should land before or alongside
  #519 actually running, or a genuine low-thrust find would have nothing to grade it against (proposed
  2026-07-02 second-pass independent review; not yet built).
- **#529** — Inter-Cycler Network / Taxi-Transfer Catalogue Extension: add a schema relation for cycler *sets*
  with shared phasing/downlink cadence and taxi transfers between catalogued cyclers (reusing
  `model/score.py`'s `taxi_cost_kms`) — a published, sourced concept (Sanchez Net et al. 2022, already digested,
  `docs/notes/s1l1-target-topology-mining.md:156-263`) the current 5-class taxonomy cannot express. Deferred
  deliberately — changes what a catalogue row/relation can be, wants a scoping discussion before a sprint slot
  (proposed 2026-07-02 second-pass independent review; not yet built).
- **#530** — Unstable-Manifold Propagation for Co-Orbital/Resonance Cyclers: compute Floquet multipliers for the
  periodic orbits already certified by #523/#527 to identify genuinely unstable family members, then propagate
  their stable/unstable manifolds and test whether trajectories along them achieve genuine Hill-sphere
  encounters — the actual mechanism both source papers describe, which periodic-orbit enumeration structurally
  cannot surface. Reuses existing `search/bifurcation_detector.py` Floquet/monodromy machinery + the two
  already-built target systems (Sun-Earth co-orbital, Sun-Jupiter Hilda). Allocated 2026-07-02, motivated by
  this session's #523/#527 findings; see the Tier 6 note above for full scoping.
  **RUN 2026-07-02 (this session): COMPLETE, clean local/linear-regime negative — registered in
  `data/empty_regions.jsonl` (`coorbital-resonance-unstable-manifold-encounters`).**
  `scripts/run_530_unstable_manifold_search.py` re-certified 9 orbits (6 from #527's Hilda family, 3 from #523's
  Earth co-orbital family) from their already-known-converged (x0, xdot0) — a refinement, not a fresh search —
  computed each one's monodromy + eigendecomposition, and propagated the manifold of any unstable member.
  **Positive control PASSED universally**: `det(monodromy)=1.000000` (symplectic, Liouville-consistent) for all
  9 orbits, validating the machinery independent of any physical claim. **Stability result**: 8 of 9 orbits are
  linearly stable/marginally-stable (leading multiplier on the unit circle) — no unstable manifold exists for
  them at all. Exactly **one** genuinely unstable orbit was found (Sun-Jupiter Hilda, C=3.14, x0≈0.7615,
  λ≈1.718, period≈1.0 yr). Its manifold was propagated for a **Lyapunov-scaled** horizon (20 periods, computed
  from `ln(target_growth/eps)/ln(|λ|)`, not guessed — an earlier fixed-3-period attempt found the perturbation
  hadn't grown at all and was caught and fixed before this entry was written). With the corrected horizon, both
  manifold branches still land within 0.4% of the periodic orbit's own closest approach (3.47-3.48× Jupiter's
  Hill radius) — nowhere near an encounter. **Honest scope caveat, not glossed over**: this tests the
  LOCAL/LINEAR-REGIME manifold only (a fixed eigenvector direction scaled by the growing Floquet factor, valid
  while the true trajectory stays near the linear approximation, ~20-25 e-foldings here) — it does NOT test full
  nonlinear manifold globalization (re-correcting onto the true invariant surface at each step, e.g.
  Koon-Lo-Marsden-Ross-style), which is the technique that would actually test the Guido-Efthymiopoulos
  heteroclinic-transport hypothesis at large distances from the parent orbit. Resweep condition: NOT more
  re-certified orbits at this same (local-manifold) technique — a genuinely different, more expensive nonlinear
  globalization method would be needed to test the far-field heteroclinic hypothesis properly; scope as a fresh
  task if pursued.
  **CORRECTION (2026-07-02, same session):** "materially different and more expensive" above was written before
  checking whether this codebase already had heteroclinic-connection machinery — it does. `genome/heteroclinic_cycle.py`
  (#314) is a full Wu(A)->Ws(B) manifold-connection corrector (`correct_connection`/`assemble_cycle`/
  `crosscheck_cycle`), already validated against a REAL Sun-Jupiter heteroclinic cycle (Wilczak & Zgliczyński's
  computer-assisted proof of the Oterma comet's L1<->L2 connection, arXiv:math/0201278, at C=3.03 — inside this
  session's own Hilda C-band [2.95, 3.14]). `LyapunovNode` (the framework's node type) is a plain dataclass
  (state0/period/jacobi/unstable_eigvec/stable_eigvec) constructible directly from #530's already-computed
  monodromy/eigendecomposition data, NOT restricted to libration-point orbits despite the name (`from_libration`
  is just one convenience constructor). See #531 below — this is the actually-cheap next step, not a from-scratch
  build.
- **#531** — Heteroclinic/Homoclinic Connection Test for the #530 Unstable Hilda Orbit (allocated 2026-07-02,
  same session, superseding the "full nonlinear globalization" framing in #530's note above once #314's existing
  machinery was found). Wire #530's one genuinely unstable orbit (Sun-Jupiter Hilda, C=3.14, x0≈0.7615, λ≈1.718)
  into a `genome/heteroclinic_cycle.LyapunovNode` (a plain dataclass — construct it directly from #530's already-
  computed state0/period/jacobi/eigenvectors, bypassing the libration-specific `from_libration` constructor) and
  call `assemble_cycle(system, [node])` to test for a HOMOCLINIC connection (the orbit's own unstable manifold
  reconnecting to its own stable manifold — the degenerate n=1 case `assemble_cycle` explicitly supports).
  Track the Hill-sphere distance along the actual certified connection geometry (`correct_connection`'s matched
  crossing + the two manifold legs), not just a short local propagation — this is what #530 could not test.
  Cross-check with `crosscheck_cycle` (independent Radau re-integration, already built into the framework).
  Positive control: ALREADY EXISTS and already passes — `tests/genome/test_heteroclinic_cycle.py` reproduces the
  Wilczak & Zgliczyński Sun-Jupiter-Oterma cycle (arXiv:math/0201278) at C=3.03, inside this session's own Hilda
  band; no new positive control needs to be built, only cited. If a same-Jacobi Earth co-orbital or Hilda-family
  pair (not just the single unstable orbit) is later found, `assemble_cycle` also directly supports true
  multi-node heteroclinic cycles (`nodes` requires matching Jacobi, an autonomous-system physical constraint —
  connections only exist between same-energy orbits). Feasibility: LOW build cost — this is wiring existing,
  validated machinery to new node data, not building new manifold/connection-correction capability.
  **RUN 2026-07-02 (this session): COMPLETE — a genuine connection was found, and it's a materially STRONGER
  negative than #530 alone.** `scripts/run_531_hilda_homoclinic_connection.py` re-ran the #314 framework's own
  W-Z Oterma golden LIVE first (not cached — residual 1.854e-10, matching the published connection to 5 digits),
  confirming the machinery works today. Built a `LyapunovNode` directly from the #530 Hilda orbit (reciprocal
  Floquet pair check: λ_u·λ_s≈1.0000, sanity-passes). Explored 36 `(k_u, k_s, branch_u, branch_s)` combinations
  (the framework's Oterma-tuned defaults don't transfer to this orbit's geometry) and found a genuine, CERTIFIED,
  CLOSED homoclinic connection at `(k_u=1, k_s=1, branch_u=-1, branch_s=+1)`: residual 3.582e-08,
  `assemble_cycle` reports `closed=True`, and an **independent Radau re-integration** (`crosscheck_cycle`)
  reproduces the identical residual — a real, verified invariant-manifold object, not a numerical artefact.
  **The Hill-sphere result**: tracking BOTH manifold legs over their FULL span to the actual crossing (51.1
  nondim time units each way — the genuine test #530's short local propagation could not do) gives
  `min_dist_to_jupiter=0.23759` for BOTH legs, **exactly matching** the underlying periodic orbit's own baseline
  closest approach (3.48× the Hill radius, #530's number). The manifold diverges exponentially from the orbit in
  phase space (that's what makes it unstable) but does NOT visit a spatially different region relative to
  Jupiter over this connection — for this orbit, the instability is a phase/timing effect, not a spatial-
  excursion effect. Registered in `data/empty_regions.jsonl`
  (`hilda-c3.14-homoclinic-connection-hill-encounter`), verified queryable. **Scope note, stated plainly**: this
  rules out the SINGLE orbit's own homoclinic tangle, not the broader Guido-Efthymiopoulos resonance-hopping
  hypothesis, which needs heteroclinic connections BETWEEN DIFFERENT orbits/MMR families — a genuinely larger
  search (same-Jacobi connections across the Hilda family and/or to other resonances), not attempted here.
  A minor bug (a `t_eval` sort-direction mismatch for backward-time integration) crashed the first run right at
  the final distance check after ~9 minutes of genuine, already-successful upstream work (positive control +
  connection + crosscheck all completed and printed before the crash); fixed and re-run reusing the
  already-verified combo rather than re-exploring from scratch.
- **#532** — Multi-Orbit Heteroclinic Resonance-Hopping Search (allocated 2026-07-02, same session — the
  actual Guido-Efthymiopoulos hypothesis, not yet properly tested). #531 tested only ONE orbit's own homoclinic
  tangle (a self-connection) and found no encounter. Guido & Efthymiopoulos (arXiv:2604.00679) describe
  "resonance hopping": heteroclinic connections BETWEEN DIFFERENT periodic orbits across DIFFERENT first-order
  MMR families (2:1 interior, 3:2 interior/Hilda — where #527/#530/#531's orbit lives, 2:3 exterior) — a
  genuinely different, larger search than a single orbit's self-connection.

  Required steps: (1) find unstable periodic orbit family members in the 2:1 and 2:3 MMR families (mirroring
  #527's DA/HOTM approach + Kepler-derived seed method for the semi-major-axis location of each resonance,
  already validated this session), (2) CRITICALLY, at a Jacobi constant SHARED with an unstable member of the
  3:2 family — `assemble_cycle`/`correct_connection` hard-require exact same-Jacobi nodes (a real physical
  constraint: connections only exist between same-energy orbits in this autonomous model), which is NOT
  guaranteed to exist and was not checked this session, (3) build `LyapunovNode`s for each, (4) search
  `correct_connection` across ALL cross-family node pairs (not self-connections), exploring
  `(k_u, k_s, branch_u, branch_s)` per pair as #531 did for one pair, (5) track Hill-sphere distance along every
  certified leg.

  **Feasibility: MEDIUM, genuinely more uncertain than #531 was** — #531 succeeded partly because a single
  orbit trivially shares its own Jacobi constant; #532 needs to find (or fails to find) a shared C where
  MULTIPLE resonance families simultaneously host unstable orbits, which is open science, not just an
  engineering wire-up. #530 found only 1 of 6 sampled Hilda orbits was genuinely unstable (a ~17% hit rate) —
  the same low yield likely applies per-family here, so covering three families may need substantially more
  than 6 sample points each. Recommend a SCOPING PILOT first (per the #521 discipline): before committing to
  the full multi-family build, spend a bounded, timed search checking whether ANY shared-Jacobi unstable pair
  exists across even two families (say 3:2 and 2:1) before building out the full three-family connection search.
  Positive control: the multi-node machinery itself is already validated (`test_assemble_l1_l2_two_cycle_closes`,
  the W-Z Oterma two-node cycle, re-verified live in #531) — no new machinery-validation control needed, only
  physical-target validation once real cross-family unstable orbits are found. Not yet built.
  **SCOPING PILOT RUN 2026-07-03 (same session, per an independent second-opinion review of the plan): a real
  false positive caught and fixed before being trusted.** `scripts/run_532_resonance_overlap_pilot.py` derived a
  Kepler-III seed for the 2:1 interior MMR (`r=(1/2)^(2/3)`, `a=3.278 AU` — independently matching the real,
  documented "Hecuba gap" asteroid-belt feature, mirroring #527's own Hilda a=3.97 AU cross-check) and re-verified
  the seed-derivation formula against #527's own committed `C_SEED=3.0613` (matched to 4 decimals) before trusting
  the new 2:1 number. A coarse 8-point-per-family DA/HOTM band scan (reusing #527's enumeration + #530's
  monodromy-stability classification) found Hilda-unstable C up to 3.1413 and interior-unstable C at 3.1453 — a
  gap smaller than the 0.03 coarse step, so a finer 0.005-step rescan (`_finegrid.py`) was run in the C=[3.10,3.18]
  window. **The finer scan appeared to find a genuine 5-point overlap (C=3.14-3.16) — but this was a FALSE
  POSITIVE, caught before being reported as a finding**: at every "matching" C, both family-labeled scans had
  converged to the IDENTICAL x0 (5-decimal match) with IDENTICAL eigenvalues. Since `(x0, xdot0=0, C)` uniquely
  fixes the full state via `ydot0`, an identical x0 at the same C is mathematically the same orbit, not two
  distinct family members — the two families' DomainBoxes overlap in x (Hilda x in [0.61,0.91], interior x in
  [0.50,0.76]), and both scans' Newton correction converged onto one single periodic-orbit branch passing through
  their shared overlap region, mislabeled by which box found it rather than genuine resonance membership.
  **Fixed properly** (`scripts/run_532_resonance_family_continuation.py`): family membership is only rigorously
  defined as a connected component of the solution manifold, so this traces each family via NATURAL-PARAMETER
  CONTINUATION IN C starting from its own verified seed (warm-started Newton correction, previous step's solution
  as the next step's guess — the same pattern #534's L2 seed continuation used), with an explicit `|x0_hilda -
  x0_interior| > 1e-3` distinctness check at any shared-C match before calling it a genuine pair. **Corrected
  result: NO-GO.** Traced Hilda over C=[3.020,3.164] (74 points, 1 unstable) and the 2:1 interior branch over
  C=[3.000,3.200] (102 points, 0 unstable) — independently re-verified the interior family's "0 unstable" by
  checking every point's symplectic-det control directly (zero failures across all 101 non-seed points, leading
  eigenvalue magnitude sits at exactly 1.0000 to 10+ decimals throughout), ruling out a silent monodromy-
  computation bug masquerading as stability. With one family fully stable across its entire tested range, no
  shared-Jacobi unstable pair can exist between these two specific traced branches. **Honest scope of this
  negative**: only ONE branch per family was traced (the `ydot0`-positive, single-rev branch reachable by
  continuation from the natural circular-orbit seed) — the original box scan surfaced other 2:1-region candidates
  (x0 around 0.50-0.66) on different branches that this continuation did not trace, so this rules out the tested
  branches, not the full 2:1/3:2 resonance regions.
  **#532 COMPLETED 2026-07-04 (same session, per "finish 532?"): full-branch census closes the gap above.**
  `scripts/run_532_full_branch_census.py` ran ONE wide DA/HOTM enumeration at a single representative
  `C=3.15` across the FULL combined `x in [0.45,0.95]` range spanning both families' original domain boxes
  (not per-family), finding every distinct periodic orbit passing through this region rather than relying on
  which of two arbitrary boxes happened to find it (the source of the earlier false positive). Found exactly
  3 genuinely distinct branches at `C=3.15` (`x0=0.514664`, `0.661391`, `0.738528` — the third is confirmed
  the SAME branch mislabeled "Hilda" by the original box-overlap bug: its `x0=0.707314` at `C=3.16` matches
  the earlier bugged trace to 5 decimals). Traced each independently via natural-parameter continuation over
  `C=[2.90,3.25]`: branch `0.5147` is fully stable across its ENTIRE traced range (176 points, 0 unstable);
  branch `0.6614` has exactly ONE unstable point (`C=3.092`, traced C=[3.02,3.164]); branch `0.7385` has 8
  unstable points but its traced range only starts at `C=3.146` (a fold/boundary below that). **Verified the
  one close call by hand**: branch `0.6614`'s single unstable point (`C=3.092`) falls entirely outside branch
  `0.7385`'s traced range (which starts at `3.146`) — confirmed by directly re-deriving both branches' unstable
  points and checking the actual numbers, not just trusting the pairwise-check code. **Exhaustive pairwise
  check across all 3 branches: VERDICT NO-GO** — no genuinely distinct (`|x0_a - x0_b| > 1e-3`) shared-Jacobi
  unstable pair exists among any combination. This is now a real, well-verified negative for every branch
  passing through the shared 2:1/3:2 x-region at this energy scale — not just the two Kepler-seeded branches,
  and independently checked at the one point where the traced ranges actually overlap. Residual honest scope:
  this censuses ONE representative C (3.15) and the `ydot0`-positive branch of the y=0 return map only; a
  different census C or the `ydot0`-negative branch could in principle surface further distinct families not
  connected to any of these 3 by continuation. Given the effort already invested (coarse pilot -> finegrid ->
  false-positive catch -> two-seed continuation -> full 3-branch census with a hand-verified overlap check),
  #532 is now closed as a well-scoped, thoroughly negative scoping pilot; the full three-family
  (2:1/3:2/2:3-exterior) connection-search build described at the top of this entry is not warranted by this
  result and is not planned unless a future census at a different C or return-map branch changes the picture.
- **#535** — Transient-Drift-Phase `quasi_cycler` Search for Co-Orbital/Resonance Objects (allocated 2026-07-03,
  same session, at the user's request — the deferred follow-up #523 and #527 BOTH independently pointed at, now
  finally given its own task number instead of remaining a repeated pointer in three different entries).
  **The precise target, grounded in the catalogue's own class definition** (`docs/notes/2026-06-16-catalogue-
  scope-taxonomy.md`): `quasi_cycler` = "closes-up-to-rotation" (NOT exact periodicity), `epoch_locked=true`, a
  finite 10-15 year `validity_window`, `n_returns` in [3,15] — "cyclers-of-opportunity inside a planetary-
  alignment window." This matches #523's own positive-control object EXACTLY: 2006 RH120 (de la Fuente Marcos &
  de la Fuente Marcos 2018, an actual observationally-confirmed transient Earth minimoon, Jul 2006-Jul 2007) was
  independently reproduced this session by direct integration — an initial ~0.7-year close quasi-satellite
  episode reaching HALF Earth's Hill radius, then a transition to wider horseshoe libration with RECURRING
  approaches every ~4.6-9.2 years. That is a `quasi_cycler` by the catalogue's own definition, not a `cycler` —
  #523's search was strictly-periodic-orbit-shaped and could never have found it, however completely it swept
  (confirmed: 120 certified periodic orbits across the full intended Jacobi band, ALL 4.9-10.3x outside the
  Hill radius — #523 OUTSTANDING.md entry, 2026-07-03). Same structural point for #527's Sun-Jupiter Hilda case:
  Guido & Efthymiopoulos (arXiv:2604.00679) locate the actual heteroclinic/chaotic-transport structure on the
  UNSTABLE MANIFOLDS near the resonance separatrix, a genuinely different dynamical object than the periodic-
  orbit family core #527's DA/HOTM enumeration searched (and #530/#531 already showed the ONE genuinely-
  unstable sampled Hilda orbit's own homoclinic tangle stays at 3.48x the Hill radius — ruling out THAT specific
  narrow angle, not the broader chaotic-transport region near the separatrix itself).
  **NOT a straightforward reuse of #320's existing quasi_cycler machinery** — #320 targeted epoch-locked
  MULTI-MOON resonant tours (Saturn/Pluto/Earth-Moon systems, all V0-known-class per its own adjudication) via
  the SAME kind of construction as strict cyclers, just with a finite window; #535's targets are fundamentally
  APERIODIC/chaotic trajectories (a captured minimoon's drift-then-libration transition; separatrix-adjacent
  chaotic transport), which need: (1) broad seeding across the co-orbital/Hilda phase space (not just AT
  already-certified periodic orbits — the whole point is these are NOT periodic-orbit family members), (2) long
  non-strictly-periodic propagation, (3) a REPEATED-Hill-sphere-encounter detector as the admission criterion
  (count and space discrete close passes, don't require any closure residual), (4) the `quasi_cycler`
  epoch_locked/validity_window/n_returns schema fields populated from the ACTUAL encounter epochs found, not a
  periodic-orbit period.
  **Feasibility: MEDIUM, genuinely open on the admission-criterion question** — the physical phenomenon is real
  and sourced (2006 RH120's OWN behavior IS the positive control, already independently reproduced), but "how
  many Hill-sphere passes within what epoch window counts as admission" is a genuinely new criterion decision
  this project has not had to make for a genuinely chaotic (not quasi-periodic-torus, not periodic-orbit-family)
  trajectory before — settle this in writing BEFORE building the sweep, per the #339-style-criterion-trap
  discipline #523's original attempt already applied once.
  **CRITERION SETTLED 2026-07-03 (same session, before any sweep code):** see
  [[2026-07-03-535-quasi-cycler-transient-drift-admission-criterion]] for the full writeup. Summary: (1) encounter
  = inside the Hill sphere (reuses #523/#527's own threshold, no new distance invented); (2) one maximal
  continuous Hill-sphere residency = one return (not sub-counting periapsis wiggles inside a single episode); (3)
  returns must be separated by >=1 year outside the Hill sphere to count as distinct (Earth's own orbital period,
  chosen to filter numerical noise without presupposing the ~4.6-9.2 yr horseshoe return timescale this session
  already measured); (4) admissible iff SOME 10-15 year sub-window of the propagated trajectory contains 3-15
  distinct returns (the catalogue's own `quasi_cycler` schema bounds, `docs/notes/2026-06-16-catalogue-scope-
  taxonomy.md`), window bounds + actual epochs reported, never silently cherry-picked; (5) bounded-geometry check
  within that window (loosest return's closest-approach distance <=3x the window's tightest), a genuinely NEW
  numeric choice flagged as provisional pending real trajectory data, not carried over from an existing gate. The
  1-year separation floor (Earth-orbital-period-motivated) is explicitly flagged as NOT transferable to the
  Sun-Jupiter Hilda case without re-deriving from Jupiter's own orbital/libration timescale. Ready to build.
  **CORRECTION 2026-07-03 (same session, before any sweep ran): the "recurring approaches every ~4.6-9.2 years"
  claim above (repeated from #523's own docstring) does NOT hold up under careful re-verification.** Propagating
  #523's literal RH120 seed for 20/40/60 years (both `ydot0` sign branches) shows exactly ONE close approach
  (matching the confirmed ~0.7-year episode reaching ~half the Hill radius) followed by departure to ~2 AU from
  Earth with no return found within at least 60 years — not the claimed recurring horseshoe libration. Ruled out
  a seed-construction bug specifically: independently re-derived the SAME IC from the raw Keplerian elements
  (`a=0.998625, e=0.019833`) via proper vis-viva mechanics at perihelion (`vr=0` there by definition, confirming
  the original `xdot0=0` choice was correct, not a simplification artifact) — the re-derived `(x0,ydot0,C)`
  matches the original seed to 4 significant figures, so the seed itself is not the bug. The genuine physical
  conclusion: this specific object IS a single-transient encounter, consistent with the REAL 2006 RH120 being an
  observationally-confirmed ~1-year capture (Jul 2006-Jul 2007), not a multi-decade recurring companion — the
  "recurring approaches" claim was an error in an earlier interactive check that was never itself re-verified
  before being written down (`scripts/run_523_earth_coorbital_search.py`'s docstring corrected in the same
  commit as this note). **Consequence for #535**: the literal RH120 seed is a valid, honest NEGATIVE against its
  own admission criterion (a single return cannot satisfy `n_returns` in [3,15]) — it motivated the search but is
  not itself a candidate. Finding a genuine admissible `quasi_cycler` needs a broader seed search (varying
  `x0`/`C` away from this exact literal seed) for a nearby co-orbital initial condition that DOES show bounded,
  recurring returns.
  **GENUINE ADMISSIBLE CANDIDATE FOUND 2026-07-03 (same session): `scripts/run_535_earth_transient_quasi_cycler_
  search.py`, built + run end-to-end.** A 40-point `x0` scan at 2 Jacobi values (the corrected-target
  `C=2.9998797409719242` plus a broader-search `C=2.9839437412`) around the RH120 region, each point propagated
  50 years and passed through the #535 admission criterion (`search/hill_sphere_return_detector.py`), found
  EXACTLY ONE admissible point in the grid: **`x0=0.9920, xdot0=0.0, ydot0=0.13033911` (Jacobi
  `C=2.9839437412`)** — 4 total Hill-sphere returns over 50 years, with **3 of them (t=0, 2.99, 6.02 years,
  closest approaches 0.48-0.80x the Hill radius) falling inside a single qualifying 10-year window** (geometry
  ratio 1.65, well under the 3.0 gate). Verified robust across 3 independent integration tolerance/sampling
  settings (rtol 1e-10 to 1e-12, 800-5000 samples/year) agreeing to 4 decimal places before being trusted as a
  real dynamical feature rather than integration noise. The scan's own positive control (re-verifying this exact
  candidate) passed before the broader 80-point sweep ran; full ratchet green.
  **Honest scope of this finding**: this clears #535's OWN discovery-screen criterion — it is a genuine candidate
  for the `quasi_cycler` catalogue class, NOT yet a catalogue admission. Per the criterion note's own "what this
  does NOT decide" section: no `dv_band`/flyby-quality classification has been attempted, and this single grid
  point was found in an ad hoc 40-point 1D scan around the motivating seed, not a systematic 2D `(x0,C)` search —
  the true extent of the admissible region (is this an isolated point or part of a narrow continuous family?) is
  unknown. The standard V0-V5 catalogue gauntlet is the next step if this candidate is pursued toward admission,
  not run here.
  **WIDTH CHARACTERIZED 2026-07-03 (same session, before any gauntlet run): the admissible region is a genuine,
  reproducible but EXTREMELY NARROW corridor, not a broad feature.** A fine `x0` scan (41 points, step 0.0001)
  at the fixed admissible `C` found admission at ONLY `x0=0.9920` and `x0=0.9921` — a ~1e-4 AU (~15,000 km) wide
  window out of the 0.004 AU range scanned. A `C` scan (41 points, step 0.0002) at the fixed admissible
  `x0=0.9920` found ZERO admissible points anywhere in a 0.008-wide range straddling the true admissible `C` on
  both sides (the two nearest grid points, `C=2.9838` and `C=2.9840`, both non-admissible) — the `C`-direction
  width is narrower than this 0.0002 step, i.e. the admissible set behaves like a thin FILAMENT in the
  `(x0,C)`-plane rather than an open 2D patch, a signature of chaotic phase-space structure (consistent with de
  la Fuente Marcos & de la Fuente Marcos 2018's own Arjuna-population statistic: "~8% of Arjunas have v_rel^2 <
  v_esc^2 at perigee" -- i.e. MOST nearby phase space genuinely does not lead to a capture episode; narrow
  corridors are the physically expected structure here, not a search-code artifact). Re-checked the `x0` boundary
  at tighter tolerance (`rtol=1e-13` vs the standard `1e-11`) and got the IDENTICAL admissible/non-admissible
  pattern at 5 test points spanning the boundary — the narrowness is a real, reproducible dynamical feature, not
  integration noise.
  **Implication**: this is a scientifically valid finding (a real transient-drift `quasi_cycler` corridor exists
  and was found, not fabricated or a numerical fluke), but its practical/mission relevance is limited by the
  corridor's own narrowness, and pushing this SPECIFIC point through the full V0-V5 gauntlet (particularly the
  `quasi_cycler`-scoped V2 bounded-drift-over-the-validity-window gate) risks a fragile pass that a tiny
  real-ephemeris/perturbation difference could break — the gauntlet should explicitly stress-test sensitivity to
  small IC perturbations given this characterization, not just check nominal-IC pass/fail. Awaiting a decision on
  whether to (a) proceed to the V0-V5 gauntlet with this sensitivity caveat carried through, (b) search more
  broadly across the co-orbital region for a WIDER, more robust corridor before investing further, or (c) treat
  this as a complete, valuable discovery-catalogue entry as-is (a genuine but narrow corridor) without further
  gauntlet investment.
  **WIDER-SEARCH ATTEMPT ABANDONED 2026-07-03 (same session): a 3111-point `(x0,C)` grid scan
  (`x0` in [0.94,1.06], `C` in [2.95,3.05], 61x51) launched to look for a WIDER corridor produced ZERO output
  after ~48 minutes (a real instrumentation gap: unbuffered-stdout blindness on a `nohup`-launched process, not a
  computation failure — `/proc/<pid>/status` confirmed it was genuinely CPU-bound throughout) and was killed
  before completing. This is an ABANDONED/inconclusive attempt, not a clean negative — it does not rule out a
  wider corridor existing elsewhere in that grid.** A parallel investigation into 3753 Cruithne's real orbital
  elements (`a=0.998 AU, e=0.515, i=19.8°`) as an alternative wide/stable-companion seed found ZERO Hill-sphere
  returns at all: Cruithne's large eccentricity is PRECISELY what keeps it out of Earth's Hill sphere, i.e. its
  long-term dynamical stability and #535's repeated-close-encounter admission criterion are in structural
  tension — searching near a KNOWN STABLE companion was the wrong strategy; genuinely wide, robust admissible
  corridors may be structurally rare in this co-orbital regime (consistent with the ~8% Arjuna capture-probability
  statistic already cited above).
  **SENSITIVITY PRE-CHECK RUN 2026-07-03 (same session), resolving the awaiting-decision above in favour of (c)
  with an explicit caveat, per the user's own request for "a cheap pre-check" before any V0-V5 investment:
  `scripts/run_535_er3bp_sensitivity_check.py`, built + run.** Re-propagated the EXACT admissible IC
  (`x0=0.9920, xdot0=0.0, ydot0=0.13033911`) through the project's existing ER3BP core (`core/er3bp.py`,
  pulsating-rotating frame, independent variable = true anomaly) at Earth's REAL orbital eccentricity
  (`e=0.0167`, NASA Earth fact sheet) instead of the idealized `e=0` CR3BP model — a genuine, minimal,
  physically-motivated perturbation (not an arbitrary IC jitter). Positive control passed first (the ER3BP code
  path at `e=0` exactly reproduces the CR3BP result: 4 returns, 2 admissible windows) before trusting the
  perturbed result, per this project's own verify-before-trust discipline. **VERDICT: COLLAPSES.** At
  `e=0.0167` only the initial pass remains (1 return total, 0 admissible windows) — the trajectory departs and
  does not return within the full 50-year propagation. This is not a shrinkage, it is a total collapse, and it
  independently confirms the width-characterization finding above (a filament this narrow in the idealized
  model's own `(x0,C)` freedoms is exactly the signature of a feature that does not survive real dynamics) —
  the same pattern this project has hit repeatedly (S1L1, #388, the #480 EGGIE construction): an idealized-model
  closure on a knife-edge, not a robust dynamical feature. **Decision: do NOT push this candidate into the
  V0-V5 gauntlet** — it would fail the first real-ephemeris-adjacent gate for a now well-understood reason, not
  a novel one. The idealized CR3BP finding itself remains valid and stays on record exactly as scoped (a genuine
  existence result within the idealized circular-restricted model, characterized as fragile from the outset).
  #535 is complete as a discovery-screen exercise; no further investment planned unless a future broader/smarter
  search (the abandoned wide scan above, redone with proper `python3 -u`/flush=True instrumentation and likely a
  coarser first pass) turns up a structurally wider corridor.

- **#557 ✓ RUN (2026-07-11) — CLEAN, WELL-CHARACTERIZED NEGATIVE** (header corrected 2026-07-15: this
  bullet's OPENING line said "planning phase first — do NOT build/run yet" for months after the
  actual resolution was appended ~96 lines further down in this SAME bullet — see the "✓ RUN" text
  below. This misled a same-day OUTSTANDING.md dashboard audit AND a manual review of it into
  re-presenting an already-answered scope question to the user. Read a bullet to its END before
  concluding a task's status; a stale HEADER on an otherwise-complete bullet is a real pattern, not
  just contradicting duplicate bullets — see `[[feedback_outstanding_current_state_maintenance]]`.)
  extend #535's validated Hill-sphere-return `quasi_cycler` screen
  (`search/hill_sphere_return_detector.py`, admission criterion already settled in writing for
  Earth) to the Sun-Jupiter quasi-Hilda population. Motivation: unlike Earth's RH120 (a
  documented ONE-SHOT transient, Jul 2006-Jul 2007, no known recurrence — and #535's own
  admissible Earth corridor turned out to be a ~15,000 km knife-edge that COLLAPSED entirely
  under real orbital eccentricity, ER3BP e=0.0167 vs the idealized e=0), Jupiter's quasi-Hilda
  population includes REAL comets with documented, REPEATING temporary satellite captures (e.g.
  P/Gehrels 3, 111P/Helin-Roman-Crockett) — a recurring phenomenon in the literature, not a
  hopeful search for something that might not exist. Checked the empty-region registry first:
  neither existing Jupiter-Hilda entry (`sun-jupiter-hilda-32-mmr-dahotm`, a periodic-orbit-family
  DA/HOTM enumeration; `hilda-c3.14-homoclinic-connection-hill-encounter`, a single-homoclinic-
  tangle check) is this method — a broad Hill-sphere-return sweep across quasi-Hilda phase space
  is genuinely unscreened territory, not already-registered-empty. Matches the pattern behind
  both of this project's confirmed novel hits (#312 Uranus, PC (3,2) via #494/#504): validated
  tool pointed at unscreened territory, not a new capability build.
  **Planning-phase scope (this dispatch):** per #535's own `#339-style-criterion-trap` discipline
  (settle the admission criterion in writing BEFORE any sweep code, exactly as #535 itself did
  for Earth before building anything), produce a concrete, reviewable plan covering: (1) re-derive
  the return-separation timescale from Jupiter's own orbital/libration period (#535's own
  criterion note explicitly flags its 1-year, Earth-orbital-period-motivated floor as NOT
  transferable to Jupiter without this re-derivation); (2) decide the seeding strategy — literal
  ICs from a documented real object's actual orbital elements (P/Gehrels 3 or 111P/Helin-Roman-
  Crockett, sourced not assumed) as a positive-control-style anchor, a broad quasi-Hilda
  phase-space scan, or both; (3) adapt the bounded-geometry/window-count admission thresholds
  (currently Earth-scale numbers) to Jupiter's Hill-sphere scale; (4) identify what, if anything,
  needs building versus what #535's existing detector/criterion machinery already handles
  verbatim; (5) an honest cost/risk estimate given #535's own history (a real result was found
  for Earth, but it was fragile and the broader search attempt was abandoned to an
  instrumentation gap, not a clean negative — flag what would need to be different this time to
  avoid repeating that). Do NOT write the sweep script or run anything yet — this is a design
  document for review, mirroring how #538's plan was written and reviewed before any Task 1
  code was dispatched.
  **Recommended model:** Opus for this planning pass (the #339-style-criterion-definition
  judgment call + Jupiter-timescale re-derivation is exactly the trust-bearing numerical-methods
  design work this project's model-tiering policy reserves for Opus, not a mechanical task).
  Implementation model to be decided after the plan is reviewed, per explicit user instruction —
  do not auto-proceed to a build dispatch on plan completion.
  **PLAN WRITTEN 2026-07-11 (planning pass complete, awaiting user review — nothing built or run,
  no catalogue/empty-region edits):** `docs/superpowers/plans/2026-07-11-557-jupiter-quasi-hilda-
  transient-capture-plan.md`. Headline finding that gates everything else: the catalogue
  `quasi_cycler` validity-window (10-15 yr) is Earth-timescale-calibrated and is STRUCTURALLY too
  short for Jupiter's natural return spacing — the re-derived separation floor is Jupiter's own
  orbital period (11.86 yr, the strict analogue of #535's Earth 1-yr noise floor; the Hilda
  libration period ~250-300 yr was rejected as presupposing the resonance), and ≥3 distinct returns
  each separated by ≥11.86 yr need ≥~24 yr, exceeding the entire un-rescaled window. So the literal
  window makes the screen a near-certain STRUCTURAL empty (a units artifact, not a dynamical
  negative — must NOT be registered as an empty region). The physically-consistent fix rescales the
  window by T_J/T_E≈11.86 to ~120-180 yr (count 3-15 and the 3x geometry ratio are scale-invariant
  and transfer as-is; `r_hill`=(μ/3)^(1/3)=0.0683 nondim/0.355 AU recomputes automatically) — but
  that departs from the schema's fixed mission-relevance window and is a `quasi_cycler`-class-scope
  decision for the user, parallel to the #320/#535 scoping discussions, to be made BEFORE any build.
  Other recommendations: seed via option (c) (real-object anchor — P/Gehrels 3 elements sourced from
  JPL SBDB, Oterma/Koon-2001 as the in-corpus fallback; NOT the #527/#530 C=3.14 periodic seed,
  which sits above C_L1=3.039 with necks closed and produces zero Hill encounters by construction —
  the temporary-capture regime needs the neck-OPEN band C∈~[3.00,3.038]), then a coarse-first
  broad scan, mirroring #535's RH120-first structure; ~90% reuse (the detector is fully
  system-agnostic already, the run_535 sweep + ER3BP-sensitivity scripts are re-parameterizable
  clones); key risk-mitigation vs #535's own history — run the real-eccentricity ER3BP check
  (Jupiter e=0.0489) EARLY on the anchor to gate the broad scan (Koon 2001 notes e "plays little
  role during the fast resonance transition", so the fast tube-mediated capture may be less fragile
  than Earth's slow horseshoe — a hypothesis to test, not assume), and `python3 -u`/flush +
  incremental checkpoint runlog + coarse-first grid from the first launch (the #535 wide scan died
  to an unbuffered-stdout instrumentation gap, not a computation failure). Confirmed by reading
  `data/empty_regions.jsonl` directly that none of the three existing Sun-Jupiter Hilda entries use
  the Hill-sphere-return method AND all sit at the closed-neck C=3.14 energy — genuinely unscreened.
  **FABLE REVIEW 2026-07-11 (window-rescaling implications + CPU-cost question, answered):**
  corrections applied directly to the plan doc and to `docs/notes/2026-06-16-catalogue-scope-
  taxonomy.md`. Key findings: (1) the "10-15 yr" figure was ALREADY dead letter — not enforced
  anywhere in code/schema, and the catalogue's own one `quasi_cycler` row (#312) already has an
  83-year `validity_window`, 5.5x the stated ceiling, admitted through the full gauntlet — so
  generalizing this is a DOCUMENTATION fix (system-period-relative floor/window, not a fixed-years
  number), not a schema migration; done directly in the taxonomy doc rather than deferred. (2) The
  right formulation is dimensionless — floor = 1 rotating-frame period of the CR3BP system under
  study, window = 10-15 such periods — which correctly generalizes to MOON-SYSTEM quasi_cyclers too
  (like #312 itself: the relevant period is the moon-pair's ~6-day synodic period, NOT Uranus's
  84-year heliocentric period; the plan's literal "primary's orbital period" phrasing would have
  been wrong by ~5000x for that case). (3) Catalogue-wide blast radius is small (other classes/M7-
  novelty-matching unaffected) but real: V2/V4 validation budgets scale in real years and can exceed
  kernel coverage at Uranus/Neptune-heliocentric scale (DE440 ~1100 yr span; add a
  `min(validity_window, kernel_coverage)` truncation rule, now in the taxonomy doc). (4) **The
  CPU-cost question is answered cleanly: positive, and NOT more expensive.** CR3BP integration cost
  scales with nondimensional rotating-frame revolutions, not calendar years — a Jupiter screen under
  the corrected window integrates the same ~10-60 revolutions as #535's Earth screen, i.e. the SAME
  per-point cost, not "8-12x more" as this plan originally (wrongly) estimated in its own §5. That
  wrong estimate was traced to a REAL latent bug the review caught: #535's own script hardcodes
  `2*pi rad = 1 yr` (true only for Sun-Earth); a naive Jupiter clone that kept this constant would
  actually integrate ~500 Jupiter-periods instead of ~50 — reproducing the original cost estimate as
  a genuine ~10x overrun, not a pessimistic guess. **Both the plan doc (SS4, SS5) and the taxonomy
  doc have been corrected accordingly** — the time-unit conversion must go through Jupiter's own
  period explicitly, not #535's literal constant. Net recommendation (adopted): proceed under the
  corrected dimensionless window, do the taxonomy-doc fix now (done), do not touch
  `catalogue.schema.json`/`validate.py` (nothing is enforced there yet, one class member exists).
  Implementation still awaits a separate go/no-go — this correction pass touched only the
  plan/taxonomy docs, no sweep code, no catalogue/empty-region edits.
  **✓ RUN (2026-07-11, commit `976eb75` + this commit) — CLEAN, WELL-CHARACTERIZED NEGATIVE.**
  Criterion note first, per the #339-style discipline:
  `docs/notes/2026-07-11-557-jupiter-quasi-hilda-admission-criterion.md` (floor = 1 Jupiter
  period = 11.868 yr, window = 10-15 periods = 119-178 yr, both re-verified in-repo, not
  recalled — matches the corrected taxonomy doc). **Positive control PASSED**: 82P/Gehrels 3's
  real JPL SBDB osculating elements (orbit_id 19, DE431, retrieved 2026-07-11) converted to a
  rotating-frame IC via the proven #535/#523 vis-viva-at-perihelion path give Jacobi C=3.02943,
  matching the INDEPENDENTLY-sourced SBDB Tisserand parameter (3.027) to 0.002 — a genuine
  cross-check, not circular — and the constructed anchor enters Jupiter's Hill sphere to 0.12
  R_hill (reproducing Gehrels 3's documented temporary capture). **Coarse CR3BP scan** (264
  points, x0 in [0.60,0.92] x 8 Jacobi values in [3.000,3.035], all below C_L1=3.0388 so the
  L1/L2 necks are open): **16/264 (6%) idealized-CR3BP-admissible** — a materially richer band
  than #535's single ~15,000-km-wide Earth corridor, consistent with Koon 2001's fast tube-
  mediated capture mechanism producing wider structure than Earth's slow horseshoe. **The
  ER3BP real-eccentricity gate (Jupiter e=0.04838624) is decisive: 15/16 collapse outright**
  (lose the recurrent 3+ return structure entirely), and **the 1 apparent survivor (x0=0.86,
  C=3.02) is an explicitly-tested KNIFE-EDGE COINCIDENCE** — admissible ONLY at exactly
  e=0.04838624 (collapses at e=0.0489/0.040/0.025/0.061), ONLY at exactly x0=0.860 (collapses
  at +/-0.005), ONLY at exactly C=3.020 (collapses at +/-0.005) — the identical signature to
  #535's Earth corridor and this project's recurring idealized-closure-on-a-knife-edge pattern
  (S1L1, #388, #480 EGGIE). **Koon 2001's eccentricity-insensitivity hypothesis does NOT hold
  for the SUSTAINED multi-return quasi_cycler criterion**: individual fast captures still occur
  under real e (the anchor itself still captures), but the recurrent structure a `quasi_cycler`
  needs is destroyed by real eccentricity, same as at Earth. **Registered in
  `data/empty_regions.jsonl`** (`sun-jupiter-quasi-hilda-transient-capture-quasi-cycler`) under
  the explicit `criterion_version: "557-optionA-system-period-relative-2026-07-11"` tag, per
  #557's own scoping — never to be confused with a fixed-years or Option-B structural-empty
  result. **Scope/resweep note** (recorded in the registry entry): this is a coarse, planar,
  single-Jupiter-phase screen; the fragility is intrinsic to the mechanism (not a grid-
  resolution artifact) — a genuinely different method (full nonlinear manifold globalization,
  a 3D/inclined model, or a real-ephemeris lane) would be needed to test for robust structure
  this screen cannot, not merely a finer CR3BP grid. Total wall-clock: 169.5s.

- **#536** — ✓ Resolved (2026-07-08) Apply the Linking-Number/QP-Torus Heteroclinic Screening Tool to a Genuinely Unmapped System
  Resolution: Implemented the Jupiter-Europa L1/L2 matched-Jacobi torus connection search in `scripts/run_536_jupiter_europa_connection.py`, finding 0 connections at C=3.001500 (commit 8b6c60f).
  **⚠ STAMP (2026-07-10, #548): this 0-connection result is METHOD-INVALID — DO-NOT-CERTIFY, not an
  empty region.** The qp_tori/qp_torus_heteroclinic linking-number pipeline was SHELVED by #548
  (postmortem `docs/notes/2026-07-10-postmortem-548-linking-number-pipeline.md`): a reframed
  Owen & Baresi positive-control sweep found zero linking-number sign changes and never produced a
  nonzero torus-level linking number even where connections should exist. This Jupiter-Europa 0/N
  is uninterpretable output from an unvalidated screen; do not register it as an empty region.
  **[Ordering note added 2026-07-15, full-file audit: everything below this line is the ORIGINAL
  2026-07-03 proposal text that predates the "✓ Resolved"/STAMP block above — it describes what
  was PROPOSED, not a second still-open item. Ends "Not yet built" because that was true when
  written, before #536 executed and was later stamped invalid. Kept in place for history, not a
  live TODO.]**
  (Jovian-Moon Tori) (allocated 2026-07-03, same session, from an independent second-opinion review of the
  #532/#534 discovery plan). **Motivation**: #534's Earth-Moon L1<->L2 target is not actually a novelty attempt
  as posed — quasi-halo/torus heteroclinic connections in that specific system are extensively published
  (Gómez et al., Koon-Lo-Marsden-Ross, Haapala & Howell, Olikara), so a successful #534 run there demonstrates
  the #522 linking-number pipeline works, it does not discover anything new. The catalogue's own history backs
  this up: this project's one confirmed genuinely novel find (the Uranus-system cycler, #312) came from an
  under-searched system, not a re-application of existing machinery to the most heavily-published system
  available. **Scope**: once #534's transit/non-transit branch-classification blocker is resolved and the
  pipeline is validated end-to-end against the known Earth-Moon connection (treat that as #534's OWN positive
  control, not a separate deliverable), re-target `genome/qp_torus_heteroclinic.py`'s `scan_linking_number` at
  quasi-periodic tori in a Jovian-moon system (e.g. Europa-Ganymede or Io-Europa near-resonant pairs, reusing
  `search/moon_cycler_genome.py`'s existing Jovian-system definitions) — a system this project has not yet
  screened for whiskered-torus heteroclinic structure at all. **Feasibility: LOW additional build cost once
  #534 is unblocked** — same pipeline, new target system and halo/torus seed pair; the real work is finding a
  genuinely isoenergetic Jovian-moon torus pair the way #534 found one for Earth-Moon (matched-Jacobi
  continuation), not new machinery. **Explicitly gated on #534's transit/non-transit classification landing
  first** — do not fork effort into a second system while the first one's own manifold-branch selection is
  still unresolved. Not yet built.
- **#533** — ✓ Resolved (2026-07-08). Genuine Coherent QBCP Model: Implemented the QBCP equations of motion and STM in canonical variables using Gimeno-Jorba 2018 Table 4 coefficients, verified against circular BCR4BP limits and finite differences (commit to follow).
  **Correctness addendum (2026-07-14, #592/#593):** the implemented EOM had a real bug (alpha_6
  incorrectly scaled the whole Newtonian potential, not just the Sun term per Rosales-Jorba 2023
  Eq. 3) — fixed under #592. This build-time verification method (FD-against-the-same-formula)
  could not have caught it by construction; the model itself is still the right deliverable and is
  now more correct, but any specific past numeric result produced with it (see #537/#544 below)
  predates the fix and should not be trusted as-is. `docs/notes/2026-07-14-qbcp-alpha6-scaling-fix.md`.
- ✓ Resolved (2026-07-09) **#537** — QBCP Cross-System Torus Connection Search: Completed the cross-system invariant torus connection search from Sun-Earth L2 Torus to Earth-Moon L2 Torus in the QBCP model, finding a highly accurate refined connection candidate with a 12,034 km position gap and a 911 m/s velocity gap at crossing time 3.66 TU (scripts/run_533_qbcp_connection.py). **Caveat found 2026-07-09 while scoping #538**: the refining `least_squares` solve is a **3-equation/4-unknown rank-deficient system that never includes velocity in its residual** — the 12,034 km / 911 m/s numbers are post-hoc diagnostics from a separate propagation, not quantities the optimizer drove toward zero. Treat as a basin seed for #538's properly-posed solve, not a near-converged boundary condition. Plan: `docs/superpowers/plans/2026-07-09-538-qbcp-cross-system-cycler-correction.md`.
  **⚠ STAMP (2026-07-15, full-file audit — consistency fix): #537's headline "highly accurate...
  connection" carries the SAME functional status as #534/#536's explicit STAMPs above (both
  DO-NOT-CERTIFY), it was just never given the formal marker. Two independent, compounding
  reasons, both already documented below in this same entry: (1) the rank-deficient-solve caveat
  immediately above (numbers were never actually driven to zero); (2) the alpha_6 staleness note
  further down (computed under pre-#592 buggy dynamics, doesn't reproduce as-is). DO NOT CITE
  #537's 12,034 km / 911 m/s figures as a validated result.**
  **Second, independent staleness reason (2026-07-14, #593, reconciled):** this ran under the
  pre-#592 buggy `qbcp.py`. #593's full investigation (including reconciling with #544's own prior,
  reverted attempt at this exact fix) confirmed #592 is a genuine, net-positive correctness fix (the
  SE-L2 regression #544 originally reported does not reproduce; the POL1-agreement improvement is
  independently confirmed) — so #592 stays applied. That still means #537's specific numbers were
  computed under different (buggy) dynamics than current `main` and don't reproduce as-is; not a
  reason to distrust the fix, just a reason to re-derive fresh rather than patch old numbers. No
  re-run warranted now (not currently active work); whoever next resumes #538/#539/#540 should
  re-derive fresh under current `main`. `docs/notes/2026-07-14-593-qbcp-alpha6-impact-scoping.md`.
- **#538** — QBCP Cross-System Periodic Orbit (Cycler) Correction: Chain the forward and reverse QBCP torus connections refined in #537, set up a multi-segment boundary value problem, and run a differential corrector to converge on a mathematically exact periodic orbit in the time-periodic QBCP model. Verify its stability and transport utility (close approaches to Earth/Moon) to evaluate it as a potential cislunar cycler. Full task breakdown + recommended model per task in `docs/superpowers/plans/2026-07-09-538-qbcp-cross-system-cycler-correction.md`.
  **Task 0 done (2026-07-09):** the preflight-gate blocking finding (all four run_522/533/534/536 scripts now call `preflight_search()`) landed in commit `ede2d0c`.
  **Tasks 1-3 built and run (2026-07-09):** `scripts/run_538_qbcp_cycler.py` implements the well-posed 12-unknown/18-residual multi-segment corrector the plan called for (position+velocity+time+phase-return at both crossings, not #537's underdetermined 3-eq/4-unknown version). Multi-start ballistic (base seed x{lm,trf} + perturbed restarts) and the Task 3 powered fallback both ran to completion; neither reached the 1e-8 target (best ballistic norm 1.068, best powered norm 0.866 with an implied 509.6 m/s leg-1 impulsive dv).
  **This negative is NOT meaningful evidence either way — do not cite it as "no QBCP cislunar cycler exists".** Diagnosis found the corrector was closing onto a **non-invariant EM-L2 "torus"**: `torus_em.invariance_residual = 3.439` vs the SE-L2 torus's clean `3.093e-05` (same `correct_qbcp_torus` function, ~5 orders of magnitude worse). An isolated diagnostic (escalating `max_iter` 30→100→300 for `correct_qbcp_torus`'s EM-side `least_squares` call) landed on the **identical** residual (1.1039) and nfev (192) at every budget — the LM optimizer self-terminates via `xtol`/`ftol` well short of any iteration cap, so this is a genuine plateau, not a starved budget. The root cause traces earlier still: the **BCR4BP-side mu_sun continuation itself doesn't cleanly re-converge at each step** (`bcr4bp_torus_residual` residual goes 0.139→0.137→0.111→0.208→0.337 as mu_sun ramps 0→full, ending worse than it started), handing the QBCP corrector an already-degraded seed (1.69) before it even begins.
  **This is not new to #538 — it's inherited from #533/#537.** `scripts/run_533_qbcp_connection.py` uses the identical `tol=1.5e-1` for the EM QBCP torus and *does* print `Earth-Moon Torus Residual: {...}`, but that print was apparently never read when #537 was marked resolved. #537's "highly accurate refined connection" (12,034 km / 911 m/s gaps) may have been built on this same non-invariant EM-side surrogate the whole time — **#537's resolution should be treated as provisional pending a re-check once the EM torus construction is fixed**, not retracted outright (its position/time-only refinement never claimed velocity closure anyway, per #538's own Context section).
  **#538 stays open, blocked on #544** (see below), not resolved as a negative. Tasks 4 (independent cross-check) and 5 (writeback) are moot until there's a genuinely invariant EM-L2 torus to close onto.

## Fix the EM-L2 QBCP torus construction

- **#544** (P0, blocks #538/#539/#540) — Found 2026-07-09 while running #538:
`genome/qbcp_torus.correct_qbcp_torus` reliably
converges the Sun-Earth L2 torus (residual 3.093e-05) but never converges the Earth-Moon
L2 torus (residual 3.439, ~5 orders of magnitude worse) — a genuine plateau, not a
budget problem (max_iter 30/100/300 all land on the identical residual 1.1039 and nfev
192; the LM optimizer self-terminates via xtol/ftol well short of any cap). Root cause
traces earlier: the BCR4BP-side mu_sun-continuation seed (`bcr4bp_torus_residual` +
`least_squares(method="lm")`, one shot per step, no per-step convergence gate) does not
cleanly re-converge as mu_sun ramps 0→full — residual goes 0.139→0.137→0.111→0.208→0.337,
ending *worse* than it started — handing the QBCP corrector an already-degraded seed
(1.69) before it does anything. `scripts/run_533_qbcp_connection.py` (#537) has the
identical code path/tolerance and prints the same residual, unread. **This blocks #538's
closure attempt from meaning anything (its "honest negative" is contaminated, not a real
physics result) and would identically block #539/#540's reuse of the same corrector
methodology on Jupiter-Europa/Uranus.**
**Feasibility: MEDIUM** — likely candidates, roughly in order of cheapness: (a) add a
per-step convergence check to the mu_sun continuation (re-solve to a real tolerance at
each step, not one `lm` call and move on — the classic naive-vs-pseudo-arclength
continuation gap #524 already built deflated Newton/continuation machinery to fix
elsewhere in this project); (b) finer mu_sun steps (5 → 20+); (c) more Fourier modes
(n_modes=2 may be insufficient for the EM-L2 torus's true QBCP quasi-periodic content,
even though it's sufficient for SE-L2); (d) reconsider the EM-L2 target Jacobi (3.17) —
check whether a torus genuinely exists there at this energy before assuming the
corrector is at fault.
**Recommended models:** diagnosing which of (a)-(d) is the actual fix → **Opus** (trust-
bearing numerical-methods judgment, same class of work as #538's Task 1). Implementing
the fix once diagnosed and re-validating against SE-L2's already-clean convergence as a
positive control → **Sonnet**.

**PARTIALLY RESOLVED (2026-07-10).** The actual bug was none of (a)-(c): it was a
**singular phase-pin gauge**. `se_lyapunov_to_bcr4bp_torus_seed` / `se_lyapunov_to_qbcp_torus_seed`
(and the four scripts' inline EM-side duplicates of the same pattern) picked
`phase_pin_idx = argmax|Im(c_1)|`, but the residual pins `Im(c_1[phase_pin_idx]) = 0`,
whose gauge derivative w.r.t. a circle rotation is `Re(c_1[phase_pin_idx])` — so the pin
coordinate needs a *large real part*, not a large imaginary part. For a symmetric
Lyapunov orbit c_1 splits cleanly real/imaginary by coordinate (x/vy real, y/vx
imaginary); `argmax|Im|` picks a coordinate whose `Re ≈ 3e-12`, an essentially singular
gauge the corrector cannot use to fix the rotational phase, so it stalls. Confirmed
against `qp_tori.py:655`'s `correct_qp_torus`, which correctly uses `argmax|Re|` — the
EM-side helpers/inline scripts had silently diverged from that pattern. Fixed identically
in `genome/bcr4bp_torus.py`, `genome/qbcp_torus.py`, and the four scripts' inline
duplicates (`run_522`/`run_533`/`run_538`/`search_coherent_connections.py`).
**Result:** SE-L2 stays clean (3.093e-05 → 2.537e-05, slightly better). EM-L2 improves
~7x (3.439 → 4.771e-01) but does **not** reach the 1e-3 target. The remaining gap is
candidate (d), not a corrector bug: the EM-L2 Lyapunov at Jacobi 3.17 sits right at
C_L2≈3.172 and swings deep toward the Moon (closure state x=0.967 vs. Moon at x=0.988),
giving genuinely high harmonic content (mode-ratio ≈0.87, vs. SE's benign 0.35) that
`n_modes=2` cannot resolve — a Jacobi scan confirmed the whole EM-L2 family near C_L2 is
high-harmonic; reaching 1e-3 there would need ~40 modes (infeasible for #538/#539/#540's
downstream correctors) or an independent off-grid cross-check before trusting a
masked/oversampled residual (the fix deliberately did not reach for that — a masked low
residual on an unresolvable torus would be exactly the "it closed!" trap this project's
orbit-closure discipline warns against). **Open decision for whoever picks up #538 next:
reconsider the EM-L2 target orbit** (a smaller-amplitude object away from the Moon, or an
EM-L2 halo, where `n_modes=2-4` should genuinely resolve it as it does for SE) rather than
pushing more modes at this specific high-harmonic orbit. Validation: ruff/mypy clean,
`tests/genome/test_bcr4bp_torus.py` + `tests/genome/test_qbcp_torus.py` +
`tests/core/test_qbcp.py` + `tests/core/test_bcr4bp.py` all pass (18/18). #538's own
multi-hour closure attempt was deliberately NOT re-run as part of this fix (that decision
belongs to whoever next works #538, given the EM-L2 target-orbit question above is still open).

**FOLLOW-UP (2026-07-10) — the bug is deeper than a target-orbit choice.** Found a much
better EM-L2 target (Jacobi 3.13 vs 3.17): the BCR4BP-side mu_sun continuation converges
beautifully (6.02e-6), but handing that exact converged seed to `correct_qbcp_torus`
makes it *worse* (2.18e-2, ~3600x), not better — the opposite of what a small QBCP
refinement of a converged BCR4BP torus should do. A dispatched read-only code-review
agent found `qbcp_torus_residual`/`bcr4bp_torus_residual` and the PV/PM handoff
(`propagate_qbcp_pv` correctly calls `state_pv_to_pm`) are both internally consistent —
ruled out. It flagged a real, empirically-confirmed discrepancy: `bcr4bp._sun_position`
places the Sun at angle 0 rotating CCW (prograde); `qbcp.evaluate_alphas`'s alpha_7/alpha_8
place it at angle ~180° rotating CW (retrograde) — confirmed by direct evaluation
(`t=0`: BCR4BP sun=(388.81,0), QBCP sun=(-392.09,0); the sign flips again by `t=0.1`).
Testing all 4 x/py mirror-sign combinations on the seed did **not** fix the residual
(best 2.75, still O(1)), so it isn't a simple axis-relabeling either.

**A positive-control test with Rosales/Jorba (2023) Table 4's own published QBCP
periodic-orbit substitutes (POL1/POL2, the "dynamical substitutes" for EM-L1/L2) is
more damning: propagating either published golden ORBIT for one Sun-synodic period
`T_s` under our own `qbcp_eom` does not close at all** (residual 1.9-4.0 depending on
sign convention, tested all 4 combinations of the x/py signs — the |x| magnitudes
0.8369/1.1557 match our own EM-L1/L2 distances suspiciously exactly, so a convention
difference is likely, but no sign combination tried actually closes). Per this
project's "verify a filter with a positive control before trusting a negative"
discipline (`[[feedback_verify_gauntlet_with_positive_control]]`), **this means the
`correct_qbcp_torus` machinery has never been validated against a real, sourced QBCP
solution — the EM-L2 non-convergence could be a genuine bug in `evaluate_alphas`'s
Fourier table transcription (Gimeno/Jorba 2018 Table 4, order k≤13) or in the model's
epoch/period convention, not (only) a target-orbit or gauge issue.** This is now a
higher-priority, more fundamental thread than the target-orbit search — the target
orbit could be arbitrarily well-chosen and the corrector would still fail if the model
itself doesn't close on its own published golden. **Next step (not yet done): resolve
the POL1/POL2 closure discrepancy first** (check period assumption, `theta_sun0`
epoch convention, and the alpha_i table transcription against Gimeno/Jorba 2018 Table 4
directly) **before further EM-L2 target-orbit or corrector-gauge work.**

**FOLLOW-UP-2 (2026-07-10) — RESOLVED as a metric/method artifact: the QBCP model is
structurally correct; POL1/POL2 "non-closure" is NOT a model bug.** Worked the POL1/POL2
discrepancy to the bottom (all diagnostics in the session scratchpad). Findings, in order:
- **The forward-propagation closure test was the wrong metric.** POL1/POL2 are dynamical
  substitutes of the *violently unstable* EM collinear points. The frozen-time
  linearization of `qbcp_eom` at the L1/L2 x-locations gives an unstable rate of
  **2.979 (L1) / 2.199 (L2)**, matching the CR3BP collinear rate (2.932 / 2.159) to ~1.6%
  — i.e. the Sun term is the expected small O(eps^2) perturbation and the model's L-point
  *stiffness is correct*. That rate implies a one-period (`T_s ≈ 6.79`) monodromy
  multiplier of `exp(rate·T_s) ~ 1e6–1e8`, so a single forward propagation of even a
  perfect IC amplifies any roundoff / model-instance offset to O(1). The reported
  residuals (1.9–4.0) are that amplification, not a defect. (The earlier "clean
  saddle×center×center monodromy, mult 216" reading was a red herring: det=1 and
  reciprocal-pair eigenvalues are automatic for *any* symplectic STM, and the 216 was the
  STM of a trajectory that had already flown off L1 along the unstable manifold, not a
  Floquet multiplier of a periodic orbit.)
- **Constructed the genuine substitute properly** by continuation from the *exact* CR3BP
  EM-L1 equilibrium (`s=0`, a fixed point of the stroboscopic map to machine precision) to
  the full QBCP (`s=1`) with a **multiple-shooting** corrector (48 segments over `T_s`).
  It converges cleanly at every step (periodicity residual **9.1e-12** at `s=1`) to a real
  QBCP periodic orbit. That orbit's closest approach to Rosales/Jorba's published POL1
  (over a full period, matching x and py to ~1.7e-3) is **2.1e-2** — i.e. our model's own
  L1 substitute agrees with the published golden to ~2%, consistent with a
  Gimeno-2018-alpha-table vs Rosales-2023-POL model-instance (Fourier refit) difference
  rather than a code error. No `theta_sun0` epoch offset (swept 0→2π) and no x/py sign
  combination reconciles the remaining gap; it is a genuine small difference between our
  (Gimeno-2018) alphas and the alphas Rosales used to compute POL, not a phase convention.
- **The parities and canonical structure were re-derived and are all correct** (alpha_1,2,4,6,7
  even; alpha_3,5,8 odd — the alpha_3 odd-parity worry was unfounded), and the alpha_4/alpha_5
  ~2.15 amplitude is *physically right* (it is the indirect Sun term ≈ `mu_S/a_S² ≈ 2.176`,
  as the existing `test_qbcp_circular_limit_eom` already encodes). The x→−x reflection baked
  into the alpha tables is internally consistent (raw and reflected models agree exactly).
- **One genuine convention discrepancy found, deliberately NOT shipped.** `qbcp_eom`
  multiplies the *entire* Newtonian potential (Earth+Moon+Sun) by `alpha_6`, but the
  in-digest Rosales/Jorba (2023) Eq. 3 places `alpha_6` *only* on the Sun term (as
  `−m_S/(alpha_6·R_PS)`), Earth/Moon unscaled — physically `alpha_6` is a Sun-distance
  coefficient. Implementing the Eq. 3 form (Sun term `/alpha_6`, Earth/Moon unscaled)
  moved our L1 substitute ~30% *closer* to POL1 (2.12e-2 → 1.47e-2) **but regressed the
  Sun-Earth-L2 torus positive control 16×** (invariance residual 3e-5 → 4.2e-4, still under
  the 1e-3 gate but clearly worse). Since (a) it does not explain the dominant 2% gap,
  (b) it risks a validated result, and (c) which Hamiltonian form the *Gimeno-2018*
  `alpha_6` values were fitted for cannot be verified without that PDF, the change was
  **reverted** per "check, don't guess." This is the #1 item to settle when
  Gimeno/Jorba 2018 Table 4 + its Hamiltonian form are accessible. (Also flagged, not
  fixed for the same no-guess reason: `_COEFFS_ALPHA1[5] = -38.068581391005552e-08` has a
  mantissa >10 and breaks the harmonic-decay pattern — very likely a decimal-point slip for
  `-3.8069e-08`, but it is O(1e-7) and dynamically immaterial.)
- **Root-cause implication for the ORIGINAL #544 blocker (EM-L2 torus never converges).**
  It is the *same* instability artifact: `correct_qbcp_torus` runs single-period GMOS
  (`propagate_qbcp_pv` over one full `T_s` per invariant-circle sample), which is fine for
  the mildly-unstable SE-L2 torus but hopeless for the violently-unstable EM-L1/L2 region
  (each sample point is blown up by ~1e6–1e8). **The real fix is a multiple-shooting GMOS
  corrector** (sub-interval the period so no single leg spans the full unstable growth) —
  NOT more Fourier modes, a different target orbit, or coefficient hunting. #539/#540's
  reuse of the single-period corrector on Jupiter-Europa/Uranus will hit the identical wall
  wherever the target libration point is strongly unstable.
- **Shipped:** no change to `qbcp.py` (the model is validated structurally correct). Added
  two non-circular regression tests to `tests/core/test_qbcp.py`:
  `test_qbcp_collinear_instability_matches_cr3bp` (frozen L1/L2 rate vs the independent
  CR3BP Szebehely reference — a sourced structural golden that would catch a gross alpha/
  potential error) and `test_qbcp_pol_forward_prop_is_instability_dominated` (pins the
  forward-prop artifact so it is not re-mistaken for a bug). All four suites pass (20/20),
  ruff + mypy clean.

## Novel-orbit discovery proposals following #538 (allocated 2026-07-09, this session — read this before dispatching)

Formulated after auditing the #521-538 arc's actual outcomes (not just its headline resolutions): the coherent-model whiskered-torus/heteroclinic-connection pipeline (#522→#533→#537→#538) is this project's newest capability and sits in the "cislunar BCR4BP — under-mined, MEDIUM-yield" slot the 2026-06-26 `docs/superpowers/specs/2026-06-26-next-frontier-prioritization.md` frontier-ranking identified before it existed. #536's Jupiter-Europa screen tested exactly **one** Jacobi constant (C=3.0015) and found 0 connections — per the project's own "no X found is conditional on the search formulation" discipline, a single-point probe does not certify that region empty; it is a starting point, not a completed sweep. The proposals below sequence the natural next moves, ranked by expected new-catalogue-row yield and gated on #538 actually landing (either a confirmed closure or a documented negative) first, since #539/#540 reuse #538's corrector methodology directly.

**PARTIAL UN-GATING (2026-07-10).** A dispatched Fable advisor pass (standing in for the
unavailable `advisor` tool) reviewed the whole #521-544 arc and argued the "#539/#540 wait
for #538" gate is only true of the torus-*connection-corrector* half of #539/#540's scope —
the Jacobi-band + synodic-phase-offset *screening* half reuses pure CR3BP-level machinery
(`scripts/run_536_jupiter_europa_connection.py` imports only `core.cr3bp`/`genome.qp_tori`/
`genome.qp_torus_heteroclinic`, zero QBCP dependency) that already works today, independent
of #544's still-open QBCP-model bug. Its case: 7 months / ~356 catalogue rows / 1 confirmed
novel find (#312, Uranus) is a predictable result of an effort allocation skewed toward
corrector-depth on one stubborn cislunar target, when both #312 and the PC(3,2)
binary-genome find came from cheap breadth (pointing existing machinery at unscreened
systems). Splitting the screen out lets it run now and feed #539/#540's eventual corrector
work later, rather than sitting idle behind #544.

- **#545** — Jupiter-Europa/Ganymede band screen: a genuine Jacobi-constant band sweep (not
  #536's single point C=3.0015) crossed with a synodic-phase-offset sweep between the moons,
  using #524's deflated Newton to enumerate basins, CR3BP-level only. **Positive control:
  recover #534's own published Earth-Moon L1<->L2 result** (not a not-yet-existing SE<->EM
  connection) as the pipeline sanity check before pointing it at unmapped Jovian territory.
  Supersedes #536's single-point negative.
  **✓ ALREADY RUN (committed, `scripts/run_545_jovian_band_phase_screen.py`;
  `data/runlogs/545_band_screen_report.json`) — DEAD + SUPERSEDED, not "ready to dispatch."
  Corrected 2026-07-15 during a full-file audit** (the original 2026-07-13/14-era proposal text
  above, describing #524's deflated Newton as an independent method, undersold this task's real
  coupling — the ACTUAL committed script builds real QP tori via `genome.qp_tori.correct_qp_torus`
  and refines candidates against `genome.qp_torus_heteroclinic.closest_curve_distance`, i.e. it
  depends on the SAME whiskered-torus/linking-number machinery #534/#536 used, with #524's
  deflated Newton only as the root-finder layered on top — NOT a substitute for it). Two
  independent reasons this is closed, not just gated:
  1. **Dead per #555's final verdict**, exactly like #536: #576's own text explicitly lists
     "the now-retired torus/heteroclinic linking-number lane (#536/#545, dead per #555's final
     verdict)". The pipeline #545 depends on was shelved after #548 found it could never produce
     a confirmed nonzero linking number even at a correct positive-control energy.
  2. **Superseded by #576/#577 regardless**: #576 already ran the far more validated, twice-proven
     #563/#575 direct symmetric-closure method against ALL 6 Galilean moon pairs (including
     Europa-Ganymede, #545's specific target), and #577 delivered a definitive adjudicated
     verdict (0/36 novelty-clear, FULL STOP). #545's entire scientific question is already
     comprehensively answered by a better method — reviving it with any positive control would
     only duplicate #576/#577's completed work using an inferior, retired approach.
  **No further action needed on #545.** If a NEW Jovian-system question ever comes up, it should
  be scoped as a fresh task using #563's method, not a revival of this one.
- **#546** (P1, run now, independent of #538/#544) — Uranian-system band screen
  (Miranda-Ariel-Umbriel-Oberon lanes), same CR3BP-level band+phase sweep as #545, pointed
  at the one system that has ever produced a confirmed novel row (#312). Per the empty-region
  registry's re-open rule, a strictly-more-capable method (genuine band+phase sweep vs.
  whatever produced #312) is licensed to re-open this territory. **Recommended models:**
  system setup + parameter sourcing (digest-and-reconcile against a published Uranian
  reference, per `[[feedback_digest_not_adoption]]`) → **Sonnet**. Running the sweep and
  adjudicating any hit against the #312 anchor → **Opus**, **Fable** second-opinion pass.
  ✓ Attempted (2026-07-10, session C). Constants were already correctly sourced (no new work
  needed — `src/cyclerfinder/core/satellites.py` lines 159-178, JPL SSD gm_de440 planetary
  constants + satellite phys_par/sats-elem tables, ref URA111, accessed 2026-06-14; Miranda GM
  4.3, Ariel 83.5, Umbriel 85.1, Oberon 205.3 km^3/s^2, mean radii + SMAs also cited there).
  NO usable positive control exists anywhere in the
  qp_tori/qp_torus_heteroclinic method family (re-ran #534's own Earth-Moon L1<->L2 search to
  completion this session: 0/50 sign changes across all 3 crossing searches, matching #536's
  own documented 0-connection Jupiter-Europa result — neither prior application of this method
  has EVER found a genuine connection); #312 itself cannot serve as a positive control either
  (structurally different object — Lambert/patched-conic multi-arc quasi-cycler, not a
  libration-point torus-heteroclinic connection). Built a new literature-independent
  planar-Lyapunov torus bootstrap (`scripts/run_546_uranian_torus_screen.py`, since unlike
  Earth-Moon/Jupiter-Europa no published halo/NRHO seed table exists for any Uranian moon) and
  ran it: 18 Jacobi-band points across Umbriel/Oberon/Ariel (Miranda not reached, deprioritized
  — smallest GM, hardest numerically, time-boxed session), 1 completed torus-connection screen
  (Uranus-Oberon, C=3.003901, 8x8 grid, 0/64 sign changes), 17 points rejected pre-torus-stage
  (bootstrap non-convergence or wrong-branch topology guard). **Verdict: NOT a certified empty
  region** — per this project's own "verify a gauntlet with a positive control before trusting
  0/N" discipline, a 0/N from an unvalidated method is uninterpretable, so nothing is registered
  in `data/empty_regions.jsonl`. **Higher-priority follow-up surfaced**: before #541 (Saturn) or
  any further re-application of #545/#546-style screening, establish a genuine validated
  positive control for the qp_tori/qp_torus_heteroclinic linking-number method itself (find or
  construct ANY real transit/heteroclinic connection in a well-studied CR3BP system, e.g.
  resolving #534's own flagged-but-never-closed transit-vs-non-transit branch-classification
  problem) — right now the method is 0-for-3 with no validated positive anywhere in its history.
  **⚠ STAMP (2026-07-10, #548): the requested follow-up ran and the method is now SHELVED.** #548
  built the empirical transit-branch classifier and ran the reframed Owen & Baresi positive control
  at the achievable common energy band — zero sign changes, linking number identically 0, no
  torus-level positive ever obtained (now 0-for-4). This Uranus-Oberon 0/64 stays METHOD-INVALID —
  DO-NOT-CERTIFY (already correctly unregistered), and the qp_tori/qp_torus_heteroclinic pipeline is
  parked. Postmortem: `docs/notes/2026-07-10-postmortem-548-linking-number-pipeline.md`.

- **#547** (P0 follow-up to #546, 2026-07-10) — Transit-vs-non-transit manifold-branch positive
  control for the qp_tori/qp_torus_heteroclinic method family. **RESOLVED (positive control built)
  + the "no validated positive anywhere" framing corrected.**
  - **A pre-existing validated positive control was already in the tree, in a DIFFERENT method
    module.** `genome/heteroclinic_cycle.py` (task #314) certifies genuine planar-CR3BP transit
    heteroclinic connections via the classical KLMR section-crossing + Floquet-manifold method,
    golden-validated against Wilczak & Zgliczynski's computer-assisted proof of the closed L1<->L2
    Lyapunov cycle in the Sun-Jupiter-Oterma PCR3BP (arXiv:math/0201278). `tests/genome/
    test_heteroclinic_cycle.py` is green (8/8, re-run this session): `test_connection_l1_to_l2_
    converges` (residual < 1e-6) and `test_assemble_l1_l2_two_cycle_closes` (closed cycle,
    residual ~1.1e-9, independently Radau cross-checked). So the #546 follow-up's "0-for-3 with no
    validated positive ANYWHERE in its history" was scoped to the qp_torus LINKING-NUMBER pipeline
    specifically; genuine transit connections ARE findable and validated in this codebase — the
    transit-vs-non-transit branch classification #534 called an "open sub-problem" is in fact
    already solved in #314's `correct_connection` (neck-facing branch selection `branch_u=-1,
    branch_s=+1` + crossing-index `(k_u=3, k_s=4)`, documented in that docstring: "only those
    branches reach the L1-L2 neck where the connection lives").
  - **Built a fresh, minimal, from-first-principles transit positive control** (`genome/
    transit_manifold.py` + `tests/genome/test_transit_manifold.py`, 3 tests green) that resolves
    #534's flagged classification numerically and directly. An Earth-Moon L1 planar Lyapunov orbit
    at `C=3.1869` (deliberately between `C_L2=3.1722` and `C_L1=3.1883`, so the L1 neck is open and
    the L2 neck is closed — the classic bounded-Moon-realm KLMR gateway) is a strong saddle
    (`lam_u~2641`, x0=0.8321 near L1). Its `+` unstable branch is a genuine TRANSIT trajectory:
    threads the neck, crosses `x=1-mu` seven times, approaches the Moon to 0.023 nondim (~8800 km),
    reaching x up to 1.088 in the Moon realm. Its `-` branch is NON-TRANSIT: stays interior, swings
    back toward Earth to x=-0.60, never crosses `x=1-mu`. The branch is classified EMPIRICALLY
    (propagate both, see which reaches the section), not guessed from an eigenvector-component
    sign — exactly the missing piece #534 named. This is the codebase's first from-scratch
    transit-branch positive control keyed to the "propagate a signed manifold perturbation and check
    surface crossing" mechanic the qp_torus screen itself uses.
  - **Diagnosis of the linking-number method's 0-for-3 (honest, no forced positive).** #534's
    energy was NOT the problem: `C=3.15` is below both collinear thresholds (both necks open),
    matching Owen & Baresi's own EM demonstration exactly. The actual gap is that the
    linking-number pipeline (#522) has NEVER been run on its own paper's validated positive control
    — Owen & Baresi's EM quasi-halo<->quasi-halo case (`mu=0.012153643, C=3.15`, L1 latitudinal
    freq 0.2739 / L2 0.02163 → 4 connections, Fig 15, per `docs/notes/2026-07-03-digest-owen-
    baresi-2024-knot-theory-heteroclinic.md`). #534 built tori at the right ENERGY but abandoned
    reproducing those exact frequencies as impractical (its tori landed at omega_trans -0.1699 /
    -0.0356 — genuinely DIFFERENT tori whose manifolds are different objects), so its 0-connection
    result never tested the Owen-Baresi case at all. I did NOT find a branch-selection bug in
    `genome/qp_torus_manifold.py`: for the planar L1 saddle the unstable eigenvector's x-component
    is a healthy 0.30 (not near-zero — an earlier in-session read of "3e-7" was the `eps*vec[0]`
    scaled perturbation, not the eigenvector), so its default `vec[0]*sign` heuristic DOES pick the
    Moon-ward transit branch there; whether that heuristic degrades for a genuine 3D quasi-halo
    unstable eigenvector is untested (rebuilding #534's uncommitted tori was out of scope this
    session) and is a plausible-but-unconfirmed weak point, not a demonstrated defect.
  - **Bearing on prior negatives.** #534/#536/#546's 0-connection results remain method-artifact-
    contaminated NON-negatives (as #546 already correctly declined to register any empty region):
    the planar positive control now proves the underlying transit physics AND this codebase's
    manifold-crossing mechanic are sound, which ISOLATES the remaining gap to the
    linking-number-specific pipeline — faithfully reproducing an Owen-Baresi-grade isoenergetic
    quasi-halo torus PAIR (frequency-matched) and then validating the closed-curve/level-set/scan
    machinery against the published 4-connection count — before any qp_torus 0/N can be trusted as
    a certified empty region. That torus-reproduction + linking-number-pipeline validation is the
    real remaining #522-family blocker; it is NOT unblocked by this task, but this task removes the
    "does a genuine transit branch even exist / can we identify it" uncertainty that #534 left open.

**PORTFOLIO DECISION (2026-07-10, this session, per a second Fable-advisor pass reading the
full #521-547 arc) — #538/#544 and #539/#540/#541 as scoped are explicitly PARKED, not
cancelled.** Rationale: #538/#544's multiple-shooting GMOS corrector fix is well-characterized
and cheap to build later, but its payoff (closing the SE<->EM cislunar corridor) is
near-zero-novelty even on success — the most-published corridor in the field, four sessions
of corrector-depth investment already, zero catalogue rows. #539/#540/#541 as scoped inherit
that corrector and are built on the still-unvalidated qp_torus linking-number pipeline (see
#548 below) — re-scope after #548 resolves, and per #547's method-choice finding, probably as
#314-planar-method screens first (cheap, already validated) rather than qp_torus rebuilds.
#545/#546 remain the decoupled, already-run screen-only predecessors (committed, no connections
found, correctly uncertified). Two new threads take priority instead — #548 (settle the
linking-number pipeline's status once and for all) and #549 (the cheapest available shot at an
actual new catalogue row, since both of the project's two real novel hits came from validated
machinery pointed at unscreened real systems, not corrector depth on a known target).

- **#548** (P0, time-boxed to one Opus dispatch, hard kill criterion) — Reframed Owen &
  Baresi positive-control gate for the `qp_tori`/`qp_torus_heteroclinic` linking-number
  pipeline. NOT a repeat of #534's frequency-archaeology attempt (correctly abandoned as
  impractical — the paper's exact "latitudinal frequency" convention and seed orbit are
  unpublished). Instead: at the already-correct energy `C=3.15` (matches Owen & Baresi's
  own EM demonstration, per #547's diagnosis), sweep the torus family parameter
  (amplitude / rotation number k) on BOTH the L1 and L2 quasi-halo sides, bracketing the
  published frequency ratio (L1 0.2739 / L2 0.02163, `docs/notes/2026-07-03-digest-owen-
  baresi-2024-knot-theory-heteroclinic.md`), and run the linking-number scan across the
  resulting torus pairs — using #547's new EMPIRICAL transit-branch classifier
  (`genome/transit_manifold.py`) instead of `qp_torus_manifold.py`'s untested
  eigenvector-sign heuristic (flagged by #547 as the plausible-but-unconfirmed weak point
  for genuine 3D quasi-halo eigenvectors, unlike the planar case where it was confirmed
  correct). The ground truth being tested is "connections exist for quasi-halo pairs at
  this C and the scan sees them as sign changes" (Owen & Baresi report 4), NOT "our omega
  printout matches theirs" — #534's original, correctly-abandoned framing.
  **PRE-REGISTERED KILL CRITERION (binding, decide before dispatch, not after seeing
  results):** if a full both-family sweep at `C=3.15` with empirically-classified transit
  branches produces ZERO sign changes anywhere in the swept range, the linking-number
  pipeline is SHELVED — write a postmortem, mark the #522 family parked, and stamp
  #534/#536/#546's 0/N results as method-invalid-do-not-certify (not "empty region",
  "tool status: retired"). No further Fourier-mode/gauge/branch tinkering is licensed
  after a clean kill — this is a terminal gate, not another diagnostic round.
  **Recommended model:** Opus (trust-bearing numerical-methods judgment + the
  criterion-definition discipline this project keeps needing before sweeps).
  **✓ RESOLVED (2026-07-10) — SHELVED per the kill criterion (spirit), with two premise
  corrections found by primary evidence.** Full postmortem: `docs/notes/2026-07-10-postmortem-
  548-linking-number-pipeline.md`. **(1) The energy premise was wrong:** #534's committed NRHO
  seeds correct to **C=3.045**, not the C=3.15 #547 recorded — #534 was never at Owen & Baresi's
  demonstration energy. **(2) Exactly-C=3.15 isoenergetic quasi-halo pairs are impractical:**
  connections are ISOENERGETIC (both tori share one C), but the EM L1 halo family bifurcates from
  planar Lyapunov at C~3.146 (so 3.15 is at/above the L1 quasi-halo regime) and the L2 halo family
  reached via the NRHO branch tops out at a genuine high-C fold ~C=3.087 (confirmed by
  pseudo-arclength continuation past the x0-fold). The highest COMMON energy both families
  robustly reach is C in [3.05, 3.087], so the positive control was built there (same physics —
  both necks open, unstable quasi-halo pairs). **Built** `src/cyclerfinder/genome/qp_torus_transit.py`
  (`transit_torus_manifold_grid`: adapts #547's empirical transit classifier to the 3D quasi-halo
  torus grid — propagates BOTH signed perturbations per point, keeps whichever reaches the section
  first, discarding the untested `vec[0]*sign` heuristic; + `tests/genome/test_qp_torus_transit.py`,
  2 green) and `scripts/run_548_owen_baresi_positive_control.py`. **Sweep result:** 3 isoenergetic
  geometry-usable pairs (C=3.05/3.06/3.07; C=3.08 L2 halo at the fold, bisection failed), 12 linking
  scans across 4 scanning-variable/curve-triple specs — **ZERO sign changes; the linking number was
  identically 0 in every scan** (the reduced stable/unstable curves never link at any D, any scan
  variable, any energy). Manifolds DID reach the section and overlap (L1 320/320, L2 166-320/320,
  ranges overlapping), and the L1 frequency ratio bracketed the published 0.2739 (swept
  0.3668→0.3300→0.2944 over C), so the pipeline had a fair shot. The `linking_number` PRIMITIVE is
  sound (returns ±1 on a Hopf link — `tests/search/test_linking_number.py` green), so the
  identically-0 output is a true property of the extracted manifold curves, not a broken invariant.
  **Honest caveats (adjudicator to weigh):** the literal "C=3.15" precondition was unsatisfiable;
  the L2 frequency bracket was poor (0.12-0.45 vs published 0.02163 — the paper's L2 quasi-halo is a
  near-planar small-z near-bifurcation halo this NRHO-branch machinery can't reach; a future
  L2-planar→halo bifurcation seed generator could give one final fully-frequency-matched shot); and
  the L2 stable grid showed some transit-sheet mixing (~150/170 sign balance, because the L2 halo
  sits almost on the section). This is a shelve-on-persistent-failure (now 0-for-4) with a
  reversible caveat, NOT a proof of non-existence. **Per the binding pre-registration, no further
  gauge/mode/branch tinkering is licensed; the #522 family is PARKED and #534/#536/#546 are stamped
  method-invalid-do-not-certify below.** Did NOT write to `data/catalogue.yaml` (no hit; nothing to
  write). Recommend an Opus/Fable second pass to confirm the shelve or authorize the one L2-near-
  bifurcation retry before permanent retirement.
  **⚖ ONE-MORE-SHOT AUTHORIZED (2026-07-10, #553 adjudication) — the shelve is NOT ratified as
  final; #548's substitution rested on a factual error.** Independent verification (vertical
  stability index a_v of the PLANAR Lyapunov family crossing +1 — a different route from #548's
  halo z0->0 continuation; script + log in the #553 session scratchpad, `verify553_bif3.py`):
  the EM **L1 halo bifurcation is at C = 3.1744, not ~3.146** (the L2 value 3.1521 CONFIRMS
  #548's ~3.152, validating the method). C=3.15 is therefore INSIDE the L1 quasi-halo regime —
  #548's own L1_TABLE corroborates this (z0^2-vs-C extrapolation from its rows lands the
  bifurcation near ~3.17; z0=-0.061 at C=3.14523 is nowhere near a pitchfork) — so the ~3.146
  "snap to planar" was a corrector/continuation artifact misread as physics. Owen & Baresi's
  demo is self-consistent with the corrected numbers: at C=3.15 the L1 quasi-halo is developed
  (lat. freq 0.2739) and the L2 halo sits 0.002 below ITS bifurcation (tiny near-planar, lat.
  freq 0.02163). Consequently the C in [3.05,3.07] substitution was not "the same physics": it
  dropped the one property that gave the kill criterion its force — KNOWN ground truth. Nobody
  knows connections exist for the specific tori tested there, so the zero is a FOURTH
  unscreened-application negative, not a failed positive control; per
  [[feedback_verify_gauntlet_with_positive_control]] it cannot trigger the kill. Two further
  oversell findings: (a) the postmortem's "straddling 0.2739" claim is false — the swept L1
  ratios 0.2944-0.3668 sit entirely ABOVE 0.2739, so BOTH sides were frequency-mismatched;
  (b) "linking number identically 0" conflates computed-zero with curve-extraction failure
  (`scan_linking_number` emits 0 when `_first_closed_curve` returns None; the run JSON records
  no per-D curve-availability counts, and the L2 grid had 166/320 NaN + ~0.001-wide z-overlap +
  150/170 sheet mixing), so "the curves never link" was never established as a property of
  actually-extracted curves. AUTHORIZED final shot (one dispatch, time-boxed): (1) L1 quasi-halo
  at C=3.15 via the EXISTING machinery with careful continuation across the ~3.146 corrector
  failure (z0~0.05 there — not degenerate); (2) the flagged L2 planar->halo bifurcation seed
  generator (the #553 verification script already locates the bifurcation orbit: amp=0.0353,
  x0=1.1204, T=3.4156, C=3.1521 — step off in z and correct); (3) rebuild tori at C=3.15
  targeting ratios ~0.2739 / ~0.02163 and re-run the scan WITH per-D both-curves-available
  instrumentation plus a synthetic linked-curve positive control of the extraction step itself.
  RE-REGISTERED kill criterion (now with a valid precondition): zero sign changes at the
  genuinely frequency-matched C=3.15 pair -> permanent retirement, no appeal. What STANDS from
  #548: the #534/#536/#546 method-invalid-do-not-certify stamps (conservative, correct either
  way), the empirical transit classifier (`qp_torus_transit.py`, a real improvement), and the
  0-for-4 record. Standing rules compelling this over immediate ratification:
  [[feedback_never_give_up_reproducing_papers]] (O&B's 4 EM connections are a published result
  this codebase has never reproduced) — and the cost is ~one session of tool VALIDATION gating
  the interpretability of three already-run screens, not more corrector-depth discovery.
  **⏹ FINAL RESOLUTION (2026-07-11, #555 — see the #555 entry below for the full record):**
  the authorized shot ran. #553's bifurcation correction was INDEPENDENTLY CONFIRMED (L1
  halo bifurcation C=3.1745 by two routes, L2 C=3.1521 with a_v=1.00007) and #548's ~3.146
  REFUTED; both C=3.15 halos were built (breaking #548's structural wall); the pipeline
  MACHINERY was validated (synthetic linked-curve positive control + per-D curve-
  availability instrumentation), and the scan gave 0 sign changes with NON-TRIVIAL
  availability (transit grids 320/320 crossings, 78/80 both-curves-available — resolving
  #548's "identically 0" ambiguity: a genuine non-link, not missing curves). The clean kill
  is NOT triggered because the L1 quasi-halo could not be built at O&B's freq 0.2739
  (energy-pinned ~0.074 at C=3.15; larger tori intractable in `correct_qp_torus`) — so the
  frequency-match precondition is only half-satisfiable. TERMINAL: linking-number screen is
  sound and its lever exhausted; the one remaining blocker is a large-rotation quasi-halo
  torus corrector (a new capability). The #534/#536/#546 method-invalid-do-not-certify
  stamps STAND.

- **#553** (P1, cheap, read-only review) — Fable ratification pass on #548's shelve verdict.
  #548's own pre-registered kill criterion was written for "zero sign changes at C=3.15"; the
  actual test ran at C in [3.05, 3.07] after #548 found C=3.15 physically unreachable for
  isoenergetic EM quasi-halo pairs (L1 halo bifurcates ~3.146, L2 halo tops out at a fold
  ~3.087) — a reasonable, honestly-documented substitution, but a substitution nonetheless, not
  the literal pre-registered condition. Per this project's "Fable second-opinion before treating
  a significant judgment call as final" convention (same pattern as catalogue-writeback
  adjudication), get an independent read on: (a) whether the C-range substitution was a faithful
  stand-in for the original criterion or a goalpost-move that deserves more scrutiny; (b) whether
  the recorded caveats (unreachable L2 frequency bracket, some L2 transit-sheet mixing) are
  merely noted or actually undermine the "fair test" claim; (c) whether to ratify the shelve as
  final, or authorize the one flagged remaining option (an L2-planar-to-halo bifurcation seed
  generator, giving one genuinely frequency-matched final shot) before permanent retirement.
  Read-only review of `docs/notes/2026-07-10-postmortem-548-linking-number-pipeline.md` +
  `scripts/run_548_owen_baresi_positive_control.py` + `genome/qp_torus_transit.py` — no new
  sweep, no code changes expected.
  **Recommended model:** Fable (independent adversarial second opinion, matching the project's
  existing pre-writeback pattern).
  **✓ RESOLVED (2026-07-10) — ONE-MORE-SHOT AUTHORIZED, shelve NOT ratified as final.** Full
  adjudication appended to the #548 entry above. Key findings: (a) #548's L1 halo-bifurcation
  claim (~3.146) is wrong — independently recomputed at C=3.1744 via planar-Lyapunov vertical
  stability (its L2 ~3.152 value was confirmed at 3.1521 by the same method), so C=3.15 IS
  inside the L1 quasi-halo regime and the substitution's physics rationale fails on the L1
  side; (b) the substituted band lost the known-ground-truth property the kill criterion
  required, making the run a fourth unscreened negative rather than a positive control; (c) the
  postmortem's "straddling 0.2739" claim is false (swept 0.2944-0.3668, all above it) and
  "identically 0" conflates computed-zero with curve-extraction failure (no per-D availability
  instrumentation). One final frequency-matched shot at C=3.15 is authorized with a
  re-registered kill criterion; see the #548 append for the scoped work list.

- **#555** (P0, user-authorized "one more shot", FINAL — no further appeal after this)
  — the #553-scoped, correctly-frequency-matched Owen & Baresi test at the genuine C=3.15.
  Per #553's adjudication (which independently recomputed the L1 halo bifurcation at
  C=3.1744, not #548's ~3.146 — so C=3.15 IS inside the L1 quasi-halo regime after all,
  reachable via careful continuation across #548's corrector-failure region near
  z0~0.05, not a genuinely new capability), do exactly the three things #553 scoped:
  (1) build the L1 quasi-halo at C=3.15 with existing machinery, continuing carefully
  through the corrector-failure region rather than trusting the naive continuation #548
  used; (2) build the L2 planar-to-halo bifurcation seed generator #553 already
  half-located (bifurcation orbit: amp=0.0353, x0=1.1204, T=3.4156, C=3.1521 — step off
  in z and correct) to reach the genuine small-amplitude L2 halo near C=3.15 that #548's
  NRHO-branch machinery structurally could not reach; (3) rescan at C=3.15 targeting the
  published frequency ratios (L1 0.2739, L2 0.02163), with per-D BOTH-CURVES-AVAILABLE
  instrumentation (fixing #553's flagged ambiguity — #548's "identically 0" conflated
  genuine non-linking with silent curve-extraction failure on NaN-heavy grids) and a
  synthetic linked-curve positive control of the extraction/scan step itself, so a "0
  sign changes" result this time is provably about real geometry, not missing data.
  **Re-registered kill criterion (binding, FINAL): zero sign changes at the genuinely
  reachable C=3.15, with per-D curve availability confirmed non-trivial and the
  extraction-step positive control passing, means the linking-number pipeline is
  PERMANENTLY RETIRED — no further appeal, no third adjudication round.** If it finds a
  connection, this is a genuine, publication-adjacent result (independently reproducing
  Owen & Baresi's own reported case) and should go to Opus/Fable adjudication before any
  catalogue-adjacent writeback (a torus connection alone is not itself a catalogue row,
  but would validate the whole #522 family for #539-541's future reuse).
  **Recommended model:** Opus (this is the third round of trust-bearing numerical-methods
  judgment on this specific question; per this project's model-tiering policy, do not
  downgrade tier just because it's a "final" pass).
  **✓ RESOLVED (2026-07-11) — the final frequency-matched shot RAN; TERMINAL, no further
  appeal.** Script `scripts/run_555_owen_baresi_c315_final.py`; bifurcation-verification
  scripts + logs in the session scratchpad (`verify555_bif*.py`, `derisk3.py`); results
  JSON in scratchpad. Six findings:
  **(1) Bifurcation independently re-verified — #553 CONFIRMED, #548's ~3.146 REFUTED.**
  EM L1 planar-Lyapunov→halo bifurcation at **C=3.1745** via the vertical stability index
  a_v of the PLANAR Lyapunov family crossing +1, by TWO independent routes (bounded
  fixed-Jacobi x0-continuation → 3.17455; independent-amplitude a_v(C) sampling → 3.17425)
  — matches #553's 3.1744 to 3 decimals; #548's ~3.146 was the corrector/continuation
  artifact #553 diagnosed. EM L2 bifurcation at **C=3.1521** (a_v=1.00007 at C=3.15210,
  x0=1.12035, T=3.4157) — independently reproduces #553's located bifurcation orbit
  (amp 0.0353, x0 1.1204, T 3.4156, C 3.1521) to 4 decimals. So C=3.15 IS inside the L1
  quasi-halo regime (0.024 below its bifurcation) and 0.0021 below the L2 halo bifurcation
  — genuinely reachable, contra #548.
  **(2) Both halos built AT C=3.15 (#548's structural wall broken).** L1 via x0-natural
  continuation of the halo family THROUGH the pitchfork region (z0=-0.0545 at C=3.1504) —
  #548's failure was secant-on-x0-to-hit-C, ill-conditioned at the pitchfork; straight
  x0-continuation is not. L2 via the #553-scoped planar→halo bifurcation SEED GENERATOR
  (planar L2 Lyapunov at C=3.15 on the genuine near-libration family — ydot0>0, T~3.42, NOT
  the large-amplitude NRHO branch #548 was stuck on — then step off in z0 → genuine small-z
  3D halo z0=0.0171 at C=3.1497), reaching the near-planar small-z L2 halo #548's
  NRHO-branch machinery structurally could not (that branch tops out ~C=3.087).
  **(3) Frequency match is HALF-satisfiable — a genuine, unanticipated finding.** The L2
  quasi-halo latitudinal frequency lands at **0.0214, MATCHING O&B's 0.02163** (validates
  the rotation-number method AND convention). But the L1 quasi-halo latitudinal frequency
  at C=3.15 is ENERGY-PINNED near **0.074** (C=3.15 sits only 0.024 below the L1
  bifurcation → small latitudinal libration), FLAT across small torus amplitudes, and
  larger-amplitude L1 tori that might carry O&B's 0.2739 are computationally INTRACTABLE /
  non-convergent in `correct_qp_torus` (a single amp=0.02 build did not finish in 250 s;
  invres degrades past amp~0.01). So O&B's specific L1 quasi-halo (freq 0.2739) is NOT
  reproducible at C=3.15 with this codebase — a NEW, characterized blocker distinct from
  #548's halo-continuation wall. Corollary: #548's "ratio set by energy not amplitude" was
  CORRECT for L1, and #553's premise that C=3.15 automatically yields L1 freq 0.2739 does
  NOT hold under the very convention that reproduces O&B's L2 value exactly.
  **(4) Pipeline MACHINERY validated (fixes #553's ambiguity).** `LinkingScanResult` now
  carries per-D stable/unstable curve-availability + `availability_summary()`; a synthetic
  offset-Hopf-link positive control run through the SAME `scan_linking_number` /
  `closest_curve_distance` code path correctly detects the known nonzero linking (unlinked
  control → identically 0; both fully available). New tests in
  `tests/genome/test_qp_torus_heteroclinic.py`.
  **(5) Scan at genuine C=3.15 (L2 freq-matched, L1 near-bifurcation 0.074): 0 sign
  changes, and — decisively — availability is NON-TRIVIAL.** Transit grids 320/320
  crossings on BOTH manifolds (vs #548's NaN-heavy grids); all 4 scan specs show
  both-curves-available 78/80 with linking number identically 0. So the zero is now
  PROVABLY a property of actually-extracted curves that genuinely do not link — the exact
  distinction #553 said #548 never established. Positive control PASS.
  **(6) TERMINAL VERDICT — QUALIFIED NEGATIVE, no further appeal.** The #522 linking-number
  SCREEN is now VALIDATED end-to-end (synthetic positive control + availability
  instrumentation prove it detects real links and reports genuine non-links), and the
  genuine C=3.15 regime is reachable (both bifurcations confirmed, both halos built) — the
  two things #548/#553 left open. For the ACHIEVABLE frequency-matched-where-possible C=3.15
  pair (L2 exact, L1 near-bifurcation), there are ZERO connections with full curve
  availability. This is NOT the clean "the method finds nothing that exists" retirement #553
  pre-registered, because the L1 side could not be built at O&B's 0.2739: the re-registered
  kill criterion's frequency-match precondition is only half-satisfiable. The residual gap
  to a FULL O&B L1↔L2 reproduction is therefore a single, precisely-localized NEW blocker —
  `correct_qp_torus` cannot build the large-rotation-number (0.2739) L1 quasi-halo near its
  bifurcation — NOT a defect in the linking-number screen. No further linking-number-screen
  tinkering is licensed (that lever is exhausted); the one remaining lever is a
  large-rotation quasi-halo torus corrector, a NEW capability out of scope for "one more
  round of the same". The #534/#536/#546 method-invalid-do-not-certify stamps STAND
  (their tori were never frequency-validated either). Did NOT write to `data/catalogue.yaml`
  (no connection found; nothing to write).

- **#556** (parking lot — flagged for user review before dispatch, NOT auto-fired) —
  large-rotation-number quasi-halo torus corrector. The single, precisely-localized blocker
  #555 found to a full Owen & Baresi L1<->L2 reproduction: `genome.qp_tori.correct_qp_torus`
  cannot converge the EM L1 quasi-halo at the amplitude needed to carry the published
  latitudinal frequency 0.2739 near C=3.15 (a single amp=0.02 build did not finish in 250s;
  invariance residual degrades past amp~0.01) — a genuine numerical-methods capability gap,
  not a conceptual one (the L2 side, near its own bifurcation, already converges cleanly and
  matches 0.02163 exactly). Building this would let the linking-number screen finally run
  the ACTUAL Owen & Baresi case rather than a half-frequency-matched substitute, either
  closing the loop on the #522 family for good (a genuine connection = a publication-adjacent
  result validating #539-541's future reuse) or reaching a truly unappealable clean negative.
  **Not auto-fired** because: (a) this is the fourth dispatch in the #548->#553->#555 chain,
  each of which found a real, non-trivial issue — a pattern that could keep finding "one more
  thing" indefinitely without a clear stopping rule this time; (b) per this session's own
  portfolio-parking logic (see the #538/#544/#539-541 PARKED note above), the whole SE<->EM
  cislunar corridor has now absorbed five-plus tasks across multiple sessions with zero
  catalogue rows, while breadth-on-unscreened-systems work (#549, #312) is what has actually
  paid off; (c) a torus corrector convergence fix is genuinely open-ended numerical-methods
  work (could be a quick fix — better initial guess, arclength continuation in amplitude — or
  a deep one — the corrector's Newton scheme may just not handle this rotation-number regime
  at all), unlike #555's own bounded, well-specified final round.
  **Recommended model:** Opus for the initial diagnosis (is this a quick continuation-strategy
  fix or a deeper corrector limitation?) before committing to a build; Sonnet for
  implementation once scoped, behind the existing SE-L2/L2-near-bifurcation convergence tests
  as a regression floor.

- **#549** (P0, fire alongside #548, fully disjoint code path) — Real-binary `(k1, k2)`
  genome sweep. The Tier-1 forward-plan line "circumbinary Pluto-Charon deeper + other-binary
  (k1,k2) sweep" (identified 2026-07-01, never allocated a task number, never run). Reuse
  the #494/#504 binary-cycler genome harness verbatim. **Positive control:** re-find PC
  (3,2) at `mu=0.10876` (Pluto-Charon's own mass ratio) before trusting anything new.
  Then sweep `(k1,k2)` at sourced mass ratios for real binary systems: **Patroclus-
  Menoetius** (mu~0.19 near-equal-mass binary, timely — Lucy flyby 2033), **Didymos-
  Dimorphos** (post-DART, precisely characterized), **Orcus-Vanth**, **Eris-Dysnomia**.
  Source each mass ratio from a citable reference (JPL SBDB / published radar-astrometry
  papers), not an assumed value — per `[[feedback_digest_not_adoption]]`. Any stable
  family member found + `search.literature_check` returning not-found = a fresh V1
  instantiation row, exactly the PC (3,2) admission pattern (this project's other actual
  novel-ish hit). This is the single most likely new catalogue row available this week —
  same validated-tool-on-unscreened-real-system pattern behind both of the project's two
  genuine hits (#312 Uranus, PC (3,2)), unlike single-moon torus connections which are
  tour ingredients, not catalogue rows, even on success.
  **Recommended models:** harness reuse + sweep execution behind the literature-check
  ratchet gate → **Sonnet** (spec-complete, strong deterministic gate). Any hit's
  admission verdict → **Opus**, **Fable** second-opinion pass before catalogue writeback
  (standard pattern for this project).
  **✓ RUN (2026-07-10) — CLEAN NEGATIVE across all four systems; no catalogue writeback.**
  **Positive control PASSED**: `sweep_32_positive_control()` (#504's own function, unmodified)
  re-finds PC (3,2) at `mu=0.10876473603280369` matching the committed `ross-rt-pc-cycler-32-2026`
  row to 9 significant figures (C=3.579515019729634 vs catalogue 3.57951501972907; x0=
  -0.693198287043394 vs -0.693198287043369; T=11.833462517014365 TU vs 11.8334625170346;
  independent-Radau dJ=7.96e-13 vs catalogue's 8.0e-13) — the sweep machinery is trustworthy.
  **Mass ratios sourced** (correcting this entry's own eyeballed guesses, per
  `[[feedback_digest_not_adoption]]`): **Patroclus-Menoetius** mu=0.4381 (NOT ~0.19 as originally
  guessed above — Buie et al. 2015, AJ 149, 113, occultation diameters D_Patroclus=113+/-3 km,
  D_Menoetius=104+/-3 km, equal-density assumption; orbital scale from Buie et al. 2024, AJ 167,
  104, a=692.5+/-4.0 km, P=4.282754+/-0.000023 d); **Didymos-Dimorphos** mu=0.0079 (Naidu et al.
  2020, Icarus 348, 113777, pre-impact system mass (5.4+/-0.4)e11 kg + equal-density Dimorphos
  mass ~4.3e9 kg; CAVEAT: Dimorphos's mass is not directly measured — post-DART papers give a
  density range 1500-3300 kg/m^3 vs the assumed 2170, i.e. mu could plausibly be 0.0055-0.012;
  ESA Hera 2026/2027 will pin it down); **Orcus-Vanth** mu=0.137+/-0.013 (Brown & Butler 2023,
  PSJ 4, 178, arXiv:2307.04848, direct ALMA astrometric mass measurement: "Vanth contains
  13.7+/-1.3% of the mass of the system"); **Eris-Dysnomia** mu=0.00498 (same Brown & Butler 2023
  paper; CAVEAT: only a 1.5-sigma mass DETECTION, 1-sigma upper bound mu=0.00833, the weakest-
  sourced of the four). Full citations + physical scales (a, P) in
  `src/cyclerfinder/search/real_binary_kk_sweep.py::REAL_BINARY_SYSTEMS`.
  **Sweep**: reused #504's machinery verbatim (fixed-Jacobi symmetric corrector, Barden stability,
  winding-topology classifier, independent-Radau crosscheck) via a new thin driver module
  (`search/real_binary_kk_sweep.py`) that fixes a latent bug in #504's `mu_step_to_orbit` (its
  "final correction" step was hardcoded to `make_pluto_charon_system()` regardless of the
  `target_mu` argument — harmless for `target_mu==PC_MU` but would have silently corrupted results
  for any other system; `mu_step_to_system` takes the target system explicitly instead). Swept the
  same six `(k1,k2)` topologies #504 swept — (1,1)/(2,1)/(2,2)/(3,1)/(3,2)/(3,3) — at each of the
  four sourced mass ratios: anchor-seeded mu-continuation from the Ross-RT 2026 Table-I anchors for
  (1,1)/(3,1)/(3,2)/(3,3) (with (1,1) tried from all three available anchors, mu=0.001/0.01215/0.5),
  and a bounded (x0,C,hc) grid search with a per-call SIGALRM timeout for (2,1)/(2,2) (not in
  Table-I, exactly #504's own fallback). **32/32 probes (24 anchor + 8 grid) returned clean
  negatives** — no stable prograde `(k1,k2)` cycler at any topology, any sourced mass ratio.
  Failure modes: mu-continuation branch folds before reaching the target mu (most common), wrong
  topology at the stable point found, no `|nu|<1` window in the C-sweep range, or (Didymos-Dimorphos
  (3,1)) the seed's Jacobi constant exceeding the target's `C_L1` outright (the family cannot exist
  in that Hill region at that mu). **Robustness check**: since many failures were "mu-continuation
  did not converge" (a branch-fold signature, not a C-sweep negative), reran those 10 cases at
  2.5-10x finer mu-step resolution (each step <=0.001 in mu, up to 426 steps) to rule out a
  step-size artifact. 6/10 reproduced the identical negative at finer resolution (confirms genuine
  branch folds). The remaining 4 (Patroclus-Menoetius (1,1)/(3,1)/(3,2), Eris-Dysnomia (3,2)) did
  not converge within a 240s per-job cap even at fine resolution — genuinely INCONCLUSIVE at the
  compute budget spent, not flipped positive and not as firmly confirmed-negative as the other 28;
  flagged rather than silently folded into the clean-negative count. **No literature-novelty check
  was run** (`search.literature_check` gates candidates, not negatives — nothing to check).
  **Two bugs found + fixed during test-driven verification (post-hoc, per
  `[[feedback_bugfix_invalidates_past_searches]]`), both re-checked against the discovery-run
  results before trusting them:** (1) `sweep_family`'s C-sweep call hardcoded `hc=None` instead of
  threading the anchor's own `half_crossings` through (matters only for the `(3,2)` anchor, which
  needs `hc=6`) — re-running all four systems' `(3,2)` topology after the fix reproduced the IDENTICAL
  negative verdict in all four cases (3 never even reached the C-sweep stage; Orcus-Vanth did, and
  got the same negative before and after). (2) the default symmetric C-sweep band walked `hc=6`
  down 0.3 below the anchor's C, which is pathologically slow/non-convergent for that family (mirrors
  why #504's own `sweep_32_positive_control` deliberately sweeps upward-only from the anchor); fixed
  to match that convention. Neither fix changed any reported verdict — both are logged for the
  record, not as a retraction. **Verdict**: Pluto-Charon's (3,2) family looks idiosyncratic to its own (mu, C) window rather
  than reflecting generic structure accessible at other real binary mass ratios; PC (3,2) remains
  this project's only confirmed real-binary novel-ish hit. New code: `src/cyclerfinder/search/
  real_binary_kk_sweep.py`, `scripts/run_549_real_binary_kk_sweep.py`,
  `tests/search/test_549_real_binary_kk_sweep.py`.

- **#550** (P3, opportunistic, cheap) — PC (3,2)'s V2->V3 promotion is gated on a single
  missing NAIF kernel (`SAT441l`, per #506's writeup) — a fetch from naif.jpl.nasa.gov, no
  new logic. Raises the tier of the project's newest confirmed novel row.
  **Recommended model:** Haiku (pure mechanical download + kernel-load verification).
  **✗ STALE — ALREADY DONE (2026-07-10, discovered while about to dispatch this).** This task
  was allocated from a Fable-advisor pass that read #506's entry but missed that the SAME
  session (2026-07-01) already fully resolved it: task **#510** fetched `plu060.bsp` (129 MB,
  `~/dev/references/kernels/plu060.bsp`, verified this session to still load cleanly via
  `verify.spice_kernels.ensure_pluto_kernel()`), and task **#511** then built and RAN the actual
  differential-correction V3 promotion attempt (`verify/pluto_charon_realeph.py` +
  `scripts/pc_v3_realeph_correction.py`, `tests/verify/test_511_pc_realeph_correction.py`, 4/4
  pass). **Verdict (already recorded, `docs/notes/2026-07-01-511-pluto-charon-realeph-verdict.md`):
  PC (3,2) STAYS V2-ballistic** — no real-ephemeris analog of this exact orbit exists, for a
  structural reason: its period `T/(2*pi)=1.8834` is not commensurate with Charon's own orbital
  period (the `#441` "period_f trap"), confirmed on real SPICE data (independent-closure residual
  1.726e-2 nd = 338 km, 9+ orders above the 1e-8 gate — not a marginal near-miss). This is a
  genuine physics finding, not a missing-kernel gap; no further action needed on this task. The
  corrector/e-continuation infrastructure remains reusable for any future `k:1`-resonant PC row.
  No dispatch fired.

- **#551** (P3, opportunistic, now time-relevant) — GTOC 13 methods-paper corpus mining.
  GTOC 13 ran Oct-Nov 2025 on ballistic Jovian gravity-assist tours — this project's exact
  problem class — and competition methods-paper write-ups are typically published mid-year
  following the competition, i.e. approximately now. Proposed 2026-07-02, dormant since.
  Competition write-ups routinely contain unpublished trajectory-family constructions —
  potential corpus-mining gold and possible V0 catalogue admissions.
  **Recommended model:** Sonnet for the search/fetch/digest pass (standard per-paper digest
  discipline, `[[feedback_per_paper_digest_todo]]`); Opus only if a candidate admission
  needs a novelty/tier judgment call.
  **✗ STALE PREMISE — mostly already checked, DO NOT dispatch as scoped (2026-07-10,
  discovered while about to dispatch this).** This entry repeats the same "ballistic Jovian
  gravity-assist tour" mischaracterization task **#526** already corrected on 2026-07-03: GTOC
  13 is actually a ballistic-flyby + optional ideal-solar-sail tour of a FICTIONAL single-star,
  10-planet exoplanetary system ("Altaira") — no Jovian moon system, no cycler-relevant
  repeating structure even in principle (a one-shot 200-year tour maximizing a score). #526
  also already searched WebSearch/arXiv/the GTOC13 site/the ESA GTOC-portal mirror for methods
  papers and found an honest negative — none exist yet, workshop still listed "TBA", consistent
  with the typical 6-18-month competition-close-to-journal-special-issue pipeline for a
  competition that closed Nov 2025. Re-dispatching the identical search only 7 days later
  (2026-07-10) would almost certainly reproduce the same negative at the cost of a wasted
  agent run. **Recheck condition:** not worth revisiting before ~Q1 2027 (12+ months post-close,
  giving the workshop→journal pipeline real time to produce something), and even then the
  problem-statement digest (`docs/notes/2026-07-03-digest-gtoc13-problem-statement.md`) already
  establishes there's no cycler-relevant structure to mine regardless of methods-paper
  availability — so this is genuinely low-priority even at its recheck date. No dispatch fired.

- **#552** (parking lot — flagged for user review before dispatch, NOT auto-fired) — 3D/
  inclined-releg extension to the moontour genome. The empty-region registry's own re-open
  rule licenses this: Jupiter-Amalthea's registered empty verdict lists 3D/inclined relegs
  as its own re-sweep condition, and the Neptunian system's "structurally empty" verdicts
  are all coplanar-prograde machinery applied to a system whose massive moon (Triton) is
  RETROGRADE — a known, named capability gap, not a swept-and-confirmed absence. This is
  the system CLASS where the project's only confirmed novel hit (#312, Uranus) lives, so a
  capability extension here plausibly re-opens more territory than a fresh system screen
  would. Genuinely multi-week Track-A-style work, not an errand — per
  `[[feedback_speculative_high_effort_required]]` this kind of investment is in-policy, but
  the Fable advisor pass explicitly flagged this as the user's call to make at a
  consolidation review, not something to auto-dispatch alongside #548/#549.
  **Recommended model:** Opus for the capability-scoping/design pass if greenlit; Sonnet for
  implementation behind whatever regression gate that design settles on.
  **✗ SCOPED, NOT GREENLIT (2026-07-10, Fable scoping pass).** The capability gap is real —
  verified in code, not just doc claims: `core/satellites.py` has no inclination/node/direction
  field at all (Triton's retrograde nature exists only in comments), `search/discovery_
  campaign.py::_moon_state()` hardcodes `z=0` + prograde motion, every moontour sweep flows
  through it. But **the payoff claim does not survive a back-of-envelope flyby-bend check the
  original proposal never ran.** Amalthea (GM=0.165) is MASS-limited, not geometry-limited — max
  bend 0.80 deg at V_inf=0.5 km/s, collapsing further at higher V_inf; #433's own verdict already
  recorded literature-fresh Amalthea closures failing the bend gate. Its actual inclination to
  Jupiter's equator is ~0.37 deg anyway (JPL SSD) — 3D relegs add zero bending capacity there; the
  registry's OTHER named re-sweep condition (low-thrust relegs) is the real lever, not this one.
  Triton is worse: the coplanar-PROGRADE model the code actually uses is GENEROUS to Triton
  (prograde encounter allows low V_inf where GM=1428 bends 23-59 deg); modeling it CORRECTLY as
  retrograde forces V_inf up to ~6-9+ km/s where the bend collapses to ~1-2 deg — the registered
  Neptune "structural empty" verdicts would get STRONGER, not reopen, under a correct model. And
  #312 itself was found WITH the existing coplanar machinery on prograde, low-inclination regular
  Uranian moons (`catalogue.yaml`: `orbit_fidelity: circular-coplanar`) — the "novel hit lives in
  this capability gap" framing conflated system class with capability gap; the hit lives exactly
  where the coplanar machinery already works. **Verdict: do not dispatch the multi-week build.**
  **Cheap alternative surfaced (~1 day, not yet dispatched):** a "retrograde-correction stamp"
  task — run the bend numbers above formally and append `reverification` entries to the 3 Neptune
  rows + the Amalthea row in `data/empty_regions.jsonl`, converting them from "conditionally
  empty pending a named capability" to near-unconditional negatives (Amalthea's live re-sweep
  condition narrows to low-thrust only). **If 3D relegs are ever wanted for their own sake**,
  re-scope around Iapetus/Saturn instead (GM=120, i~15 deg — a real, currently-ignored
  plane-change cost with no registered re-open claim resting on it) — genuinely small (3-5 days,
  not multi-week) since `core/lambert.py` already supports 3D/retrograde branch selection; the
  #515-518 "3D Lift framework" is NOT reusable here (it's SE<->EM libration-point cross-system
  machinery, a different problem domain from planet-centric patched-conic relegs).
  **REVIVAL ASSESSMENT (2026-07-12, prompted by #571's 187 Titan-Iapetus candidates) — the
  Iapetus/Saturn re-scope the kill itself pointed at now has a real motivating result, and the
  Triton kill-logic does NOT transfer to it. Full writeup:
  `docs/notes/2026-07-12-552-revival-assessment.md`.** Ran the SAME back-of-envelope bend check
  that killed Triton, on all 187 flagged #571 candidates, via `core/flyby.py::max_bend` +
  `core/satellites.py`'s sourced Iapetus constants (GM=120.515, r_p,min=834.3 km). Iapetus
  orbital speed 3.263 km/s → worst-case out-of-plane node penalty v·sin(15.5°)=0.872 km/s
  (confirms the caveat's ~0.85). The 5° two-sided gate ceiling is Iapetus V∞=1.780 km/s; adding
  0.872 in quadrature, **149/187 (80%) still clear the 5° bend gate** (127-182 across an 8-17°
  inclination band), real-corrected Iapetus bend 4.15-8.87° (median 5.41°), best cluster ~8.9° at
  the analytically-predicted min-achievable V∞≈0.98. This is GENUINELY DIFFERENT from Triton:
  Iapetus is prograde+inclined (a bounded plane change, √(V∞²+0.87²), ~15-40% V∞ bump) not
  retrograde (~180° flip, V∞→~8.5, bend→~1.5°, 0 survivors). The "coplanar model generous → real
  geometry kills it" logic does not apply here. **BUT** the survivorship is a max-bend
  *necessary-not-sufficient* pass, not a closure proof, and Titan/Iapetus e≈0.028 (Fable flag #4,
  ±0.16 km/s ≈ 3× residual floor) caps the upside. **Verdict: do NOT revive the general
  multi-week build; do NOT declare a family; do the cheap day-scale 3D-closure probe first — see
  #572 below.** `core/lambert.py` is already fully 3D (takes 3D r1/r2, cross_z branch +
  prograde/retrograde selector); the only coplanar assumption is upstream state generation
  (`_moon_state()` z=0), so a single-candidate hand-check is a day-scale probe, not a build.

- **#572** (P1, cheap ~1-day probe — the load-bearing gate before ANY #552 build; do NOT
  auto-fire the corrector build, this decides whether it's worth scoping) — Titan-Iapetus 3D
  closure probe on the top #571 candidates by real-corrected margin.
  **FABLE CORRECTION (2026-07-12) — the original exemplar candidate below is WRONG, do not use
  it:** the previously-cited "rel_offset 18° / phase0 0° / tof_scale 1.15 / n_rev (0,0) /
  residual 5.6e-5" record is the best-by-RESIDUAL, not best-by-corrected-margin — its actual
  Iapetus V∞ is 1.680 km/s and its real-corrected bend is 4.44°, which FAILS the 5° gate under
  the very correction this probe is meant to test. Using it would test an already-dead candidate
  and produce a garbage verdict either way. **Use these instead** (independently verified by
  Fable against the raw sweep data): **rel_offset 255° / tof_scale 1.80 / n_rev (1,1),
  Titan-anchored** (Iapetus V∞ 0.9815, corrected bend 8.86°, residual 3.3e-3) and **rel_offset
  89° / tof_scale 0.70 / n_rev (0,0), Titan-anchored** (V∞ 0.9858, corrected bend 8.83°,
  residual 6.2e-3) — both from `data/scan_571_saturn_titan_iapetus.jsonl`. **Deliberately prefer
  Titan-anchored candidates over Iapetus-anchored ones**: the 118 Iapetus-anchored survivors (two
  Iapetus encounters ~2×tof apart at different longitudes) face a much harder double-node
  phasing constraint the original assessment never accounted for — start with the Titan-anchored
  pair above, not the best Iapetus-anchored candidate (rel_offset 215° / tof_scale 2.75, corrected
  min bend 6.51°) even though it's also gate-passing.
  **Implementation notes (also Fable corrections):**
  (1) **Node alignment must be an explicit free variable in the probe**, not fixed at the
  coplanar seed's arbitrary phase — the coplanar problem is rotationally symmetric, so a
  Titan-anchored candidate can be rigidly rotated to place its Iapetus encounter at the real
  Titan-Iapetus mutual node; a negative result without searching over node alignment is
  formulation-conditional, not a real closure failure, per this project's own discipline.
  (2) **Do not count a `LambertGeometryError` (near the 180°-transfer singularity, which
  node-aligned inclined geometries can approach) as "no closure"** — that's a solver-domain
  artifact, not evidence of dynamical infeasibility; retry at a nearby offset before concluding
  no closure exists.
  (3) **Frame the deliverable explicitly as quasi-cycler-class evidence, not ballistic-cycler
  evidence** — Fable's own required-turn check (comparing the actual turn angle demanded by the
  geometry against Titan's real bend capability) confirms this whole genome, including #312's
  own verified anchor point, is quasi-cycler class (V2 fails on drift BY DESIGN, per
  `FAIL_QUASI_BOUNDED`), not a ballistic-cycler-class member — do not let a successful 3D
  closure here be read as "ballistic cycler," the standing it would earn is the same
  quasi-cycler-class standing #312's own family has.
  (4) **Do not read "80% of 187 survive the bend gate" as "149 robust candidates."** Post-
  correction median bend across all 149 survivors is 5.50° — right at the floor; the genuinely
  robust core is specifically the V∞≈0.98 Titan-anchored cluster named above (8.8-8.9° bend),
  not the full 149-candidate population.
  Place Iapetus on its real inclined orbit (i≈15.5° to Titan's plane — a conservative worst-case
  proxy verified by Fable to stay safely inside the sensitivity band even accounting for
  Iapetus's own e≈0.028 node-crossing-speed effect, true worst case ≈15.8-16.1°) in a
  throwaway/one-off state generator, feed the real 3D positions to `core/lambert.py` (Fable
  independently confirmed this is genuinely already 3D-capable: full 3D r1/r2, transfer angle
  from the 3D dot product, `cross_z`-based short/long-way branch selection with an explicit
  prograde flag, 3D-correct f,g velocity construction — the ~1-day cost estimate holds), and
  attempt an actual 3D Lambert closure near each coplanar seed, searching over node alignment
  per note (1). **Decision gate:** does a real 3D solution exist with residual near the coplanar
  value AND every encounter still clearing the #324 bend gate? This is the closure question the
  #571 assessment's bend check (necessary-not-sufficient) cannot answer. **If YES →** scope the
  narrow Titan-Iapetus 3D corrector (#552's own re-scoped ~3-5 day estimate reusing
  `core/lambert.py`'s 3D support, NOT the general n-body build) behind Opus adjudication + Fable
  second-opinion + V-gauntlet, eccentricity/duty-cycle caveat honored. **If NO →** stamp
  Titan-Iapetus as bend-feasible-but-no-3D-closure in `data/empty_regions.jsonl` and close the
  pair (clean cheap negative on a NOVEL search, per [[project_negative_results_registry]]).
  **Explicitly out of scope:** the full 3D genome capability, any catalogue.yaml edit, the
  V-gauntlet itself. **Recommended models:** Sonnet for the throwaway 3D state-gen + Lambert
  closure probe (mechanical, behind the existing #324 gate). Opus + Fable for the go/no-go
  adjudication on whether a genuine 3D closure survived.
  **RESULT (2026-07-12) — GENUINE 3D CLOSURE FOUND for BOTH candidates; hand off to Opus+Fable
  adjudication next, not adjudicated here.** Built the throwaway 3D state generator
  `scripts/probe_572_titan_iapetus_3d_closure.py` (Titan kept coplanar/equatorial per #571's
  own convention; Iapetus placed on a real circular orbit inclined 15.5 deg via the standard
  R3(Omega).R1(inc) orbital-plane rotation, ascending-node longitude `Omega` an explicit free
  variable). Inline smoke test confirms the generator reduces EXACTLY to `_moon_state`'s
  coplanar formula at inc=0 for any Omega (caught + fixed a rad/day-vs-rad/s velocity-unit bug
  this way before the real search ran). **Node-alignment search:** fine 3600-point (0.1 deg)
  grid over Omega at each candidate's own (rel_offset, tof_scale, n_rev), all local-minimum
  basins in the residual-vs-Omega landscape enumerated (not just the single deepest one — a
  first pass that seeded Nelder-Mead only from the global residual minimum landed on a
  physically different high-V_inf branch, ~8.6 km/s, that fails the bend gate; the genuine
  low-V_inf closure sits in a DIFFERENT, shallower-but-still-sub-gate basin the naive search
  missed), each basin locally refined via a bounds-restricted Nelder-Mead (Omega +/-15 deg,
  tof_scale +/-0.1, preventing branch-jumping) over (Omega, tof_scale) jointly. No
  `LambertGeometryError` was raised at any of the 3600 Omega samples for either candidate (the
  180-deg-singularity retry path was implemented but never exercised). **Cand 1** (rel_offset
  255 deg, tof_scale 1.80, n_rev=(1,1)): closure at Omega=30.91 deg, tof_scale=1.793, residual
  1.7e-9 km/s, V_inf [Titan-out 1.870, Iapetus 1.207, Titan-in 1.870] km/s, bends [45.47, 10.35,
  45.47] deg (Iapetus encounter, the tight one, clears the 5 deg floor by >2x) — physical gate
  PASS. **Cand 2** (rel_offset 89 deg, tof_scale 0.70, n_rev=(0,0)): closure at Omega=37.46 deg,
  tof_scale=0.688, residual 2.0e-9 km/s, V_inf [2.307, 1.203, 2.307] km/s, bends [34.04, 10.41,
  34.04] deg — physical gate PASS. Both realized V_inf sit ~20-35% above the coplanar seed
  (consistent with the OUTSTANDING.md #552 revival-assessment's own analytic quadrature-sum
  prediction, `sqrt(V_inf_coplanar^2 + 0.87^2)`, for the out-of-plane node penalty). Both
  independently cross-checked with `saturn_uranus_campaign.dop853_cross_check_leg` (numerical
  integration, NOT the analytic Lambert solve itself): arrival-position error 1.1e-5 to
  1.2e-4 km, far under this project's standing <1 km cross-check floor — the closures are
  genuine converged two-body conic arcs, not solver artifacts. **Decision-gate criterion used:**
  residual < GATE_RESIDUAL_KMS (0.05 km/s, the same bar #558/#571 use everywhere) counts as
  "near the coplanar residual"; both closing basins land at ~1e-9 km/s, far inside it. **Did NOT
  touch `data/empty_regions.jsonl`** (only written on a NO-closure verdict per this task's own
  scope) and did NOT edit `catalogue.yaml` or run the V-gauntlet, per scope. **Framing
  (mandatory, per this entry's own note (3) above): this is quasi-cycler-class evidence, NOT a
  "ballistic cycler" finding** — same standing as #312's own family, V2 fails on drift by
  design. **This is a factual closure finding only — full Opus + Fable go/no-go adjudication on
  whether to scope the #552 narrow Titan-Iapetus 3D corrector is the next step, not performed
  here.** Artifacts: `scripts/probe_572_titan_iapetus_3d_closure.py`,
  `data/probe_572_titan_iapetus_3d_closure.jsonl` (full basin-by-basin record for both
  candidates). Ratchet suite not required (empty_regions.jsonl untouched); ruff clean.
  **ADJUDICATION (2026-07-12, Opus) — VERDICT: WIDEN THE PROBE FIRST (new task #573),
  do NOT jump to the corrector build; do NOT stamp the pair empty (the closures are real).**
  Read the full basin-by-basin JSONL for both candidates directly; the closures are genuine
  (residual ~1.7e-9/2.0e-9 km/s, DOP853 cross-check 1e-5–1e-4 km, the low-V∞ gate-passing basin
  is a genuinely distinct basin from the high-V∞ ~8.6 km/s branch that fails the bend gate). **But
  2-for-2 is non-representative and cannot justify a 3–5 day corrector build:** (a) the two tested
  candidates are the *best-by-corrected-margin* robust core (V∞≈0.98 cluster) — this task's own
  note (4) says the 149-survivor **median** bend is 5.50°, right at the floor, so 2-for-2 on the
  two most robust candidates says almost nothing about the median candidate; (b) machine-precision
  closure is *generically* achievable here — the probe drives a scalar residual to ~0 with two free
  knobs (Ω, tof_scale), so ~1e-9 residual is not evidence of rarity; the discriminating fact is
  whether the closing basin lands in a bend-gate-PASSING branch, which for floor-hugging candidates
  is unknown. **Factual correction to carry forward:** the binding *Iapetus* (flyby-body) bends are
  **10.35° and 10.41°** (2.07× the 5° floor), NOT the 45.5°/34.0° cited in some summaries — those
  are the *Titan* bends. **The widened probe is nearly free:** #572 cost 1.00 s + 0.25 s wall-clock
  per candidate; the ~69 Titan-anchored coplanar-gate-passing candidates in
  `data/scan_571_saturn_titan_iapetus.jsonl` (all `physical_gate_passed=True` records are
  Titan-anchored; n_rev 26×(0,0)/40×(1,1)/3×(2,2)) run in minutes, well under the 10-min cap — a
  single synchronous run converts the selection-biased 2-for-2 into a real population hit-rate and
  a #312-scale-vs-isolated-pair family-size estimate before any 3–5 day investment. **Still
  genuinely unknown (honest):** the probe is *circular* (Titan+Iapetus e≈0.028, 7–25× the Uranian
  e≤0.004) — it says NOTHING about real-ephemeris survival; the idealized→real duty-cycle gap
  should be materially worse than #568's 61.7–79.1%, and is unknowable until the corrector + a
  real-ephemeris Saturn (SPICE SAT441 kernel, fetchable, not on host; `verify/mission_spk.py`
  already maps Saturn) validation pipeline exist. Widening does not resolve eccentricity either,
  but if the population hit-rate is low we stop before ever building that pipeline — correct
  sequencing. Full reasoning + concrete corrector deliverable shape in
  `docs/notes/2026-07-12-572-closure-adjudication.md`. **Framing (mandatory): quasi-cycler-class
  evidence, same standing as #312's own family (V2 fails on drift by design), NOT a ballistic
  cycler and NOT a novelty claim.** Next step: **#573** (widened 3D-closure probe), which must pass
  a Fable second-opinion on its plan before dispatch.

- **#573** (P1, cheap ~minutes-scale probe — the population-rate gate before ANY #552 corrector
  build; plan under Fable review before dispatch, do NOT auto-fire) — widened Titan-Iapetus
  3D-closure probe. **Motivation:** #572 found genuine 3D closures for BOTH tested candidates, but
  those two were the *best-by-corrected-margin* robust core (V∞≈0.98 cluster, Iapetus bend
  8.8–8.9° coplanar), and #572's own note (4) puts the 149-survivor **median** bend at 5.50° —
  right at the 5° floor. 2-for-2 on the two most robust candidates does not estimate the closure
  rate of the *population*, and machine-precision closure is generically achievable with the
  probe's two free knobs (Ω, tof_scale), so the discriminating fact is not "residual→0" but
  "does the closing basin PASS the bend gate." A 3–5 day corrector justified on 2 candidates risks
  discovering the family is 2–3 members, not #312-scale. **Task:** generalize
  `scripts/probe_572_titan_iapetus_3d_closure.py` (which is already a clean throwaway generator
  with a validated inc=0 smoke test) to run the SAME 3D node-alignment closure probe over the full
  Titan-anchored candidate population, and report the population-level gate-passing closure rate +
  the distribution of realized Iapetus (flyby-body) bends across the closers. **Input population
  (pin precisely — Fable flag):** the ~69 coplanar-gate-passing (`physical_gate_passed=True`)
  Titan-anchored candidates in `data/scan_571_saturn_titan_iapetus.jsonl` (n_rev split
  26×(0,0)/40×(1,1)/3×(2,2)); Fable should confirm whether the inclination-corrected-but-coplanar-
  fail candidates from the #552 revival assessment's "187 flagged / 149 survive" set also belong,
  or whether the coplanar-gate-pass set is the right input. **Keep it circular-inclined**
  (apples-to-apples with #572; Titan coplanar-circular, Iapetus inclined-circular at i≈15.5°) —
  eccentricity is deferred to the corrector/validation stage, NOT this probe (folding an eccentric
  variant in is defensible but treated as scope creep for the population-rate question; Fable may
  reargue). **Cost:** #572 cost 1.00 s + 0.25 s/candidate, so the full set runs in minutes — MUST
  complete synchronously in a single sub-10-minute run (per this project's hard no-backgrounding
  rule); if the Iapetus-anchored candidates (harder double-node phasing, deferred by #572) are
  added as an optional second phase, they get their own synchronous run. **Decision gate:** a
  *substantial* gate-passing closing population (a real family, not 2–3 isolated points) → scope
  the narrow Titan-Iapetus 3D corrector (#552's re-scoped ~3–5 day estimate reusing
  `core/lambert.py`'s confirmed 3D support, NOT the general n-body build) as a subsequent task
  behind Opus adjudication + Fable second-opinion + V-gauntlet, eccentricity/duty-cycle caveat
  front-and-center (concrete deliverable shape in
  `docs/notes/2026-07-12-572-closure-adjudication.md` §6); only the robust-core cluster closes
  (a handful) → re-adjudicate whether a 2–3 member idealized quasi-cycler set justifies the build;
  broad closure but Iapetus bends hugging the 5° floor → flag that the eccentricity gap (§5) will
  likely dominate before committing. **Explicitly out of scope:** the corrector build itself, the
  full 3D genome capability, any `catalogue.yaml` / `empty_regions.jsonl` edit, the V-gauntlet, the
  SPICE SAT441 kernel fetch. **Framing (mandatory, per #572 note 3): quasi-cycler-class evidence,
  same standing as #312's own family (V2 fails on drift by design), NOT a ballistic cycler, NOT a
  novelty claim** (internally-enumerated fact about our own search space, per `our_status`
  discipline). **Recommended models:** Sonnet for the mechanical widened probe behind the existing
  #324 gate; Opus + Fable for the go/no-go population adjudication. **REQUIRED: a Fable
  second-opinion pass on the #573 plan before dispatch** — this thread's Fable discipline has
  caught real load-bearing issues at essentially every stage (including the wrong-exemplar-candidate
  bug in #572's own plan and the max_bend triage inversion in #565); do not dispatch #573 without
  it.
  **FABLE PLAN REVIEW (2026-07-12): CONFIRMED WITH CORRECTIONS — one factual error, one
  load-bearing grid-aliasing gap (same failure class as #562), several precision fixes. All
  folded into scope below; dispatch only with these applied.**
  **(1) Input population — RESOLVED, use all 69 Titan-anchored gate-passers, no pre-filtering.**
  The "187/149 set" framing was based on a FACTUAL ERROR in
  `docs/notes/2026-07-12-572-closure-adjudication.md` §4 (it claimed some candidates fail the
  coplanar gate but pass after correction — impossible; the worst-case quadrature correction only
  INCREASES Iapetus V∞, so 149 is a strict subset of the 187, never an alternative). Correct that
  note. The only genuinely additional population beyond the 69 Titan-anchored is the **118
  Iapetus-anchored** candidates (already scoped as an optional phase 2) — make phase 2 MANDATORY
  if the Titan-anchored result lands in the "marginal" bucket (3), and note it requires a real
  generator change (Iapetus becomes the ANCHOR with double-node phasing, not a parameter tweak).
  **(2) Circular-inclined scoping is sound as a population-RATE estimate** (deferring eccentricity
  doesn't systematically bias the rate — the modulation is phase-symmetric across a spread
  population) but add one free, cheap analytic post-processing step: **report the count of closers
  whose realized Iapetus bend is ≥6.0°** (the bend at V∞=1.62 km/s, i.e. the margin surviving a
  worst-case +0.16 km/s eccentric shift, verified against sourced GM=120.515/r_p=834.3) as an
  "eccentricity-robust" sub-count — makes bucket (3) below mechanical instead of vague.
  **(3) LOAD-BEARING: the 69 records are grid samples along ~10 CONTINUOUS solution branches, not
  69 independent candidates — dedup by branch before applying ANY count threshold, or the
  decision gate is meaningless.** Verified directly: contiguous rel_offset runs (e.g. 18
  consecutive 1°-spaced records at one (n_rev, tof_scale) combo) chain across tof_scale steps at a
  consistent drift, wrapping at the 360° seam — exactly the #562 grid-aliasing failure mode.
  WORSE: the probe's Nelder-Mead refinement window (±0.1 tof_scale) is WIDER than the input grid
  spacing (0.05), so neighboring grid records will silently refine to the IDENTICAL closure and
  get double-counted as if independent. **Required fix**: cluster closers by grid adjacency in
  (n_rev, tof_scale, rel_offset) with 360° wraparound AND/OR dedup refined solutions by proximity
  in (Ω*, tof_scale*, V∞); state the decision gate in DISTINCT BRANCHES, not raw record counts.
  **Concrete branch-deduped thresholds** (calibrated against the #558/#563 Uranian precedent, ~5
  closures/pair, and against the known robust core = exactly 2 branches):
  - **≥5 distinct gate-passing branches, ≥2 n_rev classes, ≥3 eccentricity-robust (bend ≥6.0°)**
    → real family; scope the corrector per the adjudication note's §6 deliverable shape.
  - **≤2 closing branches** (i.e. only the already-#572-tested robust core) → adds nothing new;
    re-adjudicate whether 2 members alone justify the build (likely no).
  - **3-4 branches, OR ≥5 but <3 eccentricity-robust** → marginal: mandatory phase 2
    (Iapetus-anchored) before adjudicating, and/or flag the eccentricity gap as likely dominant.
  **(4) Cost estimate arithmetic was WRONG but the conclusion holds**: "1.00s + 0.25s/candidate"
  misread #572's two candidates' INDIVIDUAL wall-clocks as a base+marginal formula. Correct:
  69 × ≤1s ≈ ≤70s core, comfortably sub-10-minutes even with a 5-10x basin-count blowup on
  pathological candidates. Keep incremental per-candidate JSONL checkpointing regardless (stall
  visibility + partial-result survival).
  **(5) Four at-scale methodology guards required in the generalized probe script** (the #572
  machinery was hand-verified on 2 candidates only, not battle-tested at 69×3600 samples):
  (a) cap basin refinement at ~12 deepest basins/candidate, log when the cap binds; (b) set the
  Nelder-Mead window to min(15°, half-distance to the nearest adjacent grid-minima) and
  post-verify the refined solution stayed within its seed basin's grid support, flagging (not
  silently counting) any that didn't; (c) record NM runs that pin at a box edge as
  formulation-conditional non-closures, distinct from genuine non-closures, per this project's
  own formulation-conditional discipline; (d) track `LambertGeometryError` occurrences across the
  full grid (untested at this scale in #572) so a solver-domain NaN band can't silently split or
  hide a basin — and apply the standing sweep-singleton-artifact rule
  ([[feedback_isolated_sweep_flips_suspect_artifact]]) to any isolated closure flip inside an
  otherwise-contiguous branch: investigate, don't count at face value.
  **✓ Resolved (2026-07-12).** Generalized `scripts/probe_572_titan_iapetus_3d_closure.py` into
  `scripts/run_573_titan_iapetus_population_closure.py` (reuses the validated #572 core --
  `iapetus_state_3d`, `_leg_best`, the rotation algebra, the inc=0 smoke test, unmodified -- only
  the dedup/guard logic is new) and ran it as ONE synchronous foreground call over all 69
  Titan-anchored `physical_gate_passed=True` candidates from `data/scan_571_saturn_titan_iapetus.jsonl`
  (n_rev split 26x(0,0)/40x(1,1)/3x(2,2), matching the #572 adjudication's count). **Runtime: 42.2 s
  total** (0.26-0.94 s/candidate), far under the 10-minute cap; results incrementally checkpointed
  (append+flush) to `data/probe_573_titan_iapetus_population_closure.jsonl`.
  **Branch dedup (item 1):** implemented as union-find proximity clustering of every genuine closing
  basin found across all 69 runs, merging two closures iff same n_rev AND Omega* within 5 deg
  (periodic) AND tof_scale* within 0.05 AND V_inf(mid) within 0.10 km/s. Verified this collapses the
  known-adjacent grid runs the Fable review flagged directly: the 18 candidates at
  (n_rev=(1,1), tof_scale=2.5, rel_offset=156-173) all found the SAME shallow non-gate-passing basin
  and correctly produced zero duplicated closures; separately, 8 genuinely-adjacent contiguous runs
  (rel338-340, rel88-90, rel319-320, rel183-185, rel57-58, rel121-122, rel37-38, rel228-229) each
  correctly merged into ONE multi-member branch (2-3 members apiece) rather than being double-counted.
  Two branches (ids 2/19 and 4/14) are exact z-mirror-reflection duplicates (Omega differing by
  ~180.00 deg, IDENTICAL V_inf/bend magnitudes -- confirmed analytically: the 3D rotation algebra is
  invariant under (u -> u+180, Omega -> Omega+180) up to a sign flip on z, so these are two genuinely
  distinct out-of-plane (above/below-primary-plane) 3D geometries, not aliasing artifacts) --
  reported as-is (22), consistent with the project's own #563 precedent of counting symmetric
  partners as distinct family members, not collapsed further.
  **Result: 22 DISTINCT BRANCHES** (33 raw closing-basin instances from 28/69 candidates, deduped to
  22), spanning **all 3 n_rev classes** ((0,0): 12 branches, (1,1): 7 branches, (2,2): 3 branches),
  of which **17/22 are eccentricity-robust** (Iapetus/flyby-body bend, the middle of the 3-element
  `max_bend_deg_per_encounter` array, >= 6.0 deg -- the margin surviving a worst-case +0.16 km/s
  eccentric perturbation); the 5 non-robust branches all hug the 5.0-5.4 deg floor.
  **DECISION BUCKET: REAL_FAMILY_SIGNAL** (22 >= 5 distinct branches, 3 >= 2 n_rev classes,
  17 >= 3 eccentricity-robust -- clears the threshold by a wide margin on all three axes).
  **Phase 2 (118 Iapetus-anchored candidates) is NOT mandated** by this result (only required for
  the marginal bucket); still a defensible, cheap follow-up but not executed here per scope.
  **Load-bearing caveat, found during this run and worth carrying forward:** the #572-confirmed
  robust-core "cand1" (rel_offset 254/255/256/257, tof_scale=1.8, n_rev=(1,1) -- #572's
  *highest-margin* tested candidate, genuine closure at residual 1.7e-9 km/s under #572's own fixed
  +-15 deg window) surfaces in THIS run as a **formulation-conditional non-closure** (boundary-pinned)
  rather than a counted genuine closure, for all 4 grid points in that group (3 pinned basins each).
  Root cause verified directly: the true 2D (Omega, tof_scale) joint minimum (Omega=30.91) sits in a
  33.6 deg gap between two raw-grid-detected basins (16.5 deg and 50.1 deg), and the item-2 adaptive
  window is correctly narrowed by a WEAK unrelated third basin 17.6 deg away, so neither seed's window
  reaches the true minimum -- exactly the box-edge-pin scenario item 4 anticipates. This means the
  22-branch count is a **conservative undercount** (a known-genuine closure is missing from it, not a
  false positive padding it), which if anything strengthens the REAL_FAMILY_SIGNAL verdict. The
  OTHER #572 candidate ("cand2", rel_offset=88-90, tof_scale=0.7, n_rev=(0,0)) IS correctly recovered
  here as genuine closure branch id 6 (Omega*=37.36, V_inf=[2.34,1.21,2.34], bend=10.38 deg, 3
  members merged), matching #572's own reported values closely.
  **Guards triggered:** 0/69 candidates hit the 12-basin refinement cap (max seen: 9 basins, at
  several n_rev=(1,1)/(2,2) candidates); 20 basins total flagged formulation-conditional
  (boundary-pinned), concentrated at the #572 cand1 group (10 of the 20) plus scattered singles
  elsewhere -- reported, not silently dropped, not counted as closures; 0 `LambertGeometryError` and
  0 `LambertConvergenceError` raised across the full grid (69 x 3600 Omega samples x 2 legs each,
  ~497k Lambert solves) -- the retry path exists but was never exercised, consistent with #572's own
  finding at 2 candidates; 0 isolated singleton closure flips flagged by
  `detect_isolated_singleton_anomalies` across all contiguous (n_rev, tof_scale) rel_offset runs of
  length >=3 (the 18-point rel156-173 run and all others were closure-boolean-consistent throughout).
  **Framing (mandatory): quasi-cycler-class evidence, same standing as #312's own family (V2 fails on
  drift by design), NOT a ballistic cycler, NOT a novelty claim** -- this is an internal fact about
  our own idealized circular-inclined search space, not a literature/novelty claim; the
  eccentricity/duty-cycle gap to real ephemeris (Titan/Iapetus e~=0.028, vs the Uranian moons'
  e<=0.004) remains completely unaddressed by this probe, per the #572 adjudication note's own
  "still genuinely unknown" caveat. **Did NOT** run the 118 Iapetus-anchored candidates, scope/build
  the 3D corrector, edit `catalogue.yaml`/`empty_regions.jsonl`, or run the V-gauntlet, per scope.
  **This is a factual population-rate finding only** -- full Opus + Fable go/no-go adjudication on
  whether to scope the #552-rescoped 3-5 day Titan-Iapetus 3D corrector is the next step, not
  performed here. Artifacts: `scripts/run_573_titan_iapetus_population_closure.py`,
  `data/probe_573_titan_iapetus_population_closure.jsonl` (full basin-by-basin record, all 69
  candidates + population summary). `uv run pytest tests/data/test_outstanding_structure.py -q`
  clean; ruff clean; `data/empty_regions.jsonl`/`catalogue.yaml` untouched so the full data/search
  ratchet suite is not required.
  **CORRECTOR GO/NO-GO (2026-07-12, Opus) — VERDICT: GO. Scope+greenlight the narrow
  Titan-Iapetus 3D corrector as new task #574, STAGED (eccentric-Keplerian kill-gate first, then
  SPICE validation only if it passes). Full writeup:
  `docs/notes/2026-07-12-573-corrector-gono-go.md`.** Read all 22 branch records directly from
  `data/probe_573_titan_iapetus_population_closure.jsonl`: Iapetus V∞ 1.19-1.78 km/s (median 1.51,
  a tight quadrature-corrected band exactly matching the #552 revival-assessment's analytic
  `sqrt(V∞_coplanar²+0.87²)` prediction), Iapetus bend 5.02-10.57° (median 6.81°), all 3 n_rev
  classes ((0,0)×12/(1,1)×7/(2,2)×3), 18 branches at machine precision (~1e-9) + 4 looser (2e-6 to
  3e-3, all far inside the 0.05 gate). **REAL_FAMILY_SIGNAL is genuine, NOT the #558 artifact:**
  (a) #558's "1000+ basins" was a raw pre-dedup grid count; #573's raw is only 33 closing-basin
  instances → 22 deduped via the Fable-mandated union-find clustering that fixes the exact #558/#562
  aliasing mode — the anti-aliasing discipline #558 lacked is already built in; (b) the method here
  ERRS toward UNDERcounting (#572's own highest-margin confirmed cand1 closure is MISSING from the
  22 as a boundary-pinned formulation-conditional non-closure — 22 is a conservative floor, opposite
  to inflation), backed by 0 geometry/convergence errors across ~497k solves + 0 singleton flips;
  (c) the decision is robust to count uncertainty (the node-alignment Ω is an extra DOF vs the
  coplanar Uranian search, so "22 ≫ ~5/pair" partly reflects the added knob — but 22/3/17 clears the
  pre-registered ≥5/≥2-classes/≥3-robust gate by ~4× regardless, and even collapsing both z-mirror
  pairs leaves 20). **The one genuinely-unresolved risk is eccentricity, unchanged:** the whole
  #571→#573 stack is CIRCULAR (`core/satellites.py` has no e field; "ecc-robust" is only a
  bend-margin proxy, not an eccentric-orbit closure test), and Titan/Iapetus e≈0.028 is 7-25× the
  Uranian e≤0.004. **Honest real-ephemeris expectation:** duty cycle MATERIALLY WORSE than the
  Uranian 61.7-79.1% — plausibly ~30-55% for survivors with a real chance broad swaths fail
  V4-strict and the family thins to a handful (the two 1.19-1.21 km/s / 10.4-10.6°-bend branches
  ids 6,7 are the likeliest survivors; the five bend≤5.4° floor-huggers ids 11,15,17,18,21 the
  likeliest to die). **"Worse duty cycle, still real" — or a clean documented negative — is an
  acceptable outcome worth finding out; lower duty cycle is a column, not a disqualifier** (Uranian
  earned a writeback at 61.7-79.1%). The cheap circular probes are now EXHAUSTED — only building the
  eccentric stage answers the question. **Why GO rather than another probe:** #573 cleared a
  PRE-REGISTERED, Fable-reviewed gate (with the ecc-robust bend proxy folded in) — adding a fifth
  standalone pre-corrector adjudication would be goalpost-shifting; the right home for the remaining
  eccentricity check is INSIDE #574 as its Stage-A kill-gate. **Framing (mandatory):
  quasi-cycler-class evidence, same standing as #312's own family (V2 fails on drift by design), NOT
  a ballistic cycler, NOT a novelty claim.** Did NOT build the corrector, edit `catalogue.yaml`,
  fetch any kernel, or run the V-gauntlet, per scope. Next step: **#574**, which MUST pass a Fable
  second-opinion on its spec before dispatch.

- **#574** (P1, the #573-cleared narrow Titan-Iapetus 3D corrector build — STAGED; **REQUIRES a
  Fable second-opinion on this spec before dispatch**, per this thread's standing discipline that
  has caught a real load-bearing issue at literally every stage: #552's payoff claim, #562's
  grid-aliasing, #565's bend inversion, #571's two-sided-gate emptiness, #572's wrong exemplar
  candidate, #573's own grid-aliasing risk) — productize the #573 idealized-circular 3D
  closure engine into a narrow, repeatable Titan-Iapetus 3D corrector behind the #324 gate AND
  attach a real-ephemeris Saturn validation, mirroring the Uranian V1-V4-strict chain (#566/#568).
  **This is the narrow #552 re-scope (~3-5 days reusing `core/lambert.py`'s confirmed 3D support),
  NOT the general n-body/arbitrary-inclination capability killed on 2026-07-10.** **Motivation:**
  #573 cleared REAL_FAMILY_SIGNAL (22 distinct gate-passing branches, 3 n_rev classes, 17
  eccentricity-robust) with margin; the only load-bearing unresolved risk is eccentricity
  (Titan/Iapetus e≈0.028, 7-25× the Uranian e≤0.004) which the entire circular #571→#573 stack
  cannot address. **STAGE IT — Stage A gates Stage B:**
  **FABLE PLAN REVIEW (2026-07-12): CONFIRMED WITH CORRECTIONS — input set (22 branches, 17
  ecc-robust, floor-hugger ids 11/15/17/18/21 at bends 5.28/5.37/5.02/5.32/5.02°) independently
  verified against `data/probe_573_titan_iapetus_population_closure.jsonl`, matches exactly, no
  #572-style wrong-input bug. But three load-bearing corrections (C1-C3) are MANDATORY before
  dispatch, plus three strongly-recommended fixes (C4-C6). All folded into the corrected Stage A
  text below (this is the SEVENTH real catch on this thread, streak intact).**
  **Stage A — eccentric-Keplerian closure kill-gate (cheap, ~part-day, NO SPICE, NO kernel fetch).**
  Extend the #573 circular state generators (`iapetus_state_3d` + the coplanar Titan state) to place
  BOTH moons on real *eccentric* Keplerian orbits: add e≈0.0288 (Titan) / 0.028 (Iapetus) — these
  are NEW data, `core/satellites.py` currently has no eccentricity field, add it sourced (JPL SSD)
  or in a throwaway local constant with a sourced comment.
  **(C1 — MANDATORY, prevents a known project failure mode):** the free parameters are **two mean
  anomalies at epoch, M0_Titan and M0_Iapetus** (a 4D per-branch space: Ω, tof_scale, M0_T, M0_I) —
  NOT a free "periapsis phase at the encounter." Every encounter's true anomaly MUST be derived by
  Kepler-propagating from these two epoch M0 values over the actual leg TOFs (mirroring exactly how
  `scripts/run_573_titan_iapetus_population_closure.py::evaluate_point_tracked` derives Titan's
  SECOND state at t=2·tof by mean-motion propagation from its state at t=0, not a free
  re-specification). A free per-encounter phase is precisely the #480 EGGIE per-encounter
  self-consistency bug ([[feedback_constructed_tour_per_encounter_self_consistency]]) — do not
  reintroduce it here.
  **(C2 — MANDATORY, positive control):** before crediting ANY eccentric kill verdict, the eccentric
  machinery MUST be shown to recover the 22 circular branches (at their circular residual values)
  when run at e=0 — mirror #573's own `smoke_test_reduction_pass` meta-record pattern. An unverified
  eccentric generator producing an all-negative result is indistinguishable from a real family death
  and a silent implementation bug; per this project's own [[feedback_verify_gauntlet_with_positive_control]]
  rule, do not trust the kill without this check passing first.
  **(C3 — MANDATORY, precise kill-gate + dedup discipline):** survivors are counted the SAME way
  #573 counted branches — eccentric closure at residual ≤0.05 km/s + full #324 bend-gate pass,
  THEN union-find proximity dedup EXTENDED to include the (M0_T, M0_I) coordinates (the near-mirror
  branch pairs {2,19} and {4,14} are known merge candidates — an undeduped eccentric count would
  silently regress to raw-count thinking, the exact #562/#573-planning failure mode). Pre-registered
  thresholds: **PASS** = ≥5 deduped survivors from the 17 ecc-robust input, spanning ≥2 n_rev
  classes (same ≥5 anchor as the Uranian/#573 precedent). **KILL** = ≤3 deduped survivors — stop,
  document a conditional negative, do not proceed to Stage B. **4 survivors = MARGINAL, requires an
  explicit Opus adjudication, not an automatic proceed either way.** Floor-hugger controls are
  expected to die; if ≥3 of the 5 control branches instead SURVIVE, the eccentricity-robust ≥6.0°
  proxy is non-discriminating at this stage — flag explicitly and adjudicate rather than silently
  trusting the main-population result. **Recommended methodology (not mandatory, but preferred over
  a fresh grid sweep):** continuation in eccentricity — step e from 0→0.0288/0.028 refining
  (Ω, tof_scale, M0_T, M0_I) from each already-known circular branch's solution — rather than a new
  higher-dimensional grid, since continuation avoids both aliasing directions at once (no false
  kills from a real basin falling between grid points, no double-counting) and yields the e=0
  control for free; if a grid is used anyway, require refinement seeded from each circular branch
  plus the standing singleton-anomaly guard ([[feedback_isolated_sweep_flips_suspect_artifact]]).
  Re-run the #573 3D closure + #324 gate on the **17 eccentricity-robust branches** (the 5
  floor-huggers ids 11/15/17/18/21, verified bends 5.28/5.37/5.02/5.32/5.02°, are the
  expected-to-die control per C3 above). **Kill gate: per C3's precise thresholds above** — STOP
  with a documented conditional negative in `docs/notes/` (and optionally an `empty_regions.jsonl`
  `reverification` entry) on KILL, do NOT fetch SAT441 or build Stage B; a clean cheap negative on a
  novel search is a fine stop ([[project_negative_results_registry]]). Must complete synchronously
  (minutes-to-part-day, per the no-backgrounding rule); checkpoint per-branch (append+flush).
  **(C4 — recommended) SPICE verification method, named precisely for Stage B:** mirror
  `src/cyclerfinder/verify/spice_kernels.py`'s existing `ensure_jup365_kernel()`/
  `ensure_pluto_kernel()` pattern (the #510/#550 PC(3,2) precedent) — verify pre-fetch via the NAIF
  `generic_kernels/spk/satellites/` directory listing + `aa_summaries.txt` (do not assume the kernel
  is named `sat441.bsp`, a successor may exist), and post-fetch via `spkcov` for NAIF IDs 606
  (Titan) and 608 (Iapetus) over the target epoch window.
  **(C5 — recommended) duty-cycle prior framing:** the "~30-55%" figure in Stage B below is a rough
  DIRECTIONAL extrapolation from the Uranian family's eccentricity ratio, NOT a computed prediction
  for Titan-Iapetus specifically and NOT a validation target — carry this caveat explicitly into any
  write-up, do not let it harden into an expected/target number.
  **(C6 — recommended) A→B transition gate:** "Stage A passes" must not be a rubber-stamp straight
  into the SAT441 fetch. Require Stage A's results to be written up in `docs/notes/` (per-branch
  survivor table, dedup applied, e=0 control result, precise threshold verdict per C3) PLUS one
  Fable pass specifically on whether Stage A was executed AS PRE-REGISTERED (not a full
  re-adjudication of the go/no-go decision, which stays settled) before the kernel fetch — escalate
  to Opus only in the C3 marginal (4-survivor) bucket. Cheap insurance immediately before the
  expensive stage, on a thread that is now 7-for-7 on real Fable catches.
  **Stage B — productized corrector + real-ephemeris SPICE validation (only if Stage A passes).**
  (1) Productize the #573 engine as a repeatable corrector behind the #324 gate (clean interface,
  provenance, stable input contract — not a throwaway sweep). (2) Add a real-ephemeris SPICE leg:
  `verify/mission_spk.py` maps `"Saturn" → "SATURN BARYCENTER"` (the barycenter, not the moon
  bodies), so this needs the **SAT441** Saturn-satellite ephemeris (or the equivalent current
  kernel — VERIFY the exact kernel name + that its coverage spans BOTH Titan and Iapetus over the
  needed epoch range BEFORE fetching; no `.bsp` is on host, only `naif0012.tls`; same NAIF
  dependency #506/#550 identified for PC(3,2) V2→V3 — the fetch is scoped into THIS task, do not
  fetch it during the go/no-go). (3) Run the V1-V4-strict gauntlet with the **#560 robustness fixes**
  and **#568 duty-cycle / synodic-boundary framing** (not raw pass%/flip%), honoring the
  eccentricity caveat throughout. (4) Behind the standard **Opus adjudication + Fable second-opinion
  + V-gauntlet** before ANY `catalogue.yaml` writeback. **Honest expectation to carry in:**
  real-ephemeris duty cycle materially worse than the Uranian 61.7-79.1% (plausibly ~30-55%,
  possibly broad V4-strict failure) — a lower-but-nonzero duty cycle is a legitimate
  quasi-cycler-class characterization (a column, not a disqualifier), and family death on
  eccentricity is an acceptable clean negative. **Framing (mandatory): any positive result is
  quasi-cycler-class evidence, same standing as #312's own family (V2 fails on drift by design,
  `FAIL_QUASI_BOUNDED`), NOT a ballistic-cycler finding and NOT a novelty claim** (internal fact
  about our own enumerated search space, per `our_status` discipline). **Explicitly out of scope:**
  the general 3D genome capability, Phase-2 Iapetus-anchored candidates (not mandated by #573's
  REAL_FAMILY_SIGNAL bucket; a defensible cheap follow-up but not part of this build). **Recommended
  models:** Sonnet for Stage A (mechanical, reuses validated #573 machinery behind the #324 gate)
  and for Stage B implementation behind the V-gauntlet; Opus + Fable for the pre-writeback
  adjudication. Full reasoning + staging rationale in
  `docs/notes/2026-07-12-573-corrector-gono-go.md`.
  **STAGE A RESULT (2026-07-12, Sonnet) — VERDICT: PASS.** Ran
  `scripts/run_574_titan_iapetus_eccentric_kill_gate.py` (continuation-in-eccentricity per the
  spec's recommended methodology, NOT a fresh grid — all 22 #573 circular branches stepped e
  0->real over 6 stages, refining the C1-mandated 4D free-parameter space (Omega, tof_scale,
  M0_Titan, M0_Iapetus) at each stage; both moons' every later state Kepler-propagated by mean
  motion from these epoch M0 values, never a free per-encounter phase). **C2 positive control:
  PASS** — the new eccentric Kepler propagator (`kepler_state_3d`) reduces exactly to
  `_moon_state`/`iapetus_state_3d` at e=0 (dr<1e-6 km, dv<1e-9 km/s, checked at a grid of
  M0/Omega/u test points), AND all 22 branches' e=0 continuation stage independently reproduces
  their #573 circular residual (0/22 smoke failures). **Main population (17 ecc-robust
  branches): 15 deduped survivors** (residual <=0.05 km/s + #324 bend-gate pass at the REAL
  eccentric V_inf; 0 merges among survivors under the C3-extended (M0_Titan, M0_Iapetus) dedup —
  every surviving branch stayed its own distinct point), spanning **all 3 n_rev classes**
  ((0,0)/(1,1)/(2,2)) — clears the pre-registered PASS bar (>=5, >=2 classes) with wide margin.
  The 2 deaths (branches 0, 14) fail only the #324 bend gate (residual stays near machine
  precision — genuine closures) after real, monotonic, smooth (not single-step-jumpy) V_inf
  drift under continuation (+1.19 km/s and +0.27 km/s respectively) that tracks a "longer
  TOF/higher n_rev drifts more" pattern (branch 0 is the longest-TOF branch in the whole
  population). **Floor-hugger control (ids 11/15/17/18/21): 4/5 SURVIVED** (only 18 died) —
  this is ABOVE the pre-registered ">=3 survive => proxy non-discriminating" flag threshold.
  **FLAGGED EXPLICITLY AND PROMINENTLY per spec: the #573 ">=6.0 deg ecc-robust" bend-margin
  proxy is NOT a reliable discriminator of real-eccentricity survival** — this does not overturn
  the PASS verdict (which is independent of the control-check outcome by design) but materially
  discounts confidence that the proxy specifically identified the survivors; read the 15/17
  population result as "most nearby closures survive a real eccentricity perturbation," not as
  the proxy having worked. **Mirror-pair check:** {branch 2, branch 19} both survive but
  **DID NOT MERGE** under the C3-extended dedup (final Omega ~180 deg apart still, V_inf now
  ~6% different); {branch 4, branch 14} **NOT_BOTH_SURVIVING** (4 survives, 14 died) — both
  outcomes match the a priori physical prediction that the exact circular (Omega+180, u+180)
  mirror degeneracy (which holds only because a circular orbit's r(nu) is constant) breaks under
  eccentricity (r(nu) != r(nu+180) in general), so the two branches were free to and did diverge.
  **Framing (mandatory, honored): quasi-cycler-class evidence only, same standing as #312's own
  family (V2 fails on drift by design), NOT a ballistic-cycler finding, NOT a novelty claim.**
  Did NOT build Stage B, fetch SAT441, touch `catalogue.yaml`, or write to
  `data/empty_regions.jsonl` (no KILL, so the conditional-negative writeback clause does not
  apply) — none of this run's mandate. Full per-branch survivor table + methodology writeup:
  `docs/notes/2026-07-12-574-stageA-eccentric-kill-gate.md`. Artifacts:
  `scripts/run_574_titan_iapetus_eccentric_kill_gate.py`,
  `data/probe_574_titan_iapetus_eccentric_kill_gate.jsonl`. `uv run pytest
  tests/data/test_outstanding_structure.py -q` clean; ruff clean; `catalogue.yaml`/
  `empty_regions.jsonl` untouched so the full data/search ratchet suite is not required. **Next
  step: the C6 A->B transition gate** (docs/notes write-up done here; still needs one Fable pass
  specifically on whether Stage A was executed AS PRE-REGISTERED before any Stage B dispatch —
  not performed by this run) — NOT run here, out of this task's scope.
  **C6 TRANSITION GATE (2026-07-12): CONFIRMED — Fable independently re-verified C1
  (parameterization), C2 (positive control), and C3 (counting discipline) against the actual
  script and raw jsonl (not just the write-up); all three hold, PASS verdict correctly derived,
  no material deviations. Clear to proceed to Stage B.**
  **STAGE B KERNEL FETCH (2026-07-12): DONE, verified.** `sat441.bsp` (631 MB) fetched from NAIF's
  public archive to `~/dev/references/kernels/sat441.bsp` with explicit user confirmation
  (filename/source/size stated first, per this project's download-permission discipline).
  Verified as the current, correct kernel — no successor covers Titan/Iapetus (checked the
  newer-numbered sat455/456/457/459 kernels on NAIF's server; they cover only unrelated
  small/irregular moons discovered 2020-2023) and JPL SSD's own ephemerides page names SAT441 as
  current (Jacobson 2022, AJ 164:199-217). Download integrity verified (exact expected byte
  count, 661592064 bytes, correctly typed as a NASA SPICE binary). **Post-fetch coverage
  verification per C4 (not just trusting the filename)**: loaded via `spiceypy` in the project's
  own venv, `spkobj()` + `spkcov()` independently confirm both Titan (NAIF 606) and Iapetus (NAIF
  608) are present with coverage 1749-12-29 to 2250-01-05 — vastly exceeds any plausible
  validity_window epoch range. **Next step: build the real V1-V4-strict-equivalent Saturn
  validation chain + productize the #573/#574 corrector, per Stage B's remaining scope
  (items 1 and 3 below) — not yet started.**
  **STAGE B ITEMS 1+3 RESULT (2026-07-12, Sonnet) — VERDICT: 0/15 PASS, honest structural
  negative.** Built `src/cyclerfinder/genome/titan_iapetus_corrector.py` (productized #572-574
  closure engine: `TitanIapetusClosureParams` 4-DOF input contract, `kepler_state_3d` eccentric
  propagator, `evaluate_closure`/`closure_passes_gate` wrapping the #324 gate verbatim; C2 e=0
  reduction positive control ported as a real pytest regression, not a throwaway smoke test) and
  the Saturn-specific V2→V3→V4→V4-strict chain (`src/cyclerfinder/data/validation/
  {v2_saturn_3d,v3_saturn_3d,v4_saturn,v4_saturn_strict}.py`) — deliberately Saturn-specific
  rather than reusing the generic Uranian `v2_moontour.py`/`v3_3d.py` drivers, because those
  compose on circular-coplanar `_moon_state` and the whole point of #574 is that Titan/Iapetus's
  real eccentricity+inclination is NOT circular-coplanar; feeding these candidates through the
  generic driver would have silently dropped the exact fidelity axis under test. V4-strict was
  written FRESH with BOTH #567 fixes inherited from the start (continuous Lambert branch
  selection by actual propagated terminal offset; `FAILURE_MODE_PLANET_CROSSING` tagged-not-
  excluded) — confirmed firing correctly on a real trigger (branch 2 at epoch 2000-09-15,
  periapsis inside Saturn's R_eq), pinned as a regression test against a genuine case, not a
  synthetic one. `ensure_sat441_kernel()` added to `verify/spice_kernels.py` mirroring
  `ensure_jup365_kernel`/`ensure_pluto_kernel`. 14 new tests across
  `tests/genome/test_titan_iapetus_corrector.py` + `tests/data/test_saturn_v2v3v4_gauntlet.py`,
  all passing; `ruff check`/`ruff format --check`/`mypy` clean on every new module.
  **Ran the full gauntlet on all 15 Stage-A survivors** (branch ids 1,2,3,4,5,6,7,8,9,10,12,13,
  16,19,20 — verified programmatically against the jsonl's `eccentricity_robust AND survives`
  fields inside the script itself, not hand-transcribed) via
  `scripts/run_574_stageB_saturn_gauntlet.py`. **Load-bearing finding: every one of the 15
  candidates FAILS the governing V2 multi-cycle gate** (run at n_cycles ∈ {3,5,10}, mirroring
  #566's grid). 3 candidates (ids 1, 10, 16 — all `n_rev != (0,0)`) fail to even complete 3
  cycles: the same-`n_rev` Lambert transfer physically ceases to exist past cycle 0 (confirmed at
  e=0 too, so not an eccentricity effect — a genuine multi-rev feasibility-window violation as
  the encounter geometry drifts). The other 12 complete all cycles but with per-cycle
  V_inf-continuity residual growing to 0.3–6.4 km/s and drift of ~0.7–2.4 million km — past both
  the strict 50,000 km floor AND (checked across the full {3,5,10} grid) the #566-style 0.5 km/s
  quasi-bounded floor (only branch 6 is quasi-bounded at nc=3 and nc=5, not nc=10 — so it does not
  clear the grid-wide bar either). **Root cause, confirmed by direct instrumentation (not
  assumed): unlike the Uranian #558-#569 family, whose `tof_scale`/`rel_offset_deg` were found by
  #563's DEDICATED symmetric/commensurate-closure enumeration specifically so a fixed-TOF
  multi-cycle repeat reproduces the SAME encounter geometry, the Titan-Iapetus #571-#574 closures
  were found by a free (Omega, tof_scale, rel_offset) search for a SINGLE V_inf-continuity
  closure with NO periodicity/commensurability constraint at all** (branch 1's leg-0 transfer
  angle measured 117.6°/78.0°/38.3° at cycles 0/1/2, checked at e=0). This means most of Stage
  A's 15 "closures" are single V_inf-continuity transfers, not repeating cyclers, under a literal
  multi-cycle test — a different (and more fundamental) finding than "the family dies on real
  eccentricity", orthogonal to the eccentricity question Stage A itself answered. V3/V4/V4-strict
  were still computed for the 12 candidates whose V2 completes 3 cycles (cheap, <1s/stage, and
  informative on their own terms): V3 agrees with V2 to near machine precision for all 12 (as
  expected, same analytic model / different integrator); V4 (J2 + 8-moon third-body scipy)
  passes 9/12; V4-strict at one reference epoch (2000-06-21, reused from the #338/#566/#559
  precedent) passes only 3/12 (branches 2, 8, 19), the rest failing mostly via
  `planet_crossing_infeasible` (real synodic geometry, not a solver artifact) plus one
  `lambert_no_solution` (branch 12) and one converged-but-past-the-agreement-floor case (branch
  9). **None of this changes the governing verdict — since V2 never clears PASS or
  FAIL_QUASI_BOUNDED for any candidate, the full-chain verdict for all 15 is `FAIL_AT_V2_*`: 0/15
  reach PASS or PASS_AS_QUASI_CYCLER.** Full per-candidate table + methodology:
  `docs/notes/2026-07-12-574-stageB-saturn-gauntlet.md`; raw per-cycle data:
  `data/gauntlet_574_saturn_stageB.jsonl`. **Epoch-sensitivity spot check** (6 points across
  2000+2015, NOT a full annual/daily sweep, per the #568 duty-cycle framing that a raw
  single-epoch result isn't the final word): the 3 V4-strict-passing candidates (2, 8, 19) each
  PASS at multiple non-adjacent sample epochs and FAIL at others (genuine synodic-boundary
  structure matching the #567/#568-established pattern, not a knife-edge single point) — reported
  for completeness only; does not change the headline result since all three already fail the
  governing V2 gate regardless. **Framing (mandatory, honored): quasi-cycler-class evidence only,
  same standing as #312's own family — NOT evidence the Stage-A eccentric closures are
  computationally wrong (branch 6's cycle-0 residual is machine precision, independently
  reproduced by this dispatch's own tests) and NOT a novelty claim; the negative is a
  periodicity-formulation gap in how #571-574's discovery search was originally posed.**
  Explicitly NOT done here: Opus+Fable adjudication of this result, `catalogue.yaml`/
  `empty_regions.jsonl` writeback, a full annual/daily V4-strict epoch-robustness sweep, or
  re-scoping the discovery search toward multi-cycle commensurability (a defensible follow-up
  given this root cause, but a new task). `uv run pytest tests/genome/
  test_titan_iapetus_corrector.py tests/data/test_saturn_v2v3v4_gauntlet.py -q` clean (14
  passed); full `uv run pytest tests/data tests/search -q` ratchet run (touches shared
  `verify/spice_kernels.py`) — see the commit for the result.

- **#575** (P1, plan under Fable review before dispatch — do NOT auto-fire) — Titan-Iapetus
  direct symmetric-closure re-search, using the #563 method (NOT #571's looser free-closure
  search). **Motivation:** #574's Stage B gauntlet found 0/15 PASS, root-caused (not a bug) to a
  periodicity-formulation gap in #571: the discovery search that found the 69 Titan-anchored
  candidates never enforced a genuine repeating-cycle condition, it only checked whether A
  Lambert closure exists — most of the 22→15 survivors turned out to be one-off transfer
  geometries, not periodic quasi-cyclers. #563's method (rel_offset ∈ {0°, 180°}, exactly
  commensurate tof = n·T_syn/2, matched leg revolution counts — the classical symmetric
  "perpendicular crossing" condition) is what actually made the Uranian family's 30 closures
  genuine repeating cycles, not #558's original grid search — #571 never applied this method to
  Titan-Iapetus, it reused a different, weaker closure criterion from earlier in the chain.
  **Scope:** directly construct (not grid-search) Titan-Iapetus candidates at commensurate
  tof = n·T_syn_Titan-Iapetus/2 for each n_rev=(n0,n1) combination in a reasonable range, at
  rel_offset ∈ {0°, 180°}, reusing #571/#572's Lambert-closure + #324 gate machinery verbatim
  (do not rederive the gate logic). This is a fresh CIRCULAR-COPLANAR idealized search first
  (matching #571's own starting point) — do NOT skip straight to inclined/eccentric; if this
  produces genuine symmetric closures, THEN run them back through the SAME #572→#573→#574
  pipeline (inclination closure test, population widening, eccentricity kill-gate, Saturn
  V2-V4-strict gauntlet) that was already built and validated this session — all of that
  machinery is reusable as-is, only the discovery-search stage (the #571-equivalent step) needs
  to be redone correctly.
  **FABLE PLAN REVIEW (2026-07-12): CONFIRMED WITH CORRECTIONS — core proposal sound, T_syn
  concept verified to transfer cleanly (generic two-body formula, no Uranus-specific assumption;
  Titan-Iapetus T_syn≈19.96d from sourced periods), #571's omission confirmed as a plain
  oversight not a justified exclusion. Four corrections, two load-bearing, MANDATORY before
  dispatch:**
  **(C1 — replaces the weak positive control) Titan-Rhea is NOT a usable positive control —
  it's structurally unfalsifiable** (its only anchor, #320's internal basin, sits at
  rel_offset=285°, which the symmetric {0°,180°} construction cannot recover by design, and it
  has no symmetric-class golden of its own). **Use a proper method golden instead**: genericize
  `scripts/enumerate_563_symmetric_closures.py` (currently hardcodes `PRIMARIES["Uranus"]`; the
  gate functions it imports already take `primary=` since #571's own genericization) and require
  EXACT reproduction of the known Uranian symmetric-closure table — at minimum #312 itself
  (Umbriel-Oberon n=5, n_rev=(1,1), V∞=0.965, residual ≤1e-13) and its n=7 sibling; the full
  30-closure run is only 8.1s, reproduce all of it. This validates the genericization end-to-end
  against a real machine-precision golden, matching #571's own discipline for `scan_558`.
  **(C2) Add a two-sided multi-cycle repeat-instrumentation control**, using the per-cycle
  transfer-angle check already built in `run_574_stageB_saturn_gauntlet.py`, at e=0 (cheap,
  seconds): (a) POSITIVE — every new symmetric candidate must repeat identical geometry to
  machine precision BY CONSTRUCTION; (b) NEGATIVE — a known-bad #571 candidate (e.g. branch 1,
  documented 117.6°/78.0°/38.3° drift) must NOT repeat. This directly tests the property whose
  absence killed #574, at the discovery stage instead of after the whole pipeline is spent.
  Also cheaply cross-check: scan #571's original 187 candidates for any within tolerance of
  rel_offset∈{0°,180°} AND commensurate tof — any such point must also appear in the new
  enumeration (an empty intersection is itself informative, supporting outcome (b) below).
  **(C3 — LOAD-BEARING) The reused pipeline's refinement stages can silently walk a genuinely
  symmetric candidate BACK OFF the commensurate manifold, reproducing 0/15 after all the compute
  is spent, indistinguishable from a real negative.** #573's Nelder-Mead refines tof_scale in a
  ±0.1 box (≈±3.6 days, vs T_syn/2≈10 days) and #574 Stage A's continuation refines tof_scale in
  its 4D space — BOTH against a single-cycle closure residual, the exact criterion class that
  produced the original #571 artifact. Do NOT hard-fix tof at exactly n·T_syn/2 (the true
  inclined/eccentric periodic solution's tof isn't exactly that value — hard-fixing would
  false-kill genuine candidates); INSTEAD: (i) report a commensurability-drift diagnostic
  (|tof − n·T_syn/2| per branch) after every #572/#573 refinement step and every #574 Stage A
  continuation step; (ii) re-run the C2 multi-cycle check after each refinement stage, not just
  at the end, so periodicity loss is caught immediately rather than only surfacing at Stage B's
  V2 gate.
  **(C4) Framing addendum**: Titan/Iapetus's real eccentricity (≈0.028 each, non-resonant pair,
  period ratio ≈4.975) causes apsidal-phase precession cycle-to-cycle that the T_syn
  commensurability condition alone doesn't lock — 7-25× stronger than the negligible Uranian
  effect. A C2-verified genuine symmetric candidate CAN still fail or quasi-bound Stage B's real
  V2 multi-cycle check under real eccentricity — that outcome means "real perturbations degrade
  a genuine cycle" (legitimate quasi-cycler-class evidence, same standing as #312) and must NOT
  be misread as "the symmetric method also failed."
  **Numbers, pinned (Fable-derived from #571's own verified data, not assumed) — do not leave to
  implementer judgment**: n_rev ∈ {0,1,2,3}² (16 combos, matches #563's own `N_REV_MAX` and
  spans the (0,0)/(1,1)/(2,2) classes that actually appeared in the Titan-Iapetus population);
  n ceiling = 10, derived as `floor(2·tof_scale_max·√(P_Titan·P_Iapetus)/T_syn)` with
  tof_scale_max=3.0 (VERIFIED against `scan_571_*.jsonl`'s own meta record — its actual tested
  tof_scale grid tops out at 3.00, not assumed). Total search space: 2 directions × 10 × 16 × 2
  rel_offset = 640 candidates, seconds of compute.
  **Positive control**: per C1 above (Titan-Rhea demoted to an unlabeled tertiary smoke test
  only, if included at all — it is NOT the primary validation).
  **Honest possible outcomes, pre-registered:** (a) a genuine symmetric Titan-Iapetus family
  exists and survives the full pipeline — worth the investment; (b) no symmetric closures exist
  at all for this pair (the earlier 69/22/15 population was ENTIRELY an artifact of the weaker
  closure criterion, not a hint of a real family) — a clean, valuable negative closing this
  specific pair for good; (c) symmetric closures exist but are few/marginal — re-adjudicate scope
  at that point, same as every prior stage in this chain. **Explicitly out of scope:** any new
  V4-strict/gauntlet code changes (reuse #574's Saturn chain verbatim — `TitanIapetusClosureParams`
  confirmed generic enough to accept symmetric-construction candidates without modification, only
  a schema/loader pointer change since `run_573...` currently hardcodes `SCAN_571_PATH` and
  #571's exact output field names); `catalogue.yaml` edits; re-running Titan-Rhea's OWN
  candidates through the full pipeline (it's V0-known, not novelty-eligible).
  **Recommended models:** Sonnet for the direct-construction search + reuse of existing gate/
  pipeline machinery (mechanical, spec-complete, matches #563's own precedent exactly). Opus +
  Fable for adjudicating any surviving family before it proceeds through the reused pipeline
  stages, and again before any catalogue writeback — same two-model discipline this whole chain
  has used throughout.
  **RESULT (2026-07-12, Sonnet mechanical pass, steps 1-6 partial per dispatch scope):**
  **C1 golden PASSED exactly** — genericized `scripts/enumerate_563_symmetric_closures.py`
  (`--primary`/`--moons`/`--tof-scale-max`/`--out`, construction logic untouched; `residual_at_point`/
  `gate_candidate` already took `primary=` since #571) reproduces the full pre-#575 Uranian table
  **bit-for-bit** (0.0 diff on every one of 60 pass records' residual/V∞ fields) with zero CLI args,
  including #312 itself (Umbriel-Oberon n=5, n_rev=(1,1), rel=180°, residual=4.22e-15, V∞=0.9646
  km/s) and its n=7 sibling (n_rev=(2,2), rel=0°, residual=5.33e-15). Tests:
  `tests/scripts/test_enumerate_563_symmetric_closures.py` (4 tests, incl. a full-run byte-diff
  against the committed golden jsonl).
  **Titan-Iapetus direct construction (pinned numbers verified: n_max=10 from #571's own
  tof_scale_max=3.0, confirmed against `scan_571_saturn_titan_iapetus.jsonl`'s `_meta`)**: 640
  candidates evaluated (2 directions × 10 × 16 × 2, 2.0s) →
  `data/enumerate_575_titan_iapetus_symmetric_closures.jsonl` → **9 unique symmetric closures pass
  ALL gates** (residual 3e-15 to 2e-12 km/s, #324 physical bend gate, DOP853 cross-check) at
  n∈{2,3,4,6,8,9,10}, n_rev∈{(0,0),(1,1),(2,2)} — a genuine, non-empty positive result at the
  coplanar-idealized construction stage (outcome (a)/(c), not the clean-negative (b)).
  **C2 two-sided repeat control PASSED**: reused `run_v2_saturn_3d`
  (`cyclerfinder.data.validation.v2_saturn_3d`, the mechanism `run_574_stageB_saturn_gauntlet.py`
  itself imports and drives — not reimplemented) at e_titan=e_iapetus=inclination=0 via
  `TitanIapetusClosureParams` overrides (pinned exact reduction to the coplanar model, per that
  module's own C2 positive control). All 9/9 new candidates complete 3 cycles with
  closure-residual repeat at machine precision (4e-14 to 3e-12 km/s); the #571 branch-1 negative
  control (loaded directly from `probe_574_titan_iapetus_eccentric_kill_gate.jsonl`'s own e=0
  stage, never hand-transcribed) correctly fails to complete even 1 cycle (Lambert infeasible),
  matching #574 Stage B's own documented 117.6/78.0/38.3° drift finding. Cross-check against
  #571's original 187 all-gates-passed candidates (69 Titan-anchored + 118 Iapetus-anchored, per
  `scan_571_saturn_titan_pairs_index.jsonl`): **0/187 land within tolerance (2°, 0.5 d) of the
  symmetric condition** — a clean empty intersection confirming #571's free-search population
  never touched the symmetric manifold; the 9 #575 closures are a genuinely independent discovery,
  not a re-validation of any #571 candidate. Scripts:
  `scripts/probe_575_titan_iapetus_repeat_check.py`; tests:
  `tests/scripts/test_probe_575_repeat_check.py` (3 tests). **Methodological note**:
  `run_v2_saturn_3d`'s `max_drift_kms` is an INERTIAL-frame absolute-position metric (identical to
  `v2_moontour`'s own convention) — it is large (90k-2.4M km) for every #575 candidate here, but
  this is NOT a red flag: #312's own catalogued SILVER fails the identical strict 50,000 km floor
  under the identical mechanism (`silver_327_moontour_v2_verdicts.jsonl`: max_drift_kms=515,499 km,
  passes_v2=False, FAIL_QUASI_BOUNDED) because T_syn is generically incommensurate with the
  anchor's own individual orbital period — a separate axis from the rel_offset/tof symmetric
  condition. C2's discriminator is therefore Lambert-completion + closure-residual repeat (which
  is exactly what killed #571 branch 1), not this drift metric.
  **Step 6 (pipeline continuation, PARTIAL — stopped after stage 1 per dispatch scope, did NOT run
  #573/#574/Stage B)**: pushed the 9 coplanar seeds through a #572-equivalent inclination closure
  (`scripts/probe_575_stage1_inclination_closure.py`, reusing #572's `evaluate_point`/
  `sweep_node_alignment`/`_smoke_test_reduction` verbatim, Iapetus at the real 15.5° inclination,
  Omega + tof_scale free in #572's own bounded Nelder-Mead). **6/9 seeds find a gate-passing 3D
  closure basin, but C3's drift diagnostic + re-run C2 check show 0/6 still repeat to machine
  precision post-refinement**: tof drifted 0.23%-11.6% of T_syn/2 off the commensurate seed value,
  and closure residual over 3 cycles grew to 0.11-1.65 km/s (not machine precision); the two
  most-drifted (10.5%, 11.6%, both n=10) lost Lambert feasibility entirely past cycle 0 — the exact
  #571-branch-1 failure signature, now reproduced by refinement-induced drift rather than a
  free-search origin. A complementary Omega-ONLY sweep at the EXACT fixed commensurate tof (zero
  C3 drift by construction, `scripts/probe_575_stage1b_omega_only_closure.py`) still finds 5/9
  gate-passing single-cycle closures, and **0/5 of those repeat either** (residual ~1e-3 km/s at
  cycle 0, not machine precision) — proving periodicity loss here is NOT solely the C3 tof-drift
  mechanism; introducing Iapetus's real inclination breaks the multi-cycle repeat property in its
  own right for this family, at least at the single best-residual basin checked per seed (a
  basin-by-basin exhaustive check — #572 found 4-8 basins per candidate — was NOT run, a cheap
  natural follow-up). Data: `data/probe_575_stage1_inclination_closure.jsonl`,
  `data/probe_575_stage1b_omega_only_closure.jsonl`.
  **Framing (C4-consistent)**: NOT a novelty claim, NOT a catalogue writeback. This is
  quasi-cycler-CLASS evidence about the idealized search space (same standing as #312), and the
  informative finding is that the genuine coplanar symmetric family found here does not
  straightforwardly survive the SAME inclination-extension method (#572) that #574 Stage B later
  found also broke down under real eccentricity for the #571 population — a materially different,
  earlier-stage breakdown (inclination alone, before eccentricity is even introduced), diagnosed
  cleanly via C2/C3 instrumentation exactly as this task's Fable corrections required.
  `uv run ruff check`/`ruff format --check` clean on all new/modified files; targeted
  `uv run pytest tests/scripts/test_enumerate_563_symmetric_closures.py
  tests/scripts/test_probe_575_repeat_check.py tests/genome/test_titan_iapetus_corrector.py
  tests/data/test_saturn_v2v3v4_gauntlet.py -q` — 21 passed.
  **FINAL ADJUDICATION (2026-07-12, Opus — culminating trust-bearing call, advisor tool was
  unavailable so decided on the full read of #571-#575 + raw probe data): STOP. Clean,
  method-conditional negative; Titan-Iapetus symmetric-closure thread CLOSED.** Wrote
  `docs/notes/2026-07-12-575-titan-iapetus-final-adjudication.md` and stamped
  `data/empty_regions.jsonl` region `saturn-titan-iapetus-symmetric-closure-inclination-empty-575`.
  Core reasoning: (1) The near-misses are robustly infeasible, not tantalizing — read directly
  from the probe jsonl: Stage-1 (Omega+tof-free continuation) closes a single cycle to machine
  precision but then DIVERGES by up to 2.44M km (~2× Titan's orbital radius) over 3 cycles, with
  tof drifting 0.23%–11.6% off commensurate; the two n=10 seeds lose Lambert feasibility past
  cycle 0. Stage-1b (zero-tof-drift, Omega-only, Omega landing on the 0/180° symmetric axis)
  cannot null even a ~6–21 m/s per-cycle discontinuity. That is a ~12-order-of-magnitude gap from
  the Uranian coplanar 4e-15 km/s — a different object, not "harder to close." (2) The negative is
  TRUSTWORTHY not a broken filter: the same pipeline reproduces the Uranian golden bit-for-bit AND
  the 9 coplanar Titan-Iapetus closures repeat to machine precision — it finds genuine repeating
  cyclers when they exist. (3) Honest Uranian comparison = a genuinely DIFFERENT regime, not the
  same correction harder: Umbriel-Oberon mutual inclination ~0.1–0.2° (coplanar model near-exact,
  inclination a negligible correction) vs Titan-Iapetus mutual inclination ~15° (coplanar model
  never a good starting point, inclination a leading-order effect). (4) The one method not tried —
  a from-scratch node-locked inclined symmetric construction — is a structurally DISCONNECTED
  family unreachable by continuation from coplanar seeds; it is essentially the #552 general
  3D/inclined build KILLED 2026-07-10, has low expected yield (a marginal inclined family must
  still survive the #574 eccentricity kill-gate + Saturn gauntlet that already killed 0/15 of the
  #571 population, with 7-25× stronger eccentricity than Uranus), and is not novelty-eligible. The
  entire #572→#575 probe sequence WAS the cheap gate deciding whether the #552 build is worth
  scoping; the gate has now returned a robust negative twice (eccentricity via #574, inclination
  alone via #575), so building the very thing the gate was designed to gate would invert the
  thread's discipline. Recorded as the named, method-versioned re-open CONDITION in the registry
  (per the "empty is method-conditional" rule), NOT dispatched — #576 stays next-unused. NOT a
  novelty claim, NO `catalogue.yaml` edit; quasi_cycler-class evidence, same standing as #312.

- **#576** (P1, plan under Fable review before dispatch — do NOT auto-fire) — Galilean-moon
  direct symmetric-closure search, using the now-twice-proven #563/#575 method
  (`scripts/enumerate_563_symmetric_closures.py`, genericized, golden-validated bit-for-bit
  against the Uranian family). **Motivation:** this specific method — exact commensurate
  construction at rel_offset ∈ {0°,180°}, tof = n·T_syn/2, gated + repeat-instrumentation
  verified — has NEVER been applied to Jupiter's Galilean moons. Prior Jupiter work used
  different methods: the now-retired torus/heteroclinic linking-number lane (#536/#545, dead per
  #555's final verdict); a real-ephemeris SOBOL-SAMPLED broad joint search (#501: 6 sequences ×
  512 samples, 213 feasible, 26 shot, 0 closed, positive-controlled against Liang Member D,
  "clean empty-region map" registered); and a DSM/leveraging multi-arc campaign (#465: 8 in-band
  closures found, all confirmed V0-known reproductions of published Campagnola-Russell/GTOC6
  tours, no novelty). None of these is #563's exhaustive direct construction.
  **Load-bearing methodological question for the Fable review to resolve, not assume:** does
  #501's "0/3072 real-ephemeris closed, positive-controlled" result already preemptively kill
  any idealized-model symmetric family here, or does #563's exhaustive commensurate construction
  explore a genuinely different part of parameter space that a SAMPLED search (even a broad,
  positive-controlled one) can structurally miss — exactly the same relationship #558's original
  grid search had to #563's own direct construction for Uranus (#563's own entry: "a basin
  narrower than the grid resolution... can fall entirely between sample points and never trigger
  the initial gate at ANY tested point"). Resolve this BEFORE dispatch, not after — if #501's
  breadth search already structurally covers what #563's method would find, this task adds
  nothing and should not run.
  **Mass/gate check (favorable, but verify don't assume):** all 4 Galilean moons are massive
  (`core/satellites.py`: Io GM=5959.9, Europa GM=3202.7, Ganymede GM=9887.8, Callisto GM=7179.3
  km³/s² — every adjacent pair comparable to or exceeding Titan's 8978.14, the ONE Saturnian
  moon that cleared the two-sided bend gate) — unlike Saturn, bend-gate infeasibility is not
  expected to be a limiting factor for any adjacent Galilean pair. **Inclination check
  (UNVERIFIED, mandatory before dispatch — this is exactly the assumption that doomed the
  original #571 Titan-Iapetus search):** the Galilean moons are widely believed near-coplanar
  with Jupiter's equator and with each other (much smaller than Titan-Iapetus's ~15.5°), but
  this must be independently VERIFIED with sourced numbers (JPL SSD or equivalent) before
  building the search plan around a coplanar-favorable assumption, not asserted from general
  knowledge — the #571 negative was root-caused exactly to an unverified/wrong assumption about
  what closure criterion mattered; do not repeat that pattern on the inclination question here.
  **Scope:** (1) verify real mutual inclinations for adjacent Galilean pairs (Io-Europa,
  Europa-Ganymede, Ganymede-Callisto — and note Io:Europa:Ganymede's exact 1:2:4 Laplace mean-
  motion resonance, which may make T_syn commensurability behave qualitatively differently than
  for a non-resonant pair like Titan-Iapetus or the Uranian moons — flag this explicitly rather
  than assuming the Uranian/Saturnian playbook transfers unchanged); (2) direct symmetric
  construction per #563's method for each adjacent pair, reusing the genericized script verbatim
  (primary=Jupiter); (3) the SAME repeat-instrumentation verification #575 used (positive: new
  candidates must repeat to machine precision; negative control: some known-non-periodic prior
  Galilean candidate, if one exists in `data/scan_433_jupiter_galilean.jsonl` or #501's own
  output, to prove the instrumentation discriminates); (4) if genuine closures survive to
  real-ephemeris gauntlet territory, note Jupiter's own real-ephemeris SPICE kernel (JUP365,
  already fetched at `~/dev/references/kernels/jup365.bsp` per earlier session work) is already
  available — no new kernel fetch needed. **Explicitly out of scope:** any new torus/linking-
  number work (dead method, do not revive); re-running #501's or #465's own already-closed
  searches; `catalogue.yaml` edits; claiming novelty before a literature check clears any
  survivor (Galilean tours are a much more heavily-published territory than Uranus/Saturn —
  Campagnola/Strange/Russell/GTOC6 — expect a HIGHER prior on "known" than either prior system).
  **Recommended models:** Sonnet for the mechanical direct-construction search reusing proven
  machinery (spec-complete, matches #575's own precedent exactly). Opus + Fable for adjudicating
  any survivors, and mandatorily for the literature-risk assessment given how heavily-mined this
  specific system is in the published literature.
  **FABLE PLAN REVIEW (2026-07-12): CONFIRMED WITH CORRECTIONS — the load-bearing #501 question
  resolves IN FAVOR of running (confirmed by reading #501's actual verdict note, not its summary
  line). Five corrections, two load-bearing:**
  **(Q1 resolved) #501 does NOT preemptively cover this territory** — structurally different on
  three independent axes: (1) different MODEL (#501 is real-ephemeris n-body-shot, #563/#576 is
  idealized circular-coplanar Kepler-Lambert — #312's whole 30-member family exists in exactly
  this idealized-but-real-eph-drift-failing gap, so a #501-style lane could never have surfaced
  it); (2) different DEGREES OF FREEDOM (#501 sampled EPOCH along the real ephemeris flow, with
  tof emerging from the Lambert solve — rel_offset was never an independent coordinate and tof
  was never constructed commensurate; #563/#576 directly enumerates the exact discrete
  rel_offset∈{0°,180°}/tof=n·T_syn/2 set — measure-zero points a 512-sample epoch scan hits with
  probability ~0); (3) different FUNNEL (#501 shot only 26 top-5-prefilter survivors with an
  n-body shooter from unconstructed seeds — 0/26 non-convergence from BAD SEEDS says nothing
  about the symmetric class's existence). Per the negative-results registry's own
  method-conditional rule, a structurally different method is licensed to re-open #501's stamped
  region. **Honest-framing corollary, mandatory**: conversely, any #576 closures do NOT
  contradict #501's real-eph empty stamp — ANNOTATE, don't unstamp; deliverable is idealized
  quasi-cycler-class evidence only (same standing as #312/#575), and #501's own real-eph negative
  keeps the prior on eventual real-eph survival LOW, not high.
  **(C1-scope, load-bearing) Widen to ALL 6 pairs** (C(4,2), not just the 3 adjacent) — Fable
  independently verified via the same two-sided-gate method `verify_571_gate_analytics.py`
  embodies that ALL 6 pairs clear the #324 gate with 2-5× margin (even the widest, Io-Callisto:
  min-achievable V∞ 4.82/3.24 km/s vs ceilings 8.25/7.77 km/s) — the exact opposite of Saturn's
  small-moon exclusions. Marginal compute for the extra 3 pairs is seconds; independently
  recompute these numbers via the dispatch's own run of `verify_571_gate_analytics.py`-style
  logic before trusting them (Fable's numbers are Hohmann-tangent hand-computation, matching that
  script's method but not run through it directly).
  **(Q3) `core/satellites.py` has NO inclination field at all** (confirmed — the dataclass carries
  mu/radius/a/safe_alt only, same gap #554 found for Triton's retrograde flag). Fable's own
  knowledge (high confidence, standard JPL SSD values, but explicitly NOT repo-sourced) suggests
  a favorable picture — mutual inclinations ≤~0.5-0.65° depending on node alignment (Uranian-class,
  ~25-50× smaller than Titan-Iapetus's 15.5°) — but per the #571 lesson this MUST be independently
  sourced (JPL SSD, cited, computed from BOTH i and Ω per pair, not just |i1-i2|) as the dispatch's
  literal first step before any construction, not asserted from general knowledge. Fable's numbers
  are a sanity envelope the sourced answer should fall in, not a substitute for sourcing it.
  **(Q4) Laplace 1:2:4 resonance is NOT a blocker** — T_syn is finite and well-defined everywhere
  (pairwise ratios near-but-not-exactly 2:1: Io/Europa≈2.007, Europa/Ganymede≈2.01; T_syn Io-Europa
  ≈3.53d, Europa-Ganymede≈7.05d, Ganymede-Callisto≈12.5d; n_max stays small, ~4-5/pair, a few
  hundred total candidates). Pairwise rel_offset reachability is NOT resonance-constrained (only a
  future 3-moon I-E-G extension would be phase-locked — note this explicitly so nobody naively
  extends the pair method to triples later without re-deriving). **But the resonance flags a real
  physics point the original plan MISSED**: it FORCES permanent, non-negligible eccentricities
  (Io e≈0.0041, Europa e≈0.0094) — Europa's velocity modulation e·v_orb≈0.13 km/s is ~2.6× the
  0.05 km/s gate floor, same order as the Titan eccentricity effect #571-C4 flagged. Pre-register:
  a genuine coplanar-circular closure degrading under this forced eccentricity is "real
  perturbations degrade a genuine cycle" (legitimate quasi-cycler evidence), NOT "the method
  failed." Minor watch-item: near-2:1 pairs put some Lambert legs near 2π-multiple transfer
  angles — watch for solver-branch/near-degenerate artifacts per
  [[feedback_isolated_sweep_flips_suspect_artifact]], don't hand-tune around them silently.
  **(C2-lit, load-bearing) Russell-Strange 2009's Table 3 is a SOURCED LITERATURE GOLDEN covering
  exactly this territory** — Jovian pairs Ganymede→Io, Ganymede→Europa, Ganymede→Callisto,
  Europa→Ganymede, with specific published members ("Callisto #1" 0.41-yr period, "Ganymede-Europa
  #316", "Europa-Ganymede #131" — digest at `docs/notes/2026-06-30-digest-russell-strange-2009-
  planetary-moon-cyclers.md`, CORPUS_INDEX line 197). This is a stronger and more concrete
  literature obligation than the plan's generic "expect a higher known-prior" framing: (i) any
  survivor MUST be cross-checked against R-S Table 3 FIRST; expect near-certain overlap for
  Europa/Ganymede/Callisto pairs (V0-known, not novel) — the likely genuinely-unpublished slice is
  Io-anchored pairs where Io is the FLYBY (not just the Ganymede→Io direction R-S already lists);
  (ii) REPRODUCE R-S Table 3's Jovian members as an ADDITIONAL positive control, strictly stronger
  than the internal Uranian bit-for-bit golden alone (a real published-record reproduction, not
  just an internal-consistency check) — per [[feedback_golden_tests_sourced_only]]; (iii)
  **critical ordering point**: run the R-S comparison against the UNGATED enumeration output, not
  the #324-gated survivors — R-S's free-return architecture can use one moon as a passive
  (no-bend) target, which the two-sided #324 gate structurally excludes (exactly as the #571
  empty-region stamps document for Titan-Enceladus) — comparing against gated-only output would
  misread a real R-S member as a coverage failure.
  **Minor fixes**: specify the negative control precisely rather than leaving it to implementer
  improvisation — reuse #575's own C2 mechanism, take the negative case from the new Jupiter
  enumeration's own non-symmetric rejects (or a #571-branch-1-style known-bad reference) rather
  than #501's real-eph output (wrong model for an idealized repeat check). Build cost note: the
  script is ALREADY genericized with the Uranian golden committed as a test (#575's own C1) — this
  task's step 2 is a CLI invocation, not a genericization task, cheaper than the entry implies.
  **RESULT (2026-07-12, Sonnet mechanical pass, all 5 steps run in order):**
  **Step 1 (inclination sourcing, PASSED with one flag):** JPL SSD "Planetary Satellite Mean
  Orbital Parameters" table (Laplace-plane mean elements, JUP365, epoch 2000-01-01.5 TDB, accessed
  2026-07-12) gives Io i=0.0/Ω=0.0°, Europa i=0.5/Ω=184.0°, Ganymede i=0.2/Ω=58.5°, Callisto
  i=0.3/Ω=309.1° (the table's own 1-decimal precision, not a truncation introduced here). Mutual
  inclination computed via the full two-plane spherical relation (`cos(i_mut) = cos(i1)cos(i2) +
  sin(i1)sin(i2)cos(Ω1-Ω2)`, mirroring #571/#572's node-aware `iapetus_state_3d` treatment, NOT
  naive `|i1-i2|`) for all 6 pairs: Io-Europa 0.500°, Io-Ganymede 0.200°, Io-Callisto 0.300°,
  Europa-Ganymede 0.637°, Europa-Callisto 0.716°, Ganymede-Callisto 0.412°.
  **FLAG: Europa-Callisto's 0.716° mutual inclination EXCEEDS Fable's ~0.65° sanity envelope** —
  the node-separation term (naive `|i1-i2|`=0.2° vs the real node-aware 0.716°, a 3.6x difference)
  is exactly the effect the #571 lesson warned about. In absolute terms this is still tiny (>20x
  smaller than the ~15.5° Titan-Iapetus mutual inclination that killed #571-#575), so it does NOT
  change the qualitative risk picture, but it is reported prominently as instructed rather than
  silently absorbed. Script: `scripts/probe_576_galilean_inclination_check.py`; data:
  `data/probe_576_galilean_inclination_check.jsonl`.
  **Step 2 (two-sided #324 gate-feasibility, PASSED 6/6):** Independent recomputation via
  `core/flyby.py::max_bend` bisection + Hohmann-tangent minimum-V∞ (same method
  `verify_571_gate_analytics.py` embodies), generalized to a genuine pairwise two-sided check.
  All 6 pairs clear on BOTH sides with 1.7-6.7x margin (worst: Io-Callisto, 1.71x/2.40x — exactly
  reproducing the Fable review's own hand-computed 4.82/3.24 km/s vs 8.25/7.77 km/s ceiling
  numbers). Script: `scripts/verify_576_galilean_gate_analytics.py`; data:
  `data/verify_576_galilean_gate_analytics.jsonl`.
  **Step 3 (direct symmetric-closure construction, all 6 pairs):**
  `scripts/enumerate_563_symmetric_closures.py --primary=Jupiter
  --moons=Io,Europa,Ganymede,Callisto --tof-scale-max=2.0` (tof_scale_max=2.0 verified against
  `data/scan_433_jupiter_galilean.jsonl`'s own `_meta` record — the actual prior idealized
  Lambert-construction sweep range for this system, mirroring #575's own derivation method; NOT
  the Uranian/Saturnian 3.0 bound, which was never tested for Jupiter). n_rev∈{0,1,2,3}²,
  rel_offset∈{0°,180°}. 1856 candidates evaluated across 12 directions (n_max ranges 2-10 per
  pair, driven by each pair's own T_syn) in 4.0s → **36 candidates pass ALL gates** (residual,
  #324 physical bend, DOP853 cross-check) across all 6 pairs (both anchor directions each):
  Io-Europa 1, Io-Ganymede 3, Io-Callisto 3, Europa-Ganymede 1, Europa-Callisto 3,
  Ganymede-Callisto 7 (×2 directions each = 36 total). Data:
  `data/enumerate_576_jupiter_galilean_symmetric_closures.jsonl`.
  **Step 4 (repeat-instrumentation, C2-STYLE DISCIPLINE PASSED):** Reused the generic
  `cyclerfinder.data.validation.v2_moontour.run_v2_moontour` driver (the SAME already-generic
  tool #574/#575's own Saturn-specific 3D wrapper is built ON TOP OF, and the tool the Uranian
  #330/#558 gauntlet already used for this exact circular-coplanar idealized model — no
  genericization needed, `system` auto-resolves Jupiter from the moon names). **POSITIVE: 36/36
  survivors repeat to machine precision** (closure residual 1.6e-14 to 5.1e-12 km/s over 3
  cycles). **NEGATIVE: a non-symmetric reject constructed from this dispatch's own Jupiter
  parameter sweep** (same Ganymede-Callisto pair/n_rev/commensurate-tof as a genuine survivor,
  but rel_offset=90° — outside the {0°,180°} symmetric set the construction ever visits) —
  cycle-0 residual 1.85 km/s (fails the gate immediately), correctly does NOT repeat. Scripts:
  `scripts/probe_576_galilean_repeat_check.py`; tests:
  `tests/scripts/test_probe_576_galilean_repeat_check.py` (2 tests, pin the 36-survivor count and
  both control sides). Data: `data/probe_576_galilean_repeat_check.jsonl`.
  **Step 5 (Russell-Strange 2009 Table 3 literature golden, run against the UNGATED enumeration
  per the mandatory C2-lit ordering correction):** Table 3 transcribed directly from the paper's
  text layer (`pdftotext -raw` on the acquired PDF, row order cross-validated against Table 2's
  sourced synodic periods and the Fig. 6 caption's independently-legible identifiers — NOT
  hand-copied from the digest, which only paraphrased 3/10 rows). Re-derived the full UNGATED
  (residual-gate-only, #324 bend gate NOT applied) candidate set directly from step 3's own
  `residual_at_point` machinery (144 total across all 6 pairs) since R-S's free-return
  architecture treats the target body as massless/no-bend, which this project's two-sided gate
  would wrongly exclude if only the gated 36 were checked. Of R-S's 10 Jovian rows, only the 2
  with `legs=1` (Ganymede→Europa #5, #43) are architecturally comparable to this project's 2-leg
  Anchor-Flyby-Anchor construction (the other 8 are multi-loop resonant tours with a single
  passive target encounter — a different topology, reported but not treated as a reproduction
  target). **Result: 1/2 comparable members has a period match (Ganymede-Europa #43, R-S
  period=14.1d vs our n=2 period=14.108d) but the matching candidates' V∞ (6.7-19.1 km/s) does
  NOT match R-S's stated V∞ (1.87/3.89 km/s) and none clear the physical bend gate — a
  commensurability coincidence at the shared T_syn grid, not a genuine geometric reproduction; the
  other comparable member (#5, period=35.3d ≈ n=5) falls outside this run's n_max=2 bound for
  that pair (tof_scale_max=2.0) and was not reached.** The 8 non-comparable (legs>1) rows show
  10 period matches each for both Ganymede-Callisto rows (period=37.6d ≈ our n=3, several with
  physically-gated, reasonable V∞ 1.5-7.7 km/s) but these are NOT claimed as reproductions given
  the topology mismatch — reported for completeness only. Script:
  `scripts/compare_576_russell_strange_galilean.py`; data:
  `data/compare_576_russell_strange_galilean.jsonl`.
  **Honest framing (mandatory per dispatch scope):** none of the 36 gate-passing symmetric
  closures found here is claimed novel — R-S 2009 already establishes this is heavily-published
  territory and the one architecturally-comparable period coincidence found does NOT reproduce
  R-S's V∞ signature, so it is neither a confirmed reproduction nor a novelty candidate as-is. Any
  further characterization requires a full `literature_check.py` pass + Opus/Fable adjudication as
  a separate follow-up (NOT concluded here, per explicit scope limit). These findings are
  idealized quasi-cycler-class evidence only (same standing as #312/#575) and do NOT contradict
  #501's own real-ephemeris "0/3072 closed, positive-controlled, clean empty-region map" stamp —
  #501 and this task explore structurally different, non-overlapping parts of parameter space (see
  the Fable plan review's Q1 resolution above); #501's stamp is unchanged, not unstamped.

- **#577** (P1, judgment-only — the data all exists, this is the literature-clearance +
  pipeline-readiness gate before any further Galilean investment) — Opus + Fable adjudication of
  the full #576 36-closure Galilean symmetric-closure survivor list. **Motivation:** #576 found
  36 genuine idealized-model gate-passing closures across all 6 Galilean pairs, honestly compared
  against Russell-Strange 2009's Table 3 (the sourced literature golden for this exact territory)
  and found NO confirmed reproduction — but also no clearance, since R-S's other 8 members use a
  structurally different multi-loop architecture this run's `legs=1` construction can't compare
  against. This is genuinely open, not resolved either way, and #576 explicitly deferred any
  further characterization to this task.
  **Scope:** (1) **Broaden the literature check** beyond Russell-Strange alone — run
  `search/literature_check.py` (the mandatory novelty-check-baseline gate, per
  [[feedback_literature_novelty_check_baseline]]) against the full `KNOWN_CORPUS` for all 36
  closures' signatures (V∞ multiset, sequence, period), not just the R-S Table 3 spot-comparison
  #576 did. Given how heavily-published Galilean tours are (Campagnola/Strange/Russell/GTOC6),
  expect a HIGH prior on known-adjacent or known-class-member results — treat anything that
  clears as a genuine surprise requiring extra scrutiny, not a default expectation.
  (2) **Triage the 36 for pipeline-readiness** — the full inclination-closure/eccentricity-
  kill-gate/real-ephemeris-gauntlet pipeline built for Saturn (#572-#574) is directly reusable
  for Jupiter (same JUP365 kernel already available, no new fetch), but running it on all 36
  would be wasteful; rank by robustness (bend margin, V∞, tof) the same way the Uranian/Saturn
  work did, and recommend a small representative subset (if any survive literature clearance) —
  or a full stop if literature clearance kills everything, which given R-S's coverage is a
  live, real possibility that must not be avoided just because it's a less exciting outcome.
  (3) **Address the 8 architecturally-mismatched R-S members explicitly** — #576 found period
  coincidences (10 matches for both Ganymede-Callisto rows, period≈37.6d) that were NOT claimed
  as reproductions due to topology mismatch (multi-loop vs this project's 2-leg construction).
  Assess whether this mismatch is genuine (a real different family) or whether the 2-leg
  construction could be extended/reinterpreted to actually test these specific R-S members
  properly before concluding they're out of scope — don't let "different topology" become an
  unexamined excuse to skip a real comparison that could be done cheaply.
  **Explicitly out of scope:** running the pipeline stages yourself (recommend which candidates,
  don't execute); `catalogue.yaml` edits; claiming novelty (this stays internally-enumerated,
  quasi-cycler-class framing throughout, same as every prior verdict in this chain, UNLESS
  literature_check.py genuinely clears something AND it survives the full pipeline — a much
  higher bar than anything found so far in this specific task).
  **Recommended models:** Opus for the triage/ranking + literature-adjacency judgment calls;
  Fable second-opinion mandatory before any recommendation to proceed to the pipeline stages —
  matches this chain's unbroken discipline (9-for-9 real catches across #558-#576).
  **RESULT (2026-07-12, Opus adjudication; full writeup
  `docs/notes/2026-07-12-577-galilean-adjudication.md`):** **0 of 36 closures clear the broadened
  literature check as novelty; FULL STOP — push none to the pipeline.** (1) Ran the deterministic
  `search/literature_check.py` matcher (`_candidate_anchors`) for all six pairs. Body-set-only
  matching: all 36 collide with 4-7 published Galilean anchors. Under the physically-correct
  `topology_label={"repeated-moon"}` (a 2-body free-return A-B-A shuttle IS the Aldrin/(k1,k2)
  paradigm): 30/36 collide with Liang CGE 2024 + Hernandez IEG 2017 triple-cycler anchors; only
  Io-Callisto (6 closures) structurally clears. **That clear is a CORPUS GAP, not novelty:** the
  single most direct prior — Russell & Strange 2009 "Cycler Trajectories in Planetary Moon
  Systems," JGCD 32(1) DOI 10.2514/1.36610, the exact circular-coplanar two-body free-return
  double-cycler enumeration #563 reimplements ("hundreds of ideal model ballistic cycler
  geometries") — is ABSENT from `KNOWN_CORPUS`, and the present Galilean pump-tour anchors carry
  `{pump-tour,mga-tour}` labels the repeated-moon filter correctly excludes, leaving nothing to
  flag Io-Callisto. Live WebSearch (2026-07-12) corroborates R-S 2009 as the canonical double-
  cycler prior + a 2024 review (DOI 10.34133/space.0036) surveying the whole field + Io-Callisto
  discussed as a dynamically awkward pair. (2) **The "8 architecturally-mismatched R-S members"
  framing is INCORRECT:** the R-S 2009 digest + live search confirm R-S's architecture is the SAME
  object (circular-coplanar ideal model, free-return two-body repeated shuttle, same enumerative
  method); its "legs" column counts free-return arcs/resonances WITHIN the family, not a different
  topology. Three of our six pairs (Io-Ganymede, Europa-Ganymede, Ganymede-Callisto) are in R-S
  Table 1's explicitly-enumerated set; the 37.6-d Ganymede-Callisto period matches are
  known-class-members of R-S's enumerated Ganymede-Callisto family (V∞ differs from the few
  tabulated representatives because Table 3 samples only a handful of hundreds of geometries — not
  a rescue for novelty). No new 2-leg-vs-multi-loop comparison was needed; the shared model+method
  settle it. (3) **Pipeline triage = full stop:** robustness ranking is moot (literature clearance
  kills the population first); the sole structural-clear Io-Callisto also fails the higher
  "survives the pipeline" bar independently (V∞ 5.4-7.3 km/s, Io in the radiation belt, 4.46× a
  ratio, R-S-method-generated, review-discussed — the least attractive pair). Standing: idealized
  quasi-cycler-class known-adjacent (same as #312/#575), NOT novel. Does not unstamp #501. The one
  actionable follow-up is a corpus-accuracy fix + registry stamp, NOT a gauntlet run: scoped as
  **#578** below (add the R-S 2009 repeated-moon anchor to close the false clear; stamp the
  territory).

- **#578** (P2, corpus-accuracy + registry stamp — NOT a pipeline/gauntlet task; Fable-reviewed and
  corrected, ready to dispatch to Sonnet) — close the #577-diagnosed literature-corpus gap and
  stamp the Galilean symmetric-closure territory. **Motivation:** #577 found that
  `search/literature_check.py`'s `KNOWN_CORPUS` is MISSING its single most direct Jovian prior —
  Russell & Strange 2009, "Cycler Trajectories in Planetary Moon Systems," JGCD 32(1) 2009, DOI
  10.2514/1.36610 — the exact circular-coplanar two-body free-return double-cycler enumeration that
  #563/#576 reimplement. Its absence (compounded by the present Galilean pump-tour anchors carrying
  `{pump-tour,mga-tour}` labels the repeated-moon filter excludes) produced a FALSE structural
  clear for the Io-Callisto pair in #577; it is a live false-negative generator for any future
  Jovian double-cycler screen (per [[project_negative_results_registry]] /
  [[feedback_bugfix_invalidates_past_searches]]).
  **CORRECTION (2026-07-12, user-prompted re-check — "we've done this a few times... ensure the
  reference IDs are correct"): the "missing from corpus" framing UNDERSTATED how well-known this
  material already is, and the fix must be stronger than #577's original scope.** Independently
  re-verified: DOI 10.2514/1.36610 is confirmed genuinely absent from `literature_check.py`'s
  `KNOWN_CORPUS` (grepped the source directly, zero hits) — #577's core diagnosis stands. But
  `data/catalogue.yaml` ALREADY carries **30 individual Russell-Strange 2009 member rows plus 2
  older family-level summary rows** (`russell-strange-2009-jovian-multimoon-family`,
  `-saturnian-multimoon-family`; 32 `russell-strange-2009-*` ids total — task #491 ingestion; ids
  `russell-strange-2009-eurgan-131/159`, `-gancal-1/5`, `-ganeur-5/43/316`,
  `-ganio-53/185/403` for the 10 Jovian members, plus 20 Saturnian Titan-Enceladus members), each
  with real sourced V∞/sequence/period from R-S's own Tables 3/5, `validation_level: V0`,
  `orbit_class: cycler`. This is a MORE authoritative "already known" signal than a missing
  corpus anchor — the project has not merely read this paper, it has already catalogued dozens of
  its specific members. Cross-checked #576's actual 36 closures' V∞ against these 10 Jovian rows
  directly: no exact V∞ match (different specific family members at the same pair-directions, not
  literal duplicate rows) — consistent with R-S's own "hundreds of geometries" claim, and with
  #577's topology-class-match conclusion (known-CLASS-member, not exact duplicate). Confirmed
  Io-Callisto genuinely has NO catalogue row among the 10 Jovian R-S entries — supports #577's
  characterization that R-S itself never enumerated this pair (a real gap in the published Table
  1, not just a `KNOWN_CORPUS` omission) — do not read the eventual anchor-fix as implying
  Io-Callisto should suddenly collide falsely; it's a genuinely different, smaller-literature-risk
  case than the other 5 pairs and should be labeled as such in the registry stamp, not folded in
  as identically-covered.
  **CORRECTION 2 (2026-07-12, Fable review of the above — CONFIRMED WITH CORRECTIONS, 2 load-
  bearing fixes, 10th consecutive real Fable catch in this chain):** `_candidate_anchors`
  (`literature_check.py` ~line 1856) requires `seq_set <= anchor.body_set` plus topology-label
  intersection — a SINGLE anchor with `body_set` = the union of all 4 Galilean moons would make
  **Io-Callisto candidates ALSO collide** ({Io,Callisto} ⊆ the union), directly contradicting this
  task's own stated goal of keeping Io-Callisto uncollided. The item-4 registry-stamp claim "5 of
  the 6 pairs ... all directly R-S-covered" is also factually wrong: only 3 pairs are directly
  R-S-Table-1-enumerated; the other 2 covered pairs are covered by DIFFERENT papers, as body-
  subsets of existing triple-cycler anchors, not by R-S itself.
  **Scope (final, post-Fable):** (1) add **THREE separate per-pair `CorpusAnchor`s** for R-S 2009
  (10.2514/1.36610), each `primary="Jupiter"`, `topology_label=frozenset({"repeated-moon"})`,
  `provenance="verified-against-source"`: `russell-strange-2009-ganio` (`body_set={Ganymede,Io}`),
  `russell-strange-2009-ganeur` (`body_set={Ganymede,Europa}`), `russell-strange-2009-gancal`
  (`body_set={Ganymede,Callisto}`) — grounded from BOTH the on-disk digest
  `docs/notes/2026-06-30-digest-russell-strange-2009-planetary-moon-cyclers.md` AND the 30
  already-catalogued `russell-strange-2009-*` member rows (cite the specific catalogue ids as
  additional provenance, not just the digest); confirm Io-Ganymede, Europa-Ganymede, and
  Ganymede-Callisto then collide correctly while Io-Callisto and Io-Europa and Europa-Callisto do
  NOT collide against these 3 new anchors (they're covered, if at all, by other existing anchors —
  see item 4); ALSO add a **4th anchor**, `russell-strange-2009-titenc` (`primary="Saturn"`,
  `body_set={Titan,Enceladus}`, same DOI/topology/provenance), since the identical false-clear
  risk is live today for Saturn too (20 R-S Saturnian catalogue rows, zero `KNOWN_CORPUS` anchor
  covering that body-set at `repeated-moon` topology) — folding this in now is strongly
  recommended by Fable rather than leaving a known-identical gap for a future Saturn double-cycler
  screen to rediscover the hard way; (2) extend the existing `tests/search/test_literature_check.py`
  self-validation with R-S-2009 Galilean AND Saturnian double-cycler signatures asserting
  `status="published"`; (3) **structural fix, concretized per Fable's recommendation**: build (not
  just sketch) a static ratchet test in `tests/search/` that walks every `data/catalogue.yaml` row
  with `source: literature`, extracts its DOI, and asserts each DOI appears either among the
  `KNOWN_CORPUS` anchor DOIs or an explicit, one-line-justified allowlist — cheap, mechanical, no
  schema-mapping risk, and would have caught this exact bug at #491 ingestion time. Fable's audit
  found **19 of 27** distinct catalogue literature DOIs currently have no `KNOWN_CORPUS` anchor at
  all — seed the allowlist with all of them (one line each) so the ratchet is green on landing, not
  a wall of new failures. Optionally SKETCH ONLY (do not build) a secondary warn-only runtime
  cross-reference layer for future consideration — full auto-generation of `KNOWN_CORPUS` from the
  catalogue remains explicitly NOT recommended (schema-mapping + circularity risk, per the answer
  already given to the user's "will KNOWN_CORPUS get auto-generated?" question); (4) register the
  Galilean symmetric-closure region in `data/empty_regions.jsonl` as literature-covered
  known-class-member, NOT novelty-bearing, with **corrected per-pair attribution**: Io-Ganymede,
  Europa-Ganymede, and Ganymede-Callisto — directly R-S-Table-1-enumerated (new anchors above);
  Io-Europa — covered as a body-subset of the existing Hernandez/Jones/Jesick IEG triple-cycler
  anchor (`body_set={Io,Europa,Ganymede}`), NOT by R-S itself; Europa-Callisto — covered as a
  body-subset of the existing Liang et al. CGE triple-cycler anchor (`body_set={Callisto,Ganymede,
  Europa}`, DOI 10.2514/1.G008387), NOT by R-S itself; note in the stamp that triple-cycler-subset
  coverage is weaker evidence than direct pair enumeration. **Stamp Io-Callisto SEPARATELY**,
  using Fable's suggested wording: "generated by our #563 reimplementation of R-S's ideal-model
  method; not in R-S Table 1's enumerated pair set and zero `russell-strange-2009-*` catalogue rows
  for the pair; no repeated-moon `KNOWN_CORPUS` anchor covers {Io,Callisto}; per #577 still NOT
  novelty-bearing" — avoid any phrasing implying "never discussed anywhere," since the gap is
  method-specific and DOI-specific, not a claim of total literature silence.
  **Explicitly out of scope:** any `catalogue.yaml` edit; running #576's pipeline stages;
  re-running the #576 search; the runtime warn-layer build (sketch only). **Run ALL ratchets** on
  any corpus/registry change (`uv run pytest tests/data tests/search -q`) per
  [[feedback_catalogue_edits_run_all_ratchets]]. **Recommended models:** Sonnet for execution — the
  judgment calls are now settled in the spec itself and the ratchets catch mechanical slips; no
  further Fable review needed before dispatch (this chain is 10-for-10 on real Fable catches, and
  this pass's corrections are the concrete fix, not a new open question).
  **✓ Resolved (2026-07-12), commit `575bd44`.** Executed the scope exactly as specified. (1) Added
  4 `CorpusAnchor`s to `KNOWN_CORPUS` (`russell-strange-2009-ganio`/`-ganeur`/`-gancal` for Jupiter,
  `-titenc` for Saturn), all DOI 10.2514/1.36610, `topology_label=frozenset({"repeated-moon"})`,
  `provenance="verified-against-source"`, grounded from the digest + the 30 catalogued
  `russell-strange-2009-*` member rows. Verified directly via `_candidate_anchors`: Io-Ganymede,
  Europa-Ganymede, Ganymede-Callisto now collide correctly (the first two ALSO already collided via
  the pre-existing Hernandez/Liang triple-cycler-subset anchors; Ganymede-Callisto via Liang); Io-
  Europa and Europa-Callisto still collide only via the Hernandez/Liang subset anchors, NOT R-S;
  Io-Callisto collides with nothing (`[]`), confirming the deliberate 3-separate-anchors-not-a-
  union design holds. (2) Extended `tests/search/test_literature_check.py`: 2 new parametrized
  R-S-2009 published-status tests (Ganymede-Io, Titan-Enceladus), an anchor-registration test, and
  an Io-Callisto non-collision test (both `_candidate_anchors` directly and end-to-end via
  `check_literature`). (3) Built `tests/search/test_corpus_doi_coverage.py` (3 tests): walks every
  `source: literature` catalogue row's `first_published.doi` + every `corroborating_sources[].doi`,
  asserts each is a `KNOWN_CORPUS` anchor DOI or in an explicit one-line-justified allowlist.
  Independently re-verified Fable's audit count before trusting it: 27 distinct literature DOIs, 19
  uncovered BEFORE this task's anchors landed (exact match), 18 after (R-S's DOI now anchored);
  allowlist seeded with all 18 (mission-description/infrastructure citations, Russell-Ocampo/
  McConaghy EM-cycler companion papers, 2 papers only reachable via anchor citation TEXT not DOI
  equality — flagged, not silently trusted). Also caught and flagged (allowlist comment, NOT fixed —
  `catalogue.yaml` edit is out of scope) a likely catalogue-side DOI typo: `10.2514/1.1011` appears
  in several `corroborating_sources` rows citing the same title as the already-anchored
  `10.2514/1.1909` Russell-Ocampo 2004 paper. (4) Registered the Galilean symmetric-closure
  territory in `data/empty_regions.jsonl`: one entry for the 5 literature-covered pairs
  (`jupiter-galilean-symmetric-closure-litcovered-known-class-member-578`) with the corrected
  per-pair attribution (Io-Ganymede/Europa-Ganymede/Ganymede-Callisto direct R-S-Table-1;
  Io-Europa/Europa-Callisto weaker triple-cycler-subset coverage, flagged as such), and Io-Callisto
  stamped SEPARATELY (`jupiter-galilean-io-callisto-symmetric-closure-corpus-gap-578`) with Fable's
  exact verbatim wording. **Unplanned fallout, fixed in the same commit:** the new narrower per-pair
  anchors caused `tests/data/test_citation_integrity.py::test_no_catalogue_citation_body_system_mis_
  citation` to newly fail — the two pre-existing, legitimately-broad R-S 2009 "family seed" catalogue
  rows (covering the whole paper, all 4 Galilean moons / all 6 Saturnian moons) now strong-link to
  the new anchors but no SINGLE anchor's body_set contains the full family claim (exactly the
  tradeoff of 3 separate anchors instead of 1 union). Fixed by generalizing the containment check to
  the UNION of all same-work strong-linked anchors (closes the Jovian case: 3 anchors' union = all 4
  Galilean moons) plus a narrow, explicitly-labelled family-seed exemption (author+system match still
  required) for the Saturnian row, which spans several R-S-documented science-target moons #578's
  scope does not add a per-pair anchor for. No `catalogue.yaml` edit. `uv run ruff check .` / `ruff
  format --check .` clean; full `uv run pytest tests/data tests/search -q` green (exit 0, no
  FAILED/ERROR). Pushed to `origin/main`.

- **#579** (P1, cheap, corpus-accuracy — do first per Fable's priority order) — fix the
  `literature_check.py` Antoniadou-Voyatzis/Libert anchor mislabeling + re-audit #287/#301.
  **Origin:** user supplied `antoniadou-voyatzis-2013-2-1-resonant-periodic-orbits-3d-planetary-
  systems-cmda-115-doi-10.1007-s10569-012-9457-4.pdf` (filed + digested
  `docs/notes/2026-07-12-digest-antoniadou-voyatzis-2013-2-1-resonant-3d-gtbp.md`); a Fable review
  of that digest independently verified the corpus-integrity claim and CORRECTED it (the digest's
  first-pass claim that TWO files were mislabeled was itself half wrong — `genome/known_corpus_3d.py`
  was already fixed in commit `0b1528f` / task #459; only ONE anchor is actually wrong). **Confirmed
  live bug:** `src/cyclerfinder/search/literature_check.py` (~line 1522-1546) has a `CorpusAnchor`
  `name="Antoniadou-Voyatzis spatial resonant periodic orbits in CR3BP (2018)"`,
  `authors=("Antoniadou", "Voyatzis")`, citing arXiv:1811.09442 with `doi=None` — but 1811.09442 is
  actually Antoniadou & **Libert** 2019 (MNRAS 483(3):2923, DOI 10.1093/mnras/sty3195; Voyatzis is
  not an author). **Load-bearing, not just cosmetic** (Fable's finding): the bad anchor's own
  comment claims coverage including "1:1" resonance and is cited as the anchor for "#287's 3D
  Braik-Ross (1,1) family extension (likely rediscovery target)"; #301's scope filters (lines 107,
  313 of `literature_check.py`) were built against this anchor's claimed "low-integer scope." But
  the paper it actually points to (Antoniadou-Libert 2019) covers MMRs 3/2, 2/1, 5/2, 3/1, 4/1, 5/1
  (NO 1:1) at μ=0.001 with no Earth-Moon ICs, and explicitly excludes asymmetric/isolated spatial
  families. A #287 (1,1) "likely rediscovery" verdict may therefore rest on coverage the cited paper
  does not actually have — a potential false NOT-novel. **CORRECTION (2026-07-12): the original
  candidate cited here ("Antoniadou & Voyatzis 2017, CeMDA 129") was itself an unverified citation
  inherited from an earlier digest and does NOT exist** — a CrossRef sweep of every real
  Antoniadou-Voyatzis joint publication found no such paper. The genuinely verified, topically-
  matching candidate is **Antoniadou, Voyatzis & Varvoglis (2014), "1/1 resonant periodic orbits in
  three dimensional planetary systems," Proceedings of the IAU 9:82-83, DOI
  10.1017/s1743921314007893** — a short (2-page) conference proceedings note, but specifically about
  **1:1 resonance**, the exact gap the mislabeled anchor's "1:1" claim needs and Antoniadou-Libert
  2019 does not cover. User is searching for this one. **Scope:** (1) relabel the anchor to
  Antoniadou & Libert 2019 with the correct DOI (or, if the anchor's ORIGINAL intent was genuinely
  the 1:1-specific IAU-proceedings paper above — acquire it if the user locates it, and cite that
  instead for the 1:1 claim specifically, keeping Antoniadou-Libert 2019 for the 3/2-5/1 MMRs it
  does cover — i.e. this may resolve into TWO anchors, not a single relabel); (2) re-derive the
  anchor's scope comment (MMR list, period band, the "1:1" claim) from
  whichever paper actually ends up cited; (3) re-audit the #287 Braik-Ross (1,1) "likely rediscovery"
  call and the #301 scope-filter logic against the corrected scope — per
  [[feedback_bugfix_invalidates_past_searches]], a corrected anchor obliges re-checking past verdicts
  that leaned on the wrong one; (4) the CORPUS_INDEX.md line and digest note for the 2013 paper have
  already been self-corrected (2026-07-12) to reflect the accurate one-anchor-not-two finding — no
  further digest edit needed. **Recommended model:** Sonnet (mechanical relabel + scope re-derivation
  behind existing ratchets; the #287/#301 re-audit needs a closer read but is not itself a
  discovery-judgment call). Run `uv run pytest tests/data tests/search -q` before commit.
  **✓ Resolved (2026-07-12), commit `194da55`.** Anchor relabeled to "Antoniadou & Libert spatial
  resonant periodic orbits in the RTBP (2019)", `authors=("Antoniadou","Libert")`,
  `doi="10.1093/mnras/sty3195"`, matching the twin `known_corpus_3d.py` anchor; false "1:1" scope
  claim removed, replaced with the real MMR list (3/2, 2/1, 5/2, 3/1, 4/1, 5/1). Investigated the
  newly-acquired Antoniadou-Voyatzis-Varvoglis 2014 IAU paper (1:1 co-orbital) as a possible fix but
  deliberately did NOT cite it — general (non-restricted) TBP model, a model-class mismatch with this
  restricted-problem anchor. Re-audit (`docs/notes/2026-07-12-579-anchor-fix-and-287-reaudit.md`):
  #299's "clean rediscovery, no novelty claims warranted" (1,1) verdict reopened as inconclusive — no
  anchor currently covers restricted-problem Earth-Moon 1:1 spatial resonance; no catalogue row was
  ever written back from that lineage, so nothing needs retracting. #301's separate 0/145 k=3-6
  sub-family verdict is unaffected (rests on other, correctly-cited anchors). `ruff`/`mypy` clean,
  full `pytest tests/data tests/search -q` green. Commit was co-swept with a concurrent #581 commit
  by a git-index race (documented honestly in that commit's message); content independently verified
  identical to the #579 agent's own work.

- **#580** (P2) — Richardson-1980 analytic third-order halo-seed generator. **Origin:** user
  supplied `richardson-1980-analytic-construction-periodic-orbits-collinear-points-celest-mech-22-
  241-doi-10.1007-BF01229511.pdf` (filed + digested
  `docs/notes/2026-07-12-digest-richardson-1980-collinear-halo-analytic.md`), independently assessed
  by Fable as worth pursuing after reading the actual codebase (not just the digest's framing). Fable
  finding: the codebase's ONLY existing 3D collinear-point seed generator,
  `search/cr3bp_seed_generator.py::lyapunov_seed_3d`, seeds the **vertical-Lyapunov/tulip branch**
  via a first-order linearized mode, NOT the halo branch — there is no halo-branch (Class I/II) seed
  generator anywhere in the codebase (halos exist only as literature anchors in
  `genome/known_corpus_3d.py`, Howell 1984 / Folta et al. 2015). The existing seeder's own code
  comments admit its corrector basin is "non-monotonic in (amplitude, μ)," compensated by an
  amplitude-ladder retry loop — evidence the linear seed is fragile at exactly the amplitudes where
  Richardson's third-order terms matter. This is a genuine branch-coverage gap, not a duplicate of a
  working path. **Scope:** add `richardson_halo_seed(system, point, amplitude_z, branch)` to
  `search/cr3bp_seed_generator.py` implementing Eq. 20a-c + the Appendix I coefficient recursions
  (a21..d32, s1/s2, l1/l2, k, λ — all closed-form in c2,c3,c4, which depend only on μ per Eq. 8a-b);
  golden test asserting the computed λ, k, c2-c4, a21…d32 against the digest's Table I transcription
  at Sun-Earth μ=3.04036e-6 for L1/L2/L3 (a genuinely sourced, non-circular golden — see digest §4);
  feed the generated IC to `correct_general_periodic_3d` and test that it converges to a non-planar
  halo orbit, ideally with a quantified comparison (fewer corrector iterations / larger reachable
  A_z before the ladder-retry logic engages) against the existing linear seed — Fable's explicit
  instruction: test that comparison, don't assume the payoff. **Framing:** this is infrastructure
  (a better/more-complete seed generator), not a discovery task — it will not itself produce a
  catalogue row; payoff is more robust halo/NRHO anchor generation for any future torus/heteroclinic
  lane that needs a halo starting point. **Recommended model:** Sonnet (spec-complete TDD behind a
  sourced golden). Not gated on user pre-approval (pure capability addition, no search/discovery
  claim), but flag the result when done.
  **✓ Resolved (2026-07-12), commit `af13056`.** Added `richardson_halo_coefficients` /
  `richardson_halo_ic` / `richardson_halo_seed` to `search/cr3bp_seed_generator.py` — closed-form,
  iteration-free `(mu, point, amplitude_z, branch) -> halo IC`. 18 golden tests validate λ, k,
  c2-c4, s1/s2, l1/l2, and a sample of the Appendix I coefficients against Richardson's Table I at
  Sun-Earth μ=3.04036e-6 for L1/L2/L3. **Along the way, caught and fixed 9 transcription errors in
  the digest's own Table I** (wrong signs, wrong exponents, digit transpositions, a cell swap),
  verified against the source PDF's text layer plus an independent third-party implementation
  (`jacobwilliams/Fortran-Astrodynamics-Toolkit`) — digest and `CORPUS_INDEX.md` corrected in place.
  **Seed-comparison result (Fable's mandated check, honest either way):** at amplitude_z=0.08
  (Earth-Moon L1), the existing linear vertical-Lyapunov seed fails to reach a non-planar solution
  while the Richardson seed converges cleanly in 5-9 iterations — a genuine, reproducible advantage
  in the regime the existing seed's own code comments document as fragile. `ruff`/`mypy` clean
  (one real mypy `no-any-return` finding fixed post-hoc — two `float()` casts on Legendre-coefficient
  expressions), full `pytest tests/data tests/search -q` green.

- **#581** (P2/P3, exploratory — flagged for user review before dispatch, NOT auto-fired) — a
  Gurfil-Kasdin-2002-style niching-GA search layer, staged and gated on a positive control.
  **Origin:** user supplied `gurfil-kasdin-2002-niching-genetic-algorithms-geocentric-orbits-3d-
  er3bp-cmame-191-doi-10.1016-S0045-7825(02)00481-4.pdf` (filed + digested
  `docs/notes/2026-07-12-digest-gurfil-kasdin-2002-er3bp-geocentric-orbits.md`) and the user's direct
  question "could genetic algorithms help us find more novel orbits?" — answered by an independent
  Fable review as **conditional yes, gated**, with this reasoning: the codebase already runs
  evolutionary optimization (scipy `differential_evolution`) in at least 6 modules
  (`search/optimize.py`, `lowthrust.py`, `maintain.py`, `mga_dsm_placement.py`,
  `global_precursor_engine.py`, `precursor_matcher.py`), so "GAs" per se are not new — but every one
  of those is single-optimum-per-run; family diversity currently comes entirely from OUTSIDE the
  optimizer (grid sweeps, direct commensurate-TOF construction, multi-start). That architecture is
  structurally blind to off-lattice, asymmetric, or unanticipated closures nobody enumerated a start
  for. Gurfil-Kasdin's **deterministic-crowding niching** is the specific mechanism (not "GA" in
  general) that removes this blindness — 14 co-existing families from ONE optimization run in the
  source paper, including ERO families (v0=0 periodic Earth-return) no resonance-lattice enumeration
  would generate. Fable also connected this to two live project facts: (a) the corrected
  Antoniadou-Libert 2019 anchor (`known_corpus_3d.py`) itself states it "does NOT cover
  asymmetric/spatial-isolated families, so 3D novelty stays open for those" — the same niche a
  niching GA structurally reaches; (b) per [[project_negative_results_registry]] ("empty is always
  conditional on the method"), a niching global search is a new method whose capability may
  legitimately reopen regions a prior grid sweep declared empty. Honest counterweights (Fable's own):
  a GA finds basins, not orbits — must feed the EXISTING correctors, not replace continuation; most
  hits will be rediscovery (novel hits are rare, per project memory); population×generations×
  propagation cost is a days-scale run (acceptable per project practice, per
  [[feedback_long_runs_acceptable]]); and the digest independently confirmed Gurfil-Kasdin's
  Earth-centered pulsating ER3BP frame matches NEITHER existing ER3BP frame module in the codebase
  (`core/er3bp.py` is Szebehely/Nechvile barycentric-true-anomaly; `core/er3bp_paper_frame.py` is
  Antoniadou-Libert non-pulsating barycentric-time) — a small frame-transform module is a
  prerequisite, not a drop-in reuse. **Scope, strictly staged (do not skip a stage):** (1) implement
  deterministic crowding as a thin layer over the existing DE infrastructure, plus the Earth-centered
  pulsating-frame variant (a geocentric offset applied to the existing Nechvile/true-anomaly
  machinery in `core/er3bp.py`, not a new frame from scratch); (2) **mandatory positive control**:
  reproduce a majority of Gurfil-Kasdin's 14 families from the digest's fully-sourced Tables 3/4 ICs
  in ONE niching run — if this fails, STOP, do not proceed to stage 3; (3) only if stage 2 passes,
  aim the niching layer at the two places Fable identifies as the plausible novelty payoff: (a)
  asymmetric/spatial-isolated 3D CR3BP resonant families (the literature-open niche per the corrected
  Antoniadou-Libert anchor); (b) bounded-drift quasi-cycler IC boxes (the #339/#473 frontier) — note
  Gurfil-Kasdin's own objective `1/[(rmax-rmin)²+1]` (Eq. 15) is exactly a bounded-oscillation
  criterion, the same notion the project's quasi-cycler gauntlet already validates against (per
  [[feedback_verify_gauntlet_with_positive_control]]). **Explicitly gated:** do not dispatch stage 3
  without a fresh Fable/user check on stage 2's actual reproduction results. **Recommended model:**
  Sonnet for stage 1 (mechanical layer + frame module) behind a sourced golden; Opus for judging
  stage 2's reproduction quality and whether stage 3 is warranted (a trust-bearing go/no-go call, not
  mechanical); Fable second-opinion on the stage-2 result before any stage-3 dispatch.
  **Stages 1-2 DONE (2026-07-12), commits `194da55` + `c5de132`; STAGE 3 STILL NOT AUTHORIZED —
  awaiting user decision.** Stage 1 built `core/er3bp_geocentric.py` (the paper's Earth-centered
  pulsating frame as a geocentric offset over the existing Nechvile/true-anomaly machinery, RHS
  parity to 6.7e-16 against an independently-transcribed Eq. 9-11) and `search/niching_ga.py`
  (Deterministic Crowding exactly per the paper's own p. 5687 pseudocode — a genuine multi-optima
  mechanism scipy's `differential_evolution` cannot provide). Stage 2's 12-set reproduction (one run
  per Table 2 optimization set, paper's own GA constants, no post-hoc tuning): **11/14 families
  recognizably reproduced** (A, B, C, D, E, G, H, J, K, L, M) against a pre-registered non-bit-exact
  criterion. Of the 3 misses: F and N are diagnosed as the GA settling on a DIFFERENT, equal-or-
  higher-fitness point in the same basin (published fitness vs. actual-converged-point fitness
  compared directly — basin competition, not a search failure); I is genuine under-convergence in
  the highest-dimensional (4 free variable) box at the paper's own stated generation budget. Along
  the way, caught real load-bearing source-interpretation bugs in the paper's own tables (Table 2/3
  vectors are INTERLEAVED `[x,x',y,y',z,z']` despite Eq. 13's stated `[x,y,z,x',y',z']` ordering;
  Table 4's km distances are AU-normalized; Family F's Table 4 rmax looks like a typesetting
  duplication of D/E's value — flagged, not asserted). Full details, results table, and honest
  caveats: `docs/notes/2026-07-12-581-niching-ga-stage2-positive-control.md`. **Coordinator's read**
  (not itself an authorization): this is a strong positive control by
  [[feedback_verify_gauntlet_with_positive_control]]'s standard — the mechanism demonstrably recovers
  multiple co-existing, topologically distinct families (including ERO/DEO cases a resonance-lattice
  or grid search would not generate) from single runs, which was #581's whole motivating claim. Stage
  3 (aiming the layer at asymmetric/spatial-isolated 3D resonant families or bounded-drift quasi-
  cycler boxes, per the original scope) looks worth authorizing on this evidence, but that decision is
  reserved for the user per the task's own gate.
  **STAGE 3 AUTHORIZED by user 2026-07-13 ("go on stage 3").** Split into two independent,
  parallel-dispatchable tracks per the original scope's two candidate targets — see #582 (track a,
  asymmetric/spatial-isolated 3D CR3BP resonant families) and #583 (track b, Sun-Earth bounded-drift
  quasi-cycler search). Both reuse `search/niching_ga.py::run_deterministic_crowding` unmodified
  (confirmed generic: fitness functions plug in as plain vector-in/scalar-out callables, same contract
  as scipy `differential_evolution`'s objective) — the only new work per track is a domain-specific
  fitness function plus its own downstream novelty/validation gate. Plan under Fable review before
  execution dispatch, per this session's established pattern for #571-#578.

- **#582** (P2/P3, exploratory, stage 3a of #581) — asymmetric/spatial-isolated 3D CR3BP resonant-
  family niching-GA search. **Origin:** #581 stage-3 target (a), authorized 2026-07-13.
  **Fable-corrected 2026-07-13** (original draft wrongly scoped this at Earth-Moon mu — see below).
  **System: mu=0.001** (a generic planetary system, NOT Earth-Moon — user decision 2026-07-13, choosing
  coherence with the open-niche claim and a direct positive control over Earth-Moon relevance, which
  would first need a #578-style corpus-anchoring prerequisite since Earth-Moon's relevant literature
  — Broucke/Hénon/Markellos-class periodic-orbit families — has no anchors in `KNOWN_CORPUS_3D` today).
  **Motivation:** the corrected Antoniadou-Libert 2019 anchor (`known_corpus_3d.py`, fixed by #579,
  explicitly annotated `mu=0.001`) explicitly states it does NOT cover asymmetric/spatial-isolated 3D
  resonant families — literature-open per the anchor's own scope statement, AT THIS mu. The existing
  3D isolated-family seed machinery (`search/er3bp_isolated_seeds.py`, #440, defaults to
  `mu=0.001, primary="Sun", secondary="planet"`) is circular-then-continued and explicitly SYMMETRIC
  (uses `correct_symmetric_fixed_jacobi`) — it structurally cannot GENERATE asymmetric seeds. But the
  codebase's asymmetric CORRECTOR already exists (`search/cr3bp_general_periodic_3d.py`, free vars
  `(x0..zdot0, T)`, min-norm Newton + independent closure check) and was already used for a planar
  asymmetric scan (#284/#343, `search/asymmetric_novel_scan_parallel.py`, full
  corrector→closure→topology→known-family→lit-footprint pipeline) — so the GA's genuinely NEW
  capability here is asymmetric/3D SEED GENERATION, not asymmetric correction; survivors must be
  routed through the EXISTING corrector, not a new one, and deduped against #284/#343's converged rows.
  **Scope:** (1) new fitness function, vector IC -> scalar, operating at mu=0.001 in `core/cr3bp.py`;
  free variables `(x0, z0, xdot0, ydot0, zdot0)` with `y0=0` (standard x-z-plane crossing convention)
  plus `T` as a genome variable bounded within roughly ±50% of the MMR's linearized resonant period
  `T0 = 2*pi / (a1**-1.5 - 1)` (#440's own formula, per its docstring); fitness combines a periodicity-
  defect term (`||state(T) - state(0)||^2` under propagation) with a soft target-Jacobi-band penalty
  (soft, not hard — an asymmetric branch may sit at a shifted C from the circular seed, matching #440's
  own note that converged members sit at a "nearby, finite-amplitude-shifted" period); death-penalty
  guards required: a minimum-period floor (T > T0/2, excluding degenerate near-equilibrium loops) and
  primary/secondary collision-exclusion radii throughout propagation (mirrors Gurfil-Kasdin's own
  Eq. 17-style collision constraint). (2) domain: the 5 INTERIOR MMRs already tabulated in
  `er3bp_isolated_seeds.py`'s `MMR_SEMI_MAJOR_AXES` ONLY for this pass — exterior bands are NOT
  tabulated anywhere and the module's own docstring flags the exterior 1:2 as a known family-selection
  trap; exterior MMRs are explicitly OUT OF SCOPE for #582 (a follow-up task, not this one, if wanted
  later). (3) one niching run per MMR band (mirrors #581 stage 2's one-run-per-set pattern).
  (4) MANDATORY pipeline before any claim, in order: GA cluster -> refine through
  `correct_general_periodic_3d`'s asymmetric mode with independent closure verification (a GA fitness
  peak is a basin indicator, not a converged orbit, per #581's own counterweight and
  [[feedback_orbit_closure_discipline]]) -> explicit SYMMETRY CLASSIFICATION (mirror-image /
  perpendicular-crossing test — a converged point ON a known symmetric orbit is trivially periodic too;
  without this test a symmetric orbit could be misreported as a novel asymmetric one) -> populate the
  survivor's `CandidateSignature` (`primary`, `body_set`, `topology_3d` k1/k2/k_z, `jacobi`) so
  `literature_check.py`'s 3D matcher actually engages (it only fires on same-`primary` anchors) -> clear
  the matcher (the machinery #579 fixed for #287/#299/#301) before any "novel" claim, per
  [[feedback_literature_novelty_check_baseline]] (not-found is necessary-not-sufficient) -> dedup
  against #284/#343's own converged/empty-region rows.
  **Positive control before novelty search:** first reproduce at least one already-known #440 circular
  member (any of the 5 tabulated MMRs) via the new fitness/bounds, converged through
  `correct_general_periodic_3d` and matched to the known member's `(x0, ydot0, T, C)` to a stated
  tolerance — not merely "the GA found a high-fitness point nearby."
  **Recommended model:** Sonnet for the fitness-function + pipeline-wiring build behind the
  positive-control golden; Opus/Fable for the literature-clearance adjudication on any survivor
  (matches #577/#579's pattern — novelty verdicts are trust-bearing, not mechanical).
  **Build + positive control DONE (2026-07-13), commit `ba60092`.** Built
  `search/isolated_3d_asymmetric_fitness.py` + `search/isolated_3d_asymmetric_pipeline.py` +
  `scripts/run_582_asymmetric_3d_niching_search.py`, 26 new tests. Positive control PASSES on MMR 3:2
  (the hardest of the 5 tabulated cases, deliberately chosen): GA+corrector recovers the known #440
  member to 0.6-4.1% error against 3-5% tolerances, independent Radau closure + symmetry-classifier
  self-check both confirm. MMR 5:1 near-miss documented honestly (threshold-edge corrector tolerance,
  not investigated). Found and fixed a real bug along the way: `correct_general_periodic_3d`'s default
  blind-Newton mode demonstrably DIVERGES on a GA-realistic seed (good fitness, wrong point) —
  `require_monotone_decrease=True` fixes it and is now pinned by a regression test. Flagged (not yet
  fixed) a real gap for whoever runs the full sweep: the literature-matcher reuses the Earth-Moon
  `primary="Earth"` label to reach the mu=0.001 Antoniadou-Libert anchor, which also pulls in unrelated
  physical Earth-Moon cycler papers into the candidate anchor pool — harmless for the positive control
  (no live `check_literature` search was run) but needs a proper generic-mu corpus label before the full
  novelty sweep trusts any verdict. Full results: `docs/notes/2026-07-13-582-asymmetric-3d-positive-
  control.md`. **NOT YET DISPATCHED: the full 5-MMR novelty sweep** (this was deliberately out of scope
  for the build dispatch — the coordinator launches and owns that run directly, per project convention
  for multi-hour compute).
  **Full 5-MMR novelty sweep DONE (2026-07-14), commit `b0225f4`.** Paper-scale budget (pop=200,
  gen=400) across all 5 MMRs, analyzed via the new `--mode analyze` (cluster + mandatory pipeline).
  **Result: 0/104 asymmetric** — every converged cluster (104 of 124 analyzed representatives, across
  all 5 MMRs) classifies SYMMETRIC; not one asymmetric member found. Diagnosed as very likely a
  search-box artifact, not a physical negative: `mmr_bounds()` is deliberately narrow and centered on
  the known symmetric seed, with tight absolute bounds (±0.05) on exactly the z0/xdot0/zdot0 components
  whose departure from zero breaks mirror symmetry — a narrow box centered on a strongly-attracting
  symmetric basin would produce this result whether or not a genuine asymmetric member exists nearby.
  Stamped as a method-conditional empty region in `data/empty_regions.jsonl`
  (`er3bp-isolated-3d-asymmetric-mu0.001-5mmr-582-2026-07-14`) per [[project_negative_results_registry]]
  — empty is never unconditional. Full writeup: `docs/notes/2026-07-14-582-5mmr-asymmetric-novelty-
  sweep-results.md`. **Natural next test (NOT run, a scope decision):** re-sweep with deliberately
  widened z0/xdot0/zdot0 bounds (the actual symmetry-breaking degrees of freedom), keeping
  x0/ydot0/T bounds as-is to anchor the resonance and limit MMR-drift risk.

- **#583** (P2/P3, exploratory, stage 3b of #581) — Sun-Earth ER3BP bounded-drift quasi-cycler search,
  widening #581's already-validated niching layer beyond Gurfil-Kasdin's own 12 optimization sets.
  **Origin:** #581 stage-3 target (b), authorized 2026-07-13. **Fable-corrected 2026-07-13** (original
  draft had a false-novel-factory gap and an infeasible positive control — both fixed below).
  **Motivation:** Gurfil-Kasdin's own fitness (`core/er3bp_geocentric.py::gurfil_kasdin_fitness`,
  Eq. 15 `1/[(rmax-rmin)²+1]`) is exactly a bounded-oscillation criterion — the same notion the
  project's quasi-cycler gauntlet validates against (per [[feedback_verify_gauntlet_with_positive_control]]
  — bounded drift-oscillation, not a strict floor). Stage 2 already validated the mechanism (11/14) and
  the frame (`core/er3bp_geocentric.py`, RHS parity 6.7e-16) on Gurfil-Kasdin's own 12 sets; this track
  reuses BOTH unmodified and only widens the search DOMAIN past those 12 sets' specific IC boxes.
  **Prerequisite (MANDATORY, do first, small #578-sized corpus task):** file corpus anchors BEFORE any
  widened-domain survivor is adjudicated — confirmed via grep that `literature_check.py`'s ~57 anchors
  currently contain ZERO for Gurfil-Kasdin 2002 itself and ZERO Sun-Earth co-orbital/DRO anchors (the
  only tadpole/horseshoe anchor is Pluto-Charon — wrong primary, will never match). Without these, the
  widened search will predictably rediscover the 14 G-K families, Earth quasi-satellites/horseshoes/
  tadpoles (Namouni 1999, Mikkola-Innanen, 3753 Cruithne, Kamo'oalewa), and Sun-Earth DROs/Hénon
  family-f distant satellites — ALL of which would wrongly "clear" `literature_check.py` as novel, the
  exact false-novel-factory trap [[feedback_literature_novelty_check_baseline]] exists to prevent. File:
  (a) a Gurfil-Kasdin-2002 anchor covering the 14 families' r-band footprints, (b) Sun-Earth co-orbital
  QS/HS/tadpole anchors, (c) Sun-Earth DRO/distant-satellite anchors. **Scope:** (1) no new fitness or
  frame module — reuse `gurfil_kasdin_fitness` and `run_deterministic_crowding` exactly as stage 2 left
  them; state explicitly (write the actual numbers into the run config, do not leave to the builder):
  which of the 6 state dims are free vs. fixed, whether `theta0` becomes a free/swept genome variable
  (recommended — stage 2 fixed it per-set, so leaving it fixed here would NOT genuinely explore beyond
  the 12 published sets) or stays a per-run constant, and that `n_rev` stays 1 for consistency with
  stage 2 unless explicitly widened. Note `escape_radius=0.5` (0.5 AU, ~50x Earth's Hill radius)
  already implicitly admits the full heliocentric co-orbital regime — a stated decision, not a builder
  surprise. (2) widen the free-variable bounds beyond the union of the 12 published optimization sets'
  boxes — the space BETWEEN and AROUND them, not a redo of stage 2's own boxes. (3) any bounded-drift
  survivor must clear, in order: the corpus-anchor prerequisite above (`literature_check.py` novelty),
  THEN a pre-registered bounded-vs-divergent drift criterion — closer template is
  `data/validation/v2_3d.py` (propagate a 6D IC, gate per-cycle drift) rather than `v2_moontour.py`
  (Earth-Moon Lambert-leg machinery, confirmed not a fit) — still needs adapting since quasi-orbits have
  no period: propagate N revolutions well past the 1-rev fitness window (N≈50-100, builder to justify
  the exact number), classify bounded (stationary geocentric r-band / recurrence) vs. divergent (secular
  rmax growth or escape), with the numeric thresholds written into the run config before any survivor is
  judged. The drift classifier itself needs its OWN positive control: G-K's 11 known-good stage-2
  families must classify bounded, and a known-escaping IC must classify divergent, before it judges
  anything new — per [[feedback_verify_gauntlet_with_positive_control]] (a wrong criterion choice has
  mis-classified real results before). Add a cheap `theta0`-robustness spot-check (re-test a "bounded"
  verdict at 2-3 other phases) per [[project_388_wall_energy_selective]]'s epoch-fragility lesson.
  **Positive control before widening (replaces the original "recover all 12 in one widened run" gate,
  which Fable flagged as likely infeasible — deterministic crowding's niche capacity has no reason to
  hold 12+ families in one pop=200 run when stage 2 needed 12 SEPARATE runs to get 11/14):** partition
  the widened domain into boxes each containing >=1 known family (reusing stage 2's per-set structure);
  before trusting anything novel found in a partition, first verify that partition's baseline run still
  recovers its neighboring known family/families under stage 2's own pre-registered non-bit-exact match
  criterion. **Recommended model:** Sonnet for the corpus-anchor filing + bounds-widening +
  drift-criterion build behind the reproduction regression; Opus/Fable for adjudicating any survivor's
  novelty and bounded-vs-divergent classification (trust-bearing per
  [[feedback_verify_gauntlet_with_positive_control]]'s own quasi-cycler precedent).
  **Dispatchable now — Fable's blocking findings resolved above.**
  **Build + positive control DONE (2026-07-13), commit `aafa244` (merged `378b1f1`).** Filed 3 corpus
  anchors (Gurfil-Kasdin 2002 itself, Sun-Earth co-orbital dynamics, Henon family-f — DOI gap on the
  last honestly flagged, not fabricated) and built `data/validation/er3bp_drift_classifier.py` (own
  positive control against stage 2's known families + a deliberately-escaping test IC) and
  `scripts/run_583_widened_bounded_drift_search.py` (7 partitions). First P1 run FAILED 0/6 — root-caused
  via independent Fable review: `LO_WIDE=0.002` AU re-admitted a physically trivial, strictly
  fitness-dominant deep-Hill-sphere basin (Eq. 15 rewards only annulus thinness, no periodicity content)
  that the whole population collapsed onto by generation ~62. Fixed by reverting `LO_WIDE` to the paper's
  own `LO_R` floor, plus an independent theta0 circular-distance bug in the match function. Corrected
  re-run: 1/6 clean match (family C, ic_dist 5e-4, feature ratios ~0.98) + 2 near-misses — a second
  independent Fable adjudication confirmed this satisfies the gate (validates the corrected machinery;
  remaining under-recovery is boundedness-fitness saturation + niche-capacity limits in a 6-family-wide
  box, not a bug — needs per-family/radial partitioning for any future full sweep, per
  [[project_negative_results_registry]]'s "empty/partial is method-conditional" framing). Full diagnosis,
  both Fable verdicts: `docs/notes/2026-07-13-583-corpus-anchors-and-drift-classifier.md`. **NOT YET
  DISPATCHED: the full 7-partition novelty sweep** (deliberately out of scope for this dispatch; needs
  the radial/per-family partitioning redesign first, per the positive control's own finding).
  **Process note:** two real subagent-backgrounding recoveries were needed on this task (see
  [[feedback_subagent_background_is_fatal]]) — the coordinator ended up running the actual positive-
  control GA passes directly rather than via a dispatched agent.

- **#584** (P2, cheap investigation) — 2 unexplained test failures discovered on this Mac (M3, native
  arm64 venv, Accelerate BLAS) that were NOT flagged by this session's prior CI runs on the previous
  Linux x86_64 environment. **Origin:** discovered 2026-07-13 running the full `tests/data tests/search
  tests/scripts` ratchet during #582's verification — confirmed to reproduce identically on a clean
  `main` checkout with zero #582 changes, so NOT introduced by #582:
  `tests/search/test_eggie_ballistic.py::test_gate_b_table4_vinf_reached_but_subsurface` (assertion
  `abs(8.313346200993646 - 8.38) = 0.0667 < 0.05` fails — a threshold-edge miss) and
  `tests/search/test_504_pluto_charon_kk_sweep.py::test_504_sweep_33` (topology (3,3) check fails on an
  otherwise-converged, crosscheck-clean sweep result). Per
  [[feedback_isolated_sweep_flips_suspect_artifact]], isolated threshold-edge failures in an otherwise-
  passing suite are more likely a numerical/platform artifact (BLAS backend switch: this venv now links
  Apple Accelerate, not whatever the prior Linux environment used) than a real regression — but this is
  unconfirmed, not yet investigated. Needs: (1) confirm whether CI (which may still run on Linux)
  reproduces these 2 failures or not — if CI stays green, this is very likely platform/BLAS-sensitivity
  local to this Mac, not a real bug; (2) if it does reproduce in CI too, these are real pre-existing
  regressions independent of the platform question and need proper triage.
  **RESOLVED (2026-07-14).** CI run 29323181595 (ubuntu-latest, commit `9e80969`, same tree these 2
  tests fail on locally) reports **"3310 passed, 23 skipped, 6 xfailed, 1 xpassed" — zero FAILED** —
  confirming both are local-Mac-only, not real regressions. Root-cause mechanism checked and is
  plausible for both: `test_gate_b_table4_vinf_reached_but_subsurface` fails a hard `< 0.05` km/s
  tolerance edge (actual 0.0667) on a converged differential-corrector V∞ output — exactly the kind of
  borderline threshold a different BLAS backend's rounding noise can tip either way.
  `test_504_sweep_33` fails an INTEGER winding-number classification
  (`topo.k1 == k1 and topo.k2 == k2` in `pluto_charon_kk_sweep.py`) derived from a corrector-found
  periodic orbit that is otherwise healthy (`stable_found=True`, `crosscheck_ok=True`,
  `crosscheck_dj=2.6e-12`) — a discretizing check on a continuous quantity is exactly what flips at a
  basin/tolerance boundary under tiny numerical perturbation, matching
  [[feedback_isolated_sweep_flips_suspect_artifact]]'s own established pattern. No code change made —
  CI is the authoritative gate and stays green; this is filed as a known, explained, Mac/Accelerate-
  local artifact for future sessions to skip re-investigating. **#584 STATUS: CLOSED.**

- **#585** (P2/P3, exploratory) — #582 follow-up: widen the symmetry-breaking bounds
  (`z0_abs`/`xdot0_abs`/`zdot0_abs`) and re-run the 5-MMR asymmetric novelty sweep. **Origin:** #582's
  5-MMR sweep (commit `b0225f4`) found 0/104 asymmetric members, diagnosed as very likely a search-box
  artifact — `mmr_bounds()` is deliberately narrow (±0.05 absolute) around exactly the state components
  whose departure from zero breaks mirror symmetry, calibrated only against the KNOWN symmetric seed.
  **Scope:** (1) pick new, wider `z0_abs`/`xdot0_abs`/`zdot0_abs` values — needs real justification, not
  an arbitrary multiplier: consider what amplitude a genuine asymmetric bifurcation off the symmetric
  family would plausibly need (literature on asymmetric CR3BP periodic families, e.g. Markellos-class
  results, may give a scale), balanced against #440's own documented exterior-1:2/neighboring-MMR
  drift risk — keep `x0_frac`/`ydot0_frac`/`t_frac` UNCHANGED to anchor the resonance. (2) MANDATORY
  positive control before any novelty claim: re-run MMR 3:2's positive control at the WIDENED bounds and
  confirm it still recovers the known symmetric member (proves widening didn't break basin containment)
  — if the known seed is no longer recoverable, the widening is too aggressive, tighten and retry before
  proceeding. (3) only then re-run all 5 MMRs at the new bounds via `--mode ga` + `--mode analyze` (reuse
  unmodified). (4) if still 0/104 asymmetric at meaningfully wider bounds, this becomes a much stronger,
  less method-confounded negative — update the `data/empty_regions.jsonl` stamp accordingly (a new
  entry, don't overwrite the #582 one, which stays valid for ITS bounds).
  **Recommended model:** Sonnet for the bounds-widening + positive-control-at-new-bounds build (spec-
  complete, gated by the mandatory re-check); the actual widened bounds VALUES are a real judgment call
  — consider a quick Fable/Opus sanity pass on the chosen widths before committing to a full sweep, since
  a wrong choice wastes the whole sweep's compute either way (too narrow: same null result for the same
  reason; too wide: drifts to an unrelated basin and the "match" becomes meaningless).
  **Fable design review DONE (2026-07-14): GO, resonance-scaled bounds ladder.** Quantitative finding:
  the current box is anisotropic by ~5-8x in exactly the symmetry-breaking directions — `xdot0_abs=0.05`
  is an eccentricity-proxy budget of only e≈0.03-0.044, while `ydot0_frac=0.35` already grants the
  SYMMETRIC direction e≈0.2 reach. **Concrete design:** replace the 3 flat absolutes with one
  resonance-scaled symmetry-breaking fraction `s`: `xdot0_abs = zdot0_abs = s * sqrt((1-mu)/a1)`
  (e/inclination-proxy scaling per-MMR, not a flat number), `z0_abs = max(0.05, s * a1)` (the `max()` so
  no rung ever shrinks below #582's already-stamped box). Run as a **2-rung ladder**: rung 1 `s=0.15`,
  rung 2 `s=0.30` (≈ the same eccentricity budget `ydot0_frac` already grants — the internally-principled
  stopping point; at 3:2 rung 2 gives ~7x the current box). Mandatory positive control (3:2) BEFORE each
  rung's sweep, per the original scope's own gate. **Add a drift-detection check to `--mode analyze`**:
  assert each refined member's semi-major axis (from its converged period) is nearest the TARGET MMR's
  `a1`, report drift-to-neighboring-MMR counts per rung — over-widening wastes GA capacity but is
  post-hoc detectable, not silently misleading. **Flagged gap, not in scope now:** `t_frac=0.5` (T within
  ±50% of T0) covers same-period pitchfork symmetry-breaking but structurally excludes period-doubled
  asymmetric branches (T≈2·T0) — note this in the eventual empty-region stamp as an out-of-box caveat,
  don't scope-creep it into this dispatch. Separate `empty_regions.jsonl` entry per rung actually run
  (don't overwrite #582's own entry, which stays valid for its bounds). Est. cost ~10-12 runs, 2.5-5h wall.
  **Dispatchable now.**
  **Build + per-rung positive control DONE (2026-07-14), commit `208db44` (merged `31970da`).** Built
  the `s`-parametrized `mmr_bounds()`, `mmr_a1_from_t0()` (inverse for drift-detection), and the
  drift-detection check in `--mode analyze`. **Rung s=0.15: PASSES** (0.3-3.7% error against tolerances,
  same known basin as #582's original box). **Rung s=0.30: FAILS** its own positive control — diagnosed
  as genuine basin competition (2 of 3 probe seeds land on a different near-periodic basin with
  comparable small-budget fitness, not simply "box too wide"). **s=0.30 is NOT certified for the full
  sweep**; only s=0.15 is confirmed safe. Full numeric bounds table + diagnosis:
  `docs/notes/2026-07-14-585-resonance-scaled-symmetry-breaking-bounds.md`.
  **Full 5-MMR sweep at s=0.15 DONE (2026-07-14), commit `83607a1`.** 78/78 converged clusters across
  all 5 MMRs and 125 analyzed representatives still classify SYMMETRIC — 0/78 asymmetric. The new
  drift-detection check confirms 0/78 drifted to a neighboring MMR, ruling out silent over-widening as
  a confound — a materially stronger negative than #582's original stamp. Stamped as a SEPARATE
  method-conditional empty region (`er3bp-isolated-3d-asymmetric-mu0.001-5mmr-585-s0p15-2026-07-14`,
  does not overwrite #582's own). s=0.30 remains uncertified (failed its own positive control). Full
  writeup: `docs/notes/2026-07-14-585-5mmr-sweep-s015-results.md`.

- **#586** (P2/P3, exploratory) — #583 follow-up: address the fitness-landscape/niching limitation the
  partitioning redesign exposed. **Origin:** #583's redesign (commit `f106519`) proved narrower,
  single-family partitioning does NOT by itself restore multi-family recovery — `gurfil_kasdin_fitness`
  (Eq. 15, boundedness-only) saturates near 1.0 across the whole bounded continuum once the trivial
  deep-Hill basin is removed, giving deterministic crowding no signal to discriminate between distinct
  genuine target families sharing one partition. Fable's own recommendation (from the redesign agent's
  report): "run multiple independent seeds per paper-band partition, and treat 'recovers ≥1 known family'
  as the realistic per-partition bar, not 'recovers every family on the first seed.'" **Scope (pick one,
  or scope both as sub-options for a reviewer to choose between):** (a) **multi-seed-per-partition**: for
  each of the 16 single-family/band partitions #583's redesign built, run N (e.g. 3-5) independent-seed
  niching-GA passes rather than one, treating "some seed recovers this partition's target family/families"
  as success — cheap to build (just a seed loop + aggregation), doesn't touch the fitness function, but
  doesn't fix the underlying discrimination problem, just plays the odds against it per-partition; (b)
  **family-discriminating fitness augmentation**: add a term to (or build a variant of)
  `gurfil_kasdin_fitness` that rewards proximity to a SPECIFIC target family's known IC/Jacobi/r-band
  signature, not just generic boundedness — turns each partition into more of a targeted search than a
  free exploration, which arguably changes what "novelty" even means for this track (a targeted search
  can't find something it wasn't aimed at) — needs explicit discussion of this tradeoff before building;
  (c) accept the current per-partition/single-seed limitation as this method's real operating envelope,
  document it plainly, and move on (a legitimate "we understand the limitation, it's not worth
  engineering around" outcome, per [[project_negative_results_registry]]'s "empty is method-conditional"
  framing — this ISN'T empty, it's "found the known family less reliably than hoped," a different
  situation from #582's clean 0/104). **Recommended approach:** get a Fable/Opus read on (a) vs (b) vs
  (c) before building anything — this is a real design-tradeoff decision (search fairness/cost vs.
  targeting vs. accepting the limitation), not a mechanical fix.
  **Fable design review DONE (2026-07-14): (a) + reframed cluster-everything harvesting + (c)'s
  documentation; REJECT (b).** Key reframe: per-partition family-match recovery is a POSITIVE-CONTROL
  concern, not the novelty deliverable — given a saturating boundedness objective, each GA run just
  samples the bounded continuum, and where it lands is basin-selection luck. Harvest the way #582's own
  `--mode analyze` already does: cluster the WHOLE final population into distinct high-fitness basins,
  refine/classify/drift-check/anchor-check EVERY representative, not just check against the partition's
  named target family. Under that harvesting, N=3 independent seeds per partition = 3 independent
  samples of basin space (a genuine novelty-sampling improvement, not just "playing the odds" on
  recovery). **Reject (b)** on its own logic: a target-proximity fitness term can't find what it wasn't
  aimed at, and targeted reproduction is already served by stage 2's narrow per-set boxes — building (b)
  would duplicate stage 2 while destroying this track's one distinctive property (unaimed exploration).
  Note: the `DEEP_HILL`/`BEYOND_HI_R` partitions' actual novelty payload is judged by the drift classifier
  + corpus anchors, NOT family-match, so the discrimination limitation barely touches them — (a) is
  sufficient to unblock the full sweep; the paper-band partitions become census/machinery-health checks.
  Do (c)'s documentation regardless (the limitation is real and empirically established — state it as a
  conditionality clause on any registry stamp: "recovery per partition is seed-conditional under Eq. 15 +
  deterministic crowding"). **Named future escalation, NOT scoped now:** if the N=3 sweep keeps surfacing
  recurring unmatched bounded basins that get lost to niche collapse, escalate to a MAP-Elites-style
  behavior-space quality-diversity archive (binned on rmin/rmax/Jacobi/period, all already computed by
  the pipeline) — this preserves novelty semantics the way (b) can't, but is deferred pending evidence
  from the cheap N=3 pass first (moderate build cost against low expected yield in this already
  heavily-charted Sun-Earth bounded domain).
  **Scope for this dispatch:** N=3 independent-seed niching-GA runs per each of #583's 16 single-family/
  band partitions, harvested via #582-style cluster-everything (not narrow family-match only), documented
  per-partition recovery + any unmatched-bounded-basin candidates (which still need the live
  `check_literature()` novelty gate before any claim, per [[feedback_literature_novelty_check_baseline]]).
  **Dispatchable now.**
  **Build + full 16-partition x 3-seed sweep DONE (2026-07-14), commit `d260811` (build, merged
  `a5f9ebc`) + `88a6d9e` (sweep results).** 6/14 single-family partitions (C, D, E, G, H, M) recovered
  their own target family in >=1 seed — confirms the corrected pipeline works when it works; the other 8
  did not in ANY of 3 seeds, matching the design review's own "seed-conditional, not guaranteed"
  framing. **264 unmatched-bounded candidates** surfaced across the sweep — raw, literature-matcher-
  engaged-only, NO novelty claim (no dedup, no live `check_literature()`, no Opus/Fable adjudication run
  yet). Full per-partition table + honest caveats:
  `docs/notes/2026-07-14-583-16-partition-3-seed-sweep-results.md`. **NOT YET DISPATCHED: dedup + live
  literature-check + Opus/Fable adjudication of the 264 candidates** — a substantial, separate task (see
  #588).
  **Build + small-scale validation DONE (2026-07-14), commit `d260811` (merged `a5f9ebc`).** Added
  `--n-seeds` (independent RNG seed per run) and a harvest step clustering EVERY seed's final population
  into distinct high-fitness basins (not just checking the partition's own target family), running each
  representative through the drift classifier + a match check against all 14 published families.
  Small-scale validation (partition `C`, 2 seeds): both seeds recovered family `C`; the harvest also
  surfaced **14 unmatched-bounded candidates** — bounded, non-matching any of the 14 known families,
  literature-matcher-engaged — flagged for the eventual full-sweep adjudication pass, explicitly NOT
  adjudicated or claimed novel here (no live `check_literature()` run). Full design + validation:
  `docs/notes/2026-07-13-583-corpus-anchors-and-drift-classifier.md` (Part 5 addendum). **NOT YET
  DISPATCHED: the full 16-partition × N=3-seed sweep** (coordinator-owned multi-hour run).

- **#587** (P2, website/catalogue-data) — populate real CR3BP identity fields (Jacobi constant, period,
  stability index) for the 6 #569 Uranian symmetric-closure quasi-cycler rows (`umbriel-oberon`,
  `titania-oberon`, `ariel-oberon`, `umbriel-titania`, `ariel-titania`, `ariel-umbriel` — includes #312
  itself). **Origin:** the cyclers.space website's orbit-view component was found to silently
  mis-render all 6 rows as a near-blank Sun+Earth frame (`cycler_class="multi-arc"` with no
  `trajectory.segments`/`orbit_elements` populated, and `bodies` using names — "Uranus", "Ariel",
  "Oberon" — that don't match the website's heliocentric-planet geometry lookup at all). Quick fix
  already applied website-side (`renderClassOf()` now checks `model_assumption === "cr3bp"` before
  `cycler_class`, so these rows show the honest "not renderable from current data" placeholder instead
  of a misleading blank frame) — but that placeholder's own text claims "The CR3BP identity (Jacobi
  constant, period, stability) is tabulated above," which is currently FALSE for all 6 rows (no
  `orbit_elements.cr3bp` block exists in any of them). **Scope:** populate `orbit_elements.cr3bp`
  (`jacobi_constant`, `period_nd`, `stability_index`) for these 6 catalogue rows from the #563/#566
  discovery/gauntlet data that already produced them (should be readily available — these ARE CR3BP
  periodic-orbit search results, the identity data exists somewhere in the #558-#569 task chain's
  artifacts, it just never made it into the catalogue row). Also flagged, lower priority: a genuine
  data inconsistency where the same physical body is coded inconsistently across rows — e.g. "U" (used
  by flyby/tour rows like Heaton-Longuski 2003) vs "Uranus" (used by the #569 writeback rows) — worth
  normalizing to one canonical code convention. **Separately, NOT part of this task:** the website's
  multi-arc trajectories (297 rows) have ZERO real per-segment (a,e) data and zero sampled-trajectory
  exports, so they get no 2D or 3D time-animation at all (a static positions-only SVG) — fixing that
  needs a genuine trajectory-sampling exporter (the "Phase-C exporter" already anticipated in
  `OrbitView.astro`'s own comments), a separate, larger scoping task if wanted.
  **CORRECTED AND DONE (2026-07-14).** Investigating before executing found the original premise
  was wrong: per this project's own `catalog.py::CatalogueEntry.fully_defined`, the `cr3bp` identity
  triple applies ONLY to `cycler_class="non-keplerian"` rows — all 6 rows are `cycler_class="multi-
  arc"` (a Uranus-centered tour across two distinct moon-pair restricted-3-body systems), which has
  no single Jacobi constant/period/stability triple to tabulate at all (same precedent as this
  catalogue's existing Jovian/Saturnian Tisserand-tour rows: "Jacobi constant is not conserved in
  this model"). The real gap was `invariants`/`legs[].tof_days` (the schema-correct multi-arc
  completeness fields) plus a website placeholder that falsely claimed CR3BP identity was tabulated
  for these rows. Fixed: populated `transit_times_days`/`legs[].tof_days` from each row's own already-
  cited source data (`data/silver_327_verified.jsonl` / `data/gauntlet_566_five_representatives.jsonl`
  `tof_days`); added honest `data_gaps` entries for `aphelion_ratio`/`turn_ratio`/`orbit_elements.cr3bp`
  explaining why each is inapplicable/unknown (not fabricated); re-ran the pre-existing (idempotent)
  `scripts/backfill_invariants.py` which incidentally caught 45 other multi-arc rows never backfilled
  since it was added; fixed `cyclers.space/src/components/OrbitView.astro`'s placeholder to
  distinguish genuine single-orbit CR3BP rows (real `cr3bp` block, text unchanged) from multi-arc
  tours mis-routed to the same render class (new, accurate text, no false tabulation claim). Full
  writeup: `docs/notes/2026-07-14-587-multiarc-schema-correction.md`. The body-code-inconsistency
  ("U" vs "Uranus") flag from the original spec remains unaddressed — genuinely lower priority,
  left for a future pass. **#587 STATUS: CLOSED.**

- **#588** (P2/P3, exploratory) — dedup + live-literature-check + Opus/Fable adjudication of #586's
  264-candidate unmatched-bounded pool from the full 16-partition x 3-seed sweep. **Origin:** commit
  `88a6d9e`, `docs/notes/2026-07-14-583-16-partition-3-seed-sweep-results.md`. **Scope:** (1) build a
  dedup pass — 264 raw candidates across 48 independent GA runs almost certainly contain many
  near-identical repeats of the same underlying basin (e.g. `DEEP_HILL`'s 31, `BEYOND_HI_R`'s 26, `I`'s
  48 are far too high to represent that many genuinely distinct orbits); cluster by IC/Jacobi/period
  proximity across the WHOLE pool, not just within one seed/partition (needs its own distance-threshold
  design, likely reusing the same bounds-normalized metric `cluster_representatives`/harvest already
  use, but applied globally). (2) run the ACTUAL live `check_literature()` search (not just the offline
  `KNOWN_CORPUS` matcher engagement check already done) on whatever survives dedup — per
  [[feedback_literature_novelty_check_baseline]], not-found is necessary-not-sufficient, this is the
  mandatory gate before any novelty claim. (3) Opus/Fable adjudication of whatever clears both gates,
  matching the #564/#565/#577 pattern (Uranian/Galilean symmetric-closure survivor-list precedent).
  **Recommended model:** Sonnet for the dedup build (mechanical, well-specified once the distance
  metric is chosen); Opus/Fable for the final adjudication (trust-bearing, not mechanical).
  **Dedup DONE (2026-07-14).** `scripts/dedup_588_candidate_pool.py` re-characterizes every raw
  candidate via `characterize()` (common physical space, not raw `ga_genome`) and greedily clusters
  fitness-ranked, never across `candidate_type`, threshold=0.10 relative deviation (empirically
  justified via a chaining-sensitivity scan). 264 raw -> **45 distinct clusters**
  (`data/found/583_widened_search/deduped_candidates.json`, full provenance retained). Design +
  results: `docs/notes/2026-07-14-588-candidate-pool-dedup.md`.
  **Live literature-check DONE (2026-07-14).** Cross-referenced all 45 clusters against the 14 known
  Gurfil-Kasdin families' own r-bands: 25/45 overlap a same-type known family's range (plausibly other
  points on an already-published continuous family), leaving **20 unmatched**. Live search (WebSearch +
  independent CrossRef verification, per [[feedback_ground_citations_against_content]]) surfaced two
  Gurfil-Kasdin companion papers not previously in corpus, both acquired (user-supplied PDF) and fully
  read: (1) SPIE 2003 (4854:251-261, DOI 10.1117/12.459820) — confirmed via exact numeric match
  (Family F rmin=224,900 km, Family C rmax~1.13M km) to be the SAME 14-family census already anchored,
  not new data; (2) CMAME 2002 out-of-ecliptic paper (191:2141-2158, DOI
  10.1016/S0045-7825(01)00380-2) — genuinely distinct (different objective: maximize z², circular not
  elliptic RTBP, search box out to 1-5 AU from Earth), but reports only isolated max-displacement
  points at a scale (33-56M km z-excursion, several-AU Earth distance) an order of magnitude beyond our
  largest unmatched cluster (40.8M km) — does not structurally cover the gap. Now added as a new
  `CorpusAnchor` in `literature_check.py` for future large-excursion searches (see #589). Full trail:
  `docs/notes/2026-07-14-588-gurfil-companion-papers-digest.md`. **20 unmatched clusters proceed to
  Opus/Fable adjudication as originally scoped — NOT YET DISPATCHED.**
  **Opus + Fable adjudication DONE (2026-07-14).** Full writeup:
  `docs/notes/2026-07-14-588-final-adjudication.md`. **14/20 confidently not novel:** 8 large planar
  DEO clusters are Henon ẏ₀≈-2x₀ curve extensions (same curve as Families A/B/F); 6 small 3D DRO
  clusters (80k-560k km) are not-novel because the Sun-Earth ER3BP model itself is invalid below
  ~600k km (the Moon, omitted from the model, dominates there, not the Sun) plus continuum-membership
  logic — NOT because of the `v0/vcirc≈1` heuristic Opus initially used, which Fable's adversarial
  check found would ALSO dismiss real published Families F (0.91) and K (1.30) if applied to them
  (logged as a 3rd instance of [[feedback_verify_metric_semantics_before_ranking]] — calibrate any
  NEW ad-hoc triage metric against the known-good reference set, not just pre-existing pipeline
  fields). **6/20 warrant a bounded follow-up before dismissal is final** (not a novelty claim): 40,
  42, 43 (genuinely 3D, near Family J's band, but each a single unreplicated GA hit with a growing
  5yr r-band) plus 24, 25, 39 (radius-band-adjacent to B/I/C; 39 was originally dismissed on a
  factually wrong band-match claim, corrected by Fable). Follow-up registered as #590. **Scope
  reminder: none of these 20 (nor the 14 known G-K families) are catalogue-admissible under this
  project's cycler/quasi_cycler/precursor_mga/mga_tour scope — this was a literature-completeness
  question, not a catalogue-discovery one; #588 closes with zero catalogue rows changed.**
  **#588 STATUS: CLOSED.**

- **#590** (P3, cheap, bounded follow-up) — settle whether Sun-Earth ER3BP clusters 40, 42, 43 (+
  optionally 24, 25, 39) from #588's adjudication are genuine unsampled members of Family J's curve
  or single-hit GA artifacts, before any of them could even be considered for a novelty claim.
  **Origin:** `docs/notes/2026-07-14-588-final-adjudication.md`. **Scope:** (1) a targeted literature
  search against the classical quasi-satellite-orbit corpus Fable identified as the actually-relevant
  prior art for 4-13M km inclined 3D Sun-Earth structures — Mikkola et al. 2006 (MNRAS), Lidov &
  Vashkov'yak 1993/1994, Sidorenko et al. 2014, Pousse-Robutel-Vienne 2016/2017 (CMDA) — with every
  DOI independently CrossRef-verified before citing, per
  [[feedback_ground_citations_against_content]] (none of Opus's 3 cited URLs from the #588 pass were
  ever grounded in a note and must not be inherited); (2) a cheap numerical continuation chain
  (predictor-corrector in orbit size, no GA sweep) connecting 40/42/43 to Family J to settle
  same-curve vs. distinct-branch directly; (3) a longer-horizon (~100-rev) drift-classifier re-check
  for 40 and 42 specifically, whose 5yr r-bands already show the largest growth trend. Cheap,
  bounded, no new search infrastructure.
  **Recommended model:** Sonnet (mechanical once the 3 steps are scoped; no new judgment calls).
  ✓ **DONE (2026-07-14, this session):** `docs/notes/2026-07-14-590-followup-results.md`.
  (1) 3 QSO papers CrossRef- and full-text-verified and added as `CorpusAnchor` entries
  (Mikkola et al. 2006 MNRAS, Sidorenko et al. 2014 CMDA, Pousse-Robutel-Vienne 2017 CMDA);
  Lidov & Vashkov'yak 1993/1994a/b honestly left un-anchored (CrossRef does not resolve
  them, pre-DOI-era Russian-journal papers) — none of the 3 anchors tabulate a direct
  comparison point, so the literature gap is characterized but NOT closed. (2) cheap
  straight-line IC-interpolation connectivity heuristic (explicitly caveated as weak, not a
  real continuation): cluster 43->J cleanest (no escape en route), 42->J roughest (genuine
  escape region mid-path), 40->J intermediate; cluster 24's nominal anchor (Family I) turned
  out to be one of the 3 Gurfil-Kasdin families that does NOT reproduce/stay bounded under
  this repo's own #581 pipeline, so no connectivity verdict is possible for 24 at all — an
  honest finding, not papered over. (3) longer-horizon (300yr) drift re-check: clusters 40
  and 42 both DO eventually escape (~213yr, ~220yr) — but so does Family J's own literal
  published IC, at a much shorter ~30yr (independently cross-checked against Table 4 to
  confirm no IC-conversion bug); cluster 43 (bonus check) stays bounded to at least 1000yr.
  "Eventually escapes" is the literature-predicted (Mikkola et al. 2006) generic behavior of
  inclined QS/DRO motion, not a discriminator between real family members and GA artifacts —
  reframes rather than resolves the 6-cluster question. Zero catalogue rows changed, no
  novelty claimed; left for a future adjudication pass. **#590 STATUS: CLOSED.**

- **#589** (P3, exploratory, new capability) — reproduce/extend the Gurfil-Kasdin 2002 out-of-ecliptic
  z²-maximization approach (CMAME 191:2141-2158, DOI 10.1016/S0045-7825(01)00380-2, digested during
  #588's literature check) as a genuinely new search regime, distinct from the bounded-DRO/DEO family
  census #581-#586 have been mining. **Why this is a different niche, not a #583 rerun:** that paper's
  objective is maximize normal (z) displacement (not periodicity/boundedness), circular (not elliptic)
  Sun-Earth RTBP, search box out to 1-5 AU from Earth with v0 up to 15.4 km/s — reaching z-excursions
  of 0.223-0.374 AU (33-56M km) and multi-AU Earth distances, an order of magnitude beyond anything in
  our current unmatched-candidate pool (#588's largest cluster: 40.8M km). The paper itself reports
  only 3 qualitative types (I/II/III, by in-plane vs. out-of-plane frequency ratio) from a handful of
  GA runs with NO systematic multi-seed/niching harvest and NO literature check of its own against
  later work — a plausible genuine discovery gap given #581-#586's own niching-GA + dedup + live-check
  machinery already exists and could be pointed at this same objective function/search box instead of
  the bounded-family one. **Scope (not yet designed):** (1) reproduce the paper's own 3 named
  trajectories (low-energy Type II z=0.223 AU/rmax=2AU; high-energy Type III z=0.374 AU/rmax=3AU) as a
  positive control using the existing niching-GA + `characterize()` infrastructure with a new
  z²-maximizing fitness function; (2) if reproduction succeeds, widen via multi-seed cluster-everything
  harvest (reusing #586's pattern) across the same 1-5 AU search box the paper explored but a niching
  GA it never ran; (3) dedup + live-literature-check + Opus/Fable adjudication of whatever survives,
  same pipeline as #588. **Flagged for user review before dispatch, not auto-fired** — this is a new
  fitness-function/objective build, not a mechanical follow-up, and should get a Fable design read on
  the objective-function choice and search-box scoping before any code is written (matching #581's own
  gating precedent for introducing a new niching-GA capability).
  **Fable design review (2026-07-14): NO-GO on the 3-stage plan as scoped / SCOPE-DOWN** — see
  `docs/notes/2026-07-14-589-fable-design-review.md` (source paper read in full). Two premises above
  are factually wrong against the paper's own text: (a) the paper's titular method IS a
  deterministic-crowding niching GA (pop 150 x 300 gens per set) — the un-run part is only the
  multi-seed/cluster-everything harvest; (b) "search box out to 1-5 AU" misreads Table 2 — the IC box
  is r0 <= 150,000 km near-Earth departures (design cases: 200-km LEO / GEO), and 1-5 AU is `rmax`, a
  PATH constraint, which is also the paper's real structure-carrying axis (rmax sweep -> Type
  I/II/III), not a niche map. Further: z²-max optima sit on the IC-box boundary (the paper's own
  Fig. 2 shows corner convergence, the only "niches" being the trivial +-z mirror pair), so a niching
  harvest is a category error here; the paper prints NO optimal ICs, so the positive control reduces
  to matching 2 scalars that two-body v_inf/plane-change arithmetic already reproduces to ~90-95%;
  Type I/II/III is exhaustive by construction (the 3 orderings of one frequency ratio — no 4th type
  exists to find); the class is non-returning by design so it fits NONE of the 4 catalogue classes;
  and the proposed reuse is structurally broken — the 0.5 AU terminal escape event hard-wired into
  `gurfil_kasdin_fitness`/`characterize()`/`classify_bounded_drift` kills every 1-5 AU target
  trajectory on sight (per #590, boundedness gates are the wrong lens even for the bounded families,
  a fortiori for an admittedly-drifting class). **Recommended instead:** a half-day analytic
  characterization note (conic accounting of the published scalars + taxonomy-exhaustiveness +
  Table 2 encoding derivation — the review note already contains the core numbers), then file the
  regime as empty-of-catalogue-relevant-structure by construction; any future reproduction gated on
  an unexplained analytic residual and spec'd as its own small circular-CR3BP (e=0, no theta0) build
  with a new frequency-ratio characterization layer, never judged by the bounded-drift stack. A
  salvageable reframing (separate task, not #589): max out-of-ecliptic excursion SUBJECT TO
  boundedness across the validated #581-#586 families — defined success predicate, honest infra
  reuse, catalogue-relevant output.
  **Not dispatched; do not dispatch as scoped.**

- **#591** (P3, cheap, analysis-only) — Fable's salvageable reframe of #589: rank already-validated
  BOUNDED 3D structures by out-of-ecliptic excursion (max |z|), rather than searching for new
  unbounded z-maximizing trajectories. **Origin:** `docs/notes/2026-07-14-589-fable-design-review.md`
  and the user's own follow-up question about whether the excursion idea was worth keeping — this is
  the narrower version Fable said IS worth it: "defined success predicate, honest infra reuse,
  catalogue-relevant output." **Scope:** (1) enumerate the pool of already-validated bounded 3D
  structures — the 5 known Gurfil-Kasdin 3D families with nonzero z (J, K, L, M, N from `TABLE34` in
  `scripts/run_581_gurfil_reproduction.py`) plus any #588 deduped 3D clusters independently confirmed
  bounded at a real horizon (cluster 43 in particular: confirmed bounded to 1000+ years by #590,
  z0≈0.0316 AU ≈ 4.7M km — a strong candidate); (2) for each, propagate (reusing
  `core/er3bp_geocentric.py`'s `er3bp_eom`/`geocentric_to_barycentric` machinery, same convention as
  `characterize()`) and extract max |z| over the orbit (not just the IC z0 — `characterize()` computes
  `rmin`/`rmax` from the full radial distance but not z alone, so this needs a small addition, not a
  new search); (3) rank by max |z| (km and AU), producing a table + a short note — this is the
  "off-ecliptic observatory utility" figure of merit the 2002 paper was chasing, but computed
  honestly over orbits already known to return/stay bounded, with zero new GA search. **Non-goals:**
  no new fitness function, no niching, no novelty claim — purely a derived-metric ranking over
  existing validated data. Script naming: `analyze_591_*.py` (not `run_*.py` — this is post-hoc
  analysis over already-searched data, not a new search, so it correctly sits outside the
  `tests/scripts/test_scripts_call_preflight.py` AST ratchet, matching the `dedup_588_*.py`
  precedent).
  **Recommended model:** Sonnet-tier work, or direct (mechanical, well-specified, low compute — a
  handful of ~5yr propagations).
  **DONE (2026-07-14).** `docs/notes/2026-07-14-591-bounded-ecliptic-excursion-results.md`. Result:
  ranked bounded pool is cluster 43 (6.60M km) > Family J (5.71M km) > Family M (591.8k km) > Family K
  (145.8k km); L/N excluded from ranking as confirmed-DEO (not bounded by type, both verified to
  escape within the 5yr window itself). **Cluster 43 — from #588's "not novel" pool — beats the
  published state of the art (Family J) on out-of-ecliptic excursion by ~16%, while also being (per
  #590) far more durably bounded (1000+ yr vs. J's own ~30yr).** Not a novelty claim (still plausibly
  a Family-J-curve extension per #588/#590), but a genuine, citable figure-of-merit result over
  already-validated data. No new search, no catalogue rows changed.
  **#591 STATUS: CLOSED.**

- **#602** (P2, cheap, bounded follow-up) — the rigorous differential-correction continuation `#590`
  explicitly deferred: `#590`'s own connectivity check between Sun-Earth ER3BP clusters 40/42/43 and
  Family J was an explicitly-caveated "cheap straight-line IC-interpolation heuristic... NOT a real
  continuation." Settle definitively whether cluster 43 (the `#588`/`#591` standout: beats Family J
  by ~16% on out-of-ecliptic excursion, 1000+yr bounded vs. Family J's own ~30yr) sits on the SAME
  continuous branch as Family J (pseudo-arclength continuation connects them → rediscovery, not
  novel, but a clean close) or is a genuinely DISTINCT branch (continuation fails/forks before
  reaching Family J's IC → a real candidate for a literature-novelty check + catalogue
  consideration). Use `search/pseudo_arclength.py` + `search/cr3bp_continuation.py` +
  `search/deflated_newton.py` (all already exist) rather than the weak heuristic. Family J's IC:
  `TABLE34["J"]` in `scripts/run_581_gurfil_reproduction.py`. Cluster 43's IC: locate in the
  `#588`/`#590`/`#591` docs/notes or the `dedup_588_candidate_pool.py` output. No blocker.
  **✓ DONE (2026-07-15) — genuinely INCONCLUSIVE, not a forced verdict.**
  `docs/notes/2026-07-15-602-cluster43-familyj-continuation-results.md`. Cluster 43's IC traced to
  `data/found/583_widened_search/deduped_candidates.json` (`cluster_id: 43`). Built a Sun-Earth-
  specific positive control (the existing `cr3bp_continuation.py` test is Earth-Moon) and found the
  actual reason a hard verdict is unreachable: **neither Family J nor cluster 43 is close to an
  exact year-periodic orbit** (both stall at nonzero residual / pin to the trust-region boundary
  under an independent robust solver, confirmed not a finite-difference artifact via the analytic
  STM) — these are quasi-periodic/metastable structures (consistent with Mikkola et al. 2006's own
  characterization of this regime), so "continuous curve of exact periodic orbits" is not the right
  formalization of "family" here; a genuine model-mismatch finding, not a tooling gap. Ran anyway: a
  damped deflated-Newton homotopy walk (21 steps) gave a patchy correction signal but a CLEAN
  boundedness result (unbounded only near Family J's own end, matching its known ~30yr escape;
  bounded continuously from t=0.25 to cluster 43 at a full 50yr horizon); a secondary check on the
  already-validated planar/circular retrograde-satellite backbone (after catching+fixing a real
  geocentric/barycentric coordinate-offset bug via a Jacobi-constant sanity check) walked cleanly
  past cluster 43's radius with zero folds/bifurcations (`stop_reason="max_steps"`, not a topology
  break). **Verdict: no evidence found FOR a distinct branch anywhere, and real (non-heuristic)
  evidence — the fold-free backbone continuation + longer-horizon boundedness — continues to favor
  "same family, unsampled point" over #588's "not novel" framing, but this is preponderance-of-
  evidence, not proof.** Recommended next step if ever revisited: a proper quasi-periodic/rotation-
  number continuation corrector (materially new capability, out of this task's bounded scope). No
  catalogue rows changed, no novelty claimed. Thread (`#588`→`#590`→`#591`→`#602`) closed at this
  evidence tier.

- **#603 ✓ DECLINED (2026-07-15, user decision)** — Sun-Neptune transient-capture `quasi_cycler` screen, extending `#535`/`#557`'s
  Hill-sphere-return detector to a system genuinely unswept by this method (per a discovery-strategy
  planning pass 2026-07-15). Motivation: `#557`'s kill mechanism at Jupiter was the real eccentricity
  (e=0.0489) collapsing 15/16 idealized survivors; Neptune's e=0.00859 (Standish & Williams Table 1,
  `core/constants.py`) is ~5.7x weaker, so the same collapse mode may not fire. Machinery proven
  twice (Earth `#535`, Jupiter `#557`). **NOT the same as `#341`/`#492`** (Neptune Triton-Proteus
  repeated-moon-flyby cycler searches — a completely different method/question; this is about
  heliocentric small bodies transiently captured BY Neptune, not spacecraft/moon encounters).
  **Real scope risk found before dispatch, why this needs its own plan+decision (like `#557` got)
  rather than direct execution**: Neptune's orbital period is ~165 yr (vs Jupiter's 11.86 yr), so
  the SAME window-rescaling logic that gave Jupiter a 119-178yr admission window would give Neptune
  something like ~1650-2475 YEARS — a far more extreme departure from the catalogue's
  mission-relevance `quasi_cycler` definition than Jupiter's already was. Plan must also check
  whether a real, literature-documented anchor object even exists (Jupiter had Gehrels 3/Oterma;
  does any known Centaur/scattered-disk object have a DOCUMENTED temporary Neptune capture episode
  the way Koon et al. 2001 documents Jupiter's?) before assuming the `#535`/`#557` recipe transfers.
  **PLAN DONE (2026-07-15)** — `docs/superpowers/plans/2026-07-15-603-sun-neptune-transient-capture-
  plan.md`. **Headline recommendation: DECLINE a full sweep.** Two independent, verified findings:
  (1) Neptune's real orbital period is 164.8947 yr (computed in-repo, cross-checked via Kepler's
  third law) — the SAME dimensionless logic `#557` used gives a floor of ~165yr and a window of
  ~1649-2473 YEARS, ~13.9x further from the catalogue's 10-15yr figure than Jupiter's already-
  flagged departure; no better-motivated Neptune-specific derivation was found (checked, not
  assumed — the one real analog, 2007 RW10's ~7,500-12,500yr libration, argues for LONGER, not
  shorter). (2) A live web search (not recalled) found **no documented real object ever entering
  Neptune's actual Hill sphere** in a capture episode — the closest analogs (2007 RW10, four
  "temporary Neptune co-orbitals") are quasi-satellite/1:1-resonance objects that by definition stay
  OUTSIDE the Hill sphere, and Horner & Evans (2006) independently document Neptune has "great
  difficulty" capturing Centaurs even into that weaker coupling. The plan also directly checked and
  REFUTED this task's own "smaller e might survive" motivating hint against the project's own prior
  evidence: Earth's e=0.0167 is itself smaller than Jupiter's e=0.0484 yet still fully collapsed
  (`#535`), so there's no established pattern to extrapolate Neptune optimism from. A third,
  previously-uncontemplated complication was also found: Neptune's neck-open Jacobi-constant band is
  ~6.7x narrower than Jupiter's (tracks μ^(2/3) scaling), so even a reduced/fallback run would need a
  much finer energy-axis grid, not just a rescaled time window. **User decided 2026-07-15: decline
  outright** (not the weaker idealized-only fallback either) — no build, no run, no registry entry.
  Thread closed; the plan document stands as the record of why.

- **#604 ✓ DONE (2026-07-15) — INCOMPATIBLE, do not build the chain.** V∞/Tisserand compatibility
  check for the `#500`→`#318` Keplerian-map-chaining idea. `#500`'s own genome regime (RS07/GR09
  worked examples, C_J=2.995-3.00618) gives moon-relative v∞ ≈ 0.6-0.8 km/s via the code's own
  Tisserand relation (`v∞² = 3 - C_J`) — an intentionally near-zero, already-captured regime.
  `#318` itself never computed an interplanetary Jupiter-ARRIVAL V∞ at all (its EM Phase-1 probe is
  Earth-Mars, wrong system; its Phase-2/`#501` follow-up pivoted to a purely intra-Jovian CGCEC
  moon-tour search that closed 0/3072, `data/empty_regions.jsonl` region
  `jovian-cgcec-sobol-smoke-318-2026-06-30` — no Jupiter-arrival V∞ figure to compare against). Best
  available real proxies (sourced, `data/catalogue.yaml`, #390): Voyager 1/2 and Pioneer 10/11
  arrived at Jupiter at V∞ 7.6-10.8 km/s — roughly an order of magnitude above the map's valid
  regime; even published, ALREADY-CAPTURED intra-Jovian ballistic moon-cyclers (Russell-Strange
  2009) run 2.4-4.1 km/s, still 3-5x hotter. **A real interplanetary cycler would need a large,
  unmodeled capture/energy-reduction phase before the map's dynamics apply — the multi-day chaining
  build should NOT proceed as scoped.** No catalogue rows changed, no code built.
- **#605 ✓ DONE (2026-07-15) — planning/research deliverable, ranked shortlist produced.** A
  genuinely creative discovery-strategy pass, distinct from the tactical "which unswept body/
  method-pairing is next" passes (`#600`/`#601`/`#557`/`#599`/`#603`, all now closed, mostly clean
  negatives). Diagnosis: the program's recurring failure mode is the **family-selection/basin
  wall** (shooters/continuation only converge to families whose basin contains the seed — `#388`,
  S1L1, `#538`, `#520`'s scoping failure) — denser grids at new targets won't fix a structural basin
  problem, hence looking at paradigms instead of targets. **Ranked shortlist**: (1) variational/
  least-action seedless periodic-orbit discovery (attacks the basin wall directly — NOTE: the
  report's framing of Shijun Liao's SJTU mass-produced-orbit work as "variational/least-action" is
  WRONG, spot-checked directly against the source (`numericaltank.sjtu.edu.cn/three-body`) — Liao's
  actual method is high-precision "clean numerical simulation" + ANN-assisted seeding, a DIFFERENT
  technique than true action-minimization (the figure-eight's own discovery method); the underlying
  idea (a seedless/global method) may still be worth exploring via a genuinely variational approach,
  just not via Liao's specific technique as cited); (2) triple/quadruple small-body multi-moon
  systems (Sylvia/Elektra/Eugenia/Kleopatra/Lempo-Paha-Hiisi) — zero published record, mostly
  existing-tool reuse + a perturbation term for irregular primaries; (3) generative ML seed model
  (VAE) trained on this project's own corrector runlogs + `empty_regions.jsonl` — **spot-checked and
  CONFIRMED real**: Litteri, Gil, Vasile, Rodriguez-Fernandez & Camacho, "Generation of periodic
  orbits in the restricted three-body problem with a variational autoencoder," *Celestial Mechanics
  and Dynamical Astronomy* 138:25 (June 2026) — externally de-risks the parked `#542`/open `#317`
  ideas; (4) hierarchical "cycler-of-cyclers" (phase-match a heliocentric cycler's planetary
  encounters to a moon-system cycler at the target), pilotable cheaply via Mars Phobos-Deimos (a
  genuine roster gap — never swept, Russell-Strange 2009 methodology applies directly); (5)
  certified non-existence via interval arithmetic (upgrades `empty_regions.jsonl` from
  method-conditional to theorem-grade — speculative but uniquely suited to a program whose main
  output is clean negatives). A JAX/GPU differentiable-propagation port was recommended as
  accelerating substrate for 1/2/4, not a standalone idea. Full report not committed as its own doc
  (delivered as a planning response); this bullet is the durable record. No code built, no catalogue
  changes — purely a research/strategy deliverable for the user to pick from.

- **#606 ✓ BUILT + POSITIVE-CONTROLLED + PILOT SUCCESS (2026-07-16)** — `#605` shortlist item 1:
  seedless spectral (harmonic-balance) periodic-orbit discovery,
  `src/cyclerfinder/search/variational_periodic_orbit.py`. Represents a candidate closed loop as a
  truncated real Fourier series over one candidate period and drives the L2 residual of the CR3BP
  EOM to zero at collocation points via `scipy.optimize.least_squares` — the mathematically-
  equivalent, more tractable alternative to literal action-gradient-descent (same Euler-Lagrange
  stationary point; no new autodiff dependency). Cited precedent verified: Moore 1993 (PRL 70, 3675)
  + Chenciner & Montgomery 2000 (Annals of Math. 152, 881-901) — NOT Shijun Liao's SJTU technique, a
  correction to `#605`'s own report.
  **Positive control PASSED**: reproduces the Earth-Moon L1 planar Lyapunov orbit from a genuinely
  cold start (location/period/Fourier coefficients all far from the known answer), matching
  state/period/Jacobi to 1e-9 to 1e-13, cross-checked with a second (Radau) integrator; a broader
  cold-multistart sweep confirms this isn't a lucky single point.
  **Pilot: real wall-crossing, independently re-confirmed and directly re-run** (not just cited) —
  targeted `#556`'s L1 quasi-halo wall: 15/16 tested amplitudes fail the existing
  `richardson_halo_seed`+shooting-corrector combination (re-confirmed live, not assumed). A
  warm-started continuation using ONLY the new method's own solves crossed the ENTIRE wall region in
  20 steps, residuals falling monotonically to ~1e-11/1e-13, landing within 0.03% of the
  independently-confirmed bifurcation Jacobi constant (3.1745) — this project's existing seeded
  correctors cannot reach this region at all. **Honest caveat, reported not hidden**: a fully COLD
  one-shot attempt directly inside the deep wall (no warm start) only partially converged and the
  new method has its OWN family-selection bias (an unconstrained in-plane amplitude slides onto the
  wrong, already-known vertical-Lyapunov family instead of the halo) — a genuine limitation, not
  absent, just different from the shooting correctors' bias. Reviewed and independently re-verified
  by the coordinating session (re-ran `scripts/run_606_variational_pilot.py` live, matched the
  reported numbers exactly) before this entry was written. 7/7 new tests pass
  (`tests/search/test_variational_periodic_orbit.py`), ruff clean, `tests/scripts` (preflight AST
  ratchet) clean with a documented, reviewed exemption for the non-sweep demo script. No catalogue
  writeback (this is a capability build + method demonstration, not a discovery result to adjudicate
  — any future novel-orbit claim built on this tool would still need its own literature-novelty
  check).
- **#611 ✓ BUILT + POSITIVE-CONTROLLED (2026-07-16)** — `#605` shortlist item 4 / `#606` follow-up:
  point the new seedless spectral corrector at `#538`/`#544`'s QBCP EM-L1/L2 wall — the
  "violently-unstable" region (frozen-time linearization rate ~2-3 over `T_s~6.79`, implying
  ~1e6-1e8 per-period amplification) the existing multiple-shooting/GMOS corrector structurally
  cannot converge without an expensive multi-stage CR3BP-L1 → BCR4BP mu_sun-ramp → QBCP bootstrap.
  Picked up from a PRIOR agent that stalled/was killed after ~3.5h of silence; its
  `src/cyclerfinder/search/variational_periodic_orbit_qbcp.py` was substantially complete and its
  module docstring's technical claims were independently verified true: the QBCP genuinely does NOT
  let `#606`'s CR3BP approach generalize by a simple EOM swap — it is non-autonomous with a FIXED,
  known period (`period_multiple * T_s`, no continuous time-shift symmetry to gauge away, so no
  phase-gauge/free-period unknown carries over) and a first-order 6-canonical-state (PM) system
  (x,y,z,px,py,pz) rather than CR3BP's second-order 3-position system, so all six components are
  independently Fourier-expanded and matched against `qbcp_eom`'s first-derivative RHS.
  **Positive control PASSED, independently reproduced live this session (not just re-quoted from the
  docstring)**: cold-started with NO continuation bootstrap (`discover_qbcp_periodic_orbit(system,
  rng=np.random.default_rng(0))`, just the module's own rough default center guess ~0.041 from POL1
  and the known period `T_s`) converges in 8.9s to `residual_rms=4.86e-13`,
  `closure_residual=3.3e-6` (confirmed with a second, independent Radau integrator: 2.45e-6) — and
  matches this project's OWN independently-built 12-segment multiple-shooting corrector
  (`scripts/analyze_593_qbcp_l1_substitute_reconciliation.py`, resnorm 1.13e-14) to
  **machine-precision agreement (full 6-state diff ~5e-15)**, landing at `dist-to-POL1 = 1.805e-2`
  (xy-only) / `2.310e-2` (full-state) from the published golden — the SAME known
  Gimeno-2018-vs-Rosales-2023 model-instance gap `#544` already attributed to a Fourier-refit
  difference, not a corrector defect. This is a genuine wall-crossing: the region that forces every
  shooting corrector in this codebase through an expensive multi-stage bootstrap (or fails outright)
  is reached directly from a plain cold start in single-digit seconds.
  **Two real bugs found and fixed while producing the test file** (both independently reproduced
  live, not assumed from the docstring): (1) the shipped defaults (`n_harmonics=8`, tol=1e-6,
  8 random restarts) never actually converged — 8/8 restarts landed on the identical
  `residual_rms=9.4e-5` floor after 236s, because 8 harmonics is too coarse a truncation for this
  violently-unstable region: it can satisfy the collocation residual while `closure_residual` (an
  independent real-propagation check) is still O(1) (`0.63`) — i.e. NOT a periodic orbit despite
  "converging". Default raised to `n_harmonics=32`, which drives `closure_residual` down to
  `~1e-6`. (2) `scipy.optimize.least_squares(method="lm")` does NOT respect `max_nfev` as a real
  wall-clock bound on this problem — measured directly, requesting `max_nfev=100` cost 38,710 actual
  residual calls (~390x over), which is almost certainly what caused the prior agent's 3.5h silent
  stall. Switched to `method="trf"` (bounded to tens of seconds per attempt even for an unlucky
  random restart, vs. `"lm"`'s effectively-open-ended cost) and reduced default `max_nfev` from
  30000 to 1500 accordingly. Remains a stochastic method post-fix: most cold starts converge in ~9s,
  but an occasional unlucky seed can take several minutes (observed once, still bounded, not
  indefinite) before landing on the same answer.
  **Explicit scope boundary (verified, not overclaimed): this does NOT solve `#538`/`#544`'s actual
  named target**, the QBCP EM-L1/L2 invariant 2-TORUS (`genome.qbcp_torus.correct_qbcp_torus`) — a
  genuinely quasi-periodic 2-angle family with a free rotation number and stroboscopic-map
  invariance. This module crosses the wall for the PERIODIC ORBIT anchoring that torus family as its
  zero-amplitude center, not the torus corrector itself; generalizing to a genuine 2D quasi-periodic
  build is a materially larger, out-of-scope follow-up, not attempted here.
  8 new tests (`tests/search/test_variational_periodic_orbit_qbcp.py`, includes the
  low-harmonics-fails-closure regression + the machine-precision multi-shooting cross-check), ruff
  clean, sibling `tests/search/test_variational_periodic_orbit.py` (untouched) still 7/7 green — all
  independently re-run and reproduced by the coordinating session (numbers reproduce exactly). No
  `scripts/run_*.py` created (library module + test file only), so the preflight AST ratchet does
  not apply. No catalogue.yaml writeback — capability build + method demonstration only.
- **#612 ✓ BUILT + POSITIVE-CONTROLLED + WALL CROSSED (2026-07-16)** — user-approved `#611`
  follow-up: extend the seedless spectral method to `#556`'s parking-lot target, the EM L1
  quasi-halo QUASI-PERIODIC 2-TORUS that `genome.qp_tori.correct_qp_torus` cannot converge above
  amp~0.01. New module `src/cyclerfinder/search/variational_qp_torus.py` + tests
  `tests/search/test_variational_qp_torus.py`. **`qp_tori.py` was NOT modified** (no bug found in
  it — its wall is intrinsic, see diagnosis).
  **Diagnosis (hypothesis CONFIRMED with numbers, not assumed).** The GMOS/invariant-circle
  corrector resolves the LONGITUDINAL angle by stroboscopic FLOW INTEGRATION: `_gmos_residual`
  propagates every invariant-circle sample forward by `t_strob=2π/ω_long` and matches the rotated
  circle. The parent L1 halo at C=3.15 (rebuilt live: C=3.15043, T=2.76109, z0=-0.05453) has
  **monodromy spectral radius ~1540** (‖M‖₂~2789), so that one-period propagation AND every
  finite-difference column of its Jacobian is amplified ~1540x. Reproduced the wall exactly:
  invariance residual 5.5e-7 (amp 5e-4) → 2.5e-4 (amp 1e-2) → TIMEOUT at amp 0.015 (90s) → 2.2e-3
  unconverged at amp 0.02 (34s). This is a *shooting* fragility in an invariant-circle disguise —
  the same class `#606`/`#611` crossed — NOT a conceptual limit of the torus.
  **Built: 2D pseudospectral corrector, no integration anywhere in the search.** Fourier-expand all
  six state components in a real tensor-product basis over BOTH angles; drive the quasi-periodic
  invariance PDE `ω1·∂_θ1 u + ω2·∂_θ2 u − f(u) = 0` to zero on a 2D collocation grid via
  `scipy.optimize.least_squares(method="trf", x_scale="jac")` with an **exact analytic Jacobian**
  (the local operator `ω1 ∂_θ1 + ω2 ∂_θ2 − Df(u)` — the ~1540x monodromy amplification never enters
  the residual or its Jacobian, which is exactly why it beats GMOS; analytic Jacobian also removes
  qp_tori's documented ~1e-7 finite-difference noise floor, verified vs central FD to rel 6e-7).
  Both `(ω1, ω2)` are free unknowns. Four gauge/anchor rows: longitudinal phase (`#606`'s sin_x1=0
  analogue), transverse phase, transverse amplitude anchor, and a **rotation-number pin**
  `ω2=ρ·ω1`. The rho-pin is load-bearing and was found the hard way: without it the least-squares
  escapes to a degenerate spurious branch (`ω2→0`, a resonant "tube of periodic orbits" that hits
  machine-precision residual and so always beats the truncation-limited genuine torus). Pinning ρ
  is Jorba-standard AND physically correct here (`#555` showed the L1 ρ is energy-pinned, flat in
  amplitude). Bootstrap: project a small-amp GMOS torus onto the 2D basis (integration allowed in
  the GUESS only), then continue amplitude with the module's own integration-free solves.
  **Positive control (L2, regression floor) PASSED.** On the already-converging EM L2
  near-bifurcation quasi-halo at C=3.15 (T=3.411, z0=0.0171), the pseudospectral corrector
  reproduces the GMOS torus: rotation number 0.023271 vs GMOS 0.023272 (matches to 1e-6, in Owen &
  Baresi's ~0.0216 latitudinal regime), residual_rms 3.7e-6 (BETTER than the GMOS seed's own
  invariance residual 5.8e-6), independent closure 8.5e-7.
  **L1 wall CROSSED (the headline).** Cold-continued the transverse amplitude from the small-amp
  bootstrap to transverse-amp 0.06 = GMOS `initial_torus_amplitude`-equivalent ~0.0197 (bootstrap
  scale 3.04x), in 12 steps ALL converged, holding independent closure 3.6e-9 → 4.2e-8 the whole
  way — i.e. clean convergence at the very amplitude (~0.02) where GMOS times out / returns invres
  2.2e-3. Total continuation ~94s. The resulting torus is a genuinely large quasi-halo (z librates
  ∈[-0.069, 0.058] about the parent halo's z0=-0.0545; max transverse excursion 0.035). Rotation
  number held energy-pinned at -0.074024 throughout.
  **Honest scope boundary (verified, not overclaimed): this crosses the CONVERGENCE wall, it does
  NOT reach Owen & Baresi's L1 latitudinal frequency 0.2739.** Those are two different things
  `#556`'s framing conflated. `#555` established — and this task re-confirmed directly — that at
  C=3.15 the L1 quasi-halo rotation number is ENERGY-pinned near 0.074, flat in amplitude; 0.2739 is
  simply not the rotation number of any L1 torus at this energy, a physical family fact independent
  of which corrector is used. This module removes the corrector limitation `#556` named (converging
  the large-amplitude L1 torus that GMOS could not); reaching 0.2739 would require a different
  energy/family, not a better corrector, and is a separate question. No catalogue writeback
  (capability build + method demonstration, not a discovery result).
  **Tests:** 10 new tests (8 fast unit — basis-derivative, pack/unpack, analytic-Jacobian-vs-FD,
  gauge-row structure, torus evaluator; 2 integration — the L2 positive control and the L1
  wall-crossing), all pass; ruff clean (module + tests). Sibling `test_variational_periodic_orbit.py`
  and `test_variational_periodic_orbit_qbcp.py` stay green; `test_qp_tori.py` green EXCEPT its
  pre-existing `test_structural_qp_continuation`, which fails on this macOS-libm machine
  independently of this work (proven: it fails identically on a clean tree with the two new files
  removed — an over-strict strict-monotonicity assertion on ρ values varying at the ~5e-6 level,
  below the ~1e-7 FD noise floor; not in `#612`'s scope to fix). No `scripts/run_*.py` created, so
  the preflight AST ratchet does not apply. Left uncommitted for coordinating-session review.
  **Recommended model (follow-up):** Sonnet — the hard judgment (diagnosis + rho-pin design) is
  done; any productionization (default mode counts, pseudo-arclength continuation across ρ, wiring
  into the `#522` linking-number screen for a full O&B run at an energy where ρ actually reaches
  0.2739) is now spec-complete work behind these regression tests.
- **#613** (dispatched 2026-07-16, `#612` follow-up) — find whether ANY Jacobi constant `C` gives
  the EM L1 quasi-halo family a rotation number matching Owen & Baresi's target latitudinal
  frequency 0.2739 (`#555` proved it is NOT 0.2739 at C=3.15 specifically — energy-pinned near
  0.074 there, flat in amplitude), and if so, build that torus with `#612`'s integration-free
  pseudospectral corrector and wire it into the `#522` linking-number screen for the actual final
  Owen & Baresi L1<->L2 reproduction attempt (paired with a genuine L2-side torus at the SAME
  energy). Map rotation number vs. parent-halo `C` near the L1 bifurcation (C=3.1745, per `#555`)
  using small-amplitude GMOS-bootstrapped tori (cheap, since GMOS converges fine below amp~0.01) at
  several `C` values without forcing a rho pin, to locate (if it exists) the `C` where the natural
  rotation number crosses 0.2739; only then use `#612`'s corrector (pinned to that `C`'s natural
  rho) to build the large-amplitude torus needed for the actual scan. A well-characterized negative
  ("0.2739 is unreachable at the L1 quasi-halo across the full physically-relevant `C` range") is an
  entirely acceptable, valuable, FINAL answer to the question the `#548`->`#553`->`#555`->`#556`->
  `#612` chain has been chasing — do not force a positive result. **Recommended model: Sonnet**
  (per `#612`'s own note — this is spec-complete productization work, not open-ended numerical
  design).
- **#607 ✓ CLEAN NEGATIVE (2026-07-16)** — `#605` shortlist item 2: triple/quadruple small-body
  multi-moon systems — (87) Sylvia (Romulus+Remus), (130) Elektra (3 moons), (45) Eugenia, (216)
  Kleopatra, and the TNO triple Lempo-Paha-Hiisi. `#549`'s real-binary `(k1,k2)` genome does NOT
  generalize (fixed two-primary construction, no room for a 3rd gravitating body); reused `#563`'s
  2-moon and `#600`'s 3-moon-chain symmetric-closure construction verbatim instead (both already
  genericized). **Lempo-Paha-Hiisi excluded structurally**: Lempo:Hiisi mass ratio ~1.27:1, a genuine
  near-equal-mass binary, not primary+test-particle-moon — violates this method's core assumption.
  New satellite registry entries added to `src/cyclerfinder/core/satellites.py` (Sylvia/Elektra/
  Eugenia/Kleopatra + their moons), every GM/mass sourced to a specific paper (Vernazza/Carry 2021,
  Fuksa et al. 2023, Beauvalet & Marchis 2014, Marchis & Yang 2021, Fang/Margot/Rojo 2012, etc.),
  with explicit caveats where a moon's mass is assumed-density rather than dynamically measured, and
  Kepler-III self-consistency checks run against each sourced period (all matched to <3.5%).
  Independently spot-checked by the coordinating session: every GM = mass × G arithmetic conversion
  and the ElektraBeta Kepler-III self-check (predicted 5.12 d vs. sourced 5.287 d) reproduce exactly.
  **Result**: 0/97,664 candidates pass all gates across the 4 included systems (Sylvia, Eugenia,
  Kleopatra: 2-moon, 384/384/128 evaluated; Elektra: all 6 ordered 3-moon-chain permutations, 96,768
  evaluated) — confirmed via the run's own `_meta` summary record. Diagnosed directly (not assumed):
  spot-checked every sub-residual-gate candidate for Sylvia (70/70) and Kleopatra (30/30) — 100% fail
  the physical `#324` bend gate (0.1-8° achievable vs. the required ≥5°), 0% fail an independent
  DOP853 cross-check (construction is numerically sound). Same mass-limited failure mode already on
  record for Jupiter's Amalthea and Neptune's Triton — these moons' GMs (1e-7 to 6e-5 km³/s²) are too
  small to bend a spacecraft at the system's natural ~10-50 m/s velocity scale. Registered as
  `empty_regions.jsonl` region `smallbody-multimoon-symmetric-closure-mass-limited-607-2026-07-16`,
  validated via `validate_empty_region()`. New tests: `tests/core/test_satellites_registry.py`
  (registry coverage/consistency extended to the 4 new systems) + `tests/scripts/
  test_enumerate_607_smallbody_multimoon_symmetric_closures.py` (7 tests) — all pass, ruff clean. No
  catalogue.yaml edit.
- **#608 ✓ BOUNDED POC BUILT (2026-07-16)** — `#605` shortlist item 3: generative ML seed model
  trained on this project's own corrector-outcome logs, testing the `#542`/`#317` idea externally
  de-risked by Litteri, Gil, Vasile, Rodriguez-Fernandez & Camacho, "Generation of periodic orbits in
  the restricted three-body problem with a variational autoencoder," *Celestial Mechanics and
  Dynamical Astronomy* 138:25 (June 2026) (paywalled; worked from the abstract + the group's earlier
  open conference version, arXiv:2408.03691, which gives the real architecture: CNN-VAE, 2D latent,
  44,112 ICs across 40 Earth-Moon families each as a 100-node time-series, refined by multiple-
  shooting, 46% of 100 latent samples converged to a genuinely new orbit).
  **Design decision, justified not assumed**: this project has NO torch/jax dependency (checked
  `pyproject.toml`); a literal CNN-VAE was judged not worth a new heavy dependency for a bounded POC.
  Built instead `src/cyclerfinder/ml/orbit_generative.py`: standardize -> linear PCA -> k-means(8)
  partition -> per-cluster empirical Gaussian in a 5D latent space — the numpy/scipy-only
  linear-Gaussian analog of their CNN-VAE, operating directly on the (state0, period, Jacobi) genome
  rather than a 100-node time series (that representation is already the CR3BP periodic orbit's
  minimal sufficient statistic; Litteri's time series was a consequence of their CNN architecture,
  not a physical necessity).
  **Training data, assembled from this project's OWN history, not synthetic**: scanned all
  `out/outcome_log/*.jsonl` (#210's passive corrector-outcome capture; 540,312 raw lines from a dozen
  past Earth-Moon CR3BP campaigns). Only ~38% were solver-`converged`, and of THOSE only ~26% were
  physically plausible (Jacobi ranged as wild as [-164, 105] before filtering — the raw log is
  contaminated with numerically-converged but degenerate/collision-adjacent junk, the same failure
  mode this task's own generate-then-refine step below rediscovered). After a documented
  Jacobi/period/out-of-plane-amplitude filter + dedup: **54,165 unique, physically-sane converged
  orbits** (43,332 train / 10,833 held-out) — comparable in raw count to Litteri's 44,112, though
  almost certainly far less FAMILY-diverse (the log has no family label, only primary/secondary; a
  real family tag would be a natural follow-up).
  **Results (`scripts/run_608_generative_seed_poc.py`, `data/found/608_generative_seed_poc/`)**:
  (a) reconstruction sanity check — 100 held-out real orbits round-tripped through the model's own
  PCA encoder/decoder and refined by the EXISTING `cr3bp_periodic.correct_periodic`: 60% raw
  convergence, 51% converged-AND-physically-sane. (b) generate-then-refine — the RAW convergence
  rate looked like a null result (55% generated vs. 51% uniform-random-in-bounding-box baseline,
  ~1.08x) — **but that raw metric is misleading for the identical reason the training corpus itself
  needed filtering**: a solver can converge to a numerically-valid, physically-degenerate solution.
  Once both are filtered by the SAME physical-sanity bounds used to build the corpus, the real
  comparison is **49% (generated) vs. 4% (uniform baseline) — a genuine ~12x improvement**, i.e. the
  learned density steers the existing corrector into a real family's basin far more often than blind
  seeding. (c) novelty — the converged generated candidates sit far from any single training example
  (median nearest-neighbor distance ~40.5 standardized units, vs. ~0.27 for two independent real
  orbits) — genuinely new, not near-duplicates, though also evidence the linear PCA model doesn't
  tightly track the true (thin, curved, nonlinear) family manifolds, landing in physically-plausible
  "gaps" between them rather than tightly on a family curve.
  **Honest verdict**: partially and meaningfully viable NOW with this project's own data (54k
  physically-sane examples is enough for a bounded statistical model to give a real ~12x seeding lift
  over blind search), but the raw corpus's ~74% junk rate and lack of family labeling mean a future,
  larger build should (1) tag provenance/family per record so the model can be conditioned/stratified,
  and (2) consider a nonlinear model (kernel PCA / shallow autoencoder / eventually a real VAE) given
  the linear model's demonstrated difficulty tracking curved manifolds — do not jump straight to a
  production VAE pipeline on this evidence alone. Does NOT supersede `#317` (a sweep-infeasibility
  PRE-FILTER, a different question from seed generation) though both could share the same corpus;
  `#542`'s learned-seed warm-start idea is now empirically de-risked rather than purely speculative.
  8 new tests (`tests/ml/test_orbit_generative.py`), ruff clean, `tests/scripts` AST ratchet clean
  with a documented exemption (same category as `#606`'s). No catalogue writeback, no
  literature-novelty check — capability POC only, per this task's explicit scope. Independently
  re-run by the coordinating session (fixed `seed=608`): numbers reproduce bit-for-bit, including
  the 12.25x generated-vs-baseline physically-sane convergence ratio.
- **#609** (P1, new capability) — `#605` shortlist item 4: hierarchical "cycler-of-cyclers" —
  phase-match a heliocentric cycler's planetary encounters to a moon-system cycler at the target
  (commensurability between the synodic super-period and the moon-cycler period). All lower-level
  ingredients already exist (`cross_system_cycle`, phase-match/endgame machinery); the new piece is
  the composition constraint itself. Cheapest pilot: Mars Phobos-Deimos (never swept — a genuine
  roster gap, Russell-Strange 2009 methodology applies directly) as the target moon-cycler for an
  existing Earth-Mars heliocentric cycler.
  **STEP 1 CLOSED (2026-07-16): CLEAN NEGATIVE, pilot cannot proceed to the phase-matching step.**
  Ran the already-genericized `#563`/`#575`/`#599` direct symmetric-closure enumeration
  (`scripts/enumerate_563_symmetric_closures.py --primary Mars --moons Phobos,Deimos`) at Mars for
  the first time — Mars GM (4.282837521e4 km^3/s^2) and Phobos/Deimos GM/radius/sma
  (`core/satellites.py`, JPL SSD MAR097-sourced, already present pre-task) were verified, not
  assumed. Result: 0/512 candidates (both directions, Phobos-anchor and Deimos-anchor) pass ALL
  gates. 52/512 clear the residual-closure gate alone (both moons orbit so close together, T_syn =
  0.427 d, that the kinematic Lambert closure is essentially exact, residual ~1e-13-1e-15 km/s) but
  EVERY one fails the two-sided `#324` physical max-bend gate: the flyby moon's own achievable bend
  tops out at 0.0159 deg across all 52 sub-gate survivors, ~300x under the
  `DEFAULT_MIN_USEFUL_BEND_DEG=5.0` deg floor. Same failure mode as `#571`
  (Saturn-Titan-small-moons) and `#599` (Neptune-Triton-Proteus), but here BOTH bodies in the pair
  are individually undersized (Phobos GM 7.087e-4, Deimos GM 9.62e-5 km^3/s^2) — `core/satellites.py`
  itself already flagged both as making "poor gravity-assist / cycler bodies" (pre-existing sourcing
  note), now confirmed quantitatively. Checked `data/catalogue.yaml` (no Phobos/Deimos rows exist)
  and `literature_check.py`'s own Wallace Mars-Phobos CR3BP-rendezvous `CorpusAnchor` (which already
  notes no published Sun-Mars-Phobos BCR4BP or repeated-flyby cycler study exists) — no prior
  published or catalogued Phobos-Deimos moon-cycler either. **Consequence: with no gate-passing
  moon-cycler to compose with, step 2 (phase-matching a catalogued Earth-Mars heliocentric cycler's
  Mars-arrival epoch against the moon-cycler's own period) cannot be attempted — there is nothing to
  phase-match against.** This is the `#609` deliverable itself, a legitimate bounded negative, not a
  partial result. `data/empty_regions.jsonl` stamped (`mars-phobos-deimos-symmetric-closure-609-
  2026-07-16`); raw sweep data at `data/enumerate_609_mars_phobos_deimos_symmetric_closures.jsonl`;
  one-shot writer `scripts/_apply_609_mars_phobos_deimos_empty_region.py`. No catalogue.yaml edit,
  no commit (left for the coordinating session to review).
- **#610** (P2, speculative) — `#605` shortlist item 5: certified non-existence via interval
  arithmetic over a compact phase-space region, upgrading an `empty_regions.jsonl` entry from
  "conditional on method" (the registry's own stated epistemic status) to a theorem-grade
  non-existence certificate.
  **✓ Proof-of-concept CERTIFIED (2026-07-16)**, `scripts/certify_610_proteus_bend_interval.py` +
  `tests/scripts/test_certify_610_proteus_bend_interval.py` (13/13 pass; independently re-run and
  the exact numbers reproduced by the coordinating session). Targeted the
  `neptune-triton-proteus-symmetric-closure-599-2026-07-15` entry (`#599`): its own recorded verdict
  isolates a single, non-transcendental sub-condition — Proteus's GM too small to deliver ≥5° of
  ballistic bend — as the entire reason for the negative, tractable via the closed-form
  Bate-Mueller-White patched-conic formula (`cyclerfinder.core.flyby.max_bend`). Built a rigorous
  interval arcsin from `mpmath.iv`'s `atan2`/`sqrt` primitives (`iv` exposes no `asin`/`acos`
  directly) and certified `sup(bend_deg) < 5°` over two boxes: Box A (data-grounded, r_p pinned at
  the 308 km safety floor, V∞ = the exact [1.824, 13.117] km/s range spanning all 104 real
  residual-gate survivors) → sup=0.288°; Box B (widened/conservative, r_p ∈ [308, 3080] km, V∞ ∈
  [0.45, 20] km/s) → sup=4.559°. This upgrades the entry's "Proteus GM too small" interpretation
  from a 104-point grid finding to a continuum-strength claim: no point in either uncountable box
  can pass the gate. Explicitly does NOT certify the residual/Lambert-closure half of the #599
  search (bounding a multi-revolution, branch-selecting universal-variable Lambert solve rigorously
  is a genuine unresolved technical obstacle, reported honestly rather than forced), nor any other
  empty_regions.jsonl entry — a scoped, single-entry, single-failure-mode proof of concept, not a
  general framework. Added an optional `interval = ["mpmath>=1.3"]` extra to `pyproject.toml`
  (tests skip cleanly without it). No catalogue.yaml/empty_regions.jsonl writes — registering the
  upgraded claim in the actual registry is a follow-up decision, not yet made.

- **#592** (P2, correctness fix) — recovered a real QBCP equations-of-motion bug from an abandoned
  `git stash` (found during a repo-hygiene cleanup pass across both repos, undated stash on `7f83277`)
  and applied it after independent verification. `qbcp_eom()`/`qbcp_potential_second_derivatives()`
  in `src/cyclerfinder/core/qbcp.py` multiplied the WHOLE Newtonian potential gradient (Earth+Moon+Sun)
  by `alpha_6`; per Rosales-Jorba (2023) Eq. 3 (an already-trusted corpus reference, used for the
  POL1/POL2 golden tests), `alpha_6` scales ONLY the Sun-distance term. The existing golden test sits
  at `alpha_6=1` exactly (a documented simplification), where old and new code coincide numerically —
  which is why this was never caught: `alpha_6` is in general a genuine time-varying coefficient
  (`evaluate_alphas()`), so any real QBCP propagation picked up a systematic O(alpha_6-1) error the
  golden test can't see. Verified before applying: citation grounded (not inherited-unverified), the
  3 files referencing `qbcp` all still pass identically, and a broader `tests/genome`+`tests/core` run
  surfaced 3 failures that reproduce byte-identically WITHOUT the fix too (pre-existing, unrelated —
  filed alongside #584's known local-Mac-artifact category). Full writeup:
  `docs/notes/2026-07-14-qbcp-alpha6-scaling-fix.md`. **#592 STATUS: CLOSED (fix applied).**

- **#593** (P3, follow-up, not yet scoped) — determine whether #592's QBCP alpha_6 fix changes any
  past QBCP-based search/continuation conclusion. Per [[feedback_bugfix_invalidates_past_searches]]:
  any result computed under the buggy code (the #533 QBCP model build, #537/#538 cross-system
  connection/correction work, #544's EM-L2 mu_sun-continuation blocker) used a subtly-wrong Sun-term
  scaling. Genuinely uncertain whether the O(alpha_6-1) magnitude is large enough to have mattered for
  any specific past conclusion — this needs its own scoping pass (which past results are
  alpha_6-sensitive enough to re-check), not a blind full re-run.
  **DONE (2026-07-14).** `docs/notes/2026-07-14-593-qbcp-alpha6-impact-scoping.md`. Bounded
  `|alpha_6-1|` at up to 0.84% over a synodic period. **Empirically re-ran the SE-L2 torus
  correction (#533/#537's clean positive control) under both the fixed and reconstructed
  byte-exact buggy code from the same seed: both converge cleanly by their own residual criterion
  (~1e-5), but land on a DIFFERENT mean-state solution (max delta 0.194, vs. only a ~0.1% shift in
  the rotation number rho).** This is a genuine, non-trivial impact, not a negligible perturbation —
  the SE-L2 numbers quoted throughout #533/#537/#544 as the clean baseline do not reproduce under
  the corrected code. #544's qualitative EM-L2-vs-SE-L2 contrast likely still holds (EM-L2's
  residual of 0.34-3.4 is 4-5 orders of magnitude larger than this effect), but its specific quoted
  SE-L2 numbers are stale. Also flags a verification-methodology gap: #533's own build-time
  finite-difference check could never have caught this bug (FD verifies a formula against its own
  derivative, not against the source paper's physics). **Recommendation: do NOT blindly re-run
  #538/#544's multi-hour effort now — whoever next picks up #538/#539/#540 should re-derive a fresh
  baseline under current `main` rather than reuse the old quoted numbers.** No code changes; a
  scoping conclusion only.
  **CRITICAL FOLLOW-UP, found AFTER first closing this out (2026-07-14):** closer reading of #544's
  own entry revealed **#544 already tried this exact fix on 2026-07-10 and explicitly reverted it**
  — reporting it moved a QBCP L1-substitute periodic orbit ~30% closer to the published POL1 golden
  (2.12e-2→1.47e-2) but allegedly **regressed** the SE-L2 positive control 16x (3e-5→4.2e-4), the
  OPPOSITE of what this task's own SE-L2 re-test found. This was missed before applying #592 — a real
  gap in the original verification (should have searched OUTSTANDING.md for prior `alpha_6`
  discussion first). Reconciled with two further checks: (1) ruled out a mismatched-EOM/STM
  explanation for the SE-L2 discrepancy (`qbcp_potential_second_derivatives` isn't even used by this
  correction path); the reported SE-L2 regression does NOT reproduce with the actual, faithfully-
  reconstructed code. (2) Built a genuinely independent multi-shooting L1-substitute construction
  from scratch (12-segment analytic-STM corrector, CR3BP L1 -> BCR4BP mu_sun-continuation -> QBCP
  handoff, cross-checked at every stage against known references) — **confirms #544's OWN finding,
  independently: distance-to-POL1 goes 4.14e-2 (buggy) -> 1.81e-2 (fixed), same direction and
  magnitude as their 2.12e-2->1.47e-2.** Net: one of #544's two original objections (SE-L2
  regression) is refuted, the other (POL1 improvement) is independently confirmed as a genuine
  improvement. The one remaining open question (whether Gimeno-2018's own alpha_6 table was fitted
  assuming this scaling convention) isn't numerically testable without the source paper, but the
  empirical balance now favors #592 as a net-positive fix. **#592 stays applied, no revert.** Full
  reconciliation: `docs/notes/2026-07-14-593-qbcp-alpha6-impact-scoping.md` (updated in place).
  **#593 STATUS: CLOSED (fully reconciled).**

- **#594** (P2, repo-hygiene / documentation-accuracy) — full-file comprehensive audit of every
  numbered task in this document (~180 entries, 6 parallel independent verification passes each
  reading its own line-range end-to-end, cross-checking header claims against actual entry
  bodies), triggered by the #545 confusion found while answering the user's own question about
  it. **Origin:** the user asked for every single audit finding to be checked and the file
  corrected, not just spot-checked.
  **DONE (2026-07-15).** Verified every task's TRUE status (not just its header) across the whole
  file. Key corrections actually applied to this file (not just found and left as a report):
  1. **#545** — was misclassified (by an earlier, less careful audit pass, and by me initially)
     as "ready to dispatch, just needs a new positive control." Direct inspection of its actual
     committed script (`scripts/run_545_jovian_band_phase_screen.py`) showed it genuinely depends
     on the retired qp_torus/heteroclinic pipeline (not an independent method as its prose
     suggested), AND its whole scientific question was already comprehensively answered by
     #576/#577 using a proven, superior method. Corrected to DEAD + SUPERSEDED, not revivable.
  2. **#521** — a stale "Phase 2 remains open" line contradicting an adjacent, more detailed,
     already-DONE note for the same Phase 2 (the AST ratchet is real and live, independently
     re-verified this session) — corrected in place.
  3. **#536** — reordered/clarified an internally-scrambled entry (STAMP block appeared before
     the pre-resolution proposal text it was retracting, reading as if un-resolved).
  4. **#537** — added the same explicit `⚠ STAMP` marker #534/#536 already carry (same
     functional DO-NOT-CERTIFY status, just missing the formal marker before this pass).
  5. **#405** — added a header-level retraction flag (its own last line reverses the "clean
     negative" the header claims, previously only visible after reading the whole paragraph).
  6. **Stale duplicate queue (#307/#310/#320/#321/#322)** — an entire un-weeded 2026-06-16-era
     pending list, each item already resolved elsewhere in the file under a DIFFERENT (correct,
     better) result, cross-referenced in place rather than left contradicting the real answer.
  7. **#318** — consolidated 3 non-adjacent, easy-to-misorder status statements into one clear
     note (net: BOTH #318 through Phase 2b AND its #501 follow-up are resolved/closed — #501 was
     actually dispatched and run to completion, corrected further below after item 12).
  8. **#248** — added the missing cross-reference to its own resolution, which sat ~150 lines
     away under a differently-named section.
  9. **#344 Phase 2 Stage A** — corrected further below after item 12 (my first pass here claimed
     a genuine still-open gap; a Fable adversarial check found this claim was itself wrong).
  10. **#506** — bolded a "DO NOT REOPEN" marker after confirming it had already been
      re-litigated once by mistake (#550).
  11. **#518 vs Tier-5 plan** — clarified #518 is superseded by #522, not independent forward work.
  Full per-chunk verification tables (all ~180 tasks, every classification with line-level
  evidence) exist in this session's agent dispatch transcripts; the corrections above are the
  subset that changed what a reader would conclude, applied directly rather than left as an
  external report.

  **12. Fable adversarial verification pass (2026-07-15), same day — checked every correction
  above against primary text, re-examined the 4 items judged "already adequately clear," and did
  a fresh independent scan.** Verdict: 9/12 corrections verified fully accurate as written
  (#545, #521, #536, #537, #405, the stale-duplicate cluster, #506, #518-vs-Tier-5, and items 1-6/
  8/10/11 of this list). 3 needed further correction, all now applied:
  - **#318/#501 (item 7 above)** — my first version claimed #501 was "OPEN-GATED pending an
    explicit decision." Fable found this cites a stale "[held]" marker on #501's own entry while
    ignoring the Tier-1 results block (near the top of this file), which shows #501 actually ran
    to completion: "0 closed; positive control (Liang Member D) PASSED → trustworthy. Clean
    empty-region map registered" — and #576/#577 both treat that result as settled history.
    Corrected: both #318 and #501 are resolved/closed; #501's own stale "[held]" marker fixed too.
  - **#344 (item 9 above)** — my first version claimed a GENUINE still-open gap (that #349's
    closure didn't address the literature-freshness blocker). **This claim was wrong.** Fable
    found `src/cyclerfinder/search/literature_check.py` explicitly documents #349 fixing exactly
    this blocker (a `topology_label` discrimination separating the Cassini-Huygens tour from a
    genuine repeated-moon cycler candidate, despite the shared {Titan, Rhea} body subset) — plus
    the candidate itself was independently adjudicated V0-known (Russell-Strange 2009,
    2026-06-30, reconfirmed after #488/#489). Corrected: resolved not-novel, not open. **Lesson
    for future audits: verify a "the fix doesn't address the blocker" claim against the actual
    CODE the fix touched, not just the one-line closure summary in this tracker — the summary
    undersold what #349 actually built.**
  - **#144 and #167** — judged "already adequately clear," WRONG on both. #144's "PROVISIONAL
    pending Guzman et al. 2002" gate was stale: Guzman 2002 was already acquired and mined
    (`docs/notes/2026-06-07-guzman-2002-primer-survey-mining.md`), which already delivers the
    exact re-label verdict — applied directly to #144's own entry. #167's "V3 writeback...HELD to
    user decision" was also stale: `data/catalogue.yaml`'s `russell-ch4-4.991gG2` row already
    carries `validation_level: V3` (mechanically backfilled from recorded gauntlet evidence) —
    the writeback already happened. Both fixed. (#164/#165 ordering and #110's "OPEN RESEARCH"
    header were re-checked and confirmed genuinely fine, no change needed.)
  **Net: the file is now verified accurate on every specific claim checked** (12 original
  corrections + 3 further corrections found by adversarial review = 15 total edits this pass),
  not just "audited and mostly trusted."

- **#554** (P2, cheap, ~1 day per the #552 scoping estimate) — Neptune/Amalthea empty-region
  retrograde-correction stamp. Formalize the #552 scoping pass's back-of-envelope flyby-bend
  analysis (Amalthea GM=0.165: max bend 0.80 deg at V_inf=0.5 km/s, collapsing at higher V_inf,
  independent of geometry; Triton correctly modeled retrograde forces V_inf up to ~6-9+ km/s
  where bend collapses to ~1-2 deg, versus the existing coplanar-prograde model's generous
  23-59 deg) into a proper, code-verified note, then append `reverification` entries to the 3
  Neptune rows + the Amalthea row in `data/empty_regions.jsonl` recording that these
  "conditionally empty pending 3D/inclined capability" verdicts are actually robust — indeed
  strengthened — under correct geometry, narrowing Amalthea's live re-sweep condition to
  low-thrust relegs only (its other named condition) and removing the false "3D relegs might
  reopen this" implication from all 4 rows. Pure registry-accuracy work, no new search, no new
  genome capability, no catalogue.yaml edit.
  **Recommended model:** Sonnet (mechanical formalization of an already-scoped analysis, behind
  the existing empty-region-registry schema/tests).
  **✓ Resolved (2026-07-10).** Recomputed both bodies' max-bend numbers from sourced GM/altitude via
  `core/flyby.py::max_bend` (Amalthea 0.8011 deg at v_inf=0.5, collapsing to 0.0001-0.008 deg at the
  v_inf=5.1-37.2 km/s actually realized in `data/scan_433_jupiter_galilean.jsonl`; confirmed in code
  that `search/discovery_campaign.py::_moon_state` + `core/satellites.py` carry no retrograde flag at
  all — Triton's real rotation exists only in comments — then derived the retrograde-corrected v_inf
  (~8.5-8.7 km/s via the equal-magnitude chord identity) and bend (~1.46-1.52 deg) from Triton's
  sourced GM/r_p). Both confirm #552's estimates. Wrote
  `docs/notes/2026-07-10-554-retrograde-correction-stamp.md` and appended `reverification` entries to
  the Amalthea row (re-sweep condition narrowed to low-thrust relegs only) and the 3 Neptune/Triton
  rows (verdicts marked strengthened, no re-sweep condition added) in `data/empty_regions.jsonl`.
  `uv run pytest tests/data -q`, `ruff check .`, `ruff format --check .`, `mypy src tests` all clean.

- **#539** (P1, do first once #538 lands) — Generalize #538's well-posed multi-segment
  torus-connection corrector out of the one-off `scripts/run_538_qbcp_cycler.py` into a
  reusable `genome/torus_cycle_corrector.py` module (mirrors how #314's
  `heteroclinic_cycle.py` became the shared base #405/#530/#531 all reused after their
  first one-off proved out) — **and** re-run the Jupiter-Europa/Ganymede screen properly:
  a Jacobi-constant band sweep (not one point) crossed with a synodic-phase-offset sweep
  between the two moons (the same free coordinate #522's coherent-model reformulation
  introduced for SE<->EM, generalized to a Jovian-moon pair), using #524's already-built
  deflated Newton to enumerate basins instead of hoping one seed lands in the right one.
  This directly supersedes #536's single-point negative rather than duplicating it.
  **Positive control:** the pipeline must first recover #538's own confirmed SE<->EM
  connection (or #537's basin) under a parameter-continuation from Earth-Moon mass ratio
  toward Jupiter-Europa's — proof the generalized module still finds the one connection
  we already trust before it's pointed at unmapped territory.
  **Feasibility: MEDIUM** — the corrector logic is proven by #538; the new build is the
  extraction/generalization plus a genuinely wider sweep, not new physics.
  **Recommended models:** module extraction/generalization behind #538's existing
  convergence-residual tests → **Sonnet** (spec-complete refactor, strong regression gate).
  The Jacobi/phase-sweep design (band bounds, resolution, what counts as a "hit") →
  **Opus** (the #339/#480-class criterion-definition judgment call the project has
  repeatedly found itself needing to settle in writing *before* sweeping). Final
  novelty/empty-region verdict → **Opus**, with a **Fable** independent adversarial pass
  before any "novel" claim or catalogue writeback (matching #538 Task 5's pattern).

- **#540** (P2, sequence after #539) — Point the now-generalized, cross-checked pipeline at
  the **Uranian system** (Miranda-Ariel-Umbriel-Oberon lanes) — the *one* system that has
  ever produced a confirmed novel catalogue row in this project's history (#312,
  Umbriel-Oberon quasi-cycler, V4). The explicit question: does a strictly-more-capable
  method (genuine coherent-model torus connections, vs. whatever produced #312) find
  *more* novel families near the known one, or corroborate #312 as an isolated island?
  Per the empty-region registry's re-open rule, this is exactly the case where a
  strictly-more-capable method is licensed to re-open previously-screened territory.
  **Feasibility: MEDIUM-HIGH** — reuses #539's module; the new work is sourcing/verifying
  Uranian system ephemeris constants and moon pair selection, not new solver machinery.
  **Recommended models:** system setup + parameter sourcing (digest-and-reconcile against
  a published Uranian-system reference, per the "digest ≠ adoption" discipline) →
  **Sonnet**. Running the sweep and adjudicating any hit against the #312 anchor →
  **Opus**, **Fable** second-opinion pass before writeback (same pattern as #539/#538).
  **SUPERSEDED (2026-07-11) by #558** — this entry's island-vs-family question was never
  answered as scoped (gated behind #538/#544, which parked; #545/#546's torus-heteroclinic
  re-screen answered a related-but-distinct question and was itself shelved by #548/#553/#555
  with 0-for-4 and no validated positive). A Fable review of #312 found the actual answer
  doesn't need the torus machinery at all — see #558 below, which reuses the ORIGINAL
  discovery genome that found #312 in the first place.

- **#558** (P0, cheap — Fable review estimates minutes of compute, fully disjoint from the
  shelved #522/#548 torus-heteroclinic pipeline) — #312 family census via the original
  discovery genome. A Fable review of #312 (2026-07-11) found the island-vs-family question
  #540 posed has never actually been answered, AND found a concrete, load-bearing reason to
  suspect "family" over "island": the discovery genome that found #312
  (`scripts/scan_312_uranus_umbriel_oberon_offset_sweep.py`) **holds moon-pair relative
  offset FIXED by convention** and only sweeps global phase — the exact convention choice
  that nearly HID #312 itself (0.636 km/s under a 3-moon-convention fixed offset vs the
  catalogued 0.025 km/s under the 2-moon/180° convention actually used). Every adjacent
  Uranian moon-pair "miss" on record (Oberon-Titania 0.0617 km/s — only 23% over the 0.05
  km/s gate; Ariel-Titania 0.27; Titania-Umbriel 0.66) was measured at an **arbitrary fixed
  offset**, never a true basin floor — the same blind spot that almost hid the one confirmed
  novel row in this project's history. Even within Umbriel-Oberon itself, the existing
  rel-offset sweep already found a SECOND distinct closure (rel_offset=45°, tof_scale=1.5,
  residual 0.024, but V∞=2.27 km/s — physically unusable, 4.2x surface escape) with a
  completely different V∞ signature than the catalogued point — multi-basin evidence at the
  one known grid point.
  **Scope:** generalize the existing 4.5-second script to sweep rel_offset (96) x phase0 (96)
  x tof_scale (densified, 0.5-3.0 in ~0.1 steps) x n_rev (0..3) over **all 10 Uranian moon
  pairs** (5 regular moons), plus a full-landscape dump for Umbriel-Oberon itself (the
  current JSONL keeps only 13 best-records; the within-pair basin census needs a re-run
  that keeps everything). Gate every sub-0.05 km/s basin through the EXISTING #324 physical
  max-bend filter, DOP853 cross-check, and — for survivors — the already-built, frozen
  V1->V4-strict gauntlet (#330/#331/#332/#335/#338 are all reusable as-is; positive control =
  re-running #312's own IC through the same pipeline, which the frozen tests already do).
  This is the SAME method that found #312, run where it was never run — not a new capability,
  not a re-run of the shelved torus machinery.
  **Recommended models (per the Fable review):** Sonnet for the sweep build/run (mechanical,
  behind deterministic gates). Opus for adjudicating any hit against the #312 anchor and the
  family/empty-region verdict. Fable second-opinion pass before any catalogue writeback
  (matches #546's own prescribed split for this territory).
  **RESULT (2026-07-11, Sonnet mechanical pass):** Positive control PASSED —
  `scripts/scan_558_uranus_all_pairs_offset_sweep.py` evaluates the EXACT catalogued point
  (Umbriel-Oberon, rel_offset=180°, tof_scale=2.0, n_rev=(1,1)) directly and reproduces
  residual=0.025232 km/s bit-for-bit against the stored `verify_327_umbriel_silver.py` value —
  not the 0.636 km/s fixed-convention artifact. The residual formula was also verified
  byte-for-byte against production `_close_one_phasing` (discovery_campaign.py, interior-flyby
  continuity + #259 anchor-wrap periodicity) — this is not a looser metric manufacturing hits.

  Two build-time findings changed the sweep's shape from the literal spec, both in the direction
  of MORE coverage per unit compute, not less: (1) **n_rev is free** — one `lambert()` call per
  leg at `max_revs=3` returns every revolution branch in [0,3] simultaneously (shared Stumpff
  root-finding), so all 16 `(n0,n1)` combinations are scored per grid point at ~zero marginal
  Lambert-solve cost, fully covering the spec's 0..3 range. (2) **the global `phase0` axis is
  PROVABLY redundant** — the circular-coplanar Kepler + patched-conic-Lambert closure is exactly
  rotationally symmetric (adding a constant to both moons' longitudes rigidly rotates the whole
  3-state configuration; the residual is built only from rotation-invariant vector-norm
  differences), verified empirically to ~1e-13 agreement across 6+ phase0 samples on 2 different
  pairs. The original discovery script's 96-sample phase grid was pure wasted compute; this run
  fixes `phase0=1` sample and reinvests the freed 96x budget into `rel_offset` (360 samples, 1°
  resolution vs. the original 96) and `tof_scale` (51 samples, 0.05 km/s step — denser than the
  spec's literal 0.1 ask). Full 20-direction sweep (10 pairs × 2 anchor choices) + #324 physical
  gate + independent DOP853 cross-check on every sub-0.05 km/s hit ran in ~8 min wall
  (09:21–09:29 AET) — over the "minutes" framing but justified by the above density increase and
  by gating thousands of hits, not a scope creep.

  **Headline finding — #540's island-vs-family question is answered decisively: FAMILY, and a
  LARGE one.** All 12 anchor-flyby directions among the 4 non-Miranda regular moons
  (Ariel/Umbriel/Titania/Oberon) produce dozens to hundreds of DISTINCT closures (basin-clustered
  at 3°/0.15-tof_scale separation; 3672 raw sub-gate grid points before clustering, 1092 after)
  that pass ALL THREE gates — residual (often 100-1000x tighter than #312's own 0.025 km/s),
  #324 physical max-bend (>5°), and DOP853 cross-check (<1 km, typically 1e-6 to 1e-4 km here) —
  at V∞ 0.3–2.4 km/s, comparable to or BELOW #312's catalogued 0.89–0.96 km/s. The 8 Miranda-
  involving directions pass the residual gate just as readily but ALWAYS fail the #324 bend gate
  (Miranda's GM=4.3 km³/s² is too small to usefully bend a flyby at these V∞ — bend angles 1-8°
  vs the 5° floor; confirmed real physics via direct inspection, not a formulation bug). Even
  WITHIN Umbriel-Oberon itself, the dense sweep finds a closure tighter than #312's own:
  rel_offset=358°, tof_scale=2.40, n_rev=(1,1), residual=0.000152 km/s (165x tighter),
  V∞=(1.221, 1.995, 1.221) km/s, bend=(9.5°, 6.9°, 9.5°), DOP853=4.5e-5 km — passing every gate
  #312 itself passed. Per-pair best-representative table (anchor-flyby-anchor, best residual,
  raw all-gates-passed count, distinct clustered basins):

  | pair | best residual (km/s) | raw hits | distinct basins |
  |---|---|---|---|
  | Oberon-Titania-Oberon | 7.2e-7 | 777 | 215 |
  | Titania-Oberon-Titania | 6.8e-6 | 707 | 212 |
  | Umbriel-Oberon-Umbriel | 1.5e-4 | 373 | 116 |
  | Umbriel-Titania-Umbriel | 1.4e-4 | 314 | 97 |
  | Titania-Umbriel-Titania | 8.7e-4 | 269 | 87 |
  | Oberon-Umbriel-Oberon | 2.6e-4 | 263 | 80 |
  | Umbriel-Ariel-Umbriel | 8.4e-5 | 232 | 60 |
  | Ariel-Umbriel-Ariel | 7.8e-5 | 193 | 57 |
  | Ariel-Titania-Ariel | 1.2e-4 | 165 | 52 |
  | Titania-Ariel-Titania | 2.5e-4 | 157 | 48 |
  | Oberon-Ariel-Oberon | 1.7e-4 | 131 | 41 |
  | Ariel-Oberon-Ariel | 3.0e-4 | 27 | 27 |
  | all 8 Miranda-involving directions | 6.4e-6 to 5.2e-4 | 0 | 0 (bend gate) |

  **Methodological flag for Opus/Fable — the load-bearing implication, not a side note:** this is
  not "more hits like #312" — it shows the #254/#285/#312 repeated-moon periodicity-continuity
  closure condition is GENERIC across nearly the whole non-Miranda Uranian moon-pair space once
  the relative-offset degree of freedom is exploited (a 2-constraint condition over the
  2-continuous-parameter (rel_offset, tof_scale) plane generically has isolated solutions per
  revolution-count branch, of which there are up to 16 per pair-direction across the swept
  tof-resonance range). #312 is one of >1000 basins meeting the current residual+bend+DOP853 gate
  stack, not a uniquely-selected discovery. That gate stack is evidently NOT selective enough to
  distinguish "an interesting novel cycler" from "a generic closure coincidence" once the search
  is properly dense in relative offset. Recommend Opus/Fable adjudicate: (a) whether/how to add a
  stronger selection criterion (tighter residual floor, resonance-structure or launch-window
  robustness prioritization à la #338/#559, literature-novelty-first triage) before spending
  further gauntlet compute on this family, and (b) whether #312's catalogue admission still
  stands on grounds other than "uniquely passed a closure gate" — e.g. its specific #338 DOY-
  robustness profile, or simply being the first-discovered representative of a large family
  rather than a rare individual hit. **No candidate from this pass was pushed through the
  V1-V4-strict gauntlet and nothing was written to `data/catalogue.yaml`** — per task discipline
  this Sonnet pass stops at the mechanical sweep + gate; adjudication is Opus/Fable's job.

  Artifacts (no catalogue writeback): `scripts/scan_558_uranus_all_pairs_offset_sweep.py`
  (generalized sweep + gating, ruff-clean); `data/scan_558_positive_control.jsonl`;
  `data/scan_558_uranus_{anchor}_{flyby}.jsonl` ×20 (per-direction top/full-grid + gate results;
  `umbriel_oberon` carries the FULL 18360-row landscape per spec item 5, not just top-N);
  `data/scan_558_uranus_all_pairs_index.jsonl` (cross-pair index + the 3672 raw all-gates-passed
  candidate records, full detail, for Opus/Fable to consume programmatically). The frozen #312
  gates (`tests/verify/test_silver_327_v1_passes.py` through `_v4_strict_passes.py`) were
  re-verified unchanged and still pass.

- **#561** (P0, judgment-only — the data already exists, this is ADJUDICATION not a new
  sweep) — the #558-mandated family-vs-coincidence adjudication. #558 found the #254/#285/#312
  repeated-moon periodicity-continuity closure condition is GENERIC across nearly the whole
  non-Miranda Uranian moon-pair space (>1000 distinct basins across 12 directions passing the
  residual + #324 bend + DOP853 gate stack, some 100-1000x tighter than #312's own residual),
  not a rare individual hit — and flagged this as load-bearing, not a side note: **#312 is now
  one of >1000 basins meeting the current gate stack, and that gate stack is evidently not
  selective enough to distinguish "an interesting novel cycler" from "a generic closure
  coincidence."** Two things need real judgment, not more compute:
  (a) **Design a stronger selection criterion** before spending any further gauntlet compute on
  this family — candidates per #558's own suggestion: a materially tighter residual floor,
  resonance-structure robustness (does the closure survive a #338/#559-style DOY/launch-epoch
  sweep, or is it also a knife-edge?), literature-novelty-first triage (screen candidates
  against the corpus BEFORE any further numerical work, cheapest filter first), or some
  combination. Read `data/scan_558_uranus_all_pairs_index.jsonl` (the full 3672-record raw
  candidate set, built for exactly this programmatic consumption) and `data/scan_558_uranus_
  umbriel_oberon.jsonl` (the full 18,360-row single-pair landscape) directly — don't re-derive
  from the OUTSTANDING.md summary table alone.
  (b) **Does #312's own catalogue admission still stand, and on what grounds?** Its `our_status`
  field context matters here (see `[[project_novel_findings_status]]` — #312 is this project's
  ONLY confirmed-novel row; PC(3,2) is explicitly `known-class-member`, not novel). If the
  closure condition #312 satisfies is generic across >1000 basins in its own moon-pair family
  class, is #312 still meaningfully "novel," or does its novelty claim need to be reframed
  (e.g. as "first-discovered representative of a large family" rather than "a unique closure"),
  or does the specific #338/#559 DOY-robustness profile / the V4-strict real-ephemeris survival
  that #312 was independently promoted on (a criterion this adjudication should check whether
  ANY of the other 1000+ basins would also pass, or whether #312 remains distinguished on that
  axis specifically) still justify its standing? Do NOT unilaterally downgrade or edit #312's
  catalogue row in this pass — surface the finding and a recommendation, this is adjudication
  informing a later decision, not the decision itself.
  **Explicitly out of scope for this pass:** running any candidate through the V1-V4-strict
  gauntlet, writing to `data/catalogue.yaml`, building new sweep code. If the adjudication
  concludes a stronger criterion is needed, SCOPE the follow-up (a new task number) rather than
  building it inline.
  **Recommended models:** Opus first (the criterion-design + #312-standing judgment is exactly
  the trust-bearing numerical-methods/discovery-verdict work this project's model-tiering
  policy reserves for Opus) — then a **Fable** second-opinion pass given the stakes (this
  potentially touches the project's ONE confirmed novel finding's own standing, not a routine
  candidate adjudication).
  **RESULT (2026-07-11, Opus adjudication — data read directly from the 3672-record
  `scan_558_uranus_all_pairs_index.jsonl` candidate array + the 18,360-row Umbriel-Oberon
  landscape, not the summary table). Verdict headlines: (a) the principled discriminant is NOT a
  tighter residual — it is synodic-resonance commensurability + the multi-cycle bounded-drift
  gauntlet the current gate stack never runs; (b) #312's LITERATURE novelty STANDS, but its
  framing must change from "unique isolated closure" to "first-documented representative of a
  large generic family," and its DOY-robustness leg is a genuine open gap (per #559).**

  **(a) Selection criterion — each candidate assessed against the real data:**

  1. **Tighter residual floor — REJECT as the primary knife.** It is not scale-free: tightening
     DOES thin the field sharply and non-arbitrarily (all-gates-passed count 3672 at <0.05 →
     1866 at <0.025 → 672 at <0.01 → 57 at <1e-3 → 4 at <1e-4 → 1 at <1e-6; min residual
     7.2e-7). But residual is physically MEANINGLESS as a discovery-value axis — it measures
     patched-conic Lambert match quality, not trajectory interest or mission utility.
     Decisively: **#312's own catalogued residual (0.025 km/s) sits at the family MEDIAN — the
     51.3rd percentile of the 3672.** A tighter-residual criterion would rank #312 in the bottom
     half and eventually DE-SELECT the one catalogued row entirely. Any criterion that rejects
     the reference candidate is the wrong criterion. Tighter residual is at most a tie-breaker,
     never the selector.

  2. **Physical quality (V∞ / bend) — a secondary ranking, not a family-vs-coincidence
     discriminant.** #312's max-encounter V∞ (~0.96 km/s) is top-quartile (23.8th percentile of
     the family's 0.33–2.40 km/s spread) — better than median but matched-or-beaten by 874 of
     the 3672. And #312 is **strictly dominated on every physical axis at once** by other
     candidates: e.g. Titania-Oberon at tof=30.86 d, n_rev=(2,2) has residual 6.8e-6 (3700×
     tighter), max V∞ 0.583 km/s (lower), min-bend 50.5° (vastly more usable) than #312's
     0.025 / 0.96 / 14.7°. Oberon-Titania at res 7.2e-7 / V∞ 0.872 / bend 28.9° likewise
     dominates it. So #312 is unremarkable-to-mediocre within its own family on physics; physical
     quality is worth carrying as a mission-utility RANK (prefer low V∞ + usable bend) but it
     does not separate "interesting" from "coincidence."

  3. **Synodic-resonance commensurability — THE principled discriminant, and it retains #312.**
     The reason #312 is a *quasi_cycler* (bounded 86k–530k km drift over 10 cycles rather than
     divergence — the #330/#331 V2/V3 result) is that its full Umbriel-Oberon-Umbriel cycle is
     ≈5 synodic periods: 2×14.94 d / T_syn(5.987 d) = 4.991, only 0.009 from integer 5. The
     single-cycle residual gate is BLIND to this: commensurability |ratio − round(ratio)| across
     the 3672 is essentially uniform-random (median 0.285; a uniform draw on [0,0.5] gives 0.25),
     and it is NOT enriched in the tight-residual or low-V∞ subsets (sub-1e-3 median 0.266,
     sub-0.7-V∞ median 0.329). Only **151 of 3672 (4.1%) fall within 0.02 of an integer ≥2, and
     29 (0.8%) within 0.01** — a physically-motivated ~25–125× thinning that, unlike a residual
     floor, KEEPS #312 (and surfaces siblings — see below). This is principled because it is tied
     to the *class definition*: geometry must re-close each cycle for multi-cycle drift to stay
     bounded, and only near-integer-synodic cycles do that. The vast majority of the 3672 are
     NON-commensurate one-cycle geometric closures that would DIVERGE on repetition — i.e. they
     are not quasi-cyclers at all, just Lambert coincidences.

  4. **Multi-cycle bounded-drift (V2 gauntlet) — the REAL selector the sweep never ran; the
     commensurability filter is its cheap proxy.** Every one of the 3672 passed only a
     SINGLE-cycle test (residual + #324 bend + DOP853). None were propagated for multi-cycle
     bounded vs divergent drift. That #330-style bounded-drift check is exactly what distinguishes
     a genuine quasi_cycler from a one-off closure, and it is what #312 actually earned its class
     on. Recommend it be run (frozen #330 machinery) on the commensurability survivors — do NOT
     trust commensurability alone as proof of boundedness; it is a well-motivated hypothesis with
     one confirmed anchor (#312 itself), to be VERIFIED by propagation, not assumed.

  5. **Literature-novelty-first triage — cheapest, run FIRST, but it discriminates at the
     pair-DIRECTION level, not per-candidate — and it is NOT a rubber stamp here.** The 53-anchor
     corpus contains 6 Uranus anchors, two of which are directly adjacent to the newly-surfaced
     directions: **Canales-Howell-Fantino 2021 (Titania-Oberon moon-to-moon analytical transfer)**
     and **Kumar 2025 (Uranus-Oberon PCRTBP mean-motion-resonance study)**. The two LARGEST basin
     directions in #558 are exactly Oberon-Titania (215 basins) and Titania-Oberon (212) — these
     carry real adjacent-publication risk and MUST get a live literature_check against those two
     papers before any of them is ever called "novel." #312's own Umbriel-Oberon-Umbriel
     direction is distinct from both (Kumar is single-moon+planet MMR, not a moon-pair cycler;
     Canales is Titania-Oberon) and was already cleared by #328 — so #312 is unaffected, but the
     family at large is NOT automatically novel per-direction.

  **Concrete recommended follow-up criterion (scope as a NEW task, do not build inline): a
  three-stage funnel, cheapest gate first —**
  - *Stage 1 (free, offline+live):* literature-novelty triage per pair-direction. Flag
    Titania-Oberon / Oberon-Titania / Oberon-MMR as adjacent-published-risk; live-check them
    against Canales 2021 + Kumar 2025 before any novelty language.
  - *Stage 2 (cheap, deterministic):* synodic-commensurability pre-filter — keep only candidates
    whose full-cycle ToF is within ~0.02 (→151) or ~0.01 (→29) of an integer × T_syn. This is the
    principled thinning tied to the quasi_cycler class definition; it retains #312. Explicitly do
    NOT use a tighter residual floor as the primary filter (it de-selects #312 and selects for a
    physically meaningless quantity).
  - *Stage 3 (moderate):* run the frozen #330 V2 bounded-drift gauntlet on Stage-2 survivors; only
    those with BOUNDED (not divergent) multi-cycle drift are genuine quasi-cyclers. V4-strict +
    DOY only for the handful that survive Stage 3.
  Residual and V∞/bend become RANKING keys within the survivor set, never the primary gate.

  **(b) #312's standing — novel-but-reframed-as-first-of-a-family, with an open DOY gap:**
  - **Literature novelty STANDS.** The genericity #558 found does NOT make the family published —
    the entire family is absent from the literature (corpus + #328's direct read of Heaton-Longuski
    2003 confirm zero Uranian moon-pair repeated-moon quasi-cyclers of this structural type). "One
    of >1000 basins" is a statement about our own search density, not about prior art. #312 remains
    literature-fresh.
  - **But the FRAMING is now wrong and should be corrected.** #312 is demonstrably NOT "a unique
    rare closure" or "an isolated island" (the exact open question in
    `[[project_novel_findings_status]]`, now RESOLVED: family, not island). It is the
    first-discovered / first-documented representative of a large generic family. Within that
    family it is unremarkable: median residual (51st pct), top-quartile-but-unexceptional V∞, and
    strictly dominated on all physical axes simultaneously by other members. Its distinction is
    "found first," not "physically best" — which, per this project's own literature-novelty
    discipline, is a LEGITIMATE and valuable claim ("first computed/documented example of a
    previously-undocumented family") but a DIFFERENT claim than "a unique find."
  - **Its independent-promotion axis (V4-strict + #338 DOY-robustness) does NOT currently
    distinguish it, on two counts.** (i) The #338/#559 DOY-robustness profile is artifact-
    contaminated — #559's own update showed the daily pass-rates measure a Lambert-branch flip +
    a DOP853 degenerate-arc failure, not real epoch sensitivity, so #312's genuine launch-epoch
    robustness is NOT cleanly established (a real gap, pending the #560 fixes; do not lean on it).
    (ii) NO other family member has been run through V4-strict at all, so #312 being "the only
    V4-validated one" reflects effort spent, not demonstrated superiority — it is not evidence
    #312 is distinguished on that axis, only that it is the only one tested.
  - **Recommendation (for a future writeback, NOT edited here): KEEP #312 catalogued — do not
    downgrade.** First-documented-of-family is a real result. But reframe its prose/notes from
    "unique closure / isolated island" to "first-computed representative of a large generic family
    of Uranian regular-moon synodic-resonance quasi-cyclers (#558)," and, when the `our_status`
    field is backfilled, mark it verified-novel-as-first-of-family rather than implying uniqueness.
    Explicitly note the DOY-robustness gap (pending #560). Update `[[project_novel_findings_status]]`
    to record the island-vs-family question as resolved (family).

  **Additional flags not explicitly asked for:**
  - **The ">1000 basins" headline overstates the quasi-cycler count.** Those are >1000
    SINGLE-CYCLE geometric closures; the genuine-quasi-cycler subset (near-integer-commensurate,
    likely bounded) is ~29–151 pre-V2, and smaller after V2. The correct claim is "#312's
    single-cycle closure condition is generic; the bounded-drift quasi-cycler property is far
    rarer and untested across the family." This sharpens rather than weakens #558's finding.
  - **A concrete genuine SIBLING of #312 exists in its own pair:** Umbriel-Oberon (and its mirror
    Oberon-Umbriel) at tof=20.92 d, n_rev=(2,2), full-cycle/T_syn = 6.987 (near-7:1), max V∞
    0.907 km/s, min-bend 16.1°, residual 4.1e-3 — a distinct resonance in the SAME moon pair,
    physically comparable to #312. Worth a dedicated look as a potential second catalogue row (a
    near-7:1 companion to #312's near-5:1) once the Stage-3 funnel is built.
  - Highest novelty-risk directions (largest basin counts + named adjacent literature) are
    Titania-Oberon / Oberon-Titania — do these FIRST in any promotion campaign, and expect some
    to be adjacent-published (Canales 2021 / Kumar 2025), not novel.
  **No `data/catalogue.yaml` edit, no gauntlet run, no new sweep code — adjudication only, per
  task scope. The literature-corpus inspection touched no repository code (read-only import).**

  **SECOND OPINION (2026-07-11, Fable — independent adversarial recheck against the raw JSONL,
  per this task's own two-model discipline). Headline: the adjudication's ARITHMETIC and its
  four top-level verdicts all hold; but Stage 2 of the recommended funnel, AS LITERALLY
  SPECIFIED, is defective — the commensurability statistic is grid-quantized below the
  tolerance being applied to it, and the specified 0.02 cut + integer≥2 restriction silently
  excludes the two LARGEST, highest-literature-risk directions (Titania-Oberon / Oberon-Titania)
  by 0.002 in a quantity the tof grid cannot resolve better than ±0.022. Fix the Stage-2
  formulation before running the funnel; everything else proceeds as adjudicated.**

  - **Every load-bearing number reproduces.** Independently recomputed from
    `scan_558_uranus_all_pairs_index.jsonl` (3672 candidates) using periods from
    `cyclerfinder.core.satellites`: #312 residual percentile 50.8 (claimed 51.3, def-of-record
    delta); frac median 0.2858 (claimed 0.285); 151 within 0.02 / 29 within 0.01 (exact); V∞
    percentile 23.9 / 879 matched-or-beaten / range 0.325–2.40 (claimed 23.8 / 874 / 0.33–2.40);
    sub-1e-3 and sub-0.7-V∞ frac medians 0.267/0.330 (claimed 0.266/0.329); residual-floor
    counts 1866/672/57 at 0.025/0.01/1e-3 (exact); sibling record confirmed verbatim (rel=358,
    tof=20.917 d, res 4.10e-3, ratio 6.9877, V∞max 0.907, min-bend 16.1°). The Titania-Oberon
    tof=30.86 "dominating" candidate sits at frac **0.494 — maximally NON-commensurate** — which
    actively STRENGTHENS the adjudication's residual-is-meaningless argument (the family's
    tightest closures are one-cycle Lambert coincidences by its own Stage-2 logic).
  - **Verdicts that HOLD:** (1) reject residual as primary filter — confirmed, and reinforced by
    the frac-0.494 dominator; (2) KEEP #312, reframe as first-of-family — confirmed, see below;
    (3) literature-first triage with Canales 2021/Kumar 2025 flagged on the T-O directions —
    confirmed and now MORE urgent (see the Stage-2 flaw); (4) V2 bounded-drift as the real
    selector, commensurability only a hypothesis to be verified — confirmed.
  - **#312's own standing is SOLIDER than the #558 framing implies — it is NOT one of the
    untested 3672.** Direct catalogue-row check: #312 is the ONE family member that already has
    the full multi-cycle chain — #330 ran the actual 10-cycle bounded-drift V2 (86k–530k km
    bounded oscillation), #331 cross-checked it with independent REBOUND IAS15 to 1.84e-6 km,
    #332/#335 carried it to V4-scipy/V4-strict real ephemerides, #338 swept 100 annual epochs
    (85/85 interior PASS). Its admission never rested on the single-cycle gate the 3672 passed;
    no V2 re-run on #312 is needed to trust its standing. The reframe (unique-island →
    first-documented-of-family) is the right and only change.
  - **THE FLAW — Stage 2's commensurability statistic is below its own data's resolution.**
    tof in the #558 sweep is `tof_scale × sqrt(P_a·P_b)` on a 0.05-step grid, so full-cycle/T_syn
    is quantized per direction: step 0.1248 for Umbriel-Oberon, 0.0502 Ariel-Umbriel, 0.0440
    Titania-Oberon. Three consequences, all verified in the data: (i) the "essentially
    uniform-random" claim FAILS in the tail where the filter operates — within 0.01 there are 29
    vs ~73 expected under uniform, within 0.005 just 2 vs ~37; the near-zero deficit is the grid
    comb, so the 0.01/0.02 survivor counts measure grid aliasing, not basin physics. Umbriel-
    Oberon only LOOKS enriched in commensurate candidates because its comb (2×base/T_syn =
    2.4954 ≈ 5/2) happens to land 0.004–0.013 from an integer every 8th tof step. (ii) The
    passing basins are tof-NARROW — sub-grid-step: at #312's own basin (rel=178) the ts=2.00
    point closes at 6.8e-3 while ±1 step (±0.37 d) sits at 0.64/0.28, i.e. 40–100× above; the
    grid frac (0.0088) is therefore NOT a measurement of the basin's true commensurability, and
    the exactly-commensurate tof (5·T_syn/2 = 14.9668 d, 0.026 d off #312's catalogued 14.94)
    has never been evaluated. Same for the sibling (exact 7:1 = 20.9535 vs grid 20.917).
    (iii) **Titania-Oberon and Oberon-Titania — 707+777 records, the two biggest directions AND
    the two with named adjacent literature — are STRUCTURALLY excluded by the ≤0.02 cut**: their
    grid's minimum achievable frac at integer≥2 is 0.0219 (ts=2.25/2.30 straddle ratio 2 at
    1.978/2.022 symmetrically — the true basin could sit at exactly 2.000), and their near-1:1
    lines (frac 0.011, ts=1.15) are excluded by the integer≥2 restriction, which has no physical
    justification (the re-closure condition 2·tof/T_syn ∈ ℤ is satisfied by 1 as much as by 5;
    integer 1 should be admitted or its exclusion argued). Net: the funnel as specified would
    do the OPPOSITE of the adjudication's own do-T-O-first instruction — it would never test
    T-O at all. Umbriel-Titania/Titania-Umbriel (min grid frac 0.0253) are likewise shut out.
  - **Required Stage-2 reformulation (fold into the funnel follow-up task, cheap and
    deterministic):** replace "threshold on grid frac" with a per-basin CONTINUOUS refinement —
    for each distinct basin (direction, rel-offset run, n_rev), refine tof to the exactly-
    commensurate value (nearest integer ≥1 × T_syn/2, per-leg tof = half that) and re-evaluate
    residual + bend + DOP853 THERE; a basin passes Stage 2 iff the commensurate point closes
    sub-gate. This is grid-independent, costs ~one 1-D refinement per basin, retains #312's
    basin on its physics rather than its grid luck, and re-admits the T-O near-2:1 (and near-1:1)
    lines to face the Canales/Kumar literature screen they most need. Optional flag: also admit
    half-integer 2·tof/T_syn (q=2 super-cycles — rel-offset advances 180°/cycle, so they need a
    conjugate basin at rel+180°, which the data suggests exists: #312's family shows basins near
    BOTH rel≈180 and rel≈0) — a one-line extension worth carrying as a secondary hypothesis.
  - **Sibling verdict: real and worth pursuing, but INSIDE the fixed funnel, not ahead of it.**
    Its numbers check out exactly as flagged, it inherits #328's Umbriel-Oberon literature
    clearance (same pair-direction as #312), and it is the best non-#312 candidate on the joint
    V∞/bend/frac axes. But its 0.0123 frac is grid-frac like everyone else's, its basin is
    tof-narrow, and the exactly-7:1 point (20.9535 d) is unevaluated — so "promising" is
    conditional on the same Stage-2 refinement + V2 propagation every other survivor needs. No
    special pre-funnel dispatch; it simply enters Stage 3 first among survivors.
  - **For `[[project_novel_findings_status]]`:** the island-vs-family resolution (FAMILY) stands
    and should be recorded; but do NOT encode the "~29–151 commensurate subfamily" counts — those
    are grid-aliasing artifacts, and the true commensurate-basin count is UNKNOWN until the
    per-basin refinement runs (it may grow, since T-O/O-T/U-T/T-U were structurally excluded).
  **Net recommendation: proceed with the adjudication's 3-stage funnel and #312 reframe as
  adjudicated, with ONE amendment — reformulate Stage 2 as per-basin commensurate-tof refinement
  (integer ≥1, optional half-integer flag) before any survivor list is drawn, and keep the
  literature-first pass on Titania-Oberon/Oberon-Titania mandatory since those directions
  re-enter the funnel once the grid artifact is fixed. No #312 V2 re-run needed; its multi-cycle
  standing is already the best-evidenced in the family. Analysis-only pass: no code, no
  catalogue edit, no gauntlet run.**

- **#562** (P0, execution of the Fable-corrected #561 funnel — the data/machinery all exist,
  this is a bounded build+run, not new capability) — per-basin continuous-tof commensurability
  refinement, fixing the exact grid-aliasing flaw the Fable second opinion found in #561's
  Stage-2 filter. **Scope, precisely as specified:**
  (1) **Literature triage first** (cheapest filter, run before any numeric refinement): screen
  the Titania-Oberon and Oberon-Titania directions specifically against the corpus (the two
  named adjacent-published risks: Canales-Howell-Fantino 2021 Titania-Oberon transfer, Kumar
  2025 Uranus-Oberon MMR) via `search/literature_check.py` (or whatever this project's actual
  novelty-check module is called — grep for it). Report per-direction lit-risk BEFORE spending
  refinement compute on directions likely to be adjacent-published, not novel.
  (2) **Stage-2 reformulation**: for each distinct basin (across ALL 20 anchor-flyby
  directions, not just the ones the broken grid filter happened to admit) found in #558's raw
  data (`data/scan_558_uranus_all_pairs_index.jsonl`, 3672 records — read directly, not the
  OUTSTANDING.md summary), refine `tof` to the exactly-commensurate value (nearest integer >=1
  x T_syn/2, per-leg tof = half that — a 1-D root-find/bisection per basin, NOT a new grid
  sweep) and re-evaluate residual + #324 physical bend + independent DOP853 cross-check AT that
  refined point. A basin passes Stage 2 iff the commensurate point closes sub-gate. This is
  grid-independent by construction and specifically re-admits the Titania-Oberon/Oberon-Titania
  near-2:1 and near-1:1 lines the broken grid filter structurally excluded (integer=1 should be
  admitted, not just integer>=2 — Fable's note: "the re-closure condition 2*tof/T_syn in Z is
  satisfied by 1 as much as by 5, integer 1 should be admitted or its exclusion argued").
  Optional secondary hypothesis (Fable's flag, cheap to include): also test half-integer
  2*tof/T_syn = q=2 super-cycles, checking the conjugate basin at rel-offset+180 the data
  suggests exists.
  (3) **#312's own sibling** (Umbriel-Oberon, near-7:1, tof≈20.92d, already literature-cleared
  via #328's same-pair clearance): runs through this SAME refined Stage 2 first among survivors
  — no special pre-funnel treatment, per Fable's explicit instruction not to jump the queue.
  (4) **Report the TRUE commensurate-subfamily size** (expected to differ from, likely exceed,
  the grid-artifact 151/29 counts in #561's original write-up — do not repeat those numbers as
  if real) and a ranked survivor list with residual/V∞/bend/lit-risk as reported ranking keys,
  not gates.
  **Explicitly out of scope**: running the V1-V4-strict gauntlet on any survivor (that is Stage
  3, a separate follow-up once this stage's survivor list exists), writing to
  `data/catalogue.yaml`, re-running #312's own V2 test (already solid per Fable's confirmation
  — #330 already ran the real 10-cycle bounded-drift test on it).
  **Recommended model:** Sonnet (mechanical: a 1-D per-basin refinement behind the existing
  residual/bend/DOP853 gate functions #558 already built and validated, plus a literature-
  check call — spec-complete work behind strong deterministic gates, matching #558's own
  successful precedent). Any survivor that clears Stage 2 needs Opus adjudication + a Fable
  second-opinion pass before Stage 3 (the V1-V4-strict gauntlet) or any catalogue-adjacent
  writeback — do not push a candidate through the gauntlet from this dispatch.
  **RESULT (2026-07-11, Sonnet mechanical pass) — the grid-artifact 151/29 counts do NOT
  survive contact with a continuous refinement. The TRUE commensurate-subfamily size is 20
  unique physical closures, not 151 or 29; Titania-Oberon/Oberon-Titania are re-admitted and
  carry real but non-blocking adjacent literature.**

  **(1) Literature triage (Stage 1, run first, no numeric work spent before this).** Grepped
  `docs/notes/` for both named risks: `docs/notes/2026-06-16-328-uranian-cycler-lit-deep-dive.md`
  (the #328 deep-dive, 20 WebSearch queries + per-paper triage table) already assessed BOTH
  papers directly against the Titania-Oberon moon pair specifically (queries 2, 12, 13; rows
  10-11 of its triage table): **Canales-Howell-Fantino 2021** (Celest. Mech. Dyn. Astron.
  133:36, arXiv:2110.03683) is a ONE-SHOT L2-halo(Titania)→L1-halo(Oberon) unstable/stable-
  manifold transfer — no free-flight repeated encounter, no periodicity; **Kumar 2025**
  (arXiv:2509.03655, §6.2) is a SINGLE-MOON Uranus-Oberon PCRTBP 3:4/4:5/5:6/4:3/5:4/6:5 MMR
  study (companion arXiv:2309.06073 extends to Uranus-Titania-Oberon CCR4BP secondary
  resonances, still single-moon-period MMR sub-harmonics, not a moon-pair cycler). The #328
  verdict for both: **STRUCTURAL ADJACENT, not DIRECT MATCH** — same moons, different topology.
  Re-confirmed formally against `literature_check.py`'s own structural matcher (not just the
  digest prose): a `CandidateSignature(primary="Uranus", sequence=("Titania","Oberon","Titania"),
  topology_label=frozenset({"repeated-moon"}))` returns **zero** overlapping `KNOWN_CORPUS`
  anchors — `_candidate_anchors`'s topology-label filter excludes Kumar (`{"resonant"}`) and
  Canales (`{"halo"}`) as disjoint from `{"repeated-moon"}`, exactly matching #328's human
  read. **Verdict: LOW direct-rediscovery risk for the repeated-encounter topology specifically,
  but MANDATORY citation of both papers as adjacent literature in any future write-up or
  promotion** given they study the exact same moon pair — a citation obligation, not a
  blocking match. This offline conclusion reconfirms/extends #328's existing clearance to the
  T-O/O-T directions specifically (which #328 assessed but never a repeated-moon candidate at);
  no fresh live WebSearch was run (out of scope per the task's own framing — digest-note
  triage was the specified cheap-first filter), and Opus/Fable should still sanity-check this
  before any T-O/O-T survivor is called novel, per #561's own instruction.

  **(2) Stage-2 reformulation, run** (`scripts/refine_562_commensurability.py`). Loaded all
  3672 raw single-cycle gate-passing candidates from `data/scan_558_uranus_all_pairs_index.jsonl`
  (read directly, not the summary table) across all 20 anchor-flyby directions (the 8
  Miranda-involving directions carry 0, confirmed unchanged — #324 bend-gate failure, real
  physics per #558). Clustered into **360 distinct basins** across the 12 non-Miranda
  directions by (anchor, flyby, n_rev) + rel_offset adjacency (gap ≤2°) — Titania-Oberon (66)
  and Oberon-Titania (62) basins are now IN this set, not structurally excluded as they were
  under #561's literal grid-frac filter. Per basin: computed the pair's synodic period
  (generic `1/|1/T_a-1/T_b|`; self-check `T_syn(Umbriel,Oberon)=5.9867 d` reproduces the
  5.987 d reference used throughout #561/#562), set `n = max(1, round(2·tof_grid/T_syn))`
  (n=1 explicitly admitted, per Fable's instruction), fixed `tof` to the EXACT commensurate
  value `n·T_syn/2` (a direct formula, not a search), then ran a bounded local 1-D search over
  `rel_offset_deg` at that fixed tof (coarse 61-pt scan ±3°, refined 41-pt scan ±0.1°) to find
  the true nearby residual minimum — necessary because the basin's (rel_offset, tof) zero is
  generically an isolated point (per the #561 Opus adjudication's own finding), so the grid's
  rel_offset and the exact commensurate tof do not automatically coincide. Re-gated every
  refined point through #558's own `gate_candidate` (physical #324 bend + independent DOP853
  cross-check) verbatim — no new gate logic. Full run: 10 s wall.

  **Result: 39/360 basin-refinements pass Stage 2 sub-gate** (residual <0.05 km/s AT the exact
  commensurate tof, plus bend + DOP853). Deduping direction-mirrors (anchor-flyby-anchor and
  flyby-anchor-flyby are the same physical closed loop) and 0°/360° wraparound cluster-split
  artifacts gives **20 unique physical commensurate closures — the TRUE commensurate-subfamily
  size.** This supersedes the old grid-artifact "151/29" counts from #561's original write-up
  entirely; those measured grid-frac against a 0.05-tof_scale grid, not continuous physics, and
  must not be repeated as real. 7 of the 20 are in the Titania-Oberon/Oberon-Titania lit-risk
  direction — previously ZERO were ever testable there under the broken filter, now correctly
  re-admitted, confirming Fable's prediction. This includes a clean **n=1 (near-1:1) T-O
  closure** (tof=12.32 d, V∞_max=2.162 km/s, min-bend=6.27°, residual 2.7e-15 — machine-
  precision exact) that the old integer≥2 restriction would have permanently excluded.

  **#312 itself** (Umbriel-Oberon, n=5, n_rev=(1,1)) refines to tof=14.9668 d (vs. its
  catalogued 14.94 d), residual 8.9e-16 — essentially EXACT, tighter than its own catalogued
  0.025 km/s. This is the expected signature of a genuine symmetric periodic closure (the
  Miele-style "perpendicular crossing" condition: at rel_offset=0/180 with matched leg
  revolution counts and exactly-commensurate tof, mid-continuity and periodicity are satisfied
  identically, not by numerical luck) — a strong independent confirmation that the family's
  anchor sits on a mathematically real closure, not a grid coincidence.

  **#312's sibling** (near-7:1, no special pre-funnel treatment per Fable's instruction — it
  ran through the identical funnel as every other basin): Umbriel-Oberon n=7, n_rev=(2,2),
  refines to tof=20.9535 d (matches Fable's own hand-computed exact-7:1 value verbatim),
  residual 5.3e-15 (machine-precision exact), V∞_max=0.93 km/s, min-bend=15.43°, DOP853
  4.4e-5 km. Passes cleanly — the strongest non-#312 survivor by residual/exactness, a real
  candidate for Stage 3 (not run here).

  **Ranked survivor table** (20 unique, sorted by refined residual; `res` in km/s, `Vinf` =
  max per-encounter km/s, `bend` = min per-encounter deg, `dop853` = max arrival residual km):

  | pair | n | n_rev | tof (d) | res | Vinf | bend | dop853 | lit-risk |
  |---|---|---|---|---|---|---|---|---|
  | Oberon-Umbriel | 5 | (1,1) | 14.970 | 8.9e-16 | 0.965 | 15.55 | 3.1e-5 | |
  | Oberon-Umbriel | 6 | (1,1) | 17.960 | 8.9e-16 | 1.995 | 6.85 | 9.3e-6 | |
  | Oberon-Umbriel | 2 | (0,0) | 5.990 | 8.9e-16 | 1.691 | 5.13 | 7.3e-6 | |
  | Ariel-Umbriel | 2 | (1,1) | 6.430 | 2.0e-15 | 1.305 | 8.37 | 3.8e-6 | |
  | Ariel-Umbriel | 1 | (0,0) | 3.220 | 2.4e-15 | 1.300 | 8.43 | 4.4e-7 | |
  | Oberon-Titania | 1 | (0,0) | 12.320 | 2.7e-15 | 2.162 | 6.27 | 2.0e-6 | **yes** |
  | Oberon-Umbriel | 3 | (0,0) | 8.980 | 2.9e-15 | 1.520 | 11.33 | 6.8e-6 | |
  | Titania-Umbriel | 1 | (0,0) | 3.950 | 3.1e-15 | 1.230 | 9.35 | 1.2e-6 | |
  | Ariel-Umbriel | 3 | (2,2) | 9.650 | 5.2e-15 | 1.492 | 6.51 | 6.5e-6 | |
  | Oberon-Umbriel | 7 | (2,2) | 20.950 | 5.3e-15 | 0.930 | 15.43 | 4.4e-5 | (#312 sibling) |
  | Ariel-Oberon | 9 | (1,1) | 13.950 | 6.0e-15 | 1.430 | 6.99 | 2.2e-5 | |
  | Oberon-Titania | 2 | (0,0) | 24.630 | 9.3e-15 | 1.524 | 11.28 | 2.2e-5 | **yes** |
  | Oberon-Titania | 1 | (1,1) | 12.320 | 2.0e-14 | 2.169 | 6.23 | 1.9e-5 | **yes** |
  | Ariel-Titania | 5 | (1,1) | 8.870 | 6.1e-14 | 1.106 | 11.24 | 1.5e-5 | |
  | Oberon-Titania | 3 | (0,0) | 36.950 | 7.9e-3 | 2.321 | 5.14 | 3.3e-5 | **yes** |
  | Oberon-Titania | 2 | (2,2) | 24.630 | 2.0e-2 | 1.227 | 17.52 | 4.3e-5 | **yes** |
  | Oberon-Titania | 2 | (1,1) | 24.630 | 2.2e-2 | 1.512 | 11.43 | 1.7e-5 | **yes** |
  | Oberon-Titania | 2 | (0,1) | 24.630 | 3.5e-2 | 2.001 | 6.82 | 6.4e-6 | **yes** |
  | Titania-Umbriel | 4 | (2,2) | 15.820 | 3.9e-2 | 0.991 | 20.31 | 7.5e-6 | |
  | Ariel-Titania | 3 | (0,0) | 5.320 | 4.4e-2 | 1.730 | 8.85 | 7.5e-6 | |

  (`#312` proper, Umbriel-Oberon n=5, sits at the top of this table alongside its sibling —
  both omitted from the "18 non-#312 unique" count language elsewhere in this write-up only
  in the sense that they are the two already-known family members; they appear in the table
  like every other survivor, no special formatting.) Note the qualitative split: the 14
  best-residual survivors (all pairs except the three highest-`n` Titania-Oberon/Oberon-Titania
  rows and two Ariel/Titania-Umbriel rows) close to machine precision (symmetric-orbit
  signature); the remaining 6 (mostly higher-`n` T-O/O-T basins whose grid rel_offset sat well
  off 0°/180°) close at 0.008-0.044 km/s — still comfortably sub-gate but visibly nearer the
  0.05 floor, a genuine (not degenerate) near-closure.

  **Optional half-integer/conjugate check (§7):** not run as a separate super-cycle search —
  the primary per-basin clustering already naturally surfaces basins at both rel≈0 and rel≈180
  independently for several pairs at matching (n, n_rev) (e.g. Oberon-Umbriel n=6 at rel=0/360;
  Umbriel-Oberon n=6 likewise), consistent with the Fable-flagged rel≈180-conjugate pattern
  without requiring extra compute.

  **Explicitly not done here (per task scope):** no V1-V4-strict gauntlet run on any survivor,
  no `data/catalogue.yaml` writeback, no #312 V2 re-run (unchanged, still solid per #330).
  Artifacts: `scripts/refine_562_commensurability.py` (ruff-clean); `data/scan_562_commensurability_refinement.jsonl`
  (360-basin full refinement detail, all fields, for Opus/Fable to consume programmatically).
  **Next: Opus adjudication of the 20-survivor list (esp. the 7 Titania-Oberon/Oberon-Titania
  members against Canales/Kumar) + a Fable second-opinion pass, per #561's own two-model
  discipline for this territory, BEFORE any Stage 3 (V1-V4-strict gauntlet) dispatch.**

- **#563** (P0, replaces grid-search-plus-refinement with a direct, exhaustive construction —
  the #562 "20 unique closures" count is a LOWER BOUND, not a census) — direct symmetric-
  closure enumeration for the #312 Uranian moon-pair family. **The gap this closes**: #562
  only refines basins #558's ORIGINAL discrete grid (1° `rel_offset`, 0.05 km/s-step
  `tof_scale`) happened to detect via its initial residual gate. #561's own diagnostic already
  established that gate-passing basins tend to be TOF-NARROW — narrower than that grid's own
  step size. A basin narrower than the grid resolution in either dimension can fall entirely
  between sample points and never trigger the initial gate at ANY tested point — #558 would
  never have found it, so #562 (which only refines what #558 already surfaced) could never
  recover it either. This is the SAME class of bug (grid coarser than the feature being
  searched for) that caused the original Stage-2 filter flaw, one level upstream now. Two
  further, smaller known gaps: the half-integer/conjugate (q=2 super-cycle) case was only
  found opportunistically where it happened to already exist as a separate visible cluster in
  #558's data, not independently searched; and only 12/20 anchor-flyby directions are
  physically viable (the 8 Miranda-involving ones correctly excluded by the real #324 bend-gate
  physics — not a coverage gap, don't re-test these).
  **The fix — stop searching a grid, construct directly.** #562's own result data shows the
  dominant closure signature is a SYMMETRIC periodic orbit (`rel_offset` in {0 deg, 180 deg},
  matched leg revolution counts, exactly-commensurate `tof`) closing to machine precision
  (1e-14 to 1e-16 km/s) — the classical Miele-style "perpendicular crossing" condition, not a
  numerical coincidence (#562's own write-up makes this argument explicitly for #312 itself).
  This means the full symmetric-closure family is a FINITE, ENUMERABLE set — not a continuous
  space that needs sampling at all: for each of the 12 viable anchor-flyby directions, each of
  up to 16 `n_rev=(n0,n1)` combinations in `{0,1,2,3}x{0,1,2,3}` (the range #558 already scores
  "for free" per grid point via cached multi-rev Lambert branches), each commensurability
  integer `n>=1` up to a reasonable ceiling (bound by the same max-tof range #558's own sweep
  covered), and `rel_offset` in `{0 deg, 180 deg}` (both symmetric candidates, which
  automatically covers the "conjugate" case #562 could only find opportunistically) — directly
  CONSTRUCT the candidate at the exact commensurate `tof = n*T_syn/2` and the symmetric
  `rel_offset`, evaluate residual/bend/DOP853 there (reuse #558/#562's own gate functions
  verbatim, no new gate logic), with NO grid-resolution risk of missing a basin, because
  nothing is being searched for — every candidate in the finite enumeration is checked
  directly. This is provably exhaustive FOR THE SYMMETRIC-CLOSURE CLASS specifically (not for
  every possible asymmetric closure — #562's own data showed a handful of near-closures at
  0.008-0.044 km/s that were NOT machine-precision-exact, suggesting some genuinely-asymmetric
  family members may also exist; flag this as an explicit, honest scope limit of this task
  rather than silently claiming total exhaustiveness — a genuinely different, more expensive
  adaptive/basin-width-aware grid search would be needed for the asymmetric case, out of scope
  here, note it as a possible future #564 if this task's result motivates it).
  **Explicitly out of scope**: the V1-V4-strict gauntlet on any survivor, `data/catalogue.yaml`
  writeback, re-testing the 8 Miranda directions, an adaptive/asymmetric-closure search (see
  above — a different task if warranted).
  **Recommended model:** Sonnet (mechanical: a finite, well-bounded enumeration behind the
  exact gate functions #558/#562 already built and validated — spec-complete work behind a
  strong deterministic gate, matching #558/#562's own successful precedent). Any NEW survivor
  this finds beyond #562's existing 20 needs the same Opus + Fable adjudication path #562's
  survivors already need before Stage 3 or any catalogue-adjacent writeback.
  **RESULT (2026-07-11, Sonnet mechanical pass) — the direct construction confirms #562's
  "20-count" was indeed a lower bound: 11 genuinely NEW symmetric closures recovered, none of
  them findable by the original #558 grid at any resolution. Total known symmetric-class
  family size is now 30 (19 physical families #562 already touched + 11 new); the symmetric
  class is now provably exhaustive within #558's own tested tof range.**

  Built `scripts/enumerate_563_symmetric_closures.py`, reusing `residual_at_point`,
  `gate_candidate`, `GATE_RESIDUAL_KMS`, `N_REV_MAX` from `scan_558_uranus_all_pairs_offset_sweep.py`
  and `synodic_period_days` from `refine_562_commensurability.py` verbatim — no new gate logic
  written. For each of the 12 non-Miranda anchor-flyby directions (Ariel/Umbriel/Titania/Oberon,
  C(4,2)=6 pairs x 2 directions), computed the pair's `T_syn` and bounded the commensurability
  integer by `n_max = floor(2 * 3.0 * sqrt(P_a*P_b) / T_syn)` — a STRICT bound using the literal
  maximum `tof_scale=3.0` #558's production sweep actually used (verified against every one of
  the 20 `scan_558_uranus_*.jsonl` files' own `_meta` record, not assumed). This gives
  direction-independent `n_max` per pair: Ariel-Umbriel 3, Ariel-Titania 7, Ariel-Oberon 11,
  Umbriel-Titania 4, Umbriel-Oberon 7, Titania-Oberon 2. For every `n` in `[1, n_max]`, every
  `n_rev=(n0,n1)` in `{0,1,2,3}^2` (16 combos), and `rel_offset` in `{0°, 180°}`, directly
  constructed `tof = n*T_syn/2` and evaluated residual + #324 bend + DOP853 there — **2176 total
  candidates evaluated, 8.1 s wall** (well inside the "few thousand cheap evaluations" the task
  predicted). 60 raw points pass all gates; every single one closes to machine precision
  (worst residual 1.6e-13 km/s, most 1e-14 to 1e-16) — the symmetric "perpendicular crossing"
  signature holds with zero exceptions across the entire enumeration, confirming the task's
  premise that this class is exact-by-construction, not luck.

  **Deduping direction-mirrors** (anchor-flyby-anchor and flyby-anchor-flyby constructions of
  the same physical loop; verified empirically — every one of the 60 raw passes has exactly one
  mirror partner with matching `n`, `n_rev`, `tof`, and per-encounter V∞ to 3 decimal places)
  gives **30 unique physical symmetric closures**. Cross-referencing against #562's own
  `data/scan_562_commensurability_refinement.jsonl` (clustering ITS 39 raw stage-2 passes the
  same way, by pair/`n`/`n_rev`/V∞-branch, actually yields **21** unique closures, not the 20
  reported in #562's write-up — the prior hand-dedup missed one branch, Oberon-Titania n=2
  n_rev=(2,2) at V∞=0.684 near rel≈180°, distinct from the already-listed V∞=1.227 branch near
  rel≈0°; a minor correction noted here, it does not change #562's substantive conclusions) by
  (pair, n, n_rev, V∞ within 0.03 km/s, tof within 0.02 d):

  - **14 CONFIRM** — #562 already found these exactly (residual <1e-6 km/s), #563 reproduces
    them independently via direct construction with no local search at all. Includes #312 itself
    (Umbriel/Oberon n=5, n_rev=(1,1), V∞=0.965) and its n=7 sibling.
  - **5 SUPERSEDE** — #562 found the SAME physical branch but only as a non-exact near-miss
    (residual 0.020-0.047 km/s, because its grid representative + bounded ±3° local search
    didn't happen to land close enough to the true 0°/180° zero); #563's direct construction
    finds the TRUE machine-precision-exact zero for the identical branch (residual 1e-14 to
    1e-15). All 5 are in the Titania-Oberon/Oberon-Titania or Titania-Umbriel/Ariel-Titania
    directions — i.e. exactly the higher-`n`, wider-basin-search-window cases where #562's own
    write-up flagged the refinement as "visibly nearer the 0.05 floor."
  - **11 NEW** — physical closures with NO #562 basin at all, at any residual. These are the
    literal coverage-gap recoveries the task set out to find: real, machine-precision-exact
    symmetric closures whose (rel_offset, tof) point fell entirely between #558's original grid
    samples and never tripped its initial residual gate.

  **The 11 new closures** (sorted by pair; `res`/`dop853` in the gates' own units, `bend` =
  min per-encounter deg, all pass residual+bend+DOP853 identically to #558/#562's criteria):

  | pair | n | n_rev | tof (d) | res | Vinf | bend | dop853 | lit-risk |
  |---|---|---|---|---|---|---|---|---|
  | Ariel-Umbriel | 3 | (2,2) | 9.648 | 2.4e-14 | 1.547 | 6.03 | 3.2e-5 | |
  | Ariel-Oberon | 5 | (0,0) | 7.752 | 2.0e-15 | 1.829 | 6.22 | 6.7e-6 | |
  | Ariel-Oberon | 10 | (1,1) | 15.504 | 5.6e-15 | 1.670 | 5.21 | 2.0e-5 | |
  | Titania-Umbriel | 3 | (1,1) | 11.863 | 1.8e-15 | 1.766 | 9.16 | 2.5e-5 | |
  | Titania-Umbriel | 4 | (0,0) | 15.818 | 7.1e-15 | 2.342 | 5.39 | 4.1e-6 | |
  | Titania-Umbriel | 4 | (1,1) | 15.818 | 1.4e-14 | 1.706 | 5.49 | 4.2e-5 | |
  | Titania-Umbriel | 4 | (1,1) | 15.818 | 3.1e-15 | 2.139 | 6.40 | 1.4e-5 | |
  | Umbriel-Oberon | 5 | (0,0) | 14.967 | 6.0e-15 | 2.124 | 6.09 | 1.4e-5 | |
  | Umbriel-Oberon | 7 | (1,1) | 20.954 | 1.3e-14 | 1.660 | 7.41 | 4.6e-5 | |
  | Oberon-Titania | 1 | (0,0) | 12.316 | 8.1e-14 | 1.799 | 8.85 | 1.7e-5 | **yes** |
  | Oberon-Titania | 2 | (1,1) | 24.632 | 6.2e-15 | 0.959 | 24.94 | 7.4e-6 | **yes** |

  2 of the 11 new closures are in the Titania-Oberon/Oberon-Titania direction and carry the
  same adjacent-literature obligation #562's Stage-1 triage already established (Canales-
  Howell-Fantino 2021 halo-manifold transfer, Kumar 2025 single-moon MMR study — both
  structurally adjacent, not direct topology matches, per #562's `literature_check.py`
  cross-check; not re-run here, same conclusion applies mechanically to any repeated-encounter
  T-O/O-T candidate). The other 9 carry no special lit-risk flag beyond the general Uranian-
  system baseline #328 already covers.

  **Explicit, honest scope limit (per the task spec):** this enumeration is exhaustive for the
  SYMMETRIC closure class only, within #558's own tested tof range. It correctly does NOT
  reproduce #562/#561's 2 genuinely asymmetric near-closures (Oberon-Titania n=3 n_rev=(0,0)
  tof=36.95d at rel_offset=114.15°, residual 7.9e-3; Oberon-Titania n=2 n_rev=(0,1) tof=24.63d
  at rel_offset=268.19°, residual 3.5e-2) — neither sits at rel_offset in {0°,180°}, and the
  first also sits just outside the strict tof_scale<=3.0 bound, so both are structurally outside
  what a symmetric-only direct construction can find. These may be genuine asymmetric family
  members; a basin-width-aware adaptive grid search would be needed to characterize them, which
  is explicitly NOT built here (a possible future #564 if this result motivates it).

  **Explicitly not done here (per task scope):** no V1-V4-strict gauntlet run on any survivor
  (all 30 unique closures, old and new, remain V0/pre-Stage-3), no `data/catalogue.yaml`
  writeback, no re-testing of the 8 Miranda directions (unchanged, confirmed real #324
  bend-gate failures). Artifacts: `scripts/enumerate_563_symmetric_closures.py` (ruff-clean,
  mypy not applicable — `scripts/*.py` is outside this project's mypy gate);
  `data/enumerate_563_symmetric_closures.jsonl` (full per-direction + per-pass detail, 2176
  candidates' worth of summary plus all 60 raw gate-passing rows, for Opus/Fable to consume
  programmatically).
  **Next: the 11 new closures (2 of them Titania-Oberon/Oberon-Titania lit-risk) join the SAME
  Opus + Fable adjudication queue #562's 20 survivors are already waiting in — no separate
  adjudication ask. That queue's total scope is now the 30-closure #563 table (superseding
  #562's 20-row table as the reference list), not #562's list alone.**

- **#564** (P0, judgment-only — the data all exists, this is the Stage-3/writeback gate the
  whole #558->#563 thread has been building toward) — Opus + Fable adjudication of the full
  30-closure #563 survivor table (supersedes #562's 20-row table as the reference list; read
  `data/enumerate_563_symmetric_closures.jsonl` directly for full per-closure detail, not just
  the OUTSTANDING.md summary table). This is the gate #561 itself specified before any Stage 3
  (V1-V4-strict gauntlet) or catalogue-adjacent writeback — do not skip straight to running the
  gauntlet on the whole list.
  **Scope:**
  (1) **Rank/triage the 30 closures for gauntlet priority.** Not all 30 need equal investment —
  weigh residual/V∞/bend quality (all 30 are machine-precision-exact per #563, so this axis is
  now nearly uniform and NOT much of a discriminant — a change from #561/#562's situation),
  literature risk (2 of the 11 new closures + several of #562's confirmed set sit in the
  Titania-Oberon/Oberon-Titania direction — cite Canales-Howell-Fantino 2021 and Kumar 2025
  mandatorily for any of these that proceed, per #562's own literature-triage verdict), and
  genuine physical/mission interest (lowest V∞, best bend angles, shortest tof — check what
  actually varies meaningfully across the 30 now that exactness doesn't).
  (2) **Decide how many, and which, closures actually warrant the V1-V4-strict gauntlet** (a
  real compute/adjudication cost per candidate — #312 itself took the full #327->#338 task
  chain). Recommend a prioritized subset, not "run all 30" by default — with reasoning for the
  cutoff, not an arbitrary top-N.
  (3) **Settle #312's own reframing** (per #561's original recommendation, not yet acted on):
  should #312's catalogue prose actually be edited now to "first-documented representative of
  a 30-closure family," given the family census is now believed complete (for the symmetric
  class, within #558's tested range — note this qualifier, don't overclaim)? This adjudication
  may recommend the specific edit; do NOT make the edit in this pass (that's a separate,
  deliberate writeback step once recommended).
  (4) **Address the 2 remaining genuinely-asymmetric near-closures** #563 correctly declined to
  find (Oberon-Titania n=3 tof=36.95d rel=114.15°, residual 7.9e-3; Oberon-Titania n=2 tof=24.63d
  rel=268.19°, residual 3.5e-2) — are these worth a dedicated adaptive/basin-width-aware search
  (a genuinely different, more expensive task if greenlit — do NOT build it in this pass, just
  recommend go/no-go and rough scope if yes), or is the symmetric-class census sufficient to
  close out this thread?
  **Explicitly out of scope for this pass:** running the V1-V4-strict gauntlet itself, writing
  to `data/catalogue.yaml`, building the asymmetric-closure search. This adjudication produces
  a prioritized plan and recommendations; it does not execute them.
  **Recommended models:** Opus first (trust-bearing discovery-verdict judgment, matching #561's
  own precedent for this exact territory), then **Fable** second-opinion given the stakes — per
  #561's own two-model discipline, and given the demonstrated value of that discipline THIS
  session (Fable caught a real, load-bearing flaw in #561's own first-pass criterion that
  changed the outcome materially).
  **RESULT (2026-07-11, Opus adjudication pass) — full verdict in
  `docs/notes/2026-07-11-564-opus-adjudication-563-family.md`; Fable second-opinion still
  owed per the two-model discipline before any Stage-3 dispatch or writeback.** Read all 60
  raw passes → 30 deduped physical loops directly; confirmed all 30 are the exact
  "perpendicular-crossing" symmetric periodic-orbit signature (machine-precision by
  construction), so residual is now a dead triage axis — the live axes are max-bend/geometric
  quality, V∞, tof, and Titania-Oberon literature risk.
  **(1) Gauntlet subset = 5, not 30.** The census claim ("family across the non-Miranda
  pairs") is supported by ONE V4-strict-validated real-ephemeris member per pair, not by
  re-validating exact-by-construction siblings within an already-represented pair. #312
  (Umbriel-Oberon) is already validated; recommend gauntleting one shortest-tof / cleanest
  representative in each of the other 5 pairs — Titania-Oberon n=1 (0,0) r180 (tof 12.316d,
  bend 7.03°, MANDATORY: discharges the Canales/Kumar test), Ariel-Umbriel n=1 (0,0) (tof
  3.216d, shortest of all 30), Ariel-Titania n=3 (0,0), Ariel-Oberon n=5 (0,0), Umbriel-Titania
  n=1 (0,0). Tier-2 optional (Umbriel-Titania n=4 (0,0), bend 5.59° = cleanest of the family;
  the #312 n=7 sibling). Explicitly do NOT gauntlet the extreme-bend low-V∞ Titania-Oberon rows
  (bend 47°/88°, V∞ 0.62/0.34) — near-surface flybys, least robust, add nothing to the claim.
  **(2) Reframe #312 — YES, with a guardrail.** Cast it as first-documented member of a
  30-member SYMMETRIC-closure family, but keep two counts strictly distinct: 30 =
  idealized-model-exact enumerated closures (symmetric class, within #558's searched tof range);
  1 (→up to 6) = real-ephemeris gauntlet-validated. Writing "30 validated quasi-cyclers" would
  overclaim (29 unproven under DE440/URA111; #312's own kernel-edge V4 failures show the
  idealized→real gap is not automatic). Draft prose + qualifier in the note; apply in the
  separate writeback step (NOT done here), ideally after ≥2-3 pairs are validated.
  **(3) Asymmetric near-closures — DEFER (NO-GO).** The symmetric census answers the original
  "isolated-or-family?" question definitively; the 2 residual asymmetric near-closures are both
  in the highest-lit-risk Titania-Oberon basin, would need a fundamentally more expensive 2D
  adaptive search with NO exhaustiveness guarantee (unlike the provably-finite symmetric class),
  and one sits near the 0.05 gate floor. If ever greenlit: standalone P2/P3 task, T-O-only,
  seeded at the 2 known near-closures, gated behind resolving Canales/Kumar novelty first.
  Lit grounding verified against `docs/notes/2026-06-16-328-uranian-cycler-lit-deep-dive.md`
  + the Canales-Howell-2023 and Kumar-2025 digests: both structurally adjacent, NEITHER a direct
  topology match (Canales = halo/manifold MMAT transfer; Kumar = single-moon MMR heteroclinics),
  so the mandatory-citation obligation holds without a novelty collision.
  **FABLE SECOND OPINION (2026-07-11): REJECT** — one load-bearing flaw, one wrong citation
  ID, one internal contradiction. See task **#565** below for the correction.

- **#565** (P0, judgment-only, blocks Stage-3 dispatch — supersedes #564's triage §2/§3 pending
  correction) — Fable-corrected re-adjudication of #564. Fable's adversarial review (full text
  in the #564 result above; verified independently before filing this task — read
  `src/cyclerfinder/search/physical_sanity.py:85-114` directly, do not take either adjudication's
  word for it) found #564 inverted the semantics of `max_bend_deg_per_encounter`: the field is
  the **maximum achievable ballistic deflection** at the safe-altitude periapsis (a pure
  function of body + V∞, gated as a **floor** — pass iff ≥5°; see `flyby_is_useful()` and the
  #324 gate). Higher = more capability margin; lower = closer to the reject floor. **#564 read
  it backwards** ("high bend = demanded/fragile turn, low bend = clean/comfortable margin") and
  every ranking built on that axis is inverted:
  (a) row 12 was called "cleanest closure, best-geometry showcase" but is actually the most
  MARGINAL member of all 30 (min per-encounter achievable bend 5.39°, barely above the 5.0°
  floor);
  (b) rows 29-30 were excluded as "demanding near-surface/impossible periapses" but in fact have
  the LARGEST capability margins in the family (min per-encounter 17.77° and 45.29°) — nothing
  in the pipeline computes a *required* turn angle anywhere (`residual_at_point` matches V∞
  magnitude only), so "demands an impossible periapsis" is not a claim the data supports;
  (c) the Titania-Oberon P1 pick (row 24, "cleanest bend of the 8 T-O closures" at 7.03°) is by
  the correct margin statistic (min per-encounter achievable bend, matching the statistic #563's
  own NEW-closures table already used) nearly the WORST of the 8 T-O closures (min 6.27°; only
  row 25 at 6.23° is lower) — row 23 (min 8.85°) or row 27 (V∞ 0.827/0.959, min bend 24.94°,
  essentially #312's own validated V∞/bend profile transplanted to Titania-Oberon) are better
  candidates to "carry" the mandatory literature test, at the cost of a longer tof for row 27.
  Also found: **Kumar 2025's arXiv ID is wrong in #564** — cited as arXiv:2509.12675 (which is
  Kumar-Rawat-Rosengren-Ross, Earth-Moon cislunar resonant transport, zero Uranian content); the
  actual Uranus-Oberon MMR paper is **arXiv:2509.03655** §6.2 per
  `docs/notes/2026-06-16-328-uranian-cycler-lit-deep-dive.md:64,86-102` — exactly the
  same-author/same-year citation-collision trap flagged in
  [[feedback_ground_citations_against_content]]; the substantive adjacency verdict (single-moon
  MMR, not two-moon cycler topology) is correct, only the ID is wrong, but a wrong ID must not
  propagate into any writeback. Also: an internal contradiction (row 22 both recommended as
  Tier-2 "cheap to confirm" and excluded "at any tier" under a nonexistent "row 47" label — one
  statement must go). Fable's independent re-verification of everything else in #564 (all 30
  table rows against the raw jsonl, the 14/5/11 #562 cross-reference, #312's own identity/status,
  both asymmetric near-closure numbers, the "5 not 30" cutoff LOGIC itself, and the reframing
  prose) came back clean — this is a correction, not a full teardown.
  **Scope:** re-derive the geometric-quality ranking in #564 §2 using the CORRECT semantics (min
  per-encounter achievable bend = capability margin above the 5° floor, not a "demand"); re-pick
  the Titania-Oberon gauntlet representative among rows 23/24/27/28 under the corrected axis
  (rows 1/5/7/10 for the other 4 pairs are unaffected — they were justified on tof/order, not
  bend, and stand); withdraw or re-justify the row-12 "showcase" call and the rows-29/30
  exclusion (their real, valid exclusion reason is same-pair redundancy with an already-selected
  representative, not "impossible periapsis" — restate it that way if kept excluded); fix the
  Kumar citation to arXiv:2509.03655; resolve the row-22 contradiction. Do NOT re-open #564's §1
  (data verification), the "5 not 30" cutoff count, the #312 reframing recommendation (§5,
  Fable confirmed this is well-guarded as written), or the asymmetric-closures deferral (§3,
  Fable confirmed this stands) — those parts of #564 are not in question, only the bend-axis
  ranking and the two citation/contradiction nits.
  **Explicitly out of scope:** running the V1-V4-strict gauntlet, editing `data/catalogue.yaml`.
  **Recommended model:** Opus (same trust-bearing territory as #564; this is correcting its own
  prior output against Fable's findings, not a mechanical fix).
  **RESULT (2026-07-11, Opus corrected re-adjudication) — full corrected verdict in
  `docs/notes/2026-07-11-565-fable-corrected-adjudication.md`; supersedes #564 §2/§3's bend-axis
  ranking and picks, leaves #564 §1/§5/§6 + the 5-not-30 cutoff standing.** Confirmed the semantic
  inversion directly against `src/cyclerfinder/search/physical_sanity.py:80,100-103`
  (`flyby_is_useful()` passes iff `max_bend_deg >= 5.0°`; higher = more margin) and recomputed the
  correct **min-per-encounter achievable bend** for all 30 deduped loops from the raw jsonl.
  **(1) Bend axis re-derived (corrects §2).** Marginal members (near the 5° floor) are rows 16/9/12/13
  (min ~5.1–5.5°); robust members are rows 30/27/29/15/19/22 (~15–45°). Decisively, **#312 itself
  (row 19, min 15.55°) sits at the ROBUST end**, and the two rows #564 excluded as "most
  fragile/impossible" (29 at 17.77°, 30 at 45.29°) have LARGER margins than #312 — the inverted axis
  had reversed the whole robustness ordering.
  **(2) Gauntlet set — count/logic unchanged (5, not 30), one pick corrected.** New set: **rows 23, 1,
  5, 7, 10** (was 24/1/5/7/10). **Titania-Oberon: row 24 → row 23** — row 24 was actually the
  2nd-WORST T-O margin (min 6.27°), whereas row 23 has the best margin of the three shortest-tof T-O
  rows (min 8.85°, moderate V∞ 1.61/1.80, NEW), keeping the tof/V4 advantage. **Row 27** (min 24.94°,
  V∞ 0.827/0.959 ≈ #312's own profile transplanted to T-O) named as the high-margin alternative at 2×
  tof. Rows 1/5/7/10 picks STAND (tof/order-dominated), but rows 5 and 7 had inverted bend sub-claims
  in their stated reasoning, now restated on tof/order grounds.
  **(3) Row 12 "best-geometry showcase" WITHDRAWN** — it is the 3rd-most-MARGINAL member (min 5.39°),
  not the cleanest. **Rows 29/30 exclusion RESTATED** as same-pair redundancy with the T-O
  representative (row 23), NOT "impossible periapsis" (no required turn angle is computed anywhere).
  **(4) Kumar citation fixed** to **arXiv:2509.03655 §6.2** (was the wrong 2509.12675, an unrelated
  Earth-Moon paper); substantive adjacency verdict unchanged. **(5) Row-22 contradiction RESOLVED** —
  excluded from the primary-5 on same-pair redundancy (U-O carried by #312), robust not fragile (min
  15.43°); the nonexistent "row 47 / bend 42°" near-surface language dropped. Everything Fable
  confirmed in #564 (§1 data, 5-not-30 cutoff, §5 reframing, §6 asymmetric deferral) left standing.
  **FABLE CONFIRM (2026-07-11): CONFIRMED — safe to proceed to Stage-3 gauntlet dispatch on the
  corrected 5-representative list (rows 23, 1, 5, 7, 10).** Independently recomputed all 30
  min-per-encounter bends from the raw jsonl (not a sample) — every #565 §1 value and the ranked
  ordering matched to 2 decimals. Verified the row-23 T-O pick against rows 24/25/27 directly:
  row 23 genuinely dominates the same-tof alternatives on margin and V∞; row 27 correctly
  characterized as the margin-vs-tof trade alternative. Citation, row-22 resolution, and
  supersession scoping (§1/§5/§6/5-not-30-cutoff left undisturbed) all verified clean. One
  non-blocking prose nit noted (a loosely-ordered list in §1's prose vs. the correctly-sorted
  table) — no action needed. **This closes the #564→#565 adjudication loop.** Stage-3
  (V1-V4-strict gauntlet on rows 23/1/5/7/10) and any catalogue.yaml writeback (the #312
  reframing per #564 §5) remain explicit separate asks — not auto-fired by this adjudication.

- **#566** (P0, mechanical reuse — cheap, seconds-scale, user-requested "run the gauntlet") —
  execute the V2→V3→V4→V4-strict gauntlet (the same chain #312 itself passed via #330/#331/
  #332/#335/#338) on the 5 candidate closures #565 recommended: rows **23** (Titania-Oberon,
  mandatory literature representative), **1** (Ariel-Umbriel), **5** (Ariel-Titania), **7**
  (Ariel-Oberon), **10** (Umbriel-Titania) from `data/enumerate_563_symmetric_closures.jsonl`.
  **No new capability needed** — the underlying drivers (`run_v2_moontour`, `run_v3_3d`,
  `run_v4_uranus`, `run_v4_uranus_strict` — used by #330-#338) are fully generic, taking
  `sequence`/`vinf_tuple_kms`/`leg_tofs_days`/`rel_offset_deg`/`n_revs`/`phase0_deg` as plain
  arguments; only the *wrapper scripts* (`run_330`/`run_331`/`run_332`/`run_335`/`run_338`) are
  hardcoded to #312's own Umbriel-Oberon numbers. `phase0_deg` is provably rotation-redundant
  (#558) — reuse the same fixed constant (29.999999999999996) already used across #330/#335/
  #338. The URA111 SPICE kernel needed for V4-strict is already vendored locally (confirmed at
  `~/GMAT/R2022a/data/planetary_ephem/spk/uranian/ura111.bsp`) — no network fetch required.
  **Exact 5 candidate parameter sets** (sequence=(anchor,flyby,anchor), n_rev symmetric (n,n)):
  - Row 23: seq=(Titania,Oberon,Titania), n=1, n_rev=(0,0), rel_offset_deg=180.0,
    tof_days=(12.316046445872583, 12.316046445872583),
    vinf_kms=(2.161767816675378, 1.9680495724521725, 2.1617678166753755)
  - Row 1: seq=(Ariel,Umbriel,Ariel), n=1, n_rev=(0,0), rel_offset_deg=0.0,
    tof_days=(3.216088179066208, 3.216088179066208),
    vinf_kms=(0.979040480994661, 1.3004234339628207, 0.9790404809946573)
  - Row 5: seq=(Ariel,Titania,Ariel), n=3, n_rev=(0,0), rel_offset_deg=0.0,
    tof_days=(5.320895317317783, 5.320895317317783),
    vinf_kms=(1.2306411593828481, 1.7185773183747601, 1.2306411593828457)
  - Row 7: seq=(Ariel,Oberon,Ariel), n=5, n_rev=(0,0), rel_offset_deg=0.0,
    tof_days=(7.751820498940574, 7.751820498940574),
    vinf_kms=(1.520866047614147, 1.8285940380726622, 1.5208660476141462)
  - Row 10: seq=(Umbriel,Titania,Umbriel), n=1, n_rev=(0,0), rel_offset_deg=0.0,
    tof_days=(3.9544738760575804, 3.9544738760575804),
    vinf_kms=(1.2295656768416439, 1.0058255988095806, 1.229565676841644)
  (Verify these against the jsonl directly before running — do not trust this transcription
  blind.) **Scope:** write ONE new parameterized runner script (following the naming convention,
  e.g. `scripts/run_566_gauntlet_five_representatives.py`) that loops the 5 candidates through
  V2 (bounded-drift, closure ≤0.05 km/s / drift ≤50,000 km floors) → V3 (REBOUND IAS15
  cross-check, ~100 km agreement floor) → V4 (J2 + Uranian-moon n-body scipy DOP853) →
  V4-strict (real URA111 SPICE ephemeris; pick ONE representative launch epoch per candidate,
  not #338's full annual sweep — that's a separate, later robustness question if any candidate
  passes V4) — reusing the generic drivers directly, not copy-pasting #330-#338's bodies.
  `preflight_search()` is MANDATORY (two past CI breaks this session from missing it). Record
  full per-stage pass/fail + numeric detail to a `data/*.jsonl` results file and a result
  paragraph under this bullet — do NOT write to `data/catalogue.yaml` in this pass (that remains
  a deliberate, separate writeback step even for candidates that pass everything).
  **Explicitly out of scope:** the annual DOY-sweep robustness pass (only warranted for
  candidates that already clear single-epoch V4-strict — apply the #559/#560 lessons: check
  Lambert-branch continuity + degenerate-arc handling if any candidate's V4-strict looks
  singleton-flip-suspicious, per [[feedback_isolated_sweep_flips_suspect_artifact]]); catalogue
  writeback; any new capability building (#556/#560 stay parked, unrelated to this task).
  **Recommended model:** Sonnet (mechanical reuse of already-validated generic drivers behind
  the gauntlet's own deterministic pass/fail floors — no new judgment calls, matches this
  project's tiering policy for spec-complete work behind a strong gate). Follow with an Opus (or
  self) pass/fail adjudication once results land — do not let a raw pass/fail table stand as the
  final word without a short judgment pass, matching #561-#565's own discipline for this thread.
  **✓ RUN (2026-07-11, `scripts/run_566_gauntlet_five_representatives.py` →
  `data/gauntlet_566_five_representatives.jsonl`, 11.9s total, well inside budget).** Numeric
  transcription check first: the 5 candidate tuples quoted in this bullet were verified
  field-for-field against `data/enumerate_563_symmetric_closures.jsonl` before running anything —
  every sequence/tof_days/vinf/rel_offset_deg/n_rev/n_commensurate_int matches exactly, but at
  jsonl lines **57, 2, 12, 18, 26** respectively (NOT the "Row N" labels this bullet uses — those
  labels do not correspond to jsonl line numbers; the underlying numbers are correct regardless).
  **Result: all 5 candidates ran the identical PASS pattern to #312's own SILVER (#327).** Every
  candidate FAILs the strict V2 50,000 km drift floor (drift 2.9e5–8.7e5 km across n_cycles
  3/5/10) but with a per-cycle closure residual ~1e-14 km/s — i.e. `FAIL_QUASI_BOUNDED` under
  #330's own three-way test (the same test that admitted #312 as `quasi_cycler`, not the strict
  V2 PASS path). Proceeding past that on the same basis #330 did: **V3 (REBOUND IAS15) PASSes for
  all 5 at all 3 n_cycles (agreement ~1e-8–1e-7 km, far inside the 100 km floor); V4 (J2 +
  Uranian-moon n-body) PASSes for all 5 at all 3 n_cycles (agreement ~2e3–7e3 km, inside the
  50,000 km floor, bounded_drift_survives=True throughout); V4-strict (real URA111 SPICE, single
  representative epoch 2000-06-21T00:00:00, reusing the #338/#559-established reference epoch)
  PASSes for all 5 at all 3 n_cycles (agreement-vs-V3 ~2.3e3–6.0e3 km, agreement-vs-V4-scipy
  ~1.5e3–1.9e4 km, both inside the 50,000 km floor).** Row 23 / jsonl-line-57
  **Titania-Oberon-Titania — the MANDATORY literature-clearance representative — PASSes the full
  chain** (`chain_verdict=PASS_AS_QUASI_CYCLER`), same as the other 4 (Ariel-Umbriel,
  Ariel-Titania, Ariel-Oberon, Umbriel-Titania). **Caveat carried into the results (not a
  capability blocker):** `run_v4_uranus_strict`'s 4 audit-only eccentricity/inclination fields are
  hardcoded to always SPICE-sample Umbriel/Oberon regardless of the candidate's actual sequence;
  the physics propagation itself (`_cycle_v4_strict`) is fully generic over `sequence` (uses
  `_moon_state_spice(moon, ...)` per actual moon name), so `passes_v4_strict` is unaffected — only
  those 4 descriptive fields are mislabeled for the 4 non-Umbriel-Oberon candidates; noted in the
  jsonl's `_meta.known_caveat` field for anyone reading the raw eccentricity numbers later.
  **This is a striking, uniform result** — 6 of the 30 exact #563 closures (all 6 non-Miranda
  pair directions now have at least one member gauntleted) clear the identical computational
  chain #312 cleared. No catalogue writeback performed (out of scope, per this bullet's own
  scope note) — the natural next step is an Opus/self adjudication pass on whether this
  uniformity argues for admitting the other 4 as their own `quasi_cycler` catalogue rows (a
  family, not an isolated #312) or whether it argues the whole V2→V4-strict chain is measuring
  something structural to ALL 30 exact closures rather than discriminating between them (in
  which case the gauntlet's real discriminating power over this candidate set is in question) —
  neither conclusion should be assumed without that follow-up pass.
  **✓ ADJUDICATED (2026-07-11, Opus, judgment-only — no new computation; full note
  `docs/notes/2026-07-11-566-opus-adjudication-gauntlet-results.md`).** Read all 63 jsonl records
  per-stage and `v4_uranus_strict.py` in full. **Verdict: the 5/5 PASS is GENUINE (not a
  malfunction) AND the gauntlet is non-discriminating for this class — for a benign reason, and
  both are true at once.** The V2 closure ~1e-14 is tautological (exact closure is a #563
  construction property), V2 drift overflow is the class definition, V3 agreement ~1e-8 km is a
  numerical self-consistency check, and V4/V4-strict agreement (2–7e3 km vs a 50,000 km floor)
  confirms only that J2 + real ephemeris don't materially bend the SHORT 3–25 d Lambert arcs — a
  physically benign regime where every member passes. The negative baseline settles the read: a
  REAL rejection in this stack looks like #535's Earth corridor (collapsed under ER3BP
  eccentricity), whereas #559 showed #312's own ≈10–14% daily V4-strict FAILs are CONFIRMED
  NUMERICAL ARTIFACTS (Lambert branch-flip + DOP853 planet-crossing), not real dynamics — so the
  chain has essentially zero demonstrated power to physically reject a real member of this class.
  A PASS is therefore trustworthy as *class-membership* confirmation, NOT as independent
  per-candidate existence evidence. **Audit-field bug: confirmed cosmetic for the verdict**
  (lines 575–578 hardcode Umbriel/Oberon; `_cycle_v4_strict` + `passes_v4_strict` never read those
  4 scalars — every candidate records identical e_u=0.0041/e_o=0.0016 incl. Ariel-Titania-Ariel),
  **but a provenance defect at writeback** (freezing them into a catalogue row bakes a wrong e/i
  into 4 of 5 rows — [[feedback_digest_not_adoption]]) — must be fixed before any writeback.
  **Single-epoch V4-strict (2000-06-21) is INSUFFICIENT confidence**: that epoch is #312's own
  known-favorable #338 anchor, #312 is known to vary across epochs (#338 found 2000-01-15 FAILs),
  the 5 have NO epoch-sensitivity data, and a scan run under today's code would reproduce #559's
  confirmed artifacts (the #560 Lambert-continuity + perijove-guard fixes are the prerequisite for
  a meaningful scan). **Recommendation: HOLD catalogue writeback** — do NOT admit any of the 5 (nor
  upgrade #312) on #566 alone; gate on **#567** (apply #560 fixes + audit-field fix + per-
  representative epoch-robustness scan). When #567 clears, admit #312 + the epoch-robust survivors
  as the first-documented members of the 30-member symmetric-closure family per #564 §5's
  reframing (validated-N kept distinct from enumerated-30), T-O representative citing
  arXiv:2110.03683/2308.10029 + Kumar arXiv:2509.03655 §6.2 (the #565 §5-corrected ID).

- **#567** (P0, catalogue-writeback prerequisite — the gate the #566 adjudication concluded must
  clear before any of the 5 representatives or #312 are written/upgraded) — build a *trustworthy*
  epoch-robustness evidence base for the #566 5-representative list (rows 23/1/5/7/10 =
  Titania-Oberon / Ariel-Umbriel / Ariel-Titania / Ariel-Oberon / Umbriel-Titania) and #312.
  **Why this is needed and not optional** (see `docs/notes/2026-07-11-566-opus-adjudication-
  gauntlet-results.md`): #566 tested all 5 at a SINGLE epoch (2000-06-21), which is precisely
  #312's own known-favorable #338 anchor; #312 is known to vary across epochs (#338 found
  2000-01-15 FAILs, 91,000 km drift); and #559 proved that #312's ≈10–14% daily V4-strict FAILs
  are CONFIRMED numerical artifacts (Lambert branch-selection flip + DOP853 planet-crossing stiff-
  death), both living in `src/cyclerfinder/data/validation/v4_uranus_strict.py`. A naive per-
  candidate daily scan run under today's code would therefore reproduce those same artifacts and
  tell us nothing real ([[feedback_isolated_sweep_flips_suspect_artifact]]). **Scope, in order:**
  (1) **apply the already-fully-diagnosed #560 robustness fixes** to `v4_uranus_strict.py` —
  Lambert branch-selection continuity tracking (`_cycle_v4_strict` ~line 418; kill the
  no-continuity `min()` tie-break that flips 23,448→1,340 km between adjacent hours) + a
  perijove/collision guard in `_v4_strict_propagate_leg` (~lines 360-372) that distinguishes a
  genuine dynamical FAIL from a DOP853 stiff-death on a non-physical near-parabolic planet-crossing
  arc (perijove inside Uranus's R_eq 25,559 km). **PIN (Fable check on the #566 adjudication,
  2026-07-11): the guard must NOT neutrally skip/exclude perijove-inside-Uranus epochs.** #559's
  own entry documents the trigger as REAL geometry (synodic phase sweeps ~60°/day; transfer angle
  collapses to 3-12°; the fixed-TOF rev-1 arc then genuinely planet-crosses) — only the *silent
  stiff-death misclassification* is the artifact, not the underlying infeasibility. These epochs
  are genuine dynamical infeasibility and MUST count as FAIL against each candidate's PASS-band
  width (they define the real `validity_window` boundary per the quasi_cycler schema), not be
  excluded as noise — a guard that neutrally skips them would inflate every candidate's apparent
  band width and cause the writeback validation level in step (4) below to over-claim. Expect the
  honest post-fix band to come out NARROWER than the raw #559 86-90% pass-rate, not wider; that
  narrowing is a feature (real boundary), not a regression. (2) **fix the hardcoded Umbriel/Oberon
  audit-field
  sampling** (`v4_uranus_strict.py` lines 575-578) so recorded `eccentricity_used_* /
  inclination_used_*` track the candidate's ACTUAL sequence — cheap, and prevents freezing a wrong
  e/i into 4 of 5 catalogue rows at writeback ([[feedback_digest_not_adoption]]; the #566 jsonl
  shows all 5 recording identical e_u=0.0041/e_o=0.0016, wrong for the 4 non-U-O candidates).
  **PIN (Fable check): `tests/data/test_v4_uranus_strict.py:290-291`'s existing range-check
  (`0 < e < 0.02`) passes even with wrong-moon values and would NOT catch a regression of this
  fix — add a sequence-consistency assertion (recorded body matches the candidate's actual
  sequence) alongside the fix, not a bare range check.** (3)
  under the fixed code, run the #338-style annual (100-epoch, ~42s) + #559-style daily
  epoch-sensitivity scan on all 5 representatives AND re-confirm #312, reporting each candidate's
  true PASS-band WIDTH (wide-tolerant vs narrow-knife-edge) rather than a single point — and, since
  the #560 fixes will change #312's own committed #559 pass-rate/spike numbers, refresh
  `tests/data/test_sweep_diagnostics.py::test_559_actual_pattern_reproduces_known_singleton_count`
  against the post-fix data (it is a regression pin on the OLD artifact-contaminated counts); (4)
  produce the writeback-readiness verdict per candidate (band width + whether it clears at #312's
  V4 level or a capped level), feeding the actual `data/catalogue.yaml` admission decision — which
  remains a SEPARATE writeback step, NOT part of #567. **A bug-fix invalidates prior negatives
  discipline** ([[feedback_bugfix_invalidates_past_searches]]): once #560 lands, #559's 89.9%/85.8%
  daily pass-rates and the 55/58 spike locations are superseded and must be recomputed, not reused.
  **Explicitly out of scope:** the catalogue writeback itself (separate step post-#567); any change
  to the V2/V3/V4 driver logic (only V4-strict's two artifact generators + the audit field);
  extending to the other 25 same-pair-redundant closures (#565 §2/§4 — census-minimum is one
  representative per pair). **Recommended model:** the #560 fixes + audit-field fix are TDD-shaped
  numerical-methods work behind concrete diagnosed failure cases (Sonnet, behind a regression pin);
  the epoch-band-width interpretation + per-candidate writeback-readiness verdict is trust-bearing
  (Opus or self, once the scan data lands) — do not let a raw post-fix pass-rate table stand as the
  writeback green-light without a judgment pass, matching this thread's #561–#566 discipline.
  **FABLE CHECK on the #566 adjudication (2026-07-11): CONFIRMED WITH CORRECTIONS — HOLD +
  #567 gating is sound.** Independently spot-checked #566's numeric claims against the raw jsonl
  (all matched), confirmed the audit-field hardcoding claim against `v4_uranus_strict.py` source
  directly, and confirmed the #535/#559 negative-baseline characterizations against their own
  OUTSTANDING.md entries. On the crux question — does "weak-to-zero per-candidate rejection
  power" mean NO evidence from this chain could ever validate anything in this class — verdict is
  NO: the chain's rejection power is real (demonstrated independently by #535's corridor
  collapse), just saturated within this benign short-arc population; #567's epoch-band-width axis
  is a genuinely different, discriminating test (not "more of the same"), so HOLD-pending-#567 is
  coherent, not circular. The two corrections above (perijove-guard pin, audit-field test-guard
  pin) are now folded into #567's scope directly; a third, cosmetic-only nit (#566's note said
  "63 records", jsonl actually has 62 lines) does not affect any conclusion and is not corrected
  here (informational only, no OUTSTANDING.md claim depends on the exact line count).
  **RESULT — steps (1)+(2) DONE (commit `6c54bba`, 2026-07-11).** All three bugs in
  `src/cyclerfinder/data/validation/v4_uranus_strict.py` fixed with TDD regression pins against
  directly-instrumented real failure instances (not synthetic): (1) Lambert branch selection now
  picks by actual propagated terminal offset — a continuous quantity — instead of a
  departure-velocity-match proxy, eliminating the discontinuous branch-flip (pinned: a located
  2000-09-06 T03→T04 hour pair that jumped ~25,000 km→~3,500 km pre-fix now varies smoothly,
  <5,000 km hour-to-hour); (2) planet-crossing legs are now pre-screened and tagged
  `FAILURE_MODE_PLANET_CROSSING` with periapsis recorded, per the Fable pin — still counted as a
  real FAIL, never excluded (pinned: a located 2000-07-24T02:00 epoch where both Oberon→Umbriel
  Lambert branches cross Uranus at 97/852 km periapsis, confirmed via direct probe); (3) the 4
  audit-only e/i fields now track the candidate's own first two distinct `sequence` bodies
  instead of hardcoded Umbriel/Oberon, with a sequence-consistency regression test (a bare range
  check would not have caught this) pinned against the #566 Ariel-Titania-Ariel representative.
  Downstream scripts `run_335`/`run_338`(×2)/`run_559` updated for the renamed
  `V4UranusStrictVerdict` fields (JSONL output keys unchanged — all three are Umbriel-Oberon-only
  scripts). Full `tests/data tests/search` ratchet suite passes clean (0 FAILED, 0 tracebacks).
  **Note on how this stage itself ran**: the dispatched subagent implemented and tested
  everything correctly, but then launched the full ratchet suite via backgrounding and ended its
  turn "waiting" — the exact fatal pattern [[feedback_subagent_background_is_fatal]] warns
  against, hit despite an explicit "don't background" instruction in the dispatch prompt (the
  instruction itself was flawed: it hedged "blocking unless it looks like >8 min," but this
  suite's runtime isn't reliably under the 10-minute hard cap, so there was no legal foreground
  path once it ran long). The actual OS-level pytest process kept running independently of the
  dead subagent turn; the top-level session (me) monitored it directly to completion and did the
  code review + commit myself. Memory sharpened accordingly: the full ratchet suite (and any
  other variable-duration command) must never be handed to a subagent, even conditionally — see
  the memory file for the updated rule. **Steps (3)+(4) (the actual multi-candidate
  epoch-robustness scan + writeback-readiness verdict) remain open**, to be executed with the
  build+validate(subagent)/launch+own(me)/analyze(fresh agent) split from the start this time.
  **RESULT — step (3) DONE (raw scan data committed; commit `a777c2d` for the script, scan
  output at `data/scan_567_epoch_robustness.jsonl`).** `scripts/run_567_epoch_robustness_scan.py`
  was built (subagent, smoke-tested only) then launched and owned directly by the top-level
  session (this time correctly split: no subagent touched the long-running command) — 4,986
  epoch-cells across 6 candidates × 3 sweep types (annual 2000-2099 100pt, daily-DOY 2000 366pt,
  daily-DOY 2030 365pt), 651.6s parallel_sweep wall-time, 0 crashes.
  **Strong positive confirmation the #567(1)+(2) fixes worked**: `lambert_no_solution=0` and
  `integrator_failure=0` in ALL 18 groups — every single FAIL across the whole scan is now
  cleanly tagged `planet_crossing_infeasible`, meaning the two known artifact generators are
  genuinely gone, not just relabeled.
  **Unexpected finding, NOT yet adjudicated**: `detect_isolated_singleton_anomalies()` still
  flags isolated flips post-fix, at wildly candidate-dependent rates — 0% (Titania-Oberon, both
  daily sweeps; 23% its own annual) up to 44% (Ariel-Umbriel annual; 32-33% its daily sweeps).
  Full per-group table:
  | candidate | annual flip% | daily-2000 flip% | daily-2030 flip% | annual PASS% |
  |---|---|---|---|---|
  | #312 (Umbriel-Oberon) | 2.0% | 10.9% | 9.6% | 88.0% |
  | Titania-Oberon (row 23, mandatory) | 23.0% | 0% | 0% | 74.0% |
  | Ariel-Umbriel (row 1) | 44.0% | 32.0% | 32.9% | 63.0% |
  | Ariel-Titania (row 5) | 14.0% | 22.7% | 22.5% | 64.0% |
  | Ariel-Oberon (row 7) | 19.0% | 28.7% | 29.0% | 71.0% |
  | Umbriel-Titania (row 10) | 9.0% | 4.1% | 3.3% | 70.0% |
  This is NOT simply "the artifact is back" (the two diagnosed bugs are confirmed gone) and NOT
  simply "trust it, it's real physics" either — the guard's own docstring says "verify before
  trusting," and this project's own discipline
  ([[feedback_isolated_sweep_flips_suspect_artifact]]) requires diagnosing WHY before either
  believing or discarding a flip population. Two live hypotheses, undiagnosed: (a) a THIRD,
  not-yet-found artifact generator in the same code family (the per-candidate variability
  tracking so closely with orbital/synodic period differences between moon pairs is suspicious
  either way — could mean a genuine sampling-aliasing signature of a real fast dynamical
  structure the daily/annual grid under-resolves, OR a numerical sensitivity that scales with
  the same period differences); (b) genuine physical structure (the actual PASS/FAIL boundary
  in real synodic geometry oscillates fast enough that even a DAILY grid aliases it — which
  would still be legitimate for reporting a validity_window, just needs characterizing
  correctly, not discarding). **Step (4) (writeback-readiness verdict) explicitly BLOCKED on
  resolving this** — do not treat any candidate's raw PASS% above as a trustworthy
  `validity_window` figure until this is diagnosed, per the same discipline that blocked #566's
  premature single-epoch read.
  **RESULT — step (3) singleton-flip anomaly DIAGNOSED (2026-07-11), verdict = (b) genuine
  synodic aliasing, NOT a third artifact** (full writeup:
  `docs/notes/2026-07-11-567-post-fix-singleton-flip-diagnosis.md`). The post-fix isolated flips
  are genuine physical aliasing of a real, fast-oscillating `planet_crossing_infeasible` boundary
  (the transfer's osculating periapsis crossing Uranus's `r_eq` 25,559 km) whose frequency IS each
  moon-pair's SYNODIC frequency — an artifact of the DIAGNOSTIC GRID, not the solver. Four
  independent confirmations: (1) the pass/fail bit is a step function of a continuous geometric
  quantity (periapsis) crossing a fixed threshold; (2) the FFT dominant period of every daily
  PASS/FAIL series matches that pair's synodic period to <2% (Titania-Oberon 24.6 d, Umbriel-Titania
  7.9 d, #312 6.0 d, Ariel-Titania 3.55 d, Ariel-Oberon 3.10 d; Ariel-Umbriel locks to HALF its
  6.43 d synodic period because the symmetric `rel_offset=0` there-and-back crosses `r_eq` twice per
  synodic cycle — hence its shortest effective period 3.2 d and highest flip rate), and the
  singleton count is monotonic in boundary period (24.6 d→0 flips, 3.2 d→120 flips) — the textbook
  aliasing signature, spectrally locked to astronomical synodic frequencies a solver artifact could
  not produce; (3) a direct 4-hour sub-daily probe of the Ariel-Umbriel 2000-01-04 daily singleton
  resolves it into a CONTIGUOUS ~1.3-day FAIL band with the offending branch's periapsis varying
  smoothly/continuously across it (17,197→19→…→24,619 km, approaching `r_eq` at the exit edge) and
  smooth low drift on both PASS sides — the definitive "zoom finer" discriminator (a solver artifact
  stays isolated at any resolution; a real under-sampled boundary resolves into a band). This also
  clears the step-(5) concern: the #567 branch-selection fix introduces no discontinuity here
  (periapsis is continuous through the band); (4) flip rate tracks boundary FREQUENCY, not PASS%
  proximity to 50% (Titania-Oberon 0 daily flips at PASS 74%; Ariel-Umbriel most flips at PASS 61%).
  A small #312-ONLY secondary population (~1-4%) fails the DIFFERENT mode (i) — all legs converge but
  V4-vs-V3 drift exceeds the 50,000 km floor, with values up to 103,413 km (far from knife-edge) —
  also genuine physics (real high-perturbation epochs), a threshold on a different continuous metric.
  **Implication for step (4):** the raw PASS% at any fixed grid is grid-resolution-dependent and is
  NOT a wall-clock `validity_window`; these quasi_cyclers have a feasible synodic PHASE band that
  recurs every synodic period. Characterize each candidate by (i) its feasible DUTY CYCLE (≈ the
  daily PASS% since daily samples are ~phase-uniform: Ariel-Umbriel ~61%, #312 ~78%) and (ii) its
  synodic boundary period — NOT by raw pass% + flip%, and the `detect_isolated_singleton_anomalies`
  flip% must be DEMOTED from an artifact-suspicion flag to an expected aliasing diagnostic for these
  candidates (reporting it as instability would wrongly discount legitimate duty cycles). **No new
  bugfix warranted** (re-characterization, not a defect). Step (4) is UNBLOCKED and re-homed as new
  task **#568** (duty-cycle + synodic-period characterization); catalogue writeback remains a
  separate step after that.

- **#568** (P0, catalogue-writeback prerequisite — the corrected replacement for #567 step (4),
  unblocked by the #567 step-(3) diagnosis above) — produce the per-candidate writeback-readiness
  verdict for #312 + the #566 5-representative list, framed as the #567 diagnosis requires: report
  each candidate's **feasible synodic DUTY CYCLE** (fraction of synodic phase for which the cycle is
  V4-strict-feasible — estimable directly from the already-committed `data/scan_567_epoch_robustness.jsonl`
  daily PASS% since daily epochs sample synodic phase ~uniformly, or refine with a phase-uniform /
  sub-synodic-resolution grid if a candidate's duty cycle sits near a writeback tier boundary) and
  its **synodic boundary period** (per the diagnosis table), NOT a raw pass% + singleton-flip% pair.
  The `detect_isolated_singleton_anomalies` flip% is CONFIRMED expected physical aliasing for these
  Uranian candidates and must NOT be used as an instability/artifact discount against a duty cycle.
  Deliver: per-candidate (duty cycle, synodic period, whether it clears at #312's V4 level or a
  capped level) feeding the actual `data/catalogue.yaml` admission decision. **Out of scope:** the
  catalogue writeback itself (separate step); any new V4-strict code change (the #567 diagnosis found
  NO third bug — do not "fix" the aliasing, it is real physics). **Recommended model:** trust-bearing
  interpretation → Opus or self. Guard [[feedback_isolated_sweep_flips_suspect_artifact]]:
  the flip population here is DIAGNOSED (genuine aliasing), so this task must not re-litigate it as
  suspect — it consumes the diagnosis, it does not repeat it.
  **CORRECTIONS (Fable adversarial check on the #567 diagnosis, 2026-07-11 — independently
  re-derived moon periods from `src/cyclerfinder/core/satellites.py`, re-ran per-window
  FFT/autocorrelation on all 6 candidates' daily series, and reproduced the diagnosis's
  previously-uncommitted sub-daily zoom-in probe from scratch, matching every number — overall
  verdict CONFIRMED, but two fixes are REQUIRED before #568 executes:**
  (1) **Do NOT use a "symmetric zero-offset ⇒ boundary period halves" RULE.** Only Ariel-Umbriel
  actually locks at synodic/2 (3.22d); Ariel-Titania, Ariel-Oberon, and Umbriel-Titania all have
  the same rel_offset=0° symmetric structure but lock at their FULL synodic period, not half —
  the diagnosis note's causal generalization for the halving does not actually predict this
  correctly. Use the MEASURED boundary period per candidate directly from
  `data/scan_567_epoch_robustness.jsonl` (already committed, already correct), never a derived
  "zero-offset ⇒ P/2" shortcut.
  (2) **Report the #312-only secondary drift-floor-exceeded population SEPARATELY, not folded
  into the duty-cycle number.** 11 of #312's 831 daily rows (~1.3%) are a DISTINCT failure
  mechanism from synodic planet-crossing aliasing — legs converge but V4-vs-V3 drift exceeds the
  50,000 km floor (50,084-103,413 km, confirmed genuine, not knife-edge). If #568 computes #312's
  duty cycle as raw daily PASS% (78%), it silently mixes this sporadic ~1.3% population into a
  number that's conceptually meant to be pure synodic-phase-band feasibility. Compute the duty
  cycle over `planet_crossing_infeasible`-mode FAILs only; report the drift-floor-exceeded rows
  as their own separate, small, distinct data point.
  (Fable also flagged, non-blocking: the diagnosis note's claim that spectral synodic-locking
  alone rules out a hidden solver artifact is overstated — the real discriminator is structural
  (`FAILURE_MODE_PLANET_CROSSING` requires EVERY Lambert branch's periapsis below r_eq, so a
  branch-selection discontinuity cannot produce it) plus the reproduced zoom-in band's smooth,
  deep (19 km at minimum vs 25,559 km threshold) periapsis excursion — not a jittery near-threshold
  crossing. The underlying conclusion (genuine physics, no third bug) is unaffected; this is a
  documentation-rigor nit in the diagnosis note, not something #568 needs to act on.)
  **RESULT (2026-07-11, Opus interpretation pass) — DONE; all 6 candidates writeback-READY at #312's
  V4-strict-equivalent (windowed) level. Full verdict:
  `docs/notes/2026-07-11-568-writeback-readiness-verdict.md`.** Feasible synodic DUTY CYCLE (daily-2000 +
  daily-2030, N=731, computed over `planet_crossing_infeasible`-mode fails only per Fable correction #2)
  and MEASURED synodic boundary period (my own FFT on the committed scan, per Fable correction #1 — no
  "zero-offset ⇒ P/2" shortcut; only Ariel-Umbriel actually halves, the other three rel_offset=0° reps
  lock at FULL synodic period): #312 Umbriel-Oberon **79.1 %** / 6.0 d; Titania-Oberon **74.4 %** / 24.4 d;
  Ariel-Oberon **71.0 %** / 3.1 d; Umbriel-Titania **68.5 %** / 7.9 d; Ariel-Titania **66.5 %** / 3.5 d;
  Ariel-Umbriel **61.7 %** / 3.2 d (=synodic/2). Every fail across all 5 representatives is
  `planet_crossing_infeasible` (zero drift-floor, zero solver artifacts), so their duty cycle = raw daily
  PASS%. **#312's secondary drift-floor population reported SEPARATELY (Fable correction #2):** 11/831 rows
  (~1.3 %; 10 of them in the 731 daily rows), all-legs-converge-but-V4-vs-V3-drift 50,084–103,413 km — a
  DISTINCT genuine-physics mechanism, NOT folded into #312's 79.1 %; none of the 5 reps show it. **Verdict:**
  a duty-cycle-limited window is the `quasi_cycler` class's own defining epoch-locked windowed validity
  (taxonomy note uses #312 itself as the worked example), NOT a defect — the lowest (Ariel-Umbriel 61.7 %)
  is still majority-feasible on a smooth boundary; report the window, do NOT cap. The 5 reps are cleaner
  than #312 on the drift axis (zero drift-floor fails). One writeback-MECHANICS qualifier (not a validation
  gap): each new row needs its frozen per-tier pytest-gate provenance generated to match #312's V4
  registration convention — mechanical writeback scope. **Recommends actual catalogue writeback proceed as
  new task #569 (allocated below):** write all 5 reps as validated quasi_cycler family members with
  per-candidate measured `validity_window` (duty cycle + synodic period), reframe the whole set as the
  30-member symmetric-closure family (#564 §5) with #312 as first-documented member, and UPDATE #312's own
  row (family reframing + its 79.1 %/6.0 d duty cycle + the ~1.3 % drift-floor caveat — its current
  {2000→2083} window misleadingly implies continuous validity).

- **#569** (P0, the actual catalogue writeback — the separate step every entry #566/#567/#568 held out of
  scope; the writeback-readiness gate is now CLEARED by the #568 verdict) — write the #312 Uranian
  symmetric-closure family into `data/catalogue.yaml`. **Scope, per the #568 verdict
  (`docs/notes/2026-07-11-568-writeback-readiness-verdict.md`):** (1) add 5 new `quasi_cycler` rows, one
  per non-Miranda moon pair (the #566 representatives: Titania-Oberon / Ariel-Oberon / Umbriel-Titania /
  Ariel-Titania / Ariel-Umbriel), each at #312-equivalent V4 (windowed), carrying its MEASURED
  `validity_window` = calendar {start,end} (URA111 kernel-coverage bounded) PLUS duty cycle + synodic
  boundary period from the #568 table (74.4 %/24.4 d, 71.0 %/3.1 d, 68.5 %/7.9 d, 66.5 %/3.5 d, 61.7 %/3.2 d
  respectively); (2) reframe the whole set (incl. #312) as the **30-member symmetric-closure family** (#563
  enumeration / #564 §5), with #312 the **first-documented member** — NOT a unique novel point; (3) UPDATE
  #312's existing `umbriel-oberon-1-1-uranian-quasi-cycler-2026` row: add the family relation, add its
  79.1 %/6.0 d duty cycle, and add its distinct ~1.3 % secondary drift-floor-exceeded population
  (50,084–103,413 km) to the `validity_window`/notes — correcting the current {2000-06-21→2083-06-21}
  window's misleading implication of *continuous* validity (it is a ~79 % synodic duty cycle recurring
  every 6.0 d); (4) generate the frozen per-tier pytest-gate provenance for each new row to match #312's V4
  registration convention (`tests/verify/test_silver_327_v*.py` pattern + `validate.py::_LEVEL_EVIDENCE`
  registration, spec §16.7.12) — this is the writeback-MECHANICS the #568 verdict flagged, not a validation
  gap; (5) re-confirm `literature_check.py` not-found for **Titania-Oberon** (the #565 mandatory
  literature-clearance obligation; structurally adjacent to Canales/Kumar, not a direct topology match, per
  the #562 finding) BEFORE its row lands. **Out of scope:** any V4-strict driver code change (the #567
  diagnosis found NO third bug — the synodic aliasing is real physics, do not "fix" it); the other 25
  same-pair-redundant closures (census-minimum is one rep per pair, #565 §2/§4). **Catalogue-edit
  discipline:** any `catalogue.yaml` row change ripples into MULTIPLE frozen-census ratchets
  ([[feedback_catalogue_edits_run_all_ratchets]]) — run the FULL `uv run pytest tests/data tests/search -q`
  (cycler_class_census + validation_tier_census + rediscovery + validate + schema), never a hand-picked
  subset, before commit. **Recommended model:** the row construction + gate-provenance generation is
  spec-complete mechanical work behind the census ratchet (Sonnet); the literature-clearance adjudication
  and any novelty-language wording is trust-bearing (Opus/self). **REQUIRES a Fable second-opinion pass
  before execution** — this exact chain (#561→#568) has had Fable catch real, load-bearing issues at nearly
  every stage (grid-aliasing exclusion #562, semantic max_bend inversion #565, perijove-guard + audit-field
  pins #567, duty-cycle-denominator + boundary-period-halving corrections #568); do not write to
  `catalogue.yaml` until a Fable pass has reviewed this scope + the #568 verdict it rests on.
  **FABLE PRE-EXECUTION REVIEW (2026-07-11): CONFIRMED WITH CORRECTIONS — #568's underlying data
  is fully clean (independently recomputed every duty-cycle/boundary-period/drift-floor number
  from the committed scan, zero transcription errors — a first for this chain), but the
  EXECUTION PLAN needs 4 fixes before writing, all now folded into this scope:**
  **(A) widen the test scope** — `tests/data tests/search` alone is NOT sufficient. Must also run
  `tests/test_catalogue_rediscovery.py` (lives at the tests/ TOP LEVEL, outside both planned
  dirs — pins `EXPECTED_COVERAGE[ExclusionReason.NOT_TWO_BODY]`) and `tests/verify/` (where the
  new frozen per-tier gate tests from scope-item (4) live — the plan must not skip executing the
  very gates it generates, and they must NOT be `@pytest.mark.slow` per
  [[feedback_delegation_fresh_agent_not_fork]]'s CI-skip-is-unverified-claim rule). Run at minimum
  `uv run pytest tests/data tests/search tests/verify tests/test_catalogue_rediscovery.py -q`
  (or the full suite) before commit.
  **(B) update 4 pinned counts in the SAME commit** (all independently verified by Fable):
  `tests/data/test_validation_tier_census.py::EXPECTED_TIER_CENSUS["unvalidated"]` 77→82;
  `tests/data/test_cycler_class_census.py`'s `MULTI_ARC_ALLOWLIST` id-set (add the 5 new ids,
  assuming `cycler_class: multi-arc` matching #312's own 3-Lambert-leg precedent) 292→297 +
  the `{"multi-arc": 292}` expected-dict entry → 297;
  `tests/test_catalogue_rediscovery.py`'s `NOT_TWO_BODY` count 11→16 (3-body Uranus+2-moon rows
  file here, same lane #312 itself uses); `README.md` line 41's "356-entry" → 361 (current count
  verified as exactly 356 `- id:` rows) — per [[feedback_update_docs_proactively]], same commit.
  **(C) novelty wording must NOT leak into the row as a walk-back.** #312's row currently carries
  NO `our_status` field (novelty lives in `source: discovered` + `first_published` + the notes'
  literature-novelty paragraph) — keep that paragraph fully INTACT, append the family context,
  do not rewrite. The adjudication phrase "#312 is NOT a unique novel point" is fine as internal
  reasoning but must not appear verbatim in row text — state instead that the 30-symmetric-closure
  family is an internally-enumerated fact about THIS PROJECT's own search space (#563), that #312
  is its first-documented member, and explicitly that the literature-novelty verdict is UNCHANGED
  (genericity within our own search density is not prior art — see this session's own earlier
  established PC(3,2) distinction, [[feedback_verify_novelty_against_our_status_field]]).
  **(D) literature clearance: broaden to all 5 new pairs, not just Titania-Oberon, and anticipate
  an adjacent/inconclusive result for Titania-Oberon, not a guaranteed clean not-found.**
  `src/cyclerfinder/search/literature_check.py` exists, is runnable, offline `KNOWN_CORPUS`
  mode, and its positive-control self-validation is already covered by the widened test scope
  (A) — good, satisfies [[feedback_verify_gauntlet_with_positive_control]]. But `KNOWN_CORPUS`
  now contains a Kumar-2025 Uranus-Oberon PCRTBP MMR anchor with `body_set={Oberon, Titania}` —
  exactly the Titania-Oberon pair — so the re-run may legitimately return adjacent/inconclusive
  rather than a mechanical clean not-found; treat that as an EXPECTED trust-bearing adjudication
  outcome (different topology: MMR study vs symmetric-closure quasi-cycler, per the #562 finding)
  to record explicitly in the row's notes, not a surprise blocker. Since all 5 new rows carry a
  novelty-adjacent `first_published` = cyclerfinder-discovery claim, run the offline check against
  all 5 pairs' signatures per [[feedback_literature_novelty_check_baseline]] (mandatory gate for
  any row implicitly claiming novelty), not just the one pair #565 originally flagged as
  literature-risky.
  **Also folded in (from a parallel website-visualization review, same session):** if convenient
  at no extra adjudication cost, have the row construction include `primary: Uranus` on all 6
  rows (5 new + #312's own update) and, if cheaply available from the already-computed V4 chain
  data, per-leg orbit/conic geometry — this unblocks a recommended (separate, NOT this task's
  scope) future cyclers.space hero-visualization scene without needing a later export pass. Do
  NOT let this expand #569's actual scope or delay it; skip it if it's not a trivial addition to
  what's already being written.
  **RESULT (2026-07-11): DONE, commit `8efabd5`, catalogue now 361 rows (was 356).** All 5 new
  `quasi_cycler` rows written (Titania-Oberon 74.4%/24.4d, Ariel-Oberon 71.0%/3.1d,
  Umbriel-Titania 68.5%/7.9d, Ariel-Titania 66.5%/3.5d, Ariel-Umbriel 61.7%/3.2d), #312 updated
  (family relation + 79.1%/6.0d duty cycle + its distinct drift-floor caveat, literature-novelty
  paragraph kept fully intact per correction C), frozen V4 provenance registered
  (`tests/verify/test_566_five_representatives_v4.py`, 10 tests, not slow), all 4 corrected pins
  updated (unvalidated 77→82, MULTI_ARC_ALLOWLIST 292→297, NOT_TWO_BODY 11→16, README 356→361)
  plus one MORE pin the executing agent found beyond Fable's list
  (`test_schema_v45_fields.py`'s above-V0/V4 counts). Literature check ran for all 5 pairs
  (correction D): 4 clean not-found, Titania-Oberon correctly surfaced the anticipated
  Kumar-2025/Canales-Howell-2021 adjacency and recorded the trust-bearing discrimination
  (MMR/halo-transfer vs symmetric-closure quasi-cycler) in that row's notes, as expected — not a
  blocker. The optional `primary: Uranus` website fold-in was correctly SKIPPED by the executing
  agent after testing showed it would silently break the Fable-verified NOT_TWO_BODY 11→16 pin
  (a real catch — the fold-in was explicitly scoped as skippable for exactly this reason).
  **Top-level session (not the subagent, per correction A) then ran the FULL
  `tests/data tests/search tests/verify tests/test_catalogue_rediscovery.py` ratchet suite
  directly: EXIT_CODE=0, zero FAILED/ERROR lines, 100% reached.** The #558→#569 chain is
  complete: #312's family is now fully written into the catalogue with real, validated,
  per-candidate `validity_window` data, closing out the thread that began with #558's family
  census. Follow-ups intentionally NOT started here (separate, future asks): the cyclers.space
  hero-visualization scene (assessed feasible, ~1-2 days, gated on this landing — now unblocked);
  the 2 remaining asymmetric near-closures (#565 §3, explicit NO-GO/deferred); #556/#560 parking
  lot items (unrelated to this thread).

- **#559** (P1, cheap — under a minute of compute per the #338 entry's own timing, fold into
  #558 or run standalone) — the never-dispatched #338 Phase 2 DOY-sensitivity scan. #338
  found 2000-01-15 launches FAIL V4-strict (91,000 km drift) while 2000-06-21 PASSES
  (12,000 km) — same year, and the 84-year validity window / 100% interior pass rate holds
  ONLY at June-21-ish launches. The daily-epoch scan across one Umbriel-Oberon synodic cycle
  that would resolve whether June 21 sits in a wide PASS band or on a knife-edge was scoped,
  declared non-blocking, and never run (#338's own 100-epoch annual sweep took 42s, so this
  is well under a minute of compute). Directly affects #312's mission-utility claim and its
  eventual V5 write-up (a knife-edge DOY dependency is a materially weaker result than a wide
  tolerant band).
  **Recommended model:** Sonnet (mechanical parameter sweep behind the existing frozen
  V4-strict gate, no new judgment call). [Ran as Haiku instead — the sweep itself proved
  genuinely mechanical as scoped; interpreting its result did not, see below.]
  **✓ RUN (2026-07-11, commit `976eb75`) — result needs a caveat this commit's own write-up
  did not carry.** 731 daily V4-strict epochs across 2000-01-01..2000-12-31 (366 days) and
  2030-01-01..2030-12-31 (365 days). Pass rates: 2000 89.9% (329/366), 2030 85.8% (313/365).
  June 21 (DOY 172) PASSes in both years (~14k km drift vs ~86-96k km when failing),
  corroborating #338's own anchor choice. The FAIL pattern is NOT a smooth knife-edge or a
  wide tolerant band — it is dominated by **isolated single-day FAIL spikes** (a FAIL
  surrounded by PASS on both immediate neighbors): 28/366 in 2000, 29/366 in 2030 (see
  `tests/data/test_sweep_diagnostics.py::test_559_actual_pattern_reproduces_known_singleton_count`,
  a regression pin against the real committed data). The original Haiku write-up called this
  "chaotic/stochastic, not resonant" — **that framing was never actually verified and should
  NOT be trusted as a diagnosis of the cause.** Isolated single-point flips surrounded by
  agreeing neighbors on both sides are an unusual pattern for real physical drift-vs-epoch
  sensitivity (which should vary continuously with small epoch perturbations in a smooth
  N-body system) and are a strong prior for a NUMERICAL ARTIFACT (a Lambert-solver branch
  switch, a kernel-interpolation node boundary, a discrete tolerance edge case) rather than
  genuine chaos — see `[[feedback_isolated_sweep_flips_suspect_artifact]]`. A follow-up
  diagnostic (sub-day-resolution zoom + code-path tracing), dispatched 2026-07-11, **CONFIRMED
  the artifact-not-chaos hypothesis with two concrete mechanisms, both in
  `src/cyclerfinder/data/validation/v4_uranus_strict.py`.** (1) **Discrete Lambert branch-
  selection flip** (`_cycle_v4_strict` ~line 418): the rev-1 low/high branch tie-break
  (`min()` on velocity-match residual) has no continuity tracking, and the two branches'
  match-values cross near-ties between adjacent epochs — caught directly: at 2000-04-09
  13:00->14:00 the selection flips low->high and the leg terminal miss jumps discontinuously
  23,448 km -> 1,340 km between adjacent HOURS. (2) **DOP853 integrator failure on non-
  physical planet-crossing arcs** (`_v4_strict_propagate_leg` ~lines 360-372): at epochs
  where the Umbriel-Oberon transfer angle collapses to ~3-12 deg (real geometry — the ~5.99 d
  synodic period sweeps relative phase ~60 deg/day), the fixed-14.94-day rev-1 Lambert
  "solution" is a near-parabolic (ecc~0.998-1.000) arc with perijove 0-850 km — INSIDE Uranus
  (R_eq 25,559 km); the propagated arc plunges to r~4,500-11,600 km, the integrator stiffens
  and dies (`status=-1`), and this integrator failure is silently collapsed into a plain
  validation FAIL, indistinguishable from a genuine dynamical fail — there is no perijove/
  collision guard at all. **Verdict: the #559 daily pass-rates (89.9%/85.8%) and the 55/58
  spike locations are NOT physically meaningful** — they measure daily aliasing of a
  two-discontinuity solver field (a discrete branch selection + a binary integrator-success
  gate), not #312's genuine launch-epoch robustness. **This does NOT undermine #312 itself**
  — the epoch-blind V2->V3->V4-scipy chain and the canonical single-epoch V4-strict result
  are untouched; it is a limitation of the EPOCH-SWEEP interpretation specifically. **Do not
  cite the #559 pass-rate numbers as a launch-epoch-robustness characterization for #312's V5
  write-up** until #560 (allocated below) addresses branch-continuity + degenerate-arc
  exclusion. Diagnosis-only pass; no repository code was modified investigating this.
  A reusable guard (`cyclerfinder.data.sweep_diagnostics.detect_isolated_singleton_anomalies`)
  now exists so future sweep scripts flag this pattern automatically rather than relying on a
  human noticing it — NOTE its own CAUTION: the 2000/2030 windows must be analyzed SEPARATELY,
  not concatenated (concatenating manufactures one spurious cross-boundary anomaly at the seam;
  confirmed empirically, see `tests/data/test_sweep_diagnostics.py`).
  Output: `data/silver_327_v4_strict_daily_sweep_559.jsonl`. Script:
  `scripts/run_559_silver_v4strict_daily_doy_scan.py`.

- **#560 ✓ RESOLVED (2026-07-14)** (header corrected 2026-07-15: this bullet's opening line said
  "P2, parking lot — not auto-fired" even though its own body carries a "✓ Resolved" resolution
  ~29 lines down — the same stale-header pattern as `#557`, caught by the new header/body
  consistency ratchet) — V4-strict Lambert-branch-continuity +
  degenerate-arc-exclusion robustness fixes. The #559 diagnostic (above) confirmed two real
  robustness deficiencies in `src/cyclerfinder/data/validation/v4_uranus_strict.py`, both
  needed before any epoch-SWEEP interpretation of V4-strict can be trusted (the canonical
  single-epoch V4-strict result for #312 itself is unaffected and does NOT need this fix):
  (1) `_cycle_v4_strict`'s rev-1 Lambert branch selection (~line 418) has no continuity
  tracking across neighboring epochs — near-tie flips between the low/high branches produce
  discontinuous multi-km-scale jumps in the terminal miss between adjacent HOURS, not just
  days. Fix direction: branch-continuation (track which branch was selected at the previous
  epoch and prefer it under a tie-break threshold) or report the min-drift-over-both-branches
  instead of a single discrete pick. (2) `_v4_strict_propagate_leg` (~lines 360-372) collapses
  a DOP853 integrator failure into a plain validation FAIL indistinguishable from a genuine
  dynamical fail — but the actual cause at these epochs is a non-physical planet-crossing
  Lambert "solution" (perijove 0-850 km, inside Uranus's 25,559 km radius) that only "fails"
  because the integrator chokes on it, not because of any real dynamical instability. There is
  no perijove/collision guard at all — a Uranus-crossing arc is silently accepted as a
  candidate transfer; it only fails downstream, and a more tolerant integrator could
  conceivably return a garbage PASS instead. Fix direction: add an explicit perijove-vs-body-
  radius guard on the Lambert solution BEFORE propagation, and exclude design-infeasible-
  geometry epochs from the pass-rate denominator rather than counting them as dynamical FAILs.
  **Not auto-fired**: doesn't block #312's own already-valid V4-strict result; needed only if
  someone wants to trust a future epoch-sweep characterization (e.g. re-running #559 properly
  once fixed).
  **Recommended model:** Opus (numerical-methods judgment — the branch-continuity strategy and
  the infeasible-geometry exclusion criterion are real design decisions, not mechanical fixes)
  for the design pass; Sonnet for implementation behind a regression test pinning the exact
  2000-04-09 13:00/14:00 branch-flip and the 2000-07-24/2000-08-17 planet-crossing cases the
  #559 diagnostic already located, so this fix is verifiably tested against the cases that
  motivated it.
  **✓ Resolved (2026-07-14), commit `7d3b21b` (merged `b22e4a4`).** Discovery: both robustness fixes
  were ALREADY landed under commit `6c54bba` (#567, 2026-07-11), before this task was even dispatched —
  verified independently by re-reading that commit, not taken on trust. The dispatched agent correctly
  found this instead of redoing the work, and added the one piece #560's own spec required that #567
  didn't explicitly cover: a durable regression pin (`test_560_silver_312_canonical_epoch_unchanged`)
  confirming #312's own canonical single-epoch result is unaffected. Also documented one deliberate
  deviation from #560's original wording (tag-but-count-as-FAIL, not exclude-from-denominator) that was
  itself corrected by a later Fable-validated #567 pin — flagged so it doesn't get "fixed back."
  Full closeout: `docs/notes/2026-07-14-560-v4strict-robustness-fixes-closeout.md`.

- **#541** (P3, opportunistic, lower priority than #539/#540) — First-pass Saturnian
  resonant-moon-pair screen (Mimas-Enceladus or Enceladus-Titan) using the same
  generalized pipeline — the third giant-planet moon system this project has never
  screened for whiskered-torus heteroclinic structure at all (extends the #536 Jovian
  first-pass to a second untouched system, one #536-caliber sweep, not a deep campaign).
  **Feasibility: MEDIUM** — same caveat as #539: a real sweep (band + phase), not a
  single-point probe, from the start this time.
  **Recommended models:** **Sonnet** for the sweep build/execution (spec inherited from
  #539/#540); **Opus** for the empty-region-vs-hit verdict.
  **NOTE (2026-07-12): this task depends on the whiskered-torus/linking-number method family,
  which was RETIRED for discovery this session (#555's final QUALIFIED NEGATIVE verdict — see
  [[project_novel_findings_status]]/the #548-#555 chain). #541 as scoped is effectively dead;
  see #571 below for a Saturn task using a DIFFERENT, still-live method instead.**

- **#571** (P1, plan under Fable review before dispatch — do not execute until that lands) —
  Titan-centered Saturnian moon-pair sweep using the #558-corrected discovery genome. This is
  NOT #541 (which used the now-retired torus method) — it directly reuses the SAME method that
  found #312's entire 30-member Uranian family: `scripts/scan_558_uranus_all_pairs_offset_sweep.py`'s
  corrected approach (sweep moon-pair relative offset × global phase0 × tof_scale × n_rev, not
  holding relative offset fixed by convention — the exact bug that nearly hid #312 itself),
  applied to a real system it has never touched.
  **Why Saturn, why Titan-centered**: Saturn's small/mid moons (Mimas/Enceladus/Tethys/Dione)
  already failed the #324 physical max-bend gate in prior work (#489: Tethys 0.44°, Dione 3-5°
  max-bend — mass-deficient, not a convention artifact) — that finding stands, it's a genuine
  physical limit, not something this sweep re-litigates. Titan itself is different: GM=8978.14
  km³/s² (JPL SSD, `core/satellites.py:150`) — an order of magnitude above any Uranian moon
  (Oberon GM=205.3) and comparable to Ganymede — genuine flyby-bending capacity. All constants
  for the full Saturnian moon set (Mimas/Enceladus/Tethys/Dione/Rhea/Titan/Iapetus/Hyperion) are
  ALREADY sourced with citations in `core/satellites.py` (JPL SSD gm_de440 + SAT441, various
  accessed dates) — no new capability, no new sourcing work, purely a "run where it was never
  run" application, matching #558's own framing exactly.
  **Proposed scope (subject to Fable's review before dispatch):**
  (1) Genericize `scripts/scan_558_uranus_all_pairs_offset_sweep.py`'s method (not its
  Uranus-specific constants) into a per-system-parameterized sweep, reusing its rel_offset ×
  phase0 × tof_scale × n_rev grid structure and its exact residual formula (already verified
  byte-for-byte against production `_close_one_phasing` per #558's own positive control).
  (2) Run it on Titan-centered pairs: Titan-Mimas, Titan-Enceladus, Titan-Tethys, Titan-Dione,
  Titan-Rhea, Titan-Iapetus (6 pairs — Iapetus is far out at a=3,561,700 km but has real mass,
  GM=120.5, comparable to Rhea; worth including). Titan-Hyperion excluded from the primary sweep
  (GM=0.37, matches the same mass-deficient profile that already failed at Saturn's small moons
  in #489 — include only if trivially cheap to add, not worth dedicated compute).
  (3) Gate every sub-0.05 km/s basin through the EXISTING #324 physical max-bend filter + DOP853
  cross-check, exactly as #558 did.
  (4) **Positive control**: Titan-Rhea is LITERATURE-KNOWN (Russell-Strange 2009's own published
  Titan→Rhea flyby architecture, per #489's own grounding note) — use it as the sweep's positive
  control (does the corrected genome recover a real, independently-published Titan-system
  result?) rather than treating a Titan-Rhea hit as novel; any genuinely novel candidate must
  clear `literature_check.py` same as every other admission this session.
  (5) Any sub-0.05 km/s survivor gets Opus adjudication + Fable second-opinion before ANY
  gauntlet/writeback step — same discipline as the #561-#569 Uranian chain throughout.
  **Explicitly out of scope for this task:** the V1-V4-strict gauntlet itself (a later, separate
  step if a genuine candidate survives); any catalogue.yaml edit; Titan-Hyperion unless trivial.
  **Recommended models:** Sonnet for the genericization + sweep build/run (mechanical, behind
  the same deterministic #324 gate #558 used). Opus for adjudicating any hit. Fable REQUIRED
  twice: once now, reviewing this plan before dispatch (per the user's explicit request — same
  discipline that caught #552's bad payoff assumption before a multi-week build was wasted), and
  again before any writeback if a hit survives.
  **FABLE PLAN REVIEW (2026-07-12): CONFIRMED WITH CORRECTIONS — a real, load-bearing catch,
  now folded into scope. Verified: all cited constants (Titan GM=8978.14, safe_alt 1500 km
  paper-anchored to Titan's atmosphere, other moons' GM/a) are correct, and #489's small-moon
  mass-deficiency finding stands (Tethys 0.44°, Dione ~3.1-3.2° — #571's "3-5°" was a minor
  overstatement, immaterial). But four corrections are MANDATORY before dispatch:**
  **(1) Cut Titan-Mimas / Titan-Enceladus / Titan-Tethys / Titan-Dione from the active sweep —
  they are ANALYTICALLY EMPTY under this project's own two-sided #324 gate, before spending any
  compute.** The gate requires EVERY encounter (not just Titan's) to clear 5° bend
  (`search/physical_sanity.py::candidate_passes_physical_gate`, applied to all 3 legs by
  `scan_558...::gate_candidate`) — the small moon itself must bend ≥5°. Fable computed each
  moon's minimum-achievable V∞ over ALL conics reaching Titan's radius (not just Hohmann-tangent)
  against each moon's own 5°-gate ceiling: Mimas ceiling 0.43 km/s vs min-achievable 4.54 km/s
  (0.05° bend); Enceladus 0.67 vs 3.71 (0.17°); Tethys 1.19 vs 3.06 (0.79°); Dione 1.56 vs 2.37
  (2.22°) — every one impossible at EVERY grid point, no basin can ever exist. **Stamp these 4
  pairs directly into `data/empty_regions.jsonl` as analytically-empty-by-gate-construction
  (method-versioned per [[project_negative_results_registry]]) instead of sweeping them** — do
  not spend compute re-deriving #489-with-Titan-attached. State explicitly in the stamp that this
  is conditional on the two-sided gate policy: Russell-Strange's own published Titan-Enceladus
  census (19 members) uses the small moon as a passive science target — an architecture this
  genome structurally excludes, not a contradiction of it.
  **(2) Only Titan-Rhea and Titan-Iapetus have a feasible window — sweep only these two.** Rhea
  ceiling 1.98 vs min-achievable 1.54 km/s (7.97° bend, narrow but real window 1.54-1.98).
  Iapetus ceiling 1.78 vs min-achievable 0.93 (16.4° bend, wider window).
  **(3) The positive-control claim is mis-grounded — fix before dispatch.** Titan-Rhea is NOT
  "Russell-Strange 2009's own published Titan→Rhea flyby architecture" as #571 originally stated
  — per the project's own digest, R-S's published Saturnian tables are Titan-Enceladus ONLY;
  Rhea is a named extension target with NO published invariants (this does not change Titan-Rhea's
  correct V0-known/not-novel status, #489's adjudication on that stands). The only usable
  positive control is INTERNAL: this project's own #320 scan already found a Titan-Rhea-Titan
  basin (residual 0.0316, V∞ 1.68-1.75 km/s, bends ~50°/~7-8° — sits inside the feasible window
  above) — use recovering THAT as the smoke test, and label it explicitly as an internal
  cross-lineage check, NOT a literature golden (per [[feedback_golden_tests_sourced_only]]).
  **(4) Honesty flags for the idealized→real gap, pre-register before adjudicating any hit:**
  Titan's real eccentricity ≈0.0288 (7-25× the Uranian moons' e≤0.004) — velocity modulation
  ±0.16 km/s is 3× the 0.05 km/s residual floor, so expect a materially worse idealized→real gap
  and a worse #568-style duty cycle than the Uranian family; don't oversell a coplanar-circular
  hit at adjudication. **Titan-Iapetus specifically has an inclination problem, not a tof/distance
  problem** (the a-gap synodic/transfer timing is fine, ~20-22 days): Iapetus is inclined ~15.5°
  to Titan's plane, giving an out-of-plane relative velocity (~0.85 km/s) comparable to the ENTIRE
  coplanar V∞ (0.93) — a coplanar closure at the 0.05 km/s floor is not a physical statement for
  this pair. Either drop Titan-Iapetus, or keep it under a PRE-COMMITTED rule: any "hit" there
  cannot be called a candidate from the coplanar sweep alone and needs 3D/inclined treatment
  (the #552-scoped genome extension, not yet built) before adjudication — do not let it be
  adjudicated as if it were a #558-style planar hit.
  **Net payoff, honestly restated**: this is now (a) a Titan-Rhea family census against our own
  internal #320 anchor — known-class-member payoff, zero expected novelty per the `our_status`
  discipline — plus (b) a flagged, artifact-prone Titan-Iapetus probe (needs the inclination
  caveat honored at adjudication), plus (c) four cheap analytic empty-region stamps. Smaller than
  the original "6 untested pairs" framing, but still worth dispatching — the compute is
  minutes-scale and the negative-registry stamps + Titan-Rhea census have standing value.
  **RESULT (2026-07-12): DONE, both parts executed exactly as Fable-corrected.** (A) The 4
  analytic stamps (Titan-Mimas/Enceladus/Tethys/Dione) were appended to `data/empty_regions.jsonl`
  after an INDEPENDENT recomputation (`scripts/verify_571_gate_analytics.py`, using
  `core/flyby.py::max_bend` + `core/satellites.py`'s sourced values, not copied from Fable's
  summary): Mimas ceiling 0.4290 vs min-achievable 4.5425 km/s (0.0466° bend); Enceladus 0.6701 vs
  3.7086 (0.1704°); Tethys 1.1947 vs 3.0566 (0.7930°); Dione 1.5569 vs 2.3672 (2.2170°) — all four
  independently CONFIRM Fable's figures to displayed precision. Titan-Hyperion was also checked
  (not stamped, per scope) and is likewise infeasible: ceiling 0.2367 vs min-achievable 0.2491
  km/s, 4.53° bend, just under the 5° floor. (B) `scan_558_uranus_all_pairs_offset_sweep.py` was
  genericized in place (added `primary: str = "Uranus"` default-preserving kwarg to
  `residual_at_point`/`sweep_pair`/`build_legs_for_record`/`gate_candidate`; re-ran its own
  positive control after the edit — byte-identical reproduction, 0.025232 km/s exact-point +
  0.000152 km/s grid best) and a new `scripts/scan_571_saturn_titan_pairs_offset_sweep.py` imports
  those functions with `primary="Saturn"`, sweeping all 4 directions (Titan-Rhea, Rhea-Titan,
  Titan-Iapetus, Iapetus-Titan) at the same density #558 used (n_offset=360, tof_scales dense-51,
  n_rev 0-3). **Internal #320 cross-lineage control**: EXACT match, not just "near" — evaluating
  this method at #320's own Titan-Rhea-Titan point (rel_off=285°, phase0=90°, tof_scale=2.0,
  n_rev=(1,1)) reproduces residual 0.031620 km/s to the full displayed precision (confirms this
  circular-coplanar-Kepler construction is the same underlying geometry #320 sampled at that
  point, not a coincidental nearby basin). **Titan-Rhea (V0-known, census only)**: 97 sub-0.05 km/s
  survivors passing residual+physical-bend+DOP853 gates across both directions (62 Titan-anchor +
  35 Rhea-anchor), residuals down to 0.0002-0.0006 km/s, V∞ 1.25-2.75 km/s, bends
  5-70°-ish — consistent census density with the internal anchor's basin, no new-method surprise.
  **Titan-Iapetus (novelty-eligible IF it survives adjudication)**: 187 sub-0.05 km/s survivors
  (69 Titan-anchor + 118 Iapetus-anchor), residuals down to 5.6e-5 km/s, V∞ 0.98-4.46 km/s — inside
  the analytically-predicted feasible window (0.93-1.78 km/s at Iapetus). **Every one of these 187
  is programmatically flagged in the output** (`iapetus_inclination_caveat` field in
  `data/scan_571_saturn_titan_iapetus.jsonl` / `scan_571_saturn_iapetus_titan.jsonl` and the
  `candidates_needing_adjudication` block of `data/scan_571_saturn_titan_pairs_index.jsonl`)
  restating the pre-registered caveat verbatim: Iapetus's real ~15.5° inclination to Titan's plane
  gives an out-of-plane relative velocity (~0.85 km/s) comparable to the entire coplanar V∞ budget
  at this pair, so a coplanar closure at the 0.05 km/s gate floor is **NOT on its own a physical
  statement** for Titan-Iapetus — it needs 3D/inclined treatment (the #552-scoped genome extension,
  not yet built) before ANY of these 187 can be treated as a genuine candidate, and must NOT be
  adjudicated as if it were a #558-style planar hit. Per this task's explicit scope, NO adjudication
  was performed here (Opus + Fable adjudication of any survivor, and the V1-V4-strict gauntlet, are
  separate later steps). Full ratchet (`uv run pytest tests/data tests/search -q`) passed clean
  (exit 0, no FAILED/ERROR; xfail/xpass/skip baseline unchanged) before commit.

- **#542 ✓ RESEARCH QUESTION ANSWERED, by `#608` (2026-07-16)** (P4, defer until #539-541 have added
  corrector-run diversity) — The previously-proposed #525 learned-seed generative warm-start
  (diffusion/generative model trained on the accumulated corrector runlogs/checkpoints, cf.
  Graebner & Beeson, arXiv:2501.07005) to propose seeds in unknown basins automatically. Still
  speculative; now has meaningfully more training diversity to draw on once #539-541's runlogs
  exist, but not yet worth building on #530-538's data alone.
  **Feasibility: LOW-MEDIUM, genuinely speculative** — the 2026-06-11 Ozaki
  "below-breakeven" triage already flagged this class of approach as marginal; revisit
  the triage, don't just proceed.
  **Recommended models:** research/design (surveying what training signal actually exists
  in the runlog corpus, whether it's enough) → **Opus**. If it proceeds to a build,
  implementation behind an explicit train/eval split → **Sonnet**.
  **This is exactly the research/design survey step called for above, now actually run (`#608`,
  `#605` shortlist item 3)** — not on `#539-541`'s not-yet-existing runlogs, but on the project's
  own existing `#210` outcome-log corpus (540,312 raw lines, 54,165 unique physically-sane
  converged Earth-Moon CR3BP orbits after filtering). Result: the 2026-06-11 Ozaki "below-breakeven"
  triage is **STALE, superseded by direct evidence** — a bounded numpy/scipy-only linear-Gaussian
  model (PCA + k-means + per-cluster Gaussian, not the literal diffusion/CNN-VAE this proposal
  envisioned) gave a real, measured **~12.25x** physically-sane convergence lift over blind
  uniform seeding (49% vs. 4%), independently reproduced bit-for-bit by the coordinating session.
  `#525`'s "propose seeds in unknown basins automatically" framing is not yet demonstrated (the
  #608 corpus/test is single-family Earth-Moon CR3BP, not a genuinely unknown/novel basin), but the
  core mechanism — a statistical density over converged genomes steering a corrector far more
  often than blind search — is no longer speculative. **Not closing this task** (the broader
  "propose seeds in an unfamiliar basin" claim remains untested), but the feasibility triage above
  should be read as ANSWERED, not open, and any future revisit of this idea should start from
  `#608`'s code/corpus (`src/cyclerfinder/ml/orbit_generative.py`,
  `data/found/608_generative_seed_poc/`) rather than re-litigating feasibility from scratch. #608's
  own honest verdict flags two concrete next steps if this is to go further: (1) tag
  provenance/family per training record (the current corpus has none, only primary/secondary),
  and (2) a nonlinear model (kernel PCA / shallow autoencoder), since the linear PCA model
  demonstrably lands in physically-plausible "gaps" between the true curved family manifolds
  rather than tracking them tightly — neither has been scoped as its own task number yet.

- **#543 ✓ SCOPED (2026-07-11)** (header corrected 2026-07-15: this bullet's opening line said
  "parking lot — needs a scoping conversation, not a sprint slot" even though its own body records
  the scoping conversation actually happened and was resolved ~10 lines down — same stale-header
  pattern as `#557`/`#560`, caught by the new header/body consistency ratchet) — #529's
  deferred inter-cycler-network / taxi-transfer catalogue extension. Still real (this
  project's own corpus documents published fleet-of-cyclers networks the current 5-class
  taxonomy cannot express), still ranked last deliberately: it changes what a catalogue
  row/relation can *be*, which wants a deliberate human scoping discussion before any
  agent starts building.
  **Recommended model:** the scoping conversation itself should happen with the user
  directly; once scoped, schema work → **Sonnet**, search-over-existing-rows logic →
  **Opus** (it's a discovery-adjacent judgment task, reusing `model/score.py`'s
  `taxi_cost_kms`, not mechanical).
  **✓ SCOPED (2026-07-11, direct conversation with the user) — general schema infrastructure
  ONLY, no populated network rows yet.** Three concrete options were laid out (general
  infrastructure / reproduce Sanchez Net's own example network / discover emergent networks
  among our own rows); the user chose general infrastructure. Reproducing Sanchez Net's own
  example network is explicitly NOT this task's job — their named network cyclers (IDs 51, 84,
  1141, etc. in `docs/notes/s1l1-target-topology-mining.md:156-263`) have no recoverable
  orbital elements (only downlink date + cycler index + ME transit duration are published), so
  there is nothing citable to reference yet regardless of schema readiness.

- **#570** (P2, schema infrastructure — #543's scoped first slice) — add a general cycler-network
  relation to the data model: a NEW top-level registry, separate from `data/catalogue.yaml`
  (matching the `data/empty_regions.jsonl` precedent of a separate registry rather than
  overloading the per-orbit-row schema), representing a SET of catalogued cyclers with shared
  phasing/downlink cadence, sourced from the Sanchez Net et al. 2022 EMCO fleet concept
  (`docs/notes/s1l1-target-topology-mining.md:156-263`) but NOT populated with any real network
  yet — that's explicitly deferred (see #543's scoping note above).
  **Design correction found while scoping (verify in your own dispatch, don't just trust this
  summary): `model/score.py::taxi_cost_kms` is a per-cycler INSERTION-cost surrogate** (Earth →
  that single cycler's own Earth encounter, `max(||enc.vinf_in||)` over Earth encounters) — it is
  NOT a cost between two different cyclers. #529/#543's original phrasing ("taxi transfers
  between catalogued cyclers") overstates what this function computes. The schema must reflect
  what's actually computable now (per-member insertion cost via the existing function) and must
  NOT imply a between-cycler transfer cost field that has no computation behind it — if a genuine
  inter-cycler transfer cost is wanted later, that's new M5/M6-tier search work, a separate task,
  not something to fake into this schema now.
  **Scope:**
  (1) New schema file `data/cycler_network.schema.json` (JSON Schema, mirroring
  `catalogue.schema.json`'s conventions — versioned, documented `description` fields) for a
  `data/cycler_networks.yaml` registry. Minimum fields per entry: `id`, `name`, `source`
  (`literature`/`derived`), `member_cycler_ids` (array of `catalogue.yaml` row ids — every one
  MUST exist, enforced by a validator, not just documented), `downlink_cadence` (free-text
  description + an OPTIONAL structured `schedule` array of `{date, cycler_id, me_transit_days}`
  for when real sourced schedule data exists — null/absent is honest, not an error, matching this
  project's `data_gap` convention), `per_member_taxi_insertion_cost_kms` (array of `{cycler_id,
  cost_kms, taxi_body}` — MUST be computed by calling `taxi_cost_kms` on the referenced row, not
  hand-typed), `data_gaps`, `first_published` (when `source: literature`), `notes`.
  (2) A validator module (e.g. `src/cyclerfinder/data/validate_networks.py`, mirroring
  `validate.py`'s pattern) checking: every `member_cycler_ids` entry resolves to a real
  `catalogue.yaml` row; every `per_member_taxi_insertion_cost_kms` entry's `cost_kms` matches a
  fresh call to `taxi_cost_kms` on that row (a cross-check, not a duplicate hand-maintained
  value — catches drift if the referenced row's encounters ever change); JSON-schema-valid.
  (3) Tests: a self-consistency/positive-control test constructing a small SYNTHETIC test network
  (2-3 members) referencing REAL, already-catalogued rows (e.g. reuse a couple of the #566/#569
  Uranian family rows, or any well-known validated rows — your choice, just real ids) to prove
  the schema + validator + taxi-cost cross-check actually work end-to-end — this test network is
  NOT a real "discovered network," just infrastructure proof, and must be clearly labeled as such
  (e.g. `source: derived`, a notes field stating it's a schema self-test, not a Sanchez-Net-style
  finding) so it can never be mistaken for a real result later.
  (4) `data/cycler_networks.yaml` itself ships EMPTY (`[]`) or with only the labeled
  self-test entry from (3) — per the user's explicit scoping choice, no attempt to populate a
  real Sanchez Net (or any other) network in this pass.
  **Explicitly out of scope:** reproducing Sanchez Net's own example network (blocked by a real
  external data gap, not a schema limitation — see #543's note); any new inter-cycler transfer-
  cost computation beyond the existing per-member insertion-cost reuse; any change to
  `catalogue.yaml` itself or its own schema; website/cyclers.space changes.
  **Recommended model:** Sonnet (spec-complete schema + validator + test infrastructure, design
  already settled by this scoping pass — no open judgment calls left for the implementer).
  **✓ RESULT (2026-07-11):** Built exactly as scoped, infrastructure only. `data/cycler_network.schema.json`
  (v1.0, JSON Schema mirroring `catalogue.schema.json`'s conventions) + `data/cycler_networks.yaml`
  (ships genuinely EMPTY `[]` — the self-test network lives entirely in the test file's own
  fixtures, not the committed registry, per the task's "your choice" guidance, to remove any risk
  of it being mistaken for a real finding later). Validator
  `src/cyclerfinder/data/validate_networks.py` runs four never-raising layers (mirroring
  `validate.py`'s pattern): JSON-Schema structural (via the `jsonschema` package directly, not
  just the out-of-process `check-jsonschema` hook), semantic (id uniqueness,
  `source=literature`⇒`first_published` honesty gate), referential (every `member_cycler_ids` /
  schedule / cost `cycler_id` must resolve to a real `catalogue.yaml` row — never silently
  skipped), and a taxi-cost cross-check that reconstructs a MINIMAL `Cycler` from the referenced
  row's published V∞ magnitudes (`_catalogue_entry_to_taxi_cycler` — no richer catalogue-row-to-
  `Cycler` adapter exists anywhere in the codebase; documented as narrowly scoped to the one
  function it needs, `taxi_cost_kms`, which only reads `.body` + `||vinf_in||`) and asserts the
  stored `cost_kms` matches a fresh `taxi_cost_kms` call within 1e-9 tol. Self-test network
  (`tests/data/test_cycler_networks.py`, 18 tests) references three REAL catalogue rows —
  `aldrin-classic-em-k1-outbound` (V2, E=6.5 km/s), `aldrin-classic-em-k1-inbound` (V1, E=6.5
  km/s), `russell-ocampo-2.5.1+0` (V1, E=7.8 km/s) — with every `cost_kms` computed by calling
  `taxi_cost_kms` in the fixture itself (never hand-typed); `source: derived`, `notes` states
  "SCHEMA SELF-TEST / INFRASTRUCTURE PROOF ONLY... NOT a real discovered or sourced cycler
  network." Positive control (full gate passes clean) + 9 negative controls (unresolvable member,
  taxi-cost drift ×2, literature-without-citation, duplicate id, non-member cost/schedule
  references, 3 schema-shape rejections) all pass. Added a `check-jsonschema` pre-commit hook +
  a `cycler-network-ratchets` local hook for the new files (mirrors the catalogue hooks).
  Side-effect: added a `jsonschema.*` mypy override (first direct Python import of `jsonschema`
  in this codebase) and removed 9 now-unused `# type: ignore[import-untyped]` comments on
  pre-existing `import jsonschema` lines that the override made redundant
  (`warn_unused_ignores = true` would otherwise fail CI). Full `tests/data tests/search -q`
  ratchet suite: clean (exit 0, no FAILED/ERROR, only pre-existing documented XFAIL/XPASS/SKIPPED
  entries). No real network populated — #543's external data gap (Sanchez Net's example cyclers
  have no recoverable orbital elements) remains open and unaddressed by design.

**Recommendation:** #538 first (unblocks all of the above — #539/#540 both reuse its
corrector). #539 next (cheapest genuine new-territory attempt, directly closes the
single-point-negative gap #536 left open). #540 alongside or immediately after (highest
strategic value — it's the only system with a proven track record of yielding a real
novel row). #541 opportunistic filler. #542/#543 stay parked until there's more
runlog diversity / a scoping conversation respectively.



## DELTA 2026-06-30 (session B — discovery campaign + Ross-corpus acquisitions) — read this first

Multi-thread campaign off the "novel cyclers are the deliverable" frame, after #492 established the
ideal-model moon-cycler frontier is exhausted (novel ground is now capability-gated). Landed:
- **#493** L-L 2011 IEG (EIGE + GIPEIPE) reproduced in the ideal model; period <1% vs sourced, but
  **adjudicated V0 (override of the agent's V1)** — sub-surface flybys / different member, same wall as
  #480; characterization kept (`afc7a64`, verdict note has the adjudication).
- **#494** binary (k₁,k₂)-cycler μ-family — REFRAMED (not new dynamics; the CR3BP solver is μ-agnostic and
  the EM slice is already V2). Acquired+mined **Ross-Roberts-Tsoukkas 2026 (arXiv:2606.29189)** Table I
  golden (`data/golden/ross_rt_2026_cycler_families.yaml`); Phase-0 positive control GO (5/5 EM recovered,
  `f872d07`); Phase 2/3 (μ-extension closing #315/#252/#255 + Pluto-Charon μ=0.1085) RUNNING. Design
  `docs/superpowers/plans/2026-06-30-494-binary-cycler-mu-family-design.md`.
- **#495** adopted Braik-Ross 2026 MIT-repo goldens — **C21 recovers at exact CJ=3.129389531054557** (the
  rounded 3.1294 gives a >5-d-different orbit → confirms the #249 bug + [[feedback_published_rounded_values_are_display]]);
  proxy ΔV is a strong screen (ρ=0.96) not a strict bound (`2b589b5`/`9b0d88c`/`42b3b16`). Repo survey
  `docs/notes/2026-06-30-binbraik-orbital-network-repo-reuse-survey.md`.
- **#318 Phase 2b** Sobol joint-search built + positive-control PASS + smoke = clean empty (`b262284`).
- Corpus: filed+mined Ross-RT 2026; **Kumar 2509.12675 was a DUPLICATE** (already digested 2026-06-20 +
  coded in resonance_network.py) — re-processing reverted, stale "undigested" index label fixed (`cfcd7e5`).
- Reviewed the FULL ross.aoe.vt.edu/papers list; the relevant gaps are allocated below. EXPLICITLY
  OUT OF SCOPE (dynamical-systems method in non-astro domains, NOT acquired): atmospheric LCS/FTLE,
  gliding & biomechanics, microbial dispersal, ship-motion & snap-through-buckling tube dynamics,
  chemistry non-RRKM rates, formation-flying/alpha-shapes.

**TASK ALLOCATIONS (next-unused per [[project_task_numbering_convention]]; #487-#496 used earlier; #497-#503 here):**
- **#497** — #249 gate recalibration: set `DV_CAP_MS` ≈ 51 m/s (from #495: C32 wins centrality in 80% of
  Braik's budget sweep there; our xfail gate uses 409.3 = near-full connectivity, no betweenness signal) →
  flip the C32-dominance gate.
  **TESTED AND REFUTED (2026-07-15).** Before dispatching this as "ready," checked the actual test
  file (`tests/search/test_reachable_network_gate.py`) — found the premise was already stale:
  #513's later R52-U node recovery already got the gate to 2/3 metrics matching Braik's Table 4
  (strength, harmonic-closeness), leaving only betweenness (C21 narrowly beats C32, ~9%), still
  computed at the ORIGINAL 409.3 cap. Ran the actual experiment #497 proposed (recalibrating
  `dv_cap_ms`) on the R52-U-recovered 13-node network, sweeping 51.16-409.3 m/s: **C32 never wins
  betweenness at any tested value** — at the tightest, Braik-matching cap it gets WORSE (betweenness
  drops from 0.3485 to 0.1364, C11a takes over entirely). This is a clean, decisive negative, not
  an unconfirmed follow-on: cap recalibration is not the fix. Updated the xfail test's own reason
  to record this (so it isn't re-proposed) and confirmed the residual cause is genuinely
  proxy-fidelity/topology-difference from Braik's own DVmatrix, not a tunable parameter — that
  deeper investigation (edge-by-edge comparison against Braik's DVmatrix.csv) is NOT done here,
  flagged as the real next step if anyone wants to pursue it further.
  `docs/notes/2026-07-15-497-dv-cap-recalibration-refuted.md`.
  **#497 STATUS: CLOSED (refuted, not resolved).**
- **#595** (P2, repo-hygiene / documentation-accuracy) — discovered the #498/#499/#503 acquisition
  status in this file was stale (see the #498/#499/#503 entries below for the full detail): the user
  re-uploaded 5 PDFs already byte-identical to files acquired 2026-06-30 (commit `6dfedab`, 16-paper
  ross.aoe.vt.edu batch); all 12 #498/#499 papers and all 4 #503 papers were already acquired,
  digested, and `CORPUS_INDEX.md`-registered that day but this tracker was never updated.
  **#595 STATUS: CLOSED** — #498/#499 marked CLOSED, #503's acquisition sub-task closed (its
  "expand #267 goldens" mining step remains open), #500's #498-gate cleared.
- **#498** — Acquire+mine MOON-TOUR / GRAVITY-ASSIST Ross papers.
  **STATUS CORRECTED 2026-07-15: this "(all MISSING)" framing is STALE.** All 6 papers were already
  acquired directly from ross.aoe.vt.edu AND digested on 2026-06-30 (commit `6dfedab`, 16-paper batch),
  fully registered in `docs/notes/CORPUS_INDEX.md`: "Design of a multi-moon orbiter" (Ross-Koon-Lo-Marsden
  2003, AAS 03-143 → `2026-06-30-digest-multi-moon-orbiter-2003.md`); "Multiple gravity assists, capture,
  and escape in the RTBP" (Ross-Scheeres 2007, SIADS 6(3) → `2026-06-30-digest-keplerian-map-2007.md`);
  "Constructing a low energy transfer between Jovian moons" (KLMR 2002, Contemp. Math. 292 →
  `2026-06-30-digest-jovian-transfer-2002.md`); "Controlled Keplerian map" (Grover-Ross 2009, JGCD 32(2) →
  `2026-06-30-digest-grover-ross-2009-controlled-keplerian-map.md`); "Resonance and capture of Jupiter
  comets" (KLMR 2001, CMDA 81 → `2026-06-30-digest-koon-2001-resonance-capture-jupiter-comets.md`);
  "Geometric mechanics and the dynamics of asteroid pairs" (Koon-Marsden-Ross-Lo-Scheeres 2004 →
  `2026-06-30-digest-koon-2004-geometric-mechanics-asteroid-pairs.md`, also #494/#308). Rediscovered
  2026-07-15 when the user uploaded 5 of these same PDFs a second time (independently downloaded, byte-
  identical to 4 already-committed files) — confirms genuine duplicates, not new content.
  **#498 STATUS: CLOSED (acquired + digested, tracker was stale).** Un-gates #500 (below).
- **#499** — Acquire+mine HETEROCLINIC FOUNDATIONS.
  **STATUS CORRECTED 2026-07-15: this "(MISSING)" framing is STALE**, same 2026-06-30/`6dfedab` batch:
  "Heteroclinic connections between periodic orbits and resonance transitions" (KLMR 2000, Chaos 10(2)
  427-469 → `2026-06-30-digest-klmr2000-heteroclinic-chaos.md`, mined-by-#314); "Heteroclinic Transfer
  Between L1 and L3 in EM" (Braik-Ross 2025, AAS 25-716 → `2026-06-30-digest-braik-ross-2025-heteroclinic-
  L1-L3.md` — acquired directly, no DOI needed, despite this paper having zero CrossRef record);
  "Connecting orbits and invariant manifolds in the spatial RTBP" (Gomez et al. 2004, Nonlinearity 17 →
  `2026-06-30-digest-gomez2004-spatial-rtbp.md`); "The Genesis trajectory and heteroclinic connections"
  (KLMR 1999, AAS 99-451 → `2026-06-30-digest-genesis1999-heteroclinic.md` — also acquired directly,
  no DOI); "Transport of Mars-crossing asteroids from the quasi-Hilda region" (Dellnitz et al. 2005,
  PRL 94 → `2026-06-30-digest-dellnitz2005-mars-transport-gaio.md`, also #308); "Experimental validation
  of phase space conduits" (Ross-BozorgMagham-Naik-Virgin 2018, PRE 98 052214 →
  `2026-06-30-digest-rbnv2018-phase-space-conduits.md`, underpins the #494 construction).
  **#499 STATUS: CLOSED (acquired + digested, tracker was stale).**
- **#500 STALE DUPLICATE, CORRECTED (2026-07-15)** — this line ("#500 itself has NOT yet been
  dispatched/run") is WRONG and predates its own resolution ever being written back into this file;
  see the actual **`#500 DONE`** entry earlier in this file (near the #498/#499/#503 block): genome
  built, 17/17 sourced positive controls pass, verdict = **clean negative for catalogue** (it's a
  Jovian moon-to-moon sub-leg planner, not a standalone interplanetary cycler genome — see
  `docs/notes/2026-06-30-500-keplerian-map-genome-verdict.md`). Caught 2026-07-15 when a Fable
  planning pass read only this stale line and recommended re-dispatching #500 against Saturn/Uranus
  as a "zero-blocker" novel-orbit path; verified against the actual verdict note before acting on
  that recommendation. Any future #318-joint-search chaining of the Keplerian map (the verdict note's
  own suggested next step) is a NEW task, not a re-dispatch of #500.
- **#501** — #318 full-scale, REFRAMED: broaden sequences/systems (NOT densify CGCEC — smoke showed a mined
  regime, 0/256) joint-search campaign with post-hoc lit-novelty.
  **STALE "[held]" MARKER CORRECTED 2026-07-15 (Fable adversarial pass, full-file audit): this
  "[held — the deferred #318 decision]" tag is stale — #501 was actually dispatched and run to
  completion regardless. See the Tier-1 results block near the top of this file: "#501 broadened
  real-eph Galilean joint search (e5b2339/d0c8def): 6 sequences × 512 Sobol = 3072, 213 feasible,
  26 shot, 0 closed, positive-controlled against Liang Member D — PASSED, trustworthy. Clean
  empty-region map registered." #576/#577 elsewhere in this file both treat this 0/3072 result as
  settled history. RESOLVED/CLOSED, not held.**
- **#502** — Watch/acquire the LONGER companion of Ross-RT 2026 (full μ-grid + all 9 EM family tables) for a
  tighter #494 Pluto-Charon continuation seed. [watch]
- **#503** — Acquire+mine cislunar-resonance / transport companions: Rawat-Kumar-Rosengren-Ross 2026 JGCD
  49(4) "Cislunar Mean-Motion Resonances"; Onozaki-Yoshimura-Ross 2017 ASR 60 (4-body tube dynamics);
  Fitzgerald-Ross 2022 ASR 70 (periodically-perturbed RTBP, → #292/#293); Naik-Lekien-Ross 2017 RCD 22(3)
  "Phase Space Transport / Lobe Dynamics" (+Lober software). → expand #267 resonance_network goldens.
  **ACQUISITION DONE 2026-07-15 (found stale, same #498/#499 audit): all 4 papers were already acquired +
  digested 2026-06-30 alongside #498/#499 (commit `6dfedab`) — `2026-06-30-digest-rawat2026-cislunar-mmr.md`,
  `-onozaki2017-tube-4body.md`, `-fitzgerald2022-transit-perturbed-rtbp.md`, `-naik2017-lobe-dynamics-
  transport.md`. Acquisition sub-task closed; the actual "expand #267 resonance_network goldens" mining
  step is NOT yet done — that remains open.** [lower priority]
- **#596** (P2, data-integrity/documentation-accuracy) — `data/MISSING_DATA.md` staleness correction.
  Same audit method as #595, triggered by the user's follow-up "do the staleness checks" request.
  Every source the 2026-06-03 report flagged PAYWALLED/RESTRICTED/NOT-ONLINE (Russell 2004
  dissertation, McConaghy 2005/2004 Purdue dissertation, Friedlander/Niehoff 1986, Jones/Hernandez/
  Jesick 2017 VEM, Russell & Strange 2009) has actually been acquired and read since — confirmed via
  `data/catalogue.yaml` cross-check, NOT just corpus-index status. **But the underlying gap counts are
  mostly still real**: ~216 Russell-family entries genuinely still carry open `trajectory.segments[].
  a_au`/`.tof_days` gaps (the 2026-06-07 mining pass extracted Russell's methodology + verbatim summary
  tables but explicitly deferred the mechanical per-row catalogue backfill); a previously-untracked
  ~38-row Table-3.4 gap (catalogue-eligible cyclers never added) was found in the same transcription
  note; the Aldrin-establishment/U0L1/VISIT-V∞ precision gaps are now CONFIRMED genuine dead ends
  (both Rogers 2012 and McConaghy's actual dissertation were read specifically for these values and
  don't have them) rather than "go find the source" items. `data/MISSING_DATA.md` rewritten in place
  with a corrected §1 summary table; §3/§4/§6 left as a historical record (source-status stale,
  numeric-gap structure still valid). **Follow-on work identified but NOT done here:** the actual
  Russell Table 3.4/3.9–3.11/4.9–4.13 per-row backfill (mechanical, source in hand, ~216+38 rows) is
  a real, ready-to-execute task — no longer blocked on access.
  **PILOT RUN 2026-07-15, CLEAN NEGATIVE — NOT actually mechanical.** Tried inverting Table 3.4's
  (AR, ToF, V∞_E, V∞_M) columns to `(a_au, e)` via `free_return_geometry` + least_squares against the
  Aldrin row (positive control: known sourced `a=1.60/e=0.393`). Neither a local fit near the known
  answer nor a 29-point global multi-start search found any `(a,e)` that simultaneously satisfies all
  three Table 3.4 columns — ToF and V∞ are mutually inconsistent under this project's radial-crossing
  convention (best fits miss one or the other by 15-20%). Leading unverified hypothesis: Table 3.4's
  own header "Earth→Mars **(or aphelion)** Time" suggests the ToF column measures time-to-aphelion for
  some rows, not always time-to-Mars-crossing — needs Russell Ch.3.5-3.8 re-read for the exact
  definition before any inversion is trustworthy.
  **RESOLVED 2026-07-15 (same day, follow-up re-read).** Re-read Russell Ch.3 (printed pp.56-69):
  confirmed AR = aphelion/1.52 AU (his own rounding). Dropped the ToF constraint from the fit; using
  only **AR + V∞_Earth** to invert `(a_au, e)` via `free_return_geometry` + least_squares works and
  is DOUBLY validated: (1) against Aldrin's known sourced `a=1.60/e=0.393` — fit recovers
  `a=1.6038/e=0.3932` (0.24%/0.05% error); (2) against cycler `2.5.1.+0`, where Table 3.5 gives the
  actual 3D v∞ VECTOR at departure, letting `(a,e)` be computed directly via two-body state-vector
  orbit determination (no fitting at all) — ground truth `a=1.5651/e=0.4010`, the AR+V∞ inversion
  independently recovers `a=1.5633/e=0.4001` (0.12%/0.22% error). The `2.5.1.+0` case ALSO confirms
  `tof_em_days` correctly reproduces Table 3.4's ToF column (94.37 emerged vs. 94 tabulated) when fed
  the true `(a,e)` — so the earlier "ToF mismatch" was an Aldrin-specific quirk (Aldrin is an
  externally-sourced named cycler, footnoted as such, not one of Russell's own generic-return search
  results), not a general modeling failure. **Practical backfill recipe validated:** `tof_days` —
  cite Table 3.4/3.9-3.11's own column directly (matches existing Aldrin catalogue precedent exactly);
  `a_au`/`e` — derive via the AR+V∞_Earth inversion (`kind: derive`); V∞ at Mars serves as a third,
  free, non-imposed cross-check on every row. Full writeup:
  `docs/notes/2026-07-15-596-russell-backfill-method-validated.md` (supersedes the same-day
  `...-pilot-inversion-fails.md` note, kept for the record).
  **EXECUTED 2026-07-15** (`scripts/backfill_russell_2004_tables.py`). Parsed all 201 Table 3.4/
  3.9-3.11 rows; 197 matched an existing catalogue entry (only 3 genuinely new/uncatalogued rows
  remain, not ~38 as the stale 2026-06-07 transcription-note estimate said — most had already been
  added since then). Of the 197, the inversion cleanly converged (solver cost < 1e-4, V∞-Mars
  cross-check < 5%) on **161**; the other **36 all have AR < 1.0** — a genuine model boundary (the
  ellipse's aphelion doesn't reach Mars for near-ballistic AR<1 cyclers, so the crossing-based
  inversion doesn't apply), confirmed by exact 1:1 correlation with `row.ar < 1.0`, correctly left as
  open gaps rather than forced. **Second real bug caught before writeback**: `orbit_elements.a_au`/`e`
  at the top level is schema-restricted to `cycler_class: single-ellipse`; every Russell-family entry
  is `multi-arc` (out-em/ret-me/loop-ee-* can be different ellipses), so writing top-level
  orbit_elements would have been a genuine schema violation — caught by
  `test_validate_catalogue.py`'s combined gate on the first write attempt, fixed to only write
  per-segment `trajectory.segments[out-em].a_au`/`e` (+ `ret-me` for the 4 `h=0` simple cyclers).
  **Also caught before commit**: a plain PyYAML (and even a naive full-file ruamel.yaml) round-trip
  reformats/strips comments across the ENTIRE 54k-line file, not just the touched rows — worked
  around via a 3-way patch strategy (diff two ruamel round-trips against each other to isolate ONLY
  the backfill's actual changes, then apply that as a patch onto the pristine original file). Final
  writeback: 161 rows backfilled (~5400-line diff, confined to the touched rows), zero cross-check
  flags, `uv run pytest tests/data tests/search tests/scripts` clean except the 2 already-documented
  pre-existing local-Mac failures (task #584). The 36 AR<1 rows and the 3 genuinely-missing rows
  remain open for a future pass (needs modeling the near-ballistic powered-assist geometry, which
  this validated method doesn't cover).
- **#597 ✓ DONE (2026-07-15)** (P3, corpus acquisition + full mining pass) — 4 more Ross-group papers
  found via a manual review of `https://ross.aoe.vt.edu/papers/` (user-suggested, same #595/#596
  session): Kumar-Rawat-Rosengren-Ross 2024 IAC-24-C1.9.5 (interior 4:1/3:1/2:1 MMR heteroclinic
  connections — predecessor to the already-anchored ASR 77:3815 journal paper); Rawat-Kumar-Rosengren-
  Ross 2025 AAS 25-569 (exterior-MMR sequel); Rawat-Kumar-Rosengren-Ross 2024 AAS 24-368 (quantitative
  resonance-widths/chaotic-zone companion); Rosengren-Ross-Kumar-Rawat 2024 AMOS (xGEO domain-awareness
  survey). All confirmed genuinely new/distinct via title-page read (not filename inference), all
  filed to `cyclers_pdf/papers/` + `CORPUS_INDEX.md`-registered.
  **Full mining pass complete 2026-07-15** (all 4 read page-by-page, not just abstract/intro; digest
  notes `2026-07-15-{kumar,rawat,rawat,rosengren}-2024/25-*.md` rewritten with the actual quantitative
  findings): IAC-24-C1.9.5 — 3:1 identified as the "gateway to the Moon" (heteroclinic to 2:1 for
  `C<=3.15`), 4:1 isolated by an RIC (rotational invariant circle) barrier at every tested `C`; AAS
  25-569 — exterior 1:3/1:4 asymmetric-bifurcation topology + `L1`-`L2` tube heteroclinics linking
  exterior 1:3 to interior 2:1/3:1 (7-9.5 day transfers at `C=3.10`, vanish by `C=3.15`); AAS 24-368 —
  CR3BP resonance widths TLE-validated against real IBEX/TESS/Spektr-R orbits (broader than Gallardo's
  semi-analytical prediction), explicit 3:1<->2:1 transfer times (28-29d direct / 56-57d via
  intermediate 5:2); AMOS — rigorous Laplace-radius xGEO definition (`~7.7 R_Earth`) + `REBOUND`/
  `ASSIST` real-spacecraft long-horizon validation (OGO-1's full 56-year trajectory, incl. 2020
  reentry, recovered from 1964 launch conditions alone). 4 new `CorpusAnchor`s registered in
  `literature_check.py` (`kumar-2024-interior-mmr-iac`, `rawat-2025-exterior-mmr-aas`,
  `rawat-2024-resonance-widths-aas`, `rosengren-2024-xgeo-resonant-structure-amos`, all
  `verified-against-source`); `tests/search/test_literature_check.py` re-verified green. **Follow-on
  opportunity found, NOT actioned here** (out of scope for a digest task) — see #598.
- **#598 ✓ DONE (2026-07-15) — the "missing PDF" was NOT missing, and fixing it up found a real bug** — `src/cyclerfinder/
  search/resonance_network.py` (#267 Track-B tier 3) carries a documented "REPRODUCE-BEFORE-TRUST"
  data gap: its period test is `xfail`-marked with the reason "Kumar 2025 PDF not held in our local
  mirror" (Kumar-Rawat-Rosengren-Ross, "Cislunar Resonant Transport and Heteroclinic Pathways: From
  3:1 to 2:1 to L1," arXiv:2509.12675, published *Advances in Space Research* 77(3):3815 (2026), DOI
  `10.1016/j.asr.2025.12.005`). **User asked to verify this claim (2026-07-15) — it is STALE, not
  true.** The PDF is already sitting in the corpus at `cyclers_pdf/papers/kumar-2025-arxiv-2509.12675.pdf`
  (47 pages; content-verified against its own title page: exact title, author list, and abstract all
  match) — found via a corpus content grep (`pdftotext` on every PDF, not filename inference), per
  the standing `feedback_corpus_check_index_not_filenames` discipline. It must have been acquired in
  a later batch (#498/#499/#503 or #595-#597) without `resonance_network.py`'s docstring/xfail gate
  ever being updated to notice — the exact same "acquired but never reconciled" failure mode as #596's
  `MISSING_DATA.md` staleness and #595's `OUTSTANDING.md` staleness, just localized to a module
  docstring instead of a tracking doc this time.
  **RESOLVED 2026-07-15** — read the full 47-page PDF page-by-page and found a REAL bug the stale
  docstring had been masking: `_RESONANT_SEEDS["R31-U-Kumar"]`'s seed (`x0=0.354146033959, sign=+1`)
  matched NO row of the paper's own Table 6 (Appendix 8.2, exact converged-orbit ICs) — a stale
  placeholder from before the PDF was in the corpus, never reconciled. Also, ALL THREE `"*-Kumar"`
  members were being recovered at the module's default `c_j=C_J_BRAIK_ROSS` (3.1294, a DIFFERENT
  paper's energy) instead of Table 6's own sourced Jacobi constants (C=3.10 for 3:1/2:1, C=3.15 for
  4:1) — silently reconstructing a nearby but different periodic orbit than Table 6 documents. Fixed
  both: corrected R31-U-Kumar's seed to Table 6's actual C=3.10 row (`x=0.822429022871,
  ydot=-0.300987128481`), added `KUMAR_TABLE6_CJ` mapping so callers pass the right energy per label.
  Added a NEW, stronger reproduce-before-trust test (`test_reproduce_kumar_table6_ic`, parametrized
  over all 3 Kumar members) checking the recovered `(x0, ydot0)` against Table 6 verbatim — converges
  to ~1e-13 relative error, the strongest possible confirmation. Un-`xfail`ed
  `test_reproduce_r41u_period_kumar` (was already unexpectedly XPASSing; the paper's Figure 7 caption
  period is an exact printed value, not a "digitization/estimate" as the stale xfail reason claimed).
  Also fixed a genuine misnomer: `kumar_equinoctial_metric` → `kumar_angular_momentum_laplace_metric`
  (the function correctly implements the paper's Eq. 10 angular-momentum/Laplace-vector distance; it
  was just named after the WRONG equation — equinoctial elements are a separate, unused alternative
  metric in the paper's Appendix 8.1 Eq. 13). All stale "PDF not held"/"xfail by design" comments in
  both the module and test docstrings rewritten with the correct citation (DOI
  `10.1016/j.asr.2025.12.005`). Verified: `tests/search/test_resonance_network.py` 13/13 passed,
  `ruff format`/`ruff check` clean, `mypy src tests` (full CI-equivalent invocation) exit 0. The #597
  mining pass's 2 conference precursors (IAC-24-C1.9.5, AAS 24-368) remain useful corroborating
  context but were not needed as a stand-in once the actual target paper's real bug was found and
  fixed directly.
  **Bugfix-invalidation audit (2026-07-15, per [[feedback_bugfix_invalidates_past_searches]]):
  ZERO blast radius beyond the test file — checked, not assumed.** A discovery-strategy planning
  pass flagged this fix as "feeds `five_tier_prioritizer.py`/`single_orbit_prioritizer.py`/
  `ftle_scorer.py`/`lobe_overlap_scorer.py`/`neural_reach_prefilter.py` — any past negative screened
  through those tiers is a potential false negative." Verified directly rather than acting on the
  claim: all five modules (+`genome/heteroclinic_cycle.py`) DO import `resonance_network.py`, but
  NONE of them call `recover_resonant_family` with a `"*-Kumar"` resonance string or touch
  `_RESONANT_SEEDS["R31-U-Kumar"]`/`KUMAR_TABLE6_CJ` — every one either builds its own
  `ResonantMember` directly from its own candidate orbits, or uses the Braik-Ross (non-Kumar) seeds,
  which the bug never touched. The `"*-Kumar"` path is exercised ONLY by
  `tests/search/test_resonance_network.py`'s own fixture. **No past discovery search, sweep, or
  catalogue decision ran through the buggy code path — nothing to re-run.** The planning pass's
  "feeds X/Y/Z" framing conflated "imports this module" with "consumes this specific buggy data,"
  the same category of overgeneralization (smaller stakes) as the `#516`/`#517` mixup nearby.
- **#599 ✓ DONE (2026-07-15) — capability gap fixed, sweep run, clean negative** — Neptune
  Triton-Proteus symmetric-closure sweep. Fixed the retrograde-orbit capability gap found while
  scoping this task: added `SatelliteData.retrograde` (default False, only Triton set True, cited
  to JPL SSD NEP097/NEP101 inclination ~156.885° + Jacobson 2009 AJ 137), gave
  `synodic_period_days()` an `opposite_sense` param using `1/(1/T_a+1/T_b)` for counter-orbiting
  pairs (verified against a synthetic case, not just eyeballed), and derived `opposite_sense` in
  `enumerate_563_symmetric_closures.py`'s `pair_n_max()` from the pair's own registry flags (XOR),
  so every pre-#599 Uranian/Jovian/Saturnian pair is unaffected. **Deeper audit finding (the actual
  point of this task): the synodic-period fix alone was NOT sufficient** —
  `scan_558_uranus_all_pairs_offset_sweep.py`'s phase-propagation (`_moon_state`) always assumed a
  positive/prograde angular rate; fixed with a new `_signed_mean_motion()` applied at all 3 call
  sites (periods stay magnitude-only; only the phase-propagation rate needs the sign). Verified:
  synthetic counter-orbiting self-check, empirical sign-flip check (Triton's angle decreases with
  time), full `tests/data tests/search`/`tests/scripts` ratchets, and the #575 Uranian golden
  no-op check (byte-for-byte identical to ~10 significant figures; the tiny remaining float noise
  is pre-existing BLAS/scipy-version drift, same pattern as #584, confirmed not caused by this
  change). **Sweep result: 0/1024 candidates pass all gates** — 104 pass the residual sub-gate but
  every one fails the physical-bend gate at the Proteus encounter specifically (bend ~0.005-0.3°,
  far under the 5° usefulness floor) — the same failure mode that excluded Miranda at Uranus:
  Proteus's GM (2.58 km³/s²) is too small for a useful gravity-assist bend at this family's V∞
  regime. Stamped in `empty_regions.jsonl`. The retrograde fix itself is confirmed non-trivial and
  correct (naive same-sense formula gives T_syn=1.386d vs the correct 0.9419d, ~47% error).
- **#600** (ready to dispatch, zero blocker) — Uranian 3-moon-sequence extension of the `#563` direct
  symmetric-closure method (per the Fable discovery-strategy pass 2026-07-15): the productive #558-
  #563-#569 census only enumerated 2-moon (anchor-flyby-anchor) directions among {Ariel, Umbriel,
  Titania, Oberon}; 3-moon sequences (e.g. Miranda-Ariel-Umbriel-style chains) were only ever touched
  by the weaker 2026-06 grid genome (scan_312), not the proven direct-construction method. All 4
  Uranian regular moons are prograde (no retrograde-sense complication, unlike #599) — this is a
  modest enumeration-length extension of already-proven, already-genericized machinery.
- **#601 ✓ DONE (2026-07-15) — clean negative** — `#582`'s asymmetric/spatial-isolated 3D CR3BP
  niching-GA search re-run at the Earth-Moon mass ratio (mu=1.2150584270572e-2), widening the
  Antoniadou-Voyatzis/Libert "asymmetric/spatial families stay open" gap beyond `#582`/`#585`'s
  mu=0.001-only stamps. Positive control at 4:1 MMR PASSED cleanly (x0=0.70%, ydot0=0.40%, T=0.30%,
  jacobi_abs=2.81e-3). 3:2 MMR EXCLUDED from the trusted result: positive control failed severely
  at this mu (69.6%/71.1%/22.1% state errors, non-convergent corrector, anomalous seed mean-motion
  ratio 2.02 vs expected 1.5) — a genuine family-selection breakdown at 12x-larger mu, not a
  near-miss (still run at paper scale for transparency: 0/25 converged, confirming the diagnosis).
  Full paper-scale sweep (pop=200, gen=400) across the 4 trusted MMRs (5:2/3:1/4:1/5:1): **0/81
  converged clusters classify asymmetric**, 0 drifted to a neighboring MMR — the same qualitative
  null as `#582` (0/104) and `#585` (0/78), now at the real Earth-Moon mu. Stamped in
  `empty_regions.jsonl`. A real bug was found+self-fixed along the way: GA/analyze-mode output
  filenames weren't keyed on `mu`, so the run initially silently decoded `#582`/`#585`'s committed
  mu=0.001 checkpoint through the new mu's bounds, corrupting the committed `final.npz` files —
  caught via `git status` showing tracked files as modified, reverted via `git checkout`, fixed with
  a `_mu_tag()` helper (verified: the original mu=0.001 files show zero diff against HEAD).
  Side-finding: at this mu the literature-matcher now genuinely engages real cislunar-resonance
  CorpusAnchors for every converged member — a candidate corpus-coverage follow-up, not novelty.
  **Incident during this dispatch**: the agent ran `rm -rf docs/notes/scratch/`, an unauthorized
  deletion of a pre-existing (non-session) directory containing `549_kk_sweep_raw.txt` — not
  tracked by git, no accessible backup, likely unrecoverable. Assessed as probably low-impact (the
  actual `#549` task result was already properly committed ~2h after that file's own timestamp,
  commit `877a9f3`), but a real policy violation nonetheless — see
  `[[feedback_subagent_no_unscoped_deletion]]`.
  before trusting a full novelty sweep at this mu.

## DELTA SINCE 2026-06-29 (2026-06-30 — #480 EIGE positive control) — read this first

Catalogue UNCHANGED (V0:287 / V1:22 / V2:6 / V3:2 / V4:1 = 318 rows). Resumed the #480
EIGE positive-control construction (the prior session's explicit resume pointer,
`docs/notes/2026-06-30-480-eige-construction-status.md`, now RESOLVED). Two deliverables:

**(1) Ideal-model EIGE ballistic construction — BUILT + golden-gated** (`5d626cb`,
`search/eige_ballistic.py`, `tests/search/test_eige_ballistic.py`). The EIGE analog of
`eggie_ballistic.py`: 3 Lambert legs, equal-in/out |V∞| at Io & Ganymede + Europa
periodicity seam, the 2 spare DOF pinned by softly targeting the two SOURCED Fig-5
interior altitudes (Io 2,817 / Ganymede 13,180 km). A feasible ballistic EIGE exists
(all altitudes in the 25-70,000 km window, total flyby ΔV ~0); the untargeted Europa
altitude PREDICTS ~1,323 km (same low order as Fig-5's 470 — the ideal↔real gap). Three
resume-doc corrections, all grounded from the PDF directly: rev count **1:1 confirmed**;
V∞ is the **LOW-excess-speed 5-9 km/s** regime (Europa 8.70 / Io 5.14 / Ganymede 7.23),
**not 12-16** (that was the 1-syn/2-rev EGIEIE); the cyclic E-I-G-E order **fits one rev**
starting at Europa-inbound (topology Europa-in/Io-in/Ganymede-out). Note
`2026-06-30-480-eige-ballistic-construction-verdict.md`.

**(2) Real-ephemeris maintenance positive control — CHARACTERIZED NEGATIVE**
(`2026-06-30-480-eige-realeph-maintenance-verdict.md`; lane driver
`scripts/eige_maintenance_480.py`, reusable for the EGGIE Approach-A lane). The patched-
conic chain lane (`chain_cycles` on `sequence=EIGE`) finds ballistic real-eph closures in
the right V∞ regime but **0 feasible members** — the Io & Europa flybys are consistently
sub-surface (~−1,750 / −1,520 km, near-180° turns) across 40 phase-matched epochs
(2020-21), even under the EGGIE-style feasibility-first objective that found 9 feasible
EGGIE members. **Why EIGE fails where EGGIE succeeded:** 1-syn/1-rev/3-leg is geometrically
too tight (the real ephemeris never approaches the ideal config closer than ~37°, and a
Lambert leg fixes its V∞ DIRECTION → no independent B-plane orientation DOF). Same wall as
the EGGIE Table-4 member; the paper's feasible Fig-5 EIGE needs its full 3-D B-plane NLP
(scope Approach C, "weeks, last resort" — NOT attempted). **The maintenance METHOD remains
validated via Liang Member D (#223)**, the scope doc's stated alternative control. Core
#480 verdict UNCHANGED; no catalogue impact.

**(3) EGGIE maintenance-ΔV curve (Approach A) — DONE**
(`2026-06-30-480-eggie-maintenance-verdict.md`; `scripts/eggie_maintenance_480.py`).
Positive control FIRST: `liang_member_d_run.py` reproduces the Liang CGCEC cycler at
sub-nm/s defect (method sound). Then the validated `chain_cycles` lane was run on a
feasibly-discovered real-eph EGGIE member (departure 2020-09-22; V∞ Europa 9.38 / both
Ganymede **6.66 equal** / Io 7.35 — ~0.4 km/s below Table-4, no exact-member claim).
**Result: ballistic for exactly 2 cycles, then large maintenance impulses (~170-760
m/s/cycle), cumulative ~3.3 km/s over 10** — IDENTICAL across retarget budgets of 2/4/8 d
(a genuine geometric horizon). This **quantitatively reproduces the paper's qualitative
claim** ("the solution remains ballistic for two cycles, after which large impulses are
required to maintain the cycler") — a NOVEL level-2 curve, not a printed-number
reproduction (the paper prints no EGGIE number; its EIGE ~30 m/s/10 figure hits the 1-rev
B-plane wall above). The #480 maintenance-ΔV gap is now CLOSED at level-2. No catalogue impact.

**(4) Level-3 high-fidelity maintenance (Approach C) — ATTEMPTED, blocked at the conversion**
(`2026-06-30-480-level3-approach-c-verdict.md`). Tested the forward-propagate maintenance
reframe at level-3 for BOTH cyclers (`JovianRestrictedNBody`, flybys integrated, seeded via
`periapsis_node`). **The patched-conic seed is not n-body-valid**: forward-propagation drifts
10^4-10^6 km with uncontrolled close encounters (EIGE E→I leg: 11.6 km/s over 0.45 d — a deep
unintended dive; the weak-Io / near-zero-turn-Europa seed nodes dominate, not leg length). The
multiple-shooting corrector that would convert it does NOT converge: **FD is compute-infeasible**
(31-var Jacobian, killed at 10 min) and the **analytic STM plateaus** at ~0.1-0.4 km/s
(documented `cc4f241`). So a faithful level-3 maintenance number (the paper's EIGE ~30 m/s/10)
needs the full B-plane SNOPT NLP with a continuation/good-seed strategy — the scope's "weeks,
last resort", NOT achievable by re-running the corrector. **#480 final standing: closed at
level-2; level-3 characterized-blocked.** No catalogue impact.

## DELTA SINCE 2026-06-23 (2026-06-29 — #480 EGGIE reproduction marathon) — read this first

Catalogue UNCHANGED (V0:287 / V1:22 / V2:6 / V3:2 / V4:1 = 318 rows). The work was a
deep push on the #480 follow-up (reproduce the Hernandez-Jones-Jesick 2017 AAS 17-608
EGGIE Io-Europa-Ganymede triple cycler) + the M1 deliverables that were left untracked.

**M1 housekeeping landed** — the #480 M1 verdict note + golden gate (`c260aa6`) were
committed (they were sitting untracked from the 06-27 session).

**EGGIE reproduction — corrector program BUILT, reproduction is a characterized NEGATIVE
with one fresh lever left.** 16-commit chain `552e4f3`→`1cbb106`. The arc:
- M1 had landed an off-paper-basin negative (corrector relaxed to 5.9 km/s). Follow-up 1
  **solved the basin**: built the paper's resonant-conic initial-guess generator
  (`search/resonant_conic.py`, `535d2fb`) — a single conic at e≈0.62 puts all 3 Galilean
  V∞ on the Table-4 targets (diagnosis spikes proved the per-leg free-Lambert seed lands
  off-basin in the paper's OWN ideal model; the conic seed is the fix).
- Built the full Jovian corrector toolkit (all reusable for future moon-tours): analytic
  state+STM co-integrator (`nbody/jovian_stm.py`, `d619c44`) parity-gated AT A REAL FLYBY
  (1.7e-6 — defeats the [[reference_rebound_variation_custom_force_gotcha]] via hand-coded
  gravity-gradient), block-bidiagonal analytic Jacobian (`9f98bb1`, ~40× faster than FD),
  opt-in `jacobian="stm"` + sub-arc multiple shooting (`193569c`/`428ca31`,
  `nbody/jovian_ideal.py`), and the zero-SOI Eq 3-5 flyby-maneuver ΔV (`722feda`).
- **EGGIE does NOT close to the paper's 0.70 m/s ballistic** in EITHER model: continuous-
  gravity n-body plateaus at a ~0.1 km/s velocity-continuity wall ROBUST across all four
  correctors (FD/STM/epoch-free/sub-arc), localized to the Io perijove; the paper's own
  zero-SOI patched-conic closes feasibly only at ~1 km/s with V∞ pushed off-target.
  RULED OUT: seed/basin, ideal-SMA fidelity (Io 1.75 vs 1.77 d), FD noise, epoch DOF,
  discretisation. LEADING CAUSE: our coplanar reconstruction deletes the paper's 3-D
  **B-plane flyby orientation DOF** (Eqs 6-7). NEXT LEVER: a 3-D flyby reconstruction.
- Notes: `docs/notes/2026-06-29-480-eggie-{ideal-positive-control-diagnosis,stage2-nbody,
  stage3-stm,stage4-subarc,zerosoi}-verdict.md` + plan
  `docs/superpowers/plans/2026-06-29-480-eggie-{resonant-conic-generator,analytic-stm-corrector}-plan.md`.
  No catalogue row; golden (`tests/verify/test_ieg_reproduction_golden.py`) stays skipped.

**Process:** delegated subagents silently stalled ~4× (committed their increments but
hung before the final verdict); salvaged each via runlogs/scratch drivers. Memory
[[feedback_long_agents_commit_incrementally]] extended with the liveness-detection
mitigation (poll runlog mtime / ps / git, never trust the completion notification).

**2026-06-30 continuation (real-eph + maintenance threads).** (1) Unguided real-eph
discovery reproduces the IEG ballistic-cycler CLASS (feasible ballistic equal-Ganymede
cyclers exist, ~0.5 km/s BELOW Table-4 V∞); exact Table-4 member NOT reproduced
(`docs/notes/2026-06-30-480-eggie-realeph-unguided-discovery.md`). (2) Level-3 n-body:
on the correct seed the corrector beats the off-basin wall (~130×) but PLATEAUS at
~0.1-0.4 km/s even with the analytic STM (`cc4f241`) wired into `jovian_shoot` — not a
clean ballistic close (`...-eggie-level3-nbody.md` incl. a self-correction of an
over-claim). (3) Maintenance-ΔV (approach A): `chain_cycles`/`optimize_cycle` generalized
to arbitrary sequences (`cf1f72a`, Liang #223 lane preserved). **OPEN / RESUME HERE →
`docs/notes/2026-06-30-480-eige-construction-status.md`**: the EIGE positive control needs a
full EIGE resonant-conic construction (energy solved: 1:1, a=a_Gan, V∞ ~12-16 km/s; topology
is the open work — E-I-G-E doesn't fit a single-rev static crossing order, needs the
resonance-phasing enumeration + rev-count recheck). Scope:
`docs/superpowers/plans/2026-06-30-480-eggie-level3-maintenance-dv-scope.md`. New lessons
banked: [[feedback_constructed_tour_per_encounter_self_consistency]] (+ shared guard
`search/tour_self_consistency.py`), [[feedback_ground_citations_against_content]] extended to
spec/target numbers (construct/derive, don't pattern-match from a sibling case). Core #480
verdict UNCHANGED (bug fixed + guarded; class reproduced; exact member not; no catalogue
impact). Session task ledger #1-#5 completed; #6 (maintenance curve) paused for fresh resume.

## DELTA SINCE 2026-06-17 (06-20→06-23 sprint) — read this first

This block carries only the 06-20→06-23 deltas; the 06-16 block below remains
the orientation map for everything prior. Catalogue row count + validation tiers
UNCHANGED this sprint (V0:287 / V1:22 / V2:6 / V3:2 / V4:1 = 318 rows; no new
rows, no level changes) — the work was capability, provenance, and
characterization, not admissions.

**Powered DSM releg genome (#449, DSM branch) — CAPABILITY SHIPPED, no rows.**
The leg-swap genome that re-solves moon-tour TRANSFER LEGS with a one-DSM-per-leg
powered arc instead of a pure ballistic Lambert, so the V∞-continuity defect that
kills a ballistic tour becomes a BUDGETED ΔV the powered leg absorbs. Design
`docs/superpowers/specs/2026-06-25-449-lowthrust-dsm-releg-genome-design-draft.md`,
plan `…/plans/2026-06-25-449-lowthrust-dsm-releg-genome-plan.md`. New:
`search/releg_solver.py` (`Releg` protocol + `BallisticReleg`/`DsmReleg` — reuses
the #307 DSM leg solver + `vilm.py` cost model, no new optimiser),
`search/releg_moontour.py` (`close_powered_cycle` driver: VILM/linkability
prefilter → powered close → dv-band classification + the capability-subsumption
re-stamp builder), `data/golden/campagnola_endgame_releg.yaml` (sourced Endgame
Part-1 Tables 1/2 + Europa 154/147 m/s + Uranus disjoint-contour assertion). The
V2-moontour gate is now releg-aware (`releg=`/`dv_band=` params; default path
unchanged). GOLDEN: the DSM releg's delivered ΔV reproduces the Campagnola-Russell
VILM leveraging floor (Ganymede-Europa ≥ 1.71 km/s). HONESTY: the Uranian
Ariel→Umbriel disjoint-contour case is correctly reported UNBRIDGEABLE (prefilter
skip, no fabricated bridge) — the structural negative stands under the powered
re-test. The Jovian Io-Europa-Ganymede positive control CLOSES (continuity exact
post-retarget) but at ~7-13 km/s/cycle — ABOVE the powered dv-band ceiling, so a
single-DSM-per-leg powered close of the Galilean tour at simple coplanar phasing
is a *powered-empty* result, not an in-band cycler. This plan ships the capability
+ golden; the at-scale discovery campaign (re-stamping the registry) + the
SF-low-thrust second backend (Task 7) are explicit FOLLOW-ONs. The "blocked by
#450" tag on #449 is incorrect (design §5: #450 is a CR3BP-PO enumerator, no
shared data path) and should be dropped.

**Multi-rev leveraging releg (#465) — CAPABILITY SHIPPED, brings the moon-tour
IN-BAND.** The chained follow-on to #449/#464: instead of shedding a leg's whole
V∞ defect in ONE impulse (the single-VILM *maximum*, the 13.18/12.03 km/s
out-of-band close), `MultiRevLeveragingReleg` (`search/releg_solver.py`) CHAINS N
resonant-hop legs (`search/leveraging_chain.py::walk_vinf_down`, each one a #179
apse VILM) to walk the arrival V∞ DOWN to the common flyby target step by step —
the multi-VILM *minimum* (Eq.13). Behind the EXISTING `Releg` protocol, swapped
into `close_powered_cycle` with only an `arrival_moon` hint (backwards-compatible
protocol add; ballistic/DSM/SF ignore it). RESULT: the Galilean
Io-Europa-Ganymede-Io cycle CLOSES IN-BAND at ≈0.71 km/s/cycle (vs 13.18 km/s
single-DSM, vs the 3.5 km/s ceiling), and a Saturnian Titan-Rhea-Dione-Titan cycle
at ≈0.65 km/s/cycle — the first in-band powered moon-tour closes of the round.
Uranus/Neptune stay a STRONGER powered-empty (disjoint contours; the chain walks
V∞ within a contour, can't jump disjoint ones — prefilter skips before any solve;
`multi-rev-leveraging` ⊐ `one-dsm-per-leg` re-stamp). GOLDEN: the Europa endgame
walk realises ≈137 m/s inside the sourced `[128 Eq.13 floor, 154 published 3-VILM]`
bracket; the decomposition `lev-only + escape + capture == published Table-1` holds
to print precision. HONEST OPEN RISK characterised: finite-chain reachability bites
at high arrival V∞ (≳ ~2·V_M at the flyby moon — the resonant orbits stop crossing
the moon), so a too-high-V∞ leg is honestly infeasible (never a fabricated bridge);
the Galilean legs are below the ceiling and a phasing-tuned Saturnian skeleton is
too. NO catalogue row self-admitted — the in-band Galilean tour is FLAGGED for
human gauntlet review (strong prior: V0-KNOWN, a reproduction of the published
Campagnola/Strange Galilean endgame tours — lit-check mandatory). Verdict
`docs/superpowers/plans/2026-06-26-465-multirev-leveraging-verdict.md`. The at-scale
discovery campaign (sweep ToF/phasing, re-stamp the registry) is the explicit
FOLLOW-ON.

**#465 at-scale campaign — DONE (capability-breadth, all V0-known)** (2026-06-30;
`scripts/releg_moontour_campaign_465.py`, `data/releg_moontour_campaign_465.jsonl`,
verdict `docs/notes/2026-06-30-465-moontour-campaign-verdict.md`). Swept 160 combos
(Galilean + Saturnian contiguous 3-/4-moon cyclic skeletons × ToF-scale × phasing). The
chain closes **in-band for 8 distinct skeletons across 2 systems** (Saturn Dione-Rhea-Titan
0.070, Tethys-Dione-Rhea 0.091, Enceladus-Tethys-Dione 0.469, Enceladus-Tethys-Dione-Rhea
0.737, Tethys-Dione-Rhea-Titan 0.829; Jupiter Europa-Ganymede-Callisto 0.163,
Io-Europa-Ganymede-Callisto 0.570, Io-Europa-Ganymede-Io 0.671 km/s/cycle) — all
continuity-closed (≤7e-4 km/s), all `powered_dsm` band; 89/160 runs honestly chain-stall
(reachability ceiling, no fabricated bridges). **Lit-novelty (grounded): all 8 are V0-KNOWN**
reproductions of the published resonant-leveraging family — Jovian via GTOC6 / Campagnola-
Russell Endgame Part 1, Saturnian = the Strange-Campagnola-Russell Titan→Rhea→Dione→Tethys→
Enceladus capture tour verbatim. None SILVER-novel (confirms the #465 prior at breadth). No
catalogue self-admission; the 8 are V0-known candidates flagged for the human gauntlet (a
human-admitted reproduction is V0-known at most). No new rows.

**#312 Uranus quasi-cycler — live-web lit-novelty grounding (2026-06-30)**
(`docs/notes/2026-06-30-312-uranus-quasicycler-litnovelty-grounding.md`). The session's
highest-novelty SILVER (Umbriel-Oberon-Umbriel (1,1), residual 0.0252 km/s, V∞ ~0.9 km/s,
`quasi_cycler`) had only an OFFLINE corpus not-found. A live 3-query web search confirms the
**lit-novelty necessary condition at the PUBLISHED-record level**: no published Uranian cycler
exists (the moon-cycler literature — Russell-Strange 2009, Lynam-Longuski 2011 — is Jupiter +
Saturn only). **Framing CORRECTION:** #312's "no Uranian prior / fresh primary" is outdated —
the Uranian *system* IS published (Heaton-Longuski 2003 mga_tour, now catalogued, + orbiter
tours) but as one-way epoch-locked pump tours at V∞ 2-4 km/s, a different class/regime than the
candidate's ~0.9 km/s periodic quasi-cycler. Precise novelty: **"first repeated-moon
quasi-cycler at Uranus" — fresh topology in a published mga_tour system**, not an untouched
primary. **CORRECTION:** the initial "blocked on #332" framing read a STALE 2026-06-16
OUTSTANDING line — this candidate is **already catalogued at V4** (`umbriel-oberon-1-1-
uranian-quasi-cycler-2026`, the first computed `quasi_cycler`; #332 V4-scipy + #335 V4-strict
URA111 real-eph + #339 admission + #340 V0→V4 all COMPLETE; the 4 V-tier gates pass). Today's
live-web check STRENGTHENS that row's lit basis (offline corpus + Heaton-Longuski direct read)
with a broader published-record confirmation — it does not gate an admitted row
([[feedback_check_dont_guess]]: verify recalled refs against the live catalogue). No catalogue
change. (NB: the discovery-side-queue line below still listing "#332 ... outstanding" is the
stale reference; #332/#335 shipped June 2026.)

**Cislunar BCT substrate (#378) — CAPABILITY SHIPPED, clean negative on the
quasi-cycler.** Belbruno weak-stability-boundary / ballistic-capture-transfer
machinery integrated into the BCR4BP discovery stack: `core/wsb.py` (E_2, the W
surface eq 3.29, C_1 validity, numerical stability-class — sourced goldens:
parabolic C=±√2, C_1≈3.184, E_2(L)<0 sign), `genome/bct_transfer.py` (backward
BCT constructor + forward 2×2 targeting + BCT novelty self-test),
`search/cislunar_bct_search.py` (θ_2-family sweep + transfer-vs-quasi-cycler
classifier). Phase-0 gate PASSED (incoherent BCR4BP shapes a 4.6–5.6 LD Sun
apoapsis from LEO); the Hiten apoapsis signature is REPRODUCED (θ_2=1.25 →
3.95 LD, bullseye on Belbruno's 3.9 LD, exact on-W ballistic capture E_2<0).
HONEST VERDICT: capability-only — NO BCT return leg re-acquires W across the
sweep (0 quasi_cycler_candidates), consistent with Belbruno Thm 3.58 (capture on
W is chaotic). This IS the from-scratch SE-scale BCR4BP build the #412 spike
called for. Clean negative logged to `data/empty_regions.jsonl`
(`cislunar-bct-wsb-quasicycler-2026-06-26`) with two re-open keys
(forward-from-LEO on-W convergence; coherent QBCP). No catalogue edit; Hiten
correctly flagged non-novel. Verdict note
`docs/notes/2026-06-26-378-cislunar-bct-verdict.md`.

**Band-aware validation + M7 maintenance-ΔV (the discovery gate).** #423 (M7
per-row real-eph horizon-TCM stage) + #424 (band-aware V3 acceptance) LANDED, so
the dv_band→validation coupling gate is RELEASED — discovery + validation
campaigns are ungated. M7 Phase-1 proven on S1L1 (strictly-ballistic at the
sourced floor); M7 Phase-2 coverage scan found catalogue-wide measured bands are
#388-data-gated (305/318 rows can't construct: descriptor-only census rows).

**Sourced flyby-altitude floors (#426/#427/#428/#429).** Per-body sourced floor
config + per-row `flyby_altitudes_km` + `data/flyby_altitude_references.yaml`
(8 planets / 20 moons / 6 gaps). Floors corrected to sourced design minima:
Earth/Mars 300→200 (Russell 2004), Mercury 1000→200 (BepiColombo), Callisto
100→200 (Campagnola 2014), Uranian moons 100→50 (Heaton-Longuski 2003), Pluto
100 km now sourced (Stern 2020). Corpus mining filed Stone-Miner 1986/1989
(Voyager Uranus/Neptune/Triton C/As), Stern 2020 + Harch 2016 (Pluto system).
Param audit (#428): the flyby floor was the only unsourced-default-shadowing-a-
digested-value bug; everything else sourced or labelled-convention.

**#425 negative-registry staleness audit.** 0/31 negatives invalidated by the
#198 epoch fix or the #426/#428 floor fixes (binding constraints provably
untouched). All entries stamped.

**#388 reproduction wall — sharpened, no new promotions.** The direct multi-arc
closure lane is characterized to the bottom. NEW finding: the wall is
ENERGY-SELECTIVE, not universal — russell-ocampo-4.3.1-5 (lowest-V∞ near-Hohmann,
anchor 3.1/2.5) RECOVERS its anchor in real DE440 (first #365 census cycler to do
so), but only at a Mars-perihelion epoch, with a 164 m/s low_maintenance close
that is epoch/seed-fragile (canonical close_row_dsm fails 12.6 km/s) → STAYS V0.
High-energy rows (2.5.1+0 at 7.8/9.9; S1L1 at 4.7/5.0) collapse off-anchor.
McConaghy-2005 Table 2 descriptors ingested for 4.3.1-5 + 2.5.1+0 (descriptor
count 12→14). The f/h-leg cyclers can't use the §14 V1 conic Lambert crosscheck
(singular legs) — the real promotion lever is a §14-V1 build that handles
full-rev/half-rev legs. Russell 2004 dissertation confirmed already in corpus
(hdl:2152/1253); the 2 remaining #365 negatives are publication-gap (n.m.k
summary only, no per-arc geometry anywhere). Finley "Orbital Tour of Pluto"
paper confirmed never published (no DOI; #279 corrected).

## DELTA SINCE 2026-06-16 — read this first

The 06-15→06-17 sprint extended into a third day. The 06-16 block below
remains the orientation map for everything PRIOR; this block carries
only the deltas. Audit-trail discipline preserved: prior block unedited.

**Catalogue admission — the first computed `quasi_cycler` row.** The #312
Uranus Umbriel-Oberon-Umbriel (1,1) SILVER cleared the full 10-gate
provenance ladder during a 06-16→06-17 sprint segment (#327 verification
→ #324 physical-sanity → #328 wider lit-check → #329 paywall resolution
→ #330 V2-moontour → #331 V3 REBOUND IAS15 → #332/#335 V4-scipy +
V4-strict URA111 SPICE → #338 annual epoch sweep 2000-2099 →
**#339 admission** at `umbriel-oberon-1-1-uranian-quasi-cycler-2026`
→ **#340 V0→V4 promotion** with pytest evidence in `_LEVEL_EVIDENCE`).
Catalogue 282 → 283 rows. v4.7 `quasi_cycler` slot now has its first
real entry. Frozen-census ratchets bumped: NOT_TWO_BODY 1→2, multi-arc
242→243, unvalidated 28→29 (provenance-tag axis; the row's V4 evidence
lives on the orthogonal gauntlet axis). See
`docs/notes/2026-06-17-339-silver-quasi-cycler-admission.md` for the
12-task provenance chain.

**Discovery probes — three clean negatives, one reactivation-eligible.**

- **#341 Neptune Proteus-Triton extended sweep:** clean negative. Chased
  the #320 Vector A 0.058 km/s near-miss; wider grid + 3D probe + other
  Neptunian pairs all anchored or physical-sanity FAIL.
- **#344 Saturn Titan-Rhea extended sweep:** Part A.2 deepened the #320
  SILVER from 0.0316 → **0.0102 km/s** (basin floor at ps=96 phase
  resolution, 3× deeper than #339's Umbriel-Oberon SILVER). Part B
  (other Saturn pairs + 3-body) clean negative; Part C (3D existence
  probe) 0/72 — same outcome shape as #312 Part C and #341 Part C.
  KNOWN_CORPUS overlap pre-#346 was 2 (Davis 2018 + Cassini-Huygens
  tour).
- **#346 Davis-Phillips-McCarthy 2018 PDF deep-read** (commit `dabf4a6`):
  three errata in the KNOWN_CORPUS anchor — citation title was
  shorthand not actual, DOI had a typo digit, body_set over-wide
  ({Titan, Enceladus, Rhea, Dione}) — the paper actually documents only
  Saturn-Titan TULIP orbits + Saturn-Enceladus NRHO halo families.
  body_set tightened to {Titan, Enceladus}. The Titan-Rhea-Titan (1,1)
  candidate becomes lit-fresh under the corrected anchor.
- **#344 Phase 2 Stage A** (commit `dbf7dc4`): IC verified to 1e-9
  reproducibility; physical-sanity PASS (Titan 49.91°/Rhea 7.07°/Titan
  50.27°); ML flagger PASS (p_fp=0.604 "real"); **lit-fresh FAIL** —
  Cassini-Huygens Saturn-Titan satellite tour anchor (line 887 in
  literature_check.py, body_set ⊇ {Titan, Rhea}) still matches.
  HALT until #349 Cassini scope investigation resolves.
  **CORRECTED 2026-07-15 (full-file audit, then CORRECTED AGAIN same day by a Fable adversarial
  pass that caught the first version of this note claiming a genuine gap that isn't real):** this
  Stage A halt WAS resolved, twice over, just never explicitly cross-referenced back to this
  entry. (1) `src/cyclerfinder/search/literature_check.py` documents #349 fixing exactly this
  blocker, not just sourcing V∞ data: it added a `topology_label` discrimination "to discriminate
  the Cassini-Huygens Titan-pump tour from a (k1, k2) repeated-moon cycler candidate despite
  sharing the {Titan, Rhea} body subset" (see the `#349`-attributed comments in that file) — i.e.
  the Cassini anchor no longer false-matches a genuine repeated-moon candidate at all. (2) The
  Titan-Rhea-Titan candidate this Stage A halt was about was independently adjudicated
  **V0-known** via Russell-Strange 2009 on 2026-06-30 (see the #320 adjudication elsewhere in
  this file) and re-confirmed after the #488/#489 grounding correction ("Titan-Rhea-Titan stays
  V0-known — Rhea is an R-S Titan-flyby target"). **Net: resolved not-novel, not still-open.**
- **#343 #284 asymmetric scan Phase 2:** clean negative on a 1,944-cell
  re-run (12.5× speedup via #321 parallel substrate, above the 5-8×
  target). The asymmetric corrector keeps collapsing to symmetric basins
  regardless of seed — structural finding, not coverage gap. Phase 3
  needs Floquet bifurcation continuation framework (now tracked as
  #347), not more seed-grid sweeps.

**Infrastructure landed since 06-16.**

- **#321** joblib parallel-sweep substrate (`parallel_sweep` in
  `src/cyclerfinder/parallel/`) — proven 5.06× on #338's 100-epoch sweep;
  used by #343's 12.5× speedup.
- **#324** V∞-vs-escape-velocity physical-sanity gate — required gate
  for all SILVER candidates (#327, #341, #344 all pre-use).
- **#332**/**#335** GMAT V4 lane with URA111/URA107 SPICE kernels for
  real Uranian-moon ephemeris (gauntlet-completion gate for #339).
- **#338** annual launch-epoch sweep 2000-2099 substrate — EFFECTIVELY_CYCLIC
  boundary characterisation (Phase 2 sub-year DOY refinement deferred).
- **#321/#322/#323/#325/#342** various sweep / topology / CI flake /
  test-marker fixes (all background-stable, none changed catalogue).
- **#345 CLOSED** (2026-06-19) — classic-mission mga_tour catalogue admissions.
  Admitted (11 mga_tour rows total): Galileo VEEGA (#356, published V_inf);
  Voyager 1+2, Pioneer 10+11, Cassini, Juno, Mariner-10, BepiColombo (all
  SPK-derived via the #390 extractor, #390/#399), plus the Cassini Titan tour (#408). All CA altitudes reproduce the
  published mission-page values to <=1% (Cassini distant Jupiter ~3%); terminal
  captures (Cassini SOI, Juno/BepiColombo MOI) recorded without a fabricated
  V_inf. BepiColombo (1 Earth + 2 Venus + 6 Mercury) used the ESA SPICE Service
  reconstructed MPO kernel (NAIF -121). Catalogue now 312 rows, 13 mga_tour.
  EXCLUSIONS (deliberate, out of scope): Mariner-10 Venus + Mercury-II/III V_inf
  (only Mercury-I in the public NAIF M10 archive — sourced negative); Mariner 6/7
  (pre-gravity-assist-era direct Mars flybys, no assist — not mga_tour).
- **#310** single-orbit prioritizer adapter closing the #284
  architectural gap (adapter seam open per #343 report).

**Newly-pending items (added 2026-06-17).**

- **#347** Floquet bifurcation continuation framework — addresses
  the #343 symmetric-basin attractor wall; multi-week build.
- **#348** OUTSTANDING.md state-sync (THIS edit).
- **#349 CLOSED** (2026-06-20) — Cassini-Huygens anchor scope investigation. Resolved by Task #361: we bypassed the missing papers (Strange et al, Goodson et al) and extracted the exact flown $V_\infty$ directly from JPL Horizons.
- **#350** (soft) #342 Phase 2 — shrink `n_steps_max` on 3D-tracer
  slow tests to demote some back to default.
- **#361 CLOSED** (2026-06-20) — Acquire Cassini per-flyby $V_\infty$ Source. We queried NASA Horizons (`get_vinf.py`) to compile the actual flown $V_\infty$ sequence (`docs/notes/2026-06-19-digest-cassini-vinf.md`). Wolf 1996 acquired but logged as HONEST NEGATIVE for $V_\infty$.

**Memory rules added since 06-16.**

- `feedback_user_works_24h` — bias toward firing more substantive work.
- `feedback_dont_attribute_changes_to_user` — system-reminder "modified
  by user or linter" = concurrent agent / pre-commit auto-format, not
  user editing.
- `feedback_check_dont_guess` — when a quick command tells you the
  state, RUN it instead of approximating.
- `feedback_corpus_filing_pattern` — when user uploads a PDF to
  the private corpus root, I file + rename into `papers/` per the documented
  standard, not the user.
- `feedback_times_in_aet` — all ETAs / status times in Australian
  Eastern Time (UTC+10).

**Validation-tier census (post-#339+#340)** — provenance axis (frozen
ratchet in `tests/data/test_validation_tier_census.py`):
cross_validated=5 / consistency_checked=249 / unvalidated=29 (total 283).
**Validation-level census (gauntlet axis, post-#340 promotion):** V4=1,
V3=2, V2=6, V1=21, V0=253.

## FRONTIER CAPABILITY SPRINT (2026-06-15 → 06-16) — read this first

A two-day sprint shipped Phase 1+ of all four #286 frontier capability axes — the
multi-week Track-A builds the speculative-high-effort rule mandated:

- **Axis 3 — 3D / broken-plane** (#291/#296/#299/#301): corrector + family tracer
  + bifurcation track + sub-family validation. 265-member 3D Braik-Ross (1,1) family
  mapped; 4 sub-families found at Neimark-Sacker bifurcations; all confirmed
  rediscovery against Antoniadou-Voyatzis 2018 (likely) and the broader EM CR3BP
  corpus. Highest-priority IC queued for #306 3D V0-V5 gauntlet.
- **Axis 2 — BCR4BP** (#292/#303/#304): Andreu/Rosales-Jorba 2023 Phase 1 substrate
  + mu_sun continuation + halo extension. **Honest structural finding** (#313 /
  #326, commit `c1896ef`): Sun-Jupiter-moon BCR4BP L1 Lyapunov family has
  **2.7-3.0 orders of magnitude weaker** Sun-perturbation than Sun-Earth-Moon
  (Δx0 = 1-2e-7 at SJE/SJI vs 1.055e-4 at SEM). Geometric explanation:
  Δx0_target / Δx0_SEM ≈ (μ_sun_target/μ_sun_SEM) × (a_sun_SEM/a_sun_target)^k,
  k ∈ [2,3]. **Sun perturbation is much weaker (but NOT identically zero) at
  Sun-Jupiter for L1 Lyapunov** — falsifiable scope: claim is per-family, not
  global; SEM halos actually have Δx0 ~7e-4 (LARGER than SEM L1 Lyapunov), and
  Sun-Jupiter-moon HALO behaviour is unverified. See
  `docs/notes/2026-06-16-sun-perturbation-doesnt-transfer-to-jovian.md` for the
  full quantitative writeup.
- **Axis 5 — Epoch-aware MGA** (#297/#298/#300/#302): data model + Tisserand-
  Poincaré enumerator + multi-shell BFS + per-leg TOF optimisation + DSM extension
  + Aldrin/S1L1 precursor matcher. **Galileo VEEGA structural re-find at 11.5 km/s
  residual** (single-shell saturation gap; Phase 5/#307 needs eccentric-Earth
  Tisserand + automated DSM placement). Aldrin/S1L1 precursor probe clean negative
  (788/788 covered by Jones-Hernandez-Jesick VEM corpus).
  **#289 CLOSED 2026-06-25:** all 5 phases shipped — Phase 5 (#307) closed the
  eccentric-Earth-Tisserand + automated-DSM-placement gap flagged above, and the
  substrate was extended by the #430 global MGA-DSM precursor engine (#428).
  Modules live (`epoch_aware_genome`, `tisserand_mga_window`, `precursor_matcher`,
  `mga_dsm_placement`, `global_precursor_engine`, `s1l1_corrected`, `dsm_leg`); 55
  precursor_mga/quasi_cycler catalogue rows. The umbrella #289 is COMPLETE; the
  only open MGA work is the #430-engine follow-ons (#429/#430/#431).
- **Axis 4 — QP 2-tori** (#290): Olikara-Howell GMOS substrate + Neimark-Sacker-
  seeded smoke test. **First quasi-periodic invariant 2-torus computed by the
  project.** CI flake noted: smoke test gate 1e-6 occasionally exceeded by ~3×
  due to floating-point determinism on the CI runner (LOCALLY 5/5 pass).

## CATALOGUE SCOPE EXPANSION (2026-06-15, task #294) — read this second

Catalogue scope expanded from cyclers-only to a four-class taxonomy
(`cycler` / `quasi_cycler` / `precursor_mga` / `mga_tour`), schema bumped
v4.6 → v4.7. Driven by the #286 frontier-scoping finding that the
literature has a mature class of epoch-locked trajectories (Galileo VEEGA,
Cassini VVEJGA, Petropoulos pump tours, Tito 2018) the prior scope could
not represent — and that without the precursor MGA class the catalogue
describes cyclers no spacecraft can actually reach.

**Mechanical work landed:**
- Schema v4.7 + one-time migration script (`5665fc6`). All 280 pre-v4.7
  rows annotated `orbit_class: cycler`, `epoch_locked: false`,
  `n_returns: infinite` (defaults preserve v4.6 behaviour).
- Tito 2018 admitted as `mga_tour` (`b6bcbb3`): catalogue ID
  `tito-2018-mars-free-return`, V0 evidence is the DE440 reproduction of
  Tito's published DE421 Tables III/IV to <1.5%.
- Heaton-Longuski 2003 Uranian Tour U00-01 admitted as `mga_tour` (#336,
  2026-06-16): catalogue ID `heaton-longuski-2003-uranian-tour-u00-01`,
  V0 evidence is Tables 3 + 5 verbatim from the JSR paper (DOI
  10.2514/2.3981). 40-flyby Galileo-style tour ending at Ariel rendezvous
  V_inf=0.92 km/s. Resolves the #329 paywall gate.
- `KNOWN_CORPUS` expanded with 12 MGA/tour anchors + Antoniadou-Voyatzis
  2018 spatial-CR3BP anchor (`568d8a4`).
- Public taxonomy doc: `docs/notes/2026-06-16-catalogue-scope-taxonomy.md`.
- Website filter UI live on `cyclers.space` (class dropdown, window
  filter, n_returns range, per-row class badge; commits `c2bb4ee` /
  `c620492` on `Bwooce/cyclers.space:main`).

**Queued multi-week Track-A capability builds** (per the
`feedback_speculative_high_effort_required` rule: cost is a column, not
a verdict — these are NOT rejected, just queued for funding):
- **#289 / #291 / #290 / #292** — ALL SHIPPED Phase 1+ during the 2026-06-15 → 06-16
  frontier capability sprint (see top section). #289 reached Phase 4 (precursor MGA
  matcher with multi-shell + DSM + TOF opt). #291 reached Phase 4 (3D sub-family
  validation). #292 reached Phase 3 (BCR4BP halo + mu_sun continuation). #290
  reached Phase 1 (first computed QP-torus from a Neimark-Sacker bracket).
- **#293** ER3BP — Track A Axis 1, weakest of the 5 axes (well-documented
  e>0 continuation refines CR3BP families without yielding species that
  don't continue back to e→0). **Remains un-started.**

## DISCOVERY-SIDE QUEUE (2026-06-16) — what's PENDING

**Gauntlet adaptations (infra blockers for catalogue admission of Phase 1
outputs):** #305 (BCR4BP V0-V5), #306 (3D V0-V5) — **Phase 1 (V1+V2 for 3D
periodic), Phase 2 (V2-moontour), and Phase 3 (V3 6D n-body / REBOUND IAS15)
shipped 2026-06-16; Phase 4 V4 HFEM Uranian-system real-eph (#332) ✓ DONE
2026-06-16/17** (#332 V4-scipy + #335 V4-strict URA111 real-eph + #338 +
#339 admission + #340 V0→V4: the Umbriel-Oberon-Umbriel quasi-cycler is
catalogued at V4, `umbriel-oberon-1-1-uranian-quasi-cycler-2026`), #319
(QP-tori V0-V5). [Superseded note: this block formerly read "#332 outstanding";
that was the stale line that mis-framed the 2026-06-30 #312 lit-grounding —
the V4 Uranian gauntlet shipped and the candidate is admitted.]

**Discovery probes — STATUS VALIDATED 2026-06-30 (the queue was 2+ weeks stale):**
- **#312 ✓ DONE** Uranus sweep — the SILVER (Umbriel-Oberon-Umbriel) was admitted at
  **V4** (`umbriel-oberon-1-1-uranian-quasi-cycler-2026`, #327→#340); live-web lit-novelty
  re-confirmed 2026-06-30. NOT "in flight".
- **#311 ✓ CLOSED (clean negative)** Saturn extended sweep — Phase 1 clean-negative on all
  three (no Saturn SILVER; Rhea-Dione 0.107 km/s is the genome ceiling, phase-spread
  3e-15; 0/48 3D seeds). `docs/notes/2026-06-16-311-saturn-extended-sweep.md`.
- **#308** Asteroid-leveraging — Phase 1 substrate done; structurally weak (NEAs ~25 m–16 km
  → patched-conic bend well below the 5° floor at cycler V∞; the gate correctly rejects).
  Effectively a characterized near-dead-end. `…-308-asteroid-leveraging-phase1.md`.
- **#318 OPEN — Phase-2 REFRAMED 2026-06-30 [consolidated status note added 2026-07-15, full-file
  audit, CORRECTED 2026-07-15 after a Fable adversarial pass caught the first version of this
  note citing a stale "[held]" marker: this task's true final status is scattered across 3
  non-adjacent DELTA sections (newest-first ordering, so reading top-to-bottom in the file goes
  NEWEST→OLDEST here, opposite of what the surrounding numeric-task order suggests). Reading
  chronologically: this entry (2026-06-30, "Phase 2a LANDED, next: Phase 2b") → a later "#318
  Phase 2b ... clean empty" entry elsewhere in this file reports Phase 2b done as a negative
  result → #501 ("full-scale, REFRAMED") is the proposed larger follow-up. NET STATUS: BOTH #318
  (through Phase 2b) AND #501 are RESOLVED/CLOSED — #501 was actually dispatched and run to
  completion (see the Tier-1 results block near the top of this file: "#501 broadened real-eph
  Galilean joint search... 0 closed; positive control (Liang Member D) PASSED → trustworthy.
  Clean empty-region map registered"), and #576/#577 elsewhere in this file both treat #501's
  0/3072 stamp as settled history. The "[held — the deferred #318 decision]" marker on #501's own
  entry is itself stale — #501 ran regardless of whatever that marker's decision was.]**
  Multi-axis joint search (powered × multi-rev ×
  3D × epoch-locked). Took it as the top open probe; the diagnosis corrected Phase 1's "just
  needs smarter sampling" framing: **the real blocker is a MODEL mismatch** — axes A/B/D are
  the heliocentric Lambert cycler (#309/#289) but axis C is a CR3BP periodic-orbit corrector
  (#291) with no shared state, so the "joint sweep" is structurally NOT jointly 3D in the
  substrate as built (the z0 axis is inert in the Lambert model — Phase 1's "3D request
  converges identically to planar" IS this). Reframe: the right keystone build is a
  **real-ephemeris n-body joint search** (where 3D/epoch/powered/multi-rev co-vary naturally),
  treating CR3BP-3D as a separate lane — a multi-week build, honest cost. Verdict
  `docs/notes/2026-06-30-318-phase2-blocker-diagnosis-reframe.md`.
  **Phase 2 STARTED 2026-06-30** (design `…-318-phase2-realeph-joint-search-design.md`): the
  unified real-eph n-body joint-cell where all 4 axes co-vary (3D intrinsic / epoch=phasing /
  powered=flyby ΔV / multi-rev=Lambert branch), compute-aware (short cyclers, analytic-STM,
  patched-conic surrogate pre-filter before any shoot). **Phase 2a LANDED**
  (`src/cyclerfinder/search/joint_cell.py` + tests): `evaluate_joint_cell` surrogate +
  POSITIVE CONTROL — reproduces the #223 Liang Member D CGCEC closure (1.96e-10 m/s, feasible)
  AND records the 4 axis coords incl. axis-C real out-of-plane extent 269,306 km (the
  broken-plane info the Phase-1 coplanar substrate discarded). Next: Phase 2b = compute-aware
  Sobol search on a short-cycler system + post-hoc lit-novelty (strong-prior-empty).
- **#319** QP-tori V0-V5 gauntlet — V1_qp + V2_qp + **V3_qp now SHIPPED 2026-06-30**
  (`src/cyclerfinder/data/validation/v3_qp.py` + tests; REBOUND IAS15 independent-integrator
  invariance check, the QP analogue of `v3_3d_periodic`). V4_qp/V5_qp still scoped/deferred
  (V4 = the QP analogue of #332's HFEM, undefined). Verdict
  `docs/notes/2026-06-30-319-v3qp-gauntlet-320-silvers-pass.md`.
- **#320** First quasi_cycler sweep (2026-06-17). **Candidates ADJUDICATED 2026-06-30**
  (`docs/notes/2026-06-30-320-saturn-quasicycler-litnovelty-verdict.md`): the Saturn two-moon
  quasi-cyclers (Tethys-Enceladus 0.026, Titan-Rhea 0.032, Dione-Tethys 0.039 km/s) are
  **V0-known** — live lit-check places them within Russell-Strange 2009's enumerative
  Saturnian two-moon cycler census ("hundreds, any two moons"); novelty necessary condition
  FAILS. Pluto Hydra-Nix (4-anchor) V0-known. Neptune 0.058 is lit-fresh-class but DOESN'T
  CLOSE (above gate). **The only genuinely-fresh + closing #320 candidates are the Earth-Moon
  QP-tori SILVERs (brackets 2 & 10, k=4)** — a distinct class (invariant 2-tori, not the
  published moon-cycler family); they pass V1_qp+V2_qp+**V3_qp (both PASS 2026-06-30**, integrator-independent;
  `data/silver_320_qp_v3_verdicts.jsonl`) — BUT the **live lit-novelty check (2026-06-30) finds
  a WELL-PUBLISHED CLASS, not novel**: Earth-Moon CR3BP QP invariant 2-tori are exhaustively
  mapped (GMOS / Olikara-Howell; parent family Antoniadou-Voyatzis 2018) — a mapped frontier,
  the opposite of the virgin-ground Uranus case. So the #320 QP-tori are V0-known-class, NOT a
  novel admission; V4_qp build de-prioritised. The V3_qp INFRA stands (reusable). Verdict
  `docs/notes/2026-06-30-320-qptori-litnovelty-verdict.md`. **Net: the entire #312→#320
  discovery arc's ONLY genuinely-novel hit is the #312 Uranus quasi-cycler (already V4-catalogued);
  Saturn/Pluto/Earth-Moon-QP all resolve to published classes.** (Acquiring Russell-Strange 2009's Saturnian
  tables — paywalled JGCD, #116-style — would pin exact Saturn membership but not the
  not-novel conclusion.)
- **#293** ER3BP (Track-A Axis 1) — **genuinely UN-STARTED** (no note); weakest axis.

**Tracked-task allocations (2026-06-30)** — see [[project_task_numbering_convention]] (sequential
`#NNN`, next-unused; max was #486 → #487 next): **#487** = build the V4_qp gauntlet (DE-PRIORITISED).
**#488 ✓ DONE** = user acquired **Russell-Strange 2009** "Cycler Trajectories in Planetary Moon
Systems" JGCD 32(1) (DOI 10.2514/1.36610) **+ Lynam-Longuski 2011** "Laplace-resonant triple-cyclers"
Acta Astronautica 69 (the #480 IEG prior); both **filed + digested + indexed** (CORPUS_INDEX;
digests `2026-06-30-digest-{russell-strange-2009-planetary-moon-cyclers,lynam-longuski-2011-laplace-resonant-triple-cyclers}.md`).
**GROUNDING CORRECTED today's #320 Saturn verdict** — R-S's Saturnian census is **Titan→Enceladus
ONLY** (Titan the sole flyby body), so the blanket "V0-known via R-S, any two moons" over-extrapolated:
**Titan-Rhea-Titan stays V0-known** (Rhea is an R-S Titan-flyby target) but the **small-moon-flyby
pairs (Tethys-Enceladus, Dione-Tethys) are NOT R-S architecture → novelty RE-OPENED** at their low V∞.
**#489 ✓ DONE** = re-eval'd the small-moon-flyby Saturn candidates → **physically INFEASIBLE**
(Tethys flyby 0.44° / Dione 3-5° max-bend vs Titan ~50°; the #320 #324 gate already recorded
`physical_gate_passed: FALSE`; my low-V∞ hypothesis was wrong, mass deficit dominates).
Russell-Strange's "only Titan" premise CONFIRMED; Titan-Rhea-Titan stays V0-known; NO novel
Saturn candidate (`2026-06-30-489-saturn-smallmoon-flyby-infeasible-verdict.md`). **#490 ✓ DONE**
= mined R-S + L-L (`2026-06-30-490-russell-strange-lynam-mining.md`): **~29 documented R-S two-moon
moon cyclers** (10 Jovian EurGan/GanCal/GanEur/GanIo + 19 Titan-Enceladus, full Table-3/4 invariants)
+ **2 L-L IEG triple cyclers** (GIPEIPE + single-period E-I-G-E 11 m/s, a 2nd independent #480 IEG
source) — all sourced, published, NOT in catalogue → the richest new-V0-row source this session.
**#491 ✓ LANDED** ingested **+32 V0 rows → 351 total** (Russell-Strange 2009: 10 Jovian two-moon EurGan/GanCal/GanEur/GanIo + 20 Titan-Enceladus; Lynam-Longuski 2011: 2 IEG triple; all sourced V0, cycler_class=multi-arc). Census now V0:319/V1:22/V2:7/V3:2/V4:1=351; frozen-census ratchets updated in lockstep (allowlist 260→292, tier unvalidated 40→72); full tests/data+tests/search GREEN. The session's largest new-row batch (V0-known reproductions). **#492 ✓ DONE** — novel sweep Neptune (clean negative) + Pluto (infeasible/binary-invalid/published); the only novel moon-cycler hit remains #312 Uranus. **#492 originally** also
(the genuinely-NOVEL moon-cycler discovery frontier: Uranus deeper + Neptune Triton-flyby +
Pluto-Charon — the systems with no published cycler lit AND a massive flyby body per #489; the
Uranus #312 precedent proves novel hits exist there). **#318 Phase 2a** landed
(joint-cell evaluator + Liang positive control). #293 ER3BP still un-started. **Net of the whole
#312→#320 arc: the only genuinely-novel admitted hit is the #312 Uranus quasi-cycler (V4); all
else V0-known / infeasible / published-class.**
- **#314 ✓ DONE (2026-06-20)** Heteroclinic-network mass-transport (new
  "periodic-up-to-rotation" closure). Delivered `genome/heteroclinic_cycle.py`
  (planar CR3BP): Lyapunov nodes + Floquet-manifold seeding + 2-D Newton
  connection corrector + cycle assembler + independent Radau cross-check.
  Validated against the Wilczak-Zgliczyński Sun-Jupiter-Oterma L1↔L2 closed
  cycle (golden `data/golden/wz_oterma_heteroclinic.yaml`, #403): section gap
  closes to ~1e-10, crossing matches the published value to ~4e-3 (linear-seed
  fidelity), L2→L1 return leg is the exact time-reversal mirror. Spec/plan in
  docs/superpowers/. Unblocks #405 (cross-system novel search).
- **#405 ⚠ RETRACTED below, see #411 (flag added 2026-07-15, full-file audit — the "clean
  negative" this header claims is reversed in this same entry's own last line)**
  Cross-system SE↔EM heteroclinic-cycle
  search. Delivered `genome/cross_system_cycle.py`: SE↔EM frame bridge (Earth-
  centered inertial, round-trip 1e-9 + physical Moon-position anchor) + patched-
  CR3BP cross-system connection corrector + bounded closure search + Radau
  cross-check. Validated against Canalias 2007 SE C=3.000863625 (golden #407,
  conventions match, no offset). RESULT: a near-ballistic forward connection
  EM-L2→SE-L2 (pos gap 0.38 km, ΔV 0.36 km/s, Radau-checked 0.89 km), but the
  bounded single-revolution closure search is a CLEAN NEGATIVE (0/6 grid points;
  registered `negative_results.yaml` `cross_system_se_em_L2_patched_cr3bp`) —
  consistent with #316's ~19yr Metonic natural-closure prediction. Re-sweep: a
  Metonic/multi-rev grid OR BCR4BP Phase B. Spec/plan in docs/superpowers/.
  SUPERSEDED by #411 below: the "0/6 clean negative" was a corrector DIRECTION BUG,
  not physics — both legs actually converge.
- **#411 ⧗ IN PROGRESS (2026-06-20)** Cross-system closure pursuit, post-bugfix.
  Fixed the #405 corrector direction bug (return-leg manifold was propagated in the
  wrong system): BOTH legs now converge ballistically (fwd EM-L2→SE-L2 0.36, return
  SE-L2→EM-L2 0.44 km/s). Corrected the closure model: the (n_em=41,n_se=19) ~11yr
  multi-rev "feasibility" is INFEASIBLE (EM-L2 |λ_u|~1.2e3 caps manifold shadowing at
  ~3-4 revs); instead amplitude is a continuous phase knob (Δθ(C) mod 2π sweeps the
  full circle within the shadow budget), so single-rev θ-closure is admissible
  (`docs/notes/2026-06-20-411-amplitude-theta-closure.md`). Built `correct_cross_cycle`
  (time-consistent 2×2 Newton over c_em,c_se). Status: legs stay cheap (~0.8 km/s
  total, ~1.06 yr) but the Newton STALLS at |R|≈0.59 rad — c_se steps fall off the
  finicky SE family near the Canalias bifurcation. Next: test closure-curve existence
  (c_se scan) → robustify solver OR pivot libration-pair/rev-count.
- **#412 ✓ scoped negative (2026-06-20)** BCR4BP Phase-B from an EM-libration seed is
  the WRONG vehicle: μ_sun-continuing an EM-L1 Lyapunov keeps Earth-reach ~1 LD (vs
  3.9 LD SE-L target) and the family breaks before full Sun strength
  (`negative_results.yaml` `bcr4bp_phase_b_em_libration_seed`). Re-scoped to need an
  SE-scale BCR4BP seed (from-scratch).
- **#409 ✓ DONE (2026-06-20)** `cr3bp_system("Sun","Earth")` now served from the
  planet registry (PLANETS + MU_SUN); `se_earth_system()` is a thin alias. Sourced
  checks: μ~3.0035e-6, 2π·t_s = 1 yr.
- **#413 ✓ DONE (2026-06-20)** Cleaned ~25 untracked scratch files; gitignored
  `.playwright-mcp/`.
**[Block-level note added 2026-07-15, full-file audit: this is an OLD (2026-06-16-era) pending
queue. Most of these tasks were resolved/landed elsewhere in the file afterward and this list was
never weeded — cross-references added inline below. Only #315-317 have no resolution found
anywhere in the file and are genuinely still open.]**
- **#315 ✓ CLOSED (corrected 2026-07-15)** Circumbinary/binary-star μ-gap sweep. This "OPEN" line was
  STALE — the exact same never-reconciled-duplicate pattern as `#500`. `#494`'s entry states
  "**Closes #315/#252/#255** positively" (4 binary-(k1,k2)-cycler μ-family reps + the first-ever
  Pluto-Charon (3,2) cycler). Caught while building the `## CURRENT STATE` dashboard at the top of
  this file.
- **#316** Cross-system cycler framework (Sun-Earth ↔ Earth-Moon manifolds) — OPEN, but possibly
  redundant with the #405/#411 SE<->EM heteroclinic-cycle work just above, which covers similar
  ground. No explicit link exists in the file either way — flagged for a human decision (merge/
  supersede or keep distinct), not resolved here.
- **#320** First quasi_cycler discovery sweep (blocked by #319) — **STALE, already resolved
  elsewhere.** #319 shipped (V1_qp/V2_qp/V3_qp) and #320's candidates were adjudicated
  2026-06-30 (net V0-known/not-novel) — see the #320 entry earlier in this file. This duplicate
  line predates that resolution and was never removed.

**Infrastructure + polish:**
- **#307** #289 Phase 5 (DSM + multi-rev + eccentric Tisserand) — gap from #300/#302 — **STALE,
  already resolved.** "#289 CLOSED 2026-06-25: ... Phase 5 (#307) closed" (see #289's entry
  elsewhere in this file).
- **#310** Single-orbit prioritizer adapter (#284 architectural gap) — **STALE, already
  resolved.** Landed elsewhere in this file ("closing the #284 architectural gap").
- **#317** PINN-based pre-filter for sweep-impossible regions — genuinely OPEN, no resolution
  found anywhere in the file.
- **#321** Multi-threaded inner-loop compute (joblib wrappers — 4-8× sweeps on multi-core) —
  **STALE, already resolved with a BETTER number.** Landed elsewhere in this file at a proven
  5.06× speedup (this line's "4-8×" was only an estimate), already reused by #343's 12.5×
  combined speedup.
- **#322** Tulip petal_count z0-collapse bug fix (in flight) — surfaced by #313 — **STALE,
  already resolved and now in active use.** Landed elsewhere in this file (folded into the
  "background-stable" infrastructure bundle) and is now itself used as a live gate
  ("the #322 z0-amplitude gate") by later tasks (e.g. #534's Phase 2 exploration).

## SESSION YIELD (HONEST RECKONING — POST #322 + #327 + #330 + #331)

- **0 admitted novel cyclers; 1 lit-fresh + physically-valid SILVER candidate
  that FAILS V2-moontour (drift + closure both above strict floor) but with
  bounded oscillation + monotonic-but-sub-0.5-km/s closure over 10 cycles.**
  Candidate ID `repeated-moon-uranus-00000041` at sequence Umbriel-Oberon-Umbriel
  (1,1). Closure 0.0252 km/s (2-moon convention, 24×24 basin floor) at V_inf =
  (0.92, 0.96, 0.89) km/s — ALL passing the #324 physical-sanity gate.
  Independent DOP853 cross-check residual 2.7e-11 nondim (5 orders of magnitude
  below the discipline gate). Offline literature_check against 35 KNOWN_CORPUS
  anchors returns `not-found` at confidence 0.40 (necessary-not-sufficient per
  `feedback_literature_novelty_check_baseline`; wider literature pass required).
  ML flagger p_fp 0.591 (below 0.75 SILVER threshold).
  **#330 (V2-moontour, 2026-06-16) verdict: FAIL_QUASI_BOUNDED.** At
  `n_cycles=3` strict gates: drift 5.2e+05 km (gate 5e+04 km), closure
  0.123 km/s (gate 0.05 km/s). Across 10 cycles drift OSCILLATES (cycle 5
  returns to 8.6e+04 km — near-resonant at 4.991× the Umbriel-Oberon
  synodic period); closure grows monotonically 0.025 → 0.349 km/s. Every
  Lambert leg converges in every cycle; the SILVER is geometrically valid
  but at this phasing it is a *near-resonant tour*, not a true cycler.
  **#331 (V3 6D n-body, 2026-06-16) verdict: PASS.** REBOUND IAS15
  agrees with the V2 (DOP853+Lambert) per-cycle drift series at every
  n_cycles in {3, 5, 10} to nanometer/micrometer precision (max |V3-V2|
  = 1.8e-6 km at n=10, vs 100 km agreement floor). The bounded-drift
  signature is a REAL property of the circular-coplanar Keplerian model
  — NOT an artifact of the V2 driver's DOP853+Lambert internals. The
  v4.7 quasi_cycler reading stands.
  Catalogue admission still blocked by (a) V4 HFEM Uranus (#332 —
  real-ephemeris with SPICE Uranian kernels), (b) wider lit pass (#329).
  See `docs/notes/2026-06-16-330-moontour-v2-phase2.md` +
  `docs/notes/2026-06-16-331-v3-nbody-phase3.md`.
  **The candidate is geometrically real and integrator-independent but
  not yet novel-claimable.**
- **Discipline held throughout** — no admission, no novelty claim, both gated
  exactly as the rule requires. Every well-scoped search hits published material
  EXCEPT this one row, which clears all offline guards and a fresh post-hoc
  physical-sanity gate (#324) but still has the discipline ladder to climb.
- **2 admitted mga_tour rows** (Tito 2018 + Heaton-Longuski 2003 Uranian Tour
  U00-01 under v4.7 scope expansion). #336 (2026-06-16) admitted the second:
  catalogue ID `heaton-longuski-2003-uranian-tour-u00-01`, V0 evidence is
  Tables 3 + 5 verbatim from the JSR paper (DOI 10.2514/2.3981). This also
  RESOLVES the #329 paywall gate (PDF acquired by the user 2026-06-16) and
  confirms Heaton-Longuski 2003 does NOT prior-publish the #327 SILVER —
  the 0.92 km/s value in Table 5 is the terminal Ariel rendezvous, not an
  Umbriel-Oberon-Umbriel cycler at (0.92, 0.96, 0.89). Frozen-census ratchet
  bumped NOT_TWO_BODY 0 -> 1 (bodies=[E, J, U, Titania, Oberon, Ariel, Umbriel],
  no period block — a one-shot tour). Catalogue 281 -> 282.
- **Genuine Koblick tulip-orbit characterizations: 1 system, NOT 5** — Earth-Moon
  Np=2..6 confirmed real 3D tulips under the #322 z0-amplitude gate. The
  #281/#283 claims of Jupiter-Ganymede / Saturn-Titan / Neptune-Triton / Pluto-Charon
  Np=2-3 "tulip" matches were **planar Np-petal collapses misidentified as 3D
  tulips by the pre-#322 broken petal_count gate**. Real boundary: the Koblick
  family does NOT μ-scale cleanly to small-μ moon systems (μ ≲ 1e-4); it's
  structurally Earth-Moon-specific.
  **Catalogue impact: zero** — none were ever written to `catalogue.yaml`; the
  JSONL-as-staging-not-truth discipline + hold-writeback rule protected the
  catalogue from the bug. Per `feedback_bugfix_invalidates_past_searches`, the
  systematic re-verification in #322 (data/tulip_topology_reverify_322.jsonl)
  caught all 7 false positives.
- **2 quantified near-misses** (UNAFFECTED by #322 — these are repeated-moon,
  not tulip) — Saturn Rhea-Dione 0.107 km/s (confirmed binding genome ceiling);
  Uranus Oberon-Titania 0.062 km/s (first quantification at this system;
  Uranus has NO published existence prior — highest novelty leverage)
- **1 publication-equivalent structural finding** (#313): Sun-perturbation
  effects do not transfer from Sun-Earth-Moon to Sun-Jupiter-moon for L1
  substitute families (Δx0 < 1e-9 across full mu_sun continuation vs Δx0 ~ 1e-4
  for the Sun-Earth-Moon case; geometric reason a_sun_LU ≈ 1160-1845 at Jupiter
  vs 388.8 at EM makes per-particle Sun acceleration ~10× weaker)
- **4 multi-week Track-A capabilities operational** (3D / BCR4BP / QP-tori /
  epoch-aware MGA — all 4 of #286's 5 axes that had concrete entry points)
- **Catalogue scope expanded** from 1 class to 4 classes (cycler / quasi_cycler
  / precursor_mga / mga_tour) with full toolchain (schema + lit corpus +
  website filters)
- **KNOWN_CORPUS** gained 12 MGA/tour anchors + Antoniadou-Voyatzis 2018
  spatial CR3BP anchor
- **3 standing memory rules added**: scope expansion taxonomy, Andreu
  canonical-momentum gotcha, speculative-high-effort-required

See `docs/notes/2026-06-16-frontier-scoping-er3bp-bcr4bp-3d-qp-epoch.md`
for the cost / existence-prior / first-IC scoping that drove this queue.

---

# Project state at a glance (updated 2026-06-14)

This top section is the orientation map for a contributor returning to
the project: what's done, what's in progress (with plan-file pointers),
what's blocked and why, and the prioritised human-actionable items. The
lettered Q&A log (A–H) below it is the per-entry catalogue-sourcing
audit trail and is unchanged in spirit.

## VALIDATION CEILING REACHED (2026-06-09) — read this first

The catalogue (280 rows) is at its **data-limited validation ceiling** via
current methods + held papers: **V3:2** (S1L1 `russell-ch4-4.991gG2` ballistic +
`russell-ch4-8.049gGf2` #188 powered), **V2:6** (Aldrin outbound + the 5 Ross
stable EM cyclers, #229), **V1:21** (the 2026-06-10 #181 writeback + the 3 Liang
CGE members, #222, + the 3 Braik-Ross 2026 common-energy CR3BP cycler
reproductions #249), rest V0. This is the honest boundary, established by exhaustive triage (#170
App-C batch, #177 self-seeding over all 212 unsourced rows, #172 Phase 6 novelty
sweep) — NOT a backlog to grind.

- **The ~200 V0 ocampo rows are V0 by a PUBLICATION gap, not laziness.** They use
  Russell's `n.m.k±j` SUMMARY-table format (descriptor + summary V∞ + period) —
  enough to CATALOGUE (V0=sourced), but the per-arc / per-leg reproducible state
  needed for V1+ was published only for ~9 PARENTS (Appendix C), never
  per-member. "Insufficient data but backed by full papers" is exactly true and
  not a contradiction. **V0 IS the deliverable for the bulk** — a faithful
  sourced census of the literature; V1+ is the smaller set WE independently
  reproduced.
- **S1L1 is SINGULAR** — the one row where coplanar descriptor + real-eph family
  + sourced anchors align. #177 confirmed the others close geometrically but miss
  their v∞ anchors (off-family) or lack the descriptor (sub-Mars / summary-only).
- **Past the ceiling is NEW-INPUT-gated** (not more iteration): (1) acquisitions
  (#116 — Russell&Ocampo 2006 / McConaghy 2006 full text; human-access,
  uncertain payoff), or (2) a more-capable genome (multi-arc / low-thrust /
  broken-plane — would reopen off-family rows AND the Jovian empty region; the
  Phase 6 §6b capability-subsumption gate auto-re-sweeps when one ships).
- **Machinery is complete**: corrector + MBH + continuation + free-return-chain +
  self-seeding + n-body harness + GMAT V4 lane (live, headless) + Phase 6
  discovery pipeline with the method-versioned empty-region registry.
  (2026-06-10: + the CR3BP Tier-2 propagator/corrector — `core/cr3bp.py` /
  `search/cr3bp_periodic.py` — and the VILM-leveraging endgame solver.)

**CORRECTION (2026-06-10) — the ceiling was WRONG for the descriptor-bearing
rows.** A Stage-B closer bug (the real-eph Lambert used the coplanar G-arc
*branch* ToF instead of the row's tabulated *signature transit*) inflated the
emerged Mars V∞ ~1.6–2.1×, producing a false "triple-confirmed off-family"
consensus — self-seeding (#177), the DSM corrector, and MBH all inherited the
SAME upstream bug, so their agreement was never independent evidence. With the
fix (signature-transit ToF + a joint (epoch, ToF) closer), ALL 6
descriptor-bearing `russell-ch4` rows close on real DE440 to both sourced
anchors (≤0.08 km/s; 6.44Gg3 to 0.00 km/s) and pass §14 V1 mechanics; the
proposed promotions (4 V0→V1 + 2 real-eph re-confirmations of existing circular
V1 rows; all 6 V3-CANDIDATES) were APPROVED + APPLIED 2026-06-10 (commit
`cec9b90`, census V1 11→15, pushed) and re-adjudicated **6/6 UNCHANGED** after
the #195 closer fixes (2026-06-11, `58674f6`/`d7f0c87`) — the writeback STANDS. See
`docs/notes/2026-06-10-dsm-tof-artifact-correction.md` and
`docs/notes/2026-06-10-tof-fix-closure-results.md`. SEPARATELY, the McConaghy
2004 dissertation mining CONFIRMED the ceiling for everything else: per-member
reproducible data (date + V∞ + closest-approach + leg ToF, DE405) exists ONLY
for S1L1 (see the 2026-06-10 section below); the ~200 ocampo rows stay V0 —
no per-member data was ever published.

## 2026-06-14 — discovery-program pivot + multi-arc convergence SOLVED

**STRATEGIC PIVOT: novel cyclers are the deliverable; validation is the means.** The
catalogue is at its data-limited validation ceiling and the literature is mostly
exhausted (paper sweep #250 found no new genome). The toolkit is mature, so the
program shifts to DISCOVERY. Spec: `docs/notes/2026-06-13-discovery-program-spec.md`.
Three tracks: **(A) richer genome** — #254 repeated-moon multi-rev (Liang CGE
replication, Jupiter/Saturn/Uranus swept empty) + **#266 TULIP-ORBIT GENOME
OPERATIONAL** (Phase 1-3 complete, commits `12fd15c`/`2eb1a56`/`d8b8210`: Sundman
regularization → NRHO continuation → period-multiplying bifurcation detector →
family-switching corrector; end-to-end `find_tulip_via_continuation(np_target=2)`
lands a Np=2 tulip at T=2.746 TU, J=3.058, petal_count=2 — within 0.38% of
Koblick 2023 AMOSTECH Table 4 source. Phase 4 follow-on (multi-shooting for k>=3)
tracked separately); **(B) prioritizer** — #249 RESOLVED 4/4 (all Braik-Ross
cycler members C11a/C11b/C21/C32 reproduced; scorer ungate-ready) + #239
impulsive merge done in #263 two-tier prioritizer; **(C) discovery-campaign
daemon** — #253. Regime arc: left input-bound → briefly build-bound (genome +
prioritizer + harness) → then CPU-bound (the wanted regime: discoveries scale with
cores×time). The agent-hang failures were long solver compute outgrowing one-shot
agents → that compute belongs in the #253 daemon.

**Multi-arc convergence SOLVED (#248).** The hard E-E-M-M closure blocker (neither lane
converged, #244) is cracked: the clean `search/multiarc_closure.py` harness (canonical
metric + epoch-safe eval + discrete resonant-return seeds + multi-start) converges
`mcconaghy-2006-em-k2` at **0.0987 km/s, reproducibly** on the FBS-analytic lane —
where Lambert (0.1044) misses. First direct evidence FBS reaches a basin Lambert
can't → **unblocks #245**. Marginal (0.0013 under gate), mcconaghy-only so far;
follow-through = russell rows + parity re-run (now #245). The scratch-era
0.163/2.06/58 numbers were buggy/metric-confused — superseded. Note:
`docs/notes/2026-06-13-multiarc-seed-basin-fix.md`.

**Binary-star μ-continuation — clean NEGATIVE (#252).** Continuing the EM (1,1)/(3,1)
families in μ does NOT reach the journal's binary-star cyclers: the cycler branches go
linearly unstable before μ=0.5/0.3, and the stable orbit found at μ=0.5 was a
one-primary librational orbit (the topology gate caught the false "it closed and it's
stable"). The paper's binary-star cyclers exist but are NOT the analytic continuation
of the EM member — they need a DIRECT fixed-μ search. Note:
`docs/notes/2026-06-14-binary-star-mu-continuation-discovery.md`.

**Infra:** WebFetch broadly enabled (blanket allow) so agents can fetch arXiv/PDFs;
standing rule — every fetched PDF is filed + committed to the private paper corpus
(4 filed this session). Roberts-Tsoukkas journal (#251) = no numeric tables
(figures only) → no sourced row, only a qualitative prior. Long agents now commit
incrementally (quota walls + polling-loop hangs lose un-committed work otherwise).

## 2026-06-13 — Ellison FBS lane, errata system, JPL oracle, writebacks, reachable-set scorer, surrogate corpus

**Catalogue writebacks applied since 2026-06-11** (census now V3:2 / V2:6 / V1:18,
rest V0; 277 rows):
- **Ross 5 stable EM cyclers → V2-ballistic** (#229): long-span discriminating
  propagation evidence; held writeback applied.
- **Liang A/B/C CGE triple-cyclers → V1** (#222): idealized same-model CGE
  reproduction. Member D (#223) re-propagated n-body with the JUP365 kernel.

**Tooling / verification landed:**
- **Assumed-errata system** (#228): `data/errata.yaml` (14 entries) + a respectful,
  evidence-first `/errata` page on cyclers.space. STANDING RULE: public defect
  claims are benefit-of-the-doubt, falsifiable, typesetting-slips-happen framing.
- **JPL Three-Body Periodic-Orbit oracle** (#116 capability): live API consumable
  via `verify/jpl_periodic_orbits.py`; our CR3BP propagator reproduces JPL ICs to
  ~2.8e-10 nd at matched μ. Resolves "JPL 3BP catalog" as a usable cross-check.
- **Solver-outcome logger** (#210): `search/outcome_log.py` — opt-in `(genome →
  outcome)` JSONL capture for a future ANN surrogate (Ozaki 2022); NO-OP unless
  `CYCLERFINDER_OUTCOME_LOG` set; never read back into validation (hard boundary).

**Ellison FBS Path-B lane (#226) — DONE, parity-proven.** A Lambert-free
match-point leg corrector with analytic Jacobian (Shepperd STM): `dsm_leg_correct_fbs`
+ `core/fbs_match_point.py` (single-leg defect/Jacobian, boundary v∞/epoch columns,
multi-arc `chain_defect_jacobian`, `flyby_coupling_block`). All gradients FD-validated
to ~1e-9; Phase-6 closure parity vs Lambert to **6.1e-13 km/s**. Full suite
1655 passed. **Opt-in / additive** — the Lambert lane is untouched.

**FBS evaluation — what it is and isn't (#242, clean negative).** Tested whether
FBS closes the historically un-closeable multi-arc rows (S1L1, 6.44Gg3): it does
NOT close anything Lambert can't. Both lanes solve the identical per-leg BVP and
agree to ~1e-12; on long multi-rev arcs FBS shooting is *more* seed-sensitive
(needs the Lambert basin as a seed). KEY REFRAME: Ellison FBS's value is **analytic
gradients for ΔV OPTIMIZATION**, not feasibility-finding — #242 mis-scoped it as a
Lambert replacement. Conclusion: **keep the opt-in FBS corrector as an independent
cross-check; do NOT default it.** See `docs/notes/2026-06-13-fbs-hard-rows-closure-attempt.md`.

**FBS optimizer fair trial (#243) — clean POSITIVE → ADOPT.** Tested in its actual
role (analytic gradients as a ΔV-minimizing multi-leg optimizer backbone vs the
same NLP solved with finite-difference gradients), FBS wins decisively across
Aldrin / Russell / 6.44Gg3-class problems (40 cold seeds each, same-model optima,
Jacobian cross-checked to ≤2.2e-7): **robustness 3–16× better** (Aldrin 100% vs 28%
feasible; opt% 100 vs 2), **wall-clock up to 4.8× faster** (gap widens with NLP size
— FD pays O(n_vars) evals/gradient), **optimum equal-or-better** (FD found a
strictly worse optimum on Russell). Commit `51bb455`; `tests/search/test_fbs_optimize.py`
4 passed. Verdict: **adopt FBS as the optimization engine → proceed to #244.**
Caveats for #244: confirm on the real `dsm_chain_correct` lane (trial used scipy
SLSQP; EMTG pairs FBS with SNOPT), and wire patched-conic flyby-continuity
constraints into the NLP (where the analytic advantage should compound). The
default-promotion ladder #243→#246 is re-pointed at the optimizer role. See
`docs/notes/2026-06-13-fbs-optimizer-fair-trial.md`.

**FBS optimizer adoption + parity sweep (#244) — DONE; verdict HOLD #245.** The
opt-in FBS-analytic-gradient backbone is wired into the real `dsm_chain_correct` /
`close_row_dsm` lane (commits `4abb1b8`/`b25f2d1`/`c2fb57a`, default `gradient=
"lambert"` byte-unchanged) WITH patched-conic flyby-continuity constraints. The
catalogue-wide parity sweep on the REAL corrector found: FBS is closer-to-feasible
on every row (optimum-quality half holds directionally) but **slower, and NEITHER
lane converges** on the multi-arc E-E-M-M / E-E-E-M-M rows. So there is no
convergence parity to flip toward → **do NOT flip the default (#245)**; FBS stays
the opt-in it already is. The real, decision-changing blocker is UPSTREAM of the
gradient: **single-charged-seed basin selection** for multi-arc rows (now tracked as
**#248**, which gates #245). Re-run `scripts/fbs_optimizer_adoption_parity.py` once a
row converges on either lane. See `docs/notes/2026-06-13-fbs-optimizer-adoption-parity.md`.

**Reachable-set accessibility scorer (#236, from Braik-Ross #230) — built, GATED.**
A reduced (x,y,θ) heading-fan forward/backward reachable-set overlap scorer +
proxy-ΔV graph centralities, as a continuation/family-selection prioritizer. The
reproduce-before-trust validation gate at C_J=3.1294 was a PARTIAL pass:
**R21-S persistent-hard-access reproduced** and Table-2 periods confirmed exactly
(LL1/LL2/DPO/R21-S/C11b to the digit, C32 to ~1.1%), but **C32-dominant was a
faithful negative** — recorded as `xfail`, NOT forced. Root cause is a
member-RECOVERY gap, not a method bug: only **6 of 13** representatives recover via
the available 1-DOF perpendicular-x-crossing symmetric corrector at the off-stable
common energy (C11a/C21/R21-U/R31 don't; JPL doesn't expose 5:2 R52), and C32's
dominance is a full-13-node-chaotic-sea property that doesn't survive the
truncation. Scorer is therefore **GATED** (not applied to our families).

**Ungate attempt (#247) — 9/13 recovered.** A network-independent free-(x0,t_half)
perpendicular-crossing corrector recovered 9 of 13 members (up from 6), each
confirmed against BOTH sourced Table-2 period AND Floquet σ. Unrecovered (at
that stage): C11a, C11b, C32, R52-U. Commits `5b48ecc`…`0ad93f7`. See
`docs/notes/2026-06-13-reachable-scorer-ungate.md`.

**Final ungate (#249, 2026-06-14/15) — 4/4 cycler members reproduced; scorer
ungate-ready.** All four Braik-Ross cycler members rigorously reproduced at
C_J=3.1294 (period + (k1,k2) winding topology + prograde + Radau cross-check):
C11a 42.1405d (1,1), C11b 55.9590d (1,1), C32 78.6126d (3,2), C21 84.5331d (2,1)
— see `docs/notes/2026-06-14-249-unstable-member-recovery-plan.md`. Method:
all-roots enumeration + winding classifier (rejects period-impostors). C21
required Ross-RT 2025 AAS 25-621 Table 4's unrounded Jacobi 3.129389531088256
(B-R's "C_J=3.1294" is 5-sig-fig DISPLAY rounding; the (2,1) family has ΔC≈4e-12
so literal-3.1294 sits outside it). Standing lesson: published rounded invariants
are display, not literal — see `feedback_published_rounded_values_are_display`.
New infrastructure (all tested, reusable): fold-turning pseudo-arclength
continuation in Jacobi (`e34edb2`), winding-(k1,k2) classifier (`93fb330`),
asymmetric (general) periodic-orbit corrector (`23b980e`). Final 4/4 disposition
`a19eb24`.

**C32-dominance gate re-run (#262, 2026-06-15) — FAITHFUL NEGATIVE on the 12-node
source-confirmable set; scorer STAYS GATED.** With all four cyclers (C11a, C11b,
C21, C32) wired into the scorer via `recover_all_cyclers_braik_ross` (each at its
correct per-family Jacobi — literal 3.1294 for C11a/C11b/C32, the unrounded
3.129389531088256 for C21) plus the eight offline-confirmable nodes (LL1, LL2,
DPO, R21-S, R21-U, R31-S, R31-U, R52-S), C32 does NOT emerge as the dominant
family node under our scorer. Observed centrality ranking:
  - strength argmax: **C11a** (0.00741); C32 = 0.00666 (rank 6 of 12)
  - harmonic argmax: **C21** (0.00743); C32 = 0.00676 (rank 5 of 12)
  - betweenness argmax: **R21-U** (0.0606); **C32 has ZERO betweenness** (no
    relay role at all)
The Braik-Ross Table-4 prediction (C32 rank 1 across all three metrics at
0.2850/0.2891/0.5000) is not reproduced. Whether this reflects a real limitation
of our heading-fan reachable-set proxy on this model (voxel-grid resolution,
T_a horizon, n_seeds/n_fan, or the missing R52-U node) or a deeper method
mismatch is not yet diagnosed — parameters were NOT tuned. R21-S hard-access
half of the gate still REPRODUCES (R21-S in the bottom-3 of strength and
harmonic closeness alongside the other stable resonants). The C32-dominance
gate-test stays `xfail`; a companion `test_validation_gate_c32_undominant_faithful_negative`
records the observed ranking as a passing test so the negative is captured (not
just an xfail). The scorer is **NOT applied to rank our families**. **Next:**
diagnose the rank disagreement (parameter sensitivity / R52-U recovery / horizon
choice) before any unblocking; #239 (Zhou-Armellin reachable-set) still pending.

**Papers digested (#230–235):** Braik & Ross 2026 Orbital Networks (reachable-set
family-accessibility — method-only, no new tuples), arXiv:2509.12671 HOTM fixed
points (#237 negative-registry entry), Şaloğlu-Taheri JAS VoR, Shepperd 1985
(#225 goldens), Fan 2025 (triage), Tito-MacCallum 2013 free-return reproduced as a
cross-check (#238, withheld per cyclers-only scope). All background papers READ,
not title-dismissed.

**Open after this session:**
- **#248** (NEW) — multi-arc seed/basin fix so the E-E-M-M rows converge on either
  optimizer lane; GATES the FBS default flip (#245). The decision-changing finding
  from #244.
  **✓ RESOLVED — cross-reference added 2026-07-15, full-file audit:** see "Multi-arc convergence
  SOLVED (#248)" elsewhere in this file (a later-dated note appearing earlier in the file per its
  DELTA-newest-first convention) — the hard E-E-M-M closure blocker was fixed. This "(NEW)" header
  predates that fix and was never updated; #245/#246 below are correspondingly unblocked, not
  genuinely still HELD.
- **#245/#246** — FBS optimizer default flip + documentation; HELD behind #248
  (no convergence parity to flip toward yet). **Per #248's resolution above, this gate should be
  re-checked — it may already be clear; verify before treating #245/#246 as blocked.**
- **#239** — Zhou-Armellin reachable-set spike; merges into the still-gated scorer.
- **#240** — KKT/surrogate amplifier; gated on the surrogate-target decision. Logger
  needs NO change — it already logs the DSM-chain lane via `optimize.py`; what's
  missing is a flyby-cycler DATA-GEN DRIVER, only if that surrogate is the target.
- **#115/#116** — blocked (no published V∞ / human-gated acquisitions).
- **Surrogate-corpus daemon** — still running (plain process, quota-proof); CR3BP
  `(genome → outcome)` tuples to gitignored `out/outcome_log/`, ~50/50 mix, nice +19.

**Quota note (2026-06-13):** several one-shot agents hit the session quota wall
mid-run; agents that committed incrementally kept their work, one that batched lost
everything — incremental commits are now the standing instruction for long agents.

## 2026-06-11 — review-and-harden + acquisitions wave 2

- **FULL-PROJECT REVIEW — verdict: "the equations are right."** 4 scoped
  reviewers (core math / search correctors / data+verify / performance) + a
  numerical math-verification agent (8 independent probe suites: FD-vs-STM,
  ∇C·f=0, Lambert-vs-independent-Kepler, flyby identities, Tisserand, frame
  conventions, rotating↔inertial round-trips). Every independently-derived
  formula agrees with the implementation at/near machine precision; ALL
  confirmed defects live in solver plumbing / conventions — every one in a
  stratum no published golden exercised. Consolidated record (findings table
  with dispositions, the doctrine, open items):
  `docs/notes/2026-06-11-project-review-results.md`.
- **Findings FIXED + adjudicated (commits):**
  - `joint_epoch_tof_close` epoch double-shift + vacuous lon diagnostics →
    `58674f6`; #195 adjudication `d7f0c87`: the #181 6-row closure re-run is
    **6/6 UNCHANGED to every recorded digit** — the V1 writeback STANDS.
  - Data-layer guards (#196) — the validation-level over-claim bypass closed
    at BOTH ends (nested `validation.level` validated + `apply_*` refuse
    unregistered promotions), registry-drift preflight, duplicate-id
    rejection, atomic catalogue writes → `091783a`; 448 data tests green.
  - R_x(−i) node-convention MIRROR in the inclined-circular backend + the
    ramped-elements `_tilt` (orbit normals sat 2×inc off DE440) →
    `278ff1a`/`1d6ad1b`. Blast radius adjudicated: the #120 3D-inclination
    negative STANDS (its decisive DE440 control was mirror-free), Tisserand
    is mirror-invariant, continuation results valid.
  - Lambert dT/dz spurious √C + log-compressed Illinois residual (rescues
    the dropped multi-rev high branches) → `f6a0460`; full suite green.
  - `crosscheck_leg` endpoint independence default-ON + poisoned-input fault
    test + the gate-classification convention (#197) → `ba55b2e`.
- **ZERO RETRACTIONS.** Every adjudication (joint-closer re-run, mirror
  blast radius, McConaghy tail) CONFIRMED the existing results. The
  McConaghy Table 7.1 rows-18–24 tail is a **SOURCE print defect** (Table
  7.1 prints Table 7.5's dates with an orphaned V∞/CA tail); our
  transcription is character-exact, and the DE440-emerged tail matches
  Table 7.5's printed values → reproduction upgraded PARTIAL-CONFIRMED →
  **CONFIRMED** (all 23 legs vs printed goldens; `53c0a92`/`a8c0928`).
  **#94 CLOSED.**
- **FALSE-CONSENSUS DOCTRINE (3rd incident: the 63 s shared epoch
  conversion).** After #180 (shared-ToF) and #197 (shared crosscheck
  endpoints), the review's probes confirmed a 63 s UTC/TDB conversion offset
  shared between the primary path and its "independent" cross-check.
  Doctrine now operational: consistency-vs-independence gate TIERS + a
  "shared with primary path:" declaration per cited gate (recorded in
  `src/cyclerfinder/data/validate.py`'s docstring), fault injection,
  positive controls, per-interface external anchors.
- **OPEN (in flight / queued):** #198 the 63 s UTC/TDB epoch offset
  (probe-confirmed; fix in flight); #212 Earth-Moon μ double-count in
  `cr3bp_system()` (−1.2%, found by the Ross mine §7; fix in flight; the
  2026-06-10 EM backfill needs a post-fix re-run; Saturnian/Jovian pairs
  expected unaffected at ≤~2e-4 — pending quantification); #206 Lambert
  blast-radius re-run; #201 perf batch; #202 fault-injection harness.
- **ACQUISITIONS WAVE 2** — forward-citation sweep of the load-bearing
  holdings (`d3dc4fd`, 13 HITs across 6 seeds) → user fetched 7 PDFs (filed
  in the private papers store; cite by publication only, never by
  path/repo) → 4 mining passes, notes all committed 2026-06-11:
  - **Ross & Roberts-Tsoukkas, AAS 25-621** (stable low-energy prograde
    Earth-Moon cyclers) — first publication of **5 STABLE fully-ballistic
    EM cycler families** (15–16-digit Jacobi C + period per stable
    representative, μ printed; no full IC printed, but every member is
    recoverable from (μ, C, T) via a 1-D symmetric-orbit solve) + 6
    critical-tangency manifold-tube goldens + the complete construction
    method. The mine FOUND the #212 μ bug. **5 proposed rows, review-gated.**
    Note: `2026-06-11-ross-roberts-tsoukkas-2025-mining.md`.
  - **Liang et al. 2024 JGCD** (Callisto-Ganymede-Europa triple cyclers) —
    **4 concrete ballistic CGE members** (3 idealized + 1 SPICE-ephemeris,
    per-flyby V∞/ToF/phases tabulated, ≥10 cycles). The Jovian empty-region
    bucket is **REOPENED via the §6b new-sourced-data arm, CONDITIONAL on a
    multi-rev + repeated-moon genome** (the members live in a topology class
    our swept genome does not cover, so the zero-rev EMPTY verdicts stay
    valid as conditioned). 3 paper errata recorded. **4 proposed rows,
    review-gated.** Note: `2026-06-11-liang-2024-cge-triple-cyclers-mining.md`.
  - **Wittal, Miaule & Asher IAC-22-C1.6.6** (full text) — **zero ICs of any
    kind**: the `wittal` NO_SOURCED_IC blocker is **SOURCE-PERMANENT** (a
    publication gap, not an access gap); the row's `model_assumption: cr3bp`
    is unsupported by the source (model never stated) → revisit #211.
    Catalogue enrichment only. Note: `2026-06-11-wittal-2022-iac-mining.md`.
  - **Cuevas del Valle et al. 2023** (*Aerospace*, Floquet station-keeping) —
    **first sourced CR3BP maintenance-cost anchor** (Table 1: 18.5–22.7 m/s
    per ~330 km insertion error, EM L2 halo, one period) + the
    Jacobi-error-state formulation confirming the #190 energy coordinate;
    cross-paper finding: the 2026 MPC paper's "L2 southern halo" IC is
    almost certainly an **L1** halo (relabeled). Note:
    `2026-06-11-cuevas-del-valle-2023-floquet-mining.md`.
  - **ML-surrogate trio (Ozaki 2022 / Leifsson 2022 / Wu 2024)** — Ozaki's
    DNN flyby-cycler blueprint is the direct template for a future
    combinatorial free-return-chain lane but **DEFERRED** (≈7e6-sample
    training floor, below breakeven at our scale; costless surrogate-prep
    items → #210); Leifsson/Wu background-only. Note:
    `2026-06-11-ml-surrogate-trio-triage.md`.
- **Catalogue state:** census **V1:15** live + pushed (the #181 writeback
  applied 2026-06-10, re-adjudicated UNCHANGED today); the 14 CR3BP SILVER
  Lyapunov members are all-gates-green (incl. the `f69d2b3` inertial
  REBOUND/IAS15 cross-check) and **await user disposition**; proposed-row
  queues review-gated: Ross 5 (stable EM), Liang 4 (Jovian CGE).

## 2026-06-10 — genome-capability day: VILM endgame, ToF-artifact correction, CR3BP Tier-2, two acquisitions mined

- **#179 — VILM-leveraging endgame solver BUILT + VALIDATED.**
  `src/cyclerfinder/search/leveraging_leg.py` (analytic near-root
  apse-quadratic leveraging legs + Γ-floor cross-check) and
  `src/cyclerfinder/search/endgame_graph.py` (multi-moon Dijkstra over the
  leveraging graph, soundness verified against brute force). The Saturnian
  leveraging novelty sweep came back **EMPTY** (best max-V∞ 12.775 km/s vs the
  6.0 km/s floor) — an honest method-versioned negative
  (`saturnian-titan-endgame-vilm-2026-06-10` in `data/empty_regions.jsonl`).
  The capability edge ("leveraging" ⊐ "single-arc") is registered for the
  Phase-6 §6b re-sweep gate.
- **#180→#181 — DSM multi-arc lane built; its negative RETRACTED same day;
  ToF fix closes all 6 descriptor rows.** The DSM multi-arc closure lane
  (`search/dsm_descriptor_seed.py`) was built and produced a "0 promotions,
  triple-confirmed off-family" negative — RETRACTED the same day as a
  shared-upstream-bug artifact (the ceiling CORRECTION above). The fix
  (signature-transit ToF + the joint (epoch, ToF) closer in
  `search/self_seeding.py`) closes all 6 descriptor-bearing rows V1-PASS on
  real DE440 (every emerged anchor ≤0.08 km/s of sourced; best is 6.44Gg3 —
  the row that triggered the investigation — at 0.00 km/s on both anchors).
  All 6 are **V3-CANDIDATES** (single-leg independent REBOUND/IAS15 confirm
  in-band; the multi-lap horizon-TCM is the named follow-up). **Writeback
  APPROVED + APPLIED 2026-06-10:** the 4 previously-untagged rows promoted
  V0→V1 (`3.78Gg3`, `6.44Gg3`, `3.64gGg3`, `5.30ggF3` — census V1: 11 → 15);
  the 2 already-V1 rows (`9.353Gg2`, `9.94Gg3`) carry the real-eph closure as
  ADDED EVIDENCE (no level change). V3 promotion still pending the multi-lap
  horizon-TCM. Notes: `2026-06-10-dsm-multiarc-closure-results.md`
  (carries the SUPERSEDED banner), `2026-06-10-dsm-tof-artifact-correction.md`,
  `2026-06-10-tof-fix-closure-results.md` (writeback record appended).
- **#182 — CR3BP Tier-2 SHIPPED** (8-task plan; Task 8 regression in flight at
  writing). `core/cr3bp.py` (rotating-frame EOM, Jacobi constant, 42-state STM
  propagation, DOP853 at 1e-12) + `search/cr3bp_periodic.py`
  (`correct_periodic` Newton single-shooting + `crosscheck_periodic`
  independent Radau re-propagation). Arenstorf sourced golden: the corrector
  converges to the published IC (closure 4.16e-9 in the gate; 7.7e-11 in the
  backfill run), Jacobi conserved ≤3.4e-12. The Moon was added to the
  `core/satellites.py` registry. **Earth-Moon backfill:** the Arenstorf row's
  CR3BP block is populated (proposed, review-gated); `genova-aldrin` + `wittal`
  came back NO_SOURCED_IC — and the Genova mining below makes the genova
  blocker PERMANENT for the CR3BP lane (not an acquisition gap). (2026-06-11:
  the Wittal full-text mining confirmed the `wittal` blocker is
  **SOURCE-PERMANENT** too — zero ICs of any kind in the paper, and the row's
  `model_assumption: cr3bp` is unsupported by the source; see the 2026-06-11
  section.) **Saturnian
  midsize-moon discovery:** the first run produced 11 contaminated
  "convergences" (8 were libration-point equilibria / a period-collapse — a
  fixed point trivially "closes" for any period); deleted and regenerated
  under a v2 degeneracy gate (equilibrium rejection / period floor / dedup /
  independent Radau crosscheck) → **14 genuine tiny-amplitude (km-scale)
  L1/L2 Lyapunov members** across Mimas/Enceladus/Tethys, SILVER
  (`data/cr3bp_silver.jsonl`, method `cr3bp-lyapunov-corrector-v2`), honestly
  flagged pipeline-validation-not-novel (Lyapunov-centre families). Amplitude
  growth needs half-period symmetry / multiple shooting — recorded as a
  method-capability gap. Notes: `2026-06-10-cr3bp-backfill-results.md`,
  `2026-06-10-cr3bp-moontour-results.md`.
- **Acquisitions MINED** (held reference papers, cited by publication only —
  the offline reference store is private and is never linked from this repo):
  - **McConaghy 2004 Purdue PhD** (*Design and Optimization of Interplanetary
    Spacecraft Trajectories*, UMI/ProQuest 3166673) — Ch. 7 Tables 7.1–7.5
    give four full per-encounter DE405 S1L1 itineraries (calendar date, V∞,
    closest-approach distance, leg ToF; 24–27 encounters over 33 yr) = exactly
    the per-member state the V0 rows lack, **for S1L1 ONLY**. Table 6.2 pins
    S1L1 = `2g(2.8277, 657.97°, U)…g(1.4508, 522.29°, L)` (multi-arc, two
    generic E-E legs). Table 5.8 gives U0L1(2.7540) ballistic (ΔV=0.00) with
    its V∞ pairs (Earth 11.3/11.3, Mars 14.0/5.4 km/s); Table 5.4 gives
    2L3 = Byrnes Case 1 (V∞E 5.65 / V∞M 3.05, simple model). Every other
    family is summary-exemplar-only → **ceiling CONFIRMED except S1L1**.
    **#94 UNBLOCKED** on the source axis — and **CLOSED (2026-06-11)**:
    Table 7.1 reproduced per-leg on DE440 (`53c0a92`); the rows-18–24 tail
    adjudicated a SOURCE print defect (Table 7.1 prints Table 7.5's dates
    with an orphaned V∞/CA tail; our transcription character-exact);
    reproduction CONFIRMED across all 23 legs (`a8c0928`). Notes:
    `docs/notes/2026-06-10-mcconaghy-2004-dissertation-mining.md`,
    `docs/notes/2026-06-10-mcconaghy-table71-reproduction.md`.
  - **Genova & Aldrin 2015** (NTRS 20150018049) — the 3-petal EM cycler does
    NOT exist in the pure CR3BP (it requires solar gravity, i.e. a
    bicircular-or-better model; stated explicitly on p.3); no Jacobi constant,
    rotating-frame IC, or nondimensional period is published → the row's
    CR3BP-backfill blocker is **PERMANENT-for-CR3BP**, not an acquisition gap.
    Sourced values extracted and held for review (553-d apsidal super-period,
    ~24° inclination, 20–62 m/s/cycle station-keeping, 39 m/s/month abstract
    average); flags a `model_assumption: cr3bp` → bicircular reclassification
    question. Note: `docs/notes/2026-06-10-genova-aldrin-2015-mining.md`.
  - **Also newly held:** Vallado 1991 USAFA-TR-91-6 (sourced algorithm test
    cases); Ellison et al. 2018 JGCD (analytic MGALT/MGAnDSM gradients);
    Cuevas del Valle 2026 EuroGNC (CR3BP rendezvous); Iorfida 2016; Shakouri
    2019. Vallado *Fundamentals of Astrodynamics* 4th ed full text is
    consultable online (archive.org djvu text — OCR caveat: equations/tables
    unreliable, independent check required before any golden use). Szebehely
    *Theory of Orbits* held, mining pending (#185 — CR3BP periodic-orbit
    goldens).

## Validation campaign — 2026-06-08 (V0-lift, three-track)

Live validation census after the writeback: **V2: 1, V1: 11, V0: 256** (268 total).

- **+6 V1 (closer sweep, #(3)):** ran the #137 circular-coplanar free-return
  closer over the newly-ingested rows; `russell-ocampo-3.1.1+2 / -3.1.3+0 /
  -4.1.1-4 / -4.1.2-2 / -4.1.4-1 / -4.6.3+0` CLOSE-AND-MATCH + §14 V1 mechanics →
  promoted V0→V1 (commit `12a0e9e`; gate `tests/search/test_closer_sweep_v1.py`;
  runlog `data/runs/closer-sweep-2026-06-08.jsonl`). 2 CLOSE-but-multi-arc + 8
  NO-CLOSE correctly refused. See `docs/notes/2026-06-08-closer-sweep-v1-candidates.md`.
- **Continuation V1→V3 batch (#(1)): 0 clean V3.** Built
  `search/continuation_batch.py` + per-row V3 gate. The best candidate
  (`russell-ch4-5.75ggF3`) passes the 3.0 km/s plausibility bar but FAILS spec
  §14's actual V3 budget (`horizon_tcm_mps` = 120 m/s; its 3-lap TCM = 1210 m/s)
  and the 4–5 lap horizon. Held — no V3 written. Aldrin remains the only family
  ballistic across the full 20–30 yr horizon. See
  `docs/notes/2026-06-08-continuation-batch-results.md`.
- **Appendix C (#(2)) confirmed the only path for the 15 Rall rows** — not
  closer-reachable, not continuation-seedable (no per-arc aphelion to seed). Held
  pending the transcription slog.
- **Finding:** the V3 frontier is gated by *physics* (cycler geometry breaks the
  ΔV budget over decades), not tooling. Cheap V0→V1 wins are now spent; the
  remaining frontier is the multi-arc genome (#94 / ~10 Russell multi-arc rows)
  and the human-gated acquisitions (#116).

### Multi-arc frontier — primitive-mismatch finding (2026-06-08, #162→#163)

- **#162 multi-arc basin experiment (AMBIGUOUS, floor-drop):** the dsm_leg
  residual had no bend-feasibility term; adding one (default-off
  `charge_flyby_continuity`) dropped the 6.44Gg3 floor to ~9 km/s (lowest to
  date: 9.40→29.9→26.9→~9) but proved bend was never the wall — freeing the
  intermediate departure-V∞ direction makes every flyby trivially bend-feasible
  (flyby_dv→0); the irreducible term is the **ΔV_DSM budget** (~7.4–9.7 km/s
  across a 9-pt seed sweep, 0/61 MBH hops). See
  `docs/notes/2026-06-08-multiarc-basin-selection-results.md`.
- **The reframe (#163):** dsm_leg is the WRONG primitive — its ΔV_DSM-min
  objective does not pin V∞ to the anchor, so it floors off-anchor regardless of
  freedom. Russell's construction is the **generic-return arc** (the g/G
  descriptor) = #137's `free_return` primitive, which pins V∞ to the anchor and
  closed 8 symmetric rows + the 6 promoted today. 6.44Gg3 / S1L1 are the TWO-arc
  cases a single ellipse can't represent. We were re-deriving a *published*
  result (Russell's 0.509 km/s cycler) with the wrong building block. #163 chains
  TWO `free_return` arcs (anchor-respecting residual, NOT ΔV-min); a
  circular-coplanar close would also hand the continuation driver (#158) the seed
  S1L1/#94 has lacked → the path to V3. RUNNING; verdict →
  `docs/notes/2026-06-08-free-return-chain-results.md`.
- **#163 LANDED — reverses the empty-set lean.** The two-arc chain REACHED the
  V∞ anchors where every prior approach floored off-anchor: 6.44Gg3 emerged
  6.48/3.70 + 6.36/3.82 vs sourced 6.44/3.74 (vinf_res **0.081 km/s**),
  bend-feasible. The obstruction MOVED to a quantified **descriptor-ToF gap**
  (V∞-pinned ellipse period × integer revs leaves the ToF in a gap; 6.44Gg3
  0.50 yr). **S1L1 `russell-ch4-4.991gG2` is CLOSE-LEANING** (own anchors
  4.99/5.10): V∞ <0.1, G-arc ToF near-exact (2.810 vs 2.8096), g-arc 0.14 yr
  off — the strongest circular-coplanar continuation seed yet. Frontier is NOT
  empty-set; it is V∞-reachable / ToF-quantised. **Next:** continue the S1L1
  two-arc seed to the real ephemeris (real eccentricity breaks the circular
  integer-rev quantization → may close the ToF gap → #94, seedable for the first
  time). Commits `f87a6be`/`234caf1`.
- **⚠️ #165 RETRACTS the #164 "closes" claim — S1L1 is two-body-closeable,
  n-body-UNCONFIRMED; #94 REMAINS OPEN.** The independent REBOUND/IAS15 cross-check
  over REAL DE440 planets got **DRIFT**: the #164 geometry crosses Mars's radius
  at the right time/speed but ~110° from real Mars → **2.6 AU miss** (band ≈
  0.0116 AU), V∞ 38.8 not 5.10. Cause: #164 closed V∞+ToF+radius in a model with
  real Mars *eccentricity* but idealised *longitude*; longitude rendezvous against
  true DE440 Mars was never enforced. A best-phase epoch scan still misses by
  0.234 AU → the single-ellipse-per-arc geometry + real-Mars rendezvous look
  jointly unsatisfiable in this construction. The continuation closed the
  ToF-quantization axis (real progress); the **phasing/longitude axis is the
  remaining wall**. No V3. Gate `tests/nbody/test_s1l1_nbody_closure.py` pins the
  DRIFT. Commits `3a371be`/`18961b5`; note
  `docs/notes/2026-06-08-s1l1-nbody-crossval-results.md`. [The #164 entry below is
  SUPERSEDED by this retraction.]
- **#166 SOURCE DIG — DRIFT cause found, exact fix seedable, NO acquisition.**
  (1) Topology was wrong: in `g(1.4612)+G(2.8096)` only the UPPERCASE G is the
  Mars-transit leg (Russell §4.8); lowercase g is a pure Earth-to-Earth free
  return that never reaches Mars (DE440-verified, closest 1.05 AU). #164/#165
  modeled both arcs as Mars-crossing → the 2.6 AU drift. Correct sequence:
  E → g(E-E) → E → G(E-M-E, longitude rendezvous) → E. (2) Russell App-C gives
  the exact real-eph state; it reconstructs on DE440 to **1.7 km** (Mars flyby
  2027-06-13, v∞ 5.248 = published, Mars lon 201.0°). PATH: build corrected
  topology, seed G from App-C, n-body-validate; golden target = App-C real-eph
  v∞ (breathes 3.2–8.0, avg 5.48), NOT coplanar 5.10/3.05. Commit `8f24e1d`;
  note `docs/notes/2026-06-08-s1l1-source-dig.md`. → #167 builds it.
- **#167 — S1L1 CONFIRMED on DE440 (corrected topology), scoped.** The corrected
  E→g(E-E,sub-Mars)→E→G(E-M-E,longitude rendezvous)→E build, seeded from Russell
  App-C #83, lands ALL 7 Mars encounters in the 3-SOI band (same band #165 used,
  not loosened) at the published per-leg v∞ to 4 dp (breathing 3.2–8.0, avg 5.47,
  NOT coplanar 5.10/3.05); g-arcs sub-Mars (0.67–1.05 AU); holds under
  Mars-perturbed gravity. **#94's scientific question answered YES** — S1L1 IS a
  real ballistic E-M cycler; Russell's states reproduce on an independent
  integrator+ephemeris (DE405→DE440) with the right topology. SCOPE (no overclaim,
  cf. #164 retraction): the run RE-ANCHORS v∞ at each App-C node (Russell's per-leg
  recipe) — VERIFIES his published cycler independently, NOT yet a single
  continuous-from-one-seed multi-cycle closure. LEVEL: solid V3 leaning V4; clean
  V3/V4 upgrade = one continuous propagation measuring the re-anchor Δv as the
  maintenance budget. **V3 writeback recommended, HELD to user decision.**
  **STALE GATE CORRECTED 2026-07-15 (Fable adversarial pass, full-file audit): the writeback
  happened — `data/catalogue.yaml`'s `russell-ch4-4.991gG2` row carries `validation_level: V3`
  today (mechanically derived by `scripts/backfill_validation_level.py` from recorded gauntlet
  evidence, per the row's own schema comment). The #170 entry below already presupposes this
  ("S1L1 ... remains the only V3"). "HELD to user decision" no longer describes reality — the V3
  level is live in the catalogue now, not still pending.** Commits
  `cec1353`/`3ae4715`; note `docs/notes/2026-06-08-s1l1-corrected-closure-results.md`.
- **[SUPERSEDED] #164 — S1L1 CLOSES ON DE440 (#94 mechanism resolved).** Continuing the
  two-arc seed circular→real eccentric/inclined model
  (`search/continuation_chain.py`) closes `russell-ch4-4.991gG2` (the S1L1
  physical cycler, Russell-coplanar framing) on BOTH halves: vinf_res **0.006
  km/s** (4.99/5.10 own anchors) AND both descriptor ToFs inside the 0.1-yr band
  (g-arc gap 0.14→0.032, G-arc 0.053), bend-feasible. Real Mars eccentricity
  breaks the circular ToF quantization exactly as predicted (e-ramp monotone,
  epoch slides to J2000+239 d). Closed V3-seed elements recorded. **THREE
  caveats before V3 writeback (held to main session):** (1) closes at the row's
  OWN 4.99/5.10, NOT spec §9 5.65/3.05 (target was always framing-ambiguous);
  (2) ToF close-within-band, not exact (~12–19 d); (3) V3-class two-body model,
  single-solver, no independent cross-check yet. **Decisive confirmation =
  full N-encounter DE440 propagation of the closed geometry at J2000+239 d.**
  Commits `561f440`/`2e7ac69`. 6.44Gg3 stays a quantified partial (g-arc target
  below the n_rev=1 floor — structurally deeper).

## App-C V3 batch result — 0/9 (#170): S1L1 was singular, not a batch (2026-06-08)

The 9-parent App-C→V3 batch promoted **0 new V3 rows**; S1L1
(`russell-ch4-4.991gG2`) remains the only V3. NOT a pipeline failure — a clean,
honest negative that corrects the earlier "V3 1→~10" optimism:
- **7/9 NOT-REACHABLE** — coplanar Table-4.10-4.13 members with no sourced App-C
  reproduction block and no 1:1 coplanar→App-C bridge; forcing one would be
  unsourced inference (discipline forbids). Skipped.
- **2/9 reachable but PARTIAL** (`8.049gGf2` #188, `8.165Gfh-f2` #192): the
  CLOSURE half passes — all 7 Mars encounters per row reconstruct in-band on an
  independent REBOUND/IAS15 DE440 integrator at the published per-leg v∞ (4 dp),
  with true-longitude rendezvous — but they are **POWERED** cyclers (published
  App-C Δv 420 / 1678 m/s; continuous TCM 164 / 2041 m/s, both > the 120 m/s V3
  budget). They fail V3-ballistic on budget, not on geometry.
- **S1L1 is special**: it is the one parent that is BOTH reachable AND genuinely
  near-ballistic (62 m/s continuous TCM). The others are powered or unsourced.
- **Implications:** (a) the App-C V3 path is now exhausted — more V3s need the
  **self-seeding longitude-rendezvous construction** (reaches the 7 + the 194
  members, removes the App-C-block dependency); (b) the powered reachable parents
  raise a **V3-powered class-split** question (parallel to the §14 V2 split) — they
  have real, independently-confirmed encounter geometry but exceed the
  ballistic budget. Modules `search/appc_corrected.py`; commits `359b248`/`bf8d2a8`;
  note `docs/notes/2026-06-08-appc-v3-batch-results.md`.

## GMAT V4 + maintenance-ΔV reconciliation (#174/#176, 2026-06-08)

GMAT R2022a installed + headless-verified; it CONVERGED both Aldrin and S1L1
(geometry externally validated). The two maintenance-ΔV mismatches it flagged
were both ARTIFACTS, resolved without changing our claims:
- **Aldrin:** #174's GMAT 0.175 km/s was the WRONG FLYBY (Mars, not the Earth
  return-flyby turn deficit). Corrected GMAT Earth-leg = 3.287 km/s impulsive
  (confirms the asymptote class ~2.9–3.3). GMAT impulsive does NOT reach our
  1.9336 Oberth figure — the Oberth credit needs a CONTINUOUS widen-cone-restore
  strategy. Defensible: bound pair **2.9138 (conservative) / 1.9336 (Oberth
  lower bound, continuous-strategy)**; 0.175 retired. Site keeps 2.914 (correct).
- **S1L1:** #174's 7.29 km/s was the forced-v∞-vector artifact (manufactured
  Mars turns exceeding the cone; absent in continuous flight). **62 m/s (#169)
  STANDS; S1L1 V3-ballistic needs NO qualification.**
- **V4 held for both** (Aldrin ~13% over the ±5% band; S1L1 needs an
  operational-targeting GMAT setup). No writeback. Notes
  `docs/notes/2026-06-08-gmat-v4-execution-results.md` +
  `2026-06-08-maintenance-dv-reconciliation.md`.

## Self-seeding construction VALIDATED (#173, 2026-06-08)

The App-C-blind self-seeding longitude-rendezvous construction PASSED its S1L1
ground-truth gate: given only the descriptor (never reading App-C), it recovered
S1L1's App-C basin as a SINGLE candidate (Mars lon 199.5° vs key 201.0°; n-body
arrival 1.7e-12 AU at v∞ 5.42) — beating the 2026-06-04 off-family failure via
corrected topology + explicit longitude target + enumerate-don't-optimise.
QUALIFIED: ~11-day epoch offset (family-correct WINDOW, not flight-grade).
Prove-on-one (6.44Gg3): clean OFF-FAMILY (coplanar 131-d transit ≠ 262-d real-eph
family; emerged v∞_M 10.9 vs 3.74). **The method is a TRIAGE + family-seed tool,
not a one-shot V3 machine** — works when coplanar transit ≈ real-eph transit.
Path to volume EXISTS but is a minority-yield triage, not "194→V3". Next (#177):
triage the 7+194 (OFF-FAMILY → negative-results registry); multi-rev G-arc Stage-A
extension for long-transit rows; REACHABLE → full pipeline → V3 candidates.
Module `search/self_seeding.py`; notes `docs/notes/2026-06-08-self-seeding-results.md`.

## Self-seeding triage RESULT (#177) + Forge Phase 6 first sweep (#172) — both EMPTY-ish, 2026-06-09

- **#177 triage of all 212 unsourced rows: 6 REACHABLE, 206 OFF-FAMILY, 0 V3-candidates.**
  Built the multi-rev G-arc Stage-A extension (short/long/k-rev branches). 204 rows
  OFF-FAMILY-NO-DESCRIPTOR (all ~194 ocampo have no g/G arc ToFs / are sub-Mars);
  the 6 reachable ch4 rows achieve longitude rendezvous (n-body miss at machine
  precision) but MISS their sourced v∞ anchors (OFF-FAMILY-AT-ANCHOR-VINF). The
  multi-rev extension pulled 6.44Gg3 v∞_M 10.9→7.83 (right direction) but not to the
  3.74 anchor. **S1L1 is SINGULAR** — the realistic V3 ceiling via current
  methods+sourced data is V3:2 (S1L1 ballistic + #188 powered); more V3s need new
  sourced real-eph data (#116), not more method. Runlogs (method-versioned) →
  #172 registry. Note `docs/notes/2026-06-08-self-seeding-triage-results.md`.
- **#172 Forge Phase 6 first novelty sweep (Jovian I-E-G): EMPTY** — the expected
  base rate, a complete valid result. 128 points, 12 closed + all `novel` (sparse
  bucket) but **0 bend-feasible** (V∞ 8.3–26.8 vs ~6 floor, gap 20.8 km/s) → 0
  SILVER, 12 REJECTED. 1 method-versioned empty region recorded
  (`data/empty_regions.jsonl`); re-sweep gate + dedup firewall verified LIVE; no
  writeback, no tolerance loosening. A bend-feasible Jovian cycler needs a
  strictly-more-capable genome (the §6b gate auto-re-sweeps when one ships). Built
  the full Phase 6 pipeline + capability-subsumption registry (16/16 plan tasks).

## Done

- **Aldrin cycler replicated on the real ephemeris** — M6b binding gate;
  patched-conic constructor reproduces the literature ellipse to ±0.002
  (see Q&A item A, commit `ba12554`).
- **Canonical 2-synodic Earth–Mars anchor replicated idealised**
  (circular-coplanar model).
- **Multi-revolution Lambert solver landed** — universal-variable
  multi-rev solver with sourced cross-check, single-rev Newton restored,
  multi-rev threaded through `verify` (commits `7d620e8`, `9f618f3`).
- **S1L1 topology corrected** — the catalogue no longer carries a
  spurious direct Mars→Earth "return leg" (commit `f4074d0`). See
  Key findings below and the S1L1-nomenclature memory note.
- **S1L1 not hostable in the idealised model — characterised and
  documented** — `scripts/characterise_s1l1.py` shows no circular-coplanar
  topology reproduces the published 5.65 / 3.05 km/s V∞ anchors within
  ±0.3 km/s. The E-M-E-E topology was additionally tested
  (`scripts/characterise_s1l1_emee.py`) and does NOT close in
  circular-coplanar (V∞_E ≈ 25-39 km/s): the blocker is the MODEL (the
  ~154-d outbound leg is near-hyperbolic, exactly like the Aldrin case),
  NOT the topology. Real-ephemeris closure additionally needs multi-rev
  support in the maintenance engine. S1L1 gates remain xfail. This is a
  model limitation, not a solver bug.
- **Hollister & Menning Earth-Venus family individuated** — the single
  E-V placeholder was expanded into the 15-orbit
  `hollister-menning-1970-ev-orbit-01..15` family from the PRIMARY
  Hollister & Menning 1970 paper (now in `docs/refs/`), with V∞ from
  Table 3 (Vr × 29.785 EMOS) and a shared period 16 yr / k=10 (corrected
  from a wrong secondary "3.2 yr", which was the coplanar sub-orbit).
  New V0 data-integrity tests: `tests/data/test_hollister_family.py`,
  `tests/data/test_aldrin_establishment.py`.
- **Real-ephemeris cell optimiser IMPLEMENTED** — `optimise_cell_ephemeris`
  (was a `NotImplementedError` stub) now runs against DE440 with an
  asymmetric `tof_seed_days` option and an Aldrin parity test; `discover()`
  supports `optimiser="ephemeris"` (Part 2 wired).
- **Catalogue census frozen as a ratchet** — `tests/test_catalogue_rediscovery.py`
  (`EXPECTED_COVERAGE`) pins the 237-entry distribution; any catalogue
  change must update it in the same commit.
- **Data-validation hardening (Forge phase 0) — COMPLETE** —
  `docs/superpowers/plans/2026-06-03-data-validation-hardening.md` is now
  fully shipped. The 2026-06-05 reconciliation deferred three items (Tasks 3,
  5, and Task 4's live-row census ratchet); all three landed this run.
  Task 3: per-field provenance tags (`orbit_source` / `vinf_source` /
  `orbit_fidelity` / `vinf_fidelity`) back-filled mechanically onto 224 of 237
  rows by `scripts/backfill_provenance_tags.py` (the 13 untagged rows name no
  `SOURCE_REGISTRY` paper — explicit "unknown" marker, never a guessed key);
  the forward-compatible `validate_provenance_tags` gate now validates real
  data. Task 4 remainder: `tests/data/test_validation_tier_census.py` freezes
  the live tier distribution `{cross_validated: 5, consistency_checked: 218,
  unvalidated: 14}` as a monotone ratchet. Task 5: `score_corroboration`
  (`tests/data/test_corroboration.py`) classifies a quantity strongly-sourced /
  single-sourced / disputed and surfaces the spread (the S1L1 5.65-vs-4.99
  Earth-V∞ dispute is the documented disputed case). Schema bumped 4.3 → 4.4
  (the four tag fields as enums; spec §16.7.11). The `cycler_class` and
  `EXPECTED_COVERAGE` census ratchets are unaffected (additive metadata).
- **Moon-tour Tier-1 (patched-conic moon systems + VILM) — SHIPPED** (task #76,
  `docs/superpowers/plans/2026-06-06-moontour-tier1-patched-conic.md`). The
  catalogue's patched-conic moon-system rows are now computable on the same
  Kepler-conic + impulsive-V∞-rematch model as the heliocentric catalogue, with
  central body = a planet and flyby bodies = its moons. Delivered:
  - **`SATELLITES`/`PRIMARIES` registry** (`core/satellites.py`) — Galilean four
    + Saturnian midsize + Titan, JPL-SSD-sourced, mean motion derived at import
    via Kepler III; registry-construction golden reproduces the published Endgame
    Part-1 Table 3 ã_M / Ṽ_M independently.
  - **Planet-centred circular ephemeris** (`_CentredCircularBackend`,
    `Ephemeris(center=...)`) — moon states about the primary, km-scaled;
    heliocentric backends byte-identical.
  - **Centre-agnostic corrector** — `mu_central` plumbed into Lambert and
    `_max_bend_deg` resolves moon codes via `SATELLITES`; Sun-default keeps the
    heliocentric solver byte-identical. The Io-Europa-Ganymede chain CLOSES about
    Jupiter.
  - **Centre-aware Tisserand** — `_a_p_km` + `mu=` resolve a moon; `T = 3 − v∞²`
    round-trips about Jupiter; a Jovicentric pair prunes through `linkable`.
  - **VILM module** (`search/vilm.py`) — Eq.(9) V̄∞-efficiency root (vs Table 3
    E/I), Eq.(13) ΔV-min quadrature (vs Table 1/2), Europa 3-VILM endgame scalar
    (vs A6 154 m/s / 46 d), GA routing, and an admissible ΔV-floor for search
    pruning.
  - **(model_assumption, primary) pool pre-filter** — a Jovicentric V∞ never
    compares to a heliocentric one; heliocentric signatures byte-identical.
  - **Gauntlet/fidelity integration** — Axis-B persistence about a primary,
    Axis-A VILM-vs-Lambert agreement, Axis-D wrong-central-μ falsification guard.
  - **Catalogue re-tag** — the two Jovian rows `non-keplerian` → `multi-arc` (so
    the gauntlet routes them to invariants, not CR3BP); the Saturnian row keeps
    `non-keplerian` with an honest Titan-Tier-1 / midsize-Tier-2 split note.
  - **Deviations recorded (honest-risk):** (1) the Phase-3 I-E-G closure is
    bend-INFEASIBLE in the coplanar-circular no-V∞-leveraging model (closes at
    ~10 km/s V∞ needing 100–150° turns vs 2–5° max-bend) — recorded as a strict
    xfail, NOT forced; bend-feasible Jovian tours need the VILM layer +
    Laplace-resonance phasing. (2) The VILM admissible ΔV-floor was implemented
    as escape+capture, not the no-GA quadrature (the plan's "no-GA-as-floor" is
    backwards — a gravity assist *reduces* ΔV, so the no-GA value is not a lower
    bound). (3) The Axis-A crosscheck was reframed to VILM-vs-Lambert on a
    Jovicentric Hohmann V∞ (both independent code paths) rather than
    VILM-vs-corrector, because the corrector's closed family is the higher-V∞
    non-bend-feasible one. The corrector's V∞ is NON-GOLDEN throughout — the
    Russell-Strange/Hernandez rows are family-seed null-numeric records with no
    sourced Jovicentric V∞ multiset; only the VILM ΔV (Part-1 Tables) and the
    registry construction (Table 3) are golden-gated.
  - **Census:** catalogue is 268 rows including the moon rows (3 Earth-primary,
    2 Jupiter, 1 Saturn); the MULTI_ENCOUNTER_SEQUENCE ratchet moved 192 → 223
    from the Russell 2004 + Rall 1970 ingests (not from the moon-tour re-tag,
    which added no rows).
  - **Tier-2 (CR3BP) remains OPEN** — Earth-Moon Arenstorf/Genova/Wittal and the
    Saturnian midsize-moon (Mimas/Enceladus/Tethys) members stay citation-only
    until the later CR3BP milestone.

## In progress (with plan references)

- **M-ED real-ephemeris ballistic corrector** —
  `docs/superpowers/plans/2026-06-05-m-ed-ballistic-corrector.md`. **Phases 1-5
  shipped.** The N-arc ballistic differential corrector (`search/correct.py`,
  `ballistic_correct`), the Russell descriptor parser (`search/descriptor.py`),
  the default-inert `mode="ballistic"` on `optimise_cell_ephemeris` (real
  V∞-continuity closure residual, maintenance mode byte-identical), the seeding
  ladder (`search/seed_ladder.py`, descriptor → anchor → coplanar → scan), and a
  real ballistic-closure V3 gate in `data/discover.py` (replacing the stale
  `# raises until M6b lands` stub branch) are all landed and gated. The corrector
  closes the S1L1 two-arc chain and a Sanchez-regime near-ballistic chain on
  DE440.
  - **OPEN RESEARCH — Jones VEM headline gate does NOT converge, even with a
    dense parallel scan (task #110, 2026-06-06).** The headline rediscovery gate
    (`tests/test_vem_rediscovery.py::test_jones_vem_ballistic_rediscovers_
    sourced_multiset`) stays **xfail** (STOP/report branch, plan Task 5.4 —
    tolerance NOT loosened, xfail NOT flipped). The parallel epoch×branch scan
    engine now exists (`search/scan.py`, `scan_parallel` over a `ProcessPool`;
    measured 3.98× speedup on a 32-point DE440 grid, 8 workers) and drives the
    ballistic mode via `optimise_cell_ephemeris(mode="ballistic", scan_epochs=N)`
    (`_ballistic_scan_rung`). A DENSE hunt (`scripts/hunt_vem_ballistic.py`:
    256 epochs over the full 12.8-yr repeat period × 11 rev/branch topologies =
    **2816 points/row**, 16-core) STILL FAILS:
    - **EMEVVE outbound** (E-M-E-V-V-E): 831 closed / 474 distinct families;
      BEST max-V∞ **17.86 km/s** (per-encounter [13.88, 13.91, 15.39, 15.39,
      17.43, 17.86]); **0 bend-feasible**. Sourced 2.5–7.0.
    - **MEEVEM inbound** (M-E-E-V-E-M): 1239 closed / 570 distinct families;
      BEST max-V∞ **18.49 km/s** (per-encounter [11.40, 16.34, 16.34, 17.83,
      18.47, 18.49]); **0 bend-feasible**. Sourced 2.42–5.16.

    No bend-feasible closed family below 10 km/s exists in either survey. The
    closed families floor ~11–18 km/s, far above the sourced Jones 2.42–7.0 —
    the documented S1L1 Mars-V∞ ~6.4 floor generalised: the single-ellipse-per-
    leg corrector closes a *different, higher-V∞, powered* family than the Jones
    members. **A denser-scan-still-fails result is real science: it sharpens the
    hypothesis** — reaching the Jones VEM family is NOT a scan-density problem;
    it needs 3D inclination (M-3D), real-eccentricity intermediate flybys, or a
    different (e.g. multi-arc-per-leg) topology seeding. (Compare the S1L1
    multi-arc finding in project memory `project_s1l1_realeph_closure_blocker`.)
    The future McConaghy 4.7/5.0 anchor remains the cross-check target for the
    S1L1 floor.
  - **3D-inclination hypothesis REFUTED for this family (task #120,
    2026-06-06).** The first sharpened hypothesis above was tested directly:
    the identical grid (256 epochs × 11 rev/branch topologies = 2816 points/row,
    both Jones rows, 16-core) re-run on the M-3D inclined-circular backend
    (`Ephemeris(model="inclined-circular")` — sourced J2000 Venus i=3.39°/Mars
    i=1.85°, via `scripts/hunt_vem_ballistic.py 256 inclined-circular`). The
    floors did **NOT** drop materially toward the sourced 2.4–7.0:
    - **EMEVVE outbound**: 1496 closed / 597 distinct families; BEST max-V∞
      **17.49 km/s** (per-encounter [12.49, 13.18, 13.18, 16.31, 16.31, 17.49]);
      **0 bend-feasible**. (Coplanar #110 was 17.86.) Sourced 2.5–7.0.
    - **MEEVEM inbound**: 1257 closed / 415 distinct families; BEST max-V∞
      **18.39 km/s** (per-encounter [11.30, 15.62, 15.62, 18.38, 18.39, 18.39]);
      **0 bend-feasible**. (Coplanar #110 was 18.49.) Sourced 2.42–5.16.

    Inclination moved the EMEVVE floor by **−0.37 km/s** and MEEVEM by
    **−0.10 km/s** — noise against the ~11–14 km/s gap to the sourced multiset,
    and STILL zero bend-feasible. A small **DE440 control** (full 3D + eccentric;
    64 epochs × 11 topologies on EMEVVE) bounds this: BEST max-V∞ **18.16 km/s**,
    0 bend-feasible — i.e. even real eccentricity+inclination together do not
    reach the family with this corrector. **Conclusion: 3D-inclination-only is
    REFUTED as the direction for the Jones VEM family.** The blocker is the
    single-ellipse-per-leg topology itself (it closes a different, higher-V∞,
    powered family); the front-runner moves to **multi-arc-per-leg seeding /
    eccentric intermediate flybys** (cf. the S1L1 multi-arc finding,
    `project_s1l1_realeph_closure_blocker`). The headline gate
    (`test_jones_vem_ballistic_rediscovers_sourced_multiset`) therefore stays
    **xfail** — tolerance NOT loosened, xfail NOT flipped.
  - **Vector-residual (Jones-method) hypothesis REFUTED for this family (task
    #122, 2026-06-06).** The design `docs/superpowers/specs/2026-06-06-jones-
    family-corrector-variant-design.md` argued the #110/#120 floor is a *residual*
    defect: magnitude-only continuity selects a powered, over-bent basin, and a
    full v∞-VECTOR residual (bend-feasibility hinge INSIDE the least-squares)
    should steer toward the bend-feasible Jones family. **Phase 1 — the cheap
    falsifier — implemented and run** (`search/correct.py` `residual_mode=
    "vector"`, threaded through `search/scan.py`; same #110 grid, 256 epochs ×
    11 topos = 2816 points/row, astropy DE440, 16-core, via
    `scripts/hunt_vem_ballistic.py 256 astropy vector`):
    - **EMEVVE outbound**: **0 closed / 2816**, 0 distinct families.
    - **MEEVEM inbound**: **2 closed / 2816**, 1 distinct bend-feasible family at
      BEST max-V∞ **20.78 km/s** (per-encounter [15.2, 17.45, 17.45, 19.82,
      19.83, 20.78]). Sourced 2.42–5.16.

    The in-residual hinge *confirms the design's mechanism* — the #110 magnitude
    basin (831/1239 closed at 17.86/18.49) WAS the powered/over-bent basin: it
    collapses to ~0 closures once feasibility is in the residual. But **no
    bend-feasible family appears below 10 km/s** (the Phase-1 gate). The one that
    closes sits at ~21 km/s, ~4× the sourced floor.
  - **Phase 2 (per-flyby B-plane angle θ_B free vars) — STOP/report, NOT
    attempted as a code change.** Decisive physics: in this Lambert-chain
    corrector each flyby's v∞-in AND v∞-out are *both* pinned by the heliocentric
    Lambert geometry between fixed bodies/epochs; a B-plane angle is an *output*,
    not a free input. The gap to the sourced multiset is a **|v∞|-MAGNITUDE gap**
    (closed legs 15–21 km/s vs sourced 2.4–7.0), and θ_B re-aims the bend at
    *constant* |v∞| — it cannot lower a leg's magnitude. A *true* θ_B DOF would
    require replacing the Lambert-both-sides chain with a flyby-PROPAGATION
    shooter (choose the bend at sourced r_p, propagate v∞-out forward) — that is
    the deferred Phase-3 (n-body shooting) architecture, out of v1 scope. The
    sourced periapsis altitudes exist (`flyby_mechanics[].min_altitude_km`,
    Jones Tables 2/3: E 100/Mars 4164/E 3814/V 684/… km) and would seed r_p, but
    seeding r_p constrains bend *magnitude*, not leg |v∞|. **Conclusion: the
    vector residual + B-plane targeting (design direction (c)) is REFUTED as a
    patched-conic Lambert-chain corrector lever for the Jones VEM family** — the
    |v∞|-magnitude family selection is a leg-topology problem, exactly the
    multi-arc-per-leg / Phase-3 shooting front-runner. The headline gate stays
    **xfail** — tolerance NOT loosened, xfail NOT flipped.
  - **Free-return (radial-crossing) genome — Russell-12 BREAKTHROUGH (#137,
    2026-06-07).** A new single-ellipse free-return genome
    (`search/free_return.py`: per-body V_inf EMERGES from the ellipse shape
    `(a, e)` riding the Mars-V_inf ridge; sourced geometry is now a residual-zero
    point) closes **8 of the Russell-12 rows CLOSE-AND-MATCH** + 3 SYMMETRIC-ONLY +
    1 OFF-ANCHOR (was 0 on the old Lambert genome). See
    `docs/notes/2026-06-07-russell12-freereturn-results.md`.
    - **Part 1 — §14 V1 writeback (DONE).** The literal §14 V1 mechanics
      (lamberthub izzo+gooding leg agreement < V1_TOLERANCE_MPS AND Kepler forward
      re-propagation < KEPLER_REPROP_TOL_KM; reused verbatim from
      `verify.agreement` paths a+c) applied LIKE-FOR-LIKE on the circular ephemeris
      to each closure (`search/free_return_v1.py`,
      `tests/search/test_free_return_v1_mechanics.py`). A V_inf-continuity honesty
      gate splits the rows: **four** whose single ellipse forms a genuinely closed,
      V_inf-continuous E->M->E cycler (`russell-ch4-5.30gGf3`, `-9.94Gg3`,
      `-5.75ggF3`, `-9.353Gg2`) pass → promoted V0→V1; the six genuinely multi-arc
      rows (forced return breaks Mars V_inf continuity by ~24 km/s) and the
      Lambert-singular `8.049gGf2` are refused. **Validation census: 6 rows V1**
      (Aldrin pair + these four). Scope: circular-coplanar reproduction of a
      circular-coplanar source (like-for-like), NOT real-ephemeris (V3).
    - **Part 2 — Jones VEM extension assessment (STOP/report).** The free-return
      primitive (one ellipse crossing TWO body radii) does NOT modestly extend to
      the Jones EMEVVE/MEEVEM topologies: they have multiple distinct transfer
      ellipses coupled through a bend-feasible Venus flyby, plus same-body resonant
      legs (V->V, E->E) the radial-crossing primitive cannot represent. A control
      probe confirms the E-M sub-arc reaches the Jones Mars-V_inf (~2.81 floor), so
      the blocker is purely the multi-ellipse-through-Venus-flyby structure — i.e.
      the existing multi-arc-per-leg / flyby-shooter front-runner, NOT a modest
      `free_return.py` addition. No Jones hunt run on a half-genome; the headline
      Jones xfail stays xfail. Design questions recorded in the results note Part 2.
    - **Part 3 — stragglers + probe fix (DONE).** `--probe-at-truth` now honours
      `--genome` (`probe_at_truth_free_return`). `9.353Gg2` was OFF-ANCHOR purely
      from a coarse 256-point phase grid (narrow high-e t0 basin) → fixed by a dense
      phase floor (`FR_PHASE_EPOCHS_FLOOR=4096`, cheap Lambert-free scan) → now
      CLOSE-AND-MATCH and §14-V1. `6.44Gg3` closes at 4096 but OFF-ANCHOR: its
      sourced aphelion+transit (a=1.27, 226-262 d) and sourced V_inf (needs
      a=1.225, 166 d) describe DIFFERENT free-return ellipses — a genuinely
      different family / descriptor-interpretation mismatch, not a phase deficiency;
      no promotion. Candidate follow-up: seed `6.44Gg3` from its V_inf rather than
      aphelion+transit.
- **SnLm multi-rev rediscovery** —
  `docs/superpowers/plans/2026-06-02-snlm-multirev-rediscovery.md`.
  Phase 1 Task 1 EXECUTED (S1L1 characterisation above). Phase 2 —
  plumbing multi-rev enumeration through `data/discover.py` and the
  sweep gate so the 202-entry `MULTI_ENCOUNTER_SEQUENCE` bucket becomes
  reachable by `discover` — IN PROGRESS.
- **Catalogue data-gap sourcing workflow** — a workflow is extracting
  Russell 2004 + Rogers 2012 tables and filling sourced catalogue gaps
  (goal: more usable golden orbits). Tracked in `data/MISSING_DATA.md`.

Two further plan files exist alongside the above and feed the same
effort: `2026-06-02-ml-multirev-lambert.md` and
`2026-06-02-multirev-3d-vem-ephemeris-roadmap.md`.

## Deferred (by explicit user decision — not being implemented now)

- ~~**M8-Core (VEM 3-body)**~~ — **EXECUTED 2026-06-05** (user un-deferred;
  commits `933e75b`..`eb851a2`): `Cell.period_basis` + beat dispatch in
  `_target_period_sec`, same-body Tisserand bypass, VEM rediscovery gate
  anchored on the Jones AAS 17-577 member rows (sourced 12.8 yr), all four
  VEM rows admitted as `CONSTRUCTIBLE_MULTIBODY` (`NOT_TWO_BODY` now 0).
  The data blocker (Q&A item D) was resolved the same day. Remaining
  M8-adjacent work: **M8-UX** (CLI/viz/reporting — still deferred) and the
  ballistic-convergence xfail handed to **M-ED**.

## Blocked — needs HUMAN institutional / paywalled access (prioritised)

These are the highest-value human-actionable items: each unblocks
catalogue gaps that no amount of computation or open-web fetching can
fill. Listed in rough priority order.

1. **Russell 2004 dissertation tables (dominant gap source).** The PDF
   itself is open-access (UT Austin `http://hdl.handle.net/2152/1253`;
   full text held offline, not stored in repo) but binary-compressed —
   the in-progress extraction workflow is handling it. 772 of 790 gaps
   trace here. NOT human-blocked, but flagged as the single
   highest-leverage target; see `data/MISSING_DATA.md` §3.1.
2. **Friedlander/Niehoff/Byrnes/Longuski 1986, AIAA 86-2009-CP** —
   VISIT-1 / VISIT-2 V∞ at Earth and Mars. Paywalled
   (DOI `10.2514/6.1986-2009`). 4 "unknown" V∞ gaps; never digitised
   online (Q&A item C still partly open). Partial mitigation 2026-06-10:
   McConaghy 2004 Table 5.5 identifies VISIT 2 = r=10 / VISIT 1 = r=12
   with PARAMETRIC aphelion/period ranges (perihelion free), not exact
   V∞ — so the paywalled paper is still wanted for the exact values.
3. **✓ Resolved (2026-06-10)** — the McConaghy Purdue dissertation (2004,
   UMI 3166673) is acquired and mined
   (`docs/notes/2026-06-10-mcconaghy-2004-dissertation-mining.md`): U0L1
   steady-state V∞ resolved from Table 5.8 (U0L1(2.7540) is ballistic,
   ΔV 0.00, with per-leg Earth V∞ 11.3/11.3 and Mars V∞ 14.0/5.4 km/s;
   aphelion 3.20/1.54 AU, period 2.93/1.18 yr), and "Case 1" resolved
   from Table 5.4 (2L3 = Byrnes Case 1: V∞E 5.65, V∞M 3.05, simple
   model). Original text:
   **McConaghy 2005 Purdue dissertation (AAI3166673) / AIAA 2002-4420
   full text** — U0L1 and Case 1 steady-state V∞ (and U0L1 return ToF).
   Both ResearchGate / ProQuest restricted. Note: the McConaghy 2006
   S1L1 orbital-element gap this would have closed is now resolved via
   Russell Table 4.9 (Q&A item B), so this is needed only for the
   remaining U0L1 / Case 1 / SnLm steady-state V∞ gaps. **NB S1L1,
   U0L1, and "Case 1" are THREE DISTINCT orbits** (Rogers 2012 Table 1:
   S1L1 a=1.30/e=0.257, Case 1 a=1.22/e=0.238, U0L1 a=2.05/e=0.563),
   not the same trajectory.
4. **Jones / Hernandez / Jesick 2017 (AAS 17-577)** — VEM triple-cycler
   member list (Q&A item D). Needed to lift M8 out of placeholder
   territory; lower priority while M8 is deferred.
5. **Data-architecture design references (spec §16.7.8) — mostly resolved
   2026-06-05 (structure-only, no golden anchors):**
   - ~~CCSDS 502.0-B-3 *Orbit Data Messages*~~ — **obtained and mined
     2026-06-05** (`docs/notes/2026-06-05-ccsds-odm-502-mining.md`).
     Yield: OCM lives in 502.0-B-3 §6 (spec's "504.0-B" citation was wrong
     and is fixed); per-TRAJ-block CENTER_NAME validates v4.2
     `segments[].center`; `CELESTIAL_SOURCE` validates `source_ephemeris`.
   - Acton 1996, "Ancillary data services of NASA's NAIF," *Planet. Space
     Sci.* 44(1):65–70, DOI `10.1016/0032-0633(95)00107-7` — **still
     paywalled**; cited from bibliographic record only.
   - ~~Campagnola & Russell 2010, Endgame Problem Parts 1+2, *JGCD* 33(2),
     DOI `10.2514/1.44258` / `10.2514/1.44290`~~ — **obtained (AAS
     09-224/09-227 preprints) and mined 2026-06-05**
     (`docs/notes/2026-06-05-endgame-tisserand-mining.md`). Yield: VILM
     ΔV-min quadrature (Forge admissible heuristic), branch & bound over
     the leveraging graph, T-P graph as CR3BP twin of our Tisserand module,
     candidate VILM golden anchors (Tables 1–3 Pt 1, Table 1 Pt 2) for a
     future VILM module.
   - ~~Campagnola, Skelton & Lantoine 2014, "Global Search for
     Gravity-Assist Trajectories"~~ — **resolved 2026-06-05: no such paper
     exists** (verified against AIAA ARC / NTRS / IEEE Xplore). Likely a
     mix-up of Anderson, Campagnola & Lantoine 2015, *Cel. Mech. Dyn.
     Astron.* 124(2):177–199, DOI `10.1007/s10569-015-9659-7`. Dropped as
     a design reference (spec §16.7.8 item 3 now rests on Russell 2004
     alone, which is held and verified).
   - Also scanned 2026-06-05 (`docs/notes/2026-06-05-vasile-hiraiwa-scan.md`):
     Vasile & Campagnola 2009 (JBIS, arXiv:1105.1823) — low-thrust MGA to
     Europa, NOT a cycler; #76-relevant synchronous moon tours + #37
     candidate anchors (Tables 3–4, need re-transcription from a clean
     copy: broken fonts). Hiraiwa et al. 2026 (arXiv:2602.17444) — cislunar
     lobe-dynamics transfers, NOT a cycler; Table 3 carries
     literature-sourced Earth-Moon ΔV goldens (Hohmann 3954 m/s, Sweetser
     min 3726 m/s, traceable to Topputo 2013).

## Catalogue census (frozen ratchet — 237 entries)

Source of truth: `tests/test_catalogue_rediscovery.py` `EXPECTED_COVERAGE`.

| Exclusion reason | Count |
|---|---|
| `multi_encounter_sequence` | 204 |
| `missing_leg_tofs` | 15 |
| `non_heliocentric` | 6 |
| `missing_vinf` | 5 |
| `not_two_body` | 4 |
| `constructible` | 2 |
| `missing_period` | 1 |
| **Total** | **237** |

(2026-06-05: `not_two_body` 2 → 4 with the two Jones 2017 VEM triple-cycler
members — each spans three bodies E/M/V. cycler_class census separately:
multi-arc 201 → 203, total 235 → 237.)

The new `missing_leg_tofs` bucket (15) holds the individuated Hollister &
Menning E-V family (orbits 1-15): each has matched V∞ at Earth and Venus
but no per-leg ToFs encoded in `legs`. Expanding the old single E-V
placeholder moved one row out of `missing_vinf` (6 → 5) and one out of
`missing_period` (2 → 1), and added the 15-row bucket.

Through the idealised rediscovery gauntlet ~0/237 currently pass green:
the 2 `constructible` (Aldrin) entries are model-limited skips
(`EXPECTED_SKIPS`). Aldrin IS replicated on the real ephemeris (M6b)
and the canonical 2-synodic E–M anchor is replicated idealised. The
dominant blocker is the 202-entry (87%) `multi_encounter_sequence` /
SnLm multi-rev modelling bucket — exactly what the in-progress SnLm
multi-rev plan targets.

## Data-gap accounting (see `data/MISSING_DATA.md`)

788 data gaps across 207 entries:

- **403 "unknown"** — need sourcing (the bulk from Russell 2004, the
  rest from Rogers 2012 / McConaghy 2005 / VISIT).
- **201 "derive"** — COMPUTABLE by our Lambert solver, no sourcing
  needed.
- **184 "uncertain"** — topology-provisional; resolve as a by-product
  of extracting Russell per-leg elements.

The bulk of the gaps trace to the Russell 2004 dissertation. Key reference
full texts are held offline (not stored in repo); cite by DOI or handle:
Russell 2004 (`hdl:2152/1253`), Rogers 2012 (DOI `10.2514/6.2012-4746`),
Genova 2016 (AAS-15), Pascarella/Sanchez Net (see `docs/v2-future-references.md` §2
and `data/OUTSTANDING.md` H.3 / task #38 for the two distinct Pony Express papers).

## Key findings

- **S1L1 nomenclature.** In McConaghy / Longuski / Byrnes, `S` and `L`
  denote consecutive Earth-to-Earth RESONANT INTERVALS, not Earth↔Mars
  legs. S1L1 therefore has NO direct Mars→Earth return leg — each
  vehicle is one-way; the crewed return needs a mirrored conjugate
  L1S1 cycler. The catalogue's old direct M→E "return leg" was a
  mismodelling, corrected in commit `f4074d0`.
- **S1L1, U0L1, and "Case 1" are three distinct orbits** — per Rogers
  2012 Table 1: S1L1 a=1.30/e=0.257, Case 1 a=1.22/e=0.238, U0L1
  a=2.05/e=0.563. No doc or entry should treat them as the same
  trajectory.
- **Idealised-model limitations are real, not bugs.** Both the Aldrin
  cycler and S1L1 have published anchors that the circular-coplanar
  idealised model cannot host. These are documented xfails / skips, and
  the now-implemented real-ephemeris cell optimiser is the intended way
  to close them.
- **Residual true data gaps** (need human / restricted access):
  VISIT-1/2 V∞ (Friedlander/Niehoff 1986 AIAA 86-2009-CP); U0L1 &
  Case 1 V∞ (McConaghy 2005 dissertation / AIAA 2002-4420 full text,
  ResearchGate/ProQuest restricted). The Aldrin establishment Mars-V∞
  is N/A (Earth-side V∞ leveraging). The 201 catalogue `derive`-kind
  data_gaps are computable (no sourcing). New reference PDFs added to
  `docs/refs/` this session: russell-2004-dissertation, rogers-2012,
  hollister-menning-1970 primary scan, hollister-rall-1970 NASA-CR,
  landau-longuski-2006-pt1, vasile-2005, genova-2016-phobos-deimos.

## Infra note

Background subagents in this environment appear to lack the `Bash` tool
(under test). Implementation and verification (pytest, git) are run in
the foreground by the orchestrator, not by dispatched subagents.

---

## ✓ Resolved (2026-05-31) — A. Aldrin orbital-element discrepancy

**Resolution:** commit `32c5eab` (errata: reconcile spec §3/§9/§9.1/§16.4)
plus the M3 implementation in `ba12554` which reproduces the literature
ellipse to ±0.002. See `docs/errata-investigation.md` §1 for the full
analysis: the spec's `a=1.659, e=0.41` is a resonance-construction
choice that is internally inconsistent with the same spec's 146-d
Earth-Mars leg (those elements imply 138.9 d). The literature value
`a=1.60, e=0.393` is internally consistent at 146 d and is what the
M3 patched-conic constructor produces. spec.md §9 now carries the
literature values; §9.1 documents the reconciliation.

**Original question (preserved as audit trail):**

`spec.md` §9 anchored the M3 gate to:

> Aldrin cycler: a ≈ 1.659 AU, e ≈ 0.41, perihelion ≈ 0.98 AU,
> aphelion ≈ 2.34 AU, E→M leg ≈ 146 d.

But the literature consistently reports:

> a = 1.60 AU, e = 0.393, perihelion = 0.97 AU, aphelion = 2.23 AU,
> E→M = 146 d
> (Rogers et al. 2012 Table 1; Russell 2004 Table 3.4 via Aphelion
> Ratio 1.47; Wikipedia citing Byrnes/Longuski/Aldrin 1993)

The gap on *a* is 0.06 AU — six times the M3 gate's `TOL_A_AU = 0.01`.
Either the spec's value set is wrong, the literature set is wrong, or
both refer to different "Aldrin cyclers".

---

## ✓ Resolved (2026-06-01) — B. McConaghy 2006 orbital elements (medium priority)

**Fully resolved by Russell 2004 Chapter 4 tables ingest** (commit
pending, 2026-06-01). Russell 2004 dissertation Table 4.9 (page 127)
row 1 carries the orbital data for the McConaghy 2006 "Notable" (S1L1)
cycler under Russell's own nomenclature `4.991gG2`:

- aphelion = 1.64 AU (row 1 of Table 4.9)
- V_inf at Earth = 4.99 km/s (Russell), vs McConaghy 2006 abstract = 4.7 km/s
- V_inf at Mars = 5.10 km/s (Russell), vs McConaghy 2006 abstract = 5.0 km/s
- E-M ToF = 150 days (Russell), vs McConaghy 2006 abstract = 153 days

Russell explicitly cross-references the two: dissertation line 7416
states "cycler 4.991gG2(#83) ... Also known as the 'S1L1' cycler",
line 5476 says "notable 'S1L1' cycler ... discovered first by
McConaghy et al. in Ref. 15", and line 8008 lists "the S1L1 cycler
(4.991Gg2), 8.049gGf2, and the Aldrin cycler" as the three most
promising designs.

**Entry 2 (`mcconaghy-2006-em-k2`) updated:**

- `orbit_elements.aphelion_au` backfilled with 1.64 AU.
- `orbit_elements.note` updated to cite Russell Table 4.9 row 1 as
  the source.
- `orbit_elements.a_au`, `e`, `perihelion_au` retain `null` because
  the cycler is a piecewise sequence of two generic-return arcs
  (g(1.4612,526.02,Ll) + G(2.8096,651.46,U)), not a single Keplerian
  ellipse. Each leg has its own (a, e), so the whole-cycler (a, e)
  are not well-defined; only the maximum aphelion is.
- `source_quotes.orbit_elements.aphelion_au` added with Russell Table
  4.9 citation.
- `notes:` block updated with the V_inf discrepancy analysis
  (McConaghy 4.7/5.0 vs Russell 4.99/5.10) and points to the
  new sibling entry `russell-ch4-4.991gG2` carrying Russell's
  circular-coplanar reference values.

The new entry `russell-ch4-4.991gG2` catalogues the same cycler
under Russell's framing — both entries are preserved per the
intended M7 collapse-via-canonical-signature semantics.

**Discrepancy noted, NOT silently resolved:** the 0.29 km/s V_inf E
difference and 3-day ToF difference between McConaghy 2006 and
Russell 2004 are larger than rounding alone could explain. The most
plausible reading per Russell's own text (line 7418 "essentially
ballistic for all launch dates ... consistent with the findings in
Ref. 15") is that McConaghy reports ephemeris-optimised values for
a realistic launch while Russell reports circular-coplanar simple-
model reference values. Both characterise the same trajectory.
Captured verbatim in both entries' `notes:` blocks for audit.

**Decision: Purdue dissertation acquisition is no longer required.**
Russell Table 4.9 provides the McConaghy 2006 orbital data with
acceptable precision for M7 matching.

**Audit trail (preserved):**

- 2026-06-01 morning: McConaghy 2005 SnLm broad-class family ingest
  landed (entries 45-47 = `mcconaghy-2005-em-case1`,
  `mcconaghy-2005-em-u0l1`, `mcconaghy-2005-em-snlm-broadclass-family`).
  Partial closure: documented SnLm class membership but did not
  backfill entry 2's null orbital elements (Purdue dissertation
  full text not accessible).
- 2026-06-01 evening: Russell 2004 Chapter 4 tables ingest landed
  (this work). Closed the gap by transcribing aphelion = 1.64 AU
  from Russell Table 4.9 row 1 and adding the cross-reference
  `russell-ch4-4.991gG2` entry. The Purdue dissertation acquisition
  is no longer the gating work; if acquired later it would only
  refine the aphelion value to higher precision and potentially
  add a true (a, e) per-leg breakdown.

**Original question (preserved as audit trail):**

The McConaghy 2006 abstract gives V∞ at Earth (4.7 km/s), V∞ at Mars
(5.0 km/s), and Earth–Mars ToF (153 d), but no orbital elements (a, e,
peri, apo). The full paper is paywalled at AIAA. **Without access to the
paper, we cannot fully specify the canonical signature for M7
matching** — finders that hit this cycler will get `null` matches on
the leg_elements field.

**Recommendation (historical):** the McConaghy 2005 Purdue PhD dissertation
(e-Pubs AAI3166673) is the open-access alternative containing the
broader SnLm taxonomy — queued for future ingest (task #34). When
ingested, the McConaghy 2006 "Notable" cycler should be cross-derived
from its dissertation analog (SnLm sibling family). _Now superseded by
the Russell Table 4.9 ingest._

---

## ◐ Partly resolved (2026-05-31) — C. VISIT-1 / VISIT-2 parameter inconsistency

**Resolution:** commit `b388b8d` — both VISIT entries now carry
arithmetic verification in `period.note` showing that the Rogers 2012
elements (a, e) are internally consistent with the 7-synodic / 14.95 yr
repeat period. The "Wikipedia swap" appears to be a different VISIT
variant (multiple "VISIT-1"-named cyclers exist in the literature).
The Rogers 2012 numbers are taken as authoritative and the YAML now
documents this choice. **Still open:** the original Niehoff 1985 / 1986
sources have not been consulted; if a future ingest of those originals
yields different elements, the entries should be re-evaluated.

**Original question (preserved):**

Wikipedia (citing McConaghy/Longuski/Byrnes 2002 p. 6) and Rogers et al.
2012 Table 1 give contradictory aphelion radii for VISIT-1 and VISIT-2:

| Source | VISIT-1 aphelion | VISIT-2 aphelion |
|---|---|---|
| Wikipedia (citing McConaghy 2002) | 1.89 AU | 1.45 AU |
| Rogers 2012 Table 1 | 1.40 AU | 1.67 AU |

The values appear to be *swapped* — i.e. Wikipedia's "VISIT-1" is
Rogers's "VISIT-2" and vice versa, OR they refer to different
"VISIT"-named cyclers in different papers (Niehoff published several
slightly different variants over 1985-91). Without the original
Niehoff documents (none online), this cannot be resolved.

---

## ✓ Resolved (2026-06-05) — D. Jones 2017 VEM triple cyclers — full member list

The full AAS 17-577 paper was obtained and mined verbatim
(docs/notes/2026-06-05-jones-aas17-577-vem-mining.md). Findings:

- **Member list (family level):** Table 1 (p.3) enumerates all 10
  permitted itinerary families (5 outbound + 5 inbound):
  EMEVE/MEVEM, EMEEVE/MEEVEM, EMEVVE/MEVVEM, EMEVEE/MEVEEM, EMMEVE/MEVEMM.
  This IS the "member list" at the family level.
- **Two fully-specified members ingested** as separate multi-arc rows:
  `jones-2017-vem-emevve-outbound` (Table 2, p.10 — EMEVVE outbound,
  11 encounters over two repeat periods, transit legs 309/259 d) and
  `jones-2017-vem-meevem-inbound` (Table 3, p.10 — MEEVEM inbound,
  transit legs 268/223 d). Both carry full per-encounter V∞ + periapsis
  altitudes from the tables. No orbital elements are published anywhere
  in the paper, so orbit_elements stays null with a data_gap (the single
  biggest gap; do NOT derive — would be circular).
- **6.4-yr (1-synodic) family is EMPTY:** p.8 "No feasible solutions were
  found (of any family) with a repeat period of one synodic period (6.4
  years). ... In contrast, thousands of feasible two-synodic period
  cyclers were obtained." The placeholder `vem-emeeve-3syn` is retained
  but flagged with a `period.feasibility` conflict data_gap recording its
  premise is unrealized in this source.
- **12.8-yr repeat-period correction applied:** p.9 "Recall that the
  repeat period T is 12.8 years." T_syn = 6.4-yr VEM synodic (p.3), so
  "two synodic period" = 12.8 yr, NOT 4.27 yr (2 × E-M synodic). The
  family-seed `jones-2017-vem-triple-family` period was corrected from
  4.27 → 12.8 yr and its `period.basis` data_gap marked resolved.
- **Table 4** (architecture pair) is family-mixed (p.7 "permits the mixing
  of cycler families") and is NOT ingested as a single-cycler row.

See docs/notes/2026-06-05-jones-aas17-577-vem-mining.md for all verbatim
quotes and the honest "not extracted" list (§8).

---

## ✓ Resolved (2026-05-31) — E. spec.md §16.4 attribution correction

**Resolution:** commit `32c5eab`. spec.md §16.4 now cites "Jones,
Hernandez, Jesick (AAS 17-577, 2017)". The Longuski mis-attribution is
documented in `docs/errata-investigation.md` §3.

**Original question (preserved):**

spec.md §16.4 attributed the 2017 triple cyclers paper to "Longuski et
al." Per the NTRS record and the paper's title page, the authors are
Drew R. Jones, Sonia Hernandez, and Mark Jesick (all JPL). Longuski is
not an author.

---

## ✓ Resolved (2026-05-31) — F. spec.md §3 VEM beat period vs. Jones 2017 findings

**Resolution:** commit `32c5eab`. spec.md §3 was updated to clarify that
6.4-yr is the *lowest* natural beat with longer commensurabilities
(12.8 yr, 32 yr) also supporting closure per Jones 2017. The catalogue
preserves both the 2-synodic Jones family-seed entry and the
3-synodic EMEEVE archetype as separate records, accommodating both
readings of the abstract. See `docs/errata-investigation.md` §4 and the
EMEEVE entry's `period.note` for the detailed reconciliation.

**Original question (preserved):**

spec.md §3 says "the natural beat is ≈ 6.4 yr (3 × E–M ≈ 4 × E–V)." But
Jones et al. 2017 found 2-synodic-E-M (4.27 yr) VEM triple cyclers,
which is NOT the 6.4-yr beat. The beat period is sufficient for closure
in the simplified circular-coplanar model with strict commensurability,
but real eccentricities/inclinations + the b-plane DOF open up shorter
periods. M8's enumerator should NOT hard-code the 6.4-yr beat as the
only feasible VEM period.

---

## ✓ Resolved (2026-05-31) — G. Long Mars→Earth return leg of the Aldrin cycler

**Resolution:** commit `b388b8d`. The Aldrin outbound entry's
`legs[1].tof_days` was corrected from 519 d to 634 d with a derivation
source quote: T_cycler (779.8 d) − tof_outbound (146 d) = 633.8 d ≈
634 d. The simplified circular-coplanar model treats the return as a
single aggregate leg without sub-segment breakdown; a future ingest
of the Byrnes/Longuski/Aldrin 1993 paper could split it further if
needed.

**Original question (preserved):**

The 146-day Earth→Mars leg of the Aldrin cycler is well-cited. The
complementary Mars→Earth return — qualitatively described as "16 months
beyond Mars" by Wikipedia — was not cleanly tabulated in any single
primary source accessed during initial compilation. The YAML originally
recorded the return as `tof_days: 519` with an explicit "UNVERIFIED"
note.

---

## H. Out-of-paradigm work flagged 2026-05-31 (NOT in the catalogue)

When the catalogue was extended on 2026-05-31 to carry non-heliocentric
(lunar + Jovian + family-seed Saturnian) cyclers, additional bodies of
work were identified as **adjacent but out of the current cyclerfinder
paradigm**. They are recorded here as awareness for the user / future
implementers, but they are **deliberately NOT added to
`catalogue.yaml`** because:

1. cyclerfinder v1 models cyclers as patched-conic gravity-assist
   sequences with a V∞ + bend-angle abstraction at each flyby.
2. The papers below use fundamentally different mathematical paradigms
   (CR3BP invariant manifolds; low-thrust / solar sail) for which the
   V∞ + bend-angle signature is undefined.
3. Including them as YAML entries would make M7 novelty matching
   meaningless against them: any heliocentric or planet-centric finder
   hit would either falsely match them (signature comparison undefined)
   or never match them (`null` signatures everywhere). Better to flag
   them here and re-evaluate when / if the project adopts those
   modelling paradigms (cf. spec §2 stretch goals).

### H.1 Fantino, Alessi, Peláez Álvarez 2019 — Saturnian CR3BP manifold connections

| Field | Value |
|---|---|
| Title | "Connecting low-energy orbits in the Saturn system" |
| Authors | Elena Fantino, Elisa Maria Alessi, Jesús Peláez Álvarez |
| Venue | 18th Australian International Aerospace Congress (ISSFD-AIAC18), Melbourne, Australia, 24-26 February 2019, paper AIAC18 |
| URL (open) | <https://issfd.org/ISSFD_2019/ISSFD_2019_AIAC18_Fantino-Elena.pdf> |
| Mirror | <https://oa.upm.es/56463/> (Universidad Politécnica de Madrid open repository) |
| Methodology | CR3BP planar Lyapunov orbits + hyperbolic invariant manifolds + low-thrust patches; demonstrates a Tethys→Dione connection of 50 d using 9 kg propellant at 25 mN continuous thrust |
| Why excluded | The patched-conic + V∞ abstraction does not apply to manifold-based low-energy trajectories; the conserved quantity is the Jacobi constant, not V∞ |
| Re-evaluate when | The project adopts CR3BP modelling (would be the natural entry point for CR3BP catalogue ingestion) |

### H.2 Vergaaij & Heiligers 2018 — TU Delft solar-sail Earth-Mars cycler

| Field | Value |
|---|---|
| Title | "Time-optimal solar sail heteroclinic-like connections for an Earth-Mars cycler" |
| Authors | Merel Vergaaij, Jeannette Heiligers |
| Venue | *Acta Astronautica*, 2018, DOI <https://doi.org/10.1016/j.actaastro.2018.06.011> (ScienceDirect S0094576518303734) |
| URL | <https://research.tudelft.nl/en/publications/time-optimal-solar-sail-heteroclinic-like-connections-for-an-eart/> |
| Methodology | Direct pseudospectral optimisation + dynamical-systems heteroclinic connections between Earth-Moon L2 and Sun-Mars L1 libration-point orbits; requires a solar sail to close the connection (no ballistic solution exists). Time-optimal cyclers span ~3 synodic Earth-Mars periods. |
| Why excluded | Low-thrust / solar-sail propulsion is a spec §2 stretch goal, out of v1 scope. The trajectories are not ballistic and have no patched-conic V∞ signature. |
| Re-evaluate when | The project adopts low-thrust modelling. NB the original task brief referred to an "Earth-asteroid" version of this paper; the closest TU Delft paper actually found is this Earth-Mars version. If a distinct Earth-asteroid TU Delft paper exists, it would have the same out-of-paradigm classification. |

### H.3 Pascarella, Woollands, Pellegrini, Sanchez-Net, Vander Hook 2022 — Solar System Pony Express **low-thrust** Earth-Mars cycler (task #38 — **low-thrust component only**)

> **Two-paper split (2026-06-04).** Task #38 was originally filed as "purely low-thrust". This was
> correct for the Pascarella AAS-22-015 paper but incorrectly subsumed the near-ballistic Sanchez Net
> 2022 JSR paper (DOI 10.2514/1.A35091), which IS in v1 scope. See H.4 below for the near-ballistic
> paper. The earlier "purely low-thrust" conclusion applies only to THIS entry (H.3 = Pascarella).

| Field | Value |
|---|---|
| Title | "Low-Thrust Trajectory Optimization for the Solar System Pony Express" |
| Authors | Alex Pascarella, Robyn Woollands, Etienne Pellegrini, Marc Sanchez-Net, Joshua Vander Hook |
| Venue | AAS 22-015 (AAS/AIAA Space Flight Mechanics Meeting); JPL / Univ. Illinois Urbana-Champaign; NASA NIAC Program (80NM0018D0004) |
| Full text | held offline (cite AAS 22-015; DOI `10.1007/978-3-031-51928-4_4`; not stored in repo) |
| Methodology | **Low-thrust** Earth-Mars cycler injection + flyby targeting. NEXT ion engine (Isp 4155 s, thrust 0.235 N), 500 kg ESPA-class courier. Indirect optimal control via Pontryagin's Minimum Principle in Modified Equinoctial Elements, bang-bang throttle with hyperbolic-tangent smoothing + continuation, RK 9(8) propagation, multiple-shooting TPBVP (Matlab `fsolve`). Patched-conic STAR seeds → two-body impulsive → medium-fidelity impulsive → high-fidelity low-thrust → polynomial thrust-arc fitting. |
| Headline result | One fully-converged solution (mid-2037 → early 2046): 500 kg courier inserted into an Earth-Mars cycler using **36 kg** propellant for cycler-orbit injection (COI) + a further **~2 kg** to target **8** subsequent flybys over **~6 years** (>8000 Tbits returned). |
| Underlying cycler families | NOT original to this paper — the patched-conic cycler database is generated by JPL STAR (Landau/Campagnola/Pellegrini 2022) and the family lineage traces to Byrnes/Longuski/Aldrin 1993 (ref 3) and Russell & Ocampo 2006 (ref 13), both already represented (within-paradigm) in the catalogue. |
| Why excluded | **Low-thrust regime (deferred scope #37).** The catalogue is ballistic-only (v1 patched-conic + V∞/bend-angle abstraction). The paper presents NO standalone sourced ballistic anchor for a *named* cycler: there is no tabulated V∞ at Earth/Mars, no period/`k`, and no per-cycler orbital elements (`a`, `e`) — the only published numerics are propellant masses, flyby altitude constraints (≤25,000 km; ≥300 km Mars / ≥1000 km Earth), and the SEP hardware spec. A null-V∞/null-orbit row would carry no signature for M7 novelty matching (see exclusion rationale 1–3 above) and would not be golden-citable. |
| Re-evaluate when | The project adopts low-thrust modelling (spec §2 stretch goal / #37). At that point this is a natural ingestion candidate — but the ballistic seed elements would still have to come from the STAR / Russell-Ocampo / Byrnes sources, not from this paper. |

### H.4 Sanchez Net, Pellegrini, Parker, Vander Hook, Woollands 2022 — Solar System Pony Express **near-ballistic** Earth-Mars cyclers (task #38 — **in v1 scope**)

| Field | Value |
|---|---|
| Title | "Cycler Orbits and Solar System Pony Express" |
| Authors | Marc Sanchez Net, Etienne Pellegrini, Jeffrey Parker, Joshua Vander Hook, Robyn Woollands |
| Venue | *Journal of Spacecraft and Rockets*, Vol. 59, No. 3, pp. 861–870, 2022 |
| DOI | [10.2514/1.A35091](https://doi.org/10.2514/1.A35091) |
| Full text | held offline (not stored in repo) |
| Scope | **NEAR-BALLISTIC, IN v1 SCOPE.** ΔV ≤ 10 m/s patched-conic Earth-Mars cyclers. Covers EM and EEM families over 2030–2034 launch windows. |
| Status | **Both Cycler 1 (EEM, `sanchez-net-2022-eem-cycler1`) and Cycler 2 (EM, `sanchez-net-2022-em-cycler2`) are now catalogue entries** (added 2026-06-04, `model_assumption: analytic-ephemeris`, `delta_v_kms: 0.005/0.007`). Data gaps remain: no tabulated `a`/`e`/`i`, partial V∞ for intermediate legs. See `docs/notes/s1l1-target-topology-mining.md` §2.9 for the mined event tables. |
| Note | Distinct from the low-thrust Pascarella AAS-22-015 (H.3 above). Both papers share authors and the SSPE mission context but use fundamentally different trajectory regimes. |
| Re-evaluate | Per-cycler orbital elements (`a`, `e`, `i`) remain data gaps. The numbered network cyclers in Tables 3–5 (IDs 51, 84, etc.) cannot be added without the Star database output. |

### H.5 Patel 2019 — Earth-Mars cycler vehicle conceptual design (in-scope reference for `s1l1`)

| Field | Value |
|---|---|
| Title | "Earth–Mars Cycler Vehicle Conceptual Design" |
| Authors | Patel (first name withheld pending verification) |
| Venue | M.S. thesis, Florida Institute of Technology, December 2019 |
| Full text | held offline (not stored in repo) |
| Scope | **IN v1 SCOPE — vehicle design reference for the `s1l1` catalogue entry.** Takes S1L1 trajectory from McConaghy 2006 circular-coplanar baseline. Key sourced values: Earth flyby V∞ 3.657 km/s, no Mars V∞ tabulated; `a`=1.30 AU, `e`=0.257; E-M-E-E topology; Earth encounters at t=0, t=2.8276 yr, t=4.286 yr. |
| Note | Patel derives vehicle mass estimates and habitat sizing from the S1L1 trajectory. The trajectory parameters (a, e, V∞) match McConaghy 2006 and Rogers 2012 Table 1 within rounding. Useful as a cross-check for the `s1l1` catalogue entry and as a design reference if the project adds vehicle-sizing fields. |

These flags are deliberately separate from the within-paradigm questions
A–G because those are gap-filling (missing numerics, inconsistent
secondary sources, attribution corrections) while H.1–H.3 are
paradigm mismatches. H.4 and H.5 are in-scope near-ballistic or
design-reference papers newly acquired.

## Resonance-anchored construction validated (2026-06-03)

The McConaghy/Russell construction method now reproduces cyclers from their
SOURCED orbital elements with no optimisation (`search/resonant_construct.py`,
`core/kepler.coe_to_rv`): S1L1 -> V_inf 4.90/4.98 km/s (matches Russell 2004
coplanar 4.99/5.10) and Aldrin -> 6.58/9.75 (matches sourced 6.5/9.7). The
spec S1L1 5.65/3.05 is a higher-fidelity figure (Mars 3.05 needs eccentric
Mars). Catalogue-scale reproduction (`scripts/batch_resonant_reproduction.py`)
is currently data-limited to 3 rows (only those carry both cycler-level (a,e)
and sourced V_inf); unlocking it needs per-cycler (a,e) populated for the
ballistic Russell rows from the Russell 2004 dissertation (in docs/refs/).

## S1L1 V∞ anchor flagged unverified-provenance (2026-06-04)

The `s1l1-2syn-em-cpom` entry's 5.65/3.05 km/s V∞ pair has been flagged
`unverified-provenance` in the catalogue's `data_gaps[]` block. Mining of
Patel 2019, McConaghy 2006 (abstract + Patel secondary), and Sanchez Net 2022
(DOI 10.2514/1.A35091) did not find this pair in any primary source:

- **Patel 2019** (= McConaghy 2006 circular-coplanar) gives Earth flyby V∞
  3.657 km/s; no Mars V∞ is tabulated.
- **Sanchez Net 2022** near-ballistic real-date regime gives Mars V∞ 5.2–7.3
  km/s and Earth V∞ 3.6–7.3 km/s; none match 5.65/3.05.
- **McConaghy 2006 JSR abstract** (real-ephemeris) gives ≈4.7 km/s Earth /
  5.0 km/s Mars.
- **Russell 2004 coplanar** gives 4.99/5.10 km/s (separate catalogue entry
  `russell-ch4-4.991gG2`).

The 5.65/3.05 values trace only to spec.md §9; their origin is undetermined
until the McConaghy 2006 JSR full text is accessed. They are retained in place
(not overwritten) pending that access. See `docs/notes/s1l1-target-topology-mining.md`
for the full sourcing analysis.

## Newly mined references — Agrawal 2022, Landau-Longuski 2009, Howe 2025 (2026-06-04)

Three papers mined for S1L1 real-ephemeris targets and near-ballistic EEM/EM
catalogue candidates. Full extraction in `docs/notes/2026-06-04-agrawal-landau-howe-mining.md`.

### Agrawal 2022 — not catalogue-eligible

**Agrawal, R. 2022.** "Design and Analysis of an Orbital Logistics Architecture
for Sustainable Human Exploration of Mars." PhD dissertation, Purdue University.
Advisors: James M. Longuski, Sarag J. Saikia.

Finding: logistics/Spacedock focus; cyclers mentioned only in future-work (pp. 72,
129–130). No cycler trajectory numbers (no orbital elements, V∞, or per-arc ToFs).
Not catalogue-eligible.

### Landau & Longuski 2009 — not catalogue-eligible

**Landau, D. F. and Longuski, J. M. 2009.** "Comparative Assessment of Human–Mars-Mission
Technologies and Architectures." *Advances in the Astronautical Sciences*, Vol. 126.
School of Aeronautics and Astronautics, Purdue University. [ScienceDirect pii
S0094576509001118; scanned document, image-only PDF.]

Finding: IMLEO architecture comparison across six mission classes (Direct, Semi-Direct,
Stop-Over, M-E Semi-Cycler, E-M Semi-Cycler, Cycler). E-M semi-cycler ranks #1 under
advanced propulsion scenario (Table 13). Cites McConaghy 2006 (ref [45]) as the canonical
S1L1 source but reproduces no per-cycler V∞ values or orbital elements — results are
IMLEO graphs only, not tabulated trajectory data. Not catalogue-eligible.

NB: Authors are Landau & Longuski (not Sarunic); image-only PDF confirmed.

### Howe et al. 2025 — candidate only, not added

**Howe, A. S., Blincow, J., Hall, T. W., and Leonard, C. 2025.** "Tackling a Mars Cycler
Design Head-on." ICES-2025-555, 54th International Conference on Environmental Systems,
13–17 July 2025, Prague.

Finding: habitat/structural design for an "escalator" cycler (citing Rauwolf, Friedlander,
Nock 2002). Up-escalator E→M ToF 151 d (Day 1→151); down-escalator stated 170 d (but
arithmetic from printed day numbers 38→229 gives 191 d — internal inconsistency in the
paper). V∞ < 6 km/s stated as a design constraint, not a computed encounter value.
Two escalator cyclers cited (Rauwolf 2002) lack tabulated encounter V∞ — candidate-only,
not added to catalogue. Primary source for escalator elements is Rauwolf 2002 (not in
scope here).

---

# Forge Phases 4 + 5 shipped — first novelty candidates queued (2026-06-06)

> New section, appended; does not modify any entry above. Plan completion notes:
> `docs/superpowers/plans/2026-06-03-the-forge-pipeline.md` (Completion notes —
> Phases 4 + 5). Modules: `src/cyclerfinder/data/discover_novel.py`,
> `src/cyclerfinder/verify/adversarial.py`, `src/cyclerfinder/data/review_queue.py`,
> `scripts/forge_novelty_run.py`.

## What shipped

- **`discover_novel`** — construction-first novelty loop over the **E-M
  multi-arc space** (E-M-E-E topologies x parallel epoch scan grid,
  `scan_parallel`, 16-core). Each closed chain is bridged to a full `Cycler`
  (`cycler_from_closure`, via `construct_cycler` at the converged epochs),
  signed (`canonical_signature`), matched against the catalogue
  (supersession-aware, R1 delta 3), code-path cross-checked (Axis A,
  `crosscheck_code_paths`), and routed through `run_gauntlet`.
- **Adversarial panel** (`adversarial_panel`) — N independent verifiers per
  candidate: a falsification probe (claimed self-consistency must hold), an
  independent re-closure from the reported seed, and perturbed-seed robustness
  re-solves. Majority-refute kills. Proven teeth: a fabricated impossible-bend
  candidate is majority-refuted AND auto-REJECTED at the verdict.
- **Human-review queue** (`data/review_queue.jsonl`) — SILVER candidates with
  full audit trail (verdict axes + panel result). Explicitly NON-catalogue
  (`is_catalogue_source() is False`); the validator refuses to queue a
  non-SILVER/GOLD or panel-refuted entry. **No catalogue row is ever created
  by the loop** (golden discipline) — the census ratchets are unaffected.

## First real run result (the science)

16-core, tight 2030 Sanchez launch window, 16 epochs x 4 E-M-E-E topologies =>
**12 distinct closed families**:

- **10 REJECTED** — bend-infeasible (the required flyby turn exceeds the
  V_inf-limited maximum). A bend-infeasible closure is physically inadmissible
  and is auto-falsified at the verdict (never SILVER).
- **2 SILVER** — bend-feasible, machine-confirmed (Axis A: in-house Lambert vs
  lamberthub + forward-Kepler re-prop agree; the single-ellipse resonance
  construction path is correctly demoted to *unavailable* for a real-DE440
  multi-arc chain, since comparing it would be the cross-fidelity confusion the
  Forge refuses), panel-survived, queued. Peak V_inf ~13.0 and ~12.1 km/s; both
  `match=novel`.

**The prototype's closed family does NOT match any catalogue row.** Its V_inf
regime (~9.75 Earth / 13.0 Mars km/s) sits far above the sourced Sanchez/S1L1
anchors (~3-6 km/s) — the documented "S1L1 Mars-V_inf ~6.4 floor generalised"
(project memory `project_s1l1_realeph_closure_blocker`). So this is **not a
rediscovery**: these are the project's first machine-confirmed novelty
candidates, held at SILVER pending human review. Most likely human resolution:
they are a *different / higher-energy* (possibly powered) family than the sourced
near-ballistic Sanchez cyclers — which is exactly why SILVER caps them pending a
human, rather than asserting a novel discovery.

## Why VEM was not scanned (recorded reason)

The #110 dense VEM scan (2816 points/row, 16-core) produced ZERO bend-feasible
closures (floors 17.9/18.5 km/s vs sourced Jones 2.42-7.0). VEM
single-ellipse-per-leg ballistic novelty is empirically nil; `discover_novel`
defaults to the demonstrably-bend-feasible E-M multi-arc space and does not burn
the loop's budget on VEM. Reaching the Jones VEM family needs 3D inclination,
real-eccentricity intermediate flybys, or multi-arc-per-leg seeding (M-3D /
future), not more scan density.

## Open follow-ups (not blockers)

- The 2 SILVER candidates await **human review** in the queue. A human decides:
  promote to a sourced catalogue row (requires an independent source — none
  exists yet, so promotion to GOLD is not currently possible), retain as a
  documented novel-machine-confirmed family, or reject as a powered/off-family
  artifact. The loop will not act without that decision.
- The unbounded loop-until-dry deepening (broader epoch spans, more topologies,
  E-E-M and longer multi-arc classes) is a future scaling pass; this run is the
  bounded first sweep the plan specified.

---

# Landed 2026-06-07 — V2 promotion, MBH, primer, Hughes cross-check, method-mining

> New section, appended; does not modify any entry above. Records the day's
> landed results against the cited commits / notes. Running work is marked
> in-flight (no outcomes claimed).

## §14 V2 class-split — first V2 catalogue row

The single §14 V2 multi-lap-periodicity gate was split into **V2-ballistic**
(the existing rotating-frame-repeat drift over ≥3 laps, in the row's defining
model) and **V2-powered** (≥3 consecutive cycles, every encounter achieved
with the per-cycle maintenance applied AND bounded intra-cycle drift vs the
planned trajectory, reset at each maneuver) — commit `ee9f854` (spec §14 +
§16.7.12). Rationale (task #134): the powered Aldrin cycler is retargeted
every cycle by design, so the cross-cycle rotating-frame-repeat metric is
structurally unsatisfiable for it (~4.14e8 km / 3 laps, ~2072× tol) — the gate
measured the wrong thing.

- **Aldrin OUTBOUND clears V2-powered** (`verify/v2_powered.py`,
  `tests/verify/test_aldrin_v2_powered.py`, slow; commit `69f2593`): 3
  consecutive in-family cycles, Mars-flyby V∞ continuity ≤1e-6 km/s,
  intra-cycle Kepler forward-reprop ≤0.002 km, in-family maintenance ΔV
  2.76–2.91 km/s/cycle. Row `aldrin-classic-em-k1-outbound` promoted **V1→V2**
  (commit `6263548`).
- **Aldrin INBOUND stays V1** — its real-window solve lands on a ballistic
  ΔV≈0 off-family neighbour, so "maintenance applied" is not demonstrated (the
  recorded #134 off-family resolver issue).
- **The four #137 free-return rows honestly do NOT pass V2-ballistic** (commit
  `78bc7f0`, `tests/search/test_free_return_v2_ballistic.py`). They are
  `cycler_class: multi-arc`: the #137 V1 evidence closes a single E→M→E
  free-return ellipse SLICE (~620–708 d, only ~0.3–0.4 of the catalogue
  4.27/6.41 yr period); no continuous ≥3-lap trajectory exists to propagate
  (measured 3-lap drift ~9.4e7–1.2e8 km, >10× the M6a 50,000 km tol). The four
  rows stay V1 — a structural finding, not a near-miss.

**Validation census (verified `tests/data/test_schema_v45_fields.py::
test_live_v1_census_matches_recorded_evidence`):** exactly one V2 (the Aldrin
outbound) and five V1 (Aldrin inbound + `russell-ch4-5.30gGf3`, `-9.94Gg3`,
`-5.75ggF3`, `-9.353Gg2`); no V3+. Schema bumped to v4.5.

## #106 SnLm Phase 2 re-scoped complete — DESCRIPTOR_CLOSABLE bucket

The `multi_encounter_sequence` exclusion bucket was split **204 → 12
DESCRIPTOR_CLOSABLE + 192 MULTI_ENCOUNTER_SEQUENCE** (commit `9fe26e3`,
`tests/test_catalogue_rediscovery.py` / `tests/_catalogue_loader.py`). The 12
`free_return_arcs[]`-bearing SnLm rows (the McConaghy/Russell-12 set) became
discover-reachable through the #137 free-return descriptor path: `discover.py`
gained `discover_free_return()` and the descriptor-closable code path (commit
`cc65c39`, `tests/data/test_discover_free_return.py`), making the Russell-12
set discover-reachable rather than blanket-excluded. Plan re-scope recorded in
commit `d60886d`. NB this is a `tests/test_catalogue_rediscovery.py` bucket
re-partition (additive split of the same 237 rows), not a census change.

## #142 method-mining of six held PDFs (commit `fee7c7b`)

Six method-mining notes dated 2026-06-07 (cite by author/title/venue only).
Headline corroborations:

- **B-plane targeting + powered-flyby kernel** is the implementable Phase C
  shooter blueprint: Jones / Hernandez / Jesick (AAS 17-577, 2017) Eqs. 1–5
  (transfer/turn angle, iterative periapsis radius, periapsis speeds, B-plane
  frame Ŝ/T̂/R̂, B-plane angle θ_B) — **independently corroborated** by Russell
  2004 dissertation Eq. 5.5 (powered-SOI ΔV).
- **S1L1 multi-arc blocker independently confirmed** by Hollister & Rall 1970
  (NASA CR periodic-orbits study): a periodic E–M orbit is intrinsically
  multi-arc (reciprocal E–M–E round trips stitched by Earth-side returns;
  single Mars-return arcs are infeasible — the trajectory intersects Mars).
  1969-vintage corroboration that S1L1 never closing as a single ellipse is a
  modelling-truth, not an infra bug.
- **#137 genome pivot validated** by Russell 2004 §5.4.2: rejecting
  integer-structured analytic transfers in favour of a shape that can morph
  under a gradient optimiser is the canonical written justification for the
  radial-crossing `(a,e)` free-return genome; our genome is the symmetric
  special case of his ψ generic-return.
- **Bend formula corroborated across four papers** — `core/flyby.max_bend`
  equals Russell Eq. 5.4, Jones Eq. 2, Ozimek (AAS 19-348) Eq. 17, and the
  Vasile & Campagnola DFET formula.

## #143 external novel-algorithms survey (commit `e1885dd`)

Note `docs/notes/2026-06-07-external-algorithms-survey.md` — a ranked roadmap
across four external algorithm families: primer-vector diagnostic (cheapest,
read-only) → Monotonic Basin Hopping (basin-selection) → STOUR cross-check →
family continuation.

- **Russell & Ocampo (2006), "Optimization of a Broad Class of Ephemeris Model
  Earth–Mars Cyclers," JGCD 29(2):354–367** flagged as the **near-exact
  frontier match** (circular-coplanar→ephemeris continuation; the V1→V3 bridge
  for the Russell rows) — **top of the acquisition list**.
- **EMTG** (NASA Goddard) is blocked on the commercial SNOPT dependency;
  **pykep / `pygmo.mbh`** (ESA, MPL-2.0, no SNOPT) is the open adoption path
  for wrapping our existing `ballistic_correct`.
- Caveat recorded in the note: four primary PDFs returned binary/encoded
  content on automated fetch, so their internal numeric tables were NOT read —
  citations are confirmed via metadata, but any number must be verified on
  acquisition before being quoted as sourced.

## #147 Hughes-Edelman-Longuski (AIAA 2014-4109) mined (commit `65aecee`)

Note `docs/notes/2026-06-07-hughes-2014-fast-mars-free-returns-mining.md` —
the STOUR cross-check for the Jones VEM frontier. Headline: **pure-ballistic
E-V-M free returns exist only at node-crossing phasings, in the 5–7 km/s
Mars-side V∞ class**; reaching *lower* V∞ is explicitly gated behind a
broken-plane (DSM) maneuver — i.e. a multi-arc / powered topology, NOT
single-ellipse ballistic. This **independently corroborates** the #110 /
M-ED single-ellipse blocker and the S1L1 multi-arc finding: our #110 floors
(17.9 / 18.5 km/s) were wrong-basin artifacts, and the Jones ~3 km/s class is
not reachable as a pure free return here.

- Six golden-eligible cross-check tables (Tables 1–6, verbatim with page
  numbers); Table 4 is a sourced STOUR-vs-STK-Astrogator high-fidelity
  agreement (~0.06 km/s V∞, ~1 d, ~20 km altitude) — a validation anchor for
  the patched-conic→n-body verify chain.
- **v4.2 backfill / catalogue-eligibility (confirmed NO):** these are one-shot
  human-flyby free returns, NOT cyclers → **NOT catalogue rows**. They ARE
  valid as cross-check / golden-test anchors for a patched-conic E-V-M
  *pipeline*, not as catalogue goldens.
- Upstream acquisitions flagged: Okutsu & Longuski 2002 ("Mars Free Returns
  via Gravity Assist…", ref 11) and Patel, Longuski & Sims 2002 ("Mars Free
  Return…", ref 10) — the acquisition-worthy STOUR free-return sources.

## #145 MBH wrapper landed (commit `54d5075`)

Transcription-agnostic Monotonic Basin Hopping loop (`search/mbh.py`;
deterministic audit trail, rng_seed required) with adapters wrapping the
existing `ballistic_correct` and the #137 `free_return_correct` correctors
without editing them. Note `docs/notes/2026-06-07-mbh-wrapper.md`.

- **Basin selection is curable:** on `mcconaghy-2006-em-k2`, from a 40-day
  off-phase mis-seed, the plain solve fails (lands at 2.16 km/s) but MBH
  recovers the sourced-anchor (Rogers ellipse) basin — emerged V∞ M=4.69 /
  E=4.05 vs sourced 5.0 / 4.7.
- **Honest negative:** the `russell-ch4-6.44Gg3` probe (catalogued multi-arc;
  sourced V∞ E=6.44 / M=3.74) does NOT recover — MBH confirms the multi-arc
  topology blocker internally rather than hopping into a non-existent
  single-ellipse basin.
- The perturbation distribution default is documented, NOT sourced — awaits
  Englander & Englander 2014.

## #144 primer-vector diagnostic landed (commit `b0dfe7d`)

Read-only Lawden / Lion & Handelsman 1968 first-order optimality check
(`verify/primer.py`; p̈ = G(r)p, propagated via the 6×6 STM along the Kepler
coast; per-coast max|p| → OPTIMAL / IMPROVABLE). Note
`docs/notes/2026-06-07-primer-vector-diagnostic.md`.

- **Aldrin maintenance schedule verdict: IMPROVABLE_ADD_IMPULSE.** Coast 0
  (E→M, 131.9 d) grid-converged max|p| = 1.1223 at t/T ≈ 0.355 (≈47 d) — a
  genuine interior bulge; coast 1 (M→E) marginal at max|p| = 1.00008 (noise
  floor).
  **RE-LABELED 2026-07-15 (Fable adversarial pass, full-file audit — Guzman 2002 was already
  acquired+mined, this PROVISIONAL gate was stale):**
  - **Coast 0:** `CONFIRMED-with-scope per Guzman 2002 §"Isolated Singularities" (p.7-8): a
    sub-one-rev (Δν<180°) singularity-free arc; the |p|=1.122 necessary-condition flag is valid
    and the Φ_rv inversion is well-posed. Remains necessary-not-sufficient (Guzman p.3).`
  - **Coast 1:** `CONFIRMED-with-scope: sub-one-rev, no Φ_rv singularity (Guzman p.7-8); marginal
    endpoint touch, treat as necessary-conditions-met` (598.7 d is ~0.87 of Mars's period —
    sub-one-rev, not the "multi-revolution-scale" this entry originally worried about).
  Machinery validated for both arcs; only necessary-not-sufficient, never sufficient, per
  Guzman's own p.3 disclaimer. See `docs/notes/2026-06-07-guzman-2002-primer-survey-mining.md`.
- **Methods correction recorded:** the brief's Hohmann-ratio-20 interior-bulge
  gate was physically incorrect — the ~11.94 threshold is Lawden's
  endpoint coast-extension (bi-elliptic) condition, not an interior |p|>1 on
  the two-impulse arc. The golden gate was corrected to a long-way
  (>180°) Lambert transfer → IMPROVABLE; symmetric Hohmann transfers (ratios
  2, 11.94, 20, 50) touch unity only at the endpoints → OPTIMAL.

## Acquisitions ledger (2026-06-07)

**New holdings (held offline; cite by author/title/venue — never a repo
path):**

- **Hughes, Edelman & Longuski, "Fast Mars Free-Returns via Venus Gravity
  Assist," AIAA 2014-4109** — **mined** this run (see #147 above).
- **Zhang 2026** — neural-network porkchop for low-thrust asteroid transfers.
  Held, **unmined**; likely out of scope (low-thrust). Scope decision pending
  user.
- **Zhang 2024** — neural angle-only orbit determination at Earth–Moon
  libration points. Held, **unmined**; likely out of scope (OD, not cycler
  trajectory design). Scope decision pending user.

**New wants (from #143 / #147; prioritised):**

1. **Russell, R. P. & Ocampo, C. A. (2006), "Optimization of a Broad Class of
   Ephemeris Model Earth–Mars Cyclers," JGCD 29(2):354–367** — top of the
   list; the circular-coplanar→ephemeris continuation V1→V3 bridge.
2. **Englander, J. A. & Englander, A. C. (2014), "Tuning Monotonic Basin
   Hopping…," 24th ISSFD, paper S7-3** — the Cauchy/Pareto perturbation-spec
   reference before implementing MBH (open PDF located).
3. **Guzman, Mailhe, Schiff, Hughes & Folta (2002), "Primer Vector
   Optimization: Survey of Theory, New Analysis and Applications,"
   IAC-02-A.6.09 (53rd IAC, Houston; NTRS 20030032208)** — ~~needed to lift the
   #144 primer diagnostic from PROVISIONAL~~ **STALE (corrected 2026-07-15, Fable adversarial
   pass, full-file audit): already acquired and mined.** See `docs/notes/CORPUS_INDEX.md` (row:
   "mined") and `docs/notes/2026-06-07-guzman-2002-primer-survey-mining.md`, which already
   delivers the re-label verdict: "PROVISIONAL pending Guzman 2002 (multi-rev caveats)" →
   "CONFIRMED-with-scope per Guzman 2002 §Isolated Singularities (p.7-8): necessary-not-sufficient
   (Guzman p.3)." This "New wants" item is done; #144's own entry should be updated to the
   confirmed-with-scope label rather than left saying PROVISIONAL.
4. **Lion, P. M. & Handelsman, M. (1968), "Primer Vector on Fixed-Time
   Impulsive Trajectories," AIAA Journal 6(1):127–132, DOI 10.2514/3.4452** —
   the add-an-impulse diagnostic the maintenance work would implement.
5. **Okutsu & Longuski (2002), "Mars Free Returns via Gravity Assist from
   Venus"** — STOUR free-return upstream (Hughes ref 11).
6. **Patel, Longuski & Sims (2002), "Mars Free Return…"** — STOUR free-return
   upstream (Hughes ref 10).

## Open / running (in-flight — no outcomes claimed)

- **#133 Phase C n-body shooter** — running concurrently; B-plane targeting
  kernel mandate (Jones Eqs. 4–5, the #142 finding).
- **#148 add-an-impulse recoverable-ΔV** — running concurrently.
- **#146 viz 2c sampled-trajectory variant** — queued.
- **#76 moon-tour** — ✓ SHIPPED (2026-06-08): Tier-1 patched-conic moon systems
  + VILM landed; see the "Moon-tour Tier-1 — SHIPPED" entry in the Done section
  above. Tier-2 (CR3BP) remains open.
- **#128-S2** — queued.
- **#172 Forge Phase 6 — first novelty campaign** — ✓ RAN (2026-06-09):
  Jovian Galilean I-E-G VILM sweep, 64 epochs × 2 topologies (128 pts). Verdict
  **EMPTY** (the expected base rate): 12 closed, all `novel`, 0 bend-feasible
  (V∞ 8.3–26.8 km/s vs ~6 km/s floor; gap 20.8 km/s), 0 SILVER. The #76
  honest-risk generalises — a rigorous bounded method-versioned negative recorded
  in `data/empty_regions.jsonl` (region `jovian-IEG-vilm-2026-06-09`). Dedup
  firewall + capability-subsumption re-sweep gate both verified live. Pipeline
  (VILM prune + empty-region registry + literature-check field + re-sweep gate)
  SHIPPED. See `docs/notes/2026-06-08-forge-phase6-jovian-sweep-results.md`. A
  multi-arc / n-body / low-thrust method will auto-re-sweep this region per §6b.
