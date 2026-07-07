# Outstanding questions — cyclerfinder catalogue

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

# Project state at a glance (updated 2026-06-30)

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
- **Tier 3 (validation + consolidation):** PC (3,2) V3 upgrade gated on SAT441l kernel (#506 scope done);
  #487 V4_qp gauntlet (de-prioritised).
- **Tier 5 (3D Dynamics & Multi-Rev Search):** #515 3D lift framework for cross-system cycles; #516 Multi-Revolution 3D Patched Search (n_em, n_se > 1); #517 Asymmetric/Mixed Libration Pairs in 3D; #518 3D BCR4BP Continuation; #520 Comprehensive 3D sweep (needs a positive control before its negative is trustworthy — see DELTA above).
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

**TASK ALLOCATIONS (next-unused per [[project_task_numbering_convention]]; #512-#514 committed; #515-#518 for session C working-tree; #519 for low-thrust proposal; #520 for the comprehensive sweep; #521-#526 for the 2026-07-02 review's gate + novel-orbit proposals; #527-#529 for the same-day second-pass review; #530 for the #523/#527-motivated unstable-manifold follow-up; #531 for the #314-reuse heteroclinic-connection follow-up; #532 for the multi-orbit resonance-hopping follow-up; #533 for the genuine QBCP model build; #534 for the #522-split single-system torus connection search; #535 for the transient-drift-phase quasi_cycler search; #536 for the Fable-review-motivated Jovian-moon-tori heteroclinic screening follow-up; #537 next-unused):**
- **#512** — (n_em, n_se) Resonance Sweep: Run sweep driver and build analytic wrap table for #411 cross-system cycle. (Resolved)
- **#513** — R52-U Recovery: Recover R52-U from sourced Braik-Ross initial conditions to partially flip the C32-dominance gate. (Resolved)
- **#514** — NAIF Kernel-Freshness Checker: Build monthly workflow and document NAIF kernel freshness. (Resolved)
- **#515** — 3D Lift Framework: Shipped cross-system cycle corrector over fixed out-of-plane amplitudes, physically scaling system lengths and aligning Floquet signs (uncommitted working-tree).
- **#516** — Multi-Revolution 3D Patched Search: Search with n_em, n_se > 1 to bypass the single-revolution phase-closure wall (uncommitted working-tree).
- **#517** — Asymmetric/Mixed Libration Pairs in 3D: Scan for EM-L1 <-> SE-L2 and EM-L2 <-> SE-L1 crossings (uncommitted working-tree).
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
  mechanically skipped. `uv run pytest tests/data tests/search -q` stays green. Building the actual
  `preflight_search()` gate + AST ratchet (Phase 2) remains open.
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
- **#536** — Apply the Linking-Number/QP-Torus Heteroclinic Screening Tool to a Genuinely Unmapped System
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
- **#533** — Genuine Coherent QBCP Model (allocated 2026-07-03, same session — path (a) from #522's scoping pass,
  formally split out as its own task at the user's request rather than left as inline prose under #522). Closes
  the real prerequisite gap #522 found: `core/bcr4bp.py`'s own docstring (lines 14-30) states outright that
  Andreu's coherent QBCP needs 8 Fourier tables `alpha_i(theta_S)` that are "NOT in the in-repo digest" — a
  literature-data-acquisition gap, not a missing function. Concrete first step: acquire + digest the Andreu
  QBCP paper(s) (the alpha_i Fourier-table source; not yet identified/acquired this session — check
  Gomez-Masdemont-Simo-era QBCP literature and any existing CORPUS_INDEX near-misses first) and extract the 8
  tables. Then build the QBCP equations-of-motion module (self-consistent Sun trajectory, not the circular
  approximation `core/bcr4bp.py` already has) with a positive control against a published QBCP invariant-object
  reproduction (e.g. a QBCP halo/Lissajous family or torus already in the literature). Real multi-week
  acquisition+build, matching the cost of the other Track-A capability axes. Explicitly an ALTERNATIVE to
  #522's cheaper, unvalidated path (b) (anchoring `theta` to the Sun's synodic angle already tracked in the
  existing EM-frame `core/bcr4bp.py`) — the two are not both required; #522's cross-system SE<->EM torus
  connection can proceed via EITHER #533 (if built) OR path (b) (if it validates), whichever the user/future
  session picks. Not yet started — no acquisition, no digest, no code.


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
  flip the C32-dominance gate. [active follow-on]
- **#498** — Acquire+mine MOON-TOUR / GRAVITY-ASSIST Ross papers (all MISSING): "Design of a multi-moon
  orbiter" (Ross-Koon-Lo-Marsden 2003, AAS 03-143); "Multiple gravity assists, capture, and escape in the
  RTBP" (Ross-Scheeres 2007, SIADS 6(3), +control-map software); "Constructing a low energy transfer between
  Jovian moons" (KLMR 2002, Contemp. Math. 292); "Controlled Keplerian map" (Grover-Ross 2009, JGCD 32(2));
  "Resonance and capture of Jupiter comets" (KLMR 2001, CMDA 81); "Geometric mechanics and the dynamics of
  asteroid pairs" (Koon-Marsden-Ross-Lo-Scheeres 2004 — binary dynamics, also #494/#308). Feeds
  #318/#465/#494/#500.
- **#499** — Acquire+mine HETEROCLINIC FOUNDATIONS (MISSING): "Heteroclinic connections between periodic
  orbits and resonance transitions" (KLMR 2000, Chaos 10(2) 427-469 — ref [25] of Ross-RT); "Heteroclinic
  Transfer Between L1 and L3 in EM" (Braik-Ross 2025, AAS 25-716); "Connecting orbits and invariant
  manifolds in the spatial RTBP" (Gomez et al. 2004, Nonlinearity 17 — ref [26]); "The Genesis trajectory
  and heteroclinic connections" (KLMR 1999, AAS 99-451); "Transport of Mars-crossing asteroids from the
  quasi-Hilda region" (Dellnitz et al. 2005, PRL 94 — set-oriented GAIO, also #308); "Experimental
  validation of phase space conduits" (Ross-BozorgMagham-Naik-Virgin 2018, PRE 98 052214 — the ΔC
  tube-cross-section scaling = ref [28], underpins the #494 construction). Feeds #314/#405/#411/#496.
- **#500** — SPECULATIVE genome: evaluate the controlled Keplerian map (Ross-Scheeres 2007 / Grover-Ross
  2009) as a moon-tour gravity-assist cycler genome — complements #465 multi-rev leveraging + #318 joint
  search. Gated on #498.
- **#501** — #318 full-scale, REFRAMED: broaden sequences/systems (NOT densify CGCEC — smoke showed a mined
  regime, 0/256) joint-search campaign with post-hoc lit-novelty. [held — the deferred #318 decision]
- **#502** — Watch/acquire the LONGER companion of Ross-RT 2026 (full μ-grid + all 9 EM family tables) for a
  tighter #494 Pluto-Charon continuation seed. [watch]
- **#503** — Acquire+mine cislunar-resonance / transport companions: Rawat-Kumar-Rosengren-Ross 2026 JGCD
  49(4) "Cislunar Mean-Motion Resonances"; Onozaki-Yoshimura-Ross 2017 ASR 60 (4-body tube dynamics);
  Fitzgerald-Ross 2022 ASR 70 (periodically-perturbed RTBP, → #292/#293); Naik-Lekien-Ross 2017 RCD 22(3)
  "Phase Space Transport / Lobe Dynamics" (+Lober software). → expand #267 resonance_network goldens.
  [lower priority]

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
- **#318 OPEN — Phase-2 REFRAMED 2026-06-30** Multi-axis joint search (powered × multi-rev ×
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
- **#405 ✓ Phase A DONE (2026-06-20)** Cross-system SE↔EM heteroclinic-cycle
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
- **#315** Circumbinary/binary-star μ-gap sweep
- **#316** Cross-system cycler framework (Sun-Earth ↔ Earth-Moon manifolds)
- **#320** First quasi_cycler discovery sweep (blocked by #319)

**Infrastructure + polish:**
- **#307** #289 Phase 5 (DSM + multi-rev + eccentric Tisserand) — gap from #300/#302
- **#310** Single-orbit prioritizer adapter (#284 architectural gap)
- **#317** PINN-based pre-filter for sweep-impossible regions
- **#321** Multi-threaded inner-loop compute (joblib wrappers — 4-8× sweeps on multi-core)
- **#322** Tulip petal_count z0-collapse bug fix (in flight) — surfaced by #313

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
- **#245/#246** — FBS optimizer default flip + documentation; HELD behind #248
  (no convergence parity to flip toward yet).
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
  maintenance budget. **V3 writeback recommended, HELD to user decision.** Commits
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
  floor). **PROVISIONAL pending Guzman et al. 2002** (where the linearised
  theory fails on long multi-rev arcs).
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
   IAC-02-A.6.09 (53rd IAC, Houston; NTRS 20030032208)** — needed to lift the
   #144 primer diagnostic from PROVISIONAL.
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
