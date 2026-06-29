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
B-plane wall above). The #480 maintenance-ΔV gap is now CLOSED at level-2; only the
level-3 B-plane NLP (Approach C, "weeks, last resort") remains. No catalogue impact.

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
- `feedback_cyclers_pdf_filing_pattern` — when user uploads a PDF to
  cyclers_pdf root, I file + rename into `papers/` per the documented
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
shipped 2026-06-16; Phase 4 V4 HFEM Uranian-system real-eph (#332)
outstanding**, #319 (QP-tori V0-V5). Without #332, even literature-fresh
candidates cannot pass V4.

**Discovery probes (the strategic-answer ideas, ranked by novelty leverage):**
- **#312** Uranus extended sweep (in flight) — 0.062 km/s near-miss is the
  session's highest-novelty-leverage probe; 3D extension is the discovery shot
- **#311** Saturn extended sweep (low-thrust + 3D)
- **#308** Asteroid-leveraging cycler search (fresh ground per #302's
  structural conclusion)
- **#318** Multi-axis joint search (powered × multi-rev × 3D × epoch-locked) —
  strategic-answer keystone; no single published paper has done this joint sweep
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
standing rule — every fetched PDF is filed + committed to the private `cyclers_pdf`
repo (4 filed this session). Roberts-Tsoukkas journal (#251) = no numeric tables
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
