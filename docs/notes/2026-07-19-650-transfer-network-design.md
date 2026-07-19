# #650 design read — inter-cycler transfer-compatibility network (2026-07-19)

Analysis-only design specification for the `#650` follow-up implementation task (Sonnet-tier,
mechanical). Defines exactly: eligibility, the compatibility gates, the phase-window model, the
edge cost metric, and the on-disk artifact. Every threshold below is either sourced or explicitly
labelled **convention** (per `[[feedback_digest_not_adoption]]`). No code was written in this
pass; no catalogue rows were touched.

Inputs actually read (not summarized from memory): `#650`/`#645`/`#570`/`#604` bullets in
`data/OUTSTANDING.md`; `docs/notes/2026-07-18-645-creative-strategy-review.md` (item 5);
`data/catalogue.schema.json` (v5.1) + a programmatic survey of all 361 `data/catalogue.yaml`
rows; `src/cyclerfinder/search/tisserand.py` (API), `src/cyclerfinder/model/score.py`
(`taxi_cost_kms`), `src/cyclerfinder/core/constants.py` (PLANETS: `mu_km3_s2`, `radius_km`,
`safe_alt_km`), `src/cyclerfinder/core/satellites.py` (SATELLITES: same fields + `sma_km`),
`data/flyby_altitude_references.yaml`, `data/cycler_network.schema.json` (#570).

## 0. What machinery actually exists (corrections to the task framing)

- **`#604` built NO code.** Its bullet says so explicitly ("No catalogue rows changed, no code
  built") — it was an analysis that compared V∞ regimes using the existing Tisserand relation.
  The reusable machinery is `search/tisserand.py`: `vinf_to_tisserand(body, vinf_kms)` /
  `tisserand_to_vinf` (heliocentric, PLANETS-keyed) and the `linkable`/`linkable_3d` contour
  intersectors. Key reduction: **at a FIXED body, the map v∞ ↔ Tisserand parameter is a
  monotone bijection** (T_P = 3 − v∞²·a_P/μ_S relation), so "same Tisserand w.r.t. X" is
  EXACTLY "same v∞ magnitude at X". The #604-class gate for a shared body therefore reduces to
  comparing the two rows' catalogued `vinf_kms` at that body — no contour machinery needed.
  (`linkable()` is for finding *intermediate linking orbits between different bodies* — out of
  scope here; noted as a possible future extension for pairs that share no body.)
- **`#570`'s schema is for fleet NETWORKS, not pairwise edges** (`member_cycler_ids`,
  `downlink_cadence`, per-member insertion cost), and its own scoping explicitly excluded any
  between-cycler transfer-cost field ("not something to fake into this schema"). #650's output
  is a DERIVED analysis artifact, recomputable from the catalogue — it belongs under
  `data/found/` (the `#317` precedent), NOT in the curated `data/cycler_networks.yaml`
  registry. If adjudicated-real networks ever emerge from the graph, populating
  `cycler_networks.yaml` is a separate later task.
- **`model/score.py::taxi_cost_kms` is Earth→one-cycler insertion cost** (max ‖vinf‖ at Earth),
  not a between-cycler cost — reconfirmed; not reused as the edge metric, but its
  magnitude-only philosophy (catalogue stores no V∞ vectors) carries over.

## 1. Catalogue data-availability findings (measured 2026-07-19, 361 rows)

- Eligible classes `cycler`+`quasi_cycler`: **328 rows** (322 + 6). `precursor_mga` (20),
  `mga_tour` (12), `resonant_po` (1) excluded — see §2.
- `vinf_kms_at_encounters`: present on 351/361 rows, but 59 entries have `vinf_kms: null`;
  after dropping null-vinf entries and the abstract `P1`/`P2` binary-genome bodies,
  **291 eligible rows have ≥1 usable (body, vinf) encounter** — these are the graph's nodes.
- Body vocabulary (usable entries): `E` 289, `M` 250, `V` 22, `Moon` 11, plus moon names
  (Io/Europa/Ganymede/Callisto/Titan/Enceladus/Ariel/Umbriel/Titania/Oberon/Charon) and
  `Pluto`. All resolvable: PLANETS by `code` then by `name`, SATELLITES by name.
- Candidate pairs sharing ≥1 usable body: **31,906** (E 31,626; M 27,966; V 136; Titan 190;
  Enceladus 190; Ganymede 78; Europa 28; Callisto 10; Io 3; each Uranian moon 3). All-closed-
  form per-pair work → the sweep is minutes of compute, not hours.
- **Timing/epoch data — the decisive finding.** The catalogue stores NO absolute encounter
  epochs or phases for ANY of the 322 `cycler` rows: `epoch_locked=false` is a class
  INVARIANT (cyclers "close at any epoch"), `validity_window`/`launch_epoch` are null by
  schema design. The 6 `quasi_cycler` rows (all Uranian moon-pair, #569) carry
  `validity_window {start, end, synodic_period_days, synodic_duty_cycle_pct,
  synodic_boundary_period_days}` — real windows and cadences, but still **no per-encounter
  timestamps**. Periods: 240 rows E-M basis, 15 E-V, 3 VEM-syn; **70 eligible rows have no
  usable `period.years` at all** (mostly moon-tour rows whose sources tabulate members in
  paywalled proceedings; the 6 QC rows carry cadence in `validity_window` instead). Intra-cycle
  encounter timelines (`trajectory.segments[].tof_days`) exist for many rows, but only
  **2 rows** have the full per-segment a/e/tof needed to *derive* an intrinsic encounter phase
  from geometry.
- **Consequence stated upfront (this is #650's honest bottleneck):** a *deterministic*
  phase-alignment verdict between two catalogue rows is IMPOSSIBLE from catalogue data for
  every cycler×cycler pair — the relative phase is simply not recorded (and for the
  epoch-free class, not recordable). The phase-window check must therefore be a
  **statistical model over an explicitly-unknown relative phase** (§5), not a yes/no check.
  This is a data-availability limit, not a physics limit — consistent with
  `[[project_validation_ceiling]]` (past it is new-input-gated).

## 2. Eligibility (Q1 — what rows participate)

A row is a **node** iff ALL of:
1. `orbit_class` ∈ {`cycler`, `quasi_cycler`}. Rationale: an edge represents a steady-state →
   steady-state handoff. `precursor_mga`/`mga_tour` are one-shot (`n_returns=1`) epoch-locked
   paths — a spacecraft "on" one is not in a repeating schedule; precursor→cycler coupling is
   already modeled by `inserts_into`. `resonant_po` has no encounters by definition.
2. ≥1 entry in `vinf_kms_at_encounters` with non-null `vinf_kms` and `body` ∉ {`P1`, `P2`}
   (abstract binary-genome placeholders — "P1" in two different synthetic systems is not a
   shared physical body).
3. No exclusion for `superseded_by`: superseded rows stay in the sweep but every edge record
   carries both rows' `superseded_by` status so the graph consumer can filter.

**Candidate pair** = unordered node pair (i < j by id) sharing ≥1 usable encounter body.
**Edge key** = `(id_a, id_b, body)` — a pair sharing both E and M produces two edge records
(the graph is a multigraph; a real itinerary would pick the cheaper).

Directionality/sense: the catalogue stores V∞ *magnitudes only* — no vectors, no
inbound/outbound encounter geometry (row-level `sense` is an E-M leg label, not encounter
geometry). So no geometric in/out gate is computable, and none is needed under the §3 handoff
model (a periapsis burn reshapes any incoming hyperbola into any outgoing one, up to the bend
limit — recorded as metadata, §3). Both rows' `sense` values are echoed on the edge for the
consumer.

## 3. Gate 2 — V∞/Tisserand compatibility and the ΔV_hop cost (Q2, Q4)

Physical model (the one real design decision, made here): the transfer is a **same-body
powered-flyby handoff** — the taxi departs cycler A on A's approach hyperbola at shared body X
and performs a single impulsive burn at periapsis to leave on cycler B's departure hyperbola.
This is the Oberth-optimal way to change v∞ magnitude, and it is the ONLY cost computable from
magnitude-only catalogue data.

For edge `(A, B, X)`:
- μ_X, R_X, h_safe(X) from PLANETS (`mu_km3_s2`, `radius_km`, `safe_alt_km`) or SATELLITES
  (same fields). Periapsis floor **r_p = R_X + h_safe(X)** (the same floor convention the #426/
  #427 flyby machinery uses; `data/flyby_altitude_references.yaml` is provenance for it).
- For each usable encounter pair (a ∈ A's X-encounters, b ∈ B's X-encounters):
  `dv = | sqrt(v∞_a² + 2μ_X/r_p) − sqrt(v∞_b² + 2μ_X/r_p) |`  (km/s)
- **`dv_hop_kms` = min over (a, b)** of `dv` (best-case handoff; also record `dv_hop_max_kms`
  and which encounter pair achieved the min).
- `delta_vinf_kms = min |v∞_a − v∞_b|` over the same pairs.
- Metadata: per-row max ballistic bend at the floor,
  `bend_max_deg = 2·asin(1 / (1 + r_p·v∞²/μ_X))` for each row's min-v∞ at X, recorded so a
  consumer can see how much direction freedom the flyby offers; and for heliocentric bodies
  (X ∈ {V, E, M}) both rows' `tisserand = vinf_to_tisserand(X, v∞)` values (null for moons —
  `tisserand.py` is heliocentric-keyed; the v∞ comparison carries the identical information).

**Semantics — lower bound, and the right polarity.** With V∞ *directions* absent from the
catalogue, `dv_hop_kms` is a LOWER BOUND on the true single-impulse transfer cost (magnitude
change is the floor of vector change; direction mismatch beyond `bend_max` costs more). So:
an edge classified EXPENSIVE is trustworthily expensive (the negative is clean), while a CHEAP
edge is a *candidate* requiring vector-level follow-up before any claim — the correct polarity
for this program (`[[feedback_orbit_closure_discipline]]`: negatives trusted, positives
gated). Every edge record carries `direction_data: "absent"`.

**Bands (all conventions, labelled as such in code and output):**
- `B0_ballistic_compatible`: `delta_vinf_kms ≤ 0.1` — same Tisserand contour at X to within
  catalogue precision (published vinf values are rounded at the 0.01–0.1 km/s scale, per
  `[[feedback_published_rounded_values_are_display]]`); a pure-flyby handoff is possible *in
  principle*, pending direction data.
- `B1_cheap`: `dv_hop_kms ≤ 0.5` (order of a typical DSM budget — convention).
- `B2_moderate`: `dv_hop_kms ≤ 2.0` (convention).
- `B3_expensive`: `> 2.0`.
Bands CLASSIFY, they do not filter: **every candidate edge is recorded** (census completeness —
the "no cheap edges" negative needs the full distribution, not a survivor list).

**Do NOT rank or combine with any other field into a scalar** — `dv_hop_kms` is the edge
weight; phase feasibility is a separate labelled dimension (§5). A weighted ΔV+wait scalar
would be an unlabelled convention and was considered and rejected.

## 4. Gate 1 — shared body, body identity

Body identity = normalized name via a single resolution function: try PLANETS by `code`
(`E`, `M`, `V`, `Pl`…), then PLANETS by `name` (`Pluto`), then SATELLITES by name (`Moon`,
`Io`, … `Charon`); raise on anything unresolvable (after the §2 P1/P2 exclusion nothing should
be). Catalogue convention (verified on rows): `vinf_kms` at body X is X-relative hyperbolic
excess speed regardless of the row's system, so same-body comparison is apples-to-apples even
across regimes (e.g. a heliocentric E-M cycler's Earth v∞ vs a cislunar row's Moon v∞ are
never compared — different bodies; an E-M row's 6.5 km/s Earth v∞ vs a hypothetical cislunar
row's Earth-relative v∞ would be, correctly, and the ΔV gate would price the huge gap).

## 5. Gate 3 — the phase-window model (Q3)

**What is computable.** Node timing data = repeat period T (from `period.years`·365.25 d; for
the 6 QC rows, `validity_window.synodic_period_days`; rows with neither → `phase_status:
"no_period_data"`, edge still carries all §3 fields) + intra-cycle encounter offsets φ_i for
the shared body (cumulative `trajectory.segments[].tof_days` walk when the segment chain is
complete and consistent with the vinf encounter list, else the n encounters at X are spread
uniformly over T; record which via `phase_timeline: "segments_derived" | "uniform_assumed"`).
The RELATIVE phase δ0 between two rows' schedules is **not recorded anywhere** (§1) — so it is
treated as an explicit uniform unknown, and the phase check outputs *statistics over δ0*, not
a verdict.

**Encounter-coincidence model.** A-encounters at X: t = φ_i^A + j·T_A (j = 0, 1, …);
B-encounters: t = δ0 + φ_k^B + m·T_B. A handoff opportunity = a pair of encounters within a
window w (|t_A − t_B| ≤ w). Window convention (labelled):
- `w_handoff = 2·r_SOI(X) / min(v∞_A, v∞_B)` — both vehicles inside X's sphere of influence
  simultaneously, the strict same-flyby rendezvous condition.
  `r_SOI = a_X · (μ_X/μ_parent)^(2/5)` — a_X from PLANETS `sma_au` (parent μ_Sun) or
  SATELLITES `sma_km` (parent planet's μ). Illustrative: Earth ≈ 9.2e5 km → w ≈ 3.3 d at
  v∞ 6.5; Umbriel ≈ 3.1e3 km → w ≈ 1.9 h at v∞ 0.9.
- Sensitivity grid: w ∈ {w_handoff, 1 d, 10 d, 30 d} (relaxed values proxy for a taxi flying
  a short phasing/connecting arc — NOT modeled in this pass, so relaxed-w columns are
  sensitivity/reporting only, never used in the cheap-edge classification).

**Algorithm (deterministic, no judgment calls).** For each edge needing the phase model:
sample δ0 at 720 uniform points in [0, T_B). For each δ0 and each intra-cycle offset pair
(i, k), the coincidence condition is |(φ_i^A − φ_k^B − δ0) + j·T_A − m·T_B| ≤ w; find the
earliest such time ≤ horizon **H = 100 yr** (convention; `[[feedback_long_runs_acceptable]]`
says don't cap for monitorability, but H here is a *model* horizon — coincidences beyond a
century are not mission-meaningful). Vectorize over j with numpy (centered-mod distance to
the nearest B-encounter); a brute-force O(events²) cross-check on one small constructed case
is a required unit test. Outputs per edge per w:
- `p_align` — fraction of δ0 samples achieving ≥1 coincidence within H;
- `median_wait_years` / `p90_wait_years` — over the succeeding δ0 samples;
- `phase_status` ∈ {`recurrent` (p_align ≥ 0.99 at w_handoff — incommensurate-period
  behaviour: alignment is a waiting-time question, not existence), `phase_locked`
  (p_align < 0.99 — commensurate/near-commensurate periods freeze the offset: aligned-or-
  never, and WHICH is unrecorded → honest label `indeterminate`, with p_align as the prior
  probability), `no_period_data`, `not_computed_dv_gated`}.
- Note for the implementer: this single statistical algorithm handles commensurate and
  incommensurate pairs uniformly — no number theory, no case split. Same-basis E-M pairs
  (T_A = k_A·T_syn, T_B = k_B·T_syn) will naturally come out `phase_locked` with p_align ≈
  (coincidences per superperiod)·2w/T_super; cross-system pairs come out `recurrent` with
  E[wait] ≈ T_A·T_B/(2·w·n_A·n_B).
- **Compute gating (deterministic rule):** run the phase model only for edges with
  `dv_hop_kms ≤ 1.0` (phase feasibility only matters where ΔV is plausibly cheap; B1 + margin)
  OR body ∉ {V, E, M} (all moon-system edges — only ~500, and they are the physically
  interesting subset, see §7). Everything else: `phase_status: "not_computed_dv_gated"`.
- **Epoch-locked special case:** when BOTH rows are `epoch_locked` with a `validity_window`
  (only QC×QC pairs — ≤15 pairs today, all Uranian): first compute the window intersection
  [max(starts), min(ends)]. Empty → `phase_status: "epoch_disjoint"` (a REAL deterministic
  negative, the only kind this data permits). Non-empty → run the statistical model with H =
  the overlap duration, status `epoch_window_overlap`, and echo both rows'
  `synodic_duty_cycle_pct` (feasibility duty cycle further discounts alignment; multiply
  p_align by both duty cycles for the reported `p_align_duty_adjusted`).

**Physical meaning of the statistics (state in the artifact README):** for `recurrent` pairs
the relative phase is partly a *design freedom* (which synodic cycle each cycler is
established on), so p_align/median_wait also read as "co-phasing headroom" for a designed
fleet — the statistic is meaningful in both the wait interpretation and the design one. For
`phase_locked` same-basis pairs the offset is intrinsic trajectory geometry that the catalogue
does not record; deciding it would require re-deriving each trajectory's synodic encounter
phase from full conic geometry, which only 2/328 rows carry data for — out of scope, recorded
as the follow-up's honest limit.

## 6. Cheap-edge classification and the headline census (Q4 cont.)

`cheap_edge = (band ∈ {B0, B1}) AND (phase_status = "recurrent" with median_wait_years ≤ 20
at w_handoff, OR phase_status = "epoch_window_overlap" with p_align_duty_adjusted ≥ 0.5)`.
The 20-yr wait ceiling ≈ mission-lifetime scale — convention, labelled. `phase_locked`
(indeterminate) edges are NEVER counted cheap regardless of p_align — the honest census
separates "cheap and realizable" from "cheap iff an unrecorded phase happens to align"
(reported as its own count: `cheap_dv_phase_indeterminate`).

## 7. Expected outcome, honest risk assessment (Q5)

Pre-computation calibration (closed-form, from the surveyed data):
- **ΔV gate will NOT be the bottleneck for heliocentric pairs**: the ~200 Russell-Ocampo E-M
  rows cluster in published v∞ families, so thousands of B0/B1 same-body edges are guaranteed
  (e.g. `aldrin-classic-em-k1-outbound` × `-inbound`, both E at 6.5 → dv_hop = 0). A permissive
  metric is not the risk here.
- **Phase is the bottleneck, exactly as #645 predicted — but with structure #645 didn't
  name:** (a) same-basis E-M pairs are `phase_locked`: alignment is intrinsic-but-unrecorded →
  indeterminate, a DATA gap, not measure-zero physics; (b) cross-system pairs are `recurrent`
  but wait-limited — E[wait] ≈ T_A·T_B/(2w): two 2.135-yr cyclers at Earth (w≈3.3 d) → ~250 yr.
  So the heliocentric graph's honest headline is likely "many ΔV-cheap edges, ~none
  phase-realizable within mission time, and same-family pairs indeterminate on missing phase
  data."
- **The moon-tour sub-catalogue is where genuinely cheap edges plausibly exist**: short
  periods (days-weeks) beat the waiting-time formula — e.g. two Uranian moon-pair rows sharing
  Umbriel: w ≈ 1.9 h, T ≈ 6-13 d → E[wait] ≈ 1-3 yr, with v∞ mismatches of only ~0.03-0.07
  km/s (B0/B1). The Titan/Enceladus (190-pair) and Galilean blocks are the volume subset —
  though ~45 moon-tour rows lack periods (paywalled member tables), which will surface as
  `no_period_data` and is itself a citable data-gap finding.
- **Verdict on #645's "no cheap edges anywhere" framing**: for the heliocentric majority,
  likely correct and now *quantifiable* (waits of order centuries at the strict window);
  for the moon-system minority, plausibly WRONG — a small phase-feasible cheap-edge
  subnetwork among the Uranian/Saturnian rows is a live possibility. Either outcome is
  citable. The clean-negative trustworthiness demanded by the task holds because dv_hop is a
  lower bound (expensive edges are truly expensive) and the phase statistics are δ0-uniform
  (no optimistic phase assumption anywhere).
- **Scope adjustment for the implementation task (the data-availability finding, stated
  upfront per the task's own request):** the phase check is statistical-only for all
  cycler×cycler pairs; deterministic phase verdicts are possible ONLY for QC×QC
  window-intersection (≤15 pairs) and are otherwise new-input-gated. Do not let the
  implementer discover this mid-build.

## 8. On-disk artifact (Q6)

`data/found/650_transfer_network/` (the `#317` `data/found/<task>_.../` precedent):
- `edges.jsonl` — one record per (id_a, id_b, body): all §3-§5 fields:
  `{id_a, id_b, body, dv_hop_kms, dv_hop_max_kms, delta_vinf_kms, band, vinf_a_kms,
  vinf_b_kms, tisserand_a, tisserand_b, bend_max_a_deg, bend_max_b_deg, r_p_km,
  direction_data, sense_a, sense_b, superseded_a, superseded_b, phase: {status,
  timeline_a, timeline_b, period_a_days, period_b_days, windows: [{w_days, p_align,
  median_wait_years, p90_wait_years}], p_align_duty_adjusted?}, cheap_edge}`.
- `summary.json` — census: node/pair/edge counts, per-band histograms per body,
  cheap-edge count + list, `cheap_dv_phase_indeterminate` count, connected components +
  degree-ranked hubs of the cheap-edge subgraph, counts of every `phase_status`, data-gap
  tallies (`no_period_data`, null-vinf rows), and the parameter conventions echoed verbatim.
- NOT `data/cycler_networks.yaml` (curated registry, different semantics — §0); no
  `catalogue.yaml` writes; no `empty_regions.jsonl` write (that registry is for swept dynamic
  regions; a "no cheap edges" outcome belongs in `data/negative_results.yaml` IF the
  coordinating session decides to register it — adjudication is not the sweep task's call).

## 9. Implementation plan (for the Sonnet dispatch)

- New module `src/cyclerfinder/data/transfer_network.py` (pure closed-form: eligibility,
  body resolution, dv_hop, bend, r_SOI, phase statistics; no integrators, no network access).
- Script `scripts/run_650_transfer_network.py` writing the §8 artifact. NOTE: the
  `tests/scripts/test_scripts_call_preflight.py` AST ratchet will need the same exemption
  category as `#317`/`#606`/`#608` (no region_id/n_points sweep to preflight) — handle it,
  and run `tests/scripts` in the verification set
  (`[[feedback_verify_scope_must_include_tests_scripts]]`).
- Tests `tests/data/test_transfer_network.py`, constructed fixtures per house style:
  dv_hop formula against a hand-computed value; band boundaries; r_SOI formula vs a published
  Earth SOI (~9.2e5 km); phase statistics vs brute-force on a small constructed pair
  (commensurate AND incommensurate case); epoch-window intersection incl. the disjoint case.
- **Positive controls (mandatory before trusting any negative,
  `[[feedback_verify_gauntlet_with_positive_control]]`):** (1) `aldrin-classic-em-k1-outbound`
  × `aldrin-classic-em-k1-inbound` MUST produce an E-edge with dv_hop = 0, band B0;
  (2) the Uranian QC rows sharing a moon MUST produce moon edges with non-empty epoch-window
  intersections (their #569 validity windows overlap by construction). A sweep returning zero
  B0 edges anywhere is broken, not a finding.
- Runtime: ~32k pairs closed-form + phase model on the dv/moon-gated subset — minutes; no
  background run, no checkpointing needed.
- Ruff + `uv run pytest tests/data tests/search tests/scripts -q` before commit; check exit
  code not a summary line (`[[project_pytest_no_final_summary_line]]`).
- No judgment calls remain: every threshold, formula, gate, status vocabulary, and output
  field is fixed above; anything genuinely ambiguous discovered mid-build (e.g. a body name
  this survey missed) should be surfaced, not resolved silently.
