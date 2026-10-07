# #864 — Project feasibility and future review: will we find more novel cyclers, and how?

> **Notice (2026-10-07, `#911`):** the six Uranian (1,1) quasi-cycler rows, including `#312`'s
> `umbriel-oberon-1-1-uranian-quasi-cycler-2026`, were WITHDRAWN from the catalogue on 2026-10-04
> (`#888`: a demanded-turn check found they are not ballistic trajectories; upheld by `#937`), and
> `umbriel-1-2-torus-homoclinic-uranus-2026` was WITHDRAWN on 2026-10-03 (`#882`: its connection is not a
> trajectory of its own model). `europa-3-4-crnbp-torus-jupiter-2026` is `known-class-member`. Any
> statement below that calls these rows novel, validated or catalogued is superseded. This note is a
> dated record and is otherwise not rewritten.

**Date**: 2026-09-05 (AET). **Task**: `#864` (user-requested, analysis-only: no code, no catalogue
writeback, no dispatches from this task itself). **Question as asked**: "Do a project review of the
feasibility and future of this project; will we find more novel cyclers? how will we do it?"

**Method**: a multi-agent review run as two workflows plus coordinator verification. Phase 1: six
read-only evidence readers over disjoint bodies (discovery track record; wall inventory and the fate
of every prior strategy pass; literature density; capability map and compute; catalogue/deliverable
state; open-task ledger) — 371 tool uses, all six returned. Phase 2: six lens strategists
(exploit-the-template, chain orbits, heliocentric return, redefine-the-target, consolidate-and-
publish, program manager), three independent judges (evidence/history, dynamics/novelty,
owner-value), one synthesizer, two adversarial refuter lenses per load-bearing claim, one revision
pass. The first refutation pass was cut short by a session quota (2 of 16 refuters and the revision
completed); the coordinator then verified every load-bearing claim directly (Sec. 6) and the
refutation/revision pass was re-run once the quota reset. Every "VERIFIED" below was read from a
file, re-run, or re-integrated by the coordinator; "INFERRED" is an agent judgment not
independently re-derived.

---

## 1. Answer

**Will we find more novel cyclers?** Probably yes, but only one to three more encounter-bearing
catalogue-class rows over the next quarter, only from three specific places, and none will reach
the spec's "validated" (V5) bar without external review. The moon-system template that produced
every hit so far is spent for encounter-bearing cyclers except one unsearched cell. Novel
*dynamical objects* (tori, homoclinic connections) will keep arriving at roughly one in five pairs
if the existing pipelines are pointed at untouched pairs — but they are not cyclers, and the owner
should decide explicitly whether they count.

**Where.** In expected-value order:

1. **The heliocentric Venus-Earth-Mars lane the spec named on day one and the program abandoned on
   2026-06-07 on a misdiagnosis.** This review's single most consequential finding (Sec. 5): the
   Jones-Hernandez-Jesick 2017 (AAS 17-577) VEM anchor, xfail for three months, reproduces at
   patched-conic level in 0.6 s once each Lambert leg is allowed 3-8 revolutions. Every June attack
   (#110/#120/#122/#133) enumerated at most one revolution on one leg. So the "seeding/basin wall"
   as applied to VEM was a topology-enumeration gap, those negatives are void for this family, and
   Jones's own stated scope (1-2 synodic periods, at most six flybys per cycle) leaves the 3- and
   5-synodic (19.2/32-yr) and >6-flyby ballistic VEM classes unsearched by anyone. Caveats from the
   refutation pass (Sec. 7a.5): the probe reproduces Jones's stage-1 broad-search chain, not a
   converged member (two interior mismatches exceed his 200 m/s tolerance); the heliocentric ceiling
   this quarter is V3, not spec-V4, because the installed GMAT is a Linux binary and its lane never
   produced a V4 row; and whether a member of a class Jones excluded on practicality grounds counts as
   novel is an owner policy decision under the project's own #577 precedent.
2. **The Russell-Strange 2009 one-working-node class** (#819 made concrete) at its three bend-capable
   working nodes: Uranus (Titania/Oberon flyby body, Miranda passive target — marginal at the #818 2%
   parasitic-turn gate and ~4.3 deg inclined; Ariel/Umbriel inadmissible as passive targets), Saturn
   (Titan flyby body, icy moons as passive targets — a census, since Russell-Strange own Titan-
   Enceladus), and Neptune (Triton flyby body, Proteus passive — Sec. 7a.4). The method is published
   and positive-controllable against 32 catalogue rows (TitEnc#235), R-S never touched Uranus or
   Neptune, and #817 ruled the class admissible and unswept. But single-moon Oberon 4:5 free-returns
   ARE published (Kumar-Anderson AAS 24-288, ~0.74x Oberon SOI), and under #817's own reading of the
   #577 standard an R-S-class member at a new pair is known-class, so novelty here is ~1 in 10 unless
   the owner's novelty policy says otherwise. Hard-gated on first reading the unfiled UOP-era Uranian
   tour papers (Landau 2025, Ellison 2025, PSJ 2026, AAS 25-668).
3. **Planet-moon CR3BP symmetric periodic orbits at Neptune-Triton — already closed during this
   review.** The Miceli-Bosanac 4:5 saddle passes 6,629 km from Triton (0.55x its 11,966 km sphere of
   influence; coordinator integration, Sec. 6), and feeding #781's two on-axis homoclinic seeds to the
   existing 1-DOF C-pinned corrector closed TWO symmetric periodic orbits (periods 2.26 and 2.86 base
   periods, own periselenes 0.58x and 0.46x SOI, Radau-confirmed; two refuters and the coordinator
   independently, Sec. 7a.1). Under the project's own #811/#855 per-member periselene-vs-SOI rule they
   are `cycler`-class rows today, with no schema change. Novelty is thin (theorem-generic accumulation
   orbits; "first computed at Neptune-Triton") and pending an ESM-family collision check and two
   unread Neptune papers; ceiling V1. A second, cheap one-working-node cell at Neptune (Triton flyby
   body, Proteus passive target: #599's 104 residual survivors were never re-gated under the #818
   passive rule) was also found by the refuters (Sec. 7a.4).

**Probabilities (after the refutation pass, Sec. 7a).** ~70% for at least one novelty-claimable
cycler-class row at any tier in 12 weeks — but almost all of that is the two V1 Neptune-Triton orbits
the refuters already closed (Sec. 7a.1), whose novelty is thin; ~25-30% for one at V3+ (the
heliocentric lane tops out at V3 this quarter because the GMAT spec-V4 route was found not to exist
on this machine, and the Uranian R-S novelty collapsed to ~1 in 10 under the project's own #577/#817
standard). Every prior hit was narrowed under adversarial review — as happened again inside this
review; expect "first computed in an unsearched class of a known concept", not a new mechanism. A
single owner decision — the NOVELTY POLICY on author-excluded classes and known architectures at new
systems (Sec. 10) — now decides whether the remaining discovery items can produce anything the owner
may call new, as opposed to census rows.

**How.** Run the template — a published, positive-controlled method pointed for the first time at
an under-published, bend-capable cell — *enumeration-first rather than corrector-first*, under the
pre-registration discipline #858 already discovered (positive control at the end target, numeric
kill criterion, cost tied to a measured unit cost, writeback path that exists today), and put the
writeback home and registry stamps in place before dispatch. The 14-day stop since 2026-08-22 is a
real signal that the *campaign menu* was exhausted, not the project's value.

---

## 2. Vitals (VERIFIED)

- First commit 2026-05-31; 14 weeks old. 2155 commits: Jun 1345 / Jul 644 / Aug 150 / Sep 0. Last
  commit 2026-08-22 06:48 AET; zero commits, notes or task registrations Aug 23 - Sep 4. New task
  ids per ISO week fell from ~100-160 (June) to ~25-33 (Aug 4-22) to 0.
- 864 task numbers; ~16 analysis-only strategy passes since the 2026-06-13 discovery spec, ~60
  shortlist items, essentially all executed. 591 dated notes; ~265 corpus papers (private repo).
- 147k LOC src (155 search modules), 113k LOC tests (3991 test functions), mypy --strict, ruff,
  pre-commit ratchets. Machine: Apple M3 8-core / 16 GB, no fp64 GPU (#697); June-era unit costs came
  from a 16-core Linux box. Self-hosted CI runner on the same Mac since 2026-08-09.

## 3. What the project has actually found (VERIFIED from rows and notes)

- Catalogue 399 rows: literature 368 (92%), this-project known-class 23, discovered 8 (2%).
  validation_level: null 240 / V0 99 / V1 44 / V2 8 / V3 2 / V4 6 / **V5 0**. By the spec's own trust
  ladder (docs/spec.md sec 14) only V5 plus catalogue and literature miss may be called a discovery.
- The 8 project-originated rows: six Uranian moon-pair (1,1) symmetric-closure quasi_cyclers
  (#312/#563/#569; in-house V4-strict via Python+SPICE URA111, not the spec's independent-codebase
  V4; 62-79% synodic duty cycles 2000-2083); `umbriel-1-2-torus-homoclinic-uranus-2026` (CCR4BP, V1,
  never within 181,000 km of a moon, ~1%-duty real-eph near-miss window recurring 10/10 epochs);
  `europa-3-4-crnbp-torus-jupiter-2026` (N=5 CRNBP, V1, bare torus, no manifold/stability/transfer
  computed). Two of the three "novel findings" are not cyclers; schema v5.3/v5.4 were minted to
  admit them. `our_status` is absent on all 8. The six flagship rows lack a `primary:` key and so
  default to heliocentric under data/README's rule (data defect, VERIFIED).
- Unwritten-back findings: #781 Neptune-Triton 4:5 homoclinics (the only novel result of the
  #753-#786 arc; Neptune-Triton has zero catalogue rows; the #856 registry requires a catalogued
  endpoint); #767 Saturn-Titan homoclinics; #782's Saturn-Titan 3:4<->6:5 chain periodic orbit.
- Hit rates: 3-4 of ~22 method families produced a novel-claimed arc; within productive families
  1-in-4 to 1-in-8 systems; roughly one written arc per month Jun-Jul, zero in Aug and Sep. Every
  method built as a *new method* (ER3BP isolated families, Floquet branches, 3D lifts, binaries,
  cross-system SE-EM, WSB, GAIO, DA/HOTM, Sobol+n-body, precursor MGA, asymmetric closures, seedless
  spectral, deflation, generative seeds, homotopy, W-Z proofs, CLV extraction) returned empty or
  known-class. **The template every hit followed**: a published, positive-controlled method applied
  for the first time to an under-published, mass-rich, commensurate outer-planet moon system or
  pair (Uranus x2, Jupiter x1, Neptune x1); each method's hit came on its first such application;
  zero hits in Earth-Mars, Earth-Moon or Sun-planet regimes, where ~60 of the 96 registered
  negatives sit. Every hit was overclaimed at first framing and narrowed by adversarial review.
- Negative-results registry: 96 method-versioned records, but the last stamp is 2026-08-10 and the
  CCR4BP negatives (#695/#696/#703/#716), the connection-lane negatives (#759/#780/#783/#786), the
  #861 Oberon gate and #791 Stage 0 are unstamped (Jun 62 / Jul 32 / Aug 1). Stamps by month tell
  the same story as commits.

## 4. Walls, and which are physics

Eighteen walls are catalogued in the Phase-1 output (walls-frontier reader). Physics walls that
should not be re-attacked: W5 mass-limited bend gate (small moons, asteroid moons, Phobos/Deimos,
Proteus, Amalthea, Miranda-as-flyby-body; several interval-certified); W12 Galilean hub energy
(V-inf 4-5 km/s links all 12 ordered Galilean pairs at once, so energy filters cannot prune tour
sequences); W15 incommensurate forcing (only Laplace-locked Jovian and near-2:1 Uranian pairs admit
a stroboscopic N>=4 reduction); W16 the equal-leg asymmetric-closure continuum; the SE<->EM core of
W2. Method walls with cheap, named, untried steps: W9 seeding-topology misidentification on systems
without a published IC table (#862 crossing-index sweep untried; the #861 addendum showed the true
Oberon 4:5 saddle is reachable from hc=3); W11 chain closure (the #687 multiple-shooting corrector has
no Jacobi row — VERIFIED in code; the Jacobi-pinned extension #858 called "the single highest-
leverage build for the whole combinatorial program" was never registered); W13 writeback has no
home for connection/chain findings — **but see Sec. 5.2: the existing `cycler`/`resonant_po` classes
plus the #855 per-member rule already provide one**; W18 CPU headroom never taken (memoisation,
numba kernels). W1 (family-selection/basin, the #388 lane) stands for shooters seeded out of basin
but — Sec. 5.1 — was misapplied to VEM.

Recurring failure patterns (all VERIFIED with examples): unanchored scale-ups die on seeding, not
compute (#776, #859, #861); cost models off by orders of magnitude in both directions (#859 measured
0.5-0.7 CPU-h vs 50-100 estimated; #790 "minutes per closure" vs 5 dispatches / 8 days for one
chain); headlines flip under review in both directions; closed work gets re-proposed (#792 duplicated
#680; this review's own exploit-template strategist re-proposed the #600 near-miss refinement that
#663 had already run and killed on bend — caught by all three judges; root cause: the #600 stamp
carries no #663 cross-reference); discovery outputs lack a schema home; analysis passes have
outnumbered rows.

## 5. This review's own findings

### 5.1 The Jones VEM anchor was never unreachable — the June attacks never enumerated its topology

The spec's headline novelty region ("VEM ... where genuinely new orbits may exist") has as its only
published anchor Jones-Hernandez-Jesick 2017, AAS 17-577 (Tables 2/3: EMEVVE outbound and MEEVEM
inbound 12.8-yr triple cyclers, per-encounter V-inf 2.4-7.0 km/s, flyby altitudes). The gate
`tests/test_vem_rediscovery.py::test_jones_vem_ballistic_rediscovers_sourced_multiset` has been xfail
since 2026-06-06 after four attacks: #110 (dense 2816-point epoch x topology scan, floors 17.9/18.5
km/s, 0 bend-feasible), #120 (inclination, refuted), #122 (vector residual, refuted), #133 (n-body
flyby-propagation shooter seeded from a conic near-miss survey that only found the 29-33 km/s basin).
The June diagnosis was "leg-topology / seeding problem: multi-arc-per-leg seeding is the front-
runner", and no task after #142 touched it.

The heliocentric strategist wrote a 0.6-second probe (scratchpad only, not committed:
`probe_jones_truth_multirev.py`) that, at the PUBLISHED encounter dates, solves every prograde Lambert
branch with 0..8 revolutions on each leg against DE440 and picks the per-leg branch chain minimising
the summed interior |v-inf| mismatch. Published V-inf and altitudes are *outputs compared*, not
inputs. **Coordinator re-run 2026-09-05 14:11 AET (VERIFIED)**:

| Case | Summed interior mismatch (9 flybys) | Max per-encounter mismatch | Altitude check |
|---|---|---|---|
| EMEVVE (Table 2), MAX_REVS=8 | 1.136 km/s | 0.229 km/s (2035-05-22 E) | Venus 684 vs 684 km; all 9 interior above surface |
| MEEVEM (Table 3), MAX_REVS=8 | 0.803 km/s | 0.210 km/s (2035-11-17 E) | Mars 251 vs 249 km; all above surface |
| EMEVVE, MAX_REVS=1 (June-style control) | 42.6 km/s | 14.1 km/s | five negative altitudes |
| MEEVEM, MAX_REVS=1 (control) | 80.8 km/s | 15.8 km/s | five negative altitudes |

Chosen chains: EMEVVE n0s,n1h,n5h,n8l,n1h,n0s,n1h,n1h,n7h,n3l; MEEVEM n0s,n1l,n4l,n5h,n1h,n0s,n3h,
n0s,n6h,n1h. The long legs need 3-8 revolutions. The identification is specific: MAX_REVS=6 breaks
EMEVVE (legs 3 and 8 need n8l/n7h); MAX_REVS=10 returns the identical chains (refuter re-run).

Why June missed it (VERIFIED by reading the code): `scripts/hunt_vem_ballistic.py::_topologies`
enumerates all-direct plus ONE 1-rev leg at a time (11 topologies); the gate test hard-codes
`per_leg_revs=(0,)*n_legs`; `nbody/shooter.py::near_miss_survey` defaults `branch_topologies` to zero
revs. A multi-rev Lambert leg IS a single ellipse, so the June "single-ellipse-per-leg ... needs
multi-arc-per-leg" diagnosis was a misread. Per `[[feedback_bugfix_invalidates_past_searches]]`, the
#110/#120/#122/#133 VEM negatives are void for this family. No empty-region stamp mentions Venus or
VEM (0/96, VERIFIED). Jones p.3 (PDF read by a refuter): "attention is restricted to families with 1
or 2 synodic periods in a cycle, and a maximum of six flybys per cycle"; p.8: none at one synodic
period, "thousands" at two. Citation indexes show no later VEM sweep (OpenAlex 1 citer, Semantic
Scholar 2, both non-VEM; INFERRED for conference-paper under-coverage).

Caveat carried honestly: a 0.01-0.23 km/s conic mismatch is what Jones's own stage-1 tolerated
(<=200 m/s) before SNOPT n-body correction; the anchor is reproduced at conic level, and the n-body/
GMAT rungs (V3/V4) still have to run. Unintended intermediate flybys on 5-8-rev legs must be filtered
(no such filter exists in src).

### 5.2 Chain orbits and resonant POs already have a catalogue home

The brief (and the #858 review) said no orbit_class exists for a resonance-transition chain orbit.
Both strategists and all three judges found, and the coordinator VERIFIED, that #811 (Vaquero EM
writeback) and #855 (C=3.13 adjudication) set orbit_class PER MEMBER by the orbit's own periselene vs
the moon SOI, with two integrators: inside the SOI -> `cycler`, outside -> `resonant_po`. Coordinator
integration of the Miceli-Bosanac ESM4 4:5 saddle (x0=-0.969056422016, ydot0=-0.126767119783,
C=2.987089791658, T=30.398418802755, mu=2.0895e-4, DOP853 rtol=atol=1e-12): return error 2.4e-9,
closest Triton approach **6,628.6 km = 0.554x SOI (11,966 km)**. So W13 closes for #781 at zero CPU
via a literature-sourced `cycler` row for the 4:5 saddle carrying the homoclinics as provenance, and a
closed Neptune-Triton chain inheriting those passes is a cycler-class row without a schema change.
The #782 Saturn-Titan chain (~20,100 km from Titan, 0.46x SOI per two agents) is likewise a `cycler`
row — a known-reproduction of Vaquero's family.

### 5.3 Corrections to the Phase-1 brief found by the strategists/judges (all VERIFIED)

- "Powered/low-thrust cyclers never swept" was wrong for moon systems: #464 (Sims-Flanagan
  low-thrust releg, 12.03 vs 13.18 km/s/cycle, 3.4x over the ceiling) and #465 (160 combos; Uranus/
  Neptune stamped powered-empty) ran. But those Uranus/Neptune stamps probed V-inf 4-15 km/s only,
  against a Uranian natural scale of 0.2-1.3 km/s — they are band-conditional and must not block a
  low-V-inf Uranian search (amend in the stamp).
- "#695 Io-Europa near-miss never re-run after #702" was misleading: #702 cross-checked every saved
  #695 candidate and all fail the off-torus threshold at 562-823 km independent of the bug; only an
  energy/family continuation is a new experiment.
- "#388-lane negatives not yet re-run after #849" is stale for the DSM lane: #849 re-ran
  `close_row_dsm` under both posings with zero verdict flips.
- "#587 cannot be found" — #587 is registered and CLOSED; the defect is a `data_gaps.todo_ref`
  pointing at a closed task while the gap persists.
- The #600 3-moon near-miss is dead (#663 found the exact closure, bend ~0.83 deg); the stamp lacks
  the cross-reference.

## 6. Load-bearing claims and their verification status

| Id | Claim | Status |
|---|---|---|
| LB1 | Multi-rev Lambert chain reproduces both Jones cyclers at published dates | VERIFIED by coordinator re-run + one refuter re-run (MAX_REVS 6/8/10 robustness) |
| LB2 | Every June VEM attack enumerated <=1 revolution on <=1 leg | VERIFIED by coordinator (three files read) |
| LB3 | Jones restricted scope to 1-2 synodic / <=6 flybys; no later VEM sweep published | Scope VERIFIED from the PDF by a refuter; absence of later work INFERRED (citation indexes + searches) |
| LB4 | Neptune-Triton 4:5 saddle periselene inside Triton's SOI -> `cycler` class by #811/#855 | VERIFIED by coordinator integration (6,628.6 km, 0.554x SOI) and rule read |
| LB5 | Symmetric chains accumulate on on-axis homoclinics; #781's seeds close with the 1-DOF corrector | HOLDS as a VERIFIED COMPUTATION (two refuters + coordinator: two orbits closed, Sec. 7a.1) — but the registered recipe (hc=13/15) is WRONG and returns 0/6 |
| LB6 | No published periodic Uranian moon cycler; R-S Table 1 excludes Uranus; UOP papers are one-way | REFUTED IN PART (0.7): single-moon Oberon 4:5 free-returns are published (Kumar-Anderson AAS 24-288, ~0.74x SOI); only the TWO-moon R-S class is unpublished; R-S scope VERIFIED from the PDF; UOP one-way VERIFIED for PSJ 2026 text, INFERRED for Landau/Ellison bodies |
| LB7 | GMAT B-plane lane is the only spec-V4 route and converged on Aldrin/S1L1 | REFUTED (0.85 both lenses; coordinator VERIFIED): binary is Linux x86-64 ELF; June figures retired by #176; both rows FAILED the V4 predicate; no spec-V4 row exists; heliocentric ceiling this quarter is V3 |
| LB8 | Every cycler hit followed the template; ~15 new-method families returned zero rows; Uranus closed on three axes | REFUTED IN PART (0.65): sub-claims hold, but a second unswept one-working-node cell exists at Neptune (Triton working node, Proteus passive; #599's 104 survivors never re-gated under #818 — VERIFIED from the stamp), and the inductive base is one arc |

## 7. Ranked roadmap (synthesis, adjudicated by three judges; revisions from the refutation pass in
Sec. 7a)

| # | Type | Item | Effort (days) | Kill criterion | Writeback / ceiling |
|---|---|---|---|---|---|
| 1 | enabler | Week-1 integrity bundle: `primary: Uranus` on 6 rows; `our_status: candidate-novel` + `discovery_run` block on 8 rows; stamp the ~10 unstamped negatives with reopen conditions; cross-ref #663 on the #600 stamp; band caveat on #465 stamps; record the VEM topology-gap invalidation; register G1; mark #789 SHELVED; fix #790 text | 3-4 | none (cap ledger hygiene at 2 days) | existing fields |
| 2 | reproduction | Jones VEM anchor to V3/V4: seed `ballistic_correct` at published dates with the probe's rev/branch chains, fix the gate to take rev/branch as a sourced-derived input, n-body shoot, GMAT B-plane targeting; flip the xfail only on a genuine pass | 2-4 | conic seed does not close <0.1 km/s on either row after 1 day; n-body cannot reach <=200 m/s in 2 CPU-h/member | existing `jones-2017-*` rows; spec-V4 |
| 3 | discovery-cycler | Jones-style itinerary-growth enumerator (near-Hohmann seeds, per-leg ToF grids, all revs 0-8 x low/high, v-inf continuity <=100-200 m/s, altitude window, unintended-flyby filter), positive-controlled by re-finding Tables 2/3 from Hohmann seeds and the 6.4-yr negative; then sweep the 19.2/32-yr and >6-flyby classes over 3 opportunities | 8-11 | control fails to re-find either 2022 member; sweep yields zero chains within 200 m/s -> stamp | `cycler`, multi-arc, analytic-ephemeris, epoch_locked; spec-V4 via GMAT |
| 4 | discovery-cycler | Neptune-Triton: seed `correct_symmetric_fixed_jacobi` (via `attempt_chain_closure_symmetric`) with #781's on-axis seeds (x=1.16933872, hc=13; x=-1.38561105, hc=15; k+-2) and #767's unused Saturn-Titan on-axis seed; write back the 4:5 saddle as a literature `cycler` row carrying #781's homoclinics, plus the Saturn-Titan 3:4 and #782 chain rows | 1.5 | all seeds fail, drift >1e-2, fail Radau, or close onto a multiple cover -> stamp | `cycler` by #855; V1 |
| 5 | enabler | Close the corpus gap: file and digest the UOP-era Uranian papers (PSJ 2026 10.3847/PSJ/ae680c, Landau 2025, Ellison 2025, AAS 25-668) plus index entries for Baresi-Owen 2026, Acta 2026 Saturnian tour, Bellome 2023, the 2024 Jovian review, Brown et al. 2024; topology-labelled anchors into `literature_check.py`; fortnightly watch | 3-4 (+0.25/fortnight) | none; if a UOP paper has a repeating two-moon Uranian cycle, item 6 becomes a reproduction and the preprint section is cut | CORPUS_INDEX + anchors (the documented-literature-review half of V5) |
| 6 | discovery-cycler | Uranian Russell-Strange one-working-node campaign (#819): `uranus_system()` in `moon_cycler_genome.py`, R-S itinerary enumerator, gate on TitEnc#235 (lifts 32 R-S rows V0->V1), then Titania/Oberon flyby body with Miranda's deflection modelled, then Ariel; V2->V3->V4-strict URA111 lane + daily duty scans | 8-14 | TitEnc#235 not reproduced in 5 build-days; zero itineraries <=8 pass bend + #818 self-consistency -> method-conditional stamp | `quasi_cycler`, primary Uranus; in-house V4-windowed |
| 7 | enabler | Campaign charter + lint (deliverable class; positive control at end target or UNANCHORED; non-cycler-vocabulary literature precheck; registry query; measured unit cost or GUESS; numeric kill; writeback path that exists today; stamping plan); `data/found_index.jsonl` + registry query; #795 CLI + detached supervisor with lock file and heartbeat, proven by a 1-h SIGKILL/duplicate-launch dummy run | 3-3.5 | dummy run drops or duplicates a cell -> no real campaign runs | data/campaigns/, scripts/ops/ |
| 8 | enabler | G1 Jacobi-pinned multiple-shooting corrector (algebraic elimination of one velocity component at the energy node; weighted-row variant as cross-check; arclength wrapper), control ladder Arenstorf -> re-close #782 -> Earth-Moon 3:1<->2:1 asymmetric from #840's legs; symmetric on-axis chain pilot over published alphabets; asymmetric tail only if the EM control passes | G1 2-3; pilot 3-5; tail 5-8 (flex) | fails to re-close #782 -> bug; EM control fails after seeding variants -> asymmetric chain is the wall, stamp and close #790's general form | rows per item 4; V1 |
| 9 | consolidation | Honest spec-V4 for the six Uranian rows via a ~200-line tudatpy driver (Uranus J2 + five moons from URA111), or GMAT body-authoring fallback, or an explicit "V4-internal" relabel | 3-5 | Tudat cannot ingest URA111 in 2 days and GMAT authoring exceeds 3 more -> relabel | validation_level + evidence note |
| 10 | discovery-cycler (flex) | Titan-Iapetus node-locked inclined quasi-cycler (the #575 stamp's own named, never-built reopen condition) with a low-maintenance fallback; fund only if item 3 or 6 kills early | 4-6 | 0 node-locked roots AND best TCM >300 m/s per 7 cycles -> replace the #575 stamp | `quasi_cycler`, primary Saturn; in-house V4 |
| 11 | publication | LICENSE + CITATION version/DOI + Zenodo release of catalogue/registry/connections/errata; one-sentence spec sec 14 amendment permitting a candidate-framed preprint; arXiv preprint (Uranian family, negative-results registry method, errata, Jones reproduction, any new rows); clean-clone reproduction package; outreach to the Kumar/Anderson, Russell-Strange and UOP groups; notify the 7 errata authors | 15-22 (weeks 9-12) | owner declines a public license -> stop at the DOI step | the only route to V5 |
| 12 | discovery-object (idle-time) | CCR4BP energy-continuation of the Io-Europa (16.7 km) and Umbriel-Titania-secondary (43 km) near-misses; untouched Uranian pairs Ariel-Umbriel, Ariel-Titania, Umbriel-Oberon (Titania-Oberon excluded: Kumar's own pair); #862 -> Oberon re-gate (>=4/6 at a 2-h budget) before any atlas alphabet; every negative stamped the same day | 8-15 spread | per pair: no guard-passing candidate at seed energy and across the C-loop -> stamp; lane: 0 hits after 3 pairs -> close | `torus_homoclinic` v5.3; V1; not cyclers |

Owner-days total 55-82 against ~60-72 available in 12 weeks; items 10, the G1 asymmetric tail, and
item 12 are the slack.

**Twelve-week calendar (decision gates).** Week 1 (Sep 7-13): items 1 and 4, Jones conic seed. Gate A
Sep 13: 6 rows carry `primary: Uranus`, 8 carry `our_status`, registry >=106 with #663 on the #600
stamp, Neptune-Triton chain verdict recorded, Jones conic chain closed or bound recorded. Week 2: item 2
n-body/GMAT, item 3 build starts, item 5 digests. Gate B Sep 20: Jones at V3/V4 or item 3 NOT
chartered; UOP verdict written (item 6 GO or REPRODUCTION); charter lint live. Weeks 3-4: item 3
positive control, item 7, item 6 design + TitEnc#235 gate, G1. Gate C Oct 4: enumerator control
passes -> sweep dispatched; TitEnc#235 reproduced -> Uranus GO; supervisor proven; G1 EM-control
verdict. Weeks 5-6: VEM sweep (detached, checkpointed), Uranus enumeration + URA111 lane. Gate D Oct
18: survivor tallies, Uranus V4-windowed or stamp, pilot hit-rate table, Oberon re-gate verdict.
Weeks 7-8: adjudication (Fable pass on every hit, class-appropriate literature checks, GMAT V4 on top
VEM survivors), rows written or regions stamped. Gate E Nov 1: rows-vs-stamps tally by class;
preprint scope decided. Weeks 9-10: item 9, item 11a, site sync. Gate F Nov 15: DOI resolves, V4
label honest, site count equals catalogue count. Weeks 11-12: preprint, reproduction package,
outreach, errata notifications. Gate G Nov 29: next-quarter decision. Miss a gate by more than a
week and the following wave is not dispatched; slack goes to consolidation, never to another
analysis pass.

### 7a. Revisions from the refutation pass (16 refuters, 2 lenses per claim; coordinator re-verified)

Five of the eight claims survived intact with wording corrections; LB5 became a verified
computation; LB6 and LB8 were refuted in part; LB7 was refuted outright by both lenses. The
consequences, each re-checked by the coordinator before adoption:

1. **Neptune-Triton chain orbits are CLOSED, not predicted (LB5).** Both refuters ran the shot; the
   coordinator re-ran it (2026-09-05 14:53 AET, `attempt_chain_closure_symmetric`, the #782 path,
   C=2.987089791658). PRIMARY seed x=1.16933872 with `ydot0_sign=-1`, `half_crossings=5`: converged in
   4 Newton iterations, crossing residual 2.86e-14, x0=1.1693356584, T=68.7485 (2.26 base periods,
   not a multiple cover), |lambda|=1.66e5, Radau full-period closure 7.6e-9, dC 4e-14, closest Triton
   approach 6,897 km = **0.58x SOI**. SECONDARY seed x=-1.38561105 with `ydot0_sign=+1`, `hc=4`:
   residual 3.30e-14, x0=-1.3853113158, T=86.898 (2.86 base periods), |lambda|=126.6, Radau closure
   6.8e-11, closest approach 5,455 km = **0.46x SOI**. Every other index tried (hc 7-17 including the
   literal manifold indices 13/15 and k+-2) did NOT converge in 8 iterations. **The recipe registered
   in #868 (hc=13/15, k+-2) returns 0/6 and would have stamped a live lane empty** — a false negative
   caught only because the refuters ran the computation instead of arguing about it. Corrected
   recipe: pick the `ydot0` sign under which the seed shadows the base orbit, take the FIRST y=0
   crossing with |xdot|<1e-2 near a base-orbit perpendicular point, one member per seed. Both orbits
   are `cycler`-class by their OWN periselene (the #855 rule requires the closed orbit's, not the
   parent's). They are theorem-generic (symmetric POs accumulating on a symmetric homoclinic), so the
   defensible novelty is "first computed symmetric periodic orbits accumulating on a Neptune-Triton
   homoclinic"; they must be checked against Miceli-Bosanac's ESM resonant families (T/2pi = 10.94
   and 13.83 sit near integers) and two unfiled Neptune papers (Campagnola et al. ISSFD 2014;
   Miceli-Bosanac AIAA 2024-1280) before any writeback. Ceiling V1 (unstable CR3BP, no Neptune SPICE
   lane). Held for adjudication per orbit-closure discipline; NOT written back by this task.
2. **The GMAT spec-V4 route does not exist today (LB7 refuted, 0.85 both lenses; VERIFIED).** The
   installed `~/GMAT/R2022a/bin/GmatConsole` is a Linux x86-64 ELF binary and cannot run on this
   arm64 Mac. The June convergences were retired by #176 (Aldrin targeted the wrong flyby; S1L1's 7.29
   km/s was a re-anchoring artifact); both rows FAILED the two-part V4 predicate and V4 was held OPEN
   (Aldrin is V2, S1L1 V3 today). The lane is a per-flyby SOI-to-SOI turn check that never propagates
   a heliocentric leg. **Heliocentric ceiling this quarter is V3** (REBOUND/IAS15 + DE440), not
   spec-V4; a spec-V4 verdict needs a GMAT re-host (Rosetta 2 or an x86_64 container, 0.5-1 day) plus
   a never-designed continuous multi-cycle chain lane (GUESS 3-5 days). Item 9 is reordered: GMAT
   re-host -> run the existing untested `gmat_v4_uranus_generate.py` -> Tudat -> explicit
   "V4-internal". The label "independent codebase, same SPK" is the honest one.
3. **Published single-moon Uranian free-returns exist (LB6 refuted in part, 0.7).** Kumar-Anderson
   AAS 24-288's Oberon 4:5 family passes ~7,177 km from Oberon (0.0123 nondim in #861's own
   mu-continuation data) = ~0.74x Oberon's SOI, so it is a literature `cycler`-class row under the
   #855 rule and should enter the catalogue as such; the UOP inclination-cranking sequence publishes
   repeated Titania-Titania legs. What remains unpublished is the TWO-moon (flyby body -> distinct
   passive target) Russell-Strange class. Ariel/Umbriel are inadmissible passive targets (35-56%
   parasitic turn); Miranda is marginal (2.6-3.0% at 100 km vs the #818 2% gate) and ~4.3 deg
   inclined; and #817 itself says an R-S-class member at a new pair is not novelty-bearing under the
   #577 standard. **Item 6's novelty expectation falls from ~1 in 3 to ~1 in 10**; it survives for its
   policy-independent value (32 R-S rows V0->V1 via the TitEnc#235 control, a Saturn Titan-flyby
   passive-target census, the `uranus_system` genome).
4. **A missed one-working-node cell at Neptune (LB8 refuted in part, 0.65; VERIFIED from the stamp).**
   #599's stamp records 104 of 1024 candidates passing the residual gate and failing ONLY at
   Proteus's bend (0.005-0.3 deg) while Triton's own bend was 1.4-30+ deg; `reopen_condition` is
   None; the #818 passive-node gate (which Proteus passes trivially at 0.1-0.4% parasitic turn) was
   built after #599 and never applied. Stage 0 (re-gate the 104 survivors) costs seconds; Stage 1
   needs the node-locked inclined construction (Triton is retrograde, ~157 deg to Proteus) that the
   #575 stamp names as never built — the same build Titan-Iapetus needs, so item 10 becomes a shared
   one-working-node node-locked construction for two cells. The one-working-node class is therefore
   open at THREE working nodes (Uranus, Saturn-Titan passive targets, Neptune-Triton), not one, and
   LB8's inductive base is honestly ONE arc (#312->#569) whose positive control post-dated the hit.
5. **LB1/LB2/LB3 hold, tightened.** The probe reproduces Jones's STAGE-1 broad-search chain, not a
   converged member: two interior mismatches (0.229, 0.210 km/s) exceed Jones's own <=200 m/s
   stage-1 tolerance; the largest deviation from a published V-inf is 0.27 km/s; 16 of 18 interior
   altitudes agree within 5.5% but two are off by 16-23%; a +-0.5 d date-rounding test left EMEVVE at
   0.275 km/s. P(conic gate flips) is revised from ~85% to ~60-65%. Per-leg revs are (0,1,5,8,1) /
   (0,1,4,5,1) — the gap is JOINT multi-rev assignment across 4 of 5 legs, not one long leg. Jones's
   "thousands" are stage-1 chains (2 members tabulated fully ballistic). Two of 96 registry stamps do
   mention Venus, as an intermediate body in precursor-MGA insertions — the correct statement is "no
   VEM-cycler or Earth-Venus cycler stamp exists". The catalogue's `k=3` (E-M synodic basis;
   `vem-emeeve-3syn` = 6.4 yr = Jones's EMPTY class) collides with "3-synodic" in Jones's T_syn units
   (19.2 yr): any registration must state k in T_syn units and years. **Under the project's own
   #577/#578 precedent a member of a class the source authors EXCLUDED on stated practicality grounds
   is known-class unless the owner writes a novelty policy** distinguishing "author-declared scope
   exclusion" from "not enumerated" — this single decision now gates the novelty (not the value) of
   items 3, 6, 7 and 10.

**Revised headline probabilities.** P(>=1 novelty-claimable cycler-class row at ANY tier in 12 weeks)
~70% (up from ~50%, almost all of it the two V1 Neptune-Triton orbits already in hand, whose novelty is
thin). P(>=1 at V3+) ~25-30% (down from 35-45%: spec-V4 unreachable this quarter, Uranian novelty
collapsed under the project's own standard, heliocentric novelty policy-gated).

**Revised ranking (final synthesis).** 1 integrity bundle (#865, expanded with the corrections above)
-> 2 Neptune-Triton rows with the verified recipe + row-first writeback of the 4:5 saddle, the
Oberon 4:5 family, the Saturn-Titan 3:4 and the #782 chain (#868 corrected) -> 3 Jones reproduction
to V3 (#866) -> 4 corpus gap + the NOVELTY POLICY decision (#869) -> 5 charter/lint/supervisor (#871,
now with REQUIRED `seed_diagnostics_run` and `environment_check` fields, the two failure modes this
review's own refutation found) -> 6 Jones enumerator and excluded classes, policy-gated (#867) -> 7
Triton-working / Proteus-passive Stage-0 re-gate (#874, seconds) -> 8 G1 corrector + symmetric pilot
(#872; the Neptune-Triton cell is done) -> 9 one-working-node campaign at Saturn and Uranus (#870
reframed, census unless policy GO) -> 10 shared node-locked inclined construction for Titan-Iapetus
and Triton-Proteus (flex) -> 11 make the V4 label honest, GMAT re-host first (#873a) -> 12 publication
package (#873b/c; wording "first two-moon Uranian cycler family", "first computed symmetric POs on a
Neptune-Triton homoclinic", "independent codebase, same SPK") -> 13 dynamical-object bets, idle-time.

**What is still missing after this review** (from the revision pass, agreed by the coordinator): the
ESM-family collision check for the two closed orbits; `ballistic_correct` has not been run at Jones's
dates with the multi-rev chain; Landau 2025 / Ellison 2025 / Campagnola 2014 / AIAA 2024-1280 /
Hollister 1969 / Minovitch 1967 unread; the Oberon 4:5 SOI entry checked for one member only; the
30-pair bend-gate table exists only as a summary; GMAT re-host and tudatpy installability on arm64
untested; the novelty policy unwritten; the 240 null-tier rows' semantics undetermined; no live
`literature_check.py` run on any proposed candidate; no measured CCR4BP pair-shot cost on the M3.

## 8. Do-not-do list (each item grounded in a closed task or a physics wall)

- Do not re-run the #600 Uranian 3-moon near-miss refinement (#663 did it; bend fails at ~0.83 deg).
- Do not re-run #695 Io-Europa verbatim under the #702 fix (#702's cross-check: 562-823 km off-torus
  independent of the bug); only an energy/family continuation is new.
- Do not re-run #110's epoch x 11-topology scan, #133's near-miss survey, or any VEM corrector scan
  that enumerates zero or one revolution per leg; do not loosen `VEM_VINF_TOL_KMS` or flip the xfail
  without a genuine converged pass; do not pose the heliocentric lane with the #388 Russell-psi
  generator + homotopy (the anchor reproduces from the published topology in 0.6 s).
- Do not run any more #563-class symmetric-closure sweeps anywhere (exhaustive at Uranus; asymmetric
  degenerate/empty; Saturn #571/#575/#655; Jupiter #577 0/36 literature-clear; small moons mass-
  limited, several interval-certified). The strategists' bend-gate table over all 30 moon pairs
  confirms no other pair passes.
- Do not build CCR4BP on Titania-Oberon (Kumar's own pair), on any Saturn pair with the inner moon as
  base (mu_pert 40-8000x below the JEG reference), on any Neptune pair, or on Miranda-Ariel.
- Do not run #790 as registered; do not re-attempt single-shooting or artificial-homotopy chain
  closure (#773/#775 are sharp negatives); do not append a lone unweighted Jacobi row; do not build
  the ML seed proposer before a zero-ML baseline exists; do not claim Jupiter-Europa 3:4<->5:6 or
  Saturn-Titan 3:4<->6:5 chains as novel (Anderson-Lo 2011; Vaquero 2013); Earth-Moon chains are
  controls only.
- Do not un-shelve the Resonant Atlas as a discovery campaign before #862 passes >=4/6 on Oberon and a
  writeback decision exists; do not fund #791 Stage 1, a Pluto-Charon resonant-manifold lane, or a
  Uranus N=5 CRNBP.
- Do not launch a low-thrust/SEP or powered "novel cycler" sweep (schema already admits powered bands;
  powered moon tours swept #464/#465; heliocentric powered Earth-Mars literature-owned; thrust
  replacing a missing bend is a Merrill-class forced periodic trajectory, not a cycler).
- Do not build another new method family aimed at W1-W3; do not re-sweep Earth-Mars, Earth-Moon or
  Sun-planet regimes for novelty.
- Do not run another analysis-only strategy pass until Gates A and B have run; do not dispatch any
  campaign without a lint-passing charter; do not mint a new orbit_class or schema version after the
  fact to admit a finding.
- Do not relabel any row V5 or set `verified-novel` on internal review; do not present the six
  Uranian rows as spec-V4 externally until item 9 lands or the label reads "V4-internal"; do not treat
  a Zenodo DOI as the publication; do not de-publish V0/null rows from cyclers.space (amend the spec
  to the labelled-census policy the site already follows); do not stamp the 240 null rows V0 unless
  the V0 check is shown to have run.
- Do not call any torus, homoclinic, heteroclinic, resonant PO or inter-moon transfer a cycler in the
  README, OUTSTANDING or a preprint; do not put a chain orbit into `manifold_connections.yaml`; do not
  give a subagent >10 minutes of compute or let it self-background.

## 9. Operating rules the review recommends adopting

1. No campaign dispatch without a lint-passing charter (item 7 fields).
2. Stamp every negative into `data/empty_regions.jsonl` the day it is decided, with reopen_condition;
   cross-reference superseding tasks on old stamps.
3. Enumeration-first, not corrector-first, for encounter-bearing searches: enumerate ALL
   revolution/branch topologies before any corrector runs.
4. A bug fix or a topology-enumeration gap invalidates past negatives: record the invalidation and
   re-scan (VEM: #110/#120/#122/#133 void).
5. "It closed!" and "0/N" are both hypotheses until independently re-derived.
6. Distinguish class in every artefact: encounter-bearing cycler/quasi_cycler vs dynamical object;
   CR3BP periodic orbits follow the #811/#855 per-member periselene-vs-SOI rule.
7. Owner launches any run >10 minutes under the supervisor, detached, n_workers=4 while CI shares the
   Mac; subagents build/validate with <=10 minutes of compute.
8. Every cost estimate cites a measured unit cost or says GUESS; measure the first cell of any new lane
   before dispatching the second.
9. Miss a gate by >1 week and the following wave is not dispatched; slack goes to consolidation.
10. Read the corpus before claiming novelty for anything Uranian or Neptune-Triton; keep the
    fortnightly literature watch dated; anchors go into `literature_check.py` code, not just notes.

## 10. Open decisions for the owner

- **NOVELTY POLICY (#875; the load-bearing one, take before Gate B)**: (i) is a member of a class the
  source authors EXCLUDED on stated practicality grounds (Jones 2017: 3-/5-synodic VEM, >6 flybys)
  `candidate-novel` or known-class? (ii) is a known architecture (Russell-Strange one-working-node) at
  a system the authors never treated (Uranus, Neptune) `candidate-novel` or known-class, per #817's
  reading of #577? (iii) are theorem-generic accumulation orbits on a homoclinic (the two Neptune-
  Triton orbits) `candidate-novel`, and what `our_status` do they and the two torus rows take? Write
  the answer into spec sec 16.4/16.5 as a one-paragraph rule. Every novelty probability in roadmap
  items 6, 7, 9 and 10 is conditional on it; their census/V-tier value is not.
- **GMAT re-host route (#873a)**: the installed R2022a is a Linux x86-64 build; choose Rosetta 2 (not
  installed) or an x86_64 container (colima/docker), 0.5-1 day either way, before any spec-V4 claim.
- **Mission scope**: do novel dynamical objects count, or only encounter-bearing cyclers? Decides
  item 13's owner-time and the README's language for the two torus rows.
- **Powered-cycler admissibility rule** (spec amendment): Russell-Ocampo <300 m/s per 7 cycles =
  cycler; powered_dsm to 3.5 km/s/cycle = powered cycler; thrust-shaped forced periodic trajectories
  out of scope. Resolves the Merrill-digest vs schema v4.8 contradiction; needed for item 10 path B.
- **Spec sec 14 one-sentence amendment** permitting a candidate-framed preprint (today "submitted for
  publication" is reserved for V5, so V5 can never be reached).
- **Writeback route for zero-row systems**: row-first (recommended; item 4b executes it) vs relaxing
  the connection registry's admission rule.
- **our_status vocabulary** for the two torus rows: `candidate-novel` vs a class-explicit value.
- **License split** for the Zenodo release (MIT code / CC-BY-4.0 data).
- **Null validation_level semantics** on 240 rows.
- **Portfolio flex**: fund item 10, the G1 asymmetric tail, item 12 idle time?
- **Heliocentric allocation**: accept the re-cut from the program-manager strategist's 10% cap to a
  primary lane (~25% of owner time) on the strength of the verified probe (all three judges
  recommend it).
- **Chain V-tier above V1** (deferred): accept a project-documented station-keeping budget as spec
  sec 14's "documented operational dV" for unstable CR3BP rows?

## 11. Process notes and uncertainties

- Agent counts: Phase 1 six agents (1.22M tokens); Phase 2 twenty-seven agents (2.94M tokens), of
  which 14 refuters and the revision pass failed on a session quota at first attempt and were
  re-run after the reset. One strategist proposal (the #600 refinement) was a fatal re-proposal of
  closed work and was caught by all three judges — the same failure class the registry-stamp gap
  explains.
- The probe script lives only in the session scratchpad; item 2 should commit a cleaned version as
  `scripts/probe_jones_multirev_truth.py` under `tests/scripts` preflight.
- INFERRED items that matter: absence of post-2017 VEM sweeps (conference proceedings under-indexed);
  UOP-era Uranian papers unread; the Neptune-Triton chain closure is a prediction; n-body survival of
  a 19.2-yr VEM chain is untested; per-pair CCR4BP cost on the M3 is unmeasured.
- Heliocentric multi-planet space beyond Earth-Mars (15 Hollister-Menning Earth-Venus rows, 4 Jones
  rows, all validation null) has never been swept with the mature toolkit; #287 (VEM heliocentric
  re-scan) was registered in June and never built. Item 3 replaces it.
- No file other than this note and `data/OUTSTANDING.md` was modified by this task.
