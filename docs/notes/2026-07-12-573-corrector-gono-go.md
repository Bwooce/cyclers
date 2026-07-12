# #573 corrector go/no-go — does REAL_FAMILY_SIGNAL justify scoping the narrow #552 Titan–Iapetus 3D corrector?

**Date:** 2026-07-12
**Task:** the final go/no-go on the chain #552 → #571 → #572 → #573. Decide whether to
scope+greenlight the narrow Titan–Iapetus 3D corrector build (#552's own re-scoped ~3–5 day
estimate, reusing `core/lambert.py`'s confirmed 3D support — NOT the general
n-body/arbitrary-inclination capability that was killed on 2026-07-10).
**Scope constraint:** read/analysis/writing only. Do NOT build the corrector, do NOT edit
`catalogue.yaml`, do NOT fetch any SPICE kernel, do NOT claim novelty (quasi-cycler-class,
internally-enumerated-search-space framing throughout). Analysis reuses
`data/probe_573_titan_iapetus_population_closure.jsonl` (read directly, all 22 branch records),
`scripts/run_573_titan_iapetus_population_closure.py`, `verify/mission_spk.py`,
`core/satellites.py`, and the full #552/#571/#572/#573 chain in `data/OUTSTANDING.md`.

---

## 1. Decision in one line

**GO — scope and greenlight the narrow Titan–Iapetus 3D corrector as new task #574, but STAGED:
its first milestone is a cheap eccentric-Keplerian closure kill-gate (no SPICE) that must pass
before the expensive real-ephemeris SPICE SAT441 validation pipeline is built.** This honors the
pre-registered #573 decision contract (which #573 cleared with margin — not scoping now would be
goalpost-shifting against a Fable-reviewed gate) while preserving the cheap-check-first discipline
by making the build's own Stage 1 the cheapest possible attack on the one load-bearing unresolved
risk (eccentricity). Require a Fable second-opinion on the #574 spec before dispatch, per this
thread's standing discipline.

## 2. The #573 data, read directly — the 22 branches (not just the aggregate counts)

Reading every branch record in `data/probe_573_titan_iapetus_population_closure.jsonl`:

| id | nmem | n_rev | Ω° | tof | resid km/s | V∞ Titan | V∞ Iap | bend_Iap° | ecc-rob |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 1 | (2,2) | 180.9 | 2.80 | 5.5e-10 | 3.425 | 1.547 | 6.53 | ✓ |
| 1 | 3 | (1,1) | 357.9 | 1.65 | 7.0e-10 | 3.478 | 1.378 | 8.11 | ✓ |
| 2 | 1 | (0,0) | 186.8 | 0.56 | 7.8e-10 | 4.025 | 1.380 | 8.09 | ✓ |
| 3 | 2 | (0,0) | 351.7 | 0.89 | 7.9e-10 | 1.787 | 1.590 | 6.20 | ✓ |
| 4 | 1 | (0,0) | 186.8 | 0.56 | 8.0e-10 | 4.331 | 1.574 | 6.32 | ✓ |
| 5 | 1 | (0,0) | 333.8 | 0.51 | 9.5e-10 | 3.962 | 1.546 | 6.53 | ✓ |
| 6 | 3 | (0,0) | 37.4 | 0.69 | 1.0e-09 | 2.341 | 1.205 | 10.38 | ✓ |
| 7 | 2 | (0,0) | 345.0 | 0.69 | 1.1e-09 | 2.025 | 1.193 | 10.57 | ✓ |
| 8 | 1 | (0,0) | 51.1 | 0.50 | 1.1e-09 | 4.169 | 1.617 | 6.00 | ✓ |
| 9 | 2 | (1,1) | 3.7 | 1.91 | 1.2e-09 | 3.562 | 1.420 | 7.66 | ✓ |
| 10 | 1 | (2,2) | 0.9 | 2.80 | 1.5e-09 | 3.772 | 1.335 | 8.60 | ✓ |
| 11 | 3 | (1,1) | 72.9 | 2.27 | 1.5e-09 | 2.642 | 1.730 | 5.28 | ✗ |
| 12 | 1 | (1,1) | 347.3 | 1.68 | 1.6e-09 | 3.835 | 1.341 | 8.53 | ✓ |
| 13 | 2 | (0,0) | 321.8 | 0.75 | 1.6e-09 | 3.633 | 1.470 | 7.18 | ✓ |
| 14 | 1 | (0,0) | 6.8 | 0.56 | 2.1e-09 | 4.331 | 1.574 | 6.32 | ✓ |
| 15 | 1 | (1,1) | 45.6 | 2.52 | 3.2e-09 | 1.473 | 1.714 | 5.37 | ✗ |
| 16 | 2 | (2,2) | 345.9 | 2.58 | 4.0e-09 | 4.220 | 1.482 | 7.08 | ✓ |
| 17 | 1 | (1,1) | 358.8 | 2.00 | 4.1e-09 | 2.254 | 1.775 | 5.02 | ✗ |
| 18 | 1 | (0,0) | 5.4 | 1.01 | 2.1e-06 | 3.296 | 1.722 | 5.32 | ✗ |
| 19 | 1 | (0,0) | 6.8 | 0.56 | 5.0e-06 | 4.025 | 1.380 | 8.09 | ✓ |
| 20 | 1 | (0,0) | 350.1 | 0.78 | 2.4e-05 | 3.520 | 1.450 | 7.37 | ✓ |
| 21 | 1 | (1,1) | 165.6 | 1.68 | 2.9e-03 | 3.344 | 1.775 | 5.02 | ✗ |

**Aggregate feel:** Iapetus V∞ 1.19–1.78 km/s (median 1.51) — a tight, physically-coherent band,
exactly the `sqrt(V∞_coplanar² + 0.87²)` quadrature-corrected range the #552 revival assessment
predicted analytically. Iapetus bend 5.02–10.57° (median 6.81°), Titan V∞ 1.47–4.33, tof_scale
0.50–2.80 spanning all three n_rev classes ((0,0)×12, (1,1)×7, (2,2)×3). Branches 18–21 have
notably looser residuals (2e-6 to 3e-3 km/s) but all sit far inside the 0.05 km/s gate; the other
18 are at machine precision (~1e-9). Branches {2,19} and {4,14} are the z-mirror-reflection pairs
the #573 note documents (identical V∞/bend, Ω differing ~180°) — kept distinct per the #563
symmetric-partner precedent. Even collapsing both mirror pairs leaves **20 branches**, still ~4×
the ≥5 threshold.

## 3. Is REAL_FAMILY_SIGNAL strong enough, or is it a smooth-problem artifact (the #558 concern)?

The sharpest version of the worry: #558's *original uncorrected* Uranian sweep found "1000+
basins" before the real family turned out to be a precisely-enumerated 30 — could Titan–Iapetus's
22 be the same grid-aliasing inflation on a smooth, generically-solvable problem? **My answer: no,
and the reasons are concrete.**

**(a) This is NOT the #558 situation — the anti-aliasing discipline that #558 initially lacked is
already built into #573.** #558's "1000+" was a *raw pre-dedup* grid-basin count. #573's raw count
is only **33 closing-basin instances** (from 28/69 candidates), deduped via Fable-mandated
union-find proximity clustering (same n_rev ∧ Ω* within 5° ∧ tof_scale* within 0.05 ∧ V∞ within
0.10 km/s) to **22**. The exact failure mode that produced #558's 1000+ (grid samples along a
continuous branch double-counted; the ±0.1 refinement window wider than the 0.05 grid spacing) was
identified *by name* in the #573 Fable plan review and explicitly fixed and verified (the 18-point
rel156–173 run collapsed to zero duplicate closures; eight adjacent-run groups each merged to one
branch). The raw count here (33) is not remotely "1000+"; the machinery already applied the
correction #558 had to learn.

**(b) The method here ERRS toward UNDERcounting, opposite to inflation.** The single most
reassuring fact in the whole run: #572's *highest-margin* confirmed closure (cand1, rel254–257,
tof 1.8, n_rev (1,1), a genuine residual-1.7e-9 closure under #572's own fixed ±15° window)
surfaces in #573 as a **formulation-conditional non-closure** (boundary-pinned) because its true 2D
minimum (Ω=30.91°) sits in a 33.6° gap between two raw-grid basins that the adaptive window
correctly refused to branch-jump across. A *known-genuine* closure is therefore **missing** from
the 22 — so 22 is a conservative floor, not padded aliasing. An artifact-inflation failure would
show the opposite signature (spurious extra closures), which is absent: 0 geometry errors, 0
convergence errors across ~497k Lambert solves, 0 isolated singleton-flip anomalies.

**(c) The decision is robust to residual count uncertainty.** One legitimate difference from the
Uranian precedent tempers a naive "22 ≫ ~5/pair, so it's 4× richer" read: the #573 probe carries an
**extra free dimension** (the node-alignment angle Ω) that the coplanar Uranian symmetric-closure
enumeration (#563) did not — more knobs can surface more basins, so part of the 22-vs-~5/pair gap
is the added DOF, not intrinsically richer dynamics. **But the go/no-go does not depend on the
exact number.** Whether the true independent-branch count is ~14 or ~22, it clears the
pre-registered ≥5 / ≥2-classes / ≥3-robust threshold by a wide margin on all three axes (22/3/17),
spans every n_rev class, and — per (b) — is a floor. The signal survives any defensible tightening
of the dedup.

**Verdict on the artifact question:** REAL_FAMILY_SIGNAL is a genuine idealized-circular family
signal, comparable in standing to the Uranian family that *did* earn a catalogue writeback
(#569), not a smooth-problem counting artifact. What it is **not** is a real-ephemeris result —
see §5.

## 4. What the corrector actually needs to prove, and how much already exists

`scripts/run_573_titan_iapetus_population_closure.py` (reusing `probe_572`'s `iapetus_state_3d`,
`sweep_node_alignment`, `_leg_best`, the R3(Ω)·R1(inc) rotation algebra, and the inc=0 smoke test)
is **already most of an idealized-circular 3D closure engine** — it takes coplanar seeds, places
Iapetus on a real inclined orbit, searches node alignment + tof_scale, enumerates/dedups basins,
and gates each against the #324 physical gate. So the "corrector" is emphatically **not** a
from-scratch build. What is genuinely still missing for a real, repeatable, V-gauntlet-ready tool:

1. **Eccentricity.** The entire #571→#573 stack is *circular* — Titan z=0 coplanar-circular,
   Iapetus inclined-**circular**. `core/satellites.py` carries GM/R/a but **no eccentricity field
   at all** (verified). "eccentricity_robust" in #573 is only a bend-margin proxy (Iapetus bend
   ≥6.0°, i.e. surviving a worst-case ±0.16 km/s shift *at the bend gate*), **not** an actual
   eccentric-orbit closure test. This is the load-bearing gap.
2. **Productization behind the #324 gate.** The #573 script is a throwaway sweep, not a repeatable
   corrector with a clean interface, checkpointing/provenance, and a stable input contract.
3. **A real-ephemeris (SPICE) leg.** `verify/mission_spk.py` maps `"Saturn" → "SATURN
   BARYCENTER"` (verified) — the *barycenter*, not the moon bodies. A real Titan/Iapetus validation
   needs the **SAT441** satellite ephemeris (or equivalent). **No `.bsp` of any kind is on this
   host** (verified: only `naif0012.tls` LSK in `verify/kernels/`). SAT441 is the same fetchable
   NAIF dependency #506/#550 already identified for PC(3,2) V2→V3.
4. **The V1–V4-strict gauntlet** with the #560 robustness fixes and #568 duty-cycle/synodic-boundary
   framing (not raw pass%/flip%), honoring the eccentricity caveat throughout.

The original ~3–5 day estimate still looks right for items 2–4 *as a Saturn analogue of the Uranian
V1–V4 chain*, precisely because item 1's engine already exists. But item 1 (eccentricity) is both
the cheapest to attack AND the most likely to kill the whole thing — which is why #574 must **stage
it first** (§6).

## 5. Honest real-ephemeris expectation — the eccentricity picture

The Uranian family (#568) achieved **61.7–79.1% real-ephemeris duty cycles at e ≤ 0.004**.
Titan e ≈ 0.0288 and Iapetus e ≈ 0.028 — **7–25× larger**. Per the #552 revival assessment this is
~±0.16 km/s velocity modulation at each encounter, ~3× the 0.05 km/s residual floor and comparable
to the entire ~1.2 km/s Iapetus-encounter tolerance budget. The 17/22 bend ≥6.0° "ecc-robust"
sub-count shows the *flyby-bend feasibility* survives that shift — but that is necessary,
**not sufficient**, for *closure* survival: eccentricity perturbs the whole transfer geometry and
the moons' relative phasing, not just the bend margin.

**My honest expectation:** the real-ephemeris duty-cycle picture will be **materially worse than
the Uranian 61.7–79.1%** — I would not be surprised by ~30–55% for the surviving branches, with a
real (non-negligible) probability that broad swaths fail a V4-strict-equivalent gate entirely and
the family thins to a handful of eccentricity-tolerant members. The two ~1.19–1.21 km/s
low-V∞/high-bend (10.4–10.6°) branches (ids 6, 7) are the most likely survivors; the five bend-≤5.4°
floor-huggers (ids 11, 15, 17, 18, 21) are the most likely to die.

**Does that change the go decision? No — "worse duty cycle, still real" is an acceptable and
worth-finding-out outcome, and a clean documented negative is equally acceptable.** The Uranian
family earned a writeback at 61.7–79.1%; a Titan–Iapetus family that lands lower but non-zero is a
legitimate quasi-cycler-class characterization (lower duty cycle is a *column*, not a
disqualifier), and a family that dies on eccentricity is a clean, publishable-internally negative
on a novel search — exactly what `[[project_negative_results_registry]]` and
`[[feedback_never_give_up_reproducing_papers]]` (this is a *novel* search, so a clean negative is a
fine stop) intend. The cheap circular probes are now **exhausted** — they are all circular by
construction and cannot answer the eccentricity question. The only way to learn the answer is to
build the eccentric stage. What we must NOT do is pay for the SAT441 fetch + full SPICE V-gauntlet
*before* the cheapest eccentric check has had its chance to kill it — hence §6's staging.

## 6. The #574 scope — staged, eccentricity-first

I allocate **#574** (see the OUTSTANDING.md TASK ALLOCATIONS ledger update). Full spec is in the
#574 OUTSTANDING.md entry; the shape:

- **Stage A — eccentric-Keplerian closure kill-gate (cheap, ~part-day, NO SPICE, NO kernel
  fetch).** Extend the #573 circular generators to place Titan AND Iapetus on real *eccentric*
  Keplerian orbits (add e ≈ 0.0288 / 0.028; solve Kepler's equation for true anomaly; sweep the
  moons' periapsis phase at the encounter as the new free parameter), and re-run the #573 3D
  closure + #324 gate on the **17 eccentricity-robust branches** (the non-robust five are expected
  to die and are the control). **Kill gate:** if the eccentric velocity modulation cannot be
  absorbed by the (Ω, tof_scale, phase) knobs for a substantial subset — i.e. closures collapse or
  the surviving population is only 2–3 branches — **stop here** with a documented conditional
  negative; do NOT fetch SAT441 or build the SPICE pipeline. This is the cheapest possible attack
  on the single load-bearing risk, decoupled from the SPICE cost.
- **Stage B — productized corrector + real-ephemeris SPICE validation (only if Stage A passes).**
  Productize the #573 engine behind the #324 gate; add the SAT441-backed real-ephemeris leg
  (`verify/mission_spk.py` Saturn mapping + the kernel fetch *scoped as part of #574*, not done
  now); run the V1–V4-strict gauntlet with #560 fixes and #568 duty-cycle framing; honor the
  eccentricity caveat throughout; behind the standard Opus adjudication + Fable second-opinion +
  V-gauntlet before ANY `catalogue.yaml` writeback.

**Framing (mandatory, unchanged from the whole chain):** any positive result is
**quasi-cycler-class evidence, same standing as #312's own family** (V2 fails on drift by design,
`FAIL_QUASI_BOUNDED`) — **not** a ballistic-cycler finding and **not** a novelty claim (an internal
fact about our own idealized enumerated search space, per `our_status` discipline).

## 7. Why not the two narrower alternatives

- **"Hold off / run yet another cheap probe first."** Rejected as goalpost-shifting. #573's
  decision gate was **pre-registered and Fable-reviewed** specifically to be the go signal, with a
  pre-registered eccentricity proxy (bend ≥6.0°) folded in; #573 cleared both (22/3/17). Adding a
  *new* pre-corrector adjudication cycle after a pre-registered gate passes is the exact
  discipline-erosion the project warns against. The right home for the remaining eccentricity check
  is **inside** #574 as its Stage-A gate, not as a fifth standalone probe+adjudication.
- **"Corrector for only the ecc-robust subset."** Effectively adopted — Stage A runs on exactly the
  17 ecc-robust branches, with the 5 floor-huggers as the expected-to-die control. Restricting the
  *whole build* to that subset up front would discard the control and the two mirror-pair members
  prematurely; staging captures the same cost saving without pre-judging.

## 8. Uncertainty flags for the Fable second-opinion pass on #574

1. **Is greenlighting the build the right call vs. one more cheap circular probe?** My §7 argument
   is that #573 cleared a pre-registered gate and the only remaining question (eccentricity) is
   *not* answerable by any circular probe — so the next increment MUST introduce eccentricity, and
   the cheapest place for it is Stage A of #574. If Fable thinks the eccentric-Keplerian check
   deserves to be its own standalone task (#574 = probe only, corrector = #575), that's a
   defensible re-slice — the substance (eccentric check before SPICE) is identical either way; I
   folded it into one staged task to avoid a redundant adjudication cycle.
2. **Eccentric-Keplerian ≠ real ephemeris.** Stage A uses eccentric *Keplerian* orbits (two-body,
   analytic, no third-body/SPICE). A Stage-A *pass* does not prove real-ephemeris survival (Stage B
   still must), but a Stage-A *fail* is a genuine cheap kill. It is a one-sided filter — correct and
   intended, but Fable should confirm the framing isn't oversold.
3. **Duty-cycle estimate (~30–55%).** This is a judgment extrapolation from the e-ratio (7–25×) and
   the ±0.16 km/s vs 0.05 km/s floor, NOT a computed number. It is directional, not a prediction to
   hold the build to; the honest statement is "materially worse than 61.7–79.1%, plausibly still
   non-zero, worth finding out."
4. **SAT441 sufficiency / kernel identity.** I cite SAT441 as the fetchable Saturn-satellite
   ephemeris (per the #572 adjudication and #550's PC(3,2) precedent) but did not verify the exact
   kernel name/coverage covers both Titan and Iapetus over the needed span — Stage B must confirm
   the precise kernel before fetching; that verification is scoped into #574, not done here.
5. **Should Phase-2 Iapetus-anchored (118 candidates) fold in?** #573 did not mandate it
   (REAL_FAMILY_SIGNAL, not the marginal bucket). I did not add it to #574 — the Titan-anchored 22
   already clear every threshold and adding the harder double-node phasing is scope creep for the
   corrector. Defensible either way; Fable may reargue.
