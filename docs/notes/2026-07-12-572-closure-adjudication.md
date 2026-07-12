# #572 closure adjudication — do the Titan–Iapetus 3D closures justify scoping the narrow #552 corrector?

**Date:** 2026-07-12
**Task:** adjudicate the #572 probe result (genuine 3D closures for both tested
Titan–Iapetus candidates) and decide the next step: scope+build the narrow
Titan–Iapetus 3D corrector (#552's re-scoped ~3–5 day estimate), widen the probe
first (cheap intermediate), or don't proceed.
**Scope constraint:** read/analysis/writing only. No corrector built, no new
sweep run, no `catalogue.yaml` / `empty_regions.jsonl` edit. Analysis reuses the
existing `data/probe_572_titan_iapetus_3d_closure.jsonl`, `data/scan_571_saturn_titan_iapetus.jsonl`,
and the #572 probe source.

---

## 1. Recommendation in one line

**WIDEN THE PROBE FIRST (new task #573), then scope the corrector only if the
widened probe shows a substantial gate-passing closing population.** Do NOT jump
straight to the 3–5 day corrector build off a 2-candidate sample. Do NOT stamp the
pair empty — the closures are real. The widened probe is nearly free (~1 s/candidate)
and answers the one question 2-for-2 on the best-selected candidates cannot: what
fraction of the *population* actually admits a gate-passing 3D closure.

## 2. What the probe data actually shows (read directly, both records)

Both tested candidates found a genuine 3D closure that clears both gates:

| | cand1 (rel255/tof1.80/n11) | cand2 (rel89/tof0.70/n00) |
|---|---|---|
| closing Ω (node) | 30.91° | 37.46° |
| refined tof_scale | 1.793 | 0.688 |
| residual | 1.69e-9 km/s | 2.01e-9 km/s |
| V∞ [Titan, Iapetus, Titan] | [1.870, 1.207, 1.870] | [2.307, 1.203, 2.307] |
| bends [Titan, **Iapetus**, Titan] | [45.47, **10.35**, 45.47] | [34.04, **10.41**, 34.04] |
| n closing basins / n basins | 2 / 8 | 1 / 4 |

The closures are real converged conic arcs — the probe's own DOP853 cross-check
put arrival error at 1.1e-5 to 1.2e-4 km, far under the 1 km floor. I re-read the
full basin enumeration for both candidates; the verdict records are internally
consistent with the basin lists, and the low-V∞ gate-passing basin is genuinely a
*different* basin from the high-V∞ (~8.6–8.7 km/s) global-residual basin that
fails the bend gate. The all-basin dump confirms the high-V∞ basins (Ω≈176°)
sit at ~3.25° bend — below the 5° floor — exactly as the probe's design notes warned.

**One factual correction to the hand-off summary.** The dispatch framing cited
"Iapetus bend 45.5° and 34.0°." Those are the **Titan** bends. The **Iapetus**
(flyby-body, i.e. the *binding*) bends are **10.35° and 10.41°** — clearing the 5°
floor by ~2.07×, not the ~7–9× that "45.5°/34.0°" would imply. Still a clean pass,
and it matches the #571 assessment's own analytic prediction (coplanar Iapetus V∞
≈0.98 → real-corrected bend ≈8.9°; the realized 3D V∞ rose to ~1.2 km/s so the bend
landed a little higher at ~10.4°). But the real Iapetus margin is 2×, not 9× — this
is the correct number to carry forward, and it is the flyby-body margin that governs
whether *floor-hugging* candidates survive.

## 3. Why 2-for-2 is genuinely non-representative (the load-bearing caveat)

Two independent reasons this sample cannot be read as a population closure rate:

**(a) Selection bias — these are the best-margin candidates by construction.**
The #572 entry (Fable correction) chose exactly the V∞≈0.98 Titan-anchored cluster,
the *most robust* candidates in the population (coplanar corrected bend 8.8–8.9°).
The #572 entry's own note (4) is explicit: post-correction **median** bend across the
149 bend-gate survivors is **5.50° — right at the floor**; the robust core is *only*
this V∞≈0.98 cluster. So 2-for-2 on the two most robust candidates says almost
nothing about the median candidate, which sits at the floor and could flip either
side of it after the real 3D correction.

**(b) Machine-precision closure is generically achievable and is NOT the
discriminating fact.** The probe drives a scalar residual to ~0 using **two free
knobs** (Ω and tof_scale). A 2-parameter fit to a 1-D residual is expected to reach
~0 whenever *any* basin exists — the 1.7e-9 / 2.0e-9 residuals are not evidence of
rarity. The genuinely discriminating content is whether the zero-residual solution
lands in a **low-V∞, bend-gate-passing** branch rather than the high-V∞ (~8.6 km/s,
3.25° bend) branch that also closes to machine precision. For the robust core it
does, with 2× margin. For a floor-hugging candidate the closing basin's Iapetus bend
could land just under 5° — and we have zero data on how often that happens.

**Consequence:** a 3–5 day corrector justified on 2-for-2 risks discovering the
"family" is 2–3 members, not a 30-member #312-scale family. The cost of finding that
out the expensive way is days; the cost of finding it out the cheap way is minutes.

## 4. The widened probe is nearly free — this is the decisive cost argument

The #572 probe cost **1.00 s (cand1) and 0.25 s (cand2)** wall-clock — the full
3600-point Ω grid plus per-basin Nelder-Mead refinement, per candidate (measured,
in the JSONL `elapsed_s` fields). The relevant population is modest:
`data/scan_571_saturn_titan_iapetus.jsonl` carries **69 coplanar-gate-passing
Titan-anchored candidates** (all `physical_gate_passed=True` records are Titan-
anchored; n_rev split 26×(0,0), 40×(1,1), 3×(2,2)). Even the broader "187 flagged"
inclination-corrected set (the assessment's §2 population, which includes candidates
that fail the *coplanar* gate but pass after correction) runs in minutes at ~1 s each.

So the entire population-level closure rate can be measured **synchronously in a
single sub-10-minute run** — this is not a multi-day sweep. It is the same cheap-
necessary-check-before-capital discipline that produced #572 itself: #572 was the
cheap closure check gating the corrector; the widened probe is the cheap population
check gating the *decision that the corrector is worth 3–5 days*. Running it first
is strictly dominant: it costs minutes and it converts a selection-biased anecdote
into a hit-rate + family-size estimate.

Note on population reconciliation: the raw scan file shows 69 Titan-anchored coplanar-
gate-passers; the "187 flagged / 149 survive / 118 Iapetus-anchored" figures in the
#552 revival assessment are the *inclination-corrected* re-check across a broader
candidate set, not the coplanar-gate-pass count. #573 must define its input population
precisely from the scan file (see the task entry). The Iapetus-anchored candidates
face the harder double-node phasing constraint the #572 entry deliberately deferred —
#573 should start with the Titan-anchored set the existing probe machinery already
handles, and treat Iapetus-anchored as an optional second phase.

## 5. What is still genuinely unknown — the eccentricity gap

The probe used **circular** orbits (Titan coplanar-circular per #571 convention,
Iapetus inclined-circular). It therefore says **nothing** about whether these
closures survive real ephemeris. The pre-registered caveat is material:

- **Titan e ≈ 0.0288** (7–25× the Uranian moons' e≤0.004), **Iapetus e ≈ 0.028**.
  Per the #552 revival assessment §5.2 this is ~±0.16 km/s velocity modulation —
  ~3× the 0.05 km/s residual floor and comparable to the entire ~1.2 km/s Iapetus
  encounter V∞ budget's tolerance.
- The Uranian family (#568) achieved 61.7–79.1% real-ephemeris duty cycles at
  e≤0.004. With eccentricity an order of magnitude larger, the idealized→real
  duty-cycle gap for Titan–Iapetus should be **materially worse** — the realistic
  upside is capped, and it is entirely possible a perfect circular 3D closure still
  fails a V4-strict-equivalent real-ephemeris validation.

**This is genuinely unknowable until the corrector plus a real-ephemeris Saturn
validation pipeline exist.** Neither the current probe nor the widened probe resolves
it — both are circular. That is *acceptable* and even *desirable* sequencing: the
widened probe is the cheap filter that tells us whether the population is worth ever
building the real-ephemeris pipeline for. If the widened population hit-rate is low,
we stop before touching eccentricity at all. If it is high, eccentricity becomes the
central risk the corrector + validation stage must confront head-on (and honestly may
still kill it — that is a legitimate later outcome, not a reason to skip the cheap
step now).

## 6. Concrete shape of the eventual corrector deliverable (for scoping, NOT built here)

If #573 clears, the narrow #552 corrector's deliverable would be a **Saturn-system
real-ephemeris validation pipeline for the Titan–Iapetus 3D closures, mirroring the
Uranian V1–V4-strict chain (#566/#568)** — not the general n-body/arbitrary-inclination
capability. Concretely:

1. A narrow 3D corrector that takes a coplanar-seed Titan–Iapetus candidate, places
   Iapetus on its real inclined orbit, and refines to closure over (Ω, tof_scale) —
   essentially productizing the #572 throwaway probe behind the existing #324 gate.
2. A real-ephemeris leg (SPICE). **Feasibility:** `verify/mission_spk.py` already
   maps `"Saturn" → "SATURN BARYCENTER"`; the Saturn satellite ephemeris kernel
   **SAT441** is a known, small NAIF fetch from naif.jpl.nasa.gov (it is *not* on
   this host — same dependency #506/#550 already identified for PC(3,2) V2→V3). So the
   real-ephemeris deliverable is feasible but carries a one-time kernel-fetch dependency
   that #573's result should justify before incurring.
3. The V1–V4-strict gauntlet run with the **#560 robustness fixes** and the **#568
   duty-cycle / synodic-boundary framing** (not raw pass%/flip%), honoring the
   eccentricity caveat throughout, behind the standard Opus adjudication + Fable
   second-opinion + V-gauntlet before any writeback.
4. Framing discipline (per #572 note 3): any positive result is **quasi-cycler-class
   evidence** (V2 fails on drift by design, `FAIL_QUASI_BOUNDED`), the same standing
   #312's own family holds — **not** a ballistic-cycler finding, and **not** a novelty
   claim (internally-enumerated fact about our own search space, per `our_status`
   discipline).

I am **not** allocating a task number for the corrector here — it is gated on #573's
population result and should be scoped concretely only once we know the family size.

## 7. Decision gate for #573 → corrector

- **If #573 finds a substantial gate-passing closing population** (a real family,
  not 2–3 isolated points): scope the narrow corrector (§6) as a subsequent task,
  behind Opus adjudication + Fable second opinion, with the eccentricity/duty-cycle
  caveat front-and-center.
- **If #573 finds only the robust-core cluster closes** (a handful): re-adjudicate —
  a 2–3 member idealized quasi-cycler set may not justify a 3–5 day real-ephemeris
  build; consider a lighter characterization or a conditional stamp.
- **If #573 finds broad closure but at Iapetus bends hugging the 5° floor**: the
  eccentricity gap (§5) will likely dominate and the corrector is high-risk — flag
  loudly before committing.

## 8. Uncertainty flags for the Fable second-opinion pass on #573

1. **Population definition.** I read 69 coplanar-gate-passing Titan-anchored
   candidates in the scan file vs. the assessment's "187 flagged / 149 survive."
   These count different things (coplanar-gate-pass vs. inclination-corrected
   survivorship). #573 must pin its exact input set — my recommendation is the full
   Titan-anchored candidate set the probe machinery already handles, but Fable should
   confirm whether the inclination-corrected-but-coplanar-fail candidates also belong.
2. **Is 2-for-2 really non-representative, or am I over-discounting it?** My §3
   argument rests on (a) explicit best-margin selection and (b) 2-DOF machine-precision
   closure being generic. If either is wrong — e.g. if the closing basin's bend is
   somehow *not* generically near the coplanar-corrected prediction — the widen-first
   case weakens and one could argue for proceeding. This is the judgment call most
   worth a second opinion.
3. **Iapetus bend margin.** I corrected the hand-off's "45.5°/34.0°" to the true
   flyby-body bends 10.35°/10.41° (2× floor). Fable should sanity-check that the
   encounter ordering (Titan, Iapetus, Titan) → bend index 1 is Iapetus is correct
   against `candidate_passes_physical_gate`'s per-encounter output.
4. **Should #573 stay circular, or add an eccentric variant?** I deliberately scoped
   #573 as circular-inclined (apples-to-apples with #572) to keep it a minutes-scale
   run and defer eccentricity to the corrector/validation stage. Fable may argue a
   cheap eccentric spot-check on the robust core is worth folding in; I judged it
   scope creep for the population-rate question, but it is defensible either way.
