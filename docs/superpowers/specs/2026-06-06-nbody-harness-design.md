# N-body harness — DESIGN DRAFT (no code)

**Status:** design draft (no code), 2026-06-06. Task #136. Sole write target is
this file. Built **once** for three consumers, not three times.

**The three consumers (the harness is sized to all three at once):**

1. **SILVER validation rungs 2–3** (#131) — propagate the Forge's two closed
   patched-conic E-M-E-E SILVER candidates through full n-body, measure closure
   degradation and the correction ΔV needed to restore periodicity (the Jones
   <200 m/s analogue). Gate: small correction ⇒ robust; large ⇒ patched-conic
   artifact.
2. **The Jones shooter** (#133, the deferred Phase 3 of
   `2026-06-06-jones-family-corrector-variant-design.md`) — near-miss conic seed
   → multiple-shooting differential correction in n-body → ballistic VEM
   trajectory matched against the sourced Jones multisets.
3. **Aldrin V4** (#134, spec §14 V4, `docs/spec.md:398`) — independent codebase +
   ephemeris reproduces the Aldrin trajectory and **our computed 2.9138 km/s**
   per-synodic maintenance ΔV (and powered periodic solution) within tolerance.

**Substrate read:** `core/ephemeris.py` (time/frame convention + the trap),
`search/correct.py` (the conic seed/solution object — `BallisticClosureResult`),
`scripts/correct_s1l1_twoarc.py` (the prototype corrector), `verify/propagate.py`
(the existing Kepler multi-lap drift machinery the harness mirrors),
`verify/gauntlet.py` (the verdict tiers the rungs feed),
`verify/spice_kernels.py` (#129 — kernel access + the same-kernel discipline),
`docs/notes/2026-06-06-performance-profile.md` (ephemeris cost reality),
`docs/notes/2026-06-06-silver-candidates-russell-diligence.md` (the candidates).

---

## 0. The non-negotiable: the time/frame conversion trap (read FIRST)

Every consumer compares an n-body propagation against **our** layer (the conic
seed in `BallisticClosureResult`, the catalogue V∞, the 2.9138 km/s). A
half-minute time-scale or frame slip silently corrupts every comparison and would
manufacture a "closure degradation" that is pure bookkeeping error. The harness
MUST replicate `core/ephemeris.py`'s conventions exactly:

- **Epoch convention.** `t_sec = 0` is **2000-01-01T12:00:00 TDB** (J2000), and
  the astropy backend builds `Time(_J2000_EPOCH.timestamp() + t_sec,
  format="unix", scale="tdb")` (`ephemeris.py:262-263`). The `t_sec` axis is a
  **TDB seconds-since-J2000** axis, NOT UTC and NOT TT. The "~64.184 s with
  drift" trap is exactly the TDB↔TT/UTC offset: an n-body tool that ingests our
  `t_sec` as UTC (or as ephemeris-time-without-leap-seconds) places every planet
  ~minute-of-arc off and corrupts the V∞ comparison. SPICE-native tools (Tudat,
  GMAT, spiceypy) must convert our `t_sec` → ET via the LSK (#129's
  `ensure_leapseconds_kernel`) and treat our axis as TDB throughout; the
  TDB↔ET identity (within ~ms) is the safe path, TT/UTC is the trap.
- **Frame.** Our state is **heliocentric J2000-ecliptic** (`+x` to the J2000
  vernal equinox, `+z` the ecliptic north pole), obtained by rotating ICRS
  (equatorial) by the J2000 mean obliquity −23.4392911° about +x
  (`ephemeris.py:269-275`). Any n-body tool returns barycentric or ICRS-equatorial
  states by default. The harness MUST apply (a) Sun-subtraction to heliocentric
  and (b) the **same obliquity rotation** before any vector is compared to a
  `BallisticClosureResult` V∞ or a catalogue value. Mars/Venus inclination lives
  in this rotation; getting it wrong reads as a fake out-of-plane V∞ component.
- **Same-kernel discipline (#129).** The cross-check's whole value is an
  *independent reader* over the *identical data* (DE440). The in-pipeline rungs
  (consumer 1, 2) should drive rebound/Tudat from the **same DE440 BSP** astropy
  cached (`spice_kernels.astropy_de440_bsp_path()`), so a disagreement is a
  reader/frame/time bug, not a different ephemeris. (The fully-independent V4
  stamp is the *deliberate* exception — see §4.)

A **golden two-body reduction test** (§5) and an **Earth-self-propagation anchor**
(§5) exist precisely to catch a violation of any bullet above before it touches a
science number.

---

## 1. Tooling trade-offs and recommendation

Three realistic options; the recommendation is a **deliberate split**.

### Option A — `rebound` (pip) + ephemeris-driven planets

- **What:** REBOUND with the **IAS15** adaptive high-order integrator — a
  genuinely independent integrator codebase (not our `core/kepler.py`, not
  scipy). Planets supplied as **point masses whose positions come from DE440**
  (the `assist` companion, or a thin ephemeris-force callback reading the same
  BSP). The spacecraft is a massless test particle; planets are perturbers.
- **For:** pip-installable, light deps, IAS15 is a publication-grade integrator,
  energy/momentum diagnostics built in. Independent of our Lambert/Kepler stack
  so it is a real cross-check for the in-pipeline rungs. Forking the worker
  inherits the DE440 mmap for free (perf profile §3 — fork-share is ~0 ms).
- **Against:** "ephemeris-driven point masses" is a *restricted* n-body (planets
  on rails, no spacecraft back-reaction — correct here, since the spacecraft is
  massless, so this is a feature not a bug). REBOUND alone has no SPICE time/
  frame machinery — the §0 conversions are **our** responsibility in the glue.
  `assist` adds a dep but does the DE440 force eval + frames correctly; weigh
  `rebound`+`assist` vs hand-rolling the ephemeris-force callback.
- **Independence:** strong for an *integrator*; ephemeris is shared-DE440 by
  design (good for rungs, see §4).

### Option B — Tudat / `tudatpy` (or pykep)

- **What:** SPICE-native ephemerides, built-in propagators, variational
  equations (state-transition matrices) and a built-in multiple-shooting /
  differential-correction toolchain.
- **For:** SPICE time + frames are handled natively (the §0 trap is largely
  solved by the library), and it ships the **STM/variational machinery the Jones
  shooter (consumer 2) actually needs** — multiple-shooting wants ∂(final
  state)/∂(node state) Jacobians, which Tudat produces analytically rather than
  by finite difference. pykep gives Lambert + Taylor propagation in the same
  family.
- **Against:** heavier install (Tudat is conda-centric, large binary deps —
  awkward in our uv-managed venv and in CI). Some independence is lost for the
  rungs if Tudat also reads DE440 via SPICE (still an independent integrator +
  reader, which is what matters).
- **Independence:** independent integrator + independent reader (CSPICE), shared
  DE440 data — same independence class as REBOUND for the rungs, *plus* the STM.

### Option C — NASA GMAT (scripted batch)

- **What:** a fully independent flight-dynamics application; scripted `.script`
  batch runs, parsed output.
- **For:** the **formal V4 stamp** (spec §14 explicitly names GMAT). Maximally
  independent: separate application, separate ephemeris, separate everything —
  exactly what "independent codebase + ephemeris reproduces … within tol" wants
  for consumer 3.
- **Against:** awkward in CI (GUI-era app, scripted batch is brittle, slow,
  hard to pin/version in a uv project). Not a library — cannot be called inside
  `least_squares`/the shooter loop. Wrong tool for the rungs and the shooter;
  right tool for the stamp.
- **Independence:** maximal (the gold standard for V4).

### Recommendation — SPLIT (honest weigh-up)

**One in-pipeline Python library for consumers 1 + 2, GMAT scripts for the
consumer-3 V4 stamp only.**

- **In-pipeline library: REBOUND/IAS15 as the baseline, with Tudat as the
  shooter-grade upgrade if/when the Jones shooter needs analytic STMs.** Reasons:
  the rungs (consumer 1) need only forward propagation + a correction solve —
  REBOUND is the lightest genuinely-independent integrator and installs cleanly
  in our venv/CI. The shooter (consumer 2) needs Jacobians; if finite-difference
  STMs over REBOUND prove too slow/noisy (likely, given the ephemeris cost in the
  perf profile), promote that consumer to Tudat for analytic variational
  equations. Design the harness's propagator behind a **single interface** (§3)
  so REBOUND↔Tudat is a backend swap, not a rewrite — mirroring how
  `core/ephemeris.py` hides `circular`/`astropy`/`inclined` behind one `state()`.
- **GMAT: V4 stamp only.** A small `scripts/gmat_v4_*.script` + a parser, run
  **out of CI** (a manual/batched lane, like the perf-profile reproduction
  scripts), producing a one-line "GMAT reproduces 2.9138 km/s within tol" record
  for the Aldrin row. Do **not** force GMAT into the rung/shooter loops.

Why not a single tool: GMAT can't be called inside a corrector loop; REBOUND
alone lacks analytic STMs for the shooter; Tudat alone is a CI/venv burden for
the cheap rung consumer. The split puts each tool where its cost is justified.

---

## 2. Force model

**Restricted n-body, heliocentric, ephemeris-driven point-mass perturbers.**

- **Bodies:** Sun (central) + the perturbers each consumer's trajectory actually
  passes near, planets as DE440-driven point masses (positions on rails from the
  shared BSP; spacecraft massless ⇒ no back-reaction, so "on rails" is exact).
- **Per consumer:**
  - **Rungs / SILVER E-M-E-E (consumer 1):** Sun + **Earth + Mars**, mandatory;
    **Jupiter** included as a perturber. Justification: the candidates are
    heliocentric E-M chains with multi-rev Earth-Earth loop legs spanning ~4–12 yr
    (`silver-candidates` note); over multi-year heliocentric arcs Jupiter's
    secular tug is the largest third-body term after the encounter bodies. The
    rung *metric* is correction ΔV at the ~tens-to-hundreds-of-m/s level (the
    Jones <200 m/s analogue), so a body that contributes ≳ a few m/s over a leg
    matters — Jupiter does at multi-year baselines; Venus/inner planets do not
    for an E-M chain. **Include Earth, Mars, Jupiter; Sun central.**
  - **Jones VEM shooter (consumer 2):** Sun + **Venus + Earth + Mars**,
    mandatory (they are the flyby bodies). **Jupiter:** include — the Jones
    repeat period is 12.8 yr (jones design §0.1) and the gate tolerance is
    0.5 km/s with a sourced ballistic target Jones removed to <200 m/s with
    SNOPT; at that precision Jupiter's perturbation over a 12.8-yr arc is not
    obviously below 100 m/s and must be carried to match Jones's own n-body
    model. **Does Jupiter matter at the 100 m/s level?** Yes, plausibly, over
    12.8 yr — so include it; the timestep-convergence + body-inclusion
    sensitivity test (§5) quantifies it and lets us drop it only with evidence.
  - **Aldrin V4 (consumer 3):** match whatever GMAT/Tudat's default high-fidelity
    EM-cycler model uses (Sun + planets through Jupiter at least); the V4 point is
    an *independent* model, so use the external tool's standard planetary set
    rather than minimising — reproduction within tol is the goal, not a minimal
    model.
- **Explicitly OUT (justify):**
  - **SRP** — out. The catalogue and all targets are *ballistic/gravity-only*
    trajectories with no spacecraft area/mass model; SRP is mission-design
    hardware detail absent from every source value. Including it would compare
    our gravity-only seed against a gravity+SRP propagation and fabricate a
    spurious correction ΔV. (A consumer that later models a specific spacecraft
    can add it; not in this harness.)
  - **Relativity (GR / PPN)** — out. The heliocentric V∞ comparison is at the
    km/s level with a 0.1–0.5 km/s tolerance; the GR perihelion-precession
    perturbation on a multi-year heliocentric arc is orders of magnitude below
    that and below the DE440 ephemeris frame/time error budget itself. Carrying
    it adds dependency and risk for no measurable change. (DE440's *planet*
    positions already embed GR; we only omit GR on the *spacecraft* dynamics.)
  - **Non-spherical gravity / moons** — out; heliocentric point-mass perturbers
    only. Flyby bending is handled at the patch points by the conic/B-plane
    layer, not by integrating through a planet's sphere of influence (consistent
    with the patched-conic seed contract).

The standing rule: **include a body iff it can move the consumer's metric by more
than that consumer's tolerance over that consumer's baseline; prove inclusion/
exclusion with the §5 sensitivity test rather than asserting it.**

---

## 3. Shooting architecture

A **multiple-shooting differential corrector over encounter nodes**, sharing the
seed contract with `search/correct.py` so the conic solution flows straight in.

### Node / variable structure

- **Nodes** = the encounters (the patch points). This matches both the conic
  model (jones design §1.3: "the encounters *are* the patch points", legs are
  single Lambert arcs) and `correct.py`'s `b{i}_in/b{i}_out` node vocabulary.
- **Free variables (per node):** node heliocentric state (or, equivalently, the
  outgoing V∞ vector + epoch), and the inter-node times of flight. This
  generalises `correct.py`'s `x = [t0, *tof]` (with one slack leg pinned by the
  sourced period) to a full multiple-shooting vector `x = [{node states}, {node
  epochs}, {ToFs}]`. The slack-leg / period-pin convention carries over (a leg
  ToF reconstructed as `period − Σ(free legs)`, `correct.py:_reconstruct_tofs`).
- **Defects (continuity constraints):** propagate each leg forward in n-body from
  node i; the defect is `(propagated state at node i+1) − (node i+1 state)`. Drive
  defects → 0. This is the n-body analogue of `correct.py`'s magnitude-continuity
  residual — but here it is **full state continuity in real dynamics**, which is
  strictly stronger (it is what makes the trajectory *truly* ballistic, the thing
  Jones's SNOPT step delivers).
- **Flyby constraints at each node:** periapsis radius `r_p ≥ r_p_safe` (reuse
  `PLANETS[body].safe_alt_km`, `correct.py:_max_bend_deg`); the bend the node
  imposes must be a real gravity-assist turn (the B-plane-targeted powered flyby
  of the jones design §1.1, Eqs.1–2). The **vector-residual / B-plane lessons
  from the jones corrector design carry in directly** — feasibility lives *inside*
  the solve, not post-hoc.
- **Periodicity:** the wrap node (home Earth) must map to the start node modulo
  the repeat period — the n-body version of `correct.py`'s closure term
  (`b{last}_in` ↔ `b0_out`).

### How the conic seed maps in

`BallisticClosureResult` carries `t0_sec`, `tof_days`, `vinf_per_encounter_kms`,
and (via `_vinf_nodes`) the per-encounter V∞ vectors. The harness seeds:

- node epochs from `t0_sec + cumulative(tof_days * 86400)` — **on the TDB
  J2000 axis** (§0), handed to the n-body tool with the time conversion applied
  once at the boundary;
- node states from `ephem.state(body, epoch)` (planet) + the seeded outgoing V∞
  (spacecraft) — exactly the `v_sc = v_planet + vinf` reconstruction
  `verify/propagate.py:434` and `correct_s1l1_twoarc.py` already do;
- ToFs and the slack/period pin directly from the result.

So the conic corrector's output object **is** the shooter's initial guess; no new
seed format is invented (golden-discipline-safe: the seed is sourced/computed
upstream, the shooter only refines it).

### Convergence / divergence handling

- **Solver:** Newton / Levenberg–Marquardt on the defect+constraint residual
  (mirror `correct.py`'s `least_squares(method="lm")` choice for continuity with
  the conic layer). Jacobians: finite-difference over REBOUND (baseline) or
  analytic **state-transition matrices** if on Tudat (the reason consumer 2 may
  need Tudat — finite-diff STMs over an expensive ephemeris propagation are slow
  and noisy; the perf profile shows propagation, not the optimiser, is the cost).
- **Divergence is a first-class outcome, surfaced not crashed.** Follow
  `correct.py:359-379` and `verify/propagate.py:733` (`KeplerConvergenceError` →
  an honest non-converged report). A diverged shoot for a rung candidate is
  itself the verdict signal (large/unbounded correction ⇒ patched-conic
  artifact, feeds REJECTED-style evidence), not an exception to swallow.
- **Multiple-shooting over single-shooting** precisely because long multi-year
  multi-rev legs (12.8-yr Jones, multi-rev E-E loops) make single-shooting
  ill-conditioned; node-splitting at the encounters keeps each propagation arc
  short enough to stay in the linear regime for the Newton step.

### The "correction ΔV" convention (the SILVER rung metric)

This is the headline number for consumer 1 and must be defined unambiguously:

- **Correction ΔV** = the sum of the impulsive ΔV that, applied **at the encounter
  nodes** (or as TCMs along the legs — pick one and state it), restores full
  n-body periodicity from the conic seed. Concretely: after the shooter converges
  the defects to zero, the correction ΔV is the **change in the spacecraft
  velocity discontinuity at each node** between the raw conic seed (which is
  *not* n-body-continuous) and the corrected ballistic solution — i.e. the ΔV the
  flybys could *not* absorb gravitationally and that a real mission would burn.
- This is the **direct Jones <200 m/s analogue** (jones design §2 (d), spec §14
  V4 "maintenance ΔV"): small correction ⇒ the conic seed was a faithful shadow
  of a real ballistic trajectory ⇒ robust SILVER candidate; large correction ⇒
  the seed lives only in patched-conic land ⇒ artifact.
- **Accounting must be sign- and node-explicit and reconcilable** with the
  existing `maintain.py` per-synodic maintenance ΔV (the 2.9138 km/s the V4
  consumer reproduces) — the rung ΔV and the maintenance ΔV are the *same kind of
  quantity* (ΔV to keep the cycle periodic) and should use one convention across
  consumers 1 and 3 so the gauntlet can compare them.

---

## 4. Independence accounting (§14 V4)

Spec §14 V4 demands **"independent codebase + ephemeris"**. Catalogue the exact
combinations and what each consumer is entitled to:

| Combination | Codebase indep? | Ephemeris indep? | Satisfies V4? | Used by |
|---|---|---|---|---|
| **GMAT** scripted, GMAT's own ephemeris | yes (separate app) | yes (separate eph) | **YES — the V4 stamp** | consumer 3 (Aldrin) |
| **Tudat/CSPICE** + a *separately-downloaded* DE440 | yes (separate integrator + CSPICE reader) | yes (own kernel copy) | **YES (weaker than GMAT but qualifies)** | consumer 3 alt, consumer 2 |
| **REBOUND/IAS15** + **shared astropy DE440 BSP** | yes (independent integrator) | **NO — same BSP** | **NO** (same data) — but valid as an in-pipeline *cross-check rung*, not the V4 stamp | consumers 1 + 2 |
| our `core/kepler.py` + astropy backend | no (our code) | no (our layer) | no | (the thing being checked) |

Rules:

- **The V4 stamp (consumer 3) requires independent ephemeris too** — so it uses
  GMAT (gold) or Tudat with its *own* kernel copy. REBOUND-on-shared-DE440 does
  **not** earn the V4 badge, by design.
- **The in-pipeline rungs (consumers 1, 2) deliberately SHARE DE440** with our
  astropy backend (#129 same-kernel discipline): there, an independent
  *integrator* over *identical data* is exactly the right tool — a disagreement
  is a reader/frame/time bug, which is what the rung wants to catch. They are
  cross-check rungs (Axis-A-like, feeding the gauntlet), **not** V4 stamps, and
  the design must label them as such so a SILVER candidate is never mislabelled
  V4-passed on a shared-ephemeris run.
- **What the rungs may share with our layer:** the DE440 BSP, the J2000/TDB axis,
  the obliquity rotation (so the comparison is apples-to-apples). What they may
  **not** share: the integrator/propagator (must be REBOUND/Tudat, not our
  `kepler.py`) and the Lambert/closure code (must re-derive continuity in real
  dynamics, not re-run `correct.py`).

---

## 5. Validation of the harness itself (golden discipline)

The harness is *new code that produces science numbers*; per
`feedback_golden_tests_sourced_only`, its self-tests must not be circular. Five
checks, cheapest first:

1. **Two-body limit reduction (the keystone golden test).** Configure the n-body
   harness with **Sun only** (no perturbers) and propagate a state; it MUST equal
   our `core/kepler.py::propagate` to tight tolerance (both are then the same
   two-body problem). This is a *cross-implementation* golden test — neither side
   is the "sourced" side because both must agree with the analytic Kepler
   solution — and it directly catches a frame/time/units slip in the glue (§0).
2. **Energy / integrals conservation.** With Sun-only, the Keplerian energy and
   angular momentum are conserved; with perturbers, the Jacobi-like integral /
   total energy drift over a full propagation must stay below an integrator-
   accuracy floor (IAS15 should give ~machine precision energy conservation on
   the two-body sub-problem — a direct integrator-health check).
3. **Timestep / tolerance convergence.** Tighten the IAS15 accuracy parameter (or
   Tudat tolerance) until the final state and the reported correction ΔV stop
   moving below the consumer's tolerance; report the converged setting per
   consumer. This same sweep **doubles as the §2 body-inclusion sensitivity test**
   (re-run with/without Jupiter; if the metric moves < tolerance, exclusion is
   justified — evidence, not assertion).
4. **External ephemeris anchor (the independence check).** Propagate **Earth
   itself** as a test particle from its DE440 state and compare to the DE440
   ephemeris position over a year — OR, more directly, assert the harness's
   planet-state ingestion (after the §0 conversions) reproduces
   `Ephemeris("astropy").state(body, t_sec)` to numerical precision at sampled
   epochs. This is the anchor that proves the time/frame conversion (§0) is
   correct *before* any spacecraft number is trusted. The sourced side is DE440
   itself.
5. **Aldrin reproduction as an end-to-end golden anchor.** The harness, run on
   the Aldrin cell, must reproduce **our** 2.9138 km/s maintenance ΔV (consumer
   3) within tol — and, on the *independent* (GMAT/own-kernel) path, this is
   simultaneously the V4 result. (Golden caveat: 2.9138 is *our computed* value,
   so the in-pipeline reproduction is a consistency check, while the
   independent-tool reproduction is the genuine external validation. State which
   is which in the test.)

Mirror `verify/propagate.py`'s discipline: a frozen result dataclass (epochs,
defect norms, per-node correction ΔV, total correction ΔV, converged/diverged,
integrator settings used) so downstream gauntlet wiring fills fields rather than
reshaping, and divergence returns an honest non-converged record.

---

## 6. Phasing (recommended order, with reasons)

1. **Phase A — substrate + the §0 conversions + the §5 self-tests, on the
   simplest consumer (SILVER rungs, consumer 1).** Do this first because: (a) the
   rungs need only forward propagation + a correction solve (no STM, no GMAT),
   the lightest build; (b) it forces the time/frame anchor (§5.4) and the
   two-body golden (§5.1) to be correct *before* anything harder depends on them;
   (c) it produces an immediately useful result (rung 2–3 verdicts on the two
   held SILVER candidates) that the gauntlet can consume. REBOUND/IAS15 +
   shared-DE440.
2. **Phase B — multiple-shooting + Jones near-miss survey integration (consumer
   2).** Builds on Phase A's propagator behind the §3 interface; adds the
   node/defect/Jacobian machinery and consumes the Phase-3a near-miss conic
   survey output as seeds. This is where a **Tudat promotion** is decided (analytic
   STMs if finite-diff over REBOUND is too slow/noisy). Second because it is the
   hardest dynamics and gated on the conic-corrector Phase 1/2 landing
   (jones design §4) — no point shooting before there is a Jones-family seed.
3. **Phase C — GMAT V4 stamp (consumer 3), last.** A separate
   `scripts/gmat_v4_*` lane + parser, run **out of CI**, producing the formal V4
   record for Aldrin (and, optionally, for a Phase-B-corrected Jones/SILVER
   candidate that earned it). Last because: GMAT is the most operationally awkward
   (CI-hostile, scripted-app), it is a *stamp* not a loop component, and the
   in-pipeline Aldrin consistency check (§5.5) already de-risks the number before
   the external stamp is sought.

Rationale summary: **rungs → shooter → GMAT** = simplest-and-foundational →
hardest-dynamics → most-awkward-but-only-a-stamp. Each phase's self-tests gate the
next; the §0 conversions are proven in Phase A and reused unchanged downstream.

---

## 7. Open questions (verbatim, for the owner / plan author)

1. **REBOUND vs Tudat for consumer 2 (the shooter):** do we commit to REBOUND
   baseline + finite-difference STMs and only promote to Tudat if it proves too
   slow/noisy, or pay the Tudat install/CI cost up front to get analytic
   variational equations from the start? (Affects the venv/CI footprint and the
   Phase-B schedule.)

2. **`assist` dependency:** for the REBOUND path, do we take the `assist`
   companion (clean DE440 force eval + frames, extra dep) or hand-roll the
   ephemeris-force callback against the shared BSP (no extra dep, but we own the
   §0 conversions and the interpolation entirely)?

3. **Correction-ΔV convention (consumer 1):** node-impulse accounting vs
   distributed-TCM accounting — which is the canonical SILVER rung metric, and
   must it be made numerically identical to `maintain.py`'s per-synodic
   maintenance ΔV (the 2.9138 km/s convention) so the gauntlet compares like with
   like across consumers 1 and 3?

4. **Jupiter at the 100 m/s level (consumer 2):** is the §5.3 sensitivity test
   alone sufficient to settle Jupiter inclusion for the 12.8-yr Jones arc, or do
   we want to match Jones/AAS-17-577's *exact* body set if it is stated in the
   source (to keep the comparison model-faithful rather than tolerance-justified)?

5. **V4 tolerance for 2.9138 km/s:** what is "within tol" for the Aldrin
   maintenance-ΔV reproduction — an absolute m/s band, a percentage, or tied to
   the spread between our value and the published-but-ungapped Aldrin (the
   catalogue records the classic-Aldrin maintenance ΔV magnitude is *not*
   published, only a turn-angle test)? This sets the V4 pass/fail line and must be
   sourced or explicitly chosen by a human, not back-fit.

6. **Shared-vs-independent ephemeris for the rungs:** confirm the rungs run on the
   shared astropy DE440 BSP (cross-check rung, not V4) — and that a SILVER
   candidate's rung result is **never** recorded as a V4 pass. Is there appetite
   to *also* run the rungs on an independently-downloaded kernel to upgrade a
   strong rung result toward V4 without GMAT?

7. **GMAT in/out of CI:** confirm GMAT stays a manual/batched out-of-CI lane
   (like the perf-profile reproduction scripts) rather than a gating check —
   acceptable for a once-per-promising-candidate V4 stamp?

---

## Approval (2026-06-06)

User-approved with all recommendations accepted. The seven open questions
resolve as: (Q1) REBOUND-first with finite-diff STMs; promote to Tudat only
on demonstrated need (slow/noisy STM evidence, not anticipation); (Q2)
hand-roll the ephemeris force callback against the shared DE440 BSP — no
`assist` dependency until proven needed; (Q3) node-impulse correction-ΔV
convention, COMPARABLE to maintain.py's per-synodic value (same physical
meaning, documented mapping) but not forced numerically identical; (Q4) the
§5.3 sensitivity test governs Jupiter/body-set choices unless AAS 17-577
states its body set (check the mining note at execution); (Q5) the V4
tolerance for the 2.9138 km/s reproduction is a human-declared percentage
chosen BEFORE the run — 5% — never back-fit; (Q6) rungs run on the shared
DE440 as cross-check-only (never recorded as V4); the independent-kernel
upgrade is a later nice-to-have; (Q7) GMAT stays a manual, out-of-CI lane —
never a gating check. The split tooling (REBOUND in-pipeline / Tudat
conditional / GMAT stamp) and the planets-on-rails force model are adopted
as drafted. Phasing: rungs → shooter → GMAT, with §0 time/frame conversion
proofs FIRST. Next step: implementation plan.
