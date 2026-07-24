# #707 — CCR4BP catalogue schema design proposal for `#701`'s Umbriel-Titania torus/homoclinic connection

**Date**: 2026-07-24
**Scope**: design proposal only. No `data/catalogue.schema.json` or `data/catalogue.yaml`
change is made by this note. This is the last design step before a writeback decision on
`#701`/`#704`/`#705`'s CCR4BP homoclinic-connection result, dispatched in parallel with `#706`
(fresh literature re-check). Follows the `#684` precedent: research the conflict, present a
concrete recommendation, let the coordinating session/user decide before touching the schema.

## 1. What object is actually being represented

`#701` found a genuine (ghost-guard-verified: independent Radau/DOP853 integrator agreement to
~1e-7 km, off-torus distance ~1900-4900 km, quasi-Jacobi gap ~1e-14) **homoclinic connection**
in the idealized CCR4BP (Uranus + Umbriel, Titania as a periodic circular-orbit perturber): a
spacecraft riding a quasi-periodic invariant **torus** around a spacecraft:Umbriel 1:2-exterior
resonant periodic orbit can depart along the torus's **unstable manifold**, flow for
`t_u ~ 19.0` nondim time units (`~12.53` days), and return along the **stable manifold** of the
**same torus** after `t_s ~ 19.09` TU (`~12.59` days) — closing the idealized position/velocity
gap to machine precision (`corrected_pos_gap_km = 5.66e-10`, `corrected_vel_gap_km_s =
5.11e-14`).

`#704`/`#705` then showed this idealized connection does **not** survive verbatim under real
DE440-class ephemeris (headline 2030-01-01 gap: 174,709 km / 1.458 km/s) — but a **narrow
synodic-phase window** exists near each synodic cycle where the mismatch collapses to a
**comparable near-miss** (best 2030 point: 84.5 km / 5.9 m/s), and this narrow-window structure
**recurs at all 10 tested epochs** spanning 2000-2083 (best points 5.4-142.8 km / 1-13 m/s,
closest approach to either moon always far outside its physical radius — no collision-course
degeneracy). This is evidence of a **real, recurring, correctable near-connection**, not a
one-off numerical coincidence.

Structurally this is **unlike every other row in the catalogue**:

- It is not a discrete multi-encounter Lambert-arc tour (`#312`'s own family, every classical
  ballistic cycler): there is **no flyby of either moon**. Closest approach to Umbriel across all
  10 tested epochs is ~181,200-181,800 km (Umbriel's radius is ~584 km) and to Titania
  534,600-828,600 km (Titania's radius ~788 km) — both moons stay dynamically irrelevant to the
  trajectory's own geometry except as gravitational forcing terms. Titania is a perturber, never
  an encounter target.
- It is not a single strictly-periodic orbit either (the existing `cycler_class: non-keplerian`
  rows — Arenstorf, Genova-Aldrin, the Ross-RT/Braik-Ross CR3BP families, resonant_po rows): the
  underlying object is a **2-D invariant torus** (quasi-periodic, Titania-forced), and the
  connection itself is **two distinct manifold flow segments** (unstable departure + stable
  arrival), not one continuously-propagated periodic loop.
- It is epoch-locked in a real sense (the comparable near-miss window is narrow and recurs on the
  ~7.9-day Umbriel-Titania synodic period, exactly like `#312`'s own duty-cycle character) — so it
  is **not** eligible for the `#684`/schema-v5.2 epoch-free CR3BP KAM-corridor carve-out (that
  carve-out is specifically for rows with **no** real-ephemeris window; `#701`/`#704`/`#705`'s
  entire contribution is establishing that a real-ephemeris window *does* exist here).

## 2. Q1 — Does `orbit_class: quasi_cycler` fit?

**No. Recommend a new `orbit_class` value: `torus_homoclinic`.**

The schema's own description of `quasi_cycler` is explicit about what it was built for:
"closes-up-to-rotation, epoch_locked, 3-15 returns inside a planetary-alignment window (cyclers
of opportunity)". Every existing `quasi_cycler` row (the six `#312`/`#569` Uranian moon-pair
rows) is a **discrete alternating body-encounter tour** — `n_returns` there literally counts
consecutive repetitions of the same multi-leg flyby cycle, and the class's entire mission-
actionable meaning is "revisits these two named moons in this order, repeatedly, within this
calendar window."

`#701`'s object shares only the shallow, non-defining property of `quasi_cycler` (a bounded,
epoch-locked, real-ephemeris-characterized structure) but fails the two properties that actually
carry the class's meaning:

1. **No named-body encounter at all.** `quasi_cycler`'s "cyclers of opportunity" framing and its
   `sequence_canonical`/`vinf_kms_at_encounters` machinery presuppose a spacecraft that visits
   physical bodies. This object visits neither Umbriel nor Titania in any physically meaningful
   sense (per the closest-approach data above) — the "encounter" being characterized is
   Torus→(unstable manifold)→(stable manifold)→**the same Torus**, an entirely intangible,
   non-physical rendezvous.
2. **A qualitatively different repeat structure.** `quasi_cycler`'s repeats are *consecutive
   flights of the identical multi-leg tour*. `#701`'s "recurrence" is *the same one-shot
   torus-departure-and-return opportunity reappearing once per ~7.9-day synodic cycle* — there is
   no evidence (nor is it claimed) that a spacecraft could chain multiple homoclinic excursions
   end-to-end the way a `quasi_cycler` chains multiple flyby cycles.

I considered whether this could be handled the way `#684` handled its own genuine schema
conflict — a narrow, additive, conditional carve-out inside an existing enum value, rather than a
new value — since that is this schema's own established preference for minimizing enum
proliferation. I don't think it transfers here: `#684`'s carve-out worked because the object
being carved out (an epoch-free KAM corridor) was still **the same kind of thing** as existing
`quasi_cycler` rows (a bounded quasi-periodic neighborhood of an existing periodic orbit) — only
one axis (epoch-locking) differed. Here, the mission-actionable meaning itself (encounter-based
tour vs. non-encounter torus-homoclinic transfer) is different, and forcing it into
`quasi_cycler` would corrupt that enum's meaning for every future consumer that (reasonably)
assumes a `quasi_cycler` row has a real multi-body `sequence_canonical` and per-encounter
`vinf_kms_at_encounters`. A new value costs nothing (additive, no census-count change, no
existing-row impact) and keeps `quasi_cycler`'s crisp meaning intact.

I also considered `resonant_po` (no encounter utility, not epoch-locked, `n_returns='infinite'`)
since the underlying torus's base periodic orbit genuinely has that flavor — but `resonant_po`'s
own defining property is explicitly **"NO demonstrated transport utility"**, and `#701`'s object
is precisely the opposite: the whole point is a *transport* connection (depart the torus, arrive
back at it via a different route) — `resonant_po` would systematically hide the interesting part
of this result.

**Genuine uncertainty flagged honestly:** naming is a judgment call, not physics. I propose
`torus_homoclinic` (parallels the class of dynamical object — a homoclinic connection to a
torus — directly). `manifold_connection` or `torus_transfer` (both suggested in the dispatch)
are equally defensible; I'd resist `torus_transfer` specifically because "transfer" in this
catalogue's existing vocabulary (`inserts_into`, `precursor_mga`) already means "one-time MGA
insertion into a *different* steady-state cycler," which is not this — the connection returns to
the **same** structure it left, which is the single most important semantic fact about a
homoclinic (vs. heteroclinic) connection and should be visible in the name.

## 3. Q2 — Does `cycler_class: multi-arc` fit?

**No. Recommend `cycler_class: non-keplerian` (an existing value, not a new one), with the
distinguishing structure carried by the new `ccr4bp_provenance` block rather than a new
`cycler_class` enum value.**

`multi-arc` (per its usage precedent — `#312`'s own row, the Jovian/Saturnian Tisserand tours) is
for a chain of **distinct Kepler/Lambert arcs between different bodies**, each with its own
`center`/`tof_days_bounds`/V∞-continuity node. `#701`'s object has two flow segments too (unstable
departure, stable arrival), but they are not arcs **between different bodies** — they are the two
halves of **one** manifold structure attached to **one** torus. There is no midpoint encounter
node (no flyby, no V∞-continuity handoff to a second body) the way `#312`'s Umbriel→Titania and
Titania→Umbriel legs have. Forcing this into `multi-arc` would misrepresent it as a two-body
patched-conic tour when it is a single torus's own excursion-and-return.

`non-keplerian` is the schema's existing bucket for "a single rotating-frame dynamical structure,
not a chain of Kepler ellipses, not a chain of Lambert arcs" — exactly what this object's
*overall* identity is (one torus, one connection). It already satisfies the schema's structural
invariant (the `allOf` rule forcing `orbit_elements.a_au`/`e` to null for `multi-arc`/
`non-keplerian` — this object plainly has neither).

The one thing `non-keplerian` rows conventionally carry that doesn't cleanly fit is
`orbit_elements.cr3bp`'s single `{jacobi_constant, period_nd, stability_index, mass_ratio,
libration_point, family}` tuple, which describes **one strictly-periodic orbit's** identity. I
recommend populating that tuple with the **base resonant periodic orbit the torus is built
around** (Umbriel:spacecraft 1:2-exterior resonance; `period_nd` = the torus's own stroboscopic
period, `mass_ratio` = `mu`, `family` = a descriptive string), with `jacobi_constant` and
`stability_index` left null + a `data_gaps` entry (the CCR4BP has **no exact conserved quantity**
— see `quasi_jacobi_gap`'s own docstring — so a literal Jacobi constant doesn't exist for this
system the way it does for pure CR3BP; reporting the approximate `quasi_jacobi_gap` instead, in
`ccr4bp_provenance`, is the honest substitute). This reuses the existing field for what it can
honestly carry (the base orbit identity) rather than stretching it to describe the torus+
connection, which is what the new `ccr4bp_provenance` block is for.

I don't think a **new** `cycler_class` value is needed on top of this — the schema's own
`cycler_class` field is described as "structural kind of orbit" at a coarse level (three buckets
only, unchanged since schema v4), and the real distinguishing detail (torus vs. strict periodic
orbit, homoclinic manifold structure) is better carried by the additive `ccr4bp_provenance` block
than by inventing a fourth top-level structural bucket for what is, after all, still "not a
Kepler ellipse, not a discrete multi-body arc chain" at the `cycler_class` grain.

## 4. Q3 — What does `sequence_canonical` mean here?

**Recommend: null, with a `data_gaps` entry of `kind: "not-applicable"`** — following the *exact*
precedent `#312`'s own Umbriel-Titania-Umbriel row already set for `invariants.aphelion_ratio`
and `orbit_elements.cr3bp` (both marked `not-applicable` with an explanatory note, rather than
fabricating a value or silently omitting the field).

`sequence_canonical` (a free-form string field, not schema-constrained — `additionalProperties:
true` at the row level) is used exclusively, across all 381 rows that carry it, as a literal
alternating-body-flyby chain (`"E-M"`, `"Umbriel-Titania-Umbriel"`, `"E-E-M-M"`). There is no
alternating body sequence here — the "sequence" is Torus→unstable-manifold→stable-manifold→same
Torus, a single self-referential loop with no intermediate body encounter to name. Inventing a
pseudo-sequence string (e.g. `"Torus-Torus"` or `"Umbriel(1:2)-self"`) risks exactly the failure
mode this project's own memory flags repeatedly (grounding claims against content, not
inheriting/fabricating labels): any future downstream consumer that treats `sequence_canonical`
as a literal flyby chain (a reasonable assumption, given every existing usage) would silently
misread it.

I recommend leaving `sequence_canonical` null/absent with:
```yaml
data_gaps:
  - path: "sequence_canonical"
    kind: "not-applicable"
    note: >
      This is a torus-homoclinic connection (departs the SAME quasi-periodic torus along its
      unstable manifold, returns to it along its stable manifold) — there is no alternating
      named-body flyby sequence for this field's convention to describe (neither Umbriel nor
      Titania is closely encountered; see ccr4bp_provenance.real_ephemeris_evidence for the
      closest-approach sanity data). Left null rather than fabricated, per the #312 row's own
      precedent for invariants.aphelion_ratio / orbit_elements.cr3bp.
```

If the coordinating session wants a purely human-readable label (not a machine-meaningful
"sequence"), a separate free-text `name`/`notes` description (already present on every row) is
the right place for it, not `sequence_canonical`.

## 5. Q4 — What do `n_returns`/`validity_window` mean here?

This is the least clean-cut of the five questions; I'll state a recommendation and flag the
genuine uncertainty.

**`n_returns`: recommend `1`, not `'infinite'` and not a repeat-count like `#312`'s `10`.**
`#701`'s object is fundamentally a **single** departure-and-return use per opportunity — there is
no evidence (and it isn't claimed) that a spacecraft could chain multiple homoclinic excursions
back-to-back the way `#312`'s row chains 10 consecutive multi-arc cycles. Using `10` here (to
mirror `#312`'s convention, since `#705` tested 10 epochs) would be a **category error**: `#312`'s
`10` counts consecutive repetitions of the *same continuously-flown* tour within one calendar
span; `#705`'s `10` counts independently-tested, ~9-year-spaced **anchor epochs** used to
characterize recurrence — a sampling-density concept, not a repeat-count. Conflating the two
would silently misstate what was actually demonstrated. `n_returns: 1` is the honest claim: *one*
opportunity, characterized at *many* epochs.

**`validity_window`: recommend the coarse `{start: 2000, end: 2083}` calendar bounds (matching
the tested-epoch span), but explicitly NOT populating the existing
`synodic_duty_cycle_pct`/`synodic_boundary_period_days` sub-fields — flag this as a genuine
`data_gaps` entry instead.** Reasoning: those three sub-fields (schema v5.1) were formalized
around `#568`'s methodology — a **dense, continuous, daily-sampled** multi-year scan (`#312`'s own
row: N=731 daily samples across 2000+2030) reducing to a single feasible-fraction percentage
against **one binary pass/fail criterion** (V4-strict planet-crossing infeasibility — a hard
physical constraint). `#701`'s evidence is structurally different in **two** ways:

1. **Sparse epochs, not a continuous scan.** Only 10 discrete anchor epochs (~9-year spacing)
   were tested, each with its own dense 300-point *synodic-phase* scan — this characterizes the
   *within-cycle* duty structure well at each anchor, but says nothing about whether the
   near-miss quality drifts *between* the tested epochs (a slow secular effect with a period
   shorter than ~9 years, out of phase with all 10 samples, could in principle exist undetected —
   `#705`'s own stated honest caveat). Claiming a `synodic_duty_cycle_pct` computed only from
   these 10 anchors would silently overstate the density of the underlying evidence relative to
   `#312`'s own daily-sampled provenance for the *same* field.
2. **A continuous mismatch magnitude, not a binary constraint.** `#312`'s duty cycle is
   feasible/infeasible against a hard physical wall (does periapsis clear Uranus's radius).
   `#701`'s evidence (`#704`'s `duty_cycle_dense_synodic_scan`) is instead a **threshold ladder**
   over an inherently continuous position-mismatch magnitude (`below_500km` / `below_1000km` /
   `below_2000km` / ... / `below_50000km`, each with its own `fraction_of_period`) — because
   "feasible" here means "correctable within some correction-ΔV budget," and there is no single
   natural threshold the way there is a hard planetary radius. Forcing this into one
   `synodic_duty_cycle_pct` number would require arbitrarily picking one budget threshold and
   hiding the real multi-threshold structure that is actually the more informative evidence.

I recommend the ladder structure live in the new `ccr4bp_provenance` block (which is
purpose-built for exactly this kind of connection-specific evidence) rather than distorting the
existing `validity_window` sub-fields to fit a shape they weren't built for, and flag the
epoch-density gap explicitly:
```yaml
data_gaps:
  - path: "validity_window.synodic_duty_cycle_pct"
    kind: "unknown"
    note: >
      Only 10 discrete real-ephemeris anchor epochs (2000, 2009, ..., 2083; ~9-year spacing)
      have been tested (#705), each densely scanned across its own single ~7.9-day synodic
      cycle (300 points, #704's methodology) — NOT a continuous multi-year daily scan across
      the whole calendar span the way #312's own synodic_duty_cycle_pct was measured. A
      denser continuous scan (daily sampling across the full 2000-2083 span, matching #312's
      own evidence density) is a natural next step, not yet run. The within-cycle duty
      structure at each of the 10 tested anchors is instead recorded as a threshold ladder in
      ccr4bp_provenance.real_ephemeris_evidence (fraction of synodic phase below each of
      several position-mismatch bands), since "feasible" for this connection is a continuous
      correction-budget threshold, not #312's binary hard-physical-wall criterion.
    todo_ref: "#707"
```

**A further, adjacent open question this raises (not explicitly asked, but relevant to the
coordinating session): `validation_level`.** This object doesn't map cleanly onto the existing
V0-V5 ladder either. It has cleared the equivalent of solver cross-check (`#701`'s ghost guard:
independent Radau vs. DOP853 integrator agreement to ~1e-7 km) — a **V1**-equivalent claim. It has
**not** cleared V2 (there is no multi-lap periodic structure to test "≥3 continuous laps" against
— the object is a one-shot transfer, so V2's own criterion doesn't have an analogous instance to
apply to), and it has not cleared V3 (V3 requires a phase-matched real window with a
budget-bounded horizon TCM; `#704`/`#705` found a *near*-miss requiring correction but never ran a
budget-bounded correction burn to confirm the fix-up cost is closable under any stated budget —
that would be the natural next gate). I'd recommend `validation_level: V1` if a level is recorded
at all, with an explicit note that the V-ladder's V2/V3 criteria don't have a literal analogue for
a one-shot torus-homoclinic transfer and a bespoke gate (a documented TCM-budget-bounded
correction at one or more of the 10 recurring near-miss epochs) is the natural way to earn a
higher tier. This is a genuine open design question I'm flagging, not answering — it's adjacent
to but not the same question as `n_returns`/`validity_window`.

## 6. Q5 — The additive `ccr4bp_provenance` block

Following the `bcr4bp_provenance` precedent exactly (additive, optional, nullable,
`additionalProperties: true`, no schema-version-bump-forcing structural change — this would be
schema v5.3, purely additive like v5.0/v5.1/v5.2 before it).

One additional finding while designing this: **`model_assumption`'s enum has no CCR4BP value.**
The existing enum (`circular-coplanar`, `circular-inclined`, `analytic-ephemeris`, `cr3bp`,
`bicircular`, `null`) added `bicircular` specifically for BCR4BP; CCR4BP (Sun-absent, a THIRD body
on a circular orbit forcing the base two-body-restricted system — a structurally different
4-body model than BCR4BP's Sun-perturbation) needs its own value for the same reason. I recommend
adding `"ccr4bp"` to that enum alongside this change — using `"cr3bp"` for this row would misstate
the model (the whole point of `#701`'s result is the Titania-forcing term, which `cr3bp` doesn't
have), and would risk a cross-fidelity comparison bug analogous to the S1L1 precedent the
schema's own `orbit_fidelity` docstring warns about. `orbit_fidelity`/`vinf_fidelity`'s existing
tiers (`circular-coplanar` for the idealized CCR4BP numbers, `real-de440` for the `#704`/`#705`
real-ephemeris numbers) don't need a new value — those are fidelity TIERS, not model IDENTITY,
and CCR4BP's idealized outputs are still circular-coplanar-tier the same way CR3BP's are.

### Schema snippet (proposed, NOT applied)

```jsonc
"model_assumption": {
  "type": ["string", "null"],
  "enum": ["circular-coplanar", "circular-inclined", "analytic-ephemeris", "cr3bp", "bicircular", "ccr4bp", null],
  "description": "... (existing v4.6 text) ... Schema v5.3 (task #707) additively widens the enum with ccr4bp (Concentric Circular Restricted 4-Body Problem -- a base 3-body restricted system, e.g. Uranus-Umbriel, forced by a THIRD moon on its own circular orbit, e.g. Titania -- structurally distinct from bicircular's Sun-perturbation term). Reused unmodified for every existing row."
},

"ccr4bp_provenance": {
  "type": ["object", "null"],
  "description": "Schema v5.3 (task #707): ADDITIVE, OPTIONAL, NULLABLE provenance block for a CCR4BP (Concentric Circular Restricted 4-Body Problem) torus-homoclinic-connection row (orbit_class='torus_homoclinic'), admitted per the #689-#705 pipeline (core.ccr4bp / search.variational_ccr4bp_torus / search.ccr4bp_manifold_globalize / search.ccr4bp_heteroclinic_search). Mirrors bcr4bp_provenance's own additive pattern exactly. Absent/null for every other row (changes NO existing census count). Carries the CCR4BP system parameters, the idealized-model connection's own geometric identity, and the real-ephemeris recurrence evidence that the planar-CR3BP / real-eph row fields cannot express. Provenance/audit only -- not a promotion gate (validation_level carries that).",
  "additionalProperties": true,
  "properties": {
    "mu": {
      "type": ["number", "null"],
      "description": "Base moon's CR3BP mass ratio, GM_moon / (GM_primary_system + GM_moon) -- e.g. Umbriel's mass ratio in the Uranus-Umbriel two-body reduction. mu=0 has no special meaning here (unlike bcr4bp's mu_sun=0 CR3BP-exact reduction) since mu is the BASE system's own mass ratio, not the perturbation strength. Null/absent if unrecorded."
    },
    "mu_gan": {
      "type": ["number", "null"],
      "description": "Perturbing (third) moon's nondimensional mass in the base system's mass unit, GM_perturber / (GM_primary_system + GM_base_moon) -- e.g. Titania's mass ratio in the Uranus-Umbriel mass unit. Field name kept consistent with core.ccr4bp.CCR4BPSystem's own mu_gan attribute (historically named after Ganymede, the perturber in the first-built Jupiter-Europa-Ganymede system) even when the actual perturbing body is not Ganymede -- avoids a translation-layer name mismatch between this provenance block and the code that produced it. mu_gan=0 recovers CR3BP exactly (structural sanity check, per core.ccr4bp's own docstring). Null/absent if unrecorded."
    },
    "a_gan": {
      "type": ["number", "null"],
      "exclusiveMinimum": 0,
      "description": "Perturbing moon's semi-major axis in base-moon-SMA units (e.g. Titania SMA / Umbriel SMA, ~1.640 for Uranus-Umbriel-Titania). Same mu_gan-derived naming convention. Null/absent if unrecorded."
    },
    "omega_gan": {
      "type": ["number", "null"],
      "description": "Perturbing moon's synodic angular rate in the base system's rotating frame (n_perturber/n_base - 1; negative when the perturber is the slower-moving outer body, per core.ccr4bp.two_body_synodic_rate). Same mu_gan-derived naming convention. Null/absent if unrecorded."
    },
    "base_resonance": {
      "type": ["string", "null"],
      "description": "Descriptive tag for the base periodic orbit the torus is built around, e.g. 'spacecraft:Umbriel=1:2-exterior' -- echoed into orbit_elements.cr3bp.family as well; kept here too so ccr4bp_provenance is self-contained. Null/absent if unrecorded."
    },
    "torus_period_tu": {
      "type": ["number", "null"],
      "exclusiveMinimum": 0,
      "description": "Base torus's own stroboscopic-map period, nondimensional time units (TU) of the base CR3BP system. Echoed into orbit_elements.cr3bp.period_nd. Null/absent if unrecorded."
    },
    "torus_rho_strob": {
      "type": ["number", "null"],
      "description": "Stroboscopic-map rotation number (rho) of the invariant torus -- the quasi-periodic winding rate transverse to the base periodic orbit. Null/absent if unrecorded."
    },
    "torus_closure_residual": {
      "type": ["number", "null"],
      "minimum": 0,
      "description": "The GMOS-style (or equivalent) torus corrector's own closure residual at convergence, nondimensional -- the torus construction's own validity gate, analogous to corridor_measurement.closure_residual. Null/absent if unrecorded."
    },
    "connection": {
      "type": ["object", "null"],
      "description": "The idealized-CCR4BP homoclinic/heteroclinic connection's own geometric identity, as produced by search.ccr4bp_heteroclinic_search.RefinedConnection + GhostGuardReport.",
      "additionalProperties": true,
      "properties": {
        "theta2_u": {"type": ["number", "null"], "description": "Unstable-manifold departure phase (torus's second angle coordinate), radians."},
        "t_u_tu": {"type": ["number", "null"], "minimum": 0, "description": "Unstable-branch elapsed flow time, nondimensional TU."},
        "t_u_days": {"type": ["number", "null"], "minimum": 0, "description": "Same, converted to days via this row's own physical time unit."},
        "theta2_s": {"type": ["number", "null"], "description": "Stable-manifold arrival phase, radians."},
        "t_s_tu": {"type": ["number", "null"], "minimum": 0, "description": "Stable-branch elapsed flow time, nondimensional TU."},
        "t_s_days": {"type": ["number", "null"], "minimum": 0, "description": "Same, converted to days."},
        "idealized_pos_gap_km": {"type": ["number", "null"], "minimum": 0, "description": "Position mismatch at the refined connection's local optimum, idealized CCR4BP model, km."},
        "idealized_vel_gap_km_s": {"type": ["number", "null"], "minimum": 0, "description": "Velocity mismatch at the refined connection's local optimum, idealized CCR4BP model, km/s."},
        "off_torus_km": {"type": ["number", "null"], "minimum": 0, "description": "Ghost-guard off-torus sanity distance -- how far the departure state is from the unperturbed torus at the same elapsed time (rules out a trivial near-departure pseudo-match), km."},
        "quasi_jacobi_gap": {"type": ["number", "null"], "description": "Base-CR3BP Jacobi-constant gap between the two connection endpoints (search.ccr4bp_heteroclinic_search.quasi_jacobi_gap) -- an APPROXIMATE physical-plausibility check, not an exact conserved quantity (the CCR4BP has none). Reported, not gated."},
        "integrator_delta_km": {"type": ["number", "null"], "minimum": 0, "description": "Ghost-guard independent-integrator (Radau vs. DOP853) consistency check: the km difference between the two integrators' position-gap estimate for the same connection. Small = genuine physics, not a chaos-amplified numerical artifact."},
        "residual_norm": {"type": ["number", "null"], "minimum": 0, "description": "The 4-unknown least-squares refinement's own converged residual norm (nondimensional)."}
      }
    },
    "real_ephemeris_evidence": {
      "type": ["object", "null"],
      "description": "Real-ephemeris consistency + epoch-recurrence evidence for this connection (the #704/#705 findings). Distinct in kind from validity_window's synodic_duty_cycle_pct trio (v5.1) -- see the #707 design note Q4 for why this connection's evidence is a continuous correction-magnitude ladder at sparse discrete epochs, not a dense binary-feasibility scan.",
      "additionalProperties": true,
      "properties": {
        "force_model": {"type": ["string", "null"], "description": "Real-ephemeris force-model identity/kernel used for reconstruction, e.g. 'DE440 + URA111, n-body'."},
        "comparable_pos_gap_threshold_km": {"type": ["number", "null"], "minimum": 0, "description": "Pre-declared 'comparable to the reference near-miss' position-gap gate used to judge epoch recurrence (e.g. 10x the reference epoch's own best point)."},
        "comparable_vel_gap_threshold_km_s": {"type": ["number", "null"], "minimum": 0, "description": "Companion velocity-gap gate."},
        "tested_epochs": {
          "type": ["array", "null"],
          "description": "One entry per independently-tested real-ephemeris anchor epoch.",
          "items": {
            "type": "object",
            "additionalProperties": true,
            "properties": {
              "base_epoch_utc": {"type": "string", "description": "ISO-8601 anchor epoch this synodic-phase scan was centered on."},
              "local_minimum_epoch_utc": {"type": ["string", "null"], "description": "ISO-8601 epoch of the best (locally-minimized) near-miss point found within this anchor's synodic-phase scan."},
              "local_minimum_pos_gap_km": {"type": ["number", "null"], "minimum": 0},
              "local_minimum_vel_gap_km_s": {"type": ["number", "null"], "minimum": 0},
              "closest_approach_km": {
                "type": ["object", "null"],
                "description": "Closest approach distance to each dynamically-relevant body during this local-minimum trajectory, km -- the collision-course sanity check (both moons' physical radii are far smaller than every observed value here).",
                "additionalProperties": {"type": ["number", "null"], "minimum": 0}
              },
              "comparable_to_reference": {"type": ["boolean", "null"], "description": "Whether this epoch's local minimum passed the comparable_pos_gap_threshold_km/comparable_vel_gap_threshold_km_s gate."}
            }
          }
        },
        "n_epochs_tested": {"type": ["integer", "null"], "minimum": 0},
        "n_epochs_comparable": {"type": ["integer", "null"], "minimum": 0},
        "fraction_epochs_comparable": {"type": ["number", "null"], "minimum": 0, "maximum": 1, "description": "n_epochs_comparable / n_epochs_tested -- the epoch-recurrence verdict fraction."},
        "synodic_duty_cycle_bands_pct": {
          "type": ["object", "null"],
          "description": "Fraction of ONE synodic cycle (at a representative anchor epoch) for which the position mismatch is below each of several thresholds -- a continuous correction-budget ladder, NOT a single binary feasible/infeasible fraction (see validity_window.synodic_duty_cycle_pct's own docstring for the contrast). Keys are threshold labels (e.g. 'below_500km'), values are {fraction_of_period, n_distinct_windows}.",
          "additionalProperties": {
            "type": "object",
            "properties": {
              "fraction_of_period": {"type": ["number", "null"], "minimum": 0, "maximum": 1},
              "n_distinct_windows": {"type": ["integer", "null"], "minimum": 0}
            }
          }
        },
        "resolution_caveat": {"type": ["string", "null"], "description": "Honest statement of the evidence's own resolution limits (e.g. epoch spacing, points-per-cycle) -- mandatory whenever tested_epochs is non-null, mirroring dv_band_source's non-null-requires-source pattern."}
      }
    },
    "method": {
      "type": ["string", "null"],
      "description": "Corrector/search-pipeline identity/provenance tag, e.g. 'ccr4bp-manifold-globalize+heteroclinic-search+ghost-guard' (search.variational_ccr4bp_torus + search.ccr4bp_manifold_globalize + search.ccr4bp_heteroclinic_search). Null/absent if unrecorded."
    }
  }
}
```

## 7. Worked example — what `#701`'s row would look like (DRAFT, not committed)

```yaml
- id: umbriel-1-2-torus-homoclinic-uranus-2026
  name: "Umbriel 1:2-exterior resonant torus homoclinic connection (Titania-forced CCR4BP)"
  source: discovered
  trajectory_regime: ballistic
  model_assumption: ccr4bp   # PROPOSED new enum value, see Sec 6
  cycler_class: non-keplerian
  orbit_class: torus_homoclinic   # PROPOSED new enum value, see Sec 2
  epoch_locked: true
  n_returns: 1
  validity_window:
    start: "2000-01-01T00:00:00Z"
    end: "2083-01-01T00:00:00Z"
    # synodic_duty_cycle_pct / synodic_boundary_period_days deliberately NOT populated -- see
    # data_gaps below and Sec 4. synodic_period_days IS populated (well-defined, matches the
    # torus's own base period exactly -- see cross-check note in source_quotes).
    synodic_period_days: 7.908977583921671
  validation_level: V1   # see Sec 4's open validation_level discussion; V2/V3 have no literal
                          # analogue yet for a one-shot torus-homoclinic transfer
  source_ephemeris: "URA111 (JPL/NAIF Uranian satellite SPK kernel, 1900-2099 coverage) + DE440"
  orbit_source: derived
  vinf_source: derived
  orbit_fidelity: circular-coplanar
  vinf_fidelity: real-de440   # the #704/#705 real-ephemeris recurrence check is the row's
                               # headline evidence; the idealized CCR4BP numbers are
                               # circular-coplanar-tier (see ccr4bp_provenance.connection)
  bodies: ["Uranus", "Umbriel", "Titania"]   # Titania modeled (perturber), never encountered
  sequence_canonical: null   # see Sec 3 -- data_gaps entry below
  sense: "n/a"
  orbit_elements:
    reference_frame: rotating-synodic
    center: Uranus
    a_au: null
    e: null
    cr3bp:
      jacobi_constant: null   # CCR4BP has no exact conserved quantity -- see data_gaps
      period_nd: 11.991104994404195
      stability_index: null
      mass_ratio: 1.4685981867903303e-05   # mu (Uranus-Umbriel)
      libration_point: null
      family: "umbriel-spacecraft-1:2-exterior-resonant-torus"
      lunit_km: 265986.0
  ccr4bp_provenance:
    mu: 1.4685981867903303e-05
    mu_gan: 3.915686587341081e-05
    a_gan: 1.6403043769220937
    omega_gan: -0.5239871813408119
    base_resonance: "spacecraft:Umbriel=1:2-exterior"
    torus_period_tu: 11.991104994404195
    torus_rho_strob: 5.995568015306847
    torus_closure_residual: 0.00014270615699557066
    connection:
      theta2_u: 3.045143535866834
      t_u_tu: 19.002858810413148
      t_u_days: 12.533722657930381
      theta2_s: 3.2367817357763884
      t_s_tu: 19.092874545266827
      t_s_days: 12.593094369669227
      idealized_pos_gap_km: 5.661540346264633e-10
      idealized_vel_gap_km_s: 5.1129191191902014e-14
      off_torus_km: 1927.9247765134512
      quasi_jacobi_gap: -1.2434497875801753e-14
      integrator_delta_km: 1.0512238869207226e-07
      residual_norm: 1.1159187446079244e-14
    real_ephemeris_evidence:
      force_model: "DE440 + URA111 n-body (see #704/#705 driver scripts for exact integrator)"
      comparable_pos_gap_threshold_km: 844.6019822482435   # 10x #704's own 2030 best point
      comparable_vel_gap_threshold_km_s: 0.05904706092575937
      tested_epochs:
        - base_epoch_utc: "2000-01-01T00:00:00"
          local_minimum_epoch_utc: "2000-01-01T11:52:25"
          local_minimum_pos_gap_km: 43.32017412945477
          local_minimum_vel_gap_km_s: 0.013072019096617437
          closest_approach_km: {Umbriel: 181621.49263441944, Titania: 810695.1962968499}
          comparable_to_reference: true
        - base_epoch_utc: "2009-03-23T03:59:59"
          local_minimum_epoch_utc: "2009-03-24T00:02:35"
          local_minimum_pos_gap_km: 60.98124123293917
          local_minimum_vel_gap_km_s: 0.008812134375914481
          closest_approach_km: {Umbriel: 181623.5, Titania: 825874.5}
          comparable_to_reference: true
        # ... (remaining 8 tested epochs: 2018, 2027, 2036, 2046, 2055, 2064, 2073, 2083 --
        #      all 10/10 comparable; full per-epoch data in
        #      data/found/705_ccr4bp_epoch_robustness_scan/result.json's epoch_scans[])
      n_epochs_tested: 10
      n_epochs_comparable: 10
      fraction_epochs_comparable: 1.0
      synodic_duty_cycle_bands_pct:
        below_500km: {fraction_of_period: 0.0033333333333333335, n_distinct_windows: 1}
        below_1000km: {fraction_of_period: 0.0033333333333333335, n_distinct_windows: 1}
        below_2000km: {fraction_of_period: 0.01, n_distinct_windows: 3}
        below_5000km: {fraction_of_period: 0.01, n_distinct_windows: 3}
        below_10000km: {fraction_of_period: 0.03333333333333333, n_distinct_windows: 4}
        below_20000km: {fraction_of_period: 0.06333333333333334, n_distinct_windows: 5}
        below_50000km: {fraction_of_period: 0.17333333333333334, n_distinct_windows: 6}
      resolution_caveat: >
        Only 10 discrete anchor epochs tested (~9-year spacing, 2000-2083), each with a 300-point
        dense synodic-phase scan (#704/#705 methodology) -- not a continuous multi-year daily
        scan. A slow secular effect with a period shorter than ~9 years and out of phase with all
        10 samples could in principle still exist undetected, though 10/10 hits with no
        near-misses among the 10 makes this unlikely.
    method: "ccr4bp-manifold-globalize+heteroclinic-search+ghost-guard"
  data_gaps:
    - path: "sequence_canonical"
      kind: "not-applicable"
      note: >
        Torus-homoclinic connection (departs the SAME torus's unstable manifold, returns via its
        stable manifold) -- no alternating named-body flyby sequence exists to describe. Neither
        Umbriel nor Titania is closely encountered (see ccr4bp_provenance.real_ephemeris_evidence
        closest_approach_km, always >181,000 km / >534,000 km respectively, far outside either
        moon's physical radius). Left null rather than fabricated, per the #312 row's own
        precedent for invariants.aphelion_ratio / orbit_elements.cr3bp.
      todo_ref: "#707"
    - path: "orbit_elements.cr3bp.jacobi_constant"
      kind: "not-applicable"
      note: >
        The CCR4BP (time-periodic Titania forcing) has NO exact conserved quantity -- see
        search.ccr4bp_heteroclinic_search.quasi_jacobi_gap's own docstring. The approximate
        substitute (ccr4bp_provenance.connection.quasi_jacobi_gap) is reported there instead.
      todo_ref: "#707"
    - path: "invariants.aphelion_ratio"
      kind: "not-applicable"
      note: "Uranus-centered torus-homoclinic connection, not a single heliocentric ellipse."
      todo_ref: "#707"
    - path: "validity_window.synodic_duty_cycle_pct"
      kind: "unknown"
      note: >
        Only 10 discrete real-ephemeris anchor epochs tested (#705), each densely scanned across
        its own single synodic cycle (#704) -- NOT a continuous multi-year daily scan the way
        #312's own synodic_duty_cycle_pct was measured. Recorded instead as a per-anchor
        threshold ladder in ccr4bp_provenance.real_ephemeris_evidence.synodic_duty_cycle_bands_pct
        (a continuous correction-budget concept, not #312's binary hard-physical-wall criterion).
        A denser continuous scan matching #312's own evidence density is a natural next step.
      todo_ref: "#707"
  first_published:
    authors: ["cyclerfinder discovery campaign"]
    year: 2026
    title: "Umbriel 1:2-exterior resonant torus homoclinic connection -- CCR4BP (Uranus-Umbriel-Titania)"
    venue: "cyclerfinder project; task chain #689 -> #693 -> #694 -> #699 -> #701 -> #702 -> #704 -> #705 -> #706 -> #707"
  priority_date: "2026-07-24"
  notes: |
    DRAFT ROW -- illustrative only, produced by #707's schema design proposal, NOT a writeback.
    Pending #706's literature-clearance re-check and coordinating-session/user sign-off on the
    orbit_class/cycler_class/ccr4bp_provenance schema questions this note raises.
```

## 8. Summary of recommendations

| Question | Recommendation | Confidence |
|---|---|---|
| Q1 `orbit_class` | New value `torus_homoclinic` (not `quasi_cycler`, not `resonant_po`) | High on "not quasi_cycler"; medium on exact naming |
| Q2 `cycler_class` | Existing `non-keplerian` (not `multi-arc`, no new enum value needed) | High |
| Q3 `sequence_canonical` | `null` + `data_gaps` `not-applicable`, following `#312`'s own precedent | High |
| Q4 `n_returns`/`validity_window` | `n_returns: 1`; calendar `validity_window` bounds only, NOT the `synodic_duty_cycle_pct` trio (evidence shape mismatch) — flag as `data_gaps` | Medium — genuine judgment call, flagged honestly |
| Q5 `ccr4bp_provenance` | New additive block per Sec 6 snippet, PLUS a `model_assumption` enum addition (`"ccr4bp"`) discovered as a necessary companion change | High on structure; the `model_assumption` addition is a new finding worth separate confirmation |

No schema or catalogue file has been modified by this task. All of the above is a proposal for
the coordinating session and user to review, decide, and (if approved) implement as an explicit
schema-bump + writeback task, exactly as `#684`'s own schema question was handled.
