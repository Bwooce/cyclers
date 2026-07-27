# #735 — N=5 CRNBP torus catalogue schema design proposal (`#720`/`#724`/`#729` object)

**Date**: 2026-07-27
**Scope**: design proposal only. No `data/catalogue.schema.json` or `data/catalogue.yaml`
change is made by this note. This is V0-V5 vetting-chain step 3 (the `#707` analogue) for the
`#714`->`#717`->`#720`->`#721`/`#722`/`#724`->`#726`->`#729` N=5 CRNBP torus. Follows `#707`'s
own precedent exactly: research the conflict, present concrete options WITH tradeoffs, recommend,
let the coordinating session/user decide before touching the schema. Sources read in full for
this design: `docs/notes/2026-07-27-724-final-confirmation-n5-torus-novelty.md` (mandatory claim
language), the `#726`/`#729` reports in `data/OUTSTANDING.md` + `data/found/
729_crnbp_epoch_torus_robustness_scan/result.json` (the actual evidence structure this design
must hold), `docs/notes/2026-07-24-707-ccr4bp-catalogue-schema-design.md` (process template),
`data/catalogue.schema.json` (current v5.3 in full), the live
`umbriel-1-2-torus-homoclinic-uranus-2026` row, `src/cyclerfinder/core/crnbp.py`, and
`src/cyclerfinder/search/variational_crnbp_torus.py`.

## 1. What object is actually being represented

A **bare quasi-periodic invariant 2-torus** — no connection, no manifold computation of any
kind — in a genuine N=5 restricted model (spacecraft + Jupiter + Europa + Io + Ganymede,
`core.crnbp`, the Laplace-locked Jupiter-Europa frame). `#724`'s **exact, mandatory claim
language**, quoted verbatim (this design must be able to carry it without loosening it):

> First computed quasi-periodic invariant-torus substitute of a mean-motion-
> resonant periodic orbit in a Laplace-locked Jupiter-Io-Europa-Ganymede
> restricted five-body model: the planar Jupiter-Europa 3:4 resonant orbit of
> Kumar et al. 2021 (exterior to Europa, a ≈ 1.211 Europa SMA), continued to
> the physical Io mass with Ganymede at its physical (non-rate-idealized)
> synodic rate and Io exactly Laplace-slaved (ω_Io = -2·ω_Gan), at the
> physical libration-center phase (Φ_L = 180°), via 2D pseudospectral Fourier
> collocation (n1=2, n2=20; on-grid residual RMS 1.23e-4, off-grid invariance
> ~2.1e-3).

With `#724`'s mandatory qualifiers, all of which the draft row in §9 carries verbatim:
- The MODEL class is Baresi, Owen & Scheeres's (AAS 23-201, 2023; ISSFD 2024) — cite both; the
  Tri-Circular/Laplace-locked idea is theirs, and they computed the first N=5 tori (L1/L2
  Lyapunov-family substitutes), Floquet stability, manifolds, and Europa<->Ganymede transfers.
- Exactly TWO novelty axes: (i) orbit family substituted (mean-motion-resonant, Kumar-class, vs
  L1/L2 Lyapunov planar); (ii) rate model (physical Ganymede synodic rate vs both rates
  idealized to exact rationals). Method is a difference but NOT independently novelty-bearing;
  phase is NOT a differentiator (TCP already used 180°).
- "Exterior Jupiter-Europa 3:4" or plain "Jupiter-Europa 3:4"; NEVER "interior".
- **No stability, manifold, or transfer computation exists for OUR object** (TCP has all three
  for theirs) — the row must not imply otherwise.
- Never "first N=5 CRNBP torus" / "first torus in the Laplace-locked model".
- `search/literature_check.py`'s gate still applies at actual writeback time (not run by this
  design task).

Real-ephemeris standing (`#726`/`#729`, `data/found/729_crnbp_epoch_torus_robustness_scan/
result.json`): the idealized torus does **NOT** survive real SPICE ephemeris generically
(`#726` headline at torus point θ1=θ2=0, naive 2030-01-01 epoch, one Ganymede-synodic forcing
period = 12.478 TU / 7.05 d: pos_gap 3.68e5 km / vel_gap 8.19 km/s — comparable to Ganymede's
orbital radius; verified genuine vs an idealized-fed noise floor of ~11% of that). But `#729`'s
epoch(+torus-point) scan (10 epochs 2000-2083, 300-point dense synodic-phase scan + bisection
refine per epoch, mirroring `#705`'s methodology, PLUS a systematic 5-point torus-point axis)
found a **recurring narrow near-miss window**: at torus points (0,0) and (π,π) the per-epoch
refined local minima are 503-1950 km and 790-1719 km respectively at ALL 10 epochs (200-2500x
tighter than the generic collapse); (π/2,π/2) also narrow 10/10 (2410-4042 km); (π,0) and (0,π)
genuinely collapse 10/10 (2.6e4-1.3e6 km). **The narrow/collapse dichotomy is
torus-point-dependent but PERFECTLY epoch-stable** — `#729`'s own strongest evidence that this
is real dynamical structure, not noise. This is the same recurring-window pattern that
ultimately supported the Umbriel-Titania writeback, but the *measurement is structurally
different*: the gap here is a **model-shadowing mismatch** (idealized CRNBP forward flow of a
chosen torus point over one forcing period vs real-SPICE propagation of the same departure
state), NOT a connection-quality metric — there is no departure/target manifold pair, and no
`comparable_to_reference` gate against a prior single-epoch connection (the pre-declared gate is
`#729`'s own absolute 5000 km narrow-near-miss threshold). §7's evidence block is shaped around
exactly this difference.

**What this object is NOT** (the load-bearing distinction for the whole design): it is not
`#708`'s object one system over. The Umbriel-Titania row is a torus-HOMOCLINIC-CONNECTION — a
computed, ghost-guard-verified manifold intersection (unstable departure + stable return to the
SAME torus, closed to machine precision) — and the `torus_homoclinic` class, its `n_returns=1`
semantics, and the `ccr4bp_provenance.connection{}` payload are all built around that connection.
For the N=5 object **no manifold has been globalized and no connection has been attempted**
(`#714` shortlist item 3, explicitly gated and not dispatched; `#720`'s own discipline section:
"No whisker/manifold globalization, no heteroclinic search"). The object is the torus itself,
plus real-ephemeris shadowing evidence for it.

## 2. Q1 — Does `orbit_class: torus_homoclinic` fit?

**No — and this is the highest-confidence answer in this note.**

The schema's own v5.3 definition is unambiguous about what the class IS: "a homoclinic
connection to a quasi-periodic invariant torus in a CCR4BP (departs the torus's unstable
manifold, returns via its stable manifold to the SAME torus) — epoch_locked=true, n_returns=1
(ONE departure-and-return opportunity per recurring real-ephemeris synodic window...)". Three
independent disqualifiers:

1. **The defining object is the connection, not the torus.** Every distinguishing property in
   the class description (departs/returns, the two manifold segments, the `n_returns=1`
   "departure-and-return opportunity" semantics) describes a computed manifold intersection.
   None has been computed here. Labelling this row `torus_homoclinic` would be a factual
   overclaim of exactly the kind `#724`'s mandatory-qualifier list bars ("No stability,
   manifold, or transfer computation exists for OUR object yet — do not imply otherwise").
2. **The class's payload contract would be violated.** The one existing `torus_homoclinic` row
   carries a fully-populated `ccr4bp_provenance.connection{}` (phases, flow times, machine-
   precision gap closure, ghost-guard fields). A consumer that reads `orbit_class:
   torus_homoclinic` may reasonably assume that payload exists. For this row every `connection`
   field would be null — the class's defining evidence, absent.
3. **`n_returns=1` has no referent.** There is no departure-and-return opportunity to count.
   (§5 discusses what `n_returns` honestly means for a bare torus.)

Could the class be *widened* ("`torus_homoclinic`, connection optional")? Rejected for the same
reason `#707` rejected widening `quasi_cycler`: it would corrupt the enum's crisp meaning for
its existing row and for every future consumer, to save one enum value that costs nothing to
add. If/when `#714` item 3 (N=5 manifold/connection work) is eventually run AND succeeds, the
natural outcome is a *separate* connection row (or a supersedes/superseded_by pair) — not a
retroactive blurring of this class.

## 3. Q2 — Does `orbit_class: resonant_po` fit?

**No — but it is the strongest existing-value candidate and deserves the honest case for it
first.** Interestingly, `#707` rejected `resonant_po` for the Umbriel object because that
object's whole point was demonstrated transport ("`resonant_po`'s defining property is
explicitly 'NO demonstrated transport utility', and `#701`'s object is precisely the
opposite"). The N=5 torus **fails that objection's premise**: it genuinely has no demonstrated
transport utility — no connection, no manifold, no encounter. On the transport axis,
`resonant_po` fits this object exactly as well as it fits `em-cycler-21-3d-spatial-2026`. It is
also literally built around a mean-motion-**resonant** orbit (the Kumar-class Europa 3:4), and
its catalogue purpose (carry a computed dynamical object for the record, not as a usable
cycler) matches.

But it fails on three grounds, two of them hard:

1. **Encodability: `resonant_po` presupposes a strictly-periodic orbit, and a torus cannot be
   honestly encoded under its field conventions.** The class description says "periodic orbit"
   throughout, and the existing member's identity payload is `orbit_elements.cr3bp`'s
   strictly-periodic tuple: a single `state_nd` initial condition, a single `period_nd`, a
   `jacobi_constant`, a Floquet `stability_index`. For this object: there is no single IC (the
   object is a 2-parameter family of states, a Fourier coefficient array `coeffs[4][2n1+1]
   [2n2+1]`); "period" is only the stroboscopic forcing period (the motion itself is
   quasi-periodic with irrational rotation number ~0.4965); the CRNBP has **no conserved
   quantity at all** (time-periodic forcing — same reason the Umbriel row nulls
   `jacobi_constant`); and no Floquet stability has been computed. The JSON Schema would
   mechanically pass (those fields are nullable) — but *every field the class exists to carry*
   would be null or reinterpreted, which is the definition of a class misfit. This is the
   dispatch's own "check the schema's field requirements" question, answered: legal, but only
   vacuously.
2. **The epoch invariant actively misstates the evidence.** The v4.7/v4.9 ratchet
   (`tests/data/test_schema_v47_orbit_class.py`) hard-requires `resonant_po` =>
   `epoch_locked=false` / `n_returns='infinite'` — "reachable at any epoch and repeats
   indefinitely". That is true of the *idealized* torus, but this row's headline real-world
   evidence is `#726`'s **generic collapse** plus `#729`'s **narrow recurring windows**: under
   real ephemeris the object is emphatically NOT reachable at any epoch/phase — model-shadowing
   holds only at specific torus points within ~0.7% of the synodic cycle (below-5000 km duty
   fraction at the 2000 anchor: 0.0067). Encoding it `epoch_locked=false` would erase the
   single most important thing `#726`/`#729` established, exactly the failure `#707` Q4 was
   designed to avoid. (One could carve out an exception as v5.2 did for `quasi_cycler` — but
   v5.2's carve-out *relaxed* toward epoch-freedom for an object with NO real-ephemeris window;
   this would be the reverse carve-out, *tightening* toward epoch-locking, for a class whose
   only other member is epoch-free. Two members, two opposite invariants: the class would mean
   nothing.)
3. **Softer: purpose mismatch.** All 22 `resonant_po`-adjacent rows are
   `our_status: known-class-member` corroboration entries of published classes. This row is a
   claimed-novel (narrow claim, `#724`) discovery with a real-ephemeris evidence chain — the
   class's "carried for known-class corroboration" framing undersells it, though this alone
   would not disqualify.

## 4. Also considered: `quasi_cycler` (v5.2 epoch-free KAM-corridor carve-out)

Not in the dispatch's three questions, but it is the *other* place the catalogue already keeps
bare quasi-periodic tori, so it must be explicitly ruled out. The v5.2 carve-out
(`corridor_measurement`, `#682`/`#684`) admits epoch-free CR3BP KAM tori — but its own
diagnostic condition is `model_assumption='cr3bp'` AND `validity_window` null (no
real-ephemeris window), and its meaning is "a measured corridor AROUND an already-catalogued,
linearly-stable periodic orbit — a characterization, NOT a novel discovery (KAM theory
guarantees existence)". This object fails every prong: not CR3BP (the forcing IS the object's
point); it HAS a real-ephemeris window (the `#729` evidence is the row's core); it is not a
neighborhood characterization of an existing catalogued row; and KAM persistence at the
physical `mu_Io` was NOT guaranteed a priori (that is what the `#720` continuation
established). The non-carve-out `quasi_cycler` branch fails for `#707`'s original reasons (no
named-body encounter sequence, no repeat-flown tour). Rejected.

## 5. Recommendation: a new `orbit_class` value

**Recommend a new enum value. Preferred name: `quasi_periodic_torus`.** Reasoning mirror of
`#707` Sec 2: a new value is additive, costs no census change and no existing-row impact, and
every existing value would either overclaim (`torus_homoclinic`), misstate the epoch evidence
(`resonant_po`), or presuppose machinery this object lacks (`quasi_cycler`). Naming options,
honestly weighed:

| Candidate | For | Against |
|---|---|---|
| `quasi_periodic_torus` (recommended) | Names the dynamical object class directly, model-agnostic — also fits a future CCR4BP/QBCP/BCR4BP bare-torus row (e.g. a known-reproduction of Kumar et al. 2021's own N=4 tori, or of TCP's Lyapunov substitutes) without another enum bump | Longest name; must be carefully distinguished in its description from the v5.2 KAM-corridor `quasi_cycler` subclass (done in the proposed description below) |
| `crnbp_torus` | Unmistakably THIS object | Encodes model identity into the object-class axis — the schema already separates those (model lives in `model_assumption`); a future N=4 bare-torus row would force yet another value |
| `resonant_torus` | Captures the novelty axis (resonant-family substitute) | Excludes future Lyapunov-family torus substitutes (TCP's own class) for no gain; the resonance identity already lives in `orbit_elements.cr3bp.family` / provenance `base_resonance` |

Proposed enum-description text (additive, would be schema v5.4):

> `quasi_periodic_torus` (schema v5.4, task #735) = a bare quasi-periodic invariant torus — a
> dynamical substitute of a periodic orbit under time-periodic forcing (CCR4BP/CRNBP/QBCP
> class) — with NO computed manifold connection. Distinct from `torus_homoclinic` (which IS a
> computed connection; a `torus_homoclinic` row's `connection{}` payload is populated, a
> `quasi_periodic_torus` row has none) and from the v5.2 epoch-free CR3BP KAM-corridor
> `quasi_cycler` subclass (which characterizes the corridor around an already-catalogued
> stable CR3BP periodic orbit, is `model_assumption='cr3bp'`-scoped, and has no real-ephemeris
> window — whereas a `quasi_periodic_torus` row's forcing model is the object's point and its
> real-ephemeris shadowing evidence, where present, is the row's headline). Distinct from
> `resonant_po`: the underlying object is a 2-D torus with irrational rotation number, not a
> strictly-periodic orbit, and its real-ephemeris character is window-structured, not
> reachable-at-any-epoch.

### `epoch_locked` / `n_returns` invariant for the new class (genuine judgment call)

Two defensible options — the same tension `#707` Q4 hit, sharpened because the idealized object
here really does wind forever:

- **Option A (recommended): `epoch_locked=true`, `n_returns=1`.** Semantics, stated in the enum
  description and enforced by extending the existing Python ratchet: *one forcing period of
  demonstrated real-ephemeris model-shadowing per recurring narrow window*. This has a real
  referent — `#726`/`#729`'s check propagates exactly ONE Ganymede-synodic forcing period
  (12.478 TU / 7.05 d) and measures the return-to-model gap; the demonstrated unit of validity
  is one period per window, recurring. It mirrors `torus_homoclinic`'s own gate shape, keeps
  the "the real-ephemeris window IS the row's evidence" principle from `#707` Q4 (the object
  generically collapses off-window — `epoch_locked=false` would misstate this, per §3 point 2),
  and avoids fabricating an infinite real-world persistence claim nobody has tested.
- **Option B: `epoch_locked=false`, `n_returns='infinite'`.** Honors the *idealized* torus's
  own character (it exists at any epoch in the model, like `resonant_po`/the KAM carve-out),
  pushing all real-ephemeris structure into the provenance block. Rejected as the
  recommendation because the top-level fields are what census consumers filter on, and
  "epoch-free" is the opposite of what `#726` measured — but listed because it is the
  intellectually cleanest reading of "what is the object *an sich*".

Medium confidence; this is the design's largest judgment call and is explicitly flagged for the
user, exactly as `#707` flagged its Q4.

### `cycler_class`

**Existing `non-keplerian`, no new value** — identical reasoning to `#707` Sec 3, and strictly
easier here (one torus, not even a two-segment connection; plainly "a single rotating-frame
dynamical structure, not a chain of Kepler/Lambert arcs"). Satisfies the `allOf` null-`a_au`/`e`
invariant. `orbit_elements.cr3bp` carries the base resonant orbit's identity (family string,
stroboscopic `period_nd`, `mass_ratio=mu`), with `jacobi_constant`/`stability_index` null +
`data_gaps` (no conserved quantity in the CRNBP; no Floquet computed — mandatory honesty per
`#724`).

### `model_assumption`

**New enum value `"crnbp"` needed** (companion change, exactly parallel to `#707`'s discovery
that `ccr4bp` was needed). Using `"ccr4bp"` would misstate the model — the entire first novelty
axis is the second simultaneous perturber. Proposed description: "Circular Restricted N-Body
Problem (Negri & Prado 2022 / Gilliam 2024 lineage; `core.crnbp`): a base CR3BP forced by TWO
OR MORE additional bodies on concentric circular coplanar orbits. N is not encoded in the enum
value — it is `2 + 1 + len(crnbp_provenance.perturbers)`; the one admitting row (N=5,
Jupiter-Europa base + Io + Ganymede) is Laplace-locked to a single forcing clock (see
`crnbp_provenance`)." This future-proofs the enum against N=6+ without another bump.
`orbit_fidelity: circular-coplanar` for the idealized numbers / `vinf_fidelity` discussion
carries over from `#707` unchanged (fidelity TIERS, not model identity).

## 6. The provenance block: `crnbp_provenance` (new), not an extended `ccr4bp_provenance`

Three options considered:

- **P1 — extend `ccr4bp_provenance` with `mu_io`/`a_io`/`omega_io`/`theta_io0`/`theta_gan0`
  optional fields.** Pros: one block, one consumer path. Cons: the block's own v5.3 description
  scopes it to `orbit_class='torus_homoclinic'` rows and its centerpiece is `connection{}` —
  every N=5 row would carry a permanently-null defining payload; the `mu_gan`-naming convention
  ("the perturber", singular) breaks at two perturbers; and each future N demands another
  hand-named field triple. The `#707` naming principle ("field name kept consistent with the
  code's own attribute to avoid a translation-layer mismatch") now argues the OTHER way: the
  N=5 code's own structure is `CRNBPSystem(mu, perturbers=(CRNBPPerturber(mu, a, omega,
  theta0), ...))` — a *list*, not named scalars.
- **P2 (recommended) — a new additive `crnbp_provenance` block with an N-body perturber LIST**,
  mirroring `core.crnbp.CRNBPSystem`/`CRNBPPerturber` 1:1. Pros: code-to-catalogue mapping is
  exact; future-proof for any N (including, if ever wanted, re-expressing an N=4 system as a
  one-element list — `ccr4bp_provenance` stays frozen for its one existing row); the torus and
  real-ephemeris sub-blocks can be shaped for what was actually measured instead of inheriting
  connection-shaped fields. Cons: a second parallel block — a cross-row consumer comparing the
  Umbriel row and this row must handle both shapes (mitigated: the `real_ephemeris_evidence`
  per-epoch item shape is kept deliberately field-compatible, §7).
- **P3 — generalize `ccr4bp_provenance` in place (migrate the Umbriel row to a list form).**
  Rejected outright: violates the project's additive-only schema discipline, forces a
  writeback-touching migration of a shipped, user-approved row, and buys nothing P2 doesn't.

### Proposed schema snippet (NOT applied; would be schema v5.4)

```jsonc
"crnbp_provenance": {
  "type": ["object", "null"],
  "description": "Schema v5.4 (task #735): ADDITIVE, OPTIONAL, NULLABLE provenance block for a CRNBP (Circular Restricted N-Body Problem, N>=5) row (orbit_class='quasi_periodic_torus'), admitted per the #714-#729 pipeline (core.crnbp / search.variational_crnbp_torus / search.crnbp_real_ephemeris_consistency). Mirrors ccr4bp_provenance's additive pattern but generalizes the perturber parameters to a LIST (matching core.crnbp.CRNBPSystem.perturbers exactly), so N=6+ needs no schema change. Absent/null for every other row (changes NO existing census count). Provenance/audit only -- not a promotion gate.",
  "additionalProperties": true,
  "properties": {
    "mu": {
      "type": ["number", "null"],
      "description": "BASE system's CR3BP mass ratio, GM_secondary / (GM_primary_system + GM_secondary) -- e.g. Europa's ratio in the Jupiter-Europa reduction (core.crnbp.CRNBPSystem.mu)."
    },
    "base_pair": {
      "type": ["string", "null"],
      "description": "Human-readable base-system tag, e.g. 'Jupiter-Europa'. The rotating frame, DU, and TU are this pair's."
    },
    "perturbers": {
      "type": ["array", "null"],
      "description": "One entry per extra perturbing body, IN THE SAME ORDER as core.crnbp.CRNBPSystem.perturbers in the producing run (order is load-bearing for reproduction). Each mirrors core.crnbp.CRNBPPerturber field-for-field, plus provenance strings.",
      "items": {
        "type": "object",
        "additionalProperties": true,
        "required": ["body", "mu", "a", "omega"],
        "properties": {
          "body": {"type": "string", "description": "Perturbing body name, e.g. 'Io'."},
          "mu": {"type": "number", "description": "Nondimensional mass GM_j/(GM_primary_sys+GM_secondary), base-system mass unit (CRNBPPerturber.mu)."},
          "a": {"type": "number", "exclusiveMinimum": 0, "description": "Orbit radius, base-secondary-SMA units (CRNBPPerturber.a)."},
          "omega": {"type": "number", "description": "Synodic angular rate in the base rotating frame, n_j/n2 - 1 (CRNBPPerturber.omega)."},
          "theta0": {"type": ["number", "null"], "description": "Synodic phase at t=0, rad (CRNBPPerturber.theta0)."},
          "rate_provenance": {"type": ["string", "null"], "description": "How omega was set: registry-derived physical rate, resonance-projected (e.g. 'laplace-projected: omega_io = -2*omega_gan exactly'), or idealized -- the honesty axis #724's rate-model novelty claim rests on."}
        }
      }
    },
    "resonance_lock": {
      "type": ["string", "null"],
      "description": "The multi-perturber commensurability making the forcing single-clock periodic, with its phase constraint -- e.g. 'Galilean Laplace lock: omega_io = -2*omega_gan exact (observed-period ratio -2.000000001); phi_L = theta_io0 + 2*theta_gan0 = pi (physical libration center, Sinclair 1975)'. Null if the perturbers are independently periodic (no such row yet)."
    },
    "torus": {
      "type": ["object", "null"],
      "description": "The quasi-periodic torus's own identity and convergence evidence (search.variational_crnbp_torus.CRNBPTorusVariationalResult). The analogue of ccr4bp_provenance's torus_* scalars plus the seed/continuation lineage that IS this object's positive control (#714: no literature-grade N=5 dynamical control exists; the control is continuation from the validated N=4 limit).",
      "additionalProperties": true,
      "properties": {
        "base_resonance": {"type": ["string", "null"], "description": "Base periodic orbit tag, e.g. 'spacecraft:Europa=3:4-exterior (Kumar et al. 2021 class)'. Echoed into orbit_elements.cr3bp.family."},
        "period_tu": {"type": ["number", "null"], "exclusiveMinimum": 0, "description": "Stroboscopic forcing period (theta1 clock, = period_multiple Ganymede-synodic cycles), nondim TU. Echoed into orbit_elements.cr3bp.period_nd."},
        "rotation_number": {"type": ["number", "null"], "description": "Torus rotation number (omega2/omega1 winding, irrational-by-construction)."},
        "rho_strob": {"type": ["number", "null"], "description": "Stroboscopic-map rotation number, matching ccr4bp_provenance.torus_rho_strob's convention."},
        "n1": {"type": ["integer", "null"]}, "n2": {"type": ["integer", "null"]},
        "residual_rms": {"type": ["number", "null"], "minimum": 0, "description": "On-grid pseudospectral collocation residual RMS at convergence."},
        "closure_residual": {"type": ["number", "null"], "minimum": 0, "description": "Independent nonlinear flow-closure residual (search.variational_crnbp_torus._independent_closure: full coupling-term-included crnbp_eom DOP853 flow vs the Fourier surface -- methodologically independent of the algebraic solve). MANDATORY to quote alongside residual_rms per #724 ('off-grid invariance ~2.1e-3 alongside the on-grid 1.23e-4'), since it is the truncation-floor-limited number."},
        "seed_lineage": {"type": ["string", "null"], "description": "Full continuation chain from the validated anchor, e.g. 'Kumar-2021 JE 3:4 exterior resonant PO (CR3BP) -> #690 CCR4BP torus (n1=2,n2=20) -> mu_io=0 CRNBP seed -> 8-step mu_Io continuation to physical 4.7043400305e-05'."},
        "continuation": {
          "type": ["object", "null"],
          "additionalProperties": true,
          "properties": {
            "parameter": {"type": ["string", "null"], "description": "e.g. 'mu_io'"},
            "n_steps": {"type": ["integer", "null"], "minimum": 0},
            "final_value": {"type": ["number", "null"]},
            "monotonicity_note": {"type": ["string", "null"]}
          }
        }
      }
    },
    "real_ephemeris_evidence": { /* see Sec 7 -- torus-point-resolved shadowing evidence */ },
    "method": {
      "type": ["string", "null"],
      "description": "Pipeline identity tag, e.g. 'variational-crnbp-torus(mu_io-continuation)+crnbp-real-ephemeris-consistency (search.variational_crnbp_torus + search.crnbp_real_ephemeris_consistency)'."
    }
  }
}
```

## 7. `real_ephemeris_evidence` for `#729`'s evidence shape

The `ccr4bp_provenance.real_ephemeris_evidence` shape cannot be reused verbatim, for three
structural reasons grounded in what `#726`/`#729` actually measured:

1. **The gap is a model-shadowing mismatch, not connection quality.** There is no
   departure/arrival manifold pair; "target" is the idealized model's own forward flow of a
   chosen torus point over one forcing period, compared against real-SPICE propagation of the
   same departure state (`#726`'s object-type adaptation).
2. **There is a torus-point axis.** `#729`'s headline structural finding — narrow at
   (0,0)/(π,π)/(π/2,π/2), collapse at (π,0)/(0,π), perfectly epoch-stable — is
   two-dimensional evidence. A flat `tested_epochs[]` array (the `#705` shape) would either
   flatten it away or force one row per torus point. The block must be per-torus-point,
   per-epoch.
3. **The recurrence gate is different.** `#705` gated on "comparable to `#704`'s reference
   connection" (a relative 10x threshold). `#729` has no reference connection; its pre-declared
   gate is an absolute `narrow_near_miss_threshold_km = 5000.0`, and the CONTRAST that makes
   "narrow" meaningful is `#726`'s generic-collapse magnitude (3.68e5 km) and noise-floor
   control — both must be recorded or the narrow claim is uninterpretable.

Proposed shape (the per-epoch item keeps `ccr4bp_provenance`'s field names —
`base_epoch_utc`, `local_minimum_epoch_utc`, `local_minimum_pos_gap_km`,
`local_minimum_vel_gap_km_s`, `closest_approach_km{}` — deliberately, so a cross-row consumer's
per-epoch reader works on both blocks; only the verdict boolean differs,
`narrow_near_miss` vs `comparable_to_reference`, because the gates genuinely differ):

```jsonc
"real_ephemeris_evidence": {
  "type": ["object", "null"],
  "description": "Real-ephemeris MODEL-SHADOWING evidence (the #726/#729 findings): idealized-CRNBP forward flow of a chosen torus point over one forcing period vs real-SPICE propagation of the same departure state -- a position/velocity mismatch at the comparison time, NOT a connection-quality metric (no manifold pair exists for this object class). Torus-point-RESOLVED: #729's narrow/collapse dichotomy is torus-point-dependent but epoch-stable, and that 2-D structure is the evidence.",
  "additionalProperties": true,
  "properties": {
    "check_type": {"type": ["string", "null"], "description": "e.g. 'model-shadowing over one Ganymede-synodic forcing period (search.crnbp_real_ephemeris_consistency.check_torus_survives_real_ephemeris)'."},
    "force_model": {"type": ["string", "null"]},
    "forcing_period_days": {"type": ["number", "null"], "exclusiveMinimum": 0, "description": "The real forcing (Ganymede-synodic) period the check propagates over, days."},
    "narrow_near_miss_threshold_km": {"type": ["number", "null"], "minimum": 0, "description": "#729's pre-declared ABSOLUTE narrow-window gate (5000 km) -- unlike ccr4bp_provenance's relative comparable_*_threshold (no reference connection exists here)."},
    "positive_control": {"type": ["object", "null"], "additionalProperties": true, "description": "The checker's own validity gate: {pos_gap_km, vel_gap_km_s, note} for a case where model and reality MUST agree (#726: 1.75e-2 km / 1.13e-7 km/s)."},
    "generic_collapse_reference": {"type": ["object", "null"], "additionalProperties": true, "description": "The off-window magnitude that makes 'narrow' meaningful: #726's headline naive-epoch gap {pos_gap_km, vel_gap_km_s, epoch_utc, noise_floor_note}."},
    "torus_points": {
      "type": ["array", "null"],
      "description": "One entry per tested torus point (theta1_0, theta2_0) -- the #729 second scan axis.",
      "items": {
        "type": "object",
        "additionalProperties": true,
        "properties": {
          "theta1_0": {"type": "number"}, "theta2_0": {"type": "number"},
          "n_synodic_points_per_epoch": {"type": ["integer", "null"]},
          "n_epochs_narrow": {"type": ["integer", "null"]},
          "n_epochs_tested": {"type": ["integer", "null"]},
          "per_epoch": {
            "type": ["array", "null"],
            "items": {
              "type": "object",
              "additionalProperties": true,
              "properties": {
                "base_epoch_utc": {"type": "string"},
                "local_minimum_epoch_utc": {"type": ["string", "null"]},
                "local_minimum_pos_gap_km": {"type": ["number", "null"], "minimum": 0},
                "local_minimum_vel_gap_km_s": {"type": ["number", "null"], "minimum": 0},
                "closest_approach_km": {"type": ["object", "null"], "additionalProperties": {"type": ["number", "null"], "minimum": 0}},
                "narrow_near_miss": {"type": ["boolean", "null"]}
              }
            }
          }
        }
      }
    },
    "epoch_stability_note": {"type": ["string", "null"], "description": "The torus-point-dependent-but-epoch-stable dichotomy statement (#729's key structural finding)."},
    "synodic_duty_cycle_bands_pct": { /* identical shape to ccr4bp_provenance's -- threshold-label keys -> {fraction_of_period, n_distinct_windows}; measured at ONE representative (torus_point, anchor) */ },
    "resolution_caveat": {"type": ["string", "null"], "description": "Mandatory whenever torus_points is non-null (same rule as ccr4bp_provenance)."}
  }
}
```

`validity_window` handling follows `#707` Q4 verbatim: coarse calendar bounds
`{2000-01-01, 2083-01-01}` + `synodic_period_days: 7.053933330456041` (the real Ganymede-synodic
forcing period, `#729`'s own `real_ganymede_synodic_period_days`), the `synodic_duty_cycle_pct`
trio deliberately NOT populated (10 sparse anchors, not a continuous daily scan), with the same
`data_gaps` entry.

## 8. `validation_level` (flagged, not settled — same open question as `#707` Sec 5)

The Umbriel row earned V1 on the ghost guard's independent Radau-vs-DOP853 integrator
cross-check. This object's nearest equivalent is weaker: the torus corrector's independent
closure check (`_independent_closure`) IS methodologically independent of the algebraic solve
(a full nonlinear DOP853 flow through the coupling-term-included `crnbp_eom` vs the spectral
surface — different method class), and `#724` §2 reproduced the whole pipeline bit-for-bit
(4e-11 relative) in the foreground; but **no second-integrator (e.g. Radau) cross-check has
been run for the torus**, and `#726`'s checker also propagates DOP853 only. Recommendation:
**V0 if a strict same-evidence-class reading of the ladder is applied; V1 defensible** if the
coordinating session accepts spectral-solve-vs-nonlinear-flow as the "independent solver
cross-check" spirit — in which case a cheap Radau re-closure at writeback time (minutes) would
make it airtight and is recommended regardless. Genuinely a user decision; the draft row below
carries V0 with the upgrade path documented, erring on the honest side.

`our_status`: left absent, matching the Umbriel row's own explicit precedent ("left absent...
rather than asserted, since no independent verification gauntlet analogous to spec §16.4's
ladder has been run for this object class"). The `#724` narrow-novelty language lives in
`notes` verbatim; `search/literature_check.py` must still be run by the writeback task.

## 9. Worked example — DRAFT row (illustrative only, NOT committed)

All values marked SOURCED trace to `core.crnbp.jupiter_europa_io_ganymede_default()` (printed
this task), `data/found/729_crnbp_epoch_torus_robustness_scan/result.json`, or `#724` §2's
foreground reproduction. Values the writeback task must extract from a fresh
`scripts/verify_724_rerun_continuation.py` run (marked TBD) are the final result object's
`rho_strob` and exact `closure_residual`/`rotation_number` at full precision.

```yaml
- id: europa-3-4-crnbp-torus-jupiter-2026
  name: "Jupiter-Europa 3:4 exterior resonant torus substitute, Laplace-locked N=5 CRNBP (Io+Ganymede forced)"
  source: discovered
  trajectory_regime: ballistic
  model_assumption: crnbp   # PROPOSED new enum value (Sec 5) -- ccr4bp would misstate the model (the Io term is novelty axis 1)
  cycler_class: non-keplerian
  orbit_class: quasi_periodic_torus   # PROPOSED new enum value (Sec 5) -- NOT torus_homoclinic (no connection computed, Sec 2), NOT resonant_po (torus not PO; epoch invariant misfit, Sec 3)
  epoch_locked: true    # Option A (Sec 5): real ephemeris shows GENERIC COLLAPSE off-window (#726); the recurring narrow window IS the evidence
  n_returns: 1          # one forcing period of demonstrated real-ephemeris model-shadowing per recurring narrow window (Sec 5 Option A) -- NOT a tour repeat count
  delta_v_kms:
  v_infinity_leveraging_dv_kms:
  fleet_size:
  flyby_mechanics:      # no flyby of any moon -- see data_gaps
  primary: "Jupiter"
  bodies: ["Jupiter", "Europa", "Io", "Ganymede"]   # Io and Ganymede are gravitational PERTURBERS, never closely encountered (closest approaches ~2.1e5 / ~4.7e5 km, see crnbp_provenance.real_ephemeris_evidence)
  sequence_canonical:   # not-applicable -- see data_gaps (no alternating named-body flyby sequence exists)
  sense: "n/a"
  period:
    pair:
    k:
    years:
    note: >
      No beat-period/encounter-pair concept applies to a bare quasi-periodic torus. The
      stroboscopic forcing period (theta1 clock) is carried in
      crnbp_provenance.torus.period_tu; the motion itself is quasi-periodic (irrational
      rotation number ~0.4965), so no single orbital period exists to state here.
  validity_window:
    start: "2000-01-01T00:00:00Z"
    end: "2083-01-01T00:00:00Z"
    # synodic_duty_cycle_pct / synodic_boundary_period_days deliberately NOT populated
    # (10 sparse anchors, not a continuous daily scan) -- see data_gaps, per #707 Q4's precedent.
    synodic_period_days: 7.053933330456041   # SOURCED: 729 result.json:real_ganymede_synodic_period_days (the real forcing period)
  validation_level: V0   # Sec 8: no second-integrator cross-check run for the torus yet; the
                          # independent nonlinear-flow closure check + #724's bit-for-bit foreground
                          # reproduction are documented in crnbp_provenance.torus; a cheap Radau
                          # re-closure at writeback would support V1 (user decision).
  source_ephemeris: "jup365.bsp (JPL/NAIF Galilean satellite SPK) via core.ephemeris; naif0012.tls"
  orbit_source: derived
  vinf_source: derived
  orbit_fidelity: circular-coplanar   # the idealized-CRNBP torus numbers (crnbp_provenance.torus)
  vinf_fidelity: real-de440   # the #726/#729 real-SPICE shadowing evidence is the row's headline real-world claim (no vinf_kms_at_encounters exists to pair it with -- see data_gaps)
  orbit_elements:
    reference_frame: rotating-synodic
    center: Jupiter
    a_au:
    e:
    note: "Jupiter-Europa rotating-frame quasi-periodic torus; Keplerian elements inapplicable. The cr3bp tuple carries the BASE resonant orbit's identity; the torus's own identity is in crnbp_provenance."
    cr3bp:
      jacobi_constant:   # CRNBP (two time-periodic forcings) has NO exact conserved quantity -- see data_gaps
      period_nd: 12.47771183806834   # SOURCED: 729 result.json:torus_period_tu -- stroboscopic forcing period, nondim TU
      stability_index:   # NOT computed -- no Floquet/stability work exists for this object (#724 mandatory qualifier) -- see data_gaps
      mass_ratio: 2.528017724591319e-05   # SOURCED: core.crnbp.jupiter_europa_io_ganymede_default().mu (Jupiter-Europa reduction)
      libration_point:
      family: "europa-spacecraft-3:4-exterior-resonant-torus (Kumar et al. 2021 seed class; EXTERIOR, never 'interior' -- #724 correction)"
      lunit_km: 671100.0   # SOURCED: 729 result.json:l_km_europa (Europa SMA)
      tunit_s: 48843.87840180734   # DERIVED: l_km / v_unit_europa_km_s (13.739695166696011, 729 result.json); one-line unit conversion, not a golden
  crnbp_provenance:   # PROPOSED schema v5.4 block (Sec 6)
    mu: 2.528017724591319e-05          # SOURCED: jupiter_europa_io_ganymede_default().mu
    base_pair: "Jupiter-Europa"
    perturbers:   # order matches core.crnbp.CRNBPSystem.perturbers in the producing run: (Io, Ganymede)
      - body: "Io"
        mu: 4.70434003054117e-05       # SOURCED: jupiter_europa_io_ganymede_default().perturbers[0].mu
        a: 0.6285203397407242          # SOURCED: same, .a (Io SMA / Europa SMA)
        omega: 1.0071053713566571      # SOURCED: same, .omega
        theta0: 3.141592653589793      # SOURCED: same, .theta0 (pi -- physical Laplace libration center, #723)
        rate_provenance: "laplace-projected: omega_io = -2*omega_gan EXACT (observed sidereal-period ratio -2.000000001, ~1.1e-9 relative residual; registry-rate ratio -1.9996 attributed to SMA rounding -- core.crnbp module docstring)"
      - body: "Ganymede"
        mu: 7.804763238533231e-05      # SOURCED: same, .perturbers[1].mu
        a: 1.59499329459097            # SOURCED: same, .a
        omega: -0.5035526856783286     # SOURCED: same, .omega
        theta0: 0.0
        rate_provenance: "registry two_body_synodic_rate -- PHYSICAL, non-rate-idealized (ephemeris-period value -0.5036473892, 1.9e-4 relative, registry-SMA rounding; both ~7.3e-3 from TCP's idealized -0.5) -- novelty axis 2 (#724)"
    resonance_lock: "Galilean Laplace lock: omega_io = -2*omega_gan exact; phi_L = theta_io0 + 2*theta_gan0 = pi exactly (physical libration center, Sinclair 1975 / Murray & Dermott; TCP Table 1 uses the same phase -- phase is NOT a novelty axis, #724)"
    torus:
      base_resonance: "spacecraft:Europa=3:4-exterior (Kumar et al. 2021, arXiv:2109.14815; a ~= 1.2114 Europa SMA, BETWEEN the moons' orbits)"
      period_tu: 12.47771183806834     # SOURCED: 729 result.json:torus_period_tu
      rotation_number: 0.496468269     # SOURCED: #724 SS2 step-8 rot (full precision TBD from the writeback re-run)
      rho_strob:                       # TBD: extract from the final CRNBPTorusVariationalResult.rho_strob at writeback
      n1: 2
      n2: 20
      residual_rms: 1.2343143649e-04   # SOURCED: #724 SS2 step 8 (matches #723 to 4e-11 relative); 729 result.json:torus_residual_rms
      closure_residual: 2.321e-3       # SOURCED: #724 SS2 -- inside the 5e-3 gate; n2=20 TRUNCATION FLOOR, same as the #690 baseline. MANDATORY to quote alongside residual_rms (~2.1e-3 off-grid invariance vs 1.23e-4 on-grid, #721 SS2 / #724)
      seed_lineage: "Kumar-2021 JE 3:4 exterior resonant PO (CR3BP, perp residual 7.2e-13) -> #690 CCR4BP torus n1=2/n2=20 (residual 1.2210263312e-04) -> mu_io=0 CRNBP seed at theta_io0=pi (1.2209943365e-04, byte-exact N=4 reduction gate) -> 8-step mu_Io continuation to physical 4.7043400305e-05 (self-generated positive control per #714: no literature-grade N=5 dynamical control exists)"
      continuation:
        parameter: "mu_io"
        n_steps: 8
        final_value: 4.7043400305e-05
        monotonicity_note: "residual monotone across steps (one -9.5e-10 re-solve dip at step 1, #723); reproduced bit-for-bit in the foreground by #724 (scripts/verify_724_rerun_continuation.py)"
    real_ephemeris_evidence:
      check_type: "model-shadowing over one Ganymede-synodic forcing period: idealized-CRNBP forward flow (propagate_crnbp) of a chosen torus point vs real-SPICE propagation of the same departure state (search.crnbp_real_ephemeris_consistency.check_torus_survives_real_ephemeris, #726; NOT a connection-quality metric)"
      force_model: "Jupiter point-mass (system-GM) central term + Europa + Io + Ganymede REAL SPICE (jup365.bsp) third-body perturbations, Jupiter-centred J2000 -- verbatim per 729 result.json:force_model_notes / module docstring"
      forcing_period_days: 7.053933330456041
      narrow_near_miss_threshold_km: 5000.0   # SOURCED: 729 result.json:summary.narrow_near_miss_threshold_km (pre-declared ABSOLUTE gate; no reference connection exists for a relative gate)
      positive_control:
        pos_gap_km: 1.75e-2
        vel_gap_km_s: 1.13e-7
        note: "#726's mandatory positive control over an 8.06e5 km / 1.8e5 s arc -- consistent with the accepted ~2.07e-4-relative system-GM-folding approximation"
      generic_collapse_reference:
        epoch_utc: "2030-01-01T00:00:00"
        pos_gap_km: 3.68e5
        vel_gap_km_s: 8.19
        note: "#726 headline at (0,0), naive epoch/phase -- comparable to Ganymede's orbital radius (1.07e6 km); verified genuine (idealized-fed control ~11% of the real gap; 0.02x-window gap 501 km still ~10x the 50 km noise floor)"
      torus_points:
        - theta1_0: 0.0
          theta2_0: 0.0
          n_synodic_points_per_epoch: 300
          n_epochs_narrow: 10
          n_epochs_tested: 10
          per_epoch:   # SOURCED verbatim: 729 result.json:primary_scans[] (bisection-refined local minima) -- ALL 10 anchors, per #708's own no-compression precedent
            - base_epoch_utc: "2000-01-01T00:00:00"
              local_minimum_epoch_utc: "2000-01-07T20:16:37"
              local_minimum_pos_gap_km: 1715.5883486989114
              local_minimum_vel_gap_km_s: 0.03106051633661316
              closest_approach_km: {Europa: 418849.08077241696, Io: 210372.12654338437, Ganymede: 515065.93786844524}
              narrow_near_miss: true
            - base_epoch_utc: "2009-03-23T03:59:59"
              local_minimum_epoch_utc: "2009-03-29T09:33:53"
              local_minimum_pos_gap_km: 655.3534456605724
              local_minimum_vel_gap_km_s: 0.023423867614985116
              narrow_near_miss: true
            # ... remaining 8 anchors (2018-2083): local minima 502.9-1950.4 km, vel gaps
            #     0.0095-0.0343 km/s, all narrow_near_miss: true; full per-epoch data in
            #     data/found/729_crnbp_epoch_torus_robustness_scan/result.json:primary_scans[]
        - theta1_0: 3.141592653589793
          theta2_0: 3.141592653589793
          n_synodic_points_per_epoch: 60
          n_epochs_narrow: 10
          n_epochs_tested: 10
          # per-epoch minima 790.2189007573971 - 1719.2544692912043 km, all narrow (secondary_scans[])
        - theta1_0: 1.5707963267948966
          theta2_0: 1.5707963267948966
          n_synodic_points_per_epoch: 60
          n_epochs_narrow: 10
          n_epochs_tested: 10
          # per-epoch minima 2409.93 - 4042.35 km, all narrow
        - theta1_0: 3.141592653589793
          theta2_0: 0.0
          n_epochs_narrow: 0
          n_epochs_tested: 10
          # GENUINE GENERIC COLLAPSE: minima 25531.6 - 1306778.1 km, 0/10 narrow
        - theta1_0: 0.0
          theta2_0: 3.141592653589793
          n_epochs_narrow: 0
          n_epochs_tested: 10
          # GENUINE GENERIC COLLAPSE: minima 33164.1 - 39984.2 km, 0/10 narrow
      epoch_stability_note: >
        The narrow/collapse dichotomy is torus-point-dependent but PERFECTLY epoch-stable
        (every tested point is 10/10 narrow or 10/10 collapse across 2000-2083) -- itself
        strong evidence of real dynamical structure, not noise (#729's key structural
        finding, extending beyond #705's single-axis scope).
      synodic_duty_cycle_bands_pct:   # SOURCED: 729 result.json:primary_scans[0].duty_cycle_dense_synodic_scan -- ONE representative (torus point (0,0), 2000-01-01 anchor), 300-pt dense scan (pre-bisection; the refined minimum at this anchor is 1715.6 km, below the dense-scan floor of 3793.4 km)
        below_500km: {fraction_of_period: 0.0, n_distinct_windows: 0}
        below_1000km: {fraction_of_period: 0.0, n_distinct_windows: 0}
        below_2000km: {fraction_of_period: 0.0, n_distinct_windows: 0}
        below_5000km: {fraction_of_period: 0.006666666666666667, n_distinct_windows: 2}
        below_10000km: {fraction_of_period: 0.013333333333333334, n_distinct_windows: 4}
        below_20000km: {fraction_of_period: 0.013333333333333334, n_distinct_windows: 4}
        below_50000km: {fraction_of_period: 0.04666666666666667, n_distinct_windows: 5}
      resolution_caveat: >
        10 discrete anchor epochs (~9.2-year spacing, 2000-2083), 300 synodic-phase points
        per primary anchor (60 per secondary combo) + 28-iteration bisection refine -- NOT a
        continuous multi-year daily scan. The duty-cycle bands are from the pre-bisection
        dense scan at ONE representative (torus-point, anchor) pair. Torus-point sampling is
        5 points, not a dense (theta1, theta2) map; the narrow/collapse basin boundary on
        the torus is uncharacterized. A slow secular effect out of phase with all 10 anchors
        could in principle exist undetected.
    method: "variational-crnbp-torus(mu_io-continuation)+crnbp-real-ephemeris-consistency (src/cyclerfinder/search/{variational_crnbp_torus,crnbp_real_ephemeris_consistency}.py; scripts/{verify_724_rerun_continuation,run_729_epoch_torus_robustness_scan}.py)"
  data_gaps:
  - path: "sequence_canonical"
    kind: "not-applicable"
    note: >
      Bare quasi-periodic torus -- no alternating named-body flyby sequence exists. No moon
      is closely encountered (closest approaches at the #729 local minima: Io ~2.1e5 km,
      Europa ~4.2e5 km, Ganymede ~5.2e5 km -- far outside every SOI). Left null rather than
      fabricated, per the #312/#708 rows' precedent.
    todo_ref: "#735"
  - path: "orbit_elements.cr3bp.jacobi_constant"
    kind: "not-applicable"
    note: "The CRNBP (two simultaneous time-periodic forcings) has NO exact conserved quantity; unlike the #708 row there is not even a quasi-Jacobi connection-endpoint gap to report (no connection exists)."
    todo_ref: "#735"
  - path: "orbit_elements.cr3bp.stability_index"
    kind: "unknown"
    note: "No Floquet/stability, manifold, or transfer computation exists for this object (#724 mandatory qualifier; TCP has all three for THEIR Lyapunov-substitute tori -- this row must not imply parity). #714 shortlist item 3 (manifold/connection work) is gated and not dispatched."
    todo_ref: "#735"
  - path: "validity_window.synodic_duty_cycle_pct"
    kind: "unknown"
    note: >
      10 sparse anchors with per-anchor dense synodic scans, not a continuous daily scan
      (same evidence-shape mismatch as the #708 row, per #707 Q4). Within-cycle structure is
      recorded as the threshold ladder in crnbp_provenance.real_ephemeris_evidence.
    todo_ref: "#735"
  - path: "vinf_kms_at_encounters"
    kind: "not-applicable"
    note: "No named-body encounter exists. vinf_fidelity (real-de440) records the fidelity tier of the row's headline real-SPICE shadowing evidence, per the #708 row's convention."
    todo_ref: "#735"
  legs: []
  corroborating_sources:
  - authors: ["Baresi, N.", "Owen, D.", "Scheeres, D. J."]
    year: 2023
    title: "Exploiting the Laplace Resonance for Designing Trajectories in the Jupiter-Io-Europa-Ganymede System"
    venue: "AAS/AIAA Astrodynamics Specialist Conference, Big Sky, MT, AAS 23-201"   # paper number per the PDF's own running header (#722 corrected #721's 'AAS 23-257'); no DOI
    note: "MODEL-CLASS source (mandatory citation, #724): the Tri-Circular/Laplace-locked N=5 idea is theirs; first N=5 tori (L1/L2 Lyapunov-family substitutes), Floquet stability, manifolds. Our object differs on exactly two axes: orbit family substituted (mean-motion-resonant vs Lyapunov) and rate model (physical Ganymede synodic rate vs both rates idealized to exact rationals -- their own 2*pi*k periodicity construction requires the idealization)."
  - authors: ["Owen, D.", "Baresi, N.", "Scheeres, D. J."]
    year: 2024
    title: "Transfer Trajectory Design in the Jupiter-Io-Europa-Ganymede Tri-circular Problem"
    venue: "29th International Symposium on Space Flight Dynamics (ISSFD), Darmstadt, Germany"   # no DOI; acquired via Surrey Open Research (#722)
    note: "MODEL-CLASS source (mandatory citation, #724): planar TCP QPO continuation (70th/75th Lyapunov members) + Europa<->Ganymede transfers; its own SSVI future-work text confirms beyond-Lyapunov substitutes were NOT computed."
  - authors: ["Kumar, B.", "Anderson, R. L.", "de la Llave, R.", "Gunter, B."]
    year: 2021
    title: "Computation and Analysis of Jupiter-Europa and Jupiter-Ganymede Resonant Orbits in the Planar Concentric Circular Restricted 4-Body Problem"
    venue: "AAS/AIAA Astrodynamics Specialist Conference, AAS 21-651; arXiv:2109.14815"
    note: "SEED source: the planar Jupiter-Europa 3:4 EXTERIOR resonant orbit/torus class this object substitutes (their own abstract's 'exterior Jupiter-Europa' terminology -- never 'interior', #724 correction). Their model is N=4 CCR4BP; no N=5/multi-perturber extension attempted (confirmed by #727's full-text read of their 2023 Acta Astronautica follow-on)."
  first_published:
    authors: ["cyclerfinder discovery campaign"]
    year: 2026
    title: "Quasi-periodic invariant-torus substitute of the Jupiter-Europa 3:4 exterior resonant orbit in a Laplace-locked restricted five-body model"
    venue: "cyclerfinder project; task chain #714 -> #717 -> #720 -> #721 -> #722 -> #723 -> #724 -> #726 -> #729 -> #735"
    doi:
  priority_date: "2026-07-27"
  notes: |
    DRAFT ROW -- illustrative only, produced by #735's schema design proposal, NOT a
    writeback. Pending user sign-off on the orbit_class/model_assumption/crnbp_provenance
    schema questions this note raises, AND a fresh search/literature_check.py clearance run
    at writeback time (#724's novelty verdict is conditional on the literature searched to
    date).

    [At writeback, the #724 SS4 headline claim paragraph and its mandatory qualifiers go
    here VERBATIM -- see docs/notes/2026-07-27-724-final-confirmation-n5-torus-novelty.md
    SS4. Never "interior"; never "first N=5 CRNBP torus"; TCP cited as the model class with
    their prior tori/stability/manifolds/transfers; exactly two novelty axes; no
    stability/manifold/transfer computation exists for THIS object.]

    Real-ephemeris standing (#726/#729): the idealized torus does NOT survive real SPICE
    ephemeris generically (headline collapse 3.68e5 km / 8.19 km/s over one forcing
    period), but a narrow near-miss window (refined minima 503-1950 km at torus point
    (0,0), 790-1719 km at (pi,pi)) recurs at ALL 10 tested anchor epochs spanning
    2000-2083; the narrow/collapse dichotomy across the 5 tested torus points is perfectly
    epoch-stable -- evidence of real dynamical structure. our_status left absent per the
    #312/#708 precedent.
```

## 10. Summary of recommendations

| Question | Recommendation | Confidence |
|---|---|---|
| `orbit_class` | New value `quasi_periodic_torus` (NOT `torus_homoclinic` — no connection computed, would overclaim; NOT `resonant_po` — torus-vs-PO encodability + the epoch invariant would misstate `#726`'s generic collapse; NOT `quasi_cycler`/KAM-carve-out — cr3bp-scoped, windowless, characterization-not-discovery) | High on all three rejections; medium on the exact name (`crnbp_torus`/`resonant_torus` listed with tradeoffs) |
| New-class invariant | `epoch_locked=true` / `n_returns=1` ("one forcing period of demonstrated real-ephemeris shadowing per recurring window") — Option B (`false`/`'infinite'`, idealized-object reading) presented as the alternative | Medium — the design's largest judgment call, flagged for the user exactly as `#707` flagged its Q4 |
| `cycler_class` | Existing `non-keplerian`, no new value | High |
| `model_assumption` | New enum value `"crnbp"` (N carried by the perturber-list length, future-proof for N=6+) | High |
| Provenance block | New additive `crnbp_provenance` with an N-body `perturbers[]` LIST mirroring `core.crnbp.CRNBPSystem` 1:1 (option P2), `ccr4bp_provenance` untouched; per-epoch item field names kept compatible with `ccr4bp_provenance`'s for cross-row consumers | High on P2-vs-P1/P3; the sub-block field inventory is open to editorial trimming |
| `real_ephemeris_evidence` shape | Torus-point-RESOLVED (per-point per-epoch arrays), absolute `narrow_near_miss_threshold_km` gate + `generic_collapse_reference` + `positive_control` recorded (the contrast IS the evidence), `epoch_stability_note` for the `#729` dichotomy finding | High — directly shaped by `data/found/729_.../result.json`'s actual structure |
| `validation_level` | V0 (honest strict reading); V1 defensible on the spectral-vs-nonlinear-flow independence argument; a cheap Radau re-closure at writeback recommended either way | Medium — genuine user decision, mirroring `#707` Sec 5 |
| `sequence_canonical` / `period` / `jacobi_constant` / `stability_index` / `vinf_kms_at_encounters` | Null + `data_gaps` entries per the `#312`/`#708` no-fabrication precedent (stability gap doubles as the mandatory "no stability/manifold work exists" `#724` qualifier) | High |

No schema or catalogue file has been modified by this task. All of the above is a proposal for
the coordinating session and user to review, decide, and (if approved) implement as an explicit
schema-bump (v5.4) + writeback task, exactly as `#707`'s design was handled by `#708`.
