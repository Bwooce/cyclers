# #770 — writeback design proposal for `#766`'s C=3.0041 homoclinic self-connection vs `europa-3-4-crnbp-torus-jupiter-2026`

**Date**: 2026-07-31
**Scope**: design proposal only. No `data/catalogue.yaml` or `data/catalogue.schema.json`
change is made by this note. This is the `#707`/`#735` "present before build" step for the
question `#766`'s own closing section explicitly deferred: now that a genuine homoclinic
self-connection exists at the torus row's own cited seed energy (C=3.0041), should the
catalogue row change — `orbit_class`, `validation_level`, and/or new provenance fields?
User-approval-gated per this project's own schema-design convention; nothing below is
implemented. Sources read in full for this design:
`docs/notes/2026-07-29-761-torus-seed-continuation-tractability.md`,
`docs/notes/2026-07-29-766-torus-seed-homoclinic-connection-c30041.md`,
`docs/notes/2026-07-27-735-n5-crnbp-torus-catalogue-schema-design.md` (the original
`quasi_periodic_torus` design), `docs/notes/2026-07-24-707-ccr4bp-catalogue-schema-design.md`
(the sibling `torus_homoclinic` design), `data/catalogue.schema.json` v5.4 in full (the
`orbit_class` enum text, `ccr4bp_provenance.connection{}` field-by-field, `crnbp_provenance`),
the live `europa-3-4-crnbp-torus-jupiter-2026` row in full (including its `#738` V1 basis and
`data_gaps`), `docs/spec.md` §14 (the V0-V5 gauntlet definitions), and
`src/cyclerfinder/search/jovian_resonant_{families,connections}.py`'s sourced constants.

## 1. What the new evidence actually is — and, critically, what object it attaches to

`#761` + `#766` together established, at self-consistency grade:

1. Kumar et al. 2021's "arbitrarily chosen" seed point (C=3.0041, the energy the catalogued
   row's `seed_lineage` cites) is a genuine interior member of the SAME continuous 3:4-LO
   family as the confirmed Anderson-Lo orbit — 24-step continuation, bidirectional
   (reverse re-lands to 2.4e-15), independently corroborated by Kumar's own published
   22,052 km closest approach (reproduced to 0.073%). Endpoint member:
   `x0=-1.3852484456241640`, `ydot0=0.5988394002678391`, `period=25.312119648766764`,
   `C=3.0041` exactly, real saddle `|lambda|=54.589750588953734`.
2. A genuine, transversal homoclinic self-connection (`Wu ∩ Ws` of that orbit) EXISTS at
   C=3.0041: primary hit `(branch_u=+1, branch_s=+1, k_u=6, k_s=6)`,
   `tau_u=10.72913431175392`, `tau_s=14.582984371438714`, crossing
   `(x, xdot)=(-1.4220714951697728, -2.016580027963677e-10)` (on the symmetry axis, the
   same structural type as Anderson & Lo's own published Table-2 point), Newton residual
   `1.97e-10`, ghost distance `0.0368` (37x the `GHOST_GUARD_DELTA=1e-3` guard),
   independent Radau cross-check `2.42e-8` (inside the `<=1e-6` mandate),
   forward/backward re-approach `2.25e-5`/`5.64e-8`, plus an independently-converged
   mirror pair `(5,6)`/`(6,5)` (residuals `3.33e-10`/`6.80e-10`, 67x ghost margins,
   reflection-symmetric to `<1e-6`).
3. Honesty framing carried by `#766` in every place a number appears: **self-consistency
   only** — Kumar 2021 publishes no connection state at this energy, so nothing here is a
   reproduction, and `#754`'s own `(k_u=3,k_s=3)` C_flyby branch does NOT continue to this
   energy (folds at `C≈2.99941`); this is a different `(branch,k)` combination.

**The load-bearing fact for this whole design** (from `#766` §9, unchanged from
`#750`/`#761`): this connection attaches to the **CR3BP seed periodic orbit** — computed in
the **planar CR3BP** at the **papers' own mass ratio** (`ANDERSON_LO_MU =
2.5266448850435e-5`, bit-identical in Anderson & Lo 2011 and Kumar 2021) — **not to the
catalogued N=5 torus**. Three separations, all still open:

- **Model**: the row's object lives in the CRNBP (`model_assumption: crnbp`, Io+Ganymede
  forced); the connection lives in the unforced planar CR3BP.
- **Mass ratio**: the row's `mu = 2.528017724591319e-05` (project DE440-derived); the
  connection's `mu = 2.5266448850435e-5` (the papers'; 0.054% different).
- **Object**: the row's torus's own NUMERICAL seed is still the project's two-body-limit
  proxy (`C=2.9040` at project `mu`, ~409,500 km Europa approach), not Kumar's `C=3.0041`
  orbit; whether the proxy seed's family connects onto this branch remains open (`#750`'s
  attempt broke near `C≈2.926` and nothing since re-answers it).

No manifold of the TORUS has been globalized and no connection of the TORUS has been
attempted — `#724`'s mandatory qualifier ("No stability, manifold, or transfer computation
exists for OUR object") remains literally true. What `#761`/`#766` strengthened is the
**seed lineage**: the cited ancestor orbit is now family-confirmed, saddle-confirmed, and
connection-bearing at its own cited energy.

## 2. Q1 — Does `orbit_class` change to `torus_homoclinic`?

**No — keep `quasi_periodic_torus`. This is the highest-confidence answer in this note,
and the reason is object identity, not evidentiary strength.**

### 2.1 The `#735` condition, quoted exactly, is still true

`#735` §1 (the original design's load-bearing distinction):

> **What this object is NOT** (the load-bearing distinction for the whole design): it is
> not `#708`'s object one system over. The Umbriel-Titania row is a
> torus-HOMOCLINIC-CONNECTION — a computed, ghost-guard-verified manifold intersection
> (unstable departure + stable return to the SAME torus, closed to machine precision) —
> and the `torus_homoclinic` class, its `n_returns=1` semantics, and the
> `ccr4bp_provenance.connection{}` payload are all built around that connection. For the
> N=5 object **no manifold has been globalized and no connection has been attempted**

and `#735` §2's first disqualifier:

> **The defining object is the connection, not the torus.** Every distinguishing property
> in the class description (departs/returns, the two manifold segments, the `n_returns=1`
> "departure-and-return opportunity" semantics) describes a computed manifold
> intersection. None has been computed here.

The schema v5.4 enum text is object-scoped the same way: `quasi_periodic_torus` = "a bare
quasi-periodic invariant torus ... with NO computed manifold connection ... no manifold
has been globalized and no connection attempted **for this object**."

The dispatch framing "that condition is no longer true" holds only under a loose reading
("no connection exists anywhere in this row's lineage"). Under the class's own
object-scoped wording the condition is **still true**: `#766`'s connection is of the seed
periodic orbit, in a different model, at a different mass ratio. The torus's manifolds
remain uncomputed. Reclassifying would assert precisely the thing `#724`'s mandatory
qualifier bars, with evidence about a *different object*.

### 2.2 `torus_homoclinic`'s own schema semantics, field by field

The v5.3 class definition: "a homoclinic connection **to a quasi-periodic invariant
torus** in a **CCR4BP** (departs **the torus's** unstable manifold, returns via **its**
stable manifold to the SAME torus) — epoch_locked=true, n_returns=1 (ONE
departure-and-return opportunity **per recurring real-ephemeris synodic window**...)".
Three structural mismatches before any field is even reached: wrong connected object
(periodic orbit, not torus), wrong model (CR3BP, not a forced 4/5-body model), and the
class invariant presupposes real-epoch-anchored recurrence evidence — `#766` has **zero
real-ephemeris content** (nothing analogous to `#704`/`#705`'s 10-epoch scan was run).

`ccr4bp_provenance.connection{}` field-by-field against `#766`'s actual data:

| Schema field | `#766` referent? |
|---|---|
| `theta2_u` / `theta2_s` (torus second-angle departure/arrival phases) | **None.** The connection is parameterized by `tau_u`/`tau_s` (manifold-leg parameters along a periodic orbit) + crossing indices `k_u`/`k_s` + branch signs — a periodic-orbit manifold parameterization, not torus angles. |
| `t_u_tu` / `t_s_tu` / `*_days` | Exists (`t_u = |t_s| = 134.269` TU ≈ 5.31 periods) — the one clean mapping. |
| `idealized_pos_gap_km` / `idealized_vel_gap_km_s` | Different construction: `#766`'s figure of merit is a nondim section-plane Newton residual (`1.97e-10`), not a full-state km/km-s gap pair at a patch point. Convertible, but not the same measurement. |
| `off_torus_km` (ghost guard: distance from the unperturbed torus) | Different guard: `#766`'s `ghost_distance` (`0.0368` nondim) is distance of the section crossing from the ORBIT's own trivial section point. Same spirit, different object and units. |
| `quasi_jacobi_gap` | **No referent — inverted.** The field exists because the CCR4BP has NO exact conserved quantity; the CR3BP HAS one, and `#766`'s connection conserves the exact Jacobi C=3.0041 by construction. |
| `integrator_delta_km` (independent Radau/DOP853 check) | Analogue exists: the `2.42e-8` Radau cross-check. Comparable evidence class. |
| `residual_norm` | Analogue exists (Newton residual). |

So roughly half the payload has no referent and the block's own description scopes it to a
CCR4BP torus-connection pipeline (`search.ccr4bp_*`), which is not the producing code here
(`search.jovian_resonant_connections`).

### 2.3 The evidentiary-bar sub-question, answered honestly

The dispatch asks whether `torus_homoclinic` presupposes "a literal V1-class
formally-gated connection with an independence gate" that `#766`'s self-consistency-only
result lacks. Checked precisely: **no** — the sibling Umbriel-Titania connection is itself
a novel discovery with no published state gated against; its V1 basis is the ghost guard's
independent-integrator agreement (`integrator_delta_km = 1.05e-7`), i.e. the same
*self-consistency-with-independence-gate* class of evidence `#766` has (Radau `2.42e-8`,
ghost margin 37x, forward/backward re-approach, mirror-pair corroboration). On the
mechanical-verification axis `#766` is genuinely comparable. What it lacks relative to the
sibling is the **real-ephemeris recurrence evidence** (`#704`/`#705`'s 10/10 epochs) that
the class invariant is built around — and, decisively, the right *object*. So even
granting the evidence bar, reclassification fails on identity; and even granting identity,
it would fail on the missing real-ephemeris axis. The rejection is overdetermined.

### 2.4 If a torus connection is ever computed, the path is a new row, not reclassification

`#735` §2 already answered the hypothetical future: "If/when `#714` item 3 (N=5
manifold/connection work) is eventually run AND succeeds, the natural outcome is a
*separate* connection row (or a supersedes/superseded_by pair) — not a retroactive
blurring of this class." That remains the right mechanism, and `#766`'s result is not that
event (it is CR3BP seed-lineage work, not N=5 manifold work).

## 3. Q2 — Does `validation_level` change? (spec §14 walked explicitly)

**No — V1 stands, unchanged.** Two independent reasons, either sufficient:

**(i) Wrong object.** `validation_level` gates the ROW's object (the N=5 torus).
Spec §14's V2-ballistic explicitly evaluates "in the row's defining model"; `#766`'s
evidence is about a different object in a different model at a different `mu`. It cannot
promote (or demote) this row's level any more than the sibling row's evidence could.

**(ii) Even in principle, a homoclinic SELF-connection does not reach V2+.** Walking the
gates for a hypothetical future where an equivalent connection existed for the torus
itself, in the row's own model:

- **V2 (≥3 continuous laps, bounded rotating-frame drift)**: a homoclinic point defines a
  single bi-asymptotic excursion — depart, flow ~5.3 periods, re-approach asymptotically.
  There are no "laps"; the drift-over-laps instrument has no referent. The transversality
  `#766` established does imply (Smale-Birkhoff) a horseshoe and hence infinitely many
  nearby periodic orbits that shadow multi-excursion itineraries — but none has been
  computed, and any such multi-excursion periodic orbit would be a NEW object (exactly the
  Anderson-Lo 3:4↔5:6 / Vaquero "resonant chain" construction, cf. `#768`), earning its
  own row and its own gauntlet, not a promotion of this one. This mirrors the sibling
  row's own standing note ("V2/V3 have no literal analogue for a one-shot
  torus-homoclinic transfer") and the bare-torus finding of `#735` §8 — the homoclinic
  structure changes the *in-principle reachability of chained laps* (via the horseshoe)
  but does not itself constitute them.
- **V3 (ephemeris realisation, budget-bounded TCM)**: no real-ephemeris work of any kind
  exists for the connection — not even the sparse-anchor recurrence scan the sibling has.
  Fails at the threshold.
- **V4 (independent codebase + ephemeris)**: nothing run.
- **V5 (novelty + expert review)**: orthogonal to mechanical validation; and note the
  connection itself should NOT be spun as novelty-bearing — Anderson & Lo 2011 publish
  homoclinic connections of this same family (at other energies); a connection at
  C=3.0041 specifically is unpublished but is an extension along a published family, and
  `#766` deliberately never claims otherwise. The project's novel-findings register
  (three entries) is unchanged by this result.

Conclusion: the row's V1 (per `#738`'s Radau cross-check of the torus's own closure claim)
is the correct, unchanged level. The new evidence is provenance, not promotion.

## 4. Q3 — What additive fields would carry `#766`'s connection data?

Three options considered:

- **P1 — reuse/extend `ccr4bp_provenance.connection{}`** (populate it on this row, or add
  the missing fields to it). Rejected: that block is scoped by its own description to
  CCR4BP `torus_homoclinic` rows and the `#689`-`#708` pipeline; half its fields have no
  referent here (§2.2); and populating a `connection{}` on a `quasi_periodic_torus` row
  would blur exactly the class boundary the v5.4 enum text draws ("its connection{}
  payload is populated for that class and has no referent here").
- **P2 (recommended) — a new OPTIONAL sub-block under `crnbp_provenance.torus`, named
  `seed_orbit_homoclinic`**, following the `#738` `radau_cross_check` precedent exactly:
  `#738` added its evidence block row-side under `crnbp_provenance.torus` WITHOUT a schema
  bump (legal — `crnbp_provenance` and `torus{}` both carry `additionalProperties: true`;
  verified: `radau_cross_check` appears nowhere in `catalogue.schema.json`). The block
  lives next to `seed_lineage`, which is exactly the thing it evidences. A later v5.1-style
  formalization pass (naming the fields in the schema after the shape has settled) remains
  available and is the project's own established pattern for ad-hoc-then-formalized fields.
- **P3 — formalize as schema v5.5 now** (same block, but with named/typed schema fields in
  this same change). Defensible — but heavier than the evidence requires: unlike
  `#707`/`#735` this introduces no new class, no new enum value, no new top-level block,
  and only one row will ever carry it in its current shape. Presented as the alternative
  if the user prefers schema-named fields from day one.

Naming, honestly weighed:

| Candidate | For | Against |
|---|---|---|
| `seed_orbit_homoclinic` (recommended) | Names the object the evidence attaches to (the SEED ORBIT, not the torus) — the single most important semantic fact; sits naturally beside `seed_lineage` | Slightly long |
| `seed_lineage_connection` | Ties to the existing `seed_lineage` field name | "connection" unqualified invites confusion with `ccr4bp_provenance.connection{}`'s torus-connection semantics |
| `cr3bp_seed_homoclinic` | Model visible in the name | Encodes model identity into a field name (the `#735` naming principle argues model identity belongs in a `model:` field inside the block, as proposed below) |

### Proposed row-side block (DRAFT — the exact writeback shape; values from `#761`/`#766`'s notes, full-precision re-extraction from the producing functions at writeback time per the `#735` TBD convention)

```yaml
  crnbp_provenance:
    torus:
      # ... existing fields unchanged (seed_lineage, radau_cross_check, ...) ...
      seed_orbit_homoclinic:   # #770 design / #77x writeback: SUPPLEMENTARY, NON-GATING
        # seed-lineage evidence. This is a connection of the SEED ORBIT in ITS OWN base
        # model -- NOT a connection of this row's N=5 torus, whose manifolds remain
        # uncomputed (#724 mandatory qualifier unchanged; orbit_class stays
        # quasi_periodic_torus per #770's design note). Not a promotion gate
        # (validation_level V1 basis is torus.radau_cross_check, unchanged).
        object: >
          Planar-CR3BP homoclinic self-connection (Wu ∩ Ws) of the 3:4-LO family member
          at the seed's own cited energy C=3.0041 -- Kumar et al. 2021's "arbitrarily
          chosen" seed point, confirmed by #761 to be a genuine interior saddle on the
          SAME continuous family as the confirmed Anderson-Lo 3:4-LO.
        model: "planar CR3BP, mu = 2.5266448850435e-5 (ANDERSON_LO_MU -- the papers' own bit-identical value; NOT this row's project mu 2.528017724591319e-05, 0.054% apart)"
        orbit:   # SOURCED: continue_34lo_to_kumar_c() endpoint, re-derived fresh by #766
          x0: -1.3852484456241640          # run-to-run float noise vs #761's -1.3852484456241585 documented in #766 Sec 1
          ydot0: 0.5988394002678391
          period_nd: 25.312119648766764
          jacobi_constant: 3.0041          # EXACT (the CR3BP has a conserved quantity, unlike this row's CRNBP)
          lambda_max: 54.589750588953734   # real saddle; ~19x weaker instability than the C_flyby member's 1036
          family_membership: "#761: 24-step continuation from the confirmed C_flyby member, bidirectional (reverse re-lands to 2.4e-15); closest approach 22,035.8 km vs Kumar 2021's published 22,052 km (0.073%)"
        connection:   # SOURCED: #766 primary hit (find_homoclinic, rank_by_residual=True)
          branch_u: 1
          branch_s: 1
          k_u: 6
          k_s: 6
          tau_u: 10.72913431175392
          tau_s: 14.582984371438714
          t_u_tu: 134.269                  # TBD full precision at writeback (~5.31 orbital periods; |t_s| equal)
          crossing_x: -1.4220714951697728
          crossing_xdot: -2.016580027963677e-10   # on the symmetry axis -- same structural type as Anderson & Lo's own Table-2 point
          newton_residual: 1.9723783076056544e-10
          ghost_distance: 0.03682304954560878     # nondim section-plane distance from the orbit's own trivial point
          ghost_guard_delta: 1.0e-3               # 37x margin -- real, not razor-thin (#766 Sec 4 discusses why smaller than #754's 145x)
          radau_crosscheck: 2.42e-8               # TBD full precision; DOP853-vs-Radau via crosscheck_cycle, <=1e-6 mandate
          backward_reapproach: 5.638e-8           # homoclinic_reapproach_check: re-traces the unstable leg to its epsilon seed
          forward_reapproach: 2.252e-5            # and the stable leg forward (roundoff amplified ~|lambda|^5.3; see #766 Sec 5)
          mirror_pair_note: "(k_u,k_s)=(5,6)/(6,5) independently converged (residuals 3.33e-10/6.80e-10, ghost 0.0670 = 67x each, reflection-symmetric to <1e-6) -- transversality at this energy is corroborated, not an isolated index fluke"
        evidence_class: >
          SELF-CONSISTENCY ONLY: no published state exists at this energy to gate against
          (Kumar 2021 reports no connection for its own seed orbit) -- never a
          reproduction claim, and NOT novelty-bearing (Anderson & Lo 2011 publish
          homoclinic connections of this same family at other energies; this extends a
          published family to an unpublished energy).
        caveats: >
          #754's own (k_u=3,k_s=3) C_flyby connection branch does NOT continue to this
          energy (fold/tangency at C~2.99941; #766 Sec 3) -- this is a genuinely
          different (branch,k) combination. The torus's own NUMERICAL seed remains the
          project's two-body-limit proxy (C=2.9040 at project mu); whether the proxy
          seed's family connects onto this branch is still open (#750/#761/#766's
          unchanged honest gap) -- this block certifies the CITED lineage's energy, not
          the numerical seed's.
        method: "search.jovian_resonant_families.continue_34lo_to_kumar_c + search.jovian_resonant_connections.{build_34lo_kumar_c_node, find_homoclinic(rank_by_residual=True), homoclinic_reapproach_check}; tests/search/test_jovian_resonant_connections.py (27/27, not slow-marked); docs/notes/2026-07-29-76{1,6}-*.md"
```

### Companion row edits at writeback (no schema change, no reclassification)

- **`notes`**: append a short "Seed-lineage standing (#761/#766)" paragraph mirroring the
  block (family-confirmed + saddle + self-consistency-grade homoclinic at C=3.0041 in the
  seed's own CR3BP; torus manifolds still uncomputed; not a promotion, not novelty-bearing),
  and extend the "Discovery + verification chain" line with `-> #750 -> #761 -> #766 ->
  #770 (design) -> #77x (writeback)`.
- **`data_gaps` `orbit_elements.cr3bp.stability_index`**: keep `kind: unknown` and the
  `#724` qualifier verbatim (still true for the torus), appending one sentence noting the
  seed orbit's own CR3BP saddle character is now established (`#761`), so the gap is
  specifically the N=5 object's Floquet/manifold work, not total ignorance of the lineage.
- **`data_gaps` `orbit_elements.cr3bp.jacobi_constant`**: its "not even a quasi-Jacobi
  connection-endpoint gap to report (no connection exists for this object)" note stays
  literally true (object-scoped); optionally add a pointer to
  `torus.seed_orbit_homoclinic.orbit.jacobi_constant` for the seed's own exact C.
- **Nothing else changes**: `orbit_class`, `epoch_locked`, `n_returns`,
  `validation_level`, `model_assumption`, census counts, and every ratchet-frozen field
  are untouched. All catalogue ratchets (`uv run pytest tests/data tests/search -q`, per
  `[[feedback_catalogue_edits_run_all_ratchets]]`) plus full mypy must run at writeback.

## 5. Q4 — Recommendation

**Option (b): keep `orbit_class: quasi_periodic_torus`, add the OPTIONAL
`crnbp_provenance.torus.seed_orbit_homoclinic` block (row-side, `#738`-precedent, no
schema bump) documenting `#766`'s connection as supplementary, non-gating seed-lineage
evidence. `validation_level` stays V1.**

Options considered:

- **(a) Reclassify to `torus_homoclinic` + connection data — REJECTED.** Overclaims on
  object identity (the connection is of the CR3BP seed periodic orbit at the papers' `mu`,
  not of the catalogued N=5 torus), on model, and on the class's real-ephemeris-anchored
  invariant (no epoch evidence exists for this connection at all). It would also flatly
  contradict `#724`'s mandatory qualifier, which remains true. Notably this is NOT because
  `#766`'s evidence is weak — its independence-gated self-consistency package is of the
  same class as the sibling row's own V1 basis (§2.3); it is because the evidence is about
  a different object. Tradeoff acknowledged per the dispatch: rejecting (a) is not
  underclaiming, because (b) still records the full result where it belongs.
- **(b) RECOMMENDED — see above.** Tradeoffs: the row grows another bespoke sub-block
  (mitigated: it sits under the already-bespoke `torus{}` next to `seed_lineage` and
  `radau_cross_check`, and is self-documenting); a schema-side consumer won't find named
  fields for it until a later v5.1-style formalization (mitigated: same is already true of
  `radau_cross_check`, and `additionalProperties: true` makes it legal; P3 offered if the
  user prefers formalizing now).
- **(c1) Do nothing to the row — REJECTED as underclaiming.** The row's `seed_lineage` is
  the catalogue's own citation of Kumar's seed point; `#761`/`#766` materially upgrade
  that lineage (family-confirmed, saddle-confirmed, connection-bearing at the cited
  energy) with committed, CI-running evidence. The catalogue is the project's record of
  record; leaving this only in `docs/notes` recreates the docs/code divergence
  `[[feedback_digest_not_adoption]]` warns about.
- **(c2) A separate catalogue row for the connection itself — REJECTED for now.** The
  connection is a planar-CR3BP object in a mined envelope, at best an unpublished-energy
  member of a published class (Anderson & Lo's own Jupiter-Europa 3:4 homoclinic
  connections) — it would enter, if at all, as a known-class-adjacent corroboration row
  needing its own `literature_check.py` gate, and nothing forces that. `#735` §2's
  "separate connection row" path stays reserved for its intended trigger: an actual N=5
  torus-manifold connection (`#714` item 3), which this is not.

## 6. Summary of recommendations

| Question | Recommendation | Confidence |
|---|---|---|
| Q1 `orbit_class` | UNCHANGED `quasi_periodic_torus` — `#735`'s object-scoped "no connection for THIS object" condition is still true; `torus_homoclinic` fails on connected-object (PO vs torus), model (CR3BP vs CCR4BP-class), missing real-ephemeris invariant basis, and ~half of `connection{}`'s fields having no referent | High |
| Q2 `validation_level` | UNCHANGED V1 (`#738` basis) — wrong object for promotion; and even in principle a homoclinic SELF-connection reaches none of V2/V3/V4/V5 without computing a new (chained/multi-excursion) object that would be its own row | High |
| Q3 new fields | OPTIONAL `crnbp_provenance.torus.seed_orbit_homoclinic{}` per §4's draft — row-side, `#738` `radau_cross_check` precedent, NO schema bump (P3 schema-v5.5 formalization offered as the alternative); plus the §4 companion `notes`/`data_gaps` touch-ups | High on content; medium on P2-vs-P3 packaging and the exact block name (alternatives listed) |
| Q4 overall | **Option (b)** — record without reclassifying; (a) overclaims object identity, (c1) underclaims a genuine committed result, (c2) has no trigger yet | High |

No schema or catalogue file has been modified by this task. All of the above is a proposal
for the coordinating session and USER to review, decide, and (if approved) implement as a
small writeback task — exactly as `#707`'s design was handled by `#708` and `#735`'s by
`#736`.
