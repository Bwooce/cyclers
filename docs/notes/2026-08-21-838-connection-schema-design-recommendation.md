# #838 — cross-row manifold-connection evidence: schema design recommendation

**Date**: 2026-08-21
**Scope**: design proposal only. No `data/catalogue.yaml`, `data/catalogue.schema.json`, or any
other data/schema file is modified by this note. This follows the `#707` -> v5.3 and `#735` ->
v5.4 process precedent exactly (both of those schema bumps were made by a *separate*,
user-approved implementation task — `#708` and `#736` — after a design note like this one, never
by the task that first wanted the change). Sources read in full for this design: `#838`'s own
registration in `data/OUTSTANDING.md`; `docs/notes/2026-08-11-822-vaquero-em-free-transfer.md`;
`docs/notes/2026-08-12-828-vaquero-connection-tier-adjudication.md`;
`docs/notes/2026-08-21-827-kumar-table5-reproduction.md` (its `#854` section);
`docs/notes/2026-07-24-707-ccr4bp-catalogue-schema-design.md` and
`docs/notes/2026-07-27-735-n5-crnbp-torus-catalogue-schema-design.md` (process templates); the
live `ccr4bp_provenance.connection` block (`umbriel-1-2-torus-homoclinic-uranus-2026`) and
`crnbp_provenance.torus.seed_orbit_homoclinic.connection` block
(`europa-3-4-crnbp-torus-jupiter-2026`) in `data/catalogue.yaml`; the two `#828` `ADDED
EVIDENCE` comment annotations (catalogue lines ~58044 and ~58146); `data/cycler_network.schema.json`
(v1.0, `#570`) + `src/cyclerfinder/data/validate_networks.py`;
`scripts/migrate_catalogue_scope_2026-06-15.py` (migration-script pattern); and
`data/found/822_vaquero_em_free_transfer/results.json` +
`src/cyclerfinder/search/vaquero_em_cycler_connections.py::OVERLAP_GRID_ICS` (the actual numbers
the worked example in §5 uses — nothing is re-derived, every value below is copied from an
existing artifact).

---

## 1. The question, restated precisely

Should the catalogue gain a structural way to record a **verified manifold connection between two
orbits where at least one endpoint is a catalogued row and the connection is not intrinsic to
that row's own object**?

The two existing connection payloads do not cover this case, by their own design:

- `ccr4bp_provenance.connection` (schema v5.3, `#707`/`#708`): the row *is* the connection — a
  torus-homoclinic object whose `orbit_class` (`torus_homoclinic`) presupposes it.
- `crnbp_provenance.torus.seed_orbit_homoclinic.connection` (added under `#770`/`#779` inside the
  v5.4 block): a homoclinic of the row's own *seed orbit*, in the seed's own base model —
  explicitly "SUPPLEMENTARY, NON-GATING seed-lineage evidence" of that one row.

Both are annotated "Provenance/audit only — not a promotion gate", and both are **intrinsic**: the
connection belongs to the single row that carries it.

The motivating case is different in kind, and — per `#838`'s registration, which is explicit on
this — it is **half-catalogued**, not symmetric: `#822`'s two row-touching heteroclinics each have
exactly ONE endpoint that is a catalogue row (`vaquero-21-c266-em-cycler-2013` as Wu origin at
C=2.66; `vaquero-31-c254-em-cycler-2013` as Ws destination at C=2.54), while the counterpart node
in each (3:1@C=2.66, 2:1@C=2.54) is an uncatalogued family member. No case with BOTH endpoints
catalogued exists in the corpus today. The design must make the half-catalogued connection
expressible first, and must not assume the two-row case.

`#828` handled `#822`'s two instances deliberately WITHOUT a schema change: comment-only `ADDED
EVIDENCE (…, NO LEVEL CHANGE)` annotations on the two rows' `validation_level` comments
(`yaml.safe_load` pre/post identical), and registered `#838` for this design instead. `#827`/`#854`
has since surfaced a **third** instance on one of the same rows: a Wu(3:1)->Ws(2:1) connection at
C=2.54, digit-grade matched (7.6e-07) to a state Kumar et al. (2026) themselves published, at
Kumar's own printed mu — a different evidence kind on the same node, not covered by `#828`'s
existing annotation.

## 2. Recommendation (read this first)

**Yes to a structural home, but NOT as a catalogue field and NOT as a catalogue row: adopt a new
dedicated top-level registry — `data/manifold_connections.yaml` validated by a new
`data/manifold_connection.schema.json` (schema v1.0, its own version line, like
`data/cycler_network.schema.json` v1.0) — with a referential-integrity semantic gate mirroring
`validate_networks.py`, and admission requiring at least one endpoint to resolve to a catalogue
row id.** Endpoints are an ordered pair of tagged objects, each either a `row_ref` (with mandatory
identity evidence) or an inline `uncatalogued` descriptor, so the half-catalogued case is the
base case, not a degenerate one. `data/catalogue.yaml` and `data/catalogue.schema.json` are
untouched in this slice — **no v5.5 bump and no migration script are required** — and the
existing rows keep only short comment pointers into the registry.

## 3. Reasoning, point by point

### 3.1 Is a structural field/object needed at all, or is prose still enough?

Prose was the right call for `#828` (two instances, no home, design not yet approved), but the
evidence that it is already straining is concrete:

- The two `#828` annotations are each a **single-line YAML comment of ~2,500+ characters**
  (catalogue lines ~58044/~58146) carrying full numeric payloads as free text — `k=(40,27)`,
  `branches=(-1,+1)`, Newton residual `3.337e-10`, 4-state gap `4.484e-06`, Radau gap,
  ghost margins, transit times, node-identity deltas to `3.1e-15`. None of it is visible to
  `yaml.safe_load`, none of it is queryable, and none of it is cross-checkable by any ratchet.
- The pattern has recurred **three times in ten days** on this one family pair alone, and two
  more instances are already registered: `#854` needs a *second* annotation on the *same*
  C=2.54 row (a different evidence kind — Kumar's published digits at Kumar's mu — which
  `#828`'s annotation explicitly does not cover); `#840` would *amend both* existing annotations
  (one-way -> round-trip at each band edge); `#839`, if it lands a heteroclinic at C=3.13,
  produces a fourth connection touching a third row. Every one of these, under prose-only,
  means re-editing giant comment strings that no parser sees.
- The strongest forward-looking argument is `#839` itself: its registered question is whether
  "demonstrated transport utility" flips `vaquero-31-c313`'s `orbit_class` from `resonant_po` to
  `cycler`. "Which catalogued orbits have verified transport connections, of what grade?" is
  *exactly* a query over connection records. Free-text comments cannot answer it; a registry can.
- Finally, the **uncatalogued endpoint has no structured home anywhere today**. The 2:1@C=2.54
  node's identity `(x0=1.0905…, ydot0=-0.8232…, T=5.9412…, |lambda|=3.1498)` lives only in a
  Python constant (`OVERLAP_GRID_ICS`) and `results.json`. Half of every half-catalogued
  connection is currently invisible to the data layer entirely.

So: structural, yes. The comment convention stays available for genuinely one-off marginalia, but
cross-row connection evidence is now a recurring, growing, query-relevant class.

### 3.2 Field on endpoint row(s), own catalogue row, or other shape?

**Neither of the first two. A separate registry file.**

*Against a field on the endpoint row(s):*
- **Ownership/mirroring is unsolvable cleanly.** With one catalogued endpoint the field's host row
  is forced (fine), but the moment a second endpoint is catalogued (a live possibility — `#839`'s
  C=3.13 counterpart 2:1 node could be admitted someday, or `#840`+future admissions close the
  loop) you must either duplicate the payload on both rows (drift risk between two hand-kept
  copies of the same evidence, in a 60k-line file) or pick one arbitrarily (the other row's reader
  never sees it). `#838`'s registration warns against designing only for the symmetric case; a
  row-field design is the shape that *degrades* when symmetry arrives.
- **It blurs the intrinsic/extrinsic line the schema already draws.** `#828`'s adjudication leans
  on exactly this distinction ("both existing payloads are INTRINSIC to their own row's object").
  A third per-row connection block that is *not* intrinsic would make the convention "a
  connection block on a row describes the row's own object" false, silently, for every future
  consumer.
- Catalogue rows are already enormous; per `#854` and `#840` a single row can accumulate
  *multiple* connection-evidence records. A list-of-connections field on rows reinvents the
  registry, worse, inside the biggest file in the repo.

*Against its own catalogue row:*
- A connection is **not an orbit**, and every obligation a catalogue row carries is defined for
  one object's own trajectory: `orbit_class`/`cycler_class`, `sequence_canonical`,
  `validation_level` (per `#828` Sec. 2, spec §14's ladder measures "one object's own trajectory
  under increasing model fidelity" — a connection "moves sideways on the ladder, not up", so a
  connection row could never honestly fill the field), the census counts, and every frozen-census
  ratchet (`[[feedback_catalogue_edits_run_all_ratchets]]`). A connection pseudo-row would be
  null/not-applicable on nearly every axis and would contaminate the census the ratchets freeze.
- The v5.3 `torus_homoclinic` precedent is not a counterexample: that row is admitted as an
  *orbit-class object* (a torus with its own real-ephemeris window evidence) whose identity
  happens to be a homoclinic structure — the row is still one object.

*For a dedicated registry:* this is not a new invention — it is the repo's **third** use of an
established pattern. `data/cycler_network.schema.json` (v1.0, `#570`) states the principle in its
own description: a cross-row relation is "a NEW top-level relation, separate from
data/catalogue.yaml (matching the data/empty_regions.jsonl precedent of a dedicated registry
rather than overloading the per-orbit-row schema)". A network (a *set* of rows sharing cadence)
and a connection (a *directed pair* of orbits, at least one a row) are the same species of thing:
extrinsic relations over catalogue rows. The registry pattern also comes with a proven
enforcement shape: `validate_networks.py`'s three-stage gate (JSON-Schema structural / semantic
cross-field / referential — every `row_ref` must resolve to a live catalogue id) transfers
directly as a new `validate_connections.py`.

Discoverability from the row side is preserved the cheap way: the touched rows keep (or gain) a
**short** comment pointer — "cross-row connection evidence: see data/manifold_connections.yaml id
`<id>`" — which is comment-only and parse-identical, the same zero-risk edit class `#828` used.

### 3.3 The half-catalogued case, concretely

Each entry carries an ordered `endpoints` pair (index 0 = unstable-manifold origin, index 1 =
stable-manifold destination; entries are **directed**). Each endpoint is one of:

- `row_ref: <catalogue id>` — plus a **mandatory** `identity_evidence` string stating how the
  connection's node was matched to the row's recorded numbers (for `#822`/`#828` this is real
  data: re-converged node reproduces `state_nd` to 3.1e-15 / 4.4e-15, `period_nd` to ~1e-13,
  etc.), and an optional `model_note` for cases where the entry's model constants differ from
  the row's recorded ones (needed for `#854`: Kumar's printed mu `1.2150584270572e-2` vs the
  row's recorded `1.215058439469525e-2`, ~1.0e-8 relative — see §6).
- `uncatalogued: {…}` — an inline descriptor with the node's own identity: family label, `x0`,
  `ydot0` (or full `state_nd`), `period_nd`, `jacobi_constant`, stability, and a `derivation`
  provenance string. This is what makes the half-catalogued case first-class: when there is no
  second row to point to, the second endpoint is *described in place*, from numbers that already
  exist in the producing artifact — never a dangling reference, never a fabricated row.

Semantic-gate rule: **at least one endpoint must be a `row_ref`** (that is the admission
criterion — it is what makes an entry catalogue-adjacent). Connections with zero catalogued
endpoints (e.g. `#822`'s other 11 grid points, the `#754`/`#759`/`#767`/`#781` systems with no
rows at all) stay where they live today: `data/found/` + `docs/notes/`. The registry is an
evidence *index over the catalogue*, not a second archive of all connection research.

### 3.4 Consistency with the two intrinsic connection payloads

There is genuine tension here and it should be named rather than papered over: after this
proposal the file set contains **three** places a `connection:` mapping can appear. Two
mitigations make that acceptable, and one alternative makes it clearly the lesser evil:

- **Same vocabulary, deliberately.** The registry's `connection`/`evidence` sub-blocks reuse the
  field names the intrinsic payloads and the producing code already use (`branch_u/branch_s`,
  `k_u/k_s`, `tau_u/tau_s`, `newton_residual`, `ghost_distance_*`, `radau_*`,
  `backward/forward_reapproach`, `t_u/t_s`, `epsilon`) — the same names
  `seed_orbit_homoclinic.connection` uses and `results.json` records. A consumer that can read
  one can read all three.
- **Same governing sentence, verbatim.** The registry schema's description carries the exact
  clause both intrinsic blocks already carry: *"Provenance/audit only — not a promotion gate
  (validation_level carries that)."* `#828` Sec. 2's ruling (a connection is a two-object,
  same-fidelity transport statement; the §14 ladder measures one object under increasing
  fidelity) is cited in the schema text so no future adjudicator re-litigates it from scratch.
- **The scope rule is one sentence:** *intrinsic to one row's own object -> inside that row's
  provenance block; extrinsic (between/touching rows) -> the registry.* The alternative — one
  unified home, i.e. migrating the two intrinsic payloads out into the registry — is strictly
  worse: the v5.3 row's `orbit_class: torus_homoclinic` presupposes its connection (the row IS
  the connection), the v5.4 seed-homoclinic is seed-lineage evidence meaningless apart from its
  row, and moving either would be a breaking rework of two user-approved schema bumps for zero
  representational gain.

So: three payload sites, one shared vocabulary, one explicit scope rule — not three incompatible
representations.

### 3.5 Migration and retroactive scope

- **No catalogue migration.** The registry is a new file with its own v1.0 schema (the
  `cycler_network.schema.json` versioning convention — registries version independently of the
  catalogue's v4.x/v5.x line). `data/catalogue.schema.json` is not touched, so there is **no
  v5.5 bump** and no `migrate_catalogue_scope_*`-style line-by-line script (that pattern exists
  for row-field backfills; nothing here adds a row field). If the user *also* wants a
  machine-readable back-pointer field on rows (`connection_refs: [<registry ids>]`), that IS an
  additive v5.5 catalogue bump — nullable, no migration script needed since only touched rows
  gain it, plus a symmetric-consistency check in the semantic gate. I recommend **deferring**
  that: with ~3 touched rows, comment pointers suffice, and the field can be added later without
  reshaping the registry.
- **Retroactive: yes, for `#822`'s two row-touching connections — populate at adoption.** Every
  number already exists in `data/found/822_vaquero_em_free_transfer/results.json`,
  `OVERLAP_GRID_ICS`, and the `#828` note; population is transcription plus the referential
  gate, zero recomputation (a hard constraint of this design, honored). This also converts the
  two giant `#828` comments into candidates for shrinking to short pointers (comment-only,
  parse-identical edit — optional, user's call, see §6).
- **`#854`'s finding: entry drafted but gated on `#854`'s own adjudication.** `#854` is a
  registered, undispached adjudication; this design should not pre-empt its verdict. What the
  design *does* change is the shape of `#854`'s deliverable: if adopted, `#854` writes a registry
  entry (plus a one-line comment pointer) instead of a second 2,500-character comment on the
  same row. Likewise `#840`'s future band-edge reverse demonstrations become new directed
  entries (or a populated `reverse_of` cross-link), not comment amendments.

## 4. Proposed registry schema, in outline

`data/manifold_connection.schema.json`, v1.0 (task `#838`), validating
`data/manifold_connections.yaml` (a YAML list, shipping with the two `#822` entries below).
Top-level description carries: the intrinsic-vs-extrinsic scope rule (§3.4), the "provenance/audit
only — not a promotion gate" sentence verbatim, and the ≥1-`row_ref` admission criterion.
Per-entry required fields: `id`, `kind` (`heteroclinic`|`homoclinic`), `model` (`type` enum
`cr3bp`|`ccr4bp`|`crnbp` + `system` + `mass_ratio`), `endpoints` (exactly 2, ordered
unstable-origin -> stable-destination, each `row_ref`+`identity_evidence` or `uncatalogued{…}`),
`connection` (geometry: branches, k's, taus, crossing, residual, epsilon), `evidence` (the
verification battery), `evidence_class` (free text: self-consistency vs digit-grade-reproduction
vs published-state-matched — the `#822`-vs-`#827` distinction `#854` turns on), `provenance`
(task refs, data path, module, commit, notes). Optional: `reverse_of`/`round_trip_note` (the
`#828` asymmetry honesty, `#840`'s hook), `dv_kms` (0.0 for a true heteroclinic).
Gate: `src/cyclerfinder/data/validate_connections.py`, three stages mirroring
`validate_networks.py` (structural / semantic: ordered-endpoint roles, ≥1 row_ref,
identity_evidence mandatory with row_ref / referential: every `row_ref` resolves to a live
catalogue id), wired into `tests/data/`.

## 5. Worked example — `#822`'s two real row-touching connections (DRAFT, not committed)

All values below are transcribed from `data/found/822_vaquero_em_free_transfer/results.json`
(sweep rows C=2.54, C=2.66), `vaquero_em_cycler_connections.OVERLAP_GRID_ICS`, and the `#828`
note — nothing re-derived.

```yaml
# data/manifold_connections.yaml  (schema: data/manifold_connection.schema.json v1.0, #838)
- id: em-vaquero-hetero-wu21c254-ws31c254-2026
  kind: heteroclinic
  model:
    type: cr3bp
    system: "Earth-Moon planar CR3BP"
    mass_ratio: 0.01215058439469525   # results.json:mu -- identical to both rows' recorded mass_ratio
  jacobi_constant: 2.54
  endpoints:
    # index 0 = unstable-manifold ORIGIN, index 1 = stable-manifold DESTINATION (directed)
    - uncatalogued:
        family: "Vaquero 2013 Sec 4.4.7 Earth-Moon 2:1 Periodic Cycler family member"
        x0: 1.0905363960533268          # OVERLAP_GRID_ICS[(2, 2.54)]
        ydot0: -0.8231863180949408
        period_nd: 5.941227735609639
        jacobi_constant: 2.54
        lambda_max: 3.149761131186353   # unstable (negative real saddle; see #822 note)
        derivation: >
          Re-converged at call time by build_vaquero_overlap_node (the same #799
          fixed-Jacobi corrector, tol=1e-12, half_crossings=3) from the #799-archived
          guess vendored in OVERLAP_GRID_ICS; |lambda| staleness-gated to 1e-6 relative.
          NOT a catalogued row -- #811 catalogued only her four printed-TOF family
          endpoints; this member's C sits at the 3:1 row's edge, not the 2:1 family's.
    - row_ref: vaquero-31-c254-em-cycler-2013
      identity_evidence: >
        #828 independent re-run: the connection module's re-converged 3:1 node reproduces
        this row's recorded state_nd to 3.1e-15, period_nd to 1.1e-13, jacobi_constant to
        5.3e-15, stability_index to 1.4e-08 relative (rtol=atol=1e-13 STM tolerance) --
        the transfer genuinely arrives at THIS catalogued orbit.
  connection:      # field names shared with crnbp seed_orbit_homoclinic.connection / results.json
    branch_u: -1
    branch_s: 1
    k_u: 40
    k_s: 27
    tau_u: 3.8298582936439947
    tau_s: 3.7380443596924104
    epsilon: 1.0e-4
    crossing_x: 0.23322507758762506
    crossing_xdot: -1.785596430399187
    newton_residual: 3.336901283531185e-10
  evidence:
    full_state_gap: 4.483722190016192e-06
    ydot_signs_match: true
    ghost_distance_from: 0.24275809881390115   # 243x the 1e-3 guard
    ghost_distance_to: 1.7999326727722387
    radau_gap: 1.8053763121228701e-06
    backward_reapproach: 1.0216784882076572e-09
    forward_reapproach: 3.8053502861397056e-04
    t_u_nd: 38.852089547100014
    t_s_nd: -27.627817271177918
    jacobi_drift: [2.7355895326763857e-12, 1.9741097645464833e-11]
    dv_kms: 0.0
  evidence_class: >
    SELF-CONSISTENCY ONLY (this project's full #822 verification battery, independently
    re-run bit-for-bit by #828): Vaquero prints no transfer state at any C -- never a
    digit-grade reproduction claim. Not novelty-bearing (her own pp.171-172 existence
    assertion; literature gate verdict `published`, #822 Sec 6).
  round_trip_note: >
    NOT demonstrated at this C: #822's reverse Wu(3:1)->Ws(2:1) demo is at C=2.60 only
    (#840, registered, would close this at the band edges).
  provenance:
    task_refs: ["#822", "#828"]
    data: "data/found/822_vaquero_em_free_transfer/results.json (sweep C=2.54)"
    module: "src/cyclerfinder/search/vaquero_em_cycler_connections.py"
    commit: "13ae76c0"
    notes:
      - "docs/notes/2026-08-11-822-vaquero-em-free-transfer.md"
      - "docs/notes/2026-08-12-828-vaquero-connection-tier-adjudication.md"

- id: em-vaquero-hetero-wu21c266-ws31c266-2026
  kind: heteroclinic
  model: {type: cr3bp, system: "Earth-Moon planar CR3BP", mass_ratio: 0.01215058439469525}
  jacobi_constant: 2.66
  endpoints:
    - row_ref: vaquero-21-c266-em-cycler-2013     # here the CATALOGUED endpoint is the ORIGIN --
      identity_evidence: >                         # the ordered-pair shape carries the asymmetry
        #828: re-converged 2:1 node reproduces this row's state_nd to 4.4e-15, period_nd
        to 5.2e-14, jacobi_constant to 8.0e-15, stability_index to 4.6e-13 relative.
    - uncatalogued:
        family: "Vaquero 2013 Sec 4.4.7 Earth-Moon 3:1 Periodic Cycler family member"
        x0: 0.8919375112041409          # OVERLAP_GRID_ICS[(3, 2.66)]
        ydot0: -0.7577709445440712
        period_nd: 6.283952207405823
        jacobi_constant: 2.66
        lambda_max: 12.783603246229166
        derivation: "build_vaquero_overlap_node, as above; not a catalogued row"
  connection:
    branch_u: 1
    branch_s: 1
    k_u: 37
    k_s: 32
    epsilon: 1.0e-4
    crossing_x: 0.869164          # #822 note Sec 4 table; full precision in results.json
    crossing_xdot: -0.575466
    newton_residual: 4.511e-10
  evidence:
    full_state_gap: 6.464e-06
    radau_gap: 3.515e-06
    ghost_distance_from: 0.5986
    ghost_distance_to: 0.5759
    t_u_nd: 41.5741
    t_s_nd: -34.6830
    dv_kms: 0.0
  evidence_class: "SELF-CONSISTENCY ONLY -- as the C=2.54 entry."
  round_trip_note: "NOT demonstrated at this C (reverse demo at C=2.60 only; #840)."
  provenance:
    task_refs: ["#822", "#828"]
    data: "data/found/822_vaquero_em_free_transfer/results.json (sweep C=2.66)"
    module: "src/cyclerfinder/search/vaquero_em_cycler_connections.py"
    commit: "13ae76c0"
    notes:
      - "docs/notes/2026-08-11-822-vaquero-em-free-transfer.md"
      - "docs/notes/2026-08-12-828-vaquero-connection-tier-adjudication.md"
```

The pending `#854` entry (drafted only when `#854` adjudicates in favor) would be a third item:
Kumar Table-5 C=2.54, `kind: heteroclinic`, `model.mass_ratio: 1.2150584270572e-2` (Kumar's own
printed mu — NOT the rows'), endpoint 0 = `row_ref: vaquero-31-c254-em-cycler-2013` with a
`model_note` stating the ~1.0e-8-relative mu difference and how node identity was established
across it, endpoint 1 = uncatalogued 2:1 member, `evidence_class` = "DIGIT-GRADE REPRODUCTION of
a published state (match_distance 7.611e-07 vs Kumar et al. 2026 Table 5; runner-up 1.307)" —
which is exactly the evidence-kind distinction the free-text comments currently cannot express
queryably.

## 6. What I am NOT certain about — reviewer sanity checks before approving

1. **The admission criterion (≥1 catalogued endpoint).** It deliberately excludes `#822`'s other
   11 grid points and all no-row systems (`#754`/`#759`/`#767`/`#781`). Alternative: admit every
   verified connection project-wide. I recommend against (the registry would become a second
   `data/found/`), but this is a scope judgment the user should confirm.
2. **The `#854` mu question.** Whether an endpoint at Kumar's mu (1.0e-8 relative from the row's
   recorded `mass_ratio`) may carry a `row_ref` at all, or must be `uncatalogued` with a soft
   `corresponds_to_row` hint. I lean `row_ref` + mandatory `model_note` (the family-member
   identity is the load-bearing claim and `#827` established it), but the gate cannot verify
   identity across differing mu mechanically — this is exactly the kind of call `#854`'s own
   adjudication should make, and the schema should not force it prematurely.
3. **Whether to shrink the two existing `#828` comments to short pointers** once their content is
   registry-resident (comment-only, parse-identical edit) or leave them verbatim as the
   historical record. Either is safe; duplication-drift argues for shrinking, audit-trail
   conservatism for leaving.
4. **Deferral of the row-side `connection_refs` field (the would-be v5.5).** If the user expects
   row-side machine queries soon (e.g. `#839`'s `orbit_class` adjudication tooling), adding it
   now alongside the registry may be cheaper than a second proposal later.
5. **Directionality convention.** Ordered endpoints (Wu origin first) + separate entries per
   direction, with `reverse_of` links, vs. one entry with a `directions:` list. I chose separate
   directed entries because `#822`/`#840` demonstrate directions at *different C values* and with
   *different evidence*, so one-entry-per-demonstration matches the actual evidence granularity —
   but the `#840` band-edge runs will stress-test this choice first.
6. **Naming.** `manifold_connections` vs `connections` (too generic next to
   `cycler_networks.yaml` / `inserts_into`) vs `heteroclinics` (excludes homoclinics). Ids
   follow the row-id spirit (`em-vaquero-hetero-…-2026`) but no convention exists yet for
   relation ids; the user may prefer another shape.

## 7. Summary

| Question | Recommendation | Confidence |
|---|---|---|
| Structural home needed? | Yes — prose is at 3 recurrences in 10 days, ~2.5KB unqueryable single-line comments, with `#854`/`#840`/`#839` queued against the same rows and an `orbit_class` adjudication (`#839`) that will want to query connections | High |
| Shape | Dedicated top-level registry `data/manifold_connections.yaml` + `data/manifold_connection.schema.json` v1.0 + `validate_connections.py` (the `#570` `cycler_network` / `empty_regions` precedent), NOT a row field, NOT a catalogue row | High on the two rejections; medium-high on registry vs row-field |
| Half-catalogued case | Ordered endpoint pair; each endpoint `row_ref`+mandatory `identity_evidence` OR inline `uncatalogued` descriptor; gate requires ≥1 `row_ref`; no symmetry assumed | High |
| Consistency | Intrinsic payloads untouched; shared field vocabulary; shared "not a promotion gate" clause; scope rule stated in all three descriptions; three-sites tension acknowledged, unification rejected as breaking | High |
| Migration/retro | No catalogue bump (no v5.5), no migration script; retroactively populate `#822`'s two row-touching entries at adoption (pure transcription); `#854` entry gated on `#854`'s verdict; `connection_refs` row field deferred | High on no-bump; medium on the deferral |

No schema, catalogue, or registry file has been modified by this task. All of the above is a
proposal for the user to review, decide, and (if approved) implement as an explicit
registry-creation + population task, exactly as `#707`'s and `#735`'s designs were handled by
`#708` and `#736`.
