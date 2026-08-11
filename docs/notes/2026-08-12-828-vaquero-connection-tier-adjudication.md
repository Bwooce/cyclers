# `#828`: do `#822`'s Vaquero free transfers upgrade any `#811` row's validation level?

**Task:** `#828`, registered 2026-08-11 (found during `#822`), dispatched 2026-08-12. Adjudicate
whether `#822`'s 13/13-grid-point verified Wu(2:1)→Ws(3:1) heteroclinic free transfers (+ the
reverse direction at C=2.60) upgrade any of `#811`'s six `vaquero-*` catalogue rows' validation
levels, per this project's own `#388` promotion discipline — and if warranted, perform the
writeback under full catalogue-edit ratchet discipline. The sibling of `#826` (`#820`'s
Russell-row analogue of the same adjudication pattern).

---

## Verdict (read this first)

**ZERO TIER PROMOTIONS on all six rows. All six stay V1.** Two of the six rows earn a
comment-only `ADDED EVIDENCE (… NO LEVEL CHANGE)` annotation; four earn nothing at all, because
`#822` computed nothing that involves them.

**`#822`'s own registration framing — "the natural upgrade path toward V2/V3-class connection
evidence for the `#811` rows" — is WRONG, and is retracted here.** It was written before the
connection was computed (the dispatch correctly flagged it as a hypothesis, not a conclusion).
A heteroclinic connection is not a tier path at all, for either family. The live V2 path for
this family pair is `#823`, and `#823`'s two targets are **disjoint** from the two rows `#822`
actually touches.

`#822`'s connection evidence itself is genuine, admissible and fully reproducible — that was
independently confirmed here, from scratch (Sec. 3). The no-promotion verdict is *not* a doubt
about `#822`'s numbers. It is that the numbers are evidence of a **different kind** from what
the §14 ladder measures.

---

## 1. Which rows are even in scope

`#811`'s six rows, with the Jacobi constant each sits at, against `#822`'s overlap band
`C ∈ [2.54, 2.66]`:

| row | family | C | in `#822`'s band? |
|---|---|---|---|
| `vaquero-21-c198-em-resonant-po-2013` | 2:1 | 1.98 | **no** |
| `vaquero-21-c246-em-cycler-2013` | 2:1 | 2.46 | **no** |
| `vaquero-21-c247-em-cycler-2013` | 2:1 | 2.47 | **no** |
| `vaquero-21-c266-em-cycler-2013` | 2:1 | 2.66 | **yes** — Wu origin at the band's high edge |
| `vaquero-31-c254-em-cycler-2013` | 3:1 | 2.54 | **yes** — Ws destination at the band's low edge |
| `vaquero-31-c313-em-resonant-po-2013` | 3:1 | 3.13 | **no** |

**Four of the six rows are outside the overlap band entirely.** `#822` computed no connection
involving those orbits, so for them the promotion question does not even arise — there is zero
evidence of any kind, not weak evidence. This is the single largest fact of the adjudication and
`#822`'s "the `#811` rows" framing (plural, undifferentiated) obscured it.

Of `#822`'s 13 same-C connections, each pairs a 2:1 node with a 3:1 node **at the same C**. Only
two of those 26 node instances are catalogued rows: the 2:1 node at C=2.66 and the 3:1 node at
C=2.54. The counterpart node in each of those two connections (the 3:1 node at C=2.66, the 2:1
node at C=2.54) is *not* a catalogued row. So neither of `#822`'s two row-touching connections
has both of its endpoints in the catalogue.

**Node identity is exact, not assumed.** `vcc.OVERLAP_GRID_ICS[(2, 2.66)]` and `[(3, 2.54)]` are
byte-identical to those rows' recorded `state_nd` / `period_nd` / `stability_index`; and the
*re-converged* node (`build_vaquero_overlap_node`, not the vendored guess) reproduces them too —
see Sec. 3.

## 2. The gates, clause by clause

The authoritative gate definitions read for this adjudication: spec §14's table + its V2/V3
class-split notes; `src/cyclerfinder/data/validate.py::validate_validation_level` +
`_LEVEL_EVIDENCE`; `src/cyclerfinder/data/provenance.py::classify_validation`; and the per-class
gauntlet modules in `src/cyclerfinder/data/validation/` (there is no planar-CR3BP-specific
module; `v2_3d.py` carries the CR3BP-class V2-ballistic criteria and is the closest written
gate for these rows).

**The discriminating constraint, stated once.** Every rung of §14 is a statement about **one
object's own trajectory under increasing model fidelity** — V0 internal consistency, V1
independent solver re-derivation, V2 multi-lap persistence, V3 real ephemeris, V4 external code.
A heteroclinic connection is a statement about **two objects, in the same model, at the same
fidelity**. It moves sideways on the ladder, not up. This is not a novel reading invented for
this task: the catalogue schema already says it in its own words for the only two connection
payloads it carries — `ccr4bp_provenance` (schema v5.3, `#708`) and `crnbp_provenance` (v5.4,
`#736`) are both annotated *"Provenance/audit only — not a promotion gate (validation_level
carries that)"*.

Checked rung by rung against what `#822` actually computed, for the two in-band rows:

- **V3-ballistic / V3-powered — structurally impossible.** Both require ephemeris realisation
  (phase-matching to a real launch window, an ephemeris-mode horizon TCM over 3–5 laps, and for
  V3-powered, encounters independently confirmed on the real ephemeris). `#822` is pure planar
  CR3BP throughout; the module touches no ephemeris, no epoch, and no ΔV budget. Nothing in the
  V3 criteria is even addressed, let alone met.
- **V2-ballistic — structurally inapplicable AND mechanically ineligible.** `v2_3d.py`'s written
  criteria (spec §14, `V2_N_CYCLES_MIN = 3`, module-constant and explicitly "NOT test-tunable"):
  propagate **the row's own IC** for ≥3 consecutive periods in the row's defining model and
  assert the state stays within a bounded distance (50,000 km floor) of the original IC at each
  return. `#822` propagates *manifold legs* — trajectories seeded `x0 + ε·v` off the orbit and
  integrated ~4.5–9 periods precisely so they **asymptotically depart** it. That is the exact
  opposite of the bounded-drift persistence V2 measures; a longer, cleaner manifold leg is
  *worse* V2 evidence, not better. Independently, both in-band rows are **unstable**
  (|λ| = 4.428 and 11.26), and each row's own existing `_LEVEL_EVIDENCE` entry already records
  "NOT V2 (an UNSTABLE orbit cannot satisfy V2-ballistic bounded-drift-over->=3-laps)". Two
  independent grounds, either sufficient.
- **V2-powered — inapplicable.** These are ballistic CR3BP periodic orbits with no documented
  per-cycle maintenance maneuver to execute; there is no powered budget to hold them to.
- **V1 — already held, and `#822` adds nothing to it.** `#822`'s node re-derivation calls
  `cp.correct_symmetric_fixed_jacobi` — **the same corrector `#799` used**, invoked explicitly
  as a staleness guard (it checks the recomputed |λ| against `#799`'s archived value to 1e-6
  relative). Same-corrector re-convergence is not an independent solver cross-check. And
  `#822`'s independent-Radau cross-check is applied to the **manifold legs' section crossings**,
  not to either node's own period-map closure — so it is not a second-integrator confirmation of
  the orbit's own periodicity either. `#822` supplies no new evidence at the rung the rows
  already hold.
- **V0** is the internal-consistency floor, already implied.

**The forcing device (what actually settles it).** `validate_validation_level` is a pure
membership check against `_LEVEL_EVIDENCE`: a promotion *is* an evidence string somebody has to
write. Drafting one for `vaquero-21-c266` at V2 makes the failure immediate — there is no lap
count to state, no drift figure to state, and the orbit is unstable. If the evidence entry
cannot be written honestly, there is no promotion. (This is the same instrument `#826` used, and
the same "when in doubt, V0" discipline `validate_validation_level`'s own docstring states.)

## 3. Independent re-verification of `#822`'s numbers (golden discipline)

Per the dispatch and `[[feedback_orbit_closure_discipline]]`, `#822`'s "verified" framing was not
taken at face value. `find_free_transfer` was re-run from scratch at three C values — the
mandated primary C=2.60 plus the two that actually matter here, the band edges C=2.66 and
C=2.54, at `#822`'s own recorded per-C settings.

**All three reproduce `#822` bit-for-bit** (`|Δ| = 0.000e+00` on both crossing components at all
three C), including the seed/refine/converge counts and the whole verification battery:

| C | k | branches | residual | 4-state gap | Radau gap | ghost | t_u, t_s (nd) | passed |
|---|---|---|---|---|---|---|---|---|
| 2.60 | (45,30) | (+1,+1) | 8.037e-11 | 5.704e-06 | 4.097e-06 | 0.5202/0.4809 | +43.4534, −30.8432 | yes |
| 2.66 | (37,32) | (+1,+1) | 4.511e-10 | 6.464e-06 | 3.515e-06 | 0.5986/0.5759 | +41.5741, −34.6830 | yes |
| 2.54 | (40,27) | (−1,+1) | 3.337e-10 | 4.484e-06 | 1.805e-06 | 0.2428/1.7999 | +38.8521, −27.6278 | yes |

`n_seeds/n_refined/n_converged` also matched exactly: 331/2/1, 31/1/1, 225/6/1. `#822`'s result
is real and reproducible.

**Node identity, independently checked** — the re-converged node vs. the catalogue row:

| row | Δ`state_nd` (x0, ẏ0) | Δ`period_nd` | Δ`jacobi_constant` | Δ`stability_index` |
|---|---|---|---|---|
| `vaquero-21-c266-em-cycler-2013` | 0.0, 4.4e-15 | 5.2e-14 | 8.0e-15 | 4.6e-13 rel |
| `vaquero-31-c254-em-cycler-2013` | 0.0, 3.1e-15 | 1.1e-13 | 5.3e-15 | 1.4e-08 rel |

The `stability_index` figures are at `#799`'s own `rtol=atol=1e-13` STM tolerance. A margin worth
recording rather than glossing (`[[feedback_verify_automated_ghost_guard_booleans]]`): the first
pass used `barden_stability`'s **default** `rtol=atol=1e-12` and c254's Barden ν came out
1.26e-04 relative off the recorded value — chased down rather than waved away, and it is purely
STM-integration tolerance on a |λ|=11.26 orbit (1e-12 → 1.3e-04, 1e-13 → 1.4e-08; 1e-14 is
*worse*, 1.4e-06, because SciPy clamps `rtol` at 2.2e-14 and an over-tight `atol` adds noise).
Not a discrepancy in the recorded value.

**This is no-downgrade evidence, not promotion evidence** — the
`[[feedback_bugfix_invalidates_past_searches]]` cuts-both-ways check `#826` also ran. `#822`'s
1e-6-relative |λ| staleness gate passing at (2, 2.66) and (3, 2.54), plus the table above,
re-confirms those two rows' recorded orbit identity from a code path built eleven tasks later.
Neither row degrades.

## 4. Is promotion even the right *category*? (the dispatch's question 3)

No — and the catalogue has no field for the right category either.

**How prior connection work was recorded, which is the best available precedent:**

- `#767` (Saturn-Titan homoclinic), `#781` (Neptune-Triton homoclinic), `#759` (Jupiter-Europa
  heteroclinic): **zero catalogue references each.** Weak precedent, though — those systems have
  no catalogue rows at all, so there was nothing to annotate.
- `#786` (Earth-Moon Class-1, targeting `casoliva-7-3b`/`7-3c`, which *are* catalogued rows):
  no writeback — but it was a clean negative, so again weak precedent.
- `#766`/`#779` (Jupiter-Europa 3:4-LO homoclinic) is the one **positive** precedent with a
  catalogue row in play, and it is instructive: the connection was recorded as a nested
  `crnbp_provenance.torus.seed_orbit_homoclinic.connection{}` **provenance block**, with the
  row's `validation_level` untouched — and the schema text for it says so explicitly.

Both existing connection payloads (`#708`'s `ccr4bp_provenance.connection`, `#736`'s
`crnbp_provenance…connection`) are **intrinsic** to their row's own object: the row *is* the
connection object, or the connection belongs to the row's own seed orbit. `#822`'s connection is
different in kind — it runs **between two different rows' orbits**. There is no schema field for
cross-row connection evidence, and inventing one here would be wrong on process as well as on
merit: both existing connection blocks arrived via a **separate, user-approved design-proposal
task** (`#707` → v5.3, `#735` → v5.4) before their schema bump. **No schema change was made
here.** Registered instead as `#838`.

What *was* done is the narrow, established, zero-risk thing: the `ADDED EVIDENCE (…, no level
change)` comment convention already in this catalogue (`#181`, on `russell-ch4-9.353Gg2` and
`russell-ch4-3.78Gg3`), applied to the two in-band rows' `validation_level` comments only. No
field value, no enum, no census count, no `_LEVEL_EVIDENCE` entry changes — verified by parsing
the catalogue before and after and asserting the parsed data is **identical** (`git diff --stat`:
2 insertions, 2 deletions, both comment text).

`#826` declined even this on its own rows; the difference is substantive, not a softer standard
here. `#826`'s evidence was inadmissible (bend-infeasible closures), off-family
(`5.30ggF3`), or non-independent (`mcconaghy`), so citing it would have misled. `#822`'s
evidence is admissible, independently reproduced here, and literally about these two orbits —
it is worth a reader knowing, provided it is labelled for what it is.

The annotations state the asymmetry honestly: `#822` demonstrated the reverse direction
Wu(3:1)→Ws(2:1) **only at C=2.60**, so neither in-band row has both directions demonstrated at
its own C, and **no round-trip transport claim is available for either row**.

## 5. What the real V2 path is

`#823` — the already-registered multi-lap bounded-drift V2-ballistic candidacy run — is the only
currently-visible tier path for this family pair, and it is untouched by `#822`. Its targets are
the two **linearly stable** rows, `vaquero-21-c198` (ν = 0.4852) and `vaquero-21-c246`
(ν = −0.9274), whose `_LEVEL_EVIDENCE` entries record "NOT V2 (no multi-lap bounded-drift run
has been performed this task)" — a *gap*, not an impossibility, unlike the four unstable rows'
"an UNSTABLE orbit cannot satisfy V2-ballistic", which is mechanical ineligibility.

Note the disjointness: `#823`'s two targets (C=1.98, 2.46) and `#822`'s two in-band rows
(C=2.66, 2.54) share no member. **The free transfer only ever touches unstable members** — it
has to, since Vaquero's own p.172 statement is an *unstable-to-unstable* one and `#799` found
the 2:1 family's stability transition at C≈2.46/2.47, below the overlap band. The two facts are
the same fact: the property that makes an orbit connectable (a saddle) is the property that
disqualifies it from V2-ballistic. Worth stating plainly, because it means no amount of further
connection work can ever promote an in-band member.

## 6. Follow-ups registered

- `#838` — schema design proposal (user approval required before any bump, per the `#707`/`#735`
  precedent): should the catalogue carry a **cross-row connection-evidence** object or field, for
  a verified manifold connection running *between two different rows' orbits*? Neither existing
  connection payload fits (both are intrinsic to a single row's own object). **The motivating
  case must be stated precisely**: `#822` is the *first* case in this corpus of a verified
  connection with **at least one** endpoint that is a catalogued row distinct from the other
  endpoint — *not* "the first connection between two catalogued rows", since neither of `#822`'s
  two row-touching connections has both endpoints catalogued (Sec. 1), and `#754`/`#759`'s
  Jupiter-Europa 3:4↔5:6 pair is no precedent at all, that system having no catalogue rows.
  The design therefore has to decide whether a *half-catalogued* connection is expressible, not
  assume the symmetric two-row case. Includes the prior question of whether such an object
  should be its own row rather than a field.
- `#839` — `vaquero-31-c313-em-resonant-po-2013` sits at **C=3.13, inside** Kumar-Rawat-
  Rosengren-Ross 2026's published `C_J ∈ [3.00, 3.15]` heteroclinic band, which is `#827`'s
  digit-grade reproduction target. It is also one of the two `resonant_po` rows — classed that
  way for "**no demonstrated transport utility**". If `#827` lands a heteroclinic touching that
  member, an `orbit_class` question (`resonant_po` → `cycler`) goes live. That is a
  *different kind* of writeback from a tier bump, and it is the only place this family's
  connection work could actually move a catalogue field. Gated on `#827`.
- `#840` — `#822`'s reverse direction Wu(3:1)→Ws(2:1) was demonstrated at C=2.60 only. Running
  it at the two band edges (C=2.54, 2.66) would give each in-band catalogued row a demonstrated
  round trip at its own C. Cheap (the machinery exists and the time-reversal symmetry guarantees
  existence), and it would strengthen the two annotations from one-way to round-trip — though,
  per Sec. 2, still not a tier promotion.

## 7. Verification run

- Independent re-run script and log: scratchpad (`reverify_828.py`), results summarised in
  Sec. 3. No `#822` artefact was modified; `data/found/822_vaquero_em_free_transfer/results.json`
  was read only.
- Catalogue edit is **comment-only**: `yaml.safe_load` of the pre-edit and post-edit files
  compares **equal**.
- Full `uv run pytest tests/data tests/search -q` ratchet (never a subset, per
  `[[feedback_catalogue_edits_run_all_ratchets]]`): exit 0, zero FAILED/ERROR.
- `uv run ruff check .`, `uv run ruff format --check .`, full `uv run mypy src tests`: clean.
