# Schema v4 Upgrade & Backfill Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax. Work on `main` — do NOT branch. Golden-discipline throughout: every backfilled value must trace to a published source, never to a value our own code computed. Spec: `docs/spec.md` §16.7 (committed 9cdc6e6).

**Goal:** Make the catalogue *structurally honest* about each cycler's kind (single-ellipse / multi-arc / non-keplerian), remove the heliocentric+2-body biases that block moon-tour and n-body entries, and backfill the new fields from the sources we already hold — without breaking the canonical signature or any existing consumer.

**Architecture:** Three layers, in order. (1) **Loader** — parse the v4 fields onto `CatalogueEntry`, dispatch `fully_defined` by class, add a physical-invariant gate. (2) **Classification sweep** — tag every row's `cycler_class` conservatively (default single-ellipse; the 6 non-heliocentric → non-keplerian; multi-arc only on positive source evidence). (3) **Source-gated data backfill** — lift `invariants{}` / `cr3bp{}` / `period.basis` from the prose `notes` and `docs/refs/` PDFs into typed, testable fields. Validation dispatch and site rendering come last.

**Tech Stack:** Python 3.11 (uv), pytest, ruff, mypy; `data/catalogue.yaml`; loaders `src/cyclerfinder/data/catalog.py` + `tests/_catalogue_loader.py`; site `../cyclers.space` (Astro/TS).

---

## Why prose `notes` ARE a data source here (design note)

This plan treats the existing prose `notes` in `catalogue.yaml` as **first-class reference material**, mined for backfill. That needs saying because §16.7 called invariants "buried in notes" as if prose were the problem. It isn't.

The objection was never to the *provenance* of prose-sourced data — a number a human transcribed from Russell into a `note` is exactly as sourced as one in a typed field. The objection is to prose as the *validation surface*:

1. **A test cannot assert on free text.** `aphelion_ratio: 1.44` can be checked against a source by `test_russell_invariants`; "aphelion ratio about 1.44" cannot, without brittle regex that itself becomes a bug source.
2. **It isn't queryable or dispatchable.** The loader, the matcher, the site, and the validation gauntlet can't branch on a sentence. They can branch on a field.
3. **Ambiguity.** Prose hides units, which-arc, and rounding; a typed field with a `source` forces those to be made explicit at ingest, when the human still has the paper open.

So prose `notes` are a fine **capture** location and a poor **validation** location. The fix is not to distrust them — it's to **lift** them into typed fields *with their source attached*. That lift is Phases 3–5 of this plan, and the prose note is one of the sanctioned inputs to it.

---

## Data sources we hold (inputs to the backfill)

- **Prose `notes` in `data/catalogue.yaml`** — 1303 lines mention aphelion/turn/transit/jacobi/stability. Primary mining target for invariants already transcribed.
- **`docs/refs/` PDFs** (10): `russell-2004-dissertation.pdf` (multi-arc invariants + which sequences are *generic-return*, i.e. genuinely multi-arc), `rogers-2012-vinf-leveraging-cyclers` (S1L1-family a/e), `hollister-menning-1970` + `hollister-rall-1970` (Earth–Venus), `Low-Excess-Speed-Triple-Cyclers-of-Venus-Earth-and-Mars` (VEM n-body → `period.basis`), `genova-aldrin-2015-earth-moon-cycler` + `genova-2016-phobos-deimos-PADME` (planet-centric), `russell-strange-2009` Jovian/Saturnian multimoon (non-keplerian), `landau-longuski-2006`, `vasile-…-2005`.
- **Already-extracted tables** in catalogue rows (Russell 3.x/4.x, Rogers Table 1/4) — reuse, don't re-derive.
- **JPL three-body periodic-orbit catalog** (https://ssd-api.jpl.nasa.gov/doc/periodic_orbits.html) — *reference shape* for `cr3bp{}`; only a citable value source where a non-keplerian row corresponds to a catalogued family.

**Hard rule:** if a value isn't in a source we hold, it stays a `data_gaps[]` entry — never invented, never our-computed.

---

## Classification rules (the crux — do not shortcut)

`cycler_class` is **not** mechanically derivable from "null `a_au`". 199 rows have multi-seg + null `a_au`, but most are *single-ellipse* cyclers (McConaghy SnLm, Hollister) whose `a` is simply *unfilled* — a data gap, not a multi-arc orbit. Classify by this precedence:

1. **`non-keplerian`** — `primary != "Sun"` (rotating-frame/CR3BP). Mechanical. Exactly the 6 known rows: `arenstorf-em-figure8-1963`, `genova-aldrin-2015-em-3petal-cycler`, `wittal-2022-em-cycler-family`, `hernandez-2017-jovian-ieg-triple-family`, `russell-strange-2009-jovian-multimoon-family`, `russell-strange-2009-saturnian-multimoon-family`.
2. **`multi-arc`** — ONLY with positive source evidence that the orbit is *generic-return* (a different ellipse per leg, no single repeating `(a,e)`). This is a per-row allowlist built by reading the source, NOT a heuristic. Seed members: the Russell-Ocampo generic-return rows (`russell-ocampo-*` whose source sequence is E-E-M-M generic-return).
3. **`single-ellipse`** — the default for everything else, including rows with a currently-null `a_au` (those are `single-ellipse` with an `a`-shaped `data_gaps[]` entry, not multi-arc).

A row that looks multi-seg but is a single repeating ellipse (the flybys re-match V∞ to *the same* orbit) is `single-ellipse`. When in doubt, default `single-ellipse` and open a `data_gaps[]` asking for source confirmation — never guess `multi-arc`.

---

## Phase 1 — Loader & gate (no data change)

### Task 1.1: Parse v4 fields onto `CatalogueEntry`

**Files:**
- Modify: `src/cyclerfinder/data/catalog.py` (`CatalogueEntry` dataclass + its constructor/`from_dict` site)
- Test: `tests/data/test_catalog_v4_fields.py` (create)

- [ ] **Step 1: Write the failing test**

```python
# tests/data/test_catalog_v4_fields.py
from cyclerfinder.data.catalog import load_catalogue  # or the existing loader entrypoint

def test_cycler_class_defaults_to_single_ellipse():
    entries = {e.id: e for e in load_catalogue()}
    # a row with no cycler_class key reads as single-ellipse
    e = entries["s1l1-2syn-em-cpom"]
    assert e.cycler_class == "single-ellipse"

def test_v4_fields_present_and_defaulted():
    entries = {e.id: e for e in load_catalogue()}
    e = entries["s1l1-2syn-em-cpom"]
    assert e.orbit_elements_reference_frame == "heliocentric-inertial"
    assert e.invariants is None          # only set for multi-arc
    assert e.cr3bp is None               # only set for non-keplerian
    assert e.period_basis is None        # only set for n-body
```

- [ ] **Step 2: Run it, confirm it fails** (`AttributeError: cycler_class`). Run: `uv run pytest tests/data/test_catalog_v4_fields.py -q`

- [ ] **Step 3: Add fields + parsing.** On the dataclass add, with v3-compatible defaults:

```python
    cycler_class: str = "single-ellipse"        # single-ellipse | multi-arc | non-keplerian
    orbit_elements_reference_frame: str = "heliocentric-inertial"
    orbit_elements_center: str = "Sun"
    invariants: dict[str, Any] | None = None    # multi-arc only
    cr3bp: dict[str, Any] | None = None          # non-keplerian only
    period_basis: tuple[dict[str, Any], ...] | None = None  # n-body only
```

In the row→entry construction, read them from the dict (the `raw` round-trip already preserves unknown keys, so writeback is safe):

```python
    oe = row.get("orbit_elements") or {}
    cycler_class=row.get("cycler_class", "single-ellipse"),
    orbit_elements_reference_frame=oe.get("reference_frame", "heliocentric-inertial"),
    orbit_elements_center=oe.get("center", "Sun"),
    invariants=row.get("invariants"),
    cr3bp=oe.get("cr3bp"),
    period_basis=tuple((row.get("period") or {}).get("basis") or ()) or None,
```

- [ ] **Step 4: Run test → PASS.** Then full suite stays green: `uv run pytest -q -m "not slow"`
- [ ] **Step 5: Commit** — `data/catalog: parse schema-v4 fields (cycler_class, frame, invariants, cr3bp, period.basis)`

### Task 1.2: Dispatch `fully_defined` by class

**Files:** Modify `src/cyclerfinder/data/catalog.py` (`fully_defined` property); Test `tests/data/test_catalog_v4_fields.py`

- [ ] **Step 1: Failing test** — a `multi-arc` row with null `a_au` but full `invariants{}` + per-segment data + V∞ + no `data_gaps` must read `fully_defined == True`; a `non-keplerian` row is `fully_defined` iff its `cr3bp{}` identity tuple (`jacobi_constant`,`period_nd`,`stability_index`) is present.

```python
def test_fully_defined_multi_arc_uses_invariants_not_ae(make_entry):
    e = make_entry(cycler_class="multi-arc", a_au=None, e=None,
                   invariants={"aphelion_ratio": 1.44}, ...)  # full segments + vinf, no gaps
    assert e.fully_defined is True

def test_fully_defined_non_keplerian_uses_cr3bp(make_entry):
    e = make_entry(cycler_class="non-keplerian", a_au=None, e=None,
                   cr3bp={"jacobi_constant": 3.0, "period_nd": 6.2, "stability_index": 1.0}, ...)
    assert e.fully_defined is True
```

- [ ] **Step 2: Run → fail** (current code returns False on null `a_au`).
- [ ] **Step 3: Implement dispatch:**

```python
    @property
    def fully_defined(self) -> bool:
        if self.data_gaps:
            return False
        if self.cycler_class == "non-keplerian":
            cr = self.cr3bp or {}
            return all(cr.get(k) is not None for k in ("jacobi_constant", "period_nd", "stability_index"))
        if self.cycler_class == "multi-arc":
            core_ok = bool(self.invariants) and any(v is not None for v in (self.invariants or {}).values())
        else:  # single-ellipse
            core_ok = self.orbit_elements_a_au is not None and self.orbit_elements_e is not None
        if not core_ok:
            return False
        if not self.vinf_kms_at_encounters or any(v is None for _, v in self.vinf_kms_at_encounters):
            return False
        return bool(self.legs_tof_days) and all(t is not None for t in self.legs_tof_days)
```

- [ ] **Step 4: Run → PASS** + full suite green.
- [ ] **Step 5: Commit** — `data/catalog: dispatch fully_defined by cycler_class`

### Task 1.3: Physical-invariant gate (folds in #73 loader validation)

**Files:** Create `src/cyclerfinder/data/validate.py` (or extend existing validation in `catalog.py`); Test `tests/data/test_schema_invariants.py`

- [ ] **Step 1: Failing test** — a row with `cycler_class` in {multi-arc, non-keplerian} AND non-null top-level `orbit_elements.a_au` must raise/collect a `SchemaError`; a clean catalogue passes.

```python
def test_multi_arc_must_not_have_top_level_a():
    bad = {"id": "x", "cycler_class": "multi-arc", "orbit_elements": {"a_au": 1.6}}
    errs = validate_schema_invariants([bad])
    assert any("a_au" in m for m in errs)

def test_current_catalogue_passes_invariants():
    import yaml
    rows = yaml.safe_load(open("data/catalogue.yaml"))
    assert validate_schema_invariants(rows) == []
```

- [ ] **Step 2: Run → fail** (no validator yet). Note: this test passing on the real catalogue is a *ratchet* — Phase 2 must keep it green.
- [ ] **Step 3: Implement** `validate_schema_invariants(rows) -> list[str]`: for each row, if `cycler_class in {"multi-arc","non-keplerian"}` and `orbit_elements.a_au`/`e` is non-null → error; if `non-keplerian` and `primary == "Sun"` → error; if `period.basis` present, each entry has `pair`+`k`. Return collected messages (don't raise — report all).
- [ ] **Step 4: Run → PASS** + suite green.
- [ ] **Step 5: Commit** — `data/validate: schema-v4 physical-invariant gate (closes #73)`

---

## Phase 2 — Classification sweep (data: `cycler_class` only)

### Task 2.1: Tag every row's `cycler_class`

**Files:** Modify `data/catalogue.yaml` (add `cycler_class:` to all 233 rows); Create `scripts/classify_cycler_class.py` (one-shot, auditable); Test `tests/data/test_cycler_class_census.py`

- [ ] **Step 1: Failing census test** (the ratchet):

```python
def test_cycler_class_census():
    import yaml
    rows = yaml.safe_load(open("data/catalogue.yaml"))
    by_class = Counter(r.get("cycler_class", "MISSING") for r in rows)
    assert by_class["MISSING"] == 0                       # every row tagged
    assert by_class["non-keplerian"] == 6                 # the 6 non-heliocentric rows
    # multi-arc set is an explicit, source-justified allowlist (frozen here):
    multi = {r["id"] for r in rows if r.get("cycler_class") == "multi-arc"}
    assert multi == EXPECTED_MULTI_ARC                    # defined from source reading, see Task 2.2
```

- [ ] **Step 2: Run → fail.**
- [ ] **Step 3: Apply tags.** `scripts/classify_cycler_class.py` writes `cycler_class` per the precedence rules above: `primary != Sun` → non-keplerian; id in the source-built `MULTI_ARC_ALLOWLIST` → multi-arc; else single-ellipse. Run it, review the diff by eye (this is a data commit — confirm no single-ellipse cycler got mis-tagged multi-arc).
- [ ] **Step 4: Run census + `validate_schema_invariants` (Task 1.3) → PASS.** The invariant gate will FAIL here if any newly-tagged multi-arc/non-keplerian row still carries a top-level `a_au` — if so, that `a_au` was either wrong (remove, add `data_gaps`) or the row is actually single-ellipse (retag). Resolve each, don't suppress.
- [ ] **Step 5: Commit** — `data/catalogue: tag cycler_class on all 233 rows (conservative, source-gated)`

### Task 2.2: Build & document the multi-arc allowlist from source

**Files:** Create `docs/notes/multi-arc-classification.md`; this *precedes* Task 2.1's `EXPECTED_MULTI_ARC`.

- [ ] **Step 1:** For each candidate Russell-Ocampo `*` row, read `russell-2004-dissertation.pdf` (use Read on the PDF; if a page won't extract, log it as a gap — do NOT guess). Record per row: is the source sequence *generic-return* (distinct per-leg orbit → multi-arc) or a single repeating ellipse (→ single-ellipse)? Cite the table/section.
- [ ] **Step 2:** Write the verdict table to `docs/notes/multi-arc-classification.md` with citation per row. The set of multi-arc verdicts becomes `EXPECTED_MULTI_ARC` / `MULTI_ARC_ALLOWLIST`.
- [ ] **Step 3: Commit** — `docs: multi-arc classification verdicts with Russell-2004 citations`

> Ordering note: do Task 2.2 *before* 2.1's data write (2.1 consumes the allowlist). They are split so the source-reading judgement is reviewable on its own.

---

## Phase 3 — `invariants{}` backfill for multi-arc (source-gated)

### Task 3.1: Lift published invariants from notes + Russell into typed fields

**Files:** Modify `data/catalogue.yaml` (multi-arc rows); Test `tests/data/test_russell_invariants.py`

- [ ] **Step 1: Failing golden test** — for each multi-arc row with a published aphelion ratio, assert the typed `invariants.aphelion_ratio` equals the sourced value (EXPECTED side from the paper/notes, never our compute):

```python
RUSSELL_INVARIANTS = {  # value : source citation, transcribed from the paper
    "russell-ocampo-2.1.1+2-case2": {"aphelion_ratio": 1.444},  # Russell 2004, Table X
}
def test_multi_arc_invariants_match_source():
    rows = {r["id"]: r for r in yaml.safe_load(open("data/catalogue.yaml"))}
    for rid, exp in RUSSELL_INVARIANTS.items():
        inv = rows[rid].get("invariants") or {}
        for k, v in exp.items():
            assert inv.get(k) == pytest.approx(v), rid
```

- [ ] **Step 2: Run → fail.**
- [ ] **Step 3: Backfill.** For each multi-arc row: move the aphelion-ratio/transit-time/turn-ratio value from its prose `note` into `invariants{}`, attach `source`. Where the note's value is ambiguous (units/which-arc), open the PDF and resolve; if unresolvable, leave the field null and add a `data_gaps[]` entry. Also relocate the now-redundant prose (or keep as human context, but the typed field is authoritative). The existing top-level `aphelion_au` that was shoehorned in stays only if it's a true cycle-level aphelion; otherwise move it under `invariants`.
- [ ] **Step 4: Run golden test + invariant gate + full suite → PASS.**
- [ ] **Step 5: Commit** — `data/catalogue: lift Russell multi-arc invariants from notes to typed fields`

---

## Phase 4 — `cr3bp{}` backfill for non-keplerian (source-gated)

### Task 4.1: Populate the CR3BP identity tuple where a source gives it

**Files:** Modify `data/catalogue.yaml` (the 6 non-keplerian rows); Test `tests/data/test_cr3bp_identity.py`

- [ ] **Step 1: Failing test** — for each non-keplerian row where a source publishes it, `orbit_elements.cr3bp` carries the sourced `jacobi_constant`/`period_nd`/`stability_index` (+ `mass_ratio`, `family`); EXPECTED values cite the source (e.g. `genova-aldrin-2015`, `russell-strange-2009`, or the JPL three-body catalog family).
- [ ] **Step 2: Run → fail.**
- [ ] **Step 3: Backfill from the PDFs we hold.** Realistic outcome: some rows (Arenstorf, Genova EM cycler) have direct source values; the Jovian/Saturnian multimoon families may only give qualitative data → populate what's sourced, mark the rest `data_gaps[]` with `source_hint`. Do NOT pull a value from the JPL API for a row unless that row demonstrably *is* that catalogued family. The honest "all-null + note" is replaced by "sourced fields + explicit gaps", not by fabricated identity.
- [ ] **Step 4: Run → PASS** + suite green.
- [ ] **Step 5: Commit** — `data/catalogue: backfill cr3bp identity for non-keplerian rows (sourced + gaps)`

---

## Phase 5 — `period.basis` for n-body (VEM) rows

### Task 5.1: Add beat-period basis to multi-pair cyclers

**Files:** Modify `data/catalogue.yaml` (VEM / triple-cycler rows); Test `tests/data/test_period_basis.py`

- [ ] **Step 1: Failing test** — VEM triple-cycler rows (from `Low-Excess-Speed-Triple-Cyclers…pdf`) carry `period.basis` as a list of `{pair,k}` whose source-stated resonances are present; 2-body rows are unchanged (flat `pair/k`, no `basis`).
- [ ] **Step 2: Run → fail.**
- [ ] **Step 3: Backfill** the basis for n-body rows from the source; leave 2-body rows alone (the flat form remains valid). Confirm `validate_schema_invariants` accepts both forms.
- [ ] **Step 4: Run → PASS** + suite green.
- [ ] **Step 5: Commit** — `data/catalogue: add period.basis to n-body (VEM) cyclers`

---

## Phase 6 — Consumers: validation dispatch + site

### Task 6.1: Validation dispatch by class

**Files:** Modify the rediscovery/validation harness (`tests/test_catalogue_rediscovery.py` and/or the data-validation-hardening entrypoint); Test alongside.

- [ ] **Step 1: Failing test** — the validator selects anchors by `cycler_class` per §16.7.5: single-ellipse → V∞+`(a,e)`+period; multi-arc → V∞+period+`invariants`; non-keplerian → `cr3bp` tuple (no Kepler check). Assert a multi-arc row is checked against its `invariants`, and that the resonance constructor is NOT invoked on a multi-arc row.
- [ ] **Step 2–4:** Implement the dispatch (a thin `anchors_for(entry)` switch), run → PASS, suite green.
- [ ] **Step 5: Commit** — `validate: dispatch expected-output anchors by cycler_class`

### Task 6.2: Site renders per-class columns

**Files:** `../cyclers.space/src/lib/catalogue.ts` (+ `types.ts`, the catalogue table component)

- [ ] **Step 1:** Extend the TS `CyclerEntry` type with the v4 fields; `isFullyDefined` mirrors the class dispatch from Task 1.2 (multi-arc → invariants, non-keplerian → cr3bp). 
- [ ] **Step 2:** Table shows the right identity column per class (single-ellipse: a/e; multi-arc: aphelion ratio + V∞; non-keplerian: Jacobi/period/stability) instead of a row of em-dashes. `npm run build` succeeds (sync step pulls canonical catalogue).
- [ ] **Step 3: Commit (cyclers.space repo)** — `catalogue: render schema-v4 per-class identity columns`

---

## Self-Review

- **Spec coverage:** §16.7.1 (cycler_class) → P1.1/P2; §16.7.2 (frame/units) → P1.1; §16.7.3 (period.basis) → P1.1/P5; §16.7.4 (invariants/cr3bp) → P3/P4; §16.7.5 (dispatch) → P6.1; physical-invariant gate → P1.3. All covered.
- **Golden-discipline:** every EXPECTED value in Phases 3–5 cites a source we hold; unresolved values become `data_gaps[]`, never invented. The prose-`note` lift is a *source* lift, not a re-derivation.
- **Non-breaking:** all loader fields default to v3 behaviour; signature untouched (v4 fields are non-signature, same as v2); the invariant gate is a ratchet that starts green after P2 and stays green.
- **Conservatism enforced:** multi-arc requires a documented source verdict (P2.2) — the null-`a` heuristic is explicitly rejected, so single-ellipse cyclers with unfilled `a` are not mis-tagged.
- **Risk:** Phase 4 yield is uncertain — the Jovian/Saturnian multimoon families may have little sourced CR3BP identity; the plan absorbs this as explicit gaps rather than fabrication. Phase 2.2 PDF extraction may hit un-extractable pages (log as gaps, don't guess).
- **Ordering:** P1 (code) is independently valuable and ships first; P2.2 precedes P2.1; data phases (3–5) are independent of each other and can be parallelised across subagents; P6 depends on P1–P5.
