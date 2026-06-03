# Data‑Validation Hardening Implementation Plan

> **For agentic workers:** REQUIRED SUB‑SKILL: superpowers:subagent‑driven‑development or executing‑plans. Checkbox steps. Work on `main` — do NOT branch (project rule). **Run AFTER the `reproduction-unlocks` workflow lands** (Tasks 3+ depend on the Russell `(a,e)` fills).

**Goal:** Turn "we reproduced N rows" into a *defensible, machine‑checked* claim by (1) structuring provenance + model fidelity, (2) classifying every reproduction into cross‑validated / consistency‑checked / unvalidated tiers, and (3) adding source‑independent physical‑consistency invariants, multi‑source corroboration scoring, and falsification guards.

**Why:** Today provenance is free‑text `note:` fields — a human can see "Rogers (a,e) vs Russell V∞" independence, but no test can *enforce* it, so a circular check (same‑source a,e and V∞) is indistinguishable from a real cross‑check. And the S1L1 5.65 vs 4.99 episode showed **model‑fidelity mismatch** is a distinct, unguarded failure mode. This plan makes both machine‑checkable.

**Tech Stack:** Python 3.11, pyyaml, pytest + uv + ruff + mypy. Touches `data/catalogue.yaml`, `tests/_catalogue_loader.py`, `src/cyclerfinder/search/resonant_construct.py`, a new `data/provenance.py` (source registry).

**Scope guard (YAGNI):** structured tags go ONLY on the derivation‑relevant fields — `orbit_elements.{a_au,e,perihelion_au,aphelion_au}`, `vinf_kms_at_encounters[].vinf_kms`, `trajectory…tof_days`, `period`. Not every field on every row.

---

## Task 1: Source + fidelity registry (no catalogue edits yet)

**Files:** Create `src/cyclerfinder/data/provenance.py`; Test `tests/data/test_provenance_registry.py`.

- [ ] Define a `SOURCE_REGISTRY: dict[str, str]` mapping short keys → full citations: `rogers-2012-t1`, `russell-2004-t34`, `russell-2004-t39_311`, `russell-2004-t49_413`, `mcconaghy-2002`, `mcconaghy-2006`, `spec-9`, `hollister-1970-t3`, `friedlander-1986`, plus pseudo‑sources `derived` and `computed`.
- [ ] Define `Fidelity = Literal["circular-coplanar", "analytic-ephemeris", "real-de440"]`.
- [ ] Test: registry keys are unique, citations non‑empty; `derived`/`computed` present.
- [ ] Commit `data/provenance: source + fidelity registry`.

## Task 2: Physical‑consistency invariants — source‑INDEPENDENT (run on ALL rows now)

**Files:** Test `tests/data/test_physical_invariants.py`. *(No dependency on the workflow — can run first.)*

- [ ] Parametrise over every catalogue row; for each non‑null field assert the physics that must hold regardless of source — these catch transcription/unit errors with **no second source needed**:

```python
# per row, where the fields exist:
#  - orbit_elements: a ≈ (peri+apo)/2 (rel 1e-2); e ≈ (apo-peri)/(apo+peri); 0 ≤ e < 1; peri ≤ a ≤ apo
#  - V∞ ≥ 0 and < 15 km/s (E-M/E-V regime sanity)
#  - period.years > 0; period.years ≈ period.k * synodic_period(pair)/yr within 5%
#  - for a ballistic 2-body cycler with (a,e): peri ≤ r_inner_body AND apo ≥ r_outer_body
#    (else the orbit can't encounter both — a construct_resonant_cycler precondition)
#  - leg tof_days (when all present) sum ≈ period within tolerance
```

- [ ] These are golden in the strict sense (pure physics, no fitted values). Expect some rows to FAIL — that surfaces real data errors; fix the data or xfail with a documented reason. Commit.

## Task 3: Back‑fill `source` + `fidelity` tags on derivation‑relevant fields

**Files:** `data/catalogue.yaml`; `tests/_catalogue_loader.py` (parse the tags). *(AFTER the workflow — it adds Russell `(a,e)`.)*

- [ ] For each derivation‑relevant field that already cites a source in its `note:`, add `source: <key>` and `fidelity: <key>` keys (back‑fill from the existing prose — cheap, mechanical, per‑row). Leave untagged where genuinely unknown.
- [ ] Extend the loader to read these into `CatalogueEntry` (e.g. `vinf_sources: dict[body,str]`, `orbit_source: str|None`, `*_fidelity`).
- [ ] Test: loader round‑trips the tags; untagged fields default to `None`. Commit.

## Task 4: Tiered validation classifier + provenance‑aware reproduction gate

**Files:** `src/cyclerfinder/data/provenance.py` (classifier); Test `tests/data/test_resonant_reproduction.py`.

- [ ] `classify_validation(orbit_source, vinf_source, *, same_fidelity) -> Tier` where `Tier ∈ {CROSS_VALIDATED, CONSISTENCY_CHECKED, UNVALIDATED}`:
  - CROSS_VALIDATED iff `orbit_source` and `vinf_source` are both set, **different**, and same fidelity.
  - CONSISTENCY_CHECKED iff both set but same source (catches transcription only).
  - UNVALIDATED otherwise.
- [ ] Reproduction gate: for each row with `(a,e)`+V∞, run `construct_resonant_cycler`, assert computed V∞ matches sourced V∞ within the fidelity‑appropriate tolerance, AND assert the row's claimed tier is actually achievable (fails loudly if a row claims CROSS_VALIDATED but shares a source/fidelity). Over‑determine: also check ToF and period where independently sourced.
- [ ] **Validation‑tier census ratchet** (mirror `EXPECTED_COVERAGE`): freeze `{CROSS_VALIDATED: N, CONSISTENCY_CHECKED: M, UNVALIDATED: K}` so the strong set is visible and monotone. Commit.

## Task 5: Multi‑source corroboration scoring

**Files:** `data/catalogue.yaml` (optional `corroborations:` list on key quantities); Test `tests/data/test_corroboration.py`.

- [ ] Where ≥2 independent sources give a quantity (e.g. S1L1 V∞: Russell 4.99, McConaghy 4.7, spec 5.65), record them and a derived **agreement spread**. Classify: `strongly-sourced` (≥2 independent within tol), `single-sourced`, or `disputed` (spread > tol → fidelity flag).
- [ ] Test: the disputed set is exactly the known fidelity‑mismatch quantities (e.g. S1L1 5.65/3.05 vs 4.99/5.10) — documents, not hides, the spread. Commit.

## Task 6: Falsification + independent‑tool guards (prove the checks have teeth)

**Files:** `tests/data/test_validation_falsification.py`.

- [ ] **Negative gate:** feed deliberately‑mismatched `(a,e)`/V∞ to the reproduction check and assert it FAILS — proves the gate isn't a no‑op.
- [ ] **Independent‑code cross‑check:** verify `construct_resonant_cycler`'s V∞ for an orbit agrees with an independent path (the existing `lambert`/lamberthub crosscheck or the optimiser) on the same geometry — reduces *code‑bug* risk, complementing the *data* cross‑check. Commit.

## Task 7: Wire into loader schema validation (closes pending #73)

**Files:** `tests/_catalogue_loader.py` / a `validate_catalogue` entry point; Test.

- [ ] Any `source`/`fidelity` tag must be a registry key; any tagged value must satisfy the Task‑2 physical ranges; a row claiming a validation tier must back it with the right provenance. One `validate_catalogue()` that fails CI on a malformed entry (the long‑pending schema‑validation task).

## Task 8: Full‑suite + lint + type gate; refresh `data/OUTSTANDING.md`

- [ ] `uv run pytest -q`; ruff + mypy. Update OUTSTANDING with the new tier census (cross‑validated vs consistency‑checked counts). Commit.

---

## Further data‑validation improvements (proposed; pick per appetite)

- **Fidelity‑matched comparison everywhere** (Task 1/3/4): the single highest‑value addition — the S1L1 5.65‑vs‑4.99 bug was purely a fidelity mismatch. Never compare across fidelities; tag and enforce.
- **Derivation audit trail:** for `derived` values, store `derived_from: [fields]` + `method:` so every derivation is reproducible and re‑checkable from the row alone.
- **Source‑snapshot drift detection:** cache the extracted source tables (e.g. the Hollister/Rogers transcriptions) as fixtures; a test asserts the catalogue values still match the snapshot, catching silent edits to golden numbers.
- **Unit‑safety types:** consider lightweight unit tagging (km/s, AU, days) on ingest to kill the classic unit‑confusion bug class.
- **Coverage of the *reproduction* gauntlet** (sibling to the exclusion census): freeze how many rows are V0‑reproduced / real‑eph‑reproduced / anchor‑only, as a ratchet, so reproduction coverage can't silently regress.
- **Adversarial re‑extraction:** for high‑value rows, a second independent agent re‑reads the source scan and must agree digit‑for‑digit before the value is promoted to CROSS_VALIDATED (defends against single‑reader OCR error — the orbit‑15 checksum catch generalised).

## Self‑Review
Task 2 (physical invariants) is independent of the workflow and the highest‑value quick win (source‑free error detection). Tasks 3–7 deliver the provenance/tier/fidelity machinery. Falsification (Task 6) ensures the gates have teeth. Everything is golden‑clean: physical invariants need no source; tier checks *enforce* independence; fidelity tags prevent the cross‑fidelity bug; nothing loosens tolerances.
