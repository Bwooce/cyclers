# Forge Phase 6 — first novelty campaign: Jovian moon-system VILM sweep (task #172)

> **For agentic workers:** REQUIRED SUB-SKILL: `superpowers:executing-plans` or
> `superpowers:subagent-driven-development`. Checkbox steps; strict TDD (write
> failing test → run **red** → minimal impl → run **green** → commit). Work on
> `main` — **do NOT branch** (project rule). uv-managed venv (no pip). Lint/type
> gate before **every** commit: `uv run ruff check .` ·
> `uv run ruff format --check .` · `uv run mypy src tests`. Fast suite:
> `uv run pytest -m "not slow"`; corrector/ephemeris/DE440 geometry tests are
> `slow`. `export PATH="$HOME/.local/bin:$PATH"`.
>
> This plan is the executable expansion of the design note
> `docs/notes/2026-06-08-forge-phase6-discovery-design.md` (read it first,
> especially §2 the pipeline, §3 the novelty gate, §5 the family choice + honest
> yield, §6 the empty-region format). Where this plan and the design note diverge,
> the design note wins.
>
> **The Forge (Phases 0–5) is COMPLETE** — `discover_novel`, the adversarial
> panel, the human-review queue, the dedup matcher, and the n-body rung are all
> shipped (`docs/superpowers/plans/2026-06-03-the-forge-pipeline.md`, Completion
> notes). This plan does NOT rebuild them. It (a) wires an incremental VILM/
> Tisserand prune ahead of the moon-system scan, (b) adds the literature-check
> field to the review queue, (c) adds the empty-region report artefact, and
> (d) runs + reports the first forward sweep of the #76 moon-system space.

---

## Goal

Run the first FORWARD novelty-discovery campaign — finding cyclers NOT in the
literature — over the single most under-explored family with the best novelty
odds and lowest cost: the **planet-centric Jovian moon-system VILM space** (#76
substrate, never swept forward). Deliver the campaign wiring (VILM prune,
literature-check field, empty-region report), run the cheapest-decisive first
sweep (Galilean I-E-G, VILM-gated), and report the outcome honestly — expecting
mostly negatives, a rigorous bounded empty-set being a success.

**Out of scope:** the heliocentric E-M space (mined; Phase 4 already swept it),
the VEM space (#110/#120/#122 empirically refuted — 0 bend-feasible), CR3BP
Tier-2 moon systems, and any catalogue writeback (novelty never auto-promotes).

---

## Design facts verified live (these OVERRIDE any stale wording)

Verified 2026-06-08 against the live tree:

1. **`discover_novel` drives the E-M multi-arc space by default**
   (`src/cyclerfinder/data/discover_novel.py:em_multiarc_topologies`). It is
   topology-parametric: `discover_novel(topologies=...)` accepts any
   `TopologySpec` set. The moon-system sweep is a *new topology set + a centred
   ephemeris*, NOT a new loop. `scan_parallel` already accepts
   `ephem_model="astropy"`; the centred-circular path is `Ephemeris(model=
   "circular", center="Jupiter")` (#76 Phase 2). **VERIFY-FIRST** the live
   `build_epoch_branch_grid` / `scan_parallel` signatures support a centred
   ephemeris before wiring — they were written for heliocentric scans.
2. **The VILM module exists** (`src/cyclerfinder/search/vilm.py`):
   `vilm_dv_floor(moon_a, moon_b)`, `vilm_dv_min(moon_a, moon_b, via=...)`,
   `min_vinf_for_vilm(moon)`, `classify_vilm_leg(...)`. These are the prune-gate
   primitives. **Note the #76 deviation:** the admissible ΔV-floor is
   escape+capture, NOT the no-GA quadrature (a GA *reduces* ΔV, so no-GA is not a
   lower bound). Use `vilm_dv_floor` as the documented floor.
3. **The corrector is centre-agnostic** (`search/correct.py:ballistic_correct`
   takes `mu_central`; `_max_bend_deg` resolves moon codes via `SATELLITES`). The
   I-E-G chain CLOSES about Jupiter but is **bend-infeasible in the
   no-leveraging model** (#76 Phase 3 honest-risk: ~10 km/s V∞, 100-150° turns vs
   2-5° max-bend, recorded as strict xfail). The open question this campaign
   tests: does VILM gating + Laplace-resonance phasing surface a bend-feasible
   tour?
4. **The dedup bucket is `(model_assumption, primary, bodies, k)`** — a
   Jovicentric V∞ never compares to a heliocentric V∞ (#76 Phase 6,
   `signature_bucket_key`). The Jovian catalogue rows are family-seed null-numeric
   (no sourced V∞ multiset), so a closed Jovicentric cycler will read `novel`.
5. **The review queue is `data/review_queue.py`** (`append_review_entry`,
   `validate_review_entry`, `is_catalogue_source()` → `False` by contract). Adding
   a `literature_check` field is additive to the JSONL payload; it does NOT touch
   catalogue schema.
6. **CONCURRENCY:** sibling agents (#170 writes catalogue/validate/nbody/search;
   #171 writes GMAT docs) may be editing `src/`/`tests/`/`data/`. Every task
   touching a shared file (`discover_novel.py`, `review_queue.py`, `vilm.py`)
   carries a **VERIFY-FIRST** re-read of the live signature before editing. Use
   explicit-pathspec commits; never reset past a sibling commit.

---

## Phasing (independently shippable)

| Phase | Theme | Tasks | depends on |
|---|---|---|---|
| **1** | Incremental VILM/Tisserand prune gate for moon-pair legs | 1.0–1.2 (3) | #76 vilm/tisserand |
| **2** | Moon-system topology set + centred-scan wiring into the novelty loop | 2.0–2.2 (3) | Phase 1 |
| **3** | Literature-check field on the review queue + promotion guard | 3.0–3.1 (2) | review_queue |
| **4** | Empty-region report + method-capability re-sweep gate | 4.0a–4.2 (5) | Phase 2 |
| **5** | Run the first sweep + honest report | 5.0–5.2 (3) | Phases 1–4 |

Phase 1 (prune) and Phase 3 (queue field) are pure and shippable alone. Phase 4
carries both the empty-region artefact and the method-capability re-sweep gate
(4.0a/4.0b are pure and shippable alone). Phase 5 is the campaign run.
**Total: 16 tasks across 5 phases.**

---

## Phase 1 — Incremental VILM/Tisserand prune gate

The Ceriotti incremental-pruning import (design note §2; mining note
`docs/notes/2026-06-07-ceriotti-2010-mga-global-opt-mining.md` §7): promote the
post-hoc per-leg feasibility test into a *level gate* that prunes a moon-pair leg
BEFORE the next leg's grid is built. The per-leg criterion is NOT the ΔV — it is
the geometric/physics feasibility (Tisserand `linkable` + `_max_bend_deg` + the
VILM ΔV-floor). The objective is a sum of non-negative per-leg terms (Bellman), so
the gate never discards a feasible candidate.

### Task 1.0 — `moon_leg_admissible(moon_a, moon_b, vinf_kms, budget_kms)` predicate

**Files:** create `src/cyclerfinder/search/moon_prune.py`; test
`tests/search/test_moon_prune.py`.

A pure predicate combining the three sourced gates for ONE moon-pair leg:
- `vilm_dv_floor(moon_a, moon_b) <= budget_kms` (the VILM leveraging ΔV-floor;
  #76 deviation: escape+capture floor),
- `linkable(moon_a, moon_b, mu=PRIMARIES[primary])` truthy (Jovicentric contour
  intersection; VERIFY-FIRST the live `linkable` signature — it may take vinf
  ranges),
- `_max_bend_deg(vinf_kms, moon_b) > 0.0` (the flyby can turn at all).

Returns `(admissible: bool, reason: str)` so the prune is *recorded*, not silent
(the empty-region report consumes `reason`).

#### Failing test
```python
"""Phase 6 Phase 1: moon-pair leg admissibility prune (plan Phase 1 Task 1.0)."""
from __future__ import annotations

from cyclerfinder.core.satellites import PRIMARIES
from cyclerfinder.search.moon_prune import moon_leg_admissible


def test_europa_ganymede_leg_returns_reasoned_verdict() -> None:
    ok, reason = moon_leg_admissible(
        "Europa", "Ganymede", vinf_kms=4.0, budget_kms=10.0,
        primary="Jupiter",
    )
    assert isinstance(ok, bool)
    assert reason  # non-empty: the prune must record WHY


def test_zero_budget_prunes_on_vilm_floor() -> None:
    ok, reason = moon_leg_admissible(
        "Europa", "Ganymede", vinf_kms=4.0, budget_kms=0.0,
        primary="Jupiter",
    )
    assert ok is False
    assert "vilm" in reason.lower() or "floor" in reason.lower()
```

Run → **red** → impl the predicate (VERIFY-FIRST the `linkable` / `vilm_dv_floor`
/ `_max_bend_deg` live signatures) → **green**. Commit:
```
search/moon_prune: VILM+Tisserand+bend leg-admissibility predicate (Forge Phase 6)
```

### Task 1.1 — `prune_topology_legs(spec, vinf_seed, budget)` — incremental box gate

**Files:** `search/moon_prune.py`; test `tests/search/test_moon_prune.py` (extend).

Walk a `TopologySpec`'s sequence leg by leg; for each consecutive moon pair call
`moon_leg_admissible`. Return `(survives: bool, per_leg_reasons: tuple[str, ...])`.
A topology survives only if EVERY leg is admissible (back-pruning: one dead leg
kills the prefix). This is the Ceriotti level-gate applied to the topology before
any corrector run — the cheap pre-filter that keeps the scan off dead sequences.

#### Failing test
```python
def test_prune_topology_records_per_leg_reasons() -> None:
    from cyclerfinder.data.discover_novel import TopologySpec
    from cyclerfinder.search.moon_prune import prune_topology_legs
    spec = TopologySpec(
        sequence=("Io", "Europa", "Ganymede", "Io"),
        per_leg_revs=(0, 0, 0), per_leg_branch=("single",) * 3,
        period_k=1, period_sec=1.0e6, tof_seed_days=(3.5, 7.2), slack_leg=2,
    )
    survives, reasons = prune_topology_legs(
        spec, vinf_seed_kms=4.0, budget_kms=10.0, primary="Jupiter",
    )
    assert isinstance(survives, bool)
    assert len(reasons) == len(spec.sequence) - 1   # one per leg
```

Run → **red** → impl → **green**. Commit:
```
search/moon_prune: incremental per-leg box gate over a topology (Forge Phase 6)
```

### Task 1.2 — golden: the prune never discards the #76 closing I-E-G family

**Files:** test `tests/search/test_moon_prune.py` (extend).

GOLDEN-DISCIPLINE guard (no fabricated value): the I-E-G chain that #76 Phase 3
*demonstrably closes* about Jupiter must survive the prune at a budget at/above
its VILM floor — i.e. the prune does not throw away the one family we know
closes. This is the Bellman "never discard a feasible candidate" invariant made a
test. The EXPECTED side is structural (survives=True), not a numeric anchor.

```python
def test_prune_keeps_the_known_closing_ieg_family() -> None:
    from cyclerfinder.data.discover_novel import TopologySpec
    from cyclerfinder.search.moon_prune import prune_topology_legs
    spec = TopologySpec(
        sequence=("Io", "Europa", "Ganymede", "Io"),
        per_leg_revs=(0, 0, 0), per_leg_branch=("single",) * 3,
        period_k=1, period_sec=(1.769 + 3.551 + 7.155) * 86400.0,
        tof_seed_days=(3.5, 7.2), slack_leg=2,
    )
    # Budget set generously above the VILM floor: the known-closing family must
    # NOT be pruned (Bellman: never discard a feasible candidate).
    survives, _ = prune_topology_legs(
        spec, vinf_seed_kms=4.0, budget_kms=50.0, primary="Jupiter",
    )
    assert survives is True
```

> If this fails, the gate is over-pruning — FIX the gate, never relax the test to
> pass (the over-pruning would invalidate every later empty-set result).

Run → **red**/**green** → Commit:
```
test: prune keeps the known-closing I-E-G family (Bellman feasibility guard) (Forge Phase 6)
```

---

## Phase 2 — Moon-system topology set + centred-scan wiring

### Task 2.0 — `jovian_galilean_topologies()` topology set

**Files:** `src/cyclerfinder/data/discover_novel.py` (additive — a sibling of
`em_multiarc_topologies`); test `tests/data/test_discover_novel_moon.py`.

> **VERIFY-FIRST (concurrency #170):** re-read `discover_novel.py`'s live
> `TopologySpec` + `em_multiarc_topologies` before adding. Add ONLY a new
> function; do not modify the existing E-M topology set.

Return the Galilean topology set the campaign sweeps: I-E-G(-I) at Laplace-
resonance ToF seeds (Io 1.769 d, Europa 3.551 d, Ganymede 7.155 d synodic
spacings), over the converging Lambert branches (low/high × n_revs ∈ {0,1}). The
slack leg absorbs the period. `period_sec` seeds from the resonance multiple
(seed only — NON-GOLDEN, the corrector refines it).

#### Failing test
```python
"""Phase 6 Phase 2: Galilean topology set for the novelty loop (plan Task 2.0)."""
from __future__ import annotations

from cyclerfinder.data.discover_novel import jovian_galilean_topologies


def test_galilean_topologies_are_jovian_moon_sequences() -> None:
    specs = jovian_galilean_topologies()
    assert specs
    for s in specs:
        assert set(s.sequence) <= {"Io", "Europa", "Ganymede", "Callisto"}
        assert s.sequence[0] == s.sequence[-1]   # closed tour
```

Run → **red** → impl → **green**. Commit:
```
data/discover_novel: Galilean Jovian topology set (Forge Phase 6 first campaign)
```

### Task 2.1 — `discover_novel_moon(...)` — centred sweep entry, VILM-pruned

**Files:** `src/cyclerfinder/data/discover_novel.py`; test
`tests/data/test_discover_novel_moon.py` (extend).

> **VERIFY-FIRST:** re-read the live `discover_novel` body + `build_epoch_branch_
> grid` + `scan_parallel` signatures. Confirm a centred ephemeris
> (`Ephemeris(model="circular", center="Jupiter")`) and `mu_central` thread
> through the scan. If they do NOT (the scan was heliocentric-only), this task's
> scope grows to plumb `center`/`mu_central` through `scan_parallel` — STOP and
> report the gap rather than forcing it.

A sibling of `discover_novel` for the moon space: defaults to
`jovian_galilean_topologies()` and `center="Jupiter"`, applies
`prune_topology_legs` (Phase 1) to each topology BEFORE building its scan grid
(skipping pruned topologies, recording their reasons for the Phase 4 report), runs
the centred scan, and yields `NoveltyFinding` per closure — reusing
`evaluate_closure` verbatim (the bridge → signature → match → agreement →
gauntlet pipeline is centre-agnostic; the signature carries `primary="Jupiter"`
via `model_assumption`/bucket).

#### Failing test (fast — assert wiring, not a DE440 closure)
```python
def test_discover_novel_moon_prunes_then_scans(monkeypatch) -> None:
    # Assert the loop applies the prune and only scans survivors; stub the scan so
    # the test stays fast (the real DE440 sweep is the Phase 5 slow run).
    import cyclerfinder.data.discover_novel as dn
    seen_topos = []
    monkeypatch.setattr(dn, "scan_parallel", lambda grid, **kw: [])
    monkeypatch.setattr(
        dn, "build_epoch_branch_grid",
        lambda **kw: seen_topos.append(kw["sequence"]) or object(),
    )
    list(dn.discover_novel_moon(base_t0_sec=0.0, n_epochs=2, budget_kms=50.0))
    # At least the known-closing I-E-G family survived the prune and was scanned.
    assert any(set(s) <= {"Io", "Europa", "Ganymede", "Callisto"} for s in seen_topos)
```

Run → **red** → impl (reuse `evaluate_closure`; do NOT duplicate the per-candidate
pipeline) → **green**. Commit:
```
data/discover_novel: VILM-pruned centred Jovian novelty sweep (Forge Phase 6)
```

### Task 2.2 — sanity-e2e: a Jovian closure routes through the full pipeline (slow)

**Files:** test `tests/data/test_discover_novel_moon.py` (extend, `@pytest.mark.slow`).

The Phase-4-style sanity gate, ported to the moon centre: assert a *closed*
Jovian chain flows bridge → signature(`primary=Jupiter` bucket) → match(`novel`,
since the bucket is null-numeric) → Axis-A → gauntlet, and that a bend-INFEASIBLE
closure routes REJECTED (not SILVER) — exactly the #76 honest-risk family. This
proves the firewall holds on a non-heliocentric centre. **NON-GOLDEN** for the V∞
value (our computation).

> Expectation per the design note §5: the #76 I-E-G closure is bend-infeasible, so
> the realistic assertion is **`verdict.tier == REJECTED` for the no-leveraging
> closure** AND `match_outcome == "novel"`. If the VILM-gated sweep *does* surface
> a bend-feasible SILVER, assert it routes to the queue (do NOT auto-promote). If
> nothing closes at all, mark xfail with the empty-set reason — do NOT loosen tol.

Run `-m slow` → record the actual outcome in the task log (closed? bend-feasible?
tier?). Commit:
```
test: Jovian closure routes bridge->signature->match->gauntlet (non-golden, slow) (Forge Phase 6)
```

---

## Phase 3 — Literature-check field on the review queue

### Task 3.0 — `literature_check` field on `ReviewQueueEntry`

**Files:** `src/cyclerfinder/data/review_queue.py`; test
`tests/data/test_review_queue.py` (extend).

> **VERIFY-FIRST (concurrency):** re-read the live `ReviewQueueEntry` +
> `validate_review_entry` + `_normalise` before editing.

Add an optional `literature_check: dict | None = None` to the queue payload (a
non-catalogue artefact — additive, no catalogue-schema impact). Shape:
`{checked: bool, reviewer: str|None, date: str|None, sources_searched:
list[str], result: str|None}`. `validate_review_entry` accepts `None` (the
machine writes `None`; a human fills it). The default for a freshly-queued SILVER
candidate is `None` = "not yet checked".

#### Failing test
```python
def test_review_entry_carries_literature_check_default_none() -> None:
    from cyclerfinder.data.review_queue import ReviewQueueEntry, validate_review_entry
    e = ReviewQueueEntry(...)  # adapt to the live constructor
    assert e.literature_check is None
    validate_review_entry(e)   # None must validate (machine-written default)
```

Run → **red** → impl → **green**. Commit:
```
data/review_queue: literature_check field (V5 documented-review record) (Forge Phase 6)
```

### Task 3.1 — promotion guard: SILVER-novel ineligible until literature checked

**Files:** `src/cyclerfinder/data/review_queue.py` (a pure predicate
`is_promotion_eligible(entry) -> bool`); test `tests/data/test_review_queue.py`
(extend).

A pure predicate (NOT an auto-promoter — `is_catalogue_source()` stays `False`):
a SILVER-novel entry is promotion-eligible only when
`literature_check is not None and literature_check["checked"] is True and
literature_check["result"] == "no-match"`. This encodes spec §613's fourth
condition ("documented literature review returned nothing") as a recorded gate.
A literature *hit* makes it ineligible-as-novel (the human downgrades it to
`known-reproduction` per spec §612).

#### Failing test
```python
def test_unchecked_silver_is_not_promotion_eligible() -> None:
    from cyclerfinder.data.review_queue import ReviewQueueEntry, is_promotion_eligible
    e = ReviewQueueEntry(...)        # literature_check=None
    assert is_promotion_eligible(e) is False

def test_clean_literature_check_makes_eligible() -> None:
    from cyclerfinder.data.review_queue import ReviewQueueEntry, is_promotion_eligible
    e = ReviewQueueEntry(..., literature_check={
        "checked": True, "reviewer": "human", "date": "2026-06-08",
        "sources_searched": ["ADS", "Jones AAS 17-577"], "result": "no-match"})
    assert is_promotion_eligible(e) is True
```

Run → **red** → impl → **green**. Commit:
```
data/review_queue: is_promotion_eligible gate (no promotion without literature check) (Forge Phase 6)
```

---

## Phase 4 — Empty-region report artefact

> Phase 4 carries SIX tasks: the empty-region artefact (4.0/4.1/4.2) AND the
> method-capability re-sweep gate (4.0a the partial-order registry; 4.0b the
> `should_sweep` gate), per design note §6a/§6b. The gate is what makes
> `empty_regions.jsonl` a *re-sweepable* record rather than a permanent foreclosure
> — without it, a recorded "empty" silently locks out every future, more-capable
> method (the #163-reopens-#137 lesson). **Total tasks now: 16 across 5 phases.**

### Task 4.0a — method-capability partial order (`MethodCapability` + `subsumes`)

**Files:** create `src/cyclerfinder/data/method_capability.py`; test
`tests/data/test_method_capability.py`.

The capability registry (design note §6a). A small enum/registry of capability
tags and a `subsumes(a, b) -> bool` partial-order predicate. `MethodCapability` is
a frozen dataclass: `genome: str`, `corrector: str`, `capability_tags:
frozenset[str]`, `git_sha: str`. The partial order is defined over the tags by an
explicit edge set (the §6b ordering), e.g.:

- `multi-arc` ⊐ `single-arc`,
- `n-body` ⊐ `patched-conic`,
- `powered`/`low-thrust` ⊐ `ballistic`,
- `one-dsm-per-leg` ⊐ `single-arc`,
- `broken-plane` ⊐ `coplanar`.

`subsumes(a, b)` is True iff every tag of `b` is reached by `a` under the partial
order (`b`'s envelope ⊆ `a`'s). It MUST be reflexive (`subsumes(a, a) is True`)
and MUST return `False` for incomparable methods (neither contains the other).

> **NON-GOLDEN / sourced-discipline:** the edge set is a *design decision*
> transcribed from §6b, not a computed value. Keep the edges in one named
> constant (`_CAPABILITY_EDGES`) so the partial order is auditable in one place.

#### Failing test
```python
"""Phase 6 Phase 4: method-capability partial order (plan Task 4.0a)."""
from __future__ import annotations

from cyclerfinder.data.method_capability import MethodCapability, subsumes

SINGLE = MethodCapability(
    genome="single-ellipse free-return", corrector="ballistic_correct",
    capability_tags=frozenset({"ballistic", "patched-conic", "single-arc",
                               "coplanar"}), git_sha="aaa")
MULTI = MethodCapability(
    genome="two-arc free-return chain", corrector="ballistic_correct",
    capability_tags=frozenset({"ballistic", "patched-conic", "multi-arc",
                               "coplanar"}), git_sha="bbb")
BROKEN = MethodCapability(
    genome="single-ellipse free-return", corrector="ballistic_correct",
    capability_tags=frozenset({"ballistic", "patched-conic", "single-arc",
                               "broken-plane"}), git_sha="ccc")


def test_subsumes_is_reflexive() -> None:
    assert subsumes(SINGLE, SINGLE) is True


def test_multi_arc_subsumes_single_ellipse() -> None:
    assert subsumes(MULTI, SINGLE) is True       # #163 ⊐ #137
    assert subsumes(SINGLE, MULTI) is False       # weaker does not subsume


def test_incomparable_methods_do_not_subsume() -> None:
    # coplanar-multi-arc vs broken-plane-single-arc: neither envelope ⊆ the other
    assert subsumes(MULTI, BROKEN) is False
    assert subsumes(BROKEN, MULTI) is False
```

Run → **red** → impl the partial order (one named edge constant) → **green**.
Commit:
```
data/method_capability: capability partial order + subsumes predicate (Forge Phase 6)
```

### Task 4.0b — the re-sweep gate `should_sweep(region, method, registry) -> bool`

**Files:** `src/cyclerfinder/data/method_capability.py` (or
`src/cyclerfinder/data/empty_regions.py` once 4.0 lands — keep it next to the
`load_empty_regions` reader); test `tests/data/test_method_capability.py`
(extend, or `tests/data/test_empty_regions.py`).

The core gate (design note §6b). Given the proposed `region` (its
region_id / box), the proposed `method: MethodCapability`, and the loaded
`empty_regions.jsonl` registry, decide whether to sweep:

```text
should_sweep(region, method, registry):
    priors = [r for r in registry if r covers the same region/box]
    if any(subsumes(prior.method_capability, method) for prior in priors):
        return False        # a prior, ≥-capable method already emptied it → SKIP
    return True             # no prior subsumes the proposed method → RE-SWEEP
                            # (new/more-capable OR incomparable both re-sweep)
```

The skip criterion is capability-subsumption, NOT region-match alone: a prior
empty over the same box only skips if its method **subsumes** the proposed one.

#### Failing test
```python
def test_weaker_method_skips_region_a_stronger_method_emptied() -> None:
    # Prior: a STRONGER (multi-arc) sweep emptied region R.
    # Proposed: a WEAKER (single-ellipse) method → nothing new to learn → SKIP.
    from cyclerfinder.data.method_capability import should_sweep
    registry = [_empty_record(region_id="R", method=MULTI)]
    assert should_sweep(region_id="R", method=SINGLE, registry=registry) is False


def test_stronger_method_re_sweeps_region_a_weaker_method_emptied() -> None:
    # Prior: a WEAKER (single-ellipse) sweep emptied R; #163-reopens-#137.
    # Proposed: a STRONGER (multi-arc) method → re-sweep the reopened ground.
    from cyclerfinder.data.method_capability import should_sweep
    registry = [_empty_record(region_id="R", method=SINGLE)]
    assert should_sweep(region_id="R", method=MULTI, registry=registry) is True


def test_incomparable_method_re_sweeps() -> None:
    # Prior emptied R with coplanar-multi-arc; proposed is broken-plane (incomparable).
    from cyclerfinder.data.method_capability import should_sweep
    registry = [_empty_record(region_id="R", method=MULTI)]
    assert should_sweep(region_id="R", method=BROKEN, registry=registry) is True


def test_no_prior_for_region_re_sweeps() -> None:
    from cyclerfinder.data.method_capability import should_sweep
    assert should_sweep(region_id="UNSWEPT", method=SINGLE, registry=[]) is True
```

(`_empty_record` is a small test helper building an `EmptyRegionReport`-shaped
object carrying `region_id` + `method_capability`.)

Run → **red** → impl → **green**. Commit:
```
data/method_capability: should_sweep capability-subsumption re-sweep gate (Forge Phase 6)
```

### Task 4.0 — `EmptyRegionReport` dataclass + serialiser

**Files:** create `src/cyclerfinder/data/empty_regions.py`; test
`tests/data/test_empty_regions.py`.

The first-class negative (design note §6). A frozen dataclass with the required
fields: `region_id`, `family`, `centre`, `topologies`, `method_capability`
(the §6a descriptor — see Task 4.0a), `search_extent`
(n_epochs / span_days / n_topologies / points_total / ephem_model / center),
`prune_gates`, `result` (closed / distinct_families / bend_feasible /
best_max_vinf_kms / vinf_floor_target_kms / gap_kms), `verdict`,
`interpretation`, `source_anchors`, `run` (date / host / cores / git_sha /
wall_s). Plus `append_empty_region(report, path)` → `data/empty_regions.jsonl`
(JSONL, mirroring `append_review_entry`).

> **The bar for a negative to count (test it):** a report is INVALID if
> `search_extent.points_total` is absent/zero (an unbounded negative is a
> silently-dropped negative), `prune_gates` is empty (can't tell if the empty
> set is an over-pruning artefact), or `method_capability` is absent/empty (an
> *unconditional* "empty" claim — "empty" is never unconditional, design note
> §6a). `validate_empty_region` enforces all three.

#### Failing test
```python
"""Phase 6 Phase 4: empty-region report is bounded + reproducible (plan Task 4.0)."""
from __future__ import annotations

import pytest

from cyclerfinder.data.empty_regions import EmptyRegionReport, validate_empty_region


def test_empty_region_requires_bounded_search_extent() -> None:
    r = EmptyRegionReport(region_id="x", family="f", centre="Jupiter",
        topologies=(), search_extent={"points_total": 0}, prune_gates=("vilm",),
        result={}, verdict="EMPTY", interpretation="", source_anchors="", run={})
    with pytest.raises(ValueError):
        validate_empty_region(r)   # points_total == 0 -> unbounded -> invalid


def test_empty_region_requires_prune_gates() -> None:
    r = EmptyRegionReport(region_id="x", family="f", centre="Jupiter",
        topologies=(), search_extent={"points_total": 2816}, prune_gates=(),
        result={}, verdict="EMPTY", interpretation="", source_anchors="", run={})
    with pytest.raises(ValueError):
        validate_empty_region(r)   # no prune gates -> can't bound over-pruning
```

Run → **red** → impl → **green**. Commit:
```
data/empty_regions: bounded+reproducible empty-region report artefact (Forge Phase 6)
```

### Task 4.1 — `append_empty_region` round-trips JSONL

**Files:** `src/cyclerfinder/data/empty_regions.py`; test
`tests/data/test_empty_regions.py` (extend).

Mirror `append_review_entry`/`load_review_queue`: append one valid report,
re-load, assert field equality. `is_catalogue_source()`-equivalent contract: the
empty-region log NEVER feeds the catalogue (it is a negative-result audit trail).

Run → **red**/**green** → Commit:
```
data/empty_regions: JSONL append + load round-trip (Forge Phase 6)
```

### Task 4.2 — orchestrator emits an empty-region report when a sweep is barren

**Files:** `scripts/forge_phase6_moon_run.py` (new — sibling of
`scripts/forge_novelty_run.py`); test
`tests/scripts/test_forge_phase6_moon_run.py` (fast, stubbed scan).

> **VERIFY-FIRST:** read `scripts/forge_novelty_run.py` for the fan-out + panel +
> queue pattern; reuse it, do not reinvent. Swap the topology set for
> `jovian_galilean_topologies` + `discover_novel_moon`.

The orchestrator: BEFORE sweeping, call `should_sweep(region, method, registry)`
(Task 4.0b) against the loaded `empty_regions.jsonl` — skip the region iff a prior
≥-capable method already emptied it (the proposed method here is the
single-ellipse no-leveraging `MethodCapability`; record the skip reason if
skipped). If swept: run `discover_novel_moon`, fan SILVER survivors through the
adversarial panel (`verify/adversarial.py`) + queue them (with
`literature_check=None`), and — if the sweep yields zero promotable candidates —
emit an `EmptyRegionReport` carrying the `method_capability` descriptor + the
actual search extent + best-achieved V∞ + gap. Assert (stubbed) that (a) a barren
sweep writes a valid method-versioned empty-region report, and (b) re-running with
the *same* method against that record skips (should_sweep → False).

Run → **red**/**green** → Commit:
```
scripts/forge_phase6_moon_run: Jovian VILM sweep orchestrator + re-sweep gate + empty-region emit (Forge Phase 6)
```

---

## Phase 5 — Run the first sweep + honest report

### Task 5.0 — cheapest-decisive first sweep (slow, real DE440)

**Files:** none (a run); record to a results note
`docs/notes/2026-06-08-forge-phase6-jovian-sweep-results.md`.

Run `scripts/forge_phase6_moon_run.py` on the Galilean I-E-G VILM-gated space:
modest extent first (e.g. 64 epochs × the Galilean topology set, `center=Jupiter`,
budget at/above the VILM floor) — the cheapest grid that decisively answers "does
the VILM-gated sweep surface a bend-feasible Jovian tour the #76 no-leveraging
closure could not?". Record: closed count, distinct families, bend-feasible
count, best max-V∞, the V∞ floor gap, and every SILVER survivor's panel verdict.

> **HONEST-RISK (binding):** if nothing closes bend-feasible (the expected
> outcome per design note §5), that is a SUCCESS — emit the empty-region report
> (carrying the single-ellipse no-leveraging `method_capability` descriptor, so a
> later multi-arc/n-body/low-thrust method re-sweeps per §6b) and STOP; do NOT
> loosen `tol_kms`, the budget, or the bend cap to manufacture a survivor. If a
> SILVER survives, run the n-body ARTIFACT rung (Task 5.1) before any human-facing
> claim.

### Task 5.1 — independent n-body cross-check of any SILVER survivor (slow)

**Files:** results note (extend).

For EACH SILVER survivor, run the n-body rung (the SILVER-ARTIFACT template,
design note §4: REBOUND/IAS15 about Jupiter over the shared DE440 BSP, body set
justified by a sensitivity test — which moons + the Sun move the rung metric by
more than tolerance). Record the rung verdict (CONFIRM / ARTIFACT / diverge).
**No promotion** — the candidate stays SILVER, `promoted=False`,
`literature_check=None`; the human decides.

> If there are ZERO SILVER survivors, this task is a no-op — record "no SILVER
> survivors; n-body rung not exercised" and proceed to the report.

### Task 5.2 — honest campaign report

**Files:** results note (finalise); `data/OUTSTANDING.md` (append a Phase 6
campaign entry — VERIFY-FIRST, append-only, do not edit sibling entries).

Write up: the family swept + why (best odds, lowest cost, fresh substrate), the
search extent, the outcome (almost certainly mostly negatives / empty), the
empty-region report path, any SILVER + its n-body rung, and the sharpened
hypothesis (the #110-style "a thorough-sweep-finds-nothing result is real
science"). State the honest yield plainly: a discovery was NOT expected on the
first sweep, and a rigorous bounded negative is the success criterion. Commit:
```
docs: Forge Phase 6 first campaign — Jovian VILM sweep results (task #172)
```

---

## Self-review

- **Reuses Phases 0–5.** No new finding loop, panel, queue, dedup, or n-body rung
  — Phase 6 is a topology set + a VILM prune + a literature-check field + a
  method-versioned empty-region artefact (with a capability-subsumption re-sweep
  gate) + one run.
- **"Empty" is never unconditional.** Every empty-region record carries a
  method-capability descriptor; the `should_sweep` gate skips ONLY when a prior
  ≥-capable method emptied the region, and re-sweeps for any more-capable OR
  incomparable method — the #163-reopens-#137 lesson encoded as a gate, so a
  recorded negative never permanently forecloses a future capability.
- **Golden-clean.** The prune gates are sourced (VILM Endgame Tables); no computed
  V∞ is a golden EXPECTED; novelty caps at SILVER; the empty-set guard forbids
  unbounded negatives; the literature gate forbids promotion without a documented
  review.
- **Honest odds, built in.** The plan *expects* mostly negatives and makes the
  empty region a first-class deliverable; the slow-run tasks carry an explicit
  "do NOT loosen tolerance to manufacture a survivor" rule.
- **Risk.** (1) `scan_parallel` may be heliocentric-only — Task 2.1 STOP/reports
  rather than forcing it. (2) A false "novel" from the null-numeric Jovian bucket
  — mitigated by the literature gate (Phase 3) + the n-body rung (Task 5.1).
  (3) The moon space may also be empty of bend-feasible novelty (the #76
  honest-risk generalises) — mitigated by the empty-region report making that a
  rigorous bounded negative, not a dead end.
```
