# SnLm Multi-Rev Rediscovery Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Re-model the S1L1 cycler with its correct Earth-to-Earth resonant topology and rediscover it as a *computed* V0→V1 cycler (Phase 1), then turn on multi-rev enumeration so the whole SnLm family (the 202-entry `MULTI_ENCOUNTER_SEQUENCE` bucket) becomes reachable by `discover` (Phase 2).

**Architecture:** The Lambert solver, `construct_cycler`, `enumerate_cells`/`feasible_cells`, and `optimise_cell_idealized` *already* support multi-rev legs (`n_revs >= 1`, `branch in {"low","high"}`). What is missing is: (1) the correct S1L1 encounter topology, which is unknown/provisional and must be **derived** by searching candidate structures that close to the published V∞ anchors (5.65 km/s Earth, 3.05 km/s Mars) — never invented; (2) the rediscovery harness/loader wiring that currently hard-codes direct (0-rev) 3-encounter cells; (3) plumbing `n_max`/`branch_set` through `discover` for the family sweep. EXPECTED values must trace to a published source or be honestly labelled as *computed* (not golden).

**Tech Stack:** Python 3.11, numpy-only production Lambert path, lamberthub (dev) for V1 crosscheck, pytest + xdist, uv-managed venv (no pip), ruff + mypy via pre-commit.

**Provenance contract (binding):** S1L1's S1/L1 interval ToFs are **not published**. This plan treats them as *computed* by our own multi-rev Lambert solver constrained to the published 5.65/3.05 km/s V∞ and the 2-synodic period, validated as a V0→V1 rediscovery + lamberthub V1 crosscheck — **NOT** as a golden gate with a sourced EXPECTED value. The published 5.65/3.05 V∞ anchors and the ~154-d outbound ToF ARE sourced and remain the only golden inputs. See `[[s1l1-nomenclature]]` and the `feedback_golden_tests_sourced_only` memory.

---

## Background: what already exists (read before starting)

- `src/cyclerfinder/search/sequence.py:57` — `Cell(bodies, sequence, period_k, per_leg_revs, per_leg_branch)`. The `id` property encodes revs/branches (e.g. `EM|E-M-E-E|k2|r010|blsl`).
- `src/cyclerfinder/search/sequence.py:145` — `enumerate_cells(body_set, l_max, k_max, n_max, branch_set=("single",))`. Already yields multi-rev cells when `n_max >= 1` and `branch_set` includes `"low"`/`"high"`.
- `src/cyclerfinder/search/sequence.py:340` — `feasible_cells(body_set, l_max, k_max, n_max, vinf_cap, ephem=None, branch_set=("single",))`.
- `src/cyclerfinder/search/construct.py:40` — `construct_cycler(sequence, encounter_times_sec, ephem, mu_sun=..., max_revs_per_leg=None, branch_per_leg=None)`. Full multi-rev support; raises `ValueError` if no Lambert solution matches the requested branch.
- `src/cyclerfinder/search/optimize.py:960` — `optimise_cell_idealized(cell, ephem, *, vinf_cap, n_starts=5, seed=0, use_de=True, rp_factors=None, target_period_sec=None, warm_starts=None) -> OptimisationResult`. The optimiser does NOT mutate the cell's discrete structure — only the interior encounter epochs. Returns `OptimisationResult(cell, best_cycler, best_score, closure_residual_kms, optimiser_history, converged, constraints_satisfied)`; `best_score` is a `Score` with `.max_vinf_kms` (`src/cyclerfinder/model/score.py:113`).
- `src/cyclerfinder/search/optimize.py:229` — `_target_period_sec(cell)` = `period_k * synodic_period_days(bodies[0], bodies[1]) * SECONDS_PER_DAY` (2-body correct).
- `tests/_catalogue_loader.py` — `classify_row`, `ExclusionReason`, `CatalogueEntry`, `load_constructible_entries`. S1L1 is currently `MULTI_ENCOUNTER_SEQUENCE` (excluded) because `_is_two_body_alternation` rejects its `sequence_canonical`.
- `tests/test_catalogue_rediscovery.py` — the parametrised V0→V1 gauntlet; `_build_cell_from_entry` hard-codes `per_leg_revs=(0,0)`, `per_leg_branch=("single","single")`. `EXPECTED_COVERAGE` (line 267) is a frozen census ratchet: `MULTI_ENCOUNTER_SEQUENCE: 202, NON_HELIOCENTRIC: 6, MISSING_VINF: 5, CONSTRUCTIBLE: 2, NOT_TWO_BODY: 2, MISSING_PERIOD: 2`.
- `data/catalogue.yaml:452` — the `s1l1-2syn-em-cpom` entry. Its `trajectory.segments` currently carries `out-em` (E→M, 154 d, 0-rev), a mis-modelled `ret-me` (M→E, null), and `loop-ee` (E→E, null). `data_gaps` already flag the topology as provisional.

---

## PHASE 1 — Re-model S1L1 as a computed SnLm exemplar

### Task 1: Characterise the S1L1 topology (exploratory derivation)

**This task DERIVES the unknown topology. Its output is a recorded finding that later tasks pin.** It is a one-shot characterisation script, not a permanent test.

**Files:**
- Create: `scripts/characterise_s1l1.py`

- [ ] **Step 1: Write the characterisation script**

```python
"""One-shot: find which 2-synodic E-M multi-rev topology closes to S1L1's
published V-infinity anchors (5.65 km/s Earth, 3.05 km/s Mars).

NOT a test. Records the winning (sequence, per_leg_revs, per_leg_branch)
so the permanent rediscovery test (Task 2) and the catalogue writeback
(Task 3) can pin it. The S1/L1 interval ToFs it prints are COMPUTED, not
sourced — see the plan's provenance contract.
"""

from __future__ import annotations

import itertools

import numpy as np

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.optimize import optimise_cell_idealized
from cyclerfinder.search.sequence import Cell

# Sourced anchors (golden inputs): spec.md §9 / catalogue s1l1-2syn-em-cpom.
VINF_E_TARGET = 5.65
VINF_M_TARGET = 3.05
VINF_TOL = 0.3            # same band as the rediscovery gauntlet
VINF_CAP = max(VINF_E_TARGET, VINF_M_TARGET) + 2.5

# Candidate topologies for a 2-synodic cycler with one Mars flyby (target
# of opportunity) on the outbound arc plus Earth-to-Earth resonant
# interval(s). Mars appears once; Earth brackets the resonant intervals.
CANDIDATE_SEQUENCES = [
    ("E", "M", "E"),        # single Earth-Earth interval bracketing Mars
    ("E", "M", "E", "E"),   # S1 (E->M->E) + L1 (E->E)
]
N_MAX = 2                   # try 0,1,2 revs per leg
BRANCHES = ("low", "high")

eph = Ephemeris(model="circular")


def vinf_by_body(result) -> dict[str, float]:
    out: dict[str, float] = {}
    for enc in result.best_cycler.encounters:
        m = max(float(np.linalg.norm(enc.vinf_in)), float(np.linalg.norm(enc.vinf_out)))
        out[enc.body] = max(out.get(enc.body, 0.0), m)
    return out


def main() -> None:
    hits = []
    for seq in CANDIDATE_SEQUENCES:
        n_legs = len(seq) - 1
        for revs in itertools.product(range(N_MAX + 1), repeat=n_legs):
            # branch choices: 0-rev -> "single"; >=1 -> low/high
            per_leg_branch_choices = [
                ("single",) if r == 0 else BRANCHES for r in revs
            ]
            for branches in itertools.product(*per_leg_branch_choices):
                cell = Cell(
                    bodies=("E", "M"),
                    sequence=seq,
                    period_k=2,
                    per_leg_revs=tuple(revs),
                    per_leg_branch=tuple(branches),
                )
                try:
                    res = optimise_cell_idealized(
                        cell, eph, vinf_cap=VINF_CAP, n_starts=5, seed=0, use_de=True
                    )
                except Exception as exc:  # noqa: BLE001 - exploratory script
                    print(f"  {cell.id}: raised {type(exc).__name__}")
                    continue
                v = vinf_by_body(res)
                ok = (
                    res.constraints_satisfied
                    and abs(v.get("E", 1e9) - VINF_E_TARGET) < VINF_TOL
                    and abs(v.get("M", 1e9) - VINF_M_TARGET) < VINF_TOL
                )
                tag = "HIT" if ok else "   "
                print(
                    f"{tag} {cell.id}: E={v.get('E', float('nan')):.3f} "
                    f"M={v.get('M', float('nan')):.3f} "
                    f"resid={res.closure_residual_kms:.4f} "
                    f"feasible={res.constraints_satisfied}"
                )
                if ok:
                    tofs_days = [
                        (leg.t_arrive - leg.t_depart) / 86400.0
                        for leg in res.best_cycler.legs
                    ]
                    hits.append((cell.id, tuple(revs), tuple(branches), tofs_days))
    print("\n=== HITS ===")
    for cid, revs, branches, tofs in hits:
        print(f"{cid}  revs={revs} branches={branches} leg_tofs_days={tofs}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run the characterisation**

Run: `uv run python scripts/characterise_s1l1.py`
Expected: a table of candidate topologies with their rediscovered E/M V∞. **Record every line tagged `HIT`** — these are the topologies whose computed V∞ matches the sourced 5.65/3.05 anchors within ±0.3 km/s. If there are zero hits, STOP and report: the circular-coplanar model may not host S1L1 (as with Aldrin — see `tests/test_catalogue_rediscovery.py:99` `EXPECTED_SKIPS`), in which case Phase 1's validation gate becomes a documented xfail with the honest "model limitation" reason rather than a passing test, and you proceed to Phase 2.

- [ ] **Step 3: Record the winning topology**

Pick the HIT with the lowest `closure_residual_kms`. Write its `(sequence, per_leg_revs, per_leg_branch)` and computed `leg_tofs_days` into a comment block at the top of `scripts/characterise_s1l1.py` headed `# FINDING (2026-06-02):`. This is the single value later tasks depend on; capturing it in-repo keeps the plan reproducible.

- [ ] **Step 4: Commit**

```bash
git add scripts/characterise_s1l1.py
git commit -m "search: S1L1 topology characterisation script + recorded finding"
```

---

### Task 2: Permanent S1L1 computed-rediscovery test

**Files:**
- Create: `tests/search/test_s1l1_rediscovery.py`

> Replace `WINNING_SEQUENCE`, `WINNING_REVS`, `WINNING_BRANCHES` below with the tuple recorded in Task 1 Step 3. The numbers in the asserts are the **sourced** anchors (5.65/3.05) — those do not change.

- [ ] **Step 1: Write the failing test**

```python
"""S1L1 computed V0->V1 rediscovery (NOT a golden gate).

The S1/L1 Earth-to-Earth resonant-interval ToFs are computed by our own
multi-rev Lambert solver constrained to the published 5.65/3.05 km/s
V-infinity anchors and the 2-synodic period; only the V-infinity anchors
are sourced. See docs/superpowers/plans/2026-06-02-snlm-multirev-rediscovery.md.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.optimize import optimise_cell_idealized
from cyclerfinder.search.sequence import Cell

# Topology recorded by scripts/characterise_s1l1.py (Task 1).
WINNING_SEQUENCE = ("E", "M", "E", "E")     # <-- pin from Task 1 finding
WINNING_REVS = (0, 0, 1)                      # <-- pin from Task 1 finding
WINNING_BRANCHES = ("single", "single", "low")  # <-- pin from Task 1 finding

VINF_E_TARGET = 5.65   # sourced: spec.md §9 / catalogue s1l1-2syn-em-cpom
VINF_M_TARGET = 3.05   # sourced
VINF_TOL = 0.3
VINF_CAP = max(VINF_E_TARGET, VINF_M_TARGET) + 2.5


def _vinf_by_body(result) -> dict[str, float]:
    out: dict[str, float] = {}
    for enc in result.best_cycler.encounters:
        m = max(float(np.linalg.norm(enc.vinf_in)), float(np.linalg.norm(enc.vinf_out)))
        out[enc.body] = max(out.get(enc.body, 0.0), m)
    return out


@pytest.mark.slow
def test_s1l1_rediscovers_published_vinf_anchors() -> None:
    cell = Cell(
        bodies=("E", "M"),
        sequence=WINNING_SEQUENCE,
        period_k=2,
        per_leg_revs=WINNING_REVS,
        per_leg_branch=WINNING_BRANCHES,
    )
    eph = Ephemeris(model="circular")
    result = optimise_cell_idealized(
        cell, eph, vinf_cap=VINF_CAP, n_starts=5, seed=0, use_de=True
    )
    assert result.constraints_satisfied
    v = _vinf_by_body(result)
    assert abs(v["E"] - VINF_E_TARGET) < VINF_TOL
    assert abs(v["M"] - VINF_M_TARGET) < VINF_TOL
```

- [ ] **Step 2: Run the test to verify it passes (topology already found)**

Run: `uv run pytest tests/search/test_s1l1_rediscovery.py -v`
Expected: PASS (the topology was selected in Task 1 precisely because it closes). If it FAILS, the pinned tuple is wrong — re-check the Task 1 finding.

- [ ] **Step 3: Commit**

```bash
git add tests/search/test_s1l1_rediscovery.py
git commit -m "tests: S1L1 computed multi-rev V1 rediscovery against sourced 5.65/3.05 anchors"
```

---

### Task 3: Write the derived topology back into the catalogue

**Files:**
- Modify: `data/catalogue.yaml` (the `s1l1-2syn-em-cpom` `trajectory.segments`, `maneuvers`, and `data_gaps`)

- [ ] **Step 1: Replace the segment skeleton with the derived topology**

Edit `data/catalogue.yaml` `s1l1-2syn-em-cpom`: replace the three provisional segments (`out-em`, `ret-me`, `loop-ee`) with segments matching `WINNING_SEQUENCE`/`WINNING_REVS`/`WINNING_BRANCHES`. Set each multi-rev segment's `tof_days`/`n_revs`/`branch` to the **computed** values from the Task 1 finding, and mark each with `provenance: computed` plus a note pointing at the rediscovery test. The outbound `out-em` keeps its sourced `tof_days: 154`. Keep `maneuvers` consistent with the new boundaries (ballistic flybys, `dv_kms: 0.0`).

Example shape (substitute the Task 1 numbers):

```yaml
    segments:
      - id: "out-em"
        from: "E"
        to: "M"
        traj_type: "keplerian-arc"
        tof_days: 154            # sourced (Spreen 2020 / Mistree-Kornfeld 2019)
        n_revs: 0
        branch: "single"
      - id: "s1-me"
        from: "M"
        to: "E"
        traj_type: "keplerian-arc"
        tof_days: <COMPUTED>     # provenance: computed (multi-rev Lambert)
        n_revs: <COMPUTED>
        branch: "<COMPUTED>"
        note: "Second half of the S1 Earth-to-Earth resonant interval after the Mars flyby. ToF computed by our multi-rev Lambert solver constrained to the published 5.65/3.05 km/s V-infinity; NOT a sourced value. See tests/search/test_s1l1_rediscovery.py."
      - id: "l1-ee"
        from: "E"
        to: "E"
        traj_type: "keplerian-arc"
        tof_days: <COMPUTED>     # provenance: computed
        n_revs: <COMPUTED>
        branch: "<COMPUTED>"
        note: "The L1 long Earth-to-Earth resonant interval. ToF computed, not sourced (same constraint as s1-me)."
```

- [ ] **Step 2: Update `data_gaps` to reflect resolution-by-computation**

Change the `trajectory.segments` topology gap `kind` from `"uncertain"` to `"derive"` and note that the topology is now derived (computed) and validated by `tests/search/test_s1l1_rediscovery.py`, while remaining unsourced (a published S1/L1 breakdown would upgrade it from computed to sourced).

- [ ] **Step 3: Run the loader + verify tests to confirm the YAML still parses**

Run: `uv run pytest tests/data/test_catalogue_loader.py tests/verify/test_real_closure.py -q`
Expected: PASS (no schema breakage; S1L1 still skipped/xfail in real-closure until Phase 1 Task 4).

- [ ] **Step 4: Commit**

```bash
git add data/catalogue.yaml
git commit -m "catalogue: replace S1L1 provisional skeleton with derived E->M + S1/L1 resonant topology (computed)"
```

---

### Task 4: Retire the S1L1 real-closure xfail to a computed-closure assertion

**Files:**
- Modify: `tests/verify/test_real_closure.py` (the Gate-2 test at line ~266)
- Modify: `src/cyclerfinder/verify/real_closure.py` (`EXPECTED_SKIPS["s1l1-2syn-em-cpom"]` at line ~222)

- [ ] **Step 1: Decide the gate based on Task 1 outcome**

If Task 1 produced HITs (S1L1 closes in circular-coplanar): rewrite the Gate-2 test to assert closure on the now-complete derived topology and remove the xfail marker. If Task 1 produced NO hits: keep the xfail but rewrite its reason to the honest "circular-coplanar model limitation" wording (mirroring the Aldrin `EXPECTED_SKIPS` reason at `tests/test_catalogue_rediscovery.py:99`), and note the idealised rediscovery test (Task 2) is the binding computed gate instead.

- [ ] **Step 2: Run the affected tests**

Run: `uv run pytest tests/verify/test_real_closure.py -q`
Expected: PASS (either a newly-green Gate-2 or a correctly-reasoned xfail).

- [ ] **Step 3: Commit**

```bash
git add tests/verify/test_real_closure.py src/cyclerfinder/verify/real_closure.py
git commit -m "verify: S1L1 closes on derived topology (or honest model-limit xfail); retire stale data-gap reason"
```

---

### Task 5: Update the coverage-census ratchet

**Files:**
- Modify: `tests/test_catalogue_rediscovery.py:267` (`EXPECTED_COVERAGE`)
- Possibly modify: `tests/_catalogue_loader.py` (`classify_row`) if S1L1 is to be admitted to the gauntlet

- [ ] **Step 1: Decide whether S1L1 enters the parametrised gauntlet**

The parametrised gauntlet (`_build_cell_from_entry`) only builds direct 3-encounter cells. Admitting S1L1 there would require teaching the loader/harness multi-rev cell construction. For Phase 1, keep S1L1 covered by its dedicated test (Task 2) and leave it classified `MULTI_ENCOUNTER_SEQUENCE` — but verify the census count is unchanged (S1L1 was already in that bucket).

- [ ] **Step 2: Run the census ratchet**

Run: `uv run pytest tests/test_catalogue_rediscovery.py -k census -v`
Expected: PASS unchanged (`MULTI_ENCOUNTER_SEQUENCE: 202` etc.). If the catalogue edit in Task 3 changed `sequence_canonical`, update `EXPECTED_COVERAGE` in the same commit and confirm the shift is intended.

- [ ] **Step 3: Commit (only if a ratchet update was needed)**

```bash
git add tests/test_catalogue_rediscovery.py tests/_catalogue_loader.py
git commit -m "tests: keep coverage census consistent after S1L1 re-model"
```

---

## PHASE 2 — SnLm family sweep through `discover`

### Task 6: Plumb `n_max` / `branch_set` through `discover`

**Files:**
- Modify: `src/cyclerfinder/data/discover.py` (`discover(...)` signature + its `feasible_cells(...)` call)
- Test: `tests/data/test_discover.py`

- [ ] **Step 1: Write the failing test**

```python
def test_discover_accepts_multirev_params(tmp_path):
    from cyclerfinder.data.discover import discover
    # A bounded multi-rev sweep must not raise and must surface at least
    # one multi-rev cell among the cells it optimises.
    ledger = tmp_path / "ledger.jsonl"
    results = list(
        discover(
            bodies=("E", "M"),
            k_synodic=2,
            vinf_cap=8.0,
            ledger_path=ledger,
            l_max=4,
            n_max=1,
            branch_set=("single", "low"),
            max_cells=3,
        )
    )
    assert results, "expected at least one solved cell"
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/data/test_discover.py::test_discover_accepts_multirev_params -v`
Expected: FAIL — `discover()` does not yet accept `l_max`/`n_max`/`branch_set`/`max_cells` (TypeError).

- [ ] **Step 3: Add the params to `discover`**

Read the current `discover(...)` signature in `src/cyclerfinder/data/discover.py` and thread `l_max`, `n_max`, `branch_set`, `max_cells` into its `feasible_cells(...)` call (which already accepts `n_max`/`branch_set`). Bound the cell stream with `itertools.islice(..., max_cells)` when `max_cells` is set. Default `n_max=0`, `branch_set=("single",)`, `max_cells=None` to keep existing callers bitwise-identical.

- [ ] **Step 4: Run to verify it passes**

Run: `uv run pytest tests/data/test_discover.py::test_discover_accepts_multirev_params -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/cyclerfinder/data/discover.py tests/data/test_discover.py
git commit -m "data/discover: thread multi-rev enumeration (n_max, branch_set, max_cells) into the sweep"
```

---

### Task 7: SnLm family sweep gate (sourced-anchor matched)

**Files:**
- Test: `tests/data/test_discover.py`

- [ ] **Step 1: Write the failing test**

```python
@pytest.mark.slow
def test_snlm_sweep_rediscovers_a_sourced_anchor(tmp_path):
    """A bounded 2-synodic E-M multi-rev sweep must surface at least one
    closed, constraint-satisfying cycler whose Earth V-infinity matches a
    sourced SnLm anchor (5.65 km/s) within the gauntlet tolerance."""
    import numpy as np
    from cyclerfinder.data.discover import discover

    ledger = tmp_path / "ledger.jsonl"
    matched = False
    for opt_result, _match, _level in discover(
        bodies=("E", "M"),
        k_synodic=2,
        vinf_cap=8.0,
        ledger_path=ledger,
        l_max=4,
        n_max=1,
        branch_set=("single", "low"),
        max_cells=24,
    ):
        if not opt_result.constraints_satisfied:
            continue
        for enc in opt_result.best_cycler.encounters:
            if enc.body == "E":
                v = max(float(np.linalg.norm(enc.vinf_in)),
                        float(np.linalg.norm(enc.vinf_out)))
                if abs(v - 5.65) < 0.3:
                    matched = True
    assert matched, "no swept multi-rev cell matched the 5.65 km/s SnLm Earth anchor"
```

- [ ] **Step 2: Run to verify it fails or passes**

Run: `uv run pytest tests/data/test_discover.py::test_snlm_sweep_rediscovers_a_sourced_anchor -v`
Expected: PASS if the sweep reaches the S1L1 basin within `max_cells=24`; if FAIL, raise `max_cells`/`n_max` or widen `branch_set` until the known S1L1 topology (Task 1) is inside the swept set, then re-run. Do NOT loosen the 0.3 km/s tolerance — that is the sourced-anchor bound.

- [ ] **Step 3: Commit**

```bash
git add tests/data/test_discover.py
git commit -m "tests: SnLm multi-rev sweep rediscovers the 5.65 km/s sourced Earth anchor"
```

---

### Task 8: Full-suite regression + lint + type gate

- [ ] **Step 1: Run the full suite**

Run: `uv run pytest -q`
Expected: all green (plus the documented xfails). Record the pass/skip/xfail counts.

- [ ] **Step 2: Lint + type**

Run: `uv run ruff check . && uv run ruff format --check . && uv run mypy src`
Expected: all pass.

- [ ] **Step 3: Commit any fixups, then stop for review**

```bash
git add -A && git commit -m "snlm: full-suite regression green after multi-rev rediscovery + sweep"
```

---

## Self-Review

**Spec coverage:** Phase 1 (S1L1 re-model + computed rediscovery) = Tasks 1–5; Phase 2 (family sweep) = Tasks 6–8. The chosen design forks (Both-sequenced; Derive+V1-rediscovery) are honoured: S1/L1 ToFs are computed and labelled, only V∞ anchors are sourced.

**Placeholder scan:** The only intentionally execution-determined values are the `WINNING_*` tuples (Task 1 derives them; Tasks 2–3 pin them) — this is a genuine data dependency with a fully-specified derivation procedure, not a hand-waved placeholder. Every code step ships real code.

**Type consistency:** `Cell(bodies, sequence, period_k, per_leg_revs, per_leg_branch)`, `optimise_cell_idealized(..., vinf_cap=, n_starts=, seed=, use_de=)`, `OptimisationResult.best_cycler/.constraints_satisfied/.closure_residual_kms`, and `Score.max_vinf_kms` match the source signatures read at plan time.

**Risk:** Task 1 may yield zero hits (S1L1 not hostable in circular-coplanar, as with Aldrin). That branch is handled explicitly in Task 1 Step 2 and Task 4 Step 1 — the work degrades to an honest model-limitation xfail rather than a forced/fabricated pass.
