# DSM Multi-Arc Closure Lane Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Wire the existing Takao η-DSM genome (`dsm_chain_correct`) into a catalogue-row closure-and-validation lane that re-derives the ~8 descriptor-bearing off-family Russell rows on-family as two-arc trajectories and proposes their off-V0 promotion.

**Architecture:** A descriptor→DSM-seed adapter feeds the existing `dsm_chain_correct` (real-eph DE440); a per-row closer emerges V∞; a gate chain (single-arc-degenerate reject → V1 anchor-match + izzo/gooding crosscheck → V3 n-body) decides promotion; a batch driver runs the 8 rows and proposes `_LEVEL_EVIDENCE` entries with NO catalogue writeback (held for review). Mirrors the #170 App-C batch.

**Tech Stack:** Python 3.11 (uv), numpy, scipy, pytest. Reuses `cyclerfinder.search.dsm_leg` (`dsm_chain_correct`, `sequence_keyed_bounds`, `dsm_chain_decision_vector`), `cyclerfinder.search.self_seeding` (`g_arc_branches`), `cyclerfinder.search.free_return_chain` (`single_arc_degenerate`), `cyclerfinder.verify.crosscheck`, the REBOUND n-body rung, and `cyclerfinder.data.validate` (`_LEVEL_EVIDENCE`).

**Spec:** `docs/superpowers/specs/2026-06-10-dsm-multiarc-closure-lane-design.md`

---

## Conventions every task must follow (project rules — read first)

- **uv only.** First line of every shell block: `export PATH="$HOME/.local/bin:$PATH"`. Use `uv run ...`. Never bare `python`, never `pip`.
- **mypy is strict (also on tests).** Two recurring gotchas: (1) every test function needs an explicit `-> None`; (2) wrap `**`/`np.*` expressions returned from a `-> float` function in `float(...)`. **Test imports MUST use the direct-module form** `import cyclerfinder.search.dsm_leg as dl` (NOT `from cyclerfinder.search import dsm_leg as dl` — the latter trips mypy `attr-defined` and the pre-commit hook rejects the commit). `uv run ruff check --fix <file>` clears SIM300 Yoda conditions. The pre-commit hook runs ruff+mypy+jsonschema and REJECTS on failure — apply these up front.
- **Commit messages:** `subsystem: description`. **No `Co-Authored-By` / AI-attribution trailers.**
- **Never push. Never branch** — work on `main`. Explicit pathspecs in `git add`; run `git status --short` before each commit and confirm only your files are staged (ignore `.playwright-mcp/`).
- **Golden tests:** the EXPECTED side traces to a published source (the Russell-table V∞ cell / the row's stored anchor), NEVER a value our own code computed.
- **No tolerance loosening, no catalogue writeback** in this plan. A row that won't close on-family is a recorded NEGATIVE — success. Writeback is a separate post-review step (Task 7 documents it; it does not run here).

## Background the implementer needs

The 8 in-scope rows are the descriptor-bearing off-family `russell-ch4` rows from the #177 triage (`docs/notes/2026-06-08-self-seeding-triage-results.md`): the 6 REACHABLE (`9.353Gg2`, `9.94Gg3`, `5.30ggF3`, `3.78Gg3`, and two more in that note's table) + the 2 OFF-FAMILY-NO-CLOSE (`5.30gGf3`, `5.75ggF3`). The 204 `russell-ocampo` rows have NO per-arc descriptor and are OUT OF SCOPE (the adapter returns `None` for them).

Key reused APIs (verified):
- `self_seeding.g_arc_branches(aphelion_au, g_tof_years, big_g_tof_years, vinf_e_anchor, vinf_m_anchor, *, max_g_revs=1, mu=MU_SUN_KM3_S2) -> list[GArcShape]` — descriptor → all Mars-crossing arc branch shapes (`[0]` is the base short-way). A `GArcShape` carries `a`, `e`, `tof_g_days`, `branch`, `g_revs` (read the dataclass at `self_seeding.py:92`).
- `dsm_leg.sequence_keyed_bounds(*, sequence, t0_window_sec, vinf_out0_bounds_kms=(1.0,5.1), charge_flyby_continuity=False) -> DsmBounds`.
- `dsm_leg.dsm_chain_decision_vector(*, t0_sec, vinf_out0_kms, alpha0, beta0, tof_days_per_leg, eta_per_leg, alpha_int_per_leg=(), beta_int_per_leg=()) -> np.ndarray` — layout `[t0, vinf_out0, alpha0, beta0, *tof, *eta, *alpha_int, *beta_int]`.
- `dsm_leg.dsm_chain_correct(x0, *, sequence, ephem, bounds=None, mu=MU_SUN_KM3_S2, rendezvous=False, tol_kms=0.1, max_nfev=200, max_revs=0, rev_branch_per_leg=None, charge_flyby_continuity=False) -> DsmChainResult`. Use `charge_flyby_continuity=True` (the #162 vector-residual mode — the only one that rewards the bend-feasible low-V∞ basin). `DsmChainResult` carries `converged`, `max_residual_kms`, `dv_dsm_per_leg_kms`, `vinf_in_kms`, `vinf_out_kms` (per-body EMERGED V∞ — the evidence), `eta_per_leg`, `tof_days_per_leg`, `t0_sec`.
- `free_return_chain.single_arc_degenerate(aphelion_au, arc_tof_years, vinf_e_anchor, vinf_m_anchor, *, mu=MU_SUN_KM3_S2, tol_kms=0.5) -> FreeReturnChainResult` — the degenerate single-ellipse check.

---

## File Structure

- **Create** `src/cyclerfinder/search/dsm_descriptor_seed.py` — `DsmChainSeed`, `seed_dsm_chain_from_descriptor`, `DsmClosureResult`, `close_row_dsm`. One responsibility: turn a catalogue row into a DSM closure attempt + result. (Closer lives with the seed because they share the row→seed→correct flow and are small.)
- **Create** `src/cyclerfinder/search/dsm_closure_gate.py` — `dsm_closure_verdict(row, result, ephem) -> DsmVerdict` (the single-arc-degenerate guard + V1 gate; V3 is invoked by the batch, slow).
- **Create** `scripts/dsm_closure_batch.py` — the 8-row batch driver (runlog + results note + proposed evidence; no writeback).
- **Create** tests: `tests/search/test_dsm_descriptor_seed.py`, `tests/search/test_dsm_closure_gate.py`.
- **Modify** none of the reused modules (pure additions).

---

## Task 1: Descriptor → DSM seed adapter

**Files:**
- Create: `src/cyclerfinder/search/dsm_descriptor_seed.py`
- Test: `tests/search/test_dsm_descriptor_seed.py`

First READ a sample row's fields: `cd /home/bruce/dev/cyclers && uv run python -c "from cyclerfinder.data.catalog import load_catalogue; r=[x for x in load_catalogue() if x['id']=='russell-ch4-9.353Gg2'][0]; import json; print(json.dumps(r, indent=2, default=str))"` to see the exact descriptor / orbit_elements / V∞-anchor / bodies fields available. Also read `scripts/triage_self_seeding.py` to see how #177 extracted `aphelion_au`, the g/G ToFs, and the V∞ E/M anchors from a row — mirror that extraction here (do not invent field names).

- [ ] **Step 1: Write the failing test**

```python
# tests/search/test_dsm_descriptor_seed.py
"""Descriptor -> DSM seed adapter (plan 2026-06-10, Component 1)."""
from __future__ import annotations

import numpy as np

import cyclerfinder.search.dsm_descriptor_seed as dds
from cyclerfinder.data.catalog import load_catalogue


def _row(row_id: str) -> dict:
    return next(r for r in load_catalogue() if r["id"] == row_id)


def test_seed_built_for_reachable_descriptor_row() -> None:
    seed = dds.seed_dsm_chain_from_descriptor(_row("russell-ch4-9.353Gg2"))
    assert seed is not None
    # Sequence is the row's encounter chain; decision vector matches the layout
    # [t0, vinf_out0, alpha0, beta0, *tof, *eta] with eta seeded ballistic (0).
    assert len(seed.sequence) >= 2
    assert seed.x0.shape[0] == 4 + 2 * (len(seed.sequence) - 1)
    n_legs = len(seed.sequence) - 1
    eta = seed.x0[4 + n_legs : 4 + 2 * n_legs]
    assert np.allclose(eta, 0.0)
    assert seed.vinf_anchor_kms > 0.0


def test_no_descriptor_row_returns_none() -> None:
    # An ocampo row has the n.m.k summary format, no per-arc g/G descriptor.
    ocampo = next(r for r in load_catalogue() if r["id"].startswith("russell-ocampo"))
    assert dds.seed_dsm_chain_from_descriptor(ocampo) is None
```

- [ ] **Step 2: Run test to verify it fails**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run pytest tests/search/test_dsm_descriptor_seed.py::test_seed_built_for_reachable_descriptor_row -v
```
Expected: FAIL — `ModuleNotFoundError` / `AttributeError: ... 'seed_dsm_chain_from_descriptor'`.

- [ ] **Step 3: Write minimal implementation**

```python
# src/cyclerfinder/search/dsm_descriptor_seed.py
"""Descriptor -> Takao-DSM closure seed + per-row closer (spec 2026-06-10).

Bridges the two multi-arc mechanisms: parse a catalogue row's 2-arc g/G
free-return descriptor, use :func:`cyclerfinder.search.self_seeding.g_arc_branches`
to get the coplanar arc shape + the transit branch matching the tabulated transit,
and assemble the decision vector + bounds that the Takao η-DSM corrector
(:func:`cyclerfinder.search.dsm_leg.dsm_chain_correct`) consumes. Pure.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import MU_SUN_KM3_S2
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search import dsm_leg, self_seeding

DAY_S = 86400.0


@dataclass(frozen=True)
class DsmChainSeed:
    sequence: tuple[str, ...]
    x0: NDArray[np.float64]            # [t0, vinf_out0, alpha0, beta0, *tof, *eta]
    bounds: dsm_leg.DsmBounds
    arc_a_au: float                    # coplanar descriptor arc shape
    arc_e: float
    transit_branch: str
    vinf_anchor_kms: float             # the row's sourced Russell-table V∞ cell


def _descriptor_params(row: dict) -> tuple[float, float, float, float, float, tuple[str, ...]] | None:
    """Extract (aphelion_au, g_tof_yr, big_g_tof_yr, vinf_e, vinf_m, sequence) from
    the row's g/G descriptor, or None if the row has no per-arc descriptor.

    MIRROR scripts/triage_self_seeding.py's extraction — use the SAME field paths it
    reads (descriptor string / segments / orbit_elements / vinf anchors). Return
    None when any required descriptor field is absent (the ocampo rows).
    """
    ...  # implement by mirroring triage_self_seeding.py (read it first)


def seed_dsm_chain_from_descriptor(row: dict) -> DsmChainSeed | None:
    params = _descriptor_params(row)
    if params is None:
        return None
    aphelion_au, g_tof_yr, big_g_tof_yr, vinf_e, vinf_m, sequence = params
    branches = self_seeding.g_arc_branches(
        aphelion_au, g_tof_yr, big_g_tof_yr, vinf_e, vinf_m
    )
    if not branches:
        return None
    arc = branches[0]  # base short-way shape; the gate may retry others
    n_legs = len(sequence) - 1
    # Seed ToFs: split the descriptor's total transit across legs (the arc branch
    # ToF on the Mars-transit leg, synodic-derived on the others). Use the branch's
    # tof_g_days for the transit leg, the row's beat for the slack leg.
    tof_seed_days = tuple(float(arc.tof_g_days) for _ in range(n_legs))
    eta_seed = tuple(0.0 for _ in range(n_legs))  # start ballistic
    t0_seed_sec = 0.0
    bounds = dsm_leg.sequence_keyed_bounds(
        sequence=sequence,
        t0_window_sec=(-arc.tof_g_days * DAY_S, arc.tof_g_days * DAY_S),
        vinf_out0_bounds_kms=(max(0.5, vinf_e - 2.0), vinf_e + 2.0),
        charge_flyby_continuity=True,
    )
    x0 = dsm_leg.dsm_chain_decision_vector(
        t0_sec=t0_seed_sec,
        vinf_out0_kms=vinf_e,
        alpha0=0.0,
        beta0=0.0,
        tof_days_per_leg=tof_seed_days,
        eta_per_leg=eta_seed,
        alpha_int_per_leg=tuple(0.0 for _ in range(n_legs - 1)),
        beta_int_per_leg=tuple(0.0 for _ in range(n_legs - 1)),
    )
    return DsmChainSeed(
        sequence=sequence,
        x0=x0,
        bounds=bounds,
        arc_a_au=float(arc.a),
        arc_e=float(arc.e),
        transit_branch=str(arc.branch),
        vinf_anchor_kms=float(vinf_m),
    )
```

Implement `_descriptor_params` by reading `scripts/triage_self_seeding.py` and reusing its exact row-field extraction. Do NOT guess field names — if the extraction path is unclear after reading that script and a sample row, report NEEDS_CONTEXT.

- [ ] **Step 4: Run test to verify it passes**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run pytest tests/search/test_dsm_descriptor_seed.py -q
```
Expected: PASS (2 passed). If the REACHABLE row's `g_arc_branches` returns empty or `_descriptor_params` can't find the fields, report BLOCKED/NEEDS_CONTEXT with the row dump — do not fabricate a descriptor.

- [ ] **Step 5: Commit**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run ruff check src/cyclerfinder/search/dsm_descriptor_seed.py tests/search/test_dsm_descriptor_seed.py
uv run ruff format src/cyclerfinder/search/dsm_descriptor_seed.py tests/search/test_dsm_descriptor_seed.py
git add src/cyclerfinder/search/dsm_descriptor_seed.py tests/search/test_dsm_descriptor_seed.py
git status --short
git commit -m "search/dsm_descriptor_seed: descriptor -> Takao-DSM closure seed adapter"
```

---

## Task 2: Per-row DSM closer

**Files:**
- Modify: `src/cyclerfinder/search/dsm_descriptor_seed.py`
- Test: `tests/search/test_dsm_descriptor_seed.py`

- [ ] **Step 1: Write the failing test**

```python
# append to tests/search/test_dsm_descriptor_seed.py
from cyclerfinder.core.ephemeris import Ephemeris


def test_close_reachable_row_emerges_vinf_near_anchor() -> None:
    # A REACHABLE descriptor row closes with the DSM genome on the real ephemeris,
    # and its EMERGED Mars V∞ lands within tolerance of the row's sourced anchor.
    row = _row("russell-ch4-9.353Gg2")
    ephem = Ephemeris(model="de440")
    res = dds.close_row_dsm(row, ephem)
    # converged is by the corrector's own residual criterion; if it converges the
    # emerged V∞ must match the sourced anchor (golden — anchor from the row, not
    # computed here). A non-converged row is a recorded negative (also valid).
    if res.converged:
        assert res.anchor_match
        assert min(res.dv_dsm_kms) >= 0.0
```

- [ ] **Step 2: Run test to verify it fails**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run pytest tests/search/test_dsm_descriptor_seed.py::test_close_reachable_row_emerges_vinf_near_anchor -v
```
Expected: FAIL — `AttributeError: ... 'close_row_dsm'`.

- [ ] **Step 3: Write minimal implementation**

```python
# add to imports
from cyclerfinder.data import validate as _validate  # for V1_TOLERANCE if exported
```
```python
# append to src/cyclerfinder/search/dsm_descriptor_seed.py

V1_TOLERANCE_KMS = 0.5  # spec §14 V1 anchor-match band (km/s); confirm vs the
# value used in data/validate.py / verify and align to it (do not invent a looser one).


@dataclass(frozen=True)
class DsmClosureResult:
    converged: bool
    max_residual_kms: float            # CONSTRAINED
    dv_dsm_kms: tuple[float, ...]      # CONSTRAINED — per-leg interior impulse
    vinf_per_encounter_kms: tuple[float, ...]  # EMERGED
    vinf_anchor_kms: float
    anchor_match: bool
    hyperbolic_impossible: bool
    seed: DsmChainSeed | None


def close_row_dsm(row: dict, ephem: Ephemeris, *, tol_kms: float = 0.1) -> DsmClosureResult:
    seed = seed_dsm_chain_from_descriptor(row)
    if seed is None:
        return DsmClosureResult(
            converged=False, max_residual_kms=float("nan"), dv_dsm_kms=(),
            vinf_per_encounter_kms=(), vinf_anchor_kms=float("nan"),
            anchor_match=False, hyperbolic_impossible=False, seed=None,
        )
    res = dsm_leg.dsm_chain_correct(
        seed.x0, sequence=seed.sequence, ephem=ephem, bounds=seed.bounds,
        tol_kms=tol_kms, charge_flyby_continuity=True,
    )
    vinf_out = tuple(float(v) for v in res.vinf_out_kms)
    # Match the Mars-encounter emerged V∞ to the sourced anchor.
    anchor = seed.vinf_anchor_kms
    best = min((abs(v - anchor) for v in vinf_out), default=float("inf"))
    hyper = any(v > 71.9 for v in vinf_out)  # heliocentric elliptic ceiling flag
    return DsmClosureResult(
        converged=bool(res.converged),
        max_residual_kms=float(res.max_residual_kms),
        dv_dsm_kms=tuple(float(d) for d in res.dv_dsm_per_leg_kms),
        vinf_per_encounter_kms=vinf_out,
        vinf_anchor_kms=anchor,
        anchor_match=best <= V1_TOLERANCE_KMS,
        hyperbolic_impossible=hyper,
        seed=seed,
    )
```

- [ ] **Step 4: Run test to verify it passes**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run pytest tests/search/test_dsm_descriptor_seed.py -q
```
Expected: PASS (3 passed). The closer may or may not converge for this row — both branches of the assert pass. Note the actual outcome in your report.

- [ ] **Step 5: Commit**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run ruff check src/cyclerfinder/search/dsm_descriptor_seed.py tests/search/test_dsm_descriptor_seed.py
uv run ruff format src/cyclerfinder/search/dsm_descriptor_seed.py tests/search/test_dsm_descriptor_seed.py
git add src/cyclerfinder/search/dsm_descriptor_seed.py tests/search/test_dsm_descriptor_seed.py
git status --short
git commit -m "search/dsm_descriptor_seed: per-row DSM closer (real-eph, emerged V∞ vs anchor)"
```

---

## Task 3: Closure gate — single-arc-degenerate guard + V1

**Files:**
- Create: `src/cyclerfinder/search/dsm_closure_gate.py`
- Test: `tests/search/test_dsm_closure_gate.py`

First READ `src/cyclerfinder/verify/crosscheck.py` for the V1 mechanism (lamberthub izzo2015 + gooding1990 per-leg agreement + Kepler reprop) and how #170 / `appc_corrected` called it, so the V1 gate reuses it rather than re-implementing.

- [ ] **Step 1: Write the failing test**

```python
# tests/search/test_dsm_closure_gate.py
"""DSM closure gate: degenerate guard + V1 (plan 2026-06-10, Component 3)."""
from __future__ import annotations

import cyclerfinder.search.dsm_closure_gate as gate
import cyclerfinder.search.free_return_chain as frc


def test_single_arc_degenerate_row_is_rejected() -> None:
    # If a plain single-ellipse already closes the row to within the degenerate
    # tolerance, the DSM lane must REJECT it (not a genuine multi-arc row).
    # Use a descriptor whose two arcs coincide (single_arc_degenerate converges).
    assert gate.is_single_arc_degenerate(
        aphelion_au=2.5, arc_tof_years=0.7, vinf_e_anchor=4.0, vinf_m_anchor=3.0
    ) in (True, False)  # smoke: callable + returns bool; real assert below


def test_degenerate_helper_matches_free_return_chain() -> None:
    # The guard delegates to free_return_chain.single_arc_degenerate and reads its
    # converged flag — no independent re-derivation.
    res = frc.single_arc_degenerate(2.5, 0.7, 4.0, 3.0)
    assert gate.is_single_arc_degenerate(
        aphelion_au=2.5, arc_tof_years=0.7, vinf_e_anchor=4.0, vinf_m_anchor=3.0
    ) == bool(res.converged)
```

- [ ] **Step 2: Run test to verify it fails**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run pytest tests/search/test_dsm_closure_gate.py -q
```
Expected: FAIL — `ModuleNotFoundError` / `AttributeError: ... 'is_single_arc_degenerate'`.

- [ ] **Step 3: Write minimal implementation**

```python
# src/cyclerfinder/search/dsm_closure_gate.py
"""DSM closure gate (spec 2026-06-10, Component 3).

Decides whether a DSM closure earns a validation-level proposal: (1) reject if a
single ellipse already closes the row (not genuinely multi-arc); (2) V1 = emerged
V∞ matches the sourced anchor AND the lamberthub two-method + Kepler-reprop
crosscheck passes. V3 (n-body) is invoked by the batch driver (slow). Pure except
for the crosscheck call. No writeback.
"""
from __future__ import annotations

from dataclasses import dataclass

from cyclerfinder.core.constants import MU_SUN_KM3_S2
from cyclerfinder.search import free_return_chain


def is_single_arc_degenerate(
    *, aphelion_au: float, arc_tof_years: float, vinf_e_anchor: float,
    vinf_m_anchor: float, tol_kms: float = 0.5,
) -> bool:
    """True iff a single free-return ellipse already closes this descriptor — in
    which case the DSM result must be rejected (no false multi-arc claim)."""
    res = free_return_chain.single_arc_degenerate(
        aphelion_au, arc_tof_years, vinf_e_anchor, vinf_m_anchor,
        mu=MU_SUN_KM3_S2, tol_kms=tol_kms,
    )
    return bool(res.converged)


@dataclass(frozen=True)
class DsmVerdict:
    level: str                  # "V0" | "V1" (V3 added by the batch on n-body pass)
    accepted: bool              # passed the degenerate guard + V1
    reason: str
```

Then add a `v1_verdict(row, result)` that returns a `DsmVerdict`: `accepted=False level="V0"` if the degenerate guard fires or `not result.converged` or `not result.anchor_match`; else run the `verify/crosscheck` V1 predicate on the closed legs and return `level="V1"` if it passes. Add the test asserting a converged+anchor-matched+non-degenerate row reaches `level="V1"` and a degenerate one stays `"V0"` (construct the inputs from a real row via Task 2's `close_row_dsm`).

- [ ] **Step 4: Run test to verify it passes**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run pytest tests/search/test_dsm_closure_gate.py -q
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run ruff check src/cyclerfinder/search/dsm_closure_gate.py tests/search/test_dsm_closure_gate.py
uv run ruff format src/cyclerfinder/search/dsm_closure_gate.py tests/search/test_dsm_closure_gate.py
git add src/cyclerfinder/search/dsm_closure_gate.py tests/search/test_dsm_closure_gate.py
git status --short
git commit -m "search/dsm_closure_gate: single-arc-degenerate guard + V1 anchor/crosscheck gate"
```

---

## Task 4: Batch driver over the 8 rows

**Files:**
- Create: `scripts/dsm_closure_batch.py`

READ `scripts/triage_self_seeding.py` (the #177 batch shape) and a `scripts/`-level batch like the App-C one for the runlog/report idioms.

- [ ] **Step 1: Write the driver**

```python
# scripts/dsm_closure_batch.py
"""DSM multi-arc closure batch over the 8 descriptor-bearing off-family rows.

Per row: seed -> close_row_dsm (real-eph) -> degenerate guard -> V1 gate; writes a
JSONL runlog and a results note with the PROPOSED _LEVEL_EVIDENCE text for passes.
NO catalogue writeback (held for session review). V3/n-body is a separate slow
step (Task 5). Run:

    export PATH="$HOME/.local/bin:$PATH"
    uv run python scripts/dsm_closure_batch.py --report /tmp/dsm_closure.txt
"""
from __future__ import annotations

import argparse
import json
import sys

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.data.catalog import load_catalogue
from cyclerfinder.search import dsm_closure_gate, dsm_descriptor_seed

ROWS = (
    "russell-ch4-9.353Gg2", "russell-ch4-9.94Gg3", "russell-ch4-5.30ggF3",
    "russell-ch4-3.78Gg3",  # + the remaining 2 REACHABLE ids from the #177 note
    "russell-ch4-5.30gGf3", "russell-ch4-5.75ggF3",  # the 2 NO-CLOSE
)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--report", default=None)
    p.add_argument("--runlog", default="data/runs/dsm-closure.jsonl")
    args = p.parse_args(argv)

    catalogue = {r["id"]: r for r in load_catalogue()}
    ephem = Ephemeris(model="de440")
    lines: list[str] = []
    with open(args.runlog, "w", encoding="utf-8") as fh:
        for rid in ROWS:
            row = catalogue.get(rid)
            if row is None:
                lines.append(f"{rid}: MISSING from catalogue")
                continue
            res = dsm_descriptor_seed.close_row_dsm(row, ephem)
            verdict = dsm_closure_gate.v1_verdict(row, res)
            rec = {
                "id": rid, "converged": res.converged,
                "max_residual_kms": res.max_residual_kms,
                "emerged_vinf": list(res.vinf_per_encounter_kms),
                "anchor_kms": res.vinf_anchor_kms, "anchor_match": res.anchor_match,
                "level": verdict.level, "accepted": verdict.accepted,
                "reason": verdict.reason,
            }
            fh.write(json.dumps(rec) + "\n")
            lines.append(
                f"{rid}: conv={res.converged} anchor_match={res.anchor_match} "
                f"-> {verdict.level} ({verdict.reason})"
            )

    report = "\n".join(lines)
    print(report)
    if args.report:
        from pathlib import Path
        Path(args.report).write_text(report + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```
Fill the two missing REACHABLE row ids from the #177 triage note's table before running.

- [ ] **Step 2: Run the batch**

```bash
export PATH="$HOME/.local/bin:$PATH"
date -Iseconds
uv run python scripts/dsm_closure_batch.py --report /tmp/dsm_closure.txt
date -Iseconds
```
Expected: a per-row line for all 8. Some converge to V1, some stay V0 (recorded negatives). Record the output VERBATIM in your report. **Do not** tune anything to force a pass.

- [ ] **Step 3: Write the results note**

Create `docs/notes/2026-06-10-dsm-closure-batch-results.md` with the verbatim table (row, converged, anchor_match, level, reason) and, for each V1 pass, the PROPOSED `_LEVEL_EVIDENCE` text (sourced-anchor + crosscheck wording, mirroring existing entries in `data/validate.py`). State explicitly: NO writeback performed; the 204 ocampo rows are out of scope (no descriptor).

- [ ] **Step 4: Commit**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run ruff check scripts/dsm_closure_batch.py
uv run ruff format scripts/dsm_closure_batch.py
git add scripts/dsm_closure_batch.py data/runs/dsm-closure.jsonl docs/notes/2026-06-10-dsm-closure-batch-results.md
git status --short
git commit -m "scripts: DSM multi-arc closure batch over the 8 descriptor-bearing rows (no writeback)"
```

---

## Task 5: V3 n-body confirmation for any V1 pass (slow, gated)

**Files:**
- Modify: `src/cyclerfinder/search/dsm_closure_gate.py`
- Test: `tests/search/test_dsm_closure_gate.py`

Only matters if Task 4 produced a V1 pass. Implement the hook + a fast contract test regardless; the actual propagation is `@pytest.mark.slow` and run manually.

- [ ] **Step 1: Write the failing contract test**

```python
# append to tests/search/test_dsm_closure_gate.py
import cyclerfinder.search.dsm_descriptor_seed as dds
from cyclerfinder.core.ephemeris import Ephemeris


def test_dsm_route_to_nbody_request_shape() -> None:
    row = next(
        r for r in __import__("cyclerfinder.data.catalog", fromlist=["load_catalogue"]).load_catalogue()
        if r["id"] == "russell-ch4-9.353Gg2"
    )
    res = dds.close_row_dsm(row, Ephemeris(model="de440"))
    req = gate.dsm_closure_to_nbody_request(row, res)
    assert req["id"] == "russell-ch4-9.353Gg2"
    assert "sequence" in req and "t0_sec" in req
```

- [ ] **Step 2: Run to verify it fails**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run pytest tests/search/test_dsm_closure_gate.py::test_dsm_route_to_nbody_request_shape -q
```
Expected: FAIL — `AttributeError: ... 'dsm_closure_to_nbody_request'`.

- [ ] **Step 3: Write minimal implementation**

Add `dsm_closure_to_nbody_request(row, result) -> dict` to `dsm_closure_gate.py` returning the dict the existing heliocentric REBOUND rung consumes (READ the #134/#170 n-body call site — `cyclerfinder.nbody.rung` or the harness used by `appc_corrected`/the App-C batch — to match keys: sequence, t0_sec, per-leg ToFs, per-leg ΔV, emerged V∞). Do NOT run propagation in this function.

- [ ] **Step 4: Run to verify it passes**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run pytest tests/search/test_dsm_closure_gate.py -q
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run ruff check src/cyclerfinder/search/dsm_closure_gate.py tests/search/test_dsm_closure_gate.py
uv run ruff format src/cyclerfinder/search/dsm_closure_gate.py tests/search/test_dsm_closure_gate.py
git add src/cyclerfinder/search/dsm_closure_gate.py tests/search/test_dsm_closure_gate.py
git status --short
git commit -m "search/dsm_closure_gate: DSM closure -> n-body (V3) request hook"
```

---

## Task 6: Full regression + lint + type gate

**Files:** none (verification only).

- [ ] **Step 1: Full suite (excluding slow)**

```bash
export PATH="$HOME/.local/bin:$PATH"
date -Iseconds
uv run pytest -q -m "not slow"
date -Iseconds
```
Expected: all PASS (no regressions).

- [ ] **Step 2: Lint + format + types**

```bash
export PATH="$HOME/.local/bin:$PATH"
uv run ruff check .
uv run ruff format --check .
uv run mypy src/cyclerfinder/search/dsm_descriptor_seed.py src/cyclerfinder/search/dsm_closure_gate.py
```
Expected: all clean. Fix anything introduced.

- [ ] **Step 3: Final commit (only if Step 2 needed fixes)**

```bash
export PATH="$HOME/.local/bin:$PATH"
git add -p
git status --short
git commit -m "search/dsm: lint + type fixes for the DSM closure lane"
```

---

## Task 7: Writeback (DOCUMENTED ONLY — runs after session review, not in this plan)

This task is **not executed by the implementer.** It records the post-review procedure so the controller can do it once the batch results are reviewed:

1. For each row the batch marked `level="V1"` (and any later n-body-confirmed `V3`), add a `_LEVEL_EVIDENCE[(row_id, level)]` entry in `src/cyclerfinder/data/validate.py` with the sourced-anchor + crosscheck wording (mirror existing entries).
2. Bump `validation_level` for that row in `data/catalogue.yaml`.
3. Run `uv run pytest -q` (the validate-invariants + over-claim-guard tests must pass) and commit as `data/catalogue: promote <row> V0->V1 (DSM multi-arc, sourced-anchor + crosscheck)`.
4. The site auto-syncs from GitHub `main` at deploy; a redeploy (`gh workflow run deploy.yml` in cyclers.space) surfaces the new tiers.

**No writeback happens during the build.** The build's deliverable is the batch results note + proposed evidence, held for review.

---

## Post-implementation (controller, not a task)

- Summarise: how many of the 8 rows closed on-family with the DSM genome (V1), how many stayed V0 (recorded negatives), and the per-row emerged-V∞-vs-anchor numbers.
- Honest framing: even a full pass is ≤8 rows; the 204 ocampo rows remain descriptor-gated (publication gap, #116-adjacent). This lifts the descriptor-bearing tail, not the bulk.
- Deferred: real-eph inclination/broken-plane on these rows; the writeback (Task 7) after review.
