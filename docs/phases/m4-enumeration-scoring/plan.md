# M4 — Enumeration + scoring

**Spec reference:** spec.md §4 (architecture — `search/sequence.py`, `model/score.py`), §5 step 2 (Tisserand enumeration) and step 4 (score & filter), §6 (interfaces sketch), §8 (M4 milestone definition — "enumerator rejects energetically infeasible sequences"), §10 (search-landscape risks), §12(d) (hard constraints, not soft regularisers), §13.1–§13.8 (the **central design document for M4** — the structural-cell decomposition, iterative deepening, Tisserand pruning, deduplication, work queue, prioritisation, execution recipe), §16 (catalogue schema — referenced by the `taxi_cost` metric and by the canonicalisation we deliberately **defer to M7**).

**Purpose:** stand up the **discrete combinatorial layer** of the cycler finder. M0–M3 produced the primitives (constants, ephemeris, Lambert, Kepler, flyby, Tisserand, resonance, cycler/leg/encounter dataclasses, patched-conic constructor) and proved them against the Aldrin published numbers. M4 now produces (a) the **cell enumerator** — generators that walk the `(body_set, sequence, k, revs, branches)` discrete structure under iterative-deepening caps, with Tisserand pruning as the gate — and (b) the **scoring layer** — a `Score` dataclass and `rank()` routine that turns built `Cycler` instances into a ranked, hard-constraint-filtered shortlist. **No optimisation, no ledger, no parallelism** — those are M5 and M7.

**Gate (definition of done):**
1. `tests/search/test_sequence.py::test_tisserand_pruning_rejects_low_vinf_em` asserts an `("E","M")` cell at `vinf_cap = 2.0 km/s` is rejected by `tisserand_feasible` (the bodies are not linkable that low).
2. `tests/search/test_sequence.py::test_tisserand_pruning_accepts_aldrin_vinf_em` asserts the same cell at `vinf_cap = 8.0 km/s` passes (the Aldrin family lives in this neighbourhood).
3. `tests/search/test_sequence.py::test_enumeration_count_em_l4_k2` asserts that `enumerate_cells(body_set=("E","M"), L_max=4, k_max=2, N_max=0, branch_set=("single",))` produces exactly the documented count (see §4.4 below for the derivation).
4. `tests/model/test_score.py::test_score_aldrin_passes_hard_constraints` asserts the Aldrin seed (rebuilt via M3's `build_aldrin_seed`) scores as `hard_constraints_pass=True` under sensible defaults, and that `composite_score` is finite and reproducible.
5. `uv run pytest`, `uv run ruff check .`, `uv run ruff format --check .`, `uv run mypy src tests` all clean. CI green on the M4 commit.

---

## 1. What this milestone delivers

Two new source modules and their test files.

### 1.1 `src/cyclerfinder/search/sequence.py` — the cell enumerator

The `Cell` frozen dataclass (the atomic unit of search per spec §13.1) plus three generator-style functions:

- `enumerate_cells(body_set, L_max, k_max, N_max, branch_set=("single",))` — yields *all* combinatorially possible `Cell` instances under the caps. **No pruning happens here.**
- `tisserand_feasible(cell, vinf_cap, ephem)` — the spec §13.3 pruning predicate. True iff every consecutive body pair in the cell's sequence is linkable at some V∞ ≤ `vinf_cap`. **This is the M4 gate.**
- `feasible_cells(body_set, L_max, k_max, N_max, vinf_cap, ephem)` — convenience: `(c for c in enumerate_cells(...) if tisserand_feasible(c, vinf_cap, ephem))`.
- `deepening_frontier(body_set, ephem, vinf_cap)` — optional generator over the iterative-deepening frontier (raises L_max, k_max, N_max progressively per spec §13.2). M4 implements the generator only. The work-queue + ledger (§13.6, §13.8) and parallel execution are explicitly **deferred to M7**.

### 1.2 `src/cyclerfinder/model/score.py` — ranking metrics

The `Score` frozen dataclass plus four routines:

- `score(cycler, ephem, vinf_cap, rp_factors=None)` — single-cycler reduction to a `Score`.
- `taxi_cost_kms(cycler)` — surrogate hyperbolic-rendezvous cost.
- `composite_score(s, weights)` — weighted-sum sortable scalar.
- `rank(cyclers, ephem, vinf_cap, n_keep=20, weights=None)` — apply hard-constraint filter then composite-score sort, return top N as `list[tuple[Score, Cycler]]`.

### 1.3 Two test files

- `tests/search/test_sequence.py` — enumeration counts, Tisserand pruning gate (the M4 binding gate), cell ID stability, frontier monotonicity.
- `tests/model/test_score.py` — Score dataclass invariants, hard-constraints semantics, taxi-cost surrogate sanity, ranking ordering.

### 1.4 Explicit non-goals (M4 boundaries)

These appear in adjacent milestones — **do not stub or partially implement them in M4**:

| Out of M4 | Where it lands | Why deferred |
|---|---|---|
| Searching encounter times within a cell to minimise residual / ΔV | **M5** `search/optimize.py` | M4 produces *structural* candidates only; the inner timing search is the next milestone's contract. |
| The work-queue, append-only ledger, atomic claim/skip protocol per §13.6 + §13.8 | **M7** | M7 owns the catalogue + finder loop; building the ledger here would couple M4 to a persistence layer it does not need to gate. |
| Parallel execution of cells across workers | **M7** (with the ledger) | Same reason — parallelism is a property of the runner, not the enumerator. |
| Canonical signature computation (spec §16.2 lexicographically-minimal rotation) | **M7** | M4 produces sequences in *some* rotation. M7's catalogue is the right place to canonicalise, where rotation-invariance is load-bearing. Coupling M4 to it would mean producing the same cell under multiple rotations *here*, slowing the enumeration. |
| Real-ephemeris (astropy) feasibility | **M6** | The `linkable` predicate is coplanar-only per M2 plan §3.2; that's correct for M4's enumeration of structural candidates. |
| CLI, viz, output JSON formatting | **M8** | M4's output is a Python iterator + a ranked list, not a serialised artefact. |

---

## 2. File tree after M4

```
cyclers/
├── … (M0/M1/M2/M3 layout preserved unchanged)
├── src/cyclerfinder/
│   ├── core/                          # unchanged
│   ├── search/
│   │   ├── __init__.py                # unchanged
│   │   ├── tisserand.py               # M2 — unchanged; sequence.py imports linkable()
│   │   ├── resonance.py               # M2 — unchanged
│   │   ├── construct.py               # M3 — unchanged
│   │   └── sequence.py                # NEW (M4)
│   └── model/
│       ├── __init__.py                # NEW re-export adds Score
│       ├── cycler.py                  # M3 — unchanged
│       └── score.py                   # NEW (M4)
└── tests/
    ├── … (M0/M1/M2/M3 tests preserved)
    ├── search/
    │   ├── __init__.py                # NEW (first test under tests/search/)
    │   └── test_sequence.py           # NEW (M4 — includes the gate tests)
    └── model/
        └── test_score.py              # NEW (M4)
```

Subpackages `verify/`, `data/`, `viz/` remain uncreated — they belong to M6/M7/M8. `search/optimize.py` stays absent until M5. **No changes to any M0–M3 source file.** M4 is purely additive at the module level; the only edit to a pre-existing file is `model/__init__.py` to re-export `Score` alongside `Cycler`, `Leg`, `Encounter` already there from M3.

---

## 3. Module designs

### 3.1 `search/sequence.py`

#### 3.1.1 The `Cell` dataclass (spec §13.1, §13.8)

```python
from __future__ import annotations
from dataclasses import dataclass

@dataclass(frozen=True)
class Cell:
    """Atomic unit of search per spec §13.1.

    Discrete structural specification of a cycler candidate. Continuous
    DOF (encounter epochs, phases) are bounded sub-problems *inside* a
    cell, optimised in M5 — not held on the dataclass.
    """

    bodies: tuple[str, ...]            # canonical body set, e.g. ("V","E","M")
    sequence: tuple[str, ...]          # flyby sequence, e.g. ("E","V","M","E","M")
    period_k: int                      # period in synodic multiples (k ≥ 1)
    per_leg_revs: tuple[int, ...]      # heliocentric revolutions per leg
    per_leg_branch: tuple[str, ...]    # "single" | "low" | "high" per leg

    @property
    def id(self) -> str:
        """Deterministic sortable identifier per spec §13.8.

        Format: ``{bodyset}|{sequence}|k{K}|r{revs}|b{branches}``

        Examples
        --------
        ``VEM|E-V-M-E-M-E|k3|r00101|blllll``  (spec §13.8 worked example)
        ``EM|E-M-E|k2|r00|bss``               (a 2-synodic E-M-E cell, M4 native)
        """
        bodyset = "".join(self.bodies)
        sequence = "-".join(self.sequence)
        revs = "".join(str(r) for r in self.per_leg_revs)
        # Branch letter: 'l'/'h' for low/high multi-rev; 's' for single (0-rev).
        branch_map = {"single": "s", "low": "l", "high": "h"}
        branches = "".join(branch_map[b] for b in self.per_leg_branch)
        return f"{bodyset}|{sequence}|k{self.period_k}|r{revs}|b{branches}"
```

**Design notes (binding):**

- **`Cell.id` is a `@property`, never a stored field.** Always derivable from the other fields, so the dataclass cannot drift out of sync with its identity. Two cells with equal contents have equal ids; matching ids implies equal contents (the format is lossless in both directions for the recognised branch vocabulary).
- **`bodies` is the canonical body set, not the sequence's planet count.** For a VEM cycler, `bodies = ("V","E","M")` is fixed; `sequence` varies. We do **not** sort `bodies` automatically here — the caller supplies the canonical order (alphabetical or `V<E<M` per spec convention; M4 expects callers to pass `("V","E","M")` in that order and tests assert this).
- **Per-leg revs / branches length matches `len(sequence)` for closed-loop cells.** For an open sequence with N encounters and N−1 legs, the tuples have length N−1; for closed-loop representation where `sequence[-1] == sequence[0]` and the cell describes a complete cycle, the tuples have length N−1 too (the last leg of length N−1 closes the loop). M4 documents this convention and tests it.
- **No canonicalisation of `sequence`.** Per spec §16.2 the canonical sequence is the lexicographically-minimal rotation; M4 deliberately produces sequences in *some* rotation (the natural one from the enumeration loop). M7's catalogue canonicalises at signature time. **This is the explicit responsibility split — do not move it into M4.**
- **Frozen + hashable.** Tuples are hashable; `Cell` is therefore hashable and usable as a `dict` key or `set` member — which the M7 ledger will need but M4 does not. We pay nothing to provide it now.

#### 3.1.2 `enumerate_cells` — the unpruned generator

```python
from collections.abc import Iterator

def enumerate_cells(
    body_set: tuple[str, ...],
    L_max: int,
    k_max: int,
    N_max: int,
    branch_set: tuple[str, ...] = ("single",),
) -> Iterator[Cell]:
    """Yield all combinatorially possible cells under the given caps.

    Per spec §13.2: under fixed caps (``L_max``, ``k_max``, ``N_max``) the
    cell set is finite and fully enumerable. This generator walks that
    set in lexicographically-stable order:

      for L in range(2, L_max + 1):            # min length 2 = one leg
        for k in range(1, k_max + 1):
          for sequence in body_set^L (each adjacent pair distinct):
            for per_leg_revs in {0..N_max}^(L-1):
              for per_leg_branch in branch_set^(L-1):
                yield Cell(...)

    No pruning. ``feasible_cells`` is the Tisserand-pruned variant.

    Parameters
    ----------
    body_set:
        Canonical body codes the sequence can visit. Cell.bodies is set
        to this tuple verbatim.
    L_max:
        Maximum number of encounters in the sequence (≥ 2).
    k_max:
        Maximum period in synodic multiples (≥ 1).
    N_max:
        Maximum heliocentric revolutions per leg (≥ 0). ``N_max = 0``
        restricts to direct legs (the M1/M3 only-supported Lambert
        regime); ``N_max ≥ 1`` requires the M4 caller to also expand
        ``branch_set`` to include ``"low"``/``"high"``.
    branch_set:
        Allowed branches per leg. Default ``("single",)`` matches the
        M1 Lambert solver's current capability. Callers exploring
        multi-rev branches once those land in Lambert will pass
        ``("single", "low", "high")``.

    Yields
    ------
    Cell
        Every combinatorial cell. Caller responsible for pruning.
    """
```

**Design notes:**

- **Sequence adjacency constraint.** Two *consecutive* bodies in a flyby sequence must differ — a flyby of body X immediately followed by another flyby of body X (with no intervening heliocentric arc) is degenerate. The enumerator skips such candidates at generation time, not as a post-filter, to keep the yielded count honest.
- **The yielded count is bounded above by `(L_max−1) · k_max · |body_set| · (|body_set|−1)^(L−1) · (N_max+1)^(L−1) · |branch_set|^(L−1)`.** For the gate test (§4.4) with `body_set=("E","M"), L_max=4, k_max=2, N_max=0, branch_set=("single",)`, the count works out to the documented number.
- **`per_leg_revs` consistency with `branch_set`.** A 0-rev leg must have `branch="single"`; a ≥1-rev leg must have `branch ∈ {"low","high"}`. The enumerator enforces this — when generating `per_leg_revs[i] == 0` it pins `per_leg_branch[i] = "single"`, and for `per_leg_revs[i] ≥ 1` it cycles through `branch_set ∩ {"low","high"}`. If `branch_set` does not include `"low"`/`"high"` but `N_max ≥ 1`, the multi-rev branches yield nothing — documented behaviour.
- **Iterator, not list.** `Iterator[Cell]` per typing; the body-set / sequence count explodes combinatorially with `L_max` and `|body_set|`, so a list-return would blow memory. Tests that need a count materialise locally via `list(...)`.

#### 3.1.3 `tisserand_feasible` — the M4 binding gate (spec §13.3)

```python
def tisserand_feasible(
    cell: Cell,
    vinf_cap: float,
    ephem: Ephemeris | None = None,
) -> bool:
    """Tisserand pruning per spec §13.3.

    True iff for every consecutive body pair ``(cell.sequence[i],
    cell.sequence[i+1])`` there exists a ``V∞ ≤ vinf_cap`` at which the
    two bodies are linkable (i.e. their constant-V∞ contours intersect
    in (a, e) space).

    This is the M4 enumerator's energetic feasibility predicate. The
    vast majority of cells die here, so compute is spent only on
    viable structures (spec §13.3).

    Parameters
    ----------
    cell:
        The cell to test.
    vinf_cap:
        Common V∞ ceiling, km/s. Each consecutive pair must be linkable
        somewhere in (0, vinf_cap].
    ephem:
        Currently unused (coplanar Tisserand needs only constants, not
        a planet state). Accepted as a parameter for API symmetry with
        ``feasible_cells`` and forward-compatibility with M6's
        ephemeris-aware variant.

    Returns
    -------
    bool
        ``True`` iff every consecutive pair has at least one linkable V∞
        below the cap.

    Notes
    -----
    Sampling resolution is fixed at 24 V∞ values evenly spaced in
    ``(0.5, vinf_cap]`` (the lower bound avoids the always-false
    asymptote near V∞ → 0). This is coarser than
    :func:`cyclerfinder.search.tisserand.linkable_region` would give
    but sufficient for the gate: a real linkable pair has a *band* of
    linkable V∞, not a discrete value. Tests pin this resolution.
    """
```

**Design notes (binding):**

- **Implementation is a straight loop over `tisserand.linkable(body_a, body_b, vinf)`** for `vinf` in a fixed 24-sample grid on `(0.5, vinf_cap]`. The grid is internal; if a future tuning is needed, change it in one place and rerun the gate test.
- **Coplanar-only.** Inherits the M2 `linkable` predicate's coplanar restriction. The M2 plan §3.2 hand-off explicitly preserved this; the M4 plan re-states it because the same gate must hold for the VEM enumeration in M8.
- **Never raises.** Mirrors `tisserand.linkable`'s contract — failures (NaN, no-bracket) return False, not exceptions. The enumerator must be able to call this in a tight loop without try/except scaffolding.
- **`ephem` parameter unused in M4** but accepted so the signature is stable when M6 introduces 3D Tisserand pruning (where the ephemeris is actually consulted for the body's instantaneous state). Marked `Ephemeris | None = None` to keep callers simple.

#### 3.1.4 `feasible_cells` — convenience composition

```python
def feasible_cells(
    body_set: tuple[str, ...],
    L_max: int,
    k_max: int,
    N_max: int,
    vinf_cap: float,
    ephem: Ephemeris | None = None,
    branch_set: tuple[str, ...] = ("single",),
) -> Iterator[Cell]:
    """Cells from :func:`enumerate_cells` filtered by :func:`tisserand_feasible`.

    The standard M4 entry point for any caller that wants the
    Tisserand-pruned cell stream (i.e. most callers — the unpruned
    ``enumerate_cells`` is exposed mostly for tests and for M7's ledger
    which counts pruned-vs-searched separately).
    """
    return (
        cell
        for cell in enumerate_cells(body_set, L_max, k_max, N_max, branch_set)
        if tisserand_feasible(cell, vinf_cap, ephem)
    )
```

#### 3.1.5 `deepening_frontier` — the §13.2 iterative-deepening generator

```python
def deepening_frontier(
    body_set: tuple[str, ...],
    ephem: Ephemeris | None = None,
    *,
    vinf_cap: float,
    L_step: int = 1,
    k_step: int = 1,
    N_step: int = 1,
    L_initial: int = 3,
    k_initial: int = 1,
    N_initial: int = 0,
    branch_set: tuple[str, ...] = ("single",),
) -> Iterator[Cell]:
    """Iterative-deepening frontier per spec §13.2.

    Yields cells in monotonically increasing complexity by raising
    the caps stepwise. After exhausting all feasible cells at
    ``(L_initial, k_initial, N_initial)``, raises the leading cap and
    continues, yielding only the *newly-added* cells at each tier.
    The infinite outer loop is the spec §13.2 "run forever" frontier;
    in M4 this generator is lazy and the caller stops it (typical
    pattern: ``itertools.islice``).

    M4 implements the **generator only**. The work-queue + ledger
    bookkeeping (spec §13.6, §13.8) — atomic claim/skip across workers,
    persistent dedup, prioritisation queues — are M7's responsibility.
    Calling this in M4 yields a finite-effort cell stream suitable for
    test exploration and for M5 to consume during optimisation
    development; it is **not** the production frontier.
    """
```

**Design notes:**

- **The newly-added set per tier.** When the caps grow from `(L, k, N)` to `(L+L_step, k, N)`, the new cells are those with `L < length(sequence) ≤ L+L_step` (the prior tier already covered shorter sequences). The generator computes this by passing `L_max = new_L` to `feasible_cells` and filtering out cells with `len(sequence) ≤ L`. Same idea for `k` and `N`. Implementation uses a small bookkeeping `set` of previously-yielded cell ids; M4 keeps this in-memory only.
- **No prioritisation in M4.** Spec §13.7 calls for VEM-first / long-period-first ordering; M4 does not implement this. The work-queue ordering is M7's concern (the ledger + queue is what makes prioritisation meaningful — without persistence, a heuristic order is a one-shot choice). M4's `deepening_frontier` walks `(L, k, N)` in their natural product order.
- **Why ship it at all in M4.** The spec §13.2 frontier *is* the central design idea. Shipping a single-process generator now means M5 can develop the optimiser against a real cell stream (instead of hand-crafted cells) without waiting for M7's ledger. The contract is: it yields cells, eventually visits all of them under monotonic cap growth, and is replaced — not extended — by M7's work-queue.

#### 3.1.6 API summary

| Symbol | Purpose | Notes |
|---|---|---|
| `Cell` | Frozen dataclass: the atomic search unit (§13.1) | `id` is a property, not a stored field |
| `Cell.id` | Spec §13.8 deterministic sortable identifier | `bodyset|sequence|kK|rrevs|bbranches` |
| `enumerate_cells(body_set, L_max, k_max, N_max, branch_set)` | All combinatorial cells under caps | Iterator; no pruning |
| `tisserand_feasible(cell, vinf_cap, ephem=None)` | Spec §13.3 pruning predicate | **M4 binding gate** |
| `feasible_cells(...)` | Pruned cell stream | `(c for c in enum if tiss_feas(c, vinf_cap))` |
| `deepening_frontier(body_set, ...)` | Spec §13.2 iterative-deepening generator | M7 replaces with ledger-backed work queue |

---

### 3.2 `model/score.py`

#### 3.2.1 The `Score` dataclass (spec §12(d) hard constraints + §5 step 4 metrics)

```python
from __future__ import annotations
from dataclasses import dataclass

@dataclass(frozen=True)
class Score:
    """Per-cycler ranking metrics, derived once and immutable.

    Mirrors spec §5 step 4 (rank seeds by total maintenance ΔV, max V∞,
    radial span, period error, taxi cost) and spec §12(d) (hard
    inequalities V∞ ≤ cap, r_p ≥ r_p_min, bend ≤ max at every flyby).

    All fields are plain floats / tuples for cheap serialisation later
    (the M7 catalogue record consumes this directly).
    """

    total_maintenance_dv_kms: float       # sum of flyby_dv at each encounter, km/s
    max_vinf_kms: float                   # max ||vinf|| over all encounters, km/s
    radial_span_au: tuple[float, float]   # (min perihelion, max aphelion) over legs, AU
    period_error_yr: float                # |cycler.period - target_period|, years
    taxi_cost_kms: float                  # surrogate hyperbolic-rendezvous V∞, km/s
    hard_constraints_pass: bool           # spec §12(d): V∞≤cap ∧ r_p≥r_p_min ∧ bend≤max
```

**Field semantics (binding):**

| Field | Definition | Source method |
|---|---|---|
| `total_maintenance_dv_kms` | `sum(flyby_dv_for(enc.body, enc.vinf_in, enc.vinf_out) for enc in cycler.encounters)` — the spec §6 powered-flyby decomposition summed across encounters. Replaces the M3 `Cycler.maintenance_dv()`'s naive velocity-discontinuity sum with the bend-and-magnitude-decomposed one. | `core.flyby.flyby_dv_for` |
| `max_vinf_kms` | `cycler.max_vinf()` directly from M3. | `Cycler.max_vinf` |
| `radial_span_au` | `cycler.radial_span()` directly from M3 — `(min perihelion, max aphelion)` across legs, in AU. | `Cycler.radial_span` |
| `period_error_yr` | `abs(cycler.period - target_period) / SECONDS_PER_YEAR`. The target is the `k * synodic_period` (or multi-body beat) the cell was constructed against; `score()` takes a `target_period_sec` argument so the caller can pass it explicitly. | computed in `score()` |
| `taxi_cost_kms` | Surrogate per §3.2.3 below. | `taxi_cost_kms()` |
| `hard_constraints_pass` | True iff at every encounter: `||vinf_in|| ≤ vinf_cap` AND `||vinf_out|| ≤ vinf_cap` AND the bend angle `is_ballistic_feasible(vin, vout, mu_planet, rp_min)` returns True. Per spec §12(d) hard inequalities; this is the v1.1 amendment binding M4 forward. | `is_ballistic_feasible`, magnitude check |

#### 3.2.2 `score()` — single-cycler reduction

```python
def score(
    cycler: Cycler,
    ephem: Ephemeris,
    vinf_cap: float,
    target_period_sec: float | None = None,
    rp_factors: dict[str, float] | None = None,
) -> Score:
    """Reduce a built cycler to a :class:`Score` record.

    Parameters
    ----------
    cycler:
        A ``Cycler`` produced by :func:`cyclerfinder.search.construct.construct_cycler`.
    ephem:
        Ephemeris used for sanity checks (not currently consumed but
        accepted for symmetry with future M6 ephemeris-mode scoring).
    vinf_cap:
        Hard ceiling on V∞ at every encounter, km/s. Spec §12(d).
    target_period_sec:
        The synodic-anchored target period the cell was built against.
        If ``None``, ``period_error_yr = 0.0`` (the caller has no target
        — typical for a hand-built test cycler). M4 callers from
        ``rank()`` always pass this.
    rp_factors:
        Per-body multiplier on ``SAFE_PERIHELION_KM[code]`` for the
        feasibility check. Default = ``1.0`` everywhere (use the
        constants-module default safe altitude). Lets callers tighten
        the flyby radius without editing the constants table.

    Returns
    -------
    Score
        Frozen, fully-populated.
    """
```

#### 3.2.3 `taxi_cost_kms()` — the surrogate (a real design decision)

Spec §2 mentions the "taxi" hyperbolic-rendezvous cost as one of the candidate metrics. The exact taxi-mission cost is a separate trajectory-design problem (it depends on the taxi's departure body, parking orbit, propulsion stage); we do **not** solve it here. We adopt a defensible **surrogate** that uses only the cycler's intrinsic quantities and is monotone in the underlying physical cost:

```python
def taxi_cost_kms(cycler: Cycler) -> float:
    """Surrogate hyperbolic-rendezvous cost for a taxi mission.

    Definition (M4):
        ``taxi_cost = max(||enc.vinf_in||)`` over the encounters whose
        body is in ``{"E"}`` (Earth-launched taxis must match Earth's
        V∞ at the cycler encounter to rendezvous from a low-Earth
        parking orbit + hyperbolic injection).

    Rationale: a taxi launched from Earth pays roughly the hyperbolic
    excess speed of the cycler at the Earth encounter as the
    propulsive component of its injection ΔV (plus a fixed escape
    component that is constant across candidates and therefore drops
    out of the *ranking*). This surrogate is **monotone in the actual
    cost**, which is what ranking needs.

    Caveats (documented to anchor future work):

    * Returns 0.0 if the cycler has no Earth encounter (e.g. a
      hypothetical Venus-Mars-only cycler). The caller is responsible
      for using a different surrogate in that case.
    * Ignores the taxi return leg (a real human-mission analysis cares
      about both directions; M4's purpose is candidate ranking, not
      mission costing).
    * Does **not** include the taxi's transfer-from-Earth-orbit ΔV,
      because that's body-of-departure dependent. The proper full
      taxi-trajectory cost is M5/M6+ territory if it ever becomes a
      ranking target.

    The choice of ``max(||vinf_in||) over Earth encounters`` (rather
    than ``mean``) reflects that the taxi must close the worst
    encounter, not the average — the bottleneck dominates.
    """
```

**Why this matters:** the spec lists "taxi cost" as a ranking metric. A poorly-chosen surrogate either ranks nothing useful (returns a constant) or ranks the wrong axis (correlates more with V∞_max than with the underlying mission cost). The max-over-Earth-V∞ surrogate is the minimal correct choice: it's nonzero, monotone in the real cost, and uses only what the `Cycler` already has. If M5/M8 ever needs a sharper number, this function is the single touch-point.

#### 3.2.4 `composite_score()` — single-number ordering

```python
DEFAULT_WEIGHTS: dict[str, float] = {
    "total_maintenance_dv_kms": 1.0,   # km/s; the primary cost
    "max_vinf_kms":             0.1,   # km/s; secondary penalty (spec §10 degeneracy guard)
    "period_error_yr":         10.0,   # yr; pin the cell's anchor period
    "taxi_cost_kms":            0.5,   # km/s; tertiary mission-relevance weight
}

def composite_score(s: Score, weights: dict[str, float] | None = None) -> float:
    """Weighted-sum scalar suitable for ascending sort (lower = better).

    Default weights documented above. Callers can override per-axis
    (e.g. a VEM-only campaign that downweights ``period_error_yr``
    because the 6.4-yr beat is itself approximate).

    Returns +inf if the cycler fails hard constraints, so any sort that
    mixes feasible and infeasible cyclers puts infeasibles last
    automatically.
    """
```

**Design notes:**

- **Weights are in the same units as their axes**, so the weighted sum is a unit-mixed scalar — purely ordinal, never reported as a "cost." Documented.
- **`hard_constraints_pass` is enforced here too**: failure returns `+inf` so the rank routine's sort places infeasibles last without an extra filter. The `rank` routine *also* filters them out; both layers exist as defence-in-depth.
- **Default weights are tunable, not load-bearing.** Tests pin a single deterministic ordering on a hand-constructed pair (one clearly better than the other), but the weights themselves are not the gate. M5/M8 will revisit.

#### 3.2.5 `rank()` — the spec §5 step 4 top-N reducer

```python
def rank(
    cyclers: list[Cycler],
    ephem: Ephemeris,
    vinf_cap: float,
    n_keep: int = 20,
    target_period_sec: float | None = None,
    weights: dict[str, float] | None = None,
) -> list[tuple[Score, Cycler]]:
    """Score every cycler, filter by ``hard_constraints_pass``, sort by
    composite ascending, return top ``n_keep`` as ``(Score, Cycler)`` pairs.

    Mirrors spec §5 step 4: "rank seeds by total maintenance ΔV, max
    V∞, radial span (must reach all target bodies), period error,
    taxi cost. Keep top-N." The composite ordering is documented in
    ``composite_score``.

    Empty input ⇒ empty output. Fewer feasible cyclers than ``n_keep``
    ⇒ output as long as the feasible set.
    """
```

**Design notes:**

- **Stable sort.** Ties on composite score break by cell id alphabetically (the M4 enumerator's `Cell.id` is sortable; we attach it post-hoc via the cycler's structure if we have one, or simply by the tuple's index in the input). The point is reproducibility — two M4 runs on the same input return identical output.
- **No early termination on quality.** `rank` is not the M5 optimiser; it reports what it gets. If `cyclers` is short, the top-N is short.
- **Empty input is fine, not an error.** The M5 optimiser may legitimately produce no feasible candidate; `rank([], ...)` returns `[]`.

#### 3.2.6 API summary

| Symbol | Purpose | Notes |
|---|---|---|
| `Score` | Frozen ranking metrics | Six fields per spec §5 step 4 + §12(d) |
| `score(cycler, ephem, vinf_cap, target_period_sec, rp_factors)` | Single-cycler reduction | Replaces M3 `Cycler.maintenance_dv()` for ranking purposes |
| `taxi_cost_kms(cycler)` | Surrogate per §3.2.3 | `max ||vinf_in|| over Earth encounters` |
| `composite_score(s, weights)` | Sortable scalar | Defaults in `DEFAULT_WEIGHTS`; returns `+inf` if hard constraints fail |
| `rank(cyclers, ephem, vinf_cap, n_keep, target_period_sec, weights)` | Top-N reducer | Empty input ⇒ empty output |
| `DEFAULT_WEIGHTS` | Module-level constant | Tunable; reported in M8 docs |

---

### 3.3 Imports / dependency graph after M4

```
constants.py    (M0)  ← root, no deps
ephemeris.py    (M1)  ← constants
lambert.py      (M1)  ← constants
kepler.py       (M1)  ← constants
flyby.py        (M2)  ← constants
frames.py       (M3)  ← constants
tisserand.py    (M2)  ← constants, flyby
resonance.py    (M2)  ← constants
model/cycler.py (M3)  ← constants, frames
construct.py    (M3)  ← constants, ephemeris, lambert, model/cycler
search/sequence.py (M4) ← search/tisserand                            [new]
model/score.py     (M4) ← constants, core/flyby, model/cycler         [new]
```

No cycles. `sequence.py` depends *only* on `tisserand.py` (it reuses `linkable()`) — it does **not** import `construct.py` (M4 enumerator does not build cyclers; that's M5's job). `score.py` depends on `flyby.py` and `cycler.py` — it scores what M3's constructor produced, with no awareness of how the cycler was constructed or enumerated.

This independence is deliberate: M4's two modules can be tested in isolation, and M5's optimiser will be the first thing that wires them together (cell → construct → score loop).

---

## 4. Tests + gate

Tests live under `tests/search/` (new directory) and `tests/model/`. Tolerances and constants are named module-level at the top of each test file.

### 4.1 `tests/search/test_sequence.py`

#### Cell dataclass

| Test | Assertion |
|---|---|
| `test_cell_is_frozen` | Assigning to a `Cell` field raises `FrozenInstanceError`. |
| `test_cell_id_format` | A hand-built `Cell(bodies=("V","E","M"), sequence=("E","V","M","E","M","E"), period_k=3, per_leg_revs=(0,0,1,0,1), per_leg_branch=("single","single","low","single","low"))` has `id == "VEM|E-V-M-E-M-E|k3|r00101|bssls s".replace(" ","")` — exact match to spec §13.8 example (with the branch alphabet mapping documented above). |
| `test_cell_id_uniqueness` | Two distinct cells from `enumerate_cells(("E","M"), 3, 1, 0)` have distinct ids; ids are stable across runs (a set of all ids has the expected size). |
| `test_cell_hashable` | A `Cell` is usable as a `dict` key (M7 ledger dependency). |

#### Enumeration counts (this is the §4.4 gate documentation)

| Test | Assertion |
|---|---|
| `test_enumeration_count_em_l2_k1` | `len(list(enumerate_cells(("E","M"), L_max=2, k_max=1, N_max=0, branch_set=("single",))))` equals **2**: two sequences of length 2 with distinct adjacent bodies (`E-M`, `M-E`), 1 leg each, 1 k value, 1 rev option (0), 1 branch ⇒ 2 × 1 × 1 × 1 = 2. |
| `test_enumeration_count_em_l4_k2` (**gate**) | `len(list(enumerate_cells(("E","M"), L_max=4, k_max=2, N_max=0, branch_set=("single",))))` equals **the documented number** computed in §4.4 below. |
| `test_enumeration_excludes_consecutive_same_body` | No yielded cell has `sequence[i] == sequence[i+1]` for any `i`. |
| `test_enumeration_iterator_not_list` | `type(enumerate_cells(...)).__name__` is a generator (`isinstance(it, Iterator)` is True; `isinstance(it, list)` is False). Catches accidental list-materialisation regressions that would blow memory at higher caps. |

#### Tisserand pruning gate (the **binding** spec §8 gate)

| Test | Assertion |
|---|---|
| `test_tisserand_pruning_rejects_low_vinf_em` (**gate**) | `tisserand_feasible(Cell(bodies=("E","M"), sequence=("E","M"), period_k=1, per_leg_revs=(0,), per_leg_branch=("single",)), vinf_cap=2.0)` is **False** — E–M is not linkable at any V∞ ≤ 2 km/s (this V∞ is well below the Aldrin family's ~5.5 km/s). |
| `test_tisserand_pruning_accepts_aldrin_vinf_em` (**gate**) | Same cell at `vinf_cap=8.0` is **True** — the Aldrin family is reachable in this V∞ band (M2 confirmed `linkable("E","M", 5.5) is True`). |
| `test_tisserand_pruning_propagates_through_sequence` | A 3-encounter cell `E-M-V` at `vinf_cap=8.0` requires *both* `(E,M)` and `(M,V)` to be linkable; a cell where `(M,V)` is not linkable returns False even if `(E,M)` is. Construct with a low cap (`vinf_cap=1.0`) to force false. |
| `test_tisserand_feasible_never_raises` | Parametrised over impossible / degenerate inputs (`vinf_cap=0.0`, `vinf_cap=-1.0`, empty sequence simulated by a single-encounter cell rejected at the dataclass level) — every call returns a `bool`, never raises. |

#### `feasible_cells` and `deepening_frontier`

| Test | Assertion |
|---|---|
| `test_feasible_cells_subset_of_enumerate` | `set(feasible_cells(("E","M"), 3, 1, 0, vinf_cap=8.0)) ⊆ set(enumerate_cells(("E","M"), 3, 1, 0))` (as `Cell.id` sets). |
| `test_feasible_cells_strict_subset_at_low_cap` | At `vinf_cap=1.0`, the feasible set is strictly smaller than the unpruned set (the gate is doing work). |
| `test_deepening_frontier_yields_in_increasing_complexity` | Take 20 cells from `deepening_frontier(("E","M"), vinf_cap=8.0)`; assert `len(cell.sequence)` is monotonically non-decreasing within each tier and never *decreases* across tiers. |
| `test_deepening_frontier_no_repeats` | The first 50 cells from the frontier have distinct ids. |

### 4.2 `tests/model/test_score.py`

#### Score dataclass

| Test | Assertion |
|---|---|
| `test_score_is_frozen` | Assignment raises `FrozenInstanceError`. |
| `test_score_field_types` | Type annotations match — `total_maintenance_dv_kms` is float, `radial_span_au` is `tuple[float, float]`, etc. (Light reflection check; mypy strict gives most of this for free.) |

#### Hard constraints (spec §12(d))

| Test | Assertion |
|---|---|
| `test_hard_constraints_pass_ballistic` | Hand-built `Cycler` with `vinf_in == vinf_out` at every encounter and `||vinf|| < vinf_cap` ⇒ `hard_constraints_pass is True`. |
| `test_hard_constraints_fail_on_vinf_cap` | Hand-build a `Cycler` with `||vinf|| = 12 km/s` at one encounter; `score(..., vinf_cap=7.0).hard_constraints_pass is False`. |
| `test_hard_constraints_fail_on_overbent_pair` | Hand-build a `Cycler` whose `vinf_in`/`vinf_out` at Mars exceed the max bend angle (M2 anchor: ~24° at 7 km/s); `hard_constraints_pass is False`. |
| `test_score_aldrin_passes_hard_constraints` (**gate**) | `cyc = build_aldrin_seed(ephem)`; `s = score(cyc, ephem, vinf_cap=12.0)`; `s.hard_constraints_pass is True` (the Aldrin V∞ values from M3 hand-off — 6.5 / 9.7 km/s — fit under a generous 12 km/s cap and within Mars/Earth max-bend at those V∞). Also `s.max_vinf_kms == pytest.approx(9.7, abs=0.1)`. |

#### Taxi cost

| Test | Assertion |
|---|---|
| `test_taxi_cost_earth_only` | A cycler with Earth encounters at V∞ = 5 and 7 km/s and no other bodies → `taxi_cost_kms == 7.0`. |
| `test_taxi_cost_zero_when_no_earth` | A hypothetical V-M-only cycler → `taxi_cost_kms == 0.0`. |
| `test_taxi_cost_aldrin` | `taxi_cost_kms(build_aldrin_seed(ephem))` equals the Aldrin Earth V∞ from M3 hand-off (~6.5 km/s) within 0.1. |

#### Composite + rank

| Test | Assertion |
|---|---|
| `test_composite_finite_and_reproducible` | Two calls on the same `Score` return bitwise-identical floats. |
| `test_composite_infinite_on_hard_constraint_fail` | A `Score` with `hard_constraints_pass=False` ⇒ `composite_score(s) == math.inf`. |
| `test_rank_orders_low_dv_first` | Three hand-built cyclers with ΔV = 0.1, 0.5, 1.0 km/s and other axes equal ⇒ `rank(...)` returns them in that order. |
| `test_rank_filters_infeasible` | A list of three cyclers (one infeasible by V∞ cap) ⇒ `rank(..., n_keep=3)` returns 2 items. |
| `test_rank_empty_input` | `rank([], ephem, vinf_cap=7.0) == []`. |
| `test_rank_n_keep_truncation` | A list of 10 feasible cyclers + `n_keep=3` ⇒ output length 3. |

### 4.3 Tolerance summary

| Layer | Quantity | Tolerance |
|---|---|---|
| Cell id format | exact equality | exact string match |
| Enumeration counts | exact equality | integer count |
| Tisserand pruning at `vinf_cap=2.0` for E–M | bool | strict `False` |
| Tisserand pruning at `vinf_cap=8.0` for E–M | bool | strict `True` |
| Aldrin score max_vinf_kms | km/s | ±0.1 |
| Aldrin score taxi_cost_kms | km/s | ±0.1 |
| Composite reproducibility | float | bitwise (no `pytest.approx`) |

### 4.4 Documented enumeration count for the gate

`enumerate_cells(body_set=("E","M"), L_max=4, k_max=2, N_max=0, branch_set=("single",))` produces:

- For `L=2` (2 encounters, 1 leg): adjacency-distinct sequences over `("E","M")` of length 2 = `{E-M, M-E}` → 2 sequences. × 2 k values × 1 rev option (0) × 1 branch = **4 cells**.
- For `L=3` (3 encounters, 2 legs): sequences of length 3 with no consecutive repeats = `{E-M-E, M-E-M}` → 2 sequences. × 2 k × 1^2 revs × 1^2 branches = **4 cells**.
- For `L=4` (4 encounters, 3 legs): `{E-M-E-M, M-E-M-E}` → 2 sequences. × 2 k × 1^3 × 1^3 = **4 cells**.

Total: **12 cells**. The gate test pins this number; if the enumeration algorithm changes (e.g. introducing a different adjacency rule, or counting closed-loop cells separately) the number changes and the test fails informatively.

For a body set of size `|B|`, length `L`, the count of adjacency-distinct sequences is `|B| × (|B|−1)^(L−1)`. For `|B|=2, L∈{2,3,4}` this gives `2, 2, 2` respectively — matching the breakdown above.

---

## 5. Risks

| # | Risk | Likelihood | Impact | Mitigation in M4 |
|---|---|---|---|---|
| 1 | **`linkable` false-negative bites the pruning gate** | medium | **high** — would discard real cyclers before they're built | Inherited M2 risk; M4 chooses a 24-sample grid on `(0.5, vinf_cap]` (denser than `linkable_region`'s default 50 across a wider range, so per-band density is comparable) and tests against the Aldrin anchor at 8 km/s. If M4 tests catch a regression in `tisserand.linkable`, fix it in M2's module, not by softening M4's gate. |
| 2 | **Enumeration count drift** | low | medium — gate test fails | The §4.4 derivation is exact; if the adjacency rule changes intentionally, update §4.4 and the test in the same commit. Otherwise the failing test is informing of an unintended algorithm change. |
| 3 | **`Cell.id` format collides with M7's signature** | low | medium — would force a rewrite of the ledger key scheme | `Cell.id` is a *cell* identifier (discrete structure); spec §16.2's signature is a *cycler* identifier (built result with V∞ multiset + leg elements). They have different sources of truth and live at different levels; M7 will store both per record. Documented in `Cell.id`'s docstring. |
| 4 | **Multi-rev branch enumeration combinatorial explosion** | medium | medium | At `branch_set=("single","low","high")` and `N_max=2`, the per-leg multiplier becomes 5 (1 single-rev + 2 each for 1- and 2-rev). For `L=6` that's `5^5 = 3125`× the single-branch count. M4's iterative-deepening generator is lazy; tests do not enumerate at high caps. Document the explosion in `enumerate_cells`'s docstring with a worked example. |
| 5 | **`taxi_cost_kms` surrogate misleads ranking** | medium | medium — would prioritise the wrong axis | Documented in §3.2.3 as a surrogate with explicit caveats; tested against the Aldrin anchor (where the Earth V∞ ≈ 6.5 km/s and the surrogate should equal that). M5/M6 has the right context to revise if needed. |
| 6 | **Hard-constraint check needs `mu_planet` and `rp_min` per body** | low | low | Sourced from `constants.PLANETS` and `constants.SAFE_PERIHELION_KM` (the same source-of-truth M2's `flyby_dv_for` already uses). `rp_factors` argument lets callers tighten the cap without editing constants. |
| 7 | **`score()` collides with `Cycler.maintenance_dv()` semantics** | medium | medium — silent disagreement | M3's `Cycler.maintenance_dv()` returns the naive velocity-discontinuity sum (no bend feasibility); M4's `Score.total_maintenance_dv_kms` uses `flyby_dv_for` (with the bend-and-magnitude decomposition). They will disagree on infeasible flybys. Documented in `score()`'s docstring; M3's method is **not** removed (M3 tests depend on it). |
| 8 | **VEM beat period is not a single `k * synodic`** | medium | low | `target_period_sec` in `score()` is a free parameter — the caller passes whatever target the cell was built against (single-pair k·synodic for E-M cells, the §16.1 beat for VEM). M4 doesn't compute this internally; that's the caller's job. |
| 9 | **`deepening_frontier` in-memory dedup grows unbounded** | low (in M4) | low (in M4) | M4's frontier is for test exploration only. Tests bound iteration to `itertools.islice(..., 50)`. M7's ledger replaces the in-memory dedup with persistent storage; the M4 frontier is *not* the production runner. |
| 10 | **Score weights treated as load-bearing** | medium | medium — would lock the optimiser's objective to one M4 author's intuition | `DEFAULT_WEIGHTS` is documented as tunable; tests pin one *ordering* of hand-built scores, not the absolute composite values. M5/M8 will revisit. |

---

## 6. Dependency additions

**None.** M4 uses only `numpy` (M0), `scipy` (M1 — already a runtime dep, consumed transitively by `tisserand.linkable`), the standard library, and the in-house M0–M3 modules. No edits to `pyproject.toml`; no `uv.lock` regeneration needed. Confirmed by `grep dependencies pyproject.toml`: current runtime deps are `numpy>=2.0` and `scipy>=1.13`, both already present.

---

## 7. Order of work

The `todo.md` mirrors this with checkboxes.

1. **Re-read predecessor docs.** Confirm M3's hand-off note (Aldrin V∞ values, `closure_residual` semantics, multi-rev branch usage observed). Confirm M2's `linkable` signature is `linkable(body_a, body_b, vinf_kms, tol_au=..., tol_e=...)`. Confirm M3's `Cycler` exposes `max_vinf()`, `radial_span()`, `encounters[i].vinf_in/out`.
2. **Write `search/sequence.py` skeleton.** Just the `Cell` dataclass + signatures of the four functions, with `NotImplementedError` bodies. Run `mypy --strict` — confirm the type signatures are clean before any logic lands.
3. **Implement `Cell.id` and test it** against the spec §13.8 worked example. `test_cell_id_format` should pass with no other code yet.
4. **Implement `enumerate_cells`.** Walk the nested loop, enforce adjacency distinctness, yield `Cell`. Land `test_enumeration_count_em_l2_k1`, `test_enumeration_count_em_l4_k2`, `test_enumeration_excludes_consecutive_same_body`, `test_enumeration_iterator_not_list`.
5. **Implement `tisserand_feasible`.** Sample `(0.5, vinf_cap]` at 24 points; call `tisserand.linkable` per pair; return True at first linkable V∞ found per pair, False if any pair has no linkable V∞ in the band. Land `test_tisserand_pruning_rejects_low_vinf_em` and `test_tisserand_pruning_accepts_aldrin_vinf_em` (the **binding gate**).
6. **Implement `feasible_cells` and `deepening_frontier`.** Land their tests including the subset assertion and the no-repeats / monotonicity properties.
7. **Write `model/score.py` skeleton.** `Score` dataclass + signatures + `NotImplementedError` bodies. mypy strict pass.
8. **Implement `taxi_cost_kms`** and its tests (Aldrin anchor, Earth-only, no-Earth cases).
9. **Implement `score()`.** Sum `flyby_dv_for` per encounter; compute `radial_span()` and `max_vinf()` from `Cycler`; compute `period_error_yr` from the optional target; evaluate hard constraints. Land `test_hard_constraints_pass_ballistic`, `test_hard_constraints_fail_on_vinf_cap`, `test_hard_constraints_fail_on_overbent_pair`, `test_score_aldrin_passes_hard_constraints` (**gate**).
10. **Implement `composite_score` and `rank`.** Land the ordering, filter, truncation, empty-input tests.
11. **Add `Score` to `model/__init__.py`** re-exports alongside `Cycler`, `Leg`, `Encounter`.
12. **Run the full local quality gate:** `uv run pytest`, `uv run ruff check .`, `uv run ruff format --check .`, `uv run mypy src tests`.
13. **Commit** as `m4: cell enumeration with Tisserand pruning + per-cycler scoring`. Push; confirm CI green.
14. **Update `docs/overview.md`** §4 milestone table: M4 → completed; M5 → planned.
15. **Hand-off note** appended to `todo.md` under `## Hand-off to M5`.

The order is "Cell + enumerator → Tisserand gate → score → rank" deliberately: each new item depends only on what came before it, and the M4 gate (the Tisserand pruning + the Aldrin score) lands as soon as the supporting infrastructure is in place.

---

## 8. Exit checklist (the gate, restated)

Before declaring M4 done:

- [ ] `uv run pytest tests/search/test_sequence.py` green; the four gate tests in §4.1 all pass at the documented tolerances.
- [ ] `uv run pytest tests/model/test_score.py` green; `test_score_aldrin_passes_hard_constraints` passes.
- [ ] `uv run pytest` green overall (no regression of M0–M3 tests).
- [ ] `uv run ruff check .` clean.
- [ ] `uv run ruff format --check .` clean.
- [ ] `uv run mypy src tests` clean under `strict = true` — including the new `Cell` and `Score` dataclasses and the `Iterator[Cell]` generators.
- [ ] CI green on the M4 commit.
- [ ] `docs/overview.md` updated: M4 status = `completed`; M5 row marked `planned`.
- [ ] `## Hand-off to M5` section appended to `phases/m4-enumeration-scoring/todo.md` covering:
  - The exact enumeration count produced for the L=4, k=2, single-branch E-M case (and any deviation from §4.4's 12-cell prediction with rationale).
  - The observed Aldrin `Score` values (max_vinf, taxi_cost, composite under default weights) for M5 to use as a regression anchor when the optimiser starts running cells through `construct → score → rank`.
  - Whether `deepening_frontier`'s in-memory dedup proved adequate at the test caps, or whether M7's ledger work needs to happen sooner.
  - Any behavioural surprises in `tisserand.linkable` discovered during the 24-sample-grid sweep that M2 might want to know about.

(Writing the M5 plan doc is the first task of M5, not an M4 exit criterion.)
