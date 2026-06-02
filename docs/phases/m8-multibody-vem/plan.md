# M8 — Multi-body (3-body) VEM cyclers: beat-period dispatch, N≥3 enumeration, rediscovery gate

> Written with the `superpowers:writing-plans` conventions. This is a
> **task-level TDD plan** (write failing test → run red → minimal impl → run
> green → commit). It is the detailed expansion of the **M-N** milestone in
> `docs/superpowers/plans/2026-06-02-multirev-3d-vem-ephemeris-roadmap.md`
> (which calls it "N-encounter + VEM enumeration & rediscovery, deps: M-L").
> Read that roadmap's "Verified current state" section first — it is the
> authoritative survey of what is already built.

---

## Goal

Make heliocentric **Venus–Earth–Mars (VEM)** triple cyclers reachable by the
existing optimiser/enumerator stack, and admit the two blocked VEM catalogue
rows (`jones-2017-vem-triple-family`, `vem-emeeve-3syn`, both currently
`ExclusionReason.NOT_TWO_BODY`, task `#76`-adjacent) into the rediscovery
census. The single concrete code blocker is named in
`src/cyclerfinder/search/optimize.py:229` `_target_period_sec` — its docstring
says the 2-body single-pair period formula "is no longer correct" for
`len(cell.bodies) >= 3` and points at
`cyclerfinder.search.resonance.multi_body_beat_days` /
`beat_period_days` as the extension point. This plan dispatches there.

**Honesty boundary (carried from spec §3, §11.3, §17).** VEM *strict* ballistic
closure in the circular-coplanar model is genuinely open research, not a
solved reproduction. Spec line 199: *"frame M8 as 'search + report best
candidates,' not 'guaranteed novel cycler.'"* This plan therefore delivers
**(a)** the period/enumeration/feasibility plumbing for N≥3 and **(b)** a
VEM rediscovery *gate that is validated only against sourced anchors that
actually exist in the catalogue* — which for both VEM rows is **period and
sequence only** (every `vinf_kms` and every `orbit_elements` field on both
rows is `null`; see `data/catalogue.yaml:1677-1703` and `1788-1806`). The
gate does **not** assert a converged ballistic V∞ signature, because no such
sourced number exists to assert against, and a value our own optimiser
computes can never be the EXPECTED side of a golden test (project memory:
*golden tests use sourced expected values only*). The optimiser-convergence
case is wired as a documented `xfail` to be flipped by **M-ED**
(real-ephemeris discovery), exactly as the roadmap's M-N test-gate row states.

---

## Architecture

The N≥3 compute stack is **already arbitrary-N** (verified read-only
2026-06-02, recorded in the roadmap doc):

- `construct_cycler` (`search/construct.py:40`) validates only `n>=2` and
  monotonic times — body-agnostic, multi-rev-capable via
  `max_revs_per_leg`/`branch_per_leg`.
- `Cycler`/`Encounter`/`Leg` (`model/cycler.py`) hold `len(legs)==len(bodies)-1`;
  closure uses `encounters[0]`/`[-1]` — N-agnostic.
- `enumerate_cells` (`search/sequence.py:145`) already loops
  `for length in range(2, l_max + 1)` over any `body_set`.
- `optimise_cell_idealized` (`search/optimize.py:960`) free-variable count is
  `n_interior = len(cell.sequence) - 2` (`optimize.py:1053`) — N-agnostic.

Only **three** things are hardcoded-to-2-body or coplanar-limited, and M8
fixes exactly those, in dependency order:

1. **Period.** `_target_period_sec` (`optimize.py:229-247`) uses
   `cell.bodies[0]`/`cell.bodies[1]` only → wrong total period T for VEM.
   Fix: dispatch to the beat-period resonance functions for `len(bodies) >= 3`,
   keep the 2-body fast path bit-identical.
2. **Feasibility.** `tisserand_feasible` (`sequence.py:265`) is honest but
   **coplanar (i=0) only** — it inherits `tisserand.linkable`'s restriction
   (`sequence.py:308-309`). It already iterates over *every* consecutive
   sequence pair, so it works structurally for VEM sequences today; M8 only
   adds tests pinning that behaviour and documents the coplanar caveat at the
   VEM-cell boundary. **No lift of the coplanar assumption here** — that is
   M-3D, deliberately out of M8 scope (see §6).
3. **Loader admission.** `tests/_catalogue_loader.py` classifies
   `len(bodies) != 2` as `NOT_TWO_BODY` (`_catalogue_loader.py:174-176`) and
   `_is_two_body_alternation` (`:134-152`) rejects N-encounter sequences. M8
   adds a *VEM-aware* classification path so the two VEM rows become
   `CONSTRUCTIBLE_MULTIBODY` (a new reason, see §5) carrying the sourced
   period anchor, and the frozen census ratchet
   (`test_catalogue_rediscovery.py:267` `EXPECTED_COVERAGE`) is updated in the
   same commit.

The beat machinery itself is **already built and tested**
(`search/resonance.py:140` `multi_body_beat_days`, `:191` `beat_period_days`;
`tests/search/test_resonance.py:44-64` pins `(4,3)` ⇒ 6.406 yr for `["V","E","M"]`).
M8 *consumes* it; it does not reimplement it.

**Critical-path note.** The roadmap marks M-N "ready once M-L (multi-rev
Lambert) lands" because VEM closure sequences (e.g. EMEEVE) contain same-era
loop legs that need `n_revs >= 1`. This plan's §2 (period) and §3 (Tisserand)
and §4 (the period-handling wiring) are **independent of M-L** and can be
executed now. §5 (loader admission as period+sequence-only) is also M-L-free.
The only M-L-dependent piece is converging an actual VEM cycler geometry,
which this plan keeps as `xfail` and hands to M-ED. Tasks are ordered so the
M-L-free work lands first.

---

## Tech stack

Python 3.11, numpy, scipy (`differential_evolution` + SLSQP via the existing
`optimise_cell_idealized` path), pytest. uv-managed venv (no pip). Lint/type
gate before every commit: `uv run ruff check .`, `uv run ruff format --check .`,
`uv run mypy src tests` (the M7 plan's gate item 11 is the standing bar).
Run the fast suite with `uv run pytest -m "not slow"`; the parametrised
rediscovery suite is `slow`.

---

## Spec references

- §1 (line 12) — VEM is a primary target ("under-explored VEM triple-cycler space").
- §3 (lines 41, 44) — VEM beat ≈ 6.4 yr (3×E–M ≈ 4×E–V); Venus is the strong
  steerer (~61° at 7 km/s), Mars weak (~22°). Sourced bend numbers reconciled in §9.1.
- §4 (line 61) — `resonance.py` slot: "synodic periods; k-synodic candidates;
  multi-body beat finder" (already shipped at M2).
- §5 step 1 (line 87) — resonance sets the candidate total periods T.
- §8 (line 152) — **M8 milestone definition:** "VEM campaign + UX: run the
  enumerator on `[Venus, Earth, Mars]` at the 6.4-yr beat; `cli`, `viz`,
  reporting, docs." (CLI/viz are a *separate* M8 sub-stream; see §6 scope note.)
- §9 (line 160) — gate anchor: "VEM beat 3×E–M ≈ 4×E–V ≈ 6.40 yr." This is the
  only sourced VEM number this plan asserts against.
- §11.3 / §17 (line 199) — "frame M8 as 'search + report best candidates,' not
  'guaranteed novel cycler.' Any closing VEM solution must pass GMAT
  verification and catalogue check before any novelty claim."
- §13.1, §13.7 (lines 354, 384) — order the frontier "VEM and other multi-body
  cells before pure Earth–Mars."
- §13.8 (line 363) — the worked cell-id example is itself a VEM cell:
  `VEM|E-V-M-E-M-E|k3|r00101|blllll`.

---

## 1. What this milestone delivers

M8 is **additive** to one source file and one test-support file, plus new
tests:

- `src/cyclerfinder/search/optimize.py` — `_target_period_sec` gains a
  multi-body dispatch branch (§2). No signature change; the 2-body path is
  byte-for-byte preserved.
- `tests/_catalogue_loader.py` — gains a VEM-aware classification path and a
  new `ExclusionReason.CONSTRUCTIBLE_MULTIBODY` (§5).
- New tests: `tests/search/test_optimize_multibody.py` (§2, §4),
  `tests/search/test_sequence_multibody.py` (§3),
  `tests/test_vem_rediscovery.py` (§4 gate, §5).
- `tests/test_catalogue_rediscovery.py::EXPECTED_COVERAGE` — census ratchet
  update (§5).

No new production module. No edits to `resonance.py`, `construct.py`,
`sequence.py` production bodies in the M-L-free portion (the Tisserand work in
§3 is **tests only** — it pins existing behaviour and documents the coplanar
caveat; it does not change `tisserand_feasible`).

---

## 2. `_target_period_sec` multi-body dispatch (M-L-free, first)

The blocker. `optimize.py:229-247` currently:

```python
def _target_period_sec(cell: Cell) -> float:
    if len(cell.bodies) < 2:
        raise ValueError(...)
    body_a = cell.bodies[0]
    body_b = cell.bodies[1]
    t_syn_days = synodic_period_days(body_a, body_b)
    return t_syn_days * cell.period_k * SECONDS_PER_DAY
```

For VEM `cell.bodies == ("V","E","M")`, `cell.period_k == 3`, this returns
`3 * T_syn(V,E)` ≈ `3 * 0.799 yr` ≈ 2.4 yr — **wrong**; the EMEEVE archetype's
sourced period is 6.41 yr (`catalogue.yaml:1782`), the 3×E-M = 4×E-V beat.

### Design decision (settle before coding)

For `len(bodies) >= 3`, the *total cycler period* is `period_k` multiples of
the **beat period** of the body set, where the beat is the smallest tuple from
`multi_body_beat_days(bodies, k_max=...)` fed to `beat_period_days`. But note
the catalogue's `period.k` for VEM is anchored on a **pair** (`period.pair:
"E-M"`, `period.k: 3`, `years: 6.41` — `catalogue.yaml:1780-1782`), not on a
beat tuple. Two consistent readings exist:

- **Reading A (beat-multiple):** `T = period_k_beat * beat_period_days(bodies, top_tuple)`
  with `period_k_beat = 1` for the natural beat. This is the resonance-pure
  interpretation matching `resonance.py`'s own docstring (line 7).
- **Reading B (anchor-pair-multiple):** `T = cell.period_k * T_syn(period_pair)`
  where `period_pair` is the catalogue's anchor pair `(E,M)`. For the EMEEVE
  row `3 * T_syn(E,M) = 3 * 2.135 = 6.405 yr` — the same number, because the
  beat *is defined* as the commensurability where `3*T_syn(E,M) ≈ 4*T_syn(E,V)`.

These agree numerically at the natural beat (that is the whole point of a
beat). **Choose Reading A** for `_target_period_sec` (a `Cell` carries no
"anchor pair" field — only `bodies` and `period_k`), and interpret
`cell.period_k` as the *beat multiple*. The catalogue loader (§5) is
responsible for translating a catalogue row's anchor-pair `period.k=3` into the
`Cell.period_k=1` beat-multiple it builds — this keeps `_target_period_sec`
dependent only on `Cell` fields and keeps the catalogue's published
anchor-pair convention intact at the loader boundary. Document this mapping in
both places.

> **Open question for the human (do not guess):** confirm Reading A and the
> "loader converts anchor-pair-k to beat-multiple-k" split before coding. The
> alternative (add a `period_pair`/`period_basis` field to `Cell`) is a wider
> change touching M4's frozen dataclass and the `Cell.id` format
> (`sequence.py:105`); flag for review rather than do it unilaterally.

### Task 2.1 — failing test: VEM 3-body period resolves to the 6.41-yr beat

Write `tests/search/test_optimize_multibody.py`:

```python
"""M8 multi-body period dispatch tests (spec §3, §9).

The only sourced VEM number we assert against is the beat period
3 x E-M ~ 4 x E-V ~ 6.40 yr (spec §9 line 160). Cross-checked against the
already-passing resonance gate tests/search/test_resonance.py:57-64
(top tuple (4,3) -> 6.406 yr) — same physics, exercised here through the
optimiser's period resolver.
"""
from __future__ import annotations

import pytest

from cyclerfinder.core.constants import DAYS_PER_JULIAN_YEAR, SECONDS_PER_DAY
from cyclerfinder.search.optimize import _target_period_sec
from cyclerfinder.search.sequence import Cell


def _vem_cell(period_k: int = 1) -> Cell:
    # Beat-multiple convention (Reading A): period_k counts beat periods.
    seq = ("E", "V", "M", "E", "M", "E")  # spec §13.8 worked VEM sequence
    n_legs = len(seq) - 1
    return Cell(
        bodies=("V", "E", "M"),
        sequence=seq,
        period_k=period_k,
        per_leg_revs=(0,) * n_legs,
        per_leg_branch=("single",) * n_legs,
    )


def test_vem_target_period_is_natural_beat_6_406yr() -> None:
    """3-body VEM cell at period_k=1 resolves to the 6.40-yr beat (spec §9)."""
    t_sec = _target_period_sec(_vem_cell(period_k=1))
    t_yr = t_sec / SECONDS_PER_DAY / DAYS_PER_JULIAN_YEAR
    assert t_yr == pytest.approx(6.406, abs=0.01), f"got {t_yr:.4f} yr"


def test_vem_target_period_scales_with_period_k() -> None:
    """period_k=2 doubles the beat (the 12.8-yr branch, spec §3 line 41)."""
    t1 = _target_period_sec(_vem_cell(period_k=1))
    t2 = _target_period_sec(_vem_cell(period_k=2))
    assert t2 == pytest.approx(2.0 * t1, rel=1e-12)
```

Run: `uv run pytest tests/search/test_optimize_multibody.py -q` → **red**
(`_target_period_sec` returns ~2.4 yr from the V-E pair).

### Task 2.2 — minimal impl: dispatch on body count

Edit `optimize.py`. Add the import at the top of the resonance import block
(currently `from cyclerfinder.search.resonance import synodic_period_days`):

```python
from cyclerfinder.search.resonance import (
    beat_period_days,
    multi_body_beat_days,
    synodic_period_days,
)
```

Replace the body of `_target_period_sec`:

```python
def _target_period_sec(cell: Cell) -> float:
    """Resolve the target heliocentric period for ``cell``, seconds.

    For a 2-body cell, ``period_k * T_syn(bodies[0], bodies[1])`` (the M5
    native case, preserved byte-for-byte). For ``len(bodies) >= 3`` (M8's
    VEM territory) the single-pair formula is wrong — the period is
    ``period_k`` multiples of the body set's *beat period*
    (``3*T_syn(E,M) ~ 4*T_syn(E,V) ~ 6.40 yr`` for ``["V","E","M"]``).
    Dispatch to :func:`~cyclerfinder.search.resonance.multi_body_beat_days`
    and :func:`~cyclerfinder.search.resonance.beat_period_days`.

    ``period_k`` here is the *beat multiple* (Reading A, M8 plan §2): the
    catalogue loader converts an anchor-pair ``period.k`` into a
    beat-multiple before building the cell.
    """
    n = len(cell.bodies)
    if n < 2:
        raise ValueError(
            f"cell.bodies must have at least 2 entries; got {cell.bodies!r}",
        )
    if n == 2:
        t_syn_days = synodic_period_days(cell.bodies[0], cell.bodies[1])
        return t_syn_days * cell.period_k * SECONDS_PER_DAY
    bodies = list(cell.bodies)
    tuples = multi_body_beat_days(bodies)
    if not tuples:
        raise ValueError(
            f"no integer beat commensurability found for bodies={bodies!r} "
            f"within the resonance.multi_body_beat_days default tolerance; "
            f"this body set has no natural cycler beat in the searched k range",
        )
    beat_days = beat_period_days(bodies, tuples[0])
    return beat_days * cell.period_k * SECONDS_PER_DAY
```

Run: `uv run pytest tests/search/test_optimize_multibody.py -q` → **green**.

### Task 2.3 — guard test: 2-body path unchanged

Add to `test_optimize_multibody.py`:

```python
def test_two_body_period_path_unchanged() -> None:
    """The 2-body fast path is byte-identical to the pre-M8 single-pair formula."""
    from cyclerfinder.search.resonance import synodic_period_days
    cell = Cell(
        bodies=("E", "M"),
        sequence=("E", "M", "E"),
        period_k=2,
        per_leg_revs=(0, 0),
        per_leg_branch=("single", "single"),
    )
    expected = synodic_period_days("E", "M") * 2 * SECONDS_PER_DAY
    assert _target_period_sec(cell) == expected  # exact, not approx
```

Run the **whole existing optimise suite** to confirm no regression:
`uv run pytest tests/search/test_optimize.py tests/test_catalogue_rediscovery.py -q`
→ green. Then `uv run ruff check . && uv run ruff format --check . && uv run mypy src tests`.

### Task 2.4 — commit

```
search/optimize: multi-body beat-period dispatch in _target_period_sec

len(bodies)>=3 now resolves the total period from resonance.beat_period_days
(VEM 3xE-M~4xE-V~6.40yr) instead of the wrong single-pair formula. 2-body
path preserved byte-for-byte. Closes the optimize.py:229 M8 extension point.
```

---

## 3. Tisserand feasibility for 3-body sequences (M-L-free; tests + doc only)

`tisserand_feasible` (`sequence.py:265`) already loops over **every**
consecutive sequence pair (`sequence.py:317` `for i in range(len(cell.sequence) - 1)`)
and returns `True` iff each pair is `linkable` somewhere in
`(0.5, vinf_cap]`. It is therefore structurally correct for VEM sequences
*today*. The one honest limitation, already documented in the function
(`sequence.py:308-309`: "**Coplanar (i=0) only** — inherits the M2 `linkable`
predicate's restriction"), is that Venus's 3.4° and Mars's 1.85° inclinations
are ignored. M8 does **not** fix that (it is M-3D). M8 only:

1. Pins that VEM sequences pass/fail the *coplanar* Tisserand gate as expected,
   so a future M-3D change is a reviewed diff against a known baseline.
2. Adds a module-level note at the VEM call site (§4) that the feasibility is
   coplanar-only and over-permissive for the inclined real geometry — a
   candidate surviving the coplanar Tisserand gate is *necessary but not
   sufficient* for real-ephemeris closure (the M-ED filter is the sufficiency
   test).

### Task 3.1 — failing test: VEM consecutive pairs are coplanar-linkable

Write `tests/search/test_sequence_multibody.py`:

```python
"""M8 Tisserand feasibility for >=3-body sequences (spec §13.3).

Honest scope: tisserand_feasible is COPLANAR (i=0) only — it inherits
tisserand.linkable's restriction (sequence.py:308). These tests pin the
coplanar baseline so the M-3D inclination lift is a reviewed diff, not a
silent behaviour change. They do NOT assert real-geometry feasibility.
"""
from __future__ import annotations

from cyclerfinder.search.sequence import Cell, tisserand_feasible


def _emeeve_cell() -> Cell:
    seq = ("E", "M", "E", "E", "V", "E")  # vem-emeeve-3syn sequence_canonical
    # NOTE: adjacency E-E at index 2 is the same-body loop leg; enumerate_cells
    # forbids it but a catalogue-derived cell can carry it. tisserand_feasible
    # must handle a same-body pair gracefully (linkable returns False / a body
    # is trivially "linkable to itself" depends on the predicate — see Task 3.2).
    n_legs = len(seq) - 1
    return Cell(
        bodies=("V", "E", "M"),
        sequence=seq,
        period_k=1,
        per_leg_revs=(0,) * n_legs,
        per_leg_branch=("single",) * n_legs,
    )


def _vem_simple_cell() -> Cell:
    seq = ("E", "V", "M", "E")  # no same-body adjacency
    n_legs = len(seq) - 1
    return Cell(
        bodies=("V", "E", "M"),
        sequence=seq,
        period_k=1,
        per_leg_revs=(0,) * n_legs,
        per_leg_branch=("single",) * n_legs,
    )


def test_vem_simple_sequence_coplanar_linkable_at_7kms() -> None:
    """E-V-M-E: each consecutive pair (E-V, V-M, M-E) is coplanar-linkable
    somewhere in (0.5, 7.0] km/s (spec §13.3)."""
    assert tisserand_feasible(_vem_simple_cell(), vinf_cap=7.0) is True
```

Run: `uv run pytest tests/search/test_sequence_multibody.py::test_vem_simple_sequence_coplanar_linkable_at_7kms -q`.
This may pass immediately (the production function already supports it) — that
is acceptable for a *characterisation* test. If it is red, investigate before
proceeding (it would indicate `linkable` cannot bracket a V-M or M-E pair, a
real finding to surface). Either way, **record the observed result in the
test as the pinned baseline.**

### Task 3.2 — same-body-pair behaviour at the EMEEVE loop leg

The EMEEVE catalogue sequence contains an `E-E` adjacency (`catalogue.yaml:1777`),
which `enumerate_cells` forbids (`sequence.py:226`) but a catalogue-derived
cell carries. Add a test that documents how `tisserand_feasible` treats it:

```python
def test_emeeve_same_body_loop_leg_behaviour() -> None:
    """Pin how tisserand_feasible treats the E-E loop-leg adjacency.

    A same-body 'leg' is a multi-rev phasing loop, not a transfer; the
    coplanar Tisserand predicate is not meaningful for it. This test PINS
    the current behaviour (whatever it is) so M-N's structural-inference
    step (which builds the EMEEVE cell with the loop leg as n_revs>=1) knows
    whether it must special-case same-body legs before calling the gate.
    """
    cell = _emeeve_cell()
    result = tisserand_feasible(cell, vinf_cap=7.0)
    # Record the observed baseline. If False because linkable(E,E,..) cannot
    # bracket (E to itself has no synodic contour intersection band), the
    # structural-inference step in a later M-N/M-ED task must skip same-body
    # pairs in the feasibility check. Document the chosen value here:
    assert result in (True, False)  # replace with the observed constant + a comment
```

> **Action when writing:** replace the `in (True, False)` placeholder with the
> *actual observed* boolean and a one-line comment explaining it. Do **not**
> leave a tautological assert in the committed test — that is a no-op gate.
> This is a characterisation finding, not a placeholder step.

Run, observe, pin, then:
`uv run ruff check . && uv run ruff format --check . && uv run mypy src tests`.

### Task 3.3 — commit

```
search/sequence: pin coplanar Tisserand baseline for VEM 3-body sequences

Characterisation tests for tisserand_feasible on E-V-M-E and the EMEEVE
loop-leg sequence. No production change: the gate is already N-agnostic and
honestly coplanar-only (sequence.py:308). Baselines the i=0 behaviour so the
M-3D inclination lift surfaces as a reviewed diff.
```

---

## 4. VEM rediscovery gate validated against sourced anchors only

This is the milestone's binding gate. **What is sourced for the two VEM rows**
(verified against `data/catalogue.yaml`):

| row | sourced anchors | NOT sourced (all `null`) |
|---|---|---|
| `jones-2017-vem-triple-family` (`:1654`) | `bodies=["V","E","M"]`, `period.pair=E-M`, `period.k=2`, `period.years=4.27`, sequence is a *placeholder* (`:1664` "canonical placeholder; actual sequences vary") | every `vinf_kms` (`:1679-1690`), every `orbit_elements` field (`:1692-1703`), `legs: []` |
| `vem-emeeve-3syn` (`:1767`) | `bodies=["E","M","E","E","V","E"]`, `sequence_canonical="E-M-E-E-V-E"`, `period.pair=E-M`, `period.k=3`, `period.years=6.41` (`:1782`) | every `vinf_kms` (`:1789-1794`), every `orbit_elements` field (`:1796-1806`), `legs: []` |

**Therefore the only golden-test-legal EXPECTED values for a VEM gate are
period and sequence.** A V∞ assertion would have to invent a number our own
code produced — forbidden by the golden-test discipline. The gate asserts the
*period resolver and cell construction*, and the optimiser-convergence case is
`xfail` pending M-ED.

### Task 4.1 — gate test: EMEEVE period anchor round-trips through the resolver

Write `tests/test_vem_rediscovery.py`:

```python
"""M8 VEM rediscovery gate (spec §8 M8 anchor, §9 line 160).

GOLDEN DISCIPLINE: both VEM catalogue rows have null V_inf and null
orbit_elements (catalogue.yaml:1677-1703, 1788-1806). The ONLY sourced
anchors are period.years and sequence_canonical. This gate asserts those and
NOTHING our own optimiser computes. The converged-geometry case is xfail
pending M-ED (real-ephemeris discovery) — see plan §4.3.
"""
from __future__ import annotations

import pytest
import yaml

from cyclerfinder.core.constants import DAYS_PER_JULIAN_YEAR, SECONDS_PER_DAY
from cyclerfinder.search.optimize import _target_period_sec
from cyclerfinder.search.sequence import Cell
from tests._catalogue_loader import CATALOGUE_PATH


def _row(entry_id: str) -> dict:
    raw = yaml.safe_load(CATALOGUE_PATH.read_text())
    for row in raw:
        if row["id"] == entry_id:
            return row
    raise AssertionError(f"catalogue row {entry_id!r} not found")


def test_emeeve_sourced_period_matches_beat_resolver() -> None:
    """The EMEEVE row's SOURCED period (6.41 yr, catalogue.yaml:1782) matches
    what _target_period_sec resolves for its VEM body set at the natural beat.

    EXPECTED side = the catalogue's published period.years (sourced from
    Jones 2017 via the row). The resolver output is the side under test.
    """
    row = _row("vem-emeeve-3syn")
    sourced_years = float(row["period"]["years"])  # 6.41, Jones 2017
    assert sourced_years == pytest.approx(6.41, abs=0.005)  # fixture sanity

    # Build the beat-multiple cell (period_k=1 = the natural beat; the row's
    # anchor-pair k=3 maps to beat-multiple 1 per plan §2 / §5).
    seq = tuple(row["sequence_canonical"].split("-"))  # E,M,E,E,V,E
    n_legs = len(seq) - 1
    cell = Cell(
        bodies=("V", "E", "M"),
        sequence=seq,
        period_k=1,
        per_leg_revs=(0,) * n_legs,
        per_leg_branch=("single",) * n_legs,
    )
    resolved_years = (
        _target_period_sec(cell) / SECONDS_PER_DAY / DAYS_PER_JULIAN_YEAR
    )
    # Sourced 6.41 yr vs the 3xE-M=4xE-V beat 6.406 yr: agree within the
    # inter-source rounding (Jones rounds to 6.4; our constants give 6.406).
    assert resolved_years == pytest.approx(sourced_years, abs=0.02)
```

Run after §2 is merged → **green** (§2 makes the resolver correct).

### Task 4.2 — xfail: full VEM ballistic convergence (hands to M-ED)

Add to the same file:

```python
@pytest.mark.xfail(
    reason="VEM strict ballistic closure is open research (spec §17 line 199); "
    "no sourced V_inf anchor exists to assert against and the circular-coplanar "
    "optimiser is not expected to converge a ballistic VEM cycler. Flipped by "
    "M-ED (real-ephemeris discovery). See roadmap M-N test-gate row.",
    strict=False,
)
@pytest.mark.slow
def test_emeeve_idealized_optimiser_converges_feasible() -> None:
    """ASPIRATIONAL: the idealized optimiser finds a feasible (constraints-
    satisfied) interior solution for the EMEEVE VEM cell. Expected to xfail in
    the circular-coplanar model; documents the M-ED handoff target.

    This asserts ONLY result.constraints_satisfied (a feasibility predicate,
    not a sourced number) — it never asserts a computed V_inf as golden.
    """
    from cyclerfinder.core.ephemeris import Ephemeris
    from cyclerfinder.search.optimize import optimise_cell_idealized

    # EMEEVE with loop leg requires multi-rev (M-L). Until M-L lands this also
    # cannot construct; the xfail covers both the M-L and the convergence gaps.
    seq = ("E", "M", "E", "E", "V", "E")
    n_legs = len(seq) - 1
    cell = Cell(
        bodies=("V", "E", "M"),
        sequence=seq,
        period_k=1,
        per_leg_revs=(0, 0, 1, 0, 0),  # the E-E loop leg is multi-rev
        per_leg_branch=("single", "single", "low", "single", "single"),
    )
    result = optimise_cell_idealized(
        cell, Ephemeris(model="circular"), vinf_cap=7.0, seed=0,
    )
    assert result.constraints_satisfied
```

Run: `uv run pytest tests/test_vem_rediscovery.py -q` → the xfail registers as
`xfailed` (or `xpassed` if M-L+convergence somehow already work, which would
itself be a reviewable finding). `strict=False` so an unexpected xpass does not
break CI but is visible in the report.

### Task 4.3 — commit

```
test: VEM rediscovery gate against sourced period anchor + M-ED xfail

Asserts the EMEEVE row's sourced 6.41-yr period round-trips through the new
beat resolver. The full ballistic-convergence case is an xfail handed to
M-ED — no computed V_inf is ever asserted as golden (both VEM rows have null
V_inf; only period/sequence are sourced).
```

---

## 5. Loader admission: VEM rows from `NOT_TWO_BODY` to a categorised multibody class

Today both VEM rows are `ExclusionReason.NOT_TWO_BODY` (`_catalogue_loader.py:174-176`)
and the census ratchet freezes that at `2`
(`test_catalogue_rediscovery.py:272`). M8 promotes them to a *categorised,
non-silent* admission so the coverage-audit invariant (task #55: "no entry
vanishes silently") holds and the VEM gate can be parametrised over them.

**Important.** The VEM rows are **not** `CONSTRUCTIBLE` in the v1 sense — they
have no sourced V∞ to rediscover, and `jones-2017-vem-triple-family`'s
sequence is an explicit placeholder (`catalogue.yaml:1664`). Admitting them as
`CONSTRUCTIBLE` would feed the V∞-asserting rediscovery suite
(`test_catalogue_rediscovery.py`) garbage targets. Instead add a distinct
reason `CONSTRUCTIBLE_MULTIBODY` that means "3-body, period+sequence sourced,
admitted to the *VEM* gate (§4) but not to the 2-body V∞ rediscovery gate."

### Task 5.1 — failing test: census expects the new multibody class

Edit `tests/test_catalogue_rediscovery.py`. Update `EXPECTED_COVERAGE`
(`:267`) in the SAME commit as the loader change (the ratchet rule):

```python
EXPECTED_COVERAGE: dict[ExclusionReason, int] = {
    ExclusionReason.MULTI_ENCOUNTER_SEQUENCE: 202,
    ExclusionReason.NON_HELIOCENTRIC: 6,
    ExclusionReason.MISSING_VINF: 5,
    ExclusionReason.CONSTRUCTIBLE: 2,
    ExclusionReason.CONSTRUCTIBLE_MULTIBODY: 2,  # M8: the two VEM rows
    ExclusionReason.MISSING_PERIOD: 2,
    # NOT_TWO_BODY now 0 — both VEM rows promoted to CONSTRUCTIBLE_MULTIBODY.
}
```

Run: `uv run pytest tests/test_catalogue_rediscovery.py -k census -q` → **red**
(`AttributeError: CONSTRUCTIBLE_MULTIBODY` / count mismatch).

### Task 5.2 — minimal impl: new reason + VEM classification path

Edit `tests/_catalogue_loader.py`. Add the enum member after `CONSTRUCTIBLE`:

```python
    CONSTRUCTIBLE_MULTIBODY = "constructible_multibody"
    """3-body (VEM) row with sourced period + sequence but null V_inf /
    orbit_elements. Admitted to the M8 VEM gate (tests/test_vem_rediscovery.py),
    which asserts only period/sequence — NOT to the 2-body V_inf rediscovery
    gauntlet (no sourced V_inf to rediscover). See M8 plan §5."""
```

In `classify_row`, replace the unconditional `NOT_TWO_BODY` branch
(`:174-176`) with a 3-body VEM sub-classification. Keep the order: heliocentric
and ballistic checks first (so a non-Sun or low-thrust VEM row would still fall
out earlier):

```python
    bodies = row.get("bodies") or []
    if len(bodies) != 2:
        # 3-body VEM rows with a sourced period are admitted to the M8 VEM
        # gate as a categorised (non-silent) class; everything else stays
        # NOT_TWO_BODY. Period must be present (otherwise MISSING_PERIOD-like
        # rows would masquerade as constructible-multibody).
        period = row.get("period") or {}
        bodyset = set(bodies)
        is_vem = len(bodies) >= 3 and bodyset <= {"V", "E", "M"}
        if is_vem and period.get("years") is not None and period.get("k") is not None:
            entry = CatalogueEntry(
                id=row["id"],
                name=row.get("name", row["id"]),
                bodies=tuple(bodies),
                sequence_canonical=row.get("sequence_canonical") or "",
                period_k=int(period["k"]),
                period_years=float(period["years"]),
                vinf_targets_kms=(),   # sourced: none — VEM rows have null V_inf
                leg_tofs_days=(),      # sourced: none — VEM rows have legs: []
            )
            return ExclusionReason.CONSTRUCTIBLE_MULTIBODY, entry
        return ExclusionReason.NOT_TWO_BODY, None
```

> **Note on the anchor-pair→beat-multiple mapping (plan §2 design decision):**
> the `CatalogueEntry.period_k` stored here is the *catalogue's anchor-pair k*
> (3 for EMEEVE), NOT the beat-multiple. The VEM gate (§4) and any future VEM
> cell builder must convert: anchor-pair `k=3` on `(E,M)` ⇒ beat-multiple `1`
> (because the natural beat *is* 3×E-M). Document this at both call sites. Do
> not silently reinterpret `period_k`.

Run: `uv run pytest tests/test_catalogue_rediscovery.py -k census -q` → **green**.

### Task 5.3 — guard: VEM rows do NOT leak into the V∞ rediscovery gauntlet

`load_constructible_entries` (`_catalogue_loader.py:224`) returns only
`CONSTRUCTIBLE` (not `CONSTRUCTIBLE_MULTIBODY`), so the V∞-asserting
`test_catalogue_entry_rediscovers` parametrisation is unchanged. Pin it:

```python
def test_vem_rows_excluded_from_vinf_gauntlet() -> None:
    """The VEM rows are CONSTRUCTIBLE_MULTIBODY, not CONSTRUCTIBLE: they must
    NOT appear in the 2-body V_inf rediscovery parametrisation (they have no
    sourced V_inf). Guards against a future loosening that would feed the
    gauntlet null/garbage V_inf targets."""
    from tests._catalogue_loader import load_constructible_entries
    ids = {e.id for e in load_constructible_entries()}
    assert "vem-emeeve-3syn" not in ids
    assert "jones-2017-vem-triple-family" not in ids
```

Run the full rediscovery suite (not slow) + lint/type:
`uv run pytest tests/test_catalogue_rediscovery.py -m "not slow" -q`
then `uv run ruff check . && uv run ruff format --check . && uv run mypy src tests`.

### Task 5.4 — wire the VEM gate to the new loader class (optional, same commit)

Optionally add a parametrised version of the §4 gate that pulls the
`CONSTRUCTIBLE_MULTIBODY` entries from the loader and asserts each row's
sourced period round-trips the resolver (generalises Task 4.1 across both VEM
rows). Skip `jones-2017-vem-triple-family` for the *sequence* leg of the
assertion (its sequence is a placeholder, `catalogue.yaml:1664`) — assert only
its period for that row, with a documented `pytest.mark` reason.

### Task 5.5 — commit

```
tests/_catalogue_loader: admit VEM rows as CONSTRUCTIBLE_MULTIBODY

Both VEM rows move out of NOT_TWO_BODY into a categorised multibody class
carrying their sourced period/sequence (V_inf stays empty — the rows have
none). The 2-body V_inf rediscovery gauntlet is unchanged (guarded). Census
ratchet updated in the same commit.
```

---

## 6. Out of scope (explicit carve-outs)

These are real M8-adjacent work items deliberately **not** in this plan, each
with where it actually lands. Do not stub them here.

| Carve-out | Where it lands | Rationale |
|---|---|---|
| **Planet-centric / non-heliocentric moon tours** (the 6 `NON_HELIOCENTRIC` rows, task **#76** — lunar / Jovian tours) | **Separate phase** (call it M9-moontour or keep under task #76). | Orthogonal physics: planet-centric `primary != "Sun"`, an incommensurable signature, a different ephemeris frame, and a different `linkable`/Tisserand basis. The roadmap (scope table line 19) explicitly says "Non-heliocentric moon tours stay under separate task #76 — planet-centric, orthogonal." Tangling it with VEM would conflate two unrelated dynamical regimes. **No code in this plan touches the `NON_HELIOCENTRIC` rows** (the census keeps them at 6). |
| **Multi-rev Lambert solver** (the `lambert(max_revs=N)` math) | **M-L** (roadmap, "ready for task-level planning: YES"). | Foundational and independent. VEM *closure-sequence* convergence needs it (the EMEEVE loop leg), which is why this plan keeps that case as `xfail`. The period/feasibility/loader work here is M-L-free. |
| **Full-3D inclination** (Venus 3.4°, Mars 1.85°) | **M-3D** (roadmap, "needs design first"). | The coplanar Tisserand limitation (§3) is honestly documented, not fixed. Real VEM closure almost certainly needs 3D (Venus's plane change is where the cheap steering comes from, spec §12.1 step 3). |
| **Real-ephemeris VEM discovery / convergence** (`optimise_cell_ephemeris`, TCM budget, phase-match) | **M-ED** (roadmap, "needs design first"). | This is where a VEM candidate is actually found and reported. The §4.2 `xfail` is its handoff target. |
| **CLI `--bodies V,E,M --period beat` + viz + reporting** (spec §8 line 152, §6 line 207) | **M8-UX sub-stream** (this plan is the M8 *search-core* sub-stream). | The spec bundles "VEM campaign + UX" under one milestone label, but the UX (cli/viz/docs) is a distinct deliverable with its own gate and is best planned separately once the search core lands. Flagged so the human can decide whether to split the milestone label. |
| **Structural inference** (catalogue `E-M-E-E-V-E @ k` → `Cell` with correct per-leg `n_revs`/`branch`) | **M-N/M-ED** (roadmap line 94). | Needs M-L to construct the loop leg. This plan hand-builds the cell in tests; the automated inference belongs after M-L. |

---

## 7. Definition of done

1. `tests/search/test_optimize_multibody.py` — VEM period resolves to the
   6.406-yr beat; scales with `period_k`; 2-body path byte-identical. **Green.**
2. `tests/search/test_sequence_multibody.py` — VEM coplanar Tisserand baseline
   pinned (with the observed same-body-pair constant documented, not a
   tautology). **Green.**
3. `tests/test_vem_rediscovery.py` — sourced EMEEVE period round-trips the
   resolver (green); full ballistic convergence registers `xfail` (handed to
   M-ED).
4. `tests/_catalogue_loader.py` — both VEM rows classify as
   `CONSTRUCTIBLE_MULTIBODY`; `NOT_TWO_BODY` count is 0; VEM rows excluded from
   the V∞ gauntlet (guarded). Census ratchet updated in the same commit.
5. No assertion anywhere uses a V∞ or `(a,e)` value our own code computed as
   the EXPECTED side — only sourced period/sequence and feasibility predicates.
6. `uv run pytest -m "not slow"`, `uv run ruff check .`,
   `uv run ruff format --check .`, `uv run mypy src tests` all clean. CI green.
7. The roadmap doc's M-N "test gates" row is satisfied: the `NOT_TWO_BODY`
   exclusion count drops, and at least one VEM triple has a rediscovery case
   (period-anchored green + convergence xfail).

---

## 8. Open design questions for the human (resolve before execution)

1. **Period convention (§2 Reading A vs B).** Confirm `_target_period_sec`
   interprets `Cell.period_k` as a *beat multiple* and the loader converts the
   catalogue's anchor-pair `period.k` (3 for EMEEVE) into beat-multiple 1.
   Alternative: add a `period_basis`/`period_pair` field to the frozen `Cell`
   dataclass + `Cell.id` format — a wider M4 change flagged but not taken
   unilaterally.
2. **`CONSTRUCTIBLE_MULTIBODY` vs a wider loader refactor.** Is a new
   `ExclusionReason` the right shape, or should the loader grow a separate
   `load_multibody_entries()` surface? The plan takes the minimal-diff path
   (new reason, V∞ gauntlet untouched); confirm that is preferred over a
   loader API change.
3. **Milestone label split.** Spec §8 bundles "VEM campaign + UX" as one M8.
   This plan delivers only the *search core*. Confirm whether CLI/viz/docs
   should be a separate phase (M8-UX) or folded back in here.
4. **EMEEVE vs the §13.8 worked sequence.** The catalogue's sourced EMEEVE
   sequence is `E-M-E-E-V-E` (`catalogue.yaml:1777`); the spec §13.8 worked
   cell-id example is `E-V-M-E-M-E`. They are different itineraries. Confirm
   the gate should anchor on the *catalogue's sourced* sequence (it must, for
   golden discipline) and the §13.8 string is illustrative only.
5. **Beat-tuple selection robustness.** `_target_period_sec` takes
   `multi_body_beat_days(bodies)[0]` (the lowest-mismatch tuple). For VEM at
   the default `k_max=6`/`tol_frac=0.02` this is unambiguously `(4,3)`
   (resonance gate `test_resonance.py:54`). Confirm we do not need to thread a
   non-default `k_max`/`tol_frac` for higher-period VEM branches (12.8 yr / 32
   yr, spec §3 line 41) in this milestone — the plan defers those to M-ED.
