# M8-Core — Multi-body (N≥3) VEM cyclers: beat-period dispatch, N≥3 enumeration, rediscovery gate (search-core only)

> Written with the `superpowers:writing-plans` conventions. This is a
> **task-level TDD plan** (write failing test → run red → minimal impl → run
> green → commit). It is the **M8-Core** (search-core) detailed expansion of
> the **M-N** milestone in
> `docs/superpowers/plans/2026-06-02-multirev-3d-vem-ephemeris-roadmap.md`
> (which calls it "N-encounter + VEM enumeration & rediscovery, deps: M-L").
> The user-facing CLI/viz/docs are carved into a separate **M8-UX** milestone
> (see §6). Read that roadmap's "Verified current state" section first — it is
> the authoritative survey of what is already built.

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
   sequence pair, so it works structurally for VEM sequences today; M8 adds a
   same-body-pair bypass (so the EMEEVE E-E loop leg is not falsely
   evaluated by the meaningless `linkable(X,X,..)` predicate), plus tests
   pinning that behaviour and documenting the coplanar caveat at the VEM-cell
   boundary. **No lift of the coplanar assumption here** — that is M-3D,
   deliberately out of M8 scope (see §6).
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
  reporting, docs." This plan delivers the **M8-Core** search half; `cli`/`viz`/
  reporting/docs are the separate **M8-UX** milestone (see §6 scope note).
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

M8-Core is **additive** to three source files and one test-support file, plus
new tests:

- `src/cyclerfinder/search/sequence.py` — `Cell` gains an optional
  `period_basis` field + `id` token (§2 Task 2.0), and `tisserand_feasible`
  gains a same-body-pair bypass (§3 Task 3.2). Existing 2-body cell ids and
  distinct-pair feasibility behaviour are unchanged.
- `src/cyclerfinder/search/optimize.py` — `_target_period_sec` gains a
  `period_basis`/multi-body dispatch (§2). No signature change; the 2-body
  basis-None path is byte-for-byte preserved.
- `tests/_catalogue_loader.py` — gains an N-agnostic multibody classification
  path and a new `ExclusionReason.CONSTRUCTIBLE_MULTIBODY` (§5).
- New tests: `tests/search/test_sequence_cell_basis.py` (§2 Task 2.0),
  `tests/search/test_optimize_multibody.py` (§2, §4),
  `tests/search/test_sequence_multibody.py` (§3),
  `tests/test_vem_rediscovery.py` (§4 gate, §5).
- `tests/test_catalogue_rediscovery.py::EXPECTED_COVERAGE` — census ratchet
  update (§5).

No new production module. No edits to `resonance.py` or `construct.py`. The
`sequence.py` change in the M-L-free portion is limited to the small
same-body-pair bypass in `tisserand_feasible` (§3 Task 3.2) plus the
`Cell.period_basis` field (§2 Task 2.0); the rest of §3 pins existing
behaviour and documents the coplanar caveat.

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
sourced period is 6.41 yr (`catalogue.yaml:1782`), the 3×E-M = 4×E-V beat. The
fix carries the catalogue's anchor pair on the Cell (`period_basis`, see the
design decision below) so the resolver can reconstruct `3 * T_syn(E,M) ≈
6.41 yr` without rewriting `period_k`.

### Design decision — traceable period resolution (settled)

**Do not silently mutate `period.k`.** The earlier draft proposed having the
loader convert the catalogue's anchor-pair `period.k=3` into a beat-multiple
`Cell.period_k=1`. That is an anti-pattern: a silent semantic shift across the
loader boundary destroys traceability. If the catalogue defines `k=3` relative
to E-M and the Cell drops it to `k=1`, a Cell printed from `vem-emeeve-3syn`
showing `period_k=1` looks like a bug versus the YAML — the Cell no longer
represents what the catalogue/user requested.

**Decision: carry the basis on the Cell instead.** Add an optional
`period_basis: tuple[str, str] | None` field to the frozen `Cell` dataclass.
The Cell then represents *exactly* what the catalogue/user requested: the
anchor pair plus the (unmodified) `period_k`.

- The loader (§5) maps the catalogue row's `period.pair` (e.g. `"E-M"`) into
  `Cell.period_basis = ("E", "M")` and keeps `period_k` as the sourced value
  (3 for EMEEVE).
- `_target_period_sec` computes the total period by dispatch:
  - **`period_basis` set:** `synodic_period_days(*cell.period_basis) *
    cell.period_k * SECONDS_PER_DAY`. For EMEEVE this is `T_syn(E,M) * 3 ≈
    2.135 * 3 = 6.405 yr` — the anchor-pair reading, with `period_k`
    untouched and fully traceable to the YAML.
  - **`period_basis is None and len(bodies) >= 3`:** fall back to the natural
    beat via `multi_body_beat_days(bodies)[0]` fed to `beat_period_days`
    (the resonance-pure default for cells with no declared basis).
  - **`period_basis is None and len(bodies) == 2`:** the existing 2-body path
    (`synodic_period_days(bodies[0], bodies[1]) * period_k`) — byte-identical
    to pre-M8.

Both the anchor-pair path and the beat-fallback agree numerically at the
natural beat (`3*T_syn(E,M) ≈ 4*T_syn(E,V) ≈ 6.40 yr`); that is the whole
point of a beat. The difference is *traceability*: with `period_basis` the Cell
faithfully echoes the catalogue, and `_target_period_sec` derives the right
number without anyone rewriting `period_k`.

**`Cell` hashability / id stability.** `Cell` is frozen and hashable and its
`id` (`sequence.py:105`) is used widely (cell ledgers, dedupe, logs). Update
`Cell.id` to incorporate `period_basis` *only when set* — e.g. append a
`|p<AB>` token (`|pEM` for `("E","M")`). When `period_basis is None` the id is
unchanged, so every existing 2-body cell id stays stable (no ledger churn,
no dedupe breakage). Task 2.0 below adds the dataclass field + id change and a
dedicated test; Task 2.2 adds the `_target_period_sec` dispatch.

### Task 2.0 — `Cell.period_basis` field + id token (M-L-free, do first)

The dataclass change underpins the rest of §2. `Cell` is frozen
(`sequence.py:57`) with fields `bodies, sequence, period_k, per_leg_revs,
per_leg_branch` (`sequence.py:99-103`) and an `id` property
(`sequence.py:105-125`).

#### Failing test — `tests/search/test_sequence_cell_basis.py`

```python
"""M8: Cell.period_basis carries the catalogue's anchor pair without
mutating period_k (plan §2 design decision)."""
from __future__ import annotations

from cyclerfinder.search.sequence import Cell


def _emeeve_basis_cell() -> Cell:
    seq = ("E", "M", "E", "E", "V", "E")
    n_legs = len(seq) - 1
    return Cell(
        bodies=("V", "E", "M"),
        sequence=seq,
        period_k=3,                  # sourced anchor-pair k, NOT rewritten
        per_leg_revs=(0,) * n_legs,
        per_leg_branch=("single",) * n_legs,
        period_basis=("E", "M"),
    )


def test_period_basis_defaults_to_none() -> None:
    """A 2-body cell built the old way has period_basis is None."""
    cell = Cell(
        bodies=("E", "M"),
        sequence=("E", "M", "E"),
        period_k=2,
        per_leg_revs=(0, 0),
        per_leg_branch=("single", "single"),
    )
    assert cell.period_basis is None


def test_period_basis_preserved_and_period_k_untouched() -> None:
    cell = _emeeve_basis_cell()
    assert cell.period_basis == ("E", "M")
    assert cell.period_k == 3  # traceable to catalogue.yaml:1782


def test_id_unchanged_when_basis_none() -> None:
    """Existing 2-body cell ids stay byte-identical (no ledger churn)."""
    cell = Cell(
        bodies=("E", "M"),
        sequence=("E", "M", "E"),
        period_k=2,
        per_leg_revs=(0, 0),
        per_leg_branch=("single", "single"),
    )
    assert cell.id == "EM|E-M-E|k2|r00|bss"


def test_id_appends_basis_token_when_set() -> None:
    """When period_basis is set, the id gains a |p<AB> token so a basis-bearing
    cell is distinguishable (and the YAML-traceable k3 stays visible)."""
    cell = _emeeve_basis_cell()
    assert cell.id == "VEM|E-M-E-E-V-E|k3|r00000|bsssss|pEM"
```

Run: `uv run pytest tests/search/test_sequence_cell_basis.py -q` → **red**
(`Cell` has no `period_basis`).

#### Minimal impl — `sequence.py`

Add the field (keyword-only-by-default-value so all existing positional
construction sites keep working unchanged):

```python
    bodies: tuple[str, ...]
    sequence: tuple[str, ...]
    period_k: int
    per_leg_revs: tuple[int, ...]
    per_leg_branch: tuple[str, ...]
    period_basis: tuple[str, str] | None = None
```

Append the token to `id` only when set:

```python
        base = f"{bodyset}|{sequence}|k{self.period_k}|r{revs}|b{branches}"
        if self.period_basis is not None:
            base += f"|p{''.join(self.period_basis)}"
        return base
```

Run: `uv run pytest tests/search/test_sequence_cell_basis.py -q` → **green**,
then the full sequence suite to confirm no id regressions:
`uv run pytest tests/search/test_sequence.py -q`, then lint/type.

Commit:

```
search/sequence: add Cell.period_basis (anchor pair) without mutating period_k

Frozen Cell gains an optional period_basis tuple so a catalogue-derived cell
echoes its anchor pair and sourced period_k verbatim (traceability). Cell.id
appends a |p<AB> token only when set; ids for existing basis-None cells are
byte-identical. Underpins the _target_period_sec dispatch (plan §2).
```

### Task 2.1 — failing test: VEM 3-body period resolves via anchor-pair basis

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


def _vem_cell(period_k: int = 3, basis: tuple[str, str] | None = ("E", "M")) -> Cell:
    # Anchor-pair convention (plan §2): the catalogue's sourced EMEEVE sequence
    # E-M-E-E-V-E (catalogue.yaml:1777) with its sourced anchor pair (E,M) and
    # sourced k=3. period_k is NOT rewritten — the basis tells the resolver how
    # to interpret it. The spec §13.8 E-V-M-E-M-E string is illustrative only.
    seq = ("E", "M", "E", "E", "V", "E")
    n_legs = len(seq) - 1
    return Cell(
        bodies=("V", "E", "M"),
        sequence=seq,
        period_k=period_k,
        per_leg_revs=(0,) * n_legs,
        per_leg_branch=("single",) * n_legs,
        period_basis=basis,
    )


def test_vem_anchor_pair_period_is_6_41yr() -> None:
    """EMEEVE cell with basis (E,M) and k=3 resolves to 3*T_syn(E,M) ~ 6.41 yr,
    the sourced beat (spec §9, catalogue.yaml:1782)."""
    t_sec = _target_period_sec(_vem_cell(period_k=3, basis=("E", "M")))
    t_yr = t_sec / SECONDS_PER_DAY / DAYS_PER_JULIAN_YEAR
    assert t_yr == pytest.approx(6.41, abs=0.02), f"got {t_yr:.4f} yr"


def test_vem_no_basis_falls_back_to_natural_beat() -> None:
    """With period_basis=None and >=3 bodies, the resolver falls back to the
    natural beat multi_body_beat_days(...)[0] ~ 6.406 yr (spec §9)."""
    t_sec = _target_period_sec(_vem_cell(period_k=1, basis=None))
    t_yr = t_sec / SECONDS_PER_DAY / DAYS_PER_JULIAN_YEAR
    assert t_yr == pytest.approx(6.406, abs=0.01), f"got {t_yr:.4f} yr"


def test_vem_anchor_period_scales_with_period_k() -> None:
    """period_k scales the anchor-pair period linearly (the higher branches,
    spec §3 line 41); doubling k doubles T."""
    t1 = _target_period_sec(_vem_cell(period_k=3, basis=("E", "M")))
    t2 = _target_period_sec(_vem_cell(period_k=6, basis=("E", "M")))
    assert t2 == pytest.approx(2.0 * t1, rel=1e-12)
```

Run: `uv run pytest tests/search/test_optimize_multibody.py -q` → **red**
(`_target_period_sec` returns ~2.4 yr from the V-E pair / ignores the basis).

### Task 2.2 — minimal impl: dispatch on period_basis, then body count

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

    Dispatch (M8 plan §2):

    * ``cell.period_basis`` set — use the catalogue's *anchor pair*:
      ``T_syn(*period_basis) * period_k``. ``period_k`` is the sourced
      catalogue value, never rewritten, so a Cell stays traceable to its
      YAML row (EMEEVE: ``T_syn(E,M) * 3 ~ 6.41 yr``).
    * ``period_basis is None`` and ``len(bodies) >= 3`` — fall back to the
      body set's *natural beat* via
      :func:`~cyclerfinder.search.resonance.multi_body_beat_days` /
      :func:`~cyclerfinder.search.resonance.beat_period_days`
      (``3*T_syn(E,M) ~ 4*T_syn(E,V) ~ 6.40 yr`` for ``["V","E","M"]``).
    * ``period_basis is None`` and ``len(bodies) == 2`` — the M5 native
      single-pair formula, preserved byte-for-byte.
    """
    if cell.period_basis is not None:
        a, b = cell.period_basis
        t_syn_days = synodic_period_days(a, b)
        return t_syn_days * cell.period_k * SECONDS_PER_DAY

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
search/optimize: period_basis + beat dispatch in _target_period_sec

period_basis-bearing cells resolve via the catalogue anchor pair
(T_syn(*basis)*period_k, EMEEVE 3xE-M~6.41yr) with period_k untouched.
Basis-None >=3-body cells fall back to resonance.beat_period_days; the
2-body path is preserved byte-for-byte. Closes the optimize.py:229 M8
extension point without mutating period_k.
```

---

## 3. Tisserand feasibility for 3-body sequences (M-L-free; same-body bypass + tests)

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
3. Adds a **same-body-pair bypass** to `tisserand_feasible` (Task 3.2). A
   same-body "leg" (e.g. `E->E` in the EMEEVE loop) is a multi-rev
   phasing/resonant loop, *not* a transfer, so the coplanar `linkable()`
   predicate cannot meaningfully evaluate it. This is a small, justified
   production change (a `continue` skip), not the blind-pinning the earlier
   draft proposed.

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
    # NOTE: adjacency E-E at index 3 is the same-body loop leg; enumerate_cells
    # forbids it but a catalogue-derived cell can carry it. tisserand_feasible
    # bypasses same-body pairs (Task 3.2) because linkable(X,X,..) is trivially
    # and meaninglessly True; loop legs are validated by M-L/M-ED, not Tisserand.
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

### Task 3.2 — bypass same-body pairs in `tisserand_feasible`

The EMEEVE catalogue sequence contains an `E-E` adjacency (`catalogue.yaml:1777`),
which `enumerate_cells` forbids (`sequence.py:226`) but a catalogue-derived
cell carries. The earlier draft proposed merely *pinning* whatever
`linkable(E,E,..)` returns. That is wrong — a same-body "leg" is a multi-rev
phasing/resonant loop, not a transfer, and a blind-pinned coplanar result is
either a meaningless False (which would falsely reject the whole EMEEVE cell)
or a meaningless trivial True. We resolve this by *skipping* same-body pairs.

#### Step 1 — mathematically verify how `linkable` handles identical bodies (finding)

Read `linkable` (`tisserand.py:302-440`). For `body_a == body_b`:

- `vinf_to_tisserand(body_a, vinf) == vinf_to_tisserand(body_b, vinf)`, so
  `t_pa == t_pb` and `_a_branches_at_e(...)` produces **identical** branch
  lists: `branches_a[i] == branches_b[i]` at every `e` sample
  (`tisserand.py:359-364`).
- The first inner-loop check (`tisserand.py:372-375`) computes
  `g_lo = a_branch_lo - b_branch_lo`. For the *same* branch value at the same
  sample this is exactly `0.0`, so `abs(g_lo) <= tol_au` is **True** and
  `linkable` returns **True immediately** at `i = 0`.

**Finding:** `linkable(X, X, vinf)` is *trivially* True (a contour is identical
to itself) for any vinf above the floor — there is no divide-by-zero (the
`brentq` refinement is never reached), but the answer is **physically
meaningless**: it asserts "X is reachable from X at fixed (a,e)", which says
nothing about a phasing loop. So pinning it would bake a meaningless True into
the baseline. Document this finding in the test docstring.

#### Step 2 — production change: skip same-body pairs

Edit `tisserand_feasible` (`sequence.py:317-326`). Inside the
`for i in range(len(cell.sequence) - 1)` loop, before computing `body_a`/
`body_b`'s linkability:

```python
        for i in range(len(cell.sequence) - 1):
            body_a = cell.sequence[i]
            body_b = cell.sequence[i + 1]
            if body_a == body_b:
                # Same-body adjacency (e.g. the EMEEVE E-E loop leg) is a
                # multi-rev phasing/resonant loop, not a transfer. The coplanar
                # linkable() predicate is meaningless here (it returns a trivial
                # True for X-to-X). Skip it; loop legs are validated by the
                # multi-rev Lambert/closure path (M-L/M-ED), not by Tisserand.
                continue
            found = False
            ...
```

The distinct-pair behaviour is unchanged byte-for-byte.

#### Step 3 — test

```python
def test_emeeve_loop_leg_bypass_returns_true() -> None:
    """An EMEEVE-style cell with an E-E loop leg returns True under the
    same-body bypass — the loop leg does NOT falsely trigger rejection.

    Finding (plan §3 Task 3.2 step 1): linkable(E, E, vinf) is trivially True
    (a contour equals itself) and physically meaningless, so we skip same-body
    pairs rather than consult linkable. Every distinct pair (E-M, M-E, E-V,
    V-E) is coplanar-linkable in (0.5, 7.0], so the cell is feasible.
    """
    assert tisserand_feasible(_emeeve_cell(), vinf_cap=7.0) is True


def test_distinct_pair_feasibility_unchanged() -> None:
    """The bypass does not alter distinct-pair behaviour: a plain E-V-M-E cell
    still passes exactly as before."""
    assert tisserand_feasible(_vem_simple_cell(), vinf_cap=7.0) is True
```

Run: `uv run pytest tests/search/test_sequence_multibody.py -q` → green, then
the full sequence/enumeration suite to confirm distinct-pair regressions are
clean: `uv run pytest tests/search/test_sequence.py -q`, then
`uv run ruff check . && uv run ruff format --check . && uv run mypy src tests`.

### Task 3.3 — commit

```
search/sequence: bypass same-body pairs in tisserand_feasible

Same-body adjacencies (the EMEEVE E-E loop leg) are multi-rev phasing loops,
not transfers; linkable(X,X,..) is trivially and meaninglessly True, so the
gate now skips them (loop legs are validated by M-L/M-ED). Distinct-pair
behaviour is byte-identical. Tests pin the EMEEVE bypass and the coplanar
E-V-M-E baseline so the M-3D inclination lift surfaces as a reviewed diff.
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
pending M-ED (real-ephemeris discovery) — see plan §4 Task 4.2.
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


def test_emeeve_sourced_period_matches_anchor_resolver() -> None:
    """The EMEEVE row's SOURCED period (6.41 yr, catalogue.yaml:1782) matches
    what _target_period_sec resolves from the row's anchor pair + sourced k.

    EXPECTED side = the catalogue's published period.years (sourced from
    Jones 2017 via the row). The resolver output is the side under test.
    """
    row = _row("vem-emeeve-3syn")
    sourced_years = float(row["period"]["years"])  # 6.41, Jones 2017
    assert sourced_years == pytest.approx(6.41, abs=0.005)  # fixture sanity

    # Build the cell exactly as the loader (§5) would: anchor pair from
    # period.pair, sourced period_k UNCHANGED (no silent k rewrite, plan §2).
    seq = tuple(row["sequence_canonical"].split("-"))  # E,M,E,E,V,E
    basis = tuple(row["period"]["pair"].split("-"))    # ("E","M")
    n_legs = len(seq) - 1
    cell = Cell(
        bodies=("V", "E", "M"),
        sequence=seq,
        period_k=int(row["period"]["k"]),  # 3, traceable to the YAML
        per_leg_revs=(0,) * n_legs,
        per_leg_branch=("single",) * n_legs,
        period_basis=basis,
    )
    resolved_years = (
        _target_period_sec(cell) / SECONDS_PER_DAY / DAYS_PER_JULIAN_YEAR
    )
    # Sourced 6.41 yr vs 3*T_syn(E,M) ~ 6.405 yr: agree within inter-source
    # rounding (Jones rounds to 6.4; our constants give 6.405).
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
    # Built as the loader would: anchor pair (E,M), sourced k=3, no k rewrite.
    seq = ("E", "M", "E", "E", "V", "E")
    n_legs = len(seq) - 1
    cell = Cell(
        bodies=("V", "E", "M"),
        sequence=seq,
        period_k=3,
        per_leg_revs=(0, 0, 1, 0, 0),  # the E-E loop leg is multi-rev
        per_leg_branch=("single", "single", "low", "single", "single"),
        period_basis=("E", "M"),
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

Asserts the EMEEVE row's sourced 6.41-yr period round-trips through the
anchor-pair resolver (period_basis=(E,M), sourced k=3, no k rewrite). The
full ballistic-convergence case is an xfail handed to M-ED — no computed
V_inf is ever asserted as golden (both VEM rows have null V_inf; only
period/sequence are sourced).
```

---

## 5. Loader admission: ≥3-body rows from `NOT_TWO_BODY` to a categorised multibody class

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
reason `CONSTRUCTIBLE_MULTIBODY` that means "≥3-body, period+sequence sourced,
admitted to the *VEM* gate (§4) but not to the 2-body V∞ rediscovery gate."

**N-agnostic.** The classification does **not** hardcode a `{"V","E","M"}`
bodyset — it admits any `len(bodies) >= 3` row with a valid period block. The
two current VEM rows are the only ≥3-body rows in the catalogue today, but a
future Jovian or Saturnian ≥3-body row would be admitted to the same
categorised class without a loader edit (keeping `linkable`/Tisserand
heliocentric assumptions aside — those remain a separate gate). This avoids a
VEM-specific special case that would have to be revisited for every new body
set.

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

### Task 5.2 — minimal impl: new reason + N-agnostic multibody classification

Edit `tests/_catalogue_loader.py`. Add the enum member after `CONSTRUCTIBLE`:

```python
    CONSTRUCTIBLE_MULTIBODY = "constructible_multibody"
    """>=3-body row with sourced period + sequence but null V_inf /
    orbit_elements. Admitted to the M8 VEM gate (tests/test_vem_rediscovery.py),
    which asserts only period/sequence — NOT to the 2-body V_inf rediscovery
    gauntlet (no sourced V_inf to rediscover). See M8-Core plan §5."""
```

In `classify_row`, replace the unconditional `NOT_TWO_BODY` branch
(`:174-176`) with an **N-agnostic** sub-classification — do **not** hardcode a
VEM bodyset. Any `len(bodies) >= 3` row with a valid period block (both
`years` and `k` present) is admitted; everything else stays `NOT_TWO_BODY`.
Keep the order: heliocentric and ballistic checks first (so a non-Sun or
low-thrust ≥3-body row would still fall out earlier):

```python
    bodies = row.get("bodies") or []
    if len(bodies) != 2:
        period = row.get("period") or {}
        if len(bodies) >= 3 and period.get("years") is not None and period.get("k") is not None:
            entry = CatalogueEntry(
                id=row["id"],
                name=row.get("name", row["id"]),
                bodies=tuple(bodies),
                sequence_canonical=row.get("sequence_canonical") or "",
                period_k=int(period["k"]),
                period_years=float(period["years"]),
                vinf_targets_kms=(),
                leg_tofs_days=(),
            )
            return ExclusionReason.CONSTRUCTIBLE_MULTIBODY, entry
        return ExclusionReason.NOT_TWO_BODY, None
```

> **Note on traceable period (plan §2 design decision):** the
> `CatalogueEntry.period_k` stored here is the *catalogue's sourced k* (3 for
> EMEEVE) on the row's anchor pair (`period.pair`), and is **never rewritten**.
> The VEM gate (§4) and any future VEM cell builder construct the `Cell` with
> `period_basis = period.pair` and `period_k` verbatim, and `_target_period_sec`
> derives `T_syn(*period_basis) * period_k` from those. There is no
> "beat-multiple" reinterpretation of `period_k` anywhere — that was rejected
> (plan §2) precisely because it destroyed traceability.

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
| **CLI `--bodies V,E,M --period beat` + viz + reporting + docs** (spec §8 line 152, §6 line 207) | **M8-UX** (a separate milestone; see the subsection below). | Decided: split out. The spec bundles "VEM campaign + UX" under one M8 label, but the UX is a distinct deliverable with its own gate. |
| **Structural inference** (catalogue `E-M-E-E-V-E @ k` → `Cell` with correct per-leg `n_revs`/`branch`) | **M-N/M-ED** (roadmap line 94). | Needs M-L to construct the loop leg. This plan hand-builds the cell in tests; the automated inference belongs after M-L. |

### Out of scope → M8-UX (separate milestone)

This plan is **M8-Core** (search-core only): the routing, the Tisserand
feasibility gate, the VEM rediscovery gate, and the loader admission. The
user-facing surface is carved into a separate **M8-UX** milestone:

- **CLI** — `--bodies V,E,M --period beat` (and anchor-pair flags) wiring the
  enumerator/optimiser to a command-line entry point.
- **Visualisation** — porkchop / trajectory / beat-diagram plots for VEM cells.
- **Reporting + docs** — the campaign report format and user documentation for
  running a VEM search.

**Rationale.** Tying UI to the search core blocks merging functional code: the
search-core changes (`period_basis`, the beat/anchor resolver, the same-body
bypass, the loader class) are reviewable and shippable on their own and gated
by golden/feasibility tests. The CLI/viz/docs have a different gate (end-to-end
behaviour, rendering) and a different review surface. Coupling them would hold
correct, tested search-core code hostage to UX iteration. M8-UX depends on
M8-Core landing first.

---

## 7. Definition of done

0. `tests/search/test_sequence_cell_basis.py` — `Cell.period_basis` carries the
   anchor pair without mutating `period_k`; `id` token appended only when set;
   existing 2-body ids byte-identical. **Green.**
1. `tests/search/test_optimize_multibody.py` — EMEEVE anchor-pair period
   resolves to 6.41 yr; basis-None ≥3-body falls back to the 6.406-yr natural
   beat; scales with `period_k`; 2-body basis-None path byte-identical. **Green.**
2. `tests/search/test_sequence_multibody.py` — the EMEEVE E-E loop leg is
   bypassed (returns True, not falsely rejected); distinct-pair coplanar
   Tisserand baseline unchanged. **Green.**
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

## 8. Resolved design decisions

These were the open questions in the original draft. Each is now **decided**;
the plan above implements the decision. Recorded here so reviewers see the
rationale without re-litigating.

1. **Period convention — REJECT silent Reading A.** Do **not** reinterpret
   `Cell.period_k` as a "beat multiple" and do **not** have the loader rewrite
   the catalogue's anchor-pair `period.k`. That destroys traceability (a Cell
   from `vem-emeeve-3syn` showing `period_k=1` looks like a bug vs the YAML's
   `k=3`). **Decision:** add the optional `period_basis: tuple[str, str] | None`
   field to the frozen `Cell` (§2 Task 2.0). The Cell echoes the anchor pair and
   the sourced `period_k` verbatim; `_target_period_sec` dispatches on
   `period_basis` (anchor-pair `T_syn(*basis)*k`), falling back to the natural
   beat only when `period_basis is None and len(bodies) >= 3`, and preserving
   the 2-body path byte-for-byte. `Cell.id` gains a `|p<AB>` token only when
   `period_basis` is set, so existing 2-body ids are stable.
2. **Loader shape — KEEP the new `CONSTRUCTIBLE_MULTIBODY` exclusion reason.**
   Minimal-diff: it satisfies the census ratchet and keeps the ≥3-body rows out
   of the 2-body V∞ rediscovery gauntlet (Task 5.3 guard) without a wider
   loader-API change. The classification is N-agnostic (`len(bodies) >= 3` +
   valid period block), not a hardcoded VEM bodyset (§5 Task 5.2).
3. **Milestone split — SPLIT into M8-Core and M8-UX.** This plan is **M8-Core**
   (search-core: routing, Tisserand, VEM rediscovery gate, loader admission).
   CLI, visualisation, reporting, and docs are carved into a separate **M8-UX**
   milestone (§6 "Out of scope → M8-UX"). Tying UI to the search core would
   block merging correct, tested functional code.
4. **EMEEVE sequence — anchor on the catalogue string `E-M-E-E-V-E`.** The gate
   anchors on the catalogue's sourced `sequence_canonical = "E-M-E-E-V-E"`
   (`catalogue.yaml:1777`) for golden discipline. The spec §13.8 worked cell-id
   `E-V-M-E-M-E` is **illustrative only** and is not used as a test anchor.
5. **Beat-tuple robustness — accept `multi_body_beat_days(...)[0]` for
   M8-Core.** For the basis-None fallback the lowest-mismatch tuple is
   unambiguously `(4,3)` ⇒ ~6.406 yr at the default `k_max=6`/`tol_frac=0.02`
   (resonance gate `test_resonance.py:54`), reliable for the natural ~6.4-yr
   beat. Threading non-default `k_max`/`tol_frac` for higher-period VEM branches
   (12.8 yr / 32 yr, spec §3 line 41) is **deferred to M-ED**. Note: the primary
   EMEEVE gate uses the anchor-pair path (decision 1), so the beat fallback is
   only exercised for cells with no declared basis.
