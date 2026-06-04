"""M7 ledger-backed discovery runner tests — plan §3.5 / §4.

Exercises :func:`cyclerfinder.data.discover.discover`:

* runner mechanics (records every attempted cell; resume skips them;
  V3 stays off by default);
* the end-to-end known-rediscovery yield (**xfail** under the M5
  optimiser regression, task #54).

The mechanics tests use a tiny, fast enumeration (``k_synodic=1``,
``l_max=2``, ``use_de=False``, ``n_starts=1``) so they exercise the
ledger/skip/record plumbing without the full M5 search cost.
"""

from __future__ import annotations

from pathlib import Path

import pytest

import cyclerfinder.data.discover as discover_mod
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.data.discover import discover
from cyclerfinder.data.ledger import Ledger, LedgerLoader
from cyclerfinder.search.optimize import optimise_cell_idealized
from cyclerfinder.search.sequence import feasible_cells

_VALID_TERMINAL = {"solved", "pruned", "failed", "searched"}


def _tiny_run(ledger_path: Path, ephem: Ephemeris) -> list[tuple[object, object, str]]:
    return list(
        discover(
            ("E", "M"),
            1,
            7.0,
            str(ledger_path),
            ephem=ephem,
            l_max=2,
            use_de=False,
            n_starts=1,
            seed=0,
        )
    )


def test_discover_writes_ledger_for_every_cell_attempted(tmp_path: Path) -> None:
    """Every enumerated cell gets a terminal ledger entry."""
    ephem = Ephemeris(model="circular")
    led = tmp_path / "ledger.jsonl"
    _tiny_run(led, ephem)

    ledger = Ledger(str(led))
    cells = list(feasible_cells(("E", "M"), l_max=2, k_max=1, n_max=0, vinf_cap=7.0, ephem=ephem))
    assert cells, "enumeration produced no feasible cells"
    for cell in cells:
        assert ledger.has(cell.id), f"cell {cell.id} not recorded"
        assert ledger.get(cell.id).status in _VALID_TERMINAL


def test_discover_resumes_from_existing_ledger(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A second run skips already-recorded cells (no optimiser calls)."""
    ephem = Ephemeris(model="circular")
    led = tmp_path / "ledger.jsonl"
    _tiny_run(led, ephem)
    n_after_first = len(Ledger(str(led)))
    assert n_after_first > 0

    calls: list[int] = []
    real = optimise_cell_idealized

    def spy(*args: object, **kwargs: object) -> object:
        calls.append(1)
        return real(*args, **kwargs)  # type: ignore[arg-type]

    monkeypatch.setattr(discover_mod, "optimise_cell_idealized", spy)
    new_yields = _tiny_run(led, ephem)

    assert calls == [], "optimiser was called on resume despite recorded cells"
    assert new_yields == [], "resume yielded new results for already-done cells"
    assert len(Ledger(str(led))) == n_after_first


def test_discover_records_signature_hash_in_ledger(tmp_path: Path) -> None:
    """Every solved ledger entry carries a non-empty signature hash tuple."""
    ephem = Ephemeris(model="circular")
    led = tmp_path / "ledger.jsonl"
    _tiny_run(led, ephem)
    for entry in LedgerLoader(str(led)):
        if entry.status == "solved":
            assert entry.signature_hashes != ()


def test_discover_skips_v3_when_disabled(tmp_path: Path) -> None:
    """``enable_v3=False`` (default) never triggers the M5 stub's NotImplementedError."""
    ephem = Ephemeris(model="circular")
    led = tmp_path / "ledger.jsonl"
    # Must not raise.
    _tiny_run(led, ephem)


@pytest.mark.xfail(
    strict=False,
    reason=(
        "Pre-existing M5 optimiser regression (task #54): the E-M k=2 "
        "optimiser returns no constraint-satisfying cycler, so discover "
        "yields no 'known' match for the 5.65 km/s Russell cycler. Same "
        "root cause as test_2syn_em_rediscovers_5_65_kms_earth. Flip "
        "strict=True once #54 lands. "
        "NOTE (2026-06-04): the 5.65 km/s anchor is unverified-provenance "
        "(catalogue data_gap vinf_kms_at_encounters, s1l1-2syn-em-cpom): "
        "traces only to spec.md §9; unconfirmed in Patel 2019 / McConaghy "
        "2006 / Sanchez Net 2022 — see docs/notes/s1l1-target-topology-mining.md."
    ),
)
def test_discover_em_k2_yields_known_for_2syn(tmp_path: Path) -> None:
    """End-to-end: the ~5.65 km/s E-M k=2 result matches ``s1l1-2syn-em-cpom``."""
    led = tmp_path / "ledger.jsonl"
    yields = list(discover(("E", "M"), 2, 7.0, str(led), seed=0))
    known = [
        (result, mr)
        for (result, mr, _level) in yields
        if abs(result.best_cycler.max_vinf() - 5.65) < 0.3
    ]
    assert known, "no E-M k=2 result near 5.65 km/s was found (task #54)"
    _result, match_result = known[0]
    assert match_result.outcome == "known"
    assert match_result.entry is not None
    assert match_result.entry.id == "s1l1-2syn-em-cpom"


def test_discover_accepts_multirev_params(tmp_path: Path) -> None:
    """The plumbing contract: ``discover`` accepts the multi-rev sweep
    params (``n_max`` / ``branch_set`` / ``max_cells``), runs without
    raising, and bounds the feasible-cell stream by ``max_cells``.

    This asserts only the plumbing — *not* that a cell solves. Whether
    any E-M k=2 cell yields a constraint-satisfying cycler is the M5
    optimiser regression's concern (task #54), exercised separately by
    :func:`test_snlm_sweep_rediscovers_a_sourced_anchor` (xfail).
    """
    ledger = tmp_path / "ledger.jsonl"
    # Must not raise; with the bounded sweep only solved cells are yielded
    # (may be empty under the M5 regression), but every *enumerated* cell
    # within the bound is recorded to the ledger.
    list(
        discover(
            bodies=("E", "M"),
            k_synodic=2,
            vinf_cap=8.0,
            ledger_path=str(ledger),
            l_max=4,
            n_max=1,
            branch_set=("single", "low"),
            max_cells=3,
        )
    )
    recorded = Ledger(str(ledger))
    assert len(recorded) == 3, "max_cells must bound the swept (recorded) cell stream"
    # At least one recorded cell must be multi-rev (n_revs >= 1), proving
    # branch_set/n_max actually widened the enumeration.
    cells = list(
        feasible_cells(
            ("E", "M"),
            l_max=4,
            k_max=2,
            n_max=1,
            vinf_cap=8.0,
            ephem=Ephemeris(model="circular"),
            branch_set=("single", "low"),
        )
    )[:3]
    assert any("r1" in c.id or "bl" in c.id for c in cells), (
        "expected a multi-rev cell in the bounded sweep"
    )


def test_discover_accepts_ephemeris_optimiser(tmp_path: Path) -> None:
    """Part 2 plumbing: ``discover(optimiser='ephemeris')`` routes cells to the
    real-ephemeris ``optimise_cell_ephemeris`` without raising. Without
    V-infinity targets no launch epoch resolves, so every cell is recorded
    ("searched") but none is yielded — the wiring contract, not a closure claim.
    """
    ledger = tmp_path / "ledger.jsonl"
    out = list(
        discover(
            bodies=("E", "M"),
            k_synodic=1,
            vinf_cap=8.0,
            ledger_path=str(ledger),
            ephem=Ephemeris(model="circular"),
            l_max=2,
            n_max=0,
            optimiser="ephemeris",
            n_starts=1,
            seed=0,
        )
    )
    assert out == [], "ephemeris mode without V-inf targets resolves no epoch -> nothing solved"
    cells = list(
        feasible_cells(
            ("E", "M"), l_max=2, k_max=1, n_max=0, vinf_cap=8.0, ephem=Ephemeris(model="circular")
        )
    )
    assert cells
    recorded = Ledger(str(ledger))
    for c in cells:
        assert recorded.has(c.id), f"cell {c.id} not recorded in ephemeris mode"


@pytest.mark.slow
@pytest.mark.xfail(
    strict=False,
    reason=(
        "Circular-coplanar model limitation, confirmed exhaustively: the "
        "S1L1 / SnLm family is not hostable in the idealised model (same as "
        "Aldrin). Both the bare multi-rev sweep AND the correct E-M-E-E "
        "topology with the sourced [154, 379, 1030] d S1/L1 seed land at "
        "V_inf_E ~25-39 km/s (see scripts/characterise_s1l1_emee.py) — the "
        "154-d E->M leg is near-hyperbolic in circular-coplanar, so the "
        "blocker is the MODEL, not the topology. discover() now also offers "
        "optimiser='ephemeris' (real-DE440), but closing S1L1/SnLm there "
        "additionally needs multi-rev support in the maintenance engine "
        "behind optimise_cell_ephemeris (currently single-rev only). The "
        "0.3 km/s tolerance is the sourced-anchor bound and is NOT loosened. "
        "Flips to a pass once the real-eph optimiser handles the L1 multi-rev "
        "leg; flip strict=True then. "
        "NOTE (2026-06-04): the 5.65 km/s anchor is unverified-provenance "
        "(catalogue data_gap vinf_kms_at_encounters, s1l1-2syn-em-cpom): "
        "traces only to spec.md §9; unconfirmed in Patel 2019 / McConaghy "
        "2006 / Sanchez Net 2022 — see docs/notes/s1l1-target-topology-mining.md."
    ),
)
def test_snlm_sweep_rediscovers_a_sourced_anchor(tmp_path: Path) -> None:
    """A bounded 2-synodic E-M multi-rev sweep must surface at least one
    closed, constraint-satisfying cycler whose Earth V-infinity matches a
    sourced SnLm anchor (5.65 km/s) within the gauntlet tolerance."""
    import numpy as np

    ledger = tmp_path / "ledger.jsonl"
    matched = False
    for opt_result, _match, _level in discover(
        bodies=("E", "M"),
        k_synodic=2,
        vinf_cap=8.0,
        ledger_path=str(ledger),
        l_max=4,
        n_max=1,
        branch_set=("single", "low"),
        max_cells=24,
    ):
        if not opt_result.constraints_satisfied:
            continue
        for enc in opt_result.best_cycler.encounters:
            if enc.body == "E":
                v = max(
                    float(np.linalg.norm(enc.vinf_in)),
                    float(np.linalg.norm(enc.vinf_out)),
                )
                if abs(v - 5.65) < 0.3:
                    matched = True
    assert matched, "no swept multi-rev cell matched the 5.65 km/s SnLm Earth anchor"
