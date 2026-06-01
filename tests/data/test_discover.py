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
    real = discover_mod.optimise_cell_idealized

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
        "strict=True once #54 lands."
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
