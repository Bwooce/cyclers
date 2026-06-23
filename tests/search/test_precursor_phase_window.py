"""#307 Task 3: cycler-cadence target phase-window auto-derivation.

The precursor matcher's alignment score (``_epoch_alignment_score``) already
accepts a target phase window; Task 3 auto-derives it from the cycler row's
published ``validity_window`` when the caller does not override. The score itself
stays informational/ungated — these tests pin the derivation only.
"""

from __future__ import annotations

from types import SimpleNamespace

from cyclerfinder.data.catalog import CATALOGUE_PATH, load_catalog
from cyclerfinder.search.precursor_matcher import _phase_window_for_entry


def test_override_always_wins() -> None:
    entry = SimpleNamespace(raw={"validity_window": {"start": "A", "end": "B"}})
    assert _phase_window_for_entry(entry, ("x", "y")) == ("x", "y")  # type: ignore[arg-type]


def test_auto_derive_from_validity_window() -> None:
    entry = SimpleNamespace(
        raw={"validity_window": {"start": "2023-01-03T00:00:00Z", "end": "2023-06-13T00:00:00Z"}}
    )
    assert _phase_window_for_entry(entry, None) == (  # type: ignore[arg-type]
        "2023-01-03T00:00:00Z",
        "2023-06-13T00:00:00Z",
    )


def test_none_when_no_window() -> None:
    assert _phase_window_for_entry(SimpleNamespace(raw={}), None) is None  # type: ignore[arg-type]
    assert _phase_window_for_entry(SimpleNamespace(raw=None), None) is None  # type: ignore[arg-type]
    # partial window (missing end) is not usable
    assert (
        _phase_window_for_entry(SimpleNamespace(raw={"validity_window": {"start": "A"}}), None)  # type: ignore[arg-type]
        is None
    )


def test_real_catalogue_row_with_validity_window() -> None:
    """A real row that publishes a validity_window auto-derives a (start, end)."""
    cat = load_catalog(CATALOGUE_PATH)
    entry = cat.by_id["aldrin-4-3-2-establishment"]
    window = _phase_window_for_entry(entry, None)
    assert window is not None
    assert len(window) == 2 and all(isinstance(s, str) and s for s in window)
