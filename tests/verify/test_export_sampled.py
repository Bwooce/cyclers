"""Unit tests for the viz-2c sampled-trajectory exporter (task #146).

These tests exercise the exporter's MECHANICS on constructed inputs — decimation
error bounds, the km->AU unit conversion, and the emitted JSON shape — without the
slow real-DE440 propagation (that is the production solver the exporter calls; it
is covered by ``test_aldrin_v2_powered.py``). No golden values: every expectation
is either a property of the construction (a known analytic arc) or a structural
invariant of the site's ``SampledTrajectory`` contract.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.constants import AU_KM
from scripts.export_sampled_trajectories import (
    _strictly_increasing,
    build_payload,
    decimate_curvature,
)


def _dense_circle(n: int, radius_km: float = AU_KM) -> np.ndarray:
    """A dense ``(n, 7)`` sample array tracing a planar circle of given radius.

    The chord error of any decimation of a circular polyline has a known
    closed form (sagitta), so this is a constructed reference, not a guess.
    Columns: ``[t_sec, x, y, z, vx, vy, vz]``; velocities are left zero (the
    decimator only reads positions).
    """
    out = np.zeros((n, 7), dtype=np.float64)
    for i in range(n):
        theta = 2.0 * np.pi * i / (n - 1)
        out[i, 0] = float(i)  # strictly increasing time
        out[i, 1] = radius_km * np.cos(theta)
        out[i, 2] = radius_km * np.sin(theta)
    return out


# ---------------------------------------------------------------------------
# decimate_curvature (Douglas-Peucker, seam-aware)
# ---------------------------------------------------------------------------

_TIGHT_TOL_KM = 1.0e3  # ~6.7e-6 AU — tight enough to keep many circle points


def test_decimate_keeps_endpoints_and_caps_count() -> None:
    dense = _dense_circle(1000)
    decimated, _err, _meta = decimate_curvature(dense, max_points=100, tol_km=_TIGHT_TOL_KM)
    assert decimated.shape[0] <= 100
    # First and last propagated states are always retained.
    assert np.array_equal(decimated[0], dense[0])
    assert np.array_equal(decimated[-1], dense[-1])


def test_decimate_error_within_tolerance() -> None:
    """Douglas-Peucker guarantees the decimated polyline stays within tol of every
    dropped point: the measured max chord error must not exceed the (escalated)
    tolerance actually used."""
    dense = _dense_circle(1500)
    decimated, err_km, meta = decimate_curvature(dense, max_points=500, tol_km=_TIGHT_TOL_KM)
    assert err_km <= meta["tol_km_used"]
    assert decimated.shape[0] <= 500


def test_decimate_looser_tolerance_keeps_fewer_points() -> None:
    """A looser chord tolerance keeps fewer points (curvature-adaptive: only the
    bends survive)."""
    dense = _dense_circle(2000)
    tight, _e1, _m1 = decimate_curvature(dense, max_points=2000, tol_km=1.0e3)
    loose, _e2, _m2 = decimate_curvature(dense, max_points=2000, tol_km=1.0e5)
    assert loose.shape[0] < tight.shape[0]


def test_decimate_detects_and_preserves_seams() -> None:
    """A constructed seam jump (craft snapping to a planet) is detected, counted,
    excluded from the chord-error metric, and its endpoints are kept so the seam
    survives in the output."""
    run_a = _dense_circle(400)
    run_b = _dense_circle(400)
    # Shift run_b far away in time and space to create one big seam between them.
    run_b[:, 0] += 1000.0  # keep times strictly increasing
    run_b[:, 1] += 3.0 * AU_KM  # a >0.3 AU jump at the join
    dense = np.vstack([run_a, run_b])

    _dec, err_km, meta = decimate_curvature(dense, max_points=500, tol_km=_TIGHT_TOL_KM)
    assert meta["seam_count"] == 1
    assert meta["max_seam_jump_au"] > 0.3
    # The seam jump (3 AU) must NOT pollute the within-run chord error.
    assert err_km / AU_KM < 0.3


def test_decimate_rejects_tiny_budget() -> None:
    with pytest.raises(ValueError, match="max_points must be >= 2"):
        decimate_curvature(_dense_circle(100), max_points=1)


# ---------------------------------------------------------------------------
# _strictly_increasing (seam-duplicate-epoch handling)
# ---------------------------------------------------------------------------


def test_strictly_increasing_nudges_seam_duplicates() -> None:
    """An equal-time seam pair (arriving leg endpoint + next leg start at the same
    encounter instant) is nudged into strict monotonicity, the site's bracketing
    precondition, with a sub-millisecond bump that preserves order."""
    times = [0.0, 10.0, 10.0, 20.0, 19.0]  # one duplicate + one out-of-order
    out = _strictly_increasing(times, eps=1.0e-3)
    assert all(out[i] > out[i - 1] for i in range(1, len(out)))
    # Untouched entries are unchanged; only the violators move, minimally.
    assert out[0] == 0.0 and out[1] == 10.0 and out[3] == 20.0
    assert out[2] == pytest.approx(10.0 + 1e-3)


def test_strictly_increasing_noop_when_already_monotone() -> None:
    times = [0.0, 1.0, 2.5, 100.0]
    assert _strictly_increasing(times) == times


# ---------------------------------------------------------------------------
# build_payload (mechanics via monkeypatched propagation — no slow DE440)
# ---------------------------------------------------------------------------


class _StubCycler:
    period = 780.0 * 86400.0


def test_build_payload_shape_and_units(monkeypatch: pytest.MonkeyPatch) -> None:
    """The emitted payload matches the site's SampledTrajectory contract and the
    km->AU conversion is exact (positions are dense_km / AU_KM, times verbatim).
    """
    import scripts.export_sampled_trajectories as mod

    dense = _dense_circle(800)
    monkeypatch.setattr(mod, "propagate_aldrin_powered", lambda **kw: (dense, _StubCycler()))

    payload = build_payload(entry_id="aldrin-classic-em-k1-outbound", n_laps=3, max_points=300)

    # Contract fields the site reads.
    assert payload["kind"] == "sampled"
    assert payload["frame"] == "eclipJ2000"
    assert isinstance(payload["fidelity"], str) and payload["fidelity"]
    assert isinstance(payload["provenance"], str) and payload["provenance"]
    assert "DE440" in payload["fidelity"]

    times = payload["timesSec"]
    pos = payload["positionsAU"]
    assert len(times) == len(pos)
    assert len(times) <= 300
    # Parallel arrays of the right inner shape.
    assert all(len(p) == 3 for p in pos)
    # Strictly increasing times (precondition of the site's bracketing search).
    assert all(times[i] > times[i - 1] for i in range(1, len(times)))

    # km->AU is the ONLY position transform: endpoints are always kept, so the
    # first/last emitted AU positions must equal the dense endpoints / AU_KM.
    assert pos[0][0] == pytest.approx(dense[0, 1] / AU_KM, rel=1e-15)
    assert pos[0][1] == pytest.approx(dense[0, 2] / AU_KM, rel=1e-15)
    assert pos[-1][0] == pytest.approx(dense[-1, 1] / AU_KM, rel=1e-15)
    assert pos[-1][1] == pytest.approx(dense[-1, 2] / AU_KM, rel=1e-15)
    # Times are emitted verbatim (seconds since J2000), no conversion.
    assert times[0] == pytest.approx(dense[0, 0])
    assert times[-1] == pytest.approx(dense[-1, 0])


def test_build_payload_meta_reports_bounded_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import scripts.export_sampled_trajectories as mod

    dense = _dense_circle(1200)
    monkeypatch.setattr(mod, "propagate_aldrin_powered", lambda **kw: (dense, _StubCycler()))
    payload = build_payload(entry_id="aldrin-classic-em-k1-outbound", n_laps=3, max_points=500)
    meta = payload["_meta"]
    assert meta["dense_points"] == 1200
    assert meta["emitted_points"] <= 500
    # The error reported in _meta equals the decimator's own measurement.
    _dec, err_km, _m = decimate_curvature(dense, 500)
    assert meta["max_chord_error_km"] == pytest.approx(err_km)
    assert meta["max_chord_error_au"] == pytest.approx(err_km / AU_KM)
    assert meta["cycler_period_days"] == pytest.approx(780.0)
    assert meta["decimation"].startswith("douglas-peucker")


def test_build_payload_rejects_unknown_id() -> None:
    with pytest.raises(ValueError, match="unsupported entry_id"):
        build_payload(entry_id="some-other-row")
