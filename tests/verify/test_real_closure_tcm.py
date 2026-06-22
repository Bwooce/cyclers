"""M7 wiring (#423) — real_closure horizon-TCM populates and feeds the #424 verdict.

The M6b ``RealClosureResult.horizon_tcm_mps`` was hard-coded ``0.0``. M7 wires
:func:`cyclerfinder.verify.real_closure._compute_horizon_tcm` (the position-targeted
continuous maintenance chain) into it under the opt-in ``compute_tcm`` flag. These
tests exercise the helper directly on a real-ephemeris E-M-E cycler (avoiding the
launch-window machinery) and confirm the result drives the #424
``v3_class_split_verdict`` programmatically — retiring the manual #175 convention.
"""

from __future__ import annotations

import numpy as np
import pytest

rebound = pytest.importorskip("rebound")

from cyclerfinder.core.ephemeris import Ephemeris  # noqa: E402
from cyclerfinder.model.cycler import Cycler, Encounter  # noqa: E402
from cyclerfinder.verify.dv_band_acceptance import v3_class_split_verdict  # noqa: E402
from cyclerfinder.verify.real_closure import _compute_horizon_tcm  # noqa: E402

_DAY_S = 86400.0


def _eme_cycler(ephem: Ephemeris) -> Cycler:
    """A minimal real-ephemeris E-M-E cycler (3 encounters); vinf vectors unused here."""
    t0 = 27.0 * 365.25 * _DAY_S
    epochs = [t0, t0 + 210.0 * _DAY_S, t0 + 480.0 * _DAY_S]
    bodies = ["E", "M", "E"]
    zero = np.zeros(3)
    encounters = []
    for body, t in zip(bodies, epochs, strict=True):
        r, v = (np.asarray(x, dtype=np.float64) for x in ephem.state(body, t))
        encounters.append(Encounter(body=body, t=t, r=r, v_planet=v, vinf_in=zero, vinf_out=zero))
    return Cycler(bodies=bodies, period=epochs[-1] - epochs[0], encounters=encounters, legs=[])


@pytest.mark.slow
def test_compute_horizon_tcm_populates_and_feeds_v3_verdict() -> None:
    """`_compute_horizon_tcm` returns a finite N-cycle TCM that the #424 verdict consumes.

    The 2-cycle E-M-E walk (E M E M E, including the inter-cycle home flyby) converges
    in Sun-cruise mode, giving a finite horizon TCM and a length-2 per-cycle tuple. The
    value is then fed to v3_class_split_verdict (#424) as a ballistic row — the verdict
    fires programmatically (PASS/FAIL against the 120 m/s/7cyc bar), proving the wiring
    that the manual #175 convention used to do by hand."""
    ephem = Ephemeris("astropy")
    cyc = _eme_cycler(ephem)

    horizon_mps, per_cycle = _compute_horizon_tcm(cyc, 2, ephem)

    assert np.isfinite(horizon_mps), "Sun-cruise E-M-E should converge to a finite TCM"
    assert horizon_mps >= 0.0
    assert len(per_cycle) == 2
    assert per_cycle[0] == pytest.approx(horizon_mps / 2)

    # #424 fires programmatically on the measured TCM (ballistic class, 120 m/s bar).
    verdict = v3_class_split_verdict(horizon_mps, n_cycles=2, dv_band="essentially_ballistic")
    assert verdict.cls == "ballistic"
    assert isinstance(verdict.passed, bool)
    assert verdict.bar_mps == pytest.approx(120.0 * 2 / 7)
