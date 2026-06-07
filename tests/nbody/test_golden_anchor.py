"""N-body GOLDEN GATE 3: planet-state ingestion == Ephemeris('astropy') (design §5.4).

GOLDEN: EXPECTED = DE440 itself (Ephemeris('astropy').state). This is the §0
anchor: it proves the time/frame conversion before any spacecraft number is
trusted. The sourced side is the ephemeris.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.nbody.forces import ingest_planet_state


@pytest.mark.slow
@pytest.mark.parametrize("body", ["E", "M", "J"])
@pytest.mark.parametrize("t_sec", [0.0, 365.25 * 86400.0, 10 * 365.25 * 86400.0])
def test_ingested_planet_state_matches_ephemeris(body: str, t_sec: float) -> None:
    ephem = Ephemeris("astropy")
    r_ref, v_ref = ephem.state(body, t_sec)
    r_got, v_got = ingest_planet_state(body, t_sec, ephem)
    assert np.allclose(r_got, r_ref, atol=1e-6)  # numerical precision
    assert np.allclose(v_got, v_ref, atol=1e-9)
