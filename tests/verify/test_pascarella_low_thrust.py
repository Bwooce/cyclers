"""Positive control validation for the Pascarella et al. 2022 (AAS 22-015) low-thrust cycler."""

from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path

import pytest
import yaml  # type: ignore[import-untyped]

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.verify import construct_real_ephemeris_cycler, verify_low_thrust_feasibility

_J2000_EPOCH = datetime(2000, 1, 1, 12, 0, 0, tzinfo=UTC)


def dt_to_t_sec(dt: datetime) -> float:
    """Convert UTC datetime to seconds since J2000."""
    return (dt - _J2000_EPOCH).total_seconds()


@pytest.fixture(scope="module")
def astropy_ephem() -> Ephemeris:
    return Ephemeris(model="astropy")


def test_pascarella_low_thrust_reproduction(astropy_ephem: Ephemeris) -> None:
    # 1. Load the sanchez-net-2022-eem-cycler1 entry from catalogue.yaml
    path = Path(__file__).resolve().parents[2] / "data" / "catalogue.yaml"
    assert path.exists()
    with path.open() as fh:
        raw = yaml.safe_load(fh)
    entry = next(e for e in raw if e.get("id") == "sanchez-net-2022-eem-cycler1")

    # 2. Construct the real ephemeris cycler on DE440 at the launch epoch (2032-05-11)
    launch_dt = datetime(2032, 5, 11, tzinfo=UTC)
    t_start = dt_to_t_sec(launch_dt)

    cycler = construct_real_ephemeris_cycler(entry, astropy_ephem, t_start)
    assert cycler is not None

    # 3. Verify low-thrust feasibility using the Pascarella NEXT engine parameters:
    # - wet mass: 500 kg
    # - Isp: 4155 s
    # - thrust: 0.235 N (0.000235 kN)
    # - segment count: 20
    # - targeting maintenance ΔV: ~163 m/s (from 2 kg propellant spent over the Targeting phase)
    # Since the catalogue row lists delta_v_kms = 0.005 (5 m/s nominal), check both
    # the catalogue's nominal 5 m/s and the paper's 163 m/s targeting ΔV.
    isp_s = 4155.0
    tmax_kn = 0.000235
    m0_kg = 500.0

    # Check nominal 5 m/s (0.005 km/s) is feasible:
    feasible_nominal = verify_low_thrust_feasibility(
        cycler,
        maintenance_dv_kms=0.005,
        isp_s=isp_s,
        tmax_kn=tmax_kn,
        m0_kg=m0_kg,
        n_segments=20,
    )
    assert feasible_nominal is True

    # Check paper's total targeting ΔV of 163 m/s (0.163 km/s) is also feasible:
    feasible_targeting = verify_low_thrust_feasibility(
        cycler,
        maintenance_dv_kms=0.163,
        isp_s=isp_s,
        tmax_kn=tmax_kn,
        m0_kg=m0_kg,
        n_segments=20,
    )
    assert feasible_targeting is True
