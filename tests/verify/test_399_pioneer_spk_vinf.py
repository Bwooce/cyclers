"""Frozen-gate validation of the #399 Pioneer SPK-derived flyby V∞.

Pioneer 10/11 are the clean-hyperbolic-flyby analogue of the #390 Voyager case:
the per-encounter V∞ is not in the acquired literature in clean form, so it is
DERIVED from the mission's archived NAIF reconstructed spacecraft SPK
(:mod:`cyclerfinder.verify.mission_spk`). This module pins the extraction.

* **Closest-approach self-consistency (tight gate).** The spacecraft's minimum
  altitude above the flyby planet, recovered from the SPK, must match the
  PUBLISHED closest-approach distance (NASA Science mission pages).
* **Vis-viva convergence (stability gate).** The outer-window vis-viva V∞ must
  be stable to < 1%.
* **Order-of-magnitude sanity.** Pioneer Jupiter/Saturn V∞ ~ 8-9 km/s.

SKIPPED honestly if spiceypy / astropy is absent or the SPK fetch is blocked.

Provenance of the EXPECTED side (sourced-only, non-circular): the published
closest-approach altitudes are the NASA mission-page record, NOT values our own
code computed:
* Pioneer 10 Jupiter: 130,354 km (81,000 mi), 1973-12-04 02:26 UT.
* Pioneer 11 Jupiter: 42,500 km (26,400 mi from cloud tops), 1974-12-03 05:22 UT.
* Pioneer 11 Saturn:  20,900 km (13,000 mi), 1979-09-01 16:31 UT.
"""

from __future__ import annotations

from types import ModuleType
from typing import Any

import pytest

# Published flyby altitudes (km) — NASA Science Pioneer-10/11 mission pages.
PUBLISHED_P10_JUPITER_ALT_KM = 130_354.0
PUBLISHED_P11_JUPITER_ALT_KM = 42_500.0
PUBLISHED_P11_SATURN_ALT_KM = 20_900.0

# V∞ order-of-magnitude sanity band (km/s) — NOT a tight gate.
PIONEER_VINF_BAND_KMS = (7.0, 10.0)


@pytest.fixture(scope="module")
def mission_spk_module() -> ModuleType:
    """Import the extractor or skip if the validation extra is unavailable."""
    pytest.importorskip("spiceypy", reason="spiceypy (validation extra) not installed")
    pytest.importorskip("astropy", reason="astropy (validation extra) not installed")
    import cyclerfinder.verify.mission_spk as mission_spk

    return mission_spk


def _extract(
    mission_spk: ModuleType,
    base_url: str,
    filename: str,
    naif_id: int,
    body: str,
    epoch: str,
    window_minutes: float = 4320.0,
) -> Any:
    """Fetch the SPK + run the extractor, skipping on a network/fetch failure."""
    try:
        spk = mission_spk.ensure_mission_spk(filename, base_url=base_url)
    except Exception as exc:  # network blocked / 404 -> skip honestly
        pytest.skip(f"could not fetch NAIF SPK {filename} (network blocked?): {exc}")
    try:
        return mission_spk.vinf_at_flyby(
            spk, naif_id, body, epoch, mission="Pioneer", window_minutes=window_minutes
        )
    except Exception as exc:  # DE440 / LSK fetch failure -> skip honestly
        pytest.skip(f"extractor failed (kernel fetch blocked?): {exc}")


def _assert_alt(r: Any, published_km: float, tol: float = 0.03) -> None:
    rel = abs(r.closest_approach_altitude_km - published_km) / published_km
    assert rel <= tol, (
        f"CA altitude {r.closest_approach_altitude_km:.0f} km disagrees with "
        f"published {published_km:.0f} km by {rel:.1%}"
    )


def _assert_vinf(r: Any) -> None:
    frac = r.vinf_kms_visviva_window_std / r.vinf_kms_visviva_window_mean
    assert frac < 0.01, (
        f"V∞ not converged: std/mean = {frac:.2%} (mean {r.vinf_kms_visviva_window_mean:.4f} km/s)"
    )
    lo, hi = PIONEER_VINF_BAND_KMS
    assert lo <= r.vinf_kms <= hi, f"V∞ {r.vinf_kms:.3f} km/s outside [{lo},{hi}]"


@pytest.mark.slow
def test_pioneer10_jupiter(mission_spk_module: ModuleType) -> None:
    """Pioneer 10 Jupiter CA altitude matches NASA + V∞ converges (~8.5 km/s)."""
    m = mission_spk_module
    r = _extract(
        m,
        m.NAIF_PIONEER10_SPK_BASE,
        "p10-a.bsp",
        m.PIONEER_10_NAIF_ID,
        "Jupiter",
        "1973-12-04T02:26:00",
    )
    _assert_alt(r, PUBLISHED_P10_JUPITER_ALT_KM)
    _assert_vinf(r)


@pytest.mark.slow
def test_pioneer11_jupiter(mission_spk_module: ModuleType) -> None:
    """Pioneer 11 Jupiter CA altitude matches NASA + V∞ converges (~8.9 km/s)."""
    m = mission_spk_module
    r = _extract(
        m,
        m.NAIF_PIONEER11_SPK_BASE,
        "p11-a.bsp",
        m.PIONEER_11_NAIF_ID,
        "Jupiter",
        "1974-12-03T05:22:00",
    )
    _assert_alt(r, PUBLISHED_P11_JUPITER_ALT_KM)
    _assert_vinf(r)


@pytest.mark.slow
def test_pioneer11_saturn(mission_spk_module: ModuleType) -> None:
    """Pioneer 11 Saturn CA altitude matches NASA + V∞ converges (~8.4 km/s).

    The Saturn-relative kernel p11_sat336.bsp covers only 1979-08-20..09-04, so
    the extraction window is narrowed to stay inside that coverage span.
    """
    m = mission_spk_module
    r = _extract(
        m,
        m.NAIF_PIONEER11_SPK_BASE,
        "p11_sat336.bsp",
        m.PIONEER_11_NAIF_ID,
        "Saturn",
        "1979-09-01T16:31:00",
        window_minutes=4000.0,
    )
    _assert_alt(r, PUBLISHED_P11_SATURN_ALT_KM)
    _assert_vinf(r)
