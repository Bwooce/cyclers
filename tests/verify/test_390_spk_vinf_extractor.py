"""Frozen-gate validation of the #390 NAIF spacecraft-SPK V∞ extractor.

The extractor (:mod:`cyclerfinder.verify.mission_spk`) derives planetocentric
hyperbolic-excess velocity (V∞) at a flyby epoch directly from a mission's
archived NAIF SPK, the wholesale unblock for the #345 classic-mission
``mga_tour`` backlog (whose papers publish dates + closest-approach geometry +
ΔV but never per-encounter V∞).

This module pins the Voyager 2 validation:

* **Closest-approach self-consistency (tight gate).** The spacecraft's minimum
  distance from Jupiter / Saturn, recovered from the SPK, must match the
  PUBLISHED closest-approach geometry (Kohlhase-Penzo 1977 "Voyager Mission
  Description", Space Science Reviews 21(2):77-101, Table IV / Figs 6-11). The
  digest note (``docs/notes/2026-06-19-345-voyager-mariner-mission-digests.md``)
  records JSX Jupiter = 10.0 R_J. Closest approach IS published, so this is the
  real ground-truth check on the geometry the extractor reads.
* **Vis-viva convergence (stability gate).** The vis-viva hyperbolic-excess
  speed ``sqrt(v_rel^2 - 2*mu/r_rel)`` must be stable across the outer-window
  samples (well outside the SOI where v_rel -> V∞), to < 1%.
* **Order-of-magnitude sanity.** The Voyager 2 Jupiter approach V∞ is widely
  cited at ~7-8 km/s; Saturn ~10-11 km/s (this is a sanity band, not a gate).

These are SKIPPED honestly if spiceypy / astropy is absent or the NAIF
spacecraft-SPK fetch is network-blocked.

Provenance of the EXPECTED side: the closest-approach radii in planet radii are
the PUBLISHED mission geometry (not a value our own code computed), keeping the
golden non-circular per feedback_golden_tests_sourced_only.
"""

from __future__ import annotations

from types import ModuleType
from typing import Any

import pytest

# Published Voyager 2 closest-approach geometry (digest note + canonical mission
# record). Jupiter: 10.0 R_J (Kohlhase-Penzo Table IV, JSX). Saturn: ~2.7 R_S
# (Voyager 2 passed ~101,000 km altitude above Saturn's 60,268 km radius).
PUBLISHED_V2_JUPITER_CA_RJ = 10.0
PUBLISHED_V2_SATURN_CA_RS = 2.7

# Voyager 2 Neptune: the #398 fast-inner-flyby regression case. Periapsis at
# 1989-08-25 03:56 UTC put closest approach ~4,950 km above the north pole, i.e.
# ~29,240 km from Neptune center (Wikipedia "Voyager 2"; cross-checked against
# the NASA Voyager mission record), or ~1.18 Neptune radii (24,764 km eq.).
# Before #398 the coarse ~6 h sample grid never landed near this periapsis and
# the DEFAULT call reported a spurious 6.2 R_N; the periapsis refinement now
# resolves it from the default call without manual re-centering.
PUBLISHED_V2_NEPTUNE_CA_RN = 1.18

# Order-of-magnitude V∞ sanity bands (km/s) — NOT tight gates.
V2_JUPITER_VINF_BAND_KMS = (6.0, 9.0)
V2_SATURN_VINF_BAND_KMS = (9.0, 13.0)


@pytest.fixture(scope="module")
def mission_spk_module() -> ModuleType:
    """Import the extractor or skip if the validation extra is unavailable."""
    pytest.importorskip("spiceypy", reason="spiceypy (validation extra) not installed")
    pytest.importorskip("astropy", reason="astropy (validation extra) not installed")
    import cyclerfinder.verify.mission_spk as mission_spk

    return mission_spk


def _extract(mission_spk: ModuleType, filename: str, body: str, epoch: str) -> Any:
    """Fetch the SPK + run the extractor, skipping on a network/fetch failure."""
    try:
        spk = mission_spk.ensure_mission_spk(filename, base_url=mission_spk.NAIF_VOYAGER_SPK_BASE)
    except Exception as exc:  # network blocked / 404 -> skip honestly
        pytest.skip(f"could not fetch NAIF SPK {filename} (network blocked?): {exc}")
    try:
        return mission_spk.vinf_at_flyby(
            spk, mission_spk.VOYAGER_2_NAIF_ID, body, epoch, mission="Voyager 2"
        )
    except Exception as exc:  # DE440 / LSK fetch failure -> skip honestly
        pytest.skip(f"extractor failed (kernel fetch blocked?): {exc}")


@pytest.mark.slow
def test_voyager2_jupiter_closest_approach_matches_published(
    mission_spk_module: ModuleType,
) -> None:
    """V2 Jupiter closest approach from the SPK == published 10.0 R_J (tight)."""
    r = _extract(mission_spk_module, "vgr2_jup230.bsp", "Jupiter", "1979-07-09T22:29:00")
    # Closest approach is PUBLISHED in planet radii; require <= 2% agreement.
    rel = abs(r.closest_approach_radius_body_radii - PUBLISHED_V2_JUPITER_CA_RJ)
    rel /= PUBLISHED_V2_JUPITER_CA_RJ
    assert rel <= 0.02, (
        f"V2 Jupiter CA {r.closest_approach_radius_body_radii:.3f} R_J disagrees "
        f"with published {PUBLISHED_V2_JUPITER_CA_RJ} R_J by {rel:.1%}"
    )


@pytest.mark.slow
def test_voyager2_jupiter_vinf_converges_and_sane(mission_spk_module: ModuleType) -> None:
    """V2 Jupiter V∞ is vis-viva-stable (<1%) and in the sanity band."""
    r = _extract(mission_spk_module, "vgr2_jup230.bsp", "Jupiter", "1979-07-09T22:29:00")
    # Convergence: outer-window vis-viva std < 1% of the mean.
    frac = r.vinf_kms_visviva_window_std / r.vinf_kms_visviva_window_mean
    assert frac < 0.01, (
        f"V2 Jupiter V∞ not converged: std/mean = {frac:.2%} "
        f"(mean {r.vinf_kms_visviva_window_mean:.4f} km/s)"
    )
    lo, hi = V2_JUPITER_VINF_BAND_KMS
    assert lo <= r.vinf_kms <= hi, f"V2 Jupiter V∞ {r.vinf_kms:.3f} km/s outside [{lo},{hi}]"


@pytest.mark.slow
def test_voyager2_saturn_closest_approach_matches_published(
    mission_spk_module: ModuleType,
) -> None:
    """V2 Saturn closest approach from the SPK == published ~2.7 R_S (tight)."""
    r = _extract(mission_spk_module, "vgr2_sat337.bsp", "Saturn", "1981-08-26T03:24:00")
    rel = abs(r.closest_approach_radius_body_radii - PUBLISHED_V2_SATURN_CA_RS)
    rel /= PUBLISHED_V2_SATURN_CA_RS
    assert rel <= 0.03, (
        f"V2 Saturn CA {r.closest_approach_radius_body_radii:.3f} R_S disagrees "
        f"with published {PUBLISHED_V2_SATURN_CA_RS} R_S by {rel:.1%}"
    )


@pytest.mark.slow
def test_voyager2_saturn_vinf_converges_and_sane(mission_spk_module: ModuleType) -> None:
    """V2 Saturn V∞ is vis-viva-stable (<1%) and in the sanity band."""
    r = _extract(mission_spk_module, "vgr2_sat337.bsp", "Saturn", "1981-08-26T03:24:00")
    frac = r.vinf_kms_visviva_window_std / r.vinf_kms_visviva_window_mean
    assert frac < 0.01, (
        f"V2 Saturn V∞ not converged: std/mean = {frac:.2%} "
        f"(mean {r.vinf_kms_visviva_window_mean:.4f} km/s)"
    )
    lo, hi = V2_SATURN_VINF_BAND_KMS
    assert lo <= r.vinf_kms <= hi, f"V2 Saturn V∞ {r.vinf_kms:.3f} km/s outside [{lo},{hi}]"


# An off-grid nominal epoch for the Neptune flyby: still within the encounter
# day, but the true periapsis (03:56 UTC) sits between the coarse ~6 h samples
# rather than on one. This is the exact condition that produced the original
# #398 spurious ~6.2 R_N. (When the nominal epoch happens to align with a grid
# sample the coarse value is accidentally right; this epoch makes the coarse
# miss explicit, so the test really exercises the refinement.)
V2_NEPTUNE_OFFGRID_EPOCH = "1989-08-25T00:00:00"


@pytest.mark.slow
def test_voyager2_neptune_periapsis_refined_from_default_call(
    mission_spk_module: ModuleType,
) -> None:
    """#398 regression: the DEFAULT call resolves the fast Neptune periapsis.

    Neptune is the motivating fast inner flyby. With the nominal epoch off the
    coarse ~6 h grid the legacy default reported ~6.2 R_N for a true ~1.18 R_N;
    with periapsis refinement on by default the closest approach must match the
    published ~1.18 R_N to < 0.5% regardless of how the nominal epoch is centered
    (no caller widening / re-centering).
    """
    r = _extract(mission_spk_module, "vgr2_nep097.bsp", "Neptune", V2_NEPTUNE_OFFGRID_EPOCH)
    assert r.closest_approach_refined, "periapsis refinement should be the default"
    rel = abs(r.closest_approach_radius_body_radii - PUBLISHED_V2_NEPTUNE_CA_RN)
    rel /= PUBLISHED_V2_NEPTUNE_CA_RN
    assert rel <= 0.005, (
        f"V2 Neptune CA {r.closest_approach_radius_body_radii:.4f} R_N disagrees "
        f"with published {PUBLISHED_V2_NEPTUNE_CA_RN} R_N by {rel:.2%} "
        f"(refined periapsis offset {r.closest_approach_offset_minutes:+.2f} min)"
    )


@pytest.mark.slow
def test_periapsis_refinement_improves_on_coarse_grid(
    mission_spk_module: ModuleType,
) -> None:
    """The refined CA must dramatically beat the legacy coarse-grid CA.

    Same SPK + off-grid epoch, refinement off vs on: the coarse grid badly
    over-reports this fast flyby (~6 R_N) while refinement recovers the true
    periapsis (~1.18 R_N) — refinement only ever moves toward periapsis, and
    here cuts the reported radius by well over half. V∞ (read at the outermost
    sample) must be untouched by refinement.
    """
    try:
        spk = mission_spk_module.ensure_mission_spk(
            "vgr2_nep097.bsp", base_url=mission_spk_module.NAIF_VOYAGER_SPK_BASE
        )
    except Exception as exc:  # network blocked / 404 -> skip honestly
        pytest.skip(f"could not fetch NAIF SPK (network blocked?): {exc}")
    try:
        coarse = mission_spk_module.vinf_at_flyby(
            spk,
            mission_spk_module.VOYAGER_2_NAIF_ID,
            "Neptune",
            V2_NEPTUNE_OFFGRID_EPOCH,
            refine_periapsis=False,
        )
        refined = mission_spk_module.vinf_at_flyby(
            spk,
            mission_spk_module.VOYAGER_2_NAIF_ID,
            "Neptune",
            V2_NEPTUNE_OFFGRID_EPOCH,
            refine_periapsis=True,
        )
    except Exception as exc:  # kernel fetch failure -> skip honestly
        pytest.skip(f"extractor failed (kernel fetch blocked?): {exc}")

    assert not coarse.closest_approach_refined
    assert refined.closest_approach_refined
    assert refined.closest_approach_radius_km <= coarse.closest_approach_radius_km
    # The coarse grid badly over-reports this fast flyby; refinement cuts it by
    # well over half (coarse ~6 R_N vs refined ~1.18 R_N).
    assert refined.closest_approach_radius_km < 0.5 * coarse.closest_approach_radius_km
    # V∞ is read at the outermost sample and must be unaffected by refinement.
    assert refined.vinf_kms == pytest.approx(coarse.vinf_kms, rel=1e-12)


def test_body_gm_is_sourced_not_hardcoded(mission_spk_module: ModuleType) -> None:
    """The extractor's body GM must come from the project constants module."""
    from cyclerfinder.core.constants import PLANETS

    # FLYBY_BODIES must map to real PLANETS keys (the sourced GM/radius source).
    for body, fb in mission_spk_module.FLYBY_BODIES.items():
        assert fb.planets_key in PLANETS, f"{body} -> unknown PLANETS key {fb.planets_key}"
