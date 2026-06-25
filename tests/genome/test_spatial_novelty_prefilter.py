"""Tests for the fast spatial-novelty pre-filter (#444).

The pre-filter must flag the C21-style case — an out-of-plane extension of a
KNOWN planar catalogued cycler — as published-mechanism BEFORE any gauntlet,
while leaving genuinely-parentless 3D candidates novelty-open.
"""

from __future__ import annotations

from cyclerfinder.genome.spatial_novelty_prefilter import (
    SPATIAL_BIFURCATION_DOI,
    classify_spatial_extension,
)


def test_c21_style_out_of_plane_extension_of_known_cycler_is_published_mechanism() -> None:
    """The exact #444 case: lifting a catalogued planar cycler to z!=0.

    The C21 candidate is k_z=10 and was lifted from ross-rt-em-cycler-21 (a
    catalogued published planar cycler). The pre-filter must flag it
    published-mechanism with the Antoniadou & Libert citation.
    """
    v = classify_spatial_extension(k_z=10, planar_root_is_known_periodic=True)
    assert v.published_mechanism is True
    assert v.is_out_of_plane is True
    assert v.doi == SPATIAL_BIFURCATION_DOI
    assert v.citation  # non-empty, hands the campaign a ready reproduction citation
    assert "vertical-bifurcation" in v.reason or "vertical-critical" in v.reason


def test_planar_candidate_is_not_flagged() -> None:
    """A planar member (k_z==0) is the root itself, not a spatial bifurcation."""
    v = classify_spatial_extension(k_z=0, planar_root_is_known_periodic=True)
    assert v.published_mechanism is False
    assert v.is_out_of_plane is False
    assert v.citation == ""


def test_parentless_3d_candidate_stays_novelty_open() -> None:
    """An out-of-plane candidate with NO known planar root is NOT downgraded.

    This is the genuine novel-3D frontier (asymmetric / isolated spatial
    families with no vertical-critical-orbit parent) — the pre-filter must let
    it through to the full gauntlet + literature check.
    """
    v = classify_spatial_extension(k_z=10, planar_root_is_known_periodic=False)
    assert v.published_mechanism is False
    assert v.is_out_of_plane is True
    assert v.citation == ""
    assert "OPEN" in v.reason
