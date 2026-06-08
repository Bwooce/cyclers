"""Tier-1 Phase 6: matcher pool pre-filters by (model_assumption, primary)
(plan Phase 6; spec §16.2 / docs/spec.md). Heliocentric rows (primary
absent->Sun) never share a bucket with Jovicentric rows."""

from __future__ import annotations

import yaml  # type: ignore[import-untyped]

from cyclerfinder.data.catalog import signature_bucket_key
from tests._catalogue_loader import CATALOGUE_PATH


def test_missing_primary_is_sun_bucket() -> None:
    assert signature_bucket_key({"model_assumption": "circular-coplanar"})[1] == "Sun"


def test_jovian_and_heliocentric_are_different_buckets() -> None:
    helio = signature_bucket_key({"model_assumption": "circular-coplanar"})
    jovian = signature_bucket_key({"model_assumption": "circular-coplanar", "primary": "Jupiter"})
    assert helio != jovian


def test_all_heliocentric_rows_bucket_as_sun() -> None:
    for row in yaml.safe_load(CATALOGUE_PATH.read_text()):
        if row.get("primary") in (None, "Sun"):
            assert signature_bucket_key(row)[1] == "Sun"
