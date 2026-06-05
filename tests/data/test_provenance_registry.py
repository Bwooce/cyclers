"""Task 1: source + fidelity registry tests.

Asserts the registry vocabulary is well-formed: unique keys, non-empty
citations, the two pseudo-sources present, the three fidelity tiers
recognised, and the pseudo-source / independence distinction.
"""

from __future__ import annotations

import pytest

from cyclerfinder.data.provenance import (
    SOURCE_REGISTRY,
    is_fidelity,
    is_independent_source,
    is_registry_key,
)

_EXPECTED_KEYS = {
    "rogers-2012-t1",
    "russell-2004-t34",
    "russell-2004-t39_311",
    "russell-2004-t49_413",
    "mcconaghy-2002",
    "mcconaghy-2006",
    "spec-9",
    "hollister-1970-t3",
    "friedlander-1986",
    "derived",
    "computed",
}


def test_registry_contains_all_planned_keys() -> None:
    """Every key the plan (Task 1) names is present."""
    missing = _EXPECTED_KEYS - set(SOURCE_REGISTRY)
    assert not missing, f"registry missing planned keys: {sorted(missing)}"


def test_registry_keys_unique_and_citations_nonempty() -> None:
    """Citations are all non-empty strings (dict keys are unique by type)."""
    for key, citation in SOURCE_REGISTRY.items():
        assert isinstance(citation, str), f"{key}: citation not a str"
        assert citation.strip(), f"{key}: citation is empty"


def test_pseudo_sources_present() -> None:
    """The two pseudo-sources are in the registry."""
    assert "derived" in SOURCE_REGISTRY
    assert "computed" in SOURCE_REGISTRY


def test_pseudo_sources_are_not_independent() -> None:
    """derived / computed are registry keys but NOT independent citations."""
    assert is_registry_key("derived")
    assert is_registry_key("computed")
    assert not is_independent_source("derived")
    assert not is_independent_source("computed")


def test_real_citations_are_independent() -> None:
    """A real citation key is an independent source."""
    assert is_independent_source("rogers-2012-t1")
    assert is_independent_source("russell-2004-t34")


def test_unknown_key_is_not_registry_or_independent() -> None:
    """An unknown key is neither a registry key nor independent."""
    assert not is_registry_key("not-a-real-source")
    assert not is_independent_source("not-a-real-source")


@pytest.mark.parametrize(
    "value",
    ["circular-coplanar", "analytic-ephemeris", "real-de440"],
)
def test_fidelity_tiers_recognised(value: str) -> None:
    """Each of the three fidelity tiers is recognised."""
    assert is_fidelity(value)


@pytest.mark.parametrize("value", ["", "de440", "coplanar", "real"])
def test_non_fidelity_strings_rejected(value: str) -> None:
    """Strings that are not exact fidelity tiers are rejected."""
    assert not is_fidelity(value)
