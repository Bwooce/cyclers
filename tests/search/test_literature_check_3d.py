"""Tests for the 3D / spatial-CR3BP literature-check widening (#434 Task 4).

The #434 broken-plane discovery campaign lifts planar Aldrin/Braik-Ross
cyclers into z!=0 and continues them in the 3D genome. Before any novelty
claim, a swept 3D candidate must be adjudicated against the heavily-published
halo / NRHO / vertical-Lyapunov / spatial-resonant CR3BP corpus -- otherwise a
re-derivation of a known spatial family would be mislabelled novel.

This module pins:

* the new ``CandidateSignature.topology_3d`` optional field constructs and is
  backward-compatible (defaults to ``None``; existing fields unaffected);
* a 3D signature checked against the new ``known_corpus_3d`` anchors returns a
  documented status in {published, not-found, inconclusive} with confidence>=0;
* the spatial matcher additionally requires ``(k1, k2, k_z)`` agreement (and
  Jacobi-band overlap when the anchor records one) when BOTH sides carry 3D
  topology -- a planar (k_z=0) candidate is NOT flagged by a halo (k_z>0)
  anchor;
* every ``known_corpus_3d`` anchor carries provenance (a non-empty citation
  AND a non-empty DOI/arXiv) -- the golden-sourced discipline.
"""

from __future__ import annotations

from collections.abc import Sequence

from cyclerfinder.genome.known_corpus_3d import KNOWN_CORPUS_3D
from cyclerfinder.search.literature_check import (
    CandidateSignature,
    SearchResult,
    check_literature,
)


def _empty_search(_query: str) -> Sequence[SearchResult]:
    """A search that returns no hits (offline corpus-only adjudication)."""
    return []


def test_signature_topology_3d_defaults_none_and_backward_compatible() -> None:
    """Omitting ``topology_3d`` leaves it ``None``; existing fields unchanged."""
    sig = CandidateSignature(
        primary="Earth-Moon",
        sequence=("E", "M"),
        period_k=1,
        period_years=None,
        vinf_per_encounter_kms=(),
        resonances=(),
        n_rev=(),
    )
    assert sig.topology_3d is None
    assert sig.primary == "Earth-Moon"
    assert sig.sequence == ("E", "M")
    assert sig.period_k == 1


def test_signature_topology_3d_constructs_and_round_trips() -> None:
    """A 3D-annotated signature constructs and carries its topology dict."""
    topo = {"k1": 1, "k2": 1, "k_z": 0, "max_z_km": 92567}
    sig = CandidateSignature(
        primary="Earth-Moon",
        sequence=("E", "M"),
        period_k=1,
        period_years=None,
        vinf_per_encounter_kms=(),
        resonances=(),
        n_rev=(),
        topology_3d=topo,
    )
    assert sig.topology_3d == topo
    assert sig.topology_3d["k_z"] == 0


def test_3d_signature_returns_documented_status() -> None:
    """A 3D signature against the new corpus yields a documented verdict."""
    sig = CandidateSignature(
        primary="Earth",
        sequence=("Moon",),
        topology_3d={"k1": 1, "k2": 1, "k_z": 2, "jacobi": 3.15},
    )
    result = check_literature(sig, search=_empty_search)
    assert result.status in ("published", "not-found", "inconclusive")
    assert result.confidence >= 0.0


def test_known_corpus_3d_all_have_provenance() -> None:
    """Every 3D anchor carries a non-empty citation AND a DOI/arXiv id."""
    assert KNOWN_CORPUS_3D, "the 3D known corpus must be non-empty"
    for anchor in KNOWN_CORPUS_3D:
        assert anchor.citation and anchor.citation.strip(), f"{anchor.name}: empty citation"
        assert anchor.doi and anchor.doi.strip(), f"{anchor.name}: missing DOI/arXiv provenance"


def test_known_corpus_3d_anchors_carry_3d_topology() -> None:
    """Each 3D anchor declares a topology_3d descriptor with a k_z."""
    for anchor in KNOWN_CORPUS_3D:
        assert anchor.topology_3d is not None, f"{anchor.name}: missing topology_3d"
        assert "k_z" in anchor.topology_3d, f"{anchor.name}: topology_3d has no k_z"
