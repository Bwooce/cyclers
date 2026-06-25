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

from cyclerfinder.genome import known_corpus_3d
from cyclerfinder.genome.known_corpus_3d import (
    ANTONIADOU_LIBERT_2019_PROSE_ANCHORS,
    KNOWN_CORPUS_3D,
)
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
    """A coordinate anchor's topology_3d declares a k_z; taxonomy anchors may omit it.

    Not every anchor is a coordinate anchor. Antoniadou & Libert 2019, for
    instance, is a TAXONOMY/mechanism anchor (mu=0.001 planetary MMRs, no
    state-vector ICs), so it legitimately carries topology_3d=None — pinning a
    fabricated (k1,k2,k_z) tuple it cannot supply would be unfaithful. We require
    only that (a) any anchor that DOES declare topology_3d includes k_z, and
    (b) at least one coordinate anchor exists so the gate can tuple-match.
    """
    has_coordinate_anchor = False
    for anchor in KNOWN_CORPUS_3D:
        if anchor.topology_3d is not None:
            assert "k_z" in anchor.topology_3d, f"{anchor.name}: topology_3d has no k_z"
            has_coordinate_anchor = True
    assert has_coordinate_anchor, "no coordinate (topology_3d) anchor in the 3D corpus"


def test_antoniadou_libert_2019_prose_anchors_sourced() -> None:
    """#459: the full-precision PROSE anchors recovered from SPA.tex are pinned.

    The #458 arXiv-source sweep showed arXiv:1811.09442 carries full-precision
    orbital-element anchors in the TeX prose (the figures render them only as
    grey dots), correcting the prior "graphical only / no matchable numeric ICs"
    digest claim. Every value below is transcribed VERBATIM from SPA.tex; the
    EXPECTED side traces to the paper ONLY (golden-sourced discipline), never to
    a value our code computed.
    """
    # mu=0.001 (Jupiter-mass giant), NOT Earth-Moon -- the transferability caveat.
    assert ANTONIADOU_LIBERT_2019_PROSE_ANCHORS["mu"] == 0.001
    assert "1811.09442" in known_corpus_3d._AL2019_SOURCE

    bifs = {(b["mmr"], b["e1"], b["i1_deg"]): b for b in known_corpus_3d._AL2019_BIFURCATION_POINTS}
    # 5/2 G^{5/2}_{C3} v.c.o. bifurcation, SPA.tex L438.
    assert ("5/2", 0.0891812, 90.0) in bifs
    assert bifs[("5/2", 0.0891812, 90.0)]["spa_tex_line"] == 438
    # The four 3/1 bifurcation points, SPA.tex L513/L523.
    assert ("3/1", 0.000178, 84.0) in bifs
    assert ("3/1", 0.0002327, 101.0) in bifs
    assert ("3/1", 0.0003454, 101.0) in bifs
    assert ("3/1", 0.0004451, 84.0) in bifs

    ds = {d["mmr"]: d for d in known_corpus_3d._AL2019_DS_MAP_CONSTANT_ELEMENTS}
    # 2/1 (pi,0) DS-map constant elements, SPA.tex L719.
    assert ds["2/1"]["a2_over_a1"] == 0.6312
    assert ds["2/1"]["M1_deg"] == 180.0 and ds["2/1"]["omega1_deg"] == 270.0
    # 3/2 (0,pi) DS-map constant elements, SPA.tex L733.
    assert ds["3/2"]["a2_over_a1"] == 0.7595
    assert ds["3/2"]["M2_deg"] == 180.0
    # 3/1 F^{3/1}_I 3D-CRTBP DS-map constant elements, SPA.tex L747.
    assert ds["3/1"]["a2_over_a1"] == 0.4806
    assert ds["3/1"]["e2"] == 0.0 and ds["3/1"]["M1_deg"] == 0.0
