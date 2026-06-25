"""Fast up-front novelty pre-filter for 3D (out-of-plane) periodic orbits.

Motivation (#444). The C21 stable 3D candidate cost a full gauntlet
(Floquet + V1 + V2 + a JPL full-family sweep) BEFORE a literature check closed
it as a member of a *published* class. The expensive verification ran first and
the cheap, decisive check ran last. This module is the cheap check, meant to run
FIRST: it flags 3D candidates that are members of an already-published class so a
campaign can route them to the reproduction ledger instead of burning a
discovery gauntlet on them.

The sourced fact it keys on (the vertical-bifurcation MECHANISM)
----------------------------------------------------------------
Antoniadou & Libert (2019, MNRAS 483(3):2923, DOI 10.1093/mnras/sty3195;
already a :data:`cyclerfinder.genome.known_corpus_3d.KNOWN_CORPUS_3D` anchor)
and the Antoniadou-Voyatzis lineage establish that **spatial (out-of-plane,
k_z>0) periodic families emanate from the vertical-critical orbits of planar
periodic families.** Lifting a *known* planar periodic cycler out of plane
therefore lands, by construction, on its spatial-bifurcation family — which is
the published mechanism, not new structure.

So the pre-filter keys on a fact the campaign already knows — *was the planar
root a known periodic family?* — NOT on a fabricated ``(k1,k2,k_z)`` fingerprint
the source cannot supply (the Antoniadou anchor is principled-``None`` for
exactly this reason; see its note in ``known_corpus_3d``). This keeps the check
honest: it asserts only the published *mechanism*, never an unsourced winding.

What it is / is NOT
-------------------
* It IS a novelty *downgrade* flag: "treat this as published-class; do not claim
  a discovery without showing the specific spatial family is absent from the
  literature for this system." Necessary-not-sufficient, like every literature
  gate (a specific spatial member could still be genuinely new — e.g. an
  ASYMMETRIC or ISOLATED family with no planar/vertical-critical limit; those
  are exactly what Antoniadou & Libert do NOT cover, so they stay novelty-open).
* It is NOT a closure/stability check and NOT a catalogue gate. It changes
  *routing* (reproduction vs discovery), saving the gauntlet cost.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

# The published mechanism's canonical citation (matches the KNOWN_CORPUS_3D
# Antoniadou & Libert anchor verbatim DOI). Stored as data so a verdict can hand
# the campaign a ready citation for the reproduction ledger.
SPATIAL_BIFURCATION_CITATION: Final[str] = (
    "Antoniadou, K. I. & Libert, A.-S., 'Spatial resonant periodic orbits in "
    "the restricted three-body problem,' MNRAS 483(3):2923-2940 (2019); DOI "
    "10.1093/mnras/sty3195; arXiv:1811.09442 (with the Antoniadou-Voyatzis "
    "2013/2014 lineage): spatial (k_z>0) families emanate from the "
    "vertical-critical orbits of planar periodic families."
)
SPATIAL_BIFURCATION_DOI: Final[str] = "10.1093/mnras/sty3195"


@dataclass(frozen=True)
class SpatialNoveltyVerdict:
    """Verdict of the fast spatial-novelty pre-filter.

    Attributes
    ----------
    published_mechanism:
        ``True`` iff the candidate is the out-of-plane extension of a KNOWN
        planar periodic family — i.e. a member of the published
        vertical-bifurcation mechanism. When ``True`` the campaign should route
        the candidate to the reproduction ledger and NOT spend a discovery
        gauntlet claiming novelty.
    is_out_of_plane:
        ``True`` iff ``k_z > 0`` (genuinely 3D). A planar (k_z==0) candidate is
        never flagged by this mechanism (the planar root itself, not a spatial
        bifurcation).
    citation:
        The published-mechanism citation when ``published_mechanism`` is True;
        empty string otherwise.
    doi:
        DOI of that citation, or empty string.
    reason:
        Human-readable explanation for the audit trail.
    """

    published_mechanism: bool
    is_out_of_plane: bool
    citation: str
    doi: str
    reason: str


def classify_spatial_extension(
    *,
    k_z: int,
    planar_root_is_known_periodic: bool,
) -> SpatialNoveltyVerdict:
    """Flag whether a 3D candidate is a published spatial-bifurcation member.

    Run this BEFORE the V0-V5 gauntlet on any out-of-plane candidate produced by
    lifting a planar family into ``z != 0``.

    Parameters
    ----------
    k_z:
        Equatorial-plane crossing count (the out-of-plane winding component).
        ``k_z == 0`` is a planar member; ``k_z > 0`` is genuinely 3D.
    planar_root_is_known_periodic:
        Whether the planar orbit this candidate was lifted from is a KNOWN
        periodic family (e.g. a catalogued cycler row, or any published planar
        periodic family). The caller knows this: a 3D broken-plane campaign that
        lifts ``ross-rt-em-cycler-21`` passes ``True``; a from-scratch direct 3D
        search with no known planar parent passes ``False``.

    Returns
    -------
    SpatialNoveltyVerdict
        ``published_mechanism=True`` iff ``k_z > 0`` AND the planar root is a
        known periodic family.

    Notes
    -----
    A ``published_mechanism=True`` verdict is NECESSARY-NOT-SUFFICIENT for
    "published": it asserts the candidate's *class* (a spatial bifurcation of a
    known planar family) is published, which is enough to downgrade novelty and
    route to reproduction. A genuinely novel 3D orbit must therefore have NO
    known planar periodic root (asymmetric / isolated spatial families that do
    not emanate from a vertical-critical orbit) — exactly the frontier
    Antoniadou & Libert leave open.
    """
    is_out_of_plane = int(k_z) > 0
    if is_out_of_plane and planar_root_is_known_periodic:
        return SpatialNoveltyVerdict(
            published_mechanism=True,
            is_out_of_plane=True,
            citation=SPATIAL_BIFURCATION_CITATION,
            doi=SPATIAL_BIFURCATION_DOI,
            reason=(
                "Out-of-plane (k_z>0) extension of a KNOWN planar periodic "
                "family: a member of the published vertical-bifurcation "
                "mechanism (spatial families emanate from vertical-critical "
                "orbits of planar families). Route to the reproduction ledger; "
                "do not claim a discovery without showing the specific spatial "
                "family is absent from the literature for this system."
            ),
        )
    if not is_out_of_plane:
        reason = (
            "Planar candidate (k_z==0): not a spatial-bifurcation member; this gate does not apply."
        )
    else:
        reason = (
            "Out-of-plane (k_z>0) but NO known planar periodic root: the spatial "
            "bifurcation mechanism does not cover it (no vertical-critical-orbit "
            "parent). Novelty stays OPEN — proceed to the gauntlet + literature check."
        )
    return SpatialNoveltyVerdict(
        published_mechanism=False,
        is_out_of_plane=is_out_of_plane,
        citation="",
        doi="",
        reason=reason,
    )


__all__ = [
    "SPATIAL_BIFURCATION_CITATION",
    "SPATIAL_BIFURCATION_DOI",
    "SpatialNoveltyVerdict",
    "classify_spatial_extension",
]
