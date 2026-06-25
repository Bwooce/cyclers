"""Sourced spatial-CR3BP family anchors for the 3D literature-check (#434).

The #434 broken-plane discovery campaign lifts planar Aldrin/Braik-Ross
cyclers into ``z != 0`` and continues them in the 3D genome. The out-of-plane
families it can rediscover -- halos, NRHOs, vertical-Lyapunov / butterfly,
spatial-resonant orbits -- are exhaustively published. This module pins those
published spatial-CR3BP families as :class:`~cyclerfinder.search.literature_check.CorpusAnchor`
records so a swept 3D candidate's structural fingerprint can be adjudicated
against the corpus BEFORE any novelty claim (the literature gate is mandatory).

DISCIPLINE (golden-sourced):

* **Every anchor carries provenance.** No anchor exists without a real
  ``citation`` AND a non-empty ``doi`` (DOI or arXiv id). This is enforced by
  ``tests/search/test_literature_check_3d.py``.
* **No fabricated ICs.** Anchors carry citation + topology (``(k1, k2, k_z)``)
  + an optional Jacobi band ONLY. We do NOT invent published initial conditions
  or numerical state vectors here -- a 3D candidate is matched on its
  ``(k1, k2, k_z)`` winding fingerprint and Jacobi level, not on a fabricated
  exact orbit.
* **k_z is the out-of-plane signature.** Vertical-Lyapunov / halo / NRHO /
  butterfly families are ``k_z > 0`` (they cross the equatorial plane). The
  planar (1,1) Braik-Ross root is ``k_z = 0``; its #287 3D extension is the
  novel-frontier target the campaign is testing for, NOT one of these anchors.

The matcher (``literature_check._candidate_anchors`` / ``_spatial_topology_matches``)
folds these anchors in only when the candidate signature carries a
``topology_3d`` dict, and requires the ``(k1, k2, k_z)`` tuple to agree (+
Jacobi-band overlap when recorded). A planar candidate never touches this
corpus, preserving the historical planar literature-check behaviour exactly.

Jacobi bands: the Earth-Moon halo/NRHO families span roughly C ~ 2.9-3.18 in
the standard mu = 0.01215 normalisation (L1/L2 collinear C ~ 3.0-3.2; deep
NRHOs run lower). Where a sourced band is not pinned to a single published
number it is left ``None`` and the ``(k1, k2, k_z)`` tuple alone decides;
recorded bands are deliberately generous (family-extent, not a single member).
"""

from __future__ import annotations

from typing import Final

from cyclerfinder.search.literature_check import CorpusAnchor

# ---------------------------------------------------------------------------
# Spatial-CR3BP known corpus.
#
# These are PUBLICATION facts (author / venue / DOI / arXiv) + the published
# family's out-of-plane topology class, NOT values our code computed.
# ---------------------------------------------------------------------------
KNOWN_CORPUS_3D: tuple[CorpusAnchor, ...] = (
    CorpusAnchor(
        # Howell 1984 -- the foundational three-dimensional periodic halo
        # family in the CR3BP (collinear L1/L2/L3). The defining halo
        # reference; every Earth-Moon L1/L2 halo a 3D lift could land on
        # traces to this continuation family. Halos cross the equatorial
        # plane (vertical structure) -> k_z > 0. The canonical Earth-Moon
        # halo continuation spans roughly C ~ 3.0-3.18 (collinear-point
        # energy band); recorded as a generous family-extent band.
        name="Howell three-dimensional periodic halo families (1984)",
        primary="Earth",
        body_set=frozenset({"Moon"}),
        topology_label=frozenset({"halo"}),
        topology_3d={"k1": 1, "k2": 0, "k_z": 2},
        jacobi_band=(2.90, 3.20),
        authors=("Howell",),
        keywords=(
            "three-dimensional periodic halo orbit CR3BP",
            "halo family collinear libration point",
            "vertical out-of-plane periodic orbit Earth-Moon",
        ),
        citation=(
            "Howell, K. C., 'Three-dimensional, periodic, 'halo' orbits,' "
            "Celestial Mechanics 32(1):53-71 (1984), DOI 10.1007/BF01358403. "
            "The foundational CR3BP halo continuation family at the collinear "
            "libration points; the canonical reference for any Earth-Moon "
            "L1/L2 halo a 3D broken-plane lift could rediscover."
        ),
        doi="10.1007/BF01358403",
    ),
    CorpusAnchor(
        # Folta-Bosanac-Guzzetti-Howell 2015 -- the Earth-Moon trajectory
        # design reference catalog: the published cislunar L1/L2 halo + NRHO
        # family catalogue (the Gateway-era L2 Southern NRHO lineage). NRHOs
        # are the deep, near-rectilinear tail of the L1/L2 halo family; still
        # k_z > 0 (vertical/out-of-plane). Deep NRHOs reach lower Jacobi than
        # the shallow halos, so the recorded band runs lower.
        #
        # NOTE: the #434 task brief named this anchor "Folta-Bosanac-Cox-Howell
        # 2017"; web verification (2026-06-24) shows the cislunar halo/NRHO
        # reference catalog is Folta, Bosanac, GUZZETTI & Howell, Acta
        # Astronautica 110 (2015), DOI 10.1016/j.actaastro.2014.07.037 -- Cox
        # is a co-author on the SEPARATE Bosanac-Cox-Howell-Folta Lunar IceCube
        # 2018 paper, not this catalogue. Author tuple corrected to the
        # verified roster; DOI confirmed against ScienceDirect / CU Experts.
        name="Folta-Bosanac-Guzzetti-Howell Earth-Moon halo/NRHO reference catalog (2015)",
        primary="Earth",
        body_set=frozenset({"Moon"}),
        topology_label=frozenset({"halo", "nrho"}),
        topology_3d={"k1": 1, "k2": 0, "k_z": 2},
        jacobi_band=(2.95, 3.18),
        authors=("Folta", "Bosanac", "Guzzetti", "Howell"),
        keywords=(
            "cislunar L1 L2 halo orbit access",
            "near-rectilinear halo orbit NRHO",
            "Earth-Moon libration point orbit transfer design",
        ),
        citation=(
            "Folta, D. C., Bosanac, N., Guzzetti, D. & Howell, K. C., 'An "
            "Earth-Moon system trajectory design reference catalog,' Acta "
            "Astronautica 110:341-353 (2015), DOI "
            "10.1016/j.actaastro.2014.07.037 -- the cislunar L1/L2 halo + "
            "near-rectilinear halo (NRHO) reference family; the published "
            "spatial-CR3BP corpus any Earth-Moon halo/NRHO 3D lift must be "
            "compared against. (The #434 brief's 'Folta-Bosanac-Cox-Howell "
            "2017' label conflated this catalogue with the Bosanac-Cox-Howell-"
            "Folta Lunar IceCube 2018 paper; roster + DOI corrected here.)"
        ),
        doi="10.1016/j.actaastro.2014.07.037",
    ),
    CorpusAnchor(
        # Antoniadou & Libert 2019 -- spatial resonant periodic orbits in the
        # RTBP (digest 2026-06-25; mu=0.001 Jupiter-mass planetary MMRs 3/2, 2/1,
        # 5/2, 3/1, 4/1, 5/1, NOT Earth-Moon). TAXONOMY/MECHANISM anchor for the
        # MATCHER (topology_3d=None): its matcher value is the bifurcation
        # MECHANISM -- spatial families (k_z>0) emanate from vertical-critical-
        # orbits of planar families -- not a Cartesian (k1,k2,k_z) Earth-Moon
        # tuple this mu cannot supply. It does NOT cover asymmetric/spatial-
        # isolated families, so 3D novelty stays open for those.
        #
        # #459 source-recovery (from the #458 arXiv-source sweep): the EARLIER
        # claim "tabulates NO state-vector ICs" was over-broad. The paper carries
        # NO Cartesian (x0,z0,ydot0) IC table, but the SPA.tex PROSE DOES carry
        # full-precision, reproducible ORBITAL-ELEMENT anchors that the figures
        # only render as grey dots -- transcribed verbatim below in
        # ANTONIADOU_LIBERT_2019_PROSE_ANCHORS with their SPA.tex line numbers.
        # These are mu=0.001 exoplanetary numbers (NOT Earth-Moon coordinates),
        # so they stay out of the (k1,k2,k_z) matcher (topology_3d remains None);
        # they are pinned as sourced taxonomy data so the corpus no longer
        # records a false "graphical-only" gap.
        name="Antoniadou & Libert spatial resonant periodic orbits in the RTBP (2019)",
        primary="Earth",
        body_set=frozenset({"Moon"}),
        topology_label=frozenset({"resonant"}),
        topology_3d=None,
        jacobi_band=None,
        authors=("Antoniadou", "Libert"),
        keywords=(
            "spatial resonant periodic orbit restricted three-body problem",
            "vertical critical orbit spatial resonance bifurcation",
            "out-of-plane resonant family taxonomy",
        ),
        citation=(
            "Antoniadou, K. I. & Libert, A.-S., 'Spatial resonant periodic "
            "orbits in the restricted three-body problem,' MNRAS 483(3):"
            "2923-2940 (2019); DOI 10.1093/mnras/sty3195; arXiv:1811.09442. "
            "Taxonomy/mechanism anchor (mu=0.001 planetary MMRs; spatial "
            "families born at vertical-critical-orbits) -- supplies no "
            "Earth-Moon coordinate ICs; does not cover asymmetric/isolated "
            "spatial families."
        ),
        doi="10.1093/mnras/sty3195",
    ),
)

# ---------------------------------------------------------------------------
# #459 -- sourced full-precision PROSE anchors from Antoniadou & Libert 2019
# (arXiv:1811.09442, SPA.tex). Recovered by the #458 arXiv-source sweep: these
# numbers live in the TeX prose/captions; the published figures render them only
# as grey/magenta/green dots, which is why the original digest mis-recorded the
# paper as "graphical only / no matchable numeric ICs".
#
# DISCIPLINE: every value below is transcribed VERBATIM from SPA.tex with its
# source line number. These are mu=0.001 (Jupiter-mass-giant) exoplanetary
# orbital elements, NOT Earth-Moon Cartesian state -- they do NOT transfer to an
# Earth-Moon (k1,k2,k_z) tuple and are deliberately kept OUT of the literature
# matcher (the CorpusAnchor above keeps topology_3d=None). They are pinned here
# as sourced taxonomy data so the corpus records that full-precision anchors DO
# exist, replacing the prior false "no matchable numeric ICs" gap.
#
# Semantics:
#   * e1, i1  = inner massless body's eccentricity / inclination (deg).
#   * A "bifurcation point" is a 3D-CRTBP -> 3D-ERTBP bifurcation (Method II) on
#     the named spatial family, given as (e1, i1).
#   * A "DS-map constant element set" is the full set of orbital elements held
#     FIXED while sweeping a dynamical-stability map: a2/a1 (semi-major-axis
#     ratio) plus the six angles per body (deg) and i2. This is a COMPLETE,
#     reproducible element specification (only the swept plane, e.g. (e1,e2) or
#     (e1,i1), is free), not a curve endpoint.
# ---------------------------------------------------------------------------

_AL2019_SOURCE: Final = (
    "Antoniadou & Libert 2019, MNRAS 483(3):2923-2940; "
    "DOI 10.1093/mnras/sty3195; arXiv:1811.09442 (SPA.tex)"
)

# 3D-CRTBP -> 3D-ERTBP bifurcation points, given as (e1, i1_deg) in the prose.
_AL2019_BIFURCATION_POINTS: tuple[dict[str, object], ...] = (
    {
        # SPA.tex L438: G^{5/2}_{C3} x-symmetric family, Method II
        "mmr": "5/2",
        "family": "G^{5/2}_{C3}",
        "e1": 0.0891812,
        "i1_deg": 90.0,
        "spa_tex_line": 438,
        "note": "v.c.o. bifurcation from 3D-CRTBP to 3D-ERTBP on the x-symmetric family",
    },
    # SPA.tex L513 (caption) / L523 (prose): four 3/1 MMR bifurcation points
    # on the spatial families emanating from circular family C1.
    {
        "mmr": "3/1",
        "family": "F_{C1}^{3/1} / G_{C1}^{3/1} / F'_{C1}^{3/1} / G'_{C1}^{3/1}",
        "e1": 0.000178,
        "i1_deg": 84.0,
        "spa_tex_line": 523,
        "note": "first of four 3/1 3D-CRTBP->3D-ERTBP bifurcation points",
    },
    {
        "mmr": "3/1",
        "family": "F_{C1}^{3/1} / G_{C1}^{3/1} / F'_{C1}^{3/1} / G'_{C1}^{3/1}",
        "e1": 0.0002327,
        "i1_deg": 101.0,
        "spa_tex_line": 523,
        "note": "second of four 3/1 3D-CRTBP->3D-ERTBP bifurcation points",
    },
    {
        "mmr": "3/1",
        "family": "F_{C1}^{3/1} / G_{C1}^{3/1} / F'_{C1}^{3/1} / G'_{C1}^{3/1}",
        "e1": 0.0003454,
        "i1_deg": 101.0,
        "spa_tex_line": 523,
        "note": "third of four 3/1 3D-CRTBP->3D-ERTBP bifurcation points",
    },
    {
        "mmr": "3/1",
        "family": "F_{C1}^{3/1} / G_{C1}^{3/1} / F'_{C1}^{3/1} / G'_{C1}^{3/1}",
        "e1": 0.0004451,
        "i1_deg": 84.0,
        "spa_tex_line": 523,
        "note": "fourth of four 3/1 3D-CRTBP->3D-ERTBP bifurcation points",
    },
)

# DS-map constant orbital-element sets (angles in degrees). Each is the full
# FIXED element specification for the named DS-map; the swept plane is noted.
_AL2019_DS_MAP_CONSTANT_ELEMENTS: tuple[dict[str, object], ...] = (
    {
        # SPA.tex L719: 2/1 MMR config (pi,0), 2D-ERTBP DS-map; swept (e1,e2)
        "mmr": "2/1",
        "config": "(pi,0)",
        "swept_plane": "(e1,e2)",
        "a2_over_a1": 0.6312,
        "omega1_deg": 270.0,
        "Omega1_deg": 270.0,
        "M1_deg": 180.0,
        "i2_deg": 0.0,
        "omega2_deg": 270.0,
        "Omega2_deg": 90.0,
        "M2_deg": 0.0,
        "spa_tex_line": 719,
    },
    {
        # SPA.tex L733: 3/2 MMR config (0,pi), 2D-ERTBP DS-map; swept (e1,e2)
        "mmr": "3/2",
        "config": "(0,pi)",
        "swept_plane": "(e1,e2)",
        "a2_over_a1": 0.7595,
        "omega1_deg": 90.0,
        "Omega1_deg": 270.0,
        "M1_deg": 0.0,
        "i2_deg": 0.0,
        "omega2_deg": 90.0,
        "Omega2_deg": 90.0,
        "M2_deg": 180.0,
        "spa_tex_line": 733,
    },
    {
        # SPA.tex L747: 3/1 MMR F^{3/1}_I 3D-CRTBP DS-map; swept (e1,i1).
        # xz-symmetric stable periodic orbit; e2=i2=omega2=Omega2=M2=0.
        "mmr": "3/1",
        "config": "F^{3/1}_I (3D-CRTBP)",
        "swept_plane": "(e1,i1)",
        "a2_over_a1": 0.4806,
        "omega1_deg": 90.0,
        "Omega1_deg": 90.0,
        "M1_deg": 0.0,
        "e2": 0.0,
        "i2_deg": 0.0,
        "omega2_deg": 0.0,
        "Omega2_deg": 0.0,
        "M2_deg": 0.0,
        "spa_tex_line": 747,
    },
)

ANTONIADOU_LIBERT_2019_PROSE_ANCHORS: dict[str, object] = {
    "source": _AL2019_SOURCE,
    "mu": 0.001,  # m2/(m0+m2), Jupiter-mass giant (SPA.tex sec 2.1); NOT Earth-Moon
    "bifurcation_points": _AL2019_BIFURCATION_POINTS,
    "ds_map_constant_elements": _AL2019_DS_MAP_CONSTANT_ELEMENTS,
}

__all__ = ["ANTONIADOU_LIBERT_2019_PROSE_ANCHORS", "KNOWN_CORPUS_3D"]
