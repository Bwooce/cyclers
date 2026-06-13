"""SOURCED CR3BP goldens from Szebehely 1967, *Theory of Orbits* (task #241).

INDEPENDENT, non-JPL cross-check anchors for our CR3BP core. Every EXPECTED
value below is transcribed from a printed table in:

    V. Szebehely, *Theory of Orbits: The Restricted Problem of Three Bodies*,
    Academic Press, 1967

cited per-test by book page / appendix / table. These are published numbers,
NOT values our own pipeline produced, so they satisfy the project's golden
discipline (the EXPECTED side traces to a source). Mining note:
docs/notes/2026-06-13-szebehely-1967-theory-of-orbits-mining.md (#185).

Two convention conversions are applied (see the mining note "Convention trap"):

1. FRAME MIRROR. Szebehely's standard reference system (Fig. 9.1a, p.449) puts
   the larger mass 1-mu at P1(+mu,0) and the smaller mass mu at P2(mu-1,0). Our
   ``core/cr3bp.py`` puts the larger mass 1-mu at (-mu,0) and the smaller at
   (1-mu,0) -- mirror-flipped in x. The Jacobi constant is mirror-invariant
   (depends on x^2, r1, r2), so C values transfer directly; collinear-point
   x-coordinates are negated (x -> -x) before feeding our routines.

2. TWO JACOBI CONSTANTS. Szebehely uses Omega_bar = r^2/2 + (1-mu)/r1 + mu/r2
   (-> C_bar, Wintner/Birkhoff constant) and the "standard" Omega = Omega_bar +
   mu(1-mu)/2 (-> C). Our ``jacobi_constant()`` evaluated at v=0 equals 2*Omega_bar
   = C_bar (no mu(1-mu) term). The Appendix I-IV tables print the STANDARD C, so
   we convert with C_bar = C - mu*(1-mu). Table III (p.457) prints BOTH C_bar and
   C and is used directly without conversion.
"""

from __future__ import annotations

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp


def _jc_at_rest(x: float, y: float, mu: float) -> float:
    """Our Jacobi constant at a point with zero velocity (= C_bar = 2*Omega_bar)."""
    return cr3bp.jacobi_constant(np.array([x, y, 0.0, 0.0, 0.0, 0.0]), mu)


# ---------------------------------------------------------------------------
# Triangular points L4/L5 -- exact algebraic anchors (location is EXACT).
# ---------------------------------------------------------------------------


def _l4_location(mu: float) -> tuple[float, float]:
    """L4 in OUR frame: equilateral apex of the primaries at (-mu,0),(1-mu,0).

    The two primaries are unit-distance apart, so the equilateral apex sits at
    x = (-mu + (1-mu))/2 = 1/2 - mu, y = +sqrt(3)/2 (r1 = r2 = 1). (For mu=1/2
    this gives x=0, matching Szebehely Table III after the frame mirror.)
    """
    return (0.5 - mu, float(np.sqrt(3.0) / 2.0))


def test_triangular_jacobi_copenhagen_mu_half() -> None:
    """L4/L5 Jacobi constant for the Copenhagen mass mu=1/2.

    Szebehely 1967, Table III (book p.457): L4 at x=0, y=+sqrt(3)/2 and L5 at
    x=0, y=-sqrt(3)/2 with C_bar = 2.7500 and C = 3.0000. For mu=1/2 the
    equilateral apex sits at x=0 exactly (location is exact, r1=r2=1), so this
    is a tight golden.

    Table III prints C_bar directly, so no convention conversion is needed; the
    frame mirror leaves x=0 unchanged and y is mirror-invariant.
    """
    mu = 0.5
    x4, y4 = _l4_location(mu)
    assert x4 == 0.0  # mu=1/2 apex on the y-axis in either frame
    for y in (y4, -y4):
        c_bar = _jc_at_rest(0.0, y, mu)
        assert c_bar == pytest.approx(2.7500, abs=5e-5)


def test_triangular_c_equals_3_all_mu() -> None:
    """Standard C(L4,5) = 3 exactly, independent of mu (Szebehely p.451).

    Szebehely states the algebraic identity C(L4,5) = 3 for ALL mu (the
    triangular points always lie at r1=r2=1, where the standard pseudo-potential
    takes the same value). Our jacobi_constant at v=0 gives C_bar; adding the
    convention term mu(1-mu) must recover the standard C = 3.0000 for every mu.
    Verified across the Earth-Moon mu and several sample mu in (0,1/2].
    """
    for mu in (0.01215068, 0.012, 0.1, 0.2, 0.3, 0.5):
        x4, y4 = _l4_location(mu)
        c_bar = _jc_at_rest(x4, y4, mu)
        c_standard = c_bar + mu * (1.0 - mu)
        assert c_standard == pytest.approx(3.0, abs=1e-12)


# ---------------------------------------------------------------------------
# Collinear points -- Copenhagen mu=1/2, Szebehely Table III (book p.457).
# x printed to 4 sig figs, so agreement is location-limited to ~4 figs.
# ---------------------------------------------------------------------------


def test_collinear_l2_copenhagen_mu_half_exact() -> None:
    """Inner collinear point at the origin for mu=1/2 (Szebehely Table III, p.457).

    For mu=1/2 the point between the equal masses sits at x=0 exactly (his frame;
    mirror leaves it at 0). Printed C_bar = 4.0000, C = 4.2500. r1=r2=1/2.
    Location is exact, so this is a tight golden.
    """
    mu = 0.5
    c_bar = _jc_at_rest(0.0, 0.0, mu)
    assert c_bar == pytest.approx(4.0000, abs=5e-5)
    assert c_bar + mu * (1.0 - mu) == pytest.approx(4.2500, abs=5e-5)


def test_collinear_l1_l3_copenhagen_mu_half() -> None:
    """Outer collinear points for mu=1/2 (Szebehely Table III, p.457).

    His frame: x=-1.1984 (1st) and x=+1.1984 (3rd), C_bar=3.4568, C=3.7068.
    Mirror x -> -x for our frame (the pair is symmetric for mu=1/2, so both map
    onto +-1.1984). x is printed to 4 sig figs -> the golden is location-limited;
    we assert agreement to the 4 printed figures.
    """
    mu = 0.5
    for x_szeb in (-1.1984, +1.1984):
        x_ours = -x_szeb  # frame mirror
        c_bar = _jc_at_rest(x_ours, 0.0, mu)
        # x is only good to 4 sig figs; |dC/dx| ~ O(1) near these points, so a
        # 1e-4 error in x propagates to ~1e-4 in C. Assert to the printed figures.
        assert c_bar == pytest.approx(3.4568, abs=2e-3)
        assert c_bar + mu * (1.0 - mu) == pytest.approx(3.7068, abs=2e-3)


# ---------------------------------------------------------------------------
# Collinear points -- high-precision Earth-Moon-range appendices.
# Columns are (mu, x_i [his frame], C_i [STANDARD C], Omega_xx, Omega_yy).
# x printed to ~13 sig figs; C is the standard constant -> convert C_bar = C - mu(1-mu).
# ---------------------------------------------------------------------------

# Szebehely 1967, Appendices I.D / II.D / III.D (book pp.216/220/224).
# (mu, x_i_his_frame, C_i_standard). The mu=0.0120 and mu=0.0123 rows are
# internally consistent to full double precision: recomputing C from the printed
# x in Szebehely's own frame reproduces the printed C to ~1e-11, and our code
# matches it to ~1e-11 -- these are tight goldens.
_APP_COLLINEAR_TIGHT = [
    # Appendix I.D, 1st collinear point (book p.216)
    (0.0120, -1.15510_01298, 3.18282_40063),
    (0.0123, -1.15625_37037, 3.18548_46255),
    # Appendix II.D, 2nd collinear point (book p.220)
    (0.0120, -0.83765_86648, 3.19880_45659),
    # Appendix III.D, 3rd collinear point (book p.224)
    (0.0120, +1.00499_99054, 3.02385_26541),
]

# The mu=0.0121 rows are INTERNALLY INCONSISTENT: recomputing the standard C from
# the printed x in Szebehely's own frame (C = x^2 + 2(1-mu)/r1 + 2mu/r2 + mu(1-mu))
# does NOT reproduce the printed C -- it disagrees by ~1e-7..3e-7, and our code
# agrees with that recompute (not with the printed C). Since the mu=0.0120 and
# mu=0.0123 rows are exact to ~1e-11, the printed (x, C) pair at mu=0.0121 cannot
# both be right; the most likely explanation is a last-digits transcription slip
# in the mining note (PDF re-check flagged in
# docs/notes/2026-06-13-szebehely-goldens.md). The printed C is still a published
# number, so we keep these rows as a LOOSE cross-check (agreement to ~3e-7) rather
# than fabricating tight agreement. The expected value is NEVER adjusted toward
# ours: it stays at Szebehely's printed C, and the tolerance documents the gap.
_APP_COLLINEAR_LOOSE = [
    (0.0121, -1.15548_72863, 3.18371_40528),  # I.D  (printed-vs-recompute ~3.0e-7)
    (0.0121, -0.83726_43231, 3.19982_77931),  # II.D (printed-vs-recompute ~1.1e-7)
    (0.0121, +1.00504_15697, 3.02405_03851),  # III.D (printed-vs-recompute ~2.0e-7)
]


@pytest.mark.parametrize(("mu", "x_szeb", "c_standard"), _APP_COLLINEAR_TIGHT)
def test_collinear_jacobi_earth_moon_range_tight(
    mu: float, x_szeb: float, c_standard: float
) -> None:
    """Jacobi constant at the printed collinear-point location matches Szebehely.

    Szebehely 1967, Appendices I.D / II.D / III.D (book pp.216/220/224): for each
    listed mu the appendix prints the collinear-point location x_i (~13 sig figs,
    his frame) and the STANDARD Jacobi constant C_i there. Feeding the
    mirror-mapped location (x -> -x) to our jacobi_constant at rest must reproduce
    C_bar = C_i - mu*(1-mu) to high precision.

    This is the tightest collinear cross-check: the location is published to full
    precision, so any disagreement is a real convention/formula defect, not a
    rounding artefact. Observed agreement on these rows is ~1e-11.
    """
    x_ours = -x_szeb  # frame mirror
    c_bar_expected = c_standard - mu * (1.0 - mu)
    c_bar_ours = _jc_at_rest(x_ours, 0.0, mu)
    assert c_bar_ours == pytest.approx(c_bar_expected, abs=1e-9)


@pytest.mark.parametrize(("mu", "x_szeb", "c_standard"), _APP_COLLINEAR_LOOSE)
def test_collinear_jacobi_earth_moon_range_loose(
    mu: float, x_szeb: float, c_standard: float
) -> None:
    """mu=0.0121 collinear rows -- loose cross-check (see _APP_COLLINEAR_LOOSE note).

    The printed (x, C) pair at mu=0.0121 is internally inconsistent at ~1e-7
    (suspected note transcription slip); the expected value is kept at the printed
    standard C and only the tolerance is relaxed. Documented in
    docs/notes/2026-06-13-szebehely-goldens.md.
    """
    x_ours = -x_szeb  # frame mirror
    c_bar_expected = c_standard - mu * (1.0 - mu)
    c_bar_ours = _jc_at_rest(x_ours, 0.0, mu)
    assert c_bar_ours == pytest.approx(c_bar_expected, abs=5e-7)


# ---------------------------------------------------------------------------
# Zero-velocity-curve smoke test -- Appendix IV C(x) table (book pp.226-229).
# ---------------------------------------------------------------------------
#
# DROPPED. The only Appendix IV sample captured in the mining note
# (mu=0.1, x=-1.0 -> C = 4.97029_60396) could not be reproduced in EITHER frame:
# computing the standard C = x^2 + 2(1-mu)/r1 + 2mu/r2 + mu(1-mu) directly in
# Szebehely's own frame at x=-1.0, mu=0.1 gives C = 4.7263636..., not 4.97029...,
# and our code matches 4.7263636... after the mirror. The 4.97029 sample is
# therefore a suspected note transcription error (the algebra forces 4.72636),
# and a golden cannot be anchored on an unverifiable number. Re-mine the actual
# Appendix IV print before wiring a C(x) smoke test. Documented in
# docs/notes/2026-06-13-szebehely-goldens.md.
