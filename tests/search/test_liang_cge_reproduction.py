"""Golden reproduction tests: Liang et al. 2024 idealized CGE triple cyclers.

EXPECTED side: the printed Tables 2-7 of Liang, Yang, Li, Bai & Qin,
"Callisto-Ganymede-Europa Triple Cyclers", *JGCD* (Engineering Note), 2024,
DOI 10.2514/1.G008387 — initial moon phases (Tables 2/4/6) and per-flyby
transit time / V-infinity (Tables 3/5/7) for Members A, B, C, transcribed in
:mod:`cyclerfinder.search.cge_scaffold` (re-verified character-by-character
in the 2026-06-12 transcription rescan, mining note).

ACTUAL side: our same-model reproduction (circular-coplanar two-body about
Jupiter at the Table 1 mean motions, fixed-time Lambert legs at the printed
ToFs) — :func:`cyclerfinder.search.cge_scaffold.reproduce_member`.

TOLERANCE RATIONALE: see :func:`cge_scaffold.vinf_print_tolerance_kms`.
Summary — the binding precision limit is the 4-decimal print quantization of
the Table 1 mean motions (half-ULP 5e-5 rad/day), which over flyby epochs up
to ~119 days bounds any reproduction's V-infinity agreement at
2 * v_moon * 5e-5 * t (0.017..0.14 km/s per flyby), plus a 1e-3 km/s floor
for the unprinted mu_Jupiter and the phase/ToF/V-infinity print quantization.
The wrong-anchor / wrong-Lambert-branch alternatives miss by >= 0.18 km/s
(asserted via selection margins), so the tolerance separates "reproduces
within input precision" from "does not reproduce" with >= 4x margin at every
flyby. Observed worst residuals (2026-06-13): A 1.52e-2, B 1.38e-2,
C 4.82e-2 km/s — each 2.8x..10x INSIDE its per-flyby tolerance. Do not widen
these tolerances; a member that fails them is a documented negative.
"""

from __future__ import annotations

from functools import cache

import pytest

from cyclerfinder.search.cge_scaffold import (
    LIANG_MEMBERS,
    MemberReproduction,
    derived_radius_km,
    reproduce_member,
    vinf_print_tolerance_kms,
)

MEMBERS = ("A", "B", "C")


@cache
def _repro(member: str) -> MemberReproduction:
    return reproduce_member(member)


def test_perijove_matches_printed_anchor() -> None:
    """r_Eu - 10000 km must hit the paper's "about 660988 km" (p. 14).

    Our derivation (registry Jupiter GM + Table 1 Europa mean motion) gives
    660993.5 km — 5.5 km from the printed rounding. 50 km tolerance covers
    the full plausible mu_Jupiter spread (r scales as mu^(1/3): the
    planet-vs-system GM ambiguity moves r_Europa by ~45 km).
    """
    r_p = derived_radius_km("Europa") - 10000.0
    assert abs(r_p - 660988.0) < 50.0


@pytest.mark.parametrize("member", MEMBERS)
def test_phase_convention_self_consistency(member: str) -> None:
    """Callisto's printed phase puts it AT the conic crossing at t_c0.

    Validates the geometry convention (phases at perijove departure t = 0;
    structure flag 1 => crossing at E = +pi/2). Tolerance: Table 2/4/6
    phase print quantization (<= 5e-5 rad) + Table 1 mean-motion half-ULP
    over t_c0 ~ 2.4 d (1.2e-4 rad) + eccentricity sensitivity to the
    unprinted mu_Jupiter => 2e-4 rad. Observed: A/B 6.5e-6, C 5.3e-5 rad.
    """
    rep = _repro(member)
    assert rep.callisto_angle_residual_rad < 2.0e-4


@pytest.mark.parametrize("member", MEMBERS)
def test_legs_are_multirev_as_published(member: str) -> None:
    """Every transfer leg resolves to a one-revolution Lambert solution.

    The paper forces "all transfer arcs ... more than one revolution"
    (p. 9): n_revs = 1 (between 1 and 2 revolutions) for all four legs.
    The high/low/high/low branch pattern is the identification regression
    (selection is unambiguous; see the margin test).
    """
    rep = _repro(member)
    assert [leg.n_revs for leg in rep.legs] == [1, 1, 1, 1]
    assert [leg.branch for leg in rep.legs] == ["high", "low", "high", "low"]


@pytest.mark.parametrize("member", MEMBERS)
def test_leg_selection_unambiguous(member: str) -> None:
    """The chosen Lambert branch beats the runner-up by >= 0.18 km/s.

    Guards the identification step: per-leg branch selection minimises the
    distance to the printed V-infinity pair, which is only meaningful if
    the alternatives are far outside the print-precision tolerance
    (<= 0.14 km/s). 0.18 km/s ~ 1.3x the largest per-flyby tolerance.
    """
    rep = _repro(member)
    for leg in rep.legs:
        assert leg.selection_margin_kms > 0.18, f"leg {leg.index} ambiguous: {leg}"


@pytest.mark.parametrize("member", MEMBERS)
def test_vinf_matches_printed_tables(member: str) -> None:
    """GOLDEN: per-flyby |V_inf| (in and out) matches Tables 3/5/7.

    Per-flyby tolerance from vinf_print_tolerance_kms (module docstring
    rationale). This asserts EVERY V-infinity the tables print, on both
    the inbound and outbound side of every flyby — not a subset.
    """
    rep = _repro(member)
    for fb in rep.flybys:
        tol = vinf_print_tolerance_kms(fb.epoch_days, fb.moon, rep.radii_km)
        assert fb.residual_in_kms < tol, f"flyby {fb.index} ({fb.moon}) in: {fb}"
        if fb.residual_out_kms is not None:
            assert fb.residual_out_kms < tol, f"flyby {fb.index} ({fb.moon}) out: {fb}"


@pytest.mark.parametrize("member", MEMBERS)
def test_vinf_continuity_at_flybys(member: str) -> None:
    """GOLDEN: ballistic flyby => inbound and outbound |V_inf| agree.

    The members are published as ballistic (residual defect Delta-v below
    1e-8 m/s, p. 13); in an exact reconstruction the magnitudes would be
    identical. Same print-precision tolerance as the table match (the
    continuity gap is a difference of two quantities with that error).
    """
    rep = _repro(member)
    for fb in rep.flybys:
        if fb.continuity_kms is not None:
            tol = vinf_print_tolerance_kms(fb.epoch_days, fb.moon, rep.radii_km)
            assert fb.continuity_kms < tol, f"flyby {fb.index} ({fb.moon}): {fb}"


@pytest.mark.parametrize("member", MEMBERS)
def test_initial_conic_vinf(member: str) -> None:
    """The initial-conic Callisto V_inf reproduces the tables' flyby-0 value.

    This side never touches a Lambert solve: it is pure (a, e, r_p) conic
    geometry against Callisto's circular velocity at the crossing — the
    cleanest single anchor. Observed: A/B 1.6e-5, C 4.1e-4 km/s. Tolerance
    1e-3 km/s (the mining-note section 5 guidance level: mu_Jupiter +
    print quantization, no mean-motion time accumulation on this term).
    """
    rep = _repro(member)
    printed = LIANG_MEMBERS[member].vinf_printed_kms[0]
    assert abs(rep.conic_vinf_kms - printed) < 1.0e-3


@pytest.mark.parametrize("member", MEMBERS)
def test_cycle_tof_in_published_range(member: str) -> None:
    """Sum of printed leg ToFs lands in the paper's per-cycle envelope.

    "All cycle ToFs ~ 100 d", spread <= 0.35 d (p. 18); Member A's Fig. 5d
    range is 99.86-100.14 d. Pure consistency check on the printed inputs.
    """
    rep = _repro(member)
    assert abs(rep.cycle_tof_days - 100.0) <= 0.35
    if member == "A":
        assert 99.86 <= rep.cycle_tof_days <= 100.14


@pytest.mark.parametrize("member", MEMBERS)
def test_defect_altitude_convention_ballpark(member: str) -> None:
    """Coarse structural check of the paper's no-Jupiter altitude fiction.

    For the two tight flybys of each cycle (#2 Callisto, #3 Europa: turn
    angles O(0.5-1 rad)), our defect_altitude_km reproduces the printed
    "Flyby Altitude" within ~10% (observed ratios 0.93..1.05). Asserted at
    [0.75, 1.33] — a coarse turn-angle-structure check ONLY: the printed
    altitudes are declared Delta-v-equivalent fictions (paper p. 16), are
    NOT catalogue anchors, and the near-zero-turn flybys (#0, #1), whose
    fictitious altitudes diverge, are deliberately not compared.
    """
    rep = _repro(member)
    spec = LIANG_MEMBERS[member]
    for idx in (2, 3):
        fb = rep.flybys[idx]
        assert fb.defect_altitude_km is not None
        ratio = fb.defect_altitude_km / spec.altitudes_printed_km[idx]
        assert 0.75 < ratio < 1.33, f"flyby {idx} ({fb.moon}) altitude ratio {ratio:.3f}"


@pytest.mark.parametrize("member", MEMBERS)
def test_worst_residual_recorded_level(member: str) -> None:
    """Regression pin: the worst V-infinity residual stays at its 2026-06-13
    level (A 1.52e-2, B 1.38e-2, C 4.82e-2 km/s) and never silently grows.

    This is OUR OWN computed value (not a published anchor): it exists so a
    future change to the Lambert solver / registry constants that degrades
    the reproduction inside the print-precision tolerance is still noticed.
    Bound = observed * 1.05 (constants-stable headroom only).
    """
    observed = {"A": 1.52e-2, "B": 1.38e-2, "C": 4.82e-2}
    rep = _repro(member)
    assert rep.max_vinf_residual_kms < observed[member] * 1.05


def test_member_b_shares_first_two_flybys_with_a() -> None:
    """Members A and B print identical first two flybys (Table 5 caveat:
    "First two rows identical to Member A — the trajectories diverge at
    the Europa branch"). Our reconstructions must agree there too."""
    rep_a = _repro("A")
    rep_b = _repro("B")
    for idx in (0, 1):
        assert rep_a.flybys[idx].vinf_in_kms == pytest.approx(
            rep_b.flybys[idx].vinf_in_kms, abs=1e-12
        )
        assert rep_a.flybys[idx].vinf_out_kms == pytest.approx(
            rep_b.flybys[idx].vinf_out_kms, abs=1e-12
        )
