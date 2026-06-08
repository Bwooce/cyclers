"""Tier-1 Phase 3 (slow): ballistic_correct runs about Jupiter on an
Io-Europa-Ganymede chain.

NON-GOLDEN: this exercises the centre-agnostic corrector infrastructure (Lambert
solved with mu_central=mu_Jupiter; the centred moon ephemeris; the SATELLITES
bend lookup). The V_inf value is OUR computation, not a sourced anchor (the
Phase-7 gauntlet does the sourced comparison).

HONEST-RISK FINDING (moon-tour Tier-1, task #76; see plan honest-risk register):
the chain CLOSES about Jupiter (V_inf-magnitude continuity + periodicity residual
-> 0, converged), proving the corrector is centre-correct. But it is NOT
bend-feasible in this coplanar-circular no-V_inf-leveraging model: at the closed
geometry the per-encounter V_inf lands at ~9-12 km/s (the moons start collinear
on perfect circles with only t0/ToF free — no relative-phase or powered-leg
freedom), forcing 100-150 deg required turns while the small moons can bend only
2-5 deg at that V_inf. A real Io-Europa-Ganymede resonant cycler needs the
V_inf-leveraging (VILM) layer + Laplace-resonance phasing (Phase 5), which the
pure ballistic patched-conic corrector cannot represent. We do NOT loosen tol_kms
or fabricate a seed to force a green (golden / honest-risk rule); bend-feasibility
is recorded as a documented xfail until the VILM layer feeds seeds. A grid sweep
over period (12.5-42 d), t0, and ToFs in both magnitude and bend-aware vector
residual modes found no bend-feasible converged solution.
"""

from __future__ import annotations

import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.satellites import PRIMARIES
from cyclerfinder.search.correct import BallisticClosureResult, ballistic_correct


# A seed that converges (residual ~1e-13) about Jupiter — found by the seed sweep
# documented in the module docstring (t0=0.6 d, free ToFs (4.0, 3.4) d, slack
# leg = the Ganymede->Io leg, period = sum of the three moon periods).
def _solve_seed() -> BallisticClosureResult:
    """Run the converging Jovicentric I-E-G seed (shared by both slow tests)."""
    ephem = Ephemeris(model="circular", center="Jupiter")
    return ballistic_correct(
        sequence=("Io", "Europa", "Ganymede", "Io"),
        per_leg_revs=(0, 0, 0),
        per_leg_branch=("single", "single", "single"),
        t0_seed_sec=0.6 * 86400.0,
        tof_seed_days=(4.0, 3.4),
        period_sec=(1.769 + 3.551 + 7.155) * 86400.0,
        vinf_cap=20.0,
        slack_leg=2,
        ephem=ephem,
        mu_central=PRIMARIES["Jupiter"],
    )


@pytest.mark.slow
def test_jovian_ieg_chain_closes_about_jupiter() -> None:
    """The corrector CLOSES the chain about Jupiter — the centre-agnostic gate."""
    r = _solve_seed()
    # Closure (magnitude continuity + periodicity) is the centre-agnostic gate:
    # it proves mu_central + the centred ephemeris produce a self-consistent
    # Jovicentric chain.
    assert r.converged
    # The computed V_inf is Jovicentric (real conic excess about Jupiter), not a
    # heliocentric figure — sanity-bound, non-golden.
    assert all(0.0 < v < 20.0 for v in r.vinf_per_encounter_kms)


@pytest.mark.slow
@pytest.mark.xfail(
    strict=True,
    reason=(
        "HONEST-RISK (task #76): the coplanar-circular no-V_inf-leveraging "
        "patched-conic model closes the I-E-G chain about Jupiter but NOT "
        "bend-feasibly (required 100-150 deg turns vs 2-5 deg max-bend at the "
        "~10 km/s closed V_inf). Bend-feasible Jovian moon tours need the VILM "
        "layer (Phase 5) + Laplace-resonance phasing. We do not loosen tol or "
        "fabricate a seed. Flip to strict pass once the VILM layer feeds "
        "bend-feasible seeds. See module docstring + plan honest-risk register."
    ),
)
def test_jovian_ieg_chain_is_bend_feasible() -> None:
    r = _solve_seed()
    assert r.converged
    assert r.bend_feasible
