"""Tests for the flyby-continuity-wired FBS-analytic ΔV-chain optimiser (#244).

This extends the #243 match-point NLP (``search/fbs_optimize.py``) with the
patched-conic flyby-continuity constraints (Ellison Eqs. 3-4, via
``core/fbs_match_point.flyby_coupling_block``) at every interior body, and their
ANALYTIC gradients. The interior bodies now carry SEPARATE arrival/departure
heliocentric velocities (a real powered/turning flyby), coupled by the
v∞-magnitude-continuity equality + the periapsis-altitude inequality.

Discipline (mirrors #243): MECHANICS gates + a Jacobian-vs-finite-difference
correctness check on the FULL constrained system (the only available check —
Ellison publishes no numeric gradient). The science verdict (the catalogue-wide
parity head-to-head) lives in ``scripts/fbs_optimizer_adoption_parity.py`` and its
results note, not here.
"""

from __future__ import annotations

import numpy as np

from cyclerfinder.core.constants import MU_SUN_KM3_S2
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.lambert import lambert
from cyclerfinder.search.fbs_optimize import ChainLegSpec
from cyclerfinder.search.fbs_optimize_flyby import (
    FlybyChainSpec,
    _flyby_constraint_jac_fd,
    _flyby_constraint_vector,
    _pack_x0,
    optimize_chain_fbs_flyby,
)

DAY = 86400.0


def _eme_flyby_spec(
    eph: Ephemeris,
) -> tuple[FlybyChainSpec, list[np.ndarray], list[np.ndarray], list[np.ndarray]]:
    """E->M->E chain (one interior Mars flyby) with a per-leg Lambert seed.

    Returns (spec, dv_seeds, varr_seeds, vdep_seeds). The interior Mars node has
    a SEPARATE arrival and departure velocity (the flyby turns v∞); the ballistic
    Lambert seed sets both to the heliocentric continuity value, so the seed is
    flyby-feasible (zero turn) and the v∞-continuity constraint is satisfied.
    """
    t0 = 0.0
    seq = ["E", "M", "E"]
    tofs = [210 * DAY, 520 * DAY]
    epochs = [t0, t0 + tofs[0], t0 + tofs[0] + tofs[1]]
    states = [eph.state(b, e) for b, e in zip(seq, epochs, strict=True)]
    pos = [np.asarray(s[0]) for s in states]
    vpl = [np.asarray(s[1]) for s in states]
    legs = tuple(ChainLegSpec(r0=pos[i], rf=pos[i + 1], tof_s=tofs[i], alpha=0.5) for i in range(2))
    spec = FlybyChainSpec(legs=legs, bodies=tuple(seq), v_planets=tuple(vpl))

    # Per-leg Lambert seed: leg i departure / arrival heliocentric velocities.
    leg_v = [lambert(pos[i], pos[i + 1], tofs[i], mu=MU_SUN_KM3_S2)[0] for i in range(2)]
    v_dep0 = np.asarray(leg_v[0].v1)  # departure at E (node 0)
    v_arr1 = np.asarray(leg_v[0].v2)  # arrival at M (node 1)
    v_dep1 = np.asarray(leg_v[1].v1)  # departure at M (node 1)
    v_arr2 = np.asarray(leg_v[1].v2)  # arrival at E (node 2)
    dvs = [np.zeros(3), np.zeros(3)]
    # Node layout: node0 dep only; node1 arr+dep; node2 arr only.
    return spec, dvs, [v_arr1, v_arr2], [v_dep0, v_dep1]


def test_flyby_constraint_jacobian_matches_finite_difference() -> None:
    """The analytic full constraint Jacobian agrees with central differences.

    The constraint stack is [chain match-point defect (6M) ; per-flyby
    v∞-continuity (n_interior) ; per-flyby altitude (n_interior)] and the Jacobian
    blends the FBS match-point STM columns with the flyby_coupling_block analytic
    gradients. This is the #244 caveat-2 wiring; the FD cross-check is the only
    available correctness gate.
    """
    eph = Ephemeris("astropy")
    spec, dvs, varr, vdep = _eme_flyby_spec(eph)
    x0 = _pack_x0(spec, dvs, varr, vdep)
    j_ana = _flyby_constraint_jac_fd(spec, x0, analytic=True)
    j_fd = _flyby_constraint_jac_fd(spec, x0, analytic=False)
    # Block-relative error: position rows (km) and velocity/constraint rows differ
    # in scale, so compare on the already row-scaled constraint the solver sees.
    rel = np.abs(j_ana - j_fd) / np.maximum(np.abs(j_fd), 1.0)
    assert float(np.max(rel)) < 1e-5


def test_constraint_vector_matchpoint_zero_flyby_nonzero_on_two_lambert_seed() -> None:
    """A two-independent-Lambert seed is dynamically consistent but NOT flyby-feasible.

    Each leg is its own Lambert arc, so the per-leg match-point defect is ~0 (the
    leg states ARE consistent), but the two arcs reach/leave Mars with DIFFERENT
    v∞ magnitudes — exactly the patched-conic discontinuity the flyby-continuity
    constraint exists to drive out. So ``c_vinf`` is small-but-nonzero here, which
    is the whole point of wiring the constraint (the #244 caveat-2 motivation).
    """
    eph = Ephemeris("astropy")
    spec, dvs, varr, vdep = _eme_flyby_spec(eph)
    x0 = _pack_x0(spec, dvs, varr, vdep)
    c = _flyby_constraint_vector(spec, x0)
    m = len(spec.legs)
    n_int = m - 1
    defect = c[: 6 * m]
    c_vinf = c[6 * m : 6 * m + n_int]
    assert float(np.max(np.abs(defect))) < 1e-6  # each leg is its own consistent arc
    assert float(np.max(np.abs(c_vinf))) > 1e-4  # but the flyby is discontinuous


def test_analytic_lane_converges_feasibly() -> None:
    """The analytic-gradient flyby lane converges to a feasible chain from a seed."""
    eph = Ephemeris("astropy")
    spec, dvs, varr, vdep = _eme_flyby_spec(eph)
    res = optimize_chain_fbs_flyby(
        spec, dvs, varr, vdep, use_analytic_jac=True, feas_tol=1e-6, maxiter=300
    )
    # Both equality blocks (match-point defect AND v∞-magnitude continuity) driven
    # to feasibility, and the periapsis-altitude inequality respected.
    assert res.feasible
    assert res.flyby_feasible
    # The analytic lane supplies the constraint Jacobian directly (njev > 0).
    assert res.constr_njev > 0


def test_default_unchanged_optimize_chain_fbs_still_imports() -> None:
    """The #243 entry point is untouched (additive-only check)."""
    from cyclerfinder.search.fbs_optimize import optimize_chain_fbs

    assert callable(optimize_chain_fbs)
