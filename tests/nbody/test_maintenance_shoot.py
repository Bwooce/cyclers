"""M7 leg-targeting solver (#423) — golden + Mars-perturbed smoke.

GOLDEN GATE: in Sun-only mode (``bodies=()``) the n-body position-targeting Newton
:func:`~cyclerfinder.nbody.maintenance_shoot.target_leg` must reproduce the INDEPENDENT
two-body :func:`cyclerfinder.core.lambert.lambert` departure velocity for the same
boundary-value leg. Lambert is the sourced/independent cross-check (a different solver,
analytic universal-variable), so agreement validates the targeting kernel without
circularity.
"""

from __future__ import annotations

import numpy as np
import pytest

rebound = pytest.importorskip("rebound")

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS  # noqa: E402
from cyclerfinder.core.ephemeris import Ephemeris  # noqa: E402
from cyclerfinder.core.flyby import flyby_dv  # noqa: E402
from cyclerfinder.core.lambert import lambert  # noqa: E402
from cyclerfinder.nbody.maintenance_shoot import (  # noqa: E402
    continuous_maintenance_chain,
    target_leg,
)
from cyclerfinder.nbody.propagator import RestrictedNBody  # noqa: E402

_DAY_S = 86400.0


def test_target_leg_sun_only_reproduces_lambert() -> None:
    """Sun-only targeting Newton recovers the two-body Lambert departure velocity.

    Real Earth->Mars leg over DE440: take r0 = Earth(t0), r_target = Mars(t1) for a
    ~210 d transfer, solve the two-body Lambert for the departure velocity, then run
    target_leg in Sun-only mode seeded with a deliberately PERTURBED guess. It must
    converge back to the Lambert velocity (km/s) and land sub-km on the target."""
    ephem = Ephemeris("astropy")
    prop = RestrictedNBody("rebound")

    t0 = 27.0 * 365.25 * _DAY_S  # ~2027, an arbitrary in-ephemeris epoch
    tof = 210.0 * _DAY_S
    t1 = t0 + tof
    r0, _ = ephem.state("E", t0)
    r_target, _ = ephem.state("M", t1)
    r0 = np.asarray(r0, dtype=np.float64)
    r_target = np.asarray(r_target, dtype=np.float64)

    sols = lambert(r0, r_target, tof, mu=MU_SUN_KM3_S2)
    v_lambert = np.asarray(sols[0].v1, dtype=np.float64)  # single-rev departure velocity

    # Seed with a perturbed guess (off by ~0.3 km/s in each component) so the Newton
    # iteration has real work to do — it must not just echo the seed.
    v_guess = v_lambert + np.array([0.3, -0.3, 0.2])
    res = target_leg(
        prop,
        r0,
        t0,
        t1,
        r_target,
        v_guess,
        bodies=(),  # Sun-only two-body — the Lambert regime
        ephem=ephem,
        tol_km=1.0,
    )

    assert res.converged, f"Sun-only targeting did not converge (miss {res.miss_km:.1f} km)"
    assert res.miss_km < 1.0, f"arrival miss {res.miss_km:.3f} km exceeds 1 km"
    # Recovers the independent Lambert departure velocity (mm/s-level agreement).
    assert res.v_dep_km_s == pytest.approx(v_lambert, abs=1e-5), (
        f"target_leg v_dep {res.v_dep_km_s} != lambert {v_lambert}"
    )


@pytest.mark.slow
def test_target_leg_mars_perturbed_converges_in_band() -> None:
    """Mars-perturbed targeting converges and lands far inside the Mars 3-SOI band.

    Same leg, now with Mars as a continuous perturber. The Mars-perturbed departure
    velocity differs from the two-body Lambert value (the perturbation is real work),
    but the solve must still drive the arrival miss sub-km — i.e. the STM Jacobian is
    correct under the perturber (the REBOUND-variation gravity-gradient gotcha is
    handled in the propagator)."""
    ephem = Ephemeris("astropy")
    prop = RestrictedNBody("rebound")

    t0 = 27.0 * 365.25 * _DAY_S
    tof = 210.0 * _DAY_S
    t1 = t0 + tof
    r0, _ = ephem.state("E", t0)
    r_target, _ = ephem.state("M", t1)
    r0 = np.asarray(r0, dtype=np.float64)
    r_target = np.asarray(r_target, dtype=np.float64)

    v_lambert = np.asarray(lambert(r0, r_target, tof, mu=MU_SUN_KM3_S2)[0].v1, dtype=np.float64)

    res = target_leg(prop, r0, t0, t1, r_target, v_lambert, bodies=("M",), ephem=ephem, tol_km=1.0)

    mars_soi_au = PLANETS["M"].sma_au * (PLANETS["M"].mu_km3_s2 / MU_SUN_KM3_S2) ** 0.4
    assert res.converged, f"Mars-perturbed targeting did not converge (miss {res.miss_km:.1f} km)"
    assert res.miss_km < 1.0, f"arrival miss {res.miss_km:.3f} km exceeds 1 km"
    assert res.miss_km / AU_KM < 3.0 * mars_soi_au  # trivially true; documents the band


def test_maintenance_chain_sun_only_matches_lambert_accounting() -> None:
    """Sun-only 3-node E-M-E chain: the flyby ΔV equals the independent Lambert accounting.

    Build a real E->M->E node sequence from DE440. In Sun-only mode each targeted leg
    is the two-body Lambert (Task-1 golden), so the interior Mars-flyby maintenance ΔV
    must equal ``flyby_dv_for("M", vinf_in, vinf_out)`` computed directly from the two
    independent Lambert legs' velocities. This validates the chain's per-node ΔV
    accounting end-to-end against an independent solver (no circularity)."""
    ephem = Ephemeris("astropy")
    prop = RestrictedNBody("rebound")

    t0 = 27.0 * 365.25 * _DAY_S
    t1 = t0 + 210.0 * _DAY_S
    t2 = t1 + 360.0 * _DAY_S
    epochs = [t0, t1, t2]
    chain_bodies = ["E", "M", "E"]

    r_e0, _ = (np.asarray(x, dtype=np.float64) for x in ephem.state("E", t0))
    r_m1, v_m1 = (np.asarray(x, dtype=np.float64) for x in ephem.state("M", t1))
    r_e2, _ = (np.asarray(x, dtype=np.float64) for x in ephem.state("E", t2))

    # Independent two-body Lambert accounting for the interior Mars flyby, using the
    # SAME min-flyby-altitude floor the chain defaults to (Russell 200 km).
    leg0 = lambert(r_e0, r_m1, t1 - t0, mu=MU_SUN_KM3_S2)[0]  # E->M
    leg1 = lambert(r_m1, r_e2, t2 - t1, mu=MU_SUN_KM3_S2)[0]  # M->E
    vinf_in_mars = np.asarray(leg0.v2, dtype=np.float64) - v_m1  # arrival at Mars
    vinf_out_mars = np.asarray(leg1.v1, dtype=np.float64) - v_m1  # departure from Mars
    rp_min_mars = PLANETS["M"].radius_eq_km + PLANETS["M"].safe_alt_km  # per-body sourced floor
    expected_dv_kms = flyby_dv(vinf_in_mars, vinf_out_mars, PLANETS["M"].mu_km3_s2, rp_min_mars)

    chain = continuous_maintenance_chain(
        epochs, chain_bodies, ephem, prop, cruise_perturbers=(), n_cycles=1
    )

    assert not chain.diverged, "Sun-only E-M-E chain should fully converge"
    assert chain.n_legs_converged == 2
    assert len(chain.nodes) == 3
    assert chain.nodes[0].is_flyby is False  # departure node, no flyby ΔV
    assert chain.nodes[2].is_flyby is False  # final arrival node, no flyby ΔV
    mars_node = chain.nodes[1]
    assert mars_node.is_flyby and mars_node.body == "M"
    # The chain's flyby ΔV reproduces the independent Lambert-based accounting.
    assert mars_node.dv_kms == pytest.approx(expected_dv_kms, rel=1e-3, abs=1e-4)
    assert chain.horizon_tcm_mps == pytest.approx(expected_dv_kms * 1000.0, rel=1e-3, abs=0.1)


@pytest.mark.slow
def test_maintenance_chain_naive_endpoint_perturber_diverges_but_exclusion_fixes_it() -> None:
    """Endpoint-exclusion is what makes the perturbed chain work (the artifact + its fix).

    PINNED FINDING (the patched-conic handoff the M7 plan flagged as hard-part #1). With
    ``exclude_endpoint_bodies=False`` the M->E return leg is propagated with Mars as a
    continuous perturber while STARTING at the Mars node — the spacecraft begins at the
    Mars centre, inside the softened core, and the targeting cannot converge. That is a
    modelling ARTIFACT, not a fuel cost. With the default ``exclude_endpoint_bodies=True``
    Mars is dropped from each leg whose endpoint it is (the flyby is patched-conic at the
    node), so the SAME ``cruise_perturbers=("M",)`` request converges — proving the
    general perturber rule removes the artifact while keeping every genuine (non-endpoint)
    third-body perturbation."""
    ephem = Ephemeris("astropy")
    prop = RestrictedNBody("rebound")

    t0 = 27.0 * 365.25 * _DAY_S
    epochs = [t0, t0 + 210.0 * _DAY_S, t0 + 570.0 * _DAY_S]

    # Naive: Mars integrated continuously from the Mars node -> divergence artifact.
    naive = continuous_maintenance_chain(
        epochs,
        ["E", "M", "E"],
        ephem,
        prop,
        cruise_perturbers=("M",),
        exclude_endpoint_bodies=False,
    )
    assert naive.diverged, "expected the naive Mars-from-centre departure leg to diverge"
    assert naive.horizon_tcm_mps == float("inf")  # unmeasurable -> honest, not forced

    # The fix: endpoint exclusion drops Mars on its own legs -> converges (== Sun-cruise).
    fixed = continuous_maintenance_chain(
        epochs, ["E", "M", "E"], ephem, prop, cruise_perturbers=("M",), exclude_endpoint_bodies=True
    )
    assert not fixed.diverged, "endpoint exclusion should remove the handoff artifact"
    assert fixed.horizon_tcm_mps < float("inf")


@pytest.mark.slow
def test_maintenance_chain_nonendpoint_perturber_is_captured() -> None:
    """A genuine third-body (non-endpoint) perturber is captured, not excluded.

    The general-system case (moon tours / Venus-flyby chains have real cross-cruise
    perturbers). On the E->M->E sequence, Jupiter is NEVER a leg endpoint, so with
    ``cruise_perturbers=("J",)`` it perturbs every cruise leg. The chain must still
    converge (Jupiter is far, no centre-of-body artifact) and its horizon TCM must
    differ from the Sun-only value — i.e. real third-body perturbation moves the
    measurement, confirming the model is not silently Sun-only for non-Earth-Mars
    systems."""
    ephem = Ephemeris("astropy")
    prop = RestrictedNBody("rebound")

    t0 = 27.0 * 365.25 * _DAY_S
    epochs = [t0, t0 + 210.0 * _DAY_S, t0 + 570.0 * _DAY_S]

    sun_only = continuous_maintenance_chain(
        epochs, ["E", "M", "E"], ephem, prop, cruise_perturbers=()
    )
    jovian = continuous_maintenance_chain(
        epochs, ["E", "M", "E"], ephem, prop, cruise_perturbers=("J",)
    )

    assert not sun_only.diverged and not jovian.diverged
    assert jovian.horizon_tcm_mps < float("inf")
    # Jupiter genuinely perturbs the cruise -> a different (here, tiny but non-zero) shift.
    assert jovian.horizon_tcm_mps != sun_only.horizon_tcm_mps


@pytest.mark.slow
def test_maintenance_chain_s1l1_reproduces_published_strictly_ballistic() -> None:
    """M7 on the real multi-rev S1L1 cycler reproduces the PUBLISHED strictly-ballistic tier.

    The decisive whole-stack reproduction: walk S1L1's full 22-node App-C sequence with the
    position-targeted chain, SEEDED PER LEG by the sourced departure V_inf (so the multi-rev
    resonant E-E legs stay in the right basin — a single-rev re-guess otherwise inflates the
    TCM ~4 orders of magnitude), using Russell 2004's SOURCED 200 km flyby floor (p.165:
    r_p,min 6578.0 km Earth / 3598.5 km Mars). Under that floor the chain measures ~0 m/s =
    strictly ballistic, REPRODUCING the published claim (Russell-Ocampo: S1L1 strictly
    ballistic). At the old unsourced 300 km floor a single marginal flyby (node 19, required
    95.0 deg bend vs 94.5 deg max) spuriously charged ~40 m/s — the floor, not the cycler,
    was the cost. See docs/notes/2026-06-23-m7-phase1-results.md. First clean multi-leg
    cycler reproduction to its published maintenance tier."""
    from cyclerfinder.search import s1l1_corrected as s1l1

    ephem = Ephemeris("astropy")
    prop = RestrictedNBody("rebound")

    node_epochs = [(s1l1.APPC_EPOCH_DAYS + ts) * s1l1.DAY_S for (_n, _b, ts, _v) in s1l1.APPC_LEGS]
    node_bodies = [b for (_n, b, _ts, _v) in s1l1.APPC_LEGS]
    # Per-leg seed = sourced departure V_inf at the departing node (rev-correct basin).
    leg_seeds = []
    for i in range(len(s1l1.APPC_LEGS) - 1):
        _n, body_i, _ts, vinf_i = s1l1.APPC_LEGS[i]
        _r, v_pl = (np.asarray(x, dtype=np.float64) for x in ephem.state(body_i, node_epochs[i]))
        leg_seeds.append(v_pl + s1l1._vinf_vec(vinf_i))

    chain = continuous_maintenance_chain(
        node_epochs,
        node_bodies,
        ephem,
        prop,
        cruise_perturbers=(),  # Earth-Mars: Sun-cruise (endpoint exclusion would do the same)
        leg_v_guess=leg_seeds,
        n_cycles=7,
    )

    assert not chain.diverged, "S1L1 sourced-seeded chain should fully converge"
    assert chain.n_legs_converged == chain.n_legs_total
    # Reproduces the published STRICTLY BALLISTIC tier under Russell's 200 km floor (~0 m/s).
    assert chain.horizon_tcm_mps < 1.0, (
        f"S1L1 M7 horizon {chain.horizon_tcm_mps:.2f} m/s should be ~0 (strictly ballistic) "
        f"at the sourced 200 km flyby floor"
    )
    # Per-node minimal-thrust flyby heights (#427): the one binding flyby (node 19) sits
    # AT the 200 km floor; every other flyby flies far higher (needs less bend). This is
    # the structured per-row tuple the catalogue/website consume.
    flybys = [nd for nd in chain.nodes if nd.is_flyby]
    assert all(nd.flyby_alt_km >= 200.0 - 1e-6 for nd in flybys)  # never below the floor
    assert min(nd.flyby_alt_km for nd in flybys) == pytest.approx(200.0, abs=2.0)  # node 19
    assert max(nd.flyby_alt_km for nd in flybys) > 1000.0  # gentle flybys fly high
