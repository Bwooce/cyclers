"""Perturbed-parity gate for the analytic Jovian state+STM co-integrator (#480).

The decisive test is at a NEAR-FLYBY leg (the EGGIE leg starting at the interior
Ganymede periapsis, ~6e3 km from the moon, deep inside its SOI) — NOT a deep-space
leg. A Jupiter-only parity test passes even with the moon gravity-gradient bug
(``memory/reference_rebound_variation_custom_force_gotcha``), so the gate MUST be
where the moon term dominates the gradient.

Three checks at the near-flyby leg:
1. First-order match: |Phi (eps d) - (f(x0+eps d) - f(x0))| / |f(x0+eps d) - f(x0)|
   <= ~1e-5 at small eps, over several random unit perturbations d.
2. Second-order consistency: the un-modelled residual |f(x0+eps d) - f(x0)
   - Phi (eps d)| / |eps d| shrinks ~linearly as eps -> 0 (it is O(eps), i.e. the
   raw residual is O(eps^2)) — confirms Phi is the true Jacobian, not a fit.
3. State parity: rf, vf from the analytic co-integrator match REBOUND's
   :meth:`JovianRestrictedNBody.propagate` endpoint to ~1e-4 rel over the same leg,
   so the STM is for the SAME dynamics the corrector residual uses.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from cyclerfinder.core.satellites import SATELLITES
from cyclerfinder.nbody.jovian import (
    MU_JUPITER_KM3_S2,
    JovianRailsCache,
    JovianRestrictedNBody,
)
from cyclerfinder.nbody.jovian_ideal import (
    EGGIE_MOONS,
    build_eggie_periapsis_seed,
    ideal_ephemeris_from_guess,
)
from cyclerfinder.nbody.jovian_stm import propagate_with_stm
from cyclerfinder.search.resonant_conic import eggie_refined_guess

# Per-component FD probe scale (km for position, km/s for velocity) — the v->r STM
# block scales like the time-of-flight (~1e6 s), so velocity columns need a far
# smaller probe than position columns (mirrors tests/nbody/test_propagator_stm.py).
_PROBE_SCALE = np.array([1.0, 1.0, 1.0, 1e-3, 1e-3, 1e-3])


@pytest.fixture(scope="module")
def flyby_leg() -> dict[str, Any]:
    """The EGGIE leg that starts closest to a Galilean moon (interior flyby leg).

    Returns the chosen leg's endpoints/epochs, the ideal ephemeris, and the
    closest-approach distance so the test can assert the gate is genuinely a
    near-flyby (within the moon SOI), per the memory rule.
    """
    guess = eggie_refined_guess()
    ephem = ideal_ephemeris_from_guess(guess)
    seed = build_eggie_periapsis_seed(guess)
    n = len(seed.sequence)

    # Pick the interior leg whose START node is closest to its moon (the strongest
    # moon gravity-gradient at leg start). Interior nodes are 1..n-2.
    best_i = None
    best_d = np.inf
    for i in range(1, n - 1):
        moon = seed.sequence[i]
        r_m, _ = ephem.state(moon, seed.epochs[i])
        d = float(np.linalg.norm(seed.node_states[i][:3] - r_m))
        if d < best_d:
            best_d, best_i = d, i
    assert best_i is not None

    moon = seed.sequence[best_i]
    sat = SATELLITES[moon]
    soi = sat.sma_km * (sat.mu_km3_s2 / (3.0 * MU_JUPITER_KM3_S2)) ** (1.0 / 3.0)
    return {
        "ephem": ephem,
        "r0": np.asarray(seed.node_states[best_i][:3], dtype=np.float64),
        "v0": np.asarray(seed.node_states[best_i][3:], dtype=np.float64),
        "t0": float(seed.epochs[best_i]),
        "t1": float(seed.epochs[best_i + 1]),
        "moon": moon,
        "approach_km": best_d,
        "soi_km": soi,
    }


def _f(r0: np.ndarray, v0: np.ndarray, leg: dict[str, Any]) -> np.ndarray:
    rf, vf, _ = propagate_with_stm(
        r0, v0, leg["t0"], leg["t1"], ephem=leg["ephem"], moons=EGGIE_MOONS
    )
    return np.concatenate([rf, vf])


def test_gate_leg_is_a_real_flyby(flyby_leg: dict[str, Any]) -> None:
    """Sanity: the gate leg starts well inside a moon's SOI (moon gradient dominates).

    At the chosen leg start the moon term must materially shape the gradient;
    otherwise the gate cannot detect the gravity-gradient bug.
    """
    assert flyby_leg["approach_km"] < flyby_leg["soi_km"], (
        f"gate leg start {flyby_leg['approach_km']:.0f} km is outside "
        f"{flyby_leg['moon']} SOI {flyby_leg['soi_km']:.0f} km — not a flyby gate"
    )
    # The moon gravity gradient (mu_m/d^3) must beat Jupiter's (mu_J/r^3) at start.
    sat = SATELLITES[flyby_leg["moon"]]
    d = flyby_leg["approach_km"]
    r = float(np.linalg.norm(flyby_leg["r0"]))
    moon_grad = sat.mu_km3_s2 / d**3
    jup_grad = MU_JUPITER_KM3_S2 / r**3
    assert moon_grad > jup_grad, (
        f"moon gradient {moon_grad:.2e} <= Jupiter {jup_grad:.2e}; gate too weak"
    )


def test_stm_perturbed_parity_at_flyby(flyby_leg: dict[str, Any]) -> None:
    """First-order match + second-order shrink at the near-flyby leg."""
    leg = flyby_leg
    rf0, vf0, phi = propagate_with_stm(
        leg["r0"], leg["v0"], leg["t0"], leg["t1"], ephem=leg["ephem"], moons=EGGIE_MOONS
    )
    f0 = np.concatenate([rf0, vf0])
    x0 = np.concatenate([leg["r0"], leg["v0"]])

    rng = np.random.default_rng(480)
    epsilons = [1e-1, 1e-2, 1e-3]
    first_order_rels: list[float] = []
    for trial in range(3):
        d = rng.standard_normal(6)
        d /= np.linalg.norm(d)
        dx_unit = _PROBE_SCALE * d  # column-scaled perturbation direction

        norm_resids: list[float] = []  # |raw residual| / |eps*dx_unit|
        for eps in epsilons:
            dx = eps * dx_unit
            f_pert = _f(x0[:3] + dx[:3], x0[3:] + dx[3:], leg)
            delta = f_pert - f0
            predicted = phi @ dx
            raw_resid = float(np.linalg.norm(delta - predicted))
            norm_resids.append(raw_resid / float(np.linalg.norm(dx)))
            if eps == epsilons[-1]:  # first-order rel at the smallest eps
                first_order_rels.append(raw_resid / float(np.linalg.norm(delta)))

        # Second-order: residual/|dx| ~ O(eps); a 10x eps drop should shrink it by
        # ~10x (allow 0.2 slack for integration-tolerance floor at the smallest eps).
        shrink_10 = norm_resids[1] / norm_resids[0]  # eps 1e-2 vs 1e-1
        shrink_100 = norm_resids[2] / norm_resids[1]  # eps 1e-3 vs 1e-2
        print(
            f"trial {trial}: norm_resids(eps=1e-1,1e-2,1e-3)="
            f"{norm_resids[0]:.3e},{norm_resids[1]:.3e},{norm_resids[2]:.3e} "
            f"shrink_10={shrink_10:.3f} shrink_100={shrink_100:.3f}"
        )
        assert shrink_10 < 0.2, f"trial {trial}: residual/|dx| not shrinking (Phi wrong)"

    fo = float(np.max(first_order_rels))
    print(f"max first-order rel error at eps=1e-3: {fo:.3e}")
    assert fo < 1e-5, f"first-order STM match {fo:.3e} > 1e-5 at the flyby leg"


def test_stm_state_matches_rebound_endpoint(flyby_leg: dict[str, Any]) -> None:
    """The analytic co-integrator's endpoint matches REBOUND over the same leg."""
    leg = flyby_leg
    rf, vf, _ = propagate_with_stm(
        leg["r0"], leg["v0"], leg["t0"], leg["t1"], ephem=leg["ephem"], moons=EGGIE_MOONS
    )
    cache = JovianRailsCache(EGGIE_MOONS, leg["ephem"], leg["t0"], leg["t1"])
    arc = JovianRestrictedNBody().propagate(
        leg["r0"], leg["v0"], leg["t0"], leg["t1"], moons=EGGIE_MOONS, cache=cache
    )
    assert arc.converged
    rel_r = float(np.linalg.norm(rf - arc.r_km) / np.linalg.norm(arc.r_km))
    rel_v = float(np.linalg.norm(vf - arc.v_km_s) / np.linalg.norm(arc.v_km_s))
    print(f"REBOUND-endpoint parity: rel_r={rel_r:.3e} rel_v={rel_v:.3e}")
    assert rel_r < 1e-4, f"position parity vs REBOUND {rel_r:.3e} > 1e-4"
    assert rel_v < 1e-4, f"velocity parity vs REBOUND {rel_v:.3e} > 1e-4"
