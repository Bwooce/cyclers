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
    jovian_defect_residual,
)
from cyclerfinder.nbody.jovian_ideal import (
    EGGIE_MOONS,
    build_eggie_periapsis_seed,
    build_subarc_seed,
    ideal_ephemeris_from_guess,
    subarc_defect_residual,
)
from cyclerfinder.nbody.jovian_stm import (
    jovian_stm_jacobian,
    propagate_with_stm,
    subarc_stm_jacobian,
)
from cyclerfinder.nbody.shooter import (
    _fd_jacobian,
    _seed_with_states,
    _serial_columns,
    _states_to_x,
    _x_to_states,
)
from cyclerfinder.search.resonant_conic import EggieGuess, eggie_refined_guess

# Per-component FD probe scale (km for position, km/s for velocity) — the v->r STM
# block scales like the time-of-flight (~1e6 s), so velocity columns need a far
# smaller probe than position columns (mirrors tests/nbody/test_propagator_stm.py).
_PROBE_SCALE = np.array([1.0, 1.0, 1.0, 1e-3, 1e-3, 1e-3])


@pytest.fixture(scope="module")
def eggie_guess() -> EggieGuess:
    """The refined in-basin EGGIE conic guess (built once per module — ~30 s)."""
    return eggie_refined_guess()


@pytest.fixture(scope="module")
def flyby_leg(eggie_guess: EggieGuess) -> dict[str, Any]:
    """The EGGIE leg that starts closest to a Galilean moon (interior flyby leg).

    Returns the chosen leg's endpoints/epochs, the ideal ephemeris, and the
    closest-approach distance so the test can assert the gate is genuinely a
    near-flyby (within the moon SOI), per the memory rule.
    """
    guess = eggie_guess
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


def test_stm_jacobian_matches_fd_on_eggie_seed(eggie_guess: EggieGuess) -> None:
    """Stage-2 gate: the analytic block-bidiagonal Jacobian matches the FD oracle.

    :func:`jovian_stm_jacobian` (assembled from per-leg analytic STMs) must agree
    with a finite-difference Jacobian of :func:`jovian_defect_residual` on the EGGIE
    periapsis seed — per nonzero 6x6 block — to well within the plan's ~1e-4 rel.
    FD stays the parity oracle. This is the Jacobian the Stage-3 ``jacobian="stm"``
    corrector relies on (and ~40x cheaper to build than the FD oracle).
    """
    ephem = ideal_ephemeris_from_guess(eggie_guess)
    seed = build_eggie_periapsis_seed(eggie_guess)
    moons = EGGIE_MOONS
    n = len(seed.sequence)
    cache = JovianRailsCache(moons, ephem, min(seed.epochs), max(seed.epochs))  # type: ignore[arg-type]

    def residual_of_x(x: np.ndarray) -> np.ndarray:
        trial = _seed_with_states(seed, _x_to_states(x, n))
        return jovian_defect_residual(
            trial,
            ephem=ephem,  # type: ignore[arg-type]
            cache=cache,
            moons=moons,
            accuracy=1e-11,
        )

    x0 = _states_to_x(seed.node_states)
    f0 = residual_of_x(x0)
    fd = _fd_jacobian(residual_of_x, x0, f0, column_eval=_serial_columns)
    stm = jovian_stm_jacobian(seed, x0, ephem=ephem, moons=moons)
    assert stm.shape == fd.shape

    overall = float(np.linalg.norm(stm - fd) / np.linalg.norm(fd))
    print(f"STM vs FD overall rel = {overall:.3e}")
    assert overall < 1e-4, f"overall STM-vs-FD rel {overall:.3e} > 1e-4"

    # Per nonzero 6x6 block (the binding parity criterion in the plan).
    n_leg = (n - 1) * 6
    n_hinge = max(0, n - 2)
    wrap0 = n_leg + n_hinge

    def block_rel(rows: slice, cols: slice) -> float:
        b = fd[rows, cols]
        bn = float(np.linalg.norm(b))
        if bn == 0.0:
            return 0.0
        return float(np.linalg.norm(stm[rows, cols] - b) / bn)

    worst = 0.0
    for i in range(n - 1):
        rows = slice(i * 6, (i + 1) * 6)
        worst = max(worst, block_rel(rows, slice(i * 6, (i + 1) * 6)))
        worst = max(worst, block_rel(rows, slice((i + 1) * 6, (i + 2) * 6)))
    wrap = slice(wrap0, wrap0 + 6)
    worst = max(worst, block_rel(wrap, slice(0, 6)))
    worst = max(worst, block_rel(wrap, slice((n - 1) * 6, n * 6)))
    print(f"worst nonzero-block rel = {worst:.3e}")
    assert worst < 1e-4, f"worst nonzero-block STM-vs-FD rel {worst:.3e} > 1e-4"


@pytest.mark.parametrize("n_subarcs", [2, 3])
def test_subarc_stm_jacobian_matches_fd(eggie_guess: EggieGuess, n_subarcs: int) -> None:
    """Stage-4 gate: the sub-arc block-bidiagonal Jacobian matches the FD oracle.

    With ``n_subarcs - 1`` interior continuity nodes per leg, :func:`subarc_stm_jacobian`
    (one analytic STM per sub-arc) must agree per nonzero 6x6 block with a
    finite-difference Jacobian of :func:`subarc_defect_residual` on the sub-arc seed —
    to well within the plan's ~1e-4 rel. FD stays the parity oracle. This is the
    Jacobian the Stage-4 ``n_subarcs>1`` corrector relies on.
    """
    ephem = ideal_ephemeris_from_guess(eggie_guess)
    seed = build_eggie_periapsis_seed(eggie_guess)
    moons = EGGIE_MOONS
    sub = build_subarc_seed(seed, n_subarcs, ephem=ephem, moons=moons)
    m = len(sub.node_states)
    n_enc = len(sub.encounter_idx)
    cache = JovianRailsCache(moons, ephem, min(sub.epochs), max(sub.epochs))  # type: ignore[arg-type]

    def residual_of_x(x: np.ndarray) -> np.ndarray:
        return subarc_defect_residual(
            sub,
            _x_to_states(x, m),
            ephem=ephem,  # type: ignore[arg-type]
            cache=cache,
            moons=moons,
            accuracy=1e-11,
        )

    x0 = _states_to_x(sub.node_states)
    f0 = residual_of_x(x0)
    fd = _fd_jacobian(residual_of_x, x0, f0, column_eval=_serial_columns)
    stm = subarc_stm_jacobian(sub, x0, ephem=ephem, moons=moons)
    assert stm.shape == fd.shape

    overall = float(np.linalg.norm(stm - fd) / np.linalg.norm(fd))
    print(f"n_subarcs={n_subarcs} sub-arc STM vs FD overall rel = {overall:.3e}")
    assert overall < 1e-4, f"overall sub-arc STM-vs-FD rel {overall:.3e} > 1e-4"

    # Per nonzero 6x6 block (the binding parity criterion).
    n_leg = (m - 1) * 6
    n_hinge = max(0, n_enc - 2)
    wrap0 = n_leg + n_hinge

    def block_rel(rows: slice, cols: slice) -> float:
        b = fd[rows, cols]
        bn = float(np.linalg.norm(b))
        if bn == 0.0:
            return 0.0
        return float(np.linalg.norm(stm[rows, cols] - b) / bn)

    worst = 0.0
    for j in range(m - 1):  # sub-arc continuity blocks (Phi_j and -R_W coupling)
        rows = slice(j * 6, (j + 1) * 6)
        worst = max(worst, block_rel(rows, slice(j * 6, (j + 1) * 6)))
        worst = max(worst, block_rel(rows, slice((j + 1) * 6, (j + 2) * 6)))
    wrap = slice(wrap0, wrap0 + 6)
    i0 = sub.encounter_idx[0]
    i_last = sub.encounter_idx[-1]
    worst = max(worst, block_rel(wrap, slice(i0 * 6, (i0 + 1) * 6)))
    worst = max(worst, block_rel(wrap, slice(i_last * 6, (i_last + 1) * 6)))
    print(f"n_subarcs={n_subarcs} worst nonzero-block rel = {worst:.3e}")
    assert worst < 1e-4, f"worst sub-arc nonzero-block STM-vs-FD rel {worst:.3e} > 1e-4"
