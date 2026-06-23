"""#430 Unified global MGA-DSM precursor search engine.

Replaces the precursor matcher's local-optimiser path with a global
differential_evolution search over (launch epoch, per-leg TOFs, per-leg DSM),
seeded with eccentric-body Tisserand-Poincaré candidates, ranking survivors by
dv_band / total ΔV. See docs/superpowers/specs/2026-06-23-global-precursor-engine-design.md.
"""

from __future__ import annotations

import numpy as np

from cyclerfinder.core.constants import AU_KM
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.tisserand_mga_window import MGAChainCandidate, find_mga_chains


def eccentric_tp_linkable_radius_au(body: str, t_sec: float, ephemeris: Ephemeris) -> float:
    """Body's ACTUAL heliocentric radius (AU) at ``t_sec``, for the eccentric
    Tisserand-Poincaré graph (Campagnola-Russell 2009 Part B). The T-P contour
    drawn at the real encounter radius shifts/widens the linkable set vs the
    mean-``a`` circular form. Reduces to ``sma_au`` on the circular backend."""
    r_km, _v = ephemeris.state(body, t_sec)
    return float(np.linalg.norm(np.asarray(r_km, dtype=np.float64))) / AU_KM


def eccentric_tp_seeds(
    *,
    first_body: str,
    seed_vinf_kms: float,
    launch_window: tuple[str, str],
    ephemeris: Ephemeris,
    intermediate_bodies: tuple[str, ...] = ("V", "E"),
    max_legs: int = 3,
    vinf_grid_kms: tuple[float, ...] = (4.0, 5.0, 6.0, 7.0, 8.0),
    tof_box_days_per_leg: tuple[float, float] = (80.0, 500.0),
    epoch_step_days: float = 60.0,
    vinf_terminal_tol_kms: float = 0.8,
) -> list[MGAChainCandidate]:
    """Enumerate Earth-launched MGA chains terminating at ``first_body`` near
    ``seed_vinf_kms``, ranked with eccentric-body (real-radius) Tisserand-Poincaré
    linkability. Returns the DE init population (MGAChainCandidate list).

    Notes
    -----
    The real :func:`find_mga_chains` signature (tisserand_mga_window.py:615)
    takes ``launch_window`` and ``planet_set`` as positional args, is pure
    geometry (it does NOT accept an ``ephemeris`` argument — the ephemeris is
    consumed later by the Phase-1 closure/validation functions), and returns a
    lazy ``Iterator[MGAChainCandidate]``. The ``ephemeris`` parameter here is
    retained for the eccentric-radius linkability (see
    :func:`eccentric_tp_linkable_radius_au`) and for downstream DE seeding; it
    is not forwarded into the pure-geometry enumerator.
    """
    candidates = find_mga_chains(
        launch_window,
        tuple(dict.fromkeys((first_body, *intermediate_bodies))),
        max_legs=max_legs,
        vinf_grid_kms=vinf_grid_kms,
        tof_box_days_per_leg=tof_box_days_per_leg,
        epoch_step_days=epoch_step_days,
    )
    out: list[MGAChainCandidate] = []
    for c in candidates:
        if c.sequence[-1] != first_body:
            continue
        if abs(c.vinf_tuple_kms[-1] - seed_vinf_kms) > vinf_terminal_tol_kms:
            continue
        out.append(c)
    return out
