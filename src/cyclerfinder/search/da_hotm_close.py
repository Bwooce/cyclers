"""Corrector micro-multistart closer for DA/HOTM enumerator candidates (#450).

The Taylor-map enumerator emits a COARSE fixed-point candidate that lands in the
corrector's neighbourhood (~3e-4 for strongly-unstable multi-rev orbits like
P5g'; the FD-Taylor floor, see the #450 decision note). The existing asymmetric
corrector ``correct_general_periodic`` has a tight isotropic convergence basin
(~1e-5 for P5g', condition ~3600) but follows the stable manifold via its STM
Newton, so a small cluster of starts around the coarse candidate reliably lands a
member inside the basin.

This module is the thin closer: a corrector micro-multistart over an expanding
cluster of starts around the candidate. It does NOT modify the corrector (additive
seam, design draft §4) -- it only chooses starting points.

Pure: numpy + ``search.cr3bp_general_periodic``.
"""

from __future__ import annotations

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.cr3bp_general_periodic import (
    GeneralPeriodicOrbit,
    correct_general_periodic,
)


def close_candidate(
    system: cr3bp.CR3BPSystem,
    x0: float,
    xdot0: float,
    *,
    c_target: float,
    period_guess: float,
    half_crossings: int,
    ydot0_sign: float = 1.0,
    tol: float = 1e-11,
    radii: tuple[float, ...] = (0.0, 1.2e-5, 4e-5, 1e-4),
    cluster: int = 7,
    dx_cap: float = 0.02,
    max_iter: int = 80,
) -> GeneralPeriodicOrbit | None:
    """Close a coarse enumerator candidate via a corrector micro-multistart.

    Tries the corrector at the candidate itself, then over expanding
    ``cluster x cluster`` grids of starts at each radius in ``radii``, returning
    the FIRST converged orbit (residual <= ``tol``). Returns ``None`` if no start
    closes -- a legitimate negative (the coarse candidate was not a true fixed
    point, or sits outside the corrector's manifold reach).

    The EXISTING ``correct_general_periodic`` is used unchanged.
    """
    for radius in radii:
        if radius == 0.0:
            offsets = [(0.0, 0.0)]
        else:
            grid = np.linspace(-radius, radius, cluster)
            offsets = [(float(sx), float(sy)) for sx in grid for sy in grid]
        for sx, sy in offsets:
            orbit = correct_general_periodic(
                system,
                x0 + sx,
                xdot0 + sy,
                c_target,
                period_guess,
                half_crossings=half_crossings,
                ydot0_sign=ydot0_sign,
                tol=tol,
                dx_cap=dx_cap,
                max_iter=max_iter,
            )
            if orbit.converged and orbit.residual <= tol:
                return orbit
    return None


__all__ = ["close_candidate"]
