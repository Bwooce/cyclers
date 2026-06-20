"""Seed adapter: Russell N-arc parent seed -> full-state shooter seed (#388).

Bridges the :class:`~cyclerfinder.search.narc_continuation.NarcSeed` (a per-leg
Russell generic-return parent description, with sourced ToFs and the parent's
resonant rev/branch shape) into the full-state multiple-shooting
:class:`~cyclerfinder.nbody.shooter.ShootingSeed` consumed by
:func:`~cyclerfinder.nbody.shooter.shoot`.

The bridge solves each leg's conic Lambert at the seed epochs (via
:func:`cyclerfinder.search.correct._vinf_nodes`) to obtain the per-encounter V_inf
vectors, then maps those into per-node full Cartesian states with
:func:`cyclerfinder.nbody.shooter.seed_from_conic`. The slack leg is the longest
seed leg (most slack to absorb the period pin), matching the
``correct.py`` / ``narc_continuation`` convention.

Pure aside from the ephemeris read. A degenerate conic seed (a Lambert pathology
at these epochs) re-raises as :class:`ValueError`, which the batch driver catches.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

import cyclerfinder.search.correct as correct
from cyclerfinder.core.lambert import LambertError
from cyclerfinder.nbody import shooter

if TYPE_CHECKING:
    from cyclerfinder.core.ephemeris import Ephemeris
    from cyclerfinder.search.narc_continuation import NarcSeed


def russell_shooting_seed(
    seed: NarcSeed,
    *,
    t0_sec: float,
    ephem: Ephemeris,
) -> shooter.ShootingSeed:
    """Build a full-state :class:`ShootingSeed` from a Russell :class:`NarcSeed`.

    The slack leg is the longest seed leg; the free ToFs (slack excluded) feed
    :func:`correct._vinf_nodes` (which re-inserts the slack leg internally), and
    the FULL per-leg ToF list feeds :func:`shooter.seed_from_conic`.

    Raises
    ------
    ValueError
        If the conic seed is degenerate at these epochs (Lambert pathology).
    """
    tofs = list(seed.tof_seed_days)
    slack_leg = int(np.argmax(tofs))
    free_tofs = [t for i, t in enumerate(tofs) if i != slack_leg]
    period_days = seed.period_sec / correct.DAY_S

    try:
        nodes = correct._vinf_nodes(
            sequence=seed.sequence,
            per_leg_revs=seed.per_leg_revs,
            per_leg_branch=seed.per_leg_branch,
            t0_sec=t0_sec,
            free_tof_days=free_tofs,
            slack_leg=slack_leg,
            period_days=period_days,
            ephem=ephem,
        )
    except (LambertError, ValueError) as exc:
        raise ValueError(f"degenerate conic seed for Russell parent: {exc}") from exc

    return shooter.seed_from_conic(
        sequence=seed.sequence,
        vinf_nodes=nodes,
        t0_sec=t0_sec,
        tofs_days=tofs,
        slack_leg=slack_leg,
        period_days=period_days,
        ephem=ephem,
    )


__all__ = ["russell_shooting_seed"]
