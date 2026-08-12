"""#830 — turn a multi-arc ballistic closure into a propagatable :class:`Cycler`.

Why this exists
---------------
``tests/search/test_free_return_v2_ballistic.py`` declined V2-ballistic promotion
for the Russell multi-arc rows on a STRUCTURAL ground:

    "the single ellipse does NOT represent the Earth-to-Earth resonant phasing
     intervals ... There is therefore no continuous >=3-lap trajectory to
     propagate for these objects — the V2-ballistic gate is structurally
     inapplicable to a single-arc slice of a multi-arc cycler."

#820's re-posed genome tiles the FULL cycler period (designated Mars-transit leg
plus every E-E phasing loop), so a continuous multi-arc trajectory now exists —
but nothing could feed it to the §12 / §14-V2 machinery, which consumes a
:class:`~cyclerfinder.model.cycler.Cycler`. The only builder in the tree
(:func:`cyclerfinder.search.free_return_v1.build_free_return_cycler`) emits the
3-encounter single-ellipse form. This module supplies the missing adapter.

The output is the OPEN form documented on :class:`Cycler`: ``n+1`` encounters and
``n`` legs, where the last encounter is the periodic wrap of the first. Per the
:class:`~cyclerfinder.model.cycler.Encounter` boundary convention the wrap
endpoints carry ``vinf_in == vinf_out`` so
:meth:`~cyclerfinder.model.cycler.Cycler.maintenance_dv` does not charge a
spurious discontinuity there.

Nothing is imposed: the encounter epochs, planet states and per-leg Lambert
solutions are re-derived from the converged ``(t0, tof_days)`` with the same
genome the corrector ran, so the reconstruction is the closure itself.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from cyclerfinder.core.constants import MU_SUN_KM3_S2
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.lambert import lambert
from cyclerfinder.model.cycler import Cycler, Encounter, Leg, SenseT
from cyclerfinder.search.correct import BallisticClosureResult, _pick

DAY_S = 86400.0


def build_multiarc_cycler(
    result: BallisticClosureResult,
    *,
    sequence: tuple[str, ...],
    per_leg_revs: tuple[int, ...],
    per_leg_branch: tuple[str, ...],
    period_sec: float,
    ephem: Ephemeris,
    sense: SenseT = "n/a",
    mu: float = MU_SUN_KM3_S2,
) -> Cycler:
    """Reconstruct an N-arc :class:`Cycler` from a converged ballistic closure.

    Parameters
    ----------
    result:
        A converged :class:`~cyclerfinder.search.correct.BallisticClosureResult`;
        its ``t0_sec`` and FULL ``tof_days`` (slack leg already re-inserted) fix
        the encounter epochs.
    sequence, per_leg_revs, per_leg_branch:
        The genome the corrector ran — the same Lambert ``(n_revs, branch)``
        selection must be reproduced or the arcs are not the closure's arcs.
    period_sec:
        The cycler period the closure was constrained to (``Cycler.period``).
    sense:
        Spec §16.2 direction tag; pass the catalogue row's own value when known.

    Raises
    ------
    ValueError
        If the ToF count does not match the sequence.
    """
    tofs: Sequence[float] = result.tof_days
    n_legs = len(sequence) - 1
    if len(tofs) != n_legs:
        raise ValueError(f"expected {n_legs} leg ToFs for {sequence}, got {len(tofs)}")

    epochs = [float(result.t0_sec)]
    for tof in tofs:
        epochs.append(epochs[-1] + float(tof) * DAY_S)

    states = [ephem.state(body, t) for body, t in zip(sequence, epochs, strict=True)]
    r_km = [np.asarray(r, dtype=np.float64) for r, _ in states]
    v_planet = [np.asarray(v, dtype=np.float64) for _, v in states]

    v_depart: list[np.ndarray] = []
    v_arrive: list[np.ndarray] = []
    for i in range(n_legs):
        sols = lambert(
            r_km[i], r_km[i + 1], float(tofs[i]) * DAY_S, mu=mu, max_revs=per_leg_revs[i]
        )
        sol = _pick(sols, per_leg_revs[i], per_leg_branch[i])
        v_depart.append(np.asarray(sol.v1, dtype=np.float64))
        v_arrive.append(np.asarray(sol.v2, dtype=np.float64))

    encounters: list[Encounter] = []
    for i, body in enumerate(sequence):
        # Boundary convention (Encounter docstring): the open sequence's first
        # and last encounters carry vinf_in == vinf_out. They are the SAME
        # physical flyby one period apart, so charging their difference would be
        # charging the frame rotation, not a manoeuvre.
        out = v_depart[i] - v_planet[i] if i < n_legs else v_arrive[i - 1] - v_planet[i]
        inn = v_arrive[i - 1] - v_planet[i] if i > 0 else v_depart[i] - v_planet[i]
        encounters.append(Encounter(body, epochs[i], r_km[i], v_planet[i], inn, out))

    legs = [
        Leg(
            sequence[i],
            sequence[i + 1],
            epochs[i],
            epochs[i + 1],
            v_depart[i],
            v_arrive[i],
            per_leg_revs[i],
            per_leg_branch[i],
        )
        for i in range(n_legs)
    ]
    return Cycler(list(sequence), float(period_sec), encounters, legs, sense=sense)
