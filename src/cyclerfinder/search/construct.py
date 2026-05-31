"""Patched-conic cycler construction.

Given a flyby sequence and a list of encounter epochs, produces a
:class:`~cyclerfinder.model.Cycler` by Lambert-solving each heliocentric
leg against an :class:`~cyclerfinder.core.ephemeris.Ephemeris`.

Scope (M3)
----------

**No optimisation.** The constructor consumes encounter times as inputs.
Searching encounter times to minimise closure residual is an M5
deliverable; searching sequences is M4. M3 only needs the deterministic
"given a schedule, compute the trajectory" forward map.

The constructor selects a single Lambert solution per leg from
``(n_revs, branch)`` arguments — it does not multiplex over branches.
This is the M3-appropriate behaviour for the Aldrin reproduction and
the 2-synodic E-M-E sanity check; M4 will exercise multi-branch
enumeration once a search is added.

References
----------

* Spec §5 step 3 (seed construction), §9 (Aldrin validation anchors).
* Plan: ``docs/phases/m3-model-construct/plan.md`` §3.3.
"""

from __future__ import annotations

from math import pi

import numpy as np

from cyclerfinder.core.constants import MU_SUN_KM3_S2, PLANETS, SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.lambert import lambert
from cyclerfinder.model import Cycler, Encounter, Leg


def construct_cycler(
    sequence: list[str],
    encounter_times_sec: list[float],
    ephem: Ephemeris,
    mu_sun: float = MU_SUN_KM3_S2,
    max_revs_per_leg: list[int] | None = None,
    branch_per_leg: list[str] | None = None,
) -> Cycler:
    """Build a :class:`Cycler` from a sequence and encounter epochs.

    Parameters
    ----------
    sequence:
        Ordered one-letter body codes (e.g. ``["E", "M", "E"]``). Each must
        be a key in :data:`cyclerfinder.core.constants.PLANETS`.
    encounter_times_sec:
        Encounter epochs in seconds from ``t = 0``. Strictly increasing;
        ``len == len(sequence) >= 2``.
    ephem:
        Planet-state provider. M3 uses ``model="circular"``.
    mu_sun:
        Solar gravitational parameter (km^3/s^2). Defaults to
        :data:`MU_SUN_KM3_S2`.
    max_revs_per_leg:
        Per-leg maximum revolution count. Defaults to ``[0] * (n - 1)``,
        i.e. direct (0-rev) Lambert solutions for every leg.
    branch_per_leg:
        Per-leg branch selection label in ``{"single", "low", "high"}``.
        Defaults to ``["single"] * (n - 1)``. For ``n_revs = 0`` only the
        ``"single"`` branch exists; the M1 Lambert solver returns it under
        that label.

    Returns
    -------
    Cycler
        ``period = encounter_times_sec[-1] - encounter_times_sec[0]``.
        Encounter ``vinf_in`` / ``vinf_out`` follow the open-sequence
        convention: at the first encounter ``vinf_in = vinf_out`` (no
        preceding leg), at the last encounter ``vinf_out = vinf_in`` (no
        succeeding leg); intermediate encounters carry the actual transfer
        :math:`V_\\infty` vectors from the two adjacent legs.

    Raises
    ------
    ValueError
        On length mismatch, non-monotonic times, unknown body code, fewer
        than two encounters, or if no Lambert solution matches the
        requested branch.
    """
    n = len(sequence)
    if n < 2:
        raise ValueError(f"sequence must have at least 2 encounters, got {n}")
    if len(encounter_times_sec) != n:
        raise ValueError(
            f"len(encounter_times_sec)={len(encounter_times_sec)} must equal len(sequence)={n}"
        )
    # Strict monotonicity of epochs.
    for i in range(n - 1):
        if encounter_times_sec[i + 1] <= encounter_times_sec[i]:
            raise ValueError(
                f"encounter_times_sec must be strictly increasing; "
                f"violation at index {i}: "
                f"{encounter_times_sec[i]} -> {encounter_times_sec[i + 1]}"
            )
    # Body codes.
    for i, body in enumerate(sequence):
        if body not in PLANETS:
            raise ValueError(
                f"unknown body code {body!r} at sequence index {i}; "
                f"valid codes: {sorted(PLANETS.keys())}"
            )

    n_legs = n - 1
    if max_revs_per_leg is None:
        max_revs_per_leg = [0] * n_legs
    if branch_per_leg is None:
        branch_per_leg = ["single"] * n_legs
    if len(max_revs_per_leg) != n_legs:
        raise ValueError(
            f"len(max_revs_per_leg)={len(max_revs_per_leg)} must equal n_legs={n_legs}"
        )
    if len(branch_per_leg) != n_legs:
        raise ValueError(f"len(branch_per_leg)={len(branch_per_leg)} must equal n_legs={n_legs}")

    # Resolve each encounter's planet state up front.
    planet_states: list[tuple[np.ndarray, np.ndarray]] = [
        ephem.state(sequence[i], encounter_times_sec[i]) for i in range(n)
    ]

    # Solve each leg.
    leg_vels: list[tuple[np.ndarray, np.ndarray, int, str]] = []  # (v_dep, v_arr, n_revs, branch)
    for j in range(n_legs):
        r_from, _ = planet_states[j]
        r_to, _ = planet_states[j + 1]
        tof = encounter_times_sec[j + 1] - encounter_times_sec[j]
        sols = lambert(
            r_from,
            r_to,
            tof,
            mu=mu_sun,
            prograde=True,
            max_revs=max_revs_per_leg[j],
        )
        requested_branch = branch_per_leg[j]
        # Pick the matching solution.
        chosen = None
        for sol in sols:
            if sol.branch == requested_branch:
                chosen = sol
                break
        if chosen is None:
            available = [(s.n_revs, s.branch) for s in sols]
            raise ValueError(
                f"no Lambert solution with branch={requested_branch!r} on leg {j} "
                f"({sequence[j]}->{sequence[j + 1]}); available={available}"
            )
        leg_vels.append((chosen.v1, chosen.v2, chosen.n_revs, chosen.branch))

    # Build Encounter list.
    encounters: list[Encounter] = []
    for i in range(n):
        r_p, v_p = planet_states[i]
        if i == 0:
            v_sc_dep = leg_vels[0][0]
            vinf_out = v_sc_dep - v_p
            vinf_in = vinf_out  # boundary convention: no preceding leg
        elif i == n - 1:
            v_sc_arr = leg_vels[-1][1]
            vinf_in = v_sc_arr - v_p
            vinf_out = vinf_in  # boundary convention: no succeeding leg
        else:
            v_sc_arr = leg_vels[i - 1][1]
            v_sc_dep = leg_vels[i][0]
            vinf_in = v_sc_arr - v_p
            vinf_out = v_sc_dep - v_p
        encounters.append(
            Encounter(
                body=sequence[i],
                t=encounter_times_sec[i],
                r=r_p,
                v_planet=v_p,
                vinf_in=vinf_in,
                vinf_out=vinf_out,
            )
        )

    # Build Leg list.
    legs: list[Leg] = []
    for j in range(n_legs):
        v_dep, v_arr, n_revs, branch = leg_vels[j]
        legs.append(
            Leg(
                from_body=sequence[j],
                to_body=sequence[j + 1],
                t_depart=encounter_times_sec[j],
                t_arrive=encounter_times_sec[j + 1],
                v_depart=v_dep,
                v_arrive=v_arr,
                n_revs=n_revs,
                branch=branch,
            )
        )

    period = encounter_times_sec[-1] - encounter_times_sec[0]
    return Cycler(bodies=list(sequence), period=period, encounters=encounters, legs=legs)


def build_aldrin_seed(
    ephem: Ephemeris,
    t_start_sec: float | None = None,
    transfer_angle_deg: float = 132.0,
    em_tof_days: float = 146.0,
) -> Cycler:
    """Build the canonical Aldrin cycler as a 2-encounter Earth -> Mars slice.

    The Aldrin Earth -> Mars leg is a heliocentric arc that **passes
    through perihelion** — the spacecraft departs Earth descending toward
    a 0.97 AU perihelion, sweeps back outward, and arrives at Mars's
    orbit (1.524 AU) on the other side of the Sun from where it started.
    The heliocentric angle between the Earth-departure and Mars-arrival
    radius vectors is therefore **~132°**, not the ~80° one would naively
    compute as the true-anomaly span on the Aldrin ellipse (a = 1.60 AU,
    e = 0.393) from r = 1.0 AU to r = 1.524 AU. The naive 80° transfer
    is the *short-way outbound* leg from peri-going-outward to apo-going-
    outward; the actual Aldrin leg goes the other way (through peri),
    which is the longer 132° arc.

    The corresponding ``t_dep`` is found by inverting the phase equation:

    .. math::

        \\lambda_M(t_{dep} + tof) - \\lambda_E(t_{dep}) &= \\Delta\\lambda \\\\
        n_M (t_{dep} + tof) - n_E t_{dep} &= \\Delta\\lambda \\\\
        t_{dep} &= \\frac{\\Delta\\lambda - n_M \\, tof}{n_M - n_E}

    For Δλ = 132°, n_M·146d ≈ 76.5°, this gives ``t_dep ≈ -120 d``. (The
    negative ``t_dep`` is benign — the M1 ephemeris has no calendar epoch,
    so any time origin is arbitrary.)

    A sweep of ``transfer_angle_deg`` from 60° to 180° (held out of the
    test surface) confirms that 132° hits the literature Aldrin elements
    (a, e, peri, apo, V∞_E, V∞_M) almost exactly: (1.60, 0.393, 0.97,
    2.23, 6.53, 9.74). Spec.md's (a = 1.659, e = 0.41) corresponds to
    roughly Δλ = 135°.

    Parameters
    ----------
    ephem:
        Planet-state provider; M3 uses ``Ephemeris(model="circular")``.
    t_start_sec:
        Earth encounter epoch in seconds. ``None`` (the default) computes
        the value that produces the requested transfer angle; pass an
        explicit value to override.
    transfer_angle_deg:
        Heliocentric longitude span Mars(arrival) minus Earth(departure),
        deg. Default 132° per the literature Aldrin elements.
    em_tof_days:
        Earth -> Mars leg time of flight, days. Default 146 (Aldrin
        published value).

    Returns
    -------
    Cycler
        ``bodies = ["E", "M"]``, 2 encounters, 1 leg, period
        ``em_tof_days * SECONDS_PER_DAY``.
    """
    tof_sec = em_tof_days * SECONDS_PER_DAY
    if t_start_sec is None:
        n_e_rad_per_s = PLANETS["E"].mean_motion_deg_day * (pi / 180.0) / SECONDS_PER_DAY
        n_m_rad_per_s = PLANETS["M"].mean_motion_deg_day * (pi / 180.0) / SECONDS_PER_DAY
        transfer_angle_rad = transfer_angle_deg * (pi / 180.0)
        t_start_sec = (transfer_angle_rad - n_m_rad_per_s * tof_sec) / (
            n_m_rad_per_s - n_e_rad_per_s
        )
    times = [t_start_sec, t_start_sec + tof_sec]
    return construct_cycler(
        sequence=["E", "M"],
        encounter_times_sec=times,
        ephem=ephem,
    )


__all__ = ["build_aldrin_seed", "construct_cycler"]
