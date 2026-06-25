"""#449 Task 5: the Releg swap seam at the discovery moontour gate.

``RepeatedMoonTarget._close_one_phasing`` (``discovery_campaign``) already takes
an injected ``lambert`` callable (line 463) — the seam the releg genome swaps. A
``BallisticReleg``-backed adapter routed through that seam must reproduce today's
ballistic result BIT-FOR-BIT (the regression lock; risk-table row "the two swap
sites drift out of sync"). The default path stays ``core.lambert.lambert``, so
every existing discovery test is untouched.
"""

from __future__ import annotations

import math

from cyclerfinder.core.lambert import lambert
from cyclerfinder.search.discovery_campaign import RepeatedMoonTarget
from cyclerfinder.search.releg_solver import BallisticReleg, ballistic_lambert_adapter


def _fixed_phasing_inputs() -> tuple[
    RepeatedMoonTarget,
    tuple[str, ...],
    list[int],
    int,
    dict[str, float],
    float,
    dict[str, tuple[float, float]],
    float,
]:
    """A fixed (sequence, n_rev, phasing, tof-scale) Jovian skeleton."""
    target = RepeatedMoonTarget(primary="Jupiter", moons=("Io", "Europa", "Ganymede"))
    consts = target._moon_consts()
    from cyclerfinder.core.satellites import PRIMARIES

    mu = PRIMARIES["Jupiter"]
    seq = ("Io", "Europa", "Ganymede", "Io")
    nrevs = [1, 0, 1]
    n_legs = len(seq) - 1
    theta0 = {m: 0.3 * (j + 1) for j, m in enumerate(sorted(consts))}
    tof_scale = 1.5
    return target, seq, nrevs, n_legs, theta0, tof_scale, consts, mu


def test_ballistic_releg_swap_matches_baseline() -> None:
    """A BallisticReleg adapter at ``_close_one_phasing`` == the lambert baseline.

    Run the closure with (a) the current ``core.lambert.lambert`` callable and
    (b) the ``BallisticReleg``-backed adapter; assert identical
    ``(feasible, worst, vinf, tofs)`` bit-for-bit.
    """
    target, seq, nrevs, n_legs, theta0, tof_scale, consts, mu = _fixed_phasing_inputs()

    baseline = target._close_one_phasing(seq, nrevs, n_legs, theta0, tof_scale, consts, mu, lambert)
    adapter = ballistic_lambert_adapter(BallisticReleg())
    swapped = target._close_one_phasing(seq, nrevs, n_legs, theta0, tof_scale, consts, mu, adapter)

    feasible_b, worst_b, vinf_b, tofs_b = baseline
    feasible_s, worst_s, vinf_s, tofs_s = swapped
    assert feasible_s == feasible_b
    assert worst_s == worst_b  # bit-for-bit
    assert vinf_s == vinf_b
    assert tofs_s == tofs_b


def test_swap_seam_feasible_skeleton() -> None:
    """Sanity: the fixed skeleton actually closes (so the equality is non-trivial)."""
    target, seq, nrevs, n_legs, theta0, tof_scale, consts, mu = _fixed_phasing_inputs()
    feasible, worst, _vinf, tofs = target._close_one_phasing(
        seq, nrevs, n_legs, theta0, tof_scale, consts, mu, lambert
    )
    assert feasible is True
    assert math.isfinite(worst)
    assert len(tofs) == n_legs
