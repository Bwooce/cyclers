"""#318 Phase 2b — Sobol sampler for the Jupiter Galilean short-cycler joint search.

Generates a low-discrepancy Sobol sample over the four #318 joint-search axes for a
closed Galilean-moon tour, returning :class:`~cyclerfinder.search.joint_cell.JointCell`
objects ready for :func:`~cyclerfinder.search.joint_cell.evaluate_joint_cell`.

Dimension layout for ``n_legs = len(sequence) - 1`` legs::

    dim 0:              epoch (continuous → JD in epoch_window → TDB ISO)
    dim 1 .. n_legs:    ToF per leg (continuous → days in tof_seed_range)
    dim n_legs+1 .. 2*n_legs:  n_revs per leg (rounded integer in n_revs_range)
    dim 2*n_legs+1 .. 3*n_legs: branch bit (< 0.5 → "high", ≥ 0.5 → "low")

Total dims: ``1 + 3 * n_legs`` (13 for CGCEC's 4 legs). Use ``n_samples`` that is a
power of 2 for best Sobol uniformity (e.g. 256 = 2^8 for the smoke run).
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.search.joint_cell import JointCell

__all__ = ["make_sobol_cells"]


def make_sobol_cells(
    n_samples: int,
    sequence: list[str],
    epoch_window: tuple[str, str],
    n_revs_range: tuple[int, int],
    tof_seed_range: tuple[float, float],
    powered_min_alt_km: float = 50.0,
    seed: int = 0,
) -> list[JointCell]:
    """Generate Sobol-sampled :class:`JointCell` objects for a Galilean-moon joint search.

    Parameters
    ----------
    n_samples:
        Number of cells (power of 2 gives best Sobol uniformity; 256 for smoke).
    sequence:
        Closed encounter sequence, ``sequence[0] == sequence[-1]``
        (e.g. ``list(jovian.CGCEC)``).
    epoch_window:
        ``(start_iso, end_iso)`` TDB ISO dates bounding the departure epoch.
    n_revs_range:
        ``(n_min, n_max)`` inclusive integer range for the per-leg Lambert rev count.
    tof_seed_range:
        ``(tof_min_days, tof_max_days)`` for the per-leg ToF seed.
    powered_min_alt_km:
        Flyby altitude floor forwarded to each cell's powered-maneuver model.
    seed:
        Sobol scramble RNG seed for reproducibility.

    Returns
    -------
    list[JointCell]
        ``n_samples`` cells with ``primary="Jupiter"``.
    """
    from astropy.time import Time
    from scipy.stats.qmc import Sobol

    seq = list(sequence)
    n_legs = len(seq) - 1
    if n_legs < 1:
        raise ValueError(f"sequence must have ≥ 2 elements; got {len(seq)}")
    if seq[0] != seq[-1]:
        raise ValueError(f"sequence must be closed (seq[0]={seq[0]!r} ≠ seq[-1]={seq[-1]!r})")

    n_dims = 1 + 3 * n_legs
    sampler = Sobol(d=n_dims, scramble=True, seed=seed)
    samples: NDArray[np.float64] = sampler.random(n_samples)

    jd_lo = float(Time(epoch_window[0], scale="tdb").jd)
    jd_hi = float(Time(epoch_window[1], scale="tdb").jd)
    tof_lo = float(tof_seed_range[0])
    tof_hi = float(tof_seed_range[1])
    n_rev_lo = int(n_revs_range[0])
    n_rev_hi = int(n_revs_range[1])

    cells: list[JointCell] = []
    for row in samples:
        jd = jd_lo + row[0] * (jd_hi - jd_lo)
        epoch_iso: str = Time(jd, format="jd", scale="tdb").isot

        tof_days = tuple(float(tof_lo + row[1 + k] * (tof_hi - tof_lo)) for k in range(n_legs))

        n_revs = tuple(
            int(
                np.clip(
                    int(np.floor(row[1 + n_legs + k] * (n_rev_hi - n_rev_lo + 1))) + n_rev_lo,
                    n_rev_lo,
                    n_rev_hi,
                )
            )
            for k in range(n_legs)
        )

        branches = tuple("low" if row[1 + 2 * n_legs + k] >= 0.5 else "high" for k in range(n_legs))

        cells.append(
            JointCell(
                primary="Jupiter",
                sequence=tuple(seq),
                epoch_iso=epoch_iso,
                n_revs=n_revs,
                branches=branches,
                tof_seed_days=tof_days,
                powered_min_alt_km=powered_min_alt_km,
            )
        )

    return cells
