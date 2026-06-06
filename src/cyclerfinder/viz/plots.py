# AP_FLAKE8_CLEAN
"""Plotting primitives — beat diagram, porkchop, trajectory.

Every entry point lazy-imports matplotlib via :func:`_require_mpl` (Agg backend,
no display) and writes a PNG to ``out_path``. matplotlib stays the optional
``[viz]`` extra; importing this module does not import matplotlib.
"""

from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

from cyclerfinder.core.constants import MU_SUN_KM3_S2

if TYPE_CHECKING:  # pragma: no cover - typing only
    from cyclerfinder.search.optimize import OptimisationResult


def _require_mpl() -> Any:
    """Lazy-import matplotlib (Agg) or raise :class:`MissingVizExtra`."""
    from cyclerfinder.viz import MissingVizExtra

    try:
        import matplotlib
    except ImportError as exc:  # pragma: no cover - exercised via monkeypatch
        raise MissingVizExtra(
            "matplotlib is required for plotting; install the viz extra: "
            "uv sync --extra viz   (or  pip install '.[viz]')"
        ) from exc
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def _iso_to_unix_sec(iso_date: str) -> float:
    """Convert an ISO date (``YYYY-MM-DD``) to POSIX seconds (UTC)."""
    dt = datetime.fromisoformat(iso_date)
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=UTC)
    return dt.timestamp()


# ---------------------------------------------------------------------------
# beat diagram
# ---------------------------------------------------------------------------


def beat_diagram(bodies: list[str], out_path: Path | str) -> None:
    """Plot pairwise synodic phase over the multi-body beat period.

    For each non-reference body the pair's synodic phase (fraction of a synodic
    cycle) is drawn across the beat period; the beat alignment is marked.
    """
    from cyclerfinder.search.resonance import (
        beat_period_days,
        multi_body_beat_days,
        synodic_period_days,
    )

    plt = _require_mpl()

    tuples = multi_body_beat_days(bodies)
    fig, ax = plt.subplots(figsize=(8, 4))
    if tuples:
        k_tuple = tuples[0]
        beat_days = beat_period_days(bodies, k_tuple)
    else:
        # No commensurability within tolerance — still draw the synodic curves
        # over an arbitrary 6-year window so the figure is informative.
        beat_days = 6.0 * 365.25
    t = np.linspace(0.0, beat_days, 400)

    ref = bodies[len(bodies) // 2] if len(bodies) >= 3 else bodies[0]
    for body in bodies:
        if body == ref:
            continue
        syn = synodic_period_days(body, ref)
        phase = (t / syn) % 1.0
        ax.plot(t / 365.25, phase, label=f"{body}-{ref} synodic phase")
    ax.axvline(beat_days / 365.25, color="k", linestyle="--", label="beat alignment")
    ax.set_xlabel("years")
    ax.set_ylabel("synodic phase (fraction)")
    ax.set_title(f"Beat diagram: {'-'.join(bodies)}")
    ax.legend(loc="upper right", fontsize="small")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


# ---------------------------------------------------------------------------
# porkchop
# ---------------------------------------------------------------------------


def porkchop(
    body_from: str,
    body_to: str,
    *,
    epoch_range: tuple[str, str],
    tof_range: tuple[float, float],
    out_path: Path | str,
    ephem: Any = None,
    n_epoch: int = 30,
    n_tof: int = 30,
) -> None:
    """Epoch x ToF grid of the departure C3 (||v_depart - v_planet||^2), as a contour.

    Uses the circular ephemeris by default (fast, deterministic). The astropy
    path is reachable by passing ``ephem=Ephemeris(model="astropy")``.
    """
    from cyclerfinder.core.ephemeris import Ephemeris
    from cyclerfinder.core.lambert import lambert

    plt = _require_mpl()
    ephem = ephem or Ephemeris(model="circular")

    t0 = _iso_to_unix_sec(epoch_range[0])
    t1 = _iso_to_unix_sec(epoch_range[1])
    epochs = np.linspace(t0, t1, n_epoch)
    tofs_days = np.linspace(tof_range[0], tof_range[1], n_tof)
    day_sec = 86400.0

    grid = np.full((n_tof, n_epoch), np.nan)
    for i, dep in enumerate(epochs):
        r1, v1 = ephem.state(body_from, float(dep))
        for j, tof_d in enumerate(tofs_days):
            tof = float(tof_d) * day_sec
            r2, _v2 = ephem.state(body_to, float(dep) + tof)
            try:
                sols = lambert(np.asarray(r1), np.asarray(r2), tof)
            except Exception:  # infeasible geometry leaves a blank cell
                continue
            if not sols:
                continue
            v_depart = np.asarray(sols[0].v1)
            c3 = float(np.linalg.norm(v_depart - np.asarray(v1)) ** 2)
            grid[j, i] = c3

    fig, ax = plt.subplots(figsize=(8, 5))
    years = (epochs - t0) / (365.25 * day_sec)
    masked = np.ma.masked_invalid(grid)
    cs = ax.contourf(years, tofs_days, masked, levels=20)
    fig.colorbar(cs, ax=ax, label="departure C3 (km^2/s^2)")
    ax.set_xlabel(f"years from {epoch_range[0]}")
    ax.set_ylabel("time of flight (days)")
    ax.set_title(f"Porkchop {body_from}->{body_to}")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


# ---------------------------------------------------------------------------
# trajectory
# ---------------------------------------------------------------------------


def trajectory(result: OptimisationResult, out_path: Path | str) -> None:
    """Heliocentric XY plot of a built cycler's legs plus planet circles."""
    from cyclerfinder.core.kepler import KeplerConvergenceError, propagate

    plt = _require_mpl()
    cycler = result.best_cycler

    fig, ax = plt.subplots(figsize=(6, 6))

    # Planet orbit circles (radius = encounter heliocentric distance).
    seen: set[str] = set()
    for enc in cycler.encounters:
        if enc.body in seen:
            continue
        seen.add(enc.body)
        radius = float(np.linalg.norm(np.asarray(enc.r)))
        theta = np.linspace(0, 2 * np.pi, 200)
        ax.plot(radius * np.cos(theta), radius * np.sin(theta), ":", alpha=0.4)

    # Each leg, sampled by two-body propagation from the departure state.
    for leg in cycler.legs:
        # Find the departure encounter (matching epoch) for r0.
        dep_enc = min(cycler.encounters, key=lambda e: abs(e.t - leg.t_depart))
        r0 = np.asarray(dep_enc.r, dtype=np.float64)
        v0 = np.asarray(leg.v_depart, dtype=np.float64)
        dt_total = leg.t_arrive - leg.t_depart
        samples = np.linspace(0.0, dt_total, 80)
        xs: list[float] = []
        ys: list[float] = []
        for dt in samples:
            try:
                r, _v = propagate(r0, v0, float(dt), mu=MU_SUN_KM3_S2)
            except KeplerConvergenceError:  # degenerate sample on a non-closed leg
                continue
            xs.append(float(r[0]))
            ys.append(float(r[1]))
        if xs:
            ax.plot(xs, ys, "-", label=f"{leg.from_body}->{leg.to_body}")

    # Encounter markers.
    for enc in cycler.encounters:
        r = np.asarray(enc.r)
        ax.plot(float(r[0]), float(r[1]), "o")

    ax.plot(0, 0, "y*", markersize=12)  # the Sun
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("x (km)")
    ax.set_ylabel("y (km)")
    ax.set_title(f"Trajectory: {result.cell.id}")
    ax.legend(loc="upper right", fontsize="small")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
