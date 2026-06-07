"""Export numerically-propagated cycler trajectories as the site's sampled JSON.

This is the missing producer half of viz-2c (task #146). The site
(``cyclers.space``) already consumes a ``SampledTrajectory`` shape
(``src/lib/three-types.ts``): parallel ``timesSec`` (seconds) / ``positionsAU``
(AU, frame ``"eclipJ2000"``) arrays plus honesty strings ``fidelity`` /
``provenance``. Until now the ONLY producer was the synthetic dev fixture; this
script emits real, propagated geometry instead.

Pipeline (the cleanest real source we have)
-------------------------------------------
1. :func:`cyclerfinder.search.maintain.optimise_aldrin_maintenance_dv` solves the
   in-family powered Aldrin E->M->E cycler at its real DE440 launch window (the
   #134 machinery; the same call the §14 V2-powered gate runs per cycle). It
   returns a :class:`~cyclerfinder.model.Cycler` whose encounters sit on the real
   planet positions.
2. :func:`cyclerfinder.verify.propagate.multi_lap_propagation` propagates that
   cycler continuously across ``n_laps`` consecutive laps, reconstructing each
   lap's geometry from the encoded leg template rotated to the lap epoch. It
   returns dense ``(t, x, y, z, vx, vy, vz)`` rows — heliocentric INERTIAL
   ECLIPTIC-J2000, km and km/s, with ``t`` in TDB seconds since J2000.
3. We decimate the dense samples (curvature-adaptive Douglas-Peucker, seam-aware;
   see :func:`decimate_curvature`), convert km -> AU (the ONLY unit change;
   positions are already in the site's ``eclipJ2000`` frame so NO rotation is
   applied — see "Frame note" below), and emit the ``SampledTrajectory`` JSON the
   site fetches lazily on "View in 3D".

Frame note (the convert.py trap, NOT triggered here)
----------------------------------------------------
``cyclerfinder.nbody.convert`` exists to rotate ICRS-equatorial vectors (what raw
SPICE / external n-body tools return) into the J2000-ecliptic frame, and to map
``t_sec`` to SPICE ET. We do NOT use it here: the project's own
``Ephemeris("astropy")`` backend already returns heliocentric J2000-ECLIPTIC
states (``core/ephemeris.py``: "+z along the ecliptic north pole"), and
``multi_lap_propagation`` Kepler-propagates in that same frame. So the propagated
positions are already in the site's ``eclipJ2000`` frame and need only km->AU.
The time axis is already TDB seconds since J2000 — identical to the site clock's
"days since J2000 = timesSec / 86400" convention (and to ``convert``'s documented
``t_sec = 0`` = J2000 + 64.184 s ET identity), so timesSec is emitted verbatim.

Decimation honesty (curvature-adaptive)
---------------------------------------
We DECIMATE (curvature-adaptive Douglas-Peucker to a target tolerance, capped at
``--max-points``), never resample/smooth: every emitted point is an actual
propagated state, so the polyline introduces no geometry the integrator did not
produce (the site's clock then LINEARLY interpolates between them — it too invents
no curvature). Douglas-Peucker keeps points where the arc bends (perihelion,
flybys) and drops them on the near-straight stretches, so a tight chord error is
achieved with far fewer points than a uniform stride would need on this
eccentric, multi-lap path. The script reports the max chord error of the
decimated polyline vs the dense reference (the largest perpendicular distance from
a dropped dense point to the decimated segment that spans it), so the fidelity
claim is bounded and measured, not asserted.

Per-leg seams (a property of THIS source, reported honestly)
------------------------------------------------------------
``multi_lap_propagation`` reconstructs each leg independently from the REAL planet
position at the lap-shifted departure epoch (plus the rotated encoded V-infinity),
rather than numerically stitching one continuous integration. The legs therefore
meet C0-discontinuously at the encounters/lap boundaries: the craft "snaps" to the
planet's true position at each encounter. These seam jumps are STRUCTURE, not
decimation error, so the decimation tolerance is enforced WITHIN continuous runs
(split at seams) and the seams are reported separately (``seam_count`` /
``max_seam_jump_au`` in ``_meta``). The fidelity string names this honestly
("kepler-legs ... per-leg reconstruction") so the curve is never read as a single
smooth integrated orbit.

Usage
-----
    uv run python scripts/export_sampled_trajectories.py            # write JSON
    uv run python scripts/export_sampled_trajectories.py --check    # verify on disk

Writes ``data/sampled/<id>.json`` in this repo (the canonical source, mirroring
``scripts/emit-planet-elements.py``); the site copies it into its ``public/data``
so it is fetched at ``/data/sampled/<id>.json`` only when the 3D view opens.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import AU_KM
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.model.cycler import Cycler
from cyclerfinder.search.maintain import optimise_aldrin_maintenance_dv
from cyclerfinder.verify.propagate import multi_lap_propagation

# The Aldrin classic E-M cycler literature/priority epoch (matches the catalogue
# `aldrin-classic-em-k1-outbound` priority_date and the §14 V2-powered gate). The
# real-window launch-phase seed that places the solve in-family on DE440.
_ALDRIN_PRIORITY_DATE = datetime(1985, 10, 28, tzinfo=UTC)

# Dense propagation resolution per lap before decimation. High enough that the
# decimated polyline's chord error is dominated by the decimation stride, not by
# the dense grid itself (the dense grid IS the reference we measure error against).
_DENSE_SAMPLES_PER_LAP = 600

# Default exhibit: the V2-powered outbound row, the first real sampled exhibit.
_DEFAULT_N_LAPS = 3
_DEFAULT_MAX_POINTS = 500

# A consecutive-sample position jump larger than this (AU) is a per-leg
# reconstruction SEAM (the craft snapping to the planet at an encounter / lap
# boundary), not a smooth arc step. The decimator splits continuous runs at
# seams so the chord-error tolerance is not polluted by structural jumps.
_SEAM_JUMP_AU = 0.30
_SEAM_JUMP_KM = _SEAM_JUMP_AU * AU_KM

# Curvature-adaptive (Douglas-Peucker) chord tolerance, km. ~7.5e-4 AU — well
# under the site's documented 5e-3 AU coincidence bound — so the decimated
# polyline is visually coincident with the dense reference at the drawn line
# widths. Points are dropped only where the arc stays within this band of the
# chord; the point cap (`--max-points`) is the hard ceiling if the tolerance
# would keep more.
_DP_TOLERANCE_KM = 7.5e-4 * AU_KM


def _git_commit() -> str:
    """Short HEAD commit SHA of this repo, for the provenance string.

    Falls back to ``"unknown"`` if git is unavailable (e.g. an exported tree) so
    the export never fails on provenance alone.
    """
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=Path(__file__).resolve().parent.parent,
            capture_output=True,
            text=True,
            check=True,
        )
        return out.stdout.strip() or "unknown"
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


def propagate_aldrin_powered(
    *,
    n_laps: int = _DEFAULT_N_LAPS,
    dense_samples_per_lap: int = _DENSE_SAMPLES_PER_LAP,
) -> tuple[NDArray[np.float64], Cycler]:
    """Propagate the in-family powered Aldrin cycler over ``n_laps`` laps (dense).

    Solves the real-DE440 in-family maintenance cycle at the Aldrin launch window,
    then propagates it continuously. Returns the dense ``(N, 7)`` sample array
    ``[t_sec, x_km, y_km, z_km, vx, vy, vz]`` (heliocentric inertial
    ecliptic-J2000) and the optimised cycler.

    This is SLOW (real DE440 BVP solve, minutes) — it is the same machinery the
    §14 V2-powered gate runs, deliberately reusing the production solver rather
    than a faster surrogate so the exported geometry is the real exhibit.
    """
    ephem = Ephemeris(model="astropy")
    res = optimise_aldrin_maintenance_dv(ephem, real_window_priority_date=_ALDRIN_PRIORITY_DATE)
    cycler = res.cycler
    t_start = float(cycler.encounters[0].t)
    mlp = multi_lap_propagation(
        cycler,
        ephem,
        n_laps,
        t_start=t_start,
        n_samples_per_lap=dense_samples_per_lap,
    )
    return np.asarray(mlp["samples"], dtype=np.float64), cycler


def _perp_error(a: NDArray[np.float64], b: NDArray[np.float64], p: NDArray[np.float64]) -> float:
    """Perpendicular distance (clamped to the segment) of point ``p`` to ``a-b``."""
    ab = b - a
    ab_len2 = float(np.dot(ab, ab))
    if ab_len2 == 0.0:
        return float(np.linalg.norm(p - a))
    tproj = min(1.0, max(0.0, float(np.dot(p - a, ab)) / ab_len2))
    return float(np.linalg.norm(p - (a + tproj * ab)))


def _seam_split_indices(pos: NDArray[np.float64], seam_jump_km: float) -> list[int]:
    """Indices where a consecutive-sample jump exceeds ``seam_jump_km`` — the
    per-leg reconstruction seams. Returned as the START index of each big jump
    (so the run boundary is between ``i`` and ``i + 1``)."""
    return [
        i
        for i in range(pos.shape[0] - 1)
        if float(np.linalg.norm(pos[i + 1] - pos[i])) > seam_jump_km
    ]


def _douglas_peucker(pos: NDArray[np.float64], lo: int, hi: int, tol_km: float) -> list[int]:
    """Douglas-Peucker on ``pos[lo..hi]`` (inclusive). Returns kept indices.

    Iterative (explicit stack) to avoid recursion limits on long dense runs.
    Always keeps the two endpoints; recursively keeps the farthest point of any
    sub-span whose max perpendicular error exceeds ``tol_km``.
    """
    if hi <= lo + 1:
        return [lo, hi] if hi > lo else [lo]
    keep = {lo, hi}
    stack = [(lo, hi)]
    while stack:
        a_i, b_i = stack.pop()
        if b_i <= a_i + 1:
            continue
        a, b = pos[a_i], pos[b_i]
        worst_err = -1.0
        worst_j = -1
        for j in range(a_i + 1, b_i):
            err = _perp_error(a, b, pos[j])
            if err > worst_err:
                worst_err, worst_j = err, j
        if worst_err > tol_km and worst_j != -1:
            keep.add(worst_j)
            stack.append((a_i, worst_j))
            stack.append((worst_j, b_i))
    return sorted(keep)


def decimate_curvature(
    dense: NDArray[np.float64],
    max_points: int,
    *,
    tol_km: float = _DP_TOLERANCE_KM,
    seam_jump_km: float = _SEAM_JUMP_KM,
) -> tuple[NDArray[np.float64], float, dict[str, Any]]:
    """Curvature-adaptive (Douglas-Peucker) decimation of dense ``(N, 7)`` rows.

    Splits the dense series into continuous runs at the per-leg reconstruction
    seams (consecutive jumps > ``seam_jump_km``), runs Douglas-Peucker on each run
    with chord tolerance ``tol_km``, then concatenates the kept rows preserving
    time order (seam endpoints are always kept, so the seams survive intact). If
    the kept set still exceeds ``max_points``, the tolerance is doubled and the
    pass repeated (so the point cap is the hard ceiling).

    Returns ``(decimated_rows, max_chord_error_km, meta)`` where ``meta`` carries
    the seam diagnostics. The chord error is measured WITHIN runs only (seam jumps
    are structure, not decimation error).

    Decimation, not resampling: every kept row is an actual propagated state.
    """
    if max_points < 2:
        raise ValueError(f"max_points must be >= 2; got {max_points}")
    n = dense.shape[0]
    pos = dense[:, 1:4]
    seams = _seam_split_indices(pos, seam_jump_km)
    # Run boundaries: [0 .. s0], [s0+1 .. s1], ... [sk+1 .. n-1].
    bounds: list[tuple[int, int]] = []
    start = 0
    for s in seams:
        bounds.append((start, s))
        start = s + 1
    bounds.append((start, n - 1))

    cur_tol = tol_km
    keep: list[int] = []
    for _ in range(12):  # bounded tolerance escalation to honour the point cap
        keep = []
        for lo, hi in bounds:
            run_keep = _douglas_peucker(pos, lo, hi, cur_tol)
            if keep and run_keep and run_keep[0] == keep[-1]:
                run_keep = run_keep[1:]
            keep.extend(run_keep)
        if len(keep) <= max_points:
            break
        cur_tol *= 2.0

    # Hard ceiling: if escalation still overshoots (e.g. a near-perfect circle
    # where DP keeps points in big jumps), uniformly thin the kept indices down
    # to max_points, always retaining the seam endpoints and the global ends so
    # no seam is erased. Still pure decimation — every survivor is a real state.
    if len(keep) > max_points:
        protected = {0, n - 1}
        for s in seams:
            protected.add(s)
            protected.add(s + 1)
        thin_stride = int(np.ceil(len(keep) / max_points))
        thinned = [idx for k, idx in enumerate(keep) if k % thin_stride == 0]
        thinned_set = set(thinned) | protected
        keep = sorted(thinned_set)

    decimated = dense[keep, :]

    # Max within-run perpendicular error of any DROPPED point vs the kept segment
    # that brackets it. Seams are skipped (a kept segment never spans a seam).
    keep_set = set(keep)
    seam_set = set(seams)
    max_err = 0.0
    kept_sorted = keep
    for s_i in range(len(kept_sorted) - 1):
        i0 = kept_sorted[s_i]
        i1 = kept_sorted[s_i + 1]
        if any(j in seam_set for j in range(i0, i1)):
            continue  # this kept segment straddles a seam: not a decimation span
        a, b = pos[i0], pos[i1]
        for j in range(i0 + 1, i1):
            if j in keep_set:
                continue
            err = _perp_error(a, b, pos[j])
            if err > max_err:
                max_err = err

    seam_jumps_km = [float(np.linalg.norm(pos[s + 1] - pos[s])) for s in seams]
    meta = {
        "tol_km_used": cur_tol,
        "seam_count": len(seams),
        "max_seam_jump_km": max(seam_jumps_km) if seam_jumps_km else 0.0,
        "max_seam_jump_au": (max(seam_jumps_km) / AU_KM) if seam_jumps_km else 0.0,
    }
    return decimated, max_err, meta


def _strictly_increasing(times: list[float], eps: float = 1.0e-3) -> list[float]:
    """Enforce strictly-increasing times (the site's bracketing-search
    precondition) by nudging any non-increasing time forward by ``eps`` seconds.

    The per-leg seams put the arriving-leg endpoint and the next leg's start at
    the SAME encounter instant; nudging the second by 1 ms preserves BOTH sampled
    positions (so the seam jump still renders) while satisfying the contract. The
    nudge is ~1e-3 s against a multi-year span — far below any visible time grain.
    """
    out: list[float] = []
    prev = -float("inf")
    for t in times:
        if t <= prev:
            t = prev + eps
        out.append(t)
        prev = t
    return out


def build_payload(
    *,
    entry_id: str,
    n_laps: int = _DEFAULT_N_LAPS,
    max_points: int = _DEFAULT_MAX_POINTS,
    dense_samples_per_lap: int = _DENSE_SAMPLES_PER_LAP,
) -> dict[str, Any]:
    """Assemble the ``SampledTrajectory`` payload for ``entry_id`` (Aldrin only).

    Currently only the V2-powered ``aldrin-classic-em-k1-outbound`` row is
    supported (the first real exhibit). Returns the exact JSON object shape the
    site's ``SampledTrajectory`` interface consumes.
    """
    if entry_id != "aldrin-classic-em-k1-outbound":
        raise ValueError(
            f"unsupported entry_id {entry_id!r}; only "
            "'aldrin-classic-em-k1-outbound' is wired today"
        )

    dense, cycler = propagate_aldrin_powered(
        n_laps=n_laps, dense_samples_per_lap=dense_samples_per_lap
    )
    decimated, chord_err_km, dmeta = decimate_curvature(dense, max_points)

    times_sec = _strictly_increasing([float(r[0]) for r in decimated])
    positions_au = [
        [float(r[1]) / AU_KM, float(r[2]) / AU_KM, float(r[3]) / AU_KM] for r in decimated
    ]
    chord_err_au = chord_err_km / AU_KM
    seam_count = int(dmeta["seam_count"])
    commit = _git_commit()

    return {
        "kind": "sampled",
        "timesSec": times_sec,
        "positionsAU": positions_au,
        "frame": "eclipJ2000",
        "fidelity": (
            f"kepler-legs + powered-flyby retarget over {n_laps} laps "
            "(DE440 encounters; per-leg reconstruction, C0 seams at flybys)"
        ),
        "provenance": (
            "cyclerfinder.verify.propagate.multi_lap_propagation of the in-family "
            "powered Aldrin E-M-E cycler "
            "(search.maintain.optimise_aldrin_maintenance_dv, real DE440 window "
            f"@1985-10-28); commit {commit}; {len(decimated)} pts (Douglas-Peucker) "
            f"from {dense.shape[0]} dense, max chord error {chord_err_au:.2e} AU "
            f"({chord_err_km:.0f} km), {seam_count} per-leg seams"
        ),
        # Diagnostics (ignored by the site's SampledTrajectory consumer; recorded
        # so the export's own honesty numbers travel with the file).
        "_meta": {
            "entry_id": entry_id,
            "n_laps": n_laps,
            "dense_points": int(dense.shape[0]),
            "emitted_points": len(decimated),
            "max_chord_error_km": chord_err_km,
            "max_chord_error_au": chord_err_au,
            "decimation": "douglas-peucker (seam-aware)",
            "dp_tolerance_km_used": float(dmeta["tol_km_used"]),
            "seam_count": seam_count,
            "max_seam_jump_km": float(dmeta["max_seam_jump_km"]),
            "max_seam_jump_au": float(dmeta["max_seam_jump_au"]),
            "cycler_period_days": float(cycler.period) / 86400.0,
            "commit": commit,
        },
    }


def _serialize(payload: dict[str, Any]) -> str:
    return json.dumps(payload, indent=2) + "\n"


def _geometry_view(payload: dict[str, Any]) -> dict[str, Any]:
    """A copy of ``payload`` with the volatile commit-SHA provenance normalised,
    for a deterministic ``--check`` that verifies the GEOMETRY is reproducible (not
    the moving HEAD the file was last emitted at). The commit field legitimately
    tracks when the data was produced, so it must not gate reproducibility —
    re-running on a later commit produces identical geometry with a different SHA,
    which is not "stale" in any meaningful sense (and would otherwise force the
    slow DE440 solve on every check).
    """
    import copy
    import re

    view = copy.deepcopy(payload)
    # Provenance carries "; commit <sha>;" — replace whatever SHA is present.
    view["provenance"] = re.sub(r"commit [0-9a-f]+", "commit <commit>", view["provenance"])
    if "_meta" in view:
        view["_meta"]["commit"] = "<commit>"
    return view


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else "")
    parser.add_argument(
        "--id",
        default="aldrin-classic-em-k1-outbound",
        help="catalogue entry id to export (only the Aldrin outbound is wired)",
    )
    parser.add_argument("--n-laps", type=int, default=_DEFAULT_N_LAPS)
    parser.add_argument("--max-points", type=int, default=_DEFAULT_MAX_POINTS)
    parser.add_argument(
        "--check",
        action="store_true",
        help="verify on-disk file is up to date (exit 1 if stale/missing)",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parent.parent
    out_path = repo_root / "data" / "sampled" / f"{args.id}.json"

    if args.check:
        if not out_path.exists():
            sys.stderr.write(f"export-sampled: {out_path} missing — run without --check.\n")
            raise SystemExit(1)
        payload = build_payload(entry_id=args.id, n_laps=args.n_laps, max_points=args.max_points)
        on_disk = json.loads(out_path.read_text())
        # Compare GEOMETRY (commit SHA normalised): the data is reproducible iff the
        # samples match, independent of which HEAD it was last emitted at.
        if _geometry_view(on_disk) != _geometry_view(payload):
            sys.stderr.write(f"export-sampled: {out_path} is stale — re-run to regenerate.\n")
            raise SystemExit(1)
        sys.stderr.write(f"export-sampled: {out_path} up to date (geometry matches).\n")
        return

    payload = build_payload(entry_id=args.id, n_laps=args.n_laps, max_points=args.max_points)
    text = _serialize(payload)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text)
    meta = payload["_meta"]
    sys.stderr.write(
        f"export-sampled: wrote {out_path} ({len(text)} bytes); "
        f"{meta['emitted_points']} pts from {meta['dense_points']} dense, "
        f"max chord error {meta['max_chord_error_au']:.2e} AU "
        f"({meta['max_chord_error_km']:.0f} km).\n"
    )


if __name__ == "__main__":
    main()
