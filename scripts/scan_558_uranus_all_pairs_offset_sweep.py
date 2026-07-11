"""#558 -- #312 family census via the original discovery genome.

Generalizes ``scripts/scan_312_uranus_umbriel_oberon_offset_sweep.py``'s
``sweep_offset()`` (relative-phase-offset sweep) from the single Umbriel-
Oberon pair to **all 10 Uranian regular-moon pairs**
(Miranda/Ariel/Umbriel/Titania/Oberon, C(5,2)=10), at BOTH anchor directions
per pair (A-B-A and B-A-B -- these are physically distinct: the anchor moon
is the one the tour repeatedly returns to).

Background (see ``data/OUTSTANDING.md`` #558, and #540 immediately above it):
the discovery genome that found #312 (the ONE confirmed novel catalogue row
in this project's history) seeds each moon's initial longitude via a FIXED
convention (``2*pi*j/len(moons)``, ``j`` = sorted-index) and never varies the
RELATIVE offset between the two moons in a pair -- only the GLOBAL absolute
phase. The addendum sweep (``scan_312_uranus_umbriel_oberon_offset_sweep.py``)
first swept this relative offset explicitly for Umbriel-Oberon and found the
fixed-convention residual (0.636 km/s, 3-moon convention) was 25x worse than
the true basin floor (0.025 km/s, the 2-moon/180-degree convention actually
catalogued). Every OTHER Uranian pair's "miss" on record was measured ONLY at
the fixed-convention grid -- this script asks whether any of them ALSO hide a
sub-gate basin once the relative offset is swept properly.

Efficiency note (n_rev dimension): a naive n_rev in [0..3] sweep would need a
fresh ``lambert()`` call per (rel_offset, phase0, tof_scale, n_rev) tuple --
16x more Lambert solves than the (n_rev fixed) original script. Instead, each
leg is solved ONCE per grid point with ``max_revs=3`` (which returns ALL
n_revs in [0,3] x both branches in a single call -- the Stumpff root-finding
this shares across revolution counts is the expensive part, not the
bookkeeping), and every (n_rev_leg0, n_rev_leg1) combination is then scored
from the already-computed solution set at ~zero marginal cost. This keeps
the n_rev axis from blowing up runtime while still covering the full 0..3
range per leg per the #558 spec.

Efficiency note (phase0 dimension is PROVABLY REDUNDANT): the circular-
coplanar Kepler + patched-conic-Lambert closure this genome uses is exactly
rotationally symmetric about the primary's polar axis. Adding a constant
``phase0`` to BOTH moons' initial longitudes rigidly rotates the entire
3-state configuration (``r0,v0,r1,v1,r2,v2``) by the same angle; the Lambert
problem is rotation-covariant (the solved ``v1,v2`` rotate the same way);
and the residual is built ONLY from norms of vector differences, which
rotations preserve exactly. So, for FIXED ``rel_offset`` and ``tof_scale``,
the residual (and every V_inf magnitude) is mathematically independent of
``phase0`` -- verified empirically here to ~1e-13 relative agreement across
6+ phase0 samples spanning 0-359 deg, on 2 different pairs (Umbriel-Oberon
AND Miranda-Titania/Ariel-Oberon), not a coincidence specific to one pair.
This means the original script's 96-sample ``n_phase`` grid was pure wasted
compute -- a 96x tax for zero information. This script exploits that: the
default ``n_phase=1`` (a single representative sample), and the freed-up
budget is spent on much denser ``rel_offset``/``tof_scale`` grids instead
(see ``main()`` CLI defaults). ``n_phase`` remains a real, settable
parameter (never silently dropped) so this claim stays auditable/falsifiable.

Positive control (run FIRST, must pass before trusting anything else): at
Umbriel-Oberon (anchor=Umbriel), rel_offset near 180 deg, tof_scale=2.0,
n_rev=(1,1), this script must reproduce residual ~0.025 km/s -- NOT the
0.636 km/s fixed-convention artifact the original discovery run recorded.

Discipline: NO catalogue writeback. NO novelty claims here -- this is the
Sonnet-tier mechanical sweep; Opus + Fable adjudicate any survivor.

Run as::

    uv run python scripts/scan_558_uranus_all_pairs_offset_sweep.py
"""

from __future__ import annotations

import itertools
import json
import math
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from cyclerfinder.core.lambert import (  # noqa: E402
    LambertConvergenceError,
    LambertGeometryError,
    lambert,
)
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES  # noqa: E402
from cyclerfinder.search.discovery_campaign import (  # noqa: E402
    DAY_S,
    _mean_motion_rad_day,
    _moon_state,
)
from cyclerfinder.search.five_tier_prioritizer import PatchedConicLeg  # noqa: E402
from cyclerfinder.search.physical_sanity import (  # noqa: E402
    DEFAULT_MIN_USEFUL_BEND_DEG,
    candidate_passes_physical_gate,
)
from cyclerfinder.search.saturn_uranus_campaign import dop853_cross_check_leg  # noqa: E402

MOONS: tuple[str, ...] = ("Miranda", "Ariel", "Umbriel", "Titania", "Oberon")
GATE_RESIDUAL_KMS = 0.05
N_REV_MAX = 3  # per-leg n_rev in [0, N_REV_MAX], per #558 spec ("0..3")

DATA_DIR = ROOT / "data"


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT
        ).strip()
    except Exception:
        return "unknown"


@dataclass(frozen=True)
class LegOption:
    n_revs: int
    vinf_out: float  # |v1_lambert - v_departure_body| at the leg's start
    vinf_in: float  # |v2_lambert - v_arrival_body| at the leg's end
    v1: np.ndarray
    v2: np.ndarray


def _leg_options(
    r_a: np.ndarray,
    v_a: np.ndarray,
    r_b: np.ndarray,
    v_b: np.ndarray,
    tof_s: float,
    mu: float,
    n_rev_max: int,
) -> dict[int, LegOption]:
    """All achievable ``n_revs in [0, n_rev_max]`` solutions for one leg.

    One ``lambert()`` call with ``max_revs=n_rev_max`` returns every feasible
    revolution count in a single solve (the root-finding is shared); this
    lets the caller score every ``(n_rev_leg0, n_rev_leg1)`` combination
    without re-solving Lambert per combination (the efficiency note above).
    """
    try:
        sols = lambert(r_a, r_b, tof_s, mu=mu, max_revs=n_rev_max)
    except (LambertGeometryError, LambertConvergenceError):
        return {}
    by_nrev: dict[int, list[Any]] = {}
    for s in sols:
        by_nrev.setdefault(s.n_revs, []).append(s)
    out: dict[int, LegOption] = {}
    for nrev, cands in by_nrev.items():
        best = min(cands, key=lambda s: float(np.linalg.norm(s.v1 - v_a)))
        out[nrev] = LegOption(
            n_revs=nrev,
            vinf_out=float(np.linalg.norm(best.v1 - v_a)),
            vinf_in=float(np.linalg.norm(best.v2 - v_b)),
            v1=best.v1,
            v2=best.v2,
        )
    return out


def residual_at_point(
    anchor: str,
    flyby: str,
    *,
    rel_offset_deg: float,
    tof_scale: float,
    n_rev: tuple[int, int],
    phase0_deg: float = 0.0,
    primary: str = "Uranus",
) -> dict[str, Any] | None:
    """Residual + V_inf at ONE explicit (rel_offset, tof_scale, n_rev) point.

    Used for the #312 reproduction check: this evaluates the EXACT known
    catalogued point directly (rel_offset=180 deg, tof_scale=2.0, n_rev=(1,1))
    rather than relying on a grid sample landing near it, so the positive
    control is a faithful reproduction check, not a coincidental nearby hit.

    ``primary`` (#571 genericization): the central body whose ``PRIMARIES``
    GM and whose moons' ``SATELLITES`` entries govern the two-body Kepler
    states below. Defaults to ``"Uranus"`` so every #558-era caller (which
    never passes this argument) is byte-for-byte unaffected -- the residual
    formula and gate logic below are untouched, only the body lookup is
    parameterized.
    """
    mu = PRIMARIES[primary]
    sat_a = SATELLITES[anchor]
    sat_b = SATELLITES[flyby]
    n_a = _mean_motion_rad_day(mu, sat_a.sma_km)
    n_b = _mean_motion_rad_day(mu, sat_b.sma_km)
    p_a = 2.0 * math.pi / n_a
    p_b = 2.0 * math.pi / n_b

    phase0 = math.radians(phase0_deg)
    rel_off = math.radians(rel_offset_deg)
    tof = tof_scale * math.sqrt(p_a * p_b)
    tof_s = tof * DAY_S
    theta_a = phase0
    theta_b = phase0 + rel_off
    r0, v0 = _moon_state(theta_a, n_a, 0.0, sat_a.sma_km, mu)
    r1, v1 = _moon_state(theta_b, n_b, tof, sat_b.sma_km, mu)
    r2, v2 = _moon_state(theta_a, n_a, 2.0 * tof, sat_a.sma_km, mu)

    n0, n1 = n_rev
    leg0 = _leg_options(r0, v0, r1, v1, tof_s, mu, max(n0, 1))
    leg1 = _leg_options(r1, v1, r2, v2, tof_s, mu, max(n1, 1))
    if n0 not in leg0 or n1 not in leg1:
        return None
    opt0, opt1 = leg0[n0], leg1[n1]
    r_mid = abs(opt0.vinf_in - opt1.vinf_out)
    r_periodic = abs(opt0.vinf_out - opt1.vinf_in)
    worst = max(r_mid, r_periodic)
    return {
        "anchor": anchor,
        "flyby": flyby,
        "rel_offset_deg": rel_offset_deg,
        "phase0_deg": phase0_deg,
        "tof_scale": tof_scale,
        "tof_days": tof,
        "n_rev": [n0, n1],
        "residual_kms": worst,
        "vinf_in": [0.0, opt0.vinf_in, opt1.vinf_in],
        "vinf_out": [opt0.vinf_out, opt1.vinf_out, 0.0],
    }


def sweep_pair(
    anchor: str,
    flyby: str,
    *,
    n_phase: int = 96,
    n_offset: int = 96,
    tof_scales: tuple[float, ...] = (0.5, 1.0, 1.5, 2.0),
    n_rev_max: int = N_REV_MAX,
    keep_top: int = 200,
    max_gate_passing_records: int = 20_000,
    primary: str = "Uranus",
) -> dict[str, Any]:
    """Sweep relative phase offset (+ global phase, tof_scale, n_rev) for the
    closed length-3 cycle ``anchor - flyby - anchor``.

    Memory/runtime note: at each ``(rel_offset, phase0, tof_scale)`` grid
    point, EVERY ``n_rev in [0, n_rev_max]`` combination for both legs is
    evaluated (from the SAME pair of ``lambert()`` calls -- see
    ``_leg_options``), but only the single BEST (lowest-residual) n_rev
    combination at that point is retained as the point's record. This keeps
    the record count equal to the (rel_offset, phase0, tof_scale) grid size
    (not 16x larger from the n_rev cross product) while still guaranteeing
    the n_rev axis is searched, not silently dropped, at every point. A
    streaming bounded max-heap (size ``keep_top``) tracks the best records
    without ever materializing the full grid in memory, EXCEPT when
    ``keep_top`` is set large enough to exceed the grid size (the explicit
    full-landscape-dump mode for Umbriel-Oberon, #558 spec item 5).

    Returns a summary dict with the best overall record, the best
    gate-passing (< 0.05 km/s) records (up to ``max_gate_passing_records``),
    and the top ``keep_top`` records by residual (for multi-basin census --
    NOT just the single best point per pair).
    """
    import heapq

    mu = PRIMARIES[primary]
    sat_a = SATELLITES[anchor]
    sat_b = SATELLITES[flyby]
    sma_a, sma_b = sat_a.sma_km, sat_b.sma_km
    n_a = _mean_motion_rad_day(mu, sma_a)
    n_b = _mean_motion_rad_day(mu, sma_b)
    p_a = 2.0 * math.pi / n_a
    p_b = 2.0 * math.pi / n_b

    best_overall = math.inf
    best_rec: dict[str, Any] = {}
    n_grid_points = 0
    n_feasible = 0
    n_gate_passing = 0
    gate_passing_records: list[dict[str, Any]] = []

    # Bounded max-heap of the ``keep_top`` smallest-residual records seen so
    # far: entries are ``(-residual, tie_break_counter, rec)`` so the heap
    # root (smallest tuple) is always the WORST of the currently-kept set.
    heap: list[tuple[float, int, dict[str, Any]]] = []
    counter = 0

    for i in range(n_offset):
        rel_off = 2.0 * math.pi * i / n_offset
        for j in range(n_phase):
            phase0 = 2.0 * math.pi * j / n_phase
            theta_a = phase0
            theta_b = phase0 + rel_off
            for ts in tof_scales:
                n_grid_points += 1
                tof = ts * math.sqrt(p_a * p_b)
                tof_s = tof * DAY_S
                r0, v0 = _moon_state(theta_a, n_a, 0.0, sma_a, mu)
                r1, v1 = _moon_state(theta_b, n_b, tof, sma_b, mu)
                r2, v2 = _moon_state(theta_a, n_a, 2.0 * tof, sma_a, mu)

                leg0 = _leg_options(r0, v0, r1, v1, tof_s, mu, n_rev_max)
                leg1 = _leg_options(r1, v1, r2, v2, tof_s, mu, n_rev_max)
                if not leg0 or not leg1:
                    continue
                n_feasible += 1

                # Best n_rev combination at THIS (rel_offset, phase0,
                # tof_scale) point -- see the memory/runtime note above.
                point_best_worst = math.inf
                point_best_rec: dict[str, Any] | None = None
                for n0, opt0 in leg0.items():
                    for n1, opt1 in leg1.items():
                        # Continuity at the middle flyby (leg0 arrival vs
                        # leg1 departure) + periodicity (leg0 departure vs
                        # leg1 arrival, i.e. the anchor's V_inf must repeat
                        # cycle-to-cycle) -- same residual definition as
                        # scan_312_uranus_umbriel_oberon_offset_sweep.py.
                        r_mid = abs(opt0.vinf_in - opt1.vinf_out)
                        r_periodic = abs(opt0.vinf_out - opt1.vinf_in)
                        worst = max(r_mid, r_periodic)
                        if worst < point_best_worst:
                            point_best_worst = worst
                            point_best_rec = {
                                "anchor": anchor,
                                "flyby": flyby,
                                "rel_offset_deg": math.degrees(rel_off),
                                "phase0_deg": math.degrees(phase0),
                                "tof_scale": ts,
                                "tof_days": tof,
                                "n_rev": [n0, n1],
                                "residual_kms": worst,
                                "vinf_in": [0.0, opt0.vinf_in, opt1.vinf_in],
                                "vinf_out": [opt0.vinf_out, opt1.vinf_out, 0.0],
                            }
                if point_best_rec is None:
                    continue

                if point_best_worst < best_overall:
                    best_overall = point_best_worst
                    best_rec = dict(point_best_rec)

                if point_best_worst < GATE_RESIDUAL_KMS:
                    n_gate_passing += 1
                    if len(gate_passing_records) < max_gate_passing_records:
                        gate_passing_records.append(point_best_rec)

                counter += 1
                key = -point_best_worst
                if len(heap) < keep_top:
                    heapq.heappush(heap, (key, counter, point_best_rec))
                elif key > heap[0][0]:
                    heapq.heapreplace(heap, (key, counter, point_best_rec))

    top_records = [rec for _key, _cnt, rec in heap]
    top_records.sort(key=lambda r: r["residual_kms"])

    return {
        "anchor": anchor,
        "flyby": flyby,
        "sequence": [anchor, flyby, anchor],
        "n_phase": n_phase,
        "n_offset": n_offset,
        "tof_scales": list(tof_scales),
        "n_rev_max": n_rev_max,
        "n_grid_points": n_grid_points,
        "n_feasible_grid_points": n_feasible,
        "n_records_total": n_feasible,
        "n_gate_passing_records": n_gate_passing,
        "best_overall_residual_kms": (best_overall if math.isfinite(best_overall) else None),
        "best_overall_record": best_rec,
        "top_records": top_records,
        "gate_passing_records": gate_passing_records,
    }


def build_legs_for_record(
    anchor: str, flyby: str, rec: dict[str, Any], *, primary: str = "Uranus"
) -> list[PatchedConicLeg]:
    """Reconstruct SI-units PatchedConicLeg objects for a specific record.

    Used ONLY for candidates that already passed the residual + physical
    gates, to run the independent DOP853 cross-check -- mirrors
    ``verify_327_umbriel_silver.py``'s pattern but built directly from the
    sweep's own (rel_offset, phase0, tof_scale, n_rev) rather than through
    the production ``close()`` convention (which cannot express an arbitrary
    relative offset).

    ``primary`` (#571 genericization, default ``"Uranus"`` for exact
    backward compatibility with every #558/#562/#563 caller).
    """
    mu = PRIMARIES[primary]
    sat_a = SATELLITES[anchor]
    sat_b = SATELLITES[flyby]
    sma_a, sma_b = sat_a.sma_km, sat_b.sma_km
    n_a = _mean_motion_rad_day(mu, sma_a)
    n_b = _mean_motion_rad_day(mu, sma_b)

    rel_off = math.radians(rec["rel_offset_deg"])
    phase0 = math.radians(rec["phase0_deg"])
    tof = rec["tof_days"]
    n0, n1 = rec["n_rev"]

    theta_a = phase0
    theta_b = phase0 + rel_off
    r0, v0 = _moon_state(theta_a, n_a, 0.0, sma_a, mu)
    r1, v1 = _moon_state(theta_b, n_b, tof, sma_b, mu)
    r2, _v2 = _moon_state(theta_a, n_a, 2.0 * tof, sma_a, mu)

    tof_s = tof * DAY_S
    sols0 = lambert(r0, r1, tof_s, mu=mu, max_revs=max(0, n0))
    best0 = min(
        (s for s in sols0 if s.n_revs == n0), key=lambda s: float(np.linalg.norm(s.v1 - v0))
    )
    sols1 = lambert(r1, r2, tof_s, mu=mu, max_revs=max(0, n1))
    best1 = min(
        (s for s in sols1 if s.n_revs == n1), key=lambda s: float(np.linalg.norm(s.v1 - v1))
    )

    km_m = 1000.0
    mu_m3_s2 = mu * (km_m**3)
    return [
        PatchedConicLeg(
            label_from=anchor,
            label_to=flyby,
            r1_m=r0 * km_m,
            v1_m_s=best0.v1 * km_m,
            r2_m=r1 * km_m,
            v2_m_s=best0.v2 * km_m,
            dt_s=tof_s,
            mu_m3_s2=mu_m3_s2,
        ),
        PatchedConicLeg(
            label_from=flyby,
            label_to=anchor,
            r1_m=r1 * km_m,
            v1_m_s=best1.v1 * km_m,
            r2_m=r2 * km_m,
            v2_m_s=best1.v2 * km_m,
            dt_s=tof_s,
            mu_m3_s2=mu_m3_s2,
        ),
    ]


def encounter_vinfs_kms(rec: dict[str, Any]) -> tuple[float, float, float]:
    """Per-encounter V_inf magnitude, matching rerun_324_physical_gate.py's
    offset-sweep-style extraction (max of the in/out asymptote at each node).
    """
    vin = rec["vinf_in"]
    vout = rec["vinf_out"]
    return tuple(max(abs(vin[k]), abs(vout[k])) for k in range(3))  # type: ignore[return-value]


def gate_candidate(
    anchor: str, flyby: str, rec: dict[str, Any], *, primary: str = "Uranus"
) -> dict[str, Any]:
    """Run the #324 physical-bend gate + DOP853 cross-check on one record.

    ``primary`` (#571 genericization, default ``"Uranus"``): forwarded only
    to :func:`build_legs_for_record` for the Kepler-state reconstruction --
    the gate logic itself (:func:`candidate_passes_physical_gate`) is
    body-agnostic (it resolves each body's own GM/radius/safe_alt from the
    ``PLANETS``/``SATELLITES`` registries via the body NAME in ``seq``, not
    via ``primary``), so this is a pure pass-through, not a gate change.
    """
    seq = (anchor, flyby, anchor)
    vinfs = encounter_vinfs_kms(rec)
    physical_pass, verdicts = candidate_passes_physical_gate(
        seq, vinfs, min_useful_bend_deg=DEFAULT_MIN_USEFUL_BEND_DEG
    )
    legs = build_legs_for_record(anchor, flyby, rec, primary=primary)
    cross_checks = [dop853_cross_check_leg(leg, rtol=1e-12, atol=1e-12) for leg in legs]
    max_dr_km = max(float(cc["dr_arrival_km"]) for cc in cross_checks)
    independent_pass = max_dr_km < 1.0  # < 1 km, matches verify_327's threshold

    return {
        "anchor": anchor,
        "flyby": flyby,
        "record": rec,
        "vinf_per_encounter_kms": list(vinfs),
        "physical_gate_passed": physical_pass,
        "max_bend_deg_per_encounter": [v.max_bend_deg for v in verdicts],
        "dop853_cross_check": {
            "max_dr_arrival_km": max_dr_km,
            "per_leg": cross_checks,
            "passed": independent_pass,
        },
        "all_gates_passed": bool(physical_pass and independent_pass),
    }


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--n-phase",
        type=int,
        default=1,
        help="global-phase samples. Default 1: phase0 is PROVABLY redundant "
        "given fixed rel_offset/tof_scale (see module docstring) -- the "
        "freed-up budget goes to --n-offset/--tof-scales density instead.",
    )
    parser.add_argument(
        "--n-offset",
        type=int,
        default=360,
        help="relative-offset samples over the full 360 deg (default: 1 deg resolution)",
    )
    parser.add_argument(
        "--tof-scales",
        type=str,
        default="dense",
        help="'dense' (0.5-3.0 step 0.05, 51 values) or 'coarse-dense' (step "
        "0.1, 26 values, the literal #558 spec density) or 'coarse' "
        "(0.5,1.0,1.5,2.0, the original script's grid) or a comma-separated list",
    )
    parser.add_argument("--n-rev-max", type=int, default=N_REV_MAX)
    parser.add_argument("--keep-top", type=int, default=200)
    parser.add_argument(
        "--positive-control-only",
        action="store_true",
        help="Only run the Umbriel-Oberon positive control and exit.",
    )
    parser.add_argument(
        "--pairs",
        type=str,
        default="all",
        help="'all' or comma-separated 'Anchor:Flyby' pairs (single direction each)",
    )
    parser.add_argument(
        "--full-dump-pair",
        type=str,
        default="",
        help="'Anchor:Flyby' -- dump the FULL grid (not just top-N) for this pair+direction",
    )
    args = parser.parse_args()

    if args.tof_scales == "dense":
        tof_scales = tuple(round(0.5 + 0.05 * k, 2) for k in range(51))  # 0.5..3.0 step 0.05
    elif args.tof_scales == "coarse-dense":
        tof_scales = tuple(round(0.5 + 0.1 * k, 2) for k in range(26))  # 0.5..3.0 step 0.1
    elif args.tof_scales == "coarse":
        tof_scales = (0.5, 1.0, 1.5, 2.0)
    else:
        tof_scales = tuple(float(x) for x in args.tof_scales.split(","))

    sha = _git_sha()
    print(f"[558] Uranus all-pair relative-offset sweep -- sha={sha}", flush=True)
    print(
        f"[558] n_phase={args.n_phase} n_offset={args.n_offset} "
        f"tof_scales({len(tof_scales)})={tof_scales[:3]}...{tof_scales[-1]} "
        f"n_rev_max={args.n_rev_max}",
        flush=True,
    )

    # --- Positive control -----------------------------------------------
    # First, the FAITHFUL reproduction check: evaluate the EXACT catalogued
    # point directly (not a grid sample landing nearby).
    print(
        "[558] POSITIVE CONTROL (exact #312 point): rel_off=180, tof_scale=2.0, n_rev=(1,1)...",
        flush=True,
    )
    ref = residual_at_point("Umbriel", "Oberon", rel_offset_deg=180.0, tof_scale=2.0, n_rev=(1, 1))
    ref_ok = ref is not None and ref["residual_kms"] < 0.03
    if ref is not None:
        print(
            f"[558]   exact-point residual = {ref['residual_kms']:.6f} km/s "
            f"(reference catalogued value 0.025232 km/s) -- reproduces: {ref_ok}",
            flush=True,
        )
    else:
        print("[558]   exact-point evaluation INFEASIBLE (no Lambert solution)", flush=True)

    print("[558] POSITIVE CONTROL: Umbriel-Oberon-Umbriel (anchor=Umbriel)...", flush=True)
    t0 = time.time()
    pc = sweep_pair(
        "Umbriel",
        "Oberon",
        n_phase=args.n_phase,
        n_offset=args.n_offset,
        tof_scales=tof_scales,
        n_rev_max=args.n_rev_max,
        keep_top=args.keep_top,
    )
    pc_elapsed = time.time() - t0
    print(
        f"[558] positive control: best residual = {pc['best_overall_residual_kms']:.6f} km/s "
        f"at rel_off={pc['best_overall_record']['rel_offset_deg']:.2f} deg, "
        f"tof_scale={pc['best_overall_record']['tof_scale']:.2f}, "
        f"n_rev={pc['best_overall_record']['n_rev']}  (elapsed {pc_elapsed:.1f}s)",
        flush=True,
    )
    pc_ok = bool(
        ref_ok
        and pc["best_overall_residual_kms"] is not None
        and pc["best_overall_residual_kms"] < 0.03
    )
    print(
        f"[558] positive control reproduces ~0.025 km/s basin floor AT THE CATALOGUED "
        f"POINT, and the grid search finds no worse overall basin: {pc_ok}",
        flush=True,
    )

    pc_out = DATA_DIR / "scan_558_positive_control.jsonl"
    with pc_out.open("w", encoding="utf-8") as fh:
        fh.write(
            json.dumps(
                {
                    "_meta": True,
                    "task": "#558 positive control -- Umbriel-Oberon reproduction",
                    "git_sha": sha,
                    "elapsed_s": pc_elapsed,
                    "positive_control_pass": pc_ok,
                    "reference_2moon_residual_kms": 0.025232272564609692,
                    "reference_3moon_residual_kms": 0.635986,
                }
            )
            + "\n"
        )
        if ref is not None:
            fh.write(json.dumps({"kind": "exact_reference_point", **ref}) + "\n")
        fh.write(json.dumps({"kind": "best_overall", **pc["best_overall_record"]}) + "\n")
        for rec in pc["top_records"][:20]:
            fh.write(json.dumps({"kind": "top20", **rec}) + "\n")
    print(f"[558] positive control written to {pc_out}", flush=True)

    if args.positive_control_only:
        return 0 if pc_ok else 1
    if not pc_ok:
        print("[558] ABORT: positive control FAILED -- not proceeding to full sweep.", flush=True)
        return 1

    # --- Full 10-pair x 2-direction sweep --------------------------------
    if args.pairs == "all":
        pair_dirs: list[tuple[str, str]] = []
        for a, b in itertools.combinations(MOONS, 2):
            pair_dirs.append((a, b))
            pair_dirs.append((b, a))
    else:
        pair_dirs = []
        for tok in args.pairs.split(","):
            a, b = tok.split(":")
            pair_dirs.append((a, b))

    full_dump_target: tuple[str, str] | None = None
    if args.full_dump_pair:
        a, b = args.full_dump_pair.split(":")
        full_dump_target = (a, b)

    index_rows: list[dict[str, Any]] = []
    candidates_needing_adjudication: list[dict[str, Any]] = []

    for anchor, flyby in pair_dirs:
        slug = f"{anchor.lower()}_{flyby.lower()}"
        out_path = DATA_DIR / f"scan_558_uranus_{slug}.jsonl"
        t0 = time.time()
        dump_full = full_dump_target == (anchor, flyby)
        keep_top = 10_000_000 if dump_full else args.keep_top
        result = sweep_pair(
            anchor,
            flyby,
            n_phase=args.n_phase,
            n_offset=args.n_offset,
            tof_scales=tof_scales,
            n_rev_max=args.n_rev_max,
            keep_top=keep_top,
        )
        elapsed = time.time() - t0
        best = result["best_overall_residual_kms"]
        print(
            f"[558] {anchor}-{flyby}-{anchor}: best={best} n_records={result['n_records_total']} "
            f"n_gate_passing={result['n_gate_passing_records']}  ({elapsed:.1f}s)",
            flush=True,
        )

        # Gate every sub-0.05 record through physical + DOP853.
        gated: list[dict[str, Any]] = []
        for rec in result["gate_passing_records"]:
            g = gate_candidate(anchor, flyby, rec)
            gated.append(g)
            if g["all_gates_passed"]:
                candidates_needing_adjudication.append(g)

        with out_path.open("w", encoding="utf-8") as fh:
            fh.write(
                json.dumps(
                    {
                        "_meta": True,
                        "task": "#558 Uranus all-pair offset sweep",
                        "anchor": anchor,
                        "flyby": flyby,
                        "sequence": [anchor, flyby, anchor],
                        "n_phase": result["n_phase"],
                        "n_offset": result["n_offset"],
                        "tof_scales": result["tof_scales"],
                        "n_rev_max": result["n_rev_max"],
                        "elapsed_s": elapsed,
                        "git_sha": sha,
                    }
                )
                + "\n"
            )
            fh.write(json.dumps({"kind": "best_overall", **result["best_overall_record"]}) + "\n")
            # Full-landscape dump for the requested pair (e.g. Umbriel-Oberon
            # itself, per #558 spec item 5): every record, not just top-N
            # (keep_top was set to a huge number above, so top_records IS
            # the full grid here).
            kind = "full_grid" if dump_full else "top_record"
            for rec in result["top_records"]:
                fh.write(json.dumps({"kind": kind, **rec}) + "\n")
            for g in gated:
                fh.write(json.dumps({"kind": "gate_result", **g}) + "\n")

        index_rows.append(
            {
                "anchor": anchor,
                "flyby": flyby,
                "path": str(out_path.relative_to(ROOT)),
                "best_overall_residual_kms": best,
                "n_gate_passing_records": result["n_gate_passing_records"],
                "n_all_gates_passed": sum(1 for g in gated if g["all_gates_passed"]),
                "elapsed_s": elapsed,
            }
        )

    index_out = DATA_DIR / "scan_558_uranus_all_pairs_index.jsonl"
    with index_out.open("w", encoding="utf-8") as fh:
        fh.write(
            json.dumps(
                {
                    "_meta": True,
                    "task": "#558 Uranus all-pair offset sweep -- index",
                    "git_sha": sha,
                    "n_pair_directions": len(pair_dirs),
                    "positive_control_pass": pc_ok,
                }
            )
            + "\n"
        )
        for row in index_rows:
            fh.write(json.dumps(row) + "\n")
        fh.write(
            json.dumps(
                {
                    "_meta": True,
                    "kind": "candidates_needing_adjudication",
                    "count": len(candidates_needing_adjudication),
                    "candidates": candidates_needing_adjudication,
                }
            )
            + "\n"
        )
    print(f"[558] DONE -- index: {index_out}", flush=True)
    print(
        f"[558] candidates passing ALL gates (residual+physical+DOP853): "
        f"{len(candidates_needing_adjudication)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
