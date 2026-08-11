"""#125 Part 2 — descriptor-seeded Russell-12 closure campaign.

For each of the 12 ``free_return_arcs[]``-bearing rows (search/descriptor.py's
12-row gate), this driver:

  1. parses the row's descriptor + trajectory segments into the corrector
     genome (sequence, per-leg ``(n_revs, branch)``, per-leg ToF seed);
  2. drives :func:`cyclerfinder.search.correct.ballistic_correct` over a modest
     epoch grid (16-core) in BOTH residual modes (magnitude AND vector);
  3. compares every closed solution against the row's SOURCED anchors
     (V∞ multiset, invariants turn_ratio / aphelion_ratio / transit_times,
     period);
  4. classifies the per-row outcome:
       - CLOSE-AND-MATCH  (closes within documented tolerances vs sourced anchors)
       - CLOSE-OFF-ANCHOR (closes but does not match a sourced anchor)
       - NO-CLOSE         (corrector did not reach the residual floor)

GOLDEN DISCIPLINE (project memory feedback_golden_tests_sourced_only): the
EXPECTED side of every match check is the row's SOURCED catalogue anchor
(``vinf_kms_at_encounters`` from Russell 2004 / McConaghy, ``invariants`` from
Russell tables, ``period.years``). Nothing our own code computes is ever used as
an EXPECTED. Tolerances are the proposed campaign defaults (V∞ 0.5 km/s, AR/TR
0.05, transits 5 d), documented and applied uniformly — never loosened per-row.

Topology mapping (the descriptor->corrector genome)
---------------------------------------------------
All 12 rows are multi-arc chains of Earth-Earth legs (Russell 2004 §4.8). The
DESIGNATED leg — the UPPERCASE letter in the row's own descriptor string, NOT
always ``arcs[0]`` (#794/#820: e.g. ``arcs[1]`` (G) for 5.30gGf3, ``arcs[2]``
(F) for the ggF/gfF rows) — is the Mars-transit leg; the remaining (lowercase)
arcs, in cyclic itinerary order after it, are the Earth-Earth phasing loops.

The genome splits the designated leg at the Mars encounter: leg 0 is E->M with
the SOURCED taxi transit t_out (segment ``out-em``); leg 1 is M->E with the
REMAINDER of the designated leg (its printed full ToF minus t_out). The
pre-#820 posing used the sourced t_in as leg 1 instead — but t_in is the
inbound-crossing taxi transit (ridden on another cycle's leg), not a per-cycle
leg, so the beyond-Mars portion of the designated conic was unaccounted and
the seeds could never tile the period. Leg 1's ``(n_revs, branch)`` is not
printed; the small candidate set implied by the printed transfer angle (g/G)
or resonance (f/F) is discriminated by residual-at-truth
(:func:`select_leg1_topology`). Each E->E loop leg seeds ``(ToF, n_revs,
branch)`` directly from the row's #794-written-back ``loop-ee-*`` segments
(primary-source-derived; M:N per the #794 convention fix — this script
previously duplicated the parse REVERSED, see the #813 blast-radius note).
The longest E->E loop is eliminated as the period slack leg (spec §2.1(a)),
exactly as ``test_correct_s1l1.py`` pins it.

Like-for-like variant (#135)
-----------------------------
The #125 Part-2 run drove the corrector on ``Ephemeris('astropy')`` (real
DE440) and compared closures to the rows' CIRCULAR-COPLANAR sourced anchors —
cross-fidelity confounded the basin question with the model mismatch. Pass
``--model circular`` to re-run the SAME genome derivation, tolerances and
anchors on ``Ephemeris('circular')``: now coplanar-vs-coplanar, against rows
that are by-construction solutions OF the circular-coplanar model.

The ``--probe-at-truth`` flag adds the decisive diagnostic: per row, seed the
corrector EXACTLY at the row's own sourced ToF geometry (transit ToFs from the
segments, E-E loop ToFs from the descriptor arcs) at the best-phase epoch
(t0 chosen as the phase minimising the residual evaluated AT the sourced
geometry), then check whether the corrector STAYS at truth (residual->0,
ToFs unchanged) or WALKS AWAY to the degenerate basin. Staying => seeding is
the whole story; walking => the residual/solver itself is implicated.

Usage::

    uv run python scripts/campaign_russell12.py [--epochs N] [--workers W] \\
        [--model circular|astropy] [--probe-at-truth] [--out FILE]
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
from concurrent.futures import ProcessPoolExecutor
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
import yaml  # type: ignore[import-untyped]

from cyclerfinder.core.constants import DAYS_PER_JULIAN_YEAR
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.data.runlog import RunLog, RunRecord, default_runlog_path
from cyclerfinder.search.correct import _residuals, ballistic_correct
from cyclerfinder.search.descriptor import arc_tof_seed_days
from cyclerfinder.search.free_return import _residuals as _fr_residuals
from cyclerfinder.search.free_return import (
    free_return_correct,
    seed_ae_from_aphelion_transit,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
CATALOGUE_PATH = REPO_ROOT / "data" / "catalogue.yaml"
DEFAULT_RUNLOG_DIR = REPO_ROOT / "data" / "runs"
DAY_S = 86400.0
J2000 = datetime(2000, 1, 1, 12, 0, 0, tzinfo=UTC)

# Proposed campaign match tolerances vs the row's SOURCED anchors.
TOL_VINF_KMS = 0.5
TOL_RATIO = 0.05  # aphelion_ratio, turn_ratio (dimensionless)
TOL_TRANSIT_DAYS = 5.0
TOL_PERIOD_YEARS = 0.05
CORRECTOR_TOL_KMS = 0.1  # convergence floor (prototype threshold)
VINF_CAP_KMS = 12.0  # generous; these rows run 4.6-10.8 km/s at Mars
# Dense phase-scan floor for the (cheap, Lambert-free) free-return t0 search.
# 256 misses the narrow basin of deep-aphelion high-e rows (#137 Part 3).
FR_PHASE_EPOCHS_FLOOR = 4096
# t0 phase points for the leg-1 (n_revs, branch) residual-at-truth selection
# in run_row (#820); the probe path reuses its own --phase-epochs instead.
LEG1_SELECT_PHASE_EPOCHS = 256

RUSSELL12_IDS = (
    "mcconaghy-2006-em-k2",
    "russell-ch4-4.991gG2",
    "russell-ch4-8.049gGf2",
    "russell-ch4-9.353Gg2",
    "russell-ch4-3.64gGg3",
    "russell-ch4-3.78Gg3",
    "russell-ch4-5.30gGf3",
    "russell-ch4-9.94Gg3",
    "russell-ch4-3.66gfF3",
    "russell-ch4-5.30ggF3",
    "russell-ch4-5.75ggF3",
    "russell-ch4-6.44Gg3",
)


def _t_sec(dt: datetime) -> float:
    return (dt - J2000).total_seconds()


def _t0_center(row: dict[str, Any]) -> float:
    """Epoch-scan centre: the row's priority date, else 2030-01-01."""
    priority = row.get("priority_date")
    if priority:
        return _t_sec(datetime.fromisoformat(str(priority)).replace(tzinfo=UTC))
    return _t_sec(datetime(2030, 1, 1, tzinfo=UTC))


# Non-Lambert branch labels used by #794's written-back loop segments, mapped
# onto the corrector's Lambert vocabulary (single/low/high). A "resonant"
# full-rev loop (N revs in M years) is posed on the multi-rev low branch,
# matching search/descriptor.py::arc_to_leg_topology.
_SEGMENT_BRANCH_MAP = {"resonant": "low", "half-rev": "single"}


def _designated_arc_index(arcs: list[dict[str, Any]]) -> int:
    """Index of the DESIGNATED arc: the UPPERCASE letter in Russell's own
    leg-descriptor notation (Russell 2004 §4.8 pp.125-127: "The transit times
    and Mars v-inf are calculated using the designated transit leg, as
    indicated by an uppercase descriptor letter"; established for this
    catalogue by #794). NOT always ``arcs[0]``: it is ``arcs[1]`` (G) for
    russell-ch4-5.30gGf3 and ``arcs[2]`` (F) for the ggF/gfF-pattern rows."""
    ups = [i for i, a in enumerate(arcs) if str(a.get("raw_descriptor") or "")[:1].isupper()]
    if len(ups) != 1:
        raise ValueError(
            f"expected exactly one designated (uppercase) arc, found {len(ups)} in "
            f"{[a.get('raw_descriptor') for a in arcs]}"
        )
    return ups[0]


def _g_arc_transfer_angle_deg(raw_descriptor: str) -> float:
    """The printed transfer angle theta (deg) of a generic g/G(t_f, theta, eps)
    descriptor (McConaghy/Russell/Longuski 2005 Eqs (1)-(12))."""
    m = re.match(r"[gG]\(\s*[^,]+,\s*([0-9.eE+-]+)\s*,", raw_descriptor)
    if not m:
        raise ValueError(f"cannot parse transfer angle from {raw_descriptor!r}")
    return float(m.group(1))


def _leg1_topology_candidates(designated: dict[str, Any]) -> list[tuple[int, str]]:
    """Candidate ``(n_revs, branch)`` for the M->E designated-remainder leg.

    full-rev (F) designated: the conic closes after N spacecraft revs in M
    Earth years and the E->M taxi transit is a fraction of one conic period,
    so the remainder sweeps exactly N-1 complete revs.

    generic (G) designated: the full E-E leg sweeps the printed transfer angle
    theta; the E->M piece sweeps an unknown theta_out in (0, 360), so the
    remainder's complete-rev count is floor(theta/360) or one less. The
    low/high Lambert branch of a multi-rev remainder is likewise not printed.
    All candidates are enumerated; :func:`select_leg1_topology` discriminates
    by residual-at-truth (the same criterion #794 used to select Lambert
    roots: the sourced geometry is the discriminator, never a tolerance)."""
    arc_type = designated["arc_type"]
    if arc_type == "full-rev":
        n_total = int(str(designated["resonance"]).split(":")[1])
        n_candidates = [max(0, n_total - 1)]
    elif arc_type == "generic":
        theta = _g_arc_transfer_angle_deg(str(designated["raw_descriptor"]))
        n_full = int(theta // 360.0)
        n_candidates = sorted({n_full, max(0, n_full - 1)}, reverse=True)
    else:
        raise ValueError(f"unsupported designated arc_type {arc_type!r}")
    out: list[tuple[int, str]] = []
    for n in n_candidates:
        if n == 0:
            out.append((0, "single"))
        else:
            out.append((n, "low"))
            out.append((n, "high"))
    return out


def build_genome(row: dict[str, Any]) -> dict[str, Any]:
    """Map a Russell row -> corrector genome (#820 re-posed per #794 semantics).

    Returns sequence, per_leg_revs, per_leg_branch, the designated-leg split
    ToFs (E->M = sourced t_out; M->E = designated-leg remainder), the E-E loop
    seeds from the #794-written-back ``loop-ee-*`` segments, the slack leg
    index, the leg-1 topology candidate set, and the target period in seconds.
    """
    segs = row["trajectory"]["segments"]
    seq = [segs[0]["from"]] + [s["to"] for s in segs]
    n_legs = len(seq) - 1

    # Leg 0 (E->M): the sourced taxi transit t_out (segment out-em).
    tof_em = float(segs[0]["tof_days"])
    loop_segs = segs[2:]  # the #794-written-back loop-ee-* segments, in order

    arcs = row["free_return_arcs"]
    d = _designated_arc_index(arcs)
    designated = arcs[d]
    # Full ToF of the designated E-E leg: printed tof_years for g/G; M Earth
    # years for f/F (the #794-fixed M:N convention, via search/descriptor.py).
    designated_tof_days = arc_tof_seed_days(
        designated["arc_type"],
        tof_years=designated.get("tof_years"),
        resonance=designated.get("resonance"),
    )
    # Leg 1 (M->E): the REMAINDER of the designated conic after the Mars
    # encounter at t_out. The sourced t_in (segment ret-me) is the
    # inbound-crossing taxi transit, ridden on another cycle's leg — it is NOT
    # a per-cycle leg ToF, and using it here (the pre-#820 posing) left the
    # beyond-Mars time unaccounted so the seeds could never tile the period.
    tof_me = designated_tof_days - tof_em
    if tof_me <= 0:
        raise ValueError(f"{row['id']}: designated leg {designated_tof_days:.1f} d <= t_out")

    # Loop arcs = the non-designated arcs in CYCLIC itinerary order after the
    # designated leg (#794: the lowercase legs, in cyclic order after it, ARE
    # the loop-ee-N segments). Cross-checked against the written-back segments;
    # the segments carry the seeds.
    loop_arcs = list(arcs[d + 1 :]) + list(arcs[:d])
    if len(loop_arcs) != len(loop_segs):
        raise ValueError(
            f"{row['id']}: {len(loop_arcs)} non-designated arcs vs "
            f"{len(loop_segs)} loop segments — descriptor/segment mapping out of sync"
        )

    leg1_candidates = _leg1_topology_candidates(designated)
    per_leg_revs: list[int] = [0, leg1_candidates[0][0]]
    per_leg_branch: list[str] = ["single", leg1_candidates[0][1]]
    ee_seeds: list[float] = []
    for arc, seg in zip(loop_arcs, loop_segs, strict=True):
        expected = arc_tof_seed_days(
            arc["arc_type"], tof_years=arc.get("tof_years"), resonance=arc.get("resonance")
        )
        if abs(expected - float(seg["tof_days"])) > 1.0:
            raise ValueError(
                f"{row['id']}: loop arc {arc.get('raw_descriptor')} ToF {expected:.1f} d "
                f"does not match segment {seg['id']} tof_days={seg['tof_days']} — "
                "descriptor/segment cyclic-order mapping out of sync"
            )
        per_leg_revs.append(int(seg["n_revs"]))
        per_leg_branch.append(_SEGMENT_BRANCH_MAP.get(str(seg["branch"]), str(seg["branch"])))
        ee_seeds.append(float(seg["tof_days"]))

    # Target period: k * T_syn(E,M). Use the sourced catalogue years directly
    # (these rows are 2-body E-M; period.years is the sourced repeat period).
    period_sec = float(row["period"]["years"]) * DAYS_PER_JULIAN_YEAR * DAY_S

    # Slack leg = the longest E-E loop (most slack to absorb the period).
    # Free legs = E->M, M->E, plus the non-slack E-E loops. With the #820
    # posing the seeds tile the period (designated + loops = k synodic), so
    # the reconstructed slack ToF now agrees with its own sourced seed.
    all_seeds = [tof_em, tof_me, *ee_seeds]
    # slack must be an E-E loop (index >= 2); pick the longest among them.
    ee_indices = list(range(2, n_legs))
    slack_leg = max(ee_indices, key=lambda i: all_seeds[i])
    free_tof = [all_seeds[i] for i in range(n_legs) if i != slack_leg]

    return {
        "sequence": tuple(seq),
        "per_leg_revs": tuple(per_leg_revs),
        "per_leg_branch": tuple(per_leg_branch),
        "free_tof": free_tof,
        "slack_leg": slack_leg,
        "period_sec": period_sec,
        "all_seeds": all_seeds,
        "designated_index": d,
        "designated_raw": str(designated.get("raw_descriptor")),
        "designated_tof_days": designated_tof_days,
        "leg1_candidates": leg1_candidates,
    }


def _genome_with_leg1(genome: dict[str, Any], n_revs: int, branch: str) -> dict[str, Any]:
    """Copy of ``genome`` with leg 1's (n_revs, branch) replaced."""
    g = dict(genome)
    revs = list(genome["per_leg_revs"])
    branches = list(genome["per_leg_branch"])
    revs[1] = n_revs
    branches[1] = branch
    g["per_leg_revs"] = tuple(revs)
    g["per_leg_branch"] = tuple(branches)
    return g


def select_leg1_topology(
    genome: dict[str, Any], *, model: str, phase_epochs: int, t0_center: float
) -> dict[str, Any]:
    """Discriminate leg 1's ``(n_revs, branch)`` by residual-at-truth.

    For each candidate, scan t0 over one period with the free ToFs PINNED at
    the truth seed and record the minimum residual; the candidate whose truth
    residual is lowest is the topology under which the sourced geometry is
    best representable (no tolerance is loosened — the losing candidates are
    reported alongside)."""
    truth_free = _truth_seed(genome)
    period_sec = genome["period_sec"]
    offsets = np.linspace(0.0, period_sec, phase_epochs, endpoint=False)
    candidates_report: list[dict[str, Any]] = []
    best: tuple[float, float, dict[str, Any], tuple[int, str]] | None = None
    for n_revs, branch in genome["leg1_candidates"]:
        gv = _genome_with_leg1(genome, n_revs, branch)
        c_res, c_t0 = float("inf"), t0_center
        for off in offsets:
            t0 = t0_center + float(off)
            try:
                r = _residual_at(gv, t0, truth_free, model, "magnitude")
            except Exception:
                continue
            if r < c_res:
                c_res, c_t0 = r, t0
        candidates_report.append(
            {"n_revs": n_revs, "branch": branch, "truth_residual_kms": round(c_res, 4)}
        )
        if best is None or c_res < best[0]:
            best = (c_res, c_t0, gv, (n_revs, branch))
    assert best is not None
    return {
        "genome": best[2],
        "leg1": {"n_revs": best[3][0], "branch": best[3][1]},
        "best_t0_sec": best[1],
        "best_truth_residual_kms": best[0],
        "candidates": candidates_report,
    }


def _run_one_epoch(args: tuple) -> dict[str, Any]:
    """Worker: run the corrector at one epoch in one residual mode."""
    (genome, t0_sec, residual_mode, model) = args
    ephem = Ephemeris(model)
    try:
        r = ballistic_correct(
            sequence=genome["sequence"],
            per_leg_revs=genome["per_leg_revs"],
            per_leg_branch=genome["per_leg_branch"],
            t0_seed_sec=float(t0_sec),
            tof_seed_days=genome["free_tof"],
            period_sec=genome["period_sec"],
            ephem=ephem,
            vinf_cap=VINF_CAP_KMS,
            slack_leg=genome["slack_leg"],
            tol_kms=CORRECTOR_TOL_KMS,
            residual_mode=residual_mode,
        )
    except Exception as exc:  # surface, don't crash the sweep
        return {
            "t0_sec": float(t0_sec),
            "mode": residual_mode,
            "converged": False,
            "error": repr(exc),
        }
    return {
        "t0_sec": float(t0_sec),
        "mode": residual_mode,
        "converged": bool(r.converged),
        "max_residual_kms": float(r.max_residual_kms),
        "vinf_per_encounter_kms": list(r.vinf_per_encounter_kms),
        "bend_feasible": bool(r.bend_feasible),
        "vinf_cap_ok": bool(r.vinf_cap_ok),
        "tof_days": list(r.tof_days),
    }


def _classify(row: dict[str, Any], genome: dict[str, Any], best: dict[str, Any] | None) -> dict:
    """Compare a closed solution to the row's SOURCED anchors -> outcome."""
    rid = row["id"]
    sourced_vinf = {e["body"]: e["vinf_kms"] for e in (row.get("vinf_kms_at_encounters") or [])}
    inv = row.get("invariants") or {}
    sourced_tr = inv.get("turn_ratio")
    sourced_ar = inv.get("aphelion_ratio")
    sourced_transits = inv.get("transit_times_days")
    sourced_period = float(row["period"]["years"])

    if best is None or not best.get("converged"):
        return {
            "id": rid,
            "outcome": "NO-CLOSE",
            "detail": "corrector did not reach the residual floor at any scanned epoch/mode",
        }

    seq = genome["sequence"]
    vinf = best["vinf_per_encounter_kms"]
    # Per-body achieved V∞ (max over encounters of that body — the binding one).
    achieved_vinf: dict[str, float] = {}
    for body, vi in zip(seq, vinf, strict=False):
        achieved_vinf.setdefault(body, vi)
        achieved_vinf[body] = max(achieved_vinf[body], vi)

    checks: list[str] = []
    vinf_match = True
    for body, src in sourced_vinf.items():
        ach = achieved_vinf.get(body)
        if ach is None:
            vinf_match = False
            checks.append(f"vinf[{body}]: no achieved value")
            continue
        ok = abs(ach - src) <= TOL_VINF_KMS
        vinf_match = vinf_match and ok
        flag = "OK" if ok else "X"
        checks.append(f"vinf[{body}]: src={src:.2f} ach={ach:.2f} d={abs(ach - src):.2f} {flag}")

    # Period: the corrector pins the period as a constraint, so it always matches;
    # report it for completeness.
    achieved_period = sum(best["tof_days"]) / DAYS_PER_JULIAN_YEAR
    period_ok = abs(achieved_period - sourced_period) <= TOL_PERIOD_YEARS
    checks.append(
        f"period: src={sourced_period:.2f} ach={achieved_period:.2f} {'OK' if period_ok else 'X'}"
    )

    # Invariants: t_out (leg 0) + the designated-leg total ToF (legs 0+1).
    # Under the #820 posing leg 1 is the beyond-Mars REMAINDER of the
    # designated conic, so the sourced t_in (the inbound-crossing taxi
    # transit, ridden on another cycle's leg) has no leg to compare against;
    # the designated leg's printed full ToF (Russell Tables 4.9-4.13, via the
    # descriptor) is the sourced anchor checked instead.
    transit_ok = None
    if sourced_transits is not None:
        ach_t_out = float(best["tof_days"][0])
        ach_designated = float(best["tof_days"][0]) + float(best["tof_days"][1])
        src_t_out = float(sourced_transits[0])
        src_designated = float(genome["designated_tof_days"])
        t_out_ok = abs(ach_t_out - src_t_out) <= TOL_TRANSIT_DAYS
        desig_ok = abs(ach_designated - src_designated) <= TOL_TRANSIT_DAYS
        transit_ok = t_out_ok and desig_ok
        checks.append(
            f"t_out: src={src_t_out:.0f} ach={ach_t_out:.1f} {'OK' if t_out_ok else 'X'}; "
            f"designated-leg ToF: src={src_designated:.1f} ach={ach_designated:.1f} "
            f"{'OK' if desig_ok else 'X'}"
        )

    matched = vinf_match and period_ok and (transit_ok in (None, True))
    outcome = "CLOSE-AND-MATCH" if matched else "CLOSE-OFF-ANCHOR"
    return {
        "id": rid,
        "outcome": outcome,
        "mode": best["mode"],
        "max_residual_kms": best["max_residual_kms"],
        "bend_feasible": best["bend_feasible"],
        "vinf_cap_ok": best["vinf_cap_ok"],
        "achieved_vinf": {k: round(v, 3) for k, v in achieved_vinf.items()},
        "sourced_vinf": sourced_vinf,
        "sourced_turn_ratio": sourced_tr,
        "sourced_aphelion_ratio": sourced_ar,
        "checks": checks,
    }


def _truth_seed(genome: dict[str, Any]) -> list[float]:
    """The row's OWN sourced ToF geometry as a free-leg seed (slack eliminated).

    ``all_seeds`` are [E->M (sourced t_out), M->E (designated-leg remainder),
    *E-E loops (#794 written-back segments)] — the by-construction
    circular-coplanar geometry, which now tiles the period (#820).
    Drop the slack leg (reconstructed by the corrector as period - sum(free))."""
    all_seeds = genome["all_seeds"]
    slack = genome["slack_leg"]
    return [all_seeds[i] for i in range(len(all_seeds)) if i != slack]


def _residual_at(
    genome: dict[str, Any], t0_sec: float, free_tof: list[float], model: str, residual_mode: str
) -> float:
    """Max residual evaluated AT a fixed genome (no solve)."""
    ephem = Ephemeris(model)
    x = np.array([t0_sec, *free_tof], dtype=np.float64)
    res = _residuals(
        x,
        sequence=genome["sequence"],
        per_leg_revs=genome["per_leg_revs"],
        per_leg_branch=genome["per_leg_branch"],
        slack_leg=genome["slack_leg"],
        period_days=genome["period_sec"] / DAY_S,
        ephem=ephem,
        residual_mode=residual_mode,
    )
    return max(abs(r) for r in res)


def probe_at_truth(row: dict[str, Any], *, phase_epochs: int, model: str) -> dict:
    """Seed-at-truth probe: does the corrector STAY at the sourced geometry?

    1. Fix the free-leg ToFs at the row's OWN sourced values (``_truth_seed``).
    2. Scan t0 over one target period; pick the phase minimising the residual
       evaluated AT that truth geometry (magnitude mode — the continuity floor).
    3. Run the full corrector seeded EXACTLY there (truth ToFs, best-phase t0).
    4. Report whether it stayed (residual->0, ToFs ~ truth) or walked away.
    """
    # Leg-1 topology selection doubles as the truth phase scan: for every
    # (n_revs, branch) candidate it scans t0 with the ToFs PINNED at truth;
    # the winner's scan IS the residual-at-truth landscape minimum.
    sel = select_leg1_topology(
        build_genome(row), model=model, phase_epochs=phase_epochs, t0_center=_t0_center(row)
    )
    genome = sel["genome"]
    truth_free = _truth_seed(genome)
    period_sec = genome["period_sec"]
    best_t0 = float(sel["best_t0_sec"])
    best_truth_res = float(sel["best_truth_residual_kms"])

    # Seed the corrector EXACTLY at (best-phase t0, truth ToFs) and let it run.
    ephem = Ephemeris(model)
    solved = ballistic_correct(
        sequence=genome["sequence"],
        per_leg_revs=genome["per_leg_revs"],
        per_leg_branch=genome["per_leg_branch"],
        t0_seed_sec=best_t0,
        tof_seed_days=truth_free,
        period_sec=period_sec,
        ephem=ephem,
        vinf_cap=VINF_CAP_KMS,
        slack_leg=genome["slack_leg"],
        tol_kms=CORRECTOR_TOL_KMS,
        residual_mode="magnitude",
    )

    # Did it stay? Compare solved free-leg ToFs to the truth seed (the slack leg
    # is reconstructed, so compare the FREE legs the solver actually moves).
    solved_full = list(solved.tof_days)
    slack = genome["slack_leg"]
    solved_free = [solved_full[i] for i in range(len(solved_full)) if i != slack]
    tof_drift_days = max(
        (abs(a - b) for a, b in zip(solved_free, truth_free, strict=False)), default=0.0
    )
    t0_drift_days = abs(solved.t0_sec - best_t0) / DAY_S
    stayed = (
        solved.converged
        and tof_drift_days <= TOL_TRANSIT_DAYS
        and best_truth_res <= CORRECTOR_TOL_KMS
    )
    return {
        "id": row["id"],
        "model": model,
        "leg1": sel["leg1"],
        "leg1_candidates": sel["candidates"],
        "best_phase_truth_residual_kms": round(best_truth_res, 4),
        "truth_residual_below_floor": bool(best_truth_res <= CORRECTOR_TOL_KMS),
        "solved_converged": bool(solved.converged),
        "solved_max_residual_kms": round(float(solved.max_residual_kms), 4),
        "tof_drift_days": round(tof_drift_days, 2),
        "t0_drift_days": round(t0_drift_days, 2),
        "truth_free_tof_days": [round(v, 1) for v in truth_free],
        "solved_free_tof_days": [round(v, 1) for v in solved_free],
        "solved_vinf_per_encounter_kms": [round(v, 3) for v in solved.vinf_per_encounter_kms],
        "verdict": "STAYED-AT-TRUTH" if stayed else "WALKED-AWAY",
    }


def probe_at_truth_free_return(row: dict[str, Any], *, phase_epochs: int, model: str) -> dict:
    """Seed-at-truth probe for the #137 free-return genome.

    The free-return genome's "truth" seed is the SOURCED ellipse ``(a, e)`` derived
    from the sourced aphelion + outbound transit (:func:`_seed_ae_from_aphelion_transit`).
    This probe: (1) fixes ``(a, e)`` at that sourced seed, (2) scans t0 over one
    period for the phase minimising the free-return residual AT that geometry,
    (3) runs :func:`free_return_correct` seeded EXACTLY there, and (4) reports
    whether it STAYED (residual->0, ``(a, e)`` unmoved) — the end-to-end visible
    truth-residual≈0 the acceptance test pins. Honours ``--genome free-return``.
    """
    rid = row["id"]
    aphelion = row["orbit_elements"].get("aphelion_au")
    transit = (row.get("invariants") or {}).get("transit_times_days")
    if aphelion is None or not transit:
        return {"id": rid, "verdict": "NO-SEED", "detail": "missing sourced aphelion or transit"}

    period_sec = float(row["period"]["years"]) * DAYS_PER_JULIAN_YEAR * DAY_S
    a_seed, e_seed = _seed_ae_from_aphelion_transit(float(aphelion), float(transit[0]))
    ephem = Ephemeris(model)

    # Match run_row_free_return's dense phase floor (#137 Part 3): the truth-residual
    # landscape of deep-aphelion high-e rows is narrow in t0.
    n_phase = max(phase_epochs, FR_PHASE_EPOCHS_FLOOR)
    best_t0, best_truth_res = 0.0, float("inf")
    for frac in np.linspace(0.0, 1.0, n_phase, endpoint=False):
        t0 = float(frac) * period_sec
        try:
            res = _fr_residuals(
                np.array([a_seed, e_seed, t0]),
                period_days=period_sec / DAY_S,
                ephem=ephem,
                bodies=("E", "M"),
                mu=132712440018.0,
            )
        except Exception:
            continue
        m = max(abs(r) for r in res)
        if m < best_truth_res:
            best_truth_res, best_t0 = m, t0

    solved = free_return_correct(
        t0_seed_sec=best_t0,
        a_seed_au=a_seed,
        e_seed=e_seed,
        period_sec=period_sec,
        ephem=ephem,
        tol_kms=CORRECTOR_TOL_KMS,
    )
    # (a, e) drift relative to the sourced seed (the free variables the solver moves).
    ae_drift = max(abs(solved.a_au - a_seed), abs(solved.e - e_seed))
    stayed = (
        solved.converged
        and best_truth_res <= CORRECTOR_TOL_KMS
        and ae_drift <= 0.05  # AU / dimensionless: the seed stays on its ridge
    )
    return {
        "id": rid,
        "model": model,
        "genome": "free-return",
        "best_phase_truth_residual_kms": round(best_truth_res, 4),
        "truth_residual_below_floor": bool(best_truth_res <= CORRECTOR_TOL_KMS),
        "solved_converged": bool(solved.converged),
        "solved_max_residual_kms": round(float(solved.max_residual_kms), 4),
        "tof_drift_days": round(abs(solved.transfer_tof_days - float(transit[0])), 2),
        "ae_drift": round(ae_drift, 4),
        "seed_a_au": round(a_seed, 4),
        "seed_e": round(e_seed, 4),
        "solved_a_au": round(solved.a_au, 4),
        "solved_e": round(solved.e, 4),
        "solved_vinf_kms": {k: round(v, 3) for k, v in solved.vinf_kms.items()},
        "verdict": "STAYED-AT-TRUTH" if stayed else "WALKED-AWAY",
    }


def run_row(row: dict[str, Any], *, epochs: int, workers: int, model: str) -> dict:
    t0_center = _t0_center(row)
    # Select leg 1's (n_revs, branch) by residual-at-truth before the epoch
    # grid (the sourced geometry discriminates the unprinted topology).
    sel = select_leg1_topology(
        build_genome(row),
        model=model,
        phase_epochs=LEG1_SELECT_PHASE_EPOCHS,
        t0_center=t0_center,
    )
    genome = sel["genome"]

    # Epoch grid over one target period centred on the priority date.
    period_sec = genome["period_sec"]
    offsets = np.linspace(-0.5 * period_sec, 0.5 * period_sec, epochs)
    tasks: list[tuple] = []
    for off in offsets:
        for mode in ("magnitude", "vector"):
            tasks.append((genome, t0_center + float(off), mode, model))

    results: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for res in ex.map(_run_one_epoch, tasks):
            results.append(res)

    # Best closed result: converged, lowest residual, then lowest max-V∞.
    closed = [r for r in results if r.get("converged")]

    def _key(r: dict[str, Any]) -> tuple:
        return (r["max_residual_kms"], max(r["vinf_per_encounter_kms"], default=1e9))

    best = min(closed, key=_key) if closed else None
    n_closed_mag = sum(1 for r in closed if r["mode"] == "magnitude")
    n_closed_vec = sum(1 for r in closed if r["mode"] == "vector")

    verdict = _classify(row, genome, best)
    verdict["model"] = model
    verdict["n_closed_magnitude"] = n_closed_mag
    verdict["n_closed_vector"] = n_closed_vec
    verdict["genome"] = {
        "sequence": "-".join(genome["sequence"]),
        "per_leg_revs": list(genome["per_leg_revs"]),
        "per_leg_branch": list(genome["per_leg_branch"]),
        "slack_leg": genome["slack_leg"],
        "seeds": [round(s, 1) for s in genome["all_seeds"]],
        "designated": genome["designated_raw"],
        "leg1_selection": {
            "selected": sel["leg1"],
            "truth_residual_kms": round(float(sel["best_truth_residual_kms"]), 4),
            "candidates": sel["candidates"],
        },
    }
    return verdict


# ---------------------------------------------------------------------------
# #137 free-return (radial-crossing) genome path.
#
# The seed (a, e) is derived from TWO SOURCED anchors -- the aphelion (Russell
# Table 4.9-4.13) and the outbound transit ToF (segments) -- so aphelion and
# transit are CONSTRAINTS (imposed), not evidence. The per-body V_inf EMERGES
# from the converged ellipse and is the EVIDENCE compared (non-circularly)
# against the independently sourced V_inf anchors. See the constraint-vs-evidence
# table in docs/notes/2026-06-06-russell12-likeforlike.md (#137 section).
# ---------------------------------------------------------------------------


def _seed_ae_from_aphelion_transit(aphelion_au: float, transit_days: float) -> tuple[float, float]:
    """Thin wrapper kept for the campaign's existing call sites.

    The derivation now lives in :func:`cyclerfinder.search.free_return
    .seed_ae_from_aphelion_transit` (promoted in task #106 so the discover path and
    this campaign share one implementation — the physics is not duplicated).
    """
    return seed_ae_from_aphelion_transit(aphelion_au, transit_days)


def run_row_free_return(row: dict[str, Any], *, phase_epochs: int, model: str) -> dict:
    """Free-return genome closure for one row (#137).

    Derives the constrained seed ``(a, e)`` from sourced aphelion + transit,
    scans t0 for the best phase, runs :func:`free_return_correct`, then compares
    the EMERGED V_inf to the sourced anchor (the only evidence). Symmetric single
    ellipse only -- the asymmetric (different-per-leg transit) rows are flagged.
    """
    rid = row["id"]
    aphelion = row["orbit_elements"].get("aphelion_au")
    transit = (row.get("invariants") or {}).get("transit_times_days")
    sourced_vinf = {e["body"]: e["vinf_kms"] for e in (row.get("vinf_kms_at_encounters") or [])}
    period_sec = float(row["period"]["years"]) * DAYS_PER_JULIAN_YEAR * DAY_S

    if aphelion is None or not transit:
        return {"id": rid, "outcome": "NO-SEED", "detail": "missing sourced aphelion or transit"}

    asymmetric = len(transit) >= 2 and abs(float(transit[0]) - float(transit[1])) > 1.0
    a_seed, e_seed = _seed_ae_from_aphelion_transit(float(aphelion), float(transit[0]))

    ephem = Ephemeris(model)
    # Deep-aphelion / high-e rows (e.g. russell-ch4-9.353Gg2, e~0.43) have a narrow
    # t0 residual basin that a coarse 256-point grid steps over, leaving the
    # corrector to drift ~0.6 km/s off-anchor (#137 Part 3 straggler diagnosis).
    # The free-return phase scan is pure residual evals (no Lambert), so a dense
    # floor is cheap and makes the best phase robust across all rows.
    n_phase = max(phase_epochs, FR_PHASE_EPOCHS_FLOOR)
    best_t0, best_res = 0.0, float("inf")
    for frac in np.linspace(0.0, 1.0, n_phase, endpoint=False):
        t0 = float(frac) * period_sec
        try:
            res = _fr_residuals(
                np.array([a_seed, e_seed, t0]),
                period_days=period_sec / DAY_S,
                ephem=ephem,
                bodies=("E", "M"),
                mu=132712440018.0,
            )
        except Exception:
            continue
        m = max(abs(r) for r in res)
        if m < best_res:
            best_res, best_t0 = m, t0

    r = free_return_correct(
        t0_seed_sec=best_t0,
        a_seed_au=a_seed,
        e_seed=e_seed,
        period_sec=period_sec,
        ephem=ephem,
        tol_kms=CORRECTOR_TOL_KMS,
    )

    # EVIDENCE: emerged V_inf vs sourced anchor (V_inf was never imposed).
    checks: list[str] = []
    vinf_match = True
    for body, src in sourced_vinf.items():
        ach = r.vinf_kms.get(body)
        if ach is None:
            vinf_match = False
            checks.append(f"vinf[{body}]: no achieved value")
            continue
        ok = abs(ach - src) <= TOL_VINF_KMS
        vinf_match = vinf_match and ok
        checks.append(
            f"vinf[{body}]: src={src:.2f} ach={ach:.2f} d={abs(ach - src):.2f} "
            f"{'OK' if ok else 'X'}"
        )

    if not r.converged:
        outcome = "NO-CLOSE"
    elif asymmetric:
        # A symmetric single ellipse cannot represent a different-per-leg transit;
        # report honestly rather than over-claim.
        outcome = "CLOSE-MATCH-SYMMETRIC-ONLY" if vinf_match else "CLOSE-OFF-ANCHOR"
    else:
        outcome = "CLOSE-AND-MATCH" if vinf_match else "CLOSE-OFF-ANCHOR"

    return {
        "id": rid,
        "outcome": outcome,
        "model": model,
        "genome": "free-return",
        "converged": bool(r.converged),
        "max_residual_kms": round(float(r.max_residual_kms), 4),
        "constrained": {  # imposed inputs (NOT evidence)
            "aphelion_au": float(aphelion),
            "transit_days": float(transit[0]),
            "asymmetric_transit": bool(asymmetric),
        },
        "derived": {  # emerged (evidence / comparable)
            "a_au": round(r.a_au, 4),
            "e": round(r.e, 4),
            "vinf_kms": {k: round(v, 3) for k, v in r.vinf_kms.items()},
            "transfer_tof_days": round(r.transfer_tof_days, 1),
            "ee_interval_days": round(r.ee_interval_days, 1),
        },
        "sourced_vinf": sourced_vinf,
        "best_phase_residual_kms": round(best_res, 4),
        "checks": checks,
    }


def _git_short_head() -> str:
    """``git rev-parse --short HEAD`` captured BY THE CALLER (runlog stays pure).

    Returns the empty string outside a git checkout / on any git error so the
    campaign never crashes on the persistence side-effect.
    """
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return ""
    return out.stdout.strip()


def _verdict_to_record(v: dict[str, Any], genome: str, code_version: str) -> RunRecord:
    """Map one campaign verdict dict -> a :class:`RunRecord` (no info dropped).

    Works for BOTH genome paths: the lambert verdict (``run_row`` -> ``_classify``)
    and the free-return verdict (``run_row_free_return``). The fields each path
    happens to carry land in the matching optional slots; absent ones stay empty.
    """
    achieved = dict(v.get("achieved_vinf") or {})
    if not achieved:
        # Free-return path nests the emerged V∞ under ``derived``.
        achieved = dict((v.get("derived") or {}).get("vinf_kms") or {})

    sourced_anchors: dict[str, Any] = {}
    for key in ("sourced_turn_ratio", "sourced_aphelion_ratio"):
        if v.get(key) is not None:
            sourced_anchors[key] = v[key]
    if v.get("constrained") is not None:
        sourced_anchors["constrained"] = v["constrained"]

    seed: dict[str, Any] = {}
    if v.get("genome") and isinstance(v["genome"], dict):
        seed = dict(v["genome"])  # lambert genome descriptor (sequence/revs/seeds)
    if v.get("derived") is not None:
        seed["derived"] = v["derived"]  # free-return emerged (a, e), tof

    solver_audit: dict[str, Any] = {"checks": v.get("checks", [])}
    for key in ("n_closed_magnitude", "n_closed_vector", "mode", "bend_feasible", "vinf_cap_ok"):
        if key in v:
            solver_audit[key] = v[key]
    if "best_phase_residual_kms" in v:
        solver_audit["best_phase_residual_kms"] = v["best_phase_residual_kms"]

    return RunRecord(
        row_id=str(v.get("id")),
        genome=genome,
        outcome=str(v.get("outcome")),
        model=str(v.get("model")),
        code_version=code_version,
        achieved_vinf_kms=achieved,
        sourced_vinf_kms=dict(v.get("sourced_vinf") or {}),
        sourced_anchors=sourced_anchors,
        seed=seed,
        residual_kms=v.get("max_residual_kms"),
        solver_audit=solver_audit,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--epochs", type=int, default=32)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument(
        "--model", type=str, default="astropy", choices=("astropy", "circular", "inclined-circular")
    )
    ap.add_argument(
        "--probe-at-truth",
        action="store_true",
        help="run the seed-at-truth diagnostic probe per row",
    )
    ap.add_argument(
        "--phase-epochs",
        type=int,
        default=256,
        help="t0 phase-scan points for the seed-at-truth probe",
    )
    ap.add_argument(
        "--genome",
        type=str,
        default="lambert",
        choices=("lambert", "free-return"),
        help="lambert = #125/#135 free-Lambert genome (default, byte-unchanged); "
        "free-return = #137 radial-crossing genome (seed (a,e) from sourced "
        "aphelion+transit, V_inf emerges as evidence)",
    )
    ap.add_argument(
        "--rows",
        type=str,
        default=None,
        help="comma-separated subset of RUSSELL12_IDS to run (default: all 12); "
        "added for the #813 blast-radius re-run of the 4 non-1:1 f/F rows",
    )
    ap.add_argument("--out", type=str, default=None)
    ap.add_argument(
        "--runlog-dir",
        type=str,
        default=str(DEFAULT_RUNLOG_DIR),
        help="directory for the per-invocation JSON-lines runlog (default ON, "
        "writes data/runs/russell12-<timestamp>.jsonl); pass empty to disable",
    )
    ap.add_argument(
        "--runlog-timestamp",
        type=str,
        default=None,
        help="filename timestamp tag for the runlog (default: UTC now)",
    )
    args = ap.parse_args()

    rows = yaml.safe_load(CATALOGUE_PATH.read_text())
    byid = {r["id"]: r for r in rows}

    run_ids: tuple[str, ...] = RUSSELL12_IDS
    if args.rows:
        requested = tuple(s.strip() for s in args.rows.split(",") if s.strip())
        unknown = [r for r in requested if r not in RUSSELL12_IDS]
        if unknown:
            raise SystemExit(f"--rows ids not in RUSSELL12_IDS: {unknown}")
        run_ids = requested

    print(
        f"genome={args.genome} model={args.model} epochs={args.epochs} "
        f"workers={args.workers} probe_at_truth={args.probe_at_truth}",
        flush=True,
    )

    verdicts: list[dict] = []
    for rid in run_ids:
        print(f"=== {rid} ===", flush=True)
        if args.genome == "free-return":
            v = run_row_free_return(byid[rid], phase_epochs=args.phase_epochs, model=args.model)
            verdicts.append(v)
            print(
                f"  {v['outcome']}  res={v.get('max_residual_kms')}",
                flush=True,
            )
            for c in v.get("checks", []):
                print(f"    {c}", flush=True)
            continue
        v = run_row(byid[rid], epochs=args.epochs, workers=args.workers, model=args.model)
        verdicts.append(v)
        print(
            f"  {v['outcome']}  closed(mag/vec)={v['n_closed_magnitude']}/{v['n_closed_vector']}",
            flush=True,
        )
        for c in v.get("checks", []):
            print(f"    {c}", flush=True)

    probes: list[dict] = []
    if args.probe_at_truth:
        print(f"\n=== SEED-AT-TRUTH PROBE (genome={args.genome}) ===", flush=True)
        for rid in run_ids:
            # Honour --genome: the free-return probe seeds the SOURCED ellipse (a, e)
            # and shows the truth-residual≈0 end-to-end; the lambert probe seeds the
            # sourced ToF geometry (the #135 diagnostic). Previously the probe always
            # ran the lambert genome regardless of --genome (the results-note caveat).
            if args.genome == "free-return":
                p = probe_at_truth_free_return(
                    byid[rid], phase_epochs=args.phase_epochs, model=args.model
                )
            else:
                p = probe_at_truth(byid[rid], phase_epochs=args.phase_epochs, model=args.model)
            probes.append(p)
            print(
                f"{rid:24s} {p['verdict']:16s} "
                f"truth_res={p.get('best_phase_truth_residual_kms', float('nan')):.3f} "
                f"solved_res={p.get('solved_max_residual_kms', float('nan')):.3f} "
                f"tof_drift={p.get('tof_drift_days')}d",
                flush=True,
            )

    # Summary table.
    print("\n=== SUMMARY ===")
    counts: dict[str, int] = {}
    for v in verdicts:
        counts[v["outcome"]] = counts.get(v["outcome"], 0) + 1
        print(f"{v['id']:30s} {v['outcome']}")
    print("\ncounts:", counts)
    if probes:
        pcounts: dict[str, int] = {}
        for p in probes:
            pcounts[p["verdict"]] = pcounts.get(p["verdict"], 0) + 1
        print("probe counts:", pcounts)

    if args.out:
        Path(args.out).write_text(
            json.dumps({"model": args.model, "verdicts": verdicts, "probes": probes}, indent=2)
        )
        print(f"wrote {args.out}")

    # Persist the per-row (seed -> basin) labels the campaign just computed
    # (additive side-effect; does not change the scientific run or stdout above).
    if args.runlog_dir:
        code_version = _git_short_head()
        runlog_path = default_runlog_path(
            args.runlog_dir, f"russell12-{args.model}", timestamp=args.runlog_timestamp
        )
        runlog = RunLog(runlog_path)
        runlog.extend([_verdict_to_record(v, args.genome, code_version) for v in verdicts])
        print(f"wrote runlog {runlog_path} ({len(verdicts)} records)")


if __name__ == "__main__":
    from cyclerfinder.search.outcome_log import enable_default_outcome_log

    enable_default_outcome_log("campaign_russell12")
    main()
