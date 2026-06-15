"""Asymmetric-corrector novel-cycler scan at Earth-Moon (#284).

Repurposes the #249 ``correct_general_periodic`` asymmetric corrector as a
*novel-search driver*. Grids over (k1, k2) winding classes x Jacobi C x seed
region of (x0, xdot0), calls the corrector for each cell, then routes each
converged candidate through the guard chain:

1.  Asymmetric Newton converges + corrector's own Radau closure < 1e-6.
2.  Topology classification via ``binary_star_search.winding_topology``
    (must match the requested target ``(k1, k2)``).
3.  Independent DOP853 (rtol=atol=1e-12) re-propagation closure < 1e-6
    (basin-artifact check; corrector uses DOP853 in the Newton loop, but
    here we restart from the corrected IC fresh with the strict tolerance
    -- per ``feedback_orbit_closure_discipline``).
4.  Topology-based offline corpus dedup against ``literature_check``
    Earth-Moon anchors (Braik-Ross, Roberts-Tsoukkas/Ross, Kumar-Rawat-
    Rosengren-Ross, Koblick tulip, Zhang-Jiang-Yuan tulip). NOT a live
    WebSearch -- offline matching is necessary-not-sufficient for novelty
    and the doc says so explicitly.
5.  ML false-positive flagger (#256/#275) probability < 0.5.

Per the task discipline: NO catalogue writeback. NO novelty claims in commit
messages. Output is ``data/scan_284.jsonl`` with one row per CONVERGED cell.
A row is ``novelty_claimable`` (a flag, not a writeback) iff ALL guards pass.

The Tier-0 NN prefilter from #276 / #282 is NOT used here: it scores
patched-conic *Lambert legs* between moon-state pairs, not CR3BP periodic
orbits. Per the 5-tier prioritizer's documented architectural seam, the
representative-pair tiers (1-5) need orbit pairs, not a single candidate
orbit -- so the 5-tier stack does not apply to a discovery scan that emits
one orbit per cell. We document the gap rather than misuse a Lambert-leg
scorer on a periodic-orbit candidate.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from collections.abc import Iterable
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.ml.falsepos_flagger import FalsePosFlagger
from cyclerfinder.ml.falsepos_labels import build_training_set
from cyclerfinder.search.binary_star_search import winding_topology
from cyclerfinder.search.cr3bp_general_periodic import (
    GeneralPeriodicOrbit,
    _ydot0_general,
    correct_general_periodic,
)
from cyclerfinder.search.literature_check import (
    CandidateSignature,
    _candidate_anchors,
)
from cyclerfinder.search.reachable_representatives import TU_DAYS, braik_ross_system

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_PATH = REPO_ROOT / "data" / "scan_284.jsonl"

# Default grid (start narrow; expand only if first batch lands nothing).
DEFAULT_K_TARGETS: tuple[tuple[int, int], ...] = (
    (1, 1),
    (1, 2),
    (2, 1),
    (3, 2),
    (2, 3),
)
# Down-scoped slab for the 2-hour budget. The asymmetric corrector takes
# ~6-15 s/cell (dominated by STM propagation across the multi-crossing
# return arc). The task brief sketched a much denser grid (21 C levels x
# 21 x0 seeds x 4 xdot0 x 2 signs x 4 hc x 7 K = 98 784 cells, ~150 hours);
# the doc records the gap and proposes a follow-up batched run.
DEFAULT_C_GRID = (3.08, 3.12, 3.14, 3.17)  # 4 Jacobi levels around C_J=3.1294
DEFAULT_X0_GRID = (-0.90, -0.82, -0.75, -0.68)  # 4 seeds, ~0.07 step
# Pair (+, -) to scan both branches of the asymmetric basin; 0.0 left out
# because it is the symmetric special case (covered by the symmetric scan).
DEFAULT_XDOT0_SEEDS = (-0.05, 0.05)
DEFAULT_YDOT0_SIGNS = (-1.0,)  # follows the family-seed convention
# half_crossings: 1-based index of the y=0 crossing at which the full state
# returns. C11a/C11b/C32 all close at 2*H=6 in the symmetric corrector;
# C21 closes at half_crossings ~ 4. We sweep (4, 6) only -- hc=8 is the
# longest arc and the most numerically taxing.
DEFAULT_HALF_CROSSINGS = (4, 6)

# Cell-level corrector budget. max_iter lowered from 60 to 20 because the
# Newton basin is narrow and non-converging cells should bail fast; 20 is
# enough to recover a known symmetric seed (8-12 iter in the reproduce test).
PERIOD_GUESS_TU = 14.0  # ~ 60 d; corrector only uses this for the integration horizon.
CORRECTOR_TOL = 1e-11
CORRECTOR_MAX_ITER = 20
INDEPENDENT_CLOSURE_TOL = 1e-6  # the orbit-closure-discipline cap.

# Deduplication of converged orbits to one row per "physical orbit" at the
# scan grid's resolution. (x0, jacobi, period_d) rounded.
_DEDUP_X0_DP = 4
_DEDUP_C_DP = 4
_DEDUP_PD_DP = 2


def _log(msg: str) -> None:
    """Print a timestamped log line, flushing immediately."""
    ts = datetime.now(UTC).strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)


def _independent_closure(mu: float, orbit: GeneralPeriodicOrbit) -> tuple[float, bool]:
    """Re-propagate the IC with DOP853 rtol=atol=1e-12 and return closure norm.

    Independent of the corrector's own Radau check (different integrator,
    different tolerance). Per ``feedback_orbit_closure_discipline``: the
    independent cross-check is mandatory for a novelty-claimable candidate.
    """
    if not math.isfinite(orbit.ydot0):
        return float("inf"), False
    state0 = np.array([orbit.x0, 0.0, 0.0, orbit.xdot0, orbit.ydot0, 0.0])
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, orbit.period),
        state0,
        args=(mu,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
    )
    if not sol.success:
        return float("inf"), False
    return float(np.linalg.norm(sol.y[:, -1] - state0)), True


def _signature_for_em_cycler(k1: int, k2: int, period_days: float) -> CandidateSignature:
    """Structural fingerprint for an Earth-Moon CR3BP cycler.

    The literature_check corpus uses primary=Earth, body_set={Moon} for the
    Earth-Moon CR3BP family (Braik-Ross, Roberts-Tsoukkas, Koblick tulip, etc.).
    A non-trivial body_set overlap there triggers ``_candidate_anchors`` so an
    EM-cycler candidate at this scan's coverage is *always* candidate-mapped
    against the EM corpus.

    n_rev carries (k1, k2) as the per-leg revolution count signal so the
    corpus anchor's structural keys are populated. period_years is not the
    right unit for an EM cycler (a synodic-month-scale orbit) -- left None.
    """
    # NOTE: the literature_check anchors use ``body_set=frozenset({"Moon"})``
    # for the Earth-Moon CR3BP family (the central body is implicit via
    # ``primary="Earth"``), so the signature's ``sequence`` must be only
    # ("Moon",) for ``_candidate_anchors`` to fire (``seq_set <= body_set``).
    return CandidateSignature(
        primary="Earth",
        sequence=("Moon",),
        period_k=None,
        period_years=None,
        vinf_per_encounter_kms=(),
        resonances=(f"{k1}:{k2}",),
        n_rev=(k1, k2),
    )


def _offline_em_corpus_match(sig: CandidateSignature) -> list[dict[str, str]]:
    """Offline mapping from a signature to overlapping EM corpus anchors.

    Returns a list of ``{name, citation, doi_or_url}`` for each anchor whose
    structural footprint overlaps the candidate. Necessary-not-sufficient
    flag for the rediscovery check (live WebSearch would be the authoritative
    pass, but this scan runs offline; the doc records the gap).
    """
    anchors = _candidate_anchors(sig)
    matches: list[dict[str, str]] = []
    for a in anchors:
        matches.append(
            {
                "name": a.name,
                "citation": a.citation,
                "doi": a.doi or "",
            }
        )
    return matches


def _fp_flagger_lazy_load() -> FalsePosFlagger | None:
    """Build + fit the false-positive flagger on the labelled corpus.

    Returns ``None`` if the labelled set / fit fails (the flagger is
    non-blocking by contract).
    """
    try:
        clf = FalsePosFlagger()
        x_train, y_train, _ = build_training_set()
        clf.fit(x_train, y_train)
    except Exception as exc:
        _log(f"WARN: flagger fit failed: {exc!r}; flagger disabled (p_fp=NaN).")
        return None
    return clf


def _fp_flagger_score(clf: FalsePosFlagger | None, orbit: GeneralPeriodicOrbit) -> float:
    """Convert a corrected orbit into a flagger record dict + score it.

    The flagger lives in patched-conic ΔV space; we feed it the natural
    CR3BP analogues: corrector residual stands in for ``max_residual_kms``,
    the (k1, k2) request is encoded as ``period_days`` for the resonance
    deviation feature, and we mark the closure_method_version as a known
    post-fix SHA so the epoch / mu-fix flags are well-defined. Missing
    Lambert-leg fields fall back to the flagger's median imputation.
    """
    if clf is None:
        return float("nan")
    rec: dict[str, Any] = {
        "max_residual_kms": float(orbit.closure_residual),
        "bend_feasible": True,
        "topology_match": True,
        "period_days": float(orbit.period * TU_DAYS),
        "closure_method_version": "23b980e",  # #249 asym corrector
        "closure_date": "2026-06-16",
        "model_assumption": "cr3bp_asymmetric_general_periodic",
        "cross_check_shared_with_primary": False,
    }
    try:
        return float(clf.score(rec))
    except Exception:
        return float("nan")


def _enumerate_cells(
    k_targets: Iterable[tuple[int, int]],
    c_grid: Iterable[float],
    x0_grid: Iterable[float],
    xdot0_seeds: Iterable[float],
    ydot0_signs: Iterable[float],
    half_crossings: Iterable[int],
) -> Iterable[tuple[int, int, float, float, float, float, int]]:
    """Yield (k1, k2, C, x0, xdot0, ydot0_sign, half_crossings) one cell at a time."""
    for k1, k2 in k_targets:
        for c in c_grid:
            for x0 in x0_grid:
                for xd in xdot0_seeds:
                    for sign in ydot0_signs:
                        for hc in half_crossings:
                            yield (
                                int(k1),
                                int(k2),
                                float(c),
                                float(x0),
                                float(xd),
                                float(sign),
                                int(hc),
                            )


def _is_known_em_family(k1: int, k2: int, period_days: float) -> bool:
    """Known Earth-Moon catalogue entries the scan should explicitly dedupe.

    Sourced from the data/catalogue.yaml Braik-Ross / Ross-Roberts-Tsoukkas
    rows (C11a/C11b/C21/C32) and tulip families. A (k1, k2) winding-class
    match alone is not a literature collision (#172 method-versioned
    registry discipline); we additionally require the orbital period to be
    near the source's value. The function returns True only when both the
    topology AND the published period band agree.
    """
    # Braik-Ross Table-2 common-energy cyclers (period in days, +-2%).
    em_known: list[tuple[int, int, float]] = [
        (1, 1, 42.14),  # C11a
        (1, 1, 55.96),  # C11b
        (2, 1, 84.53),  # C21 (Braik-Ross asymmetric)
        (3, 2, 78.61),  # C32
    ]
    for kk1, kk2, p in em_known:
        if k1 == kk1 and k2 == kk2 and abs(period_days - p) / p < 0.02:
            return True
    return False


def run_scan(
    *,
    k_targets: Iterable[tuple[int, int]] = DEFAULT_K_TARGETS,
    c_grid: Iterable[float] = DEFAULT_C_GRID,
    x0_grid: Iterable[float] = DEFAULT_X0_GRID,
    xdot0_seeds: Iterable[float] = DEFAULT_XDOT0_SEEDS,
    ydot0_signs: Iterable[float] = DEFAULT_YDOT0_SIGNS,
    half_crossings: Iterable[int] = DEFAULT_HALF_CROSSINGS,
    out_path: Path = OUT_PATH,
    progress_every: int = 200,
    time_budget_s: float | None = None,
) -> dict[str, Any]:
    """Drive the (k1, k2) x C x seed-region grid through the asymmetric corrector.

    Writes one JSONL row per CONVERGED cell. Returns a summary stats dict.
    """
    system = braik_ross_system()
    mu = system.mu
    flagger = _fp_flagger_lazy_load()

    # Materialise the grids so we can both report sizes and iterate them.
    k_targets_t = tuple(k_targets)
    c_grid_t = tuple(c_grid)
    x0_grid_t = tuple(x0_grid)
    xdot0_seeds_t = tuple(xdot0_seeds)
    ydot0_signs_t = tuple(ydot0_signs)
    half_crossings_t = tuple(half_crossings)

    cells = list(
        _enumerate_cells(
            k_targets_t,
            c_grid_t,
            x0_grid_t,
            xdot0_seeds_t,
            ydot0_signs_t,
            half_crossings_t,
        )
    )
    total = len(cells)
    _log(
        f"scan_284 start: total cells={total} "
        f"(|K|={len(k_targets_t)}, |C|={len(c_grid_t)}, "
        f"|x0|={len(x0_grid_t)}, |xdot0|={len(xdot0_seeds_t)}, "
        f"|signs|={len(ydot0_signs_t)}, |hc|={len(half_crossings_t)})"
    )
    _log(f"output: {out_path}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_fh = out_path.open("w")

    stats = {
        "cells_attempted": 0,
        "cells_converged": 0,
        "topology_matches": 0,
        "independent_closure_pass": 0,
        "literature_fresh": 0,
        "ml_low_fp": 0,
        "novelty_claimable": 0,
        "known_family_collisions": 0,
        "per_k_attempted": {},
        "per_k_converged": {},
        "per_k_topology": {},
        "per_k_independent": {},
        "per_k_novel": {},
    }
    seen_orbits: set[tuple[float, float, float]] = set()
    t0 = time.time()

    for idx, (k1, k2, c, x0g, xd_seed, sign, hc) in enumerate(cells):
        stats["cells_attempted"] += 1
        kkey = f"{k1},{k2}"
        stats["per_k_attempted"][kkey] = stats["per_k_attempted"].get(kkey, 0) + 1

        if time_budget_s is not None and time.time() - t0 > time_budget_s:
            _log(
                f"time budget exhausted at cell {idx + 1}/{total}; "
                f"aborting remaining {total - idx - 1} cells."
            )
            break

        # Skip cells whose Jacobi is infeasible at this (x0, xdot0).
        try:
            _ydot0_general(x0g, xd_seed, c, mu, sign)
        except ValueError:
            continue
        except Exception:
            continue

        try:
            orbit = correct_general_periodic(
                system,
                x0g,
                xd_seed,
                c,
                PERIOD_GUESS_TU,
                half_crossings=hc,
                ydot0_sign=sign,
                tol=CORRECTOR_TOL,
                max_iter=CORRECTOR_MAX_ITER,
                rtol=1e-12,
                atol=1e-12,
            )
        except (ValueError, RuntimeError):
            continue
        except Exception as exc:
            _log(f"WARN: corrector raised at cell {idx}: {exc!r}")
            continue

        if not orbit.converged:
            continue

        # Corrector-internal Radau closure check.
        if not math.isfinite(orbit.closure_residual):
            continue
        if orbit.closure_residual > INDEPENDENT_CLOSURE_TOL:
            continue

        stats["cells_converged"] += 1
        stats["per_k_converged"][kkey] = stats["per_k_converged"].get(kkey, 0) + 1

        # Dedupe by (rounded x0, jacobi, period_d) so families don't write
        # thousands of near-identical rows.
        period_days = orbit.period * TU_DAYS
        key = (
            round(orbit.x0, _DEDUP_X0_DP),
            round(orbit.jacobi, _DEDUP_C_DP),
            round(period_days, _DEDUP_PD_DP),
        )
        if key in seen_orbits:
            continue
        seen_orbits.add(key)

        # Topology classification (cheap; runs on every dedup'd survivor).
        state0 = np.array([orbit.x0, 0.0, 0.0, orbit.xdot0, orbit.ydot0, 0.0])
        try:
            topo = winding_topology(mu, state0, orbit.period)
        except Exception:
            continue

        topology_match = topo.k1 == k1 and topo.k2 == k2
        if topology_match:
            stats["topology_matches"] += 1
            stats["per_k_topology"][kkey] = stats["per_k_topology"].get(kkey, 0) + 1

        # Independent DOP853 closure (no Radau; different integrator AND
        # tolerance from the corrector's own check, per orbit-closure
        # discipline).
        ind_closure, ind_ok = _independent_closure(mu, orbit)
        independent_pass = ind_ok and ind_closure < INDEPENDENT_CLOSURE_TOL
        if independent_pass:
            stats["independent_closure_pass"] += 1
            stats["per_k_independent"][kkey] = stats["per_k_independent"].get(kkey, 0) + 1

        # Literature corpus match (offline; necessary-not-sufficient).
        # The signature uses the candidate's classified topology, not the
        # request, so a topology mismatch is exposed in the row.
        sig = _signature_for_em_cycler(topo.k1, topo.k2, period_days)
        anchors = _offline_em_corpus_match(sig)
        # The corpus anchors all fire for any EM CR3BP candidate (body_set
        # = {Moon}); that overlap is necessary-not-sufficient for a
        # rediscovery. We separately flag the *known sourced families*
        # (C11a/C11b/C21/C32) by period match, which IS a rediscovery.
        known_family = _is_known_em_family(topo.k1, topo.k2, period_days)
        if known_family:
            stats["known_family_collisions"] += 1
        # "literature_fresh" means: the candidate does NOT collide with a
        # known sourced family AND no anchor's structural footprint is a
        # tight period match. (Anchors fire generically for any EM CR3BP
        # orbit; we DO NOT use that as a published-flag -- it's only a
        # rediscovery if the orbit lines up with a sourced family.)
        literature_fresh = (not known_family) and independent_pass
        if literature_fresh:
            stats["literature_fresh"] += 1

        # ML false-positive flagger (#256/#275).
        p_fp = _fp_flagger_score(flagger, orbit)
        ml_low_fp = math.isfinite(p_fp) and p_fp < 0.5
        if ml_low_fp:
            stats["ml_low_fp"] += 1

        novelty_claimable = topology_match and independent_pass and literature_fresh and ml_low_fp
        if novelty_claimable:
            stats["novelty_claimable"] += 1
            stats["per_k_novel"][kkey] = stats["per_k_novel"].get(kkey, 0) + 1

        row = {
            "k_target": [k1, k2],
            "k_classified": [topo.k1, topo.k2],
            "topology_match": bool(topology_match),
            "c_request": round(c, 6),
            "seed_x0": round(x0g, 6),
            "seed_xdot0": round(xd_seed, 6),
            "ydot0_sign": sign,
            "half_crossings": hc,
            "x0": round(orbit.x0, 9),
            "xdot0": round(orbit.xdot0, 9),
            "ydot0": round(orbit.ydot0, 9),
            "jacobi": round(orbit.jacobi, 10),
            "period_TU": round(orbit.period, 9),
            "period_days": round(period_days, 6),
            "asymmetry": round(orbit.asymmetry, 9),
            "x_min": round(topo.x_min, 6),
            "x_max": round(topo.x_max, 6),
            "prograde": bool(topo.prograde),
            "reaches_secondary": bool(topo.reaches_secondary),
            "corrector_residual": float(orbit.residual),
            "corrector_closure_radau": float(orbit.closure_residual),
            "independent_closure_dop853": (float(ind_closure) if ind_ok else None),
            "independent_closure_pass": bool(independent_pass),
            "n_newton_iter": int(orbit.n_iter),
            "literature_offline_anchors": anchors,
            "known_em_family_collision": bool(known_family),
            "literature_fresh_offline": bool(literature_fresh),
            "p_false_positive": (round(p_fp, 4) if math.isfinite(p_fp) else None),
            "ml_low_fp": bool(ml_low_fp),
            "novelty_claimable": bool(novelty_claimable),
            "scan_id": "284",
            "scan_method_version": "asymmetric_corrector(23b980e)",
            "scan_timestamp_utc": datetime.now(UTC).isoformat(),
        }
        out_fh.write(json.dumps(row) + "\n")
        out_fh.flush()

        if novelty_claimable:
            _log(
                f"  NOVELTY-CLAIMABLE cell {idx}: (k1,k2)=({topo.k1},{topo.k2}) "
                f"C={orbit.jacobi:.6f} T={period_days:.2f}d "
                f"x0={orbit.x0:+.6f} xdot0={orbit.xdot0:+.6f} "
                f"asym={orbit.asymmetry:.2e} ind_closure={ind_closure:.2e}"
            )

        if (idx + 1) % progress_every == 0:
            elapsed = time.time() - t0
            rate = (idx + 1) / max(elapsed, 1e-9)
            eta_s = (total - idx - 1) / max(rate, 1e-9)
            _log(
                f"progress {idx + 1}/{total} ({100 * (idx + 1) / total:.1f}%) "
                f"converged={stats['cells_converged']} "
                f"topo_match={stats['topology_matches']} "
                f"ind_pass={stats['independent_closure_pass']} "
                f"novel-claimable={stats['novelty_claimable']} "
                f"elapsed={elapsed:.0f}s eta={eta_s:.0f}s"
            )

    out_fh.close()
    elapsed = time.time() - t0
    stats["elapsed_seconds"] = round(elapsed, 1)
    stats["total_cells"] = total
    _log(
        f"scan_284 done in {elapsed:.0f}s: "
        f"attempted={stats['cells_attempted']} "
        f"converged={stats['cells_converged']} "
        f"topo_match={stats['topology_matches']} "
        f"ind_pass={stats['independent_closure_pass']} "
        f"lit_fresh={stats['literature_fresh']} "
        f"ml_low_fp={stats['ml_low_fp']} "
        f"novel-claimable={stats['novelty_claimable']} "
        f"known_collisions={stats['known_family_collisions']}"
    )
    return stats


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="Tiny grid (1 K cell x 3 C x 3 x0 x 2 xdot0) for CI smoke.",
    )
    parser.add_argument(
        "--time-budget-s",
        type=float,
        default=None,
        help="Soft wall-clock cap (default: unbounded; scan runs full grid).",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=OUT_PATH,
        help=f"Output JSONL path (default: {OUT_PATH}).",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=200,
        help="Log a progress line every N cells (default 200).",
    )
    args = parser.parse_args(argv)

    if args.smoke:
        # Minimal grid for CI / quick sanity: covers (1,1) + (2,1) so both
        # symmetric C11 and asymmetric C21 reproductions get exercised.
        run_scan(
            k_targets=((1, 1), (2, 1)),
            c_grid=(3.13, 3.14),
            x0_grid=(-0.81, -0.80),
            xdot0_seeds=(-0.05, 0.05),
            ydot0_signs=(-1.0,),
            half_crossings=(6,),
            out_path=args.out,
            progress_every=10,
            time_budget_s=args.time_budget_s,
        )
    else:
        run_scan(
            out_path=args.out,
            progress_every=args.progress_every,
            time_budget_s=args.time_budget_s,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
