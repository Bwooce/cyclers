"""#567 step (3) -- generalized epoch-robustness scan across 6 candidates.

Background
----------
#566 gauntlet-tested 5 new #563 symmetric-closure representatives (plus the
already-catalogued #312 Umbriel-Oberon-Umbriel SILVER) through
V2->V3->V4->V4-strict at a SINGLE launch epoch (2000-06-21, #312's own
known-favorable #338 anchor). The #566 Opus adjudication (see
``docs/notes/2026-07-11-566-opus-adjudication-gauntlet-results.md``) held
catalogue writeback pending a trustworthy per-candidate epoch-robustness
scan, because:

* #312 is KNOWN to vary across epochs (#338 found 2000-01-15 FAIL,
  2000-06-21 PASS, same year).
* #559 proved #312's ~10-14% daily V4-strict FAILs under the PRE-FIX code
  were CONFIRMED numerical artifacts -- a Lambert branch-selection flip and
  a DOP853 stiff-death on genuinely planet-crossing arcs silently
  misclassified as an unexplained integrator failure.
* #567 steps (1)+(2) (commit ``6c54bba``) fixed both artifact generators in
  ``src/cyclerfinder/data/validation/v4_uranus_strict.py`` plus the
  hardcoded Umbriel/Oberon audit-field bug, and added a ``failure_mode``
  field (``FAILURE_MODE_CONVERGED`` / ``FAILURE_MODE_LAMBERT_NO_SOLUTION``
  / ``FAILURE_MODE_PLANET_CROSSING`` / ``FAILURE_MODE_INTEGRATOR_FAILURE``)
  to every per-cycle verdict so a genuine dynamical FAIL (planet-crossing)
  is never silently conflated with an unexplained solver FAIL.

This script is step (3): run the #338-style annual (100-epoch) sweep AND
the #559-style daily-DOY sweep (two year-long windows: 2000, 2030) under
the NOW-FIXED V4-strict driver, on all 6 candidates -- #312 itself plus the
5 #566 representatives -- and report each candidate's true PASS-band WIDTH
with an explicit ``failure_mode`` breakdown, not a bare pass/fail count.

Candidate provenance
---------------------
* #312 (Umbriel-Oberon-Umbriel): SILVER constants pulled verbatim from
  ``scripts/run_338_silver_v4strict_annual_sweep.py`` /
  ``scripts/run_559_silver_v4strict_daily_doy_scan.py`` (which themselves
  cite ``data/silver_327_verified.jsonl``, READ-ONLY here).
* The 5 #566 representatives (Titania-Oberon, Ariel-Umbriel, Ariel-Titania,
  Ariel-Oberon, Umbriel-Titania): parameter sets copied verbatim from
  ``scripts/run_566_gauntlet_five_representatives.py``'s ``CANDIDATES``
  tuple (itself independently verified against
  ``data/enumerate_563_symmetric_closures.jsonl`` -- not re-verified here,
  see that script's docstring for the line-by-line provenance).

Parallel substrate
-------------------
Reuses ``cyclerfinder.parallel.parallel_sweep`` (the joblib/loky substrate
first proven on this exact sweep shape by the #321 parallel demo,
``scripts/run_338_parallel_demo.py``: 100-epoch V4-strict annual sweep,
43.6s serial -> 8.6s parallel at 8 workers, 5.06x speedup, byte-for-byte
equivalence). ALL candidates x ALL sweep-type x ALL epoch cells are
assembled into ONE flat cell list and run through a SINGLE ``parallel_sweep``
call (one pool spin-up, not one per candidate).

Isolated-singleton guard
-------------------------
Per candidate, per sweep type, per CONTIGUOUS window (the daily sweep's two
year-long windows, 2000 and 2030, are never concatenated before the guard
runs -- see the CAUTION in
``src/cyclerfinder/data/sweep_diagnostics.py::detect_isolated_singleton_anomalies``,
which documents that concatenating non-contiguous windows manufactures a
spurious anomaly right at the seam), calls
:func:`cyclerfinder.data.sweep_diagnostics.detect_isolated_singleton_anomalies`
on the ``passes_v4_strict`` boolean series. With the #567 (1)+(2) fixes
landed, we expect FEWER (ideally zero) isolated flips than #559's raw
pre-fix 55/58 spike counts -- that would be a good sign (the two known
artifact generators are fixed), reported either way, not assumed.

Discipline anchors
-------------------
* READ-ONLY on ``src/cyclerfinder/data/validation/v[2-4]*.py`` -- no driver
  changes, only calls the existing (now-fixed) generic surfaces.
* NO catalogue writeback -- ``data/catalogue.yaml`` is untouched. Step (4)
  (the writeback-readiness verdict) is a SEPARATE, later, judgment pass per
  #567's own scope note.
* NO ``--no-verify``.

Run as::

    uv run python scripts/run_567_epoch_robustness_scan.py           # full (6 candidates)
    uv run python scripts/run_567_epoch_robustness_scan.py --smoke   # fast validation path
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from cyclerfinder.data.method_capability import MethodCapability  # noqa: E402
from cyclerfinder.data.preflight import preflight_search  # noqa: E402
from cyclerfinder.data.sweep_diagnostics import (  # noqa: E402
    detect_isolated_singleton_anomalies,
    singleton_anomaly_summary,
)
from cyclerfinder.data.validation.v2_moontour import run_v2_moontour  # noqa: E402
from cyclerfinder.data.validation.v3_3d import V3Verdict3D, run_v3_3d  # noqa: E402
from cyclerfinder.data.validation.v4_uranus import (  # noqa: E402
    URANIAN_PERTURBER_MOONS,
    URANUS_J2,
    URANUS_R_EQ_KM,
    V4_AGREEMENT_FLOOR_KMS,
    V4_N_CYCLES_MIN,
    V4UranusVerdict,
    run_v4_uranus,
)
from cyclerfinder.data.validation.v4_uranus_strict import (  # noqa: E402
    DEFAULT_LSK_PATH,
    DEFAULT_PCK_PATH,
    DEFAULT_URA_PATH,
    FAILURE_MODE_CONVERGED,
    FAILURE_MODE_INTEGRATOR_FAILURE,
    FAILURE_MODE_LAMBERT_NO_SOLUTION,
    FAILURE_MODE_PLANET_CROSSING,
    V4UranusStrictVerdict,
    run_v4_uranus_strict,
)
from cyclerfinder.parallel import ParallelSweepConfig, parallel_sweep  # noqa: E402

# --------------------------------------------------------------------------- #
# Candidate parameter sets
# --------------------------------------------------------------------------- #

# #558-established rotation-redundant constant, reused verbatim from
# #330/#335/#338/#559/#566 -- not re-derived here.
PHASE0_DEG = 29.999999999999996

# V4-strict requires n_cycles >= V4_N_CYCLES_MIN (== 3). Stay at the floor,
# matching #338/#559/#566's own choice -- runtime scales linearly with this.
N_CYCLES = V4_N_CYCLES_MIN


@dataclass(frozen=True)
class Candidate:
    """One of the 6 #567 epoch-robustness scan candidates."""

    candidate_id: str
    label: str
    sequence: tuple[str, ...]
    vinf_kms: tuple[float, ...]
    tof_days: tuple[float, ...]
    rel_offset_deg: float
    n_revs: tuple[int, ...]


CANDIDATES: tuple[Candidate, ...] = (
    # #312 itself -- constants pulled verbatim from
    # scripts/run_338_silver_v4strict_annual_sweep.py /
    # scripts/run_559_silver_v4strict_daily_doy_scan.py.
    Candidate(
        candidate_id="repeated-moon-uranus-00000041",
        label="#312 Umbriel-Oberon-Umbriel (SILVER, catalogued)",
        sequence=("Umbriel", "Oberon", "Umbriel"),
        vinf_kms=(
            0.9199258810725036,
            0.9604309791298091,
            0.8946936085078939,
        ),
        tof_days=(14.940560615336594, 14.940560615336594),
        rel_offset_deg=180.0,
        n_revs=(1, 1),
    ),
    # The 5 #566 representatives -- constants copied verbatim from
    # scripts/run_566_gauntlet_five_representatives.py's CANDIDATES tuple.
    Candidate(
        candidate_id="enum563-line57-titania-oberon-titania",
        label="Titania-Oberon-Titania (MANDATORY, #565 literature-clearance obligation)",
        sequence=("Titania", "Oberon", "Titania"),
        vinf_kms=(2.161767816675378, 1.9680495724521725, 2.1617678166753755),
        tof_days=(12.316046445872583, 12.316046445872583),
        rel_offset_deg=180.0,
        n_revs=(0, 0),
    ),
    Candidate(
        candidate_id="enum563-line2-ariel-umbriel-ariel",
        label="Ariel-Umbriel-Ariel",
        sequence=("Ariel", "Umbriel", "Ariel"),
        vinf_kms=(0.979040480994661, 1.3004234339628207, 0.9790404809946573),
        tof_days=(3.216088179066208, 3.216088179066208),
        rel_offset_deg=0.0,
        n_revs=(0, 0),
    ),
    Candidate(
        candidate_id="enum563-line12-ariel-titania-ariel",
        label="Ariel-Titania-Ariel",
        sequence=("Ariel", "Titania", "Ariel"),
        vinf_kms=(1.2306411593828481, 1.7185773183747601, 1.2306411593828457),
        tof_days=(5.320895317317783, 5.320895317317783),
        rel_offset_deg=0.0,
        n_revs=(0, 0),
    ),
    Candidate(
        candidate_id="enum563-line18-ariel-oberon-ariel",
        label="Ariel-Oberon-Ariel",
        sequence=("Ariel", "Oberon", "Ariel"),
        vinf_kms=(1.520866047614147, 1.8285940380726622, 1.5208660476141462),
        tof_days=(7.751820498940574, 7.751820498940574),
        rel_offset_deg=0.0,
        n_revs=(0, 0),
    ),
    Candidate(
        candidate_id="enum563-line26-umbriel-titania-umbriel",
        label="Umbriel-Titania-Umbriel",
        sequence=("Umbriel", "Titania", "Umbriel"),
        vinf_kms=(1.2295656768416439, 1.0058255988095806, 1.229565676841644),
        tof_days=(3.9544738760575804, 3.9544738760575804),
        rel_offset_deg=0.0,
        n_revs=(0, 0),
    ),
)

# --------------------------------------------------------------------------- #
# Epoch generation -- FULL (real scan) vs SMOKE (fast validation path)
# --------------------------------------------------------------------------- #

# Full annual sweep: one epoch per year, Y-06-21 -- matches #338 exactly.
FULL_ANNUAL_YEARS: tuple[int, ...] = tuple(range(2000, 2100))
ANNUAL_DOY_LABEL = "06-21T00:00:00"

# Full daily sweep: every calendar day across two year-long windows --
# matches #559 exactly (2000 brackets the known FAIL/PASS pair, 2030 is an
# independent decade-scale check).
FULL_DAILY_WINDOW_YEARS: tuple[int, ...] = (2000, 2030)

# Smoke path: 1 candidate (#312, the reference with known FAIL/PASS anchor
# points), a handful of epochs, no daily-window iteration -- validates the
# script runs end-to-end in well under a minute.
SMOKE_ANNUAL_YEARS: tuple[int, ...] = (2000, 2050, 2099)
# #338's own known FAIL (2000-01-15) / PASS (2000-06-21) anchor points, plus
# one more day -- a meaningful smoke check, not arbitrary padding.
SMOKE_DAILY_EPOCHS: tuple[str, ...] = (
    "2000-01-15T00:00:00",
    "2000-06-21T00:00:00",
    "2000-12-31T00:00:00",
)

OUT_JSONL = ROOT / "data" / "scan_567_epoch_robustness.jsonl"

_REGION_ID = "uranian-567-six-candidate-epoch-robustness-annual-daily-2026-07-11"
_METHOD = MethodCapability(
    genome=(
        "Annual (#338-style, 100 epoch/yr) + daily-DOY (#559-style, two "
        "year-long windows) launch-epoch robustness scan across all 6 "
        "#566/#312 candidates, under the #567(1)+(2)-fixed V4-strict "
        "driver -- read-only validation reporting, no genome/corrector change"
    ),
    corrector=(
        "existing V2->V3->V4-scipy->V4-strict chain "
        "(validation/v2_moontour.py, v3_3d.py, v4_uranus.py, "
        "v4_uranus_strict.py) POST #567 steps (1)+(2) bugfix (commit "
        "6c54bba): continuous Lambert branch selection by propagated "
        "terminal offset (kills the discontinuous departure-velocity-match "
        "tie-break), pre-screened + tagged FAILURE_MODE_PLANET_CROSSING "
        "(never silently excluded), and per-sequence (not hardcoded "
        "Umbriel/Oberon) audit e/i fields"
    ),
    capability_tags=frozenset(
        {
            "cr3bp",
            "real-ephemeris",
            "v4-strict",
            "uranian",
            "validation-reporting",
            "moontour",
            "epoch-robustness",
            "failure-mode-breakdown",
        }
    ),
    git_sha="working-tree",
)


def _annual_epochs(years: tuple[int, ...]) -> tuple[str, ...]:
    return tuple(f"{y:04d}-{ANNUAL_DOY_LABEL}" for y in years)


def _daily_epochs_for_window(year: int) -> tuple[str, ...]:
    start = datetime(year, 1, 1)
    end = datetime(year + 1, 1, 1)
    out: list[str] = []
    cur = start
    while cur < end:
        out.append(cur.strftime("%Y-%m-%dT00:00:00"))
        cur += timedelta(days=1)
    return tuple(out)


# --------------------------------------------------------------------------- #
# Parallel cell payload + closure (top-level for pickle safety)
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class _Cell:
    """One (candidate, sweep_type, window, epoch) cell.

    The V3 + V4-scipy verdicts are epoch-blind (circular-coplanar Kepler
    moons) so they are built ONCE per candidate outside the cell loop and
    passed into every cell for that candidate. Both are frozen dataclasses
    of floats/tuples, so they pickle cleanly under loky.
    """

    candidate_id: str
    label: str
    sequence: tuple[str, ...]
    vinf_kms: tuple[float, ...]
    tof_days: tuple[float, ...]
    rel_offset_deg: float
    n_revs: tuple[int, ...]
    sweep_type: str  # "annual" | "daily"
    window_label: str  # "2000-2099" for annual, "2000" / "2030" for daily
    epoch_utc: str
    v3: V3Verdict3D
    v4_scipy: V4UranusVerdict


def _aggregate_failure_mode(verdict: V4UranusStrictVerdict) -> str:
    """Roll a verdict's per-cycle ``failure_mode`` list up to one label.

    If every cycle's Lambert legs converged, the aggregate label is
    ``FAILURE_MODE_CONVERGED`` -- this is correct even when
    ``passes_v4_strict`` is False (e.g. every leg solved but the
    agreement-vs-V3 floor was blown, or bounded-drift failed): the FAIL is
    real dynamical drift, not a solver artifact, and the four-way taxonomy
    exists precisely to keep that distinction visible. Otherwise, returns
    the failure_mode of the FIRST non-converged cycle encountered (cycle
    order), which is the dominant/earliest cause of the epoch's FAIL.
    """
    for c in verdict.per_cycle:
        if c.failure_mode != FAILURE_MODE_CONVERGED:
            return c.failure_mode
    return FAILURE_MODE_CONVERGED


def _run_one_cell(cell: _Cell) -> dict[str, Any]:
    """Cell closure: run V4-strict at the cell's candidate/epoch, return a row dict."""
    t_run = time.time()
    try:
        v4s = run_v4_uranus_strict(
            cell.candidate_id,
            cell.sequence,
            cell.vinf_kms,
            cell.tof_days,
            cell.rel_offset_deg,
            cell.epoch_utc,
            None,
            v3_verdict=cell.v3,
            v4_scipy_verdict=cell.v4_scipy,
            n_cycles=N_CYCLES,
            n_revs=cell.n_revs,
            notes=(
                f"#567 epoch-robustness scan, candidate={cell.candidate_id}, "
                f"sweep_type={cell.sweep_type}, window={cell.window_label}, "
                f"epoch={cell.epoch_utc}, n_cycles={N_CYCLES}"
            ),
        )
        return {
            "kind": "epoch_robustness_row",
            "candidate_id": cell.candidate_id,
            "label": cell.label,
            "sweep_type": cell.sweep_type,
            "window_label": cell.window_label,
            "launch_epoch_utc": cell.epoch_utc,
            "passes_v4_strict": bool(v4s.passes_v4_strict),
            "bounded_drift_survives": bool(v4s.bounded_drift_survives),
            "epoch_failure_mode": _aggregate_failure_mode(v4s),
            "per_cycle_failure_modes": [c.failure_mode for c in v4s.per_cycle],
            "n_cycles_propagated": int(v4s.n_cycles_propagated),
            "n_cycles_requested": N_CYCLES,
            "drift_agreement_kms_vs_v3": float(v4s.drift_agreement_kms_vs_v3),
            "drift_agreement_kms_vs_v4_scipy": float(v4s.drift_agreement_kms_vs_v4_scipy),
            "audit_body1_name": v4s.audit_body1_name,
            "audit_body2_name": v4s.audit_body2_name,
            "eccentricity_used_e_body1": float(v4s.eccentricity_used_e_body1),
            "eccentricity_used_e_body2": float(v4s.eccentricity_used_e_body2),
            "inclination_used_deg_body1": float(v4s.inclination_used_deg_body1),
            "inclination_used_deg_body2": float(v4s.inclination_used_deg_body2),
            "wall_clock_s": float(time.time() - t_run),
        }
    except Exception as exc:
        return {
            "kind": "epoch_robustness_row",
            "candidate_id": cell.candidate_id,
            "label": cell.label,
            "sweep_type": cell.sweep_type,
            "window_label": cell.window_label,
            "launch_epoch_utc": cell.epoch_utc,
            "passes_v4_strict": False,
            "bounded_drift_survives": False,
            "epoch_failure_mode": FAILURE_MODE_INTEGRATOR_FAILURE,
            "per_cycle_failure_modes": [],
            "n_cycles_propagated": 0,
            "n_cycles_requested": N_CYCLES,
            "drift_agreement_kms_vs_v3": float("inf"),
            "drift_agreement_kms_vs_v4_scipy": float("inf"),
            "error": f"{type(exc).__name__}: {exc}",
            "wall_clock_s": float(time.time() - t_run),
        }


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT
        ).strip()
    except Exception:
        return "unknown"


def _build_v3_v4scipy(cand: Candidate, n_cycles: int) -> tuple[V3Verdict3D, V4UranusVerdict]:
    """V2 -> V3 -> V4-scipy chain (epoch-blind). Built once per candidate."""
    v2 = run_v2_moontour(
        cand.candidate_id,
        cand.sequence,
        cand.vinf_kms,
        cand.tof_days,
        cand.rel_offset_deg,
        None,
        n_cycles=n_cycles,
        n_revs=cand.n_revs,
        phase0_deg=PHASE0_DEG,
        notes=f"#567 epoch-robustness scan input chain, {cand.label}, n_cycles={n_cycles}",
    )
    v3 = run_v3_3d(
        cand.candidate_id,
        cand.sequence,
        cand.vinf_kms,
        cand.tof_days,
        cand.rel_offset_deg,
        None,
        v2_verdict=v2,
        n_cycles=n_cycles,
        n_revs=cand.n_revs,
        phase0_deg=PHASE0_DEG,
        notes=f"#567 epoch-robustness scan input chain, {cand.label}, n_cycles={n_cycles}",
    )
    v4 = run_v4_uranus(
        cand.candidate_id,
        cand.sequence,
        cand.vinf_kms,
        cand.tof_days,
        cand.rel_offset_deg,
        None,
        v3_verdict=v3,
        n_cycles=n_cycles,
        n_revs=cand.n_revs,
        phase0_deg=PHASE0_DEG,
        notes=f"#567 epoch-robustness scan input chain, {cand.label}, n_cycles={n_cycles}",
    )
    return v3, v4


@dataclass(frozen=True)
class _GroupSpec:
    """One contiguous (candidate, sweep_type, window) group of epochs."""

    candidate: Candidate
    sweep_type: str
    window_label: str
    epochs: tuple[str, ...]


def _build_group_specs(candidates: tuple[Candidate, ...], *, smoke: bool) -> list[_GroupSpec]:
    specs: list[_GroupSpec] = []
    for cand in candidates:
        if smoke:
            annual = _annual_epochs(SMOKE_ANNUAL_YEARS)
            specs.append(_GroupSpec(cand, "annual", "smoke-3yr", annual))
            specs.append(_GroupSpec(cand, "daily", "smoke-3day", SMOKE_DAILY_EPOCHS))
        else:
            annual = _annual_epochs(FULL_ANNUAL_YEARS)
            specs.append(
                _GroupSpec(
                    cand,
                    "annual",
                    f"{FULL_ANNUAL_YEARS[0]}-{FULL_ANNUAL_YEARS[-1]}",
                    annual,
                )
            )
            for year in FULL_DAILY_WINDOW_YEARS:
                specs.append(_GroupSpec(cand, "daily", str(year), _daily_epochs_for_window(year)))
    return specs


def _group_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    n = len(rows)
    n_pass = sum(1 for r in rows if r["passes_v4_strict"])
    n_fail = n - n_pass
    fail_mode_counts = {
        FAILURE_MODE_CONVERGED: 0,
        FAILURE_MODE_LAMBERT_NO_SOLUTION: 0,
        FAILURE_MODE_PLANET_CROSSING: 0,
        FAILURE_MODE_INTEGRATOR_FAILURE: 0,
    }
    for r in rows:
        if not r["passes_v4_strict"]:
            fail_mode_counts[r["epoch_failure_mode"]] = (
                fail_mode_counts.get(r["epoch_failure_mode"], 0) + 1
            )
    passes = [bool(r["passes_v4_strict"]) for r in rows]
    labels = [str(r["launch_epoch_utc"]) for r in rows]
    anomalies = detect_isolated_singleton_anomalies(passes, labels)
    return {
        "n_epochs": n,
        "n_pass": n_pass,
        "n_fail": n_fail,
        "pass_fraction": n_pass / n if n else 0.0,
        "fail_mode_counts": fail_mode_counts,
        "n_isolated_singleton_anomalies": len(anomalies),
        "isolated_singleton_anomaly_epochs": [a.label for a in anomalies],
        "singleton_anomaly_summary": singleton_anomaly_summary(passes, labels),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--smoke",
        action="store_true",
        help=(
            "fast validation path: 1 candidate (#312), 3 annual epochs + 3 "
            "daily epochs (6 total) -- confirms the script runs end-to-end, "
            "not a real scan"
        ),
    )
    ap.add_argument(
        "--n-workers",
        type=int,
        default=-1,
        help="parallel workers for parallel_sweep (-1 = all cores; default: -1)",
    )
    args = ap.parse_args()

    candidates = CANDIDATES[:1] if args.smoke else CANDIDATES
    group_specs = _build_group_specs(candidates, smoke=args.smoke)
    n_points = sum(len(g.epochs) for g in group_specs)

    preflight_search(
        task_no=567,
        region_id=_REGION_ID if not args.smoke else f"{_REGION_ID}-smoke",
        method=_METHOD,
        script_path=Path(__file__),
        n_points=n_points,
        override_reason=(
            "read-only validation-reporting epoch-robustness scan on 6 already "
            "gauntlet-tested (#566) / already-catalogued (#312) candidates, "
            "using the existing V2->V3->V4-scipy->V4-strict pipeline post the "
            "#567(1)+(2) artifact-generator bugfix -- not a discovery sweep. "
            "Per-point cost is directly measured by the #338 parallel demo "
            "(43.6s serial / 8.6s parallel for 100 epochs, i.e. ~0.09s/epoch "
            "at 8-way loky parallelism), so a multi-thousand-point total here "
            "is seconds-to-minutes scale, not the #520 unbudgeted-12-hour "
            "pattern this gate exists to catch."
        ),
    )

    sha = _git_sha()
    t0 = time.time()
    mode_label = "SMOKE" if args.smoke else "FULL"
    print(f"[#567] epoch-robustness scan ({mode_label}) -- sha={sha}", flush=True)
    print(f"[#567] candidates = {[c.candidate_id for c in candidates]}", flush=True)
    print(f"[#567] groups = {len(group_specs)}, total epoch-cells = {n_points}", flush=True)

    for p in (DEFAULT_LSK_PATH, DEFAULT_PCK_PATH, DEFAULT_URA_PATH):
        if not p.exists():
            print(f"[#567] FATAL: SPICE kernel missing: {p}", file=sys.stderr)
            return 1

    # Build the V2->V3->V4-scipy chain ONCE per candidate (epoch-blind).
    chains: dict[str, tuple[V3Verdict3D, V4UranusVerdict]] = {}
    print(
        f"[#567] building V2->V3->V4-scipy chains for {len(candidates)} candidate(s)...", flush=True
    )
    for cand in candidates:
        t_chain = time.time()
        chains[cand.candidate_id] = _build_v3_v4scipy(cand, N_CYCLES)
        print(
            f"[#567]   {cand.candidate_id}: chain ready (elapsed {time.time() - t_chain:.2f}s)",
            flush=True,
        )

    # Assemble the FLAT cell list (single pool, single parallel_sweep call)
    # while remembering each group's slice boundaries for post-processing.
    cells: list[_Cell] = []
    group_offsets: list[tuple[_GroupSpec, int, int]] = []  # (spec, start, end)
    for spec in group_specs:
        v3, v4_scipy = chains[spec.candidate.candidate_id]
        start = len(cells)
        for epoch in spec.epochs:
            cells.append(
                _Cell(
                    candidate_id=spec.candidate.candidate_id,
                    label=spec.candidate.label,
                    sequence=spec.candidate.sequence,
                    vinf_kms=spec.candidate.vinf_kms,
                    tof_days=spec.candidate.tof_days,
                    rel_offset_deg=spec.candidate.rel_offset_deg,
                    n_revs=spec.candidate.n_revs,
                    sweep_type=spec.sweep_type,
                    window_label=spec.window_label,
                    epoch_utc=epoch,
                    v3=v3,
                    v4_scipy=v4_scipy,
                )
            )
        group_offsets.append((spec, start, len(cells)))

    print(
        f"[#567] running {len(cells)} cells through parallel_sweep (n_workers={args.n_workers})...",
        flush=True,
    )
    cfg = ParallelSweepConfig(n_workers=args.n_workers, backend="loky", verbose=0)
    t_par = time.time()
    result = parallel_sweep(cells, _run_one_cell, config=cfg)
    par_elapsed = time.time() - t_par
    print(
        f"[#567]   parallel_sweep done: {par_elapsed:.1f}s, "
        f"n_succeeded={result.n_succeeded}, n_failed={result.n_failed}",
        flush=True,
    )
    if result.n_failed:
        print(f"[#567]   NOTE: {result.notes}", flush=True)

    all_rows: list[dict[str, Any]] = [
        r
        if r is not None
        else {
            "kind": "epoch_robustness_row",
            "error": "cell returned None (joblib-level failure)",
            "passes_v4_strict": False,
            "epoch_failure_mode": FAILURE_MODE_INTEGRATOR_FAILURE,
        }
        for r in result.results
    ]

    # Per-group summaries.
    group_summaries: list[dict[str, Any]] = []
    for spec, start, end in group_offsets:
        rows = all_rows[start:end]
        summary = _group_summary(rows)
        summary.update(
            {
                "candidate_id": spec.candidate.candidate_id,
                "label": spec.candidate.label,
                "sweep_type": spec.sweep_type,
                "window_label": spec.window_label,
            }
        )
        group_summaries.append(summary)
        print(
            f"[#567] {spec.candidate.candidate_id} / {spec.sweep_type}/{spec.window_label}: "
            f"PASS={summary['n_pass']}/{summary['n_epochs']} "
            f"({summary['pass_fraction']:.1%}) | fail_modes={summary['fail_mode_counts']} | "
            f"{summary['singleton_anomaly_summary']}",
            flush=True,
        )

    # Per-candidate roll-up across its sweep-type groups.
    candidate_rollups: list[dict[str, Any]] = []
    for cand in candidates:
        cand_groups = [g for g in group_summaries if g["candidate_id"] == cand.candidate_id]
        n_epochs_total = sum(g["n_epochs"] for g in cand_groups)
        n_pass_total = sum(g["n_pass"] for g in cand_groups)
        n_anomalies_total = sum(g["n_isolated_singleton_anomalies"] for g in cand_groups)
        candidate_rollups.append(
            {
                "candidate_id": cand.candidate_id,
                "label": cand.label,
                "n_epochs_total": n_epochs_total,
                "n_pass_total": n_pass_total,
                "pass_fraction_total": (n_pass_total / n_epochs_total if n_epochs_total else 0.0),
                "n_isolated_singleton_anomalies_total": n_anomalies_total,
                "groups": [
                    {
                        "sweep_type": g["sweep_type"],
                        "window_label": g["window_label"],
                        "n_epochs": g["n_epochs"],
                        "n_pass": g["n_pass"],
                        "pass_fraction": g["pass_fraction"],
                        "fail_mode_counts": g["fail_mode_counts"],
                        "n_isolated_singleton_anomalies": g["n_isolated_singleton_anomalies"],
                    }
                    for g in cand_groups
                ],
            }
        )

    # --------------------------------------------------------------------- #
    # Write output
    # --------------------------------------------------------------------- #
    out_rows: list[dict[str, Any]] = []
    out_rows.append(
        {
            "_meta": True,
            "task": "#567 step (3) generalized epoch-robustness scan across 6 candidates",
            "mode": mode_label,
            "candidates": [
                {
                    "candidate_id": c.candidate_id,
                    "label": c.label,
                    "sequence": list(c.sequence),
                    "vinf_kms": list(c.vinf_kms),
                    "tof_days": list(c.tof_days),
                    "rel_offset_deg": c.rel_offset_deg,
                    "n_revs": list(c.n_revs),
                }
                for c in candidates
            ],
            "phase0_deg": PHASE0_DEG,
            "n_cycles": N_CYCLES,
            "driver_floors": {
                "v4_agreement_floor_kms": V4_AGREEMENT_FLOOR_KMS,
                "v4_n_cycles_min": V4_N_CYCLES_MIN,
            },
            "v4_strict_model": {
                "integrator": "scipy DOP853",
                "uranus_j2": URANUS_J2,
                "uranus_r_eq_km": URANUS_R_EQ_KM,
                "perturber_moons": list(URANIAN_PERTURBER_MOONS),
                "spice_kernels_used": [
                    str(DEFAULT_LSK_PATH),
                    str(DEFAULT_PCK_PATH),
                    str(DEFAULT_URA_PATH),
                ],
                "post_567_bugfix": (
                    "commit 6c54bba (continuity + planet-crossing tag + audit fields)"
                ),
            },
            "n_epoch_cells_total": len(cells),
            "n_groups": len(group_specs),
            "parallel_config": {"n_workers": args.n_workers, "backend": "loky"},
            "parallel_wall_seconds": par_elapsed,
            "git_sha": sha,
        }
    )
    out_rows.extend(all_rows)
    out_rows.extend({"kind": "group_summary", **g} for g in group_summaries)
    out_rows.append(
        {
            "_meta": True,
            "kind": "headline",
            "mode": mode_label,
            "candidate_rollups": candidate_rollups,
            "writeback_to_catalogue": False,
            "next_step": (
                "#567 step (4): writeback-readiness verdict per candidate "
                "(band width + whether it clears at #312's V4 level or a "
                "capped level) -- a SEPARATE judgment pass, not part of this script"
            ),
            "elapsed_s": time.time() - t0,
        }
    )

    OUT_JSONL.parent.mkdir(parents=True, exist_ok=True)
    with OUT_JSONL.open("w", encoding="utf-8") as fh:
        for row in out_rows:
            fh.write(json.dumps(row) + "\n")

    print(f"\n[#567] wrote {OUT_JSONL}", flush=True)
    print(f"[#567] total elapsed {time.time() - t0:.1f}s", flush=True)
    for r in candidate_rollups:
        print(
            f"[#567] SUMMARY {r['candidate_id']}: "
            f"{r['n_pass_total']}/{r['n_epochs_total']} PASS "
            f"({r['pass_fraction_total']:.1%}), "
            f"{r['n_isolated_singleton_anomalies_total']} isolated singleton anomalies",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
