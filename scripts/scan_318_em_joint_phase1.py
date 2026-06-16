"""#318 Phase 1 - multi-axis joint-search EM probe.

A small Phase-1 probe that walks the four cycler-discovery axes (powered
maintenance, multi-rev Lambert, 3D / broken-plane, epoch-locked validity
window) jointly over the Earth-Mars Aldrin tour. The point is not to
discover anything -- it is to demonstrate the *capability* substrate at
:mod:`cyclerfinder.search.multi_axis_search` and surface the cell-by-cell
candidate count + best survivor for each axis corner, so a Phase 2 follow-up
knows which cells warrant smarter sampling.

The probe runs joint_axis_search at the Aldrin canonical tour with:

  * Axis A (powered budgets): (0.0, 0.05, 0.1, 0.2, 0.5, 1.5) km/s.
    Records the requested budget on every cell; the powered driver still
    computes the optimum, the cell tag is the joint-axis index.
  * Axis B (n_revs per leg): ((0, 1), (0, 1)).
    Direct + 1-rev branch on each leg; the four-cell Cartesian product
    ((0,0), (0,1), (1,0), (1,1)) is the multi-rev sweep.
  * Axis C (z0 amplitudes nondim): (0.0, 1e-3, 1e-2, 5e-2).
    Planar baseline + three 3D requests; Phase 1 records the amplitude,
    Phase 2 will drive the #291 corrector inline.
  * Axis D (launch epochs): (None, '2030-01-01').
    Epoch-blind baseline + one real-DE440 cell (will fall back to the
    powered ephemeris if astropy isn't available; the cell tag is
    recorded regardless).

Total cells: 6 powered x 4 revs x 4 z0 x 2 epoch = 192. Phase 1's
Cartesian-product is the *capability* sweep; Phase 2 will adopt smarter
sampling (Sobol, surrogate-driven).

Expected outcome: most cells rediscover Aldrin (the joint-zero corner is
the only one with a definite-source target); some cells will fail to close
(the optimiser may not find the multi-rev branch, the 3D request is
record-only so it never actually closes a 3D orbit). The interesting cells
are at the JOINT corner -- (powered + multi-rev + 3D + epoch-locked) --
where individual axes don't add but combined they might.

NO catalogue writeback. NO novelty claims.

Run as::

    uv run python scripts/scan_318_em_joint_phase1.py

Outputs ``data/scan_318_em_joint_phase1.jsonl`` -- one leading ``_meta``
row, then per-candidate rows, then a trailing summary row.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from cyclerfinder.core.ephemeris import Ephemeris  # noqa: E402
from cyclerfinder.search.maintain import _default_t0_guess  # noqa: E402
from cyclerfinder.search.multi_axis_search import (  # noqa: E402
    JointAxisCandidate,
    joint_axis_search,
)


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT
        ).strip()
    except Exception:
        return "unknown"


# ----------------------------------------------------------------------------
# Aldrin reference parameters
# ----------------------------------------------------------------------------

ALDRIN_EM_TOF_DAYS = 146.0
ALDRIN_ME_TOF_DAYS = 634.0
ALDRIN_EARTH_FLYBY_ALT_KM = 200.0
ALDRIN_TOF_JITTER = (20.0, 60.0)


# ----------------------------------------------------------------------------
# Phase 1 joint-axis grids
# ----------------------------------------------------------------------------

# Axis A: powered budgets (km/s). Records the cell's requested budget;
# Phase 1 does NOT gate on it.
POWERED_BUDGETS_KMS: tuple[float, ...] = (0.0, 0.05, 0.1, 0.2, 0.5, 1.5)

# Axis B: per-leg revs grid -- direct + 1-rev branch on each leg. Phase 1
# Cartesian product expands to ((0,0), (0,1), (1,0), (1,1)).
N_REVS_GRID_PER_LEG: tuple[tuple[int, ...], ...] = ((0, 1), (0, 1))

# Axis C: z0 amplitudes (non-dim). Planar baseline + three 3D requests.
Z0_AMPLITUDES_NONDIM: tuple[float, ...] = (0.0, 1e-3, 1e-2, 5e-2)

# Axis D: launch epochs. None = epoch-blind (analytic-circular backend);
# a UTC string activates the real-ephemeris closure witness.
LAUNCH_EPOCH_GRID: tuple[str | None, ...] = (None, "2030-01-01")


def _ephem_circular() -> Ephemeris:
    return Ephemeris("circular")


def _try_ephem_astropy() -> Ephemeris | None:
    """Best-effort astropy backend init; falls back to None on environment lacks."""
    try:
        return Ephemeris("astropy")
    except Exception:
        return None


def main() -> None:
    out_path = ROOT / "data" / "scan_318_em_joint_phase1.jsonl"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    eph = _ephem_circular()
    epoch_ephem = _try_ephem_astropy() or eph

    # Meta header row.
    meta = {
        "_meta": {
            "task": "#318 Phase 1 - multi-axis joint-search EM probe",
            "git_sha": _git_sha(),
            "driver": "cyclerfinder.search.multi_axis_search.joint_axis_search",
            "primary_sequence": ["E", "M", "E"],
            "k_synodic": 1,
            "axis_a_powered_budgets_kms": list(POWERED_BUDGETS_KMS),
            "axis_b_n_revs_grid_per_leg": [list(g) for g in N_REVS_GRID_PER_LEG],
            "axis_c_z0_amplitudes_nondim": list(Z0_AMPLITUDES_NONDIM),
            "axis_d_launch_epoch_grid": list(LAUNCH_EPOCH_GRID),
            "cells_total": (
                len(POWERED_BUDGETS_KMS)
                * len(N_REVS_GRID_PER_LEG[0])
                * len(N_REVS_GRID_PER_LEG[1])
                * len(Z0_AMPLITUDES_NONDIM)
                * len(LAUNCH_EPOCH_GRID)
            ),
            "epoch_locked_ephem_kind": (
                "astropy" if epoch_ephem is not eph else "circular_fallback"
            ),
            "discipline": [
                "Phase 1 substrate only",
                "No catalogue writeback",
                "No novelty claims",
                "The four axis modules are unmodified (composition only)",
                "Sourced golden anchor: aldrin-classic-em-k1-outbound "
                "at (powered=0, n_revs=0, z0=0, epoch=None)",
            ],
        },
    }

    candidates: list[JointAxisCandidate] = []
    with out_path.open("w") as f:
        f.write(json.dumps(meta) + "\n")
        for cand in joint_axis_search(
            primary_sequence=("E", "M", "E"),
            k_synodic=1,
            ephem=eph,
            powered_budgets_kms=POWERED_BUDGETS_KMS,
            n_revs_grid_per_leg=N_REVS_GRID_PER_LEG,
            z0_amplitudes_nondim=Z0_AMPLITUDES_NONDIM,
            launch_epoch_grid=LAUNCH_EPOCH_GRID,
            leg_tof_guesses_days=(ALDRIN_EM_TOF_DAYS, ALDRIN_ME_TOF_DAYS),
            leg_tof_bounds_days=((100.0, 250.0), (400.0, 900.0)),
            tof_jitter_half_days=ALDRIN_TOF_JITTER,
            synodic_pair=("E", "M"),
            closure_body="E",
            closure_flyby_alt_km=ALDRIN_EARTH_FLYBY_ALT_KM,
            n_starts=4,
            seed=0,
            epoch_locked_ephem=epoch_ephem,
            t0_guess_sec_for_epoch_blind=_default_t0_guess(ALDRIN_EM_TOF_DAYS),
        ):
            f.write(json.dumps(cand.as_dict()) + "\n")
            candidates.append(cand)
            print(
                f"cell={{powered={cand.powered_budget_kms_requested}, "
                f"n_revs={cand.n_revs_per_leg}, z0={cand.z0_amplitude_nondim:.0e}, "
                f"epoch={cand.launch_epoch_utc}}}  "
                f"-> dv={cand.powered_maintenance_dv_kms_per_synodic:.4f} km/s, "
                f"vinf={tuple(round(v, 2) for v in cand.vinf_tuple_kms)}",
                flush=True,
            )

        # Per-axis-corner summary.
        summary: dict = {
            "_summary": {
                "candidates_total": len(candidates),
                "cells_total": meta["_meta"]["cells_total"],
                "survivors_by_axis_a_powered": {},
                "survivors_by_axis_b_n_revs": {},
                "survivors_by_axis_c_z0": {},
                "survivors_by_axis_d_epoch": {},
                "best_per_corner": {},
            },
        }
        for c in candidates:
            ka = str(c.powered_budget_kms_requested)
            kb = str(list(c.n_revs_per_leg))
            kc = f"{c.z0_amplitude_nondim:.3e}"
            kd = str(c.launch_epoch_utc)
            for key, bucket in [
                (ka, "survivors_by_axis_a_powered"),
                (kb, "survivors_by_axis_b_n_revs"),
                (kc, "survivors_by_axis_c_z0"),
                (kd, "survivors_by_axis_d_epoch"),
            ]:
                summary["_summary"][bucket].setdefault(key, 0)
                summary["_summary"][bucket][key] += 1

        # Best per-axis-corner: lowest closure_residual on each axis bucket.
        for axis_label, key_fn in [
            ("axis_a_powered", lambda c: str(c.powered_budget_kms_requested)),
            ("axis_b_n_revs", lambda c: str(list(c.n_revs_per_leg))),
            ("axis_c_z0", lambda c: f"{c.z0_amplitude_nondim:.3e}"),
            ("axis_d_epoch", lambda c: str(c.launch_epoch_utc)),
        ]:
            best: dict = {}
            for c in candidates:
                k = key_fn(c)
                if k not in best or c.closure_residual_kms < best[k]["closure_residual_kms"]:
                    best[k] = {
                        "closure_residual_kms": c.closure_residual_kms,
                        "n_revs_per_leg": list(c.n_revs_per_leg),
                        "powered_budget_kms_requested": c.powered_budget_kms_requested,
                        "z0_amplitude_nondim": c.z0_amplitude_nondim,
                        "launch_epoch_utc": c.launch_epoch_utc,
                        "vinf_tuple_kms": list(c.vinf_tuple_kms),
                    }
            summary["_summary"]["best_per_corner"][axis_label] = best

        f.write(json.dumps(summary) + "\n")

    print()
    print(f"wrote {len(candidates)} candidate rows to {out_path}")
    print(f"cells_total = {meta['_meta']['cells_total']}")
    print(f"survivors = {len(candidates)}")


if __name__ == "__main__":
    main()
