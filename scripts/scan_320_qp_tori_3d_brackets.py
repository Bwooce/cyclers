"""#320 Vector A -- QP-tori at the 12 #299 Neimark-Sacker brackets.

Phase 1 of #320 (first systematic quasi_cycler discovery sweep). For each of
the 12 Neimark-Sacker brackets surfaced by #299's 3D bifurcation tracker over
the #296 Earth-Moon family, this script:

  1. Pulls the bracket's parent state + period from
     ``data/family_296_3d_em_11.jsonl`` at ``step_a`` (the closer-to-crossing
     family member);
  2. Runs the Olikara-Howell GMOS corrector (:func:`correct_qp_torus`) to seek
     a converged 2-torus born at the bifurcation point;
  3. If the torus converged, runs the V1_qp + V2_qp gauntlet
     (:func:`run_v1_qp` / :func:`run_v2_qp`);
  4. Records the frequency ratio, irrationality check
     (:func:`is_practically_irrational`), and a JSONL row.

The output ``data/scan_320_qp_tori_3d_brackets.jsonl`` carries:
  * one leading ``_meta`` row
  * one row per bracket
  * one trailing ``summary`` row

NO catalogue writeback. SILVER survivors (V1_qp PASS + V2_qp PASS + irrational
frequency ratio) feed a #306-style follow-up task; this Phase 1 outputs JSONL
+ doc only.

Discipline:
  * READ-ONLY on substrate modules (``qp_tori.py``, ``v1_qp.py``, ``v2_qp.py``).
  * Sourced golden discipline: parent states + Floquet pairs come from #299's
    bracket inventory (commit ``c83b6f9`` / family ``data/family_296_3d_em_11
    .jsonl`` with seed Antoniadou-Voyatzis 2018, see header ``seed`` field).

Run as::

    uv run python scripts/scan_320_qp_tori_3d_brackets.py
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from cyclerfinder.core import cr3bp  # noqa: E402
from cyclerfinder.data.validation.v1_qp import run_v1_qp  # noqa: E402
from cyclerfinder.data.validation.v2_qp import run_v2_qp  # noqa: E402
from cyclerfinder.genome.qp_tori import (  # noqa: E402
    correct_qp_torus,
    is_practically_irrational,
)

SUBFAMILIES_FILE = ROOT / "data" / "family_296_3d_subfamilies_299.jsonl"
PARENT_FAMILY_FILE = ROOT / "data" / "family_296_3d_em_11.jsonl"
OUT_PATH = ROOT / "data" / "scan_320_qp_tori_3d_brackets.jsonl"


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT
        ).strip()
    except Exception:
        return "unknown"


def _em_system() -> cr3bp.CR3BPSystem:
    """Earth-Moon CR3BP system. Same factory used by the qp_tori tests and the
    #296/#299 family files."""
    return cr3bp.cr3bp_system("Earth", "Moon")


def _load_bracket_inventory() -> list[dict[str, Any]]:
    """Return the ``bracket_inventory`` block from the #299 subfamilies header
    (12 Neimark-Sacker brackets surfaced by the 3D bifurcation tracker)."""
    if not SUBFAMILIES_FILE.exists():
        raise FileNotFoundError(f"missing {SUBFAMILIES_FILE} -- #299 must have run")
    with SUBFAMILIES_FILE.open() as f:
        for line in f:
            obj = json.loads(line)
            if obj.get("type") == "header":
                inv = obj.get("bracket_inventory") or []
                return list(inv)
    raise RuntimeError(f"no header row in {SUBFAMILIES_FILE}")


def _load_parent_member(step_index: int) -> dict[str, Any]:
    """Return the #296 3D Earth-Moon family member at ``step_index`` (the bracket's
    closer-to-crossing family member -- ``step_a`` in the bracket spec)."""
    if not PARENT_FAMILY_FILE.exists():
        raise FileNotFoundError(f"missing {PARENT_FAMILY_FILE}")
    with PARENT_FAMILY_FILE.open() as f:
        for line in f:
            obj = json.loads(line)
            if obj.get("type") == "header":
                continue
            if obj.get("step_index") == step_index:
                return dict(obj)
    raise RuntimeError(f"step_index={step_index} not found in {PARENT_FAMILY_FILE}")


def _per_bracket(
    bracket_idx: int,
    bracket: dict[str, Any],
    parent: dict[str, Any],
    system: cr3bp.CR3BPSystem,
    *,
    n_trans: int = 2,
    initial_torus_amplitude: float = 5e-4,
    tol: float = 1e-8,
    max_iter: int = 40,
    independent_tol: float = 1e-3,
) -> dict[str, Any]:
    """Compute one bracket: corrector + V1_qp + V2_qp + irrationality verdict.

    The (n_trans=2, amplitude=5e-4) defaults match the smoke test's calibration
    (see ``tests/genome/test_qp_tori.py::test_sourced_neimark_sacker_smoke``).
    """
    parent_state = np.asarray(parent["state_nd"], dtype=np.float64)
    parent_period = float(parent["T_TU"])
    k = int(bracket["k"])
    lam_a = complex(bracket["eig_a_re"], bracket["eig_a_im"])
    lam_b = complex(bracket["eig_b_re"], bracket["eig_b_im"])

    row: dict[str, Any] = {
        "bracket_idx": bracket_idx,
        "k": k,
        "step_a": int(bracket["step_a"]),
        "step_b": int(bracket["step_b"]),
        "T_a_TU": float(bracket["T_a"]),
        "C_a": float(bracket["C_a"]),
        "z0_a": float(bracket["z0_a"]),
        "d_a": float(bracket["d_a"]),
        "lam_a": {"re": float(np.real(lam_a)), "im": float(np.imag(lam_a))},
        "lam_b": {"re": float(np.real(lam_b)), "im": float(np.imag(lam_b))},
        "parent_state_nd": list(parent_state),
        "parent_T_TU": parent_period,
        "n_trans": n_trans,
        "initial_torus_amplitude": initial_torus_amplitude,
    }

    t0 = time.time()
    try:
        torus = correct_qp_torus(
            system,
            parent_state,
            parent_period,
            (lam_a, lam_b),
            k=k,
            n_long=16,
            n_trans=n_trans,
            initial_torus_amplitude=initial_torus_amplitude,
            tol=tol,
            max_iter=max_iter,
            independent_tol=independent_tol,
            notes=f"scan_320_vectorA_bracket_{bracket_idx}",
        )
    except Exception as e:
        row.update(
            {
                "corrector_status": "error",
                "corrector_error": repr(e),
                "elapsed_s": time.time() - t0,
            }
        )
        return row

    row.update(
        {
            "corrector_status": "ok",
            "invariance_residual": float(torus.invariance_residual),
            "independent_closure_residual": float(torus.independent_closure_residual),
            "corrector_converged": bool(torus.converged),
            "omega_long": float(torus.omega_long),
            "omega_trans": float(torus.omega_trans),
            "rho": float(torus.rho),
            "t_strob": float(torus.t_strob),
            "n_iter": int(torus.n_iter),
            "freq_ratio": float(torus.omega_trans / torus.omega_long)
            if torus.omega_long != 0.0
            else float("nan"),
            "freq_ratio_minus_inv_k": float((torus.omega_trans / torus.omega_long) - 1.0 / k)
            if torus.omega_long != 0.0
            else float("nan"),
            "n_modes": int(torus.n_modes),
        }
    )

    # Irrationality of the frequency ratio (a phase-locked rational ratio means
    # the "torus" collapsed to the parent's periodic orbit / one of its k-fold
    # covers).
    ratio = torus.omega_trans / torus.omega_long
    row["irrational_freq_ratio"] = bool(
        is_practically_irrational(ratio, max_denominator=10, tol=1e-3)
    )

    # Only proceed to V1/V2 gauntlets if the corrector converged (otherwise the
    # gauntlet just re-reports an infinity).
    if not math.isfinite(torus.invariance_residual) or torus.invariance_residual > 1e-4:
        row["v1_qp_status"] = "skipped_no_convergence"
        row["v2_qp_status"] = "skipped_no_convergence"
        row["silver_status"] = "FAIL_corrector"
        row["elapsed_s"] = time.time() - t0
        return row

    try:
        v1 = run_v1_qp(
            candidate_id=f"vectorA_bracket_{bracket_idx}",
            torus=torus,
            notes="scan_320_vectorA",
        )
    except Exception as e:
        row["v1_qp_status"] = "error"
        row["v1_qp_error"] = repr(e)
        row["silver_status"] = "FAIL_v1_error"
        row["elapsed_s"] = time.time() - t0
        return row
    row.update(
        {
            "v1_qp_status": "ok",
            "v1_fourier_norm": float(v1.invariance_residual_fourier_norm),
            "v1_independent_nondim": float(v1.independent_invariance_residual_nondim),
            "v1_independent_km": float(v1.independent_residual_km),
            "v1_passes": bool(v1.passes_v1_qp),
        }
    )

    try:
        v2 = run_v2_qp(
            candidate_id=f"vectorA_bracket_{bracket_idx}",
            torus=torus,
            notes="scan_320_vectorA",
        )
    except Exception as e:
        row["v2_qp_status"] = "error"
        row["v2_qp_error"] = repr(e)
        row["silver_status"] = "FAIL_v2_error"
        row["elapsed_s"] = time.time() - t0
        return row
    row.update(
        {
            "v2_qp_status": "ok",
            "v2_max_drift_nondim": float(v2.max_invariance_drift),
            "v2_max_drift_km": float(v2.max_invariance_drift_km),
            "v2_per_cycle": list(v2.per_cycle_invariance_residual),
            "v2_passes": bool(v2.passes_v2_qp),
        }
    )

    silver = bool(v1.passes_v1_qp and v2.passes_v2_qp and row["irrational_freq_ratio"])
    row["silver_status"] = "SILVER" if silver else "BRONZE_or_lower"

    row["elapsed_s"] = time.time() - t0
    return row


def main() -> int:
    sha = _git_sha()
    print(f"[320-A] QP-tori at #299 3D brackets -- sha={sha}", flush=True)
    print(f"[320-A] out={OUT_PATH}", flush=True)
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    inventory = _load_bracket_inventory()
    n_brackets = len(inventory)
    print(f"[320-A] bracket inventory: {n_brackets} Neimark-Sacker brackets", flush=True)

    system = _em_system()
    rows: list[dict[str, Any]] = []
    t0 = time.time()

    with OUT_PATH.open("w", encoding="utf-8") as fh:
        meta = {
            "_meta": True,
            "task": "#320 Vector A -- QP-tori at #299 3D Earth-Moon Neimark-Sacker brackets",
            "primary": "Earth",
            "secondary": "Moon",
            "system_mu_source": "cr3bp_system('Earth', 'Moon') -- JPL SSD gm_de440",
            "subfamilies_file": str(SUBFAMILIES_FILE.relative_to(ROOT)),
            "parent_family_file": str(PARENT_FAMILY_FILE.relative_to(ROOT)),
            "n_brackets": n_brackets,
            "corrector_defaults": {
                "n_long": 16,
                "n_trans": 2,
                "initial_torus_amplitude": 5e-4,
                "tol": 1e-8,
                "max_iter": 40,
                "independent_tol": 1e-3,
            },
            "v1_qp_floors": {
                "fourier_norm": 1.0e-5,
                "independent_nondim": 1.0e-4,
            },
            "v2_qp_floors": {
                "drift_nondim": 5.0e-2,
                "n_cycles": 3,
            },
            "git_sha": sha,
        }
        fh.write(json.dumps(meta) + "\n")
        fh.flush()

        for i, br in enumerate(inventory):
            print(
                f"[320-A] bracket {i}/{n_brackets} k={br['k']} step_a={br['step_a']}"
                f" z0_a={br['z0_a']:.5f} T_a={br['T_a']:.5f}",
                flush=True,
            )
            parent = _load_parent_member(int(br["step_a"]))
            row = _per_bracket(i, br, parent, system)
            rows.append(row)
            fh.write(json.dumps(row) + "\n")
            fh.flush()
            print(
                f"[320-A]   -> {row.get('silver_status', '?')}"
                f" inv={row.get('invariance_residual', float('nan')):.3e}"
                f" v1_passes={row.get('v1_passes')}"
                f" v2_passes={row.get('v2_passes')}"
                f" irrational={row.get('irrational_freq_ratio')}"
                f" elapsed={row['elapsed_s']:.1f}s",
                flush=True,
            )

        # Summary.
        n_converged = sum(1 for r in rows if r.get("corrector_status") == "ok")
        n_corrector_converged_strict = sum(1 for r in rows if r.get("corrector_converged"))
        n_v1_pass = sum(1 for r in rows if r.get("v1_passes"))
        n_v2_pass = sum(1 for r in rows if r.get("v2_passes"))
        n_irrational = sum(1 for r in rows if r.get("irrational_freq_ratio"))
        n_silver = sum(1 for r in rows if r.get("silver_status") == "SILVER")
        elapsed = time.time() - t0
        summary = {
            "_meta": True,
            "kind": "summary",
            "n_brackets": n_brackets,
            "n_corrector_ran": n_converged,
            "n_corrector_strict_converged": n_corrector_converged_strict,
            "n_v1_qp_pass": n_v1_pass,
            "n_v2_qp_pass": n_v2_pass,
            "n_irrational_freq_ratio": n_irrational,
            "n_silver": n_silver,
            "elapsed_s": elapsed,
            "git_sha": sha,
            "discipline": (
                "SILVER status here is V1_qp PASS + V2_qp PASS + irrational "
                "frequency ratio. NOT a catalogue admission; novelty check "
                "and #306-style follow-up are required."
            ),
        }
        fh.write(json.dumps(summary) + "\n")
        fh.flush()

    print(
        f"[320-A] DONE -- {n_silver}/{n_brackets} SILVER ({n_v1_pass} V1_pass, "
        f"{n_v2_pass} V2_pass, {n_irrational} irrational) in {elapsed:.1f}s",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
