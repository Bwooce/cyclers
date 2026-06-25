"""#305 BCR4BP V0-V5 gauntlet probe — characterize #303/#304 families + POL1.

Runs the built BCR4BP V0->V3 gauntlet (src/cyclerfinder/data/validation/
v{0,1,2,3}_bcr4bp.py) on:

  * the sourced Andreu / Rosales-Jorba POL1 L1 substitute (Sun-commensurate,
    n=1) — the golden the gauntlet was validated against;
  * the #303 planar L1 Lyapunov mu_sun-continuation family
    (data/bcr4bp_l1_family_303.jsonl);
  * the #304 L1 southern halo mu_sun-continuation family
    (data/bcr4bp_halo_family_304.jsonl).

REPORTING ONLY. This script does NOT write the catalogue and does NOT
self-admit any row. A V0-V3-passing family member is a KNOWN-REPRODUCTION
candidate (first_published=Andreu, our_status=known-reproduction) to be
FLAGGED for human review with the proposed row + attribution — never
self-admitted (per the #305 task brief and the orbit-closure discipline).

Run: ``uv run python scripts/run_305_bcr4bp_gauntlet_probe.py``
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import cyclerfinder.core.bcr4bp as bcr4bp
from cyclerfinder.data.validation.v0_bcr4bp import run_v0_bcr4bp
from cyclerfinder.data.validation.v1_bcr4bp import run_v1_bcr4bp
from cyclerfinder.data.validation.v2_bcr4bp import run_v2_bcr4bp
from cyclerfinder.data.validation.v3_bcr4bp import run_v3_bcr4bp
from cyclerfinder.genome.bcr4bp_genome import (
    IDX_XDOT,
    IDX_Y,
    IDX_YDOT,
    IDX_ZDOT,
    BCR4BPPeriodicOrbit,
    correct_bcr4bp_periodic,
)

_DATA = Path(__file__).resolve().parent.parent / "data"
_POL1_X = -0.8369141677649317
_POL1_PY = -0.8391311559808445
_POL1_SEED = np.array([_POL1_X, 0.0, 0.0, 0.0, _POL1_PY - _POL1_X, 0.0], dtype=np.float64)


def _close_pol1() -> BCR4BPPeriodicOrbit:
    sys_bcr = bcr4bp.andreu_default()
    period_fixed = bcr4bp.sun_commensurate_period(sys_bcr.omega_sun_nondim, n=1)
    return correct_bcr4bp_periodic(
        sys_bcr,
        _POL1_SEED,
        period_fixed,
        sun_commensurate_n=1,
        free_vars=(0, IDX_YDOT),
        residual_indices=(IDX_Y, IDX_XDOT, IDX_ZDOT),
        is_half_period_residual=True,
        tol=1e-10,
        independent_tol=1e-6,
        state_step_cap=0.2,
        require_monotone_decrease=False,
        max_iter=80,
    )


def _member_orbit(member: dict, header: dict) -> BCR4BPPeriodicOrbit:
    """Reconstruct a BCR4BPPeriodicOrbit from a family jsonl member row."""
    sys_bcr = bcr4bp.BCR4BPSystem(
        mu=header["mu_em"],
        mu_sun=member["mu_sun_value"],
        a_sun_nondim=header["a_sun_nondim"],
        omega_sun_nondim=header["omega_sun_nondim"],
    )
    ic = np.array(
        [member["x0"], member["y0"], member["z0"], member["vx0"], member["vy0"], member["vz0"]],
        dtype=np.float64,
    )
    return BCR4BPPeriodicOrbit(
        state_initial=ic,
        period_nondim=member["T_TU"],
        sun_commensurate_n=member["sun_commensurate_n"],
        sun_phase_drift=member["sun_phase_drift"],
        converged=True,
        corrector_residual=member["corrector_residual"],
        independent_closure_residual=member["independent_closure_residual"],
        n_iter=1,
        system=sys_bcr,
        free_vars=(0, IDX_YDOT, 6),
        residual_indices=(IDX_Y, IDX_XDOT, IDX_ZDOT),
        is_half_period_residual=True,
        notes="reconstructed from family jsonl",
    )


def _run_gauntlet(cid: str, orbit: BCR4BPPeriodicOrbit) -> dict:
    v0 = run_v0_bcr4bp(cid, orbit)
    v1 = run_v1_bcr4bp(cid, orbit)
    v2 = run_v2_bcr4bp(cid, orbit, require_commensurate=True)
    v3_pass: bool | None = None
    if v2.passes_v2_bcr4bp:
        v3 = run_v3_bcr4bp(cid, orbit, v2_verdict=v2, n_cycles=3)
        v3_pass = v3.passes_v3_bcr4bp
    return {
        "id": cid,
        "sun_phase_drift": round(float(orbit.sun_phase_drift), 4),
        "quasi_periodic": v0.quasi_periodic,
        "V0": v0.passes_v0_bcr4bp,
        "V1": v1.passes_v1_bcr4bp,
        "V1_indep_km": round(float(v1.independent_closure_kms), 4),
        "V2": v2.passes_v2_bcr4bp,
        "V2_reason": v2.reason,
        "V3": v3_pass,
    }


def _load_family(path: Path) -> tuple[dict, list[dict]]:
    rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    header = next(r for r in rows if r.get("row_type") == "header")
    members = [r for r in rows if r.get("row_type") == "member"]
    return header, members


def main() -> None:
    results: list[dict] = []

    print("== POL1 (sourced Andreu/Rosales-Jorba L1 substitute, n=1 commensurate) ==")
    pol1 = _close_pol1()
    r = _run_gauntlet("andreu-pol1-bcr4bp", pol1)
    results.append(r)
    print(json.dumps(r, indent=1))

    for tag, fname in (
        ("303-L1", "bcr4bp_l1_family_303.jsonl"),
        ("304-halo", "bcr4bp_halo_family_304.jsonl"),
    ):
        path = _DATA / fname
        if not path.exists():
            print(f"\n== {tag}: {fname} MISSING ==")
            continue
        header, members = _load_family(path)
        # Probe the Andreu-mu_sun endpoint (last member) — the published model.
        last = members[-1]
        cid = f"{tag}-musun-andreu-step{last['step_idx']}"
        orbit = _member_orbit(last, header)
        print(f"\n== {tag} (Andreu mu_sun endpoint, step {last['step_idx']}) ==")
        rr = _run_gauntlet(cid, orbit)
        results.append(rr)
        print(json.dumps(rr, indent=1))

    print("\n== SUMMARY (V0-V3; V4-SEM real-eph + V5 deferred per design draft) ==")
    for r in results:
        v3 = r["V3"]
        v3s = "n/a" if v3 is None else ("PASS" if v3 else "FAIL")
        reason = r["V2_reason"] or "ok"
        print(
            f"  {r['id']:<40} drift={r['sun_phase_drift']:>8} qp={r['quasi_periodic']!s:<5} "
            f"V0={r['V0']!s:<5} V1={r['V1']!s:<5} V2={r['V2']!s:<5} ({reason}) V3={v3s}"
        )
    print(
        "\nNOTE: any V0-V3 PASS is a KNOWN-REPRODUCTION candidate (first_published=Andreu);\n"
        "flag for HUMAN REVIEW with the proposed row + attribution. NEVER self-admit."
    )


if __name__ == "__main__":
    main()
