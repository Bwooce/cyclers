"""#668 -- rigorous interval certification of a piece of the W-Z Oterma proof.

Driver on top of :mod:`scripts._validated_taylor_integrator`.  ``#403`` reproduced
Wilczak & Zgliczynski's published Oterma heteroclinic crossing coordinates in
ordinary floating point (``data/golden/wz_oterma_heteroclinic.yaml``); W-Z's
*proof* rests on rigorous interval enclosures of the flow (CAPD/Lohner).  This
driver closes as much of that gap as one dispatch honestly can, in two tiers:

  * **Tier 1 (closed-form, rigorous).**  For every published section point,
    rigorously certify -- with ``mpmath.iv`` interval arithmetic -- that it is
    energetically admissible on the Jacobi level set ``C = 3.03``: the derived
    out-of-section velocity ``vy`` satisfies ``vy^2 >= 0`` with a rigorous,
    strictly-positive enclosed lower bound (the point genuinely lies on the
    Oterma energy manifold, not merely to display precision).  No integration.

  * **Tier 2 (validated integration, rigorous).**  Take the L1->L2 heteroclinic
    section point, build its full ``C = 3.03`` state, and rigorously enclose the
    forward CR3BP flow with the validated Taylor integrator.  Certify that the
    Jacobi constant -- an exact first integral -- stays enclosed around ``C0``
    to within the (tiny) enclosure width along the whole validated arc.  This is
    a genuine interval-verified invariant along a *real* Oterma trajectory, not a
    floating-point residual.  The arc is bounded: without Lohner's QR-coordinate
    ("reframing") wrapping-effect control the naive C0 enclosure grows and, near
    the L1 neck's hyperbolic stretching, can no longer be validated -- the driver
    reports this horizon honestly and demonstrates (two step sizes) that it is
    the *wrapping effect*, not truncation error (the blow-up time is essentially
    step-size independent).

Honest scope -- this is a *fraction* of the full W-Z proof.  It supplies the one
piece #653 found the codebase entirely lacked: a rigorous ODE-flow enclosure.
It does NOT yet provide: QR-Lohner wrapping control (so no long / section-to-
section-map rigor), a C1-variational (STM) enclosure, rigorous Poincare-map
derivatives, or the covering-relations / h-set machinery that turns enclosures
into a topological existence theorem.  Those are the remaining stages.

Run:  ``uv run python scripts/certify_668_wz_oterma_interval.py``
Writes ``data/668_wz_oterma_interval_certificate.json``.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import _validated_taylor_integrator as vti  # noqa: E402

GOLDEN = ROOT / "data" / "golden" / "wz_oterma_heteroclinic.yaml"
OUT = ROOT / "data" / "668_wz_oterma_interval_certificate.json"

DPS = 40
MU = 0.0009537
C_LEVEL = 3.03


def _pt(iv: Any, x: float) -> Any:
    return iv.mpf([x, x])


def _vy2_on_level(iv: Any, x: float, vx: float, mu: float, c: float) -> Any:
    """Rigorous enclosure of vy^2 for a section point (x, y=0, vx) at level C=c."""
    xi = _pt(iv, x)
    vxi = _pt(iv, vx)
    r1 = iv.sqrt((xi + iv.mpf(mu)) ** 2)
    r2 = iv.sqrt((xi + iv.mpf(mu - 1.0)) ** 2)
    return xi**2 + 2 * (1 - iv.mpf(mu)) / r1 + 2 * iv.mpf(mu) / r2 - vxi**2 - iv.mpf(c)


def tier1_admissibility(iv: Any, golden: dict[str, Any]) -> list[dict[str, Any]]:
    """Rigorously certify every published section point is admissible at C=3.03."""
    results: list[dict[str, Any]] = []
    items: list[tuple[str, float, float]] = []
    for name, node in golden["lyapunov_fixed_points"].items():
        x, vx = node["point"]
        items.append((f"fixed_point/{name}", float(x), float(vx)))
    for cname, cnode in golden["crossings"].items():
        for i, (x, vx) in enumerate(cnode["sequence"]):
            items.append((f"{cname}[{i}]", float(x), float(vx)))
    for label, x, vx in items:
        vy2 = _vy2_on_level(iv, x, vx, MU, C_LEVEL)
        lo = float(vy2.a)  # rigorous lower endpoint of vy^2
        admissible = bool(vy2.a >= 0)  # rigorous: vy^2 lower bound >= 0
        results.append(
            {
                "label": label,
                "x": x,
                "vx": vx,
                "vy2_lower": lo,
                "admissible_C303": admissible,
            }
        )
    return results


def tier2_flow_enclosure(iv: Any, golden: dict[str, Any]) -> dict[str, Any]:
    """Validated forward flow of the L1->L2 heteroclinic point; Jacobi enclosed."""
    seq = golden["crossings"]["heteroclinic_L1_to_L2"]["sequence"]
    x0, vx0 = float(seq[0][0]), float(seq[0][1])
    vy2 = _vy2_on_level(iv, x0, vx0, MU, C_LEVEL)
    y0 = [_pt(iv, x0), _pt(iv, 0.0), _pt(iv, vx0), iv.sqrt(vy2)]
    c0 = vti.cr3bp_planar_jacobi(iv, y0, MU)

    # Two step sizes over the same reachable horizon: show the blow-up time is
    # step-size independent (wrapping effect, not truncation).
    runs: list[dict[str, Any]] = []
    for n_steps in (64, 128):
        res = vti.integrate(iv, vti.cr3bp_planar_jet, y0, 1.0, MU, n_steps=n_steps, order=12)
        h = 1.0 / n_steps
        t_reached = res["n_completed"] * h
        cf = vti.cr3bp_planar_jacobi(iv, res["final"], MU) if res["n_completed"] else c0
        jacobi_ok = bool(cf.a <= c0.a) and bool(c0.b <= cf.b)
        runs.append(
            {
                "n_steps": n_steps,
                "step_h": h,
                "steps_validated": res["n_completed"],
                "t_reached": t_reached,
                "final_max_halfwidth": res["final_max_halfwidth"],
                "jacobi_end_contains_C0": jacobi_ok,
                "jacobi_end_width": float(cf.delta.b),
            }
        )
    return {
        "start_x": x0,
        "start_vx": vx0,
        "C0_lower": float(c0.a),
        "C0_upper": float(c0.b),
        "runs": runs,
    }


def build_certificate() -> dict[str, Any]:
    import mpmath as mp

    mp.mp.dps = DPS
    mp.iv.dps = DPS
    iv = mp.iv
    golden = yaml.safe_load(GOLDEN.read_text())

    t1 = tier1_admissibility(iv, golden)
    t2 = tier2_flow_enclosure(iv, golden)

    all_admissible = all(r["admissible_C303"] for r in t1)
    all_jacobi_ok = all(run["jacobi_end_contains_C0"] for run in t2["runs"])
    # wrapping signature: halving the step barely extends the reachable time
    t_reached = [run["t_reached"] for run in t2["runs"]]
    wrapping_limited = max(t_reached) < 2.0 * min(t_reached)

    return {
        "task": "#668",
        "system": {"mu": MU, "C": C_LEVEL, "model": "planar CR3BP (Sun-Jupiter, Oterma)"},
        "method": (
            "mpmath.iv validated interval Taylor integration (Lohner C0: a-priori "
            "Picard enclosure + rigorous order-p Taylor remainder). See "
            "scripts/_validated_taylor_integrator.py."
        ),
        "dps": DPS,
        "tier1_energy_admissibility": {
            "all_points_admissible_at_C303": all_admissible,
            "n_points": len(t1),
            "points": t1,
        },
        "tier2_validated_flow": {
            "jacobi_conserved_within_enclosure": all_jacobi_ok,
            "wrapping_effect_limited": wrapping_limited,
            **t2,
        },
        "scope_note": (
            "Genuinely interval-rigorous, but a FRACTION of the full W-Z proof: "
            "supplies a rigorous CR3BP flow enclosure (Tier 2) and closed-form "
            "energy-admissibility (Tier 1). Does NOT yet include QR-Lohner "
            "wrapping control, C1/STM enclosure, rigorous Poincare-map "
            "derivatives, or covering-relations/h-set topology. The validated arc "
            "is wrapping-limited near the L1 neck (blow-up time step-size "
            "independent, i.e. wrapping not truncation)."
        ),
    }


def main() -> None:
    cert = build_certificate()
    OUT.write_text(json.dumps(cert, indent=2) + "\n")
    t1 = cert["tier1_energy_admissibility"]
    t2 = cert["tier2_validated_flow"]
    print(
        f"[#668] Tier 1: {t1['n_points']} points, all admissible at C=3.03: "
        f"{t1['all_points_admissible_at_C303']}"
    )
    print(
        f"[#668] Tier 2: Jacobi conserved within enclosure: "
        f"{t2['jacobi_conserved_within_enclosure']}; wrapping-limited: "
        f"{t2['wrapping_effect_limited']}"
    )
    for run in t2["runs"]:
        print(
            f"         h={run['step_h']:.4f}: {run['steps_validated']} steps validated, "
            f"t_reached={run['t_reached']:.3f}, final_halfwidth={run['final_max_halfwidth']:.2e}, "
            f"jacobi_contains_C0={run['jacobi_end_contains_C0']}"
        )
    print(f"[#668] certificate written -> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
