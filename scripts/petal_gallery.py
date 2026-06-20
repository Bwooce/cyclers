"""#406 petal-plot gallery: every CR3BP cycler in the catalogue as a high-accuracy
rotating-frame rosette.

The non-keplerian (CR3BP) catalogue rows carry a sourced/derived (mu, state_nd, T)
in ``orbit_elements.cr3bp``. In the rotating frame these orbits are NOT simple
ellipses — a (k1,k2) resonant cycler that makes several radial excursions per cycle
traces a multi-lobed rosette (the "petal"/"flower" / star shapes seen in the 3D
viz). This tool propagates each at high accuracy (DOP853, tight tol) and renders a
gallery, annotating each panel with honest closure + Jacobi-drift numbers so the plot
doubles as visual QA of the exotic orbits.

Run: uv run python scripts/petal_gallery.py [--out docs/figures/petal_gallery.png]
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp

_CATALOGUE = Path("data/catalogue.yaml")


def _short_family(family: str | None, row_id: str) -> str:
    """Trim the sourced family string to its leading resonance label."""
    if not family:
        return row_id
    # "(1,1) prograde Earth-Moon cycler — long-period branch (C11b)" -> first clause
    head = family.split(" prograde")[0].split(" Earth-Moon")[0]
    tag = ""
    if "(C11" in family:
        tag = " " + family[family.index("(C11") :].split(")")[0] + ")"
    return head + tag


def _collect_cr3bp_rows() -> list[dict]:
    rows = yaml.safe_load(_CATALOGUE.read_text())
    out: list[dict] = []
    for row in rows:
        oe = (row.get("orbit_elements") or {}).get("cr3bp")
        if not oe or not oe.get("state_nd"):
            continue
        if oe.get("mass_ratio") is None or oe.get("period_nd") is None:
            continue
        out.append(
            {
                "id": row["id"],
                "family": _short_family(oe.get("family"), row["id"]),
                "mu": float(oe["mass_ratio"]),
                "state_nd": [float(v) for v in oe["state_nd"]],
                "period_nd": float(oe["period_nd"]),
            }
        )
    return out


def _propagate(mu: float, state_nd: list[float], period_nd: float):
    s0 = np.asarray(state_nd, dtype=float)
    c0 = cr3bp.jacobi_constant(s0, mu)
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, period_nd),
        s0,
        args=(mu,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        t_eval=np.linspace(0.0, period_nd, 4000),
    )
    sf = sol.y[:, -1]
    closure_nd = float(np.linalg.norm(sf[:6] - s0))
    drift = abs(cr3bp.jacobi_constant(sf, mu) - c0)
    return sol.y[0], sol.y[1], closure_nd, drift


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="docs/figures/petal_gallery.png")
    args = ap.parse_args()

    rows = _collect_cr3bp_rows()
    if not rows:
        print("no CR3BP rows with state_nd found")
        return
    n = len(rows)
    cols = 4
    nrows = math.ceil(n / cols)
    fig, axes = plt.subplots(nrows, cols, figsize=(3.4 * cols, 3.4 * nrows))
    axes = np.atleast_1d(axes).ravel()

    for ax, row in zip(axes, rows, strict=False):
        x, y, closure_nd, drift = _propagate(row["mu"], row["state_nd"], row["period_nd"])
        mu = row["mu"]
        ax.plot(x, y, lw=0.9, color="#7fa8e0")
        ax.plot(-mu, 0.0, "o", color="#2a6299", ms=7)  # Earth
        ax.plot(1.0 - mu, 0.0, "o", color="#999999", ms=4)  # Moon
        ax.set_aspect("equal")
        ax.set_title(row["family"], fontsize=8)
        t_d = row["period_nd"] * 375190.0 / 86400.0  # ~EM TU -> days
        ax.text(
            0.5,
            -0.02,
            f"T={row['period_nd']:.2f} TU (~{t_d:.0f} d)\n"
            f"closure {closure_nd * 384400:.0f} km · ΔC {drift:.1e}",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=6,
            color="#555",
        )
        ax.set_xticks([])
        ax.set_yticks([])
    for ax in axes[n:]:
        ax.axis("off")

    fig.suptitle(
        "Catalogue CR3BP cyclers — Earth-Moon rotating-frame rosettes (DOP853, high accuracy)",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=120)
    print(f"saved {out}  ({n} rosettes)")
    for row in rows:
        print(f"  {row['id']}: {row['family']}")


if __name__ == "__main__":
    main()
