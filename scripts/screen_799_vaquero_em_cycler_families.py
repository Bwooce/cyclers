"""#799: reproduce Vaquero 2013's Earth-Moon 2:1 / 3:1 periodic-cycler families.

Direct CR3BP family continuation across Vaquero's own printed Jacobi ranges
(Sec. 4.4.7, p. 171: 2:1 C in [1.98, 2.66]; 3:1 C in [2.54, 3.13]), seeded
from the Tisserand-matched two-body apoapsis construction in
``cyclerfinder.search.vaquero_em_cyclers``. Writes the full per-member
record to ``data/found/799_vaquero_em_cycler_families/results.json`` and
prints a summary including the four printed-endpoint TOF comparisons (the
only digit-grade golden values Vaquero prints for these families).

Named ``screen_*`` (not ``run_*``): a fixed two-family published-result
reproduction with no region_id/n_points sweep-region concept to preflight
(same category as ``scripts/screen_716_*``), not a catalogue-region
discovery sweep.

Foreground, single-process. Runtime: a few minutes.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path

import cyclerfinder.search.vaquero_em_cyclers as vec

OUT_DIR = (
    Path(__file__).resolve().parent.parent / "data" / "found" / ("799_vaquero_em_cycler_families")
)

#: Seed energies: mid-range grid points chosen so the printed endpoints land
#: exactly on the d_jacobi=0.01 continuation grid.
SEEDS = [
    ("2:1", 2, 2.30, vec.VAQUERO_C_RANGE_21, vec.VAQUERO_TOF_DAYS_21),
    ("3:1", 3, 2.80, vec.VAQUERO_C_RANGE_31, vec.VAQUERO_TOF_DAYS_31),
]


def main() -> None:
    system = vec.earth_moon_system()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out: dict = {
        "task": "#799",
        "generated_utc": datetime.now(UTC).isoformat(),
        "mu": system.mu,
        "l_km": system.l_km,
        "t_s": system.t_s,
        "families": [],
    }
    for label, p, c_seed, c_range, tof_prints in SEEDS:
        print(f"\n=== {label} family: seed C={c_seed}, target C in {c_range} ===", flush=True)
        rep = vec.reproduce_vaquero_family(
            system, p, c_seed=c_seed, c_range=c_range, label=f"vaquero_{p}to1"
        )
        print(
            f"  branches: down stop={rep.branch_down.stop_reason} "
            f"({len(rep.branch_down.members)} members), "
            f"up stop={rep.branch_up.stop_reason} "
            f"({len(rep.branch_up.members)} members); merged={len(rep.members)}",
            flush=True,
        )
        fam: dict = {
            "label": label,
            "p": p,
            "c_seed": c_seed,
            "c_range": list(c_range),
            "seed_geometry": asdict(rep.seed_geom),
            "stop_down": str(rep.branch_down.stop_reason),
            "stop_up": str(rep.branch_up.stop_reason),
            "n_members": len(rep.members),
            "members": [],
            "endpoint_tof_checks": [],
        }
        for m, r in zip(rep.members, rep.reports, strict=True):
            fam["members"].append(
                {
                    "jacobi": m.jacobi,
                    "x0": m.x0,
                    "ydot0": m.ydot0,
                    "period": m.period,
                    "nu_planar": r.nu_planar,
                    "abs_lambda": r.abs_lambda,
                    "nu_out_of_plane": r.nu_out_of_plane,
                    "stable_planar": r.stable_planar,
                    "perigee_km": r.perigee_km,
                    "apogee_km": r.apogee_km,
                    "moon_min_km": r.moon_min_km,
                    "a_two_body_lu": r.a_two_body_lu,
                    "tof_earth_moon_days": r.tof_earth_moon_days,
                    "perigee_in_band": r.perigee_in_band,
                    "tof_within_ceiling": r.tof_within_ceiling,
                    "moon_min_vs_l1_dist": r.moon_min_vs_l1_dist,
                    "moon_min_vs_l2_dist": r.moon_min_vs_l2_dist,
                    "crossing_residual": m.crossing_residual,
                    "radau_djacobi": m.radau_djacobi,
                }
            )
        # Printed-endpoint TOF golden checks.
        for c_print, tof_print in tof_prints.items():
            nearest = min(rep.reports, key=lambda r: abs(r.jacobi - c_print))
            reached = abs(nearest.jacobi - c_print) < 5e-3
            fam["endpoint_tof_checks"].append(
                {
                    "c_printed": c_print,
                    "tof_printed_days": tof_print,
                    "c_nearest_member": nearest.jacobi,
                    "endpoint_reached": reached,
                    "tof_ours_days": nearest.tof_earth_moon_days,
                    "tof_rel_err": abs(nearest.tof_earth_moon_days - tof_print) / tof_print,
                }
            )
            print(
                f"  endpoint C={c_print}: printed TOF={tof_print} d, "
                f"ours={nearest.tof_earth_moon_days:.3f} d at C={nearest.jacobi:.4f} "
                f"(reached={reached}, rel_err="
                f"{abs(nearest.tof_earth_moon_days - tof_print) / tof_print:.3%})",
                flush=True,
            )
        n_stable = sum(1 for r in rep.reports if r.stable_planar)
        n_band = sum(1 for r in rep.reports if r.perigee_in_band)
        c_vals = [r.jacobi for r in rep.reports]
        print(
            f"  members: C span [{min(c_vals):.4f}, {max(c_vals):.4f}], "
            f"{n_stable}/{len(rep.reports)} planar-stable, "
            f"{n_band}/{len(rep.reports)} perigee-in-LEO-GEO-band",
            flush=True,
        )
        worst_res = max(m.crossing_residual for m in rep.members)
        worst_dj = max(m.radau_djacobi for m in rep.members)
        print(f"  worst crossing residual {worst_res:.2e}, worst Radau dC {worst_dj:.2e}")
        assert not math.isnan(worst_res)
        out["families"].append(fam)

    path = OUT_DIR / "results.json"
    path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
