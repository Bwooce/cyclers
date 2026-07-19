"""#655 -- Rhea-Titan 3D-closure inclination-extension probe.

Matches #575's own C2 "inclination-extension" method EXACTLY, reusing
``scripts/probe_572_titan_iapetus_3d_closure.py``'s ``evaluate_point`` /
``sweep_node_alignment`` node-alignment-search machinery verbatim (only the
module-level ``ANCHOR``/``FLYBY``/``INCLINATION_DEG``/``CANDIDATES``
constants change -- the functions themselves reference these globals, so no
function-body edits were made; the algebra, smoke test, gate reuse, and
basin-refinement discipline are byte-identical to #572).

Why this check is needed (per #655's dispatch instruction): the base #563
symmetric-closure construction is CIRCULAR-COPLANAR. #575 found that a
genuine, machine-precision, repeat-to-machine-precision coplanar Titan-
Iapetus symmetric family (9 closures) completely fails to survive extension
to Iapetus's REAL ~15.5 deg mutual inclination (0/6, 0/9 basins found do not
repeat as multi-cycle cyclers -- see #575's own empty_regions.jsonl stamp).
The #655 base run found 3 gate-passing coplanar symmetric closures for
Rhea-Titan (``data/enumerate_655_saturn_rhea_titan_symmetric_closures.jsonl``,
Rhea-anchored direction) -- before calling these "genuine repeating cyclers"
this probe checks whether they survive the SAME kind of inclination
extension, at Rhea-Titan's OWN real (much smaller) mutual inclination.

Inclination estimate (documented conservative UPPER BOUND, core/satellites.py
carries no inclination field for any Saturn moon -- same situation #572 was
in for Iapetus): JPL SSD "Planetary Satellite Mean Elements"
(ssd.jpl.nasa.gov/sats/elem/, SAT441 ephemeris, mean elements referred to
each moon's own local Laplace plane, accessed 2026-07-19) lists Rhea's
inclination to its local Laplace plane as ~0.345 deg and Titan's as
~0.348-0.35 deg (rounded to ~0.3 deg by the page's own summary table). Since
both moons orbit close enough to Saturn that their local Laplace planes are
"effectively identical to the Saturn equator" (interior-to-Titan moons; same
source), a fair worst-case bound on the Rhea-Titan MUTUAL inclination is the
spherical-triangle-inequality upper limit i_Rhea + i_Titan =~ 0.345+0.349 =~
0.69 deg -- rounded UP to 0.7 deg here, a small but non-zero conservative
estimate, roughly 22x SMALLER than the 15.5 deg used for Titan-Iapetus. This
is the single largest qualitative difference between the two pairs and the
entire reason #654 flagged Rhea-Titan (and Dione-Rhea, Tethys-Dione) as
plausibly Uranus-class rather than Titan-Iapetus-class.

Discipline: NO catalogue writeback, NO new gate logic (reuses #572's
``candidate_passes_physical_gate``/``GATE_RESIDUAL_KMS`` verbatim), read-only
on ``data/enumerate_655_saturn_rhea_titan_symmetric_closures.jsonl`` (the 3
candidates are LOADED from that file, not hand-transcribed, per project
discipline).

Run as::

    uv run python scripts/probe_655_rhea_titan_3d_closure.py
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import probe_572_titan_iapetus_3d_closure as p572  # noqa: E402
from scan_558_uranus_all_pairs_offset_sweep import GATE_RESIDUAL_KMS  # noqa: E402

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES  # noqa: E402
from cyclerfinder.search.discovery_campaign import _mean_motion_rad_day  # noqa: E402
from cyclerfinder.search.physical_sanity import (  # noqa: E402
    DEFAULT_MIN_USEFUL_BEND_DEG,
    candidate_passes_physical_gate,
)

DATA_DIR = ROOT / "data"
ENUM_655_PATH = DATA_DIR / "enumerate_655_saturn_rhea_titan_symmetric_closures.jsonl"
OUT_PATH = DATA_DIR / "probe_655_rhea_titan_3d_closure.jsonl"

# #572's module-level globals -- overridden here (this IS the genericization
# mechanism: evaluate_point/sweep_node_alignment reference these names, not
# hardcoded literals) rather than a function-signature rewrite.
p572.PRIMARY = "Saturn"
p572.ANCHOR = "Rhea"
p572.FLYBY = "Titan"
p572.INCLINATION_DEG = 0.7  # conservative upper-bound estimate; see module docstring


def load_655_candidates() -> list[dict[str, Any]]:
    """The 3 unique Rhea-anchored #655 symmetric-closure survivors, loaded
    directly from the committed enumeration output (never hand-transcribed)."""
    cands: list[dict[str, Any]] = []
    with ENUM_655_PATH.open(encoding="utf-8") as fh:
        for line in fh:
            d = json.loads(line)
            if d.get("kind") != "pass" or d.get("anchor") != "Rhea":
                continue
            label = f"n{d['n_commensurate_int']}_nrev{d['n_rev']}_rel{d['rel_offset_deg']:.0f}"
            cands.append(
                {
                    "label": label,
                    "rel_offset_deg": d["rel_offset_deg"],
                    "tof_scale": d["tof_days"] / _sqrt_papb(),
                    "n_rev": tuple(d["n_rev"]),
                    "coplanar_residual_kms": d["residual_kms"],
                    "coplanar_iapetus_vinf_kms": max(d["vinf_per_encounter_kms"]),
                    "coplanar_bend_deg": min(d["max_bend_deg_per_encounter"]),
                }
            )
    return cands


def _sqrt_papb() -> float:
    import math

    mu = PRIMARIES["Saturn"]
    sat_a = SATELLITES["Rhea"]
    sat_b = SATELLITES["Titan"]
    n_a = _mean_motion_rad_day(mu, sat_a.sma_km)
    n_b = _mean_motion_rad_day(mu, sat_b.sma_km)
    p_a = 2.0 * math.pi / n_a
    p_b = 2.0 * math.pi / n_b
    return math.sqrt(p_a * p_b)


def main() -> int:
    sha = p572._git_sha()
    print(f"[655] Rhea-Titan 3D-closure inclination-extension probe -- sha={sha}", flush=True)
    print(f"[655] inclination_deg={p572.INCLINATION_DEG} (conservative upper bound)", flush=True)

    print("[655] smoke test: iapetus_state_3d reduces to _moon_state at inc=0 ...", flush=True)
    smoke_ok = p572._smoke_test_reduction()
    print(f"[655] smoke test PASS: {smoke_ok}", flush=True)
    if not smoke_ok:
        print("[655] ABORTING -- 3D state generator does not reduce correctly.", flush=True)
        return 1

    candidates = load_655_candidates()
    assert len(candidates) == 3, f"expected 3 Rhea-anchored #655 survivors, got {len(candidates)}"

    out_records: list[dict[str, Any]] = [
        {
            "_meta": True,
            "task": "#655 Rhea-Titan 3D-closure inclination-extension probe",
            "git_sha": sha,
            "inclination_deg": p572.INCLINATION_DEG,
            "inclination_provenance": (
                "conservative upper bound = i_Rhea + i_Titan to Saturn's local Laplace plane "
                "(JPL SSD sats/elem SAT441, ~0.345+0.349=~0.69deg, rounded up to 0.7deg); "
                "core/satellites.py carries no inclination field, same situation #572 was in"
            ),
            "gate_residual_kms": GATE_RESIDUAL_KMS,
            "min_useful_bend_deg": DEFAULT_MIN_USEFUL_BEND_DEG,
            "smoke_test_reduction_pass": smoke_ok,
            "n_candidates": len(candidates),
        }
    ]

    verdicts: list[dict[str, Any]] = []
    seq = (p572.ANCHOR, p572.FLYBY, p572.ANCHOR)
    for cand in candidates:
        print(f"[655] --- {cand['label']} ---", flush=True)
        t0 = time.time()
        sweep = p572.sweep_node_alignment(cand, n_omega=3600)
        elapsed = time.time() - t0
        sweep["elapsed_s"] = elapsed
        print(
            f"[655]   n_feasible_omega={sweep['n_feasible_omega_points']}/{sweep['n_omega_grid']}  "
            f"n_basins={sweep['n_basins_found']}  "
            f"geom_errors={sweep['n_geometry_errors_encountered']} "
            f"(resolved_by_retry={sweep['n_geometry_errors_resolved_by_retry']})  "
            f"({elapsed:.1f}s)",
            flush=True,
        )

        basin_evals: list[dict[str, Any]] = []
        for b in sweep["basins"]:
            residual_near = b["residual_kms"] < GATE_RESIDUAL_KMS
            gate_pass, gate_verdicts = candidate_passes_physical_gate(
                seq, tuple(b["vinf_kms"]), min_useful_bend_deg=DEFAULT_MIN_USEFUL_BEND_DEG
            )
            bends = [v.max_bend_deg for v in gate_verdicts]
            closure = bool(residual_near and gate_pass)
            basin_evals.append(
                {
                    "omega_deg": b["omega_deg"],
                    "tof_scale": b["tof_scale"],
                    "residual_kms": b["residual_kms"],
                    "vinf_kms": b["vinf_kms"],
                    "residual_near_coplanar": residual_near,
                    "physical_gate_pass": gate_pass,
                    "max_bend_deg_per_encounter": bends,
                    "closure": closure,
                }
            )
            vinf_str = [f"{v:.3f}" for v in b["vinf_kms"]]
            bends_str = [f"{x:.2f}" for x in bends]
            print(
                f"[655]     basin omega={b['omega_deg']:7.3f} tof_scale={b['tof_scale']:.4f} "
                f"residual={b['residual_kms']:.6f}  vinf={vinf_str}  "
                f"bends={bends_str}  "
                f"near_gate={residual_near} phys_gate={gate_pass} CLOSURE={closure}",
                flush=True,
            )

        closing_basins = [be for be in basin_evals if be["closure"]]
        closure = len(closing_basins) > 0
        if not sweep["basins"]:
            print("[655]   NO feasible 3D point found at ANY node alignment tried.", flush=True)
            verdict = {
                "label": cand["label"],
                "closure_found": False,
                "reason": "no_feasible_lambert_point_at_any_omega",
            }
        else:
            best_closing = (
                min(closing_basins, key=lambda b: b["residual_kms"]) if closing_basins else None
            )
            print(f"[655]   >>> CLOSURE VERDICT for {cand['label']}: {closure} <<<", flush=True)
            verdict = {
                "label": cand["label"],
                "closure_found": closure,
                "n_basins_evaluated": len(basin_evals),
                "n_closing_basins": len(closing_basins),
                "best_closing_basin": best_closing,
                "all_basins": basin_evals,
            }
        verdicts.append(verdict)
        out_records.append({"kind": "candidate_result", **sweep})
        out_records.append({"kind": "candidate_verdict", **verdict})

    with OUT_PATH.open("w", encoding="utf-8") as fh:
        for rec in out_records:
            fh.write(json.dumps(rec, default=str) + "\n")
    print(f"[655] results written to {OUT_PATH}", flush=True)

    print("[655] === SUMMARY ===", flush=True)
    for v in verdicts:
        print(f"[655]   {v['label']}: closure_found={v['closure_found']}", flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
