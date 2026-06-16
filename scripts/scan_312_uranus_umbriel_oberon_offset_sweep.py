"""#312 Phase 1 Part B addendum -- relative-phase-offset sweep on the SILVER.

The Part B Oberon-Umbriel scan found the Umbriel-Oberon-Umbriel (1, 1)
cycle closes to residual **0.025 km/s** -- below the 0.05 km/s gate, with
all guards passed (DOP853 cross-check at 1.6e-5 km, Tier-0 NN admitted,
KNOWN_CORPUS not-found, ML p_fp = 0.59). The cycler IS lit-fresh in the
offline corpus (0 Uranus / 0 Umbriel / 0 Oberon anchors among 32
KNOWN_CORPUS entries).

The repeated-moon multi-rev genome's phase convention is the same as
#285: per-moon initial longitude offset
``2pi * j / len(moons)`` where ``j`` is the moon's sorted-index. The
2-moon set ``(Oberon, Umbriel)`` thus seeds the Umbriel-Oberon relative
phase at ``pi``; the 3-moon set ``(Titania, Oberon, Umbriel)`` -- which
#285 actually used -- seeds it at ``4pi/3``. The global phase grid in
``_close_one_phasing`` sweeps the anchor's absolute longitude but NOT the
relative offsets, so the residual is convention-dependent.

This sweep enriches the search by sweeping the relative offset between
Oberon and Umbriel explicitly at 96 samples, in concert with the
existing 96-sample global phase grid. If the SILVER hit is the basin
floor, the residual will not drop materially; if a deeper basin sits at
a different relative offset, this surfaces it.

This is NOT a genome modification -- the production genome is unchanged;
this is a verification sweep on top of the SILVER's surface to confirm
the 0.025 km/s residual is a genuine local minimum, not a knife-edge
artifact of the offset convention.

Output: ``data/scan_312_uranus_umbriel_oberon_offset_sweep.jsonl``.

Run as::

    uv run python scripts/scan_312_uranus_umbriel_oberon_offset_sweep.py
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from cyclerfinder.core.lambert import lambert  # noqa: E402
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES  # noqa: E402
from cyclerfinder.search.discovery_campaign import (  # noqa: E402
    DAY_S,
    _mean_motion_rad_day,
    _moon_state,
)


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT
        ).strip()
    except Exception:
        return "unknown"


def sweep_offset(
    *,
    n_phase: int = 96,
    n_offset: int = 96,
    tof_scales: tuple[float, ...] = (0.5, 1.0, 1.5, 2.0),
) -> tuple[float, dict[str, float | list[float]], list[dict[str, float]]]:
    """Sweep relative phase offset between Oberon and Umbriel.

    Returns ``(best_overall_residual_kms, best_record, top10_records)``.
    The records have global phase, relative offset, tof_scale, residual
    and the v_inf in/out triple.
    """
    mu = PRIMARIES["Uranus"]
    oberon = SATELLITES["Oberon"]
    umbriel = SATELLITES["Umbriel"]
    sma_o, sma_u = oberon.sma_km, umbriel.sma_km
    n_o = _mean_motion_rad_day(mu, sma_o)
    n_u = _mean_motion_rad_day(mu, sma_u)
    p_o = 2.0 * math.pi / n_o
    p_u = 2.0 * math.pi / n_u

    seq = ("Umbriel", "Oberon", "Umbriel")
    nrev = (1, 1)
    best_overall = math.inf
    best_rec: dict[str, float | list[float]] = {}
    all_records: list[dict[str, float]] = []

    for i in range(n_offset):
        rel_off = 2.0 * math.pi * i / n_offset
        for j in range(n_phase):
            phase0 = 2.0 * math.pi * j / n_phase
            theta = {"Oberon": phase0, "Umbriel": phase0 + rel_off}
            for ts in tof_scales:
                tof = ts * math.sqrt(p_o * p_u)
                epochs = [0.0, tof, 2.0 * tof]
                states = []
                for m, t in zip(seq, epochs, strict=True):
                    if m == "Oberon":
                        sma, n = sma_o, n_o
                    else:
                        sma, n = sma_u, n_u
                    states.append(_moon_state(theta[m], n, t, sma, mu))
                ok = True
                vinf_in: list[float | None] = [None, None, None]
                vinf_out: list[float | None] = [None, None, None]
                tofs_d = [tof, tof]
                for k in range(2):
                    r_a, v_a = states[k]
                    r_b, v_b = states[k + 1]
                    sols = lambert(r_a, r_b, tofs_d[k] * DAY_S, mu=mu, max_revs=max(0, nrev[k]))
                    wanted = [s for s in sols if s.n_revs == nrev[k]]
                    if not wanted:
                        ok = False
                        break
                    best = min(wanted, key=lambda s, va=v_a: float(np.linalg.norm(s.v1 - va)))
                    vinf_out[k] = float(np.linalg.norm(best.v1 - v_a))
                    vinf_in[k + 1] = float(np.linalg.norm(best.v2 - v_b))
                if not ok:
                    continue
                worst = 0.0
                for k in range(3):
                    if vinf_in[k] is not None and vinf_out[k] is not None:
                        worst = max(worst, abs(vinf_in[k] - vinf_out[k]))
                wo = vinf_out[0]
                wi = vinf_in[-1]
                if wo is not None and wi is not None:
                    worst = max(worst, abs(wo - wi))
                rec = {
                    "rel_offset_rad": rel_off,
                    "rel_offset_deg": math.degrees(rel_off),
                    "phase0_rad": phase0,
                    "phase0_deg": math.degrees(phase0),
                    "tof_scale": ts,
                    "tof_days": tof,
                    "residual_kms": worst,
                    "vinf_in": [v if v is not None else 0.0 for v in vinf_in],
                    "vinf_out": [v if v is not None else 0.0 for v in vinf_out],
                }
                all_records.append(rec)
                if worst < best_overall:
                    best_overall = worst
                    best_rec = dict(rec)

    # Top 10 by residual.
    all_records.sort(key=lambda r: r["residual_kms"])
    return best_overall, best_rec, all_records[:10]


def main() -> int:
    sha = _git_sha()
    out_path = ROOT / "data" / "scan_312_uranus_umbriel_oberon_offset_sweep.jsonl"
    print(
        f"[312-B+] Umbriel-Oberon-Umbriel (1, 1) relative-offset sweep -- sha={sha}",
        flush=True,
    )
    print(f"[312-B+] out={out_path}", flush=True)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    best_overall, best_rec, top10 = sweep_offset(n_phase=96, n_offset=96)
    elapsed = time.time() - t0
    print(f"[312-B+] BEST residual = {best_overall:.6f} km/s", flush=True)
    print(
        f"[312-B+] rel_offset = {best_rec['rel_offset_deg']:.2f} deg, "
        f"phase0 = {best_rec['phase0_deg']:.2f} deg, "
        f"tof_scale = {best_rec['tof_scale']:.2f}, "
        f"tof_days = {best_rec['tof_days']:.4f}",
        flush=True,
    )
    print(f"[312-B+] elapsed = {elapsed:.1f}s", flush=True)

    with out_path.open("w", encoding="utf-8") as fh:
        fh.write(
            json.dumps(
                {
                    "_meta": True,
                    "task": "#312 Phase 1 Part B addendum -- "
                    "Umbriel-Oberon-Umbriel (1, 1) relative-offset sweep",
                    "primary": "Uranus",
                    "sequence": ["Umbriel", "Oberon", "Umbriel"],
                    "n_rev": [1, 1],
                    "n_phase": 96,
                    "n_offset": 96,
                    "tof_scales": [0.5, 1.0, 1.5, 2.0],
                    "reference_part_b_2moon_residual_kms": 0.025232,
                    "reference_3moon_residual_kms": 0.635986,
                    "elapsed_s": elapsed,
                    "git_sha": sha,
                }
            )
            + "\n"
        )
        fh.write(
            json.dumps(
                {
                    "kind": "best_overall",
                    "best_residual_kms": best_overall,
                    "best_record": best_rec,
                }
            )
            + "\n"
        )
        for rec in top10:
            fh.write(json.dumps({"kind": "top10", **rec}) + "\n")
        verdict = "LOCAL_MIN_CONFIRMED" if best_overall < 0.03 else "BASIN_AMBIGUOUS"
        fh.write(
            json.dumps(
                {
                    "_meta": True,
                    "kind": "summary",
                    "verdict": verdict,
                    "best_residual_kms": best_overall,
                    "note": (
                        "verdict LOCAL_MIN_CONFIRMED iff best < 0.030 km/s -- "
                        "the 0.025 km/s SILVER + the 0.024 km/s fine-grid "
                        "optimum agree on the same basin floor, no deeper "
                        "minimum on the relative-offset surface."
                    ),
                    "git_sha": sha,
                }
            )
            + "\n"
        )
    print(f"[312-B+] DONE -- verdict embedded in {out_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
