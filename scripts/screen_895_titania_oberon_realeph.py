"""#895 staged driver: ballistic Titania-Oberon arcs in a URA111 force model.

Pre-registration: ``docs/notes/2026-10-04-895-titania-oberon-realeph.md`` section 1.

This is a targeted correction around one known seed (the `#890` orbit) at five pre-registered
epochs, not a sweep of a region, so it does not call ``preflight_search`` (the ratchet in
``tests/scripts/test_scripts_call_preflight.py`` covers ``scripts/run_*.py``; the `#890` driver
``screen_890_*.py`` is the same kind of script).

Stages (each resumable, each invocation bounded by ``--budget`` seconds of wall time):

* ``control``  force-model positive control P1 and the table interpolation check;
* ``p2``       the `#890` model orbit rebuilt in this frame (P2);
* ``arc``      homotopy for one epoch and cycle count (checkpointed after every step);
* ``verify``   criteria (c) to (f) on a converged arc;
* ``sunoff``   the Sun's effect (one Newton solve with the Sun off, from a converged arc).

Run with ``uv run python scripts/screen_895_titania_oberon_realeph.py <stage> [options]``.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np

import cyclerfinder.search.titania_oberon_realeph_895 as m

OUT = Path(__file__).resolve().parents[1] / "data" / "found" / "895_titania_oberon_realeph"
CACHE = Path(os.environ.get("CF895_CACHE", "/tmp/cf895_cache"))
REFINE_890 = OUT.parent / "890_titania_oberon_candidate" / "refine.json"

#: Pre-registered target dates (section 1.4), 00:00 TDB.
EPOCHS: dict[str, str] = {
    "E1": "2030-01-12 00:00:00 TDB",
    "E2": "2031-06-13 00:00:00 TDB",
    "E3": "2035-07-01 00:00:00 TDB",
    "E4": "2040-07-01 00:00:00 TDB",
    "E5": "2045-07-01 00:00:00 TDB",
}
FIT_MARGIN_D = 30.0
TABLE_MARGIN_D = 60.0
MAX_CYCLES = 12
NOMINAL_CYCLE_D = 123.2


def log(msg: str) -> None:
    print(f"[{datetime.now(UTC).strftime('%H:%M:%S')}Z] {msg}", flush=True)


def git_sha() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def dump(name: str, obj: dict[str, Any]) -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    obj = {"_meta": {"task": "#895", "git_sha": git_sha(), "stage": name}, **obj}
    path = OUT / name
    path.write_text(json.dumps(obj, indent=1, default=float) + "\n")
    log(f"wrote {path}")
    return path


# ------------------------------------------------------------------------------------------- #
# Epoch set-up: table, circles, conjunction
# ------------------------------------------------------------------------------------------- #


def epoch_setup(tag: str) -> tuple[m.EphemerisTable, m.Circles, float]:
    t_ref = m.et_of(EPOCHS[tag])
    t0 = -TABLE_MARGIN_D * m.DAY_S
    t1 = (MAX_CYCLES * NOMINAL_CYCLE_D + TABLE_MARGIN_D + 20.0) * m.DAY_S
    table = m.build_table(t_ref, t0, t1, cache_dir=CACHE)
    circles = m.fit_circles(
        table, -FIT_MARGIN_D * m.DAY_S, (MAX_CYCLES * NOMINAL_CYCLE_D + FIT_MARGIN_D) * m.DAY_S
    )
    t_c = circles.conjunction_near(0.0)
    return table, circles, t_c


def stage_control(args: argparse.Namespace) -> None:
    rows: list[dict[str, Any]] = []
    interp: list[dict[str, Any]] = []
    epochs_info: dict[str, Any] = {}
    variants: dict[str, dict[str, Any]] = {
        "full": {},
        "no_sun": {"sun": False},
        "no_j4": {"j4": 0.0},
        "french_j2": {"j2": m.J2_FRENCH_2024},
        "ura111_pole": {"pole": m.pole_vector(*m.POLE_URA111_RA_DEC_DEG)},
        "no_self_mass": {"self_mass": False},
        "no_zonal": {"j2": 0.0, "j4": 0.0},
        "former_repo_j2": {"j2": m.J2_FORMER_REPO, "j4": 0.0},
    }
    sp = m.load_kernels()
    for tag in EPOCHS:
        table, circles, t_c = epoch_setup(tag)
        ang = math.degrees(math.acos(float(np.clip(-np.dot(circles.ez, m.POLE_IAU), -1, 1))))
        epochs_info[tag] = {
            "t_ref_et": table.t_ref_et,
            "t_c_s": t_c,
            "t_c_tdb": m.tdb_string(table.t_ref_et + t_c),
            "circles": circles.to_dict(),
            "angle_ez_from_minus_iau_pole_deg": ang,
        }
        log(f"{tag}: conjunction {epochs_info[tag]['t_c_tdb']}, ez vs -pole {ang:.4f} deg")
        # interpolation check at mid-knots over the first 40 days
        for k, naif in enumerate(m.NAIF_IDS):
            worst = 0.0
            model = m.full_model(table)
            for i in range(200):
                t = table.t0 + (i * 17 + 0.5) * table.h
                ref = np.asarray(sp.spkgeo(naif, table.t_ref_et + t, "J2000", 799)[0])
                p, _ = model.body_states(t)
                worst = max(worst, float(np.linalg.norm(p[k] - ref[:3])))
            interp.append(
                {"epoch": tag, "body": m.BODY_NAMES[k], "max_mid_knot_err_m": worst * 1e3}
            )
        for name, kw in variants.items():
            if args.quick and name not in ("full", "no_zonal", "former_repo_j2"):
                continue
            for k in m.MOONS:
                err = m.moon_positive_control(table, k, t_c, (30.0, 123.0), **kw)
                rows.append(
                    {"epoch": tag, "variant": name, "moon": m.BODY_NAMES[k], "err30_km": err[0],
                     "err123_km": err[1]}
                )  # fmt: skip
            log(
                f"{tag} {name}: "
                + ", ".join(
                    f"{r['moon']} {r['err30_km']:.2f}/{r['err123_km']:.1f}" for r in rows[-5:]
                )
            )
    full = [r for r in rows if r["variant"] == "full"]
    p1_pass = all(r["err30_km"] <= 10.0 and r["err123_km"] <= 40.0 for r in full)
    disc = {}
    for v in ("no_zonal", "former_repo_j2"):
        sel = [r for r in rows if r["variant"] == v and r["moon"] in ("Miranda", "Ariel")]
        disc[v] = all(r["err30_km"] > 10.0 for r in sel)
    dump(
        "control.json",
        {
            "criterion": "P1: every moon, every epoch: <= 10 km at 30 d and <= 40 km at 123 d "
            "(variant 'full'); discrimination: Miranda and Ariel > 10 km at 30 d for "
            "'no_zonal' and 'former_repo_j2' at every epoch",
            "P1_pass": p1_pass,
            "discrimination_pass": disc,
            "epochs": epochs_info,
            "interpolation": interp,
            "rows": rows,
        },
    )
    log(f"P1 pass {p1_pass}; discrimination {disc}")


def stage_p2(args: argparse.Namespace) -> None:
    c = m.registry_890()
    s890 = np.asarray(json.loads(REFINE_890.read_text())["refined_state"])
    y0 = m.state_890_to_inertial(c, s890)
    t_cyc = c.cycle_s
    out: dict[str, Any] = {"cycle_d": t_cyc / m.DAY_S}
    for label, model in (
        ("model_890", m.model_890(c)),
        ("oberon_about_uranus", m.model_890(c, oberon_about_barycentre=False)),
    ):
        res: dict[str, Any] = {}
        for rtol in (1e-12, 1e-13):
            half = m.propagate(model, 0.0, t_cyc / 2, y0, rtol=rtol, atol_r=1e-10, atol_v=1e-16)
            rot = m.inertial_to_rot_890(c, t_cyc / 2, half.state)
            full = m.propagate(
                model, 0.0, t_cyc, y0, rtol=rtol, atol_r=1e-10, atol_v=1e-16, record=True,
                max_steps=40000, allow_incomplete=True,
            )  # fmt: skip
            r: dict[str, Any] = {
                "half_cycle_y_m": rot[1] * c.a_t * 1e3,
                "half_cycle_vx_mm_s": rot[2] * c.a_t * c.n_t * 1e6,
                "completed_cycle": full.complete,
            }
            if full.complete:
                end = m.inertial_to_rot_890(c, t_cyc, full.state)
                r["return_m"] = float(np.linalg.norm(end[:2] - s890[:2]) * c.a_t * 1e3)
                r["return_cm_s"] = float(np.linalg.norm(end[2:] - s890[2:]) * c.a_t * c.n_t * 1e5)
            assert full.rec_t is not None and full.rec_y is not None
            tr = m.hermite_traj(full.rec_t, full.rec_y, model)
            for k in (m.I_TITANIA, m.I_OBERON):

                def bst(
                    t: float, k: int = k, model: m.ForceModel = model
                ) -> tuple[np.ndarray, np.ndarray]:
                    pp, vv = model.body_states(t)
                    return pp[k], vv[k]

                enc = m.find_minima(tr, bst, 0.0, float(full.rec_t[-1]), body=k, max_dist=30000.0)
                r[m.BODY_NAMES[k]] = [{"t_d": e.t / m.DAY_S, **e.osculating()} for e in enc]
            res[f"rtol_{rtol:g}"] = r
        out[label] = res
        for sol_name, meth, rt in (
            ("scipy_DOP853_1e-13", "DOP853", 1e-13),
            ("scipy_Radau_1e-12", "Radau", 1e-12),
            ("scipy_LSODA_1e-13", "LSODA", 1e-13),
        ):
            if label != "model_890":
                continue
            sol = m.propagate_scipy(model, 0.0, t_cyc, y0, method=meth, rtol=rt)
            end = m.inertial_to_rot_890(c, t_cyc, sol.y[:, -1])
            res[sol_name] = {
                "return_m": float(np.linalg.norm(end[:2] - s890[:2]) * c.a_t * 1e3),
                "return_cm_s": float(np.linalg.norm(end[2:] - s890[2:]) * c.a_t * c.n_t * 1e5),
            }
        log(f"{label}: {json.dumps(res, default=float)[:600]}")
    dump("p2.json", out)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("stage", choices=["control", "p2"])
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    t = time.time()
    {"control": stage_control, "p2": stage_p2}[args.stage](args)
    log(f"done in {time.time() - t:.1f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
