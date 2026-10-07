"""#998: Pluto-Charon one-working-node cyclers with a passive small moon.

Charon is the only massive flyby body; Styx, Nix, Kerberos or Hydra is a massless target (the
Russell & Strange 2007/2009 one-body architecture). The generator, the seeds and the gate are the
#942/#943 ones (scripts/run_942_enumerate.py, search/two_working_body*.py); this driver only adds
the four Pluto cells and the #320 recall control.

Ideal model: circular, coplanar, primary mass = Pluto+Charon system GM (975.5 km^3/s^2, the
registry PRIMARIES["Pluto"]), each body on a circle whose radius follows from its period by
Kepler III (the R-S convention). Periods are the mean sidereal periods fitted to plu060.bsp
(3-year arc from 2030-01-01, positions relative to the Pluto-system barycentre); see
docs/notes/2026-10-07-998-pluto-small-moon-cyclers.md sec. 1.

Cells (a = Charon):
  ps  Styx massless      pn  Nix massless      pk  Kerberos massless      ph  Hydra massless

Subcommands:
  enumerate ARGS     the run_942 driver with these cells and this task's preflight
  controlgen --out F VenMar#45 through the generator (code-path control)
  control320 --out F re-run the #320 Pluto sweep (scripts/scan_320_epoch_aware_moon_systems.py,
                     its own code path, literature step stubbed out) into F and compare with the
                     committed data/scan_320_epoch_aware_pluto.jsonl

Usage:
  uv run python scripts/run_998_enumerate.py enumerate --cell pn --k 1,2 --out DIR \
      --max-returns 3,0 --transfer-revs 0,1 --generic-revs 1,2 [--shard i/n] [--count-only]
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
import types
from pathlib import Path
from typing import Any

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES, _sat
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.search.two_working_body import CircularSystem, FlybyBody

REPO = Path(__file__).resolve().parents[1]
DAY = 86400.0

#: Mean sidereal periods (days) fitted to plu060.bsp, J2000 frame, relative to the Pluto-system
#: barycentre (NAIF 9), 3000 samples over 2030-01-01 + 3 yr. Literature (Showalter & Hamilton
#: 2015; Brozovic et al. 2015) agrees to better than 2e-4 d.
PERIODS_D = {
    "Charon": 6.38722,
    "Styx": 20.16195,
    "Nix": 24.85472,
    "Kerberos": 32.16798,
    "Hydra": 38.20202,
}
#: Target constants (km, km^3/s^2): Nix and Hydra from the registry (satellites.py); Styx and
#: Kerberos are approximate (not in the registry). The targets are massless (no turn), so these
#: enter only the SOI diagnostic; non-zero GM keeps that diagnostic from collapsing to 0.
TARGET_GM = {"Styx": 0.0007, "Nix": 0.0015, "Kerberos": 0.0011, "Hydra": 0.0020}
RADII_KM = {"Styx": 5.0, "Nix": 18.0, "Kerberos": 6.0, "Hydra": 18.5}
CHARON_GM, CHARON_RADIUS, CHARON_FLOOR = 106.1, 606.0, 100.0  # registry satellites.py
CELLS = {"ps": "Styx", "pn": "Nix", "pk": "Kerberos", "ph": "Hydra"}


def pluto_system(moon: str) -> CircularSystem:
    mu = PRIMARIES["Pluto"]
    if moon not in SATELLITES:
        # Styx and Kerberos are not in the registry; the SOI diagnostic of the generator looks
        # every encounter body up there. Add them to the in-process dict only (no file edit),
        # with the fitted Kepler a and the approximate target constants above.
        per = PERIODS_D[moon] * DAY
        a_fit = (mu * (per / (2.0 * math.pi)) ** 2) ** (1.0 / 3.0)
        SATELLITES[moon] = _sat(moon, "Pluto", TARGET_GM[moon], RADII_KM[moon], a_fit, 10.0)
    bodies = {}
    for m in ("Charon", moon):
        per = PERIODS_D[m] * DAY
        a = (mu * (per / (2.0 * math.pi)) ** 2) ** (1.0 / 3.0)
        bodies[m] = (a, per, 0.0)
    over = {
        "Charon": FlybyBody("Charon", CHARON_GM, CHARON_RADIUS, CHARON_FLOOR),
        moon: FlybyBody(moon, TARGET_GM[moon], RADII_KM[moon], 10.0),
    }
    return CircularSystem(mu, bodies, frozenset({moon}), flyby_overrides=over)


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _preflight_998(**kwargs: Any) -> None:
    preflight_search(task_no=998, script_path=Path(__file__), **kwargs)


def cmd_enumerate(argv: list[str]) -> None:
    enum = _load("run_942_enumerate", REPO / "scripts" / "run_942_enumerate.py")
    original = enum.cell_system

    def cell_system(cell: str) -> tuple[CircularSystem, str, str]:
        if cell in CELLS:
            return pluto_system(CELLS[cell]), "Charon", CELLS[cell]
        result: tuple[CircularSystem, str, str] = original(cell)
        return result

    enum.cell_system = cell_system  # main() looks the cell up through the module global
    sys.argv = [sys.argv[0], *argv]
    enum.main(preflight=_preflight_998, x1=False)


def cmd_control320(argv: list[str]) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    # The scan imports search.literature_check at module level; this task does not use that module
    # (the literature step is deferred to #972), so the import is satisfied by an inert stand-in.
    # Only the residual, V_inf and physical-gate fields are compared; the anchor-overlap field is
    # not reproduced.
    stub = types.ModuleType("cyclerfinder.search.literature_check")
    stub.CandidateSignature = lambda **kw: kw  # type: ignore[attr-defined]
    stub._candidate_anchors = lambda sig: []  # type: ignore[attr-defined]
    real = sys.modules.get("cyclerfinder.search.literature_check")
    sys.modules["cyclerfinder.search.literature_check"] = stub
    try:
        scan = _load("scan_320", REPO / "scripts" / "scan_320_epoch_aware_moon_systems.py")
    finally:
        if real is not None:
            sys.modules["cyclerfinder.search.literature_check"] = real
        else:
            sys.modules.pop("cyclerfinder.search.literature_check", None)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    _preflight_998(
        region_id="pluto-charon-hydra-nix-320-recall-control",
        method=_method("#320 Vector B repeated-moon length-3 closed cycles, recall"),
        n_points=90,
        timing_pilot_seconds_per_point=2.0,
    )
    scan._per_system_sweep(
        primary="Pluto",
        moons=("Charon", "Hydra", "Nix"),
        nrev_grid=(0, 1, 2, 3),
        out_path=args.out,
        sha="998-recall",
    )
    old = _rows(REPO / "data" / "scan_320_epoch_aware_pluto.jsonl")
    new = _rows(args.out)
    report = []
    for key, o in old.items():
        n = new.get(key)
        if n is None:
            report.append({"key": list(key), "status": "MISSING in re-run"})
            continue
        report.append(
            {
                "key": list(key),
                "old_residual": o["residual_kms"],
                "new_residual": n["residual_kms"],
                "old_vinf": o["vinf_per_encounter_kms"],
                "new_vinf": n["vinf_per_encounter_kms"],
                "old_physical": o["physical_gate_passed"],
                "new_physical": n["physical_gate_passed"],
                "old_silver": o["is_silver"],
                "new_silver": n["is_silver"],
                "max_abs_diff": max(
                    abs(o["residual_kms"] - n["residual_kms"]),
                    *(
                        abs(a - b)
                        for a, b in zip(
                            o["vinf_per_encounter_kms"],
                            n["vinf_per_encounter_kms"],
                            strict=True,
                        )
                    ),
                ),
            }
        )
    extra = [list(k) for k in new if k not in old]
    out = args.out.with_suffix(".compare.json")
    out.write_text(json.dumps({"rows": report, "extra_in_rerun": extra}, indent=1))
    print(f"compared {len(old)} committed rows; wrote {out}", flush=True)


def cmd_controlgen(argv: list[str]) -> None:
    """Code-path control of the generator in this checkout: VenMar#45 (R-S 2007, Venus flyby,
    Mars massless) must come out of the vm cell as an exact zero at the stored V_inf."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    from cyclerfinder.search.two_working_body_enum import assess, solve_structure

    enum = _load("run_942_enumerate", REPO / "scripts" / "run_942_enumerate.py")
    system, a, b = enum.cell_system("vm")
    key = "k2|LV>M/0s|LM>V/0s"
    _, cyc = enum.parse_cycle_key(key, system, a, b)
    zeros = solve_structure(
        system,
        cyc,
        phase_period_s=system.synodic_s(a, b),
        n_phase=36,
        n_split=12,
        n_refine=40,
    )
    stored = {"V": 8.219573173685143, "M": 12.96475955374864}  # data/942_cell_vm_gauntlet.json
    out = []
    for z in zeros:
        ass = assess(system, z)
        out.append(
            {
                "residual_kms": z.residual_kms,
                "status": ass.status,
                "vinf_kms": ass.vinf_kms,
                "published_R-S_2007_Table_3": {"V": 8.22, "M": 12.96, "period_d": 667.8},
                "stored_942": stored,
            }
        )
    args.out.write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps(out, default=float), flush=True)


def _rows(path: Path) -> dict[tuple[Any, ...], dict[str, Any]]:
    rows = {}
    for line in path.read_text().splitlines():
        d = json.loads(line)
        if d.get("_meta"):
            continue
        rows[(tuple(d["sequence"]), tuple(d["n_rev"]))] = d
    return rows


def _method(genome: str) -> Any:
    from cyclerfinder.data.method_capability import MethodCapability

    return MethodCapability(
        genome=genome,
        corrector="closure residual (V_inf continuity) with phase and tof-scale grid",
        capability_tags=frozenset({"ballistic", "coplanar", "patched-conic", "circular"}),
        git_sha="working-tree",
    )


def main() -> None:
    if len(sys.argv) < 2 or sys.argv[1] not in ("enumerate", "control320", "controlgen"):
        raise SystemExit(__doc__)
    cmd, argv = sys.argv[1], sys.argv[2:]
    commands = {
        "enumerate": cmd_enumerate,
        "control320": cmd_control320,
        "controlgen": cmd_controlgen,
    }
    commands[cmd](argv)


if __name__ == "__main__":
    main()
