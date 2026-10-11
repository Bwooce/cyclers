"""#1022d: re-run of the Pluto sweeps under the stored Table 11 mean motions, with a control first.

Pre-registration: docs/notes/2026-10-11-1022d-fitted-periods-prereg.md sec. 4. Old result files are
never touched; every output goes to data/1022d_rerun/.

  control   scripts/scan_320_epoch_aware_moon_systems.py _per_system_sweep for Saturn and Neptune
            (moons the change does not touch): the rows must equal the committed
            data/scan_320_epoch_aware_{saturn,neptune}.jsonl bit for bit, apart from git_sha,
            elapsed_s and the literature-anchor fields (the literature module is replaced by an
            inert stand-in, as in run_998_enumerate.py control320)
  pluto320  the same sweep for Pluto -> data/1022d_rerun/scan_320_epoch_aware_pluto.jsonl
  pluto492  scripts/scan_492_neptune_pluto.py's Pluto scan -> data/1022d_rerun/scan_492_pluto.jsonl

Usage: uv run python scripts/run_1022d_rerun.py control [saturn|neptune]
       uv run python scripts/run_1022d_rerun.py pluto320
"""

from __future__ import annotations

import importlib.util
import json
import sys
import types
from pathlib import Path
from typing import Any

from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "data" / "1022d_rerun"
TARGETS = {
    "saturn": ("Saturn", ("Enceladus", "Tethys", "Dione", "Rhea", "Titan"), 300),
    "neptune": ("Neptune", ("Triton", "Proteus"), 30),
    "pluto": ("Pluto", ("Charon", "Hydra", "Nix"), 90),
}
IGNORED = {"git_sha", "elapsed_s", "n_lit_corpus_anchors_overlap", "lit_corpus_anchor_names"}
IGNORED |= {"is_lit_fresh_class"}


def _load_scan() -> Any:
    stub = types.ModuleType("cyclerfinder.search.literature_check")
    stub.CandidateSignature = lambda **kw: kw  # type: ignore[attr-defined]
    stub._candidate_anchors = lambda sig: []  # type: ignore[attr-defined]
    real = sys.modules.get("cyclerfinder.search.literature_check")
    sys.modules["cyclerfinder.search.literature_check"] = stub
    spec = importlib.util.spec_from_file_location(
        "scan_320_1022d", REPO / "scripts" / "scan_320_epoch_aware_moon_systems.py"
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["scan_320_1022d"] = mod
    spec.loader.exec_module(mod)
    return mod, stub, real


def _sweep(tag: str, old_path: bool = False) -> Path:
    primary, moons, n_points = TARGETS[tag]
    preflight_search(
        task_no=1022,
        region_id=f"{primary.lower()}-repeated-moon-320-sweep-1022d-rerun",
        method=MethodCapability(
            genome="#320 Vector B repeated-moon length-3 closed cycles under registry orbit_circle",
            corrector="closure residual (V_inf continuity) with phase and tof-scale grid",
            capability_tags=frozenset({"ballistic", "coplanar", "patched-conic", "circular"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=n_points,
        timing_pilot_seconds_per_point=2.0,
    )
    scan, _stub, real = _load_scan()
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / f"scan_320_epoch_aware_{tag}{'_oldpath' if old_path else ''}.jsonl"
    if old_path:
        # the pre-#1022d expression, on THIS machine: (sma_km, sqrt(mu / sma^3) * 86400)
        from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
        from cyclerfinder.search.discovery_campaign import _mean_motion_rad_day

        def _old(moon: str) -> tuple[float, float]:
            sat = SATELLITES[moon]
            return sat.sma_km, _mean_motion_rad_day(PRIMARIES[primary], sat.sma_km)

        scan.orbit_circle = _old
    try:
        scan._per_system_sweep(
            primary=primary,
            moons=moons,
            nrev_grid=(0, 1, 2, 3),
            out_path=out,
            sha="1022d-rerun",
        )
    finally:
        if real is not None:
            sys.modules["cyclerfinder.search.literature_check"] = real
        else:
            sys.modules.pop("cyclerfinder.search.literature_check", None)
    return out


def _rows(path: Path) -> list[dict[str, Any]]:
    rows = []
    for line in path.read_text().splitlines():
        d = json.loads(line)
        rows.append({k: v for k, v in d.items() if k not in IGNORED})
    return rows


def _close(x: Any, y: Any, rel: float = 1e-9) -> bool:
    if isinstance(x, float) and isinstance(y, float):
        return abs(x - y) <= rel * max(abs(x), abs(y), 1e-300)
    if isinstance(x, list) and isinstance(y, list):
        return len(x) == len(y) and all(_close(p, q, rel) for p, q in zip(x, y, strict=True))
    if isinstance(x, dict) and isinstance(y, dict):
        return x.keys() == y.keys() and all(_close(x[k], y[k], rel) for k in x)
    return bool(x == y)


def compare(tag: str, new: Path, old_path_file: Path) -> dict[str, Any]:
    """Two comparisons. (1) the real control: the new code against the pre-change expression run on
    THIS machine, bit for bit. (2) against the committed file, to 1e-9 relative (the committed file
    was written on another machine; last-digit floating-point differences and one tie between
    equal-residual cells are expected and are reported, not hidden)."""
    committed = _rows(REPO / "data" / f"scan_320_epoch_aware_{tag}.jsonl")
    a, b, c = _rows(new), _rows(old_path_file), committed
    tie_rows = [
        i
        for i, (x, y) in enumerate(zip(c, a, strict=False))
        if not _close(x, y)
        and _close(
            {k: v for k, v in x.items() if k != "phase0_deg"},
            {k: v for k, v in y.items() if k != "phase0_deg"},
        )
    ]
    return {
        "tag": tag,
        "rows": len(a),
        "new_vs_oldpath_bitwise_identical": a == b,
        "new_vs_committed_rows": len(c),
        "new_vs_committed_identical_rows": sum(1 for x, y in zip(c, a, strict=False) if x == y),
        "new_vs_committed_close_rows_1e-9": sum(
            1 for x, y in zip(c, a, strict=False) if _close(x, y)
        ),
        "new_vs_committed_phase_tie_rows": tie_rows,
    }


def main() -> None:
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "control":
        tags = sys.argv[2:] or ["neptune", "saturn"]
        res = []
        for tag in tags:
            out = _sweep(tag)
            old_out = _sweep(tag, old_path=True)
            r = compare(tag, out, old_out)
            print(json.dumps(r), flush=True)
            res.append(r)
            (OUT / f"control_{tag}.json").write_text(json.dumps(r, indent=1))
        if not all(r["new_vs_oldpath_bitwise_identical"] for r in res):
            raise SystemExit("CONTROL FAILED: an untouched system changed")
    elif cmd == "pluto492":
        from cyclerfinder.search.saturn_uranus_campaign import run_prioritized_scan

        preflight_search(
            task_no=1022,
            region_id="pluto-repeated-moon-492-scan-1022d-rerun",
            method=MethodCapability(
                genome="#285 prioritized repeated-moon scan under registry orbit_circle",
                corrector="closure residual (V_inf continuity), 12 phase samples",
                capability_tags=frozenset({"ballistic", "coplanar", "patched-conic", "circular"}),
                git_sha="working-tree",
            ),
            script_path=Path(__file__),
            n_points=216,
            timing_pilot_seconds_per_point=0.1,
        )
        OUT.mkdir(parents=True, exist_ok=True)
        summary = run_prioritized_scan(
            primary="Pluto",
            moons=["Charon", "Nix", "Hydra"],
            seq_lengths=(3,),
            n_rev_grid=(0, 1, 2, 3, 4, 5),
            out_path=OUT / "scan_492_pluto.jsonl",
            git_sha="1022d-rerun",
        )
        print(json.dumps(summary.as_dict(), default=str)[:1500], flush=True)
    elif cmd == "pluto320":
        out = _sweep("pluto")
        print(f"wrote {out}", flush=True)
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
