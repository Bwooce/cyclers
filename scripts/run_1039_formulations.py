"""#1039: open-chain end formulations (pre-registration:
docs/notes/2026-10-08-1039-open-chain-formulation.md).

Formulation (b), direction-only end pins, on the positive control first: the one-cycle EGGIE (#1043
seed, the paper's ideal model), reusing ``scripts/run_1043_eggie.py`` with
``run_968_rungb.END_MODE = "direction"`` and its own data directory.

    uv run python scripts/run_1039_formulations.py --form b --target eggie --stage sigma|down|ias15

Checkpoints in data/1039_<form>_<target>/; live log in .../live/ (gitignored).
"""

from __future__ import annotations

import argparse
import importlib.util
import shutil
import sys
from pathlib import Path
from typing import Any

from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search

ROOT = Path(__file__).resolve().parents[1]


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def main() -> None:
    preflight_search(
        task_no=1039,
        region_id="open-chain-end-formulation-controls",
        method=MethodCapability(
            genome="one-cycle open chains (EGGIE ideal model; GanEur#316 jup365)",
            corrector="forward-backward shooting, DOP853 + STM; sigma continuation; end forms",
            capability_tags=frozenset({"ballistic", "n-body"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=1,
    )
    ap = argparse.ArgumentParser()
    ap.add_argument("--form", choices=["b"], required=True)
    ap.add_argument("--target", choices=["eggie"], required=True)
    ap.add_argument("--stage", choices=["sigma", "down", "ias15"], required=True)
    args = ap.parse_args()
    out = ROOT / "data" / f"1039_{args.form}_{args.target}"
    out.mkdir(parents=True, exist_ok=True)
    w = _load("run_1043_eggie", ROOT / "scripts" / "run_1043_eggie.py")
    w.OUT = out
    w.Q.OUT = out
    w.Q.R.OUT = out
    w.Q.R.END_MODE = "direction"
    seed = out / "pc_chain_n1.json"
    if not seed.exists():
        shutil.copy(ROOT / "data" / "1043_eggie" / "pc_chain_n1.json", seed)
    if args.stage == "sigma":
        w.Q.R.stage_sigma(1)
    elif args.stage == "down":
        w.Q.R.stage_down(1)
    else:
        w.Q.stage_ias15(1)


if __name__ == "__main__":
    main()
