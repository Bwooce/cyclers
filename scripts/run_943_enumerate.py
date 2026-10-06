"""#943 (X1) entry point of the two-working-body enumeration: the Jovian cells.

Cells gc (Ganymede-Callisto), ge (Ganymede-Europa) and the one-body recall cells gc1, ge1.
Same driver, flags and outputs as scripts/run_942_enumerate.py (which holds the code);
only the preflight task number differs, so the task-number guard sees #943.

Usage: uv run python scripts/run_943_enumerate.py --cell gc --k 1,2,3 ... (see run_942)
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
from typing import Any

from cyclerfinder.data.preflight import preflight_search

_DRIVER = Path(__file__).resolve().parent / "run_942_enumerate.py"


def _preflight_943(**kwargs: Any) -> None:
    preflight_search(task_no=943, script_path=Path(__file__), **kwargs)


def main() -> None:
    spec = importlib.util.spec_from_file_location("run_942_enumerate", _DRIVER)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.main(preflight=_preflight_943, x1=True)


if __name__ == "__main__":
    main()
