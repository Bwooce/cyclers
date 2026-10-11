"""#1038: gauntlet_942 and scan_320 must import without loading search.literature_check.

A task that must not load that module (it is under rework) reuses the DOP853 re-fly. The check
runs in a subprocess in which the module is blocked, so any module-level import fails loudly.
"""

from __future__ import annotations

import ast
import subprocess
import sys
import textwrap
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]

_PROBE = textwrap.dedent(
    """
    import importlib.util, sys
    sys.modules["cyclerfinder.search.literature_check"] = None  # any import now raises
    for name, rel in (
        ("gauntlet_942", "scripts/gauntlet_942.py"),
        ("scan_320", "scripts/scan_320_epoch_aware_moon_systems.py"),
    ):
        spec = importlib.util.spec_from_file_location(name, rel)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
    print("OK")
    """
)


def test_scripts_import_with_literature_check_blocked() -> None:
    proc = subprocess.run(
        [sys.executable, "-c", _PROBE],
        cwd=REPO,
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert proc.stdout.strip().endswith("OK")


def test_refly_module_imports_nothing_from_literature_check() -> None:
    tree = ast.parse((REPO / "src/cyclerfinder/search/two_working_body_refly.py").read_text())
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            assert "literature_check" not in (node.module or "")
            assert all(a.name != "literature_check" for a in node.names)
        elif isinstance(node, ast.Import):
            assert all("literature_check" not in a.name for a in node.names)


def test_gauntlet_942_reexports_the_refly_functions() -> None:
    import importlib.util

    from cyclerfinder.search import two_working_body_refly as refly

    spec = importlib.util.spec_from_file_location(
        "gauntlet_942_1038", REPO / "scripts" / "gauntlet_942.py"
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert mod.cross_check is refly.cross_check
    assert mod.fly is refly.fly
