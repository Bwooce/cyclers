"""N-body harness Phase A: package + optional-dep wiring (plan Phase A)."""

from __future__ import annotations

import importlib

import pytest


def test_nbody_package_imports() -> None:
    mod = importlib.import_module("cyclerfinder.nbody")
    assert mod is not None


def test_rebound_is_in_validation_extra() -> None:
    import tomllib
    from pathlib import Path

    root = Path(__file__).resolve().parents[2]
    data = tomllib.loads((root / "pyproject.toml").read_text())
    extras = data["project"]["optional-dependencies"]
    assert any(d.startswith("rebound") for d in extras["validation"]), (
        "rebound must join the existing 'validation' extra (design Q2: shared-"
        "DE440 cross-check class), not a new extra"
    )


def test_rebound_skips_cleanly_when_absent() -> None:
    """The fast suite must not hard-require rebound (mirror spiceypy skip)."""
    rebound = pytest.importorskip("rebound")
    assert rebound is not None
