"""Shared pytest fixtures for the cyclerfinder test tree.

The live ``data/catalogue.yaml`` is parsed by 30+ test modules — many of them
several times per module (per-test ``yaml.safe_load`` re-reads). On a sharded
run that is the same multi-hundred-row YAML parsed dozens of times per worker.
The :func:`catalogue_rows` fixture parses it ONCE per worker (session scope) and
hands the result to consumers.

Consumers must treat the returned list as READ-ONLY (it is shared across the
session); the wired-in re-parsers only read, so this is behaviour-preserving.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
import yaml  # type: ignore[import-untyped]

CATALOGUE_PATH = Path(__file__).resolve().parent.parent / "data" / "catalogue.yaml"


@pytest.fixture(scope="session")
def catalogue_rows() -> list[dict[str, Any]]:
    """Parsed ``data/catalogue.yaml`` rows, parsed once per worker (READ-ONLY)."""
    rows: list[dict[str, Any]] = yaml.safe_load(CATALOGUE_PATH.read_text())
    return rows
