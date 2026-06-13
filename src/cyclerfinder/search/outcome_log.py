"""Passive solver-outcome logger (#210) — costless training-data capture.

Every time a solver decides an outcome (converged or not, cost, validation-relevant
scalars), :func:`log_outcome` appends ONE ``(inputs -> outcome)`` JSON line to the
path named by the environment variable ``CYCLERFINDER_OUTCOME_LOG``. If that variable
is unset, the call is a NO-OP — zero cost, no file, normal runs entirely unaffected.

The accumulated JSONL is the ``(genome -> outcome)`` corpus a future ANN cycler-search
surrogate (Ozaki et al. 2022, arXiv:2111.11858, whose architecture needs ~7e6 samples)
would later train on. Capturing it now is a zero-cost byproduct of normal solves.

HARD BOUNDARY — this is TRAINING-DATA CAPTURE ONLY.
--------------------------------------------------
The log is NEVER read back to validate a cycler. An ANN / surrogate output can never be
a golden: the V0-V5 validation gauntlet still governs every catalogue claim. This module
only *writes*; nothing in the validation path ever *reads* it. Keeping that boundary
explicit is the whole point of putting the warning here.

Design constraints (#210)
-------------------------
* Opt-in: a NO-OP unless ``CYCLERFINDER_OUTCOME_LOG`` is set.
* Near-zero cost when off (one env lookup, then return).
* Append-only JSONL; the parent directory is created on demand.
* NEVER raises into the caller — capturing training data must never break a solve.
  Any failure is swallowed (logged to the stdlib ``logging`` channel and ignored).
* Records are flat and JSON-serializable; numpy scalars/arrays are coerced to
  float / list via :func:`to_jsonable`.
* Each record carries a monotonic counter (deterministic-friendly) rather than
  depending on wall-clock; a wall-clock stamp is included only as a convenience
  side field.

Pure-ish: depends only on the stdlib (os / json / logging / time / itertools) plus
numpy for the coercion helper.
"""

from __future__ import annotations

import itertools
import json
import logging
import os
import time
from pathlib import Path
from typing import Any

import numpy as np

__all__ = [
    "ENV_VAR",
    "SCHEMA_VERSION",
    "enable_default_outcome_log",
    "log_outcome",
    "to_jsonable",
]

# Bump when the record envelope (the keys this module stamps) changes shape.
SCHEMA_VERSION = 1

# The single switch: set this env var to a writable path to turn capture on.
ENV_VAR = "CYCLERFINDER_OUTCOME_LOG"

# Gitignored default location used by campaign DRIVER scripts (not by tests or
# library imports) when they opt in via enable_default_outcome_log(). ``out/`` is
# already in .gitignore, so the accruing corpus is never committed.
_DEFAULT_LOG_DIR = "out/outcome_log"

_log = logging.getLogger(__name__)


def enable_default_outcome_log(tag: str) -> str | None:
    """Opt a DRIVER SCRIPT into outcome capture without anyone remembering the env var.

    If ``CYCLERFINDER_OUTCOME_LOG`` is already set (CI, a manual export, a test),
    leave it untouched and return it. Otherwise point it at the gitignored
    ``out/outcome_log/<tag>.jsonl`` so running a real campaign accrues the corpus
    automatically. Call this ONLY from a campaign driver's ``__main__`` — never on
    library import (that would silently turn capture on everywhere, including tests).
    Returns the active log path.
    """
    existing = os.environ.get(ENV_VAR)
    if existing:
        return existing
    path = f"{_DEFAULT_LOG_DIR}/{tag}.jsonl"
    os.environ[ENV_VAR] = path
    return path


# Process-local monotonic record counter (deterministic-friendly: it does not
# depend on wall-clock and is reproducible within a single process run).
_counter = itertools.count()


def to_jsonable(obj: Any) -> Any:
    """Coerce ``obj`` into a JSON-serializable structure.

    Handles the types our solvers emit: numpy scalars -> python float/int/bool,
    numpy arrays -> nested lists, dict / list / tuple recursed elementwise, and
    plain JSON scalars passed through. Anything else falls back to ``str(obj)``
    so a stray non-serializable value can never break the append.
    """
    if obj is None or isinstance(obj, (bool, int, float, str)):
        return obj
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, dict):
        return {str(k): to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [to_jsonable(v) for v in obj]
    # Last resort: stringify (e.g. enums, dataclasses, complex) so we never raise.
    try:
        if isinstance(obj, complex):
            return {"real": obj.real, "imag": obj.imag}
        return str(obj)
    except Exception:  # pragma: no cover - defensive
        return repr(obj)


def log_outcome(
    *,
    solver: str,
    inputs: dict[str, Any],
    outcome: dict[str, Any],
    meta: dict[str, Any] | None = None,
) -> None:
    """Append ONE ``(inputs -> outcome)`` record to the outcome log.

    NO-OP unless the environment variable named by :data:`ENV_VAR`
    (``CYCLERFINDER_OUTCOME_LOG``) is set to a writable path. When set, one JSON
    object is appended (newline-terminated) to that file; the parent directory is
    created if needed.

    Parameters
    ----------
    solver:
        Name of the solver / call-site producing this outcome (e.g.
        ``"cr3bp.correct_periodic"``).
    inputs:
        The decision inputs / genome for this solve (coerced to JSON).
    outcome:
        The decided outcome scalars (converged flag, cost, residuals, …).
    meta:
        Optional free-form context (caller-supplied, deterministic-friendly).

    This function NEVER raises into the caller: any I/O or serialization error is
    swallowed and logged to the stdlib ``logging`` channel. Capturing training
    data must never break a solve.
    """
    path = os.environ.get(ENV_VAR)
    if not path:
        # Disabled: zero-cost path. No file, no error, no side effect.
        return
    try:
        record = {
            "schema_version": SCHEMA_VERSION,
            "counter": next(_counter),
            "wall_time": time.time(),  # convenience side field, NOT the ordering key
            "solver": str(solver),
            "inputs": to_jsonable(inputs),
            "outcome": to_jsonable(outcome),
            "meta": to_jsonable(meta) if meta is not None else None,
        }
        line = json.dumps(record, separators=(",", ":"), sort_keys=True)
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("a", encoding="utf-8") as fh:
            fh.write(line + "\n")
    except Exception as exc:  # never propagate into the caller
        _log.warning("outcome_log: failed to append record (%s): %s", solver, exc)
