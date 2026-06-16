"""#309 Phase 1 — low-thrust powered cycler discovery sweep (Earth-Mars).

A minimum-viable discovery probe into the LOW-THRUST regime our existing
Sims-Flanagan machinery (``src/cyclerfinder/core/sims_flanagan.py``,
``src/cyclerfinder/search/lowthrust*.py``) has covered as machinery but has
never been pointed at cycler discovery.

The driver lives in
:mod:`cyclerfinder.search.low_thrust_cycler_search` -- this script is the CLI
that sweeps it over a small (epoch x ToF-shape x revs) E-M-E grid centred on
the Aldrin phase, with the same closure-gate + literature-check + propellant
accounting the search driver wires.

NO catalogue writeback. NO novelty claims. A candidate row is a CANDIDATE for
the V0-V5 gauntlet, not a discovery -- the doc note records the verdict.

Run as::

    uv run python scripts/scan_309_low_thrust_em.py

Outputs ``data/scan_309_low_thrust_em.jsonl`` -- one leading ``_meta`` row,
then per-candidate rows, then a trailing summary.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from cyclerfinder.core.ephemeris import Ephemeris  # noqa: E402
from cyclerfinder.search.low_thrust_cycler_search import (  # noqa: E402
    LowThrustCyclerCandidate,
    sweep_low_thrust_cyclers,
)
from cyclerfinder.search.maintain import _default_t0_guess  # noqa: E402


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT
        ).strip()
    except Exception:
        return "unknown"


# --------------------------------------------------------------------------
# Sweep grid -- Aldrin-anchored E-M-E single-synodic sweep.
# --------------------------------------------------------------------------

# Aldrin reference parameters (McConaghy/Longuski/Byrnes 2002 Table 4 row
# "1L1"; catalogue ``aldrin-classic-em-k1-outbound``).
ALDRIN_EM_TOF_DAYS = 146.0
ALDRIN_ME_TOF_DAYS = 634.0
ALDRIN_EARTH_FLYBY_ALT_KM = 200.0
ALDRIN_TOF_JITTER = (20.0, 60.0)

# A small ToF-shape grid around the Aldrin centre, plus an inbound-twin
# (long outbound / short inbound). Each shape feeds the optimiser as the
# initial guess; the optimiser then drives Lambert closure + maintenance ΔV.
TOF_SHAPES_DAYS: tuple[tuple[float, float], ...] = (
    (ALDRIN_EM_TOF_DAYS, ALDRIN_ME_TOF_DAYS),  # Aldrin outbound centre
    (ALDRIN_EM_TOF_DAYS + 25.0, ALDRIN_ME_TOF_DAYS - 25.0),  # +slack outbound
    (ALDRIN_EM_TOF_DAYS - 25.0, ALDRIN_ME_TOF_DAYS + 25.0),  # -slack outbound
    (634.0, 146.0),  # inbound twin (down-escalator)
    (200.0, 580.0),  # mid-arc variant
)

# Per-leg revolution grid -- direct-only first (Aldrin's regime), 0/1/2.
REVS_GRID: tuple[tuple[int, int], ...] = (
    (0, 0),
    (0, 1),
    (1, 0),
)


def _augment_with_provenance(
    rows: list[LowThrustCyclerCandidate], sha: str
) -> list[dict[str, object]]:
    """Convert candidate rows to JSON dicts and stamp provenance."""
    out: list[dict[str, object]] = []
    for cand in rows:
        payload = cand.as_dict()
        payload["_meta"] = {"git_sha": sha, "scan": "309_em"}
        out.append(payload)
    return out


def main() -> int:
    sha = _git_sha()
    out_path = ROOT / "data" / "scan_309_low_thrust_em.jsonl"
    print(f"[309-EM] starting -- sha={sha} -- out={out_path}", flush=True)

    eph = Ephemeris("circular")
    aldrin_t0 = _default_t0_guess(ALDRIN_EM_TOF_DAYS)

    meta_row = {
        "_meta": True,
        "kind": "config",
        "scan": "309_em",
        "git_sha": sha,
        "sequence": ["E", "M", "E"],
        "k_synodic": 1,
        "tof_shapes_days": [list(s) for s in TOF_SHAPES_DAYS],
        "revs_grid": [list(r) for r in REVS_GRID],
        "aldrin_t0_guess_sec": aldrin_t0,
        "earth_flyby_alt_km": ALDRIN_EARTH_FLYBY_ALT_KM,
        "tof_jitter_half_days": list(ALDRIN_TOF_JITTER),
        "discipline": (
            "NO catalogue writeback. NO novelty claims. Candidates feed gauntlet (V0-V5)."
        ),
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as out:
        out.write(json.dumps(meta_row) + "\n")

    print(
        f"[309-EM] sweeping {len(TOF_SHAPES_DAYS)} ToF shapes x {len(REVS_GRID)} rev "
        f"grids = {len(TOF_SHAPES_DAYS) * len(REVS_GRID)} cells",
        flush=True,
    )
    rows = sweep_low_thrust_cyclers(
        sequence=("E", "M", "E"),
        k_synodic=1,
        ephem=eph,
        t0_epochs_sec=(aldrin_t0,),
        leg_tof_shapes_days=TOF_SHAPES_DAYS,
        per_leg_revs_grid=REVS_GRID,
        closure_body="E",
        closure_flyby_alt_km=ALDRIN_EARTH_FLYBY_ALT_KM,
        tof_jitter_half_days=ALDRIN_TOF_JITTER,
        synodic_pair=("E", "M"),
        search=None,  # offline; literature_check -> inconclusive (correct default)
        n_starts=4,
        seed=0,
    )
    print(f"[309-EM] sweep returned {len(rows)} converged candidates", flush=True)

    # Filter buckets used in the doc verdict.
    sf_feasible = sum(1 for r in rows if r.sims_flanagan_feasible)
    novelty_claimable = sum(1 for r in rows if r.novelty_claimable)
    cross_check_ok = sum(1 for r in rows if r.independent_cross_check_residual_kms < 1.0e-3)

    payloads = _augment_with_provenance(rows, sha)
    with out_path.open("a", encoding="utf-8") as out:
        for payload in payloads:
            out.write(json.dumps(payload) + "\n")

        summary = {
            "_meta": True,
            "kind": "summary",
            "scan": "309_em",
            "git_sha": sha,
            "enumerated_cells": len(TOF_SHAPES_DAYS) * len(REVS_GRID),
            "converged_candidates": len(rows),
            "sims_flanagan_feasible": sf_feasible,
            "novelty_claimable_offline": novelty_claimable,
            "independent_cross_check_within_1e-3_kms": cross_check_ok,
            "notes": (
                "Offline run: literature_check returns 'inconclusive' for "
                "every candidate. A real novelty pass MUST inject the live "
                "WebSearch -- see search.literature_check.check_literature."
            ),
        }
        out.write(json.dumps(summary) + "\n")

    print(
        f"[309-EM] DONE -- enumerated={len(TOF_SHAPES_DAYS) * len(REVS_GRID)}, "
        f"converged={len(rows)}, sf_feasible={sf_feasible}, "
        f"novelty_claimable_offline={novelty_claimable}, "
        f"cross_check_ok={cross_check_ok}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
