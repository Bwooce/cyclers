"""#309 Phase 1 -- low-thrust powered cycler discovery sweep (Venus-Earth-Mars).

Companion to ``scripts/scan_309_low_thrust_em.py``. Same machinery, different
encounter tour: a closed V-E-M-V triple (Venus -> Earth -> Mars -> Venus). The
Jones-Hernandez-Jesick / Hughes-Edelman-Longuski VEM corpus (the published
ballistic VEM extensions; see ``literature_check.KNOWN_CORPUS``) covers
ballistic VEM densely; LOW-THRUST VEM is sparser in the published record, so
a clean closure here is a less-saturated candidate for the gauntlet (#274) --
emphasis on candidate.

NO catalogue writeback. NO novelty claims. Run as::

    uv run python scripts/scan_309_low_thrust_vem.py

Outputs ``data/scan_309_low_thrust_vem.jsonl`` -- one leading ``_meta`` row,
then per-candidate rows, then a trailing summary.
"""

from __future__ import annotations

import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from cyclerfinder.core.ephemeris import Ephemeris  # noqa: E402
from cyclerfinder.search.low_thrust_cycler_search import (  # noqa: E402
    LowThrustCyclerCandidate,
    datetime_to_sec_since_j2000,
    sweep_low_thrust_cyclers,
)


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT
        ).strip()
    except Exception:
        return "unknown"


# --------------------------------------------------------------------------
# Sweep grid -- V-E-M-V single-synodic.
# --------------------------------------------------------------------------

# Tour: V -> E -> M -> V. Three legs; the synodic anchor is the
# Venus/Earth pair (583.9 d synodic period); per-leg revs default to (0,0,0).
# Initial ToF shapes mirror the geometry: Venus->Earth ~150 d (short inner
# hop), Earth->Mars ~250-350 d (Hohmann-like outbound), Mars->Venus ~250-350
# d (return crossover). The optimiser then drives Lambert closure.
TOF_SHAPES_DAYS: tuple[tuple[float, float, float], ...] = (
    (150.0, 290.0, 145.0),  # short-inner / Hohmann-out / short-return
    (180.0, 250.0, 155.0),
    (120.0, 320.0, 140.0),
    (200.0, 200.0, 180.0),  # balanced
    (160.0, 280.0, 200.0),
)

# Per-leg revolution grid -- direct-only first; one multi-rev cell each.
REVS_GRID: tuple[tuple[int, int, int], ...] = (
    (0, 0, 0),
    (0, 1, 0),
    (1, 0, 0),
)

# Try a few launch-window phases through the EM synodic window so the
# optimiser is not pinned to a single phase. The circular ephemeris has no
# absolute calendar; we sample three uniformly-spaced t0 epochs around J2030.
T0_EPOCHS_SEC: tuple[float, ...] = tuple(
    datetime_to_sec_since_j2000(datetime(2030, m, 1)) for m in (1, 7)
)


def _augment_with_provenance(
    rows: list[LowThrustCyclerCandidate], sha: str
) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for cand in rows:
        payload = cand.as_dict()
        payload["_meta"] = {"git_sha": sha, "scan": "309_vem"}
        out.append(payload)
    return out


def main() -> int:
    sha = _git_sha()
    out_path = ROOT / "data" / "scan_309_low_thrust_vem.jsonl"
    print(f"[309-VEM] starting -- sha={sha} -- out={out_path}", flush=True)

    eph = Ephemeris("circular")

    meta_row = {
        "_meta": True,
        "kind": "config",
        "scan": "309_vem",
        "git_sha": sha,
        "sequence": ["V", "E", "M", "V"],
        "k_synodic": 1,
        "tof_shapes_days": [list(s) for s in TOF_SHAPES_DAYS],
        "revs_grid": [list(r) for r in REVS_GRID],
        "t0_epochs_sec": list(T0_EPOCHS_SEC),
        "synodic_pair": ["V", "E"],
        "discipline": (
            "NO catalogue writeback. NO novelty claims. "
            "Candidates feed gauntlet (V0-V5). Low-thrust VEM is sparser "
            "in the published record than ballistic VEM -- a clean closure "
            "here is a candidate for the gauntlet, not a discovery."
        ),
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as out:
        out.write(json.dumps(meta_row) + "\n")

    n_cells = len(T0_EPOCHS_SEC) * len(TOF_SHAPES_DAYS) * len(REVS_GRID)
    print(
        f"[309-VEM] sweeping {len(T0_EPOCHS_SEC)} t0 x {len(TOF_SHAPES_DAYS)} ToF "
        f"shapes x {len(REVS_GRID)} rev grids = {n_cells} cells",
        flush=True,
    )
    rows = sweep_low_thrust_cyclers(
        sequence=("V", "E", "M", "V"),
        k_synodic=1,
        ephem=eph,
        t0_epochs_sec=T0_EPOCHS_SEC,
        leg_tof_shapes_days=TOF_SHAPES_DAYS,
        per_leg_revs_grid=REVS_GRID,
        closure_body="V",  # close the tour at Venus.
        synodic_pair=("V", "E"),
        search=None,  # offline; literature_check -> inconclusive (correct default)
        n_starts=4,
        seed=0,
    )
    print(f"[309-VEM] sweep returned {len(rows)} converged candidates", flush=True)

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
            "scan": "309_vem",
            "git_sha": sha,
            "enumerated_cells": n_cells,
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
        f"[309-VEM] DONE -- enumerated={n_cells}, converged={len(rows)}, "
        f"sf_feasible={sf_feasible}, novelty_claimable_offline={novelty_claimable}, "
        f"cross_check_ok={cross_check_ok}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
