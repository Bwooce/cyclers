"""Hollister & Menning 1970 Table 3 as cycles for the two-working-body corrector (#942).

Source: W. M. Hollister and M. D. Menning, "Periodic Swing-By Orbits between
Earth and Venus", J. Spacecraft and Rockets 7(10), 1970, pp. 1193-1198; Table 3
(pp. 1196-1198). The transcription is ``data/sources/hollister-menning-1970-table3.yaml``.

Every orbit is 5 blocks of [E, E, V, V, V] plus the closing E (16 yr, 5844 d,
inclined-elliptic model; p.1194). The E pair is a full-revolution return
(dates 365/366 d apart) or a symmetric return (about 490 d); each of the two
Venus returns is full-revolution (225 d) or symmetric (about 330 d).

Read errors and print errors
----------------------------
* The YAML was corrected to the print on 2026-10-05 (30 cells, including
  orbit 1 row 12 planet E and orbit 6 row 3 date 995; see
  ``docs/notes/2026-10-05-hollister-menning-1970-table3-recheck.md`` and the
  Menning 1968 thesis). Values are read as stored.
* PRINT ERRORS in the paper (the YAML matches the print): dates that break the
  225-d Venus full-revolution step, namely orbit 2 "5715" (5590 + 225 = 5815),
  orbit 4 "2573" (2358 + 225 = 2583), orbit 5 "5662" (5877 - 225 = 5652),
  orbit 6 "2585, 2810" after 2135 (a two-full-revolution block; 2810 falls after
  the next Earth date 2716), orbit 8 "4870" (5645 + 225 = 5870). The block TYPE
  is unambiguous in each case; dates reached through full-revolution legs are
  computed from the block start, so these printed dates are never used.
* Orbit 15's date span is 5843 d; kept as printed.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import yaml  # type: ignore[import-untyped]

from cyclerfinder.core.constants import SECONDS_PER_DAY
from cyclerfinder.search.two_working_body import (
    Cycle,
    LambertLeg,
    Leg,
    ResonantLeg,
    System,
    date_residual,
)

TABLE3_PATH = (
    Path(__file__).resolve().parents[3] / "data" / "sources" / "hollister-menning-1970-table3.yaml"
)

#: Paper footnote: "dates repeating after 16 yr (add 5844 days)".
PERIOD_DAYS = 5844.0

#: (orbit, block) -> Venus return types where a print error hides the type
#: from the date gaps (see module docstring).
_VENUS_TYPE_OVERRIDE: dict[tuple[int, int], str] = {(2, 4): "FS", (6, 1): "FF", (8, 4): "FF"}


#: Print errors inferred from the block structure (the page shows the left value).
#: Orbit 5 block 1 prints V 977, 1242, 1467: a 265-d first gap is too short for a
#: symmetric return at the printed V_r (the Lambert residual at the printed dates
#: is 8.9 km/s) and the turn pattern (34.0, 34.0, 38.7) is a two-full-revolution
#: block's; 1017 = 1242 - 225.
#: The other entries restore the 225-d Venus full-revolution step at the print
#: errors listed in the module docstring (they affect only date comparisons).
_PRINT_ERROR_FIXES: dict[tuple[int, int], tuple[str, float, str, float]] = {
    (5, 2): ("V", 977.0, "V", 1017.0),
    (2, 23): ("V", 5715.0, "V", 5815.0),
    (4, 9): ("V", 2573.0, "V", 2583.0),
    (5, 22): ("V", 5662.0, "V", 5652.0),
    (6, 8): ("V", 2585.0, "V", 2360.0),
    (6, 9): ("V", 2810.0, "V", 2585.0),
    (8, 23): ("V", 4870.0, "V", 5870.0),
}


@dataclass(frozen=True)
class Row:
    planet: str
    date: float
    vr_emos: float
    theta_deg: float
    rmin_radii: float


def load_table3(path: Path = TABLE3_PATH) -> dict[int, list[Row]]:
    """Table 3 rows per orbit, with the inferred print-error date fixes applied."""
    data = yaml.safe_load(path.read_text())
    out: dict[int, list[Row]] = {}
    for o in data["orbits"]:
        rows = [
            Row(
                e["planet"],
                float(e["date_j2440000"]),
                float(e["vr_emos"]),
                float(e["theta_deg"]),
                float(e["rmin_radii"]),
            )
            for e in o["encounters"]
        ]
        for (orb, i), (pl_old, d_old, pl_new, d_new) in _PRINT_ERROR_FIXES.items():
            if o["orbit"] == orb and (rows[i].planet, rows[i].date) == (pl_old, d_old):
                r = rows[i]
                rows[i] = Row(pl_new, d_new, r.vr_emos, r.theta_deg, r.rmin_radii)
        out[int(o["orbit"])] = rows
    return out


def _venus_gap_type(orbit: int, block: int, gap: float) -> str:
    """``F`` for a full-revolution gap (225 d; 215 d occurs as a print error),
    ``S`` for a symmetric return (between one and two Venus periods)."""
    if 205.0 <= gap <= 235.0:
        return "F"
    if 240.0 <= gap <= 450.0:
        return "S"
    raise ValueError(
        f"orbit {orbit} block {block}: Venus gap {gap} d is a print error (override it)"
    )


def block_types(orbit: int, rows: list[Row]) -> list[tuple[str, str]]:
    """``(earth_type, venus_types)`` per block, e.g. ``("FR", "FS")``."""
    out = []
    for k in range(5):
        i = 5 * k
        eg = rows[i + 1].date - rows[i].date
        et = "FR" if eg < 400 else "SY"
        if (orbit, k) in _VENUS_TYPE_OVERRIDE:
            vt = _VENUS_TYPE_OVERRIDE[(orbit, k)]
        else:
            g1 = rows[i + 3].date - rows[i + 2].date
            g2 = rows[i + 4].date - rows[i + 3].date
            vt = "".join(_venus_gap_type(orbit, k, g) for g in (g1, g2))
        out.append((et, vt))
    return out


def build_cycle(
    system: System, orbit: int, rows: list[Row], sy_branches: tuple[str, ...]
) -> tuple[Cycle, np.ndarray]:
    """The orbit's leg chain and the printed-date seed (seconds) of each Lambert leg start.

    ``sy_branches`` gives the multi-revolution branch (``"low"``/``"high"``) of
    each symmetric return in chain order. Dates reached through full-revolution
    legs are the block start plus whole periods (never a printed date).
    """
    legs: list[Leg] = []
    starts: list[float] = []
    pe, pv = system.period_s("E"), system.period_s("V")
    day = SECONDS_PER_DAY
    sy = iter(sy_branches)
    for k, (et, vt) in enumerate(block_types(orbit, rows)):
        i = 5 * k
        t = rows[i].date * day
        if et == "FR":
            legs.append(ResonantLeg("E"))
            t += pe
        else:
            legs.append(LambertLeg("E", "E", 1, next(sy)))
            starts.append(t)
            t = rows[i + 1].date * day
        legs.append(LambertLeg("E", "V"))
        starts.append(t)
        t = rows[i + 2].date * day
        for j, c in enumerate(vt):
            if c == "F":
                legs.append(ResonantLeg("V"))
                t += pv
            else:
                legs.append(LambertLeg("V", "V", 1, next(sy)))
                starts.append(t)
                t = rows[i + 3 + j].date * day
        legs.append(LambertLeg("V", "E"))
        starts.append(t)
    return Cycle(tuple(legs), PERIOD_DAYS * day), np.asarray(starts)


def n_symmetric(orbit: int, rows: list[Row]) -> int:
    return sum((et == "SY") + vt.count("S") for et, vt in block_types(orbit, rows))


def pick_branches(system: System, orbit: int, rows: list[Row]) -> tuple[tuple[str, ...], float]:
    """Pre-registered branch rule: the symmetric-return branch combination that
    minimises the sum of squared junction residuals AT THE PRINTED DATES."""
    best: tuple[tuple[str, ...], float] | None = None
    for combo in itertools.product(("low", "high"), repeat=n_symmetric(orbit, rows)):
        cyc, x = build_cycle(system, orbit, rows, combo)
        r = date_residual(system, cyc, x)
        if r is None:
            continue
        s = float(np.sum(r**2))
        if best is None or s < best[1]:
            best = (combo, s)
    if best is None:
        raise ValueError(
            f"orbit {orbit}: no branch combination is Lambert-feasible at the printed dates"
        )
    return best
