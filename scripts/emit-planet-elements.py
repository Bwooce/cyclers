"""Emit the sourced J2000 planet orbital elements as JSON for the site viz.

Usage:
    uv run python scripts/emit-planet-elements.py            # write data/planet-elements.json
    uv run python scripts/emit-planet-elements.py --check    # verify on-disk file matches

Single source of truth is :mod:`cyclerfinder.core.constants` (PLANETS) plus the
sourced inclined elements in :mod:`cyclerfinder.core.ephemeris`
(``inclined_planets``). This script is the only place the two are joined into the
6-element osculating set ``(a, e, i, Ω, ϖ, L0)`` the time-true visualization on
cyclers.space consumes. The site never hand-copies these numbers; it fetches the
emitted ``data/planet-elements.json`` (mirroring the catalogue sync), so every
value traces back to constants.py and cannot drift.

All angles are degrees, J2000 ecliptic frame; ``a`` is AU. Citation:
Standish & Williams, "Approximate Positions of the Planets", JPL Solar System
Dynamics, Table 1 (1800-2050 AD) — the same table that supplies every element.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from cyclerfinder.core.constants import PLANETS
from cyclerfinder.core.ephemeris import inclined_planets

_CITATION = "J2000 osculating ellipse (Standish & Williams Table 1)"


def build() -> dict[str, object]:
    """Assemble the planet-elements payload from the upstream single source."""
    inclined = inclined_planets()  # real i/Ω for inclined bodies; Earth stays 0.0
    bodies: list[dict[str, object]] = []
    for code, p in PLANETS.items():
        inc = inclined[code]
        bodies.append(
            {
                "code": code,
                "name": p.name,
                "a_au": p.sma_au,
                "e": p.ecc,
                "i_deg": inc.inc_deg,
                "lan_deg": inc.lan_deg,
                "varpi_deg": p.varpi_deg,
                "L0_deg": p.L0_deg,
                "mean_motion_deg_day": p.mean_motion_deg_day,
            }
        )
    return {
        "epoch": "J2000",
        "frame": "ecliptic",
        "citation": _CITATION,
        "bodies": bodies,
    }


def _serialize(payload: dict[str, object]) -> str:
    return json.dumps(payload, indent=2) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else "")
    parser.add_argument(
        "--check", action="store_true", help="verify on-disk file is up to date (exit 1 if stale)"
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parent.parent
    out_path = repo_root / "data" / "planet-elements.json"
    text = _serialize(build())

    if args.check:
        if not out_path.exists():
            sys.stderr.write(f"emit-planet-elements: {out_path} missing — run without --check.\n")
            raise SystemExit(1)
        if out_path.read_text() != text:
            sys.stderr.write(f"emit-planet-elements: {out_path} is stale — re-run to regenerate.\n")
            raise SystemExit(1)
        sys.stderr.write(f"emit-planet-elements: {out_path} up to date.\n")
        return

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text)
    sys.stderr.write(
        f"emit-planet-elements: wrote {out_path} ({len(text)} bytes, {len(PLANETS)} bodies).\n"
    )


if __name__ == "__main__":
    main()
