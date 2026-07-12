"""#576 step 1 -- source real Galilean-moon inclinations and compute MUTUAL
inclination for all 6 pairs (Io, Europa, Ganymede, Callisto).

Mandatory-first-step per the #576 Fable plan review: ``core/satellites.py``
carries NO inclination field at all (confirmed -- the dataclass carries
mu/radius/a/safe_alt only, same gap #554 found for Triton's retrograde flag).
Fable's own knowledge (high confidence, standard JPL SSD values) suggested a
favorable picture (mutual inclinations <~0.5-0.65 deg) but explicitly flagged
this as NOT repo-sourced and NOT a substitute for independently sourcing real
numbers -- exactly the #571 Titan-Iapetus lesson (an unverified/assumed
inclination premise was the actual root cause of that thread's negative, not
a method failure).

Source
------
JPL Solar System Dynamics "Planetary Satellite Mean Orbital Parameters"
table (https://ssd.jpl.nasa.gov/sats/elem/sep.html), Jupiter Galilean
satellites (JUP365 ephemeris fit), epoch 2000-01-01.5 TDB, mean orbital
elements referred to each moon's own LOCAL LAPLACE PLANE (the standard
reference frame for regular-satellite mean elements -- the same convention
`core/satellites.py`'s own sourcing note uses for the Uranian moons, and the
same frame #572 used for Iapetus's 15.5 deg). Fetched 2026-07-12. This is the
JPL SSD table's own published precision (one decimal place for both
inclination and node for these four bodies -- NOT a truncation introduced
here; the table lists no higher-precision figures for i/node for the
Galilean satellites):

    Io:       i = 0.0 deg,  node (Omega) =   0.0 deg  (P = 1.769 d)
    Europa:   i = 0.5 deg,  node (Omega) = 184.0 deg  (P = 3.551 d)
    Ganymede: i = 0.2 deg,  node (Omega) =  58.5 deg  (P = 7.155 d)
    Callisto: i = 0.3 deg,  node (Omega) = 309.1 deg  (P = 16.690 d)

(Io's node is formally undefined at i=0.0 deg -- listed as 0.0 by convention;
this does not affect the mutual-inclination formula below since sin(i)=0
kills the node-dependent term for any pair involving Io.)

Mutual inclination formula
---------------------------
Per the #576 mandate ("compute the MUTUAL inclination for each of the 6
pairs from BOTH i and the ascending node Omega, not just |i1-i2| -- check how
#571's own inclination handling did this, mirror that approach"): #571/#572
placed Iapetus on an explicit 3D circular orbit with its OWN ascending-node
longitude Omega relative to Titan's plane (`iapetus_state_3d`, R3(Omega).R1(inc)),
i.e. it treated inclination as a relative-plane-tilt geometry problem, not a
naive scalar subtraction. The Galilean analogue of that same physical
question -- "how tilted is moon B's orbital plane relative to moon A's
orbital plane" -- is the standard two-plane spherical-triangle relation
(e.g. Murray & Dermott, Solar System Dynamics, eq. 2.14):

    cos(i_mutual) = cos(i1) cos(i2) + sin(i1) sin(i2) cos(Omega1 - Omega2)

This reduces to the correct i1 (or i2) when the other body is exactly
coplanar (i=0, as Io is here) and is NOT the naive |i1 - i2| the #571 lesson
warns against -- |i1-i2| only recovers the mutual inclination in the special
case Omega1 == Omega2 (aligned nodes); for arbitrary node separation the true
mutual inclination can be LARGER than |i1-i2| (up to i1+i2), which is exactly
why the naive-subtraction shortcut is unsafe.

Run as::

    uv run python scripts/probe_576_galilean_inclination_check.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

DATA_DIR = ROOT / "data"
OUT_PATH = DATA_DIR / "probe_576_galilean_inclination_check.jsonl"

GALILEAN_MOONS: tuple[str, ...] = ("Io", "Europa", "Ganymede", "Callisto")

# JPL SSD "Planetary Satellite Mean Orbital Parameters" table
# (https://ssd.jpl.nasa.gov/sats/elem/sep.html), Jupiter Galilean satellites,
# JUP365 ephemeris fit, epoch 2000-01-01.5 TDB, mean elements referred to each
# moon's own local Laplace plane. Fetched 2026-07-12. Precision as published
# (1 decimal place for i/node -- the table's own limit for these bodies, not
# a truncation introduced here).
INCLINATION_DEG: dict[str, float] = {
    "Io": 0.0,
    "Europa": 0.5,
    "Ganymede": 0.2,
    "Callisto": 0.3,
}
NODE_OMEGA_DEG: dict[str, float] = {
    "Io": 0.0,  # formally undefined at i=0; harmless (sin(i)=0 kills the node term)
    "Europa": 184.0,
    "Ganymede": 58.5,
    "Callisto": 309.1,
}

# Fable's sanity envelope (explicitly NOT repo-sourced, a check target only,
# per the #576 dispatch instructions): mutual inclinations expected roughly
# <= 0.5-0.65 deg for any pair. Flag prominently (not silently) if exceeded.
FABLE_ENVELOPE_MAX_DEG = 0.65


def mutual_inclination_deg(body_a: str, body_b: str) -> float:
    """cos(i_mut) = cos(i1)cos(i2) + sin(i1)sin(i2)cos(Omega1-Omega2)."""
    i1 = math.radians(INCLINATION_DEG[body_a])
    i2 = math.radians(INCLINATION_DEG[body_b])
    d_omega = math.radians(NODE_OMEGA_DEG[body_a] - NODE_OMEGA_DEG[body_b])
    cos_i_mut = math.cos(i1) * math.cos(i2) + math.sin(i1) * math.sin(i2) * math.cos(d_omega)
    # Clamp for float safety (cos_i_mut can drift a hair past +-1 at i=0 pairs).
    cos_i_mut = max(-1.0, min(1.0, cos_i_mut))
    return math.degrees(math.acos(cos_i_mut))


def naive_abs_diff_deg(body_a: str, body_b: str) -> float:
    """The UNSAFE shortcut the #571 lesson warns against -- |i1 - i2|, ignoring
    node separation. Reported alongside the real mutual inclination so any
    divergence between the two is visible, not silently used."""
    return abs(INCLINATION_DEG[body_a] - INCLINATION_DEG[body_b])


def main() -> int:
    import itertools

    pairs = list(itertools.combinations(GALILEAN_MOONS, 2))
    print("[576-incl] JPL SSD Laplace-plane mean elements (epoch 2000-01-01.5 TDB, JUP365):")
    for m in GALILEAN_MOONS:
        print(
            f"[576-incl]   {m:10s} i={INCLINATION_DEG[m]:.1f} deg  "
            f"Omega={NODE_OMEGA_DEG[m]:.1f} deg"
        )

    records: list[dict[str, Any]] = []
    any_exceeds_envelope = False
    for a, b in pairs:
        i_mut = mutual_inclination_deg(a, b)
        naive = naive_abs_diff_deg(a, b)
        exceeds = i_mut > FABLE_ENVELOPE_MAX_DEG
        any_exceeds_envelope = any_exceeds_envelope or exceeds
        rec = {
            "pair": [a, b],
            "i_a_deg": INCLINATION_DEG[a],
            "i_b_deg": INCLINATION_DEG[b],
            "omega_a_deg": NODE_OMEGA_DEG[a],
            "omega_b_deg": NODE_OMEGA_DEG[b],
            "mutual_inclination_deg": i_mut,
            "naive_abs_diff_deg": naive,
            "exceeds_fable_envelope": exceeds,
        }
        records.append(rec)
        flag = "  <-- EXCEEDS Fable's ~0.65 deg sanity envelope" if exceeds else ""
        print(
            f"[576-incl] {a:10s}-{b:10s}: mutual={i_mut:.4f} deg "
            f"(naive |di|={naive:.4f} deg){flag}",
            flush=True,
        )

    max_mutual = max(r["mutual_inclination_deg"] for r in records)
    print(
        f"[576-incl] DONE: max mutual inclination across 6 pairs = {max_mutual:.4f} deg "
        f"({'EXCEEDS' if any_exceeds_envelope else 'within'} the ~0.65 deg envelope). "
        "For comparison, Titan-Iapetus (the pair that killed #571-#575) is ~15.5 deg -- "
        f"{'still' if max_mutual < 15.5 else 'NOT'} 20x+ smaller at the largest Galilean pair.",
        flush=True,
    )

    with OUT_PATH.open("w", encoding="utf-8") as fh:
        fh.write(
            json.dumps(
                {
                    "_meta": True,
                    "task": "#576 step 1 -- Galilean mutual-inclination sourcing",
                    "source": "JPL SSD Planetary Satellite Mean Orbital Parameters "
                    "(https://ssd.jpl.nasa.gov/sats/elem/sep.html), Jupiter Galilean "
                    "satellites, JUP365 ephemeris fit, epoch 2000-01-01.5 TDB, mean "
                    "elements referred to each moon's local Laplace plane. Accessed "
                    "2026-07-12.",
                    "inclination_deg": INCLINATION_DEG,
                    "node_omega_deg": NODE_OMEGA_DEG,
                    "fable_envelope_max_deg": FABLE_ENVELOPE_MAX_DEG,
                    "max_mutual_inclination_deg": max_mutual,
                    "any_exceeds_fable_envelope": any_exceeds_envelope,
                    "titan_iapetus_comparison_deg": 15.5,
                }
            )
            + "\n"
        )
        for rec in records:
            fh.write(json.dumps({"kind": "pair_result", **rec}) + "\n")
    print(f"[576-incl] written: {OUT_PATH}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
