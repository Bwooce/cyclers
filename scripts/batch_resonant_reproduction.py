"""Batch V0 reproduction via resonance-anchored construction.

For every 2-body catalogue row that carries BOTH a cycler-level orbit
(a,e — directly or via perihelion/aphelion) AND sourced per-body V_inf,
construct the cycler from its sourced orbit and compare the COMPUTED
coplanar V_inf to the catalogue's sourced V_inf. A match is a V0
reproduction (our patched-conic recovers the published energy from the
published geometry).

NOT a test. Diagnostic / reproducibility artifact for the family-targeted
reproduction plan (docs/superpowers/plans/2026-06-03-family-targeted-reproduction.md).

# FINDING (2026-06-03):
#   tried=3, FULL=2, miss=1.
#   - aldrin-classic-em-k1-outbound/inbound: FULL — construction gives
#     V_E=6.58 / V_M=9.75 vs sourced 6.5 / 9.7 (coplanar reproduction).
#   - s1l1-2syn-em-cpom: "miss" only vs its STORED spec V_inf (5.65/3.05):
#     construction gives 4.90/4.98, which MATCHES the Russell 2004 coplanar
#     anchors (4.99/5.10) — see tests/search/test_resonant_construct.py. The
#     5.65/3.05 is a higher-fidelity (real-ephemeris; Mars 3.05 needs eccentric
#     Mars) figure, not the coplanar signature.
#   DATA LIMITATION: only 3 of 233 rows carry both a cycler-level (a,e) and
#   sourced V_inf. The ~200 SnLm/Russell rows have V_inf but no cycler-level
#   (a,e); rows with peri/apo (VISIT, Case1, U0L1, Hollister) lack V_inf.
#   => Catalogue-scale reproduction needs per-cycler (a,e) populated for the
#   ballistic Russell rows from the Russell 2004 dissertation tables (in
#   docs/refs/), then this batch validates them. That is the Task-4 unlock.
"""

from __future__ import annotations

import yaml  # type: ignore[import-untyped]

from cyclerfinder.search.resonant_construct import construct_resonant_cycler

_TOL_KMS = 0.5  # absorbs the documented inter-source coplanar V_inf spread


def _orbit_ae(oe: dict) -> tuple[float, float, str] | None:
    a, e = oe.get("a_au"), oe.get("e")
    if a is not None and e is not None:
        return float(a), float(e), "a/e"
    peri, apo = oe.get("perihelion_au"), oe.get("aphelion_au")
    if peri is not None and apo is not None:
        peri, apo = float(peri), float(apo)
        return (peri + apo) / 2.0, (apo - peri) / (apo + peri), "peri/apo"
    return None


def main() -> None:
    from pathlib import Path

    cat = Path(__file__).resolve().parent.parent / "data" / "catalogue.yaml"
    rows = yaml.safe_load(cat.read_text())
    tried = full = part = 0
    for r in rows:
        oe = r.get("orbit_elements") or {}
        bodies = r.get("bodies") or []
        vs = r.get("vinf_kms_at_encounters") or []
        got = _orbit_ae(oe)
        if got is None or len(bodies) != 2:
            continue
        tgt = {v["body"]: v.get("vinf_kms") for v in vs}
        if any(tgt.get(b) is None for b in bodies):
            continue
        a, e, src = got
        tried += 1
        try:
            res = construct_resonant_cycler(a_au=a, e=e, bodies=tuple(bodies))
        except Exception as exc:
            print(f"  [{src}] {r['id']}: ERR {type(exc).__name__}")
            continue
        diffs = {b: abs(res.vinf_kms[b] - float(tgt[b])) for b in bodies}
        tag = (
            "FULL"
            if all(d < _TOL_KMS for d in diffs.values())
            else ("PART" if any(d < _TOL_KMS for d in diffs.values()) else "miss")
        )
        full += tag == "FULL"
        part += tag == "PART"
        got_vs = {b: (round(res.vinf_kms[b], 2), round(float(tgt[b]), 2)) for b in bodies}
        print(f"  [{src}] {r['id']}: {tag} {got_vs}")
    print(f"tried={tried} FULL={full} PART={part}")


if __name__ == "__main__":
    main()
