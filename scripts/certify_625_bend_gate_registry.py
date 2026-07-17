"""#625 -- certify the #324 bend sub-gate over every entry sharing #610's failure mode.

Generalizes ``scripts/certify_610_proteus_bend_interval.py``'s single-entry
proof of concept (rigorous ``mpmath.iv`` interval-arithmetic certification
that the Bate-Mueller-White max-bend formula stays below the ``#324`` 5 deg
gate over a box, via ``scripts/_bend_gate_interval_cert.py``'s now-shared
:func:`certify_bend_gate_over_box`) to every OTHER ``data/empty_regions.jsonl``
entry the ``#623``/``#625`` dispatch note names as sharing the identical
undersized-moon-GM failure mode:

* ``#571``'s 4 "analytically-empty" Saturn entries (Titan-Mimas/-Enceladus/
  -Tethys/-Dione) -- 2 named in the dispatch note, but the actual registry
  has 4 (checked directly per the task's own "check once you're in the
  data" instruction).
* ``#609`` (Mars Phobos-Deimos).
* ``#607`` (Sylvia/Eugenia/Kleopatra/Elektra small-body multi-moon systems).

Two certification methods, chosen per entry's own original methodology
--------------------------------------------------------------------------
**Hohmann-floor (analytic, #571 entries)**: ``#571`` never ran a numerical
sweep -- it proved the gate analytically via the Hohmann-transfer minimum-
energy-ellipse V_inf floor (the provably lowest V_inf ANY conic connecting
two circular orbits can deliver at the departure radius -- see
``scripts/verify_571_gate_analytics.py::hohmann_vinf_at_r1``, reused here
verbatim, not re-derived). Since ``max_bend`` is monotonically DECREASING in
both r_p and V_inf, certifying at this floor V_inf (and the registry's own
safety-altitude r_p floor) bounds the ENTIRE achievable continuum, not just a
discrete grid -- the strongest available certificate.

**Data-grounded (numerical, #607/#609 entries)**: ``#607``/``#609`` DID run
real numerical sweeps (finite direct construction, not a search) -- unlike
the Hohmann-floor case, the sweep's OWN commensurate-tof/rel_offset
constraints mean the achievable V_inf at each moon is NOT free to reach the
unconstrained Hohmann floor (checked directly below: for several #607
moons, the Hohmann floor's OWN best-case bend actually EXCEEDS 5 deg -- see
"Genuine finding" below). The certifiable claim here is therefore the
#610-style Box-A one: given the REAL V_inf range realized among every
residual-sub-gate survivor (reproduced from the entries' own production
enumeration code, not copied from prose), does the bend stay under the gate
across that entire realized range. This is a narrower, but still
continuum-strength (not just the finitely-many actually-evaluated points),
claim than the Hohmann-floor one.

Genuine finding -- NOT forced through (task #625 step 8's explicit allowance)
--------------------------------------------------------------------------
Reproducing the real per-body V_inf range for #607's Sylvia system
(Romulus/Remus) shows Romulus's OWN bend reaches 5.15 deg and Remus's
reaches 7.12 deg at their own lowest realized survivor V_inf -- BOTH above
the 5 deg gate. Neither moon's bend sub-gate is universally unsatisfiable
over the real survivor range, so Sylvia's system-level entry does NOT admit
a clean single-body certification the way Proteus/#609/#571's small moons
do. This is reported honestly (``certified: false`` on both bodies, with an
explanatory note) rather than forced -- it does NOT contradict the #607
entry's own core empirical finding (0/384 Sylvia candidates pass ALL gates
simultaneously; residual+DOP853+the OTHER body's bend jointly still rule
every actual point out), only the DISPATCH NOTE'S blanket "same failure
mode" framing for this one system specifically. Eugenia, Kleopatra, and
Elektra (via ElektraGamma/ElektraDelta) DO certify cleanly this way -- see
the printed table and the written registry fields.

Writeback (task #625 step 6)
--------------------------------------------------------------------------
Adds ONE new field, ``bend_gate_certified_interval``, to each targeted
``data/empty_regions.jsonl`` line -- via direct JSON-line editing (parse the
one target line, add the key, re-dump), NOT via
:class:`cyclerfinder.data.empty_regions.EmptyRegionReport` (whose fixed
dataclass fields would silently DROP other entries' already-present extra
keys, e.g. ``reverification``/``legacy_negative_results_yaml_id``, if
round-tripped through ``_to_payload``/``_from_payload``). This mirrors how
those prior extra fields were themselves added. Schema check: confirmed
against ``tests/data/test_empty_regions.py::test_real_repo_registry_loads_whole_file``
-- the loader only reads specific known keys via ``payload["key"]``/
``payload.get(...)`` and does not reject unknown top-level keys, so no
dataclass/schema change is needed (verified, not assumed).

Run as::

    uv run --extra interval python scripts/certify_625_bend_gate_registry.py [--write]

``--write`` actually edits ``data/empty_regions.jsonl``; without it, the
script only prints the certification table (dry run, the default).
"""

from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
# Also on sys.path so `scripts.X` (dotted) resolves standalone, matching how
# tests/scripts/*.py already address these modules -- needed because a couple
# of sibling scripts (this one included) are reachable BOTH ways depending on
# caller, and mypy needs a single, consistent module identity to avoid a
# "Source file found twice under different module names" error.
sys.path.insert(0, str(ROOT))

from scan_558_uranus_all_pairs_offset_sweep import (  # noqa: E402
    GATE_RESIDUAL_KMS,
    N_REV_MAX,
    residual_at_point,
)
from verify_571_gate_analytics import hohmann_vinf_at_r1  # noqa: E402

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES  # noqa: E402
from cyclerfinder.search.physical_sanity import DEFAULT_MIN_USEFUL_BEND_DEG  # noqa: E402
from scripts._bend_gate_interval_cert import HAVE_MPMATH, certify_bend_gate_over_box  # noqa: E402
from scripts.enumerate_563_symmetric_closures import N_REV_VALUES as N_REV_563  # noqa: E402
from scripts.enumerate_563_symmetric_closures import (  # noqa: E402
    REL_OFFSETS_DEG as REL_OFFSETS_563,
)
from scripts.enumerate_563_symmetric_closures import pair_n_max  # noqa: E402
from scripts.enumerate_600_3moon_symmetric_closures import (  # noqa: E402
    REL_OFFSETS_DEG as REL_OFFSETS_600,
)
from scripts.enumerate_600_3moon_symmetric_closures import (  # noqa: E402
    residual_from_options,
    three_leg_options,
)

if HAVE_MPMATH:
    import mpmath as mp

DATA_PATH = ROOT / "data" / "empty_regions.jsonl"

N_REV_600 = tuple(range(N_REV_MAX + 1))

# Conservative, generous upper bounds for the certification box -- per
# `bend_deg_interval`'s own docstring, sup(bend) is realized at the box's
# LOWER corner (r_p floor, V_inf floor), so these upper bounds never affect
# the certified result; they're just wide enough to make each box look
# obviously conservative on inspection, not tuned to any target.
_RP_WIDEN_FACTOR = 10.0
_VINF_WIDE_PLANETARY_KMS = 50.0
_VINF_WIDE_SMALLBODY_KMS = 1.0


def _rp_safe(body: str) -> float:
    s = SATELLITES[body]
    return s.radius_eq_km + s.safe_alt_km


# --------------------------------------------------------------------------
# Data-grounded real-survivor V_inf range reproduction (#607/#609) -- mirrors
# #610's own `_proteus_subgate_vinf_range`, generalized to (a) return BOTH
# bodies' per-encounter V_inf lists (not just one), and (b) a 3-moon variant
# for Elektra. Recomputes from the entries' own production enumeration code
# (never copies a prose number), exactly like #610's precedent.
# --------------------------------------------------------------------------


def subgate_vinf_ranges_2moon(anchor: str, flyby: str, *, primary: str) -> dict[str, Any]:
    """Real per-body V_inf lists among residual-sub-gate survivors, one direction.

    One direction suffices for a 2-moon pair (spot-checked: the reverse
    direction reproduces byte-for-byte identical per-body V_inf lists, since
    swapping anchor/flyby only relabels which body is visited first, not the
    physical geometry swept).
    """
    t_syn, p_a, p_b, n_max = pair_n_max(anchor, flyby, primary=primary)
    sqrt_papb = (p_a * p_b) ** 0.5
    vinf: dict[str, list[float]] = {anchor: [], flyby: []}
    n_total = 0
    n_sub = 0
    for n in range(1, n_max + 1):
        tof_days = n * t_syn / 2.0
        tof_scale = tof_days / sqrt_papb
        for n0, n1 in itertools.product(N_REV_563, N_REV_563):
            for rel in REL_OFFSETS_563:
                n_total += 1
                pt = residual_at_point(
                    anchor,
                    flyby,
                    rel_offset_deg=rel,
                    tof_scale=tof_scale,
                    n_rev=(n0, n1),
                    primary=primary,
                )
                if pt is None or pt["residual_kms"] >= GATE_RESIDUAL_KMS:
                    continue
                n_sub += 1
                vin, vout = pt["vinf_in"], pt["vinf_out"]
                vinf[anchor].append(max(abs(vin[0]), abs(vout[0])))
                vinf[flyby].append(max(abs(vin[1]), abs(vout[1])))
                vinf[anchor].append(max(abs(vin[2]), abs(vout[2])))
    return {"n_total": n_total, "n_sub": n_sub, "vinf": vinf}


def subgate_vinf_ranges_3moon(moons: tuple[str, str, str], *, primary: str) -> dict[str, Any]:
    """Real per-moon V_inf lists among residual-sub-gate survivors, aggregated
    over ALL 6 orderings of a 3-moon chain (each moon plays a different role
    -- anchor/flyby1/flyby2 -- in different orderings, and n_sub genuinely
    differs by ordering per the #607 run, so all 6 must be swept, unlike the
    2-moon case)."""
    vinf: dict[str, list[float]] = {m: [] for m in moons}
    n_total = 0
    n_sub = 0
    for anchor, flyby1, flyby2 in itertools.permutations(moons, 3):
        t_syn_a, _, _, n_max_a = pair_n_max(anchor, flyby1, primary=primary)
        t_syn_b, _, _, n_max_b = pair_n_max(flyby1, flyby2, primary=primary)
        t_syn_c, _, _, n_max_c = pair_n_max(flyby2, anchor, primary=primary)
        for n1 in range(1, n_max_a + 1):
            tof1 = n1 * t_syn_a / 2.0
            for n2 in range(1, n_max_b + 1):
                tof2 = n2 * t_syn_b / 2.0
                for n3 in range(1, n_max_c + 1):
                    tof3 = n3 * t_syn_c / 2.0
                    for ro1, ro2 in itertools.product(REL_OFFSETS_600, REL_OFFSETS_600):
                        opts_a, opts_b, opts_c = three_leg_options(
                            anchor,
                            flyby1,
                            flyby2,
                            rel_offset1_deg=ro1,
                            rel_offset2_deg=ro2,
                            tof1_days=tof1,
                            tof2_days=tof2,
                            tof3_days=tof3,
                            primary=primary,
                        )
                        for na, nb, nc in itertools.product(N_REV_600, N_REV_600, N_REV_600):
                            n_total += 1
                            pt = residual_from_options(opts_a, opts_b, opts_c, (na, nb, nc))
                            if pt is None or pt["residual_kms"] >= GATE_RESIDUAL_KMS:
                                continue
                            n_sub += 1
                            vin, vout = pt["vinf_in"], pt["vinf_out"]
                            vs = [max(abs(vin[k]), abs(vout[k])) for k in range(4)]
                            vinf[anchor].append(vs[0])
                            vinf[flyby1].append(vs[1])
                            vinf[flyby2].append(vs[2])
                            vinf[anchor].append(vs[3])
    return {"n_total": n_total, "n_sub": n_sub, "vinf": vinf}


# --------------------------------------------------------------------------
# Per-target certification
# --------------------------------------------------------------------------


def certify_hohmann_floor(
    iv: Any, moon: str, other: str, primary: str, *, label: str
) -> dict[str, Any]:
    """#571-style: analytic Hohmann-transfer V_inf floor over the FULL
    achievable continuum (not grid-restricted) -- the strongest available
    certificate, reused verbatim from #571's own method."""
    sat = SATELLITES[moon]
    other_sat = SATELLITES[other]
    mu_primary = PRIMARIES[primary]
    vinf_floor = hohmann_vinf_at_r1(sat.sma_km, other_sat.sma_km, mu_primary)
    rp_safe = _rp_safe(moon)
    result = certify_bend_gate_over_box(
        iv,
        gm_moon_km3_s2=sat.mu_km3_s2,
        rp_lo_km=rp_safe,
        rp_hi_km=_RP_WIDEN_FACTOR * rp_safe,
        vinf_lo_kms=vinf_floor,
        vinf_hi_kms=_VINF_WIDE_PLANETARY_KMS,
        gate_deg=DEFAULT_MIN_USEFUL_BEND_DEG,
        label=label,
    )
    result["method_detail"] = "hohmann-floor (analytic, full continuum)"
    result["vinf_floor_source"] = (
        "hohmann_vinf_at_r1 (verify_571_gate_analytics.py, reused verbatim)"
    )
    return result


def certify_data_grounded(
    iv: Any, moon: str, vinf_list: list[float], *, label: str, vinf_hi_kms: float
) -> dict[str, Any]:
    """#610-Box-A-style: real survivor V_inf range from production code."""
    sat = SATELLITES[moon]
    rp_safe = _rp_safe(moon)
    vinf_min, vinf_max = min(vinf_list), max(vinf_list)
    result = certify_bend_gate_over_box(
        iv,
        gm_moon_km3_s2=sat.mu_km3_s2,
        rp_lo_km=rp_safe,
        rp_hi_km=rp_safe,
        vinf_lo_kms=vinf_min,
        vinf_hi_kms=max(vinf_hi_kms, vinf_max),
        gate_deg=DEFAULT_MIN_USEFUL_BEND_DEG,
        label=label,
    )
    result["method_detail"] = "data-grounded (real residual-sub-gate survivor V_inf range)"
    result["n_survivor_encounters"] = len(vinf_list)
    result["vinf_survivor_range_kms"] = [vinf_min, vinf_max]
    return result


def main(argv: list[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write",
        action="store_true",
        help="Actually edit data/empty_regions.jsonl (default: dry run).",
    )
    args = parser.parse_args(argv)

    if not HAVE_MPMATH:
        print(
            "[625] mpmath not installed -- run `uv run --extra interval python "
            "scripts/certify_625_bend_gate_registry.py`.",
            flush=True,
        )
        return 1

    mp.mp.dps = 50
    iv = mp.iv
    iv.dps = 50

    # region_id -> list of per-body certification results.
    by_region: dict[str, list[dict[str, Any]]] = {}

    # ---- #571: 4 Saturn Titan-pair "analytically-empty" entries ----------
    titan_pairs = [
        ("saturn-titan-mimas-analytically-empty-571", "Mimas"),
        ("saturn-titan-enceladus-analytically-empty-571", "Enceladus"),
        ("saturn-titan-tethys-analytically-empty-571", "Tethys"),
        ("saturn-titan-dione-analytically-empty-571", "Dione"),
    ]
    for region_id, moon in titan_pairs:
        res = certify_hohmann_floor(iv, moon, "Titan", "Saturn", label=moon)
        by_region.setdefault(region_id, []).append(res)
        print(
            f"[625] {region_id}: {moon} sup(bend)={res['sup_bend_deg']:.4f} deg "
            f"certified={res['certified']} ({res['method_detail']})",
            flush=True,
        )

    # ---- #609: Mars Phobos-Deimos -----------------------------------------
    mars_region = "mars-phobos-deimos-symmetric-closure-609-2026-07-16"
    mars_data = subgate_vinf_ranges_2moon("Phobos", "Deimos", primary="Mars")
    print(
        f"[625] {mars_region}: reproduced {mars_data['n_total']} evaluated, "
        f"{mars_data['n_sub']} sub-gate survivors",
        flush=True,
    )
    for moon in ("Phobos", "Deimos"):
        res = certify_data_grounded(
            iv, moon, mars_data["vinf"][moon], label=moon, vinf_hi_kms=_VINF_WIDE_PLANETARY_KMS
        )
        by_region.setdefault(mars_region, []).append(res)
        print(
            f"[625] {mars_region}: {moon} sup(bend)={res['sup_bend_deg']:.4f} deg "
            f"certified={res['certified']} ({res['method_detail']})",
            flush=True,
        )

    # ---- #607: 4 small-body multi-moon systems ----------------------------
    smallbody_region = "smallbody-multimoon-symmetric-closure-mass-limited-607-2026-07-16"

    two_moon_systems = [
        ("Sylvia", ("Romulus", "Remus")),
        ("Eugenia", ("PetitPrince", "EugeniaS2")),
        ("Kleopatra", ("AlexHelios", "CleoSelene")),
    ]
    for primary, pair in two_moon_systems:
        anchor, flyby = pair
        data = subgate_vinf_ranges_2moon(anchor, flyby, primary=primary)
        print(
            f"[625] {smallbody_region} ({primary}): reproduced {data['n_total']} evaluated, "
            f"{data['n_sub']} sub-gate survivors",
            flush=True,
        )
        for moon in pair:
            res = certify_data_grounded(
                iv,
                moon,
                data["vinf"][moon],
                label=f"{primary}:{moon}",
                vinf_hi_kms=_VINF_WIDE_SMALLBODY_KMS,
            )
            by_region.setdefault(smallbody_region, []).append(res)
            print(
                f"[625] {smallbody_region}: {primary}:{moon} "
                f"sup(bend)={res['sup_bend_deg']:.4f} deg "
                f"certified={res['certified']} ({res['method_detail']})",
                flush=True,
            )

    elektra_moons = ("ElektraBeta", "ElektraGamma", "ElektraDelta")
    elektra_data = subgate_vinf_ranges_3moon(elektra_moons, primary="Elektra")
    print(
        f"[625] {smallbody_region} (Elektra): reproduced {elektra_data['n_total']} evaluated, "
        f"{elektra_data['n_sub']} sub-gate survivors",
        flush=True,
    )
    for moon in elektra_moons:
        res = certify_data_grounded(
            iv,
            moon,
            elektra_data["vinf"][moon],
            label=f"Elektra:{moon}",
            vinf_hi_kms=_VINF_WIDE_SMALLBODY_KMS,
        )
        by_region.setdefault(smallbody_region, []).append(res)
        print(
            f"[625] {smallbody_region}: Elektra:{moon} sup(bend)={res['sup_bend_deg']:.4f} deg "
            f"certified={res['certified']} ({res['method_detail']})",
            flush=True,
        )

    # ---- Summary ------------------------------------------------------------
    n_certified = sum(1 for results in by_region.values() for r in results if r["certified"])
    n_total_bodies = sum(len(results) for results in by_region.values())
    print(
        f"[625] SUMMARY: {n_certified}/{n_total_bodies} per-body certifications succeeded "
        f"across {len(by_region)} empty_regions.jsonl entries.",
        flush=True,
    )
    for region_id, results in by_region.items():
        any_cert = any(r["certified"] for r in results)
        all_cert = all(r["certified"] for r in results)
        status = (
            "ALL bodies certified"
            if all_cert
            else ("PARTIAL -- see per-body detail" if any_cert else "NO body certified")
        )
        print(f"[625]   {region_id}: {status}", flush=True)
        for r in results:
            print(
                f"[625]     {r['label']}: sup={r['sup_bend_deg']:.4f} deg "
                f"certified={r['certified']}",
                flush=True,
            )

    if args.write:
        _write_results(by_region)
    else:
        print("[625] dry run -- pass --write to actually edit data/empty_regions.jsonl", flush=True)

    return 0


def _write_results(by_region: dict[str, list[dict[str, Any]]]) -> None:
    """Add `bend_gate_certified_interval` to each targeted line, in place.

    Direct JSON-line editing (see module docstring for why this does NOT go
    through EmptyRegionReport). Only the targeted lines are touched; every
    other line (and every other key on the targeted lines) is preserved
    byte-for-byte.
    """
    lines = DATA_PATH.read_text(encoding="utf-8").splitlines()
    touched: set[str] = set()
    out_lines: list[str] = []
    for line in lines:
        stripped = line.strip()
        if not stripped:
            out_lines.append(line)
            continue
        payload = json.loads(stripped)
        region_id = payload.get("region_id")
        if region_id in by_region:
            payload["bend_gate_certified_interval"] = {
                "task": 625,
                "generalizes": 610,
                "per_body": by_region[region_id],
                "scope_note": (
                    "Certifies only that no (r_p, V_inf) in each body's recorded box can "
                    "clear the #324 physical max-bend sub-gate. Does NOT certify the "
                    "residual/Lambert-closure half of the search is exhaustive. See "
                    "scripts/certify_625_bend_gate_registry.py and scripts/"
                    "_bend_gate_interval_cert.py."
                ),
            }
            out_lines.append(json.dumps(payload, ensure_ascii=True))
            touched.add(region_id)
        else:
            out_lines.append(line)
    missing = set(by_region) - touched
    if missing:
        raise SystemExit(f"[625] ERROR: region_id(s) not found in {DATA_PATH}: {missing}")
    DATA_PATH.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    print(
        f"[625] wrote bend_gate_certified_interval to {len(touched)} entries in {DATA_PATH}",
        flush=True,
    )


if __name__ == "__main__":
    raise SystemExit(main())
