"""#888 -- demanded-turn gate: positive controls, regression, and re-screen.

Stages (each writes JSONL under ``data/found/888_turn_gate/``)::

    uv run python scripts/screen_888_turn_gate_rescreen.py --stage controls
    uv run python scripts/screen_888_turn_gate_rescreen.py --stage regression
    uv run python scripts/screen_888_turn_gate_rescreen.py --stage rescreen
    uv run python scripts/screen_888_turn_gate_rescreen.py --stage ungated
    uv run python scripts/screen_888_turn_gate_rescreen.py --stage extended

controls    McConaghy-Longuski-Byrnes 2002 Table 4 (19 Earth-Mars nPr cyclers) and
            Russell-Strange 2009 Tables 3-6 (10 generic-leg moon cyclers) through
            the gate, with the published numbers alongside.
regression  The six withdrawn Uranian (1,1) rows (data/withdrawn/).
rescreen    Every stored gate-passing symmetric closure of the two-moon
            enumerations (#563 Uranus, #576 Jupiter, #575 and #655 Saturn, #599
            Neptune, #609 Mars): rebuilt, magnitudes checked against the stored
            record, then gated at the project floor and at 50 km.
ungated     The same enumerations re-run over each file's own ranges (n_rev 0-3,
            rel_offset 0/180, tof = n T_syn / 2 up to its tof_scale bound) with
            ONLY the closure residual gate (0.05 km/s) and the demanded-turn gate:
            the #324 capacity gate, which pruned before, is dropped. Miranda is
            added at Uranus.
extended    Wider ranges in the same ideal model: n_rev 0-6 per leg, both
            multi-revolution Lambert branches per leg, tof_scale up to 6.

Closure definition and branch rule:
:func:`cyclerfinder.verify.turn_gate_closures.symmetric_closure`.
"""

from __future__ import annotations

import argparse
import datetime as dt
import itertools
import json
import math
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import yaml  # type: ignore[import-untyped]  # noqa: E402

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES  # noqa: E402
from cyclerfinder.search.discovery_campaign import _mean_motion_rad_day  # noqa: E402
from cyclerfinder.verify.turn_gate import demanded_turn_gate  # noqa: E402
from cyclerfinder.verify.turn_gate_closures import (  # noqa: E402
    RS_GENERIC_CYCLERS,
    RS_TITAN_MIN_ALT_KM,
    SymmetricClosure,
    mcconaghy_npr_cycler,
    russell_strange_generic_cycler,
    symmetric_closure,
)

OUT_DIR = ROOT / "data" / "found" / "888_turn_gate"
GATE_RESIDUAL_KMS = 0.05  # the #558 closure gate, unchanged
UNIFORM_FLOOR_KM = 50.0

ENUMERATIONS: tuple[str, ...] = (
    "enumerate_563_symmetric_closures.jsonl",
    "enumerate_576_jupiter_galilean_symmetric_closures.jsonl",
    "enumerate_575_titan_iapetus_symmetric_closures.jsonl",
    "enumerate_655_saturn_rhea_titan_symmetric_closures.jsonl",
    "enumerate_655_saturn_dione_rhea_symmetric_closures.jsonl",
    "enumerate_655_saturn_enceladus_tethys_symmetric_closures.jsonl",
    "enumerate_655_saturn_tethys_dione_symmetric_closures.jsonl",
    "enumerate_599_neptune_triton_proteus_symmetric_closures.jsonl",
    "enumerate_609_mars_phobos_deimos_symmetric_closures.jsonl",
)

#: McConaghy, Longuski and Byrnes (AIAA 2002-4420) Table 4, p.6:
#: (n, P, r, aphelion AU, V_inf Earth km/s, required turn deg, max turn deg).
MCCONAGHY_TABLE4: tuple[tuple[int, str, int, float, float, float, float], ...] = (
    (1, "L", 1, 2.23, 6.54, 84, 72),
    (2, "L", 2, 2.33, 10.06, 134, 44),
    (2, "L", 3, 1.51, 5.65, 135, 82),
    (3, "L", 4, 1.89, 11.78, 167, 35),
    (3, "L", 5, 1.45, 7.61, 167, 62),
    (3, "S", 5, 1.52, 12.27, 167, 33),
    (4, "S", 5, 1.82, 11.23, 167, 38),
    (4, "S", 6, 1.53, 8.51, 167, 54),
    (5, "S", 4, 2.49, 10.62, 134, 41),
    (5, "S", 5, 2.09, 9.08, 134, 50),
    (5, "S", 6, 1.79, 7.51, 135, 62),
    (5, "S", 7, 1.54, 5.86, 135, 79),
    (5, "S", 8, 1.34, 4.11, 136, 103),
    (6, "S", 4, 2.81, 7.93, 83, 59),
    (6, "S", 5, 2.37, 6.94, 84, 68),
    (6, "S", 6, 2.04, 5.96, 84, 78),
    (6, "S", 7, 1.78, 4.99, 85, 90),
    (6, "S", 8, 1.57, 4.02, 85, 104),
    (6, "S", 9, 1.40, 3.04, 86, 120),
)
MCCONAGHY_BALLISTIC = {"6S7", "6S8", "6S9"}  # Table 4 footnote e


def log(msg: str) -> None:
    print(f"[{dt.datetime.now().isoformat(timespec='seconds')}] {msg}", flush=True)


def git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT
        ).strip()
    except Exception:
        return "unknown"


def write_jsonl(name: str, meta: dict[str, Any], rows: list[dict[str, Any]]) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / name
    with path.open("w", encoding="utf-8") as fh:
        fh.write(json.dumps({"_meta": True, "git_sha": git_sha(), **meta}) + "\n")
        for r in rows:
            fh.write(json.dumps(r) + "\n")
    log(f"wrote {path.relative_to(ROOT)} ({len(rows)} rows)")
    return path


# ---------------------------------------------------------------------------
# controls
# ---------------------------------------------------------------------------


def stage_controls() -> None:
    rows: list[dict[str, Any]] = []
    for n, p, r, ra, vinf, req, mx in MCCONAGHY_TABLE4:
        c = mcconaghy_npr_cycler(n, p, r)  # type: ignore[arg-type]
        assert c is not None, (n, p, r)
        rep = demanded_turn_gate([c.encounter])
        e = rep.encounters[0]
        name = f"{n}{p}{r}"
        rows.append(
            {
                "kind": "mcconaghy_2002_table4",
                "cycler": name,
                "published": {
                    "aphelion_au": ra,
                    "vinf_earth_kms": vinf,
                    "required_turn_deg": req,
                    "max_turn_deg": mx,
                    "ballistic": name in MCCONAGHY_BALLISTIC,
                },
                "rebuilt": {"aphelion_au": c.aphelion_au, "vinf_earth_kms": c.vinf_earth_kms},
                "gate": rep.as_dict(),
                "verdict_agrees": rep.turn_feasible == (name in MCCONAGHY_BALLISTIC),
            }
        )
    for row in RS_GENERIC_CYCLERS:
        floor = RS_TITAN_MIN_ALT_KM if row.flyby_body == "Titan" else 0.0
        rb = russell_strange_generic_cycler(row, alt_floor_km=floor)
        rep = demanded_turn_gate(rb.encounters)
        proj_floor = SATELLITES[row.flyby_body].safe_alt_km
        rep_proj = demanded_turn_gate(rb.encounters, alt_floor_km=proj_floor)
        rows.append(
            {
                "kind": "russell_strange_2009",
                "id": row.rs_id,
                "published": {
                    "vinf_kms": row.vinf_flyby_kms,
                    "period_days": row.period_days,
                    "min_flyby_alt_km": row.min_flyby_alt_km,
                    "min_dist_primary_km": row.min_dist_primary_km,
                    "max_dist_primary_km": row.max_dist_primary_km,
                },
                "rebuilt": {
                    "leg_vinf_kms": list(rb.leg_vinf_kms),
                    "period_days": rb.period_days,
                    "min_dist_primary_km": min(rb.leg_rp_km),
                    "max_dist_primary_km": max(rb.leg_ra_km),
                    "min_required_alt_km": rep.min_required_alt_km,
                    "leg_branch": list(rb.leg_branch),
                },
                "gate_at_source_floor_km": floor,
                "gate": rep.as_dict(),
                "gate_at_project_floor_km": proj_floor,
                "gate_project_floor": rep_proj.as_dict(),
            }
        )
    write_jsonl("controls.jsonl", {"task": "#888 controls"}, rows)
    for r in rows:
        if r["kind"] == "mcconaghy_2002_table4":
            e = r["gate"]["encounters"][0]
            pub = r["published"]
            log(
                f"  {r['cycler']:4s} turn {e['demanded_turn_deg']:6.2f} "
                f"(pub {pub['required_turn_deg']}) "
                f"max {e['available_bend_deg']:6.2f} (pub {pub['max_turn_deg']}) "
                f"feasible={r['gate']['turn_feasible']} agrees={r['verdict_agrees']}"
            )
        else:
            proj_ok = r["gate_project_floor"]["turn_feasible"]
            log(
                f"  {r['id']:10s} min req alt {r['rebuilt']['min_required_alt_km']:8.0f} "
                f"(pub {r['published']['min_flyby_alt_km']}) "
                f"feasible@src={r['gate']['turn_feasible']} "
                f"feasible@project({r['gate_at_project_floor_km']} km)={proj_ok}"
            )


# ---------------------------------------------------------------------------
# regression
# ---------------------------------------------------------------------------


def rebuild_withdrawn(path: Path) -> tuple[dict[str, Any], SymmetricClosure, float]:
    """Rebuild a withdrawn row; rel_offset and n_rev chosen by V-inf match to the row.

    The rows do not carry rel_offset as a field, and the #312 row's ``legs``
    block lists n_revs 0 while its name and #885's reproduction say (1, 1); so
    rel_offset in {0, 180} and n_rev in {0, 1}^2 are tried and the combination
    whose rebuilt magnitudes reproduce ``vinf_kms_at_encounters`` is kept.
    """
    row = yaml.safe_load(path.read_text())[0]
    seq = row["sequence_canonical"].split("-")
    tof = float(row["legs"][0]["tof_days"])
    want = [float(e["vinf_kms"]) for e in row["vinf_kms_at_encounters"]]
    best: tuple[float, SymmetricClosure] | None = None
    for rel in (0.0, 180.0):
        for nr in ((0, 0), (1, 1), (0, 1), (1, 0)):
            c = symmetric_closure(
                "Uranus", seq[0], seq[1], tof_days=tof, n_rev=nr, rel_offset_deg=rel
            )
            if c is None:
                continue
            err = max(abs(a - b) for a, b in zip(c.stored_convention_vinf, want, strict=True))
            if best is None or err < best[0]:
                best = (err, c)
    assert best is not None
    return row, best[1], best[0]


def stage_regression() -> None:
    rows: list[dict[str, Any]] = []
    for path in sorted((ROOT / "data" / "withdrawn").glob("*-1-1-uranian-quasi-cycler-2026.yaml")):
        row, c, err = rebuild_withdrawn(path)
        rep = demanded_turn_gate(c.encounters)
        rep50 = demanded_turn_gate(c.encounters, alt_floor_km=UNIFORM_FLOOR_KM)
        rows.append(
            {
                "kind": "withdrawn_row",
                "id": row["id"],
                "rel_offset_deg": c.rel_offset_deg,
                "n_rev": list(c.n_rev),
                "tof_days": c.tof_days,
                "vinf_reproduction_err_kms": err,
                "gate_project_floor": rep.as_dict(),
                "gate_50km": rep50.as_dict(),
                "wrap_local_check_deg": c.wrap_local_check_deg,
            }
        )
        log(
            f"  {row['id']}: err {err:.1e} "
            + " ".join(
                f"{e.body} {e.demanded_turn_deg:.1f}/{e.available_bend_deg:.1f} ({e.ratio:.2f}x)"
                for e in rep.encounters
            )
        )
    write_jsonl("regression_withdrawn.jsonl", {"task": "#888 regression"}, rows)


# ---------------------------------------------------------------------------
# helpers shared by rescreen / ungated / extended
# ---------------------------------------------------------------------------


def gate_closure(c: SymmetricClosure) -> dict[str, Any]:
    rep = demanded_turn_gate(c.encounters)
    rep50 = demanded_turn_gate(c.encounters, alt_floor_km=UNIFORM_FLOOR_KM)
    return {
        "primary": c.primary,
        "anchor": c.anchor,
        "flyby": c.flyby,
        "tof_days": c.tof_days,
        "n_rev": list(c.n_rev),
        "rel_offset_deg": c.rel_offset_deg,
        "residual_kms": c.residual_kms,
        "vinf_kms": list(c.stored_convention_vinf),
        "wrap_local_check_deg": c.wrap_local_check_deg,
        "gate_project_floor": rep.as_dict(),
        "gate_50km": rep50.as_dict(),
        "pass_project_floor": rep.turn_feasible,
        "pass_50km": rep50.turn_feasible,
        "worst_ratio_project_floor": rep.worst_ratio,
    }


def physical_key(g: dict[str, Any]) -> tuple[Any, ...]:
    """Same periodic chain seen from either end: moon set, leg time, magnitudes."""
    a, b = g["anchor"], g["flyby"]
    va, vb = round(g["vinf_kms"][0], 6), round(g["vinf_kms"][1], 6)
    return (g["primary"], round(g["tof_days"], 6), tuple(sorted([(a, va), (b, vb)])))


def summarise(name: str, gated: list[dict[str, Any]]) -> dict[str, Any]:
    keys_all = {physical_key(g) for g in gated}
    keys_proj = {physical_key(g) for g in gated if g["pass_project_floor"]}
    keys_50 = {physical_key(g) for g in gated if g["pass_50km"]}
    s = {
        "kind": "summary",
        "set": name,
        "records": len(gated),
        "distinct_chains": len(keys_all),
        "records_pass_project_floor": sum(1 for g in gated if g["pass_project_floor"]),
        "records_pass_50km": sum(1 for g in gated if g["pass_50km"]),
        "distinct_pass_project_floor": len(keys_proj),
        "distinct_pass_50km": len(keys_50),
        "min_worst_ratio_project_floor": min(
            (g["worst_ratio_project_floor"] for g in gated), default=None
        ),
    }
    log(
        f"  {name}: {s['records']} records ({s['distinct_chains']} distinct); "
        f"pass@floor {s['records_pass_project_floor']} "
        f"({s['distinct_pass_project_floor']} distinct); "
        f"pass@50km {s['records_pass_50km']} ({s['distinct_pass_50km']}); "
        f"min worst ratio {s['min_worst_ratio_project_floor']}"
    )
    return s


# ---------------------------------------------------------------------------
# rescreen
# ---------------------------------------------------------------------------


def stage_rescreen() -> None:
    rows: list[dict[str, Any]] = []
    for fname in ENUMERATIONS:
        lines = [
            json.loads(x) for x in (ROOT / "data" / fname).read_text().splitlines() if x.strip()
        ]
        meta = lines[0]
        primary = meta["primary"]
        passes = [x for x in lines if x.get("kind") == "pass"]
        gated: list[dict[str, Any]] = []
        for p in passes:
            c = symmetric_closure(
                primary,
                p["anchor"],
                p["flyby"],
                tof_days=p["tof_days"],
                n_rev=(p["n_rev"][0], p["n_rev"][1]),
                rel_offset_deg=p["rel_offset_deg"],
            )
            assert c is not None, p
            err = max(
                abs(a - b)
                for a, b in zip(c.stored_convention_vinf, p["vinf_per_encounter_kms"], strict=True)
            )
            g = gate_closure(c)
            g["source_file"] = fname
            g["stored_vinf_kms"] = p["vinf_per_encounter_kms"]
            g["vinf_reproduction_err_kms"] = err
            gated.append(g)
        for g in gated:
            rows.append({"kind": "closure", **g})
        if gated:
            worst_err = max(g["vinf_reproduction_err_kms"] for g in gated)
            log(f"  {fname}: worst stored-vs-rebuilt V-inf error {worst_err:.2e} km/s")
        rows.append({**summarise(fname, gated), "primary": primary})
    write_jsonl("rescreen_stored.jsonl", {"task": "#888 rescreen of stored passes"}, rows)


# ---------------------------------------------------------------------------
# ungated / extended enumeration
# ---------------------------------------------------------------------------


def synodic_days(primary: str, a: str, b: str) -> float:
    mu = PRIMARIES[primary]
    sa, sb = SATELLITES[a], SATELLITES[b]
    pa = 2 * math.pi / _mean_motion_rad_day(mu, sa.sma_km)
    pb = 2 * math.pi / _mean_motion_rad_day(mu, sb.sma_km)
    if sa.retrograde != sb.retrograde:
        return 1.0 / (1.0 / pa + 1.0 / pb)
    return 1.0 / abs(1.0 / pa - 1.0 / pb)


def enumerate_pair(
    primary: str,
    anchor: str,
    flyby: str,
    *,
    tof_scale_max: float,
    n_rev_max: int,
    both_branches: bool,
) -> tuple[int, list[dict[str, Any]]]:
    mu = PRIMARIES[primary]
    pa = 2 * math.pi / _mean_motion_rad_day(mu, SATELLITES[anchor].sma_km)
    pb = 2 * math.pi / _mean_motion_rad_day(mu, SATELLITES[flyby].sma_km)
    t_syn = synodic_days(primary, anchor, flyby)
    n_max = math.floor(2.0 * tof_scale_max * math.sqrt(pa * pb) / t_syn)
    evaluated = 0
    out: list[dict[str, Any]] = []
    for n in range(1, n_max + 1):
        tof = n * t_syn / 2.0
        for n0, n1 in itertools.product(range(n_rev_max + 1), repeat=2):
            opts0 = ["low", "high"] if (both_branches and n0 > 0) else ["closest"]
            opts1 = ["low", "high"] if (both_branches and n1 > 0) else ["closest"]
            for rel in (0.0, 180.0):
                for b0, b1 in itertools.product(opts0, opts1):
                    br = (b0, b1, b0)
                    evaluated += 1
                    c = symmetric_closure(
                        primary,
                        anchor,
                        flyby,
                        tof_days=tof,
                        n_rev=(n0, n1),
                        rel_offset_deg=rel,
                        branches=br,  # type: ignore[arg-type]
                    )
                    if c is None or c.residual_kms >= GATE_RESIDUAL_KMS:
                        continue
                    g = gate_closure(c)
                    g["n_commensurate_int"] = n
                    g["branches"] = list(br)
                    out.append(g)
    return evaluated, out


def _pair_job(args: tuple[str, str, str, float, int, bool]) -> tuple[int, list[dict[str, Any]]]:
    primary, a, b, tsm, n_rev_max, both = args
    return enumerate_pair(primary, a, b, tof_scale_max=tsm, n_rev_max=n_rev_max, both_branches=both)


def run_enumeration(
    stage: str,
    *,
    n_rev_max: int,
    scale: float | None,
    both_branches: bool,
    workers: int = 4,
    only_primary: str | None = None,
) -> None:
    systems: list[tuple[str, list[str], float, str]] = []
    for fname in ENUMERATIONS:
        meta = json.loads((ROOT / "data" / fname).read_text().splitlines()[0])
        moons = list(meta["moons"])
        if meta["primary"] == "Uranus":
            moons = ["Miranda", *moons]
        if meta["primary"] == "Neptune":
            # The stored #599 construction pairs a retrograde (clockwise) Triton
            # position with a prograde circular velocity (discovery_campaign.
            # _moon_state); its V-infinity vectors are not physical. Skipped.
            continue
        if only_primary is not None and meta["primary"] != only_primary:
            continue
        systems.append((meta["primary"], moons, float(meta["tof_scale_max_bound"]), fname))
    rows: list[dict[str, Any]] = []
    t0 = time.time()
    total_eval = 0
    for primary, moons, tsm, fname in systems:
        tof_scale_max = tsm if scale is None else scale
        gated_all: list[dict[str, Any]] = []
        pairs = list(itertools.permutations(moons, 2))
        with ProcessPoolExecutor(max_workers=workers) as pool:
            results = pool.map(
                _pair_job,
                [(primary, a, b, tof_scale_max, n_rev_max, both_branches) for a, b in pairs],
            )
            for (a, b), (ev, gated) in zip(pairs, results, strict=True):
                total_eval += ev
                gated_all.extend(gated)
                npass = sum(1 for g in gated if g["pass_project_floor"])
                log(
                    f"  {primary} {a}-{b}-{a}: evaluated {ev}, "
                    f"residual-gate closures {len(gated)}, "
                    f"turn-gate passes (floor) {npass} [{time.time() - t0:.0f}s]"
                )
        for g in gated_all:
            g["source_meta"] = fname
            rows.append({"kind": "closure", **g})
        rows.append(
            {
                **summarise(f"{stage}:{fname}", gated_all),
                "primary": primary,
                "moons": moons,
                "tof_scale_max": tof_scale_max,
                "n_rev_max": n_rev_max,
                "both_branches": both_branches,
            }
        )
    suffix = "" if only_primary is None else f"_{only_primary.lower()}"
    write_jsonl(
        f"{stage}{suffix}.jsonl",
        {
            "task": f"#888 {stage} enumeration",
            "residual_gate_kms": GATE_RESIDUAL_KMS,
            "n_rev_max": n_rev_max,
            "tof_scale_max_override": scale,
            "both_branches": both_branches,
            "only_primary": only_primary,
            "rel_offsets_deg": [0.0, 180.0],
            "total_evaluated": total_eval,
            "elapsed_s": time.time() - t0,
        },
        rows,
    )


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--stage",
        required=True,
        choices=["controls", "regression", "rescreen", "ungated", "extended"],
    )
    ap.add_argument("--extended-scale", type=float, default=6.0)
    ap.add_argument("--extended-nrev", type=int, default=6)
    ap.add_argument("--workers", type=int, default=4, help="at most 4 on the shared machine")
    ap.add_argument(
        "--primary", type=str, default=None, help="restrict ungated/extended to one primary"
    )
    args = ap.parse_args(argv)
    log(f"stage {args.stage}: start")
    if args.stage == "controls":
        stage_controls()
    elif args.stage == "regression":
        stage_regression()
    elif args.stage == "rescreen":
        stage_rescreen()
    elif args.stage == "ungated":
        run_enumeration(
            "ungated",
            n_rev_max=3,
            scale=None,
            both_branches=False,
            workers=min(args.workers, 4),
            only_primary=args.primary,
        )
    else:
        run_enumeration(
            "extended",
            n_rev_max=args.extended_nrev,
            scale=args.extended_scale,
            both_branches=True,
            workers=min(args.workers, 4),
            only_primary=args.primary,
        )
    log(f"stage {args.stage}: done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
