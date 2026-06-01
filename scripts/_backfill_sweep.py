"""Throwaway sweep: inventory legs[] -> trajectory.segments backfill readiness.

Classifies every catalogue row as migrated / legacy / neither, and ranks the
legacy rows by how much real per-leg data they carry. Used to generate the
schema-v3 backfill worklist (task #66). Not part of the package.

    uv run --with pyyaml python scripts/_backfill_sweep.py
"""

from __future__ import annotations

from pathlib import Path

import yaml

CATALOGUE = Path(__file__).resolve().parent.parent / "data" / "seed_cyclers.yaml"


def classify(row: dict) -> dict:
    """Return a per-row classification record."""
    rid = row.get("id") or "?"
    name = row.get("name") or rid
    traj = row.get("trajectory") or {}
    segments = traj.get("segments") or []
    legs = row.get("legs") or []
    gaps = row.get("data_gaps") or []

    if segments:
        category = "migrated"
        unit = segments
    elif legs:
        category = "legacy"
        unit = legs
    else:
        category = "neither"
        unit = []

    n_legs = len(unit)
    n_tof = sum(1 for leg in unit if leg.get("tof_days") is not None)
    n_ae = sum(1 for leg in unit if leg.get("a_au") is not None and leg.get("e") is not None)
    orbit = row.get("orbit_elements") or {}
    has_orbit_ae = orbit.get("a_au") is not None and orbit.get("e") is not None

    return {
        "id": rid,
        "name": name,
        "category": category,
        "n_legs": n_legs,
        "n_tof": n_tof,
        "n_ae": n_ae,
        "has_orbit_ae": has_orbit_ae,
        "n_gaps": len(gaps),
        "n_encounters": n_encounters(row),
        "model": row.get("model_assumption") or "circular-coplanar",
    }


def n_encounters(row: dict) -> int:
    """Encounter count implied by sequence_canonical (deduped closing body).

    A k-synodic cycle's segment topology has one transit leg per encounter
    (closed loop), so the number of dash-separated body codes is the number
    of legs a *complete* trajectory.segments decomposition would carry.
    """
    seq = row.get("sequence_canonical") or ""
    if not seq:
        return 0
    parts = seq.split("-")
    return len(parts)


def readiness(rec: dict) -> str:
    """Rank a legacy row: full / topology-gap / partial / sparse.

    - full:         legs cover the whole encounter topology, every leg has
                    tof_days, and a single orbit-level (a, e) applies.
                    Migrates losslessly like the Aldrin/s1l1 exemplar.
    - topology-gap: every present leg has tof_days, but legs[] under-counts
                    the encounters implied by sequence_canonical (the return /
                    intermediate-loop segments are not tabulated). Needs the
                    s1l1-style data_gaps[] markers, NOT a naive rewrite.
    - partial:      some legs lack tof_days.
    - sparse:       no legs, or no per-leg tof at all.
    """
    if rec["n_legs"] == 0:
        return "sparse"
    all_tof = rec["n_tof"] == rec["n_legs"]
    if not all_tof:
        return "partial" if rec["n_tof"] > 0 else "sparse"
    # All legs have tof. Does the leg list cover the encounter topology?
    if rec["n_legs"] >= rec["n_encounters"] and rec["has_orbit_ae"]:
        return "full"
    return "topology-gap"


def main() -> None:
    rows = yaml.safe_load(CATALOGUE.read_text())
    recs = [classify(r) for r in rows]

    by_cat: dict[str, list[dict]] = {"migrated": [], "legacy": [], "neither": []}
    for rec in recs:
        by_cat[rec["category"]].append(rec)

    print(f"total entries: {len(recs)}")
    for cat in ("migrated", "legacy", "neither"):
        print(f"  {cat}: {len(by_cat[cat])}")

    ranks = ("full", "topology-gap", "partial", "sparse")
    legacy_rank: dict[str, list[dict]] = {r: [] for r in ranks}
    for rec in by_cat["legacy"]:
        legacy_rank[readiness(rec)].append(rec)

    print("\nlegacy breakdown:")
    for rank in ranks:
        print(f"  {rank}: {len(legacy_rank[rank])}")

    print("\n=== FULL (legs cover topology, all tof, single orbit a/e) ===")
    for rec in legacy_rank["full"]:
        print(
            f"  {rec['id']:<32} legs={rec['n_legs']} enc={rec['n_encounters']} "
            f"tof={rec['n_tof']} orbit_ae={rec['has_orbit_ae']}"
        )

    print("\n=== TOPOLOGY-GAP (all tof present, legs < encounters) ===")
    print(f"  ({len(legacy_rank['topology-gap'])} entries; sample:)")
    for rec in legacy_rank["topology-gap"][:10]:
        print(
            f"  {rec['id']:<32} legs={rec['n_legs']} enc={rec['n_encounters']} "
            f"orbit_ae={rec['has_orbit_ae']}"
        )

    print("\n=== PARTIAL (some but not all legs have tof_days) ===")
    for rec in legacy_rank["partial"]:
        print(f"  {rec['id']:<32} legs={rec['n_legs']} tof={rec['n_tof']}")

    print("\n=== SPARSE (no per-leg tof or no legs) ===")
    for rec in legacy_rank["sparse"]:
        print(f"  {rec['id']:<32} legs={rec['n_legs']} tof={rec['n_tof']}")

    print("\n=== MIGRATED ===")
    for rec in by_cat["migrated"]:
        print(f"  {rec['id']:<32} segments={rec['n_legs']} gaps={rec['n_gaps']}")

    print("\n=== NEITHER (citation-only / family-seed) ===")
    for rec in by_cat["neither"]:
        print(f"  {rec['id']:<32} orbit_ae={rec['has_orbit_ae']}")


if __name__ == "__main__":
    main()
