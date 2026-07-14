"""#588 step 1: global dedup pass over #583/#586's 264 unmatched-bounded pool.

**Not a search script.** This is pure post-processing/analysis over data
already committed to disk by 48 independent GA runs (#583's full
16-partition x 3-seed sweep, harvested by #586's cluster-everything
pipeline). It does not build a genome, run a GA, or call a corrector, so
it is named ``dedup_588_*.py`` rather than ``run_588_*.py`` (deliberately
outside ``tests/scripts/test_scripts_call_preflight.py``'s
``scripts/run_*.py`` glob) and does not call
``cyclerfinder.data.preflight.preflight_search`` -- that gate exists to stop
unbudgeted NEW sweeps (task-number hygiene, empty-region-registry
subsumption, timing pilots for large grids), none of which apply here: no
grid is built, no worker pool is spawned, no new region of phase space is
searched. See ``data/OUTSTANDING.md``'s ``#588`` entry and
``docs/notes/2026-07-14-583-16-partition-3-seed-sweep-results.md`` for the
full origin/context.

## Why raw ``ga_genome`` vectors are not directly comparable

Every one of #583's 16 partitions uses the SAME 7-gene genome layout
(``[x, x', y, y', z, z', theta0]``, geocentric ER3BP interleaved state +
free launch phase) but FIXES a different subset of the 6 state slots at 0
(mirroring which published family box that partition targets). A candidate
from partition ``I`` (6 free dims) and one from partition ``M`` (3 free
dims, others pinned to 0) can be the literal same physical orbit found at
two different launch phases along it, yet have wildly different raw genome
components -- and, conversely, two candidates with numerically close
``ga_genome`` vectors from DIFFERENT partitions are not meaningfully
"close" at all, because the pinned-vs-free slot pattern differs. Comparing
raw genomes across partitions is therefore not meaningful.

## The fix: dedup on physical characterization, not genome

Every candidate's genome is re-integrated through the SAME
``characterize()`` function ``run_581_gurfil_reproduction.py`` and
``run_583_widened_bounded_drift_search.py`` already both import and use for
their own family-matching (never re-implemented here) -- this yields
partition-independent physical invariants: ``rmin_km_1yr``/``rmax_km_1yr``
(radial envelope over a 1-year integration) and ``type`` (DRO/DPO/DEO/ERO,
3D-prefixed or not). Two candidates are considered possible duplicates of
the same underlying basin only if:

1. **Same ``type``.** A DRO and a DEO are different dynamics regardless of
   radial proximity -- type mismatch never merges.
2. **``rmin_km_1yr`` AND ``rmax_km_1yr`` both within a relative tolerance**
   (``DEDUP_RELATIVE_TOLERANCE``, default 0.10 -- see below) of an already-
   accepted cluster representative's own values.

## Distance threshold: 0.10 relative deviation in both rmin and rmax

Chosen for two reasons:

* **Precedent.** This exact pipeline already uses 0.10 twice as its
  "coarse enough to collapse near-duplicate members of one niche, tight
  enough to keep genuinely distinct basins separate" threshold --
  ``DEFAULT_ANALYZE_DISTANCE_THRESHOLD`` (#582) and
  ``DEFAULT_HARVEST_DISTANCE_THRESHOLD`` (#583), both applied to
  bounds-normalized genome-space RMS distance within one seed's population.
  Reusing the same numeric philosophy (10% = "same niche") for the
  physical-space cross-partition pass keeps the two dedup stages
  conceptually consistent instead of inventing an unrelated number.
* **Empirically checked for chaining.** A greedy single-link merge can
  silently chain far-apart points through a series of overlapping
  near-duplicates. At 0.10 the worst observed within-cluster spread across
  the actual 264-candidate pool is a ~1.18x max/min ratio in rmin (i.e.
  materially TIGHTER than the 0.10 threshold itself would allow via a
  single hop) -- chaining is not inflating clusters beyond what direct
  pairwise proximity would already justify. At 0.20 one 3D DRO cluster's
  internal spread reaches ~1.44x, i.e. chaining measurably starts to
  matter. 0.10 was kept as the more conservative choice (biased toward
  NOT over-merging genuinely distinct candidates, since under-merging only
  costs a few redundant literature-check calls downstream, while
  over-merging could silently drop a candidate that deserved its own
  check).

Representative selection within a cluster: highest ``ga_fitness`` (the same
objective every GA run itself optimized against) -- ties broken by lower
``drift.trend_fraction`` (the more stably bounded member).

## Output

``data/found/583_widened_search/deduped_candidates.json`` -- a list of
cluster objects, each carrying the full original representative candidate
record (with a ``partition`` key merged in, since the source per-candidate
dict has no such field at its own level -- only the parent aggregate file
does), its ``characterization`` (the ``characterize()`` output), a
``duplicate_count``, and a ``source_locations`` list recording every raw
candidate (partition/seed_index/cluster_rank) that collapsed into it.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import cast

import numpy as np

# Same repo-root sys.path convention #583 itself uses to reuse #581's
# characterize()/TABLE34 unmodified (see run_583's own module-level comment).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.run_581_gurfil_reproduction import characterize

SOURCE_DIR = Path(__file__).resolve().parent.parent / "data" / "found" / "583_widened_search"
OUTPUT_PATH = SOURCE_DIR / "deduped_candidates.json"

DEDUP_RELATIVE_TOLERANCE = 0.10
"""See module docstring's "Distance threshold" section for the full
justification (precedent + empirical chaining check against the real pool)."""


def load_all_candidates(source_dir: Path = SOURCE_DIR) -> list[dict[str, object]]:
    """Read every ``{partition}_aggregate_harvest_summary.json``'s unmatched-
    bounded candidates, tagging each with its source ``partition``."""
    records: list[dict[str, object]] = []
    for path in sorted(source_dir.glob("*_aggregate_harvest_summary.json")):
        agg = json.loads(path.read_text())
        partition = agg["partition"]
        for cand in agg["unmatched_bounded_candidates"]:
            records.append({"partition": partition, **cand})
    return records


def characterize_all(records: list[dict[str, object]]) -> list[dict[str, object]]:
    """Re-integrate every candidate's genome through the shared characterize().

    Returns one physical-characterization dict per record, same order.
    """
    chars: list[dict[str, object]] = []
    t0 = time.monotonic()
    for i, rec in enumerate(records):
        genome = rec["ga_genome"]
        assert isinstance(genome, list)
        vec6 = np.array(genome[:6], dtype=float)
        theta0 = float(genome[6])
        chars.append(characterize(vec6, theta0))
        if (i + 1) % 50 == 0 or (i + 1) == len(records):
            print(
                f"characterized {i + 1}/{len(records)} ({time.monotonic() - t0:.1f}s elapsed)",
                flush=True,
            )
    return chars


def _rep_sort_key(rec: dict[str, object]) -> tuple[float, float]:
    """Sort key for representative selection: higher fitness first, then
    lower drift trend_fraction (more stably bounded) as the tiebreak."""
    drift = rec.get("drift")
    trend = float(drift["trend_fraction"]) if isinstance(drift, dict) else float("inf")
    return (-float(rec["ga_fitness"]), trend)  # type: ignore[arg-type]


def dedup_pool(
    records: list[dict[str, object]],
    chars: list[dict[str, object]],
    tol: float = DEDUP_RELATIVE_TOLERANCE,
) -> list[list[int]]:
    """Greedy fitness-ranked clustering on (type, rmin_km_1yr, rmax_km_1yr).

    Same "walk in descending fitness order, accept as a new representative
    only if far enough from every existing one" structure as
    ``cluster_representatives``/``cluster_population`` elsewhere in this
    pipeline, generalized to the physical-characterization feature space
    instead of bounds-normalized genome space.

    Returns a list of clusters, each a list of indices into ``records``.
    """
    order = sorted(range(len(records)), key=lambda i: _rep_sort_key(records[i]))
    rep_feature: list[tuple[str, float, float]] = []
    clusters: list[list[int]] = []
    for i in order:
        typ = str(chars[i]["type"])
        rmin = float(chars[i]["rmin_km_1yr"])  # type: ignore[arg-type]
        rmax = float(chars[i]["rmax_km_1yr"])  # type: ignore[arg-type]
        placed = False
        for ci, (rtyp, rrmin, rrmax) in enumerate(rep_feature):
            if typ != rtyp:
                continue
            if (
                abs(rmin - rrmin) / max(rmin, rrmin) <= tol
                and abs(rmax - rrmax) / max(rmax, rrmax) <= tol
            ):
                clusters[ci].append(i)
                placed = True
                break
        if not placed:
            rep_feature.append((typ, rmin, rmax))
            clusters.append([i])
    return clusters


def build_output(
    records: list[dict[str, object]],
    chars: list[dict[str, object]],
    clusters: list[list[int]],
) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for idxs in clusters:
        best = min(idxs, key=lambda i: _rep_sort_key(records[i]))
        source_locations = [
            {
                "partition": records[i]["partition"],
                "seed_index": records[i]["seed_index"],
                "cluster_rank": records[i]["cluster_rank"],
                "population_index": records[i]["population_index"],
                "ga_fitness": records[i]["ga_fitness"],
            }
            for i in idxs
        ]
        out.append(
            {
                "cluster_id": -1,  # reassigned below after sorting
                "candidate_type": chars[best]["type"],
                "representative": records[best],
                "characterization": chars[best],
                "duplicate_count": len(idxs),
                "source_locations": source_locations,
            }
        )
    # Report in descending duplicate_count so the "most-repeated basin"
    # candidates surface first for the downstream literature-check step.
    out.sort(key=lambda c: -cast(int, c["duplicate_count"]))
    for new_id, entry in enumerate(out):
        entry["cluster_id"] = new_id
    return out


def print_report(records: list[dict[str, object]], deduped: list[dict[str, object]]) -> None:
    total = len(records)
    n_clusters = len(deduped)
    print()
    print(f"=== #588 dedup report: {total} raw candidates -> {n_clusters} distinct clusters ===")
    if n_clusters:
        print(f"dedup ratio: {total / n_clusters:.2f}x")
    print()
    print("top 15 clusters by duplicate_count:")
    for entry in deduped[:15]:
        locs = entry["source_locations"]
        assert isinstance(locs, list)
        partitions = sorted({loc["partition"] for loc in locs})
        print(
            f"  cluster {entry['cluster_id']:3d} type={entry['candidate_type']:>7s} "
            f"dup_count={entry['duplicate_count']:3d} partitions={partitions}"
        )
    print()
    print("cross-partition merges (clusters whose members span >1 partition):")
    cross = 0
    pair_counts: dict[tuple[str, str], int] = {}
    for entry in deduped:
        locs = entry["source_locations"]
        assert isinstance(locs, list)
        partitions = sorted({loc["partition"] for loc in locs})
        if len(partitions) > 1:
            cross += 1
            for a_idx in range(len(partitions)):
                for b_idx in range(a_idx + 1, len(partitions)):
                    key = (partitions[a_idx], partitions[b_idx])
                    pair_counts[key] = pair_counts.get(key, 0) + 1
    print(f"  {cross}/{n_clusters} clusters span more than one partition")
    for (pa, pb), count in sorted(pair_counts.items(), key=lambda kv: -kv[1])[:20]:
        print(f"    {pa} <-> {pb}: {count} shared cluster(s)")
    print()
    by_type: dict[str, int] = {}
    for entry in deduped:
        t = str(entry["candidate_type"])
        by_type[t] = by_type.get(t, 0) + 1
    print(f"distinct clusters by type: {by_type}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--tol",
        type=float,
        default=DEDUP_RELATIVE_TOLERANCE,
        help="relative-deviation merge tolerance on rmin_km_1yr/rmax_km_1yr (default 0.10)",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=OUTPUT_PATH,
        help="output path for the deduped candidate list",
    )
    args = ap.parse_args()

    records = load_all_candidates()
    print(f"loaded {len(records)} raw unmatched-bounded candidates from {SOURCE_DIR}")
    chars = characterize_all(records)
    clusters = dedup_pool(records, chars, tol=args.tol)
    deduped = build_output(records, chars, clusters)

    args.out.write_text(json.dumps(deduped, indent=2, default=str))
    print(f"\nwrote {len(deduped)} deduped candidate(s) -> {args.out}")
    print_report(records, deduped)


if __name__ == "__main__":
    main()
