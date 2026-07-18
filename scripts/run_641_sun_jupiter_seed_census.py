#!/usr/bin/env python3
"""#641: Sun-Jupiter CR3BP periodic-orbit family CENSUS using #628's productionized
generative seed model.

`#628` built a reusable API (``cyclerfinder.ml.seed_generation.
generate_and_refine_seeds``) on top of `#608`'s statistical generative seed
model (trained on this project's own Earth-Moon CR3BP corrector-outcome
corpus). `#624` independently verified the model's convergence-lift advantage
genuinely TRANSFERS to mass ratios it never trained on -- most strongly at
mu=0.001 (30x over blind seeding), more weakly at Sun-Earth mu~3e-6 (3.5x).
Neither task pointed this capability at a genuine, previously-unswept
discovery target; both were evaluation-only. This script is that first real
application: a CENSUS (not a targeted search for one specific family) of
whatever real, physically-sane periodic-orbit families the model's seeds
converge to at Sun-Jupiter, mu ~ 9.5388e-4 -- chosen deliberately because it
sits close to `#624`'s STRONGEST validated cross-mu anchor (mu=0.001,
delta_log10_mu~1.085), not the weakest, and because Sun-Jupiter CR3BP itself
was never the direct target of any prior task in this project's #600-#640
arc.

**This is explicitly NOT ``scripts/search_campaign_daemon.py`` Phase B** --
that script is Earth-Moon-only (in-distribution, not a transfer test) and its
``correct_periodic`` calls auto-log into the model's OWN training corpus by
default (``CYCLERFINDER_OUTCOME_LOG``), which would create exactly the
distribution-feedback-loop hazard `#634` already documented. This script
never sets that env var and never touches ``out/outcome_log/``.

**No catalogue writeback** -- raw CR3BP periodic orbits are not cyclers in
this project's own ``cycler``/``quasi_cycler``/``precursor_mga`` taxonomy.
This is a survey/census script; any adjudication of a genuinely novel
survivor happens in ``data/OUTSTANDING.md``, not here.

**preflight_search() exemption**: this is a FIXED-N capability census
(generate N seeds, refine, cluster, literature-check) with no
``region_id``/``n_points`` sweep-region concept to preflight -- the same
category as `#608`/`#614`/`#317`/`#624` before it (see
``tests/scripts/test_scripts_call_preflight.py``'s ``_LEGACY_EXEMPT`` entry
for this file).

**Pilot-run finding (this task's own discovery, not anticipated by `#628`/
`#624`)**: a large fraction of "converged AND physically-sane" seeds at this
mu are NOT genuine periodic orbits at all -- they are the trivial L4/L5
equilateral (and L1/L2/L3 collinear) Lagrange EQUILIBRIUM points. A fixed
point trivially satisfies ``correct_periodic``'s periodicity residual for
*any* period guess (nothing moves), and ``is_physically_sane`` has no
velocity-magnitude check, so it does not reject them. This is a genuine
convergence-target artifact of applying the corrector at a genuinely new mu
that this task ran into, not a documented `#608`/`#624`/`#628` caveat --
:func:`is_degenerate_equilibrium` filters these out before clustering.
"""

from __future__ import annotations

import json
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from cyclerfinder.core.cr3bp import cr3bp_system
from cyclerfinder.ml.seed_generation import (
    GeneratedSeed,
    SeedGenerationReport,
    expected_lift_for_mu,
    generate_and_refine_seeds,
)
from cyclerfinder.search.literature_check import (
    KNOWN_CORPUS,
    CandidateSignature,
    SearchResult,
    check_literature,
)

# `GeneratedSeed` is re-exported (not just used internally) because
# `tests/scripts/test_run_641_sun_jupiter_seed_census.py` constructs instances
# of it directly to exercise the clustering logic without a real corrector
# run -- same pattern `seed_generation.py` uses for its own `cr3bp` re-export
# (see that module's comment). `__all__` makes this explicit for mypy's
# strict-mode implicit-reexport check.
__all__ = ["GeneratedSeed"]

_REPO_ROOT = Path(__file__).resolve().parent.parent
_OUT_DIR = _REPO_ROOT / "data" / "found" / "641_sun_jupiter_seed_census"

PILOT_N = 100
FULL_N = 1000
SEED = 641

# Observed in this task's own pilot run: genuine converged periodic orbits had
# |v0| ~ 0.13 (nondimensional); trivial Lagrange-equilibrium "convergences"
# had |v0| in the 1e-12 .. 1e-17 range (Newton-iteration numerical noise
# around an exact zero). 1e-6 sits ~6 orders of magnitude below the smallest
# genuine orbit seen and ~6 orders above the largest equilibrium noise floor
# -- a wide, well-separated margin, not a knife-edge threshold.
EQUILIBRIUM_VNORM_THRESHOLD = 1e-6

# Clustering tolerances for "distinct family" grouping (task step 4: "Jacobi
# constant + period + qualitative geometry"). Jacobi is rounded to 3 decimals
# (this model/corrector's converged Jacobi values repeat to ~1e-13 precision
# within a single basin, so 3 decimals safely separates distinct basins
# without over-fragmenting one), period to 1 decimal (period bounds span
# 0.5-15.0, so 1 decimal gives ~145 distinguishable bins across the range --
# coarse enough that small solver-noise differences in period do not split a
# single family into many spurious clusters).
JACOBI_ROUND = 3
PERIOD_ROUND = 1
# "Near zero" deadband for the qualitative y0/z0 sign tags below -- values
# smaller than this in absolute terms are tagged "0" (on-axis / in-plane)
# rather than "+"/"-", matching the deadband `heuristic_family_tag`
# (`cyclerfinder.ml.orbit_generative`) uses for its own z0 tag.
GEOMETRY_DEADBAND = 1e-3


def is_degenerate_equilibrium(
    state0: NDArray[np.float64] | list[float],
    *,
    vnorm_threshold: float = EQUILIBRIUM_VNORM_THRESHOLD,
) -> bool:
    """Is this converged ``state0`` a trivial Lagrange-point equilibrium rather
    than a genuine periodic orbit?

    A fixed point of the CR3BP rotating-frame equations (velocity identically
    zero) trivially satisfies any periodicity residual for an arbitrary
    period guess -- the corrector reports ``converged=True`` and the guess's
    OWN period survives unchanged, because nothing ever moves. See this
    module's docstring for how this was discovered (this task's own pilot
    run, not a documented #608/#624/#628 caveat).
    """
    v = np.asarray(state0, dtype=np.float64)[3:6]
    return bool(np.linalg.norm(v) < vnorm_threshold)


def lagrange_point_label(state0: NDArray[np.float64] | list[float], mu: float) -> str:
    """Best-effort classical label (L1-L5) for a degenerate-equilibrium seed.

    Only meaningful when :func:`is_degenerate_equilibrium` is True for this
    seed -- a genuine orbit passing near a Lagrange point in position is NOT
    relabelled by this function (callers must gate on the equilibrium check
    first). Primary sits at ``x = -mu``, secondary at ``x = 1 - mu``.
    """
    x0, y0 = float(state0[0]), float(state0[1])
    if abs(y0) > GEOMETRY_DEADBAND:
        return "L4" if y0 > 0 else "L5"
    # Collinear ordering along the x-axis: L3 < primary(-mu) < L1 < secondary(1-mu) < L2.
    if x0 > 1.0 - mu:
        return "L2"
    if x0 > -mu:
        return "L1"
    return "L3"


def _geometry_tag(state0: NDArray[np.float64]) -> tuple[str, str]:
    y0, z0 = float(state0[1]), float(state0[2])
    y_tag = "y0" if abs(y0) < GEOMETRY_DEADBAND else ("y+" if y0 > 0 else "y-")
    z_tag = "z0" if abs(z0) < GEOMETRY_DEADBAND else ("z+" if z0 > 0 else "z-")
    return y_tag, z_tag


@dataclass(frozen=True)
class FamilyCluster:
    """One distinct candidate family: seeds sharing (rounded Jacobi, rounded
    period, qualitative y0/z0-sign geometry).
    """

    key: tuple[float, float, str, str]
    members: list[GeneratedSeed]
    representative: GeneratedSeed


def cluster_genuine_orbits(
    seeds: list[GeneratedSeed],
    *,
    jacobi_round: int = JACOBI_ROUND,
    period_round: int = PERIOD_ROUND,
) -> list[FamilyCluster]:
    """Cluster CONVERGED, PHYSICALLY-SANE, NON-degenerate seeds into distinct
    candidate families by (Jacobi constant, period, qualitative geometry).

    Callers must pre-filter to ``converged and physically_sane and not
    is_degenerate_equilibrium(...)`` seeds -- this function does not
    re-check either condition. Deliberately simple/inspectable per this
    task's "light clustering pass" scope (task step 4): two genuinely
    different families landing in the same rounded bin would get merged,
    and one family straddling a bin boundary could get split -- the same
    honest limitation `heuristic_family_tag` documents for its own coarser
    binning. Returned list is sorted largest-cluster-first.
    """
    groups: dict[tuple[float, float, str, str], list[GeneratedSeed]] = {}
    for s in seeds:
        assert s.state0 is not None and s.jacobi is not None and s.period is not None
        y_tag, z_tag = _geometry_tag(s.state0)
        key = (round(s.jacobi, jacobi_round), round(s.period, period_round), y_tag, z_tag)
        groups.setdefault(key, []).append(s)

    clusters = []
    for key, members in groups.items():
        rep = min(members, key=lambda m: m.residual if m.residual is not None else float("inf"))
        clusters.append(FamilyCluster(key=key, members=members, representative=rep))
    return sorted(clusters, key=lambda c: -len(c.members))


def _offline_corpus_hits() -> list[SearchResult]:
    """Same construction `run_436_direct_er3bp.py` uses: turn this project's own
    curated cycler-literature corpus (``KNOWN_CORPUS``) into fake "search hits"
    so :func:`check_literature` can run offline/deterministically in a script,
    per this project's `#436` precedent for a raw (non-cycler) periodic-orbit
    candidate.
    """
    hits: list[SearchResult] = []
    for anchor in KNOWN_CORPUS:
        bodies_str = " ".join(sorted(anchor.body_set))
        snippet = (
            f"{anchor.citation}. Authors: {', '.join(anchor.authors)}. "
            f"Keywords: {', '.join(anchor.keywords)}. Bodies: {bodies_str}. "
            f"Primary: {anchor.primary}. Subject: cycler trajectory mission design."
        )
        hits.append(
            SearchResult(
                title=anchor.name + " (cycler trajectory)",
                url=anchor.doi or f"https://example.org/{anchor.name.replace(' ', '_')}",
                snippet=snippet,
            )
        )
    return hits


_OFFLINE_HITS = _offline_corpus_hits()


def _tokenise(s: str) -> set[str]:
    return {t for t in "".join(c.lower() if c.isalnum() else " " for c in s).split() if len(t) > 2}


def offline_search(query: str) -> list[SearchResult]:
    q_terms = _tokenise(query)
    out: list[tuple[int, SearchResult]] = []
    for r in _OFFLINE_HITS:
        text_terms = _tokenise(r.title + " " + r.snippet)
        overlap = len(q_terms & text_terms)
        if overlap >= 2:
            out.append((overlap, r))
    out.sort(key=lambda t: t[0], reverse=True)
    return [r for _, r in out]


def cluster_topology_label(rep: GeneratedSeed) -> frozenset[str]:
    """Best-effort ``topology_label`` for a cluster's representative -- "halo"
    if it has a non-negligible out-of-plane component, otherwise unrestricted
    (empty set) since none of ``CandidateSignature``'s standard labels
    (repeated-moon/pump-tour/mga-tour/tulip/nrho/resonant/binary-coorbital)
    honestly describe a raw planar Lyapunov/DRO-type orbit -- inventing one
    would bias the anchor match rather than leaving it open, per this
    dataclass's own "empty means topology unrestricted" documented default.
    """
    assert rep.state0 is not None
    if abs(float(rep.state0[2])) > GEOMETRY_DEADBAND:
        return frozenset({"halo"})
    return frozenset()


def literature_check_clusters(clusters: list[FamilyCluster]) -> list[dict[str, object]]:
    """Mandatory literature-novelty check (per
    ``[[feedback_literature_novelty_check_baseline]]``): run
    ``search/literature_check.py``'s :func:`check_literature` against each
    distinct cluster's representative member.

    IMPORTANT HONEST CAVEAT (discovered running this task): ``KNOWN_CORPUS``
    -- and therefore every query :func:`build_queries` generates -- is scoped
    to published CYCLER trajectories (every query is suffixed "cycler
    trajectory"/"cycler resonance"/etc.). Raw CR3BP periodic-orbit families
    like Lyapunov/halo/DRO are NOT cyclers and essentially never described
    with that vocabulary in their own literature, so a "not-found" verdict
    here reflects a QUERY-VOCABULARY MISMATCH, not evidence these families
    are unpublished -- do not read a not-found from this function as support
    for novelty. See this script's own module docstring / the task's
    OUTSTANDING.md bullet for the supplementary manual grounding (JPL SSD's
    public Three-Body Periodic Orbits catalog, ``ssd-api.jpl.nasa.gov``,
    explicitly indexes Sun-Jupiter planar Lyapunov / halo / vertical / DRO /
    long-period families) that this mismatch makes necessary.
    """
    results: list[dict[str, object]] = []
    for c in clusters:
        rep = c.representative
        sig = CandidateSignature(
            primary="Sun",
            sequence=("Sun", "Jupiter"),
            topology_label=cluster_topology_label(rep),
        )
        lit = check_literature(sig, search=offline_search)
        results.append(
            {
                "cluster_key": list(c.key),
                "status": lit.status,
                "citation": lit.citation,
                "confidence": lit.confidence,
                "notes": lit.notes,
            }
        )
    return results


def _print(msg: str) -> None:
    print(f"[{time.strftime('%Y-%m-%dT%H:%M:%S')}] {msg}", flush=True)


def _seed_to_dict(s: GeneratedSeed) -> dict[str, object]:
    return {
        "converged": s.converged,
        "physically_sane": s.physically_sane,
        "state0": s.state0.tolist() if s.state0 is not None else None,
        "period": s.period,
        "jacobi": s.jacobi,
        "residual": s.residual,
        "stability_index": s.stability_index,
        "is_equilibrium": (is_degenerate_equilibrium(s.state0) if s.state0 is not None else None),
    }


def _summarize_batch(label: str, report: SeedGenerationReport, elapsed: float) -> dict[str, object]:
    sane_seeds = [s for s in report.seeds if s.converged and s.physically_sane]
    equilibria = [s for s in sane_seeds if is_degenerate_equilibrium(s.state0)]  # type: ignore[arg-type]
    genuine = [s for s in sane_seeds if not is_degenerate_equilibrium(s.state0)]  # type: ignore[arg-type]
    _print(
        f"  [{label}] n_attempted={report.n_attempted} n_converged={report.n_converged} "
        f"n_physically_sane={report.n_physically_sane} "
        f"(rate={report.converged_and_physically_sane_rate:.1%}) "
        f"in {elapsed:.1f}s ({elapsed / max(report.n_attempted, 1):.3f}s/seed)"
    )
    _print(
        f"  [{label}] of {len(sane_seeds)} sane seeds: {len(equilibria)} degenerate "
        f"Lagrange-equilibrium 'convergences', {len(genuine)} genuine periodic orbits"
    )
    return {
        "label": label,
        "n_requested": report.n_requested,
        "n_attempted": report.n_attempted,
        "n_converged": report.n_converged,
        "n_physically_sane": report.n_physically_sane,
        "converged_and_physically_sane_rate": report.converged_and_physically_sane_rate,
        "n_degenerate_equilibrium": len(equilibria),
        "n_genuine_orbit": len(genuine),
        "elapsed_s": elapsed,
    }


def main() -> None:
    t_start = time.time()
    _print("#641 Sun-Jupiter generative-seed census starting")

    system = cr3bp_system("Sun", "Jupiter")
    _print(f"  Sun-Jupiter mu = {system.mu!r} (l_km={system.l_km:.3e}, t_s={system.t_s:.3e})")

    lift = expected_lift_for_mu(system.mu)
    _print(
        f"  expected_lift_for_mu: estimated_lift={lift.estimated_lift:.2f}x, "
        f"delta_log10_mu={lift.delta_log10_mu:.3f}, "
        f"beyond_validated_range={lift.beyond_validated_range}"
    )
    _print(f"  method: {lift.method}")

    # --- Pilot phase: gauge per-seed cost before committing to the full batch ---
    _print(f"\n=== pilot batch: N={PILOT_N} ===")
    t0 = time.time()
    pilot_report = generate_and_refine_seeds(
        PILOT_N, primary="Sun", secondary="Jupiter", rng=np.random.default_rng(SEED * 10)
    )
    pilot_elapsed = time.time() - t0
    pilot_summary = _summarize_batch("pilot", pilot_report, pilot_elapsed)

    # --- Full census batch ---
    _print(f"\n=== full census batch: N={FULL_N} ===")
    t0 = time.time()
    full_report = generate_and_refine_seeds(
        FULL_N, primary="Sun", secondary="Jupiter", rng=np.random.default_rng(SEED)
    )
    full_elapsed = time.time() - t0
    full_summary = _summarize_batch("full", full_report, full_elapsed)

    sane_seeds = [s for s in full_report.seeds if s.converged and s.physically_sane]
    genuine = [
        s for s in sane_seeds if s.state0 is not None and not is_degenerate_equilibrium(s.state0)
    ]
    equilibria = [
        s for s in sane_seeds if s.state0 is not None and is_degenerate_equilibrium(s.state0)
    ]

    eq_labels: dict[str, int] = {}
    for s in equilibria:
        assert s.state0 is not None
        label = lagrange_point_label(s.state0, system.mu)
        eq_labels[label] = eq_labels.get(label, 0) + 1
    _print(f"  equilibrium breakdown (full batch): {eq_labels}")

    clusters = cluster_genuine_orbits(genuine)
    _print(f"\n=== clustering: {len(clusters)} distinct genuine-orbit family cluster(s) ===")
    for i, c in enumerate(clusters):
        rep = c.representative
        assert rep.state0 is not None
        _print(
            f"  cluster {i}: key={c.key} n_members={len(c.members)} "
            f"rep: jacobi={rep.jacobi:.6f} period={rep.period:.6f} "
            f"state0={np.array2string(rep.state0, precision=6)} "
            f"stability={rep.stability_index}"
        )

    _print("\n=== mandatory literature-novelty check (per cluster) ===")
    lit_results = literature_check_clusters(clusters)
    for r in lit_results:
        _print(f"  cluster {r['cluster_key']}: status={r['status']} citation={r['citation']!r}")
    _print(
        "  NOTE: KNOWN_CORPUS/build_queries is scoped to published CYCLER "
        "trajectories (query vocabulary always includes 'cycler') -- a "
        "not-found here reflects that vocabulary mismatch for raw CR3BP "
        "periodic-orbit families, NOT evidence of novelty. See this script's "
        "module docstring / the #641 OUTSTANDING.md bullet for the "
        "supplementary manual grounding against JPL SSD's public Three-Body "
        "Periodic Orbits catalog."
    )

    elapsed = time.time() - t_start
    _OUT_DIR.mkdir(parents=True, exist_ok=True)

    summary = {
        "task": "#641",
        "system": {"primary": "Sun", "secondary": "Jupiter", "mu": system.mu},
        "lift_estimate": {
            "estimated_lift": lift.estimated_lift,
            "delta_log10_mu": lift.delta_log10_mu,
            "beyond_validated_range": lift.beyond_validated_range,
            "method": lift.method,
            "caveat": lift.caveat,
        },
        "pilot": pilot_summary,
        "full": full_summary,
        "n_clusters": len(clusters),
        "clusters": [
            {
                "key": list(c.key),
                "n_members": len(c.members),
                "representative": _seed_to_dict(c.representative),
            }
            for c in clusters
        ],
        "equilibrium_breakdown": eq_labels,
        "literature_check": lit_results,
        "seed": SEED,
        "elapsed_s": elapsed,
    }
    summary_path = _OUT_DIR / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=False) + "\n")

    detail_path = _OUT_DIR / "full_batch_seeds.jsonl"
    with detail_path.open("w") as fh:
        for s in full_report.seeds:
            fh.write(json.dumps(_seed_to_dict(s)) + "\n")

    _print(f"\ndone in {elapsed:.1f}s")
    _print(f"  wrote {summary_path}")
    _print(f"  wrote {detail_path}")
    _print(
        "\nHONEST HEADLINE: this is a CENSUS, not a targeted search for a specific "
        "family. No catalogue writeback performed. Literature novelty check runs "
        "separately (see this task's own scripts/run_641 companion analysis / "
        "OUTSTANDING.md #641 bullet for the verdict)."
    )


if __name__ == "__main__":
    main()
