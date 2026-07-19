"""#667 -- mine JPL SSD's Three-Body Periodic Orbit catalog as a DISCOVERY
INPUT (not just #647's novelty-check gate), per `#661` shortlist item 5.

See `#667`'s own `data/OUTSTANDING.md` bullet and `#661`'s own honest framing
for the full case. In short: `#647` built a GATE (does JPL already catalog
THIS ONE candidate). This script mines the same catalog the other direction
-- bulk-retrieving whole families in systems this project's OWN search
methods have never directly targeted, then classifying every retrieved
member's GEOMETRY (does it make a recurrent, close, physically-valid approach
to the secondary body) using :mod:`cyclerfinder.search.jpl_family_census`
(#667's own new module, built for this task).

TARGET SYSTEM SELECTION (honest, checked against this project's own history,
not assumed)
-----------------------------------------------------------------------------
Of JPL's 7 indexed systems, ``sun-earth``/``earth-moon`` are extensively
covered by this project's other methods (the entire #312/#563/#606/#650
lineage lives there); ``jupiter-europa``/``saturn-titan`` have some existing
coverage (Europa via the moon-tour genome; Titan via #504-adjacent and the
#647 golden itself). ``saturn-enceladus``, ``mars-phobos``, ``sun-mars`` are
the genuine roster gaps -- confirmed directly: `#609` swept Mars Phobos-
Deimos with a DIFFERENT method (patched-conic symmetric-closure + #324 bend
gate) and found it EMPTY for that method, but never touched JPL's own
CR3BP-periodic-orbit catalog; no task in this project's history has ever
queried ``saturn-enceladus`` or ``sun-mars`` at all.

FAMILY SELECTION PER SYSTEM (documented judgment, not "query all 12
blindly")
-----------------------------------------------------------------------------
All 3 target systems have TINY-to-modest secondary mass ratios (mu ~1.6e-8
for Phobos, ~1.9e-7 for Enceladus, ~3.2e-7 for Mars-about-Sun) -- confirmed
live during this task's own exploratory probes (mu/Hill-radius/L1-L2-distance
values reported inline below): for the two
MOON systems, L1/L2 sit at almost exactly the secondary's own Hill-radius
scale (Phobos: L1/L2 ~16.6 km vs. Hill radius 16.6 km vs. body radius
11.3 km -- L1/L2 essentially AT the body's surface scale; Enceladus similar,
~950 km vs. body radius 252 km), so libration-point families (halo,
lyapunov) are NOT a priori remote from the secondary the way they would be
at, say, Earth-Moon -- worth checking, not assuming. For ``sun-mars``,
L1/L2 sit ~990,000 km from Mars (radius 3389.5 km) -- three orders of
magnitude further out -- so a full halo/lyapunov family sweep is not
warranted a priori; a single targeted sanity probe at the family's
largest-amplitude (lowest-Jacobi) end is run instead, and the decision to
skip a full sweep is only taken if that probe confirms remoteness.

Chosen per system:

* ``dro``, ``dpo`` (ALL 3 systems) -- secondary-centered by construction,
  the strongest a priori candidates.
* ``halo`` (libr 1/2, branch N/S) and ``lyapunov`` (libr 1/2) for
  ``mars-phobos``/``saturn-enceladus`` only (tiny-mu L1/L2-near-secondary
  case above) -- a single targeted probe for ``sun-mars``.
* ``resonant`` -- EXPLICITLY SKIPPED. Live-probed during this task's own
  build (see the CORRECTION note this task added to
  ``verify/jpl_periodic_orbits.py``'s module docstring): ``resonant``
  requires a ``branch`` encoding a resonance ratio ("12" for 1:2, etc.),
  and the valid branch-code set is not documented and appears large/sparse
  (5 of 6 naive guesses were valid, each returning 1300-3200 members) --
  exhaustively enumerating it is a combinatorially larger, separate
  discovery task of its own, disproportionate to this census's effort
  budget given dro/dpo/halo/lyapunov already cover the "recurrent close
  secondary approach" niche this task is after.
* ``vertical``/``axial``/``longp``/``short`` -- SKIPPED for the same
  effort-budget reason (same L1/L2-family category as halo/lyapunov,
  already probed; not exhaustively re-querying every remaining variant).

GEOMETRIC CLASSIFICATION
-----------------------------------------------------------------------------
Every sampled member is propagated one full period
(``jpl_family_census.propagate_min_distances_km``, generalizing
``real_binary_kk_sweep.min_body_clearance_km`` to a full 6-D state) and
classified (``jpl_family_census.classify_secondary_approach``) against TWO
criteria already used elsewhere in this project (not invented here): a
zero-margin physical non-crash floor (`#660`'s convention) and a
Hill-fraction "genuinely close" ceiling borrowed by analogy from
``genome.hill_screen``'s existing 0.3 PASS band. See that module's own
docstring for the honest caveat on the analogy.

Large families (up to ~2400 cataloged members) are SUBSAMPLED (up to
``N_SAMPLE`` members, evenly spaced by sorted Jacobi constant) rather than
exhaustively propagated -- a representative-extent census, not an
exhaustive one; family geometry varies smoothly with Jacobi/amplitude so a
sparse sample already captures whether ANY part of the family gets close.

NO catalogue.yaml writeback. NO novelty claim -- being numerically
catalogued at JPL is not the same as being characterized as a cycler (see
`#661`'s own honest framing); a geometrically-qualifying candidate here is
reported for adjudication, never written back directly.
"""

from __future__ import annotations

import datetime
import json
import pathlib
import sys
from dataclasses import asdict

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))

from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import PreflightBlockedError, preflight_search
from cyclerfinder.search.jpl_family_census import (
    classify_secondary_approach,
    fetch_family_window,
    hill_radius_km,
)
from cyclerfinder.verify import jpl_periodic_orbits as jpo

_REPO = pathlib.Path(__file__).resolve().parent.parent
_OUT_DIR = _REPO / "data" / "found" / "667_jpl_family_census"

_PrefetchKey = tuple[str, str, int | None, str | None]
_PrefetchVal = tuple[jpo.JplSystemConstants, list[jpo.JplOrbit]]

# Generous but FINITE server-side range filters (never an unfiltered whole-
# family fetch, per #647's own server-respect convention) -- wide enough to
# cover every family's real extent at these systems (confirmed live: every
# family queried below returns rows entirely inside jacobi in [2.9, 3.1] and
# period in [0.01, 8] TU) without hardcoding a family-specific tight window.
_JACOBI_MIN = 1.5
_JACOBI_MAX = 3.3
_PERIOD_MIN = 0.001
_PERIOD_MAX = 200.0

# Measured live during this task's own build (21-orbit timing probe on
# mars-phobos dro): ~0.156 s/orbit end-to-end (fetch already cached + DOP853
# propagate + classify). Rounded up for headroom.
_MEASURED_SEC_PER_POINT = 0.20

N_SAMPLE = 100

# (system, family, libr, branch, n_sample) -- see module docstring for the
# per-system/family selection reasoning.
_COMBOS: tuple[tuple[str, str, int | None, str | None, int], ...] = (
    # mars-phobos: dro/dpo (secondary-centered) + halo/lyapunov (tiny-mu,
    # L1/L2 near-secondary-scale, worth checking).
    ("mars-phobos", "dro", None, None, N_SAMPLE),
    ("mars-phobos", "dpo", None, None, N_SAMPLE),
    ("mars-phobos", "halo", 1, "N", N_SAMPLE),
    ("mars-phobos", "halo", 1, "S", N_SAMPLE),
    ("mars-phobos", "halo", 2, "N", N_SAMPLE),
    ("mars-phobos", "halo", 2, "S", N_SAMPLE),
    ("mars-phobos", "lyapunov", 1, None, N_SAMPLE),
    ("mars-phobos", "lyapunov", 2, None, N_SAMPLE),
    # saturn-enceladus: same reasoning as mars-phobos.
    ("saturn-enceladus", "dro", None, None, N_SAMPLE),
    ("saturn-enceladus", "dpo", None, None, N_SAMPLE),
    ("saturn-enceladus", "halo", 1, "N", N_SAMPLE),
    ("saturn-enceladus", "halo", 1, "S", N_SAMPLE),
    ("saturn-enceladus", "halo", 2, "N", N_SAMPLE),
    ("saturn-enceladus", "halo", 2, "S", N_SAMPLE),
    ("saturn-enceladus", "lyapunov", 1, None, N_SAMPLE),
    ("saturn-enceladus", "lyapunov", 2, None, N_SAMPLE),
    # sun-mars: dro/dpo full treatment; halo is a SMALL targeted sanity
    # probe only (see module docstring -- L1/L2 sit ~990,000 km from Mars,
    # three orders of magnitude beyond Mars's own radius, so a full sweep
    # is not a priori warranted; 20 samples at the family's own extent is
    # enough to confirm or refute that expectation before committing to
    # more).
    ("sun-mars", "dro", None, None, N_SAMPLE),
    ("sun-mars", "dpo", None, None, N_SAMPLE),
    ("sun-mars", "halo", 1, "N", 20),
)


def _ts() -> str:
    return datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")


def _subsample_by_jacobi(orbits: list[jpo.JplOrbit], n_sample: int) -> list[jpo.JplOrbit]:
    if len(orbits) <= n_sample:
        return orbits
    ordered = sorted(orbits, key=lambda o: o.jacobi)
    idx = [round(i * (len(ordered) - 1) / (n_sample - 1)) for i in range(n_sample)]
    seen: set[int] = set()
    out = []
    for i in idx:
        if i not in seen:
            seen.add(i)
            out.append(ordered[i])
    return out


def main() -> None:
    print(f"[{_ts()}] #667 JPL periodic-orbit catalog mining census starting.")
    _OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Pre-count how many propagations this run will actually do, for an
    # honest preflight n_points + measured timing pilot.
    print(f"[{_ts()}] Pre-fetching (cached) to size the run ...")
    prefetch: dict[_PrefetchKey, _PrefetchVal] = {}
    total_points = 0
    for system, family, libr, branch, n_sample in _COMBOS:
        constants, orbits = fetch_family_window(
            system,
            family,
            libr=libr,
            branch=branch,
            jacobi_min=_JACOBI_MIN,
            jacobi_max=_JACOBI_MAX,
            period_min=_PERIOD_MIN,
            period_max=_PERIOD_MAX,
        )
        prefetch[(system, family, libr, branch)] = (constants, orbits)
        n_take = min(n_sample, len(orbits))
        total_points += n_take
        print(
            f"[{_ts()}]   {system:18s} {family:9s} libr={libr!s:4s} branch={branch!s:4s} "
            f"n_cataloged={len(orbits):5d} n_sampled={n_take:4d}"
        )

    preflight_search(
        task_no=667,
        region_id="667-jpl-periodic-orbit-catalog-mining-census-2026-07-20",
        method=MethodCapability(
            genome=(
                "JPL SSD Three-Body Periodic Orbit catalog bulk mining "
                "(jpl_family_census.fetch_family_window) + geometric recurrent-"
                "close-secondary-approach classification "
                "(jpl_family_census.classify_secondary_approach)"
            ),
            corrector=(
                "none -- JPL's own published ICs propagated one period via DOP853 (no correction)"
            ),
            capability_tags=frozenset(
                {
                    "cr3bp",
                    "jpl-catalog-mined",
                    "secondary-centered-geometry",
                    "ballistic",
                    "coplanar-and-3d",
                }
            ),
            git_sha="working-tree",
        ),
        script_path=pathlib.Path(__file__),
        n_points=total_points,
        timing_pilot_seconds_per_point=_MEASURED_SEC_PER_POINT,
        # #667 IS registered in data/OUTSTANDING.md (its own bullet, "- **#667
        # (registered 2026-07-19, not yet dispatched)** ..."), but
        # preflight.py's own _TASK_ALLOCATION_RE regex requires the literal
        # CONTIGUOUS "**#NNN**" (nothing between the digits and the closing
        # "**"), which every post-#645 bullet header style (a parenthetical
        # registration date, a "✓ DONE (date, model)" annotation, etc.)
        # breaks -- confirmed directly: the regex currently recognizes a
        # max task number of #645 across the WHOLE file. This is a checker
        # false-negative (the real hygiene condition -- "is this task
        # actually registered" -- is satisfied), not a genuine gap, so the
        # audited override escape hatch is the right tool here rather than
        # reformatting the bullet just to game the regex.
        override_reason=(
            "#667 is registered in data/OUTSTANDING.md's own bullet, but "
            "preflight.py's _TASK_ALLOCATION_RE only matches a literal "
            "contiguous '**#NNN**' with nothing in between; every post-#645 "
            "bullet header (date/model/status annotations inside the bold) "
            "breaks that regex -- confirmed the checker's recognized task "
            "numbers currently top out at #645. Checker false-negative, not "
            "a real registration gap."
        ),
    )

    print(f"[{_ts()}] Preflight clear. Total propagations this run: {total_points}")

    all_verdicts_path = _OUT_DIR / "verdicts.jsonl"
    summary: list[dict[str, object]] = []
    n_done = 0

    with all_verdicts_path.open("w", encoding="utf-8") as fh:
        for system, family, libr, branch, n_sample in _COMBOS:
            constants, orbits = prefetch[(system, family, libr, branch)]
            sample = _subsample_by_jacobi(orbits, n_sample)
            r_hill = hill_radius_km(constants.mu, constants.lunit_km)
            print(
                f"[{_ts()}] Classifying {system}/{family} libr={libr} branch={branch}: "
                f"{len(sample)} sampled of {len(orbits)} cataloged "
                f"(mu={constants.mu:.4e}, hill_radius_km={r_hill:.2f}, "
                f"radius_secondary_km={constants.radius_secondary_km})"
            )
            n_close = 0
            n_invalid = 0
            n_stable = 0
            min_hill_fraction = float("inf")
            for orbit in sample:
                v = classify_secondary_approach(
                    orbit, constants, system=system, family=family, libr=libr, branch=branch
                )
                fh.write(json.dumps(asdict(v)) + "\n")
                n_done += 1
                if not v.physically_valid:
                    n_invalid += 1
                if v.is_close_approach:
                    n_close += 1
                if v.stability <= 1.0 + 1e-6:
                    n_stable += 1
                min_hill_fraction = min(min_hill_fraction, v.hill_fraction)

            print(
                f"[{_ts()}]   -> {n_close}/{len(sample)} close-approach-qualifying, "
                f"{n_invalid}/{len(sample)} physically-invalid (should be 0), "
                f"{n_stable}/{len(sample)} linearly-stable (stability<=1), "
                f"min_hill_fraction_seen={min_hill_fraction:.4f} "
                f"[running total propagated: {n_done}/{total_points}]"
            )
            summary.append(
                {
                    "system": system,
                    "family": family,
                    "libr": libr,
                    "branch": branch,
                    "n_cataloged": len(orbits),
                    "n_sampled": len(sample),
                    "n_close_approach_qualifying": n_close,
                    "n_physically_invalid": n_invalid,
                    "n_linearly_stable_stability_le_1": n_stable,
                    "min_hill_fraction_seen": min_hill_fraction,
                    "hill_radius_km": r_hill,
                    "radius_secondary_km": constants.radius_secondary_km,
                    "mu": constants.mu,
                }
            )

    summary_path = _OUT_DIR / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))

    print(f"\n[{_ts()}] Census complete. {n_done} orbits propagated + classified.")
    print(f"[{_ts()}] Per-orbit verdicts: {all_verdicts_path}")
    print(f"[{_ts()}] Per-combo summary:  {summary_path}")
    print(f"[{_ts()}] ---- SUMMARY TABLE ----")
    for row in summary:
        print(
            f"  {row['system']:18s} {row['family']:9s} libr={row['libr']!s:4s} "
            f"branch={row['branch']!s:4s} "
            f"close={row['n_close_approach_qualifying']}/{row['n_sampled']} "
            f"(of {row['n_cataloged']} cataloged) "
            f"min_hill_frac={row['min_hill_fraction_seen']:.4f}"
        )


if __name__ == "__main__":
    try:
        main()
    except PreflightBlockedError as exc:
        print(f"[{_ts()}] BLOCKED by preflight_search:\n{exc}")
        sys.exit(1)
