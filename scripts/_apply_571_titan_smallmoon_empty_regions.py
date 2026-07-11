"""Apply the #571 Titan-Mimas/Enceladus/Tethys/Dione analytically-empty stamps.

One-shot writer (prefixed `_` so it is treated as scratch and left untouched
by the harness): appends 4 EmptyRegionReport records to
data/empty_regions.jsonl for the four Titan-small-moon pairs the 2026-07-12
Fable plan review of #571 (``data/OUTSTANDING.md``) found are analytically
empty under this project's own two-sided #324 physical max-bend gate --
BEFORE spending any sweep compute, not after a sweep found nothing.

Numbers used here are RECOMPUTED independently by
``scripts/verify_571_gate_analytics.py`` (not copied from OUTSTANDING.md's
summary), using ``core/flyby.py::max_bend`` + ``core/satellites.py``'s
sourced GM/radius/sma/safe_alt values. They match Fable's figures to the
displayed precision.

Gate logic (``search/physical_sanity.py::candidate_passes_physical_gate``,
as applied by ``scan_558_uranus_all_pairs_offset_sweep.py::gate_candidate``
to every 3-encounter (anchor, flyby, anchor) sequence): EVERY encounter must
clear ``DEFAULT_MIN_USEFUL_BEND_DEG=5.0`` deg of ballistic bend -- the small
moon's OWN encounter (index 1 of the 3), not just Titan's, must clear 5 deg.
For each of these four small moons, the moon's own minimum-achievable
``V_inf`` (over ANY conic reaching Titan's orbital radius, realized by the
Hohmann-transfer tangential minimum-energy ellipse -- the single most
favorable-to-feasibility geometry, so if THIS fails, every less-favorable
Lambert-arc grid point the actual sweep would test fails too) exceeds the
``V_inf`` at which the moon's own 5 deg bend ceiling sits. No basin can ever
exist at any grid point in any (rel_offset, phase0, tof_scale, n_rev) sweep
of this genome for these four pairs.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from cyclerfinder.core.flyby import max_bend  # noqa: E402
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES  # noqa: E402
from cyclerfinder.data.empty_regions import (  # noqa: E402
    DEFAULT_EMPTY_REGIONS_PATH,
    EmptyRegionReport,
    append_empty_region,
)
from cyclerfinder.data.method_capability import MethodCapability  # noqa: E402
from cyclerfinder.search.physical_sanity import DEFAULT_MIN_USEFUL_BEND_DEG  # noqa: E402

git_sha = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"]).decode().strip()

MU_SATURN = PRIMARIES["Saturn"]
TITAN = SATELLITES["Titan"]


def _bend_ceiling_vinf(body: str, target_bend_deg: float) -> float:
    import math

    sat = SATELLITES[body]
    mu = sat.mu_km3_s2
    rp = sat.radius_eq_km + sat.safe_alt_km
    target_rad = math.radians(target_bend_deg)
    lo, hi = 1e-6, 50.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if max_bend(mu, rp, mid) > target_rad:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _hohmann_vinf_at_r1(r1_km: float, r2_km: float, mu_primary: float) -> float:
    import math

    a_t = 0.5 * (r1_km + r2_km)
    v_circ = math.sqrt(mu_primary / r1_km)
    v_transfer = math.sqrt(mu_primary * (2.0 / r1_km - 1.0 / a_t))
    return abs(v_transfer - v_circ)


def _bend_at_vinf(body: str, vinf: float) -> float:
    import math

    sat = SATELLITES[body]
    return math.degrees(max_bend(sat.mu_km3_s2, sat.radius_eq_km + sat.safe_alt_km, vinf))


SMALL_MOONS = ("Mimas", "Enceladus", "Tethys", "Dione")

method = MethodCapability(
    genome="#558-corrected discovery genome (circular-coplanar Kepler moon states + "
    "patched-conic Lambert, rel_offset x phase0 x tof_scale x n_rev grid, "
    "scripts/scan_558_uranus_all_pairs_offset_sweep.py, genericized to a "
    "primary= parameter for #571) -- NOT SWEPT for these 4 pairs, gate "
    "analytically evaluated instead (see module docstring).",
    corrector="none (no sweep run; this is an analytic pre-sweep gate evaluation, not a "
    "correction/closure)",
    capability_tags=frozenset(
        {
            "ballistic",
            "coplanar",
            "patched-conic",
            "single-arc",
            "analytic-gate-precheck",
        }
    ),
    git_sha=git_sha,
)

results = {}
for moon in SMALL_MOONS:
    sat = SATELLITES[moon]
    ceiling = _bend_ceiling_vinf(moon, DEFAULT_MIN_USEFUL_BEND_DEG)
    min_vinf = _hohmann_vinf_at_r1(sat.sma_km, TITAN.sma_km, MU_SATURN)
    bend_at_min = _bend_at_vinf(moon, min_vinf)
    results[moon] = {
        "gm_km3_s2": sat.mu_km3_s2,
        "sma_km": sat.sma_km,
        "safe_alt_km": sat.safe_alt_km,
        "bend_ceiling_vinf_kms": ceiling,
        "min_achievable_vinf_kms": min_vinf,
        "bend_at_min_achievable_deg": bend_at_min,
        "feasible": bool(min_vinf <= ceiling),
    }
    assert not results[moon]["feasible"], f"{moon} unexpectedly feasible -- do not stamp!"

for moon in SMALL_MOONS:
    r = results[moon]
    region_id = f"saturn-titan-{moon.lower()}-analytically-empty-571"
    report = EmptyRegionReport(
        region_id=region_id,
        family=f"Titan-{moon}-Titan (#558-corrected genome, small-moon-flyby gate-infeasible pair)",
        centre="Saturn",
        topologies=(
            {
                "sequence": ["Titan", moon, "Titan"],
                "per_leg_revs": [0, 1, 2, 3],
                "period_k": 1,
            },
        ),
        method_capability=method,
        search_extent={
            "n_epochs": 0,
            "span_days": 0.0,
            "n_topologies": 1,
            # Bounded-but-not-executed: the analytic pre-check covers the
            # FULL grid this genome's sweep would ever explore -- every
            # (rel_offset, phase0, tof_scale, n_rev) point shares the SAME
            # moon-radius pair, so the same infeasibility holds everywhere.
            # points_total is reported as the #558-style production grid
            # size that a real sweep of this pair WOULD have covered (so the
            # negative is bounded/auditable per validate_empty_region), not
            # as "1" (which would misrepresent this as a single-point check).
            "points_total": 360 * 51,  # n_offset=360 x tof_scales(dense)=51, n_phase=1
            "ephem_model": "patched-conic (circular coplanar Kepler)",
            "center": "Saturn",
            "n_offset": 360,
            "n_phase": 1,
            "tof_scales_dense_count": 51,
            "n_rev_range": [0, 3],
            "moon_gm_km3_s2": r["gm_km3_s2"],
            "moon_sma_km": r["sma_km"],
            "moon_safe_alt_km": r["safe_alt_km"],
            "titan_gm_km3_s2": TITAN.mu_km3_s2,
            "titan_sma_km": TITAN.sma_km,
        },
        prune_gates=(
            "two-sided #324 physical max-bend gate "
            "(search/physical_sanity.py::candidate_passes_physical_gate), EVERY "
            "encounter incl. the small moon's own must clear "
            f"DEFAULT_MIN_USEFUL_BEND_DEG={DEFAULT_MIN_USEFUL_BEND_DEG} deg",
            "moon's own min-achievable V_inf (Hohmann-transfer tangential minimum-energy "
            "bound, the single MOST FAVORABLE-to-feasibility conic connecting the moon's "
            "orbit radius to Titan's) vs the moon's own 5 deg bend ceiling",
        ),
        result={
            "bend_ceiling_vinf_kms": round(r["bend_ceiling_vinf_kms"], 4),
            "min_achievable_vinf_kms": round(r["min_achievable_vinf_kms"], 4),
            "bend_at_min_achievable_vinf_deg": round(r["bend_at_min_achievable_deg"], 4),
            "gap_factor": round(r["min_achievable_vinf_kms"] / r["bend_ceiling_vinf_kms"], 2),
        },
        verdict=f"ANALYTICALLY EMPTY -- Titan-{moon} fails the two-sided #324 physical "
        "max-bend gate at EVERY possible grid point (the moon's own encounter can never "
        "clear 5 deg bend at any V_inf reachable by a Lambert arc to Titan's radius). "
        "No compute spent sweeping.",
        interpretation=(
            f"{moon}'s own 5-deg-bend V_inf ceiling is {r['bend_ceiling_vinf_kms']:.4f} km/s "
            f"(GM={r['gm_km3_s2']} km^3/s^2, safe periapsis radius = radius_eq + "
            f"{r['safe_alt_km']} km safe_alt, both from core/satellites.py, JPL SSD-sourced). "
            f"The MINIMUM achievable V_inf at {moon}, over ANY conic connecting {moon}'s "
            f"orbit (a={r['sma_km']:.1f} km) to Titan's orbit (a={TITAN.sma_km:.1f} km) about "
            f"Saturn (GM={MU_SATURN} km^3/s^2), is {r['min_achievable_vinf_kms']:.4f} km/s "
            "-- realized by the Hohmann-transfer tangential minimum-energy ellipse, the "
            "single most favorable-to-feasibility geometry (any non-tangential Lambert arc "
            "the actual sweep would evaluate has a LARGER velocity mismatch with the moon's "
            "local circular velocity, so this is a strict floor, not an approximation). At "
            f"this best-case V_inf, the achievable bend is only "
            f"{r['bend_at_min_achievable_deg']:.4f} deg -- "
            f"{r['min_achievable_vinf_kms'] / r['bend_ceiling_vinf_kms']:.1f}x above the "
            f"ceiling V_inf, so nowhere near the 5 deg floor. Because the #558-corrected "
            "genome's residual/gate structure evaluates the SAME two body-radius pair "
            "(Titan's orbit, this moon's orbit) at every (rel_offset, phase0, tof_scale, "
            "n_rev) grid point -- only the RELATIVE GEOMETRY (which V_inf, at which angle, "
            "the sweep's dimensions vary) changes, never the pair of orbit radii the Lambert "
            "arc must connect -- this Hohmann floor bounds every grid point uniformly: no "
            "point in the full (rel_offset x phase0 x tof_scale x n_rev) grid can achieve a "
            f"lower V_inf at {moon} than this analytic floor, so none can clear the gate "
            "either. This is a genuine physical limit (mass-deficiency), consistent with "
            "#489's prior finding that these same four moons already fail an equivalent "
            "gate for OTHER pairings (Tethys 0.44 deg, Dione ~3.1-3.2 deg max-bend at their "
            "#489-era V_inf regimes) -- this #571 stamp extends that same physical "
            "conclusion to the specific Titan-paired case, not a re-litigation of it.\n\n"
            "CAVEAT (conditional on this project's own gate POLICY, not a universal "
            "physical claim): this is conditional on the TWO-SIDED gate this project's "
            "genome enforces (every encounter, including the small moon's own, must bend "
            "usefully). Russell & Strange (2009)'s own published Titan-Enceladus census "
            "(19 members) uses Enceladus as a PASSIVE SCIENCE TARGET, not as a "
            "trajectory-bending flyby node -- an architecture where the small moon's "
            "encounter does NOT need to clear a bend-usefulness floor at all, because Titan "
            "alone does the trajectory work and the small-moon pass is opportunistic "
            "science, not propulsion. That is a STRUCTURALLY DIFFERENT, uncovered case for "
            "this genome (which by construction requires every node in the cycle to "
            "contribute useful bending), not a contradiction of this empty-region finding."
        ),
        source_anchors="core/flyby.py::max_bend (ballistic bend-angle formula, Bate-Mueller-"
        "White Sec 6.4, this project's existing implementation); core/satellites.py sourced "
        "GM/radius/sma/safe_alt for Titan and the four small moons (JPL SSD gm_de440 + "
        "SAT441 planetary/satellite constants, cited in-file); "
        "search/physical_sanity.py::DEFAULT_MIN_USEFUL_BEND_DEG=5.0 deg (this project's "
        "existing #324 gate floor, unchanged); numbers independently recomputed by "
        "scripts/verify_571_gate_analytics.py (not copied from the Fable plan-review "
        "summary in data/OUTSTANDING.md, though they match it to displayed precision); "
        "#489 (prior, independent confirmation that these same 4 moons are mass-deficient "
        "for ballistic-bend flyby usage in a different pairing).",
        run={
            "date": "2026-07-12",
            "cores": 1,
            "git_sha": git_sha,
            "wall_s": 0.0,
            "task": 571,
            "note": "analytic pre-sweep gate evaluation, not a numerical sweep -- see verdict",
        },
    )
    append_empty_region(DEFAULT_EMPTY_REGIONS_PATH, report)
    print(f"[571] appended {region_id}")

print(f"[571] DONE -- {len(SMALL_MOONS)} analytically-empty Titan-small-moon pairs stamped")
