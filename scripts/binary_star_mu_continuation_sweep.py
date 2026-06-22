"""Task #315 — binary-star cycler μ-continuation sweep.

Track the Earth-Moon (3,1) and (1,1) prograde cycler families up the
mass-parameter (μ) branch via pseudo-arclength continuation, to test whether the
stable prograde binary-star cyclers depicted in Roberts-Tsoukkas & Ross 2026
(μ = 0.3, 0.5) lie on the *continuous* branch of the Earth-Moon families
(complements #255, which seeded the figure ICs directly).

When a branch folds before reaching its target μ, records a first-class,
method-versioned negative via :func:`append_empty_region` — capturing the
**maximum μ the branch actually reaches** (the physically meaningful number),
not the arclength endpoint where the walk happens to underflow.
"""

import logging
import subprocess

import numpy as np

from cyclerfinder.core.cr3bp import CR3BPSystem
from cyclerfinder.data.empty_regions import (
    DEFAULT_EMPTY_REGIONS_PATH,
    EmptyRegionReport,
    append_empty_region,
)
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.search.binary_star_search import winding_topology
from cyclerfinder.search.cr3bp_periodic import _xaxis_crossings, correct_symmetric_fixed_jacobi
from cyclerfinder.search.mu_continuation import continue_in_mu, scan_c_family_at_mu

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

EARTH_MOON_MU = 1.2150584270572e-2

# Continuation parameters — recorded with every negative so the empty-region is
# conditional on the method, per the negative-results-registry convention.
CONT_PARAMS = {
    "ds0": 8e-3,
    "ds_max": 1e-2,
    "ds_min": 1e-6,
    "period_jump_frac": 0.3,
    "max_steps": 5000,
}

SEEDS = [
    {
        "id": "ross-rt-em-cycler-31-2025",
        "k1": 3,
        "k2": 1,
        "c_start": 3.161784147013429,
        "period_start": 14.78849241668140,
        "x0_guess": -0.3209891696,
        "ydot0_sign": -1.0,
        "mu_target": 0.3,
    },
    {
        "id": "ross-rt-em-cycler-11-2025",
        "k1": 1,
        "k2": 1,
        "c_start": 3.151175879508174,
        "period_start": 10.29206921007976,
        "x0_guess": -0.7682140805,
        "ydot0_sign": -1.0,
        "mu_target": 0.5,
    },
]


def _git_sha():
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


def build_fold_report(seed_info, *, mu_max, mu_min, n_members, stop_reason, git_sha):
    """A first-class empty-region report for a branch that folds before its target.

    Pure (no integration) so the verified numbers from a completed run can be
    re-serialised without re-walking the (expensive) continuation.
    """
    k1, k2 = seed_info["k1"], seed_info["k2"]
    target = seed_info["mu_target"]
    return EmptyRegionReport(
        region_id=f"binary-star-mu-continuation-em-{k1}{k2}-to-mu{target}",
        family="binary-star cycler (μ-continuation from Earth-Moon)",
        centre="Earth-Moon → binary-star μ",
        topologies=(
            {"k1": k1, "k2": k2, "prograde": True, "period_tu": seed_info["period_start"]},
        ),
        method_capability=MethodCapability(
            genome=f"CR3BP symmetric prograde cycler, winding ({k1},{k2})",
            corrector="pseudo-arclength continuation + Radau full-period closure",
            capability_tags=frozenset(
                {"mu-continuation", "pseudo-arclength", "radau-closure", "cr3bp", "planar"}
            ),
            git_sha=git_sha,
        ),
        search_extent={
            "points_total": int(n_members),
            "mu_start": EARTH_MOON_MU,
            "mu_target": target,
            "mu_max_reached": mu_max,
            "mu_min_reached": mu_min,
            "ephem_model": "cr3bp",
            **CONT_PARAMS,
        },
        prune_gates=(
            f"period_jump_frac={CONT_PARAMS['period_jump_frac']} (branch-switch/fold guard)",
            f"step_underflow at ds_min={CONT_PARAMS['ds_min']}",
        ),
        result={
            "target_reached": False,
            "mu_max_reached": mu_max,
            "mu_min_reached": mu_min,
            "n_members": int(n_members),
            "stop_reason": stop_reason,
        },
        verdict="EMPTY — family folds before reaching target μ",
        interpretation=(
            f"The Earth-Moon ({k1},{k2}) prograde cycler family, pseudo-arclength-continued in "
            f"μ from the seed (μ={EARTH_MOON_MU:.6f}), climbs to μ_max={mu_max:.6f} then folds; "
            f"it does NOT connect to the Roberts-Tsoukkas/Ross 2026 binary-star target μ={target} "
            f"on the continuous CR3BP branch. Arclength walk ran to μ={mu_min:.6f} before "
            f"{stop_reason}. Confirms the #255 direct-seed negative structurally."
        ),
        source_anchors="roberts-tsoukkas-ross-2026 (AAS 25-621 / journal); complements #255",
        run={"task": 315, "script": "scripts/binary_star_mu_continuation_sweep.py"},
    )


def sweep_seed(seed_info, git_sha):
    """Continue one seed toward its target μ; return a fold report dict or None."""
    logging.info("=== %s -> mu = %s ===", seed_info["id"], seed_info["mu_target"])
    system = CR3BPSystem(
        mu=EARTH_MOON_MU, primary="Earth", secondary="Moon", l_km=384400.0, t_s=375699.8
    )
    seed = correct_symmetric_fixed_jacobi(
        system=system,
        x0_guess=seed_info["x0_guess"],
        jacobi=seed_info["c_start"],
        period_guess=seed_info["period_start"],
        ydot0_sign=seed_info["ydot0_sign"],
        half_crossings=None,
    )
    if not seed.converged:
        logging.error("Seed %s failed to converge; skipping.", seed_info["id"])
        return None

    state0 = np.array([seed.x0, 0.0, 0.0, 0.0, seed.ydot0, 0.0])
    times, _ = _xaxis_crossings(
        system, state0, seed.t_half * 1.1, with_stm=False, rtol=1e-12, atol=1e-12
    )
    idx = int(np.argmin(np.abs(times - seed.t_half))) + 1
    logging.info("Seed recovered: x0=%.6f period=%.4f half_crossings=%d", seed.x0, seed.period, idx)

    branch = continue_in_mu(
        seed=seed,
        mu_start=EARTH_MOON_MU,
        half_crossings=idx,
        ydot0_sign=seed_info["ydot0_sign"],
        mu_target=seed_info["mu_target"],
        ds0=CONT_PARAMS["ds0"],
        ds_max=CONT_PARAMS["ds_max"],
        ds_min=CONT_PARAMS["ds_min"],
        max_steps=CONT_PARAMS["max_steps"],
        period_jump_frac=CONT_PARAMS["period_jump_frac"],
    )
    if not branch.members:
        logging.error("No members on branch for %s.", seed_info["id"])
        return None

    mus = np.array([m.mu for m in branch.members])
    mu_max, mu_min = float(mus.max()), float(mus.min())
    logging.info(
        "Branch: n=%d stop=%s mu_max=%.6f mu_min=%.6f", len(mus), branch.stop_reason, mu_max, mu_min
    )

    if mu_max < seed_info["mu_target"] - 1e-6:
        return build_fold_report(
            seed_info,
            mu_max=mu_max,
            mu_min=mu_min,
            n_members=len(mus),
            stop_reason=branch.stop_reason,
            git_sha=git_sha,
        )

    # Target μ reached — check topology + stability (the positive path).
    final = max(branch.members, key=lambda m: m.mu)
    topo = winding_topology(final.mu, final.state0, final.period)
    logging.info("At mu_max: topology k1=%d k2=%d prograde=%s", topo.k1, topo.k2, topo.prograde)
    if topo.k1 == seed_info["k1"] and topo.k2 == seed_info["k2"] and topo.prograde:
        if final.stable:
            logging.info("STABLE member at target mu (nu=%s) — novelty check required.", final.nu)
        else:
            logging.info("Unstable at target mu; scanning C-family...")
            members = scan_c_family_at_mu(
                mu=seed_info["mu_target"],
                x0_guess=final.x0,
                c_center=final.jacobi,
                period_guess=final.period,
                half_crossings=idx,
                ydot0_sign=seed_info["ydot0_sign"],
                dc=5e-4,
                n_each=20,
            )
            logging.info(
                "C-scan stable members: %d (novelty check if > 0).", sum(m.stable for m in members)
            )
    else:
        logging.info("Reached target mu but topology changed to (%d,%d).", topo.k1, topo.k2)
    return None


def main():
    git_sha = _git_sha()
    reports = [r for seed in SEEDS if (r := sweep_seed(seed, git_sha)) is not None]
    for report in reports:
        append_empty_region(DEFAULT_EMPTY_REGIONS_PATH, report)
    logging.info(
        "Appended %d empty-region report(s) to %s", len(reports), DEFAULT_EMPTY_REGIONS_PATH
    )


if __name__ == "__main__":
    main()
