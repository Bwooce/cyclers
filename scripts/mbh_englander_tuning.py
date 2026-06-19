"""
Multi-dimensional MBH tuning script over the Englander & Englander 2014 parameters.
Evaluates Cauchy vs Pareto across 3 benchmark free-return genomes over multiple RNG seeds.
"""

from __future__ import annotations

import math

import numpy as np

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.free_return import _residuals
from cyclerfinder.search.mbh import make_free_return_step, mbh

DAY_S = 86400.0

# 3 Benchmark Free-Return Genomes
# [a_au, e, period_sec]
BENCHMARKS = {
    "mcconaghy-2006-em-k2": {
        "a_au": 1.30,
        "e": 0.257,
        "period_years": 4.27,
    },
    "niehoff-visit1": {
        "a_au": 1.17,
        "e": 0.193,
        "period_years": 14.95,
    },
    "niehoff-visit2": {
        "a_au": 1.52,
        "e": 0.34,
        "period_years": 14.95,
    },
}

RNG_SEEDS = list(range(1, 9))  # 8 seeds
N_HOPS = 100  # Give it a bit more hops to show Pareto long tails

# Unfrozen scale configuration for [a, e, t0]
# We use relative scale for a, e and absolute scale for t0
RELATIVE_SCALE = np.array([0.05, 0.05, 0.0])  # 5% relative for a and e
ABSOLUTE_SCALE = np.array([np.nan, np.nan, 8.0 * DAY_S])  # 8 days absolute for t0

# Mis-seed offsets
MISSEED_OFFSET_A = 0.05
MISSEED_OFFSET_E = 0.05
MISSEED_OFFSET_T0 = 40.0 * DAY_S


def _best_phase_t0(a: float, e: float, period_sec: float, ephem: Ephemeris) -> float:
    best_t0, best_res = 0.0, float("inf")
    for frac in np.linspace(0.0, 1.0, 500, endpoint=False):
        t0 = float(frac) * period_sec
        res = _residuals(
            np.array([a, e, t0]),
            period_days=period_sec / DAY_S,
            ephem=ephem,
            bodies=("E", "M"),
            mu=132712440018.0,
        )
        m = max(abs(r) for r in res)
        if m < best_res:
            best_res, best_t0 = m, t0
    return best_t0


def _hops_to_recovery(result) -> int | None:
    for i, b in enumerate(result.best_history):
        if np.isfinite(b) and b < 0.1:
            return i
    return None


def main() -> None:
    print("MBH Englander Tuning Benchmark")
    print("=" * 100)

    ephem = Ephemeris("circular")

    # Grid: cauchy vs pareto with various alphas
    configs = [
        {"dist": "cauchy", "alpha": float("nan")},
        {"dist": "pareto", "alpha": 1.01},
        {"dist": "pareto", "alpha": 1.05},
        {"dist": "pareto", "alpha": 1.08},
        {"dist": "pareto", "alpha": 1.12},
    ]

    header = f"{'Genome':<25} {'Dist':<10} {'Alpha':<8} {'Successes/8':<15} {'Mean Hops':<12}"
    print(header)
    print("-" * len(header))

    for name, params in BENCHMARKS.items():
        period_sec = params["period_years"] * 365.25 * DAY_S
        a_true = params["a_au"]
        e_true = params["e"]

        step_fn = make_free_return_step(period_sec=period_sec, ephem=ephem)
        best_t0 = _best_phase_t0(a_true, e_true, period_sec, ephem)

        misseed = np.array(
            [
                params["a_au"] + MISSEED_OFFSET_A,
                params["e"] + MISSEED_OFFSET_E,
                best_t0 - MISSEED_OFFSET_T0,
            ]
        )

        for config in configs:
            successes = 0
            hops_list = []

            for seed in RNG_SEEDS:
                result = mbh(
                    step_fn,
                    misseed,
                    n_hops=N_HOPS,
                    rng_seed=seed,
                    perturbation=config["dist"],
                    perturbation_alpha=config["alpha"] if not math.isnan(config["alpha"]) else 1.08,
                    perturbation_scale=RELATIVE_SCALE,
                    perturbation_absolute_scale=ABSOLUTE_SCALE,
                )

                hr = _hops_to_recovery(result)
                if hr is not None:
                    successes += 1
                    hops_list.append(hr)

            mean_hops = np.mean(hops_list) if hops_list else float("nan")
            print(
                f"{name:<25} {config['dist']:<10} {config['alpha']:<8.2f} "
                f"{successes:>11}     {mean_hops:>9.1f}"
            )

    print(
        "\nEvaluating Algorithm 1 (restart_bounds=15) on hardest genome "
        "(mcconaghy-2006-em-k2, pareto 1.08):"
    )
    # Evaluate Algorithm 1 on one of the genomes
    name = "mcconaghy-2006-em-k2"
    params = BENCHMARKS[name]
    period_sec = params["period_years"] * 365.25 * DAY_S
    step_fn = make_free_return_step(period_sec=period_sec, ephem=ephem)
    best_t0 = _best_phase_t0(params["a_au"], params["e"], period_sec, ephem)
    misseed = np.array(
        [
            params["a_au"] + MISSEED_OFFSET_A,
            params["e"] + MISSEED_OFFSET_E,
            best_t0 - MISSEED_OFFSET_T0,
        ]
    )

    restart_bounds = (np.array([1.0, 0.01, 0.0]), np.array([2.0, 0.6, period_sec]))

    for restart in [None, 15]:
        successes = 0
        hops_list = []
        for seed in RNG_SEEDS:
            result = mbh(
                step_fn,
                misseed,
                n_hops=N_HOPS,
                rng_seed=seed,
                perturbation="pareto",
                perturbation_alpha=1.08,
                perturbation_scale=RELATIVE_SCALE,
                perturbation_absolute_scale=ABSOLUTE_SCALE,
                stop_after_stall=restart,
                restart_bounds=restart_bounds if restart is not None else None,
            )
            hr = _hops_to_recovery(result)
            if hr is not None:
                successes += 1
                hops_list.append(hr)
        mean_hops = np.mean(hops_list) if hops_list else float("nan")
        print(
            f"restart_bounds={restart!s:<4} -> successes: {successes}/8, mean_hops: {mean_hops:.1f}"
        )


if __name__ == "__main__":
    main()
