"""One-off local sweep validating the perturbation distribution + scale choice on
OUR Gate-2 free-return recovery (task #156, Englander & Englander 2014 upgrade).

NOT a pytest gate -- a recorded one-off whose output is pasted into the addendum
of docs/notes/2026-06-07-mbh-wrapper.md. Re-run with:

    uv run python scripts/mbh_pareto_sweep.py

Setup mirrors tests/search/test_mbh.py::test_free_return_mbh_recovers_sourced_
basin_from_misseed exactly: row mcconaghy-2006-em-k2, circular model, 40-day
off-phase mis-seed of the #137 free-return genome, rng_seed=6, a/e frozen, t0
perturbed ABSOLUTELY. The Gate-2 test pins gaussian + 8-day t0 step; here we
sweep {uniform, cauchy, pareto} x a small t0-step grid (default/10 .. default,
where the Gate-2 "default" absolute t0 step is 8 days) and, for pareto, the
excursion is alpha (held at the module default 1.08) with the t0 step as the
per-gene multiplier. We report hops-to-first-recovery and final feasibility per
cell.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import yaml  # type: ignore[import-untyped]

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.free_return import _residuals
from cyclerfinder.search.mbh import make_free_return_step, mbh

REPO_ROOT = Path(__file__).resolve().parents[1]
DAY_S = 86400.0

# SOURCED constraint side (Rogers 2012 Table 1) + V_inf anchor (McConaghy 2006).
_S1L1_A_AU = 1.30
_S1L1_E = 0.257
_SRC_VINF_E = 4.7
_SRC_VINF_M = 5.0

DISTRIBUTIONS = ("uniform", "cauchy", "pareto")
# Gate-2 default absolute t0 step is 8 days; grid spans default/10 .. default.
T0_STEP_DAYS = (0.8, 2.5, 8.0)
RNG_SEED = 6
N_HOPS = 60


def _row(rid: str) -> dict[str, Any]:
    rows = yaml.safe_load((REPO_ROOT / "data" / "catalogue.yaml").read_text())
    return next(r for r in rows if r["id"] == rid)


def _best_phase_t0(ephem: Any, period_sec: float) -> float:
    best_t0, best_res = 0.0, float("inf")
    for frac in np.linspace(0.0, 1.0, 2000, endpoint=False):
        t0 = float(frac) * period_sec
        res = _residuals(
            np.array([_S1L1_A_AU, _S1L1_E, t0]),
            period_days=period_sec / DAY_S,
            ephem=ephem,
            bodies=("E", "M"),
            mu=132712440018.0,
        )
        m = max(abs(r) for r in res)
        if m < best_res:
            best_res, best_t0 = m, t0
    return best_t0


def _hops_to_recovery(result: Any) -> int | None:
    """First hop index whose best-history reaches a converged residual."""
    for i, b in enumerate(result.best_history):
        if np.isfinite(b) and b < 0.1:
            return i
    return None


def main() -> None:
    row = _row("mcconaghy-2006-em-k2")
    period_sec = float(row["period"]["years"]) * 365.25 * DAY_S
    ephem = Ephemeris("circular")
    step = make_free_return_step(period_sec=period_sec, ephem=ephem)
    best_t0 = _best_phase_t0(ephem, period_sec)
    misseed = np.array([_S1L1_A_AU, _S1L1_E, best_t0 - 40.0 * DAY_S])

    print(f"row=mcconaghy-2006-em-k2  rng_seed={RNG_SEED}  n_hops={N_HOPS}  40d mis-seed")
    print("sourced anchor: V_inf E=4.7 M=5.0; ellipse a=1.30 e=0.257")
    print()
    header = (
        f"{'dist':>8} {'t0_step_d':>10} {'recovered':>10} {'hops2conv':>9} "
        f"{'obj_kms':>10} {'a_rec':>7} {'e_rec':>7} {'vinfE':>6} {'vinfM':>6} {'acc/att':>8}"
    )
    print(header)
    print("-" * len(header))

    for dist in DISTRIBUTIONS:
        for step_d in T0_STEP_DAYS:
            result = mbh(
                step,
                misseed,
                n_hops=N_HOPS,
                perturbation=dist,
                perturbation_scale=None,
                perturbation_absolute_scale=[0.0, 0.0, step_d * DAY_S],
                rng_seed=RNG_SEED,
            )
            converged = result.best_feasible and result.best_objective < 0.1
            h2r = _hops_to_recovery(result)
            a_rec = result.best_info.get("a_au", float("nan"))
            e_rec = result.best_info.get("e", float("nan"))
            vinf = result.best_info.get("vinf_kms", {})
            ve = vinf.get("E", float("nan"))
            vm = vinf.get("M", float("nan"))
            # SOURCED-basin recovery: not just converged, but in the sourced
            # ellipse + V_inf anchor (the Gate-2 pass criterion).
            on_anchor = (
                converged
                and abs(a_rec - _S1L1_A_AU) <= 0.03
                and abs(e_rec - _S1L1_E) <= 0.03
                and abs(vm - _SRC_VINF_M) <= 0.5
                and abs(ve - _SRC_VINF_E) <= 1.0
            )
            label = "ANCHOR" if on_anchor else ("conv" if converged else "no")
            print(
                f"{dist:>8} {step_d:>10.1f} {label:>10} "
                f"{(str(h2r) if h2r is not None else '-'):>9} {result.best_objective:>10.4f} "
                f"{a_rec:>7.3f} {e_rec:>7.3f} {ve:>6.2f} {vm:>6.2f} "
                f"{result.hops_accepted}/{result.hops_attempted:>3}"
            )


if __name__ == "__main__":
    main()
