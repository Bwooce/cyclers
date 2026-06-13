#!/usr/bin/env python3
"""Per-row multi-arc closure driver (#248) — a compact canonical table.

Runs the multi-arc closure harness
(:mod:`cyclerfinder.search.multiarc_closure`) over one or more catalogue rows and
prints the canonical ``(max_residual_kms, converged)`` per row. This is the
coordinator's entry point for the bounded convergence campaign: vary ``--n-starts``
and ``--rows`` here; the harness itself is the tested module.

NO catalogue writeback — a converged optimisation is a SILVER candidate for the
V0-V5 gauntlet, never a validated cycler. Evidence only.

The default row set is ``mcconaghy-2006-em-k2`` (the closest E-E-M-M row) plus the
six russell descriptor rows the design names for the follow-up.

Examples
--------
    # closest row, single start (the #248 proof-of-life)
    uv run python scripts/multiarc_closure_run.py --rows mcconaghy-2006-em-k2 --n-starts 1

    # the full default set, more starts (the coordinator's bounded campaign)
    uv run python scripts/multiarc_closure_run.py --n-starts 25
"""

from __future__ import annotations

import argparse
import time
import warnings

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.data.catalog import load_catalog
from cyclerfinder.search.multiarc_closure import close_multiarc_row

_DEFAULT_ROWS = (
    "mcconaghy-2006-em-k2",
    "russell-ch4-6.44Gg3",
    "russell-ch4-9.353Gg2",
    "russell-ch4-3.78Gg3",
    "russell-ch4-3.64gGg3",
    "russell-ch4-9.94Gg3",
    "russell-ch4-5.30ggF3",
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--rows",
        nargs="+",
        default=list(_DEFAULT_ROWS),
        help="catalogue row ids to close (default: the #248 follow-up set)",
    )
    parser.add_argument(
        "--n-starts",
        type=int,
        default=1,
        help="multi-start budget per row (distinct discrete seeds solved)",
    )
    parser.add_argument(
        "--gradient",
        default="lambert",
        choices=("lambert", "fbs-analytic"),
        help="corrector lane (default: lambert, the convergence reference)",
    )
    parser.add_argument("--tol-kms", type=float, default=0.1, help="convergence gate (km/s)")
    parser.add_argument(
        "--n-resonant",
        type=int,
        default=4,
        help="discrete resonant returns enumerated per resonant leg",
    )
    args = parser.parse_args()

    warnings.filterwarnings("ignore")
    cat = load_catalog()
    eph = Ephemeris("astropy")

    print(
        f"{'row':<26} {'seq':<12} {'best_res_kms':>12} {'conv':>5} "
        f"{'starts':>6} {'seeds':>6} {'s':>6}"
    )
    for rid in args.rows:
        entry = cat.by_id.get(rid)
        if entry is None:
            print(f"{rid:<26} <not in catalogue>")
            continue
        t0 = time.perf_counter()
        report = close_multiarc_row(
            entry.raw,
            eph,
            n_starts=args.n_starts,
            gradient=args.gradient,
            tol_kms=args.tol_kms,
            n_resonant=args.n_resonant,
        )
        dt = time.perf_counter() - t0
        seq = "-".join(report.sequence) if report.sequence else "<no-seed>"
        print(
            f"{report.row_id:<26} {seq:<12} {report.best_max_residual_kms:>12.4f} "
            f"{report.converged!s:>5} {report.n_starts_run:>6} "
            f"{report.n_seeds_available:>6} {dt:>6.1f}"
        )


if __name__ == "__main__":
    main()
