# AP_FLAKE8_CLEAN
"""Typed argparse CLI for cyclerfinder — the user-facing M8 surface.

A thin driver over the already-shipped enumerator / optimiser / discover /
gauntlet stack. No new physics: every subcommand is a shell over an existing
function. The default body set is the M8 VEM anchor (spec §8 line 152).

Honesty boundary (spec §11.3/§17): `solve`/`report` separate *sourced* facts
(catalogue period / sequence / sourced V∞) from *computed* results (our
optimiser's V∞ / ΔV / closure residual). A computed value is never presented as
a sourced anchor.

Exit-code contract (module constants below): 0 ok · 2 usage error · 3
missing-viz-extra · 4 not-implemented (interim) · 5 no-candidates / empty-ledger.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import sys
from collections.abc import Callable, Sequence
from typing import Any

from cyclerfinder import __version__

# ---------------------------------------------------------------------------
# Exit-code contract
# ---------------------------------------------------------------------------

EXIT_OK = 0
EXIT_USAGE = 2
EXIT_MISSING_VIZ = 3
EXIT_NOT_IMPLEMENTED = 4
EXIT_NO_CANDIDATES = 5


# ---------------------------------------------------------------------------
# Per-subcommand argument builders
# ---------------------------------------------------------------------------


def _add_enumerate_args(sub: argparse.ArgumentParser) -> None:
    sub.add_argument("--bodies", default="V,E,M", help="comma-separated body codes")
    sub.add_argument("--l-max", type=int, default=6, help="max encounters in the sequence")
    sub.add_argument("--k-max", type=int, default=3, help="max period in synodic/beat multiples")
    sub.add_argument("--n-max", type=int, default=0, help="max revolutions per leg")
    sub.add_argument("--branch", default="single", help="comma-list: single|low|high")
    sub.add_argument("--vinf-cap", type=float, default=7.0, help="km/s feasibility cap")
    sub.add_argument("--feasible-only", action="store_true", help="only emit feasible cells")
    sub.add_argument("--period", default="beat", help="'beat' (basis=None) or an anchor pair 'E-M'")
    sub.add_argument("--format", default="table", choices=("table", "json", "csv"))
    sub.add_argument("--limit", type=int, default=None, help="cap rows emitted")


def _add_solve_args(sub: argparse.ArgumentParser) -> None:
    sub.add_argument("--cell-id", default=None, help="a Cell.id to optimise")
    sub.add_argument("--bodies", default=None, help="comma-separated body codes")
    sub.add_argument("--sequence", default=None, help="dash-separated flyby sequence, e.g. E-M-E")
    sub.add_argument("--k", type=int, default=None, help="period in synodic multiples")
    sub.add_argument("--period-basis", default=None, help="anchor pair 'A-B' for >=3-body cells")
    sub.add_argument("--revs", default=None, help="comma-separated per-leg revolutions")
    sub.add_argument("--branch", default=None, help="comma-separated per-leg branches")
    sub.add_argument("--fidelity", default="idealized", choices=("idealized", "ephemeris"))
    sub.add_argument("--mode", default="maintenance", choices=("maintenance", "ballistic"))
    sub.add_argument("--vinf-cap", type=float, default=7.0)
    sub.add_argument("--n-starts", type=int, default=5)
    sub.add_argument("--seed", type=int, default=0)
    sub.add_argument("--no-de", action="store_true", help="disable differential evolution pass")
    sub.add_argument("--priority-date", default=None, help="ephemeris epoch resolution, ISO date")
    sub.add_argument("--vinf-targets", default=None, help="phase-match, e.g. E=5.65,M=3.05")
    sub.add_argument("--format", default="table", choices=("table", "json"))


def _add_discover_args(sub: argparse.ArgumentParser) -> None:
    sub.add_argument("--bodies", default="V,E,M", help="comma-separated body codes")
    sub.add_argument("--k", type=int, default=3, help="period in synodic multiples")
    sub.add_argument("--vinf-cap", type=float, default=7.0)
    sub.add_argument("--ledger", required=True, help="path to the JSONL ledger to write")
    sub.add_argument("--l-max", type=int, default=4)
    sub.add_argument("--n-max", type=int, default=0)
    sub.add_argument("--branch", default="single", help="comma-list: single|low|high")
    sub.add_argument("--max-cells", type=int, default=None)
    sub.add_argument("--fidelity", default="idealized", choices=("idealized", "ephemeris"))
    sub.add_argument("--enable-v3", action="store_true", help="ballistic-closure V3 gate (M-ED)")
    sub.add_argument("--priority-date", default=None)
    sub.add_argument("--vinf-targets", default=None, help="e.g. E=5.65,M=3.05")
    sub.add_argument("--n-starts", type=int, default=5)
    sub.add_argument("--seed", type=int, default=0)
    sub.add_argument("--no-de", action="store_true")
    sub.add_argument("--format", default="table", choices=("table", "json"))


def _add_report_args(sub: argparse.ArgumentParser) -> None:
    sub.add_argument("--ledger", required=True, help="path to the JSONL ledger to read")
    sub.add_argument("--out", required=True, help="output path stem; writes .md and/or .json")
    sub.add_argument("--format", default="both", choices=("md", "json", "both"))
    sub.add_argument("--with-verdicts", action="store_true", help="attach run_gauntlet tiers")


def _add_viz_args(sub: argparse.ArgumentParser) -> None:
    viz_sub = sub.add_subparsers(dest="viz_kind")
    pork = viz_sub.add_parser("porkchop", help="epoch x ToF V-inf grid")
    pork.add_argument("--bodies", default="E,M", help="two body codes 'A,B'")
    pork.add_argument("--tof-min", type=float, default=100.0)
    pork.add_argument("--tof-max", type=float, default=400.0)
    pork.add_argument("--epoch-range", default="2032-01-01:2039-01-01")
    pork.add_argument("--fidelity", default="idealized", choices=("idealized", "ephemeris"))
    pork.add_argument("--out", required=True)

    traj = viz_sub.add_parser("trajectory", help="heliocentric XY of a built cycler")
    traj.add_argument("--cell-id", default=None)
    traj.add_argument("--bodies", default=None)
    traj.add_argument("--sequence", default=None)
    traj.add_argument("--k", type=int, default=None)
    traj.add_argument("--revs", default=None)
    traj.add_argument("--branch", default=None)
    traj.add_argument("--period-basis", default=None)
    traj.add_argument("--fidelity", default="idealized", choices=("idealized", "ephemeris"))
    traj.add_argument("--vinf-cap", type=float, default=7.0)
    traj.add_argument("--n-starts", type=int, default=5)
    traj.add_argument("--seed", type=int, default=0)
    traj.add_argument("--no-de", action="store_true")
    traj.add_argument("--out", required=True)

    beat = viz_sub.add_parser("beat", help="synodic beat-alignment diagram")
    beat.add_argument("--bodies", default="V,E,M")
    beat.add_argument("--out", required=True)


# ---------------------------------------------------------------------------
# Parser construction
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    """Build the top-level argparse parser with all five subcommands."""
    parser = argparse.ArgumentParser(
        prog="cyclerfinder",
        description="Find, rank, and verify planetary cycler trajectories.",
    )
    parser.add_argument("--version", action="version", version=f"cyclerfinder {__version__}")
    parser.add_argument("-v", "--verbose", action="count", default=0, help="increase verbosity")

    subparsers = parser.add_subparsers(dest="command", metavar="<subcommand>")
    _add_enumerate_args(subparsers.add_parser("enumerate", help="list cells + feasibility"))
    _add_solve_args(subparsers.add_parser("solve", help="optimise a single cell"))
    _add_discover_args(subparsers.add_parser("discover", help="ledger-backed discover() loop"))
    _add_report_args(subparsers.add_parser("report", help="campaign report from a ledger"))
    _add_viz_args(subparsers.add_parser("viz", help="visualisations (requires the [viz] extra)"))
    return parser


# ---------------------------------------------------------------------------
# Shared parsing / emit helpers
# ---------------------------------------------------------------------------


def _parse_bodies(spec: str, parser: argparse.ArgumentParser) -> tuple[str, ...]:
    """Parse a comma-separated body-code list, validating against constants.

    Unknown codes call ``parser.error`` (argparse exits 2).
    """
    from cyclerfinder.core.constants import SUPPORTED_BODIES

    bodies = tuple(b.strip() for b in spec.split(",") if b.strip())
    if not bodies:
        parser.error("--bodies must list at least one body code")
    unknown = [b for b in bodies if b not in SUPPORTED_BODIES]
    if unknown:
        parser.error(
            f"unknown body code(s) {','.join(unknown)}; supported: {','.join(SUPPORTED_BODIES)}"
        )
    return bodies


def _parse_period_basis(
    spec: str | None, parser: argparse.ArgumentParser
) -> tuple[str, str] | None:
    """Resolve --period / --period-basis: ``'beat'``/None → None; ``'A-B'`` → pair."""
    if spec is None or spec == "beat":
        return None
    from cyclerfinder.core.constants import SUPPORTED_BODIES

    parts = spec.split("-")
    if len(parts) != 2 or any(p not in SUPPORTED_BODIES for p in parts):
        parser.error(f"--period anchor must be 'A-B' of supported bodies; got {spec!r}")
    return (parts[0], parts[1])


def _emit(rows: list[dict[str, Any]], fmt: str, columns: Sequence[str]) -> None:
    """Print ``rows`` in the requested format to stdout."""
    if fmt == "json":
        print(json.dumps(rows, indent=2, sort_keys=True))
        return
    if fmt == "csv":
        buf = io.StringIO()
        writer = csv.DictWriter(buf, fieldnames=list(columns), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
        print(buf.getvalue(), end="")
        return
    # table — aligned columns
    if not rows:
        print("(no rows)")
        return
    widths = {c: len(c) for c in columns}
    for row in rows:
        for c in columns:
            widths[c] = max(widths[c], len(str(row.get(c, ""))))
    header = "  ".join(c.ljust(widths[c]) for c in columns)
    print(header)
    print("  ".join("-" * widths[c] for c in columns))
    for row in rows:
        print("  ".join(str(row.get(c, "")).ljust(widths[c]) for c in columns))


# ---------------------------------------------------------------------------
# enumerate handler
# ---------------------------------------------------------------------------


def _handle_enumerate(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    from cyclerfinder.search.sequence import Cell, enumerate_cells, tisserand_feasible

    bodies = _parse_bodies(args.bodies, parser)
    basis = _parse_period_basis(args.period, parser)
    branch_set = tuple(b.strip() for b in args.branch.split(",") if b.strip())

    rows: list[dict[str, Any]] = []
    for cell in enumerate_cells(bodies, args.l_max, args.k_max, args.n_max, branch_set):
        if basis is not None:
            cell = Cell(
                bodies=cell.bodies,
                sequence=cell.sequence,
                period_k=cell.period_k,
                per_leg_revs=cell.per_leg_revs,
                per_leg_branch=cell.per_leg_branch,
                period_basis=basis,
            )
        feasible = tisserand_feasible(cell, vinf_cap=args.vinf_cap)
        if args.feasible_only and not feasible:
            continue
        rows.append(
            {
                "cell_id": cell.id,
                "bodies": "".join(cell.bodies),
                "sequence": "-".join(cell.sequence),
                "period_k": cell.period_k,
                "feasible": feasible,
            }
        )
        if args.limit is not None and len(rows) >= args.limit:
            break

    _emit(rows, args.format, ("cell_id", "bodies", "sequence", "period_k", "feasible"))
    return EXIT_OK


# ---------------------------------------------------------------------------
# solve helpers + handler
# ---------------------------------------------------------------------------

_LETTER_BRANCH: dict[str, str] = {"s": "single", "l": "low", "h": "high"}


def _parse_cell_id(cell_id: str, parser: argparse.ArgumentParser) -> Any:
    """Parse a ``Cell.id`` string back into a :class:`Cell` (the inverse of ``Cell.id``).

    Validated by a round-trip assertion against ``Cell.id`` (``sequence.py``).
    """
    from cyclerfinder.search.sequence import Cell

    parts = cell_id.split("|")
    if len(parts) not in (5, 6):
        parser.error(f"malformed --cell-id {cell_id!r}")
    bodyset, sequence_tok, k_tok, revs_tok, branch_tok = parts[:5]
    basis: tuple[str, str] | None = None
    if len(parts) == 6:
        basis_tok = parts[5]
        if not basis_tok.startswith("p") or len(basis_tok) != 3:
            parser.error(f"malformed period-basis token in --cell-id {cell_id!r}")
        basis = (basis_tok[1], basis_tok[2])
    if not k_tok.startswith("k") or not revs_tok.startswith("r") or not branch_tok.startswith("b"):
        parser.error(f"malformed --cell-id {cell_id!r}")
    try:
        bodies = tuple(bodyset)
        sequence = tuple(sequence_tok.split("-"))
        period_k = int(k_tok[1:])
        per_leg_revs = tuple(int(c) for c in revs_tok[1:])
        per_leg_branch = tuple(_LETTER_BRANCH[c] for c in branch_tok[1:])
    except (ValueError, KeyError):
        parser.error(f"malformed --cell-id {cell_id!r}")
    cell = Cell(
        bodies=bodies,
        sequence=sequence,
        period_k=period_k,
        per_leg_revs=per_leg_revs,
        per_leg_branch=per_leg_branch,
        period_basis=basis,
    )
    if cell.id != cell_id:
        parser.error(f"--cell-id {cell_id!r} did not round-trip (got {cell.id!r})")
    return cell


def _parse_vinf_targets(
    spec: str | None, parser: argparse.ArgumentParser
) -> dict[str, float] | None:
    """Parse ``E=5.65,M=3.05`` into a ``{body: kms}`` mapping."""
    if spec is None:
        return None
    out: dict[str, float] = {}
    for token in spec.split(","):
        token = token.strip()
        if not token:
            continue
        if "=" not in token:
            parser.error(f"--vinf-targets entry must be BODY=KMS; got {token!r}")
        body, _, value = token.partition("=")
        try:
            out[body.strip()] = float(value)
        except ValueError:
            parser.error(f"--vinf-targets value not a float in {token!r}")
    return out


def _cell_from_args(args: argparse.Namespace, parser: argparse.ArgumentParser) -> Any:
    """Build a :class:`Cell` from explicit flags, or parse it from ``--cell-id``.

    Explicit flags take precedence when both are supplied.
    """
    from cyclerfinder.search.sequence import Cell

    has_flags = args.bodies is not None and args.sequence is not None and args.k is not None
    if has_flags:
        bodies = _parse_bodies(args.bodies, parser)
        sequence = tuple(args.sequence.split("-"))
        revs_spec = args.revs if args.revs is not None else ",".join(["0"] * (len(sequence) - 1))
        per_leg_revs = tuple(int(r) for r in revs_spec.split(",") if r != "")
        branch_spec = (
            args.branch if args.branch is not None else ",".join(["single"] * (len(sequence) - 1))
        )
        per_leg_branch = tuple(b for b in branch_spec.split(",") if b != "")
        basis = _parse_period_basis(args.period_basis, parser)
        return Cell(
            bodies=bodies,
            sequence=sequence,
            period_k=args.k,
            per_leg_revs=per_leg_revs,
            per_leg_branch=per_leg_branch,
            period_basis=basis,
        )
    if args.cell_id is not None:
        return _parse_cell_id(args.cell_id, parser)
    parser.error("solve needs --cell-id OR --bodies/--sequence/--k")


def _result_to_dict(result: Any) -> dict[str, Any]:
    """Serialise an OptimisationResult, segregating optimiser outputs under `computed`.

    Golden discipline: every value here is computed by our optimiser, never a
    sourced anchor, so it all lives under the ``computed`` key.
    """
    score = result.best_score
    return {
        "cell_id": result.cell.id,
        "converged": bool(result.converged),
        "constraints_satisfied": bool(result.constraints_satisfied),
        "computed": {
            "closure_residual_kms": float(result.closure_residual_kms),
            "max_vinf_kms": float(score.max_vinf_kms),
            "total_maintenance_dv_kms": float(score.total_maintenance_dv_kms),
            "taxi_cost_kms": float(score.taxi_cost_kms),
            "period_error_yr": float(score.period_error_yr),
            "hard_constraints_pass": bool(score.hard_constraints_pass),
        },
    }


def _handle_solve(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    from cyclerfinder.core.ephemeris import Ephemeris
    from cyclerfinder.search.optimize import (
        optimise_cell_ephemeris,
        optimise_cell_idealized,
    )

    cell = _cell_from_args(args, parser)

    if args.fidelity == "idealized":
        result = optimise_cell_idealized(
            cell,
            Ephemeris(model="circular"),
            vinf_cap=args.vinf_cap,
            n_starts=args.n_starts,
            seed=args.seed,
            use_de=not args.no_de,
        )
    else:
        vinf_targets = _parse_vinf_targets(args.vinf_targets, parser)
        result = optimise_cell_ephemeris(
            cell,
            Ephemeris(model="astropy"),
            vinf_cap=args.vinf_cap,
            priority_date_iso=args.priority_date,
            vinf_targets_kms=vinf_targets,
            n_starts=args.n_starts,
            seed=args.seed,
            mode=args.mode,
        )

    payload = _result_to_dict(result)
    if args.format == "json":
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        flat = {"cell_id": payload["cell_id"], **payload["computed"]}
        _emit([flat], "table", list(flat.keys()))
    # `solve` produced a result: success. Convergence (constraints_satisfied) is
    # carried in the payload, never conflated with the run's exit code — non-
    # convergence is honest data, not a CLI failure (spec §11.3 honesty).
    return EXIT_OK


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------


def _stub_handler(name: str) -> Callable[[argparse.Namespace, argparse.ArgumentParser], int]:
    def _handler(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
        del args, parser
        print(f"cyclerfinder {name}: not implemented yet")
        return EXIT_NOT_IMPLEMENTED

    return _handler


_HANDLERS: dict[str, Callable[[argparse.Namespace, argparse.ArgumentParser], int]] = {
    "enumerate": _handle_enumerate,
    "solve": _handle_solve,
    "discover": _stub_handler("discover"),
    "report": _stub_handler("report"),
    "viz": _stub_handler("viz"),
}


def main(argv: Sequence[str] | None = None) -> int:
    """Parse ``argv`` and dispatch to the matching subcommand handler.

    Returns the handler's int exit code. ``--version`` raises ``SystemExit(0)``
    via argparse. No subcommand → usage to stderr, return 2.
    """
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_usage(file=sys.stderr)
        return EXIT_USAGE
    return _HANDLERS[args.command](args, parser)


if __name__ == "__main__":
    raise SystemExit(main())
