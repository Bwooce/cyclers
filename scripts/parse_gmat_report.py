"""Parse a GMAT V4 ReportFile and apply the two-part V4 predicate (#171).

Consumes the text GMAT writes after a headless run of a script from
:mod:`scripts.gmat_v4_generate`, and decides V4-PASS. The predicate is two-part
(design §4):

1. **Convergence (Beeson AAS 15-278, primary):** every ``Target``/``Optimize`` block
   reported NLP convergence — defects driven below GMAT's feasibility tolerance,
   i.e. Jones AAS 17-577 §2.5's published continuity bar 1e-3 km / 1e-6 km/s. This
   is the SOURCED tolerance, not ours.
2. **Maintenance-ΔV band (self-declared ±5%, where a reference exists):** the summed
   per-flyby converged TCM magnitude reproduces OUR reference within ±5%
   (Aldrin 2.9138 km/s -> 2.768-3.060; S1L1 62 m/s -> 58.9-65.1). For the
   Mars-perturbed continuous arm there is NO prior reference number (it is the figure
   GMAT produces for the first time), so its gate is **convergence-only** and the
   produced ΔV becomes the recorded figure — you cannot band a number you do not have.

GOLDEN / HONESTY. 2.9138 km/s and 62 m/s are OUR computed values; GMAT (independent
codebase + ephemeris) is the external check, never an EXPECTED-from-source assertion.
The ±5% band is declared UP FRONT and never back-fit.

This parser does not invoke GMAT. Usage (after a manual headless GMAT run)::

    uv run python scripts/parse_gmat_report.py <report> --ref-dv 2.9138
    uv run python scripts/parse_gmat_report.py <report>            # convergence-only
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path

# A converged Target/Optimize block. GMAT prints variants of these on success.
_CONVERGED_RE = re.compile(
    r"converged|convergence achieved|targeting.*converged|optimization.*converged",
    re.IGNORECASE,
)
# An explicit non-convergence signal — any of these on a block fails the run.
_NOT_CONVERGED_RE = re.compile(
    r"did not converge|failed to converge|not converged|no convergence|exceeded.*max.*iter",
    re.IGNORECASE,
)
# A per-flyby TCM magnitude. GMAT ReportFile columns / Maneuver summaries print
# either "<name>.Magnitude = <val>" or "<name>.TotalDV = <val>" (km/s).
_TCM_RE = re.compile(
    r"(?:TCM_\w+\.(?:Magnitude|TotalDV)|Maneuver\.TotalDV)\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)"
)

DEFAULT_TOL_FRAC: float = 0.05  # the self-declared ±5% band (design §4)


@dataclass(frozen=True)
class ParsedReport:
    """What the parser extracts from a GMAT ReportFile."""

    converged: bool
    n_converged_blocks: int
    n_failed_blocks: int
    tcm_magnitudes_kms: tuple[float, ...]
    maintenance_dv_kms: float


def parse_convergence(report_text: str) -> bool:
    """True iff every Target/Optimize block converged and NONE failed (Beeson primary).

    A single explicit non-convergence line fails the whole run; absent any
    convergence evidence, returns False (no silent pass).
    """
    if _NOT_CONVERGED_RE.search(report_text):
        return False
    return bool(_CONVERGED_RE.search(report_text))


def parse_maintenance_dv(report_text: str) -> float:
    """Sum the per-flyby converged TCM magnitudes (km/s) over one cycle/synodic period."""
    return float(sum(float(m) for m in _TCM_RE.findall(report_text)))


def parse_report(report_text: str) -> ParsedReport:
    """Extract convergence + per-flyby TCM magnitudes + their sum."""
    n_ok = len(_CONVERGED_RE.findall(report_text))
    n_fail = len(_NOT_CONVERGED_RE.findall(report_text))
    tcms = tuple(float(m) for m in _TCM_RE.findall(report_text))
    return ParsedReport(
        converged=parse_convergence(report_text),
        n_converged_blocks=n_ok,
        n_failed_blocks=n_fail,
        tcm_magnitudes_kms=tcms,
        maintenance_dv_kms=float(sum(tcms)),
    )


def v4_pass(
    gmat_dv: float,
    ref_dv: float | None,
    *,
    converged: bool,
    tol_frac: float = DEFAULT_TOL_FRAC,
) -> bool:
    """The two-part V4 predicate.

    V4-PASS iff:

    * ``converged`` is True (the flyby B-plane Targets closed in GMAT's high-fidelity
      model — Beeson's "validated"); AND
    * EITHER ``ref_dv is None`` (convergence-only: no prior reference, the produced
      ``gmat_dv`` becomes the recorded figure — the Mars-perturbed continuous arm),
      OR ``|gmat_dv - ref_dv| <= tol_frac * ref_dv`` (the maintenance ΔV reproduces
      OUR reference within the self-declared ±5% band).

    ``ref_dv`` (Aldrin 2.9138 km/s; S1L1 62 m/s) is OUR computed value under GMAT's
    external check, never an EXPECTED-from-source assertion; the ±5% band is declared
    up front and never back-fit.
    """
    if not converged:
        return False
    if ref_dv is None:
        return True
    return abs(gmat_dv - ref_dv) <= tol_frac * abs(ref_dv)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Parse a GMAT V4 report and apply the V4 predicate."
    )
    parser.add_argument("report", type=Path, help="GMAT ReportFile path")
    parser.add_argument(
        "--ref-dv",
        type=float,
        default=None,
        help="OUR reference maintenance dV (km/s); omit for convergence-only.",
    )
    parser.add_argument(
        "--tol-frac", type=float, default=DEFAULT_TOL_FRAC, help="band fraction (default 0.05)"
    )
    args = parser.parse_args()

    text = args.report.read_text(encoding="utf-8")
    parsed = parse_report(text)
    verdict = v4_pass(
        parsed.maintenance_dv_kms,
        args.ref_dv,
        converged=parsed.converged,
        tol_frac=args.tol_frac,
    )
    print(
        f"converged: {parsed.converged} "
        f"({parsed.n_converged_blocks} ok, {parsed.n_failed_blocks} failed)"
    )
    print(f"per-flyby TCM (km/s): {parsed.tcm_magnitudes_kms}")
    print(f"summed maintenance dV (km/s): {parsed.maintenance_dv_kms:.6f}")
    if args.ref_dv is None:
        print("reference: convergence-only (no prior reference dV)")
    else:
        print(f"reference (OUR value): {args.ref_dv:.4f} km/s, band +-{args.tol_frac:.0%}")
    print(f"V4-PASS: {verdict}")


if __name__ == "__main__":
    main()


__all__ = [
    "DEFAULT_TOL_FRAC",
    "ParsedReport",
    "parse_convergence",
    "parse_maintenance_dv",
    "parse_report",
    "v4_pass",
]
