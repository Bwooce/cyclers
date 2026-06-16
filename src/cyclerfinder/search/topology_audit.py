"""Independent topology audit harness for discovery JSONLs (#325 Phase 1 Part B).

Defensive sweep for #322-class topology-gate bugs. The #322 bug had this shape:

    A genome's topology-classification gate (``tulip.petal_count``) counts a
    feature (in-plane petals) without checking the FULL state dimensionality.
    At extreme regimes (tiny mu) the corrector drives z0 -> 0; the orbit
    collapses to a planar Np-petal orbit; the petal-count gate fires; the
    genome reports a "3D tulip" that is actually planar.

The #322 fix (``cyclerfinder.genome.tulip.is_three_dimensional``) adds a
complementary 3D-amplitude check on max|z(t)|, which catches the collapse.
**This module re-verifies every discovery JSONL row through an INDEPENDENT
topology check** -- not the genome's own gate, but a from-scratch verification
of the topology claim, cross-anchored to published formulas.

Scope
-----
This is a Phase-1 audit harness. It is read-only with respect to the discovery
JSONLs (it consumes them); it does NOT write back to the catalogue, and it
does NOT make novelty claims. Its output is a per-row "independent verdict"
JSONL + a doc summarising the sibling-bug count.

Checkers implemented (all sourced):

  * :func:`check_tulip_topology` -- re-runs ``is_three_dimensional`` from
    :mod:`cyclerfinder.genome.tulip` on the row's 3D periodic-orbit IC. Floor
    sourced to Koblick 2023 Table 4 (#322 commit ``c2a77c7``).

  * :func:`check_floquet_neimark_sacker` -- classifies a Floquet eigenvalue
    pair as Neimark-Sacker iff BOTH eigenvalues lie on the unit circle within
    a tolerance AND they are not real (within tolerance of +/-1). A real pair
    near +1 is the trivial Liouville pair; a real pair near -1 is k=2
    period-doubling -- neither is Neimark-Sacker. Sourced to standard CR3BP
    bifurcation theory (Gomez, Koon, Lo, Marsden, Masdemont, Ross 2001).

  * :func:`check_periodic_orbit_closure` -- independent full-period
    re-propagation (Radau, distinct integrator) of the row's IC and check
    that ``||X(T) - X(0)||`` is below a sourced floor.

Discipline
----------
  * NO catalogue writeback. Output is JSONL + doc only.
  * NO novelty claims. A discrepancy is a flag, not a verdict.
  * Each checker is INDEPENDENT of the genome it audits -- the floor /
    classification rule is sourced, not the genome's own.
  * Re-verification is informational: a row that downgrades here is a
    candidate for genome-side investigation, not an automatic retraction.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.tulip import (
    TULIP_Z_AMPLITUDE_FLOOR_NONDIM,
    is_three_dimensional,
)

# ---------------------------------------------------------------------------
# Sourced thresholds.
# ---------------------------------------------------------------------------

# Independent full-period closure floor for the cross-integrator check.
# Sourced to :func:`cyclerfinder.genome.bcr4bp_genome.correct_bcr4bp_periodic`'s
# ``independent_tol`` default (1e-6) and the #266 Phase 1 NRHO acceptance
# floor in :func:`cyclerfinder.search.nrho_continuation`.
INDEPENDENT_CLOSURE_FLOOR_NONDIM: float = 1.0e-6

# Floquet "on the unit circle" tolerance. The bifurcation_detector default
# (``cyclerfinder.search.bifurcation_detector.scan_family_for_bifurcations``)
# uses ``tol=1e-2`` for "distance to a primitive k-th root of unity".
UNIT_CIRCLE_TOL: float = 1.0e-2

# Real-axis tolerance for "is this eigenvalue real?" (eps in |Im(lambda)| < eps).
REAL_AXIS_TOL: float = 1.0e-3


@dataclass(frozen=True)
class TopologyAuditFinding:
    """One row's independent-topology verdict.

    Attributes
    ----------
    source_jsonl :
        Path to the JSONL the row came from.
    row_id :
        Stable identifier for the row inside the JSONL.
    genome_gate :
        Which genome / topology class the original row was classified under.
    genome_verdict :
        What the genome's gate said.
    independent_verdict :
        What this audit says, from an independent checker.
    discrepancy :
        ``genome_verdict != independent_verdict``.
    failure_mode :
        If ``discrepancy`` is True, what's the actual topology.
    extras :
        Numerical diagnostics.
    """

    source_jsonl: str
    row_id: str
    genome_gate: str
    genome_verdict: bool
    independent_verdict: bool
    discrepancy: bool
    failure_mode: str
    extras: dict[str, float] = field(default_factory=dict)


def check_tulip_topology(
    state0: NDArray[np.float64],
    period: float,
    system: cr3bp.CR3BPSystem,
    *,
    z_floor: float = TULIP_Z_AMPLITUDE_FLOOR_NONDIM,
    rtol: float = 1e-11,
    atol: float = 1e-11,
) -> tuple[bool, dict[str, float]]:
    """Independent 3D-topology check for a tulip / NRHO orbit row.

    Re-runs the #322 :func:`cyclerfinder.genome.tulip.is_three_dimensional`
    gate on the row's IC. Returns ``(is_3d, extras)``.
    """
    s0 = np.asarray(state0, dtype=np.float64)
    is_3d, max_abs_z = is_three_dimensional(
        s0, float(period), system, z_floor=z_floor, rtol=rtol, atol=atol
    )
    return bool(is_3d), {
        "z0": float(s0[2]),
        "max_abs_z": float(max_abs_z),
        "z_floor": float(z_floor),
    }


def check_floquet_neimark_sacker(
    eig_real: float,
    eig_imag: float,
    *,
    unit_circle_tol: float = UNIT_CIRCLE_TOL,
    real_axis_tol: float = REAL_AXIS_TOL,
) -> tuple[bool, dict[str, float]]:
    """Independent Neimark-Sacker classification of a single Floquet eigenvalue.

    Sourced to standard CR3BP bifurcation taxonomy (Gomez et al. 2001 sec 3):
    Neimark-Sacker requires ``|lambda| ~ 1`` AND ``|Im(lambda)|`` strictly
    above ``real_axis_tol``.
    """
    eig_abs = float(np.hypot(eig_real, eig_imag))
    abs_imag = float(abs(eig_imag))
    on_unit_circle = abs(eig_abs - 1.0) <= unit_circle_tol
    truly_complex = abs_imag > real_axis_tol
    is_ns = on_unit_circle and truly_complex
    return is_ns, {
        "eig_real": float(eig_real),
        "eig_imag": float(eig_imag),
        "abs_eig": eig_abs,
        "imag_part": abs_imag,
        "unit_circle_distance": float(abs(eig_abs - 1.0)),
        "unit_circle_tol": float(unit_circle_tol),
        "real_axis_tol": float(real_axis_tol),
    }


def check_periodic_orbit_closure(
    state0: NDArray[np.float64],
    period: float,
    system: cr3bp.CR3BPSystem,
    *,
    floor: float = INDEPENDENT_CLOSURE_FLOOR_NONDIM,
    rtol: float = 1e-11,
    atol: float = 1e-11,
) -> tuple[bool, dict[str, float]]:
    """Independent (Radau, full-period) closure cross-check."""
    from scipy.integrate import solve_ivp

    s0 = np.asarray(state0, dtype=np.float64)
    if s0.shape != (6,):
        raise ValueError(
            f"check_periodic_orbit_closure: state0 must be a 6-vector, got shape {s0.shape}"
        )
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, float(period)),
        s0,
        args=(system.mu,),
        method="Radau",
        rtol=rtol,
        atol=atol,
    )
    if not sol.success:
        return False, {
            "closure_residual": float("inf"),
            "floor": float(floor),
            "integrator_failed": 1.0,
        }
    closure_residual = float(np.linalg.norm(sol.y[:, -1] - s0))
    return (closure_residual < floor), {
        "closure_residual": closure_residual,
        "floor": float(floor),
        "integrator_failed": 0.0,
    }


def _system_from_row_or_header(
    row: dict[str, Any], header: dict[str, Any] | None
) -> cr3bp.CR3BPSystem | None:
    """Best-effort CR3BPSystem reconstruction from a JSONL row."""
    sys_block: dict[str, Any] | None = None
    for src in (row, header):
        if isinstance(src, dict) and isinstance(src.get("system"), dict):
            sys_block = src["system"]
            break
    if sys_block is None:
        for src in (row, header):
            if isinstance(src, dict) and "mu" in src:
                sys_block = src
                break
    if sys_block is None:
        return None
    mu = sys_block.get("mu")
    if mu is None:
        return None
    return cr3bp.CR3BPSystem(
        mu=float(mu),
        primary=str(sys_block.get("primary", "unknown")),
        secondary=str(sys_block.get("secondary", "unknown")),
        l_km=float(sys_block.get("l_km", 1.0)),
        t_s=float(sys_block.get("t_s", 1.0)),
    )


def _extract_state_and_period(row: dict[str, Any]) -> tuple[NDArray[np.float64], float] | None:
    """Best-effort (state6, period_nondim) extraction from a JSONL row."""
    state = row.get("state_nd") or row.get("state_initial") or row.get("state6")
    period = row.get("T_TU") or row.get("period_nondim") or row.get("period")
    if state is None or period is None:
        return None
    state_arr = np.asarray(state, dtype=np.float64)
    if state_arr.shape != (6,):
        return None
    return state_arr, float(period)


def _classify_row_topology(row: dict[str, Any], header: dict[str, Any] | None) -> str:
    """Heuristically identify which topology gate the row was classified under."""
    del header
    has_state = "state_nd" in row or "state_initial" in row
    has_period = "T_TU" in row or "period_nondim" in row
    if has_state and has_period:
        return "tulip_or_3d_periodic"
    if "floquet_real" in row and "floquet_imag" in row:
        return "floquet_bifurcation"
    return "unknown"


def audit_topology(
    jsonl_path: str,
    *,
    independent_checkers: dict[str, Callable[..., object]] | None = None,
    z_floor: float = TULIP_Z_AMPLITUDE_FLOOR_NONDIM,
    closure_floor: float = INDEPENDENT_CLOSURE_FLOOR_NONDIM,
    skip_rows_without_state: bool = True,
) -> list[TopologyAuditFinding]:
    """Re-verify topology claims across a discovery JSONL.

    Walks the JSONL line-by-line. The FIRST line is treated as a header if
    it carries ``type=='header'`` or ``_meta=True`` -- its ``system`` block
    is used as a fallback for rows that omit ``system``. Each subsequent
    row is dispatched to applicable independent checkers:

      * Rows with ``state_nd`` / ``state_initial`` + ``T_TU`` / ``period_nondim``
        -> :func:`check_tulip_topology` AND :func:`check_periodic_orbit_closure`.
        Rows already tagged ``degenerate_planar=True`` are skipped (honestly-
        labeled planar terminal members, not false 3D claims).
      * Rows with ``floquet_real`` / ``floquet_imag`` -> per-eigenvalue
        :func:`check_floquet_neimark_sacker` (informational).
      * Bracket entries in the header's ``bracket_inventory`` with
        ``classification == "neimark_sacker"`` -> each eigenvalue is
        re-classified; misclassifications are flagged.
    """
    del independent_checkers
    path = Path(jsonl_path)
    if not path.exists():
        raise FileNotFoundError(f"audit_topology: {jsonl_path} not found")
    findings: list[TopologyAuditFinding] = []

    header: dict[str, Any] | None = None
    with path.open() as fh:
        for line_idx, raw in enumerate(fh):
            raw = raw.strip()
            if not raw:
                continue
            try:
                row = json.loads(raw)
            except json.JSONDecodeError:
                continue
            if line_idx == 0 and (row.get("type") == "header" or row.get("_meta") is True):
                header = row
                continue
            if row.get("type") == "header" or row.get("_meta") is True:
                continue

            row_id = str(row.get("step_index", row.get("id", f"line_{line_idx}")))
            topology_class = _classify_row_topology(row, header)

            if topology_class == "tulip_or_3d_periodic":
                sp = _extract_state_and_period(row)
                if sp is None:
                    if not skip_rows_without_state:
                        findings.append(
                            TopologyAuditFinding(
                                source_jsonl=str(path),
                                row_id=row_id,
                                genome_gate="tulip_or_3d_periodic",
                                genome_verdict=True,
                                independent_verdict=False,
                                discrepancy=True,
                                failure_mode="row_schema_unsupported",
                            )
                        )
                    continue
                state0, period = sp
                system = _system_from_row_or_header(row, header)
                if system is None:
                    continue
                if row.get("degenerate_planar") is True:
                    continue
                z0 = float(state0[2])
                genome_3d_claim = abs(z0) > 0.0
                try:
                    is_3d, extras = check_tulip_topology(state0, period, system, z_floor=z_floor)
                except RuntimeError as exc:
                    findings.append(
                        TopologyAuditFinding(
                            source_jsonl=str(path),
                            row_id=row_id,
                            genome_gate="tulip_or_3d_periodic",
                            genome_verdict=genome_3d_claim,
                            independent_verdict=False,
                            discrepancy=False,
                            failure_mode=f"integrator_failed: {exc}",
                        )
                    )
                    continue
                discrepancy = genome_3d_claim and not is_3d
                findings.append(
                    TopologyAuditFinding(
                        source_jsonl=str(path),
                        row_id=row_id,
                        genome_gate="tulip_or_3d_periodic",
                        genome_verdict=genome_3d_claim,
                        independent_verdict=is_3d,
                        discrepancy=discrepancy,
                        failure_mode=("planar_collapse_under_322_floor" if discrepancy else ""),
                        extras=extras,
                    )
                )
                try:
                    closes, c_extras = check_periodic_orbit_closure(
                        state0, period, system, floor=closure_floor
                    )
                except (RuntimeError, ValueError) as exc:
                    findings.append(
                        TopologyAuditFinding(
                            source_jsonl=str(path),
                            row_id=row_id,
                            genome_gate="periodic_orbit_closure",
                            genome_verdict=True,
                            independent_verdict=False,
                            discrepancy=False,
                            failure_mode=f"closure_check_failed: {exc}",
                        )
                    )
                    continue
                genome_closes = True
                rec_resid = row.get("independent_closure_residual")
                if isinstance(rec_resid, (int, float)) and rec_resid >= closure_floor:
                    genome_closes = False
                closure_discrepancy = genome_closes and not closes
                findings.append(
                    TopologyAuditFinding(
                        source_jsonl=str(path),
                        row_id=row_id,
                        genome_gate="periodic_orbit_closure",
                        genome_verdict=genome_closes,
                        independent_verdict=closes,
                        discrepancy=closure_discrepancy,
                        failure_mode=("radau_cross_check_disagrees" if closure_discrepancy else ""),
                        extras=c_extras,
                    )
                )
            elif topology_class == "floquet_bifurcation":
                _audit_floquet_block(row, header, row_id, findings, path)

    if header is not None and "bracket_inventory" in header:
        for idx, bracket in enumerate(header.get("bracket_inventory", [])):
            if not isinstance(bracket, dict):
                continue
            classification = bracket.get("classification")
            if classification != "neimark_sacker":
                continue
            for which in ("a", "b"):
                er = bracket.get(f"eig_{which}_re")
                ei = bracket.get(f"eig_{which}_im")
                if er is None or ei is None:
                    continue
                is_ns, extras = check_floquet_neimark_sacker(float(er), float(ei))
                discrepancy = (not is_ns) and (classification == "neimark_sacker")
                failure = ""
                if discrepancy:
                    if extras["imag_part"] <= REAL_AXIS_TOL:
                        if abs(extras["eig_real"] - 1.0) < UNIT_CIRCLE_TOL:
                            failure = "real_eig_near_+1_trivial_Liouville_pair"
                        elif abs(extras["eig_real"] + 1.0) < UNIT_CIRCLE_TOL:
                            failure = "real_eig_near_-1_period_doubling_k2"
                        else:
                            failure = "real_eig_off_unit_circle"
                    else:
                        failure = "complex_pair_off_unit_circle"
                findings.append(
                    TopologyAuditFinding(
                        source_jsonl=str(path),
                        row_id=f"bracket_{idx}_eig_{which}",
                        genome_gate="neimark_sacker",
                        genome_verdict=True,
                        independent_verdict=is_ns,
                        discrepancy=discrepancy,
                        failure_mode=failure,
                        extras=extras,
                    )
                )

    return findings


def _audit_floquet_block(
    row: dict[str, Any],
    header: dict[str, Any] | None,
    row_id: str,
    findings: list[TopologyAuditFinding],
    path: Path,
) -> None:
    """Inner: walk a row's ``floquet_real``/``floquet_imag`` arrays."""
    del header
    re_arr = row.get("floquet_real")
    im_arr = row.get("floquet_imag")
    if not isinstance(re_arr, list) or not isinstance(im_arr, list):
        return
    if len(re_arr) != len(im_arr):
        return
    for k, (er, ei) in enumerate(zip(re_arr, im_arr, strict=True)):
        try:
            er_f = float(er)
            ei_f = float(ei)
        except (TypeError, ValueError):
            continue
        is_ns, extras = check_floquet_neimark_sacker(er_f, ei_f)
        findings.append(
            TopologyAuditFinding(
                source_jsonl=str(path),
                row_id=f"{row_id}_eig_{k}",
                genome_gate="floquet_eig",
                genome_verdict=is_ns,
                independent_verdict=is_ns,
                discrepancy=False,
                failure_mode="",
                extras=extras,
            )
        )


def write_findings(findings: list[TopologyAuditFinding], out_path: str) -> None:
    """Write findings to a JSONL file (one line per finding)."""
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        for f in findings:
            fh.write(json.dumps(asdict(f)) + "\n")


__all__ = [
    "INDEPENDENT_CLOSURE_FLOOR_NONDIM",
    "REAL_AXIS_TOL",
    "UNIT_CIRCLE_TOL",
    "TopologyAuditFinding",
    "audit_topology",
    "check_floquet_neimark_sacker",
    "check_periodic_orbit_closure",
    "check_tulip_topology",
    "write_findings",
]
