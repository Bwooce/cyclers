"""Higher-Np tulip continuation off impactor branch (#283).

Continues #281's multi-mu tulip Np=2 characterization by attempting Np>=3 at:
  * Earth-Moon (Np in {3,4,5,6}), using the Koblick AMOSTECH Table 4 paper IC
    rows that #266 Phase 4 partially mapped. The paper rows for Np>=3 carry
    NEGATIVE ``r_min_km`` (the published IC sample sits on the *impactor*
    branch — the orbit passes through the lunar interior at the sample point).
    The negative ``r_min_km`` is GEOMETRICALLY MEANINGFUL: it indexes a valid
    bifurcation root of the family curve. The orbit IS periodic at the IC; it
    is only unphysical because it would graze inside the lunar radius. The
    contribution here is recording the Np>=3 family-member parameters as
    closed periodic orbits in the CR3BP — the family is Koblick-published; the
    impactor-vs-physical-branch distinction is downstream.
  * Other 4 systems where #281 landed Np=2 (Jupiter-Ganymede, Saturn-Titan,
    Neptune-Triton, Pluto-Charon) at Np=3, using both direct-seed and a
    multi-shoot fallback from each system's Np=2 converged state.

The deliverable is ``data/tulip_higher_np_283.jsonl``. **NO catalogue
writeback.** Tulips are not cyclers; the JSONL is the record.

Convergence + cross-check discipline:
  * Multi-shoot ``max_segment_residual < 1e-7``.
  * Forward propagate the converged IC by T_TU; closure residual < 1e-6 (the
    same threshold used in #281's cross-check column).
  * Topological gate: ``petal_count(state0, T_TU, system) == np_target``. If
    the corrector lands a periodic orbit at the WRONG petal count, record it
    as ``petal_mismatch`` — that family member is not the target.
  * Earth-Moon ONLY: cross-check ``T_TU`` against
    :data:`KOBLICK_2023_TABLE4_PAPER[Np]['tau0']`. ``|dT/T|`` reported as
    ``deviation_vs_paper_pct``; >2% flags for review.
"""

from __future__ import annotations

import json
import signal
import subprocess
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np

from cyclerfinder.core import cr3bp as cr3bp_mod
from cyclerfinder.genome.multi_shooting import multi_shoot_periodic
from cyclerfinder.genome.tulip import (
    KOBLICK_2023_TABLE4_PAPER,
    find_tulip_at_system,
    koblick_system,
    petal_count,
)

CLOSURE_TOL: float = 1e-6
MULTISHOOT_TOL: float = 1e-7
MAX_ITER: int = 120
# Per-tier wall-clock cap. Task brief: "If a particular Np's corrector hangs
# >5 min on a single system, record timeout-negative and move on."
TIER_TIMEOUT_S: int = 300


class _TimeoutError(Exception):
    pass


def _timeout_handler(_signum: int, _frame: Any) -> None:
    raise _TimeoutError("corrector timeout")


def _git_sha() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


@dataclass
class Row:
    """One (system, Np) attempt."""

    system: str
    mu: float
    l_km: float
    t_s: float
    np_target: int
    converged: bool
    path: str
    x0: float | None
    ydot0: float | None
    z0: float | None
    T_nondim: float | None
    T_seconds: float | None
    jacobi_c: float | None
    n_petals_observed: int | None
    closure_residual: float | None
    max_segment_residual: float | None
    deviation_vs_paper_pct: float | None
    notes: str
    elapsed_seconds: float
    git_sha: str

    def to_jsonl(self) -> str:
        return json.dumps(asdict(self), separators=(",", ":"), sort_keys=True)


def _attempt_direct_paper_seed(
    system: cr3bp_mod.CR3BPSystem,
    np_target: int,
) -> tuple[bool, dict | None, str]:
    """Attempt the Koblick paper IC (Np row) directly at the target mu.

    Returns ``(converged, result_dict, notes)``. ``result_dict`` has keys
    ``state0, period, max_segment_residual``.
    """
    if np_target not in KOBLICK_2023_TABLE4_PAPER:
        return False, None, f"no_paper_row_for_np={np_target}"
    row = KOBLICK_2023_TABLE4_PAPER[np_target]
    parent_state = np.array(
        [
            float(row["x0"]),
            0.0,
            float(row["z0"]),
            0.0,
            float(row["ydot0"]),
            0.0,
        ],
        dtype=np.float64,
    )
    parent_period = float(row["tau0"])
    try:
        result = multi_shoot_periodic(
            system,
            parent_state,
            parent_period,
            n_segments=np_target,
            tol=MULTISHOOT_TOL,
            max_iter=MAX_ITER,
        )
    except (RuntimeError, ValueError) as exc:
        return False, None, f"exception:{type(exc).__name__}:{exc}"
    if not result.converged:
        return False, None, f"direct_seed_no_converge:resid={result.max_segment_residual:.3e}"
    return (
        True,
        {
            "state0": result.state0.copy(),
            "period": float(result.period),
            "max_segment_residual": float(result.max_segment_residual),
        },
        "direct_seed",
    )


def _attempt_multishoot_from_np2(
    system: cr3bp_mod.CR3BPSystem,
    np_target: int,
) -> tuple[bool, dict | None, str]:
    """Compute Np=2 at this system, then attempt multi-shoot Np=target from it.

    The parent_period is scaled by ``np_target/2`` to seed the corrector near
    the right multi-petal period.
    """
    try:
        res2 = find_tulip_at_system(system, np_target=2, multi_shooting=True, tol=1e-9)
    except Exception as exc:
        return False, None, f"np2_exception:{type(exc).__name__}:{exc}"
    if res2 is None or not res2.success or res2.switched is None:
        reason = res2.reason if res2 is not None else "result_none"
        return False, None, f"np2_no_converge:{reason}"
    sw = res2.switched
    parent_state = np.array([sw.x0, 0.0, sw.z0, 0.0, sw.ydot0, 0.0], dtype=np.float64)
    parent_period = float(sw.T_TU) * (np_target / 2.0)
    try:
        result = multi_shoot_periodic(
            system,
            parent_state,
            parent_period,
            n_segments=np_target,
            tol=MULTISHOOT_TOL,
            max_iter=MAX_ITER,
        )
    except (RuntimeError, ValueError) as exc:
        return False, None, f"multishoot_exception:{type(exc).__name__}:{exc}"
    if not result.converged:
        return False, None, f"multishoot_no_converge:resid={result.max_segment_residual:.3e}"
    return (
        True,
        {
            "state0": result.state0.copy(),
            "period": float(result.period),
            "max_segment_residual": float(result.max_segment_residual),
        },
        "multi_shoot_from_np2",
    )


def _closure_residual(
    system: cr3bp_mod.CR3BPSystem,
    state0: np.ndarray,
    period: float,
) -> float:
    """Forward-propagate by T and return |state(T) - state(0)|.

    Independent cross-check: the multi-shooter enforces segment patches; this
    re-integrates the whole period in one DOP853 pass.
    """
    arc = cr3bp_mod.propagate(system, state0, period, rtol=1e-12, atol=1e-12)
    return float(np.linalg.norm(arc.state_f - state0))


def _safe_petal_count(
    state0: np.ndarray,
    period: float,
    system: cr3bp_mod.CR3BPSystem,
) -> int:
    """Return petal count or -1 if the classifier fails (e.g. deep low-perilune)."""
    try:
        return int(petal_count(state0, period, system))
    except RuntimeError:
        return -1


def _measure(
    system: cr3bp_mod.CR3BPSystem,
    sys_label: str,
    np_target: int,
    path: str,
    state0: np.ndarray,
    period: float,
    max_segment_residual: float,
    notes_extra: str,
    elapsed: float,
    sha: str,
) -> Row:
    """Pack a successful corrector outcome into a Row, including cross-checks."""
    try:
        closure_res = _closure_residual(system, state0, period)
    except Exception as exc:
        closure_res = float("nan")
        notes_extra += f";closure_check_exception:{exc}"
    n_petals = _safe_petal_count(state0, period, system)
    jacobi = float(cr3bp_mod.jacobi_constant(state0, system.mu))
    # paper-deviation cross-check (Earth-Moon only)
    dev_pct: float | None = None
    if (
        system.mu == koblick_system().mu  # same mu => same nondim normalisation
        and np_target in KOBLICK_2023_TABLE4_PAPER
    ):
        paper_t = float(KOBLICK_2023_TABLE4_PAPER[np_target]["tau0"])
        dev_pct = 100.0 * (period - paper_t) / paper_t

    # gates: closure < 1e-6 AND petal_count == np_target
    converged = bool(closure_res < CLOSURE_TOL and n_petals == np_target)
    notes_parts: list[str] = []
    if not (closure_res < CLOSURE_TOL):
        notes_parts.append(f"closure_above_tol:{closure_res:.3e}")
    if n_petals != np_target:
        notes_parts.append(f"petal_mismatch:got={n_petals}_expected={np_target}")
    if notes_extra:
        notes_parts.append(notes_extra)
    notes_str = ";".join(notes_parts) if notes_parts else "ok"

    return Row(
        system=sys_label,
        mu=float(system.mu),
        l_km=float(system.l_km),
        t_s=float(system.t_s),
        np_target=int(np_target),
        converged=converged,
        path=path,
        x0=float(state0[0]),
        ydot0=float(state0[4]),
        z0=float(state0[2]),
        T_nondim=float(period),
        T_seconds=float(period * system.t_s),
        jacobi_c=jacobi,
        n_petals_observed=int(n_petals) if n_petals >= 0 else None,
        closure_residual=float(closure_res),
        max_segment_residual=float(max_segment_residual),
        deviation_vs_paper_pct=float(dev_pct) if dev_pct is not None else None,
        notes=notes_str,
        elapsed_seconds=elapsed,
        git_sha=sha,
    )


def _negative_row(
    system: cr3bp_mod.CR3BPSystem,
    sys_label: str,
    np_target: int,
    path: str,
    notes: str,
    elapsed: float,
    sha: str,
) -> Row:
    return Row(
        system=sys_label,
        mu=float(system.mu),
        l_km=float(system.l_km),
        t_s=float(system.t_s),
        np_target=int(np_target),
        converged=False,
        path=path,
        x0=None,
        ydot0=None,
        z0=None,
        T_nondim=None,
        T_seconds=None,
        jacobi_c=None,
        n_petals_observed=None,
        closure_residual=None,
        max_segment_residual=None,
        deviation_vs_paper_pct=None,
        notes=notes,
        elapsed_seconds=elapsed,
        git_sha=sha,
    )


def _attempt_pipeline(
    system: cr3bp_mod.CR3BPSystem,
    sys_label: str,
    np_target: int,
    sha: str,
    *,
    try_direct: bool = True,
    try_multishoot: bool = True,
) -> Row:
    """Walk the escalation tiers; return the first match (or final negative)."""
    t0 = time.time()
    notes_combined: list[str] = []

    # Tier A: direct Koblick paper-row IC at this mu.
    if try_direct:
        signal.signal(signal.SIGALRM, _timeout_handler)
        signal.alarm(TIER_TIMEOUT_S)
        try:
            ok, payload, notes = _attempt_direct_paper_seed(system, np_target)
        except _TimeoutError:
            ok, payload, notes = False, None, f"direct_seed_timeout:>{TIER_TIMEOUT_S}s"
        finally:
            signal.alarm(0)
        if ok and payload is not None:
            row = _measure(
                system,
                sys_label,
                np_target,
                path="direct_seed",
                state0=payload["state0"],
                period=payload["period"],
                max_segment_residual=payload["max_segment_residual"],
                notes_extra="",
                elapsed=time.time() - t0,
                sha=sha,
            )
            if row.converged:
                return row
            # Direct-seed converged in corrector but gates fail; remember for
            # later but continue trying multi-shoot.
            notes_combined.append(f"direct_seed_gated:{row.notes}")
        else:
            notes_combined.append(notes)

    # Tier B: multi-shoot from this system's Np=2 result, with parent_period
    # scaled by Np/2.
    if try_multishoot:
        signal.signal(signal.SIGALRM, _timeout_handler)
        signal.alarm(TIER_TIMEOUT_S)
        try:
            ok, payload, notes = _attempt_multishoot_from_np2(system, np_target)
        except _TimeoutError:
            ok, payload, notes = (
                False,
                None,
                f"multishoot_timeout:>{TIER_TIMEOUT_S}s",
            )
        finally:
            signal.alarm(0)
        if ok and payload is not None:
            row = _measure(
                system,
                sys_label,
                np_target,
                path="multi_shoot_from_np2",
                state0=payload["state0"],
                period=payload["period"],
                max_segment_residual=payload["max_segment_residual"],
                notes_extra="",
                elapsed=time.time() - t0,
                sha=sha,
            )
            if row.converged:
                return row
            notes_combined.append(f"multishoot_gated:{row.notes}")
        else:
            notes_combined.append(notes)

    # Clean negative.
    return _negative_row(
        system=system,
        sys_label=sys_label,
        np_target=np_target,
        path="no_path_converged",
        notes="|".join(notes_combined) if notes_combined else "no_path_attempted",
        elapsed=time.time() - t0,
        sha=sha,
    )


def main() -> None:
    repo_root = Path(__file__).resolve().parent.parent
    out_path = repo_root / "data" / "tulip_higher_np_283.jsonl"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    sha = _git_sha()

    rows: list[Row] = []
    t_start = time.time()

    # PASS 1: Earth-Moon Np in {3, 4, 5, 6} via direct paper seed.
    em_system = koblick_system()
    for np_t in [3, 4, 5, 6]:
        print(f"[tulip-283] earth-moon np={np_t} ...", flush=True)
        row = _attempt_pipeline(em_system, "earth-moon", np_t, sha, try_multishoot=False)
        rows.append(row)
        print(
            f"[tulip-283] earth-moon np={np_t} converged={row.converged} "
            f"path={row.path} T_s={row.T_seconds} dev_pct={row.deviation_vs_paper_pct} "
            f"notes={row.notes}",
            flush=True,
        )

    # PASS 2: other systems at Np=3 via direct paper seed first, then
    # multi-shoot fallback.
    other_systems: list[tuple[str, str, str]] = [
        ("Jupiter", "Ganymede", "jupiter-ganymede"),
        ("Saturn", "Titan", "saturn-titan"),
        ("Neptune", "Triton", "neptune-triton"),
        ("Pluto", "Charon", "pluto-charon"),
    ]
    for primary, secondary, label in other_systems:
        sysm = cr3bp_mod.cr3bp_system(primary, secondary)
        print(f"[tulip-283] {label} np=3 mu={sysm.mu:.4e} ...", flush=True)
        row = _attempt_pipeline(sysm, label, 3, sha)
        rows.append(row)
        print(
            f"[tulip-283] {label} np=3 converged={row.converged} "
            f"path={row.path} T_s={row.T_seconds} notes={row.notes}",
            flush=True,
        )

    # Write JSONL.
    with out_path.open("w", encoding="utf-8") as fh:
        for r in rows:
            fh.write(r.to_jsonl() + "\n")

    summary = {
        "n_rows": len(rows),
        "n_converged": sum(1 for r in rows if r.converged),
        "wall_seconds": time.time() - t_start,
        "out_path": str(out_path),
    }
    print(f"[tulip-283] SUMMARY {summary}", flush=True)


if __name__ == "__main__":
    main()
