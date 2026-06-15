"""Multi-mu tulip Np=2 characterization sweep (#281).

Drives :func:`cyclerfinder.genome.tulip.find_tulip_at_system` across the registry
of CR3BP systems available from :data:`cyclerfinder.core.satellites.SATELLITES`,
records per-system corrector outcomes, and cross-checks at least one converged
member by direct CR3BP forward-propagation closure.

This is a discovery characterization, NOT a cycler. Tulip orbits are
periodic-orbit-family members in the CR3BP, not Earth-Moon resonant cyclers,
so the deliverable is a results JSONL at ``data/tulip_sweep_281.jsonl``. There
is NO catalogue writeback and the script never claims novelty -- the tulip
family was published by Koblick 2023; this run characterises mu-specific
instances at moon systems other than the Earth-Moon and Saturn-Titan anchors
already established at #266 / #280.

Run::

    uv run python scripts/tulip_multimu_sweep.py

The JSONL output is the deliverable. Each line carries enough provenance
(system name, mu, np_target, converged flag, x0/ydot0/z0/T_TU/T_s/jacobi_C,
closure residual, failure reason, elapsed seconds, git SHA) to reproduce the
outcome.
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

from cyclerfinder.core import cr3bp
from cyclerfinder.genome.tulip import find_tulip_at_system, petal_count

# Per-system corrector timeout in seconds. Memory rule: "if a particular
# system's corrector hangs (>3 min), record as 'did-not-converge: corrector
# timeout' and move on." Implemented via SIGALRM (Linux/POSIX only).
TIMEOUT_S: int = 180


# Candidate three-body systems. All are buildable from the existing
# (primary, secondary) registry entries in
# :mod:`cyclerfinder.core.satellites` -- no new sourcing required.
#
# Selection rationale (per the #281 task brief):
#   * Earth-Moon       : baseline anchor (Koblick 2023, confirmed at #266)
#   * Saturn-Titan     : second anchor (Phase 5 #280)
#   * Jupiter-Galileans: Io/Europa/Ganymede/Callisto, full set
#   * Saturn extras    : Enceladus / Rhea / Iapetus (skip Mimas/Dione/Tethys
#                        which are dominated by other moons in this primary)
#   * Uranus           : Titania / Oberon (the large regulars)
#   * Neptune-Triton   : retrograde + inclined, but a large body -> a
#                        legitimate mu-test even though dynamically hostile
#   * Pluto-Charon     : mu ~ 0.108, the near-binary case
#   * Mars-Phobos      : tiny mu (~ 1.7e-8) -- the small-mu robustness test
SYSTEMS: list[tuple[str, str]] = [
    ("Earth", "Moon"),
    ("Mars", "Phobos"),
    ("Jupiter", "Io"),
    ("Jupiter", "Europa"),
    ("Jupiter", "Ganymede"),
    ("Jupiter", "Callisto"),
    ("Saturn", "Enceladus"),
    ("Saturn", "Rhea"),
    ("Saturn", "Titan"),
    ("Saturn", "Iapetus"),
    ("Uranus", "Titania"),
    ("Uranus", "Oberon"),
    ("Neptune", "Triton"),
    ("Pluto", "Charon"),
]


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
class SweepRow:
    """One Np-target attempt at one system."""

    system: str
    primary: str
    secondary: str
    mu: float
    l_km: float
    t_s: float
    np_target: int
    converged: bool
    x0: float | None
    ydot0: float | None
    z0: float | None
    T_nondim: float | None
    T_seconds: float | None
    jacobi_c: float | None
    n_petals_observed: int | None
    closure_residual_propagated: float | None
    reason: str
    elapsed_seconds: float
    git_sha: str

    def to_jsonl(self) -> str:
        return json.dumps(asdict(self), separators=(",", ":"), sort_keys=True)


def _propagate_closure_residual(
    system: cr3bp.CR3BPSystem,
    x0: float,
    ydot0: float,
    z0: float,
    t_nondim: float,
) -> float:
    """Forward-propagate the (x0,0,z0,0,ydot0,0) IC by t_nondim and report |state(T) - state(0)|.

    This is the INDEPENDENT cross-check required by the orbit-closure discipline:
    the corrector inside ``find_tulip_at_system`` enforces perpendicular-crossing
    closure (planar x-z-plane crossing at t=T/2 maps back to the IC at t=T) via
    half-period STM Newton; this routine re-integrates the same IC over the FULL
    period via raw DOP853 and measures the residual. Agreement is the gate.
    """
    state0 = np.array([x0, 0.0, z0, 0.0, ydot0, 0.0], dtype=np.float64)
    arc = cr3bp.propagate(system, state0, t_nondim, rtol=1e-12, atol=1e-12)
    return float(np.linalg.norm(arc.state_f - state0))


def _attempt_one(primary: str, secondary: str, np_target: int) -> SweepRow:
    system = cr3bp.cr3bp_system(primary, secondary)
    sha = _git_sha()
    sys_label = f"{primary.lower()}-{secondary.lower()}"
    t0 = time.time()
    print(
        f"[tulip-sweep] {sys_label} np={np_target} mu={system.mu:.4e} l_km={system.l_km:.4e}",
        flush=True,
    )

    # Set a per-attempt wall-clock alarm. The #281 brief calls out >3 min as the
    # cutoff; we use 180s.
    signal.signal(signal.SIGALRM, _timeout_handler)
    signal.alarm(TIMEOUT_S)

    try:
        result = find_tulip_at_system(
            system,
            np_target=np_target,
            multi_shooting=True,
        )
    except _TimeoutError:
        signal.alarm(0)
        return SweepRow(
            system=sys_label,
            primary=primary,
            secondary=secondary,
            mu=system.mu,
            l_km=system.l_km,
            t_s=system.t_s,
            np_target=np_target,
            converged=False,
            x0=None,
            ydot0=None,
            z0=None,
            T_nondim=None,
            T_seconds=None,
            jacobi_c=None,
            n_petals_observed=None,
            closure_residual_propagated=None,
            reason="corrector_timeout",
            elapsed_seconds=time.time() - t0,
            git_sha=sha,
        )
    except Exception as exc:
        signal.alarm(0)
        return SweepRow(
            system=sys_label,
            primary=primary,
            secondary=secondary,
            mu=system.mu,
            l_km=system.l_km,
            t_s=system.t_s,
            np_target=np_target,
            converged=False,
            x0=None,
            ydot0=None,
            z0=None,
            T_nondim=None,
            T_seconds=None,
            jacobi_c=None,
            n_petals_observed=None,
            closure_residual_propagated=None,
            reason=f"exception:{type(exc).__name__}:{exc}",
            elapsed_seconds=time.time() - t0,
            git_sha=sha,
        )
    finally:
        signal.alarm(0)

    elapsed = time.time() - t0

    if result is None:
        # find_tulip_at_system returns None only when the seed itself fails to
        # converge at the new mu (no periodic orbit anywhere near the Koblick
        # seed at this mu). That IS a faithful negative.
        return SweepRow(
            system=sys_label,
            primary=primary,
            secondary=secondary,
            mu=system.mu,
            l_km=system.l_km,
            t_s=system.t_s,
            np_target=np_target,
            converged=False,
            x0=None,
            ydot0=None,
            z0=None,
            T_nondim=None,
            T_seconds=None,
            jacobi_c=None,
            n_petals_observed=None,
            closure_residual_propagated=None,
            reason="seed_no_converge",
            elapsed_seconds=elapsed,
            git_sha=sha,
        )

    # FindTulipResult: we accept Tier A direct seed matches AND Tier B
    # family-switched matches. In both cases ``result.switched`` is the
    # converged orbit and result.success carries the petal-count gate.
    switched = result.switched
    success = bool(result.success)
    if (switched is None) or (not getattr(switched, "converged", False)) or (not success):
        return SweepRow(
            system=sys_label,
            primary=primary,
            secondary=secondary,
            mu=system.mu,
            l_km=system.l_km,
            t_s=system.t_s,
            np_target=np_target,
            converged=False,
            x0=None,
            ydot0=None,
            z0=None,
            T_nondim=None,
            T_seconds=None,
            jacobi_c=None,
            n_petals_observed=None,
            closure_residual_propagated=None,
            reason=result.reason or "no_switched_member",
            elapsed_seconds=elapsed,
            git_sha=sha,
        )

    # Converged: record the IC + period, the Jacobi constant, the independent
    # petal count, and a forward-propagation closure residual cross-check.
    x0 = float(switched.x0)
    ydot0 = float(switched.ydot0)
    z0 = float(switched.z0)
    t_nd = float(switched.T_TU)
    state0 = np.array([x0, 0.0, z0, 0.0, ydot0, 0.0], dtype=np.float64)
    jacobi_c = float(cr3bp.jacobi_constant(state0, system.mu))
    try:
        n_petals = int(petal_count(state0, t_nd, system))
    except Exception:
        n_petals = -1
    try:
        closure_res = _propagate_closure_residual(system, x0, ydot0, z0, t_nd)
    except Exception as exc:
        closure_res = float("nan")
        print(f"[tulip-sweep] closure-check failed at {sys_label}: {exc}", flush=True)

    return SweepRow(
        system=sys_label,
        primary=primary,
        secondary=secondary,
        mu=system.mu,
        l_km=system.l_km,
        t_s=system.t_s,
        np_target=np_target,
        converged=True,
        x0=x0,
        ydot0=ydot0,
        z0=z0,
        T_nondim=t_nd,
        T_seconds=t_nd * system.t_s,
        jacobi_c=jacobi_c,
        n_petals_observed=n_petals,
        closure_residual_propagated=closure_res,
        reason=result.reason or "ok",
        elapsed_seconds=elapsed,
        git_sha=sha,
    )


def main() -> None:
    repo_root = Path(__file__).resolve().parent.parent
    out_path = repo_root / "data" / "tulip_sweep_281.jsonl"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[SweepRow] = []
    t_start = time.time()
    for primary, secondary in SYSTEMS:
        row = _attempt_one(primary, secondary, np_target=2)
        rows.append(row)
        print(
            f"[tulip-sweep] {row.system} np={row.np_target} converged={row.converged} "
            f"T_s={row.T_seconds} reason={row.reason} elapsed={row.elapsed_seconds:.1f}s",
            flush=True,
        )

    # Overwrite the JSONL (not append) so the file holds one canonical sweep.
    with out_path.open("w", encoding="utf-8") as fh:
        for r in rows:
            fh.write(r.to_jsonl() + "\n")

    summary: dict[str, Any] = {
        "n_systems": len(rows),
        "n_converged": sum(1 for r in rows if r.converged),
        "wall_seconds": time.time() - t_start,
        "out_path": str(out_path),
    }
    print(f"[tulip-sweep] SUMMARY {summary}", flush=True)


if __name__ == "__main__":
    main()
