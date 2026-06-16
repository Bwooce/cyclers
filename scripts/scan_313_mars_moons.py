"""Mars-Phobos / Mars-Deimos tulip + multi-rev scouts (#313 Phase 1 Part A).

Cheap parallel scout running existing tulip-orbit and multi-shoot machinery at
fresh-mu CR3BP systems Mars-Phobos (mu~1.65e-8) and Mars-Deimos (mu~2.25e-9).
#281 single-shot already attempted Mars-Phobos at Np=2 and returned
``seed_no_converge`` (the Koblick Earth-Moon NRHO seed does not land any
periodic orbit at the Mars-Phobos basin). This script re-attempts both moons
with two additional levers — multi-shooting at n_segments in {2, 3, 4} and the
Koblick Table 4 PAPER row seeds at higher Np (where the seed differs from the
default Np=1 row used by find_tulip_at_system).

Discipline (#313 brief):
  * READ-ONLY on the modules src/cyclerfinder/genome/*.py — wraps, never
    modifies.
  * NO catalogue writeback. JSONL deliverable only.
  * NO novelty claims. Tulip family at Mars-Phobos is a mu-scaling of the
    published Koblick Earth-Moon family — a lit-fresh hit at Mars-mu would be
    a candidate for a future BCR4BP/QP-tori gauntlet, not a discovery.
  * Independent CR3BP forward-propagation closure cross-check on every
    converged member (same anchor as #281's multimu sweep).

Output: ``data/scan_313_mars_phobos.jsonl`` and
``data/scan_313_mars_deimos.jsonl``. One row per (system, attempt, Np_target).

Mars-moon mu values from the JPL SSD satellite registry already in
``src/cyclerfinder/core/satellites.py`` (Phobos GM 7.087e-4 km^3/s^2,
Deimos GM 9.62e-5; both sourced from JPL SSD phys_par MAR097 + mean elements,
accessed 2026-06-14). The mu values are recomputed at script start so any
registry update flows through.
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
from cyclerfinder.genome.multi_shooting import multi_shoot_periodic
from cyclerfinder.genome.tulip import (
    KOBLICK_2023_TABLE4_PAPER,
    find_tulip_at_system,
    petal_count,
)

# Per-attempt wall-clock cap. Mars-Phobos at mu~1e-8 is structurally hostile;
# we cap generously but still bounded so the script always terminates.
TIMEOUT_S: int = 180


def _git_sha() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


class _TimeoutError(Exception):
    pass


def _timeout_handler(_signum: int, _frame: Any) -> None:
    raise _TimeoutError("attempt timeout")


@dataclass
class ScanRow:
    system: str
    primary: str
    secondary: str
    mu: float
    l_km: float
    t_s: float
    attempt: str  # "tulip_direct" | "tulip_multishoot_n{n}" | "tulip_paperseed_np{n}"
    np_target: int
    converged: bool
    x0: float | None
    y0: float | None
    z0: float | None
    xdot0: float | None
    ydot0: float | None
    zdot0: float | None
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


def _closure_residual(
    system: cr3bp.CR3BPSystem,
    state0: np.ndarray,
    t_nondim: float,
) -> float:
    """Forward-propagate state0 over t_nondim via DOP853 and report L2 closure.

    Independent cross-check: the corrector enforces perpendicular-crossing
    closure via the half-period STM Newton, this re-integrates over the FULL
    period in raw DOP853 — agreement is the orbit-closure gate.
    """
    arc = cr3bp.propagate(system, state0, t_nondim, rtol=1e-12, atol=1e-12)
    return float(np.linalg.norm(arc.state_f - state0))


def _row_template(
    sys_label: str,
    primary: str,
    secondary: str,
    system: cr3bp.CR3BPSystem,
    attempt: str,
    np_target: int,
    elapsed: float,
    git_sha: str,
) -> ScanRow:
    return ScanRow(
        system=sys_label,
        primary=primary,
        secondary=secondary,
        mu=system.mu,
        l_km=system.l_km,
        t_s=system.t_s,
        attempt=attempt,
        np_target=np_target,
        converged=False,
        x0=None,
        y0=None,
        z0=None,
        xdot0=None,
        ydot0=None,
        zdot0=None,
        T_nondim=None,
        T_seconds=None,
        jacobi_c=None,
        n_petals_observed=None,
        closure_residual_propagated=None,
        reason="",
        elapsed_seconds=elapsed,
        git_sha=git_sha,
    )


def _success_row(
    template: ScanRow,
    system: cr3bp.CR3BPSystem,
    state0: np.ndarray,
    t_nondim: float,
    reason: str,
) -> ScanRow:
    template.converged = True
    template.x0 = float(state0[0])
    template.y0 = float(state0[1])
    template.z0 = float(state0[2])
    template.xdot0 = float(state0[3])
    template.ydot0 = float(state0[4])
    template.zdot0 = float(state0[5])
    template.T_nondim = float(t_nondim)
    template.T_seconds = float(t_nondim * system.t_s)
    template.jacobi_c = float(cr3bp.jacobi_constant(state0, system.mu))
    try:
        template.n_petals_observed = int(petal_count(state0, t_nondim, system))
    except Exception:
        template.n_petals_observed = -1
    try:
        template.closure_residual_propagated = _closure_residual(system, state0, t_nondim)
    except Exception as exc:
        template.closure_residual_propagated = float("nan")
        reason = f"{reason}/closure_exception:{type(exc).__name__}"
    template.reason = reason
    return template


def _attempt_tulip_direct(primary: str, secondary: str, np_target: int, git_sha: str) -> ScanRow:
    system = cr3bp.cr3bp_system(primary, secondary)
    sys_label = f"{primary.lower()}-{secondary.lower()}"
    t0 = time.time()
    template = _row_template(
        sys_label,
        primary,
        secondary,
        system,
        attempt="tulip_direct",
        np_target=np_target,
        elapsed=0.0,
        git_sha=git_sha,
    )
    signal.signal(signal.SIGALRM, _timeout_handler)
    signal.alarm(TIMEOUT_S)
    try:
        result = find_tulip_at_system(system, np_target=np_target, multi_shooting=True)
    except _TimeoutError:
        signal.alarm(0)
        template.reason = "attempt_timeout"
        template.elapsed_seconds = time.time() - t0
        return template
    except Exception as exc:
        signal.alarm(0)
        template.reason = f"exception:{type(exc).__name__}:{exc}"
        template.elapsed_seconds = time.time() - t0
        return template
    finally:
        signal.alarm(0)
    template.elapsed_seconds = time.time() - t0
    if result is None:
        template.reason = "seed_no_converge"
        return template
    switched = result.switched
    if (switched is None) or (not getattr(switched, "converged", False)) or (not result.success):
        template.reason = result.reason or "no_switched_member"
        return template
    state0 = np.array(
        [switched.x0, 0.0, switched.z0, 0.0, switched.ydot0, 0.0],
        dtype=np.float64,
    )
    return _success_row(template, system, state0, float(switched.T_TU), result.reason or "ok")


def _attempt_tulip_paperseed(primary: str, secondary: str, np_target: int, git_sha: str) -> ScanRow:
    """Try find_tulip_at_system with the Koblick Table 4 PAPER row at np_target.

    The default find_tulip_at_system uses KOBLICK_2023_TABLE4_PAPER[1] (the Np=1
    row) as the seed and walks. Here we override with the Np=np_target paper
    row directly — it may already sit closer to the target Np orbit at the new
    mu. Only legal for np_target in KOBLICK_2023_TABLE4_PAPER keys.
    """
    if np_target not in KOBLICK_2023_TABLE4_PAPER:
        raise ValueError(f"no paper row at Np={np_target}")
    system = cr3bp.cr3bp_system(primary, secondary)
    sys_label = f"{primary.lower()}-{secondary.lower()}"
    t0 = time.time()
    template = _row_template(
        sys_label,
        primary,
        secondary,
        system,
        attempt=f"tulip_paperseed_np{np_target}",
        np_target=np_target,
        elapsed=0.0,
        git_sha=git_sha,
    )
    signal.signal(signal.SIGALRM, _timeout_handler)
    signal.alarm(TIMEOUT_S)
    try:
        result = find_tulip_at_system(
            system,
            np_target=2,  # Tier B path uses paper-seed; np_target=2 is the only Tier B value
            multi_shooting=True,
            seed_row=KOBLICK_2023_TABLE4_PAPER[np_target],
        )
    except _TimeoutError:
        signal.alarm(0)
        template.reason = "attempt_timeout"
        template.elapsed_seconds = time.time() - t0
        return template
    except Exception as exc:
        signal.alarm(0)
        template.reason = f"exception:{type(exc).__name__}:{exc}"
        template.elapsed_seconds = time.time() - t0
        return template
    finally:
        signal.alarm(0)
    template.elapsed_seconds = time.time() - t0
    if result is None:
        template.reason = "seed_no_converge"
        return template
    switched = result.switched
    if (switched is None) or (not getattr(switched, "converged", False)) or (not result.success):
        template.reason = result.reason or "no_switched_member"
        return template
    state0 = np.array(
        [switched.x0, 0.0, switched.z0, 0.0, switched.ydot0, 0.0],
        dtype=np.float64,
    )
    return _success_row(template, system, state0, float(switched.T_TU), result.reason or "ok")


def _attempt_multishoot(
    primary: str,
    secondary: str,
    np_target: int,
    n_segments: int,
    git_sha: str,
) -> ScanRow:
    """Multi-shoot Koblick paper-row seed at np_target with n_segments arcs.

    This is the direct-corrector path: the Koblick Earth-Moon paper seed is
    fed straight into the multi-shooter at the target mu without any
    intermediate continuation. Designed to catch the case where the basin of
    attraction shifts at the new mu but is still reachable by multi-shooting
    rather than family-switch continuation.
    """
    if np_target not in KOBLICK_2023_TABLE4_PAPER:
        raise ValueError(f"no paper row at Np={np_target}")
    system = cr3bp.cr3bp_system(primary, secondary)
    sys_label = f"{primary.lower()}-{secondary.lower()}"
    t0 = time.time()
    template = _row_template(
        sys_label,
        primary,
        secondary,
        system,
        attempt=f"tulip_multishoot_n{n_segments}",
        np_target=np_target,
        elapsed=0.0,
        git_sha=git_sha,
    )
    row = KOBLICK_2023_TABLE4_PAPER[np_target]
    state_seed = np.array(
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
    period_seed = float(row["tau0"])
    signal.signal(signal.SIGALRM, _timeout_handler)
    signal.alarm(TIMEOUT_S)
    try:
        result = multi_shoot_periodic(
            system,
            state_seed,
            period_seed,
            n_segments=n_segments,
            tol=1e-9,
            max_iter=80,
        )
    except _TimeoutError:
        signal.alarm(0)
        template.reason = "attempt_timeout"
        template.elapsed_seconds = time.time() - t0
        return template
    except Exception as exc:
        signal.alarm(0)
        template.reason = f"exception:{type(exc).__name__}:{exc}"
        template.elapsed_seconds = time.time() - t0
        return template
    finally:
        signal.alarm(0)
    template.elapsed_seconds = time.time() - t0
    if not result.converged:
        template.reason = (
            f"multishoot_no_converge:max_resid={result.max_segment_residual:.3e}:"
            f"n_iter={result.n_iter}"
        )
        return template
    state0 = np.asarray(result.state0, dtype=np.float64)
    return _success_row(
        template,
        system,
        state0,
        float(result.period),
        f"multishoot_n{n_segments}_ok",
    )


def _scan_one_system(primary: str, secondary: str) -> list[ScanRow]:
    git_sha = _git_sha()
    sys_label = f"{primary.lower()}-{secondary.lower()}"
    print(
        f"[scan-313] {sys_label} sha={git_sha} start at {time.strftime('%Y-%m-%dT%H:%M:%S')}",
        flush=True,
    )
    rows: list[ScanRow] = []

    # Attempt 1: tulip_direct at Np=2 (replicates #281 mars-phobos seed_no_converge).
    print(f"[scan-313] {sys_label} attempt: tulip_direct Np=2", flush=True)
    r1 = _attempt_tulip_direct(primary, secondary, np_target=2, git_sha=git_sha)
    rows.append(r1)
    print(
        f"[scan-313] {sys_label} tulip_direct Np=2 converged={r1.converged} "
        f"reason={r1.reason} elapsed={r1.elapsed_seconds:.1f}s",
        flush=True,
    )

    # Attempt 2..4: multi-shoot from Koblick Np=2 paper seed at n_segments=2,3,4.
    for n_segs in (2, 3, 4):
        print(
            f"[scan-313] {sys_label} attempt: multishoot n={n_segs} Np=2",
            flush=True,
        )
        r = _attempt_multishoot(primary, secondary, np_target=2, n_segments=n_segs, git_sha=git_sha)
        rows.append(r)
        print(
            f"[scan-313] {sys_label} multishoot n={n_segs} converged={r.converged} "
            f"reason={r.reason} elapsed={r.elapsed_seconds:.1f}s",
            flush=True,
        )

    # Attempt 5,6: tulip_paperseed at Np=3 and Np=4 (higher-petal paper rows).
    for np_target in (3, 4):
        print(
            f"[scan-313] {sys_label} attempt: tulip_paperseed Np={np_target}",
            flush=True,
        )
        r = _attempt_tulip_paperseed(primary, secondary, np_target=np_target, git_sha=git_sha)
        rows.append(r)
        print(
            f"[scan-313] {sys_label} tulip_paperseed Np={np_target} "
            f"converged={r.converged} reason={r.reason} "
            f"elapsed={r.elapsed_seconds:.1f}s",
            flush=True,
        )

    return rows


def main() -> int:
    repo_root = Path(__file__).resolve().parent.parent
    data_dir = repo_root / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    targets: list[tuple[str, str, Path]] = [
        ("Mars", "Phobos", data_dir / "scan_313_mars_phobos.jsonl"),
        ("Mars", "Deimos", data_dir / "scan_313_mars_deimos.jsonl"),
    ]

    t_start = time.time()
    for primary, secondary, out_path in targets:
        rows = _scan_one_system(primary, secondary)
        with out_path.open("w", encoding="utf-8") as fh:
            for r in rows:
                fh.write(r.to_jsonl() + "\n")
        n_conv = sum(1 for r in rows if r.converged)
        print(
            f"[scan-313] WROTE {out_path} : {len(rows)} rows, {n_conv} converged",
            flush=True,
        )

    total_elapsed = time.time() - t_start
    print(
        f"[scan-313] PART A DONE : total wall = {total_elapsed:.1f}s",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
