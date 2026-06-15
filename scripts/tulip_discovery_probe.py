"""Tulip-orbit discovery probe at alternate three-body systems (#264 Phase B).

The #266 tulip genome ships ``find_tulip_via_continuation(np_target=2)`` — a
Phase-3 reproduce-gate routine that, at the Koblick Earth-Moon system, recovers
the published Np=2 butterfly via NRHO continuation + family-switching. This
driver re-runs that pipeline at three OTHER three-body systems where tulip
orbits have NOT been characterised in our sourced corpus:

* **Jupiter-Europa** (mu ~ 2.5e-5; small, gas-giant inner moon),
* **Saturn-Titan**  (mu ~ 2.4e-4; the Davis-Phillips-McCarthy 2018 existence
  prior — they found Np=6 tulips at Titan, so a probe here has a published
  existence anchor in a different family),
* **Pluto-Charon**  (mu ~ 0.108; the binary-system case, larger than
  Earth-Moon).

This is a DISCOVERY probe. Honest expectations:

* the Koblick seed IC (``x0=1.023731, z0=-0.18, ydot0=-0.103``) is an Earth-Moon
  NRHO IC — at a system with different mu it is generally NOT a periodic orbit,
  so the symmetric corrector at the new mu may fail to converge or land on a
  different family. A non-converged seed is a faithful negative (the genome did
  not find a tulip at this mu via the published Earth-Moon seed),
* even when the seed converges, the period-doubling bifurcation at k=2 may not
  be reachable inside ``n_steps_max`` continuation steps,
* even when a k=2 bifurcation is found, the family-switching corrector may not
  converge.

NO catalogue writeback. SILVER-class hits (closure + topology match) route to
``out/outcome_log/tulip_discovery_probe.jsonl`` for human review only. The
script never modifies ``data/catalogue.yaml`` or ``data/review_queue.jsonl``.

Run::

    uv run python scripts/tulip_discovery_probe.py
"""

from __future__ import annotations

import json
import subprocess
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np

from cyclerfinder.core import cr3bp
from cyclerfinder.genome.tulip import find_tulip_via_continuation, petal_count


def _git_sha() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


@dataclass
class ProbeOutcome:
    """One probe attempt at one system."""

    primary: str
    secondary: str
    mu: float
    l_km: float
    t_s: float
    seed_converged: bool
    n_branch_members: int
    bifurcation_found: bool
    bifurcation_x0: float | None
    bifurcation_distance: float | None
    switched_converged: bool
    petals: int | None
    period_tu: float | None
    x0: float | None
    ydot0: float | None
    z0: float | None
    success: bool
    reason: str
    wall_seconds: float
    git_sha: str

    def to_jsonl(self) -> str:
        return json.dumps(asdict(self), separators=(",", ":"), sort_keys=True)


def _probe_one(primary: str, secondary: str) -> ProbeOutcome:
    """Run ``find_tulip_via_continuation(np_target=2)`` at one system."""
    system = cr3bp.cr3bp_system(primary, secondary)
    sha = _git_sha()
    t0 = time.time()
    print(
        f"[tulip-probe] {primary}-{secondary}: mu={system.mu:.4e} "
        f"l_km={system.l_km:.4e} t_s={system.t_s:.4e}",
        flush=True,
    )
    try:
        result = find_tulip_via_continuation(np_target=2, system=system)
    except Exception as exc:  # never raise -- a faithful negative is the deliverable
        return ProbeOutcome(
            primary=primary,
            secondary=secondary,
            mu=system.mu,
            l_km=system.l_km,
            t_s=system.t_s,
            seed_converged=False,
            n_branch_members=0,
            bifurcation_found=False,
            bifurcation_x0=None,
            bifurcation_distance=None,
            switched_converged=False,
            petals=None,
            period_tu=None,
            x0=None,
            ydot0=None,
            z0=None,
            success=False,
            reason=f"exception:{type(exc).__name__}:{exc}",
            wall_seconds=time.time() - t0,
            git_sha=sha,
        )

    seed_ok = result.seed.converged
    bif = result.bifurcation
    switched = result.switched

    petals_n: int | None = None
    if switched is not None and switched.converged:
        try:
            state0 = np.array(
                [
                    switched.x0,
                    0.0,
                    switched.z0,
                    0.0,
                    switched.ydot0,
                    0.0,
                ],
                dtype=np.float64,
            )
            petals_n = int(petal_count(state0, switched.T_TU, system))
        except Exception:
            petals_n = None

    out = ProbeOutcome(
        primary=primary,
        secondary=secondary,
        mu=system.mu,
        l_km=system.l_km,
        t_s=system.t_s,
        seed_converged=seed_ok,
        n_branch_members=len(result.branch_members),
        bifurcation_found=bif is not None,
        bifurcation_x0=(float(bif.members[0].state0[0]) if bif is not None else None),
        bifurcation_distance=(
            float(min(bif.dist_before, bif.dist_after)) if bif is not None else None
        ),
        switched_converged=(switched is not None and switched.converged),
        petals=petals_n,
        period_tu=(float(switched.T_TU) if switched is not None else None),
        x0=(float(switched.x0) if switched is not None else None),
        ydot0=(float(switched.ydot0) if switched is not None else None),
        z0=(float(switched.z0) if switched is not None else None),
        success=bool(result.success and petals_n == 2),
        reason=result.reason,
        wall_seconds=time.time() - t0,
        git_sha=sha,
    )
    print(
        f"[tulip-probe] {primary}-{secondary}: "
        f"seed={'OK' if seed_ok else 'FAIL'} "
        f"branch={out.n_branch_members} bif={'YES' if out.bifurcation_found else 'NO'} "
        f"switched={'OK' if out.switched_converged else 'FAIL'} "
        f"petals={petals_n} success={out.success} reason={out.reason} "
        f"wall={out.wall_seconds:.1f}s",
        flush=True,
    )
    return out


def main() -> None:
    systems: list[tuple[str, str]] = [
        ("Jupiter", "Europa"),
        ("Saturn", "Titan"),
        ("Pluto", "Charon"),
    ]
    repo_root = Path(__file__).resolve().parent.parent
    out_path = repo_root / "out" / "outcome_log" / "tulip_discovery_probe.jsonl"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    outcomes: list[ProbeOutcome] = []
    for primary, secondary in systems:
        outcomes.append(_probe_one(primary, secondary))

    # Append every outcome -- positive or negative -- so the run is auditable.
    with out_path.open("a", encoding="utf-8") as fh:
        for o in outcomes:
            fh.write(o.to_jsonl() + "\n")

    summary: dict[str, Any] = {
        "n_systems": len(outcomes),
        "n_seed_ok": sum(1 for o in outcomes if o.seed_converged),
        "n_bifurcation": sum(1 for o in outcomes if o.bifurcation_found),
        "n_switched_ok": sum(1 for o in outcomes if o.switched_converged),
        "n_success_np2": sum(1 for o in outcomes if o.success),
    }
    print(f"[tulip-probe] SUMMARY {summary}", flush=True)


if __name__ == "__main__":
    main()
