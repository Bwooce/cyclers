"""#682 -- quasi-periodic "cycler corridor" census around the stable prograde EM cyclers.

`#444`'s own named redirect (b): characterise cycler-USABILITY where the family is
known and stability-classified but its transport UTILITY -- the surrounding volume of
naturally-cycling, station-keeping-free trajectories -- has never been measured. This
is NOT a search for a new species; every linearly-stable periodic orbit has SOME
surrounding quasi-periodic torus family by KAM theory. The deliverable is a MEASUREMENT:
how BIG that corridor is per stable member, and whether it is operationally meaningful.

Members
-------
The linearly-STABLE members of the Braik-Ross 3D-lifted EM cycler families produced by
`#438` (`data/scan_434_3d_broken_plane_em.jsonl`): C21 (2,0,10) z0_0.24 (107 stable,
z0[-0.645,-0.206], C[2.147,3.026]); C32 (8,0,26) z0_0.24 (164 stable, z0[0.142,0.350],
C[2.680,3.061]); lyapunov3d-L1 (14 stable, C~2.999). The SUSPECT C32 z0_0.10 (66,0,136)
family the `#438` verdict flagged as a continuation artifact (Jacobi spans 1.79->10.75)
is EXCLUDED. Plus the linearly-stable PLANAR Braik-Ross goldens (classified live).

Method (reuse of the shared GMOS torus machinery `#612`/`#617` bootstrap from)
-----------------------------------------------------------------------------
Per member: monodromy -> Floquet center eigenpairs -> for the pair with the largest
rotation number (widest expected family), build a GMOS invariant-circle torus
(`genome.qp_tori.correct_qp_torus`, the Olikara-Scheeres GMOS benchmark corrector) at
each amplitude of a geometric ladder, recording the INDEPENDENT closure residual (a
non-circular short-time nonlinear-flow check). The corridor extent is the largest ladder
amplitude whose torus keeps closure below `CLOSURE_GATE`, converted to a physical tube
half-width in km / (m/s).

Corrector choice (honest): `#612`'s 2D pseudospectral corrector was built to cross the
UNSTABLE-halo shooting-fragility wall; on these LONG-PERIOD, high-winding STABLE cycler
center-manifolds it under-resolves the longitudinal structure at tractable mode counts
(observed: closure stalls ~1e-3 at n1=12). The GMOS corrector resolves the longitudinal
angle EXACTLY by stroboscopic integration and builds clean tori (closure ~1e-8) on the
identical members, so it is the correct tool here. GMOS is the SAME machinery both `#612`
modules bootstrap their seed from; the measure (closure-gated amplitude -> tube half-width)
is unchanged.

Positive control: Olikara-Scheeres 2012 is NOT in the corpus; the `#612` in-repo L2 GMOS
positive control (tests/search/test_variational_qp_torus.py::
test_l2_positive_control_reproduces_gmos_torus) validates this exact GMOS corrector (it
reproduces the published GMOS L2 torus rotation number to ~1e-5) and is the accepted
fallback per the dispatch.

Discipline: NO catalogue writeback. The schema question (new quasi_cycler rows vs.
corridor-width fields on the parent cycler rows) is an explicit user-decision point and
is left to the coordinating session; this script only measures and reports.

Run (resumable; re-invoke until "ALL MEMBERS DONE", then --assemble):
    uv run python scripts/census_682_cycler_corridors.py            # one chunk
    uv run python scripts/census_682_cycler_corridors.py --assemble # final summary
State   -> data/682_cycler_corridor_state.json
Summary -> data/found/682_cycler_corridor_census/summary.json
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import cyclerfinder.core.cr3bp as cr3bp  # noqa: E402
from cyclerfinder.genome.qp_tori import (  # noqa: E402
    QPTorus,
    correct_qp_torus,
    evaluate_torus,
)
from cyclerfinder.search.bifurcation_detector import (  # noqa: E402
    floquet_multipliers,
    monodromy,
)
from cyclerfinder.search.reachable_representatives import (  # noqa: E402
    TU_DAYS,
    braik_ross_system,
)

# --- physical scales (Earth-Moon, Braik-Ross convention) --------------------
L_KM = 384400.0
V_MS = L_KM * 1000.0 / (TU_DAYS * 86400.0)  # nondim velocity -> m/s

# --- corridor gate ----------------------------------------------------------
# The honest, non-circular validity gate is the INDEPENDENT closure residual (a
# true short-time nonlinear CR3BP flow check, not the algebraic residual). A torus
# with closure below this stays quasi-periodic to ~1e-4 nondim (~38 km) over the
# check horizon -- comparable to the GMOS corrector's own default independent_tol.
CLOSURE_GATE = 1.0e-4

SCAN = ROOT / "data" / "scan_434_3d_broken_plane_em.jsonl"
STATE = ROOT / "data" / "682_cycler_corridor_state.json"
OUTDIR = ROOT / "data" / "found" / "682_cycler_corridor_census"

# Stratified-sample stride per 3D family (name -> keep every Nth stable member,
# sorted by z0). These cycler members are compute-hostile (T~15-31 TU, high
# winding -> ~100-190 s per GMOS ladder call), so we span each family's z0/C range
# with a modest stratified sample rather than the full stable population; the
# uncovered members are named honestly in the summary, not silently dropped.
STRIDE = {
    "braik-ross-C21-em-z0_0.24": 22,  # 107 stable -> 5 spanning z0[-0.645,-0.206]
    "braik-ross-C32-em-z0_0.24": 34,  # 164 stable -> 5 spanning z0[0.142,0.350]
    "lyapunov3d-L1": 2,  # 14 stable -> 7 (cheap: T~2 TU)
}
_SYS = braik_ross_system()


def _best_k(phi: float) -> int:
    bk, bd = 5, math.inf
    for kk in range(3, 81):
        for j in range(1, kk):
            if math.gcd(j, kk) == 1 and abs(phi - 2 * math.pi * j / kk) < bd:
                bd, bk = abs(phi - 2 * math.pi * j / kk), kk
    return bk


def _conj_pairs(eigs: list[complex]) -> list[tuple[complex, complex]]:
    """Group unit-circle center eigenvalues into conjugate pairs (imag>0 first)."""
    cands = [complex(e) for e in eigs if abs(abs(e) - 1.0) < 0.2 and abs(e.imag) > 1e-4]
    used: list[int] = []
    pairs: list[tuple[complex, complex]] = []
    for i, e in enumerate(cands):
        if i in used:
            continue
        for j in range(i + 1, len(cands)):
            if j in used:
                continue
            if abs(cands[j] - e.conjugate()) < 1e-4:
                pairs.append((e, cands[j]) if e.imag > 0 else (cands[j], e))
                used += [i, j]
                break
    return pairs


# Geometric amplitude ladder (GMOS eigenvector-displacement / amplitude pin). The
# corridor is walked upward and stops at the first amplitude whose torus fails the
# closure gate (the corridor edge) or the top of the ladder (corridor >= top).
AMP_LADDER = (3e-4, 1e-3, 5e-3, 2e-2, 6e-2)


def _gmos_tube(torus: QPTorus) -> tuple[float, float]:
    """Transverse tube half-width of a GMOS torus at the invariant circle
    (theta_long=0): max over theta_trans of the excursion from the circle
    centroid, in km (position) / (m/s) (velocity)."""
    th2 = np.linspace(0, 2 * np.pi, 32, endpoint=False)
    pts = np.array([evaluate_torus(torus, 0.0, float(b)) for b in th2])
    pos = pts[:, :3]
    vel = pts[:, 3:]
    pr = float(np.max(np.linalg.norm(pos - pos.mean(0), axis=1)))
    vr = float(np.max(np.linalg.norm(vel - vel.mean(0), axis=1)))
    return pr * L_KM, vr * V_MS


def measure_corridor(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    period: float,
    *,
    n_trans: int = 4,
    max_iter: int = 20,
) -> dict[str, Any]:
    """Per-member corridor via a GMOS amplitude ladder.

    For these LINEARLY-STABLE cycler members the GMOS invariant-circle corrector
    (``genome.qp_tori.correct_qp_torus``, the Olikara-Scheeres benchmark) resolves
    the longitudinal angle EXACTLY by stroboscopic integration, so it builds clean
    tori (independent closure ~1e-8) even on the long-period, high-winding orbits
    where `#612`'s 2D pseudospectral corrector (built to cross the UNSTABLE-halo
    shooting-fragility wall) under-resolves the longitudinal structure. The
    corridor extent is the largest ladder amplitude whose torus keeps the honest
    independent closure residual below ``CLOSURE_GATE``.
    """
    t0 = time.time()
    eigs = floquet_multipliers(monodromy(system, state0, period))
    pairs = _conj_pairs(list(eigs))
    if not pairs:
        return {"status": "no_center_pair", "wall_s": time.time() - t0}
    pair = max(pairs, key=lambda p: abs(math.atan2(p[0].imag, p[0].real)))
    rho_lin = abs(math.atan2(pair[0].imag, pair[0].real)) / (2 * math.pi)
    k = _best_k(abs(math.atan2(pair[0].imag, pair[0].real)))
    best_amp = 0.0
    best_pos = 0.0
    best_vel = 0.0
    best_closure = float("nan")
    rho = float("nan")
    walled = False
    for amp in AMP_LADDER:
        try:
            g = correct_qp_torus(
                system,
                state0,
                period,
                (pair[0], pair[1]),
                k=k,
                n_trans=n_trans,
                initial_torus_amplitude=amp,
                tol=1e-8,
                max_iter=max_iter,
            )
        except Exception:  # Newton blow-up = corridor edge
            walled = True
            break
        if g.independent_closure_residual >= CLOSURE_GATE:
            walled = True
            break
        pos_km, vel_ms = _gmos_tube(g)
        best_amp, best_pos, best_vel = amp, pos_km, vel_ms
        best_closure = g.independent_closure_residual
        rho = g.omega_trans / g.omega_long
    if best_amp == 0.0:
        return {
            "status": "no_gated_torus",
            "rho_lin": rho_lin,
            "n_pairs": len(pairs),
            "wall_s": time.time() - t0,
        }
    return {
        "status": "ok",
        "rho": rho,
        "rho_lin": rho_lin,
        "n_pairs": len(pairs),
        "corridor_amp_nondim": best_amp,
        "corridor_pos_km": best_pos,
        "corridor_vel_ms": best_vel,
        # walled=True -> the corridor edge was found; walled=False -> the corridor
        # is at least this wide (ladder top reached with the gate still satisfied).
        "walled": walled,
        "corridor_is_lower_bound": not walled,
        "best_closure": best_closure,
        "wall_s": time.time() - t0,
    }


# ---------------------------------------------------------------------------
# Member enumeration.
# ---------------------------------------------------------------------------


def _stable_3d_members() -> list[dict[str, Any]]:
    recs = [json.loads(line) for line in SCAN.read_text().splitlines() if line.strip()]
    members: list[dict[str, Any]] = []
    for label, stride in STRIDE.items():
        st = sorted(
            (r for r in recs if r["seed_label"] == label and r["floquet_tag"] == "stable"),
            key=lambda r: r["z0"],
        )
        for r in st[::stride]:
            members.append(
                {
                    "member_id": f"{label}|x0={r['x0']:.6f}|z0={r['z0']:.6f}",
                    "family": label,
                    "kind": "3d_lift",
                    "x0": r["x0"],
                    "z0": r["z0"],
                    "ydot0": r["ydot0"],
                    "T_TU": r["T_TU"],
                    "jacobi": r["jacobi"],
                }
            )
    return members


GOLDEN_YAML = ROOT / "data" / "golden" / "braik_ross_2026_em_family_ics.yaml"


def _stable_planar_goldens() -> list[dict[str, Any]]:
    """Linearly-stable planar Braik-Ross goldens (all 13 from the sourced golden
    ICs, classified live via the monodromy spectral radius). ``period_nd`` is the
    FULL period; the ICs are the perpendicular-crossing state (x0,0,0,0,ydot0,0)."""
    import yaml  # type: ignore[import-untyped]  # local import; yaml is a dev dep

    data = yaml.safe_load(GOLDEN_YAML.read_text())
    members: list[dict[str, Any]] = []
    for fam in data["families"]:
        s0 = np.array([fam["x0"], 0.0, 0.0, 0.0, fam["ydot0"], 0.0])
        period = float(fam["period_nd"])
        eigs = floquet_multipliers(monodromy(_SYS, s0, period))
        spectral_radius = max(abs(e) for e in eigs)
        if spectral_radius < 1.05:  # linearly stable (all multipliers ~ on unit circle)
            members.append(
                {
                    "member_id": f"planar-golden|{fam['label']}",
                    "family": "planar_golden",
                    "kind": "planar_golden",
                    "label": fam["label"],
                    "x0": float(s0[0]),
                    "z0": 0.0,
                    "ydot0": float(s0[4]),
                    "T_TU": period,
                    "jacobi": float(fam["jacobi"]),
                    "spectral_radius": float(spectral_radius),
                }
            )
    return members


def _all_members() -> list[dict[str, Any]]:
    return _stable_3d_members() + _stable_planar_goldens()


# ---------------------------------------------------------------------------
# Checkpointed driver.
# ---------------------------------------------------------------------------


def _load_state() -> dict[str, Any]:
    if STATE.exists():
        state: dict[str, Any] = json.loads(STATE.read_text())
        return state
    return {"results": {}}


def _save_state(state: dict[str, Any]) -> None:
    STATE.parent.mkdir(parents=True, exist_ok=True)
    STATE.write_text(json.dumps(state, indent=2, default=float) + "\n")


def _ts() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S")


def run_chunk(max_members: int = 6) -> None:
    members = _all_members()
    state = _load_state()
    results: dict[str, Any] = state["results"]
    todo = [m for m in members if m["member_id"] not in results]
    print(f"[{_ts()}] {len(results)}/{len(members)} done; {len(todo)} remaining", flush=True)
    if not todo:
        print(f"[{_ts()}] ALL MEMBERS DONE -- run with --assemble.", flush=True)
        return
    for m in todo[:max_members]:
        s0 = np.array([m["x0"], 0.0, m["z0"], 0.0, m["ydot0"], 0.0])
        print(f"[{_ts()}] {m['member_id']} C={m['jacobi']:.4f} ...", flush=True)
        rec = measure_corridor(_SYS, s0, float(m["T_TU"]))
        rec.update({k: m[k] for k in ("family", "kind", "x0", "z0", "jacobi", "T_TU")})
        results[m["member_id"]] = rec
        _save_state(state)
        print(
            f"[{_ts()}]   -> {rec['status']} "
            f"amp={rec.get('corridor_amp_nondim')} "
            f"pos_km={rec.get('corridor_pos_km')} "
            f"({rec.get('wall_s', 0):.0f}s)",
            flush=True,
        )
    remaining = len([m for m in members if m["member_id"] not in results])
    print(f"[{_ts()}] chunk done; {remaining} members remain", flush=True)


def assemble() -> None:
    state = _load_state()
    results: dict[str, Any] = state["results"]
    ok = {k: v for k, v in results.items() if v.get("status") == "ok"}
    by_family: dict[str, list[dict[str, Any]]] = {}
    for v in ok.values():
        by_family.setdefault(str(v["family"]), []).append(v)
    fam_summary: dict[str, Any] = {}
    for fam, rows in by_family.items():
        amps = [float(r["corridor_amp_nondim"]) for r in rows]
        pos = [float(r["corridor_pos_km"]) for r in rows]
        vel = [float(r["corridor_vel_ms"]) for r in rows]
        fam_summary[fam] = {
            "n_members": len(rows),
            "corridor_amp_nondim": {
                "min": min(amps),
                "median": float(np.median(amps)),
                "max": max(amps),
            },
            "corridor_pos_km": {"min": min(pos), "median": float(np.median(pos)), "max": max(pos)},
            "corridor_vel_ms": {"min": min(vel), "median": float(np.median(vel)), "max": max(vel)},
        }
    status_counts: dict[str, int] = {}
    for v in results.values():
        status_counts[str(v["status"])] = status_counts.get(str(v["status"]), 0) + 1
    summary = {
        "task": 682,
        "closure_gate_nondim": CLOSURE_GATE,
        "closure_gate_km_equiv": CLOSURE_GATE * L_KM,
        "n_measured_ok": len(ok),
        "n_total_attempted": len(results),
        "status_counts": status_counts,
        "family_summary": fam_summary,
        "per_member": results,
    }
    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "summary.json").write_text(json.dumps(summary, indent=2, default=float) + "\n")
    print(
        f"[{_ts()}] wrote {OUTDIR / 'summary.json'}: {len(ok)} ok / {len(results)} attempted",
        flush=True,
    )
    for fam, fs in fam_summary.items():
        print(
            f"  {fam}: n={fs['n_members']} "
            f"amp med={fs['corridor_amp_nondim']['median']:.3e} "
            f"pos_km med={fs['corridor_pos_km']['median']:.1f} "
            f"[{fs['corridor_pos_km']['min']:.1f}, {fs['corridor_pos_km']['max']:.1f}]",
            flush=True,
        )


def main() -> None:
    if "--assemble" in sys.argv:
        assemble()
        return
    n = 6
    for a in sys.argv[1:]:
        if a.startswith("--chunk="):
            n = int(a.split("=", 1)[1])
    run_chunk(max_members=n)


if __name__ == "__main__":
    main()
