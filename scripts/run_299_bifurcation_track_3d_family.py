"""Bifurcation tracking + family switching on the 3D Earth-Moon (1,1) family.

Part B of #299 / #291 Phase 3.

Reads ``data/family_296_3d_em_11.jsonl`` (the Phase 2 tracer's 265 members) and
scans for period-multiplying bifurcations using the Floquet multipliers already
recorded with each member (no re-integration needed at the scan stage).

Each bracketed crossing is classified:

  * Floquet pair `e^{±iθ}` near a primitive k-th root of unity (θ = 2π j/k,
    gcd(j, k) = 1) with |λ| ≈ 1 → **period-multiplying / Neimark-Sacker**
    bifurcation; expect a kT sub-family branch.
  * Real multiplier crossing +1 → saddle-node (fold turning, no new branch).
  * Real multiplier crossing -1 → classic period-doubling (k=2).

For each viable period-multiplying bracket, the script attempts a family
switch:

  1. Take the bracketing parent member's monodromy (from the JSONL).
  2. Compute the right eigenvector of the multiplier nearest the primitive
     k-th root of unity.
  3. Project the eigenvector onto the SYMMETRIC-TULIP free-var components
     ``(z0, ydot0)`` (the IMAGINARY part of the complex eigenvector hits
     those components — the symmetry constrains y0 = xdot0 = zdot0 = 0).
  4. Perturb the parent IC along that direction with step ``+ε`` and ``-ε``.
  5. Re-correct via :func:`cr3bp_general_periodic_3d.correct_general_periodic_3d`
     with period guess ``k * T_parent`` and free vars
     ``FREE_VARS_SYMMETRIC_TULIP``.
  6. If the corrector lands on a periodic orbit with period ≈ kT, RECORD
     and continue the sub-family with the existing 3D tracer.

If the switched orbit's period collapses back to T (the corrector slid along
the parent family) or fails to converge, the bracket yields no new sub-family
— that is itself an honest signal (the bifurcation is geometrically valid in
Floquet structure but the kT branch is not accessible to the symmetric-tulip
corrector at the symmetric-IC submanifold; multi-shooting or asymmetric
free-vars would be the next escalation).

The output JSONL is:
  * Header: bracket inventory + classification.
  * Per accepted sub-family member: same shape as the Phase 2 family JSONL.

Usage::

    uv run python scripts/run_299_bifurcation_track_3d_family.py
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.bifurcation_detector import (
    floquet_multipliers,
    monodromy,
)
from cyclerfinder.search.cr3bp_3d_family_tracer import (
    Family3DMember,
    continue_general_3d_family,
)
from cyclerfinder.search.cr3bp_general_periodic_3d import (
    FREE_VARS_SYMMETRIC_TULIP,
    RESIDUAL_PERPENDICULAR_HALF_PERIOD,
    correct_general_periodic_3d,
)

FAMILY_PATH = Path("/home/bruce/dev/cyclers/data/family_296_3d_em_11.jsonl")
OUT_PATH = Path("/home/bruce/dev/cyclers/data/family_296_3d_subfamilies_299.jsonl")

EM_MU = 1.2150584270572e-2
EM_L_KM = 384400.0
EM_T_S = 375699.8


# ---------------------------------------------------------------------------
# Bracket detection (replicates scan_family_for_bifurcations on stored eigs).
# ---------------------------------------------------------------------------


def _primitive_roots(k: int) -> list[complex]:
    return [
        complex(math.cos(2 * math.pi * j / k), math.sin(2 * math.pi * j / k))
        for j in range(1, k)
        if math.gcd(j, k) == 1
    ]


def _nearest_kth_root_distance(eig: complex, k: int) -> tuple[float, complex]:
    roots = _primitive_roots(k)
    if not roots:
        return float("inf"), complex(1.0)
    dists = [abs(eig - r) for r in roots]
    idx = int(np.argmin(dists))
    return float(dists[idx]), roots[idx]


def _classify_bracket(
    eig_a: complex,
    eig_b: complex,
    k: int,
    *,
    unit_tol: float = 1e-2,
) -> str:
    """Tag the bracket by the multiplier's structural type.

    * ``period_doubling``: k=2 and the multiplier is real near -1.
    * ``saddle_node``: k>1 but the bracketing multiplier is real near +1 (i.e.
      a TRIVIAL unit multiplier shifted off-unity by numerical noise; not a
      genuine period-multiplying signal).
    * ``neimark_sacker``: complex pair on the unit circle near a primitive
      k-th root of unity for k >= 3 (or k=2 with real component near -1 AND
      imaginary near 0; classical NS for non-real k).
    * ``period_multiplying``: catch-all for k >= 3 real-root crossings (rare
      in symmetric tulip families; included for completeness).
    """
    on_unit_a = abs(abs(eig_a) - 1.0) < unit_tol
    on_unit_b = abs(abs(eig_b) - 1.0) < unit_tol
    on_unit = on_unit_a and on_unit_b
    has_imag = abs(eig_a.imag) > unit_tol or abs(eig_b.imag) > unit_tol
    if k == 2:
        if not has_imag and on_unit:
            return "period_doubling"
        return "period_doubling_complex"  # rare; complex pair crossing -1
    # k >= 3
    if on_unit and has_imag:
        return "neimark_sacker"
    if on_unit:
        return "period_multiplying"
    return "period_multiplying_offcircle"


def _scan_brackets(members: list[dict], *, k_max: int = 6, tol: float = 0.01) -> list[dict]:
    """Find adjacent-pair brackets where any multiplier crosses near a primitive k-th root."""
    per_member: list[dict[int, tuple[complex, float]]] = []
    for m in members:
        eigs = [
            complex(re, im) for re, im in zip(m["floquet_real"], m["floquet_imag"], strict=True)
        ]
        summary: dict[int, tuple[complex, float]] = {}
        for k in range(2, k_max + 1):
            best_d = float("inf")
            best_e = complex(1.0)
            for e in eigs:
                d, _ = _nearest_kth_root_distance(e, k)
                if d < best_d:
                    best_d, best_e = d, e
            summary[k] = (best_e, best_d)
        per_member.append(summary)

    brackets: list[dict] = []
    for i in range(len(members) - 1):
        m_a, m_b = members[i], members[i + 1]
        if m_b["step_index"] - m_a["step_index"] != 1:
            continue  # don't bracket across the seed-merge gap (e.g. -1 -> 0)
        sum_a, sum_b = per_member[i], per_member[i + 1]
        for k in range(2, k_max + 1):
            e_a, d_a = sum_a[k]
            e_b, d_b = sum_b[k]
            # Sign flip across tol band -> crossing.
            if (d_a - tol) * (d_b - tol) < 0.0:
                classification = _classify_bracket(e_a, e_b, k)
                brackets.append(
                    {
                        "k": k,
                        "classification": classification,
                        "step_a": int(m_a["step_index"]),
                        "step_b": int(m_b["step_index"]),
                        "T_a": float(m_a["T_TU"]),
                        "T_b": float(m_b["T_TU"]),
                        "z0_a": float(m_a["state_nd"][2]),
                        "z0_b": float(m_b["state_nd"][2]),
                        "C_a": float(m_a["jacobi_constant"]),
                        "C_b": float(m_b["jacobi_constant"]),
                        "d_a": d_a,
                        "d_b": d_b,
                        "eig_a_re": e_a.real,
                        "eig_a_im": e_a.imag,
                        "eig_b_re": e_b.real,
                        "eig_b_im": e_b.imag,
                    }
                )
    return brackets


# ---------------------------------------------------------------------------
# Family switch on a 3D symmetric tulip orbit.
# ---------------------------------------------------------------------------


def _select_period_multiplying_eigenvector(
    monodromy_matrix: np.ndarray,
    k: int,
) -> tuple[np.ndarray, complex, float] | None:
    """Pick the eigenvector of the monodromy whose eigenvalue is closest to a
    primitive k-th root of unity.

    Returns the IMAGINARY part of the complex eigenvector — for a symmetric
    tulip parent IC ``(x0, 0, z0, 0, ydot0, 0)`` the IMAGINARY part of the
    eigenvector of a primitive k-th-root-of-unity multiplier (k >= 3) hits the
    symmetric IC components ``(x, z, ydot)`` while the real part hits the
    perpendicular-zero components ``(y, xdot, zdot)``. Perturbing along the
    imaginary part keeps the perturbed IC inside the symmetric-tulip
    submanifold the corrector lives on.

    Returns ``(v, lam, dist)`` or ``None`` if no eigenvalue sits within 0.1
    of any primitive k-th root.
    """
    eigvals, eigvecs = np.linalg.eig(monodromy_matrix)
    roots = _primitive_roots(k)
    if not roots:
        return None
    best_dist = float("inf")
    best_idx = -1
    for i, lam in enumerate(eigvals):
        for r in roots:
            d = abs(complex(lam) - r)
            if d < best_dist:
                best_dist = d
                best_idx = i
    if best_idx < 0 or best_dist > 0.1:
        return None
    v_complex = eigvecs[:, best_idx]
    # For a real k=2 multiplier, the eigenvector is real; for complex pairs
    # (k >= 3), the IMAGINARY part hits the symmetric-tulip free vars.
    v = np.real(v_complex).astype(np.float64) if k == 2 else np.imag(v_complex).astype(np.float64)
    nrm = float(np.linalg.norm(v))
    if nrm < 1e-12:
        # Fall back to the other component.
        v = (np.real(v_complex) if k != 2 else np.imag(v_complex)).astype(np.float64)
        nrm = float(np.linalg.norm(v))
    if nrm < 1e-12:
        return None
    v = v / nrm
    return v, complex(eigvals[best_idx]), float(best_dist)


def _attempt_3d_family_switch(
    system: cr3bp.CR3BPSystem,
    parent: dict,
    k: int,
    *,
    eigenvector_step: float = 1e-3,
    corrector_tol: float = 1e-10,
    closure_tol: float = 1e-6,
    max_iter: int = 80,
) -> dict | None:
    """Try to land on the k-period sub-family branching off ``parent``.

    Returns the switched orbit IC + period + Jacobi + verification, or
    ``None`` if no convergence (both eigenvector signs tried; both periods
    tested for collapse-back-to-parent).
    """
    mono = np.array(parent["monodromy"], dtype=np.float64)
    pick = _select_period_multiplying_eigenvector(mono, k)
    if pick is None:
        return {"failure_reason": "no_kth_root_eigenvector"}
    v, lam, dist = pick
    state_parent = np.asarray(parent["state_nd"], dtype=np.float64)
    t_parent = float(parent["T_TU"])
    period_guess = float(k) * t_parent

    # Perturb the symmetric-tulip free vars: (x0_idx=0, z0_idx=2, ydot0_idx=4).
    # The corrector holds y0=xdot0=zdot0=0 (symmetric) and walks
    # (z0, ydot0, T); x0 is the natural parameter (fixed during correction).
    dz_v = float(v[2]) * eigenvector_step
    dyd_v = float(v[4]) * eigenvector_step

    attempts: list[dict] = []
    for sign in (+1.0, -1.0):
        state_pert = state_parent.copy()
        state_pert[2] = float(state_parent[2]) + sign * dz_v
        state_pert[4] = float(state_parent[4]) + sign * dyd_v
        # Hold y0=xdot0=zdot0=0 (symmetric tulip).
        state_pert[1] = 0.0
        state_pert[3] = 0.0
        state_pert[5] = 0.0
        try:
            result = correct_general_periodic_3d(
                system,
                state_pert,
                period_guess,
                free_vars=FREE_VARS_SYMMETRIC_TULIP,
                residual_indices=RESIDUAL_PERPENDICULAR_HALF_PERIOD,
                is_half_period_residual=True,
                tol=corrector_tol,
                max_iter=max_iter,
                independent_tol=closure_tol,
            )
        except RuntimeError:
            attempts.append({"sign": sign, "reason": "corrector_runtime"})
            continue
        if not result.converged:
            attempts.append(
                {
                    "sign": sign,
                    "reason": "not_converged",
                    "corrector_residual": float(result.corrector_residual),
                    "independent": float(result.independent_closure_residual),
                }
            )
            continue
        # Did we land back on the parent family (period ~ T_parent) or on the
        # genuine k*T branch?
        period_ratio = result.T_TU / t_parent
        accept = (period_ratio > 0.5 * k) and (period_ratio < 1.5 * k)
        # Sanity: the orbit must not have collapsed to planar.
        if result.degenerate_planar:
            attempts.append(
                {
                    "sign": sign,
                    "reason": "collapsed_to_planar",
                    "period_ratio": float(period_ratio),
                }
            )
            continue
        if not accept:
            attempts.append(
                {
                    "sign": sign,
                    "reason": "period_collapsed_to_parent",
                    "T_landed": float(result.T_TU),
                    "T_parent": float(t_parent),
                    "period_ratio": float(period_ratio),
                }
            )
            continue
        # Independent cross-check via fresh monodromy.
        try:
            mono_new = monodromy(system, result.state0, result.T_TU, rtol=1e-12, atol=1e-12)
            eigs_new = floquet_multipliers(mono_new)
        except RuntimeError:
            mono_new = None
            eigs_new = None
        return {
            "sign": sign,
            "k": k,
            "lam_parent": {"re": lam.real, "im": lam.imag},
            "lam_parent_dist": float(dist),
            "switched_state": result.state0.tolist(),
            "switched_T_TU": float(result.T_TU),
            "switched_T_days": float(result.T_TU * (system.t_s / 86400.0)),
            "switched_jacobi": float(result.jacobi),
            "switched_period_ratio": float(period_ratio),
            "switched_closure_corrector": float(result.corrector_residual),
            "switched_closure_independent": float(result.independent_closure_residual),
            "switched_floquet_abs": (
                [float(abs(e)) for e in eigs_new] if eigs_new is not None else None
            ),
        }
    return {"failure_reason": "no_sign_converged_to_kT", "attempts": attempts}


# ---------------------------------------------------------------------------
# Sub-family continuation.
# ---------------------------------------------------------------------------


def _member_to_dict(member: Family3DMember, *, system: cr3bp.CR3BPSystem) -> dict:
    """Mirror of the Phase 2 generator's member serializer."""
    orb = member.orbit
    rec = {
        "step_index": int(member.step_index),
        "arc_length": float(member.arc_length),
        "stability_tag": member.stability_tag,
        "state_nd": orb.state0.tolist(),
        "T_TU": float(orb.T_TU),
        "T_days": float(orb.T_TU * (system.t_s / 86400.0)),
        "jacobi_constant": float(orb.jacobi),
        "degenerate_planar": bool(orb.degenerate_planar),
        "corrector_residual": float(orb.corrector_residual),
        "independent_closure_residual": float(orb.independent_closure_residual),
        "n_iter": int(orb.n_iter),
    }
    if member.floquet is not None:
        rec["floquet_real"] = [float(z.real) for z in member.floquet]
        rec["floquet_imag"] = [float(z.imag) for z in member.floquet]
        rec["floquet_abs"] = [float(abs(z)) for z in member.floquet]
    return rec


def _continue_subfamily(
    system: cr3bp.CR3BPSystem,
    seed_state: np.ndarray,
    seed_period: float,
    *,
    n_steps_max: int = 30,
    step: float = 0.005,
) -> list[Family3DMember]:
    """Walk the sub-family pseudo-arclength from the switched seed (a short
    walk; the point is to confirm there IS a sub-family, not to map it).
    """
    fam = continue_general_3d_family(
        system,
        seed_state,
        seed_period,
        continuation="pseudo_arclength",
        step=step,
        n_steps_max=n_steps_max,
        direction="both",
        corrector_tol=1e-10,
        closure_tol=1e-6,
        fold_detection=True,
        monodromy_eval=True,
    )
    return fam.members


def main() -> int:
    if not FAMILY_PATH.exists():
        print(f"ERROR: {FAMILY_PATH} not found", file=sys.stderr)
        return 2
    rows: list[dict] = []
    with FAMILY_PATH.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rows.append(json.loads(line))
    members = sorted(rows[1:], key=lambda m: m["step_index"])
    print(f"Loaded {len(members)} family members from {FAMILY_PATH}", flush=True)

    print(f"[{time.strftime('%H:%M:%S')}] Scanning for bracketed period-multiplying crossings...")
    brackets = _scan_brackets(members, k_max=6, tol=0.01)
    print(f"  found {len(brackets)} brackets (k <= 6, tol=0.01)")

    # Classification histogram.
    cls_hist: dict[str, int] = {}
    for b in brackets:
        cls_hist[b["classification"]] = cls_hist.get(b["classification"], 0) + 1
    print(f"  classification histogram: {cls_hist}")

    # Pick a handful of best-quality brackets to attempt family switching on:
    # the closest crossing per k (smallest d_a + d_b), preferring k = 4, 3, 5, 6.
    brackets_by_k: dict[int, list[dict]] = {}
    for b in brackets:
        brackets_by_k.setdefault(b["k"], []).append(b)
    pick: list[dict] = []
    for k in (4, 3, 5, 6, 2):
        if k in brackets_by_k:
            pick.append(min(brackets_by_k[k], key=lambda b: b["d_a"] + b["d_b"]))
    print(f"\n  attempting family switch on {len(pick)} best brackets (one per k):")
    for b in pick:
        print(
            f"    k={b['k']}  step={b['step_a']:+4d}->{b['step_b']:+4d}  "
            f"T={b['T_a']:.4f}  z0={b['z0_a']:+.5f}  d={b['d_a']:.4f}  "
            f"cls={b['classification']}"
        )

    # System.
    system = cr3bp.CR3BPSystem(
        mu=EM_MU,
        primary="earth",
        secondary="moon",
        l_km=EM_L_KM,
        t_s=EM_T_S,
    )

    switch_results: list[dict] = []
    accepted_seeds: list[dict] = []  # successful switches
    for b in pick:
        parent = next(m for m in members if m["step_index"] == b["step_a"])
        print(
            f"\n[{time.strftime('%H:%M:%S')}] Trying family switch at "
            f"k={b['k']} step={b['step_a']:+d}..."
        )
        # Try several perturbation steps; if a tiny perturbation lands on the
        # parent (period_ratio ~ 1), a larger one is needed to escape the
        # parent's attractor basin.
        result: dict | None = None
        for eps in (1e-3, 5e-3, 1e-2, 3e-2):
            r = _attempt_3d_family_switch(
                system,
                parent,
                b["k"],
                eigenvector_step=eps,
                corrector_tol=1e-10,
                closure_tol=1e-6,
            )
            if r and "switched_state" in r:
                result = r
                result["eigenvector_step"] = eps
                break
            # No success at this eps; keep trying.
        if result is None or "switched_state" not in result:
            print(f"  no kT branch found for k={b['k']} at step {b['step_a']:+d}")
            switch_results.append({"bracket": b, "outcome": "no_branch_found", "trials": result})
            continue
        # Success: record + attempt a short sub-family walk.
        print(
            f"  ACCEPTED: T_landed={result['switched_T_TU']:.4f} TU "
            f"(ratio {result['switched_period_ratio']:.3f}x parent T)"
        )
        print(
            f"  closure: corrector={result['switched_closure_corrector']:.2e}  "
            f"independent={result['switched_closure_independent']:.2e}"
        )
        accepted_seeds.append({"bracket": b, "result": result})

        # Brief sub-family walk to confirm it IS a family, not an isolated orbit.
        seed_state = np.array(result["switched_state"], dtype=np.float64)
        try:
            sub_members = _continue_subfamily(
                system,
                seed_state,
                float(result["switched_T_TU"]),
                n_steps_max=20,
                step=0.005,
            )
        except (RuntimeError, ValueError) as exc:
            print(f"  sub-family continuation failed: {exc}")
            sub_members = []
        print(f"  sub-family walk: {len(sub_members)} members accepted")
        t_min = min((m.orbit.T_TU for m in sub_members), default=float("nan"))
        t_max = max((m.orbit.T_TU for m in sub_members), default=float("nan"))
        if sub_members:
            print(f"  sub-family T extent: {t_min:.4f} .. {t_max:.4f} TU")
        switch_results.append(
            {
                "bracket": b,
                "outcome": "accepted",
                "switched": result,
                "subfamily_n_members": len(sub_members),
                "subfamily_T_TU_extent": [t_min, t_max] if sub_members else None,
                "subfamily_members": [_member_to_dict(m, system=system) for m in sub_members],
            }
        )

    # Write outputs.
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    summary = {
        "type": "header",
        "issue": 299,
        "phase": "phase3_bifurcation_track_3d_family",
        "input": str(FAMILY_PATH),
        "n_parent_members": len(members),
        "n_brackets_total": len(brackets),
        "classification_histogram": cls_hist,
        "n_brackets_attempted": len(pick),
        "n_brackets_accepted": sum(1 for r in switch_results if r["outcome"] == "accepted"),
        "bracket_inventory": brackets,
        "discipline": (
            "Brackets are detected from the JSONL's stored Floquet multipliers "
            "(not re-integrated). Each accepted switched seed is cross-checked "
            "by (a) fresh full-period monodromy + Floquet, (b) the Phase 1 "
            "Radau independent-closure gate at 1e-6, and (c) the period ratio "
            "to the parent ~ k (collapse to parent T is REJECTED). Sub-family "
            "walks are short (n_steps_max=20) — confirmation that the seed "
            "lives on a family, not isolated."
        ),
    }
    with OUT_PATH.open("w") as f:
        f.write(json.dumps(summary) + "\n")
        for sr in switch_results:
            f.write(json.dumps(sr) + "\n")

    print()
    print(f"Wrote {OUT_PATH}  ({len(switch_results)} bracket-records + 1 header)")
    print(f"  accepted: {summary['n_brackets_accepted']}")
    print(f"  rejected: {len(switch_results) - summary['n_brackets_accepted']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
