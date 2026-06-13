"""Catalogue-wide FBS-analytic vs Lambert+FD parity sweep on the REAL corrector (#244).

DECISION-GATE evidence for #245 (the default flip). #243's fair trial showed the
FBS analytic gradients beat finite differences on an ISOLATED match-point NLP. The
two #243 caveats this run closes:

  1. Confirm the advantage holds on the REAL production chain corrector
     (:func:`cyclerfinder.search.dsm_leg.dsm_chain_correct` via
     :func:`cyclerfinder.search.dsm_descriptor_seed.close_row_dsm`), NOT a toy NLP.
  2. The patched-conic flyby-continuity constraints are wired with analytic
     gradients (``search/fbs_optimize_flyby.py``) — exercised by the FD-vs-analytic
     Jacobian cross-check below and its own test module.

The sweep
---------
Every catalogue row that has a DSM optimisation formulation (a g/G descriptor the
:func:`seed_dsm_chain_from_descriptor` lane can seed) is closed TWICE through the
identical corrector + identical charged seed — once with ``gradient="lambert"``
(the incumbent Lambert + SLSQP/least_squares finite-difference lane) and once with
``gradient="fbs-analytic"`` (each leg solved by the Ellison-2018 FBS match-point
corrector with the #226 analytic Jacobian). The ONLY difference is the
gradient/leg-evaluation backbone, so any divergence is attributable to it.

Same-model golden discipline
----------------------------
The comparison target per row is the Lambert+FD lane's OWN converged result on the
IDENTICAL dynamical model (DE440 via ``Ephemeris("astropy")``) and identical seed —
never a published/cross-model value. The analytic Jacobian is independently
cross-checked against central differences (the #243 ≤2.2e-7 discipline, repeated
here on the real flyby-wired chain). NO catalogue writeback.

Run: ``uv run python scripts/fbs_optimizer_adoption_parity.py``
"""

from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np

from cyclerfinder.core.constants import MU_SUN_KM3_S2
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.lambert import lambert
from cyclerfinder.data.catalog import load_catalog
from cyclerfinder.search.dsm_descriptor_seed import (
    close_row_dsm,
    seed_dsm_chain_from_descriptor,
)
from cyclerfinder.search.fbs_optimize import ChainLegSpec
from cyclerfinder.search.fbs_optimize_flyby import (
    FlybyChainSpec,
    _flyby_constraint_jac_fd,
    _pack_x0,
)

# Parity tolerances (km/s). The corrector's own tol_kms is 0.1; two lanes that
# reach the SAME basin agree on dV / emerged V_inf far tighter than that. A row is
# a DISAGREEMENT if the lanes differ by more than this, or one converges and the
# other does not.
_DV_PARITY_TOL_KMS = 0.5
_VINF_PARITY_TOL_KMS = 0.5
DAY = 86400.0


@dataclass
class RowParity:
    row_id: str
    sequence: tuple[str, ...]
    lam_converged: bool
    fbs_converged: bool
    lam_total_dv: float
    fbs_total_dv: float
    lam_max_res: float
    fbs_max_res: float
    lam_vinf: tuple[float, ...]
    fbs_vinf: tuple[float, ...]
    lam_wall_s: float
    fbs_wall_s: float
    dv_gap: float
    vinf_gap: float
    disagreement: bool
    disagreement_reason: str


def _vinf_gap(a: tuple[float, ...], b: tuple[float, ...]) -> float:
    """Max abs difference between aligned emerged-V_inf tuples (inf on length skew)."""
    if len(a) != len(b) or not a:
        return float("inf")
    return max(abs(x - y) for x, y in zip(a, b, strict=True))


def jacobian_crosscheck() -> float:
    """FD-vs-analytic max relative error of the flyby-wired constraint Jacobian.

    The #243 discipline (analytic gradient independently trusted) repeated on the
    real flyby-continuity-wired chain (caveat 2). Returns the max relative error.
    """
    eph = Ephemeris("astropy")
    t0 = 0.0
    seq = ["E", "M", "E"]
    tofs = [210 * DAY, 520 * DAY]
    epochs = [t0, t0 + tofs[0], t0 + tofs[0] + tofs[1]]
    states = [eph.state(b, e) for b, e in zip(seq, epochs, strict=True)]
    pos = [np.asarray(s[0]) for s in states]
    vpl = [np.asarray(s[1]) for s in states]
    legs = tuple(ChainLegSpec(r0=pos[i], rf=pos[i + 1], tof_s=tofs[i], alpha=0.5) for i in range(2))
    spec = FlybyChainSpec(legs=legs, bodies=tuple(seq), v_planets=tuple(vpl))
    leg_v = [lambert(pos[i], pos[i + 1], tofs[i], mu=MU_SUN_KM3_S2)[0] for i in range(2)]
    dvs = [np.zeros(3), np.zeros(3)]
    varr = [np.asarray(leg_v[0].v2), np.asarray(leg_v[1].v2)]
    vdep = [np.asarray(leg_v[0].v1), np.asarray(leg_v[1].v1)]
    x0 = _pack_x0(spec, dvs, varr, vdep)
    j_ana = _flyby_constraint_jac_fd(spec, x0, analytic=True)
    j_fd = _flyby_constraint_jac_fd(spec, x0, analytic=False)
    rel = np.abs(j_ana - j_fd) / np.maximum(np.abs(j_fd), 1.0)
    return float(np.max(rel))


def sweep_row(row_id: str, raw: dict, eph: Ephemeris) -> RowParity:
    """Close one row both ways and score the parity."""
    seed = seed_dsm_chain_from_descriptor(raw)
    assert seed is not None  # caller filters to seedable rows

    t = time.perf_counter()
    lam = close_row_dsm(raw, eph, gradient="lambert")
    lam_wall = time.perf_counter() - t

    t = time.perf_counter()
    fbs = close_row_dsm(raw, eph, gradient="fbs-analytic")
    fbs_wall = time.perf_counter() - t

    lam_dv = float(sum(lam.dv_dsm_kms))
    fbs_dv = float(sum(fbs.dv_dsm_kms))
    dv_gap = abs(lam_dv - fbs_dv)
    vinf_gap = _vinf_gap(lam.vinf_per_encounter_kms, fbs.vinf_per_encounter_kms)

    reasons: list[str] = []
    if lam.converged != fbs.converged:
        reasons.append(f"converged mismatch (lam={lam.converged}, fbs={fbs.converged})")
    if lam.converged and fbs.converged:
        if dv_gap > _DV_PARITY_TOL_KMS:
            reasons.append(f"dV gap {dv_gap:.3f} > {_DV_PARITY_TOL_KMS}")
        if vinf_gap > _VINF_PARITY_TOL_KMS:
            reasons.append(f"Vinf gap {vinf_gap:.3f} > {_VINF_PARITY_TOL_KMS}")
    elif not lam.converged and not fbs.converged:
        # NEITHER lane reached the tol — parity is UNDEFINED, not "clean". Flag it
        # so the aggregate verdict can never silently call this a pass. The useful
        # signal in this regime is the residual gap (which lane gets closer).
        reasons.append(
            f"neither converged (lam_res={lam.max_residual_kms:.3f}, "
            f"fbs_res={fbs.max_residual_kms:.3f})"
        )
    disagreement = bool(reasons)
    return RowParity(
        row_id=row_id,
        sequence=seed.sequence,
        lam_converged=lam.converged,
        fbs_converged=fbs.converged,
        lam_total_dv=lam_dv,
        fbs_total_dv=fbs_dv,
        lam_max_res=float(lam.max_residual_kms),
        fbs_max_res=float(fbs.max_residual_kms),
        lam_vinf=lam.vinf_per_encounter_kms,
        fbs_vinf=fbs.vinf_per_encounter_kms,
        lam_wall_s=lam_wall,
        fbs_wall_s=fbs_wall,
        dv_gap=dv_gap,
        vinf_gap=vinf_gap,
        disagreement=disagreement,
        disagreement_reason="; ".join(reasons) if reasons else "-",
    )


def main() -> None:
    print("FBS-analytic vs Lambert+FD parity sweep on the REAL DSM chain corrector (#244)")
    print("Same-model golden = the Lambert+FD lane's own converged result. NO writeback.\n")

    rel = jacobian_crosscheck()
    print(f"Jacobian FD-vs-analytic cross-check (flyby-wired chain): max rel err = {rel:.2e}")
    print(f"  -> analytic gradient {'TRUSTED' if rel < 1e-5 else 'SUSPECT'} (#243 discipline)\n")

    cat = load_catalog()
    rows = []
    for e in cat.entries:
        try:
            seed = seed_dsm_chain_from_descriptor(e.raw)
        except Exception:  # a malformed row must not abort the sweep
            seed = None
        if seed is not None:
            rows.append((e.id, e.raw))

    print(f"DSM-seedable rows: {len(rows)}\n")
    eph = Ephemeris("astropy")
    results: list[RowParity] = []
    for row_id, raw in rows:
        print(f"  sweeping {row_id} ...", flush=True)
        results.append(sweep_row(row_id, raw, eph))

    # --- parity table ---
    print("\n" + "=" * 100)
    print("PARITY TABLE (lam = Lambert+FD incumbent; fbs = FBS-analytic gradient)")
    print("=" * 100)
    hdr = (
        f"{'row':24} {'seq':10} {'lcv':4} {'fcv':4} "
        f"{'lam_res':>8} {'fbs_res':>8} {'res_x':>6} "
        f"{'lam_dV':>8} {'fbs_dV':>8} {'lam_s':>6} {'fbs_s':>6}"
    )
    print(hdr)
    print("-" * 100)
    for r in results:
        res_x = r.lam_max_res / r.fbs_max_res if r.fbs_max_res > 0 else float("inf")
        print(
            f"{r.row_id:24} {'-'.join(r.sequence):10} "
            f"{r.lam_converged!s:4.4} {r.fbs_converged!s:4.4} "
            f"{r.lam_max_res:8.3f} {r.fbs_max_res:8.3f} {res_x:6.2f} "
            f"{r.lam_total_dv:8.2f} {r.fbs_total_dv:8.2f} {r.lam_wall_s:6.1f} {r.fbs_wall_s:6.1f}"
        )

    # --- disagreement ledger ---
    print("\n" + "=" * 100)
    print("DISAGREEMENT LEDGER (rows where the two lanes diverge — for investigation)")
    print("=" * 100)
    diffs = [r for r in results if r.disagreement]
    if not diffs:
        print("  (none) — every row's FBS-analytic lane matched the Lambert+FD lane within tol.")
    else:
        for r in diffs:
            print(f"  {r.row_id}: {r.disagreement_reason}")
            print(f"      lam: conv={r.lam_converged} dV={r.lam_total_dv:.3f} vinf={r.lam_vinf}")
            print(f"      fbs: conv={r.fbs_converged} dV={r.fbs_total_dv:.3f} vinf={r.fbs_vinf}")

    # --- aggregate verdict ---
    n = len(results)
    n_diff = len(diffs)
    n_lam_cv = sum(r.lam_converged for r in results)
    n_fbs_cv = sum(r.fbs_converged for r in results)
    both_conv = [r for r in results if r.lam_converged and r.fbs_converged]
    neither_conv = [r for r in results if not r.lam_converged and not r.fbs_converged]
    print("\n" + "=" * 100)
    print("AGGREGATE")
    print("=" * 100)
    print(f"  rows swept:                 {n}")
    print(f"  lam converged:              {n_lam_cv}")
    print(f"  fbs converged:              {n_fbs_cv}")
    print(f"  both converged:             {len(both_conv)}")
    print(f"  neither converged:          {len(neither_conv)}")
    print(f"  disagreements:              {n_diff}")
    if both_conv:
        max_dv = max(r.dv_gap for r in both_conv)
        max_vinf = max(r.vinf_gap for r in both_conv)
        print(f"  max dV gap (both-conv):     {max_dv:.4f} km/s")
        print(f"  max Vinf gap (both-conv):   {max_vinf:.4f} km/s")
        lam_tot = sum(r.lam_wall_s for r in both_conv)
        fbs_tot = sum(r.fbs_wall_s for r in both_conv)
        ratio = fbs_tot / lam_tot
        print(f"  wall (both-conv): lam={lam_tot:.1f}s  fbs={fbs_tot:.1f}s  ratio={ratio:.2f}x")
    if neither_conv:
        # In the zero-convergence regime, parity over the convergence flag is
        # vacuous. The decision-relevant signal is whether the FBS-analytic lane
        # gets STRICTLY CLOSER to feasibility (lower residual) on the same seed.
        fbs_closer = sum(1 for r in neither_conv if r.fbs_max_res < r.lam_max_res)
        res_ratios = [r.lam_max_res / r.fbs_max_res for r in neither_conv if r.fbs_max_res > 0]
        med_ratio = sorted(res_ratios)[len(res_ratios) // 2] if res_ratios else float("nan")
        print(
            f"  neither-conv residual signal: fbs strictly closer on "
            f"{fbs_closer}/{len(neither_conv)} rows; median lam/fbs residual = {med_ratio:.2f}x"
        )

    # The verdict must NOT call zero-convergence "clean". Parity is only
    # established on rows where BOTH lanes actually converge.
    if len(both_conv) == 0:
        verdict = (
            "NO CONVERGENCE ON EITHER LANE — parity undefined; "
            "real-corrector blocker is seed/basin, NOT gradient backbone"
        )
    elif n_diff == 0:
        verdict = "CLEAN PARITY (all both-converged rows agree within tol)"
    else:
        verdict = f"{n_diff} DISAGREEMENT(S)"
    print(f"\n  VERDICT: {verdict}")


if __name__ == "__main__":
    main()
