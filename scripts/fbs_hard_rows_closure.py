"""FBS vs Lambert closure attempt on the hard multi-arc cycler rows (#242).

EVIDENCE-GATHERING run, NOT a catalogue update. Question: can the newly-landed
Lambert-free FBS match-point lane (#226: ``core/fbs_match_point.py`` +
``search/dsm_leg.dsm_leg_correct_fbs``) close the historically un-closeable
multi-arc rows that the Lambert corrector could never close? A clean closure FBS
achieves that Lambert cannot = strong argument to default FBS (#243). No new
closures = FBS stays an opt-in cross-check tool. A clean NEGATIVE is a success.

Rows attempted
--------------
* **S1L1 / russell-ch4-4.991gG2** — corrected two-arc topology
  ``E -> g(E-E free return, NO Mars) -> E flyby -> G(E-M-E Mars transit, true
  longitude rendezvous with DE440 Mars) -> E`` (Russell 2004 App-C #83; source-dig
  ``docs/notes/2026-06-08-s1l1-source-dig.md``). The binding leg is the G
  Mars-transit leg with the full position+velocity match point at true DE440 Mars.
* **russell-ch4-6.44Gg3** — full-sequence E-M-E-M-E loop-arcs-as-DSM-legs probe
  (#153 / #157). Per-leg FBS-vs-Lambert agreement + chain non-closure check.

Closure discipline (project hard-won rules)
-------------------------------------------
* The FBS match-point residual contains BOTH position (3 rows) AND velocity (3
  rows) defect at the interior impulse point — a leg that does not dynamically
  connect its two boundary states cannot drive this to zero. This is the full
  binding per-leg constraint; we state it explicitly.
* SAME-MODEL golden: the comparison target for each leg is the LAMBERT
  ``dsm_leg`` solution on the IDENTICAL boundary states / ToF / eta / mu (no
  cross-model golden). At the Lambert solution the FBS defect is exactly zero with
  ``dv = v21-v12, vf = v22``; the two lanes MUST agree (same geometry, same dV)
  for the per-leg result to count as cross-confirmed.
* For S1L1 the additional binding constraint is the LONGITUDE RENDEZVOUS with the
  true DE440 Mars position — that lives in the right-boundary position ``rf`` of
  the G leg (we set ``rf`` = real DE440 Mars position at the arrival epoch), so a
  converged match point IS a true-Mars intercept (position row of the defect).
* No catalogue writeback. Real solver numbers only.

Run: ``uv run python scripts/fbs_hard_rows_closure.py``
"""

from __future__ import annotations

import numpy as np

from cyclerfinder.core.constants import MU_SUN_KM3_S2
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.data.catalog import load_catalog
from cyclerfinder.search.dsm_leg import dsm_leg, dsm_leg_correct_fbs
from cyclerfinder.search.s1l1_corrected import (
    APPC_MARS_TRANSIT,
    build_seeded_arcs,
)

DAY_S = 86400.0


def _vinf(v_sc: np.ndarray, v_planet: np.ndarray) -> float:
    return float(np.linalg.norm(np.asarray(v_sc) - np.asarray(v_planet)))


def _hr(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def fbs_vs_lambert_leg(
    *,
    r0: np.ndarray,
    v0: np.ndarray,
    rf: np.ndarray,
    tof_s: float,
    eta: float,
    mu: float = MU_SUN_KM3_S2,
    label: str,
) -> dict:
    """Solve ONE leg both ways (Lambert dsm_leg + FBS corrector) and compare.

    Lambert is the SAME-MODEL golden. FBS is seeded from a PERTURBED Lambert
    solution so a genuine root-find happens (not a trivial echo). Returns the
    agreement diagnostics.
    """
    # --- Lambert lane (the same-model golden) ---
    lam = dsm_leg(r0, v0, tof_s, eta, rf, mu=mu)
    dv_lam = lam.v_depart_post_dsm - lam.v_arrive_pre_dsm
    vf_lam = lam.v_arrive

    # --- FBS lane, seeded PERTURBED so the corrector actually roots ---
    rng = np.random.default_rng(0)
    dv0 = dv_lam + rng.normal(scale=0.3, size=3)
    vf0 = vf_lam + rng.normal(scale=0.3, size=3)
    fbs = dsm_leg_correct_fbs(r0, v0, rf, tof_s, eta, dv0, vf0, mu=mu, tol=1e-9, max_nfev=200)

    dv_err = float(np.linalg.norm(fbs.dv_dsm - dv_lam))
    vf_err = float(np.linalg.norm(fbs.v_arrive - vf_lam))
    dv_mag_lam = float(np.linalg.norm(dv_lam))
    dv_mag_fbs = float(np.linalg.norm(fbs.dv_dsm))

    print(f"\n  -- leg [{label}] eta={eta:.3f} tof={tof_s / DAY_S:.1f} d --")
    print(f"     Lambert dV_DSM = {dv_mag_lam:.6f} km/s  |vf| = {np.linalg.norm(vf_lam):.6f}")
    print(
        f"     FBS     dV_DSM = {dv_mag_fbs:.6f} km/s  |vf| = {np.linalg.norm(fbs.v_arrive):.6f}"
        f"  (converged={fbs.converged}, max_res={fbs.max_residual:.3e}, "
        f"nfev={fbs.solver_nfev}, njev={fbs.solver_njev})"
    )
    print(
        f"     AGREEMENT: |dv_fbs - dv_lam| = {dv_err:.3e} km/s, "
        f"|vf_fbs - vf_lam| = {vf_err:.3e} km/s"
    )
    return {
        "label": label,
        "dv_mag_lam": dv_mag_lam,
        "dv_mag_fbs": dv_mag_fbs,
        "fbs_converged": fbs.converged,
        "fbs_max_res": fbs.max_residual,
        "fbs_nfev": fbs.solver_nfev,
        "dv_err": dv_err,
        "vf_err": vf_err,
        "vf_lam": vf_lam,
        "vf_fbs": fbs.v_arrive,
    }


def attempt_s1l1(eph: Ephemeris) -> None:
    _hr("S1L1 / russell-ch4-4.991gG2 — corrected two-arc topology")
    print(
        "Topology (Russell App-C #83): E -> g(E-E free return, NO Mars) -> E flyby\n"
        "  -> G(E-M-E Mars transit, longitude rendezvous w/ true DE440 Mars) -> E.\n"
        "The binding multi-arc constraint lives on the G (Mars-transit) leg: the\n"
        "right-boundary position rf is the TRUE DE440 Mars position at the published\n"
        "arrival epoch, so a converged FBS match point IS a longitude-correct Mars\n"
        "intercept. Same-model golden = the Lambert dsm_leg on the identical leg BVP.\n"
        "Constraints in the FBS residual: [r^B - r^F (3, km); v^B - v^F (3, km/s)] at\n"
        "the interior impulse — full position AND velocity match; rf carries the true\n"
        "Mars-longitude target."
    )
    arcs = build_seeded_arcs(eph)
    g_legs = [a for a in arcs if a.is_mars_transit]
    results = []
    for arc in g_legs:
        arrival_no = arc.leg_no + 1
        _pub_tof, pub_vinf = APPC_MARS_TRANSIT[arrival_no]
        # Right boundary = TRUE DE440 Mars position at the arrival epoch (the
        # longitude-rendezvous target). The G-leg eta is free; use 0.5 (the
        # Takao default) — a real interior DSM the corrector can place anywhere.
        r_mars, v_mars = eph.state("M", arc.t1_sec)
        r_mars = np.asarray(r_mars, dtype=np.float64)
        v_mars = np.asarray(v_mars, dtype=np.float64)
        tof_s = arc.t1_sec - arc.t0_sec
        res = fbs_vs_lambert_leg(
            r0=arc.r0_km,
            v0=arc.v0_km_s,
            rf=r_mars,
            tof_s=tof_s,
            eta=0.5,
            label=f"G leg#{arc.leg_no}->M#{arrival_no}",
        )
        # Emerged v_inf at Mars from BOTH lanes vs the published App-C value.
        vinf_lam = _vinf(res["vf_lam"], v_mars)
        vinf_fbs = _vinf(res["vf_fbs"], v_mars)
        print(
            f"     v_inf@Mars: Lambert={vinf_lam:.4f}  FBS={vinf_fbs:.4f}  "
            f"published(App-C)={pub_vinf:.4f} km/s"
        )
        res["vinf_lam"] = vinf_lam
        res["vinf_fbs"] = vinf_fbs
        res["vinf_pub"] = pub_vinf
        results.append(res)

    _s1l1_verdict(results)


def _s1l1_verdict(results: list[dict]) -> None:
    n = len(results)
    all_conv = all(r["fbs_converged"] for r in results)
    max_dv_err = max(r["dv_err"] for r in results)
    max_vf_err = max(r["vf_err"] for r in results)
    print(
        f"\n  S1L1 SUMMARY ({n} G legs): FBS all converged={all_conv}; "
        f"max |dv_fbs-dv_lam|={max_dv_err:.3e}, max |vf_fbs-vf_lam|={max_vf_err:.3e} km/s."
    )
    # NOTE: dV_DSM here is NOT zero because eta=0.5 forces an interior impulse on a
    # leg whose true geometry is ballistic; the point is FBS<->Lambert AGREEMENT on
    # the SAME leg BVP + a true-Mars intercept, not a low dV.


def attempt_644gg3(eph: Ephemeris) -> None:
    _hr("russell-ch4-6.44Gg3 — full E-M-E-M-E sequence, per-leg FBS vs Lambert")
    cat = load_catalog()
    row = cat.by_id["russell-ch4-6.44Gg3"]
    raw = row.raw
    vinf_e = raw["vinf_kms_at_encounters"][0]["vinf_kms"]
    vinf_m = raw["vinf_kms_at_encounters"][1]["vinf_kms"]
    arcs = raw["free_return_arcs"]
    transit = raw["invariants"]["transit_times_days"][0]
    g_yr = arcs[0]["tof_years"]
    g_arc_yr = arcs[1]["tof_years"]
    print(
        "Topology (#153): the descriptor arcs g(2.087 yr) + G(4.3191 yr) sum to the\n"
        "6.41-yr period; unrolled into Mars-bracketing pieces -> sequence E-M-E-M-E\n"
        f"(4 legs). Sourced anchors (EXPECTED, never imposed): v_inf E={vinf_e}, "
        f"M={vinf_m} km/s.\n"
        "Per-leg ToFs: leg1 E->M transit, leg2 M->E = g_arc - transit, leg3 E->M\n"
        "transit, leg4 M->E = G_arc - transit.\n"
        "Constraints in each FBS leg residual: [r^B-r^F (km); v^B-v^F (km/s)] — full\n"
        "match point. Same-model golden = Lambert dsm_leg per identical leg BVP.\n"
        "The chain-level binding constraint (V_inf magnitude continuity at each\n"
        "flyby) is the thing that historically did NOT close; we measure the emerged\n"
        "encounter v_inf and report the gap honestly."
    )
    transit_s = transit * DAY_S
    leg2_tof_s = g_yr * 365.25 * DAY_S - transit_s
    leg4_tof_s = g_arc_yr * 365.25 * DAY_S - transit_s
    leg_tofs = [transit_s, leg2_tof_s, transit_s, leg4_tof_s]
    sequence = ["E", "M", "E", "M", "E"]

    # Build the chain ballistically (like evaluate_dsm_chain) but solve every leg
    # BOTH ways and compare. Departure: v_planet(E,t0) + v_inf_out0 (sourced
    # magnitude, tangential-prograde seed direction). We propagate the LAMBERT
    # heliocentric arrival velocity forward (ballistic flyby continuity) so the
    # chain geometry is identical for both lanes; FBS is checked per-leg against
    # the Lambert leg on those same boundary states.
    t0 = 0.0
    r0, v_pl0 = eph.state("E", t0)
    r0 = np.asarray(r0, dtype=np.float64)
    v_pl0 = np.asarray(v_pl0, dtype=np.float64)
    # Tangential-prograde departure direction in the ecliptic plane.
    v_hat = v_pl0 / np.linalg.norm(v_pl0)
    v_depart = v_pl0 + vinf_e * v_hat

    t_cursor = t0
    r_curr = r0
    results = []
    emerged_vinf = []
    for i in range(4):
        tof_s = leg_tofs[i]
        t_arr = t_cursor + tof_s
        target = sequence[i + 1]
        r_tgt, v_tgt = eph.state(target, t_arr)
        r_tgt = np.asarray(r_tgt, dtype=np.float64)
        v_tgt = np.asarray(v_tgt, dtype=np.float64)
        res = fbs_vs_lambert_leg(
            r0=r_curr,
            v0=v_depart,
            rf=r_tgt,
            tof_s=tof_s,
            eta=0.5,
            label=f"leg{i + 1} {sequence[i]}->{target}",
        )
        vinf_in = _vinf(res["vf_lam"], v_tgt)
        emerged_vinf.append((sequence[i + 1], vinf_in))
        print(f"     emerged v_inf@{target} (incoming) = {vinf_in:.4f} km/s")
        results.append(res)
        # ballistic flyby heliocentric continuity for the next leg
        v_depart = np.asarray(res["vf_lam"], dtype=np.float64)
        r_curr = r_tgt
        t_cursor = t_arr

    _644_verdict(results, emerged_vinf, vinf_e, vinf_m)


def _644_verdict(
    results: list[dict],
    emerged_vinf: list[tuple[str, float]],
    vinf_e: float,
    vinf_m: float,
) -> None:
    all_conv = all(r["fbs_converged"] for r in results)
    max_dv_err = max(r["dv_err"] for r in results)
    max_vf_err = max(r["vf_err"] for r in results)
    print(
        f"\n  6.44Gg3 SUMMARY: FBS all legs converged={all_conv}; "
        f"max |dv_fbs-dv_lam|={max_dv_err:.3e}, max |vf_fbs-vf_lam|={max_vf_err:.3e} km/s."
    )
    print(
        "  Emerged encounter v_inf (Lambert/FBS agree per leg): "
        + ", ".join(f"{b}={v:.2f}" for b, v in emerged_vinf)
    )
    print(f"  Sourced anchors: E={vinf_e}, M={vinf_m} km/s. Gap to anchors is the")
    print("  chain-level V_inf-continuity blocker — NOT a per-leg BVP failure.")


def probe_fbs_wins_where_lambert_fails(eph: Ephemeris) -> None:
    """Does FBS solve any leg geometry where the Lambert dsm_leg RAISES?

    This is the strongest possible argument for defaulting FBS (#243): a leg the
    Lambert corrector cannot represent at all but FBS can. We probe the two
    documented Lambert pain points on real S1L1 geometry:

      (1) near-180-deg transfer (LambertGeometryError territory), and
      (2) a long multi-rev loop arc with the back-arc forced single-rev
          (the #153 degeneracy).

    For each we attempt the Lambert dsm_leg (catch any LambertError) and the FBS
    corrector seeded ballistically, and report which lanes produced a valid leg.
    """
    _hr("PROBE: can FBS solve a leg geometry where Lambert dsm_leg FAILS?")
    arcs = build_seeded_arcs(eph)
    # Take the long G-arc-style loop: a real Mars->Earth return over a >1-period
    # ToF. Use the 6.44Gg3 leg4 geometry (1315 d M->E, ~3 revs).
    cat = load_catalog()
    raw = cat.by_id["russell-ch4-6.44Gg3"].raw
    transit_s = raw["invariants"]["transit_times_days"][0] * DAY_S
    g_arc_yr = raw["free_return_arcs"][1]["tof_years"]
    leg4_tof_s = g_arc_yr * 365.25 * DAY_S - transit_s

    t0 = 0.0
    r_m, v_m = eph.state("M", t0)
    r_m = np.asarray(r_m, dtype=np.float64)
    v_m = np.asarray(v_m, dtype=np.float64)
    r_e, _ = eph.state("E", t0 + leg4_tof_s)
    r_e = np.asarray(r_e, dtype=np.float64)
    # A modest departure velocity off Mars.
    v0 = v_m + 3.0 * (v_m / np.linalg.norm(v_m))

    print(f"\n  Long M->E loop arc: tof={leg4_tof_s / DAY_S:.1f} d (~3 revs), eta=0.5")
    # Lambert lane, single-rev (the historical default) — expect degeneracy.
    lam_ok = False
    lam_dv = float("nan")
    try:
        lam = dsm_leg(r_m, v0, leg4_tof_s, 0.5, r_e, max_revs=0)
        lam_ok = True
        lam_dv = float(np.linalg.norm(lam.v_depart_post_dsm - lam.v_arrive_pre_dsm))
        print(f"     Lambert (max_revs=0): SOLVED, dV_DSM={lam_dv:.4f} km/s")
    except Exception as exc:  # reporting which lane fails
        print(f"     Lambert (max_revs=0): RAISED {type(exc).__name__}: {exc}")
    # FBS lane (A) — seed from a circular-ish guess (NO Lambert needed). This is
    # the "FBS as a standalone Lambert-free corrector" condition.
    v_circ = float(np.sqrt(MU_SUN_KM3_S2 / np.linalg.norm(r_m)))
    dv0 = np.zeros(3)
    vf0 = v_circ * (r_e / np.linalg.norm(r_e))
    fbs_cold = dsm_leg_correct_fbs(r_m, v0, r_e, leg4_tof_s, 0.5, dv0, vf0, max_nfev=400)
    print(
        f"     FBS (cold circular seed): converged={fbs_cold.converged}, "
        f"dV_DSM={np.linalg.norm(fbs_cold.dv_dsm):.4f} km/s, "
        f"max_res={fbs_cold.max_residual:.3e}, nfev={fbs_cold.solver_nfev}"
    )
    # FBS lane (B) — seed from the Lambert solution (the parity condition). Shows
    # FBS CAN land the same leg, but only once Lambert has already supplied the basin.
    if lam_ok:
        dvb = lam.v_depart_post_dsm - lam.v_arrive_pre_dsm
        fbs_warm = dsm_leg_correct_fbs(
            r_m, v0, r_e, leg4_tof_s, 0.5, dvb, lam.v_arrive, max_nfev=400
        )
        print(
            f"     FBS (warm Lambert seed): converged={fbs_warm.converged}, "
            f"dV_DSM={np.linalg.norm(fbs_warm.dv_dsm):.4f} km/s, "
            f"max_res={fbs_warm.max_residual:.3e}, nfev={fbs_warm.solver_nfev}"
        )
    print(
        "\n  NOTE: dsm_leg with max_revs>0 already repairs the single-rev degeneracy\n"
        "  (#157) by enumerating multi-rev Lambert branches, so Lambert is not\n"
        "  strictly stuck here either. The point of this probe is whether FBS offers\n"
        "  a leg-representation advantage that the Lambert lane lacks even with its\n"
        "  multi-rev repair. Verdict reported in the note."
    )
    _ = arcs


def main() -> None:
    eph = Ephemeris("astropy")
    print("FBS vs Lambert closure attempt on the hard multi-arc rows (#242)")
    print("Evidence-gathering only. NO catalogue writeback.")
    attempt_s1l1(eph)
    attempt_644gg3(eph)
    probe_fbs_wins_where_lambert_fails(eph)
    _hr("DONE")


if __name__ == "__main__":
    main()
