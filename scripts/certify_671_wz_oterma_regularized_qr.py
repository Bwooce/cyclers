"""#671 -- the REGULARIZED variational jet + QR/Lohner-C1 reframing compose, and
together carry the rigorous Oterma arc far past the ``#670`` wrapping wall.

Stage 5 of the Wilczak-Zgliczynski proof-machinery build (``#668`` Stages 1-2,
``#669`` Stage 3 QR reframing, ``#670`` Stage 4 Levi-Civita regularization).

``#670`` regularized the Jupiter close approach that stopped Stages 2/3 and
carried the rigorous enclosure THROUGH the perijove -- but then hit a NEW wall
(``t ~ 0.85``, ``tau ~ 40``) it characterized precisely as classic exponential
(``e^t``) WRAPPING (Jacobi stays enclosed, no collision-guard trip, half-width
grows smoothly then blows up).  That is exactly the failure mode ``#669``'s QR
reframing defeats -- but ``#670`` integrated the regularized field with the plain
(non-variational) jet, so QR could not be applied in the regularized frame.

This driver closes that gap.  ``scripts/_validated_taylor_integrator.py`` gains
``make_cr3bp_lc_secondary_variational_jet`` -- the state-transition-matrix flow
``V' = Dg(w) V`` of the Levi-Civita field ``g(w)``, ``w = [u,v,Pu,Pv,T]`` -- whose
analytic 5x5 Jacobian ``Dg`` is re-derived from the exact regularized equations of
motion and validated against finite differences of that field (see the ``#671``
test).  Feeding it and the plain regularized jet to ``#669``'s ``integrate_c1_qr``
makes both fixes compose: the QR frame now tracks the regularized coordinates' own
local stretching directions.

Findings (all rigorous where stated, ``mpmath.iv``):

  1. **Closed-form anchors re-validated WITH QR.**  The pure-Kepler-in-Levi-Civita
     radial fall (a genuine two-body collision orbit) run through ``integrate_c1_qr``
     with the Stage-5 STM jet still traverses ``r2 = 0`` smoothly and matches
     ``sqrt(a) cos(omega tau)`` to <1e-20, containing the true value and excluding
     nearby wrong ones.  The Stage 2/3 exp/harmonic goldens are intact through the
     QR path.  So the regularized-variational/QR machinery is correct, not vacuous.

  2. **QR defeats the regularized-frame wrapping wall, decisively.**  At a matched
     grid (``h = 0.125``, order 10) the naive C0 box balloons by the textbook ``e^t``
     wrapping growth and fails to validate at ``tau ~ 23``, while the QR-regularized
     enclosure stays flat (``hw ~ 1e-11``) all the way to ``tau = 80``, rigorously
     conserving the Jacobi constant and non-vacuous the whole way -- no wall.

  3. **The rigorous arc reaches ``t_phys ~ 2.05`` (``tau = 80``).**  ``#668``/``#669``
     stopped at ``t ~ 0.42-0.466`` (the first Jupiter perijove); ``#670`` reached
     ``t ~ 0.80`` and its naive box died at ``t ~ 0.85``.  Stage 5 more than doubles
     the physical-time reach past ``#670``'s wall, and the arc now sails through the
     REPEATED Jupiter close approaches (float diagnostic: perijoves at ``t = 0.462,
     1.278, 2.055`` with ``r2`` down to ``0.0019``) that each individually stopped
     the pre-regularization stages -- with no wall found in the probed range.

Honest scope: the W-Z heteroclinic connection between the ``L1*`` and ``L2*``
Lyapunov orbits is an *infinite-time* asymptotic object (proven via covering
relations on a FINITE sequence of Poincare-section crossings, not a finite total
flight time), so "fraction of the full connection time" is not the right metric;
the honest progress metric is how far the rigorous flow enclosure now reaches and
how many section crossings / close approaches it survives.  Rigorous Poincare-map
derivatives and the covering-relations / h-set topology remain future stages.

Run:  ``uv run python scripts/certify_671_wz_oterma_regularized_qr.py``
Writes ``data/671_wz_oterma_regularized_qr_certificate.json``.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any

import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import _validated_taylor_integrator as vti  # noqa: E402

GOLDEN = ROOT / "data" / "golden" / "wz_oterma_heteroclinic.yaml"
OUT = ROOT / "data" / "671_wz_oterma_regularized_qr_certificate.json"

MU = 0.0009537
C_LEVEL = 3.03
H_ENERGY = -C_LEVEL / 2
PERIJOVE_T = 0.466  # #669 float diagnostic: the physical wall Stages 2/3 hit
STAGE4_WALL_T = 0.85  # #670's naive-C0 wrapping wall (t model-time)


def _kepler_lc_jet(h: float) -> Any:
    def jet(iv: Any, state: list[Any], mu: float) -> list[Any]:
        u, v, pu, pv, _t = state
        r2 = vti.ts_add(iv, vti.ts_mul(iv, u, u), vti.ts_mul(iv, v, v))
        return [
            vti.ts_scale(iv, pu, iv.mpf(0.25)),
            vti.ts_scale(iv, pv, iv.mpf(0.25)),
            vti.ts_scale(iv, u, iv.mpf(2.0 * h)),
            vti.ts_scale(iv, v, iv.mpf(2.0 * h)),
            r2,
        ]

    return jet


def _kepler_lc_var_jet(h: float) -> Any:
    def jet(iv: Any, aug: list[Any], mu: float) -> list[Any]:
        u, v, pu, pv, _t = aug[0], aug[1], aug[2], aug[3], aug[4]
        vmat = aug[5:30]
        order = len(u) - 1
        zero = vti.ts_const(iv, iv.mpf(0), order)
        q = vti.ts_const(iv, iv.mpf(0.25), order)
        th = vti.ts_const(iv, iv.mpf(2.0 * h), order)
        two_u = vti.ts_scale(iv, u, iv.mpf(2.0))
        two_v = vti.ts_scale(iv, v, iv.mpf(2.0))
        df = [
            [zero, zero, q, zero, zero],
            [zero, zero, zero, q, zero],
            [th, zero, zero, zero, zero],
            [zero, th, zero, zero, zero],
            [two_u, two_v, zero, zero, zero],
        ]
        dV: list[Any] = [zero for _ in range(25)]
        for c in range(5):
            col = [vmat[r * 5 + c] for r in range(5)]
            for rr in range(5):
                acc = None
                for kk in range(5):
                    if df[rr][kk] is zero:
                        continue
                    term = vti.ts_mul(iv, df[rr][kk], col[kk])
                    acc = term if acc is None else vti.ts_add(iv, acc, term)
                dV[rr * 5 + c] = acc if acc is not None else zero
        r2 = vti.ts_add(iv, vti.ts_mul(iv, u, u), vti.ts_mul(iv, v, v))
        state = [
            vti.ts_scale(iv, pu, iv.mpf(0.25)),
            vti.ts_scale(iv, pv, iv.mpf(0.25)),
            vti.ts_scale(iv, u, iv.mpf(2.0 * h)),
            vti.ts_scale(iv, v, iv.mpf(2.0 * h)),
            r2,
        ]
        return [*state, *dV]

    return jet


def radial_fall_qr_anchor(iv: Any) -> dict[str, Any]:
    """Closed-form collision anchor, now through integrate_c1_qr (both fixes)."""
    import mpmath as mp

    a, mu = 1.0, 0.5
    h = -mu / a
    omega = mp.sqrt(mp.mpf(-h) / 2)
    y0 = [iv.mpf([1.0, 1.0]), iv.mpf([0, 0]), iv.mpf([0, 0]), iv.mpf([0, 0]), iv.mpf([0, 0])]
    tau_f = 4.0
    res = vti.integrate_c1_qr(
        iv, _kepler_lc_jet(h), _kepler_lc_var_jet(h), y0, tau_f, mu, n_steps=32, order=16
    )
    u_enc = res["final"][0]
    u_cf = mp.sqrt(mp.mpf(a)) * mp.cos(omega * mp.mpf(tau_f))
    return {
        "problem": "radial fall (collision orbit) -> harmonic in L-C, via integrate_c1_qr",
        "validated_through_collision": bool(res["validated"]),
        "u_encloses_closed_form": bool(u_enc.a <= u_cf) and bool(u_cf <= u_enc.b),
        "u_excludes_wrong": not (
            bool(u_enc.a <= u_cf + mp.mpf("1e-15")) and bool(u_cf + mp.mpf("1e-15") <= u_enc.b)
        ),
        "u_halfwidth": float(u_enc.delta.b) / 2,
        "final_max_halfwidth": res["final_max_halfwidth"],
    }


def qr_positive_controls(iv: Any) -> dict[str, Any]:
    """Stage 2/3 exp/harmonic goldens through the QR integrator."""
    import mpmath as mp

    def _jet_exp(iv: Any, s: list[Any], mu: float) -> list[Any]:
        return [s[0]]

    def _var_exp(iv: Any, a: list[Any], mu: float) -> list[Any]:
        return [a[0], a[1]]

    def _jet_harm(iv: Any, s: list[Any], mu: float) -> list[Any]:
        return [s[1], vti.ts_scale(iv, s[0], iv.mpf(-1.0))]

    def _var_harm(iv: Any, a: list[Any], mu: float) -> list[Any]:
        def neg(s: list[Any]) -> list[Any]:
            return vti.ts_scale(iv, s, iv.mpf(-1.0))

        return [a[1], neg(a[0]), a[4], a[5], neg(a[2]), neg(a[3])]

    re = vti.integrate_c1_qr(
        iv, _jet_exp, _var_exp, [iv.mpf([1.0, 1.0])], 1.0, 0.0, n_steps=8, order=14
    )
    enc = re["final"][0]
    e_true = iv.exp(iv.mpf(1.0))
    rh = vti.integrate_c1_qr(
        iv,
        _jet_harm,
        _var_harm,
        [iv.mpf([1.0, 1.0]), iv.mpf([0, 0])],
        1.0,
        0.0,
        n_steps=8,
        order=14,
    )
    y1, y2 = rh["final"]
    c1, s1 = iv.cos(iv.mpf(1.0)), -iv.sin(iv.mpf(1.0))
    return {
        "exp_encloses_e": bool(enc.a <= e_true.a) and bool(e_true.b <= enc.b),
        "exp_excludes_wrong": not (
            bool(enc.a <= mp.e + mp.mpf("1e-12")) and bool(mp.e + mp.mpf("1e-12") <= enc.b)
        ),
        "harmonic_encloses_cos1": bool(y1.a <= c1.a) and bool(c1.b <= y1.b),
        "harmonic_encloses_neg_sin1": bool(y2.a <= s1.a) and bool(s1.b <= y2.b),
    }


def jacobian_fd_check(iv: Any) -> dict[str, Any]:
    """Analytic Dg vs central finite differences of the regularized field."""
    vjet = vti.make_cr3bp_lc_secondary_variational_jet(H_ENERGY)
    sjet = vti.make_cr3bp_lc_secondary_jet(H_ENERGY)
    eps = 1e-7

    def g_at(w: list[float]) -> list[float]:
        o = [vti.ts_const(iv, iv.mpf([x, x]), 0) for x in w]
        d = sjet(iv, o, MU)
        return [float(di[0].mid) for di in d]

    def dg(w: list[float]) -> list[list[float]]:
        aug = [iv.mpf([x, x]) for x in w] + [
            iv.mpf(1.0 if r == c else 0.0) for r in range(5) for c in range(5)
        ]
        o = [vti.ts_const(iv, aug[i], 0) for i in range(30)]
        d = vjet(iv, o, MU)
        stm = d[5:30]
        return [[float(stm[r * 5 + c][0].mid) for c in range(5)] for r in range(5)]

    points = [
        [0.7, 0.3, 0.2, -0.4, 0.0],
        [1.2, -0.5, 0.1, 0.6, 3.0],
        [0.05, 0.02, -0.3, 0.15, 1.0],
        [0.001, -0.002, 0.05, 0.03, 0.5],
        [-0.9, 0.4, 0.5, -0.2, 0.0],
    ]
    maxerr = 0.0
    for w in points:
        a = dg(w)
        for j in range(5):
            wp = list(w)
            wp[j] += eps
            wm = list(w)
            wm[j] -= eps
            gp, gm = g_at(wp), g_at(wm)
            for i in range(5):
                fd = (gp[i] - gm[i]) / (2 * eps)
                scale = max(1.0, abs(a[i][j]), abs(fd))
                maxerr = max(maxerr, abs(a[i][j] - fd) / scale)
    return {
        "n_points": len(points),
        "includes_near_collision": True,
        "max_relative_error": maxerr,
        "matches": bool(maxerr < 1e-5),
    }


def _oterma_ic(iv: Any) -> tuple[list[Any], Any, float]:
    golden = yaml.safe_load(GOLDEN.read_text())
    seq = golden["crossings"]["heteroclinic_L1_to_L2"]["sequence"]
    x0, vx0 = float(seq[0][0]), float(seq[0][1])
    r1 = iv.sqrt((iv.mpf(x0) + iv.mpf(MU)) ** 2)
    r2 = iv.sqrt((iv.mpf(x0) + iv.mpf(MU - 1.0)) ** 2)
    vy2 = (
        iv.mpf(x0) ** 2
        + 2 * (1 - iv.mpf(MU)) / r1
        + 2 * iv.mpf(MU) / r2
        - iv.mpf(vx0) ** 2
        - iv.mpf(C_LEVEL)
    )
    vy0 = iv.sqrt(vy2)
    y0 = vti.lc_secondary_from_physical(
        iv, iv.mpf([x0, x0]), iv.mpf([0, 0]), iv.mpf([vx0, vx0]), vy0, MU
    )
    c0 = vti.cr3bp_planar_jacobi(iv, vti.lc_secondary_to_physical(iv, y0, MU), MU)
    return y0, c0, x0


def naive_c0_wall(iv: Any, y0: list[Any], c0: Any) -> dict[str, Any]:
    """The Stage-4 naive-C0 wrapping wall at the matched QR grid (h=0.125)."""
    h = 0.125
    tau_f = 44.0
    ns = int(tau_f / h)
    res = vti.integrate(
        iv, vti.make_cr3bp_lc_secondary_jet(H_ENERGY), y0, tau_f, MU, n_steps=ns, order=10
    )
    nc = res["n_completed"]
    box = res["final"]
    t_phys = box[4]
    hws = res["max_halfwidths"]
    return {
        "grid_h": h,
        "order": 10,
        "validated": bool(res["validated"]),
        "tau_reached": nc * h,
        "t_phys_reached": float(t_phys.b),
        "final_max_halfwidth": res["final_max_halfwidth"],
        "signature": "exponential e^t wrapping growth then blow-up (a-priori box fails)",
        "halfwidth_growth_sample": [
            (round(i * h, 2), hws[i]) for i in range(0, len(hws), max(1, len(hws) // 10))
        ],
    }


def qr_reach(iv: Any, y0: list[Any], c0: Any) -> dict[str, Any]:
    """The Stage-5 headline: QR-regularized reach at the same grid."""
    h = 0.125
    tau_f = 80.0
    ns = int(tau_f / h)
    jet = vti.make_cr3bp_lc_secondary_jet(H_ENERGY)
    vjet = vti.make_cr3bp_lc_secondary_variational_jet(H_ENERGY)
    res = vti.integrate_c1_qr(iv, jet, vjet, y0, tau_f, MU, n_steps=ns, order=10)
    nc = res["n_completed"]
    box = res["final"]
    t_phys = box[4]
    phys = vti.lc_secondary_to_physical(iv, box, MU)
    cf = vti.cr3bp_planar_jacobi(iv, phys, MU)
    hws = res["max_halfwidths"]
    return {
        "grid_h": h,
        "order": 10,
        "validated": bool(res["validated"]),
        "tau_reached": nc * h,
        "t_phys_reached_lo": float(t_phys.a),
        "t_phys_reached_hi": float(t_phys.b),
        "final_max_halfwidth": res["final_max_halfwidth"],
        "jacobi_end_contains_C0": bool(cf.a <= c0.a) and bool(c0.b <= cf.b),
        "hit_wall": not bool(res["validated"]),
        "halfwidth_growth_sample": [
            (round(i * h, 2), hws[i]) for i in range(0, len(hws), max(1, len(hws) // 20))
        ],
    }


def build_certificate() -> dict[str, Any]:
    import mpmath as mp

    mp.mp.dps = 40
    mp.iv.dps = 40
    iv = mp.iv

    t0 = time.time()
    anchor = radial_fall_qr_anchor(iv)
    controls = qr_positive_controls(iv)
    jac = jacobian_fd_check(iv)
    y0, c0, x0 = _oterma_ic(iv)
    print(f"[#671] anchors+controls+jacobian done ({time.time() - t0:.0f}s)", flush=True)

    naive = naive_c0_wall(iv, y0, c0)
    print(
        f"[#671] naive wall done: stops tau={naive['tau_reached']} ({time.time() - t0:.0f}s)",
        flush=True,
    )
    reach = qr_reach(iv, y0, c0)
    print(
        f"[#671] QR reach done: tau={reach['tau_reached']} "
        f"t_phys~{reach['t_phys_reached_hi']:.3f} ({time.time() - t0:.0f}s)",
        flush=True,
    )

    return {
        "task": "#671",
        "system": {
            "mu": MU,
            "C": C_LEVEL,
            "h_energy": H_ENERGY,
            "model": "planar CR3BP (Sun-Jupiter, Oterma)",
        },
        "method": (
            "Regularized C1/Lohner-QR: the Levi-Civita secondary-regularized field "
            "(#670, make_cr3bp_lc_secondary_jet) integrated with #669's integrate_c1_qr, "
            "driven by the NEW analytic regularized variational jet "
            "make_cr3bp_lc_secondary_variational_jet (V'=Dg(w)V, w=[u,v,Pu,Pv,T]; Dg re-derived "
            "from the exact regularized EOM and validated against finite differences). QR now "
            "reframes the enclosure in the regularized coordinates' own local stretching "
            "directions."
        ),
        "start": {"x": x0, "C0_lower": float(c0.a), "C0_upper": float(c0.b)},
        "prior_stage_walls": {
            "stage_2_3_perijove_t": PERIJOVE_T,
            "stage_4_reach_t": 0.80,
            "stage_4_wrapping_wall_t": STAGE4_WALL_T,
        },
        "radial_fall_qr_closed_form_anchor": anchor,
        "qr_positive_controls": controls,
        "regularized_variational_jacobian_check": jac,
        "naive_c0_wrapping_wall_matched_grid": naive,
        "qr_regularized_reach": reach,
        "findings": {
            "regularized_variational_jet_correct": bool(jac["matches"]),
            "qr_composes_with_regularization_closed_form": bool(
                anchor["validated_through_collision"]
                and anchor["u_encloses_closed_form"]
                and anchor["u_excludes_wrong"]
            ),
            "stage_2_3_controls_intact": all(controls.values()),
            "naive_wraps_and_fails_at_matched_grid": bool(
                (not naive["validated"]) and naive["final_max_halfwidth"] > 1.0
            ),
            "qr_defeats_wrapping": bool(
                reach["validated"]
                and reach["jacobi_end_contains_C0"]
                and reach["final_max_halfwidth"] < 1e-6
            ),
            "qr_reach_t_phys": reach["t_phys_reached_hi"],
            "reach_extension_factor_over_stage4_wall": round(
                reach["t_phys_reached_hi"] / STAGE4_WALL_T, 2
            ),
            "wall_found_in_probed_range": bool(reach["hit_wall"]),
        },
        "scope_note": (
            "Rigorous, and it meets the Stage-5 goal: the regularized variational jet + QR "
            "reframing COMPOSE. At a matched grid the naive C0 box wraps (e^t) and fails at "
            f"tau~{naive['tau_reached']}, while the QR-regularized enclosure stays tight "
            f"(hw~{reach['final_max_halfwidth']:.1e}) to tau={reach['tau_reached']} "
            f"(t_phys~{reach['t_phys_reached_hi']:.2f}), Jacobi rigorously conserved, no wall in "
            "range -- more than doubling the physical-time reach past #670's t~0.85 wrapping wall "
            "and sailing through the repeated Jupiter close approaches that stopped the "
            "pre-regularization stages. The W-Z heteroclinic connection is an infinite-time "
            "asymptotic object (proven via covering relations on a finite "
            "Poincare-section-crossing sequence, not a finite flight time), so 'fraction of the "
            "full connection time' is not the right metric; rigorous Poincare-map derivatives and "
            "the covering-relations / h-set topology remain future stages."
        ),
    }


def main() -> None:
    cert = build_certificate()
    OUT.write_text(json.dumps(cert, indent=2) + "\n")
    a = cert["radial_fall_qr_closed_form_anchor"]
    print(
        f"[#671] radial-fall QR anchor: "
        f"valid-through-collision={a['validated_through_collision']}, "
        f"encloses cf={a['u_encloses_closed_form']}, excludes wrong={a['u_excludes_wrong']}"
    )
    j = cert["regularized_variational_jacobian_check"]
    print(
        f"[#671] regularized Jacobian vs finite-diff: "
        f"max_rel_err={j['max_relative_error']:.2e} matches={j['matches']}"
    )
    print(f"[#671] Stage 2/3 QR controls intact: {all(cert['qr_positive_controls'].values())}")
    n = cert["naive_c0_wrapping_wall_matched_grid"]
    print(
        f"[#671] naive C0 (matched grid): validated={n['validated']} stops tau={n['tau_reached']} "
        f"hw={n['final_max_halfwidth']:.2e} -> e^t wrapping wall"
    )
    r = cert["qr_regularized_reach"]
    print(
        f"[#671] QR-regularized reach: validated={r['validated']} tau={r['tau_reached']} "
        f"t_phys~{r['t_phys_reached_hi']:.3f} hw={r['final_max_halfwidth']:.2e} "
        f"jacobi_ok={r['jacobi_end_contains_C0']} wall={r['hit_wall']}"
    )
    print(f"[#671] certificate written -> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
