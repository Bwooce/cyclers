"""#672 -- rigorous Poincare SECTION-CROSSING MAP with its own rigorous Jacobian.

Stage 6 of the Wilczak-Zgliczynski (W-Z) proof-machinery build (``#668`` Stages
1-2, ``#669`` Stage 3 QR reframing, ``#670`` Stage 4 Levi-Civita regularization,
``#671`` Stage 5 regularized+QR composition -- the flow enclosure now survives 3
Jupiter perijoves to ``t ~ 2.05``).

The W-Z existence proof is built on covering relations between h-sets placed along
the Oterma orbit, verified through the Poincare SECTION MAP ``P`` (not the raw
flow).  Their section is ``Theta = {y = 0}``, parameterized by ``(x, xdot)``; in
the Levi-Civita regularized coordinates ``#670``/``#671`` carry, the physical
ordinate is ``y = eta = 2 u v``, so the section condition is ``sigma(w) = 2uv = 0``.
This driver builds and certifies the necessary PRECURSOR to h-sets: a rigorous
section-crossing map with a rigorous Jacobian.  Building the actual h-sets and
proving covering relations remains the NEXT stage, explicitly not this one.

What it establishes (all ``mpmath.iv``, dps 40):

  1. **Closed-form positive control -- 2D harmonic oscillator, section {y1=0}.**
     Start ``(1, 1)``: the crossing time (``3 pi/4``), post-map state ``(0, -sqrt2)``
     AND the FULL section-map Jacobian ``DP = [[0,0],[-1/sqrt2,-1/sqrt2]]`` are all
     hand-computable in closed form (rotation STM, ``f = [y2,-y1]``, ``Dsigma =
     [1,0]``, ``DP = Dphi - f (Dsigma Dphi)/(Dsigma f)``).  The rigorous machinery
     encloses every one of them, tightly, and excludes nearby wrong values -- an
     end-to-end INDEPENDENT check that the crossing isolation + implicit-function
     derivative correction are actually correct, not merely non-crashing.

  2. **Second control -- pure-Kepler-in-Levi-Civita, section {u=0} (5D regularized).**
     ``u(tau) = cos(omega tau)`` crosses ``0`` at ``tau = pi``; post ``pu = -2``;
     ``Dsigma.f = -1/2``; ``DP`` row-``u`` vanishes (on-section).  Exercises the 5D
     regularized state and the accumulated-STM path on a case with a known answer.

  3. **Real W-Z Oterma section point.**  From the first published golden crossing
     (``x0 = 0.95229...``, on ``{y=0}``), the machinery rigorously isolates the FIRST
     forward re-crossing of ``2uv = 0``: a tight crossing-time interval, the enclosed
     physical ``(x, xdot)``, a transversality certificate (``Dsigma.f != 0``), and a
     tight, non-vacuous rigorous ``DP`` -- with the QR-accumulated flow Jacobian
     staying tight straight through the first Jupiter perijove the crossing sits just
     beyond.  (A forward shot from the ROUNDED published IC does not track the true,
     unstable heteroclinic connection, so this crossing is not expected to reproduce
     golden point 2; matching the published sequence is the job of the future h-set /
     covering-relations stage, which places h-sets to bracket the true orbit.)

Run:  ``uv run python scripts/certify_672_wz_oterma_section_map.py``
Writes ``data/672_wz_oterma_section_map_certificate.json``.
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
OUT = ROOT / "data" / "672_wz_oterma_section_map_certificate.json"

MU = 0.0009537
C_LEVEL = 3.03
H_ENERGY = -C_LEVEL / 2


# --------------------------------------------------------------------------- #
# Control 1: 2D harmonic oscillator, section {y1 = 0}.                          #
# --------------------------------------------------------------------------- #
def _harm_jet(iv: Any, s: list[Any], mu: float) -> list[Any]:
    return [s[1], vti.ts_scale(iv, s[0], iv.mpf(-1.0))]


def _harm_var(iv: Any, aug: list[Any], mu: float) -> list[Any]:
    y1, y2 = aug[0], aug[1]
    v00, v01, v10, v11 = aug[2], aug[3], aug[4], aug[5]

    def neg(s: list[Any]) -> list[Any]:
        return vti.ts_scale(iv, s, iv.mpf(-1.0))

    return [y2, neg(y1), v10, v11, neg(v00), neg(v01)]


def _sig_y1(iv: Any, st: list[Any]) -> Any:
    return st[0]


def _grad_y1(iv: Any, st: list[Any]) -> list[Any]:
    return [iv.mpf(1), iv.mpf(0)]


def harmonic_control(iv: Any) -> dict[str, Any]:
    import mpmath as mp

    y0 = [iv.mpf([1.0, 1.0]), iv.mpf([1.0, 1.0])]
    res = vti.rigorous_section_map(
        iv,
        _harm_jet,
        _harm_var,
        y0,
        0.0,
        sigma_val=_sig_y1,
        sigma_grad=_grad_y1,
        tau_max=4.0,
        n_steps=64,
        order=16,
        tau_min=0.1,
    )
    ts = res["tau_star"]
    cs = res["crossing_state"]
    dp = res["section_jacobian"]
    tstar_cf = 3 * mp.pi / 4
    y2_cf = -mp.sqrt(2)
    off = -1 / mp.sqrt(2)

    def contains(enc: Any, val: Any) -> bool:
        return bool(enc.a <= val) and bool(val <= enc.b)

    return {
        "start": [1.0, 1.0],
        "section": "{y1 = 0}",
        "tau_star_lo": float(ts.a),
        "tau_star_hi": float(ts.b),
        "tau_star_width": float(ts.delta.b),
        "tau_star_encloses_3pi_over_4": contains(ts, tstar_cf),
        "crossing_y1_halfwidth": float(cs[0].delta.b) / 2,
        "crossing_y2_encloses_neg_sqrt2": contains(cs[1], y2_cf),
        "DP_encloses_closed_form": (
            contains(dp[0][0], mp.mpf(0))
            and contains(dp[0][1], mp.mpf(0))
            and contains(dp[1][0], off)
            and contains(dp[1][1], off)
        ),
        "DP_excludes_wrong": (
            not contains(dp[1][0], off + mp.mpf("1e-12"))
            and not contains(dp[1][1], off - mp.mpf("1e-12"))
        ),
        "DP_max_halfwidth": max(float(dp[i][j].delta.b) for i in range(2) for j in range(2)) / 2,
        "dsigma_dot_f_encloses_neg_sqrt2": contains(res["dsigma_dot_f"], -mp.sqrt(2)),
        "transversal": bool(res["transversal"]),
        "closed_form": {
            "tau_star": float(tstar_cf),
            "crossing_state": [0.0, float(y2_cf)],
            "DP": [[0.0, 0.0], [float(off), float(off)]],
            "dsigma_dot_f": float(-mp.sqrt(2)),
        },
    }


# --------------------------------------------------------------------------- #
# Control 2: pure-Kepler-in-Levi-Civita, section {u = 0} (5D regularized).      #
# --------------------------------------------------------------------------- #
def _kep_jet(h: float) -> Any:
    def jet(iv: Any, s: list[Any], mu: float) -> list[Any]:
        u, v, pu, pv, _t = s
        r2 = vti.ts_add(iv, vti.ts_mul(iv, u, u), vti.ts_mul(iv, v, v))
        return [
            vti.ts_scale(iv, pu, iv.mpf(0.25)),
            vti.ts_scale(iv, pv, iv.mpf(0.25)),
            vti.ts_scale(iv, u, iv.mpf(2.0 * h)),
            vti.ts_scale(iv, v, iv.mpf(2.0 * h)),
            r2,
        ]

    return jet


def _kep_var(h: float) -> Any:
    def jet(iv: Any, aug: list[Any], mu: float) -> list[Any]:
        u, v, pu, pv, _t = aug[0], aug[1], aug[2], aug[3], aug[4]
        vmat = aug[5:30]
        order = len(u) - 1
        zero = vti.ts_const(iv, iv.mpf(0), order)
        q = vti.ts_const(iv, iv.mpf(0.25), order)
        th = vti.ts_const(iv, iv.mpf(2.0 * h), order)
        tu = vti.ts_scale(iv, u, iv.mpf(2.0))
        tv = vti.ts_scale(iv, v, iv.mpf(2.0))
        df = [
            [zero, zero, q, zero, zero],
            [zero, zero, zero, q, zero],
            [th, zero, zero, zero, zero],
            [zero, th, zero, zero, zero],
            [tu, tv, zero, zero, zero],
        ]
        dV: list[Any] = [zero for _ in range(25)]
        for c in range(5):
            col = [vmat[r * 5 + c] for r in range(5)]
            for rr in range(5):
                acc = None
                for kk in range(5):
                    if df[rr][kk] is zero:
                        continue
                    t = vti.ts_mul(iv, df[rr][kk], col[kk])
                    acc = t if acc is None else vti.ts_add(iv, acc, t)
                dV[rr * 5 + c] = acc if acc is not None else zero
        r2 = vti.ts_add(iv, vti.ts_mul(iv, u, u), vti.ts_mul(iv, v, v))
        st = [
            vti.ts_scale(iv, pu, iv.mpf(0.25)),
            vti.ts_scale(iv, pv, iv.mpf(0.25)),
            vti.ts_scale(iv, u, iv.mpf(2.0 * h)),
            vti.ts_scale(iv, v, iv.mpf(2.0 * h)),
            r2,
        ]
        return [*st, *dV]

    return jet


def _sig_u(iv: Any, st: list[Any]) -> Any:
    return st[0]


def _grad_u(iv: Any, st: list[Any]) -> list[Any]:
    return [iv.mpf(1), iv.mpf(0), iv.mpf(0), iv.mpf(0), iv.mpf(0)]


def kepler_lc_control(iv: Any) -> dict[str, Any]:
    import mpmath as mp

    h = -0.5
    y0 = [iv.mpf([1.0, 1.0]), iv.mpf([0, 0]), iv.mpf([0, 0]), iv.mpf([0, 0]), iv.mpf([0, 0])]
    res = vti.rigorous_section_map(
        iv,
        _kep_jet(h),
        _kep_var(h),
        y0,
        0.0,
        sigma_val=_sig_u,
        sigma_grad=_grad_u,
        tau_max=5.0,
        n_steps=64,
        order=16,
        tau_min=0.5,
    )
    ts = res["tau_star"]
    cs = res["crossing_state"]
    dp = res["section_jacobian"]

    def contains(enc: Any, val: Any) -> bool:
        return bool(enc.a <= val) and bool(val <= enc.b)

    return {
        "start_u": 1.0,
        "section": "{u = 0}",
        "tau_star_encloses_pi": contains(ts, mp.pi),
        "tau_star_width": float(ts.delta.b),
        "crossing_u_halfwidth": float(cs[0].delta.b) / 2,
        "crossing_pu_encloses_neg2": contains(cs[2], mp.mpf(-2)),
        "dsigma_dot_f_encloses_neg_half": contains(res["dsigma_dot_f"], mp.mpf("-0.5")),
        "DP_row_u_is_zero": all(contains(dp[0][j], mp.mpf(0)) for j in range(5)),
        "DP_row_u_max_halfwidth": max(float(dp[0][j].delta.b) for j in range(5)) / 2,
        "transversal": bool(res["transversal"]),
    }


# --------------------------------------------------------------------------- #
# Real W-Z Oterma section point.                                               #
# --------------------------------------------------------------------------- #
def _sig_2uv(iv: Any, st: list[Any]) -> Any:
    return iv.mpf(2) * st[0] * st[1]


def _grad_2uv(iv: Any, st: list[Any]) -> list[Any]:
    return [iv.mpf(2) * st[1], iv.mpf(2) * st[0], iv.mpf(0), iv.mpf(0), iv.mpf(0)]


def _oterma_ic(iv: Any) -> tuple[list[Any], float, float]:
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
    w0 = vti.lc_secondary_from_physical(
        iv, iv.mpf([x0, x0]), iv.mpf([0, 0]), iv.mpf([vx0, vx0]), vy0, MU
    )
    return w0, x0, vx0


def oterma_section_map(iv: Any) -> dict[str, Any]:
    w0, x0, vx0 = _oterma_ic(iv)
    jet = vti.make_cr3bp_lc_secondary_jet(H_ENERGY)
    vjet = vti.make_cr3bp_lc_secondary_variational_jet(H_ENERGY)
    res = vti.rigorous_section_map(
        iv,
        jet,
        vjet,
        w0,
        MU,
        sigma_val=_sig_2uv,
        sigma_grad=_grad_2uv,
        tau_max=20.0,
        n_steps=160,
        order=10,
        tau_min=1.0,
    )
    out: dict[str, Any] = {
        "start_x": x0,
        "start_vx": vx0,
        "section": "{y = 0}  (physical),  sigma = 2uv in Levi-Civita coords",
        "found": bool(res["found"]),
        "reason": res["reason"],
        "validated_to_tau": res["validated_to"],
    }
    if not res["found"]:
        return out
    ts = res["tau_star"]
    cs = res["crossing_state"]
    phys = vti.lc_secondary_to_physical(iv, cs, MU)
    dp = res["section_jacobian"]
    stm = res["stm_full"]
    out.update(
        {
            "crossing_step": res["crossing_step"],
            "tau_star_lo": float(ts.a),
            "tau_star_hi": float(ts.b),
            "tau_star_width": float(ts.delta.b),
            "physical_x_lo": float(phys[0].a),
            "physical_x_hi": float(phys[0].b),
            "physical_y_halfwidth": float(phys[1].delta.b) / 2,
            "physical_vx_lo": float(phys[2].a),
            "physical_vx_hi": float(phys[2].b),
            "physical_time_T_lo": float(cs[4].a),
            "physical_time_T_hi": float(cs[4].b),
            "dsigma_dot_f_lo": float(res["dsigma_dot_f"].a),
            "dsigma_dot_f_hi": float(res["dsigma_dot_f"].b),
            "transversal": bool(res["transversal"]),
            "stm_max_halfwidth": max(float(stm[i][j].delta.b) for i in range(5) for j in range(5))
            / 2,
            "DP_max_halfwidth": max(float(dp[i][j].delta.b) for i in range(5) for j in range(5))
            / 2,
            "note": (
                "first forward re-crossing of {y=0}; lands just past the first Jupiter "
                "perijove (T~0.462). Not expected to equal golden point 2 (rounded IC does "
                "not track the unstable connection)."
            ),
        }
    )
    return out


def build_certificate() -> dict[str, Any]:
    import mpmath as mp

    mp.mp.dps = 40
    mp.iv.dps = 40
    iv = mp.iv

    t0 = time.time()
    harm = harmonic_control(iv)
    print(f"[#672] harmonic control done ({time.time() - t0:.0f}s)", flush=True)
    kep = kepler_lc_control(iv)
    print(f"[#672] Kepler-LC control done ({time.time() - t0:.0f}s)", flush=True)
    oterma = oterma_section_map(iv)
    print(f"[#672] Oterma section map done ({time.time() - t0:.0f}s)", flush=True)

    return {
        "task": "#672",
        "stage": "Stage 6 -- rigorous Poincare section-crossing map + Jacobian",
        "system": {"mu": MU, "C": C_LEVEL, "h_energy": H_ENERGY},
        "section_definition": (
            "W-Z use Theta = {y = 0} parameterized by (x, xdot); in Levi-Civita "
            "regularized coords y = eta = 2uv, so sigma(w) = 2uv = 0. Crossing isolated in "
            "the regularized fictitious time tau (the flow is carried regularized to survive "
            "the Jupiter perijoves between crossings); physical (x, xdot) recovered by the "
            "inverse Levi-Civita map at the enclosed crossing."
        ),
        "method": (
            "Interval-Newton root of sigma(state(s))=0 on the crossing step's rigorous "
            "Taylor MODEL (poly coeffs from the tight step-start enclosure box + an "
            "order-(p+1) Lagrange remainder bounded over the whole-step a-priori enclosure); "
            "the same step certifies transversality (0 not in dsigma/ds over the step). "
            "Section-map Jacobian DP = Dphi_tau* - f (Dsigma Dphi_tau*)/(Dsigma f) with "
            "Dphi_tau* the QR-accumulated flow STM (Phi=A@B, kept tight by the Lohner-C1 "
            "reframing) composed with the crossing step's partial STM Taylor model."
        ),
        "harmonic_closed_form_control": harm,
        "kepler_lc_closed_form_control": kep,
        "oterma_real_section_map": oterma,
        "findings": {
            "harmonic_DP_matches_closed_form": bool(
                harm["tau_star_encloses_3pi_over_4"]
                and harm["crossing_y2_encloses_neg_sqrt2"]
                and harm["DP_encloses_closed_form"]
                and harm["DP_excludes_wrong"]
            ),
            "kepler_lc_control_matches_closed_form": bool(
                kep["tau_star_encloses_pi"]
                and kep["crossing_pu_encloses_neg2"]
                and kep["dsigma_dot_f_encloses_neg_half"]
                and kep["DP_row_u_is_zero"]
            ),
            "oterma_crossing_rigorously_isolated": bool(oterma.get("found")),
            "oterma_transversal": bool(oterma.get("transversal", False)),
            "oterma_DP_nonvacuous_tight": bool(
                oterma.get("found") and oterma.get("DP_max_halfwidth", 1.0) < 1e-6
            ),
        },
        "scope_note": (
            "This is the section-map machinery only (crossing isolation, crossing state, and "
            "rigorous section-map Jacobian), validated end-to-end against closed-form controls "
            "and demonstrated tight and transversal on a real W-Z Oterma section point. "
            "Constructing the actual h-sets and proving covering relations along the W-Z chain "
            "is the next stage, explicitly not this dispatch's scope."
        ),
    }


def main() -> None:
    cert = build_certificate()
    OUT.write_text(json.dumps(cert, indent=2) + "\n")
    h = cert["harmonic_closed_form_control"]
    print(
        f"[#672] harmonic control: tau*=3pi/4 enclosed={h['tau_star_encloses_3pi_over_4']} "
        f"DP matches cf={h['DP_encloses_closed_form']} excludes wrong={h['DP_excludes_wrong']} "
        f"(DP hw={h['DP_max_halfwidth']:.1e})"
    )
    k = cert["kepler_lc_closed_form_control"]
    print(
        f"[#672] Kepler-LC control: tau*=pi enclosed={k['tau_star_encloses_pi']} "
        f"pu*=-2 enclosed={k['crossing_pu_encloses_neg2']} DP row-u zero={k['DP_row_u_is_zero']}"
    )
    o = cert["oterma_real_section_map"]
    if o["found"]:
        print(
            f"[#672] Oterma section map: tau*~[{o['tau_star_lo']:.6f},{o['tau_star_hi']:.6f}] "
            f"x*~{o['physical_x_lo']:.6f} vx*~{o['physical_vx_lo']:.6f} "
            f"T~{o['physical_time_T_hi']:.4f} "
            f"transversal={o['transversal']} DP hw={o['DP_max_halfwidth']:.1e}"
        )
    else:
        print(f"[#672] Oterma section map: NOT found -- {o['reason']}")
    print(f"[#672] certificate written -> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
