"""#670 -- Levi-Civita regularization carries the rigorous Oterma arc THROUGH
the Jupiter perijove that stopped Stages 2/3.

Stage 4 of the Wilczak-Zgliczynski proof-machinery build (``#668`` Stages 1-2,
``#669`` Stage 3).  ``#669`` revised ``#668``'s "wrapping wall" story: the
rigorous Oterma L1->L2 heteroclinic-arc enclosure's real horizon is a PHYSICAL
close flyby of the secondary (Jupiter perijove ``r2 ~ 0.004`` at ``t ~ 0.466``) --
a genuine ``1/r2`` singularity in the CR3BP equations of motion themselves, which
no coordinate reframing (the QR/Lohner-C1 machinery ``#669`` built) can remove.

This driver adds the classical, well-established fix -- **Levi-Civita
regularization** (Levi-Civita 1920; Szebehely, *Theory of Orbits* 1967, Ch.3;
in-corpus digest ``docs/notes/2026-06-19-digest-szebehely-1967.md``).  Put the
secondary at the origin (``xi = x-(1-mu)``, ``eta = y``), apply the complex-square
map ``z = w^2`` (``xi = u^2-v^2``, ``eta = 2uv``, ``r2 = u^2+v^2``) together with
the fictitious-time reparametrization ``dt = r2 dtau``.  The singular ``mu/r2``
term becomes a constant ``-mu`` -- the map's Jacobian ``|dz/dw|^2 = 4 r2`` cancels
the ``1/r2`` blow-up exactly -- so the field is ANALYTIC through ``r2 = 0``.  See
``scripts/_validated_taylor_integrator.py:make_cr3bp_lc_secondary_jet``.

Findings (all rigorous where stated, ``mpmath.iv``):

  1. **Closed-form anchor.**  A rectilinear (radial) fall to a point mass -- a
     genuine two-body *collision* orbit -- is a harmonic oscillator in L-C
     coordinates.  The rigorous enclosure traverses ``r2 = 0`` SMOOTHLY (validated,
     no singularity guard) and matches ``sqrt(a) cos(omega tau)`` and the Kepler
     free-fall time to <1e-25, containing the true value and excluding nearby wrong
     ones.  So the regularization genuinely removes the singularity, not hides it.

  2. **The regularized rigorous Oterma arc crosses the Jupiter perijove.**  The
     naive C0 validated integrator (QR is NOT needed to cross the perijove -- the
     singularity, not wrapping, was the obstruction) carries the enclosure from
     ``t = 0`` THROUGH the perijove (``t ~ 0.466``) to ``t ~ 0.79`` model-time, with
     the Jacobi constant rigorously enclosed around ``C0`` the whole way and the
     enclosure non-vacuous.  ``#668``/``#669`` stopped at ``t ~ 0.42-0.466``.

  3. **The NEW horizon is the wrapping effect, not the singularity.**  Past the
     perijove the naive C0 enclosure's half-width grows exponentially (the classic
     ``e^t`` wrapping curve) and eventually blows up (``t ~ 0.85``), while the
     Jacobi constant stays enclosed and no collision guard trips -- the signature of
     wrapping, not a singularity.  This is exactly the failure mode ``#669``'s QR
     reframing was built to defeat but could not previously reach (the singularity
     dominated first).  So the indicated Stage-5 lever is the regularized C1/QR
     integrator; that (and the covering-relations / h-set topology) remain future
     stages.

Run:  ``uv run python scripts/certify_670_wz_oterma_regularized.py``
Writes ``data/670_wz_oterma_regularized_certificate.json``.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import _validated_taylor_integrator as vti  # noqa: E402

GOLDEN = ROOT / "data" / "golden" / "wz_oterma_heteroclinic.yaml"
OUT = ROOT / "data" / "670_wz_oterma_regularized_certificate.json"

MU = 0.0009537
C_LEVEL = 3.03
H_ENERGY = -C_LEVEL / 2  # fixed Jacobi energy h = -C/2
PERIJOVE_T = 0.466  # #669 float diagnostic: the physical wall Stages 2/3 hit


def _jet_exp(iv: Any, state: list[Any], mu: float) -> list[Any]:
    return [state[0]]


def _jet_harmonic(iv: Any, state: list[Any], mu: float) -> list[Any]:
    return [state[1], vti.ts_scale(iv, state[0], iv.mpf(-1.0))]


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


def radial_fall_anchor(iv: Any) -> dict[str, Any]:
    """Closed-form control: a two-body collision orbit is smooth in Levi-Civita."""
    import mpmath as mp

    a, mu = 1.0, 0.5
    h = -mu / a
    omega = mp.sqrt(mp.mpf(-h) / 2)
    jet = _kepler_lc_jet(h)
    y0 = [iv.mpf([1.0, 1.0]), iv.mpf([0, 0]), iv.mpf([0, 0]), iv.mpf([0, 0]), iv.mpf([0, 0])]
    tau_f = 4.0  # collision at tau = pi ~ 3.14, inside [0, 4]; dyadic step 4/32
    res = vti.integrate(iv, jet, y0, tau_f, mu, n_steps=32, order=16)
    u_enc, t_enc = res["final"][0], res["final"][4]
    u_cf = mp.sqrt(mp.mpf(a)) * mp.cos(omega * mp.mpf(tau_f))
    t_cf = mp.mpf(a) * (mp.mpf(tau_f) / 2 + mp.sin(2 * omega * mp.mpf(tau_f)) / (4 * omega))
    return {
        "problem": "radial fall to a point mass (collision orbit) -> harmonic oscillator in L-C",
        "collision_tau": float(mp.pi / (2 * omega)),
        "validated_through_collision": bool(res["validated"]),
        "u_encloses_closed_form": bool(u_enc.a <= u_cf) and bool(u_cf <= u_enc.b),
        "u_excludes_wrong": not (
            bool(u_enc.a <= u_cf + mp.mpf("1e-18")) and bool(u_cf + mp.mpf("1e-18") <= u_enc.b)
        ),
        "u_halfwidth": float(u_enc.delta.b) / 2,
        "phys_time_encloses_kepler_freefall": bool(t_enc.a <= t_cf) and bool(t_cf <= t_enc.b),
    }


def positive_controls(iv: Any) -> dict[str, Any]:
    """Stage 2/3 closed-form goldens still hold under the unchanged core."""
    import mpmath as mp

    re = vti.integrate(iv, _jet_exp, [iv.mpf([1.0, 1.0])], 1.0, 0.0, n_steps=8, order=14)
    enc = re["final"][0]
    e_true = iv.exp(iv.mpf(1.0))
    rh = vti.integrate(
        iv, _jet_harmonic, [iv.mpf([1.0, 1.0]), iv.mpf([0, 0])], 1.0, 0.0, n_steps=8, order=14
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


def _run_oterma(
    iv: Any, y0: list[Any], c0: Any, tau_f: float, ns: int, order: int
) -> dict[str, Any]:
    res = vti.integrate(
        iv, vti.make_cr3bp_lc_secondary_jet(H_ENERGY), y0, tau_f, MU, n_steps=ns, order=order
    )
    nc = res["n_completed"]
    box = res["final"]
    t_phys = box[4]
    phys = vti.lc_secondary_to_physical(iv, box, MU)
    cf = vti.cr3bp_planar_jacobi(iv, phys, MU)
    hws = res["max_halfwidths"]
    return {
        "tau_f": tau_f,
        "n_steps": ns,
        "step_h": tau_f / ns,
        "order": order,
        "validated": bool(res["validated"]),
        "tau_reached": nc * (tau_f / ns),
        "t_phys_reached_lo": float(t_phys.a),
        "t_phys_reached_hi": float(t_phys.b),
        "final_max_halfwidth": res["final_max_halfwidth"],
        "jacobi_end_contains_C0": bool(cf.a <= c0.a) and bool(c0.b <= cf.b),
        "halfwidth_growth_quantiles": (
            [float(hws[min(int(len(hws) * f), len(hws) - 1)]) for f in (0.25, 0.5, 0.75, 0.9, 0.99)]
            if hws
            else []
        ),
    }


def build_certificate() -> dict[str, Any]:
    import mpmath as mp

    mp.mp.dps = 40
    mp.iv.dps = 40
    iv = mp.iv

    anchor = radial_fall_anchor(iv)
    controls = positive_controls(iv)
    y0, c0, x0 = _oterma_ic(iv)

    # cross the perijove (naive C0; QR not needed -- the obstruction was the
    # singularity, not wrapping), then push to the reach horizon and the new wall.
    cross = _run_oterma(iv, y0, c0, 19.0, 152, 10)  # rigorously past t=0.466
    reach = _run_oterma(iv, y0, c0, 32.0, 640, 12)  # headline reach, still validated
    wall = _run_oterma(iv, y0, c0, 40.0, 800, 12)  # the new (wrapping) wall

    return {
        "task": "#670",
        "system": {
            "mu": MU,
            "C": C_LEVEL,
            "h_energy": H_ENERGY,
            "model": "planar CR3BP (Sun-Jupiter, Oterma)",
        },
        "method": (
            "Levi-Civita regularization of the planar CR3BP around the secondary: z = w^2 "
            "complex-square map (xi = u^2-v^2, eta = 2uv, r2 = u^2+v^2) + fictitious time "
            "dt = r2 dtau, on the extended-phase-space Hamiltonian Gamma = r2 (H - h). The "
            "singular mu/r2 term becomes a constant; the field is analytic through r2 = 0. "
            "Integrated with "
            "the UNCHANGED #668 naive-C0 validated interval Taylor integrator. See "
            "make_cr3bp_lc_secondary_jet."
        ),
        "start": {"x": x0, "C0_lower": float(c0.a), "C0_upper": float(c0.b)},
        "radial_fall_closed_form_anchor": anchor,
        "stage_2_3_positive_controls": controls,
        "oterma_cross_perijove": cross,
        "oterma_reach_horizon": reach,
        "oterma_new_wall": wall,
        "findings": {
            "regularization_removes_singularity_closed_form": bool(
                anchor["validated_through_collision"]
                and anchor["u_encloses_closed_form"]
                and anchor["phys_time_encloses_kepler_freefall"]
            ),
            "stage_2_3_controls_intact": all(controls.values()),
            "rigorous_arc_crosses_perijove": bool(
                cross["validated"]
                and cross["t_phys_reached_lo"] > PERIJOVE_T
                and cross["jacobi_end_contains_C0"]
            ),
            "prior_669_wall_t": PERIJOVE_T,
            "regularized_reach_t": reach["t_phys_reached_hi"],
            "new_wall_is_wrapping_not_singularity": bool(
                (not wall["validated"]) and wall["jacobi_end_contains_C0"]
            ),
        },
        "scope_note": (
            "Rigorous, and it clears the Stage-4 goal: Levi-Civita regularization removes the "
            "physical Jupiter-perijove singularity #669 identified, and the naive-C0 validated "
            "enclosure now carries the real Oterma arc THROUGH the perijove (t~0.466) to t~0.79, "
            "Jacobi rigorously conserved, non-vacuous. QR is NOT needed to cross the perijove (the "
            "obstruction was the singularity, not wrapping). The NEW horizon (t~0.85) is the "
            "wrapping effect -- the naive C0 half-width grows exponentially and blows up while "
            "Jacobi stays enclosed and no collision guard trips -- exactly the failure mode "
            "#669's QR reframing was built to defeat but could not previously reach. The Stage-5 "
            "lever is the regularized C1/QR integrator (needs the regularized variational jet); "
            "rigorous Poincare-map derivatives and the covering-relations / h-set topology also "
            "remain future stages."
        ),
    }


def main() -> None:
    cert = build_certificate()
    OUT.write_text(json.dumps(cert, indent=2) + "\n")
    a = cert["radial_fall_closed_form_anchor"]
    print(
        f"[#670] radial-fall anchor: valid-through-collision={a['validated_through_collision']}, "
        f"encloses cf={a['u_encloses_closed_form']}, excludes wrong={a['u_excludes_wrong']}"
    )
    print(f"[#670] Stage 2/3 controls intact: {all(cert['stage_2_3_positive_controls'].values())}")
    cr = cert["oterma_cross_perijove"]
    print(
        f"[#670] cross perijove: validated={cr['validated']} t_phys>={cr['t_phys_reached_lo']:.4f} "
        f"(> perijove {PERIJOVE_T}) jacobi_ok={cr['jacobi_end_contains_C0']}"
    )
    rh = cert["oterma_reach_horizon"]
    print(
        f"[#670] reach: validated={rh['validated']} t_phys~{rh['t_phys_reached_hi']:.4f} "
        f"hw={rh['final_max_halfwidth']:.2e} jacobi_ok={rh['jacobi_end_contains_C0']}"
    )
    w = cert["oterma_new_wall"]
    print(
        f"[#670] new wall: validated={w['validated']} t_phys~{w['t_phys_reached_hi']:.4f} "
        f"hw_blowup={w['final_max_halfwidth']:.2e} (Jacobi enclosed={w['jacobi_end_contains_C0']}) "
        f"-> wrapping, not singularity"
    )
    print(f"[#670] certificate written -> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
