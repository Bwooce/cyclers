"""#792 scoping -- closed-form structure of #680's asymmetric-closure continuum.

`#792` proposed a ~1-2 week "adaptive basin-width-aware asymmetric-closure
grid" census at Uranus. `#680` (2026-07-22, empty-regions region
``uranus-asymmetric-closure-freebeta-degenerate-manifold-2026-07-22``) had
already shown empirically that the free-``rel_offset`` (beta) V-infinity-
MAGNITUDE closure solution set is a degenerate continuous >=1-D manifold, not
a discrete set of isolated closures. This script is the reproducible artifact
for the `#792` scoping verdict (see
``docs/notes/2026-08-10-792-scoping-vs-680.md``): that manifold has an exact
CLOSED FORM, so there is nothing for any grid -- adaptive or otherwise -- to
census, and the one missing periodicity EQUALITY collapses the manifold onto
the already-catalogued symmetric (`#563`/`#569`) family.

Checks (all against `scan_558`'s own ``residual_at_point`` and the repo's own
``lambert()``, no new physics):

(A) Planar Tisserand linearity: at a circular moon of radius r, a leg with
    two-body invariants (E, h) has
        |vinf|^2 = 2E + 3*mu/r - 2*sqrt(mu)*h*r**-1.5,
    LINEAR in (E, h). Verified to ~1e-14 on real Lambert legs. Because the
    genome's closure residual matches |vinf| at BOTH radii r_A != r_B, and
    the 2x2 linear map (E, h) -> (|vinf_A|^2, |vinf_B|^2) is invertible
    (det = 4*sqrt(mu)*(r_A**-1.5 - r_B**-1.5) != 0), magnitude closure forces
    E1 == E0 AND h1 == h0 exactly: the return leg is a congruent (rotated)
    copy of the outbound conic.

(B) Both legs are Lambert arcs between the same two radii with the SAME
    duration (single ``tof`` in the genome), so congruence forces the mirror
    arc: equal transfer angles. With dnu0 = beta + n_b*tof and
    dnu1 = 2*n_a*tof - beta - n_b*tof this is ONE scalar equation
        beta == (n_a - n_b) * tof   (mod 180 deg)
    in TWO unknowns -> a 1-D line family = `#680`'s continuum, in closed
    form. Verified: at PREDICTED points, ``residual_at_point`` returns
    ~1e-13..1e-15 km/s with no refinement, in BOTH the (0,0) and (2,2) rev
    classes (the latter is the class `#680`'s deflated Newton sampled), and
    a local refine moves tof by <1e-14 d. T_syn/2 = 180/(n_a - n_b) equals
    the symmetric `#563` commensurate step exactly -- the symmetric
    enumeration is precisely the beta == 0 (mod 180) slice of this form.

(C) Along the manifold, the REQUIRED turn angle at the flyby moon and the
    moon-local vinf direction mismatch at the anchor vary SMOOTHLY: flyby
    turn feasibility is an inequality (bend is a free periapsis-choice
    parameter), so a stricter turn-enforcing formulation selects sub-arcs of
    the same curve -- it cannot isolate points.

(D) The one missing EQUALITY for true periodicity is pattern repeat:
    (n_a - n_b)*2*tof == 0 (mod 360). On the manifold this reads
    2*beta == 0 (mod 360), i.e. beta in {0, 180} -- exactly the symmetric,
    already-catalogued family. Printed as the per-cycle pattern drift
    (= 2*beta mod 360, up to sign), nonzero everywhere off the symmetric
    slice.

Discipline: no catalogue write, no registry write, no new search method.

Run as::

    uv run python scripts/check_792_manifold_closed_form.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from scan_558_uranus_all_pairs_offset_sweep import residual_at_point  # noqa: E402

from cyclerfinder.core.lambert import lambert  # noqa: E402
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES  # noqa: E402
from cyclerfinder.search.discovery_campaign import (  # noqa: E402
    DAY_S,
    _mean_motion_rad_day,
    _moon_state,
)

MU = PRIMARIES["Uranus"]
ANCHOR, FLYBY = "Ariel", "Umbriel"
SAT_A, SAT_B = SATELLITES[ANCHOR], SATELLITES[FLYBY]
N_A = _mean_motion_rad_day(MU, SAT_A.sma_km)  # rad/day
N_B = _mean_motion_rad_day(MU, SAT_B.sma_km)
N_A_DEG, N_B_DEG = math.degrees(N_A), math.degrees(N_B)
DN_DEG = N_A_DEG - N_B_DEG
TOF_SCALE_UNIT = math.sqrt((2 * math.pi / N_A) * (2 * math.pi / N_B))


def residual_norm(beta: float, tof: float, n_rev: tuple[int, int]) -> float:
    pt = residual_at_point(
        ANCHOR, FLYBY, rel_offset_deg=beta, tof_scale=tof / TOF_SCALE_UNIT, n_rev=n_rev
    )
    return math.inf if pt is None else float(pt["residual_kms"])


def refine_tof(beta: float, tof0: float, n_rev: tuple[int, int]) -> tuple[float, float]:
    """Golden-section local refine of tof around tof0; returns (tof, residual)."""
    lo, hi = tof0 - 0.05, tof0 + 0.05
    for _ in range(60):
        m1 = lo + 0.382 * (hi - lo)
        m2 = lo + 0.618 * (hi - lo)
        if residual_norm(beta, m1, n_rev) < residual_norm(beta, m2, n_rev):
            hi = m2
        else:
            lo = m1
    t = 0.5 * (lo + hi)
    return t, residual_norm(beta, t, n_rev)


def leg_vectors(beta: float, tof: float) -> dict[str, object]:
    """vinf VECTORS at the 3 encounters + per-leg (E, h), direct (0-rev) legs."""
    r0, v0 = _moon_state(0.0, N_A, 0.0, SAT_A.sma_km, MU)
    r1, v1 = _moon_state(math.radians(beta), N_B, tof, SAT_B.sma_km, MU)
    r2, v2 = _moon_state(0.0, N_A, 2.0 * tof, SAT_A.sma_km, MU)
    sol0 = lambert(r0, r1, tof * DAY_S, mu=MU, max_revs=0)[0]
    sol1 = lambert(r1, r2, tof * DAY_S, mu=MU, max_revs=0)[0]
    eh = []
    for r, v in ((r0, np.asarray(sol0.v1)), (r1, np.asarray(sol1.v1))):
        energy = 0.5 * float(np.dot(v, v)) - MU / float(np.linalg.norm(r))
        eh.append((energy, float(np.cross(r, v)[2])))
    return {
        "vinf_out_a0": np.asarray(sol0.v1) - v0,
        "vinf_in_b": np.asarray(sol0.v2) - v1,
        "vinf_out_b": np.asarray(sol1.v1) - v1,
        "vinf_in_a2": np.asarray(sol1.v2) - v2,
        "eh": eh,
    }


def angle_deg(u: np.ndarray, w: np.ndarray) -> float:
    c = float(np.dot(u, w) / (np.linalg.norm(u) * np.linalg.norm(w)))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def check_a_tisserand_linearity() -> None:
    print("(A) Tisserand linearity |vinf|^2 = 2E + 3mu/r - 2 sqrt(mu) h r^-1.5:")
    lv = leg_vectors(311.0, 2.9)
    eh = lv["eh"]
    assert isinstance(eh, list)
    for name, vec, (energy, h), r in (
        ("leg0@anchor", lv["vinf_out_a0"], eh[0], SAT_A.sma_km),
        ("leg0@flyby ", lv["vinf_in_b"], eh[0], SAT_B.sma_km),
        ("leg1@flyby ", lv["vinf_out_b"], eh[1], SAT_B.sma_km),
        ("leg1@anchor", lv["vinf_in_a2"], eh[1], SAT_A.sma_km),
    ):
        assert isinstance(vec, np.ndarray)
        lhs = float(np.dot(vec, vec))
        rhs = 2 * energy + 3 * MU / r - 2 * math.sqrt(MU) * h / r**1.5
        print(f"    {name}: |vinf|^2={lhs:.9f}  formula={rhs:.9f}  |diff|={abs(lhs - rhs):.2e}")


def check_b_manifold(n_rev: tuple[int, int], betas_k: list[tuple[float, int]]) -> None:
    print(f"(B) closed form beta == (n_a-n_b)*tof (mod 180), n_rev={n_rev} (per-point branch k):")
    print(f"    n_a-n_b = {DN_DEG:.6f} deg/day; predicted tof = (beta + 180*k)/(n_a-n_b)")
    for beta, k in betas_k:
        tof_pred = (beta + 180.0 * k) / DN_DEG
        res_pred = residual_norm(beta, tof_pred, n_rev)
        tof_root, res_root = refine_tof(beta, tof_pred, n_rev)
        drift = (2.0 * beta) % 360.0
        print(
            f"    beta={beta:7.2f}: tof_pred={tof_pred:8.4f} res(pred)={res_pred:.2e} | "
            f"refined tof={tof_root:8.4f} res={res_root:.2e} "
            f"|dtof|={abs(tof_root - tof_pred):.1e} | pattern drift={drift:6.1f} deg"
        )


def check_c_turn_angles(betas: list[float], k: int) -> None:
    print("(C) required turn / direction mismatch ALONG the manifold (smooth => no isolation):")
    for beta in betas:
        tof_root, res = refine_tof(beta, (beta + 180.0 * k) / DN_DEG, (0, 0))
        lv = leg_vectors(beta, tof_root)
        vin_b, vout_b = lv["vinf_in_b"], lv["vinf_out_b"]
        vin_a2, vout_a0 = lv["vinf_in_a2"], lv["vinf_out_a0"]
        assert isinstance(vin_b, np.ndarray) and isinstance(vout_b, np.ndarray)
        assert isinstance(vin_a2, np.ndarray) and isinstance(vout_a0, np.ndarray)
        turn_b = angle_deg(vin_b, vout_b)
        theta = N_A * 2.0 * tof_root
        c, s = math.cos(-theta), math.sin(-theta)
        rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
        turn_a = angle_deg(rot @ vin_a2, vout_a0)
        print(
            f"    beta={beta:7.2f} tof={tof_root:7.4f} res={res:.1e}: "
            f"required turn@flyby={turn_b:7.3f} deg, moon-local mismatch@anchor={turn_a:7.3f} deg"
        )


def main() -> int:
    print("#792 scoping -- closed form of #680's asymmetric V-inf-magnitude closure continuum")
    print(f"pair {ANCHOR}-{FLYBY}: n_a={N_A_DEG:.4f} deg/d, n_b={N_B_DEG:.4f} deg/d")
    print(
        f"T_syn/2 = 180/(n_a-n_b) = {180.0 / DN_DEG:.6f} d "
        "(== the #563 symmetric commensurate step and the Ariel-Umbriel #569 golden tof)"
    )
    check_a_tisserand_linearity()
    # #680 finding-4 window (beta 320..360 <-> tof 2.50..3.216 d) is branch k=-1:
    check_b_manifold((0, 0), [(b, -1) for b in (320.0, 327.0, 335.0, 344.0, 351.0, 359.0)])
    check_b_manifold((0, 0), [(b, 1) for b in (90.0, 120.0, 150.0, 179.0)])
    # the rev class #680's deflated-Newton continuum (finding 1) actually sampled; the
    # branch k carrying the (2,2) Lambert arc differs per beta (branch existence, not
    # closed-form failure -- the line family is rev-class-independent, its Lambert
    # carrier class varies along each line):
    check_b_manifold((2, 2), [(200.0, 1), (233.7, 2), (267.13, 1), (301.0, 1)])
    check_c_turn_angles([305.0, 320.0, 335.0, 344.0, 351.0, 359.0], k=-1)
    print("(D) [analytic] per-cycle A-B pattern drift on the manifold = 2*beta (mod 360);")
    print("    drift == 0  <=>  beta == 0 (mod 180)  <=>  the symmetric #563/#569 family.")
    print("    => true-periodicity (the missing EQUALITY) collapses the continuum onto the")
    print("       already-catalogued symmetric closures; turn-angle checks are INEQUALITIES")
    print("       and only select sub-arcs. Nothing for #792's adaptive grid to census.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
