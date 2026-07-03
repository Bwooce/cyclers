"""#523 Earth co-orbital (Arjuna-class) horseshoe/quasi-satellite periodic-orbit search.

Points the #450 DA/HOTM global fixed-point enumerator (already reused once for
#527) at the Sun-Earth CR3BP co-orbital region (C ~ 3, near L3/L4/L5), asking
whether a CERTIFIED periodic horseshoe/quasi-satellite orbit exists that makes
repeated close (sub-Hill-sphere) approaches to Earth -- the co-orbital-exchange
"cycler" mechanism #523 proposed, which needs no flyby bend.

ENCOUNTER CRITERION (settled in writing before this sweep runs, per the
review's own #339-style-criterion-trap warning): a trajectory counts as
making a repeated ENCOUNTER with Earth if it passes within Earth's Hill-sphere
radius (R_hill = a_Earth * (mu/3)^(1/3), the standard measure of a body's
sphere of gravitational relevance -- the same measure used for the #527 Hill
Jupiter analysis this session) during each cycle of a CERTIFIED periodic
orbit. This is deliberately a "gravitationally relevant close approach"
criterion, not a flyby-distance criterion -- #523's whole premise is that
these encounters do NOT need flyby-grade bend, unlike every prior small-body
search this project has run.

POSITIVE CONTROL: seeded from a REAL, published Earth co-orbital object's
Keplerian elements (2006 RH120, de la Fuente Marcos & de la Fuente Marcos
2018, Astron. Nachr., Table 1: a=0.998625 AU, e=0.019833, i=1.52613 deg --
notable because 2006 RH120 was an actual, observationally confirmed transient
Earth minimoon, Jul 2006-Jul 2007). Direct integration of this seed in the
Sun-Earth CR3BP independently reproduces the REAL qualitative behavior: an
initial ~0.7-year close quasi-satellite episode reaching HALF the Hill
radius.

CORRECTION (2026-07-03, #535 investigation): an earlier interactive check
claimed this seed ALSO transitions into a wider horseshoe libration with
recurring approaches every ~4.6-9.2 years. That claim does NOT hold up under
careful re-verification (long-duration integration + an independent
vis-viva re-derivation of the same IC from the raw Keplerian elements,
which reproduces this seed's (x0,ydot0,C) to 4 significant figures,
confirming the seed itself is NOT the bug): the trajectory makes exactly
ONE close approach then departs to ~2 AU from Earth and does not return
within at least 60 years, consistent with the REAL 2006 RH120 being a
genuinely single, ~1-year transient capture (not a multi-decade recurring
horseshoe companion). Only the single-encounter part of the original claim
is confirmed; the "recurring approaches" part was an error, never actually
re-checked before being written down. See #535's OUTSTANDING.md entry and
`docs/notes/2026-07-03-535-quasi-cycler-transient-drift-admission-
criterion.md` for the full investigation.

SYSTEM-SCALE NOTE (learned from #527, re-verified here rather than assumed):
Sun-Earth's characteristic time t_s ~= 58.13 days (2*pi*t_s ~= 1 year, as
expected). Direct integration of the positive-control seed showed section
returns (first y=0 crossing with matching ydot sign) ranging from ~4.4 to
~58 nondimensional time units depending on where in the horseshoe cycle the
seed starts -- much more variable than #527's clean single-resonance case.
This script uses a generous t_max=100 for the coarse enumerator and an
iterative period-guess refinement loop (up to 5 passes, chaining each pass's
own period estimate into the next) for certification, since a single-shot
correct_general_periodic call from a coarse candidate did not always
converge to tol on the first attempt in interactive testing.
"""

from __future__ import annotations

import datetime
import pathlib
import sys
import time

# Ensure the src tree is on the path when invoked as a script.
_REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "src"))

import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402

import cyclerfinder.core.cr3bp as cr3bp  # noqa: E402
from cyclerfinder.data.method_capability import MethodCapability  # noqa: E402
from cyclerfinder.data.preflight import PreflightBlockedError, preflight_search  # noqa: E402
from cyclerfinder.genome.da_hotm_backend import SamplingSectionMap  # noqa: E402
from cyclerfinder.genome.da_hotm_enumerator import DomainBox, enumerate_fixed_points  # noqa: E402
from cyclerfinder.search.cr3bp_general_periodic import (  # noqa: E402
    GeneralPeriodicOrbit,
    correct_general_periodic,
)

# 2006 RH120 sourced elements (de la Fuente Marcos & de la Fuente Marcos 2018,
# Astron. Nachr., Table 1, JPL epoch 2456800.5): a=0.998625 AU, e=0.019833,
# i=1.52613 deg. Converted (this session, interactively verified) to a
# Sun-Earth CR3BP rotating-frame seed via a planar circular-orbit-style
# coe_to_rv + frame-rotation-subtraction construction (see #527's script for
# the same method applied to Sun-Jupiter/Hilda).
SEED_X = 0.97881371
SEED_XDOT = 0.0
SEED_C = 2.9998797409719242

C_BAND = (2.9990, 2.9995, SEED_C, 3.0000, 3.0005, 3.0010)
N_RANGE = (1,)
DOMAIN_BOX = DomainBox(x_lo=0.80, x_hi=1.10, xdot_lo=-0.08, xdot_hi=0.08)
GRID = (21, 15)
# Tightened from an initial 5e-2 (2026-07-02, this session): the looser value let
# too many near-duplicate coarse candidates through, each requiring a full 5-pass
# certification -- a first attempt at this tolerance took 16+ min to cover even
# one Jacobi value (with an already-certified orbit re-certified byte-identically
# from a second nearby seed) and was killed. 1e-2 matches #527's tolerance.
RESIDUAL_TOL = 1e-2
# Post-certification dedup radius: TWO coarse seeds separated by MORE than the
# enumerator's own dedup_radius (0.01, in the coarse residual landscape) can still
# converge under Newton refinement to the SAME true fixed point (observed
# directly this session -- x0=0.83658 certified twice, byte-identically, from a
# 5e-2-tolerance run). This is a second, independent dedup on the CERTIFIED
# (x0, xdot0) output, not a substitute for the coarse one.
CERTIFIED_DEDUP_RADIUS = 1e-3
T_MAX = 100.0  # nondim; generous vs the ~4-58 range observed interactively.

_REGION_ID = "sun-earth-coorbital-horseshoe-qsat-dahotm"
_METHOD = MethodCapability(
    genome="DA/HOTM global Poincare-section fixed-point enumeration (#450)",
    corrector="correct_general_periodic (asymmetric single-shooting, analytic STM)",
    capability_tags=frozenset(
        {"cr3bp", "ballistic", "coplanar", "single-arc", "poincare-section-enumeration"}
    ),
    git_sha="working-tree",
)


def _ts() -> str:
    return datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")


def _min_dist_to_earth(
    system: cr3bp.CR3BPSystem, x0: float, xdot0: float, ydot0: float, period: float
) -> float:
    mu = float(system.mu)
    state0 = np.array([x0, 0.0, 0.0, xdot0, ydot0, 0.0], dtype=np.float64)
    t_eval = np.linspace(0.0, abs(period), 600)
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, abs(period)),
        state0,
        args=(mu,),
        method="DOP853",
        rtol=1e-11,
        atol=1e-11,
        t_eval=t_eval,
    )
    if not sol.success:
        return float("inf")
    dx = sol.y[0] - (1.0 - mu)
    dy = sol.y[1]
    return float(np.min(np.hypot(dx, dy)))


def _certify_with_refinement(
    system: cr3bp.CR3BPSystem, x0: float, xdot0: float, c_target: float, n: int
) -> GeneralPeriodicOrbit:
    """Chain up to 3 correct_general_periodic passes, each seeded from the last.

    A single-shot call from a coarse enumerator candidate did not always
    converge on the first attempt (interactive testing this session); this
    mirrors the manual refinement that DID converge -- but a genuine
    convergent seed converges in well under 30 Newton iterations (observed
    interactively), so max_iter is capped there rather than 100: a candidate
    that WON'T converge should fail fast, not spend 5 x 100 long-horizon
    iterations finding that out (an earlier version of this script did
    exactly that -- one non-converging candidate took 2+ minutes). Bails out
    after 2 consecutive passes with no residual improvement (a stalled seed,
    not making progress) instead of always spending all 3 passes.
    """
    period_guess = T_MAX
    orbit = correct_general_periodic(
        system,
        x0,
        xdot0,
        c_target,
        period_guess=period_guess,
        half_crossings=2 * n,
        ydot0_sign=1.0,
        tol=1e-11,
        t_hi_frac=1.15,
        max_iter=30,
    )
    prev_residual = float("inf")
    for _ in range(2):
        if orbit.converged and orbit.residual <= 1e-9:
            return orbit
        if orbit.residual >= prev_residual:
            break  # stalled: not improving, don't spend another pass on it
        prev_residual = orbit.residual
        x0, xdot0, period_guess = orbit.x0, orbit.xdot0, orbit.period
        orbit = correct_general_periodic(
            system,
            x0,
            xdot0,
            c_target,
            period_guess=period_guess,
            half_crossings=2 * n,
            ydot0_sign=1.0,
            tol=1e-11,
            t_hi_frac=1.15,
            max_iter=30,
        )
    return orbit


def main() -> None:
    print(f"[{_ts()}] #523 Sun-Earth co-orbital horseshoe/QS DA/HOTM search starting.")

    system = cr3bp.cr3bp_system("Sun", "Earth")
    r_hill = (float(system.mu) / 3.0) ** (1.0 / 3.0)
    n_points = len(C_BAND) * len(N_RANGE) * GRID[0] * GRID[1]

    preflight_search(
        task_no=523,
        region_id=_REGION_ID,
        method=_METHOD,
        script_path=pathlib.Path(__file__),
        n_points=n_points,
        # Measured this session via direct interactive timing on this exact
        # (system, box, grid) configuration: ~0.045 s/point for the coarse
        # enumerator pass (certification cost is separate, per-candidate).
        timing_pilot_seconds_per_point=0.05,
    )

    print(
        f"[{_ts()}] Positive control: seed from 2006 RH120's real published elements "
        f"(C={SEED_C:.6f}, x~={SEED_X}, xdot~={SEED_XDOT}) ..."
    )
    control_backend = SamplingSectionMap(system, c_target=SEED_C, ydot_sign=1.0, t_max=T_MAX)
    control_box = DomainBox(x_lo=0.85, x_hi=1.05, xdot_lo=-0.06, xdot_hi=0.06)
    control_candidates = enumerate_fixed_points(
        control_backend, control_box, 1, residual_tol=RESIDUAL_TOL, grid=(6, 6)
    )
    if not control_candidates:
        raise RuntimeError(
            "POSITIVE CONTROL FAILED: no sub-tolerance candidate recovered near the "
            "2006-RH120-derived seed. Do not trust any negative from the full scan below."
        )
    print(
        f"[{_ts()}] Positive control PASSED: {len(control_candidates)} candidate(s), "
        f"best residual={min(c.residual for c in control_candidates):.3e}"
    )

    print(
        f"[{_ts()}] Full scan: {len(C_BAND)} Jacobi values x {len(N_RANGE)} rev counts "
        f"x {GRID[0] * GRID[1]} grid points = {n_points} total."
    )

    t0 = time.time()
    n_certified = 0
    n_encounter = 0
    n_dup_skipped = 0
    certified_seen: list[tuple[float, float, float]] = []  # (c_target, x0, xdot0)
    findings: list[dict[str, float | int | bool]] = []
    for c_idx, c_target in enumerate(C_BAND):
        backend = SamplingSectionMap(system, c_target=float(c_target), ydot_sign=1.0, t_max=T_MAX)
        for n in N_RANGE:
            cands = enumerate_fixed_points(
                backend, DOMAIN_BOX, n, residual_tol=RESIDUAL_TOL, grid=GRID, dedup_radius=0.01
            )
            print(
                f"[{_ts()}] C={c_target:.4f} n={n} ({c_idx + 1}/{len(C_BAND)}): "
                f"{len(cands)} coarse candidate(s) to certify."
            )
            for i, cand in enumerate(cands):
                print(
                    f"[{_ts()}]   certifying candidate {i + 1}/{len(cands)} "
                    f"(seed x0={cand.x0:.5f} xdot0={cand.xdot0:.5f}, "
                    f"coarse residual={cand.residual:.3e}) ..."
                )
                orbit = _certify_with_refinement(system, cand.x0, cand.xdot0, float(c_target), n)
                if orbit is None or not (orbit.converged and orbit.residual <= 1e-9):
                    print(f"[{_ts()}]     did not certify.")
                    continue
                if any(
                    c == float(c_target)
                    and (orbit.x0 - x) ** 2 + (orbit.xdot0 - xd) ** 2
                    < CERTIFIED_DEDUP_RADIUS * CERTIFIED_DEDUP_RADIUS
                    for c, x, xd in certified_seen
                ):
                    n_dup_skipped += 1
                    print(f"[{_ts()}]     duplicate of an already-certified orbit, skipping.")
                    continue
                certified_seen.append((float(c_target), orbit.x0, orbit.xdot0))
                n_certified += 1
                min_dist = _min_dist_to_earth(
                    system, orbit.x0, orbit.xdot0, orbit.ydot0, orbit.period
                )
                is_encounter = min_dist < r_hill  # the written encounter criterion
                if is_encounter:
                    n_encounter += 1
                print(
                    f"[{_ts()}]     CERTIFIED C={c_target:.4f} n={n} x0={orbit.x0:.5f} "
                    f"xdot0={orbit.xdot0:.5f} period={orbit.period:.3f} "
                    f"min_dist_to_earth={min_dist:.5f} (Hill={r_hill:.5f}) "
                    f"{'[ENCOUNTER]' if is_encounter else ''}"
                )
                findings.append(
                    {
                        "c_target": float(c_target),
                        "n": n,
                        "x0": orbit.x0,
                        "xdot0": orbit.xdot0,
                        "ydot0": orbit.ydot0,
                        "period": orbit.period,
                        "residual": orbit.residual,
                        "min_dist_to_earth": min_dist,
                        "is_encounter": is_encounter,
                    }
                )

    dt = time.time() - t0
    print()
    print(
        f"[{_ts()}] Scan complete in {dt:.1f}s. Certified orbits: {n_certified} "
        f"({n_dup_skipped} duplicate(s) skipped); meeting the Hill-sphere encounter "
        f"criterion: {n_encounter}."
    )
    for f in findings:
        print(f"    {f}")


if __name__ == "__main__":
    try:
        main()
    except PreflightBlockedError as exc:
        print(f"[{_ts()}] BLOCKED by preflight_search:\n{exc}")
        sys.exit(1)
