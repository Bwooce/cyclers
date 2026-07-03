"""#532 scoping pilot: does an UNSTABLE-Jacobi-constant overlap even exist
between the Sun-Jupiter 2:1 interior and 3:2 (Hilda) mean-motion-resonance
families, before building the full cross-family heteroclinic-connection
search?

Per the Guido & Efthymiopoulos (arXiv:2604.00679) resonance-hopping
hypothesis, a connection between orbits in DIFFERENT MMR families can only
exist where both families have a genuinely UNSTABLE member at the SAME
Jacobi constant (`assemble_cycle`/`correct_connection` hard-require exact
same-Jacobi nodes -- a physical constraint of this autonomous system, not a
numerical convenience). Each family is a continuous one-parameter arc
sweeping a RANGE of C, and #530 already found only ~17% (1/6 sampled) of
Hilda-family orbits are genuinely unstable -- so the question is not "do the
families overlap in C" (they trivially can, both are wide bands) but "does
EITHER family's UNSTABLE sub-range overlap the other's", which is unknown
and must be checked before any connection-search machinery is built.

This is a SCOPING PILOT (go/no-go), not the full #532 build: it certifies
orbits across a Jacobi band for both families (reusing #527's exact DA/HOTM
enumeration + correction pipeline) and classifies each certified orbit's
stability (reusing #530's exact monodromy/leading-eigenvalue method), then
reports whether the two families' unstable-C sub-ranges intersect at all.

2:1 interior MMR seed derivation (mirrors #527's Hilda derivation exactly,
Kepler III on Jupiter's sourced SMA): heliocentric distance ratio
`r = (1/2)**(2/3)` (period ratio P_ast/P_jup = 1/2 for an interior 2:1
resonance), giving `a_2:1 = r * a_jupiter = 3.278 AU` -- independently
matching the well-documented real "Hecuba gap" location in the asteroid
belt (the 2:1 resonance's own well-known real-world signature, exactly
mirroring the Hilda-group real-world cross-check #527 already used).
Converted to a rotating-frame circular seed via `x0 = r - mu`,
`ydot0 = sqrt((1-mu)/r) - x0` (barycentric offset + rotating-frame velocity
transform at the x-axis crossing). VERIFIED against #527's own committed
seed before use: applying this exact formula to the 3:2 ratio reproduces
#527's `C_SEED=3.0613` to 4 decimal places (3.061320 computed here) --
confirming the derivation itself, not just the already-published Hilda
number, before trusting the new 2:1 seed. The 2:1 natural section-return
period (one full Jupiter rotating-frame period, `T~=6.26` nondim units,
directly measured by integration) is one-half the Hilda family's
`T~=12.6` (two Jupiter periods) -- consistent with the 2:1 vs 3:2
resonance orders, not a system-scale bug of the kind #527 found and fixed.

Positive control: the Hilda-family half of this scan reproduces #527's own
already-published positive control (a sub-tolerance candidate recovered
near the Kepler-derived Hilda seed) before either band's negative/positive
finding is trusted. Every monodromy computed is checked for the symplectic
det=1 property (#530's own control) before its stability classification is
used.
"""

from __future__ import annotations

import dataclasses
import datetime
import pathlib
import sys
import time

_REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "src"))

import numpy as np  # noqa: E402

import cyclerfinder.core.cr3bp as cr3bp  # noqa: E402
from cyclerfinder.data.method_capability import MethodCapability  # noqa: E402
from cyclerfinder.data.preflight import PreflightBlockedError, preflight_search  # noqa: E402
from cyclerfinder.genome.da_hotm_backend import SamplingSectionMap  # noqa: E402
from cyclerfinder.genome.da_hotm_enumerator import DomainBox, enumerate_fixed_points  # noqa: E402
from cyclerfinder.search.bifurcation_detector import monodromy  # noqa: E402
from cyclerfinder.search.cr3bp_general_periodic import correct_general_periodic  # noqa: E402


@dataclasses.dataclass(frozen=True)
class FamilySpec:
    label: str
    period_ratio: float  # P_asteroid / P_jupiter
    t_max: float  # section-map propagation horizon (nondim), covers the natural return
    period_guess_unit: float  # measured single-rev section-return time (nondim)
    domain_half_width_x: float
    domain_half_width_xdot: float
    c_offsets: tuple[float, ...]  # relative to the derived seed C


def _seed_from_period_ratio(mu: float, period_ratio: float) -> tuple[float, float, float]:
    """Kepler-III-derived circular seed -- returns (x0, ydot0, c_seed)."""
    r = period_ratio ** (2.0 / 3.0)
    x0 = r - mu
    v_circ = np.sqrt((1.0 - mu) / r)
    ydot0 = v_circ - x0
    state = np.array([x0, 0.0, 0.0, 0.0, ydot0, 0.0])
    c_seed = cr3bp.jacobi_constant(state, mu)
    return x0, ydot0, c_seed


HILDA = FamilySpec(
    label="3:2 Hilda",
    period_ratio=2.0 / 3.0,
    t_max=20.0,
    period_guess_unit=12.6,
    domain_half_width_x=0.15,
    domain_half_width_xdot=0.08,
    c_offsets=(-0.12, -0.09, -0.06, -0.03, 0.0, 0.03, 0.06, 0.08),
)
INTERIOR_2_1 = FamilySpec(
    label="2:1 interior",
    period_ratio=1.0 / 2.0,
    t_max=10.0,
    period_guess_unit=6.26,
    domain_half_width_x=0.13,
    domain_half_width_xdot=0.08,
    c_offsets=(-0.12, -0.09, -0.06, -0.03, 0.0, 0.03, 0.06, 0.08),
)
GRID = (21, 15)
N_RANGE = (1, 2)
RESIDUAL_TOL = 1e-2
UNSTABLE_THRESHOLD = 1.01  # #530's own margin above 1 to avoid numerical-noise false positives

_METHOD = MethodCapability(
    genome="DA/HOTM global Poincare-section fixed-point enumeration (#450), "
    "generalized across Sun-Jupiter MMR families",
    corrector="correct_general_periodic (analytic-STM single shooting) + monodromy stability "
    "classification (#530's method)",
    capability_tags=frozenset(
        {"cr3bp", "ballistic", "coplanar", "single-arc", "poincare-section-enumeration"}
    ),
    git_sha="working-tree",
)


def _ts() -> str:
    return datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")


def _classify_stability(
    system: cr3bp.CR3BPSystem, state0: np.ndarray, period: float
) -> float | None:
    """Returns the leading non-trivial |eigenvalue|, or None if all trivial."""
    mono = monodromy(system, state0, period)
    det = float(np.linalg.det(mono))
    if abs(det - 1.0) > 1e-3:
        return None  # symplectic control failed; do not trust this classification
    eigvals, _ = np.linalg.eig(mono)
    mags = np.abs(eigvals)
    non_trivial_idx = [
        i for i in range(6) if abs(eigvals[i].real - 1.0) > 1e-3 or abs(eigvals[i].imag) > 1e-6
    ]
    if not non_trivial_idx:
        return 1.0
    lead = max(non_trivial_idx, key=lambda i: mags[i])
    return float(mags[lead])


def _scan_family(system: cr3bp.CR3BPSystem, spec: FamilySpec) -> dict[float, list[dict]]:
    mu = float(system.mu)
    x0_seed, ydot0_seed, c_seed = _seed_from_period_ratio(mu, spec.period_ratio)
    c_band = tuple(round(c_seed + off, 4) for off in spec.c_offsets)
    x_lo = x0_seed - spec.domain_half_width_x
    x_hi = x0_seed + spec.domain_half_width_x
    xdot_lo = -spec.domain_half_width_xdot
    xdot_hi = spec.domain_half_width_xdot
    domain_box = DomainBox(x_lo=x_lo, x_hi=x_hi, xdot_lo=xdot_lo, xdot_hi=xdot_hi)

    print(
        f"[{_ts()}] === {spec.label}: seed x0={x0_seed:.6f} ydot0={ydot0_seed:.6f} "
        f"C_seed={c_seed:.6f}, band={c_band} ==="
    )

    results: dict[float, list[dict]] = {c: [] for c in c_band}
    for c_target in c_band:
        backend = SamplingSectionMap(
            system, c_target=float(c_target), ydot_sign=1.0, t_max=spec.t_max
        )
        for n in N_RANGE:
            cands = enumerate_fixed_points(
                backend, domain_box, n, residual_tol=RESIDUAL_TOL, grid=GRID, dedup_radius=0.01
            )
            for cand in cands:
                period_guess = spec.period_guess_unit * n
                orbit = correct_general_periodic(
                    system,
                    cand.x0,
                    cand.xdot0,
                    float(c_target),
                    period_guess=period_guess,
                    half_crossings=2 * n,
                    ydot0_sign=1.0,
                    tol=1e-11,
                )
                if not (orbit.converged and orbit.residual <= 1e-9):
                    continue
                state0 = np.array(
                    [orbit.x0, 0.0, 0.0, orbit.xdot0, orbit.ydot0, 0.0], dtype=np.float64
                )
                lead_mag = _classify_stability(system, state0, orbit.period)
                if lead_mag is None:
                    continue  # symplectic control failed; skip rather than mis-classify
                is_unstable = lead_mag > UNSTABLE_THRESHOLD
                results[c_target].append(
                    {
                        "x0": orbit.x0,
                        "xdot0": orbit.xdot0,
                        "ydot0": orbit.ydot0,
                        "period": orbit.period,
                        "n": n,
                        "lead_eig_mag": lead_mag,
                        "unstable": is_unstable,
                    }
                )
    n_certified = sum(len(v) for v in results.values())
    n_unstable = sum(1 for v in results.values() for r in v if r["unstable"])
    print(
        f"[{_ts()}] {spec.label}: {n_certified} certified orbits, "
        f"{n_unstable} unstable ({100 * n_unstable / max(n_certified, 1):.1f}%)."
    )
    for c_target, orbits in results.items():
        for r in orbits:
            flag = "UNSTABLE" if r["unstable"] else "stable"
            print(
                f"[{_ts()}]     C={c_target:.4f} n={r['n']} x0={r['x0']:.5f} "
                f"|lambda|={r['lead_eig_mag']:.4f} [{flag}]"
            )
    return results


def _positive_control(system: cr3bp.CR3BPSystem) -> None:
    """Re-verify #527's own committed Hilda positive control before trusting anything."""
    mu = float(system.mu)
    _, _, c_seed = _seed_from_period_ratio(mu, HILDA.period_ratio)
    print(
        f"[{_ts()}] Positive control: re-deriving the Hilda seed and checking it matches "
        f"#527's committed C_SEED=3.0613 ..."
    )
    if abs(c_seed - 3.0613) > 1e-3:
        raise RuntimeError(
            f"POSITIVE CONTROL FAILED: re-derived Hilda seed C={c_seed:.6f} does not match "
            f"#527's committed C_SEED=3.0613. Do not trust the 2:1 seed derivation below -- "
            f"it reuses the identical formula."
        )
    print(f"[{_ts()}] Positive control PASSED: re-derived C_seed={c_seed:.6f} matches #527.")

    control_backend = SamplingSectionMap(system, c_target=3.0613, ydot_sign=1.0, t_max=HILDA.t_max)
    control_box = DomainBox(x_lo=0.70, x_hi=0.82, xdot_lo=-0.05, xdot_hi=0.05)
    control_candidates = enumerate_fixed_points(
        control_backend, control_box, 1, residual_tol=RESIDUAL_TOL, grid=(15, 11)
    )
    if not control_candidates:
        raise RuntimeError(
            "POSITIVE CONTROL FAILED: no sub-tolerance candidate recovered near the Hilda "
            "seed -- the shared enumeration machinery itself is not reproducing #527's "
            "already-published result. Stop; do not trust the pilot below."
        )
    print(
        f"[{_ts()}] Positive control PASSED: {len(control_candidates)} candidate(s) "
        f"near the Hilda seed, best residual={min(c.residual for c in control_candidates):.3e}."
    )


def main() -> None:
    print(f"[{_ts()}] #532 resonance-overlap scoping pilot starting.")
    system = cr3bp.cr3bp_system("Sun", "Jupiter")

    n_points = sum(
        len(spec.c_offsets) * len(N_RANGE) * GRID[0] * GRID[1] for spec in (HILDA, INTERIOR_2_1)
    )
    preflight_search(
        task_no=532,
        region_id="sun-jupiter-21-32-mmr-unstable-c-overlap-pilot",
        method=_METHOD,
        script_path=pathlib.Path(__file__),
        n_points=n_points,
        # #527 measured 0.007-0.043 s/point on the identical enumerator/grid shape for the
        # Hilda family; the 2:1 family reuses the same machinery at a comparable grid size.
        timing_pilot_seconds_per_point=0.05,
    )

    _positive_control(system)

    t0 = time.time()
    hilda_results = _scan_family(system, HILDA)
    interior_results = _scan_family(system, INTERIOR_2_1)
    dt = time.time() - t0

    hilda_unstable_c = sorted(
        c for c, orbits in hilda_results.items() if any(r["unstable"] for r in orbits)
    )
    interior_unstable_c = sorted(
        c for c, orbits in interior_results.items() if any(r["unstable"] for r in orbits)
    )

    print()
    print(f"[{_ts()}] Scan complete in {dt:.1f}s.")
    print(f"[{_ts()}] Hilda (3:2) unstable C values: {hilda_unstable_c}")
    print(f"[{_ts()}] Interior (2:1) unstable C values: {interior_unstable_c}")

    if not hilda_unstable_c or not interior_unstable_c:
        print(
            f"[{_ts()}] VERDICT: at least one family shows NO unstable member in this band -- "
            f"cannot even test overlap at this sampling resolution."
        )
        return

    hilda_range = (min(hilda_unstable_c), max(hilda_unstable_c))
    interior_range = (min(interior_unstable_c), max(interior_unstable_c))
    overlap_lo = max(hilda_range[0], interior_range[0])
    overlap_hi = min(hilda_range[1], interior_range[1])
    overlaps = overlap_lo <= overlap_hi

    print(f"[{_ts()}] Hilda unstable C-range (sampled extent): {hilda_range}")
    print(f"[{_ts()}] Interior unstable C-range (sampled extent): {interior_range}")
    print(
        f"[{_ts()}] VERDICT: unstable-C ranges {'OVERLAP' if overlaps else 'DO NOT OVERLAP'} "
        f"at this sampling resolution."
    )
    if overlaps:
        print(
            f"[{_ts()}] GO: #532's full cross-family connection search is worth building -- "
            f"a shared-energy unstable pair may exist in [{overlap_lo:.4f}, {overlap_hi:.4f}]."
        )
    else:
        print(
            f"[{_ts()}] NO-GO (at this sampling resolution): no evidence of a shared-energy "
            f"unstable pair between the 2:1 and 3:2 families. A finer C-band sample could "
            f"still find one -- this is a sampling-resolution-scoped result, not a proof of "
            f"non-existence."
        )


if __name__ == "__main__":
    try:
        main()
    except PreflightBlockedError as exc:
        print(f"[{_ts()}] BLOCKED by preflight_search:\n{exc}")
        sys.exit(1)
