"""Task #665: cannonball-SRP-augmented real-binary (k1,k2) re-sweep.

Re-sweeps the SIX real binary-asteroid systems `#549`/`#657`/`#659`/`#660`
found gravity-only-negative (Patroclus-Menoetius, Didymos-Dimorphos,
Orcus-Vanth, Eris-Dysnomia, Sila-Nunam, Lempo-Hiisi -- NOT Antiope, already
closed as physically inadmissible regardless of force model per `#659`) with
a cannonball SRP perturbation added (`cyclerfinder.search.real_binary_srp`),
at the SAME six `(k1,k2)` topologies `#549`/`#657` used ((1,1) from all 3
Table-I anchors, (3,1), (3,2), (3,3) anchor-seeded; (2,1),(2,2) grid-seeded),
via the new `sweep_family_srp`/`sweep_family_grid_srp` drivers in
`real_binary_kk_sweep.py`.

Scope decisions (all documented, none silent)
------------------------------------------------
* beta = 1e-3 m^2/kg -- the TOP of the user-decided bare-rock range (1e-4 to
  1e-3 m^2/kg), i.e. the largest, most-likely-to-show-an-effect area-to-mass
  ratio in scope. C_R = 1.3 (representative bare-rock/regolith reflectivity,
  between fully-absorbing 1.0 and fully-reflecting 2.0).
* phi0 in {0, pi} -- SRP exactly along the primary-secondary axis. This is
  NOT an arbitrary "worst case" pick: `real_binary_srp`'s own docstring
  ("IMPORTANT constraint on phi0") shows any OTHER angle breaks the
  half-period symmetric-orbit corrector's mirror-symmetry assumption and
  produces a numerically-"converged" but NOT genuinely periodic result
  (confirmed directly: full-period Radau crosscheck closes to ~1e-8 at
  phi0 in {0,pi} vs ~1e-3--1e-4 off-axis). {0, pi} are also the two
  physically distinct geometries (Sun on the primary side vs the secondary
  side of the mutual orbit) -- not a cherry-pick, the only valid AND a fair
  bracket of the real geometry.
* Per-system heliocentric distance (needed for a_srp ~ 1/d_AU^2): standard
  JPL Small-Body Database osculating semi-major axes (order-of-magnitude
  sourcing is sufficient here -- the negligibility argument for 5/6 systems,
  see SRP_TO_GRAVITY_RATIO below, is robust to a factor of a few in AU).

Mandatory positive control (run first, gates trusting anything below):
`sweep_32_positive_control()` (#504's own function, UNMODIFIED) PLUS a
beta=0 regression check of the NEW `sweep_family_srp` driver against that
same PC(3,2) row (confirms the new machinery reduces exactly to the old,
already-trusted machinery at beta=0).

Usage
-----
  uv run python scripts/run_665_srp_binary_resweep.py --phase positive_control
  uv run python scripts/run_665_srp_binary_resweep.py --phase ratio_table
  uv run python scripts/run_665_srp_binary_resweep.py --phase anchors
  uv run python scripts/run_665_srp_binary_resweep.py --phase grid
"""

from __future__ import annotations

import argparse
import math
import sys
import time
from pathlib import Path

_SRC = Path(__file__).parent.parent / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import numpy as np  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import cyclerfinder.search.real_binary_srp as srp  # noqa: E402
from cyclerfinder.data.method_capability import MethodCapability  # noqa: E402
from cyclerfinder.data.preflight import preflight_search  # noqa: E402
from cyclerfinder.search.pluto_charon_kk_sweep import (  # noqa: E402
    PC_MU,
    _c_l1,
    sweep_32_positive_control,
)
from cyclerfinder.search.real_binary_kk_sweep import (  # noqa: E402
    REAL_BINARY_SYSTEMS,
    SweepResult,
    sweep_family_grid_srp,
    sweep_family_srp,
)

OUT_PATH = Path(__file__).parent.parent / "docs" / "notes" / "scratch" / "665_srp_resweep_raw.txt"

_REGION_ID = "real-binary-kk-cycler-srp-resweep-2026-07-21"
_METHOD = MethodCapability(
    genome=(
        "Real-binary (k1,k2) CR3BP cycler genome (#494/#504/#549's fixed-Jacobi symmetric "
        "corrector + Barden stability + winding-topology classifier + independent Radau "
        "crosscheck) AUGMENTED with cannonball SRP (#665, real_binary_srp.py): gravity-only "
        "seed (mu_step_to_system / grid search, unchanged) then beta-continuation + a "
        "C_srp-sweep for a stable window, at the same 6 (k1,k2) topologies #549/#657 used"
    ),
    corrector=(
        "SRP-augmented fixed-Jacobi single shooting (correct_symmetric_fixed_jacobi_srp), "
        "phi0 in {0,pi} only (mirror-symmetry constraint -- see real_binary_srp.py docstring)"
    ),
    capability_tags=frozenset(
        {"cr3bp", "binary-cycler", "k1k2-genome", "real-binary", "mu-continuation", "srp"}
    ),
    git_sha="working-tree",
)

ANCHOR_TOPOLOGIES: list[tuple[int, int, str]] = [
    (1, 1, "mu001_11"),
    (1, 1, "mu01215_11"),
    (1, 1, "mu05_11"),
    (3, 2, "mu01_32"),
    (3, 1, "mu03_31"),
    (3, 3, "mu01215_33"),
]
GRID_TOPOLOGIES: list[tuple[int, int]] = [(2, 1), (2, 2)]

SYSTEM_KEYS = [
    "patroclus-menoetius",
    "didymos-dimorphos",
    "orcus-vanth",
    "eris-dysnomia",
    "sila-nunam",
    "lempo-hiisi",
]

#: Standard JPL Small-Body Database osculating heliocentric semi-major axis
#: (AU), order-of-magnitude sourcing (see module docstring). Used only for
#: the a_srp ~ 1/d_AU^2 scaling.
SUN_DISTANCE_AU: dict[str, float] = {
    "patroclus-menoetius": 5.23,  # Jupiter L5 Trojan
    "didymos-dimorphos": 1.64,  # NEA (Apollo-class), Hera target
    "orcus-vanth": 39.2,  # Plutino (3:2 Neptune resonance)
    "eris-dysnomia": 67.8,  # scattered-disk dwarf planet
    "sila-nunam": 43.9,  # cold classical KBO
    "lempo-hiisi": 39.4,  # Plutino
}

#: Bare-rock top-of-range beta (m^2/kg) and a representative cannonball
#: reflectivity coefficient for an uncoated rocky/metallic body (1.0 fully
#: absorbing, 2.0 fully specularly reflecting -- 1.3 is a standard assumed
#: mid-range value, e.g. widely used as a default in mission SRP models).
BETA_M2_PER_KG = 1.0e-3
C_R = 1.3

PHI0_VALUES: list[float] = [0.0, math.pi]


def _fmt(r: SweepResult) -> str:
    if r.stable_found:
        return (
            f"STABLE  C={r.jacobi_mid:.7f}  x0={r.x0_mid:.9f}  "
            f"T={r.period_mid:.5f} TU ({r.period_days:.3f} d)  "
            f"nu={r.nu_mid:.2e}  topo_ok={r.topology_ok}  xcheck={r.crosscheck_ok} "
            f"(dj={r.crosscheck_dj:.2e})  clearance_ok={r.min_clearance_ok}  "
            f"method={r.method!r}"
        )
    return f"negative  method={r.method!r}  note={r.note!r}"


def _append(lines: list[str]) -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUT_PATH.open("a") as f:
        f.write("\n".join(lines) + "\n")


def _print_and_append(lines: list[str]) -> None:
    for line in lines:
        print(line, flush=True)
    _append(lines)


def phase_positive_control() -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y-%m-%dT%H:%M:%SZ")
    OUT_PATH.write_text(f"Task #665 SRP real-binary re-sweep -- run {stamp}\n")

    print("Step 1: gravity-only positive control (PC 3,2, #504's own function)...", flush=True)
    t0 = time.time()
    r_grav = sweep_32_positive_control()
    _print_and_append(
        [
            f"POSITIVE CONTROL A (gravity-only PC 3,2 @ mu={PC_MU})  {_fmt(r_grav)}  "
            f"[{time.time() - t0:.1f}s]"
        ]
    )

    print(
        "\nStep 2: SRP-driver beta=0 regression check against the SAME PC(3,2) row...",
        flush=True,
    )
    from cyclerfinder.search.pluto_charon_kk_sweep import make_pluto_charon_system

    target = make_pluto_charon_system()
    t0 = time.time()
    r_srp0 = sweep_family_srp(target, "mu01_32", beta_nd=0.0, phi0=0.0)
    _print_and_append(
        [
            f"POSITIVE CONTROL B (sweep_family_srp, beta_nd=0)  {_fmt(r_srp0)}  "
            f"[{time.time() - t0:.1f}s]"
        ]
    )
    match = (
        r_srp0.stable_found
        and r_grav.stable_found
        and abs((r_srp0.jacobi_mid or 0) - (r_grav.jacobi_mid or 0)) < 1e-8
        and abs((r_srp0.x0_mid or 0) - (r_grav.x0_mid or 0)) < 1e-6
    )
    _print_and_append([f"POSITIVE CONTROL B<->A agreement: {match}"])


def phase_ratio_table() -> None:
    """Report the a_srp/a_grav_char ratio for every system at BETA_M2_PER_KG
    -- the honest a-priori scoping check this task's own OUTSTANDING bullet
    demands before spending compute (see module docstring)."""
    lines = [f"\nSRP/gravity ratio table (beta={BETA_M2_PER_KG:.0e} m^2/kg, C_R={C_R:.2f}):"]
    for key in SYSTEM_KEYS:
        s = REAL_BINARY_SYSTEMS[key]
        system = s.to_cr3bp_system()
        d_au = SUN_DISTANCE_AU[key]
        a_srp = srp.cannonball_srp_accel_m_s2(BETA_M2_PER_KG, C_R, d_au)
        beta_nd = srp.srp_beta_nd(BETA_M2_PER_KG, C_R, d_au, system)
        l_m = system.l_km * 1000.0
        a_grav_char = l_m / (system.t_s**2)
        ratio = a_srp / a_grav_char
        lines.append(
            f"  {key:22s} d={d_au:6.2f}AU  a_srp={a_srp:.3e} m/s^2  "
            f"a_grav_char={a_grav_char:.3e} m/s^2  ratio={ratio:.3e}  beta_nd={beta_nd:.3e}"
        )
    _print_and_append(lines)


def _beta_nd_for(sys_key: str) -> float:
    s = REAL_BINARY_SYSTEMS[sys_key]
    system = s.to_cr3bp_system()
    d_au = SUN_DISTANCE_AU[sys_key]
    return srp.srp_beta_nd(BETA_M2_PER_KG, C_R, d_au, system)


def _run_anchor_job(
    sys_key: str, k1: int, k2: int, anchor_key: str, phi0: float
) -> tuple[str, int, int, str, float, SweepResult, float]:
    s = REAL_BINARY_SYSTEMS[sys_key]
    target = s.to_cr3bp_system()
    beta_nd = _beta_nd_for(sys_key)
    t0 = time.time()
    r = sweep_family_srp(
        target,
        anchor_key,
        beta_nd,
        phi0,
        radius_km_primary=s.radius_km_primary,
        radius_km_secondary=s.radius_km_secondary,
    )
    return (sys_key, k1, k2, anchor_key, phi0, r, time.time() - t0)


def _run_grid_job(
    sys_key: str, k1: int, k2: int, phi0: float
) -> tuple[str, int, int, str, float, SweepResult, float]:
    s = REAL_BINARY_SYSTEMS[sys_key]
    target = s.to_cr3bp_system()
    mu = target.mu
    beta_nd = _beta_nd_for(sys_key)
    x0_grid = np.linspace(-1.05, (1.0 - mu) - 0.05, 8)
    c_l1 = _c_l1(mu)
    c_grid = np.linspace(max(2.6, c_l1 - 1.0), c_l1 - 0.01, 6)
    hc_list = (2, 3, 4)
    t0 = time.time()
    r = sweep_family_grid_srp(
        target,
        k1,
        k2,
        beta_nd,
        phi0,
        x0_grid=x0_grid,
        c_grid=c_grid,
        hc_list=hc_list,
        period_guess=14.0,
        per_call_timeout=2,
        radius_km_primary=s.radius_km_primary,
        radius_km_secondary=s.radius_km_secondary,
    )
    return (sys_key, k1, k2, "grid_search", phi0, r, time.time() - t0)


def _emit(results: list[tuple[str, int, int, str, float, SweepResult, float]]) -> None:
    lines = ["\nResults"]
    n_gate_passing = 0
    for sys_key, k1, k2, method_key, phi0, r, elapsed in results:
        gate_pass = (
            r.stable_found
            and r.topology_ok
            and r.crosscheck_ok
            and (r.min_clearance_ok is not False)
        )
        if gate_pass:
            n_gate_passing += 1
        lines.append(
            f"[{sys_key}] ({k1},{k2}) [{method_key}, phi0={phi0:.4f}]  gate_pass={gate_pass}  "
            f"{_fmt(r)}  [{elapsed:.1f}s]"
        )
    lines.append(f"\nTOTAL gate-passing: {n_gate_passing}/{len(results)}")
    _print_and_append(lines)


def phase_anchors() -> None:
    jobs = [
        delayed(_run_anchor_job)(sys_key, k1, k2, anchor_key, phi0)
        for sys_key in SYSTEM_KEYS
        for k1, k2, anchor_key in ANCHOR_TOPOLOGIES
        for phi0 in PHI0_VALUES
    ]
    print(f"\nAnchor-seeded SRP sweep: dispatching {len(jobs)} jobs in parallel...", flush=True)
    t0 = time.time()
    results = Parallel(n_jobs=-1, verbose=10)(jobs)
    print(f"anchor jobs done in {time.time() - t0:.1f}s")
    _emit(results)


def phase_grid() -> None:
    jobs = [
        delayed(_run_grid_job)(sys_key, k1, k2, phi0)
        for sys_key in SYSTEM_KEYS
        for k1, k2 in GRID_TOPOLOGIES
        for phi0 in PHI0_VALUES
    ]
    print(f"\nGrid-seeded SRP sweep: dispatching {len(jobs)} jobs in parallel...", flush=True)
    t0 = time.time()
    results = Parallel(n_jobs=-1, verbose=10)(jobs)
    print(f"grid jobs done in {time.time() - t0:.1f}s")
    _emit(results)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--phase",
        choices=["positive_control", "ratio_table", "anchors", "grid"],
        required=True,
    )
    args = ap.parse_args()

    n_points = (
        1  # gravity-only positive control
        + 1  # SRP-driver beta=0 regression
        + len(SYSTEM_KEYS) * len(ANCHOR_TOPOLOGIES) * len(PHI0_VALUES)
        + len(SYSTEM_KEYS) * len(GRID_TOPOLOGIES) * len(PHI0_VALUES) * 8 * 6 * 3
    )
    preflight_search(
        task_no=665,
        region_id=_REGION_ID,
        method=_METHOD,
        script_path=Path(__file__),
        n_points=n_points,
        override_reason=(
            "reuses #494/#504/#549/#657's already-validated binary-cycler harness for the "
            "GRAVITY-ONLY seed step verbatim; the SRP-augmented continuation+C-sweep is bounded "
            "by construction (fixed anchor/grid topology lists, 2 phi0 values, 1 beta value), "
            "not an open-ended discovery grid needing a timing pilot; positive control (both "
            "the pre-existing gravity-only PC(3,2) AND a new SRP-driver beta=0 regression "
            "against the same row) gates trusting anything below"
        ),
    )

    if args.phase == "positive_control":
        phase_positive_control()
    elif args.phase == "ratio_table":
        phase_ratio_table()
    elif args.phase == "anchors":
        phase_anchors()
    elif args.phase == "grid":
        phase_grid()

    print(f"\n[phase={args.phase}] appended to {OUT_PATH}")


if __name__ == "__main__":
    main()
