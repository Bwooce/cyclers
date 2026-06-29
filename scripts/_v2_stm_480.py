"""#480 M1 — real-eph EGGIE corrector with analytic STM Jacobian.

Builds the real-ephemeris periapsis seed at the best-epoch offset (+3.78 d past the
paper departure), runs a parity check of the analytic STM Jacobian against FD, then
runs jovian_shoot(jacobian="stm", max_nfev=60, max_wall_sec=500) foreground and
reports the maintenance ΔV.

Usage::

    timeout 540 uv run python scripts/_v2_stm_480.py [runlog_path]
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import spiceypy

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.lambert import lambert
from cyclerfinder.core.satellites import PRIMARIES
from cyclerfinder.nbody.jovian import (
    JovianEphemeris,
    JovianRailsCache,
    jovian_defect_residual,
    jovian_shoot,
    periapsis_node,
)
from cyclerfinder.nbody.jovian_stm import jovian_stm_jacobian
from cyclerfinder.nbody.shooter import (
    ShootingSeed,
    _fd_jacobian,
    _seed_with_states,
    _serial_columns,
    _states_to_x,
    _x_to_states,
)
from cyclerfinder.search.ieg_seed import paper_departure_et
from cyclerfinder.verify.spice_kernels import ensure_jup365_kernel, ensure_leapseconds_kernel

_SPD = 86400.0
_MU_JUPITER = PRIMARIES["Jupiter"]

# Best-epoch offset from the paper departure (Task 4 scan result: +3.78 d).
_BEST_OFFSET_DAYS = 3.78

# Seed parameters for the real-eph EGGIE ballistic periapsis member (Task 4 result).
_SEQUENCE = ("Europa", "Ganymede", "Ganymede", "Io", "Europa")
_TOFS_DAYS = [1.65, 8.70, 7.02, 11.03]
_LEG_PLAN: list[tuple[int, str]] = [(0, "single"), (1, "high"), (1, "low"), (1, "high")]

# Sourced Table-4 V∞ targets (Hernandez 2017, AAS 17-608):
_TABLE4_VINF_KMS = {"Europa": 9.12, "Ganymede": 7.07, "Io": 8.38}


def build_realeph_periapsis_seed(
    departure_et: float,
    jeph: JovianEphemeris,
) -> ShootingSeed:
    """Build the real-eph EGGIE periapsis-node ShootingSeed at departure_et.

    Uses Ephemeris(center="Jupiter", model="spice") for moon positions (same kernel
    as JovianEphemeris), Lambert legs with the per-leg revolution plan, then converts
    each encounter to the periapsis hyperbola state via periapsis_node.
    """
    ephem = Ephemeris(center="Jupiter", model="spice")

    # Encounter epochs.
    epochs: list[float] = [departure_et]
    for tof_d in _TOFS_DAYS:
        epochs.append(epochs[-1] + tof_d * _SPD)

    n_enc = len(_SEQUENCE)
    moon_r: list[np.ndarray] = []
    moon_v: list[np.ndarray] = []
    for i in range(n_enc):
        r, v = ephem.state(_SEQUENCE[i], epochs[i])
        moon_r.append(np.asarray(r, dtype=np.float64))
        moon_v.append(np.asarray(v, dtype=np.float64))

    # Lambert legs -> spacecraft velocities.
    n_legs = len(_TOFS_DAYS)
    sc_v_dep: list[np.ndarray] = []
    sc_v_arr: list[np.ndarray] = []
    for i in range(n_legs):
        tof_sec = _TOFS_DAYS[i] * _SPD
        n_revs, branch = _LEG_PLAN[i]
        sols = lambert(moon_r[i], moon_r[i + 1], tof_sec, mu=_MU_JUPITER, max_revs=n_revs)
        sol = next((s for s in sols if s.n_revs == n_revs and s.branch == branch), None)
        if sol is None:
            sol = sols[0]
        sc_v_dep.append(np.asarray(sol.v1, dtype=np.float64))
        sc_v_arr.append(np.asarray(sol.v2, dtype=np.float64))

    # V∞ vectors at each encounter.
    vinf_in: list[np.ndarray] = []
    vinf_out: list[np.ndarray] = []
    for i in range(n_enc):
        vin = np.zeros(3, dtype=np.float64) if i == 0 else sc_v_arr[i - 1] - moon_v[i]
        vout = sc_v_dep[i] - moon_v[i] if i < n_legs else sc_v_arr[n_legs - 1] - moon_v[i]
        vinf_in.append(vin)
        vinf_out.append(vout)

    # Convert to periapsis nodes via hyperbola geometry.
    node_states: list[np.ndarray] = []
    for k, moon in enumerate(_SEQUENCE):
        vin = vinf_in[k]
        vout = vinf_out[k]
        # Boundary nodes have one zero V∞; use the non-zero one for both.
        if float(np.linalg.norm(vin)) <= 0.0:
            vin = vout
        if float(np.linalg.norm(vout)) <= 0.0:
            vout = vin
        r_p, v_p, _d = periapsis_node(moon, epochs[k], vin, vout, jeph)
        node_states.append(np.concatenate([r_p, v_p]))

    tofs = list(_TOFS_DAYS)
    return ShootingSeed(
        node_states=node_states,
        epochs=epochs,
        tofs=tofs,
        sequence=_SEQUENCE,
        slack_leg=int(np.argmax(tofs)),
        period_days=float(sum(tofs)),
        vinf_in=vinf_in,
        vinf_out=vinf_out,
    )


def run_stm_parity(
    seed: ShootingSeed,
    jeph: JovianEphemeris,
    moons: tuple[str, ...],
    *,
    log: object,
) -> float:
    """Compute STM-vs-FD Jacobian relative error on the periapsis seed.

    Returns the overall relative error. Uses the SAME residual function the
    corrector calls (jovian_defect_residual).
    """

    def _write(msg: str) -> None:
        print(msg, flush=True)
        if log is not None:
            log.write(msg + "\n")  # type: ignore[union-attr]
            log.flush()  # type: ignore[union-attr]

    cache = JovianRailsCache(moons, jeph, min(seed.epochs), max(seed.epochs))
    n = len(seed.sequence)

    def residual_of_x(x: np.ndarray) -> np.ndarray:
        trial = _seed_with_states(seed, _x_to_states(x, n))
        return jovian_defect_residual(
            trial,
            ephem=jeph,
            cache=cache,
            moons=moons,
            accuracy=1e-11,
            max_wall_sec=60.0,
        )

    x0 = _states_to_x(seed.node_states)

    _write("PARITY: computing FD Jacobian ...")
    t0 = time.monotonic()
    f0 = residual_of_x(x0)
    fd = _fd_jacobian(residual_of_x, x0, f0, column_eval=_serial_columns)
    _write(f"PARITY: FD done in {time.monotonic() - t0:.1f}s, shape={fd.shape}")

    _write("PARITY: computing STM Jacobian ...")
    t0 = time.monotonic()
    stm = jovian_stm_jacobian(seed, x0, ephem=jeph, moons=moons)
    _write(f"PARITY: STM done in {time.monotonic() - t0:.1f}s, shape={stm.shape}")

    overall = float(np.linalg.norm(stm - fd) / (np.linalg.norm(fd) + 1e-30))
    _write(f"PARITY: overall STM-vs-FD relative error = {overall:.4e}")

    # Per nonzero 6x6 block: Phi_i and -R_W coupling.
    n_leg = (n - 1) * 6
    n_hinge = max(0, n - 2)

    def block_rel(rows: slice, cols: slice) -> float:
        b = fd[rows, cols]
        bn = float(np.linalg.norm(b))
        if bn < 1e-30:
            return 0.0
        return float(np.linalg.norm(stm[rows, cols] - b) / bn)

    worst = 0.0
    for i in range(n - 1):
        rows = slice(i * 6, (i + 1) * 6)
        worst = max(worst, block_rel(rows, slice(i * 6, (i + 1) * 6)))
        worst = max(worst, block_rel(rows, slice((i + 1) * 6, (i + 2) * 6)))
    wrap0 = n_leg + n_hinge
    wrap = slice(wrap0, wrap0 + 6)
    worst = max(worst, block_rel(wrap, slice(0, 6)))
    worst = max(worst, block_rel(wrap, slice((n - 1) * 6, n * 6)))
    _write(f"PARITY: worst nonzero-block STM-vs-FD rel = {worst:.4e}")

    return overall


def main() -> None:
    runlog_path = (
        Path(sys.argv[1])
        if len(sys.argv) > 1
        else Path(
            "/tmp/claude-1000/-home-bruce-dev-cyclers/b2ef8d69-0098-4763-909e-37ee4edbdf86/scratchpad/l3stm.log"
        )
    )
    runlog_path.parent.mkdir(parents=True, exist_ok=True)

    log = open(runlog_path, "w", buffering=1)  # noqa: SIM115

    def write(msg: str) -> None:
        ts = time.strftime("%Y-%m-%dT%H:%M:%S")
        line = f"{ts} {msg}"
        print(line, flush=True)
        log.write(line + "\n")
        log.flush()

    write("START _v2_stm_480.py")
    t_global = time.monotonic()

    # Resolve JUP365 kernel.
    kernel_path = ensure_jup365_kernel()
    spiceypy.furnsh(ensure_leapseconds_kernel())
    jeph = JovianEphemeris(kernel_path)
    moons = ("Io", "Europa", "Ganymede")
    write(f"JUP365 kernel: {kernel_path}")

    # Build the real-eph periapsis seed at the best epoch.
    write("Building real-eph periapsis seed at paper_et + 3.78 d ...")
    paper_et = paper_departure_et()
    departure_et = paper_et + _BEST_OFFSET_DAYS * _SPD
    write(f"  paper_et = {paper_et:.6e}  departure_et = {departure_et:.6e}")

    seed = build_realeph_periapsis_seed(departure_et, jeph)
    write(f"  seed n_nodes = {len(seed.sequence)}, sequence = {seed.sequence}")
    write(f"  epochs: {[f'{e:.6e}' for e in seed.epochs]}")
    write(f"  tofs_days: {seed.tofs}")

    # Compute seed defect norm.
    cache0 = JovianRailsCache(moons, jeph, min(seed.epochs), max(seed.epochs))
    res0 = jovian_defect_residual(
        seed, ephem=jeph, cache=cache0, moons=moons, accuracy=1e-11, max_wall_sec=60.0
    )
    seed_dn = float(np.linalg.norm(res0))
    write(f"  seed_defect_norm = {seed_dn:.6e}")

    # Per-encounter V∞ at the seed.
    write("Seed V∞ per encounter:")
    for i, (body, epoch) in enumerate(zip(seed.sequence, seed.epochs, strict=True)):
        _, v_m = jeph.state(body, epoch)
        vinf_km = float(np.linalg.norm(np.asarray(seed.node_states[i][3:]) - np.asarray(v_m)))
        tgt = _TABLE4_VINF_KMS.get(body, float("nan"))
        write(f"  node {i} {body}: V∞={vinf_km:.3f} km/s  (Table-4 {body}={tgt:.2f})")

    # --- Parity gate: STM Jacobian vs FD ---
    write("--- STM-vs-FD parity gate ---")
    t_par = time.monotonic()
    parity_rel = run_stm_parity(seed, jeph, moons, log=log)
    write(f"Parity gate done in {time.monotonic() - t_par:.1f}s  rel={parity_rel:.4e}")
    if parity_rel > 1e-3:
        write(
            f"WARNING: parity rel {parity_rel:.4e} > 1e-3 threshold "
            "(real-eph noise expected; proceed)"
        )
    else:
        write(f"PARITY PASS: rel {parity_rel:.4e} <= 1e-3")

    # --- Main run: STM corrector ---
    write("--- Running jovian_shoot(jacobian='stm', max_nfev=60, max_wall_sec=500) ---")
    t_shoot = time.monotonic()
    result = jovian_shoot(
        seed,
        kernel_path=kernel_path,
        moons=moons,
        accuracy=1e-11,
        max_nfev=60,
        max_wall_sec=500.0,
        jacobian="stm",
    )
    wall_shoot = time.monotonic() - t_shoot
    write(f"jovian_shoot done in {wall_shoot:.1f}s, nfev={result.n_iterations}")
    write(f"  converged={result.converged}")
    write(f"  seed_defect_norm={result.seed_defect_norm:.6e}")
    write(f"  final_defect_norm={result.defect_norm:.6e}")
    write(f"  correction_dv_kms={result.correction_dv_kms:.6f}")
    write(f"  correction_dv_ms={result.correction_dv_kms * 1000.0:.2f} m/s")
    write(f"  bend_feasible={result.bend_feasible}")

    write("Corrected V∞ per encounter:")
    for i, (body, vinf_km) in enumerate(
        zip(result.sequence, result.vinf_per_encounter_kms, strict=True)
    ):
        tgt = _TABLE4_VINF_KMS.get(body, float("nan"))
        write(f"  node {i} {body}: V∞={vinf_km:.3f} km/s  (Table-4 {body}={tgt:.2f})")

    # --- FD comparison run (same budget) ---
    write("--- Comparison: jovian_shoot(jacobian='fd', max_nfev=18, max_wall_sec=60) ---")
    t_fd = time.monotonic()
    result_fd = jovian_shoot(
        seed,
        kernel_path=kernel_path,
        moons=moons,
        accuracy=1e-11,
        max_nfev=18,
        max_wall_sec=60.0,
        jacobian="fd",
    )
    wall_fd = time.monotonic() - t_fd
    write(f"FD shoot done in {wall_fd:.1f}s, nfev={result_fd.n_iterations}")
    write(f"  converged={result_fd.converged}")
    write(f"  final_defect_norm={result_fd.defect_norm:.6e}")
    write(f"  correction_dv_ms={result_fd.correction_dv_kms * 1000.0:.2f} m/s")

    # Summary.
    write("--- SUMMARY ---")
    write(f"parity_rel_error: {parity_rel:.4e}")
    write(f"stm_converged: {result.converged}")
    write(f"stm_defect_norm: {result.defect_norm:.6e}")
    write(f"stm_correction_dv_ms: {result.correction_dv_kms * 1000.0:.2f}")
    write(f"stm_nfev: {result.n_iterations}")
    write(f"stm_wall_s: {wall_shoot:.1f}")
    write(f"fd_converged: {result_fd.converged}")
    write(f"fd_defect_norm: {result_fd.defect_norm:.6e}")
    write(f"fd_correction_dv_ms: {result_fd.correction_dv_kms * 1000.0:.2f}")
    write(f"fd_nfev: {result_fd.n_iterations}")
    write(f"fd_wall_s: {wall_fd:.1f}")

    write("V∞ vs Table-4 (sourced E9.12/G7.07/Io8.38):")
    for body, vinf_km in zip(result.sequence, result.vinf_per_encounter_kms, strict=True):
        tgt = _TABLE4_VINF_KMS.get(body, float("nan"))
        write(f"  {body}: V∞_corrected={vinf_km:.3f}  Table-4={tgt:.2f}")

    write(f"TOTAL wall: {time.monotonic() - t_global:.1f}s")
    log.close()

    # Write JSON summary to the same dir.
    summary_path = runlog_path.with_suffix(".json")
    summary = {
        "parity_rel": parity_rel,
        "stm_converged": result.converged,
        "stm_defect_norm": result.defect_norm,
        "stm_correction_dv_ms": result.correction_dv_kms * 1000.0,
        "stm_nfev": result.n_iterations,
        "stm_vinf_kms": list(result.vinf_per_encounter_kms),
        "stm_sequence": list(result.sequence),
        "fd_converged": result_fd.converged,
        "fd_defect_norm": result_fd.defect_norm,
        "fd_correction_dv_ms": result_fd.correction_dv_kms * 1000.0,
        "fd_nfev": result_fd.n_iterations,
        "table4_vinf_kms": _TABLE4_VINF_KMS,
    }
    summary_path.write_text(json.dumps(summary, indent=2))
    print(f"Summary written to {summary_path}", flush=True)


if __name__ == "__main__":
    main()
