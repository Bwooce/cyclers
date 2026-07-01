"""Parity gate: real-eph EGGIE ballistic seed — STM Jacobian vs FD (#480).

The real-ephemeris EGGIE periapsis seed (Hernandez 2017, AAS 17-608, best-epoch
at paper_departure_et + 3.78 d) is a harder test than the ideal-model case: the
JUP365 rails introduce positional noise ~1e-2 km (spline vs SPICE direct query)
that the STM co-integrator sees but the REBOUND corrector averages via the cache.
The tolerance is therefore 1e-3 rel (vs 1e-4 for the ideal model in
``test_jovian_stm.py``).

Sources:
- Table-4 V∞ targets: Hernandez-Jones-Jesick (2017), AAS 17-608.
  Digest: docs/notes/2026-06-26-digest-hernandez-2017-ieg-triple-cyclers-aas-17-608.md
- Best-epoch offset +3.78 d: Task 4 epoch scan
  (docs/notes/2026-06-27-480-ieg-reproduction-verdict.md).
"""

from __future__ import annotations

import numpy as np
import pytest

rebound = pytest.importorskip("rebound")
pytest.importorskip("spiceypy")

from cyclerfinder.core.ephemeris import Ephemeris  # noqa: E402
from cyclerfinder.core.lambert import lambert  # noqa: E402
from cyclerfinder.core.satellites import PRIMARIES  # noqa: E402
from cyclerfinder.nbody.jovian import (  # noqa: E402
    JovianEphemeris,
    JovianRailsCache,
    jovian_defect_residual,
    periapsis_node,
)
from cyclerfinder.nbody.jovian_stm import jovian_stm_jacobian  # noqa: E402
from cyclerfinder.nbody.shooter import (  # noqa: E402
    ShootingSeed,
    _fd_jacobian,
    _seed_with_states,
    _serial_columns,
    _states_to_x,
    _x_to_states,
)
from cyclerfinder.search.ieg_seed import paper_departure_et  # noqa: E402
from cyclerfinder.verify.spice_kernels import (  # noqa: E402
    ensure_jup365_kernel,
    ensure_leapseconds_kernel,
)

_SPD = 86400.0
_MU_JUPITER = PRIMARIES["Jupiter"]
_SEQUENCE = ("Europa", "Ganymede", "Ganymede", "Io", "Europa")
_TOFS_DAYS = [1.65, 8.70, 7.02, 11.03]
_LEG_PLAN: list[tuple[int, str]] = [(0, "single"), (1, "high"), (1, "low"), (1, "high")]
_BEST_OFFSET_DAYS = 3.78

# Sourced Table-4 V∞ targets (Hernandez 2017, AAS 17-608).
_TABLE4_VINF_KMS = {"Europa": 9.12, "Ganymede": 7.07, "Io": 8.38}


def _build_realeph_periapsis_seed(departure_et: float, jeph: JovianEphemeris) -> ShootingSeed:
    """Build the real-eph EGGIE periapsis-node ShootingSeed (mirrored from _v2_stm_480.py)."""
    import spiceypy

    from cyclerfinder.verify.spice_kernels import ensure_leapseconds_kernel

    spiceypy.furnsh(ensure_leapseconds_kernel())
    ephem = Ephemeris(center="Jupiter", model="spice")

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

    vinf_in: list[np.ndarray] = []
    vinf_out: list[np.ndarray] = []
    for i in range(n_enc):
        vin = np.zeros(3, dtype=np.float64) if i == 0 else sc_v_arr[i - 1] - moon_v[i]
        vout = sc_v_dep[i] - moon_v[i] if i < n_legs else sc_v_arr[n_legs - 1] - moon_v[i]
        vinf_in.append(vin)
        vinf_out.append(vout)

    node_states: list[np.ndarray] = []
    for k, moon in enumerate(_SEQUENCE):
        vin = vinf_in[k]
        vout = vinf_out[k]
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


try:
    _KERNEL: str | None = ensure_jup365_kernel()
except Exception:  # jup365.bsp is local-only (~50 MB, absent in CI) -> skip, don't error
    _KERNEL = None

pytestmark = pytest.mark.skipif(_KERNEL is None, reason="JUP365 kernel not furnished (local-only)")


@pytest.fixture(scope="module")
def realeph_periapsis_seed() -> ShootingSeed:
    """Real-eph EGGIE periapsis seed at paper_et + 3.78 d (built once per module)."""
    kernel_path = ensure_jup365_kernel()
    ensure_leapseconds_kernel()
    jeph = JovianEphemeris(kernel_path)
    departure_et = paper_departure_et() + _BEST_OFFSET_DAYS * _SPD
    return _build_realeph_periapsis_seed(departure_et, jeph)


def test_realeph_seed_defect_finite(realeph_periapsis_seed: ShootingSeed) -> None:
    """Fast smoke: the real-eph periapsis seed has a finite residual."""
    seed = realeph_periapsis_seed
    kernel_path = ensure_jup365_kernel()
    jeph = JovianEphemeris(kernel_path)
    moons = ("Io", "Europa", "Ganymede")
    cache = JovianRailsCache(moons, jeph, min(seed.epochs), max(seed.epochs))
    res = jovian_defect_residual(
        seed, ephem=jeph, cache=cache, moons=moons, accuracy=1e-11, max_wall_sec=60.0
    )
    assert np.all(np.isfinite(res)), "seed residual is not finite"
    defect = float(np.linalg.norm(res))
    # Seed defect is in the 4e5 range (characterised in scripts/_v2_stm_480.py).
    assert 1e3 < defect < 1e7, f"unexpected seed defect {defect:.3e}"
    print(f"\nseed_defect_norm={defect:.4e}")


def test_realeph_stm_jacobian_matches_fd(realeph_periapsis_seed: ShootingSeed) -> None:
    """Parity gate: real-eph STM Jacobian matches FD to 1e-3 rel (#480 M1).

    The real-eph tolerance is 1e-3 (vs 1e-4 for the ideal model) because the JUP365
    rails cache introduces ~1e-2 km spline noise at the moon positions that the STM
    co-integrator (which queries SPICE directly) sees but the REBOUND corrector
    averages through the cache. Overall relative error 5.1e-4 was measured in
    the _v2_stm_480.py run (2026-06-30).

    Sources:
    - Residual function: jovian_defect_residual (jovian.py).
    - STM Jacobian: jovian_stm_jacobian (jovian_stm.py).
    - Parity criterion: 1e-3 rel (real-eph model; tighter 1e-4 applies to ideal in
      test_jovian_stm.py::test_stm_jacobian_matches_fd_on_eggie_seed).
    """
    seed = realeph_periapsis_seed
    kernel_path = ensure_jup365_kernel()
    jeph = JovianEphemeris(kernel_path)
    moons = ("Io", "Europa", "Ganymede")
    n = len(seed.sequence)
    cache = JovianRailsCache(moons, jeph, min(seed.epochs), max(seed.epochs))

    def residual_of_x(x: np.ndarray) -> np.ndarray:
        trial = _seed_with_states(seed, _x_to_states(x, n))
        return jovian_defect_residual(
            trial, ephem=jeph, cache=cache, moons=moons, accuracy=1e-11, max_wall_sec=60.0
        )

    x0 = _states_to_x(seed.node_states)
    f0 = residual_of_x(x0)
    fd = _fd_jacobian(residual_of_x, x0, f0, column_eval=_serial_columns)
    stm = jovian_stm_jacobian(seed, x0, ephem=jeph, moons=moons)

    assert stm.shape == fd.shape, f"shape mismatch: {stm.shape} vs {fd.shape}"

    overall = float(np.linalg.norm(stm - fd) / (np.linalg.norm(fd) + 1e-30))
    print(f"\nreal-eph STM-vs-FD overall rel = {overall:.4e}")
    assert overall < 1e-3, (
        f"real-eph STM-vs-FD overall rel {overall:.4e} > 1e-3 "
        "(real-eph noise budget; measured ~5e-4 in _v2_stm_480.py)"
    )

    # Per nonzero block (less strict at real-eph noise level).
    n_leg = (n - 1) * 6
    n_hinge = max(0, n - 2)
    wrap0 = n_leg + n_hinge

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
    wrap = slice(wrap0, wrap0 + 6)
    worst = max(worst, block_rel(wrap, slice(0, 6)))
    worst = max(worst, block_rel(wrap, slice((n - 1) * 6, n * 6)))
    print(f"real-eph worst nonzero-block rel = {worst:.4e}")
    # Worst-block tolerance is 5e-3 (real-eph noise is concentrated in some blocks).
    assert worst < 5e-3, (
        f"real-eph worst block STM-vs-FD rel {worst:.4e} > 5e-3 "
        "(measured ~3.6e-3 in _v2_stm_480.py)"
    )
