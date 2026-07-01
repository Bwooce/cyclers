"""Jupiter-central multiple-shooting corrector RUNS on the IEG seed (#480 M1 Task 4a).

The heliocentric ``shooter.shoot`` is hardwired to ``RestrictedNBody`` (central
MU_SUN, ``PLANETS`` perturbers) and so raises ``KeyError 'Io'`` + IAS15
step-collapse on the Jupiter-centred EGGIE seed. :func:`jovian_shoot` propagates
each leg with :class:`JovianRestrictedNBody` instead. These tests assert the
corrector RUNS and returns a finite, non-diverging result — NOT reproduction
convergence (the single-rev seed will not converge well; that is Task 4b's
multi-rev seed + Task 4's re-run).
"""

from __future__ import annotations

import numpy as np
import pytest

rebound = pytest.importorskip("rebound")
pytest.importorskip("spiceypy")

from cyclerfinder.nbody.jovian import (  # noqa: E402
    JovianEphemeris,
    JovianRailsCache,
    JovianRestrictedNBody,
    jovian_defect_residual,
    jovian_shoot,
)
from cyclerfinder.search.ieg_seed import ieg_eggie_seed  # noqa: E402
from cyclerfinder.verify.spice_kernels import ensure_jup365_kernel  # noqa: E402

try:
    _KERNEL: str | None = ensure_jup365_kernel()
except Exception:  # jup365.bsp is local-only (~50 MB, absent in CI) -> skip, don't error
    _KERNEL = None

pytestmark = pytest.mark.skipif(_KERNEL is None, reason="JUP365 kernel not furnished (local-only)")


def test_jovian_leg_propagates_finite() -> None:
    """Fast smoke (non-slow): one Jovian leg of the IEG seed propagates finite.

    The exact failure mode of the heliocentric shooter on this seed is a
    non-terminating / NaN leg; here a single short Jupiter-central leg must
    return finite state with no KeyError.
    """
    seed = ieg_eggie_seed(departure_et=0.0)
    kernel = ensure_jup365_kernel()
    jeph = JovianEphemeris(kernel)
    cache = JovianRailsCache(("Io", "Europa", "Ganymede"), jeph, seed.epochs[0], seed.epochs[1])
    prop = JovianRestrictedNBody()
    arc = prop.propagate(
        np.asarray(seed.node_states[0][:3]),
        np.asarray(seed.node_states[0][3:]),
        seed.epochs[0],
        seed.epochs[1],
        moons=("Io", "Europa", "Ganymede"),
        cache=cache,
        max_wall_sec=30.0,
    )
    assert np.all(np.isfinite(arc.r_km))
    assert np.all(np.isfinite(arc.v_km_s))


def test_jovian_residual_finite_on_ieg_seed() -> None:
    """The full multiple-shooting residual is finite on the IEG seed (no KeyError)."""
    seed = ieg_eggie_seed(departure_et=0.0)
    kernel = ensure_jup365_kernel()
    jeph = JovianEphemeris(kernel)
    cache = JovianRailsCache(("Io", "Europa", "Ganymede"), jeph, min(seed.epochs), max(seed.epochs))
    res = jovian_defect_residual(
        seed,
        ephem=jeph,
        cache=cache,
        moons=("Io", "Europa", "Ganymede"),
        max_wall_sec=30.0,
    )
    assert np.all(np.isfinite(res))
    assert res.size > 0


@pytest.mark.slow
def test_jovian_shoot_runs_on_ieg_seed() -> None:
    """The Jovian corrector RUNS on the IEG seed, terminates, and does not diverge.

    Asserts (Task 4a HARD GATE): finite ``defect_norm``, finite ``seed_defect_norm``,
    and the corrector did not diverge (``defect_norm <= seed_defect_norm``). Does NOT
    assert reproduction convergence — that is Task 4b (multi-rev seed) + the Task 4
    re-run. Bounded budget: small ``max_nfev`` + per-leg wall cap.
    """
    seed = ieg_eggie_seed(departure_et=0.0)
    result = jovian_shoot(
        seed,
        moons=("Io", "Europa", "Ganymede"),
        max_nfev=8,
        max_wall_sec=20.0,
    )
    assert np.isfinite(result.seed_defect_norm)
    assert np.isfinite(result.defect_norm)
    # Corrector did not diverge from the seed.
    assert result.defect_norm <= result.seed_defect_norm + 1e-6
    assert len(result.vinf_per_encounter_kms) == len(seed.sequence)
    assert all(np.isfinite(v) for v in result.vinf_per_encounter_kms)
