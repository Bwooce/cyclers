"""#140 review remediation — free-return corrector convergence + guards.

Covers the C1 / I1 fixes in :mod:`cyclerfinder.search.free_return`:

* C1 — convergence is decided by residual MAGNITUDE alone, BY DESIGN: a
  max_nfev-exhausted least_squares run (``solver_success is False``) that still
  landed a residual-good point is reported ``converged`` (the residual is the
  physics; trf's internal stopping is secondary). Pinned with a synthetic
  ``max_nfev=1`` run so the solver cannot report success but the seeded geometry
  is already a residual-zero point.
* I1 — the bodies tuple must be ordered ``(inner, outer)`` by semi-major axis;
  out-of-order bodies raise ``ValueError``.

NON-GOLDEN: every asserted number is OUR computation on the circular model, used
as a mechanics/regime fixture, never a published-anchor EXPECTED.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest
import yaml  # type: ignore[import-untyped]
from scipy.optimize import least_squares as _scipy_least_squares

import cyclerfinder.search.free_return as fr
from cyclerfinder.core.constants import DAYS_PER_JULIAN_YEAR, MU_SUN_KM3_S2
from cyclerfinder.core.ephemeris import Ephemeris

REPO_ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = REPO_ROOT / "scripts" / "campaign_russell12.py"
DAY_S = 86400.0

# Sourced S1L1 heliocentric ellipse (Rogers 2012 Table 1) — the constraint side.
_S1L1_A_AU = 1.30
_S1L1_E = 0.257


def _load_campaign() -> ModuleType:
    spec = importlib.util.spec_from_file_location("campaign_russell12", CAMPAIGN)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _row(rid: str) -> dict[str, Any]:
    rows = yaml.safe_load((REPO_ROOT / "data" / "catalogue.yaml").read_text())
    return next(r for r in rows if r["id"] == rid)


def _best_phase_t0(a_au: float, e: float, period_sec: float, *, n: int = 360) -> float:
    """Phase (t0) minimising the free-return residual at the SOURCED (a, e)."""
    ephem = Ephemeris("circular")
    best_t0, best_res = 0.0, float("inf")
    for frac in np.linspace(0.0, 1.0, n, endpoint=False):
        t0 = float(frac) * period_sec
        res = fr._residuals(
            np.array([a_au, e, t0]),
            period_days=period_sec / DAY_S,
            ephem=ephem,
            bodies=("E", "M"),
            mu=MU_SUN_KM3_S2,
        )
        m = max(abs(r) for r in res)
        if m < best_res:
            best_res, best_t0 = m, t0
    return best_t0


# ---------------------------------------------------------------------------
# C1 — residual-magnitude-only convergence (max_nfev exhausted but good)
# ---------------------------------------------------------------------------


def test_converged_is_residual_only_even_when_solver_does_not_succeed(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """C1: a max_nfev-exhausted least_squares run (``solver_success is False``)
    that nonetheless sits at a residual-zero point must report ``converged`` —
    the residual is authoritative BY DESIGN; the solver's internal stopping is
    secondary and recorded only for the audit trail.

    Construction: first run a normal (unstarved) solve to obtain the converged
    ``(a, e, t0)``; then re-seed EXACTLY there but starve the solver to
    ``max_nfev=1`` so it cannot declare success. The residual is already below the
    floor at the seed, so no refinement is needed — exactly the
    exhausted-but-residual-good case the review asked to pin."""
    period_sec = 4.27 * DAYS_PER_JULIAN_YEAR * DAY_S
    t0 = _best_phase_t0(_S1L1_A_AU, _S1L1_E, period_sec)

    # Step 1: normal solve -> the residual-zero geometry.
    converged = fr.free_return_correct(
        t0_seed_sec=t0,
        a_seed_au=_S1L1_A_AU,
        e_seed=_S1L1_E,
        period_sec=period_sec,
        ephem=Ephemeris("circular"),
        tol_kms=0.1,
    )
    assert converged.converged
    assert converged.max_residual_kms < 0.1

    # Step 2: re-seed AT the residual-zero point, but starve the solver so it
    # cannot report success.
    def _starved(fun: Any, x0: Any, **kw: Any) -> Any:
        kw = {**kw, "max_nfev": 1}
        return _scipy_least_squares(fun, x0, **kw)

    monkeypatch.setattr("cyclerfinder.search.free_return.least_squares", _starved)

    result = fr.free_return_correct(
        t0_seed_sec=converged.t0_sec,
        a_seed_au=converged.a_au,
        e_seed=converged.e,
        period_sec=period_sec,
        ephem=Ephemeris("circular"),
        tol_kms=0.1,
    )

    # The solver itself did not declare success (max_nfev exhausted)...
    assert result.solver_success is False
    assert result.solver_nfev >= 1
    # ...but the residual IS below the floor, so the physics says converged.
    assert result.max_residual_kms < 0.1
    assert result.converged is True


def test_not_converged_when_residual_above_floor() -> None:
    """The flip side: a geometry whose residual stays above the floor is NOT
    converged regardless of the solver outcome (the gate keeps teeth)."""
    period_sec = 4.27 * DAYS_PER_JULIAN_YEAR * DAY_S
    t0 = _best_phase_t0(_S1L1_A_AU, _S1L1_E, period_sec)
    result = fr.free_return_correct(
        t0_seed_sec=t0,
        a_seed_au=_S1L1_A_AU,
        e_seed=_S1L1_E,
        period_sec=period_sec,
        ephem=Ephemeris("circular"),
        tol_kms=1e-9,  # an impossibly tight floor: nothing should pass
    )
    assert result.converged is False
    assert result.max_residual_kms >= 1e-9


# ---------------------------------------------------------------------------
# I1 — bodies-order guard
# ---------------------------------------------------------------------------


def test_bodies_must_be_inner_outer_order_geometry() -> None:
    """``free_return_geometry`` rejects an (outer, inner) bodies tuple."""
    with pytest.raises(ValueError, match="inner, outer"):
        fr.free_return_geometry(1.3, 0.257, bodies=("M", "E"))


def test_bodies_must_be_inner_outer_order_corrector() -> None:
    """``free_return_correct`` rejects an out-of-order bodies tuple up front
    (before the residual swallows the ValueError)."""
    with pytest.raises(ValueError, match="inner, outer"):
        fr.free_return_correct(
            t0_seed_sec=0.0,
            a_seed_au=1.3,
            e_seed=0.257,
            period_sec=4.27 * DAYS_PER_JULIAN_YEAR * DAY_S,
            ephem=Ephemeris("circular"),
            bodies=("M", "E"),
        )


def test_bodies_inner_outer_order_accepted() -> None:
    """The correct (E, M) ordering is accepted (no raise)."""
    g = fr.free_return_geometry(1.3, 0.257, bodies=("E", "M"))
    assert g.vinf["E"] > 0.0
    assert g.vinf["M"] > 0.0


# ---------------------------------------------------------------------------
# Truth-representability beyond k=2 (the #137 free-return breakthrough)
# ---------------------------------------------------------------------------


@pytest.mark.slow
@pytest.mark.parametrize("rid", ["russell-ch4-5.30gGf3", "russell-ch4-9.94Gg3"])
def test_free_return_converges_at_sourced_ellipse_beyond_k2(rid: str) -> None:
    """#140 regression: seeded at the SOURCED aphelion/transit ellipse and best
    phase, the free-return corrector reaches residual ~ 0 on the higher-k Russell
    rows (not just the symmetric k=2 ones). Circular model, fast.

    GOLDEN DISCIPLINE: (a, e) is the sourced input (constraint); the asserted
    EVIDENCE is residual convergence on the reconstructed arc, not a rediscovered
    published value."""
    camp = _load_campaign()
    row = _row(rid)
    aphelion = row["orbit_elements"]["aphelion_au"]
    transit = row["invariants"]["transit_times_days"]
    period_sec = float(row["period"]["years"]) * DAYS_PER_JULIAN_YEAR * DAY_S
    a_seed, e_seed = camp._seed_ae_from_aphelion_transit(float(aphelion), float(transit[0]))

    t0 = _best_phase_t0(a_seed, e_seed, period_sec, n=camp.FR_PHASE_EPOCHS_FLOOR)
    result = fr.free_return_correct(
        t0_seed_sec=t0,
        a_seed_au=a_seed,
        e_seed=e_seed,
        period_sec=period_sec,
        ephem=Ephemeris("circular"),
        tol_kms=0.1,
    )
    assert result.converged, f"{rid}: residual {result.max_residual_kms} km/s"
    assert result.max_residual_kms < 0.1
