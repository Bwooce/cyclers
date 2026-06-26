"""Task #477 — VILM quadrature divergent-integral guard + bit-identity.

The leverage quadrature ``vilm._quadrature_dv_adim`` integrates ``V∞/Γ(V∞)``; the
integrand has a pole wherever ``Γ → 0`` (near V∞ ≈ 1.39 for the Exterior leg). A
V∞-band that straddles that root is a *physically infeasible* leverage step, and
``scipy.integrate.quad`` previously ground through its full subdivision budget on
it (the "integral is probably divergent, or slowly convergent" path — minutes per
moon-tour candidate). The #477 fix detects non-convergence via ``full_output`` and
returns the infeasible sentinel ``math.inf`` CHEAPLY, while keeping every
CONVERGENT (feasible) value byte-for-byte identical (the safety proof).
"""

from __future__ import annotations

import math
import time

import pytest
from scipy.integrate import quad

from cyclerfinder.search.vilm import _quadrature_dv_adim, gamma


def test_divergent_band_returns_inf_fast() -> None:
    # A band straddling the Γ pole (~1.39 exterior) makes V∞/Γ divergent. The
    # guard must return the infeasible sentinel, fast (not grind / hang).
    t0 = time.perf_counter()
    val = _quadrature_dv_adim.__wrapped__(1.2, 1.6, exterior=True)
    dt = time.perf_counter() - t0
    assert val == math.inf
    # Cheap: the bounded quad returns in well under a second (the old path took
    # seconds-to-minutes via runaway subdivision). Generous bound for CI noise.
    assert dt < 1.0


def test_convergent_band_is_bit_identical_to_plain_quad() -> None:
    # The safety proof: for every CONVERGENT (feasible) band the guarded value is
    # byte-for-byte identical to the prior plain ``quad(...)`` at the default
    # limit=50 — so no golden / feasible value moves; only divergent cases change.
    for exterior in (True, False):
        for lo, hi in ((0.3, 1.0), (0.5, 1.2), (0.8, 1.35), (1.45, 1.9)):
            ref, _ = quad(lambda v, ext=exterior: v / gamma(v, exterior=ext), lo, hi)
            got = _quadrature_dv_adim.__wrapped__(lo, hi, exterior=exterior)
            assert got == ref  # bit-identical, not approx


def test_empty_band_returns_zero() -> None:
    assert _quadrature_dv_adim.__wrapped__(1.0, 1.0, exterior=True) == 0.0
    assert _quadrature_dv_adim.__wrapped__(1.5, 1.0, exterior=True) == 0.0


def test_cache_collapses_repeated_args() -> None:
    # The #477 memoization: a repeated (lo, hi, exterior) is a VERIFIED cache hit.
    _quadrature_dv_adim.cache_clear()
    a = _quadrature_dv_adim(0.4, 0.9, exterior=True)
    b = _quadrature_dv_adim(0.4, 0.9, exterior=True)
    assert a == b == pytest.approx(_quadrature_dv_adim.__wrapped__(0.4, 0.9, exterior=True))
    assert _quadrature_dv_adim.cache_info().hits >= 1
