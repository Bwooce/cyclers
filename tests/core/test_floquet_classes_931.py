"""#931: planar symplectic Floquet classifier (``core.floquet_classes``), sourced controls only.

Sources:

* Hadjidemetriou 1975b (Celest. Mech. 12:255) eqs. 39-47 for alpha, beta, Delta, b1, b2; digest
  ``docs/notes/2026-10-04-digest-hadjidemetriou-1975b-stability-periodic-orbits-three-body.md``.
  Its Table I rows 10 and 11 need a general three-body integrator, which the repo does not yet have
  as a helper: NOT tested here (follow-up).
* Broucke 1969 (JPL TR 32-1360; AIAA J. 7:1003) printed stability regions per family, digest
  ``docs/notes/2026-10-04-digest-broucke-1969b-aiaa-stability-elliptic-periodic-orbits.md``. The
  rows below are printed starts (Tables 12-17 of the TR). Each is re-corrected here (symmetric
  half-period crossing, free x0 and ydot0; the printed 7-digit state sits within 1e-3 of the
  region boundary k = 2 for the 7 families, and large multipliers amplify the rounding for 8 and
  11), then the one-period 4x4 block is classified. The re-correction is a permitted step, not a
  value under test; the expected regions are the printed ones.
* Jorba & Olle 2004 (Nonlinearity 17:691) map Ts at K = -1, closed-form Jacobian, printed
  alpha = -3, beta = 4 + L, Delta = 1 - 4L, L_crit = 1/4, omega_crit = arccos(3/4) =
  0.72273424781342, Table 1 angles; digest ``docs/notes/2026-10-05-digest-jorba-olle-2004-
  invariant-curves-hamiltonian-hopf.md``.

Not done: the Olle, Pacha & Villanueva 2004 vertical-L4 family (Delta = -7.2e-3 at h = -1.4572,
+4.1e-3 at h = -1.4571; the critical h = -1.4571360299 is the digest agent's derived value, not a
printed one) needs a spatial CR3BP continuation that does not fit a 5 s test: a follow-up.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from numpy.typing import NDArray
from scipy.linalg import expm
from scipy.optimize import fsolve

from cyclerfinder.core.er3bp import ER3BPSystem, propagate_er3bp
from cyclerfinder.core.floquet_classes import (
    classify_invariants,
    classify_planar_monodromy,
    find_family_transitions,
    planar_block,
)

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")

_J = np.block([[np.zeros((2, 2)), np.eye(2)], [-np.eye(2), np.zeros((2, 2))]])


def _pair(k: float) -> NDArray[np.float64]:
    """2x2 block with det 1 and trace 2k: multipliers lambda + 1/lambda = 2k."""
    return np.array([[0.0, -1.0], [1.0, 2.0 * k]])


def _sympl(blocks: tuple[float, float], seed: int = 0) -> NDArray[np.float64]:
    """Symplectic 4x4 with pair indices k = (k1, k2), randomly symplectically conjugated."""
    m = np.zeros((4, 4))
    m[np.ix_((0, 2), (0, 2))] = _pair(blocks[0])
    m[np.ix_((1, 3), (1, 3))] = _pair(blocks[1])
    rng = np.random.default_rng(seed)
    sym = rng.normal(size=(4, 4))
    s = expm(0.4 * _J @ (sym + sym.T) / 2.0)
    assert np.allclose(s.T @ _J @ s, _J, atol=1e-10)
    return s @ m @ np.linalg.inv(s)


def _quartet(a: float, w: float, seed: int = 0) -> NDArray[np.float64]:
    """Symplectic 4x4 with multipliers exp(+-a +- i w), conjugated."""
    big = np.array([[a, w], [-w, a]])
    ham = np.block([[big, np.zeros((2, 2))], [np.zeros((2, 2)), -big.T]])
    m = expm(ham)
    assert np.allclose(m.T @ _J @ m, _J, atol=1e-10)
    rng = np.random.default_rng(seed)
    sym = rng.normal(size=(4, 4))
    s = expm(0.4 * _J @ (sym + sym.T) / 2.0)
    return s @ m @ np.linalg.inv(s)


@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize(
    ("ks", "regime", "region"),
    [
        ((0.3, -0.6), "stable", 1),
        ((1.25, 0.2), "saddle_centre", 6),
        ((-1.25, 0.2), "saddle_centre", 7),
        ((1.25, 3.0), "saddle_saddle", 4),
        ((-1.25, -3.0), "saddle_saddle", 5),
        ((1.25, -3.0), "saddle_saddle", 3),
    ],
)
def test_synthetic_spectra_cover_every_regime(
    ks: tuple[float, float], regime: str, region: int, seed: int
) -> None:
    c = classify_planar_monodromy(_sympl(ks, seed))
    assert (c.regime, c.region) == (regime, region)
    assert not c.boundary
    assert sorted((c.k1.real, c.k2.real)) == pytest.approx(sorted(ks), rel=1e-6)
    assert c.alpha == pytest.approx(-2.0 * (ks[0] + ks[1]), rel=1e-6)
    assert c.reciprocity_error < 1e-8


@pytest.mark.parametrize("seed", [0, 1])
def test_synthetic_complex_quartet(seed: int) -> None:
    c = classify_planar_monodromy(_quartet(0.4, 0.9, seed))
    assert (c.regime, c.region) == ("complex_quartet", 2)
    assert c.delta < 0.0
    assert not c.stable


def test_quartet_with_tiny_modulus_excess_is_not_stable() -> None:
    """Moduli exp(+-5e-4) are inside a 1e-3 modulus tolerance; Delta < 0 still says unstable."""
    m = _quartet(5e-4, 0.9)
    assert float(np.abs(np.linalg.eigvals(m)).max()) < 1.0 + 1e-3
    c = classify_planar_monodromy(m)
    assert c.regime == "complex_quartet"


def test_boundary_and_tolerance_argument() -> None:
    m = _sympl((1.0 - 2e-4, 0.2))  # pair index 2e-4 inside k = 1
    assert classify_planar_monodromy(m).regime == "boundary"
    assert classify_planar_monodromy(m).region is None
    c = classify_planar_monodromy(m, boundary_tol=1e-5)
    assert (c.regime, c.region) == ("stable", 1)


def test_six_by_six_input_takes_the_planar_block() -> None:
    m4 = _sympl((0.3, -0.6))
    m6 = np.eye(6)
    m6[np.ix_((0, 1, 3, 4), (0, 1, 3, 4))] = m4
    assert np.array_equal(planar_block(m6), m4)
    assert classify_planar_monodromy(m6).region == 1


def test_hadjidemetriou_definitions_alpha_beta_delta() -> None:
    """alpha = -trace, beta = sum of principal 2x2 minors, b = (alpha +- sqrt(Delta)) / 2."""
    c = classify_invariants(-3.0, 4.24)
    assert c.delta == pytest.approx(9.0 - 4.0 * 2.24)
    assert c.b1.real + c.b2.real == pytest.approx(-3.0)
    assert c.b1.real * c.b2.real == pytest.approx(4.24 - 2.0)  # beta = 2 + b1 b2


# --------------------------------------------------------------------------------------------
# Jorba & Olle 2004, map Ts, K = -1: exactly solvable Hamiltonian-Hopf
# --------------------------------------------------------------------------------------------


def _ts_jacobian(k: float, ell: float) -> NDArray[np.float64]:
    """Jacobian at the origin of the printed Ts (K1 = K, K2 = 0, L1 = L, L2 = -L)."""
    return np.array(
        [
            [1.0 + k + ell, k + ell, ell, ell],
            [1.0, 1.0, 0.0, 0.0],
            [-ell, -ell, 1.0 - ell, -ell],
            [0.0, 0.0, 1.0, 1.0],
        ]
    )


@pytest.mark.parametrize("ell", [0.2, 0.24, 0.249, 0.26, 0.3])
def test_ts_printed_alpha_beta_delta(ell: float) -> None:
    c = classify_planar_monodromy(_ts_jacobian(-1.0, ell))
    assert c.alpha == pytest.approx(-3.0, abs=1e-12)
    assert c.beta == pytest.approx(4.0 + ell, abs=1e-12)
    assert c.delta == pytest.approx(1.0 - 4.0 * ell, abs=1e-12)
    assert (c.regime, c.region) == (("stable", 1) if ell < 0.25 else ("complex_quartet", 2))


@pytest.mark.parametrize(
    ("ell", "angles"),
    [
        (0.24, (0.64350110879328, 0.79539883018414)),
        (0.245, (0.66752639710877, 0.77468035122454)),
        (0.249, (0.69849419144132, 0.74632549050620)),
    ],
)
def test_ts_table_1_rotation_angles(ell: float, angles: tuple[float, float]) -> None:
    """Printed Table 1 (p.697): the multiplier angles are arccos(k) of the two pair indices."""
    c = classify_planar_monodromy(_ts_jacobian(-1.0, ell))
    got = sorted(math.acos(k.real) for k in (c.k1, c.k2))
    assert got == pytest.approx(sorted(angles), abs=1e-12)


def test_ts_critical_point_is_tagged_hopf_with_direction_undetermined() -> None:
    a = _ts_jacobian(-1.0, 0.25)
    c = classify_planar_monodromy(a)
    assert c.regime == "critical"
    assert c.region is None
    assert "undetermined" in c.note
    assert math.acos(c.k1.real) == pytest.approx(math.acos(0.75), abs=1e-7)  # 0.72273424781342
    assert math.acos(0.75) == pytest.approx(0.72273424781342, abs=1e-14)
    # a single Jordan block: A - lambda I has exactly one null direction
    lam = complex(0.75, math.sqrt(1.0 - 0.75**2))
    sv = np.linalg.svd(a - lam * np.eye(4), compute_uv=False)
    assert sv[-1] < 1e-7
    assert sv[-2] > 0.1


def test_ts_family_scan_locates_delta_sign_change_at_one_quarter() -> None:
    ells = [0.2, 0.24, 0.26, 0.3]
    mons = [_ts_jacobian(-1.0, x) for x in ells]
    found = find_family_transitions(
        mons, ells, monodromy_at=lambda x: _ts_jacobian(-1.0, x), bisect_tol=1e-13
    )
    assert [(t.kind, t.index) for t in found] == [("delta", 1)]
    assert found[0].located == pytest.approx(0.25, abs=1e-12)
    assert (found[0].region_before, found[0].region_after) == (1, 2)


# --------------------------------------------------------------------------------------------
# Family scan on synthetic pair-index sweeps
# --------------------------------------------------------------------------------------------


@pytest.mark.parametrize("target", [1.0, -1.0])
def test_family_scan_finds_pair_crossing_k_plus_minus_one(target: float) -> None:
    ts = [0.0, 0.5, 1.0, 1.5, 2.0]
    t_cross = 1.2

    def at(t: float) -> NDArray[np.float64]:
        k1 = target * (0.8 + 0.2 * t / t_cross)  # passes |k| = 1 at t_cross
        return _sympl((k1, 0.1 * target), seed=3)

    found = find_family_transitions([at(t) for t in ts], ts, monodromy_at=at, bisect_tol=1e-10)
    kinds = [(f.kind, f.index) for f in found]
    assert kinds == [("k=+1" if target > 0 else "k=-1", 2)]
    assert found[0].located == pytest.approx(t_cross, abs=1e-8)


def test_family_scan_without_callback_and_with_no_change() -> None:
    mons = [_sympl((0.1 * i, -0.2)) for i in range(5)]
    assert find_family_transitions(mons) == []
    f = find_family_transitions([_sympl((0.9, 0.1)), _sympl((1.1, 0.1))])
    assert [(x.kind, x.index, x.located) for x in f] == [("k=+1", 0, None)]


# --------------------------------------------------------------------------------------------
# Broucke 1969 printed regions on re-corrected orbits (core.er3bp)
# --------------------------------------------------------------------------------------------

# (family, mu, e, x0, ydot0 printed, apoapsis start, expected region, boundary_tol)
_BROUCKE = [
    ("7P", 0.012155, 0.21, 0.0999002, 3.6218855, False, 6, 1e-3),
    ("7P", 0.012155, 0.026, 0.1459117, 3.1977458, False, 6, 5e-5),
    ("7A", 0.012155, 0.055, 0.1650050, 3.0963007, True, 1, 5e-5),
    ("7A", 0.012155, 0.35, 0.2278010, 3.0210246, True, 1, 5e-5),
    ("7A", 0.012155, 0.75, 0.2784251, 4.2150985, True, 1, 5e-5),
    ("8P", 0.5, 0.2, -0.4163137, 3.1469273, False, 6, 1e-3),
    ("8P", 0.5, 0.5, -0.4353950, 3.2408923, False, 6, 1e-3),
    ("8A", 0.5, 0.1, -0.3940499, 3.1666771, True, 4, 1e-3),
    ("8A", 0.5, 0.3, -0.3780669, 3.2891093, True, 2, 1e-3),
    ("8A", 0.5, 0.5, -0.3630355, 3.6001740, True, 2, 1e-3),
    ("8A", 0.5, 0.8, -0.3585912, 5.4899734, True, 2, 1e-3),
    ("8A", 0.5, 0.85, -0.3657577, 6.5214243, True, 1, 1e-3),
    ("11P", 0.5, 0.3, -0.1240117, 0.9964282, False, 4, 1e-3),
    ("11A", 0.5, 0.1, -0.0590284, 0.7879774, True, 6, 1e-3),
    ("11A", 0.5, 0.25, -0.0444063, 0.7358044, True, 6, 1e-3),
    ("11A", 0.5, 0.4, -0.0325438, 0.6923672, True, 3, 1e-3),
    ("11A", 0.5, 0.8, -0.0081638, 0.6195207, True, 3, 1e-3),
]


def _flow(
    sy: ER3BPSystem, f0: float, x: float, yd: float, span: float
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    _f, hist, phi = propagate_er3bp(
        np.array([x, 0.0, 0.0, 0.0, yd, 0.0]),
        (f0, f0 + span),
        sy,
        rtol=1e-13,
        atol=1e-13,
        with_stm=True,
    )
    return hist[:, -1], phi


@pytest.mark.parametrize(
    ("fam", "mu", "e", "x0", "yd0", "apo", "region", "btol"),
    _BROUCKE,
    ids=[f"{r[0]}-e{r[2]}" for r in _BROUCKE],
)
def test_broucke_printed_region_on_recorrected_orbit(
    fam: str, mu: float, e: float, x0: float, yd0: float, apo: bool, region: int, btol: float
) -> None:
    sy = ER3BPSystem(mu, e, "m1", "m2")
    f0 = math.pi if apo else 0.0

    def residual(u: NDArray[np.float64]) -> list[float]:
        end, _ = _flow(sy, f0, float(u[0]), float(u[1]), math.pi)
        return [float(end[1]), float(end[3])]

    u = fsolve(residual, [x0, yd0], xtol=1e-14)
    assert float(np.abs(u - [x0, yd0]).max()) < 1e-5  # printed state is a rounding of the orbit
    _end, phi = _flow(sy, f0, float(u[0]), float(u[1]), 2.0 * math.pi)
    c = classify_planar_monodromy(phi, boundary_tol=btol)
    assert c.region == region, (fam, e, c.regime, 2.0 * c.k1.real, 2.0 * c.k2.real)
    assert c.reciprocity_error < 1e-6


def test_7a_is_on_the_boundary_at_the_default_tolerance() -> None:
    """Broucke's 7A rows sit within 1e-3 (in k_B) of 2: the default tolerance reports boundary."""
    sy = ER3BPSystem(0.012155, 0.35, "m1", "m2")

    def residual(u: NDArray[np.float64]) -> list[float]:
        end, _ = _flow(sy, math.pi, float(u[0]), float(u[1]), math.pi)
        return [float(end[1]), float(end[3])]

    u = fsolve(residual, [0.2278010, 3.0210246], xtol=1e-14)
    _end, phi = _flow(sy, math.pi, float(u[0]), float(u[1]), 2.0 * math.pi)
    assert classify_planar_monodromy(phi).regime == "boundary"
