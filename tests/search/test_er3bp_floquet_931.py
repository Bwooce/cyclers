"""#931: ``floquet_classify`` can return "stable"; ``monodromy_eigenstructure`` is no longer fooled.

Spectra are synthetic symplectic matrices with chosen multipliers (the construction of
tests/core/test_floquet_classes_931.py), so the expected labels follow from the chosen spectra.
"""

from __future__ import annotations

import numpy as np
import pytest
from numpy.typing import NDArray
from scipy.linalg import expm

from cyclerfinder.search.er3bp_floquet import floquet_classify
from cyclerfinder.search.er3bp_periodic import monodromy_eigenstructure

_J = np.block([[np.zeros((2, 2)), np.eye(2)], [-np.eye(2), np.zeros((2, 2))]])
_IDX = (0, 1, 3, 4)


def _pair(k: float) -> NDArray[np.float64]:
    return np.array([[0.0, -1.0], [1.0, 2.0 * k]])


def _planar4(k1: float, k2: float) -> NDArray[np.float64]:
    m = np.zeros((4, 4))
    m[np.ix_((0, 2), (0, 2))] = _pair(k1)
    m[np.ix_((1, 3), (1, 3))] = _pair(k2)
    sym = np.random.default_rng(5).normal(size=(4, 4))
    s = expm(0.4 * _J @ (sym + sym.T) / 2.0)
    return s @ m @ np.linalg.inv(s)


def _quartet4(a: float, w: float) -> NDArray[np.float64]:
    big = np.array([[a, w], [-w, a]])
    return expm(np.block([[big, np.zeros((2, 2))], [np.zeros((2, 2)), -big.T]]))


def _embed(m4: NDArray[np.float64], kz: float = 0.5) -> NDArray[np.float64]:
    """Planar 6x6 monodromy: the 4x4 block plus a decoupled stable (z, z') pair."""
    m6 = np.zeros((6, 6))
    m6[np.ix_(_IDX, _IDX)] = m4
    m6[np.ix_((2, 5), (2, 5))] = _pair(kz)
    return m6


def test_stable_is_reachable() -> None:
    r = floquet_classify(_embed(_planar4(0.3, -0.6)))
    assert r.stability_tag == "stable"
    assert r.on_unit_circle


@pytest.mark.parametrize(
    ("m4", "tag"),
    [
        (_planar4(1.25, 0.2), "unstable"),  # saddle-centre
        (_planar4(1.25, -3.0), "unstable"),  # saddle-saddle
        (_quartet4(0.4, 0.9), "unstable"),  # complex quartet
        (_quartet4(5e-4, 0.9), "unstable"),  # moduli within 1e-3 of 1: a modulus test says marginal
        (_planar4(1.0 - 2e-4, 0.2), "marginal"),  # pair at k = 1 to the boundary tolerance
    ],
)
def test_planar_labels(m4: NDArray[np.float64], tag: str) -> None:
    assert floquet_classify(_embed(m4)).stability_tag == tag


def test_critical_hopf_point_is_marginal() -> None:
    ell = 0.25  # Jorba & Olle 2004 Ts, K = -1, printed critical L
    a = np.array(
        [
            [ell, -1.0 + ell, ell, ell],
            [1.0, 1.0, 0.0, 0.0],
            [-ell, -ell, 1.0 - ell, -ell],
            [0.0, 0.0, 1.0, 1.0],
        ]
    )
    assert floquet_classify(_embed(a)).stability_tag == "marginal"


def _spatial(h_diag: list[float]) -> NDArray[np.float64]:
    """Coupled symplectic 6x6: exp of J H conjugated by a random symplectic map."""
    j6 = np.block([[np.zeros((3, 3)), np.eye(3)], [-np.eye(3), np.zeros((3, 3))]])
    rng = np.random.default_rng(1)
    sym = rng.normal(size=(6, 6))
    s = expm(0.3 * j6 @ (sym + sym.T) / 2.0)
    m = expm(0.3 * j6 @ np.diag(h_diag))
    return s @ m @ np.linalg.inv(s)


def test_spatial_monodromy_uses_the_eigenvalue_rule() -> None:
    assert floquet_classify(_spatial([5.0] * 6)).stability_tag == "stable"
    assert floquet_classify(_spatial([-5.0, 5.0, 5.0, 5.0, 5.0, 5.0])).stability_tag == "unstable"


def test_eigenstructure_rejects_a_complex_quartet_as_centre() -> None:
    """Moduli exp(+-0.4) are within the old 0.5 acceptance of the circle; no centre exists."""
    with pytest.raises(ValueError):
        monodromy_eigenstructure(_quartet4(0.4, 0.9))


def test_eigenstructure_finds_the_centre_of_a_saddle_centre() -> None:
    r, w = monodromy_eigenstructure(_planar4(1.25, 0.2))
    assert r == pytest.approx(1.25 + (1.25**2 - 1.0) ** 0.5, rel=1e-6)
    assert w == pytest.approx(np.arccos(0.2), abs=1e-6)
