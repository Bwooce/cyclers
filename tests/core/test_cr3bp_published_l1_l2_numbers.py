"""#896: printed Earth-Moon L1 and L2 numbers of the circular restricted three-body problem.

Part 1, L1 linear frequencies. Jorba, A., Jorba-Cusco, M. and Rosales, J. J. (2020), "The vicinity
of the Earth-Moon L1 point in the bicircular problem", Celestial Mechanics and Dynamical Astronomy
132:11, DOI 10.1007/s10569-019-9940-2, filed in the private paper corpus as
jorba-jorba-cusco-rosales-2020-vicinity-earth-moon-l1-point-bicircular-problem-cmda-132-11-doi-10.1007-s10569-019-9940-2.pdf.
Section 3.2, page 13, compares the Floquet frequencies of the bicircular L1 replacement orbit
(tested in ``test_bcr4bp_sun_sense.py``) with "the frequencies related to the equilibrium point L1
in the RTBP (2.33438585628816 and 2.2688310655411, respectively)": the planar (in-plane elliptic)
and the vertical frequency of the linearisation at L1, in units where the Earth-Moon mean motion
is 1. Table 1 (page 3) gives the mass ratio as mu = 0.012150581623433623.

Measured with the linearisation of ``core.cr3bp`` (``cr3bp_stm_eom`` with the identity as STM) at
an L1 found to machine precision:

- at the Table 1 mu the planar and vertical frequencies are 2.3044e-9 and 2.3545e-9 BELOW the
  printed values, about 1e-9 relative. This is the "about 1e-8" agreement of the digest.
- The two misses are the same mass-ratio error. Both frequencies grow with mu (d omega/d mu = 7.802
  planar, 7.974 vertical), and both misses are removed by the same change of mu, +2.9527e-10 (the
  two estimates agree to 4e-16 in mu). At mu = 1/(1 + 81.300585) = 0.012150581918706896 the model
  gives 2.33438585628800 and 2.26883106554112, which are the printed values to their last digit
  (residuals -1.6e-13 and +2e-14).
- Table 1's double is itself exactly 1/(1 + 81.300587) (equal to the last bit). INFERRED: the
  printed RTBP frequencies were computed with an Earth/Moon mass ratio of 81.300585, two units in
  the sixth decimal from the 81.300587 behind Table 1. The paper does not say this; it is the
  single value that reproduces both printed frequencies to their printed precision. It is not a
  rounding of mu, a loose L1 root, or a different quantity: the quantities are the linearisation
  frequencies, and they match to 1e-13 once the mass ratio is changed.

Asserted: both frequencies within 1e-12 at mu = 1/(1 + 81.300585); at the Table 1 mu, both
misses within 2e-11 of the miss predicted from d omega/d mu and the mu difference; control: the
project's default Earth-Moon mu (0.01215058439469525, from ``cr3bp_system("Earth", "Moon")``)
misses both by about +1.9e-8, more than eight times the Table 1 miss.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import brentq

from cyclerfinder.core.cr3bp import cr3bp_eom, cr3bp_stm_eom, cr3bp_system

FloatArray = NDArray[np.float64]

# Jorba, Jorba-Cusco & Rosales 2020, section 3.2, p13: RTBP L1 planar and vertical frequencies.
_JORBA_L1_PLANAR = 2.33438585628816
_JORBA_L1_VERTICAL = 2.2688310655411
# Table 1, p3.
_JORBA_TABLE1_MU = 0.012150581623433623
# INFERRED (module docstring): the mass ratio behind the printed frequencies.
_MU_FROM_FREQUENCIES = 1.0 / (1.0 + 81.300585)


def _l1_x(mu: float) -> float:
    def ax(x: float) -> float:
        return float(cr3bp_eom(0.0, np.array([x, 0.0, 0.0, 0.0, 0.0, 0.0]), mu)[3])

    return float(brentq(ax, 0.5, 1.0 - mu - 1e-3, xtol=1e-15, rtol=8.9e-16))


def _l1_frequencies(mu: float) -> tuple[float, float]:
    """Planar and vertical frequency of the linearisation of ``core.cr3bp`` at L1."""
    y42 = np.concatenate([[_l1_x(mu), 0.0, 0.0, 0.0, 0.0, 0.0], np.eye(6).ravel()])
    jac = cr3bp_stm_eom(0.0, y42, mu)[6:].reshape(6, 6)
    eigs = np.linalg.eigvals(jac)
    vertical = float(np.sqrt(-jac[5, 2]))
    imaginary = [abs(e.imag) for e in eigs if abs(e.real) < 1e-9 and abs(e.imag) > 0.0]
    planar = [w for w in imaginary if abs(w - vertical) > 1e-3]
    assert len(imaginary) == 4 and len(planar) == 2  # saddle x centre x centre
    return float(planar[0]), vertical


def test_jorba_2020_l1_frequencies_at_the_mass_ratio_81_300585() -> None:
    planar, vertical = _l1_frequencies(_MU_FROM_FREQUENCIES)
    assert abs(planar - _JORBA_L1_PLANAR) < 1e-12
    assert abs(vertical - _JORBA_L1_VERTICAL) < 1e-12


def test_jorba_2020_table1_mu_misses_both_frequencies_by_one_mass_ratio_error() -> None:
    assert _JORBA_TABLE1_MU == 1.0 / (1.0 + 81.300587)
    planar, vertical = _l1_frequencies(_JORBA_TABLE1_MU)
    miss_p, miss_v = planar - _JORBA_L1_PLANAR, vertical - _JORBA_L1_VERTICAL
    assert -2.4e-9 < miss_p < -2.2e-9 and -2.45e-9 < miss_v < -2.25e-9
    h = 1e-7
    hi = _l1_frequencies(_JORBA_TABLE1_MU + h)
    lo = _l1_frequencies(_JORBA_TABLE1_MU - h)
    slope_p, slope_v = (hi[0] - lo[0]) / (2 * h), (hi[1] - lo[1]) / (2 * h)
    d_mu = _JORBA_TABLE1_MU - _MU_FROM_FREQUENCIES
    assert abs(miss_p - slope_p * d_mu) < 2e-11
    assert abs(miss_v - slope_v * d_mu) < 2e-11


def test_project_default_earth_moon_mu_is_discriminated() -> None:
    mu = cr3bp_system("Earth", "Moon").mu
    planar, vertical = _l1_frequencies(mu)
    assert planar - _JORBA_L1_PLANAR > 1.5e-8
    assert vertical - _JORBA_L1_VERTICAL > 1.5e-8
