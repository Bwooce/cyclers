"""Task 6: falsification + independent-tool guards (prove the checks have teeth).

Two guards:

* **Negative gate** — feed deliberately-mismatched ``(a,e)`` / V∞ to the
  reproduction check (compute V∞ via ``construct_resonant_cycler`` and
  compare to a *wrong* claimed value) and assert it FAILS. Proves the
  reproduction gate is not a no-op; if it passed a bogus pairing, the
  catalogue's "we reproduced this" claim would be meaningless.
* **Independent-code cross-check** — verify ``construct_resonant_cycler``'s
  V∞ for an orbit agrees with an *algebraically independent* path (a
  vis-viva + angular-momentum / flight-path-angle decomposition) on the
  same circular-coplanar geometry. Reduces *code-bug* risk in the
  constructor, complementing the *data* cross-check elsewhere.

The constructor is the circular-coplanar (lowest) fidelity model, so the
positive reproduction tolerance is loose (~0.1 km/s vs sourced); the
independent-code check is tight (same model, different algebra).
"""

from __future__ import annotations

import math

import pytest

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS
from cyclerfinder.search.resonant_construct import construct_resonant_cycler

# Sourced Aldrin classic E-M cycler (a=1.6 AU, e=0.393): Rogers et al. 2012
# / spec anchors give V_inf ~ 6.5 km/s at Earth, ~9.7 km/s at Mars.
_ALDRIN_A_AU = 1.6
_ALDRIN_E = 0.393
_ALDRIN_VINF_E = 6.5
_ALDRIN_VINF_M = 9.7

# Coplanar-fidelity reproduction tolerance (km/s): the constructor's circular-
# coplanar model differs from the sourced eccentric/inclined value at this
# scale. NOT the cross-fidelity slop that masked the S1L1 bug — it is the
# documented, model-appropriate band for THIS fidelity.
_REPRO_TOL_KMS = 0.15


def _vinf_independent(a_au: float, e: float, body: str, mu: float = MU_SUN_KM3_S2) -> float:
    """V∞ at *body* via an algebraically independent path from the constructor.

    The constructor builds full inertial state vectors and differences them.
    This path instead decomposes by vis-viva speed + angular momentum:
    at the crossing radius ``r`` the heliocentric speed is
    ``sqrt(mu (2/r - 1/a))``; its transverse component is ``h/r`` with
    ``h = sqrt(mu p)``; the planet (circular, coplanar) contributes only a
    transverse ``sqrt(mu/r)``. V∞ is the magnitude of the
    (radial, transverse-difference) vector. Same model, independent algebra.
    """
    a = a_au * AU_KM
    r = PLANETS[body].sma_au * AU_KM
    p = a * (1.0 - e * e)
    v_sc = math.sqrt(mu * (2.0 / r - 1.0 / a))
    h = math.sqrt(mu * p)
    v_t = h / r
    v_r = math.sqrt(max(v_sc * v_sc - v_t * v_t, 0.0))
    v_circ = math.sqrt(mu / r)
    return math.hypot(v_r, v_t - v_circ)


# ---------------------------------------------------------------------------
# Positive reproduction (anchor): constructor reproduces sourced Aldrin V∞
# ---------------------------------------------------------------------------


def test_constructor_reproduces_sourced_aldrin_vinf() -> None:
    """Sanity anchor: correct (a,e) yields V∞ matching the sourced multiset."""
    rc = construct_resonant_cycler(_ALDRIN_A_AU, _ALDRIN_E, ("E", "M"))
    assert rc.vinf_kms["E"] == pytest.approx(_ALDRIN_VINF_E, abs=_REPRO_TOL_KMS)
    assert rc.vinf_kms["M"] == pytest.approx(_ALDRIN_VINF_M, abs=_REPRO_TOL_KMS)


# ---------------------------------------------------------------------------
# Negative gate: deliberately-mismatched (a,e) must NOT reproduce the V∞
# ---------------------------------------------------------------------------


def test_negative_gate_wrong_orbit_fails_reproduction() -> None:
    """A different orbit must NOT reproduce the Aldrin V∞ within tolerance.

    If this passed, the reproduction check would be a no-op (it would
    'reproduce' anything), so the catalogue's reproduction claim would be
    worthless. We use a clearly different but still Mars-reaching orbit.
    """
    rc = construct_resonant_cycler(2.05, 0.563, ("E", "M"))  # mcconaghy u0l1 geometry
    e_off = abs(rc.vinf_kms["E"] - _ALDRIN_VINF_E)
    m_off = abs(rc.vinf_kms["M"] - _ALDRIN_VINF_M)
    assert e_off > _REPRO_TOL_KMS or m_off > _REPRO_TOL_KMS, (
        f"a wrong orbit spuriously reproduced Aldrin V∞ "
        f"(E off {e_off:.3f}, M off {m_off:.3f} km/s) — gate is a no-op"
    )


def test_negative_gate_perturbed_eccentricity_fails() -> None:
    """Perturbing e well beyond round-off must move the computed V∞ out of band."""
    rc = construct_resonant_cycler(_ALDRIN_A_AU, _ALDRIN_E + 0.05, ("E", "M"))
    moved = abs(rc.vinf_kms["E"] - _ALDRIN_VINF_E) > _REPRO_TOL_KMS or (
        abs(rc.vinf_kms["M"] - _ALDRIN_VINF_M) > _REPRO_TOL_KMS
    )
    assert moved, "a materially perturbed e left V∞ inside tolerance — gate too loose"


# ---------------------------------------------------------------------------
# Independent-code cross-check: two algebras, same model, must agree tightly
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("a_au", "e"),
    [
        (1.6, 0.393),  # Aldrin classic
        (2.05, 0.563),  # mcconaghy u0l1
        (1.31, 0.275),  # niehoff visit2
    ],
)
def test_constructor_vinf_matches_independent_path(a_au: float, e: float) -> None:
    """construct_resonant_cycler V∞ agrees with the vis-viva path (code crosscheck)."""
    rc = construct_resonant_cycler(a_au, e, ("E", "M"))
    for body in ("E", "M"):
        indep = _vinf_independent(a_au, e, body)
        assert rc.vinf_kms[body] == pytest.approx(indep, rel=1e-9), (
            f"constructor V∞ at {body} ({rc.vinf_kms[body]}) disagrees with "
            f"independent vis-viva path ({indep}) — suspected constructor code bug"
        )


def test_unreachable_orbit_raises() -> None:
    """An orbit that cannot reach Mars must raise, not silently return a V∞.

    This is the constructor's own reach guard — the runtime twin of the
    static reach invariant in validate_physical_invariants.
    """
    with pytest.raises(ValueError, match="does not reach"):
        # Aphelion 1.155 AU < Mars perihelion 1.381 AU.
        construct_resonant_cycler(1.05, 0.10, ("E", "M"))
