"""Fold golden — the arclength continuator walks THROUGH an ER3BP turning point
that the secant continuator cannot pass.

GOLDEN CHOICE: self-consistent capability golden (the plan's sanctioned
fallback), NOT the sourced numeric fold.

Why the fallback. Martinez-Cacho, Gil Calvo, Bombardelli & Baresi 2025 ("Planar
retrograde periodic orbits in the elliptic restricted three-body problem",
Acta Astronautica 229:430-465; digest
``docs/notes/2026-06-25-digest-planar-retrograde-ERTBP-2025.md``) report an
explicit fold bifurcation in eccentricity at ``e = 0.0324`` for the Sun-Mars
2:3:1 IDS swing QSO (sec. 6.2.1 p. 441), and their Appendix B leaves the
Sun-Mars 2:3 ER3BP steady-QSO entry BLANK ("could not be closed in the ER3BP"
- the orbit tied to that fold). The published ICs, however, are 15-digit
*curvilinear* (pulsating Nechville-type) coordinates ``(rho0, 0, 0, theta'0, 0)``
at periapsis (Appendix B), in a retrograde-QSO frame whose conversion to our
Cartesian pulsating-frame corrector (Appendix A transform Eqs A.1-A.7) is
non-trivial and not cleanly machine-transcribable from the digest. Rather than
fabricate a frame-converted IC and assert an unsourced numeric fold at the
published e, this test takes the plan's explicitly-allowed self-consistent path.

What this test asserts (the SOURCED capability claim). The capability claim -
that naive/secant natural-parameter continuation FAILS at a fold (turning point
in e where de/ds changes sign) while pseudo-arclength continuation walks
THROUGH it - is documented by both Martinez-Cacho et al. 2025 (the fold is the
worked example of their generic e-continuation failure mode; they use a
pseudo-arclength method for the Hill families) and Peng-Bai-Xu 2017 ("we tested
to directly continue ... the routine failed" at turning points; fix =
pseudo-arclength). This test reproduces that capability gap on OUR OWN folding
family, asserting the relative capability, not an unsourced numeric value.

The folding family. The Broucke 1969 (TR 32-1360, Table 12, Family 7P) Earth-Moon
floor seed continues SMOOTHLY in e only to ~0.0549 (that smooth regime is the
regression golden in ``test_er3bp_arclength_continuation.py``). Driven HARDER -
toward ``e_target = 0.30`` - this family folds almost immediately: the secant
predictor cannot leave the low-e regime (it stalls at the first turning point),
while the pseudo-arclength walker tracks the family curve through the fold and
on to the target. This is precisely the secant-vs-arclength fold gap the sources
document, and the linchpin capability for the #436 re-run (whose discriminator
flipped at exactly these folds) and the #440 isolated-family hunt.

CONVENTION (Task 1): ``period_f`` passed to the continuators is the INTEGRATION
SPAN = the HALF-period (pi for a full-2pi orbit) when
``is_half_period_residual=True``; the Broucke EM seed is built with
``period_f = pi`` (matching ``test_er3bp_arclength_continuation.py``).
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.er3bp import ER3BPSystem
from cyclerfinder.genome.er3bp_continuation import (
    continue_er3bp_family_in_e_arclength,
    continue_er3bp_family_in_e_partial,
)

# Broucke 1969 TR 32-1360 Table 12, Family 7P, Earth-Moon mu=0.0121550, Orbit 1
# (the e~0.0001 member IC). Smooth in e only to ~0.0549; folds when pushed toward
# e_target=0.30 (see module docstring).
_MU_EM = 0.0121550
_BROUCKE_EM_IC = np.array([0.1520965, 0.0, 0.0, 0.0, 3.1608994, 0.0])
# INTEGRATION span = half-period (pi) under is_half_period_residual=True.
_PERIOD_F = np.pi
# Driven hard past the smooth ~0.0549 regime: the family folds well below this.
_E_TARGET = 0.30
_TOL = 1e-10


def _em_base() -> ER3BPSystem:
    return ER3BPSystem(mu=_MU_EM, e=0.0, primary_name="Earth", secondary_name="Moon")


@pytest.mark.slow
def test_arclength_walks_through_fold_secant_cannot_pass() -> None:
    """Arclength passes a turning point in e that the secant stalls at.

    Sourced capability (Martinez-Cacho et al. 2025; Peng-Bai-Xu 2017): naive
    secant natural-parameter continuation fails at a fold while pseudo-arclength
    walks through it. Demonstrated self-consistently on the Broucke EM Family-7P
    seed driven toward e=0.30 (folds early). Asserts the RELATIVE capability -
    arclength produces strictly more members, reaches strictly further in e past
    the secant's death point, and every member closes - never an unsourced
    numeric fold eccentricity.
    """
    # SECANT (non-raising variant): reports the eccentricity at which it dies.
    secant_orbits, death_e = continue_er3bp_family_in_e_partial(
        _em_base(),
        _BROUCKE_EM_IC,
        _PERIOD_F,
        _E_TARGET,
        n_steps=40,
        is_half_period_residual=True,
        tol=_TOL,
    )

    # The secant must DIE before the target (it cannot pass the fold).
    assert death_e is not None, (
        "secant unexpectedly reached e_target; need a folding family to exercise "
        "the fold-passing capability"
    )
    secant_emax = max(o.e for o in secant_orbits)
    assert secant_emax < _E_TARGET, (
        f"secant reached e={secant_emax}, should stall below target {_E_TARGET}"
    )

    # ARCLENGTH (fold-capable): walks the family curve through the turning point.
    arclength = continue_er3bp_family_in_e_arclength(
        _em_base(),
        _BROUCKE_EM_IC,
        _PERIOD_F,
        _E_TARGET,
        ds=0.005,
        max_steps=400,
        is_half_period_residual=True,
        tol=_TOL,
    )

    arc_es = [o.e for o in arclength]
    arc_emax = max(arc_es)

    # (a) Arclength produces strictly MORE members than the stalled secant.
    assert len(arclength) > len(secant_orbits), (
        f"arclength produced {len(arclength)} members vs secant "
        f"{len(secant_orbits)}; arclength must walk further through the fold"
    )

    # (b) Arclength reaches strictly FURTHER in e than the secant's death point -
    # i.e. it walks PAST the turning point the secant stalled at. (Comfortable
    # margin: arclength reaches the target; secant dies near e=0.)
    assert arc_emax > secant_emax, (
        f"arclength e_max={arc_emax} did not exceed secant e_max={secant_emax}; "
        "it must pass the turning point"
    )
    assert arc_emax > death_e, (
        f"arclength e_max={arc_emax} did not pass the secant death_e={death_e}; "
        "it must walk through (past) the fold"
    )
    # Reaches the requested target (the family is continuable to e_target once the
    # fold is handled).
    assert arc_emax == pytest.approx(_E_TARGET, abs=1e-3), (
        f"arclength stopped at e={arc_emax}, expected to reach {_E_TARGET}"
    )

    # (c) Every arclength member is a genuinely converged periodic orbit.
    for orbit in arclength:
        assert orbit.corrector_residual < _TOL, (
            f"arclength member at e={orbit.e} did not close "
            f"(residual {orbit.corrector_residual} >= tol {_TOL})"
        )
