"""Independent-closure gate on the ER3BP arclength continuator (#441).

The #441 (Phase 2) bridge spike found the ``period_f`` trap: when the fixed
integration span is NOT a multiple of pi, the symmetric corrector can report a
member as CONVERGED (corrector residual at machine zero — it only zeroes
``(y, xdot)`` at the crossing) while the FULL orbit does not close. Only an
independent full-period re-propagation (``independent_residual``) catches this.
Continuing on such a span manufactures a family of false-positive "members".

The durable fix the spike prescribes: *gate every member on
``independent_residual``, never ``corrector_residual`` alone.* This test pins
the gate MECHANISM with numbers measured directly from the Broucke EM Family-7P
folding family (commensurate ``period_f = pi``):

* a production gate (1e-8) preserves the whole valid family (its members close
  to ~1e-10 — comfortably under 1e-8), so the gate never harms good families;
* an artificially tight gate (2e-10), set *inside* the family's natural
  independent-residual range, truncates the walk and every returned member
  satisfies the gate — proving the mechanism actually rejects members the
  ungated walk would admit.

The FULL convention reproduction (continue a #440 resonant seed on
``period_f = T/2`` and show the independent residual DIVERGES while a
commensurate span stays tiny) needs the #440 Phase-1 seed machinery and is
scheduled in the #441 Phase 2 plan, not here.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.er3bp import ER3BPSystem
from cyclerfinder.genome.er3bp_continuation import continue_er3bp_family_in_e_arclength

# Broucke 1969 TR 32-1360 Table 12 Family 7P Earth-Moon seed (same as the fold
# golden). Commensurate half-period span pi => full period 2pi => closes.
_MU_EM = 0.0121550
_BROUCKE_EM_IC = np.array([0.1520965, 0.0, 0.0, 0.0, 3.1608994, 0.0])
_PERIOD_F = np.pi
_E_TARGET = 0.20
_TOL = 1e-10

# Production gate: well above the valid family's ~1e-10 independent residual.
_PROD_GATE = 1e-8
# Tight gate set INSIDE the family's measured independent-residual range
# (members run ~9.5e-11 near e=0 up to ~3.6e-10 by e~0.18), so it must truncate.
_TIGHT_GATE = 2e-10


def _em_base() -> ER3BPSystem:
    return ER3BPSystem(mu=_MU_EM, e=0.0, primary_name="Earth", secondary_name="Moon")


@pytest.mark.slow
def test_production_gate_preserves_valid_family() -> None:
    """A 1e-8 independent gate leaves the valid commensurate family intact.

    The Broucke EM family closes to ~1e-10 in independent residual; a 1e-8 gate
    is far above that, so gating must return exactly the ungated family (no good
    member dropped).
    """
    ungated = continue_er3bp_family_in_e_arclength(
        _em_base(),
        _BROUCKE_EM_IC,
        _PERIOD_F,
        _E_TARGET,
        ds=0.01,
        max_steps=60,
        is_half_period_residual=True,
        tol=_TOL,
    )
    gated = continue_er3bp_family_in_e_arclength(
        _em_base(),
        _BROUCKE_EM_IC,
        _PERIOD_F,
        _E_TARGET,
        ds=0.01,
        max_steps=60,
        is_half_period_residual=True,
        tol=_TOL,
        independent_gate=_PROD_GATE,
    )

    # The valid family really does close well under the production gate.
    assert max(o.independent_residual for o in ungated) < _PROD_GATE, (
        "test premise broken: valid family exceeds the production gate"
    )
    # No good member is dropped.
    assert len(gated) == len(ungated)
    assert gated[-1].e == pytest.approx(ungated[-1].e, abs=1e-9)


@pytest.mark.slow
def test_tight_gate_truncates_and_every_member_satisfies_gate() -> None:
    """A gate set inside the family's residual range truncates the walk.

    Every returned member must satisfy ``independent_residual <= gate``, and the
    ungated walk must admit at least one member the gate rejects — proving the
    mechanism excludes non-closing members rather than being a no-op.
    """
    ungated = continue_er3bp_family_in_e_arclength(
        _em_base(),
        _BROUCKE_EM_IC,
        _PERIOD_F,
        _E_TARGET,
        ds=0.01,
        max_steps=60,
        is_half_period_residual=True,
        tol=_TOL,
    )
    gated = continue_er3bp_family_in_e_arclength(
        _em_base(),
        _BROUCKE_EM_IC,
        _PERIOD_F,
        _E_TARGET,
        ds=0.01,
        max_steps=60,
        is_half_period_residual=True,
        tol=_TOL,
        independent_gate=_TIGHT_GATE,
    )

    # The ungated walk admits members the tight gate rejects.
    assert max(o.independent_residual for o in ungated) > _TIGHT_GATE, (
        "test premise broken: family never exceeds the tight gate"
    )
    # Every gated member honours the gate.
    assert gated, "gate dropped the seed itself; seed must close at e=0"
    for o in gated:
        assert o.independent_residual <= _TIGHT_GATE, (
            f"gated member at e={o.e} has independent_residual "
            f"{o.independent_residual} > gate {_TIGHT_GATE}"
        )
    # And the walk is genuinely truncated relative to the ungated family.
    assert gated[-1].e < ungated[-1].e, (
        f"gated walk reached e={gated[-1].e}, not truncated below ungated e={ungated[-1].e}"
    )
