"""V2 evidence for catalogue row ``em-cycler-21-3d-spatial-2026`` (#444).

This IS the mechanical evidence the ``_LEVEL_EVIDENCE`` registry points to for the
row's ``validation_level: V2`` claim. The row is the C21 stable 3D out-of-plane
spatial extension of the planar (2,1) cycler ``ross-rt-em-cycler-21-2025`` — a
``known-class-member`` (a computed member of the published Antoniadou & Libert
2019 spatial-resonant class; spec §16.4).

These assert OUR-computed gauntlet verdicts (same-model closure + bounded drift),
NOT a sourced golden value — which is exactly what a V-gauntlet evidence test
does (the EXPECTED side is the spec §14 floor, the measured side is our verdict).
"""

from __future__ import annotations

import numpy as np

from cyclerfinder.data.validation.v1_3d import run_v1_3d
from cyclerfinder.data.validation.v2_3d import run_v2_3d
from cyclerfinder.search.reachable_representatives import braik_ross_system

# The row's exact CR3BP identity tuple (orbit_elements.cr3bp.state_nd / period_nd).
_CANDIDATE_ID = "em-cycler-21-3d-spatial-2026"
_STATE0 = np.array([0.7440212218499672, 0.0, -0.2057098355650995, 0.0, 0.35368280201143637, 0.0])
_PERIOD_NDIM = 18.167169790651315


def test_v1_same_model_3d_closure_passes() -> None:
    """V1: full-asymmetric 6D closure + independent Radau under the 1 m/s floor."""
    sysm = braik_ross_system()
    v = run_v1_3d(_CANDIDATE_ID, _STATE0, _PERIOD_NDIM, sysm)
    assert v.passes_v1 is True
    assert v.converged_corrector is True
    assert v.converged_independent is True
    # Genuinely out-of-plane (not collapsed to the planar manifold).
    assert v.degenerate_planar is False
    # 9 orders under the spec §14 V1 floor (1e-3 km/s).
    assert v.independent_closure_kms < 1.0e-3


def test_v2_long_span_bounded_drift_passes() -> None:
    """V2: 6 consecutive cycles stay bounded well under the 50,000 km floor."""
    sysm = braik_ross_system()
    v = run_v2_3d(_CANDIDATE_ID, _STATE0, _PERIOD_NDIM, sysm, n_cycles=6)
    assert v.passes_v2 is True
    assert v.n_cycles_propagated == 6
    assert v.max_drift_kms < 50_000.0
