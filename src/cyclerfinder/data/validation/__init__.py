"""V0-V5 validation gauntlet modules (#306 Phase 1 onward).

This package holds the per-tier 3D gauntlet implementations. Phase 1 ships
V1 + V2 for 3D / broken-plane CR3BP periodic orbits; later phases add V3
(independent n-body), V4 (HFEM real-eph), and V5 (mission quality).

See ``docs/spec.md`` §14 for the V0-V5 ladder definitions and the
``feedback_orbit_closure_discipline`` memory for the independence rules
each gate must observe.
"""

from cyclerfinder.data.validation.v1_3d import (
    V1_FLOOR_KMS,
    V1_FLOOR_NONDIM_DEFAULT,
    V1Verdict3D,
    run_v1_3d,
)

__all__ = [
    "V1_FLOOR_KMS",
    "V1_FLOOR_NONDIM_DEFAULT",
    "V1Verdict3D",
    "run_v1_3d",
]
