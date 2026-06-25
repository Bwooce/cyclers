"""V0-V5 validation gauntlet modules (#306 Phase 1 onward).

This package holds the per-tier 3D gauntlet implementations. Phase 1 ships
V1 + V2 for 3D / broken-plane CR3BP periodic orbits; later phases add V3
(independent n-body), V4 (HFEM real-eph), and V5 (mission quality).

See ``docs/spec.md`` §14 for the V0-V5 ladder definitions and the
``feedback_orbit_closure_discipline`` memory for the independence rules
each gate must observe.
"""

from cyclerfinder.data.validation.v0_bcr4bp import (
    V0_BCR4BP_CORRECTOR_FLOOR,
    V0_BCR4BP_PHASE_DRIFT_CONVENTION,
    V0VerdictBCR4BP,
    run_v0_bcr4bp,
)
from cyclerfinder.data.validation.v1_3d import (
    V1_FLOOR_KMS,
    V1_FLOOR_NONDIM_DEFAULT,
    V1Verdict3D,
    run_v1_3d,
)
from cyclerfinder.data.validation.v1_bcr4bp import (
    V1_BCR4BP_FLOOR_KMS,
    V1_BCR4BP_FLOOR_NONDIM,
    V1VerdictBCR4BP,
    run_v1_bcr4bp,
)
from cyclerfinder.data.validation.v2_3d import (
    V2_DRIFT_FLOOR_KMS,
    V2_N_CYCLES_MIN,
    V2Verdict3D,
    run_v2_3d,
)
from cyclerfinder.data.validation.v2_moontour import (
    V2_MOONTOUR_CLOSURE_FLOOR_KMS,
    V2_MOONTOUR_DRIFT_FLOOR_KMS,
    V2_MOONTOUR_N_CYCLES_MIN,
    MoontourCycleVerdict,
    V2MoontourVerdict,
    run_v2_moontour,
)
from cyclerfinder.data.validation.v3_3d import (
    V3_AGREEMENT_FLOOR_KMS,
    V3_N_CYCLES_MIN,
    V3CycleVerdict3D,
    V3PeriodicRegressionVerdict,
    V3Verdict3D,
    run_v3_3d,
    run_v3_periodic_regression,
)
from cyclerfinder.data.validation.v3_3d_periodic import (
    V3_PERIODIC_CLOSURE_FLOOR_NONDIM,
    V3_PERIODIC_DRIFT_AGREEMENT_FLOOR_KMS,
    V3_PERIODIC_N_CYCLES_MIN,
    V3PeriodicVerdict3D,
    run_v3_3d_periodic,
)
from cyclerfinder.data.validation.v4_uranus import (
    URANIAN_PERTURBER_MOONS,
    URANUS_J2,
    URANUS_R_EQ_KM,
    V4_AGREEMENT_FLOOR_KMS,
    V4_N_CYCLES_MIN,
    V4CycleVerdictUranus,
    V4UranusVerdict,
    run_v4_uranus,
)

__all__ = [
    "URANIAN_PERTURBER_MOONS",
    "URANUS_J2",
    "URANUS_R_EQ_KM",
    "V0_BCR4BP_CORRECTOR_FLOOR",
    "V0_BCR4BP_PHASE_DRIFT_CONVENTION",
    "V1_BCR4BP_FLOOR_KMS",
    "V1_BCR4BP_FLOOR_NONDIM",
    "V1_FLOOR_KMS",
    "V1_FLOOR_NONDIM_DEFAULT",
    "V2_DRIFT_FLOOR_KMS",
    "V2_MOONTOUR_CLOSURE_FLOOR_KMS",
    "V2_MOONTOUR_DRIFT_FLOOR_KMS",
    "V2_MOONTOUR_N_CYCLES_MIN",
    "V2_N_CYCLES_MIN",
    "V3_AGREEMENT_FLOOR_KMS",
    "V3_N_CYCLES_MIN",
    "V3_PERIODIC_CLOSURE_FLOOR_NONDIM",
    "V3_PERIODIC_DRIFT_AGREEMENT_FLOOR_KMS",
    "V3_PERIODIC_N_CYCLES_MIN",
    "V4_AGREEMENT_FLOOR_KMS",
    "V4_N_CYCLES_MIN",
    "MoontourCycleVerdict",
    "V0VerdictBCR4BP",
    "V1Verdict3D",
    "V1VerdictBCR4BP",
    "V2MoontourVerdict",
    "V2Verdict3D",
    "V3CycleVerdict3D",
    "V3PeriodicRegressionVerdict",
    "V3PeriodicVerdict3D",
    "V3Verdict3D",
    "V4CycleVerdictUranus",
    "V4UranusVerdict",
    "run_v0_bcr4bp",
    "run_v1_3d",
    "run_v1_bcr4bp",
    "run_v2_3d",
    "run_v2_moontour",
    "run_v3_3d",
    "run_v3_3d_periodic",
    "run_v3_periodic_regression",
    "run_v4_uranus",
]
