"""Task #808: wrong-topology clean-negative gate on the grid paths -- tests.

#807 found (and fixed, in `pluto_charon_kk_sweep`) a latent failure class:
a (k1,k2)-family sweep whose continuation/C-sweep loses the target branch
and instead converges a stable orbit of a DIFFERENT family was reported as
`stable_found=True, topology_ok=False` -- a confusing partial-failure state
-- instead of a clean negative for the target family. #808 closes the same
documented gap (#660 scope decision) in `real_binary_kk_sweep`'s two grid
paths, which verify only the SEED's topology and then run un-re-checked
continuations:

* `_finalize_grid_candidate` (gravity-only, used by `sweep_family_grid`),
* `sweep_family_grid_srp` (SRP: gravity-only seed, then beta-continuation
  AND a C_srp-sweep, neither re-checked) via `_topology_gated_result_srp`.

This is not hypothetical for the grid machinery: #656's PC (5,1) grid sweep
measured exactly this event (genuine (5,1) seed found by `_grid_seed_search`,
`c_sweep_find_nu_zero` with hc=None then walked onto the retrograde (4,0)
family and reported stable_found=True/topology_ok=False -- hand-diagnosed at
the time, see #656's bullet in data/OUTSTANDING.md).

Test strategy (mirrors test_660's own cheap-reconvergence precedent): rather
than re-running a ~75s discovery grid search, reconverge the known Antiope
(2,2) orbit from #659's recorded IC near-instantly, then feed it to the
finalizers under a deliberately MISLABELED target (5,1) -- a controlled,
deterministic stand-in for a C-sweep that drifted onto a different family.
The gate must convert it to a clean negative recording the RECOVERED (2,2)
topology. The correctly-labeled (2,2) path (topology gate passes, #660
clearance gate then rejects) is already pinned by
test_660_antiope_22_fails_clearance_gate_explicitly, which doubles as the
proof that #808's gate does not swallow topology-correct candidates.
"""

from __future__ import annotations

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp
from cyclerfinder.search.real_binary_kk_sweep import (
    DEFAULT_CLEARANCE_MARGIN_KM,
    REAL_BINARY_SYSTEMS,
    _finalize_grid_candidate,
    _topology_gated_result_srp,
)

# #659's recorded Antiope (2,2) converged IC (docs/notes/
# 2026-07-19-659-antiope-adjudication.md) -- same numbers test_660 reuses.
_ANTIOPE_22_X0 = -0.5742744462570041
_ANTIOPE_22_C = 3.4661023165370235
_ANTIOPE_22_T = 6.011499192614617


def _reconverge_antiope_22() -> tuple[cr3bp.CR3BPSystem, cp.SymmetricOrbit]:
    antiope = REAL_BINARY_SYSTEMS["antiope"]
    target = antiope.to_cr3bp_system()
    orbit = cp.correct_symmetric_fixed_jacobi(
        target,
        _ANTIOPE_22_X0,
        _ANTIOPE_22_C,
        _ANTIOPE_22_T,
        ydot0_sign=-1.0,
        half_crossings=None,
        tol=1e-10,
    )
    assert orbit.converged, "reconvergence from the known #659 (2,2) IC must succeed"
    return target, orbit


def test_808_grid_finalizer_gates_wrong_topology_to_clean_negative() -> None:
    """A stable orbit of a different family, labeled (5,1), must come out of
    `_finalize_grid_candidate` as a clean negative for (5,1) -- with the
    recovered (2,2) topology recorded in `note` -- NOT as the pre-#808
    confusing `stable_found=True, topology_ok=False` state."""
    target, orbit = _reconverge_antiope_22()
    antiope = REAL_BINARY_SYSTEMS["antiope"]

    result = _finalize_grid_candidate(
        5,
        1,
        target,
        orbit,
        "mislabeled_topology_gate_check",
        radius_km_primary=antiope.radius_km_primary,
        radius_km_secondary=antiope.radius_km_secondary,
        clearance_margin_km=DEFAULT_CLEARANCE_MARGIN_KM,
    )

    assert result.stable_found is False, (
        f"wrong-topology grid result must be a clean negative (note={result.note!r})"
    )
    assert result.topology_ok is False
    assert "wrong topology" in result.note
    assert "recovered (k1,k2)=(2,2)" in result.note, (
        f"note must record the RECOVERED topology, got note={result.note!r}"
    )
    # The gate runs BEFORE the #660 clearance gate: no clearance evaluation
    # should have happened on a wrong-family orbit.
    assert result.min_clearance_ok is None
    assert result.min_distance_primary_km is None
    assert result.min_distance_secondary_km is None


def test_808_srp_grid_gate_wrong_topology_clean_negative_at_beta_zero() -> None:
    """Same gate on the SRP grid path via `_topology_gated_result_srp`,
    exercised at beta_nd=0 (exact reduction to gravity-only dynamics, proven
    in test_665_real_binary_srp) so the known Antiope orbit stays a genuine
    periodic orbit of the augmented model."""
    target, orbit = _reconverge_antiope_22()

    result = _topology_gated_result_srp(
        5,
        1,
        target,
        orbit,
        0.0,  # beta_nd = 0: SRP-augmented model reduces exactly to gravity-only
        0.0,  # phi0 (on-axis, the only corrector-admissible choice per #665)
        "mislabeled_srp_topology_gate_check",
    )

    assert result.stable_found is False
    assert result.topology_ok is False
    assert "wrong topology" in result.note
    assert "recovered (k1,k2)=(2,2)" in result.note
    assert result.beta_nd == 0.0
    assert result.phi0 == 0.0


def test_808_srp_gate_passes_through_correct_topology() -> None:
    """Correctly-labeled (2,2) must pass the SRP topology gate untouched
    (topology_ok True, real diagnostics populated, no gate note) -- the gate
    only changes how a MISMATCH is reported. (The gravity-only equivalent is
    already pinned by test_660_antiope_22_fails_clearance_gate_explicitly,
    whose clearance-rejection numbers prove `_finalize_grid_candidate` still
    reaches the #660 gate for topology-correct candidates.)"""
    target, orbit = _reconverge_antiope_22()

    result = _topology_gated_result_srp(
        2,
        2,
        target,
        orbit,
        0.0,
        0.0,
        "correctly_labeled_srp_gate_check",
    )

    assert result.topology_ok is True
    assert result.stable_found is True  # #659: genuinely stable (fails only clearance)
    assert "wrong topology" not in result.note
    assert result.jacobi_mid is not None
    assert abs(result.jacobi_mid - _ANTIOPE_22_C) < 1e-9
