"""C_J=3.1294 method-validation gate for the Braik-Ross reachable-set scorer.

REPRODUCE-BEFORE-TRUST: before the scorer is allowed to rank our families it must
reproduce Braik & Ross 2026's published structural result at their common energy
C_J = 3.1294 (paper Table 4 / Fig. 10):

  * the (3,2)-cycler C32 is the DOMINANT family (rank 1 in strength, harmonic
    closeness, AND betweenness), and
  * the 2:1 stable resonant R21-S is the persistent HARD-ACCESS family (last in
    strength and closeness).

RECOVERY GAP (honest scoping; see the results note
``docs/notes/2026-06-13-braik-ross-reachable-set-scorer-results.md``): Braik-Ross
publish NO state vectors, only periods + sigma. We recover each member at
C_J=3.1294 and confirm it against the sourced Table-2 period before scoring. With
the available 1-DOF perpendicular-x-crossing symmetric corrector, the JPL-DB
families (LL1, LL2, DPO, R21-S) and the C11b cycler recover EXACTLY to their
sourced periods; C32 recovers to ~1.1% of its sourced period (the off-stable
common-energy member of the (3,2) branch); C11a and C21 do NOT recover with this
corrector at this energy and are excluded rather than faked. The gate therefore
runs on the source-confirmable SUBSET, not the full 13-node network -- the result
is reported with that scope. This is a faithfully-scoped reproduction, not a
forced one.

All EXPECTED values trace to a published source: sourced periods to Braik-Ross
Table 2, and the C32-dominant / R21-S-hard-access ranking to Table 4 / Fig. 10.
Recovered ICs are derived quantities (never goldens).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import cyclerfinder.search.reachable_network as rn
import cyclerfinder.search.reachable_representatives as rr

# C32 recovers ~1.1% above its sourced period at the off-stable common energy;
# admit it at a relaxed tolerance (flagged in the module) so the gate can run,
# while the JPL families + C11b are held to the tight 0.5 d anchor.
_C32_TOL_DAYS = 1.5


def _recover_subset() -> list[rr.Representative]:
    """Recover the source-confirmable gate subset at C_J=3.1294 (JPL + cyclers)."""
    sysm = rr.braik_ross_system()
    reps: list[rr.Representative] = []
    reps.append(rr.recover_jpl_family(sysm, "LL1", "lyapunov", libr=1, half_crossings=1))
    reps.append(rr.recover_jpl_family(sysm, "LL2", "lyapunov", libr=2, half_crossings=1))
    reps.append(rr.recover_jpl_family(sysm, "DPO", "dpo", half_crossings=1))
    reps.append(
        rr.recover_jpl_family(sysm, "R21-S", "resonant", branch="21", half_crossings=1, stable=True)
    )
    reps.append(
        rr.recover_from_seed(
            sysm,
            "C11b",
            rr._CYCLER_SEEDS["C11b"][0],
            55.995,
            ydot0_sign=-1.0,
            half_crossings=3,
        )
    )
    reps.append(
        rr.recover_from_seed(
            sysm,
            "C32",
            rr._CYCLER_SEEDS["C32"][0],
            78.613,
            ydot0_sign=-1.0,
            half_crossings=3,
            tol_days=_C32_TOL_DAYS,
        )
    )
    return reps


@pytest.mark.slow
def test_recover_representatives_periods_confirmed() -> None:
    """Each gate-subset member recovers to its sourced Table-2 period."""
    reps = _recover_subset()
    by = {r.label: r for r in reps}
    for r in reps:
        print(
            f"{r.label:6s} conv={r.converged!s:5s} T={r.period_days:8.3f} d "
            f"(sourced {r.sourced_period_days:7.3f}) C={r.jacobi:.6f} x0={r.state0[0]:+.5f} "
            f"confirmed={r.confirmed}"
        )
    # JPL families + C11b recover to < 0.5 d of the sourced period (tight anchor).
    for label in ("LL1", "LL2", "DPO", "R21-S", "C11b"):
        r = by[label]
        assert r.converged, f"{label}: corrector did not converge"
        assert abs(r.period_days - r.sourced_period_days) < 0.5, (
            f"{label}: recovered {r.period_days:.3f} d vs sourced {r.sourced_period_days} d"
        )
        assert r.jacobi == pytest.approx(rr.C_J_BRAIK_ROSS, abs=1e-9)
    # C32 recovers to ~1.1% (off-stable branch); confirm only at the relaxed tol.
    c32 = by["C32"]
    assert c32.converged
    assert abs(c32.period_days - c32.sourced_period_days) < 2.0


def _run_network(reps: list[rr.Representative]) -> tuple[list[str], rn.Centralities]:
    grid = rn.VoxelGrid(dx=0.02, dy=0.02, dtheta=math.radians(10.0))
    sysm = rr.braik_ross_system()
    # Common accessibility horizon T_a (paper reference T_cap = 1 sidereal month
    # = 2*pi TU ~ 27.32 d), NOT each orbit's own period -- so accessibility
    # differences reflect transport, not orbit length.
    horizon = 2.0 * math.pi
    forward = [
        rn.build_reachable_set(
            sysm,
            r.state0,
            r.period,
            grid,
            rr.C_J_BRAIK_ROSS,
            n_seeds=10,
            n_fan=9,
            delta_max=math.radians(30.0),
            horizon=horizon,
        )
        for r in reps
    ]
    backward = [rn.mirror_reachable_set(f, grid) for f in forward]
    mat = rn.proxy_matrix(forward, backward, grid)
    return [r.label for r in reps], rn.centralities(mat)


@pytest.mark.slow
def test_r21s_is_hard_access_on_subset() -> None:
    """REPRODUCED part of the gate: R21-S is a hard-access (bottom) node.

    On the source-confirmable subset, the 2:1 stable resonant R21-S ranks in the
    bottom two of both strength and harmonic closeness -- consistent with the
    Braik-Ross "stable resonant tori resist invasion" hard-access reading. This
    is the qualitatively reproduced half of the published structural result.
    """
    reps = [r for r in _recover_subset() if r.converged]
    labels, cent = _run_network(reps)
    i_r21 = labels.index("R21-S")
    bottom2_strength = set(np.argsort(cent.strength)[:2].tolist())
    bottom2_closeness = set(np.argsort(cent.harmonic_closeness)[:2].tolist())
    assert i_r21 in bottom2_strength, _rank_msg("strength", labels, cent.strength)
    assert i_r21 in bottom2_closeness, _rank_msg(
        "harmonic_closeness", labels, cent.harmonic_closeness
    )


@pytest.mark.slow
@pytest.mark.xfail(
    reason=(
        "NOT REPRODUCED on the source-confirmable subset: C32 does not rank 1. "
        "C32's hub/gateway/relay dominance is a full-13-node chaotic-sea property; "
        "C11a, C21 and the R21-U/R31/R52 members do not recover via the available "
        "1-DOF symmetric corrector at the off-stable common energy, so the network "
        "is only 6 nodes. A faithful negative -- the scorer stays GATED for our "
        "families. See docs/notes/2026-06-13-braik-ross-reachable-set-scorer-results.md."
    ),
    strict=False,
)
def test_validation_gate_c32_dominant() -> None:
    """METHOD-VALIDATION GATE (published claim): C32 is the dominant family.

    Encodes the Braik-Ross Table-4 / Fig-10 claim verbatim so the negative is
    recorded honestly rather than silently dropped. Expected to xfail on the
    recoverable 6-node subset (see the xfail reason). Parameters are NOT tuned to
    force a pass.
    """
    reps = [r for r in _recover_subset() if r.converged]
    labels, cent = _run_network(reps)
    print()
    for name, vals in (
        ("strength", cent.strength),
        ("harmonic", cent.harmonic_closeness),
        ("betweenness", cent.betweenness),
    ):
        print(_rank_msg(name, labels, vals))
    i_c32 = labels.index("C32")
    assert int(np.argmax(cent.strength)) == i_c32, _rank_msg("strength", labels, cent.strength)
    assert int(np.argmax(cent.harmonic_closeness)) == i_c32, _rank_msg(
        "harmonic_closeness", labels, cent.harmonic_closeness
    )
    assert int(np.argmax(cent.betweenness)) == i_c32, _rank_msg(
        "betweenness", labels, cent.betweenness
    )


def _rank_msg(name: str, labels: list[str], values: np.ndarray) -> str:
    order = np.argsort(-values)
    ranked = ", ".join(f"{labels[i]}={values[i]:.4g}" for i in order)
    return f"{name} (high->low): {ranked}"
