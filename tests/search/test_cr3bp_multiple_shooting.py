"""Multiple-shooting periodic-orbit corrector (task #687).

Positive control: the corrector MUST converge (and report a physically-sane
monodromy) on a genuine, sourced periodic orbit -- otherwise a "does not
converge" verdict on a real search seed is untrustworthy (project discipline:
verify a 0/N negative with a positive control first). The Arenstorf orbit
(Arenstorf 1963; Hairer-Norsett-Wanner "Solving ODEs I" p. 129, test B5) is the
same golden the single-shooting corrector's own tests use.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_multiple_shooting as ms

# Arenstorf golden (mu, x0, vy0, period) -- see module docstring.
MU = 0.012277471
X0, VY0, PERIOD = 0.994, -2.0015851063790825, 17.0652165601579625
_SYS = cr3bp.CR3BPSystem(mu=MU, primary="t", secondary="t", l_km=1.0, t_s=1.0)


def _arenstorf_nodes(n: int) -> tuple[list[np.ndarray], list[float]]:
    """Sample ``n`` patch points evenly in time around the Arenstorf orbit by
    forward propagation, plus the equal segment durations between them."""
    s0 = np.array([X0, 0.0, 0.0, 0.0, VY0, 0.0])
    h = PERIOD / n
    nodes = []
    s = s0.copy()
    for _ in range(n):
        nodes.append(s.copy())
        s = cr3bp.propagate(_SYS, s, h).state_f
    return nodes, [h] * n


def test_positive_control_already_periodic_reports_tiny_residual() -> None:
    # Nodes sampled from the real orbit are already continuous: ||F|| ~ 0
    # before any Newton step, and the corrector reports converged immediately.
    nodes, seg = _arenstorf_nodes(4)
    res = ms.correct_multiple_shooting(_SYS, nodes, seg, tol=1e-8)
    assert res.converged
    assert res.closure_residual < 1e-8
    assert res.period == float(np.sum(res.seg_times))


def test_positive_control_reconverges_after_node_perturbation() -> None:
    # Perturb one patch point off the orbit; the corrector must drive the
    # full-state continuity residual back to ~0 (a genuine convergence, not the
    # stall seen on a chaotic non-orbit).
    nodes, seg = _arenstorf_nodes(4)
    nodes[1] = nodes[1] + np.array([1e-3, -8e-4, 0.0, 5e-4, 3e-4, 0.0])
    res = ms.correct_multiple_shooting(_SYS, nodes, seg, tol=1e-9, max_iter=60)
    assert res.converged
    assert res.closure_residual < 1e-9
    # Independent Radau cross-check: the corrected node[0] re-closes over T.
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, res.period),
        res.nodes[0],
        args=(MU,),
        method="Radau",
        rtol=1e-11,
        atol=1e-11,
    )
    assert float(np.linalg.norm(sol.y[:, -1] - res.nodes[0])) < 1e-6


def test_positive_control_monodromy_has_trivial_unit_pair() -> None:
    # A genuine periodic orbit's monodromy has an eigenvalue pair near 1
    # (the along-flow direction). Confirms the STM product is assembled right.
    nodes, seg = _arenstorf_nodes(4)
    res = ms.correct_multiple_shooting(_SYS, nodes, seg, tol=1e-9)
    lam = ms.floquet_multipliers(_SYS, res)
    assert np.min(np.abs(lam - 1.0)) < 1e-3
    # Symplectic monodromy: reciprocal-pair structure => product of |lambda| ~ 1.
    assert abs(float(np.prod(np.abs(lam))) - 1.0) < 1e-2


def test_non_orbit_seed_does_not_converge() -> None:
    # An arbitrary non-periodic patch set must NOT be reported converged.
    rng = np.random.default_rng(687)
    base = np.array([0.6, 0.1, 0.0, 0.2, -0.3, 0.0])
    nodes = [base + 0.05 * rng.standard_normal(6) for _ in range(3)]
    seg = [1.0, 1.0, 1.0]
    res = ms.correct_multiple_shooting(_SYS, nodes, seg, tol=1e-10, max_iter=20)
    assert not res.converged


def test_length_mismatch_raises() -> None:
    nodes = [np.zeros(6), np.zeros(6)]
    try:
        ms.correct_multiple_shooting(_SYS, nodes, [1.0])
    except ValueError:
        return
    raise AssertionError("expected ValueError on nodes/seg_times length mismatch")
