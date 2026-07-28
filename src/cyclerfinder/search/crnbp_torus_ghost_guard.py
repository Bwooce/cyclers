"""#738: independent-integrator (Radau) ghost-guard cross-check for a converged
N=5 CRNBP quasi-periodic torus (:mod:`cyclerfinder.search.variational_crnbp_torus`).

Why this exists, and why it is NOT a copy of `#701`'s ``ghost_guard``
-----------------------------------------------------------------------
`#701`'s ``search.ccr4bp_heteroclinic_search.ghost_guard`` independently re-propagates a
refined heteroclinic CONNECTION's two manifold endpoints under Radau (vs the refinement's
own DOP853) to rule out chaos-amplified integrator sensitivity masquerading as a genuine
connection. There is no connection, no manifold branch, and no eigenvector direction here --
this object is a bare quasi-periodic 2-D torus. The closest existing evidence is
:func:`cyclerfinder.search.variational_crnbp_torus._independent_closure`: it already
propagates SAMPLE points on the converged torus forward through the FULL,
coupling-term-included :func:`cyclerfinder.core.crnbp.crnbp_eom` under DOP853 and checks the
result still lands on the Fourier-predicted torus surface -- methodologically independent of
the algebraic (pseudospectral collocation) solve, but using only ONE integrator. `#735` Sec 8
flagged that gap explicitly: no SECOND-INTEGRATOR cross-check exists for the torus itself, so
this row was written back at V0 pending exactly this module.

The adapted ghost-guard CONCEPT: run the identical short-time flow-consistency check under an
independent integrator (Radau -- implicit, different method family from DOP853's explicit
high-order Runge-Kutta) and confirm it agrees with the DOP853 result to a comfortable margin,
both against the torus's own Fourier-predicted target state (``*_closure_residual``, the same
quantity ``_independent_closure`` reports) and DIRECTLY against each other (``integrator_delta*``,
the sharper number: are the two integrators' own propagated end-states close, independent of
whether the torus itself is a good invariant surface).

The `#702` anchoring lesson, applied here
-------------------------------------------
`#702` found that `#701`'s original ``ghost_guard`` re-derived its own reference eigenvector
(a manifold direction's sign-continuity anchor) at a DIFFERENT phase than the quantity it was
supposed to be cross-checking against, silently comparing two different trajectories. There is
no eigenvector here (:func:`~cyclerfinder.search.variational_crnbp_torus.evaluate_torus_state`
is a deterministic, single-valued Fourier-series evaluation -- no sign ambiguity, no
phase-continuity hazard of that specific kind exists for this object class). The analogous
anchoring hazard for a TORUS cross-check is different but real: the (theta1, theta2) SAMPLE
POINTS, the elapsed time ``dt``, and the target state ``u_target`` this module compares Radau
against must be EXACTLY the same anchors ``_independent_closure``'s own DOP853 check used --
same RNG seed (``_DEFAULT_SEED = 0xC0FFEE``, `#720`/`#723`'s own value), same draw order, same
``dt_frac``/``n_samples``, same underlying ``result`` (``omega1``/``omega2``/``period``). Two
independently-drawn random sample sets would sample DIFFERENT points on the torus -- not a
cross-check of the same claim, just two unrelated measurements that could differ by chance
alone and tell you nothing about integrator (dis)agreement. :func:`radau_ghost_guard` reuses
``_independent_closure``'s own sampling code path with only the integrator ``method`` swapped,
guaranteeing byte-identical anchors (verified in
``tests/search/test_crnbp_torus_ghost_guard.py::test_dop853_side_reproduces_stored_closure_residual_exactly``,
which reproduces ``result.closure_residual`` bit-for-bit as a same-anchor sanity gate before
trusting any Radau delta reported alongside it).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.crnbp as crnbp
from cyclerfinder.core.satellites import SATELLITES
from cyclerfinder.search.variational_crnbp_torus import (
    CRNBPTorusVariationalResult,
    evaluate_torus_state,
)

# MUST match variational_crnbp_torus._independent_closure's own hard-coded seed -- this is the
# anchoring choice this module's whole docstring is about: reusing the SAME sample points the
# quantity being cross-checked (result.closure_residual) was itself computed from, not drawing
# an independent set that would silently turn this into two unrelated measurements.
_DEFAULT_SEED = 0xC0FFEE

# Europa's SMA, matching the catalogue row's own crnbp_provenance.torus... lunit_km /
# orbit_elements.cr3bp.lunit_km convention -- SOURCED from the same registry
# core.crnbp.jupiter_europa_io_ganymede_default() derives its own DU from (core.satellites, not
# re-hardcoded here).
_DEFAULT_LUNIT_KM = SATELLITES["Europa"].sma_km


@dataclass(frozen=True)
class CRNBPTorusGhostGuardSample:
    """One (theta1_0, theta2_0) anchor's per-integrator closure + cross-integrator delta.

    ``dop853_closure``/``radau_closure`` are each integrator's own
    ``||propagated_end - torus_predicted_target||`` in the SAME nondimensional 4-state
    (x, y, vx, vy) norm :func:`~cyclerfinder.search.variational_crnbp_torus._independent_closure`
    uses (so ``dop853_closure`` reproduces that function's own per-sample number exactly).
    ``integrator_delta``/``integrator_delta_pos_km`` compare the two integrators' propagated
    end-states DIRECTLY to each other (independent of how good the torus itself is) -- the
    sharper, `#701`-``ghost_guard``-style number.
    """

    theta1_0: float
    theta2_0: float
    dop853_closure: float
    radau_closure: float
    integrator_delta: float
    integrator_delta_pos_km: float


@dataclass(frozen=True)
class CRNBPTorusGhostGuardReport:
    """Independent-integrator cross-check report for one converged CRNBP torus."""

    n_samples: int
    dt_frac: float
    seed: int
    dop853_closure_residual: float
    radau_closure_residual: float
    integrator_delta_max: float
    integrator_delta_max_km: float
    reproduces_reported_closure: bool
    radau_consistent: bool
    genuine: bool
    notes: str
    samples: tuple[CRNBPTorusGhostGuardSample, ...] = field(default_factory=tuple)


def radau_ghost_guard(
    result: CRNBPTorusVariationalResult,
    *,
    dt_frac: float = 0.02,
    n_samples: int = 12,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    seed: int = _DEFAULT_SEED,
    lunit_km: float = _DEFAULT_LUNIT_KM,
    integrator_consistency_km: float = 1.0,
    closure_reproduction_reltol: float = 1e-9,
) -> CRNBPTorusGhostGuardReport:
    """Cross-check ``result``'s short-time flow-consistency claim under Radau.

    Propagates the SAME ``n_samples`` (theta1, theta2) anchor points
    :func:`~cyclerfinder.search.variational_crnbp_torus._independent_closure` draws (identical
    RNG seed and draw order -- see module docstring) forward under BOTH DOP853 (reproducing
    that function's own reported ``result.closure_residual``, verified below) and Radau (the
    independent second integrator), through the FULL coupling-term-included
    :func:`cyclerfinder.core.crnbp.crnbp_eom`.

    ``genuine`` is True only if (1) the DOP853 side reproduces ``result.closure_residual`` to
    within ``closure_reproduction_reltol`` relative (a same-anchor sanity check: if this fails,
    the sampling/propagation here has silently diverged from what the corrector itself computed,
    and NOTHING below can be trusted) AND (2) the direct Radau-vs-DOP853 position gap stays under
    ``integrator_consistency_km`` (default 1.0 km, matching the `#701`/`#708` Umbriel-Titania
    row's own ``ghost_guard`` gate exactly, for a directly comparable V1 bar).
    """
    if result.closure_residual is None or not np.isfinite(result.closure_residual):
        raise ValueError(
            "result.closure_residual is missing/non-finite -- run "
            "correct_crnbp_torus_pseudospectral (or equivalent) to populate it before "
            "cross-checking it here."
        )

    dt = dt_frac * result.period
    rng = np.random.default_rng(seed)
    samples: list[CRNBPTorusGhostGuardSample] = []
    max_dop853 = 0.0
    max_radau = 0.0
    max_delta = 0.0
    max_delta_km = 0.0
    for _ in range(n_samples):
        th1 = float(rng.uniform(0, 2 * np.pi))
        th2 = float(rng.uniform(0, 2 * np.pi))
        u0 = evaluate_torus_state(result, th1, th2)
        s6 = np.array([u0[0], u0[1], 0.0, u0[2], u0[3], 0.0], dtype=np.float64)
        t0 = th1 / result.omega1
        u_target = evaluate_torus_state(result, th1 + result.omega1 * dt, th2 + result.omega2 * dt)

        planar_by_method: dict[str, NDArray[np.float64]] = {}
        closure_by_method: dict[str, float] = {}
        for method in ("DOP853", "Radau"):
            sol = solve_ivp(
                crnbp.crnbp_eom,
                (t0, t0 + dt),
                s6,
                args=(result.system,),
                method=method,
                rtol=rtol,
                atol=atol,
            )
            if not sol.success:
                planar_by_method[method] = np.full(4, np.inf)
                closure_by_method[method] = float("inf")
                continue
            end = sol.y[:, -1]
            planar_end = np.array([end[0], end[1], end[3], end[4]], dtype=np.float64)
            planar_by_method[method] = planar_end
            closure_by_method[method] = float(np.linalg.norm(planar_end - u_target))

        delta = float(np.linalg.norm(planar_by_method["DOP853"] - planar_by_method["Radau"]))
        # Position-only (x, y) component of the SAME two end-states, physical km -- the direct
        # analogue of `#701`'s ghost_guard.integrator_delta_km, for a like-for-like V1 bar.
        delta_pos_km = (
            float(np.linalg.norm(planar_by_method["DOP853"][:2] - planar_by_method["Radau"][:2]))
            * lunit_km
        )

        max_dop853 = max(max_dop853, closure_by_method["DOP853"])
        max_radau = max(max_radau, closure_by_method["Radau"])
        max_delta = max(max_delta, delta)
        max_delta_km = max(max_delta_km, delta_pos_km)
        samples.append(
            CRNBPTorusGhostGuardSample(
                theta1_0=th1,
                theta2_0=th2,
                dop853_closure=closure_by_method["DOP853"],
                radau_closure=closure_by_method["Radau"],
                integrator_delta=delta,
                integrator_delta_pos_km=delta_pos_km,
            )
        )

    reported = float(result.closure_residual)
    reproduces = np.isfinite(max_dop853) and abs(
        max_dop853 - reported
    ) <= closure_reproduction_reltol * max(reported, 1e-300)
    radau_consistent = np.isfinite(max_delta_km) and max_delta_km < integrator_consistency_km

    notes_parts = []
    if not reproduces:
        notes_parts.append(
            f"DOP853 same-anchor reproduction FAILED: got {max_dop853:.6e}, "
            f"result.closure_residual is {reported:.6e} -- sampling/propagation here has "
            "diverged from _independent_closure's own anchors; nothing else in this report "
            "can be trusted until this is fixed"
        )
    if not radau_consistent:
        notes_parts.append(
            f"Radau/DOP853 disagree by {max_delta_km:.6e} km "
            f">= {integrator_consistency_km} km threshold"
        )
    genuine = reproduces and radau_consistent
    notes = "; ".join(notes_parts) if notes_parts else "all checks passed"

    return CRNBPTorusGhostGuardReport(
        n_samples=n_samples,
        dt_frac=dt_frac,
        seed=seed,
        dop853_closure_residual=max_dop853,
        radau_closure_residual=max_radau,
        integrator_delta_max=max_delta,
        integrator_delta_max_km=max_delta_km,
        reproduces_reported_closure=bool(reproduces),
        radau_consistent=bool(radau_consistent),
        genuine=bool(genuine),
        notes=notes,
        samples=tuple(samples),
    )
