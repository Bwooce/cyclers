"""`#783` Earth-Moon He1 homoclinic-connection stage (Barrabes-Mondelo-Olle 2009).

Follow-up to `#780` (:mod:`cyclerfinder.search.earth_moon_resonant_families`), which
recovered the He1-family L1 Lyapunov periodic orbit itself (Casoliva 2010 Table 4's own
golden anchor) but explicitly scoped OUT the homoclinic CONNECTION built from it. This
module attempts that connection, mirroring how `#767` followed `#765` (Saturn-Titan) and
`#777`/`#781` followed `#776` (Neptune-Triton) -- but via a GENUINELY DIFFERENT algorithm
class than those three: Casoliva's own Class 2 cyclers are built via Barrabes, Mondelo &
Olle's own numerical-continuation-of-homoclinic-connections method (Barrabes, E., Mondelo,
J. M. & Olle, M., "Numerical Continuation of Homoclinic Connections of Periodic Orbits in
the RTBP," *Nonlinearity* 22:2901, 2009, DOI 10.1088/0951-7715/22-12-006, filed at
``cyclers_pdf/papers/barrabes-mondelo-olle-2009-...``), NOT this project's own established
Poincare-section Newton-shooting machinery (:mod:`cyclerfinder.search.jovian_resonant_connections`'s
``correct_connection``/ghost-guard, reused verbatim by the Saturn-Titan and Neptune-Triton
connection modules). This module implements Barrabes-Mondelo-Olle's OWN system of equations
(their Eq. (6)/(7), matching Casoliva 2010's own Eq. (20) almost verbatim -- the two papers
share an author and the formulation is identical) directly, reusing only the generic
STM/eigenpair utilities (:func:`~cyclerfinder.genome.heteroclinic_cycle._planar_floquet_pair`)
from the existing codebase, per the dispatch's own explicit instruction.

HONEST HEADLINE RESULT (read before trusting any number below): this is a REPRODUCTION
ATTEMPT of Casoliva's own published He1 connection at her own single most extreme energy
point (``h = -1.45016232260699``, the "closest approaches to Earth and moon" member of the
He1 family, which her own paper reports took 5 h 20 min of CPU time with specialized
continuation code, the LONGEST of the four CPU times she reports for the whole family). This
module's own Newton correctors (both single- and multiple-shooting) make genuine, measured
progress reducing the Barrabes-Mondelo-Olle Eq. (6) matching residual, but do NOT converge to
a certified connection within this task's own compute budget -- a CLEAN NEGATIVE on the
specific reproduction target, not a fabricated pass. See "THE CONDITIONING WALL" below for
the full quantitative diagnosis, and the module's own results note
(``docs/notes/2026-08-08-783-earth-moon-homoclinic-connection.md``) for the complete
investigation log. Per this project's own "never give up reproducing papers" discipline, real
effort was expended (see that section) before concluding this -- this is reported as an
honest negative, not a shortcut.

Barrabes-Mondelo-Olle's own manifold parametrization (their Eq. (4), Casoliva 2010's
unlabelled formula between Eq. (19)/Table 4 and Table 5)
------------------------------------------------------------------------------------------
    psi(theta, xi) = phi_{theta*T/(2*pi)}(z0)
                      + xi * (Lambda)^(-theta/(2*pi)) * DPhi_{theta*T/(2*pi)}(z0) @ v0

where ``phi_t`` is the time-``t`` CR3BP flow, ``z0``/``T`` are the periodic orbit's own IC/
period, ``Lambda``/``v0`` are the (unstable or stable) monodromy eigenvalue/eigenvector, and
``theta`` is a PHASE (radians, unrestricted range -- NOT wrapped to [0, 2*pi) a priori, since
the flow's own semigroup property extends ``phi_t`` naturally to any real ``t`` including
negative). :func:`analytic_manifold_seed` implements this EXACTLY (confirmed against
Barrabes-Mondelo-Olle 2009's own printed Eq. (4), read directly off the paper's own page,
NOT reconstructed from a garbled OCR text-layer extraction of Casoliva 2010's own equivalent
formula -- see "OCR SIGN HAZARD" below for why that distinction mattered this task).
:data:`HE1_XI0` = 1e-6 is Casoliva's own printed value ("kept fixed at a small value, e.g.
10^-6") and independently Barrabes-Mondelo-Olle's own stated practice ("we have taken |xi| of
the order of 10^-6").

Barrabes-Mondelo-Olle's own matching condition and unknowns (their Eq. (6))
------------------------------------------------------------------------------------------
Unknowns: ``h, z, T, Lambda^u, v^u, Lambda^s, v^s, theta^u, T^u, theta^s, T^s`` (plus, for the
MULTIPLE-shooting version they themselves use whenever ``T^u``/``T^s`` are large -- their own
words, "in order to avoid loss of precision" -- interior patch-point states along each leg).
Equations: energy/Poincare-section/periodicity/eigenpair conditions pinning the periodic
orbit, PLUS the connection-matching condition ``phi_{T^u}(psi^u(theta^u, xi_0)) -
phi_{T^s}(psi^s(theta^s, xi_0)) = 0`` with ``T^u > 0, T^s < 0`` (BOTH propagated FORWARD in
their OWN time convention -- ``T^s`` being NEGATIVE is what makes the stable leg's own
propagation run backward; this is Barrabes-Mondelo-Olle's OWN documented sign convention,
"Here T, T^u > 0 and T^s < 0", not a guess this module made).

This module takes the periodic orbit ``(z, T, Lambda^u, v^u, Lambda^s, v^s)`` as ALREADY KNOWN
(via `#780`'s own :func:`~cyclerfinder.search.earth_moon_resonant_families.recover_he1_lyapunov`
+ :func:`~cyclerfinder.genome.heteroclinic_cycle._planar_floquet_pair`, both already converged
to ~1e-12/~1e-10 relative precision respectively -- see :func:`build_he1_connection_seed` and
:func:`eigenvector_reproduction_check`, which CONFIRM this precision directly rather than
assume it) -- NOT re-solved jointly with the connection. This is a disclosed SIMPLIFICATION
of Barrabes-Mondelo-Olle's own fully-general system (needed for their own CONTINUATION-in-
energy use case, where the orbit itself must move together with the connection at each
continuation step): since this task targets a SINGLE energy point where the orbit is already
known to near-machine precision, solving only for ``(theta^u, T^u, theta^s, T^s)`` [+ interior
patch states, multiple-shooting version] is the natural reduction, and an early diagnostic
this task ran (adding the orbit's own ``(x0, T)`` back in as free unknowns) confirmed this
simplification is NOT the source of this module's own non-convergence -- see "THE CONDITIONING
WALL" below.

THE CONDITIONING WALL (the real finding of this task)
------------------------------------------------------------------------------------------
The He1 connection's own unstable leg spans ``T^u / T ~= 4.13`` periods of an orbit whose own
unstable eigenvalue is ``Lambda^u ~= 108.60`` per period. A perturbation therefore amplifies by
``exp(4.13 * ln(108.5966)) ~= exp(19.35) ~= 2.5e8`` over that leg alone (the stable leg's own
``|T^s|/T ~= 4.08`` periods amplifies comparably backward). For single-shooting Newton to have
a usable convergence basin, the SEED must be accurate to roughly (typical CR3BP-shooting
nonlinearity scale ~1e-3) / 2.5e8 ~= 4e-12 -- five orders of magnitude tighter than this
module's own achievable seed precision (this module's own recovered orbit closes to ~7.5e-12,
its eigenvector to ~1.9e-10 relative spread in the ``M @ v_u / v_u`` ratio test -- BOTH
already near machine precision, ruling out "the frozen orbit/eigenpair is the defect" as the
cause, confirmed directly this task before concluding this).

This module's own single-shooting corrector (:func:`correct_connection_single_shoot`),
seeded directly from Casoliva's own printed Table 4 ``(theta^u, T^u, T^s)`` (project-frame,
``theta^s = -theta^u`` per Barrabes-Mondelo-Olle's own footnote/symmetry statement), reduces
the matching residual from ``0.394`` (nondim, ~150,000 km position-scale) at the initial guess
to a PLATEAU around ``0.03-0.034`` (~13,000 km) after 30-75 Newton iterations, using an
STM-CHAINED semi-analytic Jacobian (NOT finite-difference across the long chaotic legs, which
would be pure numerical noise at this amplification -- only the SHORT, bounded-to-one-period
``theta``-sensitivity uses finite differences, chained with the exact per-leg STM for the long
propagation). This module's OWN multiple-shooting corrector
(:func:`correct_connection_multi_shoot`, splitting each leg into ``n_u``/``n_s`` equal
sub-intervals with continuity unknowns at each interior patch point, a block-structured
Newton system, minimum-norm least-squares solve per Barrabes-Mondelo-Olle's own documented
practice for their "over-determined and rank-deficient" system) reaches the SAME ``~0.03``
floor at both ``n_u=n_s=8`` and larger segment counts this task tested -- floor-invariant
across single-shooting, multiple-shooting, and solver strategy (``np.linalg.solve`` vs.
``np.linalg.lstsq``), which is itself the diagnostic signature of a genuine hard limit of
this task's own numerical stack, not a tuning artifact.

Barrabes-Mondelo-Olle's OWN paper independently confirms this territory is hard: they
describe using multiple shooting "in order to avoid loss of precision" specifically because
"the integration times T^u, T^s may become large," a SECOND-ORDER variational-equation
integrator (this module's own STM machinery is first-order only), an RKF7(8) integrator at
``1e-14`` tolerance (tighter than this module's own ``1e-12``/``1e-13``), an ADAPTIVE
segmentation criterion (``||DPhi_t(z_i)||_inf < M``, ``M`` "typically of the order of tens or
hundreds," recomputed every continuation step -- this module uses a FIXED, non-adaptive
segment count, a disclosed simplification), and QR-decomposition-with-column-pivoting
minimum-norm least squares specifically because their own system is "over-determined and rank
deficient." This module's own multiple-shooting implementation is a genuine, faithful
reduction of their method (same matching condition, same manifold parametrization, same
"multiple shooting because the legs are long" motivation, same min-norm-LS solve strategy)
but does not replicate every one of these refinements (adaptive segmentation, second-order
variationals, QR-pivoted rank handling) -- any of which could plausibly close the remaining
gap, and none of which was tractable to build and validate within this task's own compute
budget on top of everything else. Casoliva's own paper separately reports her own Fig. 8
connections (four OTHER, presumably easier He1-family energies) took 121, 196, 397, and 1619 s
of CPU time -- "the continuation procedure slows down as energy increases due to the
decreasing perigee distance" -- and states her own Fig. 9 connection (THIS module's own target,
``h=-1.450162``, "the closest approaches to Earth and moon" of "all the computed connections
of the He1 family") took 5 h 20 min SEPARATELY, i.e. this is Casoliva's own hardest reported
point, beyond (not merely the largest of) the four Fig. 8 timings -- using a FULL
predictor-corrector CONTINUATION from lower, easier energies. This module instead attempted a
COLD start directly at this hardest point, which is a materially harder ask than what
Casoliva's own procedure actually did.

Per this module's own honesty discipline (see the top of the docstring and the results
note): Table 5/6 pericenter/apocenter orbital elements and the LEO-rendezvous delta-v figure
are NOT computed from either corrector's own unconverged trajectory -- doing so on a
~13,000 km-scale-mismatched trajectory would produce numbers that superficially resemble
Table 5/6's own tightly-clustered magnitudes (e.g. Table 6's own Earth-relative apogees all
sit near 354,000-358,000 km) without being a genuine reproduction, exactly the "it closed!"
failure mode this project's own orbit-closure discipline warns against.

OCR SIGN HAZARD (Table 4's own ``V^u``, a SECOND, WORSE instance of the same hazard `#780`
already found on ``p_y``)
------------------------------------------------------------------------------------------
`#780`'s own module already found and fixed one OCR sign hazard on Table 4 (the ``p_y``
component of ``Y``, printed without its own minus sign in this project's ``pdftotext``
extractions). THIS task found a SECOND, worse instance on the SAME table's ``V^u`` row:
BOTH the ``pdftotext -raw`` AND ``pdftotext -layout`` text-layer extractions render ALL FOUR
of ``V^u``'s own components as if UNSIGNED (``0.04936117474608325, 0.5669396006141868,
0.8204920465924677, 0.05418270168269977``), because the same colon-for-minus-sign ligature
substitution `#780` already documented for Table 3 collapses EVERY negative sign the SAME way
regardless of which specific numbers happen to be negative. Reading the PDF's own page 13 (JGCD
p. 1635) DIRECTLY AS AN IMAGE (not through any text-layer extraction) resolved this
unambiguously: Casoliva's own printed ``V^u = (0.04936117474608325, -0.5669396006141868,
0.8204920465924677, -0.05418270168269977)`` -- components 2 and 4 (0-indexed 1 and 3) are
negative, components 1 and 3 positive. An earlier working attempt this task used a WRONG
sign pattern (guessed from the ``p_y`` precedent, before the direct image read) and got a
connection-residual discriminator test that was ~5x WORSE than the corrected version -- caught
and fixed before any constant was committed. :data:`HE1_VU_RAW_PXPY` below is vendored from
the direct image read, cross-checked (see :func:`eigenvector_reproduction_check`) against this
module's own independently-recovered monodromy eigenvector to ~1.9e-8 relative agreement (a
tight, genuine reproduction of the ONE piece of Table 4 that DOES reproduce cleanly). Table 5
and Table 6 (also newly vendored this task) were BOTH independently image-read (not
text-layer-extracted) for the same reason -- Table 6 in particular visibly carries the same
colon hazard on several of its own negative semimajor-axis and argument-of-pericenter entries
in the raw text layer.

Literature novelty gate
------------------------------------------------------------------------------------------
Not triggered -- like `#780`, this task is a reproduction ATTEMPT of Casoliva's own published
connection (which resulted in a clean negative on the reproduction target itself), explicitly
scoped away from anything novel (``search/literature_check.py`` not invoked; nothing here is
framed as a discovery).

Pure: math/numpy/scipy + :mod:`cyclerfinder.core.cr3bp`,
:mod:`cyclerfinder.search.earth_moon_resonant_families` (`#780`, the He1 periodic-orbit
anchor), :mod:`cyclerfinder.genome.heteroclinic_cycle` (reused ONLY for
``_planar_floquet_pair``, a generic STM/eigenpair utility -- the connection-matching
machinery itself, ``correct_connection``, is deliberately NOT used here, per the module's own
docstring above).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.earth_moon_resonant_families as emf
from cyclerfinder.genome.heteroclinic_cycle import _planar_floquet_pair

_PLANAR = [0, 1, 3, 4]

# ---------------------------------------------------------------------------
# Table 4 (2010, p. 1635) continuation-variable constants, corrected/re-vendored
# this task (module docstring, "OCR SIGN HAZARD"). theta_u/T_u/T_s already
# match #780's own HE1_H/HE1_PERIOD_GUESS/HE1_TARGET_LAMBDA_U (unchanged, reused
# from emf directly, not re-vendored here).
# ---------------------------------------------------------------------------

#: Barrabes-Mondelo-Olle 2009's own manifold-offset magnitude, cross-confirmed by Casoliva
#: 2010's own printed text ("kept fixed at a small value, e.g. 10^-6") -- both sources agree.
HE1_XI0 = 1e-6

#: Table 4's own printed V^u, IN CASOLIVA'S OWN (x, y, p_x, p_y) canonical-momentum basis,
#: image-read from PDF page 13 (JGCD p. 1635) directly -- see module docstring's "OCR SIGN
#: HAZARD". NOT yet flipped to project frame or converted to (x, y, vx, vy) -- see
#: :func:`_he1_vu_project_planar` for that.
HE1_VU_RAW_PXPY: tuple[float, float, float, float] = (
    0.04936117474608325,
    -0.5669396006141868,
    0.8204920465924677,
    -0.05418270168269977,
)

#: Table 4's own footnote: "V^s is obtained from V^u by changing the signs of its second and
#: third components" -- i.e. (x, y, p_x, p_y) components 2 and 3 (1-indexed), y and p_x.
HE1_VS_RAW_PXPY: tuple[float, float, float, float] = (
    HE1_VU_RAW_PXPY[0],
    -HE1_VU_RAW_PXPY[1],
    -HE1_VU_RAW_PXPY[2],
    HE1_VU_RAW_PXPY[3],
)

#: Table 4's own printed theta^u (a PHASE, radians -- module docstring).
HE1_THETA_U = 2.180794911374517

#: Table 4's own footnote: "theta^s = -theta^u".
HE1_THETA_S = -HE1_THETA_U

#: Table 4's own printed T^u (nondim time, POSITIVE -- forward propagation).
HE1_T_U = 27.6856812343605

#: Table 4's own printed T^s (nondim time, NEGATIVE -- module docstring's own "T, T^u > 0 and
#: T^s < 0" citation; this is Casoliva's/Barrabes-Mondelo-Olle's OWN sign convention, not a
#: guess this module made).
HE1_T_S = -27.33008060916397


def _pxpy_to_vxvy(dx: float, dy: float, dpx: float, dpy: float) -> NDArray[np.float64]:
    """Canonical-momentum planar-perturbation basis change, ``p_x = vx - y, p_y = vy + x``
    (project_andreu_pol_canonical_momentum convention, already used by `#780`) -- for a
    PERTURBATION vector (not a state), the SAME linear map applies (``dpx = dvx - dy``,
    ``dpy = dvy + dx``), so ``dvx = dpx + dy, dvy = dpy - dx``."""
    return np.array([dx, dy, dpx + dy, dpy - dx])


def _he1_vu_project_planar() -> NDArray[np.float64]:
    """Project-frame (this project's own primary-secondary convention) planar
    ``(x, y, vx, vy)`` perturbation vector for Table 4's own ``V^u``, via the SAME validated
    180-degree flip `#780` uses for state vectors (``(x,y,vx,vy) -> (-x,-y,-vx,-vy)`` -- a
    perturbation vector transforms identically under this rigid rotation, module docstring)."""
    return -_pxpy_to_vxvy(*HE1_VU_RAW_PXPY)


def _he1_vs_project_planar() -> NDArray[np.float64]:
    return -_pxpy_to_vxvy(*HE1_VS_RAW_PXPY)


#: Project-frame planar ``V^u``, ``V^s`` -- NOT unit-normalised (Table 4's own printed
#: ``||V^u||^2 - 1 = 0`` equation means the RAW printed vector should already be very close
#: to unit norm; :func:`eigenvector_reproduction_check` verifies this directly).
HE1_VU_PROJECT_PLANAR: NDArray[np.float64] = _he1_vu_project_planar()
HE1_VS_PROJECT_PLANAR: NDArray[np.float64] = _he1_vs_project_planar()


# ---------------------------------------------------------------------------
# Table 5 (2010, p. 1635): pericenters/apocenters of the He1 PERIODIC ORBIT itself
# (Moon-relative throughout), image-read from PDF page 13 -- module docstring.
# NOT reproduced from a corrected trajectory this task (see module docstring's
# "THE CONDITIONING WALL" -- these are vendored as citation/gate targets for a
# future task, not computed against here).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Table5Row:
    """One row of Casoliva 2010 Table 5 -- flight time (days, from periselene 1),
    Moon-relative two-body osculating pericenter/apocenter elements."""

    label: int
    t_flight_days: float
    r_km: float
    v_kms: float
    a_km: float
    e: float
    omega_deg: float


TABLE5_ROWS: tuple[Table5Row, ...] = (
    Table5Row(1, 0.000, 6331.184, 1.283, -50268.460, 1.12595, -0.000),
    Table5Row(2, 8.983, 260154.656, 0.381, -45741.513, 6.68750, 70.184),
    Table5Row(3, 14.582, 138425.392, 0.456, -35829.505, 4.86345, 0.000),
    Table5Row(4, 20.181, 260154.656, 0.381, -45741.513, 6.68750, -70.184),
)


# ---------------------------------------------------------------------------
# Table 6 (2010, p. 1636): pericenters/apocenters of the He1 CONNECTION itself
# (labels 1-5, 15-19 Moon-relative; labels 6-14 Earth-relative), image-read from
# PDF page 14 -- module docstring.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Table6Row:
    """One row of Casoliva 2010 Table 6 -- flight time (days, from periselene 1 of the
    CONNECTION part), two-body osculating pericenter/apocenter elements relative to
    ``body`` ("Moon" or "Earth")."""

    label: int
    body: str
    t_flight_days: float
    r_km: float
    v_kms: float
    a_km: float
    e: float
    omega_deg: float


TABLE6_ROWS: tuple[Table6Row, ...] = (
    Table6Row(1, "Moon", 0.000, 6579.507, 1.260, -49783.466, 1.13216, -0.760),
    Table6Row(2, "Moon", 8.565, 246058.581, 0.340, -64657.394, 4.80558, 70.381),
    Table6Row(3, "Moon", 13.730, 142669.952, 0.466, -32943.371, 5.33076, 3.139),
    Table6Row(4, "Moon", 19.935, 300407.691, 0.493, -23344.885, 13.86824, -67.513),
    Table6Row(5, "Moon", 29.939, 49917.757, 0.649, -21805.151, 3.28926, -54.055),
    Table6Row(6, "Earth", 34.250, 71035.877, 3.060, 214414.519, 0.66870, -62.458),
    Table6Row(7, "Earth", 39.973, 358249.556, 0.600, 213738.353, 0.67611, 41.848),
    Table6Row(8, "Earth", 45.666, 67869.457, 3.141, 211932.496, 0.67976, 147.477),
    Table6Row(9, "Earth", 51.236, 354118.546, 0.605, 211500.403, 0.67432, -106.097),
    Table6Row(10, "Earth", 56.816, 68946.158, 3.113, 212828.455, 0.67605, -0.000),
    Table6Row(11, "Earth", 62.396, 354118.549, 0.605, 211500.407, 0.67432, 106.097),
    Table6Row(12, "Earth", 67.965, 67869.462, 3.141, 211932.500, 0.67976, -147.477),
    Table6Row(13, "Earth", 73.659, 358249.560, 0.600, 213738.357, 0.67611, -41.848),
    Table6Row(14, "Earth", 79.382, 71035.882, 3.060, 214414.524, 0.66870, 62.458),
    Table6Row(15, "Moon", 83.693, 49917.761, 0.649, -21805.150, 3.28926, 54.055),
    Table6Row(16, "Moon", 93.697, 300407.693, 0.493, -23344.884, 13.86825, 67.513),
    Table6Row(17, "Moon", 99.902, 142669.952, 0.466, -32943.371, 5.33076, -3.139),
    Table6Row(18, "Moon", 105.067, 246058.580, 0.340, -64657.396, 4.80558, -70.381),
    Table6Row(19, "Moon", 113.632, 6579.507, 1.260, -49783.466, 1.13216, 0.760),
)

#: Text-stated (2010, p. 1635): "the minimum perigees (numbers 8 and 12) correspond to the
#: perigees of elliptical orbits with semimajor axis of 67,869 km and eccentricity of 0.68. To
#: insert a spacecraft at this perigee from a circular orbit around the Earth of radius
#: 67,869 km, a delta-v of 717.5 m/s would be required." Matches
#: :data:`~cyclerfinder.search.earth_moon_resonant_families.HE1_LEO_DV_MPS_2010` (already
#: vendored there) -- reused, NOT re-vendored here. Directly checkable from Table 6 row 8's
#: OWN printed (r, v) via vis-viva, independent of this module's own corrector: see
#: :func:`table6_row8_leo_dv_check`.


def table6_row8_leo_dv_check(mu_earth_km3s2: float = 398600.4415) -> tuple[float, float]:
    """Self-consistency check of Table 6 row 8's own printed ``(r, v)`` against the text's
    own stated 717.5 m/s LEO-rendezvous delta-v figure -- computed DIRECTLY from Casoliva's
    own printed table values (row 8: ``r=67869.457`` km, ``v=3.141`` km/s), NOT from this
    module's own (unconverged) corrected trajectory (module docstring, "THE CONDITIONING
    WALL"). Returns ``(dv_mps, v_circular_kms)``. Standard Earth GM (398600.4415 km^3/s^2,
    IAU/DE440-consistent to well within the precision this check needs) -- not itself
    sourced from Casoliva, who does not print the GM she used."""
    row8 = next(r for r in TABLE6_ROWS if r.label == 8 and r.body == "Earth")
    v_circ = math.sqrt(mu_earth_km3s2 / row8.r_km)
    dv_mps = (row8.v_kms - v_circ) * 1000.0
    return dv_mps, v_circ


# ---------------------------------------------------------------------------
# He1 connection seed: the ALREADY-CONVERGED periodic orbit + its monodromy
# eigenpair (module docstring: taken as known, not re-solved jointly with the
# connection -- a disclosed simplification of Barrabes-Mondelo-Olle's own fully
# general system).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class HE1ConnectionSeed:
    """The He1 Lyapunov periodic orbit + its planar monodromy saddle pair, as needed by
    :func:`analytic_manifold_seed`. ``state0`` is the full 6-vector IC; ``v_u``/``v_s`` are
    UNIT-NORMALISED planar 4-vectors (``_planar_floquet_pair``'s own convention)."""

    state0: NDArray[np.float64]
    period: float
    lam_u: float
    v_u: NDArray[np.float64]
    lam_s: float
    v_s: NDArray[np.float64]
    closure_residual: float


def build_he1_connection_seed(system: cr3bp.CR3BPSystem | None = None) -> HE1ConnectionSeed:
    """Recover the He1 periodic orbit (`#780`'s own
    :func:`~cyclerfinder.search.earth_moon_resonant_families.recover_he1_lyapunov`, ``tol=1e-12``)
    and its planar monodromy saddle pair (:func:`_planar_floquet_pair`). Also directly
    measures the orbit's own one-period closure residual (module docstring: confirmed
    ~7.5e-12 this task, i.e. already near machine precision -- NOT the source of this
    module's own non-convergence)."""
    sys_ = system if system is not None else emf.earth_moon_system()
    orbit = emf.recover_he1_lyapunov(sys_)
    if not orbit.converged:
        raise ValueError("build_he1_connection_seed: recover_he1_lyapunov did not converge")
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0])
    lam_u, v_u, lam_s, v_s = _planar_floquet_pair(sys_, state0, orbit.period)
    arc = cr3bp.propagate(sys_, state0, orbit.period, with_stm=False, rtol=1e-13, atol=1e-14)
    closure_residual = float(np.linalg.norm(arc.state_f - state0))
    return HE1ConnectionSeed(
        state0=state0,
        period=orbit.period,
        lam_u=lam_u,
        v_u=v_u,
        lam_s=lam_s,
        v_s=v_s,
        closure_residual=closure_residual,
    )


@dataclass(frozen=True)
class EigenvectorReproductionCheck:
    """Cross-checks this module's own independently-recovered monodromy eigenvector
    (``v_u`` from :func:`build_he1_connection_seed`) against Table 4's own printed ``V^u``
    (:data:`HE1_VU_PROJECT_PLANAR`, module docstring's OCR-corrected reading) -- the ONE
    piece of Table 4 this module DOES reproduce cleanly (module docstring)."""

    cos_angle: float
    rel_err_unit: float
    monodromy_eig_spread_rel: float
    passed: bool


#: Generous relative tolerance for the eigenvector-direction reproduction -- justified by the
#: observed match quality (this task: ~1.9e-8 relative on the unit-vector difference,
#: cos_angle within ~2e-8 of 1.0) -- 1e-2 is deliberately loose relative to that observed
#: precision, mirroring this project's own "generous relative to observed precision floor"
#: discipline elsewhere in the codebase.
EIGENVECTOR_REPRODUCTION_REL_TOL = 1e-2


def eigenvector_reproduction_check(
    system: cr3bp.CR3BPSystem | None = None,
    seed: HE1ConnectionSeed | None = None,
) -> EigenvectorReproductionCheck:
    """Two independent checks: (a) this module's own recovered ``v_u`` (unit-normalised)
    against Table 4's own printed ``V^u`` (:data:`HE1_VU_PROJECT_PLANAR`, also normalised
    here for the comparison -- Table 4's own printed vector should already be very close to
    unit norm per its own ``||V^u||^2-1=0`` equation); (b) the monodromy self-consistency
    ``M @ v_u / v_u`` componentwise ratio spread (should be tiny -- confirms ``v_u`` is
    genuinely an eigenvector of the RECOVERED orbit's own monodromy, independent of Table 4
    entirely)."""
    sys_ = system if system is not None else emf.earth_moon_system()
    sd = seed if seed is not None else build_he1_connection_seed(sys_)
    target = HE1_VU_PROJECT_PLANAR / np.linalg.norm(HE1_VU_PROJECT_PLANAR)
    recovered = sd.v_u / np.linalg.norm(sd.v_u)
    cos_angle = float(np.dot(target, recovered))
    rel_err_unit = float(np.linalg.norm(target - recovered))

    arc = cr3bp.propagate(sys_, sd.state0, sd.period, with_stm=True, rtol=1e-13, atol=1e-14)
    assert arc.stm is not None
    stm4 = arc.stm[np.ix_(_PLANAR, _PLANAR)]
    ratio = (stm4 @ sd.v_u) / sd.v_u
    spread_rel = float((ratio.max() - ratio.min()) / sd.lam_u)

    passed = rel_err_unit < EIGENVECTOR_REPRODUCTION_REL_TOL and spread_rel < 1e-4
    return EigenvectorReproductionCheck(
        cos_angle=cos_angle,
        rel_err_unit=rel_err_unit,
        monodromy_eig_spread_rel=spread_rel,
        passed=passed,
    )


# ---------------------------------------------------------------------------
# Barrabes-Mondelo-Olle Eq. (4) manifold parametrization.
# ---------------------------------------------------------------------------


def _v6(v4: NDArray[np.float64]) -> NDArray[np.float64]:
    return np.array([v4[0], v4[1], 0.0, v4[2], v4[3], 0.0])


def analytic_manifold_seed(
    system: cr3bp.CR3BPSystem,
    seed: HE1ConnectionSeed,
    theta: float,
    *,
    direction: str,
    xi: float = HE1_XI0,
    rtol: float = 1e-13,
    atol: float = 1e-14,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Barrabes-Mondelo-Olle's own Eq. (4) manifold parametrization ``psi(theta, xi)``
    (module docstring), for ``direction in {"unstable", "stable"}``. Returns ``(state6,
    stm6)`` -- the seed's own 6-state AND the 6x6 STM ``DPhi_{theta*T/(2*pi)}(z0)`` used to
    build it (returned so callers needing ``d(seed)/d(theta)`` via finite differences, or the
    onward-propagation STM chaining, do not have to re-propagate). ``theta`` is UNWRAPPED
    (may be negative or exceed ``2*pi`` -- the flow's own semigroup property handles this
    correctly, module docstring)."""
    if direction not in {"unstable", "stable"}:
        raise ValueError(f"direction must be 'unstable' or 'stable'; got {direction!r}")
    lam = seed.lam_u if direction == "unstable" else seed.lam_s
    v4 = seed.v_u if direction == "unstable" else seed.v_s
    tstep = theta * seed.period / (2.0 * math.pi)
    arc = cr3bp.propagate(system, seed.state0, tstep, with_stm=True, rtol=rtol, atol=atol)
    assert arc.stm is not None
    scale = xi * (lam ** (-theta / (2.0 * math.pi)))
    state6 = arc.state_f + scale * (arc.stm @ _v6(v4))
    return state6, arc.stm


def _d_seed_dtheta(
    system: cr3bp.CR3BPSystem,
    seed: HE1ConnectionSeed,
    theta: float,
    *,
    direction: str,
    xi: float = HE1_XI0,
    h: float = 1e-6,
) -> NDArray[np.float64]:
    """Central-difference ``d(psi)/d(theta)`` (planar 4-vector). SAFE as finite differences:
    ``theta`` only affects the seed through a propagation bounded to at most one orbital
    period (module docstring) -- not the long chaotic legs, where finite differences would be
    pure noise."""
    sp, _ = analytic_manifold_seed(system, seed, theta + h, direction=direction, xi=xi)
    sm, _ = analytic_manifold_seed(system, seed, theta - h, direction=direction, xi=xi)
    return ((sp - sm) / (2.0 * h))[_PLANAR]


# ---------------------------------------------------------------------------
# Single-shooting corrector: Newton on (theta_u, T_u, theta_s, T_s), STM-chained
# semi-analytic Jacobian (module docstring).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ConnectionShootResult:
    """Result of either corrector below. ``converged`` requires ``residual < tol`` --
    per this module's own honest headline result, this is EXPECTED to be ``False`` for the
    He1 connection at this task's own achievable compute budget (module docstring, "THE
    CONDITIONING WALL"). ``method`` is ``"single_shoot"`` or ``"multi_shoot"``."""

    method: str
    converged: bool
    residual: float
    initial_residual: float
    theta_u: float
    t_u: float
    theta_s: float
    t_s: float
    n_iter: int
    notes: str = ""


def _seg_propagate(
    system: cr3bp.CR3BPSystem, state4: NDArray[np.float64], dt: float, *, rtol: float, atol: float
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    state6 = np.array([state4[0], state4[1], 0.0, state4[2], state4[3], 0.0])
    arc = cr3bp.propagate(system, state6, dt, with_stm=True, rtol=rtol, atol=atol)
    assert arc.stm is not None
    stm4 = arc.stm[np.ix_(_PLANAR, _PLANAR)]
    return arc.state_f[_PLANAR], stm4


def connection_residual_single_shoot(
    system: cr3bp.CR3BPSystem,
    seed: HE1ConnectionSeed,
    theta_u: float,
    t_u: float,
    theta_s: float,
    t_s: float,
    *,
    xi: float = HE1_XI0,
    rtol: float = 1e-12,
    atol: float = 1e-13,
) -> NDArray[np.float64]:
    """Barrabes-Mondelo-Olle's own matching residual ``phi_{T^u}(psi^u) - phi_{T^s}(psi^s)``
    (planar 4-vector), evaluated by direct propagation (no STM -- use
    :func:`_single_shoot_residual_and_jac` for the Newton-ready version)."""
    seed_u6, _ = analytic_manifold_seed(system, seed, theta_u, direction="unstable", xi=xi)
    seed_s6, _ = analytic_manifold_seed(system, seed, theta_s, direction="stable", xi=xi)
    pu, _ = _seg_propagate(system, seed_u6[_PLANAR], t_u, rtol=rtol, atol=atol)
    ps, _ = _seg_propagate(system, seed_s6[_PLANAR], t_s, rtol=rtol, atol=atol)
    return pu - ps


def _single_shoot_residual_and_jac(
    system: cr3bp.CR3BPSystem,
    seed: HE1ConnectionSeed,
    x: NDArray[np.float64],
    *,
    xi: float,
    rtol: float,
    atol: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    theta_u, t_u, theta_s, t_s = x
    seed_u6, _ = analytic_manifold_seed(system, seed, theta_u, direction="unstable", xi=xi)
    seed_s6, _ = analytic_manifold_seed(system, seed, theta_s, direction="stable", xi=xi)
    pu, stm_u = _seg_propagate(system, seed_u6[_PLANAR], t_u, rtol=rtol, atol=atol)
    ps, stm_s = _seg_propagate(system, seed_s6[_PLANAR], t_s, rtol=rtol, atol=atol)
    resid = pu - ps

    dsu_dtheta = _d_seed_dtheta(system, seed, theta_u, direction="unstable", xi=xi)
    dss_dtheta = _d_seed_dtheta(system, seed, theta_s, direction="stable", xi=xi)
    f_pu = cr3bp.cr3bp_eom(0.0, np.array([pu[0], pu[1], 0.0, pu[2], pu[3], 0.0]), system.mu)[
        _PLANAR
    ]
    f_ps = cr3bp.cr3bp_eom(0.0, np.array([ps[0], ps[1], 0.0, ps[2], ps[3], 0.0]), system.mu)[
        _PLANAR
    ]

    jac = np.column_stack(
        [
            stm_u @ dsu_dtheta,
            f_pu,
            -(stm_s @ dss_dtheta),
            -f_ps,
        ]
    )
    return resid, jac


#: Newton convergence tolerance -- an aspirational target this module's own correctors do
#: NOT reach for He1 (module docstring); kept at the same order as the sibling connection
#: modules' own tolerances for direct comparability.
CONNECTION_TOL = 1e-7

#: (residual, jacobian) pair, or (None, None) if a candidate step hit a close encounter.
_ShootResid = tuple[NDArray[np.float64] | None, NDArray[np.float64] | None]


def correct_connection_single_shoot(
    system: cr3bp.CR3BPSystem,
    seed: HE1ConnectionSeed,
    *,
    theta_u0: float = HE1_THETA_U,
    t_u0: float = HE1_T_U,
    theta_s0: float = HE1_THETA_S,
    t_s0: float = HE1_T_S,
    xi: float = HE1_XI0,
    tol: float = CONNECTION_TOL,
    max_iter: int = 40,
    rtol: float = 1e-12,
    atol: float = 1e-13,
) -> ConnectionShootResult:
    """Single-shooting Newton on ``(theta_u, T_u, theta_s, T_s)`` with an STM-chained
    semi-analytic Jacobian (module docstring). Seeded, by default, directly from Casoliva's
    own printed Table 4 values (project-frame). A close-secondary-encounter during a
    candidate step (a genuine risk mid-search -- this connection legitimately passes within
    ~6,500 km of the Moon at its own converged solution) is treated as a rejected step (finer
    backtracking), never a crash."""
    x = np.array([theta_u0, t_u0, theta_s0, t_s0], dtype=np.float64)

    def _safe(xv: NDArray[np.float64]) -> _ShootResid:
        try:
            return _single_shoot_residual_and_jac(system, seed, xv, xi=xi, rtol=rtol, atol=atol)
        except RuntimeError:
            return None, None

    resid0, jac0 = _safe(x)
    if resid0 is None or jac0 is None:
        return ConnectionShootResult(
            "single_shoot",
            False,
            float("inf"),
            float("inf"),
            float(x[0]),
            float(x[1]),
            float(x[2]),
            float(x[3]),
            0,
            "initial guess hit a close encounter",
        )
    resid: NDArray[np.float64] = resid0
    jac: NDArray[np.float64] = jac0
    initial_residual = float(np.linalg.norm(resid))
    rn = initial_residual
    n_iter = 0
    for i in range(1, max_iter + 1):
        n_iter = i
        if rn < tol:
            break
        try:
            step = np.linalg.solve(jac, -resid)
        except np.linalg.LinAlgError:
            step, *_ = np.linalg.lstsq(jac, -resid, rcond=None)
        alpha = 1.0
        improved = False
        for _ in range(30):
            xt = x + alpha * step
            rt, jt = _safe(xt)
            if rt is not None and jt is not None and float(np.linalg.norm(rt)) < rn:
                x, resid, jac = xt, rt, jt
                rn = float(np.linalg.norm(rt))
                improved = True
                break
            alpha *= 0.5
        if not improved:
            break
    ok = rn < tol
    notes = "" if ok else "did not reach tol -- see module docstring's conditioning-wall diagnosis"
    return ConnectionShootResult(
        "single_shoot",
        ok,
        rn,
        initial_residual,
        float(x[0]),
        float(x[1]),
        float(x[2]),
        float(x[3]),
        n_iter,
        notes,
    )


# ---------------------------------------------------------------------------
# Multiple-shooting corrector: Barrabes-Mondelo-Olle's own remedy for long legs
# (module docstring). Fixed (non-adaptive) equal-time segmentation -- a disclosed
# simplification of their own adaptive ||DPhi_t||_inf < M criterion.
# ---------------------------------------------------------------------------


def correct_connection_multi_shoot(
    system: cr3bp.CR3BPSystem,
    seed: HE1ConnectionSeed,
    *,
    theta_u0: float = HE1_THETA_U,
    t_u0: float = HE1_T_U,
    theta_s0: float = HE1_THETA_S,
    t_s0: float = HE1_T_S,
    n_u: int = 8,
    n_s: int = 8,
    xi: float = HE1_XI0,
    tol: float = CONNECTION_TOL,
    max_iter: int = 40,
    rtol: float = 1e-12,
    atol: float = 1e-13,
) -> tuple[ConnectionShootResult, list[float]]:
    """Multiple-shooting Newton, splitting the unstable/stable legs into ``n_u``/``n_s``
    EQUAL-time sub-intervals with continuity unknowns at each interior patch point (module
    docstring's Barrabes-Mondelo-Olle Eq. (7) reduction). Uses ``np.linalg.lstsq``
    (minimum-norm least squares) throughout, per Barrabes-Mondelo-Olle's own documented
    practice for their own "over-determined and rank-deficient" system (module docstring) --
    this module's own dense system is square by construction (not literally over-determined),
    but ``lstsq`` is used uniformly anyway since it degrades gracefully under the severe
    ill-conditioning this problem exhibits (module docstring), unlike ``np.linalg.solve``.

    Returns ``(result, segment_stm_norms)`` -- the second element is the per-segment
    ``||STM||_inf`` at the FINAL iterate, a diagnostic for whether the fixed segmentation
    left any one segment far more ill-conditioned than its neighbours (module docstring's own
    "adaptive segmentation" discussion; Barrabes-Mondelo-Olle's own ``M`` bound, "typically
    tens to hundreds" -- reported here, not enforced, since this module's own segmentation is
    fixed, not adaptive)."""
    if n_u < 1 or n_s < 1:
        raise ValueError(f"n_u, n_s must be >= 1; got n_u={n_u}, n_s={n_s}")

    def build_initial_guess(
        theta_u: float, t_u: float, theta_s: float, t_s: float
    ) -> tuple[list[NDArray[np.float64]], list[NDArray[np.float64]]]:
        seed_u6, _ = analytic_manifold_seed(system, seed, theta_u, direction="unstable", xi=xi)
        seed_s6, _ = analytic_manifold_seed(system, seed, theta_s, direction="stable", xi=xi)
        dtu, dts = t_u / n_u, t_s / n_s
        pu = [seed_u6[_PLANAR]]
        st = pu[0]
        for _ in range(n_u - 1):
            st, _ = _seg_propagate(system, st, dtu, rtol=rtol, atol=atol)
            pu.append(st)
        ps = [seed_s6[_PLANAR]]
        st = ps[0]
        for _ in range(n_s - 1):
            st, _ = _seg_propagate(system, st, dts, rtol=rtol, atol=atol)
            ps.append(st)
        return pu, ps

    pu0, ps0 = build_initial_guess(theta_u0, t_u0, theta_s0, t_s0)
    x0 = np.concatenate([[theta_u0, t_u0, theta_s0, t_s0], *pu0[1:n_u], *ps0[1:n_s]])
    n = x0.size
    off_u = 4
    off_s = off_u + 4 * (n_u - 1)

    def residual_and_jac(
        xv: NDArray[np.float64],
    ) -> tuple[NDArray[np.float64], NDArray[np.float64], list[float]]:
        theta_u, t_u, theta_s, t_s = xv[:4]
        dtu, dts = t_u / n_u, t_s / n_s

        seed_u6, _ = analytic_manifold_seed(system, seed, theta_u, direction="unstable", xi=xi)
        seed_s6, _ = analytic_manifold_seed(system, seed, theta_s, direction="stable", xi=xi)
        dsu_dtheta = _d_seed_dtheta(system, seed, theta_u, direction="unstable", xi=xi)
        dss_dtheta = _d_seed_dtheta(system, seed, theta_s, direction="stable", xi=xi)

        r = np.zeros(n)
        jac = np.zeros((n, n))
        stm_norms: list[float] = []
        row = 0

        pts_u: list[NDArray[np.float64]] = [seed_u6[_PLANAR]] + [
            xv[off_u + 4 * i : off_u + 4 * i + 4] for i in range(n_u - 1)
        ]
        for i in range(n_u - 1):
            end, stm = _seg_propagate(system, pts_u[i], dtu, rtol=rtol, atol=atol)
            stm_norms.append(float(np.linalg.norm(stm, ord=np.inf)))
            r[row : row + 4] = end - pts_u[i + 1]
            jac[row : row + 4, off_u + 4 * i : off_u + 4 * i + 4] = -np.eye(4)
            if i == 0:
                jac[row : row + 4, 0] = stm @ dsu_dtheta
            else:
                jac[row : row + 4, off_u + 4 * (i - 1) : off_u + 4 * (i - 1) + 4] = stm
            f_end = cr3bp.cr3bp_eom(
                0.0, np.array([end[0], end[1], 0.0, end[2], end[3], 0.0]), system.mu
            )[_PLANAR]
            jac[row : row + 4, 1] += f_end / n_u
            row += 4

        pts_s: list[NDArray[np.float64]] = [seed_s6[_PLANAR]] + [
            xv[off_s + 4 * i : off_s + 4 * i + 4] for i in range(n_s - 1)
        ]
        for i in range(n_s - 1):
            end, stm = _seg_propagate(system, pts_s[i], dts, rtol=rtol, atol=atol)
            stm_norms.append(float(np.linalg.norm(stm, ord=np.inf)))
            r[row : row + 4] = end - pts_s[i + 1]
            jac[row : row + 4, off_s + 4 * i : off_s + 4 * i + 4] = -np.eye(4)
            if i == 0:
                jac[row : row + 4, 2] = stm @ dss_dtheta
            else:
                jac[row : row + 4, off_s + 4 * (i - 1) : off_s + 4 * (i - 1) + 4] = stm
            f_end = cr3bp.cr3bp_eom(
                0.0, np.array([end[0], end[1], 0.0, end[2], end[3], 0.0]), system.mu
            )[_PLANAR]
            jac[row : row + 4, 3] += f_end / n_s
            row += 4

        end_u, stm_u = _seg_propagate(system, pts_u[-1], dtu, rtol=rtol, atol=atol)
        end_s, stm_s = _seg_propagate(system, pts_s[-1], dts, rtol=rtol, atol=atol)
        stm_norms.append(float(np.linalg.norm(stm_u, ord=np.inf)))
        stm_norms.append(float(np.linalg.norm(stm_s, ord=np.inf)))
        r[row : row + 4] = end_u - end_s
        if n_u - 1 == 0:
            jac[row : row + 4, 0] += stm_u @ dsu_dtheta
        else:
            jac[row : row + 4, off_u + 4 * (n_u - 2) : off_u + 4 * (n_u - 2) + 4] = stm_u
        f_end_u = cr3bp.cr3bp_eom(
            0.0, np.array([end_u[0], end_u[1], 0.0, end_u[2], end_u[3], 0.0]), system.mu
        )[_PLANAR]
        jac[row : row + 4, 1] += f_end_u / n_u
        if n_s - 1 == 0:
            jac[row : row + 4, 2] -= stm_s @ dss_dtheta
        else:
            jac[row : row + 4, off_s + 4 * (n_s - 2) : off_s + 4 * (n_s - 2) + 4] = -stm_s
        f_end_s = cr3bp.cr3bp_eom(
            0.0, np.array([end_s[0], end_s[1], 0.0, end_s[2], end_s[3], 0.0]), system.mu
        )[_PLANAR]
        jac[row : row + 4, 3] -= f_end_s / n_s
        return r, jac, stm_norms

    def _safe(
        xv: NDArray[np.float64],
    ) -> tuple[NDArray[np.float64] | None, NDArray[np.float64] | None, list[float]]:
        try:
            return residual_and_jac(xv)
        except RuntimeError:
            return None, None, []

    x = x0
    r0, jac0, stm_norms = _safe(x)
    if r0 is None or jac0 is None:
        return (
            ConnectionShootResult(
                "multi_shoot",
                False,
                float("inf"),
                float("inf"),
                theta_u0,
                t_u0,
                theta_s0,
                t_s0,
                0,
                "initial guess hit a close encounter",
            ),
            [],
        )
    r: NDArray[np.float64] = r0
    jac: NDArray[np.float64] = jac0
    initial_residual = float(np.linalg.norm(r))
    rn = initial_residual
    n_iter = 0
    for i in range(1, max_iter + 1):
        n_iter = i
        if rn < tol:
            break
        step, *_ = np.linalg.lstsq(jac, -r, rcond=None)
        alpha = 1.0
        improved = False
        for _ in range(30):
            xt = x + alpha * step
            rt, jt, snt = _safe(xt)
            if rt is not None and jt is not None and float(np.linalg.norm(rt)) < rn:
                x, r, jac, stm_norms = xt, rt, jt, snt
                rn = float(np.linalg.norm(rt))
                improved = True
                break
            alpha *= 0.5
        if not improved:
            break
    ok = rn < tol
    notes = "" if ok else "did not reach tol -- see module docstring's conditioning-wall diagnosis"
    result = ConnectionShootResult(
        "multi_shoot",
        ok,
        rn,
        initial_residual,
        float(x[0]),
        float(x[1]),
        float(x[2]),
        float(x[3]),
        n_iter,
        notes,
    )
    return result, stm_norms


__all__ = [
    "CONNECTION_TOL",
    "EIGENVECTOR_REPRODUCTION_REL_TOL",
    "HE1_THETA_S",
    "HE1_THETA_U",
    "HE1_T_S",
    "HE1_T_U",
    "HE1_VS_PROJECT_PLANAR",
    "HE1_VS_RAW_PXPY",
    "HE1_VU_PROJECT_PLANAR",
    "HE1_VU_RAW_PXPY",
    "HE1_XI0",
    "TABLE5_ROWS",
    "TABLE6_ROWS",
    "ConnectionShootResult",
    "EigenvectorReproductionCheck",
    "HE1ConnectionSeed",
    "Table5Row",
    "Table6Row",
    "analytic_manifold_seed",
    "build_he1_connection_seed",
    "connection_residual_single_shoot",
    "correct_connection_multi_shoot",
    "correct_connection_single_shoot",
    "eigenvector_reproduction_check",
    "table6_row8_leo_dv_check",
]
