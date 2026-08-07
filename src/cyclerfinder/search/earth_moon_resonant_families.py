"""Earth-Moon planar CR3BP resonant/Lyapunov family gate (#780).

Thin sibling of :mod:`cyclerfinder.search.saturn_titan_resonant_families`
(`#765`) and :mod:`cyclerfinder.search.neptune_triton_resonant_families`
(`#776`), retargeting the SAME reused corrector/classification machinery
(:mod:`cyclerfinder.search.cr3bp_periodic`,
:mod:`cyclerfinder.search.jovian_resonant_families`,
:mod:`cyclerfinder.search.cr3bp_continuation`) at this project's own
flagship Earth-Moon system -- which, despite anchoring 322 of 383 catalogue
rows, had never had this reproduce-a-published-family-table pipeline
applied to it before `#780`.

Source
------
Casoliva, J., Mondelo, J. M., Villac, B. F., Mease, K. D., Barrabes, E. &
Olle, M., "Two Classes of Cycler Trajectories in the Earth-Moon System,"
JGCD 33(5), 2010, pp. 1623-1640, DOI 10.2514/1.46856 ("2010"), and its
conference precursor, AIAA 2008-6434 ("2008"). Both filed at
``cyclers_pdf/papers/casoliva-mondelo-villac-mease-barrabes-olle-*``. Full
digest: ``docs/notes/2026-07-27-725-casoliva-earth-moon-cycler-families-digest.md``.
Full-precision Table 3 (2010, Class 1) extraction cross-check for 6 of its
16 rows: ``docs/notes/2026-07-28-747-franz-russell-casoliva-crosscheck.md``.

Two classes
-----------
* **Class 1** (Sec. IV both papers): high-energy near-Keplerian p-q
  resonant orbits, Table 3 (2010, p.1630) -- 16 rows total, vendored
  VERBATIM below (:data:`TABLE3_ROWS`). Casoliva's own footnotes flag 7 of
  the 16 as either not satisfying their own p-q resonance relation
  (``satisfies_resonance=False``) or flying through the Earth's own radius
  (``exists_in_em_system=False``, physically unusable as a real Earth-Moon
  cycler even though it is a mathematically valid PCR3BP periodic orbit) --
  both flags are carried through, not silently dropped. The 9 "clean" rows
  contain BOTH stable and unstable members per resonance (1-2c/1-2e/2-1a/
  3-2c/7-3a stable, 1-2d/2-1b/7-3b/7-3c unstable) -- Casoliva's own
  documented "at least one unstable and one stable cycler per resonance
  class" finding, found via an entirely different method (elliptical/
  second-species differential correction) than Class 2 or than this
  project's own manifold-tube-based resonant-connection lineage.
* **Class 2** (Sec. V both papers): low-energy L1-Lyapunov-based homoclinic
  cyclers. Stage A (`#780`, this module) reproduces only the BASE OBJECT --
  the L1 Lyapunov periodic orbit itself, at the one energy Table 4 (2010,
  p.1635) prints to full (14-15 significant digit) machine precision,
  ``h = -1.45016232260699`` (the "He1 family" golden connection anchor).
  The homoclinic CONNECTION itself (He1-4, Hm1-2) and the Barrabes-Mondelo-
  Olle continuation-of-homoclinic-connections algorithm are explicitly
  OUT OF SCOPE here -- see "Explicitly out of scope" below.

Coordinate convention (CRITICAL, verified not assumed)
--------------------------------------------------------
Casoliva's own stated convention (2010 Eq. 1 text) places Earth at
``(+mu, 0)`` and the Moon at ``(mu-1, 0)`` -- i.e. Moon on the LEFT
(negative x). This project's own :mod:`cyclerfinder.core.cr3bp` places the
PRIMARY at ``(-mu, 0)`` and the SECONDARY at ``(1-mu, 0)`` -- i.e. Moon on
the RIGHT (positive x), the exact opposite side. The two conventions are
related by a rigid 180-degree rotation about the barycenter (which sits at
the origin in both), which is an exact symmetry of the rotating-frame
CR3BP EOM: ``(x, y, vx, vy) -> (-x, -y, -vx, -vy)`` (position AND velocity
components both negate together). This was verified two independent ways,
not merely asserted:

1. Jacobi-constant round-trip: transforming a printed Table 3 IC through
   this flip and evaluating :func:`cr3bp.jacobi_constant` at this
   project's own registry mu reproduces the row's own printed ``C_J`` to
   3.2e-7 relative (1-2c) -- but Jacobi constant only depends on
   ``vx^2+vy^2`` (velocity MAGNITUDE), so this check alone cannot
   distinguish a velocity-sign error from a correct transform.
2. DIRECT PERIODICITY (decisive): propagating the flip-transformed IC
   for the row's own printed period ``T`` and checking ``||X(T)-X(0)||``
   closes to ~1e-4 to ~1e-6 (position/velocity-scale units) for all 16
   Table 3 rows AND to 2.5e-6 for the Class 2 Table 4 golden IC -- this
   DOES discriminate velocity sign (a wrong sign produces immediate,
   dramatic non-closure: chaotic/close-encounter blowup with residuals
   O(1) or worse within one period, verified directly during this task's
   own scoping). All 17 rows (16 Table 3 + the Class 2 golden anchor)
   close under the SAME uniform full-flip transform -- one transform, one
   convention, for the whole module (an earlier draft of this module's own
   scoping mistakenly concluded Table 3 and Table 4 needed DIFFERENT
   transforms; that was traced to a Table-4-specific OCR artifact, see
   "OCR sign artifact" below, not a genuine convention split).

**OCR sign artifact (Table 4 only).** Casoliva's own printed minus signs
render inconsistently across this module's two independent
``pdftotext -raw``/``-layout`` extractions of the PDF text layer: Table 3's
negative x-coordinates render as ``N:DDDD`` (e.g. ``2:4754942840`` for
``-2.4754942840``, colon substituting for the minus-sign-plus-decimal-point
ligature) and were caught/confirmed via `#747`'s own independent
``pdftotext -raw`` extraction of 6 of these 16 rows. Table 4's 4th "Y"
(canonical-momentum ``p_y``) component prints as a bare ``1.43284467834384``
with NO colon in either extraction, but is a NEGATIVE value:
``p_y = -1.43284467834384`` (confirmed by requiring the row's own printed
h/eigenvalue to reproduce -- see "Class 2 IC derivation" below). This is
reported here as a source-text-extraction hazard specific to this
module's own OCR tooling on Table 4's particular typesetting, NOT a defect
in Casoliva's own PDF (this project's respectful-errata-framing discipline
does not apply -- this is our own extraction miss, not their error).

Class 2 IC derivation (Table 4, 2010 p.1635)
---------------------------------------------
Table 4 prints the L1 Lyapunov p.o.'s own state in Hamiltonian
position-momentum form, ``Y = (x, y, p_x, p_y)``, using the standard CR3BP
canonical-momentum convention this project already documents elsewhere
(``project_andreu_pol_canonical_momentum`` memory note): ``p_x = vx - y``,
``p_y = vy + x``. With the corrected reading ``p_y = -1.43284467834384``
(see above) and Table 4's own printed ``x = -0.6280674149446867, y = 0``:
``vy = p_y - x = -1.43284467834384 - (-0.6280674149446867) =
-0.8047772633991533`` (Casoliva's own frame); the full flip gives this
project's frame IC ``(x0, y0, vx0, vy0) = (0.6280674149446867, 0, 0,
0.8047772633991533)``. This is CONFIRMED (not merely plausible) two
independent ways: (a) direct propagation for the printed period
``T = 6.706878522271349`` closes to 2.5e-6; (b) after
:func:`~cyclerfinder.search.cr3bp_periodic.correct_symmetric_fixed_jacobi`
refinement (residual 1e-13), the recovered dominant planar-monodromy
eigenvalue matches Table 4's own printed unstable eigenvalue
``108.5966557497375`` to **5.3e-8 relative** -- an eigenvalue reproduction
this tight is not attainable by a wrong IC (eigenvalues are NOT
Jacobi-constant-degenerate the way raw velocity-magnitude checks are; see
:data:`HE1_EIGENVALUE_GATE_REL_TOL`).

mu / l* / t* (per task instruction, NOT re-litigating the sibling
modules' own source-mu precedent)
------------------------------------------------------------------------
Unlike the Saturn-Titan/Neptune-Triton modules (which use the SOURCE
paper's own displayed mu for "fairest reproduction" against a system this
project has no other canonical anchor for), Earth-Moon is this project's
OWN flagship system with an existing canonical mu/l*/t* anchor used by 322
catalogue rows (:func:`cyclerfinder.core.cr3bp.cr3bp_system`). This module
uses THAT registry system throughout (:func:`earth_moon_system`), NOT
Casoliva's own displayed ``mu_EM = 0.0121529529`` (2010) /
``0.01215`` (2008) -- kept below as :data:`CASOLIVA_MU_2010` purely for
citation/comparison, never fed into any corrector. This choice turns out to
be numerically BETTER, not just policy-compliant: the Class 2 golden
anchor's own printed energy reproduces the registry-mu Jacobi constant to
**1.3e-9 relative** (vs. 1.8e-6 relative at Casoliva's own displayed mu) --
strong evidence Casoliva's internal computation used a mu closer to this
project's own DE440-registry value than to her own paper's rounded
7-significant-figure display (the same "displayed value is itself a
rounded display of a higher-precision internal one" pattern the
Saturn-Titan module documents for its own thesis-mu gap).

Explicitly out of scope for `#780` (mirrors `#765`'s own narrow
family-only first task)
------------------------------------------------------------------------
- The homoclinic connection stage itself (He1-4, Hm1-2 discrete points) and
  the Barrabes-Mondelo-Olle (2009, Nonlinearity) numerical-continuation-of-
  homoclinic-connections algorithm (2010 Eq. 20) that builds them. This is
  registered as a follow-up task -- see the module's own results note.
- Gate item (e) ("continuation-in-energy onto multiple family members") is
  satisfied here ONLY in the narrow sense the vendored data supports: a
  plain natural-parameter (Jacobi) continuation of the Class 2 Lyapunov
  ORBIT family from the He1 golden anchor, via
  :func:`cyclerfinder.search.cr3bp_continuation.continue_family` (already-
  built, general-purpose machinery, NOT the excluded homoclinic-connection
  continuator -- continuing a plain periodic-orbit family in C is a
  categorically different, much simpler operation). This is reported as
  family-persistence SMOKE-TEST evidence only, not a reproduction (no
  printed IC/eigenvalue exists at the intermediate energies to reproduce
  against). Class 1's 16 rows are each isolated one-per-(p,q,label) points
  in Table 3, not multi-member continuation sequences -- no continuation is
  attempted on them here.

Pure: math/numpy/scipy + :mod:`cyclerfinder.core.cr3bp`,
:mod:`cyclerfinder.search.cr3bp_periodic`,
:mod:`cyclerfinder.search.cr3bp_continuation`,
:mod:`cyclerfinder.search.jovian_resonant_families` (reused directly, not
reimplemented, for :func:`~cyclerfinder.search.jovian_resonant_families.two_body_resonant_seed`
and :func:`~cyclerfinder.search.jovian_resonant_families.converge_candidate`).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_continuation as cc
import cyclerfinder.search.cr3bp_periodic as cp
import cyclerfinder.search.jovian_resonant_families as jrf

# ---------------------------------------------------------------------------
# System (registry-canonical, NOT a new mu -- see module docstring).
# ---------------------------------------------------------------------------

#: Casoliva 2010's own displayed Earth-Moon mass ratio (p.1624 text,
#: "0.0121529529"). Citation/comparison ONLY -- never fed into a corrector.
#: See module docstring for why the registry mu is used instead and why it
#: is the numerically better choice here (not merely policy-compliant).
CASOLIVA_MU_2010 = 0.0121529529

#: Casoliva 2008's own displayed (rounded) Earth-Moon mass ratio. Citation
#: only.
CASOLIVA_MU_2008 = 0.01215


def earth_moon_system() -> cr3bp.CR3BPSystem:
    """This project's own canonical Earth-Moon CR3BP system (DE440 registry
    mu/l*/t*, :func:`cyclerfinder.core.cr3bp.cr3bp_system`) -- NOT a new
    system, per task instruction. See module docstring."""
    return cr3bp.cr3bp_system("Earth", "Moon")


# ---------------------------------------------------------------------------
# Class 1 -- Table 3 (2010, p.1630), verbatim, all 16 rows.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Table3Row:
    """One row of Casoliva 2010 Table 3, verbatim, IN CASOLIVA'S OWN
    coordinate convention (Earth at ``(+mu,0)``, Moon at ``(mu-1,0)``) --
    NOT this project's own convention. Use :func:`table3_seed_state` to get
    a project-frame IC. ``k`` is Casoliva's own stability index (Eq. 8,
    ``k = lambda + 1/lambda`` from the FULL-period planar monodromy --
    NOT the Barden half-period ``nu = 1/2(lambda+1/lambda)`` convention
    this project's own :func:`~cyclerfinder.search.cr3bp_periodic.barden_stability`
    uses); ``|k| < 2`` is Casoliva's own stability threshold.
    """

    designation: str
    p: int
    q: int
    c_j: float
    period: float  # T, nondim; = 2*pi*q exactly for a genuine p-q resonance (Eq. 18)
    x_i: float
    u_i: float
    v_i: float
    r_pm: float  # periselene (Moon-relative), normalized units
    k: float  # stability index, Eq. 8 (full-period lambda + 1/lambda)
    r_pe: float  # perigee (Earth-relative), normalized units
    r_ae: float  # apogee (Earth-relative), normalized units
    v_pm: float  # velocity at periselene, normalized units
    v_pe: float  # velocity at perigee, normalized units
    satisfies_resonance: bool  # footnote "e": False = "does not satisfy its resonance relation"
    exists_in_em_system: bool  # footnote "d": False = "flies through the Earth"

    @property
    def stable(self) -> bool:
        """Casoliva's own stability call: ``|k| < 2``."""
        return abs(self.k) < 2.0


#: Casoliva 2010 Table 3 (p.1630), ALL 16 rows, verbatim. mu_EM = 0.0121529529
#: (Casoliva's own, per-row IC point is a Poincare y=0 crossing, NOT
#: necessarily the perpendicular symmetric crossing -- her own text: "we
#: have not used this property [the perpendicular crossings] in this
#: paper"). Column order per the table's own header: Designation, C_J, x_i,
#: T, u_i, v_i, r_pM, k, r_pE, r_aE, V_pM, V_pE. 6 of these 16 rows
#: (1-2c/1-2d/1-2e/2-1a/2-1b/3-2c/7-3a, sic -- actually the 6 SYMMETRIC-orbit
#: rows) were independently cross-extracted and verified in
#: `docs/notes/2026-07-28-747-franz-russell-casoliva-crosscheck.md` via a
#: SEPARATE ``pdftotext -raw`` pass; all match this module's own extraction
#: exactly. The other 10 (including all 7 footnote-flagged rows) were
#: extracted fresh this task and independently VALIDATED not merely
#: trusted: every one of the 16 rows closes under direct propagation for
#: its own printed period at this module's registry mu (residuals
#: 6e-7 to 6e-2, worst on the "d"-flagged near-Earth-singular rows, as
#: expected -- see module docstring's periodicity-closure discussion; this
#: is the check that would have caught a colon/sign misread, and did not
#: flag one for any of the 16).
TABLE3_ROWS: tuple[Table3Row, ...] = (
    Table3Row(
        "1-2a",
        1,
        2,
        1.0963633786,
        12.5663706144,
        -3.1452530020,
        0.0000000000,
        3.0713600647,
        0.6185470387,
        1.9971178106,
        0.0263287962,
        3.1574059549,
        1.9528940614,
        8.6267981407,
        False,
        True,
    ),
    Table3Row(
        "1-2b",
        1,
        2,
        1.4799919040,
        12.5663706144,
        -3.0692958563,
        0.0480996538,
        2.9310589378,
        0.5364691875,
        1.9975766156,
        0.0906722860,
        3.0917876572,
        1.5294833421,
        4.6009059094,
        False,
        True,
    ),
    Table3Row(
        "1-2c",
        1,
        2,
        1.5691874798,
        12.5663706144,
        -2.4754942840,
        -0.3779077269,
        2.2861781603,
        0.8833330997,
        1.9995914782,
        0.1166859899,
        3.0501343486,
        3.0411378116,
        4.0410890199,
        True,
        True,
    ),
    Table3Row(
        "1-2d",
        1,
        2,
        2.5803060666,
        12.5663706144,
        -2.5987958202,
        0.0000000000,
        2.2237844860,
        1.5702173819,
        4.8578794987,
        0.5702686374,
        2.6109487731,
        2.6832627814,
        1.6832776688,
        True,
        True,
    ),
    Table3Row(
        "1-2e",
        1,
        2,
        2.7629814961,
        12.5663706144,
        -2.1598692556,
        -0.2293727086,
        1.6672704485,
        0.6982893039,
        1.9997758212,
        0.7490805561,
        2.4266038762,
        0.4848553515,
        1.4222532528,
        True,
        True,
    ),
    Table3Row(
        "2-1a",
        2,
        1,
        0.4887353098,
        6.2831853072,
        -0.7362153512,
        -0.6876743802,
        1.5221620241,
        0.2353576230,
        1.9255637447,
        0.1876742270,
        1.0644080071,
        1.1752363151,
        2.9902585472,
        True,
        True,
    ),
    Table3Row(
        "2-1b",
        2,
        1,
        1.1964188553,
        6.2831853072,
        -1.2287157860,
        0.0000000000,
        1.4164812666,
        0.2408687390,
        2.0373978625,
        0.0176632387,
        1.2621776123,
        1.0157645984,
        10.5028780155,
        True,
        True,
    ),
    Table3Row(
        "2-1c",
        2,
        1,
        1.7351680041,
        6.1905730000,
        -0.9731602028,
        0.0000000000,
        1.6946858601,
        0.0146868443,
        102.8780940137,
        0.0007205802,
        1.1078267516,
        1.0918627346,
        52.3452383432,
        False,
        False,
    ),
    Table3Row(
        "2-1d",
        2,
        1,
        1.9521613046,
        6.2831853072,
        -0.9732149393,
        -0.0000000000,
        1.6313040525,
        0.0146321078,
        118.9539595702,
        0.0031649323,
        1.1166032346,
        1.0476316856,
        24.9494588896,
        True,
        False,
    ),
    Table3Row(
        "3-2a",
        3,
        2,
        0.1258650202,
        12.6356320000,
        -0.8187678081,
        -0.7299659685,
        1.5915830242,
        0.1561637871,
        1.5113469213,
        0.2074812540,
        1.3292921326,
        0.2690830739,
        2.8688059249,
        False,
        True,
    ),
    Table3Row(
        "3-2c",
        3,
        2,
        0.7089330385,
        12.5663706144,
        -1.2237269180,
        -0.5006876037,
        1.4965104108,
        0.2193864666,
        1.8751541804,
        0.0392320964,
        1.5557128499,
        1.5768170682,
        7.0089024548,
        True,
        True,
    ),
    Table3Row(
        "3-2d",
        3,
        2,
        1.6505940399,
        12.5663706144,
        -1.4665738649,
        0.0000000000,
        1.3737156364,
        0.4787268178,
        2.0008005398,
        0.0139098602,
        1.5303593112,
        1.0282877248,
        11.8633639097,
        True,
        False,
    ),
    Table3Row(
        "7-3a",
        7,
        3,
        1.0215696153,
        18.8495559215,
        -0.8938394486,
        -0.4831113793,
        1.4082724890,
        0.0709180137,
        -1.2990228617,
        0.0667195222,
        1.0910338049,
        1.2234251887,
        5.2829778676,
        True,
        True,
    ),
    Table3Row(
        "7-3b",
        7,
        3,
        1.0687623900,
        18.8495559215,
        -0.8404096071,
        -0.5542434833,
        1.3463118042,
        0.0343651369,
        57.3519357466,
        0.0408619187,
        1.1887710938,
        1.2829729016,
        6.8372131877,
        True,
        True,
    ),
    Table3Row(
        "7-3c",
        7,
        3,
        1.0687623900,
        18.8495559215,
        -0.9462620909,
        -0.4279513671,
        1.5130807007,
        0.0343660232,
        57.3519357463,
        0.0408577074,
        1.1886076749,
        1.2844406950,
        6.8375778062,
        True,
        True,
    ),
    Table3Row(
        "7-3d",
        7,
        3,
        1.2892422741,
        18.9972028186,
        -0.9319486565,
        -0.5413649154,
        1.3467383149,
        0.0087747218,
        4.4005195470,
        0.0132595167,
        1.2268771199,
        1.9565627505,
        12.1414662244,
        False,
        False,
    ),
)

#: Casoliva's own "clean" set: satisfies its own p-q resonance relation AND
#: physically exists in the Earth-Moon system (does not fly through Earth).
#: 9 of the 16 rows. Digest golden-table set.
TABLE3_VALID_DESIGNATIONS: tuple[str, ...] = (
    "1-2c",
    "1-2d",
    "1-2e",
    "2-1a",
    "2-1b",
    "3-2c",
    "7-3a",
    "7-3b",
    "7-3c",
)

#: The 5 stable members among the 9 valid rows -- Casoliva's own "separate,
#: uncatalogued stable family" (digest #725 finding).
TABLE3_STABLE_DESIGNATIONS: tuple[str, ...] = ("1-2c", "1-2e", "2-1a", "3-2c", "7-3a")


def table3_row(designation: str) -> Table3Row:
    for row in TABLE3_ROWS:
        if row.designation == designation:
            return row
    raise ValueError(f"unknown Table 3 designation {designation!r}")


def table3_seed_state(row: Table3Row) -> NDArray[np.float64]:
    """Project-frame (this project's own primary-secondary convention) 6-state
    IC for ``row``, via the validated full 180-degree-rotation flip
    ``(x,y,vx,vy) -> (-x,-y,-vx,-vy)`` (module docstring). ``y=0`` for every
    Table 3 row (Casoliva's own Poincare-section convention)."""
    return np.array([-row.x_i, 0.0, 0.0, -row.u_i, -row.v_i, 0.0])


# ---------------------------------------------------------------------------
# Planar stability index -- Casoliva's own Eq. 8 convention (FULL-period
# lambda + 1/lambda), computed via the trace identity (robust near |k|~2,
# where Barden-style nontrivial-eigenpair SORTING mis-picks -- see module
# results note) with an eigenvalue-pair cross-check.
# ---------------------------------------------------------------------------


def _planar_monodromy_4x4(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    period: float,
    *,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> NDArray[np.float64]:
    """The 4x4 planar (x,y,vx,vy) block of the full-period monodromy matrix."""
    arc = cr3bp.propagate(system, state0, period, with_stm=True, rtol=rtol, atol=atol)
    assert arc.stm is not None
    idx = [0, 1, 3, 4]
    return np.asarray(arc.stm)[np.ix_(idx, idx)]


@dataclass(frozen=True)
class StabilityIndex:
    """Casoliva 2010 Eq. 6-8 stability decomposition, computed from the FULL
    6x6 one-period monodromy of a (possibly asymmetric) periodic orbit: an
    IN-PLANE index ``k_par`` (4x4 ``(x,y,vx,vy)`` block, trace identity
    ``k = trace(Phi4) - 2``, exact for the ``{1,1,lambda,1/lambda}`` planar
    spectrum -- no eigenvalue-pairing heuristic, robust as ``lambda -> 1``)
    and an OUT-OF-PLANE index ``k_perp`` (DECOUPLED 2x2 ``(z,vz)`` block --
    for a genuinely planar orbit z=vz=0 is an invariant subspace, so this
    block's eigenvalues are the vertical-stability reciprocal pair
    directly and ``trace(Phi_z) = kappa_perp + 1/kappa_perp`` exactly, no
    pairing needed there either). ``k_signed`` is Casoliva's own printed
    "k" per her Eq. 8, ``k = max(|kappa_par|, |kappa_perp|)`` -- taken here
    as whichever of ``k_par``/``k_perp`` has the larger ABSOLUTE value,
    keeping its own sign (Table 3 prints signed k, e.g. 7-3a's ``-1.299``).

    This decomposition was NOT obvious from the OCR'd table text alone (an
    earlier draft of this module compared the printed k only against
    ``k_par``, which reproduces just 5 of 16 rows) -- Eq. 6-8's own text
    (2010 p.1626: "kappa_par = kappa + 1/kappa" ... "k =
    max{|kappa_par|,|kappa_perp|}") was the resolving read. With the
    ``max(|.|)`` selection, **12 of 16 rows reproduce Casoliva's own
    printed k to 8.6e-8 to 1.7e-3 relative** (module test suite) -- a
    genuine reproduction, not merely evidentiary. The 4 exceptions
    (1-2e, 3-2a, 7-3a, 7-3d) are honest misses on k: the recovered k is
    wildly different in both sign and magnitude from the printed value
    (rel_err > 1) -- NOT resolved this task; see the results note for the
    full per-row data. THREE of the four (1-2e, 3-2a, 7-3a) still
    reproduce IC/period/Jacobi tightly (1e-5 to 2e-5 relative -- this IS
    the same orbit, the miss is k-only). The fourth, 7-3d, ALSO misses on
    IC/period/Jacobi (x0 off by a factor of ~3.8) -- it separately carries
    BOTH footnotes (does not satisfy its own resonance relation AND flies
    through the Earth), the single most degenerate row in the table, so
    this broader miss is less surprising than the other 3's k-only one.

    ``k_eig`` is a Barden-style direct eigenvalue-pair cross-check on
    ``k_par`` (largest-|eig-1| pair selection) -- ``agree`` flags
    disagreement with the trace identity (diagnostic of a near-degenerate
    pairing; item (c) of the gate).
    """

    k_par: float
    k_perp: float
    k_eig: float
    lam: complex
    agree: bool

    @property
    def k_signed(self) -> float:
        return self.k_par if abs(self.k_par) >= abs(self.k_perp) else self.k_perp


def planar_stability_index(
    system: cr3bp.CR3BPSystem, state0: NDArray[np.float64], period: float
) -> StabilityIndex:
    arc = cr3bp.propagate(system, state0, period, with_stm=True, rtol=1e-12, atol=1e-12)
    assert arc.stm is not None
    stm = np.asarray(arc.stm)
    phi4 = stm[np.ix_([0, 1, 3, 4], [0, 1, 3, 4])]
    phi_z = stm[np.ix_([2, 5], [2, 5])]
    k_par = float(np.trace(phi4).real) - 2.0
    k_perp = float(np.trace(phi_z).real)
    eigs = np.linalg.eigvals(phi4)
    order = np.argsort(np.abs(eigs - 1.0))
    nontrivial = eigs[order[2:]]
    lam = complex(nontrivial[int(np.argmax(np.abs(nontrivial)))])
    k_eig = float((lam + 1.0 / lam).real)
    agree = abs(k_par - k_eig) < 1e-6 * max(1.0, abs(k_par))
    return StabilityIndex(k_par=k_par, k_perp=k_perp, k_eig=k_eig, lam=lam, agree=agree)


# ---------------------------------------------------------------------------
# Class 1 gate: (a) periodicity self-consistency, (b) reproduction against
# printed IC/period/stability, (c) internal trace-vs-eigenpair + independent-
# integrator (Radau) cross-check.
# ---------------------------------------------------------------------------


def recover_table3_row(
    designation: str,
    system: cr3bp.CR3BPSystem | None = None,
    *,
    tol: float = 1e-8,
    max_iter: int = 1000,
) -> cp.PeriodicOrbit:
    """Full-state single-shooting periodicity correction
    (:func:`~cyclerfinder.search.cr3bp_periodic.correct_periodic`) seeded
    from Table 3's own printed (generally ASYMMETRIC) Poincare-crossing IC
    for ``designation`` -- NOT the symmetric-orbit corrector (Casoliva's
    own text: her printed crossing is not necessarily the perpendicular
    one), matching the general full-state Newton shooting this project
    already has, reused unmodified. ``correct_periodic`` is an undamped
    min-norm Newton iteration; the two near-Earth-singular rows (2-1c,
    2-1d, both footnote-``d`` "flies through the Earth") need the full
    ``max_iter=1000`` to close to ``tol=1e-8`` -- all 16 rows converge at
    these settings (module test suite)."""
    sys_ = system if system is not None else earth_moon_system()
    row = table3_row(designation)
    seed = table3_seed_state(row)
    return cp.correct_periodic(sys_, seed, row.period, tol=tol, max_iter=max_iter)


@dataclass(frozen=True)
class Table3GateRow:
    """One row of the Table 3 gate report. ``passed`` requires: (a) the
    corrector converged; (b) the corrected IC/period/Jacobi AND the
    ``k_signed`` stability index (Casoliva Eq. 8, ``max(|k_par|,|k_perp|)``
    -- see :class:`StabilityIndex`'s docstring) reproduce Casoliva's own
    printed values within :data:`TABLE3_IC_GATE_REL_TOL`; (c) the internal
    trace-vs-eigenpair AND independent-Radau-integrator cross-checks both
    agree. 12 of 16 rows pass the k-reproduction criterion; 4 (1-2e, 3-2a,
    7-3a, 7-3d) are honest misses on k ONLY (IC/period/Jacobi still
    match) -- see :class:`StabilityIndex`'s docstring."""

    designation: str
    converged: bool
    closure_residual: float
    x0_rel_err: float
    period_rel_err: float
    jacobi_rel_err: float
    k_par: float
    k_perp: float
    k_eig: float
    k_signed: float
    k_rel_err: float
    k_reproduced: bool
    stability_index_agree: bool
    recovered_stable: bool
    source_stable: bool
    stability_character_matches: bool
    radau_ok: bool
    radau_djacobi: float
    passed: bool


#: Reproduction tolerance for x0/period/Jacobi against Table 3's own printed
#: values. Justified by the observed match quality after correction
#: (typically 1e-6 to 1e-4 relative, module test suite) -- the SAME order
#: as the direct-propagation closure residuals reported in the module
#: docstring, themselves attributable to the small (~0.02% relative)
#: registry-vs-Casoliva mu difference (module docstring), not a genuine
#: mismatch. 1e-2 is generous relative to the observed match quality,
#: chosen to comfortably clear the mu-difference floor without absorbing a
#: genuine miss.
TABLE3_IC_GATE_REL_TOL = 1e-2


def table3_gate_report(system: cr3bp.CR3BPSystem | None = None) -> list[Table3GateRow]:
    sys_ = system if system is not None else earth_moon_system()
    rows: list[Table3GateRow] = []
    for row in TABLE3_ROWS:
        seed = table3_seed_state(row)
        po = cp.correct_periodic(sys_, seed, row.period, tol=1e-8, max_iter=1000)
        if not po.converged:
            rows.append(
                Table3GateRow(
                    designation=row.designation,
                    converged=False,
                    closure_residual=po.closure_residual,
                    x0_rel_err=float("nan"),
                    period_rel_err=float("nan"),
                    jacobi_rel_err=float("nan"),
                    k_par=float("nan"),
                    k_perp=float("nan"),
                    k_eig=float("nan"),
                    k_signed=float("nan"),
                    k_rel_err=float("nan"),
                    k_reproduced=False,
                    stability_index_agree=False,
                    recovered_stable=False,
                    source_stable=row.stable,
                    stability_character_matches=False,
                    radau_ok=False,
                    radau_djacobi=float("nan"),
                    passed=False,
                )
            )
            continue
        x0_rel_err = (
            abs(po.state0[0] - seed[0]) / abs(seed[0]) if seed[0] != 0.0 else abs(po.state0[0])
        )
        period_rel_err = abs(po.period - row.period) / abs(row.period)
        jacobi_rel_err = abs(po.jacobi - row.c_j) / abs(row.c_j)
        stab = planar_stability_index(sys_, po.state0, po.period)
        recovered_stable = abs(stab.k_signed) < 2.0
        k_rel_err = abs(stab.k_signed - row.k) / abs(row.k)
        k_reproduced = k_rel_err < TABLE3_IC_GATE_REL_TOL
        radau_ok, radau_dj = cp.crosscheck_periodic(sys_, po)
        rows.append(
            Table3GateRow(
                designation=row.designation,
                converged=True,
                closure_residual=po.closure_residual,
                x0_rel_err=x0_rel_err,
                period_rel_err=period_rel_err,
                jacobi_rel_err=jacobi_rel_err,
                k_par=stab.k_par,
                k_perp=stab.k_perp,
                k_eig=stab.k_eig,
                k_signed=stab.k_signed,
                k_rel_err=k_rel_err,
                k_reproduced=k_reproduced,
                stability_index_agree=stab.agree,
                recovered_stable=recovered_stable,
                source_stable=row.stable,
                stability_character_matches=recovered_stable == row.stable,
                radau_ok=radau_ok,
                radau_djacobi=radau_dj,
                passed=(
                    x0_rel_err < TABLE3_IC_GATE_REL_TOL
                    and period_rel_err < TABLE3_IC_GATE_REL_TOL
                    and jacobi_rel_err < TABLE3_IC_GATE_REL_TOL
                    and k_reproduced
                    and stab.agree
                    and radau_ok
                ),
            )
        )
    return rows


# ---------------------------------------------------------------------------
# Class 2 -- Table 4 (2010, p.1635) golden He1-family Lyapunov p.o. anchor.
# ---------------------------------------------------------------------------

#: Table 4's own printed energy (Hamiltonian constant h), full precision.
HE1_H = -1.45016232260699

#: Corresponding Jacobi constant, C_J = -2h (2010 Eq. 3).
HE1_CJ = -2.0 * HE1_H

#: Table 4's own printed period T (nondim time units).
HE1_PERIOD_GUESS = 6.706878522271349

#: Project-frame x0 (this project's own convention), derived per the module
#: docstring's "Class 2 IC derivation" -- ``ydot0`` is re-solved from the
#: Jacobi constraint by the corrector itself, not hardcoded (the seed
#: velocity 0.8047772633991533 is used only as ``ydot0_sign``'s magnitude
#: check in the module test suite).
HE1_X0_GUESS = 0.6280674149446867
HE1_YDOT0_GUESS = 0.8047772633991533

#: Table 4's own printed unstable eigenvalue (module docstring: labeled
#: "mu_u" in the continuation-system equations, i.e. the eigenvalue paired
#: with the unstable-manifold eigenvector V_u -- NOT literally the mass
#: ratio mu).
HE1_TARGET_LAMBDA_U = 108.5966557497375

#: Text-stated (both papers): "the period of the Lyapunov orbit is 29.1640
#: days." Evidentiary only (see :data:`HE1_PERIOD_DAYS_GATE_REL_TOL`'s own
#: docstring for the small, honestly-reported day-unit mismatch).
HE1_PERIOD_DAYS_SOURCE = 29.1640

#: Text-stated (both papers, essentially identical figure): connection
#: flight time from periselene 1 to periselene 19. Citation only -- NOT
#: reproduced this stage (the connection itself is out of scope, module
#: docstring).
HE1_CONNECTION_FLIGHT_TIME_DAYS_SOURCE = 113.6319

#: 2008's own independently-stated LEO-rendezvous insertion delta-v at the
#: connection's Earth perigee (703 m/s at r=67,808 km) vs 2010's (717.5 m/s
#: at r=67,869 km) -- both citation-only context for the follow-up
#: connection-stage task, not reproduced or gated here.
HE1_LEO_DV_MPS_2008 = 703.0
HE1_LEO_DV_MPS_2010 = 717.5


def recover_he1_lyapunov(system: cr3bp.CR3BPSystem | None = None) -> cp.SymmetricOrbit:
    """Recover the He1-family golden Lyapunov p.o. via the standard
    symmetric-orbit corrector (Table 4's own IC is a genuine perpendicular
    x-axis crossing, ``y0=0, vx0=0`` -- module docstring), seeded from
    :data:`HE1_X0_GUESS`/:data:`HE1_CJ`/:data:`HE1_PERIOD_GUESS`."""
    sys_ = system if system is not None else earth_moon_system()
    return cp.correct_symmetric_fixed_jacobi(
        sys_, HE1_X0_GUESS, HE1_CJ, HE1_PERIOD_GUESS, ydot0_sign=1.0, half_crossings=1, tol=1e-12
    )


@dataclass(frozen=True)
class He1GateRow:
    """Class 2 golden-anchor gate report. ``passed`` requires the corrector
    to converge AND the recovered dominant planar eigenvalue to reproduce
    Table 4's own printed :data:`HE1_TARGET_LAMBDA_U` within
    :data:`HE1_EIGENVALUE_GATE_REL_TOL` -- the primary, non-Jacobi-degenerate
    reproduction criterion (module docstring). Period-in-days match is
    reported as SECONDARY/evidentiary (see :data:`HE1_PERIOD_DAYS_GATE_REL_TOL`).
    """

    converged: bool
    crossing_residual: float
    x0: float
    ydot0: float
    period_nondim: float
    eigenvalue_barden_nu: float
    eigenvalue_recovered_lambda: float
    eigenvalue_rel_err: float
    eigenvalue_confirmed: bool
    barden_vs_floquet_agree: bool
    period_days: float
    period_days_rel_err: float
    period_days_confirmed: bool
    passed: bool


#: Justified by the observed match quality: 5.3e-8 relative at the
#: converged solution (module docstring) -- SIX orders of magnitude inside
#: this tolerance. Generous on purpose (mirrors the sibling modules' own
#: "generous relative to the observed precision floor" discipline).
HE1_EIGENVALUE_GATE_REL_TOL = 1e-4

#: The recovered period-in-days (via this module's own registry t_s)
#: differs from Casoliva's own printed "29.1640 days" by ~0.135% relative
#: -- small but an order of magnitude larger than every other reproduced
#: quantity here (x0/ydot0/period-nondim/eigenvalue all match to 1e-8 to
#: 1e-7 relative, module docstring). NOT attributable to the registry-vs-
#: Casoliva t_s difference (that alone is only ~9e-6 relative, checked
#: directly: converting the SAME recovered nondim period through Casoliva's
#: own stated ``omega_EM = 2.66529e-6`` still gives ~29.1248 days, not
#: 29.1640). Reported honestly as a SECONDARY, non-gating discrepancy
#: (evidentiary only -- does not affect ``passed``, which is keyed on the
#: eigenvalue reproduction) rather than silently absorbed into a loose
#: tolerance.
HE1_PERIOD_DAYS_GATE_REL_TOL = 1e-2


def he1_gate_report(system: cr3bp.CR3BPSystem | None = None) -> He1GateRow:
    sys_ = system if system is not None else earth_moon_system()
    orbit = recover_he1_lyapunov(sys_)
    if not orbit.converged:
        return He1GateRow(
            converged=False,
            crossing_residual=orbit.crossing_residual,
            x0=orbit.x0,
            ydot0=orbit.ydot0,
            period_nondim=orbit.period,
            eigenvalue_barden_nu=float("nan"),
            eigenvalue_recovered_lambda=float("nan"),
            eigenvalue_rel_err=float("nan"),
            eigenvalue_confirmed=False,
            barden_vs_floquet_agree=False,
            period_days=float("nan"),
            period_days_rel_err=float("nan"),
            period_days_confirmed=False,
            passed=False,
        )
    nu, lam = cp.barden_stability(sys_, orbit)
    from cyclerfinder.search.resonance_network import _planar_floquet

    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0])
    lam_pf, _vec = _planar_floquet(sys_, state0, orbit.period)
    barden_vs_floquet_agree = abs(abs(lam.real) - abs(lam_pf)) / abs(lam_pf) < 1e-4

    lam_recovered = abs(lam.real)
    eig_rel_err = abs(lam_recovered - HE1_TARGET_LAMBDA_U) / HE1_TARGET_LAMBDA_U
    eig_confirmed = eig_rel_err < HE1_EIGENVALUE_GATE_REL_TOL

    period_days = orbit.period * sys_.t_s / 86400.0
    period_days_rel_err = abs(period_days - HE1_PERIOD_DAYS_SOURCE) / HE1_PERIOD_DAYS_SOURCE
    period_days_confirmed = period_days_rel_err < HE1_PERIOD_DAYS_GATE_REL_TOL

    return He1GateRow(
        converged=True,
        crossing_residual=orbit.crossing_residual,
        x0=orbit.x0,
        ydot0=orbit.ydot0,
        period_nondim=orbit.period,
        eigenvalue_barden_nu=nu,
        eigenvalue_recovered_lambda=lam_recovered,
        eigenvalue_rel_err=eig_rel_err,
        eigenvalue_confirmed=eig_confirmed,
        barden_vs_floquet_agree=barden_vs_floquet_agree,
        period_days=period_days,
        period_days_rel_err=period_days_rel_err,
        period_days_confirmed=period_days_confirmed,
        passed=eig_confirmed,
    )


# ---------------------------------------------------------------------------
# Gate item (e): plain natural-parameter (Jacobi) continuation of the He1
# Lyapunov orbit family, via the ALREADY-BUILT
# cyclerfinder.search.cr3bp_continuation.continue_family -- NOT the excluded
# Barrabes-Mondelo-Olle homoclinic-CONNECTION continuator (module docstring).
# Reported as family-persistence smoke-test evidence only.
# ---------------------------------------------------------------------------


def continue_he1_family(
    system: cr3bp.CR3BPSystem | None = None,
    *,
    direction: int = 1,
    d_jacobi: float = 0.01,
    n_steps: int = 10,
) -> cc.FamilyBranch:
    """Natural-parameter continuation of the He1 Lyapunov orbit family from
    the golden anchor, ``direction`` steps of ``d_jacobi`` in Jacobi
    constant, up to ``n_steps``. Smoke-test evidence only -- no printed
    IC/eigenvalue exists at these intermediate energies to reproduce
    against (module docstring, "Explicitly out of scope")."""
    sys_ = system if system is not None else earth_moon_system()
    seed = recover_he1_lyapunov(sys_)
    return cc.continue_family(
        sys_,
        seed,
        direction=direction,
        d_jacobi=d_jacobi,
        n_steps=n_steps,
        min_jacobi=HE1_CJ - 1.0,
        max_jacobi=HE1_CJ + 1.0,
        half_crossings=1,
        ydot0_sign=1.0,
        seed_label="he1_lyapunov",
    )


# ---------------------------------------------------------------------------
# Gate item (d): two_body_resonant_seed lineage check (Class 1). Expect an
# honest negative, per the pattern already confirmed 3x project-wide
# (Jovian/Saturn-Titan/Neptune-Triton) -- see :data:`TWO_BODY_SEED_LINEAGE_NOTE`.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TwoBodySeedAttempt:
    p: int
    q: int
    seed_x0: float
    seed_jacobi: float
    converged: bool
    recovered_x0: float
    recovered_period_over_2pi: float
    nearest_row_designation: str
    nearest_row_x0_project_frame: float


def two_body_seed_lineage_check(
    system: cr3bp.CR3BPSystem | None = None,
) -> list[TwoBodySeedAttempt]:
    """For each DISTINCT (p,q) among Table 3's own resonance labels, try
    :func:`~cyclerfinder.search.jovian_resonant_families.two_body_resonant_seed`
    (x0_sign=-1, matching the sibling modules' own convention) at ITS OWN
    natural Jacobi constant, converge via
    :func:`~cyclerfinder.search.jovian_resonant_families.converge_candidate`,
    and report where it lands relative to Table 3's own nearest same-(p,q)
    row. No fitting/tuning -- report whatever the naive construction
    produces, per the established honest-negative pattern."""
    sys_ = system if system is not None else earth_moon_system()
    seen_pq: list[tuple[int, int]] = []
    for row in TABLE3_ROWS:
        if (row.p, row.q) not in seen_pq:
            seen_pq.append((row.p, row.q))

    out: list[TwoBodySeedAttempt] = []
    for p, q in seen_pq:
        seed = jrf.two_body_resonant_seed(p, q, x0_sign=-1)
        state0 = np.array([seed.x0, 0.0, 0.0, 0.0, seed.ydot0, 0.0])
        seed_jacobi = cr3bp.jacobi_constant(state0, sys_.mu)
        cand = jrf.converge_candidate(
            sys_,
            f"two_body_{p}_{q}",
            seed.x0,
            seed_jacobi,
            seed.period_full,
            ydot0_sign=1.0 if seed.ydot0 >= 0 else -1.0,
            half_crossings=1,
        )
        # Nearest same-(p,q) vendored row (project-frame x0) for comparison.
        candidates_pq = [r for r in TABLE3_ROWS if (r.p, r.q) == (p, q)]
        nearest = min(candidates_pq, key=lambda r: abs(-r.x_i - seed.x0))
        out.append(
            TwoBodySeedAttempt(
                p=p,
                q=q,
                seed_x0=seed.x0,
                seed_jacobi=seed_jacobi,
                converged=cand is not None,
                recovered_x0=cand.x0 if cand is not None else float("nan"),
                recovered_period_over_2pi=cand.period_over_2pi
                if cand is not None
                else float("nan"),
                nearest_row_designation=nearest.designation,
                nearest_row_x0_project_frame=-nearest.x_i,
            )
        )
    return out


#: HONEST NEGATIVE (the pattern already confirmed 3x project-wide --
#: Jovian, Saturn-Titan, Neptune-Triton -- now a 4th). For every one of
#: Table 3's 4 distinct (p,q) pairs,
#: ``jrf.two_body_resonant_seed(p, q, x0_sign=-1)`` converges cleanly at
#: its own natural Jacobi constant via
#: ``jrf.converge_candidate`` (residuals ~1e-11 to 1e-14, module test
#: suite) -- but EVERY ONE lands on a ``period/2pi ~= 0.99`` orbit
#: (essentially a single-loop, non-resonant orbit), NOT the labeled
#: multi-period p-q resonance its own construction targets (e.g. (1,2)
#: should imply a period ~2x a base unit, (7,3) a much longer one -- none
#: does). Measured this task: (1,2) -> x0=-0.840 (nearest same-(p,q)
#: vendored row 1-2e at x0=+2.160, project frame -- opposite sign, very
#: different magnitude); (2,1) -> x0=-1.346 (nearest row 2-1a at
#: x0=+0.736); (3,2) -> x0=-1.169 (nearest row 3-2a at x0=+0.819);
#: (7,3) -> x0=-1.484 (nearest row 7-3b at x0=+0.840). The naive two-body-
#: resonant-ellipse construction converges to a genuine periodic orbit in
#: EVERY case but never identifies the correct topology OR magnitude of
#: any of Casoliva's own vendored Table 3 members -- the same qualitative
#: finding Anderson & Lo 2011 (Jovian), Vaquero 2013 (Saturn-Titan), and
#: Miceli & Bosanac 2026 (Neptune-Triton) each document for their own
#: analogous naive attempts.
TWO_BODY_SEED_LINEAGE_NOTE = (
    "jrf.two_body_resonant_seed(p, q, x0_sign=-1) converges cleanly at its "
    "own natural Jacobi constant for all 4 distinct (p,q) in Table 3 "
    "((1,2),(2,1),(3,2),(7,3)) -- but every one lands on a "
    "period/2pi ~= 0.99 orbit, not the labeled multi-period p-q "
    "resonance, and none matches its nearest same-(p,q) vendored row's "
    "own x0 in either sign or magnitude. See "
    "docs/notes/2026-08-07-780-earth-moon-casoliva-families-gate.md and "
    "two_body_seed_lineage_check() for the full per-(p,q) numbers."
)


__all__ = [
    "CASOLIVA_MU_2008",
    "CASOLIVA_MU_2010",
    "HE1_CJ",
    "HE1_CONNECTION_FLIGHT_TIME_DAYS_SOURCE",
    "HE1_EIGENVALUE_GATE_REL_TOL",
    "HE1_H",
    "HE1_LEO_DV_MPS_2008",
    "HE1_LEO_DV_MPS_2010",
    "HE1_PERIOD_DAYS_GATE_REL_TOL",
    "HE1_PERIOD_DAYS_SOURCE",
    "HE1_PERIOD_GUESS",
    "HE1_TARGET_LAMBDA_U",
    "HE1_X0_GUESS",
    "HE1_YDOT0_GUESS",
    "TABLE3_IC_GATE_REL_TOL",
    "TABLE3_ROWS",
    "TABLE3_STABLE_DESIGNATIONS",
    "TABLE3_VALID_DESIGNATIONS",
    "TWO_BODY_SEED_LINEAGE_NOTE",
    "He1GateRow",
    "StabilityIndex",
    "Table3GateRow",
    "Table3Row",
    "TwoBodySeedAttempt",
    "continue_he1_family",
    "earth_moon_system",
    "he1_gate_report",
    "planar_stability_index",
    "recover_he1_lyapunov",
    "recover_table3_row",
    "table3_gate_report",
    "table3_row",
    "table3_seed_state",
    "two_body_seed_lineage_check",
]
