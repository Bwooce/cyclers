"""Saturn-Titan planar CR3BP unstable resonant-orbit families (#765).

Thin sibling of :mod:`cyclerfinder.search.jovian_resonant_families` (#753),
retargeting the SAME reused corrector/classification machinery (no
reimplementation -- see the imports below) to the Saturn-Titan system, and
attempting to reproduce Table 4.1 of

    T. M. Vaquero Escribano (2013). "Spacecraft Transfer Trajectory Design
    Exploiting Resonant Orbits in Multi-Body Environments," PhD dissertation,
    Purdue University (advisor K. C. Howell), August 2013.

PDF held at
``cyclers_pdf/papers/vaquero-2013-spacecraft-transfer-trajectory-design-resonant-orbits-multibody-environments-purdue-phd.pdf``
(born-digital text layer, md5 ``fdcbf871322b87cd1dd3448059cb2596`` -- matches
the copy `#764`'s own scoping note downloaded and text-verified; re-verified
independently this task, both via ``pdftotext`` text-layer grep AND a direct
page-image (vision) read of the rendered Table 4.1 on printed p.109 (PDF page
124), per this project's corpus-document-policy precision-critical-table
rule). Journal companion (the citable venue, same content, not independently
re-checked for typos this task): Vaquero & Howell, "Transfer Design
Exploiting Resonant Orbits and Manifolds in the Saturn-Titan System,"
*J. Spacecraft & Rockets* 50(5):1069-1085 (2013), DOI ``10.2514/1.A32412``.

UNLIKE Anderson & Lo 2011 (the Jovian module's source), which published NO
initial-condition table (`#753`-`#758` had to locate every seed by blind
grid search + continuation), Vaquero 2013's own Table 4.1 ("Initial State,
Period, and Unstable Eigenvalue for Selected Periodic Orbits, C = 3.010000",
p.109) prints the seed ``x`` (km) and ``ẏ`` (km/s) directly -- per `#758`'s
own lesson ("paper-sourced seed windows beat blind grid scans"), this module
leads with those sourced seeds rather than a blind sweep. A blind
:func:`~cyclerfinder.search.jovian_resonant_families.two_body_resonant_seed`
attempt at the literal periapsis-at-Titan two-body geometry (``x0 = -1.0``)
was tried first for completeness/comparison and -- consistent with
Anderson & Lo's OWN documented finding for the analogous Jupiter-Europa
attempt (`#753` module docstring, item 1: "converges to trivial... orbit") --
does NOT even produce a valid ``ydot0`` at ``C = 3.010000`` (the Jacobi
constraint's radicand is negative there); this is an expected, DOCUMENTED
negative, not a bug (see the module test suite).

SOURCED CONSTANTS (verified directly against the PDF text layer AND the
rendered page image this task, not inherited from `#764`'s scoping note
without re-checking):

* :data:`VAQUERO_MU` = 2.3658e-4 -- the thesis's own stated Saturn-Titan
  mass ratio (p.132, in the Earth-Moon-system comparison paragraph: "the mass
  parameter is larger in the Earth-Moon system (mu ~ 0.0122) than in the
  Saturn-Titan (mu ~ 2.3658 x 10^-4) and Jupiter-Europa (mu ~ 2.5266 x
  10^-5) systems"). NOTE: unlike Anderson & Lo's UNQUALIFIED 14-digit mu,
  Vaquero's own text uses "~" (approximately) before this value -- an
  explicit signal the displayed 5-significant-figure number is itself a
  rounded DISPLAY of whatever higher-precision mu her own code used
  internally (see the "mu precision" honest finding below). This project's
  own DE440 registry Saturn-Titan mu (``core.satellites`` registry,
  GM_Titan / (GM_Saturn-system + GM_Titan) = 8978.14 / 37940185.14 =
  2.36639e-4) differs from the thesis's rounded display value by ~0.03%
  relative -- the same class of GM-vintage delta the Jovian module documents
  for its own 0.034% Jupiter-Europa gap. This module uses the THESIS's own
  displayed value (2.3658e-4), not the registry value, for the same
  fair-reproduction reason the Jovian module gives.
* :data:`VAQUERO_C` = 3.010000 -- Table 4.1's own stated Jacobi constant
  (p.109 caption, verbatim), the energy at which all four rows are
  evaluated.
* :data:`TABLE41_TARGETS` -- the four rows of Table 4.1's own
  "Unstable Eigenvalue" (lambda_u) column (p.109, verbatim): 3:4 resonant,
  6:5 resonant, L1 Lyapunov, L2 Lyapunov.
* :data:`TABLE41_DIMENSIONAL` -- Table 4.1's own (x [km], ydot [km/s],
  T [days]) columns, verbatim, for the same four rows -- the sourced seeds
  this module's corrector runs use directly (see
  :func:`_table41_seed_nondim`), NOT a blind two-body construction.

CHARACTERISTIC-QUANTITY DERIVATION (l*, t*) -- the thesis states the mass
ratio but does NOT print a numeric l*/t* table for the Saturn-Titan system
(only the general Ch.2 formulas, l* = separation of the primaries, t* =
sqrt(l*^3 / G(m1+m2))); this module derives them from this project's own
DE440 registry (Titan's own SMA about Saturn, ``core.satellites`` registry,
1,221,870 km, as l*; GM_Saturn-system + GM_Titan as the two-body GM sum,
mirroring EXACTLY how :func:`~cyclerfinder.search.jovian_resonant_families.jupiter_europa_system`
derives its own ``t_s`` for period-in-days reporting) -- and then
SELF-VALIDATES that choice: nondimensionalizing all four of Table 4.1's own
printed (x, ydot) pairs with this l*/t* and the thesis's own mu reproduces
the STATED Jacobi constant, C = 3.010000, to 1.6e-5 relative or better for
EVERY row (see :func:`test_table41_ic_reproduces_stated_jacobi_constant` in
the module test suite) -- strong, independent evidence this module's l*/t*
choice matches whatever the thesis itself used internally, not merely a
plausible-looking guess.

HONEST FINDINGS (load-bearing, read before trusting any single row):

1. **3:4, L1 confirm to near-machine precision** (eigenvalue rel_err
   1.07e-6 and 4.56e-6 respectively) -- both the dimensionless eigenvalue
   AND the dimensional (x, ydot, T) columns reproduce Table 4.1 to <0.02%
   relative on every axis. This is strong evidence :data:`VAQUERO_MU` and
   this module's derived l*/t* are correct to within the thesis's own
   5-significant-figure mu display.
2. **6:5 is a genuine, real, SMALL near-miss on eigenvalue only**
   (rel_err = 2.34e-3, just outside :data:`TABLE41_EIGENVALUE_GATE_REL_TOL`)
   despite its OWN dimensional IC and period matching Table 4.1 to <0.02%
   relative, as tight as 3:4/L1. `#769` (see
   ``docs/notes/2026-08-05-769-saturn-titan-65-eigenvalue.md``) followed up
   on this and found the mu-precision mechanism ORIGINALLY proposed here
   (an earlier draft of this finding) does NOT hold up under closer
   measurement, and identified a better-supported alternative:

   - **mu is EXCLUDED as the cause.** The measured +0.1% mu perturbation
     shifts ALL FOUR rows' eigenvalues by a comparable RELATIVE amount
     (+0.096% for 3:4, -0.109% for 6:5, -0.093% for L1, -0.096% for L2 --
     NOT a comparable absolute amount as an earlier draft of this finding
     claimed; absolute shifts instead scale with each row's own eigenvalue
     magnitude). Note 3:4 shifts in the OPPOSITE direction from 6:5/L1/L2.
     Because relative mu-sensitivity is comparable across rows, a mu
     precision floor would produce a comparable RELATIVE eigenvalue error
     in all four rows -- not the observed 1e-6-level match in three rows
     and a 2.3e-3 miss in the fourth. Confirmed directly: bisecting mu to
     hit 6:5's target eigenvalue (191.641) exactly requires mu =
     2.36068e-4, a -0.216% shift from the thesis's stated 2.3658e-4 --
     three orders of magnitude larger than the +0.1% sensitivity probe --
     and at that mu, 3:4/L1/L2 are all pushed OUT of gate to ~2.1e-3 each
     (worse than baseline). No single mu value reconciles all four rows.
   - **A single global C (Jacobi constant) offset does NOT explain the
     table either -- adversarially reviewed and REFUTED, not just
     "not adopted."** An earlier draft of this finding claimed C was "a
     much better-supported explanation" than mu, citing a bisected DeltaC
     of 1.5e-6 relative that reconciles 6:5's eigenvalue while leaving all
     four rows inside the 1e-3 gate. A Fable adversarial review (2026-08-05)
     found that framing overclaimed on two counts: (a) the "~1569 for 6:5
     vs ~3-452 for the others" sensitivity coefficients mean 6:5 is only
     ~3.5x more C-sensitive than L1/L2 specifically (1569/452, 1569/446),
     not "~350-500x" as an earlier draft of this note claimed for all three
     -- that 500x figure only holds against 3:4, whose own coefficient is
     itself ill-determined (a ~45% spread between +/- probes). (b)
     Inverting each row's OWN baseline eigenvalue error through its OWN
     measured C-coefficient gives the implied (C_true - C_printed)/C each
     row would need: 6:5 wants -1.49e-6, but L1 wants only -1.0e-8 and L2
     only -6.1e-9 -- L1/L2 pin Vaquero's effective C to the printed
     3.010000 to ~1e-8 relative, ~150x tighter than what 6:5 would need. A
     single shared C offset is excluded by the SAME logic (and by a wider
     margin) that excluded mu: at the bisected C_fit, L1/L2's own
     eigenvalue errors get 146x/242x WORSE (4.6e-6/2.7e-6 -> 6.66e-4/6.59e-4)
     -- "still inside the 1e-3 gate" is a threshold artifact, not evidence
     FOR a shared-C explanation. Display-rounding cannot rescue it either:
     a correctly-rounded ``C = 3.010000`` (6 dp) bounds any true C within
     5e-7 absolute (1.66e-7 relative) of the printed value, which at 6:5's
     own C-coefficient explains at most ~2.6e-4 of eigenvalue error --
     an order of magnitude short of the observed 2.34e-3 miss. (The
     earlier draft's "plausibility" argument -- that the required DeltaC
     sits within this module's own C self-validation floor,
     ``test_table41_ic_reproduces_stated_jacobi_constant`` -- was also a
     category error: that residual measures THIS module's own l*/t*
     round-trip noise reproducing the PRINTED C, not any evidence about
     whether Vaquero's own internal C differed from what she printed; it
     never even enters the eigenvalue corrector, which holds C exactly at
     the printed value throughout.)
   - **What the evidence actually supports**: 6:5's own eigenvalue genuinely
     is far more C-sensitive than the other three rows' (confirmed smooth
     and monotone over a wide C window, not a fold/bifurcation artifact or
     bisection noise) -- but no single shared parameter (mu, C, or a joint
     mu+C offset checked against the full sensitivity matrix) reconciles
     all four rows simultaneously, and this row's own numerics are
     independently verified rock-solid on our side (Newton residual 3.6e-12,
     eigenvalue stable to 1.8e-10 across a 100x rtol sweep, Barden vs
     Floquet agree to 7e-10). The best-supported honest account is that the
     miss is ROW-SPECIFIC and SOURCE-SIDE -- a transcription/precision
     issue specific to 6:5's own printed 191.641 (this table already
     carries one other demonstrated erratum, L2's period, see finding 3
     below) amplified through 6:5's own unusually steep C-sensitivity --
     not a coherent single-parameter correction. This is reported HONESTLY
     as a FAIL under :data:`TABLE41_EIGENVALUE_GATE_REL_TOL` -- not fudged,
     not silently loosened -- while noting it is a well-characterized,
     small, genuinely-investigated miss, not a wild or unexplained one.
3. **L2's own printed period (T = 79.7260 days) is very likely a
   transcription/typesetting error in the thesis's own Table 4.1**, flagged
   per this project's respectful-errata-framing discipline (evidence-first,
   falsifiable, benefit of the doubt -- typesetting slips happen to
   everyone). Evidence: the SAME seed that reproduces L2's own printed
   ``x`` (km) to 0.003% relative and its own printed lambda_u to 2.7e-6
   relative converges to a self-consistent period of 3.3899 nondim time
   (8.603 days) -- a factor ~9.27 different from the printed 79.7260 days,
   NOT explainable by this module's l*/t* choice (which independently
   reproduces 3:4/6:5/L1's own printed T-in-days columns to <0.02%
   relative using the exact same l*/t*). Fig. 4.6(b) (p.109) itself plots
   the L1 and L2 Lyapunov orbits as small, SIMPLE single-loop ovals of
   comparable size -- consistent with an L2 period on the same O(1)-to-O(10)
   day order of magnitude as L1's own 8.2829-day period, not an order of
   magnitude larger. This module does NOT silently substitute a "corrected"
   period into :data:`TABLE41_DIMENSIONAL` (the printed value stays exactly
   as printed, sourced-only discipline) -- it reports the period criterion
   as an honest, explained FAIL for this one row only, with the full
   evidence in :data:`TABLE41_L2_PERIOD_ERRATA_NOTE`.
4. **is_real_unstable is True (a genuine real saddle eigenvalue, not a
   complex unit-modulus pair) for all four rows** at the eigenvalues found
   here -- unlike the Jovian module's own 5:6-LI search (which needed the
   Barden-vs-``_planar_floquet`` degenerate-eigenvalue cross-check to avoid
   a false "barely unstable" read near |lambda| ~ 1), none of Table 4.1's
   four rows sit anywhere near |lambda| ~ 1 (the smallest is 6:5's 191.6),
   so that specific pitfall does not bite here -- but the cross-check is
   still performed on every row (Barden vs ``_planar_floquet`` agree to
   <1e-8 relative in every case; see the module test suite) as standing
   discipline, not because this module found a new degenerate case.

CONNECTION STAGE EXPLICITLY OUT OF SCOPE for `#765` (its own Task-B analog,
a later task): Vaquero 2013 Sec. 4.3.1 also describes -- as figures + prose
only, no state tables -- a homoclinic connection of the 3:4 resonant orbit
(Fig. 4.9), a periodic "resonant chain" cycling indefinitely between the 3:4
and 6:5 resonances (Fig. 4.10, continued in C in Fig. 4.11), and a
falsifiable published termination claim ("it is suspected that this family
of periodic resonant chains ends for a value of Jacobi constant C <
3.01400", Fig. 4.12). None of that is attempted here.

Pure: math/numpy + :mod:`cyclerfinder.core.cr3bp`,
:mod:`cyclerfinder.search.cr3bp_periodic`,
:mod:`cyclerfinder.search.jovian_resonant_families` (reused directly, not
reimplemented -- ``converge_candidate``, ``ResonantFamilyCandidate``,
``two_body_resonant_seed``, ``basin_robustness_scan``, ``survey_candidates``
all already take a ``system`` argument and are Jupiter-Europa-agnostic per
`#764`'s own confirmed finding).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.jovian_resonant_families as jrf

# ---------------------------------------------------------------------------
# Sourced constants (verified directly against the PDF text layer AND the
# rendered page image of printed p.109 / PDF page 124, #765).
# ---------------------------------------------------------------------------

#: Saturn-Titan mass ratio, Vaquero 2013 p.132 (verbatim, her own "~" prefix
#: preserved in spirit -- see module docstring's mu-precision discussion).
VAQUERO_MU = 2.3658e-4

#: Jacobi constant Table 4.1's four rows are evaluated at (p.109, verbatim).
VAQUERO_C = 3.010000

#: Table 4.1 (p.109), "Unstable Eigenvalue" (lambda_u) column -- verbatim.
TABLE41_TARGETS: dict[str, float] = {
    "3:4": 2129.81,
    "6:5": 191.641,
    "L1": 1004.72,
    "L2": 892.850,
}

#: Table 4.1 (p.109), (x [km], ydot [km/s], T [days]) columns -- verbatim.
#: NOTE (see module docstring finding 3): the "L2" row's T = 79.7260 days is
#: flagged as a likely source transcription error; kept here EXACTLY as
#: printed (sourced-only discipline) -- never silently "corrected".
TABLE41_DIMENSIONAL: dict[str, tuple[float, float, float]] = {
    "3:4": (1.25869e6, 0.477301, 66.3312),
    "6:5": (1.14214e6, 0.545759, 71.2638),
    "L1": (1.15897e6, 0.447315, 8.2829),
    "L2": (1.25231e6, 0.549329, 79.7260),
}

#: Half-crossing index each row's own perpendicular x-axis return sits at,
#: pinned explicitly (not left to auto-selection) for reproducibility -- see
#: module docstring. Determined empirically this task by direct inspection
#: of each seed's own x-axis-crossing list at ``VAQUERO_C`` (module test
#: suite covers this): 3:4 is a genuinely convoluted multi-crossing orbit
#: (matching Fig. 4.4's "flower"-like shape); 6:5 and both Lyapunov orbits
#: are simple single-loop orbits whose OWN first perpendicular return (index
#: 1) is the correct target -- auto-selection (nearest running T/2 estimate)
#: mis-picks a later, WRONG crossing for L2 specifically because its own
#: printed T (79.7260 days, see finding 3) badly overestimates the true
#: half-period, seeding auto-selection toward the wrong crossing index.
_HALF_CROSSINGS: dict[str, int] = {"3:4": 2, "6:5": 1, "L1": 1, "L2": 1}

#: Titan's own SMA about Saturn (km, JPL SSD, ``core.satellites`` registry)
#: -- used as this module's characteristic length l*, mirroring
#: :func:`~cyclerfinder.search.jovian_resonant_families.jupiter_europa_system`'s
#: own use of Europa's SMA. Self-validated against Table 4.1's own four rows
#: (see module docstring) -- NOT an independent free parameter.
_TITAN_SMA_KM = 1221870.0

#: Primary criterion tolerance: relative error on the dimensionless unstable
#: eigenvalue |lambda_u| a recovered candidate must beat to count as a
#: reproduction of a Table 4.1 row. Same class of justification as the
#: Jovian module's own ``TABLE1_GATE_REL_TOL``: the corrector's own
#: convergence floor here is 1e-12 to 1e-14 (crossing residuals, module test
#: suite) -- eight to ten orders of magnitude below this tolerance -- and
#: THREE of Table 4.1's four rows (3:4, L1, L2) reproduce to 1e-6-1e-5
#: relative at exactly the thesis's OWN stated mu, strong evidence this
#: tolerance is generous relative to the real precision floor, not a
#: fudge-to-pass threshold. 1e-3 is kept IDENTICAL to the Jovian module's
#: own value for direct comparability across the two systems' gates; unlike
#: that module, one row here (6:5) genuinely, honestly FAILS this tolerance
#: (2.34e-3, see module docstring finding 2) -- not silently loosened to
#: absorb it.
TABLE41_EIGENVALUE_GATE_REL_TOL = 1e-3

#: Secondary criterion tolerance: relative error on the dimensional period
#: (T, days) a recovered candidate must beat. Justified by the OBSERVED
#: match quality of 3:4/6:5/L1 (0.017%, 0.020%, 0.018% respectively, all
#: roughly an order of magnitude inside this tolerance) -- unlike the
#: Jovian module's abstract "does the period land on a clean 2*pi*q
#: multiple" criterion (meaningful only for its own two-body-seed
#: construction), this module compares directly against the SOURCE's own
#: printed T-in-days column, a much more direct reproduction test made
#: possible by Vaquero printing IC tables Anderson & Lo never did. L2 FAILS
#: this criterion by a factor of ~9.27 (see module docstring finding 3) --
#: an honest, explained, un-fudged FAIL, not evidence against this
#: tolerance choice.
TABLE41_PERIOD_GATE_REL_TOL = 1e-3

#: Tertiary/evidentiary tolerance on the dimensional (x, ydot) IC match --
#: NOT part of the ``passed`` dual criterion (unit-conversion-dependent, per
#: the #765 dispatch note's own framing), reported purely as supporting
#: evidence. All four rows beat this comfortably (<=0.014% relative on both
#: x and ydot -- module test suite).
TABLE41_IC_EVIDENCE_REL_TOL = 1e-3

#: Documents the L2-period anomaly (module docstring finding 3) in one
#: importable, testable place -- not just prose. ``True`` iff this module's
#: own re-derivation independently confirms the thesis's OWN printed L2 row
#: T-in-days value looks inconsistent with its own IC + eigenvalue.
TABLE41_L2_PERIOD_ERRATA_NOTE = (
    "Vaquero 2013 Table 4.1 (p.109) prints the L2 Lyapunov orbit's period as "
    "T = 79.7260 days. This module's own re-derivation -- from the SAME "
    "row's own printed x=1.25231e6 km, ydot=0.549329 km/s, at the SAME "
    "mu=2.3658e-4, C=3.010000, l*/t* this module independently validates "
    "against the OTHER three rows -- converges to a self-consistent period "
    "of 3.3899 nondim time (8.603 days), a factor ~9.27 different, while "
    "reproducing the row's OWN x (0.003% relative) and lambda_u (2.7e-6 "
    "relative) essentially exactly. Fig. 4.6(b) (p.109) plots the L2 orbit "
    "as a small, simple single-loop oval comparable in size to the L1 "
    "orbit (8.2829-day period) -- physically consistent with an ~8.6-day "
    "period, not a ~80-day one. Flagged as a likely source transcription "
    "error (respectful-errata framing: evidence-first, falsifiable, benefit "
    "of the doubt), NOT corrected in TABLE41_DIMENSIONAL (sourced-only "
    "discipline) -- the period gate for this row honestly reports FAIL."
)


def saturn_titan_system(mu: float = VAQUERO_MU) -> cr3bp.CR3BPSystem:
    """Saturn-Titan planar CR3BP system at Vaquero's own mass ratio.

    ``t_s`` is derived from the registry's Saturn-system + Titan GM
    (mirroring :func:`~cyclerfinder.search.jovian_resonant_families.jupiter_europa_system`'s
    own approximation) and ``l_km`` from Titan's own registry SMA; both are
    used ONLY to report periods/positions in physical units and are
    self-validated against Table 4.1's own four rows (see module
    docstring) -- they never enter the nondimensional dynamics or the
    primary eigenvalue gate.
    """
    from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

    gm_pair = PRIMARIES["Saturn"] + SATELLITES["Titan"].mu_km3_s2
    t_s = float(math.sqrt(_TITAN_SMA_KM**3 / gm_pair))
    return cr3bp.CR3BPSystem(
        mu=mu, primary="Saturn", secondary="Titan", l_km=_TITAN_SMA_KM, t_s=t_s
    )


def _table41_seed_nondim(label: str, system: cr3bp.CR3BPSystem) -> tuple[float, float]:
    """Nondimensionalize Table 4.1's own printed (x, T) for ``label`` using
    ``system``'s own ``l_km``/``t_s``. Returns ``(x0_guess, period_guess)``.

    ``ydot0`` is NOT independently nondimensionalized here -- the corrector
    (:func:`~cyclerfinder.search.cr3bp_periodic.correct_symmetric_fixed_jacobi`,
    invoked via :func:`~cyclerfinder.search.jovian_resonant_families.converge_candidate`)
    re-derives ``ydot0`` from ``(x0, jacobi, sign)`` each iteration; every one
    of Table 4.1's four rows has a POSITIVE printed ``ydot`` (km/s), which
    nondimensionalizes to a positive value for every row (verified in the
    module test suite), so ``ydot0_sign=+1`` is used uniformly.
    """
    if label not in TABLE41_DIMENSIONAL:
        raise ValueError(f"label must be one of {sorted(TABLE41_DIMENSIONAL)}; got {label!r}")
    x_km, _ydot_kms, t_days = TABLE41_DIMENSIONAL[label]
    x0_guess = x_km / system.l_km
    period_guess = t_days * 86400.0 / system.t_s
    return x0_guess, period_guess


def recover_table41_candidate(
    label: str, system: cr3bp.CR3BPSystem | None = None
) -> jrf.ResonantFamilyCandidate:
    """Converge+classify Table 4.1's own sourced seed for ``label``.

    Reuses :func:`~cyclerfinder.search.jovian_resonant_families.converge_candidate`
    directly (the corrector, STM/Floquet propagation, and Barden-vs-
    ``_planar_floquet`` cross-check are all system-agnostic, per `#764`'s
    own confirmed finding -- nothing here is reimplemented).

    Raises ``ValueError`` if the corrector fails to converge (should not
    happen for these sourced seeds; the module test suite covers this as a
    standing regression check).
    """
    sys_ = system if system is not None else saturn_titan_system()
    x0_guess, period_guess = _table41_seed_nondim(label, sys_)
    hc = _HALF_CROSSINGS[label]
    cand = jrf.converge_candidate(
        sys_, label, x0_guess, VAQUERO_C, period_guess, ydot0_sign=1.0, half_crossings=hc
    )
    if cand is None:
        raise ValueError(f"sourced Table 4.1 seed for {label!r} failed to converge -- regression")
    return cand


@dataclass(frozen=True)
class Table41GateRow:
    """One row of the Table 4.1 gate report.

    ``passed`` is the PRIMARY+SECONDARY dual criterion (eigenvalue AND
    period), mirroring the Jovian module's own ``GateRow.passed`` dual-
    criterion discipline -- an eigenvalue-only match is NOT sufficient. The
    dimensional (x, ydot) IC match is reported as ``ic_rel_err``/
    ``ic_confirmed`` for evidentiary purposes only (unit-conversion-
    dependent, per the #765 dispatch note's own framing) and does NOT gate
    ``passed`` -- so a row like L2 (near-exact eigenvalue AND IC match, but
    a badly mismatched printed period, see :data:`TABLE41_L2_PERIOD_ERRATA_NOTE`)
    is visible as a well-evidenced eigenvalue confirmation with an honestly
    failing period criterion, rather than collapsing into an
    undifferentiated ``passed=False``.
    """

    label: str
    target_eigenvalue: float
    recovered_eigenvalue: float
    is_real_unstable: bool
    eigenvalue_rel_err: float
    eigenvalue_confirmed: bool
    target_period_days: float
    recovered_period_days: float
    period_rel_err: float
    period_confirmed: bool
    ic_rel_err: float
    ic_confirmed: bool
    passed: bool


def gate_report(system: cr3bp.CR3BPSystem | None = None) -> list[Table41GateRow]:
    """Recover all four Table 4.1 candidates and check each against
    :data:`TABLE41_TARGETS` at :data:`TABLE41_EIGENVALUE_GATE_REL_TOL`, AND
    against the row's own printed period (days) at
    :data:`TABLE41_PERIOD_GATE_REL_TOL` (the dual criterion -- see
    :class:`Table41GateRow`'s docstring). Honest: does not fudge or hide a
    failing row on either criterion. Includes the L1/L2 Lyapunov rows as
    cheap positive controls through the SAME classification path, per the
    #765 dispatch note.
    """
    sys_ = system if system is not None else saturn_titan_system()
    rows: list[Table41GateRow] = []
    for label, target_lambda in TABLE41_TARGETS.items():
        cand = recover_table41_candidate(label, sys_)
        eig_rel_err = abs(cand.max_eigenvalue - target_lambda) / target_lambda
        eig_confirmed = eig_rel_err < TABLE41_EIGENVALUE_GATE_REL_TOL

        x_km_target, ydot_kms_target, t_days_target = TABLE41_DIMENSIONAL[label]
        period_days_recovered = cand.period * sys_.t_s / 86400.0
        period_rel_err = abs(period_days_recovered - t_days_target) / t_days_target
        period_confirmed = period_rel_err < TABLE41_PERIOD_GATE_REL_TOL

        x_km_recovered = cand.x0 * sys_.l_km
        vel_star = sys_.l_km / sys_.t_s
        ydot_kms_recovered = cand.ydot0 * vel_star
        x_rel_err = abs(x_km_recovered - x_km_target) / x_km_target
        ydot_rel_err = abs(ydot_kms_recovered - ydot_kms_target) / ydot_kms_target
        ic_rel_err = max(x_rel_err, ydot_rel_err)
        ic_confirmed = ic_rel_err < TABLE41_IC_EVIDENCE_REL_TOL

        rows.append(
            Table41GateRow(
                label=label,
                target_eigenvalue=target_lambda,
                recovered_eigenvalue=cand.max_eigenvalue,
                is_real_unstable=cand.is_real_unstable,
                eigenvalue_rel_err=eig_rel_err,
                eigenvalue_confirmed=eig_confirmed,
                target_period_days=t_days_target,
                recovered_period_days=period_days_recovered,
                period_rel_err=period_rel_err,
                period_confirmed=period_confirmed,
                ic_rel_err=ic_rel_err,
                ic_confirmed=ic_confirmed,
                passed=eig_confirmed and period_confirmed,
            )
        )
    return rows


__all__ = [
    "TABLE41_DIMENSIONAL",
    "TABLE41_EIGENVALUE_GATE_REL_TOL",
    "TABLE41_IC_EVIDENCE_REL_TOL",
    "TABLE41_L2_PERIOD_ERRATA_NOTE",
    "TABLE41_PERIOD_GATE_REL_TOL",
    "TABLE41_TARGETS",
    "VAQUERO_C",
    "VAQUERO_MU",
    "Table41GateRow",
    "gate_report",
    "recover_table41_candidate",
    "saturn_titan_system",
]
