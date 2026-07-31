"""Neptune-Triton planar CR3BP resonant/Lyapunov/DPO/LPO periodic-orbit
family confirmation gate (`#776`).

Thin sibling of :mod:`cyclerfinder.search.saturn_titan_resonant_families`
(`#765`) and, one level further back, :mod:`cyclerfinder.search.jovian_resonant_families`
(`#753`) -- reusing the SAME corrector/classification/continuation machinery
(no reimplementation) against a third system, confirming `#764`'s own finding
that the pipeline is system-agnostic.

Source (acquired, filed, and independently re-verified this task -- see
``docs/notes/2026-08-01-776-neptune-triton-resonant-families-gate.md``):

    Miceli, G.E. & Bosanac, N. (2026), "Generating Planar Trajectories for
    Neptunian System Exploration Using Motion Primitives," *The Journal of
    the Astronautical Sciences* 73:11, DOI ``10.1007/s40295-025-00545-z``
    (open access, CC-BY 4.0).

ships FOUR Springer Electronic Supplementary Material files (Online
Resources 1-4); Resources 2-4 are machine-readable text files
(``40295_2025_545_MOESM{2,3,4}_ESM.txt``) printing, per row, a 12-decimal
nondimensional planar CR3BP state ``(x0, y0, z0, xdot0, ydot0, zdot0)``, an
"Integration Time" (independently re-verified this task, and already in
`#771`'s own scoping pass, to be the row's own FULL PERIOD -- propagating
each row's IC for its stated Integration Time returns to the IC to
1e-8-1e-11 nondim), and a Jacobi constant, under a mass ratio stated to 16
digits in every file header (``Mass ratio: 2.089503183689124e-04``). UNLIKE
Vaquero 2013 (the Saturn-Titan module's own source, dimensional km/s ICs at
6 significant figures, no printed l*/t* table), this data is ALREADY
nondimensional at full double precision -- so, unlike that module, there is
NO l*/t* sensitivity anywhere in this module's primary gate; ``l_km``/``t_s``
enter ONLY for human-readable period-in-days reporting, using the paper's OWN
stated values (JAS-2026 p.5 body text, independently re-grepped this task):
``l* = 354,760 km``, ``t* ~= 8.081353e4 s``.

PDFs + all 4 ESM files filed at ``cyclers_pdf/papers/`` (see
``docs/notes/CORPUS_INDEX.md``'s own `#776` entry for md5s + the digest
note); the companion CU Boulder PhD dissertation (Miceli 2025, freely
hosted, source of Table 2.1's equilibrium-point positive control) is filed
alongside it.

SOURCED CONSTANTS (re-verified directly against the raw ESM text files and
the two PDFs' own text layers this task, not inherited from `#771`'s
scoping note without re-checking, per
`[[feedback_ground_citations_against_content]]`):

* :data:`MICELI_MU` = 2.089503183689124e-04 -- every ESM file header's own
  stated mass ratio, verbatim (16 digits; both PDFs' own body text round
  this to "mu ~= 0.00020895", consistent).
* :data:`MICELI_L_KM` = 354760.0, :data:`MICELI_T_S` = 80813.53 -- the
  JAS-2026 paper's own stated characteristic length/time (p.5: "l* =
  354,760 km" / "t* ~= 8.081353 x 10^4 s"), used ONLY for period-in-days
  reporting. This project's own DE440-registry Neptune-Triton mu
  (``core.satellites`` registry, GM_Triton / (GM_Neptune-system +
  GM_Triton) = 1428.49546 / 6837955.59604 = 2.08907e-4) differs from the
  paper's own value by ~2.08e-4 relative -- the same GM-vintage-delta class
  the Jovian/Saturn-Titan modules each document for their own systems; this
  module uses the PAPER's own value, not the registry's, for the same
  fair-reproduction reason those modules give. The registry's own Triton SMA
  (354,800 km) also differs from the paper's own ``l*`` by ~1.13e-4
  relative -- immaterial (reporting-only, not gate-load-bearing).
* :data:`ESM_GATE_ROWS` -- ten rows vendored VERBATIM from the raw ESM text
  files (re-downloaded and re-verified this task; md5s recorded in
  ``docs/notes/CORPUS_INDEX.md`` match `#771`'s own scoping-note md5s
  exactly), spanning all three data files and four orbit classes
  (resonant / Lyapunov / DPO / LPO):

  - ``"1:7"`` (ESM3 ``Res17+x+h``, C=1.806962818639) -- the JAS-2026
    Scenario 1 high-energy TARGET orbit (paper Sec. 4.1, "a 1:7 resonant
    orbit with CJ = 1.8 and a period of 41.14 days" -- reproduced below to
    41.137 days, <0.01% relative, via this module's OWN ``l*``/``t*``).
  - ``"3:2-start"`` (ESM4 ``Res32-x+h``, C=3.028835529717) -- the
    JAS-2026 Scenario 2 low-energy START orbit (Sec. 4.2, "a 3:2 resonant
    orbit").
  - ``"L2-lyapunov-target"`` (ESM4 ``L2LyapunovTarget``,
    C=3.003706759619) -- the Scenario 2 low-energy TARGET orbit.
  - ``"L1-lyapunov"`` / ``"L2-lyapunov"`` (ESM2's own worked walkthrough
    pair, C=3.013995763057 / 3.013770826075) -- cheap positive controls
    through the same classification path (mirrors the Saturn-Titan
    module's own L1/L2 inclusion).
  - ``"DPO"`` / ``"LPO"`` (ESM4, C=3.013873926461 / 3.021659863278) --
    the Distant/Low Prograde Orbit families the paper's Sec. 4.2 low-energy
    scenario also builds primitives from.
  - ``"4:5-saddle"`` (ESM4 ``Res45+x+h``, C=2.987089791658) -- the
    Scenario-1 4:5 manifold-source family `#771`'s own eigenvalue survey
    flagged at ``|lambda| ~= 105``.
  - ``"4:7-stress"`` (ESM4 ``Res47-x+h``, first of 3, C=2.997230642137)
    -- the strongest-instability row `#771`'s own survey flagged at
    ``|lambda| ~= 1.5e4``, the deliberate stress row.
  - ``"4:3-saddle"`` (ESM4 ``Res43+x+h``, 2nd of 4, C=3.016635194282) --
    the third genuine saddle-class resonant member (this module's own
    choice, per the `#776` dispatch note's discretion to fill this slot
    from `#771`'s own characterized survey ranges: the FIRST and FOURTH
    printed ``Res43+x+h`` rows are near-unit-circle, NOT saddles -- see
    :data:`FAMILY_43_NEAR_UNIT` below and the honest finding in this
    docstring's "family-mixing" note); this row is ALSO one of the two
    members used for :func:`continue_43_saddle_family`'s own continuation
    confirmation.

  This is a DELIBERATELY NARROW ~10-row vendor, per the `#776` dispatch
  note's own explicit instruction (the user's own "prove narrow first"
  choice) -- NOT the full ~50-row dataset; see the module's own
  ``__all__``-adjacent follow-up-task note and the results note for the
  registered (not-dispatched) follow-up task number for the remaining rows.

* :data:`FAMILY_23` (3 rows, ESM4 ``Res23-x+h``) and :data:`FAMILY_43_HC2`
  (2 rows, ESM4 ``Res43+x+h``, one of which duplicates ``"4:3-saddle"``
  above) -- explicitly called out by the `#776` dispatch note itself
  ("the multi-member families in ESM4 -- e.g. four 4:3 rows, three 2:3
  rows -- give continuation targets at several C_J values each") as the
  two multi-member families :func:`continue_23_family` /
  :func:`continue_43_saddle_family` continue onto. :data:`FAMILY_43_NEAR_UNIT`
  (the OTHER two ``Res43+x+h`` rows) is vendored too, purely to document the
  honest family-mixing finding below -- not part of either continuation gate.

HONEST FINDING (family-mixing, `#776`): all four printed ``Res43+x+h`` rows
are NOT samples of one single-crossing-index continuation branch. Two
(:data:`FAMILY_43_NEAR_UNIT`, C=3.012796850854 and 3.040807590609) converge
under this module's corrector at ``half_crossings=1`` to NEAR-UNIT-CIRCLE
(non-saddle) members; the other two (:data:`FAMILY_43_HC2`, C=3.016635194282
and 3.019162345122) converge at ``half_crossings=2`` to genuine saddles
(``|lambda| ~= 16-20``). :func:`continue_43_saddle_family` (seeded from the
lower-C ``half_crossings=2`` member, walking toward the higher-C one)
succeeds CLEANLY (``StopReason.JACOBI_BOUND``, 51 gauntlet-passing members,
0 rejections, landing on the printed endpoint to 1.0e-4 relative x0 /
9.4e-6 relative period). A separate attempt continuing from the LOWER-C
``half_crossings=1`` member toward the HIGHER-C ``half_crossings=1`` member
(:func:`continue_43_near_unit_family_fold_reversal`) instead hits
``StopReason.FOLD_REVERSAL`` after 95 members, well short of the target C --
the natural-parameter walk on THIS branch turns back in C before reaching
it. This is reported as an honest, well-characterized NEGATIVE (the
"4:3" label spans (at least) two topologically distinct branches at nearby
C, not one smooth curve) -- exactly the kind of genuine partial result this
project's own standing practice treats as a fine, valuable outcome, not a
gate failure to paper over.

HONEST FINDING (two-body seed lineage, `#776`): unlike the Saturn-Titan
module's own naive :func:`~cyclerfinder.search.jovian_resonant_families.two_body_resonant_seed`
attempt (which does not even produce a valid ``ydot0`` at Vaquero's own C),
BOTH of this module's own naive attempts (``two_body_resonant_seed(4, 3,
x0_sign=-1)`` and ``two_body_resonant_seed(2, 3, x0_sign=-1)``, converged
directly at each seed's own natural Jacobi constant) DO converge cleanly
(residuals ~1e-14/1e-15) -- but to a DIFFERENT resonance topology than
either seed's own ``p:q`` label: the "4:3" seed lands on a
``period/2pi ~= 4.0`` orbit (not 3.0, and not matching any vendored 4:3
row's own period), and the "2:3" seed lands on a ``period/2pi ~= 7.0``
orbit (matching the "x7"-family period range, not 3.0). This is the SAME
qualitative finding Anderson & Lo's own paper documents for the analogous
naive Jupiter-Europa attempt ("converges to trivial... orbit") and Vaquero's
own Saturn-Titan case reproduces in a different failure mode -- a third,
independent confirmation that the plain two-body resonant-ellipse
construction is not, by itself, a reliable seed for these systems' own
labeled resonant families (``x0_sign=+1`` is even worse here: DOP853's step
size collapses at the exact secondary-radius crossing, exactly the
documented hazard in :func:`~cyclerfinder.search.jovian_resonant_families.survey_candidates`'s
own docstring -- never attempted as part of the gate, only characterized
once in scratchpad and recorded here to warn future callers).
:data:`TWO_BODY_SEED_LINEAGE_NOTE` documents this in one importable,
testable place.

POSITIVE CONTROL (Miceli 2025 PhD dissertation, Table 2.1, p. ~59,
verbatim, re-verified directly against the dissertation's own text this
task): the collinear/triangular equilibrium-point x-coordinates
(:data:`DISSERTATION_TABLE21_EQUILIBRIA`) reproduce, via
:func:`~cyclerfinder.search.reachable_representatives.lagrange_collinear_x`
(L1/L2, reused directly) plus this module's own tiny L3 root-find (same
``dUbar/dx = 0`` equation, different bracket) and the closed-form L4/L5
``x = 0.5 - mu``, to <5e-7 relative on every point (see
:func:`equilibrium_gate_report`).

CONNECTION/CHAIN/RETROGRADE/CATALOGUE-WRITEBACK STAGE EXPLICITLY OUT OF
SCOPE for `#776` (its own Task-B analog, a later task, mirroring `#765`'s
own Task-A-only scope and the standing `#768`/`#773`/`#775` lesson that
chain-closure work is genuinely hard even with a good source and must never
be bundled into the family task): this module deliberately exposes no
manifold/homoclinic/heteroclinic machinery, does not attempt any of the
paper's own retrograde families (no digit-grade anchor for them -- `#771`
Sec. 1), and makes no catalogue writeback (nothing here is
catalogue-eligible until novelty-checked via
``search/literature_check.py``, a separate future task).

Pure: math/numpy/scipy + :mod:`cyclerfinder.core.cr3bp`,
:mod:`cyclerfinder.search.cr3bp_periodic`,
:mod:`cyclerfinder.search.cr3bp_continuation`,
:mod:`cyclerfinder.search.jovian_resonant_families` (reused directly --
``_classify``, ``ResonantFamilyCandidate``, ``two_body_resonant_seed`` all
already system-agnostic per `#764`'s own confirmed finding),
:mod:`cyclerfinder.search.reachable_representatives` (``lagrange_collinear_x``,
reused directly for L1/L2).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_continuation as cc
import cyclerfinder.search.cr3bp_periodic as cp
import cyclerfinder.search.jovian_resonant_families as jrf
from cyclerfinder.search.reachable_representatives import lagrange_collinear_x

# ---------------------------------------------------------------------------
# Sourced constants (verified directly against the raw ESM text files and the
# JAS-2026 paper's / Miceli 2025 dissertation's own PDF text layers, #776).
# ---------------------------------------------------------------------------

#: Neptune-Triton mass ratio -- verbatim from every ESM file's own header
#: ("Mass ratio: 2.089503183689124e-04"), 16 digits.
MICELI_MU = 2.089503183689124e-04

#: JAS-2026 paper's own stated characteristic length (p.5 body text,
#: verbatim: "l* = 354, 760 km"). Reporting-only -- never enters the
#: nondimensional gate (this data is already nondimensional).
MICELI_L_KM = 354760.0

#: JAS-2026 paper's own stated characteristic time (p.5 body text, verbatim:
#: "t* ~= 8.081353 x 10^4 s"). Reporting-only.
MICELI_T_S = 8.081353e4

#: JAS-2026 paper's own stated characteristic mass (p.5 body text, verbatim:
#: "m* ~= 1.02457 x 10^26 kg"). Not used anywhere in this module; recorded
#: for completeness/documentation only.
MICELI_M_KG = 1.02457e26


@dataclass(frozen=True)
class EsmRow:
    """One vendored row, verbatim from a Miceli & Bosanac 2026 ESM text file.

    ``x0``/``ydot0``/``period`` (the file's own "Integration Time" column,
    independently re-verified to be the full period -- see module
    docstring)/``jacobi`` are copied VERBATIM (12 decimals) from the raw
    file; ``ydot0_sign`` is derived programmatically from ``ydot0``'s own
    sign (never hand-duplicated); ``half_crossings``/``x0_bounds`` are
    this module's OWN corrector-configuration metadata (determined
    empirically this task by direct inspection of each row's own
    crossing-index list -- see module docstring), analogous to the
    Saturn-Titan module's own ``_HALF_CROSSINGS`` table.
    """

    source_file: str  # "ESM2" | "ESM3" | "ESM4"
    source_row_label: str  # the file's own row label, verbatim
    x0: float
    ydot0: float
    period: float
    jacobi: float
    half_crossings: int
    x0_bounds: tuple[float, float] = (-2.0, 2.0)

    @property
    def ydot0_sign(self) -> float:
        return 1.0 if self.ydot0 > 0.0 else -1.0


#: Ten sourced gate rows, verbatim from the raw ESM text files, spanning all
#: three data files and four orbit classes -- see module docstring for the
#: full per-row provenance/rationale. A DELIBERATE ~10-row subset (per the
#: `#776` dispatch note's own "prove narrow first" instruction), NOT the
#: full ~50-row dataset.
ESM_GATE_ROWS: dict[str, EsmRow] = {
    "1:7": EsmRow(
        "ESM3",
        "Res17+x+h",
        7.011583577132,
        -6.902215957345,
        43.981049667607,
        1.806962818639,
        half_crossings=7,
        x0_bounds=(-10.0, 10.0),
    ),
    "3:2-start": EsmRow(
        "ESM4",
        "Res32-x+h",
        -0.901451717937,
        -0.051900595659,
        12.554472960470,
        3.028835529717,
        half_crossings=1,
    ),
    "L2-lyapunov-target": EsmRow(
        "ESM4",
        "L2LyapunovTarget",
        1.056460480396,
        -0.110092617115,
        4.110480863772,
        3.003706759619,
        half_crossings=1,
    ),
    "L1-lyapunov": EsmRow(
        "ESM2",
        "L1Lyapunov",
        0.963347425337,
        -0.026829573549,
        2.968384422151,
        3.013995763057,
        half_crossings=1,
    ),
    "L2-lyapunov": EsmRow(
        "ESM2",
        "L2Lyapunov",
        1.045065938781,
        -0.024387884231,
        3.139197034268,
        3.013770826075,
        half_crossings=1,
    ),
    "DPO": EsmRow(
        "ESM4",
        "DPO",
        1.018236356947,
        0.094629509354,
        1.522062471575,
        3.013873926461,
        half_crossings=1,
    ),
    "LPO": EsmRow(
        "ESM4",
        "LPO",
        1.010469990495,
        0.130296567487,
        0.527889248185,
        3.021659863278,
        half_crossings=1,
    ),
    "4:5-saddle": EsmRow(
        "ESM4",
        "Res45+x+h",
        -0.969056422016,
        -0.126767119783,
        30.398418802755,
        2.987089791658,
        half_crossings=3,
    ),
    "4:7-stress": EsmRow(
        "ESM4",
        "Res47-x+h",
        1.787757094694,
        -1.147924599743,
        44.514684550531,
        2.997230642137,
        half_crossings=5,
        x0_bounds=(-10.0, 10.0),
    ),
    "4:3-saddle": EsmRow(
        "ESM4",
        "Res43+x+h",
        -0.597950460396,
        -0.828492541088,
        12.795970900869,
        3.016635194282,
        half_crossings=2,
    ),
}

#: The 2:3 resonant family (ESM4 "Res23-x+h", all 3 printed members,
#: uniform half_crossings=1 -- unlike the 4:3 family below, this one is a
#: single, cleanly-continuable branch, see :func:`continue_23_family`).
FAMILY_23: tuple[EsmRow, ...] = (
    EsmRow(
        "ESM4",
        "Res23-x+h",
        1.587427060598,
        -0.882468096310,
        18.840967405783,
        3.001357255816,
        half_crossings=1,
    ),
    EsmRow(
        "ESM4",
        "Res23-x+h",
        1.508425262243,
        -0.758195762838,
        18.832770049887,
        3.026732966797,
        half_crossings=1,
    ),
    EsmRow(
        "ESM4",
        "Res23-x+h",
        1.406506637065,
        -0.594708582453,
        18.794508786554,
        3.047064550081,
        half_crossings=1,
    ),
)

#: The genuine-saddle half of the 4:3 family (ESM4 "Res43+x+h", 2nd and 3rd
#: printed rows, half_crossings=2) -- the lower-C member IS
#: ``ESM_GATE_ROWS["4:3-saddle"]``. See :func:`continue_43_saddle_family`.
FAMILY_43_HC2: tuple[EsmRow, ...] = (
    ESM_GATE_ROWS["4:3-saddle"],
    EsmRow(
        "ESM4",
        "Res43+x+h",
        -0.603501155465,
        -0.812257573381,
        12.785730262322,
        3.019162345122,
        half_crossings=2,
    ),
)

#: The near-unit-circle half of the 4:3 family (1st and 4th printed rows,
#: half_crossings=1) -- documents the honest family-mixing finding (module
#: docstring); NOT part of either primary continuation gate. See
#: :func:`continue_43_near_unit_family_fold_reversal`.
FAMILY_43_NEAR_UNIT: tuple[EsmRow, ...] = (
    EsmRow(
        "ESM4",
        "Res43+x+h",
        -0.716813429616,
        -0.540022937644,
        18.817525704645,
        3.012796850854,
        half_crossings=1,
    ),
    EsmRow(
        "ESM4",
        "Res43+x+h",
        -0.792639486296,
        -0.333250768979,
        15.402078897641,
        3.040807590609,
        half_crossings=1,
    ),
)

#: Documents the two-body-resonant-seed-lineage honest negative (module
#: docstring) in one importable, testable place.
TWO_BODY_SEED_LINEAGE_NOTE = (
    "jrf.two_body_resonant_seed(4, 3, x0_sign=-1), converged directly at its "
    "own natural Jacobi constant (2.987670303046456), lands on a "
    "period/2pi ~= 4.0 orbit -- NOT the 4:3 (period/2pi ~= 3.0) topology its "
    "own label implies, and not matching any vendored 4:3 row. "
    "jrf.two_body_resonant_seed(2, 3, x0_sign=-1), converged at its own "
    "natural Jacobi constant (2.9876334118549757), lands on a "
    "period/2pi ~= 7.0 orbit (matching the unrelated 'x7'-family period "
    "range). Both converge cleanly (residuals ~1e-14/1e-15) but to the "
    "WRONG topology -- the same qualitative finding Anderson & Lo 2011 "
    "document for the analogous naive Jupiter-Europa attempt and Vaquero "
    "2013 reproduces in a different failure mode for Saturn-Titan (see "
    "cyclerfinder.search.saturn_titan_resonant_families's own "
    "test_naive_two_body_seed_does_not_converge_at_vaquero_c). "
    "x0_sign=+1 was NOT attempted as part of any gate: it places x0 exactly "
    "at the secondary's own orbital radius, the documented DOP853 "
    "step-collapse hazard in jrf.survey_candidates's own docstring "
    "(confirmed once in scratchpad this task: both (4,3,+1) and (2,3,+1) "
    "time out under a 20s wall-clock budget)."
)

#: Miceli 2025 PhD dissertation Table 2.1 (p.~59), verbatim: equilibrium
#: x-coordinates in the Neptune-Triton CR3BP. y is 0 for L1/L2/L3 and
#: +-sqrt(3)/2 for L4/L5 (not gated here -- x alone is the positive control).
DISSERTATION_TABLE21_EQUILIBRIA: dict[str, float] = {
    "L1": 0.959217,
    "L2": 1.041493,
    "L3": -1.000087,
    "L4": 0.499791,
    "L5": 0.499791,
}

#: Gate tolerance for the equilibrium-point positive control. Justified by
#: the observed match quality (<5e-7 relative on every point, module test
#: suite) -- an order of magnitude inside this tolerance, and by the
#: dissertation's own printed precision (6 decimal places).
EQUILIBRIUM_GATE_REL_TOL = 1e-5

#: (a) Periodicity self-consistency gate: absolute nondim closure residual
#: (planar (x, y, xdot, ydot) norm) a row's own raw IC, propagated for its
#: own stated Integration Time, must beat. Justified against the paper's
#: own stated "converges to 1e-10" multiple-shooting corrector tolerance
#: (Sec. 3, cited by `#771`'s own scoping note) PLUS this module's own
#: independent re-integration (a different integrator, DOP853 at
#: rtol=atol=1e-13, which can amplify residual noise along a strongly
#: unstable direction -- observed up to 9.79e-8 for the 4:7-stress row,
#: |lambda|~1.5e4, an order of magnitude inside this tolerance; every other
#: row is <=2.4e-9).
PERIODICITY_GATE_TOL = 1e-6

#: (b) Reproduction gate: relative error on x0/ydot0/period a candidate
#: re-converged by THIS module's own corrector, seeded at the row's own
#: printed IC, must beat against the row's own printed values. Justified by
#: the observed match quality (module test suite): every one of the ten
#: gate rows reproduces to <=5.8e-11 relative on every one of the three
#: axes -- four+ orders of magnitude inside this tolerance. Tighter than
#: the Saturn-Titan module's own 1e-3 (TABLE41_EIGENVALUE_GATE_REL_TOL)
#: because this data carries NO unit-conversion/rounding loss (already
#: nondimensional to 12 decimals), unlike Vaquero's own dimensional km/s
#: table.
REPRODUCTION_GATE_REL_TOL = 1e-6

#: (c) Barden-vs-planar_floquet INTERNAL cross-check gate (NOT a
#: reproduction -- no published eigenvalue target exists for this system,
#: per `#771`'s own scoping note; see module docstring). Justified by the
#: observed match quality: every row agrees to <=6.5e-7 relative -- more
#: than an order of magnitude inside this tolerance.
CROSSCHECK_GATE_REL_TOL = 1e-5


def neptune_triton_system(mu: float = MICELI_MU) -> cr3bp.CR3BPSystem:
    """Neptune-Triton planar CR3BP system at the source's own mass ratio.

    ``l_km``/``t_s`` are the JAS-2026 paper's OWN stated characteristic
    quantities (module docstring) -- used ONLY to report periods in days;
    unlike the Saturn-Titan/Jovian modules' own registry-derived
    approximations, these are printed directly in the source, not derived.
    """
    return cr3bp.CR3BPSystem(
        mu=mu, primary="Neptune", secondary="Triton", l_km=MICELI_L_KM, t_s=MICELI_T_S
    )


def periodicity_residual(row: EsmRow, system: cr3bp.CR3BPSystem | None = None) -> float:
    """Propagate ``row``'s own raw printed IC for its own stated Integration
    Time (the full period, module docstring) and return the planar
    ``(x, y, xdot, ydot)`` closure residual -- (a) of the `#776` gate,
    data self-consistency independent of this module's own corrector.
    """
    sys_ = system if system is not None else neptune_triton_system()
    state0 = np.array([row.x0, 0.0, 0.0, 0.0, row.ydot0, 0.0])
    arc = cr3bp.propagate(sys_, state0, row.period, with_stm=False, rtol=1e-13, atol=1e-13)
    return float(np.linalg.norm(arc.state_f[:4] - state0[:4]))


def _converge_esm_row(
    row: EsmRow, system: cr3bp.CR3BPSystem, *, label: str = ""
) -> jrf.ResonantFamilyCandidate | None:
    """Converge+classify ``row``'s own sourced seed with THIS module's own
    ``x0_bounds`` (some rows, e.g. "1:7"/"4:7-stress", sit outside
    :func:`~cyclerfinder.search.cr3bp_periodic.correct_symmetric_fixed_jacobi`'s
    own default ``(-2.0, 2.0)`` clamp box -- module docstring). This is a
    thin re-plumb of :func:`~cyclerfinder.search.jovian_resonant_families.converge_candidate`
    (which does not itself expose ``x0_bounds``), reusing
    :func:`~cyclerfinder.search.jovian_resonant_families._classify` and
    :class:`~cyclerfinder.search.jovian_resonant_families.ResonantFamilyCandidate`
    directly -- NOT a reimplementation of the classification logic.
    """
    orbit = cp.correct_symmetric_fixed_jacobi(
        system,
        row.x0,
        row.jacobi,
        row.period,
        ydot0_sign=row.ydot0_sign,
        half_crossings=row.half_crossings,
        tol=1e-11,
        rtol=1e-13,
        atol=1e-13,
        x0_bounds=row.x0_bounds,
    )
    if not orbit.converged:
        return None
    lam_barden, max_eig, is_real, lam_pf = jrf._classify(system, orbit)
    return jrf.ResonantFamilyCandidate(
        label=label,
        x0=orbit.x0,
        ydot0=orbit.ydot0,
        period=orbit.period,
        jacobi=orbit.jacobi,
        ydot0_sign=row.ydot0_sign,
        half_crossings=row.half_crossings,
        crossing_residual=orbit.crossing_residual,
        barden_eigenvalue=lam_barden,
        max_eigenvalue=max_eig,
        is_real_unstable=is_real,
        planar_floquet_eigenvalue=lam_pf,
        period_over_2pi=orbit.period / (2.0 * math.pi),
    )


def recover_esm_candidate(
    label: str, system: cr3bp.CR3BPSystem | None = None
) -> jrf.ResonantFamilyCandidate:
    """Converge+classify :data:`ESM_GATE_ROWS`'s own sourced seed for ``label``.

    Raises ``ValueError`` if the corrector fails to converge (should not
    happen for these sourced seeds; the module test suite covers this as a
    standing regression check).
    """
    sys_ = system if system is not None else neptune_triton_system()
    row = ESM_GATE_ROWS[label]
    cand = _converge_esm_row(row, sys_, label=label)
    if cand is None:
        raise ValueError(f"sourced ESM seed for {label!r} failed to converge -- regression")
    return cand


@dataclass(frozen=True)
class GateRow:
    """One row of the `#776` `#776` gate report -- criteria (a)/(b)/(c) of
    the dispatch note's own dual-plus criterion. ``passed`` requires ALL
    THREE (periodicity self-consistency AND reproduction AND the internal
    cross-check) -- an eigenvalue-only or IC-only match is NOT sufficient,
    mirroring the Jovian/Saturn-Titan modules' own dual-criterion discipline.

    ``max_eigenvalue``/``is_real_unstable`` are reported for evidentiary
    completeness (several rows, e.g. "1:7"/"3:2-start"/"LPO", are genuinely
    near-unit-circle, NOT saddles -- irrelevant to ``passed``, which never
    gates on saddleness; see `#771`'s own scoping note).
    """

    label: str
    source_file: str
    source_row_label: str
    periodicity_residual: float
    periodicity_confirmed: bool
    x0_rel_err: float
    ydot0_rel_err: float
    period_rel_err: float
    reproduction_confirmed: bool
    max_eigenvalue: float
    is_real_unstable: bool
    planar_floquet_eigenvalue: float
    crosscheck_rel_err: float
    crosscheck_confirmed: bool
    passed: bool


def gate_report(system: cr3bp.CR3BPSystem | None = None) -> list[GateRow]:
    """Recover all ten :data:`ESM_GATE_ROWS` candidates and check each
    against criteria (a)/(b)/(c) of the `#776` dispatch note's own gate.
    Honest: does not fudge or hide a failing row on any criterion.
    """
    sys_ = system if system is not None else neptune_triton_system()
    rows: list[GateRow] = []
    for label, row in ESM_GATE_ROWS.items():
        resid = periodicity_residual(row, sys_)
        periodicity_ok = resid < PERIODICITY_GATE_TOL

        cand = _converge_esm_row(row, sys_, label=label)
        if cand is None:
            raise ValueError(f"sourced ESM seed for {label!r} failed to converge -- regression")

        x0_rel_err = abs(cand.x0 - row.x0) / abs(row.x0) if row.x0 != 0 else abs(cand.x0)
        ydot0_rel_err = abs(cand.ydot0 - row.ydot0) / abs(row.ydot0)
        period_rel_err = abs(cand.period - row.period) / abs(row.period)
        reproduction_ok = (
            x0_rel_err < REPRODUCTION_GATE_REL_TOL
            and ydot0_rel_err < REPRODUCTION_GATE_REL_TOL
            and period_rel_err < REPRODUCTION_GATE_REL_TOL
        )

        crosscheck_rel_err = (
            abs(cand.max_eigenvalue - cand.planar_floquet_eigenvalue) / abs(cand.max_eigenvalue)
            if cand.max_eigenvalue != 0
            else abs(cand.planar_floquet_eigenvalue)
        )
        crosscheck_ok = crosscheck_rel_err < CROSSCHECK_GATE_REL_TOL

        rows.append(
            GateRow(
                label=label,
                source_file=row.source_file,
                source_row_label=row.source_row_label,
                periodicity_residual=resid,
                periodicity_confirmed=periodicity_ok,
                x0_rel_err=x0_rel_err,
                ydot0_rel_err=ydot0_rel_err,
                period_rel_err=period_rel_err,
                reproduction_confirmed=reproduction_ok,
                max_eigenvalue=cand.max_eigenvalue,
                is_real_unstable=cand.is_real_unstable,
                planar_floquet_eigenvalue=cand.planar_floquet_eigenvalue,
                crosscheck_rel_err=crosscheck_rel_err,
                crosscheck_confirmed=crosscheck_ok,
                passed=periodicity_ok and reproduction_ok and crosscheck_ok,
            )
        )
    return rows


# ---------------------------------------------------------------------------
# (d) Equilibrium-point positive control (Miceli 2025 dissertation Table 2.1).
# ---------------------------------------------------------------------------


def _lagrange_l3_x(mu: float) -> float:
    """x-coordinate of the collinear L3 point (beyond the primary, opposite
    the secondary) -- the SAME ``dUbar/dx = 0`` equation
    :func:`~cyclerfinder.search.reachable_representatives.lagrange_collinear_x`
    solves for L1/L2, reused directly
    (:func:`cyclerfinder.search.cr3bp_periodic._ubar_grad_x_at_axis`), just
    with L3's own bracket (that function itself only handles L1/L2)."""

    def f(x: float) -> float:
        return cp._ubar_grad_x_at_axis(x, mu)

    return float(brentq(f, -1.5, -mu - 1e-3))


def recovered_equilibrium_x(point: str, mu: float = MICELI_MU) -> float:
    """x-coordinate of collinear/triangular equilibrium ``point`` (L1-L5)."""
    if point in ("L1", "L2"):
        return lagrange_collinear_x(mu, point)
    if point == "L3":
        return _lagrange_l3_x(mu)
    if point in ("L4", "L5"):
        return 0.5 - mu
    raise ValueError(f"point must be one of L1-L5; got {point!r}")


@dataclass(frozen=True)
class EquilibriumGateRow:
    point: str
    target_x: float
    recovered_x: float
    rel_err: float
    confirmed: bool


def equilibrium_gate_report(mu: float = MICELI_MU) -> list[EquilibriumGateRow]:
    """Positive-control gate: this module's own equilibrium-point x against
    :data:`DISSERTATION_TABLE21_EQUILIBRIA`."""
    out: list[EquilibriumGateRow] = []
    for point, target in DISSERTATION_TABLE21_EQUILIBRIA.items():
        recovered = recovered_equilibrium_x(point, mu)
        rel = abs(recovered - target) / abs(target)
        out.append(
            EquilibriumGateRow(
                point=point,
                target_x=target,
                recovered_x=recovered,
                rel_err=rel,
                confirmed=rel < EQUILIBRIUM_GATE_REL_TOL,
            )
        )
    return out


# ---------------------------------------------------------------------------
# Continuation-in-C_J onto the two multi-member families (item (d) of the
# `#776` gate proper: seed-lineage confirmation via continuation).
# ---------------------------------------------------------------------------


def _seed_orbit(row: EsmRow, system: cr3bp.CR3BPSystem) -> cp.SymmetricOrbit:
    orbit = cp.correct_symmetric_fixed_jacobi(
        system,
        row.x0,
        row.jacobi,
        row.period,
        ydot0_sign=row.ydot0_sign,
        half_crossings=row.half_crossings,
        tol=1e-11,
        rtol=1e-13,
        atol=1e-13,
        x0_bounds=row.x0_bounds,
    )
    if not orbit.converged:
        raise ValueError(f"seed row {row.source_row_label!r} failed to converge -- regression")
    return orbit


def continue_23_family(
    system: cr3bp.CR3BPSystem | None = None, *, d_jacobi: float = 2e-4, jacobi_tol: float = 1e-7
) -> cc.FamilyBranch:
    """Continue :data:`FAMILY_23` (ESM4 "Res23-x+h", all 3 members) from
    its own lowest-C printed member toward its own highest-C one, via the
    EXISTING :func:`cyclerfinder.search.cr3bp_continuation.continue_family`
    gauntlet (closure/period-bounds/equilibrium/Jacobi-conservation/
    independent-Radau/dedup/fold-detection all inherited for free).

    FINDING (`#776`): this walk succeeds CLEANLY --
    ``StopReason.JACOBI_BOUND``, 229 gauntlet-passing members, 0 rejections
    -- and the branch's own members land on the 2nd/3rd printed rows to
    6.0e-5/5.7e-7 and 5.7e-4/3.6e-5 relative (x0/period) respectively
    (module test suite): a genuine, tight, independent seed-lineage
    confirmation of the whole 2:3 family.
    """
    sys_ = system if system is not None else neptune_triton_system()
    seed_row = FAMILY_23[0]
    seed = _seed_orbit(seed_row, sys_)
    c_target = FAMILY_23[-1].jacobi
    n_steps = int(abs(c_target - seed_row.jacobi) / abs(d_jacobi)) + 5
    return cc.continue_family(
        sys_,
        seed,
        direction=1,
        d_jacobi=d_jacobi,
        n_steps=n_steps,
        min_jacobi=seed_row.jacobi - 1e-6,
        max_jacobi=c_target + 1e-6,
        half_crossings=seed_row.half_crossings,
        ydot0_sign=seed_row.ydot0_sign,
        seed_label="2:3-family",
        jacobi_tol=jacobi_tol,
    )


def continue_43_saddle_family(
    system: cr3bp.CR3BPSystem | None = None, *, d_jacobi: float = 5e-5, jacobi_tol: float = 1e-7
) -> cc.FamilyBranch:
    """Continue :data:`FAMILY_43_HC2` (the genuine-saddle half of the 4:3
    family, ``half_crossings=2``) from its own lower-C member (=
    ``ESM_GATE_ROWS["4:3-saddle"]``) toward its own higher-C one.

    FINDING (`#776`): succeeds CLEANLY -- ``StopReason.JACOBI_BOUND``, 51
    gauntlet-passing members, 0 rejections -- landing on the printed
    endpoint to 1.0e-4 relative x0 / 9.4e-6 relative period (module test
    suite): a second, independent multi-member family confirmation.
    """
    sys_ = system if system is not None else neptune_triton_system()
    seed_row = FAMILY_43_HC2[0]
    seed = _seed_orbit(seed_row, sys_)
    c_target = FAMILY_43_HC2[-1].jacobi
    n_steps = int(abs(c_target - seed_row.jacobi) / abs(d_jacobi)) + 5
    return cc.continue_family(
        sys_,
        seed,
        direction=1,
        d_jacobi=d_jacobi,
        n_steps=n_steps,
        min_jacobi=seed_row.jacobi - 1e-6,
        max_jacobi=c_target + 1e-6,
        half_crossings=seed_row.half_crossings,
        ydot0_sign=seed_row.ydot0_sign,
        seed_label="4:3-saddle-family",
        jacobi_tol=jacobi_tol,
    )


def continue_43_near_unit_family_fold_reversal(
    system: cr3bp.CR3BPSystem | None = None, *, d_jacobi: float = 2e-4, jacobi_tol: float = 1e-7
) -> cc.FamilyBranch:
    """Continue :data:`FAMILY_43_NEAR_UNIT` (the OTHER half of the 4:3
    family, ``half_crossings=1``) from its own lower-C member toward its
    own higher-C one -- documents the honest family-mixing NEGATIVE (module
    docstring): this walk hits ``StopReason.FOLD_REVERSAL`` well short of
    the target C (module test suite), NOT a clean confirmation. Not part of
    the primary continuation gate (:func:`continue_23_family` /
    :func:`continue_43_saddle_family` already give two clean confirmations)
    -- kept as its own function purely so this genuine, well-characterized
    negative is standing, regression-tested evidence, not just prose.
    """
    sys_ = system if system is not None else neptune_triton_system()
    seed_row = FAMILY_43_NEAR_UNIT[0]
    seed = _seed_orbit(seed_row, sys_)
    c_target = FAMILY_43_NEAR_UNIT[-1].jacobi
    n_steps = int(abs(c_target - seed_row.jacobi) / abs(d_jacobi)) + 5
    return cc.continue_family(
        sys_,
        seed,
        direction=1,
        d_jacobi=d_jacobi,
        n_steps=n_steps,
        min_jacobi=seed_row.jacobi - 1e-6,
        max_jacobi=c_target + 1e-6,
        half_crossings=seed_row.half_crossings,
        ydot0_sign=seed_row.ydot0_sign,
        seed_label="4:3-near-unit-family",
        jacobi_tol=jacobi_tol,
    )


__all__ = [
    "CROSSCHECK_GATE_REL_TOL",
    "DISSERTATION_TABLE21_EQUILIBRIA",
    "EQUILIBRIUM_GATE_REL_TOL",
    "ESM_GATE_ROWS",
    "FAMILY_23",
    "FAMILY_43_HC2",
    "FAMILY_43_NEAR_UNIT",
    "MICELI_L_KM",
    "MICELI_MU",
    "MICELI_M_KG",
    "MICELI_T_S",
    "PERIODICITY_GATE_TOL",
    "REPRODUCTION_GATE_REL_TOL",
    "TWO_BODY_SEED_LINEAGE_NOTE",
    "EquilibriumGateRow",
    "EsmRow",
    "GateRow",
    "continue_23_family",
    "continue_43_near_unit_family_fold_reversal",
    "continue_43_saddle_family",
    "equilibrium_gate_report",
    "gate_report",
    "neptune_triton_system",
    "periodicity_residual",
    "recover_esm_candidate",
    "recovered_equilibrium_x",
]
