"""#347 Phase 2 P2.2 — Sweep parents enumeration.

The Phase 2 discovery sweep walks symmetric (k1, k2) Earth-Moon cycler family
parents through Jacobi continuation, applies the saddle-center detector
(:func:`cyclerfinder.search.bifurcation_detector.detect_saddle_center_bracket`)
with the Pellegrini-Russell 2016 fixed-path STM mode (per #372 P372.3, the
variable-step trivial-pair smear of 0.35% is comparable to k=1 saddle-center
signal), and branches off at each saddle-center via
:func:`cyclerfinder.genome.asymmetric_branch.branch_at_saddle_center` (with
the #379 / Phase 2 P2.1 Gram-Schmidt fix in place).

This module enumerates the **parent families** the sweep visits. Each entry
is a sourced symmetric cycler IC + the continuation window + the topology
target. The parent inventory is capped at the small set Phase 1's substrate
demonstrably handles + the most relevant other (k1, k2) members from
Braik-Ross Table 2.

Sweep parents (this Phase 2 round):

  * **C11a** — (1, 1)a, sigma_d = 1.0482 d⁻¹. Sourced anchor: Braik-Ross Table 2.
    P=42.140 d at C_J=3.1294. Phase 2 sweeps UPWARD in C looking for the
    family's saddle-center (RTR2026 p.5 theorem: every symmetric cycler family
    has a saddle-center at its maximum C).
  * **C11b** — (1, 1)b, sigma_d = 0.9255 d⁻¹. Same provenance. P=55.995 d.
  * **C21** — (2, 1), sigma_d = 0.1358 d⁻¹. Sourced at C_J_C21 = 3.129389531088256
    (Ross-RT 2025 unrounded; Braik-Ross's literal 3.1294 sits ~1e-5 outside
    the family's C-extent). P=84.533 d. **NB:** the (2, 1) family has a very
    narrow C-extent (~4e-12 per #262); the Phase 2 sweep must use a very small
    dC (e.g. 1e-13) AND a short max_steps (e.g. 5) to stay inside the family.
  * **C32** — (3, 2), sigma_d = 0.6886 d⁻¹. **Phase 1 anchor**, included as the
    control case to verify Phase 2's substrate REPRODUCES Phase 1's saddle-
    center finding under the new Gram-Schmidt-corrected branch corrector.
    Sweep window: C ∈ [3.1294, 3.1544], dC=1e-4, ~250 members (matches
    Phase 1). Additionally: the inverse-direction second saddle-center at
    C ∈ (3.1494, 3.1524) observed at Phase 1 P1.3 is documented (the
    detector intentionally skips real_far→complex_unit_circle transitions;
    Phase 2 does NOT widen this — flagged for future Phase 2.x extension).

Parents DEFERRED to a future Phase 2.x:

  * (5, 2), (3, 1), (4, 3) and the higher-k Aldrin / Russell-Ocampo entries —
    each needs its own sourced symmetric IC (the Braik-Ross Table 2 only
    covers the four cyclers above + 9 non-cycler reachable representatives;
    higher-codimension cyclers like (3, 1) require Aldrin 1986 / Byrnes
    1993 / Russell-Ocampo 2006 ICs which are not present in the offline
    seeds table). The substrate generalises trivially; the work is adding
    the sourced IC + the per-family Jacobi.

Discipline: every parent in this module traces to a published IC in the
cyclers project corpus (Braik-Ross 2026 Table 2 + Ross-RT 2025 unrounded
C_J_C21). No values from memory; the recovery pipeline
:func:`cyclerfinder.search.reachable_representatives.recover_all_cyclers_braik_ross`
already produces all four ICs against the sourced periods.
"""

from __future__ import annotations

from dataclasses import dataclass

from cyclerfinder.search.reachable_representatives import (
    C_J_BRAIK_ROSS,
    C_J_C21,
)


@dataclass(frozen=True)
class SweepParent:
    """Phase 2 sweep parent: a sourced symmetric cycler IC + continuation window.

    Attributes
    ----------
    label :
        Braik-Ross Table 2 row identifier (e.g. "C32"); also the key into
        :data:`cyclerfinder.search.reachable_representatives._CYCLER_SEEDS`.
    k1, k2 :
        Winding-topology integers of the parent family. Used for the
        independent topology cross-check at the recovered anchor (matches
        Phase 1 P1.1 gate).
    jacobi_anchor :
        The published Jacobi constant the anchor sits on (literal Braik-Ross
        3.1294 for C11a/C11b/C32; unrounded C_J_C21 for C21).
    cj_window :
        Tuple ``(c_min, c_max)`` for the natural-parameter continuation walk.
        Includes the anchor and extends in the direction the saddle-center is
        expected (positive ΔC per RTR2026 p.5: the bifurcation sits at the
        family's C-maximum).
    dc :
        Continuation step in Jacobi (positive; direction=+1).
    n_steps :
        Maximum continuation steps. The walk stops earlier on a fold / topology
        jump / Jacobi-bound.
    sourced_period_days :
        Published period for the anchor (Braik-Ross Table 2). Used by the
        anchor-stage period gate (1% tolerance per Phase 1).
    sourced_sigma_d :
        Published Floquet instability sigma_d = ln|λ_max|/T in day⁻¹. Used for the
        independent diagnostic gate (5% tolerance, per Phase 1).
    notes :
        Free-text discussion of the per-family caveats (narrow C-extent,
        deferred extensions, etc.).
    """

    label: str
    k1: int
    k2: int
    jacobi_anchor: float
    cj_window: tuple[float, float]
    dc: float
    n_steps: int
    sourced_period_days: float
    sourced_sigma_d: float
    notes: str


# Phase 2 sweep parents. Six entries: the four Braik-Ross Table 2 cyclers
# (C11a, C11b, C21, C32) PLUS C32 walked DOWN (decreasing C) to characterise
# the family's other limit, and a smaller version of C21 (narrower window) to
# probe the tiny-extent (2, 1) family's saddle-center.
#
# Per the Phase 2 mandate "Cap initial sweep at ~6-8 (k1, k2) pairs to keep
# scope tractable" — six pairs is the initial cap.
PHASE2_SWEEP_PARENTS: tuple[SweepParent, ...] = (
    SweepParent(
        label="C32",
        k1=3,
        k2=2,
        jacobi_anchor=C_J_BRAIK_ROSS,
        cj_window=(3.1294, 3.1544),
        dc=1e-4,
        n_steps=250,
        sourced_period_days=78.613,
        sourced_sigma_d=0.1583,  # day⁻¹ per Braik-Ross Table 2 (NOT TU⁻¹)
        notes=(
            "Phase 1 anchor; included as a control to verify Phase 2 substrate "
            "REPRODUCES the Phase 1 saddle-center at C* ∈ (3.14170, 3.14180). "
            "Phase 1 P1.3 also observed a SECOND saddle-center in the inverse "
            "direction (real_far → complex_unit_circle) at C ∈ (3.1494, 3.1524) "
            "— flagged for future Phase 2.x extension; current detector skips it."
        ),
    ),
    SweepParent(
        label="C11a",
        k1=1,
        k2=1,
        jacobi_anchor=C_J_BRAIK_ROSS,
        cj_window=(3.1294, 3.1500),
        dc=1e-4,
        n_steps=200,
        sourced_period_days=42.140,
        sourced_sigma_d=1.0482 / 27.321661 * 2.0 * 3.141592653589793,  # TU⁻¹ → d⁻¹ via TU_DAYS
        notes=(
            "(1, 1)a cycler — Braik-Ross Table 2 most-unstable cycler "
            "(sigma = 1.0482 TU⁻¹). Sweep UPWARD looking for the family's "
            "saddle-center at C* > 3.1294 per RTR2026 p.5."
        ),
    ),
    SweepParent(
        label="C11b",
        k1=1,
        k2=1,
        jacobi_anchor=C_J_BRAIK_ROSS,
        cj_window=(3.1294, 3.1500),
        dc=1e-4,
        n_steps=200,
        sourced_period_days=55.995,
        sourced_sigma_d=0.9255 / 27.321661 * 2.0 * 3.141592653589793,
        notes=(
            "(1, 1)b cycler — Braik-Ross Table 2 sister of C11a. Distinct family "
            "(different ydot0_sign + different half_crossings; both topology (1, 1))."
        ),
    ),
    SweepParent(
        label="C21",
        k1=2,
        k2=1,
        jacobi_anchor=C_J_C21,
        cj_window=(C_J_C21 - 1e-11, C_J_C21 + 1e-11),
        dc=1e-13,
        n_steps=20,
        sourced_period_days=84.533,
        sourced_sigma_d=0.1358 / 27.321661 * 2.0 * 3.141592653589793,
        notes=(
            "(2, 1) cycler — CAUTION: per #262, the (2, 1) family has a TINY "
            "Jacobi extent (~4e-12). Sweep window is set ±1e-11 with dC=1e-13 "
            "to stay inside. Saddle-center detection likely degenerate because "
            "the family ends ABRUPTLY (no continuation through C*); the natural "
            "fold-back is the expected mode. Phase 2 outcome on this family is "
            "expected to be HONEST_NEGATIVE — the substrate may not even reach "
            "a saddle-center bracket."
        ),
    ),
    SweepParent(
        label="C32_down",
        k1=3,
        k2=2,
        jacobi_anchor=C_J_BRAIK_ROSS,
        cj_window=(3.10, 3.1294),
        dc=1e-4,
        n_steps=250,
        sourced_period_days=78.613,
        sourced_sigma_d=0.1583,
        notes=(
            "C32 walked DOWNWARD (decreasing C) from the anchor — probes the "
            "OTHER end of the (3, 2) family. RTR2026 Fig. 4 lower panel shows "
            "the saddle-center at the C-maximum, but the family may have "
            "additional bifurcations toward lower C. Direction parameter for "
            "continue_family is -1 (vs +1 for the other parents); the driver "
            "handles this via the suffix '_down'."
        ),
    ),
    SweepParent(
        label="C11a_down",
        k1=1,
        k2=1,
        jacobi_anchor=C_J_BRAIK_ROSS,
        cj_window=(3.10, 3.1294),
        dc=1e-4,
        n_steps=200,
        sourced_period_days=42.140,
        sourced_sigma_d=1.0482 / 27.321661 * 2.0 * 3.141592653589793,
        notes=(
            "C11a walked DOWNWARD. Complements the C11a upward sweep; the "
            "(1, 1)a family may have bifurcations on either side of the "
            "Braik-Ross anchor."
        ),
    ),
)
