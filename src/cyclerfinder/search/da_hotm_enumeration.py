"""DA/HOTM enumeration driver + novelty routing (#450 Task 6).

The production driver for the global multi-rev fixed-point enumeration lane. It
loops over a Jacobi band and revolution range, enumerates fixed-point candidates
over the section domain (design draft §3 filter cascade), runs the
reproduction/known-family screen, and classifies each survivor:

    reproduction          -- a known base family (DRO / Lyapunov member, JPL /
                             catalogue / sourced)
    known-family          -- a continuation point of a known family
    novel-PO              -- a genuinely new CR3BP periodic orbit, NOT a cycler
                             (the Png' case: a DRO/Lyapunov bridge, no Earth-Earth
                             resonant transfer leg)
    novel-cycler-candidate-- new AND carries cycler transfer structure
    uncertified           -- the corrector did not close it (kept for audit)

It REUSES the existing modules unchanged: the enumerator, the
``correct_general_periodic`` corrector, and (injectably) the reproduction screen +
literature check. The reproduction screen (``is_known_family``) and the lit-check
are injected so the lane runs offline in tests and wires the live oracles in
production (mirrors ``literature_check``'s injected ``SearchFn``).

Catalogue boundary (design draft §4): a Png'-class PO never earns a catalogue row
on PO-ness alone -- ``novel-PO`` is a discovery-ledger note, not a catalogue row.

Pure: the enumerator + corrector + injected screens.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Literal

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.da_hotm_backend import SamplingSectionMap
from cyclerfinder.genome.da_hotm_enumerator import DomainBox, enumerate_fixed_points
from cyclerfinder.search.cr3bp_general_periodic import correct_general_periodic

Classification = Literal[
    "reproduction",
    "known-family",
    "novel-PO",
    "novel-cycler-candidate",
    "uncertified",
]

# Reproduction screen: (x0, xdot0, n) -> True iff this is a KNOWN family (a JPL /
# catalogue / sourced base family). Injected so tests run offline and production
# wires the live JPL-oracle / catalogue-signature dedup.
KnownFamilyFn = Callable[[float, float, int], bool]


def classify_candidate(*, is_known: bool, has_cycler_structure: bool) -> Classification:
    """Route a CERTIFIED candidate by the reproduction screen + cycler structure.

    * known family                        -> ``reproduction``
    * unknown, no cycler transfer leg     -> ``novel-PO`` (the Png' case)
    * unknown, has a cycler transfer leg  -> ``novel-cycler-candidate``
    """
    if is_known:
        return "reproduction"
    if has_cycler_structure:
        return "novel-cycler-candidate"
    return "novel-PO"


@dataclass(frozen=True)
class EnumerationResult:
    """One classified entry in the enumeration ledger.

    Carries the certified IC (or the coarse one if uncertified), the section
    residual at emission, the certified period/closure (NaN if uncertified), the
    classification, and the ``(C, n)`` provenance.
    """

    x0: float
    xdot0: float
    ydot0: float
    c_target: float
    n: int
    section_residual: float
    period: float
    closure_residual: float
    converged: bool
    classification: Classification


def run_enumeration(
    system: cr3bp.CR3BPSystem,
    *,
    c_band: Sequence[float],
    n_range: Sequence[int],
    domain_box: DomainBox,
    is_known_family: KnownFamilyFn,
    has_cycler_structure: Callable[[float, float, int], bool] | None = None,
    residual_tol: float = 1e-2,
    grid: tuple[int, int] = (41, 31),
    dedup_radius: float = 0.01,
    ydot0_sign: float = 1.0,
    close_tol: float = 1e-11,
) -> list[EnumerationResult]:
    """Run the enumeration lane over a (C, n) band and classify each candidate.

    For each ``(C, n)``: build the sampling section map, enumerate coarse fixed
    points over ``domain_box`` (section-residual gate + dedup), certify each with
    ``correct_general_periodic`` (the existing corrector, unchanged), screen for
    reproduction, and classify. Returns the ledger (possibly empty -- an empty band
    is a legitimate negative, re-stamped in the registry by Task 7).

    ``has_cycler_structure(x0, xdot0, n)`` -- optional; whether a certified orbit
    carries a cycler transfer leg (default: never, since a CR3BP PO is not a cycler
    on PO-ness alone, design draft §4).
    """
    cycler_struct = has_cycler_structure or (lambda x0, xdot0, n: False)
    results: list[EnumerationResult] = []
    for c_target in c_band:
        backend = SamplingSectionMap(system, c_target=float(c_target), ydot_sign=ydot0_sign)
        for n in n_range:
            cands = enumerate_fixed_points(
                backend,
                domain_box,
                int(n),
                residual_tol=residual_tol,
                grid=grid,
                dedup_radius=dedup_radius,
            )
            for cand in cands:
                orbit = correct_general_periodic(
                    system,
                    cand.x0,
                    cand.xdot0,
                    float(c_target),
                    # horizon: ~2n y=0 crossings; a generous nondim guess.
                    period_guess=float(2 * n) * 2.5,
                    half_crossings=2 * int(n),
                    ydot0_sign=ydot0_sign,
                    tol=close_tol,
                )
                if orbit.converged and orbit.residual <= close_tol:
                    x0, xdot0 = orbit.x0, orbit.xdot0
                    is_known = is_known_family(x0, xdot0, int(n))
                    classification: Classification = classify_candidate(
                        is_known=is_known,
                        has_cycler_structure=cycler_struct(x0, xdot0, int(n)),
                    )
                    results.append(
                        EnumerationResult(
                            x0=x0,
                            xdot0=xdot0,
                            ydot0=orbit.ydot0,
                            c_target=float(c_target),
                            n=int(n),
                            section_residual=cand.residual,
                            period=orbit.period,
                            closure_residual=orbit.closure_residual,
                            converged=True,
                            classification=classification,
                        )
                    )
                else:
                    # Coarse candidate that did not certify; kept for audit so a
                    # band's "empty" claim is bounded (it surfaced candidates).
                    is_known = is_known_family(cand.x0, cand.xdot0, int(n))
                    results.append(
                        EnumerationResult(
                            x0=cand.x0,
                            xdot0=cand.xdot0,
                            ydot0=float("nan"),
                            c_target=float(c_target),
                            n=int(n),
                            section_residual=cand.residual,
                            period=float("nan"),
                            closure_residual=float("nan"),
                            converged=False,
                            classification="reproduction" if is_known else "uncertified",
                        )
                    )
    return results


__all__ = [
    "Classification",
    "EnumerationResult",
    "KnownFamilyFn",
    "classify_candidate",
    "run_enumeration",
]
