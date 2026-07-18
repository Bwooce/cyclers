"""JPL SSD Three-Body Periodic Orbit catalog GATE for raw CR3BP candidates (#647).

``literature_check.py``'s ``KNOWN_CORPUS`` / ``check_literature`` is scoped
entirely to CYCLER-trajectory vocabulary: every generated query is suffixed
"cycler trajectory" / "cycler resonance" / etc., and a hit only counts as a
match when the literal word "cycler"/"cyclic" appears in it. That is the
right tool for adjudicating a candidate CYCLER against the named cycler
literature -- but it is the WRONG tool for adjudicating a raw CR3BP periodic
orbit (halo, DRO, Lyapunov, resonant, ...) against the published catalog of
*those*, because none of that catalog is indexed under "cycler" vocabulary.
#641 hit this by hand: 5 physically distinct Sun-Jupiter periodic-orbit
families all got the SAME spurious "published" verdict citing one cycler
paper, because the keyword search matched on generic orbital-mechanics terms
rather than on any actual family/Jacobi/period identity.

This module is the missing, DIFFERENT-IN-KIND check: given a candidate's
dynamical system, family type, Jacobi constant, period, and mass ratio, query
JPL SSD's own "Three-Body Periodic Orbits" REST API
(https://ssd-api.jpl.nasa.gov/doc/periodic_orbits.html) using ITS OWN family
vocabulary and server-side numeric range filters
(``jacobimin/jacobimax/periodmin/periodmax/stabmin/stabmax``), and return a
structured verdict a human/gauntlet can judge on family-type + invariant
closeness, not a keyword collision.

SCOPE (read before trusting a "not-covered" or "no-match" verdict):

* The live API indexes exactly 7 systems (:data:`~cyclerfinder.verify.
  jpl_periodic_orbits.SUPPORTED_SYSTEMS`) and 12 families (:data:`~
  cyclerfinder.verify.jpl_periodic_orbits.SUPPORTED_FAMILIES`). Sun-Jupiter
  is NOT one of the 7 systems -- an earlier claim that it was (adjacent to
  #641) was a hallucinated/misattributed citation, corrected here and in
  ``jpl_periodic_orbits.py``'s module docstring.
* ``status="not-covered"`` means the SYSTEM or FAMILY is outside the JPL
  catalog's scope -- an honest "this check cannot adjudicate this candidate",
  NOT a "no match found" (which would wrongly suggest the search happened and
  came up empty) and NOT a crash.
* ``status="no-match"`` means the system/family ARE covered, the query ran
  cleanly, but no cataloged member fell inside the tolerance window around
  the candidate's (Jacobi, period) -- i.e. this SPECIFIC candidate is not the
  same orbit as any member JPL currently lists (still NOT proof of novelty:
  the catalog is a finite discretization of each family, not exhaustive).
* ``status="matched"`` means a cataloged member fell inside tolerance --
  treat as a rediscovery of that JPL-catalogued family member.
* ``status="error"`` means the query itself could not be trusted (bad
  parameters, network failure) -- never conflate with "no-match".

DISCIPLINE: this is a NUMERIC-catalog gate, sibling to (not a replacement
for) ``literature_check.check_literature``'s WebSearch-driven cycler-paper
gate. A raw periodic-orbit candidate should be run through BOTH: this module
for "is this exact family member already in JPL's numeric catalog", and
``check_literature`` for "has a NAMED published cycler already claimed this
tour" when the candidate is (or extends into) an actual cycler.

CACHING: every :func:`check_jpl_family` call passes ``cache_dir`` through to
:func:`cyclerfinder.verify.jpl_periodic_orbits.query`, which reads/writes the
raw JSON payload for that exact parameter set under a gitignored ``out/``
directory (this project's existing convention for slow-changing external
astronomical data — see ``.github/workflows/kernel-freshness.yml`` for the
sibling NAIF-kernel pattern). Repeated calls with the same window never
re-hit the live API; be respectful of it during your own testing/exploration
and prefer the cache.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

from cyclerfinder.verify import jpl_periodic_orbits as jpo

#: Verdict for a single :func:`check_jpl_family` call.
FamilyMatchStatus = Literal["matched", "no-match", "not-covered", "error"]

#: Default local cache directory (gitignored ``out/``, per this project's
#: existing external-data-caching convention).
DEFAULT_CACHE_DIR: Path = Path("out") / "jpl_periodic_orbits_cache"

#: Default absolute Jacobi-constant match window (nondim CR3BP units).
DEFAULT_JACOBI_TOL = 0.01

#: Default RELATIVE period match window (fraction of the candidate's period).
DEFAULT_PERIOD_REL_TOL = 0.02

#: How JPL's own catalog should be cited when a match anchors a verdict.
JPL_CATALOG_CITATION = (
    "NASA/JPL Solar System Dynamics, 'Three-Body Periodic Orbits' catalog API "
    "(https://ssd-api.jpl.nasa.gov/doc/periodic_orbits.html), a numeric "
    "reproduction (not a paper) of published CR3BP periodic-orbit families."
)

#: Signature ``jpl_periodic_orbits.query`` matches -- injectable so tests pin
#: a deterministic fixture-backed fake instead of hitting the network.
QueryFn = Callable[..., tuple[jpo.JplSystemConstants, list[jpo.JplOrbit]]]


def normalize_system(system: str) -> str:
    """Map a free-form system label (``"Earth-Moon"``, ``"earth_moon"``, ...)
    onto JPL's own lower-case hyphenated ``sys`` vocabulary."""
    return system.strip().lower().replace("_", "-").replace(" ", "-")


def normalize_family(family: str) -> str:
    """Map a free-form family label onto JPL's own lower-case ``family`` vocabulary."""
    return family.strip().lower().replace("_", "-").replace(" ", "-")


@dataclass(frozen=True)
class JplFamilyMatch:
    """Structured verdict for one candidate's JPL-catalog family check.

    Carries enough detail for a caller to judge match QUALITY, not just a
    boolean: which family/system was queried, how many candidates fell in the
    server-side filter window, and -- for a ``matched``/``no-match`` verdict
    -- how close the best JPL member's Jacobi/period were to the candidate's.
    """

    status: FamilyMatchStatus
    system: str
    """Normalized JPL ``sys`` value queried (or the caller's raw input for a
    ``not-covered`` system verdict)."""
    family: str
    """Normalized JPL ``family`` value queried (or the caller's raw input for
    a ``not-covered`` family verdict)."""
    n_candidates: int = 0
    """How many cataloged members fell inside the server-side filter window."""
    confidence: float = 0.0
    """``[0, 1]``; 0 for not-covered/error/no-match, else a closeness score
    combining the best member's normalized Jacobi + period distance."""
    jacobi_diff: float | None = None
    """``|matched_jacobi - candidate_jacobi|`` for the closest cataloged member."""
    period_diff: float | None = None
    """``|matched_period - candidate_period|`` (nondim TU) for the closest member."""
    matched_jacobi: float | None = None
    matched_period: float | None = None
    matched_stability: float | None = None
    jpl_mu: float | None = None
    """The system mass ratio JPL's response reports (may be requested even
    on a ``no-match`` verdict, since the query still returns system constants
    when candidates exist in a wider band -- ``None`` when the query never
    ran)."""
    mu_reconciliation: dict[str, float] | None = None
    """:func:`~cyclerfinder.verify.jpl_periodic_orbits.reconcile_mu` output
    when the caller supplied its own ``mu`` -- the JPL-vs-ours mass-ratio gap
    that floors how exactly an IC can be expected to re-close."""
    query_params: dict[str, str] = field(default_factory=dict)
    citation: str = ""
    notes: str = ""

    def to_review_block(self) -> dict[str, object]:
        """Render as a machine-written detail block (mirrors
        :meth:`~cyclerfinder.search.literature_check.LiteratureCheckResult.to_review_block`'s
        shape/discipline for the review-queue)."""
        return {
            "status": self.status,
            "system": self.system,
            "family": self.family,
            "n_candidates": self.n_candidates,
            "confidence": self.confidence,
            "jacobi_diff": self.jacobi_diff,
            "period_diff": self.period_diff,
            "matched_jacobi": self.matched_jacobi,
            "matched_period": self.matched_period,
            "matched_stability": self.matched_stability,
            "jpl_mu": self.jpl_mu,
            "mu_reconciliation": self.mu_reconciliation,
            "query_params": dict(self.query_params),
            "citation": self.citation,
            "notes": self.notes,
        }


def _confidence(
    jacobi_diff: float, jacobi_tol: float, period_diff: float, period_tol: float
) -> float:
    """Closeness score in ``[0, 1]``: 1.0 is an exact Jacobi+period match, 0.0
    is at-or-beyond the tolerance window on either axis."""
    j_frac = (
        min(1.0, jacobi_diff / jacobi_tol) if jacobi_tol > 0 else (0.0 if jacobi_diff == 0 else 1.0)
    )
    p_frac = (
        min(1.0, period_diff / period_tol) if period_tol > 0 else (0.0 if period_diff == 0 else 1.0)
    )
    return max(0.0, 1.0 - 0.5 * (j_frac + p_frac))


def check_jpl_family(
    system: str,
    family: str,
    *,
    jacobi: float,
    period: float,
    mu: float | None = None,
    libr: int | None = None,
    branch: str | None = None,
    jacobi_tol: float = DEFAULT_JACOBI_TOL,
    period_rel_tol: float = DEFAULT_PERIOD_REL_TOL,
    cache_dir: Path | str | None = DEFAULT_CACHE_DIR,
    query_fn: QueryFn = jpo.query,
) -> JplFamilyMatch:
    """Check one candidate's (system, family, Jacobi, period) against JPL SSD.

    ``system`` / ``family`` are normalized to JPL's own vocabulary (see
    :func:`normalize_system` / :func:`normalize_family`) before any lookup --
    the whole point of this check versus ``literature_check``'s cycler-
    keyword matcher. Systems/families outside JPL's coverage return an honest
    ``"not-covered"`` verdict WITHOUT any network call.

    ``jacobi_tol`` is an absolute nondim-Jacobi window; ``period_rel_tol`` is
    a RELATIVE period window (fraction of ``period``) -- both are passed to
    JPL as server-side ``jacobimin/jacobimax/periodmin/periodmax`` filters, so
    only the narrow candidate window is ever fetched (never a whole family).

    ``libr``/``branch`` are forwarded when the family requires them (see
    :data:`~cyclerfinder.verify.jpl_periodic_orbits.FAMILIES_REQUIRING_LIBR` /
    ``FAMILIES_REQUIRING_BRANCH``); omitting a required one returns
    ``status="error"`` (a malformed query, not a real "no-match").
    """
    sys_norm = normalize_system(system)
    fam_norm = normalize_family(family)

    if sys_norm not in jpo.SUPPORTED_SYSTEMS:
        return JplFamilyMatch(
            status="not-covered",
            system=sys_norm,
            family=fam_norm,
            notes=(
                f"JPL SSD's Three-Body Periodic Orbits API does not catalog the "
                f"'{sys_norm}' system (supported: {sorted(jpo.SUPPORTED_SYSTEMS)}). "
                "Honest scope gap, not a 'no-match' -- this check cannot adjudicate "
                "this candidate at all."
            ),
        )
    if fam_norm not in jpo.SUPPORTED_FAMILIES:
        return JplFamilyMatch(
            status="not-covered",
            system=sys_norm,
            family=fam_norm,
            notes=(
                f"JPL SSD's Three-Body Periodic Orbits API does not catalog the "
                f"'{fam_norm}' family (supported: {sorted(jpo.SUPPORTED_FAMILIES)}). "
                "Honest scope gap, not a 'no-match'."
            ),
        )
    if fam_norm in jpo.FAMILIES_REQUIRING_LIBR and libr is None:
        return JplFamilyMatch(
            status="error",
            system=sys_norm,
            family=fam_norm,
            notes=f"family '{fam_norm}' requires a libr (libration point) but none was given.",
        )
    if fam_norm in jpo.FAMILIES_REQUIRING_BRANCH and branch is None:
        return JplFamilyMatch(
            status="error",
            system=sys_norm,
            family=fam_norm,
            notes=f"family '{fam_norm}' requires a branch (e.g. 'N'/'S') but none was given.",
        )
    if not math.isfinite(jacobi) or not math.isfinite(period) or period <= 0:
        return JplFamilyMatch(
            status="error",
            system=sys_norm,
            family=fam_norm,
            notes=(
                f"non-physical jacobi={jacobi!r} / period={period!r} "
                "(period must be finite and > 0)."
            ),
        )

    period_min = period * (1.0 - period_rel_tol)
    period_max = period * (1.0 + period_rel_tol)
    query_params = {
        "sys": sys_norm,
        "family": fam_norm,
        "jacobimin": repr(jacobi - jacobi_tol),
        "jacobimax": repr(jacobi + jacobi_tol),
        "periodmin": repr(period_min),
        "periodmax": repr(period_max),
        "periodunits": "TU",
    }
    if libr is not None:
        query_params["libr"] = str(libr)
    if branch is not None:
        query_params["branch"] = branch

    try:
        constants, orbits = query_fn(
            sys_norm,
            fam_norm,
            libr=libr,
            branch=branch,
            jacobimin=jacobi - jacobi_tol,
            jacobimax=jacobi + jacobi_tol,
            periodmin=period_min,
            periodmax=period_max,
            periodunits="TU",
            cache_dir=cache_dir,
        )
    except Exception as exc:  # a flaky/malformed query is an error, not a no-match
        return JplFamilyMatch(
            status="error",
            system=sys_norm,
            family=fam_norm,
            query_params=query_params,
            notes=f"JPL query failed: {exc!r}",
        )

    mu_rec = jpo.reconcile_mu(constants.mu, mu) if mu is not None else None

    if not orbits:
        return JplFamilyMatch(
            status="no-match",
            system=sys_norm,
            family=fam_norm,
            n_candidates=0,
            jpl_mu=constants.mu,
            mu_reconciliation=mu_rec,
            query_params=query_params,
            notes=(
                "Query ran cleanly (system+family are JPL-covered) but no cataloged "
                f"member fell within the tolerance window (jacobi_tol={jacobi_tol}, "
                f"period_rel_tol={period_rel_tol}). NOT proof of novelty -- JPL's "
                "catalog is a finite discretization of each family."
            ),
        )

    def _dist(o: jpo.JplOrbit) -> float:
        j_frac = abs(o.jacobi - jacobi) / jacobi_tol if jacobi_tol > 0 else 0.0
        p_tol = period_rel_tol * period
        p_frac = abs(o.period - period) / p_tol if p_tol > 0 else 0.0
        return j_frac + p_frac

    best = min(orbits, key=_dist)
    jacobi_diff = abs(best.jacobi - jacobi)
    period_diff = abs(best.period - period)
    period_tol_abs = period_rel_tol * period
    confidence = _confidence(jacobi_diff, jacobi_tol, period_diff, period_tol_abs)

    return JplFamilyMatch(
        status="matched",
        system=sys_norm,
        family=fam_norm,
        n_candidates=len(orbits),
        confidence=round(confidence, 3),
        jacobi_diff=jacobi_diff,
        period_diff=period_diff,
        matched_jacobi=best.jacobi,
        matched_period=best.period,
        matched_stability=best.stability,
        jpl_mu=constants.mu,
        mu_reconciliation=mu_rec,
        query_params=query_params,
        citation=JPL_CATALOG_CITATION,
        notes=(
            f"{len(orbits)} cataloged '{sys_norm}' '{fam_norm}' member(s) fell within "
            "the tolerance window; treat the closest as a rediscovery of a JPL-"
            "catalogued family member, NOT novelty-claimable on this axis."
        ),
    )


__all__ = [
    "DEFAULT_CACHE_DIR",
    "DEFAULT_JACOBI_TOL",
    "DEFAULT_PERIOD_REL_TOL",
    "JPL_CATALOG_CITATION",
    "FamilyMatchStatus",
    "JplFamilyMatch",
    "QueryFn",
    "check_jpl_family",
    "normalize_family",
    "normalize_system",
]
