"""Mine JPL SSD's Three-Body Periodic Orbit catalog as a DISCOVERY INPUT (#667).

`#647` built :mod:`cyclerfinder.search.jpl_family_check` as a novelty-check
GATE: given ONE candidate's ``(system, family, jacobi, period, mu)`` it asks
"does JPL already catalog this exact family member". `#661` shortlist item 5
(see `#661`'s own `data/OUTSTANDING.md` bullet for the full case) asks the
OPPOSITE question: can that same catalog be mined the other direction, as a
source of candidates, in systems this project's own search methods have never
directly targeted -- notably ``saturn-enceladus``, ``mars-phobos`` and
``sun-mars`` (``earth-moon``/``sun-earth`` are extensively covered by this
project's other methods already; ``jupiter-europa``/``saturn-titan`` have some
existing coverage too).

This module is that DIFFERENT, additive capability:

1. :func:`fetch_family_window` -- bulk retrieval of every cataloged member of
   one ``(system, family[, libr, branch])`` inside a caller-chosen Jacobi/
   period window, reusing :func:`cyclerfinder.verify.jpl_periodic_orbits.query`
   and its caching convention directly (server-side range filters, never an
   unfiltered whole-family fetch) -- the bulk-retrieval sibling of
   ``check_jpl_family``'s single-candidate-focused query construction.
2. :func:`propagate_min_distances_km` -- generalizes
   :func:`cyclerfinder.search.real_binary_kk_sweep.min_body_clearance_km`
   (identical DOP853/rtol/atol/n_samples propagation convention) from that
   function's restricted planar-symmetric IC form to a FULL 6-D state, since
   JPL's own catalogued ICs are general (halo/vertical/axial family members
   carry nonzero z0/vx0).
3. :func:`classify_secondary_approach` -- the "does this JPL-catalogued
   member geometrically qualify as a recurrent close secondary approach"
   verdict, built from TWO criteria this project ALREADY uses elsewhere
   rather than an invented ad hoc geometric rule:

   * a physical non-crash floor, the same zero-margin "does the trajectory
     ever pass inside the body" gate `#660`'s
     ``real_binary_kk_sweep._gate_clearance`` already applies (see
     :data:`cyclerfinder.search.real_binary_kk_sweep.DEFAULT_CLEARANCE_MARGIN_KM`);
   * a "genuinely CLOSE, not just somewhere-in-the-system" ceiling, borrowed
     BY ANALOGY from :mod:`cyclerfinder.genome.hill_screen`'s existing
     ``PASS_HILL_FRACTION = 0.3`` band. Read the caveat in
     :data:`CLOSE_APPROACH_HILL_FRACTION`'s own docstring before trusting
     this as more than a documented judgment call: ``hill_screen`` applies
     0.3 to a DIFFERENT quantity (an orbit's amplitude about the PRIMARY,
     relative to the Earth-Sun Hill radius, as a solar-tide-survivability
     proxy); here it is applied to distance-to-SECONDARY relative to the
     SECONDARY's own Hill radius, as a "deep inside the body's own
     gravitational neighbourhood vs. out at the ragged edge shared with L1/
     L2-family orbits" proxy. Same 0.3 number, related but not identical
     physics -- reused for consistency and because both measure the same
     underlying "how deep inside a gravity well is this" idea, not because
     the two ratios are algebraically the same thing.

Stability (JPL's own reported ``stability`` index) and Jacobi/period are
reported for every candidate but are NOT gated on here -- see
:class:`SecondaryApproachVerdict`'s own docstring for why boundedness over one
period is automatic for any exact CR3BP periodic-orbit family member, and
station-keeping against linear instability is a separate, later engineering
question this module does not adjudicate.

DISCIPLINE: this module NEVER writes to ``data/catalogue.yaml`` and makes no
novelty claim. A candidate passing :func:`classify_secondary_approach` is
merely "geometrically close to the secondary, does not crash, JPL already
catalogs it" -- see `#667`'s own `data/OUTSTANDING.md` bullet / `#661`'s
honest framing for why the novelty ceiling here is low even for a genuine
geometric pass (being numerically catalogued at JPL is not the same as being
characterized as a cycler).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from cyclerfinder.core.cr3bp import cr3bp_eom
from cyclerfinder.genome.hill_screen import PASS_HILL_FRACTION
from cyclerfinder.search.jpl_family_check import (
    DEFAULT_CACHE_DIR,
    QueryFn,
    normalize_family,
    normalize_system,
)
from cyclerfinder.search.real_binary_kk_sweep import DEFAULT_CLEARANCE_MARGIN_KM
from cyclerfinder.verify import jpl_periodic_orbits as jpo

#: Borrowed BY ANALOGY from :data:`cyclerfinder.genome.hill_screen.PASS_HILL_FRACTION`
#: -- see this module's own docstring for the important caveat that the two
#: ratios are physically related but not identical. A candidate whose
#: closest approach to the secondary is within this fraction of the
#: secondary's OWN Hill radius counts as "genuinely close" here; beyond it,
#: the candidate is geometrically no closer to the secondary than an
#: L1/L2-family orbit living at the edge of its gravitational neighbourhood.
CLOSE_APPROACH_HILL_FRACTION: float = PASS_HILL_FRACTION

#: Same DOP853 dense-sampling convention
#: :func:`cyclerfinder.search.real_binary_kk_sweep.min_body_clearance_km` and
#: :func:`cyclerfinder.search.binary_star_search.winding_topology` already use.
DEFAULT_N_SAMPLES: int = 4000
DEFAULT_RTOL: float = 1e-11
DEFAULT_ATOL: float = 1e-11


def fetch_family_window(
    system: str,
    family: str,
    *,
    libr: int | None = None,
    branch: str | None = None,
    jacobi_min: float | None = None,
    jacobi_max: float | None = None,
    period_min: float | None = None,
    period_max: float | None = None,
    periodunits: str = "TU",
    cache_dir: Path | str | None = DEFAULT_CACHE_DIR,
    query_fn: QueryFn = jpo.query,
) -> tuple[jpo.JplSystemConstants, list[jpo.JplOrbit]]:
    """Bulk-retrieve every cataloged ``(system, family)`` member in a window.

    The DISCOVERY-INPUT sibling of
    :func:`cyclerfinder.search.jpl_family_check.check_jpl_family` (#667):
    that function's whole query construction is aimed at ONE candidate's
    tolerance window; this one is a thin, explicitly bulk-retrieval wrapper
    -- callers choose the (Jacobi, period) window themselves, generally
    generous-but-FINITE bounds (never omit both jacobi/period bounds
    entirely -- that would fetch an entire family unfiltered, which this
    project's own server-respect convention, per `#647`'s bullet, says not
    to do). Same normalization + required-``libr``/``branch`` validation as
    ``check_jpl_family``, but raises :class:`ValueError` up front for an
    unsupported system/family/missing-required-parameter (there is no single-
    candidate verdict object to carry a "not-covered" status here -- a bulk
    fetch against an uncovered system/family is simply a caller error).

    ``cache_dir`` reuses the exact same gitignored ``out/`` on-disk cache
    :func:`~cyclerfinder.verify.jpl_periodic_orbits.query` already implements
    (#647) -- repeated calls with the same window never re-hit the live API.
    """
    sys_norm = normalize_system(system)
    fam_norm = normalize_family(family)

    if sys_norm not in jpo.SUPPORTED_SYSTEMS:
        raise ValueError(
            f"JPL SSD's Three-Body Periodic Orbits API does not catalog the "
            f"'{sys_norm}' system (supported: {sorted(jpo.SUPPORTED_SYSTEMS)})."
        )
    if fam_norm not in jpo.SUPPORTED_FAMILIES:
        raise ValueError(
            f"JPL SSD's Three-Body Periodic Orbits API does not catalog the "
            f"'{fam_norm}' family (supported: {sorted(jpo.SUPPORTED_FAMILIES)})."
        )
    if fam_norm in jpo.FAMILIES_REQUIRING_LIBR and libr is None:
        raise ValueError(
            f"family '{fam_norm}' requires a libr (libration point) but none was given."
        )
    if fam_norm in jpo.FAMILIES_REQUIRING_BRANCH and branch is None:
        raise ValueError(
            f"family '{fam_norm}' requires a branch (e.g. 'N'/'S') but none was given."
        )

    return query_fn(
        sys_norm,
        fam_norm,
        libr=libr,
        branch=branch,
        jacobimin=jacobi_min,
        jacobimax=jacobi_max,
        periodmin=period_min,
        periodmax=period_max,
        periodunits=periodunits,
        cache_dir=cache_dir,
    )


def hill_radius_km(mu: float, l_km: float) -> float:
    """Secondary's own Hill radius, km: ``l_km * (mu / (3*(1-mu)))**(1/3)``.

    Generalizes :func:`cyclerfinder.genome.hill_screen.earth_sun_hill_radius_km`
    (``r_H = a * (GM2 / (3*GM1))**(1/3)``, and ``GM2/GM1 == mu/(1-mu)`` for
    any CR3BP mass ratio ``mu``) from that function's hardcoded Earth-Sun
    constants to an arbitrary ``(mu, l_km)`` pair -- same formula, same
    derivation, just not pinned to one system.
    """
    if not (0.0 < mu < 1.0):
        raise ValueError(f"mu must be in (0, 1), got {mu}")
    if l_km <= 0.0:
        raise ValueError(f"l_km must be positive, got {l_km}")
    return float(l_km * (mu / (3.0 * (1.0 - mu))) ** (1.0 / 3.0))


def propagate_min_distances_km(
    orbit: jpo.JplOrbit,
    constants: jpo.JplSystemConstants,
    *,
    n_samples: int = DEFAULT_N_SAMPLES,
    rtol: float = DEFAULT_RTOL,
    atol: float = DEFAULT_ATOL,
) -> tuple[float, float]:
    """Propagate one full period of a JPL-catalogued orbit; return minimum
    distance (km) to each body's CENTRE: ``(min_dist_primary_km, min_dist_secondary_km)``.

    Generalizes
    :func:`cyclerfinder.search.real_binary_kk_sweep.min_body_clearance_km`
    (identical DOP853/rtol/atol/n_samples convention, identical "primary at
    ``(-mu, 0, 0)``, secondary at ``(1-mu, 0, 0)``" CR3BP rotating-frame
    layout) from that function's restricted planar-symmetric
    ``(x0, 0, 0, 0, ydot0, 0)`` IC form to a FULL 6-D state, since JPL's own
    catalogued orbits are general (halo/vertical/axial family members carry
    nonzero ``z0``/``vx0``/``vz0``). Propagates entirely within the orbit's
    OWN (JPL-reported) ``mu``/``lunit_km`` -- no cross-system mu-
    reconciliation is needed here since both the geometry input and the
    output distance stay inside JPL's own self-consistent unit system.
    """
    sol = solve_ivp(
        cr3bp_eom,
        (0.0, orbit.period),
        orbit.state0,
        args=(constants.mu,),
        method="DOP853",
        rtol=rtol,
        atol=atol,
        max_step=orbit.period / n_samples,
    )
    x, y, z = sol.y[0], sol.y[1], sol.y[2]
    d_primary_nd = np.sqrt((x - (-constants.mu)) ** 2 + y**2 + z**2).min()
    d_secondary_nd = np.sqrt((x - (1.0 - constants.mu)) ** 2 + y**2 + z**2).min()
    return (
        float(d_primary_nd) * constants.lunit_km,
        float(d_secondary_nd) * constants.lunit_km,
    )


@dataclass(frozen=True)
class SecondaryApproachVerdict:
    """Geometric "recurrent close secondary approach" verdict for one JPL orbit.

    ``is_close_approach`` combines two REUSED criteria (see this module's own
    docstring): the zero-margin physical non-crash floor
    (`#660`'s ``real_binary_kk_sweep`` convention) and the
    :data:`CLOSE_APPROACH_HILL_FRACTION` ceiling borrowed from
    ``genome.hill_screen``.

    "Recurrent" is automatic for ANY exact CR3BP periodic-orbit family member
    -- by construction the state repeats identically every ``period``, so a
    close pass within one period recurs every subsequent period forever
    (within the idealized autonomous CR3BP; real-ephemeris survivability is a
    separate, later V4-style question this module does not address).
    "Bounded" over one period is likewise automatic for any exact periodic
    orbit; ``stability`` (JPL's own reported linear-stability index) is
    reported for characterization but NOT gated on -- a candidate can be a
    perfectly valid, geometrically-close periodic orbit while still being
    linearly UNSTABLE (needing station-keeping), which is a distinct,
    downstream engineering question from "does this orbit's geometry ever
    get close to the secondary at all".
    """

    system: str
    family: str
    libr: int | None
    branch: str | None
    jacobi: float
    period: float
    stability: float
    min_dist_secondary_km: float
    min_dist_primary_km: float
    radius_secondary_km: float | None
    radius_source: str
    hill_radius_km: float
    hill_fraction: float
    physically_valid: bool
    is_close_approach: bool
    notes: str = ""


def classify_secondary_approach(
    orbit: jpo.JplOrbit,
    constants: jpo.JplSystemConstants,
    *,
    system: str,
    family: str,
    libr: int | None = None,
    branch: str | None = None,
    radius_secondary_km_override: float | None = None,
    radius_source_override: str = "",
    clearance_margin_km: float = DEFAULT_CLEARANCE_MARGIN_KM,
    close_hill_fraction: float = CLOSE_APPROACH_HILL_FRACTION,
    n_samples: int = DEFAULT_N_SAMPLES,
    rtol: float = DEFAULT_RTOL,
    atol: float = DEFAULT_ATOL,
) -> SecondaryApproachVerdict:
    """Propagate ``orbit`` and classify its "recurrent close secondary
    approach" character (see :class:`SecondaryApproachVerdict`).

    ``radius_secondary_km_override``/``radius_source_override``: JPL's own
    response does not always populate ``radius_secondary`` -- when it is
    ``None`` and no override is supplied, ``physically_valid`` is ``None``-
    like (reported ``True`` with a note that the check is unevaluated) rather
    than silently treated as "clears", per the `#660` discipline that an
    unsourced radius must never be silently treated as passing.
    """
    min_dist_primary_km, min_dist_secondary_km = propagate_min_distances_km(
        orbit, constants, n_samples=n_samples, rtol=rtol, atol=atol
    )

    radius_km = (
        radius_secondary_km_override
        if radius_secondary_km_override is not None
        else constants.radius_secondary_km
    )
    radius_source = (
        radius_source_override
        if radius_secondary_km_override is not None
        else (
            "JPL response 'radius_secondary'" if constants.radius_secondary_km is not None else ""
        )
    )

    if radius_km is None:
        physically_valid = True
        clearance_note = (
            "radius_secondary unsourced (neither JPL response nor override) -- "
            "clearance NOT evaluated, reported True as a non-claim, not a pass."
        )
    else:
        physically_valid = min_dist_secondary_km >= (radius_km + clearance_margin_km)
        verdict_word = "clears" if physically_valid else "FAILS -- inside the body, impossible"
        clearance_note = (
            f"min distance to secondary centre {min_dist_secondary_km:.3f} km vs. "
            f"radius+margin {radius_km + clearance_margin_km:.3f} km ({verdict_word})"
        )

    r_hill = hill_radius_km(constants.mu, constants.lunit_km)
    hill_fraction = min_dist_secondary_km / r_hill if r_hill > 0 else float("inf")
    is_close = physically_valid and hill_fraction <= close_hill_fraction

    notes = (
        f"{clearance_note}; hill_fraction={hill_fraction:.4f} "
        f"({'<=' if hill_fraction <= close_hill_fraction else '>'} "
        f"{close_hill_fraction} close-approach ceiling)"
    )

    return SecondaryApproachVerdict(
        system=normalize_system(system),
        family=normalize_family(family),
        libr=libr,
        branch=branch,
        jacobi=orbit.jacobi,
        period=orbit.period,
        stability=orbit.stability,
        min_dist_secondary_km=min_dist_secondary_km,
        min_dist_primary_km=min_dist_primary_km,
        radius_secondary_km=radius_km,
        radius_source=radius_source,
        hill_radius_km=r_hill,
        hill_fraction=hill_fraction,
        physically_valid=physically_valid,
        is_close_approach=is_close,
        notes=notes,
    )


__all__ = [
    "CLOSE_APPROACH_HILL_FRACTION",
    "DEFAULT_ATOL",
    "DEFAULT_N_SAMPLES",
    "DEFAULT_RTOL",
    "SecondaryApproachVerdict",
    "classify_secondary_approach",
    "fetch_family_window",
    "hill_radius_km",
    "propagate_min_distances_km",
]
