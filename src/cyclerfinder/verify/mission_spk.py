"""Fetch + cache NAIF spacecraft SPK kernels and extract planetocentric flyby V∞.

This module is the #390 *wholesale unblock* for the #345 classic-mission
``mga_tour`` backlog. Every Voyager / Mariner-10 mission-overview AND dedicated
navigation paper acquired publishes encounter dates + closest-approach geometry
+ maneuver ΔV but NOT the per-encounter hyperbolic-excess velocity (V∞) the §14
V0 ``mga_tour`` standard requires (see
``docs/notes/2026-06-19-345-voyager-mariner-mission-digests.md``). The V∞ is
instead DERIVED here from each mission's archived NAIF SPK kernel at the
already-sourced flyby epochs.

This is *validation infrastructure only* — nothing here participates in the
production construct/score/verify pipeline (seeds-not-tracks is intact). It
extends, rather than reinvents, the on-demand NAIF-fetch pattern already used by
:mod:`cyclerfinder.verify.spice_kernels` for the leapseconds kernel, and reuses
the spiceypy ``furnsh``/``spkezr``/``str2et`` machinery already exercised by
:mod:`cyclerfinder.data.validation.v4_uranus_strict`.

Discipline
----------
* **Sourced-only.** Body gravitational parameters come from the project
  constants module (:data:`cyclerfinder.core.constants.PLANETS`) — never
  hard-coded. Every extracted V∞ carries SPK-kernel filename + NAIF body ID +
  epoch provenance.
* **Derived, not published.** The V∞ values produced here are computed by our
  own code from the kernels. They are NOT auto-admitted to the catalogue; #390
  recommends rows only, and admission is a parent-reviewed follow-on (the
  "it closed!" danger signal of feedback_orbit_closure_discipline applies to
  derived numbers feeding admission).
* **Binary kernels are never committed** — fetched on demand into the astropy
  cache dir alongside the cached DE440 / LSK.

Kernel sources (documented, not committed)
------------------------------------------
Spacecraft SPKs are fetched from the NAIF public archive. The planet-relative
reconstructed SPKs (e.g. ``vgr2_jup230.bsp``) are tiny (a few hundred kB) and
cover exactly the flyby window with the spacecraft state given relative to the
planet's barycenter, which is what the planetocentric-V∞ extraction needs.
"""

from __future__ import annotations

import math
import os
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import spiceypy as spice

from cyclerfinder.core.constants import PLANETS
from cyclerfinder.verify.spice_kernels import (
    astropy_de440_bsp_path,
    ensure_leapseconds_kernel,
)

# NAIF spacecraft-SPK base URLs (public archive). Confirmed reachable + exact
# filenames verified against the NAIF directory index on 2026-06-19.
NAIF_VOYAGER_SPK_BASE = "https://naif.jpl.nasa.gov/pub/naif/VOYAGER/kernels/spk/"
NAIF_M10_SPK_BASE = "https://naif.jpl.nasa.gov/pub/naif/M10/kernels/spk/"

# NAIF body IDs for the spacecraft (CSPICE built-in name->ID; given explicitly so
# the extractor never depends on a name-resolution kernel being loaded).
VOYAGER_1_NAIF_ID = -31
VOYAGER_2_NAIF_ID = -32
MARINER_10_NAIF_ID = -76


# Map flyby-body short names to (NAIF body name spiceypy understands, PLANETS key
# for the sourced GM). spkezr targets the planet BARYCENTER name for the giant
# planets; for the planetocentric V∞ that is the correct gravitating center used
# by the reconstructed planet-relative SPKs.
@dataclass(frozen=True)
class FlybyBody:
    """A flyby target: its SPICE center name + the project GM key."""

    spice_center: str
    planets_key: str


FLYBY_BODIES: dict[str, FlybyBody] = {
    "Jupiter": FlybyBody(spice_center="JUPITER BARYCENTER", planets_key="J"),
    "Saturn": FlybyBody(spice_center="SATURN BARYCENTER", planets_key="S"),
    "Uranus": FlybyBody(spice_center="URANUS BARYCENTER", planets_key="U"),
    "Neptune": FlybyBody(spice_center="NEPTUNE BARYCENTER", planets_key="N"),
    "Venus": FlybyBody(spice_center="VENUS BARYCENTER", planets_key="V"),
    "Mercury": FlybyBody(spice_center="MERCURY", planets_key="Me"),
}


def _spice_cache_dir(cache_dir: str | os.PathLike[str] | None = None) -> Path:
    """Return the on-disk cache dir for mission SPKs (astropy cache subdir)."""
    if cache_dir is None:
        from astropy.config.paths import get_cache_dir

        cache_dir = Path(get_cache_dir()) / "cyclerfinder_spice"
    cache_path = Path(cache_dir)
    cache_path.mkdir(parents=True, exist_ok=True)
    return cache_path


def ensure_mission_spk(
    filename: str,
    *,
    base_url: str,
    cache_dir: str | os.PathLike[str] | None = None,
) -> str:
    """Return a path to a mission SPK ``filename``, fetching it once if necessary.

    Mirrors :func:`cyclerfinder.verify.spice_kernels.ensure_leapseconds_kernel`:
    the kernel is fetched into the ``cyclerfinder_spice`` astropy-cache subdir
    (never the repo). Network is only touched on the first call; subsequent calls
    reuse the cached file.

    Parameters
    ----------
    filename:
        Exact SPK filename in the NAIF directory (e.g. ``"vgr2_jup230.bsp"``).
    base_url:
        NAIF directory URL the file lives in (one of the module constants).
    """
    cache_path = _spice_cache_dir(cache_dir)
    spk_path = cache_path / filename
    if not spk_path.exists():
        url = base_url.rstrip("/") + "/" + filename
        urllib.request.urlretrieve(url, spk_path)
    return str(spk_path)


# --------------------------------------------------------------------------- #
# Planetocentric V∞ extraction
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class VinfSample:
    """One epoch's planetocentric state sample around a flyby."""

    et_seconds: float
    offset_minutes: float
    """Minutes relative to the nominal closest-approach epoch."""
    r_rel_km: float
    """Spacecraft distance from the flyby body, km."""
    v_rel_kms: float
    """Spacecraft speed relative to the flyby body, km/s."""
    vinf_visviva_kms: float
    """sqrt(max(0, v_rel^2 - 2*mu/r_rel)) — the vis-viva hyperbolic-excess
    speed at THIS epoch. Converges to the true V∞ far from the body."""


@dataclass(frozen=True)
class FlybyVinf:
    """Extracted planetocentric V∞ for one (mission, flyby body, epoch)."""

    mission: str
    sc_naif_id: int
    flyby_body: str
    epoch_utc: str
    spk_filename: str
    mu_km3_s2: float
    body_radius_km: float
    vinf_kms: float
    """Best-estimate asymptotic V∞ — the vis-viva value at the OUTERMOST
    sampled epoch (where the planet's pull is weakest and v_rel -> V∞)."""
    vinf_kms_visviva_window_mean: float
    """Mean vis-viva V∞ across the outer-window samples (>= 1 SOI radius)."""
    vinf_kms_visviva_window_std: float
    """Std of the same — the convergence/stability evidence (<1% target)."""
    closest_approach_radius_km: float
    """Minimum spacecraft-to-body distance over the sampled window, km."""
    closest_approach_altitude_km: float
    """closest_approach_radius_km - body_radius_km (cross-check vs published)."""
    closest_approach_radius_body_radii: float
    """closest_approach_radius_km / body_radius_km (Kohlhase-Penzo Table IV
    publishes closest approach in planet radii — direct self-consistency)."""
    samples: tuple[VinfSample, ...]


def vinf_at_flyby(
    mission_spk: str,
    sc_naif_id: int,
    flyby_body: str,
    epoch_utc: str,
    *,
    mission: str = "",
    window_minutes: float = 4320.0,
    n_samples: int = 25,
    extra_kernels: tuple[str, ...] = (),
) -> FlybyVinf:
    """Extract the planetocentric hyperbolic-excess velocity (V∞) at a flyby.

    Pipeline:

    1. ``furnsh`` DE440 (for the planet barycenter ephemeris) + LSK (UTC->ET)
       + the mission SPK (spacecraft state).
    2. Convert ``epoch_utc`` to ET (``str2et``).
    3. Over a symmetric ``window_minutes`` window around the epoch, sample the
       spacecraft state relative to ``flyby_body`` (``spkezr``, J2000,
       body-centered) at ``n_samples`` epochs.
    4. At each epoch compute the vis-viva hyperbolic-excess speed
       ``sqrt(max(0, v_rel^2 - 2*mu/r_rel))``. Near periapsis this is depressed
       by the deep potential; far outside the SOI it converges to the true
       asymptotic V∞.
    5. Report the V∞ as the outermost-epoch vis-viva value, plus the
       outer-window mean/std as convergence evidence, and the closest-approach
       radius (cross-check vs the published flyby geometry).

    Body GM + radius are sourced from :data:`cyclerfinder.core.constants.PLANETS`
    (never hard-coded). ``window_minutes`` defaults to 3 days each side, large
    enough that the outer samples sit well beyond the giant-planet SOI where
    v_rel -> V∞.

    Kernel isolation: ``kclear`` is called in a finally block so the SPICE pool
    is never polluted across calls.
    """
    if flyby_body not in FLYBY_BODIES:
        raise ValueError(f"unknown flyby body {flyby_body!r}; known: {sorted(FLYBY_BODIES)}")
    if n_samples < 5:
        raise ValueError(f"n_samples must be >= 5 for convergence evidence; got {n_samples}")
    if window_minutes <= 0.0:
        raise ValueError(f"window_minutes must be > 0; got {window_minutes}")

    body = FLYBY_BODIES[flyby_body]
    planet = PLANETS[body.planets_key]
    mu = float(planet.mu_km3_s2)
    body_radius_km = float(planet.radius_eq_km)

    de440 = astropy_de440_bsp_path()
    lsk = ensure_leapseconds_kernel()

    spice.kclear()
    try:
        spice.furnsh(de440)
        spice.furnsh(lsk)
        spice.furnsh(mission_spk)
        for k in extra_kernels:
            spice.furnsh(k)

        et0 = float(spice.str2et(epoch_utc))
        half = window_minutes * 60.0
        offsets_s = np.linspace(-half, half, n_samples)

        samples: list[VinfSample] = []
        for ds in offsets_s:
            et = et0 + float(ds)
            state, _lt = spice.spkezr(str(sc_naif_id), et, "J2000", "NONE", body.spice_center)
            arr = np.asarray(state, dtype=np.float64)
            r_rel = float(np.linalg.norm(arr[:3]))
            v_rel = float(np.linalg.norm(arr[3:]))
            c3 = v_rel * v_rel - 2.0 * mu / r_rel
            vinf = math.sqrt(c3) if c3 > 0.0 else 0.0
            samples.append(
                VinfSample(
                    et_seconds=et,
                    offset_minutes=float(ds) / 60.0,
                    r_rel_km=r_rel,
                    v_rel_kms=v_rel,
                    vinf_visviva_kms=vinf,
                )
            )
    finally:
        spice.kclear()

    radii = np.array([s.r_rel_km for s in samples])
    ca_idx = int(np.argmin(radii))
    ca_radius = float(radii[ca_idx])

    # Outer window = samples beyond 5x the closest-approach radius (well outside
    # the deep potential), used for the converged V∞ + stability evidence. If the
    # SPK window is short, fall back to the outer half of the samples.
    outer = [s for s in samples if s.r_rel_km >= 5.0 * ca_radius]
    if len(outer) < 3:
        ordered = sorted(samples, key=lambda s: s.r_rel_km)
        outer = ordered[-max(3, n_samples // 2) :]
    outer_vinf = np.array([s.vinf_visviva_kms for s in outer])

    # Best-estimate V∞: the value at the single OUTERMOST sample.
    outermost = max(samples, key=lambda s: s.r_rel_km)

    return FlybyVinf(
        mission=mission,
        sc_naif_id=int(sc_naif_id),
        flyby_body=flyby_body,
        epoch_utc=epoch_utc,
        spk_filename=os.path.basename(mission_spk),
        mu_km3_s2=mu,
        body_radius_km=body_radius_km,
        vinf_kms=float(outermost.vinf_visviva_kms),
        vinf_kms_visviva_window_mean=float(np.mean(outer_vinf)),
        vinf_kms_visviva_window_std=float(np.std(outer_vinf)),
        closest_approach_radius_km=ca_radius,
        closest_approach_altitude_km=ca_radius - body_radius_km,
        closest_approach_radius_body_radii=ca_radius / body_radius_km,
        samples=tuple(samples),
    )


__all__ = [
    "FLYBY_BODIES",
    "MARINER_10_NAIF_ID",
    "NAIF_M10_SPK_BASE",
    "NAIF_VOYAGER_SPK_BASE",
    "VOYAGER_1_NAIF_ID",
    "VOYAGER_2_NAIF_ID",
    "FlybyBody",
    "FlybyVinf",
    "VinfSample",
    "ensure_mission_spk",
    "vinf_at_flyby",
]
