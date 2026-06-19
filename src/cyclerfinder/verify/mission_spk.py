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
NAIF_PIONEER10_SPK_BASE = "https://naif.jpl.nasa.gov/pub/naif/PIONEER10/kernels/spk/"
NAIF_PIONEER11_SPK_BASE = "https://naif.jpl.nasa.gov/pub/naif/PIONEER11/kernels/spk/"
NAIF_JUNO_SPK_BASE = "https://naif.jpl.nasa.gov/pub/naif/JUNO/kernels/spk/"
NAIF_CASSINI_SPK_BASE = "https://naif.jpl.nasa.gov/pub/naif/CASSINI/kernels/spk/"

# NAIF body IDs for the spacecraft (CSPICE built-in name->ID; given explicitly so
# the extractor never depends on a name-resolution kernel being loaded).
VOYAGER_1_NAIF_ID = -31
VOYAGER_2_NAIF_ID = -32
MARINER_10_NAIF_ID = -76
PIONEER_10_NAIF_ID = -23  # confirmed from p10-a.bsp coverage (#399)
PIONEER_11_NAIF_ID = -24  # confirmed from p11-a.bsp coverage (#399)
JUNO_NAIF_ID = -61  # confirmed from the Juno cruise-merge SPK coverage (#399)
CASSINI_NAIF_ID = -82  # confirmed from the Cassini reconstructed-SPK coverage (#399)


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
    # Earth flybys (Cassini, Juno, Galileo cruise) use the EARTH body center
    # (399), NOT the Earth-Moon barycenter: the Moon's mass is non-negligible,
    # and PLANETS["E"].mu is Earth's GM (398600.44), which must match the center.
    "Earth": FlybyBody(spice_center="EARTH", planets_key="E"),
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
    """Minimum spacecraft-to-body distance, km. By default this is the
    automatically-REFINED periapsis (golden-section search between the bracketing
    coarse samples), not merely the best coarse sample — the coarse grid can miss
    the true periapsis by hours on a fast inner flyby (the #398 Voyager-2 Neptune
    bug, which reported 6.2 R_N for a true 1.18 R_N)."""
    closest_approach_altitude_km: float
    """closest_approach_radius_km - body_radius_km (cross-check vs published)."""
    closest_approach_radius_body_radii: float
    """closest_approach_radius_km / body_radius_km (Kohlhase-Penzo Table IV
    publishes closest approach in planet radii — direct self-consistency)."""
    closest_approach_et_seconds: float
    """ET (TDB seconds) of the refined periapsis."""
    closest_approach_offset_minutes: float
    """Refined periapsis epoch relative to the nominal ``epoch_utc``, minutes."""
    closest_approach_refined: bool
    """True if the periapsis was refined past the coarse grid (default path)."""
    samples: tuple[VinfSample, ...]


def _r_rel_at(sc_naif_id: int, et: float, center: str) -> float:
    """Spacecraft-to-center distance (km) at ET — a single ``spkezr`` position read.

    Assumes the required kernels are already ``furnsh``-ed (this is only called
    from inside :func:`vinf_at_flyby`'s furnsh/kclear scope).
    """
    state, _lt = spice.spkezr(str(sc_naif_id), et, "J2000", "NONE", center)
    arr = np.asarray(state, dtype=np.float64)
    return float(np.linalg.norm(arr[:3]))


def _refine_periapsis(
    sc_naif_id: int,
    center: str,
    et_lo: float,
    et_hi: float,
    *,
    tol_seconds: float = 1.0,
    max_iter: int = 64,
) -> tuple[float, float]:
    """Golden-section minimisation of |r_rel|(t) on ``[et_lo, et_hi]``.

    Returns ``(et_periapsis, r_min_km)``. The interval is the bracket around the
    coarse minimum (the min-distance sample and its two neighbours), within which
    |r_rel|(t) is unimodal for a hyperbolic flyby, so golden-section converges
    deterministically. Each iteration costs one ``spkezr`` position read; the
    1-second tolerance closes in <~ log_phi(window/tol) ~ a few dozen reads.

    A noisy SPK could in principle break unimodality; the caller guards the
    result against the coarse minimum and falls back to a dense fine resample if
    the refinement does not improve on it.
    """
    inv_phi = (math.sqrt(5.0) - 1.0) / 2.0  # 1/phi ~ 0.618
    a, b = et_lo, et_hi
    c = b - inv_phi * (b - a)
    d = a + inv_phi * (b - a)
    fc = _r_rel_at(sc_naif_id, c, center)
    fd = _r_rel_at(sc_naif_id, d, center)
    for _ in range(max_iter):
        if (b - a) <= tol_seconds:
            break
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - inv_phi * (b - a)
            fc = _r_rel_at(sc_naif_id, c, center)
        else:
            a, c, fc = c, d, fd
            d = a + inv_phi * (b - a)
            fd = _r_rel_at(sc_naif_id, d, center)
    et_min = c if fc < fd else d
    return et_min, min(fc, fd)


def _refine_periapsis_dense(
    sc_naif_id: int,
    center: str,
    et_lo: float,
    et_hi: float,
    *,
    n: int = 64,
) -> tuple[float, float]:
    """Robust fallback: dense uniform resample of |r_rel|(t) on the bracket.

    Used only if golden-section fails to beat the coarse minimum (e.g. a noisy /
    multi-modal SPK segment), so the CA is still resolved well past the coarse
    grid. ``n`` position reads.
    """
    ets = np.linspace(et_lo, et_hi, n)
    rs = np.array([_r_rel_at(sc_naif_id, float(et), center) for et in ets])
    i = int(np.argmin(rs))
    return float(ets[i]), float(rs[i])


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
    refine_periapsis: bool = True,
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
    5. **Refine the periapsis** (default): the coarse grid (~6 h spacing on the
       default window) can miss the true closest approach by hours — for a fast
       flyby this hugely over-reports the CA radius (the #398 Voyager-2 Neptune
       bug: 6.2 R_N reported for a true 1.18 R_N). After locating the coarse
       minimum, bracket it with its two neighbouring samples and minimise
       |r_rel|(t) by golden-section search (≈1 s tolerance, a couple of dozen
       extra ``spkezr`` reads), with a dense-resample fallback. The reported
       closest-approach radius/epoch are the refined values. Set
       ``refine_periapsis=False`` to keep the legacy coarse-grid behaviour.
    6. Report the V∞ as the outermost-epoch vis-viva value, plus the
       outer-window mean/std as convergence evidence, and the refined
       closest-approach radius (cross-check vs the published flyby geometry).

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

        radii = np.array([s.r_rel_km for s in samples])
        ca_idx = int(np.argmin(radii))
        ca_radius = float(radii[ca_idx])
        ca_et = float(samples[ca_idx].et_seconds)
        ca_refined = False

        if refine_periapsis and n_samples >= 3:
            # Bracket the coarse minimum with its two neighbours (clamped to the
            # sampled window). |r_rel|(t) is unimodal across a hyperbolic flyby,
            # so this bracket contains the true periapsis.
            lo_idx = max(0, ca_idx - 1)
            hi_idx = min(n_samples - 1, ca_idx + 1)
            et_lo = float(samples[lo_idx].et_seconds)
            et_hi = float(samples[hi_idx].et_seconds)
            if et_hi > et_lo:
                r_et, r_min = _refine_periapsis(sc_naif_id, body.spice_center, et_lo, et_hi)
                if r_min > ca_radius:  # golden-section did not beat the grid
                    r_et, r_min = _refine_periapsis_dense(
                        sc_naif_id, body.spice_center, et_lo, et_hi
                    )
                if r_min < ca_radius:
                    ca_radius = r_min
                    ca_et = r_et
                    ca_refined = True
    finally:
        spice.kclear()

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
        closest_approach_et_seconds=ca_et,
        closest_approach_offset_minutes=(ca_et - et0) / 60.0,
        closest_approach_refined=ca_refined,
        samples=tuple(samples),
    )


__all__ = [
    "CASSINI_NAIF_ID",
    "FLYBY_BODIES",
    "JUNO_NAIF_ID",
    "MARINER_10_NAIF_ID",
    "NAIF_CASSINI_SPK_BASE",
    "NAIF_JUNO_SPK_BASE",
    "NAIF_M10_SPK_BASE",
    "NAIF_PIONEER10_SPK_BASE",
    "NAIF_PIONEER11_SPK_BASE",
    "NAIF_VOYAGER_SPK_BASE",
    "PIONEER_10_NAIF_ID",
    "PIONEER_11_NAIF_ID",
    "VOYAGER_1_NAIF_ID",
    "VOYAGER_2_NAIF_ID",
    "FlybyBody",
    "FlybyVinf",
    "VinfSample",
    "ensure_mission_spk",
    "vinf_at_flyby",
]
