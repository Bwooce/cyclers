"""Heliocentric inertial planet-state provider.

The :class:`Ephemeris` class supplies ``(r_km, v_km_s)`` heliocentric inertial
state vectors for the planets defined in :mod:`cyclerfinder.core.constants`.
Two backends:

- ``"circular"`` (M1) — every planet on its mean-motion circle in the
  ecliptic plane. At ``t_sec = 0`` every planet sits at ``theta = 0``.
- ``"astropy"`` (M6-slice, 2026-06-01) — real heliocentric position/velocity
  from astropy's bundled JPL DE440 ephemeris, computed at the absolute epoch
  ``2000-01-01T12:00:00 TDB + t_sec``. ``+x`` toward the J2000 vernal equinox,
  ``+z`` along the ecliptic north pole.

The two backends share the ``Ephemeris.state(body, t_sec)`` signature. Callers
that don't care about the epoch convention can swap one for the other freely.
Callers that DO care (M6b phase-matching, real launch-window finders) must use
``"astropy"`` because only it returns true heliocentric positions on real
dates; ``"circular"`` is a per-cycler idealisation with no absolute calendar
meaning.

References
----------
Spec §4, §6, §12.1 (phase-match). Plans:
``docs/phases/m1-core-mechanics/plan.md`` §3.1 (circular backend),
``docs/phases/m4-enumeration-scoring/plan.md`` (M4 consumer notes),
and the launch-windows slice (2026-06-01) that brought the astropy backend
forward from full M6.
"""

from __future__ import annotations

from collections import OrderedDict
from collections.abc import Sequence
from dataclasses import replace
from math import cos, pi, radians, sin
from typing import Final, Protocol

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import AU_KM, PLANETS, SECONDS_PER_DAY, PlanetData

Vec3 = NDArray[np.float64]  # shape (3,), dtype float64

# Epoch convention for the astropy backend: t_sec is TDB seconds since J2000
# (2000-01-01T12:00:00 TDB, JD 2451545.0 TDB) — identical to SPICE ET and to
# nbody/convert.py's harness axis. The backend builds Time("J2000",
# scale="tdb") + TimeDelta(t_sec) directly in the TDB scale; the J2000 instant
# is spelled as the JD constant below. (The pre-2026-06-11 code built the
# epoch from a UTC datetime's POSIX timestamp relabeled scale="tdb" — the
# design §0 "~64.184 s TDB↔TT/UTC trap" nbody/convert.py warns about — which
# put t_sec=0 at J2000 TDB + 64.184 s, ~1944 km of Earth along-track motion.)
_J2000_TDB_JD: Final[float] = 2451545.0

# Mean obliquity of the ecliptic at J2000 (IAU 2006 / IERS 2010). Used by
# the astropy backend to rotate ICRS (BCRS-aligned, equatorial) into the
# J2000 ecliptic frame the rest of cyclerfinder expects (z-axis along the
# ecliptic north pole).
_J2000_OBLIQUITY_RAD: Final[float] = radians(23.4392911)

# Body-code → astropy body name map. Derived from the central PLANETS registry
# so any (sourced) body added there is automatically resolvable here: astropy's
# solar-system body names are the lowercased planet names (``"venus"``,
# ``"earth"``, ``"mars"``, ``"jupiter"`` …). Bodies whose astropy slug differs
# from ``name.lower()`` (e.g. a barycenter or a moon) get an explicit override.
_ASTROPY_NAME_OVERRIDES: Final[dict[str, str]] = {}
# DE440 (the bundled astropy ephemeris) carries the eight major planets and the
# Pluto barycenter, but NOT the main-belt / TNO small bodies added 2026-06-14
# (#260): Ceres/Eris/Makemake/Haumea/Vesta/Pallas have no DE440 slot, so they
# cannot have a real-ephemeris state and are excluded from the astropy map.
# They remain fully usable via the circular / inclined backends (the Tisserand
# screen uses those). Exclude by NAME so a later sourced small body is caught
# too. ``get_body_barycentric_posvel`` would raise for these anyway; excluding
# them here keeps the advertised astropy body set honest.
_NO_DE440_EPHEMERIS: Final[frozenset[str]] = frozenset(
    {"Ceres", "Eris", "Makemake", "Haumea", "Vesta", "Pallas"}
)
_ASTROPY_BODY_NAMES: dict[str, str] = {
    code: _ASTROPY_NAME_OVERRIDES.get(code, data.name.lower())
    for code, data in PLANETS.items()
    if data.name not in _NO_DE440_EPHEMERIS
}


# Sourced J2000 inclination / longitude-of-ascending-node, Standish & Williams,
# "Approximate Positions of the Planets", JPL Solar System Dynamics, Table 1
# (valid 1800-2050 AD) — the same source quoted in the PLANETS comments
# (constants.py). Earth defines the ecliptic, so its inc/lan stay 0.0. These
# live here (not in PLANETS) so the live coplanar table is never mutated.
_INCLINED_ELEMENTS_J2000: Final[dict[str, tuple[float, float]]] = {
    "V": (3.39467605, 76.67984255),
    "M": (1.84969142, 49.55953891),
    "Me": (7.00497902, 48.33076593),
    "J": (1.30439695, 100.47390909),
    "S": (2.48599187, 113.66242448),
    "U": (0.77263783, 74.01692503),
    "N": (1.77004347, 131.78422574),
    # Dwarf planets / planetoids (added 2026-06-14, #260). These bodies are
    # significantly inclined, so the 3D Tisserand screen needs the real i/Omega.
    # Pluto: Standish & Williams Table 2a (the only Standish table with Pluto).
    # The rest: JPL SBDB osculating elements at epoch JD 2461200.5 (i = "in",
    # Omega = "om"), accessed 2026-06-14 (same source as constants.py sma/ecc).
    "Pl": (17.14104260, 110.30167986),
    "Ce": (10.58802780183462, 80.24862682043221),
    "Er": (43.9258279471791, 36.00477044417249),
    "Mk": (29.02785603743067, 79.2948338209406),
    "Ha": (28.20847393040364, 121.7860561329425),
    "Ve": (7.143925545058711, 103.701293265032),
    "Pa": (34.93279321851542, 172.8866193357694),
}


def inclined_planets() -> dict[str, PlanetData]:
    """Return a NEW planet table with the sourced J2000 inc/Ω filled in.

    A copy of :data:`~cyclerfinder.core.constants.PLANETS` where each body's
    ``inc_deg``/``lan_deg`` are replaced with the Standish & Williams Table 1
    values (Venus i=3.39467605°, Ω=76.67984255°; Mars i=1.84969142°,
    Ω=49.55953891°). Bodies without a tabulated entry (Earth, the
    ecliptic-defining body) keep their live ``0.0`` defaults.

    The live ``PLANETS`` dict is **never mutated** — this returns a fresh dict of
    fresh :class:`PlanetData` records (via :func:`dataclasses.replace`), so the
    coplanar default stays byte-identical for every existing caller. Intended to
    be injected into :class:`_InclinedCircularBackend` (see
    :meth:`Ephemeris.inclined_circular`).
    """
    out: dict[str, PlanetData] = {}
    for code, data in PLANETS.items():
        if code in _INCLINED_ELEMENTS_J2000:
            inc_deg, lan_deg = _INCLINED_ELEMENTS_J2000[code]
            out[code] = replace(data, inc_deg=inc_deg, lan_deg=lan_deg)
        else:
            out[code] = replace(data)
    return out


class _Backend(Protocol):
    """Strategy interface — both backends conform to a single state() method."""

    def state(self, body: str, t_sec: float) -> tuple[Vec3, Vec3]:  # pragma: no cover
        ...


def _circular_inplane_state(planet: PlanetData, t_sec: float) -> tuple[Vec3, Vec3]:
    """In-plane circular state in the orbital plane (z == 0, prograde).

    The shared kernel for both the flat and inclined backends. The returned
    state lives in the orbital plane with the reference direction (where the
    body sits at ``t_sec == 0``, i.e. the ascending-node direction for an
    inclined orbit) along ``+x``.
    """
    a_km = planet.sma_au * AU_KM
    # Convert tabulated mean motion (deg/day) to rad/s.
    n_rad_s = planet.mean_motion_deg_day * (pi / 180.0) / SECONDS_PER_DAY
    theta = n_rad_s * t_sec
    cos_t = cos(theta)
    sin_t = sin(theta)
    # Circular orbital speed = a * n.
    speed = a_km * n_rad_s
    r = np.array([a_km * cos_t, a_km * sin_t, 0.0], dtype=np.float64)
    v = np.array([-speed * sin_t, speed * cos_t, 0.0], dtype=np.float64)
    return r, v


class _CircularBackend:
    """Circular, prograde planet motion in the J2000 ecliptic frame.

    Coplanar (z == 0) by default. For a planet whose ``inc_deg != 0.0`` this
    delegates to :class:`_InclinedCircularBackend`, which rotates the in-plane
    circular state into a tilted orbital plane. The ``inc_deg == 0.0`` test is
    an EXACT float comparison so the coplanar path stays byte-identical for any
    body that has not been given a (sourced) non-zero inclination.
    """

    def __init__(self) -> None:
        self._inclined = _InclinedCircularBackend()

    def state(self, body: str, t_sec: float) -> tuple[Vec3, Vec3]:
        # KeyError is the right error if the caller passes an unknown body.
        planet = PLANETS[body]
        if planet.inc_deg != 0.0:
            return self._inclined.state(body, t_sec)
        return _circular_inplane_state(planet, t_sec)


class _CentredCircularBackend:
    """Circular moon motion about a *primary* (planet), km-scaled.

    The planet-centric Tier-1 analogue of :class:`_CircularBackend`: a moon rides
    its mean-motion circle ABOUT THE PRIMARY (Jupiter/Saturn/…), returning
    planet-centred inertial ``(r, v)`` in km / km·s⁻¹. It reads
    :data:`cyclerfinder.core.satellites.SATELLITES` (not ``PLANETS``), uses
    ``sma_km`` directly (NO ``* AU_KM``) and the about-primary
    ``mean_motion_deg_day`` (Kepler III with the primary's GM). Coplanar (z == 0),
    every moon at ``theta == 0`` at ``t_sec == 0`` — the same convention as the
    heliocentric circular backend, one centre down.

    Only moons whose ``primary`` matches the backend's ``center`` are resolvable;
    any other code raises ``KeyError``.
    """

    def __init__(self, center: str) -> None:
        # Lazy import keeps core.ephemeris -> core.satellites edge one-directional
        # and avoids a cycle if satellites ever needs ephemeris helpers.
        from cyclerfinder.core.satellites import SATELLITES

        self._center = center
        self._moons = {code: sat for code, sat in SATELLITES.items() if sat.primary == center}

    def state(self, body: str, t_sec: float) -> tuple[Vec3, Vec3]:
        # KeyError is the right error if the caller passes an unknown moon (or a
        # moon of a different primary than this backend's centre).
        moon = self._moons[body]
        a_km = moon.sma_km
        n_rad_s = moon.mean_motion_deg_day * (pi / 180.0) / SECONDS_PER_DAY
        theta = n_rad_s * t_sec
        cos_t = cos(theta)
        sin_t = sin(theta)
        speed = a_km * n_rad_s
        r = np.array([a_km * cos_t, a_km * sin_t, 0.0], dtype=np.float64)
        v = np.array([-speed * sin_t, speed * cos_t, 0.0], dtype=np.float64)
        return r, v


class _InclinedCircularBackend:
    """Circular planet motion in an orbital plane inclined to the ecliptic.

    The in-plane circular state (z == 0, body along ``+x`` at ``t_sec == 0``)
    is rotated into the J2000 ecliptic frame by ``R_z(+lan) @ R_x(+inc)``
    (the Standish ascending-node convention): first tilt the plane about the
    node line (the ``+x`` reference direction) by the inclination, then rotate
    about the ecliptic ``+z`` by the longitude of the ascending node. With this
    convention the body sits on the ASCENDING node at ``t_sec == 0`` (z == 0
    there, crossing the ecliptic going north), and the orbit normal points
    along ``n_hat = (sin(lan) sin(inc), -cos(lan) sin(inc), cos(inc))``.

    For ``inc_deg == 0.0`` the rotation is the identity, so the result is
    numerically identical to the flat in-plane state.
    """

    def __init__(self, planets: dict[str, PlanetData] | None = None) -> None:
        self._planets: dict[str, PlanetData] = PLANETS if planets is None else planets
        # Per-body orbital-plane -> ecliptic rotation is constant in t_sec (it
        # depends only on inc/Ω), so build it ONCE at construction instead of
        # rebuilding 5 trig + 2 matmuls on every state() call. Matrices are the
        # CORRECT R_z(+lan) @ R_x(+inc) (post #199 sign fix) — cached verbatim.
        self._rotations: dict[str, NDArray[np.float64]] = {
            code: self._rotation(data) for code, data in self._planets.items()
        }

    def state(self, body: str, t_sec: float) -> tuple[Vec3, Vec3]:
        planet = self._planets[body]
        r_plane, v_plane = _circular_inplane_state(planet, t_sec)
        rot = self._rotations[body]
        r = rot @ r_plane
        v = rot @ v_plane
        return (
            np.asarray(r, dtype=np.float64),
            np.asarray(v, dtype=np.float64),
        )

    @staticmethod
    def _rotation(planet: PlanetData) -> NDArray[np.float64]:
        """Orbital-plane -> ecliptic rotation ``R_z(+lan) @ R_x(+inc)``."""
        inc = radians(planet.inc_deg)
        lan = radians(planet.lan_deg)
        ci, si = cos(inc), sin(inc)
        cl, sl = cos(lan), sin(lan)
        # R_x(+inc): tilt the orbital plane about the node line (+x) so the
        # orbit normal acquires a -y, +z component before the node rotation
        # (Standish ascending-node convention: the body crosses the ecliptic
        # going NORTH at the node). The previous R_x(-inc) form mirrored every
        # orbital plane about the ecliptic (normals 2*inc away from DE440).
        rx = np.array(
            [
                [1.0, 0.0, 0.0],
                [0.0, ci, -si],
                [0.0, si, ci],
            ],
            dtype=np.float64,
        )
        # R_z(+lan): rotate the node line to its ecliptic longitude.
        rz = np.array(
            [
                [cl, -sl, 0.0],
                [sl, cl, 0.0],
                [0.0, 0.0, 1.0],
            ],
            dtype=np.float64,
        )
        return rz @ rx


class _AstropyBackend:
    """Real heliocentric states from astropy's bundled JPL DE440 ephemeris.

    Returns heliocentric ``(r_km, v_km_s)`` in the **J2000 ecliptic
    frame** (z-axis along the ecliptic north pole) at the astropy
    ``Time`` corresponding to ``J2000 + t_sec`` (TDB).

    astropy's :func:`get_body_barycentric_posvel` returns BCRS
    coordinates, which are aligned with ICRS — the equatorial mean
    pole and equinox of J2000. The cyclerfinder ``circular`` backend
    and the rest of the package (frames, propagator, dynamic frame)
    work in the **ecliptic** frame, so we rotate the ICRS output by
    the J2000 obliquity (-23.4393 deg about +x) to produce ecliptic-
    aligned coordinates. Without this rotation Earth's heliocentric
    z-component oscillates by ~ ±58 million km over a year (the
    23.4393 deg equator/ecliptic tilt), which silently breaks any
    rotating-frame transform that assumes z is the ecliptic normal —
    exactly the M6a spec §10 risk.
    """

    def __init__(self) -> None:
        # Import lazily so packages that omit the runtime dep can still import
        # this module's circular path.
        from astropy.coordinates import solar_system_ephemeris
        from astropy.time import Time

        # DE440 is bundled with astropy as of 6.x+; opt in once.
        solar_system_ephemeris.set("de440")
        # t_sec=0 reference: exactly J2000 TDB (JD 2451545.0 in the TDB scale).
        # Built ONCE as a jd-format TDB Time so the t_sec offset is pure TDB
        # second arithmetic — never a POSIX/unix timestamp relabeled "tdb"
        # (that relabeling is the ~64.184 s trap; see module-level comment).
        self._j2000_tdb = Time(_J2000_TDB_JD, format="jd", scale="tdb")
        # Pre-compute the ICRS->ecliptic rotation about +x by -obliquity.
        eps = _J2000_OBLIQUITY_RAD
        self._r_icrs_to_ecl: NDArray[np.float64] = np.array(
            [
                [1.0, 0.0, 0.0],
                [0.0, cos(eps), sin(eps)],
                [0.0, -sin(eps), cos(eps)],
            ],
            dtype=np.float64,
        )

    def state(self, body: str, t_sec: float) -> tuple[Vec3, Vec3]:
        from astropy.coordinates import get_body_barycentric_posvel
        from astropy.time import TimeDelta

        if body not in _ASTROPY_BODY_NAMES:
            raise KeyError(
                f"unknown body code {body!r}; astropy backend handles {tuple(_ASTROPY_BODY_NAMES)}"
            )
        astropy_name = _ASTROPY_BODY_NAMES[body]
        # J2000(TDB) + t_sec, entirely in the TDB scale: a scale-less
        # TimeDelta added to a TDB Time is applied as TDB seconds (uniform,
        # no leap-second steps), so t_sec=0 lands EXACTLY on JD 2451545.0 TDB.
        t = self._j2000_tdb + TimeDelta(t_sec, format="sec")
        # Heliocentric posvel: subtract Sun's barycentric posvel from the
        # body's. astropy returns CartesianRepresentation in km and km/s when
        # we cast through .xyz.to(km).value.
        body_pos, body_vel = get_body_barycentric_posvel(astropy_name, t)
        sun_pos, sun_vel = get_body_barycentric_posvel("sun", t)
        r_helio_icrs = (body_pos - sun_pos).xyz.to("km").value
        v_helio_icrs = (body_vel - sun_vel).xyz.to("km/s").value
        # Rotate ICRS (equatorial) to J2000 ecliptic so the rest of
        # cyclerfinder gets z-along-ecliptic-pole positions matching
        # the circular backend's convention.
        r_helio = self._r_icrs_to_ecl @ np.asarray(r_helio_icrs, dtype=np.float64)
        v_helio = self._r_icrs_to_ecl @ np.asarray(v_helio_icrs, dtype=np.float64)
        return (
            np.asarray(r_helio, dtype=np.float64),
            np.asarray(v_helio, dtype=np.float64),
        )

    def states(self, bodies: Sequence[str], epochs: Sequence[float]) -> list[tuple[Vec3, Vec3]]:
        """Vectorised batch of :meth:`state` over parallel ``(body, epoch)`` lists.

        Collapses the per-call ``Time`` construction and units/posvel framework
        overhead into array-``Time`` calls: one ``Time`` over all DISTINCT
        epochs, the Sun posvel computed once per distinct epoch (shared across
        every body at that epoch), and each body's posvel grouped so astropy's
        Chebyshev evaluator runs once per (body, epoch-array). Per-element output
        is byte-identical to looping :meth:`state` — same DE440 states, same
        ICRS→ecliptic rotation, same subtraction order.
        """
        from astropy.coordinates import get_body_barycentric_posvel
        from astropy.time import TimeDelta

        unknown = [b for b in bodies if b not in _ASTROPY_BODY_NAMES]
        if unknown:
            raise KeyError(
                f"unknown body code(s) {tuple(unknown)!r}; astropy backend handles "
                f"{tuple(_ASTROPY_BODY_NAMES)}"
            )
        n = len(bodies)
        # Distinct epochs (preserve first-seen order) → one array Time + one Sun
        # posvel per distinct epoch. Match the scalar path's epoch arithmetic:
        # J2000(TDB) + t_sec as a TDB-scale TimeDelta (see state()).
        epoch_index: dict[float, int] = {}
        tdb_offsets_sec: list[float] = []
        for t_sec in epochs:
            if t_sec not in epoch_index:
                epoch_index[t_sec] = len(tdb_offsets_sec)
                tdb_offsets_sec.append(t_sec)
        t_arr = self._j2000_tdb + TimeDelta(
            np.asarray(tdb_offsets_sec, dtype=np.float64), format="sec"
        )
        sun_pos, sun_vel = get_body_barycentric_posvel("sun", t_arr)

        # Group element indices by astropy body so each body's Chebyshev eval is
        # a single array call over its (distinct) epochs.
        by_body: dict[str, list[int]] = {}
        for i, b in enumerate(bodies):
            by_body.setdefault(b, []).append(i)

        results: list[tuple[Vec3, Vec3] | None] = [None] * n
        rot = self._r_icrs_to_ecl
        for body, idxs in by_body.items():
            astropy_name = _ASTROPY_BODY_NAMES[body]
            uniq = list(dict.fromkeys(epochs[i] for i in idxs))
            cols = [epoch_index[t] for t in uniq]
            sub_t = t_arr[cols]
            body_pos, body_vel = get_body_barycentric_posvel(astropy_name, sub_t)
            sun_pos_sub = sun_pos[cols]
            sun_vel_sub = sun_vel[cols]
            r_icrs = (body_pos - sun_pos_sub).xyz.to("km").value
            v_icrs = (body_vel - sun_vel_sub).xyz.to("km/s").value
            # xyz is shape (3, k); rotate all columns at once.
            r_ecl = rot @ np.asarray(r_icrs, dtype=np.float64)
            v_ecl = rot @ np.asarray(v_icrs, dtype=np.float64)
            col_for_epoch = {t: c for c, t in enumerate(uniq)}
            for i in idxs:
                c = col_for_epoch[epochs[i]]
                results[i] = (
                    np.array(r_ecl[:, c], dtype=np.float64),
                    np.array(v_ecl[:, c], dtype=np.float64),
                )
        return [r for r in results if r is not None]


# Default per-instance state() cache capacity. Profile-informed: one S1L1
# astropy solve touches ~70-80 distinct (body,epoch) pairs (Target 2 in
# docs/notes/2026-06-06-performance-profile.md), and the maintenance window
# scan + DE pass revisit overlapping launch epochs across generations. 4096
# entries comfortably covers a whole solve / window scan while bounding memory:
# each entry is a (str, float) key plus two float64 (3,) arrays ≈ 200 bytes, so
# the full cache is well under ~1 MB per Ephemeris instance.
_DEFAULT_STATE_CACHE_SIZE: Final[int] = 4096


class Ephemeris:
    """Planet-state provider; selects the backend at construction time.

    Parameters
    ----------
    model:
        Backend name.

        - ``"circular"`` (default) — analytic states for the planets in
          :data:`cyclerfinder.core.constants.PLANETS`; perfect circles in the
          ecliptic, every planet at ``theta = 0`` at ``t_sec = 0``.
        - ``"astropy"`` — real heliocentric states from astropy's bundled
          JPL DE440 ephemeris. ``t_sec = 0`` corresponds to the J2000 epoch
          (2000-01-01T12:00:00 TDB).
    cache:
        Memoise exact ``(body, t_sec)`` results in a per-instance bounded LRU
        (default ``True``). The astropy backend does no caching of its own and
        ~36-42% of the ``state()`` calls in a single solve are exact duplicates
        (the finite-difference jacobian holds ``t0`` fixed; the maintenance scan
        revisits launch epochs across DE generations). The cache returns FRESH
        copies on every access, so behaviour is byte-identical to the uncached
        path (proven by the replay fixture in ``tests/core``). The cache is
        per-instance, so a ramped backend (``continuation.ramped_ephemeris``),
        the circular backend, and a DE440 instance never share entries.
    cache_size:
        Max number of distinct ``(body, t_sec)`` entries retained (LRU
        eviction). Defaults to :data:`_DEFAULT_STATE_CACHE_SIZE`.
    """

    def __init__(
        self,
        model: str = "circular",
        *,
        center: str | None = None,
        cache: bool = True,
        cache_size: int = _DEFAULT_STATE_CACHE_SIZE,
    ) -> None:
        self._cache_enabled: bool = cache
        self._cache_size: int = cache_size
        # Keyed on (body, t_sec); values are (r, v) float64 (3,) tuples kept
        # read-only internally and copied out on every access.
        self._state_cache: OrderedDict[tuple[str, float], tuple[Vec3, Vec3]] = OrderedDict()
        if center is not None:
            # Planet-centric Tier-1 moon ephemeris: a moon on its mean-motion
            # circle about ``center`` (a primary). Only the circular model is
            # supported for centred moon states in Tier-1; the heliocentric
            # backends below are entered only when ``center is None``, so they
            # stay byte-identical.
            if model != "circular":
                raise ValueError(
                    f"center={center!r} requires model='circular' "
                    f"(Tier-1 centred moon ephemeris); got model={model!r}",
                )
            self._backend_impl: _Backend = _CentredCircularBackend(center)
            self._center: str | None = center
            self._model: str = model
            return
        self._center = None
        if model == "circular":
            self._backend_impl = _CircularBackend()
        elif model == "inclined-circular":
            # Opt-in inclined-circular backend: real J2000 inc/Ω (Standish &
            # Williams), mean sma, circular (eccentricity ignored — separable
            # follow-on). Built from a COPY of PLANETS so the live coplanar
            # table is never mutated.
            self._backend_impl = _InclinedCircularBackend(inclined_planets())
        elif model == "astropy":
            self._backend_impl = _AstropyBackend()
        else:
            raise ValueError(
                f"unknown ephemeris model {model!r}; expected 'circular', "
                "'inclined-circular', or 'astropy'",
            )
        self._model = model

    @property
    def _backend(self) -> _Backend:
        return self._backend_impl

    @_backend.setter
    def _backend(self, backend: _Backend) -> None:
        # Reassigning the backend (e.g. continuation.ramped_ephemeris injects a
        # _RampedElementsBackend post-construction) invalidates any cached
        # states computed by the previous backend, preventing cross-model
        # contamination on the SAME instance.
        self._backend_impl = backend
        self._state_cache.clear()

    @classmethod
    def inclined_circular(cls) -> Ephemeris:
        """Construct an :class:`Ephemeris` on the opt-in inclined-circular backend.

        Convenience for ``Ephemeris(model="inclined-circular")``: planets ride
        their mean-motion circles in orbital planes tilted by the sourced J2000
        inclination/node (:func:`inclined_planets`). The default
        ``Ephemeris(model="circular")`` path is unaffected and byte-identical.
        """
        return cls(model="inclined-circular")

    @property
    def model(self) -> str:
        """Name of the backend in use."""
        return self._model

    @property
    def cache_enabled(self) -> bool:
        """Whether per-instance ``state()`` memoisation is active."""
        return self._cache_enabled

    def state(self, body: str, t_sec: float) -> tuple[Vec3, Vec3]:
        """Return ``(r_km, v_km_s)`` heliocentric inertial state for ``body``.

        Parameters
        ----------
        body:
            One-letter body code from
            :data:`cyclerfinder.core.constants.PLANETS` (``"V"``, ``"E"``,
            ``"M"``). Unknown codes raise :class:`KeyError`.
        t_sec:
            Time since the M1 reference epoch, seconds. May be negative.

        Returns
        -------
        (r_km, v_km_s):
            Two ``numpy.ndarray`` of shape ``(3,)`` and ``dtype=float64``, in
            the heliocentric ecliptic J2000 frame described in the module
            docstring.
        """
        if not self._cache_enabled:
            return self._backend_impl.state(body, t_sec)
        key = (body, t_sec)
        cache = self._state_cache
        hit = cache.get(key)
        if hit is None:
            r, v = self._backend_impl.state(body, t_sec)
            # Store a read-only canonical copy; never the array we hand out.
            r_store = np.array(r, dtype=np.float64)
            v_store = np.array(v, dtype=np.float64)
            r_store.flags.writeable = False
            v_store.flags.writeable = False
            cache[key] = (r_store, v_store)
            if len(cache) > self._cache_size:
                cache.popitem(last=False)  # LRU evict oldest
            return r_store.copy(), v_store.copy()
        cache.move_to_end(key)  # LRU: mark most-recently-used
        r_store, v_store = hit
        # Fresh writeable copies every access => byte-identical to uncached path.
        return r_store.copy(), v_store.copy()

    def states(self, bodies: Sequence[str], epochs: Sequence[float]) -> list[tuple[Vec3, Vec3]]:
        """Vectorised batch of :meth:`state` over parallel ``(body, epoch)`` lists.

        ``bodies[i]`` is evaluated at ``epochs[i]``; the two sequences must be
        the same length. Returns a list of ``(r_km, v_km_s)`` tuples in the same
        order, each element identical to ``state(bodies[i], epochs[i])``.

        For the astropy backend this collapses the per-call ``Time`` construction
        and ``get_body_barycentric_posvel`` framework overhead (~2x the actual
        Chebyshev math, per the profile note) into a handful of array-``Time``
        calls. Cached ``(body, epoch)`` pairs are served from the per-instance
        LRU and only the misses are batched. For the analytic backends it is a
        thin loop over :meth:`state` (no vectorisation benefit, same results).

        Callers can adopt this where they already have the full epoch grid up
        front (e.g. a phase-match window scan). Per-element behaviour is
        byte-identical to repeated :meth:`state` calls.
        """
        if len(bodies) != len(epochs):
            raise ValueError(f"bodies/epochs length mismatch: {len(bodies)} != {len(epochs)}")
        results: list[tuple[Vec3, Vec3] | None] = [None] * len(bodies)
        misses: list[int] = []
        for i, (body, t_sec) in enumerate(zip(bodies, epochs, strict=True)):
            if self._cache_enabled:
                hit = self._state_cache.get((body, t_sec))
                if hit is not None:
                    self._state_cache.move_to_end((body, t_sec))
                    r_store, v_store = hit
                    results[i] = (r_store.copy(), v_store.copy())
                    continue
            misses.append(i)

        if misses:
            batch_states = getattr(self._backend_impl, "states", None)
            if batch_states is not None:
                computed = batch_states([bodies[i] for i in misses], [epochs[i] for i in misses])
            else:
                computed = [self._backend_impl.state(bodies[i], epochs[i]) for i in misses]
            for i, (r, v) in zip(misses, computed, strict=True):
                if self._cache_enabled:
                    r_store = np.array(r, dtype=np.float64)
                    v_store = np.array(v, dtype=np.float64)
                    r_store.flags.writeable = False
                    v_store.flags.writeable = False
                    self._state_cache[(bodies[i], epochs[i])] = (r_store, v_store)
                    if len(self._state_cache) > self._cache_size:
                        self._state_cache.popitem(last=False)
                    results[i] = (r_store.copy(), v_store.copy())
                else:
                    results[i] = (
                        np.array(r, dtype=np.float64),
                        np.array(v, dtype=np.float64),
                    )
        return [r for r in results if r is not None]

    def clear_cache(self) -> None:
        """Drop all memoised ``state()`` results for this instance."""
        self._state_cache.clear()
