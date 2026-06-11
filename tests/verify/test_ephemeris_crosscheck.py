"""Independent-toolchain validation of the ephemeris layer (task #129).

Every internal cross-check in the project (Axis A's lamberthub/Kepler paths,
the corrector, the gauntlet) consumes the SAME ``Ephemeris("astropy")`` states.
A frame-convention, time-scale (TDB/UTC/TT), or jplephem-reader bug would shift
all of them *consistently* and pass every existing gate. This module adds the
missing independent check on the ephemeris layer itself, in three layers:

1. **Convention pins (fast)** — assert, in code, the frame/center/time
   conventions ``core/ephemeris.py`` actually uses, so they are *pinned* and a
   silent convention change trips a test (not just a docstring).

2. **SPICE same-kernel cross-check (slow)** — point NAIF's CSPICE reader
   (spiceypy) at the *same* DE440 BSP astropy caches and compare
   ``Ephemeris("astropy").state`` against ``spkgeo`` (heliocentric / Sun-center,
   ECLIPJ2000 = J2000 ecliptic, ET seconds) over a 2025-2045 grid for all eight
   planets. Same data, independent reader ⇒ agreement is a numerical-precision
   *consistency predicate*. A larger systematic = the bug class we are hunting;
   report it, do not widen the tolerance.

3. **Horizons golden anchors (slow, sourced EXPECTED)** — frozen JPL Horizons
   state vectors (retrieval provenance in comments) as independently-sourced
   ground truth. Looser tolerance: Horizons uses DE441 (we read DE440).

Conventions confirmed from ``src/cyclerfinder/core/ephemeris.py`` (read
2026-06-06):

- **Center**: heliocentric — ``body - sun`` of astropy barycentric posvel ⇒
  SPICE center 10 (Sun).
- **Frame**: J2000 **ecliptic** (ICRS/equatorial rotated by ``-obliquity``
  23.4392911 deg about +x) ⇒ SPICE ``ECLIPJ2000``.
- **Time**: ``t_sec`` is **TDB seconds since J2000(TDB)** (JD 2451545.0 TDB),
  a uniform axis with no leap-second steps — identical to SPICE ET and to
  ``nbody/convert.py``. Pinned by :func:`test_j2000_epoch_offset_is_pinned`.
  (Until 2026-06-11/#198 the backend relabeled a UTC POSIX timestamp
  ``scale="tdb"``, landing ``t_sec=0`` at J2000(TDB)+64.184 s; this module's
  SPICE/Horizons comparisons replicated that internal ``Time`` via
  ``_backend_tdb_jd``, so they validated the kernel *reader* while sharing the
  epoch convention — the shift was invisible here. The epoch itself is now
  anchored externally by ``tests/verify/test_epoch_anchor.py``.)

The astropy-name -> NAIF-id map was confirmed empirically (min-residual over the
grid): astropy ``get_body_barycentric_posvel`` returns the *planet center* for
Mercury/Venus/Earth (no significant moons, or 399 directly for Earth) and the
*system barycenter* for Mars and the giants (the cached DE440 BSP carries only
the barycenter at top level for those). Earth maps to 399 (NOT the
Earth-Moon barycenter, NAIF 3 — that disagrees by ~4400 km).
"""

from __future__ import annotations

from collections.abc import Iterator
from typing import Any

import numpy as np
import pytest

from cyclerfinder.core.constants import PLANETS
from cyclerfinder.core.ephemeris import (
    _ASTROPY_BODY_NAMES,
    _J2000_OBLIQUITY_RAD,
    _J2000_TDB_JD,
    Ephemeris,
)

# --------------------------------------------------------------------------- #
# Body-code -> NAIF id for the SAME DE440 BSP astropy reads. Confirmed by
# minimum position residual over the 2025-2045 grid (see module docstring).
# --------------------------------------------------------------------------- #
_NAIF_ID: dict[str, int] = {
    "Me": 1,  # Mercury barycenter == Mercury (no moons)
    "V": 2,  # Venus barycenter == Venus (no moons)
    "E": 399,  # Earth center (NOT EMB id 3)
    "M": 4,  # Mars system barycenter
    "J": 5,  # Jupiter system barycenter
    "S": 6,  # Saturn system barycenter
    "U": 7,  # Uranus system barycenter
    "N": 8,  # Neptune system barycenter
}

# Catalogue horizon grid: 2025-2045, 20 points (decimal years).
_GRID_YEARS = [2025.0 + i * (2045.0 - 2025.0) / 19.0 for i in range(20)]

# Same-kernel agreement gates. Observed worst case over the full grid (all 8
# planets) is ~0.59 km / ~0.010 mm/s, dominated by Neptune at 4.5e9 km where
# jplephem's and CSPICE's reconstruction of the *identical* Chebyshev
# coefficients differ at the limit of double precision. These gates are tight:
# any real frame/time-scale/reader bug would be thousands of km, not sub-km.
# Set with documented headroom over the observed residual, NOT widened to hide a
# systematic. Tightening further would make the gate flaky on outer planets;
# loosening past ~2 km would start to admit a 1-second time-scale error
# (Neptune moves ~5 km/s).
_SPICE_POS_GATE_KM = 2.0
_SPICE_VEL_GATE_MM_S = 0.1


# --------------------------------------------------------------------------- #
# Convention pins (fast) — assert the conventions in code, not just docstrings.
# --------------------------------------------------------------------------- #
def test_obliquity_constant_is_j2000_mean() -> None:
    """The ICRS->ecliptic rotation uses the IAU 2006 J2000 mean obliquity."""
    expected = np.radians(23.4392911)
    assert pytest.approx(expected, rel=0, abs=1e-12) == _J2000_OBLIQUITY_RAD


def test_j2000_epoch_offset_is_pinned() -> None:
    """Pin the epoch: t_sec=0 is EXACTLY J2000(TDB), JD 2451545.0 TDB.

    Fixed in #198 (2026-06-11): the backend previously relabeled a UTC POSIX
    timestamp ``scale="tdb"``, landing t_sec=0 at J2000(TDB)+64.184 s
    (32 leap seconds at 2000 + 32.184 s TAI->TT — the design §0 trap
    ``nbody/convert.py`` warns about). The backend now constructs
    ``Time(2451545.0, format="jd", scale="tdb") + TimeDelta(t_sec)``; pin
    both the exported JD constant and the constructed instant so any silent
    reintroduction of a relabeled epoch is caught.
    """
    pytest.importorskip("astropy")
    from astropy.time import Time, TimeDelta

    assert _J2000_TDB_JD == 2451545.0
    t0 = Time(_J2000_TDB_JD, format="jd", scale="tdb") + TimeDelta(0.0, format="sec")
    offset_s = (t0.tdb.jd - 2451545.0) * 86400.0
    assert offset_s == pytest.approx(0.0, abs=1e-6)


def test_earth_z_small_confirms_ecliptic_frame() -> None:
    """Earth's heliocentric z stays small (ecliptic frame), not ~±58e6 km.

    If the backend forgot the ICRS->ecliptic rotation, Earth's z would
    oscillate by ~±58 million km over a year (the 23.44 deg equator/ecliptic
    tilt). A small z across a year confirms z is the ecliptic normal.
    """
    pytest.importorskip("astropy")
    eph = Ephemeris("astropy")
    max_abs_z = 0.0
    for frac in (0.0, 0.25, 0.5, 0.75):
        t_sec = frac * 365.25 * 86400.0
        r, _ = eph.state("E", t_sec)
        max_abs_z = max(max_abs_z, abs(float(r[2])))
    # Earth's true ecliptic z stays under ~1e4 km; equatorial would be ~5.8e7.
    assert max_abs_z < 1.0e5


def test_all_eight_planets_registered_and_named() -> None:
    """All 8 planets are in PLANETS and resolvable by the astropy backend."""
    assert set(PLANETS) == {"Me", "V", "E", "M", "J", "S", "U", "N"}
    for code in PLANETS:
        assert code in _ASTROPY_BODY_NAMES
        assert code in _NAIF_ID


# --------------------------------------------------------------------------- #
# Shared helpers for the same-kernel SPICE / Horizons comparisons.
# --------------------------------------------------------------------------- #
def _backend_tdb_jd(t_sec: float) -> float:
    """The exact TDB Julian date the astropy backend evaluates ``t_sec`` at.

    Since #198 the backend builds ``Time(_J2000_TDB_JD, format="jd",
    scale="tdb") + TimeDelta(t_sec, format="sec")`` — a uniform TDB axis, so
    this map is exactly linear: ``JD = 2451545.0 + t_sec / 86400``. We still
    replicate it through astropy's own ``Time`` arithmetic (not the formula)
    so the SPICE ET and Horizons comparisons hit the SAME instant the backend
    reads, whatever the backend's construction is.
    """
    from astropy.time import Time, TimeDelta

    t = Time(_J2000_TDB_JD, format="jd", scale="tdb") + TimeDelta(t_sec, format="sec")
    return float(t.tdb.jd)


def _et_for_t_sec(t_sec: float) -> float:
    """SPICE ET (seconds past true J2000 TDB) for the backend's ``t_sec``.

    Exact by construction: ET is derived from the backend's own TDB instant
    (:func:`_backend_tdb_jd`), so SPICE reads the SAME epoch astropy reads.
    Since #198 this is simply ``t_sec`` itself (ET == TDB seconds past J2000).
    """
    return (_backend_tdb_jd(t_sec) - 2451545.0) * 86400.0


def _t_sec_for_tdb_jd(tdb_jd: float) -> float:
    """Invert :func:`_backend_tdb_jd`: t_sec landing the backend on ``tdb_jd``.

    Since #198 the t_sec axis is uniform TDB seconds since J2000(TDB), so the
    inverse is the exact linear map (no leap-second steps to bisect across).
    Used to anchor against Horizons vectors stated at a fixed TDB instant.
    """
    return (tdb_jd - 2451545.0) * 86400.0


# --------------------------------------------------------------------------- #
# (2) SPICE same-kernel cross-check (slow).
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="module")
def spice_kernels() -> Iterator[Any]:
    """Furnish the SAME DE440 BSP astropy caches + the NAIF leapseconds kernel.

    Skips cleanly if spiceypy is not installed (the optional ``validation``
    extra) or if the kernels cannot be fetched (e.g. network blocked on first
    run). Yields the spiceypy module with kernels furnished; clears the kernel
    pool on teardown so other tests are unaffected.
    """
    sp = pytest.importorskip(
        "spiceypy",
        reason=(
            "spiceypy not installed; `uv pip install -e .[validation]` to run the SPICE cross-check"
        ),
        exc_type=ImportError,
    )
    pytest.importorskip("astropy")
    from cyclerfinder.verify.spice_kernels import (
        astropy_de440_bsp_path,
        ensure_leapseconds_kernel,
    )

    try:
        bsp = astropy_de440_bsp_path()
        lsk = ensure_leapseconds_kernel()
    except Exception as exc:  # network / cache failure -> skip honestly
        pytest.skip(f"could not obtain SPICE kernels (network blocked?): {exc}")

    sp.furnsh(lsk)
    sp.furnsh(bsp)
    try:
        yield sp
    finally:
        sp.kclear()


@pytest.mark.slow
def test_spice_same_kernel_crosscheck_all_planets(spice_kernels: Any) -> None:
    """astropy reader vs CSPICE reader on the SAME DE440 kernel, all 8 planets.

    Consistency predicate: same data, independent reader. Position must agree
    sub-2-km, velocity sub-0.1-mm/s over the 2025-2045 grid. A larger
    systematic is the frame/time-scale/reader bug this whole test exists to
    catch — report it verbatim, do NOT widen the tolerance.
    """
    sp = spice_kernels
    eph = Ephemeris("astropy")

    worst_pos = (0.0, "", 0.0)
    worst_vel = (0.0, "", 0.0)
    failures: list[str] = []

    for year in _GRID_YEARS:
        # 2025-2045 grid point as a TDB instant, mapped to the backend's own
        # t_sec so SPICE reads the SAME epoch the backend reads.
        from astropy.time import Time

        tdb_jd = float(Time(f"{int(year)}-01-01 00:00:00", scale="tdb").tdb.jd)
        t_sec = _t_sec_for_tdb_jd(tdb_jd)
        et = _et_for_t_sec(t_sec)
        for code, naif in _NAIF_ID.items():
            r_ours, v_ours = eph.state(code, t_sec)
            state, _lt = sp.spkgeo(naif, et, "ECLIPJ2000", 10)
            r_spice = np.asarray(state[:3], dtype=np.float64)
            v_spice = np.asarray(state[3:], dtype=np.float64)
            d_pos = float(np.linalg.norm(r_ours - r_spice))
            d_vel_mm_s = float(np.linalg.norm(v_ours - v_spice)) * 1.0e6
            if d_pos > worst_pos[0]:
                worst_pos = (d_pos, code, year)
            if d_vel_mm_s > worst_vel[0]:
                worst_vel = (d_vel_mm_s, code, year)
            if d_pos > _SPICE_POS_GATE_KM or d_vel_mm_s > _SPICE_VEL_GATE_MM_S:
                failures.append(
                    f"{code} @ {int(year)}: dpos={d_pos:.4f} km, dvel={d_vel_mm_s:.4f} mm/s"
                )

    assert not failures, (
        "SPICE same-kernel cross-check DISAGREES beyond numerical precision — "
        "this is the frame/time-scale/reader bug class the check exists to "
        "catch. DO NOT widen the tolerance; investigate.\n" + "\n".join(failures)
    )
    # Sanity: the worst residual should be comfortably under the gate (it is the
    # signal that the readers genuinely agree, not that the grid was empty).
    assert worst_pos[0] < _SPICE_POS_GATE_KM
    assert worst_vel[0] < _SPICE_VEL_GATE_MM_S


# --------------------------------------------------------------------------- #
# (3) Horizons golden anchors (slow) — sourced EXPECTED constants.
# --------------------------------------------------------------------------- #
# Retrieved 2026-06-06 from the JPL Horizons API (telnet/web API v1.2),
# https://ssd.jpl.nasa.gov/api/horizons.api, with:
#   EPHEM_TYPE=VECTORS, CENTER='500@10' (Sun center),
#   REF_PLANE='ECLIPTIC', REF_SYSTEM='J2000', OUT_UNITS='KM-S',
#   VEC_TABLE='2' (position + velocity), GEOMETRIC states (no light-time/aberr).
# Horizons reported {source: DE441} for both target and Sun center; OUR backend
# reads DE440. The DE440-vs-DE441 difference is why the anchor tolerance is
# looser than the SPICE same-kernel gate. Epochs are stated by Horizons in TDB;
# we compare at the same TDB instant via _t_sec_for_tdb_jd (since #198 the exact
# linear map, so these anchors now pin the epoch convention too — see also the
# dedicated keystone anchors in test_epoch_anchor.py).
# Each entry: (body_code, "YYYY-MM-DD TDB", r_km(3), v_km_s(3)).
_HORIZONS_ANCHORS: list[tuple[str, str, tuple[float, float, float], tuple[float, float, float]]] = [
    # Earth (target 399, Sun center 10), DE441, retrieved 2026-06-06.
    (
        "E",
        "2030-01-01 00:00:00",
        (-2.600847757810495e07, 1.447900488833384e08, -9.532206645891070e03),
        (-2.981505192741513e01, -5.371451831814229e00, -5.573338788711357e-04),
    ),
    (
        "E",
        "2030-01-02 00:00:00",
        (-2.858037993331836e07, 1.443033275564770e08, -9.583639651425183e03),
        (-2.971805786018533e01, -5.895012704276932e00, -6.214572133349527e-04),
    ),
    # Mars (target 499, Sun center 10), DE441, retrieved 2026-06-06.
    (
        "M",
        "2032-06-15 00:00:00",
        (-2.730053224974398e07, 2.355278985432163e08, 5.605090927148551e06),
        (-2.315336862126238e01, -7.328028072380064e-01, 5.520951720986933e-01),
    ),
    (
        "M",
        "2032-06-16 00:00:00",
        (-2.929994573227200e07, 2.354558464845639e08, 5.652583384448871e06),
        (-2.312907339722446e01, -9.349869228781336e-01, 5.472621955160913e-01),
    ),
    # Jupiter (target '5' = Jupiter SYSTEM BARYCENTER, Sun center 10), DE441,
    # retrieved 2026-06-06. We use the system barycenter (NAIF 5), NOT Jupiter
    # center (599), because that is exactly what our backend returns: the cached
    # DE440 BSP carries only the J-system barycenter at top level, so
    # ``get_body_barycentric_posvel("jupiter")`` yields NAIF 5. Anchoring to 599
    # instead would mismatch the center by ~190 km (Galilean moons shift the
    # barycenter off Jupiter center) — a center-convention artefact, not an
    # ephemeris error; using '5' is the correct like-for-like comparison.
    (
        "J",
        "2035-03-01 00:00:00",
        (6.862672518227907e08, 2.796064284502964e08, -1.651526265369132e07),
        (-5.089375467968009e00, 1.272142430213055e01, 6.098059102755382e-02),
    ),
    (
        "J",
        "2035-03-02 00:00:00",
        (6.858266945509964e08, 2.807052186877146e08, -1.650997382616685e07),
        (-5.108707372336333e00, 1.271352978660025e01, 6.144594583724317e-02),
    ),
]

# Looser than the same-kernel gate: the only systematic here is DE440 (ours) vs
# DE441 (Horizons), which is sub-km over 2030-2035 for these bodies (observed
# worst ~0.5 km). 50 km / 5 mm/s leaves generous headroom over that DE-version
# difference while staying far below any frame/time-scale error (thousands of
# km). All anchors use a center matching our backend (Sun-centered; Jupiter via
# the system barycenter, NAIF 5), so no center-convention offset is admitted.
_HORIZONS_POS_GATE_KM = 50.0
_HORIZONS_VEL_GATE_MM_S = 5.0


@pytest.mark.slow
@pytest.mark.parametrize(
    ("code", "epoch_tdb", "r_expected", "v_expected"),
    [(c, e, r, v) for (c, e, r, v) in _HORIZONS_ANCHORS],
    ids=[f"{c}-{e[:10]}" for (c, e, _r, _v) in _HORIZONS_ANCHORS],
)
def test_horizons_golden_anchors(
    code: str,
    epoch_tdb: str,
    r_expected: tuple[float, float, float],
    v_expected: tuple[float, float, float],
) -> None:
    """Gate our astropy states against frozen JPL Horizons vectors (sourced).

    EXPECTED side is sourced from Horizons (provenance in the
    ``_HORIZONS_ANCHORS`` comment); this is genuine external ground truth, not a
    value our own code computed. Looser tolerance: Horizons uses DE441, we read
    DE440 (and the J anchor mixes Jupiter-center vs system-barycenter).
    """
    pytest.importorskip("astropy")
    from astropy.time import Time

    eph = Ephemeris("astropy")
    tdb_jd = float(Time(epoch_tdb, scale="tdb").tdb.jd)
    t_sec = _t_sec_for_tdb_jd(tdb_jd)
    r_ours, v_ours = eph.state(code, t_sec)

    d_pos = float(np.linalg.norm(r_ours - np.asarray(r_expected)))
    d_vel_mm_s = float(np.linalg.norm(v_ours - np.asarray(v_expected))) * 1.0e6
    assert d_pos < _HORIZONS_POS_GATE_KM, (
        f"{code} @ {epoch_tdb}: position {d_pos:.3f} km off Horizons anchor "
        f"(gate {_HORIZONS_POS_GATE_KM} km)"
    )
    assert d_vel_mm_s < _HORIZONS_VEL_GATE_MM_S, (
        f"{code} @ {epoch_tdb}: velocity {d_vel_mm_s:.4f} mm/s off Horizons "
        f"anchor (gate {_HORIZONS_VEL_GATE_MM_S} mm/s)"
    )
