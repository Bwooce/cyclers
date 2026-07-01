"""Locate / fetch the SPICE kernels for the independent ephemeris cross-check.

This module is *validation infrastructure only*. It is imported by
``tests/verify/test_ephemeris_crosscheck.py`` to point spiceypy at the EXACT
same DE440 binary kernel that astropy has already cached, plus the NAIF
leapseconds kernel SPICE needs to convert UTC/calendar time to ET. It does
**not** participate in the production trajectory pipeline (seeds-not-tracks is
intact: nothing here is consumed by construct/score/verify of cyclers).

Why same-kernel
---------------
The whole point of the SPICE cross-check is to run an *independent reader*
(NAIF's CSPICE via spiceypy) over the *identical data* (astropy's cached DE440
BSP). If the two readers disagree by more than numerical precision, that is the
exact frame-convention / time-scale / jplephem-reader bug class the cross-check
exists to catch. We therefore deliberately read astropy's own cached BSP rather
than downloading a fresh copy.

Kernel sources (documented, not committed to the repo)
------------------------------------------------------
- DE440 BSP: whatever file ``astropy.coordinates.solar_system_ephemeris`` has
  cached for ``"de440"`` (astropy downloads it from
  ``https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp``
  on first use and stores it under the astropy download cache).
- Leapseconds (LSK): ``naif0012.tls`` from
  ``https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls``
  fetched on demand into the astropy cache dir (binary kernels are never
  committed).
"""

from __future__ import annotations

import os
import urllib.request
from pathlib import Path

# NAIF generic leapseconds kernel. naif0012.tls is the current LSK (covers all
# leap seconds through the present); ET<->UTC conversions in 2025-2045 only need
# the historical leap-second table, so this file is stable.
NAIF_LSK_URL = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls"
NAIF_LSK_FILENAME = "naif0012.tls"


def astropy_de440_bsp_path() -> str:
    """Return the on-disk path of the DE440 BSP astropy has cached.

    Triggers astropy to select DE440 and returns the underlying file the
    jplephem SPK reader has open. This is the SAME binary kernel the production
    ``Ephemeris("astropy")`` backend reads, so handing it to spiceypy gives a
    like-for-like, same-data comparison.

    Raises ``RuntimeError`` if the kernel file path cannot be recovered (e.g. a
    future astropy that no longer exposes the jplephem DAF file handle).
    """
    from astropy.coordinates import solar_system_ephemeris

    solar_system_ephemeris.set("de440")
    kernel = solar_system_ephemeris.kernel
    try:
        path = kernel.daf.file.name  # jplephem SPK -> DAF -> open file handle
    except AttributeError as exc:  # pragma: no cover - defensive
        raise RuntimeError(
            "could not recover astropy's cached DE440 BSP path from the "
            f"jplephem kernel object ({type(kernel).__name__})"
        ) from exc
    if not path or not os.path.exists(path):  # pragma: no cover - defensive
        raise RuntimeError(f"astropy DE440 BSP path {path!r} does not exist on disk")
    return str(path)


NAIF_JUP365_URL = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/jup365.bsp"
NAIF_JUP365_LOCAL = Path("~/dev/references/kernels/jup365.bsp")


def ensure_jup365_kernel() -> str:
    """Return the local path of the jup365.bsp Galilean-satellite kernel.

    Resolves ``~/dev/references/kernels/jup365.bsp`` (expanduser). This kernel
    covers Io, Europa, Ganymede, and Callisto relative to Jupiter barycentre and
    requires NO leapseconds kernel when passing ET floats directly to
    ``spiceypy.spkezr``.

    Does **not** auto-download (jup365.bsp is ~50 MB). If the file is absent,
    a ``RuntimeError`` is raised with the NAIF download URL so the user knows
    exactly where to fetch it.

    Used by the Galilean cycler real-ephemeris pipeline (#480 Task 1).
    """
    path = NAIF_JUP365_LOCAL.expanduser()
    if not path.exists():
        raise RuntimeError(
            f"jup365.bsp kernel not found at {path}. Download it from: {NAIF_JUP365_URL}"
        )
    return str(path)


NAIF_PLU060_URL = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/plu060.bsp"
NAIF_PLU060_LOCAL = Path("~/dev/references/kernels/plu060.bsp")


def ensure_pluto_kernel() -> str:
    """Return the local path of the plu060.bsp Pluto-system kernel.

    Resolves ``~/dev/references/kernels/plu060.bsp`` (expanduser). This kernel
    covers Charon (901), Nix (902), Hydra (903), Kerberos (904), Styx (905),
    and Pluto (999) relative to the solar-system barycentre, 1800-2199. Same
    NAIF data source (PLU060) already cited in ``core/satellites.py`` for the
    Charon/Nix/Hydra GM values used by the catalogue's Pluto-Charon rows, so
    there is no model-mismatch risk versus the mass data already in the repo.

    Mirrors :func:`ensure_jup365_kernel` exactly: does **not** auto-download
    (129 MB). If the file is absent, a ``RuntimeError`` is raised with the
    NAIF download URL. Used by the Pluto-Charon real-ephemeris differential
    corrector (#511).
    """
    path = NAIF_PLU060_LOCAL.expanduser()
    if not path.exists():
        raise RuntimeError(
            f"plu060.bsp kernel not found at {path}. Download it from: {NAIF_PLU060_URL}"
        )
    return str(path)


def ensure_leapseconds_kernel(cache_dir: str | os.PathLike[str] | None = None) -> str:
    """Return a path to ``naif0012.tls``, fetching it once if necessary.

    The LSK is fetched into ``cache_dir`` (default: a ``cyclerfinder_spice``
    subdirectory of the astropy cache dir, so it sits alongside the cached BSP
    and is never written into the repo). Returns the local path. Network is only
    touched on the first call; subsequent calls reuse the cached file.
    """
    if cache_dir is None:
        from astropy.config.paths import get_cache_dir

        cache_dir = Path(get_cache_dir()) / "cyclerfinder_spice"
    cache_path = Path(cache_dir)
    cache_path.mkdir(parents=True, exist_ok=True)
    lsk_path = cache_path / NAIF_LSK_FILENAME
    if not lsk_path.exists():
        # NAIF_LSK_URL is a fixed https NAIF endpoint (module constant).
        urllib.request.urlretrieve(NAIF_LSK_URL, lsk_path)
    return str(lsk_path)
