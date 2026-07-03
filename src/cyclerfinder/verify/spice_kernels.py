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
  ``https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls``,
  vendored directly in this package (see ``VENDORED_LSK_PATH`` below) --
  unlike the large binary kernels above, the LSK is a ~5 KB plain-text table
  that only changes when IERS announces a new leap second (with ~6 months'
  notice), so committing it removes a genuine CI flake (naif.jpl.nasa.gov was
  observed intermittently unreachable from GitHub Actions runners at ~10% of
  recent CI runs, 2026-07-03) without the size/versioning concerns that keep
  the binary kernels out of the repo.
"""

from __future__ import annotations

import os
import time
import urllib.request
from pathlib import Path

# NAIF generic leapseconds kernel. naif0012.tls is the current LSK (covers all
# leap seconds through the present); ET<->UTC conversions in 2025-2045 only need
# the historical leap-second table, so this file is stable.
NAIF_LSK_URL = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls"
NAIF_LSK_FILENAME = "naif0012.tls"

# Vendored copy -- see the module docstring for why the LSK (unlike the large
# binary kernels) is committed directly instead of fetched over the network.
VENDORED_LSK_PATH = Path(__file__).parent / "kernels" / NAIF_LSK_FILENAME


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


def _download_with_timeout(url: str, dest: Path, *, timeout: float = 20.0) -> None:
    """Download ``url`` to ``dest`` with an EXPLICIT per-attempt socket timeout.

    ``urllib.request.urlretrieve`` -- used here previously -- takes no
    ``timeout`` argument at all, so a single call falls back to the
    process-wide socket default (frequently unset, i.e. block until the OS's
    own TCP timeout, commonly 60-130s on Linux). That defeats a retry loop's
    entire purpose: 3 "retries" against an unresponsive host can each hang
    for over a minute before the loop's own 10s backoff even runs, turning a
    fast-fail into many minutes of dead waiting (observed directly: CI run
    28636314125, 2026-07-02, took 48 minutes before finally failing this
    way). ``urlopen(..., timeout=...)`` DOES support a real per-call timeout;
    stream the response to disk manually instead of using the convenience
    wrapper that lacks the option.
    """
    with urllib.request.urlopen(url, timeout=timeout) as response, dest.open("wb") as fh:
        fh.write(response.read())


def ensure_leapseconds_kernel(cache_dir: str | os.PathLike[str] | None = None) -> str:
    """Return a path to ``naif0012.tls``, preferring the vendored copy.

    Returns :data:`VENDORED_LSK_PATH` directly if present (the normal case --
    see the module docstring for why the LSK, unlike the large binary
    kernels, is committed). Falls back to fetching into ``cache_dir``
    (default: a ``cyclerfinder_spice`` subdirectory of the astropy cache dir)
    only if the vendored file is somehow missing, so a future NAIF leap-second
    update can still be picked up by refreshing the vendored file OR simply
    deleting it to force a live re-fetch.

    Fetches with 3 retries (10s backoff, 20s per-attempt timeout -- see
    :func:`_download_with_timeout`) on a transient network failure -- the
    same class of intermittent NAIF-server timeout the CI workflow's "Warm
    DE440 kernel (retries)" step and the real-eph SPK kernel fetch
    (``curl --retry 3 --retry-delay 10``) already guard against. This
    function previously had no retry at all and a single NAIF timeout during
    a live pytest run failed CI (2026-07-02, run 28588092512) with no
    connection to any code change -- fixed here at the source rather than
    only papering over it with more CI-side caching. A LATER incident
    (2026-07-02, run 28636314125) showed the retry alone was not enough: an
    unbounded per-attempt timeout let 3 retries consume the better part of
    an hour before failing anyway -- the explicit per-attempt timeout above
    is the fix for that second bug. A THIRD incident (2026-07-03, run
    28638767219) showed the fetch itself is genuinely flaky (~10% of recent
    CI runs) even with a bounded timeout -- fixed by vendoring the file so
    the network path is no longer on the critical path at all.
    """
    if VENDORED_LSK_PATH.exists():
        return str(VENDORED_LSK_PATH)
    if cache_dir is None:
        from astropy.config.paths import get_cache_dir

        cache_dir = Path(get_cache_dir()) / "cyclerfinder_spice"
    cache_path = Path(cache_dir)
    cache_path.mkdir(parents=True, exist_ok=True)
    lsk_path = cache_path / NAIF_LSK_FILENAME
    if not lsk_path.exists():
        # NAIF_LSK_URL is a fixed https NAIF endpoint (module constant).
        last_error: Exception | None = None
        for attempt in range(3):
            try:
                _download_with_timeout(NAIF_LSK_URL, lsk_path)
                last_error = None
                break
            except OSError as exc:
                last_error = exc
                if attempt < 2:
                    time.sleep(10)
        if last_error is not None:
            raise last_error
    return str(lsk_path)
