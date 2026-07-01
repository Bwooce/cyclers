"""Check tracked NAIF satellite SPK kernels against the live NAIF archive (#514).

Why this exists
----------------
#506 (2026-07-01) assumed the Pluto-Charon satellite kernel was
``plu058.bsp``. #510's survey (same day) found that filename no longer
exists in NAIF's public archive at all -- the current kernel is
``plu060.bsp``, at 129 MB vs. the stale estimate's ~3.8 MB. That was a
one-off manual catch. This script makes the check mechanical and repeatable
so a superseded kernel filename/version is never silently assumed again.

What "tracked" means
---------------------
The :data:`KERNEL_MANIFEST` below lists every kernel the project currently
relies on (see ``docs/notes/2026-07-01-510-naif-spk-kernel-survey.md``):
Jupiter (jup365), Uranus (ura111), Pluto (plu060). Adding a new tracked
system is a dict entry, not new per-body logic -- the comparison, parsing,
and reporting are all generic over :class:`KernelSpec`.

Two independent axes are checked and reported separately (see
:func:`check_kernel`):

* **freshness** -- does NAIF's live directory listing contain a kernel in
  the same numbered series with a HIGHER version number than the one this
  manifest currently pins? This is a pure manifest-vs-remote-listing
  comparison; it does NOT require the kernel to be present on the machine
  running the check, so it works unmodified in CI (which never holds the
  multi-GB binaries). ``STALE`` here is the only thing that fails the exit
  code / drives a CI issue.
* **local presence** -- is the pinned kernel actually sitting at its
  documented local path on THIS machine? Purely informational (mirrors
  ``ensure_jup365_kernel()`` in ``src/cyclerfinder/verify/spice_kernels.py``
  -- local-path-only, no auto-download by default). Absence alone is
  reported as ``MISSING`` but never fails the exit code.

Usage
-----
    uv run python scripts/check_kernel_freshness.py
    uv run python scripts/check_kernel_freshness.py --system pluto
    uv run python scripts/check_kernel_freshness.py --fetch   # download STALE kernels

``--fetch`` is opt-in only (never runs unattended, mirrors the
auto-download-refusal convention in ``spice_kernels.py``): it downloads the
newer kernel found upstream into the SAME local directory the stale one is
pinned to, under its OWN (new) filename -- it never overwrites the old file
or edits this script's manifest for you.
"""

from __future__ import annotations

import argparse
import re
import sys
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

NAIF_SATELLITES_DIR = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/"

_USER_AGENT = "cyclerfinder-kernel-freshness-check/1"

_SIZE_UNITS = {"K": 1024, "M": 1024**2, "G": 1024**3}

# Apache mod_autoindex row, e.g.:
#   <a href="jup365.bsp">jup365.bsp</a>                     2021-03-14 08:29  1.1G
_ROW_RE = re.compile(
    r'<a href="([^"]+\.bsp)">[^<]*</a>\s+(\d{4}-\d{2}-\d{2} \d{2}:\d{2})\s+([\d.]+[KMG]?)'
)


@dataclass(frozen=True)
class KernelSpec:
    """One tracked kernel: what we pin, where it lives locally, where to check."""

    system: str
    # The filename this project currently relies on/is documented against.
    # This is the freshness-comparison baseline -- NOT a check that this
    # exact file is on disk (see local_path for that).
    pinned_filename: str
    # Where the kernel is expected on disk if fetched (informational only).
    local_path: Path
    naif_directory_url: str
    # Matches candidate filenames in the SAME numbered series as
    # pinned_filename. Group 1 must capture an integer version number used
    # for ordering. Deliberately anchored (^...$) so it does NOT match the
    # XL/extended-length navigation variants (e.g. "ura111xl-701.bsp",
    # "jup310xl.bsp") -- those are a different kernel class the project has
    # already explicitly rejected (see scripts/install_uranian_spice.sh).
    series_pattern: re.Pattern[str]


KERNEL_MANIFEST: dict[str, KernelSpec] = {
    "jupiter": KernelSpec(
        system="jupiter",
        pinned_filename="jup365.bsp",
        local_path=Path("~/dev/references/kernels/jup365.bsp"),
        naif_directory_url=NAIF_SATELLITES_DIR,
        series_pattern=re.compile(r"^jup(\d{3})\.bsp$"),
    ),
    "uranus": KernelSpec(
        system="uranus",
        pinned_filename="ura111.bsp",
        # Lives under the GMAT install tree, not ~/dev/references/kernels --
        # see scripts/install_uranian_spice.sh for why.
        local_path=Path("~/GMAT/R2022a/data/planetary_ephem/spk/uranian/ura111.bsp"),
        naif_directory_url=NAIF_SATELLITES_DIR,
        # ura184_part-N.bsp is a newer multi-part release of the same
        # classical-satellite series (184 > 111 sorts correctly as "newer").
        series_pattern=re.compile(r"^ura(\d{3})(?:_part-\d+)?\.bsp$"),
    ),
    "pluto": KernelSpec(
        system="pluto",
        pinned_filename="plu060.bsp",
        local_path=Path("~/dev/references/kernels/plu060.bsp"),
        naif_directory_url=NAIF_SATELLITES_DIR,
        series_pattern=re.compile(r"^plu(\d{3})\.bsp$"),
    ),
}


def parse_directory_listing(html: str) -> dict[str, int]:
    """Parse an Apache mod_autoindex HTML directory listing into {filename: size_bytes}.

    Only ``.bsp`` rows are captured (``.cmt``/``.mrg`` companion files are
    irrelevant to freshness). Pure string parsing -- no network I/O.
    """
    listing: dict[str, int] = {}
    for filename, _mtime, size_str in _ROW_RE.findall(html):
        listing[filename] = _parse_size(size_str)
    return listing


def _parse_size(size_str: str) -> int:
    size_str = size_str.strip()
    if size_str and size_str[-1] in _SIZE_UNITS:
        return round(float(size_str[:-1]) * _SIZE_UNITS[size_str[-1]])
    return round(float(size_str))


def fetch_directory_listing(url: str, timeout: float = 30.0) -> str:
    """Fetch a NAIF directory listing page. Network I/O -- not called from tests."""
    req = urllib.request.Request(url, headers={"User-Agent": _USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        body: bytes = resp.read()
        return body.decode("utf-8", errors="replace")


@dataclass(frozen=True)
class KernelCheckResult:
    system: str
    status: str  # "OK" | "STALE" | "MISSING"
    pinned_filename: str
    pinned_version: int
    remote_filename: str | None
    remote_version: int | None
    remote_size_bytes: int | None
    local_present: bool
    local_size_bytes: int | None


def _series_version(spec: KernelSpec, filename: str) -> int:
    match = spec.series_pattern.match(filename)
    if match is None:
        raise ValueError(
            f"{spec.system}: pinned filename {filename!r} does not match its own "
            f"series_pattern {spec.series_pattern.pattern!r}"
        )
    return int(match.group(1))


def _current_in_listing(spec: KernelSpec, listing: dict[str, int]) -> tuple[str, int, int] | None:
    """Return (filename, version, size_bytes) for the highest-version match, or None."""
    candidates = [
        (filename, int(m.group(1)), size)
        for filename, size in listing.items()
        if (m := spec.series_pattern.match(filename)) is not None
    ]
    if not candidates:
        return None
    max_version = max(version for _, version, _ in candidates)
    top = [c for c in candidates if c[1] == max_version]
    # Multi-part releases (e.g. ura184_part-1/2/3) share a version; report
    # the largest single file as the representative for display.
    top.sort(key=lambda c: c[2], reverse=True)
    filename, version, size = top[0]
    return filename, version, size


def check_kernel(spec: KernelSpec, listing: dict[str, int] | None) -> KernelCheckResult:
    """Compare a tracked kernel's pinned version against a (possibly absent) remote listing.

    ``listing=None`` means the remote fetch itself failed (network error);
    callers should treat that as a hard error, not silently pass it here.
    """
    pinned_version = _series_version(spec, spec.pinned_filename)

    remote_filename: str | None = None
    remote_version: int | None = None
    remote_size: int | None = None
    if listing is not None:
        found = _current_in_listing(spec, listing)
        if found is not None:
            remote_filename, remote_version, remote_size = found

    local_path = spec.local_path.expanduser()
    local_present = local_path.exists()
    local_size = local_path.stat().st_size if local_present else None

    if remote_version is not None and remote_version > pinned_version:
        status = "STALE"
    elif not local_present:
        status = "MISSING"
    else:
        status = "OK"

    return KernelCheckResult(
        system=spec.system,
        status=status,
        pinned_filename=spec.pinned_filename,
        pinned_version=pinned_version,
        remote_filename=remote_filename,
        remote_version=remote_version,
        remote_size_bytes=remote_size,
        local_present=local_present,
        local_size_bytes=local_size,
    )


def format_result(result: KernelCheckResult) -> str:
    if result.status == "STALE":
        assert result.remote_filename is not None and result.remote_size_bytes is not None
        delta = None
        if result.local_size_bytes is not None:
            delta = result.remote_size_bytes - result.local_size_bytes
        delta_str = f", delta {delta:+d} bytes" if delta is not None else ""
        return (
            f"{result.system}: STALE -- pinned {result.pinned_filename} "
            f"(v{result.pinned_version}), current upstream is "
            f"{result.remote_filename} (v{result.remote_version}, "
            f"{result.remote_size_bytes} bytes{delta_str})"
        )
    if result.status == "MISSING":
        return (
            f"{result.system}: MISSING -- pinned {result.pinned_filename} is current "
            f"but not present locally at the documented path"
        )
    return f"{result.system}: OK -- {result.pinned_filename} is current and present locally"


def download_kernel(spec: KernelSpec, result: KernelCheckResult) -> Path:
    """Download the flagged-stale remote kernel into the local kernel dir.

    Only called when ``--fetch`` is explicitly passed. Writes the NEW
    filename alongside (not over) the old pinned file -- this script never
    edits KERNEL_MANIFEST for you.
    """
    if result.remote_filename is None:
        raise ValueError(f"{spec.system}: no remote filename to fetch (not STALE?)")
    dest_dir = spec.local_path.expanduser().parent
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / result.remote_filename
    url = spec.naif_directory_url + result.remote_filename
    print(f"{spec.system}: fetching {url} -> {dest}")
    urllib.request.urlretrieve(url, dest)
    return dest


def run(systems: list[str], *, do_fetch: bool) -> int:
    exit_code = 0
    for system in systems:
        spec = KERNEL_MANIFEST[system]
        try:
            html = fetch_directory_listing(spec.naif_directory_url)
        except (urllib.error.URLError, OSError) as exc:
            print(f"{system}: ERROR fetching {spec.naif_directory_url} ({exc})")
            exit_code = max(exit_code, 2)
            continue
        listing = parse_directory_listing(html)
        result = check_kernel(spec, listing)
        print(format_result(result))
        if result.status == "STALE":
            exit_code = max(exit_code, 1)
            if do_fetch:
                download_kernel(spec, result)
    return exit_code


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else None)
    parser.add_argument(
        "--system",
        action="append",
        choices=sorted(KERNEL_MANIFEST),
        help="restrict to one or more tracked systems (default: all)",
    )
    parser.add_argument(
        "--fetch",
        action="store_true",
        help="download any STALE kernel's current upstream version (never runs unattended)",
    )
    args = parser.parse_args(argv)
    systems = args.system or sorted(KERNEL_MANIFEST)
    return run(systems, do_fetch=args.fetch)


if __name__ == "__main__":
    sys.exit(main())
