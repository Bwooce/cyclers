"""Tests for the NAIF kernel-freshness checker (#514).

NO network access. All HTTP is mocked: :data:`_FIXTURE_LISTING_HTML` is a
verbatim excerpt of NAIF's real Apache mod_autoindex directory listing
(``https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/``,
captured 2026-07-01 while drafting this checker) so the parsing regex is
exercised against the real row format, not an invented one.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

import scripts.check_kernel_freshness as ckf
from scripts.check_kernel_freshness import (
    KernelSpec,
    _current_in_listing,
    check_kernel,
    format_result,
    parse_directory_listing,
    run,
)

# Verbatim excerpt (2026-07-01) covering the three tracked systems plus the
# XL/extended-length variants and non-.bsp companion files that must be
# ignored or excluded.
_FIXTURE_LISTING_HTML = """\
<img src="/icons/unknown.gif" alt="[   ]"> <a href="jup347.bsp">jup347.bsp</a>                     2025-05-15 11:58  879M
<img src="/icons/unknown.gif" alt="[   ]"> <a href="jup348.bsp">jup348.bsp</a>                     2026-03-18 10:38   57M
<img src="/icons/unknown.gif" alt="[   ]"> <a href="jup349.bsp">jup349.bsp</a>                     2026-04-13 11:15   93M
<img src="/icons/unknown.gif" alt="[   ]"> <a href="jup349.mrg">jup349.mrg</a>                     2026-04-13 11:15  546
<img src="/icons/unknown.gif" alt="[   ]"> <a href="jup365.bsp">jup365.bsp</a>                     2021-03-14 08:29  1.1G
<img src="/icons/unknown.gif" alt="[   ]"> <a href="jup365.cmt">jup365.cmt</a>                     2021-03-14 08:29   68K
<img src="/icons/unknown.gif" alt="[   ]"> <a href="mar099.bsp">mar099.bsp</a>                     2021-04-13 08:00  611M
<img src="/icons/unknown.gif" alt="[   ]"> <a href="plu060.bsp">plu060.bsp</a>                     2024-04-03 07:44  129M
<img src="/icons/unknown.gif" alt="[   ]"> <a href="plu060.cmt">plu060.cmt</a>                     2024-04-03 07:44   63K
<img src="/icons/unknown.gif" alt="[   ]"> <a href="ura111xl-701.bsp">ura111xl-701.bsp</a>               2022-04-10 13:29  2.0G
<img src="/icons/unknown.gif" alt="[   ]"> <a href="ura111xl-702.bsp">ura111xl-702.bsp</a>               2022-04-10 13:29  828M
<img src="/icons/unknown.gif" alt="[   ]"> <a href="ura116xl.bsp">ura116xl.bsp</a>                   2022-04-10 13:30  659M
<img src="/icons/unknown.gif" alt="[   ]"> <a href="ura184_part-1.bsp">ura184_part-1.bsp</a>              2025-09-26 06:42  1.9G
<img src="/icons/unknown.gif" alt="[   ]"> <a href="ura184_part-2.bsp">ura184_part-2.bsp</a>              2025-09-26 06:42  1.9G
<img src="/icons/unknown.gif" alt="[   ]"> <a href="ura184_part-2.mrg">ura184_part-2.mrg</a>              2025-09-26 06:42  480
<img src="/icons/unknown.gif" alt="[   ]"> <a href="ura184_part-3.bsp">ura184_part-3.bsp</a>              2025-09-26 06:43  369M
<img src="/icons/unknown.gif" alt="[   ]"> <a href="ura184_part-3.cmt">ura184_part-3.cmt</a>              2025-09-26 06:43   63K
"""


# --- parse_directory_listing ---------------------------------------------------


def test_parse_directory_listing_extracts_bsp_filenames_and_sizes() -> None:
    listing = parse_directory_listing(_FIXTURE_LISTING_HTML)
    assert listing["jup365.bsp"] == round(1.1 * 1024**3)
    assert listing["plu060.bsp"] == 129 * 1024**2
    assert listing["ura184_part-3.bsp"] == 369 * 1024**2


def test_parse_directory_listing_ignores_non_bsp_companion_files() -> None:
    listing = parse_directory_listing(_FIXTURE_LISTING_HTML)
    assert "jup349.mrg" not in listing
    assert "jup365.cmt" not in listing
    assert "plu060.cmt" not in listing
    assert "ura184_part-2.mrg" not in listing
    assert "ura184_part-3.cmt" not in listing


def test_parse_directory_listing_captures_every_series_relevant_to_tracked_systems() -> None:
    listing = parse_directory_listing(_FIXTURE_LISTING_HTML)
    for name in ("jup347.bsp", "jup348.bsp", "jup349.bsp", "jup365.bsp", "plu060.bsp"):
        assert name in listing
    assert "mar099.bsp" in listing  # unrelated series, present but never matched below


# --- series_pattern excludes XL variants ---------------------------------------


def test_uranus_series_pattern_excludes_xl_variants() -> None:
    listing = parse_directory_listing(_FIXTURE_LISTING_HTML)
    spec = ckf.KERNEL_MANIFEST["uranus"]
    found = _current_in_listing(spec, listing)
    assert found is not None
    filename, version, _size = found
    # The XL files (ura111xl-701.bsp, ura116xl.bsp) sort higher lexically
    # in places but must never be picked -- they are a rejected kernel class.
    assert "xl" not in filename
    assert filename.startswith("ura184_part-")
    assert version == 184


def test_jupiter_series_pattern_matches_only_numbered_jup3xx() -> None:
    listing = parse_directory_listing(_FIXTURE_LISTING_HTML)
    spec = ckf.KERNEL_MANIFEST["jupiter"]
    found = _current_in_listing(spec, listing)
    assert found == ("jup365.bsp", 365, round(1.1 * 1024**3))


# --- check_kernel: freshness x local-presence matrix ----------------------------


def _spec(system: str, pinned: str, pattern: str, local_path: Path) -> KernelSpec:
    return KernelSpec(
        system=system,
        pinned_filename=pinned,
        local_path=local_path,
        naif_directory_url=ckf.NAIF_SATELLITES_DIR,
        series_pattern=re.compile(pattern),
    )


def test_check_kernel_ok_when_pinned_is_current_and_locally_present(tmp_path: Path) -> None:
    local = tmp_path / "plu060.bsp"
    local.write_bytes(b"x" * 100)
    listing = parse_directory_listing(_FIXTURE_LISTING_HTML)
    spec = _spec("pluto", "plu060.bsp", r"^plu(\d{3})\.bsp$", local)
    result = check_kernel(spec, listing)
    assert result.status == "OK"
    assert result.pinned_version == 60
    assert result.remote_version == 60
    assert result.local_present is True


def test_check_kernel_missing_when_pinned_is_current_but_not_on_disk(tmp_path: Path) -> None:
    local = tmp_path / "plu060.bsp"  # never written
    listing = parse_directory_listing(_FIXTURE_LISTING_HTML)
    spec = _spec("pluto", "plu060.bsp", r"^plu(\d{3})\.bsp$", local)
    result = check_kernel(spec, listing)
    assert result.status == "MISSING"
    assert result.local_present is False


def test_check_kernel_stale_when_remote_has_higher_version_regardless_of_local_presence(
    tmp_path: Path,
) -> None:
    """This is the #506/#510 scenario: a pinned version has been superseded upstream.

    STALE must be reported whether or not the (old) file happens to be on
    disk locally -- freshness is a manifest-vs-remote comparison, not a
    disk-presence comparison, so this also works with zero local files (CI).
    """
    listing = parse_directory_listing(_FIXTURE_LISTING_HTML)
    spec = _spec(
        "uranus", "ura111.bsp", r"^ura(\d{3})(?:_part-\d+)?\.bsp$", tmp_path / "ura111.bsp"
    )
    result = check_kernel(spec, listing)
    assert result.status == "STALE"
    assert result.pinned_version == 111
    assert result.remote_filename == "ura184_part-1.bsp"
    assert result.remote_version == 184

    # Same result even with a local file present -- freshness ignores it.
    (tmp_path / "ura111.bsp").write_bytes(b"x" * 100)
    result_with_local = check_kernel(spec, listing)
    assert result_with_local.status == "STALE"


def test_check_kernel_stale_size_delta_is_remote_minus_local(tmp_path: Path) -> None:
    local = tmp_path / "ura111.bsp"
    local.write_bytes(b"x" * 1000)
    listing = parse_directory_listing(_FIXTURE_LISTING_HTML)
    spec = _spec("uranus", "ura111.bsp", r"^ura(\d{3})(?:_part-\d+)?\.bsp$", local)
    result = check_kernel(spec, listing)
    assert result.status == "STALE"
    assert result.remote_size_bytes == round(1.9 * 1024**3)
    assert result.local_size_bytes == 1000


def test_check_kernel_no_remote_candidates_falls_back_to_local_presence(tmp_path: Path) -> None:
    """If the remote listing has nothing matching the series at all, freshness
    can't be assessed as STALE -- fall through to local-presence reporting."""
    local = tmp_path / "nep097.bsp"
    local.write_bytes(b"x")
    spec = _spec("neptune", "nep097.bsp", r"^nep(\d{3})\.bsp$", local)
    result = check_kernel(spec, listing={})  # empty remote listing
    assert result.status == "OK"
    assert result.remote_filename is None


def test_check_kernel_pinned_filename_must_match_its_own_series_pattern(tmp_path: Path) -> None:
    spec = _spec("bogus", "not-a-match.bsp", r"^plu(\d{3})\.bsp$", tmp_path / "x.bsp")
    with pytest.raises(ValueError, match="does not match its own"):
        check_kernel(spec, listing={})


# --- format_result ---------------------------------------------------------------


def test_format_result_ok() -> None:
    listing = parse_directory_listing(_FIXTURE_LISTING_HTML)
    spec = ckf.KERNEL_MANIFEST["jupiter"]
    result = check_kernel(spec, listing)
    # jup365.bsp is not present in the test env's real HOME by construction
    # here (spec uses the real ~/dev/references/kernels path) -- assert on
    # the text shape only, not on local presence which is environment-dependent.
    text = format_result(result)
    assert text.startswith("jupiter: ")
    assert "jup365.bsp" in text


def test_format_result_stale_mentions_both_filenames(tmp_path: Path) -> None:
    listing = parse_directory_listing(_FIXTURE_LISTING_HTML)
    spec = _spec(
        "uranus", "ura111.bsp", r"^ura(\d{3})(?:_part-\d+)?\.bsp$", tmp_path / "ura111.bsp"
    )
    result = check_kernel(spec, listing)
    text = format_result(result)
    assert "STALE" in text
    assert "ura111.bsp" in text
    assert "ura184_part-1.bsp" in text


def test_format_result_missing() -> None:
    listing: dict[str, int] = {"plu060.bsp": 129 * 1024**2}
    spec = _spec("pluto", "plu060.bsp", r"^plu(\d{3})\.bsp$", Path("/nonexistent/plu060.bsp"))
    result = check_kernel(spec, listing)
    text = format_result(result)
    assert text.startswith("pluto: MISSING")


# --- run(): exit-code contract (mocked network) -----------------------------------


def test_run_exit_code_zero_when_all_ok_or_missing(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    fake_manifest = {
        "pluto": _spec("pluto", "plu060.bsp", r"^plu(\d{3})\.bsp$", tmp_path / "plu060.bsp"),
    }
    monkeypatch.setattr(ckf, "KERNEL_MANIFEST", fake_manifest)
    monkeypatch.setattr(
        ckf, "fetch_directory_listing", lambda url, timeout=30.0: _FIXTURE_LISTING_HTML
    )
    exit_code = run(["pluto"], do_fetch=False)
    assert exit_code == 0


def test_run_exit_code_nonzero_when_any_tracked_kernel_is_stale(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    fake_manifest = {
        "pluto": _spec("pluto", "plu060.bsp", r"^plu(\d{3})\.bsp$", tmp_path / "plu060.bsp"),
        "uranus": _spec(
            "uranus", "ura111.bsp", r"^ura(\d{3})(?:_part-\d+)?\.bsp$", tmp_path / "ura111.bsp"
        ),
    }
    monkeypatch.setattr(ckf, "KERNEL_MANIFEST", fake_manifest)
    monkeypatch.setattr(
        ckf, "fetch_directory_listing", lambda url, timeout=30.0: _FIXTURE_LISTING_HTML
    )
    exit_code = run(["pluto", "uranus"], do_fetch=False)
    assert exit_code == 1


def test_run_works_with_zero_local_files_present_ci_mode(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """The CI scenario: no kernel is on disk anywhere, only the remote listing exists.

    Freshness must still be assessable (STALE still fires for the superseded
    Uranus pin); presence-only absence must not itself fail the run.
    """
    fake_manifest = {
        "jupiter": _spec(
            "jupiter", "jup365.bsp", r"^jup(\d{3})\.bsp$", tmp_path / "does-not-exist.bsp"
        ),
    }
    monkeypatch.setattr(ckf, "KERNEL_MANIFEST", fake_manifest)
    monkeypatch.setattr(
        ckf, "fetch_directory_listing", lambda url, timeout=30.0: _FIXTURE_LISTING_HTML
    )
    exit_code = run(["jupiter"], do_fetch=False)
    assert exit_code == 0  # jup365 is still current upstream -> MISSING, not STALE


def test_run_reports_network_error_as_nonzero_without_crashing(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def _raise(url: str, timeout: float = 30.0) -> str:
        raise OSError("simulated network failure")

    monkeypatch.setattr(ckf, "fetch_directory_listing", _raise)
    exit_code = run(["pluto"], do_fetch=False)
    assert exit_code != 0
