"""N-body Phase A: t_sec (TDB-J2000) -> ET reuses #129's LSK (plan Phase A; design §0).

GOLDEN DISCIPLINE: EXPECTED side = SPICE's own str2et over the J2000 epoch (the
LSK-defined time scale), not a value our code computed. This catches the TDB vs
UTC/TT trap (design §0: the "~64.184 s with drift" failure).

NOTE (deviation from the plan's test sketch): the plan wrote the J2000 anchor as
``str2et("2000-01-01T12:00:00 TDB")``, but the installed spiceypy rejects the ISO
"T" delimiter when a scale tag is present. The SPICE-parseable TDB form is
``"2000-01-01 12:00:00 TDB"`` (space delimiter); it returns ET = 0.0 (correct).
The trap-guard below pins that the tagless ISO form is read as UTC (+64.184 s).
"""

from __future__ import annotations

import pytest

spiceypy = pytest.importorskip("spiceypy")

from cyclerfinder.nbody.convert import t_sec_to_et  # noqa: E402
from cyclerfinder.verify.spice_kernels import ensure_leapseconds_kernel  # noqa: E402


@pytest.mark.slow
def test_t_sec_zero_is_j2000_et() -> None:
    spiceypy.furnsh(ensure_leapseconds_kernel())
    try:
        et_j2000 = spiceypy.str2et("2000-01-01 12:00:00 TDB")
        assert t_sec_to_et(0.0) == pytest.approx(et_j2000, abs=1e-3)
    finally:
        spiceypy.kclear()


@pytest.mark.slow
def test_t_sec_is_tdb_not_utc() -> None:
    """A one-year offset stays a clean TDB second count (no leap-second jump)."""
    spiceypy.furnsh(ensure_leapseconds_kernel())
    try:
        one_year = 365.25 * 86400.0
        et = spiceypy.str2et("2000-01-01 12:00:00 TDB") + one_year
        assert t_sec_to_et(one_year) == pytest.approx(et, abs=1e-3)
    finally:
        spiceypy.kclear()


@pytest.mark.slow
def test_tagless_iso_epoch_is_the_utc_trap() -> None:
    """Demonstrate the design §0 trap: tagless J2000 ISO string reads as UTC.

    t_sec_to_et(0) must be the TDB anchor (~0), NOT the +64.184 s that the naive
    UTC interpretation of the same calendar instant produces.
    """
    spiceypy.furnsh(ensure_leapseconds_kernel())
    try:
        utc_trap_et = spiceypy.str2et("2000-01-01T12:00:00")  # tagless -> UTC
        assert utc_trap_et == pytest.approx(64.184, abs=1e-2)  # the trap
        assert abs(t_sec_to_et(0.0) - utc_trap_et) == pytest.approx(64.184, abs=1e-2)
    finally:
        spiceypy.kclear()


@pytest.mark.slow
def test_t_sec_to_et_round_trip_is_linear() -> None:
    """t_sec_to_et(t) - t_sec_to_et(0) == t for non-round t (leap-second guard).

    SPICE ET is itself TDB-seconds-past-J2000, so the mapping must be the exact
    identity offset: the difference of two conversions reproduces the raw t_sec to
    relative precision. Any leap-second / TDB↔UTC contamination would break this
    linearity (the difference would pick up a non-constant jump).
    """
    base = t_sec_to_et(0.0)
    for t_sec in (1.0, 123.456, 86_400.0 + 0.5, 3.156e7 + 17.0, -5_000_000.25):
        assert t_sec_to_et(t_sec) - base == pytest.approx(t_sec, rel=1e-9)


@pytest.mark.slow
def test_t_sec_to_et_does_not_unload_other_kernels() -> None:
    """t_sec_to_et must NOT contaminate a caller's SPICE pool (kernel isolation).

    The earlier implementation called ``kclear()`` in its finally block, which
    unloads ALL furnished kernels — silently dropping any BSP a caller had loaded.
    Furnish an unrelated kernel (the DE440 BSP astropy already caches), call
    t_sec_to_et, and assert the other kernel survives via ``ktotal``.
    """
    from cyclerfinder.verify.spice_kernels import astropy_de440_bsp_path

    other_kernel = astropy_de440_bsp_path()
    spiceypy.furnsh(other_kernel)
    try:
        before = spiceypy.ktotal("ALL")
        _ = t_sec_to_et(123.0)
        after = spiceypy.ktotal("ALL")
        # The LSK furnished+unloaded inside the call nets to zero change; the
        # caller's BSP is still loaded (the old kclear() would have dropped it).
        assert after == before
        loaded = {spiceypy.kdata(i, "ALL")[0] for i in range(after)}
        assert other_kernel in loaded
    finally:
        spiceypy.kclear()
