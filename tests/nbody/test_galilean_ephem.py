from pathlib import Path

import pytest

from cyclerfinder.verify.spice_kernels import ensure_jup365_kernel

try:
    _KERNEL: str | None = ensure_jup365_kernel()
except Exception:  # jup365.bsp is local-only (~50 MB, absent in CI) -> skip, don't error
    _KERNEL = None

pytestmark = pytest.mark.skipif(_KERNEL is None, reason="JUP365 kernel not furnished (local-only)")


def test_jup365_kernel_path_exists() -> None:
    p = ensure_jup365_kernel()
    assert Path(p).is_file()
    assert p.endswith("jup365.bsp")


def test_galilean_spice_backend_matches_spkezr() -> None:
    import spiceypy as s

    from cyclerfinder.core.ephemeris import Ephemeris

    eph = Ephemeris(center="Jupiter", model="spice")
    for et in (0.0, 7.05 * 86400.0, 28.22 * 86400.0):
        r, v = eph.state("Europa", et)  # km, km/s, J2000-equatorial, Jupiter-centred
        ref, _ = s.spkezr("502", et, "J2000", "NONE", "599")
        assert max(abs(a - b) for a, b in zip(r, ref[:3], strict=True)) < 1.0
        assert max(abs(a - b) for a, b in zip(v, ref[3:], strict=True)) < 1e-6


def test_galilean_spice_io_period_sanity() -> None:
    from cyclerfinder.core.ephemeris import Ephemeris

    eph = Ephemeris(center="Jupiter", model="spice")
    r0, _ = eph.state("Io", 0.0)
    r1, _ = eph.state("Io", 1.769 * 86400.0)
    # Io's orbital period is ~1.769 days; after one period it should be near start
    assert max(abs(a - b) for a, b in zip(r0, r1, strict=True)) < 5000.0


def test_centred_spice_backend_survives_external_kclear() -> None:
    """#824 regression: a same-process ``spice.kclear()`` call elsewhere must
    not leave ``_CentredSpiceBackend`` silently unable to serve states.

    Same bug class as ``nbody.jovian.JovianEphemeris``'s own former
    ``_FURNISHED`` cache (see that module's fix and
    ``tests/nbody/test_jovian.py::
    test_jovian_ephemeris_survives_external_kclear``): the former
    ``_JUP365_FURNISHED`` module-level flag in ``core/ephemeris.py`` had no
    way to learn that ``spice.kclear()`` (called by several other modules
    in this codebase) had wiped SPICE's entire pool, so a second
    construction skipped re-furnishing and the next ``spkezr`` call raised
    ``SpiceNOLOADEDFILES``.
    """
    import spiceypy as s

    from cyclerfinder.core.ephemeris import Ephemeris

    eph1 = Ephemeris(center="Jupiter", model="spice")
    r1, _ = eph1.state("Europa", 0.0)

    s.kclear()  # simulates an unrelated same-process caller

    eph2 = Ephemeris(center="Jupiter", model="spice")  # must re-furnish, not trust a cache
    r2, _ = eph2.state("Europa", 0.0)
    assert max(abs(a - b) for a, b in zip(r1, r2, strict=True)) < 1e-9
