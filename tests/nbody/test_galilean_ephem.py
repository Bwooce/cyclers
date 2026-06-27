from pathlib import Path

from cyclerfinder.verify.spice_kernels import ensure_jup365_kernel


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
