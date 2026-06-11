"""Keystone external epoch anchor for the astropy ephemeris backend (#198).

Every other gate in the project consumes ``Ephemeris("astropy")`` states
through the SAME ``t_sec`` axis, so a coherent shift of the whole simulation
clock cancels out of every internal comparison (the false-consensus failure
mode). Even the #129 SPICE/Horizons cross-checks replicated the backend's own
epoch construction (``_backend_tdb_jd`` in ``test_ephemeris_crosscheck.py``),
so they validated the kernel *reader*, not the epoch. Exactly this happened:
until 2026-06-11 the backend built its epoch from a UTC datetime's POSIX
timestamp relabeled ``scale="tdb"``, landing ``t_sec=0`` at J2000(TDB)
+ 64.184 s — ~1944 km of Earth along-track error at every absolute epoch —
and every test passed.

This module is the per-interface anchor that breaks that consensus: the
EXPECTED side is transcribed from JPL Horizons (provenance below), and the
``t_sec`` for each anchor is computed by the *definition* of the axis —
``t_sec = (JD_TDB - 2451545.0) * 86400`` — never by replicating the backend's
internal ``Time`` construction. The tolerance (100 km) is generous for kernel
differences (DE440 ours vs DE441 Horizons: sub-km for Earth here) but far
below the ~1944 km signature of a ~64 s epoch shift, so this test FAILS on the
pre-fix code (demonstrated against the pre-fix tree on 2026-06-11) and pins
the epoch convention for good.

Horizons provenance — retrieved 2026-06-11 from the JPL Horizons API,
https://ssd.jpl.nasa.gov/api/horizons.api, with:
  COMMAND='399' (Earth), CENTER='500@10' (Sun center), EPHEM_TYPE='VECTORS',
  REF_PLANE='ECLIPTIC', REF_SYSTEM='J2000', OUT_UNITS='KM-S', VEC_TABLE='2',
  VEC_CORR='NONE' (geometric states), TLIST as the JD(TDB) noted per anchor.
Horizons reported: Target body Earth (399) {source: DE441}, Center body
Sun (10) {source: DE441}, "Reference frame : Ecliptic of J2000.0", epochs
stated in TDB.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.ephemeris import Ephemeris

# Seconds per TDB day; with the J2000 JD this DEFINES the t_sec axis.
_J2000_JD_TDB = 2451545.0
_SEC_PER_DAY = 86400.0

# (label, JD_TDB, r_km, v_km_s) — transcribed verbatim from the Horizons
# output (query details in the module docstring).
_ANCHORS: list[tuple[str, float, tuple[float, float, float], tuple[float, float, float]]] = [
    (
        # 2000-01-01 12:00:00.0000 TDB — the J2000 epoch itself, t_sec = 0.
        "J2000-TDB",
        2451545.0,
        (-2.649903367743050e07, 1.446972967925493e08, -6.111494259536266e02),
        (-2.979426007043741e01, -5.469294939770602e00, 1.817836785027449e-04),
    ),
    (
        # 2026-01-01 00:00:00.0000 TDB — a modern epoch, exercises the
        # cumulative (leap-second-free) span of the axis.
        "2026-01-01-TDB",
        2461041.5,
        (-2.607213844816194e07, 1.447746738210197e08, -8.892861905746162e03),
        (-2.978893116564825e01, -5.396270536820538e00, 4.106639513727917e-04),
    ),
]

# Position gate: DE440-vs-DE441 for Earth at these epochs is sub-km; a ~64 s
# epoch shift is ~1944 km (Earth at ~29.8 km/s). 100 km admits kernel-version
# noise with huge headroom yet catches any time-scale slip above ~3 s.
_POS_GATE_KM = 100.0
# Velocity gate: same logic — DE-version noise is < 0.01 mm/s, while a 64 s
# shift moves Earth's velocity vector by ~13 mm/s (centripetal turn).
_VEL_GATE_MM_S = 10.0


@pytest.mark.parametrize(
    ("label", "jd_tdb", "r_expected", "v_expected"),
    _ANCHORS,
    ids=[a[0] for a in _ANCHORS],
)
def test_earth_state_anchored_to_horizons_at_exact_tdb_epoch(
    label: str,
    jd_tdb: float,
    r_expected: tuple[float, float, float],
    v_expected: tuple[float, float, float],
) -> None:
    """Earth's heliocentric state at an exact TDB instant matches Horizons.

    ``t_sec`` comes from the axis DEFINITION (TDB seconds since JD 2451545.0
    TDB), not from the backend's internals — so this anchors the epoch
    convention itself, not merely the kernel reader.
    """
    pytest.importorskip("astropy")
    eph = Ephemeris("astropy")
    t_sec = (jd_tdb - _J2000_JD_TDB) * _SEC_PER_DAY
    r_ours, v_ours = eph.state("E", t_sec)

    d_pos = float(np.linalg.norm(r_ours - np.asarray(r_expected)))
    d_vel_mm_s = float(np.linalg.norm(v_ours - np.asarray(v_expected))) * 1.0e6
    assert d_pos < _POS_GATE_KM, (
        f"{label}: Earth position {d_pos:.3f} km off the Horizons anchor "
        f"(gate {_POS_GATE_KM} km). A ~64 s epoch shift shows up as ~1944 km — "
        "if d_pos is in that ballpark, the t_sec=0 epoch has drifted off "
        "J2000(TDB) again (see #198)."
    )
    assert d_vel_mm_s < _VEL_GATE_MM_S, (
        f"{label}: Earth velocity {d_vel_mm_s:.4f} mm/s off the Horizons "
        f"anchor (gate {_VEL_GATE_MM_S} mm/s)"
    )


def test_batched_states_match_scalar_at_anchor_epochs() -> None:
    """The batched ``states()`` path shares the scalar path's epoch exactly."""
    pytest.importorskip("astropy")
    eph = Ephemeris("astropy", cache=False)
    epochs = [(a[1] - _J2000_JD_TDB) * _SEC_PER_DAY for a in _ANCHORS]
    bodies = ["E"] * len(epochs)
    batched = eph.states(bodies, epochs)
    for (rb, vb), t_sec in zip(batched, epochs, strict=True):
        rs, vs = eph.state("E", t_sec)
        assert float(np.linalg.norm(rb - rs)) == 0.0
        assert float(np.linalg.norm(vb - vs)) == 0.0
