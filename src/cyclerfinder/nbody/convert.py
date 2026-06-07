"""§0 time/frame conversion glue for the n-body harness (design §0).

The harness ``t_sec`` axis is **TDB seconds since J2000** — identical to the
``Ephemeris("astropy")`` backend's convention (``core/ephemeris.py:44-47,262-263``:
``t_sec = 0`` is 2000-01-01T12:00:00 TDB, ``Time(..., scale="tdb")``). This module
is the single boundary where that axis and frame are translated for tools that do
not share ``core/ephemeris.py``'s machinery.

Two concerns, both reusing existing project machinery rather than reinventing it
(design §0, the "~64.184 s with drift" TDB↔TT/UTC trap):

* **Time** — :func:`t_sec_to_et` maps our TDB-J2000-seconds to SPICE ET for
  SPICE-native tools. SPICE ET *is* TDB seconds past J2000, so the conversion is
  the identity ``str2et("J2000 TDB") + t_sec`` — but we route through #129's LSK
  (``verify/spice_kernels.ensure_leapseconds_kernel``) so the J2000 anchor is the
  LSK-defined one and we never accidentally go through UTC/TT.
* **Frame** — :func:`icrs_eq_to_ecliptic` / :func:`ecliptic_to_icrs_eq` rotate
  between ICRS-equatorial (what any n-body / SPICE tool returns by default) and
  the heliocentric J2000-ecliptic frame the rest of cyclerfinder uses, using the
  *same* obliquity constant ``core/ephemeris._J2000_OBLIQUITY_RAD`` (imported, not
  re-literal'd). A wrong obliquity reads as a fake out-of-plane V∞ component.
"""

from __future__ import annotations

from math import cos, sin

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.ephemeris import _J2000_OBLIQUITY_RAD

Vec3 = NDArray[np.float64]

# The J2000 epoch string in the TDB scale. SPICE ET is defined as ephemeris
# seconds past this very epoch, so str2et of it is ~0 and our t_sec axis maps to
# ET by simple addition. Spelled with the explicit "TDB" tag AND a space
# delimiter, the form spiceypy's str2et accepts: dropping the "TDB" tag (or using
# the ISO "T" delimiter, which str2et rejects with a tag) makes str2et read the
# epoch as UTC and return ET = +64.184 s — the exact design §0 trap. The TDB form
# returns 0.0, which is correct.
_J2000_TDB_STR = "2000-01-01 12:00:00 TDB"


def t_sec_to_et(t_sec: float) -> float:
    """Convert harness ``t_sec`` (TDB seconds since J2000) to SPICE ET seconds.

    Loads #129's leapseconds kernel (``ensure_leapseconds_kernel``) and anchors
    on ``str2et`` of the J2000 TDB epoch, then adds ``t_sec``. Because SPICE ET is
    itself TDB-seconds-past-J2000, this is the identity to within ~ms — but it is
    computed against the LSK-defined time scale, not assumed, so a TDB↔TT/UTC slip
    cannot creep in.

    Requires ``spiceypy`` (the ``validation`` extra). Imports it lazily so the
    fast suite need not have it installed.
    """
    import spiceypy

    from cyclerfinder.verify.spice_kernels import ensure_leapseconds_kernel

    spiceypy.furnsh(ensure_leapseconds_kernel())
    try:
        et_j2000 = float(spiceypy.str2et(_J2000_TDB_STR))
    finally:
        spiceypy.kclear()
    return et_j2000 + float(t_sec)


def _r_icrs_to_ecl() -> NDArray[np.float64]:
    """ICRS-equatorial -> J2000-ecliptic rotation: R_x(-eps).

    Identical to ``core/ephemeris._AstropyBackend``'s pre-computed matrix
    (``ephemeris.py:241-250``); built from the SAME ``_J2000_OBLIQUITY_RAD``.
    """
    eps = _J2000_OBLIQUITY_RAD
    return np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, cos(eps), sin(eps)],
            [0.0, -sin(eps), cos(eps)],
        ],
        dtype=np.float64,
    )


def icrs_eq_to_ecliptic(vec3: Vec3) -> Vec3:
    """Rotate an ICRS-equatorial 3-vector into the J2000-ecliptic frame."""
    return np.asarray(_r_icrs_to_ecl() @ np.asarray(vec3, dtype=np.float64), dtype=np.float64)


def ecliptic_to_icrs_eq(vec3: Vec3) -> Vec3:
    """Rotate a J2000-ecliptic 3-vector back into the ICRS-equatorial frame."""
    return np.asarray(_r_icrs_to_ecl().T @ np.asarray(vec3, dtype=np.float64), dtype=np.float64)


__all__ = ["ecliptic_to_icrs_eq", "icrs_eq_to_ecliptic", "t_sec_to_et"]
