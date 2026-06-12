"""FALSE-CONSENSUS DEFENCE HARNESS — fault-injection gate-rejection suite (#202).

Doctrine (memory ``orbit-closure-discipline`` "FALSE-CONSENSUS DEFENCE"
addendum; ``docs/notes/2026-06-11-project-review-results.md`` §"The
false-consensus doctrine", item 3): **agreement between N checks is only worth
what they do NOT share.** Three real incidents on this project proved the cost
of an undefended shared component:

* **#180** — three "independent" methods inherited one upstream ToF bug and
  agreed on a wrong answer.
* **#197** — ``crosscheck_leg`` read its Lambert endpoints from the artifact
  under test, so a wrong upstream position fed the SAME ``(r1, r2)`` to all
  three solvers, which happily agreed.
* **#198** — a 63 s UTC/TDB epoch-conversion offset shared between the primary
  path and its "independent" cross-check cancelled out of every internal
  comparison; every test passed for months.

The defence is fault injection (doctrine item 3): for each SHARED component a
validation gate depends on, deliberately POISON it and assert a SPECIFIC gate
FIRES. A poison the gauntlet survives is a *measured, undefended correlation*
between checks — recorded here as a FINDING, not papered over with a passing
test.

Every poison below is injected in-test via ``monkeypatch`` (source is never
mutated) and each test demonstrates BOTH halves of rejection power:

* CLEAN — the gate PASSES on the unpoisoned component (the gate is not vacuous);
* POISONED — the same gate FAILS on the injected fault (the gate has teeth).

Harness index (the executable doctrine, discoverable in one place)
------------------------------------------------------------------
========================  =======================================  ==========
Shared component          Gate asserted to fire                    Incident
========================  =======================================  ==========
Epoch convention          ``test_epoch_anchor.py`` Horizons        #198
                          absolute-epoch anchor
Frame handedness          ``test_ephemeris_inclined.py``           #199
(R_x(+inc) node sign)     DE440 orbit-normal anchor
Crosscheck endpoint       ``test_crosscheck.py``                   #197
(shared ``(r1,r2)``)      ``test_crosscheck_catches_poisoned_endpoint``
Earth-Moon mu             ``test_cr3bp.py``                        #212a
(sourced mass ratio)      ``test_earth_moon_mu_physical``
Signature / transit ToF   crosscheck endpoint-independence gate    #180
========================  =======================================  ==========

Cross-links to the rest of the doctrine
---------------------------------------
* CONSISTENCY-vs-INDEPENDENCE gate classification and the "shared with primary
  path:" evidence convention live in
  :mod:`cyclerfinder.data.validate`'s module docstring (#197).
* The POSITIVE-CONTROL convention (assert a method re-finds a KNOWN solution
  through the identical pipeline config before a negative is trusted — the rule
  that would have killed #180) is documented in
  ``docs/notes/2026-06-12-false-consensus-defence-harness.md`` and given a
  reusable helper in :mod:`cyclerfinder.verify.positive_control`.
* Per-interface external anchors (doctrine item 5) ARE the gates this harness
  poisons: ``test_epoch_anchor.py`` (epoch), ``test_ephemeris_inclined.py``
  DE440 orbit-normal (frame), ``test_cr3bp.py`` (sourced mu).
"""

from __future__ import annotations

import dataclasses
import math

import numpy as np
import pytest
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.core.ephemeris as ephemeris_mod
import cyclerfinder.core.satellites as satellites
from cyclerfinder.core.constants import PlanetData
from cyclerfinder.core.ephemeris import Ephemeris, _InclinedCircularBackend
from cyclerfinder.model.cycler import Cycler
from cyclerfinder.verify.crosscheck import POSITION_CONSISTENCY_TOL_KM, crosscheck_leg
from cyclerfinder.verify.real_closure import construct_real_ephemeris_cycler

# --------------------------------------------------------------------------- #
# Epoch / frame anchor constants — re-import the gate's own EXPECTED side so the
# poison tests assert against the SAME sourced anchor the gate uses (not a value
# this harness invented).
# --------------------------------------------------------------------------- #
_J2000_JD_TDB = 2451545.0
_SEC_PER_DAY = 86400.0

# J2000(TDB) Earth heliocentric state, transcribed verbatim from JPL Horizons in
# tests/verify/test_epoch_anchor.py (provenance there). Re-stated here so the
# epoch poison is checked against the external anchor, not the backend itself.
_HORIZONS_EARTH_J2000_R_KM = (
    -2.649903367743050e07,
    1.446972967925493e08,
    -6.111494259536266e02,
)
_EPOCH_POS_GATE_KM = 100.0  # test_epoch_anchor.py _POS_GATE_KM


# =========================================================================== #
# POISON 1 — EPOCH CONVENTION (#198)
#
# Shared component: the astropy backend's J2000(TDB) reference instant. Every
# gate consumes Ephemeris("astropy") through the SAME t_sec axis, so a coherent
# clock shift cancels out of every internal comparison — exactly how the 63 s
# offset survived. The external Horizons absolute-epoch anchor is the ONLY check
# that does not share the axis.
#
# Gate asserted: tests/verify/test_epoch_anchor.py
#   test_earth_state_anchored_to_horizons_at_exact_tdb_epoch (J2000-TDB anchor).
# =========================================================================== #


def _earth_pos_offset_from_horizons_at_j2000(eph: Ephemeris) -> float:
    """Replicate the epoch-anchor gate's J2000 position residual (km).

    Mirrors test_epoch_anchor.py exactly: t_sec from the axis DEFINITION
    (``(JD_TDB - 2451545.0) * 86400``), Earth state from the backend, compared
    to the sourced Horizons J2000(TDB) vector.
    """
    t_sec = (_J2000_JD_TDB - _J2000_JD_TDB) * _SEC_PER_DAY  # == 0.0 by definition
    r_ours, _v = eph.state("E", t_sec)
    return float(np.linalg.norm(r_ours - np.asarray(_HORIZONS_EARTH_J2000_R_KM)))


def test_epoch_anchor_passes_clean() -> None:
    """CLEAN: the epoch anchor passes on the unpoisoned backend (gate not vacuous)."""
    pytest.importorskip("astropy")
    eph = Ephemeris("astropy")
    d_pos = _earth_pos_offset_from_horizons_at_j2000(eph)
    assert d_pos < _EPOCH_POS_GATE_KM


def test_epoch_anchor_fires_on_shifted_j2000_reference(monkeypatch: pytest.MonkeyPatch) -> None:
    """POISON: shift the backend's J2000(TDB) reference by ~60 s -> anchor FAILS.

    The 63 s #198 bug survived precisely because nothing tripped; this pins that
    a coherent epoch shift now lands ~1800+ km off the external Horizons anchor
    (~29.8 km/s * 60 s ~ 1788 km), far above the 100 km gate. We poison the
    module-level reference JD the backend builds its TDB epoch from — the single
    point that defines where t_sec=0 lands — then build a fresh backend.
    """
    pytest.importorskip("astropy")
    shift_sec = 60.0
    poisoned_jd = _J2000_JD_TDB + shift_sec / _SEC_PER_DAY
    monkeypatch.setattr(ephemeris_mod, "_J2000_TDB_JD", poisoned_jd)

    eph = Ephemeris("astropy")
    d_pos = _earth_pos_offset_from_horizons_at_j2000(eph)
    # The gate (d_pos < 100 km) must now FAIL: the shift shows up as ~1800 km.
    assert d_pos > _EPOCH_POS_GATE_KM
    assert d_pos == pytest.approx(29.8 * shift_sec, rel=0.2)


# =========================================================================== #
# POISON 2 — FRAME HANDEDNESS (#199)
#
# Shared component: the inclined-circular backend's orbital-plane -> ecliptic
# rotation. The Standish ascending-node convention requires R_x(+inc); the
# pre-fix R_x(-inc) mirrored every orbital plane about the ecliptic, putting
# orbit normals 2*inc off DE440 (Venus 6.789 deg, Mars 3.699 deg). The internal
# n_hat formula test shared the same mirrored rotation and hid it; only the
# independent DE440 h = r x v anchor breaks the consensus.
#
# Gate asserted: tests/core/test_ephemeris_inclined.py
#   test_inclined_orbit_normal_anchored_to_de440.
# =========================================================================== #

_VENUS_INC_DEG = 3.39467605
_VENUS_LAN_DEG = 76.67984255


def _rotation_rx_minus_inc(planet: PlanetData) -> NDArray[np.float64]:
    """Pre-#199 mirrored rotation R_z(+lan) @ R_x(-inc).

    A faithful reconstruction of the bug: identical to the fixed
    ``_InclinedCircularBackend._rotation`` but with the R_x inclination sign
    flipped back to negative (the plane-mirror defect).
    """
    inc = math.radians(planet.inc_deg)
    lan = math.radians(planet.lan_deg)
    ci, si = math.cos(inc), math.sin(inc)
    cl, sl = math.cos(lan), math.sin(lan)
    rx = np.array(
        [[1.0, 0.0, 0.0], [0.0, ci, si], [0.0, -si, ci]],  # R_x(-inc): mirrored
        dtype=np.float64,
    )
    rz = np.array([[cl, -sl, 0.0], [sl, cl, 0.0], [0.0, 0.0, 1.0]], dtype=np.float64)
    return rz @ rx


def _de440_orbit_normal(body: str, t_sec: float) -> NDArray[np.float64]:
    r, v = Ephemeris(model="astropy").state(body, t_sec)
    h = np.cross(r, v)
    return np.asarray(h / np.linalg.norm(h), dtype=np.float64)


def _venus_backend() -> _InclinedCircularBackend:
    from cyclerfinder.core.constants import PLANETS

    venus = dataclasses.replace(PLANETS["V"], inc_deg=_VENUS_INC_DEG, lan_deg=_VENUS_LAN_DEG)
    return _InclinedCircularBackend(planets={"V": venus})


def _normal_angle_to_de440_deg(backend: _InclinedCircularBackend) -> float:
    """Angle (deg) between the inclined backend's orbit normal and DE440's."""
    from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, SECONDS_PER_DAY

    venus = backend._planets["V"]
    period = 2.0 * math.pi * math.sqrt((venus.sma_au * AU_KM) ** 3 / MU_SUN_KM3_S2)
    r0, v0 = backend.state("V", 0.13 * period)
    h = np.cross(r0, v0)
    n_model = np.asarray(h / np.linalg.norm(h), dtype=np.float64)
    n_de440 = _de440_orbit_normal("V", 100.0 * SECONDS_PER_DAY)
    return float(np.degrees(np.arccos(np.clip(float(np.dot(n_model, n_de440)), -1.0, 1.0))))


def test_orbit_normal_anchor_passes_clean() -> None:
    """CLEAN: fixed R_x(+inc) backend matches DE440's orbit normal to < 0.01 deg."""
    pytest.importorskip("astropy")
    angle = _normal_angle_to_de440_deg(_venus_backend())
    assert angle < 0.01


def test_orbit_normal_anchor_fires_on_rx_minus_inc(monkeypatch: pytest.MonkeyPatch) -> None:
    """POISON: flip R_x(+inc) -> R_x(-inc) -> DE440 orbit-normal anchor FAILS.

    The mirror puts the normal ~2*inc = 6.789 deg off DE440 for Venus, hundreds
    of times the 0.01 deg gate. We patch the backend's ``_rotation`` staticmethod
    (the single point that encodes the node-sign convention).
    """
    pytest.importorskip("astropy")
    monkeypatch.setattr(_InclinedCircularBackend, "_rotation", staticmethod(_rotation_rx_minus_inc))
    angle = _normal_angle_to_de440_deg(_venus_backend())
    # Gate (angle < 0.01 deg) must FAIL; the mirror lands at ~2*inc.
    assert angle > 0.01
    assert angle == pytest.approx(2.0 * _VENUS_INC_DEG, rel=0.02)


# =========================================================================== #
# POISON 3 — CROSSCHECK ENDPOINT (#197) — REFERENCED, NOT DUPLICATED
#
# The poisoned-endpoint fault test already exists and is the TEMPLATE for this
# whole harness:
#   tests/verify/test_crosscheck.py::test_crosscheck_catches_poisoned_endpoint
#   tests/verify/test_crosscheck.py::
#       test_crosscheck_escape_hatch_reproduces_shared_endpoint_blindness
# The escape-hatch test is the explicit FINDING that endpoints-from-artifact
# (independent_endpoints=False) is an undefended shared component — kept off by
# default for exactly that reason. This harness references those tests via the
# index in the module docstring rather than re-implementing them. The signature/
# ToF poison below (#180 class) exercises the SAME endpoint-independence gate
# from a different fault, confirming it is discoverable and load-bearing here.
# =========================================================================== #


# =========================================================================== #
# POISON 4 — EARTH-MOON mu (#212a)
#
# Shared component: the registry GMs (PRIMARIES["Earth"] system GM and the Moon
# GM) that cr3bp_system() derives the Earth-Moon CR3BP mass parameter from. The
# sourced-mu gate (Ross & Roberts-Tsoukkas 2025, mu = 1.2150584270572e-2) is the
# external anchor; a registry GM drift would shift mu and t_s consistently
# everywhere downstream.
#
# Gate asserted: tests/core/test_cr3bp.py::test_earth_moon_mu_physical.
# =========================================================================== #

_SOURCED_EM_MU = 1.2150584270572e-2  # Ross & Roberts-Tsoukkas 2025 (AAS 25-621) p.3
_EM_MU_REL_GATE = 1e-6  # test_cr3bp.py gate tolerance


def test_em_mu_anchor_passes_clean() -> None:
    """CLEAN: registry-derived Earth-Moon mu matches the sourced value to 1e-6."""
    sysm = cr3bp.cr3bp_system("Earth", "Moon")
    assert sysm.mu == pytest.approx(_SOURCED_EM_MU, rel=_EM_MU_REL_GATE)


def test_em_mu_anchor_fires_on_perturbed_gm(monkeypatch: pytest.MonkeyPatch) -> None:
    """POISON: perturb the Earth-system GM by 1% -> sourced-mu gate FAILS.

    mu = GM_Moon / GM_EarthSystem, so a +1% GM_EarthSystem perturbation shifts
    mu by ~-1% (~1e4 x the 1e-6 gate). This is the registry-drift class the
    sourced anchor exists to catch (cf. the original -1.2% double-count bug). We
    patch the PRIMARIES dict entry cr3bp reads through (same dict object as
    satellites.PRIMARIES, which cr3bp imports).
    """
    original = satellites.PRIMARIES["Earth"]
    monkeypatch.setitem(satellites.PRIMARIES, "Earth", original * 1.01)

    sysm = cr3bp.cr3bp_system("Earth", "Moon")
    rel_err = abs(sysm.mu - _SOURCED_EM_MU) / _SOURCED_EM_MU
    # Gate (rel < 1e-6) must FAIL; the 1% perturbation moves mu by ~1%.
    assert rel_err > _EM_MU_REL_GATE
    assert rel_err == pytest.approx(0.01, rel=0.05)


# =========================================================================== #
# POISON 5 — SIGNATURE / TRANSIT ToF (#180 class)
#
# Shared component: a leg's transit ToF, embedded in the cycler artifact as the
# pair (encounter epochs, leg t_depart/t_arrive). The #180 trap: a ToF bug that
# relabels the claimed epoch while the stored encounter geometry is stale would
# feed every solver the same wrong (r1, r2) if endpoints were read from the
# artifact. The crosscheck endpoint-INDEPENDENCE gate (#197) re-queries r1/r2
# from the ephemeris at the leg's own epochs, so a ToF poison that moves the
# claimed arrival epoch off the stale stored position is caught.
#
# Gate asserted: the crosscheck endpoint-independence check
#   (crosscheck_leg(..., independent_endpoints=True) endpoint_mismatch_km).
# =========================================================================== #


@pytest.fixture()
def constructed_cycler() -> Cycler:
    """A real-ephemeris single-leg E->M cycler (n=1, 780 d) — the crosscheck fixture."""
    pytest.importorskip("astropy")
    entry = {
        "id": "fault-injection-tof",
        "bodies": ["E", "M"],
        "legs": [{"from": "E", "to": "M", "tof_days": 780.0, "n_revs": 1}],
        "period": {"years": 2.135},
    }
    return construct_real_ephemeris_cycler(entry, Ephemeris(model="astropy"), 0.0)


def _poison_leg_tof_days(cycler: Cycler, leg_index: int, extra_days: float) -> Cycler:
    """Shift a leg's CLAIMED arrival epoch by ``extra_days`` while freezing geometry.

    Simulates a #180-class ToF bug: the arrival epoch label on both the leg and
    its matching encounter is shifted, but the stored encounter position
    (``Encounter.r``) is left at the original geometry (the stale upstream
    value). A check that read endpoints from the artifact would not notice; the
    independent ephemeris re-query at the new epoch must.
    """
    leg = cycler.legs[leg_index]
    new_t_arrive = leg.t_arrive + extra_days * 86400.0
    poisoned_leg = dataclasses.replace(leg, t_arrive=new_t_arrive)
    new_legs = [*cycler.legs]
    new_legs[leg_index] = poisoned_leg
    # Re-label the matching arrival encounter's epoch (so _leg_endpoints still
    # matches) but KEEP its stale position r.
    new_encs = []
    for enc in cycler.encounters:
        if enc.body == leg.to_body and enc.t == leg.t_arrive:
            new_encs.append(dataclasses.replace(enc, t=new_t_arrive))
        else:
            new_encs.append(enc)
    return dataclasses.replace(cycler, legs=new_legs, encounters=new_encs)


def test_tof_poison_clean_passes(constructed_cycler: Cycler) -> None:
    """CLEAN: the unpoisoned leg's endpoints match the re-query to round-off."""
    pytest.importorskip("astropy")
    eph = Ephemeris(model="astropy")
    result = crosscheck_leg(constructed_cycler.legs[0], constructed_cycler, eph, leg_index=0)
    assert result.endpoint_mismatch_km is not None
    assert result.endpoint_mismatch_km < 1.0e-6
    assert result.passed is True


def test_tof_poison_fires_endpoint_independence_gate(constructed_cycler: Cycler) -> None:
    """POISON: shift the leg's claimed ToF by 30 d (stale geometry) -> gate FAILS.

    The independent re-query reads Mars at the NEW arrival epoch; Mars moves far
    more than the 1 km consistency tolerance over 30 days, so endpoint_mismatch
    blows past POSITION_CONSISTENCY_TOL_KM and the leg fails regardless of solver
    agreement (the defence against the #180 shared-ToF false consensus).
    """
    pytest.importorskip("astropy")
    poisoned = _poison_leg_tof_days(constructed_cycler, leg_index=0, extra_days=30.0)
    eph = Ephemeris(model="astropy")
    result = crosscheck_leg(poisoned.legs[0], poisoned, eph, leg_index=0)
    assert result.endpoint_mismatch_km is not None
    assert result.endpoint_mismatch_km > POSITION_CONSISTENCY_TOL_KM
    assert result.passed is False
