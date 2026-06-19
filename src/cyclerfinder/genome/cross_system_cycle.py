"""#405 cross-system (Sun-Earth <-> Earth-Moon) heteroclinic-cycle framework.

Task 1: the inter-system frame bridge. Maps a 6-state between the SE-rotating CR3BP
frame, a common Earth-centered inertial frame (km, km/s), and the EM-rotating CR3BP
frame, parameterized by ``theta`` (the relative phase between the Sun-Earth line and
the Earth-Moon line). Correctness is gated by a round-trip identity AND a physical
Moon-position test (see tests/genome/test_cross_system_cycle.py).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS


def se_earth_system() -> cr3bp.CR3BPSystem:
    """Build the Sun-Earth CR3BP system (Earth is the SE secondary at 1-mu_SE).

    The shared registry (``core.satellites``) carries Earth only as a *primary* and
    has no Sun entry, so ``cr3bp.cr3bp_system("Sun", "Earth")`` cannot construct the
    SE pair. We assemble it here from the same sourced constants used everywhere else:
    ``MU_SUN_KM3_S2`` (IAU 2015 nominal solar GM), Earth's planet-only GM, and Earth's
    heliocentric SMA. Convention matches ``cr3bp.cr3bp_system``: pair GM = G(m1+m2),
    ``l_km`` = secondary SMA about primary, ``t_s`` = sqrt(l_km^3 / G(m1+m2)).
    """
    earth = PLANETS["E"]
    gm_pair = MU_SUN_KM3_S2 + earth.mu_km3_s2
    mu = earth.mu_km3_s2 / gm_pair
    l_km = earth.sma_au * AU_KM
    t_s = math.sqrt(l_km**3 / gm_pair)
    return cr3bp.CR3BPSystem(mu=mu, primary="Sun", secondary="Earth", l_km=l_km, t_s=t_s)


def em_moon_system() -> cr3bp.CR3BPSystem:
    """Build the Earth-Moon CR3BP system via the shared registry (Moon is registered)."""
    return cr3bp.cr3bp_system("Earth", "Moon")


def _rot_z(angle: float) -> NDArray[np.float64]:
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], dtype=np.float64)


@dataclass(frozen=True)
class FrameBridge:
    """Transforms 6-states between SE-rot, Earth-centered inertial (km, km/s), EM-rot.

    SE-rot and EM-rot are nondimensional CR3BP rotating frames; the inertial frame is
    Earth-centered, dimensional. ``theta`` is the inertial angle of the EM-rot x-axis
    minus that of the SE-rot x-axis (the SE-rot x-axis is taken at inertial angle 0).
    """

    se: cr3bp.CR3BPSystem
    em: cr3bp.CR3BPSystem

    def se_rot_to_inertial(self, s: NDArray[np.float64], *, theta: float) -> NDArray[np.float64]:
        r = np.asarray(s[:3], float).copy()
        v = np.asarray(s[3:], float).copy()
        r = r - np.array([1.0 - self.se.mu, 0.0, 0.0])  # Earth-centered, SE-rot, nondim
        lkm, vunit = self.se.l_km, self.se.l_km / self.se.t_s
        pos_km = r * lkm
        vel_rot = v * vunit
        omega = np.array([0.0, 0.0, 1.0 / self.se.t_s])
        vel_in = vel_rot + np.cross(omega, pos_km)
        # SE-rot x-axis is taken at inertial angle 0, so no theta rotation here.
        return np.concatenate([pos_km, vel_in])

    def inertial_to_se_rot(self, x: NDArray[np.float64], *, theta: float) -> NDArray[np.float64]:
        pos_km = np.asarray(x[:3], float).copy()
        vel_in = np.asarray(x[3:], float).copy()
        omega = np.array([0.0, 0.0, 1.0 / self.se.t_s])
        vel_rot = vel_in - np.cross(omega, pos_km)
        lkm, vunit = self.se.l_km, self.se.l_km / self.se.t_s
        r = pos_km / lkm + np.array([1.0 - self.se.mu, 0.0, 0.0])
        v = vel_rot / vunit
        return np.concatenate([r, v])

    def em_rot_to_inertial(self, s: NDArray[np.float64], *, theta: float) -> NDArray[np.float64]:
        r = np.asarray(s[:3], float).copy()
        v = np.asarray(s[3:], float).copy()
        r = r - np.array([-self.em.mu, 0.0, 0.0])  # Earth-centered, EM-rot, nondim
        lkm, vunit = self.em.l_km, self.em.l_km / self.em.t_s
        pos_emrot = r * lkm
        vel_rot = v * vunit
        omega = np.array([0.0, 0.0, 1.0 / self.em.t_s])
        vel_in_emrot = vel_rot + np.cross(omega, pos_emrot)
        rot = _rot_z(theta)  # EM-rot x-axis leads SE-rot by theta in inertial frame
        return np.concatenate([rot @ pos_emrot, rot @ vel_in_emrot])

    def inertial_to_em_rot(self, x: NDArray[np.float64], *, theta: float) -> NDArray[np.float64]:
        rot_inv = _rot_z(-theta)
        pos_emrot = rot_inv @ np.asarray(x[:3], float)
        vel_in_emrot = rot_inv @ np.asarray(x[3:], float)
        omega = np.array([0.0, 0.0, 1.0 / self.em.t_s])
        vel_rot = vel_in_emrot - np.cross(omega, pos_emrot)
        lkm, vunit = self.em.l_km, self.em.l_km / self.em.t_s
        r = pos_emrot / lkm + np.array([-self.em.mu, 0.0, 0.0])
        v = vel_rot / vunit
        return np.concatenate([r, v])
