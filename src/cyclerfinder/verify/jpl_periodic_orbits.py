"""JPL Three-Body Periodic Orbit Catalog client — INDEPENDENT-ORACLE source (#116).

A thin read-only client for JPL SSD/CNEOS's public "periodic_orbits" REST API
(https://ssd-api.jpl.nasa.gov/doc/periodic_orbits.html). It returns, per CR3BP
family, the published initial conditions (x, y, z, vx, vy, vz in the rotating
non-dimensional frame), Jacobi constant, period, and stability index, together
with the system constants (mass ratio, libration points, length/time units).

WHY THIS LIVES IN ``verify`` AND NOT ``data``: the JPL catalog is an
*independent reproduction by a separate group* of CR3BP periodic orbits. It is
used here as a cross-check ORACLE for our own corrector/continuation output (a
per-interface external anchor, the false-consensus discipline) and as a source
of sourced ICs for targeted families — NOT as a bulk catalogue import. Halo,
DRO, Lyapunov, etc. are not cycler trajectories; only a genuine *cycler* family
(if one ever appears in their family list) would feed the catalogue, and then as
V0-sourced seeds.

CONVENTION RECONCILIATION (mandatory before any cross-check is trusted): JPL's
Earth-Moon mass ratio is ``1.215058560962404e-2``; ours
(:func:`cyclerfinder.core.cr3bp.cr3bp_system` "Earth","Moon") is
``0.01215058439469525`` (Ross-sourced). They differ by ~1e-7 relative. An IC
taken from one mu and propagated under the other will NOT re-close exactly; the
residual is the constants mismatch, not a finding about the orbit. Always carry
both mu values explicitly and report which was used (:func:`reconcile_mu`).

NETWORK: the live :func:`query` uses the Python standard library (urllib), no
new dependency. Tests never hit the network — they parse a saved fixture via
:func:`parse_payload`.
"""

from __future__ import annotations

import json
import urllib.parse
import urllib.request
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

API_BASE = "https://ssd-api.jpl.nasa.gov/periodic_orbits.api"

#: Our Earth-Moon mass ratio (Ross-sourced; see core.cr3bp). Kept here so a
#: cross-check can compare against JPL's reported value without importing the
#: whole CR3BP stack.
OUR_EM_MU = 0.01215058439469525


@dataclass(frozen=True)
class JplSystemConstants:
    """CR3BP system constants as reported in a JPL periodic_orbits response."""

    name: str
    mu: float
    lunit_km: float
    tunit_s: float
    radius_secondary_km: float | None
    libration_points: dict[str, tuple[float, float, float]]


@dataclass(frozen=True)
class JplOrbit:
    """One published periodic orbit: rotating-frame IC + invariants."""

    state0: NDArray[np.float64]  # (x, y, z, vx, vy, vz), nondimensional
    jacobi: float
    period: float  # nondimensional
    stability: float


def parse_payload(payload: dict[str, Any]) -> tuple[JplSystemConstants, list[JplOrbit]]:
    """Parse a periodic_orbits JSON response into typed objects.

    Pure (no network) — the unit-test entry point. ``fields`` is honoured by
    name so a future column re-order cannot silently mis-map values.
    """
    sysd = payload["system"]

    def _pt(key: str) -> tuple[float, float, float] | None:
        v = sysd.get(key)
        if v is None:
            return None
        return (float(v[0]), float(v[1]), float(v[2]))

    libration = {k: pt for k in ("L1", "L2", "L3", "L4", "L5") if (pt := _pt(k)) is not None}
    rad = sysd.get("radius_secondary")
    constants = JplSystemConstants(
        name=str(sysd["name"]),
        mu=float(sysd["mass_ratio"]),
        lunit_km=float(sysd["lunit"]),
        tunit_s=float(sysd["tunit"]),
        radius_secondary_km=None if rad is None else float(rad),
        libration_points=libration,
    )

    fields: list[str] = list(payload["fields"])
    idx = {name: i for i, name in enumerate(fields)}
    needed = ("x", "y", "z", "vx", "vy", "vz", "jacobi", "period", "stability")
    missing = [n for n in needed if n not in idx]
    if missing:
        raise ValueError(f"JPL payload missing expected fields: {missing}")

    orbits: list[JplOrbit] = []
    for row in payload["data"]:
        state0 = np.array(
            [float(row[idx[c]]) for c in ("x", "y", "z", "vx", "vy", "vz")],
            dtype=np.float64,
        )
        orbits.append(
            JplOrbit(
                state0=state0,
                jacobi=float(row[idx["jacobi"]]),
                period=float(row[idx["period"]]),
                stability=float(row[idx["stability"]]),
            )
        )
    return constants, orbits


def query(
    system: str,
    family: str,
    *,
    libr: int | None = None,
    branch: str | None = None,
    timeout: float = 60.0,
) -> tuple[JplSystemConstants, list[JplOrbit]]:
    """Fetch a family live from the JPL API. Network call (urllib).

    ``system`` e.g. "earth-moon"; ``family`` e.g. "halo", "dro", "lyapunov";
    ``libr`` the libration-point index (1..5) where the family requires it;
    ``branch`` e.g. "N"/"S" where applicable.
    """
    params: dict[str, str] = {"sys": system, "family": family}
    if libr is not None:
        params["libr"] = str(libr)
    if branch is not None:
        params["branch"] = branch
    url = f"{API_BASE}?{urllib.parse.urlencode(params)}"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        payload = json.loads(resp.read().decode("utf-8"))
    return parse_payload(payload)


def reconcile_mu(jpl_mu: float, our_mu: float = OUR_EM_MU) -> dict[str, float]:
    """Report the JPL-vs-ours mass-ratio gap (carry this into any cross-check).

    Returns absolute and relative differences; a non-zero relative difference is
    the irreducible floor on how well a JPL IC can re-close under our integrator.
    """
    abs_diff = abs(jpl_mu - our_mu)
    rel_diff = abs_diff / our_mu if our_mu else float("inf")
    return {"jpl_mu": jpl_mu, "our_mu": our_mu, "abs_diff": abs_diff, "rel_diff": rel_diff}
