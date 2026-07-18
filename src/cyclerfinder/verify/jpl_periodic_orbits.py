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

CACHING (#647): :func:`query` accepts an optional ``cache_dir`` — when given,
the raw JSON payload for a given parameter set is read from / written to a
per-query file under that directory instead of re-fetching. ``cache_dir=None``
(the default) preserves the exact prior behaviour (always live network) for
every existing call site (``search/reachable_representatives.py``,
``scripts/search_campaign_daemon.py``, ``scripts/jpl_oracle_crosscheck.py``).
New callers (:mod:`cyclerfinder.search.jpl_family_check`) pass a gitignored
``out/`` cache dir per this project's existing data-caching convention (see
``.github/workflows/kernel-freshness.yml`` for the sibling NAIF-kernel
pattern) so repeated queries/tests never hammer the live API.

SUPPORTED SYSTEMS / FAMILIES (confirmed live 2026-07-18, see
:data:`SUPPORTED_SYSTEMS` / :data:`SUPPORTED_FAMILIES`): the API currently
indexes exactly 7 systems (``sun-earth, earth-moon, sun-mars, jupiter-europa,
saturn-enceladus, saturn-titan, mars-phobos`` — Sun-Jupiter is NOT covered,
correcting an earlier #641-adjacent hallucinated claim) and 12 families
(``halo, vertical, axial, lyapunov, longp, short, butterfly, dragonfly,
resonant, dro, dpo, lpo``). ``halo/lyapunov/longp/short/axial/vertical``
additionally require a ``libr`` (libration point); ``halo`` also requires a
``branch`` (``N``/``S``) per the API's own 400-error text.
"""

from __future__ import annotations

import hashlib
import json
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

API_BASE = "https://ssd-api.jpl.nasa.gov/periodic_orbits.api"

#: Our Earth-Moon mass ratio (Ross-sourced; see core.cr3bp). Kept here so a
#: cross-check can compare against JPL's reported value without importing the
#: whole CR3BP stack.
OUR_EM_MU = 0.01215058439469525

#: Systems the live API indexes (confirmed via a live 400-error parameter
#: listing, 2026-07-18) — Sun-Jupiter is deliberately NOT in this set (#647).
SUPPORTED_SYSTEMS: frozenset[str] = frozenset(
    {
        "sun-earth",
        "earth-moon",
        "sun-mars",
        "jupiter-europa",
        "saturn-enceladus",
        "saturn-titan",
        "mars-phobos",
    }
)

#: Orbit families the live API indexes (per its documentation page, confirmed
#: 2026-07-18), independent of which systems carry which families.
SUPPORTED_FAMILIES: frozenset[str] = frozenset(
    {
        "halo",
        "vertical",
        "axial",
        "lyapunov",
        "longp",
        "short",
        "butterfly",
        "dragonfly",
        "resonant",
        "dro",
        "dpo",
        "lpo",
    }
)

#: Families that require a ``libr`` (libration-point) query parameter.
FAMILIES_REQUIRING_LIBR: frozenset[str] = frozenset(
    {"halo", "lyapunov", "longp", "short", "axial", "vertical"}
)

#: Families that require a ``branch`` (e.g. "N"/"S") query parameter.
FAMILIES_REQUIRING_BRANCH: frozenset[str] = frozenset({"halo"})


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


def _cache_path(cache_dir: Path, params: dict[str, str]) -> Path:
    """Deterministic per-query-parameter-set cache filename.

    Keyed on ALL params (sorted) so two different filter windows for the same
    system/family never collide; the ``sys``/``family`` prefix keeps the
    directory human-browsable.
    """
    key = "&".join(f"{k}={params[k]}" for k in sorted(params))
    digest = hashlib.sha1(key.encode("utf-8")).hexdigest()[:16]
    sys_slug = params.get("sys", "unknown").replace("/", "_")
    fam_slug = params.get("family", "unknown").replace("/", "_")
    return cache_dir / f"{sys_slug}__{fam_slug}__{digest}.json"


def query(
    system: str,
    family: str,
    *,
    libr: int | None = None,
    branch: str | None = None,
    jacobimin: float | None = None,
    jacobimax: float | None = None,
    periodmin: float | None = None,
    periodmax: float | None = None,
    periodunits: str | None = None,
    stabmin: float | None = None,
    stabmax: float | None = None,
    timeout: float = 60.0,
    cache_dir: Path | str | None = None,
) -> tuple[JplSystemConstants, list[JplOrbit]]:
    """Fetch a family live from the JPL API (or a local cache). Network call (urllib).

    ``system`` e.g. "earth-moon"; ``family`` e.g. "halo", "dro", "lyapunov";
    ``libr`` the libration-point index (1..5) where the family requires it;
    ``branch`` e.g. "N"/"S" where applicable. ``jacobimin``/``jacobimax``/
    ``periodmin``/``periodmax``/``stabmin``/``stabmax`` are the API's own
    server-side range filters (inclusive) — use these to fetch a narrow window
    around a candidate's invariants instead of an entire family.
    ``periodunits`` is one of ``"s"``, ``"h"``, ``"d"``, ``"TU"`` (nondim, the
    API default) and only matters when a period filter is given.

    ``cache_dir``, when given, reads/writes the raw JSON payload for this
    EXACT parameter set under that directory instead of always hitting the
    network (#647) — ``None`` (the default) preserves the historical
    always-live behaviour for every pre-existing call site.
    """
    params: dict[str, str] = {"sys": system, "family": family}
    if libr is not None:
        params["libr"] = str(libr)
    if branch is not None:
        params["branch"] = branch
    if jacobimin is not None:
        params["jacobimin"] = repr(float(jacobimin))
    if jacobimax is not None:
        params["jacobimax"] = repr(float(jacobimax))
    if periodmin is not None:
        params["periodmin"] = repr(float(periodmin))
    if periodmax is not None:
        params["periodmax"] = repr(float(periodmax))
    if periodunits is not None:
        params["periodunits"] = periodunits
    if stabmin is not None:
        params["stabmin"] = repr(float(stabmin))
    if stabmax is not None:
        params["stabmax"] = repr(float(stabmax))

    cache_path: Path | None = None
    if cache_dir is not None:
        cache_path = _cache_path(Path(cache_dir), params)
        if cache_path.exists():
            payload = json.loads(cache_path.read_text())
            return parse_payload(payload)

    url = f"{API_BASE}?{urllib.parse.urlencode(params)}"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        payload = json.loads(resp.read().decode("utf-8"))

    if cache_path is not None:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(payload))

    return parse_payload(payload)


def reconcile_mu(jpl_mu: float, our_mu: float = OUR_EM_MU) -> dict[str, float]:
    """Report the JPL-vs-ours mass-ratio gap (carry this into any cross-check).

    Returns absolute and relative differences; a non-zero relative difference is
    the irreducible floor on how well a JPL IC can re-close under our integrator.
    """
    abs_diff = abs(jpl_mu - our_mu)
    rel_diff = abs_diff / our_mu if our_mu else float("inf")
    return {"jpl_mu": jpl_mu, "our_mu": our_mu, "abs_diff": abs_diff, "rel_diff": rel_diff}
