"""EGGIE Lambert-real seed adapter for the #480 M1 Galilean cycler reproduction.

Builds a :class:`~cyclerfinder.nbody.shooter.ShootingSeed` for the
Hernandez 2017 (AAS 17-608) Io-Europa-Ganymede **EGGIE** cycler using real
Galilean-moon positions from the jup365.bsp SPICE kernel and Lambert-solve
spacecraft velocities for each leg.  This seed is the warm-start for the Task 4
multiple-shooting corrector.

**Seed nature (not a sourced golden)**
The node positions are the real moon positions at the sourced ToF epochs anchored
at ``departure_et``.  The spacecraft velocities come from a single-revolution
Lambert solve Jupiter-centred; the resulting V-inf magnitudes are a GEOMETRIC
ESTIMATE in real ephemeris, NOT the sourced 9.12 / 7.07 / 7.07 / 8.38 km/s from
Table 4 of the digest.  Matching the paper's V-inf is Task 4/6's job; this seed
only needs to be in the right basin.

**Departure epoch note**
``departure_et=0.0`` (J2000 TDB, 2000-01-01 12:00 TDB) is a placeholder anchor.
The paper's Table 4 departure date is 29-Sep-2020 ET (approximately
ET = 656 days x 86400 ~ 6.58e7 s past J2000).  Task 4 sweeps the epoch
to find the minimum-residual departure; the default 0.0 is geometric only.

**Sources**
- Hernandez, Jones & Jesick (2017), "One Class of Io-Europa-Ganymede Triple
  Cyclers," AAS 17-608, Table 4.
  Digest: docs/notes/2026-06-26-digest-hernandez-2017-ieg-triple-cyclers-aas-17-608.md
- Jupiter GM: PRIMARIES["Jupiter"] from core/satellites.py (JPL DE440,
  1.26686534e8 km^3/s^2).
- Moon positions: jup365.bsp via Ephemeris(center="Jupiter", model="spice").
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.lambert import lambert
from cyclerfinder.core.satellites import PRIMARIES
from cyclerfinder.nbody.shooter import ShootingSeed

# ---------------------------------------------------------------------------
# Sourced EGGIE invariants - Table 4, Hernandez 2017 (AAS 17-608).
# Source: docs/notes/2026-06-26-digest-hernandez-2017-ieg-triple-cyclers-aas-17-608.md
# ---------------------------------------------------------------------------

# Encounter sequence (Europa first, by the paper's convention).
EGGIE_SEQUENCE: tuple[str, ...] = ("Europa", "Ganymede", "Ganymede", "Io", "Europa")

# Inter-encounter times of flight [days], in order (4 legs for 5 encounters).
# Sourced from Table 4: 1.59, 8.60, 7.34, 10.69 d.
EGGIE_TOFS_DAYS: list[float] = [1.59, 8.60, 7.34, 10.69]

# Total period (4 synodic periods * T_syn = 7.05 d).
# Sourced: T_syn = 7.05 d, total = 28.22 d (Table 4 sum).
EGGIE_PERIOD_DAYS: float = 28.22

# Jupiter gravitational parameter [km^3/s^2] - from PRIMARIES registry (JPL DE440).
_MU_JUPITER: float = PRIMARIES["Jupiter"]

# Seconds per day.
_SPD: float = 86400.0

# Paper EGGIE departure window (Table 4 caption: 29-Sep-2020; page-9 text + Fig 4:
# 02-Oct-2020 — the digest flags the paper's own date inconsistency). We anchor on
# 02-Oct-2020 12:00 TDB; the resulting ET is ~6.55e8 s past J2000 (a prior note had
# this 10x wrong at ~6.6e7 — the sanity guard in PAPER_DEPARTURE_ET_APPROX catches it).
PAPER_DEPARTURE_UTC: str = "2020-OCT-02 12:00 TDB"
PAPER_DEPARTURE_ET_APPROX: float = 6.55e8


def paper_departure_et() -> float:
    """SPICE ET (TDB seconds past J2000) for the paper's EGGIE departure epoch.

    Resolves :data:`PAPER_DEPARTURE_UTC` via ``spiceypy.str2et`` after furnishing the
    NAIF leapseconds kernel. Returns ~6.55e8 s (NOT ~6.6e7 — the historical 10x error).
    """
    import spiceypy

    from cyclerfinder.verify.spice_kernels import ensure_leapseconds_kernel

    spiceypy.furnsh(ensure_leapseconds_kernel())
    return float(spiceypy.str2et(PAPER_DEPARTURE_UTC))


def ieg_eggie_seed(departure_et: float = 0.0) -> ShootingSeed:
    """Build the EGGIE Lambert-real seed at ``departure_et`` (TDB seconds past J2000).

    Parameters
    ----------
    departure_et:
        Departure epoch as TDB seconds past J2000 (= SPICE ET).  Default 0.0 is
        a placeholder geometric anchor; use Task 4's epoch sweep to find the
        departure nearest the paper's 29-Sep-2020 / 02-Oct-2020 date window.

    Returns
    -------
    ShootingSeed
        Node states = [moon_position | lambert_spacecraft_velocity] (km, km/s),
        Jupiter-centred J2000-equatorial frame.  All 5 encounters are populated.
        ``vinf_in``/``vinf_out`` are spacecraft_velocity - moon_velocity at each
        node for the inbound / outbound legs respectively.
    """
    # Real Galilean moon ephemeris via jup365.bsp.
    ephem = Ephemeris(center="Jupiter", model="spice")

    # --- Build encounter epochs (departure + cumulative ToF sums). ---
    epochs: list[float] = [departure_et]
    for tof_d in EGGIE_TOFS_DAYS:
        epochs.append(epochs[-1] + tof_d * _SPD)

    # --- Gather moon positions at each encounter epoch. ---
    n_enc = len(EGGIE_SEQUENCE)
    moon_r: list[NDArray[np.float64]] = []
    moon_v: list[NDArray[np.float64]] = []
    for i in range(n_enc):
        r, v = ephem.state(EGGIE_SEQUENCE[i], epochs[i])
        moon_r.append(np.asarray(r, dtype=np.float64))
        moon_v.append(np.asarray(v, dtype=np.float64))

    # --- Lambert-solve each leg to get spacecraft velocities at each node. ---
    # lambert(r1, r2, tof_sec, mu=...) -> list[LambertSolution]; use [0] (single rev).
    # sc_v_dep[i] = spacecraft velocity at START of leg i (outgoing at node i).
    # sc_v_arr[i] = spacecraft velocity at END   of leg i (arriving at node i+1).
    n_legs = len(EGGIE_TOFS_DAYS)
    sc_v_dep: list[NDArray[np.float64]] = []
    sc_v_arr: list[NDArray[np.float64]] = []

    for i in range(n_legs):
        r1 = moon_r[i]
        r2 = moon_r[i + 1]
        tof_sec = EGGIE_TOFS_DAYS[i] * _SPD
        sols = lambert(r1, r2, tof_sec, mu=_MU_JUPITER)
        # sols[0] is always the single-rev solution (n_revs=0, branch="single").
        sc_v_dep.append(np.asarray(sols[0].v1, dtype=np.float64))
        sc_v_arr.append(np.asarray(sols[0].v2, dtype=np.float64))

    # --- Assemble per-node spacecraft velocities. ---
    # Convention (mirrors seed_from_conic in shooter.py):
    #   node_i state = [moon_r[i] | sc_v_at_node_i]
    #   For interior nodes the outgoing velocity (start of next leg) is used.
    #   For the terminal node (i = n_enc-1) the arriving velocity is used.
    #
    # vinf_in[i]  = sc arriving at node i - moon_v[i]
    #             = sc_v_arr[i-1] - moon_v[i]  for i >= 1; zero vector for i=0.
    # vinf_out[i] = sc departing node i - moon_v[i]
    #             = sc_v_dep[i] - moon_v[i]  for i < n_legs; zero for terminal.

    node_states: list[NDArray[np.float64]] = []
    vinf_in: list[NDArray[np.float64]] = []
    vinf_out: list[NDArray[np.float64]] = []

    for i in range(n_enc):
        r_enc = moon_r[i]
        # Spacecraft velocity at this node (outgoing, or arriving for the terminal).
        v_sc = sc_v_dep[i] if i < n_legs else sc_v_arr[n_legs - 1]

        node_states.append(np.concatenate([r_enc, v_sc]))

        # V-inf inbound (what the spacecraft had when it arrived here).
        vin = np.zeros(3, dtype=np.float64) if i == 0 else sc_v_arr[i - 1] - moon_v[i]
        vinf_in.append(vin)

        # V-inf outbound (what the spacecraft has when it leaves here).
        vout = sc_v_dep[i] - moon_v[i] if i < n_legs else sc_v_arr[n_legs - 1] - moon_v[i]
        vinf_out.append(vout)

    # slack_leg = index of the longest ToF (the period pin, corrector convention).
    slack_leg = int(np.argmax(EGGIE_TOFS_DAYS))

    return ShootingSeed(
        node_states=node_states,
        epochs=epochs,
        tofs=list(EGGIE_TOFS_DAYS),
        sequence=EGGIE_SEQUENCE,
        slack_leg=slack_leg,
        period_days=EGGIE_PERIOD_DAYS,
        vinf_in=vinf_in,
        vinf_out=vinf_out,
    )
