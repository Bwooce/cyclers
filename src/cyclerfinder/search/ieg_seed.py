"""EGGIE Lambert-real seed adapter for the #480 M1 Galilean cycler reproduction.

Builds a :class:`~cyclerfinder.nbody.shooter.ShootingSeed` for the
Hernandez 2017 (AAS 17-608) Io-Europa-Ganymede **EGGIE** cycler using real
Galilean-moon positions from the jup365.bsp SPICE kernel and Lambert-solve
spacecraft velocities for each leg.  This seed is the warm-start for the Task 4
multiple-shooting corrector.

**Multi-revolution topology (Task 4b)**
The paper's EGGIE cycler is a 4-synodic-period / five-spacecraft-revolution tour
(Hernandez 2017, Table 1 note: 4 synodic periods → sc:syn = 4:5, i.e. 5 spacecraft
revolutions about Jupiter).  Single-revolution Lambert legs miss the correct
trajectory topology, giving V-inf that are 2-3x off the paper and landing in the wrong
basin.

Each leg is therefore built with the **multi-revolution Lambert branch** whose
spacecraft departure V-inf magnitude is CLOSEST to the sourced Table-4 value at that
encounter (selection principle: iterate N=0..3, pick the N that minimises
|v_sc_departure - moon_velocity| - paper_vinf_target).

**Chosen per-leg revolution counts (Task 4b, sourced from Table 4 V∞)**

Leg 0 - Europa->Ganymede, ToF 1.59 d:
  N=0 (single-rev).  No multi-rev solution is feasible for the 1.59-day
  time-of-flight in this geometry; the short arc is geometry-constrained to
  single-revolution.  Characterized: node-0 and node-1 V∞ remain ~3.7 and ~4.7 km/s
  from the Table-4 targets at departure_et=0 (epoch-geometry gap, not a solver gap).

Leg 1 - Ganymede->Ganymede, ToF 8.60 d:
  N=1, HIGH branch.  Vinf ~7.04/7.04 km/s vs sourced 7.07/7.07 km/s (avg_err ~0.03).
  The Ganymede-Ganymede equal-V-inf property (adjacent same-body encounters in a
  ballistic cycler always have equal V-inf, Hernandez sec 3) is recovered precisely.

Leg 2 - Ganymede->Io, ToF 7.34 d:
  N=2, HIGH branch.  Best departure V∞ match to the 7.07 km/s target (7.96 vs 7.07,
  avg_err ≈ 6.29 — the Io arrival is still ~20 km/s vs 8.38 km/s target).  This
  epoch/geometry cannot produce a Lambert arc matching the Io arrival; the paper's
  multi-body patched-conic arc through the actual encounter geometry is not
  reproducible with the Keplerian Lambert legs at departure_et=0.

Leg 3 - Io->Europa, ToF 10.69 d:
  N=1, HIGH branch.  Vinf ~7.52/8.73 km/s vs sourced 8.38/9.12 (avg_err ~0.62).

Total spacecraft revolutions: 0+1+2+1 = **4** (paper states ~5 for the 4-synodic
EGGIE; the departure_et=0 geometry gives 4, which is a better-than-single-rev
approximation of the topology — the paper epoch ET ≈ 6.55e8 s may differ).

**Seed defect characterisation (Task 4b)**
At departure_et=0 the multi-rev seed defect_norm is ~1.58e6 km/(km/s),
reduced from ~2.47e6 (single-rev) but NOT below 1e6.  The dominant residual
is Leg 2 Io-arrival (~20 km/s vs 8.38 km/s target) — a geometry gap that
multi-rev Keplerian arcs alone cannot close at this epoch.

**Departure epoch note**
``departure_et=0.0`` (J2000 TDB, 2000-01-01 12:00 TDB) is a placeholder anchor.
The paper's Table 4 departure date is 29-Sep-2020 ET (approximately
ET = 6.55e8 s past J2000).  Task 4 sweeps the epoch to find the minimum-residual
departure; the default 0.0 is geometric only.

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

# Per-leg Lambert revolution counts and branches (Task 4b multi-rev topology).
# Selected by V∞-match principle: for each leg the (n_revs, branch) pair whose
# spacecraft departure V∞ magnitude best matches the paper's sourced Table-4
# value at that node.  See module docstring for per-leg rationale.
# Sourced V∞ targets: Europa=9.12, Ganymede=7.07, Io=8.38 km/s (Table 4).
#
#  Leg 0  Europa→Ganymede   1.59 d  n=0  single  (only feasible rev for 1.59-day arc)
#  Leg 1  Ganymede→Ganymede 8.60 d  n=1  high    (7.04/7.04 km/s vs sourced 7.07/7.07)
#  Leg 2  Ganymede→Io       7.34 d  n=2  high    (7.96 km/s dep vs sourced 7.07; best available)
#  Leg 3  Io→Europa         10.69 d n=1  high    (7.52/8.73 km/s vs sourced 8.38/9.12)
#
# Total: 0+1+2+1 = 4 spacecraft revolutions (paper states ~5 for 4-synodic EGGIE;
# departure_et=0 geometry gives 4).
EGGIE_LEG_NREVS: list[int] = [0, 1, 2, 1]
EGGIE_LEG_BRANCHES: list[str] = ["single", "high", "high", "high"]

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
    """Build the EGGIE multi-revolution Lambert-real seed at ``departure_et``.

    Each leg is built with the multi-revolution Lambert branch selected by the
    V∞-match principle (see :data:`EGGIE_LEG_NREVS` / :data:`EGGIE_LEG_BRANCHES`):
    for each leg the (n_revs, branch) pair whose spacecraft departure V∞ magnitude
    best matches the paper's sourced Table-4 value at that encounter.  The chosen
    per-leg revolution counts are [0, 1, 2, 1] (total 4 revolutions).

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
    # Multi-rev branches are selected by the V∞-match principle:
    #   for each leg use the (n_revs, branch) from EGGIE_LEG_NREVS / EGGIE_LEG_BRANCHES
    #   whose spacecraft departure V∞ best matches the paper's sourced Table-4 value.
    # lambert(r1, r2, tof_sec, mu=..., max_revs=n) -> list[LambertSolution].
    # sc_v_dep[i] = spacecraft velocity at START of leg i (outgoing at node i).
    # sc_v_arr[i] = spacecraft velocity at END   of leg i (arriving at node i+1).
    n_legs = len(EGGIE_TOFS_DAYS)
    sc_v_dep: list[NDArray[np.float64]] = []
    sc_v_arr: list[NDArray[np.float64]] = []

    for i in range(n_legs):
        r1 = moon_r[i]
        r2 = moon_r[i + 1]
        tof_sec = EGGIE_TOFS_DAYS[i] * _SPD
        n_revs = EGGIE_LEG_NREVS[i]
        branch = EGGIE_LEG_BRANCHES[i]
        sols = lambert(r1, r2, tof_sec, mu=_MU_JUPITER, max_revs=n_revs)
        # Select the solution matching the chosen (n_revs, branch).
        sol = next((s for s in sols if s.n_revs == n_revs and s.branch == branch), None)
        if sol is None:
            # Fallback: if the multi-rev branch is infeasible at this epoch/geometry,
            # use the single-rev solution.  This can happen for short-ToF legs when
            # the requested revolution is below t_min; the single-rev is always returned.
            sol = sols[0]
        sc_v_dep.append(np.asarray(sol.v1, dtype=np.float64))
        sc_v_arr.append(np.asarray(sol.v2, dtype=np.float64))

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
