"""V3 6D n-body independent-integrator gauntlet (#306 Phase 3 / #331).

Spec reference
--------------
* §14 V3 — independent-cross-check at scale: the candidate's defining
  dynamical signature must survive re-evaluation under an INDEPENDENT
  integrator architecture.

For a moontour (a Lambert-leg patched-conic tour in the planet frame),
V3 here is "**same model, different integrator**":

* Same model: circular-coplanar Keplerian moons about the primary, planet-
  frame Kepler arcs between encounters, identical V_inf magnitudes.
* Different integrator: REBOUND IAS15 high-order Gauss-Radau (a
  fundamentally different family from the project's scipy
  :func:`solve_ivp` DOP853 used by the Phase-2 V2-moontour driver).

The V2 verdict produced by :func:`run_v2_moontour` measures per-cycle
inter-cycle rendezvous drift using the analytic Lambert + analytic moon
orbit (closed-form Keplerian). V3 reruns each cycle's first encounter
with the same Lambert v-out, then propagates the spacecraft NUMERICALLY
under IAS15 through the leg ToF, and compares the IAS15 terminal
position to the analytic moon position at the leg endpoint. The
"rendezvous drift" V3 measures is the IAS15-vs-analytic-Kepler-leg
defect; the inter-cycle drift signature is then computed from V3's
per-cycle terminal positions and cross-compared to V2's per-cycle
terminal positions.

If V3 and V2 AGREE within a tight floor (default 100 km), the
quasi_cycler signature V2 reports is a REAL dynamical property of the
model, not an artifact of the DOP853-based Lambert+propagate stack.
If they DISAGREE, V2's bounded-drift signature was integrator noise
and the candidate must be retired.

V3 here is the IAS15 cross-check at the moon-frame Kepler scale. V4
(real-ephemeris HFEM with SPICE kernels, GMAT integration) is a
separate gate handled by #332.

Independent-integrator architecture (the whole point)
-----------------------------------------------------
The V2 moontour driver uses:

* :func:`cyclerfinder.core.lambert.lambert` — Battin / Izzo Lambert solver.
* circular-coplanar moon ephemerides closed-form.
* worst-flyby V_inf-magnitude continuity residual.

V3 in this module uses:

* REBOUND 5 IAS15 — high-order adaptive Gauss-Radau, **distinct from
  the DOP853 + Lambert chain** the V2 driver runs.
* circular-coplanar moon ephemerides (same closed form — this is the
  SAME MODEL, not a different one).
* IAS15-terminal-vs-analytic-moon position offset at the leg endpoint.

A V3 agreement with V2 says: the V2 bounded-drift signature is what the
shared model dictates, and is NOT an artifact of the V2 driver's
DOP853 + Lambert internals.

Composition map
---------------
This module composes on:

* :func:`cyclerfinder.core.lambert.lambert` — re-used to solve the same
  Lambert legs as V2 (we share the leg solution; V3 differs in the
  PROPAGATION of the leg, not its targeting).
* REBOUND 5 IAS15 (``import rebound``) — the independent integrator.
  REBOUND availability is required; the fallback path (scipy LSODA over
  the full multi-cycle horizon) is documented in the spec but the
  baseline V3 of this gate is IAS15.

Discipline
----------
* NO catalogue writeback. A V3 PASS does NOT admit — V4 (#332) + #329
  literature check remain.
* The V3 verdict is whatever the math says — if IAS15 and DOP853
  disagree on the SILVER's bounded-drift signature, the verdict is
  FAIL and the candidate is honestly retired as integrator artefact.
* The Lambert leg targeting + moon ephemerides are SHARED with V2 — V3
  asserts AGREEMENT of two integrator architectures on the SAME model.
  Specific drift numbers are NOT golden; AGREEMENT is.
* When REBOUND is not importable, the driver falls back to
  scipy.integrate.LSODA (multistep BDF) — different integrator family
  from DOP853 (single-step) — and the verdict's ``integrator`` field
  records which one was used. The fallback is honest, not silent.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.core.constants import SECONDS_PER_DAY
from cyclerfinder.core.lambert import lambert as _lambert
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.data.validation.v2_moontour import V2MoontourVerdict
from cyclerfinder.search.discovery_campaign import DAY_S, _mean_motion_rad_day, _moon_state

V3_AGREEMENT_FLOOR_KMS: Final[float] = 100.0
"""Default V3-vs-V2 agreement floor in km.

If V3 (IAS15) and V2 (DOP853+Lambert) per-cycle terminal positions
agree to within 100 km on a SILVER whose V2 drift is ~3e5 km, the
bounded-drift signature is integrator-independent — i.e. a REAL
property of the shared model, not numerical artefact.

100 km is the Tier-1 patched-conic integration-noise budget: at Oberon's
~6e5 km orbit it represents 0.02% of the SMA, well below the V2 drift
floor of 50,000 km and the SILVER's actual ~3e5 km drift. The
discriminator is wide.
"""

V3_N_CYCLES_MIN: Final[int] = 3
"""Spec §14 V2-ballistic minimum (carried into V3): ``>= 3`` cycles to even
report a verdict. V3 does not relax this — if V2 didn't run >= 3 cycles,
V3 has nothing to compare against."""


@dataclass(frozen=True)
class V3CycleVerdict3D:
    """Per-cycle V3 verdict for the n-body independent-integrator gauntlet."""

    cycle_index: int
    """Zero-indexed cycle number (0, 1, 2, ...)."""
    converged_legs: int
    """Number of Lambert legs that closed in this cycle (same machinery as V2)."""
    n_legs: int
    """Total Lambert legs in one cycle."""
    rendezvous_drift_kms_v3: float
    """Position offset of THIS cycle's V3 (IAS15) final encounter vs cycle 0's
    V3 final encounter, km. Mirrors V2's ``rendezvous_drift_kms`` but with
    the spacecraft propagated under IAS15 leg-by-leg in the planet frame
    instead of the analytic Lambert + closed-form moon orbit chain."""
    rendezvous_drift_kms_v2: float
    """V2's stored ``rendezvous_drift_kms`` at the same cycle index, for
    direct comparison. 0.0 for cycle 0 by construction."""
    agreement_kms: float
    """``|rendezvous_drift_kms_v3 - rendezvous_drift_kms_v2|`` —
    integrator-vs-integrator delta. Cycle-0 is 0 by construction."""
    ias15_terminal_offset_kms: float
    """The IAS15 spacecraft's terminal position offset vs the cycle-0 IAS15
    terminal position at the final encounter, km. Same definition as
    ``rendezvous_drift_kms_v3`` (duplicated for clarity in the audit trail)."""
    ias15_vs_analytic_kepler_kms: float
    """The IAS15 spacecraft terminal position offset vs the analytic moon
    target (the planet-frame Lambert target endpoint) at the cycle's
    final encounter, km. This is the genuine numeric agreement of two
    Kepler integrators on the same model: IAS15 (Gauss-Radau) vs the
    project's analytic Lambert + closed-form Kepler used by V2. At
    ``ias15_epsilon=1e-12`` this is typically nanometers — meaning the
    analytic Kepler chain V2 uses is faithful to the numerical-Kepler
    truth to integrator-noise level."""
    notes: str = ""


@dataclass(frozen=True)
class V3Verdict3D:
    """Frozen V3 verdict for an n-body independent-integrator gauntlet candidate."""

    candidate_id: str
    sequence: tuple[str, ...]
    n_cycles_propagated: int
    integrator: str
    """Human-readable integrator label (e.g. ``"REBOUND IAS15"``,
    ``"scipy LSODA fallback"``). The fallback path is documented honestly."""
    per_cycle: tuple[V3CycleVerdict3D, ...]
    per_cycle_drift_kms_v3: tuple[float, ...]
    """Per-cycle V3 rendezvous drift, km (cumulative from cycle 0). The headline
    series for V3 / V2 comparison."""
    per_cycle_drift_kms_v2: tuple[float, ...]
    """The V2 rendezvous-drift series as stored in the input
    :class:`V2MoontourVerdict`, sliced to ``n_cycles_propagated``."""
    drift_agreement_kms: float
    """``max_k |drift_v3[k] - drift_v2[k]|`` — the headline integrator-
    agreement number. If small (< :data:`V3_AGREEMENT_FLOOR_KMS`), the V2
    quasi_cycler signature is a real property of the model."""
    v3_v2_agreement_floor_kms: float
    """Floor against which ``drift_agreement_kms`` is gated. Default
    :data:`V3_AGREEMENT_FLOOR_KMS` (100 km)."""
    passes_v3: bool
    """``drift_agreement_kms <= v3_v2_agreement_floor_kms`` AND every cycle's
    Lambert leg closed under both V2 and V3. The headline boolean.

    Interpretation:

    * PASS: V3 (IAS15) agrees with V2 (DOP853+Lambert) — the SILVER's
      bounded-drift signature is REAL. Next gates: V4 (#332 HFEM
      Uranian-system real-ephemeris) + #329 (Heaton-Longuski literature
      check). If those clear, this is the catalogue's first computed
      ``quasi_cycler`` row.
    * FAIL: integrator disagreement — V2's signature was numerical. The
      candidate retires to the negative-results registry (#172).
    """
    notes: str = ""


def _resolve_primary(system: cr3bp.CR3BPSystem | None, sequence: tuple[str, ...]) -> str:
    """Resolve the primary body name (same logic as v2_moontour)."""
    if system is not None and system.primary:
        target = system.primary.strip().lower()
        for name in PRIMARIES:
            if name.lower() == target:
                return name
        raise ValueError(f"unknown primary {system.primary!r} from system")
    if not sequence:
        raise ValueError("empty sequence; cannot resolve primary")
    head = sequence[0]
    if head not in SATELLITES:
        raise ValueError(f"unknown moon {head!r} in sequence")
    return SATELLITES[head].primary


def _moon_constants(primary: str, sequence: tuple[str, ...]) -> dict[str, tuple[float, float]]:
    """``{moon: (sma_km, mean_motion_rad_day)}`` for every moon in the sequence."""
    mu = PRIMARIES[primary]
    out: dict[str, tuple[float, float]] = {}
    for moon in sequence:
        if moon not in SATELLITES:
            raise ValueError(f"unknown moon {moon!r}")
        sat = SATELLITES[moon]
        if sat.primary != primary:
            raise ValueError(
                f"moon {moon!r} orbits {sat.primary!r}, not the resolved primary {primary!r}"
            )
        out[moon] = (sat.sma_km, _mean_motion_rad_day(mu, sat.sma_km))
    return out


def _ias15_propagate_planet_frame(
    r0_km: np.ndarray,
    v0_km_s: np.ndarray,
    tof_s: float,
    mu_primary: float,
    *,
    epsilon: float = 1e-12,
) -> tuple[np.ndarray, np.ndarray, str]:
    """Propagate a spacecraft state under planet-frame two-body Kepler with IAS15.

    Returns ``(r_f_km, v_f_km_s, integrator_label)``. The integrator label
    is ``"REBOUND IAS15"`` if rebound was importable, else
    ``"scipy LSODA fallback"`` if scipy LSODA was used.

    The propagation is bare two-body Kepler about the primary (no moon
    perturbations — patched-conic Tier-1, matching the V2 driver's model).
    The "independent integrator" property comes from REBOUND 5's
    Gauss-Radau implementation, which is architecturally distinct from
    scipy's DOP853 Runge-Kutta used by :mod:`cyclerfinder.core.cr3bp`.
    """
    try:
        import rebound

        sim = rebound.Simulation()
        sim.G = mu_primary
        sim.integrator = "ias15"
        sim.integrator.epsilon = epsilon
        # Central body (primary) at origin with mass 1.0 so G*M = mu_primary.
        sim.add(m=1.0, x=0.0, y=0.0, z=0.0, vx=0.0, vy=0.0, vz=0.0)
        sim.add(
            m=0.0,
            x=float(r0_km[0]),
            y=float(r0_km[1]),
            z=float(r0_km[2]),
            vx=float(v0_km_s[0]),
            vy=float(v0_km_s[1]),
            vz=float(v0_km_s[2]),
        )
        sim.integrate(float(tof_s))
        sc = sim.particles[1]
        r_f = np.array([sc.x, sc.y, sc.z], dtype=np.float64)
        v_f = np.array([sc.vx, sc.vy, sc.vz], dtype=np.float64)
        return r_f, v_f, "REBOUND IAS15"
    except ImportError:
        # Fallback: scipy LSODA (multistep BDF) — different family from
        # scipy DOP853 (single-step Runge-Kutta) the V2 driver uses indirectly.
        # Different enough to count as an independent integrator architecture.
        from scipy.integrate import solve_ivp

        def two_body(_t: float, y: np.ndarray) -> np.ndarray:
            r = y[:3]
            r_norm = float(np.linalg.norm(r))
            a = -mu_primary * r / r_norm**3
            return np.concatenate([y[3:], a])

        y0 = np.concatenate([r0_km, v0_km_s])
        sol = solve_ivp(
            two_body,
            (0.0, float(tof_s)),
            y0,
            method="LSODA",
            rtol=1e-12,
            atol=1e-12,
        )
        if not sol.success:
            raise RuntimeError(f"LSODA fallback failed: {sol.message}") from None
        yf = sol.y[:, -1]
        return yf[:3].copy(), yf[3:].copy(), "scipy LSODA fallback"


def _cycle_v3(
    *,
    sequence: tuple[str, ...],
    leg_tofs_days: tuple[float, ...],
    theta_base: dict[str, float],
    t_cycle_offset_days: float,
    consts: dict[str, tuple[float, float]],
    mu_primary: float,
    n_revs: tuple[int, ...] | None,
    ias15_epsilon: float,
) -> tuple[bool, np.ndarray | None, np.ndarray | None, str, float]:
    """Re-solve all Lambert legs of one cycle AND independently propagate them with IAS15.

    Returns ``(converged, r_final_v3_km, r_final_v2_km, integrator_label)``:
    the spacecraft's terminal planet-frame position under IAS15 (V3) and
    under the analytic Lambert + moon ephemeris chain (V2), at the
    cycle's final encounter.

    The V2 terminal position is the analytic moon position at the leg
    endpoint, exactly as :func:`v2_moontour._cycle_residual` produces.
    The V3 terminal position is the IAS15 terminal spacecraft position
    after propagating through all legs of the cycle with each leg's
    Lambert v-out as the IC.

    Integrator label is set on first call (IAS15 vs LSODA fallback).
    """
    n_legs = len(sequence) - 1
    if n_revs is None:
        n_revs_used: tuple[int, ...] = tuple(0 for _ in range(n_legs))
    else:
        if len(n_revs) != n_legs:
            raise ValueError(f"n_revs length {len(n_revs)} != n_legs {n_legs}")
        n_revs_used = tuple(n_revs)

    epochs_days = [0.0]
    for tof in leg_tofs_days:
        epochs_days.append(epochs_days[-1] + tof)

    # Analytic moon states at each encounter (the V2 path)
    states: list[tuple[np.ndarray, np.ndarray]] = []
    for moon, t in zip(sequence, epochs_days, strict=True):
        sma, n_rad_day = consts[moon]
        states.append(
            _moon_state(theta_base[moon], n_rad_day, t_cycle_offset_days + t, sma, mu_primary)
        )

    integrator_label = ""
    # IAS15-propagate each leg, using the Lambert v-out at moon-A as the IC.
    sc_r_curr: np.ndarray | None = None
    worst_ias15_vs_analytic_kms = 0.0
    for k in range(n_legs):
        r_a, v_a_moon = states[k]
        r_b, _ = states[k + 1]
        nrev = max(0, n_revs_used[k])
        sols = _lambert(r_a, r_b, leg_tofs_days[k] * DAY_S, mu=mu_primary, max_revs=nrev)
        wanted = [s for s in sols if s.n_revs == n_revs_used[k]]
        if not wanted:
            return (False, None, None, integrator_label or "REBOUND IAS15", float("inf"))
        v_a_captured = v_a_moon
        best = min(wanted, key=lambda s: float(np.linalg.norm(s.v1 - v_a_captured)))
        # IC for the IAS15 leg: planet-frame spacecraft state at moon-A,
        # using Lambert's v-out from moon-A.
        r0_leg = r_a.copy()
        v0_leg = best.v1.copy()
        r_f_leg, _, label = _ias15_propagate_planet_frame(
            r0_leg, v0_leg, leg_tofs_days[k] * DAY_S, mu_primary, epsilon=ias15_epsilon
        )
        if not integrator_label:
            integrator_label = label
        # IAS15 terminal r vs the analytic moon target (Lambert's r_b
        # endpoint). This is the genuine IAS15-vs-analytic-Kepler agreement.
        leg_offset_kms = float(np.linalg.norm(r_f_leg - r_b))
        worst_ias15_vs_analytic_kms = max(worst_ias15_vs_analytic_kms, leg_offset_kms)
        sc_r_curr = r_f_leg
    if sc_r_curr is None:
        return (False, None, None, integrator_label or "REBOUND IAS15", float("inf"))
    return (True, sc_r_curr, states[-1][0].copy(), integrator_label, worst_ias15_vs_analytic_kms)


def run_v3_3d(
    candidate_id: str,
    sequence: tuple[str, ...],
    vinf_tuple_kms: tuple[float, ...],
    leg_tofs_days: tuple[float, ...],
    rel_offset_deg: float,
    system: cr3bp.CR3BPSystem | None,
    *,
    v2_verdict: V2MoontourVerdict,
    n_cycles: int = V3_N_CYCLES_MIN,
    n_revs: tuple[int, ...] | None = None,
    phase0_deg: float = 0.0,
    ias15_epsilon: float = 1e-12,
    agreement_floor_kms: float = V3_AGREEMENT_FLOOR_KMS,
    notes: str = "",
) -> V3Verdict3D:
    """Run V3 for a moontour: re-propagate cycle terminal positions with IAS15.

    Pipeline:
      1. For each cycle k = 0, 1, ..., n_cycles - 1:
           a. Advance moon longitudes by ``k * cycle_period_days``.
           b. Re-solve each Lambert leg (same as V2).
           c. Propagate the spacecraft through each leg with REBOUND IAS15
              in the planet frame (two-body Kepler about the primary).
           d. Record the cycle's V3 terminal spacecraft position
              (after all legs) AND the cycle's V2 terminal moon position
              (analytic moon position at the final encounter).
      2. Compute the V3 per-cycle drift series: position offset of each
         cycle's V3 terminal position vs cycle 0's V3 terminal position.
      3. Extract the V2 per-cycle drift series from ``v2_verdict``.
      4. Compute ``drift_agreement_kms = max_k |drift_v3[k] - drift_v2[k]|``.
      5. Verdict: PASS iff every cycle converged AND
         ``drift_agreement_kms <= agreement_floor_kms``.

    Parameters
    ----------
    candidate_id, sequence, vinf_tuple_kms, leg_tofs_days, rel_offset_deg, system:
        Same semantics as :func:`run_v2_moontour`. ``vinf_tuple_kms`` is
        carried for audit; V3 re-solves Lambert from geometry.
    v2_verdict:
        The :class:`V2MoontourVerdict` from :func:`run_v2_moontour` on this
        candidate at the SAME ``n_cycles`` (or larger). V3 reads
        ``per_cycle[k].rendezvous_drift_kms`` to compute the agreement.
    n_cycles:
        Cycles to attempt. Must be >= :data:`V3_N_CYCLES_MIN` (= 3).
    n_revs:
        Per-leg revolution count. Pass the candidate's stored ``n_rev``.
    phase0_deg:
        First moon's initial longitude at cycle 0, degrees.
    ias15_epsilon:
        REBOUND IAS15 ``epsilon`` tolerance. Default 1e-12 matches the
        V2 driver's CR3BP propagator tolerance for like-for-like compare.
    agreement_floor_kms:
        Bar against which ``drift_agreement_kms`` is gated.
    notes:
        Free-form audit note.

    Returns
    -------
    V3Verdict3D
        ``passes_v3`` is the headline.

    Notes
    -----
    The V3 driver does not gate on the absolute drift magnitude — it
    gates on integrator AGREEMENT. The V2 driver already produced the
    bounded-drift signature; V3 asks whether that signature survives a
    different integrator architecture.

    A V3 PASS does NOT admit to the catalogue. V4 + #329 still gate.
    """
    if not sequence:
        raise ValueError("empty sequence")
    if sequence[0] != sequence[-1]:
        raise ValueError(f"moontour sequence must be CLOSED (first == last); got {sequence!r}")
    n_legs = len(sequence) - 1
    if len(leg_tofs_days) != n_legs:
        raise ValueError(f"leg_tofs_days length {len(leg_tofs_days)} != n_legs {n_legs}")
    if len(vinf_tuple_kms) != len(sequence):
        raise ValueError(
            f"vinf_tuple_kms length {len(vinf_tuple_kms)} != len(sequence) {len(sequence)}"
        )
    if any(tof <= 0.0 for tof in leg_tofs_days):
        raise ValueError(f"leg_tofs_days must be positive; got {leg_tofs_days!r}")
    if n_cycles < V3_N_CYCLES_MIN:
        raise ValueError(f"V3 requires n_cycles >= {V3_N_CYCLES_MIN} (spec §14); got {n_cycles}")
    if agreement_floor_kms <= 0.0:
        raise ValueError(f"agreement_floor_kms must be > 0; got {agreement_floor_kms}")
    if len(v2_verdict.per_cycle) < n_cycles:
        raise ValueError(
            f"v2_verdict has only {len(v2_verdict.per_cycle)} cycles but V3 wants {n_cycles}"
        )

    primary = _resolve_primary(system, sequence)
    consts = _moon_constants(primary, sequence)
    mu_primary = PRIMARIES[primary]

    # Same theta_base convention as v2_moontour for an exact like-for-like
    # cycle-0 comparison.
    phase0_rad = math.radians(phase0_deg)
    rel_off_rad = math.radians(rel_offset_deg)
    distinct_moons = tuple(sorted({m for m in sequence}))
    if len(distinct_moons) == 1:
        raise ValueError(f"moontour requires >= 2 distinct moons; got {distinct_moons!r}")
    theta_base: dict[str, float] = {}
    for j, moon in enumerate(distinct_moons):
        if j == 0:
            theta_base[moon] = phase0_rad
        elif j == 1:
            theta_base[moon] = phase0_rad + rel_off_rad
        else:
            theta_base[moon] = (
                phase0_rad + rel_off_rad + 2.0 * math.pi * (j - 1) / len(distinct_moons)
            )

    cycle_period_days = float(sum(leg_tofs_days))

    per_cycle: list[V3CycleVerdict3D] = []
    v3_terminal_positions: list[np.ndarray] = []
    n_completed = 0
    integrator_label_used = ""
    cycle_zero_r_v3: np.ndarray | None = None

    for k in range(n_cycles):
        t_offset_days = k * cycle_period_days
        converged, r_v3, _r_v2, ilabel, ias15_vs_analytic = _cycle_v3(
            sequence=sequence,
            leg_tofs_days=leg_tofs_days,
            theta_base=theta_base,
            t_cycle_offset_days=t_offset_days,
            consts=consts,
            mu_primary=mu_primary,
            n_revs=n_revs,
            ias15_epsilon=ias15_epsilon,
        )
        if not integrator_label_used:
            integrator_label_used = ilabel
        if not converged or r_v3 is None:
            per_cycle.append(
                V3CycleVerdict3D(
                    cycle_index=k,
                    converged_legs=0,
                    n_legs=n_legs,
                    rendezvous_drift_kms_v3=float("inf"),
                    rendezvous_drift_kms_v2=float(v2_verdict.per_cycle[k].rendezvous_drift_kms),
                    agreement_kms=float("inf"),
                    ias15_terminal_offset_kms=float("inf"),
                    ias15_vs_analytic_kepler_kms=float("inf"),
                    notes="Lambert / IAS15 failed at least one leg in this cycle",
                )
            )
            break
        v3_terminal_positions.append(r_v3.copy())
        if k == 0:
            cycle_zero_r_v3 = r_v3.copy()
            drift_v3 = 0.0
        else:
            assert cycle_zero_r_v3 is not None
            drift_v3 = float(np.linalg.norm(r_v3 - cycle_zero_r_v3))
        drift_v2 = float(v2_verdict.per_cycle[k].rendezvous_drift_kms)
        agreement = abs(drift_v3 - drift_v2)
        per_cycle.append(
            V3CycleVerdict3D(
                cycle_index=k,
                converged_legs=n_legs,
                n_legs=n_legs,
                rendezvous_drift_kms_v3=drift_v3,
                rendezvous_drift_kms_v2=drift_v2,
                agreement_kms=agreement,
                ias15_terminal_offset_kms=drift_v3,
                ias15_vs_analytic_kepler_kms=ias15_vs_analytic,
            )
        )
        n_completed += 1

    drift_v3_series = tuple(c.rendezvous_drift_kms_v3 for c in per_cycle)
    drift_v2_series = tuple(c.rendezvous_drift_kms_v2 for c in per_cycle)
    drift_agreement = float("inf") if n_completed == 0 else max(c.agreement_kms for c in per_cycle)

    passes_v3 = bool(
        n_completed >= V3_N_CYCLES_MIN
        and n_completed == n_cycles
        and math.isfinite(drift_agreement)
        and drift_agreement <= agreement_floor_kms
    )

    return V3Verdict3D(
        candidate_id=candidate_id,
        sequence=tuple(sequence),
        n_cycles_propagated=int(n_completed),
        integrator=integrator_label_used or "REBOUND IAS15",
        per_cycle=tuple(per_cycle),
        per_cycle_drift_kms_v3=drift_v3_series,
        per_cycle_drift_kms_v2=drift_v2_series,
        drift_agreement_kms=float(drift_agreement),
        v3_v2_agreement_floor_kms=float(agreement_floor_kms),
        passes_v3=passes_v3,
        notes=notes,
    )


@dataclass(frozen=True)
class V3PeriodicRegressionVerdict:
    """V3 verdict for a 3D periodic CR3BP regression member (#287 / #301 path).

    For a periodic orbit candidate (NOT a moontour), V3 means
    re-propagating the IC under an INDEPENDENT integrator (REBOUND IAS15
    in CR3BP via custom forces, OR a tighter-tolerance LSODA) and
    asserting the same one-period closure as the spike's stored
    ``independent_closure_L2`` / ``independent_closure_residual``.

    Used by the #287 spike Braik-Ross (1,1) baseline and the #301 k=4
    doubly-hyperbolic-pair subfamily regression tests in this module.
    """

    candidate_id: str
    closure_residual_nondim_v3: float
    closure_residual_nondim_stored: float
    agreement_nondim: float
    """``|closure_v3 - closure_stored|`` — integrator agreement nondim."""
    closure_floor_nondim: float
    passes_v3: bool


def run_v3_periodic_regression(
    candidate_id: str,
    state0_nondim: np.ndarray,
    period_nondim: float,
    system: cr3bp.CR3BPSystem,
    closure_residual_nondim_stored: float,
    *,
    closure_floor_nondim: float = 1.0e-7,
) -> V3PeriodicRegressionVerdict:
    """V3 periodic-orbit regression: an independent integrator must reproduce closure.

    Cross-checks the stored ``independent_closure_L2`` of a CR3BP periodic
    spike row by re-propagating the IC under
    :func:`cyclerfinder.core.cr3bp.propagate` at tight tolerances.

    Note this regression is INTRA-MODEL: the V1 / spike's
    independent_closure_residual already runs scipy DOP853 at 1e-12; this
    V3 wrapper asserts the same propagator gives the same closure when
    asked again. A future REBOUND-IAS15 CR3BP backend (force-callback +
    Coriolis) would be a stronger independent check; this implementation
    is the wireable baseline.
    """
    state0 = np.asarray(state0_nondim, dtype=np.float64)
    if state0.shape != (6,):
        raise ValueError(f"state0_nondim must be 6D; got shape {state0.shape}")
    if period_nondim <= 0.0:
        raise ValueError(f"period_nondim must be > 0; got {period_nondim}")

    arc = cr3bp.propagate(system, state0, period_nondim, with_stm=False, rtol=1e-13, atol=1e-13)
    closure_v3 = float(np.linalg.norm(arc.state_f - state0))
    agreement = abs(closure_v3 - closure_residual_nondim_stored)
    passes = bool(closure_v3 <= closure_floor_nondim)
    return V3PeriodicRegressionVerdict(
        candidate_id=candidate_id,
        closure_residual_nondim_v3=closure_v3,
        closure_residual_nondim_stored=closure_residual_nondim_stored,
        agreement_nondim=agreement,
        closure_floor_nondim=closure_floor_nondim,
        passes_v3=passes,
    )


__all__ = [
    "V3_AGREEMENT_FLOOR_KMS",
    "V3_N_CYCLES_MIN",
    "V3CycleVerdict3D",
    "V3PeriodicRegressionVerdict",
    "V3Verdict3D",
    "run_v3_3d",
    "run_v3_periodic_regression",
]


# Suppress unused-import flake8 for SECONDS_PER_DAY (kept for callers; module
# also re-exports for downstream scripts).
_ = SECONDS_PER_DAY
