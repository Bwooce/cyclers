"""#650: inter-cycler transfer-compatibility network — pure closed-form.

Implements the design spec in
``docs/notes/2026-07-19-650-transfer-network-design.md`` in full: node
eligibility (spec §2), body-name resolution (§4), the same-body
powered-flyby ``dv_hop_kms`` cost metric and its B0-B3 bands (§3), ballistic
bend metadata, sphere-of-influence (``r_SOI``) computation, and the
statistical phase-alignment model over an explicitly-unknown relative phase
offset (§5). Every threshold below is a design-doc convention, named as
such; see the design doc for provenance.

Pure closed-form arithmetic over already-catalogued data: no integrators, no
network access, no catalogue writes. ``scripts/run_650_transfer_network.py``
is the driver that reads ``data/catalogue.yaml``, runs the pairwise sweep,
and writes the ``data/found/650_transfer_network/`` artifact.

Physical model (design §3): a same-body *powered-flyby handoff* — the taxi
departs cycler A on A's approach hyperbola at shared body X and performs a
single impulsive burn at periapsis to leave on cycler B's departure
hyperbola. This is the Oberth-optimal way to change ``v_inf`` magnitude, and
the only cost computable from magnitude-only catalogue data (no V_inf
vectors are recorded, so ``dv_hop_kms`` is a LOWER BOUND on the true
single-impulse cost — trustworthy as a negative, not as a positive; design
§3 "Semantics").

The phase-alignment model (design §5) is statistical, not deterministic,
because the catalogue records NO absolute encounter epochs or relative
phases for any cycler x cycler pair (``epoch_locked=false`` is a class
invariant for ``cycler`` rows). The one deterministic exception is the
window-intersection check between two ``quasi_cycler`` rows that both carry
a real ``validity_window`` (design §5 "Epoch-locked special case").
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from math import sqrt
from typing import Any, Literal

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import AU_KM, DAYS_PER_JULIAN_YEAR, MU_SUN_KM3_S2, PLANETS
from cyclerfinder.core.flyby import max_bend
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.search.tisserand import vinf_to_tisserand

# ---------------------------------------------------------------------------
# Conventions (design doc §2-§6; every value here is either a physical
# formula or an explicitly-labelled convention, echoed verbatim into
# summary.json by the driver script).
# ---------------------------------------------------------------------------

#: orbit_class values eligible to be a graph node (design §2 rule 1).
ELIGIBLE_ORBIT_CLASSES: frozenset[str] = frozenset({"cycler", "quasi_cycler"})

#: Abstract binary-genome placeholder bodies excluded from encounters (§2 rule 2).
_ABSTRACT_BODIES: frozenset[str] = frozenset({"P1", "P2"})

#: Heliocentric bodies with a Tisserand-parameter identity (§3 metadata).
HELIOCENTRIC_BODIES: frozenset[str] = frozenset({"V", "E", "M"})

#: B0 band: same Tisserand contour to within catalogue rounding precision (§3, convention).
B0_DELTA_VINF_KMS: float = 0.1
#: B1 band ceiling: order of a typical DSM budget (§3, convention).
B1_DV_HOP_KMS: float = 0.5
#: B2 band ceiling (§3, convention).
B2_DV_HOP_KMS: float = 2.0

#: Compute-gating rule: only run the phase model for dv_hop <= this (§5, convention)...
DV_PHASE_GATE_KMS: float = 1.0
#: ...OR any body outside the heliocentric set (moon-system edges, §5).
_HELIOCENTRIC_GATE_BODIES: frozenset[str] = HELIOCENTRIC_BODIES

#: Statistical model horizon for non-epoch-locked pairs, years (§5, convention).
STATISTICAL_HORIZON_YEARS: float = 100.0
#: Number of uniform delta0 samples in [0, T_B) (§5, deterministic algorithm parameter).
N_DELTA0_SAMPLES: int = 720
#: Sensitivity-grid relaxed windows, days (§5) -- reporting only, never used for cheap_edge.
SENSITIVITY_WINDOWS_DAYS: tuple[float, ...] = (1.0, 10.0, 30.0)
#: p_align threshold separating "recurrent" (waiting-time) from "phase_locked" (§5).
PHASE_LOCK_P_ALIGN_THRESHOLD: float = 0.99
#: Segment-chain-consistency tolerance vs the row's own period (encounter-offset derivation).
_SEGMENT_PERIOD_TOLERANCE_FRAC: float = 0.05

#: Mission-lifetime-scale cheap-edge wait ceiling, years (§6, convention).
CHEAP_EDGE_WAIT_YEARS: float = 20.0
#: Duty-cycle-adjusted alignment floor for the epoch_window_overlap cheap-edge branch (§6).
CHEAP_EDGE_DUTY_ADJUSTED_P_ALIGN: float = 0.5

#: Every catalogue row's V_inf values are magnitude-only (no direction data recorded).
DIRECTION_DATA_STATUS: Literal["absent"] = "absent"

PhaseStatus = Literal[
    "recurrent",
    "phase_locked",
    "no_period_data",
    "not_computed_dv_gated",
    "epoch_disjoint",
    "epoch_window_overlap",
]

Band = Literal["B0_ballistic_compatible", "B1_cheap", "B2_moderate", "B3_expensive"]


# ---------------------------------------------------------------------------
# §4 -- body-name resolution
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class BodyInfo:
    """Physical-constants record for a resolved catalogue body."""

    name: str
    """Canonical full name (e.g. ``"Earth"``, ``"Titan"``) -- the edge graph's body identity."""
    code: str
    """Lookup key usable with :func:`cyclerfinder.search.tisserand.vinf_to_tisserand`
    (a PLANETS code such as ``"E"`` for heliocentric bodies, else the same as ``name``)."""
    kind: Literal["planet", "satellite"]
    mu_km3_s2: float
    radius_eq_km: float
    safe_alt_km: float
    sma_km: float
    """Semi-major axis about the immediate primary, km (heliocentric AU->km for planets;
    already km, already about-primary, for satellites)."""
    parent_mu_km3_s2: float
    """GM of the immediate primary (Sun for a planet, the primary planet's system GM
    from :data:`cyclerfinder.core.satellites.PRIMARIES` for a satellite)."""

    @property
    def r_p_km(self) -> float:
        """Minimum safe flyby periapsis radius: planet/moon radius + safe altitude."""
        return self.radius_eq_km + self.safe_alt_km


def resolve_body(raw_name: str) -> BodyInfo:
    """Resolve a catalogue body string to its physical-constants record (design §4).

    Resolution order: PLANETS by ``code`` (``"E"``, ``"M"``, ``"V"``, ``"Pl"``, ...),
    then PLANETS by full ``name`` (``"Pluto"``), then SATELLITES by name (``"Moon"``,
    ``"Io"``, ... ``"Charon"``). Raises :class:`KeyError` on anything unresolvable --
    after the caller excludes the abstract ``P1``/``P2`` placeholders, nothing in the
    eligible catalogue subset should be unresolvable (design §1 body-vocabulary survey).
    """
    planet = PLANETS.get(raw_name)
    if planet is not None:
        return BodyInfo(
            name=planet.name,
            code=planet.code,
            kind="planet",
            mu_km3_s2=planet.mu_km3_s2,
            radius_eq_km=planet.radius_eq_km,
            safe_alt_km=planet.safe_alt_km,
            sma_km=planet.sma_au * AU_KM,
            parent_mu_km3_s2=MU_SUN_KM3_S2,
        )
    for candidate in PLANETS.values():
        if candidate.name == raw_name:
            return BodyInfo(
                name=candidate.name,
                code=candidate.code,
                kind="planet",
                mu_km3_s2=candidate.mu_km3_s2,
                radius_eq_km=candidate.radius_eq_km,
                safe_alt_km=candidate.safe_alt_km,
                sma_km=candidate.sma_au * AU_KM,
                parent_mu_km3_s2=MU_SUN_KM3_S2,
            )
    satellite = SATELLITES.get(raw_name)
    if satellite is not None:
        parent_mu = PRIMARIES[satellite.primary]
        return BodyInfo(
            name=satellite.name,
            code=satellite.name,
            kind="satellite",
            mu_km3_s2=satellite.mu_km3_s2,
            radius_eq_km=satellite.radius_eq_km,
            safe_alt_km=satellite.safe_alt_km,
            sma_km=satellite.sma_km,
            parent_mu_km3_s2=parent_mu,
        )
    raise KeyError(f"transfer_network.resolve_body: unresolvable body {raw_name!r}")


# ---------------------------------------------------------------------------
# §2 -- eligibility
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Encounter:
    """One usable (non-null-vinf, non-abstract-body) catalogue encounter entry."""

    body_raw: str
    body: str
    """Canonical body name (:attr:`BodyInfo.name`) -- the graph's body identity."""
    vinf_kms: float


def usable_encounters(row: dict[str, Any]) -> list[Encounter]:
    """Non-null-vinf, non-``P1``/``P2`` encounters, body-resolved (design §2 rule 2)."""
    out: list[Encounter] = []
    for entry in row.get("vinf_kms_at_encounters") or []:
        body_raw = entry.get("body")
        vinf = entry.get("vinf_kms")
        if vinf is None or body_raw in _ABSTRACT_BODIES or body_raw is None:
            continue
        resolved = resolve_body(body_raw)
        out.append(Encounter(body_raw=body_raw, body=resolved.name, vinf_kms=float(vinf)))
    return out


def is_node(row: dict[str, Any]) -> bool:
    """True iff ``row`` is an eligible graph node (design §2)."""
    if row.get("orbit_class") not in ELIGIBLE_ORBIT_CLASSES:
        return False
    return len(usable_encounters(row)) > 0


def usable_bodies(row: dict[str, Any]) -> set[str]:
    """Set of canonical body names ``row`` has >=1 usable encounter at."""
    return {enc.body for enc in usable_encounters(row)}


# ---------------------------------------------------------------------------
# §3 -- dv_hop cost metric, bands, bend metadata
# ---------------------------------------------------------------------------


def periapsis_speed_kms(vinf_kms: float, mu_km3_s2: float, r_p_km: float) -> float:
    """Hyperbolic periapsis speed for excess speed ``vinf_kms`` at floor ``r_p_km``."""
    return sqrt(vinf_kms * vinf_kms + 2.0 * mu_km3_s2 / r_p_km)


def dv_hop_kms(vinf_a_kms: float, vinf_b_kms: float, mu_km3_s2: float, r_p_km: float) -> float:
    """Oberth-optimal same-body powered-flyby handoff cost (design §3).

    ``dv = |sqrt(v_a^2 + 2*mu/r_p) - sqrt(v_b^2 + 2*mu/r_p)|``.
    """
    return abs(
        periapsis_speed_kms(vinf_a_kms, mu_km3_s2, r_p_km)
        - periapsis_speed_kms(vinf_b_kms, mu_km3_s2, r_p_km)
    )


def classify_band(delta_vinf_kms: float, dv_hop_kms_value: float) -> Band:
    """B0-B3 band classification (design §3). Checks are mutually exclusive & exhaustive.

    B0 is checked first on ``delta_vinf_kms`` (not ``dv_hop_kms``): since the Oberth
    map ``v -> sqrt(v^2 + c)`` has slope < 1 everywhere (``c = 2*mu/r_p > 0``),
    ``dv_hop_kms <= delta_vinf_kms`` always, so every B0 edge is automatically <= 0.1
    km/s on dv_hop too -- the B0 test is deliberately the tighter, Tisserand-native
    ``delta_vinf_kms`` gate per the design's own semantics.
    """
    if delta_vinf_kms <= B0_DELTA_VINF_KMS:
        return "B0_ballistic_compatible"
    if dv_hop_kms_value <= B1_DV_HOP_KMS:
        return "B1_cheap"
    if dv_hop_kms_value <= B2_DV_HOP_KMS:
        return "B2_moderate"
    return "B3_expensive"


def bend_max_deg(vinf_kms: float, r_p_km: float, mu_km3_s2: float) -> float:
    """Max ballistic deflection (deg) at ``vinf_kms``, floor ``r_p_km`` (reuses
    :func:`cyclerfinder.core.flyby.max_bend`, which returns radians)."""
    return float(np.degrees(max_bend(mu_km3_s2, r_p_km, vinf_kms)))


def r_soi_km(body: BodyInfo) -> float:
    """Sphere-of-influence radius, km: ``r_SOI = a_X * (mu_X / mu_parent)^(2/5)`` (design §5)."""
    return float(body.sma_km * (body.mu_km3_s2 / body.parent_mu_km3_s2) ** 0.4)


def w_handoff_days(body: BodyInfo, vinf_a_kms: float, vinf_b_kms: float) -> float:
    """Same-SOI-simultaneously rendezvous window, days (design §5).

    ``w_handoff = 2 * r_SOI(X) / min(v_inf_A, v_inf_B)``, converted from seconds to days.
    """
    v_min = min(vinf_a_kms, vinf_b_kms)
    seconds = 2.0 * r_soi_km(body) / v_min
    return seconds / 86400.0


# ---------------------------------------------------------------------------
# §5 -- period + intra-cycle encounter offsets ("phase timeline")
# ---------------------------------------------------------------------------

PhaseTimeline = Literal["segments_derived", "uniform_assumed"]


def period_days(row: dict[str, Any]) -> float | None:
    """Repeat period T, days: ``period.years`` (converted), else (quasi_cycler)
    ``validity_window.synodic_period_days``, else ``None`` (design §5)."""
    period = row.get("period") or {}
    years = period.get("years")
    if years is not None:
        return float(years) * DAYS_PER_JULIAN_YEAR
    vw = row.get("validity_window") or {}
    synodic = vw.get("synodic_period_days")
    if synodic is not None:
        return float(synodic)
    return None


def encounter_offsets_days(
    row: dict[str, Any], body: str, t_days: float | None
) -> tuple[list[float], PhaseTimeline]:
    """Intra-cycle encounter offsets phi_i for ``body`` within ``row``'s own T (design §5).

    Tries a ``trajectory.segments[].tof_days`` cumulative walk first; falls back to
    spreading the ``n`` encounters at ``body`` uniformly over T when the segment chain
    is absent, incomplete, or inconsistent with the row's own encounter list (per-body
    encounter counts must match, and the summed segment ToF must land within 5% of T).
    """
    encounters = [enc for enc in usable_encounters(row) if enc.body == body]
    n = len(encounters)
    if n == 0:
        return [], "uniform_assumed"

    segments = ((row.get("trajectory") or {}).get("segments")) or []
    if segments and t_days:
        offsets_by_body: dict[str, list[float]] = {}
        cumulative = 0.0
        valid = True
        for seg in segments:
            tof = seg.get("tof_days")
            to_raw = seg.get("to")
            if tof is None or to_raw is None:
                valid = False
                break
            try:
                to_canonical = resolve_body(to_raw).name
            except KeyError:
                valid = False
                break
            cumulative += float(tof)
            offsets_by_body.setdefault(to_canonical, []).append(cumulative)
        if valid and cumulative > 0.0:
            seg_counts = {b: len(v) for b, v in offsets_by_body.items()}
            enc_counts: dict[str, int] = {}
            for enc in usable_encounters(row):
                enc_counts[enc.body] = enc_counts.get(enc.body, 0) + 1
            period_ok = abs(cumulative - t_days) <= _SEGMENT_PERIOD_TOLERANCE_FRAC * t_days
            if seg_counts == enc_counts and period_ok:
                phis = offsets_by_body.get(body, [])
                if len(phis) == n:
                    return sorted(p % t_days for p in phis), "segments_derived"

    if t_days:
        return [i * t_days / n for i in range(n)], "uniform_assumed"
    return [0.0] * n, "uniform_assumed"


# ---------------------------------------------------------------------------
# §5 -- epoch-window intersection (quasi_cycler x quasi_cycler special case)
# ---------------------------------------------------------------------------


def _parse_iso(value: str) -> datetime:
    return datetime.fromisoformat(value.replace("Z", "+00:00"))


def epoch_window_intersection(
    row_a: dict[str, Any], row_b: dict[str, Any]
) -> tuple[datetime, datetime] | None:
    """Intersection ``[max(starts), min(ends)]`` of two rows' ``validity_window``.

    Returns ``None`` when either row lacks a window, or the intersection is empty
    (disjoint windows) -- the only *deterministic* phase negative this data permits
    (design §5 "Epoch-locked special case").
    """
    vw_a = row_a.get("validity_window") or {}
    vw_b = row_b.get("validity_window") or {}
    start_a, end_a = vw_a.get("start"), vw_a.get("end")
    start_b, end_b = vw_b.get("start"), vw_b.get("end")
    if start_a is None or end_a is None or start_b is None or end_b is None:
        return None
    start = max(_parse_iso(start_a), _parse_iso(start_b))
    end = min(_parse_iso(end_a), _parse_iso(end_b))
    if start >= end:
        return None
    return start, end


# ---------------------------------------------------------------------------
# §5 -- statistical phase-alignment model
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class WindowResult:
    """Phase-alignment statistics at one handoff window ``w_days`` (design §5)."""

    w_days: float
    p_align: float
    median_wait_years: float | None
    p90_wait_years: float | None

    def to_json(self) -> dict[str, Any]:
        return {
            "w_days": self.w_days,
            "p_align": self.p_align,
            "median_wait_years": self.median_wait_years,
            "p90_wait_years": self.p90_wait_years,
        }


def _delta0_grid(t_b_days: float, n_delta0: int) -> NDArray[np.float64]:
    """Uniform delta0 sample grid in [0, T_B), the shared grid used by
    :func:`phase_alignment_stats` and (for exact test cross-checks) callers of
    :func:`brute_force_phase_alignment`."""
    return np.linspace(0.0, t_b_days, n_delta0, endpoint=False)


def _raw_wait_days(
    phi_a_days: list[float],
    t_a_days: float,
    phi_b_days: list[float],
    t_b_days: float,
    w_days: float,
    horizon_days: float,
    delta0: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Per-``delta0``-sample earliest coincidence time (days), vectorized over ``j``.

    For each ``delta0`` and each intra-cycle offset pair ``(phi_a, phi_b)``, finds the
    earliest ``t = phi_a + j*T_A`` (``j = 0, 1, ...``) within a B-encounter window
    ``w_days`` of some ``delta0 + phi_b + m*T_B`` (``m = 0, 1, ...``), up to
    ``horizon_days``. Vectorized over ``j`` (centered-mod distance to the nearest
    non-negative-``m`` B-encounter), per design §5's own implementer note.
    ``np.inf`` where no coincidence occurs within the horizon.
    """
    n_delta0 = delta0.shape[0]
    best_wait: NDArray[np.float64] = np.full(n_delta0, np.inf)

    n_j = int(horizon_days / t_a_days) + 2
    j_index = np.arange(n_j)

    for phi_a in phi_a_days:
        t_j = phi_a + j_index * t_a_days
        t_j = t_j[t_j <= horizon_days]
        if t_j.size == 0:
            continue
        for phi_b in phi_b_days:
            diff = t_j[:, None] - delta0[None, :] - phi_b  # shape (n_t_j, n_delta0)
            m = np.maximum(0.0, np.round(diff / t_b_days))
            resid = diff - m * t_b_days
            satisfies = np.abs(resid) <= w_days
            any_sat = satisfies.any(axis=0)
            first_idx = satisfies.argmax(axis=0)
            wait_for_pair = np.where(any_sat, t_j[first_idx], np.inf)
            best_wait = np.minimum(best_wait, wait_for_pair)

    return best_wait


def phase_alignment_stats(
    phi_a_days: list[float],
    t_a_days: float,
    phi_b_days: list[float],
    t_b_days: float,
    w_days: float,
    horizon_days: float,
    n_delta0: int = N_DELTA0_SAMPLES,
) -> WindowResult:
    """Statistical phase-alignment model over an explicitly-unknown relative phase.

    Samples ``delta0`` at ``n_delta0`` uniform points in ``[0, T_B)`` (see
    :func:`_delta0_grid`) and reports, via :func:`_raw_wait_days`, the fraction of
    ``delta0`` samples achieving >=1 coincidence (``p_align``) and the median/p90
    wait time (years) over the succeeding samples.
    """
    delta0 = _delta0_grid(t_b_days, n_delta0)
    best_wait = _raw_wait_days(
        phi_a_days, t_a_days, phi_b_days, t_b_days, w_days, horizon_days, delta0
    )

    aligned = np.isfinite(best_wait)
    p_align = float(np.mean(aligned))
    if aligned.any():
        aligned_waits = best_wait[aligned]
        median_wait = float(np.median(aligned_waits) / DAYS_PER_JULIAN_YEAR)
        p90_wait = float(np.percentile(aligned_waits, 90) / DAYS_PER_JULIAN_YEAR)
    else:
        median_wait = None
        p90_wait = None
    return WindowResult(
        w_days=w_days, p_align=p_align, median_wait_years=median_wait, p90_wait_years=p90_wait
    )


def brute_force_phase_alignment(
    phi_a_days: list[float],
    t_a_days: float,
    phi_b_days: list[float],
    t_b_days: float,
    w_days: float,
    horizon_days: float,
    delta0_samples: NDArray[np.float64],
) -> NDArray[np.float64]:
    """O(events^2) reference implementation (test cross-check only, design §5).

    For each given ``delta0`` sample, naively enumerates every A-encounter time and
    every B-encounter time up to ``horizon_days`` (plus a ``w_days`` margin on the B
    side, since a B-encounter just past the horizon can still be the nearest match to
    an A-encounter just before it) and returns the smallest A-time at which some
    B-encounter is within ``w_days`` (``np.inf`` if none). Deliberately independent in
    structure from :func:`_raw_wait_days` (no ``j``-vectorization, no centered-mod
    arithmetic) -- exists to validate that function against a straightforward
    nested-loop implementation, not used by the production sweep.
    """
    a_times: list[float] = []
    for phi_a in phi_a_days:
        j = 0
        while True:
            t = phi_a + j * t_a_days
            if t > horizon_days:
                break
            a_times.append(t)
            j += 1
    a_times.sort()

    out = np.full(delta0_samples.shape[0], np.inf)
    for idx, delta0 in enumerate(delta0_samples):
        b_times: list[float] = []
        for phi_b in phi_b_days:
            m = 0
            while True:
                t = delta0 + phi_b + m * t_b_days
                if t > horizon_days + w_days:
                    break
                b_times.append(t)
                m += 1
        for t_a in a_times:
            if any(abs(t_a - t_b) <= w_days for t_b in b_times):
                out[idx] = t_a
                break
    return out


# ---------------------------------------------------------------------------
# Edge assembly -- ties §2-§6 together for one (row_a, row_b, body) candidate.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class BestPair:
    dv_hop_kms: float
    dv_hop_max_kms: float
    delta_vinf_kms: float
    vinf_a_kms: float
    vinf_b_kms: float


def _best_pair(
    encounters_a: list[Encounter], encounters_b: list[Encounter], mu_km3_s2: float, r_p_km: float
) -> BestPair:
    dvs: list[float] = []
    best_dv = float("inf")
    best_delta = 0.0
    best_va = 0.0
    best_vb = 0.0
    for enc_a in encounters_a:
        for enc_b in encounters_b:
            dv = dv_hop_kms(enc_a.vinf_kms, enc_b.vinf_kms, mu_km3_s2, r_p_km)
            delta = abs(enc_a.vinf_kms - enc_b.vinf_kms)
            dvs.append(dv)
            if dv < best_dv:
                best_dv = dv
                best_delta = delta
                best_va = enc_a.vinf_kms
                best_vb = enc_b.vinf_kms
    return BestPair(
        dv_hop_kms=best_dv,
        dv_hop_max_kms=max(dvs),
        delta_vinf_kms=best_delta,
        vinf_a_kms=best_va,
        vinf_b_kms=best_vb,
    )


@dataclass(frozen=True)
class PhaseResult:
    status: PhaseStatus
    timeline_a: PhaseTimeline | None
    timeline_b: PhaseTimeline | None
    period_a_days: float | None
    period_b_days: float | None
    windows: list[WindowResult] = field(default_factory=list)
    p_align_duty_adjusted: float | None = None

    def to_json(self) -> dict[str, Any]:
        out: dict[str, Any] = {
            "status": self.status,
            "timeline_a": self.timeline_a,
            "timeline_b": self.timeline_b,
            "period_a_days": self.period_a_days,
            "period_b_days": self.period_b_days,
            "windows": [w.to_json() for w in self.windows],
        }
        if self.p_align_duty_adjusted is not None:
            out["p_align_duty_adjusted"] = self.p_align_duty_adjusted
        return out


def _phase_gate_applies(dv_hop: float, body_raw_code: str) -> bool:
    """Design §5 compute-gating rule: dv_hop <= 1.0 km/s OR body outside {V, E, M}."""
    return dv_hop <= DV_PHASE_GATE_KMS or body_raw_code not in _HELIOCENTRIC_GATE_BODIES


def compute_phase(
    row_a: dict[str, Any],
    row_b: dict[str, Any],
    body: str,
    body_raw_code: str,
    body_info: BodyInfo,
    dv_hop: float,
    vinf_a_kms: float,
    vinf_b_kms: float,
) -> PhaseResult:
    """Full §5 phase-alignment computation for one edge (all branches)."""
    if not _phase_gate_applies(dv_hop, body_raw_code):
        return PhaseResult(
            status="not_computed_dv_gated",
            timeline_a=None,
            timeline_b=None,
            period_a_days=None,
            period_b_days=None,
        )

    t_a = period_days(row_a)
    t_b = period_days(row_b)

    epoch_locked_a = bool(row_a.get("epoch_locked")) and row_a.get("validity_window") is not None
    epoch_locked_b = bool(row_b.get("epoch_locked")) and row_b.get("validity_window") is not None

    w_handoff = w_handoff_days(body_info, vinf_a_kms, vinf_b_kms)
    all_w_days = (w_handoff, *SENSITIVITY_WINDOWS_DAYS)

    if epoch_locked_a and epoch_locked_b:
        intersection = epoch_window_intersection(row_a, row_b)
        if intersection is None:
            return PhaseResult(
                status="epoch_disjoint",
                timeline_a=None,
                timeline_b=None,
                period_a_days=t_a,
                period_b_days=t_b,
            )
        start, end = intersection
        horizon_days = (end - start).total_seconds() / 86400.0
        if t_a is None or t_b is None:
            return PhaseResult(
                status="no_period_data",
                timeline_a=None,
                timeline_b=None,
                period_a_days=t_a,
                period_b_days=t_b,
            )
        phi_a, timeline_a = encounter_offsets_days(row_a, body, t_a)
        phi_b, timeline_b = encounter_offsets_days(row_b, body, t_b)
        windows = [
            phase_alignment_stats(phi_a, t_a, phi_b, t_b, w, horizon_days) for w in all_w_days
        ]
        duty_a = (row_a.get("validity_window") or {}).get("synodic_duty_cycle_pct")
        duty_b = (row_b.get("validity_window") or {}).get("synodic_duty_cycle_pct")
        p_align_duty_adjusted = None
        if duty_a is not None and duty_b is not None:
            p_align_duty_adjusted = windows[0].p_align * (duty_a / 100.0) * (duty_b / 100.0)
        return PhaseResult(
            status="epoch_window_overlap",
            timeline_a=timeline_a,
            timeline_b=timeline_b,
            period_a_days=t_a,
            period_b_days=t_b,
            windows=windows,
            p_align_duty_adjusted=p_align_duty_adjusted,
        )

    if t_a is None or t_b is None:
        return PhaseResult(
            status="no_period_data",
            timeline_a=None,
            timeline_b=None,
            period_a_days=t_a,
            period_b_days=t_b,
        )

    phi_a, timeline_a = encounter_offsets_days(row_a, body, t_a)
    phi_b, timeline_b = encounter_offsets_days(row_b, body, t_b)
    horizon_days = STATISTICAL_HORIZON_YEARS * DAYS_PER_JULIAN_YEAR
    windows = [phase_alignment_stats(phi_a, t_a, phi_b, t_b, w, horizon_days) for w in all_w_days]
    status: PhaseStatus = (
        "recurrent" if windows[0].p_align >= PHASE_LOCK_P_ALIGN_THRESHOLD else "phase_locked"
    )
    return PhaseResult(
        status=status,
        timeline_a=timeline_a,
        timeline_b=timeline_b,
        period_a_days=t_a,
        period_b_days=t_b,
        windows=windows,
    )


def is_cheap_edge(band: Band, phase: PhaseResult) -> bool:
    """Design §6 cheap-edge classification.

    ``phase_locked`` (indeterminate) edges are NEVER cheap regardless of ``p_align`` --
    the honest census separates "cheap and realizable" from "cheap iff an unrecorded
    phase happens to align" (the ``cheap_dv_phase_indeterminate`` count, computed by
    the driver script).
    """
    if band not in ("B0_ballistic_compatible", "B1_cheap"):
        return False
    if not phase.windows:
        return False
    w_handoff_result = phase.windows[0]
    if phase.status == "recurrent":
        return (
            w_handoff_result.median_wait_years is not None
            and w_handoff_result.median_wait_years <= CHEAP_EDGE_WAIT_YEARS
        )
    if phase.status == "epoch_window_overlap":
        return (
            phase.p_align_duty_adjusted is not None
            and phase.p_align_duty_adjusted >= CHEAP_EDGE_DUTY_ADJUSTED_P_ALIGN
        )
    return False


@dataclass(frozen=True)
class Edge:
    id_a: str
    id_b: str
    body: str
    dv_hop_kms: float
    dv_hop_max_kms: float
    delta_vinf_kms: float
    band: Band
    vinf_a_kms: float
    vinf_b_kms: float
    tisserand_a: float | None
    tisserand_b: float | None
    bend_max_a_deg: float
    bend_max_b_deg: float
    r_p_km: float
    sense_a: str | None
    sense_b: str | None
    superseded_a: list[str] | None
    superseded_b: list[str] | None
    phase: PhaseResult
    cheap_edge: bool

    def to_json(self) -> dict[str, Any]:
        return {
            "id_a": self.id_a,
            "id_b": self.id_b,
            "body": self.body,
            "dv_hop_kms": self.dv_hop_kms,
            "dv_hop_max_kms": self.dv_hop_max_kms,
            "delta_vinf_kms": self.delta_vinf_kms,
            "band": self.band,
            "vinf_a_kms": self.vinf_a_kms,
            "vinf_b_kms": self.vinf_b_kms,
            "tisserand_a": self.tisserand_a,
            "tisserand_b": self.tisserand_b,
            "bend_max_a_deg": self.bend_max_a_deg,
            "bend_max_b_deg": self.bend_max_b_deg,
            "r_p_km": self.r_p_km,
            "direction_data": DIRECTION_DATA_STATUS,
            "sense_a": self.sense_a,
            "sense_b": self.sense_b,
            "superseded_a": self.superseded_a,
            "superseded_b": self.superseded_b,
            "phase": self.phase.to_json(),
            "cheap_edge": self.cheap_edge,
        }


def compute_edge(row_a: dict[str, Any], row_b: dict[str, Any], body: str) -> Edge:
    """Assemble the full edge record for shared canonical ``body`` between two nodes.

    ``row_a``/``row_b`` must already be ordered ``id_a < id_b`` by the caller (design
    §2: "unordered node pair (i < j by id)").
    """
    encounters_a = [enc for enc in usable_encounters(row_a) if enc.body == body]
    encounters_b = [enc for enc in usable_encounters(row_b) if enc.body == body]
    if not encounters_a or not encounters_b:
        raise ValueError(f"compute_edge: body {body!r} not usable on both rows")
    body_raw_code = encounters_a[0].body_raw
    body_info = resolve_body(body_raw_code)

    best = _best_pair(encounters_a, encounters_b, body_info.mu_km3_s2, body_info.r_p_km)
    band = classify_band(best.delta_vinf_kms, best.dv_hop_kms)

    tisserand_a: float | None = None
    tisserand_b: float | None = None
    if body_raw_code in HELIOCENTRIC_BODIES:
        tisserand_a = vinf_to_tisserand(body_info.code, best.vinf_a_kms)
        tisserand_b = vinf_to_tisserand(body_info.code, best.vinf_b_kms)

    min_vinf_a = min(enc.vinf_kms for enc in encounters_a)
    min_vinf_b = min(enc.vinf_kms for enc in encounters_b)
    bend_a = bend_max_deg(min_vinf_a, body_info.r_p_km, body_info.mu_km3_s2)
    bend_b = bend_max_deg(min_vinf_b, body_info.r_p_km, body_info.mu_km3_s2)

    phase = compute_phase(
        row_a,
        row_b,
        body,
        body_raw_code,
        body_info,
        best.dv_hop_kms,
        best.vinf_a_kms,
        best.vinf_b_kms,
    )
    cheap = is_cheap_edge(band, phase)

    return Edge(
        id_a=row_a["id"],
        id_b=row_b["id"],
        body=body,
        dv_hop_kms=best.dv_hop_kms,
        dv_hop_max_kms=best.dv_hop_max_kms,
        delta_vinf_kms=best.delta_vinf_kms,
        band=band,
        vinf_a_kms=best.vinf_a_kms,
        vinf_b_kms=best.vinf_b_kms,
        tisserand_a=tisserand_a,
        tisserand_b=tisserand_b,
        bend_max_a_deg=bend_a,
        bend_max_b_deg=bend_b,
        r_p_km=body_info.r_p_km,
        sense_a=row_a.get("sense"),
        sense_b=row_b.get("sense"),
        superseded_a=row_a.get("superseded_by"),
        superseded_b=row_b.get("superseded_by"),
        phase=phase,
        cheap_edge=cheap,
    )


def candidate_pairs(nodes: list[dict[str, Any]]) -> list[tuple[dict[str, Any], dict[str, Any]]]:
    """Unordered node pairs (``id_a < id_b``) sharing >=1 usable encounter body."""
    ordered = sorted(nodes, key=lambda r: r["id"])
    out: list[tuple[dict[str, Any], dict[str, Any]]] = []
    bodies_by_row = [usable_bodies(r) for r in ordered]
    n = len(ordered)
    for i in range(n):
        for j in range(i + 1, n):
            if bodies_by_row[i] & bodies_by_row[j]:
                out.append((ordered[i], ordered[j]))
    return out
