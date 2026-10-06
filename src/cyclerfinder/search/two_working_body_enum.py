"""Ideal-model enumeration for the two-working-body generator (#942 R1, #943 X1).

A *structure* is a cycle template: the sequence of legs (transfers, returns)
with the cycle period fixed at ``k`` synodic periods of the two bodies. Each
structure is solved by :func:`two_working_body.correct_dates` from a grid of
seeds (overall phase x split of the free Lambert time), the converged zeros
are de-duplicated (a shift by one synodic period is the same cycler rotated),
and every zero is assessed: flybys with minimax free directions, the #888/#937
demanded-turn gate at the project floor, the near-180-degree rejection, the
independent-propagation encounter check, and the heliocentric (or
planetocentric) extent.

Return catalogue at a body (circular model):

* full-revolution ``n:m`` resonant returns (``ResonantLeg``);
* half-revolution n-pi returns (``HalfRevLeg``; both mirror solutions);
* generic same-body returns (``LambertLeg`` with ``nrev >= 1``, both branches),
  which include Hollister's "symmetric" returns.

A massless target (Russell's architecture, R1 cell (c), the one-body X1
control) is placed on a solved one-body cycler by
:func:`place_massless_target`: the target's phase is chosen so that it sits
where a leg crosses its orbit, and that leg is split there; the corrector then
re-solves with the target's vector-continuity residual.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Iterator, Sequence
from dataclasses import dataclass, field, replace

import numpy as np

from cyclerfinder.core.lambert import LambertError, lambert
from cyclerfinder.search.two_working_body import (
    CircularSystem,
    Cycle,
    CycleReport,
    Flyby,
    HalfRevLeg,
    LambertLeg,
    Leg,
    ResonantLeg,
    Solution,
    Vec,
    _blocks,
    correct_dates,
    cycle_flybys,
    date_residual,
    encounter_self_consistency,
    eval_lambert_legs,
    fixed_duration_s,
    gate_cycle,
    kepler_step,
    rmin_radii,
    sphere_of_influence_km,
)

DAY = 86400.0


# ---------------------------------------------------------------------------
# Catalogue and structures
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CatalogueSpec:
    """Which returns a body may host."""

    resonances: tuple[tuple[int, int], ...] = ((1, 1), (2, 1), (1, 2), (3, 2), (2, 3))
    half_revs: tuple[tuple[int, int, bool], ...] = (
        (1, 0, True),
        (1, 0, False),
        (3, 1, True),
        (3, 1, False),
    )
    generic_revs: tuple[int, ...] = (1,)


def return_catalogue(body: str, spec: CatalogueSpec) -> list[Leg]:
    out: list[Leg] = [ResonantLeg(body, n, m) for n, m in spec.resonances]
    out += [HalfRevLeg(body, h, k, p, 0) for h, k, p in spec.half_revs]
    out += [LambertLeg(body, body, n, b) for n in spec.generic_revs for b in ("low", "high")]
    return out


def _min_lambert_tof(system: CircularSystem, leg: LambertLeg) -> float:
    """A loose lower bound on a Lambert leg's flight time (pruning only)."""
    if leg.frm == leg.to:
        return 0.5 * leg.nrev * system.period_s(leg.frm)
    p_in = min(system.period_s(leg.frm), system.period_s(leg.to))
    return 0.08 * p_in + leg.nrev * 0.5 * p_in


def structures(
    system: CircularSystem,
    a: str,
    b: str,
    k: int,
    *,
    max_returns: dict[str, int],
    spec: dict[str, CatalogueSpec],
    transfer_revs: Sequence[int] = (0, 1),
    visits: int = 1,
) -> Iterator[Cycle]:
    """Every cycle template ``[A-block, A->B, B-block, B->A] * visits`` with period
    ``k`` A-B synodic periods whose fixed time leaves room for the Lambert legs."""
    period = k * system.synodic_s(a, b)
    cat = {c: return_catalogue(c, spec[c]) for c in (a, b)}
    blocks: dict[str, list[tuple[Leg, ...]]] = {
        c: [
            tuple(s) for n in range(max_returns[c] + 1) for s in itertools.product(cat[c], repeat=n)
        ]
        for c in (a, b)
    }
    transfers: dict[tuple[str, str], list[tuple[Leg, ...]]] = {
        (f, t): [
            (LambertLeg(f, t, n, br),)
            for n in transfer_revs
            for br in (("single",) if n == 0 else ("low", "high"))
        ]
        for f, t in ((a, b), (b, a))
    }
    unit: list[list[tuple[Leg, ...]]] = [
        blocks[a],
        transfers[(a, b)],
        blocks[b],
        transfers[(b, a)],
    ]
    for combo in itertools.product(*(unit * visits)):
        legs: list[Leg] = []
        for part in combo:
            legs.extend(part)
        fixed = sum(fixed_duration_s(system, lg) for lg in legs if lg.fixed_tof)
        lam = sum(_min_lambert_tof(system, lg) for lg in legs if isinstance(lg, LambertLeg))
        if fixed + lam >= period:
            continue
        yield Cycle(tuple(legs), period)


# ---------------------------------------------------------------------------
# Solving
# ---------------------------------------------------------------------------


def _seed_grid(
    system: CircularSystem, cycle: Cycle, phase_period_s: float, n_phase: int, n_split: int
) -> Iterator[np.ndarray]:
    """Seeds: overall phase on ``[0, phase_period)`` x the free Lambert time split
    over the legs on a simplex grid."""
    li = cycle.lambert_index
    n = len(li)
    fixed_after = []
    for j in range(n):
        i0, i1 = li[j], li[(j + 1) % n]
        between = 0.0
        q = (i0 + 1) % len(cycle.legs)
        while q != i1:
            between += fixed_duration_s(system, cycle.legs[q])
            q = (q + 1) % len(cycle.legs)
        fixed_after.append(between)
    free = cycle.period_s - sum(fixed_after)
    if free <= 0:
        return
    levels = np.arange(1, n_split + 1)
    for t0 in np.linspace(0.0, phase_period_s, n_phase, endpoint=False):
        for w in itertools.product(levels, repeat=n - 1):
            last = n_split + 1 - sum(w) if n > 1 else 1
            if n > 1 and last < 1:
                continue
            ws = np.array([*w, last], dtype=float)
            durs = free * ws / ws.sum()
            x = np.empty(n)
            t = t0
            for j in range(n):
                x[j] = t
                t += durs[j] + fixed_after[j]
            yield x


@dataclass(frozen=True)
class Zero:
    """A converged date vector of a structure."""

    cycle: Cycle
    x: np.ndarray
    residual_kms: float


def solve_structure(
    system: CircularSystem,
    cycle: Cycle,
    *,
    phase_period_s: float,
    n_phase: int = 24,
    n_split: int = 8,
    n_refine: int = 10,
    min_sep_frac: float = 0.03,
    tol_kms: float = 1e-8,
) -> list[Zero]:
    """All distinct zeros reached from the best-residual seeds of the grid.

    Seeds are taken in residual order but a seed within ``min_sep_frac`` of the
    cycle period of one already refined (dates compared modulo the phase
    period) is skipped, so the ``n_refine`` refinements go to distinct basins.
    """
    scored: list[tuple[float, np.ndarray]] = []
    for x in _seed_grid(system, cycle, phase_period_s, n_phase, n_split):
        r = date_residual(system, cycle, x)
        if r is None:
            continue
        scored.append((float(np.sum(r * r)), x))
    scored.sort(key=lambda s: s[0])
    chosen: list[np.ndarray] = []
    sep = min_sep_frac * cycle.period_s
    for _, x0 in scored:
        if len(chosen) >= n_refine:
            break
        if any(_sep(x0, c, phase_period_s) < sep for c in chosen):
            continue
        chosen.append(x0)
    zeros: list[Zero] = []
    for x0 in chosen:
        sol: Solution = correct_dates(system, cycle, x0, tol_kms=tol_kms, max_nfev=80)
        if not sol.converged:
            continue
        x = _normalise(sol.x, phase_period_s)
        if any(_same(x, z.x, phase_period_s) for z in zeros):
            continue
        zeros.append(Zero(cycle, x, sol.max_abs_residual_kms))
    return zeros


def _sep(x: np.ndarray, y: np.ndarray, period: float) -> float:
    d = (x - y + 0.5 * period) % period - 0.5 * period
    return float(np.max(np.abs(d)))


def _normalise(x: np.ndarray, period: float) -> np.ndarray:
    return np.asarray(x - math.floor(x[0] / period) * period, dtype=np.float64)


def _same(x: np.ndarray, y: np.ndarray, period: float, tol_s: float = 0.01 * DAY) -> bool:
    d = (x - y + 0.5 * period) % period - 0.5 * period
    return bool(np.max(np.abs(d - d[0])) < tol_s and abs(d[0]) < tol_s)


# ---------------------------------------------------------------------------
# Assessment
# ---------------------------------------------------------------------------


@dataclass
class Assessment:
    """Everything the results note reports for one zero."""

    zero: Zero
    flybys: list[Flyby] | None
    report: CycleReport | None
    report_hm_floor: CycleReport | None
    max_encounter_miss_km: float
    min_soi_km: float
    r_min_km: float
    r_max_km: float
    vinf_kms: dict[str, float] = field(default_factory=dict)

    @property
    def status(self) -> str:
        return "no-directions" if self.report is None else self.report.status


def _arc_extent(
    mu: float, r0: np.ndarray, v0: np.ndarray, dt: float, n_samples: int
) -> tuple[float, float]:
    """Smallest and largest distance from the centre along a two-body arc of ``dt`` seconds.

    Sampled, then refined: when an apse is passed (the radial velocity changes sign between
    samples), its exact radius ``a(1 -/+ e)`` is used."""
    energy = 0.5 * float(v0 @ v0) - mu / float(np.linalg.norm(r0))
    h = np.cross(r0, v0)
    ecc = math.sqrt(max(0.0, 1.0 + 2.0 * energy * float(h @ h) / mu**2))
    a = -mu / (2.0 * energy) if energy < 0.0 else math.inf
    prev_rdot = float(r0 @ v0)
    lo = hi = float(np.linalg.norm(r0))
    for t in np.linspace(0.0, dt, n_samples)[1:]:
        r, v = kepler_step(r0, v0, float(t), mu)
        rn = float(np.linalg.norm(r))
        lo, hi = min(lo, rn), max(hi, rn)
        rdot = float(r @ v)
        if prev_rdot < 0.0 <= rdot:
            lo = min(lo, a * (1.0 - ecc) if math.isfinite(a) else lo)
        elif prev_rdot > 0.0 >= rdot and math.isfinite(a):
            hi = max(hi, a * (1.0 + ecc))
        prev_rdot = rdot
    return lo, hi


def leg_extent(
    system: CircularSystem,
    cycle: Cycle,
    x: np.ndarray,
    n_samples: int = 400,
    flybys: list[Flyby] | None = None,
) -> tuple[float, float]:
    """Smallest and largest distance from the central body over every leg.

    Lambert legs always; fixed legs (full-rev, half-rev) only when ``flybys`` (from
    :func:`cycle_flybys`) supply their departure directions. Without them the extent covers
    the Lambert legs only (#943: GanEur#316's half-rev apoapsis is the cycler's maximum).
    """
    legs = eval_lambert_legs(system, cycle, x)
    if legs is None:
        return math.nan, math.nan
    lo, hi = math.inf, 0.0
    for ev, leg_i in zip(legs, cycle.lambert_index, strict=True):
        leg = cycle.legs[leg_i]
        assert isinstance(leg, LambertLeg)
        r0, w0 = system.state(leg.frm, ev.t_dep)
        a_lo, a_hi = _arc_extent(system.mu, r0, w0 + ev.vinf_dep, ev.t_arr - ev.t_dep, n_samples)
        lo, hi = min(lo, a_lo), max(hi, a_hi)
    if flybys:
        for blk in _blocks(system, cycle, legs):
            for leg_f, t0 in blk.fixed:
                fb = next(f for f in flybys if f.body == blk.body and abs(f.t_s - t0) < 1.0)
                r0, w0 = system.state(blk.body, t0)
                a_lo, a_hi = _arc_extent(
                    system.mu, r0, w0 + fb.vinf_out, fixed_duration_s(system, leg_f), n_samples
                )
                lo, hi = min(lo, a_lo), max(hi, a_hi)
    return lo, hi


def assess(system: CircularSystem, zero: Zero, *, hm_floor_radii: float = 1.1) -> Assessment:
    """Flybys, the gate at the project floor and at H&M's 1.1-radius floor, and checks."""
    fl = cycle_flybys(system, zero.cycle, zero.x)
    rep = gate_cycle(system, fl) if fl else None
    rep_hm = None
    if fl:
        floors = {
            f.body: (hm_floor_radii - 1.0) * system.body(f.body).radius_km
            for f in fl
            if not system.body(f.body).massless
        }
        # one floor per body: gate each body's flybys at its own H&M floor
        rep_hm = _gate_with_body_floors(system, fl, floors)
    miss = encounter_self_consistency(system, zero.cycle, zero.x)
    bodies = {leg.frm for leg in zero.cycle.legs if isinstance(leg, LambertLeg)}
    soi = min(sphere_of_influence_km(system, b) for b in bodies)
    lo, hi = leg_extent(system, zero.cycle, zero.x, flybys=fl)
    vinf: dict[str, float] = {}
    if fl:
        for f in fl:
            vinf.setdefault(f.body, f.vinf_kms)
    return Assessment(zero, fl, rep, rep_hm, miss, soi, lo, hi, vinf)


def _gate_with_body_floors(
    system: CircularSystem, fl: list[Flyby], floors: dict[str, float]
) -> CycleReport:
    reports = [gate_cycle(system, [f], alt_floor_km=floors[f.body]) for f in fl if f.body in floors]
    encs = tuple(e for r in reports for e in r.gate.encounters)
    from cyclerfinder.verify.turn_gate import TurnGateReport

    gate = TurnGateReport(
        encounters=encs, mag_tol_kms=reports[0].gate.mag_tol_kms if reports else 1e-3
    )
    return CycleReport(tuple(fl), gate, any(r.near_180 for r in reports))


def flyby_table(system: CircularSystem, a: Assessment) -> list[dict[str, float | str]]:
    """Per-flyby rows for the results note."""
    rows: list[dict[str, float | str]] = []
    if not a.flybys or a.report is None:
        return rows
    gate_iter = iter(a.report.gate.encounters)
    for f in a.flybys:
        b = system.body(f.body)
        row: dict[str, float | str] = {
            "body": f.body,
            "t_days": f.t_s / DAY,
            "vinf_kms": f.vinf_kms,
            "turn_deg": f.turn_deg,
        }
        if not b.massless:
            e = next(gate_iter)
            row |= {
                "available_deg": e.available_bend_deg,
                "ratio": e.ratio,
                "required_alt_km": e.required_alt_km,
                "rmin_radii": rmin_radii(b, f.vinf_kms, math.radians(f.turn_deg)),
                "status": e.status,
            }
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Massless target placement (Russell architecture)
# ---------------------------------------------------------------------------


def place_massless_target(
    system: CircularSystem,
    cycle: Cycle,
    x: np.ndarray,
    target: str,
    *,
    n_samples: int = 4000,
) -> list[tuple[CircularSystem, Cycle, np.ndarray]]:
    """For each crossing of the target's orbit by a Lambert leg: set the target's
    phase so it is there, split the leg at the crossing (both sub-legs re-solved
    as Lambert arcs, the branch chosen to reproduce the original conic), and
    return ``(system, cycle, seed)`` for the corrector."""
    legs = eval_lambert_legs(system, cycle, x)
    if legs is None:
        return []
    a_t, p_t, _ = system.bodies[target]
    out = []
    for j, (ev, leg_i) in enumerate(zip(legs, cycle.lambert_index, strict=True)):
        leg = cycle.legs[leg_i]
        assert isinstance(leg, LambertLeg)
        r0, w0 = system.state(leg.frm, ev.t_dep)
        v0 = w0 + ev.vinf_dep
        ts = np.linspace(0.0, ev.t_arr - ev.t_dep, n_samples)
        rr = np.array([np.linalg.norm(kepler_step(r0, v0, float(t), system.mu)[0]) for t in ts])
        for i in range(len(ts) - 1):
            if (rr[i] - a_t) * (rr[i + 1] - a_t) >= 0:
                continue
            lo, hi = float(ts[i]), float(ts[i + 1])
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                rm = float(np.linalg.norm(kepler_step(r0, v0, mid, system.mu)[0]))
                if (rr[i] - a_t) * (rm - a_t) < 0:
                    hi = mid
                else:
                    lo = mid
            tc = 0.5 * (lo + hi)
            rc, vc = kepler_step(r0, v0, tc, system.mu)
            t_abs = ev.t_dep + tc
            theta0 = math.atan2(rc[1], rc[0]) - 2.0 * math.pi * t_abs / p_t
            bodies = dict(system.bodies)
            bodies[target] = (a_t, p_t, theta0)
            sys2 = replace(system, bodies=bodies, massless=system.massless | {target}, _fb={})
            split = _split_leg(sys2, leg, ev.t_dep, t_abs, ev.t_arr, v0, vc)
            if split is None:
                continue
            l1, l2 = split
            new_legs = list(cycle.legs)
            new_legs[leg_i : leg_i + 1] = [replace(l1, to=target), replace(l2, frm=target)]
            new_cycle = Cycle(tuple(new_legs), cycle.period_s)
            seed = np.insert(np.asarray(x, dtype=float), j + 1, t_abs)
            out.append((sys2, new_cycle, seed))
    return out


def _split_leg(
    system: CircularSystem, leg: LambertLeg, t0: float, tc: float, t1: float, v0: Vec, vc: Vec
) -> tuple[LambertLeg, LambertLeg] | None:
    """Lambert labels of the two sub-arcs that reproduce the original conic."""
    r0, _ = system.state(leg.frm, t0)
    rc = kepler_step(r0, v0, tc - t0, system.mu)[0]
    r1, _ = system.state(leg.to, t1)
    best: list[LambertLeg] = []
    for ra, rb, dt, vref in ((r0, rc, tc - t0, v0), (rc, r1, t1 - tc, vc)):
        try:
            sols = lambert(ra, rb, dt, mu=system.mu, max_revs=max(leg.nrev, 1))
        except (LambertError, ValueError):
            return None
        s = min(sols, key=lambda s: float(np.linalg.norm(s.v1 - vref)))
        if float(np.linalg.norm(s.v1 - vref)) > 1e-6 * float(np.linalg.norm(vref)):
            return None
        best.append(LambertLeg(leg.frm, leg.to, s.n_revs, s.branch))
    return best[0], best[1]
