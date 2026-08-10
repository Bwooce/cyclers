"""#816 -- unequal-leg-time discrete asymmetric-closure enumeration at Uranus.

The `#792` scoping pass (``docs/notes/2026-08-10-792-scoping-vs-680.md``,
reproducible artifact ``scripts/check_792_manifold_closed_form.py``) proved in
closed form that the `#558`/`#563`/`#680` anchor-flyby-anchor V-infinity-
magnitude closure system, under the lineage's EQUAL-leg-time genome (one
``tof`` for both legs), has NO isolated asymmetric closures: matching |vinf|
at both flyby-pair radii is an invertible affine (planar-Tisserand) map of the
leg invariants (E, h), so both magnitudes matching forces the return leg to be
a CONGRUENT (rotated) copy of the outbound conic; with equal leg times
congruence forces the mirror arc, collapsing the whole system to the 1-D line
family ``beta == (n_a - n_b)*tof (mod 180 deg)`` -- and the one missing
true-periodicity EQUALITY (pattern repeat) collapses that continuum exactly
onto the already-catalogued symmetric `#563`/`#569` family.

`#816` (registered 2026-08-10, expectations LOW) is the ONE formulation in
this lineage that is both genuinely asymmetric-capable and structurally able
to have ISOLATED (discrete) roots: INDEPENDENT leg times ``tof0 != tof1`` +
full periodicity. Three unknowns ``(beta, tof0, tof1)`` against three real
equalities:

  * E-match and h-match (equivalently the two |vinf|-magnitude residuals, by
    the same Tisserand-affine argument -- rev-class/leg-time-independent), and
  * the pattern-repeat (resonance) condition
    ``(n_a - n_b)*(tof0 + tof1) == 0 (mod 360 deg)``,

so the system is generically DISCRETE. The repeat condition is used here to
ELIMINATE one unknown exactly: ``T = tof0 + tof1 = q * T_syn`` for a positive
integer resonance order ``q`` (``T_syn = 360 deg / |n_a - n_b|``), leaving a
square 2-D root-solve ``F(beta, tof0) = [r_mid, r_periodic] = 0`` per
(pair, direction, q, n_rev0, n_rev1) -- a finite, cheap, closed-form-informed
enumeration (deflated Newton seeded from grid local minima), NOT a new
adaptive grid search. On the equal-tof diagonal ``tof0 = tof1 = T/2`` the
closed form predicts roots exactly at the symmetric `#563` goldens
(``beta == q*180 (mod 360)``), which double as positive controls; fixing
``T`` transversally cuts `#680`'s along-manifold null direction, so cond(J)
at roots is expected FINITE (isolation restored) -- reported per root.

Known kill-risks stated up front (from the task registration):

  * the zero-apsidal-rotation branch is a trivial resonant Keplerian orbit
    (leg1 ballistically continues leg0's orbit; neither flyby does any work)
    -- classified out as ``trivial_ballistic_resonant``, not a discovery;
  * along the equal-tof manifold the `#792` scoping check measured required
    turn angles of 137-156 deg vs single-digit achievable bends at these
    small Uranian moons -- the unequal-tof discrete roots may ALL fail the
    same physical-bend wall. Gates applied per root, verbatim from the
    lineage: `#324` physical bend floor (every encounter's achievable bend
    >= DEFAULT_MIN_USEFUL_BEND_DEG), required-turn feasibility (required
    turn <= achievable bend at that encounter's own V_inf -- the `#565` /
    `#680` "necessary-not-sufficient" gap, closed here), and the independent
    DOP853 per-leg cross-check (< 1 km arrival residual).

A clean negative (all discrete roots trivial/symmetric/gate-failing) is the
EXPECTED, fully acceptable outcome per the task's own registration. No
catalogue writeback; no novelty claim -- any survivor is flagged for
literature_check + Opus/Fable adjudication, never self-adjudicated.

Box/conditionality (recorded for the empty-regions stamp): circular-coplanar
patched-conic model (the lineage discovery genome); moons = the `#563` census
scope (Ariel/Umbriel/Titania/Oberon; Miranda excluded there); per-leg n_rev in
[0, 3] (the `#558` spec box, one Lambert branch per n_rev chosen by the
lineage's own min-departure-|vinf| convention); leg durations in
[0.15, 3.6]*sqrt(P_a*P_b) (superset of `#680`'s extended [0.5, 3.6] box);
resonance orders q with T = q*T_syn inside the leg-duration box.

Run as (single shot, ~1.5 h serial; or chunked as below, pool-parallel)::

    uv run python scripts/scan_816_unequal_tof_discrete_roots.py
    # chunked (same code path; each writes data/found/.../roots_partial_*.json):
    uv run python scripts/scan_816_unequal_tof_discrete_roots.py \
        --directions Ariel-Umbriel,Umbriel-Ariel --out roots_partial_au.json
    ...
    uv run python scripts/scan_816_unequal_tof_discrete_roots.py --merge
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from scan_558_uranus_all_pairs_offset_sweep import (  # noqa: E402
    _leg_options,
    _signed_mean_motion,
    residual_at_point,
)

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES  # noqa: E402
from cyclerfinder.search.deflated_newton import enumerate_roots  # noqa: E402
from cyclerfinder.search.discovery_campaign import (  # noqa: E402
    DAY_S,
    _mean_motion_rad_day,
    _moon_state,
)
from cyclerfinder.search.five_tier_prioritizer import PatchedConicLeg  # noqa: E402
from cyclerfinder.search.physical_sanity import (  # noqa: E402
    DEFAULT_MIN_USEFUL_BEND_DEG,
    candidate_passes_physical_gate,
)
from cyclerfinder.search.saturn_uranus_campaign import dop853_cross_check_leg  # noqa: E402

PRIMARY = "Uranus"
MOONS: tuple[str, ...] = ("Ariel", "Umbriel", "Titania", "Oberon")  # #563 census scope
N_REV_MAX = 3  # per-leg n_rev in [0,3], the #558 spec box
TOF_LO_SCALE = 0.15  # per-leg duration floor, x sqrt(P_a*P_b)
TOF_HI_SCALE = 3.6  # per-leg duration ceiling (matches #680's extended box)
N_BETA = 180  # 2-deg beta grid for seeding
N_TOF0 = 200  # tof0 grid points per q-window for seeding
SEED_CAP_PER_COMBO = 150
SEED_RESIDUAL_MAX_KMS = 1.0
ROOT_TOL_KMS = 1e-10  # Newton convergence tol (matches #680)
REVALIDATE_TOL_KMS = 1e-9
TURN_TRIVIAL_DEG = 0.05  # below this at BOTH nodes => ballistic resonant orbit
SYMMETRIC_DTOF_DAYS = 1e-6
MIN_LEG_DAYS_FRAC = 0.05  # discard roots with a leg shorter than 0.05*sqrt(PaPb)

DATA_DIR = ROOT / "data" / "found" / "816_unequal_tof_asymmetric_roots"


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT
        ).strip()
    except Exception:
        return "unknown"


@dataclass(frozen=True)
class PairCtx:
    anchor: str
    flyby: str
    mu: float
    sma_a: float
    sma_b: float
    n_a: float  # signed mean motion, rad/day
    n_b: float
    s: float  # sqrt(P_a * P_b), days
    t_syn: float  # 360 deg / |n_a - n_b|, days


def pair_ctx(anchor: str, flyby: str) -> PairCtx:
    mu = PRIMARIES[PRIMARY]
    sat_a, sat_b = SATELLITES[anchor], SATELLITES[flyby]
    n_a_mag = _mean_motion_rad_day(mu, sat_a.sma_km)
    n_b_mag = _mean_motion_rad_day(mu, sat_b.sma_km)
    p_a, p_b = 2 * math.pi / n_a_mag, 2 * math.pi / n_b_mag
    n_a = _signed_mean_motion(sat_a, n_a_mag)
    n_b = _signed_mean_motion(sat_b, n_b_mag)
    return PairCtx(
        anchor=anchor,
        flyby=flyby,
        mu=mu,
        sma_a=sat_a.sma_km,
        sma_b=sat_b.sma_km,
        n_a=n_a,
        n_b=n_b,
        s=math.sqrt(p_a * p_b),
        t_syn=2 * math.pi / abs(n_a - n_b),
    )


def _states(
    ctx: PairCtx, beta_deg: float, tof0: float, tof1: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Moon states at the 3 encounters (anchor t=0, flyby t=tof0, anchor t=T)."""
    r0, v0 = _moon_state(0.0, ctx.n_a, 0.0, ctx.sma_a, ctx.mu)
    r1, v1 = _moon_state(math.radians(beta_deg), ctx.n_b, tof0, ctx.sma_b, ctx.mu)
    r2, v2 = _moon_state(0.0, ctx.n_a, tof0 + tof1, ctx.sma_a, ctx.mu)
    return r0, v0, r1, v1, r2, v2


def residual_vec_unequal(
    ctx: PairCtx, beta_deg: float, tof0: float, tof1: float, n0: int, n1: int
) -> NDArray[np.float64] | None:
    """Signed [r_mid, r_periodic] with INDEPENDENT leg durations.

    Identical residual definition to ``scan_558``'s (leg0 arrival vs leg1
    departure |vinf| at the flyby; leg0 departure vs leg1 arrival |vinf| at
    the anchor), with the single generalization tof0 != tof1.
    """
    if tof0 <= 1e-3 or tof1 <= 1e-3:
        return None
    r0, v0, r1, v1, r2, v2 = _states(ctx, beta_deg, tof0, tof1)
    leg0 = _leg_options(r0, v0, r1, v1, tof0 * DAY_S, ctx.mu, max(n0, 1))
    leg1 = _leg_options(r1, v1, r2, v2, tof1 * DAY_S, ctx.mu, max(n1, 1))
    if n0 not in leg0 or n1 not in leg1:
        return None
    opt0, opt1 = leg0[n0], leg1[n1]
    return np.array([opt0.vinf_in - opt1.vinf_out, opt0.vinf_out - opt1.vinf_in])


def make_residual_fn(ctx: PairCtx, t_total: float, n0: int, n1: int) -> Any:
    def f(z: NDArray[np.float64]) -> NDArray[np.float64] | None:
        beta, tof0 = float(z[0]), float(z[1])
        return residual_vec_unequal(ctx, beta, tof0, t_total - tof0, n0, n1)

    return f


def grid_stage(
    ctx: PairCtx, q: int
) -> tuple[dict[tuple[int, int], list[NDArray[np.float64]]], float, float] | None:
    """Coarse (beta, tof0) grid; per-(n0,n1) local-minima seeds for Newton."""
    t_total = q * ctx.t_syn
    lo = max(TOF_LO_SCALE * ctx.s, t_total - TOF_HI_SCALE * ctx.s)
    hi = min(TOF_HI_SCALE * ctx.s, t_total - TOF_LO_SCALE * ctx.s)
    if hi <= lo:
        return None
    betas = np.linspace(0.0, 360.0, N_BETA, endpoint=False)
    tofs = np.linspace(lo, hi, N_TOF0)
    grids: dict[tuple[int, int], NDArray[np.float64]] = {
        (n0, n1): np.full((N_BETA, N_TOF0), np.inf)
        for n0 in range(N_REV_MAX + 1)
        for n1 in range(N_REV_MAX + 1)
    }
    for i, beta in enumerate(betas):
        for j, tof0 in enumerate(tofs):
            tof1 = t_total - tof0
            r0, v0, r1, v1, r2, v2 = _states(ctx, float(beta), float(tof0), tof1)
            leg0 = _leg_options(r0, v0, r1, v1, tof0 * DAY_S, ctx.mu, N_REV_MAX)
            leg1 = _leg_options(r1, v1, r2, v2, tof1 * DAY_S, ctx.mu, N_REV_MAX)
            if not leg0 or not leg1:
                continue
            for n0, opt0 in leg0.items():
                for n1, opt1 in leg1.items():
                    grids[(n0, n1)][i, j] = math.hypot(
                        opt0.vinf_in - opt1.vinf_out, opt0.vinf_out - opt1.vinf_in
                    )
    seeds: dict[tuple[int, int], list[NDArray[np.float64]]] = {}
    for combo, g in grids.items():
        finite = np.isfinite(g)
        if not finite.any():
            continue
        cand: list[tuple[float, float, float]] = []
        for i in range(N_BETA):
            im, ip = (i - 1) % N_BETA, (i + 1) % N_BETA  # beta wraps
            for j in range(N_TOF0):
                v = g[i, j]
                if not math.isfinite(v) or v > SEED_RESIDUAL_MAX_KMS:
                    continue
                neigh = [g[im, j], g[ip, j]]
                if j > 0:
                    neigh.append(g[i, j - 1])
                if j < N_TOF0 - 1:
                    neigh.append(g[i, j + 1])
                if all(v <= nb for nb in neigh):
                    cand.append((v, float(betas[i]), float(tofs[j])))
        cand.sort()
        picked = [np.array([b, t]) for _v, b, t in cand[:SEED_CAP_PER_COMBO]]
        # Always seed the closed-form-predicted symmetric diagonal points
        # (beta == q*180 mod 360, tof0 = T/2) -- positive-control anchors.
        for b_sym in ((q * 180.0) % 360.0, (q * 180.0 + 180.0) % 360.0):
            if lo <= t_total / 2 <= hi:
                picked.append(np.array([b_sym, t_total / 2]))
        if picked:
            seeds[combo] = picked
    return seeds, lo, hi


def _angle_deg(u: np.ndarray, w: np.ndarray) -> float:
    c = float(np.dot(u, w) / (np.linalg.norm(u) * np.linalg.norm(w)))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def _energy_h(r: np.ndarray, v: np.ndarray, mu: float) -> tuple[float, float]:
    energy = 0.5 * float(np.dot(v, v)) - mu / float(np.linalg.norm(r))
    return energy, float(np.cross(r, v)[2])


def _cond_j(ctx: PairCtx, t_total: float, n0: int, n1: int, beta: float, tof0: float) -> float:
    f = make_residual_fn(ctx, t_total, n0, n1)
    f0 = f(np.array([beta, tof0]))
    if f0 is None:
        return math.inf
    jac = np.zeros((2, 2))
    for i, h in ((0, 1e-6), (1, 1e-7)):
        zp = np.array([beta, tof0])
        zp[i] += h
        fp = f(zp)
        if fp is None:
            return math.inf
        jac[:, i] = (fp - f0) / h
    return float(np.linalg.cond(jac))


def classify_root(dtof_days: float, turn_b_deg: float, turn_a_deg: float) -> str:
    """Pure classification from the stored per-root turn/tof numbers.

    ``anchor_passive`` / ``flyby_passive`` are the two one-node-working
    generalizations of the registration's pre-declared trivial branch (a node
    whose required turn is ~0 is not doing flyby work: the trajectory passes
    it ballistically, so the object's propulsive skeleton is a SINGLE-moon
    repeated-flyby cycler with a passive crossing of the other moon -- the
    Russell-Strange passive-science-target architecture, outside this
    genome's two-sided bend-usefulness semantics; see the `#571` stamp's own
    caveat). Only ``asymmetric_candidate`` (BOTH nodes turning) can count as
    a genuine dual-flyby asymmetric closure of this lineage's genome.
    """
    if dtof_days < SYMMETRIC_DTOF_DAYS:
        return "symmetric_equal_tof"
    if turn_b_deg < TURN_TRIVIAL_DEG and turn_a_deg < TURN_TRIVIAL_DEG:
        return "trivial_ballistic_resonant"
    if turn_b_deg < TURN_TRIVIAL_DEG:
        return "flyby_passive"
    if turn_a_deg < TURN_TRIVIAL_DEG:
        return "anchor_passive"
    return "asymmetric_candidate"


def analyze_root(
    ctx: PairCtx, q: int, n0: int, n1: int, beta_deg: float, tof0: float
) -> dict[str, Any] | None:
    """Classify one converged root + run the full lineage gate battery."""
    t_total = q * ctx.t_syn
    tof1 = t_total - tof0
    res = residual_vec_unequal(ctx, beta_deg, tof0, tof1, n0, n1)
    if res is None:
        return None
    res_norm = float(np.linalg.norm(res))
    if res_norm > REVALIDATE_TOL_KMS:
        return None
    r0, v0, r1, v1, r2, v2 = _states(ctx, beta_deg, tof0, tof1)
    leg0 = _leg_options(r0, v0, r1, v1, tof0 * DAY_S, ctx.mu, max(n0, 1))
    leg1 = _leg_options(r1, v1, r2, v2, tof1 * DAY_S, ctx.mu, max(n1, 1))
    if n0 not in leg0 or n1 not in leg1:
        return None
    opt0, opt1 = leg0[n0], leg1[n1]

    vout_a0 = np.asarray(opt0.v1) - v0
    vin_b = np.asarray(opt0.v2) - v1
    vout_b = np.asarray(opt1.v1) - v1
    vin_a2 = np.asarray(opt1.v2) - v2

    e0, h0 = _energy_h(r0, np.asarray(opt0.v1), ctx.mu)
    e1, h1 = _energy_h(r1, np.asarray(opt1.v1), ctx.mu)

    turn_b = _angle_deg(vin_b, vout_b)
    theta = ctx.n_a * t_total  # anchor's rotation over one full cycle (rad)
    c, s = math.cos(-theta), math.sin(-theta)
    rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    turn_a = _angle_deg(rot @ vin_a2, vout_a0)

    # Lineage per-encounter V_inf extraction (max of in/out asymptote at each
    # node, with vin=[0, vin_b, vin_a2] / vout=[vout_a0, vout_b, 0] -- same as
    # scan_558.encounter_vinfs_kms).
    vinfs = (
        abs(opt0.vinf_out),
        max(abs(opt0.vinf_in), abs(opt1.vinf_out)),
        abs(opt1.vinf_in),
    )
    seq = (ctx.anchor, ctx.flyby, ctx.anchor)
    bend_pass, verdicts = candidate_passes_physical_gate(
        seq, vinfs, min_useful_bend_deg=DEFAULT_MIN_USEFUL_BEND_DEG
    )
    max_bends = [float(v.max_bend_deg) for v in verdicts]
    achievable_b = max_bends[1]
    achievable_a = min(max_bends[0], max_bends[2])
    turn_feasible = bool(turn_b <= achievable_b and turn_a <= achievable_a)

    klass = classify_root(abs(tof0 - tof1), turn_b, turn_a)

    # Independent DOP853 cross-check only where the verdict depends on it
    # (asymmetric/anchor-passive candidates + symmetric controls), O(0.1 s)/leg.
    dop: dict[str, Any] | None = None
    if klass in ("asymmetric_candidate", "anchor_passive", "flyby_passive", "symmetric_equal_tof"):
        km_m = 1000.0
        mu_m3_s2 = ctx.mu * km_m**3
        legs = [
            PatchedConicLeg(
                label_from=ctx.anchor,
                label_to=ctx.flyby,
                r1_m=r0 * km_m,
                v1_m_s=np.asarray(opt0.v1) * km_m,
                r2_m=r1 * km_m,
                v2_m_s=np.asarray(opt0.v2) * km_m,
                dt_s=tof0 * DAY_S,
                mu_m3_s2=mu_m3_s2,
            ),
            PatchedConicLeg(
                label_from=ctx.flyby,
                label_to=ctx.anchor,
                r1_m=r1 * km_m,
                v1_m_s=np.asarray(opt1.v1) * km_m,
                r2_m=r2 * km_m,
                v2_m_s=np.asarray(opt1.v2) * km_m,
                dt_s=tof1 * DAY_S,
                mu_m3_s2=mu_m3_s2,
            ),
        ]
        checks = [dop853_cross_check_leg(leg, rtol=1e-12, atol=1e-12) for leg in legs]
        max_dr = max(float(cc["dr_arrival_km"]) for cc in checks)
        dop = {"max_dr_arrival_km": max_dr, "passed": bool(max_dr < 1.0)}

    dop_pass = bool(dop["passed"]) if dop is not None else False
    physical = bool(bend_pass and turn_feasible and dop_pass)
    # A genuine dual-flyby asymmetric survivor must have BOTH nodes doing
    # useful turn work (class asymmetric_candidate) AND pass every physical
    # gate; passive-node roots that pass the physical gates are reported
    # separately (out-of-genome architecture, flagged for adjudication).
    survivor = bool(klass == "asymmetric_candidate" and physical)
    return {
        "anchor": ctx.anchor,
        "flyby": ctx.flyby,
        "q": q,
        "t_total_days": t_total,
        "t_syn_days": ctx.t_syn,
        "n_rev": [n0, n1],
        "beta_deg": beta_deg % 360.0,
        "tof0_days": tof0,
        "tof1_days": tof1,
        "residual_kms": res_norm,
        "cond_j": _cond_j(ctx, t_total, n0, n1, beta_deg, tof0),
        "class": klass,
        "energy_leg0": e0,
        "energy_leg1": e1,
        "h_leg0": h0,
        "h_leg1": h1,
        "congruence_dE": e1 - e0,
        "congruence_dh": h1 - h0,
        "required_turn_flyby_deg": turn_b,
        "required_turn_anchor_deg": turn_a,
        "vinf_per_encounter_kms": list(vinfs),
        "max_bend_deg_per_encounter": max_bends,
        "bend_gate_passed": bool(bend_pass),
        "turn_feasible": turn_feasible,
        "dop853": dop,
        "passes_physical_gates": physical,
        "all_gates_passed": survivor,
    }


def solve_q(ctx: PairCtx, q: int) -> list[dict[str, Any]]:
    """Enumerate discrete roots for one (anchor, flyby, q) resonance order."""
    staged = grid_stage(ctx, q)
    if staged is None:
        return []
    seeds, _lo, _hi = staged
    t_total = q * ctx.t_syn
    out: list[dict[str, Any]] = []
    for (n0, n1), seed_list in sorted(seeds.items()):
        f = make_residual_fn(ctx, t_total, n0, n1)
        roots = enumerate_roots(
            f,
            seed_list,
            tol=ROOT_TOL_KMS,
            max_iter=60,
            step_cap=np.array([20.0, 0.3 * ctx.s]),
            dedup_tol=1e-4,
        )
        seen: set[tuple[float, float]] = set()
        for r in roots:
            beta, tof0 = float(r[0]) % 360.0, float(r[1])
            tof1 = t_total - tof0
            if min(tof0, tof1) < MIN_LEG_DAYS_FRAC * ctx.s:
                continue
            if max(tof0, tof1) > TOF_HI_SCALE * ctx.s + 0.1 * ctx.s:
                continue  # outside the declared leg-duration box
            key = (round(beta, 3), round(tof0, 5))
            if key in seen:
                continue
            seen.add(key)
            rec = analyze_root(ctx, q, n0, n1, beta, tof0)
            if rec is not None:
                out.append(rec)
    return out


def q_range(ctx: PairCtx) -> range:
    q_min = max(1, math.ceil(2 * TOF_LO_SCALE * ctx.s / ctx.t_syn))
    q_max = math.floor(2 * TOF_HI_SCALE * ctx.s / ctx.t_syn)
    return range(q_min, q_max + 1)


def _solve_task(args: tuple[str, str, int]) -> tuple[str, str, int, list[dict[str, Any]]]:
    """Pool worker: one (anchor, flyby, q) solve (module-level for pickling)."""
    anchor, flyby, q = args
    recs = solve_q(pair_ctx(anchor, flyby), q)
    return anchor, flyby, q, recs


def positive_controls() -> None:
    """PC1: equal-tof reduction reproduces scan_558's residual exactly.
    PC2: the pipeline recovers the catalogued #563/#569 symmetric golden."""
    print("[PC1] equal-tof faithfulness vs scan_558.residual_at_point:")
    ctx = pair_ctx("Umbriel", "Oberon")
    tof = 2.0 * ctx.s
    pt = residual_at_point("Umbriel", "Oberon", rel_offset_deg=180.0, tof_scale=2.0, n_rev=(1, 1))
    assert pt is not None
    mine = residual_vec_unequal(ctx, 180.0, tof, tof, 1, 1)
    assert mine is not None
    diff = abs(float(np.max(np.abs(mine))) - float(pt["residual_kms"]))
    print(
        f"    #312 point: scan_558 residual={pt['residual_kms']:.15f}, "
        f"unequal-tof(tof0==tof1) worst={float(np.max(np.abs(mine))):.15f}, |diff|={diff:.2e}"
    )
    assert diff < 1e-12, "PC1 FAILED: unequal-tof residual is not a faithful generalization"

    print("[PC2] recover the catalogued Ariel-Umbriel (0,0) symmetric golden (q=1):")
    ctx = pair_ctx("Ariel", "Umbriel")
    golden_tof = 3.216088179066208  # data/enumerate_563_symmetric_closures.jsonl
    roots = solve_q(ctx, 1)
    hits = [
        r
        for r in roots
        if r["n_rev"] == [0, 0]
        and abs(r["tof0_days"] - golden_tof) < 1e-5
        and abs(r["tof1_days"] - golden_tof) < 1e-5
        and (r["beta_deg"] < 1.0 or r["beta_deg"] > 359.0 or abs(r["beta_deg"] - 180.0) < 1.0)
    ]
    assert hits, "PC2 FAILED: symmetric golden not recovered as a root of the unequal-tof system"
    h = hits[0]
    print(
        f"    recovered: beta={h['beta_deg']:.6f} tof0={h['tof0_days']:.9f} "
        f"tof1={h['tof1_days']:.9f} res={h['residual_kms']:.2e} cond(J)={h['cond_j']:.2e}"
    )
    if h["cond_j"] > 1e7:
        print("    WARNING: cond(J) large at the symmetric control -- isolation suspect")


def summarize_and_write(
    all_roots: list[dict[str, Any]], out_path: Path, t_start: float, note: str
) -> None:
    # Retrofit-reclassify from the stored per-root numbers (pure function of
    # the record), so chunk files written by any code revision merge into one
    # consistent classification; recompute the missing DOP853 check for any
    # bend+turn-passing record that predates its class being DOP-checked.
    for idx, r in enumerate(all_roots):
        r["class"] = classify_root(
            abs(r["tof0_days"] - r["tof1_days"]),
            r["required_turn_flyby_deg"],
            r["required_turn_anchor_deg"],
        )
        dop = r.get("dop853")
        if dop is None and r["bend_gate_passed"] and r["turn_feasible"]:
            fresh = analyze_root(
                pair_ctx(r["anchor"], r["flyby"]),
                r["q"],
                r["n_rev"][0],
                r["n_rev"][1],
                r["beta_deg"],
                r["tof0_days"],
            )
            if fresh is not None:
                all_roots[idx] = r = fresh
                dop = r.get("dop853")
        physical = bool(
            r["bend_gate_passed"] and r["turn_feasible"] and dop is not None and dop["passed"]
        )
        r["passes_physical_gates"] = physical
        r["all_gates_passed"] = bool(r["class"] == "asymmetric_candidate" and physical)

    by_class: dict[str, int] = {}
    for r in all_roots:
        by_class[r["class"]] = by_class.get(r["class"], 0) + 1
    survivors = [r for r in all_roots if r["all_gates_passed"]]
    passive_phys = [
        r
        for r in all_roots
        if r["class"] in ("anchor_passive", "flyby_passive") and r["passes_physical_gates"]
    ]
    asym = [r for r in all_roots if r["class"] == "asymmetric_candidate"]
    max_cond = max((r["cond_j"] for r in all_roots if math.isfinite(r["cond_j"])), default=0.0)
    max_congr = max(
        (max(abs(r["congruence_dE"]), abs(r["congruence_dh"])) for r in all_roots), default=0.0
    )

    print("\n" + "=" * 78)
    print(f"TOTAL validated discrete roots: {len(all_roots)}  classes: {by_class}")
    print(f"max cond(J) over all roots: {max_cond:.2e} (isolation diagnostic)")
    print(f"max |dE|,|dh| congruence defect over all roots: {max_congr:.2e}")
    print(f"asymmetric candidates: {len(asym)}; ALL-GATES survivors: {len(survivors)}")
    if asym:
        print("\nasymmetric candidates (gate outcomes):")
        for r in sorted(asym, key=lambda x: -min(x["max_bend_deg_per_encounter"]))[:40]:
            print(
                f"  {r['anchor']}-{r['flyby']} q={r['q']} nrev={r['n_rev']} "
                f"beta={r['beta_deg']:8.3f} tof0={r['tof0_days']:8.4f} "
                f"tof1={r['tof1_days']:8.4f} req_turn(B,A)=({r['required_turn_flyby_deg']:7.2f},"
                f"{r['required_turn_anchor_deg']:7.2f}) "
                f"achievable={[f'{b:.2f}' for b in r['max_bend_deg_per_encounter']]} "
                f"bend={r['bend_gate_passed']} turnfeas={r['turn_feasible']} "
                f"gates={r['all_gates_passed']}"
            )
    if survivors:
        print("\nSURVIVORS (NOT a novelty claim -- literature_check + adjudication required):")
        for r in survivors:
            print(f"  {json.dumps(r, default=float)}")
    if passive_phys:
        print(
            f"\nPASSIVE-NODE roots passing ALL PHYSICAL gates: {len(passive_phys)} "
            "(one node requires ~zero turn -> the object's propulsive skeleton is a"
        )
        print(
            "SINGLE-moon repeated-flyby cycler with a passive crossing of the other moon;"
            " out-of-genome architecture -- FLAGGED FOR ADJUDICATION, not a survivor):"
        )
        for r in passive_phys:
            print(
                f"  {r['anchor']}-{r['flyby']} q={r['q']} nrev={r['n_rev']} "
                f"beta={r['beta_deg']:.6f} tof0={r['tof0_days']:.6f} tof1={r['tof1_days']:.6f} "
                f"req_turn(B,A)=({r['required_turn_flyby_deg']:.9f},"
                f"{r['required_turn_anchor_deg']:.9f}) deg class={r['class']}"
            )

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        json.dump(
            {
                "task": 816,
                "git_sha": _git_sha(),
                "note": note,
                "primary": PRIMARY,
                "moons": list(MOONS),
                "n_rev_max": N_REV_MAX,
                "tof_scale_box": [TOF_LO_SCALE, TOF_HI_SCALE],
                "grid": {"n_beta": N_BETA, "n_tof0": N_TOF0},
                "gates": {
                    "bend_floor_deg": DEFAULT_MIN_USEFUL_BEND_DEG,
                    "root_tol_kms": ROOT_TOL_KMS,
                    "dop853_dr_km": 1.0,
                },
                "elapsed_s": time.time() - t_start,
                "n_roots": len(all_roots),
                "classes": by_class,
                "n_survivors": len(survivors),
                "n_passive_node_physical_gate_passers": len(passive_phys),
                "roots": all_roots,
            },
            fh,
            indent=1,
            default=float,
        )
    print(f"\nwrote {out_path} ({len(all_roots)} roots, {time.time() - t_start:.0f}s)")


def merge_partials() -> int:
    """Combine roots_partial_*.json chunks into the final roots.json + summary."""
    t_start = time.time()
    partials = sorted(DATA_DIR.glob("roots_partial_*.json"))
    if not partials:
        print("no roots_partial_*.json found to merge")
        return 1
    all_roots: list[dict[str, Any]] = []
    covered: set[tuple[str, str]] = set()
    elapsed = 0.0
    for p in partials:
        with p.open() as fh:
            d = json.load(fh)
        all_roots.extend(d["roots"])
        elapsed += float(d.get("elapsed_s", 0.0))
        for r in d["roots"]:
            covered.add((r["anchor"], r["flyby"]))
        print(f"  merged {p.name}: {d['n_roots']} roots ({d.get('note', '')})")
    expected = {(a, b) for a in MOONS for b in MOONS if a != b}
    missing = expected - covered
    if missing:
        print(f"  WARNING: directions with NO roots found in any partial: {sorted(missing)}")
        print("  (check the per-chunk logs -- a direction may legitimately have 0 roots)")
    summarize_and_write(
        all_roots,
        DATA_DIR / "roots.json",
        t_start - elapsed,  # so elapsed_s reports total compute across chunks
        f"merged from {len(partials)} chunk files",
    )
    return 0


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--directions",
        default="",
        help="comma-separated Anchor-Flyby list (default: all 12 ordered pairs)",
    )
    parser.add_argument("--out", default="roots.json", help="output filename under DATA_DIR")
    parser.add_argument("--merge", action="store_true", help="merge roots_partial_*.json chunks")
    parser.add_argument("--skip-pc", action="store_true", help="skip positive controls")
    parser.add_argument("--pool", type=int, default=6, help="worker processes over (dir, q) tasks")
    args = parser.parse_args()

    if args.merge:
        return merge_partials()

    t_start = time.time()
    print("#816 unequal-leg-time (tof0 != tof1) discrete asymmetric-closure enumeration")
    print(f"primary={PRIMARY} moons={MOONS} n_rev box=[0,{N_REV_MAX}] ")
    print(
        f"leg-duration box=[{TOF_LO_SCALE},{TOF_HI_SCALE}]*sqrt(PaPb); "
        f"pattern-repeat T=q*T_syn eliminates tof1; gates: #324 bend floor "
        f">={DEFAULT_MIN_USEFUL_BEND_DEG} deg + required-turn feasibility + DOP853 <1 km\n"
    )
    if not args.skip_pc:
        positive_controls()

    if args.directions:
        directions = []
        for tok in args.directions.split(","):
            a, b = tok.strip().split("-")
            assert a in MOONS and b in MOONS and a != b, f"bad direction {tok!r}"
            directions.append((a, b))
    else:
        directions = [(a, b) for a in MOONS for b in MOONS if a != b]

    tasks: list[tuple[str, str, int]] = []
    for anchor, flyby in directions:
        ctx = pair_ctx(anchor, flyby)
        for q in q_range(ctx):
            tasks.append((anchor, flyby, q))
    print(
        f"[{time.strftime('%H:%M:%S')}] {len(directions)} directions -> {len(tasks)} (dir,q) "
        f"tasks, pool={args.pool}",
        flush=True,
    )

    all_roots: list[dict[str, Any]] = []
    if args.pool > 1:
        import multiprocessing as mp

        with mp.get_context("spawn").Pool(processes=args.pool) as pool:
            for anchor, flyby, q, recs in pool.imap_unordered(_solve_task, tasks):
                all_roots.extend(recs)
                print(
                    f"    [{time.strftime('%H:%M:%S')}] {anchor}-{flyby} q={q}: "
                    f"{len(recs)} validated roots",
                    flush=True,
                )
    else:
        for task in tasks:
            anchor, flyby, q, recs = _solve_task(task)
            all_roots.extend(recs)
            print(
                f"    [{time.strftime('%H:%M:%S')}] {anchor}-{flyby} q={q}: "
                f"{len(recs)} validated roots",
                flush=True,
            )

    summarize_and_write(
        all_roots,
        DATA_DIR / args.out,
        t_start,
        f"directions={','.join(f'{a}-{b}' for a, b in directions)}",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
