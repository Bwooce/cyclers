"""CR3BP stable-family continuation search -- the project's first novel-discovery campaign.

Spec: ``docs/superpowers/specs/2026-06-12-cr3bp-family-continuation-search-design.md``
(Phase 2). Seeds natural-parameter (Jacobi) continuation from sourced/published
members, walks each family BOTH directions in the Jacobi constant, and for every
gauntlet-passing member runs an INDEPENDENT inertial n-body cross-check (REBOUND /
IAS15 -- a different code path AND frame than the rotating-frame DOP853 corrector,
the spec's false-consensus independence gate). Each member is then classified:

  REPRODUCTION : within dedup tolerance of a published Ross / Arenstorf member --
                 a SOURCED cross-check (validates the engine), NOT a discovery.
  NOVEL-SILVER : new (not in our sourced set), passes the full gauntlet AND the
                 inertial cross-check -> routed to the review queue as SILVER.
                 NEVER auto-promoted to the catalogue.
  EMPTY        : a Jacobi range that yields no new member under this method+version
                 -> the method-versioned negative registry (data/empty_regions.jsonl).

Seeds:
  (a) the 5 Ross & Roberts-Tsoukkas 2025 stable EM (k1,k2) cycler families;
  (b) the Arenstorf 1963 figure-eight (CR3BP test-problem mu);
  (c) the Tier-2 Saturnian midsize-moon planar Lyapunov seeds (#182).

Parallelised across (seed x direction) work units on the 16-core box (the #182
process-pool pattern). NO silent caps: every Jacobi range, step, and seed is
logged. SILVER-only, NO catalogue writeback. Planar CR3BP (PCR3BP) only.

NOVELTY SCOPE CAVEAT: dedup is against our SOURCED set only (Ross's 5 published
members, the Arenstorf orbit, and the existing data/cr3bp_silver.jsonl rows). The
JPL Three-Body Periodic Orbit Catalog is NOT yet acquired (#116); any "novel"
member is therefore "not in our sourced set", not "new to the literature". This
is stated in every output.

Usage:
    uv run python scripts/cr3bp_family_search.py
    uv run python scripts/cr3bp_family_search.py --no-writeback
    uv run python scripts/cr3bp_family_search.py --workers 8
"""

from __future__ import annotations

import argparse
import contextlib
import json
import math
import subprocess
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_continuation as cc
import cyclerfinder.search.cr3bp_periodic as cp

METHOD_TAG = "cr3bp-jacobi-continuation-v1"
INERTIAL_TAG = "cr3bp-inertial-rebound-ias15-v1"

_GIT_SHA: str = ""
with contextlib.suppress(Exception):
    _GIT_SHA = (
        subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], stderr=subprocess.DEVNULL)
        .decode()
        .strip()
    )

REPO_ROOT = Path(__file__).resolve().parent.parent
REVIEW_QUEUE_PATH = REPO_ROOT / "data" / "cr3bp_continuation_review_queue.jsonl"
EMPTY_REGIONS_PATH = REPO_ROOT / "data" / "empty_regions.jsonl"
SILVER_QUEUE_PATH = REPO_ROOT / "data" / "cr3bp_silver.jsonl"
NOTE_PATH = REPO_ROOT / "docs" / "notes" / "2026-06-12-cr3bp-continuation-results.md"
RUNLOG_PATH = REPO_ROOT / "data" / "runs" / "cr3bp-continuation"

# Ross & Roberts-Tsoukkas 2025, AAS 25-621, p. 3.
ROSS_MU = 1.2150584270572e-2
ROSS_L_KM = 384400.0
ROSS_T_S = 375699.8

# Continuation walk parameters (no silent caps -- logged in the EMPTY registry).
DEFAULT_DC = 1e-4
DEFAULT_N_STEPS = 25
# Jacobi window: below C(L1)=3.18834 (the energy ceiling for these prograde
# cyclers, Ross p. 13) and above a generous floor.
JACOBI_MIN = 3.00
JACOBI_MAX = 3.18834
# Local half-window around each seed's own C (so a seed outside the global EM band,
# e.g. Arenstorf at C~2.86, still gets a symmetric Jacobi range to walk).
JACOBI_HALFWIDTH = DEFAULT_DC * DEFAULT_N_STEPS  # = n_steps * dC, one full walk each way

# A corrected seed must reproduce its published period to within this fraction or
# it is off-family (wrong crossing index / sign) and is abandoned.
SEED_PERIOD_TOL_FRAC = 0.02

# Dedup tolerances for the REPRODUCTION test (against the sourced set).
DEDUP_X0_TOL = 5e-4  # nd: matches the corrector's window-width resolution
DEDUP_C_TOL = 5e-4
DEDUP_T_TOL = 5e-3

# Inertial cross-check (REBOUND/IAS15) parameters.
N_PERIODS_INERTIAL = 5
N_SAMP_PER_PERIOD = 128
IAS15_EPSILON = 1e-9


# ---------------------------------------------------------------------------
# Seeds
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class FamilySeed:
    """A continuation seed: a sourced/published member + its branch identity."""

    label: str
    kind: str  # "ross" | "arenstorf" | "saturnian-lyapunov"
    mu: float
    l_km: float
    t_s: float
    primary: str
    secondary: str
    x0_seed: float
    jacobi: float
    period_guess: float
    ydot0_sign: float
    half_crossings: int
    note: str


# (a) Ross 5 families (Table 3 C^stable / T^stable; x0 derived seeds from #212).
ROSS_SEEDS = [
    FamilySeed(
        "ross-(1,1)",
        "ross",
        ROSS_MU,
        ROSS_L_KM,
        ROSS_T_S,
        "Earth",
        "Moon",
        -0.7682140805,
        3.151175879508174,
        10.29206921007976,
        -1.0,
        3,
        "Ross Table 3 (1,1) nu=0 midpoint",
    ),
    FamilySeed(
        "ross-(2,1)",
        "ross",
        ROSS_MU,
        ROSS_L_KM,
        ROSS_T_S,
        "Earth",
        "Moon",
        0.7237335857,
        3.129389531088256,
        19.44043166795154,
        1.0,
        4,
        "Ross Table 3 (2,1) nu=0 midpoint; razor-thin window",
    ),
    FamilySeed(
        "ross-(3,1)",
        "ross",
        ROSS_MU,
        ROSS_L_KM,
        ROSS_T_S,
        "Earth",
        "Moon",
        -0.3209891696,
        3.161784147013429,
        14.78849241668140,
        -1.0,
        3,
        "Ross Table 3 (3,1) nu=0 midpoint",
    ),
    FamilySeed(
        "ross-(3,2)",
        "ross",
        ROSS_MU,
        ROSS_L_KM,
        ROSS_T_S,
        "Earth",
        "Moon",
        -0.3210000000,
        3.182762663084288,
        17.90058010350006,
        -1.0,
        6,
        "Ross Table 3 (3,2) nu=0 midpoint; half-period = 6th x-crossing",
    ),
    FamilySeed(
        "ross-(3,3)",
        "ross",
        ROSS_MU,
        ROSS_L_KM,
        ROSS_T_S,
        "Earth",
        "Moon",
        -0.3217380626,
        3.177224018696528,
        18.14546057589189,
        -1.0,
        5,
        "Ross Table 3 (3,3) nu=0 midpoint; wide window",
    ),
]

# (b) Arenstorf 1963 figure-eight (CR3BP test-problem mu; perpendicular x-crossing IC).
ARENSTORF_MU = 0.012277471
ARENSTORF_SEED = FamilySeed(
    "arenstorf-1963",
    "arenstorf",
    ARENSTORF_MU,
    ROSS_L_KM,
    ROSS_T_S,
    "Earth",
    "Moon",
    0.994,
    cr3bp.jacobi_constant(np.array([0.994, 0.0, 0.0, 0.0, -2.0015851063790825, 0.0]), ARENSTORF_MU),
    17.0652165601579625,
    -1.0,
    1,
    "Arenstorf 1963 / Hairer B5; UNSTABLE figure-eight",
)


def _saturnian_lyapunov_seeds() -> list[FamilySeed]:
    """Build symmetric-orbit continuation seeds from the Tier-2 SILVER Lyapunovs.

    The #182 SILVER rows are perpendicular x-axis-crossing planar Lyapunov orbits
    (state0_nd = (x0, ~0, 0, ~0, vy0, 0)). We re-seed one converged member per
    (Saturn, moon) pair as a symmetric-continuation seed (the half-period crossing
    is the 1st x-axis crossing for these single-loop Lyapunovs).
    """
    if not SILVER_QUEUE_PATH.exists():
        return []
    seen: set[str] = set()
    seeds: list[FamilySeed] = []
    for ln in SILVER_QUEUE_PATH.read_text(encoding="utf-8").splitlines():
        if not ln.strip():
            continue
        rec = json.loads(ln)
        sec = str(rec["secondary"])
        if sec in seen:
            continue
        seen.add(sec)
        system = cr3bp.cr3bp_system(str(rec["primary"]), sec)
        s0 = np.asarray(rec["state0_nd"], dtype=np.float64)
        seeds.append(
            FamilySeed(
                label=f"saturn-{sec.lower()}-lyapunov",
                kind="saturnian-lyapunov",
                mu=system.mu,
                l_km=system.l_km,
                t_s=system.t_s,
                primary=str(rec["primary"]),
                secondary=sec,
                x0_seed=float(s0[0]),
                jacobi=float(rec["jacobi_constant"]),
                period_guess=float(rec["period_nd"]),
                ydot0_sign=math.copysign(1.0, float(s0[4])),
                half_crossings=1,
                note=f"Tier-2 #182 SILVER Lyapunov seed ({rec['libration_point']})",
            )
        )
    return seeds


def all_seeds() -> list[FamilySeed]:
    return [*ROSS_SEEDS, ARENSTORF_SEED, *_saturnian_lyapunov_seeds()]


# ---------------------------------------------------------------------------
# Independent inertial cross-check (REBOUND / IAS15) -- false-consensus gate
# ---------------------------------------------------------------------------
def _rot_to_inertial_t0(
    state_nd: NDArray[np.float64], l_km: float, n_rad_s: float
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    r_km = state_nd[:3] * l_km
    v_rot_km_s = state_nd[3:] * (l_km * n_rad_s)
    omega_cross_r = n_rad_s * np.array([-r_km[1], r_km[0], 0.0])
    return r_km, v_rot_km_s + omega_cross_r


def _inertial_to_rot_nd(
    r_km: NDArray[np.float64],
    v_km_s: NDArray[np.float64],
    theta: float,
    l_km: float,
    n_rad_s: float,
) -> NDArray[np.float64]:
    c, s = math.cos(theta), math.sin(theta)
    rot = np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])
    r_rot = rot @ r_km
    v_rot = rot @ v_km_s - n_rad_s * np.array([-r_rot[1], r_rot[0], 0.0])
    return np.concatenate([r_rot / l_km, v_rot / (l_km * n_rad_s)])


@dataclass(frozen=True)
class InertialVerdict:
    """Inertial-crosscheck metrics + verdict (graded against linear instability)."""

    verdict: str  # "PASS" | "CHECK-FAILED" | "INCONCLUSIVE"
    delta_per_period_nd: tuple[float, ...]
    delta1_nd: float
    abs_lambda_per_period: float  # |monodromy eigenvalue| from nu (theory)
    t_bound_nd: float
    t_dep_pred_nd: float
    jacobi_drift_bounded: float
    noise_floor_nd: float
    diverged: bool


def inertial_crosscheck(
    member: cc.BranchMember,
    seed: FamilySeed,
    *,
    n_periods: int = N_PERIODS_INERTIAL,
) -> InertialVerdict:
    """Re-propagate a member in the inertial REBOUND/IAS15 harness and grade it.

    Like-for-like CR3BP-consistency (NOT a real-ephemeris claim): the two primaries
    sit on the exact circular two-body rail of the CR3BP idealisation; the
    spacecraft is a massless test particle; IAS15 integrates the inertial n-body
    problem; the trajectory is back-transformed to the rotating frame.

    Grading is against LINEAR-INSTABILITY THEORY (the #182 method; band never
    loosened): the per-period amplification is the monodromy eigenvalue magnitude
    ``|lambda| = |nu + sqrt(nu^2 - 1)|`` from the member's Barden nu. For a STABLE
    member (|nu|<1) |lambda|=1 -- the orbit must stay BOUNDED over n_periods (R1+R2).
    For an UNSTABLE member the predicted departure time
    ``t_dep_pred = T + ln(D/delta1)/ln|lambda|`` bounds when it must leave; observed
    departure EARLIER than predicted contradicts the orbit (CHECK-FAILED).
    """
    import rebound

    mu, l_km, t_s = seed.mu, seed.l_km, seed.t_s
    n_rad_s = 1.0 / t_s
    # G=1, masses carry GM (km^3/s^2). For a generic mu we set the system GM so the
    # mean motion n = sqrt(GM_pair / l^3) matches 1/t_s exactly (self-consistent rail).
    gm_pair = (n_rad_s**2) * l_km**3
    gm2 = mu * gm_pair
    gm1 = (1.0 - mu) * gm_pair

    period_nd = member.period
    amp = max(member.amplitude_nd, 1e-9)
    state0_nd = member.state0
    j0 = cr3bp.jacobi_constant(state0_nd, mu)

    sc_r0, sc_v0 = _rot_to_inertial_t0(state0_nd, l_km, n_rad_s)
    p1_r0 = np.array([-mu * l_km, 0.0, 0.0])
    p1_v0 = np.array([0.0, -mu * l_km * n_rad_s, 0.0])
    p2_r0 = np.array([(1.0 - mu) * l_km, 0.0, 0.0])
    p2_v0 = np.array([0.0, (1.0 - mu) * l_km * n_rad_s, 0.0])

    sim = rebound.Simulation()
    sim.G = 1.0
    sim.integrator = "ias15"
    sim.integrator.epsilon = IAS15_EPSILON
    for m, r, v in ((gm1, p1_r0, p1_v0), (gm2, p2_r0, p2_v0), (0.0, sc_r0, sc_v0)):
        sim.add(
            m=m,
            x=float(r[0]),
            y=float(r[1]),
            z=float(r[2]),
            vx=float(v[0]),
            vy=float(v[1]),
            vz=float(v[2]),
        )

    t_period_s = period_nd * t_s
    n_total = n_periods * N_SAMP_PER_PERIOD
    p2_ref_nd = np.array([1.0 - mu, 0.0, 0.0])
    masses = np.array([gm1, gm2])

    times_nd: list[float] = []
    states_nd: list[NDArray[np.float64]] = []
    jacobis: list[float] = []
    rail_devs: list[float] = []
    com_devs: list[float] = []
    diverged = False
    try:
        for j in range(1, n_total + 1):
            sim.integrate(t_period_s * (j / N_SAMP_PER_PERIOD))
            ps = sim.particles
            p1 = np.array([ps[0].x, ps[0].y, ps[0].z])
            p2 = np.array([ps[1].x, ps[1].y, ps[1].z, ps[1].vx, ps[1].vy, ps[1].vz])
            sc = np.array([ps[2].x, ps[2].y, ps[2].z, ps[2].vx, ps[2].vy, ps[2].vz])
            if not (np.all(np.isfinite(sc)) and np.all(np.isfinite(p2))):
                diverged = True
                break
            theta = n_rad_s * float(sim.t)
            state_nd = _inertial_to_rot_nd(sc[:3], sc[3:], theta, l_km, n_rad_s)
            p2_nd = _inertial_to_rot_nd(p2[:3], p2[3:], theta, l_km, n_rad_s)
            com = (masses[0] * p1 + masses[1] * p2[:3]) / float(np.sum(masses))
            times_nd.append(float(sim.t) / t_s)
            states_nd.append(state_nd)
            jacobis.append(cr3bp.jacobi_constant(state_nd, mu))
            rail_devs.append(float(np.linalg.norm(p2_nd[:3] - p2_ref_nd)))
            com_devs.append(float(np.linalg.norm(com)) / l_km)
    except Exception:
        diverged = True

    floor_nd = float(max(max(rail_devs, default=0.0), max(com_devs, default=0.0), 1e-15))

    deltas: list[float] = []
    for k in range(1, n_periods + 1):
        idx = k * N_SAMP_PER_PERIOD - 1
        deltas.append(
            float(np.linalg.norm(states_nd[idx] - state0_nd))
            if idx < len(states_nd)
            else float("nan")
        )
    delta1 = deltas[0]

    # Per-period amplification from the member's Barden nu (theory, never loosened).
    nu = member.nu
    if abs(nu) <= 1.0:
        abs_lambda = 1.0
    else:
        disc = math.sqrt(nu * nu - 1.0)
        abs_lambda = max(abs(nu + disc), abs(nu - disc))

    departure_d = 3.0 * amp
    dl = [float(np.linalg.norm(s[:3] - state0_nd[:3])) for s in states_nd]
    t_bound_nd = n_periods * period_nd if not diverged else 0.0
    n_bound = len(states_nd)
    for i, d in enumerate(dl):
        if d > departure_d:
            t_bound_nd = times_nd[i]
            n_bound = i
            break

    jac_arr = np.asarray(jacobis[:n_bound]) if n_bound else np.asarray([j0])
    jacobi_drift_bounded = float(np.max(np.abs(jac_arr - j0)))

    delta1_eff = max(delta1, floor_nd) if math.isfinite(delta1) else float("inf")
    if abs_lambda > 1.0 and 0.0 < delta1_eff < departure_d:
        t_dep_pred_nd = period_nd + math.log(departure_d / delta1_eff) / math.log(abs_lambda)
    else:
        t_dep_pred_nd = float(n_periods * period_nd)

    noise_ok = floor_nd <= 0.01 * amp
    r1 = math.isfinite(delta1) and delta1 <= 0.1 * amp
    r2 = jacobi_drift_bounded <= 1e-9
    r3 = not diverged and t_bound_nd >= min(n_periods * period_nd, 0.7 * t_dep_pred_nd)
    if not noise_ok:
        verdict = "INCONCLUSIVE"
    elif r1 and r2 and r3:
        verdict = "PASS"
    else:
        verdict = "CHECK-FAILED"

    return InertialVerdict(
        verdict=verdict,
        delta_per_period_nd=tuple(deltas),
        delta1_nd=delta1,
        abs_lambda_per_period=abs_lambda,
        t_bound_nd=float(t_bound_nd),
        t_dep_pred_nd=float(t_dep_pred_nd),
        jacobi_drift_bounded=jacobi_drift_bounded,
        noise_floor_nd=floor_nd,
        diverged=diverged,
    )


# ---------------------------------------------------------------------------
# Per-seed work unit (parallelised)
# ---------------------------------------------------------------------------
@dataclass
class SeedResult:
    """Census for one seed: its two branches, classifications, and member tuples.

    NOVELTY ACCOUNTING (the #219 correction). Every seed in this campaign is a
    *published* family representative (Ross Table 3 / Arenstorf / a sourced #182
    Lyapunov). Continuing such a seed in the Jacobi constant therefore walks ALONG
    the already-published family -- each stepped member is the SAME family at a
    neighbouring energy, NOT a new discovery. The honest classification is:

      REPRODUCTION              : within dedup tol of the published seed point --
                                  the published representative itself, re-found.
      KNOWN-FAMILY-CONTINUATION : a distinct point on the SAME continuation branch
                                  as the published seed (continuous in C, no
                                  topology jump). A same-family cross-check member,
                                  NOT a discovery. (Was mislabelled NOVEL-SILVER.)
      NOVEL-SILVER              : a member on a branch that does NOT trace back to a
                                  sourced seed -- a genuinely DISTINCT new family.
                                  None occur in this campaign (all seeds are sourced).

    ``n_distinct_new_families`` is the only count that can support a discovery
    claim; in this campaign it is 0 (every branch is a sourced family continued).
    """

    label: str
    kind: str
    mu: float
    primary: str
    secondary: str
    seed_jacobi: float
    seed_period: float
    n_members: int = 0  # total kept across both directions (seed counted once)
    n_reproduction: int = 0
    n_known_family_continuation: int = 0  # same family as the sourced seed, stepped in C
    n_novel_silver: int = 0  # distinct NEW family (not a sourced-seed continuation)
    n_stable_continuation: int = 0  # |nu|<1 KNOWN-FAMILY-CONTINUATION members
    n_stable_novel: int = 0  # |nu|<1 genuinely-novel members (distinct new family)
    n_distinct_new_families: int = 0  # the only count that supports a discovery claim
    branches: list[dict[str, object]] = field(default_factory=list)
    # Members routed to the review queue (SILVER): both known-family continuations
    # (flagged as such) and any genuinely-novel members.
    review_members: list[dict[str, object]] = field(default_factory=list)
    empty: bool = False
    error: str = ""


def _is_reproduction(member: cc.BranchMember, seed: FamilySeed) -> bool:
    """True iff the member matches the seed's published member (sourced set dedup).

    The sourced set per seed is its own published (x0, C, T) -- the seed itself.
    A member within dedup tolerance of the seed is the published orbit (a sourced
    REPRODUCTION).
    """
    return (
        abs(member.x0 - seed.x0_seed) < DEDUP_X0_TOL
        and abs(member.jacobi - seed.jacobi) < DEDUP_C_TOL
        and abs(member.period - seed.period_guess) < DEDUP_T_TOL
    )


def run_seed(seed: FamilySeed) -> SeedResult:
    """Continue one seed both directions, run the inertial gate, classify members."""
    system = cr3bp.CR3BPSystem(
        mu=seed.mu, primary=seed.primary, secondary=seed.secondary, l_km=seed.l_km, t_s=seed.t_s
    )
    result = SeedResult(
        label=seed.label,
        kind=seed.kind,
        mu=seed.mu,
        primary=seed.primary,
        secondary=seed.secondary,
        seed_jacobi=seed.jacobi,
        seed_period=seed.period_guess,
    )
    try:
        corrected = cp.correct_symmetric_fixed_jacobi(
            system,
            seed.x0_seed,
            seed.jacobi,
            seed.period_guess,
            ydot0_sign=seed.ydot0_sign,
            half_crossings=seed.half_crossings,
            tol=1e-10,
        )
    except Exception as exc:
        result.error = f"seed corrector raised: {exc!r}"
        result.empty = True
        return result
    if not corrected.converged:
        result.error = "seed corrector did not converge"
        result.empty = True
        return result

    # Seed self-consistency: the corrected seed must reproduce the published period
    # to within SEED_PERIOD_TOL_FRAC. If it does not, the (x0, half_crossings, sign)
    # do NOT place the corrector on the intended family (e.g. a wrong crossing index
    # landed it on a different orbit) -- abandon rather than continue a spurious
    # branch and emit fabricated members.
    if abs(corrected.period - seed.period_guess) > SEED_PERIOD_TOL_FRAC * seed.period_guess:
        result.error = (
            f"seed off-family: corrected T={corrected.period:.6f} vs "
            f"published T={seed.period_guess:.6f} "
            f"(>{SEED_PERIOD_TOL_FRAC:.0%}); crossing index / sign do not match"
        )
        result.empty = True
        return result

    # Jacobi window centred on the seed's own C (so seeds below the global EM floor,
    # e.g. Arenstorf at C~2.86, still get a symmetric range around them).
    c_lo = max(JACOBI_MIN, seed.jacobi - JACOBI_HALFWIDTH)
    c_hi = min(JACOBI_MAX, seed.jacobi + JACOBI_HALFWIDTH)
    if seed.jacobi < JACOBI_MIN or seed.jacobi > JACOBI_MAX:
        # Seed sits outside the global EM prograde-cycler band: walk a local window.
        c_lo = seed.jacobi - JACOBI_HALFWIDTH
        c_hi = seed.jacobi + JACOBI_HALFWIDTH

    seen_keys: set[tuple[int, int]] = set()  # (round(x0/tol), round(C/tol)) dedup across dirs
    for direction in (1, -1):
        branch = cc.continue_family(
            system,
            corrected,
            direction=direction,
            d_jacobi=DEFAULT_DC,
            n_steps=DEFAULT_N_STEPS,
            min_jacobi=c_lo,
            max_jacobi=c_hi,
            half_crossings=seed.half_crossings,
            ydot0_sign=seed.ydot0_sign,
            seed_label=seed.label,
        )
        branch_members: list[dict[str, object]] = []
        for mi, member in enumerate(branch.members):
            # Skip the seed on the -1 branch (it is also the first +1 member).
            key = (round(member.x0 / 1e-9), round(member.jacobi / 1e-12))
            if mi == 0 and direction == -1 and key in seen_keys:
                continue
            seen_keys.add(key)
            result.n_members += 1

            reproduction = _is_reproduction(member, seed)
            inertial = inertial_crosscheck(member, seed)
            # NOVELTY ACCOUNTING (the #219 correction). The branch is built by
            # continuing a *sourced* seed (Ross / Arenstorf / sourced #182 Lyapunov)
            # that already passed the off-family guard, so it traces to a published
            # family by construction. A non-reproduction member on this branch is
            # therefore the SAME family at a neighbouring energy -- a KNOWN-FAMILY
            # CONTINUATION (a same-model cross-check), NOT a discovery. The only way
            # to earn NOVEL-SILVER is a member on a branch that does NOT trace back
            # to a sourced seed; none occur here (all seeds are sourced), so this
            # campaign yields 0 genuinely-distinct new families.
            if reproduction:
                cls = "REPRODUCTION"
            elif inertial.verdict != "PASS":
                # Continuation member that fails the independence gate is NOT routed.
                cls = f"CONTINUATION-REJECTED({inertial.verdict})"
            else:
                cls = "KNOWN-FAMILY-CONTINUATION"

            rec = {
                "index": mi,
                "direction": direction,
                "classification": cls,
                "family_label": seed.label,  # the sourced family this point belongs to
                "x0": member.x0,
                "ydot0": member.ydot0,
                "state0_nd": member.state0.tolist(),
                "period_nd": member.period,
                "jacobi": member.jacobi,
                "nu": member.nu,
                "abs_lambda": member.abs_lambda,
                "stable": member.stable,
                "crossing_residual": member.crossing_residual,
                "radau_djacobi": member.radau_djacobi,
                "max_speed_nd": member.max_speed_nd,
                "amplitude_nd": member.amplitude_nd,
                "inertial_verdict": inertial.verdict,
                "inertial_delta1_nd": inertial.delta1_nd,
                "inertial_jacobi_drift_bounded": inertial.jacobi_drift_bounded,
                "inertial_t_bound_nd": inertial.t_bound_nd,
                "inertial_t_dep_pred_nd": inertial.t_dep_pred_nd,
                "inertial_abs_lambda_theory": inertial.abs_lambda_per_period,
                "inertial_noise_floor_nd": inertial.noise_floor_nd,
            }
            branch_members.append(rec)
            if reproduction:
                result.n_reproduction += 1
            elif cls == "KNOWN-FAMILY-CONTINUATION":
                result.n_known_family_continuation += 1
                # Routed to the review queue as SILVER, but flagged same-family
                # (NOT a discovery -- never inflates the new-family count).
                result.review_members.append(rec)
                if member.stable:
                    result.n_stable_continuation += 1

        result.branches.append(
            {
                "direction": direction,
                "dC": branch.d_jacobi,
                "n_steps_taken": branch.n_steps_taken,
                "n_rejected": branch.n_rejected,
                "stop_reason": branch.stop_reason.value,
                "n_members": len(branch_members),
                "members": branch_members,
            }
        )

    # EMPTY-for-novelty (the registry sense): the seed surfaced NO genuinely-distinct
    # new family. Every seed here continues a sourced family, so n_distinct_new_families
    # is always 0 -- the whole campaign is EMPTY for novelty by construction.
    result.n_distinct_new_families = result.n_novel_silver
    result.empty = result.n_distinct_new_families == 0
    return result


# ---------------------------------------------------------------------------
# Disposition / reporting
# ---------------------------------------------------------------------------
def _now_iso() -> str:
    return datetime.now(UTC).isoformat()


def write_review_queue(results: list[SeedResult]) -> int:
    """Append every routed continuation member to the review queue. SILVER only.

    Members written here are KNOWN-FAMILY-CONTINUATION points (same family as the
    sourced seed, stepped in C) -- recorded as SILVER same-model cross-check members,
    explicitly flagged NOT a discovery. Any genuinely-novel member (distinct new
    family) would also be written here; none occur in this campaign.
    """
    n = 0
    REVIEW_QUEUE_PATH.parent.mkdir(parents=True, exist_ok=True)
    with REVIEW_QUEUE_PATH.open("a", encoding="utf-8") as fh:
        for res in results:
            for m in res.review_members:
                record = {
                    "candidate_id": (
                        f"cr3bp-cont-{res.secondary.lower()}-{res.label}"
                        f"-C{m['jacobi']:.8f}-T{m['period_nd']:.6f}"
                    ),
                    "verdict_tier": "SILVER",
                    "method": METHOD_TAG,
                    "git_sha": _GIT_SHA,
                    "seed_label": res.label,
                    "seed_kind": res.kind,
                    "primary": res.primary,
                    "secondary": res.secondary,
                    "mass_ratio": res.mu,
                    "model_assumption": "cr3bp",
                    "frame": "rotating barycentric planar (PCR3BP)",
                    **m,
                    "novelty_scope": (
                        f"SAME family as the sourced seed '{res.label}' (a continuation "
                        "point at a neighbouring Jacobi constant), NOT a discovery. The "
                        "sourced set (Ross 5 / Arenstorf / existing SILVER) is the family "
                        "representative; this is another member of that published family. "
                        "JPL 3-body catalog NOT acquired (#116) -- no literature-novelty claim"
                    ),
                    "gates_passed": [
                        "corrector-perpendicular-crossing(<1e-10)",
                        "period-bounds",
                        "equilibrium-rejection",
                        "dedup-distinct",
                        "jacobi-conservation(<=1e-10)",
                        "crosscheck-radau",
                        f"inertial-crosscheck({INERTIAL_TAG})",
                    ],
                    "t_added": _now_iso(),
                    "notes": (
                        "Jacobi-continuation member from a sourced seed; full gauntlet + "
                        "independent inertial REBOUND/IAS15 cross-check. SILVER pending "
                        "human review. NO catalogue writeback."
                    ),
                }
                fh.write(json.dumps(record, ensure_ascii=True) + "\n")
                n += 1
    return n


def write_empty_regions(results: list[SeedResult]) -> int:
    """Append a method-versioned EMPTY record for every seed that yielded no new family.

    EMPTY-for-novelty means no genuinely-DISTINCT new family surfaced. Reproductions
    AND known-family continuations are both EXPECTED, sourced cross-checks (a success),
    not new families -- so a seed that only reproduces/continues a sourced family is
    EMPTY for the novelty registry.
    """
    n = 0
    EMPTY_REGIONS_PATH.parent.mkdir(parents=True, exist_ok=True)
    today = datetime.now(UTC).strftime("%Y-%m-%d")
    with EMPTY_REGIONS_PATH.open("a", encoding="utf-8") as fh:
        for res in results:
            # A seed is EMPTY for novelty unless it produced a genuinely-distinct new
            # family. Reproductions and same-family continuations are sourced
            # cross-checks (a success), not new families.
            if res.n_distinct_new_families > 0:
                continue
            record = {
                "region_id": f"cr3bp-continuation-{res.label}-{today}",
                "family": (
                    f"CR3BP symmetric periodic-orbit continuation ({res.primary}/{res.secondary})"
                ),
                "centre": res.secondary,
                "method_capability": {
                    "genome": (
                        "natural-parameter (Jacobi) continuation of symmetric "
                        "perpendicular-crossing orbits"
                    ),
                    "corrector": "correct_symmetric_fixed_jacobi + barden_stability (Ross 2025)",
                    "method": METHOD_TAG,
                    "capability_tags": ["cr3bp", "planar", "symmetric", "jacobi-continuation"],
                    "git_sha": _GIT_SHA,
                },
                "search_extent": {
                    "seed_label": res.label,
                    "seed_jacobi": res.seed_jacobi,
                    "seed_period_nd": res.seed_period,
                    "dC": DEFAULT_DC,
                    "n_steps_each_direction": DEFAULT_N_STEPS,
                    "jacobi_range": [JACOBI_MIN, JACOBI_MAX],
                    "directions": [1, -1],
                    "branches": [
                        {
                            k: b[k]
                            for k in (
                                "direction",
                                "n_steps_taken",
                                "n_rejected",
                                "stop_reason",
                                "n_members",
                            )
                        }
                        for b in res.branches
                    ],
                },
                "result": {
                    "members_kept": res.n_members,
                    "reproductions": res.n_reproduction,
                    "known_family_continuation": res.n_known_family_continuation,
                    "stable_continuation": res.n_stable_continuation,
                    "distinct_new_families": res.n_distinct_new_families,
                    "novel_silver": res.n_novel_silver,
                    "stable_novel": res.n_stable_novel,
                    "error": res.error,
                },
                "verdict": (
                    "EMPTY -- no genuinely-DISTINCT new family under this method+version in "
                    "the scanned Jacobi range (only reproductions / same-family continuations)"
                ),
                "interpretation": (
                    "Continuation reproduced the seed family and/or stepped along it in the "
                    "Jacobi constant (KNOWN-FAMILY-CONTINUATION -- a same-model sourced "
                    "cross-check), without surfacing any DISTINCT new family. Continuation "
                    "points of a published family are the SAME family, not discoveries. JPL "
                    "3-body catalog not acquired (#116); novelty scope is 'not in our "
                    "sourced set AND not a continuation point of a sourced family'."
                ),
                "run": {"date": today, "git_sha": _GIT_SHA},
            }
            fh.write(json.dumps(record, ensure_ascii=True) + "\n")
            n += 1
    return n


def build_note(results: list[SeedResult], elapsed_s: float, n_review: int, n_empty: int) -> str:
    ts = datetime.now(tz=UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    total_members = sum(r.n_members for r in results)
    total_repro = sum(r.n_reproduction for r in results)
    total_cont = sum(r.n_known_family_continuation for r in results)
    total_stable_cont = sum(r.n_stable_continuation for r in results)
    total_distinct = sum(r.n_distinct_new_families for r in results)
    lines: list[str] = [
        "# CR3BP stable-family continuation search -- results (engine-validation campaign)",
        "",
        f"**Run timestamp:** {ts}   elapsed {elapsed_s:.1f}s",
        "**Script:** `scripts/cr3bp_family_search.py`",
        f"**Method:** `{METHOD_TAG}`  (git `{_GIT_SHA or '?'}`)",
        f"**Inertial gate:** `{INERTIAL_TAG}` (REBOUND/IAS15, independent frame + code path)",
        "**Spec:** `docs/superpowers/specs/2026-06-12-cr3bp-family-continuation-search-design.md`",
        "**Status:** SILVER-only. NO catalogue writeback.",
        "",
        "## Headline (honest framing -- read this first)",
        "",
        f"- **Distinct NEW families found: {total_distinct}.** This is the only count that",
        "  can support a discovery claim, and it is ZERO. Every seed in this campaign is a",
        "  *published* family representative (Ross Table 3 / Arenstorf / a sourced #182",
        "  Lyapunov). Continuing such a seed in the Jacobi constant walks ALONG the already-",
        "  published family -- each stepped member is the SAME family at a neighbouring",
        "  energy, NOT a new orbit. A clean ZERO with the engine validated is a SUCCESS, not",
        "  a failure: it is exactly what the orbit-closure discipline expects when every",
        "  seed is sourced and the JPL 3-body catalogue (#116) is not yet acquired.",
        "- **This corrects an earlier inflated count.** A prior pass labelled every non-",
        "  point-reproduction member 'NOVEL-SILVER' (peaking at 83 'novel' / 30 'stable",
        "  novel'), because the dedup compared distance to a single published *member* while",
        "  Ross prints one representative per family -- so every other point on the same",
        "  family curve counted as a discovery. The classification now keys on FAMILY",
        "  identity / branch continuity, not point-distance:",
        "",
        f"- Seeds run: **{len(results)}**",
        f"- Members kept (gauntlet-passing, across both directions): **{total_members}**",
        f"- REPRODUCTION (published representative re-found): **{total_repro}**",
        f"- KNOWN-FAMILY-CONTINUATION (same sourced family, stepped in C; NOT a discovery):"
        f" **{total_cont}**",
        f"  -- of which linearly STABLE (|nu|<1): {total_stable_cont} (additional members of",
        "  Ross's *named* stable windows -- a same-model cross-check, not new stable families).",
        f"- NOVEL-SILVER (member on a branch NOT tracing to a sourced seed = distinct new"
        f" family): **{sum(r.n_novel_silver for r in results)}**",
        f"- Review-queue rows written: {n_review}  (`{REVIEW_QUEUE_PATH.relative_to(REPO_ROOT)}`)",
        "  -- all KNOWN-FAMILY-CONTINUATION, flagged same-family (NOT discoveries).",
        f"- EMPTY-region records written: {n_empty}  (`data/empty_regions.jsonl`)",
        "",
        "## Model & scope (caveats)",
        "",
        "- **Planar CR3BP (PCR3BP) only.** The richest stable families (halos, DROs) are",
        "  3D and need a 3D corrector (deferred).",
        "- **Novelty scope = 'not in our sourced set AND not a continuation point of a",
        "  sourced family'.** A genuinely-novel member would have to sit on a branch that",
        "  does not trace back to a published seed. The JPL Three-Body Periodic Orbit",
        "  Catalog is NOT acquired (#116); even a future genuinely-novel member is 'new vs",
        "  our sourced set', NOT a literature-novelty claim.",
        "- **Reproduction / same-family continuation is the engine-validation SUCCESS.**",
        "  Re-deriving Ross's families (and stepping smoothly along them with nu varying",
        "  continuously and Jacobi conserved) validates the continuation engine against",
        "  sourced truth. That -- not a discovery count -- is what this campaign delivers.",
        "",
        "## Per-seed census",
        "",
        "| seed | system | mu | members | repro | known-fam-cont | stable-cont "
        "| new families | branches (dir: stop @ steps) |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        br = "; ".join(
            f"{b['direction']:+d}: {b['stop_reason']}@{b['n_steps_taken']}" for b in r.branches
        ) or (r.error or "n/a")
        lines.append(
            f"| {r.label} | {r.primary}/{r.secondary} | {r.mu:.6e} | {r.n_members} "
            f"| {r.n_reproduction} | {r.n_known_family_continuation} | {r.n_stable_continuation} "
            f"| {r.n_distinct_new_families} | {br} |"
        )
    lines += [
        "",
        "## Inertial cross-check (false-consensus independence gate)",
        "",
        "Each kept member is re-propagated in the inertial REBOUND/IAS15 harness (a",
        "different code path AND frame than the rotating-frame DOP853 corrector) and",
        "graded against LINEAR-INSTABILITY THEORY -- the per-period amplification is the",
        "monodromy eigenvalue magnitude `|lambda| = |nu + sqrt(nu^2-1)|` from the member's",
        "Barden nu. R1 (one-period recurrence <= 0.1A), R2 (Jacobi drift <= 1e-9 over",
        "the bounded span), R3 (observed departure no EARLIER than the linear prediction",
        "`t_dep = T + ln(3A/delta1)/ln|lambda|`). The band is NEVER loosened. A member",
        "that fails this gate is NOT routed (classified CONTINUATION-REJECTED).",
        "",
        "**Honest limit of this gate (do not overread a PASS).** The corrector lands",
        "every member to a ~1e-10 perpendicular-crossing residual, so the inertial",
        "one-period recurrence delta1 is ~1e-9..1e-8 -- and even a strongly unstable",
        "member (|nu|~10^2, |lambda|~10^2..10^3) does not amplify that tiny seed past 3A",
        "within the 5-period span, so it stays numerically bounded and R3 records",
        "departure LATER than predicted (a PASS). The inertial gate therefore confirms",
        "each member is a GENUINE periodic orbit through an independent integrator+frame",
        "(it re-closes and conserves Jacobi) -- it is the false-consensus *consistency*",
        "gate. It does NOT independently certify STABILITY: stability is the Barden nu",
        "(reported per member), which IS discriminating here (nu spans -0.0 stable to",
        "+360 wildly unstable across these branches). Read 'inertial PASS' as 'a real",
        "orbit, cross-checked', and the STABLE verdict as 'from nu', never conflated.",
        "",
    ]
    if total_stable_cont > 0:
        lines += [
            "## STABLE same-family continuation members (NOT discoveries, NOT catalogued)",
            "",
            "**These are NOT new families.** Every stable member below is a continuation",
            "member of one of Ross's *named* (k1,k2) families, at a Jacobi constant a few",
            "steps off the published nu=0 midpoint but INSIDE that family's finite stable",
            "window (Ross Table 3 gives the window widths; e.g. (3,3) spans ~2041 km in",
            "perilune). They are additional members of the published stable WINDOW -- a same-",
            "model cross-check, a SUCCESS -- not new stable families. No literature-novelty",
            "is claimed (JPL 3-body catalog not acquired, #116). NEVER catalogue.",
            "",
            "| family | C | T (nd) | nu | x0 | ydot0 | inertial |",
            "|---|---|---|---|---|---|---|",
        ]
        for r in results:
            for m in r.review_members:
                if m["stable"]:
                    lines.append(
                        f"| {r.label} | {m['jacobi']:.8f} | {m['period_nd']:.6f} | {m['nu']:+.5f} "
                        f"| {m['x0']:.8f} | {m['ydot0']:.8f} | {m['inertial_verdict']} |"
                    )
        lines.append("")
    if total_distinct == 0:
        lines += [
            "## Distinct new families",
            "",
            "**None.** By construction: every seed is a sourced family representative, so",
            "continuation can only re-derive or step along a published family. No branch in",
            "this campaign traces to anything other than a sourced seed. Recorded as method-",
            "versioned EMPTY-for-novelty (re-sweepable when the JPL catalog (#116) or a 3D",
            "corrector ships, or when an UN-sourced seed is introduced).",
            "",
        ]
    lines += [
        "## Discipline statement",
        "",
        "Sourced goldens only (Ross Table 3 / Arenstorf); same-model validation; novelty",
        "keyed on FAMILY identity / branch continuity (a continuation point of a sourced",
        "family is the SAME family, not a discovery), not point-distance; SILVER-only with",
        "NO catalogue writeback; every Jacobi range/step/stop-reason logged (no silent caps).",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description="CR3BP Jacobi-continuation family search")
    parser.add_argument(
        "--no-writeback", action="store_true", help="skip review-queue + empty-region writes"
    )
    parser.add_argument("--workers", type=int, default=None, help="process-pool size")
    parser.add_argument("--report", type=Path, default=None, help="write the runlog here")
    args = parser.parse_args()

    seeds = all_seeds()
    t0 = datetime.now(tz=UTC)
    print(f"cr3bp_family_search  {t0.isoformat()}  method={METHOD_TAG}  git={_GIT_SHA or '?'}")
    print(f"  seeds: {len(seeds)}  ({', '.join(s.label for s in seeds)})")
    print(
        f"  Jacobi range [{JACOBI_MIN}, {JACOBI_MAX}]  dC={DEFAULT_DC}  "
        f"n_steps={DEFAULT_N_STEPS}/dir"
    )

    results: list[SeedResult]
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        results = list(ex.map(run_seed, seeds))

    elapsed = (datetime.now(tz=UTC) - t0).total_seconds()
    for r in results:
        print(
            f"  {r.label:<26} members={r.n_members:<3} repro={r.n_reproduction:<3} "
            f"novel={r.n_novel_silver:<3} stable-novel={r.n_stable_novel:<3} "
            f"{'EMPTY' if r.empty else ''}{(' ERR:' + r.error) if r.error else ''}"
        )

    n_review = 0
    n_empty = 0
    if args.no_writeback:
        print("\nWriteback SKIPPED (--no-writeback).")
    else:
        n_review = write_review_queue(results)
        n_empty = write_empty_regions(results)
        print(f"\nReview-queue rows written: {n_review} -> {REVIEW_QUEUE_PATH}")
        print(f"EMPTY-region records written: {n_empty} -> {EMPTY_REGIONS_PATH}")

    NOTE_PATH.parent.mkdir(parents=True, exist_ok=True)
    NOTE_PATH.write_text(build_note(results, elapsed, n_review, n_empty), encoding="utf-8")
    print(f"Results note written: {NOTE_PATH}")

    RUNLOG_PATH.mkdir(parents=True, exist_ok=True)
    runlog_file = RUNLOG_PATH / f"run-{t0.strftime('%Y%m%dT%H%M%SZ')}.json"
    runlog_file.write_text(
        json.dumps(
            {
                "method": METHOD_TAG,
                "git_sha": _GIT_SHA,
                "timestamp": t0.isoformat(),
                "elapsed_s": elapsed,
                "jacobi_range": [JACOBI_MIN, JACOBI_MAX],
                "dC": DEFAULT_DC,
                "n_steps": DEFAULT_N_STEPS,
                "results": [asdict(r) for r in results],
            },
            ensure_ascii=True,
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"Runlog written: {runlog_file}")

    total_stable_novel = sum(r.n_stable_novel for r in results)
    print(f"\nSTABLE novel members (|nu|<1): {total_stable_novel}")
    print("NO catalogue writeback performed.")
    if args.report is not None:
        args.report.write_text(
            build_note(results, elapsed, n_review, n_empty) + "\n", encoding="utf-8"
        )


if __name__ == "__main__":
    from cyclerfinder.search.outcome_log import enable_default_outcome_log

    enable_default_outcome_log("cr3bp_family_search")
    main()
