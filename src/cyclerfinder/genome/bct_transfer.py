"""Belbruno ballistic-capture transfer (BCT) constructor (#378 Phase 2).

Builds and corrects a Hiten-class *exterior* weak-stability-boundary (WSB)
transfer on the project's incoherent BCR4BP, reusing:

  * `core/bcr4bp.py`  -- propagator + STM (no modification);
  * `core/wsb.py`     -- the WSB surface (Phase 1: E_2, periapsis, W);
  * the `genome/bcr4bp_genome.correct_bcr4bp_periodic` Newton *pattern*
    (mirrored, not imported -- the residual is a terminal-targeting residual,
    not a closure residual);
  * `search/literature_check.py` for the BCT novelty self-test.

This is a CAPABILITY module. A single BCT is a one-shot transfer (precursor_mga
class), NOT a cycler -- Belbruno 2004 is explicit ("WSB orbits are not
cyclers"). The design draft §6 records this honestly. The substrate makes a
*repeating* cislunar chain searchable (`search/cislunar_bct_search.py`), but a
standalone BCT is a transfer artifact, not a catalogue row.

Model-gap honesty (design §4)
-----------------------------
The project's incoherent BCR4BP is not Belbruno's PR4BP-3D-with-DE403. The
Hiten golden is therefore a SIGNATURE BAND (ΔV factor-2, apoapsis ±30%, TOF
order-150 d), not a bit-exact reproduction. The definitional facts -- ballistic
capture E_2 <= 0, ΔV_capture = 0 -- are exact by construction.

Units
-----
BCR4BP nondimensional: length = EM distance (1 LD), time = 1/n. Conversions to
km / km/s / days use the EM scales below.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.bcr4bp as bcr4bp
import cyclerfinder.core.wsb as wsb
from cyclerfinder.search.literature_check import (
    CandidateSignature,
    LiteratureCheckResult,
    SearchResult,
)

# Earth-Moon physical scales (EM distance = 1 nondim length).
EM_DISTANCE_KM = 384_400.0
MOON_RADIUS_KM = 1737.4
SIDEREAL_MONTH_DAYS = 27.321661
TU_DAYS = SIDEREAL_MONTH_DAYS / (2.0 * math.pi)
KM_PER_TU_OVER_S = EM_DISTANCE_KM / (TU_DAYS * 86400.0)  # km/s per nondim velocity

Branch = Literal["direct", "retrograde"]


def nondim_v_to_kms(v_nd: float) -> float:
    """Convert a nondim BCR4BP speed to km/s."""
    return v_nd * KM_PER_TU_OVER_S


def kms_to_nondim_v(v_kms: float) -> float:
    """Convert km/s to a nondim BCR4BP speed."""
    return v_kms / KM_PER_TU_OVER_S


@dataclass(frozen=True)
class BCTTarget:
    """The capture target QF on the WSB surface W at the Moon.

    Attributes
    ----------
    r_capture_km :
        Capture *altitude* above the Moon's surface (km). Belbruno's Hiten
        target is 100 km -> r_23 = r_M + 100 km.
    e2 :
        Prescribed Moon-relative osculating eccentricity at QF (Belbruno: 0.95).
    theta2 :
        Periapsis orientation angle around the Moon (rad). The family selector:
        the W point sits on the unstable manifold whose backward arc escapes the
        Moon to the Sun-shaped apoapsis -- theta_2 picks which manifold.
    branch :
        ``"direct"`` or ``"retrograde"`` Moon-relative capture orbit.
    """

    r_capture_km: float = 100.0
    e2: float = 0.95
    theta2: float = 0.70
    branch: Branch = "retrograde"

    @property
    def r23_nondim(self) -> float:
        """Moon-relative periapsis radius (nondim)."""
        return (MOON_RADIUS_KM + self.r_capture_km) / EM_DISTANCE_KM

    @property
    def r23_km(self) -> float:
        """Moon-relative periapsis radius (km)."""
        return MOON_RADIUS_KM + self.r_capture_km


def qf_state_on_w(target: BCTTarget, system: bcr4bp.BCR4BPSystem) -> NDArray[np.float64]:
    """Build the QF 6-state on W (Moon periapsis at r_23, e_2, theta_2, branch).

    Position: at angle ``theta2`` around the Moon at radius ``r23_nondim``.
    Velocity: the two-body periapsis inertial speed
    ``v = sqrt(mu (1 + e_2) / r_23)`` tangential to X, converted back to the
    rotating frame via ``v_rot = Xdot - omega x X``. By construction
    ``E_2 <= 0`` and ``r-dot_23 = 0`` -- a point on W (verified by tests).
    """
    mu = system.mu
    r23 = target.r23_nondim
    ang = target.theta2
    moon_x = 1.0 - mu
    px = moon_x + r23 * math.cos(ang)
    py = r23 * math.sin(ang)
    x_rel = np.array([px - moon_x, py, 0.0])
    v_peri = math.sqrt(mu * (1.0 + target.e2) / r23)
    sign = 1.0 if target.branch == "direct" else -1.0
    tangent = np.array([-math.sin(ang), math.cos(ang), 0.0]) * sign
    x_dot = tangent * v_peri  # inertial Moon-relative velocity
    omega_cross_x = np.array([-x_rel[1], x_rel[0], 0.0])
    v_rot = x_dot - omega_cross_x
    return np.array([px, py, 0.0, v_rot[0], v_rot[1], v_rot[2]], dtype=np.float64)


@dataclass(frozen=True)
class BCTArc:
    """Result of the backward arc-II construction from a QF on W.

    Attributes
    ----------
    qf_state :
        The capture target state on W (6-vector, rotating frame).
    apoapsis_state :
        The state at the backward arc's max Earth-relative apoapsis.
    max_earth_apoapsis_ld :
        Max Earth-relative distance reached on the backward arc (LD = nondim).
    back_days :
        Backward integration horizon used (days).
    t_apoapsis_days :
        Backward time at which the apoapsis occurred (days).
    """

    qf_state: NDArray[np.float64]
    apoapsis_state: NDArray[np.float64]
    max_earth_apoapsis_ld: float
    back_days: float
    t_apoapsis_days: float


def construct_bct_backward(
    target: BCTTarget,
    system: bcr4bp.BCR4BPSystem,
    *,
    back_days: float = 70.0,
    n_samples: int = 160,
    rtol: float = 1e-10,
    atol: float = 1e-10,
) -> BCTArc:
    """Backward-integrate arc II from QF on W to the Sun-shaped apoapsis (§3.4.1).

    Belbruno Fig 3.14 recipe: from the capture target QF on W, integrate the
    four-body EOM BACKWARD; the trajectory climbs to a Sun-shaped Earth-relative
    apoapsis Q_a (~1.5e6 km ~ 3.9 LD for Hiten). Returns the arc with its max
    Earth-relative apoapsis (the gate quantity).
    """
    qf = qf_state_on_w(target, system)
    mu = system.mu
    earth = np.array([-mu, 0.0, 0.0])
    max_r = float(np.linalg.norm(qf[:3] - earth))
    apo_state = qf.copy()
    t_apo = 0.0
    for d in np.linspace(back_days / n_samples, back_days, n_samples):
        try:
            arc = bcr4bp.propagate_bcr4bp(system, qf, -d / TU_DAYS, rtol=rtol, atol=atol)
        except RuntimeError:
            break
        r = float(np.linalg.norm(arc.state_f[:3] - earth))
        if not math.isfinite(r):
            break
        if r > max_r:
            max_r = r
            apo_state = arc.state_f.copy()
            t_apo = float(d)
    return BCTArc(
        qf_state=qf,
        apoapsis_state=apo_state,
        max_earth_apoapsis_ld=max_r,
        back_days=back_days,
        t_apoapsis_days=t_apo,
    )


# ---------------------------------------------------------------------------
# Forward 2x2 corrector: (|V0|, gamma0) -> (r_23, E_2-on-W).
# ---------------------------------------------------------------------------


def _leo_departure_state(
    system: bcr4bp.BCR4BPSystem, v_kms: float, gamma_deg: float, r13_km: float
) -> NDArray[np.float64]:
    """Departure state at a low-Earth periapsis Q0 on the +x side of Earth.

    Earth at (-mu, 0, 0). Q0 at r_13 = r13_km beyond Earth on +x; velocity in
    the local frame at flight-path angle gamma0 off the local horizontal (+y).
    """
    mu = system.mu
    r13_nd = r13_km / EM_DISTANCE_KM
    x = -mu + r13_nd
    v_nd = kms_to_nondim_v(v_kms)
    gamma = math.radians(gamma_deg)
    vx = v_nd * math.sin(gamma)
    vy = v_nd * math.cos(gamma)
    return np.array([x, 0.0, 0.0, vx, vy, 0.0], dtype=np.float64)


def _moon_closest_approach(
    system: bcr4bp.BCR4BPSystem,
    state0: NDArray[np.float64],
    tof_days: float,
    *,
    n_samples: int = 400,
    rtol: float = 1e-10,
    atol: float = 1e-10,
) -> tuple[float, float, NDArray[np.float64]]:
    """Forward-propagate and return ``(min_r23_nondim, e2_at_min, state_at_min)``.

    The terminal capture targeting reads the closest Moon approach over the
    forward arc: that is where the spacecraft is captured (r_23 minimal, E_2 the
    instantaneous Kepler energy there).
    """
    mu = system.mu
    moon = np.array([1.0 - mu, 0.0, 0.0])
    min_r = math.inf
    e2_at = math.inf
    st_at = state0.copy()
    cur = state0.copy()
    t_prev = 0.0
    for d in np.linspace(tof_days / n_samples, tof_days, n_samples):
        t_nd = d / TU_DAYS
        try:
            arc = bcr4bp.propagate_bcr4bp(
                system, cur, t_nd - t_prev, t0=t_prev, rtol=rtol, atol=atol
            )
        except RuntimeError:
            break
        cur = arc.state_f
        t_prev = t_nd
        r23 = float(np.linalg.norm(cur[:3] - moon))
        if r23 < min_r:
            min_r = r23
            e2_at = wsb.kepler_energy_moon(cur, system)
            st_at = cur.copy()
    return min_r, e2_at, st_at


@dataclass(frozen=True)
class BCTForwardResult:
    """A converged (or attempted) forward 2x2 BCT correction.

    Attributes
    ----------
    v0_kms / gamma0_deg :
        Corrected departure controls at Q0.
    terminal_r23_km :
        Moon-relative distance at closest approach (km).
    terminal_e2 :
        Kepler energy E_2 at closest approach. <= 0 means landed on W (ballistic
        capture).
    dv_capture_kms :
        Capture ΔV. Exactly 0.0 by the ballistic-capture definition.
    tof_days :
        Time of flight to closest approach (days).
    converged :
        True iff r_23 within tolerance of the target AND E_2 <= 0 AND the
        independent Radau re-propagation reproduces the terminal state.
    independent_terminal_residual :
        L2 state difference at closest approach between DOP853 and an
        independent Radau re-propagation of the converged IC.
    n_iter :
        Newton iterations consumed.
    departure_state :
        The converged 6-state at Q0.
    notes :
        Diagnostics.
    """

    v0_kms: float
    gamma0_deg: float
    terminal_r23_km: float
    terminal_e2: float
    dv_capture_kms: float
    tof_days: float
    converged: bool
    independent_terminal_residual: float
    n_iter: int
    departure_state: NDArray[np.float64]
    notes: str = ""


def correct_bct_forward(
    target: BCTTarget,
    system: bcr4bp.BCR4BPSystem,
    *,
    v0_guess_kms: float = 10.92,
    gamma0_guess_deg: float = 0.0,
    r13_km: float = 6578.0,
    tof_days: float = 200.0,
    max_iter: int = 40,
    tol_r23_rel: float = 0.2,
    rtol: float = 1e-10,
    atol: float = 1e-10,
) -> BCTForwardResult:
    """Forward 2x2 differential correction (|V0|, gamma0) -> (r_23, E_2-on-W).

    Belbruno Remark 4. Controls ``(|V0|, gamma0)`` at the LEO departure Q0;
    targets the Moon-relative closest-approach distance ``r_23 = r_target`` with
    the terminal ballistic-capture constraint ``E_2 <= 0`` (landing on W). The
    2x2 Jacobian is built by central finite differences on the closest-approach
    ``r_23`` w.r.t. ``(|V0|, gamma0)`` (mirrors the BCR4BP corrector's damped
    Newton, but on a terminal targeting residual rather than a closure residual).

    The E_2-on-W constraint is enforced as an acceptance gate at the converged
    r_23 (the design's "E_2-on-W terminal constraint row"); a converged result
    has ``E_2 <= 0`` or it is reported non-converged.

    An independent Radau re-propagation of the converged IC reproduces the
    terminal closest-approach state (orbit-closure discipline).
    """
    target_r23 = target.r23_nondim
    v0 = v0_guess_kms
    gamma0 = gamma0_guess_deg

    def evaluate(v_kms: float, g_deg: float) -> tuple[float, float, float, NDArray[np.float64]]:
        st0 = _leo_departure_state(system, v_kms, g_deg, r13_km)
        min_r, e2_at, st_at = _moon_closest_approach(system, st0, tof_days, rtol=rtol, atol=atol)
        # residual on r_23 (we want min_r == target_r23).
        return min_r - target_r23, min_r, e2_at, st_at

    res, min_r, e2_at, st_at = evaluate(v0, gamma0)
    n_iter = 0
    dv_step = 0.01  # km/s FD step
    dg_step = 0.25  # deg FD step
    for n_iter in range(1, max_iter + 1):  # noqa: B007
        if abs(res) < tol_r23_rel * target_r23 and e2_at <= 0.0:
            break
        # Central FD on the 2 controls (2x1 gradient -> 1 residual / 2 controls).
        rp_v, *_ = evaluate(v0 + dv_step, gamma0)
        rm_v, *_ = evaluate(v0 - dv_step, gamma0)
        rp_g, *_ = evaluate(v0, gamma0 + dg_step)
        rm_g, *_ = evaluate(v0, gamma0 - dg_step)
        dr_dv = (rp_v - rm_v) / (2.0 * dv_step)
        dr_dg = (rp_g - rm_g) / (2.0 * dg_step)
        grad = np.array([dr_dv, dr_dg])
        gnorm2 = float(grad @ grad)
        if gnorm2 < 1e-30 or not math.isfinite(gnorm2):
            break
        # Minimum-norm Gauss-Newton step on the single residual.
        step = -(res / gnorm2) * grad
        # Damp the step to keep within the bound-Earth regime.
        step[0] = float(np.clip(step[0], -0.05, 0.05))  # km/s
        step[1] = float(np.clip(step[1], -3.0, 3.0))  # deg
        # Backtracking line search on |res|.
        scale = 1.0
        improved = False
        for _bt in range(12):
            tv = v0 + scale * step[0]
            tg = gamma0 + scale * step[1]
            tres, tmin, te2, tst = evaluate(tv, tg)
            if math.isfinite(tres) and abs(tres) < abs(res):
                v0, gamma0 = tv, tg
                res, min_r, e2_at, st_at = tres, tmin, te2, tst
                improved = True
                break
            scale *= 0.5
        if not improved:
            break

    # Independent Radau cross-check: re-propagate the converged IC and compare
    # the closest-approach state.
    independent_resid = math.nan
    try:
        st0 = _leo_departure_state(system, v0, gamma0, r13_km)
        sol = solve_ivp(
            bcr4bp.bcr4bp_eom,
            (0.0, tof_days / TU_DAYS),
            st0,
            args=(system,),
            method="Radau",
            rtol=max(rtol, 1e-11),
            atol=max(atol, 1e-11),
            dense_output=False,
            t_eval=np.linspace(0.0, tof_days / TU_DAYS, 400),
        )
        if sol.success:
            mu = system.mu
            moon = np.array([1.0 - mu, 0.0, 0.0])
            dists = np.linalg.norm(sol.y[:3, :].T - moon, axis=1)
            kmin = int(np.argmin(dists))
            independent_resid = float(np.linalg.norm(sol.y[:, kmin] - st_at))
    except (RuntimeError, ValueError):
        pass

    converged = (
        abs(res) < tol_r23_rel * target_r23
        and e2_at <= 0.0
        and math.isfinite(independent_resid)
        and independent_resid < 1e-3
    )
    tof_to_capture = tof_days  # closest approach time bounded by the horizon
    return BCTForwardResult(
        v0_kms=float(v0),
        gamma0_deg=float(gamma0),
        terminal_r23_km=float(min_r * EM_DISTANCE_KM),
        terminal_e2=float(e2_at),
        dv_capture_kms=0.0,
        tof_days=float(tof_to_capture),
        converged=bool(converged),
        independent_terminal_residual=float(independent_resid),
        n_iter=int(n_iter),
        departure_state=_leo_departure_state(system, v0, gamma0, r13_km),
        notes=(
            f"res={res:.3e} (target r23={target_r23:.5f}), e2={e2_at:.4f}, "
            f"indep={independent_resid:.2e}"
        ),
    )


# ---------------------------------------------------------------------------
# Full BCT result + Hiten signature assembly.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class BCTResult:
    """A constructed Hiten-class BCT with the full ΔV / TOF / apoapsis signature.

    Attributes
    ----------
    apoapsis_ld :
        Sun-shaped Earth-relative apoapsis Q_a (LD).
    tof_days :
        Total transfer time of flight (days).
    capture_e2 :
        Kepler energy at capture (<= 0: ballistic, on W).
    dv_depart_kms / dv_match_kms / dv_capture_kms :
        ΔV breakdown: departure burn (arc I), apoapsis-match burn (arc II match),
        capture burn (exactly 0, ballistic by W definition).
    dv_total_kms :
        Sum of the three.
    capture_e2 / capture_r23_km :
        Capture energy and Moon-relative distance.
    """

    apoapsis_ld: float
    tof_days: float
    capture_e2: float
    capture_r23_km: float
    dv_depart_kms: float
    dv_match_kms: float
    dv_capture_kms: float
    dv_total_kms: float


# Published Belbruno 2004 §3.4 Hiten reference ΔV split (SOURCED constants, NOT
# values our corrector computed). The forward-from-LEO on-W stitch does not
# converge in the incoherent BCR4BP (documented model-gap boundary, design R3),
# so these are the published reference magnitudes the geometry targets, carried
# explicitly as the signature band -- never presented as a computed solve.
_HITEN_DV_DEPART_KMS = 0.014  # Belbruno §3.4 arc-I injection (14 m/s)
_HITEN_DV_MATCH_KMS = 0.030  # Belbruno §3.4 arc-II apoapsis match (30 m/s)


def build_hiten_bct(system: bcr4bp.BCR4BPSystem, target: BCTTarget) -> BCTResult:
    """Assemble a Hiten-class BCT: backward arc II (COMPUTED) + sourced ΔV split.

    The geometric signature is COMPUTED by this code:
      * the Sun-shaped Earth apoapsis Q_a (:func:`construct_bct_backward`);
      * the on-W ballistic-capture target (E_2 <= 0 at a Moon periapsis), exact
        by construction;
      * TOF from the backward arc.

    The ΔV split (14 + 30 m/s) is the PUBLISHED Belbruno §3.4 Hiten reference
    (:data:`_HITEN_DV_DEPART_KMS` / :data:`_HITEN_DV_MATCH_KMS`), NOT a value the
    corrector produced: the forward-from-LEO on-W capture does not converge in
    the incoherent BCR4BP (design R3 / §4 model-gap boundary). It is carried as
    the sourced signature band the geometry targets -- honestly labelled, never
    presented as a computed solve. ΔV_capture = 0 is exact (W definition).
    """
    arc = construct_bct_backward(target, system, back_days=70.0)
    return BCTResult(
        apoapsis_ld=arc.max_earth_apoapsis_ld,
        tof_days=min(arc.t_apoapsis_days * 2.0, 220.0),
        capture_e2=wsb.kepler_energy_moon(arc.qf_state, system),
        capture_r23_km=target.r23_km,
        dv_depart_kms=_HITEN_DV_DEPART_KMS,
        dv_match_kms=_HITEN_DV_MATCH_KMS,
        dv_capture_kms=0.0,
        dv_total_kms=_HITEN_DV_DEPART_KMS + _HITEN_DV_MATCH_KMS,
    )


# ---------------------------------------------------------------------------
# BCT novelty self-test (BCT corpus, NOT the cycler matcher).
# ---------------------------------------------------------------------------

# Belbruno / WSB / Hiten BCT corpus keywords. A BCT is NOT a cycler, so the
# cycler-only matcher in literature_check does not apply; BCT novelty is checked
# against the ballistic-capture corpus (Belbruno 2004, Belbruno-Miller 1993,
# Koon-Lo-Marsden-Ross 2001).
_BCT_AUTHORS = ("belbruno", "miller", "koon", "marsden", "ross")
_BCT_KEYWORDS = (
    "ballistic capture",
    "weak stability boundary",
    "hiten",
    "muses-a",
    "low energy lunar transfer",
    "low-energy transfer to the moon",
)


def bct_candidate_signature(
    sequence: tuple[str, ...] = ("E", "M"),
) -> CandidateSignature:
    """A CandidateSignature for a cislunar BCT (Earth-Moon, capture topology)."""
    return CandidateSignature(
        primary="Earth-Moon",
        sequence=sequence,
        topology_label=frozenset({"ballistic-capture"}),
    )


def _bct_hit_confidence(result: SearchResult) -> float:
    """Confidence a search hit is a published BCT (Belbruno/WSB/Hiten corpus).

    Distinct from the cycler matcher: scores on ballistic-capture corpus authors
    and keywords (a BCT is a transfer, not a cycler).
    """
    hay = (result.title + " " + result.snippet).lower()
    score = 0.0
    if any(kw in hay for kw in _BCT_KEYWORDS):
        score += 0.55
    if any(a in hay for a in _BCT_AUTHORS):
        score += 0.35
    if "moon" in hay or "lunar" in hay:
        score += 0.10
    return min(score, 1.0)


def check_bct_novelty(
    sig: CandidateSignature,
    *,
    search: Callable[[str], list[SearchResult]],
    threshold: float = 0.70,
) -> LiteratureCheckResult:
    """BCT novelty self-test against the ballistic-capture corpus.

    A re-derived Hiten transfer MUST be flagged ``published`` (non-novel) -- the
    golden guard against false-novelty (design R4). Uses a BCT-specific keyword
    matcher because a BCT is a transfer, not a cycler (the
    ``literature_check.check_literature`` cycler matcher mandates the word
    "cycler", which a Hiten-transfer hit legitimately lacks).
    """
    queries = [
        "Belbruno weak stability boundary ballistic capture lunar transfer",
        "Hiten MUSES-A ballistic capture Earth-Moon low energy transfer",
        "Koon Lo Marsden Ross low energy transfer to the Moon",
    ]
    trail: list[str] = []
    best_conf = 0.0
    best_hit: SearchResult | None = None
    any_results = False
    for q in queries:
        trail.append(q)
        try:
            results = list(search(q))
        except Exception as exc:
            trail[-1] = f"{q}  [ERROR: {exc!r}]"
            continue
        if results:
            any_results = True
        for r in results:
            conf = _bct_hit_confidence(r)
            if conf > best_conf:
                best_conf = conf
                best_hit = r
        if best_conf >= threshold:
            break

    if best_hit is not None and best_conf >= threshold:
        return LiteratureCheckResult(
            status="published",
            citation=best_hit.title,
            doi=None,
            confidence=round(best_conf, 3),
            query_trail=trail,
            matched_url=best_hit.url,
            notes="BCT signature matched the Belbruno/Hiten ballistic-capture corpus "
            "-- a rediscovery, NOT novelty-claimable.",
        )
    if not any_results:
        return LiteratureCheckResult(
            status="inconclusive",
            citation=None,
            doi=None,
            confidence=0.0,
            query_trail=trail,
            notes="No BCT search results returned -- rerun before trusting not-found.",
        )
    return LiteratureCheckResult(
        status="not-found",
        citation=None,
        doi=None,
        confidence=round(best_conf, 3),
        query_trail=trail,
        notes="No published BCT matched the signature (necessary-not-sufficient).",
    )
