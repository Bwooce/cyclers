"""CR3BP Saturnian midsize-moon periodic-orbit discovery (Task 6, plan 2026-06-10).

Seeds planar Lyapunov guesses near the L1/L2 collinear libration points for
Saturn/Mimas, Saturn/Enceladus, Saturn/Tethys; runs the STM single-shooting
corrector; applies a DEGENERACY GATE to every converged solution; and routes
survivors to:
  - data/cr3bp_silver.jsonl  -- gate-surviving, independently cross-checked orbits
                                (method-tagged SILVER, no sourced anchor; JSONL,
                                review-gated).
  - data/empty_regions.jsonl -- a method-versioned EMPTY note for every pair with
                                no surviving family.

Degeneracy gate (v2; the v1 run was contaminated by degenerate "convergences"):
  1. Equilibrium rejection -- a libration point trivially satisfies X(T)=X(0) for
     ANY period, and the min-norm Newton step happily converges onto L1/L2/L4/L5.
     Reject if max |v| over the propagated period < 1e-6 nondim or the position
     amplitude max|r(t)-r(0)| < 1e-6 nondim.
  2. Period floor -- reject period < 0.1 nondim (period collapse is another
     trivial closure).
  3. Dedup -- survivors in the same pair whose state0 positions agree within 1e-9
     are one orbit; keep the better closure_residual.
  4. Independent cross-check -- crosscheck_periodic (Radau re-propagation, Jacobi
     drift) must pass; dJacobi is recorded in the SILVER entry.

Usage:
    uv run python scripts/cr3bp_moontour_run.py
    uv run python scripts/cr3bp_moontour_run.py --report /tmp/cr3bp_moontour.txt

IMPORTANT: NO catalogue writeback.  SILVER entries require human review before
promotion.  The EMPTY note records the search extent so a future
higher-capability sweep can supersede it.
"""

from __future__ import annotations

import argparse
import contextlib
import json
import subprocess
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp

# ---------------------------------------------------------------------------
# Method version + git SHA (best-effort; empty string if not available)
# ---------------------------------------------------------------------------
METHOD_TAG = "cr3bp-lyapunov-corrector-v2"

_GIT_SHA: str = ""
with contextlib.suppress(Exception):
    _GIT_SHA = (
        subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            stderr=subprocess.DEVNULL,
        )
        .decode()
        .strip()
    )

# ---------------------------------------------------------------------------
# Output paths
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent
SILVER_QUEUE_PATH = REPO_ROOT / "data" / "cr3bp_silver.jsonl"
EMPTY_REGIONS_PATH = REPO_ROOT / "data" / "empty_regions.jsonl"
RESULTS_NOTE_PATH = REPO_ROOT / "docs" / "notes" / "2026-06-10-cr3bp-moontour-results.md"

# ---------------------------------------------------------------------------
# Saturnian midsize-moon targets
# ---------------------------------------------------------------------------
TARGETS: list[tuple[str, str]] = [
    ("Saturn", "Mimas"),
    ("Saturn", "Enceladus"),
    ("Saturn", "Tethys"),
]

# Lyapunov-family seed amplitudes, as fractions of the libration-point distance
# gamma from the secondary.  Each amplitude lands on a different family member
# (the full-state min-norm corrector does not pin energy).
#
# Basin note (measured 2026-06-10, this commit): full-period single shooting on
# these strongly unstable orbits (|eigenvalue| ~ e^{2.5 T}) has a TINY convergence
# basin -- seeds with Ax beyond ~1e-5 nondim (~0.003 gamma) diverge or get captured
# by the L4/L5 equilibria (which the degeneracy gate rejects).  Amplitude-scaled
# continuation from a converged member does not extend the basin.  The top fracs
# here deliberately sit at/past the basin edge so non-convergence and gate
# rejections are exercised and counted.
_AMPLITUDE_FRACS: list[float] = [0.0005, 0.001, 0.002, 0.003, 0.005]

# Degeneracy-gate thresholds (nondimensional).
MAX_SPEED_FLOOR_ND = 1e-6  # equilibrium gate: max |v| over the arc
AMPLITUDE_FLOOR_ND = 1e-6  # equilibrium gate: max |r(t) - r(0)| over the arc
PERIOD_FLOOR_ND = 0.1  # period-collapse gate
DEDUP_POS_TOL_ND = 1e-9  # state0 position agreement => same orbit


def _collinear_x(mu: float, point: str) -> float:
    """x-coordinate of L1 ("L1") or L2 ("L2"): Newton on dU/dx = 0, Hill-radius start."""
    r_h = float((mu / 3.0) ** (1.0 / 3.0))
    x = 1.0 - mu - r_h if point == "L1" else 1.0 - mu + r_h
    for _ in range(100):
        s1 = x + mu
        s2 = x - 1.0 + mu
        f = x - (1.0 - mu) * s1 / abs(s1) ** 3 - mu * s2 / abs(s2) ** 3
        fp = 1.0 + 2.0 * (1.0 - mu) / abs(s1) ** 3 + 2.0 * mu / abs(s2) ** 3
        dx = f / fp
        x -= dx
        if abs(dx) < 1e-15:
            break
    return x


def _lyapunov_seed(mu: float, x_lib: float, amp_x: float) -> tuple[NDArray[np.float64], float]:
    """Linear planar Lyapunov initial condition + period guess at a collinear point.

    Linearised in-plane dynamics at a collinear point (c2 = (1-mu)/r1^3 + mu/r2^3):
        xi'' - 2 eta' - (1 + 2 c2) xi = 0
        eta'' + 2 xi'  - (1 - c2) eta = 0
    The centre eigenfrequency nu satisfies
        nu^2 = (2 - c2 + sqrt(9 c2^2 - 8 c2)) / 2,
    and the periodic solution xi = Ax cos(nu t), eta = -k Ax sin(nu t) with
        k = (nu^2 + 1 + 2 c2) / (2 nu)
    gives the seed state [x_lib + Ax, 0, 0, 0, -k nu Ax, 0] and period 2 pi / nu.
    """
    s1 = abs(x_lib + mu)
    s2 = abs(x_lib - 1.0 + mu)
    c2 = (1.0 - mu) / s1**3 + mu / s2**3
    nu = float(np.sqrt((2.0 - c2 + float(np.sqrt(9.0 * c2**2 - 8.0 * c2))) / 2.0))
    k = float((nu**2 + 1.0 + 2.0 * c2) / (2.0 * nu))
    state0 = np.array([x_lib + amp_x, 0.0, 0.0, 0.0, -k * nu * amp_x, 0.0])
    return state0, float(2.0 * np.pi / nu)


@dataclass(frozen=True)
class GateMetrics:
    """Arc diagnostics used by the equilibrium gate."""

    max_speed_nd: float  # max |v| over one propagated period
    amplitude_nd: float  # max |r(t) - r(0)| over one propagated period


def degeneracy_gate(system: cr3bp.CR3BPSystem, orbit: cp.PeriodicOrbit) -> tuple[str, GateMetrics]:
    """Apply the period-floor and equilibrium gates to a converged orbit.

    Returns ("", metrics) on pass, or (reason, metrics) with reason in
    {"period_floor", "equilibrium"} on rejection.
    """
    t_eval = np.linspace(0.0, orbit.period, 256)
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, orbit.period),
        np.asarray(orbit.state0, float),
        args=(system.mu,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        t_eval=t_eval,
    )
    states = sol.y
    speeds = np.linalg.norm(states[3:6, :], axis=0)
    disp = np.linalg.norm(states[0:3, :] - states[0:3, :1], axis=0)
    metrics = GateMetrics(
        max_speed_nd=float(np.max(speeds)),
        amplitude_nd=float(np.max(disp)),
    )
    if orbit.period < PERIOD_FLOOR_ND:
        return "period_floor", metrics
    if metrics.max_speed_nd < MAX_SPEED_FLOOR_ND or metrics.amplitude_nd < AMPLITUDE_FLOOR_ND:
        return "equilibrium", metrics
    return "", metrics


@dataclass(frozen=True)
class SilverCandidate:
    """A converged orbit that survived every gate, with its provenance."""

    orbit: cp.PeriodicOrbit
    libration_point: str
    amplitude_frac: float
    metrics: GateMetrics
    crosscheck_djacobi: float


@dataclass
class PairOutcome:
    """Per-pair tallies: seeds tried, converged, rejected-by-gate, SILVER survivors."""

    primary: str
    secondary: str
    mu: float
    l_km: float
    t_s: float
    seeds_tried: int = 0
    converged: int = 0
    rejected_period_floor: int = 0
    rejected_equilibrium: int = 0
    rejected_duplicate: int = 0
    rejected_crosscheck: int = 0
    silver: list[SilverCandidate] = field(default_factory=list)

    @property
    def rejected_total(self) -> int:
        return (
            self.rejected_period_floor
            + self.rejected_equilibrium
            + self.rejected_duplicate
            + self.rejected_crosscheck
        )


def run_moon_pair(primary: str, secondary: str) -> PairOutcome:
    """Seed, correct, gate, dedup, and cross-check one Saturn/moon pair."""
    system = cr3bp.cr3bp_system(primary, secondary)
    out = PairOutcome(
        primary=primary,
        secondary=secondary,
        mu=system.mu,
        l_km=system.l_km,
        t_s=system.t_s,
    )

    # Seed + correct + per-orbit gates.
    candidates: list[tuple[cp.PeriodicOrbit, str, float, GateMetrics]] = []
    for point in ("L1", "L2"):
        x_lib = _collinear_x(system.mu, point)
        gamma = abs(x_lib - (1.0 - system.mu))
        for frac in _AMPLITUDE_FRACS:
            out.seeds_tried += 1
            s0, t_guess = _lyapunov_seed(system.mu, x_lib, frac * gamma)
            res = cp.correct_periodic(system, s0, t_guess, max_iter=60)
            if not res.converged:
                continue
            out.converged += 1
            reason, metrics = degeneracy_gate(system, res)
            if reason == "period_floor":
                out.rejected_period_floor += 1
                continue
            if reason == "equilibrium":
                out.rejected_equilibrium += 1
                continue
            candidates.append((res, point, frac, metrics))

    # Dedup: same state0 position within tolerance => one orbit; keep best residual.
    candidates.sort(key=lambda c: c[0].closure_residual)
    kept: list[tuple[cp.PeriodicOrbit, str, float, GateMetrics]] = []
    for cand in candidates:
        pos = np.asarray(cand[0].state0[:3], float)
        is_dup = any(
            float(np.linalg.norm(pos - np.asarray(k[0].state0[:3], float))) < DEDUP_POS_TOL_ND
            for k in kept
        )
        if is_dup:
            out.rejected_duplicate += 1
        else:
            kept.append(cand)

    # Independent cross-check (Radau re-propagation + Jacobi drift).
    for orbit, point, frac, metrics in kept:
        ok, dj = cp.crosscheck_periodic(system, orbit)
        if not ok:
            out.rejected_crosscheck += 1
            continue
        out.silver.append(
            SilverCandidate(
                orbit=orbit,
                libration_point=point,
                amplitude_frac=frac,
                metrics=metrics,
                crosscheck_djacobi=dj,
            )
        )
    return out


def _now_iso() -> str:
    return datetime.now(UTC).isoformat()


def append_silver(cand: SilverCandidate, primary: str, secondary: str) -> None:
    """Append a SILVER CR3BP periodic-orbit record to data/cr3bp_silver.jsonl."""
    orbit = cand.orbit
    cid = (
        f"cr3bp-lyapunov-{secondary.lower()}-{cand.libration_point.lower()}"
        f"-J{orbit.jacobi:.6f}-T{orbit.period:.4f}"
    )
    record: dict[str, object] = {
        "candidate_id": cid,
        "verdict_tier": "SILVER",
        "method": METHOD_TAG,
        "git_sha": _GIT_SHA,
        "primary": primary,
        "secondary": secondary,
        "libration_point": cand.libration_point,
        "seed_amplitude_frac": cand.amplitude_frac,
        "jacobi_constant": orbit.jacobi,
        "period_nd": orbit.period,
        "state0_nd": orbit.state0.tolist(),
        "closure_residual": orbit.closure_residual,
        "converged": orbit.converged,
        "max_speed_nd": cand.metrics.max_speed_nd,
        "amplitude_nd": cand.metrics.amplitude_nd,
        "crosscheck_djacobi": cand.crosscheck_djacobi,
        "gates_passed": [
            "equilibrium-rejection(max|v|>=1e-6,amp>=1e-6)",
            "period-floor(T>=0.1)",
            "dedup(state0-pos-1e-9)",
            "crosscheck-radau(closure<1e-8,dJ<1e-8)",
        ],
        "t_added": _now_iso(),
        "sourced_anchor": None,
        "notes": (
            "Planar Lyapunov periodic orbit near a collinear point; seeded from the "
            "linearised centre eigensolution and corrected by STM single-shooting. "
            "Degeneracy-gated and Radau cross-checked (v2; the v1 run was contaminated "
            "by equilibrium/period-collapse pseudo-convergences and was discarded). "
            "No independent source; SILVER pending human review. NO catalogue writeback."
        ),
    }
    SILVER_QUEUE_PATH.parent.mkdir(parents=True, exist_ok=True)
    with SILVER_QUEUE_PATH.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(record, ensure_ascii=True) + "\n")


def append_empty(outcome: PairOutcome) -> None:
    """Append a method-versioned EMPTY note to data/empty_regions.jsonl."""
    record: dict[str, object] = {
        "region_id": (
            f"cr3bp-lyapunov-{outcome.secondary.lower()}-{datetime.now(UTC).strftime('%Y-%m-%d')}"
        ),
        "family": (f"CR3BP rotating-frame periodic orbits ({outcome.primary}/{outcome.secondary})"),
        "centre": outcome.secondary,
        "method_capability": {
            "genome": "CR3BP planar Lyapunov near L1/L2 (STM single-shooting)",
            "corrector": "correct_periodic (tol=1e-10, max_iter=60)",
            "method": METHOD_TAG,
            "capability_tags": ["cr3bp", "planar", "lyapunov", "collinear"],
            "git_sha": _GIT_SHA,
        },
        "search_extent": {
            "n_seeds": outcome.seeds_tried,
            "seed_points": "L1, L2",
            "amplitude_fractions": _AMPLITUDE_FRACS,
            "seeding": "linearised centre eigensolution (xi=Ax cos, eta=-k Ax sin)",
        },
        "result": {
            "converged_orbits": outcome.converged,
            "rejected_period_floor": outcome.rejected_period_floor,
            "rejected_equilibrium": outcome.rejected_equilibrium,
            "rejected_duplicate": outcome.rejected_duplicate,
            "rejected_crosscheck": outcome.rejected_crosscheck,
            "silver_survivors": 0,
        },
        "verdict": "EMPTY -- no gate-surviving planar Lyapunov orbit for this pair",
        "interpretation": (
            "Either no seed converged or every convergence was degenerate "
            "(equilibrium / period collapse / duplicate / cross-check failure). "
            "A higher-capability search (halo orbits, vertical Lyapunov, "
            "continuation) may find orbits; this search only covers the planar "
            "small-amplitude Lyapunov family."
        ),
        "run": {
            "date": datetime.now(UTC).strftime("%Y-%m-%d"),
            "git_sha": _GIT_SHA,
        },
    }
    EMPTY_REGIONS_PATH.parent.mkdir(parents=True, exist_ok=True)
    with EMPTY_REGIONS_PATH.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(record, ensure_ascii=True) + "\n")


def build_runlog(outcomes: list[PairOutcome], elapsed_s: float) -> str:
    """Build the plain-text per-pair report (also embedded verbatim in the note)."""
    ts = datetime.now(tz=UTC).isoformat()
    lines: list[str] = [
        f"cr3bp_moontour_run  {ts}  method={METHOD_TAG}  git={_GIT_SHA or '?'}"
        f"  elapsed={elapsed_s:.1f}s",
        "=" * 78,
    ]
    for o in outcomes:
        lines.append(f"pair: {o.primary}/{o.secondary}  mu={o.mu:.3e}")
        lines.append(
            f"  seeds tried        : {o.seeds_tried}  (L1+L2 x amplitude fracs {_AMPLITUDE_FRACS})"
        )
        lines.append(f"  converged          : {o.converged}")
        lines.append(
            f"  rejected degenerate: {o.rejected_total}"
            f"  [equilibrium={o.rejected_equilibrium}"
            f" period_floor={o.rejected_period_floor}"
            f" duplicate={o.rejected_duplicate}"
            f" crosscheck={o.rejected_crosscheck}]"
        )
        lines.append(
            f"  crosscheck passed  : {len(o.silver)} of "
            f"{len(o.silver) + o.rejected_crosscheck} checked"
        )
        if o.silver:
            lines.append(f"  SILVER written     : {len(o.silver)}")
            for c in o.silver:
                orb = c.orbit
                lines.append(
                    f"    {c.libration_point} Ax/gamma={c.amplitude_frac:<4}"
                    f"  J={orb.jacobi:.10f}  T={orb.period:.8f}"
                    f"  resid={orb.closure_residual:.2e}"
                    f"  max|v|={c.metrics.max_speed_nd:.3e}"
                    f"  amp={c.metrics.amplitude_nd:.3e}"
                    f"  dJ={c.crosscheck_djacobi:.2e}"
                )
        else:
            lines.append("  SILVER written     : 0 -> EMPTY record appended")
        lines.append("")
    lines.append("NO catalogue writeback performed.")
    return "\n".join(lines)


def build_note(outcomes: list[PairOutcome], runlog: str) -> str:
    """Build the markdown results note (verbatim runlog + caveats)."""
    ts = datetime.now(tz=UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    total_silver = sum(len(o.silver) for o in outcomes)
    lines: list[str] = [
        "# CR3BP Saturnian Midsize-Moon Discovery Results (degeneracy-gated, v2)",
        "",
        f"**Run timestamp:** {ts}",
        "**Script:** `scripts/cr3bp_moontour_run.py`",
        f"**Method:** `{METHOD_TAG}`  (git `{_GIT_SHA or '?'}`)",
        "",
        "**Status:** NO catalogue writeback.  SILVER entries in `data/cr3bp_silver.jsonl`",
        "(review-gated).  EMPTY notes (if any) in `data/empty_regions.jsonl`.",
        "",
        "## Why v2",
        "",
        "The v1 run (`cr3bp-lyapunov-corrector-v1`, 11 entries, never committed) was",
        'contaminated: the min-norm full-state corrector "converged" onto the L1/L2/L4/L5',
        "libration points themselves (rotating-frame speed ~1e-12, closed for any period),",
        "onto a period-collapse solution (T=4.5e-5), and onto the same equilibrium twice",
        "with different reported periods.  That output was deleted and regenerated from",
        "scratch under the v2 degeneracy gate:",
        "",
        "1. **Equilibrium rejection** -- max |v| over the propagated period >= 1e-6 nondim",
        "   AND position amplitude max|r(t)-r(0)| >= 1e-6 nondim.",
        "2. **Period floor** -- period >= 0.1 nondim.",
        "3. **Dedup** -- state0 positions within 1e-9 are one orbit (best residual kept).",
        "4. **Independent cross-check** -- `crosscheck_periodic` (Radau, vs the corrector's",
        "   DOP853) must re-close within 1e-8 with Jacobi drift < 1e-8; dJ recorded.",
        "",
        "v2 also fixes the seeding: the v1 draft's linear frequency formula was wrong",
        "(it substituted `1+2c2` where the standard formula takes `c2`, with a sign slip),",
        "and its constant-Jacobi seeding placed initial velocities off the Lyapunov family.",
        "v2 seeds the linearised centre eigensolution at Newton-solved collinear points:",
        "`x0 = x_L + Ax`, `vy0 = -k nu Ax`, `T_guess = 2 pi / nu`.",
        "",
        "## Run report (verbatim)",
        "",
        "```text",
        runlog,
        "```",
        "",
        "## Per-pair summary",
        "",
        "| Pair | Seeds | Converged | Rejected (equilibrium/period/dup/xcheck) | SILVER |",
        "|---|---|---|---|---|",
    ]
    for o in outcomes:
        rej = (
            f"{o.rejected_equilibrium}/{o.rejected_period_floor}"
            f"/{o.rejected_duplicate}/{o.rejected_crosscheck}"
        )
        outcome_s = str(len(o.silver)) if o.silver else "0 (EMPTY)"
        lines.append(
            f"| {o.primary}/{o.secondary} | {o.seeds_tried} | {o.converged} | {rej} | {outcome_s} |"
        )
    lines += [
        "",
        f"**Total SILVER candidates:** {total_silver}",
        "",
        "## Honest caveat",
        "",
        "Small collinear-point (L1/L2) planar Lyapunov families are mathematically",
        "guaranteed to exist in any CR3BP (Lyapunov centre theorem) and are NOT novel",
        "discoveries in the literature sense -- such families are tabulated for many",
        "systems.  Their value here is exercising and validating the seed -> STM-corrector",
        "-> degeneracy-gate -> independent-crosscheck pipeline on new (Saturn, midsize-moon)",
        "mass ratios.  They route to review (SILVER) regardless; promotion past SILVER",
        "would require a sourced anchor.  **NO writeback to `data/catalogue.yaml`.**",
        "",
        "All survivors are SMALL-amplitude members (Ax ~ 1e-6..1e-5 nondim, i.e. of order",
        "a kilometre at these moons): the full-period single-shooting corrector's",
        "convergence basin on these strongly unstable orbits is tiny, and larger-amplitude",
        "seeds either diverge or are captured by the L4/L5 equilibria (gate-rejected).",
        "Extending the families to useful amplitudes needs a symmetry-exploiting",
        "half-period or multiple-shooting corrector -- a method-capability gap, recorded",
        "here so a future sweep can supersede this one.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="CR3BP Saturnian midsize-moon discovery (no catalogue writeback)"
    )
    parser.add_argument("--report", type=Path, default=None, help="Write plain-text runlog here")
    args = parser.parse_args()

    t0 = datetime.now(tz=UTC)
    print("cr3bp_moontour_run: Saturnian midsize-moon CR3BP periodic-orbit discovery")
    print(f"  Targets: {', '.join(f'{p}/{s}' for p, s in TARGETS)}")
    print(f"  Method:  {METHOD_TAG}  (git {_GIT_SHA or '(unknown)'})")
    print()

    outcomes: list[PairOutcome] = []
    for primary, secondary in TARGETS:
        print(f"  [{primary}/{secondary}] seeding, correcting, gating...", flush=True)
        outcome = run_moon_pair(primary, secondary)
        outcomes.append(outcome)
        print(
            f"    seeds={outcome.seeds_tried} converged={outcome.converged}"
            f" rejected={outcome.rejected_total} silver={len(outcome.silver)}"
        )
        if outcome.silver:
            for cand in outcome.silver:
                append_silver(cand, primary, secondary)
            print(f"    -> appended to {SILVER_QUEUE_PATH}")
        else:
            append_empty(outcome)
            print(f"    -> EMPTY note appended to {EMPTY_REGIONS_PATH}")
        print()

    elapsed = (datetime.now(tz=UTC) - t0).total_seconds()
    runlog = build_runlog(outcomes, elapsed)
    print(runlog)

    RESULTS_NOTE_PATH.parent.mkdir(parents=True, exist_ok=True)
    RESULTS_NOTE_PATH.write_text(build_note(outcomes, runlog), encoding="utf-8")
    print(f"\nResults note written: {RESULTS_NOTE_PATH}")

    if args.report is not None:
        Path(args.report).write_text(runlog + "\n", encoding="utf-8")
        print(f"Runlog written:       {args.report}")

    print("\nNO catalogue writeback performed.")


if __name__ == "__main__":
    main()
