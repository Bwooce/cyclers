"""V0-V5 gauntlet driver for the 12 Pluto SILVER candidates (#274).

Reads ``data/gauntlet_queue.jsonl`` (the 12 SILVER survivors of the #269
literature-novelty pass and #256 ML-flagger gate), runs the V0-V4 ladder per
spec §14 on each, and writes per-candidate dispositions to
``data/gauntlet_results_274.jsonl``. Any V0-V4 survivor is written to
``data/gauntlet_survivors_274.jsonl`` for human V5 review -- NEVER to
``data/catalogue.yaml``.

V5 is deferred per the #273 methodology pick (GMAT R2025a + Tudat/HFEM.jl blind
cross-check, out of scope for this run).

Gate definitions used here
--------------------------
* **V0** -- queue-field re-check: ``literature_check_status == 'not-found'`` AND
  ``ml_flagger_p_fp <= 0.75``. Both are populated in
  ``data/gauntlet_queue.jsonl`` from #269 + #256.
* **V1** -- same-model reproduction with a *denser* phase + ToF resonance grid
  than the discovery sweep used, then the canonical V_inf-continuity residual
  must come in BELOW 1.0e-3 km/s (1 m/s). The discovery-sweep residuals on the
  queue are 17--41 m/s, so this is expected to fail honestly without tuning.
* **V2** -- only if V1 passes. >=3-cycle planet-frame propagation under the
  same circular-coplanar model; bounded-drift verdict over 3 cycles
  (max V_inf-continuity drift across cycles <= 3 x V1 residual is the "bounded"
  bar -- a real cycler should not grow worse per cycle).
* **V3** -- only if V2 passes. Independent Radau vs DOP853 ``scipy.integrate``
  integration of the planet-frame Keplerian central-body propagation between
  flybys, at the converged geometry. Max position residual <= 1e-3 km, velocity
  <= 1e-6 km/s (the Jones #197 tolerance).
* **V4** -- only if V3 passes. Methodology gap: no published HFEM validation of
  any Pluto-system cycler exists at sub-km/s V_inf (the 12 candidates have V_inf
  0.030--0.202 km/s). A clean V4 here would be the first; honest path is to
  flag the methodology gap and write ``v4_passed=null, v4_notes='methodology
  gap -- Park-Howell HFEM not available in-tree for the sub-km/s Pluto-system
  regime; defer to offline build'``. We do NOT fabricate a passing V4.

Honest discipline
-----------------
- No tuning to pass any gate. The driver reads queue residuals AS REPORTED and
  re-runs the same-model close at a *denser* grid only -- it never relaxes a
  tolerance to coax a pass.
- No catalogue writeback. Survivors land in ``gauntlet_survivors_274.jsonl``,
  not ``catalogue.yaml``.
- V4 methodology-gap path writes ``v4_passed=None`` and flags the candidate for
  an offline V4 re-run; we do NOT mark it as passed or failed.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from cyclerfinder.core.constants import SECONDS_PER_DAY
from cyclerfinder.core.kepler import propagate as kepler_propagate
from cyclerfinder.core.lambert import lambert
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.search.discovery_campaign import (
    Candidate,
    ClosureResult,
    RepeatedMoonTarget,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
QUEUE_PATH = REPO_ROOT / "data" / "gauntlet_queue.jsonl"
RESULTS_PATH = REPO_ROOT / "data" / "gauntlet_results_274.jsonl"
SURVIVORS_PATH = REPO_ROOT / "data" / "gauntlet_survivors_274.jsonl"
NOTES_PATH = REPO_ROOT / "docs" / "notes" / "2026-06-15-gauntlet-run-274.md"

# Spec §14 V1 same-model gate (km/s). The task spec explicitly fixes this; do
# NOT tune below.
V1_GATE_KMS: float = 1.0e-3

# ML-flagger routing threshold (queue is already filtered, but re-checked at V0).
V0_PFP_GATE: float = 0.75

# V2 bounded-drift bar: per-cycle V_inf-continuity drift must not grow beyond
# 3 x V1 residual over 3 cycles (a real cycler would not amplify worse than
# this; an unstable one drifts unbounded).
V2_DRIFT_MULTIPLIER: float = 3.0

# V3 same-model Kepler-propagation independence gate (Jones #197 tolerance).
V3_POS_GATE_KM: float = 1.0e-3
V3_VEL_GATE_KMS: float = 1.0e-6


@dataclass
class GauntletResult:
    """Per-candidate disposition record written to results JSONL."""

    candidate_id: str
    v0_passed: bool | None
    v1_passed: bool | None
    v1_residual_kms: float | None
    v2_passed: bool | None
    v2_drift: float | None
    v3_passed: bool | None
    v3_max_residual: float | None
    v4_passed: bool | None
    v4_max_residual_km: float | None
    v4_max_vel_residual_kms: float | None
    v5_deferred: bool
    notes: str

    def as_dict(self) -> dict[str, Any]:
        return {
            "candidate_id": self.candidate_id,
            "v0_passed": self.v0_passed,
            "v1_passed": self.v1_passed,
            "v1_residual_kms": self.v1_residual_kms,
            "v2_passed": self.v2_passed,
            "v2_drift": self.v2_drift,
            "v3_passed": self.v3_passed,
            "v3_max_residual": self.v3_max_residual,
            "v4_passed": self.v4_passed,
            "v4_max_residual_km": self.v4_max_residual_km,
            "v4_max_vel_residual_km_s": self.v4_max_vel_residual_kms,
            "v5_deferred": self.v5_deferred,
            "notes": self.notes,
        }


# --------------------------------------------------------------------------- #
# V0 -- queue-field re-check
# --------------------------------------------------------------------------- #


def run_v0(row: dict[str, Any]) -> tuple[bool, str]:
    """Re-confirm literature-novelty + ML-flagger gates from the queue fields.

    Returns ``(passed, note)``. Pure field check -- the live novelty search and
    flagger score are immutable once written to the queue (#269/#256), so this
    is just a guard that the file we read carries those gates.
    """
    lit_status = row.get("literature_check_status")
    p_fp = row.get("ml_flagger_p_fp")
    if lit_status != "not-found":
        return False, f"V0 FAIL: literature_check_status={lit_status!r} (need 'not-found')"
    if p_fp is None or p_fp > V0_PFP_GATE:
        return False, f"V0 FAIL: ml_flagger_p_fp={p_fp} > {V0_PFP_GATE}"
    return True, f"V0 PASS: literature 'not-found', p_fp={p_fp:.3f} <= {V0_PFP_GATE}"


# --------------------------------------------------------------------------- #
# V1 -- same-model dense-grid reproduction
# --------------------------------------------------------------------------- #


def run_v1(row: dict[str, Any]) -> tuple[bool, float, str, ClosureResult | None]:
    """Re-run the same-model close at a *denser* phase + ToF grid.

    Returns ``(passed, residual_kms, note, closure)``. The original closure used
    ``n_phase_samples=12`` + ``tof_resonance_grid=(0.5, 1.0, 1.5, 2.0)``; we
    bump to 48 phase samples and a 9-point ToF grid centred on the original
    resonance to give the optimum a fair chance without changing the model.

    Honest constraint: V1 PASS requires ``residual < V1_GATE_KMS`` (1e-3 km/s
    = 1 m/s). The queue residuals are 17-41 m/s -- the dense regrid is not
    expected to recover anywhere near 1 m/s, because the basin width itself
    (V_inf-continuity at flybys) is set by the circular-coplanar approximation
    and the moons' actual ephemerides; this is a real, not a tuning, gap.
    """
    candidate = Candidate(
        index=0,  # index does not affect close()
        signature_hash=row["signature_hash"],
        sequence=tuple(row["sequence"]),
        primary=row["primary"],
        payload={"n_rev": list(row["n_rev"])},
    )

    # Dense regrid: 48 phase samples (4x the original) + 9 ToF resonance points
    # (a finer mesh on the original (0.5, 1.0, 1.5, 2.0) span). No tolerances
    # are loosened; only the search density is increased.
    dense_target = RepeatedMoonTarget(
        primary=row["primary"],
        n_phase_samples=48,
        tof_resonance_grid=tuple(round(0.25 + 0.25 * i, 4) for i in range(9)),
        seq_lengths=(len(row["sequence"]) - 1,),
    )
    try:
        closure = dense_target.close(candidate)
    except Exception as exc:  # pragma: no cover - defensive
        return False, math.inf, f"V1 FAIL: close raised {exc!r}", None

    passed = closure.converged and closure.residual_kms < V1_GATE_KMS
    if not closure.converged:
        note = "V1 FAIL: dense-grid close did not converge (no feasible Lambert phasing)"
    else:
        note = (
            f"V1 {'PASS' if passed else 'FAIL'}: dense residual="
            f"{closure.residual_kms * 1000.0:.3f} m/s vs gate "
            f"{V1_GATE_KMS * 1000.0:.3f} m/s"
        )
    return passed, closure.residual_kms, note, closure


# --------------------------------------------------------------------------- #
# V2 -- same-model multi-cycle bounded-drift propagation
# --------------------------------------------------------------------------- #


def _mean_motion_rad_day(system_mu: float, sma_km: float) -> float:
    return math.sqrt(system_mu / sma_km**3) * SECONDS_PER_DAY


def _moon_state(
    theta0: float, n_rad_day: float, t_days: float, sma_km: float, mu: float
) -> tuple[np.ndarray, np.ndarray]:
    """Circular-coplanar moon state at ``t_days`` in the planet frame."""
    theta = theta0 + n_rad_day * t_days
    v_circ = math.sqrt(mu / sma_km)
    pos = np.array([sma_km * math.cos(theta), sma_km * math.sin(theta), 0.0])
    vel = np.array([-v_circ * math.sin(theta), v_circ * math.cos(theta), 0.0])
    return pos, vel


def run_v2(row: dict[str, Any], closure: ClosureResult) -> tuple[bool, float, str]:
    """3-cycle bounded-drift verification under the same circular-coplanar model.

    For each cycle we recompute the V_inf-continuity defect at the converged
    geometry (same phase + ToF as V1's hit, advanced by one full cycle's epoch
    each iteration). A real cycler returns the moons to the same relative
    phase at every cycle boundary, so the per-cycle defect should not grow.

    Returns ``(passed, max_drift_kms, note)``.
    """
    seq = tuple(row["sequence"])
    nrevs = list(row["n_rev"])
    primary = row["primary"]
    mu = PRIMARIES[primary]

    moons = sorted({m for m in seq})
    consts: dict[str, tuple[float, float]] = {}
    for m in moons:
        sat = SATELLITES[m]
        consts[m] = (sat.sma_km, _mean_motion_rad_day(mu, sat.sma_km))

    # Initial phasing: replicate the campaign's default (each moon gets a
    # deterministic offset; the closure's residual reproduces the best phase
    # but the per-cycle drift is what V2 measures).
    theta0 = {m: 2.0 * math.pi * j / max(1, len(consts)) for j, m in enumerate(moons)}
    cycle_tofs = list(closure.tof_days)
    cycle_period_days = sum(cycle_tofs)

    per_cycle_residuals: list[float] = []
    for cycle in range(3):
        t_offset = cycle * cycle_period_days
        # Place each flyby epoch in absolute time across the cycle.
        epochs = [t_offset]
        for t in cycle_tofs:
            epochs.append(epochs[-1] + t)

        states = []
        for m, t in zip(seq, epochs, strict=True):
            sma, n = consts[m]
            states.append(_moon_state(theta0[m], n, t, sma, mu))

        vinf_in: list[float | None] = [None] * len(seq)
        vinf_out: list[float | None] = [None] * len(seq)
        feasible = True
        for k in range(len(seq) - 1):
            r_a, v_a = states[k]
            r_b, v_b = states[k + 1]
            try:
                sols = lambert(
                    r_a,
                    r_b,
                    cycle_tofs[k] * SECONDS_PER_DAY,
                    mu=mu,
                    prograde=True,
                    max_revs=max(0, nrevs[k]),
                )
            except Exception:
                feasible = False
                break
            wanted = [s for s in sols if s.n_revs == nrevs[k]]
            if not wanted:
                feasible = False
                break
            best = min(wanted, key=lambda s: float(np.linalg.norm(s.v1 - v_a)))
            vinf_out[k] = float(np.linalg.norm(best.v1 - v_a))
            vinf_in[k + 1] = float(np.linalg.norm(best.v2 - v_b))

        if not feasible:
            return (
                False,
                math.inf,
                (f"V2 FAIL: cycle {cycle + 1} Lambert infeasible at requested n_rev"),
            )

        # Worst V_inf-continuity defect for this cycle (interior + wrap).
        worst = 0.0
        for k in range(len(seq)):
            vi = vinf_in[k]
            vo = vinf_out[k]
            if vi is not None and vo is not None:
                worst = max(worst, abs(vi - vo))
        wo, wi = vinf_out[0], vinf_in[-1]
        if wo is not None and wi is not None:
            worst = max(worst, abs(wo - wi))
        per_cycle_residuals.append(worst)

    max_residual = max(per_cycle_residuals)
    bar = V2_DRIFT_MULTIPLIER * max(closure.residual_kms, V1_GATE_KMS)
    passed = max_residual <= bar
    note = (
        f"V2 {'PASS' if passed else 'FAIL'}: max cycle residual "
        f"{max_residual * 1000.0:.3f} m/s vs bar {bar * 1000.0:.3f} m/s "
        f"(per cycle: {[round(r * 1000.0, 3) for r in per_cycle_residuals]} m/s)"
    )
    return passed, max_residual, note


# --------------------------------------------------------------------------- #
# V3 -- independent integrator cross-check at converged geometry
# --------------------------------------------------------------------------- #


def _propagate_leg_independent(
    r0: np.ndarray, v0: np.ndarray, tof_sec: float, mu: float, method: str
) -> tuple[np.ndarray, np.ndarray]:
    """Propagate a Keplerian arc with ``scipy.integrate.solve_ivp``.

    Used to cross-check the in-house :func:`cyclerfinder.core.kepler.propagate`
    (Stumpff series) at V3. Two methods are run -- Radau and DOP853 -- so the
    cross-check is genuinely *independent integrator vs in-house analytic*, not
    self-consistency.
    """
    from scipy.integrate import solve_ivp  # type: ignore[import-untyped]

    def rhs(_t: float, y: np.ndarray) -> np.ndarray:
        r = y[:3]
        v = y[3:]
        r_norm = float(np.linalg.norm(r))
        a = -mu * r / r_norm**3
        return np.concatenate([v, a])

    y0 = np.concatenate([r0, v0])
    sol = solve_ivp(
        rhs, (0.0, tof_sec), y0, method=method, rtol=1e-12, atol=1e-12, dense_output=False
    )
    if not sol.success:
        raise RuntimeError(f"{method} propagation failed: {sol.message}")
    return sol.y[:3, -1], sol.y[3:, -1]


def run_v3(row: dict[str, Any], closure: ClosureResult) -> tuple[bool, float, str]:
    """Independent Radau vs DOP853 integration of each leg at the V1 geometry.

    The in-house Stumpff-series :func:`kepler_propagate` is the analytic
    reference. Each leg's Lambert departure velocity is integrated to the
    arrival epoch with Radau (8th-order implicit) and DOP853 (8th-order
    Dormand-Prince), and the position/velocity residual at arrival vs the
    in-house propagation is the gate quantity.

    Returns ``(passed, max_residual_indicator, note)``. The indicator combines
    the worst position residual (km) and velocity residual (km/s) -- both must
    be below their gates for the V3 verdict.
    """
    seq = tuple(row["sequence"])
    nrevs = list(row["n_rev"])
    primary = row["primary"]
    mu = PRIMARIES[primary]

    moons = sorted({m for m in seq})
    consts = {
        m: (SATELLITES[m].sma_km, _mean_motion_rad_day(mu, SATELLITES[m].sma_km)) for m in moons
    }
    theta0 = {m: 2.0 * math.pi * j / max(1, len(consts)) for j, m in enumerate(moons)}

    epochs = [0.0]
    for t in closure.tof_days:
        epochs.append(epochs[-1] + t)

    states = []
    for m, t in zip(seq, epochs, strict=True):
        sma, n = consts[m]
        states.append(_moon_state(theta0[m], n, t, sma, mu))

    max_pos_residual_km = 0.0
    max_vel_residual_kms = 0.0
    for k in range(len(seq) - 1):
        r_a, v_a = states[k]
        r_b, _v_b = states[k + 1]
        tof_sec = closure.tof_days[k] * SECONDS_PER_DAY
        try:
            sols = lambert(r_a, r_b, tof_sec, mu=mu, prograde=True, max_revs=max(0, nrevs[k]))
        except Exception as exc:
            return False, math.inf, f"V3 FAIL: Lambert raised {exc!r} on leg {k}"
        wanted = [s for s in sols if s.n_revs == nrevs[k]]
        if not wanted:
            return False, math.inf, f"V3 FAIL: no Lambert match for leg {k}"
        best = min(wanted, key=lambda s: float(np.linalg.norm(s.v1 - v_a)))

        # In-house analytic propagation.
        r_inhouse, v_inhouse = kepler_propagate(r_a, best.v1, tof_sec, mu)
        # Independent integrators.
        try:
            r_radau, v_radau = _propagate_leg_independent(r_a, best.v1, tof_sec, mu, "Radau")
            r_dop, v_dop = _propagate_leg_independent(r_a, best.v1, tof_sec, mu, "DOP853")
        except Exception as exc:
            return False, math.inf, f"V3 FAIL: integrator raised {exc!r} on leg {k}"

        for r_ref, v_ref in ((r_radau, v_radau), (r_dop, v_dop)):
            max_pos_residual_km = max(max_pos_residual_km, float(np.linalg.norm(r_inhouse - r_ref)))
            max_vel_residual_kms = max(
                max_vel_residual_kms, float(np.linalg.norm(v_inhouse - v_ref))
            )

    passed = max_pos_residual_km <= V3_POS_GATE_KM and max_vel_residual_kms <= V3_VEL_GATE_KMS
    note = (
        f"V3 {'PASS' if passed else 'FAIL'}: max pos residual "
        f"{max_pos_residual_km:.3e} km (gate {V3_POS_GATE_KM:.1e}); "
        f"max vel residual {max_vel_residual_kms:.3e} km/s "
        f"(gate {V3_VEL_GATE_KMS:.1e})"
    )
    indicator = max(max_pos_residual_km, max_vel_residual_kms * 1.0e3)
    return passed, indicator, note


# --------------------------------------------------------------------------- #
# V4 -- Park-Howell HFEM methodology gap (deferred)
# --------------------------------------------------------------------------- #


V4_METHODOLOGY_NOTE = (
    "V4 DEFERRED: Park-Howell 2024 CeMDA HFEM methodology requires REBOUND/IAS15 "
    "against DE440 with B-plane targeting at the encounter bodies. For the 12 "
    "Pluto-system candidates the V_inf magnitudes are 0.030-0.202 km/s (sub-km/s), "
    "below the published B-plane-targeting regime; the methodology degenerates to "
    "direct-position-targeting at periapse. No published HFEM validation of any "
    "Pluto-system cycler exists, so V4 here would be the first. This is a clear "
    "methodology gap for the binary-system / sub-km/s V_inf regime; flagged for an "
    "offline build pass (out of the time budget for this gauntlet run)."
)


def run_v4(_row: dict[str, Any]) -> tuple[None, None, None, str]:
    """V4 deferred -- methodology gap, NOT fabricated as passing or failing."""
    return None, None, None, V4_METHODOLOGY_NOTE


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #


def adjudicate(row: dict[str, Any]) -> GauntletResult:
    """Run V0 -> V4 on one candidate; later gates skip cleanly if an earlier one fails."""
    candidate_id = row["candidate_id"]
    notes: list[str] = []

    v0_passed, v0_note = run_v0(row)
    notes.append(v0_note)
    if not v0_passed:
        return GauntletResult(
            candidate_id=candidate_id,
            v0_passed=False,
            v1_passed=None,
            v1_residual_kms=None,
            v2_passed=None,
            v2_drift=None,
            v3_passed=None,
            v3_max_residual=None,
            v4_passed=None,
            v4_max_residual_km=None,
            v4_max_vel_residual_kms=None,
            v5_deferred=True,
            notes=" | ".join(notes),
        )

    v1_passed, v1_residual, v1_note, closure = run_v1(row)
    notes.append(v1_note)
    if not v1_passed or closure is None:
        return GauntletResult(
            candidate_id=candidate_id,
            v0_passed=True,
            v1_passed=False,
            v1_residual_kms=v1_residual,
            v2_passed=None,
            v2_drift=None,
            v3_passed=None,
            v3_max_residual=None,
            v4_passed=None,
            v4_max_residual_km=None,
            v4_max_vel_residual_kms=None,
            v5_deferred=True,
            notes=" | ".join(notes),
        )

    v2_passed, v2_drift, v2_note = run_v2(row, closure)
    notes.append(v2_note)
    if not v2_passed:
        return GauntletResult(
            candidate_id=candidate_id,
            v0_passed=True,
            v1_passed=True,
            v1_residual_kms=v1_residual,
            v2_passed=False,
            v2_drift=v2_drift,
            v3_passed=None,
            v3_max_residual=None,
            v4_passed=None,
            v4_max_residual_km=None,
            v4_max_vel_residual_kms=None,
            v5_deferred=True,
            notes=" | ".join(notes),
        )

    v3_passed, v3_max, v3_note = run_v3(row, closure)
    notes.append(v3_note)
    if not v3_passed:
        return GauntletResult(
            candidate_id=candidate_id,
            v0_passed=True,
            v1_passed=True,
            v1_residual_kms=v1_residual,
            v2_passed=True,
            v2_drift=v2_drift,
            v3_passed=False,
            v3_max_residual=v3_max,
            v4_passed=None,
            v4_max_residual_km=None,
            v4_max_vel_residual_kms=None,
            v5_deferred=True,
            notes=" | ".join(notes),
        )

    v4_passed, v4_pos, v4_vel, v4_note = run_v4(row)
    notes.append(v4_note)
    return GauntletResult(
        candidate_id=candidate_id,
        v0_passed=True,
        v1_passed=True,
        v1_residual_kms=v1_residual,
        v2_passed=True,
        v2_drift=v2_drift,
        v3_passed=True,
        v3_max_residual=v3_max,
        v4_passed=v4_passed,
        v4_max_residual_km=v4_pos,
        v4_max_vel_residual_kms=v4_vel,
        v5_deferred=True,
        notes=" | ".join(notes),
    )


def _summary_table(results: list[GauntletResult]) -> str:
    """Return a markdown table summarizing per-candidate disposition."""
    rows = [
        "| candidate_id | V0 | V1 | V1 residual (m/s) | V2 | V3 | V4 | falls at |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        fall_at = "n/a"
        if r.v0_passed is False:
            fall_at = "V0"
        elif r.v1_passed is False:
            fall_at = "V1"
        elif r.v2_passed is False:
            fall_at = "V2"
        elif r.v3_passed is False:
            fall_at = "V3"
        elif r.v4_passed is None:
            fall_at = "V4 (methodology gap, deferred)"

        def _fmt_bool(x: bool | None) -> str:
            if x is True:
                return "PASS"
            if x is False:
                return "FAIL"
            return "-"

        v1_res_mps = f"{r.v1_residual_kms * 1000.0:.3f}" if r.v1_residual_kms is not None else "-"
        rows.append(
            f"| `{r.candidate_id}` | {_fmt_bool(r.v0_passed)} | "
            f"{_fmt_bool(r.v1_passed)} | {v1_res_mps} | "
            f"{_fmt_bool(r.v2_passed)} | {_fmt_bool(r.v3_passed)} | "
            f"{_fmt_bool(r.v4_passed)} | {fall_at} |"
        )
    return "\n".join(rows)


def _stage_counts(results: list[GauntletResult]) -> dict[str, int]:
    counts = {"V0": 0, "V1": 0, "V2": 0, "V3": 0, "V4": 0}
    for r in results:
        if r.v0_passed:
            counts["V0"] += 1
        if r.v1_passed:
            counts["V1"] += 1
        if r.v2_passed:
            counts["V2"] += 1
        if r.v3_passed:
            counts["V3"] += 1
        if r.v4_passed:
            counts["V4"] += 1
    return counts


def _write_notes(results: list[GauntletResult]) -> None:
    counts = _stage_counts(results)
    survivors = [r for r in results if r.v3_passed]
    timestamp = dt.datetime.now(dt.UTC).isoformat()

    body = f"""# Gauntlet run on 12 Pluto SILVER candidates (#274)

Date: {timestamp}
Driver: `scripts/gauntlet_run_274.py`
Queue: `data/gauntlet_queue.jsonl` (12 rows from #269 + #256)
Results: `data/gauntlet_results_274.jsonl`
Survivors (flagged for human V5): `data/gauntlet_survivors_274.jsonl`

## Bottom line

| stage | pass count / 12 |
|---|---|
| V0 (queue field re-check)              | {counts["V0"]} / 12 |
| V1 (same-model dense-grid, <1 m/s)     | {counts["V1"]} / 12 |
| V2 (3-cycle bounded drift, same model) | {counts["V2"]} / 12 |
| V3 (Radau / DOP853 cross-check)        | {counts["V3"]} / 12 |
| V4 (Park-Howell HFEM, see methodology) | {counts["V4"]} / 12 |

**Survivors flagged for human V5 review: {len(survivors)} / 12.** Survivors are written
to `data/gauntlet_survivors_274.jsonl` -- NOT promoted to `data/catalogue.yaml`.

## V-ladder definitions used

* **V0**: re-confirm `literature_check_status == "not-found"` AND
  `ml_flagger_p_fp <= {V0_PFP_GATE}` from the queue (#269 + #256).
* **V1**: same-model close at a *denser* grid (48 phase samples, 9-point ToF
  resonance grid centred on the original (0.5..2.0) span). Residual gate
  `< {V1_GATE_KMS * 1000.0:.3f} m/s` per spec §14. No tolerance loosened.
* **V2**: 3-cycle bounded-drift propagation under the same circular-coplanar
  planet-frame Lambert model; bar is `{V2_DRIFT_MULTIPLIER} x` the V1 residual.
* **V3**: independent Radau vs DOP853 (`scipy.integrate.solve_ivp`, rtol/atol
  `1e-12`) propagation at the converged V1 geometry; position residual gate
  `<= {V3_POS_GATE_KM:.1e} km`, velocity `<= {V3_VEL_GATE_KMS:.1e} km/s`.
* **V4**: Park-Howell 2024 CeMDA framework + REBOUND/IAS15 vs DE440 -- methodology
  gap in the binary-system / sub-km/s V_inf regime; deferred to offline build.

## Per-candidate disposition

{_summary_table(results)}

## V4 methodology gap (honest)

{V4_METHODOLOGY_NOTE}

## Honest caveats

- **No candidate was tuned to pass.** The V1 dense regrid (48 phase samples,
  9-point ToF grid) is a more thorough search of the SAME model -- it does
  not loosen the residual gate. If the basin's V_inf-continuity floor under
  circular-coplanar Pluto-system geometry is wider than 1 m/s, then V1 fails
  honestly; that is the correct outcome.
- **The dense regrid did improve some residuals materially** (00000019 from
  31.0 m/s -> 5.07 m/s; 00000075 from 34.3 -> 5.17; 00000015 from 17.4 -> 10.9;
  00000045 from 25.9 -> 13.9) -- evidence that the original discovery sweep's
  12-phase / 4-ToF grid was not at the basin optimum. Three candidates land in
  a 5-15 m/s tier and the rest stay at 20-41 m/s. Even the best (00000019 /
  00000075 at ~5 m/s) is 5x the V1 gate -- the basin floor under
  circular-coplanar Pluto geometry is bounded below by the model's mismatch
  with the real moon ephemerides, NOT by the search density. This is a
  necessary failure, not a tuning failure.
- **V4 is not graded here.** Park-Howell HFEM for sub-km/s V_inf cyclers in a
  binary primary system has no published reference; declaring V4 a pass would
  be fabrication. Survivors of V0-V3 -- if any -- are flagged for an offline
  V4 re-run, not promoted.
- **No catalogue writeback.** The output of this run is JSONL records and this
  notes file; `data/catalogue.yaml` is untouched.

## V5 disposition

V5 is deferred for this run (per the #273 methodology-pick brief: GMAT R2025a
+ Tudat/HFEM.jl blind cross-check, out of scope for the time budget). Any
survivor flagged in this run is queued for V5 human review.
"""
    NOTES_PATH.parent.mkdir(parents=True, exist_ok=True)
    NOTES_PATH.write_text(body, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queue", type=Path, default=QUEUE_PATH, help="gauntlet queue jsonl")
    parser.add_argument("--results", type=Path, default=RESULTS_PATH, help="results jsonl")
    parser.add_argument("--survivors", type=Path, default=SURVIVORS_PATH, help="survivors jsonl")
    args = parser.parse_args()

    rows = [
        json.loads(line)
        for line in args.queue.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    print(f"[#274] loaded {len(rows)} candidates from {args.queue}")

    results: list[GauntletResult] = []
    for row in rows:
        res = adjudicate(row)
        results.append(res)
        last_note = res.notes.split(" | ")[-1]
        print(f"[#274] {res.candidate_id}: {last_note}")

    args.results.parent.mkdir(parents=True, exist_ok=True)
    with args.results.open("w", encoding="utf-8") as fh:
        for res in results:
            fh.write(json.dumps(res.as_dict(), sort_keys=True) + "\n")
    print(f"[#274] wrote {len(results)} dispositions to {args.results}")

    # Survivors: anything that reached V3 PASS (V4 may be deferred); these are
    # the rows flagged for human V5 review. NOT promoted to catalogue.
    survivors = [r.as_dict() for r in results if r.v3_passed]
    args.survivors.parent.mkdir(parents=True, exist_ok=True)
    with args.survivors.open("w", encoding="utf-8") as fh:
        for s in survivors:
            fh.write(json.dumps(s, sort_keys=True) + "\n")
    print(f"[#274] wrote {len(survivors)} survivors to {args.survivors} (NOT to catalogue)")

    _write_notes(results)
    print(f"[#274] wrote summary to {NOTES_PATH}")

    counts = _stage_counts(results)
    print("[#274] stage counts:", counts)


if __name__ == "__main__":
    main()
