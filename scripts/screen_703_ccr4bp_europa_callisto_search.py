"""#703 -- CCR4BP manifold globalization + heteroclinic search: Jupiter Europa-Callisto.

Genuine discovery attempt (not a capability-proof exercise -- that was #694's
JEG positive control, `scripts/screen_694_ccr4bp_heteroclinic_search.py`,
whose structure this script mirrors, alongside #695's/#696's/#701's own
discovery-attempt scripts, `scripts/screen_695_ccr4bp_io_europa.py` /
`scripts/screen_696_ccr4bp_io_ganymede_search.py` /
`scripts/screen_701_ccr4bp_umbriel_titania_search.py`). #693's/#700's
last remaining CCR4BP candidate: mu_gan=5.667e-5 (comparable order to the
already-validated JEG system's own mu_gan), novelty-cleared by #700's own
deep literature check (zero "Callisto" mentions in either Kumar-group source
PDF, five further live-search queries with no direct Europa-Callisto
CCR4BP/PCCFBP hit), but with a real, flagged scientific-motivation caveat: the
Europa:Callisto period ratio (~4.73) has NO clean low-integer
commensurability, unlike every other CCR4BP pair this project has built --
every seed orbit tried converges to essentially the SAME near-degenerate
near-circular family rather than a physically-motivated resonant orbit (see
`tests/search/test_ccr4bp_torus_europa_callisto.py`'s own module docstring
for the full re-derivation of this finding, not just a citation of #700's own
quick scratch check).

Every module used (core.ccr4bp, search.variational_ccr4bp_torus,
search.ccr4bp_whisker, search.ccr4bp_manifold_globalize,
search.ccr4bp_heteroclinic_search) is #689-#694's code, reused UNMODIFIED,
INCLUDING #702's ghost_guard ref_vec-anchoring fix (already present in the
committed module -- this task does not need its own separate anchor-bug
investigation). The only new code is the system constructor
(core.ccr4bp_europa_callisto, #703) and this driver.

Base resonant orbit choice (documented in full in
tests/search/test_ccr4bp_torus_europa_callisto.py's module docstring):
spacecraft:Europa = 4:1 (interior), following #696's own Io-Ganymede
precedent for a near-4:1-labeled system -- the naive exterior
"moon-ratio-matching" 1:5 reading does not hard-fail under Callisto's
physical forcing (unlike Io-Ganymede's own naive 1:4 choice) but has only a
thin ~1.1x-Hill-scale margin, an order of magnitude narrower than the chosen
interior reading's ~22x margin.

IMPORTANT physical-unit note -- a COINCIDENCE here, not a correction
------------------------------------------------------------------------
search.ccr4bp_heteroclinic_search hardcodes a module-level `_L_KM = 671_100.0`
(= EUROPA's SMA) and `_v_unit_km_s` hardcodes GM_Europa + GM_Jupiter -- both
baked in for the JEG positive control, NOT parameterized per-system. Every
prior non-JEG discovery-attempt script (#695 Io-Europa, #696 Io-Ganymede,
#701 Umbriel-Titania) had to independently rescale every km-denominated field
because their base moon was NOT Europa. For Europa-Callisto, the base moon
IS Europa (unchanged from JEG), so `core.ccr4bp_europa_callisto.L_KM` equals
`671_100.0` exactly and `v_unit_km_s()` equals the module's own hardcoded
value exactly -- VERIFIED explicitly below (not just assumed from the module
docstring), so `module_native_*` and `corrected_*` fields are reported side
by side for consistency with the other three discovery-attempt scripts, and
should agree to numerical precision.

Two-tier genuineness: raw gate vs ROBUST pass (mandatory self-skepticism)
---------------------------------------------------------------------------
Mirrors #701's own mandatory-skepticism discipline exactly: a candidate that
merely crosses the raw `#694` ghost-guard's `integrator_delta_km < 1.0` km
numeric cutoff without a comfortable margin is NOT trusted as genuine here --
this script computes a stricter `robust_genuine` flag
(`corrected_integrator_delta_km < _ROBUST_INTEGRATOR_MARGIN_KM = 0.1` km,
still ~5 orders of magnitude looser than JEG's own actual ~6.5e-7 km margin)
and reports `best_robust_genuine_connection_corrected` as the trustworthy
headline result, separately from the raw-gate-only
`best_genuine_connection_corrected`.

No catalogue writeback -- a discovery-attempt script, not a vetted result.

Run:  uv run python scripts/screen_703_ccr4bp_europa_callisto_search.py
Outputs -> data/found/703_ccr4bp_europa_callisto_search/result.json
"""

from __future__ import annotations

import json
import sys
import time
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import cyclerfinder.core.ccr4bp_europa_callisto as ecal  # noqa: E402
import cyclerfinder.core.cr3bp as cr3bp  # noqa: E402
import cyclerfinder.search.ccr4bp_heteroclinic_search as hs  # noqa: E402
import cyclerfinder.search.ccr4bp_manifold_globalize as mg  # noqa: E402
import cyclerfinder.search.variational_ccr4bp_torus as vt  # noqa: E402
from cyclerfinder.genome.composed_moon_map import resonance_semimajor  # noqa: E402
from cyclerfinder.search.variational_ccr4bp_torus import evaluate_torus_state  # noqa: E402

OUT_DIR = ROOT / "data" / "found" / "703_ccr4bp_europa_callisto_search"

# Correct physical units for THIS system (base moon = Europa, SAME as JEG).
_L_KM_EUROPA = ecal.L_KM  # 671,100 km, expected to equal hs._L_KM exactly
_V_UNIT_EUROPA_KM_S = ecal.v_unit_km_s()

# See #701's own module for the full derivation of this margin choice.
_ROBUST_INTEGRATOR_MARGIN_KM = 0.1


def _resonant_symmetric_orbit(
    mu: float, p_sc: int, q_moon: int, *, max_iter: int = 80, tol: float = 1e-12, cap: float = 0.05
) -> tuple[np.ndarray, float, float]:
    """Identical test-only scaffolding to #690/#691/#694/#695/#696/#701's own scripts/tests."""
    a = resonance_semimajor(p_sc, q_moon)
    period = 2.0 * np.pi * q_moon
    th = 0.5 * period
    x0 = a - mu
    vy0 = float(np.sqrt((1.0 - mu) / a)) - x0
    res = np.inf
    for k in range(max_iter):
        s0 = np.array([x0, 0.0, 0.0, 0.0, vy0, 0.0])
        y42 = np.concatenate([s0, np.eye(6).reshape(36)])
        sol = solve_ivp(
            cr3bp.cr3bp_stm_eom, (0.0, th), y42, args=(mu,), method="DOP853", rtol=1e-12, atol=1e-12
        )
        sf = sol.y[:, -1]
        phi = sf[6:].reshape(6, 6)
        g = np.array([sf[1], sf[3]])
        res = float(np.linalg.norm(g))
        if res < tol:
            break
        jac = np.array([[phi[1, 0], phi[1, 4]], [phi[3, 0], phi[3, 4]]])
        dz = np.linalg.solve(jac, -g)
        dz = dz * (0.3 if k < 8 else 1.0)
        norm = float(np.linalg.norm(dz))
        if norm > cap:
            dz = dz / norm * cap
        x0 += dz[0]
        vy0 += dz[1]
    return np.array([x0, 0.0, 0.0, 0.0, vy0, 0.0]), period, res


def _log(msg: str, t0: float) -> None:
    print(f"[{time.time() - t0:7.1f}s] {msg}", flush=True)


def _corrected_gaps(
    refined: hs.RefinedConnection,
) -> dict[str, float]:
    """Recompute pos/vel gap in km using Europa's OWN physical units -- for
    THIS system this is expected to numerically equal the module-native
    value (base moon IS Europa), computed independently here rather than
    trusted from the module docstring alone."""
    pos_gap_km = float(np.linalg.norm(refined.state_u[:3] - refined.state_s[:3])) * _L_KM_EUROPA
    vel_gap_km_s = (
        float(np.linalg.norm(refined.state_u[3:] - refined.state_s[3:])) * _V_UNIT_EUROPA_KM_S
    )
    return {"pos_gap_km": pos_gap_km, "vel_gap_km_s": vel_gap_km_s}


def _corrected_off_torus_km(
    torus_u: vt.CCR4BPTorusVariationalResult,
    refined: hs.RefinedConnection,
    theta1_section: float,
) -> float:
    """Replicate ghost_guard's off-torus distance formula with Europa's own
    km scale (see ghost_guard's own docstring for the formula this mirrors)."""
    theta1_u_at_t = theta1_section + torus_u.omega1 * refined.t_u
    theta2_u_at_t = refined.theta2_u + torus_u.omega2 * refined.t_u
    torus_pt_u = evaluate_torus_state(torus_u, theta1_u_at_t, theta2_u_at_t)
    planar_u = refined.state_u[[0, 1, 3, 4]]
    return float(np.linalg.norm(torus_pt_u - planar_u)) * _L_KM_EUROPA


def _corrected_radau_check(
    torus_u: vt.CCR4BPTorusVariationalResult,
    torus_s: vt.CCR4BPTorusVariationalResult,
    refined: hs.RefinedConnection,
    *,
    lobe_sign_u: float,
    lobe_sign_s: float,
    theta1_section: float,
    n_segments_dir: int,
    rtol: float,
    atol: float,
) -> dict[str, float]:
    """Re-run #694's own private Radau re-propagation helper (the same
    independent-integrator check ghost_guard performs internally) and rescale
    the resulting position gap with Europa's own km unit. Calls
    hs._radau_manifold_state directly -- the same private-helper access
    pattern this project's OWN test module already uses for hs._L_KM -- not
    a modification of #694's module.

    `#702`: anchor the Radau re-check with the SEED-anchored ref_vec
    refine_candidate actually threaded through the optimization (carried on
    ``refined.ref_vec_u``/``ref_vec_s``), NOT a fresh anchor re-derived at the
    FINAL converged theta2 -- the fix is already present in the module this
    script imports; this helper reuses it, mirroring #701's own updated
    `_corrected_radau_check`."""
    ref_vec_u = refined.ref_vec_u
    ref_vec_s = refined.ref_vec_s
    su_radau = hs._radau_manifold_state(
        torus_u,
        "unstable",
        theta1_section,
        refined.theta2_u,
        refined.t_u,
        eps=mg.DEFAULT_EPS,
        lobe_sign=lobe_sign_u,
        n_segments_dir=n_segments_dir,
        rtol=rtol,
        atol=atol,
        ref_vec=ref_vec_u,
    )
    ss_radau = hs._radau_manifold_state(
        torus_s,
        "stable",
        theta1_section,
        refined.theta2_s,
        refined.t_s,
        eps=mg.DEFAULT_EPS,
        lobe_sign=lobe_sign_s,
        n_segments_dir=n_segments_dir,
        rtol=rtol,
        atol=atol,
        ref_vec=ref_vec_s,
    )
    if su_radau is None or ss_radau is None:
        return {"radau_pos_gap_km": float("nan"), "integrator_delta_km": float("inf")}
    radau_pos_gap_km = float(np.linalg.norm(su_radau[:3] - ss_radau[:3])) * _L_KM_EUROPA
    corrected_dop853_pos_gap_km = (
        float(np.linalg.norm(refined.state_u[:3] - refined.state_s[:3])) * _L_KM_EUROPA
    )
    integrator_delta_km = abs(radau_pos_gap_km - corrected_dop853_pos_gap_km)
    return {"radau_pos_gap_km": radau_pos_gap_km, "integrator_delta_km": integrator_delta_km}


def main() -> dict[str, Any]:
    t0 = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    system = ecal.jupiter_europa_callisto_default()
    _log(
        f"system: mu={system.mu:.6e} mu_gan={system.mu_gan:.6e} a_gan={system.a_gan:.4f} "
        f"omega_gan={system.omega_gan:.4f}",
        t0,
    )
    # Verify the "coincidence, not a fix" unit claim explicitly before using it.
    unit_coincidence_verified = (
        abs(_L_KM_EUROPA - hs._L_KM) < 1e-9
        and abs(_V_UNIT_EUROPA_KM_S - hs._v_unit_km_s(system)) < 1e-9
    )
    _log(
        f"physical-unit coincidence check: L_KM={_L_KM_EUROPA} vs hs._L_KM={hs._L_KM}, "
        f"v_unit={_V_UNIT_EUROPA_KM_S:.6f} vs hs._v_unit_km_s={hs._v_unit_km_s(system):.6f}, "
        f"verified={unit_coincidence_verified}",
        t0,
    )
    assert unit_coincidence_verified, "expected L_KM/v_unit to coincide with hs's own JEG constants"

    s0, period, res = _resonant_symmetric_orbit(system.mu, 4, 1)
    assert res < 1e-10, f"resonant orbit did not converge: {res:.2e}"
    _log(f"base 4:1 (interior) resonant orbit converged, perp residual {res:.2e}", t0)

    torus = vt.discover_ccr4bp_torus_from_resonant_orbit(
        system,
        s0,
        period,
        n1=1,
        n2=20,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    _log(
        f"physical-mass 4:1 CCR4BP torus built: residual_rms={torus.residual_rms:.3e}, "
        f"closure_residual={torus.closure_residual:.3e}, rho_strob={torus.rho_strob:.4f}, "
        f"period={torus.period:.4f} TU",
        t0,
    )

    # ------------------------------------------------------------------
    # Stage 1: globalize both branches, both lobes.
    #
    # t_max_periods=6.0 (NOT the #694/#701 default of 2.0): a preliminary
    # cheap coarse-only diagnostic (no refine step) at t_max_periods=2.0
    # found every coarse candidate confined to a small (gap_planar~1e-3
    # nondim) neighbourhood, all of which refine to a tiny off_torus_km
    # (40-254 km, well under the 1000 km "genuinely departed" gate) --
    # every one a trivial near-departure pseudo-match. Extending the same
    # cheap coarse-only diagnostic to t_max_periods=6.0 found MUCH tighter
    # coarse candidates (gap_planar ~1e-6 nondim, three orders of magnitude
    # smaller) at elapsed times spanning 0.6-5.9 torus periods -- i.e. the
    # genuine intersection region (if any) needs a substantially longer
    # flow-time search horizon than #694's original JEG-tuned default,
    # consistent with this system's much stronger one-period amplification
    # (|lam_u|~175 vs JEG's ~6-13) requiring more periods for the manifold
    # to genuinely separate from and then re-approach the torus. This is a
    # caller-level search-parameter tuning decision (per
    # ccr4bp_heteroclinic_search's own module docstring: mesh-refinement/
    # search-horizon choices are explicitly caller-driven, not baked into
    # the module), not a modification of any #689-#694 module.
    #
    # A further targeted follow-up (refine_candidate + ghost_guard, not just
    # the cheap coarse-only diagnostic) at t_max_periods=6.0 found the
    # search's best candidate so far climbing to off_torus_km~647 km
    # (residual_norm~4.3e-7, integrator agreement ~8.8e-7 km) -- MUCH closer
    # to the 1000 km gate than the original 2.0-period search's 40-254 km,
    # though not yet across it. n_theta2=40/n_time=200 (rather than a finer
    # 60/220 grid tried in an earlier coarse-only scan) is the exact grid
    # that follow-up used and is kept here for direct continuity with that
    # result; n_candidates raised to 8 (from #694's/#701's own default 5)
    # for more thorough coverage of this longer, richer search horizon.
    n_theta2, n_time, t_max_periods, n_segments_dir = 40, 200, 6.0, 24
    tubes: dict[str, mg.ManifoldTube] = {}
    for branch in ("unstable", "stable"):
        for lobe in (1.0, -1.0):
            key = f"{branch}_{'+' if lobe > 0 else '-'}"
            tubes[key] = mg.globalize_manifold_tube(
                torus,
                branch,
                n_theta2=n_theta2,
                t_max_periods=t_max_periods,
                n_time=n_time,
                n_segments_dir=n_segments_dir,
                lobe_sign=lobe,
            )
            n_valid = int(tubes[key].valid.sum())
            _log(f"globalized {key}: {n_valid}/{n_theta2} phases valid", t0)

    # ------------------------------------------------------------------
    # Stage 2: coarse search + refine + ghost-guard, all 4 lobe-pair combos.
    # ------------------------------------------------------------------
    combo_results: list[dict[str, Any]] = []
    best: dict[str, Any] | None = None
    best_robust: dict[str, Any] | None = None
    for u_lobe in (1.0, -1.0):
        for s_lobe in (1.0, -1.0):
            u_key = f"unstable_{'+' if u_lobe > 0 else '-'}"
            s_key = f"stable_{'+' if s_lobe > 0 else '-'}"
            candidates = hs.coarse_candidates(
                tubes[u_key], tubes[s_key], n_candidates=8, t_min_frac=0.05
            )
            _log(f"combo ({u_key}, {s_key}): {len(candidates)} coarse candidates", t0)
            refined_list: list[dict[str, Any]] = []
            for cand in candidates:
                refined = hs.refine_candidate(
                    torus,
                    torus,
                    cand,
                    lobe_sign_u=u_lobe,
                    lobe_sign_s=s_lobe,
                    n_segments_dir=32,
                )
                if refined is None:
                    continue
                guard = hs.ghost_guard(
                    torus, torus, refined, lobe_sign_u=u_lobe, lobe_sign_s=s_lobe, n_segments_dir=32
                )
                corrected = _corrected_gaps(refined)
                corrected_off_torus = _corrected_off_torus_km(torus, refined, 0.0)
                corrected_radau = _corrected_radau_check(
                    torus,
                    torus,
                    refined,
                    lobe_sign_u=u_lobe,
                    lobe_sign_s=s_lobe,
                    theta1_section=0.0,
                    n_segments_dir=32,
                    rtol=1e-13,
                    atol=1e-13,
                )
                corrected_radau_consistent = corrected_radau["integrator_delta_km"] < 1.0
                corrected_genuine = corrected_radau_consistent and corrected_off_torus >= 1000.0
                robust_genuine = (
                    corrected_genuine
                    and corrected_radau["integrator_delta_km"] < _ROBUST_INTEGRATOR_MARGIN_KM
                )
                entry = {
                    "seed": asdict(cand),
                    # Module-native fields (Europa-SMA-scaled, JEG's own hardcoded constant).
                    "module_native_pos_gap_km": refined.pos_gap_km,
                    "module_native_vel_gap_km_s": refined.vel_gap_km_s,
                    "module_native_guard_off_torus_km": guard.off_torus_km,
                    "module_native_guard_integrator_delta_km": guard.integrator_delta_km,
                    "module_native_guard_genuine": guard.genuine,
                    # Independently recomputed (Europa-SMA-scaled) physical values --
                    # expected (and verified above) to equal the module-native ones.
                    "corrected_pos_gap_km": corrected["pos_gap_km"],
                    "corrected_vel_gap_km_s": corrected["vel_gap_km_s"],
                    "corrected_off_torus_km": corrected_off_torus,
                    "corrected_radau_pos_gap_km": corrected_radau["radau_pos_gap_km"],
                    "corrected_integrator_delta_km": corrected_radau["integrator_delta_km"],
                    "corrected_genuine": corrected_genuine,
                    "robust_genuine": robust_genuine,
                    "residual_norm": refined.residual_norm,
                    "converged": refined.converged,
                    "theta2_u": refined.theta2_u,
                    "t_u": refined.t_u,
                    "theta2_s": refined.theta2_s,
                    "t_s": refined.t_s,
                    "guard_quasi_jacobi_gap": guard.quasi_jacobi_gap,
                    "guard_notes": guard.notes,
                }
                refined_list.append(entry)
                _log(
                    f"  refined: pos_gap(corrected)={corrected['pos_gap_km']:.2f} km, "
                    f"vel_gap(corrected)={corrected['vel_gap_km_s'] * 1000:.3f} m/s, "
                    f"off_torus(corrected)={corrected_off_torus:.1f} km, "
                    f"integ_delta(corrected)={corrected_radau['integrator_delta_km']:.4g} km, "
                    f"genuine(module/corrected/robust)={guard.genuine}/{corrected_genuine}/{robust_genuine}",
                    t0,
                )
                is_better = best is None or corrected["pos_gap_km"] < best["corrected_pos_gap_km"]
                if corrected_genuine and is_better:
                    best = {**entry, "u_lobe": u_lobe, "s_lobe": s_lobe}
                is_better_robust = (
                    best_robust is None
                    or corrected["pos_gap_km"] < best_robust["corrected_pos_gap_km"]
                )
                if robust_genuine and is_better_robust:
                    best_robust = {**entry, "u_lobe": u_lobe, "s_lobe": s_lobe}
            combo_results.append(
                {
                    "u_lobe": u_lobe,
                    "s_lobe": s_lobe,
                    "n_coarse_candidates": len(candidates),
                    "refined": refined_list,
                }
            )

    if best is not None:
        _log(
            f"best genuine-by-raw-gate (corrected) connection: "
            f"pos_gap={best['corrected_pos_gap_km']:.3g} km, "
            f"integ_delta={best['corrected_integrator_delta_km']:.4g} km",
            t0,
        )
    else:
        _log("NO connection clears the raw ghost-guard gate in any lobe combo", t0)
    if best_robust is not None:
        _log(
            f"best ROBUSTLY-genuine (corrected, integ_delta < {_ROBUST_INTEGRATOR_MARGIN_KM} km) "
            f"connection: pos_gap={best_robust['corrected_pos_gap_km']:.3f} km, "
            f"vel_gap={best_robust['corrected_vel_gap_km_s'] * 1000:.3f} m/s",
            t0,
        )
    else:
        _log(
            "NO robustly-genuine connection found -- every raw-gate pass is marginal/suspect",
            t0,
        )

    # ------------------------------------------------------------------
    # Seed-perturbation stability re-check on the raw-gate "best" IF it is
    # not also the robust best (mirrors #701's own mandatory-skepticism check).
    # ------------------------------------------------------------------
    ghost_sensitivity_check: dict[str, Any] | None = None
    if best is not None and (best_robust is None or best is not best_robust):
        seed_cand = hs.ManifoldCandidate(**best["seed"])
        reseed_cand = hs.ManifoldCandidate(
            theta2_u=best["theta2_u"],
            t_u=best["t_u"],
            theta2_s=best["theta2_s"],
            t_s=best["t_s"],
            gap_planar=0.0,
        )
        seed_deltas: dict[str, float] = {}
        for label, cand in [
            ("original_coarse_seed", seed_cand),
            ("converged_point_as_seed", reseed_cand),
        ]:
            refined_variant = hs.refine_candidate(
                torus,
                torus,
                cand,
                lobe_sign_u=best["u_lobe"],
                lobe_sign_s=best["s_lobe"],
                n_segments_dir=32,
                rtol=1e-13,
                atol=1e-13,
            )
            if refined_variant is None:
                continue
            corrected_radau_variant = _corrected_radau_check(
                torus,
                torus,
                refined_variant,
                lobe_sign_u=best["u_lobe"],
                lobe_sign_s=best["s_lobe"],
                theta1_section=0.0,
                n_segments_dir=32,
                rtol=1e-13,
                atol=1e-13,
            )
            seed_deltas[label] = corrected_radau_variant["integrator_delta_km"]
        fragile = len(seed_deltas) >= 2 and max(seed_deltas.values()) > 10.0 * max(
            min(seed_deltas.values()), 1e-6
        )
        ghost_sensitivity_check = {
            "target": "raw-gate best (not robust)",
            "integrator_delta_km_by_seed_variant": seed_deltas,
            "fragile_under_seed_perturbation": fragile,
            "interpretation": (
                "Two DIFFERENT starting seeds for the SAME nominal candidate give WILDLY "
                "different independent-integrator (Radau) agreement -- a numerically fragile "
                "coincidence, not a robust dynamical intersection. Treated as an UNCONFIRMED "
                "suspected ghost artifact, not a genuine connection, despite nominally "
                "clearing the raw 1.0 km gate on one of the two seedings."
                if fragile
                else "Both seed variants agree closely on integrator consistency -- no fragility "
                "detected under this check."
            ),
        }
        _log(
            f"seed-perturbation stability re-check on raw-gate best: "
            f"deltas={seed_deltas}, fragile={fragile}",
            t0,
        )

    # ------------------------------------------------------------------
    # Stage 3: mesh-refinement stability check.
    # ------------------------------------------------------------------
    mesh_check: dict[str, Any] | None = None
    mesh_target = best_robust if best_robust is not None else best
    mesh_target_kind = "robust_genuine" if best_robust is not None else "raw_gate_only"
    if mesh_target is None:
        all_entries = [(c["u_lobe"], c["s_lobe"], e) for c in combo_results for e in c["refined"]]
        if all_entries:
            u_lobe, s_lobe, entry = min(all_entries, key=lambda t: t[2]["corrected_pos_gap_km"])
            mesh_target = {**entry, "u_lobe": u_lobe, "s_lobe": s_lobe}
            mesh_target_kind = "closest_overall_not_genuine"

    if mesh_target is not None:
        n_theta2_dense, n_time_dense = 120, 300
        tube_u_dense = mg.globalize_manifold_tube(
            torus,
            "unstable",
            n_theta2=n_theta2_dense,
            t_max_periods=t_max_periods,
            n_time=n_time_dense,
            n_segments_dir=n_segments_dir,
            lobe_sign=mesh_target["u_lobe"],
        )
        tube_s_dense = mg.globalize_manifold_tube(
            torus,
            "stable",
            n_theta2=n_theta2_dense,
            t_max_periods=t_max_periods,
            n_time=n_time_dense,
            n_segments_dir=n_segments_dir,
            lobe_sign=mesh_target["s_lobe"],
        )
        dense_candidates = hs.coarse_candidates(
            tube_u_dense, tube_s_dense, n_candidates=5, t_min_frac=0.15
        )
        dense_best_gap = min((c.gap_planar for c in dense_candidates), default=float("inf"))
        reref = hs.refine_candidate(
            torus,
            torus,
            hs.ManifoldCandidate(
                theta2_u=mesh_target["theta2_u"],
                t_u=mesh_target["t_u"],
                theta2_s=mesh_target["theta2_s"],
                t_s=mesh_target["t_s"],
                gap_planar=0.0,
            ),
            lobe_sign_u=mesh_target["u_lobe"],
            lobe_sign_s=mesh_target["s_lobe"],
            n_segments_dir=48,
        )
        reref_corrected_pos_gap_km = (
            float(np.linalg.norm(reref.state_u[:3] - reref.state_s[:3])) * _L_KM_EUROPA
            if reref is not None
            else None
        )
        mesh_check = {
            "target_kind": mesh_target_kind,
            "target_was_genuine": best is not None,
            "target_was_robust_genuine": best_robust is not None,
            "n_theta2_dense": n_theta2_dense,
            "n_time_dense": n_time_dense,
            "dense_grid_coarse_best_gap_planar": dense_best_gap,
            "rerefine_at_n_segments_dir_48_pos_gap_km_corrected": reref_corrected_pos_gap_km,
            "reref_matches_original": (
                abs(reref_corrected_pos_gap_km - mesh_target["corrected_pos_gap_km"]) < 5.0
                if reref_corrected_pos_gap_km is not None
                else False
            ),
        }
        _log(
            f"mesh-refinement check: dense-grid re-refine pos_gap(corrected)="
            f"{reref_corrected_pos_gap_km}",
            t0,
        )

    result = {
        "task": "#703",
        "system": (
            "Jupiter-Europa-Callisto CCR4BP, Europa 4:1 (interior) resonant torus (physical mass)"
        ),
        "physical_unit_note": (
            "search.ccr4bp_heteroclinic_search hardcodes Europa's SMA/GM (_L_KM=671100 km) for "
            "all km conversions. For THIS system the base moon IS Europa, so this is NOT a "
            "generality gap here -- module_native_* and corrected_* fields are verified (not "
            "just expected) to agree; the 'unit_coincidence_verified' assertion above confirms "
            "this at runtime before any downstream reporting."
        ),
        "unit_coincidence_verified": unit_coincidence_verified,
        "l_km_europa": _L_KM_EUROPA,
        "v_unit_europa_km_s": _V_UNIT_EUROPA_KM_S,
        "torus_residual_rms": torus.residual_rms,
        "torus_closure_residual": torus.closure_residual,
        "torus_rho_strob": torus.rho_strob,
        "torus_period_tu": torus.period,
        "globalization_params": {
            "n_theta2": n_theta2,
            "n_time": n_time,
            "t_max_periods": t_max_periods,
            "n_segments_dir": n_segments_dir,
            "eps": mg.DEFAULT_EPS,
        },
        "combo_results": combo_results,
        "best_genuine_connection_corrected": best,
        "best_robust_genuine_connection_corrected": best_robust,
        "ghost_sensitivity_check": ghost_sensitivity_check,
        "robust_integrator_margin_km": _ROBUST_INTEGRATOR_MARGIN_KM,
        "closest_candidate_overall": mesh_target,
        "mesh_refinement_check": mesh_check,
        "base_orbit_note": (
            "Base resonant orbit is spacecraft:Europa=4:1 (interior, a~0.397), following #696's "
            "own Io-Ganymede precedent for a near-4:1-labeled system -- the naive exterior "
            "'moon-ratio-matching' 1:5 reading does not hard-fail under Callisto's physical "
            "forcing but has only a thin ~1.1x-Hill-scale collision margin (vs the chosen "
            "reading's ~22x). See tests/search/test_ccr4bp_torus_europa_callisto.py for the "
            "full collision-risk re-derivation."
        ),
        "scientific_motivation_caveat": (
            "Europa:Callisto's period ratio (~4.73) has NO clean low-integer commensurability, "
            "unlike every other CCR4BP pair this project has built (JEG 2.03, Io-Europa 2.00, "
            "Io-Ganymede 4.06, Umbriel-Titania 2.10). The chosen base orbit is the "
            "least-degenerate tractable seed available, not a resonance-anchored one -- a real "
            "difference in scientific motivation from every prior build, carried forward "
            "honestly from #700's own literature-check note, not smoothed over."
        ),
    }

    out_path = OUT_DIR / "result.json"
    out_path.write_text(json.dumps(result, indent=2, default=float))
    _log(f"wrote {out_path}", t0)
    return result


if __name__ == "__main__":
    main()
