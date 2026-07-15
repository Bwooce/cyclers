"""#600 -- direct symmetric-closure enumeration, generalized to 3-MOON chains.

Extends ``scripts/enumerate_563_symmetric_closures.py``'s direct-construction
method (anchor -> flyby -> anchor, at commensurate ``tof = n*T_syn/2``, no
iterative search) to the untested case of a 3-moon chain:
``anchor -> flyby1 -> flyby2 -> anchor``. Every prior run of the #563 method
(Uranus non-Miranda moons -> the #312 hit; Saturn Titan-Iapetus -> clean
negative; Jupiter Galilean moons -> 36 gate-passing / 0 novelty-clear) has
only ever enumerated 2-moon (single anchor-flyby direction) sequences. A
longer 3-moon tour has never been tried with this proven construction+gate
method -- only the much weaker original discovery genome (grid search, with
the attendant grid-resolution risk #563 was built to replace) ever touched
3-moon Uranian chains.

Generalizing the construction
------------------------------
The 2-moon method builds 3 states -- anchor(t=0), flyby(t=tof), anchor(t=2*tof)
-- from a SINGLE commensurate leg time ``tof = n*T_syn(anchor,flyby)/2`` and a
SINGLE relative-offset parameter (0 or 180 deg), because both legs traverse
the SAME pair and a symmetric (mirror-image) periodic-orbit family lives
exactly at those points (#562's own empirical grid sweep, reproduced 8.9e-16
km/s by #563 at the exact catalogued #312 point).

A 3-moon chain has THREE distinct pairs (anchor,flyby1), (flyby1,flyby2),
(flyby2,anchor) -- there is no single shared synodic period, so the 2-moon
mirror-symmetry argument does not carry over verbatim. The natural, direct
(non-iterative) generalization used here: construct FOUR states from THREE
independently-commensurate leg times, one per pair --

    anchor(t=0),  phase 0 (global phase is provably redundant under the
                  rotational symmetry #558's own docstring establishes --
                  the whole 3-body assembly is a rigid rotation, so WLOG
                  the anchor's initial phase is fixed at 0)
    flyby1(t=tof1),  phase rel_offset_1 at t=0, tof1 = n1*T_syn(anchor,flyby1)/2
    flyby2(t=tof1+tof2),  phase rel_offset_2 at t=0, tof2 = n2*T_syn(flyby1,flyby2)/2
    anchor(t=tof1+tof2+tof3),  tof3 = n3*T_syn(flyby2,anchor)/2

-- each ``n_i`` ranging over ``1..n_max_i`` (same per-pair bound as #563:
``n_max = floor(2*tof_max/T_syn)``, ``tof_max = 3.0*sqrt(P_a*P_b)``, the
literal max ``tof_scale`` #558's own production sweep tested), and each
``rel_offset_i`` in ``{0, 180}`` deg (the same discrete symmetric convention
#563 restricts to, rather than the 0/180-degree pair being independently
re-derived here). This is a bounded, FINITE, direct construction -- nothing is
searched for, exactly as with #563 -- but is an extrapolation, not a proven
theorem: there is no published guarantee that the individual per-pair
commensurability implies exact 3-body closure (unlike the 2-moon case, where
the commensurate point IS the analytic symmetric-orbit family). The residual
gate is the real arbiter: if the extrapolation doesn't hold, candidates will
simply show large residuals and fail the gate -- a valid, informative clean
negative, not a broken method.

Chaining: 3 legs share TWO interior nodes (flyby1, flyby2) plus the anchor
periodicity condition (leg C's arrival vs leg A's departure) -- three
continuity residuals in total (vs. the 2-moon case's two: one interior +
one periodic). The worst (max) of the three is the reported residual, same
convention as #558/#563.

Gate reuse (verbatim, no new gate logic)
------------------------------------------
``GATE_RESIDUAL_KMS``, ``candidate_passes_physical_gate``, and
``dop853_cross_check_leg`` are #558's own gate machinery, imported and used
unmodified. Only the LEG-BUILDING plumbing (3 legs instead of 2, 4-node
instead of 3-node ``vinf`` bookkeeping) is new here, mirroring #558's own
``build_legs_for_record``/``gate_candidate``/``encounter_vinfs_kms`` pattern.

Sanity check (run automatically, --skip-selfcheck to disable)
-----------------------------------------------------------------
Two cheap self-consistency checks, NOT a re-derivation of the gate:
(1) rotation invariance -- the whole 4-state construction rotated by a
    constant angle (nonzero global phase) must give byte-for-byte the same
    residual, the same rotational-symmetry argument #558 already established
    for the 2-moon case, now checked for 3 legs.
(2) cross-check the batched (leg-options-computed-once, reused across all
    n_rev combos) construction against a slow, independent per-combo
    reconstruction (fresh Lambert solves every time, no caching) at a
    handful of random points -- confirms the batching optimization here
    does not silently diverge from a naive direct implementation.

Discipline: NO catalogue writeback, NO V1-V4-strict gauntlet run here.

Run as::

    uv run python scripts/enumerate_600_3moon_symmetric_closures.py
"""

from __future__ import annotations

import itertools
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

# Reuse #558's own gate machinery + leg-construction helpers verbatim.
# Reuse #563's own per-pair n_max helper verbatim.
from enumerate_563_symmetric_closures import (  # noqa: E402
    NON_MIRANDA_MOONS,
    TOF_SCALE_MAX,
    pair_n_max,
)
from scan_558_uranus_all_pairs_offset_sweep import (  # noqa: E402
    GATE_RESIDUAL_KMS,
    N_REV_MAX,
    _leg_options,
)

from cyclerfinder.core.lambert import lambert  # noqa: E402
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES  # noqa: E402
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

DATA_DIR = ROOT / "data"
OUT_PATH = DATA_DIR / "enumerate_600_3moon_symmetric_closures.jsonl"

REL_OFFSETS_DEG: tuple[float, ...] = (0.0, 180.0)
N_REV_VALUES: tuple[int, ...] = tuple(range(N_REV_MAX + 1))  # 0..3


# --------------------------------------------------------------------------
# Construction
# --------------------------------------------------------------------------


def _states_3moon(
    anchor: str,
    flyby1: str,
    flyby2: str,
    *,
    rel_offset1_deg: float,
    rel_offset2_deg: float,
    tof1_days: float,
    tof2_days: float,
    tof3_days: float,
    primary: str,
    phase0_deg: float = 0.0,
) -> tuple[
    np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray
]:
    """The 4 states (r0,v0)..(r3,v3) for one (tof-triple, offset-pair) point.

    ``phase0_deg`` is a global rotation of the whole assembly, default 0 --
    exposed only for the rotation-invariance self-check below (see module
    docstring); every real enumeration call leaves it at 0 (WLOG, per the
    same rotational-symmetry argument #558 established for the 2-moon case).
    """
    mu = PRIMARIES[primary]
    sat_a = SATELLITES[anchor]
    sat_1 = SATELLITES[flyby1]
    sat_2 = SATELLITES[flyby2]
    n_a = _mean_motion_rad_day(mu, sat_a.sma_km)
    n_1 = _mean_motion_rad_day(mu, sat_1.sma_km)
    n_2 = _mean_motion_rad_day(mu, sat_2.sma_km)

    phase0 = math.radians(phase0_deg)
    ro1 = math.radians(rel_offset1_deg)
    ro2 = math.radians(rel_offset2_deg)

    r0, v0 = _moon_state(phase0, n_a, 0.0, sat_a.sma_km, mu)
    r1, v1 = _moon_state(phase0 + ro1, n_1, tof1_days, sat_1.sma_km, mu)
    r2, v2 = _moon_state(phase0 + ro2, n_2, tof1_days + tof2_days, sat_2.sma_km, mu)
    r3, v3 = _moon_state(phase0, n_a, tof1_days + tof2_days + tof3_days, sat_a.sma_km, mu)
    return r0, v0, r1, v1, r2, v2, r3, v3


def three_leg_options(
    anchor: str,
    flyby1: str,
    flyby2: str,
    *,
    rel_offset1_deg: float,
    rel_offset2_deg: float,
    tof1_days: float,
    tof2_days: float,
    tof3_days: float,
    primary: str = "Uranus",
    n_rev_max: int = N_REV_MAX,
    phase0_deg: float = 0.0,
) -> tuple[dict[int, Any], dict[int, Any], dict[int, Any]]:
    """Solve each of the 3 legs ONCE (all achievable n_rev in one Lambert call
    each), for reuse across the full n_rev-combination loop -- mirrors #558's
    ``sweep_pair`` efficiency pattern (one solve per grid point, not per
    n_rev combo)."""
    mu = PRIMARIES[primary]
    r0, v0, r1, v1, r2, v2, r3, v3 = _states_3moon(
        anchor,
        flyby1,
        flyby2,
        rel_offset1_deg=rel_offset1_deg,
        rel_offset2_deg=rel_offset2_deg,
        tof1_days=tof1_days,
        tof2_days=tof2_days,
        tof3_days=tof3_days,
        primary=primary,
        phase0_deg=phase0_deg,
    )
    opts_a = _leg_options(r0, v0, r1, v1, tof1_days * DAY_S, mu, n_rev_max)
    opts_b = _leg_options(r1, v1, r2, v2, tof2_days * DAY_S, mu, n_rev_max)
    opts_c = _leg_options(r2, v2, r3, v3, tof3_days * DAY_S, mu, n_rev_max)
    return opts_a, opts_b, opts_c


def residual_from_options(
    opts_a: dict[int, Any],
    opts_b: dict[int, Any],
    opts_c: dict[int, Any],
    n_rev: tuple[int, int, int],
) -> dict[str, Any] | None:
    """Residual for one explicit n_rev triple, from already-solved leg options.

    Three continuity conditions (mirrors #558's two-condition ``r_mid`` /
    ``r_periodic``, extended to the two interior nodes of a 3-leg chain plus
    the anchor periodicity condition):
      - flyby1 continuity: leg A arrival V_inf vs leg B departure V_inf
      - flyby2 continuity: leg B arrival V_inf vs leg C departure V_inf
      - anchor periodicity: leg A departure V_inf vs leg C arrival V_inf
    ``worst`` = max of the three, same "worst-case residual" convention.
    """
    na, nb, nc = n_rev
    if na not in opts_a or nb not in opts_b or nc not in opts_c:
        return None
    oa, ob, oc = opts_a[na], opts_b[nb], opts_c[nc]
    r_flyby1 = abs(oa.vinf_in - ob.vinf_out)
    r_flyby2 = abs(ob.vinf_in - oc.vinf_out)
    r_periodic = abs(oa.vinf_out - oc.vinf_in)
    worst = max(r_flyby1, r_flyby2, r_periodic)
    return {
        "residual_kms": worst,
        "r_flyby1_kms": r_flyby1,
        "r_flyby2_kms": r_flyby2,
        "r_periodic_kms": r_periodic,
        "vinf_in": [0.0, oa.vinf_in, ob.vinf_in, oc.vinf_in],
        "vinf_out": [oa.vinf_out, ob.vinf_out, oc.vinf_out, 0.0],
    }


def build_legs_for_record_3moon(
    anchor: str, flyby1: str, flyby2: str, rec: dict[str, Any], *, primary: str = "Uranus"
) -> list[PatchedConicLeg]:
    """Reconstruct SI-units PatchedConicLeg objects (3 legs) for one gated
    record -- mirrors #558's ``build_legs_for_record`` verbatim pattern,
    extended from 2 legs to 3."""
    mu = PRIMARIES[primary]
    tof1, tof2, tof3 = rec["tof1_days"], rec["tof2_days"], rec["tof3_days"]
    na, nb, nc = rec["n_rev"]

    r0, v0, r1, v1, r2, v2, r3, _v3 = _states_3moon(
        anchor,
        flyby1,
        flyby2,
        rel_offset1_deg=rec["rel_offset1_deg"],
        rel_offset2_deg=rec["rel_offset2_deg"],
        tof1_days=tof1,
        tof2_days=tof2,
        tof3_days=tof3,
        primary=primary,
    )
    tof1_s, tof2_s, tof3_s = tof1 * DAY_S, tof2 * DAY_S, tof3 * DAY_S

    def _best(sols: list[Any], n_target: int, v_ref: np.ndarray) -> Any:
        return min(
            (s for s in sols if s.n_revs == n_target),
            key=lambda s: float(np.linalg.norm(s.v1 - v_ref)),
        )

    sols_a = lambert(r0, r1, tof1_s, mu=mu, max_revs=max(0, na))
    best_a = _best(sols_a, na, v0)
    sols_b = lambert(r1, r2, tof2_s, mu=mu, max_revs=max(0, nb))
    best_b = _best(sols_b, nb, v1)
    sols_c = lambert(r2, r3, tof3_s, mu=mu, max_revs=max(0, nc))
    best_c = _best(sols_c, nc, v2)

    km_m = 1000.0
    mu_m3_s2 = mu * (km_m**3)
    return [
        PatchedConicLeg(
            label_from=anchor,
            label_to=flyby1,
            r1_m=r0 * km_m,
            v1_m_s=best_a.v1 * km_m,
            r2_m=r1 * km_m,
            v2_m_s=best_a.v2 * km_m,
            dt_s=tof1_s,
            mu_m3_s2=mu_m3_s2,
        ),
        PatchedConicLeg(
            label_from=flyby1,
            label_to=flyby2,
            r1_m=r1 * km_m,
            v1_m_s=best_b.v1 * km_m,
            r2_m=r2 * km_m,
            v2_m_s=best_b.v2 * km_m,
            dt_s=tof2_s,
            mu_m3_s2=mu_m3_s2,
        ),
        PatchedConicLeg(
            label_from=flyby2,
            label_to=anchor,
            r1_m=r2 * km_m,
            v1_m_s=best_c.v1 * km_m,
            r2_m=r3 * km_m,
            v2_m_s=best_c.v2 * km_m,
            dt_s=tof3_s,
            mu_m3_s2=mu_m3_s2,
        ),
    ]


def encounter_vinfs_kms_3moon(rec: dict[str, Any]) -> tuple[float, float, float, float]:
    """Per-encounter V_inf magnitude for the 4-node sequence, matching
    #558's ``encounter_vinfs_kms`` convention (max of in/out asymptote at
    each node), extended from 3 nodes to 4."""
    vin = rec["vinf_in"]
    vout = rec["vinf_out"]
    return tuple(max(abs(vin[k]), abs(vout[k])) for k in range(4))  # type: ignore[return-value]


def gate_candidate_3moon(
    anchor: str, flyby1: str, flyby2: str, rec: dict[str, Any], *, primary: str = "Uranus"
) -> dict[str, Any]:
    """#558's ``gate_candidate``, extended verbatim in spirit to a 4-node /
    3-leg sequence: physical bend gate (body-agnostic, arbitrary sequence
    length) + independent DOP853 cross-check on all 3 legs."""
    seq = (anchor, flyby1, flyby2, anchor)
    vinfs = encounter_vinfs_kms_3moon(rec)
    physical_pass, verdicts = candidate_passes_physical_gate(
        seq, vinfs, min_useful_bend_deg=DEFAULT_MIN_USEFUL_BEND_DEG
    )
    legs = build_legs_for_record_3moon(anchor, flyby1, flyby2, rec, primary=primary)
    cross_checks = [dop853_cross_check_leg(leg, rtol=1e-12, atol=1e-12) for leg in legs]
    max_dr_km = max(float(cc["dr_arrival_km"]) for cc in cross_checks)
    independent_pass = max_dr_km < 1.0  # matches #558/verify_327's threshold

    return {
        "anchor": anchor,
        "flyby1": flyby1,
        "flyby2": flyby2,
        "record": rec,
        "vinf_per_encounter_kms": list(vinfs),
        "physical_gate_passed": physical_pass,
        "max_bend_deg_per_encounter": [v.max_bend_deg for v in verdicts],
        "dop853_cross_check": {
            "max_dr_arrival_km": max_dr_km,
            "per_leg": cross_checks,
            "passed": independent_pass,
        },
        "all_gates_passed": bool(physical_pass and independent_pass),
    }


# --------------------------------------------------------------------------
# Enumeration
# --------------------------------------------------------------------------


def enumerate_sequence(
    anchor: str,
    flyby1: str,
    flyby2: str,
    *,
    primary: str = "Uranus",
    tof_scale_max: float = TOF_SCALE_MAX,
) -> dict[str, Any]:
    """Construct + gate every symmetric-extrapolation candidate for one
    ``anchor -> flyby1 -> flyby2 -> anchor`` sequence."""
    t_syn_a, _, _, n_max_a = pair_n_max(
        anchor, flyby1, primary=primary, tof_scale_max=tof_scale_max
    )
    t_syn_b, _, _, n_max_b = pair_n_max(
        flyby1, flyby2, primary=primary, tof_scale_max=tof_scale_max
    )
    t_syn_c, _, _, n_max_c = pair_n_max(
        flyby2, anchor, primary=primary, tof_scale_max=tof_scale_max
    )

    n_evaluated = 0
    n_infeasible = 0
    n_subgate = 0
    passes: list[dict[str, Any]] = []

    for n1 in range(1, n_max_a + 1):
        tof1 = n1 * t_syn_a / 2.0
        for n2 in range(1, n_max_b + 1):
            tof2 = n2 * t_syn_b / 2.0
            for n3 in range(1, n_max_c + 1):
                tof3 = n3 * t_syn_c / 2.0
                for ro1, ro2 in itertools.product(REL_OFFSETS_DEG, REL_OFFSETS_DEG):
                    opts_a, opts_b, opts_c = three_leg_options(
                        anchor,
                        flyby1,
                        flyby2,
                        rel_offset1_deg=ro1,
                        rel_offset2_deg=ro2,
                        tof1_days=tof1,
                        tof2_days=tof2,
                        tof3_days=tof3,
                        primary=primary,
                    )
                    for na, nb, nc in itertools.product(N_REV_VALUES, N_REV_VALUES, N_REV_VALUES):
                        n_evaluated += 1
                        pt = residual_from_options(opts_a, opts_b, opts_c, (na, nb, nc))
                        if pt is None:
                            n_infeasible += 1
                            continue
                        if pt["residual_kms"] >= GATE_RESIDUAL_KMS:
                            continue
                        n_subgate += 1
                        rec = {
                            "anchor": anchor,
                            "flyby1": flyby1,
                            "flyby2": flyby2,
                            "n1_commensurate_int": n1,
                            "n2_commensurate_int": n2,
                            "n3_commensurate_int": n3,
                            "t_syn_a_days": t_syn_a,
                            "t_syn_b_days": t_syn_b,
                            "t_syn_c_days": t_syn_c,
                            "rel_offset1_deg": ro1,
                            "rel_offset2_deg": ro2,
                            "n_rev": [na, nb, nc],
                            "tof1_days": tof1,
                            "tof2_days": tof2,
                            "tof3_days": tof3,
                            "tof_total_days": tof1 + tof2 + tof3,
                            **pt,
                        }
                        gated = gate_candidate_3moon(anchor, flyby1, flyby2, rec, primary=primary)
                        if gated["all_gates_passed"]:
                            passes.append(
                                {
                                    **rec,
                                    "vinf_per_encounter_kms": gated["vinf_per_encounter_kms"],
                                    "max_bend_deg_per_encounter": gated[
                                        "max_bend_deg_per_encounter"
                                    ],
                                    "dop853_cross_check": gated["dop853_cross_check"],
                                }
                            )

    return {
        "anchor": anchor,
        "flyby1": flyby1,
        "flyby2": flyby2,
        "sequence": [anchor, flyby1, flyby2, anchor],
        "t_syn_a_days": t_syn_a,
        "t_syn_b_days": t_syn_b,
        "t_syn_c_days": t_syn_c,
        "n_max_a": n_max_a,
        "n_max_b": n_max_b,
        "n_max_c": n_max_c,
        "n_evaluated": n_evaluated,
        "n_infeasible": n_infeasible,
        "n_subgate_residual_only": n_subgate,
        "n_all_gates_passed": len(passes),
        "passes": passes,
    }


# --------------------------------------------------------------------------
# Self-checks (cheap; NOT a re-derivation of the gate -- see module docstring)
# --------------------------------------------------------------------------


def _selfcheck_rotation_invariance(primary: str = "Uranus") -> bool:
    """Rotating the whole assembly (nonzero global phase) must not change the
    residual -- the same rotational-symmetry argument #558 established for
    the 2-moon case, checked here for the 3-leg construction."""
    anchor, flyby1, flyby2 = "Ariel", "Titania", "Oberon"
    t_syn_a, _, _, _ = pair_n_max(anchor, flyby1, primary=primary)
    t_syn_b, _, _, _ = pair_n_max(flyby1, flyby2, primary=primary)
    t_syn_c, _, _, _ = pair_n_max(flyby2, anchor, primary=primary)
    tof1, tof2, tof3 = 2 * t_syn_a / 2.0, 3 * t_syn_b / 2.0, 1 * t_syn_c / 2.0
    # Find any n_rev triple that's feasible at phase0=0 first (n_rev=(0,0,0),
    # the direct/no-revolution transfer, is almost always achievable), then
    # check that SAME triple's residual is phase0-invariant.
    opts_a0, opts_b0, opts_c0 = three_leg_options(
        anchor,
        flyby1,
        flyby2,
        rel_offset1_deg=0.0,
        rel_offset2_deg=180.0,
        tof1_days=tof1,
        tof2_days=tof2,
        tof3_days=tof3,
        primary=primary,
        phase0_deg=0.0,
    )
    n_rev_triple = None
    for na, nb, nc in itertools.product(N_REV_VALUES, N_REV_VALUES, N_REV_VALUES):
        if residual_from_options(opts_a0, opts_b0, opts_c0, (na, nb, nc)) is not None:
            n_rev_triple = (na, nb, nc)
            break
    if n_rev_triple is None:
        print(
            "[600 selfcheck] rotation invariance: NO feasible n_rev triple found -- fail",
            flush=True,
        )
        return False

    residuals = []
    for phase0_deg in (0.0, 37.0, 271.0):
        opts_a, opts_b, opts_c = three_leg_options(
            anchor,
            flyby1,
            flyby2,
            rel_offset1_deg=0.0,
            rel_offset2_deg=180.0,
            tof1_days=tof1,
            tof2_days=tof2,
            tof3_days=tof3,
            primary=primary,
            phase0_deg=phase0_deg,
        )
        pt = residual_from_options(opts_a, opts_b, opts_c, n_rev_triple)
        residuals.append(pt["residual_kms"] if pt is not None else math.nan)
    ok = all(math.isfinite(r) for r in residuals) and max(residuals) - min(residuals) < 1e-9
    print(f"[600 selfcheck] rotation invariance: residuals={residuals} ok={ok}", flush=True)
    return ok


def _selfcheck_naive_agreement(primary: str = "Uranus") -> bool:
    """Cross-check the batched construction against an independent, naive
    per-combo reconstruction (fresh Lambert solves, no shared caching) at a
    handful of points -- confirms the batching optimization matches a
    straightforward direct implementation of the same construction."""
    import random

    rng = random.Random(600)
    sequences = list(itertools.permutations(NON_MIRANDA_MOONS, 3))
    ok = True
    for _ in range(6):
        anchor, flyby1, flyby2 = rng.choice(sequences)
        t_syn_a, _, _, n_max_a = pair_n_max(anchor, flyby1, primary=primary)
        t_syn_b, _, _, n_max_b = pair_n_max(flyby1, flyby2, primary=primary)
        t_syn_c, _, _, n_max_c = pair_n_max(flyby2, anchor, primary=primary)
        n1 = rng.randint(1, n_max_a)
        n2 = rng.randint(1, n_max_b)
        n3 = rng.randint(1, n_max_c)
        tof1, tof2, tof3 = n1 * t_syn_a / 2.0, n2 * t_syn_b / 2.0, n3 * t_syn_c / 2.0
        ro1 = rng.choice(REL_OFFSETS_DEG)
        ro2 = rng.choice(REL_OFFSETS_DEG)
        na, nb, nc = (rng.randint(0, N_REV_MAX) for _ in range(3))

        # Batched path.
        opts_a, opts_b, opts_c = three_leg_options(
            anchor,
            flyby1,
            flyby2,
            rel_offset1_deg=ro1,
            rel_offset2_deg=ro2,
            tof1_days=tof1,
            tof2_days=tof2,
            tof3_days=tof3,
            primary=primary,
        )
        batched = residual_from_options(opts_a, opts_b, opts_c, (na, nb, nc))

        # Naive path: fresh, independent, per-leg lambert() calls with
        # max_revs=max(n,1) each (mirrors #563's residual_at_point style).
        mu = PRIMARIES[primary]
        r0, v0, r1, v1, r2, v2, r3, v3 = _states_3moon(
            anchor,
            flyby1,
            flyby2,
            rel_offset1_deg=ro1,
            rel_offset2_deg=ro2,
            tof1_days=tof1,
            tof2_days=tof2,
            tof3_days=tof3,
            primary=primary,
        )
        naive_a = _leg_options(r0, v0, r1, v1, tof1 * DAY_S, mu, max(na, 1))
        naive_b = _leg_options(r1, v1, r2, v2, tof2 * DAY_S, mu, max(nb, 1))
        naive_c = _leg_options(r2, v2, r3, v3, tof3 * DAY_S, mu, max(nc, 1))
        naive = residual_from_options(naive_a, naive_b, naive_c, (na, nb, nc))

        if batched is None and naive is None:
            agree = True
        elif batched is None or naive is None:
            agree = False
        else:
            agree = abs(batched["residual_kms"] - naive["residual_kms"]) < 1e-9
        ok = ok and agree
        print(
            f"[600 selfcheck] {anchor}-{flyby1}-{flyby2}-{anchor} n=({n1},{n2},{n3}) "
            f"nrev=({na},{nb},{nc}) ro=({ro1},{ro2}): batched={batched} naive={naive} "
            f"agree={agree}",
            flush=True,
        )
    print(f"[600 selfcheck] naive-agreement overall ok={ok}", flush=True)
    return ok


def main(argv: list[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=str, default="Uranus")
    parser.add_argument(
        "--moons",
        type=str,
        default=",".join(NON_MIRANDA_MOONS),
        help="Comma-separated moon list; every ordered permutation of 3 among them is enumerated.",
    )
    parser.add_argument("--tof-scale-max", type=float, default=TOF_SCALE_MAX)
    parser.add_argument("--out", type=str, default="")
    parser.add_argument(
        "--skip-selfcheck", action="store_true", help="Skip the two cheap sanity checks."
    )
    args = parser.parse_args(argv)

    primary = args.primary
    moons = tuple(m.strip() for m in args.moons.split(",") if m.strip())
    tof_scale_max = args.tof_scale_max
    is_default_uranian = (
        primary == "Uranus" and moons == NON_MIRANDA_MOONS and tof_scale_max == TOF_SCALE_MAX
    )
    if args.out:
        out_path = Path(args.out)
    elif is_default_uranian:
        out_path = OUT_PATH
    else:
        raise SystemExit("--out is required when --primary/--moons/--tof-scale-max are non-default")

    if not args.skip_selfcheck:
        rot_ok = _selfcheck_rotation_invariance(primary=primary)
        naive_ok = _selfcheck_naive_agreement(primary=primary)
        if not (rot_ok and naive_ok):
            print("[600] ABORT: self-check FAILED -- not proceeding to full sweep.", flush=True)
            return 1

    t0 = time.time()
    sequences = list(itertools.permutations(moons, 3))
    print(
        f"[600] primary={primary} moons={moons} {len(sequences)} 3-moon sequences x "
        f"{len(N_REV_VALUES) ** 3} n_rev combos x {len(REL_OFFSETS_DEG) ** 2} offset combos, "
        f"tof_scale_max={tof_scale_max}",
        flush=True,
    )

    all_results: list[dict[str, Any]] = []
    total_evaluated = 0
    total_passes = 0
    for anchor, flyby1, flyby2 in sequences:
        res = enumerate_sequence(
            anchor, flyby1, flyby2, primary=primary, tof_scale_max=tof_scale_max
        )
        all_results.append(res)
        total_evaluated += res["n_evaluated"]
        total_passes += res["n_all_gates_passed"]
        print(
            f"[600] {anchor}-{flyby1}-{flyby2}-{anchor}: "
            f"n_max=({res['n_max_a']},{res['n_max_b']},{res['n_max_c']}) "
            f"evaluated={res['n_evaluated']} sub_gate={res['n_subgate_residual_only']} "
            f"all_gates_pass={res['n_all_gates_passed']}",
            flush=True,
        )

    elapsed = time.time() - t0
    print(
        f"[600] DONE: {total_evaluated} candidates directly evaluated across "
        f"{len(sequences)} sequences, {total_passes} pass ALL gates (residual+bend+DOP853) "
        f"({elapsed:.1f}s)",
        flush=True,
    )

    with out_path.open("w", encoding="utf-8") as fh:
        fh.write(
            json.dumps(
                {
                    "_meta": True,
                    "task": "#600 3-moon chain direct symmetric-closure enumeration",
                    "primary": primary,
                    "moons": list(moons),
                    "sequences": len(sequences),
                    "n_rev_combos": len(N_REV_VALUES) ** 3,
                    "rel_offset_combos": len(REL_OFFSETS_DEG) ** 2,
                    "tof_scale_max_bound": tof_scale_max,
                    "total_evaluated": total_evaluated,
                    "total_all_gates_passed": total_passes,
                    "elapsed_s": elapsed,
                    "gate_residual_kms": GATE_RESIDUAL_KMS,
                }
            )
            + "\n"
        )
        for res in all_results:
            fh.write(
                json.dumps(
                    {
                        "kind": "sequence_summary",
                        "anchor": res["anchor"],
                        "flyby1": res["flyby1"],
                        "flyby2": res["flyby2"],
                        "sequence": res["sequence"],
                        "n_max_a": res["n_max_a"],
                        "n_max_b": res["n_max_b"],
                        "n_max_c": res["n_max_c"],
                        "n_evaluated": res["n_evaluated"],
                        "n_infeasible": res["n_infeasible"],
                        "n_subgate_residual_only": res["n_subgate_residual_only"],
                        "n_all_gates_passed": res["n_all_gates_passed"],
                    }
                )
                + "\n"
            )
            for p in res["passes"]:
                fh.write(json.dumps({"kind": "pass", **p}) + "\n")
    print(f"[600] written: {out_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
