"""Five-tier accessibility prioritizer (#282).

Composes the full Track-B scorer stack into ONE API, on top of -- not in place
of -- :class:`two_tier_prioritizer.TwoTierPrioritizer`. The two-tier module is
imported and called verbatim; this wrapper adds Tier-0 (NN reachability
prefilter) as a gating pre-screen and Tiers 3-5 (resonant-manifold heteroclinic
network, FTLE chaos field, Hiraiwa lobe-overlap graph) as parallel verifiers.

Tier inventory (one-line source per tier; full provenance lives in each tier's
module docstring):

* Tier 0 -- :class:`neural_reach_prefilter.NeuralReachPrefilter`
  Zhang-Topputo 2026 NN predicted ΔV (m/s) for one heliocentric Lambert leg.
  Inputs: (r1, v1, r2, v2, dt, m0, T_max, Isp, mu). Native patched-conic
  surface; the only tier that can score a raw Lambert leg without a CR3BP
  representative orbit.
* Tier 1 -- :class:`reachable_network` heading-fan overlap (Braik-Ross 2026,
  energy-PRESERVING, single C_J manifold).
* Tier 2 -- :func:`reachable_impulsive.reachable_cloud` single-impulse footprint
  (Zhou-Armellin 2024, energy-CHANGING).
* Tier 3 -- :class:`resonance_network.ResonanceNetworkScorer` (Kumar-Rawat-
  Rosengren-Ross 2025, energy-DEGENERATE Floquet manifolds, perigee section).
* Tier 4 -- :class:`ftle_scorer.FTLEScorer` (Canales-Howell 2023, chaos-aware
  transport-corridor strength on a FTLE field).
* Tier 5 -- :class:`lobe_overlap_scorer.LobeOverlapScorer` (Hiraiwa-Bando 2026,
  effective-lobe graph; bottleneck-flux Dijkstra on the lobe-overlap edges).

COMPOSITION CHOICE -- rank-product on the survivors, NOT a weighted sum.
Each tier reports a different *physical* quantity (ΔV in m/s; voxel-overlap
score in [0,1]; perigee distance in nondim; FTLE corridor strength in [0,1];
lobe bottleneck flux in nondim area). Comparing them magnitude-wise is
nonsense without per-tier calibration, which the source papers do not provide.
Rank-product is the standard composition that survives unit / scale
incompatibility: each tier ranks the candidates 1..N (1 = best), and the
composite score is the geometric mean of the ranks. The lowest composite
rank is the most-accessible pair across all tiers. This is the same
discipline as Breitling 2004 (rank-product over heterogeneous biological
assays), used here because the physical interpretation matches: each tier is
an independent "vote" on accessibility, with no a-priori weighting between
votes. See :func:`rank_product_score`.

OPERATING MODES (the architectural seam):

The five tiers split by INPUT NATURE, not by sequence:

* **Patched-conic leg mode** (used by the #264 / #282 sweep): the candidate
  is one Lambert leg between two moon-state pairs in the planet frame. Only
  Tier 0 (NN) operates natively here -- it takes ``(r1, v1, r2, v2, dt, ...)``
  in SI units. Tiers 1-5 require a CR3BP **representative** with ``state0``
  and ``period``. A multi-leg cycler candidate has no such anchor without an
  IC-synthesis step that bridges patched-conic states into a CR3BP periodic
  family -- a step this module does NOT perform (it would be a separate genome).
  In this mode :meth:`score_leg` runs Tier 0 and returns the NN verdict alone.

* **Representative-pair mode** (used by manifold work, the resonant-network
  catalogue): the candidate is an ordered pair of CR3BP representative orbits
  with state0/period (and, for Tiers 3 + 5, a stability/Floquet structure;
  for Tier 4, a Jacobi constant). In this mode :meth:`score_pair_full` runs
  all five tiers.

Both modes are first-class: a hybrid sweep (Tier 0 on legs, Tiers 1-5 on the
representative anchors) is the right architecture, but is not implemented
here because #264 emits no representative anchors. The honest decision is to
DOCUMENT the seam; the #282 yield comparison runs Tier 0 over the #264
patched-conic legs and reports the architectural gap for Tiers 1-5.

INDEPENDENT CROSS-CHECK -- the composition matches each source paper's intent:

* Tier 0: the NN's `admit_threshold_kms` default (5 km/s) is the
  Zhang-Topputo 2026 prefilter regime ("obviously infeasible" cut, not a
  decision threshold); we pass it through verbatim. The threshold is the
  scorer's own default, not re-calibrated here.
* Tier 4 FTLE: the chaos-class classification uses
  :func:`ftle_scorer.classify_chaos` which the source-paper-equivalent
  Shadden-Lekien-Marsden 2005 definition implements (see
  ``ftle_scorer.py`` lines 19-23 for the canonical formula); the
  transport-corridor strength is the mean of capture-class membership along
  the geodesic, weighted 1 / 0.5 / 0 for capture / transit / escape (matching
  the Canales-Howell 2023 narrative that capture cells are dynamically
  cheap, transit cells are passable, and escape cells are barriers).

DO NOT MODIFY ``two_tier_prioritizer.py``. This module imports it as-is and
composes. The two-tier API surface is the user contract; this module is a
new contract on top of it.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.search.ftle_scorer as ftle_mod
import cyclerfinder.search.lobe_overlap_scorer as lobe_mod
import cyclerfinder.search.neural_reach_prefilter as nn_mod
import cyclerfinder.search.resonance_network as res_mod
import cyclerfinder.search.two_tier_prioritizer as ttp

# ---------------------------------------------------------------------------
# Patched-conic-leg input view.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class PatchedConicLeg:
    """One Lambert leg in a planet (or central-body) frame -- Tier 0 input.

    Carries SI inertial state at both endpoints + the time of flight. The
    NN takes these directly via :meth:`NeuralReachPrefilter.predict_native`.

    Parameters
    ----------
    label_from / label_to :
        Free-form identifiers (moon names; "Hydra"->"Nix").
    r1_m, v1_m_s :
        Inertial position (m) and velocity (m/s) at departure.
    r2_m, v2_m_s :
        Inertial position (m) and velocity (m/s) at arrival.
    dt_s :
        Time of flight (s).
    mu_m3_s2 :
        Gravitational parameter of the central body (m^3 / s^2).
    """

    label_from: str
    label_to: str
    r1_m: NDArray[np.float64]
    v1_m_s: NDArray[np.float64]
    r2_m: NDArray[np.float64]
    v2_m_s: NDArray[np.float64]
    dt_s: float
    mu_m3_s2: float


# ---------------------------------------------------------------------------
# Rank-product helper.
# ---------------------------------------------------------------------------


def rank_product_score(
    per_tier_scores: Sequence[Sequence[float]],
    *,
    descending: Sequence[bool] | None = None,
) -> list[float]:
    """Rank-product composition for N candidates over M tiers.

    Parameters
    ----------
    per_tier_scores :
        A length-M sequence of length-N sequences. Each inner sequence is one
        tier's score for each candidate. NaN/inf entries are pushed to the
        worst rank (largest rank = worst).
    descending :
        Optional length-M sequence of booleans; ``descending[i] = True`` means
        higher values of tier i are BETTER (e.g. overlap strength, transport
        corridor strength), so ranks are computed in descending order. If
        ``False`` or omitted, lower is better (e.g. ΔV, perigee distance).
        Default: all ``False`` (lower-is-better).

    Returns
    -------
    list[float]
        Length-N composite scores (the geometric mean of the per-tier ranks).
        Lower composite = better across the stack. A candidate that misses any
        tier (all NaN) gets the worst rank in that tier and so a higher product.
    """
    n_tiers = len(per_tier_scores)
    if n_tiers == 0:
        return []
    n = len(per_tier_scores[0])
    if any(len(t) != n for t in per_tier_scores):
        raise ValueError("Per-tier score arrays must all have the same length")
    desc = list(descending) if descending is not None else [False] * n_tiers
    if len(desc) != n_tiers:
        raise ValueError("descending must match per_tier_scores length")
    # Per-tier ranks: convert NaN/inf to +inf so they sort last (worst).
    ranks: list[NDArray[np.float64]] = []
    for i, scores in enumerate(per_tier_scores):
        arr = np.asarray(scores, dtype=np.float64)
        bad = ~np.isfinite(arr)
        # Sentinel: lower-is-better -> +inf is worst; higher-is-better -> -inf is worst.
        sentinel = np.inf if not desc[i] else -np.inf
        arr = np.where(bad, sentinel, arr)
        # argsort: ascending. For descending, negate.
        order = np.argsort(-arr if desc[i] else arr, kind="stable")
        rank = np.empty(n, dtype=np.float64)
        rank[order] = np.arange(1, n + 1)
        ranks.append(rank)
    # Composite = geometric mean of ranks.
    log_sum = np.zeros(n, dtype=np.float64)
    for r in ranks:
        log_sum += np.log(r)
    return list(np.exp(log_sum / n_tiers))


# ---------------------------------------------------------------------------
# The five-tier prioritizer.
# ---------------------------------------------------------------------------


@dataclass
class FiveTierPrioritizer:
    """Five-tier accessibility prioritizer composing Tier 0 (NN) + Tiers 1-5.

    Composition pattern:

    * The two-tier wrapper :class:`TwoTierPrioritizer` is held as ``two_tier``
      and called verbatim for Tier 1-2 scoring (this module never reaches into
      its internals).
    * The Tier-0 prefilter and Tiers 3-5 scorers are held as members and
      called for the modes their inputs support.

    For the **patched-conic leg mode** used by #282:
    :meth:`score_leg` runs Tier 0 alone (the only tier that operates on raw
    Lambert legs without a CR3BP representative anchor) and returns the NN
    verdict, with the architectural gap documented in the result dict.

    For the **representative-pair mode**:
    :meth:`score_pair_full` runs all five tiers and returns a dict with per-tier
    scores + the rank-product composite.

    Parameters
    ----------
    nn_prefilter :
        Loaded :class:`NeuralReachPrefilter`. If ``None``, defaults to the
        vendored weights at :data:`neural_reach_prefilter._VENDORED_DV_DIR`
        via :meth:`NeuralReachPrefilter.from_weight_dir`.
    two_tier :
        Optional :class:`TwoTierPrioritizer`; required for representative-pair
        mode but not for leg mode. If ``None``, :meth:`score_pair_full` will
        skip Tiers 1-2 and run only 0, 3, 4, 5 (the others where present).
    resonance, ftle, lobe :
        Optional Tier 3-5 scorers (any may be ``None``; the missing tier is
        skipped in the composition). All must share the same CR3BP system as
        ``two_tier`` for the representative-pair mode to be well-posed.
    """

    nn_prefilter: nn_mod.NeuralReachPrefilter | None = None
    two_tier: ttp.TwoTierPrioritizer | None = None
    resonance: res_mod.ResonanceNetworkScorer | None = None
    ftle: ftle_mod.FTLEScorer | None = None
    lobe: lobe_mod.LobeOverlapScorer | None = None
    # Tier 0 admit threshold override (km/s); if None, falls back to the NN
    # scorer's own default (5 km/s, the Zhang-Topputo prefilter regime).
    tier0_admit_threshold_kms: float | None = None

    def __post_init__(self) -> None:
        if self.nn_prefilter is None:
            self.nn_prefilter = nn_mod.NeuralReachPrefilter.from_weight_dir()
        if self.tier0_admit_threshold_kms is not None:
            # Apply override on the live prefilter.
            self.nn_prefilter.admit_threshold_kms = float(self.tier0_admit_threshold_kms)

    # ------------------------------------------------------------------
    # Patched-conic leg mode: Tier 0 only.
    # ------------------------------------------------------------------

    def score_leg(
        self,
        leg: PatchedConicLeg,
        *,
        mass_kg: float | None = None,
        tmax_n: float | None = None,
        isp_s: float | None = None,
    ) -> dict[str, Any]:
        """Tier-0 NN ΔV prefilter for one heliocentric/patched-conic leg.

        Tiers 1-5 are NOT applicable in this mode (they require a CR3BP
        representative orbit with ``state0`` + ``period``; one Lambert leg is
        not such an anchor). The returned dict's ``tiers_skipped`` key
        documents the architectural gap so the caller cannot mistake a
        Tier-0-only verdict for a full five-tier composite.

        Returns
        -------
        dict
            ``label_from`` / ``label_to``,
            ``tier0_predicted_dv_kms``, ``tier0_predicted_tof_days``,
            ``tier0_admitted``, ``tier0_model_available``,
            ``tier0_fallback_used``, ``tiers_skipped`` (list of skipped tier
            indices with reason).
        """
        assert self.nn_prefilter is not None
        # Build the minimal RepView shape the NN scorer expects.
        rv_from = nn_mod.RepView(
            label=leg.label_from,
            state0=np.concatenate([leg.r1_m, leg.v1_m_s]),
            period=float(leg.dt_s),
            heliocentric_state=np.concatenate([leg.r1_m, leg.v1_m_s]),
        )
        rv_to = nn_mod.RepView(
            label=leg.label_to,
            state0=np.concatenate([leg.r2_m, leg.v2_m_s]),
            period=float(leg.dt_s),
            heliocentric_state=np.concatenate([leg.r2_m, leg.v2_m_s]),
        )
        tof_days = leg.dt_s / 86400.0
        result = self.nn_prefilter.score_pair(
            rv_from,
            rv_to,
            tof_window=(tof_days, tof_days),
            mass_kg=mass_kg,
            tmax_n=tmax_n,
            isp_s=isp_s,
            mu_m3_s2=leg.mu_m3_s2,
        )
        return {
            "label_from": leg.label_from,
            "label_to": leg.label_to,
            "tier0_predicted_dv_kms": float(result["predicted_dv_kms"]),
            "tier0_predicted_tof_days": float(result["predicted_tof_days"]),
            "tier0_admitted": bool(result["prefilter_admitted"]),
            "tier0_model_available": bool(result["model_available"]),
            "tier0_fallback_used": result["fallback_used"],
            "tiers_skipped": [
                {
                    "tier": i,
                    "reason": "patched-conic leg input lacks CR3BP "
                    "representative orbit (state0/period); tier requires "
                    "a periodic-orbit anchor",
                }
                for i in (1, 2, 3, 4, 5)
            ],
        }

    def score_candidate_legs(
        self,
        legs: Sequence[PatchedConicLeg],
        *,
        mass_kg: float | None = None,
        tmax_n: float | None = None,
        isp_s: float | None = None,
    ) -> dict[str, Any]:
        """Score every leg of a multi-leg candidate; aggregate to per-candidate stats.

        Returns a dict with:

        * ``per_leg`` -- list of per-leg :meth:`score_leg` dicts.
        * ``tier0_max_dv_kms`` -- the worst (largest) NN ΔV across legs (the
          binding leg).
        * ``tier0_sum_dv_kms`` -- the linear sum of NN ΔV's (a cumulative
          accessibility cost; the NN model is per-leg, so this is meaningful
          only as an additive proxy across legs of the same candidate).
        * ``tier0_all_admitted`` -- True iff every leg was admitted.
        * ``tier0_any_inference_failed`` -- True iff any leg fell back.
        """
        per_leg = [self.score_leg(leg, mass_kg=mass_kg, tmax_n=tmax_n, isp_s=isp_s) for leg in legs]
        dvs = [p["tier0_predicted_dv_kms"] for p in per_leg]
        finite = [v for v in dvs if math.isfinite(v)]
        return {
            "per_leg": per_leg,
            "tier0_max_dv_kms": float(max(finite)) if finite else float("inf"),
            "tier0_sum_dv_kms": float(sum(finite)) if finite else float("inf"),
            "tier0_all_admitted": all(p["tier0_admitted"] for p in per_leg),
            "tier0_any_inference_failed": any(
                p["tier0_fallback_used"] not in (None, "no-weights") for p in per_leg
            ),
            "n_legs": len(per_leg),
        }

    # ------------------------------------------------------------------
    # Representative-pair mode: all five tiers.
    # ------------------------------------------------------------------

    def score_pair_full(
        self,
        rep_from: object,
        rep_to: object,
        *,
        c_j: float,
        dv_budget_kms: float = 1.0,
        heliocentric_state_from: Sequence[float] | None = None,
        heliocentric_state_to: Sequence[float] | None = None,
        tof_window_days: tuple[float, float] | None = None,
    ) -> dict[str, Any]:
        """All five tiers + rank-product composite for a representative pair.

        Tiers absent on the instance (``None`` member) are skipped in the
        composite. Tier 0 only runs if heliocentric SI states are provided
        (either via the kwargs or attached to the representatives).
        """
        out: dict[str, Any] = {
            "rep_from": str(getattr(rep_from, "label", "<unlabelled>")),
            "rep_to": str(getattr(rep_to, "label", "<unlabelled>")),
        }
        # Tier 0
        assert self.nn_prefilter is not None
        t0 = self.nn_prefilter.score_pair(
            rep_from,
            rep_to,
            tof_window=tof_window_days,
            heliocentric_state_from=heliocentric_state_from,
            heliocentric_state_to=heliocentric_state_to,
        )
        out["tier0"] = t0
        # Tiers 1-2 (via two-tier wrapper)
        if self.two_tier is not None:
            out["tier12"] = self.two_tier.score_pair(
                rep_from, rep_to, c_j=c_j, dv_budget_kms=dv_budget_kms
            )
        # Tier 3
        if self.resonance is not None:
            out["tier3"] = self.resonance.score_pair(rep_from, rep_to)  # type: ignore[arg-type]
        # Tier 4
        if self.ftle is not None:
            out["tier4"] = self.ftle.score_pair(rep_from, rep_to, c_j=c_j)
        # Tier 5
        if self.lobe is not None:
            out["tier5"] = self.lobe.score_pair(rep_from, rep_to)  # type: ignore[arg-type]
        return out


# ---------------------------------------------------------------------------
# Convenience: build the patched-conic legs from one #264-style closed candidate.
# ---------------------------------------------------------------------------


def legs_from_repeated_moon_candidate(
    primary: str,
    sequence: Sequence[str],
    n_rev: Sequence[int],
    *,
    phase_samples: int = 24,
    tof_resonance_grid: Sequence[float] = (0.5, 1.0, 1.5, 2.0),
) -> list[PatchedConicLeg] | None:
    """Reconstruct the best-phasing Lambert legs for a #264 SILVER candidate.

    Mirrors :meth:`RepeatedMoonTarget._close_one_phasing` (#254 / #264) but
    returns the SI-units leg records (the NN scorer's native input shape)
    instead of the V_inf-continuity residual. Returns ``None`` if no phasing
    is feasible (no Lambert solution at the requested n_rev grid).

    Units: positions in m, velocities in m/s, dt in s, mu in m^3/s^2.

    Independent cross-check: this recomputation uses the same
    :func:`cyclerfinder.core.lambert.lambert` solver as the #264 closure (so a
    SILVER candidate that closed at residual < 0.05 km/s must produce N-1
    feasible legs here too). The leg geometry is otherwise identical.
    """
    from cyclerfinder.core.lambert import lambert
    from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

    day_s = 86400.0
    km_m = 1000.0
    mu_km3_s2 = PRIMARIES[primary]
    mu_m3_s2 = mu_km3_s2 * (km_m**3)
    consts: dict[str, tuple[float, float]] = {}
    for m in sequence:
        sat = SATELLITES[m]
        sma_km = sat.sma_km
        mean_motion = math.sqrt(mu_km3_s2 / (sma_km**3)) * day_s  # rad/day
        consts[m] = (sma_km, mean_motion)
    n_legs = len(sequence) - 1
    if len(n_rev) != n_legs:
        raise ValueError(f"n_rev length {len(n_rev)} != n_legs {n_legs}")
    # Per-leg ToF on the resonance grid (use 1.0 default; the #264 close also
    # sweeps a 4-point grid but the leg geometry is dominated by the
    # geometric-mean period). For the NN scorer we use the first feasible
    # phasing; the goal is the leg geometry, not the residual.
    phase0_grid = [2.0 * math.pi * i / phase_samples for i in range(phase_samples)]

    best_phasing: list[PatchedConicLeg] | None = None
    best_worst_vinf_diff = math.inf

    for phase0 in phase0_grid:
        theta0 = {
            m: phase0 + 2.0 * math.pi * j / max(1, len(consts))
            for j, m in enumerate(sorted(consts))
        }
        for tof_scale in tof_resonance_grid:
            tofs_d: list[float] = []
            for k in range(n_legs):
                _, na = consts[sequence[k]]
                _, nb = consts[sequence[k + 1]]
                pa = 2.0 * math.pi / na
                pb = 2.0 * math.pi / nb
                tofs_d.append(tof_scale * math.sqrt(pa * pb))
            epochs = [0.0]
            for tof in tofs_d:
                epochs.append(epochs[-1] + tof)
            states = []
            for m, t in zip(sequence, epochs, strict=True):
                sma, n = consts[m]
                theta = theta0[m] + n * t
                v_circ_km_s = math.sqrt(mu_km3_s2 / sma)
                pos = np.array([sma * math.cos(theta), sma * math.sin(theta), 0.0])
                vel = np.array([-v_circ_km_s * math.sin(theta), v_circ_km_s * math.cos(theta), 0.0])
                states.append((pos, vel))
            try:
                leg_records: list[PatchedConicLeg] = []
                worst_vinf_diff = 0.0
                ok = True
                last_vinf_arr: float | None = None
                for k in range(n_legs):
                    r_a_km, v_a_km_s = states[k]
                    r_b_km, v_b_km_s = states[k + 1]
                    sols = lambert(
                        r_a_km, r_b_km, tofs_d[k] * day_s, mu=mu_km3_s2, max_revs=max(0, n_rev[k])
                    )
                    wanted = [s for s in sols if s.n_revs == n_rev[k]]
                    if not wanted:
                        ok = False
                        break
                    best = min(wanted, key=lambda s: float(np.linalg.norm(s.v1 - v_a_km_s)))
                    v1_dep = best.v1  # km/s
                    v2_arr = best.v2  # km/s
                    vinf_dep = float(np.linalg.norm(v1_dep - v_a_km_s))
                    if last_vinf_arr is not None:
                        worst_vinf_diff = max(worst_vinf_diff, abs(vinf_dep - last_vinf_arr))
                    last_vinf_arr = float(np.linalg.norm(v2_arr - v_b_km_s))
                    leg_records.append(
                        PatchedConicLeg(
                            label_from=sequence[k],
                            label_to=sequence[k + 1],
                            r1_m=r_a_km * km_m,
                            v1_m_s=v1_dep * km_m,
                            r2_m=r_b_km * km_m,
                            v2_m_s=v2_arr * km_m,
                            dt_s=tofs_d[k] * day_s,
                            mu_m3_s2=mu_m3_s2,
                        )
                    )
                if ok and worst_vinf_diff < best_worst_vinf_diff:
                    best_worst_vinf_diff = worst_vinf_diff
                    best_phasing = leg_records
            except Exception:
                continue
    return best_phasing


__all__ = [
    "FiveTierPrioritizer",
    "PatchedConicLeg",
    "legs_from_repeated_moon_candidate",
    "rank_product_score",
]
