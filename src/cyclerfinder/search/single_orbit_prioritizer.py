"""Single-orbit prioritizer adapter (#310 Phase 1).

Architectural-gap closure for the #282 five-tier accessibility prioritizer.

The five-tier stack
(:class:`cyclerfinder.search.five_tier_prioritizer.FiveTierPrioritizer`) is
shaped around two input modes (see its module docstring):

* **Patched-conic leg mode** -- one Lambert leg `(r1, v1, r2, v2, dt, mu)`;
  Tier 0 (NN reachability) operates natively, Tiers 1-5 are skipped.
* **Representative-pair mode** -- an ordered PAIR of CR3BP representative
  orbits with `state0` + `period` (and Floquet structure for Tier 5);
  Tier 0 + Tiers 1-5 all run.

A **single-orbit discovery candidate** (one CR3BP IC + period + Jacobi, e.g.
the kind a #284 asymmetric-corrector scan or a #287-style 3D corrector spike
emits) fits NEITHER mode. The five-tier prioritizer has no compatible input
shape for it. That gap blocks any asymmetric-corrector batch scoring, any
sub-family-member ranking, any single-orbit shortlist.

This module is the adapter that closes the gap. It is **read-only** on
``five_tier_prioritizer.py`` and the underlying Phase 1 scorer modules; it
composes them on top of two adaptation strategies and reports an explicit
verdict on which tiers scored under which strategy.

ADAPTATION STRATEGIES (per the #284 architectural-gap analysis):

1. **Surrogate-pair neighbor**. Find the nearest catalogue / negative-
   registry pair-member to the single-orbit candidate in the tuple metric
   ``(mu_ratio, jacobi_C, period_nondim)`` -- a small, sortable, defensible
   distance that captures (system identity, energy level, characteristic
   timescale) without per-tier calibration. The surrogate becomes one slot
   of a representative pair; the candidate is the other slot. The full
   five-tier composite then runs on the (candidate, surrogate) pair. This
   gives Tier-0 / Tiers 1-2 / Tier 3 a meaningful input shape, at the cost
   of being a HYBRID accessibility verdict (`candidate -> surrogate` rather
   than `candidate -> candidate's_eventual_partner`).

2. **Parallel single-orbit pipeline**. Tier 4 (FTLE) and Tier 5 (Hiraiwa
   lobe-overlap) are scoring transport-corridor / phase-space-flux fields
   that DO NOT require a partner orbit -- both operate on a single
   representative's neighbourhood. FTLE samples the field along the geodesic
   from the candidate's planar position to itself (degenerate-same-point
   path; the scorer's existing branch reports the local cell value alone).
   Lobe-overlap builds the candidate's own lobe partition and reports the
   self-graph (no partner needed). Tiers 0-3 are marked None in this mode.

DEFAULT POLICY: try strategy 1 first; fall back to strategy 2 if no
surrogate within ``surrogate_neighbor_max_distance``. Both strategies are
first-class -- the result dataclass records which tiers scored under which
strategy via the ``notes`` field.

CROSS-CHECK -- the composition does not introduce any new physics; every
score traces verbatim to its source-paper scorer through the existing
five-tier API (Tier 0 -> :class:`NeuralReachPrefilter`; Tiers 1-2 ->
:class:`TwoTierPrioritizer`; Tier 3 -> :class:`ResonanceNetworkScorer`;
Tier 4 -> :class:`FTLEScorer`; Tier 5 -> :class:`LobeOverlapScorer`). The
adapter is the GLUE; the scoring physics is unchanged.

NO catalogue writeback. NO novelty claims. The module exposes scoring; the
caller decides what to do with the verdict.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.five_tier_prioritizer as ftp


def _as_finite_float(value: Any) -> float | None:
    """Coerce a dict-typed ``value`` to ``float`` if finite, else ``None``.

    Helper for the per-tier score extraction: scorer return dicts are
    typed ``dict[str, object]`` (the scoring submodules use ``object``
    values to admit mixed numeric / label / bool entries). The adapter
    only wants finite numerics; everything else (None, NaN, label
    strings, exceptions on float-cast) maps to ``None``.
    """
    if value is None:
        return None
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(f):
        return None
    return f


# ---------------------------------------------------------------------------
# Result dataclass.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SingleOrbitScore:
    """Verdict for one single-orbit candidate scored through the adapter.

    Attributes
    ----------
    candidate_id :
        Free-form identifier carried through to the result.
    surrogate_pair_neighbor_id :
        Catalogue (or negative-registry) id used as the pair-shape surrogate
        under strategy 1, or ``None`` if no neighbor was found and the
        parallel single-orbit pipeline (strategy 2) was used instead.
    tier_0_score :
        Tier 0 (NN reachability prefilter, Zhang-Topputo 2026) -- the
        predicted patched-conic Lambert-leg ΔV in km/s. ``None`` if the
        surrogate-pair strategy was unavailable AND the candidate carries
        no heliocentric SI state to drive the NN directly.
    tier_1_score :
        Tier 1 (Braik-Ross 2026 energy-preserving heading-fan overlap).
        ``None`` if no surrogate neighbor.
    tier_2_score :
        Tier 2 (Zhou-Armellin 2024 single-impulse footprint) -- minimum ΔV
        in km/s. ``None`` if no surrogate neighbor.
    tier_3_score :
        Tier 3 (Kumar-Rawat-Rosengren-Ross 2025 resonant-manifold
        heteroclinic perigee distance). ``None`` if no surrogate neighbor.
    tier_4_score :
        Tier 4 (Canales-Howell 2023 FTLE transport-corridor strength).
        Always populated when an :class:`FTLEScorer` is provided -- the
        single-orbit case is the degenerate-same-point geodesic the scorer
        already handles.
    tier_5_score :
        Tier 5 (Hiraiwa-Bando-Sato-Hokamoto 2026 lobe-overlap bottleneck
        flux). Populated when a :class:`LobeOverlapScorer` is provided
        AND the candidate's Floquet structure was successfully computed.
    combined_rank :
        Rank-product composite over the tiers that scored. ``None`` if
        fewer than two tiers populated (rank-product on one tier is
        trivially the rank itself; the value carries no composition
        signal).
    notes :
        Free-form narrative -- which strategy was used, why fallback
        triggered, any per-tier exceptions.
    """

    candidate_id: str
    surrogate_pair_neighbor_id: str | None
    tier_0_score: float | None
    tier_1_score: float | None
    tier_2_score: float | None
    tier_3_score: float | None
    tier_4_score: float | None
    tier_5_score: float | None
    combined_rank: float | None
    notes: str = ""


# ---------------------------------------------------------------------------
# Surrogate-neighbor lookup -- (mu, Jacobi, period) tuple-distance.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class _SurrogateCandidate:
    """One pair-shape candidate from the catalogue or neg-registry.

    Carries the minimal fields the tuple-distance metric needs PLUS the
    state / period that the prioritizer will substitute into the
    representative-pair slot.
    """

    id: str
    mu: float
    jacobi: float
    period_nd: float
    state_nd: NDArray[np.float64]


def _surrogate_distance(
    candidate_mu: float,
    candidate_jacobi: float,
    candidate_period_nondim: float,
    surrogate: _SurrogateCandidate,
) -> float:
    """Tuple-distance between candidate and surrogate over (mu, C, T).

    The metric is the L2 norm of the relative differences, with a guard
    against division by zero. Mu is the system-identity coordinate (two
    candidates at different mu live in DIFFERENT CR3BP systems and are
    not meaningfully comparable as a representative pair); Jacobi is the
    energy-level coordinate (defines the same-C manifold); period is the
    characteristic-timescale coordinate (a same-mu / same-C orbit with
    wildly different period is a different family member).

    Defensible-not-sourced: no source paper publishes a single-orbit-to-
    pair surrogate-distance recipe (the question doesn't arise in the
    paper batch). The tuple choice is the smallest physically-meaningful
    descriptor; the L2-with-relative-differences norm is the standard
    multi-axis distance when the axes have different units. The
    alternative -- a weighted sum with per-axis calibration -- requires
    a sourced calibration that does not exist.
    """
    mu_rel = abs(candidate_mu - surrogate.mu) / max(abs(candidate_mu), 1e-12)
    c_rel = abs(candidate_jacobi - surrogate.jacobi) / max(abs(candidate_jacobi), 1e-12)
    t_rel = abs(candidate_period_nondim - surrogate.period_nd) / max(
        abs(candidate_period_nondim), 1e-12
    )
    return float(math.sqrt(mu_rel * mu_rel + c_rel * c_rel + t_rel * t_rel))


def _load_catalogue_surrogates(
    catalogue_path: Path | None = None,
) -> list[_SurrogateCandidate]:
    """Scan ``data/catalogue.yaml`` for rows usable as pair-shape surrogates.

    A row qualifies iff its ``orbit_elements.cr3bp`` block carries
    ``mass_ratio``, ``jacobi_constant``, ``period_nd``, and ``state_nd``
    (a non-null 6-vector). Rows missing any field are skipped silently --
    they simply don't participate in the neighbor metric.

    Read-only. Does NOT mutate the catalogue.
    """
    import yaml  # type: ignore[import-untyped]  # local import; yaml is a dev dep

    if catalogue_path is None:
        catalogue_path = Path(__file__).resolve().parents[3] / "data" / "catalogue.yaml"
    if not catalogue_path.exists():
        return []
    raw = yaml.safe_load(catalogue_path.read_text())
    out: list[_SurrogateCandidate] = []
    for row in raw or []:
        if not isinstance(row, dict):
            continue
        oe = row.get("orbit_elements") or {}
        cr = oe.get("cr3bp") or {}
        mu = cr.get("mass_ratio")
        c_j = cr.get("jacobi_constant")
        period_nd = cr.get("period_nd")
        state_nd = cr.get("state_nd")
        row_id = row.get("id")
        if (
            mu is None
            or c_j is None
            or period_nd is None
            or state_nd is None
            or not isinstance(state_nd, (list, tuple))
            or len(state_nd) != 6
            or not isinstance(row_id, str)
        ):
            continue
        try:
            state_arr = np.asarray(state_nd, dtype=np.float64)
        except (TypeError, ValueError):
            continue
        if state_arr.shape != (6,) or not np.all(np.isfinite(state_arr)):
            continue
        out.append(
            _SurrogateCandidate(
                id=row_id,
                mu=float(mu),
                jacobi=float(c_j),
                period_nd=float(period_nd),
                state_nd=state_arr,
            )
        )
    return out


def find_surrogate_neighbor(
    candidate_mu: float,
    candidate_jacobi: float,
    candidate_period_nondim: float,
    *,
    max_distance: float | None = None,
    catalogue_path: Path | None = None,
    extra_surrogates: list[_SurrogateCandidate] | None = None,
) -> tuple[_SurrogateCandidate, float] | None:
    """Find the nearest pair-shape surrogate; ``None`` if outside ``max_distance``.

    The metric is :func:`_surrogate_distance` over the catalogue rows that
    expose ``(mass_ratio, jacobi_constant, period_nd, state_nd)``. Additional
    pair-shape candidates can be supplied via ``extra_surrogates`` (e.g. from
    a negative-registry; the registry is not loaded by this module to keep
    the import surface small, but the seam is open).
    """
    catalogue = _load_catalogue_surrogates(catalogue_path)
    if extra_surrogates:
        catalogue = list(catalogue) + list(extra_surrogates)
    if not catalogue:
        return None
    best: tuple[_SurrogateCandidate, float] | None = None
    for s in catalogue:
        d = _surrogate_distance(candidate_mu, candidate_jacobi, candidate_period_nondim, s)
        if best is None or d < best[1]:
            best = (s, d)
    if best is None:
        return None
    if max_distance is not None and best[1] > max_distance:
        return None
    return best


# ---------------------------------------------------------------------------
# A minimal "RepView with everything the scorers need" duck-typed object.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class _AdapterRep:
    """Duck-typed representative for the FTLE/two-tier/Tier-0 scorers.

    Carries ``label``, ``state0``, ``period`` -- the minimum the
    :class:`FTLEScorer` and :class:`TwoTierPrioritizer` accept (both
    expressly support any object exposing those three attributes, see the
    duck-typed branch in :func:`ftle_scorer.FTLEScorer.score_pair` and
    :func:`two_tier_prioritizer._as_repview`).
    """

    label: str
    state0: NDArray[np.float64]
    period: float


# ---------------------------------------------------------------------------
# Single-orbit scoring API.
# ---------------------------------------------------------------------------


def score_single_orbit(
    candidate_state0: NDArray[np.float64],
    candidate_period_nondim: float,
    candidate_system: cr3bp.CR3BPSystem,
    *,
    candidate_id: str = "single-orbit-candidate",
    candidate_jacobi: float | None = None,
    prioritizer: ftp.FiveTierPrioritizer | None = None,
    use_surrogate_neighbor: bool = True,
    surrogate_neighbor_max_distance: float | None = None,
    catalogue_path: Path | None = None,
    extra_surrogates: list[_SurrogateCandidate] | None = None,
    skip_tiers: tuple[int, ...] = (),
    notes: str = "",
) -> SingleOrbitScore:
    """Adapt the #282 five-tier prioritizer to score a single-orbit candidate.

    Parameters
    ----------
    candidate_state0 :
        6D CR3BP IC ``(x, y, z, xdot, ydot, zdot)`` in nondimensional
        rotating-frame coordinates.
    candidate_period_nondim :
        Period in CR3BP nondimensional time units (TU).
    candidate_system :
        The CR3BP system the candidate lives in. Used to recompute the
        Jacobi (if not supplied) and as the integration system for the
        Floquet structure (Tier 5).
    candidate_id :
        Free-form identifier carried through to the result.
    candidate_jacobi :
        Optional pre-computed Jacobi constant. If ``None``, recomputed
        from ``candidate_state0`` + ``candidate_system.mu``.
    prioritizer :
        Optional :class:`FiveTierPrioritizer` to drive Tier 0 + Tiers 1-5.
        If ``None`` and Tiers 0-3 are not skipped, a default prioritizer
        is built with no tier-1/2/3/4/5 scorers attached (so only the
        Tier-0 NN runs in surrogate mode). The caller is expected to wire
        whatever scorers it wants reachable; this adapter does not build
        them by default because each (esp. tier 4 FTLE field) is
        configuration-heavy.
    use_surrogate_neighbor :
        If True (default), try to find the nearest pair-shape surrogate
        and run the representative-pair pipeline against it. If False, go
        straight to the parallel single-orbit pipeline (FTLE + lobe only).
    surrogate_neighbor_max_distance :
        Maximum tuple-distance for an admissible surrogate. ``None`` =
        accept any (the nearest catalogue row regardless of distance).
    catalogue_path :
        Override for :data:`data/catalogue.yaml`.
    extra_surrogates :
        Optional additional pair-shape candidates (e.g. from a negative
        registry the caller has loaded).
    skip_tiers :
        Tiers to skip explicitly (return ``None`` for them regardless of
        availability). Use e.g. ``(0, 1, 2, 3)`` to score Tier 4 + Tier 5
        alone (the parallel-pipeline-only verdict).
    notes :
        Extra prose appended to the result's ``notes`` field.

    Returns
    -------
    SingleOrbitScore
        Dataclass with per-tier verdicts + the rank-product composite over
        the tiers that actually scored.
    """
    skip = set(int(t) for t in skip_tiers)
    notes_parts: list[str] = []
    if notes:
        notes_parts.append(notes)

    # Compute / accept Jacobi.
    if candidate_jacobi is None:
        candidate_jacobi = cr3bp.jacobi_constant(candidate_state0, candidate_system.mu)

    if prioritizer is None:
        prioritizer = ftp.FiveTierPrioritizer()

    # --- Strategy 1: surrogate-pair neighbor (tiers 0-3) ---
    surrogate_id: str | None = None
    surrogate: _SurrogateCandidate | None = None
    surrogate_distance: float | None = None
    if use_surrogate_neighbor and not skip.issuperset({0, 1, 2, 3}):
        found = find_surrogate_neighbor(
            candidate_system.mu,
            float(candidate_jacobi),
            float(candidate_period_nondim),
            max_distance=surrogate_neighbor_max_distance,
            catalogue_path=catalogue_path,
            extra_surrogates=extra_surrogates,
        )
        if found is not None:
            surrogate, surrogate_distance = found
            surrogate_id = surrogate.id
            notes_parts.append(
                f"strategy-1 surrogate-neighbor: id={surrogate.id} "
                f"distance={surrogate_distance:.3e}"
            )
        else:
            notes_parts.append(
                "strategy-1 surrogate-neighbor: no catalogue row within "
                f"max_distance={surrogate_neighbor_max_distance}"
            )

    candidate_rep = _AdapterRep(
        label=candidate_id,
        state0=np.asarray(candidate_state0, dtype=np.float64),
        period=float(candidate_period_nondim),
    )

    # Per-tier scores; populated below.
    t0: float | None = None
    t1: float | None = None
    t2: float | None = None
    t3: float | None = None
    t4: float | None = None
    t5: float | None = None

    # ------------------------------------------------------------------
    # Strategy 1: representative-pair against surrogate.
    # ------------------------------------------------------------------
    if surrogate is not None:
        surrogate_rep = _AdapterRep(
            label=surrogate.id,
            state0=surrogate.state_nd,
            period=surrogate.period_nd,
        )
        # Tier 0 -- NN ΔV via the prioritizer's score_pair entry point.
        if 0 not in skip and prioritizer.nn_prefilter is not None:
            try:
                t0_dict = prioritizer.nn_prefilter.score_pair(candidate_rep, surrogate_rep)
                t0 = _as_finite_float(t0_dict.get("predicted_dv_kms"))
            except Exception as e:
                notes_parts.append(f"tier0 raised: {e!r}")
        # Tier 1-2 -- two-tier prioritizer.
        if (1 not in skip or 2 not in skip) and prioritizer.two_tier is not None:
            try:
                t12 = prioritizer.two_tier.score_pair(
                    candidate_rep, surrogate_rep, c_j=float(candidate_jacobi)
                )
                if 1 not in skip:
                    t1 = _as_finite_float(t12.get("tier1_heading_overlap"))
                if 2 not in skip:
                    t2 = _as_finite_float(t12.get("tier2_impulsive_min_dv_kms"))
            except Exception as e:
                notes_parts.append(f"tier12 raised: {e!r}")
        # Tier 3 -- resonance-network scorer. The two-tier scorer expects a
        # :class:`resonance_network.ResonantMember` (it indexes per-orbit
        # Floquet structure); we duck-type via ``_AdapterRep`` and rely on
        # the scorer's documented ``getattr`` fallback (the scorer ignores
        # missing Floquet keys and reports a degraded verdict). The type
        # mismatch is intentional and bracketed by a try/except.
        if 3 not in skip and prioritizer.resonance is not None:
            try:
                t3_dict = prioritizer.resonance.score_pair(
                    candidate_rep,  # type: ignore[arg-type]
                    surrogate_rep,  # type: ignore[arg-type]
                )
                # Resonance-network scoring exposes a perigee-distance scalar; the
                # exact key name depends on the scorer implementation. We probe a
                # set of documented keys and take the first scalar that finite-
                # parses. This is the only sensible value to compare for tier
                # 3 across the candidate pair (lower-is-better).
                if isinstance(t3_dict, dict):
                    for key in ("min_perigee_distance", "perigee_distance", "min_distance"):
                        candidate_t3 = _as_finite_float(t3_dict.get(key))
                        if candidate_t3 is not None:
                            t3 = candidate_t3
                            break
            except Exception as e:
                notes_parts.append(f"tier3 raised: {e!r}")

    # ------------------------------------------------------------------
    # Strategy 2: parallel single-orbit pipeline (Tier 4 + Tier 5).
    #
    # Both run even when strategy 1 succeeded -- they DO NOT require the
    # surrogate. The "fallback" framing in the docstring is only relevant
    # for tiers 0-3 (which are SKIPPED when no surrogate is found).
    # ------------------------------------------------------------------

    # Tier 4 -- FTLE on the candidate's planar position (degenerate-same-point
    # geodesic; the scorer's existing branch reports the local cell value).
    if 4 not in skip and prioritizer.ftle is not None:
        try:
            t4_dict = prioritizer.ftle.score_pair(candidate_rep, candidate_rep)
            t4 = _as_finite_float(t4_dict.get("transport_corridor_strength"))
        except Exception as e:
            notes_parts.append(f"tier4 raised: {e!r}")

    # Tier 5 -- lobe-overlap requires a ResonantMember (Floquet pair). Build
    # one synthetically from the candidate's state + period.
    if 5 not in skip and prioritizer.lobe is not None:
        try:
            # Local import: resonance_network drags in heavy CR3BP machinery;
            # don't force it on importers of this adapter when tier 5 is
            # skipped.
            import cyclerfinder.search.resonance_network as rn

            lam, vec = rn._planar_floquet(
                candidate_system,
                np.asarray(candidate_state0, dtype=np.float64),
                float(candidate_period_nondim),
            )
            # Build a synthetic ResonantMember. The scorer reads label,
            # state0, period, unstable_eigenvector; the other fields are
            # bookkeeping (period_days, sourced_period_days, confirmed) and
            # do not influence the lobe partition.
            member = rn.ResonantMember(
                label=candidate_id,
                state0=np.asarray(candidate_state0, dtype=np.float64),
                period=float(candidate_period_nondim),
                jacobi=float(candidate_jacobi),
                sourced_period_days=float("nan"),
                period_days=float(candidate_period_nondim) * rn.TU_DAYS,
                confirmed=False,
                unstable_eigenvalue=float(lam),
                unstable_eigenvector=vec,
            )
            t5_dict = prioritizer.lobe.score_pair(member, member)
            t5 = _as_finite_float(t5_dict.get("min_path_flux"))
        except Exception as e:
            notes_parts.append(f"tier5 raised: {e!r}")

    # ------------------------------------------------------------------
    # Rank-product composite -- over the tiers that scored.
    # ------------------------------------------------------------------
    scored = [
        (idx, val, descending)
        for idx, val, descending in (
            (0, t0, False),  # ΔV, lower-is-better
            (1, t1, True),  # heading-overlap, higher-is-better
            (2, t2, False),  # impulsive ΔV, lower-is-better
            (3, t3, False),  # perigee distance, lower-is-better
            (4, t4, True),  # corridor strength, higher-is-better
            (5, t5, True),  # bottleneck flux, higher-is-better
        )
        if val is not None and math.isfinite(val)
    ]
    combined_rank: float | None = None
    if len(scored) >= 2:
        # Single-candidate rank-product: every tier ranks the lone
        # candidate at rank 1, so the geometric mean is 1.0. That's the
        # honest output -- the composite is meaningful only when the
        # adapter is called across MANY candidates (the caller composes).
        # We still expose the per-tier values; the caller can re-run
        # rank-product across its candidate set.
        combined_rank = 1.0

    if surrogate_id is not None:
        notes_parts.append(
            f"strategy-1 representative-pair tiers: "
            f"t0={'ok' if t0 is not None else 'skip'} "
            f"t1={'ok' if t1 is not None else 'skip'} "
            f"t2={'ok' if t2 is not None else 'skip'} "
            f"t3={'ok' if t3 is not None else 'skip'}"
        )
    notes_parts.append(
        f"strategy-2 parallel-pipeline tiers: "
        f"t4={'ok' if t4 is not None else 'skip'} "
        f"t5={'ok' if t5 is not None else 'skip'}"
    )

    return SingleOrbitScore(
        candidate_id=candidate_id,
        surrogate_pair_neighbor_id=surrogate_id,
        tier_0_score=t0,
        tier_1_score=t1,
        tier_2_score=t2,
        tier_3_score=t3,
        tier_4_score=t4,
        tier_5_score=t5,
        combined_rank=combined_rank,
        notes=" | ".join(notes_parts),
    )


__all__ = [
    "SingleOrbitScore",
    "_AdapterRep",
    "_SurrogateCandidate",
    "find_surrogate_neighbor",
    "score_single_orbit",
]
