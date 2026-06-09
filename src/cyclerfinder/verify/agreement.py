"""Axis A — code-path agreement cross-check (the Forge, Phase 2).

For a *solved or constructed* cycler, require **>= 2 independent code
paths** to agree on the per-encounter hyperbolic-excess speeds / geometry
before the candidate is trusted. This is the "independence-at-every-tier"
governing principle of the Forge design spec
(``docs/superpowers/specs/2026-06-03-the-forge-discovery-pipeline-design.md``)
realised as a single combiner, :func:`crosscheck_code_paths`, that emits a
frozen :class:`AgreementReport` for Phase 3's gauntlet to consume.

Three paths, each genuinely independent of the others
-----------------------------------------------------
(a) **In-house Lambert vs lamberthub** — every leg re-solved with
    ``lamberthub``'s ``izzo2015`` + ``gooding1990`` and compared to the
    in-house :func:`cyclerfinder.core.lambert.lambert`. This is M7's
    spec §14 V1 check; it is *reused verbatim* from
    :mod:`cyclerfinder.verify.crosscheck` (never reimplemented here).

(b) **Resonance construction vs the optimiser / Lambert geometry** — the
    closed-form circular-coplanar V_inf from
    :func:`cyclerfinder.search.resonant_construct.construct_resonant_cycler`
    (built analytically from the cycler's own heliocentric ``(a, e)``)
    compared, per encounter, to the V_inf the cycler actually carries
    (which was produced by the *Lambert* boundary-value solver in
    :mod:`cyclerfinder.search.construct`, or by the free optimiser in
    :mod:`cyclerfinder.search.optimize`). Two completely different
    derivations of the same physical quantity: a root-find on the
    Lambert time-of-flight equation vs a closed-form vis-viva crossing.
    When a :class:`~cyclerfinder.search.sequence.Cell` and an ephemeris
    are supplied the free optimiser
    (:func:`~cyclerfinder.search.optimize.optimise_cell_idealized`) is
    *additionally* run and its converged V_inf folded in as a third
    independent witness; if it does not converge to a feasible closing
    geometry (the documented circular-coplanar scope boundary for the
    E-M family — see ``tests/search/test_optimize.py``
    ``test_2syn_em_rediscovers_5_65_kms_earth``, xfail) that witness is
    recorded as *unavailable*, not as a failure, so the construction-vs-
    construction agreement still gates.

(c) **Forward Kepler re-propagation residual** — each leg's departure
    state ``(r_planet, v_planet + vinf_out)`` is propagated forward by the
    leg's time-of-flight with :func:`cyclerfinder.core.kepler.propagate`
    (the universal-variable two-body propagator — *not* the Lambert that
    built the leg) and the resulting position is compared to the arrival
    encounter's claimed planet position. A well-constructed cycler closes
    to metres; a corrupted state diverges by millions of km.

Golden discipline
-----------------
No path asserts a computed value against a fabricated constant. Paths
(a) and (b) are *consistency predicates* — agreement between two
independent in-house code paths. Path (c) is a physical self-consistency
residual (forward propagation must return to the geometry the cycler
claims). The sourced Aldrin V_inf in the catalogue YAML may be asserted
*additionally* by a caller, because it is independently sourced; this
module does not bake it in.

Tolerances are module constants, spec-derived, and NOT test-tunable.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

import numpy as np

from cyclerfinder.core.constants import MU_SUN_KM3_S2, PLANETS
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.kepler import KeplerConvergenceError, propagate
from cyclerfinder.model.cycler import Cycler, orbit_elements_au
from cyclerfinder.search.resonant_construct import construct_resonant_cycler
from cyclerfinder.verify.crosscheck import (
    V1_TOLERANCE_MPS,
    LambertCrosscheckResult,
    crosscheck_cycler,
)

# ---------------------------------------------------------------------------
# Tolerances (spec-derived; NOT test-tunable)
# ---------------------------------------------------------------------------

VINF_AGREEMENT_TOL_KMS: Final[float] = 0.05
"""Per-encounter V_inf agreement bound between two independent in-house
derivations, km/s. The Lambert construction and the closed-form resonance
construction agree to ~1e-13 km/s when fed the same heliocentric
``(a, e)``; the optimiser, when it converges, lands on the same family to
within its SLSQP/DE polish. 0.05 km/s (50 m/s) is two orders of magnitude
tighter than the ±0.3 km/s inter-source rounding the catalogue
rediscovery harness absorbs, yet loose enough to admit the optimiser's
finite polish — it gates *code-path agreement*, not source fidelity."""

KEPLER_REPROP_TOL_KM: Final[float] = 1.0
"""Forward Kepler re-propagation position residual bound, km. A
circular-coplanar Aldrin seed re-propagates to ~2e-4 km; a 1 % velocity
perturbation diverges to ~8e6 km. 1 km sits ~4 orders of magnitude above
the clean residual and ~7 below the corrupted one — it rejects a
propagator regression or a corrupted state without tripping on
Lambert/propagator numerical noise."""


# ---------------------------------------------------------------------------
# Per-path result records
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class LamberthubPathResult:
    """Path (a): in-house Lambert vs lamberthub agreement.

    Attributes
    ----------
    available:
        ``True`` when the path ran (it always does for a cycler with
        legs and an ephemeris).
    per_leg:
        The reused :class:`~cyclerfinder.verify.crosscheck.LambertCrosscheckResult`
        records, one per leg.
    max_diff_mps:
        Worst per-leg in-house-vs-lamberthub disagreement, m/s.
    passed:
        ``available and max_diff_mps < V1_TOLERANCE_MPS``.
    """

    available: bool
    per_leg: tuple[LambertCrosscheckResult, ...]
    max_diff_mps: float
    passed: bool


@dataclass(frozen=True)
class ConstructionOptimiserPathResult:
    """Path (b): resonance construction vs the cycler's (and optionally the
    free optimiser's) V_inf.

    Attributes
    ----------
    available:
        ``True`` when the resonance construction reached every body the
        cycler visits (i.e. the cycler's heliocentric ``(a, e)`` spans
        each body's circular radius). ``False`` if the analytic crossing
        is undefined for some body.
    resonant_vinf_kms:
        Closed-form V_inf per body code from
        :func:`construct_resonant_cycler`.
    cycler_vinf_kms:
        The cycler's own per-body V_inf magnitude (``max`` of in/out).
    construction_max_diff_kms:
        Worst per-body ``|resonant - cycler|``, km/s.
    optimiser_available:
        ``True`` iff a cell + ephem were supplied AND the optimiser
        converged to a feasible closing geometry. ``False`` records the
        documented circular-coplanar scope boundary (the optimiser does
        not close the E-M family in the idealised model) — *not* a
        failure of this path.
    optimiser_vinf_kms:
        The optimiser's per-body V_inf when ``optimiser_available``, else
        an empty mapping.
    optimiser_max_diff_kms:
        Worst per-body ``|resonant - optimiser|`` km/s when available,
        else ``None``.
    max_diff_kms:
        The binding disagreement for the gate: the construction diff, and
        the optimiser diff when that witness is available.
    passed:
        ``available and max_diff_kms < VINF_AGREEMENT_TOL_KMS``.
    """

    available: bool
    resonant_vinf_kms: dict[str, float]
    cycler_vinf_kms: dict[str, float]
    construction_max_diff_kms: float
    optimiser_available: bool
    optimiser_vinf_kms: dict[str, float]
    optimiser_max_diff_kms: float | None
    max_diff_kms: float
    passed: bool


@dataclass(frozen=True)
class KeplerRepropPathResult:
    """Path (c): forward Kepler re-propagation residual.

    Attributes
    ----------
    available:
        ``True`` when every leg propagated without a
        :class:`~cyclerfinder.core.kepler.KeplerConvergenceError`.
    per_leg_residual_km:
        Position residual ``||r_propagated - r_arrival_claimed||`` per
        leg, km, in leg order.
    max_residual_km:
        Worst per-leg residual, km (``inf`` if a leg failed to
        propagate).
    passed:
        ``available and max_residual_km < KEPLER_REPROP_TOL_KM``.
    """

    available: bool
    per_leg_residual_km: tuple[float, ...]
    max_residual_km: float
    passed: bool


@dataclass(frozen=True)
class AgreementReport:
    """Frozen Axis-A verdict — the combiner output Phase 3 consumes.

    Attributes
    ----------
    lamberthub:
        Path (a) result.
    construction_optimiser:
        Path (b) result.
    kepler_reprop:
        Path (c) result.
    n_paths_available:
        How many of the three paths actually ran.
    n_paths_passed:
        How many available paths passed their tolerance.
    agreed:
        Overall verdict — ``True`` iff **at least two** paths are
        available AND every available path passed. The ">= 2 independent
        code paths agree" predicate of the Forge spec: a single passing
        path is never sufficient, and an available-but-failing path vetoes
        the verdict.
    """

    lamberthub: LamberthubPathResult
    construction_optimiser: ConstructionOptimiserPathResult
    kepler_reprop: KeplerRepropPathResult
    n_paths_available: int
    n_paths_passed: int
    agreed: bool


# ---------------------------------------------------------------------------
# Per-path implementations
# ---------------------------------------------------------------------------


def _cycler_vinf_by_body(cycler: Cycler) -> dict[str, float]:
    """Map ``body -> max(||vinf_in||, ||vinf_out||)`` over the cycler.

    ``max`` (not the in/out average) is the binding magnitude when a body
    appears at more than one encounter, matching the catalogue
    rediscovery harness's ``_vinf_magnitudes_by_body``.
    """
    out: dict[str, float] = {}
    for enc in cycler.encounters:
        m = max(
            float(np.linalg.norm(enc.vinf_in)),
            float(np.linalg.norm(enc.vinf_out)),
        )
        out[enc.body] = max(out.get(enc.body, 0.0), m)
    return out


def _run_lamberthub_path(
    cycler: Cycler,
    ephem: Ephemeris,
    *,
    mu: float,
) -> LamberthubPathResult:
    """Path (a): reuse :func:`crosscheck_cycler` (never reimplemented)."""
    if not cycler.legs:
        return LamberthubPathResult(
            available=False,
            per_leg=(),
            max_diff_mps=0.0,
            passed=False,
        )
    per_leg = crosscheck_cycler(cycler, ephem, mu=mu)
    max_diff = max((r.max_diff_mps for r in per_leg), default=0.0)
    return LamberthubPathResult(
        available=True,
        per_leg=per_leg,
        max_diff_mps=max_diff,
        passed=max_diff < V1_TOLERANCE_MPS,
    )


def _resonant_vinf_for_cycler(
    cycler: Cycler,
    *,
    mu: float,
) -> dict[str, float] | None:
    """Closed-form resonance-construction V_inf for the cycler's geometry.

    Derives the heliocentric ``(a, e)`` from the cycler's first leg
    departure state (vis-viva + eccentricity vector via
    :func:`cyclerfinder.model.cycler.orbit_elements_au`) and runs
    :func:`construct_resonant_cycler` over the *unique* body codes the
    cycler visits. Returns ``None`` if the orbit does not span some body's
    radius (the analytic crossing is then undefined) — recorded as the
    path being unavailable, not a failure.
    """
    if not cycler.legs:
        return None
    first_leg = cycler.legs[0]
    a_au, e = orbit_elements_au(cycler.encounters[0].r, first_leg.v_depart, mu=mu)
    if not (0.0 <= e < 1.0) or a_au <= 0.0:
        # Hyperbolic / degenerate leg — the elliptic crossing construction
        # does not apply.
        return None
    bodies: list[str] = []
    for enc in cycler.encounters:
        if enc.body not in bodies:
            bodies.append(enc.body)
    # The closed-form resonance crossing is a HELIOCENTRIC single-ellipse model
    # (``construct_resonant_cycler`` resolves body radii via ``PLANETS``). It does
    # not apply to a planet-centred moon tour (moon codes live in ``SATELLITES``),
    # so the path is recorded UNAVAILABLE for such a candidate rather than
    # crashing — exactly the design-note §2 "resonance path dropped" rule, here
    # for a non-heliocentric centre. lamberthub + Kepler-reprop remain the two
    # independent Axis-A witnesses.
    if any(b not in PLANETS for b in bodies):
        return None
    try:
        rc = construct_resonant_cycler(a_au, e, bodies=tuple(bodies), mu=mu)
    except ValueError:
        return None
    return dict(rc.vinf_kms)


def _run_construction_optimiser_path(
    cycler: Cycler,
    *,
    mu: float,
    optimiser_vinf: dict[str, float] | None,
) -> ConstructionOptimiserPathResult:
    """Path (b): resonance construction vs the cycler's (+ optimiser) V_inf."""
    resonant = _resonant_vinf_for_cycler(cycler, mu=mu)
    cycler_vinf = _cycler_vinf_by_body(cycler)
    if resonant is None:
        return ConstructionOptimiserPathResult(
            available=False,
            resonant_vinf_kms={},
            cycler_vinf_kms=cycler_vinf,
            construction_max_diff_kms=float("inf"),
            optimiser_available=False,
            optimiser_vinf_kms={},
            optimiser_max_diff_kms=None,
            max_diff_kms=float("inf"),
            passed=False,
        )

    construction_diff = max(
        (abs(resonant[b] - cycler_vinf[b]) for b in resonant if b in cycler_vinf),
        default=float("inf"),
    )

    optimiser_available = False
    optimiser_diff: float | None = None
    opt_vinf: dict[str, float] = {}
    if optimiser_vinf:
        opt_vinf = optimiser_vinf
        optimiser_available = True
        optimiser_diff = max(
            (abs(resonant[b] - opt_vinf[b]) for b in resonant if b in opt_vinf),
            default=float("inf"),
        )

    diffs = [construction_diff]
    if optimiser_diff is not None:
        diffs.append(optimiser_diff)
    max_diff = max(diffs)

    return ConstructionOptimiserPathResult(
        available=True,
        resonant_vinf_kms=resonant,
        cycler_vinf_kms=cycler_vinf,
        construction_max_diff_kms=construction_diff,
        optimiser_available=optimiser_available,
        optimiser_vinf_kms=opt_vinf,
        optimiser_max_diff_kms=optimiser_diff,
        max_diff_kms=max_diff,
        passed=max_diff < VINF_AGREEMENT_TOL_KMS,
    )


def _run_kepler_reprop_path(cycler: Cycler, *, mu: float) -> KeplerRepropPathResult:
    """Path (c): forward Kepler re-propagation residual per leg.

    Each leg departs from ``(encounter.r, v_planet + vinf_out)`` — i.e.
    the spacecraft heliocentric state at departure, reconstructed from the
    cycler's *encoded* geometry — and is propagated forward by the leg
    time-of-flight with :func:`cyclerfinder.core.kepler.propagate`. The
    residual is the distance from the propagated position to the arrival
    encounter's claimed planet position.
    """
    if not cycler.legs:
        return KeplerRepropPathResult(
            available=False,
            per_leg_residual_km=(),
            max_residual_km=float("inf"),
            passed=False,
        )

    # Index encounters by (body, time) so a leg's endpoints resolve to the
    # encounter records that carry the planet states and V_inf vectors.
    enc_by_key: dict[tuple[str, float], int] = {
        (enc.body, enc.t): i for i, enc in enumerate(cycler.encounters)
    }

    residuals: list[float] = []
    available = True
    for leg in cycler.legs:
        dep_i = enc_by_key.get((leg.from_body, leg.t_depart))
        arr_i = enc_by_key.get((leg.to_body, leg.t_arrive))
        if dep_i is None or arr_i is None:
            # Malformed cycler: a leg whose epochs match no encounter.
            available = False
            residuals.append(float("inf"))
            continue
        dep_enc = cycler.encounters[dep_i]
        arr_enc = cycler.encounters[arr_i]
        r0 = np.asarray(dep_enc.r, dtype=np.float64)
        v_sc = np.asarray(dep_enc.v_planet, dtype=np.float64) + np.asarray(
            dep_enc.vinf_out, dtype=np.float64
        )
        dt = leg.t_arrive - leg.t_depart
        try:
            r_prop, _v_prop = propagate(r0, v_sc, dt, mu=mu)
        except KeplerConvergenceError:
            available = False
            residuals.append(float("inf"))
            continue
        residuals.append(float(np.linalg.norm(r_prop - np.asarray(arr_enc.r, dtype=np.float64))))

    max_res = max(residuals, default=float("inf"))
    return KeplerRepropPathResult(
        available=available,
        per_leg_residual_km=tuple(residuals),
        max_residual_km=max_res,
        passed=available and max_res < KEPLER_REPROP_TOL_KM,
    )


# ---------------------------------------------------------------------------
# Combiner
# ---------------------------------------------------------------------------


def crosscheck_code_paths(
    cycler: Cycler,
    ephem: Ephemeris | None = None,
    *,
    mu: float = MU_SUN_KM3_S2,
    optimiser_vinf_kms: dict[str, float] | None = None,
) -> AgreementReport:
    """Run the available Axis-A code paths and combine into an :class:`AgreementReport`.

    Parameters
    ----------
    cycler:
        The solved / constructed cycler under cross-check.
    ephem:
        Heliocentric state provider. Required for path (a) (the lamberthub
        cross-check); ``None`` runs only paths (b) and (c), which need
        only the cycler's own encoded geometry.
    mu:
        Heliocentric gravitational parameter, km^3/s^2. Defaults to
        :data:`cyclerfinder.core.constants.MU_SUN_KM3_S2`.
    optimiser_vinf_kms:
        Optional per-body V_inf from a *free* optimiser run on the same
        cell (e.g.
        :func:`cyclerfinder.search.optimize.optimise_cell_idealized`'s
        result). When supplied it is folded into path (b) as a third
        independent witness. Omitted (``None``) leaves path (b) gating on
        the construction-vs-construction agreement alone — the correct
        behaviour for the E-M family, whose idealised optimiser does not
        converge to a feasible closing geometry (a documented scope
        boundary, not a defect). Pass it through from the optimiser's
        result *only* when ``converged and constraints_satisfied``.

    Returns
    -------
    AgreementReport
        Per-path pass/fail + numbers + the overall ``agreed`` verdict
        (>= 2 paths available and every available path passing).
    """
    lamberthub = (
        _run_lamberthub_path(cycler, ephem, mu=mu)
        if ephem is not None
        else LamberthubPathResult(available=False, per_leg=(), max_diff_mps=0.0, passed=False)
    )
    construction_optimiser = _run_construction_optimiser_path(
        cycler, mu=mu, optimiser_vinf=optimiser_vinf_kms
    )
    kepler_reprop = _run_kepler_reprop_path(cycler, mu=mu)

    paths = (lamberthub, construction_optimiser, kepler_reprop)
    n_available = sum(1 for p in paths if p.available)
    n_passed = sum(1 for p in paths if p.available and p.passed)
    # Every available path must pass, and at least two must have run.
    all_available_passed = all(p.passed for p in paths if p.available)
    agreed = n_available >= 2 and all_available_passed

    return AgreementReport(
        lamberthub=lamberthub,
        construction_optimiser=construction_optimiser,
        kepler_reprop=kepler_reprop,
        n_paths_available=n_available,
        n_paths_passed=n_passed,
        agreed=agreed,
    )


__all__ = [
    "KEPLER_REPROP_TOL_KM",
    "VINF_AGREEMENT_TOL_KMS",
    "AgreementReport",
    "ConstructionOptimiserPathResult",
    "KeplerRepropPathResult",
    "LamberthubPathResult",
    "crosscheck_code_paths",
]
