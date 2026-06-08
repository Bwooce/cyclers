"""Batch driver over the #158 continuation: V1->V3 lift for the closable set.

This generalises the single-row #158 continuation driver
(:mod:`cyclerfinder.search.continuation`) into a *batch* tool that, for each input
catalogue row, finds the DE440 basin near the row's sourced launch window, runs the
homotopy ladder out to the true ephemeris (DE440), and returns a structured per-row
result suitable for V3 evidence assessment.

What it can lift (the KEY constraint, from the #158 deep-dive)
-------------------------------------------------------------
The continuation walks UP from a circular-coplanar CLOSURE: it can only promote
rows that ALREADY close circular-coplanar. That is the #137 set — the four Russell
V1 free-return rows (whose single heliocentric ellipse forms a closed,
V_inf-continuous E->M->E cycler on the circular ephemeris) plus the Aldrin row.
Rows lacking a circular-coplanar closure (the multi-arc / Lambert-singular set, and
S1L1) CANNOT be seeded here; they need the closer sweep / Appendix C, not this
driver.

The seed (sourced -> constraint) and the basin search
-----------------------------------------------------
For the Russell rows the seed ``(a, e)`` is derived from the SOURCED aphelion +
outbound transit via :func:`free_return.seed_ae_from_aphelion_transit` (the #106
shared physics, identical to ``scripts/campaign_russell12.py`` and
``tests/search/test_free_return_v1_mechanics.py``). For Aldrin the sourced ``(a, e)``
is taken directly (Rogers 2012 Table 1, catalogue ``aldrin-classic-em-k1``).

The continuation ``t0`` seed is NOT a circular phase scan: the circular model lives
in a θ=0-at-J2000 frame whose absolute epoch does not coincide with DE440's, so a
circular-best basin can sit ~600 d off the real window (verified). Instead the
driver sweeps seed dates across a bounded window around the SOURCED launch and runs
the continuation at each; the basin SELECTED is the converged one whose EMERGED
per-body V_inf best matches the row's INDEPENDENTLY sourced anchor (V_inf is the
row's physical fingerprint — distinct DE440 basins are distinct cyclers, and only
the V_inf-matching basin is the continuation of *this* row), tie-broken toward the
window. This mirrors the Aldrin V3 test's basin resolution (a DE440 basin near, not
exactly at, the sourced window). The seed/anchor are CONSTRAINTS (sourced); the
emerged V_inf and the phase-matched ``t0`` that EMERGE are the evidence.

What emerges (the evidence)
---------------------------
Per row the result records: the per-rung final residuals + winning rung, the
true-ephemeris best-final residual (the "ballistic within tol" closure evidence),
the EMERGED per-body V_inf vs the INDEPENDENTLY sourced anchor, and the
phase-matched ``t0`` vs the sourced launch window (the §14 V3 phase-match half).
The bounded-horizon-TCM half of V3 is built per-row in
``tests/verify/test_continuation_v3_batch.py`` (it needs the maintenance machinery
and is the slow gate).

Pure: depends only on core, search/free_return, search/continuation. No catalogue
WRITE (writebacks are consolidated by the main session); a catalogue READ helper is
provided for the test surface but the physics functions take explicit seeds.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import UTC, datetime, timedelta
from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]

from cyclerfinder.core.constants import DAYS_PER_JULIAN_YEAR, MU_SUN_KM3_S2, SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search import continuation as cont
from cyclerfinder.search.free_return import seed_ae_from_aphelion_transit

DAY_S = SECONDS_PER_DAY
_J2000 = datetime(2000, 1, 1, 12, 0, 0, tzinfo=UTC)

# The launch-window proximity gate (days) — the #158 gate the Aldrin V3 test uses.
WINDOW_GATE_DAYS = 200

# The DE440-basin seed sweep around the sourced launch window. The basin a row
# continues into depends sensitively on the seed epoch (verified: a 20 d shift can
# move the emerged V_inf from 5 to 10 km/s — distinct cyclers), so we sweep seed
# dates across [-SEED_SWEEP_DAYS, +SEED_SWEEP_DAYS] of the window in SEED_STEP_DAYS
# steps and select by V_inf fingerprint. The span stays within the WINDOW_GATE so a
# selected basin can still phase-match; the step is the granularity at which basins
# resolve (finer than the ~140 d basin width seen in the sweeps).
SEED_SWEEP_DAYS = 200
SEED_STEP_DAYS = 10

# Closure tolerance (km/s) — OUR documented floor (free_return default tol_kms);
# Russell's "0 m/s" is below his unprinted SNOPT/post-processing floor.
TOL_KMS = 0.1

# Loose emerged-V_inf corroboration band (km/s) — the #158 golden tolerance: a
# coplanar single-ellipse free-return vs the real 3-D ephemeris pins the family to
# ~1.5 km/s. Basins whose V_inf is further than this from the sourced anchor are a
# DIFFERENT family, not this row's continuation.
VINF_MATCH_KMS = 1.5

DEFAULT_CATALOGUE_PATH = Path(__file__).resolve().parents[3] / "data" / "catalogue.yaml"


# ---------------------------------------------------------------------------
# Inputs: the V1->V3 lift set (circular-coplanar closable rows + Aldrin)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class BatchRowSeed:
    """The free-return seed + sourced anchors for one liftable row.

    Attributes
    ----------
    rid:
        Catalogue row id.
    a_seed_au, e_seed:
        The CONSTRAINED free-return seed (sourced aphelion+transit for the Russell
        rows; sourced (a, e) for Aldrin).
    period_sec:
        The cycler period (catalogue ``period.years``) carried through the
        free-return genome — matches the #137 / #158 test surface (it drives only
        the diagnostic Earth-Earth interval; the closure residual is period-free).
    sourced_launch:
        The SOURCED launch window (catalogue ``priority_date`` / Russell window) —
        the EXPECTED side of the V3 phase-match half, and the centre of the basin
        sweep.
    sourced_vinf:
        The INDEPENDENTLY sourced per-body V_inf anchors (Russell tables) — the
        EXPECTED fingerprint used to select the right DE440 basin and the EXPECTED
        side of the emerged-V_inf corroboration.
    """

    rid: str
    a_seed_au: float
    e_seed: float
    period_sec: float
    sourced_launch: datetime
    sourced_vinf: dict[str, float]


@dataclass(frozen=True)
class BatchRowResult:
    """Per-row continuation outcome with V3-evidence-ready fields.

    Attributes
    ----------
    seed:
        The :class:`BatchRowSeed` driven.
    seed_offset_days:
        The selected seed's offset (days) from the sourced launch window — which
        basin the sweep landed on.
    result:
        The raw :class:`continuation.ContinuationResult` (full per-rung audit) for
        the SELECTED basin.
    best_final_residual_kms:
        The true-ephemeris best-final residual (the "ballistic within tol"
        closure evidence).
    winning_nstep:
        Ladder rung that produced ``best_final``.
    ballistic_within_tol:
        ``best_final`` converged at the ephemeris step (residual < tol).
    emerged_vinf_kms:
        The EMERGED per-body V_inf on the best-final solution.
    vinf_offset_kms:
        ``|emerged - sourced|`` per body (the corroboration metric).
    vinf_matched:
        Every body within :data:`VINF_MATCH_KMS` of its sourced anchor.
    emerged_t0:
        The phase-matched launch datetime of the best-final solution.
    window_offset_days:
        ``|emerged_t0 - sourced_launch|`` in days (the §14 V3 phase-match metric).
    phase_matched:
        ``window_offset_days < WINDOW_GATE_DAYS``.
    selected:
        True if a converged, V_inf-matching basin was found within the sweep; False
        means the sweep found no basin that is this row's continuation (a break
        point — ``result`` then carries the lowest-residual basin tried, for the
        audit trail).
    """

    seed: BatchRowSeed
    seed_offset_days: float
    result: cont.ContinuationResult
    best_final_residual_kms: float
    winning_nstep: int
    ballistic_within_tol: bool
    emerged_vinf_kms: dict[str, float]
    vinf_offset_kms: dict[str, float]
    vinf_matched: bool
    emerged_t0: datetime
    window_offset_days: float
    phase_matched: bool
    selected: bool

    @property
    def reaches_v3_closure(self) -> bool:
        """The V3 phase-match half: converged ballistic, V_inf-matched, in-window."""
        return (
            self.selected and self.ballistic_within_tol and self.vinf_matched and self.phase_matched
        )


def epoch_sec(dt: datetime) -> float:
    """Seconds since the J2000 epoch (the continuation/ephemeris time base)."""
    return (dt - _J2000).total_seconds()


def run_continuation_at_seed(
    seed: BatchRowSeed,
    t0_seed_sec: float,
    *,
    ladder: tuple[int, ...] = cont.LADDER,
    final_ephemeris: Ephemeris,
    bodies: tuple[str, str] = ("E", "M"),
    mu: float = MU_SUN_KM3_S2,
    tol_kms: float = TOL_KMS,
) -> cont.ContinuationResult:
    """Run the #158 continuation ladder once at an explicit ``t0`` seed."""
    return cont.continuation_correct(
        t0_seed_sec,
        seed.a_seed_au,
        seed.e_seed,
        seed.period_sec,
        bodies=bodies,
        mu=mu,
        tol_kms=tol_kms,
        ladder=ladder,
        final_ephemeris=final_ephemeris,
    )


def _package(
    seed: BatchRowSeed,
    seed_offset_days: float,
    result: cont.ContinuationResult,
    *,
    selected: bool,
) -> BatchRowResult:
    bf = result.best_final
    emerged_vinf = dict(bf.vinf_kms)
    vinf_offset = {
        b: abs(emerged_vinf.get(b, float("nan")) - seed.sourced_vinf[b]) for b in seed.sourced_vinf
    }
    vinf_matched = bool(vinf_offset) and all(d <= VINF_MATCH_KMS for d in vinf_offset.values())
    emerged_t0 = _J2000 + timedelta(seconds=bf.t0_sec)
    window_offset_days = abs((emerged_t0 - seed.sourced_launch).days)
    return BatchRowResult(
        seed=seed,
        seed_offset_days=seed_offset_days,
        result=result,
        best_final_residual_kms=bf.max_residual_kms,
        winning_nstep=result.winning_nstep,
        ballistic_within_tol=bf.converged,
        emerged_vinf_kms=emerged_vinf,
        vinf_offset_kms=vinf_offset,
        vinf_matched=vinf_matched,
        emerged_t0=emerged_t0,
        window_offset_days=float(window_offset_days),
        phase_matched=window_offset_days < WINDOW_GATE_DAYS,
        selected=selected,
    )


def run_continuation_for_seed(
    seed: BatchRowSeed,
    *,
    ladder: tuple[int, ...] = (3,),
    final_ephemeris: Ephemeris | None = None,
    bodies: tuple[str, str] = ("E", "M"),
    mu: float = MU_SUN_KM3_S2,
    tol_kms: float = TOL_KMS,
    sweep_days: int = SEED_SWEEP_DAYS,
    step_days: int = SEED_STEP_DAYS,
) -> BatchRowResult:
    """Find the row's DE440 basin near its window and continue to the true ephemeris.

    Sweeps seed dates across ``[-sweep_days, +sweep_days]`` of the sourced launch in
    ``step_days`` steps, runs the continuation ladder at each, and SELECTS the
    converged basin whose emerged per-body V_inf best matches the sourced anchor
    (tie-broken toward the window). The default ``ladder=(3,)`` is the documented
    winning rung (cheap, deterministic, the Aldrin V3 reproduction path); pass the
    full ladder for the golden gate. If no converged V_inf-matching basin exists in
    the sweep the lowest-residual basin tried is returned with ``selected=False``
    (a recorded break point, never a silent cap).

    The seed/anchor are CONSTRAINTS (sourced); the emerged V_inf and phase-matched
    ``t0`` are the evidence (the golden-rule separation from :mod:`free_return`).
    """
    if final_ephemeris is None:
        final_ephemeris = Ephemeris("astropy")

    lw = epoch_sec(seed.sourced_launch)
    offsets = list(range(-sweep_days, sweep_days + 1, step_days))

    best_match: BatchRowResult | None = None
    best_match_key: tuple[float, float] | None = None
    best_any: BatchRowResult | None = None
    for off in offsets:
        result = run_continuation_at_seed(
            seed,
            lw + off * DAY_S,
            ladder=ladder,
            final_ephemeris=final_ephemeris,
            bodies=bodies,
            mu=mu,
            tol_kms=tol_kms,
        )
        pkg = _package(seed, float(off), result, selected=False)
        # Track the lowest-residual basin overall (the break-point fallback).
        if best_any is None or pkg.best_final_residual_kms < best_any.best_final_residual_kms:
            best_any = pkg
        # A selectable basin: converged AND V_inf-matched AND phase-matched.
        if pkg.ballistic_within_tol and pkg.vinf_matched and pkg.phase_matched:
            total_vinf_off = sum(pkg.vinf_offset_kms.values())
            key = (total_vinf_off, abs(float(off)))
            if best_match_key is None or key < best_match_key:
                best_match_key, best_match = key, pkg

    if best_match is not None:
        return _package(
            best_match.seed, best_match.seed_offset_days, best_match.result, selected=True
        )
    assert best_any is not None
    return best_any


# ---------------------------------------------------------------------------
# Catalogue -> seed (the V1->V3 lift set)
# ---------------------------------------------------------------------------


# Aldrin's sourced ellipse is given directly (a, e), not via aphelion+transit. The
# launch window is Russell's Aug-2003 (Fig.5.7). Period is the cycler resonant
# period (catalogue period.years).
_ALDRIN_RID = "aldrin-classic-em-k1-outbound"
_ALDRIN_A_AU = 1.60
_ALDRIN_E = 0.393
_ALDRIN_LAUNCH = datetime(2003, 8, 6, tzinfo=UTC)  # Russell p.176 / Fig.5.7

# The four Russell free-return rows currently at validation_level V1 (the #137 set
# whose single ellipse closes circular-coplanar with V_inf continuity). These are
# the ONLY catalogue rows the continuation can lift besides Aldrin.
RUSSELL_V1_RIDS: tuple[str, ...] = (
    "russell-ch4-5.30gGf3",
    "russell-ch4-9.94Gg3",
    "russell-ch4-5.75ggF3",
    "russell-ch4-9.353Gg2",
)

# The full V1->V3 lift set driven by the batch (Russell V1 rows + Aldrin).
LIFT_SET_RIDS: tuple[str, ...] = (*RUSSELL_V1_RIDS, _ALDRIN_RID)


def _load_rows(catalogue_path: Path | str | None = None) -> list[dict[str, Any]]:
    path = Path(catalogue_path) if catalogue_path is not None else DEFAULT_CATALOGUE_PATH
    rows = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(rows, list):
        raise ValueError(f"catalogue at {path} is not a list of rows")
    return rows


def _parse_launch(row: dict[str, Any]) -> datetime:
    """The sourced launch window from the row's ``priority_date`` (ISO date)."""
    pd = row.get("priority_date")
    if pd is None:
        raise ValueError(f"{row.get('id')}: no priority_date for the sourced launch window")
    dt = datetime.fromisoformat(str(pd))
    return dt if dt.tzinfo is not None else dt.replace(tzinfo=UTC)


def seed_from_catalogue_row(row: dict[str, Any]) -> BatchRowSeed:
    """Build a :class:`BatchRowSeed` from a catalogue row (the lift-set members).

    Russell rows: seed ``(a, e)`` from sourced aphelion+transit (#106 shared
    physics). Aldrin: sourced ``(a, e)`` taken directly. Period is the catalogue
    resonant period; launch is the row ``priority_date`` (Aldrin overrides with the
    sourced Aug-2003 window). Raises ``ValueError`` if the row lacks a
    circular-coplanar seed (no aphelion/transit and not Aldrin) — such rows are NOT
    in the lift set.
    """
    rid = str(row["id"])
    period_sec = float(row["period"]["years"]) * DAYS_PER_JULIAN_YEAR * DAY_S
    sourced_vinf = {
        e["body"]: float(e["vinf_kms"]) for e in (row.get("vinf_kms_at_encounters") or [])
    }

    if rid == _ALDRIN_RID:
        return BatchRowSeed(
            rid=rid,
            a_seed_au=_ALDRIN_A_AU,
            e_seed=_ALDRIN_E,
            period_sec=period_sec,
            sourced_launch=_ALDRIN_LAUNCH,
            sourced_vinf=sourced_vinf,
        )

    oe = row.get("orbit_elements") or {}
    aphelion = oe.get("aphelion_au")
    transit = (row.get("invariants") or {}).get("transit_times_days")
    if aphelion is None or not transit:
        raise ValueError(
            f"{rid}: no circular-coplanar seed (missing sourced aphelion or transit) — "
            f"not in the continuation lift set"
        )
    a_seed, e_seed = seed_ae_from_aphelion_transit(float(aphelion), float(transit[0]))
    return BatchRowSeed(
        rid=rid,
        a_seed_au=a_seed,
        e_seed=e_seed,
        period_sec=period_sec,
        sourced_launch=_parse_launch(row),
        sourced_vinf=sourced_vinf,
    )


def seeds_for_ids(
    rids: tuple[str, ...] = LIFT_SET_RIDS,
    *,
    catalogue_path: Path | str | None = None,
) -> list[BatchRowSeed]:
    """Build seeds for the given catalogue ids (default: the full lift set)."""
    rows = {str(r["id"]): r for r in _load_rows(catalogue_path)}
    seeds: list[BatchRowSeed] = []
    for rid in rids:
        if rid not in rows:
            raise KeyError(f"catalogue has no row {rid!r}")
        seeds.append(seed_from_catalogue_row(rows[rid]))
    return seeds


def run_batch(
    rids: tuple[str, ...] = LIFT_SET_RIDS,
    *,
    catalogue_path: Path | str | None = None,
    ladder: tuple[int, ...] = (3,),
    final_ephemeris: Ephemeris | None = None,
    sweep_days: int = SEED_SWEEP_DAYS,
    step_days: int = SEED_STEP_DAYS,
) -> list[BatchRowResult]:
    """Run the continuation for every id in the lift set; return per-row results.

    A single ``final_ephemeris`` (default a fresh ``Ephemeris('astropy')``) is
    shared across rows for efficiency. The seed for each row is derived from the
    catalogue; rows lacking a circular-coplanar seed raise (they are not liftable).
    """
    if final_ephemeris is None:
        final_ephemeris = Ephemeris("astropy")
    seeds = seeds_for_ids(rids, catalogue_path=catalogue_path)
    return [
        run_continuation_for_seed(
            s,
            ladder=ladder,
            final_ephemeris=final_ephemeris,
            sweep_days=sweep_days,
            step_days=step_days,
        )
        for s in seeds
    ]


__all__ = [
    "LIFT_SET_RIDS",
    "RUSSELL_V1_RIDS",
    "SEED_STEP_DAYS",
    "SEED_SWEEP_DAYS",
    "TOL_KMS",
    "VINF_MATCH_KMS",
    "WINDOW_GATE_DAYS",
    "BatchRowResult",
    "BatchRowSeed",
    "epoch_sec",
    "run_batch",
    "run_continuation_at_seed",
    "run_continuation_for_seed",
    "seed_from_catalogue_row",
    "seeds_for_ids",
]
