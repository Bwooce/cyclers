"""#135/#820 like-for-like seed-at-truth probe: pins the decisive diagnostic.

HISTORY. The #135 diagnosis ("the sourced geometry is NOT a residual-zero point
of the descriptor->corrector genome") was pinned here 2026-06-06 — and #820
(2026-08-11) showed it was an artifact of a MIS-POSED genome: ``build_genome``
assumed ``arcs[0]`` was the Mars leg and used the sourced t_in as the M->E leg
ToF, whereas #794's primary-source semantics put the Mars transit on the
DESIGNATED (uppercase) arc and require the M->E leg to be that arc's
beyond-Mars REMAINDER (so the seeds can tile the period). Under the corrected
posing the corrector, seeded at truth, converges to a genuine residual-zero
closure within ~1 day of the sourced geometry whose EMERGED V-infinity matches
the row's independently sourced anchors. This file now pins THAT finding on
``russell-ch4-4.991gG2`` (Russell 2004 Table 4.9 row 1, anchors E 4.99 /
M 5.10 km/s). GOLDEN DISCIPLINE: the expected side of the V-infinity check is
the SOURCED anchor; V-infinity is never imposed. See
``docs/notes/2026-08-11-820-russell12-designated-arc-reposing.md`` (and
``docs/notes/2026-06-06-russell12-likeforlike.md`` for the superseded #135
diagnosis).
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
import yaml  # type: ignore[import-untyped]

REPO_ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = REPO_ROOT / "scripts" / "campaign_russell12.py"


def _load_campaign() -> ModuleType:
    spec = importlib.util.spec_from_file_location("campaign_russell12", CAMPAIGN)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _row(rid: str) -> dict[str, Any]:
    rows = yaml.safe_load((REPO_ROOT / "data" / "catalogue.yaml").read_text())
    return next(r for r in rows if r["id"] == rid)


@pytest.mark.slow
def test_reposed_truth_seed_converges_to_sourced_anchor() -> None:
    """#820: seeded at the row's OWN sourced geometry under the designated-arc
    posing, the corrector converges to a residual-zero closure a fraction of a
    day from truth, and the EMERGED per-body V-infinity matches the row's
    independently sourced Russell anchors (E 4.99, M 5.10 km/s)."""
    mod = _load_campaign()
    row = _row("russell-ch4-4.991gG2")
    probe = mod.probe_at_truth(row, phase_epochs=64, model="circular")

    assert probe["solved_converged"]
    assert probe["solved_max_residual_kms"] <= mod.CORRECTOR_TOL_KMS
    assert probe["tof_drift_days"] <= mod.TOL_TRANSIT_DAYS
    # Sequence is E-M-E-E: encounter 0 is Earth, encounter 1 is Mars. The
    # sourced anchors are the EXPECTED side (golden discipline); the achieved
    # values EMERGE from the converged Lambert chain (observed 5.008 / 5.107).
    vinf = probe["solved_vinf_per_encounter_kms"]
    assert abs(vinf[0] - 4.99) <= mod.TOL_VINF_KMS
    assert abs(vinf[1] - 5.10) <= mod.TOL_VINF_KMS
    # Residual AT the rounded seed itself sits just above the 0.1 km/s floor —
    # period.years is printed to 2 decimals (4.27 vs the exact 4.2708), which
    # alone injects ~1.3 d of slack-leg error. Pin the regime: sub-km/s, an
    # order of magnitude below the pre-#820 mis-posed genome's ~3.4 km/s.
    assert probe["best_phase_truth_residual_kms"] < 1.0


@pytest.mark.slow
def test_genome_designated_arc_split_and_loop_seeds() -> None:
    """#820: the genome splits the DESIGNATED (uppercase) arc — ``arcs[1]``
    (G) for this gG row, NOT ``arcs[0]`` — at the Mars encounter, and seeds
    the E-E loop from the #794-written-back segment."""
    mod = _load_campaign()
    row = _row("russell-ch4-4.991gG2")
    genome = mod.build_genome(row)
    assert genome["designated_index"] == 1
    assert genome["designated_raw"].startswith("G(")
    # Leg 0 = sourced t_out; leg 1 = designated remainder (2.8096 yr - 150 d);
    # the loop (#794 segment, 533.7 d) is eliminated as the slack leg.
    truth_free = mod._truth_seed(genome)
    assert truth_free[0] == 150.0
    assert truth_free[1] == pytest.approx(2.8096 * 365.25 - 150.0, abs=0.01)
    assert genome["slack_leg"] == 2
    assert genome["all_seeds"][2] == 533.7
    # The seeds now tile the sourced period to ~1 d (pre-#820 they could not:
    # the beyond-Mars designated-leg time had no leg to live in).
    period_days = genome["period_sec"] / 86400.0
    assert abs(sum(genome["all_seeds"]) - period_days) < 2.0


# ---------------------------------------------------------------------------
# #137 ACCEPTANCE GATE — the reworked free-return (radial-crossing) genome.
#
# The free-Lambert genome above makes truth NOT a closure point (WALKED-AWAY).
# The free-return genome (cyclerfinder.search.free_return) expresses the Mars
# transfer as a free-return arc on a single heliocentric ellipse; the per-body
# V_inf and leg ToFs EMERGE from the ellipse shape. Seeding (a, e) at the SOURCED
# S1L1 ellipse (a=1.30 AU, e=0.257; Rogers 2012, the SAME physical cycler per
# docs/notes/multi-arc-classification.md §12) at the best phase yields residual
# ~ 0 -- the sourced geometry IS now representable.
#
# GOLDEN DISCIPLINE: (a, e) is the SOURCED input (constraint); the asserted
# EVIDENCE is the EMERGED V_inf compared to the INDEPENDENTLY sourced V_inf
# anchors (Russell 4.99/5.10, McConaghy 4.7/5.0). V_inf is never imposed.
# ---------------------------------------------------------------------------

# Sourced S1L1 heliocentric ellipse (Rogers 2012 Table 1) — the constraint side.
_S1L1_A_AU = 1.30
_S1L1_E = 0.257


def _best_phase_t0(corr_mod: ModuleType, period_sec: float) -> float:
    """Scan t0 over one period; return the phase minimising the free-return
    residual at the SOURCED (a, e). Mirrors the probe's best-phase selection."""
    import numpy as np

    from cyclerfinder.core.ephemeris import Ephemeris

    ephem = Ephemeris("circular")
    best_t0, best_res = 0.0, float("inf")
    for frac in np.linspace(0.0, 1.0, 360, endpoint=False):
        t0 = float(frac) * period_sec
        res = corr_mod._residuals(
            np.array([_S1L1_A_AU, _S1L1_E, t0]),
            period_days=period_sec / 86400.0,
            ephem=ephem,
            bodies=("E", "M"),
            mu=132712440018.0,
        )
        m = max(abs(r) for r in res)
        if m < best_res:
            best_res, best_t0 = m, t0
    return best_t0


@pytest.mark.slow
@pytest.mark.parametrize(
    ("rid", "src_vinf_e", "src_vinf_m"),
    [
        ("mcconaghy-2006-em-k2", 4.7, 5.0),
        ("russell-ch4-4.991gG2", 4.99, 5.10),
    ],
)
def test_free_return_genome_makes_truth_a_closure_point(
    rid: str, src_vinf_e: float, src_vinf_m: float
) -> None:
    """#137 acceptance gate: the free-return genome reaches residual ~ 0 at the
    SOURCED geometry on the symmetric k=2 rows, and the EMERGED V_inf matches the
    INDEPENDENTLY sourced anchor (within the 0.5 km/s campaign tolerance)."""
    import cyclerfinder.search.free_return as fr
    from cyclerfinder.core.ephemeris import Ephemeris

    row = _row(rid)
    period_sec = float(row["period"]["years"]) * 365.25 * 86400.0
    t0 = _best_phase_t0(fr, period_sec)

    result = fr.free_return_correct(
        t0_seed_sec=t0,
        a_seed_au=_S1L1_A_AU,
        e_seed=_S1L1_E,
        period_sec=period_sec,
        ephem=Ephemeris("circular"),
        tol_kms=0.1,
    )

    # The decisive #137 gate: truth IS a residual-zero point now.
    assert result.converged
    assert result.max_residual_kms < 0.1
    # EVIDENCE (derived, not imposed): emerged V_inf matches the sourced anchor.
    assert abs(result.vinf_kms["E"] - src_vinf_e) <= 0.5
    assert abs(result.vinf_kms["M"] - src_vinf_m) <= 0.5
    # The emerged transfer ToF lands on the sourced ~150-153 d transit.
    assert abs(result.transfer_tof_days - 153.0) <= 5.0
