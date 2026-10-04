"""#888 -- the demanded-turn gate wired into the moon-tour validation lane.

Before #888 the V2 / V3 / V4 / V4-strict moon-tour tiers compared only V-infinity
MAGNITUDES at the encounters and restarted every leg from its own Lambert
solution, so six catalogued Uranian rows passed every tier while demanding 1.8 to
28 times the bend their moons can supply. These tests pin the wiring of
:mod:`cyclerfinder.verify.turn_gate` into the lane
(:mod:`cyclerfinder.data.validation.moontour_turn`).

GOLDEN DISCIPLINE. Expected values come from published sources or identities:

* POSITIVE CONTROL (published), at the chain-builder level the lane uses:
  Russell and Strange, "Cycler Trajectories in Planetary Moon Systems", JGCD
  32(1) 2009 (filed in the private paper corpus as
  russell-strange-2009-cycler-trajectories-planetary-moon-systems-JGCD-32-doi-10.2514-1.36610.pdf),
  Tables 3-6: every generic-leg cycler that clears the project floor is
  ACCEPTED and its needed altitude reproduces "Min flyby alt." The two-moon
  lane itself cannot ingest these cyclers (their target moon is massless and
  unphased in the ideal model), so the lane-level acceptance control below is
  the nearest chain the lane CAN ingest, a project computation (#890), and is
  labelled as such.
* WRAP (published): McConaghy, Longuski and Byrnes, AIAA 2002-4420, Table 4
  p.6 and footnote e (filed as
  mcconaghy-longuski-byrnes-2002-analysis-broad-class-earth-mars-cycler-trajectories-AIAA-2002-4420.pdf):
  for an Earth-to-Earth nPr cycler the ONLY flyby is the periodicity wrap, so
  ``search.correct._bend_feasible`` with the wrap must accept 6S7/6S8/6S9 and
  reject 1L1 (Aldrin), and without it accepts everything (vacuously).
* REGRESSION: the six withdrawn Uranian rows (``data/withdrawn/``), fed to the
  lane exactly as the #566 gauntlet and the #330 V2 run fed them (inputs read
  from those runs' stored ``_meta`` and cross-checked against the withdrawn
  rows), are REJECTED by every tier with the turn named as the reason. Only the
  verdict and conservative bounds are asserted, never our earlier digits.
"""

from __future__ import annotations

import importlib
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import yaml  # type: ignore[import-untyped]

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.data.validation.moontour_turn import (
    CycleTurnRecord,
    LegVinf,
    chain_turn_report,
    gate_cycle_record,
    summarize_turn_gate,
)
from cyclerfinder.data.validation.v2_moontour import (
    V2_MOONTOUR_CLOSURE_FLOOR_KMS,
    V2_MOONTOUR_DRIFT_FLOOR_KMS,
    V2MoontourVerdict,
    run_v2_moontour,
)
from cyclerfinder.data.validation.v3_3d import V3_AGREEMENT_FLOOR_KMS, V3Verdict3D, run_v3_3d
from cyclerfinder.data.validation.v4_uranus import (
    V4_AGREEMENT_FLOOR_KMS,
    V4UranusVerdict,
    run_v4_uranus,
)
from cyclerfinder.data.validation.v4_uranus_strict import (
    DEFAULT_LSK_PATH,
    DEFAULT_PCK_PATH,
    DEFAULT_URA_PATH,
    run_v4_uranus_strict,
)
from cyclerfinder.search.correct import _bend_feasible
from cyclerfinder.verify.turn_gate_closures import (
    RS_GENERIC_CYCLERS,
    mcconaghy_npr_cycler,
    russell_strange_generic_cycler,
)

REPO = Path(__file__).resolve().parents[2]
WITHDRAWN_DIR = REPO / "data" / "withdrawn"
GAUNTLET_566 = REPO / "data" / "gauntlet_566_five_representatives.jsonl"
SILVER_330 = REPO / "data" / "silver_327_moontour_v2_verdicts.jsonl"
SILVER_CANDIDATE_ID = "repeated-moon-uranus-00000041"
CANONICAL_EPOCH_UTC = "2000-06-21T00:00:00"  # the #338 / #566 V4-strict anchor epoch

#: Lane candidate id -> withdrawn row id (the #566 REPRESENTATIVES map plus #327).
CANDIDATE_TO_ROW: dict[str, str] = {
    "enum563-line57-titania-oberon-titania": "titania-oberon-1-1-uranian-quasi-cycler-2026",
    "enum563-line18-ariel-oberon-ariel": "ariel-oberon-1-1-uranian-quasi-cycler-2026",
    "enum563-line26-umbriel-titania-umbriel": "umbriel-titania-1-1-uranian-quasi-cycler-2026",
    "enum563-line12-ariel-titania-ariel": "ariel-titania-1-1-uranian-quasi-cycler-2026",
    "enum563-line2-ariel-umbriel-ariel": "ariel-umbriel-1-1-uranian-quasi-cycler-2026",
    SILVER_CANDIDATE_ID: "umbriel-oberon-1-1-uranian-quasi-cycler-2026",
}

#: Floors wide enough that drift and magnitude closure cannot fail, so the turn
#: is the only criterion left to reject on.
GENEROUS: dict[str, Any] = {"drift_floor_kms": 1.0e9, "closure_floor_kms": 1.0e3}

_KERNELS_PRESENT = (
    DEFAULT_LSK_PATH.exists() and DEFAULT_PCK_PATH.exists() and DEFAULT_URA_PATH.exists()
)


@dataclass(frozen=True)
class LaneCase:
    candidate_id: str
    row_id: str
    sequence: tuple[str, ...]
    vinf: tuple[float, ...]
    tofs: tuple[float, ...]
    rel_offset_deg: float
    n_revs: tuple[int, ...]
    phase0_deg: float

    @property
    def args(
        self,
    ) -> tuple[str, tuple[str, ...], tuple[float, ...], tuple[float, ...], float, None]:
        return (
            self.candidate_id,
            self.sequence,
            self.vinf,
            self.tofs,
            self.rel_offset_deg,
            None,
        )

    @property
    def kw(self) -> dict[str, Any]:
        return {"n_cycles": 3, "n_revs": self.n_revs, "phase0_deg": self.phase0_deg}


def _withdrawn_cases() -> list[LaneCase]:
    meta = json.loads(GAUNTLET_566.read_text().splitlines()[0])
    cases = [
        LaneCase(
            candidate_id=c["candidate_id"],
            row_id=CANDIDATE_TO_ROW[c["candidate_id"]],
            sequence=tuple(c["sequence"]),
            vinf=tuple(c["vinf_kms"]),
            tofs=tuple(c["tof_days"]),
            rel_offset_deg=float(c["rel_offset_deg"]),
            n_revs=tuple(c["n_revs"]),
            phase0_deg=float(meta["phase0_deg"]),
        )
        for c in meta["candidates"]
    ]
    silver = json.loads(SILVER_330.read_text().splitlines()[0])["stored_silver"]
    cases.append(
        LaneCase(
            candidate_id=SILVER_CANDIDATE_ID,
            row_id=CANDIDATE_TO_ROW[SILVER_CANDIDATE_ID],
            sequence=tuple(silver["sequence"]),
            vinf=tuple(silver["vinf_per_encounter_kms"]),
            tofs=tuple(silver["tof_days"]),
            rel_offset_deg=float(silver["rel_offset_deg"]),
            n_revs=tuple(silver["n_rev"]),
            phase0_deg=float(silver["phase0_deg"]),
        )
    )
    return cases


CASES = _withdrawn_cases()
CASE_IDS = [c.row_id for c in CASES]


@pytest.fixture(scope="module")
def withdrawn_chain() -> dict[str, tuple[V2MoontourVerdict, V3Verdict3D, V4UranusVerdict]]:
    """Default-floor V2 -> V3 -> V4-scipy for each withdrawn row (about 1 s each)."""
    out: dict[str, tuple[V2MoontourVerdict, V3Verdict3D, V4UranusVerdict]] = {}
    for c in CASES:
        v2 = run_v2_moontour(*c.args, **c.kw)
        v3 = run_v3_3d(*c.args, v2_verdict=v2, **c.kw)
        v4 = run_v4_uranus(*c.args, v3_verdict=v3, **c.kw)
        out[c.row_id] = (v2, v3, v4)
    return out


# ---------------------------------------------------------------------------
# The regression inputs are the withdrawn rows
# ---------------------------------------------------------------------------


def test_six_withdrawn_rows_are_the_lane_inputs() -> None:
    rows = sorted(WITHDRAWN_DIR.glob("*-1-1-uranian-quasi-cycler-2026.yaml"))
    assert {p.stem for p in rows} == set(CANDIDATE_TO_ROW.values())
    assert len(CASES) == 6
    for c in CASES:
        row = yaml.safe_load((WITHDRAWN_DIR / f"{c.row_id}.yaml").read_text())[0]
        assert row["id"] == c.row_id
        assert tuple(row["sequence_canonical"].split("-")) == c.sequence
        assert tuple(float(leg["tof_days"]) for leg in row["legs"]) == c.tofs


# ---------------------------------------------------------------------------
# Regression: the six withdrawn rows are rejected by every tier, on the turn
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("case", CASES, ids=CASE_IDS)
def test_v2_rejects_withdrawn_row_on_the_turn_alone(case: LaneCase) -> None:
    v2 = run_v2_moontour(
        *case.args,
        **case.kw,
        drift_floor_kms=GENEROUS["drift_floor_kms"],
        closure_floor_kms=GENEROUS["closure_floor_kms"],
    )
    # Precondition: with these floors nothing but the turn can fail.
    assert v2.n_cycles_completed == 3
    assert v2.max_drift_kms <= v2.drift_floor_kms
    assert v2.max_closure_residual_kms <= v2.closure_floor_kms
    assert v2.passes_v2 is False
    assert v2.turn_feasible is False
    assert "demanded turn exceeds the available bend" in v2.turn_failure_reason
    assert v2.worst_turn_ratio > 1.5
    for cyc in v2.per_cycle:
        assert cyc.turn_feasible is False
        # One flyby at the other moon and the wrap flyby at the anchor.
        assert [e.label for e in cyc.turn_encounters] == ["e1", "wrap"]
        assert [e.body for e in cyc.turn_encounters] == [case.sequence[1], case.sequence[0]]
        for e in cyc.turn_encounters:
            assert not e.turn_feasible
            assert e.ratio > 1.5
            assert e.required_alt_km < 0.0  # periapsis inside the moon
            assert e.alt_floor_km == SATELLITES[e.body].safe_alt_km


@pytest.mark.parametrize("case", CASES, ids=CASE_IDS)
def test_v3_rejects_withdrawn_row_on_the_turn(case: LaneCase, withdrawn_chain: Any) -> None:
    _v2, v3, _v4 = withdrawn_chain[case.row_id]
    # Precondition: the integrator agreement V3 measures still holds.
    assert v3.n_cycles_propagated == 3
    assert v3.drift_agreement_kms <= V3_AGREEMENT_FLOOR_KMS
    assert v3.passes_v3 is False
    assert v3.turn_feasible is False
    assert "demanded turn exceeds the available bend" in v3.turn_failure_reason


@pytest.mark.parametrize("case", CASES, ids=CASE_IDS)
def test_v4_rejects_withdrawn_row_on_the_turn(case: LaneCase, withdrawn_chain: Any) -> None:
    _v2, _v3, v4 = withdrawn_chain[case.row_id]
    assert v4.n_cycles_propagated == 3
    assert v4.passes_v4 is False
    assert v4.turn_feasible is False
    assert v4.worst_turn_ratio > 1.5
    assert "demanded turn exceeds the available bend" in v4.turn_failure_reason
    for cyc in v4.per_cycle:
        assert cyc.turn_feasible is False
        assert [e.label for e in cyc.turn_encounters] == ["e1", "wrap"]


#: At the canonical epoch the Umbriel-Titania row no longer reaches the turn
#: gate: since #567 every Lambert branch of its first leg is planet-crossing
#: there (the stored #566 PASS predates #567), so the lane already rejects it on
#: that. One day later all three cycles converge and the turn is what rejects it.
STRICT_EPOCH: dict[str, str] = {
    "umbriel-titania-1-1-uranian-quasi-cycler-2026": "2000-06-22T00:00:00",
}


def _strict(case: LaneCase, v3: Any, v4: Any, epoch: str) -> Any:
    return run_v4_uranus_strict(
        case.candidate_id,
        case.sequence,
        case.vinf,
        case.tofs,
        case.rel_offset_deg,
        epoch,
        None,
        v3_verdict=v3,
        v4_scipy_verdict=v4,
        n_cycles=3,
        n_revs=case.n_revs,
    )


@pytest.mark.skipif(not _KERNELS_PRESENT, reason="URA111 SPICE kernel not installed")
@pytest.mark.parametrize("case", CASES, ids=CASE_IDS)
def test_v4_strict_rejects_withdrawn_row_on_the_turn(case: LaneCase, withdrawn_chain: Any) -> None:
    _v2, v3, v4 = withdrawn_chain[case.row_id]
    epoch = STRICT_EPOCH.get(case.row_id, CANONICAL_EPOCH_UTC)
    if epoch != CANONICAL_EPOCH_UTC:
        at_canonical = _strict(case, v3, v4, CANONICAL_EPOCH_UTC)
        assert at_canonical.passes_v4_strict is False
        assert at_canonical.per_cycle[0].failure_mode == "planet_crossing_infeasible"
    v4s = _strict(case, v3, v4, epoch)
    assert v4s.n_cycles_propagated == 3
    assert v4s.passes_v4_strict is False
    assert v4s.turn_feasible is False
    assert v4s.worst_turn_ratio > 1.5
    assert "demanded turn exceeds the available bend" in v4s.turn_failure_reason
    for cyc in v4s.per_cycle:
        assert [e.label for e in cyc.turn_encounters] == ["e1", "wrap"]


@pytest.mark.parametrize(
    "case",
    [c for c in CASES if c.candidate_id != SILVER_CANDIDATE_ID],
    ids=[c.row_id for c in CASES if c.candidate_id != SILVER_CANDIDATE_ID],
)
def test_withdrawn_representatives_used_to_pass_v4(case: LaneCase) -> None:
    """What changed: the frozen #566 output recorded ``passes_v4`` True for these."""
    rows = [json.loads(x) for x in GAUNTLET_566.read_text().splitlines() if x.strip()]
    v4_rows = [
        r
        for r in rows
        if r.get("kind") == "moontour_v4_verdict" and r.get("candidate_id") == case.candidate_id
    ]
    assert v4_rows and all(r["passes_v4"] is True for r in v4_rows)


# ---------------------------------------------------------------------------
# Acceptance: the wiring does not reject everything
# ---------------------------------------------------------------------------


def test_chain_builder_accepts_published_russell_strange_cyclers() -> None:
    """Published flyable cyclers pass the chain builder the lane uses.

    Each Russell-Strange rebuild is re-expressed as per-leg V-infinity vectors
    (leg ``i`` departs with encounter ``i-1``'s outgoing vector and arrives with
    encounter ``i``'s incoming one), wrapped and gated at the PROJECT floors with
    the registry constants. Rows that the project floor excludes (TitEnc#37 and
    #314: 1377 and 1218 km against the 1500 km Titan floor, #888 note section 3)
    are not used. Tolerance on the needed altitude: the registry radii and GMs
    differ from the paper's Table 2 by up to 3 km and 0.01 %, which moves the
    needed altitude by a few km; GanIo#403 keeps the 25 km of the gate's own test
    (its two branches are ambiguous at the printed precision).
    """
    accepted = 0
    for row in RS_GENERIC_CYCLERS:
        if row.rs_id in {"TitEnc#37", "TitEnc#314"}:
            continue
        rb = russell_strange_generic_cycler(row, alt_floor_km=0.0)
        k = len(rb.encounters)
        legs = [
            LegVinf(depart=rb.encounters[i - 1].vinf_out, arrive=rb.encounters[i].vinf_in)
            for i in range(k)
        ]
        seq = (row.flyby_body,) * (k + 1)
        rep = chain_turn_report(seq, legs, legs[0].depart)
        assert len(rep.encounters) == k
        assert rep.encounters[-1].label == "wrap"
        assert rep.turn_feasible, row.rs_id
        tol = 25.0 if row.rs_id == "GanIo#403" else 10.0
        assert rep.min_required_alt_km == pytest.approx(row.min_flyby_alt_km, abs=tol)
        accepted += 1
    assert accepted == 8


def _titania_oberon_890() -> LaneCase:
    """The one turn-feasible symmetric closure of the #888 extended enumeration
    (#890): Titania-Oberon-Titania, legs of 5 T_syn / 2, 5 revolutions each,
    rel_offset 0. A PROJECT computation, not a published cycler; it is the
    nearest turn-feasible chain the two-moon lane can ingest."""
    mu = PRIMARIES["Uranus"]

    def period_days(m: str) -> float:
        return 2.0 * math.pi * math.sqrt(SATELLITES[m].sma_km ** 3 / mu) / 86400.0

    t_syn = 1.0 / (1.0 / period_days("Titania") - 1.0 / period_days("Oberon"))
    tof = 5.0 * t_syn / 2.0
    return LaneCase(
        candidate_id="890-titania-oberon-5rev",
        row_id="890",
        sequence=("Titania", "Oberon", "Titania"),
        vinf=(0.27, 0.27, 0.27),  # audit only; the lane re-solves from geometry
        tofs=(tof, tof),
        rel_offset_deg=0.0,
        n_revs=(5, 5),
        phase0_deg=0.0,
    )


def test_lane_accepts_a_turn_feasible_chain() -> None:
    """Same generous drift floor as the regression above (the V2 drift is an
    inertial position offset of the closing encounter, about 7e5 km here, as
    for every symmetric closure whose cycle is not also a multiple of the
    anchor's own period); the default magnitude-closure floor is kept. Under the
    same floors the six withdrawn rows fail and this chain passes."""
    case = _titania_oberon_890()
    v2 = run_v2_moontour(*case.args, **case.kw, drift_floor_kms=GENEROUS["drift_floor_kms"])
    assert v2.turn_feasible is True
    assert v2.turn_failure_reason == ""
    assert 0.0 < v2.worst_turn_ratio < 1.0
    assert v2.max_closure_residual_kms <= V2_MOONTOUR_CLOSURE_FLOOR_KMS
    assert v2.passes_v2 is True
    # At the default 50,000 km drift floor it fails on drift, not on the turn.
    v2_default = run_v2_moontour(*case.args, **case.kw)
    assert v2_default.turn_feasible is True
    assert v2_default.max_drift_kms > V2_MOONTOUR_DRIFT_FLOOR_KMS
    assert v2_default.passes_v2 is False
    v3 = run_v3_3d(*case.args, v2_verdict=v2, **case.kw)
    assert v3.turn_feasible is True
    assert v3.passes_v3 is (v3.drift_agreement_kms <= V3_AGREEMENT_FLOOR_KMS)
    v4 = run_v4_uranus(*case.args, v3_verdict=v3, **case.kw)
    assert v4.turn_feasible is True
    # A feasible turn changes nothing: the verdict is the old conjunction.
    old = (
        v4.n_cycles_propagated == 3
        and v4.drift_agreement_kms <= V4_AGREEMENT_FLOOR_KMS
        and v4.bounded_drift_survives
    )
    assert v4.passes_v4 is old


def test_lane_rejects_a_chain_whose_only_infeasible_flyby_is_the_wrap() -> None:
    """Titania-Umbriel-Titania, legs of 3 T_syn (commensurate, so the phasing
    repeats), one revolution each, rel_offset 0: the Umbriel flyby is feasible,
    the closing Titania flyby is not. With the drift floor opened, the pre-#888
    lane passed it (magnitude closure within its default floor); a check of the
    intermediate encounters alone would still pass it. A project computation, found by scanning the
    lane's own construction. Drift floor generous as above (inertial offset of
    the closing encounter); magnitude-closure floor at its default."""
    mu = PRIMARIES["Uranus"]

    def period_days(m: str) -> float:
        return 2.0 * math.pi * math.sqrt(SATELLITES[m].sma_km ** 3 / mu) / 86400.0

    t_syn = 1.0 / abs(1.0 / period_days("Titania") - 1.0 / period_days("Umbriel"))
    tof = 3.0 * t_syn
    v2 = run_v2_moontour(
        "wrap-only",
        ("Titania", "Umbriel", "Titania"),
        (1.0, 1.0, 1.0),
        (tof, tof),
        0.0,
        None,
        n_cycles=3,
        n_revs=(1, 1),
        phase0_deg=0.0,
        drift_floor_kms=GENEROUS["drift_floor_kms"],
    )
    # Precondition: the pre-#888 criteria all pass.
    assert v2.n_cycles_completed == 3
    assert v2.max_drift_kms <= v2.drift_floor_kms
    assert v2.max_closure_residual_kms <= V2_MOONTOUR_CLOSURE_FLOOR_KMS
    for cyc in v2.per_cycle:
        mid, wrap = cyc.turn_encounters
        assert (mid.body, mid.label) == ("Umbriel", "e1")
        assert mid.turn_feasible and mid.ratio < 0.6
        assert (wrap.body, wrap.label) == ("Titania", "wrap")
        assert not wrap.turn_feasible and wrap.ratio > 5.0
    assert v2.passes_v2 is False
    assert "Titania (wrap)" in v2.turn_failure_reason


def test_chain_builder_gates_the_wrap() -> None:
    """Identity: no turn at the intermediate flyby, 90 degrees at the wrap."""
    v = np.array([0.9, 0.0, 0.0])
    rot90 = np.array([0.0, 0.9, 0.0])
    legs = [LegVinf(depart=v, arrive=v), LegVinf(depart=v, arrive=v)]
    rep = chain_turn_report(("Titania", "Oberon", "Titania"), legs, rot90)
    assert [e.label for e in rep.encounters] == ["e1", "wrap"]
    assert rep.encounters[0].demanded_turn_deg == pytest.approx(0.0, abs=1e-9)
    assert rep.encounters[1].demanded_turn_deg == pytest.approx(90.0, abs=1e-9)
    assert not rep.turn_feasible
    rec = CycleTurnRecord(legs=legs, wrap_depart=rot90)
    gate_cycle_record(rec, ("Titania", "Oberon", "Titania"))
    summary = summarize_turn_gate([rec])
    assert not summary.turn_feasible
    assert "Titania (wrap)" in summary.turn_failure_reason


def test_unflyable_or_degenerate_cycle_is_not_feasible() -> None:
    seq = ("Titania", "Oberon", "Titania")
    v = np.array([0.9, 0.0, 0.0])
    missing_wrap = CycleTurnRecord(legs=[LegVinf(v, v), LegVinf(v, v)])
    gate_cycle_record(missing_wrap, seq)
    zero = CycleTurnRecord(legs=[LegVinf(v, np.zeros(3)), LegVinf(v, v)], wrap_depart=v)
    gate_cycle_record(zero, seq)  # zero V-infinity: recorded, not raised
    for rec in (missing_wrap, zero):
        assert rec.report is None and rec.error
        assert not summarize_turn_gate([rec]).turn_feasible
    assert not summarize_turn_gate([]).turn_feasible


# ---------------------------------------------------------------------------
# Heliocentric corrector: the wrap in _bend_feasible (published control)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("n", "p", "r", "ballistic"),
    [(1, "L", 1, False), (6, "S", 7, True), (6, "S", 8, True), (6, "S", 9, True)],
)
def test_bend_feasible_checks_the_wrap(n: int, p: str, r: int, ballistic: bool) -> None:
    """McConaghy et al. 2002 Table 4 footnote e: 6S7-6S9 ballistic, 1L1 not."""
    c = mcconaghy_npr_cycler(n, p, r)  # type: ignore[arg-type]
    assert c is not None
    psi = 2.0 * math.pi * n * (15.0 / 7.0)
    cs, sn = math.cos(-psi), math.sin(-psi)
    unrot = np.array([[cs, -sn, 0.0], [sn, cs, 0.0], [0.0, 0.0, 1.0]])
    nodes = {"b0_out": unrot @ c.encounter.vinf_out, "b1_in": c.encounter.vinf_in}
    seq = ("E", "E")
    # Intermediate encounters only: an E-E chain has none, so this is vacuous.
    assert _bend_feasible(nodes, seq, None) is True
    assert _bend_feasible(nodes, seq, None, wrap_rotation_rad=psi) is ballistic


# ---------------------------------------------------------------------------
# Discovery: the #558 candidate gate (and its generic callers)
# ---------------------------------------------------------------------------


def _scan_558() -> Any:
    """The #558 script, imported the way its sibling scripts import it (bare
    name from ``scripts/``, which is also mypy's ``mypy_path``)."""
    scripts_dir = str(REPO / "scripts")
    if scripts_dir not in sys.path:
        sys.path.insert(0, scripts_dir)
    return importlib.import_module("scan_558_uranus_all_pairs_offset_sweep")


def _record_558(case: LaneCase) -> dict[str, Any]:
    s558 = _scan_558()
    mu = PRIMARIES["Uranus"]
    anchor, flyby = case.sequence[0], case.sequence[1]
    p = [
        2.0 * math.pi * math.sqrt(SATELLITES[m].sma_km ** 3 / mu) / 86400.0 for m in (anchor, flyby)
    ]
    rec = s558.residual_at_point(
        anchor,
        flyby,
        rel_offset_deg=case.rel_offset_deg,
        tof_scale=case.tofs[0] / math.sqrt(p[0] * p[1]),
        n_rev=(case.n_revs[0], case.n_revs[1]),
        phase0_deg=case.phase0_deg,
    )
    assert rec is not None
    return dict(rec)


def test_scan_558_gate_candidate_rejects_withdrawn_closure_on_the_turn() -> None:
    s558 = _scan_558()
    case = next(c for c in CASES if c.row_id.startswith("titania-oberon"))
    g = s558.gate_candidate(case.sequence[0], case.sequence[1], _record_558(case))
    assert g["physical_gate_passed"] is True  # the #324 capacity gate let it through
    assert g["turn_gate_passed"] is False
    assert g["all_gates_passed"] is False
    assert g["turn_gate"]["worst_ratio"] > 1.5


def test_scan_558_gate_candidate_passes_turn_feasible_closure() -> None:
    s558 = _scan_558()
    case = _titania_oberon_890()
    g = s558.gate_candidate(case.sequence[0], case.sequence[1], _record_558(case))
    assert g["turn_gate_passed"] is True
    assert g["turn_gate"]["worst_ratio"] < 1.0
