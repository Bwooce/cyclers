"""S1L1 idealised multi-rev rediscovery — does a 4-encounter E-M-E-E cell with a
multi-rev Earth-to-Earth resonant interval host the published V_inf anchors in
the circular-coplanar model once the STAGE 1 multi-rev plumbing is wired?

The S1/L1 labels are consecutive *Earth-to-Earth* resonant intervals (see
[[s1l1-nomenclature]]); the encounter chain is the outbound E->M plus the
resonant Earth loop, modelled here as the 4-encounter sequence E-M-E-E with a
multi-rev final (Earth-to-Earth) leg.

PROVENANCE
----------
SOURCED EXPECTED anchors: V_inf_Earth = 5.65 km/s, V_inf_Mars = 3.05 km/s
(spec.md §9 / catalogue ``s1l1-2syn-em-cpom`` source_quotes). Tolerance is an
honest +/-0.4 km/s engineering band. The launch epoch, leg ToFs and revolution
counts are COMPUTED.

PROVENANCE FLAG (2026-06-04): the 5.65/3.05 anchor pair is unverified-provenance
per catalogue data_gaps[] (see docs/notes/s1l1-target-topology-mining.md).
Source-mining found: Patel 2019 gives Earth flyby v∞ 3.657 km/s, no Mars v∞;
McConaghy 2006 abstract gives ≈4.7/5.0 km/s; Sanchez Net 2022 near-ballistic
range is Earth 3.6-5.7 / Mars 5.2-7.3 km/s. The 5.65/3.05 pair traces only to
spec.md §9 and is unconfirmed in any mined primary source. Assert values are
kept unchanged pending McConaghy 2006 JSR 43(2):456-465 full-text access.

STATUS: ``@pytest.mark.xfail(strict=False)`` — the topology characterisation
script (``scripts/characterise_s1l1.py``) found NO circular-coplanar topology
(neither E-M-E nor E-M-E-E, over per-leg revs in {0,1,2} and branches
{single,low,high}) that reproduces the 5.65/3.05 anchors. As with the Aldrin
cycler, S1L1's published anchors appear not to be hostable in the idealised
circular-coplanar model — the real-ephemeris (Mars eccentricity) path is
needed. STAGE 1 supplies the multi-rev plumbing; the topology remains open, so
this gate stays xfail until a winning per_leg_revs is pinned (then flips to a
strict passing assert)."""

import numpy as np
import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.optimize import optimise_cell_ephemeris
from cyclerfinder.search.sequence import Cell

VINF_E, VINF_M, TOL = 5.65, 3.05, 0.4  # SOURCED anchors (spec §9); computed band


def _vinf_by_body(result: object) -> dict[str, float]:
    out: dict[str, float] = {}
    for enc in result.best_cycler.encounters:  # type: ignore[attr-defined]
        m = max(
            float(np.linalg.norm(enc.vinf_in)),
            float(np.linalg.norm(enc.vinf_out)),
        )
        out[enc.body] = max(out.get(enc.body, 0.0), m)
    return out


@pytest.mark.slow
@pytest.mark.xfail(
    strict=False,
    reason=(
        "TOPOLOGY OPEN: scripts/characterise_s1l1.py found no circular-coplanar "
        "topology (E-M-E or E-M-E-E, per-leg revs in {0,1,2}, branches "
        "{single,low,high}) reproducing the 5.65/3.05 anchors. STAGE 1 provides "
        "the multi-rev plumbing (per_leg_revs/per_leg_branch now thread through "
        "optimise_cell_ephemeris -> optimise_maintenance_dv -> _build_chain), but "
        "the idealised model cannot host the anchors any more than it can host "
        "Aldrin's. Flips to a strict passing assert once a winning per_leg_revs / "
        "topology is pinned (likely on the real ephemeris). "
        "NOTE (2026-06-04): the 5.65/3.05 anchor pair is unverified-provenance "
        "(catalogue data_gap vinf_kms_at_encounters, s1l1-2syn-em-cpom): traces "
        "only to spec.md §9; unconfirmed in Patel 2019 / McConaghy 2006 / "
        "Sanchez Net 2022 — see docs/notes/s1l1-target-topology-mining.md. "
        "Assert values kept unchanged pending McConaghy 2006 full-text access."
    ),
)
def test_s1l1_4enc_multirev_closes_to_published_vinf_anchors() -> None:
    # E-M-E-E: outbound E->M, then the S1/L1 Earth-to-Earth resonant interval as
    # a multi-rev final leg. period_k=2 (2-synodic cycler).
    cell = Cell(
        bodies=("E", "M"),
        sequence=("E", "M", "E", "E"),
        period_k=2,
        per_leg_revs=(0, 0, 1),
        per_leg_branch=("single", "single", "low"),
    )
    eph = Ephemeris("circular")
    t_syn_em_days = 779.9
    period_days = 2 * t_syn_em_days
    # Asymmetric seed: ~154 d sourced outbound, a direct E->E hop, then the
    # multi-rev resonant interval absorbing the remainder.
    result = optimise_cell_ephemeris(
        cell,
        eph,
        vinf_cap=8.15,
        priority_date_iso="2002-08-05",
        vinf_targets_kms={"E": VINF_E, "M": VINF_M},
        n_starts=5,
        seed=0,
        tof_seed_days=[154.0, 380.0, period_days - 154.0 - 380.0],
    )
    assert result.constraints_satisfied
    v = _vinf_by_body(result)
    assert abs(v["E"] - VINF_E) < TOL
    assert abs(v["M"] - VINF_M) < TOL
