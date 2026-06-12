"""S1L1 real-ephemeris rediscovery — the idealised model could NOT host the
5.65/3.05 km/s anchors (see scripts/characterise_s1l1.py); this gate tests
whether the real ephemeris (Mars eccentricity) can, like the Aldrin cycler.
Anchors are sourced (spec §9); epoch and leg ToFs are computed.

PROVENANCE FLAG (2026-06-04): the 5.65/3.05 km/s anchor pair is
unverified-provenance per catalogue data_gaps[] (s1l1-2syn-em-cpom,
docs/notes/s1l1-target-topology-mining.md). Source-mining found:
Patel 2019 gives Earth flyby v∞ 3.657 km/s, no Mars v∞; McConaghy 2006
abstract gives ≈4.7/5.0 km/s; Sanchez Net 2022 near-ballistic range is
Earth 3.6-5.7 / Mars 5.2-7.3 km/s. Assert values are kept unchanged
pending McConaghy 2006 JSR 43(2):456-465 full-text access."""

import numpy as np
import pytest

from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.optimize import optimise_cell_ephemeris
from cyclerfinder.search.sequence import Cell

VINF_E, VINF_M, TOL = 5.65, 3.05, 0.4  # sourced anchors; real-eph band


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
def test_s1l1_direct_eme_cell_is_off_family() -> None:
    """The 3-encounter DIRECT E-M-E cell does NOT host S1L1 — re-anchored
    2026-06-12 (#214), was an xfail asserting the unverified coplanar 5.65/3.05
    pair.

    S1L1 is CLOSED-CONFIRMED via the CORRECTED topology (E -> g(E-E) -> E flyby
    -> G(E-M-E transit) -> E, App-C #83 seed, per-leg Mars v_inf breathing
    3.2-8.0 km/s), gated in tests/search/test_s1l1_corrected.py. The S/L labels
    are Earth-to-Earth resonant intervals, not an E-M-E direct return leg
    (see memory s1l1-nomenclature), so this construction path is structurally the
    WRONG topology and cannot reach the family. This test now asserts that
    VERIFIED negative — the direct cell lands off-family (constraints
    unsatisfied, emerged v_inf far outside the App-C 3.2-8.0 band) — rather than
    re-asserting the wrong-model 5.65/3.05 anchor (which was also
    unverified-provenance, spec.md §9 only)."""
    cell = Cell(
        bodies=("E", "M"),
        sequence=("E", "M", "E"),
        period_k=2,
        per_leg_revs=(0, 0),
        per_leg_branch=("single", "single"),
    )
    eph = Ephemeris(model="astropy")
    t_syn_em_days = 779.9
    period_days = 2 * t_syn_em_days
    result = optimise_cell_ephemeris(
        cell,
        eph,
        vinf_cap=8.15,
        priority_date_iso="2002-08-05",
        vinf_targets_kms={"E": VINF_E, "M": VINF_M},
        n_starts=5,
        seed=0,
        tof_seed_days=[154.0, period_days - 154.0],
    )
    # Verified structural truth: the direct E-M-E topology is off-family.
    assert not result.constraints_satisfied
    v = _vinf_by_body(result)
    # Emerged Mars v_inf is far above the corrected cycler's 3.2-8.0 band, not
    # near the wrong-model 3.05 anchor — the direct cell is a different basin.
    assert v["M"] > 9.0
