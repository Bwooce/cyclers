"""M-ED Phase 4: seeding ladder (plan Phase 4; spec §3)."""

from __future__ import annotations

from cyclerfinder.search.seed_ladder import resolve_seed
from cyclerfinder.search.sequence import Cell


def _s1l1_cell() -> Cell:
    return Cell(
        bodies=("E", "M"),
        sequence=("E", "M", "E", "E"),
        period_k=2,
        per_leg_revs=(0, 0, 0),
        per_leg_branch=("single", "single", "single"),
    )


def _emevve_cell() -> Cell:
    return Cell(
        bodies=("E", "M", "V"),
        sequence=("E", "M", "E", "V", "V", "E"),
        period_k=2,
        per_leg_revs=(0, 0, 0, 0, 0),
        per_leg_branch=("single", "single", "single", "single", "single"),
        period_basis=("E", "M"),
    )


def test_descriptor_rung_used_when_arcs_present() -> None:
    arcs = [
        {"arc_type": "generic", "tof_years": 1.4612, "resonance": None, "raw_descriptor": "g(...)"},
        {"arc_type": "generic", "tof_years": 2.8096, "resonance": None, "raw_descriptor": "G(...)"},
    ]
    plan = resolve_seed(_s1l1_cell(), free_return_arcs=arcs)
    assert plan.source == "descriptor"
    assert len(plan.tof_seed_days) == 2  # two E-E arcs


def test_anchor_rung_used_when_no_descriptor() -> None:
    plan = resolve_seed(
        _emevve_cell(),
        free_return_arcs=None,
        anchor_tofs=(309.0, 259.0),  # Jones EMEVVE transit legs (sourced)
        anchor_vinf={"E": 4.72, "M": 2.50},  # Jones Table 2 (sourced)
    )
    assert plan.source == "anchor"
    assert plan.tof_seed_days[0] == 309.0


def test_coplanar_rung_used_when_no_descriptor_or_anchor() -> None:
    from cyclerfinder.core.ephemeris import Ephemeris

    plan = resolve_seed(
        _s1l1_cell(),
        free_return_arcs=None,
        anchor_tofs=None,
        anchor_vinf=None,
        coplanar_tofs=(154.0, 379.0, 932.0),
    )
    assert plan.source == "coplanar"
    assert plan.tof_seed_days == (154.0, 379.0, 932.0)
    # Ephemeris-derived coplanar warm start is also accepted (no explicit tofs).
    plan2 = resolve_seed(_s1l1_cell(), ephem=Ephemeris(model="circular"))
    assert plan2.source in ("coplanar", "scan")


def test_scan_rung_is_last_resort() -> None:
    plan = resolve_seed(_s1l1_cell())
    assert plan.source == "scan"
    assert len(plan.tof_seed_days) == len(_s1l1_cell().sequence) - 1
    assert all(t > 0 for t in plan.tof_seed_days)


def test_ladder_degrades_in_priority_order() -> None:
    cell = _s1l1_cell()
    arcs = [
        {"arc_type": "generic", "tof_years": 1.4612, "resonance": None, "raw_descriptor": "g(...)"},
        {"arc_type": "generic", "tof_years": 2.8096, "resonance": None, "raw_descriptor": "G(...)"},
    ]
    # Descriptor wins even when anchor + coplanar are also supplied.
    assert (
        resolve_seed(
            cell,
            free_return_arcs=arcs,
            anchor_tofs=(100.0, 200.0, 300.0),
            anchor_vinf={"E": 5.0, "M": 5.0},
            coplanar_tofs=(1.0, 2.0, 3.0),
        ).source
        == "descriptor"
    )
    # Anchor wins over coplanar.
    assert (
        resolve_seed(
            cell,
            anchor_tofs=(100.0, 200.0, 300.0),
            anchor_vinf={"E": 5.0, "M": 5.0},
            coplanar_tofs=(1.0, 2.0, 3.0),
        ).source
        == "anchor"
    )
    # Coplanar wins over scan.
    assert resolve_seed(cell, coplanar_tofs=(1.0, 2.0, 3.0)).source == "coplanar"
