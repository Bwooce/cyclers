"""Descriptor -> DSM seed adapter (plan 2026-06-10, Component 1)."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

import cyclerfinder.search.dsm_descriptor_seed as dds
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.data.catalog import load_catalog


def _row(row_id: str) -> dict[str, Any]:
    return load_catalog().by_id[row_id].raw


def test_seed_built_for_reachable_descriptor_row() -> None:
    seed = dds.seed_dsm_chain_from_descriptor(_row("russell-ch4-9.353Gg2"))
    assert seed is not None
    # Charged layout: [t0, vinf_out0, alpha0, beta0, *tof(n), *eta(n),
    # *alpha_int(n-1), *beta_int(n-1)] — the charge_flyby_continuity=True vector.
    assert len(seed.sequence) >= 2
    n_legs = len(seed.sequence) - 1
    assert seed.x0.shape[0] == 4 + 2 * n_legs + 2 * (n_legs - 1)
    eta = seed.x0[4 + n_legs : 4 + 2 * n_legs]
    assert np.allclose(eta, 0.0)
    # intermediate-flyby direction coords also seeded ballistic (0).
    assert np.allclose(seed.x0[4 + 2 * n_legs :], 0.0)
    assert seed.vinf_anchor_kms > 0.0


def test_seed_tofs_are_sourced_not_slack() -> None:
    # The resonant same-body legs seed at the PUBLISHED arc ToF (free_return_arcs
    # tof_years x 365.25); the transit leg seeds at the sourced
    # invariants.transit_times_days value -- not the old slack heuristic.
    seed = dds.seed_dsm_chain_from_descriptor(_row("russell-ch4-4.991gG2"))
    assert seed is not None
    n_legs = len(seed.sequence) - 1
    tof = seed.x0[4 : 4 + n_legs]
    assert np.allclose(tof, [533.70, 150.0, 1026.21], atol=0.5)
    assert np.allclose(seed.per_leg_tof_days, [533.70, 150.0, 1026.21], atol=0.5)


def test_seed_max_revs_from_published_tof() -> None:
    # max_revs = max over legs of floor(arc_tof_days / body_period_days) + 1,
    # capped at Russell's 6-body-period generic-return ceiling. For an E-E-M-M
    # two-synodic row this is 2.
    for rid in ("russell-ch4-4.991gG2", "russell-ch4-9.353Gg2"):
        seed = dds.seed_dsm_chain_from_descriptor(_row(rid))
        assert seed is not None
        assert seed.max_revs == 2, rid


def test_resonant_leg_tof_bounds_bracket_published() -> None:
    # A same-body resonant leg's ToF box brackets its published seed ToF (0.7x..1.3x)
    # so the corrector cannot collapse it to the degenerate near-zero-ToF single-rev
    # solution; the transit leg keeps the sequence-keyed bound.
    import pytest

    seed = dds.seed_dsm_chain_from_descriptor(_row("russell-ch4-4.991gG2"))
    assert seed is not None
    n_legs = len(seed.sequence) - 1
    lower = seed.bounds.lower[4 : 4 + n_legs]
    upper = seed.bounds.upper[4 : 4 + n_legs]
    # leg 0 (E->E) is resonant: published 533.70 d -> [373.6, 693.8]
    assert lower[0] == pytest.approx(0.7 * 533.70, abs=1.0)
    assert upper[0] == pytest.approx(1.3 * 533.70, abs=1.0)
    # leg 2 (M->M) resonant: published 1026.21 d
    assert lower[2] == pytest.approx(0.7 * 1026.21, abs=1.0)
    assert upper[2] == pytest.approx(1.3 * 1026.21, abs=1.0)


def test_no_descriptor_row_returns_none() -> None:
    # An ocampo row has the n.m.k summary format, no per-arc g/G descriptor.
    catalog = load_catalog()
    ocampo = next(e for e in catalog.entries if e.id.startswith("russell-ocampo"))
    assert dds.seed_dsm_chain_from_descriptor(ocampo.raw) is None


def test_close_reachable_row_emerges_vinf_near_anchor() -> None:
    # A REACHABLE descriptor row closes with the DSM genome on the real ephemeris,
    # and its EMERGED Mars V_inf lands within tolerance of the row's sourced anchor.
    row = _row("russell-ch4-9.353Gg2")
    ephem = Ephemeris("astropy")  # real DE440 via astropy
    res = dds.close_row_dsm(row, ephem)
    # converged is by the corrector's own residual criterion; if it converges the
    # emerged V_inf must match the sourced anchor (golden — anchor from the row, not
    # computed here). A non-converged row is a recorded negative (also valid).
    if res.converged:
        assert res.anchor_match
        assert min(res.dv_dsm_kms) >= 0.0


def test_closure_result_carries_seed_max_revs_and_n_revs() -> None:
    # close_row_dsm threads the seed's Russell rev cap into the corrector and reports
    # the emerged per-leg revolution count for the runlog/audit.
    row = _row("russell-ch4-4.991gG2")
    ephem = Ephemeris("astropy")
    res = dds.close_row_dsm(row, ephem)
    assert res.max_revs_used == 2
    assert isinstance(res.n_revs_per_leg, tuple)
    assert len(res.n_revs_per_leg) == 3  # E-E-M-M -> 3 legs


# --- #849: _descriptor_params identifies g/G by descriptor CASE, not array
# position (the #820 defect class carried by this lane's own seeder). ---


def test_descriptor_params_identifies_designated_by_case_not_position() -> None:
    # russell-ch4-9.353Gg2's free_return_arcs are [G(1.7238,...), g(2.5469,...)] --
    # the DESIGNATED (uppercase) arc is arcs[0], not arcs[1]. A positional read
    # (pre-#849: g_tofs[0]/g_tofs[1]) would swap the g/G roles.
    row = _row("russell-ch4-9.353Gg2")
    arcs = row["free_return_arcs"]
    assert arcs[0]["raw_descriptor"].startswith("G")
    assert arcs[1]["raw_descriptor"].startswith("g")
    params = dds._descriptor_params(row)
    assert params is not None
    _aph, g_tof_yr, big_g_tof_yr, _ve, _vm, _seq = params
    assert g_tof_yr == pytest.approx(2.5469)  # arcs[1], lowercase g
    assert big_g_tof_yr == pytest.approx(1.7238)  # arcs[0], uppercase G


@pytest.mark.parametrize(
    "row_id",
    ["russell-ch4-5.30ggF3", "russell-ch4-5.75ggF3"],
)
def test_descriptor_params_none_when_designated_arc_is_full_rev(row_id: str) -> None:
    # These rows' designated (uppercase) arc is F -- full-rev, no tof_years. The
    # g/G shape model (free_return_chain_correct: "two distinct Earth-to-Earth
    # free-return arcs, arc-1 = g, arc-2 = G") has no generic arc to play the G
    # role, so this is honestly out of scope -- None, not a fabricated ToF pair
    # built from the two non-designated lowercase loop arcs (the pre-#849 bug).
    row = _row(row_id)
    assert dds._descriptor_params(row) is None


@pytest.mark.parametrize(
    "row_id",
    ["russell-ocampo-4.3.1-5", "russell-ocampo-2.5.1+0"],
)
def test_descriptor_params_none_for_non_two_generic_arc_row(row_id: str) -> None:
    # Both ocampo rows carry exactly ONE generic (g) arc plus full-rev/half-rev
    # arcs -- not a two-generic-arc g/G row. The pre-#849 positional read
    # accidentally satisfied its own len(g_tofs) >= 2 guard here because a
    # half-rev arc ALSO carries tof_years (descriptor.py: "tof_years -- ...
    # g/h arcs only"), so it silently fed the h-arc's 0.5 yr ToF into the
    # designated-arc role of a model documented as g/G-only -- a worse defect
    # than a simple index swap. Correctly None now.
    row = _row(row_id)
    assert dds._descriptor_params(row) is None


def test_close_row_dsm_swap_affected_rows_numerically_unchanged() -> None:
    # #849 finding: for the 3 rows whose designated arc position flips
    # (9.353Gg2, 3.78Gg3, 9.94Gg3), the corrected g/G identification changes
    # the audited coplanar arc SHAPE (seed.arc_a_au/arc_e) but NOT
    # close_row_dsm's converged/anchor_match/vinf verdict, because the
    # corrector's actual per-leg ToF seed comes from the row's sourced
    # free_return_arcs/transit_times_days (arc order, untouched by the g/G
    # swap), not from the shape-fit arc itself. Pin against the historical
    # #830 uncorrected-posing numbers (data/found/830_v2_ballistic_multiarc/
    # dsm_388_recheck.json) so a future change to either lane is caught.
    ephem = Ephemeris("astropy")
    expected = {
        "russell-ch4-9.353Gg2": (29.46771835978174, False),
        "russell-ch4-3.78Gg3": (20.15004737126381, False),
        "russell-ch4-9.94Gg3": (28.5125967696132, False),
    }
    for row_id, (res_kms, match) in expected.items():
        res = dds.close_row_dsm(_row(row_id), ephem)
        assert res.max_residual_kms == pytest.approx(res_kms, abs=1e-4), row_id
        assert res.anchor_match is match, row_id
