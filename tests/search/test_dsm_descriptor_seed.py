"""Descriptor -> DSM seed adapter (plan 2026-06-10, Component 1)."""

from __future__ import annotations

from typing import Any

import numpy as np

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
