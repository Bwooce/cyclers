"""Descriptor -> DSM seed adapter (plan 2026-06-10, Component 1)."""

from __future__ import annotations

from typing import Any

import numpy as np

import cyclerfinder.search.dsm_descriptor_seed as dds
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


def test_no_descriptor_row_returns_none() -> None:
    # An ocampo row has the n.m.k summary format, no per-arc g/G descriptor.
    catalog = load_catalog()
    ocampo = next(e for e in catalog.entries if e.id.startswith("russell-ocampo"))
    assert dds.seed_dsm_chain_from_descriptor(ocampo.raw) is None
