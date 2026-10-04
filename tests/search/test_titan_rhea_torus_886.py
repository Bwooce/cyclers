"""#886 part 1: the #882 invariant-circle corrector against Kumar IAC-25-C1.9.6 (Titan-Rhea).

Expected values trace to the paper (filed in the private paper corpus as
``kumar-2025-analysis-unstable-resonant-orbits-saturn-tour-design-titan-rhea-IAC-25-C1.9.6-doi-10.52202-083087-0076.pdf``):
the printed constants, the printed Tp/T range of the Titan 3:2 low/mid-e family, the six
printed secondary-resonance ratios, and the qualitative statement that tori persist away from
those ratios and fail near them. Node counts, steps and the scan itself are in
``scripts/screen_886_titan_rhea_torus_check.py`` (staged, not in the default suite); results in
``docs/notes/2026-10-04-886-titan-rhea-torus-published-check.md``.
"""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path
from types import ModuleType

import numpy as np
import pytest
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "screen_886_titan_rhea_torus_check.py"
FAMILY = ROOT / "data" / "found" / "886_titan_rhea_torus" / "family.json"

# Printed by the paper (p. 2, p. 12, p. 13).
PAPER_TP = 2.48376
PAPER_RATIO_HI = 0.1952735
PAPER_LISTED = {(4, 21), (5, 26), (6, 31), (7, 36), (8, 41), (9, 47)}


@pytest.fixture(scope="module")
def mod() -> ModuleType:
    spec = importlib.util.spec_from_file_location("screen_886", SCRIPT)
    assert spec is not None and spec.loader is not None
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


@pytest.fixture(scope="module")
def family(mod: ModuleType) -> object:
    # The stored rows are used only as Newton seeds; every value asserted is recomputed.
    return mod.Family(mod._load(FAMILY.name))


def test_model_period_is_printed_tp(mod: ModuleType) -> None:
    sysp = mod.system(mod.MU3)
    assert sysp.ganymede_synodic_period == pytest.approx(PAPER_TP, rel=1e-14)
    # Rhea inside Titan's orbit advances in Titan's frame
    assert sysp.omega_gan > 0.0
    # Kepler radius for that rate, against the printed "about 0.4315"
    assert abs(sysp.a_gan - 0.4315) < 5e-4


def test_farey_ratios_below_50_are_the_papers_six(mod: ModuleType) -> None:
    got = set(mod.farey_in_range(0.1898475, PAPER_RATIO_HI, 49))
    assert got == PAPER_LISTED


def test_family_upper_ratio_matches_paper(mod: ModuleType, family: object) -> None:
    """The interior minimum of T gives max Tp/T; the paper prints 0.1952735.

    Agreement is limited by the printed Tp (6 significant figures, relative rounding up to
    2e-6, i.e. about 4e-7 in the ratio).
    """

    def period(x: float) -> float:
        return float(family.at(x)[1])  # type: ignore[attr-defined]

    res = minimize_scalar(
        period,
        bounds=(-0.690, -0.675),
        method="bounded",
        options={"xatol": 1e-7},
    )
    ratio_max = PAPER_TP / float(res.fun)
    assert abs(ratio_max - PAPER_RATIO_HI) < 4e-7


def test_off_resonance_member_persists(mod: ModuleType, family: object) -> None:
    """Low-e member far from every listed ratio: the circle continues to the printed mu3."""
    out = mod.continue_member(-0.7157, family, levels=(151,), with_bundles=False)
    assert out["seed_residual"] <= 1e-10
    assert out["status"] == "persist"
    assert out["final_residual"] <= 1e-10
    assert out["offnode_residual"] <= 1e-8
    ratio = out["ratio"]
    assert min(abs(ratio - p / q) for p, q in PAPER_LISTED) > 1e-4


def test_at_listed_resonance_member_fails(mod: ModuleType, family: object) -> None:
    """Tp/T = 7/36 exactly (segment 1, C about 3.028): no circle reaches even a small mass.

    A coarse step floor keeps the test cheap; the staged scan uses 1/1024 and a refined rerun.
    """
    seg = family.segments()[1]  # type: ignore[attr-defined]
    x = family.solve_ratio(7 / 36, seg)  # type: ignore[attr-defined]
    assert x is not None
    out = mod.continue_member(x, family, levels=(151,), h_floor=1.0 / 64.0, with_bundles=False)
    assert abs(out["ratio"] - 7 / 36) < 1e-10
    assert out["seed_residual"] <= 1e-10  # the representation is fine; the torus is not
    assert out["status"] != "persist"
    assert out["frac_reached"] < 0.05
    assert math.isfinite(out["final_residual"])
    assert np.isfinite(out["T"])
