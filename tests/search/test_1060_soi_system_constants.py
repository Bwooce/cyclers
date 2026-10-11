"""#1060: the sphere of influence comes from the system's own constants, registry as fallback.

Expected values are computed in the test from published numbers (JPL SSD PLU060 pages for the
Pluto system, accessed 2026-10-11: Pluto+Charon system GM 975.5 km^3/s^2, Nix GM 0.0015 km^3/s^2,
a 49300 km), not from the code under test.
"""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import pytest

from cyclerfinder.search.two_working_body import (
    CircularSystem,
    FlybyBody,
    MeanElementSystem,
    sphere_of_influence_km,
)

REPO = Path(__file__).resolve().parents[2]
MU_SYSTEM = 975.5
NIX_A, NIX_GM = 49300.0, 0.0015
NIX_SOI = NIX_A * (NIX_GM / MU_SYSTEM) ** 0.4


def _system(name: str) -> CircularSystem:
    period = 2.0 * math.pi * math.sqrt(NIX_A**3 / MU_SYSTEM)
    return CircularSystem(
        MU_SYSTEM,
        {name: (NIX_A, period, 0.0)},
        frozenset({name}),
        flyby_overrides={name: FlybyBody(name, NIX_GM, 18.0, 10.0)},
    )


def test_registered_moon_soi() -> None:
    assert sphere_of_influence_km(_system("Nix"), "Nix") == pytest.approx(NIX_SOI, rel=1e-12)


def test_unregistered_moon_does_not_raise() -> None:
    assert sphere_of_influence_km(_system("NotInRegistry"), "NotInRegistry") == pytest.approx(
        NIX_SOI, rel=1e-12
    )


def test_system_constants_win_over_registry() -> None:
    """A body on a circle of a different radius gets that radius, not the registry sma."""
    sysm = _system("Nix")
    per = sysm.bodies["Nix"][1]
    sysm.bodies["Nix"] = (60000.0, per, 0.0)
    assert sphere_of_influence_km(sysm, "Nix") == pytest.approx(
        60000.0 * (NIX_GM / MU_SYSTEM) ** 0.4, rel=1e-12
    )


def test_registry_fallback_for_system_without_bodies() -> None:
    """Earth's Laplace sphere of influence about the Sun is about 9.25e5 km (textbook value)."""
    soi = sphere_of_influence_km(MeanElementSystem(), "E")
    assert soi == pytest.approx(9.25e5, rel=0.02)


def test_unknown_body_without_constants_still_raises() -> None:
    with pytest.raises(KeyError):
        sphere_of_influence_km(MeanElementSystem(), "Planet9")


def _driver() -> object:
    spec = importlib.util.spec_from_file_location(
        "run_998_enumerate_1060", REPO / "scripts" / "run_998_enumerate.py"
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_pilot_reports_counts_and_no_assessment_errors() -> None:
    drv = _driver()
    res = drv.pilot_assess("pn", max_structures=40)  # type: ignore[attr-defined]
    assert set(res) == {"structures_tried", "zeros", "assess_errors"}
    assert res["structures_tried"] > 0
    assert res["assess_errors"] == 0


def test_pilot_would_catch_a_lookup_failure(monkeypatch: pytest.MonkeyPatch) -> None:
    """With the pre-#1060 lookup (registry only) an unregistered moon makes every assessment
    raise; the pilot must count those errors."""
    import cyclerfinder.search.two_working_body_enum as enum_mod

    def broken(system: object, code: str) -> float:
        raise KeyError(code)

    monkeypatch.setattr(enum_mod, "sphere_of_influence_km", broken)
    drv = _driver()
    res = drv.pilot_assess("pn", max_structures=60)  # type: ignore[attr-defined]
    assert res["zeros"] > 0
    assert res["assess_errors"] == res["zeros"]
