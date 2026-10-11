"""#1064: scripts/realeph_998.py's EphSystem must not take the Keplerian 1:1 shortcut at Charon.

Charon moves about the Pluto-system barycentre at a speed the conic about the system GM does not
give. With is_keplerian True (the old behaviour) a full-revolution 1:1 return built by
two_working_body.resonant_circle misses Charon by about 19,000 km after one period; with the fix
(is_keplerian False for Charon, so the general a_sc path is used) it re-flies onto Charon.
Needs plu060.bsp on disk (skipped otherwise).
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from cyclerfinder.search.two_working_body import (
    ResonantLeg,
    kepler_step,
    resonant_circle,
    resonant_vec,
)

REPO = Path(__file__).resolve().parents[2]


def _load(name: str, rel: str) -> Any:
    spec = importlib.util.spec_from_file_location(name, REPO / rel)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_is_keplerian_false_for_charon_only() -> None:
    mod = _load("realeph_998_1064", "scripts/realeph_998.py")
    eph = object.__new__(mod.EphSystem)  # no kernel needed for the predicate
    assert eph.is_keplerian("Charon") is False
    for moon in ("Styx", "Nix", "Kerberos", "Hydra"):
        assert eph.is_keplerian(moon) is True


def _system() -> tuple[Any, Any]:
    mod = _load("realeph_998_1064b", "scripts/realeph_998.py")
    spice = pytest.importorskip("spiceypy")
    from cyclerfinder.verify.spice_kernels import ensure_pluto_kernel

    try:
        kernel = ensure_pluto_kernel()
    except RuntimeError:
        pytest.skip("plu060.bsp not on disk")
    spice.furnsh(str(mod.LSK))
    spice.furnsh(kernel)
    drv = _load("run_998_enumerate_1064", "scripts/run_998_enumerate.py")
    ideal = drv.pluto_system("Nix")
    et_a = spice.str2et("2030-01-01T00:00:00")
    frame = mod.mean_plane_frame(et_a, et_a + 4 * 365.25 * 86400.0)
    et0 = mod.next_conjunction(frame, "Nix", et_a)
    return mod, mod.EphSystem(ideal, et0, frame)


def _miss_km(eph: Any, vinf: float, phi: float = 0.3) -> float:
    leg = ResonantLeg("Charon", 1, 1)
    circ = resonant_circle(eph, leg, 0.0, vinf)
    assert circ is not None
    u = resonant_vec(circ, phi)
    r, w = eph.state("Charon", 0.0)
    period = eph.period_s("Charon")
    r_end, _ = kepler_step(r, w + u, period, eph.mu)
    r_charon, _ = eph.state("Charon", period)
    return float(np.linalg.norm(r_end - r_charon))


@pytest.mark.parametrize("vinf", [0.15, 0.2, 0.3])
def test_one_to_one_leg_refly_onto_charon_with_and_without_the_fix(
    vinf: float, monkeypatch: pytest.MonkeyPatch
) -> None:
    mod, eph = _system()
    assert _miss_km(eph, vinf) < 1.0  # measured 0.196 km (Charon's own eccentricity)
    monkeypatch.setattr(mod.EphSystem, "is_keplerian", lambda self, code: True)
    assert _miss_km(eph, vinf) > 10_000.0  # the shortcut: about 18,000-19,000 km
