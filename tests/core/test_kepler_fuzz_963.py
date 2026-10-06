"""#963: random-state fuzz of :func:`cyclerfinder.core.kepler.propagate` in production units.

The bracket-safeguarded Newton of #963 must converge on every arc the searches produce. A
first version passed every canonical-units test and still failed on a Jovian arc in km
(its fallback step had an absolute size of 1), so these fuzzes run in km at the two scales
the searches use: heliocentric (|r| ~ 1.5e8 km, v 5-45 km/s, |dt| up to 2 yr) and Jovian
(|r| ~ 1.5e6 km, v 1-60 km/s, |dt| up to 40 d), elliptic and hyperbolic mixed.

Measured before #963 (2026-10-06): the heliocentric fuzz failed on 107 of 200,000 states, the
Jovian on 21 of 100,000; after: none. The fast tests run the first 2,000 states of each; the
slow ones the full sets.

Accuracy is held by the reference tests in ``test_kepler.py``; here the checks are that every
arc converges and that energy, angular momentum and a backward round trip are kept to loose
bounds that catch a wrong root or garbage. The bounds are set from the measured worst cases
(2.3e-10, 2.8e-8 and 3.6e-7, all on orbits that plunge to a small fraction of r0), times
about 30.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.kepler import propagate

_HELIO = (132712440018.0, 1.5e8, 5.0, 45.0, 6.3e7, 1)  # mu, r scale, v lo, v hi, |dt| max, seed
_JOVIAN = (126686534.0, 1.5e6, 1.0, 60.0, 3.5e6, 7)


def _fuzz(
    n: int, mu: float, rscale: float, vlo: float, vhi: float, tmax: float, seed: int
) -> list[str]:
    rng = np.random.default_rng(seed)
    bad: list[str] = []
    for i in range(n):
        r0 = rng.normal(size=3) * rscale
        r0[2] *= 0.1
        v0 = rng.normal(size=3)
        v0 *= rng.uniform(vlo, vhi) / float(np.linalg.norm(v0))
        v0[2] *= 0.1
        dt = float(rng.uniform(-1.0, 1.0)) * tmax
        try:
            r, v = propagate(r0, v0, dt, mu)
            r_back, _ = propagate(r, v, -dt, mu)
        except Exception as exc:
            bad.append(f"{i}: {type(exc).__name__}")
            continue
        r0n = float(np.linalg.norm(r0))
        e0 = 0.5 * float(v0 @ v0) - mu / r0n
        e1 = 0.5 * float(v @ v) - mu / float(np.linalg.norm(r))
        h0 = np.cross(r0, v0)
        if abs(e1 - e0) > 1e-8 * max(abs(e0), mu / r0n):
            bad.append(f"{i}: energy")
        if float(np.linalg.norm(np.cross(r, v) - h0)) > 1e-6 * float(np.linalg.norm(h0)):
            bad.append(f"{i}: angular momentum")
        if float(np.linalg.norm(r_back - r0)) > 1e-5 * r0n:
            bad.append(f"{i}: round trip")
    return bad


@pytest.mark.parametrize("case", [_HELIO, _JOVIAN], ids=["heliocentric", "jovian"])
def test_kepler_fuzz_fast_963(case: tuple[float, float, float, float, float, int]) -> None:
    bad = _fuzz(2000, *case)
    assert not bad, f"{len(bad)} bad states:\n" + "\n".join(bad[:20])


@pytest.mark.slow
@pytest.mark.parametrize(
    ("case", "n"), [(_HELIO, 200_000), (_JOVIAN, 100_000)], ids=["heliocentric", "jovian"]
)
def test_kepler_fuzz_full_963(case: tuple[float, float, float, float, float, int], n: int) -> None:
    bad = _fuzz(n, *case)
    assert not bad, f"{len(bad)} bad states:\n" + "\n".join(bad[:20])
