# A units-dependent bracket fallback passed a canonical-units sweep and broke a km-unit search

- Date: 2026-10-06
- Agent: ci-keeper-opus
- Seen before: no

What happened: the first #963 fix (5c8acd2d) stepped out of an open bracket by an absolute `2 lo + 1`. In canonical units (mu = 1, chi ~ 1) that step is sensible. In km at Jupiter (chi ~ 1e4) it is tiny. The whole sweep and the 200,000-state heliocentric comparison passed, and then the #943 ge enumeration crashed on a Jovian hyperbolic arc in km. The same state converges when rescaled to canonical units.
Workaround: reproduced the arc with a wrapper around `propagate` in the running driver. Fixed the step (start from the parabolic bootstrap, which has the units of chi and the sign of dt) and added the rtsafe rule.
Suggested fix: run numerics regression tests at the scales production uses (km and km^3/s^2 at the Sun and at Jupiter), not only in canonical units. Grep new code for absolute constants (`+ 1`, `max(|x|, 1)`) that silently assume a unit. `tests/core/test_kepler_fuzz_963.py` and the scaled hyperbolic sweep in `test_kepler.py` now do this for kepler.
