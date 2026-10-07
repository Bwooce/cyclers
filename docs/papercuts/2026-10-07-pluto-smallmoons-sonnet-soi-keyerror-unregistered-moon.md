# sphere_of_influence_km raises KeyError for a moon that is not in the satellite registry, after the whole run

- Date: 2026-10-07
- Agent: pluto-smallmoons-sonnet
- Seen before: no

What happened: `two_working_body.sphere_of_influence_km` looks the body up in PLANETS and SATELLITES and raises `KeyError` otherwise. `assess` calls it for every encounter body, so for a cell with an unregistered target (Kerberos, Styx are not in `core/satellites.py`) every exact zero failed assessment with `KeyError('Kerberos')`. The run driver catches this per zero and writes status "error", so the 6532-structure Kerberos run finished with 3111 assessment errors and zero usable results, after about 35 minutes of chunks. The pilot (Nix) had both bodies registered, so it did not show it.
Workaround: the #998 driver adds the moon to the in-process `SATELLITES` dict before building the system, and the Kerberos run was discarded and redone.
Suggested fix: either let the generator take the SOI from the system's own constants (FlybyBody already carries mu; the semi-major axis is in `CircularSystem.bodies`), or make the driver run one `assess` on a pilot zero before the sweep. A pilot for every new cell should include a zero count and an assessment-error count, not only the solver timing.
