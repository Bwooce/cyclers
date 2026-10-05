# core.kepler.propagate raises KeplerConvergenceError on an ordinary single-rev heliocentric arc

- Date: 2026-10-05
- Agent: twobody-gen-opus
- Seen before: yes. #934 added a parabolic-bootstrap retry for near-parabolic orbits. This case is elliptic (alpha > 0) and still fails.

What happened: while #942 sampled leg extents, `propagate` failed with "chi=4.831545e+04, residual=0" on a 176-d Venus-to-Earth arc. Reproducer: r0=[-88039995.65452953, -62916481.70866319, 0.0] km, v0=[-19.984711843308432, -29.60574045546609, 0.0] km/s, mu=132717453059.67786, dt=14953799.40036204 s. Any caller that catches the error as "arc infeasible" turns it into a false negative.
Workaround: `two_working_body.kepler_step`, an eccentric-anomaly Lagrange f and g propagator for elliptic orbits. Where `propagate` converges, the two agree to < 1 mm and 1e-9 km/s (test in tests/search/test_two_working_body.py).
Suggested fix: add an elliptic eccentric-anomaly path, or a bracketed (bisection-safeguarded) Newton, to `propagate`; add this reproducer as a regression test; grep for callers that catch KeplerError and treat it as infeasible.

Disposition (main, 2026-10-05): Backlog: #963 (correctness bug in core.kepler.propagate; after the fix, re-check negatives that hit KeplerError).
