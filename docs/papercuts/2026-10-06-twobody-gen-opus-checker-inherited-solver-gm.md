# The "independent" re-fly took its GM from the solver it was checking

- Date: 2026-10-06
- Agent: twobody-gen-opus
- Seen before: yes, as a class (the project's "independent cross-check" rule); first time for a constant

What happened: scripts/check_942_realeph_chain.py re-flew each leg with DOP853 using `sysm.mu`, the solver's own system. That system carried the heliocentric ideal model's GM, a 1-AU/365.25-d convention 3.8e-5 above the real solar GM, even at lambda = 1. Integrator, tolerances and code path were independent, but the physics constant was shared, so the checker agreed with the solver on a wrong model. The defect showed only when a new code path rebuilt the junctions against the real system directly (#942 D1, note 6.23).
Workaround: the checker takes the GM from the real ephemeris object; the solver's GM now follows lambda (c4ff9a41).
Suggested fix: an independent checker should take every physical constant (GM, body states, periods) from the reference model, never from the object under test, and should print the constants it used.
