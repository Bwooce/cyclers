# The #975 loop tracer accepted corrector landings on other families, and two negative controls passed on it

- Date: 2026-10-12
- Agent: earthmoon-975-opus
- Seen before: yes, in the same tracer: the #975 N1 inexact-closing false positive (2026-10-11), a closing step accepted without landing on the start.

What happened: the fixed-mu pseudo-arclength tracer (`trace_loop` in `src/cyclerfinder/search/varin_closed_families.py`) accepted any converged corrector point. At a hairpin of the i_1 oval at the Earth-Moon mu, a 0.02 predictor step overshot and the corrector converged 1.02 away on the 3-fold cover of a T~ 1 orbit. The P2 up trace then ran 500 members toward a P1 collision and never closed. The same scan found jumps in the N1 and N2 b-branch traces (step 2.0 and 1.53). Those traces had "passed" N1 and N2 by not closing, so the two negative controls were reopened.
Workaround: tool amendment A2 (#975 note sec. 11.4, commit deae8e75). A step is accepted only if the corrector lands within 0.5|step| of the predictor and the tangent turns by at most 10 deg; otherwise the step is halved. Closure also needs a forward closing tangent. The same guard is in the extremum-curve follower. A regression test is built from the committed jump.
Suggested fix: every continuation in the project (the multiarc, DSM and manifold continuations too) should check corrector displacement and tangent turn, not only convergence. Before a trace is used as evidence, scan its step lengths against the step cap: a step much longer than the cap is a jump. A negative control that passes by "not closing" or "not found" must first be checked for jumps, because a jumping tracer is a false-negative generator.
