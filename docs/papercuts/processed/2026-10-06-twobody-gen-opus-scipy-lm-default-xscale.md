# Passing x_scale=1.0 to SciPy's lm silently changed the solver

- Date: 2026-10-06
- Agent: twobody-gen-opus
- Seen before: no

What happened: I added a `--shoot-method` option to the #942 chain tool and wrote `x_scale="jac" if method == "trf" else 1.0`, believing 1.0 was the lm default. Since SciPy 1.16 the lm default is `x_scale="jac"` (the repo has 1.17.1). The shoot phase then stalled at a residual of about 1e-3 to 1e-2 where it had converged before. I first read this as a structural "no root" result. The same captured problem converged in 20 evaluations when called without the kwarg.
Workaround: pass `x_scale="jac"` explicitly (commit after 1b8cfbfa). The 1-cycle sweep made with 1b8cfbfa is diagnostic only and is void.
Suggested fix: when adding an option that passes a kwarg the code used to leave at its default, check the library's current default (`inspect.signature`) and keep it. A stall right after a "solver settings only" change points at the change first.

Disposition: Promoted: lesson added to docs/team/coordination.md known hazards (check a library's CURRENT default with inspect.signature before passing a kwarg the code left at default).
