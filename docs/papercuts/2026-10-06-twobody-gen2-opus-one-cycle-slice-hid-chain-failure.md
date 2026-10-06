# A 1-cycle validation slice did not cover the 10-cycle launch, and a half-rev root choice hid behind it

- Date: 2026-10-06
- Agent: twobody-gen2-opus
- Seen before: yes. The same kind of miss as 2026-10-06-twobody-gen-opus-slice-missed-task-number-branch.md.

What happened: I validated the GanEur#316 real-ephemeris chain with `--n-cycles 1` and then sent launch commands with `--n-cycles 10`. Both launches died at once in `initial_fixed_params` (`assert res is not None`).

Cause:
- `_solve_half_rev_e` returned the first root in e when scanning up from 0.
- For a (3, 1) half-rev the body's own circle (e = 0) is also a root.
- In cycles 2 onward, the blended Ganymede radius was a little ABOVE the circular one. That moved the circle's root to e of about 0.002-0.005, so the scan took it instead of the leg's conic (e = 0.285). The result was 87-deg turn demands, or no directions at all.
- Cycle 1's radius happened to be below the circular one, so the 1-cycle slice passed.
- The same knife-edge also existed in the ideal model, where f(0) is round-off. Reassessment showed it flipped some ev zeros from "fail" to "no-directions", but no pass.

Workaround/fix:
- The root on the leg's side of the flight-time minimum is now kept (test against R-S Table 3).
- The bare assert is replaced by an error that names the block, the body, the time and |V_inf|.
- The shoot gets a nearest-geometry seed where the blend has no minimax solution.

Suggested fix: validate long-compute launches with the EXACT launch flags (same --n-cycles and epoch options; cut only the wall time, e.g. --shoot-restarts), not a shorter chain.
