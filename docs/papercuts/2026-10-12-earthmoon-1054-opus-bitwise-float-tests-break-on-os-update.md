# Exact float equality against stored solver records broke on a macOS update and blocked a catalogue window

- Date: 2026-10-12
- Agent: earthmoon-1054-opus
- Seen before: yes, #584 (two tests flipped on the new M3/Accelerate Mac, 2026-07-14)

What happened: c27609de added tests asserting `out["C"] == rec["C"]` for converged free-C shooting seeds against records in committed jsonl files. They passed at 09:25 AEDT. macOS 27.0.1 installed at the 09:37 reboot, which changed Accelerate/libm. The same code then landed 1.9e-7 and 1.6e-9 away in C, because the free-C landing point drifts along the solution family by rounding noise over a smallest singular value of about 4e-8. Window 4 was blocked while the lead bisected.
Workaround: replaced the equalities with fixed bounds derived from the arc rounding noise and the last-step singular value, with a wrong-row positive control (note 2026-10-11-1054 sec. 10.1, scripts/probe_1054_float_sensitivity.py).
Suggested fix: a test that compares an iterative solver's float output with a value stored on disk states a derived tolerance, not `==`. Bitwise comparisons are for values computed twice in the same process, or for closed-form values. When a long-stable float test fails, check `/Library/Receipts/InstallHistory.plist` and `sysctl kern.boottime` for an OS update first.
