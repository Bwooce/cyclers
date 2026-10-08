# The literature gate's "closing repeat" rule misreads a time-ordered flyby list that starts and ends at the same body

- Date: 2026-10-08
- Agent: twobody-ext-opus
- Seen before: no (gc-1 and gc-2 also start and end with Ganymede, but dropping the last encounter still leaves a C-C pair, so the result did not change)

What happened: in the #1025 literature step, gc6-0 (cycle G G C, built by `scripts/litcheck_942_943_scope.py::sig_of` as the time-ordered sequence (Ganymede, Callisto, Ganymede)) read "published" (0.85) via the Campagnola 2019 GCGC anchor. `literature_check._has_consecutive_same_body` applies #972 F3: a final encounter equal to the first is a catalogue "closing repeat" and is dropped. In a `sig_of` signature the last Ganymede is a distinct flyby, not a repeat of the first, so the drop turns G G C into an alternating (G, C). The GCGC anchor's `alternating_scope` then no longer excludes it. With the cycle started at Callisto, (C, G, G), the same member reads "inconclusive" (R-S G-C, F7) like the other 19 gc members.
Workaround: the pre-registered result is recorded as run, the rotated-start result is reported beside it as a diagnostic, and the #1025 notes call it a signature-convention artefact.
Suggested fix: let a signature declare whether its sequence includes a closing repeat (`CandidateSignature` field, default as today for catalogue strings; `sig_of` sets "no"), or make `sig_of` start every sequence at an encounter whose predecessor is a different body. Re-run the 26 controls and the #942/#943/#1025 candidates after the change.
