# The private corpus clone was 87 commits behind, so 71 indexed PDFs looked missing

- Date: 2026-10-05
- Agent: main (found by corpus-review-fable)
- Seen before: no (the user's global rule on fetching every repo already exists)

What happened: cyclers_pdf had no upstream set, so `git status` showed nothing behind; a third of the review ran on digests only.
Workaround: fast-forwarded; set upstream on cyclers_pdf and cyclers.space.
Suggested fix: session-start check that fetches all three cycler repos and reports ahead/behind.

Disposition (main, 2026-10-05): Backlog: #964 (session-start fetch/ahead-behind check for the three cycler repos). Upstreams already set 2026-10-05.
