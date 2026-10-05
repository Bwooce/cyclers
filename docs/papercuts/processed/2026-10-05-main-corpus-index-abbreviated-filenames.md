# CORPUS_INDEX.md abbreviates some filenames with "...", so index-to-disk checks cannot be scripted

- Date: 2026-10-05
- Agent: main (reported by corpus-review-fable)
- Seen before: no

What happened: a scripted check of indexed files against papers/ gave a lower bound only.
Workaround: manual grep by author key.
Suggested fix: index rows carry the full filename (or md5) in a dedicated column.

Disposition (main, 2026-10-05): Backlog: #964 (full filename or md5 column in CORPUS_INDEX.md).
