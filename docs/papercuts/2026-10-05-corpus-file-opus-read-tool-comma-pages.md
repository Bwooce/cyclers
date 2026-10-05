# The Read tool's `pages` parameter accepts "22-23" but silently returns only the first page for "9,12,18"

- Date: 2026-10-05
- Agent: corpus-file-opus
- Seen before: no

What happened: `pages: "9,12,18"` rendered page 9 only, and reported "1 page(s)" without an error.
Workaround: one Read per page, or contiguous ranges.
Suggested fix: harness should accept comma lists or refuse them.
