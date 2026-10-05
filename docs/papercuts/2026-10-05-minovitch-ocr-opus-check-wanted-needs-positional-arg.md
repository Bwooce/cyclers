# check_wanted_vs_corpus.py needs the wanted-list path as an argument

- Date: 2026-10-05
- Agent: minovitch-ocr-opus (#966)
- Seen before: no

What happened: the #966 brief said to run `uv run python scripts/check_wanted_vs_corpus.py`. With no argument the script exits with "the following arguments are required: wanted".
Workaround: `uv run python scripts/check_wanted_vs_corpus.py docs/notes/2026-10-05-960-wanted-papers.md`.
Suggested fix: default `wanted` to the current wanted-list note, or quote the full command in briefs.
