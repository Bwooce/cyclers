# A low-load timing pilot underestimated a CI-loaded enumeration by 5x, and the 10-min tool limit forced a mixed shard layout

- Date: 2026-10-07
- Agent: twobody-ext-opus
- Seen before: no

What happened: the #973 gc k = 4 pilot (24 random structures, load average 3-4) gave 0.21 s per structure. The full run, under the CI pytest suite (6 workers, load average 17-30), averaged 1.07 s. The first 8-way shards then ran 8-10 min each, past the 10-min tool-call limit (one call was moved to the background by the harness), so the rest of the cell was split 16-way. The final layout mixes shards of 8 and of 16.
Workaround: `scripts/run_973_enumerate.py check` validates coverage from the structure keys (every structure once, no duplicate zero lines), so a mixed layout is safe; the note states the layout.
Suggested fix: pilot under the load the run will see (or scale by the current load average), and size shards for the tool limit at the loaded rate. A `check`-style coverage validator belongs in `run_942_enumerate.py` for every cell.
