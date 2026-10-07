# gauntlet_942.py and scan_320 import literature_check at module level, so they cannot be reused by a task that must not load it

- Date: 2026-10-07
- Agent: pluto-smallmoons-sonnet
- Seen before: no

What happened: `scripts/gauntlet_942.py` (the DOP853 `cross_check`) and `scripts/scan_320_epoch_aware_moon_systems.py` (the #320 sweep used as a control) both import `cyclerfinder.search.literature_check` at the top. A task whose brief forbids importing that module (the file is under rework in #972) could not reuse either script's functions.
Workaround: #998 loads each script with an inert stand-in for `literature_check` placed in `sys.modules` for the duration of the import, then restores it. The anchor-overlap field of the #320 control is therefore not reproduced.
Suggested fix: move `cross_check`, `fly` and the DOP853 re-fly into a small module (for example `verify/` or `search/two_working_body_refly.py`) with no literature import, and import the literature module lazily inside the functions that use it.
