# The two-working-body pass definition and gauntlet accept conics that pass inside the central body or through a moon's SOI between encounters

- Date: 2026-10-07
- Agent: twobody-ext-opus
- Seen before: no (the k = 1-3 cells had no such gate-passer, so it never showed)

What happened: in the #973 gc k = 4 run, 8 of the 51 gate-passing physical cyclers have a leg whose minimum distance from Jupiter is below Jupiter's radius (71,492 km), down to 185 km from the centre. They pass the pre-registered #942/#943 definition (exact zero, turn gate, near-180 rule, encounter re-propagation) and the `scripts/gauntlet_942.py` gauntlet, because neither compares `r_min_km` with the primary's radius. `assess()` already computes r_min along every leg (fixed legs included since 6.30).
Second gap, same run: 14 of the 51 have a leg that passes inside a moon's Laplace SOI between the scheduled encounters (3 of them through the moon itself, 3-32 km from its centre). The patched-conic cycle is not valid there, and nothing in the gauntlet looks for it (`encounter_self_consistency` and the DOP853 re-fly check only the leg end points).
Workaround: `scripts/run_973_enumerate.py passes` finds them (interior distance minima along every leg); the #973 note classifies them post hoc as "J" (not physical) and keeps them in the table; nothing was loosened or dropped.
Suggested fix: add the unscheduled-pass check and an r_min > R_primary (plus a margin the owner sets) check to the pass definition or to the gauntlet of the two-working-body cells, with the radius from the registry; pre-register it before the next cell (gc k = 5-6 and ev k = 4-5 are about to run).
