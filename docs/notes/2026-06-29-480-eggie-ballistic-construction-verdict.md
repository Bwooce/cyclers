# #480 ballistic-construction verdict — the equal-V∞ bug is fixed, but DISCOVERY still finds the wrong resonant member (a real construction bug, not a search-tuning gap)

**Date:** 2026-06-29 (salvaged from a stalled agent run — findings in the runlog,
nothing the agent committed; this note + the analysis are the coordinator's).
**Builds on:** the forward-verify correction
(`docs/notes/2026-06-29-480-eggie-forward-verify-correction.md`).

## What the ballistic-construction search found (192 branch/rev plans)

- **BEST FEASIBLE: NONE.** No plan produced an all-flybys-above-surface ballistic close.
- **Best by distance-to-Table-4 V∞** (plan `0s1l1l2l`): V∞ = E 9.124 / **G1 7.071 /
  G2 7.071** / Io 8.378 / E 9.125 — DEAD ON Table 4, equal Ganymede recovered; ΔV
  0.172 m/s (ballistic). BUT all flyby altitudes negative (G1/G2 −2627, Io −1817,
  E −1559 km) and ToFs `[0.574, 17.517, 7.195, 12.052]` (G→G = 17.5 d ≈ 2× Table-4's
  8.60 d). → on-V∞ but a DIFFERENT (2:1-G) resonant member, sub-surface.
- **Gentle-bend (seed e=0.64):** V∞ on-target, ToFs `[5.46, 8.43, 6.80, 10.63]`
  (G→G ≈ 8.4 d, near Table-4), ΔV 2.83 m/s; **Ganymede #2 altitude flips POSITIVE
  (2432 km)** but G1/Io/Europa-seam still sub-surface; leg-0 = 5.46 d vs Table-4 1.59 d.
- The rigid conic guess itself has near-Table-4 ToFs `[1.607, 9.006, 6.895, 10.509]`
  but WRONG V∞ at Io (6.23 vs 8.38) and Europa-return (11.65 vs 9.12).

## Diagnosis — the bug is in DISCOVERY (crossing-topology), and hardcoding ToFs is NOT the fix

The machinery is sound: V∞, ΔV, altitude, and the equal-Ganymede property all compute
correctly (the search reproduces Table-4 V∞ exactly when it lands on the right energy).
What is WRONG is **which resonant trajectory the construction discovers**:

- At the rigid conic guess's near-Table-4 ToFs, the construction computes the WRONG V∞
  at **Io and the 2nd Europa** (Io 6.2 vs 8.38; Europa-return 11.6 vs 9.12). So at the
  paper's own geometry our pipeline does not reproduce the paper's encounters — a direct
  contradiction ⇒ a construction bug, not "EGGIE is hard."
- Root cause (localized): the EGGIE **crossing-topology / rev-branch assignment**
  (`resonant_conic._EGGIE_REVS = (0,1,2,4,5)` + the in/out `_EGGIE_NU_KEYS`) maps the
  Io and repeated-Europa/Ganymede encounters to conic crossings that are NOT the paper's
  member. Because those nodes' V∞ are wrong at the correct ToFs, when the optimizer then
  forces all V∞ onto Table 4 it must distort the orbit — escaping to a different member
  (G→G 17.5 d) with sub-surface flybys. The sub-surface altitudes and wrong ToFs are
  SYMPTOMS of the topology mis-assignment, not separate bugs.

**Why NOT hardcode the paper's ToFs:** plugging in `[1.59, 8.60, 7.34, 10.69]` would be
circular — feeding the paper's answer in, not reproducing it — and would leave us with a
tool that cannot DISCOVER cyclers (the whole point of the conic initial-guess method).
The legitimate use of the paper's ToFs is only as a positive-control to localize the bug
(primitives vs discovery); that control already passed for the primitives (V∞/ΔV/altitude
formulas are right), pointing the finger squarely at the crossing-topology discovery step.

## The genuine fix (next, not a re-run)

Debug `resonant_conic`'s crossing selection so the discovery FINDS the paper's EGGIE
member on its own: the Io and repeated-body encounters must map to the conic crossings
whose ToFs come out `≈ [1.59, 8.60, 7.34, 10.69]` AND whose V∞ are equal-in/out and
above-surface. Concretely: re-derive `_EGGIE_REVS` / `_EGGIE_NU_KEYS` from the paper's
member (5 s/c revs, the SHORT 1.59 d leg-0 E→G hop, the ~1:1 G→G return — not 2:1), and
enforce repeated-encounter self-consistency (the moon must actually be at the conic
crossing at the crossing time). Then the construction should discover a feasible,
on-Table-4, ballistic EGGIE with no hardcoded ToFs. Scratch drivers `scripts/_eggie_*.py`
(removed); search runlog captured the numbers above.

Process note: the agent stalled (background-detached search; no re-wake) — salvaged per
[[feedback_long_agents_commit_incrementally]]; future search work runs foreground-bounded.

## ROOT CAUSE (debug `scripts/_eggie_topo_{debug,discover}_480.py`, 2026-06-29)

Per-encounter geometry dump (rigid guess, e=0.62) — `|dpos|` = distance from the moon
to the conic crossing at the encounter time:

| encounter | rev | t (d) | \|dpos\| km | V∞ |
|---|---|---|---|---|
| Europa-depart | 0 | 0.0 | 0 | 9.08 ✓ |
| Ganymede-1 | 1 | 1.61 | 0 | 6.78 ✓ |
| **Ganymede-2** | 2 | 10.61 | **391,088** | 3.06 ✗ |
| Io | 4 | 17.51 | 0 | 8.35 ✓ |
| **Europa-2 (return)** | 5 | 28.02 | **237,764** | 13.83 ✗ |

The three DISTINCT-moon encounters are exact; both REPEATED encounters are hundreds of
thousands of km off. Three nested facts:

1. **The construction sets each moon's phase from its FIRST encounter only** (`_moon_phase_ics`)
   and never verifies repeated encounters — so the 2nd Ganymede / 2nd Europa silently
   drift off the conic. That is the bug class.
2. **The repeated revs are mis-assigned:** the resonance-consistent Ganymede repeat is
   rev-1-inbound (t≈10.1 d, the moon only dθ≈4.68° = the designed 5.2°/synodic shift, ≈86,000 km
   off), but the topology picks the rev that is ~21°/391,000 km off.
3. **Deeper — a rigid single conic cannot host EGGIE:** the moons are ~5° off every
   crossing BY DESIGN (the 5.2° shift), and there is NO good Europa-return crossing near
   t≈28 d (best Europa crossings are t=10.7 d and t=38.7 d). So the conic is only a SEED;
   the paper's Monte-Carlo Lambert search over (phase, ToFs) finds the actual member where
   the legs connect the real shifted moon positions ballistically.

**The genuine fix (no hardcoding):** seed the optimizer (`search/eggie_ballistic.py`,
committed — `build_legs`/`ballistic_residual`/`refine`, no hardcoded ToFs) with the
RESONANCE-CONSISTENT topology (correct repeat revs) + enforce a per-encounter
self-consistency invariant (every body, incl. repeats, within its SOI of its node), then
search phase+ToFs within ±10% of the moon periods. The construction must DISCOVER the
member; the conic only seeds it. Debug scratch removed.
