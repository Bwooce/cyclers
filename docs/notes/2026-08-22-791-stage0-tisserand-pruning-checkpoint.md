# `#791` Stage 0 — Tisserand/energy-graph pruning checkpoint: pruning too weak, STOP

**Update 2026-08-22 (same day): the chained/energy-consistency filter was tried next — it fails
too, and for a sharper, more decisive reason.** See "Addendum: the chained filter" below.

**Task:** `#791`, registered 2026-08-08 (Campaign 3 of the combinatorial-search survey, "moon-tour
encounter-sequence enumeration"). Per `#858`'s review (Sec. 5.4), fund Stage 0 only first: pure
Tisserand/energy-graph pruning of the closed-tour sequence alphabet over the Galilean moons, before
evaluating a single cell. Mandatory decision rule from that review: *"if 10^3 sequences survive at
length ≤ 8, the pruning is too weak to proceed."*

**Verdict: pruning is too weak.** 3269 of 3269 new candidate sequences survive (100%) — the pure
pairwise Tisserand screen provides **zero discriminating power** on this alphabet. Per `#858`'s own
stated rule this is a STOP-and-report checkpoint: do **not** auto-fire Stage 1's patched-conic
enumeration on the full, unpruned 3269-sequence set. This is not a bug in the gate (verified below)
— it is a genuine physical fact about the Galilean system.

---

## What was run

`scripts/scan_791_stage0_sequence_pruning.py`. Alphabet: the four Galilean moons (Io, Europa,
Ganymede, Callisto) — matching every prior campaign in this lane (`#318`/`#501`/
`releg_moontour.py`). Sequences: closed tours (`seq[0] == seq[-1]`), adjacency-distinct, 3-8
encounters (matching `#791`'s own registration language and `#501`'s tag convention, e.g.
"EGE" = 3 encounters / 2 legs).

- **Total enumerated closed sequences (length 3-8):** 3276.
- **Already-stamped**: 7 — `#501`'s own 6 (`scan_501_broadened_joint_search.py`'s EGE, GCG, EGCE,
  IEI, IEGI, EGCGE) **plus** the earlier `#318` smoke-test's `jovian.CGCEC`
  (Callisto-Ganymede-Callisto-Europa-Callisto), found only by checking `#318`'s own driver script
  directly rather than trusting the `#791` registration's own "7 hand-picked sequences" summary at
  face value (the `#501` module alone only lists 6 — the two together account for the "7"). All 7
  correctly excluded from the "new" count (capability-subsumption skip, not re-swept).
- **Pure Tisserand pruning gate:** for every leg, does `tisserand.linkable` find a shared V∞ in the
  project's own established Jovian probe band (`releg_moontour._VINF_PROBE_KMS` =
  4,5,6,8,10,12,15 km/s) at which the two moons' constant-V∞ contours intersect? Deliberately
  narrower than `moon_prune.moon_leg_admissible` (no VILM ΔV-floor, no bend check) — this is the
  pure energetic screen Fable's Stage 0 spec asked for, with no budget assumption baked in this
  early.
- **Positive control:** the known-closing IEGI family (Io-Europa-Ganymede-Io) survives the pure
  Tisserand screen, as it must — a real closer cannot exist if its legs are not even energetically
  linkable.

**Result: 3276/3276 total survive, 3269/3269 new (non-stamped) survive. Zero pruned.**

## Why: the Galilean alphabet is a Tisserand-complete graph at ~3-4 km/s, not a sparse one

Independently re-verified before trusting the "0 pruned" result — a screen that prunes nothing is
as suspicious as a search that closes everything ("it closed!" and "it pruned nothing!" are the
same species of danger signal). Two checks:

1. **Falsifiability**: `linkable()` genuinely discriminates in general — e.g. Io-Callisto is linked
   at V∞=4 km/s but NOT at 8, 15, 30, or 50 km/s; Europa-Ganymede is linked at 0.1-4 km/s but not
   8+. The gate is not trivially always-True.
2. **Padding artifact check**: re-ran all 6 possible Galilean pairs with the moons' raw SMA-derived
   `a` range and **no** 25% padding (`moon_prune._a_range_au`'s own padding factor, removed for
   this check). Result unchanged: every one of the 6 possible pairs (Io-Europa, Io-Ganymede,
   Io-Callisto, Europa-Ganymede, Europa-Callisto, Ganymede-Callisto) is linkable somewhere in
   2-6 km/s even with zero padding. This is not a tolerance artifact.

The four Galilean moons are mutually Tisserand-linkable at the ~3-4 km/s energy that real tour
designers (Galileo, JUICE, and the literature `#858` Sec. 5.2 already flags as saturated —
Russell-Strange 2009, Liang-Hernandez, Campagnola) have long exploited precisely because this
system IS a near-complete connectivity graph at achievable flyby energies. A sequence built only
from 4 mutually-linkable moons is, by construction, never going to fail a *pairwise* Tisserand
screen — the alphabet's own physics guarantees every leg passes. This confirms and sharpens `#858`
Sec. 5.1's "quiet flaw" warning (multiplying sequences at the same thin per-sequence density) from
a different angle: not only is per-sequence cell density thin, the sequence-level filter meant to
cut the count *before* that stage cannot discriminate at all for this specific system.

## Disposition, per `#858`/`#791`'s own decision rule

**Do not proceed to Stage 1 as scoped.** Running the full patched-conic enumeration
(`campaign_runner`, 4-32k cells × 3269 unpruned sequences ≈ 10^7 cells) would be exactly the
combinatorial explosion `#858` Sec. 5.1 warned against, with Stage 0 contributing no narrowing at
all. Options for whoever revisits this (none attempted here, all out of scope for this checkpoint):

- **The chained/graph-level energy-consistency check was tried same-day — it also fails** (see
  Addendum below): the alphabet has a universal shared-V∞ "hub" energy (4-5 km/s) at which all 12
  possible moon pairs are simultaneously linkable, so even the strictly-stronger single-V∞
  constraint cannot discriminate. Energy-based sequence filtering appears structurally unable to
  narrow this alphabet.
- **A non-energy criterion** — e.g. phasing/resonance compatibility (does the sequence's own
  period_k admit a closed synodic pattern across the moons' actual orbital periods, not just "is a
  transfer energetically possible") or ToF-window intersection across the whole chain — is the
  next candidate the Addendum identifies, not attempted here.
- **Accept Stage 1's full cost** — spend the ~300-2800 CPU-hours `#858` Sec. 5.3 already costed,
  with the explicit understanding that this is now the ENTIRE proposed cost, not a
  pruning-narrowed fraction of it.
- **Shelve**, on the same logic `#861` used for Resonant Atlas: two decisive negatives on the
  intended cheap-first gate (pairwise AND chained energy filters), redirect effort elsewhere
  (`#790`, per `#858`'s own ranking).

No recommendation is made here between these — that is a dispatch decision, not a Stage 0 finding.
**Explicitly NOT run**: any cell evaluation, `campaign_runner` dispatch, or n-body shooting. No
`data/catalogue.yaml` writeback (methodology/pruning checkpoint only).

## Addendum: the chained filter

The pairwise screen above lets each leg use a *different* V∞ from the probe band — a weaker
constraint than a real ballistic tour actually obeys. A flyby rotates V∞ direction but preserves
its *magnitude* (`releg_moontour.py`'s own "Continuity model" docstring already encodes this for
the powered lane), so an unpowered patched-conic tour should hold V∞ quasi-constant across the
whole cycle. The natural next, strictly-stronger filter: does a **single common V∞** in the same
probe band link **every** leg of the sequence simultaneously?

Rebuilt the driver (`scripts/scan_791_stage0_sequence_pruning.py`, same file) around a
precomputed per-pair-per-V∞ linkability table (12 ordered pairs × 7 probe V∞ = 84 `linkable()`
calls total, cached once) so the chained check is `intersection of each leg's linkable-V∞ set`.
Positive control (IEGI) still survives (common V∞ = {4, 5} km/s).

**Result: still 3269/3269 new sequences survive. Zero additional pruning from the stronger
constraint.**

**Why: the alphabet has a universal "hub" energy, not just per-edge connectivity.** Checked
directly whether any single V∞ links *every one* of the 12 ordered pairs simultaneously:

| V∞ (km/s) | all 12 pairs linkable? |
|---|---|
| 4.0 | **True** |
| 5.0 | **True** |
| 6.0 | False (2/12 pairs drop) |
| 8.0 | False (7/12 pairs drop) |
| 10.0+ | False (all drop) |

At V∞ = 4-5 km/s the four-moon system is **fully connected simultaneously** — not just
edge-by-edge at possibly-different energies, but at one shared energy covering the whole alphabet
at once. Any sequence built only from these 4 moons can therefore always fall back to V∞=4 or 5
km/s for every leg, so the chained constraint — despite being physically the *correct*, stricter
requirement — collapses to the same non-discriminating result as the pairwise one. This is a
sharper version of the original finding, not a weaker one: it rules out "the pairwise screen was
just too permissive" as the explanation and confirms the real cause is structural — the Galilean
system's own physical connectivity at achievable flyby energies, the same fact that makes this
territory so literature-saturated (`#858` Sec. 5.2).

**Implication:** any Tisserand/energy-magnitude-based sequence-level filter is unlikely to
discriminate meaningfully in this specific system, because the discriminating power such filters
rely on (a sparse, mostly-disconnected reachability graph) is simply absent here — Jupiter's inner
three-plus-one Galilean moons are close enough in semi-major axis and share a common natural flyby
energy that connects all of them at once. A filter that *would* discriminate would need to test
something ENERGY doesn't capture — e.g. phasing/resonance compatibility (do the moons' orbital
periods admit a closed synodic pattern for this specific sequence and period_k, not just "is a
transfer energetically possible"), or ToF-window intersection across the whole chain. Neither is
attempted here — genuinely out of scope for a same-day Stage 0 follow-up, and a different kind of
check than "Tisserand/energy-graph pruning" as specified.

## Verification

`ruff check` / `ruff format --check` clean; full `uv run mypy src tests` clean (`scripts/` is not
a mypy root per `pyproject.toml`'s `[tool.mypy]` config, matching every other driver script in this
directory). No dedicated pytest file — this is a thin driver script over already-tested library
calls (`tisserand.linkable`, `moon_prune._a_range_au`), per this project's established convention
for driver/orchestration scripts. Data: `data/found/791_stage0_sequence_pruning/results.jsonl`
(3276 rows, one per enumerated sequence, independently spot-checked above against a from-scratch
positive/negative-control re-derivation, not trusted from the driver's own summary print alone).
