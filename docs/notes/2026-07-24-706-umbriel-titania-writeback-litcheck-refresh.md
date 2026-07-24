# #706 — Uranus Umbriel-Titania literature clearance refresh (pre-writeback gate)

**Date**: 2026-07-24 (same day as `#699`). **Scope**: re-verification only, reusing `#699`'s own
methodology, refreshed for the sharper claim `#701`/`#704`/`#705` have since established (a
real-ephemeris-survivable, epoch-recurring homoclinic connection across 10 epochs spanning
2000-2083 — not just "does a CCR4BP torus/connection exist on this pair in the literature at
all," which is all `#699` checked for).

## What was re-checked

**1. `#699`'s own note** (`docs/notes/2026-07-24-699-umbriel-titania-deep-litcheck.md`), read in
full. Verdict was CLEAR after 10 live queries + full anchor re-read + corpus scan + full-text
check of the closest adjacent paper (Kumar arXiv:2509.03655).

**2. This project's own knowledge, re-verified unchanged:**
- `src/cyclerfinder/search/literature_check.py` `KNOWN_CORPUS`: same six Uranian anchors as
  `#699` found (Heaton-Longuski 2003, Sims 2014, Kumar Uranus-Oberon PCRTBP MMR `arXiv:2509.03655`,
  Canales-Howell-Fantino Titania-Oberon halo transfer `arXiv:2110.03683`, Jarmak QUEST, UOP
  Decadal). No new anchor added since `#699`; Kumar's `body_set` is still
  `frozenset({"Oberon", "Titania"})` — Umbriel absent from every anchor.
- `docs/notes/CORPUS_INDEX.md`: no Uranus/Umbriel/Titania/Oberon/CCR4BP entries beyond the
  already-known ones. Two *new-since-#699* corpus additions exist (Kumar-Anderson-de la
  Llave-Gunter 2021 Europa-Ganymede CCR4BP, arXiv:2109.14815; Kumar-Anderson-de la Llave 2023
  Ganymede secondary-resonance-overlap CCR4BP, arXiv:2309.06073) but both are explicitly Jovian
  (Europa/Ganymede), acquired 2026-07-23 for `#688`/`#686` — not Uranian, not relevant here.
- `/Users/bruce/dev/cyclers_pdf/papers/`: directory listing re-scanned; still zero
  Umbriel/Titania/Oberon-titled PDFs beyond the already-known Heaton-Longuski and Voyager-2 items.

**3. Fresh live web searches** (WebSearch confirmed available, used directly — no access
limitation). Re-ran a representative subset of `#699`'s direct-pair queries plus new queries
targeting the sharper claim:
- `"Umbriel" "Titania" CCR4BP restricted four-body problem homoclinic`
- `Umbriel Titania real ephemeris transfer trajectory design`
- `Umbriel Titania homoclinic connection epoch recurrence`
- `Uranian moon quasi-periodic torus manifold transfer feasibility ephemeris`
- `Bhanu Kumar Uranus Titania Umbriel 2026 new paper`
- `Kumar Anderson de la Llave Uranus CCR4BP 2026 arxiv`
- `Umbriel Titania 2:1 resonance mission spacecraft trajectory AAS AIAA 2026`
- `arxiv 2026 Uranus Umbriel Titania torus connection AAS conference paper`

None surfaced a CCR4BP/torus/homoclinic/manifold paper on Umbriel-Titania specifically, at any
epoch-recurrence level of specificity. New items surfaced, all checked and ruled out:
- **arXiv:2412.20326** ("Orbital maneuvers for a space probe around Titania") — WebFetched
  directly. Single-moon (Titania-only) collision-avoidance station-keeping study; no multi-moon
  problem, no homoclinic/heteroclinic content, no Umbriel. Not relevant.
- Ariel-Umbriel resonance-capture/migration papers (arXiv:2509.24631, 2403.17896/2403.17897) —
  same subfield `#699` already correctly identified as planetary-formation dynamics, not
  astrodynamics trajectory design, and not touching Titania as a transfer target.
- JWST composition paper (arXiv:2607.05600, "carbon oxides on Uranus's large moons") — planetary
  science, unrelated to dynamics/trajectory design.
- Bhanu Kumar's personal publication page (WebFetched again) — no new paper since `#699`'s check;
  most recent items are the already-known 2025 cislunar/Earth-Moon-system papers. Nothing on
  Umbriel, Titania, or Uranus CCR4BP beyond the already-anchored `arXiv:2509.03655`.

No conference-program or preprint evidence of any 2026 AAS/AIAA/journal paper specifically
targeting a real-ephemeris-survivable or epoch-recurring homoclinic connection at any Uranian
moon pair.

## Verdict: **STILL CLEAR**

Nothing new has surfaced since `#699`'s same-day clearance. The in-repo corpus and `KNOWN_CORPUS`
are unchanged with respect to Umbriel-Titania (the two new CCR4BP corpus additions since `#699`
are Jovian, not Uranian). Fresh live searches — both a repeat of `#699`'s direct-pair queries and
new queries specifically targeting the sharper real-ephemeris/epoch-recurrence claim — found no
disqualifying prior work. The one genuinely new item surfaced (arXiv:2412.20326, Titania-only
station-keeping) does not touch the Umbriel-Titania pair or any multi-moon/homoclinic content.
Safe to proceed to the schema/writeback decision from a novelty standpoint; this task makes no
recommendation on that decision itself (that's `#707`'s and the coordinating session's call).

## What was NOT done (by design, per task scope)

No schema/writeback recommendation. No catalogue or `OUTSTANDING.md` edits. No new
`KNOWN_CORPUS` anchor added (nothing new and substantive found to anchor).
