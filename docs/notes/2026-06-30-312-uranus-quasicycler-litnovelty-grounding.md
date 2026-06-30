# #312 — Uranus Umbriel-Oberon-Umbriel quasi-cycler: live-web lit-novelty grounding

**Date:** 2026-06-30. The #312 Part-B SILVER candidate (the session's highest-novelty find
of that round) had only an OFFLINE corpus check (`offline_corpus_search`, KNOWN_CORPUS
not-found, confidence 0.4). The mandatory baseline ([[feedback_literature_novelty_check_baseline]])
requires clearing it against the PUBLISHED record, not just the corpus subset. This grounds
it with a live web search. **Outcome: the lit-novelty necessary condition HOLDS at the
published-record level (strengthened), with one framing correction.**

## The candidate (re-verified from data)
`data/scan_312_uranus_oberon_umbriel.jsonl`: **Umbriel-Oberon-Umbriel (1,1)**, residual
**0.0252 km/s** (sub the 0.05 gate), V∞ per encounter **[0.92, 0.96, 0.895] km/s** (low
excess speed), verdict SILVER. A near-ballistic PERIODIC repeated-moon cycle — `quasi_cycler`
class (passes_v2 by bounded drift-oscillation, like #339/#344, not strict closure).

## Live web search (3 queries, 2026-06-30)
1. **No published Uranian cycler / quasi-cycler surfaced.** Queries on Uranus moon cyclers
   + repeated periodic flybys returned only Uranian-system *resonance-evolution* dynamics
   (Ariel-Umbriel 5:3 / 2:1 MMR history) — and the explicit note that the Uranian regular
   moons are NOT currently in any mean-motion resonance. No spacecraft repeated-flyby cycler.
2. **The moon-CYCLER literature is Jupiter + Saturn ONLY.** Russell-Strange "Planetary Moon
   Cycler Trajectories" (2009) and Lynam-Longuski "Laplace-resonant triple-cyclers" (2011)
   — the foundational moon-cycler papers — treat two-moon Galilean and Titan-Enceladus
   cyclers; the search "did not return Uranus/Neptune moon cyclers." (Both are corpus anchors.)
3. **The Uranian published work is one-way mga_tours, a DIFFERENT object.** Heaton-Longuski
   2003 (the "Galileo-style tour of the Uranian satellites", DOI 10.2514/2.3981 — ALREADY in
   our catalogue as `heaton-longuski-2003-uranian-tour-u00-01`, `mga_tour`) and the Uranus-
   orbiter tour papers use Umbriel/Oberon/Titania flybys, but: at **V∞ ~2-4 km/s** (catalogue
   entry: Titania 3.27-3.30, Oberon 2.13-2.98, Umbriel 3.64), in a **one-way epoch-locked
   pump-down tour** (40 flybys → Ariel rendezvous), NOT a periodic ballistic cycler. The
   orbiter-tour endgames do reach <1 km/s, but as a TARGET-MATCHING approach (matching Ariel's
   orbit), not a repeated Umbriel-Oberon cycle.

## Verdict
- **Lit-novelty necessary condition HOLDS (strengthened).** No published Uranian cycler or
  quasi-cycler exists at the live-web published-record level — corroborating the #312 offline
  not-found with a broader search. The candidate is a plausibly-fresh **repeated-moon
  quasi-CYCLER at Uranus**, a regime the published moon-cycler literature (Jupiter + Saturn)
  has not reached.
- **Framing CORRECTION to the #312 note.** Its "NO published Uranian prior exists … fresh
  primary" is now outdated: the Uranian *system* IS published (Heaton-Longuski mga_tour, now
  catalogued; multiple orbiter-tour studies). The precise novelty is therefore **"first
  repeated-moon quasi-cycler at Uranus" — a fresh class/topology/V∞-regime in a system whose
  one-way mga_tours are published**, NOT an untouched primary. Grounding the claim against the
  real published record narrows it honestly ([[feedback_ground_citations_against_content]]).
- **CORRECTION (verify recalled refs against current state):** this candidate is **already
  ADMITTED** — catalogue row `umbriel-oberon-1-1-uranian-quasi-cycler-2026`, the catalogue's
  first computed `quasi_cycler`, **`validation_level: V4`** (#339 admission + #340 V0→V4
  promotion; #332 V4-scipy + #335 V4-strict URA111 real-eph + #338 EFFECTIVELY_CYCLIC all
  COMPLETE). My initial "blocked on #332" framing (and the first draft of this note) read a
  STALE OUTSTANDING line from 2026-06-16 that predates #335/#339/#340 — corrected here
  ([[feedback_check_dont_guess]]: verify a recalled doc/flag against the live catalogue before
  asserting a blocker). The four V-tier frozen gates pass (`tests/verify/test_silver_327_v*`).
- **What today's live check actually adds:** the admitted row's lit basis was the offline
  41-anchor corpus + a Heaton-Longuski-2003 direct read. This adds a **broader live-web
  published-record confirmation** (no published Uranian cycler exists; moon-cycler literature
  is Jupiter+Saturn only) — strengthening, not gating, the existing V4 row. It also sharpens
  the novelty wording: "first repeated-moon quasi-cycler at Uranus" (fresh topology in a
  published mga_tour system), not "untouched fresh primary."

## Bankable
The catalogued V4 Uranian quasi-cycler's novelty survives a fresh live published-record
lit-check (necessary condition re-confirmed at a broader scope than the admission-time offline
corpus + Heaton-Longuski read). The "fresh primary" overstatement is corrected to "fresh
quasi-cycler topology in a published mga_tour system." No catalogue change (the row is already
V4); the row's lit basis is strengthened, not gated.
