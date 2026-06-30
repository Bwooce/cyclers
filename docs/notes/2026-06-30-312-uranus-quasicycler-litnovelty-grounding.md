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
- **Still necessary-not-sufficient.** Novelty + admission remain gated on the V4 Uranian-system
  real-eph gauntlet (#332, the outstanding infra blocker) and the human review. No catalogue
  writeback; the candidate stays flagged, now with a published-record-grounded novelty scope.

## Bankable
The #312 SILVER survives a live published-record lit-check as a fresh Uranian quasi-cycler
(necessary condition met); the "fresh primary" overstatement is corrected to "fresh
quasi-cycler topology in a published mga_tour system." The admission blocker is #332 (V4
Uranian HFEM), not novelty.
