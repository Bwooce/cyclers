# #699 — Uranus Umbriel-Titania deep literature check (CCR4BP novelty)

**Date**: 2026-07-24
**Scope**: research-only, no code/catalogue changes. Single question: has anyone already
published CCR4BP/whiskered-torus/quasi-periodic-torus/heteroclinic-transfer (or any other
multi-moon 4-body dynamical-systems astrodynamics) work on the Uranus Umbriel-Titania pair
specifically — the pair `#693` flagged as the best-conditioned non-Jovian CCR4BP candidate
(`mu_pert=3.92e-5`, Δi≈0, both e small, 2.10≈2:1 period ratio) but only lightly novelty-checked
(one query).

## What was checked

**1. This project's own existing knowledge.**
- `src/cyclerfinder/search/literature_check.py`'s full `KNOWN_CORPUS`, grepped case-insensitively
  for Uranus/Umbriel/Titania/Oberon/Ariel/Miranda. Six Uranian anchors exist (added under `#328`,
  lines ~1714-1850): Heaton-Longuski 2003 tour, Sims 2014 tour, Kumar Uranus-Oberon PCRTBP MMR
  (arXiv:2509.03655), Canales-Howell-Fantino Titania-Oberon halo transfer (arXiv:2110.03683),
  Jarmak QUEST, UOP Decadal. **Read every one of these in full** (not grep snippets). Critically:
  the Kumar anchor's `body_set` is `frozenset({"Oberon", "Titania"})` and the Canales-Howell
  anchor's is also `frozenset({"Titania", "Oberon"})` — **neither anchor's `body_set` includes
  Umbriel**. No existing anchor covers Umbriel at all, let alone Umbriel-Titania.
- `docs/notes/CORPUS_INDEX.md`, grepped for CCR4BP/four-body/torus/heteroclinic/Uranus/Kumar —
  no entry beyond the already-known Jovian Kumar/Aryan-Fitzgerald-adjacent digests and the six
  Uranian tour/MMR papers above.
- `docs/notes/2026-06-16-328-uranian-cycler-lit-deep-dive.md` read in full (20-query trail, the
  source digest behind the existing Kumar/Canales-Howell anchors). It explicitly classifies Kumar
  2509.03655 as "Section 6.2: Uranus-Oberon PCRTBP... extends to Uranus-Titania-Oberon CCR4BP
  secondary resonances" — **Titania appears only as the perturbing 4th body atop an Oberon-based
  resonant-orbit family**, never as the base moon paired with Umbriel. No Umbriel mention anywhere
  in this 163-line digest.
- `ls -la /Users/bruce/dev/cyclers_pdf/papers/` (218 files) — grepped for "umbriel", "titania",
  "oberon" in filenames: **zero matches for all three**. The Kumar arXiv:2509.03655 PDF itself is
  NOT in the corpus (only cited by digest/abstract, confirming `#693`'s flagged `inherited-
  unverified` provenance caveat on that anchor). No Uranian moon-pair PDF has ever been acquired
  beyond the Voyager 2 Uranus/Neptune encounter overview and the two already-anchored tour papers.

**2. Independent live search** (WebSearch + WebFetch both confirmed available and used — no
access limitation to report). Ten distinct targeted queries run, not one generic query:
1. `"Umbriel" "Titania" CCR4BP restricted four-body problem`
2. `"Umbriel" "Titania" quasi-periodic torus invariant manifold heteroclinic`
3. `Kumar Anderson de la Llave Gunter Uranus concentric circular restricted four body problem`
4. `Uranus moon system four body problem trajectory design mission design 2025 2026`
5. `Kumar arXiv 2509.03655 multi-shooting parameterization ... celestial resonant dynamics`
6. `"Umbriel" resonant orbit spacecraft trajectory AAS OR AIAA conference paper`
7. `Aryan Fitzgerald AAS 24-103 four body invariant structures chaos Jovian multi-moon ballistic transfers Uranus`
8. `Uranian satellites resonant periodic orbit survey Miranda Ariel Umbriel Titania Oberon dynamical systems 2024 2025`
9. `Umbriel Titania 2:1 near resonance mission trajectory design spacecraft`
10. `Bhanu Kumar Uranus moons publications list resonant dynamics 2026` → led to his personal
    publication page, WebFetched directly for his complete list.

**None of the ten surfaced a CCR4BP/torus/heteroclinic/manifold-transfer paper on Umbriel-Titania
specifically.** Queries 1-2 and 6, 9 (the most direct) returned zero on-point hits — only generic
RFBP methodology papers and the already-known Jovian Kumar/Aryan-Fitzgerald hits. Query 8 surfaced
only the active Ariel-Umbriel resonance-capture/migration subfield (arXiv:2509.24631, 2403.17896/
2403.17897, 2309.04786) — planetary-formation dynamics, not astrodynamics trajectory design,
exactly the distinction `#693` already correctly drew and which is reconfirmed here.

**3. Read-through, not abstract-skimming, of the two closest candidates:**
- **Kumar arXiv:2509.03655** ("Multi-shooting parameterization methods...", published in J.
  Nonlinear Science 2026, DOI 10.1007/s00332-026-10276-6): WebFetched the arXiv abstract page, the
  arXiv HTML full text (twice, targeting Section 6), and confirmed via live search of the
  predecessor AAS 2024 conference paper ("A Survey of Oberon Mean Motion Resonant Unstable Orbit
  Properties and Connections for Uranian Tours", Kumar & Anderson, AAS/AIAA 2024) that **the
  4-body extension is explicitly framed as "Oberon resonant orbits in a restricted 4-body model
  including Titania"** — Titania is the perturber on an Oberon-base-orbit family, structurally the
  opposite pairing/role from a putative Umbriel(base)-Titania(perturber) torus. Two independent
  WebFetch passes over the arXiv HTML confirmed **"Umbriel" does not appear anywhere** in the
  accessible text of this paper.
- **Aryan & Fitzgerald, AAS 24-103** ("Four Body Invariant Structures and Chaos Analysis for
  Jovian Multi-Moon Ballistic Transfers", 2024) — re-confirmed via live search (not re-verified in
  depth, out of this task's scope) that this paper's PCCFBP tori are computed for
  **Jupiter-Europa-Ganymede and Jupiter-Callisto-Ganymede only** — confirmed Jovian-only, no
  Uranian content, correctly out of scope for `#699` (belongs to the parallel `#700` Europa-
  Callisto check instead).
- **Canales-Howell-Fantino** (arXiv:2110.03683 + companion arXiv:2308.10029): already fully
  triaged by `#328`'s digest as Titania-Oberon halo-to-halo one-shot manifold transfer (not CCR4BP
  torus, not Umbriel). Re-confirmed scope by re-reading the existing anchor text; no new search
  needed since the pairing (Titania-Oberon) and moons involved are unambiguous and don't touch
  Umbriel.

## Verdict: **CLEAR**

No disqualifying prior CCR4BP-class (or any other multi-moon 4-body dynamical-systems
astrodynamics) work was found on the Umbriel-Titania pair specifically, after:
- a full re-read of every existing in-repo anchor and digest touching the Uranian system,
- a corpus-directory scan confirming zero Umbriel/Titania/Oberon-titled PDFs are held,
- ten independent, distinctly-worded live-search queries spanning direct-pair, author-driven,
  and general-survey angles,
- author-list verification (Kumar's complete publication list, fetched directly) confirming no
  paper of his covers Umbriel in any capacity,
- a full-text check (not abstract-only) of the single closest-adjacent paper (Kumar
  arXiv:2509.03655 / AAS 24 predecessor), confirming Titania appears there only as Oberon's
  perturber, never paired with Umbriel, and "Umbriel" is absent from the paper entirely.

The Uranian neighborhood's "active adjacent work" that `#693` flagged as scoop-context risk
resolves cleanly on inspection: it is either (a) Ariel-Umbriel resonance-capture/migration — a
different subfield (planetary formation, not astrodynamics) studying a different pair — or (b)
Titania-Oberon CCR4BP/halo-transfer — a different, already-published, already-anchored pair. Both
are genuinely adjacent, neither touches Umbriel-Titania. This candidate is ready to be considered
for promotion to a build task analogous to `#695`/`#696`, pending the coordinating session's own
GO/NO-GO (not this task's call).

## What was NOT done (by design, per task scope)

No CCR4BP build, no code, no catalogue writes. No action taken on Jupiter Europa-Callisto (`#700`
covers that independently). No new `KNOWN_CORPUS` anchor was added — nothing new and substantive
was acquired (no new PDF, no new paper found); adding a low-confidence anchor "just to have added
something" was explicitly against this task's own instructions, and every paper surfaced here was
already covered by an existing anchor or is out of scope (Ariel-Umbriel: different subfield, not
disqualifying; Jovian Aryan-Fitzgerald: belongs to `#700`, not `#699`). `docs/notes/CORPUS_INDEX.md`
left unchanged for the same reason.

## Sources consulted (live)

- arxiv.org/abs/2509.03655, arxiv.org/html/2509.03655(v1) — Kumar multi-shooting paper (2 fetches)
- link.springer.com/article/10.1007/s00332-026-10276-6 — paywalled (redirect to IdP login), not
  accessible; abstract content instead recovered via WebSearch snippet (6:5 Oberon MMR secondary
  resonance ratios matched `#328`'s digest numbers exactly, cross-confirming both sources)
- bhanukumar314.github.io — Kumar's personal publication list (direct fetch)
- researchgate.net/publication/383155359 — Aryan & Fitzgerald AAS 24-103 (Jovian-only, confirmed)
- Multiple arXiv/A&A/ScienceDirect hits on Ariel-Umbriel resonance-capture/migration (2509.24631,
  2403.17896, 2403.17897, 2309.04786) — confirmed different subfield, not astrodynamics
