# #700 — Deep literature check: Jupiter Europa-Callisto CCR4BP candidate (2026-07-24)

Analysis-only (no code, no catalogue writes, no corrector/EOM changes). Dispatched in parallel
with `#699` (Uranus Umbriel-Titania) to give `#693`'s fourth-ranked CCR4BP candidate — Jupiter
Europa-Callisto, `mu_pert=5.67e-5`, 4.73(4) period ratio — its own dedicated, independent
literature check, since `#693`'s own pass only piggybacked on an adjacent Ganymede-Callisto query
and explicitly flagged this pair as under-checked.

## Part 1 — Novelty

**Existing project knowledge, checked first:**
- `src/cyclerfinder/search/literature_check.py`'s `KNOWN_CORPUS` (grepped case-insensitively for
  Europa/Callisto/Ganymede/Jovian/Galilean): no CCR4BP/whiskered-torus/heteroclinic anchor for
  ANY Jovian pair currently registered — not Europa-Ganymede, not Ganymede-Callisto, not
  Europa-Callisto. (The two Kumar-group CCR4BP papers are digested in `docs/notes/` but were never
  added to `KNOWN_CORPUS` — confirmed by grepping for their arXiv IDs / AAS numbers, zero hits.)
  Existing Jovian anchors are all classical ballistic-cycler / MGA-tour papers (Liang 2024
  Callisto-Ganymede-Europa triple cyclers, Hernandez/Jones/Jesick Io-Europa-Ganymede triple
  cyclers, Russell-Strange 2009, Niehoff 1970, Campagnola 2014) — a different subfield
  (patched-conic ballistic cyclers/tours), not CCR4BP tori.
- `docs/notes/CORPUS_INDEX.md` + the two Kumar CCR4BP digest notes
  (`2026-07-23-digest-kumar-2021-europa-ganymede-ccr4bp-resonant-orbits.md`,
  `2026-07-23-digest-kumar-2023-secondary-resonance-overlap-ccr4bp.md`): scope is Europa-Ganymede
  only (base+perturber) and the Europa-induced secondary resonances inside the Ganymede 4:3 family.
  **Independently re-verified past the digest** by `pdftotext`-extracting both source PDFs
  (`kumar-anderson-delallave-gunter-2021-...-arxiv-2109.14815.pdf`,
  `kumar-anderson-delallave-2023-...-arxiv-2309.06073.pdf`) and grepping for "Callisto":
  **zero matches in either paper**, in the full text layer, not just the abstract. Neither paper
  discusses Callisto even in passing (no future-work mention either).
- `ls /Users/bruce/dev/cyclers_pdf/papers/`: no Aryan/Fitzgerald/"AAS 24-103" file present — this
  paper (found by `#693`'s live search, flagged as "not yet acquired") is still **not acquired**
  in this project's corpus. No other undigested Jovian-multi-moon-CCR4BP paper found in the
  directory listing.

**Live search (WebSearch available, five independent targeted queries run):**
1. `"Europa" "Callisto" CCR4BP restricted four-body torus` — only surfaced the known
   Europa-Ganymede (2109.14815) and secondary-resonance (2309.06073) papers, plus unrelated hits.
2. `Aryan Fitzgerald "AAS 24-103"` — confirmed full title/authors: **Suryansh Aryan (Virginia
   Tech, Space@VT GRA) & Riley M. Fitzgerald, "Four Body Invariant Structures And Chaos Analysis
   for Jovian Multi-Moon Ballistic Transfers," AAS 24-103 (2024)**.
3. `Jupiter Europa Callisto heteroclinic connection quasi-periodic invariant torus`
4. `"Jovian Multi-Moon Ballistic Transfers" ... "Europa-Ganymede" "Ganymede-Callisto" system
   PCCFBP abstract`
5. `Kumar Anderson "de la Llave" Gunter Jupiter Callisto CCR4BP torus` (checks whether the
   *originating* CCR4BP group has any other paper covering Callisto — no hit; their published
   work stops at Europa-Ganymede + the Ganymede-4:3 secondary-resonance extension)
6. `non-adjacent moon pair CCR4BP OR PCCFBP "skip" intermediate moon` — no hits describing a
   studied "skip-one-moon" class at all.

**Key finding on the Aryan & Fitzgerald paper (the one genuinely live question):** full text is
paywalled (ResearchGate 403; not on arXiv, NTRS, or any university repository found). This is an
honest, stated gap — I could not read it page-by-page. However, five independently-phrased
WebSearch queries converged, consistently and without contradiction, on the same description of
its scope: it builds quasi-periodic invariant-torus families for **two adjacent-pair PCCFBP
systems** — "Jupiter-Europa-Ganymede-spacecraft" and "Jupiter-Callisto-Ganymede-spacecraft" — each
with **Ganymede as the shared third body/perturber**, and assesses "the existence of invariant
manifold transit connections." The reported "Europa-Callisto transfer feasibility" result is
explicitly a **Ganymede-mediated two-hop chain**: an L1 unstable manifold from Callisto reaching
Ganymede in ~60 days, then an L2 stable manifold from Ganymede reaching Europa in ~74 days — i.e.
a manifold-hopping TOUR through Ganymede as an intermediate torus, not a directly-modeled
Europa-Callisto CCR4BP system that skips Ganymede entirely (which would require Europa as base,
Callisto as sole perturber — a different `mu`/`mu_pert` pair from either of the two systems this
paper actually builds). This directly answers the task's core concern: the paper's silence on a
literal Europa-Callisto PAIR (as opposed to a Europa↔Ganymede↔Callisto tour) is real, not just an
artifact of thin checking.

**Broader finding, as requested:** no evidence surfaced, in this pass or `#693`'s, of ANY published
CCR4BP/PCCFBP paper studying a genuinely non-adjacent moon pair (skipping an intermediate moon)
for any planet. Every hit found — Kumar's Europa-Ganymede, Aryan-Fitzgerald's
Europa-Ganymede/Callisto-Ganymede, the Uranus Titania-Oberon/Umbriel-Titania work `#699` is
independently checking — pairs a base moon with its immediate neighbor(s), or chains adjacent
pairs through a shared hub. This project has now looked at this question twice (`#693`, here) and
found nothing suggesting "skip-one-moon" CCR4BP is a studied class at all — plausibly because the
model's own weak-perturbation premise (comparable-order forcing) degrades badly once you skip a
comparably-massed intermediate body, though that's a physical-plausibility inference, not a
literature-sourced claim.

**Novelty verdict: CLEAR.** No genuine prior CCR4BP/whiskered-torus/heteroclinic work found for
the direct Europa-Callisto pair, across in-repo corpus, KNOWN_CORPUS, and six independent live
searches. Necessary-not-sufficient, per this project's own standing discipline — but this pair has
now been checked far more thoroughly than `#693`'s single piggybacked query.

## Part 2 — Tractability (independent of novelty)

`#693` flagged Europa-Callisto's 4.73 period ratio as the loosest of the tractable Jovian
candidates, with "no clean low-integer commensurability." Ran a quick, disposable numerical check
(not committed; reuses `cyclerfinder.core.cr3bp` production code plus the exact same test-only
`_resonant_symmetric_orbit` Newton-shooting scaffold used in
`tests/search/test_ccr4bp_torus_io_ganymede.py`, and `genome.composed_moon_map.resonance_semimajor`
for the seed semimajor axes) — script at
`/private/tmp/claude-501/-Users-bruce-dev-cyclers/e8a086b8-fae2-4e77-b340-1425b9d3c532/scratchpad/check_europa_callisto_resonance.py`,
not part of the repo.

**Constants** (from `core.satellites` JPL-SSD-sourced registry, cross-checked against `#693`'s own
table — identical): `mu` (Europa) `= 2.528018e-5`, `mu_pert` (Callisto) `= 5.666828e-5`,
`a_pert = 2.8054` (Callisto SMA in Europa units). Simple two-body Kepler check
(`a_pert^1.5 = 4.699`) matches `#693`'s real-period-sourced 4.734 to within the expected
real-vs-Keplerian discrepancy (moon-mass correction, same order as the `#690` `omega_gan` note).

**Best low-order rational approximants of 4.7** (continued-fraction convergents): 4/1 (off by
0.70), 5/1 (off by 0.30), then 14/3 (4.667) and 19/4 (4.75) as the next tier. Contrast with the
built/ranked reference systems, all much closer to their nearest integer: Io-Europa exactly 2.000,
JEG 2.03 (off by 0.03), Umbriel-Titania 2.10 (off by 0.10), Io-Ganymede 4.06 (off by 0.06).
Europa-Callisto's nearest-integer gap (0.30, even taking the closer of 4 or 5) is 3-10x larger
than any of these.

**Numeric result:** tested seed candidates at `(p_sc, q_moon) ∈ {(3,1),(4,1),(5,1),(1,4),(1,5),
(14,3),(19,4),(33,7),(52,11),(71,15)}`. Every interior candidate (radius 0.34-0.48 Europa-SMA
units) converged to residual `1e-13`-`1e-14` (tighter than tolerance) with essentially zero
eccentricity (near-circular family member — same character as Io-Ganymede's own chosen 4:1 orbit,
per that test's own docstring noting it is "much less eccentric" than JEG's literature-sourced
Europa 3:4 orbit). Clearance to Callisto's orbit (`a_pert=2.805`) is comfortable for every interior
choice: gap `2.32`-`2.46` Europa-SMA units — **larger** margin than Io-Ganymede's own chosen 4:1
orbit needed (that system's minimum spacecraft-Ganymede distance over a full torus period was
~2.14 Io-SMA units, per that test's docstring). Exactly one candidate — the naive exterior 1:5
choice (`a=2.924`) — comes out **beyond** Callisto's orbit (gap `= -0.12`), a direct structural
analog of the exact Io-Ganymede near-collision the task asked about (that system's naive 1:4
choice sat almost exactly at Ganymede's radius and threw a genuine integrator `RuntimeError`); here
it is simply avoidable by using an interior candidate instead, with no forced choice through the
danger zone.

**Tractability verdict — findable, not harder than precedent, but with an honest caveat:** a
clean, cheaply-converged, comfortably collision-clear base periodic orbit is straightforwardly
findable for Europa-Callisto — if anything with wider safety margins than Io-Ganymede's own build
needed. This is only a first-step check (base CR3BP orbit, not torus persistence under Callisto's
actual `mu_pert` forcing, which is what `#690`/`#696`'s full pipeline verifies) so it does not
prove a torus/heteroclinic search would succeed — but it directly answers the question asked: no
sign of a "much harder/messier" base-orbit-finding wall analogous to Io-Ganymede's real one.
**The genuine caveat is softer and different in kind**: because 4.7 has no clean low-integer
commensurability, every candidate `p:q` seed tested converges to essentially the *same* near-
degenerate near-circular family at whatever radius Kepler's third law assigns it — unlike JEG
(anchored to a literature-sourced eccentric orbit at a real 3:4 resonance) or Io-Europa (the
tightest resonance in the solar system), there is no genuine, physically-motivated mean-motion
resonance for a Europa-Callisto build to be "near." Any resulting torus/heteroclinic search here
would be exploring a more arbitrary corner of phase space, not a recognized dynamical structure —
a real difference in scientific motivation from the other ranked candidates, even though it is not
a numerical build blocker.

## Verdict

**CLEAR** — novelty-clean (no direct Europa-Callisto CCR4BP/PCCFBP work found after a proper
independent search, including page-by-page-verified silence in both existing Kumar-group papers
and five further live queries) **and** a findable, cleanly-converging, collision-clear base orbit
(quick numeric check shows more comfortable margins than the Io-Ganymede precedent needed). Ready
to promote to a build task analogous to `#695`/`#696`, **with one flagged caveat for the
dispatching decision**: unlike every other ranked candidate, Europa-Callisto has no clean
low-integer MMR to anchor a dynamically-meaningful seed-orbit choice — a build here would be a
generic parameter-scan exercise rather than a resonance-motivated discovery search, which may
argue for ranking it below Io-Europa/Io-Ganymede/Umbriel-Titania on scientific-interest grounds
even though nothing here blocks it on tractability or novelty grounds.

Also flagged, unresolved from `#693`, still open: the Aryan & Fitzgerald AAS 24-103 (2024) paper
remains **unacquired** in this project's corpus (paywalled, no arXiv/NTRS/university-repository
copy found) — its Europa-Ganymede and Ganymede-Callisto PCCFBP results should still be acquired
and digested before any Jovian CCR4BP build lands, per `#693`'s own recommendation, independent of
this note's Europa-Callisto-specific verdict. No new paper was acquired this session, so nothing
new is registered in `CORPUS_INDEX.md` or `KNOWN_CORPUS`.
