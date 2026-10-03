# Digest — Howett et al. (2021), "Persephone: A Pluto-system Orbiter and Kuiper Belt Explorer"

**Digested:** 2026-10-03 (`#881`). Supplied by the owner as the publisher EPUB; filed at
`cyclers_pdf/papers/howett-2021-persephone-pluto-system-orbiter-kuiper-belt-explorer-psj-2-75-doi-10.3847-PSJ-abe6aa.epub`
(md5 `fe98c5a2069499be1744eb25205e8aa7`). Read in full from the EPUB's text (one XHTML body,
about 75,000 characters); the figures were not inspected.

**Citation (CrossRef-verified):** Howett, Robbins, Holler, Hendrix, Fielhauer, Perry, et al., *The Planetary Science Journal* **2**(2):75 (2021),
DOI 10.3847/PSJ/abe6aa; arXiv:2102.08282.

**Why it was read:** the corpus anchor for this paper had never been checked against the paper
(the 2026-06-15 Pluto review says "the full PDF was not opened"). The anchor named the first
author "Howard", gave article number 56, and carried DOI 10.3847/PSJ/abf837, which resolves to an
unrelated paper (Tollefson et al., Neptune VLA/ALMA brightness temperatures). After `#880` the
anchor's scope also had to be settled, because every Pluto-system candidate reads `published`
through it.

## What it is

A NASA concept-mission study: launch 2031 on SLS Block 2 with a Centaur kick stage, radioisotope
electric propulsion, a Jupiter gravity assist (minimum altitude 17.8 Jupiter radii), a 27.6-year
cruise with one Kuiper Belt object flyby, arrival at Pluto in 2058, a 3.1-year orbital campaign of
the Pluto system, and an optional extended mission to another Kuiper Belt object. Eleven
instruments; nominal cost $3.0 billion. Almost all of the paper is science case, payload and
spacecraft.

## The trajectory content (Sec. 4.2, "Pluto Orbit Phase") — the part that matters here

Quoted, because the anchor's scope rests on it:

- "the science orbit consists of multiple periodic orbits designed in the Pluto/Charon restricted
  three-body dynamics model."
- The orbits "span very large radial distances both in and out of the satellite plane, enabling
  close encounters with Pluto and all its satellites."
- "The periodic nature of these orbits in the Pluto/Charon rotating frame, coupled with the fact
  that Pluto and Charon are tidally locked, has repeating ground tracks naturally built into the
  trajectory."
- "The periodic orbits identified fell into two primary categories: a high out-of-plane component
  to enable high-latitude global mapping, and low altitude to enable in situ sampling. Four
  distinct periodic orbits—two from each category—... were selected (Figures 13 and 14)."
- "Persephone would have the capability to transfer between periodic orbits."
- Minor satellites: the tour is defined "by seeing when the closest approaches are to each
  satellite"; coverage is best for Styx and worst for Hydra; "The encounter velocity of the minor
  satellites would be <300 m s−1."
- Charon gravity assist appears only in Sec. 4.3, as a way to DEPART the Pluto system for the
  extended mission.

Figure 13 shows the four orbits in the Pluto/Charon rotating frame with the minor satellites'
orbits and the Lagrange points; Figure 14 shows ground-track coverage (orbits 1 and 2 global,
orbits 3 and 4 low-altitude and low-latitude).

## What it does NOT contain

- No initial conditions, periods, Jacobi constants, stability indices or family names for any of
  the four orbits. They exist only as figures.
- The words "cycler", "halo" and "CR3BP" do not appear. "Resonant" appears once, about Kuiper Belt
  populations.
- No gravity-assist tour among the small moons, and no patched-conic or Lambert construction.

## Consequences for this project

1. **The anchor's original description was right and is restored.** It read "(CR3BP periodic
   orbits)"; `#881` first renamed it "(mission concept study)" on the strength of the abstract
   alone. The paper is a concept study whose science phase is built from restricted three-body
   periodic orbits.
2. **Scope: deliberately left undeclared, now for a sourced reason.** The paper's orbits are both
   planar-ish low-altitude and strongly out-of-plane, are periodic in the binary rotating frame,
   and repeatedly encounter Pluto, Charon and (opportunistically) the small moons. No single
   standard topology label covers that, and declaring a narrower one would under-cover. So any
   Pluto-system periodic or repeating-encounter candidate should keep reading `published` against
   this anchor and be compared with its Figures 13-14 by a person. This is prior art at the level
   of a CLASS: there is nothing in it to reproduce numerically.
3. **`#492`'s third ground stands in substance.** The `#880` audit called the Persephone hit on
   the nine Pluto closures an artefact and withdrew that ground. The way the hit was generated was
   mechanical, but the citation is relevant: repeated-encounter trajectories in the Pluto system
   are published, and in the binary model that `#492`'s own note says its point-primary closures
   lacked. The audit note and the `#492` addendum are corrected accordingly.
4. **The Pluto-Charon catalogue rows are unaffected.** They are already `known-class-member` of
   the Ross and Roberts-Tsoukkas class; this paper is a second, earlier class-level precedent for
   periodic three-body orbits used at Pluto-Charon and could be added to their
   `corroborating_sources` (not done here).
5. **Training labels.** `ml/falsepos_labels.py` has two synthetic records "patterned on" this
   paper. Their `_source` strings carried the wrong author and article number (fixed). Their
   numbers (for example encounter periods of 6.387 and 12.774 days, "halo, 2:1") are invented
   shapes, flagged `_mocked`; the paper gives no periods and never uses the word halo.
