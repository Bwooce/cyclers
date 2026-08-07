# Digest: Parker, Davis & Born 2010, "Chaining periodic three-body orbits in the Earth-Moon system"

**Source:** J.S. Parker, K.E. Davis & G.H. Born, *Acta Astronautica* 67(5-6):623-638 (2010).
DOI `10.1016/j.actaastro.2010.04.003`. Filed at
`cyclers_pdf/papers/parker-davis-born-2010-chaining-periodic-three-body-orbits-earth-moon-actaastro-67-623-doi-10.1016-j.actaastro.2010.04.003.pdf`
(md5 `2c1da419dbaf8885517ae739f5901bb7`), text-layer PDF (Elsevier), 16 pages, tagged.

**Acquisition context:** uploaded directly by the user 2026-08-07, identified this session as
Vaquero 2013's own reference [69] — a journal-published companion to ref [68] (Lo & Parker 2005,
"Chaining Simple Periodic Three-Body Orbits," AAS-05-380, still unacquired, no DOI exists) — flagged
by `#773`/`#775` as the untried escape hatch for `#774`'s abandoned Saturn-Titan resonant-chain
closure.

## What this paper actually contains — and whether it's the right paper

**Yes, genuinely on-point.** This paper builds complex orbit chains in the planar Earth-Moon CRTBP
via a **multiple-shooting differential corrector applied to a sequence of "patchpoints"** — abstract
`"A multiple-shooting differential corrector is used to construct complex orbit chains and complex
periodic orbits."` This is the same general technique class `#773` already tried (it reused this
project's own `cr3bp_multiple_shooting.correct_multiple_shooting`), but this paper's own contribution
is a **specific, non-uniform patchpoint selection strategy** that `#773`'s own attempt did not use:

- Patchpoints are placed at **natural dynamical waypoints** — the orthogonal x-axis crossings of each
  periodic orbit involved (states labelled A, B, E, F in the paper's own worked example) PLUS the
  x-axis crossings closest to a theoretical heteroclinic/homoclinic connection between them (states
  C, D, G, H) — not at uniformly-spaced time samples along the trajectory.
- The corrector then adjusts a full itinerary sequence (e.g. `{...,A,B,A,B,C,D,E,F,E,F,...}`) as one
  simultaneous multi-point boundary-value problem, using the near-matching but not-quite-continuous
  states at the connection points (the paper's own Table 1 example: state C is only ~306 km/0.8 m/s
  from state A, state G only ~2058 km/28.9 m/s from state E) as a warm start close enough to converge
  without difficulty.
- **This is exactly the "adaptive patch-point insertion, not uniform resegmentation" ingredient
  `#775`'s own note flagged as missing** — `#773`'s own `n_segments=8->16` test found finer UNIFORM
  segmentation did not help; this paper's technique is qualitatively different (patchpoints chosen at
  natural orbit/connection geometry, not by uniform time subdivision).

**One important caveat, checked directly rather than assumed:** this paper's own Section 2.5
("Differential correction") does NOT itself derive the corrector's algorithm in full — it states
`"This differential corrector may be characterized as an iterative algorithm, where each iteration
involves two levels of differential correction. The algorithm is described in detail in the
literature [20,21,40]"` and defers to:
- **[20]** Howell 1984 (already in this project's corpus, `howell-1984-...halo-orbits`) — standard
  single-shooting halo corrector, not itself multiple-shooting.
- **[21]** K.C. Howell & H.J. Pernicka, "Numerical determination of lissajous trajectories in the
  restricted three-body problem," *Celestial Mechanics* 41 (1988) 107-124 — plausibly the actual
  multiple-shooting derivation source. **Not in this project's corpus; not previously flagged in the
  `#730` backlog. New acquisition candidate.**
- **[40]** R. Wilson, "Derivation of differential correctors used in GENESIS mission design,"
  Technical Report JPL IOM 312.I-03-002, Jet Propulsion Laboratory, 2003 — a well-known internal JPL
  memo in the astrodynamics community for its multiple-shooting corrector derivation (associated with
  the GENESIS mission). **Not in this project's corpus; not previously flagged. New acquisition
  candidate — likely the deepest technical source, but JPL internal memos are frequently NOT publicly
  available; may require a targeted search (course materials, personal faculty pages sometimes host
  copies) rather than a standard DOI/journal lookup.**

**Practical assessment for `#774`:** this paper alone likely provides ENOUGH detail to prototype the
patchpoint-selection strategy (the worked example gives the concrete state-selection and sequence-
construction logic) even without acquiring [21]/[40] — the missing piece was the SELECTION strategy,
not the underlying multiple-shooting math (which this project's own `cr3bp_multiple_shooting` module
already implements, reused successfully by `#773`). Acquiring Howell & Pernicka 1988 and/or the Wilson
2003 memo would sharpen the corrector's own convergence behavior but is not obviously required to
attempt the technique.

## Registration

Filed in `cyclers_pdf` (commit pending in that repo). `CORPUS_INDEX.md` and `#730`'s backlog master
list updated in the same session (new item, since this paper had no prior backlog row). Two new
acquisition candidates flagged (Howell & Pernicka 1988; Wilson 2003 JPL IOM) but not acquired this
pass. **No decision made here on whether to reopen `#774`** — that closure was explicitly confirmed
abandoned-as-scoped by the user; reopening it with this new technique is a call for the coordinating
session/user, not something this filing task decides unilaterally.
