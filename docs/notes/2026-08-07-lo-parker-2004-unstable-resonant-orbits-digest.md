# Digest: Lo & Parker 2004, "Unstable Resonant Orbits near Earth and Their Applications in Planetary Missions" (AIAA 2004-5304)

**Source:** M. W. Lo & J. S. Parker, AIAA/AAS Astrodynamics Specialist Conference, 16-19 August 2004,
Providence, RI. DOI `10.2514/6.2004-5304` (CrossRef-confirmed). Filed at
`cyclers_pdf/papers/lo-parker-2004-unstable-resonant-orbits-near-earth-planetary-missions-AIAA-2004-5304.pdf`
(md5 `da6336f61e1652eb1b6f0a4d0eeed6df`), text-layer PDF (`activePDF Toolkit`), 29 pages.

**Acquisition context:** uploaded directly by the user 2026-08-07 after this session identified it as
`#730`'s backlog item 30. **Important correction to the backlog framing** (see "What this paper is NOT"
below): this is a genuine, on-point acquisition for `#780` (Earth-Moon Casoliva lane), but it is **NOT**
the paper needed to unblock `#774`'s abandoned Saturn-Titan chain closure — that is a different Lo &
Parker paper, still unacquired (see below).

## What this paper is

A planar Earth-Moon PCRTBP study cataloguing **six families (A-F) of simple periodic symmetric orbits**
(orbits piercing the rotating-frame x-axis exactly twice per orbit, each crossing orthogonal to the
axis), built via Howell's differential corrector (Ref. 16) plus a standard pseudo-arclength continuation
method (perturb one IC component, differentially correct, predict-and-step along the family curve).

**Table 1** reproduces Broucke's (Ref. 8) six-class classification scheme for these orbits, based on
which side of which primary each of the two x-axis crossings falls on (Class 1 centred at L3, Class 2 at
M1, Class 3 at M1+M2 barycenter region, etc. — full table in the PDF).

**Table 2**, the load-bearing cross-reference for `#780`:

| Family | Broucke Class | Description |
|---|---|---|
| A | Class 6 | Lyapunov orbits about LL2 |
| B | Class 4 | Lyapunov orbits about LL1 |
| C | Class 5 | Distant retrograde orbits about the Moon |
| D | Class 5 | Low prograde orbits about the Moon |
| E | Class 5 | Distant prograde orbits about the Moon |
| **F** | **Class 3** | **Periodic resonant lunar flyby orbits** |

**Family F is the direct classification-scheme predecessor of Casoliva's own Class 1** (resonant p-q
lunar-flyby cyclers) — confirms and sharpens the backlog note's own framing ("the specific classification
of planar symmetric periodic-orbit families Casoliva's own introduction cites as its own predecessor").
Family F is further subdivided into sub-branches F1/F2/F3 (distinguished by period, ~30-90 days per
Figure 7's own axis range) — worth checking against Casoliva's own p-q labelling scheme when `#780`
builds its own family module.

**No digit-grade numeric IC/period/Jacobi table is printed for Family F** (or any family) — all family
data in this paper is presented graphically (Figures 4-9, initial-condition and period curves vs. x0),
not tabulated to precision. **This paper is a classification/method reference, not a source of sourced
numeric goldens** — `#780`'s own golden numeric anchor remains Casoliva 2008/2010's own tables, exactly
as already scoped. Do not expect this paper to substitute for that.

Sections IV-VI (not deep-read this pass, flagged for later reference if `#780` needs the corrector
methodology in more detail) cover mission-design applications of these families' invariant manifolds,
focused on resonant lunar-flyby trajectories — likely useful background for `#780`'s own Class 1 work,
though Casoliva's own papers should remain the primary implementation reference.

## What this paper is NOT — a citation correction

This project's own `#773`/`#775` results notes (`docs/notes/2026-08-01-773-resonant-chain-periodicity-closure.md:203-206`,
`docs/notes/2026-08-01-775-resonant-chain-continuation-closure.md:244-249,272-273`) flag "Lo & Parker's
own paper... ref. [68] in [Vaquero's] own thesis" as the untried escape hatch for `#774`'s abandoned
Saturn-Titan resonant-chain closure (an "iterative multi-patchpoint refinement" methodology). **Checked
directly against this paper's own text: no multi-patchpoint, patch-point, or multiple-shooting
methodology appears anywhere in it** — its corrector is a standard single-shooting Howell corrector, and
its continuation method is standard pseudo-arclength family continuation. This paper cannot be the `#774`
escape hatch.

**Checked directly against Vaquero's own bibliography** (`vaquero-2013-...-purdue-phd.txt:6618-6620`):
her own reference [68] is a **different** Lo & Parker paper:

> [68] M. W. Lo and J. S. Parker. **Chaining Simple Periodic Three-Body Orbits.** In AIAA/AAS
> Astrodynamics Specialist Conference, Lake Tahoe, California, August 7-11 2005. Paper **AAS-05-380**.

("Chaining" is also the operative word in Vaquero's own text at both cited locations — line 4248: "their
generating orbits [68]. These chains are similar to the homoclinic cycles previously discussed" — this
is unambiguously about the 2005 Lake Tahoe chaining paper, not the 2004 Providence paper filed here.) A
closely related possible backup is Vaquero's own ref [69]: J. S. Parker, K. E. Davis & G. H. Born,
"Chaining Periodic Three-Body Orbits in the Earth-Moon System," vol. 67, pp. 623-638, 2010 (a journal
version, possibly DOI-findable — not checked this pass).

**`#730`'s own backlog item 30 entry was correct all along** — it explicitly scoped this 2004 paper's
relevance to Casoliva's own citation (`#780`-relevant), not to `#774`. The conflation happened later, in
a cross-session survey that (correctly) surfaced ref [68]'s relevance to `#774` but then (incorrectly)
attached it to this already-identified-but-different backlog item's DOI. **`#774`'s own actual escape
hatch, Lo & Parker 2005 "Chaining Simple Periodic Three-Body Orbits" (AAS-05-380, Lake Tahoe), remains
unacquired** — a genuinely new, distinct acquisition target, not yet in the `#730` backlog under its own
entry. No DOI found for it in this pass (AAS conference papers frequently lack one, per this project's
own established pattern for pre-2010 AAS papers); worth a dedicated acquisition attempt if `#774` is ever
revisited.

## Citation-mining pass (background/related-work section)

References worth flagging for the `#730` backlog (checked against the current backlog list; not
re-flagging anything already present):

- **Ref. 8**: Broucke, R. "Stability of Periodic Orbits in the Elliptic, Restricted Three-Body Problem."
  *AIAA Journal*, 7(6), 1969 — the origin of the six-class classification scheme used throughout. Distinct
  from the already-filed `broucke-1968-...jpl-tr-32-1168` (a different Broucke work, JPL TR not AIAA
  Journal). Worth checking if this 1969 AIAA paper is a genuine backlog gap.
- **Ref. 18**: Parker, J.S. & Chua, W.S. (unclear exact title from this citation context) — invariant
  manifold construction reference, likely already covered by existing manifold-construction machinery in
  this project; not flagged as high-priority.
- **Ref. 19/20**: Matukuma, T. and Strömgren, E. — the pre-Broucke classification schemes this paper
  explicitly supersedes; historical interest only, not flagged for acquisition.

## Registration

Filed in `cyclers_pdf` (commit pending in that repo, separate from this one). `CORPUS_INDEX.md` and
`#730`'s backlog master list item 30 status updated in the same session. No catalogue or code changes —
digest and citation-mining only. Directly informs `#780` (dispatched, in progress as of this digest);
no action needed on `#774` from this paper specifically.
