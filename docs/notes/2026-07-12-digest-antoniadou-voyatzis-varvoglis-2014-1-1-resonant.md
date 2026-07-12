# Digest: Antoniadou, Voyatzis & Varvoglis 2014 — 1/1 resonant periodic orbits in 3D planetary systems

Single-paper digest. Read 2/2 pages (full paper — a short conference proceedings
note) on 2026-07-12 AET. User located and supplied this after I identified it
as a citation candidate via CrossRef while scoping task #579. Filed to the
private corpus as
`antoniadou-voyatzis-varvoglis-2014-1-1-resonant-periodic-orbits-3d-planetary-systems-iau-proc-310-82-doi-10.1017-S1743921314007893.pdf`.

## 1. Header

- **Title (verbatim)**: *1/1 resonant periodic orbits in three dimensional
  planetary systems*
- **Authors**: Kyriaki I. Antoniadou, George Voyatzis, Harry Varvoglis —
  Aristotle University of Thessaloniki
- **Venue**: *Complex Planetary Systems*, Proceedings IAU Symposium No. 310
  (2014), pp. 82-83 (a 2-page conference proceedings note, NOT a full journal
  article)
- **DOI**: 10.1017/S1743921314007893 (printed on p. 82, sourced from the page)

## 2. What the paper actually is — and the critical model-class caveat

**This is a short proceedings summary of the SAME general (non-restricted,
two-massive-body) three-body-problem model** used in
`docs/notes/2026-07-12-digest-antoniadou-voyatzis-2013-2-1-resonant-3d-gtbp.md`
— the paper's own text states it directly: "we utilize the spatial general
TBP in a rotating frame of reference (Antoniadou and Voyatzis (2014))." It
studies **1:1 mean-motion resonance (co-orbital motion)** between two
massive planets with non-zero mutual inclination — the spatial extension of
planar co-orbital families.

**This does NOT fix the `literature_check.py` anchor's claimed scope.** The
buggy anchor (task #579, `src/cyclerfinder/search/literature_check.py`
~line 1522-1546) is titled "...in the **Restricted** Three-Body Problem" —
i.e. it claims to cover the restricted problem (one massless body, the
model class cyclerfinder actually searches). This 2014 paper, like its 2013
sibling, is the GENERAL problem (both bodies massive) — a spacecraft/
massless-third-body candidate in the codebase's actual discovery pipeline
is NOT the same dynamical system this paper's 1:1 families were computed
in. **This paper's 1:1-resonance content cannot be used to authorize a
"published, not novel" verdict for a restricted-problem (CR3BP) candidate**
— using it that way would repeat exactly the kind of model-mismatch error
this whole #579 investigation exists to catch. It is topically adjacent
(same authors, same "1/1 resonance," same "three dimensional planetary
systems" title fragment that made it a plausible citation candidate) but
NOT a valid literature anchor for the restricted-problem gate.

## 3. Content summary (both pages, in full)

- Studies bifurcation of periodic co-orbital (1:1) orbits from planar to
  spatial as a function of planetary mass ratio ρ = m2/m1. Bifurcations to
  spatial families exist only for **ρ < 0.0205** — above that, no spatial
  1:1 families bifurcate from the planar ones at all.
- **Table 1** (p. 82, fully sourced): eccentricity pairs (e1, e2) of the
  vertical-critical orbits (v.c.o.) bounding the bifurcating segments, at
  ρ = 0.01, 0.018, 0.02, 0.0205 — two v.c.o. per ρ value (4 rows × 2
  columns of (e1,e2) pairs).
- All computed spatial periodic orbits in this family are **linearly
  stable**. The families form "bridges" connecting pairs of v.c.o., with
  mutual inclination Δi reaching a ρ-dependent maximum (Fig. 1c, up to
  ~20° at very small ρ, per the log-log plot — figure-derived, not a table
  value).
- Dynamical stability maps (FLI-based, Fig. 2) show a large-eccentricity
  regular-motion region at ρ=0.01, an isolated stable "island" near the
  periodic orbit at Δi≈10°, and stability restricted to Ω2≈270°
  (i.e. Ω2-Ω1≈180°).

## 4. Relevance / disposition

Not a codebase capability fit for the same reason the 2013 general-TBP
paper isn't (see that digest's §7) — general-TBP, not spacecraft-relevant.
Its only real value to this project is **negative**: it closes off this
specific citation candidate as a valid fix for the #579 anchor bug (it's
topically named right but model-class wrong), so #579's actual fix should
proceed as a relabel to Antoniadou-Libert 2019 (the correct paper behind
the anchor's arXiv ID) WITHOUT adding a "1:1 now covered" claim from this
paper. If a genuinely restricted-problem 1:1-resonance reference is ever
needed, this is not it — that search remains open.

## 5. References cited (p. 83)

Antoniadou, K. I. & Voyatzis, G. (2014), *Ap&SS* 349, 657 — this is the
companion paper "Resonant periodic orbits in the exoplanetary systems"
(confirmed via the earlier CrossRef sweep for task #579, DOI
10.1007/s10509-013-1679-8) that this proceedings note's "spatial general
TBP in a rotating frame" and vertical-stability method both cite; NOT yet
acquired. Hadjidemetriou, Psychoyos & Voyatzis (2009), Cel. Mech. Dyn.
Astron. 104, 23. Hadjidemetriou & Voyatzis (2011), Cel. Mech. Dyn. Astron.
111, 179. Robutel & Pousse (2013), Cel. Mech. Dyn. Astron. 117, 17.
Voyatzis (2008), ApJ 675, 802 — the DFLI diagnostic source, already noted
in the 2013 paper's digest.

End of digest.
