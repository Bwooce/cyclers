# #434 3D / broken-plane cycler discovery — verdict

**Date:** 2026-06-24
**Status:** COMPLETE — capability delivered + spike-validated; structural NEGATIVE
for a *novel* out-of-plane family, CONDITIONAL on a narrow machine-loadable seed
set. No catalogue writeback.

## What was run

`scripts/scan_434_3d_broken_plane_em.py` lifts planar Earth-Moon cyclers into
z≠0 and continues them with the delivered #291 3D genome (Phase-1 corrector +
Phase-2 pseudo-arclength tracer), classifying each member's out-of-plane
topology (k1,k2,k_z) and screening against a new spatial-CR3BP literature corpus.

**Task-0 GO/NO-GO spike gate passed first:** the tracer reproduced the #287
Braik-Ross (1,1) 3D family (201 members, closure < 1e-8) before any build.

Two seed routes (the planar-collapse finding from the gate — small z0 lifts
collapse to the planar manifold — drove the design):
- **route (i)** L1 vertical-Lyapunov (`lyapunov_seed_3d`): locked, 201 members.
  L2 did not converge to a non-planar orbit (skipped, honestly logged).
- **route (ii)** Braik-Ross C11a (1,1) z0-amplitude lock from |z0|∈{0.05..0.24}:
  **locks=5, collapses=0** — every amplitude locked onto a 3D branch.

Output: `data/scan_434_3d_broken_plane_em.jsonl` (1113 members), wall 2903 s.

## Result

| metric | value |
|---|---|
| converged 3D members | 1113 |
| closure residual | median 2.5e-11, max 6.8e-9 (all < 1e-8) |
| underlying families | 2: C11a (1,1) 3D extension + L1 vertical-Lyapunov |
| (k1,k2,k_z) clusters | (1,0,6)=603, (0,0,2)=148, (1,1,8)=139, (1,2,10)=110, (1,1,10)=60, (0,1,2)=53 |
| stability | overwhelmingly hyperbolic_pair; 14 stable in L1, 23 unstable in one C11a branch |
| novel family | **none** |

The C11a z0-lifts at 0.05 / 0.15 / 0.24 all lock to z0≈−0.241, C∈[2.99,3.12] —
this **is** the #287 spike's (1,1) 3D extension, reproduced and mapped across
z-amplitude. The 0.10 / 0.20 lifts reach a smaller-amplitude branch (z0≈−0.18).
Route (i) maps the L1 vertical-Lyapunov family (C∈[3.00,3.17]). The six topology
clusters are (k1,k2,k_z) combinations of these **two known family types**.

## Reading (honest scope of the negative)

The literature gate returned `not-found` for every cluster — but this is a
**sparse-corpus false-negative, NOT a novelty signal.** The new `known_corpus_3d`
holds only 3 anchors (Howell 1984 halos, Folta-Bosanac-Guzzetti-Howell 2015
NRHO catalog, Antoniadou-Voyatzis 2018 spatial resonant families), too thin to
tuple-match these (k1,k2,k_z) signatures. By construction the families ARE known:
route (ii) seeds *from* the Braik-Ross (1,1) cycler (its 3D extension = the #287
spike, already characterized) and route (i) seeds the textbook L1 vertical mode.
Per the literature-novelty discipline ("not-found" is necessary-not-sufficient;
human judgment governs), **no novel out-of-plane family was discovered.**

The negative is **conditional on a narrow seed set**, and the condition matters:
**only the Braik-Ross C11a (1,1) planar root was machine-loadable** as a
rotating-frame IC (the catalogue rows carry heliocentric mission elements, not
CR3BP ICs — the scan hardcodes C11a from the #287 on-disk planar baseline). So
the campaign mapped the 3D extensions of *one* planar cycler + the L1 vertical
family. The realistic novel frontier — *other* planar Aldrin/Braik-Ross roots
lifted to 3D, and other systems (Sun-planet, other moons) — was **not reached**
for lack of machine-extractable seed ICs, not because it is empty.

So the precise claim: *the z0-lift + pseudo-arclength method finds no novel
out-of-plane cycler beyond the (already-known) C11a (1,1) 3D extension and the
L1 vertical-Lyapunov family, at the one planar root currently seedable.*

## Disposition

- **No catalogue writeback** (the families are known; nothing novel cleared the
  gate; matches #432/#433/#435/#436 discipline).
- **Capability delivered + validated:** 3D vertical-Lyapunov seed generator,
  z-oscillation topology (k_z), z0-lift + pseudo-arclength discovery sweep
  (spike-reproduction golden), and the spatial-CR3BP literature corpus + matcher.
- Method-versioned negative registered in `data/empty_regions.jsonl`.
- **Follow-on (the real unlock):** make the catalogue's planar Aldrin/Braik-Ross
  rotating-frame ICs machine-loadable so the sweep can lift the *full* planar
  cycler set (not just C11a) into 3D — that is where a genuinely novel
  broken-plane cycler, if one exists, would be found. The infrastructure is ready
  to consume them.
