# Digest: Owen & Baresi 2024

**Paper:** "Applications of knot theory to the detection of heteroclinic connections between
quasi-periodic orbits"
**Journal:** Astrodynamics **8**, 577-595 (2024)
**DOI:** 10.1007/s42064-024-0201-0
**Authors:** Danny Owen, Nicola Baresi (Surrey Space Centre, University of Surrey)
**Filed:** `owen-baresi-2024-knot-theory-heteroclinic-connections-quasiperiodic-orbits-astrodynamics-doi-10.1007-s42064-024-0201-0.pdf`
(conference precursor also fetched: ISSFD 2024, same content, `issfd.org/ISSFD_2024/ISSFD2024_19-5.pdf`)

Grounds task #522 — see [[2026-07-03-digest-kumar-anderson-delallave-2025-whiskered-tori-connections]]
for the companion citation and the same acquisition/grounding context.

## Core Idea

Heteroclinic connections between quasi-periodic tori require intersecting the tori's stable and
unstable manifolds — both 3D objects in the CR3BP's 5D fixed-Jacobi phase space (co-dimension 2
each, so their generic intersection is 1D — Tables 1-2, §1: *higher*-dimensional intersection than
periodic-orbit manifolds, which only intersect at isolated points, hence quasi-periodic connections
are argued to be more common). Finding this intersection directly in 4D reduced phase space is
hard. This paper's contribution: a **purely topological SCREENING method** (linking number, from
knot theory) that cheaply flags where along a 1D scan a connection exists, before any expensive
numerical refinement — avoiding both blind brute-force search and human-in-the-loop plotting.

## Method

1. **Torus generation** — reuses **GMOS** (Gómez-Mondelo-Olikara-Scheeres) quasi-periodic-torus
   parameterization *unmodified* (§2.3, citing ref 26) for its computed Floquet matrix per invariant
   circle, which gives the stable/unstable eigenvectors `δU_s(θ)`, `δU_u(θ)` directly as a byproduct.
2. **Manifold generation** — perturb `N×M` points spanning the torus by `ε·δU_{s,u}(θ)`
   (`ε=1e-6`), propagate forward (unstable) / backward (stable) to a fixed surface of section
   (`x = 1-μ`, the secondary's position) — standard manifold globalization, same pattern this
   session's #530/#531 already used for periodic orbits. Near-Moon propagation is regularized via a
   Kustaanheimo-Stiefel-type transform (Eq 7-18) purely to keep torus-map interpolation well behaved
   (removes a coordinate singularity, not a physics requirement).
3. **Torus maps** — the ensemble of propagated endpoint states, keyed by their originating `(θ_0,
   θ_1)`, is interpolated into continuous scalar "torus map" functions per state component.
4. **Level curves + reduced closed curves** — pick a scanning variable `D` (one state component);
   extract the `θ` level curve where the torus map equals a chosen `D` value; use that curve to
   interpolate 3 OTHER state components `(A,B,C)` at the same `θ` values, producing a **closed curve
   in 3D `(A,B,C)` space** — one from the stable manifold's torus map, one from the unstable.
5. **Linking number** — a topological invariant of two closed curves that can ONLY change if the
   curves pass through each other. Computed via explicit computational geometry (§3.1, Eq 19-20):
   triangulate a fan surface bounded by one curve (from its centroid), then for every line segment of
   the OTHER curve, test intersection against every triangle (three dot-product sign checks per
   segment-triangle pair) and accumulate signed crossing counts.
6. **Scan + detect** — sweep the scanning variable `D` densely; **any change in linking number
   between consecutive `D` values is a heteroclinic connection** at that `D`. This converts an
   expensive N-D search into a cheap 1D sweep of a purely topological, robust invariant.
7. **Initial guess extraction** — at a detected `D`, interpolate `(A,B,C)` from both curves, average
   nearest-neighbor matches, use the Jacobi integral to recover the 6th state coordinate — giving a
   full 6D initial guess for a differential corrector (pages beyond p.10 — the actual refinement step
   and Sun-Earth/Earth-Moon/Jupiter-Ganymede numerical results were not yet extracted; read on
   demand before/while implementing).

## Test Systems (from abstract; page-10 numbers not yet extracted)

Demonstrated on **Earth-Moon, Sun-Earth, and Jupiter-Ganymede CR3BPs** — each a SEPARATE autonomous
CR3BP (this paper does NOT itself demonstrate a genuine cross-system SE<->EM connection within one
coherent 4-body model; each demonstration is single-system, e.g. EM-internal L1<->L2-family
connections). This matters for #522: Owen & Baresi's own worked example is a same-system positive
control, not evidence the method already spans SE<->EM.

## Relevance to #522 — RECOMMENDED FIRST-BUILD PATH

This method is substantially lighter-weight than
[[2026-07-03-digest-kumar-anderson-delallave-2025-whiskered-tori-connections]]'s GPU
mesh-collision-detection + Fourier-Taylor manifold machinery: no bundle-equation quasi-Newton
solve, no Fourier-Taylor recursion, no GPU/OpenCL requirement. The heaviest reused dependency — GMOS
torus generation — **already exists in this codebase** (`genome/qp_tori.py`, confirmed by two
independent audits, 2026-07-02/03). The genuinely new pieces (torus-map interpolation, level-curve
extraction, triangulated linking-number computation) are lightweight computational geometry,
buildable and testable in a single session, with a natural closed-form-adjacent positive control
(reproduce Owen & Baresi's own single-system Earth-Moon L1<->L2 linking-number result before
attempting any cross-system claim).

**Caveat found during #522 scoping (2026-07-03, same-day frame-gap audit):** applying EITHER paper's
method to a genuine SE<->EM cross-system connection (not just single-system) requires a model where
the SE-EM relative phase is a genuine dynamical state, not a searched/patched unknown — this
codebase does not yet have that (see the #522 OUTSTANDING.md entry for the full finding). This
digest's "recommended first-build path" note applies to a SINGLE-SYSTEM validation build; the
cross-system application is gated on that separate finding.
