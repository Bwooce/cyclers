# Digest: Kumar, Anderson, de la Llave 2025

**Paper:** "Rapid GPU-Assisted Search and Parameterization-Based Refinement and Continuation of
Connections between Tori in Periodically Perturbed Planar Circular Restricted 3-Body Problems"
**Journal:** SIAM J. Appl. Dyn. Syst. **24**(1):219-258 (2025)
**Preprint:** arXiv:2109.14814v2 (2023-10-18 revision; title matches the published version exactly)
**Authors:** Bhanu Kumar (Georgia Tech), Rodney L. Anderson (JPL/Caltech), Rafael de la Llave (Georgia Tech)
**Filed:** `kumar-anderson-delallave-2025-gpu-connections-tori-perturbed-crtbp-siam-ads-24-219-arxiv-2109.14814.pdf`

Grounds task #522 (Coherent-Model Torus-Connection Search) — cross-checked directly against the
arXiv PDF (not inherited secondhand) after a prior audit flagged the citation as "not found in
corpus," which meant unacquired, not wrong: confirmed real, on-topic, and precisely matching the
title/venue/pages #522 cites.

## Core Premise

In an autonomous PCRTBP, unstable periodic orbits have 2D stable/unstable manifolds; a 2D energy
level + Poincaré section reduces connection-finding to intersecting two 1D curves in a 2D plane.
Under a TIME-PERIODIC perturbation (a third body, e.g. the Sun for an Earth-Moon system — the
BCR4BP/QBCP case), the system loses its energy integral and gains a dimension; **most unstable
periodic orbits persist not as periodic orbits but as 2D WHISKERED TORI** (Arnol'd), and it is
these tori's manifolds — not periodic-orbit manifolds — that are the generic zero-fuel connection
structures in the perturbed 5D extended phase space (x,y,px,py,θ_p).

## Method (stroboscopic-map reduction)

1. **Stroboscopic map** `F`: the time-`2π/Ω_p` map of the periodically-forced flow (`Ω_p` = forcing
   frequency). This drops the phase space from 5D to 4D and the invariant torus from 2D to a 1D
   invariant CIRCLE `K(θ): T -> R^4` satisfying `F(K(θ)) = K(θ+ω)`, `ω = 2πΩ_1/Ω_p` (§4.1, Eq 4.1-4.3).
2. **Simultaneous bundle solve**: rather than solving for `K` alone, solve jointly for `K(θ)` and
   periodic matrix-valued `P(θ), Λ(θ)` satisfying `DF(K(θ)) P(θ) = P(θ+ω) Λ(θ)` (Eq 4.4), where
   `P`'s columns are the tangent / symplectic-conjugate / stable / unstable directions at each `θ`,
   and `Λ` is block-structured with a periodic tangent-coupling term `T(θ)` plus constant Floquet
   multipliers `λ_s < 1`, `λ_u > 1` (Eq 4.5). A quasi-Newton scheme (their prior work, ref 25) solves
   this in `O(N)` storage / `O(N log N)` ops via FFT (vs `O(N^3)` for `K`-only), backed by a-posteriori
   theorems. Continuation from the unperturbed PCRTBP (`ε=0`, ordinary periodic orbit) into the
   perturbed model (`ε>0`) roughly squares the residual each Newton step.
3. **Manifold parameterization**: the torus's 2D stable/unstable manifolds are represented as a
   Fourier-Taylor series `W(θ,s) = K(θ) + Σ_{k≥1} W_k(θ) s^k` satisfying `F(W(θ,s)) = W(θ+ω, λs)`
   (Eq 4.6-4.7); `W_1` is the known stable/unstable bundle column from `P`, higher `W_k` solved
   recursively.
4. **Search**: manifolds are discretized as 2D meshes in the 4D reduced phase space; a GPU-parallel
   (Julia + OpenCL) broad-phase/narrow-phase collision-detection algorithm (adapted from computer
   graphics) finds candidate mesh intersections, giving a 5-7x CPU-vs-GPU speedup.
5. **Refinement + continuation**: candidate mesh intersections are refined to high accuracy using
   the Fourier-Taylor manifold parameterizations, then continued through 1-parameter families of
   tori (exploiting the families' Whitney differentiability to interpolate needed parameterizations).

## Validation / Test System

Demonstrated on the **Jupiter-Europa planar elliptic RTBP (PERTBP)**, NOT Sun-Earth-Moon: `mu ≈
2.527e-5`, eccentricity (perturbation strength) `ε = 0.0094` — the PERTBP models the Jupiter-Europa
pair's own orbital eccentricity as the periodic perturbation (not a third body), with perturbation
frequency `Ω_p = 1` by construction (time normalized to the primaries' own period). No Sun-Earth or
Earth-Moon-Sun numbers are given in the pages read (pp. 1-8; validation tables are later in the
38-page paper, not yet extracted — read on demand if the GPU/mesh machinery specifically is built).

## Relevance to #522

This paper's structural setting — a periodically-perturbed autonomous CR3BP — is EXACTLY the
`core/bcr4bp.py` BCR4BP already in this codebase (Earth-Moon primaries + Sun as periodic perturber,
reduces to CR3BP at `mu_S=0`). The bundle-equation solve + Fourier-Taylor manifold machinery
described here is genuinely new/unbuilt in this codebase (confirmed absent by two independent
audits, 2026-07-02/03). It is the HEAVIER of #522's two candidate first-build paths — see
[[2026-07-03-digest-owen-baresi-2024-knot-theory-heteroclinic]] for a lighter-weight alternative
(topological linking-number screening) that reuses existing GMOS torus-generation code already in
this repo (`genome/qp_tori.py`) and is the recommended first cut.
