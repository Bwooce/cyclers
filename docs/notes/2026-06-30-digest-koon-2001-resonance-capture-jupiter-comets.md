# Koon, Lo, Marsden & Ross 2001 — Resonance and Capture of Jupiter Comets (CMDA 81) — digest

Digested 2026-06-30. Focus: L1/L2 invariant manifold tubes as conduits for resonance
transition and temporary Jupiter capture; primary input for **#267 (resonance_network /
resonant transport)**.

**Source (cite exactly):**
Koon, W. S., Lo, M. W., Marsden, J. E., & Ross, S. D., "Resonance and Capture of
Jupiter Comets," *Celestial Mechanics and Dynamical Astronomy*, Vol. 81, 2001,
pp. 27–38. DOI: 10.1023/A:1013359120468. Caltech & JPL.
© 2001 Kluwer Academic Publishers.

> 12 pp., clean digital typeset. A concise summary/companion to the full treatment
> in Koon et al. 2000 (Chaos 10:427–469), which proved existence of heteroclinic
> connections and the itinerary framework. This paper focuses specifically on the
> resonance transition of Jupiter family comets (Oterma, Gehrels 3) and makes the
> connection to mean motion resonances explicit with Poincaré sections.

---

## 0. HEADLINE VERDICT

**The paper identifies the tube-intersection mechanism linking the 3:2 (interior) and
2:3 (exterior) Jupiter resonances via the L1/L2 heteroclinic chain.** For #267
(resonance_network): this is the canonical source establishing that the first Poincaré
cut of the L1 stable/unstable manifold tubes intersects at the 3:2 resonance, and
the L2 tubes intersect at the 2:3 resonance — and these are dynamically connected
through the Jupiter-capture region. This explains why spacecraft / comets with this
energy hop between first-order (|p-q|=1) resonances rapidly. Higher-order resonances
require further Poincaré cuts.

---

## 1. METHOD / MODEL

### 1.1 Dynamical Model

**Planar Circular Restricted Three-Body Problem (PCR3BP)** for Sun-Jupiter system:

    ẍ - 2ẏ = Ω_x,    ÿ + 2ẋ = Ω_y

where Ω = (x² + y²)/2 + (1-μ)/r_S + μ/r_J.

- Total mass normalized to 1: m_S = 1-μ, m_J = μ
- Sun-Jupiter distance normalized to 1
- Angular velocity normalized to 1
- Jacobi constant: C = -2E, where E = ½(ẋ² + ẏ²) - Ω(x,y)
- Energy manifolds: 3D surfaces foliated by the 4D phase space; Poincaré sections 2D

### 1.2 Phase Space Structure

**Equilibrium points L1, L2** (collinear, unstable): saddle × center linearization.
For energies just above L2's energy, Hill's region connects three regions:
- **Interior region (S)**: inside Jupiter's orbit
- **Jupiter-capture region (J)**: bubble around Jupiter
- **Exterior region (X)**: outside Jupiter's orbit

Connected via bottleneck "necks" around L1 (S↔J) and L2 (J↔X).

**Four orbit types** in each equilibrium region (Conley 1968): (1) Lyapunov periodic
orbit, (2) asymptotic orbits on stable/unstable manifolds, (3) transit orbits (cross
the neck), (4) non-transit orbits (reflect back).

**Invariant manifold tubes**: 2D stable (W^s) and unstable (W^u) manifolds of the
Lyapunov orbits form cylindrical tubes. Transit orbits pass *inside* the tube;
non-transit orbits are outside. These tubes extend globally as phase-space conduits.

### 1.3 Resonance Transition Mechanism

**Rapid transition** (< 1 Jupiter period) from exterior to interior region:
- Comet enters the L2 stable manifold tube (exterior)
- Passes through L2 neck into Jupiter capture region
- Intersects L1 stable manifold tube
- Passes through L1 neck into interior region

The key step is the intersection of L2 unstable tube (in Jupiter region, from exterior)
with L1 stable tube (in Jupiter region, toward interior). Visualized on a Poincaré
section (vertical line through Jupiter), these appear as distorted circles that overlap.

**Connection between resonances** is computed via Poincaré sections in Delaunay
variables. Tube cross-sections appear as closed curves; their intersections mark the
resonant families.

---

## 2. SOURCED GOLDENS (verbatim with page references)

| Quantity | Value | Source |
|---|---|---|
| Sun-Jupiter mass ratio | μ = 9.537×10⁻⁴ | p. 29 |
| Jupiter orbital eccentricity | 0.0483 (stated as playing "little role during fast resonance transition") | p. 29 |
| Fast transition timescale | "less than one Jupiter period" | p. 32 |
| Jupiter orbital period (implied) | ≈ 12 years (stated as reference: "Jupiter period ≈ 12 years") | p. 32 |
| Interior resonance at first L1 tube cut | **3:2 mean motion resonance** (interior region, U1 section) | p. 33–34 |
| Exterior resonance at first L2 tube cut | **2:3 mean motion resonance** (exterior region, U4 section); secondary: 1:2 | p. 34 |
| First-order resonance condition | \|p-q\| = 1 (first Poincaré cut gives only first-order resonances) | p. 35 |
| Oterma time interval covered | AD 1910–1980 (Figure 1a) | p. 28 |
| Heteroclinic connection paper | Koon et al. 2000, Chaos 10(2):427–469 | p. 28, 31 |

**Key qualitative result with quantitative implication**: "the cross-section of the
tube is widest near the 2:3 resonance" (exterior, p. 34) — the 1:2 exterior
intersection exists but is narrower, making the 2:3 the dominant exterior resonance
for this energy range. This is a selectivity result about the relative "gate width"
between resonances.

**Higher-order resonances**: "Looking at cuts beyond the first reveals transitions
between higher order resonances. In addition, higher energies have 'larger', more
dispersive tubes, which have more intersections for a given cut number" (p. 35).

---

## 3. REUSE VERDICTS

### #267 (resonance_network / resonant transport)

**STATUS: FOUNDATIONAL — MINED. This is the canonical tube-resonance-connection
reference for the #267 resonance network.**

Key results directly reusable:

1. **Tube-resonance correspondence**: L1 stable/unstable manifold tubes intersect the
   interior phase space at the 3:2 resonance (first Poincaré cut at U1). L2 tubes
   intersect the exterior phase space at the 2:3 resonance (first cut at U4). This
   is the "address" of each resonance in phase space.

2. **Resonance network topology**: the two resonances (3:2 interior, 2:3 exterior)
   are dynamically linked via tube intersections in the Jupiter-capture region. The
   intersection set Λ (in the Jupiter region) is the "gateway" connecting them. Any
   trajectory through Λ makes the Oterma-like 2:3 → 3:2 transition.

3. **Resonance selectivity by tube width**: widest tube cross-section at 2:3 in
   exterior (dominant gateway); 1:2 intersection also present but narrower (secondary
   gateway at the same energy). This gives a quantitative measure of which resonances
   are "easy to reach."

4. **Energy dependence**: at higher energies, tubes are larger → more intersections →
   more resonances reachable per Poincaré cut number. The resonance network is
   energy-stratified.

5. **Higher-order resonances via higher Poincaré cuts**: first cut → first-order
   (|p-q|=1); second cut → second-order, etc. This provides the framework for
   constructing a full resonance network by iterating the Poincaré map.

For the resonance_network implementation: the lobe dynamics (strips P(Δ_X) and
P⁻¹(Δ_S) in the Jupiter region, their intersection Λ) provide the quantitative
transport measure between resonances. The paper stops at the first cut / qualitative
description; the 2000 Chaos paper (same authors) has the full theoretical machinery.

### #500 (Controlled Keplerian map as moon-tour cycler genome)

**STATUS: BACKGROUND THEORY — INDIRECTLY RELEVANT.** The L1/L2 tube structure
described here is the same structure that the Keplerian map (Grover-Ross 2009)
approximates. The map's stable resonant islands in the (ω, a) phase space correspond
exactly to the resonant families seen in the U1/U4 Poincaré sections here. The "lanes
of fast migration" between islands in the Keplerian map are the chaotic zones between
the tube cross-sections here.

For #500: this paper explains WHY the Keplerian map's resonance islands are at the
correct resonance semimajor axes and WHY the chaotic sea mediates inter-resonance
transport. It is not an implementation source but is the theoretical underpinning.

### #494 (k1,k2 binary cyclers / circumbinary)

**STATUS: NOT RELEVANT.** This paper treats the Sun-Jupiter PCR3BP; there is no binary
primary system and no circumbinary regime. The L1/L2 structure is for a single large
primary (Sun) + one smaller body (Jupiter). Not applicable to binary-primary cyclers.

---

## 4. NOTES ON RELATIONSHIP TO OTHER PAPERS

- **Koon et al. 2000 (Chaos 10:427–469)**: the parent paper — full theoretical proof
  of heteroclinic connections, itinerary framework, numerical construction of
  transition orbits. This 2001 paper is a focused summary applying those results to
  the comet resonance problem.
- **Lo & Ross 1997 (JPL IOM 312/97)**: first paper connecting invariant manifolds to
  Oterma's orbit — predecessor.
- **Koon et al. 2004 NYAS 1017**: uses same tube/lobe dynamics framework applied to
  asteroid binaries (digested separately).
- **Grover-Ross 2009**: the Keplerian map as an analytical Poincaré map of the same
  PCR3BP, providing a fast propagator for the dynamics described here.

---

## 5. STATUS

**DIGESTED.** Not catalogue-eligible (comet dynamics paper; no repeating spacecraft
cycler). Maps to #267 (resonance_network) as foundational theory.
