# Ross, Koon, Lo, Marsden 2003 — Design of a Multi-Moon Orbiter (AAS 03-143)

Digested 2026-06-30. Ross-group Jovian moon-tour paper, filed in relation to open threads
#318 (joint moon-tour search), #465 (multi-rev leveraging), #500 (Keplerian-map genome).

**Source (cite exactly, no file path):**
Ross, S.D., Koon, W.S., Lo, M.W., Marsden, J.E., "Design of a Multi-Moon Orbiter,"
AAS Paper 03-143, 13th AAS/AIAA Space Flight Mechanics Meeting, Ponce, Puerto Rico,
9–13 February 2003.

> Text layer present and machine-readable. All figures are ASCII-degraded but no key
> tabular data is in figures only. Full paper (15 pp. + references) read in its entirety.

---

## 1. Contribution in three lines

Introduces the **Multi-Moon Orbiter (MMO)** concept: a single spacecraft orbits several
Jovian moons in sequence (Callisto → Ganymede → Europa) using **resonant gravity assists**
outside the moon's sphere of influence combined with **ballistic capture/escape** via
invariant manifold tubes. Demonstrates a tour requiring only 22 m/s deterministic ΔV from
an initial jovicentric orbit; Europa orbit insertion costs an additional ~450 m/s.

---

## 2. Model and method

**Dynamical model:** Patched Planar Circular Restricted 3-Body Problem (PCR3BP) as building
blocks for a restricted 5-body problem (R5BP: Jupiter + Europa + Ganymede + Callisto +
spacecraft). Moons assumed coplanar, circular, no mutual gravitational interaction.

**Three building blocks:**
1. **Resonant gravity assists** — multiple near-resonance periapsis passages outside the
   Hill sphere reduce (or increase) the jovicentric semimajor axis; Tisserand parameter /
   Jacobi constant approximately conserved. Illustrated by e vs. a trajectory plot
   (Fig. 7a) and Period vs. periapse plot (Fig. 7b) following Heaton–Strange–Longuski–
   Bonfiglio convention.
2. **Ballistic capture/escape** — transit orbits inside the stable/unstable manifold tubes
   of Lyapunov orbits around L1/L2 carry the spacecraft from the exterior realm into
   temporary capture around a moon and back. No fuel used during capture/escape.
3. **Small impulsive maneuvers** — ≤2 m/s each, performed at opposition (spacecraft and
   moon in opposition as seen from Jupiter) to steer the spacecraft between the natural
   resonance-lock phases that determine the sign of the gravity-assist kick.

---

## 3. Sourced goldens (exact numbers)

All numbers from the preliminary numerical example trajectory in the paper. No tables with
ICs appear; values are scattered in text and figures.

| Quantity | Value | Location |
|---|---|---|
| Total deterministic ΔV (entire tour) | **22 m/s** | Abstract, p.1, p.9 |
| Time of flight (entire tour) | **~4 years** | Abstract, p.9 |
| Individual maneuver cap (arbitrary, engine-dependent) | **2 m/s** | p.6, footnote e |
| Callisto closest approach altitude | **~1400 km** | p.9 |
| Ganymede closest approach altitude | **~2100 km** | p.9 |
| Europa orbit insertion ΔV (100 km circular) | **~450 m/s** | pp.9, 11 |
| Time inside 12 R_J before Europa insertion | **~260 days** | p.11, footnote f |
| Jupiter-Callisto mass parameter µ (implied) | **5.667 × 10⁻⁵** | (Ross-Scheeres 2007 companion paper; consistent) |

**Earth-Moon comparison trajectories (sourced to Bollt&Meiss 1995 and Schroer&Ott 1997,
reproduced in Fig.9):**
- Bollt & Meiss [1995]: ΔV = 749.6 m/s, TOF = 748 days; same initial (59669 km circular)
  and final (precessing ellipse, perilune 13970 km) orbits
- Schroer & Ott [1997]: ΔV = 748.9 m/s, TOF = 377.5 days (same orbit pair)
- Ross [2003] (present work / Fig.8b): ΔV = **860.1 m/s**, TOF = **65 days** (same pair)
- Hohmann baseline (same pair): ΔV = **1220 m/s**, TOF = **6.6 days** (Fig.9 legend)

**Orbital inclinations of Galilean moons w.r.t. local Laplace plane** (footnote c, p.5):
- Europa: 0.467°, Ganymede: 0.172°, Callisto: 0.306°
- Orbital eccentricities: Europa 0.0002, Ganymede 0.0011, Callisto 0.0074

**No tabulated ICs, Jacobi constants, or resonance ratios for the MMO tour trajectory
itself appear in this paper.** The specific (a, e, ω) state at any epoch is
figures-only in this conference paper. The companion Chaos 2000 paper (Koon et al.)
carries the numerical detail.

---

## 4. Key physical insight for our threads

**Resonant gravity assist geometry:** Pass apojove *slightly ahead of the outer moon* (ω
slightly positive in rotating frame) → perijove decreases, a decreases. Pass perijove
*slightly ahead of the inner moon* → apojove decreases, orbit circularizes toward that
moon. Maneuvers at opposition ensure near-resonance is maintained across multiple passes.
This is the same mechanism as the Keplerian map in Ross-Scheeres 2007.

**Tisserand plot (e vs. a, Fig.7a):** the spacecraft jumps between branches of constant
Jacobi constant (= Tisserand parameter in the approximation) as control switches from one
moon to another. Orbital energy decreases via resonant GA from the outer moon; the
switch point is when perijove of the jovicentric orbit grazes the inner moon's orbit (point
E in Fig.3).

**Capture-then-escape:** ballistic capture puts spacecraft into an unstable orbit around a
moon; escape via very small ΔV is the time-reversed capture. Controllable via manifold
tube entry angle.

---

## 5. Reuse verdict vs. open threads

### #318 (joint moon-tour search — joint_cell.py, joint_sobol.py)

**MODERATE — conceptual building blocks, limited algorithmic detail here.**
The MMO paper establishes the *components* of a multi-moon tour (resonant GA → capture →
escape) but the detailed manifold-tube-intersection algorithm is in the companion papers
(Koon et al. 2000 Chaos, Koon et al. 2002 Contemp. Math — see separate digests). The
Tisserand / Jacobi constant diagram (Fig.7) is directly analogous to the Tisserand plane
approach already in our codebase. The key additive point is: **the "switching control"
criterion** (outer moon's orbit grazes perijove of inner-moon's resonant orbit) is a
natural gate to embed in joint_cell's cell-connectivity logic. Currently this paper
provides conceptual validation rather than a ready algorithm to import.

### #465 (multi-rev leveraging)

**MODERATE — the resonant GA description is the physical mechanism underlying multi-rev.**
Multiple near-resonance periapsis passages with a moon is exactly the multi-rev
leveraging scenario. The paper shows this can reduce a to any desired value with only
~2 m/s steering per pass. For #465: the maneuver-at-opposition timing rule gives a
concrete criterion for when to apply steering impulses in a multi-rev arc. The companion
Ross-Scheeres 2007 (Keplerian map) provides the quantitative tool; this paper provides
the engineering motivation and example.

### #500 (Keplerian-map genome — SPECULATIVE)

**LOW — this paper motivates the genome but the Keplerian map tool lives in Ross-Scheeres
2007.** The resonant GA building block described here is the physical phenomenon the
Keplerian map models. This paper's Fig.3 (geometry of resonant GA) and Fig.7 (e vs. a
Tisserand jumps) give a clean mission-design picture of what the map genome encodes.
Worth citing as the mission application context.

### #494 (binary/circumbinary)

**NOT RELEVANT** — purely heliocentric/jovicentric single central body. No binary/
circumbinary content.

---

## 6. Corpus-index status

**Verdict: `mined`** — full paper read, all sourced numbers extracted; no IC tables to
transcribe to catalogue. No catalogue rows (jovian moon tour, not an Earth-Mars cycler
or class-in-scope object).

---

## 7. Not extractable

- Figure 1 (full tour in jovicentric inertial frame): no numeric ICs
- Figure 6 (semimajor axis vs. time): qualitative only; no epoch-tagged states
- Figure 7 (e vs. a Tisserand diagram): resonance labels visible (E=Europa, G=Ganymede,
  C=Callisto at a=1, 1.5, 3.5 in Europa units) but no tabular data
- Figure 7b (period vs. periapse): referenced to Heaton-Strange-Longuski-Bonfiglio
  convention; axis values not recoverable from text-layer extraction
- The specific ΔV components (how the 22 m/s is distributed across passes) are not
  individually listed in the paper text
