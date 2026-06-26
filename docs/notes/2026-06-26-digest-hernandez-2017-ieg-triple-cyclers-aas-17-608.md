# Digest — Hernandez, Jones & Jesick (2017), "One Class of Io-Europa-Ganymede Triple Cyclers" (AAS 17-608)

**Digested:** 2026-06-26 (read page-by-page from the acquired PDF; every value below is
sourced to a specific page/table — no figure-read estimates unless flagged).
**Task:** #482 (acquire + digest the genuine Galilean reference for #480's reproduction golden).
**System:** Jovian. **Bodies:** {Io, Europa, Ganymede}. **DOI:** none (AAS/AIAA conference paper).

## Citation (grounded)
Sonia Hernandez, Drew R. Jones, Mark Jesick, "One Class of Io-Europa-Ganymede Triple
Cyclers," **AAS 17-608**, AAS/AIAA Astrodynamics Specialist Conference, Stevenson WA
(Columbia River Gorge), Aug 20-24 2017; Advances in the Astronautical Sciences Vol. 162
(Univelt), pp. 973-984. JPL/Caltech. (Title page + © page confirm AAS 17-608 and JPL.)

**Two-paper disambiguation (the concept-collision that earlier tripped the main session):**
The SAME authors published a *sibling* paper the same year on a DIFFERENT system —
"Low Excess Speed Triple Cyclers of Venus, Earth, and Mars," **AAS 17-577** (author order
Jones/Hernandez/Jesick), which is this paper's own reference [6]. This Galilean paper is
17-608 (Hernandez-first); the VEM paper is 17-577 (Jones-first). Same year, same "triple
cycler" term, different systems. Catalogue row `hernandez-2017-jovian-ieg-triple-family`
had the AAS number unconfirmed (candidate 17-462/17-608) — **it is 17-608**; correct it.

## Method (pp. 1-7)
- **Tisserand-graph initial guess** (Fig 1): v∞ contours 1-15 km/s for Io/Europa/Ganymede;
  Ganymede reachable via Io flyby at v∞>4 km/s, via Europa flyby at v∞>2 km/s.
- **Ideal ephemeris model** (p.3): moon orbits circular + planar, periods set so a 5.2°
  inertial displacement occurs each synodic period. Given a_Io: a_Eur=((8π+Δ)/(4π+Δ))^(2/3)·a_Io,
  a_Gan=((8π+Δ)/(2π+Δ))^(2/3)·a_Io, Δ=5.2°.
- **Conic initial-guess tool** (pp.4-5): trajectory = 2-body Jovian conic classified by
  (period = integer × T_syn, eccentricity, ω). Search box R_Jup ≤ r_p ≤ a_Io and r_a ≥ a_Gan.
  Table 1 (allowed s/c revs : synodic period): syn 1→{1:1,1:2}; syn 2→{2:1,2:3,2:5};
  syn 3→{3:1,3:2,3:4,3:5,3:7}; syn 4→{4:1,4:3,4:5,4:7,4:9}.
- **Lambert legs + powered periapsis maneuver** (pp.6-7): zero-radius SOI patched conic;
  ΔV impulse at periapsis if needed for v∞ mismatch (Eqs 3-7, B-plane formulation).
- **Monte Carlo Lambert optimization** (p.7): flyby altitudes constrained **25 km–70,000 km**;
  search params = departure phase + N-1 times-of-flight; ToF variation ~10% of body period
  (Io ≈ 4 h, Ganymede ≈ 17 h).
- **Two-level real-eph conversion** (p.10): ideal solution → Lambert in real ephemeris →
  full high-fidelity optimization. Ballistic repeatability of the high-fidelity solution
  "in general will only last for a few cycles," then maintenance ΔV required.

## Grounded constants
- **1:2:4 Laplace resonance** (Ganymede:Europa:Io periods); ref [8] Sinclair 1975.
- **T_syn = 7.05 days** = 1 Ganymede orbit = 2 Europa = 4 Io orbits (p.2).
- **5.2° inertial shift per synodic period** (p.2).
- Io orbital period 1.75 d (p.4). Flyby-altitude window 25 km–70,000 km (p.7).
- Sequence naming: first letter of each flyby body, Europa chosen as first by convention
  (EGIE ≡ GIEG). #sequences for n encounters: s(n)=3^(n-1)-2^n+1 (Eq 2): 3→2, 4→12, 5→50, 6→180.

## The printed solutions = the invariants available for #480's golden
Only THREE worked examples are tabulated (the paper states "hundreds of solutions found"
but prints no family census — a publication-gap caveat: the golden validates against THESE
examples, not a census).

### Table 3 — EGIEIE, 1 synodic period, depart 03-Oct-2020 ET (p.8)
| Flyby | ToF (d) | V∞ (km/s) | Altitude (km) | ΔV (m/s) |
|---|---|---|---|---|
| Europa | – | 12.27 | 9,260 | 1.8 |
| Ganymede | 2.03 | 7.29 | 69,990 | 0.1 |
| Io | 0.85 | 15.77 | 3,176 | 4.4 |
| Europa | 0.66 | 12.15 | 2,259 | 0.6 |
| Io | 2.87 | 15.86 | 494 | 0.0 |
| Europa | 0.65 | – | – | – |
| **Total** | **7.06** | – | – | **6.90** |
High excess speeds (12-16 km/s) — "hard to fly from a navigation point of view" (p.8).

### Table 4 — EGGIE, 4 synodic periods, depart 29-Sep-2020 (p.10) — THE LOW-ΔV PRIZE
| Flyby | ToF (d) | V∞ (km/s) | Altitude (km) | ΔV (m/s) |
|---|---|---|---|---|
| Europa | – | 9.12 | 1,444 | 0.00 |
| Ganymede | 1.59 | 7.07 | 2,155 | 0.60 |
| Ganymede | 8.60 | 7.07 | 6,263 | 0.00 |
| Io | 7.34 | 8.38 | 653 | 0.10 |
| Europa | 10.69 | – | – | – |
| **Total** | **28.22** | – | – | **0.70** |
Low excess speeds (7-9 km/s) → navigation-viable; total ΔV 0.70 m/s ("can probably be
optimized to zero"). Both Ganymede flybys at equal v∞ = 7.07 km/s (general property:
adjacent same-body encounters in a ballistic cycler occur at equal v∞).
**DATE INCONSISTENCY (the PAPER'S, not ours):** Table 4 caption says **29-Sep-2020**;
the page-9 text and Figure 4 caption say **02-Oct-2020**. Flag both; resolve at reproduction.

### Real-ephemeris EIGE, 1 synodic period / 1 rev (Fig 5, pp.10-11) — the direct #480 V3/V4 analog
Computed in the REAL ephemeris via Lambert arcs. Flyby altitudes: **Io 2,817 km,
Ganymede 13,180 km, Europa 470 km.** Ideal-model equivalent is ballistic; first cycle
optimized to ballistic, but after **10 repeat cycles ΔV grows to almost 30 m/s**.
High-fidelity gravity conversion: ballistic for **2 cycles**, then large maintenance impulses.

## Relevance to #480 (the reproduction target)
- The **EIGE real-ephemeris example** (flyby alts 2817/13180/470 km, ballistic→~30 m/s over
  10 cycles) is the closest analog to #480 M1's V3/V4 real-eph reproduction — the authors'
  own two-level conic→real-eph path is exactly #480 Approach A.
- The **EGGIE 4-synodic ballistic** cycler (0.70 m/s, equal-7.07-km/s Ganymede flybys) is the
  low-ΔV ideal-model golden.
- **#480 invariant golden** (sourced) = the resonance ratios (1:2:4), the per-flyby V∞ levels
  + ToFs + altitudes + ΔV of these tables, and T_syn=7.05 d. NOT a family census (not printed).

## Corpus relationships / further references (from the reference list, p.12)
- **Sibling:** AAS 17-577 VEM triple cyclers (ref [6]) — catalogued `jones-2017-vem-triple-family`.
- **Prior Jovian triple-cycler work:** A. E. Lynam & J. M. Longuski, "Laplace-resonant
  triple-cyclers for missions to Jupiter," Acta Astronautica 69(3-4) 2011, pp.158-167,
  **doi:10.1016/j.actaastro.2011.03.011** (ref [7] — this paper's "investigated for the first
  time in Reference 7"). HAS a DOI; worth acquiring for a second IEG invariant source.
- **Foundational moon cyclers:** R. P. Russell & N. J. Strange, "Planetary Moon Cycler
  Trajectories," JGCD 32(1) 2009, pp.143-157 (ref [5]).
- **Laplace resonance source:** A. T. Sinclair, "The orbital resonance amongst the Galilean
  satellites of Jupiter," Celestial Mechanics 12 1975, pp.88-96 (ref [8]).
