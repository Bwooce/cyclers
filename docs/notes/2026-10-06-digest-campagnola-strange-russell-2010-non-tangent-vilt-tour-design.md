# Digest: Campagnola, Strange & Russell 2010, "A Fast Tour Design Method Using Non-Tangent V-Infinity Leveraging Transfers" (#960)

S. Campagnola (USC), N. J. Strange (JPL) & R. P. Russell (Georgia Tech), "A Fast Tour Design Method Using
Non-Tangent V-Infinity Leveraging Transfers", AAS 10-164, AAS/AIAA Space Flight Mechanics Meeting, 2010.
No DOI (AAS paper).
- Filed as `cyclers_pdf/papers/campagnola-strange-russell-2010-fast-tour-design-non-tangent-v-infinity-leveraging-transfers-enceladus-AAS-10-164.pdf`.
  20 pages, text layer, md5 659fab800cb5f406b6b8e7dbfd90755f.
- I read pp.1-4 and 12-20 from the text layer. I skimmed the VILT derivation (pp.4-11). The tables
  (2-7) are scrambled in the text layer, so I quote only figures stated in the prose.

## 0. Verdict

A linear approximation of the non-tangent v_inf-leveraging transfer (VILT) solution space, plus a graphical
(Tisserand-graph) design method for multi-VILT tours. Demonstrated on a Saturn-system descent to a 200 km
Enceladus orbit: 52 gravity assists at Titan, Rhea, Dione, Tethys and Enceladus, 2.7 yr, total about
445 m/s including Enceladus orbit insertion (abstract), against 3.9 km/s for EOI from a Titan-Enceladus
Hohmann.
- **Not the wanted "Strange, Campagnola & Russell, Leveraging Flybys of Low Mass Moons to Enable an
  Enceladus Orbiter" (AAS 09-435; also cited as AAS 09-208).** That paper is this one's Ref. [10]. Its
  tour is the "second solution" that Table 7 compares against. That item is not in the merged wanted
  list; it appears only in digests.
- **`#943` X1:** no collision. Each leg uses one moon only (VILTs and resonant transfers at Titan, then
  Rhea, and so on), with no repeating two-moon segment. The leg-to-leg transfers are not designed
  ("beyond the scope", p.12).
- **Use for the X1 V2-V3 continuation:** the n:m(sigma) notation, the near-flat VILT solution space and
  the linear approximation are tooling for v_inf changes between resonances (INFERRED applicability).

## 1. Content (READ)

- **VILT (pp.2-4):** a trajectory that starts and ends at the same minor body, with one small impulsive
  maneuver for a large change in v_inf.
  - Notation n:m(sigma): sigma = 0 is a ballistic resonant transfer; sigma = +/-1 is non-resonant, with
    slightly more or fewer than m revolutions (p.4).
  - "Leveraging" was coined by Longuski and first documented by Williams 1990 (footnote, p.1).
- **Low-mass moons (p.1):** leveraging increases the effective bending of a low-mass moon. It "can enable
  reaching a 1:1 resonance when, at best, only a 19:18 resonance could be reached ballistically" (citing
  Ref. [10]).
- **Result (pp.12-13, Figs. 10-18, Tables 2-7):**
  - Legs: Titan, Rhea, Dione, Tethys, Enceladus.
  - Re-solving the exact VILTs agrees with the piecewise-linear values within 3.3 % ("less than 0.02 %
    in most cases").
  - The trajectory starts from the orbit after the pericentre raise maneuver of the Titan Saturn System
    Mission (AAS 09-356).

## 2. Gate relevance

- X1: no repeating Ganymede-Callisto or Ganymede-Europa analogue; a Saturn-system, single-moon-per-leg
  descent.
- Tooling: see the verdict.

## 3. Positive controls

- None transcribed. Tables 2-7 need page-image reads before use.

## 4. Citation mining (references 1-19)

Held: Campagnola & Russell 2010 Endgame Parts 1 and 2 [12], [13] (as AAS 09-224/09-227); Vasile &
Campagnola 2009 [14]; Uphoff et al. [2] is not held (the 1976 JSR version is in the wanted list).

Not held:
1. [10] Strange, Campagnola & Russell (2009), "Leveraging Flybys of Low Mass Moons to Enable an Enceladus
   Orbiter", AAS 09-435. The low-mass-moon leveraging paper. Add to the wanted list (Tier D, X1
   tooling).
2. [8] Sims, Longuski & Staugler (1997), "V-infinity Leveraging for Interplanetary Missions:
   Multiple-Revolution Orbit Techniques", JGCD 20(3):409-415, doi 10.2514/2.4064 (printed).
3. [9] Strange & Sims (2001), "Methods for the design of V-infinity leveraging maneuvers", AAS 01-437.
4. [11] Brinkerhoff & Russell (2009), AAS 09-222; [19] Strange et al. (2009), TSSM mission design, AAS
   09-356; [1] Williams (1990) MS thesis, Purdue; [7] Hollenbeck (1975), AAS 75-087.
