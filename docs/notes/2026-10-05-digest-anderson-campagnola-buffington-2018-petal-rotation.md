# Digest: Anderson, Campagnola & Buffington 2018, "Analysis of Petal Rotation Trajectory Characteristics" (#960)

R. L. Anderson, S. Campagnola and B. B. Buffington (JPL), "Analysis of Petal Rotation Trajectory
Characteristics", J. Guidance, Control, and Dynamics 41(4):827-840 (April 2018), doi 10.2514/1.G002571.
- Crossref-confirmed in batch 2. Presented as AIAA 2014-4350.
- Filed as `cyclers_pdf/papers/anderson-campagnola-buffington-2018-petal-rotation-trajectory-characteristics-JGCD-41-4-827-doi-10.2514-1.G002571.pdf`.
  14 pages, text layer, md5 72a204557529bcae6d9cf7e959190390.
- I read the text layer in full.

## 0. Verdict

- **Petal rotations are single-moon objects.** They alternate long and short nonresonant transfers with ONE
  gravity-assist body, n:n+ / m:m- (for example 1:1+/2:2-), to rotate the line of apsides.
- In the CR3BP, "petal rotation trajectories ... are unstable periodic orbits" (abstract).
- The paper converges them by multiple shooting and continues them for Io, Europa, Ganymede, Callisto,
  Titan and Triton.
- The Jupiter-Europa 1:1+/2:2- family is identified as Henon's period-three family g3, which intersects
  family f (the DROs) twice (pp.834-835, Conclusions).
- **No collision with `#943` X1 or `#945` R2.**
  - They are single-moon periodic orbits: one moon, repeated flybys of the same moon.
  - They are not two- or three-moon cyclers.
  - The cycler link is only a citation: petals are used in "cycler trajectory concepts [13]" (Russell &
    Strange 2009).
- **No numeric periodic-orbit initial conditions are printed.** The families are figures only (Figs. 9-24).
- Usable numbers are the mass ratios (Table 1) and the body constants (Table 2).
- Relevance: these single-moon n:n+/m:m- periodic orbits are what Campagnola 2019's Callisto petals
  "shadow". They are prior art for any single-moon resonant or petal periodic orbit the project computes at
  these moons, and the period-3 g3 identification is a family-identity control.

## 1. Content (READ)

- Patched-conic petal analysis (Sect. II):
  - Resonant and nonresonant transfers. Notation m:n with m spacecraft and n moon revolutions; the sign
    shows more or less than m inertial revolutions (p.828).
  - The apse rotation per flyby pair is Eq. (19).
  - Feasibility is limited by V_inf,max for a minimum-altitude flyby (Eq. (18), Table 2).
  - Fig. 4: rotation rate per secondary revolution against dimensionless V_inf for many n:n/m:m pairs and
    systems.
  - Figs. 5-6: dimensional results for Europa, Ganymede, Callisto and Titan.
- CR3BP (Sects. III-V):
  - Patched-conic initial guesses at the flybys are corrected by multiple shooting into periodic orbits in
    the rotating frame (Figs. 8-9).
  - Examples:
    - Jupiter-Europa 1:4+/1:5- (far from Europa; patched conics match well).
    - Jupiter-Europa 1:1+/2:2- (continued in C and V_inf; it becomes period-3 orbits near family f, and is
      family g3 per Henon).
    - Saturn-Titan.
    - Io 1:1+/2:2-, Ganymede 2:2+/2:2-, Callisto 2:2+/3:3- (Fig. 23).
  - Patched conics fail at low V_inf: below about 3 km/s at Europa and about 1.5 km/s at Titan and Callisto
    (Conclusions; Fig. 24).
- **Table 1, mass ratios (p.829):**

| System | mu |
|---|---|
| Jupiter-Io | 4.705093e-5 |
| Jupiter-Europa | 2.526645e-5 |
| Jupiter-Ganymede | 7.803691e-5 |
| Jupiter-Callisto | 5.667999e-5 |
| Saturn-Titan | 2.365805e-4 |
| Neptune-Triton | 2.087757e-4 |

## 2. Citation mining (references [1]-[45])

Held:
- [13] Russell & Strange 2009.
- [18] Campagnola & Russell 2010, endgame part 1 (held as AAS 09-224).
- [41] Campagnola & Russell 2010, part 2 (held as AAS 09-227).
- [19] Anderson & Lo 2011, JAS 58:167 (held).
- [3] Wolf & Smith 1995.
- [9] Campagnola et al. 2015 Neptune-Triton.

Not held, in priority order:
1. [45] Henon, M. (2003), "New Families of Periodic Orbits in Hill's Problem of Three Bodies", CMDA
   85(3):223-246, doi 10.1023/A:1022518422926 (printed in the reference list). [44] Henon, M. (1970),
   "Numerical Exploration of the Restricted Problem. VI. Hill's Case: Non-Periodic Orbits", A&A 9(1):24-36.
   These are the family g3 sources, the family-identity control for 1:1+/2:2-.
2. [42] Stromgren (1933), Bull. Astron. 9:87-130 (family f); [43] Lam & Whiffen (2005), AAS 05-110
   (Europa DROs).
3. [21] Anderson, Campagnola & Lantoine 2016, CMDA 124:177, doi 10.1007/s10569-015-9659-7.
4. [16] Barrabes & Gomez, spatial p-q resonant orbits; [17], [26] Anderson, "Approaching Moons from
   Resonance via Invariant Manifolds"; [22]-[25] Anderson & Lo papers.
5. Mission background, low priority: [1]-[12] Cassini, Europa, Neptune and JUICE tour papers; [27] Uphoff
   et al.; [28] Strange & Sims AAS 01-437; [29] Strange, Campagnola & Russell.
