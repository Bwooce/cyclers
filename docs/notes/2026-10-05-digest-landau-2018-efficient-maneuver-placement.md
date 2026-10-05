# Digest: Landau 2018, "Efficient Maneuver Placement for Automated Trajectory Design" (#960)

D. F. Landau (JPL), "Efficient Maneuver Placement for Automated Trajectory Design", J. Guidance, Control,
and Dynamics 41:1531-1541 (2018), doi 10.2514/1.G003172.
- Crossref-confirmed in batch 2 as the journal version of AAS 15-585.
- Filed as `cyclers_pdf/papers/landau-2018-efficient-maneuver-placement-automated-trajectory-design-JGCD-41-1531-doi-10.2514-1.G003172.pdf`.
  11-page "Article in Advance", text layer, md5 845c1a7258624f04d36f857cc93e4afb.
- I read the full text layer. Fig. 4 (p.7) was read as a page image.

## 0. Verdict for `#943` X1 (literal-collision check)

- No collision.
- The paper is a maneuver-placement method: primer-vector placement of deep-space maneuvers in broad
  searches, reducing the maneuver search to one dimension.
- It contains NO cycler and NO repeated Ganymede-Callisto or Ganymede-Europa sequence.
- The name "STAR", which Campagnola et al. 2019 cite for the GCGC first guesses, does not appear in this
  paper. It describes the method, not the software.
- Its one Jovian example is a one-way Europa endgame: the "Banzai Pipeline" of McElrath, Campagnola &
  Strange 2012, built from 180-deg transfers. Fig. 4, read off the figure:
  - Callisto 23 Sep 2027, v_inf 1.240 km/s, declination -85.6 deg.
  - Callisto 1 Oct, 1.222.
  - Ganymede 7 Oct, 1.467.
  - Ganymede 10 Oct, 1.459.
  - Europa 13 Oct, 1.643.
  - DSMs: 15 Oct (130 m/s) and 3 Nov (50 m/s).
  - Europa 1 Nov, 0.978.
  - Europa 26 Nov, 0.665, declination -0.5 deg.
  - Total 870 m/s leveraging plus capture, over 64 days.
  - This is a C-C-G-G-E-E-E pump-down, not a repeating cycle.

## 1. Content (READ)

- Sect. II: primer-vector theory. The maneuver's first-order effect on any objective J is computed
  explicitly from the primer vector along an existing transfer.
- The optimal location and direction come from the primer extrema. The magnitude is a search variable.
- A Lambert-like flight-time root solve then fixes the two-conic transfer, Eqs. (26)-(44).
- Appendix: a quasi-second-order root iteration.
- Examples (Sect. IV):
  - (A) The Europa endgame above: v_inf leveraging at apoapsis, confirmed optimal (Fig. 2, near a 6:5
    resonance).
  - (B) A Saturn tour from a Titan 2:1 resonance through Rhea, Dione and Tethys to capture at Enceladus.
    The objective is the period change per flyby. The minimum-Delta-V tour has 54 flybys (Fig. 7, Pareto
    front).
  - (C) Earth-Mars broken-plane transfers under a declination-penalised launch vehicle.

## 2. Citation mining (references [1]-[16])

Held:
- [16] Campagnola, Strange & Russell 2010 (held as the AAS 09-227 or related preprint; see the
  CORPUS_INDEX endgame rows).
- [13] Battin 1999: textbook, not held.

Not held, in priority order:
1. [14] McElrath, T. P., Campagnola, S. & Strange, N. (2012), "Riding the Banzai Pipeline at Jupiter:
   Balancing Low Delta-V and Low Radiation to Reach Europa", AIAA 2012-4809, doi 10.2514/6.2012-4809
   (printed in the reference list). The Callisto/Ganymede pi-transfer endgame.
2. [15] Strange, Campagnola & Russell 2009, AAS 09-435 (also in the Campagnola 2019 list as AAS 09-208; the
   AAS number differs between the two citing papers).
3. [3] Lantukh, Russell & Campagnola (2015), JSR 52(3):697-710, doi 10.2514/1.A32918.
4. [2] Landau, Lam & Strange (2009), AAS 09-428 (SEP broad search).
5. Method background, low priority: [1] Izzo et al. 2007; [4] Vasile & Ceriotti 2010; [5] Gad &
   Abdelkhalik 2011; [6] Englander et al. 2012; [7], [8] Lawden; [9] Fernandes 1999; [10] Carter 1990;
   [11] Jezewski & Rozendaal 1968; [12] Prussing 2010.
