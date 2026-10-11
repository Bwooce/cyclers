# Digest: Campagnola and Kawakatsu 2012, "Three-Dimensional Resonant Hopping Strategies and the Jupiter Magnetospheric Orbiter", JGCD 35(1):340-344 (#960)

Stefano Campagnola (JAXA ISAS, later JPL) and Yasuhiro Kawakatsu (JAXA ISAS). Journal of Guidance, Control, and Dynamics
35(1):340-344, January-February 2012 (Engineering Note). DOI 10.2514/1.53334. Received 6 December 2010, revised 13 June
2011, accepted 13 June 2011. 5 pdf pages, text layer, no OCR. Wanted-list row "Clipper and Jovian tour background"
(Campagnola and Kawakatsu is the first item; Uphoff, Roberts and Friedman 1976 is the second).

- **Filed as** `cyclers_pdf/papers/campagnola-kawakatsu-2012-three-dimensional-resonant-hopping-strategies-jupiter-magnetospheric-orbiter-jgcd-35-340-doi-10.2514-1.53334.pdf`.
  Numbers transcription: `data/sources/campagnola-kawakatsu-2012-numbers.yaml`.
  - The filed file is the received file, unchanged, md5 64795d64336d4f308ce0c40d8828ccf7 (3,719,003 bytes). Pages are printed
    pp.340-344 (PDF page 1 = p.340, page 5 = p.344, which is references only).
  - DOI and pages were confirmed from the first-page header, the PDF subject field and the Crossref lookup of 2026-10-11.
- **No tables.** The Note has no tables. The yaml holds every number printed in the text, the ten equations and the figure labels.
- **How it was read:** all five pages on a 150 dpi render, and the right-hand column with Eqs. 4-7 on a 220 dpi crop.
  - Text numbers: the raw and the layout text extractions and the render agree on every digit.
  - The text layer drops the minus signs, exponents and fractions of every equation. The equations were read from the
    render only and then checked by derivation (sec. 3).
  - The figure labels are drawn graphics with the render as the only witness. Plot readings (Fig. 5, Fig. 6 dots) are approximate.

## 0. Verdict

**A single-moon resonant-hopping design method and a Callisto-only Jupiter Magnetospheric Orbiter (JMO) endgame. It contains
no repeated-moon, cycler or periodic-orbit statement.**
- The method puts the Tisserand parameter first: it fixes the V_inf magnitude (Eq. 2), the largest flyby bend (Eq. 1) and the
  largest reachable inclination (Eq. 6). A pump angle and a crank angle on the V_inf sphere give the resonance and the
  inclination direction.
- A branch-and-bound search over resonant hops at one moon at fixed 100 km altitude gives almost 10,000 solutions in 2 minutes.
  The worked case starts on a 12:1 Callisto resonance and ends on a 6:7 resonance with an inclination above 30 degrees. The
  example solution reaches 48 degrees in 1.5 years.
- **Only Callisto flybys are used in the design.** Ganymede enters only in Figs. 3 and 4 as a second moon to show the graph.
  There is no two-moon sequence, no moon-to-moon transfer and no return to a previous state.

## 1. What the paper does

- **Models (p.340).** A linked-conic (zero sphere-of-influence patched-conic) model, circular moon orbit, ballistic only
  (footnote). The V_inf vector turns by an angle delta (Eq. 1). Its length is conserved through the whole hopping sequence at
  one moon. Eq. 2 and Eq. 3 give |v_inf| = |v_M| sqrt(3 - T) with T the Tisserand parameter in the (r_a, r_p, i) form.
  The paper credits the formula's first use to its ref. [2] (Campagnola and Russell 2010, Endgame Part 2, held as the AAS form).
- **Pump and crank angles (p.340, Fig. 1).** On the V_inf sphere the pump angle alpha in [0, pi] alone sets the semi-major
  axis, so the period and the resonance rho = n:m (n moon revolutions to m spacecraft revolutions). Eq. 4 and Eq. 5 are the
  two directions of that map. The crank angle k in [-pi, pi] is the azimuth. The largest inclination is arcsin(|v_inf|/|v_M|)
  (Eq. 6), reached at alpha = pi/2 + i_max and k = +-pi/2. The paper credits the angles to its refs. [4, 6, 7]. Ref. [6] is
  Uphoff, Roberts and Friedman 1976, the second paper of this wanted-list row.
- **Gravity-assist map (p.340-341, Eqs. 7-9, Fig. 2).** A flyby of bend delta moves (alpha, k) within a lens-shaped region:
  alpha_out in [alpha_in - delta, alpha_in + delta], with k_out from the spherical law of cosines. Special cases are given for
  sin alpha = 0 and for alpha_in near 0 or pi. Fig. 2 draws the (alpha, k) chart for v_inf = 4 km/s at Callisto with the
  r_p = 15 R_J and r_a = 40 R_J level sets and the reach of a 100 km flyby from a 4:1 orbit at k = 30 degrees (the drawn arrow
  gives a 2:1 orbit at i = 15 degrees).
- **The 3-D Tisserand graph (p.341, Eq. 10, Figs. 3 and 4).** Orbits are points (r_a, r_p, i). A flyby is a move along the
  level set of constant T, drawn here as i_max level sets. Fig. 3 draws the i_max = 40 degrees surface for two moons
  (v_inf 7.0 and 5.3 km/s on the labels), with the closed loop reachable from one point by a 100 km flyby. The text says high
  inclination is reached at low energy and high pericenter. Fig. 4 is the i = 0 section, the planar Tisserand graph with i_max
  and resonance level sets for Ganymede and Callisto. The paper calls the graph the linked-conic equivalent of the
  Tisserand-Poincare graph of its ref. [2].
- **Resonant hopping search (p.342).** Inputs: |v_inf|, initial (alpha, k), target alpha and k, the maximum time between
  flybys t_max (so n_max = t_max / P_M), the maximum total time, and optional constraints. At each node, Eqs. 10 and 5 give
  the range of alpha_out and so of rho, and every n:m with n < n_max in that range is a branch. Eqs. 4 and 9 give the new
  state. The sign is chosen to raise the inclination, and alpha_in is the lower bound of alpha_out so that the energy falls.
  A sequence is stored when it meets the target and discarded when the time limit is broken.
- **JMO application (pp.342-343).** After insertion and the pericenter-raise burn the spacecraft is in a 200 day orbit.
  - Start: 12:1 with Callisto, r_p = 16 R_J, v_inf = 6.717 km/s, alpha = 85 degrees, k = 0 (the star in Fig. 4).
  - Target: 6:7 resonance with i > 30 degrees. t_max = 85 days (a little over five Callisto periods). Total at most
    two years. r_p > 12 R_J for radiation.
  - Result: almost 10,000 solutions in 2 minutes. Fig. 5 plots final inclination against flight time. The solutions lie
    between about 1.4 and 2.0 years and between about 42 and 50 degrees. The circled example reaches 48 degrees at 1.5 years.
  - Figs. 6-8 show that example. The pump angle climbs from 85 degrees through the resonance ladder
    12:1, 5:1, 3:1, 2:1, 3:2, 5:4, 1:1, 6:7 (the top-axis labels of Fig. 6 and the point labels of Fig. 7). Fig. 8 labels the
    flybys Ca1 to Ca12.
  - The actual JMO design loops the initial conditions, allows longer transfers and relaxes the pericenter limit. That gives
    several 100,000 solutions in under 12 hours. The details are in ref. [14] (Campagnola and Kawakatsu, JSR, "to be published").

## 2. Tables and figures

- **Tables:** none.
- **Fig. 1 (p.340):** the V_inf sphere with the moon velocity, the spacecraft velocity, i_max, alpha and k.
- **Fig. 2 (p.341):** the (alpha, k) chart at Callisto, v_inf 4 km/s. Resonance lines labelled 4:1, 2:1, 1:1 and 6:7 along the top.
- **Fig. 3 (p.342):** the (r_a, r_p, i) surface for i_max = 40 degrees, two moons, with two closed 100 km reach loops.
- **Fig. 4 (p.342):** the planar section. Two families of vertical-ish resonance lines and two i_max level-set families.
  Upper labels 1:1, 4:3, 2:1, 3:1, 4:1, 5:1, 12:1 and 6:7. Lower labels 2:1, 4:1, 6:1, 8:1, 10:1, 12:1 and 28:1.
  The star lies at r_a near 260 R_J, r_p = 16 R_J, on the 12:1 line and on the 28:1 line.
- **Fig. 5 (p.342):** final inclination against time of flight for the whole solution set.
- **Figs. 6-8 (pp.342-343):** the example in the pump-crank chart, in the (r_a, r_p, i) graph, and in a Jupiter-centred
  inertial view with twelve Callisto flybys.

## 3. Independent checks, all ours (constants are ours, standard values, not from the paper)

- **Initial orbit.** A 12:1 orbit with Callisto has a period of 200.27 days, which matches the "200 day orbit". With r_p = 16 R_J
  its apocenter is 260.0 R_J, which matches the star in Fig. 4.
- **v_inf and alpha.** At i = 0 the Tisserand relation gives v_inf = 6.711 km/s (printed 6.717, difference 0.09 percent, within
  the spread of the moon and Jupiter constants). Eq. 4 then gives alpha = 85.1 degrees (printed 85). Eq. 6 gives i_max = 54.9
  degrees, consistent with "higher than 50 degrees".
- **Fig. 2 example.** alpha = 68 degrees at v_inf = 4 km/s gives rho = 3.999, the printed 4:1.
- **Fig. 3 labels.** v_inf for i_max = 40 degrees is 6.99 km/s at Ganymede and 5.27 km/s at Callisto, the printed 7.0 and 5.3.
- **The resonance ladder of Fig. 6.** At constant v_inf = 6.711 km/s Eq. 4 gives pump angles 85.1, 90.4, 95.3, 100.5, 105.3, 108.9, 114.1
  and 118.4 degrees for 12:1, 5:1, 3:1, 2:1, 3:2, 5:4, 1:1 and 6:7. The Fig. 6 axis has the labelled lines at about 85, 90, 95, 100, 105,
  109, 114 and 118 degrees. So Fig. 6 is a plot of a constant-v_inf sequence, as Eq. 2 requires.
- **The star on two moons.** 12 Callisto periods are 27.99 Ganymede periods, so the same orbit is a 28:1 Ganymede orbit. That is
  why the star sits on the 12:1 line of one family and on the 28:1 line of the other in Fig. 4. This supports reading the upper
  family as Callisto and the lower as Ganymede. The paper does not say which is which on the figure.
- **Flyby strength.** The largest bend of one Callisto flyby at 100 km altitude and v_inf = 6.717 km/s is 6.8 degrees
  (GM 7179 km^3/s^2, radius 2410.3 km). So each flyby moves the pump angle by at most a few degrees and the 85-to-118 degree
  climb needs many flybys. Twelve are drawn (Ca1 to Ca12).
- **No self-consistency problem found.** Every printed number that can be derived from the others was reproduced.

## 4. Relevance to this project

- **Not a cycler paper.** The sequence of resonant orbits is a one-way endgame: energy falls and inclination rises, ending on a
  6:7 orbit. The paper never states that a hop sequence returns to its start, and it states nothing about repeated moons.
  Resonant hopping here means one moon, many resonances, no closure.
- **Against the Ganymede-Callisto two-working-body rows (gc rows in `data/catalogue.yaml`).**
  - The gc rows are two-moon closed cyclers with a 37.57 day period, three Ganymede-Callisto synodic periods, and moon
    V_inf of 1.8 to 3.6 km/s. The JMO hopping is one moon (Callisto) at a constant, much larger V_inf of 6.717 km/s, and its
    orbit starts at a 200 day period with r_a about 260 R_J.
  - So the paper's regime (high V_inf, large pericenter, inclination raising, 3-D) does not overlap the gc rows' regime (low
    V_inf, short period, near-planar). The paper is background for how a single-moon resonant chain is built, not a source for a row.
- **Against the Ganymede-Europa rows.** No Europa content at all. The figures involve Ganymede and Callisto only.
- **Against #973 (two-working-body generator at longer periods).** #973 runs in the ideal Russell-Strange model (circular,
  coplanar). This paper's formulas are the 3-D version of the same linked-conic machinery. The overlap is the method: a
  Tisserand level set per moon, pump and crank angles, and the resonance ratio rho = n:m from the pump angle alone (Eq. 4 and
  Eq. 5). The paper's pump-angle formula for a resonance is an independent check on any planar resonance constructor
  (the ladder of sec. 3 is a ready test: constant v_inf, alpha from rho).
  - Resonances used: 12:1, 5:1, 3:1, 2:1, 3:2, 5:4, 1:1 and 6:7 at Callisto. Inclinations: 0 degrees at the start to 48 degrees.
  - V_inf at each moon: one value per run, 6.717 km/s at Callisto. Other values appear only as figure labels (4 km/s in Fig. 2,
    7.0 and 5.3 km/s in Fig. 3).
  - The paper makes no repeated-moon or cycler statement. It does not restrict the sequence to one moon by a theorem. It
    only treats one moon per search.
- **Method relation to the held Tisserand-graph papers.** The 3-D graph generalises Strange and Longuski 2002 (held) and
  the Tisserand-Poincare graph of Campagnola and Russell 2010 (held as AAS 09-227). Anderson, Campagnola and Lantoine 2016
  (held) uses the planar resonant-orbit view.
- **Literature-gate label the content supports.** No gate anchor exists for this paper (no "Magnetospheric" or "JMO" in
  literature_check.py or the index on 2026-10-11). If one is added: topology {resonant} and {mga-tour} (a non-repeating
  multi-flyby sequence), primary Jupiter, body set {Callisto} with {Ganymede} as a comparison moon only, a single flyby moon.
  Not {repeated-moon} and not {pump-tour} in the Titan sense (the sequence is a resonance ladder at one moon, not V_inf
  leveraging with transfers). A candidate that uses two or more moons in a closed repeating sequence is not covered by it.
- **Risk of misuse.** The paper's "resonant hopping" is not a periodic orbit and its Callisto V_inf is high. The figure
  captions give no body names on the Fig. 3 and Fig. 4 curves, so they must not be read as a Ganymede-Callisto cycler map.

## 5. Citation mining (checked against CORPUS_INDEX.md and the wanted list, 2026-10-11)

- **Held:** Campagnola and Russell 2010, Endgame Parts 1 and 2 (refs 2 and 10; held as the AAS 09-224 and 09-227 forms, same
  authors and content, journal forms JGCD 33(2):463 and 476 not separately held); Strange, Russell and Buffington 2007,
  AAS 07-277 (ref 4); Strange and Longuski 2002 (ref 9); Wolf and Smith 1995 (ref 7); Vasile and Campagnola 2009 (ref 13).
- **Already on the wanted list:** Uphoff, Roberts and Friedman 1976 (ref 6, same row); Labunsky, Papkov and Sukhanov 1998
  (ref 8, in the Outer-planet and Galileo tour-design row).
- **Not in the index and not on the wanted list (new candidates):**
  - **Campagnola, S. and Kawakatsu, Y. (2012), "Jupiter Magnetospheric Orbiter: Trajectory Design in the Jovian System",
    Journal of Spacecraft and Rockets (ref 14, "to be published").** The full JMO design. Medium priority. Volume and pages are
    not confirmed here. This is the one that would give per-flyby V_inf and altitude tables.
  - **Miller, J. K. and Weeks, C. J. (2002), "Application of Tisserand's Criterion to the Design of Gravity Assist
    Trajectories", AIAA 2002-4717 (ref 1).** Low to medium priority, the Tisserand-criterion source.
  - **Kawakatsu, Y. (2009), "V-Infinity Direction Diagram and its Application to Swingby Design", 21st ISSFD, Toulouse (ref 5).**
    Low priority, conference paper.
  - Smith and Buffington 2009 (Cassini Solstice, ref 11) and Campagnola, Takashima and Kawakatsu 2011 (ISAS radiation report,
    ref 15): not cycler-relevant, not added. Land and Doig 1960 (ref 12) and Tisserand 1896 (ref 3): textbook level, not added.
- **Wanted-list action:** in the row "Clipper and Jovian tour background" the Campagnola and Kawakatsu item is marked received.
  The row stays for Uphoff, Roberts and Friedman until that paper is filed.

*Filed as `cyclers_pdf/papers/campagnola-kawakatsu-2012-three-dimensional-resonant-hopping-strategies-jupiter-magnetospheric-orbiter-jgcd-35-340-doi-10.2514-1.53334.pdf`. Numbers transcription: `data/sources/campagnola-kawakatsu-2012-numbers.yaml`.*
