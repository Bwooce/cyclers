# Digest: Boutonnet, Langevin & Erd 2024, "Designing the JUICE Trajectory" (#960, #943)

A. Boutonnet (ESA/ESOC), Y. Langevin (IAS) & C. Erd (ESA/ESTEC), "Designing the JUICE Trajectory", Space
Science Reviews 220:67 (2024), doi 10.1007/s11214-024-01093-y. Open access (CC BY 4.0).
- Filed as `cyclers_pdf/papers/boutonnet-langevin-erd-2024-designing-juice-trajectory-space-sci-rev-220-67-doi-10.1007-s11214-024-01093-y.pdf`.
  71 pages, text layer, md5 19fda97e4d4fd19b7ee64ac00c575c7b.
- Scope (book-length paper): I read the abstract, sec. 3.5 (pp.33-44, the Ganymede-approach strategies
  and the Callisto-Ganymede-Callisto round trip) and the conclusions. I searched the full text for
  cycler, repeating and Ganymede-Callisto content. Secs. 1-2 and 4-5 were skimmed. The numbers below
  are from the text layer of the pages named.

## 0. Verdict

The JUICE tour contains no cycler and no repeating Ganymede-Callisto segment. "Cycler" does not appear in
the paper. Its Ganymede-approach strategy uses a **single "Callisto-Ganymede-Callisto round trip"**
(secs. 3.5.1-3.5.2) in which both moons bend the trajectory. This is a one-shot v_inf-reduction step,
not a periodic orbit.
- **`#943` cell gc:** no literal collision with gc-1 or gc-2. It is the nearest published one-shot
  relative by v_inf (see sec. 2).

## 1. Content (READ)

- **Abstract and conclusions:**
  - launch 14 April 2023; Jupiter arrival July 2031, near equinox, so eclipses are a major constraint;
  - a Moon-Earth double swing-by;
  - a low-energy, Callisto-assisted endgame to Ganymede orbit insertion (GOI);
  - the high-inclination phase was extended to 400 days, with a maximum inclination of 33.1 deg;
  - post-launch Delta-V margin 150 m/s.
- **Sec. 3.5.1, velocity-reduction strategy (pp.33-34, Figs. 26-27):**
  - From a Callisto v_inf of 5.22 km/s, two Callisto flybys raise perijove to about 15 RJ (a grazing
    Ganymede orbit, period 27.4 d). The Ganymede v_inf is about 3.3 km/s.
  - "The most effective way to further reduce the infinite velocity to Ganymede consists in implementing
    Callisto - Ganymede - Callisto round trips": Ganymede flybys lower the apojove to a near-grazing
    Callisto encounter at 26.5 RJ, and Callisto raises the perijove back.
  - In the optimum case the Ganymede v_inf drops to 1.5 km/s and GOI to about 0.75 km/s.
  - Several round trips would step the Callisto v_inf down (5.22 to 3.5 to 2 km/s) at a cost in time
    and dose.
- **Sec. 3.5.2, the baseline's single round trip (pp.37-40, Figs. 29-32):**
  - A near-tight Callisto flyby (about 300 km) on 2 November 2033 starts a direct transfer to Ganymede,
    v_inf about 3.35 km/s (inbound) or 3.3 km/s (outbound).
  - Two departure windows from Ganymede exist "over the 50-day cycle": 7 or 8 Ganymede periods after
    arrival. Ganymede "petal" strategies rotate the line of apses.
  - The return to Callisto comes 7.6 d (second window) or 24.3 d (first window) after departure. The
    Callisto arrival v_inf is about 1.9 km/s (orbit period 10.93 d) to about 2.3 km/s (11.4 d).
  - Callisto petals follow. The round trip ends at a Callisto outbound encounter in early May 2034, with
    v_inf 1.8-2.4 km/s and right ascension 160-250 deg, which starts the Callisto-assisted endgame.
  - The total round trip plus the Callisto position adjustment takes about 6 months.

## 2. Gate relevance (`#943` cell gc)

- **gc-1** (G/C 2.397/1.807 km/s, 37.57 d): the round trip's Callisto v_inf (1.8-2.4 km/s) overlaps
  gc-1's Callisto value. Its Ganymede v_inf (about 3.3 km/s) is 0.9 km/s higher. Not periodic: the
  Callisto v_inf changes from 5.22 to 1.8-2.4 over the sequence.
- **gc-2** (3.617/3.039): Ganymede v_inf within 0.3 km/s; Callisto about 0.7-1.1 km/s lower.
- The 50-day Ganymede-Callisto window cycle (about 4 G-C synodic periods) is the same near-repeat that
  Liang 2024 uses. gc-1 and gc-2 have k = 3 (37.57 d).
- Cite JUICE as an operational precursor in which both moons bend in a C-G-C sequence; it does not
  publish a cycler.

## 3. Positive controls

- None suitable (the dates and v_inf values are approximate, from figures and prose).

## 4. Citation mining

Not mined in full. Relevant refs: Boutonnet & Schoenmaekers 2012 (AAS 12-207, JUICE mission analysis;
already cited in the Buffington 2012 digest); Campagnola's low-energy endgame (held, AAS 09-224/09-227);
Heaton et al. 2002 (Europa Orbiter tour).
