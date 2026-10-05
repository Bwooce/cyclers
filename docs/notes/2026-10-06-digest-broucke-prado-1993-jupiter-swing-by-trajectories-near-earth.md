# Digest: Broucke & Prado 1993, "Jupiter Swing-By Trajectories Passing Near the Earth" (#960)

R. A. Broucke (University of Texas at Austin) & A. F. B. A. Prado (UT Austin and INPE), "Jupiter Swing-By
Trajectories Passing Near the Earth", AAS 93-177, AAS/AIAA Spaceflight Mechanics Meeting, 1993. No DOI
(AAS paper).
- Filed as `cyclers_pdf/papers/broucke-prado-1993-jupiter-swing-by-trajectories-passing-near-earth-AAS-93-177.pdf`.
  18 pages, text layer, md5 9d629a62f7a1c5abbe127e294b9273da.
- I read pp.1-11 and 16-18 from the text layer. I skimmed Appendices A-C and did not use the Fig. 2
  letter-plots or the Table 3 values.

## 0. Verdict

A classification of single Jupiter swing-bys in the planar circular restricted problem (Sun-Jupiter,
Lemaitre-regularised). It marks which trajectories cross the Earth's orbit before and/or after the
encounter. Nothing repeats: no periodic or cycler geometry, and no repeating Jupiter-Earth path. It is a
continuation of Broucke's AIAA 88-4220 ("The Celestial Mechanics of Gravity Assist").
- **Not the conference form of the wanted "de Almeida Prado & Broucke 1995/1996" items.** Those are, by
  Crossref, "Transfer Orbits in Restricted Problem" (JGCD 18:593-598, 1995, doi 10.2514/3.21428) and
  "Transfer Orbits in the Earth-Moon System Using a Regularized Model" (JGCD 19:929-933, 1996, doi
  10.2514/3.21720). Different titles and topics, so those items stay unmarked.
- **`#942`/`#943`:** no collision; background on swing-by classification only.

## 1. Content (READ)

- **Model (pp.4-6):** planar circular restricted three-body problem, Sun-Jupiter, canonical units, with
  Lemaitre's regularisation near Jupiter.
  - Table 1 (p.5): unit of distance 778,000,000 km, unit of time 689.567 days, unit of velocity
    13.058 km/s.
  - My check: with GM_Sun = 1.327e11 km^3/s^2, sqrt(a^3/GM) = 5.957e7 s = 689.5 d, and a/t =
    13.06 km/s. Both agree.
- **Parameters (pp.3-4):** Jacobi constant J (equivalent to v_inf), the periapsis angle psi from the
  Sun-Jupiter line, and the periapsis distance Rp. Integrate forward and backward from periapsis to
  0.5 canonical units from Jupiter. Then classify the orbit before and after as direct or retrograde,
  elliptic or hyperbolic.
- **Classes (Table 2, p.8):** 16 letters A-P (before x after), plus Z for temporary capture. Upper case:
  never crosses the Earth's orbit. Lower case: crosses it in one time direction. Bold lower case: crosses
  it in both directions (a possible Earth-Jupiter-Earth path "if a proper timing condition can be
  found").
- **Symmetry (p.8):** psi and psi + 180 deg are time-reverses: I<->C, J<->G, L<->O, B<->E, N<->H, M<->D;
  A, F, K and P are unchanged.
- **Results (pp.8-9, Fig. 2):**
  - Plots for Rp = 1.1, 1.5, 2, 5, 10 and 50 Jupiter radii, with psi from 180 to 360 deg and J from
    -1.45 to 1.55.
  - They confirm energy loss for a flyby in front of Jupiter (psi from 0 to 180 deg, maximum at 90 deg)
    and energy gain for a flyby behind it (maximum at 270 deg).
- **History (Appendices A-C, pp.13-16):** d'Alembert, Laplace (sphere of influence, comet of 1770),
  H. A. Newton 1878/1893 (cometary capture by Jupiter), and von Pirquet 1928 (Jupiter-flyby vector
  addition).
  - Also Minovitch's 1961 JPL memorandum 312-130. The corpus holds it, in the Minovitch 1961 letters and
    documents file.

## 2. Gate relevance

- None for R1 or X1. It is background for any regularised close-approach code (Lemaitre) and the
  letter classification.

## 3. Positive controls

- Table 1 canonical units (checked above). The Table 3 example trajectories (types N, j, b) are not
  transcribed; read them on the page image before use.

## 4. Citation mining (references 1-42)

Held: Minovitch 1961 JPL TM 312-130 [18] (in the Minovitch 1961 letters file); Uphoff 1989 [33] is not
held.

Not held (selected):
1. [16] Broucke, R. A. (1988), "The Celestial Mechanics of Gravity Assist", AIAA 88-4220. The base paper
   of this one. Not in the wanted list; low priority.
2. [19], [20] Dowling, Kosmann, Minovitch & Ridenoure (1990, 1991), history of gravity propulsion (IAF).
3. [21] Flandro (1966), the Grand Tour paper.
4. [36]-[42] Leverrier 1847, Callandreau 1890-1902, van Woerkom 1948, Newton 1893, Everhart 1969 (AJ
   74:735): cometary capture. Background.
