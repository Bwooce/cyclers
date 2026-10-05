# Digest: Lantukh, Russell & Campagnola 2015, "V-Infinity Leveraging Boundary-Value Problem and Application in Spacecraft Trajectory Design" (#960)

D. V. Lantukh, R. P. Russell (UT Austin) & S. Campagnola (JAXA), "V-Infinity Leveraging Boundary-Value
Problem and Application in Spacecraft Trajectory Design", J. Spacecraft and Rockets 52(3):697-710
(May-June 2015), doi 10.2514/1.A32918 (wanted-list CONFIRMED; printed on every page).
- Filed as `cyclers_pdf/papers/lantukh-russell-campagnola-2015-v-infinity-leveraging-boundary-value-problem-jsr-52-3-697-doi-10.2514-1.A32918.pdf`.
  14 pages, text layer, md5 81235e761eb10d32c721e7dae2fecc54.
- I read pp.1-2 and 12-14 from the text layer. I skimmed secs. II-VI (formulation, time-of-flight
  function, implementation, tangent-VILT design space, examples).

## 0. Verdict

A Lambert-like boundary-value formulation of the v_inf-leveraging transfer: given the two boundary
positions and the time of flight, one extra continuous parameter (the maneuver location), solved by a
1-D root solve. It handles eccentric body orbits, ephemeris positions, non-tangent maneuvers and
**inter-body** leveraging (different departure and arrival bodies). Bi-elliptic transfers are a special
case of inter-body leveraging.
- **Use for the X1 V2-V3 continuation:** the most directly usable VILT tool in the corpus. It fits a
  Lambert-based outer loop with "a similar computational burden ... as that of the conventional
  multirevolution Lambert problem" (Conclusions, sec. VIII, pp.709-710) (INFERRED fit to our generator).
- **`#943` X1:** no collision. The examples are single transfers: Delta-V-EGA characterisation, VILT
  families with an eccentric flyby body (Mercury), inter-body families, and a Jupiter capture. No
  repeating moon-tour segment.

## 1. Content (READ)

- **Formulation (secs. II-IV):** the same boundary conditions as Lambert's problem plus one free
  maneuver parameter. Algorithm 2 is the numerical solution. Sec. IV explores and bounds the number of
  solutions of the time-of-flight function.
- **Examples (secs. VI-VII):**
  - Earth Delta-V-EGA.
  - Mercury (eccentric orbit) VILT families; Table 2 gives the planet elements and Table 3 the search
    constraints.
  - Inter-body families.
  - A Jupiter capture from the GTOC6 problem: infinity, Io, Ganymede, VILM, Ganymede (Fig. 17, Table 4
    moon elements at MJD 58849).
- **Table 5 (p.709, read on the page image; it matches the text layer):** v_inf1 14.1617 km/s, v_inf2
  6.16534 km/s, TOF 108.534 d, Delta-V 0.401698 km/s, xi_L -0.761599 rad, efficiency 19.9064. "There
  are no tangent VILTs that can generate the same sequence."

## 2. Gate relevance

- X1: tooling only; see the verdict.

## 3. Positive controls

- Table 5 (Jupiter capture VILT), page-image verified.
- The Table 4 Io and Ganymede elements (p.708, text layer; re-read on the image before use) give the setup: circular, a = 422,030 and 1,070,590 km, GM 5959.92
  and 9887.83 km^3/s^2, mean anomalies 4.279603 and 2.215313 rad at MJD 58849.

## 4. Citation mining

Not mined in full; the references overlap the Campagnola 2010 and Russell-Ocampo 2005 lists. New and
not held: Sims, Longuski & Staugler 1997 (JGCD 20(3):409, doi 10.2514/2.4064), as in the Campagnola
2010 digest.
