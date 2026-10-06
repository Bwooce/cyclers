# Digest: Anderson, Campagnola & Lantoine 2016, "Broad search for unstable resonant orbits in the planar circular restricted three-body problem" (#960)

R. L. Anderson (JPL), S. Campagnola (JAXA) & G. Lantoine (JPL), "Broad search for unstable resonant orbits
in the planar circular restricted three-body problem", Celestial Mechanics and Dynamical Astronomy
124:177-199 (2016), doi 10.1007/s10569-015-9659-7 (printed). Received 14 May 2014; published online
11 December 2015. Conference precursor: AAS/AIAA Astrodynamics Specialist Conference, 2014.
- Filed as `cyclers_pdf/papers/anderson-campagnola-lantoine-2016-broad-search-unstable-resonant-orbits-pcrtbp-cmda-124-177-doi-10.1007-s10569-015-9659-7.pdf`.
  23 pages, text layer, md5 e4a1896562cf3ba2204227ba232d391b. Supplied by the owner.
- I read the abstract, secs. 1-3, 5.1-5.2, 7.2 and 8 from the text layer, and skimmed secs. 4 and 6.
  Orbit data are given only in figures.

## 0. Verdict

A methods paper on computing unstable resonant periodic orbits in the planar CR3BP. It compares:
- continuation from two-body orbits;
- a single-shooting grid search;
- "flyby maps" in the CWIC model (Campagnola, Skerritt & Russell 2012, CMDA: Keplerian arcs with CR3BP
  integration only near the secondary);
- continuation.
The focus is multi-loop resonant orbits in Jupiter-Europa (3:4, 4:5, 5:6, 6:7, 7:8, ...), with some
interior resonances in Jupiter-Ganymede.
- **No cyclers, and no two-moon orbits:** single-moon resonant-orbit families only. No `#943` collision.
- **Use:** an enumeration and positive-control source for the single-moon resonant building blocks of
  Jovian tours, and for the X1 resonant legs (INFERRED).

## 1. Content (READ)

- **Table 1 (p.180), mass ratios:** Jupiter-Europa mu = 0.0000252664488504; Jupiter-Ganymede mu =
  0.0000780369094055.
- **Sec. 2.3:** unstable resonant orbits, and their role in the Europa approach (the 3:4-5:6 sequence
  of Johannesen & D'Amario 1999).
- **Sec. 3:** continuation from two-body orbits. Fig. 4 shows interior resonant orbits in
  Jupiter-Ganymede at C = 2.99, including 3:4.
- **Sec. 4:** single-shooting grid search. It finds unique orbits across many resonances, but is
  limited to orbits with few intersections.
- **Secs. 5-6:** CWIC flyby maps. A p:q orbit is integrated only between apocentres x_a and x_b near
  the secondary, with p-1 Keplerian revolutions of Tisserand constant T approximately equal to C
  (Eq. 11). A grid over (T_i, q-tilde_i) is solved for theta_f = -pi. This computes 4:5 families with
  multiple intersections.
- **Sec. 7:** continuation of CWIC solutions, and of a typical 3:4 orbit. Unstable resonant orbits can
  be reached from stable ones in particular cases. The stable and unstable orbits sit near 4:5 and 3:4
  respectively (Poincare-section discussion, sec. 7.2).

## 2. Gate relevance

- X1: background for resonant legs. It does not bear on any collision check.
- The Table 1 mass ratios are a sourced value for Jupiter-Europa and Jupiter-Ganymede mu (check against
  the code constants before use, per the digest-not-adoption rule).

## 3. Positive controls

- Table 1 mu values only. The orbit families are in figures (no initial-condition tables), so they need
  digitising before use.

## 4. Citation mining

Held: Anderson & Lo 2010, 2011; Campagnola & Russell 2010 (as the AAS 09-224/227 versions); Campagnola,
Buffington & Petropoulos 2014; Kumar, Anderson & de la Llave 2023 (cites this paper). Not held (low priority): Anderson 2012 AAS 12-136 / 2015 JGCD
38(6):1097 ("Approaching Moons from resonance via invariant manifolds"); Anderson 2013 AAS 13-493;
Campagnola, Skerritt & Russell 2012 (CMDA, "Flybys in the planar, circular, restricted, three-body
problem", the CWIC model; also cited by McElrath 2012); Campagnola et al. (Tisserand-Poincare graph);
the other Anderson & Lo invariant-manifold papers; Lo & Parker (unstable resonant orbits near Earth). Not added.
