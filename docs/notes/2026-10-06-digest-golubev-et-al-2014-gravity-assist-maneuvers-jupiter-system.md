# Digest: Golubev, Grushevskii, Koryanov & Tuchin 2014, "Gravity Assist Maneuvers of a Spacecraft in Jupiter System" (#960, #943)

Yu. F. Golubev, A. V. Grushevskii, V. V. Koryanov and A. G. Tuchin (Keldysh Institute of Applied
Mathematics), "Gravity Assist Maneuvers of a Spacecraft in Jupiter System", Journal of Computer and
Systems Sciences International 53(3):445-463 (2014), doi 10.1134/S1064230714030083. This is the
translation of Izvestiya RAN, Teoriya i Sistemy Upravleniya 2014(3):149-167.
- Filed as `cyclers_pdf/papers/golubev-grushevskii-koryanov-tuchin-2014-gravity-assist-maneuvers-spacecraft-jupiter-system-jcssi-53-445-doi-10.1134-S1064230714030083.pdf`.
  19 pages, text layer, md5 9f2a02e3b0f7f418a7277f36ffe4c779. Supplied by the owner (batch 19, item 41).
- I read secs. 5-12 and the reference list from the text layer, and skimmed secs. 1-4.
- Wanted-list rank 16 before batch 19 (Tier A, `#943` completeness). Removed from the list in batch 19.

## 0. Verdict (sent to main and twobody-gen2-opus before filing)

**No collision with `#943` gc-1 or gc-2. No repeating G-C pattern. No cycler.**
- The "crossed" ("crisscross") G-C-G maneuvers are one-shot moves. They reduce v_inf for a Ganymede
  landing (the Laplace-P mission concept). They are not a periodic G-C sequence.
- The paper gives no v_inf values, periods or leg times for any G-C-G chain. Nothing in it can be
  compared numerically with gc-1 (k3, 37.57 d, v_inf G 2.397 / C 1.807) or gc-2 (3.617 / 3.039).
- The only published two-working-body G-C cycler known to the `#943` prior-art search is still
  Campagnola 2019 GCGC (see `2026-10-06-943-gc-prior-art-search.md`).

## 1. Content (READ)

- **Mission frame (secs. 1, 5):** Laplace-P (Roscosmos), with possible cooperation with ESA JUICE. The
  phases are JOI, a period-reduction "debut" of resonant flybys, a v_inf-reduction "middlegame", GOI,
  then landing. The debut uses Lambert arcs with period resonance (multiplicity = the smallest integer
  above T_sc / T_Ganymede). A three-step correction impulse aims the flyby at the outer side of the
  moon's orbit, for a decelerating flyby. It is integrated with the satellites' gravity (NAIF ephemeris,
  KIAM ESTK package).
- **Flyby capability table (sec. 3):** dV_max and chi_mod = dV_max / V_pl for the Moon, the Galilean
  moons, Titan and the dwarf planets. The text-layer table is scrambled; read the page image before
  citing any value.
- **Middlegame (sec. 9):** periapsis raises by apocentre burns cost 50-100 m/s each. The authors prefer
  "crisscross" flybys of a second moon. They read the JUICE event table [7] as already doing this: they
  name **G5-C11-G12, G12-C24-G25 and G25-C29-G30** as v_inf-reducing combinations. They cite Galileo [31]
  for the same idea.
- **Secs. 7-8, 10:** Jacobi and Tisserand invariants ("comet invariants") and Tisserand-Poincaré
  diagrams (after Campagnola-Russell Endgame II [23]).
- **Sec. 11 and the "reflection" analogy:** a flyby is a reflection of a beam of trajectories, like
  double refraction in a crystal (Fig. 7). A "bank of virtual maneuvers" (BVM) and "pilot charts" seed
  the search. They say blind Monte Carlo cannot find crisscross maneuvers in the real Jovian system.
- **Sec. 12, results:**
  - A radiation model, with dose accumulated as the search runs (Figs. 9-12).
  - Front growth on the T-P diagram toward Ganymede's (15 R_J, 15 R_J) point (Figs. 13-14). The front
    must stay right of 26 R_J until v_inf is low enough, or Callisto flybys become impossible.
  - Second-level chains in their notation: R11-...-R1n-C12-C21-R1n+1-..., where C12 is a "reflection"
    to the second moon and C21 a "re-reflection" back. The variant
    R11-...-R1n-C12-R21-...-R2k-C21-R1n+1-... allows resonant flybys of the second moon in between.
  - **These are templates for a one-way v_inf descent, not closed periodic cycles.**
  - The VDPP interactive tool (Figs. 16-17).

## 2. Relevance

- `#943`: a negative prior-art result, recorded in the gc note. Golubev's C12-C21 notation is the
  nearest named G-C-G pattern in the Russian literature, and it is non-periodic.
- `#960`: it shows that the JUICE G-C-G triplets were read as v_inf-reduction devices by 2014. This agrees
  with the Boutonnet 2024 JUICE digest.

## 3. Citation mining (references 1-41)

Held:
- Campagnola & Russell, Endgame I/II [22, 23]: held as the AAS 09-224 / 09-227 conference forms.
- Strange, Russell & Buffington 2007, AAS 07-277 [24].
- Szebehely 1967 [28].
- Minovitch TR 32-464 [2]: NOT held; it is wanted-list rank 1.

Not held. These are added to the wanted list (Tier D, Jovian tour-design row) unless noted:
- Strange & Longuski 2002, "Graphical Method for Gravity-Assist Trajectory Design", JSR 39(1):9-16,
  doi 10.2514/2.3800 (Crossref CONFIRMED) [33]. It is cited across the project, but no corpus file
  holds it.
- Woolley 2010, "Endgame strategies for planetary moon orbiters", PhD thesis, CU Boulder [19].
- Campagnola, Skerritt & Russell 2011, "Flybys in the planar, circular, restricted, three-body
  problem", AAS 11-245 [25].
- Campagnola, Boutonnet, Schoenmaekers, Grebow, Petropoulos & Russell 2012, "Tisserand-leveraging
  transfers", AAS 12-185 [26].
- Boutonnet & Schoenmaekers 2012, "Mission analysis for the JUICE mission", AAS 12-207 [6]. The JUICE
  CReMA reports [7, 8] are ESA-internal and are not listed.
- Borovin, Golubev, Grushevskii, Koryanov & Tuchin 2013, KIAM Preprint 72/2013 [5] (Keldysh library,
  free, Russian).
- Labunsky, Papkov & Sukhanov 1998, "Multiple Gravity Assist Interplanetary Trajectories", Gordon &
  Breach [32].
- Uphoff, Roberts & Friedman 1976 [31]: already in wanted row "Clipper and Jovian tour background".

Not listed: refs. 1, 3, 4, 9-18, 20, 21, 27, 29, 30, 34-41 (textbooks, KIAM internal reports,
presentations, radiation and software references).
