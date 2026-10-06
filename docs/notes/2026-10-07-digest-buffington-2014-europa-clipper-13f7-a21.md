# Digest: Buffington 2014, "Trajectory Design for the Europa Clipper Mission Concept" (#960 batch 30; #943)

B. Buffington (JPL, Europa Clipper pre-project mission design lead), AIAA/AAS Astrodynamics Specialist
Conference, AIAA SPACE Forum, San Diego CA, 4-7 August 2014, AIAA 2014-4105, doi 10.2514/6.2014-4105.
Crossref title: "Trajectory Design Concept for the Proposed Europa Clipper Mission"; the PDF metadata title
is the same with "(Invited)". The title printed on p.1 is "Trajectory Design for the Europa Clipper Mission
Concept".
- Source file: supplied file `buffington2014.pdf` (not modified), 17 pp., text layer (Word to
  PDF; AIAA ARC download stamp), md5 d97c43417c5f5f3847b2b0dcdf1df79c.
- Filed as `cyclers_pdf/papers/buffington-2014-trajectory-design-europa-clipper-mission-concept-aiaa-2014-4105-doi-10.2514-6.2014-4105.pdf`.
- How I read it: full text layer. Table 3 (p.7) is an image. I read it at 220 dpi on the page image
  (all 59 rows; the columns used below). Tables 1, 2, 4, 5 and 6 are images. I did not transcribe them. No
  number below depends on them, except the abstract and text values (45/5/9 flybys, 20.1 deg, 164 m/s,
  2.8 Mrad), which are in the text layer.
- Wanted-list row 45 (Clipper and Jovian tour background; this DOI is marked CONFIRMED there).
- Scripts: `buffington_vca_check.py` with `buffington_vca_check_output.txt`.

## 0. Verdict

**This conference paper is the first full published description of the 13F7-A21 Clipper tour (INFERRED; the 2012 paper covers 11F5-A21): 45 Europa, 5 Ganymede and 9
Callisto flybys over 3.5 years. It is a mission-design paper. It has no cycler and no repeating
two-moon path. No collision with `#943` gc-1, gc-2 or ge-1..3.** The filer's collision verdict is
confirmed on the Table 3 image. Details:
- **The only G-C alternation** is the one-way pump-down G1 G2 G3 C1 G4 C2 (Table 3, "Transition to
  Europa Science", 22 Oct 2028 to 16 Mar 2029). V-infinity on the image: G1 6.39, G2 6.44, G3 6.43,
  C1 5.54, G4 5.24, C2 4.35 km/s. So the Ganymede V-infinity falls from 6.39 to 5.24 and the Callisto
  V-infinity from 5.54 to 4.35. Transfers are non-resonant from G3 on ("N R" in the m, n columns); G1
  and G2 are printed 8:1 and 5:1. The G2 "5" is a likely misprint for 4 (section 2). The V-infinity falls at every pass, so this is energy pumping, not a cycler.
- **The E-C switch-flip is used once:** E28 (27 May 2030), then C3, C4, C5, C6, C7, C8, C9 (1 Jun to
  7 Nov 2030), then E29 (9 Dec 2030). Callisto V-infinity 2.80-2.97 km/s. The m:n entries are 3:4,
  1:1, "pi-trans", 1:1, 2:3, 2:3, then N R to Europa. The text (p.11) says this phase uses "seven Callisto and five
  Europa flybys". The count does not reconcile exactly: the text describes E25-E30 (six), while Table 3
  groups E25-E28 (four) with the transition and E29-E30 under "pump up". There is no Ganymede flyby in it and no
  return to Callisto later.
- **Cyclers are named only as an unused option** (p.11): the illuminated-hemisphere change "can be
  accomplished by petal rotation via Europa ..., Ganymede, or Callisto, cyclers[29], a Europa
  pi-transfer, or a 'switch-flip'". Ref. 29 is Russell & Strange 2009, JGCD 32(1):143-157 (HELD). A
  note on the wording: as printed it is a list. "Cyclers" may be a separate item after "petal rotation
  via Europa, Ganymede, or Callisto", rather than "Ganymede or Callisto cyclers". Either way, no
  cycler is flown. 13F7-A21 uses the E-C switch-flip.
- Ganymede is used only in the pump-down (G0-G4) and for end-of-mission disposal (impact, baseline
  plan, p.11). So there is no G-E chain.
- Catalogue implication: none. PROPOSAL: cite this paper, with Lam et al. 2015 Table A1, as the
  published source of the 13F7-A21 G-C pump-down if a gc result ever needs a "precursor G-C-G-C
  chain" citation. **Correction for the Lam 2015 digest (section 2): the speeds quoted there as v_inf
  are closest-approach speeds, not V-infinity.**

## 1. Content (READ)

- Mission concept (pp.1-4): a multiple-flyby architecture instead of a Europa orbiter. The spacecraft
  dips into the radiation belt for each Europa pass, then spends about 70% of a roughly 14-day orbit
  outside it to downlink ("store and forward"). Model payload of 8 instruments plus gravity science.
  Global-regional coverage means at least 11 of 14 Europa panels.
- Interplanetary (p.5): Atlas V 551, VEEGA, 21-day launch period 15 Nov to 5 Dec 2021, max C3 15 km^2/s^2,
  Jupiter arrival April 2028 held fixed. An SLS Earth-Jupiter direct 2022 option is shown (Fig. 4).
- Tour naming (p.6): "13F7-A21" = designed 2013, Flyby mission, 7th trajectory, Atlas V, launch 2021.
- Tour (p.6, Table 3 p.7): 59 targeted flybys. Max inclination 20.1 deg. Deterministic Delta-V 164 m/s
  after PRM. TID 2.8 Mrad (Si) behind 100 mil Al. Five phases:
  1. Approach: G0 500 km about 14 h before JOI; JOI at 11.1 R_J; PRM on a 202-day orbit; pump-down
     with 4 Ganymede and 2 Callisto flybys (above).
  2. Anti-Jupiter campaign (24 Europa flybys, 13.3 months). COT-1 and COT-2 are six 4:1 resonant
     transfers each at V-infinity about 4 km/s, cranked in opposite senses. A non-resonant E7-E8 transfer
     lies between them. E14-E16 avoid solar conjunction. Petal rotation E17-E24 alternates 4:1+ and 5:1-
     non-resonant transfers (V-infinity about 4.1) for gravity science. Ref. 28 is the petal-rotation paper.
  3. Illuminated-hemisphere transition (5.6 months): E25-E27 crank up and pump down, then the
     E28-C3...C9-E29 switch-flip, then E29-E30 pump up.
  4. Sub-Jupiter campaign: COT-3 (five 4:1, E31-E35), non-resonant E36-E37, then COT-4 with
     alternating 4:1 and 7:2 resonances to pull closest approach off the 180-deg meridian.
  5. Disposal by Ganymede impact.
- The COT (crank-over-the-top) rules (pp.8-9): more flybys per COT with higher V-infinity or a shorter period.
  An outbound start covers the anti-Jupiter hemisphere. The crank sign sets the node and the north-south
  build-up direction.
- Section V (p.13): comparison with the 2011 tour 11F5-A21. 13F7-A21 is longer for four reasons: a
  Ganymede+Callisto pump-down, more Europa flybys, a Europa-Callisto (not Europa-Ganymede) switch-flip
  to hold TID down, and solar-conjunction handling.
- Section VI (p.14): the baseline moved to a 2022 launch on SLS. A 2022 EVEEGA backup (7.58 yr TOF) is
  noted.

## 2. Checks (scripts kept)

`buffington_vca_check.py`: the Table 3 "Velocity" column at closest approach should equal
sqrt(v_inf^2 + 2 GM/(R + alt)). I used GM and mean radii from JPL SSD. It does: all 19 rows tested
(7 pump-down, 7 switch-flip Callisto, 5 Europa) agree to 0.008 km/s or better.

This gives a finding for the held Lam, Arrieta-Camacho & Buffington 2015 (AAS 15-657) and its digest. Lam's
Table A1 (13F7-A21) column headed "V-inf" lists 8.36, 6.93, 6.98, 6.84, 5.92, 5.81, 4.67 for G0-C2. These are
exactly Buffington's closest-approach speeds, not his V-infinity values (7.98, 6.39, 6.44, 6.43, 5.54,
5.24, 4.35). I checked them from the Lam PDF text layer. So the Lam 2015 digest line "v_inf falls ...
G 6.84 to 5.81 km/s, C 5.92 to 4.67 km/s" quotes flyby speeds. The V-infinity values are G 6.43 to 5.24
and C 5.54 to 4.35 (G3 to G4, C1 to C2). The conclusion (pumping, not a cycler) is unchanged. A second
small difference: Lam lists 6C2 as "O" (outbound); Buffington Table 3 lists Callisto2 as "I". I did not
resolve which is right.

Period-ratio check (same script): m x T_moon / n against the printed post-flyby period for 11 resonant
 rows (8:1, 4:1, 7:2 at Europa, 3:4, 1:1, 2:3 at Callisto). Ten agree within 1.1%. **G2 does not**:
 printed 5:1 with period 28.6 d, but 5 Ganymede periods are 35.8 d and 4 are 28.6 d. I re-read the cell
 at 400 dpi and it is printed "5". Lam 2015 Table A1 lists G2 as 4:1 (its m:n column for 2G2 reads 4, 1).
 So the printed 5 is a likely misprint. It changes no #943 verdict.

Cross-check with the held 21F31 digest (Campagnola et al. 2024): its tour table row "13F7 ... 3.5 ...
4:1 ... 45 / 5 / 9 ... 164 ... 20.1" agrees with this paper (3.5 yr, 45/5/9 flybys, 164 m/s, 20.1 deg).

## 3. Relation to the held Clipper papers

- Buffington, Campagnola & Petropoulos 2012 (AIAA 2012-5069, HELD, digest 2026-10-05): the 11F5-A21
  predecessor. It used a Europa-Ganymede switch-flip and had no Callisto flybys. This paper says why
  13F7-A21 changed to Callisto (TID).
- Lam, Arrieta-Camacho & Buffington 2015 (AAS 15-657, HELD): 13F7-A21 plus three later tours, with the
  full Table A1. See section 2 for the V-infinity labelling.
- Anderson, Campagnola & Buffington 2018 (JGCD, HELD): the journal form of ref. 28 (petal rotation).
- Campagnola et al. 2019 (JGCD 42(12), HELD): the tour-design techniques, including the GCGC cycler as an
  option. This 2014 paper only names cyclers.
- Lam, Buffington & Campagnola 2018 (AIAA 2018-0202, HELD) and Campagnola et al. 2024 ISSFD 21F31 (HELD):
  later tours. They have the same switch-flip structure, also with no repeating G-C or G-E segment.
- Cangahuala et al. 2025 (SSR, HELD): the flight mission design.

## 4. Citation mining

- Ref. 29 Russell & Strange 2009, "Planetary Moon Cycler Trajectories", JGCD 32(1):143-157: HELD
  (`russell-strange-2009-cycler-trajectories-planetary-moon-systems-JGCD-32-doi-10.2514-1.36610.pdf`).
- Ref. 22 Buffington, Campagnola & Petropoulos 2012 AIAA 2012-5069: HELD.
- Ref. 23 Buffington, Strange & Campagnola, "Global Moon Coverage via Hyperbolic Flybys", 23rd ISSFD 2012:
  not held (`ls | grep -i buffington` shows no ISSFD file). Not on the wanted list. Low priority (COT
  coverage geometry). Not a new candidate for the cycler work.
- Ref. 28 Anderson, Campagnola & Buffington, "Analysis of Petal Rotation Trajectory Characteristics",
  AIAA/AAS 2014: the conference form of the HELD JGCD 2018 paper. Not needed.
- Ref. 27 Petropoulos, Longuski & Bonfiglio 2000, JSR 37(6):776 (VEEGA to Jupiter): not held (`ls | grep -i
  petropoulos` shows only the 2012, 2014 and 2019 co-authored Clipper files). Already on the wanted list
  (row 67, gravity-assist path-design references). Interplanetary background.
- Ref. 19 Kloster, Petropoulos & Longuski 2010, Acta Astronautica (Europa orbiter tour with Io
  assists): not held, not in the index, not on the wanted list. New candidate only if Io-assisted tour
  background matters (for example `#953`); low priority.
- Ref. 14 Johannesen & D'Amario 1999 (Europa Orbiter tour), Ref. 17 Kahn, Campagnola & Croon 2004, and
  Ref. 18 Boutonnet & Schoenmaekers 2012 (JUICE, AAS 12-207): tour background. Not held (the held
  Johannesen and Boutonnet files are other papers: Lam-Johannesen-Kowalkowski 2008 Juno, Boutonnet et al.
  2024 JUICE SSR). Not on the wanted list. Not needed.
- Refs. 30-32 (navigation and maneuver analysis 2014): not needed.

*Check scripts and outputs named above are filed beside the PDF as `cyclers_pdf/papers/<pdf stem>-<script name>`.*

*Wanted-list row numbers in this digest are the batch-29 numbering; the list was renumbered in batch 30.*
