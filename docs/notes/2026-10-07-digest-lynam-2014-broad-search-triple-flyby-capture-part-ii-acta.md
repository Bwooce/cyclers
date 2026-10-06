# Digest: Lynam 2014, "Broad-search algorithms ... Callisto-Ganymede-Io triple flyby sequences from 2024 to 2040, Part II" (Acta Astronautica) (#960, #943)

A. E. Lynam (West Virginia University), "Broad-search algorithms for the spacecraft trajectory design of
Callisto-Ganymede-Io triple flyby sequences from 2024 to 2040, Part II: Lambert pathfinding and trajectory solutions",
Acta Astronautica 94(1):253-261, January 2014. doi 10.1016/j.actaastro.2013.07.020 (printed on p.253; Crossref title,
author, volume, issue and pages agree). Received 21 March 2013, accepted 5 July 2013, online 13 July 2013. Open access (CC BY).
- File given: `799e8c75-Broad-search_algorithms_for_the_spacecraft_traject.pdf`, 9 pages, text layer (publisher PDF).
  md5 672705501ea1659ac09908463aadf07b.
- **Proposed corpus filename:**
  `cyclers_pdf/papers/lynam-2014-broad-search-callisto-ganymede-io-triple-flyby-part2-lambert-pathfinding-acta-astronautica-94-253-doi-10.1016-j.actaastro.2013.07.020.pdf`
- **Wanted list:** not listed. Part I (below) is not listed either.
- **How I read it.**
  - I read the whole text layer. I read pp.258 (Figs. 4-5, Table 1), 259 (Fig. 6, capture timeline) and 260 (Fig. 7,
    discussion) on 130-dpi renders. Every number quoted below from those pages was read on the image.
  - The paper has no table that allows an arithmetic check except the dates. No disagreement found.
  - Page numbers are the printed journal pages.
  - The #943 collision check is in `collision.md` (written first). The verdict is **no collision**.

## 0. Verdict

**A capture-trajectory paper. It is not about cyclers.**
- It solves one inbound pass: Callisto flyby, Lambert arc to Ganymede, Ganymede flyby, Lambert arc to Io, Io flyby
  (before or after perijove). The aim is to cut the Jupiter orbit insertion (JOI) cost.
- **#943 collision: no collision** with gc-1, gc-2, ge-1, ge-2 or ge-3. There is no repeating sequence. The arrival
  V_inf at Jupiter is 3.6-5.7 km/s (Table 1). On such hyperbolas (rp 4.2-6.4 R_J) my patched-conic estimate of the moon v_inf is about
  12.1-14.6 km/s at Ganymede and 11.0-12.8 km/s at Callisto (INFERRED; `checks_lynam.out` sec. 5). The candidates have 1.4-3.9 km/s.
- **What it gives the project:** the 2024-2040 C-G-I window dates (sec. 2), for the X1 lane, and a published
  V_inf-matching + Lambert method for a three-moon pass (sec. 1). It is a neighbour of the X1 multi-moon tour work,
  not a catalogue source.
- **Catalogue implication (PROPOSAL only):** none.
- **Lynam 2012 PhD (wanted row 14; the brief says row 15, but the current list at commit 52d51c6f has the thesis at row 14):**
  - This paper is post-thesis work. The thesis is dated May 2012 (preview, acceptance form), and this paper was received
    in March 2013 from West Virginia University.
  - It carries out thesis future-work item 6.1, "Broad-search Algorithm for Ballistic Interplanetary Trajectories for
    Multiple-..." (preview TOC p.vi). It extends thesis secs. 2.7 (triple-satellite-aided capture, incl. 2.7.4-2.7.5
    Callistan triples) and 3.5 (targeting triple sequences).
  - It does not cover thesis Ch. 4 (navigation) or Ch. 5 (Laplace-resonant triple cyclers).
  - The thesis sec. 1.2 maps its chapters to refs [10]-[16] (preview p.2-3). None of the three Lynam broad-search papers
    is among them (they did not yet exist).
  - So this paper does not reduce the need for the full thesis. The thesis is wanted for Ch. 5 (cyclers), and the held
    Lynam & Longuski 2011 Acta already covers that.

## 1. Method (READ, pp.254-257)

- **Inputs:** the candidate Ganymede-flyby times from Part I's heuristic pruning; Galilean ephemeris positions and
  velocities at those times; and initial guesses for the semilatus recta p(Ca,Ga), p(Ga,Io) and the transfer times
  T(Ca,Ga), T(Ga,Io), taken from Part I's interpolation structures (p.254).
- **Moon propagation:** Callisto is propagated backward and Io forward from the Ganymede-flyby epoch by mean anomaly.
  Kepler's equation is solved by a fifth-order equation-of-the-centre series (eq. 3; e_Ca = 0.0074, e_Io = 0.0041).
- **Lambert:** third-order (Chebyshev-method) p-iteration after Herrick & Liu 1959 and Bate-Mueller-White, with
  central finite-difference derivatives (p_pert = 10,000 km, tolerance 1e-6) (eqs. 3-20).
- **V_inf matching (outer loop):** Newton-Raphson on (T_Ca,Ga, T_Ga,Io) to make the Ganymede in/out V_inf magnitudes
  equal (F1) and the turn angle equal to its maximum at 300 km altitude (F2, the inequality set to an equality) (eqs. 21-30).
- **Callisto and Io flybys:** left free. Each is discretised into 32 points on the V_inf-globe "small circle" (Strange,
  Russell & Buffington 2007 terminology; Fig. 3). 64 solutions are recorded per triple flyby.
- **Interplanetary filter:** keep solutions whose incoming Jupiter-centred V_inf is near anti-parallel to Jupiter's
  heliocentric velocity (eq. 31, angle from 180 deg).
- The author notes that the p-iteration solver has a 180-deg singularity and was "probably not the best option" (p.259).

## 2. Results (READ; image-checked)

- **Windows (Figs. 4-5, p.258; text p.257):** CGIP (C, G, Io inbound, then perijove) only in 2026 and 2033.
  CGPI (C, G, perijove, Io) in 2026 (two), 2029, 2030, 2033, 2034, 2036, 2037 and 2039: nine windows.
  CGPI is about three times as common as CGIP (p.257).
- **STK integration:** 11 windows were targeted. Five failed for lack of an Earth or Mars flyby. Six are feasible.
- **Table 1 (p.258), image-read:**

| window (Ganymede flyby) | arrival V_inf, km/s | capture orbit | flybys |
|---|---|---|---|
| CGIP 2 Apr 2026 | 4.7 | 590 d | Earth |
| CGIP 14 Jun 2033 | 5.1 | 509 d | Mars |
| CGPI 2 Apr 2026 | 5.3 | 587 d | Earth |
| CGPI 2 Dec 2029 | 3.6 | 97 d | Mars |
| CGPI 15 Jan 2033 | 5.7 | 2933 d | Mars |
| CGPI 31 Dec 2034 | 5.7 | hyperbolic | Earth and Mars |

- JOI estimates (p.258): about 100 m/s to reach 200 d and about 200 m/s to reach 100 d for the first three; 225 m/s (200 d)
  and 350 m/s (100 d) for the last two.
- **Best case, CGPI December 2029 (p.259, Figs. 6-7):** Earth launch September 2024 (C3 13.8 km^2/s^2), Mars flyby May 2026,
  Jupiter December 2029; 5.2 yr flight; unoptimised low-thrust Delta-V 7.3 km/s; peak acceleration about 2e-7 km/s^2
  (1 N for 5000 kg). At Jupiter: Callisto (95 km), Ganymede 19 h later (125 km), perijove 4.2 R_J 17 h later, Io 4 h later
  (137 km), then a 97-day orbit with no JOI burn.
- **Comparison (p.260):** against the Strange et al. 2012 Callisto-perijove-Ganymede (CPG) double capture of June 2027:
  97 d vs 354 d capture period, perijove 4.2 vs 9.4 R_J. The author says the other five windows are probably worse than
  CPG double capture.

## 3. Checks

- Fig. 7's STK label reads "2 Dec 2029 06:01:03.163", which agrees with the Table 1 date for that window (image).
- No other arithmetic is possible from the printed data. No disagreements.

## 4. Citation mining (35 references, pp.260-261)

Checked with `ls cyclers_pdf/papers | grep -i` and `grep -i docs/notes/CORPUS_INDEX.md`.

| ref | work | status |
|---|---|---|
| [30] | **Lynam, A. E. (2014), "Broad-search algorithms ... Part I: Heuristic pruning of the search space", Acta Astronautica 94(1):246-252, doi 10.1016/j.actaastro.2013.07.018** (Crossref confirms title, author, pages) | not held; not on the wanted list. **New candidate (low priority).** It defines the phase-angle pruning and the interpolation structures this paper uses. L15tq says its own method is "essentially more general and faster variants" of Part I, so Part I adds history, not data. |
| [20] | Lynam, Kloster & Longuski (2011), CMDA 109:59-84 (Laplace-resonance multiple capture) | not held as a filed paper; a manuscript copy (`5d107f14-MSAC_Laplace_Resonance.pdf`) is in batch 35 with a sibling agent |
| [21] | Lynam & Longuski (2011), JGCD 34(5):1485-1494 | in batch 35 with a sibling agent (`05b40581-lynam2011_1.pdf`, doi 10.2514/1.53251) |
| [22] | Lynam & Longuski (2012), "Preliminary analysis for the navigation of multiple-satellite-aided capture sequences at Jupiter", Acta Astronautica 70:33-43 | not held; not on the wanted list. Navigation only (thesis Ch. 4). Not needed for #943. |
| [19] | Lynam, Kloster & Longuski (2009), AAS 09-424 | not held; conference form of [20]. |
| [34] | Strange, Russell & Buffington (2007), "Mapping the V-infinity globe", AAS 07-277 | HELD (`strange-russell-buffington-2007-mapping-v-infinity-globe-AAS-07-277.pdf`) |
| [18] | Strange, Landau, Hofer et al. (2012), "Solar electric propulsion gravity-assist tours for Jupiter missions", AIAA 2012-4518 | not held; not on the wanted list. The CPG comparison case. Low priority. |
| [17] | Landau, Strange & Lam (2010), "Solar electric propulsion with satellite flyby for Jovian capture", AAS/AIAA SFM | not held. Low priority. |
| [9]-[12] | Longman 1968 (RAND); Longman & Schneider 1970 (JSR 7:570); Cline 1979 (Celest. Mech. 19:405); Nock & Uphoff 1979 (AAS 79-165) | none held (no file hit for longman, cline or uphoff; the "nock" files are Mars papers). Satellite-aided capture history; no cycler content expected. |
| [13] | Johannesen & D'Amario 1999, AAS 99-330 | not held (the one "johannesen" file is Lam et al. 2008, Juno). |
| [35] | Kloster, Petropoulos & Longuski (2010), "Europa orbiter tour design with Io gravity assists", Acta Astronautica 68:931-946 | not held (no kloster hit). Tour design; low priority. |
| others | Galileo and Cassini operations papers [1]-[8]; Lambert solvers [23]-[27] (Gooding 1990 and Battin 1999 not held); Williams 1990 and Longuski & Williams 1991 [28,29]; Vallado 2007 [31] (only Vallado 1991 USAFA TR held); Bate et al. [32]; Amat et al. 2003 [33] | none needed for the project |

- New candidates: Part I (low priority). Nothing for #943.

*Wanted-list row numbers in this digest are those of the list at commit 52d51c6f.*

*Filed as `cyclers_pdf/papers/lynam-2014-broad-search-algorithms-callisto-ganymede-io-triple-flyby-capture-part-ii-acta-astronautica-94-253-doi-10.1016-j.actaastro.2013.07.020.pdf`. Check scripts, outputs and other files named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`.*

*Wanted-list row numbers in this digest are the batch-34 numbering; the list was renumbered in batch 35.*
