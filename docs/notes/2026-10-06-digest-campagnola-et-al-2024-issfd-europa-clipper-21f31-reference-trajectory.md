# Digest: Campagnola et al. 2024, "Europa Clipper Mission Analysis: Design of the 21F31 Reference Trajectory" (ISSFD 2024) (#960, #943)

S. Campagnola, B. B. Buffington, B. Anderson, E. Pellegrini, R. Restrepo, T. Lam, B. Bradley (JPL),
C. Scott, D. Ellison, K. Bokelman, M. Ozimek & G. Pradipto (APL), "Europa Clipper Mission Analysis: Design
of the 21F31 Reference Trajectory", 29th International Symposium on Space Flight Dynamics (ISSFD), ESOC
Darmstadt, 22-26 April 2024, paper 19-3. No DOI.
- **Conference counterpart of the JAS 2025 paper** (Campagnola, Buffington, Anderson, Pellegrini,
  Restrepo & Lam, "Europa Clipper Mission Design: Design of the 21F31 Reference Tour", JAS 72, doi
  10.1007/s40295-025-00527-1). The ISSFD abstract is essentially verbatim the JAS abstract. The JAS
  version is not held and stays in the wanted list, marked "conference version held".
- Filed as `cyclers_pdf/papers/campagnola-et-al-2024-europa-clipper-mission-analysis-design-21f31-reference-trajectory-issfd-2024-19-3.pdf`.
  23 pages, text layer, md5 032c1c312ed59e416550f798323455cb. Supplied by the owner.
- I searched the full text for cycler, Ganymede-Callisto and Ganymede-Europa content, and read sec. II
  (history; Fig. 7 on the page image), sec. III-A, and sec. IV's tour phases. The Table 1 flyby list is
  printed rotated; the same list was read from the held Cangahuala 2025 SSR Table 2 instead.

## 0. Verdict (sent to the lead and twobody-gen2-opus before this digest)

**No repeating Ganymede-Callisto or Ganymede-Europa segment. No collision with `#943` gc-1, gc-2, ge-1,
ge-2 or ge-3.**
- "Cycler" appears once in the paper, as an untried option for the transition to Europa Campaign 2:
  "petal rotation using Europa flybys ..., Ganymede flybys, or Callisto flybys; different types of
  cyclers; a Europa pi-transfer; or a 'switch-flip.'"
- The flown 21F31_V6 uses a Callisto petal rotation:
  - E26 targets C01 at v_inf 4.8 km/s;
  - G06 "is then utilized to leverage the V_inf at Callisto down to 3.6 km/s";
  - from C03 the non-resonant transfers 2:2-, 1:1+ and 3:3- (the last of 43 d, with a solar
    occultation) rotate the line of apsides;
  - G07 "leverage[s] up the V_inf at Callisto to 5.1 km/s".
  - So Ganymede is used twice, as a one-shot v_inf lever.
- G-E occurs only in the one-shot pump-down (G01-G03, E01-E02, G04-G05) and the G08 impact. Pump-down
  v_inf from the SSR Table 2 page image: Ganymede 8.73-7.65 km/s, Europa 6.30/6.38.

## 1. Content (READ)

- **Tour design process (sec. II):**
  - seven design cycles over a decade with the Project Science Group;
  - 21F31_V4 adopted in 2021 (cycle 6);
  - 21F31_V6 is the final tweak (cycle 7).
- **Fig. 7 (p.6, page image): history of baseline tours.**

| Tour | Gate | Launch | Arrival | Tour yr | EC1 Europa res. | E / G / C flybys | Det. Delta-V post-PRM, m/s | i_max, deg | TID, Mrad |
|---|---|---|---|---|---|---|---|---|---|
| 13F7 | MCR | 2021 | - | 3.5 | 4:1 | 45 / 5 / 9 | 164 | 20.1 | 2.82 |
| 15F10 | SRR/MDR | 2022 | - | 3.4 | 4:1 | 42 / 4 / 8 | 118 | 21.2 | 2.99 |
| 17F12 V2 | PDR | 6/4/22 | 12/23/24 | 3.7 | 4:1 | 46 / 4 / 9 | 182 | 18.9 | 2.5 |
| 19F22/F23 | CDR | 11/7/23 | 9/29/29 | 3.84 | 4:1 | 51 / 6 / 7 | 199 | 21 | 2.88 |
| 21F31 V4 | SIR | 10/10/24 | 4/10/30 | 4.27 | 6:1 | 53 / 7 / 9 | 225.2 | 7.5 | 2.97 |
| 21F31 V6 | FRR | 10/10/24 | 4/10/30 | 4.27 | 6:1 | 53 / 7 / 9 | 214.8 | 7.7 | 2.97 |

  - The figure gives no transition type per tour, so whether an earlier tour used the GCGC cycler is
    not stated here. Campagnola 2019 (held) says the GCGC idea was "inspired by" tour 15F09.
- **Sec. III-A:** resonant (n:m) and non-resonant transfers; a pi-transfer is "a special case of a
  non-res where the time-of-flight is an integer multiple of the gravity assist body period plus 1/2".
- **Sec. IV:** 21F31_V6 has 53 E, 7 G and 9 C flybys over 4.27 yr; TID 2.97 Mrad. The pump-down was
  deliberately under-used to absorb a +/-20 m/s JOI error and a 2-h outage (alternate pump-downs
  reconnect at E04). The EC1/EC2 crank-over-the-top sequences are as in the SSR paper.

## 2. Gate relevance

- `#943` gc / ge: no collision (sec. 0). The published two-working-body G-C cycler record remains
  Campagnola 2019 GCGC alone.
- Nearest values: the petal-rotation Callisto v_inf of 3.6 is 0.56 km/s above gc-2's 3.039; the Ganymede
  values (5.35/4.30, SSR) are far above both gc candidates.

## 3. Positive controls

- None for cyclers. Fig. 7 tour statistics (above).

## 4. Citation mining

The references overlap the held Clipper set (Campagnola 2019, Buffington 2012, Lam 2015,
Anderson 2018). Not held: [3] Buffington et al., "Evolution of Trajectory Design Requirements ..." (Clipper
series); [22] Arya, French, Pellegrini & Campagnola, "Stochastic Advance DeltaV99 Evolutionary ...";
[12] Buffington & Strange 2007 (patched-integrated design); [15] Buffington, Strange & Campagnola
("Global Moon Coverage via Hyperbolic Flybys", the COT reference, already noted in the Lam 2015 digest).
Low priority for X1; not added to the wanted list.
