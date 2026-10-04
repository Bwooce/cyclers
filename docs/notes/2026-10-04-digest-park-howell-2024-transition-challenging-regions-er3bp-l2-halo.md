# Digest: Park & Howell (2024), "Characterizing Transition-Challenging Regions Leveraging the Elliptic Restricted Three-Body Problem: L2 Halo Orbits"

Beom Park (Ph.D. candidate) and Kathleen C. Howell, School of Aeronautics and Astronautics, Purdue University. AIAA conference
paper, 22 pages. The venue and paper number are NOT printed anywhere in the document (no header, footer or copyright block); the
corpus filename says "AIAA-scitech", and the reference list cites a SciTech 2024 paper (Sanaga & Howell, ref. 15), so a 2024
AIAA SciTech Forum paper is INFERRED, not READ.
Filed in the private paper corpus as
park-howell-2024-characterizing-transition-challenging-regions-ER3BP-L2-halo-orbits-AIAA-scitech.pdf

Digested 2026-10-04 from the page images (all 22 pages read). Statements are marked READ (seen on the page, with page, section,
equation, figure or table), COMPUTED (our arithmetic) or INFERRED (our reading). Page numbers are the paper's printed page numbers,
which equal the PDF page numbers. Context: tasks #893 (every core model needs an identity test and a published positive control)
and #884; `core/er3bp.py`; the later paper by the same group, digested in
`2026-10-04-digest-singh-park-howell-2026-aas-26-654-l2-families-intermediary-models.md`.

## 0. What the paper is

READ (abstract, p1): "This investigation leverages the Elliptic Restricted Three-Body Problem (ER3BP) as an intermediate model
between the CR3BP and HFEM, employing numerical continuation and bifurcation analysis to characterize the challenges associated with
the transition process. Specifically, the focus is on the Earth-Moon L2 halo orbit family. A subset of the family exhibits a
proliferation of fold bifurcations with respect to the eccentricity of the model, serving as an indicator for transition-challenging
behavior as analyzed within the ER3BP."

It is a bifurcation-and-continuation study with a single printed table (Table 1, a list of p:q ranges, section 3 below). It prints
NO initial condition, no Jacobi constant, no stability index and no table of orbit states. All results are in figures (hodographs of
x0 against e or T, and a global scatter of fold eccentricities against period). It is a CR3BP-to-ER3BP-to-ephemeris-model
(HFEM, "Higher-Fidelity Ephemeris Model") transition study for the Earth-Moon L2 southern halo family; the 9:2 NRHO is named as
the motivating orbit (p1, p4).

## 1. The model exactly as printed

### 1.1 Frame (READ, p3, section II.A, eqs. 1 to 3)

Pulsating-rotating frame. Three unit vectors: x-hat from the Earth to the Moon, z-hat along the angular momentum vector of the
celestial bodies, y-hat completing the dextral triad. Origin at the Earth-Moon barycentre. Quote: "Regardless of the instantaneous
dimensional Earth-Moon distance (characteristic length), the frame adopts a consistent nondimensional unit, where the
nondimensional Earth-Moon distance is always unity." Positions: `r_E = -mu x-hat` (eq. 1), `r_M = (1 - mu) x-hat` (eq. 2),
`r_s = x x-hat + y y-hat + z z-hat` (eq. 3), with `mu = mu_M/(mu_E + mu_M)`. "Note that the Earth and Moon position vectors are
fixed within the rotating frame on the x-hat axis, separated by a unit nondimensional distance." The primaries are point masses,
with gravitational parameters "consistent with the JPL ephemerides, DE440.bsp" (ref. 17). The numerical value of mu is NOT printed
for the main computation (the text gives "approximately 0.0125 (the Earth-Moon system value)" only as one of three sample values in
section V.C.3, p17; the actual Earth-Moon value is 0.01215, so "0.0125" is a rounded display value). The frame is the project's frame
(Earth at -mu, Moon at 1 - mu).

### 1.2 CR3BP (READ, p3, section II.B, eq. 4)

`r'' = -2 z-hat x r' + grad(Omega_C)` with `Omega_C = (x^2 + y^2)/2 + Omega`, `Omega = (1 - mu)/r_Es + mu/r_Ms`; independent variable the
nondimensional time t, dot = d/dt.

### 1.3 ER3BP (READ, p3, section II.C, eq. 5)

Independent variable the true anomaly f of the Earth-Moon conic motion; eccentricity e with 0 <= e < 1; "For the Earth-Moon system,
0.055 represents a sample realistic value" (ref. 10). Equation of motion, with prime = d/df:

`d^2 r/df^2 = -2 z-hat x dr/df + grad(Omega_E)`   (eq. 5)

with the pseudo-potential printed as `Omega_E = Omega_C/(1 + e cos f) - (e cos f) z^2/(2 + 2 e cos f)`. "Note that the ER3BP reduces to
the CR3BP and f also reduces to t for e = 0." The paper notes that other formulations exist (Hiday-Johnston & Howell 1994, ref. 18)
"cast in a non-dimensional rotating frame where the Earth and Moon are not fixed".

COMPUTED: eq. 5's potential is algebraically `[ (x^2 + y^2)/2 + Omega - e z^2 cos(f)/2 ] / (1 + e cos f)`, which is the same function as
the one in the `core/er3bp.py` docstring, `omega = (1/(1 + e cos f)) [0.5 (x^2 + y^2 - e z^2 cos f) + (1 - mu)/r1 + mu/r2]`.
So the paper's model is the project's model (pulsating frame, true anomaly, primaries fixed at -mu and 1 - mu). The paper does not
print the gradient written out, nor an identity test.

### 1.4 Symmetries and sign of e (READ, p6, section III.C.1)

"From Eq. (5), note that negative eccentricity in the ER3BP dynamics result in a phase shift in by Delta f = 180 deg." Hence for an
odd p the two counterparts are linked by continuation in positive and negative e; for an even p they are not. Only e >= 0 is used.

## 2. What is computed

### 2.1 CR3BP L2 southern halo family (READ, p4, section III.A, figs. 2, 3)

Southern branch (mirror of the northern) of the L2 halo family; parameterised by period T (days) because the period is monotonic along
the family; time unit `t* ~ 375699 s` ("approximately 375699 seconds is multiplied with the non-dimensional period"). States are
parameterised by the longitudinal angle theta; theta = 0 is the apolune. The family shown is the part above the lunar radius.
Fig. 2b is the hodograph x0 (at apolune) against T from about 4 to 15 days (x0 about 1.0 to 1.17); no tabulated values.

### 2.2 Stability, rotation number and bifurcations in the CR3BP (READ, p4 to p5, eq. 6, figs. 3, 4)

The monodromy matrix M has three reciprocal pairs, with a trivial pair lambda_{1,2} = 1. Centre eigenvalues lambda_{3,4}, lambda_{5,6}
give the rotation number `rho = arccos(Re(lambda))`, rho in [0, pi] (eq. 6). When rho is a rational multiple of 2 pi
(p:q = 2 pi : rho) the halo undergoes a period-multiplying bifurcation into a higher-period family; otherwise Hopf-type bifurcation to
two-parameter QPO families (T, rho), which collapse to higher-period periodic orbits at resonance. Example (Fig. 4): a period-tripling
bifurcation at rho = 2 pi/3 (about T = 40 days for the three-stack halo, read off the plot), where a new three-lobed family appears
with constant rho. Fig. 3 plots rho of lambda_{3,4} and lambda_{5,6} over about T = 6 to 15 days (lambda_{3,4} drops to zero near
T = 10.3 days, read off the plot; no tabulated values).

### 2.3 How periodicity is defined in the elliptic problem (READ, p6, section III.C, eq. 7)

The ER3BP is periodically perturbed with period 2 pi in f, so a CR3BP periodic orbit survives as an ER3BP periodic orbit only if the
orbit period is commensurate with 2 pi. Quote: "The ER3BP POs considered in this study are characterized by the rational ratio p:q,
where p:q = 2 pi : T, with p and q being positive coprime integers." Equivalently the stroboscopic mapping time is `T = 2 pi q/p`
(n.d.) and the second rotation angle is `rho_E = T` (eq. 7); the full ER3BP orbit closes after p halo revolutions, which is q
revolutions of the primaries (full period 2 pi q in f). READ example (p8): "the stroboscopic mapping time is fixed as
T = 7 . 2 pi/13 ~ 3.383253 n.d." for the 13:7 orbit, "approximately 14.7 days". COMPUTED: 14 pi/13 = 3.383253 (agrees), and
3.383253 x 375699 s / 86400 s = 14.71 days (agrees); 2 pi/3 = 2.0944 n.d. = 9.107 days (agrees with the text's "approximately 9.1
days" for 3:1); T ~ 14.8 days = 3.403 n.d. (agrees with the text's rho_E ~ 3.4 rad).

Notation note (INFERRED): in the appendix (p19 to p20) "T = 2 pi q" and "the total period is 2 q pi" denote the period of the whole ER3BP
orbit, while in the main text T denotes the stroboscopic (single halo revolution) time 2 pi q/p. Both are used consistently within
their sections.

Two ER3BP counterparts exist for each p:q (READ, p6, Figs. 5, 6): counterpart A is started at the perpendicular crossing with
f0 = 0 deg and counterpart B at f0 = 180 deg for odd p (both from the apolune, theta = 0, of the CR3BP halo); for even p, counterpart B
starts at the perilune (theta = 180 deg) with f0 = 0 (and appears at a 90 deg offset in f0 on the theta = 0 stroboscopic map).
The trivial pair lambda_{1,2} = 1 of the CR3BP orbit bifurcates into a saddle and a centre pair for each counterpart. ER3BP periodic
orbits are isolated (island) solutions at fixed e, which is why each p:q needs its own continuation in e.

### 2.4 ER3BP quasi-periodic orbits (READ, p6 to p7)

For irrational T/(2 pi) the periodic orbits become QPOs (two-parameter family in T and e). Computed with the GMOS algorithm (Olikara &
Scheeres, ref. 28, building on Gomez & Mondelo, ref. 27), a two-point boundary value problem for the invariance condition, with 51
discretisation nodes on the invariant curve; continuation in e rather than in rho_E. Fig. 7 example at T ~ 14.8 days (rho_E ~ 3.4 rad).

### 2.5 Transition-challenging region and how it is characterised (READ, pp1, 8 to 11, 14 to 15)

Definition from the introduction (p1, p2): "transition-challenging" behaviour is where using the CR3BP orbit as the initial guess for a
differential corrector gives "significantly shorter horizon times" in the higher-fidelity ephemeris model (HFEM) "for a subset of
the Earth-Moon L2 halo family, corresponding to a range from 8.6 to 11.0 days for the periods of the CR3BP halo orbits" (Fig. 1,
recreated from the authors' earlier paper, ref. 11, which shows horizon length 0 to 20 years against period 6 to 15 days). That
subset "includes a 3:1 synodic resonant orbit" (p2).

Characterising quantity: the fold bifurcation in the ER3BP continuation in e. READ (p10): "a fold bifurcation refers to a turn in the
eccentricity, or, as demonstrated in Fig. 8(b), a local extremum in eccentricity along the hodograph." Two representative
bifurcation diagrams (Fig. 8, p9):

- Nominal case, Fig. 8(a): the ER3BP branch leaves the CR3BP family at a single bifurcation (the loss of the trivial eigenvalue pair)
  and e increases monotonically up to e = 0.055 at constant T = rho_E, with no turn in e. Example: 13:7 PO and a nearby QPO with
  T ~ 3.383250 n.d. (Fig. 9, p9), "approximately 14.7 days and is far from the predescribed transition-challenging region".
- Fold case, Fig. 8(b): the branch turns in e (fold, a unity eigenvalue pair), returns to e = 0, connects to a distinct CR3BP
  structure at a "resonance bifurcation" (rho = rho_E + 2 k pi), and that constant-rho CR3BP branch reconnects to the original halo
  at a period-multiplying (or Hopf) bifurcation. Quote (p10): "the continuous family within the ER3BP in yellow does not reach
  e = 0.055 without connecting back to e = 0. This trait implies that a successful smooth transition process may be required to
  leverage a nearby period-multiplied PO or a QPO associated with the center mode of the original CR3BP halo orbit, indicating that
  the CR3BP PO itself potentially does not serve as a suitable initial guess for the transition process."
- Example (Fig. 10, p11): 53:21 PO with T ~ 2.48956 n.d. and nearby QPO with T ~ 2.48963 n.d.; two folds per counterpart, with
  fold eccentricities of order 3 to 4 x 10^-3 (axis "x 10^-3", read off the plot), returning to e = 0.

Relation to broken bifurcations (READ, p10 to p11 and Appendix B, p19 to p20, Fig. 21): for a one-parameter model such as the ER3BP
at fixed mu, perfect transcritical or pitchfork intersections of families are non-generic and are broken into two disconnected
continuous branches; the fold is "a practical metric", not equivalent to a broken bifurcation: "Broken bifurcations do not
necessarily involve a fold in eccentricity, and the presence of a fold bifurcation alone does not confirm a broken branch."

Global finding (READ, p14, section V.A, Figs. 15, 16, Table 1): PO continuation in e up to a maximum of 0.055 was run for many
(p, q) in the ranges of Table 1 (section 3). Quote: "The intermediate region from T = 8.6 to 11 days are characterized by the
existence of fold bifurcation at e < 0.055 and the evolution back to zero eccentricity. Beyond this period range, the yellow markers
govern, i.e., this continuation scheme reaches e = 0.055 without any fold bifurcations." and "Therefore, the proliferation of ER3BP
PO fold bifurcations before reaching e = 0.055 successfully predicts challenges in HFEM transition." Marker classes: "PO reaches
e = 0 without fold" (reached e = 0.055 without a fold); "PO folds but reaches e = 0.055" (even number of folds); "PO folds and
reverts to e = 0". Exceptions quoted (p14 to p15): the 9:2 ratio at T ~ 6.07 days folds at e ~ 0.04, outside the challenging
region; near T ~ 9.1 days, which "corresponds to the 3:1 resonance", both counterparts reach e = 0.055 without a fold; between 9.5
and 10 days four ratios reach e = 0.055 with an even number of folds. QPO continuation boundaries agree broadly with the PO fold
eccentricities (grey bands in Figs. 15, 16), with slight deviations for T < 9 days that the authors attribute to multiple broken
branches and to step-size differences.

Higher eccentricity (READ, p15, Fig. 17): continuing past 0.055 from the no-fold ratios, "all ratios exhibit turns in eccentricity
before e = 0.99, or a value near unity"; the NRHO region (T < 8.6 days) has "smaller values of eccentricity at the fold
bifurcations as compared to the region defined by T > 11."

### 2.6 Numerical caveats the paper states (READ, p13 to p14, section IV.D, Figs. 13, 14)

Fold detection depends on the formulation: using the perpendicular-crossing (mirror) assumption misses non-symmetric branches; denser
steps expose more broken branches; Fig. 13 (73:24 counterpart A) shows three different branches found with three different step
sizes, and natural-parameter continuation can jump between disconnected branches. QPO continuation fails earlier than PO continuation
when the invariant curve becomes too complex (Fig. 14, 59:21).

### 2.7 Mitigation strategies (READ, pp16 to 18, section V.C, Figs. 18 to 20)

1. Continue from a nearby CR3BP higher-period PO reached after the fold (Fig. 18, 5:2 counterpart A, ending at e = 0.055 as a
   multi-lobed orbit; f0 is shifted from 0 to 180 deg at the resonance bifurcation).
2. Deliberately jump to the other broken branch by natural-parameter continuation in e at the sharp turn (Fig. 19, 19:6 counterpart A;
   quoted "good luck is needed", citing Seydel).
3. Continue in mu (Fig. 20, three mu values for 8.6 <= T <= 11 days: 0.1, "approximately 0.0125", 0.001); fold eccentricities shift
   upward for the larger mu and downward for the smaller, but the proliferation of folds near the original period range persists.
4. Continue in T (only for QPOs; POs exist at discrete p:q); mentioned for the Hill four-body model.

Conclusion (p18): "fold bifurcations within the ER3BP prove to be an effective metric for anticipating transition challenges within
the ER3BP alone, eliminating the need for extensive numerical transition tests from the CR3BP to the HFEM." The paper itself lists the
verification for other families and systems as future work.

## 3. Printed tables, transcribed digit by digit

The paper contains one table.

**Table 1 (p15), "Range for p:q values"** (the ranges of (p, q) used in the global continuation)

| Periods (days) | Range for p | Range for q |
|---|---|---|
| T < 8.6 | p <= 60 | q <= 30 |
| 8.6 <= T <= 11 | p <= 100 | q <= 50 |
| 11 < T | p <= 50 | q <= 25 |

(READ; the inequality signs are as printed.) No table prints initial conditions, periods, Jacobi constants, eccentricities,
stability indices or bifurcation values. Numbers printed in the running text are collected here for convenience (all READ):

| Quantity | Value | Where |
|---|---|---|
| Realistic Earth-Moon eccentricity used | 0.055 | p3, throughout |
| Time scale t* | approximately 375699 s | p4 |
| Transition-challenging period range | 8.6 to 11 days (11.0 on p1; Fig. 1 vertical lines near 8.6 and 11.1) | p1, p14 |
| 3:1 resonance | rho_E = 2 pi/3 rad, "approximately 9.1 days" | p6, p14 |
| 13:7 PO | T = 7 . 2 pi/13 ~ 3.383253 n.d. (about 14.7 days); nearby QPO T ~ 3.383250 n.d. | p8, p9 |
| 53:21 PO | T ~ 2.48956 n.d.; nearby QPO T ~ 2.48963 n.d. | p11 |
| QPO example | T ~ 14.8 days, rho_E ~ 3.4 rad, 51 nodes | p7 |
| 9:2 outlier | fold at e ~ 0.04 near T ~ 6.07 days | p14 |
| Continuation settings | pseudo-arclength step ds = 0.001; n = 2p + 1 segments; 51 QPO nodes | p20 (eq. 18 text) |
| mu values explored | 0.1 (about Pluto-Charon), about 0.0125 (Earth-Moon), 0.001 (about Sun-Jupiter) | p17 |
| Free-variable / constraint lengths | 6n - 2 and 6n - 3 | p20 to p21 |

## 4. Method, as far as needed to read the numbers

Multiple shooting with n = 2p + 1 equally spaced (in f) segments from f0 to f0 + q pi (half the orbit, using the perpendicular
crossing with the x-z plane for mirror symmetry). Free vector `X = [x0 z0 ydot0 s2 ... sn e]` (eq. 16, length 6n - 2), constraint
vector `F = [F_c F_p]` (continuity plus perpendicular crossing at f0 + q pi, eq. 17, length 6n - 3), so the Jacobian has a
one-dimensional null space `V_N`; pseudo-arclength continuation adds `dX . V_N - ds = 0` (eq. 18) with ds = 0.001. The augmented
7x7 monodromy matrix (eq. 10, p19) shows that at a fold (turn in e) the 6x6 monodromy matrix has a unity eigenvalue pair with the
family tangent as the eigenvector (eqs. 13 to 15). Larger n and smaller ds raise the chance of finding extra folds.

## 5. Positive controls for the project

The paper prints no state, no initial condition and no period beyond the few values in the table of section 3, and the mu of the
computation is not printed. There is therefore NO printed periodic orbit that `core/er3bp.py` can be started from and checked for
closure. Honest candidates, in order of usefulness, all weak:

1. READ resonance arithmetic (checks the period convention, not a state): the 13:7 orbit is `T = 14 pi/13 = 3.383253` n.d. per halo
   revolution and the ER3BP orbit closes after 14 pi in f. A closure test for any CR3BP halo whose period is `2 pi q/p` seeded at
   e = 0 and continued in e at fixed f0 (0 or pi) would test the same structure.
2. READ qualitative counts that a continuation can be compared with, once seeds are made by the project's own CR3BP halo corrector
   (mu 0.01215, same frame): the 3:1 halo (T = 2 pi/3 = 2.0944 n.d., 9.1 days) has both counterparts reach e = 0.055 without a fold
   (p14), and the 13:7 (T = 3.383253) likewise; the 53:21 (T = 2.48956) has a fold at e of order 3 to 4 x 10^-3 and returns to
   e = 0 (p11). These are figure-level facts (no numbers beyond the periods), so they can falsify gross errors but not small ones.
3. Axis read-offs (INFERRED from the plots, accurate to the last printed tick, about 0.0005; not printed values): the e = 0 end of
   the x0 hodograph is about 1.0640 (3:1 counterpart A, Fig. 5a), about 1.1435 (2:1, Fig. 6a), about 1.1795 (QPO at T ~ 14.8 days,
   Fig. 7a) and about 1.176 to 1.178 (13:7, Fig. 9). These could be compared with the CR3BP halo apolune x0 at the matching period
   as a plausibility check on a seed; they must not be used as goldens.

The closest printed-number positive controls for the elliptic Earth-Moon problem remain Neelakantan & Ramanan (2022) Table 8
(already in `tests/core/test_er3bp_neelakantan_2022.py`) and, for Sun-Mercury, Peng, Bai & Xu (2017) Table 2 (see
`2026-06-25-digest-peng-2017-sun-mercury-ERTBP.md`).

## 6. What this means for the project

Scope. The paper is about periodic and quasi-periodic halo orbits of the Earth-Moon elliptic problem and the ease of transitioning
them to an ephemeris model. It contains no cycler, no quasi-cycler and no resonant orbit with flybys; its resonances are
resonances between the halo revolution period and the Moon's orbital period (p:q = 2 pi : T), with no close approach to the Moon.
It falls outside the catalogue's four classes (cycler, quasi-cycler, precursor MGA, MGA tour); it is background for the elliptic
model and for transition methodology.

Overlap with open tasks.
- #893 (identity test and positive control for every core model): the paper confirms the model equation (section 1.3) but offers no
  printed orbit; it does not close the ER3BP positive-control requirement beyond what the Neelakantan table already provides.
- #884 and #891 (the sun-forced and bicircular work): the paper's central lesson is that the CR3BP-to-higher-fidelity transition
  can be anticipated from an intermediate periodically perturbed model, and that families in such models are generically broken
  rather than intersecting. That is the same mechanism the project met in its own sun-forced periodic-orbit searches; the paper
  supports treating an isolated failure to continue as a fold, with a check by jumping branch or by continuing in a second
  parameter (mu or T), before declaring a negative.

Reusable methods (INFERRED): the augmented 7x7 monodromy test for a fold (a unity eigenvalue pair with the family tangent as
eigenvector, eqs. 10 to 15); pseudo-arclength continuation in e with n = 2p + 1 segments (eqs. 16 to 18); using the nearby CR3BP
higher-period PO after a fold (mitigation 1); the resonance classification by p:q = 2 pi : T with the two counterparts A and B.

Notes on slips, stated factually (INFERRED). (a) The text on p16 says "Peng and Hao [34]" while reference 34 lists Peng, Bai and Xu
(2017). (b) References 5 and 7 are the same paper (Dei Tos & Topputo 2017). (c) The mu value 0.0125 is a rounded form of the
Earth-Moon value (0.01215). (d) The text gives the challenging range as "8.6 to 11.0 days" (p1) and "8.6 and 11 days" elsewhere,
while Fig. 1's right-hand vertical line is drawn near 11.1 days; the discrepancy is within the plotting resolution.

## 7. Peng & Xu (2015) and the corrupt held copy

The paper does NOT reproduce any table, initial condition or result from Peng & Xu (2015, Celestial Mechanics and Dynamical
Astronomy 123:279). It cites it once, at p6 (section III.C.1), as the source of "more information on the computation of two ER3BP PO
counterparts", and the reference appears as ref. 25. The Sun-Mercury 5:2 behaviour is attributed to Peng, Bai & Xu (2017, ref. 34)
at p16. So the paper fills none of the Peng & Xu (2015) gap: no number can be recovered from it.

## 8. References cited by the paper (full citations as printed; numbers are the paper's)

Multi-revolution, resonant and elliptic-problem halo orbits:
- [1] Zimovan-Spreen, E. M., Howell, K. C., and Davis, D. C., "Near rectilinear halo orbits and nearby higher-period dynamical
  structures: orbital stability and resonance properties," Celestial Mechanics and Dynamical Astronomy, Vol. 132, No. 5, 2020, pp. 1-25.
- [14] Sanaga, R. R., and Howell, K. C., "Synodic resonant halo orbits in the Hill restricted four-body problem," AAS/AIAA
  Astrodynamics Specialist Conference, Austin, Texas, 2023.
- [15] Sanaga, R. R., and Howell, K. C., "Analyzing the challenging region in the Earth-Moon L2 halo family through Hill restricted
  four-body problem dynamics," AIAA SCITECH 2024 Forum, 2024.
- [16] Henry, D. B., Rosales, J. J., Brown, G. M., and Scheeres, D. J., "Quasi-periodic orbits around Earth-Moon L1 and L2 in the Hill
  restricted four-body problem," AAS/AIAA Astrodynamics Specialist Conference, Big Sky, Montana, 2023.
- [18] Hiday-Johnston, L., and Howell, K., "Transfers between libration-point orbits in the elliptic restricted problem," Celestial
  Mechanics and Dynamical Astronomy, Vol. 58, 1994, pp. 317-337.
- [23] Ferrari, F., and Lavagna, M., "Periodic motion around libration points in the elliptic restricted three-body problem,"
  Nonlinear Dynamics, Vol. 93, No. 2, 2018, pp. 453-462.
- [24] Campagnola, S., Lo, M., and Newton, P., "Subregions of motion and elliptic halo orbits in the elliptic restricted three-body
  problem," AAS/AIAA Space Flight Mechanics Meeting, Galveston, Texas, 2008.
- [25] Peng, H., and Xu, S., "Stability of two groups of multi-revolution elliptic halo orbits in the elliptic restricted three-body
  problem," Celestial Mechanics and Dynamical Astronomy, Vol. 123, No. 3, 2015, pp. 279-303.
- [34] Peng, H., Bai, X., and Xu, S., "Continuation of periodic orbits in the Sun-Mercury elliptic restricted three-body problem,"
  Communications in Nonlinear Science and Numerical Simulation, Vol. 47, 2017, pp. 1-15.
- [26] Jorba, A., and Villanueva, J., "On the persistence of lower dimensional invariant tori under quasi-periodic perturbations,"
  Journal of Nonlinear Science, Vol. 7, No. 5, 1997, pp. 427-473.
- [27] Gomez, G., and Mondelo, J. M., "The dynamics around the collinear equilibrium points of the RTBP," Physica D: Nonlinear
  Phenomena, Vol. 157, No. 4, 2001, pp. 283-321.
- [28] Olikara, Z. P., and Scheeres, D. J., "Numerical method for computing quasi-periodic orbits and their stability in the
  restricted three-body problem," Advances in the Astronautical Sciences, Vol. 145, 2012, pp. 911-930.
- [30] Olikara, Z. P., Gomez, G., and Masdemont, J. J., "A note on dynamics about the coherent Sun-Earth-Moon collinear libration
  points," Astrodynamics Network AstroNet-II: The Final Conference, Springer, 2016, pp. 183-192.
- [32] Rosales, J. J., Jorba, A., and Jorba-Cusco, M., "Families of halo-like invariant tori around L2 in the Earth-Moon bicircular
  problem," Celestial Mechanics and Dynamical Astronomy, Vol. 133, 2021, p. 16.
- [33] Jorba-Cusco, M., Farres, A., and Jorba, A., "Two periodic models for the Earth-Moon system," Frontiers in Applied Mathematics
  and Statistics, Vol. 4, 2018, pp. 1-14.
- [31] Seydel, R., Practical Bifurcation and Stability Analysis, Vol. 5, Springer Science & Business Media, 2009.
- [29] McCarthy, B. P., "Cislunar trajectory design methodologies incorporating quasi-periodic structures with applications," Ph.D.
  Dissertation, Purdue University, West Lafayette, Indiana, 2022.

Model-transition (CR3BP to ephemeris) studies:
- [5], [7] Dei Tos, D. A., and Topputo, F., "Trajectory refinement of three-body orbits in the real solar system model," Advances
  in Space Research, Vol. 59, No. 8, 2017, pp. 2117-2132 (listed twice).
- [6] Davis, D. C., Phillips, S. M., Howell, K. C., Vutukuri, S., and McCarthy, B. P., "Stationkeeping and transfer trajectory design
  for spacecraft in cislunar space," 2017 AAS/AIAA Astrodynamics Specialist Conference, 2017, pp. 1-20.
- [8] Oguri, K., Oshima, K., Campagnola, S., Kakihara, K., Ozaki, N., Baresi, N., Kawakatsu, Y., and Funase, R., "EQUULEUS trajectory
  design," The Journal of the Astronautical Sciences, Vol. 67, No. 3, 2020, pp. 950-976.
- [9] Boudad, K., Howell, K. C., and Davis, D. C., "Analogs for Earth-Moon halo orbits and their evolving characteristics in
  higher-fidelity force models," AIAA SCITECH 2022 Forum, 2022.
- [10] Park, B., and Howell, K. C., "Leveraging intermediate dynamical models for transitioning from the circular restricted
  three-body problem to an ephemeris model," AAS/AIAA Astrodynamics Specialist Conference, Charlotte, North Carolina, 2022.
- [11] Park, B., and Howell, K. C., "Leveraging the elliptic restricted three-body problem for characterization of multi-year
  Earth-Moon L2 halos in an ephemeris model," AAS/AIAA Astrodynamics Specialist Conference, Big Sky, Montana, 2023.
- [12] Gomez, G., Masdemont, J., and Mondelo, J. M., "Solar system models with a selected set of frequencies," Astronomy &
  Astrophysics, Vol. 390, No. 2, 2002, pp. 733-749.
- [13] Lian, Y., Gomez, G., Masdemont, J. J., and Tang, G., "A note on the dynamics around the Lagrange collinear points of the
  Earth-Moon system in a complete solar system model," Celestial Mechanics and Dynamical Astronomy, Vol. 115, 2013, pp. 185-211.

Others cited: [2] Crusan et al. 2018 (Gateway); [3] McComas et al. 2018 (IMAP); [4] Clampin 2008 (JWST); [17] Park, Folkner,
Williams and Boggs 2021 (DE440, DE441), The Astronomical Journal 161(3):105; [19] Wiesel and Pohlen 1994 (canonical Floquet
theory); [20] Williams, Howell and Davis 2023; [21] Lujan and Scheeres 2022 (Earth-Moon L2 quasi-halo family); [22] Olikara 2016
(Ph.D. dissertation, Colorado).
