# Digest: Jorba & Villanueva 1997, "On the persistence of lower dimensional invariant tori under quasi-periodic perturbations"

Citation: A. Jorba and J. Villanueva, "On the Persistence of Lower Dimensional Invariant Tori under
Quasi-Periodic Perturbations", *Journal of Nonlinear Science* 7:427-473 (1997), DOI
`10.1007/s003329900036`; received 21 May 1996, revised 8 January 1997, communicated by S. Wiggins. 47 PDF
pages (journal pp. 427-473). Filed in the private paper corpus as
`jorba-villanueva-1997-persistence-lower-dimensional-invariant-tori-quasi-periodic-perturbations-jns-7-427-doi-10.1007-s003329900036.pdf`
(md5 `d1b140fb64b37bd82fb9b42655c44e91`). Good text layer; the statement of Theorem 1 (p. 440) was also read from
the page image.

Markers: PRINTED = taken from the paper (page given); DERIVED = worked out here, not printed; COMPUTED = computed
this session; PROJECT = a statement about this repository. The paper has NO tables, NO printed numerical results and no
figures of data: the two applications (Section 4) are qualitative, so there is nothing in it that can be a numerical test.
Sections 2 to 5.4 (the proof) were read for statements and hypotheses; the technical lemmas (1 to 16, pp. 445-473) were
read in outline only and are not digested line by line.

Companions: `2026-10-04-digest-rhouma-chicone-2000-continuation-periodic-orbits.md` (commensurate orbits stay periodic) and
`2026-10-04-digest-rosales-jorba-jorba-cusco-2021-bcp-halo-like-tori-l2.md` (the numerical families of 2-tori).

## 1. What the paper does

Setting (PRINTED, p. 427-431). An autonomous analytic Hamiltonian H with l degrees of freedom has an r-dimensional
invariant torus (r = 1 is a periodic orbit) with a quasi-periodic flow of intrinsic frequencies omega_hat^(0) in R^r. A
perturbation eps H_hat is added that is analytic and depends on time quasi-periodically with s basic frequencies
omega_tilde^(0) in R^s (s = 1 is a periodic forcing such as the Sun). Result: under analyticity, nonresonance and
nondegeneracy, most of the tori survive as (r + s)-dimensional tori that have the old intrinsic frequencies plus the
perturbation's frequencies ("adding the frequencies of the perturbation to the ones they already have", abstract).
The new content over earlier work is (a) the case of lower dimensional tori (r < l) including normally ELLIPTIC ones, and (b)
exponentially small bounds on the measure of the destroyed set. Normally hyperbolic initial tori are said to be easier and
not treated explicitly (p. 430, p. 441).

The hard point is the "lack of parameters" problem: for a normally elliptic torus the small divisors involve both the intrinsic
frequencies and the NORMAL frequencies, and one cannot prescribe both. The device (p. 430, 437): let the normal frequencies
move as a function of eps (result a) or of the intrinsic frequencies (result b), and discard a Cantor set of resonant parameter
values.

## 2. The Hamiltonian form and the hypotheses

Normal form (PRINTED, Eqs. 1, 3, 4, pp. 427, 434). Angles theta = (theta_hat, theta_tilde) in R^(r+s), actions I = (I_hat, I_tilde),
normal variables z = (x, y) (m pairs, r + m = l):

    H(theta, x, I, y, eps) = omega^(0)T I + (1/2) z^T B z + (1/2) I_hat^T C I_hat + H_*(theta_hat, x, I_hat, y) + eps H_hat(theta, x, I_hat, y, eps),

with omega^(0) = (omega_hat^(0), omega_tilde^(0)) in R^(r+s), H_* of order at least 3 in the sense (monomial degree |l|_1 + 2|j|_1) with
the seminormal-form conditions P1, P2 (the coefficients of (z, I_hat), (z, I_hat, I_hat) vanish; the (z, z, I_hat) and (I_hat, I_hat) coefficients are
independent of theta_hat, the first only the trivial resonant ones). Those are reached by three Lie-series steps, which need only the
Diophantine condition below. The term in I_tilde is the device for making a time-dependent perturbation autonomous.

**Theorem 1 (PRINTED, p. 440).** Hypotheses:

- (i) H_* and H_hat are analytic in (theta, x, I_hat, y) near z = 0, I_hat = 0, 2 pi-periodic in theta, for every eps in I_0 = [0, eps_0], on a domain
  independent of eps; the dependence on eps is C^2 and the eps-derivatives of H_hat are analytic on the same domain.
- (ii) B is a symmetric constant matrix such that J_m B is diagonal with DIFFERENT eigenvalues lambda = (lambda_1, ..., lambda_m, -lambda_1, ..., -lambda_m)
  (the torus' normal flow is reducible to constant coefficients; for a periodic orbit this is Floquet, p. 431).
- (iii) NDC1: C is a symmetric constant matrix with det C != 0 (the intrinsic frequencies depend on the actions; equivalently
  det d^2 H_0 / d I^2 != 0, p. 428).
- (iv) Diophantine (nonresonance) condition: for some mu_0 > 0 and gamma > r + s - 1,

      | i k^T omega^(0) + l^T lambda | >= mu_0 / |k|_1^gamma,    k in Z^(r+s) \ {0},   l in N^(2m),  |l|_1 <= 2.

  (|k|_1 = |k_1| + ... + |k_(r+s)|.) For purely elliptic normal directions lambda_j = i nu_j, so the divisors are k.omega +- nu_j and k.omega +- nu_j +- nu_l,
  and k.omega alone for l = 0. Remark (p. 441): condition (iv) holds for all frequencies and eigenvalues except a set of measure zero.

Then, "under certain generic nondegeneracy conditions" NDC2:

- NDC2 (PRINTED, p. 466, Eq. 33 and following). After one normal-form step in eps, the eigenvalues of J_m B^(1) are written
  lambda_j(phi) = lambda_j + i u_j eps + i v_j^T (omega_hat - omega_hat^(0)) + (higher order), phi = (omega_hat, eps), u_j in C, v_j in C^r.
  NDC2: for every j with Re lambda_j = 0 (elliptic), u_j != 0 and Re(v_j) is NOT in Z^r; and the same two conditions for the differences
  u_(j,l) = u_j - u_l, v_(j,l) = v_j - v_l for every j != l with Re(lambda_j - lambda_l) = 0.
  Informal reading (p. 440, "Note"): the normal frequencies must depend on eps and on the intrinsic frequencies of the basic family of tori
  (the family of Section 2.5, parametrised by omega_hat through the change I_hat -> I_hat + C^-1 (omega_hat - omega_hat^(0)), Eq. 14).

Conclusions:

- (a) There is a Cantor set I_* in I_0 such that for every eps in I_* the Hamiltonian has a reducible (r + s)-dimensional invariant torus with the
  frequency vector omega^(0) (the intrinsic frequencies are kept fixed, the torus is deformed). For every 0 < sigma < 1 and eps_bar small enough
  (depending on sigma):

      mes([0, eps_bar] \ I_*(eps_bar)) <= exp( -(1/eps_bar)^(sigma/gamma) ).

- (b) Given R_0 > 0 small enough and a FIXED eps with 0 <= eps <= R_0^(gamma + 1) [the exponent on R_0 is small and hard to read in the print: read as
  gamma + 1], there is a Cantor set W_*(eps, R_0) inside the ball V(R_0) = { |omega_hat - omega_hat^(0)| <= R_0 } such that for every omega_hat in W_*
  the Hamiltonian at that fixed eps has a reducible (r + s)-dimensional invariant torus with frequency vector (omega_hat, omega_tilde^(0)). For every
  0 < sigma < 1 and R_0 small enough:

      mes( V(R_0) \ W_*(eps, R_0) ) <= exp( -(1/R_0)^(sigma/(gamma + 1)) ).

  Both measures are Lebesgue measure; the constants (the smallness thresholds "eps_bar small enough depending on sigma") are not given.

Remarks (PRINTED, pp. 440-441). (b) at eps = 0 says: around an r-dimensional reducible torus of the unperturbed system there is an r-dimensional
Cantor family of r-dimensional reducible tori parametrised by omega_hat, destroyed fraction exponentially small in R_0. The same holds around any
(r + s)-torus obtained for some eps != 0 provided its frequencies satisfy the same Diophantine bounds. If the perturbation is autonomous (s = 0) the
results are estimates on the measure of destroyed tori of an autonomous Hamiltonian near a lower-dimensional torus.

How the measure bound arises (PRINTED, Section 2.4, p. 437-438). With normal eigenvalues moving by at most a eps, once the Diophantine constant has
dropped to mu_n <= mu_0/2, only resonant k with |k|_1 >= K(eps) = (mu_0/(2 a eps))^(1/gamma) matter: "we do not have low order resonances nearby;
we only have to eliminate higher order ones". The removed measure is exponentially small, of order exp(-1/eps_0^c) for any 0 < c < 1/gamma.

Proof method (PRINTED). A Kolmogorov-type quadratically convergent scheme: at each step a canonical change by the time-one flow of a generating function
S (Eq. 7) kills the terms a_tilde, b, c - omega_hat^(0), E and the non-diagonal part of B (Eqs. eq1 to eq5, p. 436), divisors i k.omega^(0), i k.omega^(0) + lambda_j,
i k.omega^(0) + lambda_j + lambda_l. Dependence on parameters is only Lipschitz (Section 5.1.1), so the tori depend on (eps, omega_hat) in a Lipschitz way.

## 3. Applications printed in the paper (Section 4)

**4.1 Bicircular model near L4,5 (pp. 442-444, qualitative).** The Hamiltonian in synodic coordinates is printed (p. 442):

    H = (1/2)(PX^2 + PY^2 + PZ^2) + Y PX - X PY - (1 - mu)/r_PE - mu/r_PM - m_s/r_PS - (m_s/a_s^2)(Y sin(theta) - X cos(theta)),

with theta = omega_S t, PX = Xdot - Y, PY = Ydot + X, PZ = Zdot, Sun at X_s = a_s cos(theta), Y_s = -a_s sin(theta), and r_PS^2 = (X - X_s)^2 + (Y - Y_s)^2 + Z^2.
DERIVED check against the project: the last term is the indirect term, and the acceleration from this Hamiltonian is -m_s (r - r_S)/|r - r_S|^3 - m_s r_S/a_s^3, identical
to `core/bcr4bp.py::_sun_acceleration`; the Sun position (a_s cos(theta), -a_s sin(theta)) with theta = omega_S t is the project's corrected clockwise sense (`#891`).
Parameter eps multiplies the perturbation: eps = 0 is the RTBP, eps = 1 the bicircular model with real values; numerical constants m_s, a_s, omega_S are not printed here.

Statements: for small eps the equilibrium L4,5 becomes a periodic orbit with the Sun's period; the three families of Lyapunov periodic orbits (short period, long
period, vertical) become three Cantorian families of 2-D invariant tori (adding the Sun's frequency); the 2-D Lyapunov tori (products of two periodic families)
become 3-D tori provided they are nonresonant with the perturbation; the maximal 3-D tori become 4-D (the last "already contained in [24]", Jorba & Simo 1996). "eps = 1 is
too big to apply these results; in particular eps is big enough to cause a change of stability in the periodic orbit that replaces the equilibrium point"; for eps = 1 one
must first normal-form around the (unstable) periodic orbit, after which Theorem 1 gives tori of dimensions 1, 2, 3 in its central directions. The stable region off the
plane (centred on vertical Lyapunov orbits) is only said to be suggested by numerics ([17], [34]); existence "has not been proved rigorously" there.
4.1.1: for the non-circular (ER3BP-type) models the equilibrium is replaced by a quasi-periodic solution that exists only for a Cantor set of eps.

**4.2 Halo orbits (p. 444).** Earth-Sun RTBP near L1: halo orbits (normal behaviour centre x saddle). With the real motions of Earth, Moon, Venus etc. modelled as a quasi-periodic
perturbation with r > 0 frequencies and a parameter eps in front: for small eps the halo orbits become a Cantorian family of (r + 1)-D invariant tori with the same centre x saddle
normal behaviour. For eps = 1 refer to the bicircular remarks.

Paper's own remark on checking hypotheses (p. 441, Section 4 opening): "The nondegeneracy conditions can be checked numerically, observing if the frequencies involved depend on the
parameters. As the applications ... are perturbations of families of periodic orbits, this hypothesis can be easily verified computing the variation of the period along the family
as well as the eigenvalues of the monodromy matrix. The numerical verification of the Diophantine condition is more difficult, but ... most of the orbits of the family are going to
satisfy it."

## 4. Relation to Rhouma & Chicone 2000 and to Rosales et al. 2021

**Complementary, on the frequency-ratio axis.** Let rho = omega_S / omega_orbit = T_orbit / T_S for an orbit of the autonomous CR3BP forced with period T_S (DERIVED notation; omega_orbit = 2 pi / T0).

| | Rhouma & Chicone 2000 | Jorba & Villanueva 1997 |
|---|---|---|
| Frequency condition | exactly commensurate, M T0 = N T_S (k.omega = 0 for a nonzero integer k) | k.omega != 0 for all k, with Diophantine lower bound |
| Result | periodic orbit of period N T_S, at simple zeros of the Melnikov function (Jacobi-constant work integral) | invariant (r + s) = 2-torus with frequencies (omega_orbit, omega_S), for most frequencies near, Cantor set |
| Normal-direction conditions | the other Floquet multipliers lambda must have lambda^M != 1 (kernel dimension) | normal eigenvalues reducible and different; Diophantine involving normal frequencies; NDC2 (normal frequencies move with eps and omega_hat) |
| Frequency nondegeneracy | dT0/dC != 0 (kernel of Phi(NT) - I is one dimensional) | NDC1: det C != 0, which for r = 1 is d omega_hat / d I != 0, i.e. again dT0 / dC != 0 along the family |
| Forcing | ONE period | s frequencies, any s |
| Radius in eps | none given | none given (constants unspecified) |

The same nondegeneracy, a nonvanishing derivative of the period along the family, is needed by both. They exclude each other on a given orbit: JV needs k.omega != 0 for every k (an orbit with
T0 = (N/M) T_S violates it at k = (M, -N)), RC needs the opposite. A family of CR3BP orbits parametrised by Jacobi constant C crosses a commensurate value at isolated C: at such orbits
the surviving objects are the periodic orbits of RC (isolated, at Melnikov zeros), at the irrational values the surviving objects are JV tori, and in between lie the resonance gaps, in which
JV excludes a set of frequencies whose measure is exponentially small but which is DENSE (every rational). The paper does not describe the object inside a gap; the standard picture (islands
around the Melnikov-zero periodic orbits, with width proportional to the square root of the forcing amplitude times the Melnikov amplitude over the twist) is background and not from either paper.

**Rosales et al. 2021 (CMDA 133:16; digest `2026-10-04-digest-rosales-jorba-jorba-cusco-2021-bcp-halo-like-tori-l2.md`).** Their numerics compute families of invariant curves of the stroboscopic map
F (flow over T_S) in the bicircular problem, parametrised by the rotation number rho (their Eq. 3), and report "families are Cantorian (gaps at resonances)" with gaps passed by lowering the Sun's mass
(p. 13 of that paper); "most of the tori are hyperbolic". That is the numerical manifestation of JV conclusion (a) / (b): for r = 1 the periodic orbit that replaces L2 is the unperturbed orbit, the families of
planar and halo-like tori are the 2-tori obtained by adding the Sun's frequency, the resonance gaps are the Cantor-set removals, and the commensurate members are the periodic orbits of RC (for example the
epsilon = 0 orbit at period T/2 in their digest, section 3). JV is the theorem behind their numerics at small epsilon, with the same limits: it says nothing about epsilon = 1.

## 5. What it says concretely about Earth-Moon CR3BP orbits under the bicircular Sun

Take f = CR3BP (planar or spatial, `core/cr3bp.py`), the perturbation the Sun of `core/bcr4bp.py` with mass scaled by eps (eps = 1 physical), s = 1, omega_tilde = omega_S = 0.925196 (rate in the synodic
frame; Sun period T_S = 6.7912 TU, the project's `BCR4BPSystem.sun_period_tu`). A CR3BP periodic orbit with synodic period T0 and frequency omega_hat = 2 pi / T0 is an r = 1 torus.

- DERIVED (JV, small eps). If dT0/dC != 0 at the orbit (NDC1), all Floquet multipliers other than the Jordan pair are distinct and normal-reducible (an orbit in a family whose multipliers do not collide;
  the paper's condition (ii) needs distinct eigenvalues, so multiplier collisions and bifurcation points are excluded), NDC2 holds (the normal exponents depend on eps and on omega_hat with u_j != 0, Re v_j not integer),
  and the Diophantine condition (iv) holds, then the orbit continues as an invariant 2-torus with frequencies (omega_hat, omega_S) for eps in a Cantor set, and for fixed small eps there is such a torus for all omega_hat in
  a Cantor set of full relative measure up to an exponentially small remainder. Orbits with a ratio omega_S / omega_hat that is a rational with small numerator and denominator are not covered: they are RC orbits.
- DERIVED. Hyperbolic normal directions (halos, Lyapunov orbits of L1 and L2, centre x saddle) need only the elliptic directions' conditions (NDC2 is stated only for Re lambda_j = 0); the paper says the hyperbolic case is easier and gives better results (p. 441).
- The printed claims about the bicircular problem are only the qualitative ones of section 3 (Lyapunov families of L4,5 to Cantorian 2-D tori; Lyapunov tori to 3-D; halo orbits to (r + 1)-D tori); L1 and L2 of the Earth-Moon system are not treated explicitly; they are the same theorem applied to the families of `core/cr3bp.py` periodic orbits.
- The measure statement is asymptotic (eps or R_0 "small enough depending on sigma", constants unspecified): it gives no number for eps = 1 or for the solar coupling mu_S / a_S^3 = 5.60e-3 (COMPUTED in the Rhouma digest). The paper itself says eps = 1 is too big for L4,5. Whether a given Earth-Moon torus exists at the physical Sun is a numerical question (Rosales et al. do this for a few families).
- The QBCP (coherent, `core/qbcp.py`) Sun is also periodic with the Sun's period in the synodic frame, so s = 1 with the same theorem; an ER3BP or ephemeris forcing adds frequencies (s >= 2) and the theorem still applies provided the extended Diophantine condition in Z^(r+s) holds.

## 6. Techniques applicable to the project's problems

**Common numerical test (applies to all items below), DERIVED from the hypotheses.**
For a CR3BP periodic orbit of period T0 forced with T_S, compute:
1. Frequency nondegeneracy: dT0/dC along the family (same quantity as the Rhouma screen) at the orbit.
2. Floquet multipliers of the one-lap monodromy matrix: all distinct (condition ii), and for the normal elliptic ones their arguments nu_j T0 (mod 2 pi).
3. Small-divisor table: for k = (k1, k2) with |k1| + |k2| <= K, the numbers |k1 omega_hat + k2 omega_S|, |k1 omega_hat + k2 omega_S +- nu_j|, |k1 omega_hat + k2 omega_S +- nu_j +- nu_l|; report the minimum of divisor times |k|_1^gamma for gamma = 1.1 and 2 (condition (iv) needs a constant mu_0 > 0 that does not shrink with K; with finite K only a trend can be read, not a proof). Any divisor under, say, 1e-3 at modest |k|_1 flags a resonance gap.
4. NDC2 by finite differences: normal exponents lambda_j(eps, omega_hat) from invariant-circle computations at several eps and several nearby members of the family; u_j = d lambda_j / d eps, v_j = d lambda_j / d omega_hat; require u_j != 0 and Re v_j not an integer.
5. Existence by direct computation, which is the real test: an invariant circle of the stroboscopic map F (flow over T_S) with the rotation number T_S/T0 mod 1, by the Fourier-series continuation of Rosales et al. / Jorba 2001. PROJECT: `search/pertbp_strob_889.py` has rotation-number and invariant-circle machinery for the stroboscopic map; `search/sun_forced_periodic_884.py` has the commensurate (periodic) counterpart.

**`#905` (rerun `#884` in the corrected models).** The `#884` stored orbits are the COMMENSURATE case (RC), not JV. What JV adds for `#905`: (a) a family walk in Jacobi constant that passes a commensurate member should show invariant circles on both sides; on a walk that tracks the orbit through eps from 0 to 1, a loss of the circle near a low-order rational is the resonance gap and is expected, not a bug; (b) the rule for which members can be trusted as "persisting": nonresonant with Diophantine margin (JV) or commensurate with a simple Melnikov zero (RC); anything else (small divisors but neither) is unclassified.

**`#902` (rebuild the bicircular validation tiers on a real orbit).** A new check tier per orbit: "regime" (commensurate / Diophantine / gap) from the divisor table, and, for the Diophantine regime, an invariant-circle residual, so that a persistence claim is attached to one of the two theorems. JV supplies no printed number, so the sourced controls remain the Jorba-Jorba-Cusco-Rosales 2020 L1 orbit and the Rosales 2021 families; the Hamiltonian printed on p. 442 is a source for the Sun terms and agrees with the code (section 3), which can be recorded as an independent (a)-type check of the `#891` sense and indirect term.

**`#916` (test the printed persistence conjecture of Ross & Roberts-Tsoukkas as invariant curves of the one-period map).** JV is the relevant theorem for near-commensurate STABLE (normally elliptic) members: the claim that a stable near-commensurate cycler persists as an invariant curve is the existence of a JV torus, and JV says: only for frequencies outside a measure-small set of resonance gaps, only at small forcing, and only if the normal-elliptic conditions (iv) and NDC2 hold. So the conjecture as printed cannot be true for every near-commensurate stable member at the physical forcing; a clean test is, for each stable member, to compute the divisor table and the invariant circle and report the three regimes. A member that sits in a gap is a counterexample only to the "invariant curve" form, not to persistence of the periodic orbits (RC).

**`#922` (the `#890` orbit as an invariant two-torus once one more frequency is added).** JV allows s >= 2 frequencies (the RC theorem does not). The relevant numbers (COMPUTED from the registry's printed periods, 8.705869 d Titania, 13.463237 d Oberon): Titania-Oberon synodic period 24.637 d; ratio of Oberon's anomalistic frequency to the synodic frequency 1.830; for the five-synodic-period `#890` orbit (intrinsic period about 123.2 d, `#890` note), the ratio of Oberon's anomalistic frequency to the orbit's frequency is 5 x 1.830 = 9.150, a ratio with no small-denominator rational neighbour (9.150 = 183/20, order |k|_1 = 203). Caveats stated plainly: (1) the `#890` orbit is itself commensurate with the Oberon-synodic forcing (period exactly 5 periods), so it is the RC kind of object for that forcing and is NOT an r = 1 torus with a nonresonant frequency vector; the JV statement applies to the frequency vector (omega_hat, nu, n_O) only after the commensurate relation is treated as the unperturbed periodic orbit of the extended system, and then the closure of the extended action (the Melnikov-type condition) must be checked, which the paper does not address; (2) the paper's NDC1 requires a family of tori with frequency varying with the action, which for a periodic orbit of the time-dependent two-moon model must be built in the extended phase space; I did not work this out. The reading that is safe: adding the eccentricity frequency to a commensurate periodic orbit gives, for the rotation number 9.150 mod 1 = 0.150, a candidate invariant circle of the stroboscopic map; whether it exists is the JV (a) question and needs the numerical test of item 5 above, as `#922` already plans. This answers `#907`'s torus question only in that form: the orbit is a point inside a resonance gap of the Titania-only torus family, not a piece of a nonresonant torus.

**`#926` (Sanaga & Howell 2025 follow-ups), item (b): the multiplier root-of-unity monitor.** A multiplier at a p-th root of unity is the l = 1 or l = 2 case of the divisors in condition (iv) with k = 0 (for the repeated orbit: lambda^p = 1 means a normal frequency times p is a multiple of the intrinsic one). JV shows the monitor should watch three things, not one: the rotation number against rationals (RC resonance), the normal multiplier arguments against rationals (the |k| = 0, l = 1 divisors), and the combinations k.omega +- nu_j and k.omega +- nu_j +- nu_l with k != 0 (mixed divisors). The positive control the task asks for (a CR3BP family with a known period-multiplying bifurcation) tests the first two; the mixed divisors have no published control in the corpus.

## 7. Follow-ups (task numbers not registered here)

1. Implement the divisor table (section 6, item 3) as a function of (omega_hat, omega_S, normal exponents) and run it on every stored `#884` and Jorba-Jorba-Cusco-Rosales family member; record regime and margin per row.
2. Add an NDC2 check (item 4) for the orbits where an invariant circle is computed; the first test case is the Rosales 2021 L2 planar Lyapunov family.
3. For `#916`: classify each stable near-commensurate member as periodic-orbit-persists (RC), JV torus exists (circle computed), or gap; compare with the printed conjecture.
4. For `#922`: decide the correct unperturbed object for the extended problem (the Titania CR3BP family plus Oberon circular), and compute the rotation-number 0.150 invariant circle of the stroboscopic map as the first step.
5. Acquire the nearest sources that complete the picture, if absent from the corpus: Jorba & Simo 1996 (SIAM J. Math. Anal. 27:1704, "On quasi-periodic perturbations of elliptic equilibrium points", the paper's reference [24]); the Jorba 2001 invariant-curve algorithm cited by Rosales et al. (not checked against the corpus); the Rosales 2023 QBCP paper is held and digested.
6. Record that no source states an eps-radius: neither paper gives a number at which the theorems stop applying, so "survives the physical Sun" is always a computed claim.

## 8. Reading quality

Text layer good; Theorem 1 re-read from the page image (p. 440). Uncertain digit: the exponent of R_0 in the restriction on eps in conclusion (b) is a small superscript; read as gamma + 1. The proof sections (5.2 to 5.4 and Lemmas 1 to 16) were skimmed for hypotheses and constants only.
