# Digest: Cincotta and Simo 1999, conditional entropy as a tool to explore phase space

Date: 2026-10-04 (Sydney). Reading and reasoning, plus one throwaway numerical check of the recipe (section 7.4) that was
run outside the repository and is not committed. No existing file was edited except the corpus index row.

Source: P. Cincotta and C. Simo, "Conditional entropy: a tool to explore the phase space", Celestial Mechanics and
Dynamical Astronomy 73:195-209 (1999), DOI 10.1023/A:1008355215603. Filed in the private paper corpus as
`cincotta-simo-1999-conditional-entropy-tool-explore-phase-space-cmda-73-195-doi-10.1023-A-1008355215603-image-only.pdf`
(15 pages, pp.195-209, image-only scan with no text layer). I rendered every page at 130 dpi and read all 15 page images,
including the figures. The equations below were read from the page images; the few that are fine print are flagged.

Evidence tags: READ (p.N) is read in the paper at printed page N. COMPUTED is a calculation of mine. INFERRED is my
reasoning, not stated by the authors. "Read off the figure" means an approximate reading of a plotted curve, not a printed
number.

Related digests: `docs/notes/2026-10-04-digest-cincotta-simo-2000-megno.md` (MEGNO, the same authors' later indicator) and
`docs/notes/2026-10-04-digest-froeschle-lega-gonczi-1997-fli.md` (FLI). The task `#924` entry in `data/OUTSTANDING.md`
(chaos indicators as a screen) is the consumer of all three.

## 1. What the paper is

A short methods paper in galactic-dynamics style, not orbital mechanics. It proposes an indicator built from the
information-theoretic entropy of two nearby orbits, with the ARC LENGTH along each orbit as the random variable instead of
time (READ, abstract p.195). The concept was introduced by Nunez, Cincotta and Wachlin (1996), who used it in two- and
three-dimensional systems; this paper supplies theory that the 1996 paper lacked and a worked example, the Henon-Heiles
system (READ p.196). The full theory is deferred to Cincotta and Simo "1998, preprint" (listed on p.209 and cited again on pp.207-208 for the convergence of J/2T and the information in J^(0), J^(1); not held in the corpus, and not identified further here). The abstract's claims
(READ p.195): the tool localises stable periodic and quasi-periodic orbits, chaotic ("aperiodic") motion and unstable
periodic orbits ("the source of chaotic motion"); it gives a measure of chaos "similar to that given by the largest
Lyapunov Characteristic Number" and "does not require long time integrations, just realistic physical times".

Section layout: 1 introduction (p.195); 2 set-up, definitions and analytic results (pp.196-201); 3 numerical examples on
Henon-Heiles (pp.202-207, Figs 1-7); 4 conclusions (p.208); references (p.209). There is no table.

## 2. The method, exactly as printed

### 2.1 Setting (pp.196-197)

Hamiltonian H(p, q) = p^2/2 + phi(q), with p, q in R^N and phi smooth (eq. 1). State x = (p, q) in R^(2N); the vector field
is v = (-dH/dq, dH/dp) = (-grad phi, p) (eq. 2), so x-dot = v(x) (eq. 3). The energy surface M_h = {x : H = h} is assumed
compact. An orbit arc is gamma = {x(t; x0) : 0 <= t <= T < T_*} (eq. 4), with T_* bounding the total motion time. Its length
L_gamma = integral of ds (eq. 5), where "it depends on the metric: in this case ds^2 = dx_1^2 + ... + dx_2N^2", i.e. the
Euclidean metric on the 2N phase-space coordinates (p, q) with no scaling between position and momentum.

### 2.2 Arc-length averages and the entropy of one orbit (pp.196-197)

Average of a scalar psi along the orbit by arc length: <psi>_gamma = (1/L_gamma) integral psi(chi(s)) ds (eq. 6). Because
ds = |x-dot| dt = |v| dt (eq. 7), L_gamma = T <|v|>_T and <psi>_gamma = <psi |v|>_T / <|v|>_T (eq. 8), where <.>_T is the
finite-time average and |.| the Euclidean norm. With psi = ln|v(x)| the ENTROPY OF THE ORBIT is

    S(gamma) = - <ln|v|>_gamma + ln L_gamma                                  (eq. 9)

which equals the continuous Shannon form

    S(gamma) = - integral_0^T rho_gamma(t) ln rho_gamma(t) dt,   rho_gamma(t) = |v(x(t))| / L_gamma >= 0   (eq. 10)

with integral_0^T rho_gamma dt = 1 (READ p.197). So the orbit is turned into a probability density over its own time
parameter, proportional to the SPEED |v| = |x-dot| in phase space, and S is the entropy of that density. A harmonic
oscillator has constant |v| and rho = 1/T (INFERRED from the printed v = <v> = sqrt(2h), p.201).

Binning (READ p.197): the continuous entropy S_c = -integral rho ln rho dx (eq. 11) is the limit of the discrete S_d = -sum
mu(x_i) ln mu(x_i) (eq. 12) over a partition, "just as it is done in the Riemann sums". No bin count or width is given
anywhere in the paper and none is needed: all computations use the closed forms (27) and (28) below (READ p.202). The
binning in the task description therefore does not occur as a numerical step; it is the theoretical route to the definition.

### 2.3 Conditional entropy of two sets (pp.197-198)

For a measurable space M with probability measure mu, sets A = {a_i} and B = {b_j}: S(A|B) = sum_j mu(b_j) S(A|b_j) = -sum_j
mu(b_j) sum_i mu(a_i|b_j) ln mu(a_i|b_j) = S(A, B) - S(B), with S(A, B) the joint entropy. The SYMMETRIC conditional
entropy

    I(A, B) = (1/2) [S(A|B) + S(B|A)] = S(A, B) - (1/2) [S(A) + S(B)]       (eq. 13)

is >= 0, equals zero if and only if mu(a_i) = mu(b_i) (READ p.198, citing Arnold and Avez 1989), tends to 0 when A and B are
strongly correlated and to [S(A) + S(B)]/2 when they are independent. In the continuous case a constant "S_0 = infinity" is
dropped and positivity is kept because I is a relative entropy (p.198).

### 2.4 Applied to two nearby orbits (p.198)

gamma is the reference orbit; gamma' is an orbit on the same or a nearby energy surface M_h' with x0' = x0 + delta_0 and
vector field v' = v(x'). The joint object is the curve Gamma = gamma x gamma' in M_h x M_h' defined by V = (v, v'), with
|V|^2 = |v|^2 + |v'|^2, and rho_Gamma = |V(x(t), x'(t))| / L_Gamma. Then

    I(gamma, gamma') = S(Gamma) - (1/2) [S(gamma) + S(gamma')]              (eq. 14)

The physical argument (READ pp.198-199): for x0 in the regular region Sigma_r, gamma and gamma' stay close, diverging at most
linearly, so they stay strongly correlated and I is approximately 0. For x0 in the chaotic region Sigma_c and T greater than
a critical time T_c, they diverge exponentially, decorrelate, and I > 0. T_c "can be associated to the time needed to pass
close to an hyperbolic object with non coincident separatrices". For T less than T_c, I is about 0 for both kinds of orbit.

### 2.5 Second-order expansion and the key formulas (pp.199-201)

Write x'(t) = x(t) + delta(t), v(t) = |x-dot|, delta(t) = |delta|, d(t) = |delta-dot|, xi = d/v assumed small, and
zeta the angle between x-dot and delta-dot. Up to second order in xi (READ p.199, eq. 15-17):

    rho_gamma' ~ rho_gamma (1 + a1 + a2 - k1 a1),  rho_Gamma ~ rho_gamma (1 + a1/2 + a2/4 - k1 a1/4 + a3)    (eq. 15)
    a_j(t) = Delta_j(t)/v(t) - <Delta_j>_T/<v>_T  (j = 1, 2, 3),   k1 = <Delta_1>_T/<v>_T                    (eq. 16)
    Delta_1 = Delta v = (delta-dot . x-dot)/v = d cos zeta;  Delta_2 = (1/2)(d^2/v) sin^2 zeta;  Delta_3 = d^2/(8 v)  (eq. 17)

The paper states that rho_gamma' and rho_Gamma differ from rho_gamma, at first order, "by the instantaneous fluctuation of
Delta v(t), the first variation of |v(x(t))|". After algebra (READ p.199, eq. 18):

    I(gamma, gamma') ~ <v a1^2>_T / (8 <v>_T)                                                  (eq. 18)

COMPUTED from eq. 16: v a1^2 = Delta_1^2/v - 2 k1 Delta_1 + k1^2 v, so <v a1^2>_T = <Delta_1^2/v>_T - <Delta_1>_T^2/<v>_T,
which is the closed form I used in section 7.4. If <v>_T depends only mildly on T, <v>_T ~ X (eq. 19, p.200):
I ~ (<Delta_1^2>_T - <Delta_1>_T^2)/(8 X^2), "the variance of Delta_1". Delta_1 follows from the variational equations
delta-dot^k = (dv^k/dx^l) delta^l (eq. 20) as Delta_1 ~ (1/X) Phi_l delta q^l with Phi_l = phi_k phi_kl - phi_l (eq. 21),
which are the components of the gradient of |grad phi|^2/2 - phi, a function of position only (INFERRED from the printed
definition).

Regular motion (READ p.200, eqs 22-23): with delta q^l(t) ~ delta_0 lambda_l t (linear divergence, lambda = max lambda_l),
Delta_1 ~ (delta_0 lambda t / X) Psi(q(t)) and, provided (omega . k) T is large for all integer vectors k not equal to 0,

    I(gamma, gamma') ~ (delta_0^2 lambda^2 A / (24 X^4)) T^2                                   (eq. 23)

with A the (almost T-independent) mean of Psi^2. The numerical factors 8X^2 in (19) and 24 X^4 in (23) are fine print in
the scan; I read them as printed but checked only the T^2 dependence (section 7.4). Hence I(T; x0) depends on x0 mainly
through lambda(x0) and I_1/I_0 ~ lambda_1^2/lambda_0^2, smallest "for x0 in a neighbourhood of a stable periodic orbit"
(READ p.201). Valid for lambda T << X/delta_0, so I ~ O(delta_0^2 T^2) << 1.

The second indicator, J. Using dI/dT ~ (1/(8T)) (a1^2(T) - 8 I(T)) v(T)/<v>_T (p.201) the paper defines

    J(gamma, gamma') = d log I / d log T = (T/I) dI/dT ~ ( a1^2(T) / (8 I(T)) - 1 ) v(T)/<v>_T ~ Delta_1^2(T)/<Delta_1^2>_T - 1   (eq. 24)

(READ p.201; the last form holds for v(t) ~ <v>_T ~ X and <Delta_1>_T = 0). Limits (READ p.201):

    regular: J ~ 2 for T >> T_D, independent of delta_0, T, lambda; if delta grows like t^r then J ~ 2r        (eq. 25)
    chaotic, delta ~ delta_0 e^(sigma t): I ~ (delta_0^2/X^4) e^(2 sigma T)/(sigma T),  J ~ 2 sigma T - J_0, J_0 ~ 1,  T >> T_D   (eq. 26)

COMPUTED check of eq. 26: ln I = 2 sigma T - ln(sigma T) + const, so d ln I / d ln T = T (2 sigma - 1/T) = 2 sigma T - 1, which
matches J ~ 2 sigma T - J_0 with J_0 = 1. For eq. 25: I proportional to T^2 gives d ln I / d ln T = 2. So J is the LOCAL
LOGARITHMIC SLOPE of I against T: it is not a time average of anything, it is a derivative, so it carries the oscillations
of the orbit (paper, p.205: "while I is a time-averaged quantity, J is not").

### 2.6 Relation to the largest Lyapunov characteristic number

READ (p.207): "recalling the second of (26) we see that dJ/dT ~ LCN". So for chaotic motion J grows linearly with T and the
slope of J(T) is 2 sigma; the paper fits the linear part of J^(2)(T) by least squares, calling the slope sigma, and
states "As we are computing J^(2) instead of J^(0), the factor 2 in front of the second of (26) compensates in part the
averaging procedure" (p.207). The conclusion states the aim as "the true convergence of J/2T to the LCN" (p.208), which
confirms that the estimate is LCN approx J/(2T) (equivalently half the slope of J). The paper gives no formula for the
slope fitting window, and calls the fitted value sigma_E (Fig. 6b, p.206). Compare MEGNO, where the working rule is
LCN approx J/T because that J is defined with a different weight (MEGNO digest section 2). The two indicators are
different quantities; both equal 2 for regular motion by different routes.

### 2.7 How it is computed (p.202-204)

For the numerical work the entropy of an orbit is written in closed form

    S(T; x0) = -(1/L(T; x0)) integral_0^T |v| ln|v| dt + ln L(T; x0)                           (eq. 27)

with L(T; x0) the time integral of |v| (READ p.202, "the time-integral of the vector field"), and the entropy rate

    dS/dT = rho(T; x0) (1 - S(T; x0)) - rho(T; x0) ln rho(T; x0)                               (eq. 28)

so that dI/dT = dS(Gamma)/dT - (1/2)[dS(gamma)/dT + dS(gamma')/dT] and J = T (dI/dT)/I. Three computational routes (p.202-204):

- A1: integrate both x0 and x0' = x0 + delta_0 with the full equations (3), then use (27), (28). Used by Nunez et al.
  (1996). Restriction: delta(t) = |x'(t) - x(t)| reaches a saturation value set by the size of the system.
- A2: approximate x'(t) = x(t) + delta(t) with delta from the variational equations (20), with delta_0^k = delta_0/2,
  delta_0 = 1e-6, k = 1..4 (note the two orbits are then on different, very close energy levels, p.203-204). No saturation,
  but numerically unstable for delta_0 -> 0 or delta large; I and hence J depend on delta_0, though for 1e-7 <= delta_0 <= 1e-4
  J is "almost independent" of delta_0 and I scales as delta_0^2 (p.204). The integration is stopped when I becomes negative,
  which happens for chaotic orbits once delta reaches the system size (p.205).
- A3: use the second-order formulas (18) and the first of (24) directly. Needs no explicit gamma', is "numerically stable",
  and is used for all computations of I and J in Figs 3-7 (p.204). delta_0 is then simply a scale factor set (arbitrarily)
  to 1, with I ~ O(T^2) >> 1 for regular motion and I ~ O(exp(sigma T)) >> 1 for chaotic (p.204; the printed
  "exp(sigma T)" is as I read it, eq. 26 gives exp(2 sigma T)).

Smoothing (READ p.204): J depends on the instantaneous rho(T), so it is averaged: with T_n = T_0 + n deltaT,
J^(0)(T_n) = J(T_n), J^(1)(T_n) = (1/n) sum_(k<=n) J(T_k), J^(2)(T_n) = (1/n) sum_(k<=n) J^(1)(T_k): a repeated running mean "to
lower as much as possible the effect of fast oscillations that do not affect the aperiodical changes". The sampling step
deltaT is not stated in the paper (READ absence); the MEGNO paper's later deltaT = 0.06 is for another indicator and is not
carried over.

Cost (READ p.202): in all numerical integrations T < T_* = 10^3 T_D, about 10^4 time units for Henon-Heiles (T_D < about 10).
The LCN comparison used T = 2.5e5 (p.202). The paper says the energy-surface variable is not evaluated for rho_Gamma ("We did
not plot rho_Gamma because its behaviour is almost the same than rho_gamma'", p.204). The extra cost over a plain
variational run is two or three scalar accumulators (the closed forms), the same order as MEGNO; no explicit second orbit
is needed in route A3 (INFERRED from the formulas).

## 3. Printed test setup and results

All items are for the Henon-Heiles model (Henon and Heiles 1964); the potential itself is not printed in this paper and is
the standard phi = (q1^2 + q2^2)/2 + q1^2 q2 - q2^3/3 (INFERRED; the paper cites the model only by reference, so this form
is recorded as the standard one, not as READ). The energy is h = 0.118 (READ p.202).

### 3.1 Phase-space description of the test (READ p.202, Fig. 1a on p.203)

- Surface of section (q2, p2), p1 > 0, at q1 = 0 (Fig. 1a caption), initial conditions q1_0 = 0 and (q2, p2)_0 on the line
  p2_0 = 0 with 0.3 <= q2 <= 0.62.
- For 0.3 <= q2 <~ 0.508 a large regular component; q2 about 0.305 corresponds to the stable 1-periodic orbit.
- For 0.508 <~ q2 <~ 0.57 a 5-periodic island with its stochastic layer around the separatrix, bounded by a KAM curve that
  "disappears for h = 0.119".
- For q2 >~ 0.57 a highly stochastic region, with many small islands throughout.

### 3.2 LCN reference run (Fig. 1b, p.203; text p.202)

1000 initial conditions along the q2 axis, 0.55 <= q2_0 <= 0.60 (p2_0 = 0), T = 2.5e5. Plotted as log10 (LCN) (axis from about
-1 to -5, read off the figure; the figure caption says "log (LCN)" without a base, and the base-10 reading follows from
the axis range and the other plots, INFERRED). The plot shows plateaus of order 1e-1 to 1e-2 (the large stochastic sea,
read off the figure) and floor values near 1e-5 for the regular components. It carries no information on the regular
structure (p.207).

### 3.3 Orbit evolution (Figs 2 and 3, pp.203-204)

Five orbits at p2_0 = 0, labelled by q2_0 (READ p.204, p.203 caption):

| q2_0 | Character (paper) | Behaviour |
| --- | --- | --- |
| 0.305 | near the stable 1-periodic orbit | rho_gamma'/rho_gamma = 1 to the resolution of Fig. 2 |
| 0.5 | quasi-periodic, associated with the 1-periodic orbit | rho ratio oscillates about 1 with very small amplitude |
| 0.5085 | quasi-periodic, close to the separatrix of the 5-periodic island | small periodic pulses and a drift in rho ratio (Fig. 2b, window [0.999995, 1.000015]) |
| 0.509 | inside the stochastic layer | rho ratio diverges from 1 after "several periods"; time-rate of exponential divergence slower than for q2_0 = 0.6 |
| 0.6 | in the highly stochastic sea | rho ratio diverges "after a few periods" |

Fig. 2a (A2 procedure, rho_gamma'/rho_gamma vs T, T up to about 1500; the A2 runs for the chaotic orbits end when I turns
negative, T about 1e3 or less, p.205). In Fig. 3 (A3, T up to 3000): log I against T for the five orbits (3a), and J^(2) against
T (3b). Printed statements (READ pp.204-205): "log I has a logarithmic dependence with T (in fact, I proportional to T^2) for x0
in Sigma_r while it varies nearly linear for x0 in Sigma_c"; the smallest I is nearest the stable 1-periodic orbit; for
x0 in Sigma_r and T >~ 10^2 T_D, J^(2) is about 2, and for x0 in Sigma_c with T > T_c, J^(2) depends linearly on T; for the
quasi-periodic orbit near the separatrix (q2_0 = 0.5085) J^(2) "reaches a maximum value and then it decreases
asymptotically to 2", the maximum being the first of the pulses of Fig. 2b (p.205). Fig. 3 values read off the figure:
log I runs from about 0 to 20 over T = 0 to 3000 for q2_0 = 0.6 (largest) and stays below about 6 for the regular
orbits; J^(2) for q2_0 = 0.6 passes 10 well before T = 1000 (off the top of the plot window 0 to 10). These are
approximate, in the 0.05 relative range at best.

Statement of the margin (READ p.205): "the regular and the stochastic component are separated by several orders of
magnitude (see the scale used in the I axis)" (Fig. 5 context, p.207).

### 3.4 Scans across the island (Figs 4-5, p.206)

- Fig. 4: 3000 initial conditions on the q2 axis within the 1-periodic island, 0.3 <= q2_0 <= about 0.5, T = 4000; (a) log I,
  (b) J^(2). log I has its minimum at the stable periodic orbit and rises fast approaching the separatrix; J^(2) is about 2
  across the island (read off the figure, the left plateau about 1.75 to 1.8, rising slightly to the right, with sharp
  spikes at several q2_0 values, the first near 0.38; the paper reads these spikes as thin chaotic layers around
  high-order-resonance separatrices, p.205: "'discontinuities' ... reveal the existence of thin chaotic layers").
- Fig. 5: the same for 4000 initial conditions across the 5-periodic separatrix, 0.506 <= q2_0 <= 0.514, T = 4000: (a)
  log (log I), (b) log J^(2) with "log 2 approx 0.3" marking the regular level (base-10 logs). The stochastic layer is
  about 0.5085 to about 0.511 (read off the figure), clearly separated from the regular sides.

### 3.5 Measure of chaos and LCN comparison (Fig. 6, p.206-207)

Fig. 6 uses the same 1000 initial conditions as Fig. 1b but T = 7000: (a) log J^(2); (b) log sigma_E, the slope-derived estimate
(section 2.6). READ p.207: J^(2) is about 2 in the ordered component and clearly larger in the stochastic regions, higher in
the large stochastic sea than in the stochastic layer. "A comparison of this figure with that for the LCN (Fig. 1b)
reveals a good agreement between both magnitudes." The cost comparison, READ p.207: sigma_E "was computed for T = 7000, the
LCN was for T = 2.5 x 10^5. Actually, T = 5000 is enough to get sigma_E while the total motion time used to compute the LCN
is, perhaps, not sufficient to get a good asymptotic value", and sigma_E "shows the structure of the regular component while
Fig. 1(b) does not". That is a motion-time saving of about 36 (2.5e5 / 7000) to 50 (2.5e5 / 5000), COMPUTED; the
conclusion says "two or three orders of magnitude" (p.208), which is not what these two numbers give, so the safer
statement is the one with the printed times. No numerical value of sigma_E or LCN is printed anywhere (READ absence); the
agreement is a visual one.

### 3.6 Locating the unstable periodic orbit (Fig. 7, pp.207-208)

Fig. 7a: log J^(2)_max (the maximum of J^(2) over T) against q2_0 for 4000 initial conditions near the separatrix,
0.506 <= q2_0 <= about 0.514. The regular side shows slight sensitivity to small islands, a smooth continuum that "diverges"
at the separatrix; the stochastic side shows the same continuum plus many sharp lines. READ p.207: the maximum of J^(2)
"is due to the presence, somewhere, of this unstable orbit" (the 5-periodic hyperbolic orbit). Fig. 7b: compute T_m with
J^(2)(T_m) = J^(2)_max for each regular initial condition, integrate each to T about T_m, and plot (q2, p2) when q1 = 0,
p1 >= 0, for 1300 initial conditions in 0.506 <= q2_0 <= 0.50859 (caption: "0.50859"). The points cluster near the hyperbolic
points of the 5-periodic orbit (Fig. 7b shows clusters near q2 about 0.2 to 0.5 and p2 up to about +/- 0.1, read off the
figure). No printed coordinates of the hyperbolic point.

## 4. Stated advantages and limits

Authors' conclusions (READ p.208): the conditional entropy defined through the arc length "is an efficient tool to explore
the phase space in short motion times"; "provides an easy way to find out the location of unstable periodic orbits";
"preliminary analytical results agree with that obtained by numerical simulations (at least for very simple systems)"; the
mean exponential divergence rate can be estimated "for motion times which are two or three orders of magnitude less than
that for the computation of the LCN using the standard procedures". Still to be done (p.208): "i) a more complete theory;
ii) numerical study of 3D models and iii) application to realistic dynamical systems".

Limits visible in the paper:

- Theory is second order in xi = d/v and rests on (a) <v>_T roughly independent of T, (b) regular divergence linear in T, and
  (c) a Fourier-expansion non-resonance condition on (omega . k) T (pp.199-201); the authors call it "preliminary".
- Valid only for lambda T << X/delta_0 in the regular case (p.201), and J_0 ~ 1 in the chaotic case is an unspecified
  constant (p.201).
- Route A2 breaks once delta is of order the size of the system (p.204-205); route A3 hides this because delta_0 is a free
  scale, but a chaotic orbit with e^(sigma T) growth overflows double precision after T about 700/sigma (COMPUTED: 1e308 =
  e^709).
- Everything is on one potential system, 2D, bounded, smooth, energy-conserving, along a single line of initial conditions;
  one reference LCN image and no tabulated values. No threshold is derived for "regular"; the criteria are J^(2) about 2
  against "much larger" (p.207), and the I-value separation of "several orders of magnitude" (p.207).
- The sampling step deltaT of the smoothing is not stated, and the J^(2) smoothing hides the pulses of an orbit near a
  hyperbolic orbit (p.205: "only the first peak, for small T, is significant").
- Per the MEGNO paper (p.226 there): the conditional entropy "looks noisy and depends more strongly on the time step" when
  many high-order resonances are present, and is "more sensitive to periodic orbits" (MEGNO digest section 4).
- No mention of non-autonomous Hamiltonians, rotating frames, time-dependent potentials, singular potentials or close
  encounters (READ absence).

## 5. Integration times needed (all READ)

| Quantity | Time | Page |
| --- | --- | --- |
| Reference LCN, 1000 orbits, 0.55-0.60 | T = 2.5e5 | pp.202-203 |
| Fig. 3 time series (I, J^(2)) | T up to 3000 | p.204 |
| Fig. 2a (A2 procedure) | up to about 1500, stopped when I becomes negative | pp.203, 205 |
| Figs 4-5 scans (3000 and 4000 orbits) | T = 4000 | p.205 |
| Fig. 6, J^(2) and sigma_E for the 1000 orbits | T = 7000 (T = 5000 "enough" for sigma_E) | pp.206-207 |
| T_* general bound | 10^3 T_D, about 10^4 for Henon-Heiles (T_D <~ 10) | p.202 |
| Time for regular J^(2) to reach 2 | T >~ 10^2 T_D | p.205 |
| Standard LCN, general statement | 10^5 to 10^6 T_D | p.195 |

## 6. Printed numbers usable as sourced tests

The paper has no table and no printed entropy or exponent values. What is printed and usable:

1. Model, energy, section: Henon-Heiles at h = 0.118, section q1 = 0, p1 > 0, line p2 = 0 (p.202, Fig. 1a caption p.203).
   The potential form is the standard one and must be sourced to Henon and Heiles 1964 (not held; INFERRED).
2. Orbit classification at h = 0.118, p2_0 = 0: q2_0 about 0.305 is the stable 1-periodic orbit, 0.3 to 0.508 the regular
   1-periodic island, 0.508 to 0.57 the 5-periodic island and its layer, above 0.57 the stochastic sea (p.202). Five test
   points with a printed character: 0.305, 0.5, 0.5085 (regular, near separatrix), 0.509 (in the layer), 0.6 (stochastic)
   (p.204).
3. J^(2) tends to 2 for regular motion (eq. 25; text p.205) and to a linear growth 2 sigma T for chaotic motion (eq. 26). The
   regular asymptote is independent of delta_0, T and lambda; this is a scale-free test. For delta ~ t^r, J ~ 2r (p.201).
4. I proportional to T^2 for regular motion; I roughly exponential (e^(2 sigma T)/(sigma T)) for chaotic motion (eqs 23, 26,
   p.205).
5. For 0.5085, J^(2) rises to a maximum and then falls to 2 (p.205); the maximum identifies the hyperbolic orbit (p.207).
6. The comparison with the LCN: the sigma_E map at T = 7000 agrees (visually) with the LCN map at T = 2.5e5 (pp.206-207), and
   T = 5000 is enough (p.207). The expected-side is a plot, not numbers; any test must use a tolerance defined by us and say
   so.
7. delta_0 independence: for 1e-7 <= delta_0 <= 1e-4, J is almost independent of delta_0 and I scales as delta_0^2 (p.204).
8. Chaotic-layer separation from the regular sea: "several orders of magnitude" in I (p.207); the 1-periodic island scan
   with sharp spikes at thin layers (p.205).

What is NOT printed, so cannot be a golden: any value of I, S, J_max, sigma or LCN at a named q2_0; T_m for any orbit; the
sampling step; the Henon-Heiles equations; the hyperbolic-orbit coordinates.

## 7. Reconciliation with project code, comparison with FLI and MEGNO, and techniques

### 7.1 Project code (COMPUTED by search on 2026-10-04)

A search of `src/cyclerfinder` for `megno`, `fast.lyapunov`, `fli` and `conditional.entropy` finds nothing; the only
exponent-like code is the discrete-QR `lyapunov_exponents` field of `search/ccr4bp_whisker.py` and
`search/variational_qbcp_torus.py` (Lyapunov exponents of a known torus from its monodromy, not an indicator on a grid).
The remaining "lyapunov" hits are planar-Lyapunov periodic orbits of the CR3BP libration points (a different object).
No chaos-indicator code exists; this paper changes nothing about that. The state and tangent matrix it needs already exist
in `core.cr3bp.cr3bp_stm_eom` and `core.ccr4bp.ccr4bp_stm_eom` (see the MEGNO digest section 7.1 for the wrapper).

### 7.2 Does it add anything over FLI and MEGNO

Plainly: very little, and as a first screen nothing. The comparison, same purpose (separating regular, chaotic and unstable
periodic motion on a grid of initial conditions in a project model):

| Property | Conditional entropy (this paper) | MEGNO (2000) | FLI (1997) |
| --- | --- | --- | --- |
| Quantity | J = d log I / d log T of an entropy-difference I | mean of a time-weighted log-derivative of delta | time at which the tangent-vector norm passes a threshold |
| Regular value | J about 2 (I proportional to T^2) | Jbar about 2 (printed band 1.98 to 2.035) | no absolute value; threshold times separate in rank |
| Chaotic value | J about 2 sigma T (slope 2 sigma) | Jbar slope sigma / 2 (J = sigma T) | threshold reached early |
| Needs | the variational equations, 3 accumulators (v, Delta_1, Delta_1^2/v), no extra orbit in route A3 | variational equations, 2 accumulators | tangent norm only, 0 accumulators |
| Absolute scale | yes (J about 2) | yes (Jbar about 2) | no |
| Sensitivity to the field-norm | high: uses |x-dot| = full phase-space vector field norm, with the speed of the acceleration components entering | none beyond tangent norm | tangent-vector norm only |
| Smoothing step | J^(1), J^(2) repeated means, step not given | continuous mean exists | none |
| Evidence | one 2D potential, no printed values | one 2D potential, printed band 1.98-2.035 | two maps + asteroid model, printed orders of magnitude |
| Documented noise | "noisy and depends more strongly on the time step" (MEGNO paper p.226) | smooth | smooth |

Why it adds nothing as a first or second screen (INFERRED unless tagged):

- It has the same regular value 2 as MEGNO and nearly the same cost, but its definition ties the result to the norm of the
  FULL vector field. In the project's rotating-frame models, with Coriolis terms and a gravity that spikes at a moon,
  |x-dot| = |f(x)| is dominated by the acceleration components at a flyby; the arc-length density rho = |f|/L then
  concentrates at the flyby. Assumption (a) of the theory (<v>_T ~ X, mild T-dependence, p.200) is exactly what fails.
  MEGNO and FLI have no such assumption.
- J is a logarithmic derivative, not an average, so it carries every flyby spike. The smoothing J^(2) is needed and its step
  is unspecified; MEGNO's mean is the better-behaved equivalent.
- The quantity (Delta_1) is the first variation of the SPEED of the orbit, which ignores perturbations that change the
  direction of f without changing |f|; for a perturbation that preserves speed exactly, Delta_1 = 0 to first order and
  the indicator is blind at that order. That is a property of the formula (COMPUTED from Delta_1 = f . (A delta)/|f|), not
  something the paper discusses.
- Its distinctive claim, sensitivity to unstable periodic orbits through the J^(2) maximum, is also shown by MEGNO (Fig. 7
  there, orbit B) and, for hyperbolic structure generally, by FLI.
- The MEGNO paper by the same authors lists it among the competing tools and calls it noisier (section 4 of the MEGNO digest).

What it does add, if anything: the logarithmic-slope diagnostic J = d ln I / d ln T is a cheap local exponent that can be
computed from any scalar growth record, and the power-law reading (delta ~ t^r gives J = 2r, eq. 25) distinguishes linear
growth from other subexponential power laws, which MEGNO's "a = 0, b = 2n" also does. Not new.

### 7.3 One-off check (COMPUTED; throwaway script in the session scratch area, not committed)

Setup: Henon-Heiles, h = 0.118, q1 = 0, p2 = 0, p1 > 0 from the energy; the standard potential above; DOP853 with
rtol = atol = 1e-12; a random unit delta_0 in R^4 (not the paper's delta_0), unnormalised (delta_0 is a scale); the closed form
I = (<Delta_1^2/v>_T - <Delta_1>_T^2/<v>_T) / (8 <v>_T^2) with the accumulators evaluated at 40 equally spaced T up to 3000
and J from the finite-difference slope of ln I against ln T on that 40-point grid; no J^(1) or J^(2) smoothing. Results at
T = 3000 (the J column is the mean of the last three grid values, which is noisy, so only the first decimal is meaningful):

| q2_0 | Paper character | I(3000) | J (last 3 grid values) | log10 of the tangent norm |
| --- | --- | --- | --- | --- |
| 0.305 | stable 1-periodic | 3.1e1 | 1.978, 1.978, 1.978 | 1.75 |
| 0.5 | regular | 7.0e2 | 1.99, 1.97, 1.93 | 2.46 |
| 0.5085 | regular near separatrix | 1.8e5 | 1.63, 1.30, 2.99 | 2.80 |
| 0.509 | stochastic layer | 6.3e36 | 38, 7.0, 9.4 (noisy) | 19.9 |
| 0.6 | stochastic sea | 1.8e139 | 433, 373, 347 | 71.6 |

Reading (COMPUTED, my numbers, not the paper's): the regular orbits give J about 2 (0.305 gives 1.98, matching eq. 25; 0.5
gives 1.9 to 2.0), the 0.6 orbit gives J of order 350 which matches 2 sigma T with sigma = 71.6 ln 10 / 3000 = 0.055 (2 sigma T
= 330, within the noise of the coarse grid), and the I ratio between 0.6 and the regular orbits is more than 130 orders of
magnitude. So eqs 18, 25 and 26 hold in this model. Not matching cleanly: 0.5085 (J swings between 1.3 and 3.0 at the end,
consistent with the paper's pulses and drift, but not a clean 2 at 40 sample points) and 0.509 (J noisy and an order of
magnitude below 2 sigma T = 92 for sigma = 0.0153 from the tangent norm at this T; the finite-difference slope on a
40-point grid cannot follow a quantity this spiky and the layer orbit may not have reached its asymptotic growth by
T = 3000, as the paper says layer divergence is slower). The five-point result is the closest I have to a positive control
for the code; the control is of the formula as I read it, and the formulas (19) and (23) use the printed constants unverified.
The check used a different random delta_0 from the paper's; nothing here is a published number.

## 8. Techniques applicable to the project's problems

### `#924` (chaos indicators as a screen)

Do not implement this indicator. The recommendation of the FLI digest (FLI first as the cheap screen, MEGNO to confirm)
stands. If a second absolute-scale confirmation is ever wanted beyond MEGNO, this one offers no advantage. The one
transferable idea is the diagnostic: record ln of a growth quantity against ln T and read the local slope, which tells
power-law growth from exponential growth in a single plot; this costs nothing with the accumulators MEGNO already needs and
should be added there as a plot, not as a new indicator.

### `#908` (capture sweeps, weak-stability sampling)

The sweep problem is a grid in speed, eccentricity and phase, with a time-to-escape criterion. FLI's time-to-threshold is the
closest analogue and was already recommended. This paper adds only the observation that a finite-T indicator needs a
reference case near a stable periodic orbit (minimum I, J about 2) and one in a known chaotic layer for calibration (Table in
section 3.3), which is generic practice. The arc-length entropy would break at a close pass for the speed-norm reason in
section 7.2; the sweeps live near a flyby, so not recommended there.

### `#890` and `#895` (Titania-Oberon, Uranian; methods only)

Both are non-autonomous, planetocentric flows with prescribed moons. The theory behind J requires an autonomous
energy-conserving Hamiltonian with <v>_T approximately independent of T (p.200), neither of which holds, and the project's
flybys at about 1000 km altitude produce the speed spikes noted above. For these two the Floquet analysis of the orbit
itself is what MEGNO would repeat; what is useful is classifying the NEIGHBOURHOOD, which MEGNO or FLI does with fewer
assumptions. No method from this paper is recommended; nothing in the paper concerns planetary satellites.

### `#916` and `#917` (persistence conjecture; Casoliva rows with the fold count)

`#916` tests whether stable near-commensurate cyclers persist as invariant curves under eccentricity and the Sun; the
natural diagnostics are the one-period-map invariant curve and an indicator on a grid around it. Here the paper's J, like
MEGNO, would show a regular island as J about 2 and its separatrix as a spike (Figs 4-5), which is the right qualitative
behaviour, but MEGNO has it with a printed band; use MEGNO. `#917` (elliptic problem, blocked on `#912`) has an
explicitly time-dependent flow, outside this paper's theory.

### Where it could matter at all

The only case in which the arc-length entropy is distinctive is an autonomous, smooth, bounded Hamiltonian flow where one
wants a single curve I(T) whose T^2 against exponential growth separates regular from chaotic over many orders of
magnitude of I (Fig. 5, "several orders of magnitude"). The project's smooth autonomous models (the planar CR3BP away from
the primaries, and the libration-point neighbourhoods) qualify; a scan of a known planar CR3BP family neighbourhood would
reproduce the paper's separation. That is a curiosity, not a need.

## 9. Follow-ups (no task numbers registered here)

1. Record the decision: not implemented, because MEGNO and FLI cover the use and this indicator's speed-norm definition is
   fragile near flybys; add the paper to the `#924` notes as "read, adds nothing".
2. If the MEGNO accumulators are built, also log ln|delta| and the local slope d ln|delta| / d ln T (cheap); this is the
   only idea from this paper worth carrying, and it is not novel to it.
3. If anyone wants the control anyway: reproduce Fig. 5 on Henon-Heiles at h = 0.118 with the standard Henon-Heiles potential
   (source it to Henon and Heiles 1964, acquire it) and 4000 orbits in 0.506 <= q2_0 <= 0.514 at T = 4000, using J^(2) with a
   stated sampling step; the expected side is the 2 versus "much larger" classification and the position of the layer
   (about 0.5085 to 0.511, read off the figure), not any number. A tolerance must be ours.
4. Identify and, if published, acquire the Cincotta and Simo "1998, preprint" cited here (convergence of J/2T to the LCN, and the
   information in J^(0) and J^(1)); the MEGNO digest section 4 mentions the authors' "1999, 2000" conditional-entropy
   references. Check `docs/notes/CORPUS_INDEX.md` first.
5. Acquire Nunez, Cincotta and Wachlin 1996, CMDA 64:43 (the original, with 3D examples), if the 3D (spatial CR3BP) case is
   ever wanted. Not needed for the current work.
6. Acquire Henon and Heiles 1964, AJ 69:73, so the test potential is sourced rather than recalled.
7. A possible honest test in the repository, if a MEGNO or FLI module is built, is the five-point Henon-Heiles classification
   of section 3.3 (0.305, 0.5, 0.5085 regular; 0.509 and 0.6 chaotic) at h = 0.118: the paper's own character labels are
   the expected side, independent of our code. The five labels are printed (p.204); the potential must be sourced.
8. Flag for the MEGNO digest: its statement that the conditional entropy "is more sensitive to periodic orbits" matches
   this paper's Fig. 7 (J^(2) maximum located the 5-periodic hyperbolic orbit), so nothing to correct there.

## 10. Transcription flags (illegible or doubtful digits)

- Eq. 19 prefactor 1/(8 X^2) and eq. 23 prefactor 24 X^4: fine print at 130 dpi; consistent with eq. 18 and a T^2 law, but
  not rendered at higher resolution. Not used in tests.
- Fig. 5 axis labels: the abscissa runs from 0.506 to 0.514 with ticks every 0.001; read at 130 dpi, last tick 0.514 doubtful.
- Fig. 7b caption upper limit "0.50859": read as printed.
- Fig. 4 and Fig. 5 ordinate values (e.g. J about 1.75 to 1.8 on the left of Fig. 4b, log J about 0.3 for regular motion in
  Fig. 5b): read off the figures; approximate.
- The text on p.204 prints "I ~ O(exp(sigma T))" for the chaotic case, while eq. 26 on p.201 gives e^(2 sigma T); I take
  eq. 26 as the derived one and treat the p.204 form as shorthand (it is consistent in an order-of-magnitude sense only).
- The ordering of A2 versus A3 details, for example "delta_0^k = delta_0/2, k = 1..4", is read as printed; the two orbits are
  on slightly different energy levels, as the text says.
