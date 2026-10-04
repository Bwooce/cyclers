# Digest: Froeschle, Lega and Gonczi 1997, fast Lyapunov indicators

Date: 2026-10-04 (Sydney). Reading and reasoning, plus throwaway numerical checks (section 7.4) run outside the
repository and not committed. No existing file was edited.

Source: C. Froeschle, E. Lega and R. Gonczi, "Fast Lyapunov indicators. Application to asteroidal motion", Celestial
Mechanics and Dynamical Astronomy 67:41-62 (1997), DOI 10.1023/A:1008276418601, received 10 June 1996, accepted 4 February
1997. Filed in the private paper corpus as
`froeschle-lega-gonczi-1997-fast-lyapunov-indicators-asteroidal-motion-cmda-67-41-doi-10.1023-A-1008276418601.pdf`
(22 pages, pp.41-62, text layer present; the formulas I quote were checked against a 300 dpi page image). I read every
page and the figure captions; the figures themselves are read only for the orders of magnitude stated below.

Evidence tags: READ (p.N) is read at printed page N. COMPUTED is my calculation. INFERRED is my reasoning.

An important note on naming. The paper's "fast Lyapunov indicators" are three quantities psi1, psi2, psi3, each an INVERSE
norm of tangent vectors, and the working indicator is the TIME at which psi3 falls below a threshold. The later and more
common form, FLI(t) = sup_j ||v_j(t)|| or its logarithm (the form the review digest quotes as "FLI(t) = sup_j ||v_j(t)||"),
is a variant that does not appear in this paper (READ absence). I describe the 1997 definitions and say where the later
form differs only from memory of the literature, flagged INFERRED.

## 1. What the paper is

A short methods paper from the Nice group, in three parts (READ, abstract p.41): "a very simple and fast method to separate
chaotic from regular orbits for non-integrable Hamiltonian systems", tested on the standard map and the Henon-Heiles
potential, where it "appears to be at least as sensitive as the frequency-analysis method", and then applied to three
asteroids. The motivation (pp.41-42) is the cost of computing Lyapunov exponents for thousands of asteroids, and the puzzle
of "stable chaos": 522 Helga has a Lyapunov time of only 6900 years yet is a permanent member of the belt (Milani and
Nobili 1992). The authors' remark that a short Lyapunov time and long-term stability can coexist matters for the project
(section 7.2).

## 2. The definition, exactly as printed

Idea (READ, p.42): when computing a Lyapunov exponent the tangent vectors are renormalised to avoid overflow; "The time
necessary to reach a given value, either of the length of any vector or of the angle between vectors is taken as an
indicator of stochasticity."

Given an n-dimensional basis V_n(0) = (v1(0), ..., vn(0)) of the tangent space and an initial condition P(0), "in order to
have homogeneous quantities we define the three indicators as functions of t" (p.43, eq. 1):

    psi1(t) = 1 / ||v1(t)||^n
    psi2(t) = 1 / ( product_{j=1..n} ||vj(t)|| )
    psi3(t) = 1 / ( sup_j ||vj(t)|| )^n

with v_j(t) obtained "by integration of the variational equations" (p.43). Here n is the dimension of the phase space.
Comments, READ:

- psi1 is one vector. psi2 is related to the volume spanned (the paper says it is related to the angle between vectors
  through volume conservation, with the angle "of the order of the inverse of the norm of a given vector to the power n",
  p.43).
- psi1 and psi2 depend on the initial basis; psi3 depends on it least (Fig. 4, p.45). "In the following we will take psi3
  as fast Lyapunov indicator so that the orientation of the initial V_n basis will no more affect the results" (p.45).
- The sup is over the whole basis, but "the vector which results to be the greatest one after a short time, for example
  10^2 iterations for mappings ... or about 1000 years for the case of asteroids ... still remains the supremum for longer
  times. Therefore it is enough to follow its evolution in order to get psi3" (p.43). So after a short initial phase one
  tangent vector suffices.

The usable quantity is a time (READ, pp.45, 50): N(epsilon) is the number of iterations (or the time) "necessary for psi3
to reach a threshold" epsilon, with a maximum run length N_m; "if the threshold is not reached in less than N_m iterations
we take N(epsilon) = N_m" (p.50). Thresholds used: 1e-8 (standard map, k = 0.3), 1e-10 (standard map, k = 1.3), 1e-20
(Henon-Heiles).

Behaviour, READ:

- Chaotic orbit: all indicators "drop sharply", for the standard map at k = 0.3 down to about 1e-20 in 200 iterations
  (p.44, Fig. 2a).
- Regular orbit: psi3 decays as a power law of time (Henon-Heiles, y = 0.2: "decreases with a power law ... within a time
  span of 10000", p.48), because the tangent vector grows linearly through the shear of the invariant torus. For the
  regular standard-map orbit at t = 200 psi3 is "about 1e-5, that is 15 order of magnitude greater than for the chaotic
  orbit" (p.44).
- Relation to the exponent (section 4, pp.52-56): renormalising whenever ||v|| reaches 10^m gives the Lyapunov exponent as
  lim (N m ln 10) / sum T_i, where T_i are the times between renormalisations. The "Local Lyapunov Times" (LLTs)
  tau_i = T_i / (m ln 10) have as their mean the Lyapunov time (eq. 4). N(epsilon) computed with psi1 satisfies
  N(epsilon) approx tau * m ln 10 (m = 4 in the test, epsilon = 1e-8 with n = 2), and in the strongly chaotic regime
  (k >= 2.3) the mean and dispersion of N(epsilon) and of the LLTs agree (Fig. 16, p.58). For the regular orbit the LLTs
  grow exponentially (regularly spaced on a logarithmic scale) so their mean goes to infinity and the exponent tends to
  zero (p.54).
- Chaos-meter: the authors intend later to use the value of psi3 reached in a given time "as a chaos-meter in order to
  compute the LCI in the shortest time"; for now "the value reached by psi3 for different orbits in a given time span gives
  a first indication of the relative strength of chaos among these orbits" (p.59).

## 3. How it is computed

Variational equations along the orbit; for a map the Jacobian is applied iteratively. Basis: an orthonormal basis at
t = 0; the paper rotates it by an angle theta to test dependence (Fig. 4). The integrator and tolerances for the
Henon-Heiles and asteroid runs are not stated in this paper (READ absence). Renormalisation is needed only to avoid
overflow (p.42); threshold crossing needs no renormalisation. The asteroid model is the Sun and the four giant planets
only (p.57).

## 4. Stated advantages and limits

Advantages, READ:

- Speed. For a thin chaotic layer of the standard map at k = 1.3 near a hyperbolic point of the 1/6 resonance, Laskar et
  al. (1992) report that the Lyapunov indicator needed "5 x 10^6 iterations ... to clearly detect the chaotic motion, while
  it was already visible with the frequency analysis with 20000 iterations"; psi3 reaches 1e-20 in 2000 iterations (p.50,
  Fig. 10). A cross-section of 1.223 > |y| > 1.220 at x = 0 with epsilon = 1e-10 and N_m = 2000 iterations shows all the
  features that the frequency analysis shows with 20000 iterations (Fig. 11, p.50). At a resolution of 1e-6 in y N(epsilon)
  shows chaotic layers "which are not clearly detected with the frequency analysis" (p.51).
- For Henon-Heiles at E = 0.08333, psi3 for a chaotic orbit reaches 1e-20 "within a time span of only 10^3, which
  corresponds to about one hundred of intersections with the surface of section", while the Lyapunov indicator needs about
  10^5 iterations, "which corresponds to 1000 T_L", and ten times more to be trusted (pp.48-49).
- Asteroids (section 5, pp.57-58): with the five-body model, 153 Hilda (very stable) has psi3 still above 1e-4 after 50000
  years, 944 Hidalgo (very chaotic) reaches 1e-8 in 2000 years. 522 Helga, the intermediate case: the text says its Lyapunov indicator after 50000 years is "between 1 and 2 order of
  magnitude smaller than Hidalgo" and that "the corresponding psi3 is 4 order of magnitude smaller after only 2000 years";
  against Hilda the indicator is of the same order but psi3 is "2 order of magnitude smaller even only after 2000 years",
  "a strong indication that the asteroid Helga is slowly chaotic" (p.58). The comparison with Hidalgo reads oddly as printed
  (a less chaotic asteroid should have the larger psi3), so I quote it without interpreting it (INFERRED misprint or
  different normalisation). The Lyapunov time of Helga is close to 1e4 years; its Lyapunov indicator reaches its limit after 1e5
  years and 1e6 years are needed for confidence, though Milani (1993) obtained a 20 to 30 percent estimate in 50000 years.
- Cost: psi3 needs the whole basis for a short time and then a single vector; psi1 needs one vector throughout.

Limits, READ or visible in the paper:

- "Of course our method is not able to distinguish between tori and islands (conversely to the frequency map) which are
  both invariant curves with zero LCE" (p.51). So it separates hyperbolic from non-hyperbolic, nothing finer.
- Dependence on initial direction: strong for psi1 and psi2 (Fig. 4); psi3 mostly removes it but not entirely. A residual
  dependence on the initial point of a regular orbit remains through the geometry of the torus, in the amplitude of
  the short-term variation about the mean power law (Fig. 6, p.47). For the chaotic orbit of the standard map the spread of
  N(epsilon) over 1000 points on one orbit is wide because the chaos is small and the orbit is affected by the
  hyperbolic point; still about 90 percent of chaotic N(epsilon) are more than an order of magnitude smaller than the regular
  orbit's, and the largest chaotic value is below the smallest regular one (Fig. 5, p.46).
- Thresholds and run lengths are chosen by hand (1e-8, 1e-10, 1e-20; N_m of 2000 or 1e4). The paper gives no rule.
- The link to the Lyapunov time (section 4) is claimed for strong chaos only: "for small values of the non linearity
  parameter k the situation is opposite due to the complexity of the topology (chain of islands, cantori etc.) of the
  chaotic zone and of the dynamics within it (sticking phenomena)" (pp.56-57).
- Evidence base: two textbook models with one or two degrees of freedom, and three asteroids in a smooth, nearly
  integrable five-body model. Nothing in the paper concerns close encounters, singular potentials or non-autonomous driving
  (READ absence).
- In the conclusion the authors list the method as complementary to frequency analysis, the sup-map methods and the helicity
  angle method (p.60), not a replacement.

## 5. Integration times needed (all READ)

| Case | Time to separate regular from chaotic | Page |
| --- | --- | --- |
| Standard map k = 0.3, orbit in the small chaotic zone | 200 iterations (psi3 to 1e-20); N(1e-8) histogram over 1000 starts of up to 1e4 iterations | pp.44-46 |
| Standard map k = 1.3, thin layer at the 1/6 resonance | 2000 iterations; the Lyapunov indicator needed 5e6, frequency analysis 2e4 | p.50 |
| Henon-Heiles E = 0.08333 | 1e3 time units (about 100 section crossings) | p.48 |
| Hilda, Hidalgo, Helga | 2000 years to flag chaos; 5e4 years for Hilda to remain above 1e-4 | pp.57-58 |
| Lyapunov indicator (for comparison) | about 1e5 iterations (1000 T_L) for Henon-Heiles; 1e5 to 1e6 years for Helga | pp.49, 58 |

## 6. Every worked example, in brief

- Standard map x_{i+1} = x_i + k sin(x_i + y_i), y_{i+1} = x_i + y_i (mod 2 pi) (eq. 2, p.43). k = 0.3: chaotic orbit
  at x0 = y0 = 0.001 inside the small chaotic zone about the hyperbolic fixed point at the origin; regular torus at
  x0 = 1, y0 = 0 (Figs 1-3); basis rotation test (Fig. 4); N(1e-8) histograms for a chaotic and a regular orbit (Fig. 5);
  second torus at x0 = 0.2 (Fig. 6).
- Henon-Heiles x'' = -x - 2xy, y'' = -y - x^2 + y^2 (eq. 3, p.48) at E = 0.08333 on the section x = 0, x-dot > 0: regular
  orbit y = 0.2, y-dot = 0, and a separatrix-like orbit y = 0.12 with homoclinic small chaotic zones (Figs 7-8).
- Sensitivity test against frequency analysis on the standard map at k = 1.3, orbit at x0 = 0, y0 = 1.2218 (the sign is
  lost in the text layer; the cross-section is quoted for 1.223 to 1.220 in magnitude) near the hyperbolic point of the 1/6
  chain (Figs 9-12, pp.50-54), and a Henon-Heiles cross-section at y-dot = 0 with epsilon = 1e-20, N_m = 1e4 (Fig. 13).
- Section 4 (Figs 14-16): norm history of one tangent vector renormalised at 10 (m = 1) for the torus and the chaotic orbit;
  histograms of LLTs; chi-square distance between spectra after n up to 1e7 iterations (eq. 5), converging within about
  1e6 iterations; mean and spread of LLTs and N(epsilon) against k.
- Asteroids (Figs 17-19): Hilda, Hidalgo, Helga, section 5.

## 7. Techniques applicable to the project's problems

### 7.1 Recipe

The project already integrates the tangent flow, so the indicator is almost free.

1. Integrate the project's state with the STM it already carries (`core.cr3bp.cr3bp_stm_eom`, `core.ccr4bp.ccr4bp_stm_eom`,
   a 6-state plus flattened 6 x 6 Phi; for `#895` build the tangent from `ForceModel.accel_grad`). With Phi(0) = identity the
   columns of Phi(t) are the paper's orthonormal basis v_1(t) ... v_n(t), n = 6 (or 4 for the planar problem).
2. psi3(t) = 1 / (max over the columns of ||Phi_j(t)||)^n, with the Euclidean norm of the column in non-dimensional
   variables. psi3 needs no extra ODE and no renormalisation until the norm threatens to overflow.
3. Choose a threshold epsilon and a maximum run time N_m, and report N(epsilon): the first time psi3 < epsilon, or N_m if
   never. Rank orbits by N(epsilon), or by psi3 at one fixed time.
4. Calibrate on the same model: run one case known to be regular and one known to be hyperbolic, and take epsilon between
   them. The paper gives no threshold rule (READ absence), so the calibration is the project's own.
5. A regular orbit shows psi3 proportional to t^(-n) (linear growth of the sup norm). An exponentially growing norm is the
   chaotic signature. Both can be read from the log-log slope instead of a threshold.
6. Use a grid of initial conditions, as the paper does for its cross-sections, and read the map: low N(epsilon) is hyperbolic
   or chaotic, high is regular.
7. Never read a short N(epsilon) as escape: the paper's own Helga example shows a short Lyapunov time with no macroscopic
   instability (section 7.2).

### 7.2 What it says about a hyperbolic cycler and about a quasi-cycler

A hyperbolic cycler. Every initial vector with a component along the unstable direction grows like the largest Floquet
multiplier, so the indicator signals hyperbolic behaviour from the first cycle, and does so without asking for the multiplier.
For the `#890` orbit (largest multiplier 8.4e5 per 123.162 d cycle, `data/OUTSTANDING.md` `#890`, taken per cycle, INFERRED)
a sup norm of 1e4 (epsilon = 1e-8 with n = 2) is reached after about ln(1e4)/ln(8.4e5) = 0.68 of a cycle (COMPUTED,
uniform growth assumed). The indicator therefore repeats what the Floquet analysis already gives for the orbit itself.

A quasi-cycler. The question is whether nearby motion is regular quasi-periodic or hyperbolic. The indicator can answer
that over a neighbourhood grid, with the limit that it cannot tell a torus from an island chain (p.51). It also cannot say
whether hyperbolic motion is dangerous: the paper's introduction presents "stable chaos" (Helga, Lyapunov time 6900 years,
permanent member; Milani and Nobili 1992) and the conjecture of Morbidelli and Froeschle (1996) that, in the Nekhoroshev
regime, diffusion times are exponentially long compared with Lyapunov times (pp.41-42). For the project this is a direct
caution: a positive indicator on a cycler-like or quasi-cycler orbit is expected and is not a disqualification. What
matters is relative strength across a family or epoch, which is exactly how the authors use it ("a first indication of the
relative strength of chaos among these orbits", p.59).

It cannot certify the invariant torus asked about in `#907` (a partially hyperbolic torus also has exponentially growing
tangent directions, INFERRED).

### 7.3 Close flybys

The paper has no close-encounter case (READ absence). A flyby multiplies the tangent norm by a factor g in a short time, so
the threshold crossing happens at the encounter that carries the cumulative product past 1/epsilon^(1/n) and N(epsilon) is
quantised by the encounter times (INFERRED). Across a flyby the resolved physical-time integration is enough; a threshold
test is not sensitive to the integration variable since only the norm matters, so, unlike a quantity weighted by time,
a regularised independent variable changes only the time at which the crossing is recorded (convert back to physical time).
For passes close enough to need regularisation, the tangent equations must be regularised consistently (INFERRED). The
mixed units of position and velocity in the norm change N(epsilon) by an encounter-dependent amount, so keep the same
non-dimensionalisation across a comparison.

### 7.4 The smallest reproducible positive control

The standard map k = 0.3 pair is the smallest control in the paper (a map, six lines of code, no integrator).

Printed values (READ):

1. k = 0.3, chaotic orbit x0 = y0 = 0.001: psi3 falls to about 1e-20 within 200 iterations (p.44; the figure reading is
   approximate).
2. k = 0.3, regular torus x0 = 1, y0 = 0: psi3 about 1e-5 at 200 iterations (p.44, stated in the text).
3. N(1e-8): the histogram for the regular orbit lies entirely to the right of the chaotic one, and about 90 percent of the
   chaotic N values are more than an order of magnitude smaller (p.46, Fig. 5; N_m = 1e4).
4. For the regular orbit psi3 follows a power law of time.

One-off check (COMPUTED; throwaway script in the session scratch area, not committed). Standard map eq. (2) with the
tangent map of Jacobian [[1 + k cos(x + y), k cos(x + y)], [1, 1]], orthonormal initial basis, psi3 = 1/(sup norm)^2:

| Case | psi3 at 100 | at 200 | at 1000 | N(1e-8) | Paper |
| --- | --- | --- | --- | --- | --- |
| chaotic, (0.001, 0.001) | 2.5e-9 | 3.2e-18 | 1.1e-50 | 38 | about 1e-20 at 200 |
| regular, (1, 0) | 1.5e-4 | 2.4e-5 | 1.8e-6 | 9412 | about 1e-5 at 200 |

The regular value agrees with the printed one in order of magnitude; the chaotic one reaches 1e-18 at 200 against the
printed 1e-20, so the gap between the two orbits is about 13 orders against the printed "15" (approximate figure reading on the
paper's side). The qualitative separation and the N(1e-8) ordering reproduce. Henon-Heiles at E = 0.08333, y = 0.2,
y-dot = 0 with the full 4 x 4 tangent flow, DOP853 at 1e-12: psi3 = 1.8e-3, 1.8e-7, 1.1e-11 at t = 100, 1e3, 1e4, so
psi3 falls by 4 orders per decade in time, the t^(-n) law with n = 4 for linear growth of the norm (COMPUTED; matches the
paper's "power law" in kind, not in printed numbers).

Two further paper cases did NOT reproduce with my initial conditions, and I make no claim either way: the standard map at
k = 1.3 with x0 = 0, y0 = 1.2218 (printed: 1e-20 after 2000 iterations; mine: after 81), and Henon-Heiles y = 0.12,
y-dot = 0 (printed: 1e-20 within 1e3; mine: 1.2e-8 at 1e3, no sign of chaos). The paper's initial conditions are rounded
and sit on very thin chaotic layers near a hyperbolic point, and the sign of y0 for the k = 1.3 case is lost in the text
layer, so I read this as a sensitivity of those particular points and not a failure of either side (INFERRED). They are
not usable as controls without the authors' unrounded values.

The goldens are the printed numbers above. A second control with no dynamics at all: a vector with norm linear in time must
give psi3 proportional to t^(-n); a norm exp(sigma t) must cross epsilon at t = ln(epsilon^(-1/n))/sigma. These closed
forms are independent of any project code (COMPUTED).

## 8. Comparison with MEGNO for the project, and which to implement first

Both indicators consume the same thing, the tangent flow the project already integrates, and both share the limits that
matter most here: neither separates a torus from an island chain, neither has been shown for close encounters or
non-autonomous driving, and for a hyperbolic cycler both just repeat the Floquet analysis of the orbit itself, so their use
is on grids of initial conditions or across families. They differ as follows (see
`docs/notes/2026-10-04-digest-cincotta-simo-2000-megno.md`).

- Output. MEGNO has an absolute scale: Jbar near 2 (printed band 1.98 to 2.035) for regular motion and linear growth for
  chaotic motion, so a map can be read without calibration. FLI (1997 form) has no absolute scale; the paper gives no
  threshold rule, and the indicator is meant to rank orbits against each other, so every new model needs a calibration pair.
- Cost and code. FLI is the sup of the STM column norms, with no extra differential equations, no weights, no renormalisation
  until overflow, and no long runs: the paper shows separation in a few hundred to 2000 map iterations or about 100 section
  crossings, and a time-to-threshold lets a run stop early. MEGNO needs two accumulators, renormalisation and about 500 to 1000
  characteristic periods.
- Time weighting. MEGNO's weight t makes it depend on when flybys occur (late encounters dominate); a threshold crossing
  of the norm depends only on the cumulative growth. For multi-flyby orbits the FLI behaves more simply.
- Evidence. The paper's FLI controls are on two maps and a reduced asteroid model, with printed orders of magnitude; the
  MEGNO paper has a flow with a printed numerical band that I reproduced in a one-off run.

Recommendation: implement FLI (psi3 and its time-to-threshold) first, as a screening pass over a grid, because it is a
ten-line addition on top of the STM columns the code already returns, it can stop each run early, and it is the natural
analogue of the project's present time-to-escape criterion in the stability sweeps (`#378`, `#681`, `#908`). Calibrate it
on a known regular and a known hyperbolic case in each model. Then implement MEGNO as the confirmation on the survivors,
since its absolute value 2 gives a model-independent regular band and is the stronger positive control. Build the shared
tangent-flow wrapper once. This is a judgement about effort and risk, not a measured comparison: I have not run both
indicators on one project orbit (INFERRED).
