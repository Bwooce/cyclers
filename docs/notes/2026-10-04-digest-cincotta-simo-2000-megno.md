# Digest: Cincotta and Simo 2000, simple tools to study global dynamics (the MEGNO paper)

Date: 2026-10-04 (Sydney). Reading and reasoning, plus one throwaway numerical check of the recipe (section 7.4) that was
run outside the repository and is not committed. No existing file was edited.

Source: P. M. Cincotta and C. Simo, "Simple tools to study global dynamics in non-axisymmetric galactic potentials - I",
Astron. Astrophys. Suppl. Ser. 147:205-228 (2000), DOI 10.1051/aas:2000108, received 30 July 1999, accepted 4 September
2000. Filed in the private paper corpus as
`cincotta-simo-2000-simple-tools-global-dynamics-non-axisymmetric-galactic-potentials-MEGNO-aas-147-205-doi-10.1051-aas-2000108.pdf`
(24 pages, text layer good). I read every page (pp.205-228), the figures at page-image resolution, and re-rendered the
MEGNO definitions (pp.214-215) at higher resolution to read the equations exactly.

Evidence tags: READ (p.N) is read in the paper at printed page N. COMPUTED is a calculation of mine. INFERRED is my
reasoning, not stated by the authors.

Notation. The paper writes the indicator with a calligraphic J and its time average with a barred J. The review digest
(`docs/notes/2026-10-04-digest-musielak-quarles-2014-three-body-problem-review.md`) and much later literature write Y and
<Y>. They are the same quantities: J = Y, J-bar = <Y>. I use J and Jbar below, and write MEGNO for the name.

## 1. What the paper is

A galactic-dynamics paper about a star in a 2D non-axisymmetric potential, not an orbital-mechanics paper. Three parts
(READ, abstract p.205):

1. Sections 2-3 (pp.207-213): an approximate third integral from a pendulum model that separates "loop" and "box"
   orbits, and a representation of the whole energy surface on the sphere S^2.
2. Section 4 (pp.213-219): the new indicator, MEGNO, applied to the 2D logarithmic potential.
3. Section 5 (pp.219-225): global structure of the logarithmic potential: stability of the two principal periodic
   orbits, stochastic layers, their connections.

Only section 4 and the discussion in section 6 (pp.225-226) matter for the project. The abstract's headline claims,
READ: MEGNO "is suitable to investigate the phase space structure associated to a general Hamiltonian", "is effective to
obtain a picture of the global dynamics and, also, to derive good estimates of the largest Lyapunov characteristic
number in realistic physical times", and "the whole chaotic domain appears to be always confined to narrow filaments,
with a Lyapunov time about three characteristic periods" (p.205).

## 2. The definition, exactly as printed

Lyapunov characteristic number (LCN), READ (p.214, eq. 21): sigma(gamma) = lim_{t to infinity} (1/t) ln( delta(gamma(t)) /
delta0 ), where delta(gamma(t)) = |delta(gamma(t))| is the length of an infinitesimal displacement from the orbit gamma
and delta satisfies the variational equations delta-dot = Lambda delta, with Lambda the Jacobian matrix of the vector
field. In integral form (p.214):

    sigma(gamma) = lim_{T to inf} (1/T) integral_0^T [ deltadot(gamma(t)) / delta(gamma(t)) ] dt = < deltadot / delta >,

where "deltadot = delta . deltadot / delta is the time derivative of delta(gamma(t))" (the time derivative of the NORM,
not the vector) "and < > denotes the usual time-average".

MEGNO, READ (p.214, eq. 22): "Let us define the Mean Exponential Growth factor of Nearby Orbits, J, as"

    J(gamma(T)) = (2/T) integral_0^T [ deltadot(gamma(t)) / delta(gamma(t)) ] t dt.

So the logarithmic derivative of the displacement is weighted by time t, and the weighting is what separates regular from
irregular motion.

Regular (quasi-periodic) orbit, READ (p.214, eq. 23): "J(gamma_r(T)) approx 2 (1 - ln(1 + lambda T)/(lambda T)) + O(gamma_r(T))",
where lambda > 0 is "the linear rate of divergence in a neighbourhood of gamma_r" and O(...) is an oscillating term of zero
mean. Derivation as printed: for stable regular motion, delta approx delta0 [1 + lambda t + t u(t)] with u bounded, |u| < b <
lambda. The result oscillates about 2 with amplitude |J - 2| < 4 ln((lambda + b)/(lambda - b)) for large T (p.214).
lambda is "a measure of the lack of isochronicity", precisely the largest eigenvalue of the matrix d omega / d I (omega the
frequencies, I the actions). Smaller lambda means slower convergence to 2, so convergence is slowest near a stable
periodic orbit (p.214).

Irregular orbit, READ (p.215, eq. 24): J_i = J(gamma_i(T)) approx sigma T + O-hat(gamma_i(T)), with sigma the LCN and
O-hat an oscillating term "in general neither quasiperiodic nor periodic".

Averages, READ (p.215, eq. 25; and eq. 26, p.215): J itself has no formal limit, "but J/T and the mean values Jbar_r, Jbar_i
have an asymptotic law for T to infinity":

    Jbar_r approx 2,   Jbar_i approx sigma_i T;   J_r/T approx 2/T,   J_i/T approx sigma_i.

The mean is taken from samples, Jbar(T_k) = (1/k) sum_{i=1..k} J(T_i), T_k = T0 + k DeltaT, with DeltaT approx 0.06
(eq. 26, p.215); the paper says the continuous mean T^-1 integral_0^T J dt "would provide a smoother behaviour ... which is
independent of the time step". A "unique" form follows: Jbar(T) approx a T + b, "with a = 0, b = 2 for regular,
quasiperiodic motion and a = sigma, b = b0 approx 0 for irregular one. If delta grows with some power of t, say n, as it
could happen in some degenerated cases, a = 0, b = 2 n. Only when the phase space has a hyperbolic structure, where
nearby orbits diverge exponentially with time, a is not zero and the MEGNO grows with time. This occurs for irregular,
chaotic motion and also, for instance, for unstable periodic orbits" (p.215; Giorgilli et al. 1997 is cited there for a picture of the hyperbolic structure).

Printed slip, to be treated respectfully: eq. (25) writes Jbar_i approx sigma_i T, but J_i approx sigma T from eq. (24)
makes the mean of J over [0, T] equal to sigma T / 2. The text itself says so two paragraphs later: for the least-squares
fit "we add a factor 2 in the derived slope to compensate the average introduced in J ... since for an irregular orbit J
grows nearly linear, the slope derived from Jbar would be underestimated in a factor 2" (p.215). So the working rule is
LCN approx 2 x (slope of Jbar), equivalently LCN approx J/T. The check by hand (COMPUTED): with constant d(ln delta)/dt =
sigma, J = (2/T) integral sigma t dt = sigma T and its mean is sigma T / 2. I also checked eq. (23) by hand: with
delta = 1 + lambda t, J = (2/T) integral lambda t/(1 + lambda t) dt = 2 (1 - ln(1 + lambda T)/(lambda T)). It matches.

A useful identity that the paper does not print (COMPUTED): integrating (22) by parts gives

    J(T) = 2 [ ln(|delta(T)|/|delta0|) - (1/T) integral_0^T ln(|delta(t)|/|delta0|) dt ],

that is, twice the gap between the final logarithmic growth and its time average. This gives the same two limits at once
(exponential growth ln delta = sigma t gives sigma T; linear growth gives 2 (1 - ln(1+lambda T)/(lambda T))), and shows that
MEGNO needs only the history of ln|delta|.

## 3. How it is computed

READ (p.215): "The computation of J was done using (22) for a given set of initial conditions. All the integrations were
carried out for a realistic time scale, T ~ 10^3 T_D approx 3000 where T_D is the period of the long-axis periodic orbit.
Therefore the computational effort per unit time is almost the same needed to compute the LCN but comparatively shorter
motion times are required. The renormalization of delta (if necessary), proceeds naturally from (22)." The variational
equations are solved with delta0 "along the x axis for loops and along the p_x axis for boxes with |delta0| = 1 and random
sign", with a Runge-Kutta 7/8 (DOPRI8, Prince and Dormand 1981; Hairer, Norsett and Wanner 1987), energy preserved to about
10^-13 (p.215).

The extra differential equations for the indicator are not printed in this paper (READ absence); the review digest cites
equations for them from the authors' other papers. The paper only says to evaluate (22) along the orbit.

For the sphere representation, the authors' Poincare-section variables (eqs 18-20) are specific to their 2D system and
not needed here.

The least-squares estimate of the LCN, READ (p.215): fit the slope of Jbar over "the last 85% of the time interval
(450 <= T <= 3000), just to avoid the initial transient", then double it. They call the result sigma_ls.

## 4. The authors' stated advantages and limits

### Against the Lyapunov characteristic number (pp.206, 215-216, 225-226)

- Time needed. READ (p.206): motion times of order 10^4 periods give only "lower bounds"; "~ 10^5 - 10^6 periods are
  necessary to obtain an accurate determination"; "motion times >~ 10^5 or 10^6 orbits turn out to be too large (in a
  computational sense)". MEGNO reaches useful estimates in "realistic times, ~ 10^3 periods" (p.207).
- The LCN estimate sigma(T) = ln T / T for a regular orbit converges slowly, with a relative error about T_L / T, where
  T_L = 1/sigma is the Lyapunov time. READ (p.216): "T ~ 10^3 T_D is not enough to separate a chaotic region with
  T_L ~ 10^3 T_D from the regular one". The least-squares estimate gives T_L ~ 10^6 T_D for the regular component in 10^3
  T_D, while plain sigma gives T_L ~ 10^2 T_D; "to get such long values of T_L ... by means of the computation of sigma,
  the total motion time should be T ~ 10^7 T_D" (p.216).
- The LCN "does not furnish any information about the structure of the regular component" (p.214); MEGNO does (peaks,
  valleys, resonance widths, section 4.2).
- Hyperbolic features. MEGNO "allows to identify clearly regular and irregular motion as well as stable and unstable
  periodic orbits" (p.226).

### Against the other tools (section 6, pp.225-226)

The authors list (a) Poincare sections: only for 2D, very long integrations for fine structure; (b) the LCN; (c)
spectral analysis (Binney and Spergel; Carpintero and Aguilar): classifies orbits in "relatively short motion times, less
than 10^3 periods" (p.225); (d) the frequency map analysis (FMA, Laskar): "very precise estimates of the
frequencies" for "short-to-moderate motion times (<~ 10^3 periods...)", giving the rotation numbers that label each
resonance, but developed for near-integrable systems; they find that "Only a few of the resonances marked in [Figs 10 and
12] are present in Figs 8 and 9 of PL96" (the FMA study) and that MEGNO shows the size and internal structure of a
resonance; (e) stretching numbers, helicity and twist angles (Contopoulos and Voglis), whose "main advantage ... is its efficiency to
separate regular and stochastic domains in rather short motion times, <~ 50 periods" (p.226); (f) the authors' own conditional
entropy indicator (Cincotta and Simo 1999, 2000), which is more sensitive to periodic orbits and "looks noisy and depends
more strongly on the time step" when many high-order resonances are present (p.226). Their conclusion: "the MEGNO is the
simplest way to obtain such information on the phase space" and they suggest combining it with a precise spectral
analysis for a complete picture (p.226).

The fast Lyapunov indicator (FLI) is not discussed in this paper (READ absence; it appeared in the literature after
the work cited here). The review digest's comparison with FLI comes from the review, not from this paper.

### Limits stated or visible in the paper

- No formal limit of J for a regular orbit: "for arbitrary gamma_r, the formal limit, lim J(gamma_r(T)), does not exist"
  (p.215); only the mean Jbar tends to 2. J has oscillations; for a regular orbit near an unstable periodic orbit it has
  "several local maxima of decreasing amplitude" (Orbit B in Fig. 7, p.215), the amplitude falling about like 1/T.
- Convergence to 2 is slow where the twist is small, that is, near stable periodic orbits (p.214).
- The estimate of the least-squares LCN needs a long fit window (the last 85 percent) and a factor 2 (p.215).
- Evidence base: ensembles of about 400 to 3500 initial conditions along a single line of a 2D autonomous Hamiltonian
  with a smooth, singularity-free potential, integrated to about 10^3 characteristic periods. Nothing is shown for
  non-autonomous systems, for singular or close-encounter dynamics, for three or more degrees of freedom (deferred to a
  future paper, p.226), or for using a single orbit instead of an ensemble ("will be the subject of a future paper",
  p.226).
- The statement of the asymptotic values is physical argument, not a theorem: "some analytical arguments ... are behind
  this method ... but a rigorous theory is still lacking" (said of the conditional entropy, p.226); eq. (23) rests on the
  assumed form delta = delta0 [1 + lambda t + t u(t)] (p.214).

## 5. Integration times needed (all READ)

| Quantity | Time | Page |
| --- | --- | --- |
| Single orbits in Fig. 7 (regular A, B; irregular C) | T = 3000, about 10^3 T_D with T_D about 3 | p.215 |
| Ensembles of 3500 orbits (Fig. 9), 400 orbits (Fig. 11b) | T = 3000 | pp.216-219 |
| Sufficient for the resonance structure | "shorter motion times, about 500 T_D, would be enough" | p.218 |
| Stochastic-layer widths and sigma_ls (Figs 16-17) | T = 10^4 T_x for the LCN, 3.5 x 10^4 T_x on the Poincare section | pp.222-224 |
| Lyapunov time of the main stochastic layer | T_L about 3 T_x, with sigma T_x in 0.2-0.4 | p.224 |

For comparison the standard LCN would need 10^4 to 10^6 periods (p.206). The sampling step for Jbar is DeltaT = 0.06
(p.215).

## 6. Every worked example, in brief

- Pendulum model (section 2, pp.207-212). For a 2D non-axisymmetric potential phi = Phi(m_q), m_q = R^2 + z^2/q^2, the
  Hamiltonian is expanded in the flatness parameter alpha = (1 - q^2)/q^2. Averaging over the unperturbed motion gives an
  approximate invariant K = p_theta^2/2 - (alpha/2) g cos 2 theta and a frequency Omega^2 = alpha g / 2 (eqs 5, 8). Boxes
  are |K| < Omega^2, loops K > Omega^2, the separatrix at K = Omega^2. Valid for q near 1; the authors say the pendulum
  fails as q falls (a banana sub-family appears at q = 0.7). Tested on the logarithmic potential
  phi = (p0^2/2) ln(x^2 + y^2/q^2 + rc^2), p0^2 = 2, rc = 0.1, energy h = -0.4059, q = 0.9, 0.8, 0.7 (Figs 1-4).
- Sphere S^2 representation (section 3, pp.212-214, Figs 5-6): a Hopf-type map of the energy surface section to the
  sphere, so that the periodic orbits along the axes appear at the poles.
- MEGNO, three orbits (section 4.2, Figs 7-8, pp.215-216), q = 0.9: A (x0 about 0.33) stable quasi-periodic, "saturates
  very fast from below to 2 without any significant oscillation"; B (x0 about 0.02) stable but very close to an unstable
  periodic orbit, decaying oscillations about 2; C (x0 about 0.002) in the stochastic layer, J grows nearly linearly.
  Fig. 8: the LCN estimates for A and C against the theoretical values; the printed theoretical values for the regular
  orbit are -3.18 and -2.57 in logarithm, which are log10(2/T) and log10(ln T / T) at T = 3000 (COMPUTED, both reproduce).
  They are therefore the two formulas of eq. (25) and the sigma(T) = ln T / T law evaluated at T = 3000, not independent
  data.
- Ensembles (Fig. 9, p.217): about 3500 orbits per family for q = 0.9, 0.8, 0.7, sigma against sigma_ls. The two agree in
  the gross stochastic layer; sigma_ls shows the regular structure underneath, with sigma about -2.57 over all regular
  regions (log10).
- Resonances (Figs 10-12, pp.218-219): peaks of sigma_ls and of Jbar mark unstable periodic orbits and quasi-periodic
  orbits near them; valleys mark orbits locked inside a resonance; the width of a peak or valley measures the size of
  the resonance. Fig. 11b for q = 0.7 loops: "the range in Jbar is [1.98, 2.035] and the dotted line is the level
  Jbar = 2". Text: "Along this interval Jbar is very close to 2 ... In any case 1.98 < Jbar < 2.035".
- Parameter space (Fig. 13, p.220): the physical region q0(h) <= q <= 1 with q0^2(h) = 1/2 - exp(-2h) in the scaled units
  and q0 about 0.696 at h = -0.4059.
- Stability of the long- and short-axis periodic orbits (section 5.1, Fig. 14, pp.220-222): the normal variational equations
  reduce to a Hill/Mathieu-type equation (eq. 31) and the stability index is the trace of the monodromy matrix; the
  instability tongues start at (h, beta) = (0, n^2) and widen with energy; the width of each tends to 1.
- Stochastic layers and their connections (section 5.2, Figs 15-17, pp.222-225): at h = 15 and 7 the layers around the 1:1,
  4:3, 3:2, 2:1, 5:3 resonances form one channel; the layer width and sigma are computed for q = 0.9, 0.8, 0.71 over h from
  2.1 to 6.1 with T = 10^4 T_x; log sigma = c(q) h + d(q) with c about -0.457, -0.438, -0.421 and d about -0.09 to -0.06
  (eq. 36); sigma T_x about 0.3, bounded by 0.2-0.4, independent of h. These slopes are consistent with base-10 logs
  (-0.434 is log10 e), which is also what the Fig. 8 values need.
- Appendix (pp.226-227): analytic results for the limit rc to 0 of the logarithmic potential (for the collinear problem with force
  -1/x, the natural continuation of a collision orbit crosses the origin in the same direction, unlike the elastic bounce
  of the two-body Kepler regularisation) and the Hermite-polynomial
  solutions of the normal variational equation, giving beta = n for the instability tongue edges.

## 7. Techniques applicable to the project's problems

The project currently has no MEGNO, FLI or Lyapunov-time code in `src` as far as a text search for "MEGNO" and "Lyapunov"
shows (INFERRED from the search; the only Lyapunov-exponent code is a discrete-QR routine for whisker tori in
`search/ccr4bp_whisker.py`, which is a Floquet-style computation on a known torus, not an indicator on a grid). The
review digest's statement that no chaos-indicator code exists is therefore right for MEGNO and FLI but the whisker module
should be known when a positive control for an exponent is needed.

### 7.1 Recipe in eight lines

1. State: the project's flow already integrates a state and a tangent matrix: `core.cr3bp.cr3bp_stm_eom` and
   `core.ccr4bp.ccr4bp_stm_eom` carry the 6-state plus a flattened 6 x 6 matrix Phi. Seed Phi with column 1 equal to a unit
   random vector delta0 and the other columns zero; then delta(t) is column 1 and delta-dot is column 1 of the returned
   dPhi. For `#895` use `search/titania_oberon_realeph_895.ForceModel.accel_grad` (gravity-gradient tensor) to build
   delta' = (delta-v, G delta-r).
2. Wrap the right-hand side: append two scalars y and w. xi = (delta . delta-dot)/(delta . delta) is d ln|delta|/dt;
   y' = t xi; w' = 2y/t (start at small t0 > 0 with y = w = 0).
3. J = 2y/t; Jbar = w/t. These are the paper's J and its continuous mean.
4. Integrate in physical time with an explicit 8th-order method at tight tolerance (the paper used DOPRI8 with energy
   drift about 1e-13). Renormalise delta every chunk and carry y and w: xi does not depend on |delta|.
5. Run to about 10^3 characteristic periods; the paper says about 500 is enough for classification (p.218).
6. Classify by Jbar(T): within about 0.05 of 2 means regular, quasi-periodic (the paper's printed band is 1.98 to 2.035);
   linear growth means hyperbolic or chaotic, with LCN about 2 x slope of Jbar over the last 85 percent of the run, or J/T.
7. Before trusting a map, reproduce the positive control of section 7.4.
8. Use it on a grid of initial conditions, not a single arc, and read peaks as hyperbolic orbits and valleys as
   resonance-locked orbits (p.218); report Jbar and J/T together, and the integration time.

### 7.2 What MEGNO can and cannot say about a cycler or a quasi-cycler

Cannot add anything about the cycler orbit itself. A cycler is hyperbolic, so its tangent flow grows like the largest
Floquet multiplier, J approx sigma T with sigma = ln|multiplier|/period, and Jbar grows at sigma/2 (eqs 24, 25 with the
factor-2 correction above). That repeats the Floquet analysis the project already does. As an internal consistency check
only (not a golden): the `#890` orbit has a largest multiplier of 8.4e5 over a cycle of 123.162 d (`data/OUTSTANDING.md`
`#890`, taking the multiplier to be per cycle, INFERRED), giving sigma = ln(8.4e5)/123.162 d = 0.111 per day (COMPUTED) and
Jbar leaving the paper's regular band (above 2.035) after about 37 days, about a third of a cycle (COMPUTED, neglecting the
initial transient). Hyperbolicity of that strength separates from 2 within a fraction of a cycle.

Can say something about the NEIGHBOURHOOD. The question for a quasi-cycler is whether nearby motion is regular
quasi-periodic or chaotic. A grid of initial conditions about the cycler, or about an arc, with Jbar near 2 after hundreds
of periods marks regular motion; Fig. 7 shows how a hyperbolic periodic orbit appears in the map (B: a regular orbit near an
unstable periodic orbit has decaying peaks; Figs 10-12: peaks at unstable periodic orbits and the stochastic layer, valleys
inside resonance islands). The natural targets in the project are the stability-map sweeps (`#378`, `#681`, `#908`),
where a grid exists and a time-to-escape criterion is used today.

Cannot certify a torus. A partially hyperbolic invariant torus or a whisker also has an exponentially growing tangent
direction, so Jbar would grow linearly on it (INFERRED from eq. 24 and the statement that only the absence of a hyperbolic
direction gives a = 0). It therefore does not answer the `#907` question of whether the `#895` arcs lie on an invariant
torus.

Cannot work on short arcs. The `#895` chain is seven encounters in 370 days; the paper's evidence is for about 10^3 periods,
and it needs about 500 for classification (p.218). The paper's own caveat that a single orbit, as opposed to an ensemble,
is for a future paper (p.226) applies.

Does not by itself apply to non-autonomous models. The values 2 and sigma T were argued for an autonomous Hamiltonian
close to integrable (action-angle form). The `#890` model is time-periodic and the `#895` model is driven by an
ephemeris; the paper does not say what happens (INFERRED: a time-periodic system can be treated by the stroboscopic map,
a quasi-periodically driven one needs the extended phase space; nothing in the paper confirms 2 for these).

### 7.3 Close flybys and regularisation

The paper has no singular or close-encounter case (READ absence). The identity in section 2 gives a direct view
(COMPUTED): if the growth factor of |delta| across encounter k is g_k, occurring at time t_k, then y = integral t d ln|delta|
collects t_k ln g_k from each encounter, so

    J(T) approx (2/T) sum_k t_k ln g_k + smooth part.

Each flyby therefore adds a lump to J, weighted by its time, and late encounters count more. For a cycler that repeats the
same flyby every Delta, J grows like T ln g / Delta, a linear growth with sigma = ln g / Delta, as it should. It also
means a single close encounter early in a run changes J by a fixed amount that is later diluted by the 1/T factor, while
one late in the run dominates.

Whether regularised variables are needed is a numerical question and not a MEGNO one. If the flyby is well resolved in
physical time, as the `#890` and `#895` flybys at about 1000 km altitude are, the integral is smooth enough. Use physical
time for the weight t. If an integrator is run in a regularised (Sundman) time s, then y' = t(s) d ln|delta|/ds is still
correct, but the paper's J and Jbar are defined in physical time and a quantity formed with s as the weight is a different
number (INFERRED). For genuine near-collision passes (`core/cr3bp_regularized.py`) the tangent equations also need
regularising consistently.

The norm. The paper uses the Euclidean norm of delta in its own phase-space variables. Mixing position and velocity
components gives the norm a units choice; for a non-dimensional CR3BP state one would use the non-dimensional (x, v).
The asymptotic values depend on the growth law, not on the choice of equivalent norm (INFERRED), but the numbers at finite T
and the oscillating term do.

### 7.4 Positive control the project could reproduce

The control is the 2D logarithmic potential, which is a few lines of code with an analytic Hessian and no singularity.

Model (READ, eqs 15-16, pp.210, 212): phi(x, y) = (p0^2/2) ln(x^2 + y^2/q^2 + rc^2), H = (px^2 + py^2)/2 + phi, with
p0^2 = 2, rc = 0.1 and energy h = -0.4059 (the value of Papaphilippou and Laskar 1996, p.210). Initial conditions: y0 = 0,
px0 = 0, x0 given, py0 from the energy (py0 = sqrt(2 (h - phi(x0, 0)))); delta0 along the x axis in the paper, any generic
direction in my runs.

Printed values to compare against (the expected side must be these, not my runs):

1. q = 0.7, h = -0.4059, loop orbits with x0 in 0.04 to 0.076, y0 = 0, px0 = 0, T = 3000 (about 10^3 T_D): Jbar lies in
   the range [1.98, 2.035] (Fig. 11b and text, p.219), apart from the narrow resonance structure the figure shows.
2. q = 0.9, x0 about 0.33 (Orbit A): Jbar "saturates very fast from below to 2 without any significant oscillation" (p.215).
3. q = 0.9, x0 about 0.002 (Orbit C, stochastic layer): J grows nearly linearly with T and reaches roughly 280 at T = 3000
   as read off Fig. 7 (an approximate figure reading, not a printed number); the LCN, Fig. 8, is about 10^-1 (log10
   about -1).
4. q = 0.9, x0 about 0.02 (Orbit B, near an unstable periodic orbit): Jbar close to 2 with decaying oscillations, the end
   value read from Fig. 7 about 2.03 (approximate).
5. Two further printed numbers that can be recomputed from the paper's own formulas, which checks the parameters but not the
   dynamics: P = (2h - p0^2 ln rc^2)^(1/2) approx 2.9 (p.210) gives 2.898 (COMPUTED), and q0 approx 0.696 for
   q0^2 = 1/2 - rc^2/exp(2h/p0^2) gives 0.6965 (COMPUTED).
6. The -3.18 and -2.57 values of Fig. 8 are log10(2/T) and log10(ln T / T) at T = 3000 and test nothing about the code.

A one-off check of the recipe (COMPUTED; throwaway script in the session scratch area, not committed). Setup: DOP853 with
rtol = atol = 1e-12, a random unit 4-vector delta0 (not the paper's delta0), renormalised every 25 time units, y and w
accumulated as in section 7.1, py0 > 0, T = 3000, Jbar the continuous mean. Energy drift at the end of every run was at most
1e-9.

| q | x0 | J(3000) | Jbar(3000) | J/T | Paper |
| --- | --- | --- | --- | --- | --- |
| 0.9 | 0.33 | 2.27 | 1.989 | 7.6e-4 | A: saturates to 2 (agrees) |
| 0.9 | 0.02 | 2.83 | 2.078 | 9.4e-4 | B: about 2.03 read from the figure; mine is larger, and x0 about 0.02 lies next to an unstable periodic orbit, so I do not claim agreement |
| 0.9 | 0.002 | 278 | 166 | 0.093 | C: J about 280 (agrees), Jbar about 140 read from the figure (same order only, since delta0 and the exact x0 differ) |
| 0.7 | 0.045 | 2.21 | 1.995 | 7.4e-4 | inside [1.98, 2.035] |
| 0.7 | 0.065 | 2.60 | 1.996 | 8.7e-4 | inside [1.98, 2.035] |
| 0.7 | 0.09 | 2.41 | 2.000 | 8.0e-4 | outside the printed window, regular as expected |

Orbit C gives J/T = 0.093, log10 about -1.03, in line with the printed LCN of about 10^-1 (Fig. 8). The first two printed
checks (A, and the band at q = 0.7) are matched by the recipe, which is the evidence that the accumulators and the variational
equations are right. This is a check of the recipe against the paper's behaviour, not a golden value; the goldens are the
printed numbers above. The run took 13 to 34 s per orbit in plain Python.

A second control that needs no dynamics at all, for testing the accumulator code alone (independent closed forms, COMPUTED):
feed delta(t) = 1 + lambda t and expect J = 2 (1 - ln(1 + lambda T)/(lambda T)) (eq. 23); feed delta(t) = exp(sigma t) and
expect J = sigma T and Jbar = sigma T / 2.

A third control with a known answer inside the project: a known hyperbolic periodic orbit (the `#890` orbit, or any CR3BP
orbit with a Floquet multiplier from the project's monodromy code) must give a slope of J of ln|multiplier| / period. This
is an internal consistency check, not a golden.

## 8. Summary for the coordinator

- Definition: J(T) = (2/T) integral_0^T [d|delta|/dt / |delta|] t dt; regular quasi-periodic motion gives Jbar to 2 (and J to
  2 (1 - ln(1+lambda T)/(lambda T)) with oscillation); chaotic motion and unstable periodic orbits give J approx sigma T
  and Jbar approx sigma T / 2 (the paper prints sigma T in eq. 25, a slip its own factor-2 remark corrects).
- It needs only the tangent flow the project already integrates plus two scalars.
- Advantage stated: about 10^3 periods versus 10^5 to 10^6 for the LCN; shows the structure of the regular component.
- It cannot add to the Floquet analysis of a hyperbolic cycler, cannot certify a torus, and is unproven for short arcs,
  non-autonomous models and close encounters; it suits grid maps.
- Positive control: 2D logarithmic potential, printed band 1.98 < Jbar < 2.035 for q = 0.7 loops at T = 3000, and Orbit A
  saturating to 2; reproduced by the recipe in a one-off run.
- No FLI treatment in the paper; the paper's own stated omission is the single-orbit and 3D use (deferred).
