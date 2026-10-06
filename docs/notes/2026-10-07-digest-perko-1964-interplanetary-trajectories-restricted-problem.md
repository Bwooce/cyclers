# Digest: Perko 1964, "Interplanetary Trajectories in the Restricted Three-Body Problem" (AIAA Journal) (#960)

L. M. Perko (Lockheed Missiles and Space Company, Palo Alto), "Interplanetary Trajectories in the Restricted Three-Body
Problem", AIAA Journal 2(12):2187-2192, December 1964. doi 10.2514/3.2761 (printed on every page margin).
"Presented as Preprint 64-52 at the AIAA Aerospace Sciences Meeting, New York, January 20-22, 1964; revision received
June 3, 1964" (p.2187). Perko thanks J. V. Breakwell for guidance. The matching frame is "due to J. V. Breakwell" (p.2190 footnote).
- File given: `e0fc80d5-perko1964.pdf`, 6 pages, text layer (AIAA reprint, downloaded 2013). md5 47570b51d4736683b4202765ae5498d4.
- **Proposed corpus filename:** `cyclers_pdf/papers/perko-1964-interplanetary-trajectories-restricted-three-body-problem-aiaa-j-2-12-2187-doi-10.2514-3.2761.pdf`
- **Wanted list:** row 28 holds two items: Perko's 1964 Stanford PhD thesis ("Asymptotic matching in the restricted three-body
  problem") and this AIAA J. paper. This file is the second item. **The row stays open for the thesis.**
- **How I read it.** The text layer loses most of the equations. So I viewed all 6 pages on 110-dpi renders (pp.2187-2192). I read Fig. 3 on a 300-dpi crop.
  I read eqs. (6), (9), (10), (15) and the matched constants on the page images (pp.2188, 2190-2192).
  Checks: `checks_perko1964.py`, output in `checks_perko1964.out`.

## 0. Verdict

**The first journal statement of the Breakwell-Perko matched-asymptotic patched conic, for interplanetary (O(1) angular-momentum)
arrivals.** It is planar, the secondary is on a circular orbit, and only the arrival side is treated.
- The outer solution is a perturbed heliocentric conic with true anomaly theta as the independent variable. The inner solution is a hyperbola about m1.
  They are matched where r2 = O(mu^(1/2)). The matching gives V_inf, the asymptote offset, the pericentre time and the
  axis rotation of the planet hyperbola from the heliocentric initial conditions plus definite integrals.
- **Accuracy claim (p.2190, p.2192):**
  - Each expansion is correct to **O(mu^(3/2))** where theta1 - theta = O(mu^(1/2)).
  - The dominant second-order terms are mu^2 ln(theta1 - theta)/(theta1 - theta).
  - Validity needs **h0 >> mu^(1/2) and |1 - h0| >> mu^(1/2)** (p.2190, image). h0 is the initial angular momentum about m0.
  - **No composite solution:** "Since the inner and outer solutions are written in terms of different independent variables that are not
    related in any elementary fashion, a composite solution cannot be obtained." The uniformly valid approximation is eqs. (7), (8) and (13),
    each used in its own region.
  - **There is no comparison with numerical integration.** The O(mu^(3/2)) claim is asserted ("It can be shown"), not demonstrated.
- **Numerical example: Fig. 3 only** (Mars' effect on Earth-Mars trajectories; sec. 2). No table.
- **What it gives the project:**
  - the scope conditions (h0 and |1 - h0| both >> mu^(1/2)) for any `#948`/`#944` use of first-order matching;
  - the closed-form matched constants for the arrival hyperbola (sec. 1).
- **Catalogue implication (PROPOSAL only):** none. Method paper.

## 1. Method (READ on the images)

- **Equations (1)-(4), p.2187-2188.** Units: distance d (m0-m1), time (d^3/G m0)^(1/2). In polar form about m0, the equations of
  motion are rewritten with theta as the independent variable. The dependent variables are t and u = 1/r. Expansion: u = u0 + mu u1 + ..., t = t0 + mu t1 + ....
- **Zero order:** a conic with h0, e0, a0, omega0. The initial conditions are chosen so that the unperturbed conic hits m1 at theta = theta1
  (r0(theta1) = 1, t0(theta1) = theta1).
- **Singular expansion about theta1, eq. (6), p.2188:**
  - a = e0 sin(theta1 - omega0) / [h0 (1 - h0)];
  - b = -e0 cos(theta1 - omega0) / [2 (1 - h0)^2];
  - K0 = h0^2 (1 - h0) / [(1 + a^2)^(3/2) |1 - h0|^3].
- **First order, eqs. (7)-(8), p.2189-2190.**
  - The 1/(theta1 - theta) terms cancel.
  - The singular part of u1 is then carried by the cosine integral Ci(theta1 - theta), which behaves like ln(theta1 - theta).
  - The time t1 has a (K0/h0^3) ln(theta1 - theta) singularity.
- **Outer expansion near m1, eqs. (9)-(11).** It uses phi = (theta1 - theta)/mu^(1/2), in an m1-centred non-rotating frame turned by alpha = tan^-1(1/a).
- **Inner solution, eqs. (12)-(15), p.2191.** The m1-centred hyperbola is written in the variable F. Its far-field expansion is in z = mu e^|F|.
- **Matched constants (pp.2191-2192, image):**
  - V_inf = |1 - h0| (1 + a^2)^(1/2);
  - h2 = (1 - h0) G7(theta1);
  - t_p = t0(theta1) + (mu/V_inf^3) [1 + ln(mu h0 e2 / (2 V_inf^2))] + mu G6(theta1);
  - e2 = (1 + V_inf^2 h2^2)^(1/2).
  - G6 and G7 are definite integrals of bounded functions from theta0 to theta1.
  - The asymptote offset is Delta = mu h2^2/(e2^2 - 1)^(1/2) = mu |G7(theta1)| / (1 + a^2)^(1/2).
- **Compared with Lagerstrom & Kevorkian 1963** (abstract, p.2187): they treated initial angular momentum about m0 of O(mu^(1/2))
  (some Earth-Moon trajectories), in one rotating frame with distance along the m0-m1 line as the independent variable. This paper
  treats O(1) angular momentum (interplanetary arrivals), with different variables, and says "the basic ideas and results are, however, the same".
- **Compared with Breakwell & Perko 1965** (HELD; digest `2026-10-06-digest-breakwell-perko-1965-...`):
  - The 1965 paper is general: three dimensions, eccentric and inclined planets, transition matrices, departure and arrival, and a
    flyby map. It claims O(lambda^2) through a composite form.
  - This 1964 paper is the planar, circular, arrival-only special case, with theta as the independent variable. It claims O(mu^(3/2)), and it says no composite can be formed.
  - Both give the same structure: a gross bias in t_p plus a local (mu/V_inf^3) ln(...) term, and an asymptote offset.

## 2. Numerical example: Fig. 3 (p.2192; 300-dpi crop; graph readings, approximate)

- **Case:** "a class of Earth-Mars trajectories leaving a massless Earth at perihelion, with initial conditions such that the
  unperturbed heliocentric conic intersects Mars". I take "at perihelion" to mean that the transfer conic's perihelion is at Earth (INFERRED; check (2) below supports it).
- **Abscissa:** transfer angle theta1, 60-290 deg.
- **Delta (left axis, miles):**
  - near 0 below about 140 deg;
  - peak about +2,450 mi near theta1 = 185-190 deg;
  - zero near 214 deg;
  - minimum about -2,800 mi near 248 deg;
  - back to about +1,000 mi near 262 deg, then a dashed continuation that falls toward about +300 mi at 290 deg.
- **Delta t_p (right axis, days, plotted downward):** the peak is about 0.19-0.20 d near theta1 = 180 deg. It is near 0 below 100 deg and above about 260 deg.
- **Text:** Mars makes the particle "arrive earlier". The deflection is "away from the sun for transfers less than approximately 214 deg
  and toward the sun for larger transfer angles". There is a "small interval about theta1 = pi/2" with no points, because |1 - h0| >> mu^(1/2) fails there. The curves are faired across it.
- **Mismatch:** Fig. 3 shows Delta positive again beyond about 258 deg, but the text says the deflection is toward the Sun (negative) for all
  angles above 214 deg. The printed Delta formula uses |G7|, so it gives only a magnitude, yet Fig. 3 plots a signed Delta. The sign convention is not stated. I did not resolve this. A possible cause (INFERRED): the positive bump near 262 deg is within
  about 10 deg of the second h0 = 1 point at 270 deg (sec. 3), where |1 - h0| < 0.02 and the method is near its limit.
- **Compared with Breakwell-Perko 1965 Fig. 3** (from its digest: Delta_2 peak about +2.5 x 10^3 mi; t_p bias about -0.2 d near 180 days).
  The magnitudes agree with this Fig. 3, so it is probably the same computation, plotted there against flight time (INFERRED).

## 3. Checks (`checks_perko1964.py` / `.out`)

0. **Delta formula.** With e2^2 - 1 = V_inf^2 h2^2, h2 = (1 - h0) G7 and V_inf = |1 - h0| (1 + a^2)^(1/2), the printed
   Delta = mu h2^2/(e2^2 - 1)^(1/2) reduces to mu |G7|/(1 + a^2)^(1/2), the printed right-hand side. The script checks this
   numerically. So the transcriptions of h2, V_inf and Delta agree with each other.
1. **V_inf.** The matched V_inf = |1 - h0| (1 + a^2)^(1/2), with the printed a, equals the zero-order patched-conic speed of the
   conic relative to m1 on a unit circular orbit. Three random cases agree to 1e-6. This is a structural check only: at zero order the matching recovers the patched conic.
2. **The gap at pi/2.** For a conic with perihelion at Earth (r_p = 1/1.5237) that reaches r = 1, h0^2 = (1 + e0) r_p and cos(theta1) = (h0^2 - 1)/e0.
   So h0 = 1 (where |1 - h0| >> mu^(1/2) fails) falls exactly at **theta1 = 90 deg**. This matches the text's "small interval about pi/2"
   and supports the perihelion-departure reading.
   - **The same geometry gives h0 = 1 again at theta1 = 270 deg.** The dashed part of Fig. 3 (about 260-290 deg) probably marks
     this second failure region (INFERRED; the text mentions only pi/2).
   - In both regions |1 - h0| falls below about 0.02 within +-10 deg (mu^(1/2) = 0.00057 for Mars).
3. **Size of Delta t_p.** The local term (mu/V_inf^3)[1 + ln(mu h0 e2 / 2 V_inf^2)] alone gives -0.29, -0.19 and -0.09 d for V_inf = 2.6, 3.0 and 4.0 km/s
   (h0 e2 ~ 1 assumed; the G6 term is omitted). This is order-of-magnitude only. The sign (early arrival) and the size (about 0.2 d) agree with Fig. 3.

## 4. Citation mining (5 references, p.2192)

Checked with `ls cyclers_pdf/papers | grep -i` and `grep -i docs/notes/CORPUS_INDEX.md`.

| ref | work | status |
|---|---|---|
| 1 | Lagerstrom & Kevorkian (1963), "Matched conic approximations to the two fixed force center problems", Astron. J. 68 (March 1963) | not held (no "lagerstrom" hit). Not on the wanted list: row 64 names only their J. Mécanique 2:189 paper and Lancaster 1968. **Candidate:** add it to row 64 (low priority; the method origin). |
| 2 | Lagerstrom & Kevorkian (1963), "Earth-to-moon trajectories in the restricted three body problems", J. Mécanique 2 (June 1963) | not held. **Wanted row 64.** |
| 3 | Lagerstrom & Kevorkian (1963), "Some numerical aspects of earth-to-moon trajectories in the restricted three body problem", AIAA Preprint "63-3891" as printed (August 1963; read on a 400-dpi crop, no footnote mark; I did not verify the number) | not held; not on the wanted list. By its title, the only cited item about numerical accuracy (not read). **Candidate:** add it to row 64 (low priority). |
| 4 | Whittaker (printed "K. T."; it is E. T.), A Treatise on the Analytical Dynamics of Particles and Rigid Bodies (Dover 1944) | not held; textbook |
| 5 | Moulton, An Introduction to Celestial Mechanics (1914) | not held; textbook |

- Held follow-ons (not cited here, for orientation): Breakwell & Perko 1965 (`breakwell-perko-1965-...`), Perko 1967 (`perko-1967-...`; the error estimates),
  Breakwell & Perko 1974 (`breakwell-perko-1974-...`), Perko 1974 (`perko-1974-...`), Guillaume 1975 (`guillaume-1975-...`), Kevorkian & Lancaster 1968 (`kevorkian-lancaster-1968-...`).
- **Proposal for row 28:** mark the AIAA J. paper received (Preprint 64-52; journal 2(12):2187). Keep the row open for the Stanford thesis.

*Filed as `cyclers_pdf/papers/perko-1964-interplanetary-trajectories-restricted-three-body-problem-aiaa-j-2-12-2187-doi-10.2514-3.2761.pdf`. Check scripts, outputs and other files named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`.*

*Wanted-list row numbers in this digest are the batch-30 numbering; the list was renumbered after batch 34.*
