# Digest: Prussing 2000, "A Class of Optimal Two-Impulse Rendezvous Using Multiple-Revolution Lambert Solutions" (#960 batch 30)

J. E. Prussing (University of Illinois at Urbana-Champaign), The Journal of the Astronautical Sciences 48(2-3):131-148
(April-September 2000), doi 10.1007/BF03546273. Presented as AAS 00-250 at the Battin Astrodynamics Symposium, March 2000.
- Given file: supplied file `06fa03fb-trussing2000.pdf`, 18 pp., text layer, md5
  22259627d16ad79813c9b2ee88dd9031. (The upload name has a typo, "trussing".)
- Filed as `cyclers_pdf/papers/prussing-2000-class-optimal-two-impulse-rendezvous-multiple-revolution-lambert-solutions-jas-48-131-doi-10.1007-BF03546273.pdf`.
- How I read it: the full text; page images of pp. 132-134 (Eqs. 2-8, Fig. 1-2), 141 (Table 1), 144 (Table 2), 145 (Tables 3-4).
  Equations 9-14 and the Appendix are from the text layer (garbled; used for the logic only). Every table cell was
  compared with my own recomputation (section 3), which also checks the digits I read.
- Wanted-list row 46 ("Prussing 2000 (10.1007/BF03546273)" in the methods row). Also cited as [51] in
  `2026-06-10-saloglu-2023-iso-impulse-mining.md` and as [14] in the Russell-Ocampo 2005 digest.

## 0. Verdict

**A clean statement of the multi-revolution Lambert solution count and of how to find all solutions with the classical
Lagrange form. The formulas and all four printed tables reproduce, with one cell 0.8% off and one printed-formula
convention that needs a sign fix. The project's own multi-rev Lambert agrees with it.**
- Count: if the transfer time allows N_MAX complete revolutions, there are exactly 2 N_MAX + 1 solutions: one with N = 0
  and two (low and high energy) for each N = 1..N_MAX; the two N_MAX solutions merge at t = t_min,N_MAX.
- Catalogue: no row depends on it. **Proposal only:** none. For `src/cyclerfinder/core/lambert.py`, cite it as an independent
  statement of the 2N+1 rule (the module cites only Vallado and Bate-Mueller-White; grep of `src` for "prussing" finds only the
  2010 primer-vector chapter in `verify/primer.py` and `verify/primer_refine.py`).
- Use for the project: a ready oracle for multi-rev solution counts and for the same-circle rendezvous problem (Tables 1-4).
  It is also the clean statement of "waiting is not optimal, final coast is": an optimal two-impulse transfer may use a
  shorter time and then coast with the target (section 4).

## 1. Lagrange multi-rev formulation (READ pp. 132-134, image)

- Eq. 1 Lambert's theorem: t_f = F(a, r1 + r2, c). Eq. 2: sqrt(mu) t_f = a^(3/2) [2 N pi + alpha - beta - (sin alpha - sin beta)].
- Eq. 3: sin(alpha/2) = (s/(2a))^(1/2), sin(beta/2) = ((s - c)/(2a))^(1/2), with s = (r1 + r2 + c)/2 and chord c.
  Minimum-energy ellipse a_m = s/2, where alpha = pi.
- Eq. 4: sqrt(mu) t_mN = (s/2)^(3/2) [(2N + 1) pi - beta_m + sin beta_m], Eq. 5: sin(beta_m/2) = ((s - c)/s)^(1/2).
  t_mN splits each N-curve into a lower portion (alpha = alpha_0, the principal arcsine) and an upper portion (alpha = 2 pi - alpha_0).
- Eq. 6, parabolic time: sqrt(mu) t_p = (sqrt 2 / 3)[s^(3/2) - sgn(sin theta)(s - c)^(3/2)]. N = 0 has no minimum time (it tends to t_p).
- For N > 0 there is a minimum time t_min,N at a > a_m. Eq. 7: dt_f/da = (1/2)(a/mu)^(1/2) f(a) / (sin(alpha - beta) + (sin alpha - sin beta)),
  Eq. 8: f(a) = [6 N pi + 3(alpha - beta) - (sin alpha - sin beta)] [sin(alpha - beta) + (sin alpha - sin beta)] - 8[1 - cos(alpha - beta)].
  t_min,N is at f(a) = 0 (Newton from 1.001 a_m). N_MAX = largest N with t_min,N <= T. Then solve g(a) = t_f(a) - T = 0
  by Newton for the 2 N_MAX + 1 roots, with Eq. 7 as the derivative. Eccentricity from Eq. 10 (text layer only).
- **Printed convention problem.** The paper says beta = beta_0 for 0 <= theta < pi and beta = 2 pi - beta_0 for pi <= theta < 2 pi.
  With that, N = 0 gives negative flight times for theta > pi (computed: theta = 252 deg, a = 1, both branches: -0.687 and -0.300;
  the standard beta = -beta_0 gives 0.313 and 0.700). The two forms differ by 2 pi in (alpha - beta), so the printed one shifts
  N by one. All tables below reproduce with beta = -beta_0 and none (47 rows differ) with the printed form. I suspect a
  typesetting slip, not a computing error. Use beta = -beta_0 for theta > pi.
- The paper says the formulation is "not the most efficient computationally"; it names the Battin-Vaughan algorithm as better.

## 2. Same-circle rendezvous (READ pp. 137-143)

- Two-impulse, time-fixed rendezvous of a spacecraft with a target in the same circular orbit. Canonical units, mu = 4 pi^2,
  circular period 1. Target starts "initial angle" ahead, so theta = initial angle + 360 T (mod 360).
- An initial coast never helps; a final coast can (the transfer is done early, then both coast).
- Zero transfer angle (chord 0): the optimal transfer is tangential, like Hohmann. The two impulses are equal and opposite.
  Then dv = 2|v - v_c| with transfer period T/N. Minimum-energy time t_mN = sqrt(2) N / 4 in circular periods:
  0, 0.35355, 0.70711 for N = 0, 1, 2 (printed; reproduced).
- Primer vector (Eq. 12 and conditions): the primer magnitude must stay <= 1 and reach 1 at interior impulses.
  Fig. 6 (T = 0.9) is an optimal two-impulse primer history. Fig. 7 is a non-optimal one: the excursion above 1 says
  add a midcourse impulse. A four-impulse trajectory at T = 2.3 costs 1.19 against 1.33 for the best two-impulse (text).
  The linear theory says at most four impulses; the nonlinear results here are "more irregular".
- Tables 1-4 (initial angles 180, 36, 144, 288 deg): columns T, theta, dv*, N*, N_MAX, L/H (lower or higher energy, smaller or larger a).
  dv* is the total of both impulses in units where the circular speed is 2 pi (check: T = 0.5, theta = 0, N = 1: a = 0.5^(2/3),
  dv = 2(2 pi - 4.036) = 4.49, printed 4.50).
- Observations in the text: local minima of dv* at the zero-angle times; each is lower than the one before, so no optimal T
  exists. At T = 3.6 (Table 3), N_MAX = 10 (21 trajectories) but N* = 4. The text calls Table 4 "initial angle of 228" but its
  caption and values use 288; the caption is right (image, and theta = 288 + 360 T reproduces the theta column).

## 3. Checks (COMPUTED)

Scripts `check_prussing_tables.py`, `check_prussing_eqs.py`, `check_prussing_project_lambert.py` and their `.out` files.
- Independent Lagrange solver (Eqs. 2-3, beta = -beta_0 for theta > pi), roots by bracketing, t_min,N by minimisation,
  dv from the Lagrange f and g coefficients and circular-velocity matching at both ends, for all 67 table rows.
  Results: theta, N*, N_MAX and dv* all match in 66 rows (dv within 0.6% of the print; most to 3-4 digits).
  - One cell differs: Table 3, T = 3.5: printed 0.552 (image), I get 0.557 (0.8%). N* = 3 and N_MAX = 6 agree. Unresolved.
  - At theta = 180 deg (Table 2, T = 0.4) the velocity direction is undefined (g = 0); I used theta = 180 deg - 1e-6 and
    got 2.427 against 2.43.
  - L/H at theta = 0: the paper prints H on all 13 zero-angle rows. My label (a < 1 = L) agrees on 3 and not on 10 (those with
    a = (T/N)^(2/3) < 1). The label is a convention at c = 0, where alpha_0 = beta_0 is degenerate; N*, N_MAX and dv* agree on all 13.
- Eqs. 4-8 on the Fig. 2 geometry (r1 = 1, r2 = 1.524, theta = 75 deg): t_p = 0.1976 (printed 0.197); t_mN = 0.49613,
  1.53985, 2.58358 for N = 0, 1, 2 (agrees with direct evaluation of Eq. 2 at a_m); f(a) = 0 gives a = 1.05130, t_min,1 = 1.49148
  and a = 1.03657, t_min,2 = 2.55510, equal to the direct minimiser and to Fig. 2 read by eye. Newton from 1.001 a_m converges.
- Project Lambert: `core/lambert.py` is Vallado Algorithm 5.2 (universal variable z) with multi-rev branches from a per-revolution minimum
  time (`_min_time_of_revolution`) and a derivative-free solve on each side. For r1 = r2 = 1, T = 1.4, theta = 324 deg, mu = 4 pi^2
  (Table 1 row: N_MAX = 2) it returns 5 = 2 N_MAX + 1 solutions. Their v1 match my Lagrange solutions to 1.6e-12 (canonical units).
  Its labels are "low" = larger a (the paper's H) and "high" = smaller a (the paper's L): opposite to the paper's naming.
  It raises `LambertGeometryError` at theta = 0 and 180 deg, so the paper's zero-angle optimum (the tangential case) is outside its range.
  No repository file was changed (checked with `git status`; numba cache sent to the scratch folder).

## 4. Relevance to the project

- Multi-rev legs matter for cycler arcs (multiarc, resonant tours). The 2N+1 rule and the t_min,N criterion are what the
  project's `max_revs` loop implements. Prussing gives the criterion in (a, alpha, beta) form; the project uses z. The same.

## 5. Citation-mining

- Not held: Battin 1987 (book; also Battin 1959/1999 on row 46); Lancaster, Blanchard and Devaney 1966 (JSR 3:1436);
  Sun, Vinh and Chern 1987 (JAS 35:213); Battin and Vaughan 1984 (JGCD 7:662); Gooding 1990 (Celest. Mech. 48:145);
  Loechler 1988 (MIT thesis); Ochoa 1991 thesis and Ochoa-Prussing AAS 92-194; Prussing 1979 (JGCD 2:442); Prussing and Conway
  1993 (book); Stern 1964 (MIT thesis, source of the t_min analysis); Prussing and Chiu 1986 (JGCD 9:17); Lawden 1963;
  Lion and Handelsman 1968; Prussing 1969, 1970, 1995. Only Battin is on the wanted list (row 46).
  New candidates that matter: Gooding 1990 and Battin-Vaughan 1984 (the efficient multi-rev solvers), Ochoa-Prussing 1992
  (the original of this formulation), Lancaster-Blanchard-Devaney 1966.
- HELD: Hollister and Prussing 1965 (`hollister-prussing-1965-optimum-transfer-mars-via-venus-aiaa-65-700...`), not cited here.
  Prussing 2010 (primer vector chapter) is cited in `src` but not held under that name.

*Check scripts and outputs named above are filed beside the PDF as `cyclers_pdf/papers/<pdf stem>-<script name>`.*

*Wanted-list row numbers in this digest are the batch-29 numbering; the list was renumbered in batch 30.*
