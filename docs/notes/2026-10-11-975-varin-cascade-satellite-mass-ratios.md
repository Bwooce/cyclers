# #975: the Varin closed-family cascade i_5-i_8 at the satellite mass ratios: pre-registration (2026-10-11)

Owner of this note: earthmoon-opus. Dispatch: lead, 2026-10-11 (from `#971` R14). Written before any
code or run. No catalogue edits. Nothing found here is called novel; the literature step comes after.

Sources read for this note:
- OUTSTANDING `#975`;
- `docs/notes/2026-10-07-971-fable-corpus-review-2.md` sec. 3 R14;
- the Varin-Bruno digests: `2026-10-07-digest-varin-bruno-2009-closed-families-ssr.md` (the SSR form of Varin's
  KIAM 16/2008; tables identical in all 886 cells) and `2026-10-06-digest-bruno-varin-kiam-preprints-2005-2010.md`
  (KIAM 16/2008 translation; corrected by the SSR digest on two points), plus the Bruno-Varin family digests;
- `data/sources/varin-bruno-2009-closed-families-tables.yaml` (Tables 1, 2, 3.1-3.7, printed values);
- the Voyatzis-Kotoulas digests with their R14 verdicts:
  - `2026-10-07-digest-voyatzis-kotoulas-2005-exterior-resonances-neptune.md`: only an EXTERIOR closed p = -7
    family at the Sun-Neptune mu 5.178e-5, figure only, no numbers. R14 is not pre-empted, but the wording must
    be limited to the interior side;
  - `2026-10-11-digest-kotoulas-voyatzis-2004-...`: no periodic-orbit table, no interior family, no closed family.
    No collision.

## 1. The object

- **Frame and family i.** Varin's frame is the planar CR3BP, synodic, with the origin at P1 (mass 1 - mu), P2 at
  (1, 0) and canonical momenta y. Family i is the generating family of symmetric P1-circling orbits.
- **Closed families i_k.** As mu varies, the characteristics of the symmetric family are level lines of a
  surface over (crossing, mu).
  - A closed family i_k branches off family i at a saddle point mu'_k and shrinks to one orbit at an
    elliptic point mu''_k.
  - It exists as a closed loop (an oval) of symmetric periodic orbits for mu'_k < mu < mu''_k.
  - Both ends obey m_41 = m_14 = 0 of the monodromy matrix, with the symmetry conditions x2(T/2) = y1(T/2) = 0,
    and Tr = trace(M) - 2 = 2 there.
- **Printed data (Table 1, 8 digits).** The saddle and elliptic orbits for k = 1-4: (a~, e~, T~ = T/2pi,
  C, Tr~_v) and the mu values (Table 2).

| k | mu'_k (saddle) | mu''_k (elliptic) | T~ at mu'_k | T~ at mu''_k |
|---|---|---|---|---|
| 1 | 4.13129887e-3 | 3.66863029e-2 | 1.23370552 | 1.97012501 |
| 2 | 6.61705554e-4 | 5.27272358e-3 | 2.26009772 | 2.97373803 |
| 3 | 2.15292269e-4 | 1.88241384e-3 | 3.26879728 | 3.97406324 |
| 4 | 9.54304953e-5 | 8.86552296e-4 | 4.27305773 | 4.97446586 |

- **Varin's empirical law (Table 2; "empirical asymptotics", not a theorem):** mu'_k ~ k^(-8/3) mu'_1 and
  mu''_k ~ k^(-8/3) mu''_1. Extrapolated (my arithmetic; predictions, not data):

| k | mu'_k (law) | mu''_k (law) | expected T~ at mu'_k / mu''_k (trend k + 0.27 / k + 0.97) |
|---|---|---|---|
| 5 | 5.65e-5 | 5.02e-4 | about 5.28 / 5.97 |
| 6 | 3.47e-5 | 3.08e-4 | about 6.28 / 6.97 |
| 7 | 2.31e-5 | 2.05e-4 | about 7.28 / 7.97 |
| 8 | 1.61e-5 | 1.43e-4 | about 8.28 / 8.97 |

- **Families to continue:** i_5, i_6, i_7, i_8: their saddle (birth) orbits, their elliptic (death) orbits, and
  each closed loop at every target mass ratio that lies inside its predicted (mu'_k, mu''_k).
- **Seeds:**
  - (a~, e~) for k = 5-8 come from the k = 1-4 trend of Table 1. At mu'_k, a~ goes 0.689, 0.803, 0.855,
    0.886 and e~ goes 0.728, 0.838, 0.885, 0.911. At mu''_k, a~ goes 0.814, 0.887, 0.913, 0.929 and e~ goes
    1.439, 1.245, 1.176, 1.138.
  - T~ follows the trend above, and mu follows the law.

## 2. Conventions (fixed now; they are the controls' first test)

- **State conversion (Varin to the project's barycentric rotating frame, x = x1 - mu, y = x2,
  u = dx1/dt = y1 + x2, v = dx2/dt = y2 - x1):**
  - x1(0) = a~ |2 - |e~||;
  - y2(0) = sgn(a~ e~) sqrt(|e~| (1 - mu) / |x1(0)|);
  - x2(0) = y1(0) = 0.
- **Jacobi constant:** C_J = C_printed - mu (1 - mu). This follows from H = (u^2+v^2)/2 - (x1^2+x2^2)/2
  - (1-mu)/r1 - mu/r2 + mu x1, which gives -2H = C_J - mu^2. It agrees with the yaml header, "subtract
  mu(1 - mu)".
- **Traces:** Tr = trace(M4) - 2 and Tr_v = trace of the vertical 2x2 block. Both are invariant under the
  constant linear change of coordinates. Tr~ = (1 + log2|Tr|) sgn Tr when |Tr| > 2.
- **Units at production scale:** each system's registry mu, l* (km) and t* (s)
  (`cyclerfinder.core.cr3bp.cr3bp_system`). Distances are reported in km. Floors and SOI/Hill radii come from the
  satellite and planet registries.
  - No absolute constant in the numerics may assume a unit (CLAUDE.md, 2026-10-06 lesson).
  - Every tolerance is relative to the nondimensional state, and all of them are stated below.

## 3. Target mass ratios (the lead's set) and what the law predicts at each

| system | mu (registry) | i_k predicted to exist (mu'_k < mu < mu''_k; k = 1-4 printed, k >= 5 law) | in scope |
|---|---|---|---|
| Earth-Moon | 1.2150584e-2 | i_1 only (Varin computes it there: Table 3.6). No i_5-i_8 (mu > mu''_5) | yes; i_5-i_8 predicted ABSENT (a check, sec. 6) |
| Jupiter-Ganymede | 7.8049605e-5 | i_5 to i_10 (law); not i_4 (mu'_4 = 9.54e-5 > mu) | yes: i_5-i_8 |
| Saturn-Titan | 2.3669534e-4 | i_3, i_4 (printed), i_5, i_6 (law); not i_7 (mu''_7 = 2.05e-4 < mu) | yes: i_5, i_6 (i_3, i_4 recorded) |
| Neptune-Triton | 2.0895046e-4 | i_4 (printed), i_5, i_6 (law); i_3 NO (mu'_3 = 2.15e-4 > mu); i_7 at the edge (mu''_7 = 2.05e-4, law) | yes: i_5, i_6; i_7 decided by the computed mu''_7 |
| Pluto-Charon | 1.0876474e-1 | none (mu > mu''_1) | yes; every i_k predicted ABSENT (a check, sec. 6) |

- **Gating:**
  - The lead's dispatch names these five as the target set.
  - Titania (3.5e-5) and Europa (2.5e-5) are in R14's list but NOT in this set. They are excluded unless the
    owner or the lead opens them.
  - The Sun-Neptune mu 5.178e-5 (Voyatzis-Kotoulas) is used only as a shape control, if run (sec. 5).
  - No system here needs an owner ruling to RUN. Any catalogue row needs an owner decision.
- **The predictions above are the law's.** The computed mu'_k and mu''_k decide, and a target mu within 5% of
  a computed end is reported as "at the edge".

## 4. Method

All new code goes in `src/cyclerfinder/search/varin_closed_families.py`. It reuses the project's planar CR3BP
propagator and STM (`cyclerfinder.core.cr3bp.propagate`, DOP853, rtol = atol = 1e-12) and the multiple-shooting
corrector (`cyclerfinder.search.cr3bp_multiple_shooting`) for members with |Tr| above about 1e4.

1. **Critical-point solver (Varin 2000, 2009 eq. (2)-(3)).**
   - Unknowns: x1(0), y2(0), T, mu.
   - Equations: x2(T/2) = 0, y1(T/2) = 0, m_41(T) = 0 and m_14(T) = 0. The m_ij are elements of the 4x4
     monodromy in Varin's canonical coordinates, obtained from the project's STM by the constant linear map
     of sec. 2.
   - Newton with a finite-difference Jacobian. A central difference over the four unknowns, step 1e-7
     relative, is checked against a 1e-6 step. Convergence: all four residuals below 1e-11.
   - It finds saddles (mu'_k) and elliptic points (mu''_k). Which one is identified from the Hessian sign of
     the surface (the m_41 and m_14 derivatives), and recorded.
2. **Seeding k = 5-8.**
   - The seed is (a~, e~, T~, mu) from the sec. 1 trends, at both ends.
   - If Newton does not converge from the seed, the fallback is mu-continuation.
     - Continue the k = 4 critical point's defining equations, minus the mu unknown, along mu, holding
       m_41 = 0: the "y2-extremum curve".
     - Look for the next saddle by the change of sign of m_14.
   - The fallback is recorded as used or not.
3. **Loop at a target mu.**
   - From the saddle orbit at mu'_k, follow the y2-extremum curve (x2 = y1 = m_41 = 0; 3 equations,
     4 unknowns; pseudo-arclength) in mu to the target mu. That gives one member of i_k at the target.
   - Then trace the oval at fixed mu by pseudo-arclength continuation in (x1(0), y2(0), T) under
     x2(T/2) = y1(T/2) = 0.
   - Step control: arclength 1e-4 to 1e-2 (nondim), halved on a failed corrector.
   - A member whose half-period |Tr| exceeds 1e4 is re-corrected by multiple shooting (8 arcs per T~ unit) and
     its residual recorded.
4. **Per member, recorded:**
   - C_J, T (TU and days) and Tr, Tr_v (stability: |Tr| <= 2 in-plane stable);
   - min r2 (km) against P2's floor, SOI and Hill radius; min r1 (km) against P1's floor;
   - max r1;
   - perpendicular crossings (it is symmetric by construction);
   - the inertial winding about P1.
5. **Symmetric and asymmetric branches.**
   - The loops are symmetric by construction.
   - Where Tr crosses +2 along a loop (a symmetry-breaking pitchfork candidate), one asymmetric branch is
     followed with the #1000 eta method. It is labelled "descriptive", because its control in #1000 was
     self-referential.
   - At most 50 continuation steps per branch. Only existence and the first members are recorded.
6. **Cycler-class test per member** (the project's #855 reading):
   - P2 periapsis between P2's floor (radius + safe altitude) and its SOI;
   - P1 periapsis above P1's floor;
   - one P2 pass per cycle.
   - R14 says the birth orbits pass P2 at about 1.1 Hill radii, outside the SOI. So the test is per member, not
     per loop.

## 5. Controls (run first; the hard stops are marked)

**Positive control P1 (hard stop): reproduce Varin's printed members at his mu, to his digits, before any
continuation.**
- (a) Each of the 8 Table 1 orbits, converted by sec. 2 and integrated at its printed mu with no correction,
  gives |y1(T/2)| <= 1e-8 at the x2 = 0 crossing nearest pi T~. (The digest reached 5e-9 with an independent
  script.)
- (b) C_J of each reproduces C_printed - mu(1 - mu) to 1e-8.
- (c) Tr~ = 2 to 1e-4 (the saddle points are sensitive). Tr~_v reproduces the printed value to 1e-6.
- (d) The critical-point solver, started from each printed orbit with its four unknowns perturbed by 1e-4
  (relative), converges back to the printed (x1(0), y2(0), T~, mu). Agreement to 5e-8 relative, the printed
  8-digit rounding, for all 8.
- Failure of any of (a)-(d) on any row stops the task here.

**Positive control P2 (loop tracer, hard stop):**
- At mu_M (Earth-Moon, Table 3.6), i_1 exists.
- Two independent starts must meet the same closed oval. One continues the y2-extremum curve up from Table 1's
  mu'_1 saddle to mu_M. The other continues down from mu''_1's elliptic orbit.
- Both trace a loop that closes on itself: return to within 1e-6 of the starting (x1(0), y2(0), T) after one
  circuit, every member closing to 1e-9 (multiple shooting where needed).
- The two loops coincide: every member of one lies within 1e-5 of the other, in (x1(0), y2(0)) at matched T.
- Descriptive, not gating: which of Table 3.6's reproducible critical orbits lie on the loop. Varin's tables do
  not say which rows belong to i_1.

**Negative control N1 (hard stop): a closed family must not be "found" where Varin's printed data says it
does not exist.**
- Continue i_2's y2-extremum curve from mu'_2 toward Earth-Moon's mu (1.215e-2 > mu''_2 = 5.27e-3).
- It must end at a fold at mu''_2, reproducing Table 1's elliptic orbit to 1e-6 relative. It must NOT reach
  1.215e-2.
- If the tracer returns a "closed i_2 loop" at Earth-Moon mu, the method is broken: stop.

**Negative control N2:** the same for i_1 toward Pluto-Charon's mu (0.1088 > mu''_1): it must end at mu''_1.

**Shape control S1 (optional, not gating):** at the Sun-Neptune mu 5.178e-5, the Voyatzis-Kotoulas EXTERIOR
closed family (p = -7, the 6/7-7/8 band) is figure-only. If time allows, a closed loop in that band is sought
and compared by shape. It is reported, never used as a pass or fail.

## 6. Pass criteria

- **"Closed family i_k exists at mass ratio mu"** requires all of:
  1. The critical-point solver finds the saddle mu'_k and the elliptic point mu''_k, with mu'_k < mu < mu''_k.
  2. A member reached at mu by the y2-extremum continuation from mu'_k.
  3. The fixed-mu pseudo-arclength trace from that member returns to it: |Delta(x1(0), y2(0), T)| < 1e-6
     after a circuit of at least 20 members.
  4. The trace does not join family i: no member within 1e-4 of a family-i member at the same mu. This is
     checked against a family-i continuation over the loop's x1(0) range.
  5. Every member closes to 1e-9, by single or multiple shooting.
  6. The loop is the k-th one: T~ within [k, k+1] at both critical ends, and the inertial winding count about
     P1 equals the one Varin's k-th family has at mu'_k (taken from the k = 1-4 printed orbits; checked on the
     reproduced Table 1 orbits).
- **"Absent at mu":** the critical points are found and mu lies outside (mu'_k, mu''_k), or the extremum curve
  folds before mu (as in N1). "Not found" (the solver fails at the seed and at the fallback) is reported as
  "not found by this method", never as "absent".
- **"Collides with a catalogue row":** a catalogue row collides when ALL of these hold:
  - it is a CR3BP row of the same primary-secondary pair, with mass_ratio within 1e-6 relative of the registry mu;
  - it has a state;
  - its orbit lies on the loop: same C_J within 1e-6, and the row state within 1e-6 (nondim) of a loop
    member's trajectory at that C_J;
  - or it has the same period-topology (T~, winding) within 1e-3 and C_J inside the loop's C_J range: a
    class-level overlap, reported separately.
- **Collision candidates in the catalogue today:**
  - Earth-Moon CR3BP rows (Casoliva, Vaquero, Ross & Roberts-Tsoukkas, Schwaniger, the #1047 rows), against i_1
    at mu_M;
  - Pluto-Charon CR3BP rows (ross-rt-pc-cycler-32-2026 and the 5,1 row), recorded though no i_k is predicted;
  - no single-moon CR3BP row with a state at Ganymede, Titan or Triton. The Russell-Strange and multi-moon rows
    are patched-conic or multi-moon models, so class-level only.
- **Literal-collision check against published orbits:**
  - Varin's own Tables 3.x where the mu matches (Earth-Moon, i_1);
  - Voyatzis-Kotoulas 2005 (exterior, so no interior collision expected);
  - Colombo, Franklin and Munford 1968 (two closed families at mu_J; NOT held, wanted row 56): any loop at a mu
    near mu_J is flagged as unchecked against them. No target mu is near mu_J.

## 7. What each outcome means

- **i_5 (and up) found at Ganymede, Titan or Triton, with a cycler-class member:** a candidate row set for the
  owner. The #875 reading would be (iii) with R14's caveat: "first computed; predicted by Varin & Bruno's
  empirical asymptotics". The literature step (#1045 v3 gate with topology labels, plus the literal table) comes
  first.
- **Found, but no member is cycler-class** (every P2 pass outside the SOI): a census entry in this note, no row
  proposal.
- **The computed mu'_k and mu''_k miss the k^(-8/3) law by more than 20%:** reported as a finding about the
  law. The target table of sec. 3 is then redone with the computed values.
- **Not found** (solver and fallback fail for k = 5): the method's limit, recorded as "not found by this
  method". It is entered in the negative-results registry conditional on the method, not as absence.
- **Earth-Moon:** i_1 at mu_M is Varin's own published family, so known. The check is whether any catalogue row
  lies on it (sec. 6). Its members' cycler-class status is recorded.
- **Pluto-Charon:** expected absent for all k. If something is found, it means the law or the solver is wrong
  there, and the run stops for review.

## 8. Runtimes and who runs what

- P1: 8 orbits, each with a critical-point solve (about 4 x 5 STM propagations per Newton step; T~ up to 5).
  Under 2 minutes. Mine.
- P2, N1, N2: extremum-curve continuations over mu in 1e-3 to 3.7e-2 (about 100-300 steps of T~ 1-3 orbits)
  and two loop traces at mu_M. Estimated 5-15 minutes. Mine if a pilot shows under 8 minutes per call
  (checkpointed); otherwise handed to the lead.
- k = 5-8 critical points: T~ 5-9 (30-56 TU), each Newton step about 20 STM propagations. About 1-5 minutes
  per point, 8 points. Mine, in calls under 8 minutes.
- Loops at the targets: Ganymede i_5-i_8, Titan i_5-i_6, Triton i_5-i_6 (i_7 if it exists). That is 8-9 loops
  of 200-1,000 members with T~ 5-9, many highly unstable (multiple shooting). Estimated 0.5-3 hours each, so
  4-20 hours in total. Handed to the lead as one checkpointed command, Pool(2), progress every member chunk.

## 9. Order

1. This note committed.
2. Code and tests: unit tests at production scale (the registry mu, l*, t* of each target system; distances
   in km) and the P1 control as a test.
3. P1, then P2, N1, N2, with a stop on any failure.
4. The k = 5-8 critical points.
5. The loop command to the lead.
6. The results and the literature step in a dated note.

## 10. Control results so far (2026-10-11, 20:2x-21:2x AEDT): N1 is ambiguous, stopped for a ruling

- **P1: PASS** (4fea9c2d, `data/975_varin/p1_control.json`). All 8 Table 1 orbits:
  - (a) |y1(T/2)| <= 4.8e-9;
  - (b) Jacobi with the mu(1 - mu) convention <= 6.8e-9;
  - (c) Tr~ = 2 to 1.1e-5, and Tr~_v to the printed digit;
  - (d) the critical-point solver returns to the printed values to <= 1.0e-8 relative.
  - The half-period form of the solver (determinants of the half-residual Jacobian; the same conditions
    as m_41 = m_14 = 0) also returns to the printed values. It is 2-3 times faster.
- **N1** (`data/975_varin/curves/n1_a.jsonl`, `n1_b.jsonl`, `loops/n1_b_loop*.jsonl`). The y2-extremum curve
  leaves the mu'_2 saddle in TWO directions, and both go up in mu. The pre-registration assumed one curve.
  - Branch a (y2 increasing): folds at mu = 5.2720e-3. The fold, refined by the critical-point solver, is
    Table 1's mu''_2 elliptic orbit to 1.0e-8 relative. As the control requires.
  - Branch b (y2 decreasing): continues to the Earth-Moon mu (reached at T~ 2.049, x1 0.9102, y2 0.5073).
    The N1 text "it must NOT reach 1.215e-2" is literally violated by this branch. My reading: from a saddle
    the m_41 = 0 locus has two branches; one is the oval's extremum (it folds), the other runs along the
    family-i piece. Reaching mu_M on it is expected and is not a closed i_2.
  - The pre-registered failure condition is "the tracer returns a closed i_2 loop at Earth-Moon mu". From
    branch b's endpoint the fixed-mu trace does NOT close.
    - A first version reported "closed", but with a closure of 8.2e-4. Its closing step was inexact: a
      projected step, not a landing. That run is kept as `n1_b_loop_inexact_closing.jsonl`, and the logic is
      fixed.
    - The fixed trace passes 6.3e-3 from its start at member 228 without landing on it (the landing
      refinement fails to reach 1e-6). It then continued for 927 members without returning. So no closed loop
      by the sec. 6 criterion.
- **Disclosure:** before the controls finished I ran one exploratory k = 5 saddle solve (half form) from the
  sec. 1 extrapolated seed, to time the solver. It converged in 11 s to mu = 5.0379e-5, T~ 5.2756,
  Tr = 2.0000000003. This is not used, and it is not a result until P2/N1/N2 pass.
- Runtime: curves take 4-7 minutes at k = 2. A loop trace at mu_M runs more than 900 members in 7 minutes
  without closing. The k >= 5 critical points are much faster than the sec. 8 estimate.
- **Ruling requested:** amend N1 to "the oval branch from mu'_k folds at mu''_k and reproduces Table 1's
  elliptic orbit; any other branch that reaches the test mu must not yield a closed loop by the sec. 6
  criterion" (passed on present evidence), or treat N1 as failed and stop.

### 10.1 Post-hoc amendment A1 to N1 and N2 (lead ruling (a), 2026-10-12)

N1 (and N2 alike) now reads: "the oval branch from mu'_k folds at mu''_k and reproduces Table 1's elliptic
orbit (to 1e-6 relative); any other branch that reaches the test mu must not yield a closed loop by the sec. 6
criterion".
- **Why:** a saddle of the characteristic surface emits two y2-extremum branches, both going up in mu. The
  pre-registration assumed one. The real failure condition was always a closed i_k loop at the test mu, and
  the fixed-mu trace does not produce one.
- **The labelled false positive:** the first "closed" trace, with a closure of 8.2e-4 from an inexact closing
  step, is kept as `loops/n1_b_loop_inexact_closing.jsonl`.
- **N1 under A1: PASS.**
- The exploratory k = 5 solve stays disclosed and unused until the controls pass.

### 10.2 N2 (under A1): PASS; P2 curves done, loop traces to the lead

- **N2, from the mu'_1 saddle toward the Pluto-Charon mu 0.10876:**
  - Branch a folds at mu = 3.6678e-2. Refined, it is Table 1's mu''_1 elliptic orbit to 2.5e-9 relative.
  - Branch b reaches the Pluto-Charon mu. Its fixed-mu trace does not close: 989 members, never back within
    0.34 of its start, then the corrector failed.
  - So **N2 passes** (`curves/n2_a.jsonl`, `n2_b.jsonl`, `loops/n2_b_loop.jsonl`).
- **P2 curves:**
  - `curves/p2_up.jsonl`: from the mu'_1 saddle, the oval branch (y2 increasing, the branch that folds in
    N2), reaches mu_M at x1 0.738495, y2 1.229627, T~ 1.67547.
  - `curves/p2_down.jsonl`: from the mu''_1 elliptic orbit downward, reaches mu_M at x1 0.326705,
    y2 2.190603, T~ 2.02719.
  - The two loop traces at mu_M and their comparison (`compare`: each start must lie on the other loop to
    < 1e-6 after refinement) are handed to the lead.

## 11. P2 as traced, a tracer defect, N1 and N2 reopened (earthmoon-975-opus, 2026-10-12, 13:40-14:00 AEDT)

Owner of this section: earthmoon-975-opus (took over #975 from earthmoon-opus after the session restart).
Written before any code change and before any re-run. Evidence: the committed traces under `loops/`.

### 11.1 P2 as traced by the lead (launched 09:20 AEDT)

- `loops/p2_down_loop.jsonl` closed: 129 members, closure 1.6e-13.
- `loops/p2_up_loop.jsonl` did not close: 500 members, then killed by the session restart.
- `compare`: the start of each lies on the other trace (7.8e-11 and 8.6e-13).
- **The P2 wording is "both traces close and each start lies on the other's trace".** As traced, the second
  half passes and the first half FAILS. P2 is not passed on this evidence.

### 11.2 Why the up trace left the oval: a corrector jump at member 66

- Distance of each up member to the down family, by correction onto the down family (`distance_to_loop`,
  seeded from the three nearest down members): members 0-65 lie on it to <= 2e-11. From member 66 on, the
  distance is 1.02 and grows.
- The step from member 65 to member 66 is 1.02 in p = (x1(0), y2(0), T~). The step size was 0.02. T~ jumps
  from 2.031 to 2.984.
- Member 65 (0.3282, 2.1857, 2.0308) is next to the down start (0.3267, 2.1906, 2.0272). There the oval has a
  hairpin: x1(0) is at its minimum and y2(0) at its maximum on the oval, so the tangent is nearly pure T. A
  0.02 predictor step overshoots the turn. The corrector then converges on the pseudo-arclength hyperplane to
  a solution on another family.
- That family is the 3-fold cover of a T~ ~ 1 orbit near P1: member 66 at T/3 also has zero half residuals
  (8e-12). The trace then followed it toward the P1 collision (x1 -> 0, y2 -> infinity, T~ -> 3), which is
  the steady y2 rise in the lead's log.
- **Cause: a tracer defect.** `trace_loop` accepts any corrector landing. It has no limit on the distance
  from the predictor and no limit on the tangent turn per step. It is not a bifurcation of the oval, not a
  sign convention and not a genuine second branch through the start.

### 11.3 The same defect in the N1 and N2 b-loop traces: N1 and N2 are REOPENED

The same scan (step length against the 0.02 cap, chord turn angle) on every committed trace:

| trace | jumps (member: step length) | landing |
|---|---|---|
| `p2_up_loop` | 66: 1.02 | 3-fold cover of a T~ 1 orbit (T~ 2.03 -> 2.98) |
| `p2_down_loop` | none (65: 0.0203, the y2-min tip) | - |
| `n1_b_loop` | 347: 0.066; 682: 2.00 | 682: a T~ 3.005 orbit at the same (x1, y2) as a T~ 1.006 member |
| `n2_b_loop` | 87: 1.53 | a T~ 0.500 orbit (T~ 2.00 -> 0.50) |
| `n1_b_loop_inexact_closing` | none | - |
| the six extremum curves | none (max chord turn 22 deg) | - |

- N1 and N2 passed under A1 because their b-branch traces did not close. Those traces jumped families, so
  "did not close" is not established. **N1 and N2 are reopened.** The A1 PASS lines in sec. 10.1 and 10.2 are
  withdrawn until the v2 re-traces below.
- The near-return of `n1_b_loop` at member 228 (6.3e-3 from its start) was checked: the family through
  members 200-260 passes 8.17e-4 from the N1 start at its closest, and its tangent there is antiparallel to
  the start tangent (dot -0.95). It is a different branch passing near the start, not a closure. It does not
  fail N1. The same 8.17e-4 is the "closure" of `n1_b_loop_inexact_closing`.
- The extremum curves (`curves/*.jsonl`) show no jump, but the follower has the same missing check.

### 11.4 Tool amendment A2 (pre-registered; tracer code changes after this commit)

- **Loop tracer:** a predictor-corrector step is accepted only if (a) the corrector converges (as now), (b) the
  corrected point lies within 0.5 |step| of the predictor, and (c) the tangent at the corrected point turns by
  at most 10 deg from the previous tangent. A rejected step is halved, as a failed corrector is now. The
  landing refinement at closure uses the same accepted-step test.
- **Closure test:** as now (back within 1e-6 of the start after at least 20 members), plus the closing
  tangent points the same way as the start tangent (dot > 0). This rejects a doubling-back. The analysis also
  checks that no member j coincides with member n - j.
- **Extremum-curve follower:** the same tests (b) and (c), in q = (x1, y2, T~, log mu).
- The unguarded step stays available for `distance_to_loop` and for a regression test built from committed
  up member 65 and its tangent: the unguarded 0.02 step lands about 1.0 away; the guarded step rejects it.
- If the 0.5 or 10 deg values make a trace stall (step below ds_min), the change is written here as an
  amendment before any result is read.

### 11.5 Re-runs and their expected results (pre-registered)

All new traces get new names (`*_v2`). The committed traces stay as evidence.

1. **P1:** re-run (`p1`). It does not call the tracer. It is re-run because the module changed. Expected:
   identical to `p1_control.json` (pass, same numbers to the printed digits).
2. **The six extremum curves** (`n1_a`, `n1_b`, `n2_a`, `n2_b`, `p2_up`, `p2_down` as `*_v2`) under the
   guarded follower. Expected: the same end type ("fold" or "reached mu_stop"). A "reached mu_stop" end equal
   to the committed end to 1e-8 relative. A fold end, refined by the critical-point solver, gives Table 1's
   elliptic orbit to 1e-6 relative, as in sec. 10.
3. **P2 v2:** `p2_up_loop_v2` and `p2_down_loop_v2` from the committed curve ends, same start directions as the
   lead's runs. Expected for each: it closes (closure < 1e-6, closing tangent dot > 0) after about 130-250
   members (more than 129 because the guard shortens steps at the two tips).
   - Coincidence, metric 1 (all members): each member of one trace lies within 1e-5 of the other trace's
     family by orthogonal correction (`distance_to_loop` seeded from the three nearest members).
   - Coincidence, metric 2 (the "matched T" wording): for each member whose unit tangent has |T component|
     >= 0.1, (x1(0), y2(0)) corrected at fixed T onto the other family, starting from the nearest member of
     the other trace, lies within 1e-5. Members with a smaller T component are near a T extremum, where fixed-T
     correction is ill-conditioned; they are covered by metric 1 only, and their count is reported.
   - P2 passes if both traces close and both metrics hold. Otherwise P2 fails and the task stops.
4. **N1 v2 and N2 v2:** from each b-curve end at the test mu (`curves/n1_b.jsonl` at mu_M, `curves/n2_b.jsonl`
   at mu_PC), a guarded trace in each direction (the lead's direction, x1 increasing, and the opposite), at
   most 1500 members each.
   - **FAIL (hard stop):** either direction closes by the 11.4 test.
   - **PASS:** both directions reach an open end. An open end is a collision limit: the crossing's distance to
     P1 (|x1(0)|) or to P2 (|x1(0) - 1|) falls below 0.02 and the speed |y2(0) - x1(0)| rises monotonically over
     the last 100 members.
   - **Inconclusive** (budget spent, or the corrector fails at ds_min): reported as inconclusive with the end
     state, and sent to the lead for a ruling. Not read as a pass.
5. **The cascade (sec. 4.2-4.3) stays stopped** until P2, N1 and N2 have their v2 verdicts. The exploratory
   k = 5 solve (sec. 10) stays disclosed and unused.

### 11.6 v2 results (2026-10-12, 13:49-14:40 AEDT; A2 code deae8e75; machine load 18-45 on 8 cores)

- **P1 re-run: PASS, identical.** `p1_control.json` differed from the committed file only in the "seconds"
  timings, so the committed file is kept.
- **P2 v2: PASS** (`loops/p2_up_loop_v2.jsonl`, `loops/p2_down_loop_v2.jsonl`, `p2_v2_coincidence.json`,
  check `scripts/check_975_p2_coincidence.py`):

| trace | members | closure | A2 rejections (displacement / turn / corrector) | metric 1 max | metric 2 max (n checked, n skipped) | max member half residual |
|---|---|---|---|---|---|---|
| up v2 | 175 | 8.7e-13 | 1 / 27 / 4 | 1.2e-10 | 2.8e-10 (124, 51) | 9.4e-11 |
| down v2 | 174 | 8.4e-14 | 2 / 29 / 3 | 1.2e-10 | 6.1e-11 (123, 51) | 9.5e-11 |

  - Both close with a forward closing tangent. No doubling-back pair. Member counts are inside the
    pre-registered 130-250.
  - Both independent starts (the mu'_1 saddle continued up, the mu''_1 elliptic orbit continued down) meet
    the same closed oval at the Earth-Moon mu. The P2 wording is met in both halves.
- **N1 v2 and N2 v2: INCONCLUSIVE by the 11.5 rule.** No trace closed. No trace reached an open end as
  defined (crossing within 0.02 of a primary). Each ended with "corrector failed" (step below 1e-6):

| trace | members | end state p | min r1 / min r2 along the orbit (nondim) | half residual recomputed without STM | residual at DOP853 tol 1e-10 / 1e-12 / 1e-13 / 3e-14 (y1) |
|---|---|---|---|---|---|
| N1 plus | 371 | (0.98546, -0.68719, 0.99133) | 9.3e-4 / 1.5e-2 | 1.3e-8 | 1.8e-6 / 1.3e-8 / 3.5e-9 / 1.5e-9 |
| N1 minus | 124 | (0.99416, -1.29735, 2.04247) | 2.0e-2 / 5.9e-3 (at member 114) | 4.9e-8 (member 114) | 1.8e-5 / 4.9e-8 / 2.4e-8 / 4.6e-9 |
| N2 plus | 131 | (0.93033, -1.01781, 2.01318) | 1.1e-2 / 7.0e-2 | 1.1e-7 | 2.7e-5 / 1.1e-7 / -6.1e-9 / -1.0e-8 |
| N2 minus | 62 | (0.88504, -0.71673, 0.95147) | 3.8e-4 / 1.2e-1 | 2.3e-9 | 3.6e-7 / 2.3e-9 / -3.4e-9 / -2.4e-9 |

  - At every stall the integration error of the half-period map (1e-9 to 1e-7) is above the corrector's
    1e-10 target. Each stall is a precision wall. In three of the four the orbit passes within 1e-3 to 2e-2
    of P1 or P2 mid-orbit, and the half-residual Jacobian is large (largest singular value up to 2e6).
  - These are not shown to be family ends. The 11.5 open-end test looks only at the crossing, so it cannot
    see a mid-orbit close pass.
  - The v1 traces went "further" only because they jumped to other families at or near the same places
    (n2_b_loop jumped at member 87, next to the N2-plus wall).
- **Ruling requested from the lead** before any further N run (11.5 item 4: inconclusive goes to the lead).
  The cascade stays stopped.

### 11.7 The six extremum curves under the guarded follower: unchanged (2026-10-12, 14:40-15:05 AEDT)

`curves/*_v2.jsonl`, comparison `curves_v2_compare.json`. Same starts, directions and mu_stop as the
committed curves.

| curve | end (v1 / v2) | check | result |
|---|---|---|---|
| p2_up | reached mu_stop / same | v2 end against v1 end, max relative difference | 4.4e-13 |
| p2_down | reached mu_stop / same | same | 3.2e-13 |
| n1_b | reached mu_stop / same | same | 1.2e-11 |
| n2_b | reached mu_stop / same | same | 3.9e-11 |
| n1_a | fold / fold | fold refined by the critical-point solver (half form, 3 iterations) against Table 1 mu''_2 | 1.0e-8 |
| n2_a | fold / fold | same, against Table 1 mu''_1 | 2.5e-9 |

All six meet the 11.5 item 2 expectation. A2 rejected 0-6 steps per curve, all for the turn limit. The
oval-branch parts of N1 and N2 (the folds at mu''_k) stand under A2. Only the b-branch traces at the test mu
are open (sec. 11.6).

### 11.8 Correction to 11.6: the four stalls are of two kinds (2026-10-12, 15:05-15:10 AEDT)

The 11.6 table gave min r1 and min r2 from 200,001 time samples. Near P1 those samples are 1.5e-3 to 2e-3
apart in position, so the near-P1 values were upper bounds. Event-located periapses (DOP853 rtol 3e-14,
atol 1e-16; `data/975_varin/live/v2_periapsis.log`, not tracked) for the last members:

| trace | true periapsis at the stall | trend along the last 60 members | kind |
|---|---|---|---|
| N1 plus | r1 5.27e-5 (about 20 km from Earth's centre, inside the Earth) | r1 0.26 -> 1.1e-4 -> 5.3e-5, then flat as the steps shrink | collision type |
| N2 minus | r1 1.45e-4 (inside Pluto) | r1 0.17 -> 3.0e-4 -> 1.45e-4, then flat | collision type |
| N1 minus | r2 5.85e-3, r1 1.83e-2 | r2 levels off at 5.85e-3 | sensitivity type (largest Jacobian singular value 2e6) |
| N2 plus | r1 1.06e-2, r2 6.97e-2 | r1 levels off at 1.06e-2 | sensitivity type (largest Jacobian singular value 7e4) |

- **Diagnostic, not pre-registered:** the last member of each sensitivity-type stall was re-corrected with
  the project's multiple-shooting corrector (`cr3bp_multiple_shooting.correct_multiple_shooting`; 8, 16 and
  32 nodes from one DOP853 rtol 1e-13 arc). All six runs converge in 2 iterations to a full-state closure of
  1e-13 to 9e-13, and the node moves by 2e-13 to 4e-12. So these members are periodic orbits to 1e-13. The
  tracer stalls because its single-shooting half-period residual carries the integration error times the
  sensitivity (1e-12 x 1e4-1e6), not because the family ends.
- **Independent residual check on P2 v2:** every member's half residual, recomputed with the plain
  propagation (no STM; the tracer used the STM-augmented one), is at most 9.0e-11 (up) and 9.2e-11 (down). The
  P2 PASS stands.
- **For the resume (lead ruling 14:45: (c) rejected; (a) diagnostic only; (b) regularised propagation for the
  verdict):**
  - Regularisation (Levi-Civita) is what can carry the two collision-type traces through or into their
    collision orbits.
  - The two sensitivity-type traces need a symmetric multiple-shooting corrector inside the tracer (sec. 4.3
    already names multiple shooting; its |Tr| > 1e4 threshold would need an amendment, since |Tr| is 20-600
    at these stalls). Regularisation alone does not fix them.
  - The cascade at T~ 5-9 will meet both kinds.
- Housekeeping at the pause: the scripts preflight tests pass after the check script was added. Nothing was
  running for #975 at 15:10.
