# #884 -- Which catalogued Earth-Moon cycler families survive the Sun (BCR4BP)

**Date:** 2026-10-04. **Code:** `src/cyclerfinder/search/sun_forced_periodic_884.py`, gates in
`tests/search/test_sun_forced_periodic_884.py` (13 tests, about 15 s with 4 workers), staged driver
`scripts/screen_884_sun_forced_em_cyclers.py`, outputs under `data/found/884_sun_forced_em_cyclers/`
(`controls/`, `families/`, `melnikov/`, `continue/`, `verify/`, `summary.json`; `continue_ds0.1/` and
`verify_ds0.1/` are a superseded coarse-step run kept as the record of a branch jump, see section 5).
No catalogue row is written; the coordinator decides.

**Published reference:** Brown, Peterson, Henry & Scheeres, "Structure of periodic orbit families in
the Hill restricted 4-body problem", SIAM J. Appl. Dyn. Syst. 24(1):346-375 (2025), filed in the
private paper corpus as
`brown-peterson-henry-scheeres-2025-periodic-orbit-families-hill-restricted-4-body-problem-siads-24-1-346-doi-10.1137-24M1637301-published.pdf`
(digest: `docs/notes/2026-10-03-digest-brown-2024-hr4bp-periodic-orbit-families.md`).

## 1. Short answer

- Every catalogued Earth-Moon family tested at a low-order commensurate member (a = 1 or 2, see
  section 2) has Sun-forced "dynamical equivalents" at physical Sun mass in the bicircular model. The
  survey covered 13 commensurate members of 9 families. At each member the equivalents sit exactly at
  the phases fixed by the reversing symmetry: Sun on the Earth-Moon line or perpendicular to it at the
  orbit's x-axis crossing.
- Not every phase survives:
  - For Ross-RT C21 (3:1), the two equivalents at the near-perpendicular phase fold back at
    eps = 0.389 and reconnect to a different CR3BP orbit, at C = 3.1214 against the start's 3.1252.
  - For one Ross-RT C11 member (5:2), one of its two equivalents does the same, folding at eps = 0.398.
- All forced cycler-class equivalents are linearly unstable, as their CR3BP parents already are. The
  only linearly stable forced orbits found belong to the two families whose CR3BP member is itself
  stable: R21-S at 1:2 and Casoliva 1:2(d) at 3:2. Neither is cycler-class.
- The periselene stays inside the lunar SOI for every cycler-class family, so all of them stay
  cycler-class. One exception: the Casoliva 2:1(b) low-perilune member's forced equivalents reach
  periselene 1,440-1,453 km. That is below the lunar radius (1,737 km), so they pass through the Moon.
- The coordinator's a = 3 target (C32 at 8:3) is a much weaker resonance, about 0.5 percent of the
  a = 2 Melnikov amplitude. Its continuations advance slowly. One of four branches reached eps = 1
  in the superseded coarse run but has not been re-confirmed under the guarded stepping. The rest are
  incomplete; see section 7 for the command.

## 2. Model difference, and what transfers from Brown et al.

**Their model.** The Hill restricted four-body problem (HR4BP) is coherent: the Earth and Moon move
on Hill's variation orbit, which is perturbed by the Sun. It is parameterised by m (m_SEM = 0.0808).
Time is scaled so that 2 pi = one synodic month. The equations are pi-periodic: the solar term is the
pure quadrupole tide, (3/4) m^2 [(x^2 - y^2) cos 2 tau - 2 x y sin 2 tau], plus the variation of the
primaries, which is also a 2-tau effect. So their forcing period Tg = pi is HALF a synodic month.

**Ours.** We use the standard incoherent BCR4BP (`core/bcr4bp.py`, Andreu / Rosales-Jorba constants).
The Sun is on a circle; its direct term plus the indirect term acts on the particle only, and the Moon
does not feel the Sun. The forcing period is the full synodic month:

    Tg = 2 pi / omega_S = 6.791193872 TU   (omega_S = 0.925195985520347)

Ross-RT, Braik-Ross and Casoliva rows use a slightly different TU convention; the 6.79117 TU in the
brief is that convention. Our homotopy multiplies the whole solar term by eps (mu_sun = eps * mu_sun_phys)
and holds omega_S fixed.

We did not use the QBCP (`core/qbcp.py`). It is coherent, and so closer to the HR4BP, but it has no
single scalar that switches the Sun off. Its eight alpha functions would have to be homotoped
individually.

Separately, and out of scope: the QBCP alpha_1 table in `core/qbcp.py` contains the literal
`-38.068581391005552e-08`, which looks like a transcription slip worth checking against Gimeno-Jorba
2018 Table 4.

**What transfers.**
1. *Commensurability.* The necessary condition a T* = n Tg holds for any time-periodic forcing, so it
   transfers. But Tg is different: our Tg = 2 x theirs.
2. *The phase function.* In both models the first-order bifurcation function is the solar work along
   the unperturbed orbit (Brown's Eq. 2.7):

       Mel(theta0) = -2 int_0^P v . a_sun(x(t), theta0 + omega_S t) dt

   This is the projection of the forcing onto the cokernel of M^a - I, which is grad C. Shifting the
   start point is equivalent to shifting theta0 (their Proposition 1). If T* = (n/a) Tg, then Mel is
   2 pi / a periodic in theta0, so only Sun-angle harmonics k that are multiples of a contribute.
3. *Which harmonics exist.* In the BCR4BP, the direct-minus-indirect potential is
   mu_S / a_S * sum_{l >= 2} (r / a_S)^l P_l(cos gamma), whose Sun-angle harmonics are:
   - k = 0, 2 from the quadrupole, which matches the HR4BP tide;
   - k = 1, 3 from the octupole, about r / 389 weaker; the HR4BP has no octupole;
   - k = 0, 2, 4 from l = 4, about (r / 389)^2 weaker; and so on.

   Hence:
   - a = 1 or a = 2 resonances get a quadrupole-strength first-order Melnikov function (Brown's
     a = 1 is our a = 2: our T* = Tg/2 is their T* = pi);
   - a = 3 gets octupole strength only;
   - a >= 4 is negligible at first order.
4. *Zeros at symmetric phases.* Both models have the reversing symmetry
   (x, -y, z, -vx, vy, -vz, -t, -theta). For an x-axis-symmetric orbit started at its perpendicular
   crossing, Mel is odd in theta0. So zeros at theta0 = 0 and pi are forced, and for a = 2 also at
   pi/2. That is their Proposition 3, and it transfers exactly. Non-symmetric zero locations and
   amplitudes do not transfer: the BCR4BP omits the Sun's effect on the Moon (the variation), which is
   the same order as the tide.
5. *Specific orbits do not transfer.* Their homotopy in m also changes the forcing frequency, so their
   starting CR3BP orbits sit at T* = k pi in sidereal units. Ours sit at multiples of Tg/2 = 3.3956 TU.
   The same named family member (for example "L2 halo at T = pi") is therefore a different orbit in
   the two studies. Only structure (zero counts, symmetric phases, folds and reconnections) can be
   compared.

## 3. Gates (measured)

All in `tests/search/test_sun_forced_periodic_884.py`, all passing. Expected values are identities,
independent code, or the published symmetry prediction.

| Gate | Result |
|---|---|
| eps = 0 RHS vs `core/cr3bp.cr3bp_eom`; eps = 1 RHS and STM RHS vs `core/bcr4bp` (separate code) | agree to 1e-13 / 1e-12 |
| Jacobi conserved with Sun off, 3 Tg on the L1 member | < 1e-10 |
| Map inverse (Tg forward then back, eps = 1) | < 1e-10 (5.9e-12 measured) |
| STM vs central finite differences; dX/deps vs finite eps | 1e-6 / 1e-5 relative |
| Monodromy symplectic in canonical coordinates p = v + omega x r | relative defect < 1e-10 (the Cartesian STM fails the same check, as it should) |
| Melnikov: Simpson quadrature vs eps-sensitivity ODE vs central finite-eps Jacobi change | 1e-8 and 1e-6 relative |
| Mel(theta0 + 2 pi / a) = Mel(theta0) (a = 2) | 1e-9 relative |
| Melnikov zeros at the symmetric phases (L1 Lyapunov, a = 2) | within 1e-7 of {0, pi/2} |
| Forced orbit: multiple-shooting residual; Radau closure of every segment | < 1e-10; < 1e-7 |
| Continuation monotone in eps (0 to 0.2) and reversible back to the start nodes | < 1e-6 |
| Symmetric forced orbit vs independent half-period perpendicular-crossing shooting in the forced model | 1e-9 |

In the survey (section 6), every eps = 1 orbit closes to between 4e-15 and 7e-11 (DOP853 multiple
shooting), and Radau reproduces every segment to 4e-10 or better. Reversibility was checked on 16
survivors (C11 x 8, Casoliva 2:1(b) x 8) plus C21 3:1 (two survivors and two fold-backs). Each
survivor's reversed branch lands at eps = 0 on its start orbit: within 2e-4 of the start nodes, with a
phase shift of at most 3e-5 TU and a fixed-point residual of at most 5e-11 (`verify/`).

The family selection is also guarded. The CR3BP walks accept a step only if the corrector lands
within half a step of the predictor.

## 4. Positive control against Brown et al.

Brown et al. print no initial conditions, so the control is qualitative. We used the two planar
Lyapunov members they list, placed at our commensurate periods (seeds: Braik-Ross 2026 Table 2 LL1 and
LL2). The L2 halo was set aside because at our Tg/2 = 3.3956 TU it sits next to the Lyapunov-halo
bifurcation, where an extra multiplier is near 1.

- **L1 Lyapunov, T* = Tg/2 (a = 2; their T* = pi case).**
  - *Prediction (Prop. 3):* zeros at the xz-plane crossing for Sun phases 0 and pi/2, giving two
    equivalents per pi.
  - *Measured:* exactly two zeros, at 3.1415926532 and 1.5707963267. Both continue to eps = 1, with
    closure 2e-14 and Radau 3e-14. Two off-zero phases (0.37 and 0.81 of the period) fail to start at
    any eps down to 1e-6, the analogue of their 100-point test.
  - The two equivalents differ in stability along the phase direction: the pi/2 branch has a
    unit-circle pair, the 0 branch a real pair. This is the alternation predicted by the sign of
    dMel/dtheta0.
- **L2 Lyapunov, T* = Tg (a = 1; their L4 / L2 T* = 2 pi cases have four zeros).**
  - *Measured:* four zeros. Two are exact (0 and pi). The other two are at pi/2 - 0.00142 and
    3 pi/2 + 0.00142, mirror images of each other under the reversing symmetry.
  - The 0.0014 rad offset is the octupole, which the HR4BP lacks. In the HR4BP, theta0 and
    theta0 + pi give the same object, so four zeros collapse to two objects, as they report for L4.
    In the BCR4BP, 0 and pi are distinct orbits, with Floquet magnitudes differing at the 1e-5
    level.
  - All four continue to eps = 1. Off-zero phases do not start.

**Outcome:** reproduced, in the sense in which the papers' claims transfer:
- zero counts;
- zero locations at the symmetric phases;
- only Melnikov zeros continue;
- a different model shifts only the non-symmetric zeros, and only at octupole order.

The fold-and-reconnect behaviour in section 6 is the analogue of their statement that "HR4BP families
can also connect two different CR3BP orbits". The 9:2 NRHO case is not testable this way: its
first-order Melnikov function is flat (a = 9), as Brown's own Prop. 2/3 discussion implies.

## 5. A branch jump caught and fixed

The first survey pass used eps-continuation steps up to ds = 0.1 with no jump guard, and every branch
"reached eps = 1". Rerunning one of them, Ross-RT C21 3:1 at theta0 = pi/2, with ds_max = 0.02 showed
something different: the branch folds at eps = 0.389 and returns to eps = 0. The coarse run had jumped
to another branch.

Reversibility at the same coarse step size did not catch this, because a jump can be retraced. The
fix has two parts:
- the corrector must land within half a step of the predictor;
- ds_max = 0.02.

The whole survey was rerun with the fix. The coarse run is kept in `continue_ds0.1/` and
`verify_ds0.1/` and is superseded. Under the guard, every coarse "survivor" except the C21 3:1
pi/2-pair re-converged to the same eps = 1 orbit, with the same Floquet magnitudes to 3 digits.

## 6. Survey

Each family was walked in the CR3BP (pseudo-arclength in (x0, [z0], vy0, T), symmetric half-period
shooting) from its catalogued row, collecting every crossing of T = k Tg/2 (k = 1..12) and of
T = 8/3 Tg and 5/6 Tg. Every member was re-corrected at the BCR4BP's mu = 0.0121505816.

**Family identities found on the way:**
- Braik-Ross C32 and Ross-RT C32 are ONE family. It folds in period (T 16.46-19.25 TU), so it crosses
  8/3 Tg twice: at C = 3.12875 (periselene 26,209 km) and at C = 3.18264 (4,457 km).
- Ross-RT C33 does not reach 8/3 Tg (its walk covers T 18.12-19.09).
- Braik-Ross C11a, C11b and Ross-RT C11 are one family.
- The #438 spatial C21 row and the #682 C21 3D corridor seed are one family.

"Phases" below are theta0 at the member's x-axis crossing. "Cycler" means periselene at most 66,182.9 km
(the lunar Laplace SOI used by `data/validate.py`), with lengths converted at 384,400 km.

| Family (seed row) | T*/Tg (a) | CR3BP member C, peri, max abs lambda | Zeros (equivalents) | eps = 1 reached | Forced stability (max abs lambda over P) | Peri at eps = 1 (cycler?) |
|---|---|---|---|---|---|---|
| C11 (braik-ross-c11a) | 2/1 (1) | 3.12107, 25,848 km, 3.8e5 | 0, pi/2, pi, 3pi/2 (4) | 4/4 | unstable, 4e5 to 5e5 | 24,678-26,303 (yes) |
| C11 | 3/2 (2), C = 3.15107 | 23,253 km, 1.0e3 | 0, pi/2 (2) | 2/2 | unstable, 3.8e7 / 1.4e7 | 21,855 / 23,999 (yes) |
| C11 | 3/2 (2), C = 3.09204 | 21,995 km, 1.4e5 | 0, pi/2 (2) | 2/2 | unstable, 3.2e10 / 1.0e10 | 21,671 / 22,115 (yes) |
| C11 (ross-rt-11) | 5/2 (2), C = 3.07692 | 19,954 km, 5.8e6 | 0, pi/2 (2) | 1/2: theta0 = 0 folds at eps 0.398, reconnects to a C = 3.0800 orbit | unstable, 2.5e14 | 19,726 (yes) |
| C21 planar (ross-rt-21) | 3/1 (1) | 3.12523, 26,007 km, 3.7e4 | 0, pi/2 + 4e-4, pi, 3pi/2 - 4e-4 (4) | 2/4: the pi/2 pair folds at eps 0.389, reconnects to a C = 3.1214 orbit | unstable, 9.1e5 | 24,403 (yes) |
| C21 3D corridor (#682 / #438 spatial) | 2/1 (1) | 3.09895, 42,003 km, 1.1e5 | 0, pi/2, pi, 3pi/2 (4) | 4/4 | unstable, 1.7e5 / 6.6e4 | 40,404-43,269 (yes) |
| C21 3D corridor | 5/2 (2) | 3.01251, 75,329 km, 37 | four zeros pi/4 apart; Mel amplitude 3e-7, i.e. the quadrupole term vanishes by symmetry (Brown's L1-vertical situation) | incomplete (eps 0.08-0.50) | -- | (no: 75,329) |
| C32 (braik-ross-c32 = ross-rt-32) | 5/2 (2), C = 3.15168 | 23,304 km, 4.3e3 | 0, pi/2 (2) | 2/2 | unstable, 1.1e7 / 2.5e7 | 21,741 / 23,312 (yes) |
| C32 | 5/2 (2), C = 3.18010 | 8,288 km, 80 | pi/2, pi (2) | 2/2 | unstable, 1.8e4 / 2.9e3 | 8,854 / 8,290 (yes) |
| C32 | 8/3 (3), C = 3.12875 | 26,209 km, 2.8e5 | pi/3, 2pi/3 (2); Mel amplitude 2.0e-4 vs 0.02 at 5/2 | incomplete (eps 0.18, 0.56); coarse run reached eps = 1 at 2pi/3 (unconfirmed) | -- | -- |
| C32 | 8/3 (3), C = 3.18264 | 4,457 km, 79 | pi/3, 2pi/3 (2); Mel amplitude 1.0e-4 | incomplete (eps 0.007, 0.016; slow) | -- | -- |
| C31 (ross-rt-31) | 5/2 (2) | 3.05583, 14,294 km, 6.2e3 | pi/2, pi (2) | 2/2 | unstable, 4.4e7 / 3.1e7 | 14,007 / 14,323 (yes) |
| C31 | 8/3 (3) | 3.02472, 10,208 km, 1.6e4 | 0, pi/3 (2); Mel amplitude 3.6e-5 | incomplete (eps 0.95, 0.90) | -- | -- |
| Casoliva 2:1(b) | 1/1 (1), C = -0.02092 | 13,843 km, 7.6 | pi/2 - 0.012, pi, 3pi/2 + 0.012, 2pi (4) | 4/4 | unstable, 7.8 | 13,367-13,441 (yes) |
| Casoliva 2:1(b) | 1/1 (1), C = 0.56729 | 1,810 km, 970 | pi/2, pi, 3pi/2, 2pi (4) | 4/4 | unstable, 950 | 1,440-1,453: below the lunar radius, impacts |
| R52-S planar (braik-ross) | 2/1 (1) | 3.19716, 74,124 km, 9.4 | 0, pi/2, pi, 3pi/2 (4) | 4/4 | unstable, 9.2 / 10.4 | 74,869-75,133 (no, resonant) |
| R52-S | 5/2 (2) | 3.02010, 78,421 km, 1.4e3 | pi/2, pi (2) | 2/2 | unstable, 1.7e6 / 2.7e6 | 78,616-79,747 (no) |
| R52-S | 8/3 (3) | 2.97116, 65,148 km, 1.1e3 | pi/3, 2pi/3; Mel amplitude 2.7e-5 | 2pi/3: incomplete (eps 0.75); pi/3: first step returned to eps = 0 (start ambiguous) | -- | -- |
| R21-S planar (braik-ross) | 1/2 (2) | 3.42072, 195,515 km, 1 (stable) | pi/2, pi (2) | 2/2 | pi/2: unstable 1.16; **pi: linearly stable** | ~194,000-197,000 (no) |
| Casoliva 1:2(d) | 3/2 (2) | 3.28578, 352,882 km, 1 (stable) | pi/2, pi (2) | 2/2 | pi/2: unstable 1.01; **pi: linearly stable** | ~348,000 (no) |

**Not reached or not applicable:**
- Vaquero 2:1 and R21-S at 5:6 (a = 6): the Melnikov amplitude is 5e-13 and 2e-13, i.e. flat at
  first order as the harmonic argument predicts. No continuation was attempted.
- Vaquero 3:1, Casoliva 1:2(c), 2:1(a), 3:2(c), Braik-Ross C11b's own walk, and R31-S: no
  half-integer-Tg member inside the range walked.
- Casoliva 7:3(b) and 7:3(c) have no perpendicular x-axis crossing (asymmetric) and were not treated.
- The #682 L1 Lyapunov 3D corridor family does not reach Tg/2 (walk T 1.80-3.12 TU).
- The C32 3D corridor (T 31.24-31.32) contains no target.
- Two walks (Ross-RT C21 downward and R21-S upward) exceeded the per-command limit in near-collision
  propagation. They are recorded as not completed, so C21 planar's 5/2 member, if any, is unknown.

**Pattern:**
- For a = 1 and 2, every equivalent of every cycler-class family reaches physical Sun mass, except
  the three fold-backs.
- The phase-direction stability alternates between neighbouring zeros. A forced orbit is linearly
  stable only if its CR3BP parent already was.
- The a = 3 (8:3) resonances have a Melnikov amplitude 0.3-1 percent of the same family's a = 2
  member, consistent with the octupole ratio r / a_S (r ~ 1-2 Earth-Moon units against a_S = 389).
  Their continuations move the fixed point a long way along the orbit per unit eps (C32 at
  C = 3.18264 drifts 0.3 in state by eps = 0.025), so they need very small steps.

## 7. Not done, and how to finish

- **Finish the a = 3 branches and the degenerate C21 3D 5/2 branches.** These are C32 8/3 (four
  branches), C31 8/3, C21 3D 8/3, R52 8/3, and C21 3D 5/2 (four branches). The stage is resumable:
  each invocation restarts every branch from its last converged point. Repeat until no `wall-clock`
  entries remain:

      uv run python scripts/screen_884_sun_forced_em_cyclers.py --stage continue --max-a 3 --max-steps 400 --budget-s 100

  Each invocation takes about 6 minutes with 4 workers. Estimates:
  - C31, C21 3D and R52: about 2-5 more invocations.
  - The C32 member at C = 3.18264 (periselene 4,457 km): it advanced about 0.01 in eps per
    invocation, so at the current step control it would need several hours. The faster route is to
    continue it at fixed eps in theta0, or to use theta0 as the continuation parameter.

  Then run `--stage verify` (fold-back landing and reversibility) and `--stage summary`.
- The two uncompleted family walks (section 6).
- Out-of-plane members of the planar families. Every planar family here was continued in the plane;
  the forced orbits' 6x6 monodromy includes the out-of-plane pair.
- A QBCP cross-check of the survivors. This means correcting the eps = 1 BCR4BP solutions directly in
  the QBCP. It is not done: the QBCP has no scalar homotopy.
- Literature: Komachi (ASC 2026, "ballistic cycler between EML2 and SEL2 in the BCR4BP") has not been
  obtained (title only), and `search/literature_check.py` has not been run. Same primary and model
  class as Brown et al., different orbit class. Labelling by spec 16.4 is left to the coordinator.

## 8. Where the brief was wrong or could be sharpened

- "Synodic month = 6.79117 TU" is the Ross-RT/Braik-Ross time unit. The BCR4BP's own value is
  6.791193872 TU, which is the one used here.
- The suggested targets (C32 at 8:3, Vaquero 2:1 at 5:6) are a = 3 and a = 6 resonances. In a
  quadrupole-dominated forcing these are first-order flat (a = 6) or octupole-weak (a = 3). The
  natural targets are the half-integer multiples of Tg (a = 1, 2), which is Brown's own choice
  (T = b pi, i.e. b half synodic months).
- The coordinator's "C32 brackets 8:3 at 78.61 d and 78.90 d" is in fact the C32 family crossing 8:3
  twice, because it folds in period. The 78.90 d row (Ross-RT C33) is on a different branch that does
  not reach 8:3.
- "Brown's forcing period is the synodic month": in the HR4BP, the forcing period is HALF the synodic
  month (Tg = pi with 2 pi = one synodic month).
