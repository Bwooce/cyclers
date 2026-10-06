# Digest: Chenciner and Llibre 1988, "A note on the existence of invariant punctured tori in the planar circular restricted three-body problem" (#960 batch 32)

Alain Chenciner (Universite Paris VII) and Jaume Llibre (Universitat Autonoma de Barcelona), Ergodic Theory and Dynamical
Systems 8* (the Conley memorial issue), 63-72 (1988). DOI 10.1017/s0143385700009330. 10 pp. Title and "8*, 63-72" confirmed on the page image.
- Source scan: `e1c9514c-a-note-on-the-existence-of-invariant-punctured-tori-in-the-plana-1988.pdf`, md5
  19740a056d7748e7c526de2a244e1814. It has a publisher text layer, but the theorem statements and the normal form are garbled in it.
- Proposed corpus filename:
  `chenciner-llibre-1988-note-existence-invariant-punctured-tori-planar-circular-restricted-three-body-problem-etds-8-63-doi-10.1017-S0143385700009330.pdf`
- How I read it: all theorem statements, the regularised Hamiltonian, the return-map normal form and the KAM theorem were read
  on 150 dpi page images (pp.64, 65, 69, 70). I did not rederive the proofs or the normal form. The only computations are two
  small checks (sec. 5, `check_chenciner.py`, `check_chenciner.out`).

## 0. Verdict

**This is the paper that turns transversal ejection-collision (EC) orbits into invariant punctured tori by KAM, and that proves
the EC orbits exist for every mass ratio.** It is the missing source behind the "four EC orbits for any mu" line quoted in the held
Rodriguez del Rio thesis digest (part a, line 271) and in the Llibre 1982 digest (line 213).
- It is analytic and perturbative: the Jacobi constant must be large, with no explicit bound. **No numbers** other than the
  coefficients of the return-map normal form (sec. 3).
- Hill's problem is covered by a remark (p.71, due to C. Simo).
- **Proposal only:** no catalogue row. Use it as the theory base for #899 (large-C tori and EC orbits). Never cite it for a
  specific C value or orbit.

## 1. Setting (pp.63-65)

- Planar circular restricted problem, complex coordinates `x = x1 + i x2`, `y = dx/dt + i*omega*x`, `omega = sqrt(G(mu + nu))`,
  masses `mu, nu` in R+ with `G = 1`, `mu + nu = 1` (so `omega = 1`). The zero-mass body orbits the primary of mass `nu` (as in Conley).
- Main parameter `epsilon` by the energy level `H = -1/epsilon^2`. The paper calls H "one half of" the Jacobi constant, so
  large Jacobi constant means small epsilon (my reading: C is about `2/epsilon^2`; the paper gives no such formula, and its H carries "the same constant 2 mu as Conley", so levels differ from other sources by a constant). For small
  epsilon there are three Hill regions. Epsilon = 0 is a collision.
- Levi-Civita regularisation `x = 2 z^2`, `y = w/(epsilon z-bar)`, `dt = 2 epsilon |x| dt'` gives
  `K(z, w) = {1 + 2 i epsilon (z-bar w - z w-bar)} |z|^2 + |w|^2 - nu epsilon^2 - mu epsilon^3 g(z)`, with
  `g(z) = 2|z|^2 {1/|2 z^2 + 1| - 1 + z^2 + z-bar^2} = |z|^2 {2|z|^4 + 3(z^4 + z-bar^4) + O(|z|^6)}`.
  - **Check:** I confirmed the series by numerical evaluation: the remainder over `|z|^8` tends to a constant (16.0),
    as an `O(|z|^6)` bracket requires. The printed leading terms are right.
- The level `K = 0` is a 3-sphere `|z|^2 + |w|^2 = nu epsilon^2`; the Levi-Civita map is a twofold cover of `H = -1/epsilon^2`
  away from the circle `z = 0`. Without the last term (`mu` terms) it is the rotating two-body problem: invariant tori labelled by
  angular momentum; the zero-angular-momentum torus contains `z = 0` and is made of EC orbits (Fig. 2).
- Annulus of section `A` bounded by the direct and retrograde quasi-circular periodic orbits (Conley 1963; Fig. 3). `P_epsilon` is
  the first return map on `A`. Blowing up `z = 0` (polar `z = r e^{i theta}`) gives McGehee-like variables and a collision
  torus `u^2 + v^2 = nu`, with two circles of equilibria `C+-` and their asymptotic manifolds `W+`, `W-` (Fig. 4).

## 2. Theorems, as printed

- **Theorem 1 (p.65).** If `mu*nu` nonzero and epsilon is small enough, the circle `z = 0` meets its image under `P_epsilon` in
  **exactly eight transversal points**. Remark (i): because of the twofold cover these are **four** transversal EC orbits in
  the original (McGehee) picture. Finiteness follows from analyticity of `P_epsilon` up to the boundary (remark iii).
  - Proof idea (pp.66-68): at epsilon = 0, `W+ = W-` is explicit: `u0 = 0`, `r0 = sqrt(nu)/cosh(sqrt(nu) t'')`, `v0 = -sqrt(nu) tanh(sqrt(nu) t'')`.
    The first-order shift of `W+` on the annulus `v = 0` is `u = phi(theta) = 12 mu nu^(5/2) (int_0^inf cosh^-7 x dx) epsilon^6 sin(4 theta) + O(epsilon^7)`.
    By time-reversal symmetry `W-` is the graph of `phi(-theta) = -phi(theta)`. They cross where `sin 4 theta = 0`: theta = k pi/4, k = 0..7,
    eight transversal points (Fig. 5). So "mu*nu nonzero and epsilon small" is the whole hypothesis (the amplitude is proportional to mu*nu^(5/2)).
  - Check: `int_0^inf cosh^-7 = 5 pi/32`, so the coefficient is `(15 pi/8) mu nu^(5/2)`, about 5.89 mu nu^(5/2) (my arithmetic).
- **Theorem 2 (p.65).** If `mu*nu` nonzero, in any neighbourhood of epsilon = 0 in R+ there is an interval of epsilon such that
  `z = 0` meets an **uncountable number of invariant curves** of `P_epsilon`, each in finitely many points. In each such
  interval there are at least two epsilon for which `z = 0` contains a pair `a, P_epsilon(a)` on the same invariant curve.
  - Reading: for some large Jacobi constants the circle `z = 0` (the EC orbits' footprint) is cut by invariant tori of the
    Levi-Civita 3-sphere, so the original problem has **invariant punctured tori that can contain EC orbits** (Fig. 8).
  - The last claim comes from a topological lemma (p.70) on two crossing closed curves and a family of homeomorphisms.

## 3. The KAM step and its numbers (pp.69-70)

- Truncated normal form of `P_epsilon` on `A` in coordinates `(phi, rho)` in `(R/Z) x R`, boundaries about `|rho| = 1`,
  computed by Chenciner in a course based on Conley's thesis (not in print):
  `P_epsilon(phi, rho) = (phi + 1/2 - nu epsilon^3/2 - (3/2)(1 - mu/4) nu^2 epsilon^6 rho + O(epsilon^7), rho + O(epsilon^7))`.
  The circle `z = 0` is the graph of a function `zeta` with `||zeta||_0 = O(epsilon^3)`. The paper says it did the calculation only to this order.
- Rotation number `omega = -1/2 + (nu/2) epsilon^3`, twist `gamma = (3/2)(1 - mu/4) nu^2 epsilon^6` (after `phi -> -phi`).
  Pick epsilon_1 with `|(nu/2) epsilon_1^3 - 1/2 - p/q| >= [(3/2)(1 - mu/4) nu^2 epsilon_1^6] / q^(2 + beta)` for all p/q
  (an uncountable set). On the epsilon window `[epsilon_1', epsilon_1'']` (length `O(epsilon_1^4)`) the invariant curve of the
  Siegel-Moser theorem (analytic, intersection property) sweeps from `Sigma >= 1/2` to `Sigma <= -1/2`. It must cross `z = 0`.
  - The twist is O(epsilon^6) and the perturbation O(epsilon^7), which is why the theorem applies for small epsilon.
  - The intersection property holds because `P_epsilon` preserves a measure positive on open sets.
- Arithmetic only (my computation from the printed formulas, not in the paper): `omega` is -0.4960 at epsilon = 0.2 and
  -0.4995 at 0.1 (mu = 0); `gamma` is 9.6e-5 and 1.5e-6. At the Earth-Moon mass ratio mu = 0.01215 the values change by under 3%.
  These show the tori are very thin in rotation number, a reason numerical detection needs high precision.

## 4. What it gives the project (#899)

- **Existence of EC orbits for all mu** at large C, with exactly four transversal ones for epsilon small. This sharpens Llibre 1982
  (at least two, small mu) and Lacomba-Llibre 1988 (transversal, small mu), both digested.
- **KAM tori through the EC circle.** The only held-digest source for the statement "invariant tori in the energy surface contain
  or meet EC orbits". The bound on C is unquantified, so it cannot support any specific Jacobi constant.
- **Transversality.** Theorem 1 gives it; it is also used in the Rodriguez del Rio thesis (Theorem 1 there quotes it).
- **Caution on the two constants:** "four" in the thesis equals "eight points / twofold cover" here. A project claim of
  "eight EC orbits" would be wrong by the cover.

## 5. Citation mining

- Lacomba and Llibre 1988 (ref [5], "to be published"): this batch, wanted-list row 36.
- Conley 1963, "On some new long periodic solutions of the plane restricted three body problem", CPAM 16:449-467: not held,
  not on the wanted list. This is the setting for the annulus and the Levi-Civita treatment. Proposed addition.
- Birkhoff 1935 (Ann. Sc. Norm. Super. Pisa 4:267, first memoir on the restricted problem): not held (the held Birkhoff file is
  the 1915 Palermo paper). Not on the list.
- Chenciner, "Le probleme de la lune et la theorie des systemes dynamiques", Paris VII preprint: not held.
- Devaney 1981 "Singularities in classical mechanical systems": not held (the held Devaney file is CMP 80:465). Row 36.
- McGehee 1978 ICM talk: not held. Moser 1970 (CPAM 23:609, regularisation of Kepler): not held.
- Sanders 1982 (Celest. Mech. 28:171, Melnikov and averaging): not held. Siegel and Moser 1971 (Lectures on Celestial
  Mechanics, ref [9], source of the invariant-curve theorem): not held.
- Related held works: Llibre 1982 (HELD), Llibre-Martinez-Alfaro 1985 (HELD), Rodriguez del Rio 2021 thesis (HELD).

*Filed as `cyclers_pdf/papers/chenciner-llibre-1988-invariant-punctured-tori-planar-circular-restricted-three-body-etds-8-63-doi-10.1017-s0143385700009330.pdf`. Check scripts, outputs and other files named above are filed beside it as `cyclers_pdf/papers/<pdf stem>-<file name>`.*

*Wanted-list row numbers in this digest are the batch-30 numbering; the list was renumbered after batch 34.*
