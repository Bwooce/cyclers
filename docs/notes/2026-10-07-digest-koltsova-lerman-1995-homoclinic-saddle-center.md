# Digest: Koltsova and Lerman 1995, "Periodic and Homoclinic Orbits in a Two-Parameter Unfolding of a Hamiltonian System with a Homoclinic Orbit to a Saddle-Center" (#960 batch 32)

O. Yu. Koltsova (Volga State Academy of Water Transport) and L. M. Lerman (Research Institute for Applied Mathematics and
Cybernetics, Nizhny Novgorod), International Journal of Bifurcation and Chaos 5(2):397-408 (1995), received 4 February 1994.
DOI 10.1142/s0218127495000338 (given in the task brief; the DOI is not printed on the page, the title and pages are). 12 pp.
- Source scan: `837710e0-koltsova1995.pdf`, md5 3ba1922e34773c291c3bb7c75dfa07e2. It is image-only (no text layer).
- OCR copy made here: `ocrmypdf --force-ocr -l eng` (not --redo-ocr), file `koltsova1995-ocr.pdf`, md5
  99d5b710c07b3b99e7d3805903d4827e, 12 pp. Page 2 of the OCR copy was rendered and is not blank (full text page, p.398).
  The text layer has the columns interleaved and loses some formulas; I used it for orientation only.
- Proposed corpus filename (file the OCR copy):
  `koltsova-lerman-1995-periodic-homoclinic-orbits-two-parameter-unfolding-hamiltonian-homoclinic-saddle-center-ijbc-5-397-doi-10.1142-S0218127495000338.pdf`
- How I read it: I read all 12 page images (pp.397-408) at 110 dpi. Every theorem and proposition below is stated from the image.
  I did not rederive the proofs. No number in the paper is a computed orbit, so there was nothing to integrate.

## 0. Verdict

**Pure local-to-global Hamiltonian theory. It contains no numbers and no restricted-problem model.** It matters to the project
only as the source of the "four-fold" structure of homoclinic orbits to a saddle-center (the L1 or L2 type point), which the
held Barrabes-Mondelo-Olle 2009b digest cites as agreeing with it.
- The setting is any real-analytic two-degree-of-freedom Hamiltonian with a saddle-center `p` and a homoclinic loop `Gamma`.
  The collinear points L1 and L2 of the planar CR3BP are saddle-centers, so the theory applies there. The paper says so
  (p.397: collinear libration point, doubly asymptotic orbits of Deprit-Henrard and Lidov-Vashkov'yak).
- It does not treat ejection-collision orbits or invariant tori of the CR3BP. The only CR3BP citation is Llibre,
  Martinez and Simo 1985 (Poincare homoclinic orbits near L2), which we hold.
- **Proposal only:** no catalogue row. Use as background for #899 and for reading BMO 2009b. Do not cite it for any value.

## 1. Setting (pp.398-400)

- Eigenvalues of the saddle-center are +-i*omega and +-lambda. A Lyapunov family of periodic orbits leaves `p`.
- Normal form (Moser 1958, Russmann 1964): `H_mu = h(xi, eta) = lambda*xi + omega*eta + ...`, with `xi = x1*y1`,
  `eta = (x2^2 + y2^2)/2`. After time scaling, `lambda = 1` and `gamma = omega/lambda > 0`.
- Two-parameter unfolding `H_mu`, `mu` in R^2, because a loop to a saddle-center has codimension 2 in the space of
  Hamiltonians. (Reversible systems reduce it to codimension 1.)
- Two cases for how the loop re-enters `p`: A (local axis x1 < 0) and B (x1 > 0). **Proposition 1** gives the topology of the
  level `V_h = {H = h}` near `p` union `Gamma`: in case A a torus-plus-cylinder split for h < 0, a pinched union at h = 0,
  and a connected sum for h > 0. In case B there is no returning mechanism for h < 0. Result: study h = 0 and h > 0.
- The global map `S_h` is area preserving; **Lemma 1** shows its derivative at the loop is `diag(b, 1/b)` after rotations.
  The genericity condition is that `det D(F_mu, G_mu)/D(mu1, mu2)` is nonsingular at the origin.
- Local map `T_h` is a rotation by `Delta_h(eta) = -a_h'(eta) ln(d1 d2 / a_h(eta))`. It is discontinuous at the loop.

## 2. Theorems, as printed

- **Proposition 2 / Lemma 2 (pp.400-402).** The linearisation at the loop is a linear Hamiltonian system with time-dependent
  matrix. It is either "linear integrable" (two Lagrangian cylinder foliations coincide) or "linear nonintegrable".
  If the full system is integrable, the linearisation is linear integrable. This is the paper's invariant form of the
  genericity condition (Section 4).
- **Proposition 3 (p.402).** At mu = 0, on the level H = 0 (case A), the return map `f0` on a punctured disk
  `0 < eta <= eta0` has the twist property: `d psi1/d eta = (-gamma + O(sqrt(eta))) / (lambda(phi) * eta) < 0`.
  The twist goes to infinity at the loop. KAM invariant curves are expected only when the system is near integrable.
- **Theorem 1 (p.403).** If the linearisation is linear nonintegrable, there are four countable families of one-circuit periodic
  orbits `p_n^i`, i = 1..4, accumulating to the loop in the level H = 0. With `s = gamma*(b - 1/b)`:
  if |s| > 2 all are hyperbolic; if |s| < 2, two of the four families are elliptic and two hyperbolic, and which pair swaps with the
  sign of s. The proof is a fixed-point count (implicit function theorem on `cos Delta = D1(r, phi)`).
- **Theorem 2 (p.404).** For 0 <= h <= h*, |mu1| <= mu1*, |mu2| <= mu2*, the level `H_mu = h` has a set `P_{h,mu}` of one-circuit
  periodic orbits; the limit as h, mu go to 0 is the set from Theorem 1.
- **Theorem 3 (p.405).** (1) With a linear nonintegrable linearisation there are four countable families of parameter
  values `mu_n^i -> 0` with a two-circuit homoclinic orbit to `p` (homoclinic doubling). (2) If also `b` different from `(1 + sqrt 5)/2`,
  there are eight countable families with a three-circuit homoclinic orbit (tripling). The paper proves the tripling only in
  outline: "Details of this proof will appear elsewhere" (p.406 remark, because of r ln r terms).
- **Hypothesis (p.407, not proved).** An n-circuit homoclinic orbit exists for a sequence `mu_n` for every n.
- **Theorem 4 (p.407).** In the three-parameter space (h, mu1, mu2), h > 0, there is a cone-shaped region where the
  stable and unstable manifolds of the Lyapunov orbit `l_h` meet in four transversal one-circuit homoclinic orbits. Hence
  no analytic additional integral exists for |mu| small (via Cushman 1978).

## 3. What it means for the project

- **Existence, not computation.** It shows that near a saddle-center loop the dynamics is rich but gives no recipe for finding
  loops in the CR3BP. For that the project has numerical continuation (BMO 2009b method; Llibre-Martinez-Simo).
- **Transversality.** Theorem 4 gives transversal homoclinic orbits to the Lyapunov orbit of a saddle-center whose loop
  has nonintegrable linearisation. The CR3BP version at L2 is Llibre-Martinez-Simo 1985 (held). This paper is the general
  statement. It does not give any CR3BP orbit.
- **KAM.** Proposition 3 is the only link to tori: the return map has an infinite twist at the loop, so Aubry-Mather
  theory applies, but the paper says (p.403) that for large |b| numerical simulation shows no invariant curves near the loop.
  So it gives no existence of tori.
- It is not an ejection-collision paper. Do not file it under the #899 ejection-collision count. Use it for the
  homoclinic-doubling background only.

## 4. Citation mining

All "held" checks: `ls cyclers_pdf/papers | grep -i <name>` and `grep -i <name> CORPUS_INDEX.md`.
- Llibre, Martinez and Simo 1985, JDE 58:104-156: HELD (`llibre-martinez-simo-1985-transversality-invariant-manifolds-...`).
- Lerman 1987 (Gorky Univ., Russian; English in Selecta Math. Sov. 10:297, 1991): not held, not on the wanted list. This is the
  source of the four families; Theorem 1 here restates it. A proposal: add to the wanted list.
- Koltsova 1991 (homoclinic doubling, Russian, Nizhny Novgorod): not held. Koltsova 2003 and this paper are wanted-list row 37.
- Deprit and Henrard 1965 (AJ 70:271, doubly asymptotic orbits at L-points): not held, not on the wanted list.
- Lidov and Vashkov'yak 1975 (IPM preprint 115): not held; later Lidov-Vashkov'yak papers are in row 67 only.
- Conley 1969 (JDE 5:136): not held, not on the list. Moser 1958 and Russmann 1964 (normal form): not held.
- Mielke, Holmes and O'Reilly 1992 (JDDE 4:95): not held. Grotta Ragazzo 1993 preprint: not held.
- Aubry and Le Daeron 1983, Mather 1982 (twist-map theory): not held. Cushman 1978, Churchill-Pecelli-Rod 1979: not held.
- Cross-link: the held digest `2026-10-05-digest-barrabes-mondelo-olle-2009b-horseshoe-homoclinic-orbits.md` (line 59) cites
  this paper for its four-fold homoclinic structure.

*Filed as `cyclers_pdf/papers/koltsova-lerman-1995-periodic-homoclinic-orbits-two-parameter-unfolding-hamiltonian-saddle-center-ijbc-5-397-doi-10.1142-s0218127495000338.pdf`.*

*Wanted-list row numbers in this digest are the batch-30 numbering; the list was renumbered after batch 34.*
