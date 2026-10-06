# Digest: O'Brien 1972, "Almost Periodic and Quasi-Periodic Solutions of Differential Equations" (#960 batch 26)

G. C. O'Brien, PhD thesis, Australian National University, Canberra, 1972 (supervisor W. A. Coppel).
- Filed as `cyclers_pdf/papers/obrien-1972-almost-periodic-quasi-periodic-solutions-differential-equations-phd-thesis-anu.pdf`.
  152 pages, text layer (about 31.7k words), md5 39e2ef1b7c598ed1ddf3f38b29e9bfb0. Supplied by the owner.
- I read Chapter I (the summary of the whole thesis) and the statement of Theorem 3.6.1, on the page
  image (thesis p.39), from the text layer. Chapters II, IV and V were skimmed through the Chapter I
  summary only. A short digest, as asked; general background, not a wanted-list row.

## 0. Verdict

Pure ODE theory: quasi-periodic existence and structure theory. **Only Chapter III touches Hamiltonian
systems, and none of it is celestial mechanics.** There are no restricted-problem examples and no
numbers to use. It is useful only as a citable statement of the frequency-count bound below.

- **Theorem 3.6.1 (p.39, read on the page image).** Let phi(t) be a quasi-periodic solution of the
  conservative Hamiltonian system J x' = dH/dx, with x in R^(2n) and H having a continuously
  differentiable derivative. Suppose its spatial extension Phi(t_1, ..., t_k) is 2 pi/beta_i-periodic
  and C^1 in each t_i, with beta_1..beta_k linearly independent over the integers. **Then the frequency
  basis has at most n elements: k <= n.**
  - This proves the second half of Cherry's claimed theorem (a holomorphic Hamiltonian system has at most
    n basic frequencies) under much weaker hypotheses.
  - O'Brien gives partial counter-examples that cast doubt on the first half: that a holomorphic
    conservative system on a compact invariant set has only quasi-periodic solutions. One example is a
    non-holomorphic system in R^3 with a positive integral invariant and an almost periodic but not
    quasi-periodic solution.
- **Relevance to `quasi_cycler` and the torus routes:**
  - The theorem gives a sourced upper bound. A quasi-periodic orbit of the planar CR3BP (n = 2) has at
    most 2 independent frequencies; the spatial CR3BP (n = 3) has at most 3.
  - A candidate torus whose fitted frequency set has more independent frequencies than the degrees of
    freedom is a fitting artefact, not a torus.
  - That is the whole use. Existence of tori (KAM) is not proved here.

## 1. Content (READ from Chapter I)

- **Ch. II:** the frequency basis of almost periodic solutions of almost periodic ODEs x' = psi(x, t)
  in R^p. There can be at most a finite number (< p) of frequencies beyond those of psi's frequency basis
  (Cartwright's result, with a shorter proof that also builds a family of related almost periodic
  solutions).
- **Ch. III:** quasi-periodic solutions of Hamiltonian systems and Cherry's theorem (above).
- **Ch. IV:** reducibility of holomorphic quasi-periodic linear systems x' = A x + P(phi) x, with real
  distinct eigenvalues of A and a Diophantine frequency vector |(k, omega)| > gamma |k|^-tau.
  - This is a Floquet theory obtained by accelerated (Newton) convergence, after Mitropolskii &
    Samoilenko, with explicit constants.
  - It uses the small-divisor sum estimates of Siegel, Moser and Russmann.
- **Ch. V:** the same for differentiable (not holomorphic) P, using the Moser and Russmann holomorphic
  approximation. Fewer derivatives are lost than in Mitropolskii & Samoilenko.
- **Bibliography:** three parts (general; almost periodic solutions; accelerated convergence),
  pp.108-152. It is a 1972 survey of the KAM and Moser literature. No celestial-mechanics items were
  mined.
