# Digest: Pfenniger 1985a, "Numerical study of complex instability. I. Mappings" (#960 batch 30)

D. Pfenniger (Observatoire de Geneve), Astronomy & Astrophysics 150:97-111 (1985), ADS 1985A&A...150...97P.
Received 2 January 1985, accepted 21 March 1985.
- Given file: supplied file `pf97-I-ocr.pdf`, 15 pp., md5 b9f29c97ac30066349fd9187b16a0b00 (force-OCR copy). It is pages 1-15 of the 32-page scan
  supplied file `pf97.pdf` (md5 3e2d9f1099e5951f7fc27ed13a2b4849; image-only JBIG2 scan, ADS stamp on every page).
- Filed as `cyclers_pdf/papers/pfenniger-1985a-numerical-study-complex-instability-I-mappings-aa-150-97-ads-1985AA-150-97P.pdf`.
- The first OCR copy (`ocrmypdf --redo-ocr`) lost the page images of this ADS JBIG2 scan. The filed copy was re-made with `ocrmypdf --force-ocr` (md5 b9f29c97ac30066349fd9187b16a0b00); every page was checked to render.
- How I read it: the full text, then page images.
  - Image-checked (journal pp. 98, 99, 101, 103): Eqs. 3-9, Fig. 1, Eqs. 12-23, the section 4.2.1-4.2.3 text and Fig. 5 caption.
  - Text layer only: sections 4.2.4-4.2.6, the captions of Figs. 10-16, and the conclusions.
- Wanted-list row 38 (Pfenniger 1985, A&A 150:97 and 150:112; Olle 2000 and Olle-Pfenniger 1999 stay on the row).
  Follow-up of the held Jorba-Olle 2004 digest (`2026-10-05-digest-jorba-olle-2004-invariant-curves-hamiltonian-hopf.md`,
  lines 30 and 107).

## 0. Verdict

**This is the origin of the 4D maps T_s and T_t that Jorba-Olle 2004 describe. It gives a clear numerical account of
complex instability and of what happens at its onset. It has no catalogue content.**
- It defines the transition of a periodic orbit (3 degrees of freedom, or any 4D symplectic return map) from stable to
  complex unstable: two equal eigenvalue pairs on the unit circle leave it together.
- It finds that the transition behaves like a Hopf bifurcation in a conservative system: families of isolated invariant
  curves (2D tori in phase space) appear. Paper II finds the same in a galactic potential.
- Catalogue: no hit for "Pfenniger", "complex instab" or "complex unstable" in `data/catalogue.yaml` or `src`.
  **Proposal only:** none for the catalogue. If the stability classifier is ever extended to spatial orbits, cite this
  paper with the held Broucke 1969b for the (alpha, beta) plane.

## 1. Definitions (READ pp. 98-99, image)

- Return map T of the 4D section; A = Jacobian at the fixed point, det A = 1, symplectic. Eigenvalues come in pairs
  (lambda, 1/lambda). Two pairs give four cases: (1) one quadruplet lambda, 1/lambda, conj(lambda), 1/conj(lambda) =
  complex instability; (2) two pairs on the unit circle = stability; (3) one circle pair and one real pair =
  semi-instability; (4) two real pairs = double instability.
- Eq. 3: lambda^4 + alpha lambda^3 + beta lambda^2 + alpha lambda + 1 = (lambda^2 + b1 lambda + 1)(lambda^2 + b2 lambda + 1).
  Eq. 4: b1 = -(lambda + 1/lambda), b2 = -(mu + 1/mu), Delta = (b1 - b2)^2. Eq. 5: alpha = b1 + b2, beta = 2 + b1 b2.
- Stable if |b1| <= 2, |b2| <= 2 and Delta >= 0. **Complex instability is Delta < 0**, where b1, b2 are complex.
  The (alpha, beta) plane is the Broucke diagram (Fig. 1): seven regions (stable; complex; even-odd, even-even, odd-odd;
  even and odd semi-instability). Only a stable orbit can reach any other region.
- Check: Eq. 5 gives Delta = alpha^2 - 4 beta + 8 (Eq. 15 prints the same). Eq. 8, the parabolic arc where the
  transition happens, is alpha = -4 cos(2 pi k), beta = 2 + 4 cos^2(2 pi k). Then Delta = 16 cos^2 - 8 - 16 cos^2 + 8 = 0.
  The eigenvalue arguments are +-2 pi k.

## 2. The transition (READ pp. 99-100, 103-109)

- On the arc, the two pairs are equal. A is then either diagonalisable (two uncoupled 2D systems; the eigenvalues cannot
  leave the circle) or has two 2x2 Jordan blocks (Eq. 9). Only the Jordan case lets them leave. My gloss: this is the
  Krein collision of two pairs. The paper cites Krein and Jakubovic but does not develop signatures.
- If k = m/n is rational, n-periodic orbits may bifurcate (points P_n on Fig. 1). If k is irrational, no periodic orbit
  bifurcates at that point, but n-periodic bifurcations with large n lie arbitrarily near it on the arc.
- The bifurcating object is a family of isolated invariant curves in the section. It lies on the unstable side of the
  central orbit and is stable ("normal" Hopf), or on the stable side and is unstable ("reverse" Hopf).
- Section 5: the stable curves are separated by multi-periodic orbits; the unstable ones make thin stochastic layers
  that fatten with the perturbation. The paper has no algorithm for finding unstable invariant curves.

## 3. What it implies for monodromy classification (my application)

- A spatial CR3BP orbit has a 6x6 monodromy matrix with the trivial double eigenvalue 1. The remaining 4x4 block has the
  structure above, so the Broucke (alpha, beta) test is the right classifier for spatial orbits.
- A planar orbit has one non-trivial pair only, so complex instability cannot occur. A stability index that tests only
  |b| <= 2 and omits Delta >= 0 can call a complex-unstable spatial orbit stable.
- Eigenvalues on the circle do not prove stability when they are equal (Krein 1950, cited): a double eigenvalue needs the
  Krein signature or a local test.
- The paper's own examples show 12 periodic branches (n = 3) leaving the central orbit, so a continuation tool near the
  transition meets many close periodic orbits.

## 4. The two maps (READ p. 101, image)

- T_s (Eq. 12): x1' = D[x1 + K1 sin(x1+x2) + L1 sin(x1+x2+x3+x4)], x2' = x1 + x2, x3' = D[x3 + K2 sin(x3+x4) +
  L2 sin(x1+x2+x3+x4)], x4' = x3 + x4, mod 2 pi. T_t (Eq. 13) uses tan in the L term. D < 1, = 1, > 1 contracts,
  preserves or dilates volume (used only as a numerical device). Froeschle's map is T_s with L1 = L2, D = 1.
- Eq. 15: alpha = -(4 + K1 + L1 + K2 + L2), Delta = (K1 - K2 + L1 - L2)^2 + 4 L1 L2. So Delta < 0 needs L1, L2 of opposite
  sign, and Froeschle's map cannot reach complex instability. The paper sets L1 = -L2 = L, K2 = 0, so Delta = K(K + 4L)
  and L_crit = -K/4 for -8 < K < 0.
- Script `check_pfenniger_map.py` (output `check_pfenniger_map.out`): the Jacobian of Eq. 14 has det 1 and unit-modulus
  eigenvalues below L_crit, and a quadruplet above it (moduli 0.8475 and 1.1799 at K = -1, L = 0.3). L_crit = -K/4
  holds for K = -1 (0.25), -4 (1.00), -6 (1.50) and -1.171573 (0.292893). The paper prints L_crit = 1.50 (Fig. 5
  caption, image) and 1.00, 0.292893 (Fig. 11-12 captions, text layer only).

## 5. Numerical findings (READ section 4)

- Newton search with least squares (HFTI) for n = 3 to 8 periodic orbits, 16-digit arithmetic on an Apple II and an HP-1000F.
- n = 3 (K = -6; image): 12 (2 x 2 x 3) branches leave the origin "like umbrella ribs", four distinct families. For T_s the
  stable branches become odd semi-unstable at L = 1.668 (period doubling), then stable again, then complex unstable near
  L = 1.749. The two semi-unstable branches become even-odd unstable near L = 1.628. n = 5, 7 behave like n = 3.
- n = 4 is a "strong resonance": two families, eight branches. The T_s families reach a maximum at L = 1.397 (image).
  n = 6, 8 have two families.
- Near the complex-unstable orbit, consequents stay confined for a long time (Magnenat 1982a,c saw this). For T_t the
  bifurcation is on the stable side and escape is fast. Dissipation (D slightly different from 1) catches stable curves.

## 6. Citation-mining

- Not held, not on the wanted list: Pfenniger 1984 (A&A 134:373, "Paper 0"); Magnenat 1982a,b,c; Contopoulos and Magnenat
  1985; Froeschle 1971, 1972; Gel'fand and Lidskii 1955; Krein 1950; Krein and Jakubovic 1983; Lichtenberg and Lieberman
  1983; Guckenheimer 1984; Henon 1966 and Martinet and Mayer 1975 (2D n = 3 bifurcation examples); Hadjidemetriou 1985.
  New candidates that matter: Magnenat 1982a,c and Contopoulos-Magnenat 1985 (3D complex-unstable orbits), Hadjidemetriou 1985.
- HELD: Broucke 1969 (`broucke-1969-periodic-orbits-elliptic-restricted-three-body-problem-jpl-tr-32-1360...`) and Broucke
  1969b (`broucke-1969b-stability-periodic-orbits-elliptic-restricted-three-body-problem-aiaa-j-7-1003...`);
  Hadjidemetriou 1975 (`hadjidemetriou-1975b-stability-periodic-orbits-three-body-problem-celest-mech-12-255...`).
- Poincare 1892 is not held and not listed.

*Check scripts and outputs named above are filed beside the PDF as `cyclers_pdf/papers/<pdf stem>-<script name>`.*

*Wanted-list row numbers in this digest are the batch-29 numbering; the list was renumbered in batch 30.*
