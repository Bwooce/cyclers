# Digest: Gómez & Llibre 1981, "A Note on a Conjecture of Poincaré" (#960, wanted-list row 31)

G. Gómez and J. Llibre (Universitat Autònoma de Barcelona), Celestial Mechanics 24:335-343 (1981),
doi 10.1007/BF01230393, ADS 1981CeMec..24..335G (stamped on every page image). Received 5 December 1979, accepted 7 April 1980.
- Proposed filename: `cyclers_pdf/papers/gomez-llibre-1981-note-conjecture-poincare-density-periodic-orbits-restricted-three-body-celest-mech-24-335-doi-10.1007-BF01230393-ads-scan.pdf`.
  - The file to file is the `ocrmypdf --force-ocr -O 3 -l eng` copy `1981CeMec24335G-ocr.pdf`, md5
    61a510c4cf8051b110ad01c36b9212d3 (506,507 bytes).
  - Source: image-only ADS JBIG2 scan `1981CeMec24335G.pdf`, 9 pp., md5 96882002a4bbd8b54e974b2c34677096 (130,356 bytes),
    from https://articles.adsabs.harvard.edu/pdf/1981CeMec..24..335G (first try).
  - Dark-pixel fractions at 100 dpi, original against OCR copy: p.1 0.0527/0.0558, p.5 0.0478/0.0506, p.9 0.0382/0.0402.
  - PDF page n = journal page 334 + n.
- **How it was read:** pp.335-337 and 339 on 110 dpi page images (theorem, eqs. 2-3, dr/da, q(a), Table I); pp.338,
  340-343 (twist map, the proof of Theorem 2, references) on the OCR text layer.
- No numeric table, so no YAML. Check script: `checks/check_gomez_llibre.py` / `.out` (ours).

## 0. Verdict

**A short proof of a weak form of Poincaré's conjecture that periodic orbits are dense in the planar circular restricted
problem.** Theorem (p.335): for any fixed Jacobi constant C and any epsilon > 0 there is mu_0 > 0 such that for
mu in [0, mu_0] the bounded orbits NOT in the closure of the set of periodic orbits have Lebesgue measure < epsilon.
- Method: at mu = 0 (rotating two-body problem) the pericentre section r' = 0, r'' > 0 is a union of rings on which the
  return map is a twist map r' = r, phi' = phi + f(r), f(r) = -2 pi a^(3/2) (eq. 4). For small mu the map is an
  area-preserving perturbation; KAM (Moser-Rüssmann, Arnold) gives invariant curves leaving out small measure (Theorem 1),
  and periodic points accumulate on every invariant curve (Theorem 2, proved in an appendix, Lemmas 1-3).
- It is NOT a collision-manifold paper. Wanted-list row 31 groups it with "collision-manifold and density theorems"; this
  paper supplies only the density half. The collision-manifold half of row 31 is Devaney 1981 (still open).

**For the project (PROPOSALS only):**
- It closes the Gómez-Llibre part of wanted-list row 31 (line 89). Devaney 1981 stays open in that row.
- It is the predecessor that Llibre 1982 Theorem D improves (held digest
  `2026-10-05-digest-llibre-1982-restricted-problem-small-mu.md`, sec. 6 and line 257, which marks this paper "not held").
  The two use the same mu = 0 ring tables and the same twist map: Llibre 1982 eq. 4.2 is this paper's eq. 3, and its
  f(r) = -2 pi a^(3/2) is this paper's eq. 4. Llibre 1982 sharpens "periodic orbits" to ten classes C_ij of SYMMETRIC
  periodic orbits. The 1981 statement is the weaker one.
- Relation to the held Lacomba-Llibre 1988 digest (`2026-10-07-digest-lacomba-llibre-1988-transversal-ejection-collision-orbits.md`):
  little overlap. Lacomba-Llibre is about ejection-collision orbits and the collision manifold (a torus); this paper never
  uses the collision manifold. One shared fact: this paper's collision ellipse z_c = 1/C (semimajor axis 1/C, p.338) has
  apocentre 2/C, which is the mu = 0 ejection-collision circle r = 2/C of Lacomba-Llibre Prop. A.2 (our arithmetic,
  apocentre = 2a for a rectilinear ellipse). Lacomba-Llibre does not cite this paper (no "Gomez" hit in that digest).
- `#899`: background only. It gives no orbits, no numbers and no method the project can run. Low value beyond completing
  the Llibre 1982 citation chain.
- No catalogue change.

## 1. Content

- **Sec. 1** (p.335): RTBP in the synodic frame with primaries at (mu, 0) and (mu - 1, 0) (Szebehely's convention, [8]);
  C = 2 Omega - (x'^2 + y'^2). Poincaré's 1892 conjecture (quoted from Méthodes nouvelles vol. 1 p.82) and
  Schwarzschild's 1896 version are stated.
- **Sec. 2, the two-body rotating problem** (pp.336-339). Bounded orbits are Kepler ellipses; for fixed C,
  +-b = (C/2) a^(1/2) - (1/2) a^(-1/2) (eq. 2; + direct, - retrograde). Circular orbits a = +-b meet at a1, a2, a3 for C > 3,
  at a1 and a2 = a3 = 1 for C = 3, at a1 only for C < 3. On the pericentre section the ellipse (a, b) is the circle
  r = a(1 - e) = a{1 - [1 - (1 - Ca)^2/(4a^3)]^(1/2)} (eq. 3). Zeros of dr/da come from q(a) = (1 - Ca)(1 + Ca)^2 + 32 a^3:
  one negative root for C < 3, plus two positive roots z1, z2 for 3 <= C < 32^(1/3), one positive root for C > 32^(1/3).
  z_c = 1/C is the collision ellipse. Parabolic orbits map to r = C^2/8. Table I (p.339) lists, by C range and direction,
  the inner and outer radii of the rings of elliptic orbits (bounded by parabolic, collision, circular or z2 orbits).
- **Sec. 3, proof** (pp.339-340). At mu = 0 an ellipse is periodic in the rotating frame iff a = q^(2/3), q rational,
  so periodic orbits are dense ("see [2] p.709", Birkhoff). For small mu, on rings R inside the Table I rings, T_mu is
  r' = r + mu g1, phi' = phi + f(r) + mu g2, area-preserving (the polar synodic coordinates are canonical positions,
  Szebehely p.354), |df/dr| > 0, C^5 norm small, so the Moser-Rüssmann theorem applies.
- **Appendix** (pp.340-343): a proof of Theorem 2 ("well known", Moser 1973 p.54, but the authors had "never seen a proof").
  Acknowledgement: C. Simó, C. Marchal, and C. Perelló for the idea of the proof.

## 2. Checks (`checks/check_gomez_llibre.out`, ours)

- Eq. 2 against C = 1/a +- 2 [a(1 - e^2)]^(1/2) over 1000 random ellipses: max residual 1.8e-15.
- Eq. 3 equals a(1 - e) with b from eq. 2 (max 5.3e-15); the printed dr/da agrees with a finite difference of eq. 3
  (max relative difference 2.4e-9). r -> C^2/8 as a -> infinity (C = 2: 0.4999999970 against 0.5).
- Circular direct orbits (2 a^(3/2) - C a + 1 = 0): no positive root at C = 2.9; double root a = 1 at C = 3; two roots
  (a = 0.7023, 1.4557) at C = 3.1. As stated.
- q(a) roots: C = 2.5 one negative and two complex; C = 3.0 roots 1, 1, -0.2; C = 3.1 two positive (3.968, 0.579) and one
  negative; at C = 32^(1/3) = 3.174802 the cubic term vanishes (q is quadratic, one positive root 0.5096); C = 3.3 one
  positive, two negative. All as stated on p.337.

## 3. Citation mining

References (p.343, text layer):

| Reference | Held? | Wanted list |
|---|---|---|
| [1] Arnold 1963, Russian Math. Surveys 18:9-36 (KAM proof) | NOT held (the only `arnold` file is Capinski et al. on Arnold diffusion) | no row. Not needed. |
| [2] Birkhoff 1923, Collected Mathematical Papers vol. 1 (cited at pp.706 and 709) | The Collected Papers are NOT held. Birkhoff 1915, Rend. Circ. Mat. Palermo 39:265-334, is HELD and digested; we did not check whether the cited pages reprint it. | no row |
| [3] McGehee 1969, thesis (Univ. Wisconsin), homoclinic orbits | NOT held (index has McGehee only as citations) | no row for the thesis. Low priority. |
| [4] Moser 1962, Nachr. Akad. Wiss. Göttingen (invariant curves) | NOT held | no row (row 51 is Moser 1953, a different paper). Not needed. |
| [5] Moser 1973, Stable and Random Motions (Princeton) | NOT held | no row. Book. Low priority. |
| [6] Poincaré 1892, Méthodes nouvelles vol. 1 | NOT held | no row. |
| [7] Schwarzschild 1896, Astron. Nachr. 147:17-24 | NOT held (no file, not in index) | no row. Historical; low priority. |
| [8] Szebehely 1967, Theory of Orbits | HELD (`szebehely-1967-theory-of-orbits-...-book.pdf`) | n/a |

- No new must-have lead. Wanted-list row 31 after this paper: Devaney 1981 ("Singularities in classical mechanical
  systems", doi 10.1007/978-1-4899-6696-4_7) remains.

*Filed as `cyclers_pdf/papers/gomez-llibre-1981-note-conjecture-poincare-celest-mech-24-335-doi-10.1007-BF01230393-ads-1981CeMec-24-335G.pdf`.*
