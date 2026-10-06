# Digest: Pfenniger 1985b, "Numerical study of complex instability. II. Barred galaxy bulges" (#960 batch 30)

D. Pfenniger (Observatoire de Geneve), Astronomy & Astrophysics 150:112-128 (1985), ADS 1985A&A...150..112P.
Received 21 January 1985, accepted 21 March 1985.
- Given file: supplied file `pfenniger_150_112-ocr.pdf`, 17 pp., md5 a45ff8d22c922efbb058c6afffb5487b (force-OCR copy).
- **Identity check:** pp. 16-32 of supplied file `pf97.pdf` (md5 3e2d9f1099e5951f7fc27ed13a2b4849) are this same paper.
  I compared page images of pf97.pdf pp. 16 (title page: "II. Barred galaxy bulges", A&A 150, 112-128) and 32
  (reference list, journal p. 128, ending "Pfenniger 1985a, A&A 150, 97 (Paper I)"), and pp. 23 (journal p. 119, Fig. 5),
  against the OCR text of pages 1, 8 and 17 of the given file. Page counts agree (32 = 15 + 17). File it once.
  Pages 1-15 of pf97.pdf are Paper I (separate digest).
- Filed as `cyclers_pdf/papers/pfenniger-1985b-numerical-study-complex-instability-II-barred-galaxy-bulges-aa-150-112-ads-1985AA-150-112P.pdf`.
- The first OCR copy (`ocrmypdf --redo-ocr`) lost the page images of this ADS JBIG2 scan. The filed copy was re-made with `ocrmypdf --force-ocr` (md5 a45ff8d22c922efbb058c6afffb5487b); every page was checked to render.
- How I read it: the full text layer; page images for pp. 112 (title, summary), 119 (Fig. 5 and the search method) and 128.
  Numbers below come from the text layer except where marked "image".
- Wanted-list row 38 (with Paper I). Follow-up of the Jorba-Olle 2004 digest.

## 0. Verdict

**An application of Paper I to a 3D rotating bar model of a galaxy. It shows the Hopf-like bifurcation at the end of a
stable periodic family in a real potential, and it is a worked search for the 3-periodic orbits. Galactic dynamics;
no direct use for planetary cyclers.**
- The main result: the family B_z2 (4/4/1 resonance, the box-shaped bulge orbit) turns complex unstable. At the
  transition a Hopf-like bifurcation of 2D tori occurs, as in the maps of Paper I. The unstable side feeds "semi-ergodic"
  orbits: long-lived chaotic orbits that diffuse in angular momentum and height (Section 4).
- Catalogue: nothing to add (no hit for "Pfenniger" or "complex instab" in `data/catalogue.yaml`).
- **Use for the project:** a method reference. Section 3 gives a practical recipe to find the multi-periodic orbits that
  accompany a complex-unstable transition in a Hamiltonian flow (below). That may help if a spatial CR3BP family is ever
  found to lose stability by Delta < 0.

## 1. Model and family (READ section 2)

- Rotating barred-galaxy potential with a closed-form bar (parameters GM_B, bar minor semi-axis c, bar major semi-axis a,
  corotation at the bar end). Strong bar: c = 0.6, GM_B = 0.1. Weak bar: GM_B = 0.01. (Text layer.)
- B_z2 is the stable 3D 4/4/1 resonant family. Its stability ends abruptly on the retrograde side of the characteristics,
  starting already at infinitesimal bar perturbation. This is where complex instability arises (Fig. 3).

## 2. The transition (READ section 3)

- In Fig. 3 the transition occurs where the stability indices b_i are near +1. By Paper I (b_i = -2 cos 2 pi k) this is
  k = 1/3 or 2/3, a 3-periodic resonance. So B_z2 sits near a 3-periodic bifurcation.
- The bifurcation is a two-parameter event. To follow a family one must vary two parameters, here the energy H and the
  bar minor semi-axis c. At fixed potential, only H varies, and the bifurcating object is then a family of invariant
  curves in the 4D section, not a family of periodic orbits.
- Two 3-periodic bifurcations were found: one near c = 0.9 on the stable side of B_z2 (unstable orbits, "reverse"
  Hopf-like), and one near GM_B = 0.09 on the unstable side (stable and semi-unstable orbits, "normal"). As in Paper I,
  there are 2 x 2 distinct 3-periodic orbits at once (even-semi-unstable and even-even unstable types in the reverse case).
- Image (journal p. 119, Fig. 5 caption): an unstable 3-periodic orbit with H = -0.190, x0 = 0.017230, z0 = 0.749163,
  xdot0 = 0.00261170, zdot0 = 0.00462144, c = 0.8878, GM_B = 0.1. It winds helically round the B_z2 orbit and closes after
  three turns. I read these on the image. I have not re-integrated it (the potential's closed form is in "Paper 0", not held).

## 3. Search method (READ section 3.2, image p. 119)

- Newton on the extended map, extended with model parameters: solve (A - I) dx = Delta x in the least-squares sense
  with the HFTI algorithm (Lawson and Hanson 1974), since A - I is near-singular when det A = 1 and an eigenvalue is
  near 1. The Jacobian of the extended map is built from two extra orbits with variation 1e-7. (Eqs. 5-6, image.)
- Close to the bifurcation the four eigenvalues of the bifurcating families are all near +1, and 12 periodic points
  lie close together. The paper says this explains why these families were missed in Contopoulos and Magnenat 1985.
- Invariant curves: contract the section momenta at each crossing (a volume-contracting "dissipation" that keeps H),
  then the stable tori act as attractors (Fig. 9). Confinement of consequents for more than 500 turns marks the tori (Fig. 7).

## 4. Semi-ergodic diffusion (READ section 4, text layer)

- Orbits started on the complex-unstable part of B_z2 diffuse over long times and fill the inner
  halo in a roughly exponential surface density, 2.5 times steeper than the model's (Fig. 12 caption). The paper suggests
  this populates the inner halo after bar formation. This is stellar dynamics, not relevant to cyclers.

## 5. Implications for monodromy classification

- Same as Paper I section 3. The new point: the transition was met in a continuous family of a physical model, so complex
  instability is not an artefact of the toy maps. For a CR3BP spatial family, expect it at crossings where the (alpha, beta)
  point meets the parabolic arc. The paper shows the closest rational resonance sets what bifurcates.
- A 3-periodic bifurcation from a stable orbit is not "period tripling" in the usual one-dimensional sense. Four
  distinct families appear, not two.

## 6. Citation-mining

- Not held: Pfenniger 1984a,b (Paper 0; A&A 134:373 and 141:171); Magnenat 1982a,b; Contopoulos and Magnenat 1985; Combes
  and Sanders 1981; Kormendy and Illingworth 1982, 1983; Sellwood 1980, 1981; Henon 1983 (Les Houches); Lawson and Hanson 1974.
  None is on the wanted list (row 38 names only Pfenniger 1985, Olle and Olle-Pfenniger).
- HELD: Broucke 1969 and 1969b (as in the Paper I digest). No other cited work is a cycler or restricted-problem
  periodic-orbit paper in the corpus.
- Olle 2000 and Olle-Pfenniger 1999 (row 38) are not cited here; they come from the Jorba-Olle digest.

*Wanted-list row numbers in this digest are the batch-29 numbering; the list was renumbered in batch 30.*
