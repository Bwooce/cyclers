# Digest: Jefferys 1974, "Stability regions for quasiperiodic motion in the restricted problem of three bodies" (#960 batch 21)

W. H. Jefferys (University of Texas at Austin), Astronomical Journal 79(6):710-721 (June 1974), ADS
1974AJ.....79..710J. Received 4 March 1974.
- Filed as `cyclers_pdf/papers/jefferys-1974-stability-regions-quasiperiodic-motion-restricted-problem-three-bodies-aj-79-710-ads-1974AJ-79-710J.pdf`.
  12 pages, ADS scan with a text layer, md5 3d72752095c53076d6964e9f52bc1221. Supplied by the owner. Not
  a wanted-list row; filed as general background.
- I read secs. I-VII and the references from the text layer. No numbers are taken from figures.

## 0. Verdict

**Background for the quasi_cycler class and the torus routes: a numerical atlas of where quasi-periodic
(KAM-like) motion exists in the planar CR3BP.** There are no periodic-orbit tables and no cycler content.
- **Method:**
  - Surfaces of section in configuration space: points where r-dot = 0, i.e. the pericentres and
    apocentres about the larger primary.
  - Direct apsides are plotted in the upper half-plane and retrograde ones in the lower.
  - Nine mass ratios, many Jacobi constants, over 4000 orbits.
  - The full set of 276 sections is in Jefferys 1971, "An Atlas of Surfaces of Section for the
    Restricted Problem of Three Bodies", Univ. Texas Publ. Dept. Astron. Ser. II, 3, 6; also AMRL 1034.
- **Mass ratios** (sec. VII; text layer, simple values):
  - 0.00095388 (Sun-Jupiter, 1/1048.355);
  - 0.01215067 (Earth-Moon, 1/82.3);
  - 0.09090909 (Darwin's 1/11);
  - 0.2 (Moulton); 1/3; 0.5 (Moulton, Copenhagen); 0.8;
  - 0.98784933 (lunar orbiter); 0.99904612 (Jupiter orbiter).
- **Findings (sec. V):**
  - **Stable retrograde orbits appear at much lower C than stable direct orbits, for every mass ratio.**
    At fixed C and position, a retrograde particle has lower two-body energy, so it is more tightly
    bound.
  - Stable direct orbits appear only once the zero-velocity curve starts to enclose the primary.
  - The resonance islands appear and disappear with C, with a rapid evolution just after the
    zero-velocity curve closes at the inner Lagrange point (Earth-Moon Figs. 20-25).
  - Stable Lagrangian (tadpole) regions appear for the two smallest mass ratios only.
- **Sec. VI:** Hénon-style stability diagrams (Figs. 28-36), one per mass ratio, built from the
  surfaces of section. They agree with Hénon 1965b at mu = 0.5. Note the conventions: Hénon's C is
  smaller by mu(1 - mu), and his x-axis is shifted by mu.

## 1. Relevance

- `quasi_cycler` and torus routes: a qualitative map of the regions where invariant tori exist, per mass
  ratio. It is a cross-check that a candidate quasi-periodic cycler lies in a region where tori are
  numerically observed, not a quantitative criterion.
- It agrees with Hénon 1965b and Hénon 1969 (retrograde orbits are more stable).

## 2. Citation mining

- Held: Arenstorf 1963 (now held), Hénon 1965b and 1969 (now held), Hénon & Guyot 1970, Hénon & Heiles
  1964, and Szebehely 1967.
- Not held:
  - Jefferys 1965 (AJ 70:393), 1966 (AJ 71:306), 1970a,b and the 1971 Atlas;
  - Goudas 1963 (Icarus 2:1);
  - Moser 1955 and 1962; Arnold 1963; Cushman 1970.
- Added to the wanted list: the Jefferys 1971 Atlas (the full section set) and Jefferys 1965 (AJ 70:393,
  symmetric periodic orbits). Low priority.
