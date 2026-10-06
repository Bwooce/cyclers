# Digest: Kevorkian & Lancaster 1968, "An Asymptotic Solution for a Class of Periodic Orbits of the Restricted Three-Body Problem" (#960 batch 27)

J. Kevorkian and J. E. Lancaster (University of Washington), Astronomical Journal 73(9):791-806
(November 1968), doi 10.1086/110701, ADS 1968AJ.....73..791K.
- Filed as `cyclers_pdf/papers/kevorkian-lancaster-1968-asymptotic-solution-class-periodic-orbits-restricted-three-body-aj-73-791-doi-10.1086-110701.pdf`.
  16 pages, ADS scan with an OCR text layer, md5 d8bd87028217821748f91bdc8fe99f2f (fetched by
  fetch-sonnet).
- I read the abstract and secs. I, V (part), VI and VII from the text layer. The equations are garbled
  in the layer: read them on the page images before using any.
- Wanted list: removed in batch 27.

## 0. Verdict

**The earliest matched-asymptotic construction of Earth-Moon periodic orbits that pass close to both
primaries: a one-parameter family with 1:2 commensurability in the Earth-Moon synodic frame, small mu.**
- It builds on Lagerstrom & Kevorkian 1963 (Earth-to-Moon matched expansions).
  - It uses a more general matching that keeps the Jacobi integral uniform to O(mu), and adds a solution
    for perturbed highly eccentric ellipses.
  - With the x-axis symmetry, it gives a family parameterised by the angular momentum near the Earth
    (lambda).
- **Checked numerically** (sec. VI) against Boeing Scientific Research Laboratories exact integrations
  (Goudas and Bray's program):
  - initial conditions at the crossing C2 for lambda = 0.50-1.75 (Fig. 3; three computed periodic
    orbits, 10 integrations, 47 min);
  - the asymptotic forms x' = -1.719 + 0.0486 lambda - 0.00421 lambda^2 and dy'/dt = 1.761 -
    0.1122 lambda + 0.00243 lambda^2 (eqs. 120a,b; mu = 0.0122);
  - the integral constants I_11 = -0.382, I_22 = 2.73, I_33 = -0.407 and I_44 = 1.487 at rho_0 = 0.7584
    (eq. 118).
- The second crossing passes within 1e-3 of the Moon (perilune). Elements near the Moon are compared in
  Fig. 4.
- **"The techniques ... can be trivially generalized to arbitrary commensurabilities."** This is an
  asymptotic existence-and-construction result, not a catalogue.
- **Use:** `#948` R4 / `#944` X2 attribution: the first analytic Earth-Moon "both primaries" periodic
  family. It is the asymptotic counterpart of Arenstorf's numerical remark (Arenstorf 1963, held). It
  predates Breakwell-Perko 1974 and Perko 1974 (held).

## 1. Citation mining

- Held: Arenstorf 1963a (as Amer. J. Math.; TN D-1859 is its NASA form) and 1963b, Birkhoff 1915 (wanted),
  Breakwell & Perko, Szebehely.
- Not held: Lagerstrom & Kevorkian 1963a,b and 1964; Lancaster 1968 PhD (Univ. Washington); Kevorkian &
  Brachet 1968; Barrar 1965 (AJ 70:3); Shi & Eckstein 1967 (AJ 72:685); Perko 1964 (wanted); Cole 1968.
  Added: Lagerstrom & Kevorkian 1963a and Lancaster 1968 PhD (one low-priority row).
