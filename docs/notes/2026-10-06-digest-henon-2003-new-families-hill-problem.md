# Digest: Hénon 2003, "New Families of Periodic Orbits in Hill's Problem of Three Bodies" (#960)

M. Hénon (Observatoire de la Côte d'Azur), "New Families of Periodic Orbits in Hill's Problem of Three
Bodies", Celestial Mechanics and Dynamical Astronomy 85(3):223-246 (2003), doi 10.1023/A:1022518422926.
- Filed as `cyclers_pdf/papers/henon-2003-new-families-periodic-orbits-hill-problem-three-bodies-cmda-85-223-doi-10.1023-A-1022518422926.pdf`.
  24 pages, text layer, md5 5b9283c22fea14f924e29f52dbe42217. Supplied by the owner. (The text layer
  drops the symbol Gamma for the Jacobi constant.)
- I read secs. 1-3 and 4.8 from the text layer, and Table XIII (p.243) on the page image. I skimmed the
  other family sections.

## 0. Verdict

The sequel to Hénon's 1969 Hill's-problem study (H5). It extends the symmetric periodic-orbit search to
double- and triple-periodic orbits. It finds seven new families (Ha, Hb, Hc, Hd, He, Hf, Hg) and
completes **family g3** (begun in Hénon 1970). All families end, at both ends, in orbits of growing size
as Gamma goes to minus infinity, and continue into second-species families, which are identified.
- **g3 is the family-identity control named by Anderson, Campagnola & Buffington 2018:** their
  Jupiter-Europa 1:1+/2:2- petal family is Hénon's g3. Table XIII below gives sourced g3 initial
  conditions, which I re-integrated (sec. 3).
- **`#945` R2 / `#953` R7:** no direct bearing. These are single-secondary Hill-problem orbits; R2 and R7
  are multi-moon Jovian triples. Relevance is limited to the single-moon resonant building blocks
  (INFERRED).

## 1. Content (READ)

- **Hill's equations (1), p.223:** xi'' = 2 eta' + 3 xi - xi/r^3, eta'' = -2 xi' - eta/r^3. Jacobi
  constant Gamma = 3 xi^2 + 2/r - xi'^2 - eta'^2 (2).
- **Symmetries (sec. 1.2):** Sigma (about the xi-axis, with time reversal); the eta-axis symmetry; and
  their product (about the origin). The orbits studied are Sigma-symmetric, with two perpendicular
  xi-axis crossings.
- **Method (sec. 1.3 and sec. 2):**
  - Start perpendicular at (xi_0, 0), with eta-dot > 0 from Gamma. The orbit is N-periodic if it
    crosses the xi-axis 2N times and crossing N is perpendicular.
  - A systematic (Gamma, xi_0) grid scan: -3 <= Gamma <= 6, -2.1 <= xi_0 <= 1, steps 0.005 x 0.002,
    1801 x 1551 nodes. Sign changes of xi-dot at crossing N locate the families.
  - This recovers all H5 simple-periodic families (N = 1).
- **Second-species continuation (sec. 3):** the limits are written as arc sequences ({i, -2, e} etc.,
  Hénon 1997 Table 5.1). A misprint in Hénon 1997 Table 5.1 is noted (p.~229).
- **New families (sec. 4):** Ha (maximum Gamma = 4.271428, a reflection at a critical orbit of g
  described twice), Hb, Hc, Hd, Hg, He and Hf, each with a numerical-data table (Tables VI-XII) and its
  second-species limits.
- **Family g3 (sec. 4.8, Fig. 14, Table XIII):**
  - entirely triple-periodic and symmetric about the eta-axis;
  - maximum Gamma = 3.806201; Gamma goes to minus infinity along both branches;
  - one branch tends to {+1, -1}, intersecting family f (described three times) at Gamma = -1.411618
    and 0.015388;
  - an orbit with two vertical ejections at Gamma = 2.238611;
  - the other branch tends to {+1, i, -1, e}.
- **Table XIII (p.243, page image), family g3, N = 3, xi_1 = -xi_0:**

| Gamma | xi_0 | xi_1 |
|---|---|---|
| -4 | -2.490440 | 2.490440 |
| -1.411618 | -1.199879 | 1.199879 |
| -0.7 | -0.779678 | 0.779678 |
| 0.015388 | -0.655072 | 0.655072 |
| 1 | -0.698379 | 0.698379 |
| 2.238611 | -0.752381 | 0.752381 |
| 3.806201 | -0.778960 | 0.778960 |
| 2 | -0.880924 | 0.880924 |
| 0 | -1.200457 | 1.200457 |

- **Final remarks (sec. 5):** most of these orbits are strongly unstable. Hill's problem is a first
  approximation for real systems. Asymmetric orbits are a possible extension.

## 2. Gate relevance

- The petal/g3 identity control for Anderson et al. 2018 is now sourced numerically.
- `#945`, `#953`: background only.

## 3. Positive controls (independently re-integrated, 2026-10-06)

DOP853, rtol = atol = 1e-12, eta-dot_0 from Gamma; the crossing number counts from the start:
- Gamma = 1, xi_0 = -0.698379: crossing 3 at xi = +0.698379 with xi-dot = 1e-6, t = 2.563606.
  Periodic (N = 3).
- Gamma = 3.806201, xi_0 = -0.778960: crossing 3 at xi = +0.778960 with xi-dot = -1e-6, t = 3.321094.
  Periodic.
- Gamma = -1.411618, xi_0 = -1.199879: every crossing is perpendicular at +/-1.199879 (half-period
  2.299322). This is family f traversed three times, as stated for this intersection point.
So the table reads as stated, to the 6 printed decimals.

## 4. Citation mining

Held: Hénon 1997 and 2001 (Generating Families I and II, LNP m52/m65). Not held, already in the wanted
list: Hénon 1969 (H5, A&A 1:223, "Hill's case: periodic orbits and their stability") and Hénon 1970
(A&A 9:24, "Hill's case: non-periodic orbits"). Not held, not added: Strömgren 1935 (background).
