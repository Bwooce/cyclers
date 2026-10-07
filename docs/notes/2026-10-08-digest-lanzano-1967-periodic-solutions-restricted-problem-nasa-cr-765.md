# Digest: Lanzano 1967, "Periodic Solutions for the Restricted Three-Body Problem of Celestial Mechanics" (#960; consumers #944, #948)

P. Lanzano (North American Aviation, Inc., Downey, Calif.), NASA CR-765, Washington D.C., May 1967. Final
technical report on contract NASw-1309. NTRS 19670014390, accession N67-23719. CFSTI price $3.00. No DOI:
a Crossref search (2026-10-08) returns only Lanzano's Icarus papers, not this report.
- Supplied file `6f4ae10e-19670014390.pdf`, 48 PDF pages (covers, front matter iii-vii, report pp. 1-42, back
  cover), md5 c193b5bf638a51bda9eff09a4beea4b1. NTRS scan (Acrobat Capture 3.0). The old text layer has the
  prose in letter-spaced form and drops every equation.
- OCR copy (the file to be filed), made with `ocrmypdf --force-ocr -l eng`:
  `cyclers_pdf/papers/lanzano-1967-periodic-solutions-restricted-three-body-problem-nasa-cr-765-ntrs-19670014390.pdf`, 48 pp.,
  md5 20ded21050861e56841f249fb6cee035, PDF/A-2b, 9.4 MB, made with `-O 3` added (section 6).
  Render check at 100 dpi on PDF pp. 1, 33 and 45: same pixel size as the original, mean grey difference 3.1,
  4.2 and 2.0 of 255, ink fraction equal within 0.001. The page images survive.
- Proposed filename: `cyclers_pdf/papers/lanzano-1967-periodic-solutions-restricted-three-body-problem-nasa-cr-765-ntrs-19670014390.pdf`.
- How I read it: all prose in the old text layer. Every equation and number quoted below was read on 130 dpi
  page renders (report pp. 3-6, 11-12, 14-18, 20-28, 31-33, 35-41), and p.28 again at 300 dpi.
- Not on the wanted list (grep of `2026-10-05-960-wanted-papers.md` for "lanzano": no hit). Not in the corpus
  (`ls papers | grep -i lanzano`: no hit; CORPUS_INDEX mentions Lanzano only in a Broucke 1969 citation list).

## 0. Verdict

**A short theory report with two unrelated parts. It has no numeric table and no orbit, so there is no
`lanzano-1967-tables.yaml`. It gives the project two things: (1) a period-area identity for simple closed
periodic orbits, and (2) a summary (no proof) of Lyapunov-type series for the two L4/L5 families. Neither
touches cyclers, close passes, or Earth-Moon second-species orbits.**

- **Part 1 ("This section deals with the second task of the contract", p.1), "Intrinsic properties of periodic solutions for the circular problem"
  (pp. 1-28).** Planar circular restricted problem, any mass ratio.
  1. A polar series for the zero-velocity curves near a primary (pp. 5-10; eq. (6), recursion (7) on p.8), in
     r0 = 2 mu1/(J - 3 mu2) with mu = mu2/mu1 (p.7). Valid only where r1 < 1, and with mu and r0 small, "which means
     not only that the ratio of the primary masses must be small but also that certain values of the Jacobi
     constant might have to be excluded" (p.7). No number is given.
  2. An intrinsic (frame-free) form of the orbit equation: curvature
     k = -2/V + (1/V^2)(Ux^2 + Uy^2)^(1/2) sin gamma (p.14), with gamma the angle between the tangent and the normal
     to the zero-velocity curves.
  3. **A period identity (p.16):** for a periodic orbit C of period tau that is a simple closed curve (no loops),
     lies wholly inside the allowed region, has no common point with the zero-velocity curve, encloses no primary,
     and is described counterclockwise, "2 tau + 2 pi = -I", I = the double integral over the enclosed area of
     Laplacian(ln V(x, y; J)). Orbits that enclose a primary need a limit with small circles cut out (p.16); the
     result for that case is not written down. **Checked (section 3).**
  4. "Periodic displacement of a periodic orbit" (pp. 18-28): the variational equations (9) along a known periodic
     orbit, with the Jacobi first integral (10) (isoenergetic). In tangent/normal coordinates (p, q) (method cited
     to Smart 1953 and Lanzano 1961) the normal displacement obeys a Hill equation, q'' + Q(t) q = 0, with
     Q(t) = V''/V + 2 P^2(t) + 2 - Uxx - Uyy and P(t) = 1 + (x' y'' - x'' y')/V^2 (p.19), and the tangential part
     follows from d(p/2V)/dt = (q/V) P(t). After **truncation of the periodic Q(t) to its first two Fourier
     terms** (p.20), the problem becomes the
     Mathieu equation (11), q'' + [Q0 - 2 Q1 cos 2t] q = 0. The rest is the classical small-Q1 perturbation
     theory of Mathieu functions: eqs. (12)-(16), and the Floquet relation cos(pi c) = q(pi; Q0, Q1) with
     q(0) = 1, q'(0) = 0 (p.26). **Checked (section 3): eq. (16) is right; the printed q4(pi) has an error.**
- **Scope of the "existence" claim in Part 1.** The Summary (p.v) says that the existence of periodic orbits that
  are an infinitesimal isoenergetic displacement of a given periodic orbit "has been proved by ascertaining
  periodic solutions to a Mathieu equation". The body does less: it truncates Q(t) to two terms ("Neglecting all
  the coefficients but Q0 and Q1", p.20) and then gives the condition (16) for a periodic Mathieu solution with
  rational c ("obviously, for a periodic solution, c must be a rational number", p.21). That is a formal
  first-harmonic approximation, not an existence theorem for the full variational system. No orbit and no mu is
  treated.
- **Part 2, "Periodic solutions about a Lagrangian point" (pp. 29-41).** A summary of a
  generalisation of Lanzano (1965), Icarus 4:223-241 (doi 10.1016/0019-1035(65)90001-1, not held), which holds
  the proofs ("For a rigorous proof of this statement, one should see Lanzano (1965)", p.38).
  - Families: **the two families of periodic orbits about a triangular point** (L4 or L5; the text says "a
    Lagrangian Triangular point" and does not pick one).
  - Mass range: **the smaller mass below 0.03852 of the total** (p.40), so that the linear eigenvalues are pure
    imaginary, lambda1 = i delta1, lambda2 = i delta2 (pp. 32-33). Plus a non-resonance condition: with
    |delta1| < |delta2|, **delta2/delta1 must not be an integer** (p.40; p.37 writes it as lambda2/lambda1).
  - Method: a symplectic linear change D that diagonalises the linear part (p.33), then formal power series in two
    complex variables u, v (eq. 22) with an auxiliary system u' = u alpha(uv), v' = v beta(uv) (eq. 23), solved
    order by order by the recurrence (26). Hamiltonian structure gives beta = -alpha, so uv is constant and
    u = u0 exp(i alpha1 t) (p.38). This is the Lyapunov centre theorem done by normal-form series. The report
    states that the series are valid "within the radius of convergence of their series representation" (p.40);
    it does not prove convergence here.
  - Results stated: synodic period tau = 2 pi/alpha1, alpha1 a series in rho^2 = u0 u0-bar whose leading terms are
    delta1 and delta2 for the two families (p.40). First order: an ellipse centred on the triangular point; its
    eccentricity depends only on the eigenvalues; its major axis lies in quadrants 2 and 4 if mu1 > mu2, and in 1
    and 3 if mu2 > mu1 (pp. 38-39). Second order: "families of oval shaped quartics" (p.39).
  - **Symmetry: not discussed.** The report does not say whether the orbits are symmetric about any axis.
  - The second family comes from "interchanging the relative role played by the two eigenvalues" (p.39). The
    report does not say that the second family then needs delta1/delta2 non-integer, which is automatic because
    |delta1| < |delta2| (our remark).
- **What it gives the project.**
  - `#944`: nothing on second-species or close-encounter families. The Hill-equation reduction (p.19) is the
    standard normal-displacement form; its Mathieu truncation is not a stability test we should adopt. The L4/L5 part is the textbook long- and
    short-period family pair. It is a citation for the Routh limit and the non-resonance condition, not a source
    of orbits.
  - `#948`: nothing. No Earth-Moon number, no collision orbit, no two-body continuation.
  - The period identity (Part 1, item 3) is a cheap, frame-free consistency test for a symmetric-orbit corrector
    on simple orbits that enclose no primary (for example planar Lyapunov orbits). See section 3 for the sign rule.
- **Catalogue (PROPOSAL only):** none. No catalogue row cites Lanzano (grep of `data/catalogue.yaml` for
  "lanzano": no hit), and the report has no orbit to add.

## 1. Set-up and conventions (pp. 3-5, page images)

- Planar circular restricted problem. Synodic frame (0; x, y) rotating about the centroid; P1 at (-mu2, 0),
  P2 at (mu1, 0), mu1 + mu2 = 1, distance P1P2 = 1 (p.4).
- Equations of motion (1): x'' - 2y' = Ux, y'' + 2x' = Uy, with
  U = (x^2 + y^2)/2 + sum mu_j/r_j = (1/2) sum mu_j (r_j^2 + 2/r_j) - (1/2) mu1 mu2 (eq. 2, p.3).
- Jacobi integral (4): V^2 = x'^2 + y'^2 = 2U(x, y) - constant; zero-velocity curves (5):
  sum mu_j (r_j^2 + 2/r_j) = J (p.5). In Part 1, J denotes the Jacobi constant of eq. (5) and V^2 = 2U - J
  (p.15).
- Part 2 uses a different notation (p.29, read in the text layer only: "we have been compelled ... to use a notation different from the one
  employed in the other task"): origin at the triangular point, axes parallel to the synodic axes, Hamiltonian
  form y' = J H_y with J the 4x4 symplectic matrix (pp. 30-31).

## 2. Part 2 formulas read on the page (pp. 32-33)

- Characteristic equation of the linear part: 2 lambda^2 = -1 +/- sqrt(delta), delta = 1 - 27 mu1 mu2 (p.32).
- delta1 = ((1 - sqrt(delta))/2)^(1/2), delta2 = ((1 + sqrt(delta))/2)^(1/2) (p.33); lambda3 = -lambda1,
  lambda4 = -lambda2.
- Eigenvector components c1k = gamma + 2 lambda_k, c2k = -3/4 + lambda_k^2, c3k = 3/4 + gamma lambda_k +
  lambda_k^2, c4k = gamma + (5/4) lambda_k + lambda_k^3 (p.33). gamma is not defined in the pages I read; it is
  presumably the off-diagonal second derivative of the potential at the triangular point (proportional to
  (mu1 - mu2)); I did not verify this.
- Symplectic matrix D = (C1/p1; C2/p2; C1-bar; C2-bar), p1 = i delta1 sqrt(delta) (5/2 - sqrt(delta)),
  p2 = -i delta2 sqrt(delta) (5/2 + sqrt(delta)) (p.33). Not checked.
- **Arithmetic check of 0.03852 (ours):** real delta needs 27 mu(1 - mu) <= 1, so mu <= (1 - sqrt(1 - 4/27))/2 =
  0.0385209. The printed "0.03852" (p.40) is this Routh value. The characteristic equation above is the
  standard lambda^4 + lambda^2 + (27/4) mu(1 - mu) = 0. For the Earth-Moon mu = 0.01215 (ours): delta1 = 0.2982,
  delta2 = 0.9545, ratio 3.20, so the condition holds there.

## 3. Checks (scripts in `checks/`, outputs kept)

### 3.1 Period identity, p.16

- Derivation (ours, from the p.14 curvature formula and Green's theorem): integrate k V dt over one period.
  The left side is the total turn of the tangent, +2 pi for a counterclockwise simple curve, -2 pi for a
  clockwise one. The right side is -2 tau plus a line integral that Green's theorem turns into -I
  (counterclockwise). So 2 tau + 2 pi = -I for counterclockwise orbits, as printed, and **2 tau - 2 pi = +I for
  clockwise orbits**. The report states the counterclockwise case only (p.15, "when described counterclockwise").
- Numerical test, `lanzano_period_integral_check.py` / `_output.txt`: planar Lyapunov orbits about L1 and L2 at
  mu = 0.012150 (ours; these orbits are clockwise and enclose no primary). I is computed on a 1201 x 1201 grid
  with a 5-point Laplacian of ln V inside the orbit polygon.
  Result (all three orbits clockwise, minimum V^2 inside the orbit 0.0007-0.0098, so ln V is regular):

  | orbit | tau (TU) | I (grid) | 2 tau - 2 pi | -(2 tau + 2 pi) |
  |---|---|---|---|---|
  | L1, x0 = 0.826918 | 2.71722 | -0.84530 | -0.84875 | -11.718 |
  | L1, x0 = 0.806918 | 3.08374 | -0.11476 | -0.11570 | -12.451 |
  | L2, x0 = 1.135680 | 3.38752 | +0.49286 | +0.49185 | -13.058 |

  The clockwise form 2 tau - 2 pi = I holds to the grid accuracy (0.1-0.8 per cent). The printed form, applied
  blindly to a clockwise orbit, would be wrong by about 4 pi. So the identity is right, but the orientation
  condition on p.15 is essential: the simple orbits that enclose no primary and that the project meets most
  (planar Lyapunov orbits, linear L4/L5 orbits) are clockwise in the rotating frame. I did not test a counterclockwise orbit numerically; that case rests on the derivation.

### 3.2 Mathieu results, pp. 20-28

- `lanzano_mathieu_check.py` / `_output.txt` (DOP853, rtol 1e-13):
  - **Eq. (16) is right.** I take Q0 from the printed series (16) to order Q1^6 and integrate (11) to t = pi with
    q(0) = 1, q'(0) = 0. Then q(pi) - cos(pi c) is 1e-7 at Q1 = 0.2 and falls by a factor of 256-257 each time Q1
    halves (c = 0.3, 0.5, 1.5), so the error is O(Q1^8), as it should be for a series correct to Q1^6. Eq. (16)
    is also the standard Mathieu expansion of a for non-integer order (McLachlan form).
  - Leading coefficients checked by hand (ours): C1/C0 = -Q1/(4(c+1)), C2/C0 = Q1^2/(32(c+1)(c+2)), C3/C0 =
    -Q1^3/(384(c+1)(c+2)(c+3)) (pp. 24-25), eq. (15) B(c; k) = -A/(4k(c+k)), and the product formula
    Ck/C0 = (-1)^k Q1^k Gamma(c+1)/(2^(2k) k! Gamma(c+k+1)) (p.25) all follow from c^2 - (c+2k)^2 = -4k(c+k).
    The higher-order terms in C1/C0 and C3/C0 were not checked.
  - **q2(pi) = pi sin(pi sqrt(Q0)) / (4 sqrt(Q0) (Q0 - 1)) (p.28) is right** (symbolic derivation from eq. 16,
    `lanzano_q4_sympy.py`; difference 0).
  - **The printed q4(pi) (p.28) is wrong in one coefficient.** Printed:
    q4(pi) = (15 Q0^2 - 32 Q0 + 8) pi sin(pi sqrt(Q0)) / (64 (Q0-1)^3 (Q0-4) Q0 sqrt(Q0))
             - pi^2 cos(pi sqrt(Q0)) / (32 Q0 (Q0-1)^2).
    The derivation from eq. (16) gives **-35 Q0** in place of -32 Q0; the cos term is right. Numerical test: with
    the printed q4 the residual of q(pi) scales as Q1^4 (it does not remove the Q1^4 term); the residual without
    any q4 term, divided by Q1^4, equals the derived q4 (0.04572, 0.02251, 0.003157 at Q0 = 0.3, 2.2, 6.5; the
    printed form gives 0.2552, 0.05439, 0.003295). I read "32" on a 300 dpi render of p.28; the digits are clear.
    So it is a slip in the report (typing or algebra), not an OCR artefact. It does not affect eq. (16) or any
    other result.

## 4. Text items kept as printed

- Summary (p.v) says "has been proved"; the body truncates Q(t) first (section 0).
- p.37 says "the ratio lambda2/lambda1 ... is not an integer"; p.40 says "delta2/delta1". Same condition.
- The two Icarus papers named in the Foreword as "will shortly appear" did appear (Crossref): "Contributions to
  the Elliptic Restricted Three-Body Problem", Icarus 6:114-128 (1967), doi 10.1016/0019-1035(67)90009-7; and
  "Stability of a Class of Periodic Orbits in the Restricted Three-Body Problem", Icarus 7:105-113 (1967),
  doi 10.1016/0019-1035(67)90053-x.

## 5. Citation mining (p.42, reference list read in the text layer only; held/not-held by `ls cyclers_pdf/papers | grep -i` and CORPUS_INDEX grep)

| Reference | Held? | Note |
|---|---|---|
| Lanzano 1965, "Periodic Motion about a Lagrangian Triangular Point", Icarus 4(3):223-241, doi 10.1016/0019-1035(65)90001-1 | no | The proofs for Part 2. Low priority: L4/L5 Lyapunov families are textbook. |
| Lanzano 1967, Icarus 7:105-113, doi 10.1016/0019-1035(67)90053-x ("Stability of a class of periodic orbits") | no | The full paper behind Part 1. Possibly a stability use of the Mathieu reduction. Low priority. |
| Lanzano 1967, Icarus 6:114-128, doi 10.1016/0019-1035(67)90009-7 (elliptic problem) | no | The elliptic-problem version of Part 2. Low priority. |
| Lanzano 1963, "A class of periodic solutions for the restricted three-body problem", Icarus 2:364-375, doi 10.1016/0019-1035(63)90066-6 | no | Not cited by CR-765; found in the same Crossref search. Title suggests a family result; check its scope before ranking. |
| Lanzano 1960 (XI IAC Stockholm, Jacobi integral and terminal guidance), Lanzano 1961 (JAS 8:40-47, Hill's lunar theory and satellites) | no | Not relevant to cyclers. |
| Kopal 1959 "Close Binary Systems"; Smart 1953 "Celestial Mechanics"; Struik 1950; Fluegge 1956; Meixner & Schaefke 1954 | no | Textbooks. |

Proposed wanted-list addition: none required. If `#944` wants the L4/L5 family proofs, add Lanzano 1965 at low
priority.

## 6. Filing note

A plain `--force-ocr` copy was 30.7 MB (24 times the 1.3 MB scan: ocrmypdf rasterises the 1-bit CCITT pages to
400 dpi colour). The filed copy adds `-O 3` (still `--force-ocr`) and is 9.4 MB; its render check is in the
header. Precedent copies in the corpus are 0.9-2.7 MB. The plain copy was deleted. Log: `ocr-lanzano-O3.log`.

*Filed as `cyclers_pdf/papers/lanzano-1967-periodic-solutions-restricted-three-body-problem-nasa-cr-765-ntrs-19670014390.pdf`.*
