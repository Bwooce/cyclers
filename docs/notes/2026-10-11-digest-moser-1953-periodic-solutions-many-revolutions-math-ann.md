# Digest: Moser 1953, "Periodische Lösungen des restringierten Dreikörperproblems, die sich erst nach vielen Umläufen schließen" (#960 wanted-list row 51)

J. Moser (Göttingen), Mathematische Annalen 126:325-335 (1953), doi 10.1007/BF01343166. Received 12 January 1953.
German. English title: "Periodic solutions of the restricted three-body problem that close only after many revolutions".
- Source: GDZ Göttingen, Math. Ann. vol. 126 (PPN235181684_0126), structure LOG_0033, article PDF
  `https://gdz.sub.uni-goettingen.de/download/pdf/PPN235181684_0126/LOG_0033.pdf` (free). 12 pp.: GDZ cover sheet
  (p.1, with the GDZ terms of use) + journal pp.325-335 (PDF p.n = journal p.323 + n). 600 dpi CCITT page images, no text
  layer on the scan pages. Source md5 fc9773b6dac90d10dc3ef793ad0f73e0 (1,266,629 bytes).
- The file to file is the `ocrmypdf --force-ocr -O 3 -l deu+eng` copy, md5 37dad4376db2b078d71f47a214fa974c (1,084,448
  bytes). It keeps the GDZ cover sheet (provenance and licence); drop p.1 if the corpus convention wants the article only.
  Dark-pixel fractions at 100 dpi, source/OCR: p.2 0.1062/0.1062, p.6 0.0484/0.0484, p.12 0.0316/0.0316. Page 2 image
  shows the title and "JÜRGEN MOSER in Göttingen". Language deu+eng, not the brief's eng: the text is German (as the held
  Burrau 1913 copy, OCR'd deu).
- Proposed filename: `moser-1953-periodische-loesungen-restringierten-dreikoerperproblems-nach-vielen-umlaeufen-math-ann-126-325-doi-10.1007-BF01343166.pdf`.
- **How it was read:** all 11 article pages on page images (90 dpi overview); eq. (7) and the theorem statement on a
  300 dpi crop of p.328. The OCR text was used to find things and for the German wording.
- Checks (ours): `work/moser_twist.py/.out` (eq. 7 at mu = 0); `work/moser_sympy.py/.out` (sympy: the D1 expansion on
  p.332 and theta(u, v) of eq. (14) both follow from eqs. (11) and (13) as printed).
- **Full English translation** (articles pages 325-335, not the GDZ cover):
  `moser-1953-periodische-loesungen-restringierten-dreikoerperproblems-nach-vielen-umlaeufen-math-ann-126-325-doi-10.1007-BF01343166-en-translation.tex`
  and `.pdf` (tectonic, 8 pp., no overfull boxes), after the Hénon 1965a translation conventions. Equations were
  copied from 200-300 dpi images. One print slip is marked: eq. (15) prints "a ccs psi" for a cos psi-tilde.

## 0. Verdict

**A local existence theorem: near each circular periodic orbit of the restricted problem with small mu (rotating-frame
rate omega, not resonant of order 1, 3 or 4), there are infinitely many periodic orbits that close only after p
revolutions, p any sufficiently large prime.** It is the sharper, local version of Birkhoff 1915 sec. 20 (held): Birkhoff
used the ring map and Poincaré's last geometric theorem in the large; Moser uses Birkhoff's fixed-point theorem for
area-preserving maps near an elliptic fixed point with non-zero twist ("inhaltstreue analytische Abbildungen im Kleinen"),
which places the new orbits close to the generating circular orbit. The proof's hard step is the twist coefficient c1,
computed at mu = 0: c1 = -3 (omega + 1) pi / omega^3, non-zero (eq. 7). **We reproduce eq. (7) at mu = 0** (sec. 2).

**For the project (PROPOSALS only):**
- **`#944` X2 (second-species enumeration) and `#948` R4 (near-both-masses corrector):** Moser gives nothing for either.
  His orbits are near-circular (small eccentricity), first-kind, collision-free, small-mu (no bound given), and lie in an
  arbitrarily small neighbourhood of the generating circular orbit once p is large enough. They never pass near the second mass. Cite
  Moser 1953 only for "infinitely many many-revolution periodic orbits exist near non-resonant circular orbits for small
  mu", next to Birkhoff 1915 sec. 20 and Arenstorf 1963.
- **Relation to Arenstorf 1963 (held, digested):** Moser's p.328 remark draws the line. At mu = 0 his orbits are
  ellipses with period commensurable with 2 pi, as Arenstorf's are. One could try to get such orbits by Poincaré
  continuation from the ellipses (Arenstorf's later route), but the hypotheses would need proof, and the two sets need not
  coincide: "these lie near the periodic starting solutions [the circular orbits], those near the ellipse orbits". So the
  project's existing citation of Arenstorf 1963 for eccentric commensurable ellipses stays correct; Moser covers the
  near-circular end only.
- **Excluded resonances, stated in mean motions (ours):** with n0 = omega + 1 (inertial mean motion of the circular orbit,
  Jupiter/Moon rate = 1), continuation needs omega != 1/n and the twist theorem needs omega != 3/n, 4/n (n = +-1, +-2, ...).
  For example omega = 1 (n0 = 2, the 2:1 inner circle), 3 (n0 = 4), 1/2 (n0 = 3/2), 3/2 (5/2), -1/2 (n0 = 1/2, outer)
  are all excluded. The theorem says nothing at these low-order commensurabilities (in our opinion, these are where most cycler-relevant
  orbits sit).
- No catalogue change. Wanted-list row 51: Moser 1953 now in hand (Koopman 1927 also fetched by this batch).

## 1. Content (sections I-III)

**I (pp.325-328): problem and result.**
- Planar circular restricted problem, G = 1, masses 1 - mu (P1) and mu (P2), 0 <= mu <= 1, rotating frame with rate 1,
  origin at the centre of mass, primaries on the x1 axis. Hamiltonian (eq. 1, p.326, image):
  H = (1/2)(y1^2 + y2^2) + x2 y1 - x1 y2 - mu/((x1 + mu - 1)^2 + x2^2)^(1/2) - (1 - mu)/((x1 + mu)^2 + x2^2)^(1/2),
  xdot_k = H_{y_k}, ydot_k = -H_{x_k}.
- At mu = 0 the circular orbits (eq. 2): x1 = rho cos omega t, x2 = rho sin omega t, y1 = -(omega + 1) rho sin omega t,
  y2 = (omega + 1) rho cos omega t, with (omega + 1)^2 rho^3 = 1, omega != 0, period 2 pi/omega. By Poincaré's continuation
  method they continue to small mu with the same period "if omega != 1/n for n = +-1, +-2, ...", with initial values
  xi_k(mu), eta_k(mu) as convergent power series and Jacobi value h(mu).
- The map (eqs. 3, 3', 4): start at x1 = xi1(mu) + rho x, x2 = 0, y1 = eta1(mu) + (omega + 1) rho y, H = h(mu); follow to
  the next crossing of x2 = 0 near t = 2 pi/omega; read X, Y the same way. The map X = f(x, y; mu), Y = g(x, y; mu) has
  the origin as a fixed point and is area-preserving (eq. 5: f_x g_y - f_y g_x = 1), because (1) is canonical.
- Fixed points of the p-th iterate that are not fixed points of lower iterates are periodic solutions that cross
  x2 = 0 exactly p times near the starting orbit, i.e. "close only after p revolutions"; and conversely.
- At mu = 0 the eigenvalues of the linear part are conjugate complex, lambda = exp(+- 2 pi i/omega). If lambda^3 != 1 and
  lambda^4 != 1 (which implies lambda != +-1), Birkhoff's normal form (eq. 6) exists: A = a cos chi + b sin chi + P4,
  B = -a sin chi + b cos chi + Q4, chi = c0 + c1 (a^2 + b^2), P4, Q4 of order at least 4.
- Birkhoff's fixed-point theorem (Surface transformations, Acta Math. 43, 1920, ch. III sec. 42ff): if c1 != 0, every
  neighbourhood of the origin contains fixed points of the p-th iterate for all sufficiently large p. With p prime and the
  neighbourhood small, the p images are distinct.
- Eq. (7) (p.328, 300 dpi crop): chi = 2 pi/omega - [3 (omega + 1)/omega^3] pi (a^2 + b^2), lambda = exp(+- 2 pi i/omega).
- **Theorem (p.328, set in roman type; our translation):** "Let omega != 0, 3/n, 4/n (n = +-1, +-2, ...) be fixed.
  Then for every fixed sufficiently small mu and every sufficiently large prime p there are, near the periodic starting
  solution, periodic solutions of (1) that close only after p revolutions, i.e. solutions that meet the plane x2 = 0
  exactly p times in a suitable neighbourhood of the starting point x1 = rho, x2 = 0, y1 = 0, y2 = (omega + 1) rho.
  Consequently there are infinitely many periodic solutions of (1) for every sufficiently small mu."
- At mu = 0 these are ellipses of P3 whose period is commensurable with 2 pi; they close in the rotating frame only after
  many revolutions. Poincaré continuation from these ellipses would need its hypotheses proved, and need not give the same
  orbits (see Verdict).
- The result also holds when H is replaced by any function analytic near each point with x1^2 + x2^2 != 0, 1 that has the
  mu = 0 form above.

**II (pp.329-331): the mu = 0 system and the continuation of the circular orbits** (equations read on 200 dpi images).
- At mu = 0, H and D = x2 y1 - x1 y2 are integrals (note the sign: D is minus the angular momentum). With
  U = x1/r + D y2 and V = x2/r - D y1 (a Laplace-Runge-Lenz-type pair) one has Udot = V, Vdot = -U. The transformation
  T: (x_k, y_k) -> (Phi, H, U, V) (eq. 9; Phi = arctg(x2/x1), H = (1/2)(y1^2 + y2^2) + D - 1/r) is one-to-one near every
  point of the circular orbit (Jacobian -omega (omega + 1) rho != 0, independent of t) and reduces (8) to one first-order
  equation, Phidot = -1 - D^(-3)(1 - u cos(Phi + t) - v sin(Phi + t))^2 (eq. 10), with U = u cos t + v sin t,
  V = -u sin t + v cos t.
- U^2 + V^2 = 1 - 2D^2/r + D^2 (y1^2 + y2^2) = 1 + 2 D^2 (H - D) (eq. 11). The circular orbits are exactly U = V = 0 for
  all t (eq. 12: Phi = omega t, H = h0 = -1/(2 rho) - (omega + 1) rho^2, U = V = 0).
- Continuation with period 2 pi/omega: the determinant of (Phi - 2 pi, U - u, V - v) with respect to (h, u, v) is
  Delta = Phi_h 4 (sin(pi/omega))^2; with D_h = D(3D - 2h)^(-1) = -1/omega and Phi_h = -6 pi/(omega^2 rho^2) at t = 2 pi/omega,
  Delta = -(24 pi/(omega^2 rho^2)) (sin(pi/omega))^2, non-zero for omega != 1/n (n = +-1, +-2, ...) (p.331, image). The
  initial values are then power series in mu. The case "1/omega an integer" was treated by E. Hölder (1931, critical
  circles).

**III (pp.331-335): the normal form at mu = 0.** Expansion of the map to third order in (u, v), an area-preserving
transformation given by a generating function F(x, b) (p.335), and coefficient comparison give chi(a, b) = eq. (7).
Footnote 7 (added in proof): a shorter route that also gives more than chi is in Moser's forthcoming paper "Über
periodische Lösungen kanonischer Differentialgleichungssysteme", Nachr. Akad. Wiss. Göttingen 1953.

## 2. Check of eq. (7) at mu = 0 (ours, `moser_twist.out`)

Kepler motion (G = M = 1) in axes rotating at rate 1. Moser's section coordinates x = x1/rho - 1, y = y1/((omega + 1) rho)
make the map area-preserving, so a^2 + b^2 is the area enclosed by the invariant curve divided by pi. For an ellipse
(a, e) the invariant curve is {(r(M), rdot(M))}; the map advances the mean anomaly by 2 pi n/(n - 1) per return; the
Jacobi value H = -1/(2a) - sqrt(a(1 - e^2)) is held at its circular value. Slope d chi/d(a^2 + b^2), at e = 1e-3, 3e-3,
1e-2, against eq. (7):

| omega | ours (e = 1e-3) | ours (e = 1e-2) | eq. (7): -3 (omega + 1) pi / omega^3 |
|---|---|---|---|
| 2.0 | -3.534288 | -3.533872 | -3.534292 |
| 0.7 | -46.711433 | -46.681960 | -46.711728 |
| -0.35 | 142.880819 | 142.671833 | 142.882931 |
| -0.5 | 37.698867 | 37.674622 | 37.699112 |
| -0.2 | 942.430680 | 937.790106 | 942.477796 |

Of these omega, only 0.7 and -0.35 satisfy the theorem's hypothesis (2 = 4/2, -0.5 = 1/(-2), -0.2 = 1/(-5) are
excluded); the check tests formula (7) itself, which holds for every omega != 0. The small-e limit agrees with eq. (7) to 1e-6 relative at omega = 2, and the leading term 2 pi n/(n - 1) = 2 pi + 2 pi/omega
gives Moser's c0 = 2 pi/omega (mod 2 pi). The same result follows by hand: at fixed H, dn = (3/2) n e^2/(n - 1) and
d chi/dn = -2 pi/(n - 1)^2, so d chi = -3 pi n e^2/(n - 1)^3 with a^2 + b^2 = e^2 to leading order. This checks the formula;
it does not check Moser's normal-form algebra in sec. III.

## 3. Citation mining

Checked with `ls cyclers_pdf/papers | grep -i` and `CORPUS_INDEX.md`:

| Reference | Held? | Wanted list |
|---|---|---|
| Birkhoff 1915, Rend. Circ. Mat. Palermo 39 (printed "1-70") | HELD (`birkhoff-1915-...-39-265-...`; Moser's "1-70" is probably the offprint numbering of pp.265-334; INFERRED) | n/a |
| Birkhoff 1935, 1936, "Sur le problème restreint des trois corps" I, II, Ann. Sc. Norm. Sup. Pisa (2) 4:267-306, 5:1-42 | NOT held | no row. Birkhoff's later restricted-problem papers; low priority |
| Siegel, Himmelsmechanik lectures, Göttingen 1951/52, pp.175-176 | NOT held | no row (lecture notes) |
| Birkhoff, Dynamical Systems (1927) | NOT held | no row (book) |
| Birkhoff 1920, Surface transformations, Acta Math. 43:1-119 | NOT held | no row. Source of the fixed-point theorem; low priority |
| E. Hölder 1931, "Die Verzweigungsgleichungen für die kritischen Kreise des restringierten Dreikörperproblems", Ber. Verh. Sächs. Akad. Wiss., math.-phys. Kl. 83:179-184 | NOT held | no row. The integer-omega (resonant circle) case; low priority |
| Moser 1953, Nachr. Akad. Wiss. Göttingen (footnote 7) | NOT held | no row. Shorter derivation; low priority |

- Wanted-list row 51: Moser 1953 received (this batch). Koopman 1927 also received (separate digest).
- Birkhoff 1915 digest, sec. 5: "Moser 1953 ... still on row 64 [now 51]" can be closed.

*Filed as `cyclers_pdf/papers/moser-1953-periodische-loesungen-restringierten-dreikoerperproblems-nach-vielen-umlaeufen-math-ann-126-325-doi-10.1007-BF01343166.pdf` (with the other files named above).*
