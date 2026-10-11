# Digest: Koopman 1927, "On Rejection to Infinity and Exterior Motion in the Restricted Problem of Three Bodies" (#960 wanted-list row 51)

B. O. Koopman (Harvard University), Transactions of the American Mathematical Society 29(2):287-331 (April 1927), doi
10.1090/S0002-9947-1927-1501390-9. Presented to the Society January 2, 1926; received April 1926. "In substance this paper
is the author's Doctor's thesis, presented at Harvard, 1925." 45 pp.
- **Source: Internet Archive, not AMS.** ams.org returned Cloudflare 403 pages (not bypassed; fetch-log). The article was
  taken from the Internet Archive scan of the whole issue, item `sim_american-mathematical-society-transactions_1927-04_29_2`
  (233 pp., md5 c489c82dae033e5d0ca0ba99fb9980eb, not access-restricted). **This source is not on the brief's list; owner
  to accept or refetch from AMS in a browser.** It is not a shadow library, and a 1927 US printing is in the US public
  domain since 2023.
- Pages 287-331 = issue PDF pages 49-93, cut with qpdf: md5 2e2376db5e8d0397a1cc1f568b8c84a9 (8,007,790 bytes, 45 pp.,
  400 ppi JPEG2000 greyscale). All 45 journal page numbers were checked in the text layer (287-292, 294-331) and page 293
  on its image. Page 1 image shows the title and "BERNARD OSGOOD KOOPMAN"; the issue's next article starts on PDF p.94
  (Starke ends on p.286 at PDF p.48).
- Filed copy: the 8 MB source cut itself (md5 2e2376db5e8d0397a1cc1f568b8c84a9, 45 pp.), whose archive.org text layer is usable (15,477 words). A force-OCR -O 3 copy was 40 MB (five times larger, re-rasterised JPEG2000 pages), so it was not filed.
- Proposed filename: `koopman-1927-rejection-to-infinity-exterior-motion-restricted-problem-three-bodies-trans-ams-29-287-doi-10.1090-S0002-9947-1927-1501390-9.pdf`.
- PDF page n = journal page 286 + n. **How it was read (low priority, per the brief):** the introduction, every theorem
  statement and the chapter openings on page images (110-250 dpi); proofs only skimmed in the OCR text.
- Checks (ours): `work/koopman_escape.py/.out`, `work/koopman_escape_neg.py/.out`.

## 0. Verdict

**Classical qualitative theory of escape ("rejection to infinity") and of motion OUTSIDE the zero-velocity oval that
contains both primaries.** It extends Birkhoff 1915 (held) from interior to exterior motion: for small mu and C > 32^(1/3)
there is a ring-shaped surface of section, infinitely many symmetric periodic exterior orbits classed by two integers, and
infinitely many orbits that go to, come from, or both come from and go to infinity. It has no numerical orbits.

**For the project (PROPOSALS only):**
- **`#944` X2 and `#948` R4:** nothing usable. Koopman's results are for exterior motion, small mu (no bound) and
  C > 32^(1/3) = 3.1748. They do not concern orbits near the second mass or second-species orbits. Arenstorf 1963 cites it
  as part of the classical existence background, and that is its only role here.
- **A sourced escape certificate.** Ch. I Theorems 2-3 give a sufficient condition for escape, valid for any mu: at a point
  with r0 >= 5 (distance from the barycentre, unit separation of the primaries), moving outward, and
  (dr/dt)0^2 >= 2/(r0 - 1), the half orbit goes to infinity without apsides, with an asymptote and a limiting speed. The
  Note on p.299 extends it to any masses of total 1 that stay within unit distance of the origin, moving in any way
  (Jacobi integral not needed). This is usable as a stop or label test in any CR3BP or ephemeris run that has to classify
  escape (for example ejected-flyby or capture studies) in place of an ad hoc "r > R and energy > 0" rule. Our check
  (sec. 2) is consistent with it.
- **Exterior symmetric periodic orbits** (Ch. IV Theorems 1-4) are the exterior analogue of Birkhoff's sec. 20
  many-revolution orbits: existence for every rotation number pair inside a window set by the inner boundary's advance.
  They are outer resonant orbits around both primaries at small mu. Background only; no catalogue row needs them.
- No catalogue change. Wanted-list row 51 (Moser 1953 + Koopman 1927): both now in hand.

## 1. Content

**Frame and notation** (Ch. I sec. 1, pp.288-289). Planar restricted problem; "two bodies S and J, of masses mu and
1 - mu" (note: S has mass mu), circles about the centre of gravity O at unit angular velocity, OS = 1 - mu, OJ = mu.
Fixed axes XY and rotating axes xy (Jacobi form). Jacobi integral x'^2 + y'^2 = 2 Omega - C. eta is the specific energy
of P in the fixed frame (V^2/2 + Pi, with -Pi the force function, as we read eqs. 3-4); G = r^2 dtheta/dt.

**Chapter I. Orbits extending to infinity** (pp.288-304).
- Theorem 1 (p.290): while eta <= eta0 < 0, P stays inside the rotating equipotential oval Pi = eta0, hence within a fixed
  distance of the origin.
- Theorem 2 (p.290-291): for a half orbit Ab from A with (a) r0 >= 5 (> (5 + 17^(1/2))/2), (b) 0 <= psi0 <= pi/2 or
  3pi/2 <= psi0 < 2pi (psi = angle from OP to the tangent), (c) eta >= 0 along Ab: then Ab goes to infinity without apsides,
  theta has a limit, and eta and V^2/2 tend to a finite limit eta_inf >= 0.
- Theorem 3 (p.298): Theorem 2's conclusion holds when (c) is replaced by (c') (dr/dt)0^2 >= 2/(r0 - 1) at A. "A test
  which is directly applicable" to initial conditions. p.298: in terms of the osculating eccentricity, Theorem 1 keeps the
  orbit bounded while e <= 1 - h, and (c) may be replaced by e >= 1 along Ab; footnote: related to the stability of the
  major axis.
- Theorem 4 (p.299): for any h > 0 there are half orbits to infinity with eta_inf > h and with eta_inf < h; Corollary:
  the Jacobi constant of such half orbits may be chosen arbitrarily. Note (p.299): Ch. I sec. 2 holds for any masses of
  total 1 within unit distance of the origin in any continuous motion.
- Theorem 5 (sec. 3, p.300, OCR text): any half orbit that reaches infinity has a limiting direction from the origin, a limiting
  tangent direction and a limiting specific energy (the conclusions of Theorem 2). Further theorems on the asymptote; the
  parabolic-like case eta_inf = 0 has the asymptote at infinity.

**Chapter II. Infinity as a singularity** (pp.304-310). A rational substitution reduces infinity to a singular point of a
first-order analytic system whose right members vanish there, with expansions starting at third degree. The cubic
truncation is solved exactly, but no method extends its properties to the full system (Introduction, p.287). Theorem 1
(p.308): the set of states whose half orbits reach infinity is a 4-dimensional continuum; Theorem 2: eta_inf is continuous
on that set. Theorem 3 (p.310): for any value of the Jacobi constant there are infinitely many orbits going to
infinity with any assigned eta_inf > 0; similarly for orbits coming from infinity.

**Chapter III. The surface of section** (pp.310-323).
- Birkhoff proved a ring-shaped surface of section for motion inside an oval around one body for small mu and
  C > 32^(1/3) (p.311, image). Koopman proves the same for motion OUTSIDE an oval containing both S and J, "with precisely
  the same restrictions for mu and C", by continuity from mu = 0.
- Theorem 1 (p.321, image): when C > 32^(1/3), for each 0 <= mu <= mu0 (mu0 small) there is a ring-shaped surface of
  section for exterior motion. The inner boundary is a closed stream line, a symmetric periodic orbit without double
  points; the outer boundary corresponds to infinity. T advances all points of the inner boundary by more than a constant
  omega > 0 (at mu = 0 the advance reduces to a0^(-3/2) as printed). There are points arbitrarily near the outer boundary
  that never complete a circuit of the ring. The reflection (x, y, x', y') -> (x, -y, -x', y') carries the ring into itself.
- Theorem 2 (p.321): a topological statement on the boundaries of the ring's image (proof by an inversion representation;
  footnote: Poincaré first used it, incorrectly, for interior motion).
- Sec. 4 (p.322): T has a positive area integral invariant (vanishing at the periodic-orbit boundary) and is the product
  of two involutions R and U, R being (x, y, x', y') -> (x, -y, -x', y'). (A paraphrase of Birkhoff 1915 sec. 14.)

**Chapter IV. Periodic and unstable motions** (pp.323-331).
- Sec. 1: symmetric periodic orbits. Right-angle crossings of the x axis fall into classes I-IV (lower or higher apsis,
  right or left of the oval). alpha = number of intersections of a closed stream line with the ring; beta = number of
  circuits of the ring after alpha applications of T. At mu = 0, alpha is the number of rotations of the ellipse in the
  xy plane, beta the number of revolutions of the particle in the ellipse, beta/alpha the mean angular velocity in fixed
  space (p.325).
  - Theorem 1 (p.325): 0 < 2 pi l <= k omega < 2 pi (l + 1) => symmetric periodic orbit, crossings of the same class,
    alpha = 2k, beta = 2l (or equal submultiples).
  - Theorem 2 (p.325): 0 < (2l + 1) pi <= k omega < (2l + 2) pi => classes I and III or II and IV, alpha = 2k + 1, beta = 2l.
  - Theorem 3 (p.326): 0 < (2l + 1) pi <= k omega < (2l + 3) pi => classes I and IV or II and III, alpha = 2k, beta = 2l + 1.
  - Theorem 4 (p.326): same window => classes I and II or III and IV, alpha = 2k + 1, beta = 2l + 1.
- Sec. 2 (pp.326-330): with mu so small that T increases w everywhere inside the ring: the set of points whose orbits go
  to infinity without higher apsides is a 2-dimensional continuum; there are infinitely many orbits that go to infinity,
  that come from infinity, and that both come from and go to infinity with any number of higher apsides (some symmetric);
  orbits that go to but do not come from infinity are not symmetric. An inner limiting set of bounded orbits exists; it has
  positive area if and only if the boundary periodic orbit is "absolutely stable". Theorem 5 (p.328) is an alternative on
  the escape regions; the rest of sec. 2 assumes that the two escape sets are distinct, an unproved hypothesis (p.328:
  "Only Theorem 6 (below) is true without this assumption"). Theorem 6 (p.330): T is intransitive (Birkhoff's sense);
  this does not need the unproved assumption.
- Sec. 3 (p.331): Birkhoff's theorem (infinitely many stream lines through a point of a closed curve twice) extends under
  an extra hypothesis. Corollaries: infinitely many orbits pass twice through any point of the plane not on the oval, or
  through a point of the oval; infinitely many orbits have two higher (or lower) apsides at an assigned distance.

## 2. Check of the escape criterion (ours, `koopman_escape.out`)

Inertial-frame integration (DOP853, rtol 1e-10) of 200 random states per mu (mu = 0.5, 0.0121506, 0.001) with r0 in [5, 25],
outward motion, random tangential speed, and (dr/dt)0^2 = 2/(r0 - 1) exactly (the weakest case allowed), for 2000 time
units. In all 600 runs dr/dt stays positive (minimum over all runs 0.10) and r reaches at least 290. This is consistent
with Theorem 3; it is not a proof.
Negative control (`work/koopman_escape_neg.py/.out`): 50 states per mu at r0 = 5, zero tangential speed, outward, with
(dr/dt)0^2 at half the threshold: 50/50 fall back (dr/dt < 0 within 200 time units) for each mu. So the test can fail, and
the threshold matters. (The criterion is sufficient only; states between half and the full threshold were not mapped.)

## 3. Citation mining

Checked with `ls cyclers_pdf/papers | grep -i` and `CORPUS_INDEX.md`:

| Reference | Held? | Wanted list |
|---|---|---|
| Birkhoff 1915, Rend. Circ. Mat. Palermo 39 (sec. 13, 14, 20, 21) | HELD | n/a |
| Birkhoff 1917, Dynamical systems with two degrees of freedom, Trans. AMS 18:199-300 | NOT held | no row. Low priority |
| Birkhoff 1920, Surface transformations, Acta Math. 43:1-119 | NOT held | no row (also cited by Moser 1953). Low priority |
| Birkhoff 1912, Bull. Soc. Math. France 40:305-323 (recurrent motion) | NOT held | no row |
| Poincaré 1912, Sur un théorème de géométrie, Rend. Circ. Mat. Palermo 33:375-407 | NOT held | no row (last geometric theorem). Low priority |
| Poincaré, Méthodes nouvelles vol. III | NOT held | no row (book) |

- No reference bears on second-species orbits or on cycler construction. No new wanted-list row is proposed.
- Wanted-list row 51: remove (Moser 1953 and Koopman 1927 both received).

*Filed as `cyclers_pdf/papers/koopman-1927-rejection-to-infinity-exterior-motion-restricted-problem-three-bodies-trans-ams-29-287-doi-10.1090-S0002-9947-1927-1501390-9.pdf`.*
