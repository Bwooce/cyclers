# Digest: Birkhoff 1915, "The restricted problem of three bodies" (#960 batch 29)

George D. Birkhoff (Cambridge, Mass.), Rendiconti del Circolo Matematico di Palermo 39:265-334 (1915),
doi 10.1007/BF03015982 (confirmed against Crossref on 2026-10-06: title, journal, volume 39, pages 265-334,
issued 1915). Read at the Palermo meeting of August 1914; printed 9 May 1915.
- Filed as `cyclers_pdf/papers/birkhoff-1915-restricted-problem-of-three-bodies-rend-circ-mat-palermo-39-265-doi-10.1007-BF03015982.pdf`.
  70 pp, md5 706eeced961a52769529477aaaebf011. Image scan of the 2008 Palermo/Springer retro-digitisation
  ("PageGenie", from a TIFF) with an old OCR text layer. The file has no DOI in its metadata.
- **How I read it.** I read the whole paper from the text layer. I then read these journal pages on the page
  image at 110 dpi: 291, 292, 298, 316, 326, 330, 331, 332 and 333. Every formula and statement quoted below
  from those pages was checked there. Sections 2-8, 14-15, 17-19 and 21 I read from the text layer only. I
  quote prose from them, not formulas.
- **Text layer: usable for prose, not for maths.** Page 5 (journal p.269), page 40 (p.304) and page 52
  (p.316) from `pdftotext -layout`: the words are right. Formulas and digits are not. Examples seen on
  p.298: the exceptional mean motion `n = +-k/(k-1)` comes out as `n -~- +___k/(k ~ I)`, and the radical
  in the limit `C > 32^(1/3)` comes out as `1 ~`. Digits (1/I, 0/o, 2/z) swap freely. A re-OCR of this scan is
  unlikely to fix the maths. **Verdict: no
  re-OCR needed for search and prose. Read every equation on the image.**
- Wanted-list row 64 (Birkhoff 1915 + Moser 1953 + Koopman 1927): Birkhoff 1915 is now held. (Wanted-list row numbers in this digest are the batch-28 numbering; the list was renumbered in batch 29.) Row 64 shrinks
  to Moser 1953 and Koopman 1927 (see section 5).

## 0. Verdict

**Classical theory paper. It gives theorem-level support for two claims the project cares about, and it
gives no numbers a catalogue row can use.**

(a) **Existence of periodic orbits that close after many revolutions (section 20, p.330-333).** For small mass
ratio mu and Jacobi constant above a bound, Birkhoff proves an infinite family of symmetric periodic orbits,
labelled by two integers (k, l), that close after about 2k passes of the "rotating ellipse" and 2l turns of
the rotating frame. The condition is a rotation-number window:
`k*rho < 2*l*pi < k*sigma` (same-class crossings; p.331, read on image). At mu = 0 this window reads
`a1^(3/2) < l/k < a2^(3/2)`, which is the plain statement that any rational mean-motion ratio between the
retrograde and direct circular limits is realised. Moser 1953 (title: "...die sich erst nach vielen
Umlaeufen schliessen") is the later, sharper version of this idea. Birkhoff is the origin.

(b) **Symmetric orbits and continuation in mu (sections 11, 12, 20).**
- Section 11 (p.296-298) continues the circular orbits from mu = 0 to mu > 0 by the perpendicular-crossing
  (symmetry) condition `(x, y, t) -> (x, -y, -t)`. This is Poincare/Moulton continuation done with
  Birkhoff's geometric argument. It names the exceptional values where it fails (section 2 below).
- Section 20 (p.333, read on image) states a continuation-in-mu rule for the (k, l) orbits: they "continue to
  exist, with K1 and K2 in the same order, unless, as mu varies from zero to the given value, one of the
  rotation numbers rho, sigma passes through the value 2*pi*l1/k1, 2*pi*l2/k2". That is a bifurcation
  criterion for the whole many-revolution family.

**What it does NOT give.** No numerical orbit, no table, no bound on "mu sufficiently small", no stability
analysis of the (k, l) orbits, no Jacobi-constant values for them, and no statement about orbits that
encircle both primaries or pass near both (all of Parts III-IV stays inside the zero-velocity oval about
one body J; the many-revolution theorem is for J the dominant body). The Poincare-Birkhoff fixed-point theorem is cited from Birkhoff's own 1913-14 papers and
is not proved here.

**Catalogue implication (PROPOSAL only).** None for a row. Use it as a citation: "Birkhoff 1915, section 20"
for "symmetric periodic orbits with arbitrary rational mean-motion ratio exist near the circular orbits
about the primary for small mu". A row that claims theorem support for a many-revolution inner orbit must
state that the window `rho < 2*pi*l/k < sigma` is for mu small and C above a limit, and must check the row's
mu and C against that regime.

## 1. Section-by-section content

**Part I (sections 1-5, p.269-279): equations and regularisation.**
- Section 1: the Jacobi integral and a reduction to third order using the direction angle phi of the velocity.
  Orbits become stream lines of a steady, incompressible flow in (x, y, phi) space.
- Section 2: a conformal change `x + iy = f(u + iv)`, `dt = |f'|^2 d tau` takes the equations (1') and (3')
  into the same form. This is the engine for sections 3-5.
- Section 3: the equations of normal and tangential displacement along an orbit (eq. 10). The normal
  equation is Jacobi's equation of the calculus of variations for the orbits as extremals at fixed C.
- Section 4: Levi-Civita's regularisation of collision with one primary, `x + iy = p^2 - q^2 + 2ipq`,
  rederived non-canonically. Orbits near collision are analytic cusps.
- Section 5: **a regularisation of both collisions at once** (eq. 15, `Z = 2w + ...`, eqs 20-21), with the only
  singularity left at infinity. A note added January 1914 says Thiele 1895 (Astron. Nachr. 138:1-10) had
  already given a double regularisation of the other, less simple, kind. Birkhoff says his algebraic form is
  different. Waldvogel 1967 (held) generalises this Birkhoff regularisation to three dimensions.
- Use: only for collision-safe integration near a primary. Parts III-IV use the simpler Levi-Civita form.

**Part II (sections 6-8, p.279-288): topology of the state space.**
- Zero-velocity regions for all C (cases I-VI, p.281-282): the five critical values at the five
  libration points. The text gives no numbers except `C <= 3` for case VI (motion everywhere), which is
  the known value at L4/L5.
- The manifold of states at fixed C is non-singular except at the five critical C. Case I (motion inside the
  oval about one body): a 3-sphere with antipodal points identified (equivalent to RP^3). Case III (both
  bodies inside one oval): the region between two concentric spheres, antipodal pairing on each. Cases IV-VI
  are listed with their gluings.
- Poincare had a two-to-one representation by "the space of inversion" and missed that it is two-to-one.
  Birkhoff's one-to-one sphere form fixes this (p.266-267).
- Use: none for numerics. This is the setting in which the "ring" and "disc" sections are defined.

**Part III (sections 9-15, p.288-313): the ring transformation.**
- Section 9 (mu = 0, p.288-291): the two-body problem in the rotating frame. Rotating Kepler ellipses with
  `b = V(a)` (eq. 35). Two circular orbits exist for C > 3: a retrograde one (radius a1) and a direct one
  (a2). **As C falls from infinity to 3, a1 rises from 0 to 1/4 and a2 from 0 to 1** (p.291, read on image).
  Arithmetic check: with C = 1/a + 2*sqrt(a) (direct) or 1/a - 2*sqrt(a) (retrograde), a = 1 and a = 1/4 both
  give C = 3.0000 exactly (python).
- Section 10 (p.292, image): the ring `psi = 0` in (a, theta, psi) coordinates, with `a'=a`,
  `theta' = theta - 2*pi*a^(3/2)` (eq. 38). Periodic iff `k/l = a^(-3/2)`; the mean motion k/l is the
  angular rate in fixed space in units of the frame rate. The geometric picture: orbits become positively
  tangent to concentric circles, a ring of two leaves joined at the heavy body.
- Section 11 (p.295-298): analytic continuation of the circular orbits to mu > 0. The orbit that crosses the
  x-axis perpendicularly twice, with its reflection, is a symmetric periodic orbit. **Retrograde circular
  orbits continue for every C > 3. Direct ones continue except at exceptional C** where
  `n = k/(k-1)`, i.e. mean motions 2/1, 3/2, ... (p.298, image). The largest exceptional C belongs to n = 2
  and equals `32^(1/3)`. Check: n = 2 gives a = 2^(-2/3), and C = 1/a + 2*sqrt(a) = 2^(5/3) = 3.17480
  = 32^(1/3) (python). The next two, n = 3/2 and 4/3, give C = 3.0575 and 3.0285.
- Section 12 (p.298-302): for `C >= C1 > 32^(1/3)` the continuation works for mu small and independent of C.
  The distortion is `rho' = rho(sqrt(1-mu) + mu*rho^2*f)`, `theta' = theta + mu*rho^2*g` (eq. 42, image).
  This uses the mean-anomaly parametrisation.
- Sections 13-14 (p.303-311): the ring `T` for mu > 0 as the first-return to positive tangency with the
  family of auxiliary periodic circles; `T` leaves an area integral invariant and **`T = R*U`, a product of
  two involutions, R = reflection in the x-axis** (p.310-311). This product form is the base of the
  symmetric-orbit count in section 20. The paper attributes the ring idea to Poincare and says the area
  invariant is what lets Poincare's last geometric theorem give infinitely many periodic orbits.
- Section 15 (p.311-313): a necessary and sufficient condition for such a ring to exist: along every orbit,
  exterior tangencies to the auxiliary curves must eventually exceed interior ones.

**Part IV (sections 16-21, p.314-334): periodic orbits.**
- Sections 16-18 (p.314-325): **retrograde periodic orbits exist whenever there is a closed oval of zero
  velocity about J.** The Hill-limit case is done first, for `C > 3^(4/3)` (p.316, image). `3^(4/3)` is 4.32675,
  the same Gamma where Henon 1969 (held) starts the Hill families at Gamma = 3^(4/3). The proof is by a
  continuity argument on the point where a projected orbit returns to the x-axis. For general mu it is a
  curve-intersection argument (figure 5). The result: at least one symmetric retrograde orbit with one
  retrograde circuit, tangent to the y-axis direction at two points of the x-axis, convex in Hill's case.
  No bound on mu is needed, and J may be either body: in the Hill limit J is the small body (S is sent to
  infinity, p.314, text layer). So this is the earliest existence theorem for retrograde satellites of the
  small body in Hill's problem, Henon's family f, the Hill DROs. For general mu it covers retrograde
  orbits about either body while its oval is closed.
- Section 19 (p.325-328, p.326 image): a second reduction, to a map of a **discoid** into itself, valid
  while a retrograde orbit exists and mu is small (eq. 60). Brouwer's fixed-point theorem then gives at least one
  direct periodic orbit of simple type, as long as the total number of retrograde circuits grows with time.
  Birkhoff says this is "probable" for all mu with an oval about J, and notes that the ring map T of Poincare
  fails once the direct orbit becomes unstable, whereas the discoid map does not.
- Section 20 (p.328-333): the symmetric orbits. See section 3 below.
- Section 21 (p.333-334): a corollary of the area invariant. For any closed curve Gamma in the state sphere
  (not a stream line), infinitely many orbits pass through a point of Gamma twice. This gives infinitely many
  orbits through any plane point, any point of the zero-velocity oval, and through J (collision) twice.

## 2. Results that bear on (a) many-revolution periodic orbits

- Closure condition at mu = 0 (p.292): `k/l = a^(-3/2)`, k ellipse circuits per l frame turns, line of apsides
  regressing by `2*pi*l/k` per circuit.
- For mu > 0 and C above the section 12 limit, rotation numbers rho (retrograde boundary) and sigma
  (direct boundary) of T replace the two ends of the window. At mu = 0, `rho = 2*pi*a1^(3/2)` and
  `sigma = 2*pi*a2^(3/2)` (p.330, image).
- Orbit counts (p.330-332, image): the four-class crossing scheme (lower or higher passage, at opposition
  or conjunction) gives, for the first class, k and l with `k*rho < 2*l*pi < k*sigma` an orbit with first
  integer 2k and second 2l (or equal submultiples). Three further rules (p.331-332, image) give odd
  integers: `k*rho < (2l+1)*pi < k*sigma` gives (2k, 2l+1); `k*rho < (2l+1)*pi` and `(2l+2)*pi < k*sigma`
  gives (2k-1, 2l+1); `k*rho < 2l*pi` and `(2l+1)*pi < k*sigma` gives (2k-1, 2l).
- The orbit is found as the crossing of the k-th image of a symmetric-crossing segment with another such
  segment on the ring. The proof is the Poincare-Birkhoff style intermediate-value argument, not a
  perturbation series. So it needs only that rho and sigma are the two boundary rotation numbers, with
  rho < sigma (the "twist"), and that rho and sigma lie in (0, 2*pi).
- Regime limits to remember:
  - `C >= C1 > 32^(1/3) = 3.1748`, so every direct mean motion n = k/l must exceed 2 (inner orbits only);
  - mu "sufficiently small", not quantified;
  - J must be the dominant body here (the two-body generator with a1, a2 is about J with mass 1 - mu).
- For a cycler or tour program this is a statement about inner resonant orbits near a primary, not about
  transfer orbits between two bodies. Birkhoff proves nothing about orbits that pass near both primaries.

## 3. Results that bear on (b) symmetric periodic orbits found by continuation in mu

- Symmetry facts used: equations (1) and (3) are invariant under `(x, y, t) -> (x, -y, -t)`. A symmetric
  orbit crosses the x-axis perpendicularly exactly twice (p.328, text layer). The two crossings are labelled
  by four classes. This is the structure behind every modern "symmetric periodic orbit" search in the
  synodic frame (Henon 1965, Broucke 1969).
- Section 11: the continuation of a circular orbit needs the derivative `alpha_2` of the second-crossing
  velocity to be non-zero. Birkhoff shows this derivative vanishes only if the time between crossings is a
  multiple of the half-period in fixed space (p.297). This gives the exceptional mean motions of section 0.
- Section 20 (p.333, image): continuation in mu of a (k, l) orbit with the same integers keeps the order of
  the points K1, K2 on the section "unless rho or sigma passes through `2*pi*l/k`". At mu = 0 the
  intersection is unique, so a point cannot pass through another, and extra intersections can only appear
  or vanish in pairs (even number). This is the first statement of the "orbits exist until they meet a
  resonance or fold" picture that continuation codes look for.
- Not provided: the (k, l) family for a finite, specified mu such as Earth-Moon, or for orbits about the
  light body. Only the single retrograde orbit of sections 16-18 holds for all mu.

## 4. What Arenstorf 1963 takes from it

Source for this section: the held Arenstorf digest
(`docs/notes/2026-10-06-digest-arenstorf-1963-amer-j-math-periodic-solutions-keplerian-continuation.md`,
lines 40-48 and 82) plus my reading of Birkhoff. I did not re-read Arenstorf's paper.
- The Arenstorf digest records that he uses "Birkhoff's **symmetric** periodicity condition (eqs 4-5, 13) in
  place of the general [Poincare] one". In Birkhoff 1915 this condition is the perpendicular-crossing
  construction of section 11 (pp.296-298) and the class scheme of section 20. Arenstorf applies it to
  eccentric Kepler orbits (apocentre start, 0 < eps < 1) instead of Birkhoff's circular generators.
- Arenstorf's mu = 0 generating orbit has `a = (m/k)^(2/3)`. Birkhoff's `a = (l/k)^(2/3)` (with `k/l = a^(-3/2)`)
  is the same commensurability with the roles of the symbols changed. This is consistent.
- Arenstorf's exceptional value (i), `eps = sqrt(1 - a^(-3))` for a > 1 and direct motion, and his
  collision condition (ii), are for eccentric generators. They are not Birkhoff's `n = k/(k-1)` list. I did
  not attempt to show whether the circular limit eps -> 0 of Arenstorf's (i) is Birkhoff's list. That is an
  open reconciliation, not a claim.
- Birkhoff gives existence by a topological argument with no analytic dependence on the parameter. Arenstorf
  gives holomorphic families in (eps, mu). The two are complementary: Birkhoff covers all rational windows
  inside the twist range with existence only. Arenstorf gives smooth families with explicit exceptional
  sets, and works for generators with eps up to 1.

## 5. Citation mining

All held/not-held checks done by `ls cyclers_pdf/papers` and `grep -i` on
`docs/notes/CORPUS_INDEX.md` (no CORPUS_INDEX.md exists under `cyclers_pdf`; the index is in
`docs/notes`). Wanted-list check against `docs/notes/2026-10-05-960-wanted-papers.md`.
- Poincare, Les methodes nouvelles de la mecanique celeste, vol. I (1892), vol. III (1899); and Poincare,
  "Sur un theoreme de geometrie", Rend. Circ. Mat. Palermo 33:375-407 (1912): not held. Poincare is in
  wanted row 47 only as a generic textbook item. The 1912 paper (Poincare's last geometric theorem, the
  fixed-point conjecture Birkhoff then proves) is a **new candidate**, low priority.
- Birkhoff 1913, "Proof of Poincare's geometric theorem", Trans. AMS 14:14-22, and 1914 "Demonstration du
  dernier theoreme de geometrie de Poincare", Bull. Soc. Math. France 42:1-12: not held, not on the list.
  **New candidates**, low priority. They are the source of the fixed-point theorem used here.
- Levi-Civita 1906, Acta Math. 30:305-327 ("Sur la resolution qualitative du probleme restreint des trois
  corps") and 1904, Ann. Mat. 9:1-32: not held (the index has no file with Levi-Civita as author; hits in
  the index are inside digests of other papers). Not on the list. **New candidate**, low priority. Also
  Levi-Civita 1900 Ann. Mat. 5:221-307 (instability criteria): not held, new candidate.
- Thiele 1895, Astron. Nachr. 138:1-10 (double regularisation, note p.265): not held, not on the list. New
  candidate, low priority (historical).
- Hill 1878, Amer. J. Math. 1:5-26, 129-147, 245-260 (lunar theory; the Hill problem): not held, not on the
  list. **New candidate**, medium priority for the Hill-limit baseline, since Henon 1969 (held) builds on it.
- Darwin 1897, "Periodic orbits", Acta Math. 21:99-242: not held. The wanted list has "Darwin 1911" in row 47,
  which is a different paper. Darwin 1897 is a **new candidate**. It is the first large numerical survey of
  retrograde and direct periodic orbits (equal masses and other ratios) that Birkhoff cites as evidence
  (p.267, 307).
- Moulton 1906, Trans. AMS 7:537-577 (periodic solutions of the three-body problem, lunar application) and
  Moulton 1913, Proc. 5th Int. Congr. Math. 2:182-187 (relations among families): not held, not on the list.
  **New candidates**, low to medium. The 1913 paper is the source of the idea that retrograde families
  cannot form cusps.
- Whittaker 1901-02, MNRAS 62:346-352: not held. Wanted row 47 lists "Whittaker" generically (probably
  Analytical Dynamics), so it is a **new candidate** as the specific paper, low priority.
- Brouwer 1910, Math. Ann. 69:176-180 (fixed point of a disc map): not held, not on the list. Textbook-level,
  not worth chasing.
- Moser 1953, Koopman 1927: still on row 64. This paper is held now; they remain wanted at low priority.
- Waldvogel 1967 (spatial Birkhoff regularisation): HELD (`waldvogel-1967-verallgemeinerung-birkhoff-
  regularisierung-...pdf`, digest `2026-10-04-digest-waldvogel-1967-spatial-birkhoff-regularisation.md`).
  It generalises Birkhoff's section 5 and is the practical follow-up for regularised propagation.
- Henon 1969 (HELD, digest `2026-10-06-digest-henon-1969-hill-case-periodic-orbits-stability.md`): shares the
  value `3^(4/3)` with Birkhoff section 16. Cross-link only.
- Szebehely 1967 book (HELD, `szebehely-1967-theory-of-orbits-restricted-problem-three-bodies-book.pdf`):
  covers this material in modern form. I did not search it for Birkhoff-specific pages.

## 6. Wanted-list actions (proposal)

- Row 64: remove Birkhoff 1915. Keep Moser 1953 and Koopman 1927.
- Add one low-priority row: Birkhoff 1913 (Trans. AMS 14) and 1914 (Bull. SMF 42); Poincare 1912 (Rend.
  Palermo 33:375); Levi-Civita 1906 (Acta Math. 30) and 1900 (Ann. Mat. 5); Moulton 1906 and 1913;
  Darwin 1897; Hill 1878; Thiele 1895. Note that only Hill 1878 and Darwin 1897 feed the Hill-limit baseline
  already built from Henon 1969.

## 7. Cleanup record

I rendered journal pages 291, 292, 298, 308, 316, 326, 330-333 as PNG while checking (page 308 was
rendered but not read), and made a full text dump `all.txt`.
I removed all PNG renders and the text dump before finishing. Only this digest remains in the folder.
Python checks (C values at a = 1/4, 1; 32^(1/3) = 2^(5/3); 3^(4/3)) were typed inline and not saved.
