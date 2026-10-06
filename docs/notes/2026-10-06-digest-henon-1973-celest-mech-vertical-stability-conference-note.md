# Digest: Hénon 1973, "Vertical Stability of Periodic Orbits in the Restricted Problem" (Celest. Mech. conference note) (#960 batch 29)

M. Hénon (Observatoire de Nice), Celestial Mechanics 8:269-272 (1973), doi 10.1007/BF01231427. Copyright 1973 D. Reidel.
Title as printed on the page image: "VERTICAL STABILITY OF PERIODIC ORBITS IN THE RESTRICTED PROBLEM*" (no "I. Equal masses" subtitle).
Footnote: "Presented at the Conference on Celestial Mechanics, Oberwolfach, Germany, August 27-September 2, 1972."
- Filed as `cyclers_pdf/papers/henon-1973b-vertical-stability-periodic-orbits-restricted-problem-oberwolfach-note-celest-mech-8-269-doi-10.1007-BF01231427.pdf`. Original upload, 4 pp, md5 0e53f02caf5a4463b8a31e6416196b39.
  It is a clean Springer typeset PDF. Filed in batch 29 (the earlier held Hénon 1973 file is the A&A paper, see below).

**How I read it.** I rendered all 4 pages at 200 dpi and read the images. The text is short and clean, so all formulas
and figure readings below come from the images. The two figures are line plots with no numeric tables. Figure readings
are approximate (about 0.02 in x) and are checked against the A&A tables.

## 0. Verdict

**This is the conference summary that announces the A&A 28:415 paper. It is not an erratum and not a different case.**
- Its last sentence is "A fuller account of these results will be given elsewhere." The fuller account is the held
  Hénon 1973 A&A 28:415-426 ("... I. Equal Masses"; digest `2026-10-06-digest-henon-1973-vertical-stability-equal-masses.md`).
  The A&A paper is dated "received 20 July 1973", after this Oberwolfach talk of August 1972. I infer the link from the
  matching title, the matching notation (families, m1, m2v, m1v) and the same example orbit. The note does not name
  the A&A paper.
- **Scope difference.** This note treats the six main families f, g, h, i, l, m "in Strömgren's notation" and a range of
  mu. The A&A paper (Part I) is restricted to mu = 0.5 (equal masses). So the note has two things the A&A paper
  lacks: (1) Fig. 2, the movement of the vertical critical orbits with mu, and (2) the Bray-Goudas link (section 3).
- **What it gives the project:** nothing new numerically. It gives the clean statement of the vertical-stability
  criterion and the physical argument for it (section 1), plus a pointer to the mu-dependence figure.
- **Catalogue implication (PROPOSAL only):** none. No catalogue row cites it. If a row ever needs a vertical-stability
  citation for a planar symmetric orbit, cite the A&A paper for numbers and this note for the criterion.

## 1. Content (all read on page images, pp. 269-272)

- Motivation (p. 269): planar stability is not enough. "The reduction to a plane problem is physically meaningful for
  the periodic orbit itself, but not for the study of its stability." In the linear approximation the planar and
  vertical perturbations decouple, so vertical stability can be studied on its own.
- Vertical variational equations, with rotating-frame coordinates x, y, z and velocity u, v, w (p. 269):
  d(Dz)/dt = Dw, d(Dw)/dt = Omega_zz Dz, with
  Omega_zz = -(1-mu)[(x+mu)^2 + y^2]^(-3/2) - mu[(x-1+mu)^2 + y^2]^(-3/2).
  (Dz, Dw are Hénon's delta-z and delta-w. Omega_zz is the second z-derivative of the effective potential on z = 0.)
- Monodromy (p. 270): (Dz1, Dw1) = [[a_v, b_v],[c_v, d_v]] (Dz0, Dw0), with a_v d_v - b_v c_v = 1.
  Subscript v separates these from the planar coefficients of Hénon 1965 (Ann. Astrophys. 28:992 = Hénon 1965b, held).
- From a perpendicular crossing of the x-axis, a_v = d_v. The eigenvalue equation is
  lambda^2 - 2 a_v lambda + 1 = 0.
  - |a_v| > 1: real roots, one of modulus above 1, vertically unstable.
  - |a_v| < 1: complex pair of modulus 1, "vertically stable in general".
- Vertical critical orbit: an orbit on a family where a_v crosses +-1. Such an orbit "in general ... represents the
  intersection of a family of plane periodic orbits with a family of three-dimensional periodic orbits" (p. 272).
- The planar index is a (Hénon 1965); planar stable when |a| < 1. Planar critical orbit m1.

### Figure 1 (p. 270), family m at mu = 0.5, a_v (solid) and a (dashed) against x
x is the abscissa where the orbit crosses the right-hand x-axis. Read on the image:
- a_v falls from large positive values at small x, passes a_v = 1 at m2v (x about 0.87 on the plot), reaches a minimum
  tangent to a_v = -1 at m1v (x about 1.25), then rises toward a_v about 0.75 at x = 4.
- The dashed a curve passes a = 1 at m1 (x about 0.62), has a minimum near -0.8, and merges with a_v at large x.
- Text: orbits of large size are vertically stable. "There is a definite interval on the family, between m2v and m1,
  where the orbits are stable in the plane but not in the perpendicular direction." (m2v lies at larger x than m1.)
- "At m1v, the curve is tangent to a_v = -1; this is a peculiarity of the case mu = 0.5."

**Cross-check against the A&A tables (checked).** The held A&A companion table
(`henon-1973-vertical...-tables.txt`, lines 176-178, which the A&A digest says were checked on the image) lists, for
family m: m1v x0 = 1.247872 (type A_v-D_v), m2v x0 = 0.859393 (type C_v), m1 x0 = 0.617017 (type C), taking absolute
values (the sign lost in the OCR witness was settled in that digest). My image readings of Fig. 1 (m1v 1.25, m2v
0.87, m1 0.62) agree within the plotting resolution. I did not re-read the A&A tables on the A&A page images.

### Figure 2 (p. 271), vertical critical orbits against mu (0 to 1)
- Curves m1v, m'1v and m2v (and their mirror images at negative x), with the primaries M1 and M2 as dashed lines.
- Text: "Periodic orbits are vertically unstable inside m2v and also between m1v and m'1v."
- The m1v and m'1v curves cross at mu = 0.5 (seen on the image), which is why mu = 0.5 is special (the tangency to
  a_v = -1 above).
- For families f, g, h, i, l: "instability appears first in the plane and only later in the perpendicular
  direction. Vertical instability is generally much milder than horizontal instability: a_v never reaches very high
  values." Family m is the exception (vertical instability first).

## 2. Link to three-dimensional families

"Vertical critical orbits m1v and m'1v for mu = 0.4 have been found to be the end orbits of the three-dimensional
families originating from L3 and L2, studied by Bray and Goudas (1967)." (p. 272, read on the image.)
This is the same mechanism as the project's bifurcation notes: a vertical critical orbit is a branch point to a 3D
family.

## 3. Citation mining

The note has only two references (p. 272):
- Bray, T. A. & Goudas, C. L. (1967), Adv. Astron. Astrophys. 5:71. Not held (no hit in `ls | grep -i bray\|goudas`
  and none in CORPUS_INDEX). Not on `2026-10-05-960-wanted-papers.md` (no hit for Bray or Goudas). **New candidate.**
  Reason: 3D families from L2 and L3 and their end orbits at vertical critical orbits.
- Hénon 1965 (Ann. Astrophys. 28:992) = Hénon 1965b: HELD (`henon-1965b-exploration-numerique-probleme-restreint-II-...`).
  Hénon 1965a is also held. Neither is a wanted row.
- Strömgren's family letters are used without a citation. Strömgren 1933 is HELD (`stromgren-1933-connaissance-actuelle-...`).

Unresolved: whether the A&A paper's introduction cites this note. I did not check that.
