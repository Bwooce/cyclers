# Digest: Howell 1987, "Consecutive Collision Orbits in the Limiting Case mu = 0 of the Elliptic Restricted Problem" (#960, wanted-list row 40)

K. C. Howell (Purdue), Celestial Mechanics 40:393-407 (1987), doi 10.1007/BF01235854, ADS 1987CeMec..40..393H (stamped on
every page image). Received 26 May, accepted 14 July 1987.
- Proposed filename: `cyclers_pdf/papers/howell-1987-consecutive-collision-orbits-limiting-case-mu-0-elliptic-restricted-problem-celest-mech-40-393-doi-10.1007-BF01235854-ads-scan.pdf`.
  - The file to file is the `ocrmypdf --force-ocr -O 3 -l eng` copy `1987CeMec40393H-ocr.pdf`, md5
    0c9fec71df615396aede32e7e248cb32 (515,940 bytes).
  - Source: image-only ADS JBIG2 scan `1987CeMec40393H.pdf`, 15 pp., md5 9eb7d308b53697e312d1152ca10c5273 (245,855 bytes),
    from https://articles.adsabs.harvard.edu/pdf/1987CeMec..40..393H (second try; the first returned HTTP 504).
  - Dark-pixel fractions at 100 dpi, original against OCR copy: p.1 0.0429/0.0451, p.8 0.0375/0.0391, p.15 0.0154/0.0159.
  - PDF page n = journal page 392 + n. Pages 12-13 (Figs. 8-9) are printed rotated 90 degrees.
- **How it was read:** pp.393-396, 401-403 and 406 on 110 dpi page images (model, eqs. 1-15, Fig. 6 caption, Table 1);
  Figs. 2 (p.397) read on the image for family labels; pp.398-400, 404-405, 407 (figures, references) on the text layer.
- **Companion `howell-1987-table1.yaml`**: Table 1 and the Fig. 6 caption. Check: `checks/check_howell.py` / `.out` (ours).

## 0. Verdict

**The mu = 0 elliptic-problem extension of Hénon's 1968 timing condition for consecutive collision orbits, in eccentric
anomalies, with family diagrams at e_p = 0.4 and 0.98 and the hyperbolic (escape-type) case.** It is a sister paper to the
held Gómez & Ollé 1986 (which uses M2's true anomaly). Howell says so (p.394) and gets the same C_ij existence condition
(eq. 8). All four Table 1 values and the Fig. 6 caption reproduce (sec. 2).

**For the project (PROPOSALS only):**
- It closes wanted-list row 40 (line 98) completely. The row's context: cited by Patel-Longuski-Sims 1998 as [25] (the held
  digest `2026-10-05-digest-patel-longuski-sims-1998-mars-free-return-trajectories.md` lists it as not held).
- It adds little beyond the held Gómez-Ollé 1986: the same physics in different variables, plus the e_p = 0.4 diagrams and a
  closed-form hyperbolic family with its rectilinear limit (eqs. 14-15). For any code that reproduces Gómez-Ollé, Table 1
  is a second, independent set of four goldens (same parabola, different e_p, different anomaly).
- Positive control we ran (sec. 2): the same parabola solver gives Hénon 1968's parabolic orbit at e_p = 0 and the
  Gómez-Ollé e_p = 0.5 values, so Howell, Gómez-Ollé and Hénon are mutually consistent.
- Relation to the held Bruno 1994 digest (`2026-10-06-digest-bruno-1994-restricted-3-body-problem-plane-periodic-orbits-book.md`):
  Bruno treats only the circular problem (that digest has no hit for "elliptic restricted", "e_p", Howell or Gómez-Ollé).
  Its Table IV.1 (97 Hénon arc solutions with tau/pi, eta/pi, transcribed in the held file) is therefore the e_p = 0 slice of
  Howell's diagrams, where Howell's gamma equals M2's true anomaly tau. Howell's e_p = 0.4 and 0.98 families are the
  elliptic continuation that Bruno does not cover. We did not compare individual rows.
- Relevance to cyclers: consecutive collision orbits are the mu -> 0 generating orbits of second-species orbits, the
  skeleton of resonant flyby (free-return) geometries, which is why Patel 1998 cites it. No catalogue row uses it. No
  catalogue change.

## 1. Content

**Model** (pp.394-396, read on the images). mu = 0; M1 at the origin; M2 on a direct ellipse a = 1, period 2 pi,
eccentricity e_p. M3 meets M2 at Q and P, symmetric about the line of apsides, at M2 eccentric anomalies -gamma and
+gamma. Four sign flags (p.395): s0 M3 periapsis positive/negative x; s1 M3 direct/retrograde; s2 M3 at periapsis/apoapsis
at t = 0; s3 M2 at periapsis/apoapsis at t = 0; s = s0 s2. M2: x = s3 cos E_p - e_p, y = s3 (1 - e_p^2)^(1/2) sin E_p.
- a and e of M3 from (gamma, eta) (eqs. 2-3); timing condition t_c = gamma - e_p sin gamma = a^(3/2)(eta - s2 e sin eta)
  (eq. 4 as printed) becomes the single equation (5) in (gamma, eta). Eq. 5 as printed has (gamma - s3 e_p sin gamma) on
  the right; eq. 4 as printed has no s3. Ours: the time from M2's starting apse is gamma - s3 e_p sin gamma, so eq. 5 is
  right and eq. 4 omits s3 (a misprint or a shorthand; kept visible).
- Families keep Hénon's names A_i, B_i, C_ij (p.396). Rule: A families s = -s3, B families s = +s3, C families
  s = s3 (-1)^(i+j). Figs. 2-5: (gamma/pi, eta/pi) diagrams for e_p = 0.4 and 0.98, s3 = +1 and -1, each curve marked with
  the sign triplet s0 s1 s2. The families "are clearly distinguishable even at a primary eccentricity of .98".
- Special cases (p.401): (i) the straight line eta = gamma (M3 orbit = M2 orbit, a = 1, e = e_p), shown dashed; (ii) when
  cos gamma = s3 e_p, eq. 5 does not depend on s, giving a pair of solutions symmetric about the y-axis (Fig. 6:
  gamma/pi = 2.36901, eta/pi = .761903, type -+-+); (iii) double points gamma/pi = i, eta/pi = j where C families start or end;
  slopes there by eq. 7, d gamma/d eta = +-[(2 (i/j)^(2/3) - 1 + s3 e_p (-1)^i)/(1 + s3 e_p (-1)^i)]^(1/2); existence of C_ij
  needs j <= [2/(1 - s3 e_p (-1)^i)]^(3/2) (eq. 8), the same as Gómez-Ollé.
- **Hyperbolic case** (pp.402-406): eqs. 9-12. Solutions only for s0 = s1 = -1, s3 = +1 and s0 = +1, s1 = s3 = -1; one
  family in each case (Figs. 8-9). As gamma -> 0, eta -> infinity and the orbit becomes rectilinear:
  eta = -ln[gamma^2/(2 s (e_p - s3))] (eq. 14), |a| = (1 - s3 e_p) gamma^2, e = 1 + [(1 + s3 e_p)/(1 - s3 e_p)] gamma^2/2
  (eq. 15). The other end is the parabola eta = 0 (Table 1), which is also the start of the elliptic family A0.

Table 1 (p.406, image), parabolic M3 orbits:

| e_p | s3 | gamma/pi |
|---|---|---|
| .4 | +1 | .14432 |
| .4 | -1 | .17697 |
| .98 | +1 | .04581 |
| .98 | -1 | .19031 |

## 2. Cross-checks (`checks/check_howell.out`, ours)

- **Independent model.** M3 on a parabola about M1 (GM = 1) with pericentre on the x-axis, passing through M2's position
  at M2's time from its starting apse (gamma - s3 e_p sin gamma); solve for gamma by bracketing and Brent. This does not use
  Howell's eqs. 2-5.
  - e_p = .4, s3 = +1: 0.14432350 (pericentre on -x). Printed .14432.
  - e_p = .4, s3 = -1: 0.17696715 (pericentre on +x). Printed .17697.
  - e_p = .98, s3 = +1: 0.04580529. Printed .04581.
  - e_p = .98, s3 = -1: 0.19031870. Printed .19031; this rounds to .19032, one unit off in the last digit. Kept visible.
  - Each case has exactly one root, and the pericentre side agrees with the hyperbolic-case sign rule of p.403
    (s3 = +1 needs s0 = -1).
- **Eq. 5 as printed**, solved at eta = 1e-3 with s = -s3 (the A-family rule), gives the same four gamma/pi to 5 digits
  (0.14432, 0.17697, 0.04581, 0.19032). With s = +s3 the only root is the trivial gamma -> 0 one.
- **Fig. 6 caption:** cos(2.36901 pi) = 0.400000 = s3 e_p with e_p = .4 (the s-independent case of p.401). Eq. 5 at that
  gamma has the root eta/pi = 0.761903, and both s = +1 and s = -1 give residual below 1e-5 and the same a = 1.81566.
- **Positive control and cross-check with held papers.** The same parabola solver at e_p = 0 gives gamma/pi = 0.163926,
  Hénon 1968 Table 2's parabolic row (0.16393, per the held Gómez-Ollé digest). At e_p = 0.5, converting gamma to M2's true
  anomaly from its starting apse gives tau/pi = 0.231149 (periapsis start) and 0.105578 (apoapsis start). The held
  Gómez-Ollé 1986 digest has printed 0.2318 (text) / 0.2317 (caption), with the digest's own solution 0.23115, and printed
  0.1055. So our solver agrees with that digest's solution, and Howell's p.406 remark holds: in eccentric anomaly the
  apoapsis start has the larger value, in true anomaly the periapsis start does ("the reverse of the plots shown in [9]").

## 3. Citation mining

References (p.407, text layer checked against the image):

| Reference | Held? |
|---|---|
| [1] Hénon 1968, Bull. Astron. 3:377 | HELD (`henon-1968-orbites-interplanetaires-...-bull-astron-3-377-...-french.pdf`, digest 2026-10-04) |
| [2] Perko 1974, SIAM J. Appl. Math. 27:200 | HELD (`perko-1974-...`) |
| [3] Hitzl 1977, A&A 54:47 | HELD (`hitzl-1977a-...-aa-54-47-...`) |
| [4] Perko 1976, Rocky Mtn. J. Math. 6:675 | HELD (`perko-1976b-...-rmj-1976-6-4-675.pdf`) |
| [5] Hitzl 1977, AIAA J. 15:1410 | HELD (`hitzl-1977-...-aiaa-j-15-1410-...`) |
| [6] Hitzl & Hénon 1977, Acta Astronaut. 4:1019 | HELD (`hitzl-henon-1977b-...`) |
| [7] Hitzl & Hénon 1977, Celest. Mech. 15:421 | HELD (`hitzl-henon-1977-...`) |
| [8] Bruno 1981, Celest. Mech. 24:255 | HELD (`bruno-1981-...-celest-mech-24-255-...-ads-scan.pdf`) |
| [9] Gómez & Ollé 1986, Celest. Mech. 39:33 | HELD (`gomez-olle-1986-...`, digest `2026-10-04-digest-gomez-olle-1986-elliptic-restricted-mu-0.md`) |

- All nine references are held. No new lead. The wanted list has none of them.

*Filed as `cyclers_pdf/papers/howell-1987-consecutive-collision-orbits-limiting-case-mu-0-elliptic-restricted-problem-celest-mech-40-393-doi-10.1007-BF01235854-ads-1987CeMec-40-393H.pdf`. Tables: `data/sources/howell-1987-table1.yaml`.*
