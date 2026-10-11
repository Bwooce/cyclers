# Digest: Llibre and Olle 2001, "The motion of Saturn coorbital satellites in the restricted three-body problem" (A&A 378:1087-1099)

**Digested:** 2026-10-11, owner upload, read in full from the native text layer (two extractions). Pages 5 to 11 were
rendered at 150 dpi and every number the paper prints was read on them. Two independent checks were run (below).
**Citation:** Jaume Llibre (Universitat Autonoma de Barcelona) and Merce Olle (Universitat Politecnica de Catalunya),
Astronomy & Astrophysics 378:1087-1099 (2001), DOI 10.1051/0004-6361:20011274. Received 7 March 2001, accepted
3 August 2001.
**Filed as** `cyclers_pdf/papers/llibre-olle-2001-motion-saturn-coorbital-satellites-restricted-three-body-problem-aa-378-1087-doi-10.1051-0004-6361-20011274.pdf` (13 pages, native text layer, md5
5d0787cd421c6ff80736f01b6e34cbd2). Check scripts and outputs are filed beside it as
`cyclers_pdf/papers/llibre-olle-2001-motion-saturn-coorbital-satellites-restricted-three-body-problem-aa-378-1087-doi-10.1051-0004-6361-20011274-check_llibre.py`, `-check_llibre.out`, `-check_lpoints.py` and `-check_lpoints.out`.
Value transcription: `data/sources/llibre-olle-2001-values.yaml`.

## Verdict

The paper has no numbered tables. Its numbers are in the text and the figure captions, and all are transcribed. It is a
horseshoe-orbit paper only. There are no tadpole orbits in it and no quasi-satellite orbits. It computes symmetric
horseshoe periodic orbits at the Saturn-Janus mass parameter 3.5e-9 and at other mass parameters. It prints one set of
full initial conditions (at 6.5e-9, not at the Janus value). For the Janus-Epimetheus orbit it prints only x at the two
axis crossings and a period. It gives no collision with a catalogue row, because the catalogue has no horseshoe,
tadpole or Janus-Epimetheus row, and it shares no mass parameter with the held Schanzle and Taylor sources except in
one manifold figure.

## What the paper does

- **Model.** Planar circular restricted problem with Saturn and Janus as primaries (mu = 3.5e-9) and Epimetheus as the
  massless body. A horseshoe periodic orbit is defined as a symmetric periodic orbit that surrounds only L3, L4 and L5
  and crosses the x axis exactly twice, both times orthogonally.
- **Mechanism (sec. 2).** For mu = 0 every sidereal ellipse with semi-major axis (p/q)^(2/3) is a periodic orbit in the
  rotating frame with synodic period 2 pi p. For small mu and C close to 3, pairs of such orbits (circles of radius
  r2 below 1 and r3 above 1 are the simplest) join through a close pass to the small primary, which gives one
  horseshoe orbit. This is the paper's answer to where the horseshoe orbits come from and where they sit.
- **Computation (sec. 3).** Fix C, start at x slightly above the zero-velocity-curve crossing x_o, step x by 1e-7 until
  the sign of x' at the next axis crossing changes, and refine. Levi-Civita regularisation is used near either primary.
  Stability from the half-period monodromy matrix.
  - **Janus and Epimetheus.** There are thousands of horseshoe orbits for x in [1.0002, 1.0004], few stable. One stable
    orbit is shown: x(0) = 1.0002979809, x(T/2) = 0.9996731232, printed "T = 1465.3 days", minimum angular separation
    4.42 degrees against the real 5.6, eccentricity 0.9e-4. Real period 2928.6 days.
  - **Families.** At fixed C the orbits form a continuum in x (Fig. 7), from the zero velocity curve inward, ending where
    the return point reaches the small primary and the orbit oscillates between y below and above 0. Near the L2
    Lyapunov orbit the x dependence is extremely thin (orbits turn down or up for x differing by 1e-13).
  - **Range of C.** Horseshoe orbits exist for C below 3 as well (no zero velocity curve, loops cross) but those are
    unstable. They vanish as C approaches C1 (none at C = 3.000015 for mu = 6.5e-9). Many orbits sit near alpha = -2,
    so bifurcations are dense and were not followed.
- **Other mass parameters (sec. 4).** Same behaviour for mu up to 0.01. For 3 <= C < C3 three regimes in mu, built from
  the invariant manifolds of the L3 Lyapunov orbit: horseshoe-shaped manifolds and homoclinic orbits for mu up to
  0.004 (shown at the Sun-Jupiter value 0.000953875); mixed, with collisions, for mu in [0.0045, 0.099] (shown at 0.01);
  no horseshoe periodic orbits at mu = 0.1 for C3 <= C <= C1, the manifold flow no longer horseshoe-shaped.
  The paper claims new homoclinic orbits not found in Henrard 1983.

## Values transcribed (all in the yaml)

Collinear points and Jacobi constants for mu = 3.5e-9 and 6.5e-9 and 0.01; the Janus-Epimetheus orbit; two
mu = 6.5e-9 horseshoe orbits with their generating circles and ellipses; the L2 Lyapunov orbit and the sequence of orbits
near its stable manifold (x = 1.0035181044137 and 1.0035181044138); the C = 3.00001 count of 29 orbits; two mu = 0.01
horseshoe initial x values for C = 2.9 and two for C = 3.01.

## Our checks (marked ours; nothing from them is a golden value)

- **Collinear points.** All printed x(L1), x(L2), x(L3), C1, C2, C3 for mu = 3.5e-9 reproduce. C1 and C3 for
  mu = 6.5e-9 and 0.01 reproduce (C1 is printed truncated, not rounded). The Sun-Jupiter C used for Fig. 16 is just below
  C3, as the section needs.
- **The printed mu = 6.5e-9 orbit (Fig. 9c,d).** The Jacobi constant of its initial state is 3.000000065000049, so it is
  the orbit of the C = 3.000000065 family. Integration (DOP853, tolerances 1e-13 and 1e-15) to the first return to y = 0
  gives x = 0.999705835 and y' = 4.55692027e-4 with x' of -1.3e-7. The printed y'(T/2) matches to its digits. **The
  printed x(T/2) = 0.999795834471204 does not match.** It differs in the fourth decimal, and agrees with my value if
  the third digit after the decimal point is read as 0 instead of 9. The paper's own direct generating circle
  r = 0.999705659 supports that reading. I read the digits on the page render, and they are as printed. This is a
  probable typesetting slip, not an erratum claim.
- **The printed "T = 1465.3 days".** INFERRED, not computed: it is the half period. The real full period is 2928.6 days,
  the text calls the half period "close to 1464", and Fig. 5a plots T/2 in days. The paper's symbol T means the full
  period elsewhere.
- **Not checked.** The Janus-Epimetheus orbit cannot be reproduced from what is printed, because y'(0) and C are not
  given. A shooting solve for y' that returns x(T/2) = 0.9996731232 was not done.

## Relevance to this project

- **Held horseshoe sources.** Schanzle 1967 and Taylor 1981 (digests of 2026-10-11) are at mu = 0.000953875 with printed
  ICs to 9 and 5 digits. This paper is at mu = 3.5e-9 and 6.5e-9 (and 0.01 for pictures) and prints ICs at 6.5e-9 only. The
  mass parameters differ by five orders, so no row can match another literally. The three together give the same
  family structure at three mass ratios: a continuum in x at fixed C near L3, symmetric with two orthogonal crossings,
  with collision orbits and non-horseshoe symmetric orbits interleaved.
- **Taylor's finding against this paper.** Taylor 1981 reports many distinct symmetric families h(n,m) at the Sun-Jupiter
  value. This paper reports a continuum at fixed C and does not classify families by loop counts. It cites Taylor
  only for segments of several families at mu_SJ.
- **Saturn moon rows and #523.** The catalogue and the corpus index hold no Janus, Epimetheus, horseshoe or tadpole row.
  `#523` (co-orbital-exchange cyclers) names Janus-Epimetheus and is not built. This paper is the only held numeric source
  at the Janus mass ratio, and it supports the encounter criterion question in `#523`: its horseshoe orbits pass the small
  primary at a separation set by x (minimum angular separation 4.42 degrees for the fitted orbit, which is about 0.077
  of the orbital radius, ours), well outside a Hill radius of order mu^(1/3), about 1.5e-3 (ours). Under the `#523`
  criterion of passing inside the Hill sphere, the fitted Janus-Epimetheus orbit does not encounter Janus. Only orbits
  near the end of the x continuum, where the return point reaches the small primary, would.
- **Positive control.** The printed mu = 6.5e-9 initial state is usable as a control for a corrector or monodromy code at
  a very small mass ratio: it closes to the printed y' and to my x = 0.999705835, with the one printed digit slip above.
  It is a sourced initial state, not a golden outcome, because the printed x(T/2) is wrong.
- **Mass-parameter window.** The paper's claim that horseshoe periodic orbits exist up to mu about 0.1 and not at
  mu = 0.1 (shown for one C range, mu = 0.1, C3 <= C <= C1) is a bound only as stated. It is not evidence for the
  Pluto-Charon ratio, which the paper does not treat.
- **No cycler content.** A horseshoe orbit's two near passes per period are a slow exchange, not a cycler pass. No
  catalogue change is proposed.

## Figures (described)

- **Fig. 1 (p. 1088).** Zero velocity curves for mu above 0 in the rotating plane, with L1 to L5 and the two primaries;
  regions 3 < C < C3 and C3 < C < C1 marked.
- **Fig. 2 (p. 1089).** mu = 0, C above 3: the allowed region inside radius r_i below 1 and outside r_o above 1, with the
  circular orbits of radius r1 (retrograde), r2 (direct) and r3 (retrograde).
- **Fig. 3 (p. 1089).** The (a, b) semi-axis plane for fixed C: curves of sidereal ellipses, with a = (p/q)^(2/3), that
  become periodic in the rotating frame, for C above, at and below 3.
- **Fig. 4 (p. 1089).** Characteristic curves C against x of the circular families l, i, h and m.
- **Fig. 5 (p. 1090).** For C = 3.00000006, 3.00000007 and 3.00000008, over x in [1.0002, 1.0004]: minimum angular
  separation, half period in days (1200 to 1750) and the stability parameter, one point per horseshoe orbit.
- **Fig. 6 (p. 1090).** The fitted Janus-Epimetheus orbit: a U-shape in y below 0 from x about 1 round to x about -1, a
  detail of the outer arc and a detail of the return near the small primary.
- **Fig. 7 (p. 1091).** Initial x (1.0005 to 1.0025) against minimum angular separation (0 to 9 degrees), mu = 6.5e-9,
  C = 3.000000065: a curve that peaks near 8.8 at small x and falls toward 0.
- **Fig. 8 (p. 1092).** The mu = 0 model orbit with p = 21, q = 20 over half a period: a ten-loop rosette hugging a
  circle, with and without a zero velocity curve.
- **Fig. 9 (p. 1092).** Pieces of the mu = 6.5e-9 horseshoe orbit against its generating ellipses and circles.
- **Figs. 10 and 11 (p. 1093).** An orbit at C = 2.9999999 with intersecting loops, and Levi-Civita views of the return
  point near the small primary.
- **Figs. 12 and 13 (p. 1094).** Orbits near the L2 Lyapunov orbit in Levi-Civita coordinates, and horseshoe orbits near
  the L1 and L2 Lyapunov orbits at C = 3.00001.
- **Figs. 14 and 15 (p. 1095).** mu = 0.01 families at C = 3.0199 and 3.01, and two horseshoe orbits at C = 2.9.
- **Fig. 16 (p. 1096).** Sun-Jupiter mass, C = 3.001906457: flows on the stable manifolds of the L3 Lyapunov orbit and the
  curves in the (x, x') plane, with horseshoe-shaped crescents.
- **Figs. 17 and 18 (p. 1098).** mu = 0.01 and 0.1: manifold flows, horseshoe-shaped at 0.01 and not at 0.1.

## Citation mining

Held: Szebehely 1967, Taylor 1981 (cited here) and Llibre, Martinez and Simo 1985. Schanzle 1967 is held but this paper
does not cite it.
On the wanted list: the Barrabes and Mikkola 2005 paper and Henrard 2002 in the same row as this paper.
Not held and not wanted before today (one new low-priority row):
- Dermott and Murray (1981), Icarus 48:1 and 48:12, the coorbital motion of 1980S1 and 1980S3. The paper compares with
  these. Also Dermott, Murray and Sinclair (1980), Nature 284:309.
- Yoder, Colombo, Synnott and Yoder (1983), Icarus 53:431, the Trojan-like analytic approximation applied to Janus
  and Epimetheus.
- Spirig and Waldvogel (1985) and Waldvogel and Spirig (1988), the singular-perturbation treatment of the Saturn coorbitals.
- Petit and Henon (1986), Icarus 66:536, satellite encounters in Hill's problem; Henon and Petit (1986), Celest. Mech. 38:67.
- Henrard (1983), Celest. Mech. 31:115, and Garfinkel (1977, 1978, 1980), formal horseshoe solutions.
Background: Rabe 1961 (AJ 66:500), Aksnes 1985, Cors and Hall (in preparation), Font 1990 thesis, Murray and Dermott 1999.
