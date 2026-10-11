# Digest: Jbara 2025, "Chaotic Dynamics and Zero-Velocity Structures in the Pluto-Charon CR3BP" (#960, wanted-list row 62)

A. W. Jbara (independent researcher), arXiv:2510.13479 [astro-ph.EP], v2 2025-12-10 ("Draft version December 11, 2025",
AASTeX 7.0.1). doi 10.48550/arXiv.2510.13479. The #1052 prereg note records v1 as submitted 2025-10-15 (arXiv abstract page). It is not peer reviewed
and no journal version is known.
- Proposed filename: `cyclers_pdf/papers/jbara-2025-chaotic-dynamics-zero-velocity-structures-pluto-charon-cr3bp-arxiv-2510.13479v2.pdf`.
  - arXiv PDF, 11 pp. (8 pp. of text and 3 pp. of appendices A-E), with a text layer, md5 5dda375cd9c2ec85c719a2a3f81ff56b
    (1,638,856 bytes).
- **How it was read:** the whole text layer. Pp. 2-6 (the equations, every figure with its title and caption, and every
  number) were read on 100 dpi page images. The figure titles carry mu = 0.1088; that value is not in the body text
  and is visible only on the images.
- Checks: `checks/check_pluto.py` (Jbara block) and `checks/check_jbara_mu.py` / `.out`. All values in these files are ours.

## 0. Verdict

**This is a short educational study of the planar circular restricted three-body problem (CR3BP) at the Pluto-Charon
mass ratio.** It plots zero-velocity curves (ZVCs) at the Jacobi constants of L4 and L1. It finds the Lagrange points
by Newton-Raphson and integrates a few trajectories with fixed-step RK4 (h = 5e-4).
- Two low-mu controls are shown: a tadpole and a horseshoe at mu = 0.001.
- At Pluto-Charon mu, the "tadpole" attempt escapes and hits Charon at T about 37. The "horseshoe" attempt breaks up
  after about 2 librations.
- Two runs start near L1: one stays near Pluto (neck closed) and one transits to Charon (neck open).
- It concludes that Pluto-Charon cannot hold long-lived Trojans. That is the known result that mu exceeds Routh's
  value (0.0385).
- Despite the title, it computes no chaos indicator: no Lyapunov exponent, FLI, MEGNO or Poincaré section. "Chaotic"
  means only "sensitive in a few runs".
- It has no periodic orbit, no table and no complete initial state. Each start is given only relative to L-points,
  with "tangential" or "vertical" velocity and a target C.
- It does not mention the small moons or their resonances.

**For the project (PROPOSALS only):**
- **#998 small-moon cyclers: none.** The paper's region and method do not overlap the #998 lane: it covers L-point
  neighbourhoods in the CR3BP, with no small moons, flybys or periodic closure. It is no prior publication of a #998
  type and gives no control.
- **Literature-gate anchor (#1052, row 5).** The prereg note already corrects the title and authors.
  - Its row 5 reads: "| 5 Tadpole (arXiv:2510.13479) | {binary-coorbital} unchanged | unset | arXiv abstract read
    2026-10-11: ... 'cannot host long-lived Trojan companions'. The label is right. The anchor's citation is wrong:
    its title string (...) is not the arXiv title, and its authors field is a placeholder ('Pluto-Charon', 'Charon').
    Proposed correction: authors ('Jbara',), the real title, provenance verified-against-source at abstract level.
    Not a periodic-orbit paper. |".
  - From the full text, one more point: the only working tadpole and horseshoe orbits are at mu = 0.001 (Figs. 3, 4).
    Every Pluto-Charon co-orbital attempt FAILS (Figs. 5, 6). The paper publishes no Pluto-Charon co-orbital orbit;
    it publishes the absence of a stable one.
  - The anchor's keyword "Pluto-Charon L4 L5 L3 periodic orbit" (`literature_check.py`) is not supported by the text.
  - Proposal: #1052 should decide whether a {binary-coorbital} label that makes a Pluto-Charon co-orbital signature
    "published" is right for a paper whose result is negative. One option is to keep the anchor as "studied, none
    stable" and drop the "periodic orbit" keyword. Provenance can rise to verified-against-source (full text).
    Do not edit `literature_check.py` here.
- **mu.** The body text gives mu ~ 0.109 (eq. 2) from masses 1.3e22 and 1.6e21 (kg implied, "S. Stern 2003"). Those
  masses give 0.1096 (ours). The figure titles say mu = 0.1088.
  - The printed Jacobi constants C_L4 = 2.9031 and C_L1 = 3.6210 hold together only for mu in [0.10874, 0.10878]
    (ours, `check_jbara_mu.out`). That interval holds the project's 0.1087647, the plu060 kernel's 0.1087740
    and Brozović-Jacobson's printed 0.1087759.
  - It excludes 0.109 (C_L4 2.9029, C_L1 3.6216) and 0.1096, and just excludes 0.1088 (C_L1 3.6211).
  - So the numbers were computed with a JPL-type mu of about 0.10877. The "0.109" and its cited masses are not what
    was used. The difference does not matter for any project use.
- Wanted list: row 62 is closed by this paper. No new row is proposed; the one possible lead (Sicardy 2010) is low
  priority (sec. 4).
- No catalogue row cites this paper. No catalogue change is proposed.

## 1. Content

**Model** (sec. 2, p.2; read on the image).
- Planar CR3BP, barycentric rotating frame. Unit separation, total mass and angular velocity. Pluto at (-mu, 0),
  Charon at (1 - mu, 0).
- Omega = (x^2 + y^2)/2 + (1 - mu)/r1 + mu/r2 (eq. 3), and C = 2 Omega - (xdot^2 + ydot^2) (eq. 5).
- Equations of motion: xddot = dOmega/dx + 2 ydot, yddot = dOmega/dy - 2 xdot (eqs. 9-10).
- This C convention omits the mu(1 - mu) term of Szebehely's. So C_L4 = 3 - mu + mu^2 (ours), not 3.
- Collinear points: from f(x) = 0 by Newton-Raphson, with the region-split forms in eqs. C20-C22. We checked the signs.
- Speed from a target C: |v| = sqrt(2 Omega - C_target) (eq. E29).
- A collision is r1 or r2 < 1e-6 (normalised units; not a body radius).
- Planar justification: Charon e "~ 0.005" (cited to Rhoden et al. 2015) and i < 1 deg. Brozović-Jacobson 2024 give a
  Plutocentric e of 1.61e-4, so the 0.005 in the paper is too large by a factor of 30. It does not affect the CR3BP runs.

**Numerical results, as printed** (all in normalised units; "drift in C" is the paper's energy-error figure):

| Fig. | System | Setup | Outcome | Drift in C |
|---|---|---|---|---|
| 1 | PC | ZVC at C = 2.9031 + 1e-3 (C_L4 = 2.9031) | small ovals round L4, L5 | - |
| 2a | PC | C = 3.6210 + 1e-2 (C_L1 = 3.6210) | L1 neck closed | - |
| 2b | PC | C = 3.6210 - 1e-2 | L1 neck open | - |
| 3 | mu = 0.001 | tadpole; C_L4 - 5e-4, tangential v, start (L4x - 0.04, L4y - 0.03), h 5e-4 | stable libration round L4 | 6.7e-15 |
| 4 | mu = 0.001 | horseshoe; C_L3 - 5e-5, vertical v, start (L3x - 0.027, L3y) | horseshoe | 4.9e-15 |
| 5 | PC, mu 0.1088 | as Fig. 3 | escapes; hits Charon at T ~ 37 | 5.1e-9 |
| 6 | PC, mu 0.1088 | horseshoe; C_L3 - 1e-3, vertical v, start (L3x - 0.027, L3y) | breaks up after ~2 librations | 8.9e-13 |
| 7a | PC, mu 0.1088 | C_L1 + 1e-2, vertical v, start (L1x - 0.04, L1y) | stays in the Pluto lobe | 4.9e-9 |
| 7b | PC, mu 0.1088 | C_L1 - 1e-2, same start | transits to Charon at T ~ 27-29 | 8.1e-8 |

- Integration lengths T are not printed, except the collision and transit times.
- Routh's critical value is quoted as 0.03852 (Sokol'skii 1975); ours (1 - sqrt(23/27))/2 = 0.0385209.
- Lagrange points and C at mu = 0.1088 (ours, for the record):
  - positions: L1 0.59251, L2 1.26260, L3 -1.04526;
  - Jacobi constants: C_L1 3.62111, C_L2 3.47974, C_L3 3.10828, C_L4 2.90304.
  - The paper prints no L-point positions and no C_L2 or C_L3.
- Sec. 2.6 calls tadpoles and horseshoes "well-known periodic solutions". In general they are quasi-periodic. Noted
  only.

**Discussion** (sec. 4). It states the model limits: planar, circular, a massless particle. It speculates that
Pluto-Charon manifolds may connect with Sun-Pluto ones ("interplanetary superhighways"). Nothing is computed for that.

## 2. Witnesses

- Witness C (decides): the page images of pp. 2-6. Every number in the table above was read there. The figure titles
  (mu = 0.1088, mu = 0.001) are only on the images.
- Witness A: the text layer. It agrees on every caption number. It loses the figure titles and the radicals in eqs.
  7-8 and A1-A3.

## 3. Cross-checks (ours)

- The C_L4 and C_L1 bracket on mu is given above. Each printed constant is reproduced to its 4 decimals for mu of about
  0.10877.
- Eqs. C20-C22 signs: correct for each region.
- Appendix B: the derivation of C is standard and correct.
- Fig. 2 caption (a) says "All L points are in the forbidden region". At C = C_L1 + 0.01 that holds by the C ordering
  (C_L1 is the largest). Consistent.
- **Citation defects, kept visible:**
  - "Stern, S. 2003, ARA&A 30, 185, doi 10.1146/annurev.aa.30.090192.001153". The volume and DOI are Stern 1992
    ("The Pluto-Charon system"), so the year is wrong.
  - "Llanos & Matas 2022" has no venue.
  - "Szebehely & Grebenikov 1969, Soviet Ast. 13, 364" is cited for the CR3BP definition. That attribution is
    unusual; not checked.
  - Sicardy 2010 (CMDA 107:145) is cited as "consistent with previous findings" that Pluto-Charon cannot host Trojans.
    Its Crossref title (2026-10-11) is "Stability of the triangular Lagrange points beyond Gascheau's value". We did
    not read it. From the title it studies when triangular points can be stable ABOVE the Routh/Gascheau limit, so
    "consistent with" may not describe it. Flagged, not judged.

## 4. Citation mining

12 references. Against `ls cyclers_pdf/papers | grep -i` and `CORPUS_INDEX.md`:

| Reference | Held? | Wanted list |
|---|---|---|
| Sicardy 2010, CMDA 107:145, doi 10.1007/s10569-010-9259-5 | NOT held | no row. Possible low-priority row for `#949` X4 (binary co-orbital background): triangular-point stability above Gascheau's value bears directly on heavy-secondary co-orbitals |
| Sokol'skii 1975, PMM 39:342 (Routh value) | NOT held | no row; textbook result, not needed |
| Calleja, Doedel, Humphries, Lemus-Rodríguez & Oldeman 2012, CMDA 114:77 | NOT held (two other Calleja papers are held; they are different works) | no row; CR3BP continuation background |
| Gómez, Koon, Lo et al. 2001, Adv. Astronaut. Sci. 109 (AAS 01-301) | NOT held as cited; the 2004 Nonlinearity 17 journal form (`gomez-koon-lo-marsden-masdemont-ross-2004-...`) is HELD | n/a |
| Murray & Dermott 1999, Solar System Dynamics | NOT held | no row; textbook |
| Buie et al. 2006, AJ 132:290 | NOT held | no row (also cited by Brozović-Jacobson) |
| Stern et al. 2015, Science 350:aad1815 | NOT held | no row |
| Stern 1992 (cited as 2003), ARA&A 30:185 | NOT held | no row |
| Rhoden et al. 2015, Icarus 246:11 | NOT held | no row; Jbara takes Charon's e ~ 0.005 from it (not checked against Rhoden); PLU060 gives 1.61e-4 |
| Clemence 1966; Szebehely & Grebenikov 1969; Llanos & Matas 2022 | NOT held | no row |

*Proposed filing: `cyclers_pdf/papers/jbara-2025-chaotic-dynamics-zero-velocity-structures-pluto-charon-cr3bp-arxiv-2510.13479v2.pdf`.*

*Filed as `cyclers_pdf/papers/jbara-2025-chaotic-dynamics-zero-velocity-structures-pluto-charon-cr3bp-arxiv-2510.13479v2.pdf`.*
