# Digest: Bruno 1994, *The Restricted 3-Body Problem: Plane Periodic Orbits*

- **Citation:** Bruno, A. D. (1994). *The Restricted 3-Body Problem: Plane Periodic Orbits.* de Gruyter Expositions in Mathematics 17. Walter de Gruyter, Berlin and New York. Preface by V. G. Szebehely. Translated from the Russian by B. Erdi. ISBN 3-11-013703-8. doi:10.1515/9783110901733. Russian original: *Ogranichennaya zadacha trekh tel: Ploskie periodicheskie orbity*, Nauka, Moscow, 1990. (The author's name is also spelled Brjuno / Bryuno.)
- **Pages:** 376 PDF pages. The book has xiv front pages and 362 printed pages. The book says it has 93 figures and 12 tables (p. iv). I count 13 numbered tables (see section 3).
- **PDF md5:** `f9e3773ac86ccdbdf0d049e9bdcbd67b` (`md5 -q`; 85 420 814 bytes, text layer present). Filed as `cyclers_pdf/papers/bruno-1994-restricted-3-body-problem-plane-periodic-orbits-de-gruyter-expositions-math-17-doi-10.1515-9783110901733.pdf`. Not compressed: the pages are 150 dpi anti-aliased greyscale, `ocrmypdf --optimize 2` saved 3 % and a lossy grey-JPEG rewrite still gave 62 MB, so the original is kept.
- **Page numbers:** for the main text, PDF page = printed page + 14 (printed p. 1 = PDF p. 15). In this digest "p. N" is the printed page. "PDF N" is the PDF page.
- **How I read it:**
  - Read in full (text layer): Introduction (pp. 1-9), Chapters III, IV, VI, VII, VIII, IX, and every chapter "Notes" section from Chapter III on.
  - My text filter kept prose lines and dropped most display-equation lines. So I did not read every equation. I checked these on the page image: (III.18)-(III.23), (III.40)-(III.42), the V^2 relation on p. 131, (IV.20), (VII.26)-(VII.27), (VII.45)-(VII.47), the Tr_1 and ϖ̇_1 formulas on p. 260, (VIII.50)-(VIII.51), and the Id/Ir segment definitions on p. 289. Other equation numbers below come from the text layer. The text layer garbles them (for example "(111.35)" for (III.35)), so treat them as "about here".
  - Read in part: Chapter V (opening, pp. 178-180, p. 186, pp. 190-193). Chapter X (all, but quickly).
  - Skimmed: Chapters I and II. I read the section headings and theorem statements, plus pp. 79, 84, 86-87 (natural families, trace cases, global section).
  - Appendix: all 8 tables viewed as page images at 150 dpi. Tables IV.1, IV.2, IV.3, V.1 and VII.1 also viewed as images.
  - Bibliography (pp. 347-356): read in full from the text layer.

---

## 0. Verdict

**What the book is.** It is a mathematical monograph on the planar circular restricted problem for small mass ratio μ. Its organising idea is Poincaré's: classify periodic orbits by their limits as μ → 0 ("generating solutions"). The book works out the μ = 0 problem completely in a frame centred on the big body P1 (not the barycentre). In that limit P2 has zero mass but is still a point that P3 can hit. So orbits of the two-body problem P1-P3 are cut at collisions with P2 into "arcs". The book:

1. lists all arc solutions that start and end at a collision with P2 (Hénon's symmetric families S, and Bruno's asymmetric families T_N);
2. studies where the Jacobi constant has extrema on these arcs (these are the candidates for stable second-species orbits);
3. finds the regular generating families: first kind (circular, Id and Ir), second kind (elliptic resonant, E_N^± and G_N), arcs that collide with P1 (K-families), and hyperbolic orbits;
4. gives first-order (O(μ)) perturbations of the trace and of the apse motion along these families.

It uses two tools throughout. One is the method of normal forms (Chapters I-II). The other is a "global section" Γ (z1 = 0, i.e. the apse points, r' = 0) and the symmetry plane Π, with Bruno's coordinates (ã, ẽ) on Π (eq. III.22, p. 108).

It is **not** a catalogue of families at finite μ. There are no initial conditions in rotating x, y, ẋ, ẏ at finite μ. All tables are μ = 0 generating data, or O(μ) coefficients.

**Coordinate warnings for anyone who uses the numbers.**
- The origin is at P1 (p. 95). Units: m1 + m2 = 1, P1P2 distance = 1, angular speed = 1. P2 sits at x1 = 1, x2 = 0.
- At μ = 0 the Jacobi constant is C = -2H = 2c + 1/a (p. 171), where c is the area (angular momentum) integral and a the semi-major axis of the P1-centred conic. At μ = 0, P1 is the barycentre, so this C is the standard Jacobi constant (the Tisserand form) with no conversion. Bruno's frame and the usual barycentric frame differ only at O(μ) (origin shift and indirect term). At a collision with P2 the synodic speed obeys V^2 = 3 - 2c - 1/a = 3 - C (p. 131, checked on image).
- On Π, ã = ±a and ẽ ∈ (-2, 2) with e = |1 - |ẽ||. ẽ > 0 is direct, ẽ < 0 is retrograde, |ẽ| > 1 is a pericentre point, |ẽ| < 1 an apocentre point (pp. 108-110, 123; image checked).
- Bruno's sign ε differs from Hénon's (1968) ε by the factor ε″ (p. 121).
- "Tr" is the trace of the 2x2 reduced monodromy. Stable means |Tr| < 2. On second-kind families Tr = 2 + μ Tr_1 + O(μ^2) (VII.26), so Tr_1 < 0 means stable for small μ.

**What it gives to #944 X2 (generating families and second-species orbits as μ → 0).** This is the book's main subject. It is directly useful.
- The arc families S (Hénon 1968 names A_0, A_1, ..., B_1, B_2, ..., C_12, C_23, C_24, ...) are the building blocks of second-species orbits. Their defining equations are (III.40)-(III.42) on p. 121 (image checked). Equation (III.42) links τ and η only. Equation (III.41) then gives a and e. Equation (IV.20) on p. 158 gives the closed form of all characteristics in Bruno's variables N^-1, θ, e*.
- **Table IV.1 (pp. 133-135)** gives 97 arc solutions with (ã, ẽ, τ/π, η/π). This is the best numeric control in the book for second-species generating arcs. I transcribed it in full (section 3). I checked all 97 rows against the three equations (III.40), with the signs ε, ε', ε″ taken from the row itself. 94 rows fit to better than 2e-4. 3 rows do not fit (rows 20, 62, 90); see the notes in section 3.
- Chapter IV section 3 and Chapter V give where C has extrema on the S families. Hitzl and Hénon (1977b) proved that stable second-species orbits with one close approach tend to these extremal points (p. 176). The book proves that Theorem IV.3.2 lists all of them (Chapter V). It does not tabulate them. It points to Hitzl and Hénon (1977a) for tables with τ/π, η/π < 10 (p. 193).
- The book does **not** build full multi-arc second-species periodic orbits or match them across collisions. It says this is "far from exhaustive" and open (p. 7, p. 176).

**What it gives to #948 R4 (Earth-grazing second-species cyclers built from e = 1 arcs).** Useful, but read the mapping carefully. The book has three different "e = 1" objects:
1. **Arcs that collide with P1** (rectilinear, c = 0, e = 1). These are the K-families of Chapter IX (KS_k^±, K_k^j, KI_k, and limits KL_k^j). Appendix Tables 5-8 give θ/π versus N^-1 (or Δ) for k = 1..5. These are relevant only if R4 treats the grazed body as P1.
2. **Parabolic collision orbits with P2** (a^-1 = 0, e = 1). This is family TL (pp. 219, 256-257), the limit of T_N as N = 1/p, p → ∞. Appendix Tables 3 and 4 give θ_cr versus c for TL and TL_{p+q}. Table 3 also gives the limit family GL with its O(μ) trace and apse-rate coefficients. These are relevant if R4's grazed body is P2 and the arcs are near-parabolic.
3. **Hyperbolic part of A_0** (Chapter VI, Fig. VI.11, p. 225). Hénon (1968, Table 1) tabulates it. The book has only a figure.

For R4 the most important lead in the book is a citation, not a table. **Hénon (1968), "Sur les orbites interplanétaires qui rencontrent deux fois la Terre"** is the source of all S families, and it is about orbits that meet the Earth twice. Bruno (1978a) = Celest. Mech. 24 (1981) 255-268, "On periodic flybys of the moon", is also directly on topic. Both are cited and **both are held** (Hénon 1968 Bull. Astron.; Bruno 1981 Celest. Mech. 24:255 = the book's "Bruno 1978a"); corrected by the corpus filer.

**Symmetric-family controls.**
- The book never uses Strömgren letters (a, b, c, f, g, h, i, l, m, n). Strömgren, Darwin and Breakwell are not cited anywhere (I searched body and bibliography). So the book gives **no** crosswalk to Strömgren families. See section 2.
- What it does give in its own notation:
  - **Id and Ir** (direct and retrograde circular orbits about P1): closed form (III.18)-(III.19), p. 107 (image checked). x1 + i x2 = a exp[i(n - 1)t], y1 + i y2 = i n a exp[i(n - 1)t], with n = ±a^(-3/2) (+ for Id, - for Ir). Period T = 2π/|n - 1|, trace Tr = 2 cos T. **There is no numeric table for Id or Ir.** The closed form is the control.
  - **Ir+ and Ir-** (p. 289, image checked): Ir is split only by the S families, into Ir+ (a < 1) and Ir- (a > 1). This is where the label "IR+" comes from.
  - **Id_p^+ and Id_p^-** (p. 289, image checked): Id is split at the first-order resonances N = (p+1)/p and N = (p-1)/p (VIII.50)-(VIII.51), p. 286. Id_p^+ is the segment p/(p-1) > N > (p+1)/p > 1. Id_p^- is the segment 1 ≥ p/(p+1) > N > (p-1)/p. So "ID1" means Id_1^+ (N > 2, i.e. a < 2^(-2/3) ≈ 0.630) or Id_1^- (N < 1/2). Fig. VIII.15 (p. 289) shows how the E_N^± segments join these pieces.
  - **A_0, A_1, B_1, ...**: Hénon's symmetric arc families (p. 121; Fig. III.17 on p. 122 is Hénon's diagram; Fig. III.18 on pp. 124-125 and Fig. IV.12 on pp. 150-151 show them on Π). Table IV.1 has sample points on A_0, A_1, A_2, B_1, B_2, B_3, C_12, C_23, C_24, C_34, C_35.
  - **E_N^+ (θ = 0) and E_N^- (θ = π/(p+q))**: symmetric second-kind families, N = (p+q)/p (p. 244-245). **G_N**: asymmetric second-kind families. They exist only for p + q = 1, ε' = +1, one per p (hypothesis on p. 252, from computation).
  - **L1-L5**: the Lagrange points. At μ = 0, L3 lies on E_1^- (θ = π), and L4, L5 lie on G_1 (θ = ±π/3) (p. 254). Families ℒ3 ⊂ E_1^- and ℒ4, ℒ5 ⊂ G_1 emanate from them (pp. 293-295).

**Numbers usable as test controls (page = printed page).**

| What | Where | Use |
|---|---|---|
| Closed form of Id, Ir: T = 2π/\|n-1\|, Tr = 2cos T | eq. (III.18)-(III.19), p. 107 | μ → 0 limit of first-kind families, incl. Ir+ (h?) |
| Table IV.1: 97 symmetric arcs (ã, ẽ, τ/π, η/π) | pp. 133-135 | second-species generating arcs; check (III.40)-(III.42) |
| Table VII.1: e and c where G_{1/p} meets E_{1/p}, p = 1..20 and ∞, 9-13 digits | p. 252 | high-precision generating-family bifurcation points |
| Appendix Table 1: G_{1/p}, p = 1..5 (e, c, θ_cr, θ, Tr_1, ϖ̇_1) | pp. 328-330 | asymmetric second-kind generating family; O(μ) trace |
| Appendix Table 2: T_N and E_N^± (ε'e, c, θ_cr, Tr_1^±, ϖ̇_1^±), 8 values of N | pp. 331-334 | symmetric second-kind generating families; O(μ) stability |
| Appendix Table 3: TL and GL (c, θ_cr, θ, ⟨Tr_1⟩, ⟨ϖ̇_1⟩) | p. 335 | parabolic (e = 1) collision limit; R4 lead |
| Appendix Table 4: TL_{p+q}, EL^± for p+q = 1..4 | pp. 336-339 | parabolic limit of E_N |
| Appendix Tables 5-8: K-families (arcs that collide with P1) | pp. 340-345 | e = 1 rectilinear generating arcs |
| GL at θ = π: c1 = 2.4771498, c2 = 0.29500129 (d1 = 3.0681355, d2 = 0.04351288) | p. 262 | endpoint of GL; equals Table VII.1 row p = ∞ |
| Stationary-point families: Δ = -(21/8)μ on ℒ3, (27/4)μ on ℒ4, ℒ5 | p. 294 | sign check of stability near L3, L4, L5 |

**Bottom line.** For #944 X2 this book is the main theory reference and gives one good arc table (IV.1). For #948 R4 it gives the parabolic family TL / GL tables and points to Hénon 1968 and Bruno 1978a/1981 (Celest. Mech. 24). For the Strömgren-letter controls (a ... n) it gives nothing directly: use Hénon 1997 / Bruno-Varin 2006 for that mapping, not this book.

---

## 1. Chapter by chapter

### Preface (Szebehely), Basic notation, Introduction (pp. v-xiv, 1-9)
- Szebehely says the book is the "volume on periodic orbits" he planned after *Theory of Orbits* ch. 8.
- Basic notation (pp. xi-xiv) defines all family names: Id, Ir, S, T_N, A_l, B_l, C_jk, E_N^±, G_N, EL, GL, Q_N, KS_k^±, KJ_k, KI_k, KL_k^j, and L_j. Note ε' = sgn c (direction), ε″ = +1 at pericentre and -1 at apocentre.
- The Introduction states the program. As μ → 0 there are three limit problems: (1) two bodies P1-P3; (2) Hill's problem; (3) two bodies P2-P3 (Bruno 1978a). The book treats only problem 1.
- Generating solutions can be regular or singular. Poincaré's second-species orbits have singular generating solutions (p. 5). The study of these is "far from exhaustive" (p. 7).

### Chapter I. The normal form of a Hamiltonian system (pp. 13-53) [skimmed]
- Normal form near a stationary point: linear (complex and real forms), then nonlinear normalization with resonances (Theorems 2.1-2.4).
- Lowering the number of degrees of freedom at resonance (Theorem 3.1).
- The analytic sets A (periodic) and B (conditionally periodic) on which the normalizing map converges (Theorems 4.1-4.3).
- Generalizations: small parameter, a manifold of stationary points (Theorems 5.1-5.2). These are used later for the L-point circle.
- No tables.

### Chapter II. Normalization near a cycle or a torus (pp. 54-92) [skimmed]
- Normal form of linear Hamiltonian systems with periodic coefficients; nonlinear normal form near a periodic solution (Theorems 2.1-2.2).
- For two degrees of freedom (§4.2, pp. 84-86): five cases of Tr. Tr > 2: unstable family. Tr = 2: the family has an extremum of H at that orbit. Third-order resonance (Tr = -1) and Tr = 2: unstable. Order 2 or 4 (Tr = -2 or 0): either.
- §3.3 (p. 79) defines "natural families" (a family continued until it ends). §4.4 (pp. 86-87) defines the global section Γ and "characteristics" (the curves where a family cuts Γ).
- No tables.

### Chapter III. Periodic solutions and arc solutions (pp. 95-132)
- Sets up the problem in P1-centred synodic coordinates (III.1)-(III.5). H = H0 + μR. At μ = 0, P2 is a massless point that can still be hit.
- Two-body facts: N = |n| = a^(-3/2) (III.10). Rational N = (p+q)/p (III.14). Synodic orbits for N = 3, 2, 5/3, 3/2, 4/3, 1, 1/2 are drawn in Figs. III.1-7 (pp. 100-106). A small unnumbered table on p. 99 gives the e of each drawn orbit (7/9, 1/3, 1/9, 5/9, 1).
- Circular families: (III.18)-(III.19), p. 107 (see section 0). For rational N the torus D_N^c is filled with periodic orbits of period 2πp and Tr = 2.
- Plane Π and coordinates ã, ẽ (III.22)-(III.23), p. 108. Global section Γ: z1 = 0 (III.25)-(III.26), with cylindrical coordinates (ã, θ, ẽ), p. 111-112.
- Arcs with consecutive collisions with P2 (§3). Collision equations (III.35). Two kinds: T_N (asymmetric, rational N, collision twice at the same sideral point) and S (symmetric; Hénon 1968). Equations (III.40)-(III.42), p. 121. Four domains ω1-ω4 on Π (Fig. III.18). Intersection types I, II, II*, III, III*, IV* of characteristics (§3.5). Collision velocity V^2 = 3 - C and exit angle α' = π - α for S arcs (p. 131).

### Chapter IV. Properties of arc solutions of the families S and T_N (pp. 133-177)
- §1: orbit drawings of A_0-A_2, B_1-B_3, C_12, C_23, C_24, C_34, C_35 (Figs. IV.1-11, pp. 136-149). **Table IV.1** (pp. 133-135) lists the drawn orbits. Fig. IV.12 (pp. 150-151) shows the characteristics on Π with arrows at extremal points.
- §2: new variables N^-1 = a^(3/2) and e* (IV.15). Theorem 2.1 (IV.20), p. 158: e* = ±cos{(θ - π(N^-1 - 1)k)/2}, ε″ sgn(a - 1) = (-1)^k. **Table IV.2** (p. 158) tabulates the singular curve f (IV.17).
- Theorem 2.5 (p. 163): at intersections of type I, e is transcendental. So no two such points share an algebraic e.
- §3: Jacobi constant at μ = 0. On ω', -2√2 < C < 3 (p. 168). Fig. IV.17 shows level lines. **Table IV.3** (p. 172) gives C ranges and sgn(dC/da) at type-II points. Theorem 3.2 (p. 171) lists the extremal points. Schanuel's conjecture would imply that no two "singular points" share the same C (Theorem 3.3, p. 176).
- Notes: Hénon and Guyot (1970): on second-species generating orbits Tr = ±∞, and its sign changes at extrema of C. Hitzl and Hénon (1977b): stable second-species orbits with one passage tend to extremal points (p. 176).

### Chapter V. Extrema of the Hamiltonian on families of arc solutions (pp. 178-193) [read in part]
- Main theorem: Theorem IV.3.2 lists all extremal points on the S characteristics.
- Proof uses a function G (V.3) whose level lines G = -kπ give the extremal points. G has one saddle, G_cr = 0.0153860 (p. 186). G_1 has min 0.327365 at r = 0.912871 (p. 179).
- **Table V.1** (p. 186): six segment types; sign of -κ sgn(a-1); number of extrema (even/odd, lower bound); range of k. Signs and parity only, no values.
- Result: types I and V have two extrema, III and IV have one, II and VI have none (p. 193). A computed root on P = 0: c = -0.617025016, a = 0.853408456, e = 0.744233279 (p. 191).
- Numeric support is in Bruno (1978b, Tables 1-5), which is not in this book.

### Chapter VI. Trajectories of collisions (pp. 194-232)
- Extends Chapter III to all orbits: elliptic, parabolic, hyperbolic, and rectilinear. Universal two-body solution (VI.7) by a Levi-Civita-type parameter (pp. 196-198).
- Classification (pp. 198-199): I (circles: Id, Ir), D (ellipses), P (parabolas), H (hyperbolas), C1 (collision with P1, c = 0).
- §2: all orbits that hit P2. System (VI.22). Families V^{k,κ}, W_d, W_r (Theorem 2.1, p. 212). The set C2 is a 2-parameter family, a double cover of the domain ω (Theorem 2.2, p. 219). At a^-1 = 0 the family V^{0,κ} becomes TL (p. 219).
- §3: Theorem 3.1 (p. 222): every arc in C2C2 belongs to T_N or to S (Hénon's result, reproved). A sign table on p. 226 gives sgn v1 by segment type. Families u^{l,λ} of orbits that hit both P1 and P2 (pp. 226-228, Fig. VI.12 on p. 227).
- §4: critical solutions = those that hit or come arbitrarily close to P2. Parabolic orbits are non-localizable. List of noncritical localizable solutions (pp. 229-230) leads into Chapters VII-X. Notes (pp. 231-232) restate Poincaré's first kind / second kind / second species.

### Chapter VII. Solutions of the second kind (pp. 235-265)
- §1: general averaging theory for rational N. Generating solutions satisfy ∂[R]/∂θ = 0 (VII.20). Trace Tr = Σ μ^k Tr_k, Tr_0 = 2, Tr_1 = -f″(0) T_0^2 N^-2 ∂^2[R]/∂θ^2 (VII.26); period T_1 (VII.27) (p. 241, image checked).
- §2: in Delaunay elements the symmetry gives two symmetric families for every N: E_N^+ (θ = 0) and E_N^- (θ = π/(p+q)) (p. 244-245). Other zeros are G_N.
- §3: R = R1 + R2 + R3 (direct, indirect). [R2] via Bessel functions (Theorem 3.1). (VII.45)-(VII.47), p. 249: [R] = const + A e^s {cos(p+q)θ + O(e^2)}, s = |q| (direct) or |q'| (retrograde).
- §4-5: computation (program "Pertur", Bruno 1976). An unnumbered table on p. 249 gives p_max for each p+q = 1..11. Result: G_N exist only for p + q = 1 and ε' = +1 (Hypothesis, p. 252). **Table VII.1** (p. 252) gives the G_{1/p} ∩ E_{1/p} points. Appendix Tables 1-2 give perturbations. The E_N^± that are unstable for small μ are the branches that turn toward P2 (Fig. VII.3, p. 253). This explains why asteroids avoid close approaches to Jupiter.
- §6-8: L3 on E_1^-, L4/L5 on G_1 (p. 254). §7: the limit p → ∞: parabolic families TL, GL, EL^± (Hénon's suggestion, letter of 16 Oct 1979). Appendix Tables 3-4. §8: conditionally periodic second-kind tori Q_N.
- Notes (pp. 263-265): history of errors on E_N stability (Schwarzschild, Andoyer, Message, Frangakis). Uno (1937) gave the first correct existence proof; Arenstorf (1963) and Barrar (1965) gave later ones.

### Chapter VIII. Solutions of the first kind (pp. 266-297)
- §1-2: normal form near a periodic solution with frequency ratio λ = p/q. Bifurcation pictures for q = 1, 2, 3, 4, ≥ 5 (Figs. VIII.3-13). Stability changes where the characteristic has a vertical tangent.
- §3: application. Uses Poincaré elements. First-order resonances (q̃ = 1) only for ε' = +1, N = (p+1)/p and N = (p-1)/p (VIII.50)-(VIII.51), p. 286. There Id breaks and joins E_N^±. For q̃ = 2 an instability interval J appears on Id (p. 287). q̃ = 3 and ≥ 4 cases include Ir with ε' = -1 (pp. 287-288).
- §4 (p. 288-290): definitions of Id_p^±, Ir^± and Fig. VIII.15 (see section 0 and section 2). For N > 2√2 (a < 1/2) the E_N^± do not meet the S families, so E_N is a closed natural family.
- §5: the circle of stationary points a = 1, e = 0. Three generating points: L3 (θ = π), L4, L5 (θ = ±π/3). Δ = -(21/8)μ on ℒ3 (unstable); +(27/4)μ on ℒ4, ℒ5 (stable) (p. 294). Model system (VIII.66): tadpole and horseshoe orbits (Fig. VIII.16).

### Chapter IX. Generating arc solutions from C1C1 (pp. 298-316)
- §1: uses coordinates centred at the P1-P2 barycentre shift (IX.1) and modified Delaunay elements. Proves that E_N are generating even when the orbit passes through P1 (Theorem 1.3, regularized by eccentric anomaly). Theorems 1.4-1.5: generating arcs that begin and end at P1 (e = 1, c = 0).
- §2: condition for an arc of k elementary pieces to be generating is the integral (IX.24) = 0. Symmetric families KS_k^+ (θ = 0) and KS_k^- (θ = π). Asymmetric K_k^j, K_k^0, KI_k. Critical families u^{l,λ} (dashed in Fig. IX.1). **Appendix Tables 5-6.**
- §3: limits as N^-1 → ∞ (parabolic arcs): families KL_k^j. Hypothesis: exactly k - 1 limiting families KL_k^j for k > 1, period 2 in Δ (p. 313). **Appendix Tables 7-8.**
- §4: open question of the order of the collision orbits along E_N(μ). It cites Broucke (1968, Fig. 22) as an example.
- Notes: K-families are named after Krasinskii (1963). KS_k^± were found by Levi-Civita (1903) and Moulton (1913) ("closed orbits of ejection").

### Chapter X. Perturbations of hyperbolic orbits (pp. 317-326) [read quickly]
- Theorem 1.1: a normal-form theorem for systems whose perturbation decays like |t|^-2.
- Result: in coordinates (a, θ, c, τ) every hyperbolic orbit that avoids P1 and P2 is generating. In Cartesian coordinates none is a uniform limit (p. 324-325).
- Orbits with equal in and out asymptotes satisfy two equations (X.31) and form one-parameter families. The book doubts they are useful.
- No tables.

### Appendix (pp. 327-345)
- Eight tables. See section 3 for the full list. The description on p. 327 says that more precise versions are in Bruno (1976, Tables 3-4), Bruno (1980a, Tables 1-2) and Bruno (1981, Tables 1-4).

---

## 2. Notation crosswalk

**The book gives no Strömgren, Darwin or Broucke letter names.** I searched the full body text and the bibliography. "Strömgren", "Darwin", "Breakwell" and "Copenhagen" do not appear. Broucke is cited three times (1962, 1968), but only for figures and for agreement of numbers. He is not cited for family names. Szebehely (1967) is cited for regularization and for the Lagrange points, not for family letters. So any a/b/c/f/g/h/i/l/m/n mapping must come from another source.

What the book does give:

| Bruno's name | Meaning | Other names the book gives | Page |
|---|---|---|---|
| A_0, A_1, A_2, ...; B_1, B_2, ...; C_12, C_23, C_24, ... | symmetric arcs P2 → P2 at μ = 0 (families S) | Hénon (1968) names, from his Fig. 4 (= Fig. III.17). Bruno: "Henon named A0, A1, A2, ...; B1, B2, ...; C12, C23, C24, ... we denote by the same symbol" | p. 121, Fig. III.17 p. 122 |
| ε, ε', ε″ | sign triplet of an S arc | "ε' and ε″ are analogous to ε' and ε″ in Hénon (1968) and Guillaume (1970), but our ε differs from the ε of Hénon by the factor ε″" | p. 121 |
| S | the set of symmetric arcs | Hénon (1968) "named S" | p. 115 |
| Id, Ir | circular direct / retrograde | Poincaré "first kind" | pp. 231, 290, 295 |
| E_N^±, G_N | elliptic resonant | Poincaré "second kind"; "Schwarzschild solutions" | pp. 231, 264, 290 |
| arcs of S and T_N | | Poincaré "second species" (generating) | pp. 231-232, 290 |
| Id_p^+ | segment p/(p-1) > N > (p+1)/p > 1 | ("ID1" = Id_1^±) | p. 289 |
| Id_p^- | segment 1 ≥ p/(p+1) > N > (p-1)/p | | p. 289 |
| Ir^+, Ir^- | Ir with a < 1 / a > 1 | ("IR+") | p. 289 |
| L1-L5 | Lagrange points | L3 ∈ E_1^-, L4, L5 ∈ G_1 at μ = 0 | p. xiv, p. 254 |
| ℒ3, ℒ4, ℒ5 | families from L3, L4, L5 | ℒ3 ⊂ E_1^-(μ); ℒ4, ℒ5 ⊂ G_1 (Liapunov short-period families) | pp. 293-295 |
| K-families (KS, KJ, KI, KL) | arcs P1 → P1 | named after Krasinskii (1963) | p. 316 |
| KS_k^± | symmetric P1-collision arcs | "found by Levi-Civita (1903) and by Moulton (1913)" | p. 316 |

Quotes for the key lines (from the text, checked on image):
- p. 121: "These solutions form infinitely many analytic curves, which Hénon named A0, A1, A2, ...; B1, B2, ...; C12, C23, C24, .... To each curve, in the phase space there corresponds its own family of symmetric arc solutions that we denote by the same symbol as the curve."
- p. 289: "We shall denote the segments of the family Id for p/(p-1) > N > (p+1)/p > 1 by Id_p^+, while the segments for 1 ≥ p/(p+1) > N > (p-1)/p are denoted by Id_p^-. The family Ir is divided only by the families S into altogether two parts: Ir^+ (for a < 1) and Ir^- (for a > 1)."

**Not from this book (unverified, for orientation only).** In later work (Hénon 1997; Bruno and Varin 2006, "Family h ...") the Strömgren family h is linked to the generating family Ir^+. I have not checked this against those sources. Do not cite this book for it.

---

## 3. Tables usable as controls

### 3.1 List of all numbered tables

| Table | Page (PDF) | Caption (short) | Columns | Rows |
|---|---|---|---|---|
| IV.1 | pp. 133-135 (PDF 147-149) | Elements of orbits shown in Figs. IV.1-11 | Number, Family, ã, ẽ, τ/π, η/π | 97 |
| IV.2 | p. 158 (PDF 172) | Values on the curve f defined by (IV.17) | a, e, N^-1, ẽ*, x, φ(x) | 27 |
| IV.3 | p. 172 (PDF 186) | Values of C and signs of dC/da | sgn(a-1), ε', C (inequality), sgn(dC/da) | 4 (inequalities and signs only) |
| V.1 | p. 186 (PDF 200) | Properties of segments of characteristics | Type, sign of -κ sgn(a-1), number of extrema, values of k | 6 (signs and parity only) |
| VII.1 | p. 252 (PDF 266) | e and c at intersections of G_N and E_N, N = 1/p | p, e1, e2, c1, c2 | 21 (p = 1..20, ∞) |
| App. 1 | pp. 328-330 (PDF 342-344) | The families G_{1/p} and their perturbations | e, c, θ_cr, θ, Tr_1, ϖ̇_1 | 109 (p = 1: 24, p = 2: 24, p = 3: 22, p = 4: 20, p = 5: 19) |
| App. 2 | pp. 331-334 (PDF 345-348) | The families T_N and perturbations of E_N^+ and E_N^- | ε'e, c, θ_cr, Tr_1^+, ϖ̇_1^+, Tr_1^-, ϖ̇_1^- | 144 (8 blocks × 18: N = 1, 1/2, 1/3, 1/4, 1/5, 2, 2/3, 2/5; ε'e = ±0.1..±0.9) |
| App. 3 | p. 335 (PDF 349) | The families TL, GL, and perturbations of GL | c, θ_cr, θ, ⟨Tr_1⟩, ⟨ϖ̇_1⟩ | 22 (c = 2.4 .. 0.3) |
| App. 4 | pp. 336-339 (PDF 350-353) | The families TL_{p+q} and perturbations of EL_{p+q} | c, θ_cr, ⟨Tr_1⟩^+, ⟨ϖ̇_1⟩^+, ⟨Tr_1⟩^-, ⟨ϖ̇_1⟩^- | 120 (4 blocks × 30: p+q = 1..4; c = 3.0 .. -3.0) |
| App. 5 | pp. 340-341 (PDF 354-355) | The families K_1, K_2, K_3 | N^-1, ψ̃/π, then θ/π on K_1^0/K_1^1, K_2^1, K_3^1, K_3^2, K_2^0, K_3^0 | 50 (N^-1 = 0.03 .. 1.99) |
| App. 6 | pp. 342-343 (PDF 356-357) | The families K_4, K_5 | N^-1, then θ/π on KI_4/K_4^1, K_4^0/K_4^2, K_4^3, K_5^0/K_5^1, K_5^2, K_5^2/K_5^3, KI_5/K_5^4 | 50 |
| App. 7 | p. 344 (PDF 358) | The families KL_1, KL_2, KL_5 | Δ, KL_1, KL_2, KL_5^1..KL_5^4 (θ/π) | 25 (Δ = 0.03 .. 0.99) |
| App. 8 | p. 345 (PDF 359) | The families KL_3, KL_4 | Δ, KL_3^1, KL_3^2, KL_4^1, KL_4^2, KL_4^3 (θ/π) | 25 |

That is 13 numbered tables. The book's own count is 12 (p. iv). Small unnumbered tables: p. 99 (orbit label → e for Figs. III.1-7), p. 123-126 (parity of p, q → domain ω_i of C_{p,p+q}; announced at the foot of p. 123), p. 226 (segment type → sign of v1), p. 249 (p+q → p_max).

**Text-layer warnings.** The text layer drops leading "1"s in the appendix (".21535" is 1.21535 in Table 1). It drops the "p = 2..5" block labels in Table 1. It garbles the Table 2 and 4 headers and the Table 8 title ("TL3" is KL_3). In Table IV.1 it drops leading "1"s (".95339" is 1.95339) and misreads family subscripts (rows 65-78 are B_3, not B_1). In Table 1, block p = 2, row e = 0.81 (p. 329), θ_cr is printed 0.41427 on the page image too. Its neighbours are 1.25969 and 1.57895, so this is likely a misprint for 1.41427 in the book itself. Always take values from the page image.

### 3.2 Transcription 1: Table IV.1 (second-species generating arcs, μ = 0)

Source: pp. 133-135 (PDF 147-149). Read from the 150 dpi page images, row by row. The text layer agreed with the image except for the dropped leading "1"s and the family labels, which I took from the image.

Meaning of the columns: ã, ẽ are the coordinates of the arc's symmetric point on Π (III.22): a = |ã|, e = |1 - |ẽ||, ã < 0 means the point is at θ = π (P3 opposite P2), ẽ < 0 means retrograde. 2τ is the time between the two collisions with P2, and 2η is the change of eccentric anomaly (p. 121). Units: P1 at origin, P1P2 = 1, angular speed 1.

**Check I ran (a transcription check, not an independent golden).** For each row I took the signs from the row: ε = sgn ã (p. xiii: ã = εa), ε' = sgn ẽ, ε″ = sgn(|ẽ| - 1) (ẽ = ε'(1 + ε″e), p. 201 text layer; pericentre if |ẽ| > 1, p. 123). Then I tested all three equations of (III.40), p. 121, with no absolute values and no best-fit choice: cos τ = εa(cos η - ε″e), sin τ = εε'a√(1 - e^2) sin η, τ = a^(3/2)(η - ε″e sin η). Logs: `t41_strict.log` (signed check) and `t41_check.log` (an earlier unsigned check) in this folder.
- 94 of 97 rows fit all three equations; every residual is below 2e-4 (the largest is row 5, where τ/π = 1.94 has only 2 decimals). This includes the 17 rows with integer η/π.
- 3 rows do not fit. I re-read them on the image and the print is as given:
  - Row 20 (A_1): residuals up to 0.9. In the unsigned check, the printed (τ/π, η/π) = (1.78, 1.16316) fits a = 1.36584, e = 0.30735, not the printed ã, ẽ. One pair is likely a misprint. Not resolved.
  - Row 62 (B_2): fails as printed. With η/π = 1.90883 instead of the printed 2.90883, all three equations fit to 6e-6. So this is very likely a misprint of the leading digit.
  - Row 90 (C_23): residuals about 0.02. Not resolved.
  - The corpus filer re-read rows 20, 62 and 90 on the page images (PDF 148-149). The transcription matches the print in all three, so the misfits are in the book, not in our copy.

| No. | Family | ã | ẽ | τ/π | η/π |
|---|---|---|---|---|---|
| 1 | A0 | -1.37126 | -1.87614 | 0.2164 | 0.4 |
| 2 | A0 | -1 | -1 | 0.5 | 0.5 |
| 3 | A0 | -1.06131 | -0.04349 | 0.9 | 0.51924 |
| 5 | A0 | -1.58265 | 0.62882 | 1.94 | 0.95930 |
| 6 | A0 | -1.58740 | 0.62996 | 2 | 1 |
| 7 | A0 | -1.59392 | 0.62198 | 2.08 | 1.05390 |
| 8 | A0 | -1.94110 | 0.02417 | 2.88 | 1.33449 |
| 9 | A0 | -2.06961 | -0.03887 | 3.16 | 1.31929 |
| 10 | A0 | -2.40018 | -0.39725 | 3.84 | 1.08094 |
| 11 | A0 | -2.51984 | -0.39685 | 4 | 1 |
| 13 | A1 | -2.78148 | 0.03366 | 4.82 | 1.26937 |
| 14 | A1 | -2.51984 | 0.39685 | 4 | 1 |
| 16 | A1 | -2.09757 | -0.02948 | 2.86 | 0.68125 |
| 18 | A1 | -1.58740 | -0.62996 | 2 | 1 |
| 20 | A1 | -1.01085 | -1.77607 | 1.78 | 1.16316 |
| 21 | A1 | -1.22983 | 1.54507 | 2.42 | 1.61139 |
| 22 | A1 | -1.31037 | 1.23686 | 3 | 2 |
| 23 | A1 | -1.38149 | 1.57281 | 3.54 | 2.33988 |
| 25 | A1 | -1.62785 | -1.70211 | 4.42 | 2.31488 |
| 26 | A1 | -1.84202 | -1.45712 | 5 | 2 |
| 27 | A2 | -1.84202 | 1.45712 | 5 | 2 |
| 28 | A2 | -1.75038 | 1.76829 | 4.38 | 1.68842 |
| 29 | A2 | -1.59546 | -1.94443 | 3.84 | 1.62932 |
| 30 | A2 | -1.43400 | -1.35184 | 3.24 | 1.82965 |
| 32 | A2 | -1.01887 | -1.17035 | 2.48 | 2.46532 |
| 34 | A2 | -1.11216 | 0.25002 | 3.26 | 2.54293 |
| 35 | A2 | -1.21141 | 0.82548 | 4 | 3 |
| 36 | A2 | -1.21879 | 0.77252 | 4.26 | 3.21053 |
| 37 | A2 | -1.33256 | 0.09164 | 4.82 | 3.41141 |
| 38 | A2 | -1.39895 | -0.06771 | 5.16 | 3.40105 |
| 41 | B1 | 2.08008 | 0.48075 | 3 | 1 |
| 42 | B1 | 1.95711 | 0.28446 | 2.48 | 0.73952 |
| 44 | B1 | 0.72634 | -1.54399 | 0.85511 | 1.25647 |
| 45 | B1 | 0.85656 | 1.54216 | 1.24 | 1.40004 |
| 46 | B1 | 1.29462 | 0.21700 | 1.72 | 1.40613 |
| 47 | B1 | 1.57091 | -0.05684 | 2.16 | 1.37409 |
| 48 | B1 | 1.82052 | -0.45988 | 2.68 | 1.18578 |
| 49 | B1 | 2.08008 | -0.48075 | 3 | 1 |
| 51 | B2 | 1.95339 | -1.57480 | 5.68 | 2.17714 |
| 52 | B2 | 1.83933 | -1.94297 | 5.18 | 2.33921 |
| 53 | B2 | 1.66394 | 1.69049 | 4.56 | 2.30388 |
| 54 | B2 | 1.51172 | 1.64708 | 3.44 | 1.67523 |
| 55 | B2 | 1.32171 | -1.92590 | 2.84 | 1.58467 |
| 56 | B2 | 1.20558 | -1.27253 | 2.36 | 1.71519 |
| 57 | B2 | 0.81355 | -0.57229 | 1.78677 | 2.32 |
| 58 | B2 | 0.88772 | 0.31730 | 2.22 | 2.44068 |
| 59 | B2 | 1.15655 | 1.76811 | 2.74 | 2.44361 |
| 60 | B2 | 1.29870 | -1.92374 | 3.16 | 2.41990 |
| 61 | B2 | 1.38827 | -1.46565 | 3.56 | 2.29492 |
| 62 | B2 | 1.65419 | -1.41227 | 4.14 | 2.90883 |
| 63 | B2 | 1.84468 | -1.94319 | 4.82 | 1.66135 |
| 64 | B2 | 1.99766 | 1.78826 | 5.4 | 1.71841 |
| 65 | B3 | 1.50543 | 0.12991 | 5.76 | 3.37390 |
| 66 | B3 | 1.40572 | 0.71138 | 5 | 3 |
| 67 | B3 | 1.38301 | 0.57374 | 4.6 | 2.72510 |
| 68 | B3 | 1.28840 | 0.09683 | 4.18 | 2.57972 |
| 69 | B3 | 1.22289 | -0.06513 | 3.86 | 2.56246 |
| 70 | B3 | 1.12893 | -0.82437 | 3.32 | 2.72534 |
| 71 | B3 | 0.86207 | -1.33211 | 2.74754 | 3.34 |
| 72 | B3 | 0.92722 | 1.60949 | 3.26 | 3.45889 |
| 73 | B3 | 1.12977 | 0.12185 | 3.82 | 3.45825 |
| 74 | B3 | 1.20236 | -0.08686 | 4.16 | 3.44099 |
| 75 | B3 | 1.25279 | -0.57516 | 4.52 | 3.34246 |
| 77 | B3 | 1.40572 | -0.71138 | 5 | 3 |
| 78 | B3 | 1.59279 | -0.05570 | 5.84 | 2.62895 |
| 79 | C12 | -0.76848 | 0.09495 | 0.90027 | 1.60802 |
| 80 | C12 | -0.94747 | -0.10710 | 1.14 | 1.51977 |
| 81 | C12 | -0.98748 | -0.68029 | 1.38453 | 1.51262 |
| 82 | C12 | -0.62996 | -0.41260 | 1 | 2 |
| 83 | C12 | -0.55970 | -0.10229 | 0.96210 | 2.16 |
| 84 | C12 | -0.60100 | 0.24240 | 1.06052 | 2.16 |
| 85 | C12 | -0.62996 | 0.41260 | 1 | 2 |
| 86 | C23 | -0.75484 | 1.37062 | 2.10967 | 3.16 |
| 87 | C23 | -0.66966 | -1.75109 | 1.89181 | 3.27192 |
| 88 | C23 | -0.76314 | -1.31037 | 2 | 3 |
| 89 | C23 | -0.87804 | -1.20291 | 2.21562 | 2.74 |
| 90 | C23 | -0.80487 | 1.55520 | 1.79428 | 2.63384 |
| 91 | C24 | 0.61917 | 0.32024 | 2.06194 | 4.14 |
| 92 | C24 | 0.58915 | -0.19188 | 1.94389 | 4.16863 |
| 93 | C24 | 0.72013 | -0.30183 | 2.14089 | 3.68791 |
| 94 | C24 | 0.65945 | 0.22368 | 1.89955 | 3.73166 |
| 95 | C34 | -0.80197 | 0.51279 | 3.20632 | 4.33081 |
| 96 | C34 | -0.76743 | 0.14142 | 3.12 | 4.38517 |
| 97 | C34 | -0.73847 | -0.36418 | 2.84303 | 4.31195 |
| 98 | C34 | -0.79479 | -0.72852 | 2.92404 | 4.1 |
| 99 | C34 | -0.99980 | -0.81364 | 3.44 | 3.50034 |
| 99* | C34 | -0.98357 | -0.09859 | 3.14 | 3.50590 |
| 100 | C34 | -0.90219 | 0.19958 | 2.82 | 3.54325 |
| 101 | C34 | -0.83212 | 0.73816 | 2.82894 | 3.78 |
| 102 | C35 | 0.69727 | 1.59295 | 3.12490 | 5.23849 |
| 103 | C35 | 0.65763 | -1.76053 | 2.89925 | 5.26 |
| 104 | C35 | 0.71138 | -1.40572 | 3 | 5 |
| 105 | C35 | 0.74046 | -1.37698 | 3.08124 | 4.88 |
| 106 | C35 | 0.80633 | -1.77810 | 3.16 | 4.59988 |
| 107 | C35 | 0.75706 | 1.85110 | 2.88 | 4.62305 |
| 108 | C35 | 0.72905 | 1.58305 | 2.84915 | 4.72 |
| 109 | C35 | 0.71138 | 1.40572 | 3 | 5 |

Sanity anchors in this table: row 2 (ã = ẽ = -1) is the point where every A_l passes (p. 123). Rows 6, 18 have ã = -4^(1/3) (N = 1/2). Rows 41, 49 have ã = 3^(2/3) (N = 1/3) and a(1 - e) = 1.

### 3.3 Transcription 2: Table VII.1 (G_{1/p} ∩ E_{1/p}, μ = 0)

Source: p. 252 (PDF 266). Read from the 150 dpi image, and rows 1-15 again at 300 dpi. Meaning: for N = 1/p (outer resonance, a = p^(2/3)), the asymmetric family G_{1/p} meets the symmetric family E_{1/p}^+ at two eccentricities e1, e2. c = √(a(1 - e^2)) is the area integral. The row p = ∞ is the parabolic limit (GL at θ = π, p. 262).

**Check I ran (transcription check).** c = √(p^(2/3)(1 - e^2)) for each pair (log `t71_check.log`). 39 of 40 pairs agree to ≤ 2e-9. Two do not, and the image shows them as printed: p = 6, c2 (misfit 6.4e-7; 0.322532932 would fit, so probably a transposed "932/293"); p = 12, c1 (misfit 6.0e-6; e1 gives c = 1.959287642). The p = ∞ row agrees with c = √(2d) for d1 = 3.0681355, d2 = 0.04351288 (p. 262) to 1e-8.

| p | e1 | e2 | c1 | c2 |
|---|---|---|---|---|
| 1 | (none) | 0.9175580167816 | (none) | 0.397601919 |
| 2 | 0.03652000274671 | 0.9592682766639 | 1.259080585 | 0.355923512 |
| 3 | 0.1221081097494 | 0.9717785352926 | 1.431456936 | 0.340219546 |
| 4 | 0.2011204732351 | 0.9779197645660 | 1.554964915 | 0.331736188 |
| 5 | 0.2671135721069 | 0.9816217379803 | 1.647844066 | 0.326326694 |
| 6 | 0.3218370252145 | 0.9841213933602 | 1.720440843 | 0.322532293 |
| 7 | 0.3676726526559 | 0.9859353649015 | 1.778940452 | 0.319702576 |
| 8 | 0.4065730081956 | 0.9873190391055 | 1.827236590 | 0.317497181 |
| 9 | 0.4400138013117 | 0.9884137591018 | 1.867896574 | 0.315722471 |
| 10 | 0.4690971294234 | 0.9893043817106 | 1.902681630 | 0.314258404 |
| 11 | 0.4946520818323 | 0.9900450844074 | 1.932841275 | 0.313026528 |
| 12 | 0.5173106310825 | 0.9906721638200 | 1.959293620 | 0.311973288 |
| 13 | 0.5375615689308 | 0.9912109211981 | 1.982703021 | 0.311060428 |
| 14 | 0.5557884042777 | 0.9916795592362 | 2.003608610 | 0.310259851 |
| 15 | 0.5722961272161 | 0.9920915904449 | 2.022409904 | 0.309549173 |
| 16 | 0.5873303250571 | 0.9924568397231 | 2.039427458 | 0.308918884 |
| 17 | 0.6010910213447 | 0.9927830786156 | 2.054918288 | 0.308358299 |
| 18 | 0.6137428356630 | 0.9930758721560 | 2.069091068 | 0.307870889 |
| 19 | 0.6254225451620 | 0.9933448143246 | 2.082117129 | 0.307342419 |
| 20 | 0.6362447862 | 0.9935859201040 | 2.094138573 | 0.306945540 |
| ∞ | 1 | 1 | 2.477149760 | 0.295001294 |

(p = 20, e1 is printed with only 10 decimals.)

### 3.4 Transcription 3: Appendix Table 3 (families TL, GL; parabolic limit, μ = 0 and O(μ))

Source: p. 335 (PDF 349). Read from the 150 dpi image. The text layer dropped the leading "1" of most rows in columns c and θ; the image is clean.

Meaning (pp. 256-262, 327): c = G = ε'√(2d) labels a parabola with pericentre distance d from P1. θ_cr is the angle θ of the collision family TL (parabolas that hit P2), blank where no collision is possible (d > 1, i.e. c > √2). θ is the pericentre angle on the generating family GL (the limit of G_{1/p} as p → ∞). ⟨Tr_1⟩ and ⟨ϖ̇_1⟩ are defined by (VII.56), p. 261 (image checked): ⟨Tr_1⟩ = 12π^2(p+q)^2 ⟨∂^2[R]/∂θ^2⟩ and ⟨ϖ̇_1⟩ = ⟨∂[R]/∂G⟩, where ⟨f⟩ is the coefficient of a^(-3/2) in f (p. 258). With the formulas on p. 260 this gives Tr_1 ≈ a^(5/2)⟨Tr_1⟩ and ϖ̇_1 ≈ a^(-3/2)⟨ϖ̇_1⟩ as a → ∞. So the two columns scale differently. All on ε' = +1. I did no numeric check of this table, except that GL must reach θ = π at c1 = 2.4771498 and c2 = 0.29500129 (p. 262). The table's ends (θ = 2.32197 at c = 2.4; θ = 3.01941 at c = 0.3) are consistent with that.

| c | θ_cr | θ | ⟨Tr_1⟩ | ⟨ϖ̇_1⟩ |
|---|---|---|---|---|
| 2.4 | | 2.32197 | -0.0880 | 0.02031 |
| 2.3 | | 2.00629 | -0.3782 | 0.01902 |
| 2.2 | | 1.82684 | -1.0919 | 0.01188 |
| 2.1 | | 1.70201 | -2.6495 | -0.00586 |
| 2.0 | | 1.60134 | -5.6900 | -0.03982 |
| 1.9 | | 1.51090 | -11.0040 | -0.09424 |
| 1.8 | | 1.42501 | -19.3123 | -0.16859 |
| 1.7 | | 1.34317 | -30.9403 | -0.25390 |
| 1.6 | | 1.26881 | -45.6212 | -0.33081 |
| 1.5 | | 1.20860 | -62.7387 | -0.37149 |
| 1.4 | 0.08646 | 1.17157 | -82.0389 | -0.34664 |
| 1.3 | 0.31006 | 1.16725 | -104.2218 | -0.23669 |
| 1.2 | 0.50655 | 1.20242 | -130.4601 | -0.04362 |
| 1.1 | 0.70446 | 1.27837 | -160.5341 | 0.20356 |
| 1.0 | 0.90413 | 1.39085 | -190.8398 | 0.45524 |
| 0.9 | 1.10380 | 1.53278 | -214.1978 | 0.65729 |
| 0.8 | 1.30155 | 1.69737 | -222.1931 | 0.76781 |
| 0.7 | 1.49567 | 1.87993 | -208.8556 | 0.76755 |
| 0.6 | 1.68474 | 2.07922 | -173.5282 | 0.66160 |
| 0.5 | 1.86766 | 2.29964 | -121.3409 | 0.47411 |
| 0.4 | 2.04358 | 2.56049 | -61.2264 | 0.23933 |
| 0.3 | 2.21194 | 3.01941 | -2.7368 | -0.00721 |

(⟨Tr_1⟩ < 0 on the whole of GL, so GL is stable to first order in μ, as the text says for G_N on p. 251.)

### 3.5 Control for Ir^+ (the "h" candidate) — closed form, no table

The book has no numeric table for Id or Ir. The control is (III.18)-(III.19), p. 107 (image checked):
- x1 + i x2 = a exp[i(n - 1)t], y1 + i y2 = i n a exp[i(n - 1)t], with n = -a^(-3/2) on Ir.
- Period T = 2π/|n - 1| = 2π/(1 + a^(-3/2)); trace Tr = 2 cos T.
- Ir^+ is the part with a < 1 (p. 289). These formulas hold at μ = 0 only, in the P1-centred frame.
- Appendix Table 2 rows with ε'e < 0 give the retrograde second-kind families E_N^± that end on Ir for even q (p. 289). They are not Ir itself.

**Nearest tabulated numbers: Appendix Table 2, block N = 2, rows with ε'e < 0 (p. 333, PDF 347).** Fig. VIII.15 (p. 289) shows E_2^- and E_2^+ attached to Ir^+ near ã ≈ ±2^(-2/3) ≈ ±0.630, in the lower (retrograde) half. These rows lie on E_2^±, **adjacent to, not on, Ir^+**. N = 2 means a = 0.629961. Read from the page image at 200 dpi. Columns: ε'e, c, θ_cr (θ on the collision family T_2; blank = no collision possible), Tr_1^+ and ϖ̇_1^+ (on E_2^+, θ = 0), Tr_1^- and ϖ̇_1^- (on E_2^-, θ = π/2). Spot check: c = ε'√(a(1 - e^2)) holds (e.g. e = 0.1 gives 0.7897).

| ε'e | c | θ_cr | Tr_1^+ | ϖ̇_1^+ | Tr_1^- | ϖ̇_1^- |
|---|---|---|---|---|---|---|
| -0.9 | -0.3459 | 0.5904 | -134.0022 | 0.4462 | -36.9488 | 0.3533 |
| -0.8 | -0.4762 | 0.6666 | -99.9431 | 0.2230 | -41.9896 | 0.9001 |
| -0.7 | -0.5668 | 0.8456 | -64.9806 | 0.0063 | -60.6526 | 2.2473 |
| -0.6 | -0.6349 | 1.3040 | -38.0283 | -0.1818 | -535.5500 | 25.5374 |
| -0.5 | -0.6873 | | -20.1067 | -0.3404 | 42.8576 | -4.6686 |
| -0.4 | -0.7274 | | -9.4007 | -0.4764 | 12.2821 | -2.4673 |
| -0.3 | -0.7571 | | -3.6574 | -0.5985 | 4.0087 | -1.7545 |
| -0.2 | -0.7776 | | -1.0162 | -0.7158 | 1.0409 | -1.3850 |
| -0.1 | -0.7897 | | -0.1217 | -0.8380 | 0.1221 | -1.1488 |

As e → 0 both Tr_1 values go to 0, as expected where E_2^± meet the circular family. For small retrograde e, Tr_1^+ < 0 (stable to first order) and Tr_1^- > 0 (unstable).

---

## 4. Citation mining

The bibliography has about 200 entries (pp. 347-356; my regex count is 207). Counts below are mentions in the body text (pp. 1-326), from the text layer.

**Most used, on topic (in order of weight):**
1. **Hénon, M. (1968).** "Sur les orbites interplanétaires qui rencontrent deux fois la Terre." Bull. Astron. Sér. 3, 3(3), 377-402. 15 mentions. The source of the S families A/B/C, Fig. III.17 and Tables 1-9 (τ/π, η/π, ε, ε', ε″, a, e, C; Table 1 is the hyperbolic part of A_0). **HELD** (`henon-1968-orbites-interplanetaires-...-french.pdf`; corrected by the corpus filer). Highest-value companion for #944 and #948.
2. **Bruno, A. D. (1972b).** Researches on the restricted three-body problem II. Periodic solutions and arcs for μ = 0. Preprint IPM No. 75 (Russian) = **Celest. Mech. 18 (1978) 9-50.** 7 mentions. = the held "Brjuno 1978 part II".
3. **Bruno (1973).** Researches ... III. Properties of solutions for μ = 0. Preprint IPM No. 25 = **Celest. Mech. 18 (1978) 51-101.** 4 mentions. = the held "Brjuno 1978 part III".
4. **Bruno (1978b).** Extrema of the Hamiltonian function on families of arc-solutions ... μ = 0. Preprint IPM No. 103 (Russian). 9 mentions. Holds the numeric tables behind Chapter V. NOT HELD.
5. **Bruno (1976).** Periodic solutions of the second kind ... Preprint IPM No. 95 (Russian). Holds the fuller versions of Appendix Tables 1-2 and the program listing. NOT HELD.
6. **Bruno (1980a, b, c).** Asymptotics of second-kind solutions (No. 51); trajectories of collisions (No. 148); trajectories with consecutive collisions (No. 149). All Russian preprints. NOT HELD.
7. **Bruno (1981).** Generating arc-solutions of the restricted three-body problem. Preprint IPM No. 25 (Russian). Fuller versions of Appendix Tables 5-8. See the ambiguity note below.
8. **Bruno (1978a).** On periodic flybys of the moon. Preprint IPM No. 91 = **Celest. Mech. 24 (1981) 255-268.** Directly relevant to #948. HELD (`bruno-1981-on-periodic-flybys-of-the-moon-...`).
9. **Hitzl, D. L. and Hénon, M. (1977a).** Critical generating orbits for second species periodic solutions of the restricted problem. Celest. Mech. 15, 421-452. HELD.
10. **Hitzl and Hénon (1977b).** The stability of second species periodic orbits in the restricted problem (μ = 0). Acta Astronaut. 4, 1019-1039. HELD.
11. **Hitzl (1977a)** (Perko timing condition, A&A 54, 47-55): NOT HELD. **Hitzl (1977b)** (AIAA J. 15, 1410-1418): HELD.
12. **Perko, L. M.** 1974 (SIAM J. Appl. Math. 27, 200-237), 1976 (Celest. Mech. 14, 395-427), 1977 (Celest. Mech. 16, 275-290), 1981 (Celest. Mech. 24, 155-171). 12 mentions. Bruno criticises Perko's treatment of arc multiplicities (p. 177). 1974, 1976, 1977 and 1981 all HELD.
13. **Krasinskii, G. A.** 1963 (closed collision trajectories; Byull. ITA 9, 154-163), 1973a, 1973b (Minor Planets). 14 mentions. The K-families are his. NOT HELD.
14. **Guillaume, P.** 1969 (A&A 3, 57-76), 1970 (thesis, Liège), 1974 (A&A 37, 209-218), 1975 (Celest. Mech. 11, 213-254, "Linear analysis of one type of second species solutions"). 11 mentions. 1975 HELD (`guillaume-1975a-...`); 1969, 1970, 1974 not held.
15. **Uno, T. (1937).** Jap. J. Astron. Geophys. 15, 149-191. First correct E_N existence proof. NOT HELD.
16. **Message, P. J.** 1958, 1959, 1966, 1970; **Message and Taylor 1978**; **Frangakis** 1968, 1973a, b. Asymmetric G_N history. NOT HELD.
17. **Hénon and Guyot (1970).** Stability of periodic orbits in the restricted problem. In: Giacaglia (ed.), Periodic Orbits, Stability and Resonances, Reidel, 349-374. HELD.
18. **Hénon 1965a, 1965b** (Ann. Astrophys. 28, 499-511; 992-1007). HELD. **Hénon 1969** (A&A 1, 223-238). HELD. Hénon 1970 (A&A 9, 24-36, Hill's case non-periodic). NOT HELD.
19. **Broucke, R. (1968).** Periodic orbits in the restricted three-body problem with Earth-Moon masses. NASA Tech. Rep. 32-1168, Pasadena (as the book cites it). HELD. (Also Broucke 1962 thesis, Louvain, NOT HELD.) The book calls him "M. R. Broucke".
20. **Szebehely, V. (1967).** Theory of Orbits. HELD. **Arenstorf (1963)** (Amer. J. Math. 85, 27-35). HELD. **Moulton** 1913 (Proc. London Math. Soc. 11, 367-384, ejection orbits; NOT HELD) and 1920 (Periodic Orbits; NOT HELD). Colombo, Franklin and Munford (1968) and Markellos (1974) are cited as μ ≈ 0.001 agreement checks (NOT HELD). Deprit and Henrard (1965, 1969, 1970) on L3 asymptotic orbits and the Trojan manifold (NOT HELD). Barrar (1965) (NOT HELD). Poincaré 1892-1899 and 1902a,b.

**Status of the held list against this book:**

| Held item | In this book? |
|---|---|
| Hénon 1997 (Generating Families) | not cited (published after the book) |
| Hénon 2001 | not cited (after the book) |
| Hitzl-Hénon 1977a, 1977b | cited |
| Breakwell-Perko 1974 | **not cited** |
| Brjuno 1978 parts II and III | cited, as Bruno 1972b and 1973 (Russian preprints of Celest. Mech. 18, 9-50 and 51-101) |
| Bruno 1981 | cited, but see the note below |
| Broucke 1968 | cited |
| Szebehely 1967 | cited |
| Hénon 1965a, b | cited |
| Hénon 1969 | cited |
| Hénon 1973 | **not cited** |
| Hénon-Guyot 1970 | cited |
| Bruno-Varin 2006 | not cited (after the book) |
| Perko 1974, 1976, 1977, 1981 | cited (all held) |
| Arenstorf 1963 | cited |
| Strömgren 1933 | **not cited** |

**Ambiguity (resolved by the corpus filer: the held "Bruno 1981" is the Celest. Mech. 24:255 paper, i.e. the book's 1978a).** The held "Bruno 1981" could be either of two items. In this book, "Bruno 1981" is the Russian preprint *Generating arc-solutions of the restricted three-body problem* (IPM No. 25, 1981). But the English paper in Celest. Mech. 24 (1981) 255-268 is the book's **Bruno 1978a**, *On periodic flybys of the moon*. The corpus holds the Celest. Mech. paper (`bruno-1981-on-periodic-flybys-of-the-moon-...`). It is the one relevant to #948. The IPM No. 25 preprint is not held.

**Wanted, not held (priority order for these routes):** Krasinskii 1963; Bruno 1978b and 1981 preprints (Russian; for the full numeric tables); Hitzl 1977a; Uno 1937; Moulton 1913. These are now one Tier C row of the wanted list (batch 29).

---

*Table IV.1 and Table VII.1 transcriptions are filed beside the PDF as `cyclers_pdf/papers/bruno-1994-restricted-3-body-problem-plane-periodic-orbits-de-gruyter-expositions-math-17-doi-10.1515-9783110901733-table-IV-1-transcription.txt` and `...-table-VII-1-transcription.txt`. The check scripts were scratch and are not kept.*
