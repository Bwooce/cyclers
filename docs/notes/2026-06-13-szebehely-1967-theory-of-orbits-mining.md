# Szebehely 1967, *Theory of Orbits* — periodic-orbit / golden mining (#185)

**Date:** 2026-06-13
**Source:** V. Szebehely, *Theory of Orbits: The Restricted Problem of Three Bodies*,
Academic Press, 1967 (cited as "Szebehely 1967, *Theory of Orbits*" + page).
**Task:** mine the foundational CR3BP textbook for periodic-orbit GOLDENS / SEEDS
usable in our CR3BP lane (#182 corrector, #218 continuation, #116 JPL-oracle
cross-check). Reproduce-before-use; do NOT wire here.

**Verdict: ADOPT-CANDIDATE (libration-point goldens, Appendices I-IV of Ch. 4) +
REFERENCE-ONLY for the periodic-orbit families (Ch. 8-9 are graphical/qualitative,
no tabulated reproducible state vectors).** One convention trap must be cleared
before any number is wired (see "Convention trap" below). A follow-on task is
proposed.

---

## TOC orientation

The PDF is the complete book (~668 book pages; the scan renders ~2 book pages per
tool page, so PDF page N ≈ book page 2·(N−9)+2; back-matter indexes confirm full
coverage). Chapters relevant to us:

- **Ch. 1** Description of the restricted problem — standard CR3BP derivation
  (sidereal → synodic → dimensionless). Confirms the equation set we use.
- **Ch. 4** Totality of Solutions — **Appendices I-IV (book pp. 214-230)** carry the
  high-precision numeric tables (collinear-point locations, Jacobi constant at the
  points, Hessian Ω_xx/Ω_yy, all vs. μ). THE wireable numerics.
- **Ch. 5** Motion near the equilibrium points — Appendices I-III (pp. 309-318):
  roots of the characteristic equation at collinear / triangular points (eigenvalue
  goldens — not transcribed this pass; candidate for a future stability cross-check).
- **Ch. 8** Periodic Orbits — definitions, analytic continuation, Whittaker's
  criterion. Existence theory; **no tabulated ICs**.
- **Ch. 9** Numerical Explorations — Copenhagen category (p. 455), periodic lunar
  orbits (p. 497), motion around triangular points (p. 514). The periodic-orbit
  **families are presented as figures + C(x) characteristic curves** reproduced from
  Strömgren's Copenhagen reports; numbers quoted in text are read-from-figure
  approximations, NOT reproducible 15-digit state vectors.

## Chapters / pages actually read

Front matter + TOC (book pp. v-xvi); Ch. 1.1-1.5 (pp. 7-17, equation set);
Ch. 4 references + **Appendices I-IV in full (pp. 212-230)**; Ch. 9.2-9.4
(pp. 446-461, organization of numerical results, system-comparison tables,
Copenhagen Class (a)/(c)); back-matter indexes (pp. 654-668).

---

## ADOPT-CANDIDATE goldens — collinear-point tables (Ch. 4 Appendices I-IV)

All transcribed exactly from the printed tables. **Candidate / reproduce-before-use.**
Szebehely numbers the collinear points 1st/2nd/3rd; his "C_i" is the **standard**
Jacobian constant (includes the +μ(1−μ) term — see trap). Ω_xx/Ω_yy are the
second derivatives of his standard Ω at the point.

### Appendix I.D — First collinear point, Earth-Moon range (book p. 216)
Columns: μ, x₁, C₁, (Ω_xx)₁, (Ω_yy)₁
```
μ=0.0120   x1=-1.15510 01298   C1=3.18282 40063   Ωxx=7.38673 72270   Ωyy=-2.19336 88635
μ=0.0121   x1=-1.15548 72863   C1=3.18371 40528   Ωxx=7.38282 21914   Ωyy=-2.19141 10957
μ=0.0123   x1=-1.15625 37037   C1=3.18548 46255   Ωxx=7.37506 09860   Ωyy=-2.18753 04930
```

### Appendix I.E — First collinear point, solar-system bodies (book p. 217)
The critical (terminating) μ value and the Earth-Moon entry:
```
EARTH-MOON   μ=0.01215 06683
CRITICAL VALUE  μ=0.03852 08965  ->  x1=-1.21443 88479  C1=3.55119 28604
DARWIN'S VALUE  μ=0.09090 90909  ->  x1=-1.26408 29085  C1=3.59418 16095
```

### Appendix II.D — Second collinear point, Earth-Moon range (book p. 220)
Columns: μ, x₂, C₂, (Ω_xx)₂, (Ω_yy)₂
```
μ=0.0120   x2=-0.83765 86648   C2=3.19880 45659   Ωxx=11.28425 11404   Ωyy=-4.24212 55702
μ=0.0121   x2=-0.83726 43231   C2=3.19982 77931   Ωxx=11.29152 37750   Ωyy=-4.14576 18875
```
(Note II.A end at μ→0.5 gives the clean check x₂=0.0, C₂=4.25, Ω_xx=17.0,
Ω_yy=−7.0 — an exact algebraic anchor for the μ=1/2 case.)

### Appendix III.D — Third collinear point, Earth-Moon range (book p. 224)
Columns: μ, x₃, C₃, (Ω_xx)₃, (Ω_yy)₃
```
μ=0.0120   x3=+1.00499 99054   C3=3.02385 26541   Ωxx=3.02111 60924   Ωyy=-0.01055 80462
μ=0.0121   x3=+1.00504 15697   C3=3.02405 03851   Ωxx=3.02129 30404   Ωyy=-0.01064 65202
```
(III.A at μ→0 reduces cleanly: x₃=1.0, C₃=3.0, Ω_xx=3.0, Ω_yy=0.0 — anchor.)

### Appendix IV — Jacobian constant on the x-axis, C = 2Ω(x,0) (book pp. 226-229)
Dense C(x) tables for μ ∈ {1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.1, 0.2, 0.3}, x from
−3.0 to +3.0. These are zero-velocity-curve / Hill-region anchors: at v=0 the Jacobi
constant on the x-axis equals 2Ω(x,0). Useful as a smoke test that our Ω/zero-
velocity machinery reproduces the printed C(x) curve for a given μ (sample, e.g.
μ=0.1, x=−1.0 → C=4.97029 60396).

### Standard μ=1/2 reference values (Ch. 9 Table III, book p. 457)
Exact comparison table for the Copenhagen mass (μ=1/2), columns x, C̄, C:
```
L1  x=-1.1984   C̄=3.4568   C=3.7068
L2  x= 0        C̄=4.0000   C=4.2500
L3  x=+1.1984   C̄=3.4568   C=3.7068
L4  x=0,y=√3/2  C̄=2.7500   C=3.0000
L5  x=0,y=-√3/2 C̄=2.7500   C=3.0000
```
This row is the cleanest cross-check pair: it prints **both** C̄ (our convention)
and C (Szebehely's standard) side by side, and gives the exact identity
**C̄ = C − μ(1−μ)** (here 3.7068 − 0.25 = 3.4568; 3.0000 − 0.25 = 2.7500).
Standard-convention algebraic identity worth pinning: **C(L4,5) = 3** exactly,
independent of μ (p. 451).

---

## Convention trap (must clear before wiring ANY number above)

1. **Frame mirror.** Szebehely's *standard reference system* (Fig. 9.1a, p. 449)
   puts the **larger** mass 1−μ at **P1(+μ,0)** and the **smaller** mass μ at
   **P2(μ−1,0)**. Our code (`core/cr3bp.py`) puts the **larger** mass 1−μ at
   **(−μ,0)** and the smaller μ at **(1−μ,0)** — mirror-flipped in x. The Jacobi
   constant is mirror-invariant (depends on x², r1, r2), so **C values transfer
   directly**; collinear-point **x-coordinates are sign/label-permuted** between the
   frames and must be remapped before any position comparison.

2. **Two Jacobi constants.** Szebehely uses Ω̄ = ½r² + (1−μ)/r1 + μ/r2 (→ C̄, =
   Wintner's/Birkhoff's constant) and Ω = Ω̄ + ½μ(1−μ) (→ C, his "standard").
   **Our `jacobi_constant()` = C̄** (no μ(1−μ) term). The **Appendix I-IV tables
   print C_i = the standard C**, so to compare to our code use **C̄ = C − μ(1−μ)**,
   or recompute. Table III (p. 457) is safe because it prints both.

3. **Hessian.** Ω_xx/Ω_yy in the appendices are derivatives of his standard Ω; the
   μ(1−μ) term is constant so the second derivatives are convention-independent —
   these can feed an eigenvalue/stability cross-check directly once the point
   location is remapped.

---

## REFERENCE-ONLY (foundational, not wireable)

- **Periodic-orbit families (Ch. 9.4 Copenhagen Classes (a)-(r), lunar orbits,
  triangular-point families).** Presented as figures + C(x) characteristic curves
  from Strömgren's reports; the in-text numbers (e.g. Class (a) C ≈ 3.71, 3.55,
  3.00, 2.50; collision orbit at C ≈ 2.10; minimum C ≈ 1.31; axis crossings
  x₁ ≈ 0.76, 1.37, …) are **read-from-figure approximations**, not reproducible
  state vectors. Good for qualitative family-topology orientation (which class is
  retrograde/direct, around which point, where collision orbits sit), **not** for a
  golden test. A few exact algebraic values do appear: μ=1/2 consecutive-collision
  orbit around L2 at **C=2.432913** (p. 461), and the L2-limit ratio 0.2278093…
  with start C₂=4.25.
- **Zero-velocity curves / Hill regions** (Ch. 4.6-4.7) — qualitative; the numeric
  substance is the Appendix IV C(x) tables already captured above.
- **Equations of motion** (Ch. 1) — confirms our standard CR3BP set; no new numerics.
- **Stability eigenvalue tables** (Ch. 5 Appendices I-III, pp. 309-318) — roots of
  the characteristic equation at the collinear/triangular points vs. μ. Not
  transcribed this pass; flagged as a future candidate (see follow-on).

---

## Discipline / sourcing

- Every transcribed number is a **published value** (Szebehely's tables), so a golden
  built on them satisfies "EXPECTED side traces to a source" — provided the two
  convention conversions above are applied. They are NOT values our code computed.
- Nothing wired here. The L-point appendix numbers are candidates → reproduce first.

## Proposed follow-on task

**"Szebehely L-point golden cross-check (same-model)"** — small, costless:
1. For the Earth-Moon μ (use the published μ=0.01215 06683, App I.E) and for μ=1/2,
   compute our collinear-point locations + C̄ and compare to App I-III (after the
   frame remap) and to Table III. Pass criterion: agreement to the printed ~13 sig
   figs in C̄ after the C̄ = C − μ(1−μ) conversion.
2. Optionally add the App IV C(x) zero-velocity smoke test for one μ.
3. If the Ch. 5 characteristic-root tables are wanted, a second pass can transcribe
   pp. 309-318 for an eigenvalue/stability golden.
This gives the #116 oracle lane a sourced, independent (non-JPL, 1967-vintage)
anchor for the libration-point geometry and Jacobi calibration, decoupled from any
value our own pipeline produced.

## #116 wishlist

Szebehely 1967 *Theory of Orbits* mined — **strike off the #116 wishlist.**
Outcome: ADOPT-CANDIDATE (L-point appendix goldens, pending the convention follow-on)
+ REFERENCE-ONLY for periodic-orbit families.
