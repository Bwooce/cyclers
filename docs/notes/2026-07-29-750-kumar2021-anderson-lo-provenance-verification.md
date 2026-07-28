# `#750`: digit-for-digit provenance check — Kumar et al. 2021 (AAS 21-651) vs. Anderson & Lo 2011's confirmed 3:4-LO

**Task:** `#750`. Checks whether the catalogued N=5 CRNBP torus discovery
(`europa-3-4-crnbp-torus-jupiter-2026`, `data/catalogue.yaml`, task chain
`#714`→`#736`) — built by continuing "the planar Jupiter-Europa 3:4 resonant
orbit of Kumar et al. 2021" — is, digit-for-digit, the *same specific orbit*
as this session's OWN independently-built and reviewer-confirmed 3:4-LO
orbit (`#753`→`#755`, reproducing Anderson & Lo 2011's Table 1 eigenvalue).
Read the Kumar et al. 2021 PDF directly (already in the private corpus at
`cyclers_pdf/papers/kumar-anderson-delallave-gunter-2021-europa-ganymede-
resonant-orbits-ccr4bp-AAS-21-651-arxiv-2109.14815.pdf`), the catalogue row's
own `crnbp_provenance` block, `src/cyclerfinder/search/jovian_resonant_
families.py`'s stored `3:4-LO` candidate, and this project's own #690/#724
seed-construction script (`scripts/verify_724_rerun_continuation.py`). No
code or catalogue changes — pure provenance verification.

---

## Verdict (read this first)

**NOT the same specific orbit.** Three distinct numeric objects exist under
the "Jupiter-Europa 3:4 exterior/outer resonant orbit" label, and they are
quantitatively very different — most tellingly, their closest approach to
Europa spans more than two orders of magnitude (409,500 km / 18,700–22,000 km
/ 1,657 km). Whether they are nonetheless *members of the same underlying
continuous family* (same p:q=3:4 resonance ratio, same "exterior"/"Outer"
branch, same base dynamical system — a real, non-trivial claim, just a
weaker one than "the same orbit") is **genuinely indeterminate** without a
dedicated continuation study this task did not attempt (and which this
session's own `#753`/`#755`/`#756`/`#758` chain has already shown is
fragile and non-trivial even for closely-related in-project comparisons, due
to well-documented "fractal sensitivity" in this exact parameter region).

**Practical consequence for `#761`:** this project's own catalogued torus's
seed orbit is a *third*, even more distant point than either paper's, and is
not verified to lie on Anderson & Lo's specific unstable/close-flyby branch
at all. `#754`/`#759`'s homoclinic/heteroclinic connection machinery — built
entirely around the close-Europa-flyby instability mechanism of the
confirmed 3:4-LO orbit — does **not** straightforwardly carry over to the
catalogued torus without first re-deriving/continuing the catalogued seed
onto (or near) that specific branch, which is itself an open, unverified
question, not a citation-chasing exercise.

---

## 1. The three numeric objects, side by side

| Quantity | Kumar et al. 2021's own seed (p.8, "3:4 Jupiter-Europa Resonant CCR4BP Tori") | Anderson & Lo 2011's confirmed 3:4-LO (this session's `#753`→`#755` chain) | This project's own catalogued-torus seed (`#690`/`#724`) |
|---|---|---|---|
| Jupiter-Europa `μ` | **2.5266448850435×10⁻⁵** (derived from the paper's own Table 1 `Gm_i` values — see §2 below; bit-identical to Anderson & Lo's own stated value) | 2.5266448850435×10⁻⁵ (`ANDERSON_LO_MU`, p.169 of the 2011 paper, verbatim) | 2.528017724591319×10⁻⁵ (this project's own CRNBP-reduction registry value; 0.054% relative from the papers' shared value) |
| Jacobi constant `C` | **3.0041** ("Starting from a periodic orbit with Jacobi constant value 3.0041 ... (arbitrarily chosen)", p.8) | 2.99163956830415 (`ANDERSON_LO_C_FLYBY`, p.179 — the energy of the paper's *own separate* ballistic 3:4↔5:6 flyby trajectory, NOT arbitrary) | 2.9040363618665213 (this project's own re-derivation this task, §4 below; built by fixing period = 8π exactly, not by matching either paper's `C`) |
| `x0` (perpendicular crossing) | **not published** (no IC table anywhere in the paper — confirmed by direct read, §3) | −1.4304078294961569 (opposite side of Jupiter from Europa — Anderson-Lo's "Outer" region) | +1.6101945301604663 (same side as Europa, near apoapsis of a modest ellipse) |
| Period / 2π | not published directly; implied ≈4 only in the μ₂→0 two-body limit (paper's own text, §4) | 4.085910 (2.1% off the "clean" `q=4` value — a real, converged offset, reviewer-accepted as expected for a strongly unstable family) | exactly 4.0 (period = 25.132741228718345 = 8π to machine precision, by construction) |
| Closest approach to Europa | **22,052 km → 18,721 km** as Ganymede's mass ratio μ₃ is turned on (paper's own stated numbers, p.8) | **1,657 km** (≈97 km surface altitude — the paper's own attributed instability mechanism) | **409,500 km** (this task's own numerical check, §4) |
| Max/unstable eigenvalue `\|λ\|` | not stated for the base periodic orbit (Kumar 2021 reports Floquet multipliers of the *torus*, not a Table-1-style periodic-orbit eigenvalue) | 1036.116116695996 (recovered vs. Anderson & Lo's own Table 1 target 1036.116088, rel. err 2.8×10⁻⁸) | not computed (this seed was never independently classified for stability in this project's own work — it is used purely as a continuation starting point) |

**Headline discriminator:** closest approach to Europa differs by roughly an
order of magnitude between Kumar's own reported orbit (~19,000–22,000 km) and
Anderson & Lo's confirmed 3:4-LO (~1,657 km), and by another order of
magnitude again between Kumar's orbit and this project's own catalogued-torus
seed (~409,500 km). All three carry the "3:4" resonance-ratio label and the
same "exterior"/"Outer" broad geometric classification, but none of the three
is a close numeric match to either of the other two.

---

## 2. Reading Kumar et al. 2021 directly — what the paper does and does not state

Read via `pdftotext -layout` plus direct page-image inspection (pages 7-8,
Table 1) to rule out an OCR/column-alignment misread before trusting the
numbers.

**Abstract (p.1):** "we compute tori corresponding to **exterior**
Jupiter-Europa and **interior** Jupiter-Ganymede PCRTBP resonant periodic
orbits" — confirms the catalogue row's own "EXTERIOR, never interior" note
(`#724` correction) is reading the right word for the right body pair.
"Exterior" here is a semi-major-axis/two-body-limit notion (the resonant
orbit's own period exceeds Europa's, hence `a > 1` in Europa units for
`p:q=3:4` since `q(4) > p(3)`) — this matches Anderson & Lo's own "Outer
region" (the "O" in "3:4-LO") in the broad sense that both describe an orbit
lying farther from Jupiter than the secondary moon on average, **not** two
unrelated concepts (an initial concern this task was specifically asked to
rule out). So the "exterior"/"Outer" correspondence checks out qualitatively.

**Section "3:4 Jupiter-Europa Resonant CCR4BP Tori" (p.8), the actual source
of the seed used downstream in the paper:**

> "Starting from a periodic orbit with Jacobi constant value 3.0041 and
> ω = 3.097849 (arbitrarily chosen), the Jupiter-Europa frame planar CCR4BP
> stroboscopic map invariant circles computed at various μ₃ continuation
> steps are displayed in Figure 1... In physical units, the closest approach
> to Europa decreases from 22052 km to 18721 km [as μ₃ is turned on to its
> physical value]."

This is the **entire** numeric specification Kumar 2021 gives for its base
Jupiter-Europa 3:4 periodic orbit: a Jacobi constant (3.0041) and a torus
rotation-number-like parameter ω (3.097849), both explicitly flagged by the
authors themselves as **"arbitrarily chosen"** — i.e. Kumar picked some
representative member of the 3:4-exterior family to demonstrate their
continuation *method*, not to reproduce or match any specific prior orbit.
**No `(x0, ẏ0)` initial condition, no explicit period, and no Table-1-style
eigenvalue for the base periodic orbit appear anywhere in the paper's 20
pages** — confirmed by a full read (not just the seed section): no
subsequent table, figure caption, or appendix restates this in IC form. This
independently confirms (rather than just inherits) `variational_ccr4bp_torus
.py`'s own pre-existing docstring claim that "the published paper does NOT
list initial conditions / energy / rotation number to pixel-reproduce."

**Citation check (Introduction, p.1-2, and full reference list, p.19-20):**
Kumar 2021 cites Anderson & Lo 2010 [1] and 2011 [2] only as **general
background** for "the final approach to Europa after a series of flybys of
the Galilean moons" — a broad motivational citation, not a claim of reusing
either paper's specific computed orbit or numeric result. No equation, table,
or figure in Kumar 2021 references Anderson & Lo's Table 1, their `C_flyby`
value, or their `3:4-LO` naming. `#745`'s digest already found this citation
relationship; this task confirms it is exactly as loose as that digest
described (motivational, not a data-reuse citation) — the two-hop lineage
`#745` flagged as "plausible, not confirmed" remains exactly that: plausible
common ancestry (same resonance class, same authorship overlap via Anderson,
same base system), not a demonstrated shared numeric result.

**Mass ratio cross-check (bonus finding, resolves the task's "check whether
the same ~0.02–0.03% gap applies here" question):** Table 1 (p.8) lists
`Gm_i` values for Jupiter/Europa/Ganymede whose row labels, when checked
against this project's own registry (`core/satellites.py`, DE440-derived:
Europa `μ=3202.739 km³/s²`, Ganymede `μ=9887.834 km³/s²`), are printed
**transposed** — the value in the row labeled "Europa" (9886.997 km³/s²)
matches Ganymede's real mass, and the value in the row labeled "Ganymede"
(3200.9998 km³/s²) matches Europa's real mass (the Orbital Period column, by
contrast, is correctly assigned to each body). Recomputing `μ =
Gm_Europa/(Gm_Jupiter+Gm_Europa)` using the mass-consistent (not the
row-label-consistent) value gives **μ = 2.526644885043503×10⁻⁵ — bit-identical
(relative difference 1.1×10⁻¹⁵) to Anderson & Lo 2011's own stated
μ=2.5266448850435×10⁻⁵** (self-consistently confirmed via Kumar's own stated
final μ₃=7.804102777055038×10⁻⁵, which only reproduces exactly under this
reading of the table). **So Kumar 2021 uses the literal, identical-vintage
Jupiter-Europa mass ratio as Anderson & Lo 2011** (unsurprising given Anderson
is a common co-author) — there is **no** mass-ratio-vintage gap between the
two papers at all; the ~0.02–0.05% gap this session has documented elsewhere
(`#745`'s digest) is entirely between this *project's* DE440-derived value
and the *papers'* (shared, older-vintage) value, not between the two papers.
Flagged respectfully as a minor, apparently genuine transcription erratum in
Kumar et al. 2021's own Table 1 (GM column rows swapped relative to body
labels) — not itself load-bearing for this task's verdict, and not something
this task attempts to correct or report externally.

---

## 3. Anderson & Lo 2011's confirmed 3:4-LO — exact numbers (from this session's own chain)

Pulled directly from `src/cyclerfinder/search/jovian_resonant_families.py`'s
`_TABLE1_CANDIDATE_SEEDS["3:4-LO"]` and the `#755` results note (not just the
rounded prose in that note — the note's own reported full-precision numbers):

```
mu (Anderson & Lo 2011, p.169, verbatim) = 2.5266448850435e-5
Jacobi constant (ANDERSON_LO_C_FLYBY, p.179, verbatim) = 2.99163956830415
x0 = -1.4304078294961569, ydot0_sign = +1, half_crossings in {5,6,7}
period = 25.672528919046933  (period/2pi = 4.085910, a real 2.1% offset from
  the "clean" q=4 two-body value -- reviewer-accepted as expected for a
  strongly unstable family whose period need not track 2*pi*q, per the
  paper's own Eq. 6 text, p.171)
crossing_residual = 2.07e-13
max eigenvalue (Barden, real) = 1036.116116695996
  vs. Anderson & Lo's own Table 1 target 1036.116088 -- rel. err 2.77e-8
Europa closest approach = 0.00247 nondim = 1657 km (~97 km surface altitude)
Spatial envelope: x in [-1.430, 1.258], y in [-1.384, 1.384]
  (visually matches the paper's own Fig. 16(a) "flower" orbit)
```

This orbit was reviewer-CONFIRMED (coordinating session, `#755`'s own note)
as the paper's genuine 3:4-LO family member — on the weight of the
near-machine-precision eigenvalue match plus independent shape/close-approach
corroboration, despite the period offset. It is `-1.43` on the Jupiter side
*opposite* Europa (Anderson-Lo's own "Outer region" of their Loop/no-Loop ×
Inner/Outer taxonomy), makes an extreme close flyby of Europa, and is very
strongly unstable (`|λ|≈1036`).

---

## 4. This project's own catalogued-torus seed — exact numbers (re-derived this task)

The catalogue row's `crnbp_provenance.torus.seed_lineage` field states: *"Kumar-2021
JE 3:4 exterior resonant PO (CR3BP, perp residual 7.2e-13, re-verified this
session)"*. That seed is built by `scripts/verify_724_rerun_continuation.py`'s
`resonant_symmetric_orbit(mu, 3, 4)` — **a fundamentally different corrector
strategy from `jovian_resonant_families.py`'s `correct_symmetric_fixed_jacobi`**:
it *fixes the period at exactly `2*pi*q = 8*pi` (the pure two-body Kepler
value)* and Newton-solves only for `(x0, ẏ0)` to make that a genuine
perpendicular-crossing periodic orbit at that exact period, rather than
fixing the Jacobi constant and letting the period float (as the `#753`→`#755`
Anderson-Lo-reproduction chain does). Re-ran it this task (read-only, no
files changed) to get exact digits and a genuine closest-Europa-approach
number, which the project's own prior work had never computed for this seed:

```
mu (project CRNBP registry, ccr4bp.jupiter_europa_ganymede_default()) = 2.528017724591319e-05
x0 = 1.6101945301604663
vy0 = -0.9647838966320111
period = 25.132741228718345  (= 8*pi to machine precision, by construction)
perpendicular-crossing residual = 7.208453296892155e-13
  (matches the catalogue row's own quoted "7.2e-13" exactly)
Jacobi constant = 2.9040363618665213
r_from_Jupiter range = [0.8125, 1.6103]  (a=1.2114, e=0.329 -- matches the
  catalogue's own "a ~= 1.2114 Europa SMA, BETWEEN the moons' orbits" note)
Closest approach to Europa (this task's own numerical check,
  200,000-point dense propagation over one full period) = 0.6102 nondim
  = 409,518 km
```

This orbit is a **third, distinct point**: same side as Europa (unlike
Anderson-Lo's 3:4-LO, which sits on the *opposite* side of Jupiter), a
gentle, near-two-body-limit ellipse (e=0.33, never gets closer than ~410,000
km to Europa — over 20x farther than even Kumar's own reported orbit, and
~250x farther than Anderson-Lo's 3:4-LO), at yet a third distinct Jacobi
constant (2.9040, vs. Kumar's 3.0041 and Anderson-Lo's 2.99164). It was
explicitly built as a theoretical two-body-limit stand-in for "the object
CLASS Kumar 2021 computed" (this project's own `variational_ccr4bp_torus.py`
docstring already says as much, honestly, in its "honest-partial" framing) —
**not** a claim of matching Kumar's specific `C=3.0041` member, and this task
confirms it indeed does not match either paper's orbit numerically.

**A brief, inconclusive continuation attempt (this task, not adopted as
evidence, reported only in the interest of full disclosure):** tried a
small-step natural-parameter (Jacobi-constant) continuation from this seed
upward toward both Kumar's `C=3.0041` and Anderson-Lo's `C_flyby=2.99164`.
The continuation lost convergence around `C≈2.926` (a jump to a
non-perpendicular-crossing state), and separately, re-converging the exact
same seed at Anderson-Lo's own `μ` (a mere 0.054% change from the project's
registry `μ`) caused an immediate branch jump (`x0` from 1.61 to 0.97) rather
than a smooth perturbation — both symptomatic of the same "fractal
sensitivity" this session's own `#753` module docstring already documents at
length for this exact family/energy region. This is **not** strong evidence
either way about whether the three orbits share a connected family; it is
evidence that a quick, naive continuation attempt cannot resolve the
question, consistent with (not contradicting) the "indeterminate" verdict
below.

---

## 5. Convention/unit check (task requirement — ruling out a false mismatch)

All three numbers above are already in the **same convention**: planar CR3BP
rotating-synodic frame, non-dimensional units (`l_km=671,100` Europa SMA,
`t_s=48,843.878...` s), Jupiter-Europa primary/secondary pair, `x`-axis
through the primaries, secondary at `x=1-μ`. No rotation/reflection/unit
conversion is needed to compare `x0`, period, or Jacobi constant directly —
this task confirms the comparison above is a true apples-to-apples read, not
an artifact of a hidden convention difference. (The one real convention
subtlety — Kumar 2021's `ω`/rotation-number parameter is a property of the
*quasi-periodic torus/invariant circle*, not of the base periodic orbit
itself, since a periodic orbit alone has no internal rotation number — was
checked directly against the paper's Eq. 10-11 KAM-parameterization
formalism and is **not** a like-for-like quantity with Anderson & Lo's
`|λ|` Floquet-multiplier eigenvalue; this task does not attempt to convert
between them, since doing so is out of scope and not needed for the
verdict.)

---

## 6. Why "indeterminate" (not a clean "different family") is the honest call on the *family* question

Three considerations keep the broader "same family" question open rather
than closed negative:

1. **Same base system, confirmed bit-identical mass ratio** (§2) — both
   papers are working in the literal same dynamical system, not a
   mass-ratio-vintage-confounded comparison.
2. **Same resonance ratio and same broad geometric branch** — both are
   `p:q=3:4` "exterior"/"Outer" orbits, and neither this project's nor
   Anderson & Lo's own work has found a competing `3:4-LI`/`3:4-NI` branch to
   confuse this with (Table 1 lists only one 3:4 row).
3. **A genuinely plausible physical story exists** for why a single
   continuous family, sampled at nearby-but-different Jacobi constants,
   could show the closest-approach and eigenvalue differences observed: this
   exact energy region is independently documented (this session's own
   `#753` module) as "fractally sensitive" — small Jacobi-constant changes
   can produce wildly different eigenvalues/geometry for nearby family
   members. A ~0.4% Jacobi-constant difference (Kumar's 3.0041 vs.
   Anderson-Lo's 2.99164) is entirely consistent with landing on a
   qualitatively different-looking member of the *same* continuum, exactly
   as the `#755` reviewer ruling already accepted for 3:4-LO's own 2.1%
   period offset from the naive two-body value.

But no task (this one, or `#753`/`#755`/`#756`/`#758` before it) has actually
run a continuous, bifurcation-checked continuation connecting Kumar's
`C=3.0041` point to Anderson-Lo's `C=2.99164` point, so "same family, just
sampled at different energies" remains a **plausible, unconfirmed
hypothesis**, not a demonstrated fact — and this task's own brief attempt
(§4) hit exactly the kind of branch-jump behavior that would make such a
continuation nontrivial to trust even if attempted. Calling this "different
family, full stop" would overclaim a negative this task cannot support;
calling it "same family" would overclaim a positive nobody has verified.
**Indeterminate, with the specific-orbit question (not just the family
question) answered NO with high confidence** is the accurate summary.

---

## Summary answer (for the coordinating session, re: `#761`)

- **Same specific orbit? No** — high confidence. Kumar 2021's own stated
  Jacobi constant (3.0041) and closest-Europa-approach (18,700–22,000 km)
  are quantitatively distinct from Anderson & Lo 2011's confirmed 3:4-LO
  (2.99164, ~1,657 km) by roughly an order of magnitude on the approach
  distance, and Kumar 2021 explicitly calls its own choice "arbitrarily
  chosen" rather than a reproduction of any specific prior orbit. No IC
  table exists in Kumar 2021 to check digit-for-digit in the first place.
- **Same underlying continuous family (same resonance branch, different
  energy sample)? Indeterminate** — plausible on general grounds (identical
  base system/mass-ratio, same resonance ratio and geometric branch, and a
  known fractal-sensitivity mechanism that could explain the large
  quantitative gap at a small `C` difference) but not demonstrated by any
  continuation this task or any prior task has run.
- **Does this project's own catalogued torus's seed match either paper's
  orbit? No** — it is a third, even more distant point (closest Europa
  approach ~409,500 km), built as a theoretical two-body-limit proxy for the
  object *class*, not a reproduction of Kumar 2021's specific `C=3.0041`
  member.
- **Bonus finding:** Kumar et al. 2021 uses the bit-identical Jupiter-Europa
  mass ratio as Anderson & Lo (2011) (μ=2.5266448850435×10⁻⁵) — there is no
  mass-ratio-vintage gap between the two *papers*; the previously-documented
  ~0.02–0.05% gap is only between this *project's* DE440-derived value and
  the (shared) older-vintage value the papers both use. Also flagged: Table
  1 of Kumar et al. 2021 appears to have its Europa/Ganymede `Gm_i` row
  labels transposed relative to the bodies' real masses (the Orbital Period
  column is correctly assigned) — a minor, likely genuine transcription
  erratum, not itself load-bearing for this task's verdict.

**Practical implication (context only, per the dispatching session's own
framing — not a decision made here):** `#754`/`#759`'s 3:4-LO-specific
connection-building machinery targets an extreme, close-Europa-flyby,
strongly-unstable member of this resonance class. The catalogued N=5 torus's
own underlying seed is a gentle, distant, weakly-corrected member built for
a different purpose. Porting the connection machinery over is not a
"same-object, different-citation" shortcut — it would require first
establishing (via a dedicated, bifurcation-aware continuation study, not
attempted here) whether and how the catalogued seed's family connects to the
confirmed 3:4-LO branch at all.

---

## Verification

Read-only task: no source, test, or catalogue files changed. Re-derivations
in this note (§2 mass-ratio arithmetic, §4 seed re-computation and
closest-approach check) used only existing project code
(`cyclerfinder.core.cr3bp`, `cyclerfinder.core.ccr4bp`,
`cyclerfinder.search.cr3bp_periodic`, `scripts/verify_724_rerun_continuation
.py`'s own `resonant_symmetric_orbit` function) via ad hoc read-only Python,
not new files. `uv run pytest tests/data/test_outstanding_structure.py
tests/data/test_outstanding_header_body_consistency.py -q` run before
committing the `OUTSTANDING.md` update (see that commit's own message for the
result).
