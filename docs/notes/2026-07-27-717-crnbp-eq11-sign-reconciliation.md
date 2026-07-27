# #717 Step 0: Negri & Prado (2022) Eq. 11 coupling-term sign reconciliation

**Task:** `#717`, shortlist item 1 of 3 from `#714`'s CCR5BP/CRNBP scoping pass
(`docs/notes/2026-07-27-714-ccr5bp-crnbp-discovery-strategy-pass.md`). Mandatory
Step 0, to be resolved from source PDFs before any code was written for the N=5
Jupiter-Europa-Io-Ganymede CRNBP core module (`src/cyclerfinder/core/crnbp.py`).

## The discrepancy

Two prior digests transcribed Negri & Prado (2022)'s general N-body EOM Eq. 11
(the "inner sum" mutual coupling term among extra perturbing bodies) with
**opposite signs**:

- `docs/notes/2026-07-26-digest-negri-prado-2022-crnbp.md` (§2), the digest of
  Negri & Prado's own paper:

  ```
  x'' = ... - sum_j mu_j [ (x+mu2-Rj*cos psi_j)/rj^3
                            - sum_{k!=j} mu_k (Rj*cos psi_j - Rk*cos psi_k) / rkj^3 ]
  ```

  an inner **MINUS** between the direct term and the coupling sum.

- `docs/notes/2026-07-26-digest-gilliam-bettinger-2024-crnbp-jovian.md` (§1),
  the digest of Gilliam's independently-derived thesis EOM:

  ```
  x'' = ... - sum_j mu_j [ (x+mu2-Rj*cos psi_j)/rj^3
                            + sum_{k!=j} mu_k (Rj*cos psi_j - Rk*cos psi_k) / rkj^3 ]
  ```

  an inner **PLUS**.

Both digests claimed term-for-term agreement with the general Negri-Prado
form, so at least one transcription had a sign slip.

## Resolution: read the source PDFs directly (not guessed, not split)

Sources read directly, in `/Users/bruce/dev/cyclers_pdf/papers/` (private
corpus, read-only reference for this task):

- `negri-prado-2022-circular-restricted-n-body-problem-jgcd-doi-10.2514-1.G006430-arxiv-2307.10881.pdf`
- `gilliam-2025-crnbp-multibody-systems-thesis-afit-etd-8309.pdf`

**Negri & Prado's own Eqs. 7-9** (the pre-simplified, inertial-frame form,
*before* their Eq. 11 substitution) are textually UNAMBIGUOUS — no
multi-line-bracket rendering ambiguity, unlike their own Eq. 11a/b (which IS
genuinely hard to parse correctly from a naive `pdftotext -layout` extraction,
and is the likely source of the bad digest's error). Eq. 8, transcribed
directly:

```
p_N'' = G { M1(p1-pN)/|p1-pN|^3 + M2(p2-pN)/|p2-pN|^3
            + sum_{j=3}^{N-1} Mj [ (pj-pN)/|pj-pN|^3
                                    + sum_{k=1,k!=j}^{N-1} (Mk/(M1+M2)) (pk-pj)/|pk-pj|^3 ] }
```

The direct term and the inner-sum coupling term are **ADDED** (same sign)
inside the `Mj[...]` bracket — an inner PLUS. Substituting Eq. 10's
`rho_j - rho_N = -(x+mu2-Rj*cos(psi_j), y-Rj*sin(psi_j), z)` and
`rho_k - rho_j = -(Rj*cos(psi_j)-Rk*cos(psi_k), Rj*sin(psi_j)-Rk*sin(psi_k), 0)`,
normalising by `G(M1+M2)=1`, and moving to the rotating frame reproduces
Eq. 11 with an inner PLUS inside the bracket and an overall MINUS in front of
`sum_j mu_j[...]` — i.e. **Gilliam's digest transcription is correct**; the
`2026-07-26-digest-negri-prado-2022-crnbp.md` digest has the transcription
error.

Independently confirmed twice more:

1. **Gilliam's own thesis, Eqs. 25-27** (a clean, native text-layer PDF — no
   OCR/bracket-rendering issue at all): same inner PLUS, matching the Negri &
   Prado Eq. 7-9 derivation term-for-term. The thesis states this derivation
   was done independently via Lagrangian mechanics (not copied from Negri &
   Prado), so this is a genuinely independent confirmation.
2. **A second, independent source bound into the same thesis PDF**: an
   appended reprint of Gilliam, Bettinger et al., *Icarus* 429 (2025) 116455
   ("debris propagation... using the circular restricted 3- and N-body
   problems"), Eqs. 9-11 — same inner PLUS again, a third independent
   transcription of the same physics by (mostly) the same authors in a
   different, later publication.

Three independent sources, all agreeing: **the digest of Negri & Prado's own
paper (`2026-07-26-digest-negri-prado-2022-crnbp.md`) has the error; the
Gilliam digest was correct.** The error is best explained as an artefact of
that PDF's Eq. 11a/b multi-line bracket typesetting, which genuinely does not
extract unambiguously under a naive `pdftotext -layout` read (the bracket
open/close glyphs land on different visual rows than the terms they
delimit) — exactly the class of hazard `#714` itself flagged before
dispatching this task's Step 0.

## The correct term

```
a_j = -mu_j * [ (r - r_j)/|r-r_j|^3 + r_j/a_j^3
                 + sum_{k!=j} mu_k * (r_j - r_k) / |r_j-r_k|^3 ]
```

(direct + indirect, exactly `ccr4bp.py`'s existing single-perturber term,
PLUS the new mutual coupling term — all with the SAME sign, all subtracted
from the base CR3BP acceleration.)

## A further finding, beyond the sign question (Step 1, `#717`)

While implementing this in `src/cyclerfinder/core/crnbp.py`, a structural
property of this correctly-signed term was found and verified three
independent ways (see that module's docstring for the full derivation and
proof): **the mutual coupling term, however many extra perturbers there are,
contributes EXACTLY ZERO to the spacecraft's total acceleration** — an
antisymmetric pairwise cancellation intrinsic to the formula (for any
unordered pair of extra bodies `{j,k}`, the ordered contributions `(j,k)` and
`(k,j)` are exact negatives of each other, so the sum over all pairs is
always zero, regardless of masses, positions, or N). Checked: (1) the
algebraic argument above; (2) numerically against the physical Io+Ganymede
N=5 system; (3) a from-scratch, independent re-transcription of Gilliam's
Eq. 25 (not reusing any of `crnbp.py`'s code) at both N=5 and N=6 (three
extra bodies).

This refines (does not invalidate) `#714`'s own tractability assessment,
which had characterized the coupling term as "genuinely new N>=5 physics"
that could NOT be captured by "just call[ing] `_ganymede_acceleration`
twice." That specific warning turns out to be incorrect once the term is
correctly summed: a naive superposition of independent single-perturber
(`ccr4bp`-style direct+indirect) terms is mathematically IDENTICAL to the
full, correctly-signed Eq. 11 CRNBP acceleration, for any N. The genuinely
new N=5 physics is simply that the spacecraft now feels TWO independent
periodic forcings simultaneously (Io's own, Ganymede's own — the latter
unchanged from `#690`'s existing CCR4BP work), not a mutual coupling effect.
`src/cyclerfinder/core/crnbp.py` still implements the coupling term
faithfully (rather than silently omitting it) so that this cancellation is a
TESTED structural property (`tests/core/test_crnbp.py`) rather than an
unverified assumption baked into a simplified implementation.
