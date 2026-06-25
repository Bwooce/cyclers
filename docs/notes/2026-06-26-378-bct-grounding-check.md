# #378 BCT grounding / acquisition decision gate (Task 0.1)

**Date:** 2026-06-26
**Source consulted:** `docs/notes/2026-06-17-digest-belbruno-2004.md` §3.4.1
(backwards-integration recipe), Remark 4 (forward 2×2 differential correction),
and the "BCT construction methodology (canonical Belbruno recipe)" section.

## Decision: **buildable from digest** (no acquisition required to start)

The forward 2×2 BCT recipe is fully specified at the granularity needed to map it
onto the existing BCR4BP Newton corrector
(`genome/bcr4bp_genome.correct_bcr4bp_periodic`) and the #380 augmented-quadrature
seam (`search/bvp_integral.propagate_augmented_bcr4bp` +
`correct_with_integral_constraints`). Every equation the substrate needs is in the
digest:

| Need | Digest anchor | Status |
|---|---|---|
| Kepler two-body energy `E₂ = ½|Ẋ|² − μ/r₂₃` (P₂-centred inertial) | Def 3.10 / eq 3.6 | present, closed-form |
| Ballistic capture = `E₂ ≤ 0` | Def 3.11 | present |
| WSB surface `W = J⁻¹(C) ∩ Σ ∩ σ` | Def 3.12 / eq 3.9 | present |
| Analytic `W`: `C = −r₂(±2√(μ(1+e₂)/r₂)+r₂) + μ(1−e₂)/r₂ + A(r₂,θ₂)` | Lemma 3.21 / eq 3.29 | present |
| Validity domain `C < C₁` (Earth-Moon C₁≈3.184) | Def 3.22 | present |
| Parabolic golden `C = ±√2 + O(μ)` | Lemma 3.34 / eq 3.39 | present, closed-form |
| Numerical stability-class algorithm (one-rev label) | §3.2.1 | present at algorithm granularity |
| Backward two-arc construction recipe | §3.4.1 / Fig 3.14 | present (control: QF on W, e₂≈0.95, r_M+100 km; back ~45 d to Q_a; fwd ~100 d from Q₀; match at Q_a) |
| Forward 2×2 corrector: control `(|V₀|, γ₀)` → target `(r₂₃, i_M)`, prescribe `(t₀, r₁₃, i_E)` | Remark 4 | present (control vars, target vars, Newton + 2×2 Jacobian all named) |
| Hiten signature (ΔV 44 m/s = 14+30, TOF ~150 d, Q_a ~1.5e6 km ≈ 3.9 LD, e₂≈0.95, ΔV_capture=0) | §3.4 | present, all values sourced |

## What is NOT in-repo (named acquisition fallbacks, NOT needed to start)

The convergence-engineering detail papers Belbruno cites for the forward method —
[33] Belbruno-Carrico 2000, [34] Belbruno-Humble-Coil 1997, [39]
Belbruno-Miller JGCD 16(4):770-775 1993 — are **not** in the corpus. They supply
operational targeting heuristics (initial-guess strategy, step damping), not new
governing equations. The 2×2 Newton on `(|V₀|, γ₀) → (r₂₃, i_M)` is a standard
finite-difference/STM single-shooting solve, which the existing BCR4BP corrector
already implements; the substrate is a thin assembler over it.

**Specific named acquisition IF the forward 2×2 corrector stalls (R3):**
[39] Belbruno & Miller, JGCD 16(4):770-775, 1993 — the published Hiten BCT
operational-targeting detail. OCR → digest → CORPUS_INDEX per corpus-document
policy before adoption. Not pulled now: no missing governing equation.

## Conclusion

Proceed with Phase 0.2 (apoapsis-reach feasibility spike) and Phase 1
(`core/wsb.py`) directly. The honest model gap (project incoherent BCR4BP vs
Belbruno's PR4BP-3D-with-DE403) makes the Hiten golden a **signature band**, not a
bit-exact reproduction — already documented in the design draft §4 and consistent
with the BCR4BP O(ε²) POL1 gap.
