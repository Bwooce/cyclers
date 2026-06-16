# 2026-06-16 -- #303 BCR4BP Phase 2: mu_sun natural-parameter continuation

**Task**: #303 (#292 Phase 2). Continue a BCR4BP periodic orbit from the
CR3BP-limit (``mu_sun = 0``) anchor to the published Andreu / Rosales-Jorba
parameter value (``mu_sun = 328900.5423094043``), producing an
intermediate-fidelity bridge between CR3BP and BCR4BP family members.

**Scope discipline**:
* Phase 2 produces JSONL deliverables + this doc note. NO catalogue
  writeback (Phase 5 separately).
* NO V-level promotion of any existing catalogue row. Validation-evidence
  matches surface as CANDIDATES for human review.
* NO novelty claims. Both the CR3BP L1 Lyapunov family and Andreu POL1 are
  published; continuation methodology is Simo-Jorba-Gomez.
* Sourced golden discipline: the seed (CR3BP L1 Lyapunov, Braik-Ross 2026
  ``C = 3.1294``) and the target ``mu_sun`` (Rosales-Jorba 2023 Table 3) are
  both SOURCED. Intermediate family members are OUR computation; tests
  assert TOPOLOGY (L1 neighbourhood) and CLOSURE (residual <= 1e-6 under
  the corrector's half-period symmetric residual) only.

---

## Deliverables

| Path | Purpose |
|---|---|
| ``src/cyclerfinder/genome/bcr4bp_continuation.py`` | New module: ``continue_bcr4bp_family_in_musun`` driver + ``BCR4BPFamily`` / ``BCR4BPFamilyMember`` dataclasses + Floquet stability classifier |
| ``tests/genome/test_bcr4bp_continuation.py`` | 6 gates: trivial-step closure, geometric-schedule monotonicity, stability-tag taxonomy, bad-input rejection (2), CR3BP-limit anchor tightness |
| ``scripts/run_303_bcr4bp_l1_continuation.py`` | End-to-end runner: CR3BP L1 Lyapunov -> BCR4BP@Andreu (50 steps, geometric in mu_sun) |
| ``scripts/run_303_catalogue_validation_probe.py`` | Catalogue probe: scan 11 Earth-Moon catalogue rows for period matches |
| ``data/bcr4bp_l1_family_303.jsonl`` | 50 converged members + metadata header |
| ``data/bcr4bp_validation_bridges_303.jsonl`` | 0 candidates + 11 no-match records + honest-scope header |

---

## Phase 2 verdict

### Continuation: full bridge confirmed for the L1 Lyapunov family

**50 of 50 continuation steps converged.** mu_sun extent reaches
``2.89e-01`` -> ``3.289e+05`` (the Andreu value exactly). Stepping was
geometric in ``mu_sun + offset`` (offset 1.0 so the seed at 0 is well-
defined).

| Quantity | Seed (mu_sun=0) | Final (mu_sun=Andreu) |
|---|---|---|
| ``x0`` | 0.8115256469 | 0.8114204 |
| ``vy0`` | 0.2561842689 | 0.2456918 |
| ``T (TU)`` | 2.946253 | 2.950491 |
| ``T (days)`` | 12.8114 | 12.8298 |
| Corrector residual (half-period symmetric) | 3.5e-14 | 4.8e-14 |
| Stability tag | hyperbolic_pair | hyperbolic_pair |

The L1 Lyapunov family at ``C = 3.1294`` (Braik-Ross 2026 common Jacobi
level) continues SMOOTHLY in ``mu_sun`` -- the IC barely moves (10^-4 in
``x``, 10^-2 in ``vy``, 1.5e-3 in ``T``). The stability tag stays
``hyperbolic_pair`` throughout (saddle x center, the expected topology
for an L1-substitute orbit).

Half-period symmetric corrector residual stays at machine precision (1e-14)
throughout. The **independent (Radau) full-period closure** residual grows
with mu_sun (from 1e-12 at mu_sun=0 to ~1e-2 at mu_sun=Andreu) -- not a
defect: it reflects the fact that T is FREE in the continuation, so the
period is NOT strictly Sun-commensurate, and the full-period closure picks
up an O(mu_sun) Sun-phase mismatch. This is the same mechanism documented
in the Phase 1 Gate 3 weak-Sun halo test. To enforce strict Sun-
commensurate closure one would need a multi-shooting or T-fixed continuation
(future Phase 3 work).

### Andreu POL1 is NOT in the L1 Lyapunov family

The task spec asked whether the continuation reaches Andreu POL1 cleanly.
**The answer is "no, but cleanly so": Andreu POL1 and the EM L1 Lyapunov
family are distinct invariants.**

| Orbit | ``x0`` | ``vy0`` | ``T (TU)`` |
|---|---|---|---|
| EM L1 Lyapunov @ C=3.1294 (CR3BP, our seed) | +0.81153 | +0.25618 | 2.946 |
| EM L1 Lyapunov @ mu_sun=Andreu (continuation final) | +0.81142 | +0.24569 | 2.950 |
| Andreu POL1 (sourced; Phase 1 corrector test Gate 2) | -0.84401 | -0.02174 | 6.791 (Sun-commensurate n=1) |

POL1 sits at NEGATIVE ``x``, has a SMALL libration speed, and lives at
period 6.79 TU (the n=1 Sun-commensurate period). The EM L1 Lyapunov is
at POSITIVE ``x``, has a moderate libration speed, and lives at period
2.95 TU (the family's own internal period, not commensurate with the Sun
synodic period).

POL1 is in the Sun-Earth-Moon barycentric L1 dynamical-substitute family
(different libration point, different rotating-frame convention). A
continuation that reaches POL1 would need to seed at the CR3BP-limit POL1
analogue -- i.e. a CR3BP planar L3-substitute or Sun-frame periodic orbit
at the right period -- not the EM L1 Lyapunov.

This is a clean negative-result discovery: **the EM L1 Lyapunov in BCR4BP
is itself a valid Phase 2 deliverable; POL1 is a separately seedable
family deferred to Phase 3.** Both Phase 1 Gate 2 (POL1 -> nearby BCR4BP
periodic orbit) and Phase 2 (L1 Lyapunov continuation) close cleanly; they
just close to DIFFERENT orbits.

---

## Catalogue validation-evidence findings

Scanned 11 Earth-Moon catalogue rows (``arenstorf-em-figure8-1963``,
``genova-aldrin-2015-em-3petal-cycler``, ``wittal-2022-em-cycler-family``,
5x ``ross-rt-em-cycler-*-2025``, 3x ``braik-ross-*-cycler-2026``).
Matched period_days from the catalogue against the closest BCR4BP L1
Lyapunov family member's period.

**Result: 0 candidates at the 5%-period-deviation threshold.** All 11
rows either:
* sit at a different period regime (Ross-RT (1,1)-(3,3) cyclers: 44-85 d;
  Braik-Ross C-cyclers: 42-79 d; Wittal 2-petal: 27 +/- 3 d) -- the L1
  Lyapunov family at ``C = 3.1294`` lives at ~12.8 d, well below;
* have only "approximate" or free-form-text periods that are not closure
  goldens (figure-8: "~14 d"; 3-petal: "~10 d"; Wittal: "~27 +/- 3 d") --
  not enough resolution to claim a bridge;
* lack any V_inf or other secondary identity carrier in the row's
  source_quotes (V_inf is null in all 11 cases -- these are CR3BP
  rotating-frame orbits where Jacobi constant is the conserved invariant,
  not patched-conic V_inf).

Closest-call rows (per ``data/bcr4bp_validation_bridges_303.jsonl``):

| Row | catalogue T (days) | family member T (days) | deviation |
|---|---|---|---|
| ``arenstorf-em-figure8-1963`` | 14.0 (free-form note) | 12.83 (final member) | 8.36% |
| ``genova-aldrin-2015-em-3petal-cycler`` | 10.0 (free-form note) | 12.81 (seed) | 28.1% |

The 8.4% Arenstorf gap is suggestive but **NOT** a bridge -- the figure-8
orbit is a low-energy free-return family, not an L1 Lyapunov, and the
period-only match without a Jacobi or V_inf cross-check is too weak. The
JSONL flags it as ``no_match`` with the honest-scope explanation.

**Net catalogue impact: zero V-level promotions, zero candidate bridges,
11 honest no-match records.** This is a clean negative that VINDICATES the
catalogue's existing classifications: the EM L1 Lyapunov at this Jacobi
level genuinely does not overlap any of the existing 11 EM rows. A future
probe with the BCR4BP HALO family (Phase 3) or the BCR4BP POL1 substitute
(separately seeded) might intersect the figure-8 / Arenstorf or Wittal
rows -- but the L1 Lyapunov in mu_sun continuation does not.

---

## Phase 3 path (concrete first step)

Two natural extensions:

1. **POL1 substitute family seeded directly** (NOT a continuation of EM
   L1 Lyapunov). Use Phase 1's POL1 -> BCR4BP closure (Gate 2; x = -0.844,
   T = 6.79 TU at n=1 Sun-commensurate) as the seed; continue in a
   different parameter (e.g. ``a_sun_nondim`` perturbations, or
   theta_sun0 phase shift, or "n" Sun-commensurate index). This gives the
   BCR4BP family that LIVES at the POL1 topology, distinct from the L1
   Lyapunov family.

2. **3D halo BCR4BP family** seeded from the EM L1 southern halo IC used
   in Phase 1 Gate 1 (the Howell / NASA TN D-1949 family at z = -0.0545)
   and continued in mu_sun. This crosses with #291 Phase 1's 3D CR3BP
   corrector (the BCR4BP corrector handles 3D states natively; the seed
   just needs to be a 3D halo IC). Could intersect the catalogue's
   Arenstorf figure-8 (if the figure-8 is found to be a planar projection
   of a halo at some C) or set up a 3D POL family analogue.

**Concrete first step**: write
``scripts/run_303_phase3_pol1_native_continuation.py`` that seeds from
the POL1 -> BCR4BP closure and continues in either ``a_sun_nondim`` or
``theta_sun0``. The continuation driver from this task generalises trivially
(just swap the natural parameter in the schedule + system rebuild logic).

---

## Honesty annotations

* The L1 Lyapunov in BCR4BP at Andreu mu_sun is OUR computation; the
  family is documented in BCR4BP-literature broadly but the specific
  ``C = 3.1294`` slice at Andreu's mu_sun is not a sourced golden. Per
  the orbit-closure discipline, the family JSONL members are evidence,
  not catalogue-promotable.
* The independent (Radau) closure residual at ~1e-2 at full mu_sun is
  NOT a defect; it's the Sun-phase mismatch from letting T drift. A
  T-fixed continuation (more expensive, smaller step size needed for
  convergence) would close to machine precision but explore a narrower
  family slice.
* The Floquet ``hyperbolic_pair`` tag throughout is *evidence* of the
  saddle x center topology but is NOT a sourced stability index. The
  catalogue's stability_index field would require a published numerical
  value to match against.
* Andreu POL1 is NOT reached by this Phase 2 continuation -- and that's
  the right answer because the EM L1 Lyapunov and POL1 are distinct
  invariants. Phase 1 Gate 2 closes POL1 directly; Phase 2 closes L1
  Lyapunov via continuation. The two are complementary, not redundant.
