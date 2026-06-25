# QP 2-Tori GMOS Family Continuation (#333 / #290 Phase 2) — Design Draft

**Date:** 2026-06-25
**Status:** DRAFT — scoping only, no implementation. User reviews before build.
**Lineage:** #290 Phase 1 built the single-torus GMOS genome
(`genome/qp_tori.py` — `QPTorus`, `correct_qp_torus`, `_correct_gmos`, the
Neimark-Sacker seeder, the irrationality test) plus a **natural-parameter
amplitude** continuation stub (`genome/qp_tori_continuation.py` —
`continue_qp_family`, explicitly *not* pseudo-arclength, stops at folds). #319
adapted the V0-V5 gauntlet to tori (`data/validation/v1_qp.py`,
`v2_qp.py`). #320 ran the first sweep (12 Neimark-Sacker brackets off the #296
Braik-Ross (1,1) Earth-Moon family via #299's tracker) → 2 SILVER Earth-Moon
tori, no catalogue writeback (lit-check has no QP fingerprint yet). This Phase 2
turns the *single converged torus* into a properly-continued *family of tori*.

## Goal

Given one converged 2-torus (an invariant circle under the stroboscopic map),
**continue it into a 1-parameter family**, tracking the family across folds the
way `cr3bp_3d_family_tracer.continue_general_3d_family` and
`er3bp_continuation.continue_er3bp_family_in_e_arclength` continue periodic
orbits — i.e. a real **pseudo-arclength** walker over the torus unknowns, not the
fold-blind amplitude stepper that exists today. Report-only; no catalogue
writeback in this task (matches the #320 / #430 discipline — the QP lit-check
fingerprint and a V3+_qp gauntlet are separate, deferred work).

## 1. Phase-1-built vs Phase-2-needed — the precise gap

**Phase 1 built (single torus):** `_correct_gmos(system, x0, n_trans, n_samples,
phase_pin_idx, amplitude_pin, tol, max_iter) -> (x_final, residual, n_iter)`
solves the GMOS invariance equation `phi_{t_strob}(u(θ)) = u(θ+rho)` in
Fourier-mode space for **one** torus. Unknown vector (`_pack_unknowns`):
`c_0` real (6) + `c_1..c_N` complex (12N) + `rho` (1) + `t_strob` (1) =
`6 + 12N + 2`. Closed by two pin rows (`_residual_real`): a **phase pin**
`Im(c_1[phase_pin_idx]) = 0` (kills rotation invariance) and an **amplitude pin**
`|c_1| − amplitude_pin = 0` (kills the trivial `c_n→0` collapse). Jacobian is
**finite-difference** (`scipy.optimize.least_squares`, `trf`, `diff_step=1e-5`) —
the corrector docstring names this the accuracy bottleneck (residual bottoms
~1e-7 at N=2) and explicitly defers an analytic Jacobian to "Phase 2".

**Phase 1 also built (amplitude stub):** `continue_qp_family(...)` — secant
predictor in the packed `x` vs amplitude, corrector = `_correct_gmos` with
`amplitude_pin = amp`. It holds amplitude **fixed** at each step and **cannot
traverse folds** (records `fold_detected`, stops). It relaxes the gate to
`max(tol, 1e-5)` to avoid mistaking FD-Jacobian noise for a fold.

**Phase 2 needs (the family continuator) — five concrete deltas:**

1. **Promote the continuation parameter into the unknown vector.** Today amplitude
   is a *fixed pin*; the reference walkers embed the parameter (`T` in the 3D
   tracer, `e` in ER3BP) **inside** `z` and replace the fixed-pin row with the
   arclength row. The QP unknown vector must gain a free continuation parameter
   `λ` so the solver can move *along* the family, including around folds where the
   parameter is non-monotone.
2. **Replace the amplitude-pin row with a pseudo-arclength row** `tau·(z − z_pred)
   = 0`, mirroring `_correct_pseudo_arclength` / `_correct_arclength`. The **phase
   pin stays** (rotation invariance is intrinsic to the torus parameterization,
   present at every family member, not a continuation device).
3. **SVD null-space tangent predictor** over the (now rank-deficient by exactly
   one) torus Jacobian — the analogue of `_tulip_tangent` / `_arclength_tangent`
   (last right-singular vector, oriented by `+dot` with the previous tangent).
4. **Fold detection + traversal** (sign change of the tangent's parameter
   component), recorded as a `FoldPoint`-analogue — the capability the amplitude
   stub structurally lacks.
5. **An analytic (or at least variational) Jacobian** of the GMOS residual w.r.t.
   the augmented unknown vector. The FD Jacobian noise floor (~1e-7) currently
   *masquerades as folds*; pseudo-arclength fold detection is unreliable on top of
   FD noise. This is the long-deferred Phase-2 Jacobian, now load-bearing.

## 2. Recommended continuation parameter — **Jacobi-energy (`C_J`) primary**

Choices considered: amplitude `|c_1|` (current), rotation number `rho`, stroboscopic
period `t_strob`, longitudinal frequency `omega_long`, Jacobi-analogue energy `C_J`.

**Recommendation: continue in the Jacobi constant `C_J` as the primary family
parameter, embedded in `z` as a free unknown with `C_J(c_0, ...)` appended as an
explicit constraint row** — for three reasons:

- **Physical family coordinate, golden-comparable.** The #320 tori are already
  reported by `C_J` (SILVER Bracket 2: `C_J=3.03196`; Bracket 10: `C_J=3.12624`;
  #290 smoke: `C_J=3.12785`). A `C_J`-parameterized family is directly
  cross-checkable against any future sourced Olikara/Howell torus family (those
  are published as energy/Jacobi families) and against the CR3BP family the torus
  bifurcated from (same energy coordinate).
- **Monotone-ish but fold-safe.** Energy typically varies monotonically along a
  torus family but can fold; pseudo-arclength embedding (not fixed-`C_J` natural
  stepping) handles the folds, so we get the best of both.
- **Avoids the resonance trap of `rho`.** Continuing directly in `rho` walks
  straight into Arnold tongues / phase-locking (Risk §1); `rho` should be
  **monitored**, not driven. Energy continuation lets `rho` evolve freely and we
  *detect* lock-in as a diagnostic rather than fighting it as a constraint.

`rho`-continuation is offered as an **optional secondary mode** (`param="rho"`) for
deliberately sweeping toward a resonance to study a tongue, but it is not the
default and carries the Arnold-tongue handling burden directly.

## 3. Architecture — new module `genome/qp_tori_arclength.py`

Mirrors `cr3bp_3d_family_tracer` / `er3bp_continuation` structure exactly, reusing
the Phase-1 GMOS residual machinery unchanged where possible.

**Augmented unknown vector** `z`:
```
z = [ pack_unknowns(c_0, c_{1..N}, rho, t_strob) , C_J ]      # length 6+12N+2+1
```
**Residual** `R(z)` (rows):
```
F_n            GMOS invariance, |n| ≤ N           (real+imag, masked tail)   [Phase 1 _gmos_residual]
phase pin      Im(c_1[phase_pin_idx]) = 0                                     [Phase 1, kept]
energy row     jacobi(c_0, system) − C_J = 0      (ties the free C_J to state)
arclength      tau · (z − z_pred) = 0                                         [NEW]
```
Note the amplitude pin is **dropped** — the arclength row supplies the missing
constraint, and `|c_1|` is now free to vary along the family (that *is* the
family). The energy row makes `C_J` a genuine readable coordinate rather than an
implicit function of the modes.

### New signatures

```python
@dataclass(frozen=True)
class QPTorusFamilyMember:
    torus: QPTorus                 # the converged Phase-1 representation
    jacobi: float                  # C_J family coordinate
    arclength_s: float             # accumulated pseudo-arclength
    tangent: NDArray               # unit tangent at this member
    rho: float                     # monitored rotation number
    freq_ratio: float              # omega_trans/omega_long (resonance monitor)
    is_practically_irrational: bool
    near_resonance: ResonanceFlag | None   # p:q, distance, tongue-width estimate
    fold_index: int | None         # set if this member is a fold point
    residual_norm: float
    extras: dict[str, float]

@dataclass(frozen=True)
class QPTorusFold:                 # analogue of FoldPoint / cr3bp tracer
    member_index: int
    param_at_fold: float           # C_J at the fold
    tangent_param_component: float # the component that changed sign

@dataclass(frozen=True)
class QPFamily:
    members: list[QPTorusFamilyMember]
    folds: list[QPTorusFold]
    resonance_crossings: list[ResonanceCrossing]   # p:q tongues encountered
    terminated_reason: str         # "max_steps" | "corrector_fail" |
                                   # "resonance_lock" | "mode_truncation_breach"
    seed_torus_id: str
```

```python
def _gmos_residual_and_jac(
    z: NDArray, system: CR3BPSystem, n_trans: int, n_samples: int,
    phase_pin_idx: int, *, analytic: bool = True,
) -> tuple[NDArray, NDArray]:
    """GMOS invariance + phase-pin + energy residual AND its Jacobian w.r.t. the
    augmented z (modes, rho, t_strob, C_J). Analytic block-structured Jacobian
    via variational propagation of phi_{t_strob} over the sample points
    (per-sample STM) + FFT linearity; FD fallback (analytic=False) for parity
    testing. Replaces the least_squares FD Jacobian of Phase 1."""

def _arclength_tangent(jac: NDArray, prev: NDArray | None) -> NDArray:
    """Unit null-space tangent = last right-singular vector of the SVD of the
    (rank-1-deficient) residual Jacobian; oriented by +dot with prev.
    Mirrors cr3bp_3d_family_tracer._tulip_tangent / er3bp _arclength_tangent."""

def _correct_arclength_torus(
    z_pred: NDArray, tau: NDArray, system: CR3BPSystem, *,
    n_trans: int, n_samples: int, phase_pin_idx: int,
    tol: float, max_iter: int = 60,
    mode_cap: float = 0.1, rho_cap: float = 0.05, cj_cap: float = 1e-2,
) -> NDArray | None:
    """Newton onto {R(z)=0, tau·(z−z_pred)=0} on the augmented Jacobian
    [dR/dz; tau], np.linalg.solve with lstsq fallback + per-block step caps.
    Mirrors _correct_pseudo_arclength / _correct_arclength."""

def continue_qp_family_arclength(
    seed_torus: QPTorus, *,
    param: Literal["jacobi", "rho"] = "jacobi",
    ds: float = 5e-3, max_steps: int = 200, direction: Literal["both","fwd","rev"] = "both",
    corrector_tol: float = 1e-8, n_trans: int | None = None,
    fold_detection: bool = True,
    resonance_max_denominator: int = 12, resonance_tol: float = 1e-3,
    mode_truncation_guard: float = 1e-4,    # |c_N| tail-energy ceiling
    on_step: Callable[[QPTorusFamilyMember], None] | None = None,
) -> QPFamily:
    """Pseudo-arclength continuation of a converged 2-torus into a family.
    Predictor: z_pred = z_cur + ds·tau (SVD null tangent). Corrector:
    _correct_arclength_torus. Monitors freq_ratio for p:q lock-in (Arnold
    tongues) and |c_N| for mode-truncation breach; bumps n_trans or halves ds
    on breach. Walks BOTH directions from the seed by default. Incremental
    on_step callback for detached-run progress logging."""
```

Reused unchanged: `QPTorus`, `_pack_unknowns` / `_unpack_unknowns`,
`_enforce_reality`, `evaluate_invariant_circle`, `is_practically_irrational`,
`correct_qp_torus` (to produce the seed), `CR3BPSystem.jacobi`. The legacy
`continue_qp_family` (amplitude stub) stays as a fast natural-parameter mode for
fold-free local sweeps; it is **superseded, not deleted**.

## 4. Data flow

```
converged seed torus (correct_qp_torus off a #299 NS bracket)
  → build augmented z0 (pack modes,rho,t_strob + C_J)
  → _gmos_residual_and_jac → _arclength_tangent (seed tangent, both directions)
  → loop (each direction):
        z_pred = z_cur + ds·tau
        _correct_arclength_torus → z_next  (or None → halve ds / terminate)
        unpack → QPTorus → QPTorusFamilyMember
        monitor: freq_ratio (resonance), |c_N| (truncation), fold sign-change
        on_step(member)   # incremental JSONL
  → QPFamily {members, folds, resonance_crossings, terminated_reason}
  → harvest note + data/family_333_qp_*.jsonl
  → (irrational, V1_qp+V2_qp-passing members) → candidate family rows  (NOT written here)
```

## 5. Bite-sized TDD build sequence

Each step is one red→green cycle with a sourced or self-consistent assertion.

1. **Analytic Jacobian parity.** `_gmos_residual_and_jac(analytic=True)` matches
   `analytic=False` (FD) to `<1e-6` on the #290 smoke torus. *Golden: the FD
   Jacobian itself (parity), the standard variational-derivative guard.* This
   unblocks reliable fold detection (Risk §2).
2. **Tangent + reduces-to-corrector.** `_arclength_tangent` returns a unit vector
   in the 1-D null space; a single `_correct_arclength_torus` step with `ds=0`
   reproduces `_correct_gmos` on the seed (the new corrector is a strict
   generalization).
3. **Forward family step.** From the #290 smoke torus, one arclength step
   produces a *new* converged torus with `C_J` shifted by ≈`ds·tangent_Cj`,
   `invariance_residual < corrector_tol`, and `is_practically_irrational` still
   true. *Self-consistent capability golden.*
4. **Both-directions + monotone energy.** N steps each way yield a family whose
   `C_J` is monotone away from the seed (absent a fold), and whose endpoints
   bracket the seed energy.
5. **Fold traversal.** On a family seeded near a known fold (constructed by
   continuing toward decreasing `n_trans` headroom), the walker crosses the fold
   (records a `QPTorusFold`, `C_J` non-monotone, continues past) where the legacy
   `continue_qp_family` stops. *This is the headline Phase-2 capability test.*
6. **Resonance monitor.** A `rho`-mode continuation deliberately steered toward a
   low-order `p:q` flags `near_resonance` and emits a `ResonanceCrossing` (does
   not silently produce a phase-locked "torus"); cross-checks the #320 screened
   1:4 partner bracket as the known lock-in case.
7. **Mode-truncation guard.** A family pushed until `|c_N|` exceeds
   `mode_truncation_guard` either auto-bumps `n_trans` or terminates with
   `terminated_reason="mode_truncation_breach"` (never reports a member whose tail
   energy invalidates the truncation).
8. **Each family member passes V1_qp.** Every converged member fed to
   `run_v1_qp` PASSes (Fourier-norm `<1e-5`, off-grid `<1e-4`) — the family is
   genuinely a chain of valid tori, gauntlet-consistent.
9. **Determinism.** Fixed seed + `ds` + `max_steps` → reproducible member
   ordering and `C_J` sequence.

### Golden / validation target

**No fully-tabulated published QP-torus family is harvestable** from the corpus
(confirmed: the #290 note says the project's tori "sit on a continuous family
that Olikara 2016 would have populated" but quotes no member values; the
Andreu/Rosales-Jorba QBCP digest explicitly notes halo/QP-torus ICs with rotation
numbers are *not* tabulated with full coordinates in the open papers — "would need
the thesis"). Therefore the validation target is **two-tier**:

- **Capability golden (primary, available now):** the project's own #290 smoke
  torus + the two #320 SILVER Earth-Moon tori. A correct walker (a) reproduces
  each as a family member, (b) connects the #290 smoke torus and SILVER Bracket 2
  if they lie on one family (both near `C_J≈3.03–3.13`, both k=4 Earth-Moon —
  *a hypothesis the continuation itself tests*), and (c) at the family's
  low-amplitude end limits back onto the #299 Neimark-Sacker bifurcation point
  (the torus shrinks to the parent periodic orbit as `|c_1|→0`) — this
  **bifurcation-limit check is a strong, physically-sourced consistency anchor**
  even without external member tables.
- **Sourced golden (deferred, flagged as an acquisition):** an Earth-Moon torus
  family from Olikara-Scheeres 2010 / Olikara 2016 (Purdue) / Howell-Howell 2014 /
  Henderson-Howell 2008. None are digested with member coordinates yet; this
  Phase-2 build should **emit a corpus-acquisition follow-on task** (per the
  never-give-up-reproducing-papers discipline) rather than block on it. The
  capability goldens make the build testable today; the sourced golden upgrades
  it to a reproduction when a family table is acquired.

## 6. Risks

**(1) Rotation-number resonances / Arnold tongues — THE KEY RISK.** As a torus
family is continued, `rho` drifts; whenever `omega_trans/omega_long` passes through
a low-order rational `p:q`, the invariant circle **phase-locks** and the GMOS
problem degenerates (the 2-torus collapses onto a `q`-fold cover of a periodic
orbit — a resonance / Arnold tongue). The corrector can *converge* there to a
spurious "torus" that is actually a periodic orbit (exactly the #320 screened-out
1:4 partner). **Mitigation:** (a) continue in **energy, not `rho`** (recommendation
§2) so we cross tongues transversally instead of driving into them; (b) at every
member compute `is_practically_irrational(freq_ratio, max_denominator=12)` and emit
a `ResonanceCrossing` when within `resonance_tol` of a low-order rational — flag,
log, and *step over* the tongue (locally shrink `ds` to resolve it, then resume)
rather than reporting locked members as tori; (c) never admit a member that fails
the irrationality test to the family's "valid torus" set. Tongues that are too wide
to step over cleanly terminate the run with
`terminated_reason="resonance_lock"` — an honest boundary, not a silent collapse.

**(2) Fourier-mode truncation.** `N=2` (the #320 working point) is thin; as the
torus grows along the family the high-frequency content (`|c_N|`) rises and the
fixed-`N` truncation silently invalidates the representation. **Mitigation:** the
`mode_truncation_guard` monitors tail energy `|c_N|/|c_1|`; on breach, auto-bump
`n_trans` (re-pack `z`, re-solve) or terminate with
`mode_truncation_breach`. Pairs with the analytic Jacobian — FD noise at small `N`
otherwise hides the tail growth.

**(3) Per-step re-solve cost.** Each family step is a full GMOS Newton solve;
Phase 1's FD Jacobian makes one solve already slow (FD dominates, the documented
bottleneck), and a family is hundreds of solves. **Mitigation:** the analytic
(variational) Jacobian (build step 1) is the primary lever — it both speeds each
iteration and de-noises fold detection (so it is not optional polish, it is on the
critical path). Per-block step caps + secant warm-starts (predictor seeds the
corrector close) keep iteration counts low. Long runs are acceptable
(per `feedback_long_runs_acceptable`); the `on_step` callback gives incremental
checkpointed JSONL so a multi-hour family walk is monitorable and resumable, not a
black box (`feedback_incremental_progress_reports`).

## 7. Decomposition (for writing-plans)

~9 TDD tasks tracking §5: (1) analytic GMOS Jacobian + FD parity;
(2) tangent + corrector-generalization; (3) single forward step;
(4) both-directions energy monotonicity; (5) fold traversal;
(6) resonance monitor + `ResonanceCrossing`; (7) mode-truncation guard;
(8) per-member V1_qp consistency + determinism; (9) campaign runner
`scripts/run_333_qp_family.py` (seed off #290/#320 tori, detached, incremental
JSONL → harvest note) + corpus-acquisition follow-on task for the sourced golden.
Estimated ~700–1000 LOC across `genome/qp_tori_arclength.py` + tests
(`tests/genome/test_qp_tori_arclength.py`).

## What this draft does NOT decide (for review)

- Whether to also build a V3_qp/V4_qp tier for family members (currently V0–V2_qp
  only exist; family members are validated to V2_qp here, no new tier built).
- Whether the #290 smoke torus and #320 SILVER Bracket 2 actually lie on one
  family (the continuation tests this; not assumed).
- Catalogue writeback policy for a *family* of tori (one row? a family row with a
  member table?) — deferred to a post-build taxonomy decision, gated on the QP
  lit-check fingerprint extension (Olikara/Howell anchors) which is separate work.
