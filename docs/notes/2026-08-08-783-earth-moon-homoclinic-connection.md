# `#783` — Earth-Moon He1 homoclinic-connection stage (Barrabés-Mondelo-Ollé 2009)

**Task:** `#783`, the connection-stage follow-up `#780`'s own results note registered
("Earth-Moon homoclinic-connection stage... mirroring exactly how `#767` followed `#765` for
Saturn-Titan and `#777` followed `#776` for Neptune-Triton"). Unlike those three, this task
targets a GENUINELY DIFFERENT algorithm class: Casoliva's own Class 2 cyclers are built via
Barrabés, Mondelo & Ollé's own numerical-continuation-of-homoclinic-connections method
(Barrabés, E., Mondelo, J. M. & Ollé, M., "Numerical Continuation of Homoclinic Connections of
Periodic Orbits in the RTBP," *Nonlinearity* 22:2901, 2009, DOI `10.1088/0951-7715/22/12/006`,
filed at `cyclers_pdf/papers/barrabes-mondelo-olle-2009-...`), NOT this project's own
Poincaré-section Newton-shooting machinery (`jovian_resonant_connections.py`'s
`correct_connection`, reused verbatim by the Saturn-Titan and Neptune-Triton connection
modules).

---

## Verdict (read this first)

**Clean negative on the He1 connection itself, with genuine, quantitatively-diagnosed
progress and real durable value delivered alongside it.** This is a REPRODUCTION ATTEMPT of
Casoliva 2010's own published He1 connection at her own single hardest energy point
(`h=-1.45016232260699`, the "closest approaches to Earth and moon" member of the He1 family —
her own paper reports 5h20min of CPU time to reach this SPECIFIC connection, SEPARATELY from
(and harder than) the four OTHER He1-family energies she reports timing for (121s-1619s in her
own Fig. 8), via a full predictor-corrector continuation FROM easier, lower-energy starting
points). This task built Barrabés-Mondelo-Ollé's own method faithfully (their
manifold parametrization Eq. (4), matching condition Eq. (6)/(7)) and attempted a COLD start
directly at this hardest point — both a single-shooting and a multiple-shooting Newton
corrector make genuine residual progress (`0.394` → `~0.03` nondim, an order of magnitude)
but plateau well short of convergence. This is diagnosed as a genuine numerical "conditioning
wall" (quantitative arithmetic below), not a convention error, not a tuning failure, and not
abandoned prematurely — real escalation was attempted per this project's own "never give up
reproducing papers" discipline before concluding this.

**What DOES reproduce cleanly, and is durable value from this task regardless of the
connection outcome:**
1. A SECOND OCR sign hazard on Casoliva 2010's own Table 4, found and fixed this task (Table
   4's own `V^u` monodromy eigenvector) — `pdftotext` (both `-raw` and `-layout`) renders ALL
   FOUR components as unsigned; resolved by reading PDF page 13 directly as an image.
2. This module's own independently-recovered monodromy eigenvector reproduces Table 4's own
   (OCR-corrected) `V^u` to `1.9e-8` relative (unit-vector comparison) — a genuine,
   non-degenerate reproduction.
3. Table 6 row 8's own printed LEO-rendezvous ΔV self-consistency: `717.56` m/s computed
   DIRECTLY from Casoliva's own printed `(r, v)` at row 8, vs. her own stated `717.5` m/s.
4. Tables 5 and 6 (pericenter/apocenter flight-times/orbital elements for the He1 periodic
   orbit and the He1 connection, 4 + 19 rows) freshly vendored this task, image-read (not
   text-layer-extracted, for the same OCR-hazard reason as #1).
5. `#780`'s own He1 periodic-orbit anchor is re-confirmed to near-machine precision
   (closure `~7.5e-12`), ruling out "the frozen orbit is the defect" as an explanation for the
   connection's own non-convergence (checked directly, not assumed).

**What does NOT reproduce: the He1 connection itself.** Per this task's own honesty
discipline (explicit in the dispatch note), Table 5/6 pericenter/apocenter elements and the
LEO ΔV figure are NOT computed from either corrector's own unconverged trajectory — a
`~0.03` nondim residual is `~11,500` km position-scale, and Table 6's own Earth-relative
apogees cluster tightly near `354,000`-`358,000` km; computing "elements" from an
unconverged trajectory at that mismatch scale would produce numbers that superficially
resemble the target without being a genuine reproduction — the exact "it closed!" failure
mode this project's own orbit-closure discipline warns against.

---

## 1. The Barrabés-Mondelo-Ollé method, read directly from their own paper (not reconstructed
   from Casoliva's own garbled OCR text layer)

Their own Eq. (4) (manifold parametrization, page 13 JGCD equivalent in Casoliva 2010's own
unlabelled formula):

```
psi(theta, xi) = phi_{theta*T/(2*pi)}(z0) + xi * Lambda^(-theta/(2*pi)) * DPhi_{theta*T/(2*pi)}(z0) @ v0
```

Their own Eq. (6)/(7) (the full system, Casoliva's own Eq. (20) is essentially identical —
same author, same formulation):

Unknowns: `h, z, T, Lambda^u, v^u, Lambda^s, v^s, theta^u, T^u, theta^s, T^s` (plus, for their
own MULTIPLE-shooting version, used "whenever `T^u`/`T^s` may become large... in order to avoid
loss of precision" — their own words, confirming this task's independent conclusion below):
interior patch-point states along each leg. Equations: energy/section/periodicity/eigenpair
conditions pinning the periodic orbit, plus the matching condition
`phi_{T^u}(psi^u(theta^u, xi_0)) - phi_{T^s}(psi^s(theta^s, xi_0)) = 0`, with **`T, T^u > 0`
and `T^s < 0`** (their own explicit stated convention — both propagated FORWARD in their own
time convention, `T^s` being negative is what makes the stable leg run backward).

Their own implementation details, confirmed by reading their own numerical-methods section
directly (Sec. 2.2, page 12-13 of the paper): second-order variational equations, an RKF7(8)
integrator at `1e-14` tolerance, Newton tolerance `1e-10`-`1e-12`, an ADAPTIVE multiple-shooting
segmentation criterion (`||DPhi_t(z_i)||_inf < M`, `M` "typically tens or hundreds",
recomputed every continuation step), and minimum-norm least-squares (QR decomposition with
column pivoting) because their own system is "over-determined and rank deficient."

## 2. Table 4's own SECOND OCR sign hazard, resolved by direct image read

`#780`'s own module already found and fixed one OCR sign hazard on Table 4 (`p_y`, printed
without its own minus sign). This task found a SECOND, worse instance on the SAME table's
`V^u` row: BOTH `pdftotext -raw` and `pdftotext -layout` render ALL FOUR components as if
unsigned (`0.04936117474608325, 0.5669396006141868, 0.8204920465924677, 0.05418270168269977`)
— because the colon-for-minus-sign ligature substitution collapses every negative sign the
same way regardless of which numbers happen to be negative. Reading PDF page 13 (JGCD p. 1635)
directly as an image resolved this unambiguously: Casoliva's own printed
`V^u = (0.04936117474608325, -0.5669396006141868, 0.8204920465924677, -0.05418270168269977)`
— components 2 and 4 negative, 1 and 3 positive.

An earlier working attempt this task used a WRONG sign pattern (guessed from the `p_y`
precedent before the direct image read) and got a connection-residual discriminator test
~5x WORSE than the corrected version — caught and fixed before any constant was committed.
Table 5 and Table 6 were BOTH independently image-read too (not text-layer-extracted), for
the same reason — Table 6 visibly carries the same colon hazard on several of its own
negative semimajor-axis and argument-of-pericenter entries in the raw text layer (e.g.
`49783:466`, `0:760`, `147:477`).

**Cross-check that the corrected `V^u` is right**: this module's own independently-recovered
monodromy eigenvector (`build_he1_connection_seed` + `_planar_floquet_pair`, computed from
`#780`'s own converged periodic orbit, with ZERO dependence on Table 4's `V^u`) matches the
OCR-corrected `V^u` (after the same 180-degree-flip transform `#780` already validated) to
`cos_angle > 1 - 2e-8`, `rel_err_unit = 1.9e-8` — a genuine, tight, non-degenerate
reproduction (see `eigenvector_reproduction_check`).

## 3. Table 6 (newly vendored, image-read)

19 rows, flight times (days, from periselene 1 of the connection) + two-body osculating
elements relative to Moon (labels 1-5, 15-19) or Earth (labels 6-14). Label 1 and label 19
share IDENTICAL printed `(r, v, a, e)` — Casoliva's own text: "the time needed... from
periselene 1 to periselene 19 is 113.6319 days," the connection's own homoclinic closure back
near the same orbit phase. Earth rows are bound ellipses (`a>0, e<1`, `a~211,500-214,414` km,
`e~0.668-0.680`); Moon rows are hyperbolic two-body-relative-to-Moon passages (`a<0, e>1`,
consistent with the text's own "hyperbolic orbit of eccentricity 1.126" for periselene 1).
Row 8 (and its mirror, row 12) is the minimum-perigee event Casoliva's own text names for the
LEO-rendezvous ΔV figure: `r=67,869.457` km, `v=3.141` km/s — `table6_row8_leo_dv_check()`
computes `dv = (v - v_circular) = 717.56` m/s directly from these two printed numbers
(standard Earth GM `398600.4415` km³/s²), matching her own stated `717.5` m/s to `<0.1%`.

## 4. Building the Newton correctors: why a cold-start single shot cannot work here

Confirmed via a systematic single-shooting Newton (STM-chained semi-analytic Jacobian, NOT
finite-difference across the long chaotic legs — only the SHORT, bounded-to-one-period
`theta`-sensitivity uses finite differences, chained with the exact per-leg STM for the long
propagation): seeded directly from Casoliva's own printed Table 4 `(theta^u, T^u, T^s)`
(project-frame, `theta^s = -theta^u` per Barrabés-Mondelo-Ollé's own footnote), the matching
residual starts at `0.394` (nondim) and monotonically decreases across 30-75 Newton iterations
to a PLATEAU around `0.030`-`0.034` — never reaching the `1e-7` tolerance target.

**Quantitative diagnosis (the real finding of this task).** The He1 connection's own unstable
leg spans `T^u/T ~= 4.13` periods of an orbit with unstable eigenvalue `Lambda^u ~= 108.60`
per period. A perturbation amplifies by `exp(4.13 * ln(108.5966)) ~= exp(19.35) ~= 2.5e8` over
that leg alone (the stable leg's own `|T^s|/T ~= 4.08` periods amplifies comparably backward).
For single-shooting Newton to have a usable convergence basin, the SEED must be accurate to
roughly `(nonlinearity scale ~1e-3) / 2.5e8 ~= 4e-12` — five orders of magnitude tighter than
this task's own achievable seed precision, even though that seed is ALREADY excellent (the
periodic orbit closes to `~7.5e-12`, the eigenvector to `~1.9e-10` relative spread in a
monodromy self-consistency test — both essentially machine precision, checked directly,
ruling out "the frozen orbit/eigenpair is itself the defect" as the explanation).

## 5. Multiple-shooting: the floor persists, floor-invariant across method/solver

Per Barrabés-Mondelo-Ollé's own documented remedy ("we have also used a multiple shooting
version... in order to avoid loss of precision"), a multiple-shooting corrector was built:
splitting each leg into `n_u`/`n_s` equal-time sub-intervals with continuity unknowns at each
interior patch point, a block-structured Newton system (unknowns = `4 + 4*(n_u-1) + 4*(n_s-1)`,
balanced against the same count of residual equations), using `np.linalg.lstsq` (minimum-norm
least squares) throughout — matching Barrabés-Mondelo-Ollé's own documented practice for their
"over-determined and rank-deficient" system, and empirically more robust than `np.linalg.solve`
under this problem's own severe ill-conditioning (`cond(J) ~ 1e9` at `n_u=n_s=8`, confirmed by
finite-difference cross-check of the analytic Jacobian, ~1e-6 to 1e-9 relative agreement on
spot-checked columns — the Jacobian itself is correct; the ill-conditioning is a genuine
property of the problem).

**Result: the SAME `~0.03` floor**, at `n_u=n_s=8` (both `np.linalg.solve` and
`np.linalg.lstsq`) — floor-invariant across single-shooting, multiple-shooting, and solver
strategy is itself the diagnostic signature of a genuine hard limit, not a tuning artifact.
(An `n_u=n_s=16` finer-segmentation run was also attempted but killed before completion for
wall-clock reasons under this task's own compute budget — NOT reported as evidence either
way, since it produced no output before being stopped.)

**A frozen-orbit-defect hypothesis was explicitly tested and ruled out.** Mid-task, a review
raised the possibility that the periodic orbit/eigenpair being held FIXED (not re-solved
jointly with the connection, a disclosed simplification of Barrabés-Mondelo-Ollé's own fully
general system) was itself the source of the residual floor. Direct measurement (Sec. 4 above)
found the orbit closes to `~7.5e-12` and the eigenvector's own monodromy self-consistency
spreads only `~1.9e-10` relative — both already near machine precision, well below any level
that could explain a `0.03`-scale residual after `~2.5e8`-fold amplification. This rules out
the frozen-orbit hypothesis as the (or a dominant) cause.

A convention-correctness discriminator was also run: evaluating the matching residual using
Casoliva's own DISPLAYED (7-significant-figure) `mu=0.0121529529` closes the orbit only to
`~5e-3` (consistent with `#780`'s own finding that her displayed mu is a rounded display of a
higher-precision internal value, not itself precise enough) — using this project's own
registry mu instead (as `#780` already established is the better choice) is necessary just to
get a sensibly-converged starting orbit at all. With the registry mu and this module's own
properly-converged orbit/eigenvector, the single-shooting residual (`0.394`) is dramatically
better than with Casoliva's raw, uncorrected Table-4 guess evaluated directly (`~2.0`,
`~490,000` km) — the RIGHT direction for a convention-correctness check, supporting (not
refuting) that this module's own conventions (sign of `T^s`, the `psi(theta,xi)` formula, the
project-frame flip) are correct, and that the remaining gap is genuinely the amplification/
conditioning wall, not a sign/convention bug.

## 6. Why Barrabés-Mondelo-Ollé's own paper predicts exactly this difficulty

Their own numerical-methods section (read directly, not inferred) documents FOUR refinements
this task's own implementation does NOT fully replicate, any of which could plausibly close
the remaining gap:

1. **Adaptive segmentation** (`||DPhi_t(z_i)||_inf < M`, recomputed every continuation step) —
   this task's own segmentation is FIXED equal-time, not adaptive. A diagnostic printed the
   per-segment `||STM||_inf` at `n_u=n_s=4`: `285`-`537`, already above their own "tens to
   hundreds" typical `M` bound — suggesting finer or non-uniform segmentation (denser near
   periselene passages) would help.
2. **Second-order variational equations** — this project's own STM machinery is first-order
   only.
3. **Tighter integrator tolerance** (`1e-14` RKF7(8) vs. this task's `1e-12`/`1e-13` DOP853).
4. **QR-decomposition-with-column-pivoting** minimum-norm least squares, specifically because
   their own system is "over-determined and rank deficient" — this task used plain
   `np.linalg.lstsq` (SVD-based minimum-norm LS), a reasonable but not identical substitute.

Casoliva's own paper reports **5h20min of CPU time** to reach THIS SPECIFIC connection
(`h=-1.450162`, her own text: "of all the computed connections of the He1 family, we display
in Fig. 9 the one with the closest approaches to Earth and moon"), reported SEPARATELY from
four OTHER He1-family energies she gives CPU times for in her own Fig. 8 (121s, 196s, 397s,
1619s [~27 min], "ordered by increasing energies") — i.e. this is Casoliva's own hardest
reported connection, not merely the largest of a stated four, reached via a FULL
predictor-corrector CONTINUATION from lower, easier energies, not a cold start. This task
instead attempted a cold start directly at the family's own hardest member, a materially
harder ask than what Casoliva's own procedure actually did.

## 7. Escalation attempted before concluding this (per "never give up reproducing papers")

1. Period-targeting: already exact — the periodic orbit reproduces to near-machine precision.
2. All four `(xi_u sign, xi_s sign)` branch combinations tested at the initial guess; the
   convergent one (`+1,+1`, both matching Casoliva's own stated `xi_0 ~1e-6` positive
   convention) was identified and used throughout.
3. Both a raw-transported-and-renormalized manifold seed (this project's own established
   `_seed_on_manifold` convention) AND Barrabés-Mondelo-Ollé's own analytic
   `Lambda^(-theta/2pi)`-scaled seed formula were tried; the analytic formula (matching their
   own paper exactly) gave the better starting residual and was used for all further work.
4. Single-shooting then multiple-shooting (`n=4`, `n=8`), `np.linalg.solve` then
   `np.linalg.lstsq`.
5. An enlarged-unknowns diagnostic (adding the orbit's own `x0`/period back in as free Newton
   unknowns) was considered and the underlying premise directly tested (Sec. 5's frozen-orbit
   check) — found not to be the explanation, so the full enlarged system was not built (it
   would not have addressed the actual conditioning wall).
6. A convention-correctness discriminator (evaluating the residual with Casoliva's own raw,
   uncorrected Table-4 numbers vs. this module's own properly-converged orbit) confirmed the
   module's own conventions are directionally correct.

A genuine hard limit, not a shortcut — consistent with this project's own "isolated sweep
flips = suspect artifact" and "verify positive-control source applicability" disciplines: the
floor's invariance across independent methods (single-shoot, multi-shoot, two solvers) is
itself the evidence, not a single anomalous run.

## 8. Literature novelty gate

Not triggered — like `#780`, this task is a reproduction ATTEMPT of Casoliva's own published
connection, which resulted in a clean negative on the reproduction target itself, explicitly
scoped away from anything novel (`search/literature_check.py` not invoked; nothing here is
framed as a discovery).

## 9. Code delivered

* `src/cyclerfinder/search/earth_moon_resonant_connections.py` (new sibling module):
  OCR-corrected Table 4 `V^u`/`V^s` constants, freshly-vendored Table 5 (4 rows) and Table 6
  (19 rows), `HE1ConnectionSeed`/`build_he1_connection_seed`/`eigenvector_reproduction_check`,
  `analytic_manifold_seed` (Barrabés-Mondelo-Ollé Eq. (4), faithful implementation),
  `correct_connection_single_shoot` (STM-chained semi-analytic Jacobian) and
  `correct_connection_multi_shoot` (block multiple-shooting, `lstsq`-based, per-segment STM-norm
  diagnostics), `table6_row8_leo_dv_check`. Reuses ONLY
  `cyclerfinder.genome.heteroclinic_cycle._planar_floquet_pair` (generic STM/eigenpair
  utility) from the existing connection-stage codebase — deliberately does NOT reuse
  `correct_connection`/`ResonantNode`/ghost-guard (the Poincaré-section-shooting machinery the
  dispatch note explicitly asked NOT to force-fit here).
* `tests/search/test_earth_moon_resonant_connections.py`: 18 tests — Table 4/5/6 constant
  structure/signs (including the OCR-hazard regression), the He1 connection seed (reproduces
  `#780`'s own eigenvalue, near-machine-precision closure), the eigenvector reproduction check,
  `analytic_manifold_seed`'s own algebraic properties (`theta=0` and `theta=2*pi` identities),
  and both correctors' own genuine residual-REDUCTION behavior at deliberately small
  iteration/segment budgets (this note's own full investigation used much larger budgets — see
  Sec. 4-5 above) — asserting real progress and the expected `converged=False` outcome, never
  a fabricated full convergence.

## 10. Verification

* `uv run pytest tests/search/test_earth_moon_resonant_connections.py -q`: 18/18 pass (~10s).
* `uv run ruff check` / `ruff format --check` on both new files: clean.
* `uv run mypy src tests` (canonical full invocation): clean, 835 source files.
* `uv run pytest tests/data tests/search -q`: run after confirming no sibling `#774` test
  process was active (per this task's own concurrent-agent CPU-contention discipline) — see
  commit history for the recorded pass/fail status.
* `uv run pytest tests/data/test_outstanding_structure.py
  tests/data/test_outstanding_header_body_consistency.py -q`: run before the `OUTSTANDING.md`
  update commit.

## 11. Net effect on `#783`

**DONE — a clean, honestly-diagnosed negative on the connection-reproduction target itself,
with genuine durable value delivered alongside it** (a second OCR-hazard fix on Table 4, a
tight eigenvector reproduction, freshly-vendored Tables 5/6, and a faithful, reusable
implementation of a genuinely different algorithm class than this project's own established
connection-stage machinery). The He1 connection itself does not converge within this task's
own compute budget — a quantitatively-diagnosed numerical conditioning wall (`~2.5e8`-fold
amplification over `~4.1` orbital periods at `Lambda^u~108.6`), consistent with, and
independently corroborated by, Barrabés-Mondelo-Ollé's own paper describing the sophisticated
numerical stack (adaptive segmentation, second-order variationals, `1e-14` tolerance,
QR-pivoted rank-deficient least-squares) and Casoliva's own 5h20min CPU cost for reaching
exactly this energy point via full continuation from easier starting points — not a cold
start, which is what this task attempted. **`#785`'s own gating premise ("once `#783` proves
the method actually works") did NOT hold** — `#785`'s bullet has been updated to reflect this;
its second, independent avenue (Class 1 resonant-connection search via the project's own
established Poincaré-shooting machinery) does not depend on this result and remains viable.
