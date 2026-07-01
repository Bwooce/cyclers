# Sun perturbation does not transfer to Sun-Jupiter-moon systems for the L1 Lyapunov family

**Task:** #326 (write-up of the #313 structural finding)
**Companion data:** `data/scan_313_sun_jupiter_europa.jsonl`, `data/scan_313_sun_jupiter_io.jsonl`, `data/bcr4bp_l1_family_303.jsonl`
**Companion commit:** `2ad6d01` (data + script), `5c2b741` (prior #313 verdict note)
**Status:** structural finding; clean negative on a published, falsifiable hypothesis.

---

## 1. Setup

The **circular restricted three-body problem (CR3BP)** models the motion of a massless
particle under two primaries on circular co-orbits. The **bicircular restricted
four-body problem (BCR4BP)** adds a third gravitating body (the Sun) on a
circular orbit *about the primary–secondary barycentre*, and the test particle
sees its gravity but does not perturb it (incoherent / kinematic). The Sun
contribution is parametrised by four constants — mu (the primary mass ratio),
mu_sun (the Sun mass in primary+secondary = 1 units), a_sun_nondim (the Sun's
semi-major axis in primary–secondary distance units LU), and omega_sun_nondim
(the Sun's angular frequency in the synodic frame).

For the **Sun-Earth-Moon (SEM)** system the canonical values are

| constant | value | source |
|---|---|---|
| mu_EM | 0.012150581600000 | Rosales-Jorba 2023 Table 3 |
| mu_sun | 328900.5423094043 | Rosales-Jorba 2023 Table 3 |
| a_sun_nondim | 388.8111430233511 | Rosales-Jorba 2023 Table 3 |
| omega_sun_nondim | 0.925195985520347 | Rosales-Jorba 2023 Table 3 |

This is the regime studied by Andreu (1998), Simó / Jorba / Gómez (1995),
Gimeno-Jorba (2018), Rosales-Jorba (2023). In that literature the Sun
perturbation is large enough to materially deform CR3BP families and to
create genuinely four-body invariants (the POL1 / POL2 substitute families,
the Sun-commensurate halos).

The natural question this note addresses: **does the same machinery, applied
at Sun-Jupiter-moon mass ratios and length scales, produce the same kind of
qualitative deformation of the CR3BP L1 Lyapunov family?** Jupiter's μ_sun
(Sun-vs-Jupiter system mass ratio) is much smaller (~1048 vs ~329 000), but
the satellite's LU is also much smaller, so a_sun_nondim grows substantially.
Which effect dominates is an empirical question.

## 2. Method

The substrate is the `BCR4BPSystem` of #292 (Phase 1), parameterised by
`(mu, mu_sun, a_sun_nondim, omega_sun_nondim)`. The continuation driver of
#303 (Phase 2) and the corrector of #292 are reused unchanged.

For each system:

1. Build the primary–secondary CR3BP from sourced GM values.
2. Correct a planar L1 Lyapunov at a Jacobi level just below C_L1 via
   `correct_symmetric_free_period` (the routine used in #303 for Earth-Moon).
3. Build a `BCR4BPSystem` with `mu_sun = 0` and the Sun-system constants
   (mu_sun, a_sun_nondim, omega_sun_nondim computed from GM_Sun and the
   primary–secondary geometry).
4. Re-converge the CR3BP seed via `correct_bcr4bp_periodic` at `mu_sun = 0`
   — this is the **CR3BP-limit anchor**, a structural-correctness test:
   if BCR4BP at mu_sun=0 does not reproduce the CR3BP orbit to numerical
   precision, the constants are wrong and we stop.
5. Continue the family in `mu_sun` from 0 to the system's full mu_sun via
   `continue_bcr4bp_family_in_musun` (geometric step schedule, 50 steps).
6. Record the converged IC drift, half-period corrector residual, and an
   independent full-period (Radau) closure residual.

The deliverable is per-step (x0, y0, z0, vx0, vy0, vz0, T_TU, mu_sun_value,
corrector_residual, independent_closure_residual, stability_tag) recorded
to JSONL, plus per-system anchor and continuation summary rows.

### System constants used in #313

Sourced from the existing project registries (`core/constants.py`,
`core/satellites.py`) — no new constants introduced.

* GM_Sun = 1.32712440018e11 km³/s² (IAU 2015 / JPL DE440)
* GM_Jupiter_sys = 1.26686534e8 km³/s² (JPL SSD `gm_de440`)
* Jupiter SMA = 5.20288700 AU = 7.7822e8 km (JPL DE440)
* AU = 1.49597870700e8 km (IAU 2012 Resolution B2 exact)
* Europa: SMA 671100 km, GM 3202.739 km³/s² (JPL SSD)
* Io: SMA 421800 km, GM 5959.916 km³/s² (JPL SSD)

Computed BCR4BP constants (line 1 of each scan JSONL, `row_type=constants`):

| System | mu | a_sun_nondim | mu_sun | omega_sun_nondim |
|---|---|---|---|---|
| Sun-Earth-Moon | 0.01215058 | 388.811143 | 328900.5423 | 0.925196 |
| Sun-Jupiter-Europa | 2.5281e-5 | 1159.798565 | 1047.5655 | 0.999181 |
| Sun-Jupiter-Io | 4.7045e-5 | 1845.284060 | 1047.5655 | 0.999592 |

## 3. Quantitative result

All three continuations terminated with `survival_fraction = 1.0` (50/50
members converged). The IC drift across the full continuation is dramatically
different between SEM and the Jupiter systems:

| Metric | SEM (303) | Sun-Jupiter-Europa | Sun-Jupiter-Io |
|---|---|---|---|
| Target mu_sun | 328900.542 | 1047.565 | 1047.565 |
| CR3BP-limit anchor `corrector_residual` | n/a (built at limit) | 2.54e-14 | 7.86e-14 |
| CR3BP-limit anchor independent closure | n/a | 4.30e-13 | 2.03e-12 |
| seed x0 | 0.81152565 | 0.97681122 | 0.97208813 |
| seed vy0 | 0.25618426 | 0.02395627 | 0.02408637 |
| final x0 (full mu_sun) | 0.81142014 | 0.97681143 | 0.97208824 |
| final vy0 (full mu_sun) | 0.24569193 | 0.02395267 | 0.02408502 |
| **Δx0 endpoint** | **1.055e-4** | **2.078e-7** | **1.082e-7** |
| **Δvy0 endpoint** | **1.049e-2** | **3.594e-6** | **1.348e-6** |
| Final corrector residual | 4.81e-14 | 6.74e-15 | 6.40e-11 |
| Final independent closure | 1.08e-2 | 2.03e-7 | 6.11e-8 |
| Final stability_tag | hyperbolic_pair | hyperbolic_pair | hyperbolic_pair |

**Headline ratio (Δx0):**
- SEM / Sun-Jupiter-Europa ≈ 1.055e-4 / 2.078e-7 ≈ **508×**
- SEM / Sun-Jupiter-Io ≈ 1.055e-4 / 1.082e-7 ≈ **975×**

That is **2.7–3.0 orders of magnitude** smaller Sun-perturbation effect at
Sun-Jupiter-moon than at Sun-Earth-Moon for the L1 Lyapunov family. The Δvy0
ratios tell the same story (∼3000× and ∼7700×). At Jupiter the converged
BCR4BP IC differs from the CR3BP seed by only ~10⁻⁷ in nondim units; this
is far below any practical observability threshold and is dominated by the
corrector's own numerical noise (the anchor independent closure is itself
~10⁻¹³–10⁻¹², so the continuation's ~10⁻⁷ drift is a few orders above noise
but qualitatively negligible).

The independent (Radau) full-period closure residual rises from ~10⁻¹³ at
the CR3BP-limit anchor to ~2×10⁻⁷ at full mu_sun for SJE and ~6×10⁻⁸ for SJI
— consistent with a real but tiny Sun perturbation. For SEM the same metric
rises to 1.08e-2 (per #303 — the independent residual is *reported, not used
as a hard gate* because T is a free variable in the continuation and the Sun
phase is not enforced commensurate; the corrector half-period residual stays
at machine precision throughout).

## 4. Geometric explanation

The Sun's gravitational acceleration on a test particle in the primary–
secondary synodic frame, at heliocentric distance r ≈ a_sun_nondim, is
approximately

  a_Sun_direct ~ μ_sun / a_sun_nondim²

and the **tidal** (differential) acceleration relevant to a perturbation of
a relative orbit of size O(1 LU) is

  a_Sun_tidal ~ μ_sun × LU / a_sun_nondim³ = μ_sun / a_sun_nondim³

In the rotating synodic frame, the mean Sun acceleration at the
primary–secondary barycentre is balanced by centrifugal terms, so the
*differential* (tidal) term is the dominant perturbation on the relative
motion. Computing both for our three systems:

| System | μ_sun | a_sun_nondim | μ_sun/a_sun² (direct) | μ_sun/a_sun³ (tidal) |
|---|---|---|---|---|
| SEM | 3.289e5 | 388.81 | 2.176e0 | 5.596e-3 |
| SJE | 1.048e3 | 1159.80 | 7.788e-4 | 6.715e-7 |
| SJI | 1.048e3 | 1845.28 | 3.076e-4 | 1.667e-7 |

**Ratio SEM / SJE:**
- μ_sun alone: 313.97×
- a_sun⁻²: 1/0.1124 = 8.90×
- a_sun⁻³: 1/0.0377 = 26.52×
- direct accel: 313.97 × 8.90 = **2794×**
- tidal accel: 313.97 × 26.52 = **8333×**

**Ratio SEM / SJI:**
- direct accel: **7072×**
- tidal accel: **33563×**

The observed Δx0 ratios (508× and 975×) sit *between* the linear
μ_sun ratio (314×) and the tidal-acceleration ratio (8333× / 33563×). This
makes sense: at small but finite mu_sun the leading effect on the IC scales
sub-linearly with the per-particle perturbation, because the continuation
adjusts the IC to maintain periodicity rather than letting drift accumulate
linearly. The observation falls in the order-of-magnitude regime predicted
by the geometric scaling — the Sun-Jupiter case is **dramatically weaker per
particle** than the Sun-Earth-Moon case, and that is what we measure.

**The bottom line for back-of-envelope estimation:** if you want to know
whether Sun-perturbation effects in a primary–secondary system will be
"Andreu-class large" or "machine-precision negligible", the right
single-number summary is

  ratio ≈ μ_sun_target / μ_sun_SEM × (a_sun_SEM / a_sun_target)^k

with **k ∈ [2, 3]** depending on whether the orbit is more sensitive to the
direct or the tidal term. For Sun-Jupiter-Europa this gives ratios from
1/3000 to 1/8000 relative to SEM, predicting essentially-negligible effects
on the SEM-scale Sun perturbations.

## 5. Why it matters

This is a falsifiable, sourced **negative** result on the question:
*should BCR4BP analyses developed at Sun-Earth-Moon extend to Sun-Jupiter-moon
systems?* The answer for the L1 Lyapunov family at this Jacobi level is

> **No. Δx0 ≲ 2e-7 nondim across full mu_sun continuation — the BCR4BP
> orbit is indistinguishable from its CR3BP parent at any practical
> tolerance.**

Practical consequences:

1. Future investigations of Jovian-moon dynamics need not build a BCR4BP
   stack expecting Andreu-class deformations of the L1 Lyapunov family —
   the CR3BP is already a quantitatively accurate model for this family
   at this Jacobi level.
2. The intuition "Sun perturbation matters in restricted N-body problems"
   needs to be re-evaluated per system. It is a function of *both* μ_sun
   and a_sun_nondim, with the second entering to power 2–3. The Sun is
   ~1048 in μ_sun units at Sun-Jupiter, but its semi-major axis is
   ~3–5× larger in moon LU, so the perturbation is two-to-four orders of
   magnitude weaker.
3. Conversely, this **strengthens** the case for treating Sun-Earth-Moon
   as a genuinely four-body problem when working at the L1 Lyapunov
   family — the Sun is not numerically negligible there.

## 6. Falsifiable scope — what this note does NOT establish

The claim is precisely scoped:

- **Planar L1 Lyapunov family**, at the specific Jacobi level used (just
  below C_L1) and the specific x0/vy0 root chosen by the corrector. Other
  members of the same family at different amplitudes may behave differently.
- **Sun-Jupiter-Europa** and **Sun-Jupiter-Io** specifically. Other Jovian
  moons (Ganymede a_sun_nondim ≈ 728, Callisto a_sun_nondim ≈ 414) lie
  between SEM and SJE/SJI in a_sun_nondim; their behaviour interpolates
  but should be measured.
- **Standard (incoherent) BCR4BP** per Simó / Jorba / Gómez. The **coherent
  Quasi-BiCircular Problem (QBCP)** of Andreu, with α-table corrections,
  was *honestly deferred* in #292 Phase 1 because the 2026-06-14 digest
  did not source the α tables. QBCP corrections might recover some Sun
  effect at Jupiter — this is an open question.

**The note explicitly does NOT establish anything about:**

- **3D halo families at Sun-Jupiter-moon.** Note that contrary to a
  superficial expectation, #304's halo continuation at *SEM* itself
  shows **Δx0 ≈ 7.5e-4** across the same mu_sun range — even larger than
  the L1 Lyapunov result (1.05e-4). Halos are *more* sensitive to Sun
  perturbation at SEM, not less. We have NOT scouted the Sun-Jupiter-moon
  halo case directly; the geometric a_sun scaling suggests it should still
  be small at Jupiter, but the family-dependent factor remains to be
  measured.
- **Distant retrograde orbits (DROs).** Large-amplitude orbits sample a
  larger Sun-acceleration gradient; the tidal scaling argument may need
  to be re-derived.
- **L4 / L5 triangular libration points.** Different linear stability
  structure; resonance with Sun-commensurate motions could matter.
- **Higher-eccentricity primary–secondary systems** or **eccentric Sun
  motion** (ER4BP, ER3BP-with-Sun). Those introduce time-dependent terms
  not captured by the BCR4BP machinery.
- **Sun-Saturn-moon** (e.g. Sun-Saturn-Titan). a_sun_nondim is yet larger
  than SJI (Saturn is ~9.5 AU vs Jupiter's 5.2 AU) and μ_sun_Saturn-sys is
  smaller still, so this should be an even *weaker* Sun perturbation — but
  not measured.
- **Sun-Mars-moon** (where the Mars-moon LU is tiny — Phobos ~ 9377 km,
  giving a_sun_nondim ≈ 24000+). #313 Part A (Mars moons tulip search)
  did not reach the BCR4BP step due to seeding problems, so we have no
  data here either.

Recommended Phase 2 / 3 work: map the regime boundary in (μ_sun, a_sun_nondim)
space — at what combination does the Sun perturbation switch from
"machine-precision negligible" to "Andreu-class material"? A first cut at
the simple Sun-tidal scaling argument suggests the boundary is roughly
μ_sun / a_sun_nondim³ ≳ 10⁻⁴ — SEM sits at 5.6e-3 (well above), Sun-Jupiter
sits at 6.7e-7 to 1.7e-7 (three orders of magnitude below).

## 7. Caveats and future work

1. **Single-member metric.** Δx0 and Δvy are reported for one corrected IC
   per system. The full family at multiple amplitudes was not mapped here;
   individual members at higher amplitude may have larger drifts (#304's
   halo result hints at this).
2. **Independent-closure ambiguity at full mu_sun.** The corrector enforces
   symmetric half-period closure; the full-period (Radau) check accumulates
   the Sun-phase residual because T is free and Sun commensurability is not
   enforced. The reported "final independent closure" of ~2e-7 (SJE) /
   ~6e-8 (SJI) therefore conflates true family drift with Sun-phase free
   parameter slop. This is the *same* discipline used in #303 — the half-
   period corrector residual is the binding metric, the independent residual
   is reported but not used as a gate. Sourcing: orbit-closure discipline
   memo.
3. **QBCP refinement not run.** The standard BCR4BP is known to differ from
   Andreu's coherent QBCP at O(μ_sun²); whether the QBCP corrections shift
   the Sun-Jupiter conclusion is open.
4. **No 3D continuation at Sun-Jupiter-moon.** The corrector substrate works
   in 3D but the #313 scout used the planar Lyapunov seed. A 3D halo
   continuation at SJE / SJI would close the gap with #304 and is
   recommended.

## 8. References

- **Andreu, M. A.** (1998). *The Quasi-Bicircular Problem.* PhD dissertation,
  Universitat de Barcelona. (In the private paper corpus — not publicly
  cited from this repo per project policy.)
- **Simó, C., Gómez, G., Jorba, À., Masdemont, J.** (1995). The bicircular
  model near the triangular libration points of the RTBP. In *The dynamical
  behaviour of our planetary system* (Roy & Dvorak, eds.), Kluwer.
- **Gimeno, J., Jorba, À.** (2018). On the numerical computation of QBCP
  parameters and substitute periodic orbits. *Frontiers AMS* 4: 32.
  Tables 3–4.
- **Rosales, J., Jorba, À., Jorba-Cuscó, M.** (2023). Families of
  Sun-perturbed Halo and Lyapunov orbits in the Earth-Moon system.
  *Celestial Mechanics and Dynamical Astronomy* 135: 15. **Table 3
  (primary parameter source for the SEM BCR4BP constants used here).**
- **#292 Phase 1** — BCR4BP substrate construction and CR3BP-limit anchor
  test (`docs/notes/2026-06-16-292-bcr4bp-phase1.md`).
- **#303 Phase 2** — mu_sun continuation driver, SEM L1 Lyapunov family
  bridge (`docs/notes/2026-06-16-303-bcr4bp-phase2-musun-continuation.md`,
  `data/bcr4bp_l1_family_303.jsonl`, commit `aa53e39`).
- **#304 Phase 3** — SEM halo family BCR4BP continuation
  (`docs/notes/2026-06-16-304-bcr4bp-halo-phase3.md`,
  `data/bcr4bp_halo_family_304.jsonl`).
- **#313 Phase 1 Part B** — Sun-Jupiter-Europa / Sun-Jupiter-Io scouts
  (`scripts/scan_313_sun_jupiter_moons.py`,
  `data/scan_313_sun_jupiter_europa.jsonl`,
  `data/scan_313_sun_jupiter_io.jsonl`, commit `2ad6d01`).
- **#313 Phase 1 verdict** — multi-system scouts write-up
  (`docs/notes/2026-06-16-313-multi-system-scouts.md`, commit `5c2b741`).

## 9. Data provenance and reproducibility

- All numeric values in this note trace to the JSONLs above and the
  `_build_system_constants` function in `scripts/scan_313_sun_jupiter_moons.py`.
- Re-run with `uv run python scripts/scan_313_sun_jupiter_moons.py`
  (wall-clock 18.3s for both moons on the development laptop, per commit
  `2ad6d01`).
- Sanity check on the BCR4BP machinery: `uv run pytest
  tests/genome/test_bcr4bp_genome.py tests/genome/test_bcr4bp_continuation.py`
  — 15 passed on 2026-06-16.
- No catalogue writeback. This note is a structural-finding deliverable only.
