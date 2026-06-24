# Digest: Peng, Bai & Xu 2017 — Sun-Mercury ER3BP periodic-orbit continuation

Single-paper digest. Read 15/15 pages of the PDF on 2026-06-25 AET.
Same-system (Sun-Mercury) direct check on the #435 e>0 continuation campaign.

## 1. Header

- **Title (verbatim)**: *Continuation of periodic orbits in the Sun-Mercury
  elliptic restricted three-body problem*
- **Authors (CONFIRMED from title page)**: Hao Peng (a,b), Xiaoli Bai (a,
  corresponding, `xiaoli.bai@rutgers.edu`), Shijie Xu (b). NOTE: the task
  brief said "Peng & Xu (verify authors)" — the actual author list is
  **three** authors, Peng, Bai & Xu, with **Xiaoli Bai** (not Xu) as the
  corresponding author.
  - (a) Department of Mechanical & Aerospace Engineering, Rutgers, The State
    University of New Jersey, NJ, USA
  - (b) School of Astronautics, Beihang University, Beijing, China
- **Venue**: *Communications in Nonlinear Science and Numerical Simulation*
  47 (2017) 1-15
- **DOI**: 10.1016/j.cnsns.2016.11.005
- **Received/Revised/Accepted**: 7 May 2015 / 5 November 2016 /
  8 November 2016. Available online 12 November 2016.
- **Length**: 15 pages including references.
- **Acknowledgment of note**: Hao Peng thanks **Prof. Carles Simó** for
  discussion during a visit to the Polytechnic University of Catalonia in
  Barcelona, and credits Simó with suggesting the paper's topic (p. 15).
  This places the work directly in the Barcelona-school ER3BP lineage —
  relevant given the project's existing Andreu/Rosales-Jorba and Ollé/Pacha
  cross-references.

## 2. What the paper actually is

The paper continues **resonant Halo orbits** from the Sun-Mercury **CRTBP**
into the Sun-Mercury **ERTBP** (= ER3BP) by gradually increasing the
eccentricity `e` from 0. The continued objects are called **Multi-revolution
Elliptic Halo (ME-Halo)** orbits around L1 and L2. The central finding is that
the continuation is **non-trivial**: the eccentricity does **not** vary
monotonically along the family characteristic curve, the curves exhibit folds
(turning points) and even closed loops, and **most of the continuation
attempts FAIL to reach the real Mercury eccentricity** `e ≈ 0.2056`.

This is a same-system, same-fidelity-tier counterpart to the project's #435
work, which continued CR3BP Lyapunov/DRO seeds into e>0 at Sun-Mercury and
reported "6/6 survive, no bifurcation, no novel family." Peng et al.'s result
is materially **richer and partly contradictory** in character (see §8).

## 3. The ER3BP model (pulsating frame)

- **Frame**: synodic *pulsating* frame, x-axis always pointing from `m1` to
  `m2`, rotating non-uniformly. Distance unit = instantaneous `r12`
  (pulsates), so the primaries are fixed on the x-axis in this frame (p. 2
  §2.1).
- **Independent variable**: the true anomaly `f` (not time `t`). The chain
  rule `df/dt = (1 + e cos f)^2 / (1 - e^2)^{3/2}` (Eq. 1, p. 2).
- **EOM** (Eq. 2, p. 2): `x'' - 2 y' = ω_x`, `y'' + 2 x' = ω_y`,
  `z'' = ω_z`, with the pulsating potential (Eqs. 3-4, p. 2):
  `ω(x,y,z,f) = (1 + e cos f)^{-1} Ω̃(x,y,z)`,
  `Ω̃ = ½(x²+y²) + (1-μ)/r1 + μ/r2 + ½ μ(1-μ) - ½ e cos f · z²`.
- **No Jacobi integral**: explicitly stated (p. 2, p. 3) — `ω` depends on `f`,
  so the ER3BP is a non-autonomous, **2π-periodic** system, not energy-
  reducible. Setting `e = 0` recovers the autonomous CRTBP exactly.
- **Mass ratio**: `μ_Sun-Mercury = 1.660 × 10⁻⁷` (p. 3 §3.1). (Cross-check:
  Restrepo-Russell 2018 Table 1 gives `1.6601209e-07` — consistent.)
- **Eccentricity used**: real Mercury `e_Sun-Mercury ≈ 0.2056` (p. 1, p. 4
  marks the vertical lines `e = ±0.2056` on every characteristic-curve figure).
  Continuation is carried over the wider band `e ∈ [-0.5, 0.5]` (analysis
  accurate to `|e| ≤ 0.5`; the M5N2 starter was actually pushable to
  `e ∈ (-0.9, 0.7)` but the monodromy eigenvalues blow up beyond `|e|=0.5`,
  reaching condition numbers ~10⁸, so the reliable window is `|e| ≤ 0.5`,
  p. 6).

## 4. The starter families and the "measurable constraint"

The continued objects are **NOT** Lyapunov orbits and **NOT** DROs. They are
**resonant Halo orbits** — specifically multi-revolution Halo orbits that
satisfy a rational period-commensurability with the 2π-periodic ER3BP forcing.

The key structural requirement (Eq. 5, p. 3), the **"measurable constraint"**:

```
T_C = T_E / M = 2 N π / M ,   M, N ∈ ℕ⁺
```

where `T_E` is the full period of the periodic orbit in the ERTBP and `T_C` is
the period of the original (CRTBP) starter Halo orbit. Because the ER3BP
forcing has period `2π`, a periodic orbit must close after the orbit revolves
`M` times while the primaries revolve `N` times — i.e. the starter Halo period
must be a **rational multiple of 2π**. A standard CRTBP Halo period is usually
LESS than 2π, so a resonant (multi-revolution) Halo is needed to satisfy the
commensurability. This is the "Multi-revolution" in ME-Halo.

**Resonance starters tested** (p. 4, §3.1): three pairs `M5N2`, `M7N3`, `M9N4`
(all have **odd M**), continued around **both L1 and L2**. Their CRTBP periods
are `T_C = 4π/5, 6π/7, 8π/9` respectively. Initial CRTBP Halo states are
tabulated (Table 1, p. 4):

| M | N | Location | x0       | z0       | ẏ0        |
|---|---|----------|----------|----------|-----------|
| 5 | 2 | L1       | 0.997148 | 0.004534 | 0.005621  |
| 5 | 2 | L2       | 1.002834 | 0.004548 | -0.005592 |
| 7 | 3 | L1       | 0.996848 | 0.004315 | 0.005790  |
| 7 | 3 | L2       | 1.003134 | 0.004335 | -0.005765 |
| 9 | 4 | L1       | 0.996664 | 0.004096 | 0.005803  |
| 9 | 4 | L2       | 1.003318 | 0.004125 | -0.005784 |

(These are perpendicular x-z-plane crossing states; `y0 = ẋ0 = ż0 = 0`. North
Halo taken as the starter throughout.)

## 5. Continuation method + fold handling (the crux for #436/#437)

- **Method**: a **standard tangential (pseudo-arclength) continuation** in the
  eccentricity parameter `e`, combined with a **multi-segment optimization
  (multiple-shooting) corrector** (p. 3 §3.1, references Allgower-Georg 1990
  for the standard continuation machinery, citation [20]/[21]). The orbit is
  divided into `n` segments; the unknown vector is
  `Γ₀ = [γ₀ᵀ, q₁ᵀ, q₂ᵀ, …, q_{n-1}ᵀ]ᵀ ∈ ℝ^{6n-3}` where `qᵢ = [xᵢ,yᵢ,zᵢ,ẋᵢ,ẏᵢ,żᵢ]`
  and `γ₀ = [x₀,z₀,ẏ₀]ᵀ` is the three nonzero perpendicular-crossing
  components. Continuation parameter is `e`.
- **Corrector engine**: MATLAB `fmincon` with a nonlinear constraint tolerance
  of 1e-11 and an optimization tolerance of 2e-11; integrator `ode45` at
  1e-11 rel/abs tolerance (p. 4). Continuation step auto-alternated between
  1e-5 and 2.5e-3.
- **FOLD / TURNING-POINT HANDLING — directly relevant to #436/#437**: this is
  the paper's most important methodological lesson. The authors report
  (p. 6 §3.2) that **directly continuing the Halo starter in `e` toward a
  target eccentricity FAILS** at turning points: *"when there are two
  corresponding ME-Halo orbits for a given value of e, our routine does not
  have any particular tricks to control the convergence. So it converges to
  two ME-Halo orbits randomly with different initial guess… we have tested to
  directly continue the Halo starter to eccentricities beyond the turning
  points, but the routine failed."* The **tangential (pseudo-arclength)
  continuation is what makes the fold traversable** — it parametrizes by
  arclength along the characteristic curve rather than by `e`, so it folds back
  through turning points smoothly (the "cusps" in Figs. 3-4 are **not** real
  sharp turns; the curve folds back smoothly, p. 6). This is precisely the
  failure mode the project's #436 secant + 2-var symmetric corrector hit:
  a parameter-stepped (in e) corrector cannot pass a fold; an arclength-
  stepped one can. **Peng et al. confirm that arclength continuation is the
  correct tool and that naive parameter-continuation in `e` fails at folds.**

## 6. Which families continue, how far in e, and do they DIE / BIFURCATE

This is the heart of the same-system comparison with #435. Findings, with
page refs:

- **All three starters (M5N2, M7N3, M9N4) successfully continue OFF e=0** into
  nonzero eccentricity, around both L1 and L2 (p. 4 §3.2). So all families
  **do** have a circular (`e=0`) limit — they are continued FROM the CRTBP.
- **The characteristic curves are NON-MONOTONIC in e and have folds/loops**
  (p. 5 §3.3, Figs. 2-5). For M5N2 there are **three separated characteristic
  curves A, B, C** (Fig. 2, p. 4-5):
  - **Curve A**: open curve spanning both positive and negative e; crosses
    `e=0` five times; crosses `e=±0.2056` twice. This is the curve carrying
    the orbits that reach real Mercury eccentricity.
  - **Curve B**: a **closed LOOP**, entirely at positive e, that **never
    reaches `e_Sun-Mercury = 0.2056`** (p. 5, p. 8 §3.3.2). It has four local
    extrema (two maxima, two minima) in e. It is an isolated loop —
    topologically disconnected from curve A.
  - **Curve C**: an open curve entirely at **negative e**, also does not reach
    +0.2056 (p. 6, p. 10 §3.3.2 / Fig. 12).
  - For M7N3 and M9N4 there is by contrast **only one continuous
    characteristic curve** each (p. 6). M9N4's curve folds twice near e=0
    around L1 but not around L2 (asymmetry between L1 and L2).
- **Do families DIE or fail to reach the target e?** YES. The decisive
  same-system finding: **NOT all of the continued family reaches the real
  Mercury eccentricity.** Quote (p. 14): *"one should be careful to check
  whether these ME-Halo orbits can reach the Sun-Mercury eccentricity during
  continuation, since some parts of the curve A and all the curve B in Fig. 2
  cannot."* The M5N2 family around each libration point yields **exactly four**
  ME-Halo orbits at the exact `e = 0.2056`: two on curve A (one at e>0, one at
  e<0 of the ±0.2056 lines), and two on curve C at e<0 (p. 13 §3.5). So
  **eight** real Sun-Mercury ME-Halo orbits total (four per libration point).
- **Do they BIFURCATE?** The characteristic curves change **stability
  structure** along the way (complex-instability onsets), and the **separation
  into three disconnected curves A/B/C is itself a bifurcation-of-families
  signature** — the Periapsis Group vs Apoapsis Group classification "is not
  enough to distinguish them" because for certain eccentricities there are two
  or more ME-Halo orbits, and "some discrete periodic orbits in the CRTBP are
  connected through a smooth variation along e" (Conclusions, p. 14). So the
  paper finds **fold-and-loop topology and stability transitions**, distinct
  from the project's #435 clean "all survive smoothly, no bifurcation."

### 6.1 Stability

- **All M5N2 ME-Halo orbits are HIGHLY UNSTABLE** and exhibit **complex
  instability** for certain parameter ranges (Abstract; p. 8-12). The
  monodromy has three reciprocal eigenvalue pairs `(λᵢ, 1/λᵢ)`; stability
  indices `kᵢ = λᵢ + 1/λᵢ`. Complex instability is flagged when two pairs go
  conjugate off the real axis (Eq. 8, p. 9). Figs. 9/11/13 (L1) and Fig. 16
  (L2) plot `k_{2,3}` vs continuation step.
- The largest eigenvalue can reach `~10⁸` modulus at `|e|=0.5` (p. 6),
  underscoring these are not gentle orbits — they have at least one close
  Mercury encounter per period that amplifies deviations (p. 13-14).

## 7. Tabulated ICs + eigenvalues — usable as a same-model golden

This is the most directly reusable artefact for an ER3BP continuator golden.

**Table 2 (p. 14)** — initial states `γ₀ = [x₀, z₀, ẏ₀]` of the eight exact
Sun-Mercury ME-Halo orbits with M5N2 (`y0 = ẋ0 = ż0 = 0`; `f₀` is the starting
true anomaly):

| Location | Orbit | f₀ | x0       | z0         | ẏ0          |
|----------|-------|----|----------|------------|-------------|
| L1 | 1) | 0 | 0.997862 | 0.00587395 | 0.00198001 |
| L1 | 2) | π | 0.996917 | 0.00389863 | 0.00635667 |
| L1 | 3) | π | 0.99671  | 0.00359628 | 0.00724124 |
| L1 | 4) | π | 0.996201 | 0.00296523 | 0.00870505 |
| L2 | 1) | 0 | 1.00211  | 0.00588612 | -0.00195444 |
| L2 | 2) | π | 1.00329  | 0.00360261 | -0.00722798 |
| L2 | 3) | π | 1.00308  | 0.00390416 | -0.00634580 |
| L2 | 4) | π | 1.0038   | 0.00297788 | -0.00868286 |

**Table 3 (p. 14)** — monodromy eigenvalues `λ_{1,2,3}` of those eight orbits:

| Location | Orbit | λ1        | λ2              | λ3                  |
|----------|-------|-----------|-----------------|---------------------|
| L1 | 1) | 1.092E06   | 5.0771 + 2.593i | 0.15622 - 0.079782i |
| L1 | 2) | -5.7936E05 | 3.0747 + 1.1612i| 0.28464 + 0.1075i   |
| L1 | 3) | 74,342     | 0.89813 + 0.43974i | 0.26376 + 0.96459i |
| L1 | 4) | 88,537     | 1.5723          | 1.2991              |
| L2 | 1) | 1.042E06   | 4.9962 + 2.6003i| 0.15749 + 0.081965i |
| L2 | 2) | -5.1771E05 | 3.0947 + 1.1471i| 0.2841 + 0.10531i   |
| L2 | 3) | 43,176     | -2.2681         | 0.86481 + 0.5021i   |
| L2 | 4) | 65,906     | 1.2733          | 0.9648 + 0.26299i   |

These (μ, e, f₀, IC, eigenvalues) constitute a **complete, sourced, same-model
golden** for an ER3BP multiple-shooting continuator: seed from the M5N2 CRTBP
Halo (Table 1), continue in e along curve A/C with arclength continuation, and
the converged ME-Halo at e=0.2056 must reproduce Table 2 ICs and Table 3
eigenvalues. There is no Jacobi-equivalent (the ER3BP has no Jacobi integral);
the eigenvalue triple is the stability invariant to match instead.

## 8. FOR OUR USE — verdict on #435 (corroborate vs challenge)

**Verdict: Peng et al. partially CHALLENGE the #435 result, and they
CORROBORATE the existence-of-circular-limit half of the #432/#436 question.**

Breaking it into the two sub-claims:

### 8.1 On #435 ("Sun-Mercury families survive into e, no bifurcation, no
novel family")

- **CORROBORATES** the bare survival claim: the Sun-Mercury CRTBP families
  Peng continues **do** survive into `e > 0` off the circular limit — they all
  start at `e=0` and reach nonzero e. None vanishes at the origin.
- **CHALLENGES** the "no bifurcation / smooth, featureless survival" framing.
  In the **same Sun-Mercury system**, Peng finds the characteristic curves are
  **non-monotonic in e, develop folds and turning points, split into multiple
  disconnected curves (A/B/C for M5N2), close into isolated loops (curve B),
  and undergo complex-instability transitions.** Crucially, **parts of the
  family CANNOT be continued to the real Mercury e=0.2056** — curve B never
  reaches it, and parts of curve A do not either. So "all survive smoothly to
  the physical eccentricity" is **NOT** generically true at Sun-Mercury; it
  depends entirely on which family and which branch.
- **The likely reconciliation**: #435 continued **single-revolution** Lyapunov
  and DRO seeds. Peng continues **multi-revolution resonant Halo** seeds
  (M5N2/M7N3/M9N4), which is a different and richer family class — exactly the
  class where the 2π-resonance constraint (Eq. 5) forces multi-revolution
  topology and produces the fold/loop structure. **The #435 "no novel family"
  conclusion may simply be a consequence of the seeds chosen (low-rev Lyapunov/
  DRO), not a property of the Sun-Mercury ER3BP.** The frontier #435 did not
  probe is the multi-revolution resonant-Halo class, where Peng demonstrably
  finds fold-rich, multi-branch, partly-non-reaching behaviour.

### 8.2 On #432/#436 (Antoniadou & Libert 2018 challenge: isolated high-e
elliptic families with no circular limit)

- **Peng does NOT corroborate the existence of an elliptic-ONLY family with no
  circular limit.** Every ME-Halo family in this paper is continued FROM the
  CRTBP (`e=0`); all of curves A/B/C trace back to a circular Halo starter
  (curve B is a loop, but it is a loop reached by continuation from the e=0
  starter, not an isolated e>0-only island). The paper explicitly frames its
  whole method as "resonant Halo orbits are used as continuation starters…
  with variations of the eccentricity" (Abstract; p. 3). So **Peng is silent
  on / does not support** the Antoniadou-Libert isolated-elliptic-family claim
  — it is a different phenomenon. Peng's curve B is the closest analogue (an
  isolated loop in the (ξ, e) characteristic space), but it still attaches to
  the circular family through the continuation curve A↔B↔C structure of the
  M5N2 family.
- **However**, Peng DOES corroborate the *mechanism* that makes isolated
  elliptic families plausible: he shows the same-system characteristic curves
  fold and split into disconnected branches, and he cites Ollé & Pacha (1999,
  ref [22]) — "some discrete periodic orbits in the CRTBP are connected through
  a smooth variation along e" and "similar facts were also observed by Ollé and
  Pacha between some other restricted three-body problems in limiting cases"
  (p. 5). Ollé-Pacha is the canonical "periodic orbits that bifurcate from
  limiting restricted problems" reference — i.e. orbits with no naive circular
  limit. So Peng points AT the Ollé-Pacha mechanism without himself exhibiting
  an elliptic-only family.

### 8.3 Implications for the #435 verdict and the #437 fold-aware follow-on

1. **#435 verdict should be NARROWED, not overturned.** The honest scope of
   #435 is: *low-revolution Lyapunov/DRO seeds at Sun-Mercury survive into e
   smoothly with no fold/bifurcation up to e=0.2056.* That is consistent with
   Peng (he doesn't continue those seeds). #435 should NOT be read as "the
   Sun-Mercury ER3BP has no fold-rich / bifurcating / non-reaching families" —
   Peng demonstrates the opposite for the multi-revolution resonant-Halo class
   in the SAME system. Recommend annotating the #435 result with this scope
   limit and citing Peng et al. 2017 as the same-system counterexample for the
   resonant-Halo class.

2. **#437 (fold-aware continuation) is VALIDATED as the right follow-on.**
   Peng's central methodological lesson is exactly #437's premise: **naive
   parameter-continuation in e FAILS at folds; arclength/tangential
   continuation is required.** This is precisely the #436 failure (secant +
   2-var symmetric corrector proved unreliable at turning points). Peng's
   pseudo-arclength + multiple-shooting recipe (§5 above) is the proven fix in
   this exact system. Recommend #437 adopt: (a) pseudo-arclength
   parametrization (step along curve, not along e); (b) a multi-segment
   (multiple-shooting) transcription, `Γ₀ ∈ ℝ^{6n-3}`, rather than single
   shooting; (c) explicit fold detection (sign change of the arclength
   tangent's e-component) so the continuator KNOWS it has passed a turning
   point and does not mis-classify the two coexisting branches.

3. **Use Tables 2 & 3 as the #437 acceptance golden.** The eight M5N2 exact
   Sun-Mercury ME-Halo ICs (Table 2) + eigenvalues (Table 3) are a same-model,
   same-μ, same-e sourced target. If the #437 continuator reproduces those
   eight orbits from the Table 1 starters via arclength continuation along
   curves A and C, it is validated. (Golden caveat per project rule: the
   EXPECTED side here traces to a published source, Peng Tables 2-3 — not a
   value our own code computed — so it is a legitimate golden.)

4. **Probe the multi-revolution resonant-Halo class at Sun-Mercury.** The
   actionable scientific gap #435 left open is the M5N2/M7N3/M9N4 resonant-Halo
   class. That class — not low-rev Lyapunov/DRO — is where fold/loop/multi-
   branch structure lives in this system, and (per the project's
   `feedback_speculative_high_effort_required` memory) the multi-revolution
   ER3BP regime is exactly a legitimate frontier, not an auto-reject.

## 9. KNOWN_CORPUS / literature-novelty impact

Recommend (do NOT implement — controller consolidates) a `CorpusAnchor` for
`literature_check.py` covering: Sun-Mercury (`μ=1.660e-7`), ER3BP / pulsating
frame, **Multi-revolution Elliptic Halo (ME-Halo)** orbits around L1/L2,
resonances M5N2/M7N3/M9N4, `e ∈ [-0.5, 0.5]`. A discovery candidate that is a
Sun-Mercury resonant elliptic Halo around L1/L2 at M:N ∈ {5:2, 7:3, 9:4} should
return `status="published"` citing Peng, Bai & Xu 2017, DOI
10.1016/j.cnsns.2016.11.005. This is the FIRST same-system ER3BP periodic-orbit
anchor in the corpus and directly hardens the novelty gate for #435/#437
output. Also worth digesting (referenced here, may not be in corpus): Campagnola
et al. (refs [13]/[15]) on elliptic Halo orbits in the ER3BP, and Ollé & Pacha
1999 (ref [22], Astron. Astrophys. 351:1149) on orbits bifurcating from limiting
restricted problems — the latter is the direct #432/#436 Antoniadou-Libert
mechanism reference.

## 10. Errata / vocabulary notes

- **Author attribution**: brief said "Peng & Xu" — correct attribution is
  **Peng, Bai & Xu** with **Bai** corresponding. Flag for the controller.
- **Negative eccentricity**: Peng defines an ER3BP with **negative e** as a
  formal device — an orbit in the Periapsis Group starting at `f₀=0` with
  `-e<0` is equivalent to the Apoapsis Group orbit starting at `f₀=π` with
  `+e>0` (Eq. 6, p. 4: `(-e)cos f = e cos[f + (2k+1)π]`). So the blue (e<0)
  branches in Figs. 2-5 are physically the Apoapsis-Group / phase-shifted-epoch
  representation, not a literally negative eccentricity. Any continuator
  reproducing the curves must implement this `f₀` ↔ `-e` equivalence or it will
  see "half" the characteristic curve only.
- **No methodological contradiction** with the project's CR3BP stack; the e=0
  limit of Peng's EOM is the standard autonomous CRTBP the project uses.

---

## Reference appendix: Peng et al. bibliography (selected, p. 15)

- [6] Broucke, R. *Stability of periodic orbits in the elliptic, restricted
  three-body problem*, AIAA J. 7(6):1003-9, 1969 — complex-instability origin.
- [7] Katsiaris, G. *The three-dimensional elliptic problem*, NATO ASI 1972 —
  3D ER3BP continuation precursor.
- [10] Sarris, E. *Families of symmetric-periodic orbits in the elliptic
  three-dimensional RTBP*, Astrophys. Space Sci. 162:107-22, 1989.
- [13] Campagnola, S. et al. *New techniques in astrodynamics for moon systems
  exploration*, USC 2010 — elliptic Halo / BepiColombo capture context.
- [15] Campagnola, Lo, Newton P. *Subregions of motion and elliptic Halo orbits
  in the ER3BP*, AAS/AIAA 2008 — direct elliptic-Halo-in-ER3BP reference.
- [19] Peng & Xu. *Stability of two groups of multi-revolution elliptic Halo
  orbits in the ER3BP*, Celest. Mech. Dyn. Astron. 123:279-303, 2015,
  DOI 10.1007/s10569-015-9635-2 — the authors' OWN prior ME-Halo paper; the
  "continuation failures in our early studies" this paper resolves. WORTH
  ACQUIRING for the project (same system, prior stability tables).
- [21] Allgower & Georg. *Introduction to Numerical Continuation Methods*,
  SIAM 2003 — the tangential-continuation textbook.
- [22] Ollé, M. & Pacha, J.R. *The 3D elliptic restricted three-body problem:
  periodic orbits which bifurcate from limiting restricted problems*, Astron.
  Astrophys. 351:1149-64, 1999 — the orbits-with-no-naive-circular-limit
  reference; direct #432/#436 Antoniadou-Libert mechanism citation.
- [23] Gómez, G. & Mondelo, J.M. *The dynamics around the collinear equilibrium
  points of the RTBP*, Phys. D 157:283-321, 2001.

End of digest.
