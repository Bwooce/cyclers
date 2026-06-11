# Triage: ML-surrogate trio (Ozaki 2022 / Leifsson 2022 / Wu 2024)

**Date:** 2026-06-11
**Scope:** GENOME ROADMAP triage (search-method improvements), not catalogue mining.
Three papers, one combined note. Depth: triage — full transcription only where a
published trajectory solution warranted it (Ozaki Tables 7–8).

**Verdicts at a glance:**

| Paper | Verdict |
|---|---|
| Ozaki et al. 2022 (JGCD, DNN cycler surrogate) | **USEFUL (method blueprint, deferred)** — direct template for a future combinatorial free-return-chain lane; below breakeven at our current scale. Two example cycler itineraries transcribed. |
| Leifsson et al. 2022 (ICCS, NN-uncertainty adaptive sampling) | **BACKGROUND** — the acquisition-loop *pattern* is worth keeping; the specific two-NN method is dominated by GPR/EGO at our budgets. |
| Wu, Sicard & Gadsden 2024 (ESWA, PIML review) | **BACKGROUND-ONLY** (as expected) — taxonomy extracted; no trajectory/astrodynamics-relevant citations found. |

---

## 1. Ozaki et al. 2022 — Asteroid Flyby Cycler Trajectory Design Using Deep Neural Networks

**Source:** N. Ozaki, K. Yanagida, T. Chikazawa, N. Pushparaj, N. Takeishi,
R. Hyodo, "Asteroid Flyby Cycler Trajectory Design Using Deep Neural Networks,"
*Journal of Guidance, Control, and Dynamics*, Vol. 45, No. 8 (Aug 2022),
pp. 1496–1511, DOI 10.2514/1.G006487 (read as arXiv:2111.11858v3; section/table
numbers per the arXiv version).

### 1.1 Problem formulation (Sec. II)

- **Bodies/model:** Sun-only two-body dynamics; Earth on a *circular ecliptic
  orbit*, zero-radius-SOI patched conics (Sec. II.A, Eq. 1). Asteroid positions
  from real ephemerides (JPL SBDB) inside the optimization; asteroid gravity
  ignored. Transcription: **MGA-1DSM** (one DSM per leg), refined by direct
  multiple shooting (Sec. II.C.1, III.C: 9 variables/phase, 12 equality
  constraints Eqs. 6–11, bounds Eqs. 12–17, SNOPT).
- **Cycler structure:** Earth free-return trajectories (Russell & Ocampo
  classification: full-revolution / half-revolution / generic, parameterized by
  `(m, n, type, v∞)` — Sec. II.B, Table 1, v-infinity globe Fig. 3) chained by
  Earth gravity assists; an asteroid flyby is inserted into each Earth–Earth
  leg. This is exactly the decomposition our free-return-chain machinery uses
  for Earth–Mars work, applied to Earth–asteroid targets. Note their explicit
  contrast (Sec. I): asteroid flyby cyclers do *not* require periodicity of the
  body geometry because each leg targets a different asteroid.
- **Decomposition (Sec. II.C):** the global problem is split into an
  **Earth-Asteroid-Earth (EAE) block** (two consecutive MGA-1DSM legs; inputs
  `t0, v∞,0`, free-return parameters, asteroid elements; outputs `tf, v∞,f,
  ΔV` — Figs. 4–5) and a combinatorial outer search over flyby sequences
  (greedy/beam tree search, Sec. II.C.2, III.F).

### 1.2 What the surrogate predicts and how it prunes (Sec. III)

- **Predicted quantity:** the *outcome of the EAE-block trajectory
  optimization* — not states. Outputs are differences from the underlying
  free-return reference, `Δtf = tf − tf,FR`, `Δv∞,f = v∞,f − v∞,f,FR`, and
  `ΔV`, each squashed by `tan⁻¹(Δz/χ)` with user 90-percentile scales
  (Δtf,90% = 30 d, Δv∞,f,90% = 2 km/s, ΔV90% = 1 km/s) — Eqs. 25–27. Learning
  the *residual relative to an analytic astrodynamics object* (free-return) is
  the load-bearing trick; raw Lambert ΔV alone correlates weakly with the true
  optimal ΔV (Fig. 17b) while the DNN gets ~0.1 km/s (Fig. 17a).
- **Inputs (Tables 2–3, Appendix A):** free-return info `(m, n, type, v∞,0)`;
  asteroid elements `(a, e, i, ΔΩ_t0, ω, M_t0)`; screening features — either
  Lambert ΔVs `(Δv0, Δv1, Δvtotal)` (Eqs. A.2–A.4) or closest-approach state
  difference `(δr_CA, δv_CA)` (Eqs. A.5–A.7); plus initial-guess quantities
  (T, n★, η_t★, v∞ vectors in the RSW frame).
- **Epoch independence (Sec. III.C, Eq. 18, Fig. 10):** post-process the
  database with `Ω → ΔΩ_t0 = Ω − λ⊕(t0)` and `M → M_t0` so the surrogate is
  reusable for any Earth departure epoch without retraining. Cheap and directly
  applicable to any synodic-symmetric problem.
- **Database generation (Sec. III.B–C):** screen candidates with Lambert
  (3-day flyby-epoch grid, keep ΔV < 3 km/s) or closest-approach distance
  (< 5×10⁶ km, Eqs. 2–4), then NLP-optimize the survivors.
- **KKT pseudo-asteroid amplification (Sec. III.D, Eqs. 19–24) — the
  standout idea.** Every converged optimal trajectory spawns ~10 fictitious
  asteroids that the *same* trajectory optimally flies by: position = spacecraft
  position at t1 (Eq. 19); velocity perturbed only perpendicular to the flyby
  Lagrange-multiplier λ6:8 (Eqs. 21–23), preserving all KKT conditions, with
  α_rand bounded so the pseudo-asteroid SMA stays in range (Eq. 24). One NLP
  solve → 11 labelled samples; data-generation rate rises 8.77 → 47.6
  samples/s (Table 4). This is solver-agnostic algebra on stored multipliers.
- **Architecture/training (Sec. III.E, IV.C, App. C):** plain fully-connected
  feed-forward, 5 layers × 1024 units, ELU hidden + sigmoid output,
  batch-norm (large DBs) or dropout (small), Adam lr 1e-4, minibatch 1024,
  MSE loss summed over the three outputs (Eqs. 28–31). Key empirical finding
  (Table 5): **≥ 7×10⁶ samples are needed** before validation loss drops two
  orders of magnitude (6.04e-4 vs ~2e-2); best case 7 (11.7M samples,
  pseudo-asteroids) reaches 3.28e-4 @ 4.5k epochs. Hyperparameter sensitivity
  Table 9: deeper helps (12×1024 → 2.46e-5) at lower throughput. Training ran
  on an NVIDIA Quadro GV100; data generation on an 18-core i9.
- **Pruning mechanism (Sec. III.F, IV.D):** beam search over
  `(free-return params × 15,340 SBDB asteroids)` (q ≤ 1.4 au, Q ≥ 0.8 au,
  OCC ≤ 6); ~300k children/parent, ~10k survive initial screening; the DNN
  evaluates 10,000 EAE costs in **10 s vs 1140 s** for the NLP (~114×). Whole
  search (beam width 100, depth 5): **~10 h with surrogate vs ~7 days
  without**, including database pre/post-processing. Final sequences are then
  re-optimized end-to-end by patched multiple shooting (Sec. IV.E) with ±0.2
  km/s / ±7 d tolerances; DNN single-block ΔV error ≤ ~0.1 km/s, but errors
  accumulate over the chain (their ID1-0704 noted as the large-error case).

### 1.3 Published cycler solutions (Sec. IV.F) — transcribed in full

DESTINY⁺ extended-mission scenario: Earth return after the Phaethon flyby at
**2028 MAY 05 12:13:59 TDB, v∞,0 = 2.684 km/s** (Sec. IV.A). Model caveats for
any catalogue use: patched-conic Sun-only, *circular-ecliptic Earth*, real
asteroid ephemerides; epochs TDB; no independent reproduction here (V0-grade,
publication-sourced).

**Table 7 — ID 2-8373 (nearly ballistic, long transfer; total ΔV = 65 m/s; two
B-type targets):**

| Date time, TDB | Event | v∞ (or v_rel), km/s | Perigee alt, km | ΔV, km/s |
|---|---|---|---|---|
| 2028 MAY 12 00:52:51 | Earth flyby | 2.567 | 232385 | |
| 2028 DEC 27 10:50:45 | 2017 YV8 flyby | 10.063 | | |
| 2029 MAY 12 07:05:22 | Earth flyby | 2.567 | 232386 | |
| 2030 FEB 10 00:45:50 | 2021 BA flyby | 9.453 | | |
| 2030 JUN 13 16:53:49 | Deep space maneuver #1 | | | 0.0285 |
| 2031 MAY 12 14:49:00 | Earth flyby | 2.543 | 40380 | |
| 2031 DEC 26 19:28:45 | 1989 UQ flyby | 6.501 | | |
| 2033 MAY 12 03:13:28 | Earth flyby | 2.543 | 29185 | |
| 2033 NOV 12 10:58:02 | 2017 UX5 flyby | 13.859 | | |
| 2034 MAY 17 14:53:20 | Deep space maneuver #2 | | | 0.0219 |
| 2035 MAY 12 15:40:17 | Earth flyby | 2.545 | 37295 | |
| 2036 JAN 30 23:48:10 | Deep space maneuver #3 | | | 0.0141 |
| 2036 DEC 01 07:20:19 | 1988 XB flyby | 11.423 | | |
| 2037 OCT 31 20:46:55 | Earth flyby | 2.545 | n/a | |

**Table 8 — ID 1-1798 (short transfer; total ΔV = 124 m/s; includes contact
binary 2000 WO107):**

| Date time, TDB | Event | v∞ (or v_rel), km/s | Perigee alt, km | ΔV, km/s |
|---|---|---|---|---|
| 2028 MAY 04 12:08:25 | Earth flyby | 2.691 | 11984 | |
| 2028 NOV 19 01:59:41 | Deep space maneuver #1 | | | 0.0228 |
| 2028 NOV 23 07:49:22 | 1998 XX2 flyby | 8.545 | | |
| 2029 MAY 13 22:38:40 | Deep space maneuver #2 | | | 0.0010 |
| 2029 NOV 03 09:26:10 | Earth flyby | 2.689 | 24795 | |
| 2030 MAY 17 15:23:49 | 2003 LN6 flyby | 4.180 | | |
| 2030 NOV 03 15:38:28 | Earth flyby | 2.689 | 297689 | |
| 2031 JUL 26 22:13:06 | 2016 JJ17 flyby | 8.369 | | |
| 2031 NOV 03 21:50:36 | Earth flyby | 2.689 | 1180355 | |
| 2032 FEB 05 07:41:40 | Deep space maneuver #3 | | | 0.0169 |
| 2032 JUL 29 09:25:09 | 2005 QP11 flyby | 4.054 | | |
| 2032 NOV 03 21:37:26 | Earth flyby | 2.717 | 245400 | |
| 2033 AUG 11 07:23:50 | Deep space maneuver #4 | | | 0.0835 |
| 2033 DEC 13 18:52:26 | 2000 WO107 flyby | 28.477 | | |
| 2034 APR 13 04:44:34 | Earth flyby | 2.779 | n/a | |

Not transcribed: Table 6 (DNN-vs-NLP total-ΔV comparison across example
sequences, Sec. IV.E) — referenced in the text but not rendered in the copy
read; re-pull from the published JGCD version if those comparison numbers are
ever needed.

**Catalogue scope call (flag, no writeback):** these are Earth–asteroid
free-return cyclers, not Earth–Mars. If the catalogue ever takes an
asteroid-cycler class, both rows are V0-ready (event tables above are
per-member reproducible state, unusually good for V0); v4.2 backfill fields
would be center=Sun, mixed source_ephemeris (circular-Earth patched conic +
JPL SBDB asteroids), ToF bounds from the tables. Deferred — out of current
Earth–Mars scope.

### 1.4 Honest applicability to our MBH/scan lanes (16-core CPU, no GPU)

The surrogate's economics: it amortizes ~10⁶–10⁷ NLP solves of *one fixed
parametric block* across a combinatorial outer loop (15,340 asteroids × deep
tree). Our situation differs on every axis:

- **No combinatorial outer loop today.** Earth–Mars cyclers have two bodies and
  low-dimensional genomes; MBH/scan lanes run 10³–10⁵ corrector solves, far
  below the ≥7×10⁶-sample training floor they measured (Table 5). Generating
  that training set IS the search, several times over — net negative.
- **CPU feasibility, for the record:** their base (non-amplified) generation
  ran at 5.75–8.77 samples/s on 18 cores (Table 4) → ~10⁶ solves ≈ 1.5–2 days
  on our 16 cores; with KKT amplification 47.6 samples/s → 7×10⁶ in ~2 days.
  Training a 5×1024 MLP on 7M samples CPU-only is slow (their GPU did ~200k
  samples/s/epoch) but not impossible — order days. So the blocker is economic,
  not technical.
- **When it WOULD pay:** if/when the multi-arc or low-thrust genome ships and a
  Phase 6 re-sweep becomes a *sequence* search over free-return-chain blocks
  `(m, n, type, v∞)` — structurally identical to their beam search over EAE
  blocks — this paper is the direct blueprint, including the residual-vs-
  free-return output encoding, the epoch normalization (Eq. 18), and the KKT
  amplification.
- **Minimum viable adoption now (costless):** (1) make MBH/corrector lanes
  *log* every (seed features → converged ΔV / residual / multipliers) tuple so
  a training database accretes for free; (2) the KKT pseudo-body amplification
  is genome-agnostic algebra on stored converged solutions — implementable and
  testable standalone whenever (1) exists; (3) steal the feature design:
  condition any future screen on cheap analytic objects (Lambert ΔVs,
  free-return reference) and predict the residual, never the raw map.

**Hard rule if ever adopted:** the surrogate prunes *search order only*; every
catalogued/golden value still comes from the real solver chain (consistent with
the sourced-goldens discipline — a NN output is never an EXPECTED value).

---

## 2. Leifsson et al. 2022 — Global Surrogate Modeling by Neural Network-Based Model Uncertainty

**Source:** L. Leifsson, J. Nagawkar, L. Barnet, K. Bryden, S. Koziel,
A. Pietrenko-Dabrowska, "Global Surrogate Modeling by Neural Network-Based
Model Uncertainty," *ICCS 2022* (LNCS), DOI 10.1007/978-3-031-08757-8_35.
7 pages + refs.

### 2.1 The adaptive-sampling loop (Sec. 2, Algorithm 1)

Two feed-forward NNs: a predictor `ŷ_f` fit to `(X,Y)_f`, and an uncertainty
model `ŝ²` fit to the *squared spatial error* `s²(x) = (ŷ_f(x) − y_f(x))²`
(Eq. 2) evaluated on a **separate** LHS-sampled labelled set `(X,Y)_u`. Each
cycle: refit predictor → compute errors on the uncertainty set → fit `ŝ²` →
next infill point `P = argmax ŝ²` (maximized with differential evolution) →
evaluate truth at P, append to `(X,Y)_f` → repeat until `max ŝ²` < tol or max
cycles. tanh activations, lr 1e-3, 3000 epochs, TensorFlow.

### 2.2 Results — toy benchmarks only (Sec. 3)

- **Forrester 1D** (Eq. 3): 3 initial + 10 infill points (10 uncertainty
  samples): RMSE → ~0.1; max model variance 250 → 0.01 (4 orders, Fig. 2).
- **Branin 2D** (Eq. 4): 10 initial LHS + 50 infills (100 uncertainty points):
  RMSE → 0.2; max ŝ² 3.7×10⁵ → 0.1 (6 orders, Fig. 4).

No engineering/domain data; no head-to-head against EGO/Kriging — the authors
themselves defer the state-of-the-art comparison to future work (Sec. 4).

### 2.3 Honest transferability

- **The pattern is right, the implementation isn't (for us).** The
  predict + uncertainty → argmax-infill loop is the generic acquisition-function
  skeleton we'd want for gating expensive corrector/NLP evaluations. But:
  - It is **pure exploration** (max predicted error = global model accuracy),
    not optimum-seeking (no expected-improvement / exploitation term). For a
    search lane we care about minima, not global fidelity.
  - The uncertainty NN needs its **own labelled set** — i.e., extra truth
    evaluations of the expensive function — which is precisely the budget the
    method is supposed to conserve. Against a 1140 s/solve NLP this overhead
    is material; GPR gives the variance for free from the same data.
  - At our dimensionalities (genome dims ~4–8) and budgets (10²–10⁴ samples),
    GPR's cubic cost is irrelevant; their stated motivation for NNs over GPR
    (scaling past ~10⁴ samples / high dimension, Sec. 1) only activates at
    Ozaki-scale databases.
- **Takeaway kept:** if a seed-proposal/acquisition layer is ever added to MBH
  (cheap screen deciding which basins get full corrector treatment), implement
  it as classic GPR/EGO (their refs: Jones et al. 1998 EGO; Forrester & Keane
  2009) first; switch the predictor to a NN only when the logged-solve database
  (Sec. 1.4 item 1 above) passes ~10⁴–10⁵ samples. This paper then describes
  the swap-in.

**Verdict: BACKGROUND.** Method paper validated on toys; useful as a named
pattern and a pointer into the EGO literature, not as an implementation to copy.

---

## 3. Wu, Sicard & Gadsden 2024 — PIML review (condition monitoring)

**Source:** Y. Wu, B. Sicard, S. A. Gadsden, "Physics-informed machine
learning: A comprehensive review on applications in anomaly detection and
condition monitoring," *Expert Systems With Applications* 255 (2024) 124678
(open access). 107 papers surveyed. SKIMMED ONLY (Secs. 1–3.4.1 read; 3.4.2–5
TOC/table level), per triage plan.

### 3.1 The taxonomy (Sec. 3) — kept as a map

Four integration frameworks (their enumeration, Sec. 3 intro):

1. **Physics embedded in feature space** (Sec. 3.1) — physics-in-features:
   physics-guided input feature augmentation, synthetic data from physical
   models / FEM / digital twins (3.1.1, Table 1), and transfer learning from a
   physics-based source domain (3.1.2, Table 2, Figs. 5–6). *This is what Ozaki
   does* (Lambert/free-return features + screened-NLP synthetic database).
2. **Data-enhanced refinement of physical models** (Sec. 3.2, Table 3,
   Fig. 7) — residual/hybrid modelling: ML learns the *discrepancy* between a
   physics model's prediction and truth. The natural pattern for any future
   "physics-informed seeding of correctors": analytic model proposes, ML
   predicts the correction. Noted limitation (their words): the ML learns the
   discrepancy, not the system — weak interpretability.
3. **Physics-informed regularization** (Sec. 3.3, Table 4) — physics-regularized:
   PINN-style loss terms penalizing violation of governing equations
   (3.3.1, Eqs. 6–10: Loss = λ₁L_data + λ₂L_PDE + λ₃L_BC + λ₄L_IC),
   data-driven DE solutions (3.3.2), and hybrid losses bolted onto
   CNN/autoencoder/RNN architectures (3.3.3). Known pain: harder loss
   landscape, generalization issues (Sec. 3.3 closing).
4. **Physics-guided design of architectures** (Sec. 3.4, Table 5) —
   physics-architecture: physically meaningful nodes/connections, physics-based
   activation functions, and the **Hamiltonian Neural Network** class (Shen et
   al. 2023, building on Greydanus et al. 2019) encoding energy conservation
   for dynamical systems — the only architecture here with a conceivable
   astrodynamics use (conservative-dynamics surrogate propagators).

This maps cleanly onto the planned 3-way mental model: physics-in-features =
(1); physics-regularized = (3) (+ (2) as the residual special case);
physics-architecture = (4).

### 3.2 Citations relevant to trajectory/astrodynamics surrogates

**None directly.** All 107 surveyed works are condition-monitoring/anomaly
detection (bearings, fatigue, batteries, structural health, power grids).
Checked Tables 1–5 and Secs. 3.1–3.4.1. Two tangentially worth remembering:

- **Greydanus et al. 2019, Hamiltonian Neural Networks** (via Shen et al.
  2023, Sec. 3.4.1/3.3.3) — energy-conserving learned dynamics; relevant only
  if we ever consider a learned propagator (we shouldn't, for catalogue work).
- **Yang et al. 2019, CoPhIK (physics-informed CoKriging)** (Table 1,
  Sec. 3.1 discussion) — multifidelity GP fusing low-fidelity physics-model
  output with sparse high-fidelity data, with greedy active learning; the
  multifidelity pattern matches our tier ladder (Lambert → patched-conic →
  DE440 → n-body) better than anything NN-based in this review.

**Verdict: BACKGROUND-ONLY**, as expected. Do not re-mine; the taxonomy above
is the extractable content.

---

## 4. Genome-roadmap implications (synthesis)

**Does an NN-surrogate lane earn a roadmap task now? No.** Three independent
reasons converge:

1. **Economics (Ozaki):** the measured training floor is ~7×10⁶ labelled
   solves; surrogate breakeven needs ~10⁶⁺ repeated evaluations of one fixed
   block. Our lanes run 10³–10⁵ solves over two bodies — the training set would
   cost more than the searches it accelerates.
2. **The catalogue is new-input-gated, not evaluation-gated.** The validation
   ceiling is a publication-data limit; no amount of search acceleration moves
   it.
3. **Acquisition-function value at our scale is delivered by GPR/EGO, not NNs**
   (Leifsson's own scaling argument cuts against NNs below ~10⁴ samples).

**What earns a (small, deferred) task instead — phase 1, costless prep:**

- **Log-everything instrumentation:** have MBH/corrector/scan lanes persist
  every (seed features, converged ΔV, residual, multipliers) tuple. Zero new
  compute; builds the training database retroactively for whenever the
  breakeven flips. Without this, any future surrogate lane starts from zero.
- **KKT pseudo-body amplifier (Ozaki Sec. III.D)** as a standalone utility over
  the logged solutions — ~11× database amplification by algebra, genome-
  agnostic, unit-testable against the KKT conditions.
- **Feature-design conventions now:** any screen/score we add should condition
  on cheap analytic references (Lambert ΔV, free-return parameters) and predict
  residuals (Ozaki Eq. 25 pattern), with synodic epoch normalization (Eq. 18
  pattern).

**Trigger to revisit (write into the roadmap, not a task yet):** the day a
shipped multi-arc / low-thrust genome turns Phase 6 re-sweeps into a
combinatorial *sequence* search over free-return-chain blocks (beam/tree over
`(m, n, type, v∞)` × bodies), Ozaki 2022 is the implementation blueprint
end-to-end, with Leifsson→EGO as the infill layer and Wu's framework-1/2
patterns (physics features + residual learning) as the design language. Until
then: surrogates may reorder the search queue, but never decide catalogue
truth — final numbers always come from the real solver chain.
