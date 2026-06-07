# Zhang, Shi & Zheng 2024 — neural angle-only OD for Earth-Moon libration orbits (scope-check mining)

Mined 2026-06-07, user-requested scope check. Triaged as likely out-of-scope
before the read; the verdict below confirms it.

**Source (cite exactly, no file path):**
Zhang, Z., Shi, Y., & Zheng, Z., "Improving Angle-Only Orbit Determination
Accuracy for Earth-Moon Libration Orbits Using a Neural-Network-Based Approach,"
*Remote Sensing* **2024**, 16, 3287. DOI 10.3390/rs16173287. MDPI, open access
(CC BY). Beijing Institute of Technology / Peng Cheng Laboratory.

> Clean digital MDPI typeset; equations, all 12 tables, and figures read
> unambiguously. Full read pp.1-25 (all content + references).

---

## 0. VERDICT — honest scope call (read first)

**OUT OF SCOPE. Nothing maps to the cycler catalogue or to any active design
lane; the single thin methodological touch-point is already covered by our
settled ML-surrogate stance.**

This is an orbit **determination** paper: given a spacecraft already flying a
known periodic orbit, and a stream of noisy ground-based angle-only (RA/DEC)
observations, estimate its state more accurately. It does **not** design,
enumerate, or optimise any trajectory. The "neural network" is a measurement
**outlier-detector** bolted onto a batch least-squares estimator — it never
produces a trajectory, a Δv, a sequence, or a transfer. The dynamics are
Earth-Moon cislunar (CR3BP-seeded, ephemeris-refined), not heliocentric cycler.
There is no overlap with our Earth-Mars patched-conic cycler search, and the
nearest possible consumer (#76 moon-tour) shares only the *words* "Earth-Moon"
and "libration orbit," not any reusable method, constant, or convention.

The one genuine (small) value: the paper is a clean worked example of the
**prune/propose-but-never-testify** rule our ML-surrogate note already enforces —
its NN is explicitly *advisory* (it flags suspect measurements; the classical
LSM still does the estimating), which is exactly the discipline we already
demand. It corroborates our stance; it does not extend it. Details in §4b.

---

## 1. Problem setup (pp.1, 6-11)

- **What is being determined:** the 6-D state `x = [r; v]` of a spacecraft on an
  Earth-Moon **libration point orbit (LPO)**. Training families (Fig.4, p.9;
  Tables 3-4, pp.9-10): **Southern L2 Halo (25 orbits), Northern L2 Halo (25),
  L2 Lyapunov (21)** — 69 distinct periodic orbits after de-duplication (p.11).
  The **test/evaluation** orbit (deliberately held out of training to show
  generalisation) is a **4:1 synodic resonant near-rectilinear halo orbit
  (NRHO)** near Earth-Moon L2 (abstract p.1; §4, p.18; Fig.12, p.19).
- **Observer / geometry (pp.6-7):** a single **ground-based optical station**
  measuring line-of-sight angles only — right ascension α and declination β
  (Eq.14, p.7). Station = **Eglin, FL** (lat 30.57°, lon -86.21°, alt 34.7 m;
  Table 2, p.7). Angle noise σ_RA = σ_DEC = **2 arcsec**, Gaussian. Elevation
  constraint 0-90° (spacecraft invisible when behind Earth) — this gives the
  short, sparse **visible arcs** (Figs.6, 12).
- **Dynamics model (p.6, Eq.13):** *high-fidelity ephemeris* point-mass gravity
  from Earth + Moon + Sun, with Sun/Moon positions from the **VSOP2000**
  planetary theory (refs 64-66). Integrated with **RK4(5)**, rel/abs tol 1e-12.
  Periodic orbits are first generated in the **CRTBP** (μ implicit), then
  corrected into the ephemeris model with a **two-level differential corrector**
  (ref 73, Howell-Pernicka), epoch 1 Jan 2020, propagated ~1 month (27.32 d).
  So: CRTBP for the *catalogue of shapes*, ephemeris for the *actual dynamics*.
- **Why angle-only OD is weak here (observability):** the paper does not give a
  formal observability-rank analysis. The stated weakness is empirical and
  data-driven: angle-only data are **short-arc, sparse, and noisy** (p.1, p.2),
  the dynamics are highly nonlinear (three-body), and **range is unobserved**
  (angles only) so the estimate is fragile to measurement outliers. The whole
  contribution is aimed at the *outlier* failure mode, not at a geometric
  degeneracy/observability-cone argument. (Honest note: the task framing hoped
  for a "geometric analysis of when angle-only data is degenerate." This paper
  does **not** provide one — see §4c.)

---

## 2. The NN method (pp.2-3, 11-18)

**What it learns:** NOT an initial-state estimator, NOT a dynamics surrogate,
NOT an observability weighter. It is a **binary classifier of measurement
quality** — for each incoming angle measurement it predicts "accurate" vs
"inaccurate," where *inaccurate* is *defined* (p.2, p.7) as true measurement
error > 2σ (i.e., > 4 arcsec). Two independent networks: one for RA, one for DEC.

**Architecture (Fig.7, p.13; Table 6, p.14):** a **modified bidirectional LSTM
(BiLSTM)**. Sequence in = per-epoch feature vector `s_k ∈ R^14` (Eq.16, p.8):
relative measurement error `ε̃_k ∈ R^2`, the measurement `z_k ∈ R^2`, station
ECI position `r_S ∈ R^3`, LSM estimated state `x̂ ∈ R^6`, measurement interval
`Δt_k ∈ R^1`. Output per epoch = `y_k ∈ R^4` (one-hot accurate/inaccurate for RA
and DEC; Eq.17, p.8). Wrappers around the BiLSTM cells: an initialization layer,
a pre-BiLSTM DNN, a post-BiLSTM DNN (re-weighting), and a sigmoid layer. Chosen
hyperparameters: **1 hidden layer, 32 neurons, Leaky ReLU** (Table 7, p.15 —
selected by min validation cross-entropy). The bidirectional structure is the
point: it uses *future* as well as past epochs to judge the current measurement,
which a forward-only LSTM cannot (p.2).

**What it does with the classification (the "compensation," Eqs.23-25, pp.17-18):**
once a measurement is flagged inaccurate, it is *nudged* by one STD in the
direction of its relative error sign: `α_k ← α_k − sign(ε̃_k^RA)·σ_RA` (and the
DEC analogue). Magnitude is exactly 1σ — deliberately small, to avoid harming
the (~13%) misclassified good measurements. The corrected measurement stream is
then re-fed to the LSM. The NN never touches the estimate directly.

**Baseline it improves on:** **batch least-squares (LSM)** — Gauss-Newton over
the STM-propagated residuals (Eqs.1-10, Algorithm 1, pp.2-4), threshold 1e-6.
The proposed pipeline is **BiLSTM-LSM** (Algorithm 2, p.18): run LSM → build
feature sequences → BiLSTM flags outliers → compensate → re-run LSM. A second
baseline, **LSAR** (least sum of absolute residuals, an LP-based robust
estimator, ref 77, Eqs.29-33), is also compared.

**Training-data generation (pp.11-12, 15):** discretise the 69 ephemeris
trajectories at 6 h spacing → 7590 candidate points; each point seeds one OD
scenario with randomised ToF (2-6 h), n measurements (121-361), and initial
estimate error (Table 5, p.12; sampled by Latin hypercube). After visibility
filtering, **3761 valid samples** (each a full RA/DEC sequence). Split 80/20
(3009 train / 662 val). Trained in PyTorch, 100 iterations, lr 1e-4 decayed
0.98/iter, ~2 h on an RTX 4090. A separate **3629-sample** set
(~905,018 individual measurements) is generated fresh for testing.

---

## 3. Quantitative results (pp.16-22)

**Classification accuracy (held-out test, Tables 8-9, pp.16-17):**
- RA BiLSTM overall accuracy **99.06%**; TPR accurate 99.62%, TPR inaccurate
  (i.e., outliers caught) **87.39%**, FNR inaccurate 12.61%.
- DEC BiLSTM overall **99.00%**; TPR accurate 99.60%, TPR inaccurate **86.40%**,
  FNR 13.60%.
- Plain English: it almost never falsely condemns a good measurement, and it
  catches ~87% of the genuinely bad ones (misses ~13%).

**Bidirectional vs unidirectional (Table 10, p.17):** BiLSTM beats forward-only
and backward-only LSTM on the hard (inaccurate) class by a wide margin —
inaccurate-class TPR: BiLSTM 87.4%/86.4% (RA/DEC) vs forward 83.6%/34.5% vs
backward 84.5%/64.3%. The DEC unidirectional collapse (34.5%) is the headline
argument for going bidirectional.

**OD accuracy improvement (the money table — Table 12, p.21; Fig.16, p.22),
300 Monte-Carlo runs on the held-out 4:1 NRHO:**

| Meas. duration | LSM σ_R (km) | BiLSTM-LSM σ_R (km) | LSM σ_V (m/s) | BiLSTM-LSM σ_V (m/s) |
|---|---|---|---|---|
| 2 h | 122.73 | 108.02 | 34.01 | 30.33 |
| 3 h | 57.02 | 52.21 | 15.46 | 14.15 |
| 4 h | 28.73 | 26.50 | 8.86 | 8.00 |
| 5 h | 20.47 | 18.61 | 5.62 | 5.11 |
| 6 h | 14.46 | 13.31 | 4.62 | 4.18 |

- **Reduced ratio (Eq.28, Fig.16):** position error cut **11.99% → 7.96%**
  (2 h → 6 h); velocity **10.82% → 9.56%**. Headline = **"~10%"** (abstract,
  conclusions). Improvement shrinks as more data accumulate (outliers matter
  less with more measurements).
- vs **LSAR** (Table 12): BiLSTM-LSM beats it on every duration; LSAR is in fact
  *worse than plain LSM* here, because LSAR is built for gross outliers (>10σ)
  and the paper's "inaccurate" measurements are mild (2-4σ), not outliers.

---

## 4. Reusability against the ONLY three candidate consumers

### (a) #76 moon-tour patched-conic lane — DOES NOT MAP

The only superficial contact is the shared phrase "Earth-Moon libration orbit."
Concretely:
- **Dynamics/ephemeris treatment:** their model is **VSOP2000** (an analytic
  planetary theory) for Sun/Moon, integrated in full 3-body ephemeris from a
  CRTBP seed. #76 Tier-1 uses a *circular planet-centred moon ephemeris* and a
  patched-conic/Lambert model; #76 Tier-2 uses CR3BP. **No constant, frame, or
  convention transfers.** The paper gives no μ value, no Earth-Moon distance, no
  mass parameter we could cite — it leans entirely on VSOP2000/DE-class
  ephemeris by reference (refs 64-66) and on CRTBP families from Daniel's thesis
  (ref 71) and the Howell-Pernicka two-level corrector (ref 73).
- **The L2 Halo / Lyapunov / NRHO tables (Tables 3-4, pp.9-10):** these are
  **CRTBP nondimensional initial conditions** (x0, z0, ẏ0, period T, Jacobi J)
  sourced from Daniel 2006 (ref 71) and Zhou (ref 72) — i.e., they are
  *re-printed* from prior sources, not new, and they are CRTBP-internal, not
  cycler rows. If #76 ever wants an L2-family seed catalogue, the *primary*
  sources (refs 71, 72) are the citation, not this paper.
- **Verdict:** nothing usable. Wrong dynamical regime (cislunar OD vs Earth-Mars
  cycler design), wrong problem (estimation vs design), no sourced constants.

### (b) Our validation / ML-surrogate culture — CORROBORATES, does not extend

This is the only place the paper earns its read. Its architecture is a textbook
instance of our governing rule (`docs/notes/2026-06-07-ml-surrogate-investigation.md`:
*"a surrogate may PRUNE and PROPOSE but may never constitute EVIDENCE"*):
- The BiLSTM **only flags** suspect measurements; the **classical LSM remains the
  estimator** (Algorithm 2, p.18). The NN never produces the answer — it pre-
  filters the input. That is exactly Paper 3's sound design in our ML note (ANN
  ranks; full OCP decides) and exactly our prune-but-never-testify line.
- Its evaluation discipline is clean and worth noting as a positive example:
  **train/val/test on disjoint sample sets** (3761 train+val, separate 3629
  test, p.15), and crucially **the headline OD evaluation is on a held-out orbit
  family the network never trained on** (the 4:1 NRHO, §4 p.18) — a genuine
  generalisation test, not in-distribution self-grading. Error is reported as
  MC-derived 3σ STDs over 300 runs with explicit reduced-ratio (Eq.28), and the
  classifier's failure mode is reported honestly (FNR ~13% — it *misses* bad
  measurements, the safe direction for a pre-filter).
- **One honest caveat for our culture, not a transfer:** their "compensation"
  step (Eqs.23-25) lets the NN's classification *modify the measurement data*
  before estimation. That is a half-step past pure pruning — the NN's output
  perturbs an input the estimator then trusts. They bound the risk (nudge = 1σ
  only), but by our stricter rule this would be the line we do **not** cross: we
  would prune/down-weight a flagged measurement, never *edit its value* and feed
  the edit forward as if real. Useful contrast to sharpen our own boundary.
- **Verdict:** no code or method to adopt; it is a clean external corroboration
  of the advisory-only stance we already hold, plus one boundary-sharpening
  contrast (don't let the surrogate rewrite the evidence).

### (c) Observability insight on degenerate angle-only data — STRETCH; mostly absent

The task flagged this as possibly informative-but-a-stretch. Honest finding:
**it is a stretch, and the paper does not deliver the geometric argument.** There
is no observability-Gramian/rank analysis, no "when is angle-only degenerate"
cone geometry, no condition-number study. The paper's notion of "bad evidence"
is purely statistical (measurement error > 2σ), detected by a learned sequence
model — not geometric. The only transferable *idea* is the framing that
**evidence quality is per-observation and sequence-dependent** (a measurement's
trustworthiness depends on its neighbours in time, which is why bidirectional
context helps). That loosely rhymes with our own "evidence quality" thinking but
gives no method, threshold, or geometric criterion we could lift. Do not
overclaim it.

---

## 5. Catalogue eligibility & v4.2 flags

**Catalogue eligibility: NONE.** No trajectory is designed, enumerated, or
optimised; the paper estimates the state of pre-existing periodic orbits from
tracking data. Nothing is golden-eligible. (Confirmed in one line, as scoped.)

**v4.2 backfill flags: n/a.** No cycler rows trace to this paper.
- *center:* Earth-centered equatorial inertial (ECI) for OD; CRTBP rotating
  frame (Earth-Moon) for the orbit-family seeds. Not heliocentric — no overlap.
- *tof_days_bounds:* n/a — "time of flight" here is the **observation window**
  (2-6 h), not a transfer duration.
- *source_ephemeris:* **VSOP2000** (analytic planetary theory; refs 64-66) for
  the OD dynamics; CRTBP seeds from Daniel 2006 (ref 71). No DE4xx used. Not
  relevant to any catalogue row, so nothing to record.

---

## 6. Honest "not extractable / not present" list

- No Earth-Mars content of any kind; entirely Earth-Moon cislunar OD.
- No formal observability analysis (the hoped-for degeneracy geometry is absent).
- No reusable constants (μ, distances, masses) — all dynamics are by-reference
  (VSOP2000) or re-printed CRTBP ICs from refs 71-72.
- No trajectory, Δv, sequence, or transfer is produced — the output is a state
  estimate's error statistics, not a design.
- The CRTBP L2-family ICs (Tables 3-4) are sourced from Daniel 2006 / Zhou; if
  ever wanted for #76 Tier-2, cite those primaries, not this paper.

---

## 7. Bottom line

Faithfully mined; honestly out of scope. The paper is a competent cislunar
angle-only **orbit-determination** method (BiLSTM measurement-outlier detector +
LSM), with a clean ~10% OD-error reduction on a held-out 4:1 NRHO. It designs
no trajectories, supplies no cycler-usable constants, and shares with #76 only
the words "Earth-Moon libration orbit." Its one lasting value to us is as a clean
external example of the **advisory-only / prune-not-testify** ML discipline our
surrogate note already mandates — plus a boundary-sharpening reminder not to let
a surrogate *rewrite* the evidence (their Eq.23-25 compensation is one step
further than we would go). Nothing to build, nothing to catalogue.
