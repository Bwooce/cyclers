# #698 — GPU-native N-body algorithms (MCPI / symplectic / parareal): follow-up to #697 (2026-07-23)

Analysis-only (no production code touched, no catalogue/OUTSTANDING writes), same discipline as
`#679`/`#686`/`#693`/`#697`. User question after `#697`'s negative: *"what about using the new
apple gpu libraries? and finding a new algorithm which does fit gpus well to solve n-body
problems?"* This pass answers both halves: (1) a fresh check on Apple fp64 GPU status, and
(2) whether a fixed-computational-shape algorithm — Modified Chebyshev-Picard Iteration (MCPI),
fixed-step symplectic integration, or parallel-in-time methods — sidesteps `#697`'s
adaptive-step-control blocker for THIS project's dynamics regime.

Inputs actually read/run: `#697`'s note in full (its measured numbers are load-bearing below);
targeted web checks 2026-07-23 AET (sources at end); and a small disposable Picard-Chebyshev
proof-of-concept run on this machine against the project's own dynamics class (planar Earth-Moon
CR3BP, scratchpad `mcpi_poc_698.py`, results reproduced verbatim below; not committed — this note
records everything it produced).

## Verdict, delivered first

**No — the GPU-native-algorithm route does not rescue Apple-GPU acceleration for this project,
and the reason is arithmetic, not algorithmics: even a perfectly GPU-shaped, fixed-iteration,
batch-uniform fp64 algorithm runs on an Apple GPU at the ~35–120 GFLOPS emulated-fp64 ceiling,
below the 271–471 GFLOPS this machine class's CPU already delivers on exactly the GEMM-shaped
kernels MCPI is built from.** The dispatch's premise — "GPUs can win even at the emulated-fp64
penalty for sufficiently large, uniform batches" — is quantitatively false on Apple silicon,
because the correct CPU comparison point is Accelerate/AMX GEMM, and the same batching that
would fill the GPU also fills AMX. The algorithm question and the Apple-GPU question therefore
decouple: no algorithm choice can fix a raw fp64-throughput deficit at the only precision the
project's 1e-12 gates accept.

**On the algorithms themselves, judged on their merits (i.e. for future native-fp64 GPU
hardware, per `#697` revisit trigger 1):** MCPI is real and has genuine CR3BP prior art — but my
own PoC on Earth-Moon CR3BP confirms the literature's caveat pattern: it converges cleanly on
benign and moderate (L1-region) arcs, then **fails exactly where this project's science lives**
— close approaches force a ≥40x segment-length shrink, and worse, the discrete Picard operator
exhibits **spurious fixed points**: runs that "converge" to 1e-12 fixed-point tolerance while
being wrong by 1e-3 to 2.0 in state. That is a structural "it closed!" hazard for a codebase
whose entire trust discipline is built on not believing convergence flags. Fixed-step symplectic
(Wisdom-Holman) is a category mismatch: it targets long-term near-Keplerian secular integration,
breaks down at close encounters by construction, and its selling point (bounded energy drift)
is not the currency the correctors/gates trade in (local 1e-12 trajectory + STM accuracy).
Parareal/PFASST: the dispatch's skepticism is confirmed by the literature — parareal diverges
for many chaotic systems and speedup saturates at low processor counts.

## Part 1: Apple fp64 GPU status re-check (2026-07-23) — unchanged, negative

- **No hardware fp64 on any Apple GPU generation including M5.** A March-2026 M5-GPU roofline
  analysis measures fp32 at ~3.85 TFLOPS (scalar, tuned) and contains no GPU fp64 path at all;
  its only fp64 number is the *CPU AMX*: `cblas_dgemm` 471 GFLOPS. Metal still exposes no fp64
  type; the only route remains software emulation (`metal-float64`-class, ~1/32–1/64 of fp32).
- **No framework change**: PyTorch MPS still rejects float64 outright (mid-2026 threads);
  MLX fp64-on-GPU remains an open feature request (ml-explore/mlx #799); jax-mps still documents
  "does not support float64 (unsupported on Metal)". Nothing new since `#697`'s snapshot of five
  days ago, as expected — it is a silicon gap, and the M5 generation did not close it.
- **The M5 datapoint strengthens `#697`'s conclusion**: emulated fp64 from 3.85 TFLOPS fp32 is
  ~60–120 GFLOPS, vs 471 GFLOPS measured on the same chip's CPU/AMX. The CPU-side fp64 advantage
  *widened* with the new generation, not narrowed.

## Part 2a: the arithmetic that decouples "GPU-native algorithm" from "Apple GPU"

MCPI's inner loop is precisely fixed dense linear algebra: RHS evaluation at N Chebyshev nodes,
then multiplication by precomputed (N+1)×(N+1) integration/transform matrices; batching B
trajectories turns it into GEMM of shape (N+1)×(N+1)×(6B). This is the *best possible* GPU shape
— and also the exact kernel `#697` measured on this CPU at **271–353 GFLOPS fp64** (Accelerate/
AMX, numpy+scipy both confirmed built against Accelerate). Apple-GPU emulated fp64 tops out at
**~35–90 GFLOPS (M3) / ~60–120 (M5)**. A perfectly uniform fp64 MCPI batch on the Apple GPU is
therefore a guaranteed **3–8x slowdown** vs running the identical batch through the CPU's BLAS —
before counting the cost of writing/validating an emulated-fp64 Metal linear-algebra stack that
does not currently exist. fp32 MCPI is excluded for the same reason as ever: 1.2e-7 epsilon vs
1e-12 gates. This holds regardless of batch size B, because both devices are FLOP-bound on the
same GEMM. **No fixed-shape algorithm can escape this; the blocker `#697` called "structural"
survives the algorithm change.**

(All published GPU-MCPI successes are CUDA on NVIDIA hardware with *native* fp64 — Bai &
Junkins' original CUDA implementation; Masat/Colombo/Boutonnet 2023 ran fp64 at rtol 1e-12 on a
GTX 1050 for 15.8x over sequential C on 13,509 trajectories. The literature confirms the
algorithm is GPU-ready; it says nothing about GPUs without fp64.)

## Part 2b: MCPI on this project's dynamics — literature + PoC

### What the literature actually supports

- **Genuine CR3BP prior art exists** (better than expected): Swenson, Woollands & Junkins 2017
  (J. Astronautical Sciences) applied MCPI to differential correction of **45,000+ halo orbits
  across three CR3BP systems**, up to ~10x faster than RK7/8 *with increased robustness*, using
  Chebyshev interpolation to target plane crossings instead of integrating the 42-dim
  variational equations. Bai & Junkins 2012 used MCPI for translunar-halo station-keeping.
  MCPI is not confined to weakly-perturbed Kepler catalogs.
- **But the convergence-domain caveats are structural**: the Picard domain of convergence in
  Cartesian variables is ~**1/3 of an orbit** and eccentricity-dependent; Woollands et al.
  needed **Kustaanheimo-Stiefel regularization** to extend it to ~90% of an orbit for perturbed
  Lambert problems. The production-grade "Adaptive Picard-Chebyshev" line (Woollands & Junkins,
  JGCD 2019) exists precisely because segment length and polynomial degree must be **tuned per
  trajectory segment** — i.e., the mature version of MCPI reintroduces data-dependent
  adaptivity. Masat 2023's GPU version needed 41–42 iterations, ~0.87-period segments at 200
  nodes, and Keplerian **warm-starting**.

### PoC on planar Earth-Moon CR3BP (this machine, 2026-07-23; vanilla cold-start Picard, CGL nodes, tol 1e-12, truth = DOP853 at 1e-13)

| case | dt (TU) | min Moon dist | N | iters | converged | max err vs DOP853 |
|---|---|---|---|---|---|---|
| benign Earth orbit | 0.5–2.0 | 0.5 | 64 | 19–31 | yes | 2e-12–5e-12 |
| L1 region | 0.5–1.0 | 0.13–0.15 | 64 | 17–23 | yes | 4e-13–6e-13 |
| L1 region (arc dips to 0.018 of Moon) | 2.0 | 0.018 | 64 | — | **no** | — |
| L1 region, same arc | 2.0 | 0.018 | 128 / 200 | 83 / 68 | **"yes"** | **2.8e-2 / 2.2e-3** |
| close lunar approach (r₂≈0.004) | 0.5–2.0 | 0.004 | 64–200 | — | **no** | — |
| close approach, shorter segments | 0.2–0.3 | 0.004 | 64–128 | 67–69 | **"yes"** | **1.7–2.2 (!)** |
| close approach, dt=0.10 | 0.10 | 0.004 | 128 | 46 | yes | 8.9e-4 |
| pre-approach only (never reaches periapsis) | 0.05 / 0.02 | 0.015 / 0.028 | 64 | 21 / 13 | yes | 4e-13 |

Three findings, all consistent with the literature and directly relevant to this project:

1. **MCPI works where the field is smooth** — including the L1 region, matching Swenson 2017's
   positive halo result. Benign arcs of ~2 TU (≈8.7 days) converge at N=64 in ~20–30 iterations
   to 1e-12-class agreement with DOP853.
2. **Close approaches destroy the uniform-shape premise.** Segments containing a
   r₂≈0.004 lunar passage either diverge outright (N up to 200, dt down to 0.2) or need
   dt≪0.05 — a **≥40x segment shrink** relative to the benign case. Batch members near vs far
   from the secondary then need wildly different segmentation — exactly the per-trajectory
   heterogeneity (different effective step counts) that the GPU pitch needed to eliminate. The
   known cures (KS/Levi-Civita regularization — conceptually the same trick as `#670`/`#671`'s
   `scripts/_validated_taylor_integrator.py` W-Z arc; adaptive segment/degree tuning; warm
   starts) all reintroduce data-dependent, per-trajectory control flow, and time-regularization
   additionally desynchronizes batch members in physical time.
3. **The discrete Picard operator has spurious fixed points — a structural "it closed!"
   hazard.** Multiple runs above satisfied the 1e-12 fixed-point convergence criterion while
   being wrong by 1e-3 to **2.0** (state units) against DOP853: when the polynomial basis cannot
   resolve the near-singular feature, the iteration converges to the fixed point of the
   *discretized* operator, which is far from the true solution, with no internal warning. For a
   codebase whose gates are built on distrusting convergence flags, adopting MCPI would require
   a full independent-verification methodology (every MCPI result re-checked by an adaptive
   fp64 integrator — i.e., the CPU pipeline it was meant to replace) plus re-derivation of the
   STM/covariance machinery in the Chebyshev basis. Caveat, stated honestly: this PoC is
   *vanilla* cold-start Picard without the error-feedback/quasi-linearization terms of
   production MCPI, which would improve iteration counts and enlarge the domain somewhat — but
   the convergence-domain mechanism it demonstrates is the same one the Woollands literature
   documents and engineers around, and the spurious-fixed-point failure mode is basis-resolution
   driven, which error feedback does not remove.
4. **No usable implementation exists to adopt.** No maintained GPU-capable (or even CPU) Python
   MCPI package surfaced; the ecosystem is research CUDA/C/MATLAB code, largely unreleased
   (Masat 2023 publishes no repository). This would be a from-scratch build plus from-scratch
   validation.

## Part 2c: fixed-step symplectic (Wisdom-Holman / MVS) — category mismatch

- **The GPU-shape claim is only half true.** Pure fixed-step WH is uniform per step and
  batchable — GENGA (CUDA, fp64 Kepler solver) proves the pattern on real hardware. But WH
  requires a near-Keplerian splitting (dominant central body + small perturbation); the
  project's rotating-frame CR3BP/BCR4BP/CCR4BP corrector work near libration points and
  separatrices is not that regime, and the perturbation term is not small exactly where it
  matters.
- **Close encounters are the canonical WH failure**, stated plainly in the GPU codes' own
  papers: a symplectic integrator "cannot reduce the timestep during close encounters without
  destroying the symplectic properties"; every production code (MERCURY-style hybrids, QYMSYM,
  GENGA) handles encounters by **switching integrators mid-flight** (smooth changeover to
  Bulirsch-Stoer/direct N-body) — data-dependent branching that re-imports the divergence
  problem into the batch, which is why GENGA's engineering is dominated by encounter bookkeeping.
- **Wrong accuracy currency.** Symplecticity buys bounded energy/structure drift over 1e6+
  orbits — valuable for planet-formation statistics, irrelevant to a differential corrector
  that needs |residual| ≤ 1e-12 *pointwise* plus a consistent STM after ~1–10 periods. A fixed
  step small enough for 1e-12 local accuracy through a close approach, applied globally,
  multiplies the FLOP count by orders of magnitude vs DOP853's adaptive grid — negating the
  batch win. Torus correctors, manifold globalization, and the V-gauntlets would all need their
  precision guarantees re-derived on a method that gives them nothing in return.
- The repo's Levi-Civita work (`#670`/`#671`) is the *conceptual* answer to close-approach
  stiffness and is already deployed on the CPU side where it belongs (interval-validated
  arcs); porting it into a batch-uniform GPU setting founders on the physical-time
  desynchronization noted above.

## Part 2d: Parareal/PFASST — dispatch skepticism confirmed

Brief, per the dispatch's own priority: the literature is unambiguous that vanilla parareal
degrades or **diverges for chaotic systems** (LLNL's MGRIT-for-chaos report states parareal
"diverges in many cases" and MGRIT/XBraid "struggles"; Rayleigh-Bénard and cylinder-flow studies
show degradation with increasing Reynolds number; observed speedups saturate at ~20 processors
as iteration counts climb). Recent research (finite-time convergence guarantees, weighted
proximity functions on Lorenz-96, 2026 preprints) improves *statistical/shadowing-sense*
convergence — not the 1e-12 pointwise trajectory accuracy the gates require. Also parareal
parallelizes across *time* for a *single* trajectory; this project's abundant parallelism is
across *trajectories* (sweeps), where parareal solves a problem we do not have. **Confirmed
lower priority; rejected.**

## Considered and explicitly rejected (do not re-surface without new evidence)

- **Apple-GPU MCPI (emulated fp64)**: 3–8x *slower* than the CPU's Accelerate GEMM on the
  identical kernel; strictly negative return. (Part 2a.)
- **fp32 MCPI triage on Apple GPU**: same rejection as `#697`'s fp32-triage item, now with the
  added spurious-fixed-point hazard — an fp32 Picard fixed point is even less trustworthy.
- **Fixed-step symplectic for corrector/gauntlet work**: category mismatch (Part 2c).
- **Parareal/PFASST**: diverges for chaotic dynamics; wrong parallelism axis (Part 2d).
- **CPU-side batched MCPI as a DOP853 replacement for sweeps**: the only lane with *any*
  residual merit on this hardware (Accelerate GEMM at 300+ GFLOPS on the batch, ~20–30
  iterations for smooth arcs), but rejected *now* because (a) `#697`'s two CPU recommendations
  (ephemeris memoisation ~−70% solve time; numba RHS ~3–10x) are cheaper, validated-in-place,
  and unexecuted — they come first; (b) the sweeps that dominate wall-clock are corrector
  sweeps whose interesting members pass near secondaries, i.e. the regime where MCPI's batch
  uniformity dissolves; (c) the validation build-out (independent cross-check methodology for a
  method with demonstrated silent-wrong-answer modes) is disproportionate to an unproven
  speedup over numba+DOP853.

## Revisit triggers

1. **Native-fp64 GPU hardware access (NVIDIA/AMD)** — unchanged from `#697`, but this note
   *upgrades* it: the Masat/Colombo/Boutonnet two-level Picard-Chebyshev augmentation is the
   proven starting design, and the workload it fits is far-from-secondary batch propagation
   (planetary-moon tour screening, measure-level sweeps), NOT close-approach corrector work,
   which stays on the adaptive CPU pipeline under any hardware.
2. **Apple ships hardware GPU fp64** (track ml-explore/mlx #799) — then re-run Part 2a's
   arithmetic; the convergence findings (Part 2b) still apply and still confine MCPI to smooth
   arcs with independent verification.
3. **A published MCPI variant demonstrating certified convergence through close approaches
   without per-trajectory adaptivity** (none found as of 2026-07-23) would reopen Part 2b.

## Honest caveats

- Paywalled papers (Swenson 2017, Woollands JGCD 2019, Bai & Junkins JAS 2013) were assessed
  from abstracts/indexing services, not full texts; iteration/segment numbers quoted for them
  are from those summaries. Masat 2023 was read via its arXiv HTML.
- The PoC is a deliberately minimal cold-start Picard-Chebyshev (numpy chebfit/chebint per
  iteration, not the fixed-matrix production form); its *performance* is meaningless and was
  not measured — only its convergence behavior, which matches the literature in both the
  positive (L1/halo-class arcs) and negative (close-approach) directions.
- Part 1 web checks are point-in-time (2026-07-23 AET); the M5 roofline source is a third-party
  measurement, not an Apple spec sheet, but it agrees with every framework's own documentation.

## Sources (web checks, 2026-07-23)

- M5 GPU roofline (fp32 ~3.85 TFLOPS; only fp64 datapoint is CPU/AMX dgemm 471 GFLOPS):
  [michaelstinkerings.org](https://www.michaelstinkerings.org/apple-m5-gpu-roofline-analysis/)
- PyTorch-MPS fp64 still rejected: [Apple dev forums](https://www.developer.apple.com/forums/thread/797778); MLX fp64 request still open: [ml-explore/mlx#799](https://github.com/ml-explore/mlx/issues/799); Metal fp64 emulation ~1/32–1/64: [philipturner/metal-float64](https://github.com/philipturner/metal-float64)
- MCPI foundations: [Bai & Junkins, MCPI for IVPs, JAS 2013](https://link.springer.com/article/10.1007/s40295-013-0021-6); [Bai & Junkins, translunar halo station-keeping, 2012](https://onlinelibrary.wiley.com/doi/10.1155/2012/926158)
- MCPI in CR3BP differential correction (45k halos, 3 systems, ≤10x vs RK7/8):
  [Swenson, Woollands, Junkins, JAS 2017](https://link.springer.com/article/10.1007/s40295-016-0110-4)
- Convergence domain ~1/3 orbit Cartesian → ~90% with KS regularization:
  [Woollands et al., JGCD](https://dx.doi.org/10.2514/1.G001028); [multiple-rev perturbed Lambert, JAS 2017](https://link.springer.com/article/10.1007/s40295-017-0116-6)
- Adaptive (per-segment tuned) Picard-Chebyshev: [Woollands & Junkins, JGCD 2019](https://arc.aiaa.org/doi/10.2514/1.G003318)
- GPU fp64 Picard-Chebyshev (CUDA, GTX 1050, 15.8x, 41–42 iterations, warm-started):
  [Masat, Colombo, Boutonnet, Acta Astronautica 2023 / arXiv 2301.03989](https://ar5iv.labs.arxiv.org/html/2301.03989)
- Symplectic close-encounter breakdown + GPU hybrids: [QYMSYM](https://www.sciencedirect.com/science/article/abs/pii/S1384107611000303); [GENGA, ApJ 2014](https://iopscience.iop.org/article/10.1088/0004-637X/796/1/23); [GENGA II, ApJ 2022](https://iopscience.iop.org/article/10.3847/1538-4357/ac6dd2)
- Parareal/MGRIT vs chaos: [LLNL MGRIT-for-chaos final report](https://www.osti.gov/servlets/purl/1994037); [MGRIT for chaotic systems, SISC](https://doi.org/10.1137/22m1518335); [Rayleigh-Bénard PinT performance](https://arxiv.org/pdf/2001.01609); [finite-time PinT guarantees, 2026](https://arxiv.org/html/2604.00855)
