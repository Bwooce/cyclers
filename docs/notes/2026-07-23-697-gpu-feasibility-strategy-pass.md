# #697 — Apple-GPU (Metal/MPS) feasibility strategy pass (2026-07-23)

Analysis-only (no production code touched, no catalogue/OUTSTANDING writes), mirroring
`#679`/`#686`'s format. User question: "consider if any of our orbit search is likely to be
able to be done on the apple gpus instead of cpus." Per the standing strategy-pass discipline,
the tractability verdict is delivered first and honestly.

Inputs actually read/run (not summarized from memory): a very-thorough read-only survey of
`src/cyclerfinder/{search,core,nbody,parallel}` (spot-verified for every load-bearing claim
below); `docs/notes/2026-06-06-performance-profile.md` in full; `data/OUTSTANDING.md` wall-clock
history (`#520`, `#618`, the 3.5 h `method="lm"` stall, the days-scale real-eph runs); a small
standalone grounding benchmark run on this machine 2026-07-23 (scratchpad `bench_697.py`,
disposable, results reproduced verbatim below; not committed — it contains nothing beyond what
this note records); `numpy.show_config()`/`scipy.show_config()` on the live venv; and five
targeted web checks on the mid-2026 Apple-GPU tooling landscape (PyTorch MPS, jax-metal,
jax-mps, MLX, Metal fp64 emulation — sources at the end).

## Verdict, delivered first

**No. No part of this pipeline should be ported to the Apple GPU, and the blocker is
structural, not effort: Apple GPUs have no hardware float64, every Python framework that
targets them confirms it, and this project's numerics are float64-mandatory at
`rtol=atol=1e-12`.** Additionally, on this specific machine (base M3, 8-core CPU, 16 GB) the
CPU's Accelerate/AMX units already deliver more fp64 dense-linear-algebra throughput than the
GPU could reach even via software-emulated fp64 — the GPU is not merely hard to use here, it is
*slower* at the only precision this project can accept. The genuinely available speedups are
CPU-side and already partially scoped (`2026-06-06-performance-profile.md`); they are listed in
"What to do instead."

## The hard constraint: no float64 on Apple GPUs, anywhere in the stack (verified mid-2026)

- **Metal hardware**: Apple-silicon GPUs have no fp64 ALUs. The only route is software
  emulation (double-float / extended-precision tricks, e.g. the `metal-float64` project), at
  **~1/32–1/64 of fp32 throughput** with real IEEE-compliance caveats.
- **PyTorch MPS backend**: float64 is rejected outright ("Cannot convert a MPS Tensor to
  float64 dtype as the MPS framework doesn't support float64"), still the case as of May 2026
  issue traffic. fp32/fp16/bf16 only.
- **jax-metal**: experimental, no fp64, and effectively stalled — recent JAX versions
  (0.8.x) fail on it; developer-forum threads openly ask whether it is abandoned.
- **jax-mps** (new, alpha v0.10.10, 2026-07-18; a PJRT plugin lowering StableHLO onto MLX):
  explicitly documents "does not support float64 (unsupported on Metal)". Single-device, alpha.
- **MLX** (Apple's own framework): no float64 on GPU (open feature request, ml-explore/mlx
  #799); fp64 exists only as a CPU-stream fallback. No numerical-ODE ecosystem beyond
  ML-oriented fixed-step `odeint` toys (flow-matching repos); nothing with adaptive-step
  fp64 correctness guarantees.

Why fp32 is not an option for this project's core loops: machine epsilon is ~1.2e-7, five
orders of magnitude above the `1e-12` tolerances the correctors, gates, and V-tier evidence
runs are built on. The whole orbit-closure discipline ("it closed!" is the danger signal;
independent cross-checks; interval certificates) is premised on trustworthy residuals at
1e-9..1e-12; fp32 residual floors would silently convert closure claims into noise. A
fp32-GPU result could never satisfy a gate — it could only ever be a triage layer in front of
a fp64 CPU re-run (see "considered and rejected" item 1).

## Measured baseline on this machine (M3, 8-core CPU, 16 GB unified; bench_697.py, 2026-07-23)

The venv's numpy 2.4.6 **and** scipy 1.17.1 are both built against **Accelerate** (confirmed
via `show_config()`), so the AMX/Neural-adjacent matrix units are already engaged — the "free"
Apple-silicon speedup is already captured, without any GPU work:

| Kernel | Result |
|---|---|
| fp64 GEMM n=512 / 2048 / 4096 | 353 / 301 / 271 GFLOPS |
| fp32 GEMM n=2048 / 4096 | 1022 / 857 GFLOPS |
| dense `np.linalg.eig`, n=1000 / 3000 (fp64) | 0.53 s / 10.2 s |
| `solve_ivp` DOP853 CR3BP 6-state, 2π, rtol 1e-9 / 1e-12 | 42 ms (8.3k RHS) / 69 ms (13.5k RHS) |
| `solve_ivp` DOP853 CR3BP+STM 42-state, 2π, rtol 1e-9 / 1e-12 | 149 ms (10.6k RHS) / 268 ms (17.3k RHS) |
| Per-RHS cost incl. stepper overhead | ~5 µs (6-state), ~15 µs (42-state), pure-Python/numpy callback |

The comparison that kills the GPU case quantitatively: the base-M3 GPU peaks at ~2.2–2.8
TFLOPS **fp32**; at the 1/32–1/64 emulation penalty its fp64 ceiling is **~35–90 GFLOPS** —
i.e. **3–8x slower than the CPU's measured 271–353 GFLOPS fp64 GEMM**, before counting
emulation-library integration cost and kernel-launch overhead on 6×6/42-dim problems. Even in
pure fp32 the GPU's headroom over this CPU's ~0.9–1.0 TFLOPS Accelerate GEMM is only ~2–3x.

## Workload-by-workload verdicts

1. **Adaptive `solve_ivp` propagation (DOP853/Radau, 6/42-dim, pure-Python RHS)** — the
   dominant per-trajectory pattern (~842 call sites). *Worst possible GPU shape*: sequential
   data-dependent step control, 8k–17k tiny RHS evaluations per period, 42-dim state,
   fp64-mandatory. Kernel-launch latency alone (~µs-scale per dispatch) would exceed the
   entire RHS cost. **No.**
2. **Differential correction / shooting** (`least_squares` + re-integration, 6×6–~40×40
   solves): control-flow-divergent (iteration counts are data-dependent per point — exactly
   what made `#520`'s sweep 9–13.5 h serial), tiny matrices, fp64 gates. **No.**
3. **Grid sweeps (`scan_parallel`, `parallel_sweep.py`)** — embarrassingly parallel *across*
   points, but each task is a full adaptive corrector (pattern 2), and the measured profile
   says the cost center is **`Ephemeris.state()` at 87–89% of solve time** (an
   interpolation/table workload), with process scaling already memory-bandwidth-limited (3.30x
   at 8 workers, 41% efficiency). A GPU port would have to reproduce the DE440 Chebyshev
   evaluation + the whole corrector in fp64 on hardware with no fp64. **No** — but this is
   precisely where the *CPU* headroom lives (below).
4. **Pseudospectral torus correctors** (`variational_{qp,qbcp,ccr4bp}_torus.py`) — the most
   GPU-*shaped* kernel in the repo (dense collocation Jacobians, numpy-vectorized einsum RHS
   over the grid), but small: a converged `n1=28, n2=9` torus is ~6.5k unknowns, and the dense
   TRF factorizations at that size run in seconds-to-minutes on a 300-GFLOPS fp64 CPU. Emulated
   GPU fp64 at ≤90 GFLOPS would be a slowdown; fp32 would destroy the corrector. **No.**
5. **Validated Taylor / interval integration** (`scripts/_validated_taylor_integrator.py`,
   mpmath.iv arbitrary-precision directed rounding): scalar, branch-heavy,
   arbitrary-precision, correctness-critical — the survey's prior "anti-GPU" label is
   confirmed. **No, categorically.**
6. **Set-oriented / GAIO transfer operator** (`set_oriented_transfer_operator.py`) — flagged
   in the dispatch as the most promising candidate; on inspection it is not: the cost is
   dominated by the **box-mapping** phase (n_boxes × n_samples_per_box independent *short
   adaptive integrations* — i.e. pattern 1 again), while the linear-algebra phase is one
   small-k `scipy.sparse.linalg.eigs` (or dense eig under the 2000-box cutoff, ~seconds). The
   GPU-friendly part is the part that is already cheap. **No at current scale** (revisit
   trigger below).

## What to do instead (the real, already-identified headroom is CPU-side)

Ranked by evidence strength; none of these are GPU work, and items 1–2 were already scoped in
`2026-06-06-performance-profile.md` §"Recommendations" but never executed:

1. **Ephemeris epoch-grid memoisation + batched `Time`/posvel evaluation** — `state()` is
   87–89% of a solve; one S1L1 solve makes 120 `state()` calls with only 70 distinct
   (body,epoch) pairs; the profile note estimates **~−70% solve wall-clock**. The batch path
   (`nbody/forces.py` `Ephemeris.states()`) already exists and is under-used. This directly
   attacks the workload class (#520-style corrector sweeps) that has actually burned the
   project's wall-clock hours.
2. **numba-`@njit` the force-model RHS kernels** (CR3BP/BCR4BP/QBCP/CCR4BP EOM+STM). numba
   0.65.1 is already a dependency, already used for Lambert/Kepler/Stumpff — but *not* for the
   force fields, which are pure-Python/numpy at ~5–15 µs/RHS. A jitted 42-dim RHS is
   realistically sub-µs; with DOP853's own Python stepper overhead remaining, a **~3–10x**
   per-propagation gain is plausible for the CR3BP-family models (not the astropy-bound
   real-eph path, which is item 1's territory). Cheap, incremental, zero correctness-model
   change (same integrator, same tolerances) — but still requires re-validation discipline per
   [[feedback_bugfix_invalidates_past_searches]] if any RHS is rewritten.
3. **If batch propagation ever becomes the true bottleneck**: the right tool is a compiled
   *CPU* ensemble integrator in fp64 (e.g. heyoka.py's LLVM Taylor integrator with SIMD batch
   mode — built for exactly this astrodynamics-ensemble use case, fp64-native, runs on
   arm64), not a Metal port. Noted as an option, not a recommendation — nothing currently
   justifies the integration-stack churn.
4. **Accept the parallel-efficiency ceiling knowingly**: the measured 41% efficiency at 8
   workers is memory-bandwidth contention; more workers on this box will not help, and a GPU
   shares the same unified-memory bandwidth.

## Considered and explicitly rejected (do not re-surface without new evidence)

- **fp32 GPU triage layer** (massive fixed-step fp32 batch propagation on MLX/jax-mps for
  FTLE/box-transition/periapse-map *statistics*, with fp64 CPU confirmation of anything
  interesting): the only defensible GPU lane in principle — measure-level outputs can tolerate
  fp32 trajectory divergence. Rejected on three grounds: (a) the workloads it would serve
  (set-oriented box-mapping, FTLE scoring) are not current bottlenecks — coarse pilot grids run
  in minutes on CPU; (b) the payoff ceiling on this hardware is ~2–3x even for perfectly
  GPU-shaped fp32 work (GPU ~2.2–2.8 TFLOPS vs CPU's measured ~1 TFLOPS fp32), while the port
  cost is a full force-model rewrite in a framework with alpha-maturity numerics; (c) every
  candidate framework (MLX, jax-mps) is single-device alpha with ML-, not
  dynamics-grade, validation. Cost/payoff is upside-down.
- **Emulated-fp64 Metal kernels** (`metal-float64`-style): ≤90 GFLOPS ceiling vs 271–353
  measured on CPU; a strictly negative-return engineering project.
- **PyTorch-MPS or jax-metal ports of anything**: no fp64, and jax-metal is effectively
  unmaintained (breaks on current JAX).
- **GPU for the torus correctors' dense algebra**: right shape, wrong size (~6.5k unknowns)
  and wrong precision; CPU Accelerate already handles it.
- **GPU sparse eig for the transfer operator**: the eig is seconds; the box-mapping is the
  cost, and it is pattern-1 work.
- **Buying the speedup via more CPU processes**: measured bandwidth ceiling (3.30x @ 8
  workers) says no.

## Revisit triggers (conditions under which this verdict should be re-examined)

1. **Hardware change**: access to a machine with native-fp64 GPU compute (any NVIDIA/AMD HPC
   part; fp64 at 1:2–1:32 of fp32) changes the entire calculus — batched STM propagation and
   box-mapping at 10⁵–10⁶ trajectories are then genuinely portable (CUDA has a mature
   ecosystem for exactly this). The *workload analysis* in this note transfers; the verdict
   does not.
2. **Apple ships hardware fp64 or an fp64 MLX GPU path** (track ml-explore/mlx #799).
3. **A measure-level workload becomes the actual bottleneck at ≥100x current scale** (e.g. a
   set-oriented campaign needing ≥10⁶ box-samples per grid): the fp32-triage lane above then
   deserves a fresh cost/payoff pass, still with fp64 CPU confirmation as the gate.

## Sources (web checks, 2026-07-23)

- PyTorch MPS fp64 rejection: [pytorch forums](https://discuss.pytorch.org/t/typeerror-cannot-convert-a-mps-tensor-to-float64-dtype-as-the-mps-framework-doesnt-support-float64-please-use-float32-instead/180852), [docling issue, May 2026](https://github.com/docling-project/docling/issues/3483), [Apple dev forums](https://www.developer.apple.com/forums/thread/797778)
- jax-metal status: [jax#16435 fp64 crash](https://github.com/jax-ml/jax/issues/16435), [jax#34109 recent-macOS failure](https://github.com/jax-ml/jax/issues/34109), [Apple dev forums "is anyone working on jax-metal?"](https://developer.apple.com/forums/thread/743350)
- jax-mps (MLX-backed PJRT, alpha, no fp64): [PyPI](https://pypi.org/project/jax-mps/)
- MLX fp64 feature request: [ml-explore/mlx#799](https://github.com/ml-explore/mlx/issues/799)
- Metal fp64 emulation at 1/32–1/64 throughput: [philipturner/metal-float64](https://github.com/philipturner/metal-float64)
