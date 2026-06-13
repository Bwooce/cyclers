# Liang et al. 2024 Member D (`liang-2024-cgcec-ephemeris-2033`) — n-body re-propagation lane (#223)

**Source row:** `liang-2024-cgcec-ephemeris-2033` (V0). **Mining note:**
`docs/notes/2026-06-11-liang-2024-cge-triple-cyclers-mining.md` (§4.4 + TRANSCRIPTION
RESCAN). **Idealized scaffold:** `src/cyclerfinder/search/cge_scaffold.py` (#222,
commit 70a1cc5; results `docs/notes/2026-06-13-liang-abc-reproduction.md`).
**Kernel:** JPL NAIF JUP365 (`jup365.bsp`, NAIF generic kernels; local copy
outside the repo), coverage 1600–2200, bodies 501–504 + 599 confirmed present.

## 0. ACHIEVABLE-CEILING ASSESSMENT (written before building)

**Member D's published record is epoch + sequence + kernel + four scalars.**
Departs Callisto 2033-09-25T18:04:43, returns to Callisto 2036-06-22T01:44:39
(~1000 d, CGCEC × 10 cycles); max required Δv per cycle 1.0383e-7 m/s; all
flybys above 100 km real altitude (p. 19); cycle-ToF band 99.4–100.5 d (figure
trace, Figs. 8d–e). **No state vector, no numeric per-flyby table** — per-flyby
epochs/V∞/altitudes exist only as figure traces (mining note §4.4, rescan item 9).

**Therefore EXACT reproduction of THEIR trajectory is impossible from print.**
Every quantity that would let us re-derive their member (an anchor state, a
per-flyby epoch table, even one numeric V∞) is unpublished. Matching their two
printed epochs to the second cannot be demanded of an independently constructed
trajectory either: the construction is a chained local optimization whose
solution depends on solver internals they do not print.

**The honest ceiling is SAME-FAMILY evidence, not THEIR-member confirmation:**
construct OUR OWN candidate in the same family — CGCEC topology, multi-rev
legs, seeded from the validated idealized scaffold phased to their published
departure epoch on real JUP365 geometry — drive it to ballistic closure, and
compare its qualitative signature (sequence achieved, per-cycle ToF inside
their 99.4–100.5 d figure band, flyby ordering/timing pattern, altitudes
> 100 km, residual-Δv scale) against the published traces.

**Why that is V1-class and NOT V3 (spec §14):** V3-ballistic/-powered both
require the row's encounters to be *independently confirmed on the real
ephemeris* — "the cycler's encounters" means THEIR encounters (in-band miss,
per-leg v∞ match against published values). Member D publishes no encounter
to confirm: there is no per-flyby epoch, no v∞ value, no altitude figure to
match within a tolerance. A same-family candidate of our own construction
demonstrates that the family the row describes EXISTS in the stated model at
the stated epoch with the stated qualitative signature — that is existence /
plausibility evidence for the row's claims, the V1 evidence class ("solver
cross-check" generalised: the row's claim re-derived by an independent
construction), not encounter-level confirmation. **No disagreement with the
tasked ceiling: V1-class is the correct cap.** Promoting beyond V1 would
require the authors' data (the row's `data_gaps` already say exactly this).

A clean failure to close ballistically is likewise a finding (the
idealized→real gap is what this lane measures) and would be recorded without
writeback.

## STATUS: RUN EXECUTED (#223) — 2026-06-13

The lane was run end to end. Driver `scripts/liang_member_d_run.py`; toolkit
`src/cyclerfinder/nbody/jovian.py`, now numerically exercised and covered by
`tests/nbody/test_jovian.py` (7 tests: 4 kernel-free closed-form + 3 slow
kernel-backed — spline accuracy, one real-geometry conic leg, and a
propagator energy-conservation guard). **No jovian.py bug was found** in the
run (the module's equations behaved; see §5 for the one seed-construction
caveat the result exposes). **No writeback is proposed — see the verdict.**

## 1. MODEL AS RUN

* **Kernel** — `jup365.bsp` (1.14 GB, not committed), via `CYCLERFINDER_JUP365`.
  Coverage spans the 2033-09-25 → 2036-06-22 window; bodies 501–504 + 599 read.
* **Frame / axis** — Jupiter-centered J2000 equatorial; TDB seconds since J2000
  (departure epoch 2033-09-25T18:04:43 TDB → t0 = 1064556283.000 s).
* **GMs** — registry (`core/satellites.py`): Jupiter 1.26686534e8 km³/s²;
  Galilean moon GMs (JPL SSD). Departure-epoch moon |r| (km), real JUP365:
  Io 421441, Europa 666291, Ganymede 1069523, Callisto 1896001 — Galilean
  geometry, as expected.
* **Conic chain** — `chain_cycles`: multi-rev Lambert legs (branch plan
  high/low/high/low, all 1-rev) between real moon positions; per-cycle
  Nelder-Mead local optimisation of the four interior flyby epochs (±3 d) to
  minimise the powered-flyby defect Δv; chained cycle-to-cycle with the inbound
  V∞ inherited. Seed ToF = idealized Member A first-cycle legs
  (31.8973, 18.1697, 29.9343, 19.9747 d).
* **n-body** — `JovianRestrictedNBody` (REBOUND/IAS15, ε 1e-11): massless
  spacecraft, central Jupiter point mass, four Galilean moons on JUP365 rails
  (0.02-d cubic-spline cache; pinned <1 km in the smoke test, ~1e-2 km claim)
  with the indirect term; softened at the moon SURFACE only. Per-cycle
  multiple shooting (`shoot_cycle`, `least_squares`/trf, 31 free vars =
  node-0 velocity + nodes 1–4 full states + 4 epoch offsets) targeting leg
  continuity + encounter-band hinges + the boundary V∞-magnitude pin.

Cycles run: **3** of the published 10 (each cycle is an independent
signature test; 3 is decisive at ~2.1 s for the whole conic chain and ~20 min
per shot cycle — runtime, not a modelling, limit). The conic chain is cheap;
the n-body shoot is the cost (one cycle, 60 nfev ≈ 1194 s wall).

## 2. PER-CYCLE SIGNATURE (patched-conic chain on real JUP365)

| cyc | cycle ToF (d) | in 99.4–100.5 | Σ defect Δv (m/s) | max defect Δv (m/s) | min flyby alt (km) | all alts > 100 km |
|----:|--------------:|:-------------:|------------------:|--------------------:|-------------------:|:-----------------:|
| 1   | 100.2832      | yes           | 1.961e-10         | 8.732e-11           | 1469.3             | yes               |
| 2   | 99.9569       | yes           | 4.100e-09         | 2.098e-09           | 15351.3            | yes               |
| 3   | 100.0081      | yes           | 3.675e-09         | 2.024e-09           | 6129.8             | yes               |

Per-cycle interior V∞ (km/s), inbound at C/G/C/E/C:

| cyc | C(out) | G | C | E | C(in→next) |
|----:|-------:|---:|---:|---:|-----------:|
| 1   | 5.5833 | 6.9671 | 5.6530 | 3.8693 | 5.5390 |
| 2   | 5.5390 | 6.8857 | 5.6141 | 4.2809 | 5.6962 |
| 3   | 5.6962 | 7.1757 | 5.7664 | 5.6205 | 5.9801 |

**Conic-level signature MATCHES the published family on every axis:**
CGCEC sequence constructed each cycle; ToF **3/3 inside the 99.4–100.5 d
figure band**; every flyby altitude > 100 km (min 1469 km, cycle 1 Callisto);
powered-flyby defect Δv at the **1e-10 – 1e-9 m/s** scale, i.e. *below* the
published max-required 1.0383e-7 m/s — the chain finds near-ballistic
flyby continuity on real geometry. V∞ cluster (5.5–7.2 km/s) sits in the
idealized A/B family band (A printed C/G/C/E/C = 5.67/6.99/5.67/4.67/5.87).

## 3. N-BODY RE-PROPAGATION (idealized → real continuous-gravity gap)

Seed = patched-conic periapsis nodes (`periapsis_node`) at the converged
conic epochs. **Raw seed defect** (one n-body propagation per leg, no
correction; `--diag`), |Δr| per leg / |Δv| per leg:

| cyc | leg1 Δr (km) | leg2 | leg3 | leg4 | leg1 Δv (m/s) | leg2 | leg3 | leg4 |
|----:|-------------:|-----:|-----:|-----:|--------------:|-----:|-----:|-----:|
| 1   | 49037  | 145376 | 259061 | 533705 | 308  | 696  | 4408 | 2767 |
| 2   | 221625 | 104793 | 239629 | 882329 | 1702 | 476  | 3939 | 3721 |
| 3   | 234620 | 217069 | 444320 | 896473 | 1835 | 984  | 7673 | 4898 |

The patched-conic seed is **10⁴–10⁶ km and 10²–10³ m/s per leg away from
n-body continuity** — the idealized→real gap this lane exists to measure,
quantified. (Most of the position gap is the periapsis-node vs Lambert-leg
geometry mismatch — §5 — compounded by continuous moon gravity along the arc.)

## 4. BALLISTIC CLOSURE ATTEMPT (cycle 1, `shoot_cycle`, 60 nfev, 1194 s)

| leg | seed |Δr| (km) | final |Δr| (km) | final |Δv| (m/s) | node→moon dist (km) |
|----:|---------------:|----------------:|-----------------:|--------------------:|
| 1   | 49037  | 2.770   | 216.284 | 29993 |
| 2   | 145376 | 1.035   | 135.166 | 16463 |
| 3   | 259061 | 10.900  | 148.356 | 2831  |
| 4   | 533705 | 1.545   | 0.006   | 33585 |

The multiple-shooter pulls the seed to **near-continuity in position**
(10⁴–10⁶ km → **1–11 km** on 10⁶ km legs, 5 orders of magnitude) with all
nodes inside their moons' SOI bands, **but velocity continuity stalls at
~135–216 m/s on three of four legs** within the 60-nfev budget (converged
flag = False; threshold <1 km AND <0.01 m/s). It is settling onto a *powered*
solution (~100–200 m/s/leg residual ΔV), not a ballistic one. (A
deeper-budget / homotopy-continuation shoot might drive the velocity residual
down further, but that is exactly the work the published anchors would let us
target and that print does not support — see the verdict.)

## 5. ONE SEED-CONSTRUCTION CAVEAT (not a bug)

The large *raw* seed defect (§3) is dominated by `periapsis_node` placing each
node at the moon-centered hyperbola **periapsis** (deep in the well, periapsis
*speed* `sqrt(V∞² + 2μ/r_p)`) while the conic Lambert legs connect moon
**centers** at the slower Jupiter-frame leg velocity. The two are geometrically
inconsistent by construction, so a propagation from a periapsis node does not
aim at the next periapsis — hence 10⁵-km seed misses even though the *conic*
chain is near-ballistic. This is the documented patched-conic→n-body bridge,
not an error: `shoot_cycle` is given full freedom on nodes 1–4 and absorbs the
mismatch (§4, position closes to km-level). It is recorded so a future deeper
attempt can consider a center-node seed (node at the moon center, Jupiter-frame
leg velocity) which would start far closer to n-body continuity. The module's
equations (Lambert legs, defect/turn formulas, periapsis hyperbola, energy
conservation) all checked out — the energy guard test confirms the propagator
is near-Keplerian (rel. ΔE < 1e-4 over 5 d) far from the moons.

## 6. VERDICT

**The same-family candidate reproduces the Member D family signature at the
patched-conic level on real JUP365 geometry (sequence, ToF band 3/3, all
flybys > 100 km, near-ballistic flyby defects below the published 1.04e-7 m/s)
— but does NOT close ballistically in the restricted n-body model from this
seed and budget: the n-body shoot settles onto a powered ~100–200 m/s/leg
solution.**

Per the ceiling assessment (§0): the conic-level match IS the V1-class
same-family existence evidence the row's print permits — an independent
construction re-derives the family the row describes, in the stated model, at
the stated epoch, with the stated qualitative signature. The n-body
non-closure is a clean, honest finding: it measures the idealized→real
continuous-gravity gap (large; §3) and shows our local corrector reaches
geometric (position) continuity but not ballistic (velocity) continuity within
budget. **Neither outcome supports V3** — that would require THEIR per-flyby
encounters confirmed on the real ephemeris, which Member D does not publish.

## 7. PROPOSED WRITEBACK — WITHHELD

**No V0 → V1 writeback is proposed for `liang-2024-cgcec-ephemeris-2033`.**

The conic-level signature match is genuine same-family evidence, but it rests
on the *idealized* scaffold's family (Member A leg ToFs) phased to real
geometry, and the n-body re-propagation — the step that would turn "same
family in the conic model" into "this trajectory exists in a continuous-gravity
model" — did NOT close ballistically. Promoting V0 → V1 on a candidate whose
only clean success is at the same patched-conic fidelity as the already-V1
idealized Members A/B/C (this note's §0 ceiling, and the A/B/C reproduction
note) would add no fidelity beyond what those rows already carry, while the
one *new* test this lane ran (real-ephemeris n-body closure) returned a
negative. The honest record is: family plausible at conic fidelity on real
geometry; real-ephemeris ballistic closure not demonstrated. That is below the
bar for changing the row's validation level. The row's `data_gaps` (no
per-flyby anchors) remain the binding constraint, exactly as §0 anticipated.

Re-open only with: (a) the authors' per-flyby epoch/V∞ data (would enable a
targeted, anchored closure and a real V3 test), or (b) a deeper n-body
campaign — center-node seed (§5) + continuation in moon-gravity strength +
larger shoot budget — that drives a self-constructed candidate to genuine
ballistic closure (which would still cap at V1-class same-family, never V3).
