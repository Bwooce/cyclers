# High-Order Transfer Map fixed-point search — mining note (#231)

**Date:** 2026-06-13
**Paper:** Zhou, Anoè, Armellin, Qiao, Li, "Identifying Fixed Points in the
Three-Body Problem Using a High-Order Transfer Map," arXiv:2509.12671 (v1, 16
Sep 2025).
**Triage:** V0-sourced CR3BP row candidate + independent novelty cross-check +
METHOD-NOTED (alternative discovery engine).

## Paper summary

The paper proposes an **exhaustive** periodic-orbit (fixed-point) search for the
CRTBP built on **high-order transfer maps (HOTMs)** in the **differential
algebra (DA)** framework (DACEyPy). The pipeline:

1. Construct a single-revolution HOTM (Taylor polynomial) that maps a Poincaré
   section (`y=0`, `ydot>=0`, `0<x<1-mu`) back to itself. The Jacobi constant and
   the section constraint each remove one DOF, so the planar map has only
   `(dx0, dxdot0)` inputs and the spatial map `(dx0, dxdot0, dz0, dzdot0)`.
2. **Automatic Domain Splitting (ADS)** with feasibility detection (LDB bounder)
   keeps each polynomial accurate while pruning infeasible subdomains.
3. A **two-stage polynomial optimization** (RPO / convex SCOP via MOSEK): stage 1
   finds combinable subdomain *sequences*; stage 2 refines them to fixed points.
   A single-rev HOTM is **chained** to reach multi-revolution orbits without
   re-propagating — the source of "arbitrary revolution count" scalability.

Applied to **Earth-Moon CRTBP**: planar POs to **9 revolutions**, spatial to
**4**. Recovers known DROs and Lyapunov orbits, **plus a previously undocumented
hybrid family** (`Png'`) with mixed DRO/Lyapunov character.

Paper's Earth-Moon constants (Table 1): `mu = 0.012150584269940354`,
`mu_E = 398600.43543609598`, `mu_M = 4902.8000661637961` km^3/s^2.
First planar case Jacobi constant **C_J = 3.00022**.

## The hybrid family (Q1): `Png'`

The undocumented family is denoted **`Png'`** — hybrid POs, "at least one segment
resembling a Lyapunov-type trajectory and another resembling a DRO-type segment."
The paper speculates (from Hénon's g'-family bifurcation picture) that `Png'`
may exist **only for odd n**, and frames these as dynamical *bridges* between
DRO and Lyapunov orbits — i.e. candidate mission transition structures. All
listed fixed points were post-refined by differential correction (half-period
`ydot -> ~1e-14`).

**Transcribed candidate IC — lowest-revolution member, P5g' (n=5), Table 3,
first planar case (C_J = 3.00022), planar state `[x0, 0, 0, xdot0, ydot0, 0]`:**

| field | value (nondimensional) |
|---|---|
| x0   | `0.807357887647950` |
| xdot0| `-0.0956081545978604` |
| ydot0| `0.433518861583397` |
| Period | `11.1751086919436` TU |

These are **CANDIDATE** values (transcribed from the printed table — reproduce
before any wiring). Note the IC is **not a perpendicular crossing** (`xdot0 != 0`),
so `Png'` members are *general* fixed points, not the symmetric perpendicular
class; `correct_symmetric_fixed_jacobi` / `barden_stability` do **not** apply.

Other `Png'` members printed (Table 3): P7g'-I (x0=0.916929181700578,
xdot0=-0.175717632213615, ydot0=0.528042521610021, T=11.6960585356669),
P7g'-II (x0=0.843301779064331, xdot0=0.0264518382948230, ydot0=0.433467279267690,
T=11.6716427018773), P7g' (x0=0.800533327625822, xdot0=-0.0783668113509505,
ydot0=0.441933007711387, T=15.0692523706070), P9g' (x0=0.807337935300132,
xdot0=-0.0956506138795539, ydot0=0.433522872600990, T=20.9914771396290). A
P3g' was found only in the *second* planar case (wider domain, Table 5:
x0=0.852098052983502, xdot0=-0.187721536949396, ydot0=0.368541113320107,
T=9.25206400352203, C_J=3.020052).

### Novelty cross-check vs our #218/#219 continuation campaign

`Png'` is a real **independent-group** family discovery. Was it in scope for our
campaign that found 0 distinct NEW families from sourced seeds?

- **Out of our scope as run.** Our continuation seeds from sourced *symmetric*
  perpendicular-crossing ICs and continues single families. `Png'` is (a)
  **non-symmetric** (xdot0 != 0) and (b) a **multi-revolution bifurcation** of the
  DRO g'-family reached by a global polynomial-optimization sweep over a 2-D
  Poincaré domain, not by single-family continuation. Our campaign could not have
  surfaced it because it never enumerates multi-rev fixed points over a domain;
  it refines/continues from a seed it is already handed.
- This is **consistent** with #218/#219's honest "0 new families" conclusion: the
  gap is a **method-capability** gap (global multi-rev enumeration), not a defect
  in our corrector. It is exactly the kind of region the negative-results registry
  flags as "empty *conditional on the method*."

## JPL-oracle triangulation (Q2)

Live `cyclerfinder.verify.jpl_periodic_orbits.query("earth-moon", ...)`.
**mu reconciliation:** JPL `mu = 0.01215058560962404` vs ours
`0.01215058439469525` -> **rel diff 1.0e-7** (the irreducible re-closure floor;
paper's mu is yet a third value, 1.03e-8 from ours).

| JPL family | members | nearest to C=3.00022 | nearest to T=11.175 |
|---|---|---|---|
| dro | 10998 | C=3.000114, T=1.576, x0=0.88494, vy0=0.47062 | C=1.541, T=6.305 (max) |
| lyapunov L1 | 3108 | C=3.000292, T=4.329, x0=0.76895, vy0=0.48103 | C=2.742, T=7.446 (max) |
| lyapunov L2 | 4298 | C=3.000187, T=4.539, x0=1.02279, vy0=0.80651 | C=2.873, T=8.214 (max) |

Verdict:
- **Base families confirmed in JPL.** The JPL DRO at C~3.00022 (x0=0.88494,
  vy0=0.47062, T=1.576) matches the paper's n=1 DRO (Table 3: x0=0.88500968,
  vy0=0.470630, T=1.5745) to ~1e-4 — exactly the mu/quantization gap expected.
  The JPL L1 Lyapunov at C~3.00029 (x0=0.768948, vy0=0.481028) matches the
  paper's listed Lyapunov IC (x0=0.7688974950452078, vy0=0.4811655138591737)
  that lay *outside* its search domain. Independent triangulation: **the paper's
  base orbits are real.**
- **The hybrid P5g' is NOT in the JPL catalog.** JPL's DRO family tops out at
  period ~6.3 TU and Lyapunov families at ~7.4/8.2 TU; **no member anywhere near
  T=11.175 TU exists**. JPL enumerates single-period families, not the
  multi-revolution g'-bifurcation orbits. The "previously undocumented" claim
  survives this independent triangulation — `Png'` is genuinely outside the
  standard JPL family parametrization.

## Our-corrector closure (Q3)

P5g' IC run through `search.cr3bp_periodic.correct_periodic` at **our** mu
(`cr3bp_system("Earth","Moon")`):

- C(s0) at our mu = **3.0002200008** (target 3.00022 — matches; the IC is
  energy-consistent under our mu despite the 1e-8 mu offset).
- **converged = True**, closure residual = **3.45e-12** (full-state periodicity).
- recovered **T = 11.17510869199** (paper 11.1751086919436 — agrees to ~11 sig
  figs), recovered **C = 3.0002200012**.
- recovered state0 = `[0.807357885, 2.3e-9, ~0, -0.095608147, 0.433518864, ~0]`
  (drift from the transcribed IC < 3e-9 — clean reproduction).
- Stability (full-period monodromy; Barden symmetric form N/A for this
  non-symmetric orbit): eigenvalue moduli
  `[2.78e-4, 0.436, ~1, ~1, 2.29, 3.60e+3]`; max |lambda| ~ **3598**,
  nu = 0.5(lambda+1/lambda) ~ **-1799** -> **strongly unstable** (as expected for
  a multi-rev bridge orbit; would need manifold/station-keeping to exploit).

**Closure verdict: clean same-model closure.** Our corrector independently
reproduces P5g' at our mu to ~1e-12 residual and ~11 sig figs in period. This is
a **same-model V1-class candidate** (sourced IC + independently re-closed by our
own corrector), but **HELD** (see follow-on) — it is a CR3BP periodic orbit, not
itself a cycler trajectory.

## METHOD-NOTED (Q4): DA / HOTM as a discovery engine

The HOTM+ADS+RPO pipeline is a genuinely different discovery engine from ours:
a **single precomputed single-rev polynomial map** chained for arbitrary
revolution counts, with **exhaustive** (domain-covering, not seed-local) fixed-
point enumeration over a Jacobi-constant slice. This is precisely the
multi-revolution global-enumeration capability our continuation campaign lacks
and that let an independent group surface `Png'`. **Relevance:** if the project
ever wants to *discover* (not just continue) new multi-rev CR3BP families in a
target energy band, a DA/HOTM sweep (DACEyPy + a convex solver) is the
method-of-record to consider. **No implementation expected or done here** — noted
for the negative-results registry as a capability that would subsume our current
sourced-seed continuation sweep.

## PROPOSED FOLLOW-ON

**Recorded outcome: clean closure + base-family cross-confirm, but HELD — no
row ingest.**

- P5g' closes cleanly under our corrector (res 3.45e-12, T to 11 sig figs) and
  the paper's *base* DRO/Lyapunov orbits independently triangulate against JPL.
  The transcription and our model are sound.
- **Hold rationale (not a defect):** `Png'` is a CR3BP **periodic orbit**, not a
  cycler trajectory (no Earth-Earth resonant transfer structure; it is a
  DRO/Lyapunov-bridge in the Earth-Moon rotating frame). Per the catalogue
  discipline (JPL-oracle note; validation-ceiling memory), non-cycler CR3BP
  families do not feed the cycler catalogue. The cross-check value is the
  **method/novelty triangulation**, not a new catalogue row.
- **Negative-results registry entry (proposed):** record that exhaustive
  multi-rev fixed-point enumeration over the C_J~3.00 Earth-Moon DRO band yields
  `Png'`-type bridge families that our **single-family continuation method cannot
  reach** — i.e. our "0 new families" sweep is empty *conditional on the
  continuation method*, and a DA/HOTM sweep would be the capability that subsumes
  it. Re-open this region only if such a sweep is implemented.
- **If the project later adopts a DA/HOTM discovery lane:** the five `Png'` ICs
  above (all C_J=3.00022 except P3g' at 3.020052) are ready same-model
  cross-check targets — each should re-close under `correct_periodic` exactly as
  P5g' did.

## Provenance

- arXiv:2509.12671 (Table 1 constants; Table 3 first-planar fixed points;
  Table 5 second-planar P3g'). Cite the arXiv id only.
- JPL SSD periodic_orbits API (live), via
  `src/cyclerfinder/verify/jpl_periodic_orbits.py`.
- Closure via `src/cyclerfinder/search/cr3bp_periodic.py::correct_periodic` at
  `cr3bp_system("Earth","Moon")` mu.
