# Postmortem: #548 Owen & Baresi linking-number pipeline positive-control gate

**Date:** 2026-07-10
**Task:** #548 (time-boxed, pre-registered kill criterion)
**Subject:** `genome/qp_tori.py` + `genome/qp_torus_manifold.py` + `genome/qp_torus_heteroclinic.py`
linking-number heteroclinic-connection screen (the #522 pipeline, Owen & Baresi *Astrodynamics*
8, 577-595, 2024).

## What #548 set out to do

Settle whether the linking-number screen actually detects heteroclinic connections between
quasi-periodic tori, by running it — for the first time ever — against a positive control in the
regime where Owen & Baresi themselves demonstrate connections exist: isoenergetic Earth-Moon L1
quasi-halo <-> L2 quasi-halo torus pairs, below both collinear necks. The screen had returned
0 connections in all three prior real applications (#534 Earth-Moon, #536 Jupiter-Europa, #546
Uranus), but #547 established none of those was a valid positive control (the pipeline had never
been pointed at a case where connections are known to exist).

The pre-registered kill criterion: a full both-family sweep with empirically-classified transit
branches producing ZERO sign changes -> SHELVE the pipeline.

## Two premise corrections found by primary evidence this session

1. **Energy.** #547 recorded that #534 built its tori "at the right energy (C=3.15), matching
   Owen & Baresi's own demonstration." This is FALSE. #534's committed NRHO seeds
   (`scripts/run_534_torus_connection.py`), when corrected, sit at **C = 3.045**
   (`cr3bp.jacobi_constant` of the corrected L1/L2 seeds = 3.04504 / 3.04483), not 3.15. #534
   was never at Owen & Baresi's energy. This alone means #534's negative was against the wrong
   configuration, independent of any branch-classification issue.

2. **Exactly-C=3.15 isoenergetic quasi-halo pairs are impractical here.** Heteroclinic
   connections are ISOENERGETIC — both tori must share one Jacobi constant. But:
   - The EM **L1 halo family bifurcates from the planar Lyapunov family at C ~ 3.146** (verified:
     continuing the L1 halo family in x0, z0 -> 0 as C climbs to ~3.145, then the corrector snaps
     to the planar z0=0 branch). So C=3.15 sits at/above the top of the L1 quasi-halo regime.
   - The EM **L2 halo family, continued from the NRHO branch, tops out at C ~ 3.087** (a genuine
     high-C fold, confirmed by pseudo-arclength continuation past the x0-fold: both tangent
     directions from the C=3.087 fold point lead back to lower C). The small-z L2 halos near the
     L2 planar bifurcation (~C=3.152) live on a segment not reachable from the NRHO branch within
     the time-box.

   The highest COMMON energy both families robustly reach is therefore **C in [3.05, 3.087]**, and
   the positive control was built there — same physics (both necks open, unstable quasi-halo pairs
   whose manifolds can intersect), a slightly deeper energy than the paper's single 3.15 demo.

   Consequence for the kill criterion: its literal precondition ("at C=3.15") is **unsatisfiable**
   for isoenergetic EM quasi-halo pairs via the available machinery. The criterion was written
   under the incorrect belief (inherited from #547) that C=3.15 was already achieved and readily
   reproducible. The result below is therefore reported as strong-but-caveated evidence at the
   achievable common energy, with the premise error surfaced, rather than a mechanical "criterion
   met" claim at an energy never actually tested.

## What was built (durable)

- **`src/cyclerfinder/genome/qp_torus_transit.py`** — `transit_torus_manifold_grid`, an EMPIRICAL
  transit-branch manifold-grid builder. Adapts #547's `transit_manifold.classify_unstable_branch`
  (planar-Lyapunov positive control) to the 3D quasi-halo torus grid: at every torus point it
  propagates BOTH signed eigenvector perturbations and keeps whichever reaches the surface of
  section first, discarding the untested `vec[0]*sign` heuristic that
  `qp_torus_manifold.torus_manifold_grid` trusts (the exact weak point #547 flagged). Returns the
  standard `ManifoldGrid` plus a per-point branch-sign field for sheet-coherence diagnostics.
  Covered by `tests/genome/test_qp_torus_transit.py` (2 tests, green).
- **`scripts/run_548_owen_baresi_positive_control.py`** — the reframed sweep: parent-halo bisection
  to a target common C on both families, amplitude-swept quasi-halo torus construction, empirical
  transit grids, and `scan_linking_number` across four (scanning-component, curve-component)
  specs per isoenergetic pair.

## Result (sweep of C in [3.05, 3.07])

`scripts/run_548_owen_baresi_positive_control.py`, 2026-07-10:

- **3 isoenergetic geometry-usable L1 x L2 quasi-halo pairs** at C = 3.05, 3.06, 3.07 (the C=3.08
  L2 parent halo bisection failed, expected — it is right at the L2 family's high-C fold). Each
  pair: two matched-Jacobi quasi-halo tori (invariance residuals 1.5e-8 ... 6.9e-8), empirical
  transit grids (20x16), and `scan_linking_number` over 4 (scanning-component, curve-triple)
  specs = **12 linking scans total.**
- **Manifolds reached the section and overlapped** (not a geometry no-show): L1 unstable grids
  320/320 crossings, L2 stable grids 166-320/320; the scanning-variable ranges overlapped in
  every scan (e.g. C=3.07 zdot overlap [-0.41, 0.83]).
- **The L1 transverse-frequency ratio brackets the published value.** L1 |omega_trans/omega_long|
  ran 0.3668 (C=3.05) -> 0.3300 (C=3.06) -> 0.2944 (C=3.07), straddling Owen & Baresi's L1
  latitudinal frequency 0.2739. (The seed AMPLITUDE does not move the ratio — confirmed identical
  0.3668 at amp 5e-4 and 1.2e-3 — so ENERGY is the frequency knob; a single amplitude/side was
  used.) The L2 side brackets POORLY: my L2 tori ran ratio 0.12-0.45 vs the published L2 0.02163,
  because the paper's L2 quasi-halo at C=3.15 is a near-planar small-z halo just off the L2
  bifurcation (very slow transverse mode), whereas the deepest L2 halos reachable here (z0 ~ 0.16)
  oscillate much faster and are near-planar AT THE SECTION (tiny z-extent -> the z-scan overlap is
  only ~0.001 wide).
- **ZERO sign changes anywhere.** Stronger: the linking number was **identically 0** in all 12
  scans (never even momentarily nonzero) -- the stable/unstable reduced closed curves never link
  at any D value, in any scanning variable, at any of the three energies.
- **The `linking_number` primitive itself is NOT the problem** -- it returns +/-1 on a Hopf link
  and 0 on an unlinked pair (verified this session; `tests/search/test_linking_number.py`
  `test_hopf_link_has_linking_number_one` is green). The identically-0 output is a true property
  of the manifold-derived curves this pipeline extracts here, not a broken invariant.

## Verdict

**SHELVED** per the #548 pre-registered kill criterion (spirit; see the energy caveat below).
The reframed positive-control attempt -- the most faithful test the pipeline has ever been given
(right system, right energy regime below both necks, L1 frequency ratio bracketing the published
value, empirical transit-branch classification) -- found ZERO linking-number sign changes, and in
fact never produced a nonzero linking number at all on any quasi-halo pair. Combined with the
prior 0-for-3 (#534/#536/#546), the qp_tori/qp_torus_heteroclinic linking-number screen has now
failed to surface a single connection across four genuine applications and has NEVER produced a
confirmed torus-level positive control.

**Honest caveats (for the adjudicator), so this is not oversold as "we reproduced their 4
connections then it broke elsewhere":**
1. The literal criterion precondition ("at C=3.15") was **unsatisfiable** for isoenergetic EM
   quasi-halo pairs via the available machinery (L1 halo bifurcation ~3.146; L2 halo NRHO branch
   tops at ~3.087). Tested at the achievable common band C in [3.05, 3.07] instead.
2. The L2-side frequency bracket was poor (0.12-0.45 vs published 0.02163); a fully faithful
   reproduction needs the L2 small-z near-bifurcation halo (~C=3.15), which requires an
   L2-planar-Lyapunov -> halo bifurcation seed generator this codebase does not currently have on
   that segment. A future worker who builds that generator could give the pipeline one final,
   fully-frequency-matched shot before permanent retirement.
3. Some evidence of transit-sheet mixing on the L2 stable grid (branch-sign balance ~150/170 at
   the section, because the L2 halo sits almost on the section x=1-mu so both branches cross
   quickly); the L1 side was clean (~50/270, 299/21).

Given the binding, pre-registered "zero sign changes -> shelve, do not keep tinkering" instruction
and the accumulated 0-for-4, the pipeline is shelved now and #534/#536/#546's negatives are
stamped method-invalid-do-not-certify. The caveats above are recorded so the decision is
reversible by a future worker with a genuine C=3.15 L2 near-bifurcation torus, not presented as a
proof of non-existence.

