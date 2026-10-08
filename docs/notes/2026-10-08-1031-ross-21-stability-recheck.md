# #1031: ross-rt-em-cycler-21-2025 stability re-check (the #1030 flag is retracted)

Task `#1031` (earthmoon-opus, 2026-10-08), from `#1030` (`d5fa6f1c`). Script
`scripts/run_1031_ross21_recheck.py`, output `data/1031_ross21/scan.json`. Optional patch
`data/1031_ross21/1031_state_precision.patch` (ydot0 at full precision). No catalogue edit.

## 1. Verdict

**The row's claim holds. The stored state is linearly stable in the plane, with Barden nu = +0.050
(the stored value is +0.05007). Its stable window exists in C from 3.129389531074196 to the fold at
3.12938953109233, and that fold equals the paper's C^max to 5e-15.**

The `#1030` flag (nu = -1.35 at "the stored C") was my error, and I retract it.
- The audit re-closed the row at the Jacobi constant of the ROUNDED stored state,
  3.129389531068003. It should have used the stored C^stable, 3.129389531088256.
- The difference is 2.0e-11. That is more than the whole stable window is wide in C (1.8e-11).
- The re-closed orbit was therefore a real family member, but one just outside the window.
- I fixed the audit script to use the stored C and re-ran it. The row now gives in-plane k_par =
  +0.048 (nu = +0.024) and k_perp = -1.140. Both are stable, and the row is consistent.

## 2. Method

- The paper (Ross & Roberts-Tsoukkas 2025, AAS 25-621, Table 3; mining note
  `2026-06-11-ross-roberts-tsoukkas-2025-mining.md` sec. 4) gives for (2,1):
  - C^stable = 3.129389531088256;
  - C^max = 3.129389531092325;
  - T^stable = 19.44043166795154 TU;
  - perilune width of the stable subfamily Delta_p_m = 4.23 km.
  So C^max - C^stable = 4.1e-12: the window lies against a fold of the family in C, where a
  fixed-C corrector is ill-conditioned.
- I parametrised the family by x0 (the start crossing, x0 > 0) and solved ydot0 by Newton on
  xdot at the 4th y = 0 crossing (t = 9.72, the half period). The fixed-x0 corrector comes from
  `scripts/run_997_lineage.py`.
- I scanned x0 in steps of 5e-7 over +-2e-5 around the row's x0. At each member: C, T, Barden
  nu = k_par / 2 and k_perp from the full-period monodromy, and the perilune distance.
- I bisected in x0 for nu = 0 and nu = +-1.
- Model: the row's mu = 0.012150584270572 (Ross p.3).

## 3. Results

| point | x0 | C | T (TU) | Barden nu | k_perp | perilune (km from centre) |
|---|---|---|---|---|---|---|
| nu = -1 edge (period doubling side) | 0.723732775152 | 3.129389531074196 | 19.4401829751 | -1.003 | -1.1398 | 26,069.5782 |
| nu = 0 midpoint | 0.723733540777 | 3.129389531087740 | 19.4401616795 | -0.006 | -1.1398 | 26,069.5759 |
| row state (x0 as stored, ydot0 re-solved) | 0.7237335857 | 3.1293895310882553 | 19.4401604300 | +0.050 | -1.1398 | 26,069.5758 |
| nu = +1 edge = fold in C | 0.723734308355 | 3.129389531092330 | 19.4401403299 | +0.991 | -1.1399 | 26,069.5736 |

- **|nu| < 1 holds for x0 in [0.723732775, 0.723734308]** (width 1.53e-6), that is, for C from
  3.129389531074196 (nu = -1) up to the fold in C (nu = +1).
- **k_perp is about -1.14 throughout**, so the window is vertically stable as well.
- Outside the window, nu runs linearly at about 1.3 per 1e-6 in x0: from -26 to +26 over the
  +-2e-5 scan.
- Against the paper:
  - C at the nu = +1 edge is 3.129389531092330 against C^max = 3.129389531092325 (difference
    5e-15). The fold is the tangent bifurcation, as it should be.
  - C at our nu = 0 point is 3.129389531087740 against C^stable = 3.129389531088256 (difference
    5e-13).
  - T at our nu = 0 point is 19.44016168 against T^stable = 19.44043167 (difference -2.7e-4 TU).
    This is the same offset the `#212b` adoption note recorded, and the row's own data_gaps
    entry explains it.
  - **Perilune width: unresolved.** Our perilune distance changes by only 0.0046 km across the
    window, against the paper's Delta_p_m = 4.23 km, about 900 times larger. Either the paper
    measures a different quantity (another pass, or periapsis of an osculating orbit), or its
    window is defined differently. I have not resolved this. It does not affect the
    stable/unstable verdict.
- **The row's V2 claim "Barden |nu| < 1 in the published stable window" holds**, at C^stable
  itself (nu = +0.05 at the stored x0).

## 4. Proposed patch (optional)

`data/1031_ross21/1031_state_precision.patch`:
- Replaces the stored ydot0, 0.4137707374 (10 digits), with 0.41377073737552744, re-solved at the
  stored x0. The new state's C is 3.1293895310882553, equal to C^stable to 1e-15. The comment
  explains why.
- With the 10-digit value, any check that re-closes the row at its own state's C lands outside
  the window. That is what misled `#1030`.
- x0 is unchanged, and the tests that pin x0 = 0.7237335857 (`test_cr3bp_ross_families.py`,
  `test_binary_star_search.py`) are unaffected.
- The patch applies with `git apply --check` and passes the schema check.

## 5. Changes to `#1030`

- `scripts/run_1030_stability_audit.py` now re-closes at the row's stored `jacobi_constant`.
- I re-ran it; `data/1030_stability_audit/audit.json` is regenerated. ross-21 is now k_par =
  +0.048, k_perp = -1.140, stable. ross-31 is unchanged (nu = 0.01545 reproduced).
- The `#1030` note's sec. 3 item 2 is corrected to point here.

## 6. Applied (2026-10-08)

The precision patch was applied in the lead's ratchet window of 2026-10-08, in the same commit as the `#970` row and the `#1030` wording.
