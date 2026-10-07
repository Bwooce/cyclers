# #1039: how to pose an open chain in continuous gravity (the end-pin artefact)

Status: PRE-REGISTRATION of formulation (b) (written and committed before any #1039 run).
Time box about 4 hours of calls, from 2026-10-08 10:20 AEDT. No catalogue writes.

## 1. Problem

Three results stop where the chain ends are pinned to the patched-conic V_inf VECTORS (#968
amendment 12):
- #968 rung (b): Europa impact at sigma 0.62.
- #1043: a fold at sigma 0.505, with the end Europa nodes as the growing defect.
- Any real-ephemeris standing for the gc rows.
The inference is that the pinned ends absorb the full-mass shift. #1039 tests that.

## 2. Formulation (b): direction-only end pins, free magnitudes, a period-like condition

The 6 end rows of amendment 12 are replaced by 6 rows that leave the V_inf magnitudes free:
1. The first node's inbound V_inf DIRECTION equals the patched-conic chain's (2 rows: components along
   two unit vectors normal to the target direction).
2. The last node's outbound V_inf DIRECTION likewise (2 rows).
3. Magnitude match: |V_inf in, first node| = |V_inf out, last node| (1 row; the paper's own junction
   condition, as in the date corrector's wrap).
4. Period-like: t_last - t_first equals the patched-conic chain's span (1 row).

Everything else is as #1043: moon-relative nodes, gauges, forward-backward legs, rtol 1e-13, sigma
continuation from 0.02 (factor 1.2, halving to 1e-5, step cap), the floors, the acceptance checks
of the object's own task. Jacobian rows for the end conditions are by central differences (no
propagation is involved).

## 3. The positive control comes first

- Control: the one-cycle EGGIE (#1043 seed, `data/1043_eggie/pc_chain_n1.json`, the paper's ideal
  model).
- PASS: it reaches sigma = 1 at the lane floors with the #1043 acceptance (25 km floor, no
  unscheduled Hill-radius pass), and the identity holds at sigma = 0.02 (V_inf within 0.02 km/s of
  the chain).
- Only if the control passes is formulation (b) tried on jup365 GanEur#316 (rung (b), #968
  amendment 13's single cycle).
- If the control folds again near sigma = 0.5, the artefact is not (only) the ends, and that is
  reported.

Expected outcome for the control under (b): reaches sigma = 1 with probability about 0.5. The
magnitudes are freed, but the directions stay pinned, and #1043's growing defect was the Europa end
nodes' periapsis, which also depends on direction.

## 4. Formulation (a') (pre-registered later, only if (b) fails or after it)

The N = 3 chain with Levenberg-Marquardt, or an end-pin homotopy (pins relaxed continuously from
the patched-conic values). Its exact definition is written before it runs.

## 5. Results

(pending)

## 5. Results

(Time box start: 2026-10-08 10:05 AEDT, not 10:20 as written in the header.)

### 5.1 Formulation (b) on the positive control (one-cycle EGGIE, ideal model): PASS

`data/1039_b_eggie/sigma_n1.json`.
- Identity at sigma = 0.02: V_inf 9.0669 / 7.0808 / 7.0807 / 8.1985 / 9.0669 km/s against the chain's
  9.068 / 7.082 / 7.082 / 8.207 / 9.068 (max difference 0.0085 <= 0.02).
- The continuation converged at the lane floors at every step from 0.02 to sigma = 1, in 23 points:
  no failed step, no noise-floor-limited point, no fold.
- At sigma = 1:

  | Node | V_inf (km/s) | r_p (km) | altitude (km) |
  |---|---|---|---|
  | E0 | 9.032 | 2,643 | 1,082 |
  | G1 | 7.046 | 3,659 | 1,028 |
  | G2 | 7.041 | 4,396 | 1,764 |
  | I | 7.880 | 9,547 | 7,726 |
  | E1 | 9.032 | 3,225 | 1,664 |

  All above the 25 km floor; no unscheduled pass inside any Hill radius.
- Beyond the control criterion, with #1043's own criteria:
  - IAS15 re-fly of every half-arc at sigma = 1: max 1.0e-6 km, PASS.
  - Identity going down: sigma = 0.01 max |dV_inf| 0.0044, sigma = 0.005 0.0023 km/s (converged at the
    floors), PASS.
- Reading: with the end MAGNITUDES free (directions pinned, the paper's magnitude match at the ends,
  the span held), the fold of #1043 at sigma = 0.505 disappears. The #1043 artefact WAS the vector end
  pins.
- Consequence for #1043: under formulation (b) the published one-cycle EGGIE CLOSES in continuous
  gravity (ideal model), meeting every #1043 criterion: floors, the 25 km floor, no unscheduled pass,
  identity at sigma <= 0.01 and IAS15. Full-mass V_inf: E 9.03, G 7.04-7.05, I 7.88 km/s. The lead
  records this against #1043.
- Next, as registered: formulation (b) on jup365 GanEur#316, one cycle (#968 rung (b)).
