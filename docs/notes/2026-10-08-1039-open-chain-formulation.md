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
