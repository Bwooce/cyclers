# #1023: re-run of the Jovian n-body results voided by the translation-only wrap

Status: (a) done (consistent synodic period, additive; impact list below); (b) EGGIE pre-registration
in sec. 3; (c) the real-ephemeris items in sec. 4. No catalogue writes; nothing is called novel.

The VOID list is in `docs/notes/2026-10-07-968-jovian-nbody-positive-control.md` sec. 2.1. Lead
rulings 2026-10-08: (a) fix the ideal model's synodic period additively, with an impact list; (b) go
on the EGGIE plan; (c) the real-ephemeris items are NOT RE-RUNNABLE until #1039, with retraction
lines.

## 1. Second defect: the ideal Galilean model's synodic period

`resonant_conic.ideal_moon_smas` builds the Hernandez et al. 2017 (AAS 17-608, p.3) ideal model.
Its smas are chosen so that in one synodic period every moon advances 2 pi k + Delta, with
Delta = 5.2 deg (Io 8 pi + Delta, Europa 4 pi + Delta, Ganymede 2 pi + Delta). The consistent
synodic period is therefore T_syn = (2 pi + Delta) / n_G = 7.1054 d (registry Io sma, lane mu).
Over it every moon advances exactly 5.2 deg, so the configuration rotates rigidly.
`ideal_t_syn()` returns the ideal GANYMEDE period, 7.0042 d. Over that period the moons advance
by different angles (Ganymede 0, Europa -5.13, Io -15.38 deg); over four of them by 0, -20.50 and
-61.51 deg. So the model as coded has no exactly periodic EGGIE, and the #968 note sec. 2.1
statement "no exact periodic orbit exists" holds for the CODED model only. In the consistent
model the configuration repeats rigidly (20.8 deg per 4 T_syn).

The paper prints T_syn = 7.05 d (p.2) and a Table 4 EGGIE total ToF of 28.22 d (= 4 x 7.055). This
matches neither 7.004 nor 7.105 d, which suggests the paper's own a_Io differs from the registry Io.
Recorded, not resolved.

Fix, additive (lead ruling): test first, `06355284`
(`tests/search/test_1023_ideal_t_syn_consistent.py`); then `ideal_t_syn_consistent()`, `a333579b`.
`ideal_t_syn` and its callers are UNCHANGED. mypy src tests clean; test_resonant_conic passes.

## 2. IMPACT LIST: what rests on the old 7.0042 d (the lead decides the switch and re-runs, #1040)

Numeric change if a caller switches to `ideal_t_syn_consistent`: T_syn 7.00419 -> 7.10536 d (+1.444 %).
Every resonant sma a = (mu (n_syn T_syn / (2 pi n_rev))^2)^(1/3) grows by +0.961 %.

| Item | Uses | Old -> new | What rests on it |
|---|---|---|---|
| `resonant_conic.eggie_resonant_sma` (4:5) | `ideal_t_syn` | a 909,420.4 -> 918,156.8 km; craft period 0.8 T_syn; cycle 28.017 -> 28.421 d | `eggie_initial_guess` / `eggie_refined_guess` (the #480 Stage-1 resonant-conic seed), `jovian_ideal.build_eggie_*_seed`, `ideal_eggie_shoot` (Stage 2-4), `scripts/eggie_maintenance_480.py` |
| `eige_ballistic.EIGE_RESONANT_SMA_KM` (1:1) | `ideal_t_syn` | 1,055,288.9 -> 1,065,426.6 km; cycle 7.004 -> 7.105 d | the EIGE construction (`docs/notes/2026-06-30-480-eige-ballistic-construction-verdict.md`), `scripts/eige_maintenance_480.py` (`...-eige-realeph-maintenance-verdict.md`), `tests/search/test_eige_ballistic.py` |
| `ll2011_ballistic.T_LAPLACE_S` and `GIPEIPE_SMA_KM` (1:2) | `ideal_t_syn` | 7.004 -> 7.105 d; 664,790.3 -> 671,176.7 km | the #493 LL2011 reproduction (`docs/notes/2026-06-30-493-ll2011-ieg-reproduction-verdict.md`: "period 7.004 d vs sourced 7.055 d, -0.72 %, PASS (<1 %)"; with 7.105 d it is +0.71 %, still inside 1 %), `scripts/ll2011_493_reproduce.py`, `tests/search/test_ll2011_ballistic.py` (the 1 % gate; it would still pass) |
| `tests/search/test_resonant_conic.py` lines 60-64 | `ideal_t_syn` | internal consistency (resonant period = 0.8 T_syn) | passes either way (no sourced value) |
| `eggie_ballistic` (#480 ballistic construction) | `ideal_moon_smas` only; leg ToFs free; seam = |V_inf| magnitude match | no T_syn; but its "periodic" seam is a magnitude match over a model with no rigid repeat | `docs/notes/2026-06-29-480-eggie-ballistic-construction-verdict.md`, `tests/search/test_tour_self_consistency.py` |
| `tests/verify/test_ieg_reproduction_golden.py` | the #480 pipeline | skipped; its skip reason cites the Stage-2/3/4 plateaus (already VOID, #968 2.1) | the skip text |
| catalogue rows | none | `hernandez-2017-jovian-ieg-triple-family`, `lynam-longuski-2011-ieg-single-period`, `lynam-longuski-2011-gipeipe` are V0 / no level; no validation rests on 7.004 d | notes only |
| notes | | `2026-06-27-480-ieg-reproduction-verdict.md`, `2026-06-29-480-eggie-*` (all), `2026-06-30-480-eggie-*`, `2026-06-30-480-eige-*`, `2026-06-30-480-level3-approach-c-verdict.md`, `2026-06-30-493-ll2011-ieg-reproduction-verdict.md` | the ideal-model numbers in each |
| data files | none found | no `data/` file names or contents tie to #480/#493/EIGE/LL2011 (grep) | |

Not affected: #968, #1004, #1034 (R-S model, no `resonant_conic`).

## 3. EGGIE in the consistent ideal model: PRE-REGISTRATION (before any EGGIE run)

Script `scripts/run_1023_eggie.py`; data `data/1023_eggie/`; live log `data/1023_eggie/live/`.

### 3.1 Model

- Io, Europa and Ganymede on circular coplanar orbits at the `ideal_moon_smas` radii, periods from
  Kepler III with the lane mu. The synodic period is `ideal_t_syn_consistent()` = 7.1054 d and the
  cycle is T = 4 T_syn = 28.4214 d. Over T every moon advances 20.8 deg, so the configuration
  repeats RIGIDLY, rotated by 20.8 deg.
- Initial phases (assumption, stated now): Io 0 and Europa 0 at t = 0, Ganymede 90 deg, i.e. the
  Laplace angle lambda_I - 3 lambda_E + 2 lambda_G = 180 deg as in the real system. The formula
  makes that angle constant in time (n_I - 3 n_E + 2 n_G = 0 exactly). The paper's own phases are not
  printed.

### 3.2 Step 1, the patched-conic EGGIE in this model

- The #943 date corrector (`two_working_body.correct_dates`) on the cycle E>G | G>G | G>I | I>E
  (Lambert legs; revolutions 0-2 and both branches enumerated per leg), period T, all three moons
  massive. Seeds: Table 4 dates (0, 1.59, 10.19, 17.53 d) plus phase shifts over one synodic period.
  The paper's ToFs serve only as seeds.
- Accepted roots: max residual < 1e-9 km/s and an independent Kepler re-propagation miss < 1 km.
- The root used is the one whose V_inf is nearest Table 4 (E 9.12, G 7.07, I 8.38 km/s). Its gate
  (demanded vs available turn) is reported at the paper's 25 km floor and at the project floors. A
  root that fails the 25 km gate is still continued (the paper's own EGGIE has 0.70 m/s of flyby
  Delta-V), but it is labelled.

### 3.3 Step 2, continuous gravity (corrector A)

- Forward-backward multiple shooting. Nodes E, G1, G2, I (state and epoch each, periapsis gauge to
  its own moon) and the wrap node E' = (Q x_E, t_E + T), with Q the 20.8 deg rotation. 28 unknowns,
  28 residuals. Moon-relative node coordinates, DOP853 + analytic STM, rtol 1e-13.
- Joint mass scale sigma on all three moons (GMs and softening radii x sigma): from 0.02 up to 1
  (natural continuation, factor 1.2, halving to 1e-5), and down to 0.01 and 0.005 for the identity.
  Intermediate points may be noise-floor-limited (as #968 amendment 11); verdict points may not.

### 3.4 What "closes" means, and why the published object would meet it

In this model the configuration repeats rigidly, so EXACT periodicity in the frame that rotates
20.8 deg per cycle is a valid criterion (lead ruling (b)). CLOSES =
1. at sigma = 1, every match and the wrap meet the lane floors (1e-3 km, 1e-6 km/s), gauges < 1e-9;
2. every node altitude >= 25 km (the paper's floor), no unscheduled pass inside any of the three
   moons' Hill radii;
3. identity: at the smallest sigma <= 0.01 that converges at the floors, every V_inf is within
   0.02 km/s of the step-1 patched-conic root (criterion 4c form); full-mass agreement with Table 4
   is NOT a criterion (#968 sec. 8.2);
4. IAS15 re-fly of every half-arc at sigma = 1: < 1e-2 km, 1e-7 km/s.

Why the published object would meet it: the paper's EGGIE is a 4-synodic-period cycler of its
ideal model with 0.70 m/s of flyby Delta-V ("can probably be optimized to zero"). In a model whose
configuration repeats rigidly, a ballistic continuous-gravity counterpart, if it exists, is exactly
periodic up to the rotation. Its sigma -> 0 limit is the patched-conic generating orbit (step 1).

### 3.5 Outcomes (fixed now)

- CLOSES: criteria 1-4 met. EGGIE exists in continuous gravity in the consistent ideal model.
- DOES NOT CLOSE: a fold (two solutions at one sigma, bracket; #1004 checks 1-2) or an impact
  (periapsis below the scaled radius or the 25 km floor) before sigma = 1. The diagnosed reason is
  stated.
- NUMERICAL STOP: six halvings without convergence and no fold evidence; undecided.
- No step-1 root near Table 4 (all > 0.5 km/s away): EGGIE NOT RE-RUNNABLE in this model; the
  roots found are listed.

Expected outcome: CLOSES with probability about 0.6. EGGIE's V_inf (7-9 km/s) are high, and the
patched conic is better there; a 0.94-type near-limit flyby is not expected at those speeds.

### 3.6 Result of step 1: EGGIE NOT RE-RUNNABLE in the consistent model (`data/1023_eggie/pc_roots.json`)

- The date corrector found 156 distinct exact roots (residual < 1e-9 km/s, Kepler re-fly miss < 1 km)
  of the E>G | G>G | G>I | I>E cycle over all 50 revolution/branch combinations and 12 seed phases.
- The nearest root to Table 4 is 0.652 km/s away (V_inf E 9.772, G 6.497, I 7.858). It demands turns
  of 117.5, 0.5, 129.8 and 163.5 deg, which is physically impossible at these speeds.
- NONE of the 156 roots passes the turn gate, at the paper's 25 km floor or at the project floors.
- By the rule fixed in 3.5 (all roots > 0.5 km/s from Table 4), EGGIE is NOT RE-RUNNABLE in the
  consistent ideal model. The continuous-gravity steps were not run (no generating orbit).
- Conditional on two inputs:
  1. The phase assumption: Laplace angle 180 deg, as in the real system. The paper does not print
     its phases.
  2. The model: the paper's printed T_syn = 7.05 d matches neither the coded 7.004 d nor the
     consistent 7.105 d (sec. 1), so the paper's own ideal model is not pinned down by what it
     prints.
- The #480 history fits this reading. The resonant-conic seed put the V_inf on the Table 4 values
  but was "structurally NON-ballistic" (forward-verify correction note), and the paper's own EGGIE
  needs 0.70 m/s of flyby Delta-V.
- Follow-up suggestions, not run: the same root search (i) in a T_syn = 7.05 d model and (ii) over
  the Laplace angle; for #1040.

## 4. Real-ephemeris items: NOT RE-RUNNABLE until #1039 (lead ruling (c))

| Item | Classification | Reason |
|---|---|---|
| EGGIE level-3 jup365 (`2026-06-30-480-eggie-level3-nbody.md` and its CORRECTION; `...-eggie-realeph-stm.md`) | NOT RE-RUNNABLE until #1039 | `jovian_shoot` closes a real-ephemeris cycle with a periodicity wrap; the real configuration does not repeat, so even the rotation-fixed wrap is only approximately satisfiable; the open-chain formulation is unsolved (#968 sec. 9) |
| #318 CGCEC Sobol smoke n-body stage | NOT RE-RUNNABLE until #1039 | same |
| #501 broadened joint-search n-body stage (6 sequences) | NOT RE-RUNNABLE until #1039 | same |

Registry: seven APPEND-ONLY re-stamps in `data/empty_regions.jsonl` (lines 120-126), one per void
stamp, region_id `<old>-nbody-retracted-1023`. Each keeps the patched-conic prefilter result, drops
the `n-body-shoot` / `analytic-stm` capability tags and the jovian_shoot prune gate, records the
retracted n_shot / n_close, carries the method sha, and says in `verdict` why the n-body stage is
retracted (translation-only wrap plus a non-repeating real configuration). Script
`scripts/run_1023_retractions.py` (idempotent). The original seven lines are untouched.
Checks: tests/data/test_empty_regions.py, test_empty_regions_da_hotm.py, test_method_capability.py,
the registry-reading script tests and tests/scripts/test_scripts_call_preflight.py pass.

## 5. Summary per VOID item

| Item | Outcome |
|---|---|
| EGGIE Stage 2 / 3 / 4 ideal-model plateaus | NOT RE-RUNNABLE in the consistent ideal model: no patched-conic EGGIE near Table 4, none gate-passing (3.6); in the coded model no periodic orbit exists (sec. 1) |
| EGGIE level-3 real-ephemeris | NOT RE-RUNNABLE until #1039 (sec. 4) |
| #318 n-body stage (1 stamp) | NOT RE-RUNNABLE until #1039; re-stamped as retracted |
| #501 n-body stages (6 stamps) | NOT RE-RUNNABLE until #1039; re-stamped as retracted |

## 6. EGGIE model variants (#1040 work, lead ruling): PRE-REGISTRATION (before any variant run)

Same root search as 3.2 (cycle E>G | G>G | G>I | I>E, all 50 revolution/branch combinations, 12
seed phases over one synodic period, exact roots with residual < 1e-9 km/s and re-fly miss < 1 km,
gate at 25 km and at the project floors). Reported per variant: the root count, the nearest-root
distance to Table 4 (max over the three moons of |V_inf - Table 4|), and the nearest gate-passing
root. One commit per variant, data `data/1023_eggie/pc_roots_<variant>.json`.

- (i) `tsyn705`: T_syn = 7.05 d as printed (p.2). Smas rebuilt consistently for that period: with
  Delta = 5.2 deg, n_I = (8 pi + Delta)/T_syn, n_E = (4 pi + Delta)/T_syn, n_G = (2 pi + Delta)/T_syn,
  and a = (mu / n^2)^(1/3) with the lane mu. Io is therefore NOT the registry Io (a_Io
  417,834 km instead of 421,800; computed in the run). Cycle 4 x 7.05 = 28.20 d (Table 4 total
  28.22 d). Laplace angle 180 deg.
- (ii) `laplace`: the consistent model (3.1), Laplace angle swept 0-345 deg in 15-deg steps (Ganymede's
  initial phase = (angle - lambda_I + 3 lambda_E)/2 with Io = Europa = 0). The nearest root is
  reported over all angles and per angle.
- (iii) `coded`: the OLD coded model, the `ideal_moon_smas` radii with T = 4 x `ideal_t_syn()` =
  28.017 d. Laplace angle 180 deg. The date residual matches V_inf MAGNITUDES at junctions only, so it
  is solvable although the configuration does not repeat rigidly. A root here is quasi-periodic, not
  exactly periodic.

Readings, fixed now:
- A root within 0.5 km/s of Table 4 in (iii) and in no other variant: the paper's own model is the
  Ganymede-period one, and its EGGIE is quasi-periodic in the rigid-rotation sense. This would explain
  the #480 history.
- A root within 0.5 km/s in (i) or (ii): EGGIE is reproducible in that variant, and the continuous
  steps 2-3 of sec. 3 become runnable there (a separate pre-registration).
- No root within 0.5 km/s in any variant: the Table 4 object is unreproduced under every model we can
  build; ESCALATION to the owner.

Expected outcomes, stated now: (i) a near root with probability 0.4; (ii) at some angle 0.35;
(iii) 0.35. Overall I expect at least one variant to give a root within 0.5 km/s, with probability
about 0.6.
