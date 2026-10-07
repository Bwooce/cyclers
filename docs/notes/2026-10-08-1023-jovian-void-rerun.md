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

## 3. EGGIE in the consistent ideal model: PRE-REGISTRATION

(pending)
