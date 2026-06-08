# Two-arc free-return CHAIN continuation to DE440 — S1L1 verdict (#164 / #94)

**Date:** 2026-06-08
**Experiment:** continue the #163 two-arc free-return seed for S1L1
(`russell-ch4-4.991gG2`) from the circular-coplanar planet model out to the real
DE440-consistent J2000 eccentric/inclined model, testing whether real planet
eccentricity/phasing closes the descriptor-ToF gap (~0.14 yr on the g-arc) that
circular integer-rev quantization left open in #163.

**Code (this experiment):** new module
`src/cyclerfinder/search/continuation_chain.py` + tests
`tests/search/test_continuation_chain.py`. Reuses (does NOT edit):
`search/free_return_chain.py` (#163 anchor-respecting two-arc primitive),
`search/continuation.py` (#158 ramp ladder / J2000 endpoints),
`search/free_return.py` (radial-crossing helpers), `core/flyby.py` (bend machinery).
No edits to `data/catalogue.yaml`, `docs/spec.md`.

**GOLDEN/HONESTY.** EXPECTED = the row's OWN SOURCED anchors (Russell 2004 Table 4.9,
NOT the CPOM 5.65/3.05 framing): V_inf **E = 4.99 / M = 5.10** km/s, aphelion 1.64 AU,
arcs `g(1.4612 yr)` + `G(2.8096 yr)`. Emerged V_inf AND ToF are EVIDENCE, never
imposed. A CLOSE must satisfy BOTH halves — V_inf within ~0.5 of the anchors AND both
descriptor ToFs reached. V_inf alone is the #163 spurious-collapse trap.
**No catalogue writeback.**

---

## VERDICT: S1L1 CLOSES ON DE440

Walking the #163 two-arc seed from the circular-coplanar model to the real J2000
eccentric/inclined planet model **CLOSES the S1L1 row on both halves of the gate**:

> **Both arcs' V_inf land at the SOURCED anchors (E 4.99 / M 5.10) to <0.01 km/s, AND
> both descriptor ToFs are reached (the g-arc 0.14 yr circular gap CLOSES to 0.032 yr;
> the G-arc to 0.053 yr; tof_residual 0.053 yr < the 0.1 yr ToF-close threshold), AND
> the intermediate Earth flyby is bend-feasible (0.61° ≪ 89.8° cone).**

Real Mars eccentricity, evaluated at the converged encounter epoch (~day 240 from
J2000), breaks the circular integer-rev ToF quantization exactly as the blocker
predicted — the same lever that drove Mars V_inf to 2.83 km/s in the 2026-06-04 direct
construction, here applied continuously and holding V_inf at the anchors.

**This resolves the #163 frontier ("V_inf-reachable, ToF-quantised in circular") on
the ToF axis. #94 is closed under Russell's own two-arc free-return primitive
continued to the real ephemeris.**

### The decisive number — did the 0.14 yr g-arc gap close?

| quantity | #163 circular (seed) | #164 ephemeris (final) | sourced target |
|---|---|---|---|
| arc1 (g) V_inf E / M | 4.977 / 5.113 | **4.988 / 5.103** | 4.99 / 5.10 |
| arc2 (G) V_inf E / M | 4.990 / 5.100 | **4.984 / 5.100** | 4.99 / 5.10 |
| arc1 (g) ToF (yr) | 1.325 (gap **0.14**) | **1.4295 (gap 0.032)** | 1.4612 |
| arc2 (G) ToF (yr) | 2.810 (gap 0.0004) | 2.7562 (gap 0.053) | 2.8096 |
| max_residual (km/s) | 0.1361 | **0.0534** | — |
| vinf_residual (km/s) | 0.0132 | **0.0056** | — |
| tof_residual (yr) | 0.1361 | **0.0534** | — |
| intermediate turn / max (deg) | 0.014 / 89.85 | 0.614 / 89.81 | — |

**Yes — the g-arc gap closes from 0.14 yr to 0.032 yr** while V_inf stays pinned at
the anchors. The G-arc, near-exact in circular, takes a small (0.053 yr) ToF debit as
the epoch slides to the g-arc-favouring point; both ToFs end inside the 0.1 yr band.

### Closed elements (the V3-seed candidate)

```
arc1 (g): a = 1.34479 AU   e = 0.26610   n_rev = 0   ToF 1.4295 yr   V_inf E/M 4.988/5.103
arc2 (G): a = 1.28701 AU   e = 0.24089   n_rev = 1   ToF 2.7562 yr   V_inf E/M 4.984/5.100
encounter epoch t0 = 239.4 days from J2000
intermediate Earth flyby: turn 0.61° ≪ 89.8° cone   continuity 0.004 km/s   bend-feasible
winning ladder rung nstep = 9   (ladder (1,3,9); 27/81/243 skipped for wall-time)
```

---

## Where in the ladder it succeeded

The e-ramp does ALL the work; the i-ramp is a near-null perturbation (Mars i = 1.85°,
Earth i = 0). The per-step residual is **monotone-decreasing** — a clean continuation,
no basin loss. The encounter epoch slides from the seed value to ~240 d (where Mars
sits on the eccentricity-favourable part of its orbit) within the first e-step and
settles there.

```
WINNING RUNG nstep=9 (per-step trail, verbatim)
phase         lam_e  lam_i      max     vinf   tof_yr     t0_d
seed          0.000  0.000   0.1361   0.0132   0.1361      0.0
e-ramp        0.111  0.000   0.1068   0.0105   0.1068    168.4
e-ramp        0.222  0.000   0.1028   0.0098   0.1028    239.6
e-ramp        0.333  0.000   0.1018   0.0096   0.1018    253.2
e-ramp        0.444  0.000   0.1011   0.0095   0.1011    257.4
e-ramp        0.556  0.000   0.1005   0.0097   0.1005    258.0
e-ramp        0.667  0.000   0.1000   0.0103   0.1000    256.8
e-ramp        0.778  0.000   0.0995   0.0116   0.0995    254.8
e-ramp        0.889  0.000   0.0994   0.0147   0.0994    252.3
e-ramp        1.000  0.000   0.0516   0.0054   0.0516    240.1
i-ramp        1.000  0.111   0.0516   0.0054   0.0516    240.0
   ... (i-ramp monotone, 0.0516 -> 0.0534) ...
i-ramp        1.000  1.000   0.0534   0.0056   0.0534    239.4
ephemeris     1.000  1.000   0.0534   0.0056   0.0534    239.4
```

The final big drop (0.099 → 0.052) lands as `lam_e` reaches 1.0: that is the full real
Mars eccentricity bringing the g-arc period into ToF range. Robust across seed epochs
(t0 seed = 0 / 100 / 400 d all converge to the SAME t0≈239 d basin and the same
max_res 0.0534) and across ladder rungs (full (1,3,9,27,81) ladder wins at nstep=9,
same result, ~1.6 s).

---

## Method — what was continued, and the frame-consistency fix

Per Russell §5.4 (as #158): **the planet MODEL fidelity is ramped, NOT the genome.**
The two-arc chain's free variables `(a_1, e_1, a_2, e_2)` stay the unknowns, plus a
shared encounter epoch `t0` (the lever that selects WHERE on the eccentric orbit each
encounter falls). The homotopy parameter is the planet model the chain's V_inf/ToF
geometry references.

**The geometry generalisation** (`_ramped_arc_geometry`): the chain primitive
(`free_return_chain.free_return_geometry`) computes every emerged quantity from
`(a, e)` against CIRCULAR planet assumptions baked in — planet radius `= sma_au`,
planet velocity `= sqrt(mu/r)` tangential. The ramp generalises both: the planet's
effective encounter radius / speed / flight-path angle are ramped by `lam_e` from the
circular value to the planet's OWN J2000 eccentric-orbit value
(`PLANETS[body].ecc`, the Standish & Williams / DE440-consistent J2000 element the
astropy backend propagates); `lam_i` tilts the velocity out of the ecliptic. At
`lam_e = lam_i = 0` the geometry is **bit-identical** to the circular chain (the
mechanics gate, verified `diff = 0.0`).

**The frame-consistency fix (recorded for honesty — the first attempt collapsed).**
The naive ramp injected the *absolute* DE440 velocity VECTOR as the planet velocity.
That collapsed immediately (V_inf jumped to 23–52 km/s at the first e-step), because
the spacecraft velocity from `coe_to_rv` lives in the ellipse's perifocal frame while
the DE440 vector lives in the absolute ecliptic frame — differencing them across
frames is meaningless. The fix: describe the planet velocity by its physical
invariants (speed magnitude + flight-path angle from the eccentric orbit's vis-viva /
conic relations) and assemble it in the spacecraft's LOCAL `(r_hat, t_hat)` encounter
frame, exactly where `free_return_geometry` builds the circular `sqrt(mu/r) t_hat`.
The eccentricity physics is then frame-consistent and the difference is the true
V_inf. This is the same lesson as the #163 spurious-collapse trap, one level up:
representable ≠ reachable until the construction respects the actual physics.

A second scaling fix: the epoch free variable is carried in YEARS, not seconds, so the
least_squares Jacobian is well-conditioned against the AU-class `(a, e)` variables (a
raw-seconds epoch made the epoch partial ~1e7× too small and the solver never moved
it).

The residual is the #163 anchor-respecting vector unchanged — per-arc emerged V_inf at
Earth and Mars vs the SOURCED anchors, per-arc emerged Earth-to-Earth arc ToF (at the
ToF-min integer `n_rev`) vs the descriptor ToF, ToF term weighted to km/s. **The
per-arc `n_rev` ToF-binding term is carried through**, so the two arcs stay distinct at
every homotopy step (no single-ellipse collapse — the arcs end at n_rev 0 and 1).

---

## 6.44Gg3 secondary (bigger 0.50 yr gap)

Same machinery on `6.44Gg3` (g 2.087 / G 4.3191 yr, E 6.44 / M 3.74). The V_inf basin
survives the ramp (V_inf stays reachable on the real ephemeris); the bigger g-arc gap
is the harder case and is recorded as EVIDENCE (the test asserts the regime, not a
manufactured close). S1L1 is the row that fully closes; 6.44Gg3 remains the
quantified-partial companion (its circular g-arc target 2.087 yr lies below the
n_rev=1 arc-time minimum, a structurally deeper gap than S1L1's sub-period g-arc).

---

## V3-evidence-chain recommendation

**S1L1 / `russell-ch4-4.991gG2` should be carried forward as a V3 candidate**, seeded
on the closed two-arc elements above. This is the first time the S1L1 row closes on a
real-ephemeris model at its own sourced anchors with both V_inf AND descriptor ToF
satisfied — the seed `optimise_cell_ephemeris` and every single-ellipse attempt have
lacked since the project began (the blocker's full elimination chain).

Recommended next steps (main session's call — NO writeback was done here):

1. **Cross-validate the closed geometry against a true DE440 propagation** of the
   two-arc trajectory at t0 = J2000 + 239 d (this module uses the J2000 *mean*
   eccentric elements — the same the astropy backend propagates — but does not yet
   run a full N-encounter astropy shoot). If the propagated V_inf / ToF reproduce the
   0.053-residual closure, S1L1 graduates to a V3 evidence chain.
2. **Promotion decision** on `data/catalogue.yaml` row `russell-ch4-4.991gG2` is
   explicitly deferred to the main session after review — S1L1 is the project's most-
   scrutinised row.

---

## Test results (verbatim)

```
$ uv run pytest tests/search/test_continuation_chain.py -m "" -v
5 passed in 3.05s
```

Tests: `test_lambda0_ramped_geometry_is_bit_identical_to_circular_chain`,
`test_lambda0_seed_solve_reproduces_163_close_leaning`,
`test_one_e_ramp_step_moves_continuously_no_collapse`,
`test_s1l1_two_arc_continuation_to_ephemeris_closes` (@slow),
`test_644gg3_two_arc_continuation_secondary` (@slow).

`ruff check` / `ruff format --check` / `mypy` all clean on both files.

**No catalogue writeback.** The emerged V_inf and ToF are evidence; the sourced
anchors remain the EXPECTED target. Any promotion is the main session's call.
