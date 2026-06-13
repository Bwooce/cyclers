# Discovery daemon (#253/#254) — three correctness fixes (#259)

The first repeated-moon hunt (2026-06-14, #253 daemon driving the #254 genome
over Jupiter + Saturn) produced 9 SILVER candidates that were ALL artifacts and
were purged. Three confirmed bugs in
`src/cyclerfinder/search/discovery_campaign.py` (`RepeatedMoonTarget`) had to be
fixed before a re-run is meaningful. This note records the fixes, the tests, and
the bounded smoke. No catalogue writeback; no real `data/` registry writes — the
smoke used temp paths only.

## Fix A — `--primary` moon-set bug

**Bug:** `RepeatedMoonTarget` swapped the primary GM (`PRIMARIES[primary]`) but
kept a hardcoded Galilean moon default `("Io","Europa","Ganymede","Callisto")`.
A Saturn campaign therefore enumerated JOVIAN moons with Saturn's mu — physically
meaningless.

**Fix:** new `_registry_moons_for(primary)` derives the moon set from
`core.satellites.SATELLITES` filtered to the moons that actually orbit `primary`,
returned sorted (deterministic). `moons` now defaults to `()` and is resolved in
`__post_init__` (frozen dataclass -> `object.__setattr__`); an explicit `moons=`
still overrides. An unknown primary with no registered moons raises.

**Test result (`test_fix_a_saturn_enumerates_only_saturn_moons`):**
- Saturn target moons = `(Dione, Enceladus, Mimas, Rhea, Tethys, Titan)` —
  Saturnian only, zero Jovian.
- Jupiter target moons = `(Callisto, Europa, Ganymede, Io)` — Jovian only.

## Fix B — trivial zero-rev degeneracy

**Bug:** all 4 purged Jovian SILVER hits were `n_rev=[0,0]` — the trivial
direct-transfer corner, not the multi-rev cyclers targeted (Liang CGE >=10
cycles).

**Rule chosen (documented in the enumerator):** a candidate must have **at least
one leg with `n_rev >= 1`**. The all-zero-rev tuple is skipped. Reference for
"what counts as a real repeated-moon cycler" = the Liang members themselves,
whose legs are all `n_rev = 1`; so this rule excludes the junk but keeps the
genuine class. (We exclude rather than segregate: the all-zero corner is the
prior single-ellipse zero-rev genome's territory, already swept.)

**Test results:** `test_fix_b_all_zero_rev_excluded` (with a grid that includes
0, no enumerated candidate is all-zero-rev) and
`test_fix_b_liang_style_legs_still_allowed` (all-`n_rev=1` closed cycles still
enumerate).

## Fix C — gate was V_inf-continuity-only (necessary, not sufficient)

**Investigation.** The #254 close() set the canonical residual = worst per-flyby
V_inf-MAGNITUDE continuity defect over the INTERIOR flybys of an **open path**
`[m_0, ..., m_{k-1}]`. Two gaps:

1. **Open-path topology.** The enumeration produced open paths and called them
   cyclers. A repeated-moon cycler is a CLOSED loop that returns to its anchor
   moon so the tour repeats (Liang CGE = Callisto-Ganymede-Callisto-Europa-
   Callisto — starts AND ends at Callisto). An open path has no anchor to
   re-close on.
2. **Wrap-around flyby never checked.** Even given a closed sequence, the old
   residual linked only interior flybys; it never compared the cycle-start
   outbound V_inf (`vinf_out[0]`) against the cycle-end inbound V_inf
   (`vinf_in[-1]`) — the SAME anchor flyby across the cycle boundary. So a
   trajectory whose closing flyby did not match its opening flyby could still
   score a low residual. V_inf continuity at interior flybys is
   necessary-not-sufficient for the loop to close on itself.

**Fix.**
- `_sequences()` now emits only CLOSED cycles: `seq[0] == seq[-1]`, each leg
  changes moons, >=2 distinct moons, `k >= 3`.
- `_close_one_phasing()` now folds the **anchor wrap defect**
  `|vinf_out[0] - vinf_in[-1]|` into the worst-continuity residual, so a sub-gate
  residual means every flyby — INCLUDING the anchor — is V_inf-continuous, i.e.
  the closed loop actually closes on itself in the relative (rotating) frame
  (same anchor moon, matched hyperbolic excess speed = a ballistic re-encounter).

**Why NOT a one-cycle phase-return-to-zero term.** A natural-looking stronger
criterion is "after one cycle ToF every moon returns to its start phase." This is
FALSE for genuine cyclers and would reject the Liang goldens: measured per
member, one cycle (sum of the 4 printed leg ToFs ~100 d) leaves Europa
0.96–1.23 rad (6.4e5–8.3e5 km) from its start phase — the true Liang repeat is a
multi-cycle (synodic, >=10-cycle) recurrence, not a single-cycle phase return.
Imposing single-cycle phase closure would be a WRONG binding constraint (the
orbit-closure-discipline failure mode: every binding constraint in the residual,
but never a constraint the real object does not satisfy). The closed-topology +
anchor-wrap continuity is the correct relative-frame closure for a single
ballistic cycle; the full multi-cycle phase recurrence is the coordinator's
deeper-search concern, not a build-time gate.

**Gate-before vs after.**
- Before: residual = worst V_inf-magnitude continuity over INTERIOR flybys of an
  open path. Anchor wrap dropped; open paths accepted.
- After: candidates are closed cycles (anchor-terminated); residual = worst
  V_inf-magnitude continuity over ALL flybys including the anchor wrap.

**Liang members still pass.** The Liang reproduce-before-search gate
(`moon_cycler_genome.reproduce_before_search_gate`, high-fidelity sourced-ToF
route in `cge_scaffold`) is unchanged and still passes A/B/C:

| member | worst_continuity (km/s) | worst_vinf_res (km/s) | tol (km/s) | passed |
|---|---|---|---|---|
| 111-highperijove (A) | 8.09e-3 | 1.52e-2 | 1.37e-1 | True |
| 110-highperijove (B) | 7.81e-3 | 1.37e-2 | 1.38e-1 | True |
| 111-lowperijove  (C) | 9.03e-3 | 4.82e-2 | 1.38e-1 | True |

(The campaign's coarse generic-ToF close() is NOT a Liang reproducer — its
resonance-ToF model does not match Liang's printed per-leg ToFs — so the Liang
gate correctly stays in the genome module on the sourced route; the campaign fix
is about not accepting open-path / zero-rev / wrap-broken junk.)

## Bounded smoke (temp paths; NOT the full hunt)

`RepeatedMoonTarget(seq_lengths=(3,), n_rev_grid=(0,1), n_phase_samples=6,
tof_resonance_grid=(0.5,1.0,1.5))`, `gate=0.05 km/s`, `max_candidates=12`,
`empty_registry=[]`, temp routing:

- **Jupiter:** moons `(Callisto, Europa, Ganymede, Io)`; bodies in enum = Galilean
  only; `any_all_zero_rev=False`; `all_closed_cycles=True`. enumerated 13,
  evaluated 12, closed 5, failed_close 7, **silver_routed 0**, empty_routed 1.
- **Saturn:** moons `(Dione, Enceladus, Mimas, Rhea, Tethys, Titan)` — Saturnian
  only (Fix A); `any_all_zero_rev=False`; `all_closed_cycles=True`. enumerated 13,
  evaluated 12, closed 12, failed_close 0, **silver_routed 0**, empty_routed 1.

Zero trivial SILVER hits in either smoke (the prior junk is gone); both route a
clean empty-region negative under this tiny bounded box. The real multi-thousand-
candidate hunt is the coordinator's job after this lands.

## Status

ruff check + ruff format clean; mypy clean on the module;
`tests/search/test_discovery_campaign.py` + `tests/search/test_moon_cycler_genome.py`
pass. No catalogue writeback; smoke wrote only temp paths.
