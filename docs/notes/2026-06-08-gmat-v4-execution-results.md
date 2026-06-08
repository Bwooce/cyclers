# GMAT V4 execution results — Aldrin + S1L1, B-plane flyby targeting (#171/#174)

**Type:** execution + results (the V4 lane RUN, not just built). Builds on the
design `docs/notes/2026-06-08-gmat-v4-design.md` and the plan
`docs/superpowers/plans/2026-06-08-gmat-v4-flyby-targeting.md`. **GMAT was actually
run** (the design's "GMAT not installed" blocker is RESOLVED — see
`memory/reference_gmat_install.md`).

**No catalogue / `validate.py` / `spec.md` writeback is performed here** — that is a
separate human-reviewed step (#175 owns those files). This note records the V4
verdicts and the recommended `_LEVEL_EVIDENCE` evidence chain for review.

---

## 0. Environment (the independent external check)

- **Tool:** NASA GMAT R2022a, console build (Build Date Jan 10 2023), at
  `~/GMAT/R2022a/`. Independent codebase + ephemeris (**DE405**, loaded
  automatically) — the spec §14 V4 independence requirement.
- **Headless invocation:** `cd ~/GMAT/R2022a/bin && env -u DISPLAY ./GmatConsole
  --run <script>`. Runs with no X display, writes the `ReportFile` to
  `~/GMAT/R2022a/output/`. The GMAT install, generated `.script` files, and reports
  are NOT committed (they reference the 390 MB install; out of CI by choice, Q7).
- **Lane pieces (committed, unit-tested without GMAT):**
  `src/cyclerfinder/verify/bplane.py` (Jones Eqs.1-5 kernel),
  `scripts/gmat_v4_generate.py` (generator), `scripts/parse_gmat_report.py`
  (parser + two-part V4 predicate), and their tests.

**GMAT ran clean headless** for every script generated. The lane is end-to-end
executable in-environment, not manual-only.

---

## 1. The lane (what runs)

`generate` → `GmatConsole --run` → `parse`:

1. **Generate.** `uv run python scripts/gmat_v4_generate.py aldrin --out <f.script>`
   (single Mars flyby) or `... s1l1 --out <dir>` (7 per-Mars-flyby scripts). Each
   script seeds a **Mars-relative incoming hyperbola** at the SOI from that node's
   sourced `v_inf^-`, then a GMAT `Target`/`Vary`/`Achieve` block (a
   `DifferentialCorrector`) varies a periapsis TCM and **Achieves the outgoing
   `v_inf^+` velocity** at the outbound SOI (Beeson: BC at flyby periapse, Mars
   gravity in the EOM). The Jones B-plane goal (Ŝ,T̂,R̂ / θ_B / B·R,B·T) is recorded
   as the approach geometry in the script comments; the achieved goal is the next
   leg's outgoing asymptote, which is the maintenance the flyby must deliver.
2. **Run GMAT** headless. The DC drives the defects below the Jones continuity bar
   (1e-3 km / 1e-6 km/s).
3. **Parse.** `uv run python scripts/parse_gmat_report.py <combined-log> [--ref-dv]`
   reports convergence + summed maintenance ΔV + the two-part V4 verdict.

---

## 2. Aldrin (#134) — V4 result

- **Script:** single Mars return flyby, seed `v_inf` 8.88 km/s, required turn 93.0°
  (the documented lap-0 powered member,
  `docs/notes/2026-06-07-aldrin-continuation-v3-evidence.md` §2/§4).
- **GMAT run:** **CONVERGED** — "The Targeter converged!" in 3 iterations, defects
  below 1e-6 km/s on the outgoing velocity components.
- **GMAT maintenance ΔV = 0.1753 km/s** (the converged periapsis TCM that delivers
  the full 93° return-flyby turn).
- **Reference (OUR value, under external check):** 2.9138 km/s/synodic → ±5% band
  **2.768–3.060 km/s**.
- **V4 predicate (band):** **FAIL** — 0.175 km/s is ~16× BELOW the band.
- **V4 predicate (convergence-only):** **PASS**.

### Reading (honest, not tuned)

The 2.9138 km/s reference is a **heliocentric asymptote-rotation impulse-schedule**
figure (`maintain.py` "asymptote" model: the bend rotated at the V_inf asymptote,
no Oberth credit, summed across the heliocentric legs). GMAT does the burn **at
periapsis with full Oberth credit**, deep in Mars's well, where a small impulse
produces a large asymptote rotation — so the same 93° turn costs **~16× less**
(0.175 vs 2.9138 km/s). GMAT's independent high-fidelity answer says our 2.9138 km/s
reference is **conservative (an over-estimate)**, not wrong: the cross-check is that
the flyby targeting CLOSES in real dynamics, which it does. (Our own
`core/flyby.dv_powered_flyby_periapsis` patched-conic Oberth estimate, 9.7 km/s at
this geometry, is itself an over-estimate vs GMAT's 0.175 — it computes the cone at
the SAFE periapsis only; GMAT lets the DC find the actual periapsis depth.)

This is exactly the "the reference is OUR value; GMAT is the external check" case the
design called for. The band FAIL is an honest finding, reported per the golden rule
(do not tune to hit the band).

---

## 3. S1L1 (`russell-ch4-4.991gG2`) — V4 result

- **Scripts:** 7 per-Mars-flyby scripts, one per App-C Mars encounter
  (legs 3,6,9,12,15,18,21), each seeded from the consecutive `APPC_LEGS`
  (v_inf^-, v_inf^+) pair (the sourced Russell 2004 Appendix-C block).
- **GMAT run:** **all 7 CONVERGED** (7 ok, 0 failed; ~3 s total headless).
- **Per-flyby converged maintenance ΔV (km/s):**

  | Mars node | \|v_inf in\| | \|v_inf out\| | \|Δ\|v_inf\|\| | GMAT TCM |
  |---|---|---|---|---|
  | Mars3  | 5.818 | 5.248 | 0.570 | 0.334 |
  | Mars6  | 5.764 | 7.693 | 1.930 | 1.273 |
  | Mars9  | 3.752 | 4.657 | 0.905 | 0.858 |
  | Mars12 | 5.532 | 3.198 | 2.334 | 1.671 |
  | Mars15 | 6.946 | 6.263 | 0.683 | 0.401 |
  | Mars18 | 5.121 | 8.046 | 2.925 | 2.157 |
  | Mars21 | 4.646 | 3.219 | 1.427 | 0.598 |
  | **sum** | | | | **7.292 km/s** |

- **Reference (OUR value):** 62 m/s patched-conic horizon-TCM over 7 cycles → ±5%
  band **58.9–65.1 m/s**.
- **V4 predicate (band):** **FAIL** — 7292 m/s is ~118× ABOVE the band.
- **V4 predicate (convergence-only):** **PASS** (the B-plane-targeted run closes
  with bounded, recorded per-cycle TCM = 7.29 km/s — the figure GMAT produces for
  the first time).

### Reading (honest, not tuned)

The per-flyby TCM tracks the **v_inf MAGNITUDE jump between consecutive App-C
nodes** (0.57–2.93 km/s; the App-C v_inf "BREATHES 3.2–8.0 km/s across the seven
cycles" — `s1l1_corrected.py` docstring). A flyby cannot change |v_inf|
ballistically, so each independent B-plane-targeted Mars encounter that must deliver
the **next leg's exact sourced v_inf^+** pays that magnitude change as a real TCM.

**The 62 m/s patched-conic horizon-TCM does NOT survive per-flyby high-fidelity
B-plane targeting.** The 62 m/s (#169) was measured on a CONTINUOUS single-seed
chain where the flyby patch absorbed magnitude differently; treating each App-C node
as a standalone flyby that re-anchors to its own sourced v_inf^+ surfaces the full
node-to-node |Δv_inf| as fuel. This corroborates the standing S1L1 finding that the
App-C block is a **per-leg reproduction recipe (each leg re-anchored to its own
sourced v_inf), not a single continuously-flown low-ΔV cycler** (memory:
`project_s1l1_realeph_closure_blocker` — "S1L1 is actually MULTI-ARC ...
single-ellipse never closed"). The 7.29 km/s is the honest cost of forcing the
App-C node sequence node-by-node in real Mars-flyby dynamics.

---

## 4. V4 verdicts (summary)

| Row | GMAT converged? | GMAT ΔV | Reference (OUR) | Band | V4 (band) | V4 (conv-only) |
|---|---|---|---|---|---|---|
| Aldrin #134 | YES (3 it.) | 0.175 km/s | 2.9138 km/s | 2.768–3.060 | **FAIL** (16× low) | **PASS** |
| S1L1 4.991gG2 | YES (7/7) | 7.292 km/s | 62 m/s | 58.9–65.1 m/s | **FAIL** (118× high) | **PASS** |

Both rows **converge** in GMAT's independent high-fidelity flyby targeting (the
Beeson primary predicate, Jones 1e-3 km / 1e-6 km/s). Neither reproduces OUR
reference ΔV within the self-declared ±5% band — Aldrin far below (Oberth credit at
periapsis), S1L1 far above (App-C node-to-node |v_inf| jumps). **Both band results
are reported as honest external-check findings; nothing was tuned to hit a band.**

---

## 5. Recommended V4 writeback (for #175 / the main session — NOT applied here)

The V4 gate is two-part by design: convergence AND band. Because **both rows FAIL
the band check**, a strict two-part V4 PASS is NOT supported for either row. Two
honest options for the reviewer:

1. **Convergence-only V4 (recommended record).** Both rows CONVERGE in GMAT's
   independent high-fidelity model — the trajectory + flyby targeting close. If the
   project accepts convergence-only V4 for rows whose ΔV reference measures a
   different quantity than GMAT's per-flyby TCM (the documented split for the
   Mars-perturbed arm, design §4/§6), then both qualify for a **V4-convergence**
   stamp recording GMAT's produced figure (Aldrin 0.175 km/s; S1L1 7.29 km/s) — NOT
   a band-confirmed V4.

2. **Hold V4 OPEN, record the external-check finding.** The cleaner, stricter
   reading: the two-part predicate is not satisfied, so canonical V4 (band-confirmed)
   stays OPEN. The GMAT run is recorded as evidence that (a) the flyby geometry is
   physically realizable (converges) and (b) OUR ΔV references measure a different
   quantity than GMAT's per-flyby periapsis TCM — Aldrin's 2.9138 is conservative
   (no Oberth), S1L1's 62 m/s under-counts the App-C node-to-node |Δv_inf|.

**Recommendation:** option 2 for the strict gate, with the GMAT numbers recorded as
the external-check evidence. The proposed `_LEVEL_EVIDENCE` evidence string shape
(if a convergence-only V4 is adopted, option 1) per row:

- **Aldrin** `aldrin-classic-em-k1-outbound`:
  `"V4-converged (GMAT R2022a, DE405): single Mars return-flyby B-plane Target
  converged 3 it. to Jones 1e-6 km/s; periapsis-TCM maintenance 0.175 km/s
  (#171/#174 results note). OUR ref 2.9138 km/s is a heliocentric no-Oberth
  asymptote-rotation figure; GMAT's periapsis-Oberth answer is ~16x lower. Band
  (±5% of 2.9138) NOT met — convergence-only."`
- **S1L1** `russell-ch4-4.991gG2`:
  `"V4-converged (GMAT R2022a, DE405): 7/7 per-Mars-flyby B-plane Targets converged
  to Jones 1e-6 km/s; summed per-flyby maintenance 7.29 km/s (#171/#174). OUR ref
  62 m/s patched-conic horizon-TCM (#169) NOT met (118x) — the App-C node-to-node
  |Δv_inf| (0.57-2.93 km/s) surfaces as fuel when each node is targeted
  independently; corroborates the MULTI-ARC App-C-reproduction reading.
  Convergence-only; replaces the V3 entry's 'canonical V4 (GMAT) remains OPEN' only
  if convergence-only V4 is accepted."`

The over-claim guard (`validate.py:661`) requires a recorded evidence pointer for
any V1+ level; the strings above cite GMAT version + ephemeris + convergence +
parsed ΔV + the band verdict, the shape of the existing V3 S1L1 entry
(`validate.py:615`).

---

## 6. Reproduce

```sh
export PATH="$HOME/.local/bin:$PATH"
# Aldrin
uv run python scripts/gmat_v4_generate.py aldrin --out /tmp/aldrin_v4.script
cd ~/GMAT/R2022a/bin && env -u DISPLAY ./GmatConsole --run /tmp/aldrin_v4.script > /tmp/aldrin.log 2>&1
cat /tmp/aldrin.log ~/GMAT/R2022a/output/gmat_v4_aldrin-classic-em-k1-outbound_MarsReturn.report > /tmp/aldrin_combined.log
cd -; uv run python scripts/parse_gmat_report.py /tmp/aldrin_combined.log --ref-dv 2.9138
# S1L1 (7 scripts)
uv run python scripts/gmat_v4_generate.py s1l1 --out /tmp/s1l1_v4
for s in /tmp/s1l1_v4/*.script; do (cd ~/GMAT/R2022a/bin && env -u DISPLAY ./GmatConsole --run "$s"); done
cat ~/GMAT/R2022a/output/gmat_v4_russell-ch4-4.991gG2_Mars*.report  # sum the dv columns
```

The generator/parser/kernel + their tests are committed; the GMAT install, scripts,
and reports live in `$HOME`/`/tmp` (not committed). **No CI test invokes GMAT (Q7).**

---

## 7. Honesty / golden discipline

- The only SOURCED anchors are each row's published v_inf nodes (Aldrin documented
  Mars geometry; S1L1 `APPC_LEGS`) and the Jones AAS 17-577 §2.5 continuity
  tolerance (1e-3 km / 1e-6 km/s). Aldrin 2.9138 km/s and S1L1 62 m/s are OUR
  computed values — GMAT is their **independent external check**, never an
  EXPECTED-from-source assertion.
- The ±5% band was declared up front (design §4) and is NOT back-fit. Both band
  FAILs are reported as honest results; nothing was tuned to hit a band.
- The B-plane kernel test is a self-consistency round-trip (Jones tabulates no
  worked example), not a golden.
