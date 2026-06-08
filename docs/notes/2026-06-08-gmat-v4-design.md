# GMAT V4 with B-plane flyby targeting — DESIGN NOTE (#171)

**Date:** 2026-06-08
**Type:** design pass (no production code). Deliverable = this note + the executable
plan `docs/superpowers/plans/2026-06-08-gmat-v4-flyby-targeting.md`.
**Extends, does NOT duplicate:** `docs/superpowers/plans/2026-06-06-nbody-harness.md`
Phase D (the GMAT V4 lane: `scripts/gmat_v4_aldrin.py` generator +
`scripts/parse_gmat_report.py` parser, the 5% Aldrin tolerance, Q7 "GMAT is a
manual out-of-CI lane"). This note adds the **B-plane flyby-targeting** layer that
#169 proved is the missing piece, and **generalises the single-row Aldrin lane to
any confirmed cycler** (Aldrin #134 + S1L1 `russell-ch4-4.991gG2`).

---

## 0. The central finding — GMAT is NOT installed in this environment

```
$ which gmat GMAT GMAT_console    → exit 1 (nothing found)
$ ls /opt/GMAT /opt/gmat /usr/local/GMAT ~/GMAT ~/gmat   → none
$ find / -iname "GMAT*" -type f | grep -iE "bin|console|RunGmat"   → empty
```

**GMAT is not present.** This is not a blocker for the *buildable* deliverable — it
is the design's central reality and it confirms Q7 of the nbody-harness plan. The
GMAT execution step is a **manual human step**. Everything we can build and
unit-test without GMAT is:

1. the **script generator** (`scripts/gmat_v4_*.py`) — pure string templating,
   testable with `assert "..." in text`;
2. the **report parser + V4 predicate** (`scripts/parse_gmat_report.py`) — parses a
   fixture report string, applies the declared tolerance;
3. the **B-plane targeting goal computation** (the `Ŝ,T̂,R̂`/`θ_B`/`B·R`/`B·T`
   numbers the generator writes into the GMAT `Target`/`Achieve` block) — pure
   geometry from `(v∞⁻, v∞⁺)`, testable as a self-consistency check;
4. the **manual run protocol** (run-book) — a documented human procedure with the
   GMAT-install prerequisite flagged at the top.

**We do NOT design a CI-automated GMAT lane.** No CI test may invoke GMAT (Q7). The
install + the run are the gating dependency, called out explicitly in §5 and in the
plan's go/no-go gates.

---

## 1. Why B-plane flyby targeting (the #169 finding)

#169 (`docs/notes/2026-06-08-s1l1-continuous-v4-results.md`) ran a continuous
single-seed REBOUND/IAS15 propagation of S1L1 through all 20 App-C nodes:

- **Sun-only / patched-conic flyby patch:** holds cleanly — total maintenance ΔV
  **62 m/s over 7 cycles** (8.8 m/s/cycle, 52% of the 120 m/s V3 budget). This
  STRENGTHENS the V3.
- **Mars-perturbed continuous (real flyby gravity):** every leg *leaving* a Mars
  node **DIVERGES** (>1e7 km, several non-converged); total continuous ΔV ~122 km/s.

**The ~122 km/s is a handoff ARTIFACT, not a fuel cost.** The App-C nodes are
patched-conic states that do not account for the continuous deflection through the
finite Mars encounter. A real continuous flight nulls that deflection *inside the
SOI* — with **B-plane-targeted flyby maintenance** (a small TCM that aims the
post-flyby trajectory at the next encounter's B-plane target), not by burning
122 km/s. The true high-fidelity maintenance ΔV is what a B-plane-targeted GMAT run
measures. This is the number canonical V4 needs, and it is the gap the existing
Phase D (which reproduces a *powered-periodic* Aldrin solution, no flyby
station-keeping) does not yet close.

The Beeson recipe (`docs/notes/2026-06-07-beeson-2015-emtg-gmat-toolchain-mining.md`)
names the operation precisely: the medium→high-fidelity jump is **moving the
boundary conditions from body-center to the periapse of the hyperbolic flyby** so
the flyby body's gravity enters the EOM, then **re-converging the match-point
defects with the NLP** (feasibility-first, homotopy if hard). "Validated" = NLP
convergence; the artifact signal is correction **COST**, not raw closure error.

---

## 2. The generator — what a confirmed cycler hands GMAT

`generate_gmat_v4_script(row, out_path, *, epoch_iso, mode)` (generalising
`generate_aldrin_script` from Phase D Task D.0). Inputs come straight off a
**confirmed catalogue row** (Aldrin or S1L1); outputs a GMAT `.script` following
Beeson's Algorithm-1 two-section skeleton (init + mission sequence with an embedded
`Target`/`Optimize` block). What it carries:

| Item | Source (Aldrin) | Source (S1L1) |
|---|---|---|
| **Seed epoch + departure state** | Aldrin outbound real-DE440 cycler (`maintain.py`/`bvp.py` powered-periodic lock) | App-C leg-2 Earth departure, 2026-12-15, v∞=(−2.278, 5.322, 0.574) km/s (`s1l1_corrected.continuous_chain`, `start_leg_no=2`) |
| **Flyby sequence** | E-M-E (single powered Mars flyby) | E-g(E-E,sub-Mars)-E-G(E-M-E)-E topology, the 20 App-C nodes, 7 Mars encounters |
| **Per-flyby B-plane targeting goals** | one Mars flyby B-plane target from (v∞⁻, v∞⁺) | per-Mars-flyby B-plane target from each App-C node's (v∞⁻, v∞⁺) |
| **Boundary conditions** | at Mars flyby **periapse** (not body-center) | at each Mars flyby **periapse** |
| **Force model** | GMAT's own DE + Sun..Jupiter point masses (V4: external tool's standard set) | same |
| **Reference ΔV (under external check)** | our computed **2.9138 km/s/synodic** | our computed **62 m/s / 7 cycles** patched-conic horizon-TCM (the number GMAT confirms + extends to the Mars-perturbed continuous case) |

**Generalisation contract:** the generator takes a `row` descriptor (id, sequence,
node states/epochs, per-encounter v∞ vectors, reference ΔV) and a `mode`
(`"powered-periodic"` for Aldrin's single-flyby reproduction; `"flyby-station-keep"`
for S1L1's per-flyby B-plane TCM chain). Aldrin is the simplest instance (one
flyby, one mode); S1L1 reuses the same templating with the chain unrolled (Beeson:
loops "completely expanded in script form"). Node states are read the same way the
existing harness does — from the row's sourced v∞ nodes (`correct._vinf_nodes`
b{i}_in/_out vectors for Aldrin; `APPC_LEGS` for S1L1), NOT from any self-computed
attribute.

---

## 3. The B-plane flyby-targeting block (the new kernel)

This is the piece Phase D lacks and #169 demands. The math is the Jones AAS 17-577
kernel (`docs/notes/2026-06-07-jones-aas17-577-method-deepdive.md` §2.2, Eqs.1-5),
body-centered equatorial, pole `k̂=(0,0,1)`:

- `Ŝ = v∞⁻/|v∞⁻|` (incoming asymptote direction)
- `T̂ = (Ŝ × k̂)/‖Ŝ × k̂‖`
- `R̂ = Ŝ × T̂`
- turn `δ = ∠(v∞⁻, v∞⁺)`; periapsis `r_p` solved from Eq.2; `θ_B` from Eq.5:
  `θ_B = atan2(v̂∞⁺·R̂, v̂∞⁺·T̂) − π`
- the **B-vector target** components `B·R = |B| sin θ_B`, `B·T = |B| cos θ_B`, with
  `|B| = r_p √(1 + 2μ/(r_p v∞⁻²))` (impact parameter from periapsis radius).

The generator turns these into the GMAT **Target/Vary/Achieve** structure (mirroring
Beeson's `Create Variable` / `Vary` / `Calculate` / `NonlinearConstraint`
primitives). Per Mars flyby:

```
Target 'FlybyTCM_<i>' DC                         % differential corrector
   Vary    TCM_<i>.V                             % the maintenance impulse components
   Vary    TCM_<i>.N
   Vary    TCM_<i>.B
   Maneuver TCM_<i>(Sat)
   Propagate NearMars(Sat) {Sat.Mars.Periapsis}  % BC at flyby periapse, gravity in EOM
   Achieve  Sat.MarsBPlane.BdotR = <B·R target>  % the Jones B-plane goal
   Achieve  Sat.MarsBPlane.BdotT = <B·T target>
EndTarget
```

The achieved `B·R`/`B·T` are the **goal**; the varied TCM is the **maintenance ΔV**.
Summing the converged `|TCM_i|` across the cycle is the **true high-fidelity
maintenance ΔV** — the number #169's Mars-perturbed DIVERGE said was non-trivial and
that the patched-conic 62 m/s only bounds in Russell's cruise model. For Aldrin
(single flyby, powered-periodic) the same block degenerates to one Target nulling
the Mars-encounter defect, recovering the 2.9138 km/s per-synodic value.

**Run order (Beeson §3.4):** feasibility-first (drive the B-plane/continuity defects
to ~0), then an optional optimality pass; homotopy (Sun+planets → +moon, Jones §2.5)
if the dense corrector stalls.

**Honesty:** Jones tabulates no worked (v∞⁻, v∞⁺)→(θ_B, r_p, B·R, B·T) example, so a
unit test of the B-plane geometry is a **self-consistency check** (round-trip:
target B-vector → propagate → recover the intended v∞⁺ turn), NOT a golden. The
sourced side stays the row's published v∞ nodes.

---

## 4. The parser + the V4 gate

`parse_gmat_report.py` (extending Phase D Task D.1):

- `parse_maintenance_dv(report_text) -> float` — sum of the converged per-flyby TCM
  magnitudes (`Maneuver.TotalDV` lines) over one cycle/synodic period.
- `parse_convergence(report_text) -> bool` — did **every** `Target`/`Optimize` block
  report NLP convergence (defects below GMAT's feasibility tolerance)? Beeson's
  acceptance is **convergence**, not a numeric gap — this is the primary predicate.
- `v4_pass(gmat_dv, ref_dv, *, converged, tol_frac=0.05) -> bool` — **two-part**:
  (a) `converged is True` (the trajectory + flyby targeting closed in GMAT's
  high-fidelity model — Beeson's "validated"); AND (b) `|gmat_dv − ref_dv| ≤
  tol_frac · ref_dv` (the maintenance ΔV reproduces ours within the declared band).

### The sourced acceptance criterion (declared up front, never back-fit)

- **Convergence gate (Beeson, primary):** every flyby B-plane Target / the mission
  Optimize converges — defects driven below GMAT's feasibility tolerance, Jones'
  continuity bar **1.0e-3 km position / 1.0e-6 km/s velocity** (AAS 17-577 §2.5, the
  published ephemeris-continuity tolerance). This is a SOURCED tolerance (Jones), not
  ours.
- **Maintenance-ΔV band (self-declared, per Q5 / Beeson §3.3):** the precedent gives
  no external % tolerance, so we declare **±5%** ourselves, up front:
  - **Aldrin:** ref = 2.9138 km/s → pass band **2.768–3.060 km/s**.
  - **S1L1:** ref = 62 m/s over 7 cycles (the patched-conic horizon-TCM #169
    measured) → pass band **58.9–65.1 m/s** *for the patched-conic confirmation
    arm*. The Mars-perturbed continuous arm has **no prior reference number** (it is
    the number GMAT is being asked to *produce* for the first time) — so its V4 gate
    is **convergence only** (the B-plane-targeted run closes with bounded, recorded
    per-cycle TCM), and the produced TCM becomes the recorded V4 figure, not a band
    check. This is the honest split: you can only band-check a number you already
    have.

`v4_pass` is **V4-PASS iff converged AND (band-check where a reference exists)**.
2.9138 km/s and 62 m/s are OUR computed values — GMAT is the **independent external
check** (independent codebase + ephemeris, spec §14 V4 / spec.md:399). They are the
*reference under check*, never an EXPECTED-from-source assertion.

### How a V4 writeback would be gated

A V4 pass adds an `(id, "V4")` entry to `_LEVEL_EVIDENCE` in
`src/cyclerfinder/data/validate.py` (the over-claim guard at
`validate.py:661` refuses any V1+ level not backed by a recorded evidence pointer).
The evidence string would cite: the GMAT version, the generated script, the parsed
report, the converged per-flyby TCM sum, and the pass/band result — exactly the
shape of the existing V3 S1L1 entry (`validate.py:615`). **The writeback is a
separate human-reviewed step after a real GMAT run; this design does not write it.**
The S1L1 V3 entry already states "canonical V4 (GMAT) remains OPEN" — that line is
what a V4 pass would replace.

---

## 5. Targets, sequence, scope estimate, gating dependency

**Order: Aldrin first** (#134 — simplest: single powered Mars flyby, documented
2.9138 km/s reference, the existing Phase D already targets it), **then S1L1**
(`russell-ch4-4.991gG2` — settles its open V4; the multi-flyby B-plane chain).

| Piece | Scope | Notes |
|---|---|---|
| B-plane goal computation (Ŝ,T̂,R̂, θ_B, r_p, B·R, B·T) from (v∞⁻,v∞⁺) + self-consistency test | **LOW** | pure geometry, Jones Eqs.1-5; reuses `core/flyby.py` turn/r_p; no GMAT |
| Generator: Aldrin powered-periodic script (one flyby Target) | **LOW–MED** | extends Phase D Task D.0; string templating; Beeson Algorithm-1 skeleton |
| Generator: S1L1 flyby-station-keep chain (per-Mars-flyby Target, unrolled) | **MED** | same templating, loop unrolled over 7 Mars nodes; node states from `APPC_LEGS` |
| Parser + two-part V4 predicate (convergence + band) | **LOW** | fixture-string tested; no GMAT |
| Manual run-book (install, generate, run, parse, record) | **LOW** (writing) / **HIGH** (executing) | the GMAT install + run is the human bottleneck |
| **Actually running GMAT + the V4 writeback** | **HIGH — BLOCKED** | GMAT not installed; manual; out of CI |

**Gating dependency (binding):** GMAT is **not installed** and the run is a **manual
human step**. The script generator, B-plane kernel, parser, and run-book are fully
buildable and unit-testable **without GMAT** and are the shippable deliverable of any
execution of the plan. The V4 *verdict* (and any catalogue writeback) waits on a
human installing GMAT (version pinned in the run-book) and running the generated
script. No CI test gates on GMAT (Q7).

---

## 6. Honest blockers / risks

1. **GMAT not installed (the headline).** Buildable: generator + parser + B-plane
   kernel + run-book. Not buildable here: the GMAT verdict. Mitigation: everything is
   unit-tested without GMAT; the run-book pins the install + manual procedure; the
   writeback is a downstream human step.
2. **S1L1 Mars-perturbed continuous has no prior reference ΔV.** Its V4 gate is
   convergence-only (the B-plane-targeted run closes with bounded recorded TCM); the
   produced number becomes the recorded figure, not a band check. You cannot band a
   number you do not yet have. (The patched-conic 62 m/s arm CAN be band-checked.)
3. **GMAT operational brittleness (Q7 Risk 5).** GUI-era app; scripted batch is
   brittle/slow/version-sensitive. Mitigation: manual out-of-CI lane; generator/parser
   unit-tested without GMAT so the logic is covered when the app is absent; GMAT
   version pinned in the run-book.
4. **B-plane geometry has no Jones worked example.** The B-plane unit test is a
   self-consistency round-trip, not a golden. The sourced side stays the row's
   published v∞ nodes; the GMAT run (independent codebase + ephemeris) is the genuine
   external check.
5. **Aldrin 2.9138 / S1L1 62 m/s are OUR values.** They are the reference under
   external check, never EXPECTED-from-source. Documented in every parser docstring
   and the writeback evidence string. This preserves the golden discipline
   (`feedback_golden_tests_sourced_only`): the only sourced anchors are the published
   v∞ nodes and the Jones continuity tolerance (1e-3 km / 1e-6 km/s).

---

## 7. One-line summary

GMAT is **not installed** → the V4 lane is a **manual, out-of-CI** human step; the
buildable deliverable is a **confirmed-row→GMAT generator** carrying per-flyby
**Jones B-plane targets** (Ŝ,T̂,R̂ / θ_B / B·R,B·T) at flyby periapse + a **parser**
whose V4 predicate is **NLP convergence (Jones 1e-3 km / 1e-6 km/s) AND the
maintenance ΔV within ±5% of our reference** (Aldrin 2.9138 km/s; S1L1 62 m/s
patched-conic, convergence-only for the Mars-perturbed continuous arm GMAT produces
for the first time), gated into `_LEVEL_EVIDENCE` by a separate human writeback after
a real run.
