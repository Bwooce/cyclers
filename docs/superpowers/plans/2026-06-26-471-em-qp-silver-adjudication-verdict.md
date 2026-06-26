# #471 — Earth-Moon 3D QP-tori SILVER adjudication verdict

Date: 2026-06-26
Task: #471 — adjudicate the two overlooked, already-V1_qp+V2_qp-validated
Earth-Moon 3D quasi-cycler SILVER candidates from the #320 Vector-A sweep.
Substrate (read-only): #290 (`qp_tori`), #296/#299 (3D family tracer +
bifurcation brackets), #319 (`v1_qp`/`v2_qp`), #320 (the sweep), #328/#434
(`literature_check`), #428 (flyby-altitude references), #466 (energy-walk
confirming Bracket-2 is a real family member).

## Verdict (one line)

**Both candidates are DOUBLE NEGATIVE: (1) a REPRODUCTION of the published
Earth-Moon CR3BP 2D quasi-periodic invariant-tori family (Olikara-Scheeres
2012 / Olikara 2016 / McCarthy-Howell 2019-2021), bifurcating off the
already-rediscovered Antoniadou-Voyatzis 2018 spatial resonant family; AND
(2) NON-transport LIBRATION tori — closest lunar approach 49,576-54,815 km
(28-32 lunar radii, no flyby), rotating-frame encounter speed 0.2-0.4 km/s,
only grazing the outer lunar SOI. NOT transport quasi-cyclers. NO catalogue
row. The C21 pattern (#447) repeats, with a published-family identity added.**

This is the honest Outcome (a)+(b) the brief flagged as most likely.

## The candidates (recovered ICs)

Parent periodic-orbit states (rotating-frame, nondim; the tori are GMOS
displacements of amplitude 5e-4 nd ≈ 192 km around these — far too small to
change any transport verdict):

| | k | step_a | T_a (TU) | T (days) | C_J | freq ratio | parent state (x,0,z,0,vy,0) |
|---|---|---|---|---|---|---|---|
| **Bracket 2** | 4 | +8 | 10.27935 | 44.64 | 3.03196 | -0.218 (irrational) | (-0.80549, 0, -0.23594, 0, -0.13157, 0) |
| **Bracket 10** | 4 | +110 | 9.30832 | 40.42 | 3.12624 | -0.265 (irrational) | (-0.76428, 0, -0.09460, 0, -0.27929, 0) |

Source rows: `data/scan_320_qp_tori_3d_brackets.jsonl` (idx 2 / idx 10) and
`data/family_296_3d_subfamilies_299.jsonl`. Both V1_qp PASS + V2_qp PASS +
irrational frequency ratio = SILVER under the #320 gauntlet. EM mu = 0.0121506,
l = 384,400 km, t = 375,190 s (`cr3bp_system('Earth','Moon')`, gm_de440).

## Gate 1 — LITERATURE NOVELTY → REPRODUCTION (both)

### Parent family is already-adjudicated published
The brackets are Neimark-Sacker bifurcation points ON the #296 (1,1) 3D
Earth-Moon family, which #299 (`data/lit_check_299_3d_family.jsonl`)
adjudicated as a **clean reproduction of Antoniadou-Voyatzis 2018** spatial
resonant periodic orbits in the CR3BP — all 265 members hit that anchor at
confidence 0.85. The Jacobi band (C_J 3.03 and 3.13) sits squarely inside the
published Earth-Moon halo / NRHO / vertical-family band (C ~ 3.0-3.18;
`known_corpus_3d.py`, Howell 1984 / Folta-Bosanac 2015).

### The QP-tori themselves are a published family (live WebSearch, #468 lesson)
`literature_check` has no QP-torus-shaped corpus anchor (documented gap in the
#320 note), so its offline verdict is `inconclusive` — not believable on its
own. Per `feedback_literature_novelty_check_baseline`, I grounded it with a
LIVE WebSearch. Earth-Moon CR3BP 2D quasi-periodic invariant tori — computed
by the GMOS / stroboscopic-map collocation method (= our `correct_qp_torus`
substrate), bifurcating off L1/L2 spatial periodic families (halo / vertical-
Lyapunov / quasi-halo), continued by rotation-number / frequency-ratio — are a
mature, thoroughly published field:

* **Olikara & Scheeres 2012**, "Numerical method for computing quasi-periodic
  orbits and their stability in the restricted three-body problem," Adv.
  Astronaut. Sci. 145:911-930 — the canonical stroboscopic-map 2D-tori method
  our substrate implements.
* **Olikara 2016** PhD thesis (CU Boulder), "Computation of Quasi-periodic
  Tori and Heteroclinic Connections in Astrodynamics Using Collocation."
* **McCarthy & Howell 2019 / 2021** (AAS 19-329, AAS 21-270) — quasi-periodic
  orbits in the Earth-Moon / Sun-Earth-Moon system for trajectory design.
* **"Dynamics in the Vicinity of the Stable Halo Orbits"** (J. Astronaut. Sci.
  2023, s40295-023-00379-7) — explicitly computes 2D quasi-periodic tori
  bifurcating off the Earth-Moon L2 halo with rotation-number continuation,
  vertical families of invariant 2-tori, and the halo family born at the
  stability bifurcation. This is precisely the Neimark-Sacker → 2-torus
  structure our brackets reproduce.

**Verdict: REPRODUCTION.** Matched source family: Olikara-Scheeres 2012 EM
QP-tori (parent family Antoniadou-Voyatzis 2018). Not lit-fresh.

## Gate 2 — TRANSPORT UTILITY → LIBRATION (both), the decisive gate

Propagated each parent over its full period (DOP853 rtol/atol 1e-13, 20k
samples, CA refined on a local 2k-point grid). Moon at x_nd = 0.98785, Earth
at x_nd = -0.01215. Lunar SOI = 66,183 km (0.17217 nd); Moon radius 1737 km;
sourced lunar flyby design-CA floor (Genova-Aldrin 2015 farside perilune)
3,000 km.

| | closest Earth | closest Moon | Moon SOI | v_rot @ Moon-CA | z @ Moon-CA |
|---|---|---|---|---|---|
| **Bracket 2** | 177,601 km | **54,815 km** (0.1426 nd) | 66,183 km | 0.404 km/s | 0.0196 nd |
| **Bracket 10** | 200,285 km | **49,576 km** (0.1290 nd) | 66,183 km | 0.219 km/s | 0.0683 nd |

Reading:

* Both DO dip inside the lunar SOI (0.75-0.83 × SOI) — slightly different from
  C21 (#447), which stayed entirely outside SOI at 122,628 km. But the lunar
  closest approach is **49,576-54,815 km = 28-32 lunar radii**, four orders of
  magnitude above the 3,000 km sourced lunar design-CA floor and well outside
  any flyby-altitude reference. There is no lunar flyby / gravity assist here.
* The rotating-frame encounter speed at lunar CA is **0.2-0.4 km/s** — a
  low-energy libration passage, not a gravity-assist encounter (cf. the #320
  Saturn Titan-Rhea SILVER at ~1.7 km/s, or any operational cycler).
* These are exactly what their published identity says they are: 2D
  quasi-periodic tori librating around an L1/L2 spatial periodic orbit. They
  oscillate in the cislunar region (x_nd ranging -0.80→0.90, z_nd up to 0.24)
  but never encounter either primary for transport. The torus offset (≈192 km)
  is negligible against the 50,000 km CA, so the verdict is geometry-robust.

**Verdict: LIBRATION torus, NO transport utility** — `resonant_po` /
quasi-halo-class structure, NOT a transport `quasi_cycler`. The C21 / #447
outcome, now with a confirmed published-family label.

## Gate 3 — V3 + admission proposal → N/A

V3 (independent-model cross-check) is gated on a candidate being NOVEL **and**
transport-relevant. Both candidates fail BOTH precursor gates (reproduction +
libration). Per task discipline (never self-admit; V3 only on a survivor), V3
is not run and **no catalogue admission is proposed**. Nothing flagged for the
human gauntlet from #471.

## Classification (final, per bracket)

| | Lit | Transport | V3 | Final classification |
|---|---|---|---|---|
| **Bracket 2** | REPRODUCTION (Olikara-Scheeres 2012 EM QP-tori) | LIBRATION (Moon CA 54,815 km, 0.40 km/s) | N/A | Known-reproduction libration torus — NOT a transport quasi_cycler. No row. |
| **Bracket 10** | REPRODUCTION (Olikara-Scheeres 2012 EM QP-tori) | LIBRATION (Moon CA 49,576 km, 0.22 km/s) | N/A | Known-reproduction libration torus — NOT a transport quasi_cycler. No row. |

## Disposition

* **No catalogue writeback** (no `quasi_cycler` row; data/catalogue.yaml
  untouched by #471). Nothing for the human gauntlet.
* The #320 Vector-A SILVER survivors are now fully adjudicated and closed: the
  V1_qp/V2_qp gauntlet was working correctly (it certified genuine 2-tori), but
  certified-torus ≠ transport-cycler. The transport gate and the live lit-check
  are what convert a SILVER torus into a known-class non-cycler.
* Calibration value: the #320 Vector-A QP-torus lane reproduces a published
  cislunar libration-torus family with no operational cycler in it. This is the
  expected result for tori lifted off L1/L2 spatial periodic orbits (they
  librate near the libration point; they do not transit between primaries).
  The QP-torus avenue for NEW transport quasi-cyclers in Earth-Moon is
  characterised negative at these brackets.

## Files

* `data/scan_320_qp_tori_3d_brackets.jsonl` (idx 2, idx 10 — source rows)
* `data/family_296_3d_subfamilies_299.jsonl` (parent family)
* `data/lit_check_299_3d_family.jsonl` (parent-family reproduction record)
* `data/flyby_altitude_references.yaml` (lunar design-CA floor 3,000 km)
* `docs/notes/2026-06-17-320-first-quasi-cycler-sweep.md` (the sweep)
* `docs/notes/2026-06-16-299-3d-lit-check-bifurcation.md` (parent-family lit-check)
* this verdict.

## Discipline reminders honoured

* Live WebSearch grounded the QP-torus novelty question (the offline corpus has
  no torus anchor; `not-found`/`inconclusive` is necessary-not-sufficient) —
  `feedback_literature_novelty_check_baseline` / the #468 lesson.
* Transport measured against sourced per-body floors, not assumed — the C21 /
  #447 discipline.
* NEVER self-admitted a catalogue row; V3 only on a survivor.
* data/catalogue.yaml + census ratchets untouched (sibling #470 owns those).
