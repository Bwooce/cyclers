# #480 — EIGE ballistic-cycler construction (ideal model): VERDICT

**Date:** 2026-06-30. Status: **DONE** (ideal-model side). The EIGE positive-control
construction — the resume-pointer's open work in
`2026-06-30-480-eige-construction-status.md` — is built and golden-gated. This is the
EIGE analog of the EGGIE ballistic construction (`eggie_ballistic.py`).

## What EIGE is (grounded from the PDF, not the digest alone)
Re-read AAS 17-608 pp.10-11 directly: *"Figure 5 shows a **one synodic period, one rev**
cycler with sequence **EIGE** computed in the real ephemeris using Lambert arcs. The
flybys occur at an altitude of **2,817 km at Io, 13,180 km at Ganymede, and 470 km at
Europa**. The equivalent cycler in the ideal model is ballistic … after **10 repeat
cycles, the ∆V increases to almost 30 m/s**."*

- **Rev count CONFIRMED 1:1** (one synodic period / one rev) — directly from the PDF, not
  inferred. So `a_sc = resonant_sma(1,1,T_syn) = a_Ganymede = 1,055,289 km`, `T_sc = 7.004 d`.
- **No EIGE departure epoch and no EIGE V∞ are printed** — only the three altitudes. So
  an epoch-exact member reproduction is impossible from the paper; the positive control
  is necessarily **class-level**, with the altitudes as the geometry cross-check.

## Corrections to the resume-pointer doc (construct/derive, don't pattern-match)
1. **V∞ regime is LOW excess speed (~5-9 km/s), NOT 12-16.** The construction-status doc
   asserted 12-16 km/s for EIGE by analogy to the sibling EGIEIE (Table 3, 1-syn/**2**-rev).
   Computing `conic_vinf` at `a=a_Gan` directly: the feasible ballistic EIGE sits at
   **Europa ~8.70, Io ~5.14, Ganymede ~7.23 km/s** — the navigation-viable band, same as
   EGGIE. (12-16 only appears at high e≈0.85-0.9, which is not where the feasible member is.)
2. **The topology fits one rev as a CYCLIC sequence.** The doc worried "E-before-I-before-G
   does not appear in a single rev's static order." Resolved: the tour is cyclic and starts
   at the Europa **inbound** crossing; over one full revolution from E-in the order
   E→I→G→E does appear. Settled topology: **Europa-in, Io-in, Ganymede-out, Europa-in (wrap)**.
3. **Per-encounter self-consistency is automatic for EIGE.** Io and Ganymede appear once
   (trivially consistent); the repeated Europa (depart + wrap) closes by the resonance
   (`T_sc = T_syn = 2·T_Europa` → Europa returns to its departure phase exactly). No
   repeated-encounter drift bug class here (unlike EGGIE's 2nd Ganymede).

## The construction (`src/cyclerfinder/search/eige_ballistic.py`)
3 Lambert legs between circular-coplanar moons; free = (φ_Io, φ_Gan, 3 ToFs), Europa
phase the gauge. Residuals = equal-in/out |V∞| at Io and Ganymede + Europa periodicity
seam (3 hard). That leaves a 2-DOF family; the 2 spare DOF are pinned by **softly
targeting the two SOURCED Fig-5 interior altitudes** (Io 2,817 / Ganymede 13,180 km) —
methodologically identical to `eggie_ballistic` softly targeting the sourced Table-4 V∞.

**Result (deterministic; a range of resonant-conic seeds e∈[0.61,0.635] all converge to
the same fixed point):**
- Ballistic to ~1e-13 km/s; cycle closed (seam ~1e-13); total flyby ΔV ~1e-10 m/s.
- **All flyby altitudes feasible** (in the 25-70,000 km window).
- Reaches both targeted interior altitudes (Io 2,817, Ganymede 13,180 km) ballistically.
- **Europa altitude is the free PREDICTION: ~1,323 km** (NOT targeted) — same low order as
  Fig-5's printed 470 km. The residual gap is the ideal↔real-ephemeris difference (Fig 5
  is real-ephemeris; the ideal coplanar model lacks the B-plane freedom). Non-circular.
- ToF seed (days): **[0.422, 1.118, 5.474]**, sum 7.013 ≈ 1 T_syn. Leg plan: all 3 legs
  prograde, n_revs=0, single.

Golden: `tests/search/test_eige_ballistic.py` (5 assertions, all in the default suite).

## What's next (the actual maintenance-ΔV positive control — NOT yet done)
The ideal ballistic EIGE is the SEED. The maintenance number lives in the **real
ephemeris**: feed `eige_tof_seed_days()` + `sequence=EIGE_SEQUENCE` + a 3-leg single-rev
branch plan to `nbody.jovian.chain_cycles` at a real DE440 epoch where the geometry
matches, optimise cycle-1 to ballistic, chain 10 cycles, and report the maintenance-ΔV
curve. Positive control PASSES iff cycle-1 is ~ballistic and the 10-cycle growth is the
**~30 m/s order** the paper prints. (`chain_cycles` is already generalised to arbitrary
sequence + branch plan, `cf1f72a`, and validated on the Liang CGCEC lane, #223.)

## Scope / honesty
No catalogue impact (a reproduced published tour is human-admitted V4-ceiling at most).
The EIGE V∞ values are REPORTED construction outputs (the paper prints none) — not
golden-asserted against a source. Core #480 verdict unchanged.
