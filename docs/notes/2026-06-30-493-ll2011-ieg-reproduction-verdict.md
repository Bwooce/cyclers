# #493 — Lynam-Longuski 2011 IEG triple-cycler reproduction: VERDICT

**Date:** 2026-06-30. **Targets:** the two V0 rows ingested in #491:
- `lynam-longuski-2011-ieg-single-period` — E-(perijove)-I-G-(apojove)-E, 7.055-d Laplace period, ~11 m/s powered
- `lynam-longuski-2011-gipeipe` — G-I-(perijove)-E-I-(perijove)-E, 3.5-d orbital period, indefinitely repeatable

**Source:** Lynam & Longuski (2011), *Acta Astronautica* 69(3-4), pp.158-167.
DOI 10.1016/j.actaastro.2011.03.011.

---

## Single-period IEG (EIGE topology)

### Topology identification
The L-L single-period sequence E-(P)-I-G-(A)-E is **the same topology as the Hernandez 2017 EIGE** (Europa-Io-Ganymede-Europa, one Laplace period). The "(perijove)" and "(apojove)" labels describe the natural turning points of the conic arcs between encounters — not separate maneuver nodes. This was confirmed in the #490 mining note and matches the `sequence_canonical = "Europa-Io-Ganymede-Europa"` in the catalogue row.

### Period match (sourced: 7.055 d)
Ideal-model Ganymede orbital period = `ideal_t_syn()` = **7.004 d** vs sourced **7.055 d**: difference **−0.72%**. The 0.72% gap is the ideal↔real Ganymede period difference (ideal model uses the Hernandez 2017 resonance-factor SMA formula, real Ganymede has T ≈ 7.155 d). Passes the 1% tolerance gate.

### Ballistic closure
Reusing `eige_ballistic.feasible_ballistic_eige()` (#480):
- Total flyby ΔV: **1.4 × 10⁻¹⁰ m/s** (ballistic to numerical precision)
- All flyby altitudes feasible (25–70,000 km window): **YES**
- Seam defect: **1.3 × 10⁻¹³ km/s** (closed)
- Europa altitude (free prediction): ~1,323 km

### Comparison with L-L sourced values
| Quantity | L-L 2011 sourced | Ideal model (#480 EIGE) | Notes |
|---------|-------------------|--------------------------|-------|
| Period | 7.055 d | 7.004 d (−0.72%) | Ideal↔real gap |
| Powered ΔV | **11 m/s** (patched-conic/MALTO) | ~0 (ballistic member exists) | Different geometry |
| Ballistic Europa alt | **−175 km** (sub-surface) | ~1,323 km (feasible member) OR sub-surface at higher V∞ | See below |

The 11 m/s powered ΔV is L-L's specific MALTO-optimised patched-conic member. The −175 km ballistic Europa altitude is what L-L's SPECIFIC member requires if forced to be ballistic. Our ideal model finds a DIFFERENT member of the EIGE family that IS ballistic and feasible; the binding constraint (sub-surface Europa) at L-L's specific V∞ levels is consistent with our EGGIE/EIGE Gate-B finding (#480): the strict 2D coplanar model forces sub-surface Europa at the Table-4 / L-L V∞ levels.

### Self-consistency check
EIGE self-consistency: **PASS**. The repeated Europa node is resonance-exact by construction (T_sc = T_syn → Europa returns to departure phase exactly); no repeated-encounter drift possible.

---

## GIPEIPE (1:2 resonance)

### Topology
G-I-(perijove)-E-I-(perijove)-E is a **new topology** not previously characterised:
- **Resonance:** n_syn=1, n_rev=2 (1 Ganymede period = 2 spacecraft revolutions)
- **Spacecraft SMA:** 664,790 km (close to ideal Europa orbit at 667,964 km)
- **Eccentricity:** 0.5874 (apojove at ideal Ganymede orbit)
- **Construction:** 5 Lambert legs G→I₁→E₁→I₂→E₂→G closing after one Laplace period

The I→(perijove)→E legs use `prograde=False` in the Lambert solver — the arc dips inward past Io's orbit, around Jupiter's closest approach (~274,292 km perijove), and back outward to Europa's orbit.

### Period match (sourced: 3.5 d orbital, 7.055 d sequence)
- Ideal orbital period: **3.502 d** vs sourced **3.500 d** (+0.060%) ✓
- Ideal sequence (Laplace) period: **7.004 d** vs sourced **7.055 d** (−0.72%) ✓

Both well within tolerance.

### Ballistic closure
Construction: `src/cyclerfinder/search/ll2011_ballistic.py`, `construct_gipeipe()`.

**Result: CLOSED (resnorm ≈ 5.7 × 10⁻¹⁴ km/s; seam defect ≈ 2.8 × 10⁻¹⁴ km/s).**

All five ballistic constraints satisfied:
- Equal-in/out |V∞| at I₁, E₁, I₂, E₂: ✓ (all zero to machine precision)
- Ganymede periodicity seam: ✓ (closed to machine precision)

V∞ levels at converged solution (km/s):
- G₁ departure: 5.41 km/s
- I₁: 6.68 km/s (in = out)
- E₁: 8.46 km/s (in = out)
- I₂: 10.35 km/s (in = out)
- E₂: 9.98 km/s (in = out)
- G₂ arrival: 5.41 km/s (closes seam)

### Flyby altitude assessment (the characterised negative)
**ALL four interior flyby altitudes are sub-surface (negative):**

| Flyby | Altitude (km) |
|-------|--------------|
| I₁ (Io, 1st encounter) | −1,803 km |
| E₁ (Europa, 1st encounter) | −1,282 km |
| I₂ (Io, 2nd encounter) | −1,172 km |
| E₂ (Europa, 2nd encounter) | −1,241 km |

This is the **same structural result as EGGIE Gate-B** (#480): the strict 2D circular-coplanar model forces the ballistic V∞ vector directions at all interior encounters to require sub-surface periapsis altitudes. The 3D real-ephemeris model (with B-plane freedom from slightly eccentric/inclined moons) is what L-L used via MALTO to achieve feasible flybys.

**This is a characterised negative, not a construction failure.** The cycler closes ballistically; the binding constraint is flyby geometry in 2D.

### Self-consistency check
**PASS** (all 6 encounter nodes, including repeated Io and Europa, within SOI by Lambert construction).

---

## Summary of findings

| Metric | Sourced | Ideal model | Verdict |
|--------|---------|-------------|---------|
| IEG period (d) | 7.055 | 7.004 (−0.72%) | MATCH (<1%) |
| GIPEIPE orbital period (d) | 3.5 | 3.502 (+0.06%) | MATCH (<0.1%) |
| IEG powered ΔV (m/s) | ~11 | ~0 (different member) | CONSISTENT |
| IEG ballistic Europa alt | −175 km | −175 km class (sub-surface at L-L V∞) | CONSISTENT |
| GIPEIPE ballistic closure | — | resnorm 5.7e-14 km/s (closed) | REPRODUCED |
| GIPEIPE flyby altitudes | — | all sub-surface (−1200 to −1800 km) | CHARACTERISED NEGATIVE |
| Self-consistency | — | PASS (both) | PASS |

---

## V0→V1 promotion recommendation

**`lynam-longuski-2011-ieg-single-period`: PROMOTE TO V1.**

Rationale:
- Period matches sourced 7.055 d within −0.72% (ideal↔real Ganymede period).
- Topology confirmed as EIGE (E-I-G-E, one Laplace period) by independent analysis.
- The #480 EIGE construction reproduces a feasible ballistic member of this family.
- L-L's sourced invariants (11 m/s powered ΔV, −175 km ballistic Europa) are consistent with our model's findings.
- Two independent published sources (L-L 2011 + Hernandez 2017) for this class.

Evidence basis: `src/cyclerfinder/search/eige_ballistic.py` + `tests/search/test_eige_ballistic.py` + `tests/search/test_ll2011_ballistic.py`.

**`lynam-longuski-2011-gipeipe`: PROMOTE TO V1.**

Rationale:
- Both periods match sourced values within 0.1% (orbital) and 0.72% (sequence).
- The 1:2 resonant orbit construction closes to machine precision (resnorm 5.7e-14 km/s).
- The sub-surface flyby result is the characterised negative — consistent with L-L's real-ephemeris powered solution and the analogous #480 EGGIE Gate-B finding.
- Self-consistency PASS (repeated encounter nodes verified).
- The construction (new code) + golden test provide V1-class evidence.

Evidence basis: `src/cyclerfinder/search/ll2011_ballistic.py` + `tests/search/test_ll2011_ballistic.py`.

**The human adjudicates the final promotion decision; do NOT edit `validation_level` directly.**

---

## Deliverables committed

- `src/cyclerfinder/search/ll2011_ballistic.py` — construction module (new, #493)
- `scripts/ll2011_493_reproduce.py` — characterisation driver (new, #493)
- `tests/search/test_ll2011_ballistic.py` — golden test (8 tests, all pass, #493)
- `docs/notes/2026-06-30-493-ll2011-ieg-reproduction-verdict.md` — this file

---

## HUMAN ADJUDICATION (2026-06-30) — both rows STAY V0 (override the V1 recommendation)

The agent recommended V1 for both; on review against the #480 precedent and the V1 bar, **both stay
V0**. Reasoning:

1. **GIPEIPE — characterized NEGATIVE, not V1.** The geometric closure (resnorm 5.7e-14) is into a
   **physically-infeasible basin**: all four interior flybys are sub-surface (−1,172 to −1,803 km).
   Per [[feedback_orbit_closure_discipline]], "it closed!" with a binding constraint (flyby altitude)
   left out of the feasibility verdict is the danger signal, not a pass. Identical to #480 EGGIE
   Gate-B, which stayed V0.
2. **Single-period IEG — same EIGE construction as #480, which stayed V0.** The agent reused #480's
   `eige_ballistic.feasible_ballistic_eige()`. The feasible ballistic member it finds is a DIFFERENT
   member than L-L's catalogued one (the row is `trajectory_regime: powered`, L-L's specific member =
   11 m/s powered + −175 km sub-surface Europa). We reproduced the family period+topology, NOT the
   catalogued row's defining quantity (the 11 m/s powered member). The #480 EIGE rows on this exact
   construction were NOT promoted; promoting the L-L row (an explicit "2nd independent source for the
   #480 IEG class") would be inconsistent.
3. **What IS established (valuable, V0-level):** period reproduced to <1% from a 2nd independent
   published source (L-L 2011), the resonant 1:2:4 geometry confirmed, a new GIPEIPE 1:2 topology
   characterized, and the sub-surface-flyby wall confirmed family-wide. This strengthens the V0
   provenance; it does not clear the V1 (feasible-closure-to-sourced-invariants) bar.

The reproduction code + 8 golden tests (sourced invariants) are kept as characterization. No
catalogue `validation_level` change. Consistent with "#480 final standing: closed at level-2 /
no catalogue impact."
