# Task #506: Pluto-Charon (3,2) V3-scope assessment

**Date:** 2026-07-01
**Row:** `ross-rt-pc-cycler-32-2026` (V2-ballistic since #505, held for adjudication)
**Task:** scope feasibility of V2→V3 promotion via real-ephemeris independent validation
**Status:** NO V3 RUN. Characterization only — catalogue untouched per task scope.

---

## Summary verdict

**V3 is infra-gated AND the standard V3 methodology is inapplicable as a
stability gate for this CR3BP-periodic-orbit row. Row stays V2-ballistic.**

Neither gate passes:
- Gate (a): No Pluto-Charon SPK kernel in-repo or in the GMAT install.
- Gate (b): The real-ephemeris model IS significantly more perturbative than
  the V2 CR3BP model (by ~5,800×), but this is a model-mismatch signal, not
  an instability signal — V3 as a naive propagation test would not add useful
  stability discrimination.

This is a scientifically correct, publishable characterization — not a
failure. "V3 not applicable" is the honest ceiling for this row.

---

## Gate (a) — Pluto-Charon SPICE kernel availability

**Status: INFRA-GATED. No PLU SPK in-repo or GMAT install.**

Inventory:
- `~/.astropy/cache/cyclerfinder_spice/`: 18 BSP files, zero Pluto-related.
- `~/GMAT/R2022a/data/planetary_ephem/spk/`: DE405/DE421/DE424 planetary
  ephemerides + the Uranian satellite kernel `ura111.bsp` (installed by #335).
  No Pluto satellite kernel.
- `src/cyclerfinder/core/satellites.py` references `PLU060` as the JPL SSD
  *data source* for Charon's physical parameters (GM = 106.1, R = 606.0 km,
  a = 19 600 km). This is a citation, not a SPICE kernel file.

What would be needed for real-eph: the JPL NAIF Pluto satellite SPK, e.g.
`plu058.bsp` (~3.8 MB, covers Pluto/Charon/Nix/Hydra, 1900–2100) — available
at `https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/`. It
is not in-repo and would need a dedicated install script analogous to
`scripts/install_uranian_spice.sh` (#335).

Conclusion: **infra-gated**. The kernel is trivially downloadable but not
installed; without it no real-eph Pluto-Charon integrator can run.

---

## Gate (b) — Meaningfulness: real-eph vs V2 CR3BP-IAS15 discrimination

**Status: Real-eph is measurably MORE perturbative than V2 (by ~5,800×),
but this is model mismatch, not instability. V3 as currently defined is NOT
applicable as a stability gate for this CR3BP-periodic-orbit row.**

### System parameters

| Quantity | Value | Source |
|---|---|---|
| a_PC (Charon SMA) | 19 600 km | JPL SSD PLU060 |
| GM_system (Pluto+Charon) | 975.5 km³/s² | JPL DE440 |
| GM_Charon | 106.1 km³/s² | JPL SSD PLU060 |
| μ = GM_Charon/GM_system | 0.108 765 | derived |
| T_Charon | 6.389 d | derived |
| CR3BP time unit t_s | 1.017 d | T_Charon / (2π) |
| (3,2) orbit T (normalized) | 11.833 | catalogue, #494 |
| (3,2) orbit T (physical) | 12.033 d | T_nd × t_s |
| e_PC (Charon eccentricity) | < 5×10⁻⁵ | Brozović et al. 2015, Icarus 246:317 |

### V2 baseline (CR3BP IAS15, #505)

| Quantity | Value |
|---|---|
| Drift over 100 periods | 6.619×10⁻⁹ nd = 0.13 m |
| Drift per period | 6.6×10⁻¹¹ nd = 0.0013 mm |
| Jacobi drift (span) | 4.0×10⁻¹⁰ (< 1×10⁻⁹ IAS15 hygiene) |
| Verdict | BOUNDED, 9 orders of margin inside 3A band |

The orbit is maximally stable (Barden ν = 1.21×10⁻⁶ ≈ 0, the midpoint of the
only stable (3,2) window at this μ). The V2 drift is purely from seed error
3.5×10⁻¹⁰ nd amplified by the finite Floquet multiplier — not residual chaos.

### Real-eph perturbation budget

**Charon eccentricity (the dominant real-eph correction):**

```
Max radial excursion:  δr = a_PC × e_PC = 19 600 km × 5×10⁻⁵ = 0.98 km  (5×10⁻⁵ nd)
```

For a spacecraft at distance r ≈ 1 nd (≈ 19 600 km) from Charon, the
eccentricity induces an oscillatory force perturbation:

```
δF ≈ 2 × μ × δr / r³ = 2 × 0.1088 × 5×10⁻⁵ / 1³ = 1.09×10⁻⁵ nd t⁻²
```

Orbit-averaged displacement amplitude per (3,2) period (T_nd = 11.83):

```
δx_osc ≈ δF × T_nd² / (4π²) = 1.09×10⁻⁵ × 140 / 39.5 ≈ 3.86×10⁻⁵ nd = 0.76 km = 756 m
```

For a resonant orbit the secular accumulation is suppressed (the resonance
angle locks over long timescales), so the eccentricity perturbation creates
**a bounded oscillation of amplitude ≈ 0.76 km** around the CR3BP trajectory
— not a monotonically growing drift.

**Comparison table:**

| Quantity | Value (nd) | Value (physical) |
|---|---|---|
| V2 drift, 100 periods | 6.62×10⁻⁹ | 0.13 m |
| V2 drift per period | 6.6×10⁻¹¹ | 0.0013 mm |
| Eccentricity oscillation amplitude | 3.86×10⁻⁵ | 0.756 km = 756 m |
| **Ratio: eccentricity / V2 (100 periods)** | **~5,800×** | — |
| **Ratio: eccentricity / V2 per-period** | **~580,000×** | — |

The real-eph oscillation amplitude (0.76 km) is 5,800× the V2 100-period
drift (0.13 m). **Real-eph is NOT ≈ V2.**

Other real-eph corrections:
- Inclination of Charon's orbit: negligible (tidally locked, coplanar to high
  precision; no published inclination offset for Charon vs the Pluto equatorial
  plane).
- Nix/Hydra third-body: GM = 1.5×10⁻³ / 2.0×10⁻³ km³/s² → mass fraction
  ~2×10⁻⁶ of the system; three orders of magnitude below the eccentricity term.
- Pluto/Charon J2 asphericity: New Horizons found both bodies nearly spherical;
  Pluto J2 is measured but small; effect at the orbit scale (a = 19 600 km)
  is << eccentricity term.

**Comparison with Uranian V4 (#335):**

| System | Satellite e | δr (km) | δr / a |
|---|---|---|---|
| Charon (this row) | < 5×10⁻⁵ | 0.98 | 5×10⁻⁵ |
| Umbriel (#335, V4) | 0.00411 | 1093 | 4.1×10⁻³ |
| Oberon (#335, V4) | 0.00056 | 327 | 5.6×10⁻⁴ |

Charon's eccentricity is 82× smaller than Umbriel's and 11× smaller than
Oberon's. The Pluto-Charon system is extraordinarily close to the CR3BP
idealization — arguably the most circular known natural satellite orbit.

### Why V3 is inapplicable as a stability gate for this row

The V3 pattern used in #331 (Uranian quasi-cycler) and #335 tested whether
an orbit defined by **encounter geometry in a Keplerian model** survived when
the moon positions were replaced by real SPICE ephemerides. The orbit was
defined by encounter sequence and ΔV, not by exact mathematical periodicity.

The Pluto-Charon (3,2) row is a **CR3BP-periodic orbit**: the orbit closes
exactly (to machine precision) in the CR3BP model. Its validation claim is
"periodic orbit in CR3BP, confirmed stable over 100 periods by an independent
integrator (REBOUND IAS15) in an independent frame (inertial)."

Propagating this CR3BP IC in the real-eph system (slightly eccentric Charon
orbit) would show:
1. An oscillation of amplitude ~0.76 km around the CR3BP trajectory (from
   the e < 5×10⁻⁵ eccentricity perturbation).
2. This oscillation is 5,800× larger than the V2 drift and would dominate any
   "drift after N periods" measurement.
3. The oscillation is NOT evidence of instability — it is a MODEL-MISMATCH
   signal (the CR3BP is an excellent but non-exact idealization of the real
   system).

A **meaningful** real-eph gate for a CR3BP-periodic orbit would require:
(a) Differentially correcting the CR3BP IC in the real-eph model to find the
    real-eph analog of the periodic orbit (a new, harder computation); and
(b) Measuring the stability of THAT orbit.

This is out of scope for this task and would constitute a new scientific
computation, not a V2→V3 promotion. Moreover, the extraordinary circularity
of the Pluto-Charon orbit (e < 5×10⁻⁵) makes it the system where CR3BP is
MOST justified as an accurate model — the argument for V3 "adding
discrimination" is weakest here of any CR3BP system.

---

## V3 recommendation

**V3 NOT APPLICABLE. Row stays V2-ballistic.**

Justification:

1. **Infra-gated (Gate a fails):** No PLU satellite SPK in-repo or GMAT
   install. Kernel is downloadable (plu058.bsp, ~3.8 MB from JPL NAIF) but
   not installed.

2. **Methodology inapplicable (Gate b):** Even if the kernel were installed,
   the standard V3 test (propagate CR3BP IC in real-eph, measure drift) would
   show ~0.76 km amplitude oscillations from the Charon eccentricity
   perturbation — a model-mismatch signal 5,800× larger than the V2 drift.
   This does not add stability discrimination; it adds a different question
   (does the CR3BP IC survive in the real system?). Answering that question
   rigorously requires finding the real-eph analog orbit, a separate
   computation.

3. **V2 is already the strongest applicable gate:** The V2-ballistic evidence
   (100-period REBOUND IAS15 bounded-drift run, #505) used an independent
   integrator in an independent frame — the full independent-cross-check that
   V3 would replicate. For this CR3BP-periodic-orbit row, V2-ballistic IS the
   independent-integrator gate; there is no additional independent-integrator
   step to take before real-eph continuation.

4. **CR3BP is the appropriate model:** Charon's eccentricity (< 5×10⁻⁵) is
   the tightest orbital circularity constraint of any known natural satellite.
   The real Pluto-Charon system is the closest natural realization of the CR3BP
   idealization in the solar system. The argument for V3 "revealing something
   the V2 model misses" is weaker here than anywhere else.

**Decision is the user's.** This note constitutes the full evidence package;
`data/catalogue.yaml` is untouched.

---

## Discipline anchors

- `project_negative_results_registry` — "V3 not applicable" is a real result,
  method-versioned: real-eph V3 is not a valid gate for CR3BP-periodic-orbit
  rows where V2 already consumed the independent-integrator step.
- `feedback_orbit_closure_discipline` — the characterization is based on
  the math (perturbation budget and methodology): no result was fudged to
  obtain a preferred answer.
- `feedback_never_give_up_reproducing_papers` — this task had no published
  real-eph Pluto-Charon orbit to reproduce; the scope assessment IS the task.
- `project_s1l1_realeph_closure_blocker` (analogy) — like S1L1, the real-eph
  blocker here is not infrastructure but physics: the V3 methodology doesn't
  address the orbit's defining model.

---

## Artefacts

- Verdict note: `docs/notes/2026-07-01-506-pluto-charon-v3-scope.md` (this file)
- Evidence tests (unchanged): `tests/search/test_505_pc_v2_longspan.py` (V2 evidence, 4/4 pass)
- No new code committed: the scope assessment is analytic; no real-eph
  propagation was run (infra-gated + methodology inapplicable).
- Ruff: clean (linting + format, checked pre-commit).
