# Rosengren-Ross-Kumar-Rawat 2024 — The resonant structure of xGEO, cislunar domain awareness (AMOS)

**Date filed/digested:** 2026-07-15 (light digest); **full mining pass:** 2026-07-15 (#597/#15).
**Trigger:** same `ross.aoe.vt.edu/papers/` review as the three companions above. Confirmed genuinely
new via title-page read; full 20-page body read page-by-page for this pass — turned out to have
substantially more original content than the "survey/motivation only" characterization in the light
digest.

**Citation:**
> Rosengren, A. J., Ross, S. D., Kumar, B., and Rawat, A. (2024). "The Resonant Structure of xGEO and
> Implications for Cislunar Domain Awareness." Advanced Maui Optical and Space Surveillance
> Technologies Conference (AMOS) 2024.
> No DOI found (AMOS conference paper, not indexed in CrossRef as of this filing).

Registered as `CorpusAnchor` `rosengren-2024-xgeo-resonant-structure-amos` in
`src/cyclerfinder/search/literature_check.py` (provenance `verified-against-source`).

## Revised scope (corrects the light digest's "survey only, no new result" characterization)

The light digest undersold this paper — it DOES contain original quantitative content, not just a
survey framing. Two substantive contributions:

1. A **rigorous, sourced dynamical definition of "xGEO"** via the Laplace radius (Section 2.1.2),
   not just an informal "beyond GEO" label.
2. **Real-spacecraft case studies** (Section 4.2) reproducing decades-long orbital histories of actual
   cataloged objects (TESS, Spektr-R, OGO-1, Luna 3-adjacent) purely from TLE data using `REBOUND`/
   `ASSIST` n-body propagation, cross-validated against JPL Horizons and GMAT.

## The Laplace-radius xGEO definition (Section 2.1.2, Eq. 10-13)

Standard MEO secular-resonance theory (Eq. 1-9) treats the satellite Hamiltonian as
`H = H^Kep + H^Earth-oblateness + H^Moon + H^Sun`, with Earth-oblateness (`J2`) precession dominating
lunisolar secular effects for `a <~ 5 R_Earth` (the classical MEO regime). Beyond GEO altitude, the
lunisolar secular terms become comparable to the oblateness term and eventually dominate — the paper
derives the **Laplace radius** `r_L` as the exact semimajor axis where the lunisolar apsidal+nodal
secular rate equals the oblateness rate:

```
r_L^5 = a^5 * omega_Earth / (omega_Moon + omega_Sun) ~= 7.7 R_Earth   (Eq. 13)
```

where `omega_Moon`, `omega_Sun`, `omega_Earth` are the respective secular apsidal/nodal precession-rate
coefficients (Eq. 8-9, 11-12), each `~ a^{3/2}` scaled by the perturbing body's own orbital elements.
**This gives xGEO a precise, falsifiable numerical boundary (`~7.7 R_Earth`) rather than an informal
label** — beyond `r_L`, atmospheric drag is negligible (even for grazing objects) and solar-radiation-
pressure secular effects are orders of magnitude smaller than lunisolar, even for high area-to-mass
objects (`A/m > 1 m^2/kg`), confirming xGEO is fundamentally governed by third-body (R3BP-style)
dynamics, not perturbed-two-body dynamics. Formalized as: **xGEO = restricted four-body problem
(R4BP)**, approachable either locally (perturbed-Hamiltonian, valid near-Earth) or globally (R3BP
methods, valid throughout).

## Kozai-Lidov formalism (Section 2.1.3, Eq. 14-17)

Restates the classical doubly-averaged CR4BP ("vZLK") Hamiltonian
`H_vZLK ~ -(a^2/8)(mu_Sun/a_Sun^3 + mu_Moon/a_Moon^3) [2+3e^2-(1-e^2+5e^2 sin^2(omega))sin^2(i)]`,
integrable with constants `c1 = (1-e^2)cos^2(i)`, `c2 = e^2[2/5 - (1-c1/(1-e^2))sin^2(omega)]`. Notes
the **critical inclination `i_crit ~= 39.2 deg`**: above this, stationary solutions at
`omega = pi/2, 3pi/2` exist with `omega`-libration (Kozai-Lidov cycles: coupled `e`-`i` oscillation).
Flags that xGEO objects with ecliptic inclination `> 39.2 deg` are expected to show octupole-order KL
effects (chaotic prograde<->retrograde flips) PLUS interaction with lunar orbital/precession
frequencies and MMRs — i.e. KL and MMR effects are NOT independent in this regime, an important
caveat for interpreting any single-mechanism analysis of a real xGEO object's long-term behavior.

## MMR width formula (Section 2.1.4, Eq. 18-19) and the "Kirkwood gaps of cislunar space" framing

Restates the standard MMR resonance condition
`Phi_lmphqj = (l-2p+q)*M_dot - (l-2h+j)*M_moon_dot + secular terms ~= 0`, simplifying (secular rates
`<<` orbital rates) to `a ~= a_{k:k_Moon} = (k/k_Moon)^{2/3} * a_Moon` (Eq. 19) — the same resonance-
center formula used across the whole paper cluster. Explicitly poses **"What are the Kirkwood gaps of
cislunar space?"** as a live open question motivating this whole research program, and notes ESA's
SMART-1 low-thrust lunar transfer deliberately used a **combination of coast arcs and weak lunar-
gravity-assist passages through MMRs** to raise perigee and rotate inclination/perigee-argument for
capture — a concrete existing-mission example of MMR-exploitation for trajectory design (distinct from
IBEX/TESS's passive stable-resonance parking).

## Real-spacecraft case studies (Section 4, Figures 1, 5-10) — the paper's most novel content

- **Fig. 1**: a full historical+current xGEO space-object catalog snapshot (TLE-derived, ~90+ named
  objects: Explorer 14, Vela series, OGO series, IMP series, ISEE, INTEGRAL, Chandra/CXO, XMM-Newton,
  THEMIS, Cluster II, Prognoz, ASTRON, Chang'e, Longjiang, TESS, IBEX, Spektr-R, etc.) in both
  `(a,i)` and `(a,e)` planes, marking the Laplace radius, Earth-reentry curve `a(1-e)=R_Earth`, the
  Moon's Hill-sphere-crossing curve, and MMR/secular-resonance locations. Highlights **Luna 3
  (1959 Theta 1)** — the first spacecraft to circumnavigate the Moon and return to an Earth-crossing
  elliptical orbit, which subsequently crashed after just 11 revolutions due to lunar-perturbation-
  driven eccentricity growth (a satisfactory dynamical treatment only very recently achieved, ref [30]
  = Amato-Malhotra-Sidorenko-Rosengren 2020).
- **IBEX** (2008-051A): nominal 2-year mission orbit became chaotic/unpredictable beyond 2.5 years due
  to previously-uncharacterized lunar perturbations; the extended mission was placed into a **stable
  3:1 lunar MMR** to avoid altitude/eclipse constraint violations.
  **TESS** (2018-038A): established via a lunar swing-by into a **stable 2:1 lunar MMR** (`P_Moon/2`).
  **Spektr-R** (2011-037A): "seemingly" also near a 3:1 MMR, but with high ecliptic inclination that
  additionally permits Kozai-Lidov + other secular-resonant phenomena (not fully characterized) —
  echoed from the resonance-widths companion paper's TLE-validation figure.
- **Fig. 7-8, TESS 30-day and ~2-year `REBOUND` propagations**: comparing a full Earth-Moon-Sun
  4-body model against an Earth-Moon-only 3-body model shows **only the more realistic 4-body model
  reproduces JPL Horizons ephemeris** — the CR3BP (no Sun) is explicitly shown to be inadequate for
  TESS specifically, a useful caveat on the limits of this whole paper cluster's planar-CR3BP-only
  methodology when applied to a REAL (not idealized) spacecraft.
- **Fig. 9, Spektr-R decadal (~12-year) propagation from TLE-derived initial conditions using
  `ASSIST`**: accurately captures the resonant transition visible in the real TLE history, though the
  underlying dynamical DRIVER of that transition is explicitly flagged as "yet to be fully identified
  and understood quantitatively" — an open problem, not resolved by this paper.
- **Fig. 10, OGO-1 (1964-054A) full 56-year recovery**: OGO-1's orbital lifetime was badly
  mispredicted in the 1960s; the catalog has a real ~30-year TLE gap (tracking resumed 2001). Using
  ONLY the initial 1964 launch state + `REBOUND` (full perturbation set + two-body regularization near
  close approach, no further differential correction), the paper **fully recovers OGO-1's entire
  56-year trajectory including its correctly-predicted 29 August 2020 reentry** — a striking practical
  validation of long-horizon N-body propagation for xGEO SDA, achieved via knowledge of the object's
  secular+resonant motion rather than TLE-based orbit determination.

## Relevance to this project

- Lower DIRECT relevance to this project's own cycler-search machinery than the other 3 papers (this
  one contributes no new heteroclinic/manifold method), but the **Laplace-radius xGEO definition**
  (`r_L ~ 7.7 R_Earth`, Eq. 13) is a clean, sourced, quotable boundary if this project ever needs to
  frame a cislunar cycler/tour result relative to the GEO/xGEO regime distinction (e.g. distinguishing
  "this cycler passes through the xGEO belt" from "this cycler stays circumterrestrial").
  Notable: this paper explicitly classifies **Earth-Moon cyclers** as one of the R3BP-governed xGEO
  dynamical regimes its group studies (alongside libration-point orbits and MMRs) — direct textual
  confirmation this whole author cluster considers "cycler" the same family of object this project
  searches for.
- The Kozai-Lidov `i_crit ~= 39.2 deg` threshold and the KL-MMR interaction caveat are useful
  background if this project's inclined/3D cislunar work (any future spatial extension of the interior/
  exterior MMR results in the companion papers) needs to account for high-inclination secular effects
  alongside MMR structure.
- The demonstrated `REBOUND`/`ASSIST` long-horizon N-body validation methodology (OGO-1's 56-year
  recovery from launch conditions alone) is a useful precedent if this project ever needs to validate
  a long-duration real-ephemeris cycler propagation against a genuinely independent, decades-scale
  historical case.

**Status: full mining pass complete.** No catalogue impact, no writeback (methods/results reference).
`CorpusAnchor` registered.
