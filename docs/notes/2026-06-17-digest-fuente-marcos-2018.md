# de la Fuente Marcos & de la Fuente Marcos 2018 — Geometric characterization of the Arjuna orbital domain

**Deep-read verdict note, 2026-06-17 AET.** Supersedes the file-only entry in
`2026-06-17-mars-cycler-wave-digest.md`.

## Header

- **Title:** *Geometric characterization of the Arjuna orbital domain*
- **Authors:** C. de la Fuente Marcos, R. de la Fuente Marcos (Universidad
  Complutense de Madrid)
- **Venue:** Astronomische Nachrichten, 11 August 2018; manuscript dated
  Oct 2014, accepted 2014, online 2015-02-09. arXiv:1410.4104v2 [astro-ph.EP].
  Note: the cited venue + arXiv-id labels are consistent (publisher
  metadata) but the on-page date is 11 August 2018 with a 2014 receipt
  history — this is a single paper, the survey's "2018" label is correct.
- **Pages:** 19 (excluding references continuation)

## What the paper actually is

A geometric / Monte Carlo characterisation of the **Arjuna class** of
near-Earth asteroids — small bodies (10-60 m) on Earth-like orbits
(0.985 < a < 1.013 AU, 0 < e < 0.1, 0 < i < 8.56°). These objects undergo
repeated 1:1 mean-motion resonance trappings with Earth and can become
**temporary Trojans, horseshoe librators, quasi-satellites, and even
transient natural satellites** ("minimoons"). The 1:1 commensurability with
Earth's orbit is the defining invariant.

**13 known objects** (Table 1, page 8) with their Keplerian elements,
MOIDs (Minimum Orbit Intersection Distances), and Δv-to-rendezvous values
from echo.jpl.nasa.gov/~lance/delta_v/delta_v_rendezvous.html:

| Object         | a (AU)   | e        | i (°)   | MOID (AU) | Δv (km/s) |
|----------------|----------|----------|---------|-----------|-----------|
| 2003 YN107     | 0.988697 | 0.013988 | 4.32115 | 0.004436  | 4.879     |
| 2006 JY26      | 1.010021 | 0.083042 | 1.43927 | 0.000052  | 4.364     |
| **2006 RH120** | 0.998625 | 0.019833 | 1.52613 | 0.000679  | **3.820** |
| 2008 KT        | 1.010800 | 0.084087 | 1.98426 | 0.000705  | 4.425     |
| 2008 UC202     | 1.010179 | 0.068692 | 7.45350 | 0.002890  | 5.471     |
| **2009 BD**    | 1.008614 | 0.040818 | 0.38516 | 0.003565  | **3.870** |
| 2009 SH2       | 0.991719 | 0.094175 | 6.81073 | 0.000390  | 5.070     |
| 2010 HW20      | 1.010924 | 0.050111 | 8.18503 | 0.008984  | 5.690     |
| 2012 FC71      | 0.989158 | 0.088009 | 4.94929 | 0.056850  | 4.686     |
| 2012 LA11      | 0.998946 | 0.096338 | 5.10630 | 0.007504  | 4.751     |
| 2013 BS45      | 0.993858 | 0.083801 | 0.77530 | 0.011514  | 4.083     |
| 2013 RZ53      | 0.999722 | 0.048260 | 1.50653 | 0.001805  | 4.198     |
| 2014 EK24      | 1.003690 | 0.072150 | 4.72541 | 0.033854  | 4.846     |

Δv to rendezvous ranges 3.82-5.69 km/s. **2006 RH120 (3.82 km/s) and 2009
BD (3.87 km/s)** are the cheapest — both noted explicitly as cheaper than
LEO-to-Moon (6.0 km/s). 2006 RH120 was actually a transient Earth natural
satellite Jul 2006-Jul 2007 (Kwiatkowski et al. 2009).

Key findings:

- ≈8% of Arjunas have v_rel² < v_esc² at perigee → probability of capture
  as transient minimoon is ~8%.
- Gravitational focusing factor ≥ 1000 for 0.053% of Arjunas; > 10-100 for
  6.7% → 10-1000-fold increase in impact cross-section vs Apollo/Aten
  populations.
- Estimated population: ≥ 172 objects > 30 m (conservative); likely 10³ at
  metric-sized, possibly 10⁴ if synodic-period distribution is biased
  toward long synodic periods.
- Most Arjunas escape detection: long synodic periods (43+ yr), small sizes
  (mostly H > 22 mag), solar-elongation-at-perigee biased toward
  daytime sky (<90° from Sun for ~50%). Gaia will NOT significantly
  improve detection rate. Space-based surveys needed.

## Catalogue / KNOWN_CORPUS relevance

**Directly relevant for #308 (asteroid-leveraging cycler search) and the
Adamo NEA-supply branch.** This paper provides:

1. **A sourced asteroid candidate list** (13 known Arjunas with full
   Keplerian elements at JPL epoch 2456800.5 = 2014-May-23 0h UT) — could
   be ingested as a `KNOWN_CORPUS` data anchor or as input to a
   Lambert-targeting scan from a candidate `quasi_cycler` heliocentric
   waypoint orbit.
2. **Δv-rendezvous economics** (3.82-5.69 km/s) — the canonical reference
   for "asteroid-leveraging is cheaper than LEO-to-Moon" claims that
   Adamo 2025's slide 16 alludes to ("493 catalogued NEAs orbit Sun
   between Earth and Mars"). Adamo's worked example (1996 XB27) is NOT
   in this Arjuna list (1996 XB27 is between Earth and Mars at a ~1.2 AU,
   well outside Arjuna's a < 1.013 AU constraint) — so Adamo's NEAs are
   the broader 493-NEA population, not specifically the Arjuna sub-class.
3. **Capture mechanism** — Kozai resonance + 1:1 mean-motion + close E-M
   encounters can drive transitions between Trojan / horseshoe /
   quasi-satellite / minimoon states. This is precisely the dynamical
   regime where "asteroid-leveraging cyclers" would operate (sticking-and-
   releasing from the 1:1 resonance manifold).

**No new KNOWN_CORPUS anchor for "Arjuna cyclers"** because the paper does
**not** propose a cycler architecture — it's an *observational/Monte Carlo
characterisation* of natural objects. No trajectory design, no Lambert
scan, no Δv tour planning beyond the single-target rendezvous values.

However, the paper's dynamical-class definition (1:1 Earth co-orbital,
low-e low-i, transient resonance) could be a **methodology citation**
for any future Arjuna-based `quasi_cycler` row — analogous to how Antoniadou-
Voyatzis 2018 is a spatial-CR3BP citation anchor.

### Possible KNOWN_CORPUS extension (deferred)

A note for when #308 reactivates: the de la Fuente Marcos pair have a
**3-paper series** (2013a, 2014a/b/c) cited in this paper as
``de la Fuente Marcos & de la Fuente Marcos 2013a, MNRAS, 434, L1`` and
2014 sequels — the 2013a paper appears to be the *first* Arjuna-class
classification paper. Acquisition of 2013a (MNRAS 434:L1) would be the
natural primary anchor; the 2018/2014 Astron. Nachr. paper digested here
is the *characterisation* follow-up.

## Errata vs the pre-read survey

The pre-read survey filed this as "*Arjuna-class Earth co-orbital NEOs.
Directly relevant for #308 asteroid-leveraging cycler search. Tagged for
#308 reactivation context if it ever re-fires.*" The deep-read **confirms**
this verdict with two clarifications:

1. The paper is observational/dynamical-class characterisation, NOT a
   cycler-architecture proposal. It provides target data (Table 1) and a
   dynamical-class definition, not a trajectory family.
2. The natural anchor for any future Arjuna-class catalogue work is the
   **de la Fuente Marcos & de la Fuente Marcos 2013a (MNRAS 434:L1)**
   paper, of which the 2018 Astron. Nachr. paper is the
   characterisation sequel. Acquisition wishlist update.

## Action items for parent

1. **No new KNOWN_CORPUS anchor.** This paper is reference/data, not a
   cycler architecture.
2. **Acquisition wishlist:** add the de la Fuente Marcos **2013a (MNRAS
   434:L1)** paper. It is the *original* Arjuna-class definition paper
   and would be the natural KNOWN_CORPUS primary citation if #308 ever
   produces an Arjuna-leveraging cycler.
3. **#308 reactivation data:** Table 1 (13 known Arjunas with full
   Keplerian elements + Δv-to-rendezvous + MOID) is the sourced input
   list for any Lambert scan from a candidate Sun-Earth Arjuna-leveraging
   `quasi_cycler` orbit. Cheapest targets: 2006 RH120 (Δv = 3.82 km/s,
   transient minimoon Jul 2006-Jul 2007) and 2009 BD (Δv = 3.87 km/s).
4. **Cross-link to Adamo 2025:** Adamo's "493 NEAs orbit Sun between Earth
   and Mars" claim is the broader NEA-population statistic; this paper's
   13-object Arjuna list is the *Earth-co-orbital subset* (a ≈ 1 AU, not
   between Earth and Mars). The two populations are disjoint by definition.
   No correction needed — both are valid for different NEA-supply use cases.
5. **Population estimate for forward search:** ≥ 172 objects > 30 m,
   likely ≥ 1000 metric-sized → the unobserved Arjuna population is
   substantial. Any forward Arjuna-leveraging search should be epoch-
   parametric (population evolves on 10-100 yr e-folding times) and
   should not assume the 13-object table is exhaustive.
