# #588: final adjudication of the 20 unmatched Sun-Earth ER3BP clusters

**Date:** 2026-07-14
**Pipeline:** #583/#586 sweep (264 raw candidates) -> dedup (45 clusters) -> known-family cross-check
(20 unmatched) -> live literature-check (2 Gurfil-Kasdin companion papers digested, neither closes
the gap) -> Opus adjudication -> Fable adversarial second opinion (this note consolidates both).

## Final disposition

**14 of 20 clusters: confidently not novel.**

- **8 large planar DEO clusters** (ids 1, 2, 3, 7, 20, 35, 36, 44) sit within 1% of Henon's
  foundational Sun-Earth ẏ₀≈-2x₀ retrograde-satellite curve (the same curve Gurfil-Kasdin's own
  Families A/B/F sit on), at radii (5.4M-40.8M km) beyond where the original paper's 12
  characterization sets happened to land. Curve extension, not a new family.
- **6 small 3D DRO clusters** (ids 10, 12, 14, 29, 37, 38, radii 80k-560k km): an Opus pass initially
  dismissed these via a `v0/vcirc≈1` ("Sun-negligible, trivially Keplerian") criterion. **A Fable
  adversarial check found this criterion invalid** — applied to the known Gurfil-Kasdin families
  themselves, it scores Family F at 0.91 and Family K at 1.30, i.e. it would dismiss or nearly
  dismiss published, real families. The metric is radius in disguise (a pure monotone function of r
  for any retrograde orbit) and carries no independent discriminating power. **The final verdict
  (not novel) survives, but on the correct grounds Fable substituted:** the Sun-Earth ER3BP model
  itself is invalid at these radii — several sit at or below the Moon's 384k km orbit, where the
  Moon (not the Sun) is the dominant perturbation the model omits — combined with continuum-
  membership logic (retrograde-satellite structure is continuous down to the primary; smaller-radius
  members are expected, not novel). This distinction matters for future work: **the `v0/vcirc`
  heuristic must not be reused as a triage axis** (see
  [[feedback_verify_metric_semantics_before_ranking]], third instance logged).

**6 of 20 clusters warrant a bounded follow-up before dismissal is treated as final** — NOT a
novelty claim, an open question:

- **Primary (3): clusters 40, 42, 43** — genuinely 3D (z0 up to ~3.2M km), radii 4.1M-13.4M km
  bracketing Family J's 6.9M-8.5M km band, ẏ₀/x₀ = -2.03 to -2.08 (close to J's -1.99). Leading
  hypothesis is "J extended," but each is a **single unreplicated GA hit** (duplicate_count=1, vs.
  2-25 for every dismissed cluster) with a growing 5yr r-band (trend 0.066-0.107) — "marginal GA
  artifact" is at least as plausible as "real family member" from current evidence.
- **Secondary (3): clusters 24, 25, 39** — near-planar, radius-band-adjacent to Families B/I/C.
  Cluster 39 was originally dismissed on a factually wrong claim (Opus stated it matches a DRO band;
  it actually sits in the gap between Families C and B, at r/r_Hill≈1.0, also single-hit) — Fable
  flagged this as a real citation/reasoning error and recommended 39 ride along with 24/25 rather
  than stay in the confident-dismissal bucket.

## Recommended follow-up (registered as #590, not yet dispatched)

1. **Targeted literature search for the corrected corpus.** Fable's independent search found the
   actually-relevant prior art for large-radius inclined 3D Sun-Earth structures (4-13M km ≈ 3-7
   Hill radii) is the **classical quasi-satellite orbit (QSO) literature**, not the 2025-2026
   ER3BP-bifurcation papers Opus's transcript cited (which were never grounded in a note per
   [[feedback_ground_citations_against_content]] and are topic-adjacent — planar or cis-lunar — not
   on point). Candidates to check: Mikkola et al. 2006 (MNRAS, "Stability limits for the
   quasi-satellite orbit"), Lidov & Vashkov'yak 1993/1994, Sidorenko et al. 2014, Pousse, Robutel &
   Vienne 2016/2017 (CMDA, quasi-satellite motion including the ER3BP). None of these DOIs are
   independently CrossRef-verified yet — that must happen before any is cited as a corpus anchor.
2. **Cheap numerical continuation.** A single predictor-corrector chain in orbit size connecting
   40/42/43 to Family J would settle "same curve, just unsampled" vs. "distinct bifurcated branch"
   directly — no GA sweep needed.
3. **Longer-horizon stability re-check.** Re-run the drift classifier at a longer horizon (e.g.
   100 revolutions) for 40 and 42 specifically, whose 5yr r-bands already show the largest growth
   trend among survivors.
4. Optionally phase-normalize the 6 dismissed small-radius clusters to a common osculating section
   (they currently sit at essentially arbitrary orbit phase, x0≈0) to make "curve extension" a proven
   claim rather than a plausible one — lower priority, not blocking anything.

## Scope reminder (Fable's framing correction, important for future readers)

**None of these 20 clusters — nor any of the 14 published Gurfil-Kasdin families — are
catalogue-admissible** under this project's own scope (cycler / quasi_cycler / precursor_mga /
mga_tour, see [[project_catalogue_scope_expanded_2026-06-15]]). Geocentric DROs and quasi-satellite
orbits are neither cyclers nor tours by themselves. #588's question was purely "new to the published
Sun-Earth orbit-family literature record" — a literature-completeness question, not a
catalogue-discovery one — unless a genuine link into an actual cycler/precursor structure is found
downstream (matching the precursor scope rule: precursors must `inserts_into` a cycler). Nothing here
changes the catalogue; closing #588 does not add rows.

## #588 status: CLOSED (adjudication complete). #590 registered for the bounded 6-candidate follow-up.
