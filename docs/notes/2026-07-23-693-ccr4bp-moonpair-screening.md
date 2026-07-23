# #693 — CCR4BP cross-solar-system moon-pair screening pass (2026-07-23)

Analysis-only (no code, no catalogue writes, no new corrector/EOM/manifold search), mirroring
`#679`'s and `#686`'s discovery-strategy-pass format. Purpose, in the coordinating session's own
framing: the `#688`→`#689`→`#690`→`#691` build chain produced a full CCR4BP (EOM+STM, torus
corrector, whisker diagnostic) capability, positive-controlled at Jupiter-Europa-Ganymede
specifically because that pair already has published torus results (Kumar-Anderson-de la
Llave-Gunter 2021, arXiv:2109.14815) — so a heteroclinic-connection discovery run there would NOT
be novel. Before committing to any further multi-day per-system CCR4BP build, screen cheaply
across every other same-primary moon pair in the outer solar system on two independent axes:
**tractability** (does the real geometry satisfy the model's circular/coplanar/weak-forcing
assumptions) and **novelty** (has anyone already published CCR4BP/whiskered-torus/heteroclinic
work on that specific pair).

## Inputs read

`src/cyclerfinder/core/ccr4bp.py` (full docstring + `CCR4BPSystem` fields — the model's own
stated assumptions and the Europa/Ganymede default derivation); `#689`/`#690`/`#691`'s full
`data/OUTSTANDING.md` bullets (the corrector's documented weak-forcing near-degeneracy caveat at
`mu_gan=7.8e-5`, and the eccentric-orbit-vs-circular-moon distinction); `src/cyclerfinder/core/
satellites.py` (full registry — confirmed it carries GM/radius/sma/derived mean-motion only, NO
eccentricity or inclination for any moon, so those had to be sourced externally); `search/
literature_check.py`'s full `KNOWN_CORPUS` (2852 lines, read via targeted grep + two full-section
reads) — found two directly load-bearing existing anchors (Kumar Uranus-Oberon/Titania-Oberon
CCR4BP secondary-resonance work; Canales-Howell-Fantino Titania-Oberon one-shot halo transfer,
NOT CCR4BP); `docs/notes/2026-06-16-328-uranian-cycler-lit-deep-dive.md` (the source digest behind
the Kumar Uranus anchor, read in full for what it actually established vs. inherited); `docs/
notes/2026-07-23-digest-kumar-2021-europa-ganymede-ccr4bp-resonant-orbits.md` and the companion
2023 secondary-resonance digest (grepped for scope — neither mentions Io or any pair beyond
Europa/Ganymede); `docs/notes/CORPUS_INDEX.md` (grepped for CCR4BP/whisker/heteroclinic entries;
also confirmed `#679`/`#686`'s own strategy notes are NOT registered there, so this note isn't
either, matching precedent); `ls /Users/bruce/dev/cyclers_pdf/papers/` (read access confirmed,
full directory listing scanned for Uranus/Saturn/Neptune CCR4BP titles — none found beyond the
three Jovian Kumar papers already in corpus). External sourcing: JPL SSD "Mean Orbital Elements of
Regular Planetary Satellites" (`https://ssd.jpl.nasa.gov/sats/elem/`, fetched 2026-07-23 — the
sma/e/i/period table for every Jupiter/Saturn/Uranus/Neptune regular moon below) and JPL SSD
"Planetary Satellite Physical Parameters" (`https://ssd.jpl.nasa.gov/sats/phys_par/`, fetched
2026-07-23 — GM for Neptune's small moons not carried in-repo); six live WebSearch queries for
CCR4BP/whiskered-torus/heteroclinic hits on specific candidate pairs (queries and results below).
**Honest limitation, stated explicitly per the task's own mandate:** WebSearch/WebFetch WERE
available in this sandbox (contrary to my working assumption going in) and were used for six
targeted queries; this is a real but partial literature check — not every venue is indexed, and
"no hit in six queries" is necessary-not-sufficient for novelty, exactly as `search/
literature_check.py`'s own docstring discipline states.

## 1. Model assumptions actually used as the tractability bar (grounded, not asserted)

`core/ccr4bp.py`'s own docstring states the CCR4BP idealization as: two moons on CONCENTRIC
CIRCULAR COPLANAR orbits about the primary, no mutual moon-moon force, time-periodic in the
Europa-synodic frame at the perturber's synodic rate REGARDLESS of exact commensurability (the
`#689`/`#690` build already used Europa:Ganymede's actual PHYSICAL 2.014 ratio, not an idealized
exact 2:1 — this is the load-bearing precedent for how "close to a low-integer ratio" a candidate
needs to be: close enough that a base resonant periodic orbit exists to seed the torus corrector
from, not exactly commensurate).

**Reference bar, from the built and positive-controlled system itself (not asserted, read off
`#689`'s own default derivation + JPL SSD):**

| Quantity | Europa (base) | Ganymede (perturber) |
|---|---|---|
| Real eccentricity | 0.009 | 0.001 |
| Real inclination (to Jupiter equator) | 0.5° | 0.2° |
| `mu` / `mu_gan` (perturber mass ratio) | `mu`=2.528e-5 | `mu_gan`=7.805e-5 |
| Period ratio (physical) | — | 2.0297 (JPL; `#689` used 2.014, a slightly different physical epoch/fit) |

`#690`'s own bullet documents the forcing as "very weak" even at `mu_gan=7.8e-5` — the corrector
needed minimal theta1 harmonics (n1=1) and strong gauge/rho weights to avoid the trust-region
wandering into spurious O(1) structure. **This is the single most load-bearing quantitative fact
for this screen**: `mu_gan=7.8e-5` is not a comfortable working point, it is close to the
documented edge of what the corrector handles cleanly. Any candidate pair with a perturber `mu`
an order of magnitude or more BELOW 7.8e-5 inherits that risk untested and unmitigated.

Sanity filter applied before the numeric table (per the dispatch's own criteria): both moons
prograde, both e ≲ 0.1, mutual inclination ≲ 5°.

## 2. Numeric survey — all candidate pairs

GM/e/i/P sourced from JPL SSD (`sats/elem/`, `sats/phys_par/`, both fetched 2026-07-23); GM
cross-checked against `core/satellites.py`'s in-repo registry where present (exact match on every
overlapping value — Io/Europa/Ganymede/Callisto/Mimas/Enceladus/Tethys/Dione/Titan/Hyperion/
Miranda/Ariel/Umbriel/Titania/Oberon/Triton/Proteus GMs all identical to the repo's own JPL-SSD-
sourced constants, confirming both pulls hit the same upstream table). `mu_base` = inner moon's
mass ratio to (primary+inner); `mu_pert` = outer moon's mass ratio to (primary+inner) — i.e. the
exact `core/ccr4bp.py` convention (`mu`, `mu_gan`), inner moon plays Europa's role, outer plays
Ganymede's. Period ratio = P_outer / P_inner.

### Jupiter (beyond the built Europa-Ganymede pair)

| Pair (in→out) | mu_base | mu_pert | P ratio | e_in / e_out | Δi (deg) |
|---|---|---|---|---|---|
| Europa→Ganymede (reference, built) | 2.528e-5 | 7.805e-5 | 2.030 | 0.009 / 0.001 | 0.3 |
| **Io→Europa** | 4.704e-5 | 2.528e-5 | **2.0000** | 0.004 / 0.009 | 0.5 |
| **Io→Ganymede** | 4.704e-5 | 7.805e-5 | 4.059 | 0.004 / 0.001 | 0.2 |
| Europa→Callisto | 2.528e-5 | 5.667e-5 | 4.734 | 0.009 / 0.007 | 0.2 |
| Ganymede→Callisto | 7.804e-5 | 5.667e-5 | 2.333 (≈7:3) | 0.001 / 0.007 | 0.1 |
| Io→Callisto | 4.704e-5 | 5.667e-5 | 9.469 | 0.004 / 0.007 | 0.3 |

Io:Europa's 2.0000 ratio is the tightest, most precisely librating resonance in the solar system
(the defining member of the Io-Europa-Ganymede 4:2:1 Laplace chain) — geometrically the cleanest
CCR4BP candidate surveyed, cleaner than the built Europa-Ganymede reference itself. Ganymede:
Callisto's 2.333 ≈ 7:3 is a REAL near-resonance, independently corroborated live (arXiv:2607.03505,
"The recent crossing of the 7:3 resonance between Ganymede and Callisto," surfaced by WebSearch —
not previously known to this project). All six pairs pass the e/Δi sanity filter comfortably
(every value ≤ the built reference's own Europa e=0.009 / Δi=0.3° except Callisto's e=0.007, still
under).

### Saturn

| Pair (in→out) | mu_base | mu_pert | P ratio | e_in / e_out | Δi (deg) |
|---|---|---|---|---|---|
| Mimas→Tethys | 6.60e-8 | 1.086e-6 | 2.003 (real 4:2 incl. resonance) | 0.020 / 0.001 | 0.5 |
| Enceladus→Dione | 1.90e-7 | 1.928e-6 | 1.997 (real 2:1, drives tidal heating) | 0.005 / 0.002 | 0.0 |
| Tethys→Dione | 1.086e-6 | 1.928e-6 | 1.450 | 0.001 / 0.002 | 1.1 |
| Titan→Hyperion | 2.366e-4 | 9.77e-9 | 1.334 (real, famous 4:3) | 0.029 / 0.105 | 0.3 |
| Mimas→Enceladus | 6.60e-8 | 1.90e-7 | 1.454 | 0.020 / 0.005 | 1.6 |

**Every Saturn `mu_pert` is 40x (Enceladus→Dione) to 8,000x (Titan→Hyperion) weaker than the
already-documented-as-"very weak" JEG `mu_gan=7.8e-5`.** Titan→Hyperion has a strong `mu_base`
(2.37e-4, actually stronger than Europa's) but the perturber term is 5 orders of magnitude below
Ganymede's — the model would collapse to an essentially-unperturbed Saturn-Titan CR3BP with
negligible 4-body content, plus Hyperion's real orbit is the most eccentric candidate surveyed
(e=0.105, an order above the reference bar) and Hyperion's famous chaotic ROTATION (a distinct,
attitude-not-orbit phenomenon, doesn't directly break a point-mass model but is a documented
oddity of the system). This is a quantified version of `#686`'s qualitative Saturn rejection, not
a re-assertion of it: every pair genuinely fails on `mu_pert` alone, by 1.5-4 orders of magnitude
below a forcing level already flagged as barely-tractable.

### Uranus

`i` values below are measured relative to Uranus's own equatorial (Laplace) plane, per the JPL SSD
source convention — this directly confirms the task's framing question: Uranus's ~98° axial tilt
relative to the ECLIPTIC does NOT enter this table at all, because CCR4BP only cares about the
moons' mutual geometry, and every regular Uranian moon orbits within a few degrees of Uranus's own
equator (Δi between any two of the five ≤ 0.1° except Miranda, whose own inclination to that same
equatorial plane is 4.4° — a real, independently-known dynamical anomaly of Miranda specifically,
not an artifact of the tilt).

| Pair (in→out) | mu_base | mu_pert | P ratio | e_in / e_out | Δi (deg) |
|---|---|---|---|---|---|
| Miranda→Ariel | 7.42e-7 | 1.441e-5 | 1.783 | 0.001 / 0.001 | **4.4** |
| Ariel→Umbriel | 1.441e-5 | 1.469e-5 | 1.644 (~5:3) | 0.001 / 0.004 | 0.1 |
| **Umbriel→Titania** | 1.469e-5 | **3.916e-5** | 2.101 (~2:1) | 0.004 / 0.002 | **0.0** |
| Titania→Oberon | 3.916e-5 | 3.543e-5 | 1.547 (~3:2) | 0.002 / 0.002 | 0.0 |
| Ariel→Titania | 1.441e-5 | 3.916e-5 | 3.454 | 0.001 / 0.002 | 0.1 |
| Umbriel→Oberon | 1.469e-5 | 3.543e-5 | 3.249 | 0.004 / 0.002 | 0.1 |

**Umbriel→Titania is, by these numbers, the single best-conditioned candidate outside the built
Jovian system**: `mu_pert=3.92e-5` (only 2x below the JEG reference, comfortably inside a
plausibly-tractable range), both eccentricities below the reference bar, Δi effectively zero, and
a 2.10 period ratio not much looser than JEG's own 2.03. This directly refines `#686`'s blanket
"no present-day MMR at Uranus → no CCR4BP lane" verdict: there is no EXACT commensurability at
Uranus (correct, and well-known — none of the five classical moons are presently locked), but the
model does not require one (the `#690` physical-2.014-ratio precedent), and Umbriel:Titania's
2.10 is a near-resonance with better mass conditioning than several of the JEG-adjacent Jovian
options above. Titania→Oberon looks similarly strong on paper but is DISQUALIFIED on novelty
grounds (§3).

### Neptune

| Pair | mu_base | mu_pert | P ratio | e_in / e_out | Δi (deg) |
|---|---|---|---|---|---|
| Despina→Galatea | 1.71e-8 | 2.78e-8 | 1.281 | 0.000 / 0.000 | 0.0 |
| Galatea→Larissa | 2.78e-8 | 3.73e-8 | 1.295 | 0.000 / 0.001 | 0.2 |
| Larissa→Proteus | 3.73e-8 | 3.78e-7 | 2.022 | 0.001 / 0.000 | 0.2 |

**Triton, explicitly checked and disqualified, not silently skipped:** `core/satellites.py`'s own
registry (citing Jacobson 2009, AJ 137:4322) records Triton's inclination to Neptune's Laplace
plane as ~156.885°; the independently-fetched JPL SSD `sats/elem/` table gives 157.3° (same
source lineage, 0.4° apart — consistent). `i > 90°` is retrograde by convention: Triton orbits
opposite the sense of every other regular Neptunian moon and opposite Neptune's own rotation. This
is not a marginal Δi — it is off the coplanar-prograde premise by ~157° against a ≲5° sanity
filter, i.e. Triton fails the model's basic geometric premise by more than 30x the filter's own
threshold. No CCR4BP pairing involving Triton is geometrically meaningful; the "circular coplanar
corotating" idealization simply does not describe this orbit.

The remaining Neptune pairs (Despina/Galatea/Larissa/Proteus — Nereid is omitted from
`core/satellites.py` itself, GM undetermined per the registry's own comment) are all geometrically
clean (circular, coplanar, and Larissa→Proteus even has a genuine 2.02 near-2:1 ratio) but every
`mu_pert` is 200x (Larissa→Proteus) to 2,800x (Despina→Galatea) weaker than the JEG reference —
worse than every Saturn pair. Neptune has no tractable CCR4BP lane at all: Triton fails on
geometry, everything else fails on mass.

## 3. Novelty check (per pair that cleared §2's tractability filter)

Method, per pair: (1) grep `search/literature_check.py`'s full `KNOWN_CORPUS` for a same-body-set
anchor; (2) grep `docs/notes/CORPUS_INDEX.md` + relevant digest notes; (3) `ls` the private
`cyclers_pdf/papers/` corpus for a matching title; (4) live WebSearch, since (contrary to this
task's own working assumption) WebSearch/WebFetch turned out to be available in this sandbox.
`search/literature_check.py`'s own `check_literature()` engine was NOT invoked directly — it's
built to consume a `CandidateSignature` for an actual discovered trajectory (V∞/resonance/n_rev
per encounter), which doesn't exist yet for any of these (screening precedes any search run); its
`KNOWN_CORPUS` anchors were grepped directly instead, which is the same underlying data the engine
would match against.

| Pair | In-repo anchor? | Live search result | Verdict |
|---|---|---|---|
| Io-Europa | None found | 3 targeted queries, all surfaced only the Europa-Ganymede papers (2109.14815, 2309.06073) and unrelated Galilean-capture papers; no Io-Europa CCR4BP/torus hit | **Not found** (necessary-not-sufficient) |
| Io-Ganymede | None found | 2 targeted queries, same result — no Io-Ganymede-specific hit | **Not found** |
| Europa-Callisto | None found | Surfaced only the Ganymede-Callisto hit below, not Europa-Callisto specifically | **Not found**, but only lightly checked — deserves its own dedicated query before any build |
| **Ganymede-Callisto** | None in `KNOWN_CORPUS` | Surfaced **Aryan & Fitzgerald, "Four Body Invariant Structures and Chaos Analysis for Jovian Multi-Moon Ballistic Transfers," AAS 24-103 (2024)** — a "Planar Concentric Circular Restricted Four Body Problem" (PCCFBP, same model class under a different acronym) computing quasi-periodic invariant tori for BOTH Jupiter-Europa-Ganymede AND **Jupiter-Callisto-Ganymede**, with invariant-manifold transit connections assessed for both | **PUBLISHED — disqualified.** Not previously in this project's corpus or `KNOWN_CORPUS`; flagged as a genuine new find, not yet acquired/digested/indexed (out of this screening task's scope; recommend a follow-up acquisition task) |
| Umbriel-Titania | None found | 1 targeted query; surfaced only Ariel-Umbriel resonance-CAPTURE/migration papers (different subfield) and the Titania-Oberon anchor below, no Umbriel-Titania-specific hit | **Not found** |
| **Titania-Oberon** | **Yes** — `literature_check.py`'s existing "Kumar Uranus-Oberon PCRTBP MMR study (2025)" anchor (arXiv:2509.03655), added 2026-06-16 during `#328`'s Uranian lit deep-dive | Anchor's own citation text: "Section 6.2 studies Uranus-Oberon PCRTBP ... plus heteroclinic connections; **extends to Uranus-Titania-Oberon CCR4BP secondary resonances**." The `#328` digest note gives the actual sub-harmonic ratios (25/69, 21/58, 17/47, 30/83, 13/36, 22/61, 9/25) | **PUBLISHED — disqualified.** Honest caveat: this anchor's `provenance` field is the DEFAULT `"inherited-unverified"` (no explicit override in the anchor definition) — the claim traces to a live-search abstract-level read during `#328`, not a page-by-page source read of arXiv:2509.03655 itself. Per `can_anchor_decision()`'s own discipline this citation is NOT decision-grade for a catalogue promotion, but it is more than sufficient to disqualify a NEW build target at the screening stage — the honest risk is scoop/overlap, not a false negative |
| Ariel-Umbriel | None found (CCR4BP-specific) | Surfaced real, active literature — but in a DIFFERENT subfield: planetary-formation resonance-capture/migration dynamics (arXiv:2509.24631 "Capture and escape from the 2:1 resonance between Ariel and Umbriel in a fast-migration scenario"; arXiv:2305.08794 on the 5:3 passage) | **Adjacent-field "not found," not a direct hit** — no CCR4BP/whiskered-torus/astrodynamics-trajectory-design paper surfaced, but this is a demonstrably active research pair from a different angle, which is a real scoop-context risk even though it doesn't disqualify outright |
| Saturn / Neptune pairs | None found (any) | 1 targeted Saturn query — explicitly returned "no CCR4BP/whiskered torus papers for Saturn moons" | **Not found**, but moot — none clear §2's tractability filter |

## 4. Ranked shortlist

Ranking combines §2's quantitative tractability with §3's novelty verdict; PUBLISHED pairs are
listed for completeness but ranked out of contention.

1. **Jupiter Io-Europa** — tractability: EXCELLENT (exact 2.000 ratio, the tightest resonance in
   the solar system; `mu_pert=2.53e-5`, comparable order to the base leg of the already-validated
   JEG system; e/Δi both under the reference bar). Novelty: not found, three targeted checks.
   **Top recommendation for the next per-pair build, pending user GO.**
2. **Jupiter Io-Ganymede** — tractability: EXCELLENT (`mu_pert=7.80e-5`, essentially IDENTICAL
   forcing strength to the already-validated JEG perturber term, since both are Ganymede/Jupiter-
   dominated; 4.06 ratio is a real Laplace-chain consequence; e/Δi both under bar). Novelty: not
   found, two targeted checks. **Second recommendation** — arguably even lower build risk than #1
   since the perturber-forcing regime is literally already proven to converge in `#690`.
3. **Uranus Umbriel-Titania** — tractability: the best-conditioned NON-Jovian candidate found
   (`mu_pert=3.92e-5`, only 2x below JEG; e/Δi excellent). Novelty: not found, one targeted check,
   but sits in a Uranian neighborhood (Titania-Oberon, Ariel-Umbriel) with genuinely active
   published/adjacent work — higher scoop-context risk than the Jovian options, and a second,
   more targeted lit pass is warranted before any build commitment.
4. **Jupiter Europa-Callisto** — tractability: GOOD (`mu_pert=5.67e-5`, comparable to JEG) but no
   clean low-integer commensurability (4.73 ratio). Novelty: only lightly checked (piggybacked on
   the Ganymede-Callisto query); needs its own dedicated search before ranking above #3.
5. **Jupiter Ganymede-Callisto** — tractability GOOD, but **DISQUALIFIED**: published (Aryan &
   Fitzgerald AAS 24-103, 2024) — quasi-periodic tori + manifold-transit connections already
   computed for this exact pair.
6. **Uranus Titania-Oberon** — tractability EXCELLENT by the numbers (best mass balance and
   eccentricities of any pair surveyed), but **DISQUALIFIED**: published (Kumar arXiv:2509.03655
   Section 6.2 + companion secondary-resonance work), per this project's own pre-existing
   `literature_check.py` anchor.
7. **Uranus Ariel-Umbriel** — tractability GOOD (`mu_pert=1.47e-5`, ~5x below JEG, still plausibly
   workable) but real adjacent-field literature exists (resonance-capture/migration dynamics);
   ranked below Umbriel-Titania on scoop-context risk alone, not disqualified.
8. **Everything else surveyed (all Saturn pairs, all Neptune pairs, Uranus Miranda-Ariel, Jupiter
   Io-Callisto) — REJECTED on tractability, not novelty.** Saturn: every `mu_pert` is 40-8,000x
   below the already-marginal JEG reference. Neptune: Triton geometrically disqualified (157°
   inclination, not a borderline case); the remaining small moons are 200-2,800x below the JEG
   reference. Miranda-Ariel: Δi=4.4° (near the 5° sanity-filter edge) plus a near-degenerate
   `mu_base` for the base pair itself. Io-Callisto: decent forcing but a 9.47 period ratio with no
   plausible low-integer seed orbit to build a torus from.

## 5. Considered and explicitly rejected (with reasons)

- **Any Saturn CCR4BP edition** (Mimas-Tethys, Enceladus-Dione, Tethys-Dione, Titan-Hyperion,
  Mimas-Enceladus) — quantified in §2: every `mu_pert` is 1.5-4 orders of magnitude below the
  already-flagged-as-barely-tractable JEG reference. This sharpens `#686`'s qualitative Saturn
  rejection into a numeric one; no new evidence changes the verdict.
- **Any Neptune CCR4BP edition** — Triton geometrically disqualified (157° inclination, not
  borderline); the small regular moons are mass-negligible (200-2,800x below reference). No
  tractable Neptunian lane exists at all, a stronger negative than Saturn's.
- **Uranus Miranda-Ariel** — Δi=4.4° sits close to the sanity filter's own 5° edge and the base
  pair's `mu_base=7.4e-7` is itself near-degenerate; not rejected as impossible, but clearly worse
  than Umbriel-Titania or Ariel-Umbriel on every axis.
- **Jupiter Ganymede-Callisto and Uranus Titania-Oberon** — both PUBLISHED (§3); rejected on
  novelty, not tractability (both are otherwise strong candidates by the numbers).

## 6. User decision points (flagged, not assumed)

1. **GO/NO-GO on dispatching a real per-pair build** for Io-Europa and/or Io-Ganymede (items 1-2)
   — each is its own multi-day build chain analogous to `#689`-`#691`, even though the EOM/STM/
   corrector/whisker MODULES themselves are almost certainly directly reusable (same primary,
   same code structure, just a different `mu`/`mu_gan`/`a_gan`/`omega_gan` parameter set via a new
   `..._default()` constructor function) — this screening pass did not verify that reuse claim
   with code, only with the model's documented structure.
2. **Whether to acquire, digest, and CORPUS_INDEX the newly-found Aryan & Fitzgerald AAS 24-103**
   paper — not previously in this project's corpus; found only by live search during this task.
   Acquisition/digest is out of this screening task's own scope but is a natural, cheap follow-on
   (the paper directly informs any future Jovian CCR4BP work, including whether its PCCFBP
   Europa-Ganymede results agree with `#690`'s own independently-derived torus).
3. **Whether the Kumar arXiv:2509.03655 anchor's `inherited-unverified` provenance should be
   ground-truthed** (acquire + read the actual Section 6.2 + the companion secondary-resonance
   extension) before it is trusted as a hard disqualifier for any FUTURE Uranus-Titania-Oberon
   catalogue-relevant claim — sufficient to disqualify a build target here, per `can_anchor_
  decision()`'s own discipline it is NOT yet decision-grade for anything stronger.
4. **A dedicated, deeper lit-check specifically for Europa-Callisto** (item 4) before it could be
   promoted above Umbriel-Titania — this pass only lightly checked it.

## 7. Recommended dispatch order

**Io-Europa first** (item 1: cleanest geometry of anything surveyed, including the already-built
reference system, and clears the novelty check). **Io-Ganymede second**, effectively free to
de-risk in the same pass since its perturber-forcing regime is byte-for-byte the same strength
already validated by `#690`. Uranus Umbriel-Titania (item 3) is a legitimate third candidate but
should wait on decision point 4's own confirmation and, honestly, a second more targeted lit pass
given the Uranian neighborhood's demonstrated literature activity — not because the numbers are
weak (they are the best of any non-Jovian pair) but because the scoop-context risk is measurably
higher there than in the Jovian system. No Saturn or Neptune CCR4BP build is recommended at any
priority; the mass-forcing negative there is now quantified, not just asserted.
