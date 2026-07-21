# Topputo & Belbruno 2015 — Earth-Mars Transfers With Ballistic Capture

**Deep-read digest, 2026-07-22 AET.** Acquired from the user directly (PDF attachment) for task
`#681` (Sun-Mars WSB repeating-capture quasi-cycler, `#679` shortlist item 2). Full text-layer
PDF, no OCR needed.

## Header

- **Title:** *Earth-Mars Transfers With Ballistic Capture*
- **Authors:** F. Topputo (Politecnico di Milano), E. Belbruno (Princeton)
- **Venue:** Celestial Mechanics and Dynamical Astronomy 121, 329-346 (2015)
- **DOI:** 10.1007/s10569-015-9605-8
- **Format:** Full peer-reviewed article, 18 pages, text-layer PDF (readable directly, no OCR).

## What the paper actually is

Constructs a new class of Earth-to-Mars transfers ending in **ballistic capture** — the
spacecraft first transfers from Earth to a point `x_c` several million to tens-of-millions of km
from Mars (still on a heliocentric orbit close to Mars's own orbit), then follows the natural
dynamics of the **planar elliptic restricted three-body problem** (Sun-Mars, `e_p = 0.093419`)
from `x_c` into ballistic capture at a chosen Mars periapsis radius `r_p`, with no further
maneuver needed after the single `ΔV_c` burn at `x_c`. This is the interplanetary-distance
analog of the same "stable sets" methodology Belbruno/Topputo used for lunar ballistic-capture
transfers (Hiten, GRAIL).

**Key machinery** (Sect. 4, "Mars stable sets and ballistic capture orbits"):
- A grid of initial conditions on a radial segment `l(θ)` from Mars, at periapsis of an osculating
  ellipse with fixed eccentricity `e` and initial true anomaly `f0` of the primaries.
- **n-stability**: an orbit is *n-stable* if it completes `n` full revolutions about Mars without
  ever completing one revolution about the Sun, then returns to negative Kepler energy w.r.t.
  Mars. `n`-*unstable* is the complement. Backward-stability (`-m`-stability) is defined
  analogously by backward integration.
- **Weak stability boundary** `∂W_n`: the locus where stability changes — a Cantor-like structure
  (alternating stable/unstable intervals) on each radial line, giving the characteristic "gaps" in
  the capture-cost plots (Fig. 8b, explicitly noted in the paper's own discussion).
- **Capture set** `C^n_{-1}(e, f0) = W̄_{-1}(e,f0) ∩ W_n(e,f0)` — points that are `-1`-unstable
  (escape Mars in backward time / approach in forward time) AND `n`-stable (perform ≥n forward
  revolutions about Mars once captured). These are the orbits of practical interest: candidate
  ballistic-capture states.
- Grid parameters used in the paper (their Sect. 4): `r` (radial distance, Δh=50km for h≤30,500km
  then Δh=500km up to h=250,000km), `θ` (0-360°, Δθ=1°), `e∈{0.90,0.99}`, `f0∈{0,π/4,π/2}`,
  `n∈{-1..6}` — **375,394 initial conditions** integrated per grid, via variable-order
  Adams-Bashforth-Moulton with Levi-Civita regularization near Mars close approach.

## Concrete, digit-grade reproduction targets (positive-control candidates for `#681`)

These are the specific numbers a from-scratch reimplementation should be checked against before
trusting any new Sun-Mars WSB search:

**Physical constants** (Table 4 — note the paper's own typo: both Mars rows are mislabeled `a_E`/
`e_E` in the printed table, values are for Mars):
- `μ_S = 1.32712×10^11 km³/s²`, `AU = 149,597,870.66 km`
- `μ_E = 3.98600×10^5 km³/s²`, Earth: `a=1.000000230 AU`, `e=0.016751040`
- `μ_M = 4.28280×10^4 km³/s²` (Mars), `a=1.523688399 AU`, `e=0.093418671`
- Sun-Mars CR3BP-elliptic mass parameter `μ = m_M/(m_S+m_M) = 3.2262081094×10^-7`

**Reference bitangential Hohmann transfers** (Table 5, four perihelion/aphelion combinations —
the paper's own baseline patched-conics sanity check):

| Case | Earth @ | Mars @ | ΔV1 (km/s) | ΔV2,∞ (km/s) | ΔV total (km/s) | Δt (days) |
|---|---|---|---|---|---|---|
| H1 | perihelion | perihelion | 2.179 | 3.388 | 5.568 | 234 |
| H2 | perihelion | aphelion | 3.398 | 2.090 | 5.488 | 278 |
| H3 | aphelion | perihelion | 2.414 | 3.163 | 5.577 | 239 |
| H4 | aphelion | aphelion | 3.629 | 1.881 | 5.510 | 283 |

**Ballistic-capture cost vs. periapsis radius** (Table 3, for `e=0.99` capture states,
`f0` in the first quadrant — the paper's headline result, Fig. 8b's underlying data):

| Point | r_p (km) | ΔV_c (km/s) | ΔV2 [H3-equiv Hohmann] (km/s) | Savings S (%) | Δt (x_c→r_p, days) |
|---|---|---|---|---|---|
| (A) | 49,896 | 2.033 | 2.116 | −4.0 | 434 |
| (B) | 73,896 | 2.036 | 2.267 | −11.3 | 433 |
| (C) | 91,897 | 2.039 | 2.344 | −14.9 | 432 |
| (D) | 113,897 | 2.041 | 2.414 | −18.2 | 431 |

`ΔV_c` is remarkably **flat (~2.03-2.04 km/s) across a wide range of `r_p`** — a sharp contrast
to the Hohmann `ΔV2` cost, which grows with `r_p`. This flatness is itself a distinctive,
checkable signature of the method (result A, Eq. 7 in the paper: `ΔV_c < ΔV_2` for `r_p` above a
threshold).

**Crossover periapsis radii** (Table 2 — where ballistic capture starts winning over Hohmann):

| f0 | r_p^(1) (km) | r_p^(2) (km) | ΔV_c (km/s) |
|---|---|---|---|
| 0 | 29,000 | 46,000 | 2.09 |
| π/4 | 26,000 | 40,000 | 2.03 |
| π/2 | 22,000 | 34,000 | 1.96 |

**Named worked examples** (Sect. 5, Fig. 6/7 — the two fully-constructed transfer cases):
- **Case 1**: `x_c` chosen from `C^6_{-1}(0.99, π/4)`, `x_c ≈ 1×10^6 km` from Mars (trailing
  Mars slightly). ~1 year transit `x_c → r_p`.
- **Case 2**: `x_c` chosen from `C^6_{-1}(0.99, π/2)`, `x_c ≈ 23×10^6 km` from Mars — notably far,
  yet approximately the same transit time as Case 1.
- The `N` (number of points) reported for the sample capture set `C^6_{-1}(0.99, π/4)` shown in
  Fig. 3 is **N=597** — a concrete count to check a from-scratch grid computation against
  (same `e`, `f0`, `n` — exact grid resolution matters for an exact match, but the *order of
  magnitude and qualitative Cantor-gap structure* should reproduce regardless).

## What this paper does NOT give us (relevant to `#681`'s actual object class)

This paper is entirely about **one-shot** Earth→Mars ballistic-capture transfers — a single
capture event, explicitly not a repeating capture↔escape cycle. Sect. 4.2 ("Long term behavior of
the capture orbits") integrates one capture orbit backward 50 Mars revolutions (~94 years) and
finds it does **not** return to a second ballistic capture within that span — i.e. the paper's
own worked example is evidence AGAINST easy repeatability, not for it. This is directly relevant
to `#681`'s own honestly-stated risk (Belbruno 2004 Thm 3.58's chaos argument): the one existing
long-integration data point in the literature shows a *typical* capture orbit does NOT recur
ballistically over decadal timescales. `#681` should not read this paper as making repeating
capture look easy — if anything it's a mild negative signal the task's own risk framing should
cite explicitly.

## Catalogue / KNOWN_CORPUS relevance

Not itself a `quasi_cycler`/cycler admission candidate (one-shot transfer, not a repeating
object) — this is a **methodology and positive-control source**, not a family to mine for rows.
No KNOWN_CORPUS anchor is proposed; register as a methods reference only.

## Action items for `#681`

1. Reproduce the stable-set/capture-set machinery (`W_n(e,f0)`, `C^n_{-1}(e,f0)`) for the
   Sun-Mars elliptic-restricted problem using the Table 4 constants above.
2. Positive control, in order of cheapness: (a) reproduce the four Hohmann baseline cases (Table
   5) as a basic patched-conics sanity check; (b) reproduce the flat `ΔV_c ≈ 2.03-2.04 km/s`
   trend across `r_p ∈ [49896, 113897]` km at `e=0.99` (Table 3) as the actual three-body-dynamics
   validation; (c) if pursuing exact grid reproduction, target `N=597` for `C^6_{-1}(0.99, π/4)`.
3. Before proposing a repeating-capture search, address the Sect. 4.2 finding above explicitly —
   it is direct evidence from this paper's own worked example that ballistic capture orbits do
   not obviously self-repeat, which should sharpen (not just gesture at) `#681`'s own stated risk.
