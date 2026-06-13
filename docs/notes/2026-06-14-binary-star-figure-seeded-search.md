# Binary-star figure-seeded cycler search (#255) — provisional negative

**Date:** 2026-06-14
**Task:** #255 (figure-read ICs, not continuation — the complement to #252).
**Result:** **provisional NEGATIVE** for the figure's stable prograde binary-star
cyclers under figure-read seeding. No catalogue writeback.

## What was built (validated, kept)
- `src/cyclerfinder/search/binary_star_search.py` — figure-seeded fixed-mu search
  with a **winding-number (k1,k2) topology classifier** (k1 = signed revolutions
  about P1 at (-mu,0), k2 = about P2 at (1-mu,0)). The classifier **reproduces the
  published Earth-Moon (3,1) and (1,1) Ross labels exactly and prograde** — the
  reproduce-before-search gate, locked in as a permanent test
  (`tests/search/test_binary_star_search.py`, 3 passing). Mandatory Barden
  stability + independent-Radau Jacobi cross-check on topology survivors.
- `scripts/binary_star_figure_search.py` — driver over the Fig. 3 (x0, C) windows
  at mu = 0.1 / 0.3 / 0.5. CAVEAT: the full enumeration (n_x0 x n_c x 6 hc x 2
  signs, each converged orbit getting a 4000-step winding propagation) is too slow
  for the mu=0.1 exterior box (long-period orbits) — use lean grids / the focused
  classification path below.

## The deciding result (mu = 0.3, the probe's apparent lead)
A coarse probe flagged a STABLE orbit at mu=0.3 (x0=0.5987, C=3.0, hc=3,
nu=-0.498) as a "cycler" via the crude `reaches-both-realms` proxy, which *looked*
like it contradicted #252. The rigorous winding classifier corrects this:

```
x0=0.59870277  C=3.00000000  T=11.553017  cross_res=5.8e-14  radau_dJ=1.0e-13
(k1,k2)=(1,0)  w1=-1.000  w2=0.000  prograde=False  xrange=[-1.932, 0.985]
nu=-0.49809  stable=True
FIGURE-MATCH (3,1)? NO -- different family (1,0), retrograde
```

The orbit is a genuine, stable, machine-precision periodic orbit, but it is a
**(1,0) RETROGRADE** libration about the primary — it crosses past P2 without
*winding* it (w2=0), and runs the wrong direction. It is **not** the depicted
prograde (3,1) cycler. The crude proxy was a false positive; the winding gate
(validated to reproduce the known EM labels) rejects it correctly.

## Honest verdict
- #255 at mu=0.3 **does not contradict #252** — it confirms it with a rigorous
  topology gate: the figure's stable prograde (3,1) was **not located** under the
  figure-read seeds tried (hc 1-4, x0 in [0.4, 0.75], C in [2.96, 3.14]).
- The figure prints **no numbers**, so the (x0, C) seed is a guess; absence under a
  bounded guessed-seed scan is a weak negative, not proof the cycler is absent.
- The discipline worked: the rigorous classifier caught a "discovery" the crude
  proxy would have mis-claimed. **A clean negative is a success.**
- mu=0.1 (exterior (1,3)) and mu=0.5 ((1,1) equal-mass) not yet swept (the mu=0.1
  box is the slow one). They remain available via the leaner focused path.

Pairs with [[feedback_orbit_closure_discipline]] and
[[feedback_literature_novelty_check_baseline]] (no novelty claim without the gates).
