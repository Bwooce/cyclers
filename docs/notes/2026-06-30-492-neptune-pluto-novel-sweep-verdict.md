# #492 — novel moon-cycler sweep (Neptune + Pluto): no novel candidate

**Date:** 2026-06-30. First pass of the novel moon-cycler discovery frontier (the unswept
viable systems — Uranus was done in #285/#312). Ran the #285 prioritized repeated-moon scan on
the two systems with a massive flyby body (the #489 constraint) + no obvious published cycler
lit: **Neptune (Triton-flyby)** and **Pluto (Charon-flyby)**. **Verdict: no novel candidate** —
Neptune is a clean negative; Pluto's closures are infeasible, model-invalid, or published-adjacent.
`scripts/scan_492_neptune_pluto.py` → `data/scan_492_{neptune,pluto}.jsonl`. No catalogue change.

## Neptune (Triton, Proteus) — CLEAN NEGATIVE
70 enumerated, 6 closed, **0 SILVER, 0 near-miss (< 1 km/s even)**. Triton-Proteus-Triton forms
no sub-gate cycler: the SMA ratio is large (Triton 354,800 km / Proteus 117,600 km ≈ 3.0) and
Triton is retrograde — the ideal-model conic + resonance phasing does not close near the gate.
Conditional on the available SATELLITES moons (Triton, Proteus only; Nereid/inner moons absent).

## Pluto (Charon, Nix, Hydra) — closures, but NONE novel
210 enumerated, 22 closed, **0 SILVER, 9 BRONZE** sub-gate closures (Charon-Nix-Charon 0.017–0.031,
Hydra-Nix-Hydra 0.009–0.026, Nix-Charon-Nix, Nix-Hydra-Nix, …). All disqualified, three ways:
1. **Small-moon-flyby ones (Hydra-Nix-Hydra, Nix-Hydra-Nix): physically INFEASIBLE** (#489) — Nix
   μ=0.0015, Hydra μ=0.002 km³/s² (tinier than the Saturn icy moons) cannot provide a gravity assist.
2. **Charon-flyby ones (Charon-Nix-Charon): model-INVALID** — Charon is ~12% of Pluto's mass; the
   Pluto-Charon system is a true BINARY, so the restricted-3-body **point-primary assumption the
   scan uses is violated**. A "cycler" computed in a Pluto-point-primary model is not physically
   valid (the real dynamics are circumbinary; the small moons orbit the BARYCENTER in 3:4:5:6
   resonance with Charon). This needs a circumbinary (P-type) model, not the CR3BP genome.
3. **Lit-flagged:** the scan's `literature_check` matched **Persephone** (Howard, Stern et al.,
   "Persephone: A Pluto-system Orbiter and Kuiper Belt Explorer," PSJ 2(2):56, 2021,
   arXiv:2102.08282) — verified live: Persephone designs a Pluto-system **minor-satellite orbiter
   tour** (closest-approach encounters, <300 m/s). Related architecture; the Pluto-system tour
   regime is published, not virgin ground.

## Standing — the novel moon-cycler frontier is exhausted
Across the whole arc, the **only genuinely-novel moon-cycler hit remains the #312 Uranus
Umbriel-Oberon quasi-cycler** (V4-catalogued). The viable-flyby unswept systems add nothing:
Neptune (Triton-Proteus) doesn't close; Pluto is binary-model-invalid + published-adjacent. The
#489 massive-flyby constraint + the binary-regime issue together explain why: the systems with a
massive flyby body AND a valid point-primary model AND no published lit are essentially just
Uranus, which is already mined.

**Re-open keys** (method-versioned empty, for `data/empty_regions.jsonl`): (a) a circumbinary
(P-type / Pluto-Charon) cycler genome would make the Pluto Charon-flyby closures physically
meaningful — a genuinely different model, not a re-run; (b) more Neptunian moons (Nereid + the
inner Naiad…Proteus group) if added to SATELLITES could reopen Neptune. (c) Uranus-deeper (3D /
all pairs) is the remaining in-Uranus angle but #285/#312 already swept it hard.

## Net
#492 first pass complete: no new novel candidate, no new row. The frontier's lesson: novel
moon-cycler discovery is gated by (massive flyby body) ∧ (valid point-primary model) ∧ (no
published lit) — a near-empty intersection now that Uranus is mined.
