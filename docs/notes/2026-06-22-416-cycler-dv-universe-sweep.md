# #416 — Corpus-wide cycler maintenance-ΔV sweep (close-out)

**2026-06-22.** To close #416 ("multi-source the #415 ΔV bands") with no cycler paper
left unchecked, a 5-agent parallel sweep read every cycler-relevant PDF in the corpus
(~45 papers across Earth-Mars, Earth-Moon, VEM/CGE triple, human-Mars) for deterministic
**maintenance/station-keeping/DSM** ΔV with its **basis**. Plus the user-supplied
Friedlander 1986 (the canonical early primary, previously the lone genuine acquisition
gap — see `2026-06-22-digest-friedlander-niehoff-byrnes-longuski-1986.md`).

The question under test: is the near-ballistic tier set **< 1 / < 10 / < 300 m/s total
over 7 cycles (real-ephemeris, best launch window)** corroborated by an *independent*
(non-Russell-lineage) same-basis source?

## Consolidated evidence (deterministic cycler maintenance ΔV only; taxi/transfer/injection excluded)

| Source | System | Value | Basis | Class | Model | vs tiers |
|---|---|---|---|---|---|---|
| **Russell-Ocampo 2006 JGCD** | E-M | **<1 / <10 / <300 m/s** (9/39/74 parents) | total/7-cycle, best window | ballistic→near-ballistic | real-eph | **defines the tiers** |
| Russell 2004 diss. | E-M | same 9/39/74 | total/7-cycle | same | real-eph | = 2006 (same author) |
| Friedlander-Niehoff-Byrnes-Longuski 1986 | E-M | VISIT "free orbit, navigation-only" (≈0); Escalator 230 m/s/orbit → ~1.6 km/s/15yr | qual. + per-orbit | ballistic (VISIT) + powered (Escalator) | 20-yr prop. | **corroborates ballistic floor (1986, pre-Russell); no tier census** |
| McConaghy 2006 (S1L1) | E-M | ~10 m/s (outbound, best); 0 (33-yr 4-veh) | total/30-yr (=7 two-syn) | near-ballistic | DE405 | **corroborates <10 tier (same basis; Longuski school)** |
| McConaghy 2004 | E-M | 0 (n=7, VISIT-1, ballistic S1L1); Aldrin 1.73/2.04 km/s | total/repeat | ballistic + powered | DE405 | corroborates 0-floor; powered far above |
| Byrnes-Longuski-Aldrin 1993 | E-M | 1.73 / 2.04 km/s, 3 of 7 orbits | total/15-yr | powered | real-eph | powered band |
| Rauwolf 2002 | E-M | 1.561 / 1.605 km/s, 3 of 7 | total/15-yr | powered | conic | powered band |
| Chen-McConaghy-Landau-Longuski-Aldrin 2005 | E-M | 0.58–1.05 km/s/synodic | per-synodic | powered (low-thrust) | DE405 | powered band (≫ tiers) |
| Chen 2002 | E-M | 1.85 km/s/synodic (3-veh); assumed 300 m/s TCM | per-synodic | powered | conic+real dates | powered band |
| Landau-Longuski 2006 | E-M | 0.78 km/s (cycler DSM, taxi-co-optimized) | avg/synodic | powered (chosen) | conic | powered; affirms DSM=0 ballistic versions exist |
| Howe 2025 (cites Rauwolf) | E-M | 1500 m/s/15yr (SEP escalator) | total/15-yr | powered SEP | idealized | powered (escalator) |
| **Sanchez-Net 2022 (Pony Express)** | E-M | ≤60 m/s total; ≤10 m/s/flyby (**search filter**) | total/3-cycle | near-ballistic | **conic only** | **NOT corroboration — defers real-eph as "large"** |
| Pascarella 2022 (Pony Express) | E-M | ~163 m/s (from 2 kg @ Isp 4155 s) | total/8-flyby/~6yr | low-thrust SEP | real-eph | magnitude-consistent w/ <300; not impulsive, 1 candidate |
| **Jones-Hernandez-Jesick 2017** | **VEM** | ≤200 m/s "differentially correctable to entirely ballistic"; converges ballistic | per-flyby → whole-traj | ballistic | real-eph | **independent (JPL/Jesick): corroborates floor + ~200-300 boundary in another system** |
| **Liang 2024** | **CGE** (Jovian) | max **1.04e-7 m/s/cycle**; 200 m/s epoch-sensitivity | per-cycle (10-cycle) | ballistic | hi-prec eph | **independent: corroborates <1 floor in another system** |
| Patel 2019 | S1L1 | 0 (ballistic, inherits McConaghy) | per-cycle | ballistic | real-eph (cited) | supports <1 |
| Pontani-Conway 2018 | E-M | "small or ballistic, no propellant" (no #) | — | ballistic | — | supports (qual.) |
| Genova-Aldrin 2015 | **E-Moon** | 20–62 m/s/cycle (~39 m/s/month); 4-petal 19 m/s/2mo | per-lunar-cycle | station-kept | real-eph N-body | EM-system band (not E-M) |
| Wittal 2022 | E-Moon | <50 m/s lifetime / 5 cycles (5-petal) | total/5-cycle | quasi-ballistic | perturbed N-body | EM-system band |
| Ross / Roberts-Tsoukkas 2025/2026 | E-Moon | 0 by construction (stable PCR3BP); **no m/s figure** | — | ballistic (stable) | PCR3BP | idealized; perturbed SK = future work |
| Merrill 2025 | E-Moon | 10–40 m/s/orbit (low-thrust **forcing** cost) | per-orbit | forced low-thrust | CR3BP | not natural-cycler upkeep |
| Rogers 2012 / 2015 | E-M | — (V∞-leveraging = **establishment** cost, not maintenance) | — | — | conic+STOUR | excluded (establishment) |

## Verdict (closes #416)

1. **The ballistic floor and the ~200–300 m/s "correctable-to-ballistic" boundary are
   now genuinely multi-sourced and multi-system** — independent of Russell's line:
   Friedlander 1986 (VISIT, pre-Russell), **Jones 2017 (VEM, JPL/Jesick, real-eph)**,
   **Liang 2024 (CGE, real-eph, ~1e-7 m/s/cycle)**, McConaghy 2004/2006 (S1L1). The
   recurring **200 m/s "differentially correctable to entirely ballistic"** threshold
   (Jones 2017; echoed by Liang) is an independent sourced anchor for the < 300 m/s
   tier boundary.
2. **The *exact* < 1 / < 10 / < 300 m/s tiered census (9/39/74 parents) remains
   Russell-Ocampo's own quantification** (peer-reviewed JGCD 29(2), 2006 = Russell 2004
   diss.). No other paper reproduces that specific parent-count census in the same
   basis. But it is no longer a *fragile* single source: the band STRUCTURE it encodes
   is corroborated across independent groups and three systems.
3. **Pony Express is NOT a corroboration** of the cutoffs: Sanchez-Net 2022 is
   patched-conic and explicitly defers real-ephemeris ΔV as likely "large";
   Pascarella 2022 is low-thrust fuel-optimal (one candidate, ~163 m/s), magnitude-
   consistent only.
4. **No outstanding acquisition gap remains** — Friedlander 1986 (the last genuine
   cycler-ΔV primary not in the corpus) is now acquired + digested.

**Net:** #416 closed. The bands are well-grounded; cite the near-ballistic *cutoffs* as
"Russell-Ocampo 2006 JGCD census; ballistic floor + ~200–300 m/s boundary independently
corroborated (Friedlander 1986; Jones 2017 VEM; Liang 2024 CGE; McConaghy 2004/2006)."

## Housekeeping found during the sweep
- `genova-aldrin-2015-earth-moon-cycler-AAS-15.pdf` == `...-free-return-...pdf` (byte-identical, md5 3efc6cb3) — same paper "A Free-Return Earth-Moon Cycler Orbit"; dedup'd.
- `chen-russell-ocampo-2002-...-free-return.pdf` is actually Chen, McConaghy, Okutsu & Longuski, "A Low-Thrust Version of the Aldrin Cycler" (AIAA 2002-4421) — not Russell-Ocampo, not free-return; renamed.
- `chen-russell-ocampo-2005-...-multiple-impulse-free-return.pdf` is actually Chen, McConaghy, Landau, Longuski & Aldrin, "Powered Earth-Mars Cycler with Three-Synodic-Period Repeat Time" (JSR 42(5)) — renamed.
