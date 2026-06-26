# #473 — Resonance-lock moon-tour search under a RELATIVE drift criterion

**Date:** 2026-06-26
**Status:** complete
**Predecessor:** #468 (in-band tours) → #470 (0/10 strict admit) → #473 STEP-0 (drift-shape re-adjudication, since corrected)

## The criterion (and why the absolute-km version was wrong)

A moon-tour is a bounded **quasi_cycler** (the #339 class) only if BOTH hold:

1. **Shape** — the per-cycle rendezvous-drift series oscillates-and-returns over
   ≥10 cycles (bounded, not monotonic divergence). [`scripts/_drift_shape_473.py`]
2. **Relative magnitude** — `max_over_cycles(drift_km) / SMA_of_normalising_moon ≤ R_REF`.

The prior STEP-0 used an **absolute** ~530k–800k km envelope, calibrated to #339's
530k km peak. That is **wrong**: 530k km is *Uranian-scale*. #339's turn-around moon
Oberon orbits at SMA 583,511 km, so 530k km is only **0.91 of that orbit radius**.
Carried verbatim into a larger system the absolute bar mis-scales — Jovian moons
(Ganymede SMA 1.07M km, Callisto 1.88M km) and Saturnian (Titan 1.22M, Rhea 527k)
legitimately drift *more in absolute km* purely because the system is bigger. The
absolute bar therefore either wrongly rejects everything Jovian/Saturnian or (as the
prior run did) flips everything on shape alone and admits 2.1M–3.76M km excursions
that are not rendezvous-quality.

**Fix — scale-invariant ratio.** Normalise by the SMA of the tour's **outermost
(largest-SMA) encountered moon**, the same moon #339's reference envelope is quoted
against. `R_REF = 0.91` calibrated directly from #339 (530,000 / 583,511 = 0.9083),
with a documented margin `R_MARGIN = 1.0` (worst-cycle drift up to one full
normalising-orbit radius is still rendezvous-scale). This is **non-circular**: #339
is the published/accepted V4 quasi_cycler, so "at least as good as #339 in *relative*
drift" is the principled bar. Helper: `scripts/_relative_drift_473.py`.

## Step 1 — Positive control (HARD GATE)

`scripts/_posctl_473_check.py` reproduces #339 (Umbriel-Oberon-Umbriel, 14.9406 d/leg,
n_rev=(1,1), rel_off 180°) at n_cycles=10: **PASS** — 10/10 cycles, strict V2 FALSE
(it is a quasi_cycler), drift `bounded-oscillating-and-returns` in the 86k–530k km
band. Its relative ratio is **530,000 / 583,511 (Oberon) = 0.908 ≤ R_MARGIN** by
construction.

## Step 2 — #470 re-adjudication under the RELATIVE criterion

Re-classifying the 10 stored #470 tours (`scripts/_v2_readjudicate_470.py`,
`out/readjudicate_470_summary.json`). Verdict = FLIP iff shape∈bounded AND
ratio≤1.0 AND not strict.

| Tour | norm moon (SMA km) | max drift (km) | shape | ratio | verdict |
|---|---|---|---|---|---|
| Galilean IEG  (Io-Eur-Gan-Io)      | Ganymede 1,070,400 |   832,718 | osc-and-returns | **0.778** | **FLIP → bounded quasi** |
| Galilean EGC  (Eur-Gan-Cal-Eur)    | Callisto 1,882,700 | 1,206,015 | osc-and-returns | **0.641** | **FLIP → bounded quasi** |
| Jovian IE pair (Io-Eur-Io)         | Europa     671,100 |   802,724 | oscillating     | 1.196 | reject (ratio>1) |
| Jovian GC pair (Gan-Cal-Gan)       | Callisto 1,882,700 | 2,140,770 | oscillating     | 1.137 | reject (ratio>1) |
| Galilean CGE  (Cal-Gan-Eur-Cal)    | Callisto 1,882,700 | 3,763,988 | osc-and-returns | 1.999 | reject (ratio>1) |
| Saturnian TRD (Tit-Rhe-Dio-Tit)    | Titan   1,221,870 | 2,443,520 | osc-and-returns | 2.000 | reject (ratio>1) |
| Saturnian RDT (Rhe-Dio-Tet-Rhe)    | Rhea      527,070 | 1,029,053 | osc-and-returns | 1.952 | reject (ratio>1) |
| Saturnian DTE (Dio-Tet-Enc-Dio)    | Dione     377,420 |   754,136 | osc-and-returns | 1.998 | reject (ratio>1) |
| Saturnian RD pair (Rhe-Dio-Rhe)    | Rhea      527,070 |   934,618 | osc-and-returns | 1.773 | reject (ratio>1) |
| Saturnian TE pair (Tet-Enc-Tet)    | Tethys    294,670 |   563,366 | monotonic-div   | 1.912 | reject (divergent) |

**Outcome: 2/10 flip** (Galilean IEG, EGC), not the prior 9/10. The prior
absolute-km classifier over-flipped 7 tours that are bounded in *shape* but drift
1.14–2.0 normalising-orbit-radii — far worse than #339 in relative terms — plus the
one genuinely divergent TE pair.

## Step 3 — Main resonance-lock search

`scripts/_resonance_lock_moontour_473.py`, 1716 cells (11 skeletons × 13 synodic
multiples × 12 relative-phase offsets × 1 anchor-longitude), each a full
`run_v2_moontour` at n_cycles=10 with the relative criterion. `phase0` was reduced
from 4 values to 1: the drift is `||r_final_k - r_final_0||`, invariant under a
whole-tour rotation (verified — identical max_drift to the decimal across phase0 ∈
{0,90,180,270}), so phase0 is a redundant axis. The #477 VILM ΔV-floor pre-filter
pruned 0/1716 (all #468-derived legs are within budget by construction). Run:
1716/1716 ok, 880 s, on the #474 multiprocessing prewarm substrate.

**Strict cyclers (max drift ≤ 50k floor — the prize): only the #339 POSCTL itself**,
rediscovered at mult=4.5 (max drift 19,602 km, ratio 0.034). This is apparatus
validation — the search re-finds the published strict cycler at its true lock — NOT
a new discovery. No NEW system produced a strict cycler.

**Bounded-quasi survivors (shape∈bounded AND ratio≤1.0): only 2 distinct NEW
families**, both outer-Galilean, plus the #339 sentinel:

| Skeleton | best ratio | best max drift (km) | mult | shape | #cells passing |
|---|---|---|---|---|---|
| Galilean EGC (Eur-Gan-Cal-Eur) | 0.693 | 1,303,998 | 3.0 | osc-and-returns | 132 |
| Galilean IEG (Io-Eur-Gan-Io)   | 0.766 |   820,230 | 0.5 | osc-and-returns | 132 |
| #339 POSCTL (Umb-Obe-Umb)      | 0.794 |   463,222 | 1.5 | osc-and-returns | 144 |

**8 of 10 #468 skeletons produce NO survivor at ANY resonance-lock multiple/phase**:
both inner-Jovian pairs (IE, GC), Galilean CGE, and all four Saturnian tours (TRD,
RDT, DTE, RD, TE) are either monotonic-divergence or ratio>1.0 across the entire
1716-cell grid. The per-skeleton "best ratio" line in the run log can show a *low*
ratio for these (e.g. Jovian GC pair 0.048) but that cell is shape=monotonic-
divergence — a reject; low relative drift on a divergent series is not a quasi_cycler.

The 2 surviving families are exactly the 2 that flip in Step 2 — fully consistent.
Resonance-locking the ToFs did not unlock any new family; it only sharpened the
envelope within the IEG/EGC/#339 basins.

## Step 4 — Survivors → official cross-check + literature

The search driver IS the official `run_v2_moontour` gauntlet, so the IEG and EGC
verdicts above are already official-gauntlet results (10/10 cycles, strict-V2 FALSE
→ quasi_cycler, bounded-oscillating-and-returns, relative ratio ≤ R_MARGIN).

Literature novelty: both survivors are **canonical Galilean multi-moon tours**
(Io-Europa-Ganymede and Europa-Ganymede-Callisto resonant linkages are the most
heavily-published moon-tour families — Strange/Campagnola/Russell Tisserand-graph
tours, Lynam Galilean tours, JUICE/Europa-Clipper tour design). They are **not
novel** and were already surfaced by #468/#470, not fresh from this search. No
catalogue write is warranted: they remain V0 quasi_cycler candidates whose envelope
(ratio 0.69-0.77, 0.8M-1.3M km) is real but materially worse than #339's (0.034 at
its true mult=4.5 lock); a catalogue row would need a sourced V4 reproduction, which
neither has. catalogue.yaml is unchanged.

## Verdict

Under a #339-calibrated **relative** drift criterion the #470 over-flip collapses
from 9/10 to **2/10** (Galilean IEG 0.766, EGC 0.693) and the resonance-lock search
adds **no new strict or quasi cyclers** — it only re-finds #339 as a strict cycler
(apparatus validation) and the two known Galilean quasi families; the other 8 #468
tours are genuinely off-basin, and no catalogue change results.
