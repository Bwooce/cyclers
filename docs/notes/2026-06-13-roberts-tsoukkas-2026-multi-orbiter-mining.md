# Roberts-Tsoukkas & Ross 2026 journal — FULL-TEXT transcription + reproduction attempt (#251)

**Date:** 2026-06-13
**Task:** #251 — full-text fetch + transcribe + reproduce the journal extension
of AAS 25-621. WebFetch was re-enabled for this run (blanket allow); the #250
sweep could only read the abstract.
**Source:** M. Roberts-Tsoukkas (advisor S. D. Ross), "Stable Prograde
Earth-Moon Multi-Orbiter Cyclers via Three-Body Dynamics," journal manuscript
hosted by the Virginia Space Grant Consortium (posted Apr 2026, 7 pp.). It is
the journal/undergraduate-research extension of **AAS 25-621** (Ross &
Roberts-Tsoukkas 2025, already mined and ADOPTED as the 5 `ross-rt-em-cycler-*`
rows; see `docs/notes/2026-06-11-ross-roberts-tsoukkas-2025-mining.md`).
**Writeback: NONE.** Sourced-golden discipline applies.

---

## Headline (corrects the #250 abstract-only verdict)

**WebFetch WORKED** (fetched the 1.1 MB PDF; the page→markdown step couldn't
read its compressed text streams, but the binary was saved and `pdftotext`
extracted the full 7-page text layer plus the figures were read directly).

**The full text DISCONFIRMS the #250 head-start assumption that this paper
contains transcribable higher-μ (μ, C, T) tuples and/or printed ICs.** The
journal version is a short (7-page) manuscript that **prints NO numeric data of
any kind** — no tables, no Jacobi constants, no periods, no stability indices,
no initial conditions. The binary-star-range results exist ONLY as three
rotating-frame orbit plots (Figure 3) and one schematic bifurcation diagram
(Figure 4), whose axes carry no usable goldens (Figure 4 is even plotted in
*offsets from unprinted critical values*: `C − C_cr` and `x − x_cr`,
×10⁻³–10⁻⁴ scale).

**Reproduction outcome: a clean NEGATIVE.** There is nothing to reproduce as a
sourced golden — the EXPECTED side does not exist in print. No new rows are
proposable from this paper. This is reported honestly per the
orbit-closure-discipline "a clean negative is success" rule.

This note SUPERSEDES the verdict in
`docs/notes/2026-06-13-roberts-tsoukkas-2026-multi-orbiter-journal-mining.md`
(written under the WebFetch denial), which assumed the PDF held tables to fetch.
It does not. The "USER FETCH the PDF" follow-on in that note is now DISCHARGED:
the PDF was fetched and contains no transcribable tuples.

## Bibliographic diff vs AAS 25-621 (#209)

| | AAS 25-621 (2025) | This journal version (2026) |
|---|---|---|
| Authorship | Ross & Roberts-Tsoukkas | Roberts-Tsoukkas (advisor Ross) |
| Length | 16 pp. | 7 pp. |
| Numeric tables | Tables 1–4 (15–16-digit C/T/x; 6 tangency consts) | **NONE** |
| Printed ICs | none (recoverable via 1-D solve) | **none** (no C/T either) |
| Earth-Moon families | 5: (1,1),(2,1),(3,1),(3,2),(3,3) | mentions (1,1) and (3,3) (Fig. 2); (3,2) bifurcation (Fig. 4) — all figure-only |
| Higher-μ families | none | **3 figures** at μ=0.1, 0.3, 0.5 (no numbers) |
| Method | full corrector + Barden monodromy + arc-length cont. (Eqs. 9–16) | prose summary only; equations reduced to EOM (1), U (2), C (3) |

The journal version is **strictly a numeric subset (in fact empty)** of the AAS
manuscript on the catalogue-grade columns. Everything quantitative we hold
already comes from AAS 25-621 and remains the version-of-record for our 5 rows.
The AAS Tables-3/4 C_(2,1)/C_(3,1) self-inconsistencies our two `kind: conflict`
data_gaps parked "for the 2026 journal" are **NOT resolved** here — the journal
prints no C-bound columns at all, so the conflicts stand unaddressed.

## What is genuinely NEW (qualitative only)

1. **Binary-star-range cycler families exist** (Fig. 3), at mass parameters the
   catalogue has never held:
   - **μ = 0.1**: a class **(1,3)** cycler, described as an *exterior* cycler
     (second perpendicular crossing in the exterior realm). Stable.
   - **μ = 0.3**: a class **(3,1)** cycler. Stable.
   - **μ = 0.5** (equal masses): a class **(1,1)** cycler. Stable.
     (NB: the Fig. 3 panel is labelled μ = 0.5 and the caption says "equal
     masses" — consistent; the CONCLUSION prose garbles this as "μ = 0.1" in
     two places. Treat μ = 0.5 from the figure panel as authoritative; the
     prose "(μ = 0.1)" for the equal-mass case is a typo.)
   - Binary-star range defined in-paper as **0.1 ≤ μ ≤ 0.5** (Results §).
2. **"Stable subfamily for every family" elevated to a conjecture with a
   mechanism.** The AAS version observed it empirically; the journal argues it
   from the **saddle-center bifurcation** at each family's C^max: that
   bifurcation always spawns one stable + one unstable branch, so a stable
   subfamily exists near C^max for every (k1,k2) family. Stated as "believed
   highly probable," not proven. Backed by the Global Orbit Theorem (Koon et
   al.) guaranteeing unstable cyclers at small C, forcing a creation point C* in
   between, observed to always be saddle-center type.
3. **A second structural conjecture:** at a given μ and any C, there exist at
   minimum zero and at maximum **two** symmetric cyclers (from the geometric
   method: the image P_B^{(k2+1)/2}(R) intersects {ẋ=0} in at most two points).
4. The (3,3) Earth-Moon family is now stated to have **5 stable subfamilies**
   ("the cycler family with the most stable subfamilies considered thus far") —
   consistent with the AAS note's "five distinct stability windows" for (3,3).

None of (1)–(4) carries a transcribable number.

## Reproduction — what is and isn't possible

**Not possible from this paper.** A same-model reproduction requires an EXPECTED
(μ, C, T) (and ideally σ) to target. The journal prints none. We could, in
principle, *attempt to rediscover* a (1,3) cycler at μ=0.1, a (3,1) at μ=0.3, or
a (1,1) at μ=0.5 with our adopted Ross corrector and arc-length continuation —
the method is exactly what we already hold and validated (ADOPTED 2026-06-12).
**But that would be a DISCOVERY, not a reproduction**: any C/T/σ we produced
would be our own computed output with no published value to check it against —
the EXPECTED side would be missing, violating the sourced-golden rule. Such a
result could not be a sourced V1 row; at best it would be a self-found candidate
needing the full V-ladder gauntlet, with the figure used only as a topology
sanity check (does the rediscovered orbit's rotating-frame shape match the
plotted petal/loop count).

Decision under ambiguity: **do NOT run a speculative binary-star rediscovery in
this pass.** Rationale: (a) it is out of the #251 reproduction scope (which is
"recover the paper's printed members," and there are none); (b) per the
discovery program a self-found candidate at a never-catalogued μ would need the
full gauntlet, not a quick recover-and-check; (c) it belongs to Track-C
discovery-daemon work, not a transcription task. Logged as a follow-on below.

## Proposed catalogue rows

**NONE.** No row is proposable: no sourced (μ, C, T) exists. Would-be
validation level: **N/A** — a sourced row needs a sourced tuple, and there is
no tuple. (Were we to self-discover a binary-star member, it would enter as a
self-found candidate subject to the full ladder, not as a sourced V1.)

The 5 existing `ross-rt-em-cycler-*-2025` rows (now V2 via
`docs/notes/2026-06-13-ross-v2-longspan-evidence.md`) are **unaffected** — their
goldens come from AAS 25-621, which this paper does not revise.

## Follow-on (review-gated, NOT done here)

1. **Track-C discovery candidate (not sourced):** point the adopted Ross
   corrector + arc-length continuation at the three binary-star topologies
   ((1,3)@μ=0.1 exterior, (3,1)@μ=0.3, (1,1)@μ=0.5) to self-discover stable
   members in a never-catalogued μ regime. Output would be candidate rows
   requiring the full V-ladder; the figures serve only as topology checks. This
   is genuinely new catalogue territory but is a DISCOVERY effort, separate from
   this transcription task.
2. **Prioritizer prior (Track B):** adopt the saddle-center / "stable subfamily
   near C^max for every (k1,k2) family at any μ" claim as a non-empty-region
   prior for the discovery daemon — a qualitative prior (no number needed),
   already partially captured by the #250 journal note. Predicts every (k1,k2)
   topology has a stable, ballistic member near its C^max regardless of μ.
3. **Negative-registry entry:** record "Roberts-Tsoukkas 2026 journal — no
   numeric data; binary-star families figure-only" so a future agent does not
   re-fetch expecting tables. (The PDF *was* fetched; there is nothing to
   transcribe.)

## Dedup / supersession

- Supersedes `2026-06-13-roberts-tsoukkas-2026-multi-orbiter-journal-mining.md`'s
  "ADOPT gated on user fetch" verdict: fetched, no data, → no adoption possible
  from this paper; the qualitative priors are the only takeaway.
- Does NOT supersede AAS 25-621 as our citable version-of-record for the 5
  Earth-Moon rows — the journal prints none of those numbers, so AAS 25-621
  remains the numeric source.
- Forward-citation lineage unchanged (same lab; Braik & Ross already mined).
</content>
</invoke>
