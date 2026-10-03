# #577 — Opus adjudication of the #576 Galilean symmetric-closure survivor list

**Date:** 2026-07-12. **Task:** #577 (judgment-only — data all pre-existing from #576).
**Scope:** broaden the literature check beyond the R-S Table 3 spot-comparison #576 did; resolve
the "8 architecturally-mismatched Russell-Strange members" question; triage the 36 idealized-model
closures for real-ephemeris-pipeline readiness. **Out of scope (obeyed):** running any pipeline
stage; editing `data/catalogue.yaml`; claiming novelty absent a genuine clear + strong pipeline-
survival prior.

**Inputs read directly (not via the #576 summary):**
`data/enumerate_576_jupiter_galilean_symmetric_closures.jsonl` (36 gate-passing closures + 12
direction summaries), `data/compare_576_russell_strange_galilean.jsonl` (R-S Table 3 comparison),
`data/probe_576_galilean_repeat_check.jsonl` (repeat instrumentation), and the sourced digest
`docs/notes/2026-06-30-digest-russell-strange-2009-planetary-moon-cyclers.md`.

---

## Headline verdict

**0 of 36 closures clear the broadened literature check as genuine novelty. Full stop — push none
to the inclination-closure / eccentricity-kill-gate / real-ephemeris gauntlet.** Every one of the
six Galilean pairs is published double-cycler / triple-cycler / pump-tour territory. The single
pair that survives the *offline structural matcher* under the physically-correct topology label
(Io–Callisto) does so only because of a **corpus gap** (below), not novelty, and is independently
the least pipeline-worthy pair on physical grounds.

---

## 1. Broadened literature check (mandatory novelty-baseline gate)

I ran the deterministic core of `search/literature_check.py` — `_candidate_anchors()`, the
primary + body-set-subset + topology-intersection matcher that drives `check_literature`'s anchor
set and query trail — against a `CandidateSignature` for each of the six pairs
(`sequence=(anchor, flyby, anchor)`, real per-encounter V∞, real `n_rev`). Result depends on the
candidate's `topology_label`, so I ran both settings:

**(a) Body-set-only matching (`topology_label` empty — the historical default):**
every one of the six pairs overlaps 4–7 published anchors. Even Io–Callisto (the widest,
non-adjacent pair) overlaps 4: Strange/Campagnola/Russell, Niehoff (1970), Petropoulos–Longuski
(2000), Heaton–Strange–Longuski (2002). Under body-set matching **all 36 collide with published
Galilean moon-tour literature — nothing clears.**

**(b) Physically-correct `topology_label={"repeated-moon"}`** (a 2-body free-return A–B–A shuttle
that repeats the same encounter pair *is* the Aldrin / (k1,k2) repeated-encounter paradigm, per
the project's standard vocabulary):

| Pair (closures) | Overlapping repeated-moon anchors | Verdict |
|---|---|---|
| Europa–Ganymede (2) | Liang CGE; Hernandez IEG | collides |
| Io–Europa (2) | Hernandez IEG | collides |
| Io–Ganymede (6) | Hernandez IEG | collides |
| Europa–Callisto (6) | Liang CGE | collides |
| Ganymede–Callisto (14) | Liang CGE | collides |
| **Io–Callisto (6)** | **none** | structurally clears |

So 30/36 collide with the Liang CGE and/or Hernandez IEG triple-cycler anchors; only Io–Callisto
(3 geometries × 2 anchor directions = 6 closures) clears.

**Why the Io–Callisto "clear" is an artifact, not novelty — a genuine corpus finding.** The single
most direct prior for this *exact* construction — **Russell & Strange 2009, "Cycler Trajectories
in Planetary Moon Systems," JGCD 32(1), DOI 10.2514/1.36610** — is **absent from `KNOWN_CORPUS`.**
That paper's method is verbatim what #563/#576 reimplements: "the equations for calculation of
nonresonant free-return are derived from two-body dynamics, and a broad search in an ideal model
consisting of circular and coplanar celestial body orbits is performed to find *hundreds of ideal
model ballistic cycler geometries*" (confirmed by live WebSearch 2026-07-12 and by the on-disk
digest). It is a **repeated-moon** two-body free-return shuttle spanning the Galilean moons. Because
that anchor is missing, and because the Galilean pump-tour anchors that *are* present
(Strange/Campagnola/Russell, Petropoulos, Niehoff, Heaton) carry `topology_label={pump-tour,
mga-tour}` — which the repeated-moon filter correctly excludes as a different family — no anchor is
left to flag Io–Callisto. The clear is a false negative produced by a corpus gap, exactly the
false-negative-generator pattern the negative-results discipline warns about. Adding the R-S 2009
double-cycler anchor (scoped as #578 below) makes all six pairs collide correctly.

Live WebSearch independently corroborates all of this: R-S 2009 is the canonical Galilean
double-cycler prior (Ganymede–Io, Ganymede–Europa, Ganymede–Callisto, Europa–Ganymede searched);
a 2024 review ("Review of Trajectory Design and Optimization for Jovian System Exploration",
*Space: Science & Technology*, DOI 10.34133/space.0036) surveys the whole double-cycler field; and
Io–Callisto specifically is *discussed* in the literature as a dynamically awkward pair ("Callisto
does not participate in the Laplace resonance… presents unique challenges"). The published R-S
enumerative method trivially generates the Io–Callisto ideal-model geometry — it is a
known-method reproduction, not a discovery.

Per the module's own load-bearing discipline, "not-found is NECESSARY-NOT-SUFFICIENT for novelty."
Io–Callisto's structural clear does not meet the #577 bar ("genuinely clears **and** strong reason
to believe it survives the full pipeline"): see §3.

**Net: 0/36 novelty-clearable.** Standing is idealized quasi-cycler-class evidence, known-adjacent
(same standing as #312/#575), *not* a novelty candidate.

## 2. The "8 architecturally-mismatched Russell-Strange members" question — resolved

#576 declined to call its 37.6-day Ganymede–Callisto period matches reproductions because R-S's
"multi-loop" architecture supposedly differs from this project's 2-leg construction, and marked
only the two `legs=1` rows "architecturally comparable." **That framing is incorrect, and the
mismatch is not a genuine different-family situation.** The R-S 2009 digest (read page-by-page from
the acquired PDF) and the live search agree: R-S's architecture is *"free-return cyclers that
repeatedly shuttle a spacecraft between two bodies,"* built in the **same circular-coplanar ideal
model** by the **same enumerative free-return search** that finds *"hundreds of ideal model
ballistic cycler geometries."* This is the identical object #563 constructs. R-S's Table 1 lists
the pairs they actually enumerated — **Ganymede→Io, Ganymede→Europa, Ganymede→Callisto,
Europa→Ganymede** — so three of our six pairs (Io–Ganymede, Europa–Ganymede, Ganymede–Callisto)
sit squarely inside R-S's explicitly-enumerated territory. The R-S Table 3 "legs" column counts
free-return arcs / central-body revolutions per cycle (i.e. *different resonances of the same
two-body family*), not a distinct topology; "multi-loop" is a within-family parameter, not a
different architecture.

Consequently the Ganymede–Callisto period-37.6d matches are **members of R-S's explicitly-enumerated
Ganymede–Callisto double-cycler family** (a *known-class-member* status, the same category as
PC(3,2)), not coincidences and not out-of-scope. The V∞ signatures do not match R-S's *specific
tabulated representatives* (3.24/3.34, 3.18/3.26 km/s) — but Table 3 shows only a handful of
"promising" samples of hundreds of enumerated geometries, so a V∞ mismatch against the tabulated
few is expected and does **not** rescue novelty; our closures are simply other members of the same
published family. The topology-mismatch caveat was over-cautious in the "not a reproduction"
direction; the correct reading is *known-family-member*, which only strengthens #576's bottom line
(not novel). Nothing here needed a new 2-leg-vs-multi-loop comparison to be built — the shared
ideal model and shared enumerative method settle it.

## 3. Pipeline-readiness triage — recommendation: FULL STOP

Robustness ranking (the axes prior Uranian/Saturn work used — required-bend margin, V∞ magnitude,
tof), for completeness, is dominated by the low-V∞ near-Hohmann Ganymede–Callisto and
Europa–Ganymede geometries (V∞ 1.5–3.5 km/s, though these sit near the flyby-bend ceiling, i.e.
*low* margin), with the high-V∞ Io–Callisto / Io–Ganymede geometries (V∞ 5.4–7.7 km/s, large
bend margin but physically hostile) at the other end. **The ranking is moot**: literature clearance
kills the whole population before robustness matters. No candidate is worth the gauntlet:

- The five triple-cycler-covered pairs are known repeated-moon territory (Liang CGE 2024, Hernandez
  IEG 2017) and known R-S double-cycler territory (three of them explicitly enumerated in R-S
  Table 1). Pushing them to a real-ephemeris gauntlet would at best re-derive a published family.
- Io–Callisto, the sole structural-clear, fails the higher bar on independent grounds: (i) it is
  generated by R-S 2009's published enumerative method; (ii) it is discussed in the review
  literature as dynamically awkward; (iii) its V∞ is 5.4–7.3 km/s with Io buried in Jupiter's
  radiation belt and a 4.46× semi-major-axis ratio — the least attractive Galilean pair, with no
  reason to expect it survives inclination-closure + eccentricity-kill + real-ephemeris continuation
  as a useful cycler. There is **no strong reason to believe it survives the pipeline**, so per
  #577 scope it is not novelty-claimable regardless of the structural clear.

This does **not** unstamp #501's real-ephemeris "0/3072 closed, clean empty-region map"; the two
explore structurally different parameter space (the #576 Fable review's Q1 resolution stands).

## 4. What IS worth doing (scoped as #578, non-pipeline)

The actionable output is not a gauntlet run but a **corpus-accuracy fix + registry stamp**: add the
R-S 2009 double-cycler paper (10.2514/1.36610) to `KNOWN_CORPUS` as a `repeated-moon` anchor over
the Galilean body set so the matcher correctly flags this territory (closing the Io–Callisto false
clear), and register the Galilean symmetric-closure region as literature-covered known-adjacent in
the negative-results registry. This prevents the same false-negative from recurring on any future
Jovian double-cycler screen. Requires a Fable second-opinion before dispatch (this chain is 9-for-9
on real Fable catches).

---

## Bottom line

The #576 population is real idealized-model physics and honestly instrumented, but it is
comprehensively pre-published: R-S 2009 double-cyclers (same model, same method, three of six pairs
explicitly enumerated), Liang 2024 / Hernandez 2017 triple-cyclers (the other structural coverage),
and the Strange/Petropoulos/Niehoff pump-tour corpus (all pairs, body-set level). **0/36 novel;
full stop on the pipeline; fix the corpus gap and stamp the territory (#578).**

## Addendum 2026-10-03 (`#869`): the review cited above has now been read

Section 1 leaned on "a 2024 review ... [that] surveys the whole double-cycler field" and on
Io-Callisto being "*discussed* in the literature as a dynamically awkward pair". Both came from
live-search snippets; the paper was not in hand. It now is (Yang, Hu, Bai & Li, *Space: Science &
Technology* 3:0036, digest `2026-10-03-digest-yang-2023-jovian-system-trajectory-review.md`):

- It was published on 23 May 2023, not in 2024.
- It does not survey the double-cycler field. Cyclers get one bullet paragraph (p. 5), which says
  Russell and Strange "first investigated the double-cyclers for exploring the Galilean moons.
  Specifically, the cyclers of Ganymede-Io, Ganymede-Europa, Ganymede-Callisto, and
  Europa-Ganymede are searched", then cites Lynam-Longuski and Hernandez-Jones-Jesick for the
  Io-Europa-Ganymede triple cyclers. No cycler numbers are tabulated. Liang et al. is not cited.
- It mentions no Io-Callisto or Europa-Callisto cycler at all. The quoted phrase about Callisto
  ("does not participate in the Laplace resonance ... presents unique challenges") is not in the
  paper; the nearest text is a satellite-aided CAPTURE table note, "Callisto is not in Laplace
  resonance but near-resonances can be found".

**What this changes.** Nothing for the five pairs covered by Russell-Strange's own enumeration or
by the triple cyclers. For Io-Callisto, ground (ii) of Section 3 ("discussed in the review
literature as dynamically awkward") is unsupported and is withdrawn; grounds (i) (generated by the
published enumerative method) and (iii) (physically unattractive: V-infinity 5.4-7.3 km/s, Io deep
in the radiation belt) stand. The review independently CONFIRMS that Russell and Strange searched
only four pairs, so Io-Callisto and Europa-Callisto are pairs no located source enumerates.

**Open question for the owner (not decided here).** The novelty policy decided after this note
(`#875`, spec 16.4) keeps `known-class-member` for "a computed member of a class the source authors
themselves searched -- same system, same architecture, a member they did not happen to tabulate",
and says the 0/36 ruling stands because these are same-system members of Russell-Strange's Jovian
class. Whether a PAIR the authors never searched (Io-Callisto) is "a member they did not happen to
tabulate" or an unsearched sub-class is not settled by that text. It matters only if someone
wants to pursue Io-Callisto, which ground (iii) argues against.

**Owner decision 2026-10-03:** Io-Callisto is NOT pursued. The six `#576` closures (three
geometries read from each end; legs 7.92 / 8.91 / 9.90 d; V-infinity 5.4-7.3 km/s) stay
uncatalogued. Under `#875` as decided they would be same-system members of Russell and Strange's
Jovian class, so the open question above is left unanswered rather than ruled on. Reopen only
for a mission-driven reason or if the owner later rules that an unenumerated pair in a treated
system is an unsearched sub-class (the same question arises for an untreated base moon at Uranus,
`#882` item 8).
