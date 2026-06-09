# Self-seeding longitude-rendezvous construction — DESIGN (task #173)

**Date:** 2026-06-08
**Type:** design note (no production code; reading-only on src/data). Companion
executable plan: `docs/superpowers/plans/2026-06-08-self-seeding-construction.md`.
**Status:** DESIGN — not built. The S1L1-self-seed validation gate (§5) is the
cheapest-decisive go/no-go and MUST pass before any unsourced member is attempted.

---

## 0. The problem in one paragraph

The proven S1L1 pipeline (#166→#167) closes `russell-ch4-4.991gG2` on DE440 with a
corrected topology — `E → g(E-E free return, sub-Mars) → E flyby → G(E-M-E transit,
true-longitude rendezvous with DE440 Mars) → E` — but it is **seeded from Russell
Appendix-C #83's printed per-leg real-eph state block**. That block is what supplies
the one ingredient every from-scratch attempt lacked: the *longitude-rendezvous
constraint* — the departure epoch at which real DE440 Mars actually sits at the
spacecraft's encounter point (Mars at ecliptic longitude 201.0° on 2027-06-13 for
S1L1 leg 2). The #170 batch then proved that **only 2 of the 9 russell-ch4 parents
have a sourced App-C block at all**; the other 7, plus the ~194 ocampo members, have
NO published real-eph seed. To reach those rows the printed seed must be **replaced
by a search**: given only the row's descriptor (a, e, ToF, ψ, the g/G topology),
*find* the departure epoch where real DE440 Mars is at the encounter longitude, then
run the existing #167 closure + #164 continuation + the independent n-body confirm.

**This design scopes that search and — critically — makes the historical off-family
failure (2026-06-04) a first-class possible outcome.** The honest prior is that this
MAY still land off-family even with the new ingredients. The design's job is to (a)
specify the epoch/longitude search, (b) name the on-family guard and reckon with why
it might still fail, (c) maximise reuse, and (d) put the S1L1-known-answer gate first
so a wrong method is caught for the price of one row, not 194.

---

## 1. The longitude-rendezvous epoch search

### 1.1 What the App-C seed actually supplied (so we know what to replace)

From #166, the App-C block gave, per Mars-transit (G) leg: a shared epoch + a
per-leg `time-start` (⇒ the **Earth-departure epoch**) and a `v∞`-vector. Reconstructed
on DE440 (`v_sc = v_planet(DE440) + v∞`, then Kepler-Sun to the published arrival
epoch), the spacecraft lands on real Mars to 1.7 km at the published v∞ to 4 dp, at
the same heliocentric longitude as Mars (Δlon 4e-7°). So the seed supplied **two
coupled things**: (i) the **encounter epoch** (when), and (ii) the **departure v∞
vector** (which heliocentric ellipse). The longitude rendezvous is the constraint
"real Mars is *there* at the arrival epoch." For an unsourced row we have NEITHER the
epoch nor the v∞ vector — only the coplanar descriptor (a, e of each arc, the two
ToFs, the v∞ magnitudes, ψ).

### 1.2 The search, concretely (two stages)

**Stage A — descriptor → coplanar geometry (already solved, reused).** The #163/#164
two-arc free-return chain (`free_return_chain.py` + `continuation_chain.py`) already
takes the descriptor `(a₁,e₁,a₂,e₂, n_rev, ToFs)` and produces a coplanar-then-DE440
two-arc geometry whose emerged v∞ matches the row's own anchors and whose ToFs match
the descriptor. This fixes the **shape** of the G arc (its (a,e), hence its
heliocentric encounter radius, the encounter true anomaly, and the **encounter
heliocentric longitude RELATIVE to the Earth-departure direction**). What it does NOT
fix is the **absolute epoch** — where in real time the departure happens so that real
Mars is at that encounter longitude. This is exactly the axis #164 closed (ToF
quantization) and the axis #165 found still open (longitude/phase vs the true
ephemeris).

**Stage B — the new piece: a longitude-rendezvous epoch solve.** The G arc's geometry
gives a function `lon_sc_encounter(t_depart)` = the heliocentric ecliptic longitude
at which the spacecraft crosses Mars's radius, as a function of departure epoch
(the arc shape is epoch-independent in the circular limit; on DE440 it drifts slowly
with Earth's true position at departure). We need the epoch where this equals real
Mars's longitude at the arrival epoch:

```
residual_lon(t_depart) = lon_sc_encounter(t_depart) − lon_Mars_DE440(t_depart + ToF_G)
```

This is a scalar root in one variable (`t_depart`) over **one synodic period** of the
Earth–Mars system (~2.135 yr). Two ways to do it, in increasing fidelity:

1. **Phase/longitude scan (the cheap, robust opener).** Sweep `t_depart` across one
   synodic period in coarse steps (e.g. 10 d), evaluate `residual_lon`, bracket every
   sign change, then Brent/bisection-refine each bracket to a root. This yields ALL
   the synodic-phase candidates (there are typically 1–2 per synodic period where the
   longitudes line up). It is a scan over a *known* periodic structure, NOT a free
   high-dimensional optimisation — this is the structural difference from the
   2026-06-04 failure (§2).

2. **Lambert/shooting node on true DE440 Mars (the on-family refinement).** At each
   bracketed root, refine with a Lambert solve: from real Earth at `t_depart` to real
   DE440 **Mars position** at `t_depart + ToF_G` (the TRUE 3-D position, not just the
   radius), solving for the transfer that matches the descriptor's revolution count
   and whose departure v∞ magnitude is closest to the row's Earth anchor. The Lambert
   arc *by construction* arrives at Mars's true position — longitude rendezvous is
   then automatic; the residual that remains is the **v∞-vs-anchor** mismatch and the
   **ToF-vs-descriptor** mismatch (the on-family test, §1.3). The g arc is then
   back-solved to close the cycle to Earth (it carries NO Mars constraint — the
   topology fix).

The deliverable of the search per row is: an Earth-departure epoch + departure v∞
vector that puts the spacecraft on **true DE440 Mars** at the encounter — i.e. exactly
the (epoch, v∞-vector) the App-C block printed for S1L1, now FOUND rather than read.

### 1.3 The residual (every binding constraint — the discipline rule)

Per the orbit-closure discipline ("all binding constraints in the residual; 'it
closed!' on a subset is the danger signal"), the self-seeding residual MUST contain
every one of these, and the on-family test is that ALL are small at the SAME epoch:

| term | meaning | why it must be in |
|---|---|---|
| `residual_lon` | SC encounter longitude − DE440 Mars longitude | the constraint App-C supplied; #165's 110° miss was this term omitted |
| `vinf_M − anchor_M` | emerged Mars v∞ vs the row's real-eph Mars anchor | the #163 spurious-collapse trap is a v∞-only residual; needs the *real-eph* anchor (breathing), not coplanar |
| `vinf_E − anchor_E` | emerged Earth-departure v∞ vs Earth anchor | family identity |
| `ToF_G − descriptor_G` | transit time vs descriptor | the #164 quantization axis; binds n_rev |
| `ToF_g − descriptor_g` | g-arc time vs descriptor | keeps the two arcs distinct (no single-ellipse collapse) |
| `g-arc closest-Mars ≫ band` | g arc stays sub-Mars (aphelion < ~1.4 AU) | the topology guard; #164/#165 violated this |
| intermediate Earth-flyby bend ≤ max | flyby feasibility | a closure with an infeasible turn is not a cycler |

A self-seed is ON-FAMILY only if `residual_lon → 0` **simultaneously with** v∞ at both
bodies near the row's *real-eph* anchors AND both ToFs in band AND the g arc sub-Mars
AND the flyby bend-feasible — at one epoch. Any subset closing alone is the danger
signal, not a result.

---

## 2. The on-family guard — and why it might STILL fail (honest)

### 2.1 The 2026-06-04 failure mode (what we must beat)

The blocker records: free epoch scans / `optimise_cell_ephemeris` driven to hit the
v∞ anchors **landed off-family** (e.g. V_E≈26.5, V_M≈23.6; or "best Vinf-continuity
epoch lands off-family Earth Vinf ~34"). Mechanism: the resolver minimised *summed*
leg-v∞ mismatch over a high-dimensional free search; the multi-rev leg dominated, the
front legs were mis-phased, and the minimiser settled in a degenerate basin that
satisfied the scalar metric while being physically the wrong cycler. The lesson: a
free real-eph multiple-shooting optimisation, with v∞-magnitude-only objective, has
many spurious low-residual basins and no guarantee of selecting the published family.

### 2.2 The three NEW ingredients that should beat it

1. **The corrected topology constrains the search.** The g arc carries NO Mars
   constraint and is *forced* sub-Mars (aphelion < 1.4 AU); only the G arc targets
   Mars. The 2026-06-04 search let both arcs float against Mars (or used the wrong
   E-M-E-E cell), which is most of why it wandered. Halving the Mars-constrained DOF
   removes the dominant off-family basins.

2. **The longitude constraint is EXPLICIT, not emergent.** The Lambert/shooting node
   (Stage B option 2) targets Mars's TRUE 3-D position. An off-family basin that
   crosses Mars's *radius* at the wrong longitude (the #165 artifact) is now
   *infeasible by construction* — the Lambert arc cannot arrive anywhere but at real
   Mars. The free search could trade longitude error for v∞ fit; the explicit
   longitude target forbids that trade.

3. **The scan is over a KNOWN periodic structure (the synodic period), not a free
   optimisation.** `residual_lon(t_depart)` is a smooth scalar with O(1) roots per
   synodic period. Bracketing-then-refining every sign change ENUMERATES the
   synodic-phase candidates instead of gradient-descending from one seed into the
   nearest (possibly degenerate) basin. We then test each enumerated candidate against
   the full residual (§1.3). This is the structural antidote to "the minimiser picked
   a wrong basin": we don't trust a minimiser to find the family — we enumerate the
   phase candidates and check each.

4. **The anchor-respecting residual + real-eph (breathing) target.** Reusing the
   #163/#164 anchor-respecting objective (emerged − sourced → 0), against the
   *real-eph* anchors (the breathing v∞ from the row's own model, NOT coplanar
   5.10/3.05), is what kept #137/#163 on-family where dsm_leg's ΔV-min floated off.

### 2.3 Why it might STILL go off-family (the honest risk register)

The design must NOT claim these ingredients guarantee success. Concrete failure modes
that survive the new ingredients:

- **Multiple synodic-phase roots, several plausibly on-family.** The scan may return
  2+ epochs each with small `residual_lon`; if more than one also has near-anchor v∞,
  family selection is ambiguous and we have no sourced tiebreak (S1L1 *did* have
  one — the App-C epoch). Mitigation: the S1L1 gate (§5) tells us which root the
  published seed corresponds to and whether the scan even surfaces it; for unsourced
  rows we report ALL on-family-looking roots, not a single forced pick.
- **The real-eph v∞ "breathes" — but we don't know the per-cycle envelope for an
  unsourced row.** For S1L1 the App-C table gave the per-leg breathing (3.2–8.0). For
  an unsourced row we have only the coplanar descriptor v∞. If the row's real-eph v∞
  departs far from the coplanar value (as S1L1's does: coplanar 5.10 vs real-eph
  3.2–8.0), the anchor we test against is itself uncertain — a correct closure could
  be mis-judged as off-anchor (the inverse-#164 trap). Mitigation: test v∞ against a
  *band* around the coplanar anchor wide enough to admit real-eph breathing, and
  report the emerged v∞ as evidence rather than gating hard on the coplanar number.
- **The descriptor → (a,e) map may not host a real-eph longitude solution at all.**
  #165's epoch scan over a full synodic period found a 0.234 AU best-phase floor for
  the (then wrong-topology) S1L1 arc — a warning that an arc shape can be
  longitude-unsatisfiable against the true ephemeris. A corrected-topology unsourced
  row may simply have no epoch where its G-arc shape rendezvous with real Mars at the
  anchor v∞. That is a clean EMPTY-SET negative, not a bug.
- **Multi-arc rows whose g arc is itself constrained** (the S1L1 g is benign
  sub-Mars; some ocampo members may have a g that must thread a specific phase) could
  reintroduce the coupling we removed.

**Verdict on the prior:** plausible-but-unproven. The honest expectation is that the
S1L1 self-seed gate (§5) is genuinely informative — it could fail, and a clean failure
there KILLS the method before it touches an unsourced row. We do NOT claim success
until that gate passes.

---

## 3. Reuse map — what is NEW vs reused

| component | status | source |
|---|---|---|
| descriptor → two-arc (a,e,n_rev) coplanar geometry | **REUSE** | `free_return_chain.py` (#163) |
| circular→DE440 homotopy (e-ramp breaks ToF quantization) | **REUSE** | `continuation_chain.py` (#164) |
| corrected topology (g sub-Mars, G transit, one Mars enc/cycle) | **REUSE** | `s1l1_corrected.py` topology (#167) |
| per-leg reconstruction recipe `v_sc = v_planet(DE440)+v∞` | **REUSE** | `s1l1_corrected.build_seeded_arcs` (#167) |
| g-arc sub-Mars clearance check | **REUSE** | `s1l1_corrected.g_arc_clearances` (#167) |
| independent n-body confirm (REBOUND/IAS15 over DE440, 3-SOI band) | **REUSE** (gate) | `tests/nbody/test_s1l1_corrected_nbody.py` pattern (#167) |
| continuous-from-one-seed TCM (maintenance budget) | **REUSE** | `s1l1_corrected.continuous_chain` (#169) |
| **`residual_lon(t_depart)` + synodic-period bracket-and-refine scan** | **NEW** | this design §1.2 |
| **Lambert/shooting node onto TRUE DE440 Mars position** | **NEW** (refine) | this design §1.2 (a Lambert primitive may already exist — VERIFY) |
| **on-family multi-term gate at a single epoch (§1.3)** | **NEW** (assembly) | this design §1.3 |
| **enumerate-and-test phase candidates (vs free-optimise)** | **NEW** (method) | this design §2.2 |

So the genuinely new surface is small and bounded: a one-variable longitude residual,
a synodic-period scan that brackets its roots, an optional Lambert refinement onto
true Mars, and the assembly of the full on-family gate. Everything downstream of "we
have a (epoch, v∞-vector) seed" is the existing #167 pipeline unchanged. **This is
why the S1L1 gate is cheap:** we are testing one new search stage feeding a proven tail.

---

## 4. Why S1L1 is the right known-answer test

S1L1 is the ONLY row where we have BOTH (a) the App-C printed seed (the answer key)
and (b) a CONFIRMED #167 closure to recover. So it is the unique row on which we can:

1. run the self-seeding search **without ever reading the App-C block** (pure
   descriptor → epoch/v∞), and
2. check the FOUND seed against the App-C answer key (epoch within tolerance of
   2026-12-15; Mars encounters in-band at the breathing v∞; the same 201.0° longitude
   rendezvous).

If the search recovers the App-C geometry from scratch, the method is validated on a
known answer and we have earned the right to apply it to an unsourced row. If it does
NOT — if the scan misses the App-C epoch, or lands a different (off-family) basin —
then we have learned the method is not yet trustworthy, for the price of one row, and
the App-C answer key tells us exactly HOW it failed (wrong root selected? longitude
residual never zeroes? v∞ off-anchor?). Either way the gate is decisive.

---

## 5. THE VALIDATION GATE (the near-perfect one) — S1L1 self-seed

**Acceptance target (sourced):** S1L1's OWN App-C-confirmed geometry from #167 — this
is a same-model golden (the EXPECTED side traces to Russell App-C #83, never to our
own search output):

| quantity | App-C answer key (#166/#167) | self-seed must recover |
|---|---|---|
| Earth-departure epoch (G leg 2) | 2026-12-15 | within ±1 synodic-phase tolerance (target ≤ a few days; report exact) |
| Mars-flyby epoch | 2027-06-13 | within tolerance |
| Mars longitude at flyby | 201.0° | rendezvous Δlon < ~0.5° |
| Mars v∞ (leg 2) | 5.248 (breathing 3.2–8.0 over 7 cycles) | emerged v∞ in the real-eph band |
| 7 Mars encounters | all in 3-SOI band (≈0.0116 AU) | independent n-body confirm, SAME band, NOT loosened |

**The gate is binary and decisive:**

- **PASS** — the synodic scan surfaces the App-C epoch (≤ tolerance), the on-family
  residual (§1.3) is small there, and the recovered seed drives the #167 tail to the
  same 7-encounters-in-band confirmation. ⇒ method validated; proceed to ONE unsourced
  member.
- **FAIL** — the scan misses the App-C epoch, or lands a different basin, or the
  longitude residual never zeroes at an anchor-consistent epoch. ⇒ method NOT
  validated; report the failure mode (the App-C key tells us which), do NOT touch any
  unsourced row. A clean FAIL here is a SUCCESSFUL outcome of this task (it kills a
  bad method cheaply).

**No catalogue writeback at any point in this design's scope.** S1L1 is already a
recommended V3 from #167; this task does not change its level — it tests a *method*.

---

## 6. After the gate — prove-on-one, never batch-of-194

If and only if §5 PASSES: apply to exactly ONE unsourced member (a russell-ch4 row
from the 7, e.g. a near-ballistic coplanar one, or one ocampo member). Run the full
search + #167 tail + independent n-body confirm + continuous-TCM budget. Three honest
verdicts, each a first-class result:

- **CONFIRMED** — self-found seed, encounters in-band at real-eph v∞, TCM ≤ budget ⇒
  the first unsourced row reachable without a published seed (recommend V3, held).
- **PARTIAL** — encounters in-band but TCM over budget (the #170 powered-parent
  pattern) ⇒ honest PARTIAL, no writeback.
- **OFF-FAMILY / EMPTY-SET** — no epoch hosts the longitude rendezvous at the anchor
  v∞ (the historical failure, now a *quantified* negative with the corrected topology)
  ⇒ a clean negative; record it. Per the negative-results-registry memory, an
  EMPTY-SET here is method-versioned and feeds the anti-catalogue.

Only after one member is decisively classified does a batch become a question — and
that is a separate task with its own brief. **Prove on one before any batch.**

---

## 7. Honest expected outcome (the brief's rule 4)

The single most likely outcome, stated plainly: **the S1L1 self-seed gate is a real
coin-flip, and even a PASS there may not transfer to unsourced rows.** The 2026-06-04
off-family failure is the base rate for from-scratch real-eph family selection; the
three new ingredients (corrected topology, explicit longitude target, enumerate-don't-
optimise) are *reasons to expect better*, not a proof. The design is built so that:

- a clean FAIL at §5 is a first-class, cheap, decisive negative (method killed before
  any unsourced row);
- a PASS at §5 earns exactly ONE unsourced attempt, whose OFF-FAMILY/EMPTY-SET outcome
  is equally first-class;
- nothing is claimed CONFIRMED, and nothing is written to the catalogue, on the
  finding-search's say-so — the independent n-body confirm and the same-model golden
  remain mandatory throughout.

The success criterion for *this task's method* is not "it closed an unsourced row" —
it is "we now know, decisively, whether self-seeding can reach the unsourced rows, and
the S1L1 known-answer gate told us so for the price of one row."
